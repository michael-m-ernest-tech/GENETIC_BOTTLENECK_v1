#!/usr/bin/env python3
"""Microsatellite heterozygosity and mutation-drift equilibrium analysis.

CLI:
  python3 bottleneck_analyzer.py genotypes.csv --model TPM --ne 50 \
      --mu 0.001 --replicates 500 --out-prefix results

Web UI:
  python3 bottleneck_analyzer.py --web
  Then open the local URL printed by the program.

Input CSV:
  First column: sample ID
  Remaining columns: allele pairs per locus, e.g. L01_A1,L01_A2,L02_A1,L02_A2
  Alleles must be integer repeat sizes. Empty pairs and NA/. values are treated as
  missing; a partially missing pair is excluded for that locus.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import random
import statistics
import sys
from email import policy
from email.parser import BytesParser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


MISSING = {"", "na", "nan", ".", "null", "none"}


def parse_genotypes(text: str) -> dict:
    """Parse paired allele columns and calculate per-locus Ho and He."""
    reader = csv.reader(io.StringIO(text))
    try:
        header = next(reader)
    except StopIteration:
        raise ValueError("The CSV is empty.")

    header = [x.strip() for x in header]
    if len(header) < 3 or (len(header) - 1) % 2:
        raise ValueError(
            "Expected a sample ID column followed by an even number of allele columns."
        )
    if any(not x for x in header):
        raise ValueError("Every CSV column must have a header.")

    locus_names = []
    for i in range(1, len(header), 2):
        left, right = header[i], header[i + 1]
        # For standard paired names, use the shared prefix as the locus name.
        name = left
        for suffix in ("_A1", "_Allele1", "_1"):
            if left.lower().endswith(suffix.lower()):
                name = left[: -len(suffix)]
                break
        if name == left:
            name = f"Locus{i // 2 + 1}"
        locus_names.append(name or f"Locus{i // 2 + 1}")

    allele_pairs = [[] for _ in locus_names]
    sample_ids = set()
    for line_num, row in enumerate(reader, start=2):
        if not row or not any(cell.strip() for cell in row):
            continue
        if len(row) != len(header):
            raise ValueError(
                f"Row {line_num} has {len(row)} fields; expected {len(header)}."
            )
        sample_id = row[0].strip()
        if not sample_id:
            raise ValueError(f"Row {line_num} has a blank sample ID.")
        if sample_id in sample_ids:
            raise ValueError(f"Duplicate sample ID {sample_id!r} on row {line_num}.")
        sample_ids.add(sample_id)

        for locus_i, col in enumerate(range(1, len(row), 2)):
            a, b = row[col].strip(), row[col + 1].strip()
            a_missing, b_missing = a.lower() in MISSING, b.lower() in MISSING
            if a_missing and b_missing:
                continue
            if a_missing != b_missing:
                continue  # Exclude partial genotypes for this locus.
            try:
                allele_pairs[locus_i].append((int(a), int(b)))
            except ValueError:
                raise ValueError(
                    f"Non-integer allele at row {line_num}, locus "
                    f"{locus_names[locus_i]!r}: {a!r}, {b!r}."
                )

    if not sample_ids:
        raise ValueError("No sample rows were found.")

    loci = []
    for name, genotypes in zip(locus_names, allele_pairs):
        if not genotypes:
            raise ValueError(f"No complete genotypes found for {name}.")
        counts = {}
        heterozygotes = 0
        for a, b in genotypes:
            counts[a] = counts.get(a, 0) + 1
            counts[b] = counts.get(b, 0) + 1
            heterozygotes += a != b

        n_genes = 2 * len(genotypes)
        he = 1.0 - sum((count / n_genes) ** 2 for count in counts.values())
        ho = heterozygotes / len(genotypes)
        loci.append({
            "locus": name,
            "n": len(genotypes),
            "allele_count": len(counts),
            "allele_frequencies": {
                str(allele): count / n_genes for allele, count in sorted(counts.items())
            },
            "Ho": ho,
            "He": he,
            "difference": ho - he,
        })

    return {"sample_count": len(sample_ids), "loci": loci}


def wilcoxon_greater(differences: list[float], alternative: str = "Ho - He > 0") -> dict:
    """One-sided Wilcoxon signed-rank test: H1 is median difference > 0.

    Uses an exact sign-flip distribution for <=25 nonzero paired differences,
    including tied absolute ranks. Uses a tie-corrected normal approximation above
    that size. Zero differences are omitted, as in the usual signed-rank test.
    """
    values = [float(x) for x in differences if abs(float(x)) > 1e-12]
    n = len(values)
    if n == 0:
        return {
            "alternative": alternative,
            "n_nonzero_loci": 0,
            "W_positive": 0.0,
            "p_value": 1.0,
            "method": "No nonzero differences",
            "median_difference": 0.0,
        }

    # Round magnitudes so float noise does not break genuine ties.
    ordered = sorted((round(abs(v), 12), i) for i, v in enumerate(values))
    ranks = [0.0] * n
    tie_sizes = []
    i = 0
    while i < n:
        j = i + 1
        while j < n and ordered[j][0] == ordered[i][0]:
            j += 1
        average_rank = ((i + 1) + j) / 2.0  # 1-based ranks
        tie_sizes.append(j - i)
        for k in range(i, j):
            ranks[ordered[k][1]] = average_rank
        i = j

    w_positive = sum(rank for rank, value in zip(ranks, values) if value > 0)

    if n <= 25:
        # Doubled average ranks are integers, enabling exact subset-sum counts.
        scaled = [int(round(rank * 2)) for rank in ranks]
        observed = int(round(w_positive * 2))
        total = sum(scaled)
        ways = [0] * (total + 1)
        ways[0] = 1
        reachable = 0
        for rank in scaled:
            for score in range(reachable, -1, -1):
                if ways[score]:
                    ways[score + rank] += ways[score]
            reachable += rank
        p_value = sum(ways[observed:]) / float(2 ** n)
        method = "Exact one-sided sign-flip distribution"
    else:
        mean = n * (n + 1) / 4.0
        variance = n * (n + 1) * (2 * n + 1) / 24.0
        variance -= sum(t**3 - t for t in tie_sizes) / 48.0
        if variance <= 0:
            p_value, method = 1.0, "Normal approximation (degenerate variance)"
        else:
            # Continuity correction for the upper tail.
            z = (w_positive - mean - 0.5) / math.sqrt(variance)
            p_value = 0.5 * math.erfc(z / math.sqrt(2.0))
            method = "Tie-corrected normal approximation"

    return {
        "alternative": alternative,
        "n_nonzero_loci": n,
        "W_positive": w_positive,
        "p_value": min(1.0, max(0.0, p_value)),
        "method": method,
        "median_difference": statistics.median(values),
    }


def poisson(rng: random.Random, lam: float) -> int:
    """Draw a Poisson variate using Knuth's exact method in small chunks."""
    if lam <= 0:
        return 0
    total = 0
    while lam > 0:
        part = min(lam, 20.0)
        threshold = math.exp(-part)
        product = 1.0
        count = 0
        while product > threshold:
            count += 1
            product *= rng.random()
        total += count - 1
        lam -= part
    return total


def branch_mutation_shift(
    rng: random.Random,
    branch_generations: float,
    mu: float,
    model: str,
    p_step: float,
    geometric_p: float,
) -> int:
    """Return repeat-size change along one coalescent branch."""
    events = poisson(rng, mu * branch_generations)
    shift = 0
    for _ in range(events):
        sign = -1 if rng.random() < 0.5 else 1
        if model == "SMM" or rng.random() < p_step:
            shift += sign
        else:
            # TPM multi-step sizes are at least 2, with a geometric tail.
            size = 2
            while rng.random() > geometric_p:
                size += 1
            shift += sign * size
    return shift


def simulate_equilibrium_locus(
    gene_copies: int,
    ne: int,
    mu: float,
    model: str,
    p_step: float,
    geometric_p: float,
    rng: random.Random,
) -> tuple[int, float]:
    """Simulate one neutral equilibrium sample using a diploid Kingman coalescent.

    Pairwise coalescence rate is 1/(2Ne) per generation. Mutations are placed
    on tree branches at rate mu per gene copy per generation.
    """
    # Node entries are (height in generations, child node IDs or None).
    nodes = [(0.0, None) for _ in range(gene_copies)]
    active = list(range(gene_copies))
    height = 0.0

    while len(active) > 1:
        k = len(active)
        height += rng.expovariate(k * (k - 1) / (4.0 * ne))
        i, j = sorted(rng.sample(range(k), 2))
        child_a, child_b = active[i], active[j]
        parent = len(nodes)
        nodes.append((height, (child_a, child_b)))
        active.pop(j)
        active.pop(i)
        active.append(parent)

    root = active[0]
    leaf_alleles = []
    stack = [(root, 0)]
    while stack:
        node_id, allele = stack.pop()
        node_height, children = nodes[node_id]
        if children is None:
            leaf_alleles.append(allele)
            continue
        for child in children:
            child_height, _ = nodes[child]
            branch_length = node_height - child_height
            child_allele = allele + branch_mutation_shift(
                rng, branch_length, mu, model, p_step, geometric_p
            )
            stack.append((child, child_allele))

    counts = {}
    for allele in leaf_alleles:
        counts[allele] = counts.get(allele, 0) + 1
    he = 1.0 - sum((count / gene_copies) ** 2 for count in counts.values())
    return len(counts), he


def quantile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    pos = (len(ordered) - 1) * q
    low = int(pos)
    high = min(low + 1, len(ordered) - 1)
    frac = pos - low
    return ordered[low] * (1 - frac) + ordered[high] * frac


def equilibrium_comparison(
    loci: list[dict],
    ne: int,
    mu: float,
    model: str,
    replicates: int,
    seed: int,
    p_step: float,
    geometric_p: float,
) -> list[dict]:
    """Compare each observed He to simulated equilibrium He at similar allele count.

    Exact allele-count matches are used when at least 20% of simulations match.
    Otherwise, simulations with the nearest allele count are used and the result
    labels that fallback. The Monte Carlo upper-tail p-value is P(He_sim >= He_obs);
    a small value means the observed He is high relative to this model baseline.
    """
    rng = random.Random(seed)
    results = []
    for locus in loci:
        sims = [
            simulate_equilibrium_locus(
                2 * locus["n"], ne, mu, model, p_step, geometric_p, rng
            )
            for _ in range(replicates)
        ]
        exact = [he for k, he in sims if k == locus["allele_count"]]
        if len(exact) >= max(10, int(replicates * 0.2)):
            baseline = exact
            conditioning = "Exact allele-count match"
        else:
            distance = min(abs(k - locus["allele_count"]) for k, _ in sims)
            baseline = [he for k, he in sims if abs(k - locus["allele_count"]) == distance]
            conditioning = f"Nearest simulated allele count (distance {distance})"

        upper_tail = (1 + sum(value >= locus["He"] for value in baseline)) / (
            len(baseline) + 1
        )
        results.append({
            **locus,
            "equilibrium_mean_He": statistics.mean(baseline),
            "equilibrium_median_He": statistics.median(baseline),
            "equilibrium_95_low": quantile(baseline, 0.025),
            "equilibrium_95_high": quantile(baseline, 0.975),
            "equilibrium_upper_tail_p": upper_tail,
            "equilibrium_n": len(baseline),
            "conditioning": conditioning,
        })
    return results


def analyze(
    text: str,
    *,
    model: str = "TPM",
    ne: int = 50,
    mu: float = 0.001,
    replicates: int = 500,
    seed: int = 42,
    p_step: float = 0.9,
    geometric_p: float = 0.5,
) -> dict:
    parsed = parse_genotypes(text)
    wilcoxon = wilcoxon_greater([x["difference"] for x in parsed["loci"]])
    comparisons = equilibrium_comparison(
        parsed["loci"], ne, mu, model, replicates, seed, p_step, geometric_p
    )
    # Classic bottleneck-style test: observed He vs equilibrium He (Heq) per locus.
    wilcoxon_eq = wilcoxon_greater(
        [x["He"] - x["equilibrium_mean_He"] for x in comparisons],
        alternative="He - Heq > 0",
    )
    return {
        "sample_count": parsed["sample_count"],
        "locus_count": len(comparisons),
        "wilcoxon": wilcoxon,
        "wilcoxon_he_vs_eq": wilcoxon_eq,
        "model": {
            "name": model,
            "Ne": ne,
            "mutation_rate": mu,
            "replicates": replicates,
            "seed": seed,
            "TPM_one_step_probability": p_step,
            "TPM_geometric_probability": geometric_p,
            "baseline_method": "Neutral diploid Kingman coalescent with mutations",
        },
        "per_locus": comparisons,
        "interpretation": (
            "Wilcoxon (Ho - He > 0) tests heterozygote excess relative to Hardy-Weinberg "
            "expectation. Wilcoxon (He - Heq > 0) is the bottleneck-style test: a small "
            "p-value means He is higher than expected at mutation-drift equilibrium. "
            "A small per-locus equilibrium upper-tail p-value indicates the observed He "
            "exceeds the simulated baseline at similar allele count. "
            "These results depend on the chosen Ne, mutation rate, and mutation model."
        ),
    }


def csv_report(result: dict) -> str:
    output = io.StringIO()
    fields = [
        "locus", "n", "allele_count", "Ho", "He", "difference",
        "equilibrium_mean_He", "equilibrium_median_He",
        "equilibrium_95_low", "equilibrium_95_high",
        "equilibrium_upper_tail_p", "equilibrium_n", "conditioning",
        "allele_frequencies",
    ]
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for row in result["per_locus"]:
        clean = dict(row)
        clean["allele_frequencies"] = json.dumps(clean["allele_frequencies"])
        writer.writerow({key: clean[key] for key in fields})
    return output.getvalue()


WEB_PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Microsatellite Bottleneck Analyzer</title>
<style>
:root{color-scheme:light;--ink:#182230;--muted:#617083;--line:#dce3eb;--blue:#2357d6;--pale:#f4f7fb}
*{box-sizing:border-box}body{margin:0;background:#f2f5f9;color:var(--ink);font:15px/1.5 system-ui,Segoe UI,sans-serif}
main{max-width:1180px;margin:32px auto;padding:0 20px}.hero{background:#14253c;color:white;padding:28px 30px;border-radius:18px}
h1{margin:0 0 4px;font-size:27px}.hero p{margin:0;color:#c7d3e3}.grid{display:grid;grid-template-columns:340px 1fr;gap:18px;margin-top:18px}
.card{background:white;border:1px solid var(--line);border-radius:14px;padding:20px;box-shadow:0 5px 20px #18304c0a}
h2{font-size:18px;margin:0 0 12px}.field{margin:12px 0}.field label{display:block;font-weight:650;margin-bottom:5px}
input,select{width:100%;padding:10px;border:1px solid #bdc8d5;border-radius:8px;background:white;color:var(--ink)}
button{border:0;background:var(--blue);color:white;font-weight:700;padding:11px 15px;border-radius:9px;cursor:pointer}
button.secondary{background:#eaf0fb;color:#1d438f;margin-left:6px}button:disabled{opacity:.6;cursor:wait}
.hint,.subtle{font-size:13px;color:var(--muted)}.hint{margin-top:3px}.error{color:#a52626;font-weight:600;margin-top:10px}
.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}.metric{background:var(--pale);padding:13px;border-radius:10px}.metric b{display:block;font-size:21px}.metric span{font-size:12px;color:var(--muted)}
.test{margin-top:14px;padding:14px;border-left:4px solid var(--blue);background:#f7f9fd;border-radius:8px;white-space:pre-line}
.tablewrap{overflow:auto;margin-top:16px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:right;padding:8px;border-bottom:1px solid var(--line);white-space:nowrap}th{background:#f5f7fa;position:sticky;top:0}th:first-child,td:first-child{text-align:left}
canvas{width:100%;height:230px;margin-top:12px}.hide{display:none}.wide{grid-column:1/-1}
@media(max-width:820px){.grid{grid-template-columns:1fr}.metrics{grid-template-columns:repeat(2,1fr)}main{margin:14px auto;padding:0 12px}}
</style></head>
<body><main>
<header class="hero"><h1>Microsatellite Bottleneck Analyzer</h1>
<p>Heterozygosity statistics, one-sided Wilcoxon tests, and an SMM/TPM equilibrium comparison.</p></header>
<div class="grid">
<section class="card"><h2>Input and model</h2>
<form id="form">
<div class="field"><label for="file">Genotype CSV</label><input id="file" name="file" type="file" accept=".csv,text/csv" required>
<div class="hint">Sample ID first, then paired allele columns per locus (integer repeat sizes).</div></div>
<div class="field"><label for="model">Mutation model</label><select id="model" name="model"><option>TPM</option><option>SMM</option></select></div>
<div class="field"><label for="ne">Effective population size (Ne)</label><input id="ne" name="ne" type="number" min="2" value="50"></div>
<div class="field"><label for="mu">Mutation rate per allele per generation</label><input id="mu" name="mu" type="number" min="0.00000001" step="any" value="0.001"></div>
<div class="field"><label for="reps">Equilibrium replicates per locus</label><input id="reps" name="replicates" type="number" min="50" max="10000" value="500"></div>
<div class="field"><label for="pstep">TPM one-step mutation probability</label><input id="pstep" name="p_step" type="number" min="0" max="1" step="0.01" value="0.9"></div>
<div class="field"><label for="gprob">TPM geometric probability</label><input id="gprob" name="geometric_p" type="number" min="0.01" max="1" step="0.01" value="0.5"></div>
<div class="field"><label for="seed">Random seed</label><input id="seed" name="seed" type="number" value="42"></div>
<button id="run" type="submit">Analyze CSV</button><div id="error" class="error"></div></form>
<p class="hint">TPM multi-step sizes are at least 2 repeats with a geometric tail. The equilibrium comparison is model-dependent.</p>
</section>
<section class="card"><h2>Results</h2><div id="empty" class="subtle">Choose a CSV and run the analysis.</div>
<div id="results" class="hide">
<div class="metrics" id="metrics"></div><div class="test" id="test"></div>
<div><canvas id="chart" width="900" height="260"></canvas></div>
<button class="secondary" id="dlcsv" type="button">Download locus results</button>
<button class="secondary" id="dljson" type="button">Download full JSON</button>
<p class="hint" id="modelnote"></p>
</div></section>
<section class="card wide"><h2>Per-locus results</h2><div class="tablewrap"><table id="table"></table></div></section>
</div></main>
<script>
let latest=null;
const form=document.getElementById('form');
const fmt=(x,n=4)=>Number(x).toFixed(n);
const esc=s=>String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
form.addEventListener('submit',async e=>{
 e.preventDefault();document.getElementById('error').textContent='';
 const btn=document.getElementById('run');btn.disabled=true;btn.textContent='Simulating…';
 try{
  const fd=new FormData(form);const r=await fetch('/analyze',{method:'POST',body:fd});
  const data=await r.json();if(!r.ok)throw new Error(data.error||'Analysis failed');
  latest=data;render(data);
 }catch(err){document.getElementById('error').textContent=err.message}
 finally{btn.disabled=false;btn.textContent='Analyze CSV'}
});
function render(d){
 document.getElementById('empty').classList.add('hide');
 document.getElementById('results').classList.remove('hide');
 const w=d.wilcoxon, q=d.wilcoxon_he_vs_eq;
 document.getElementById('metrics').innerHTML=
  `<div class="metric"><b>${d.sample_count}</b><span>Individuals</span></div>
   <div class="metric"><b>${d.locus_count}</b><span>Loci</span></div>
   <div class="metric"><b>${Number(w.p_value).toExponential(3)}</b><span>Wilcoxon p (Ho &gt; He)</span></div>
   <div class="metric"><b>${Number(q.p_value).toExponential(3)}</b><span>Wilcoxon p (He &gt; Heq)</span></div>`;
 document.getElementById('test').textContent=
  `Wilcoxon, ${w.alternative}: W+ = ${fmt(w.W_positive,2)}, p = ${Number(w.p_value).toExponential(4)} `+
  `(${w.method}; ${w.n_nonzero_loci} nonzero loci; median diff ${fmt(w.median_difference)}).\n`+
  `Wilcoxon, ${q.alternative}: W+ = ${fmt(q.W_positive,2)}, p = ${Number(q.p_value).toExponential(4)} `+
  `(${q.method}; ${q.n_nonzero_loci} nonzero loci; median diff ${fmt(q.median_difference)}).`;
 document.getElementById('modelnote').textContent=
  `Equilibrium model: ${d.model.name}, Ne=${d.model.Ne}, μ=${d.model.mutation_rate}, `+
  `${d.model.replicates} replicates per locus. Small equilibrium upper-tail p-values mean `+
  `observed He is high relative to simulated He at similar allele count.`;
 const cols=['locus','n','allele_count','Ho','He','equilibrium_mean_He','equilibrium_95_low',
 'equilibrium_95_high','equilibrium_upper_tail_p','conditioning'];
 const labels=['Locus','N','Alleles','Ho','He','Eq. mean He','Eq. 2.5%','Eq. 97.5%',
 'Eq. upper-tail p','Conditioning'];
 let t='<thead><tr>'+labels.map(x=>`<th>${x}</th>`).join('')+'</tr></thead><tbody>';
 d.per_locus.forEach(x=>{t+='<tr>'+cols.map(k=>{
  let v=x[k];if(['Ho','He','equilibrium_mean_He','equilibrium_95_low','equilibrium_95_high'].includes(k))v=fmt(v);
  if(k==='equilibrium_upper_tail_p')v=Number(v).toExponential(2);
  return `<td>${esc(v)}</td>`}).join('')+'</tr>'});
 document.getElementById('table').innerHTML=t+'</tbody>';
 draw(d.per_locus);
}
function draw(rows){
 const c=document.getElementById('chart'),ctx=c.getContext('2d'),W=c.width,H=c.height;
 ctx.clearRect(0,0,W,H);const left=48,right=15,top=15,bottom=35;
 const max=Math.max(.1,...rows.flatMap(x=>[x.Ho,x.He,x.equilibrium_mean_He,x.equilibrium_95_high]));
 const y=v=>H-bottom-(v/max)*(H-top-bottom), step=(W-left-right)/Math.max(1,rows.length);
 ctx.strokeStyle='#dce3eb';ctx.fillStyle='#617083';ctx.font='12px system-ui';
 for(let i=0;i<=4;i++){let v=max*i/4, yy=y(v);ctx.beginPath();ctx.moveTo(left,yy);ctx.lineTo(W-right,yy);ctx.stroke();ctx.fillText(v.toFixed(2),4,yy+4)}
 rows.forEach((x,i)=>{let cx=left+step*(i+.5);
  ctx.strokeStyle='#9aa9bb';ctx.beginPath();ctx.moveTo(cx,y(x.equilibrium_95_low));ctx.lineTo(cx,y(x.equilibrium_95_high));ctx.stroke();
  [[x.Ho,'#2357d6',-5],[x.He,'#d36929',0],[x.equilibrium_mean_He,'#16835d',5]].forEach(([v,color,dx])=>{
   ctx.fillStyle=color;ctx.beginPath();ctx.arc(cx+dx,y(v),3.7,0,Math.PI*2);ctx.fill()});
  ctx.fillStyle='#617083';ctx.fillText(x.locus,left+step*i+2,H-10);
 });
 ctx.fillStyle='#2357d6';ctx.fillText('Ho',left,12);ctx.fillStyle='#d36929';ctx.fillText('He',left+32,12);
 ctx.fillStyle='#16835d';ctx.fillText('Equilibrium mean He (bars: 95% interval)',left+65,12);
}
function download(name,text,type){const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type}));
 a.download=name;a.click();URL.revokeObjectURL(a.href)}
document.getElementById('dlcsv').onclick=()=>{if(latest)download('bottleneck_per_locus.csv',toCsv(latest),'text/csv')};
document.getElementById('dljson').onclick=()=>{if(latest)download('bottleneck_summary.json',JSON.stringify(latest,null,2),'application/json')};
function toCsv(d){
 const keys=['locus','n','allele_count','Ho','He','difference','equilibrium_mean_He',
 'equilibrium_median_He','equilibrium_95_low','equilibrium_95_high',
 'equilibrium_upper_tail_p','equilibrium_n','conditioning','allele_frequencies'];
 const q=v=>'"'+String(typeof v==='object'?JSON.stringify(v):v).replaceAll('"','""')+'"';
 return keys.join(',')+'\n'+d.per_locus.map(x=>keys.map(k=>q(x[k])).join(',')).join('\n')+'\n';
}
</script></body></html>"""


def make_handler(defaults: dict):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, fmt, *args):
            print(f"[web] {self.address_string()} - {fmt % args}")

        def send_json(self, status: int, obj: dict):
            data = json.dumps(obj).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            if urlparse(self.path).path != "/":
                self.send_error(404)
                return
            data = WEB_PAGE.encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def do_POST(self):
            if urlparse(self.path).path != "/analyze":
                self.send_error(404)
                return
            try:
                length = int(self.headers.get("Content-Length", "0"))
                if length <= 0 or length > 20_000_000:
                    raise ValueError("Upload must be between 1 byte and 20 MB.")
                body = self.rfile.read(length)
                mime = (
                    f"Content-Type: {self.headers.get('Content-Type')}\r\n"
                    "MIME-Version: 1.0\r\n\r\n"
                ).encode("ascii") + body
                message = BytesParser(policy=policy.default).parsebytes(mime)
                fields = {}
                csv_text = None
                for part in message.iter_parts():
                    name = part.get_param("name", header="content-disposition")
                    payload = part.get_payload(decode=True) or b""
                    if name == "file":
                        csv_text = payload.decode("utf-8-sig")
                    elif name:
                        fields[name] = payload.decode("utf-8", errors="replace")
                if not csv_text:
                    raise ValueError("Upload a genotype CSV file.")
                options = {
                    "model": fields.get("model", defaults["model"]).upper(),
                    "ne": int(fields.get("ne", defaults["ne"])),
                    "mu": float(fields.get("mu", defaults["mu"])),
                    "replicates": int(fields.get("replicates", defaults["replicates"])),
                    "seed": int(fields.get("seed", defaults["seed"])),
                    "p_step": float(fields.get("p_step", defaults["p_step"])),
                    "geometric_p": float(fields.get("geometric_p", defaults["geometric_p"])),
                }
                validate_options(options)
                self.send_json(200, analyze(csv_text, **options))
            except Exception as exc:
                self.send_json(400, {"error": str(exc)})

    return Handler


def validate_options(options: dict):
    if options["model"] not in {"SMM", "TPM"}:
        raise ValueError("Mutation model must be SMM or TPM.")
    if options["ne"] < 2:
        raise ValueError("Ne must be at least 2.")
    if options["mu"] <= 0 or not math.isfinite(options["mu"]):
        raise ValueError("Mutation rate must be a finite number greater than zero.")
    if not 50 <= options["replicates"] <= 10000:
        raise ValueError("Replicates must be between 50 and 10,000.")
    if not 0 <= options["p_step"] <= 1:
        raise ValueError("TPM one-step probability must be between 0 and 1.")
    if not 0 < options["geometric_p"] <= 1:
        raise ValueError("TPM geometric probability must be in (0, 1].")


def run_web(args):
    options = {
        "model": args.model,
        "ne": args.ne,
        "mu": args.mu,
        "replicates": args.replicates,
        "seed": args.seed,
        "p_step": args.p_step,
        "geometric_p": args.geometric_p,
    }
    validate_options(options)
    server = ThreadingHTTPServer((args.host, args.port), make_handler(options))
    actual_port = server.server_address[1]
    print(f"Analyzer running at http://localhost:{actual_port}/")
    print("Press Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped.")
    finally:
        server.server_close()


def build_parser():
    parser = argparse.ArgumentParser(
        description="Analyze microsatellite heterozygosity and equilibrium baselines."
    )
    parser.add_argument("csv", nargs="?", help="Input genotype CSV; omit to start web UI.")
    parser.add_argument("--web", action="store_true", help="Start the local web interface.")
    parser.add_argument("--host", default="127.0.0.1", help="Web server host (default: localhost).")
    parser.add_argument("--port", type=int, default=8000, help="Web server port.")
    parser.add_argument("--model", choices=["TPM", "SMM"], default="TPM")
    parser.add_argument("--ne", type=int, default=50, help="Diploid effective population size.")
    parser.add_argument("--mu", type=float, default=0.001, help="Mutation rate per allele per generation.")
    parser.add_argument("--replicates", type=int, default=500, help="Equilibrium simulations per locus.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--p-step", type=float, default=0.9, help="TPM one-step mutation probability.")
    parser.add_argument("--geometric-p", type=float, default=0.5, help="TPM geometric-tail probability.")
    parser.add_argument("--out-prefix", help="Write PREFIX_per_locus.csv and PREFIX_summary.json.")
    return parser


def main():
    args = build_parser().parse_args()
    if args.web or not args.csv:
        run_web(args)
        return

    options = {
        "model": args.model,
        "ne": args.ne,
        "mu": args.mu,
        "replicates": args.replicates,
        "seed": args.seed,
        "p_step": args.p_step,
        "geometric_p": args.geometric_p,
    }
    validate_options(options)
    try:
        text = Path(args.csv).read_text(encoding="utf-8-sig")
        result = analyze(text, **options)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(2)

    for key in ("wilcoxon", "wilcoxon_he_vs_eq"):
        w = result[key]
        print(f"Wilcoxon alternative: {w['alternative']}  "
              f"W+ = {w['W_positive']:.3f} | p = {w['p_value']:.6g} ({w['method']}) | "
              f"median diff = {w['median_difference']:.6f}")
    print(f"Samples: {result['sample_count']} | Loci: {result['locus_count']}")
    print(f"Equilibrium model: {args.model}, Ne={args.ne}, mu={args.mu}, reps={args.replicates}")
    print("\nLocus\tN\tK\tHo\tHe\tEq.mean.He\tEq.upper-tail.p")
    for x in result["per_locus"]:
        print(
            f"{x['locus']}\t{x['n']}\t{x['allele_count']}\t{x['Ho']:.4f}\t"
            f"{x['He']:.4f}\t{x['equilibrium_mean_He']:.4f}\t"
            f"{x['equilibrium_upper_tail_p']:.4g}"
        )

    if args.out_prefix:
        prefix = Path(args.out_prefix)
        csv_path = prefix.with_name(prefix.name + "_per_locus.csv")
        json_path = prefix.with_name(prefix.name + "_summary.json")
        csv_path.write_text(csv_report(result), encoding="utf-8", newline="")
        json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(f"\nWrote {csv_path} and {json_path}")


if __name__ == "__main__":
    main()
