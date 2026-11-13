#!/usr/bin/env python3
"""Score the few-genes against many-genes simulation (sim_genes.sbatch) against the true node ages.
Arms: exact (true branch lengths, from the simulation study sim1), few (23 genes, bootstrap SEs), few_rule (23 genes, SE rule),
many_rule (1,000 genes, SE rule). For each extant node: 95% interval coverage, error of the posterior median relative to the true
age, absolute error in Myr, and interval width. Pooled over replicates, overall and by true age.  First 10% of samples dropped.

  sim_genes_summary.py SIM2DIR SIM1DIR [OUT.tsv]
"""
import csv, glob, math, os, statistics as st, sys

D2, D1 = sys.argv[1], sys.argv[2]
OUT = sys.argv[3] if len(sys.argv) > 3 else None
ARMS = [('exact', D1, 'true branch lengths'), ('few', D2, '23 genes, bootstrap SE'), ('few_rule', D2, '23 genes, SE rule'),
        ('many_rule', D2, '1,000 genes, SE rule')]
BINS = [(0, 5), (5, 10), (10, 20), (20, 40), (40, 60), (60, 101)]


def q(x, a):
    i = a * (len(x) - 1); lo = int(i)
    return x[lo] + (x[min(lo + 1, len(x) - 1)] - x[lo]) * (i - lo)


rows, reps = [], []
for rd in sorted(glob.glob(f'{D2}/r*'), key=lambda s: int(s.rsplit('r', 1)[1])):
    r = os.path.basename(rd)
    logs = {a: f'{base}/{r}/{a}.log' for a, base, _ in ARMS}
    if not all(os.path.exists(f) and sum(1 for _ in open(f)) > 1000 for f in logs.values()):
        continue
    reps.append(r)
    truth = {l.split('\t')[0]: float(l.split('\t')[1]) for l in open(f'{rd}/truth.tsv').read().splitlines()[1:]}
    for a, f in logs.items():
        lg = list(csv.DictReader(open(f), delimiter='\t'))
        lg = lg[len(lg) // 10:]
        for p, t in truth.items():
            if not p.startswith('age_n') or p not in lg[0]:
                continue
            x = sorted(float(z[p]) for z in lg)
            lo, med, hi = q(x, .025), q(x, .5), q(x, .975)
            rows.append({'rep': r, 'arm': a, 'node': p[4:], 'true': t, 'median': med, 'lo95': lo, 'hi95': hi,
                         'in95': int(lo <= t <= hi), 'relerr': (med - t) / t, 'abserr': abs(med - t), 'width': hi - lo})
if OUT:
    with open(OUT, 'w') as fh:
        w = csv.DictWriter(fh, rows[0].keys(), delimiter='\t'); w.writeheader(); w.writerows(rows)
print(f'replicates with all four arms: {len(reps)}; nodes per arm: {sum(1 for x in rows if x["arm"] == "exact")}')
print(f'{"arm":10} {"nodes":>6} {"cov95":>6} {"median rel err":>14} {"median |err| Myr":>16} {"mean |err| Myr":>14} {"median width Myr":>16}')
for a, _, desc in ARMS:
    s = [x for x in rows if x['arm'] == a]
    print(f'{a:10} {len(s):6d} {st.fmean(x["in95"] for x in s):6.3f} {100 * st.median(x["relerr"] for x in s):13.1f}% '
          f'{st.median(x["abserr"] for x in s):16.2f} {st.fmean(x["abserr"] for x in s):14.2f} {st.median(x["width"] for x in s):16.2f}   {desc}')
print('\nby true node age (Ma): coverage / median relative error / median absolute error (Myr)')
print(f'{"arm":10} ' + ' '.join(f'{f"{lo}-{hi if hi < 100 else 100}":>20}' for lo, hi in BINS))
for a, _, _ in ARMS:
    out = []
    for lo, hi in BINS:
        s = [x for x in rows if x['arm'] == a and lo <= x['true'] < hi]
        out.append(f'{st.fmean(x["in95"] for x in s):.2f}/{100 * st.median(x["relerr"] for x in s):+.0f}%/{st.median(x["abserr"] for x in s):.2f}' if s else '-')
    print(f'{a:10} ' + ' '.join(f'{o:>20}' for o in out))
# paired: per node, is many_rule closer to the truth than few?
pair = {}
for x in rows:
    pair.setdefault((x['rep'], x['node']), {})[x['arm']] = x
both = [v for v in pair.values() if 'few' in v and 'many_rule' in v]
closer = sum(abs(v['many_rule']['median'] - v['many_rule']['true']) < abs(v['few']['median'] - v['few']['true']) for v in both)
print(f'\npaired nodes: many_rule closer to the truth than few on {closer} of {len(both)} ({100 * closer / len(both):.1f}%)')
