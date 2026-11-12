#!/usr/bin/env python3
"""Score the simulation study (sim_study.py / sim_study.sbatch): for each replicate and arm (exact, dna), coverage of the true value
by the 50% and 95% central credible intervals, relative error of the posterior median, and relative interval width, pooled over
extant node ages, fossil ages, root age and the FBD rates.  First 10% of the logged samples are discarded (burn-in is already excluded).

  sim_study_summary.py SIMDIR [OUT.tsv]       (SIMDIR holds r1, r2, ...; each with truth.tsv, exact.log, dna.log)
"""
import csv, glob, math, os, statistics, sys

D = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else None
GROUPS = ('node', 'fossil', 'root', 'speciation_rate', 'extinction_rate', 'psi', 'num_samp_anc', 'branch_rates_mean')


def group(p):
    return 'node' if p.startswith('age_n') else 'fossil' if p.startswith('FOSSIL') else p


def q(x, a):
    i = a * (len(x) - 1)
    lo = int(math.floor(i))
    return x[lo] + (x[min(lo + 1, len(x) - 1)] - x[lo]) * (i - lo)


def ess(v):
    n = len(v)
    m = statistics.fmean(v)
    var = statistics.pvariance(v)
    if var == 0:
        return float(n)
    s, k = 0.0, 1
    while k < n // 2:
        r = sum((v[i] - m) * (v[i + k] - m) for i in range(0, n - k, max(1, (n - k) // 2000))) / (len(range(0, n - k, max(1, (n - k) // 2000))) * var)
        if r < 0.05:
            break
        s += r
        k += max(1, k // 2)
    return n / (1 + 2 * s)


rows = []
for rd in sorted(glob.glob(f'{D}/r*'), key=lambda s: int(s.rsplit('r', 1)[1])):
    if not os.path.exists(f'{rd}/truth.tsv'):
        continue
    truth = dict(l.split('\t') for l in open(f'{rd}/truth.tsv').read().splitlines()[1:])
    for arm in ('exact', 'dna'):
        f = f'{rd}/{arm}.log'
        if not os.path.exists(f):
            continue
        lg = list(csv.DictReader(open(f), delimiter='\t'))
        if len(lg) < 100:
            continue
        lg = lg[len(lg) // 10:]
        for p, v in truth.items():
            col = 'root_age_tree' if p == 'root' else f'age_{p}' if p.startswith('FOSSIL') else p
            if group(p) not in GROUPS or col not in lg[0]:
                continue
            x = sorted(float(r[col]) for r in lg)
            t = float(v)
            med = q(x, 0.5)
            rows.append({'rep': os.path.basename(rd), 'arm': arm, 'param': p, 'group': group(p), 'true': t, 'median': med,
                         'lo95': q(x, 0.025), 'hi95': q(x, 0.975), 'in95': int(q(x, 0.025) <= t <= q(x, 0.975)),
                         'in50': int(q(x, 0.25) <= t <= q(x, 0.75)),
                         'relerr': (med - t) / t if t else float('nan'), 'relwidth': (q(x, 0.975) - q(x, 0.025)) / t if t else float('nan'),
                         'ess': ess([float(r[col]) for r in lg]) if group(p) in ('root', 'speciation_rate', 'num_samp_anc') else float('nan')})

if OUT:
    with open(OUT, 'w') as fh:
        w = csv.DictWriter(fh, rows[0].keys(), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
reps = {(r['rep'], r['arm']) for r in rows}
print('replicates scored:', {a: sum(1 for x in reps if x[1] == a) for a in ('exact', 'dna')})
print(f"{'arm':6} {'group':18} {'n':>6} {'cov95':>6} {'cov50':>6} {'med rel err':>11} {'mean |rel err|':>14} {'rel width':>9} {'min ESS':>8}")
for arm in ('exact', 'dna'):
    for g in GROUPS:
        s = [r for r in rows if r['arm'] == arm and r['group'] == g]
        if not s:
            continue
        e = [r['relerr'] for r in s if not math.isnan(r['relerr'])]
        es = [r['ess'] for r in s if not math.isnan(r['ess'])]
        print(f"{arm:6} {g:18} {len(s):6d} {statistics.fmean(r['in95'] for r in s):6.3f} {statistics.fmean(r['in50'] for r in s):6.3f} "
              f"{statistics.median(e) if e else float('nan'):11.4f} {statistics.fmean(abs(x) for x in e) if e else float('nan'):14.4f} "
              f"{statistics.median(r['relwidth'] for r in s if not math.isnan(r['relwidth'])):9.3f} {min(es) if es else float('nan'):8.0f}")
