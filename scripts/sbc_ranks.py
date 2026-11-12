#!/usr/bin/env python3
"""Rank test for the realistic-size SBC (sim_study.py sbc / sbc_study.sbatch). For each replicate: first 10% of samples dropped, the rest
thinned to 99 draws, rank of the true value among them (0..99, ties at random). Correct sampler and model -> ranks uniform. Reports a
chi-square test on 20 bins, the mean rank (0.5 expected) and the share of ranks in the outer 5% on each side (0.05 expected; too many =
intervals too narrow, too few = too wide; a lopsided split = bias). Also the median ESS of each parameter (draws must be near independent).

  sbc_ranks.py SBCDIR
"""
import csv, glob, math, os, random, statistics, sys

D = sys.argv[1]
K = 99
random.seed(1)
PAR = {'speciation_rate': 'speciation_rate', 'extinction_rate': 'extinction_rate', 'psi': 'psi', 'branch_rates_mean': 'branch_rates_mean',
       'root': 'root_age_tree', 'num_samp_anc': 'num_samp_anc', 'sbc_side1': None, 'sbc_side2': None, 'sbc_fossil': None}


def chisq_sf(x, df):
    a, x2 = df / 2.0, x / 2.0
    if x2 < a + 1:
        s, t, k = 0.0, 1.0 / a, 0
        while abs(t) > 1e-15 * abs(s) or k == 0:
            s += t
            k += 1
            t *= x2 / (a + k)
        return 1 - s * math.exp(-x2 + a * math.log(x2) - math.lgamma(a))
    b, c, d = x2 + 1 - a, 1e300, 1 / (x2 + 1 - a)
    h = d
    for i in range(1, 300):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1 / d if d else 1e300
        c = b + an / c
        c = c if c else 1e-300
        h *= d * c
    return math.exp(-x2 + a * math.log(x2) - math.lgamma(a)) * h


def ess(v):
    n, m = len(v), statistics.fmean(v)
    var = statistics.pvariance(v)
    if var == 0:
        return float(n)
    s = 0.0
    for k in range(1, n // 2):
        r = sum((v[i] - m) * (v[i + k] - m) for i in range(n - k)) / ((n - k) * var)
        if r < 0.05:
            break
        s += r
    return n / (1 + 2 * s)


ranks = {p: [] for p in PAR}
esses = {p: [] for p in PAR}
sizes = []
for rd in glob.glob(f'{D}/r*'):
    f = f'{rd}/exact.log'
    if not os.path.exists(f) or not os.path.exists(f'{rd}/truth.tsv'):
        continue
    truth = dict(l.split('\t') for l in open(f'{rd}/truth.tsv').read().splitlines()[1:])
    lg = list(csv.DictReader(open(f), delimiter='\t'))
    if len(lg) < 200:
        continue
    lg = lg[len(lg) // 10:]
    step = len(lg) / K
    thin = [lg[int(i * step)] for i in range(K)]
    sizes.append((int(truth['n_extant']), int(truth['n_fossils'])))
    for p, col in PAR.items():
        if p not in truth:
            continue
        if col is None:
            col = truth[p]         # log column: age_n12 (truth key age_n12) or age_FOSSIL_x (truth key FOSSIL_x)
            t = float(truth[col] if col in truth else truth[col[4:]])
        else:
            t = float(truth[p])
        if col not in thin[0]:
            continue
        v = [float(r[col]) for r in thin]
        less = sum(x < t for x in v)
        eq = sum(x == t for x in v)
        ranks[p].append(less + random.randint(0, eq))
        esses[p].append(ess([float(r[col]) for r in lg]))
print(f'replicates {len(sizes)}; extant tips median {statistics.median(s[0] for s in sizes) if sizes else 0} (range {min(s[0] for s in sizes)}-{max(s[0] for s in sizes)}), '
      f'fossils median {statistics.median(s[1] for s in sizes)} (range {min(s[1] for s in sizes)}-{max(s[1] for s in sizes)})')
print(f"{'parameter':18} {'n':>4} {'chi2 p':>8} {'mean rank':>9} {'low 5%':>7} {'high 5%':>7} {'median ESS':>10}")
for p in PAR:
    r = ranks[p]
    if not r:
        continue
    bins = [0] * 20
    for x in r:
        bins[min(19, x * 20 // (K + 1))] += 1
    e = len(r) / 20
    chi = sum((b - e) ** 2 / e for b in bins)
    print(f'{p:18} {len(r):4d} {chisq_sf(chi, 19):8.4f} {statistics.fmean(r) / K:9.3f} {sum(x < 5 for x in r) / len(r):7.3f} '
          f'{sum(x > K - 5 for x in r) / len(r):7.3f} {statistics.median(esses[p]):10.0f}')
