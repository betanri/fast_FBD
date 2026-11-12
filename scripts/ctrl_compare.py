#!/usr/bin/env python3
"""Score the controlled data test (Step 10; ctrl_arms.sbatch). Reference: full-DNA FBD in RevBayes on the 23-gene alignment
(rb_e2/e2_SET_dna_r*.log, collapse-expand move removed). Each Fast FBD arm (ffbd, runs/ctrl_SET_ARM_r{1,2}.log) is compared
with it on extant node ages (medians, matched by clade name nK, identical in every arm), the root age and the mean clock rate.
The noise floor is full DNA against itself (chains 1-3 against chains 4-6). Then the branch-length inputs are compared edge
by edge (log-log r, ratio to the 23-gene lengths, and the same for short edges).

  ctrl_compare.py SET [BASE] [--tsv OUT.tsv]   BASE = the run directory (default: current directory); OUT.tsv: per node and arm, median and 95% interval
"""
import csv, glob, math, os, statistics as st, sys

TSV = None
if '--tsv' in sys.argv:
    i = sys.argv.index('--tsv'); TSV = sys.argv[i + 1]; del sys.argv[i:i + 2]
S = sys.argv[1]
B = sys.argv[2] if len(sys.argv) > 2 else '.'
DUMP = []
ARMS = [('orig23', 'runs/e2_{S}_r{r}.log', 'original Fast FBD input of test 4a (23 genes)'),
        ('g23', 'runs/ctrl_{S}_g23_r{r}.log', '23-gene phylogram, bootstrap SE'),
        ('allgene', 'runs/ctrl_{S}_allgene_r{r}.log', 'all-gene phylogram, bootstrap SE'),
        ('allgene_rule', 'runs/ctrl_{S}_allgene_rule_r{r}.log', 'all-gene phylogram, SE rule sqrt(b/20000)'),
        ('large_tree_a', 'runs/ctrl_{S}_julian_r{r}.log', 'genome-scale phylogram of a larger tree, pruned to the clade, SE rule'),
        ('large_tree_b', 'runs/ctrl_{S}_c01_r{r}.log', 'a second larger-tree phylogram, pruned, SE rule')]


def load(files):
    out = {}
    for f in files:
        rows = list(csv.DictReader(open(f), delimiter='\t'))
        rows = rows[len(rows) // 10:]
        for k in rows[0]:
            if k.startswith('age_n') or k in ('root_age_tree', 'branch_rates_mean'):
                out.setdefault(k, []).extend(float(x[k]) for x in rows if x[k] != '')
    return out


def q(v, a):
    v = sorted(v); i = a * (len(v) - 1); lo = int(i)
    return v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (i - lo)


def pearson(x, y):
    mx, my = st.fmean(x), st.fmean(y)
    sx = math.sqrt(sum((a - mx) ** 2 for a in x)); sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def compare(ref, arm, label):
    nodes = sorted((k for k in ref if k.startswith('age_n') and k in arm), key=lambda k: int(k[5:]))
    DUMP.extend((label, k[4:], q(ref[k], .5), q(ref[k], .025), q(ref[k], .975), q(arm[k], .5), q(arm[k], .025), q(arm[k], .975)) for k in nodes)
    mr = [q(ref[k], 0.5) for k in nodes]; ma = [q(arm[k], 0.5) for k in nodes]
    pct = [100 * abs(a - r) / r for a, r in zip(ma, mr)]
    young = [p for p, r in zip(pct, mr) if r < 20]
    inside = sum(q(ref[k], 0.025) <= a <= q(ref[k], 0.975) for k, a in zip(nodes, ma))
    wid = [(q(arm[k], 0.975) - q(arm[k], 0.025)) / (q(ref[k], 0.975) - q(ref[k], 0.025)) for k in nodes]
    print(f'{label:14} {len(nodes):5d} {pearson(mr, ma):6.4f} {st.median(pct):7.2f} {max(pct):7.2f} '
          f'{(st.median(young) if young else float("nan")):9.2f} {st.fmean(a - r for a, r in zip(ma, mr)):+9.2f} '
          f'{inside:3d}/{len(nodes):<3d} {st.median(wid):6.2f} {q(arm["root_age_tree"], .5):7.2f} '
          f'{q(arm["branch_rates_mean"], .5) * 1000:8.3f}')


chains = sorted(glob.glob(f'{B}/rb_e2/e2_{S}_dna_r*.log'))
ref = load(chains)
print(f'{S}: full-DNA reference, {len(chains)} RevBayes chains, {len(ref["root_age_tree"])} samples after 10% burn-in each; '
      f'root {q(ref["root_age_tree"], .5):.2f} Ma, mean clock rate x1000 {q(ref["branch_rates_mean"], .5) * 1000:.3f}')
print(f'{"arm":14} {"nodes":>5} {"r":>6} {"med %":>7} {"max %":>7} {"<20Ma %":>9} {"mean dMyr":>9} {"in95":>7} {"width":>6} '
      f'{"root":>7} {"rate e3":>8}')
h = len(chains) // 2
compare(load(chains[:h]), load(chains[h:]), 'fullDNA noise')
for name, pat, desc in ARMS:
    fs = [f'{B}/{pat.format(S=S, r=r)}' for r in (1, 2)]
    if all(os.path.exists(f) for f in fs):
        compare(ref, load(fs), name)
print('  (r, median and largest % difference of node-age medians against full DNA; <20Ma: median % on nodes under 20 Ma;'
      ' in95: arm median inside the full-DNA 95% interval; width: interval width ratio, median over nodes)')

print(f'\n{S}: branch-length inputs, edge by edge against the 23-gene ML lengths (g23)')
tabs = {'g23': f'{B}/ctrl/{S}/g23/terms/branch_table.tsv', 'allgene': f'{B}/ctrl/{S}/allgene/terms/branch_table.tsv',
        'large_tree_a': f'{B}/ctrl/{S}/julian/tab/branch_table.tsv', 'large_tree_b': f'{B}/ctrl/{S}/c01/tab/branch_table.tsv'}
bl = {k: {int(r['term']): (float(r['b_ml']), float(r['se_used'])) for r in csv.DictReader(open(f), delimiter='\t')}
      for k, f in tabs.items() if os.path.exists(f)}
g = bl['g23']
short = [j for j in g if g[j][0] < 0.01]
print(f'{"input":9} {"tree length":>11} {"log r":>6} {"ratio med":>9} {"short log r":>11} {"short ratio":>11} {"SE/b med":>8}   '
      f'(short = {len(short)} edges under 0.01 subs/site in g23)')
for k, v in bl.items():
    j = sorted(v)
    lx = [math.log(g[i][0]) for i in j]; ly = [math.log(v[i][0]) for i in j]
    sx = [math.log(g[i][0]) for i in short]; sy = [math.log(v[i][0]) for i in short]
    print(f'{k:9} {sum(v[i][0] for i in j):11.4f} {pearson(lx, ly):6.3f} {st.median(v[i][0] / g[i][0] for i in j):9.3f} '
          f'{(pearson(sx, sy) if len(short) > 2 else float("nan")):11.3f} '
          f'{(st.median(v[i][0] / g[i][0] for i in short) if short else float("nan")):11.3f} '
          f'{st.median(v[i][1] / v[i][0] for i in j):8.3f}')

if TSV:
    with open(TSV, 'w') as fh:
        fh.write('set\tarm\tnode\tdna_median\tdna_lo95\tdna_hi95\tarm_median\tarm_lo95\tarm_hi95\n')
        for r in DUMP:
            fh.write(f'{S}\t' + '\t'.join(f'{x:.4f}' if isinstance(x, float) else x for x in r) + '\n')
