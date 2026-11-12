#!/usr/bin/env python3
"""Compare the fixed-rate prior runs of ffbd, MrBayes and RevBayes (make_fixed3.py): extant node ages, fossil ages and the number of
sampled ancestors. MrBayes trees (.t) are read directly (sampled ancestor = fossil tip with zero branch length).

  compare_fixed3.py SET LOCALDIR [OUT.tsv]   OUT.tsv: per node and program, median and 95% interval; LOCALDIR holds v3_SET_r{1,2}.log (ffbd), v3rb_SET_r{1,2}.log (RevBayes), v3mb_SET_r{1,2}.t (MrBayes)
"""
import csv, re, statistics, sys

S, D = sys.argv[1], sys.argv[2]
ent = open(f'fast_fbd_test/fish_root/data/{S}_root/extant_nodes.Rev').read()
clades = {m.group(1): frozenset(re.findall(r'"([^"]+)"', m.group(2))) for m in re.finditer(r'^clade_(n\d+) = clade\((.*)\)$', ent, re.M)}
BURN = 0.1


def logcols(f):
    r = list(csv.DictReader(open(f), delimiter='\t'))
    r = r[int(len(r) * BURN):]
    out = {k: [float(x[f'age_{k}']) for x in r] for k in clades}
    out['num_samp_anc'] = [float(x['num_samp_anc']) for x in r]
    return out


def parse(s, i=0):
    """newick -> (children, name, length); returns node, index"""
    ch = []
    if s[i] == '(':
        i += 1
        while True:
            c, i = parse(s, i)
            ch.append(c)
            if s[i] == ',':
                i += 1
                continue
            i += 1
            break
    m = re.match(r'[^,():;\[]*', s[i:])
    name = m.group(0)
    i += len(name)
    ln = 0.0
    if i < len(s) and s[i] == ':':
        m = re.match(r'[-0-9.eE+]+', s[i + 1:])
        ln = float(m.group(0))
        i += 1 + len(m.group(0))
    return (ch, name, ln), i


def mb(f):
    txt = open(f).read()
    tr = dict(re.findall(r'^\s+(\d+) ([^,;\s]+)[,;]', txt.split('translate', 1)[1].split(';', 1)[0] + ';', re.M))
    trees = re.findall(r'^\s*tree gen\.(\d+) = (?:\[[^\]]*\]\s*)*(.*);', txt, re.M)
    trees = trees[int(len(trees) * BURN):]
    out = {k: [] for k in clades}
    out['num_samp_anc'] = []
    for g, t in trees:
        root, _ = parse(t)
        depth, tips, nodes = {}, {}, []
        def rec(n, d):
            ch, name, ln = n
            d2 = d + ln
            if not ch:
                tips[tr[name]] = (d2, ln)
                return {tr[name]}
            s = set().union(*[rec(c, d2) for c in ch])
            nodes.append((frozenset(s), d2))
            return s
        rec(root, 0.0)
        H = max(d for nm, (d, ln) in tips.items() if not nm.startswith('FOSSIL'))
        byext = {}
        for s, d in nodes:
            e = frozenset(x for x in s if not x.startswith('FOSSIL'))
            if e and (e not in byext or d > byext[e]):   # nodes above the MRCA with the same extant set only add fossils: keep the youngest (largest depth)
                byext[e] = d
        for k, c in clades.items():
            out[k].append(H - byext[c])
        out['num_samp_anc'].append(sum(1 for nm, (d, ln) in tips.items() if nm.startswith('FOSSIL') and ln == 0.0))
    return out


runs = {}
for r in (1, 2):
    runs[f'ffbd_r{r}'] = logcols(f'{D}/v3_{S}_r{r}.log')
    runs[f'rb_r{r}'] = logcols(f'{D}/v3rb_{S}_r{r}.log')
    runs[f'mb_r{r}'] = mb(f'{D}/v3mb_{S}_r{r}.t')
pool = {p: {k: runs[f'{p}_r1'][k] + runs[f'{p}_r2'][k] for k in runs[f'{p}_r1']} for p in ('ffbd', 'mb', 'rb')}
med = lambda v: statistics.median(v)
pct = lambda a, b: 100 * abs(a - b) / b
print(f'{S}: samples per chain', {k: len(v['n1']) for k, v in runs.items()})
for a, b in [('ffbd', 'mb'), ('rb', 'mb'), ('ffbd', 'rb')]:
    d = [pct(med(pool[a][k]), med(pool[b][k])) for k in clades]
    print(f'  node-age medians, {a} vs {b}: median {statistics.median(d):.2f}%  max {max(d):.2f}%')
for p in ('ffbd', 'mb', 'rb'):
    d = [pct(med(runs[f'{p}_r1'][k]), med(runs[f'{p}_r2'][k])) for k in clades]
    print(f'  {p} chain 1 vs chain 2: median {statistics.median(d):.2f}%  max {max(d):.2f}%')
print('  mean sampled ancestors: ' + '  '.join(f'{p} {statistics.fmean(pool[p]["num_samp_anc"]):.3f}' for p in ('ffbd', 'mb', 'rb')))
worst = sorted(clades, key=lambda k: -pct(med(pool['rb'][k]), med(pool['mb'][k])))[:5]
print('  largest RevBayes-MrBayes differences (median age ffbd / MrBayes / RevBayes):')
for k in worst:
    print(f'    {k} ({len(clades[k])} tips): {med(pool["ffbd"][k]):.2f} / {med(pool["mb"][k]):.2f} / {med(pool["rb"][k]):.2f}')

if len(sys.argv) > 3:
    qq = lambda v, a: sorted(v)[int(a * (len(v) - 1))]
    with open(sys.argv[3], 'w') as fh:
        fh.write('set\tnode\ttips\tprogram\tmedian\tlo95\thi95\n')
        for k in clades:
            for p in ('ffbd', 'mb', 'rb'):
                v = pool[p][k]
                fh.write(f'{S}\t{k}\t{len(clades[k])}\t{p}\t{med(v):.4f}\t{qq(v, .025):.4f}\t{qq(v, .975):.4f}\n')
        for p in ('ffbd', 'mb', 'rb'):
            for c in range(0, 1 + int(max(pool[p]['num_samp_anc']))):
                fh.write(f'{S}\tnum_samp_anc={c}\t\t{p}\t{sum(1 for x in pool[p]["num_samp_anc"] if x == c) / len(pool[p]["num_samp_anc"]):.5f}\t\t\n')
