#!/usr/bin/env python3
"""Test 5b: MrBayes 3.2.7 full-DNA FBD (make_mb_dna.py; IGR clock) against RevBayes full-DNA FBD (rb_e2, 6 chains, 23 genes) and
ffbd given the 23-gene phylogram (runs/ctrl_SET_g23_r*.log), on the extant node ages of A06 and A12. MrBayes tree lengths are
rate x time, so each node age is (tree height - depth) / clockrate of that sample. The MrBayes prior-only run (same priors, no data)
shows how much of each age comes from the DNA. First 10% of every chain dropped.

  compare_mbdna.py SET CLADES.Rev [BASE] [--tsv OUT]   BASE = directory holding rb_e2/, mbdna/ and runs/ (default .)
"""
import csv, glob, re, statistics as st, sys
TSV = None
if '--tsv' in sys.argv:
    i = sys.argv.index('--tsv'); TSV = sys.argv[i + 1]; del sys.argv[i:i + 2]
S, CL = sys.argv[1], sys.argv[2]
B = sys.argv[3] if len(sys.argv) > 3 else '.'
ent = open(CL).read()
clades = {}
for m in re.finditer(r'^clade_(n\d+) = clade\((.*)\)$', ent, re.M):
    s = set()
    for p in m.group(2).split(','):
        p = p.strip(); s |= clades[p[6:]] if p.startswith('clade_') else {p.strip('"')}
    clades[m.group(1)] = frozenset(s)
BURN = 0.1


def logs(files):
    out = {k: [] for k in clades}
    for f in files:
        r = list(csv.DictReader(open(f), delimiter='\t')); r = r[int(len(r) * BURN):]
        for k in clades: out[k] += [float(x[f'age_{k}']) for x in r]
    return out


def parse(s, i=0):
    ch = []
    if s[i] == '(':
        i += 1
        while True:
            c, i = parse(s, i); ch.append(c)
            if s[i] == ',': i += 1; continue
            i += 1; break
    m = re.match(r'[^,():;]*', s[i:]); name = m.group(0); i += len(name); ln = 0.0
    if i < len(s) and s[i] == ':':
        m = re.match(r'[-0-9.eE+]+', s[i + 1:]); ln = float(m.group(0)); i += 1 + len(m.group(0))
    return (ch, name, ln), i


def mb(files):
    out = {k: [] for k in clades}
    for f in files:
        txt = open(f).read()
        tr = dict(re.findall(r'^\s+(\d+) ([^,;\s]+)[,;]', txt.split('translate', 1)[1].split(';', 1)[0] + ';', re.M))
        trees = re.findall(r'^\s*tree gen\.\d+ .*?\[&clockrate=([-0-9.eE+]+)\]\s*(.*);', txt, re.M)
        trees = trees[int(len(trees) * BURN):]
        for cr, t in trees:
            cr = float(cr); t = re.sub(r'\[[^\]]*\]', '', t)
            root, _ = parse(t); tips, nodes = {}, []
            def rec(n, d):
                ch, name, ln = n; d2 = d + ln
                if not ch:
                    tips[tr[name]] = d2; return {tr[name]}
                s = set().union(*[rec(c, d2) for c in ch]); nodes.append((frozenset(s), d2)); return s
            rec(root, 0.0)
            H = max(d for nm, d in tips.items() if not nm.startswith('FOSSIL'))
            byext = {}
            for s, d in nodes:
                e = frozenset(x for x in s if not x.startswith('FOSSIL'))
                if e and (e not in byext or d > byext[e]): byext[e] = d
            for k, c in clades.items(): out[k].append((H - byext[c]) / cr)
    return out


def q(v, a):
    v = sorted(v); return v[min(len(v) - 1, int(a * len(v)))]


rb = sorted(glob.glob(f'{B}/rb_e2/e2_{S}_dna_r*.log'))
A = {'RevBayes DNA (6 chains)': logs(rb), 'RevBayes DNA chains 1-3': logs(rb[:3]), 'RevBayes DNA chains 4-6': logs(rb[3:]),
     'MrBayes DNA': mb([f'{B}/mbdna/mbdna_{S}_r1.t', f'{B}/mbdna/mbdna_{S}_r2.t']),
     'MrBayes DNA chain 1': mb([f'{B}/mbdna/mbdna_{S}_r1.t']), 'MrBayes DNA chain 2': mb([f'{B}/mbdna/mbdna_{S}_r2.t']),
     'MrBayes prior only': mb([f'{B}/mbdna/mbprior_{S}_r1.t', f'{B}/mbdna/mbprior_{S}_r2.t']),
     'ffbd 23-gene phylogram': logs([f'{B}/runs/ctrl_{S}_g23_r1.log', f'{B}/runs/ctrl_{S}_g23_r2.log'])}
M = {a: {k: st.median(v) for k, v in d.items()} for a, d in A.items()}
I = {a: {k: (q(v, .025), q(v, .975)) for k, v in d.items()} for a, d in A.items()}
print(f'{S}: {len(clades)} extant nodes; samples', {a: len(next(iter(d.values()))) for a, d in A.items()})
print(f'{"comparison":58} {"r":>7} {"median %":>9} {"max %":>7} {"mean shift Myr":>15} {"inside 95% of ref":>18}')
def cmp(x, y, lab):
    ks = list(clades); a = [M[x][k] for k in ks]; b = [M[y][k] for k in ks]
    ma, mb_ = st.fmean(a), st.fmean(b)
    r = sum((p - ma) * (o - mb_) for p, o in zip(a, b)) / (sum((p - ma) ** 2 for p in a) * sum((o - mb_) ** 2 for o in b)) ** .5
    d = [100 * abs(p - o) / o for p, o in zip(a, b)]
    ins = sum(I[y][k][0] <= M[x][k] <= I[y][k][1] for k in ks)
    print(f'{lab:58} {r:7.4f} {st.median(d):9.2f} {max(d):7.2f} {st.fmean(p - o for p, o in zip(a, b)):+15.2f} {ins:>12}/{len(ks)}')
cmp('MrBayes DNA chain 1', 'MrBayes DNA chain 2', 'noise: MrBayes DNA chain 1 vs chain 2')
cmp('RevBayes DNA chains 1-3', 'RevBayes DNA chains 4-6', 'noise: RevBayes DNA chains 1-3 vs 4-6')
cmp('MrBayes DNA', 'RevBayes DNA (6 chains)', 'MrBayes DNA vs RevBayes DNA')
cmp('ffbd 23-gene phylogram', 'MrBayes DNA', 'ffbd 23-gene phylogram vs MrBayes DNA')
cmp('ffbd 23-gene phylogram', 'RevBayes DNA (6 chains)', 'ffbd 23-gene phylogram vs RevBayes DNA')
cmp('MrBayes prior only', 'MrBayes DNA', 'MrBayes prior only vs MrBayes DNA')
root = max(clades, key=lambda k: len(clades[k]))
print('root (Ma):', {a: round(M[a][root], 1) for a in ('RevBayes DNA (6 chains)', 'MrBayes DNA', 'ffbd 23-gene phylogram', 'MrBayes prior only')})
if TSV:
    with open(TSV, 'w') as fh:
        fh.write('node\tntips\tprogram\tmedian\tlo95\thi95\n')
        for a in M:
            for k in clades: fh.write(f'{k}\t{len(clades[k])}\t{a}\t{M[a][k]:.4f}\t{I[a][k][0]:.4f}\t{I[a][k][1]:.4f}\n')
