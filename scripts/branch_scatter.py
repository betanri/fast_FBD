#!/usr/bin/env python3
"""Branch lengths in time (Myr): full-DNA FBD against Fast FBD, one panel per data set.
Each side is summarized as an extant-only time tree (fossils pruned, conditional median node ages; FBD/TIME_TREES/extant_time_tree.py).
A branch is compared when both trees contain the same child clade with the same parent clade; its length is parent age minus child age.
Top row: node ages (linear). Bottom row: branch lengths (log). --ages-only: node ages only.
Usage: branch_scatter.py OUT.png OUT.tsv  NAME FAST.extant.tre FULL.extant.tre  [NAME FAST FULL ...]
FULL may also be a NEXUS tree with fossils (e.g. an original MCC tree): tips older than the youngest tip by >0.01 Myr are pruned first."""
import math, statistics as st, sys
import dendropy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def branches(path):
    t = dendropy.Tree.get(path=path, schema='nexus', preserve_underscores=True, rooting='force-rooted')
    t.calc_node_ages(is_force_max_age=True)
    old = [l.taxon.label for l in t.leaf_node_iter() if l.age > 0.01]
    if old:
        t.prune_taxa_with_labels(old)
        t.suppress_unifurcations()
        t.calc_node_ages(is_force_max_age=True)
    out = {}
    for nd in t.postorder_node_iter():
        nd.tips = frozenset([nd.taxon.label]) if nd.is_leaf() else frozenset().union(*(c.tips for c in nd.child_nodes()))
    for nd in t.preorder_node_iter():
        if nd.parent_node is not None:
            out[(nd.tips, nd.parent_node.tips)] = nd.parent_node.age - nd.age
    return out, {nd.tips: nd.age for nd in t.internal_nodes()}


def branches_from_logs(spec):
    """spec = 'ffbdlog:INPUT.txt:LOG1,LOG2': node ages from the ffbd log files (every logged sample, 10% burn-in per chain),
    clades from the extant nodes of the ffbd input file. Uses far more samples than the saved trees."""
    _, inp, logs = spec.split(':', 2)
    L = open(inp).read().split('\n')
    i = L.index(next(l for l in L if l.startswith('tips ')))
    nt = int(L[i].split()[1]); name = {}
    for l in L[i + 1:i + 1 + nt]:
        p = l.split(); name[int(p[0])] = (p[1], p[2])
    j = L.index(next(l for l in L if l.startswith('initial_tree ')))
    nn = int(L[j].split()[1]); kids = {}
    for l in L[j + 1:j + 1 + nn]:
        p = l.split(); v, a, b = int(p[0]), int(p[2]), int(p[3])
        if a >= 0: kids[v] = (a, b)
    k = L.index(next(l for l in L if l.startswith('extant_nodes ')))
    ext = {int(l.split()[0]): l.split()[1] for l in L[k + 1:k + 1 + int(L[k].split()[1])]}
    tipset = {}
    def tips(v):
        if v not in tipset:
            tipset[v] = frozenset([name[v][0]]) if v not in kids else (tips(kids[v][0]) | tips(kids[v][1]))
        return tipset[v]
    ext_tips = {v: frozenset(x for x in tips(v) if name_type[x] == 'E') for v in ext} if False else None
    etype = {nm: t for nm, t in name.values()}
    clade = {nm: frozenset(x for x in tips(v) if etype[x] == 'E') for v, nm in ext.items()}
    vals = {nm: [] for nm in ext.values()}
    for lf in logs.split(','):
        rows = [l.rstrip('\n').split('\t') for l in open(lf)]
        h, rows = rows[0], rows[1:]
        rows = rows[len(rows) // 10:]
        for nm in vals:
            c = h.index('age_' + nm)
            vals[nm] += [float(r[c]) for r in rows if len(r) == len(h)]
    age = {clade[nm]: st.median(v) for nm, v in vals.items()}
    out = {}
    for cl, a in age.items():  # parent = smallest strictly larger extant clade
        par = min((c for c in age if cl < c), key=len, default=None)
        if par is not None:
            out[(cl, par)] = age[par] - a
        if True:
            pass
    for nm_, t in etype.items():  # terminal branches: tip age 0
        if t == 'E':
            s = frozenset([nm_]); par = min((c for c in age if s < c), key=len)
            out[(s, par)] = age[par]
    return out, age


def corr(x, y):
    mx, my = st.mean(x), st.mean(y)
    return sum((p - mx) * (q - my) for p, q in zip(x, y)) / math.sqrt(sum((p - mx) ** 2 for p in x) * sum((q - my) ** 2 for q in y))


args = sys.argv[1:]
AGES_ONLY = args[0] == '--ages-only'
if AGES_ONLY:
    args = args[1:]
png, tsv, sets = args[0], args[1], [args[i:i + 3] for i in range(2, len(args), 3)]
fig, axes = plt.subplots(1 if AGES_ONLY else 2, len(sets), figsize=(4.6 * len(sets), 4.6 if AGES_ONLY else 9.0), squeeze=False)
with open(tsv, 'w') as fh:
    fh.write('set\tquantity\tn\tr\tmedian_abs_pct\tshare_within_10pct\tmean_signed_Myr\n')
    for j, (name, fast, full) in enumerate(sets):
        (a, aa), (b, ba) = (branches_from_logs(fast) if fast.startswith('ffbdlog:') else branches(fast)), branches(full)
        # row 1: node ages (linear), every clade in both trees
        k = [c for c in aa if c in ba and ba[c] > 0]
        x, y = [ba[c] for c in k], [aa[c] for c in k]
        r = corr(x, y); pc = [abs(q - p) / p * 100 for p, q in zip(x, y)]; ms = sum(q - p for p, q in zip(x, y)) / len(k)
        fh.write(f'{name}\tnode_age\t{len(k)}\t{r:.4f}\t{st.median(pc):.2f}\t{sum(v <= 10 for v in pc) / len(pc):.3f}\t{ms:.3f}\n')
        ax = axes[0][j]; hi = max(x + y) * 1.05
        ax.plot([0, hi], [0, hi], ls='--', lw=1, color='#888888', zorder=1)
        ax.scatter(x, y, s=10 if len(x) > 200 else 22, alpha=0.55 if len(x) > 200 else 0.8, color='#2a78d6', edgecolors='none', zorder=2)
        ax.set_xlim(0, hi); ax.set_ylim(0, hi)
        ax.set_title(f'{name}: node ages\n{len(k)} nodes, r = {r:.3f}, Fast {ms:+.1f} Myr on average', fontsize=10)
        ax.set_xlabel('Full FBD node age (Ma)'); ax.set_ylabel('Fast FBD (ffbd) node age (Ma)')
        print(name, 'node ages', len(k), 'r', round(r, 4), 'median |%|', round(st.median(pc), 2), 'mean signed', round(ms, 2), flush=True)
        if AGES_ONLY:
            axes[0][j].grid(True, color='#e6e6e6', lw=0.6); axes[0][j].set_axisbelow(True)
            continue
        # row 2: branch lengths (log)
        keys = [q for q in a if q in b and a[q] > 1e-3 and b[q] > 1e-3]
        x, y = [b[q] for q in keys], [a[q] for q in keys]
        rl = corr([math.log(v) for v in x], [math.log(v) for v in y]); pc = [abs(q - p) / p * 100 for p, q in zip(x, y)]
        fh.write(f'{name}\tbranch_length\t{len(keys)}\t{rl:.4f}\t{st.median(pc):.2f}\t{sum(v <= 10 for v in pc) / len(pc):.3f}\t{sum(q - p for p, q in zip(x, y)) / len(keys):.3f}\n')
        ax = axes[1][j]; lo, hi = min(x + y) * 0.8, max(x + y) * 1.25
        ax.plot([lo, hi], [lo, hi], ls='--', lw=1, color='#888888', zorder=1)
        ax.scatter(x, y, s=10 if len(x) > 200 else 22, alpha=0.55 if len(x) > 200 else 0.8, color='#2a78d6', edgecolors='none', zorder=2)
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_title(f'{name}: branch lengths\n{len(keys)} branches, r(log) = {rl:.3f}, median diff {st.median(pc):.1f}%', fontsize=10)
        ax.set_xlabel('Full FBD branch length (Myr)'); ax.set_ylabel('Fast FBD (ffbd) branch length (Myr)')
        print(name, 'branches', len(keys), 'r(log)', round(rl, 4), 'median |%|', round(st.median(pc), 2), flush=True)
        for row in (0, 1):
            axes[row][j].grid(True, which='major', color='#e6e6e6', lw=0.6); axes[row][j].set_axisbelow(True)
fig.tight_layout()
fig.savefig(png, dpi=170)
print('wrote', png, tsv)
