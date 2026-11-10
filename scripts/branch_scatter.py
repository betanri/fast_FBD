#!/usr/bin/env python3
"""Branch lengths in time (Myr): full-DNA FBD against Fast FBD, one panel per data set.
Each side is summarized as an extant-only time tree (fossils pruned, conditional median node ages; FBD/TIME_TREES/extant_time_tree.py).
A branch is compared when both trees contain the same child clade with the same parent clade; its length is parent age minus child age.
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
    return out


args = sys.argv[1:]
png, tsv, sets = args[0], args[1], [args[i:i + 3] for i in range(2, len(args), 3)]
fig, axes = plt.subplots(1, len(sets), figsize=(4.6 * len(sets), 4.6))
axes = axes if len(sets) > 1 else [axes]
with open(tsv, 'w') as fh:
    fh.write('set\tn_branches\tr_log\tmedian_abs_pct\tshare_within_10pct\tmedian_full_Myr\tmedian_fast_Myr\n')
    for ax, (name, fast, full) in zip(axes, sets):
        a, b = branches(fast), branches(full)
        keys = [k for k in a if k in b and a[k] > 1e-3 and b[k] > 1e-3]
        x = [b[k] for k in keys]; y = [a[k] for k in keys]
        lx, ly = [math.log(v) for v in x], [math.log(v) for v in y]
        mx, my = st.mean(lx), st.mean(ly)
        r = sum((p - mx) * (q - my) for p, q in zip(lx, ly)) / math.sqrt(sum((p - mx) ** 2 for p in lx) * sum((q - my) ** 2 for q in ly))
        pc = [abs(q - p) / p * 100 for p, q in zip(x, y)]
        w10 = sum(v <= 10 for v in pc) / len(pc)
        fh.write(f'{name}\t{len(keys)}\t{r:.4f}\t{st.median(pc):.2f}\t{w10:.3f}\t{st.median(x):.3f}\t{st.median(y):.3f}\n')
        lo, hi = min(x + y) * 0.8, max(x + y) * 1.25
        ax.plot([lo, hi], [lo, hi], ls='--', lw=1, color='#888888', zorder=1)
        ax.scatter(x, y, s=10 if len(x) > 200 else 22, alpha=0.55 if len(x) > 200 else 0.8, color='#2a78d6', edgecolors='none', zorder=2)
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_title(f'{name}\n{len(keys)} branches, r(log) = {r:.3f}, median diff {st.median(pc):.1f}%', fontsize=10)
        ax.set_xlabel('Full FBD branch length (Myr)'); ax.set_ylabel('Fast FBD (ffbd) branch length (Myr)')
        ax.grid(True, which='major', color='#e6e6e6', lw=0.6); ax.set_axisbelow(True)
        print(name, 'branches', len(keys), 'r(log)', round(r, 4), 'median |%|', round(st.median(pc), 2), 'within 10%', round(w10, 3), flush=True)
fig.tight_layout()
fig.savefig(png, dpi=180)
print('wrote', png, tsv)
