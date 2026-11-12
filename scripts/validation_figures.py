#!/usr/bin/env python3
"""Validation figures from the per-node and per-replicate tables written by the scoring scripts.

  validation_figures.py DATADIR OUTDIR        (the repository copies: validation_figures.py validation validation)

DATADIR holds (any missing input skips its figure):
  three_A06.tsv, three_A12.tsv    compare_fixed3.py SET LOCALDIR OUT.tsv      (test 3a: ffbd, MrBayes, RevBayes; no data)
  ctrl_A06.tsv, ctrl_A12.tsv      ctrl_compare.py SET --tsv OUT.tsv           (test 4e: phylograms against 23-gene DNA)
  sbc_ranks.tsv                   sbc_ranks.py SBCDIR OUT.tsv                 (test 2b: simulation-based calibration)
  sim_study.tsv                   sim_study_summary.py SIMDIR OUT.tsv         (test 4c: simulation study)
"""
import csv, math, os, statistics as st, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

D, O = sys.argv[1], sys.argv[2]
os.makedirs(O, exist_ok=True)
BLUE, ORANGE, AQUA, GREY, INK, MUTED = '#2a78d6', '#eb6834', '#1baf7a', '#9a9993', '#0b0b0b', '#52514e'
NAMES = {'A06': 'Gobiaria', 'A12': 'Protacanthopterygii'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#bdbcb6', 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})


def rows(f):
    for p in (f'{D}/{f}', f'{D}/plotdata_{f}'):       # the copies in validation/ carry a plotdata_ prefix
        if os.path.exists(p):
            return list(csv.DictReader(open(p), delimiter='\t'))
    return None


def logx(ax):
    from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 15, 20, 30, 40, 50, 70, 100, 150, 200, 300, 500]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    ax.xaxis.set_minor_formatter(NullFormatter())


def grid(ax):
    ax.grid(True, color='#ececea', lw=0.6)
    ax.set_axisbelow(True)


def diag(ax, lo, hi):
    ax.plot([lo, hi], [lo, hi], ls='--', color=GREY, lw=1, zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)


# ---- test 3a: three programs, same fixed-rate prior, no data
three = {S: rows(f'three_{S}.tsv') for S in NAMES}
if all(three.values()):
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.4))
    for j, S in enumerate(NAMES):
        r = three[S]
        node = {}
        for x in r:
            if x['node'][1:].isdigit():
                node.setdefault(x['node'], {})[x['program']] = float(x['median'])
        node = {k: v for k, v in node.items() if v['mb'] > 0}
        mb = [v['mb'] for v in node.values()]
        ax = axes[0][j]
        grid(ax)
        ax.axhline(0, color=GREY, lw=1, ls='--', zorder=1)
        for p, lab, c, z in [('rb', 'RevBayes 1.3.2', ORANGE, 3), ('ffbd', 'ffbd', BLUE, 4)]:
            d = [100 * (v[p] - v['mb']) / v['mb'] for v in node.values()]
            ax.scatter(mb, d, s=36, color=c, edgecolor='white', lw=1.2, zorder=z,
                       label=f'{lab}: median |diff| {st.median(abs(x) for x in d):.1f}%, largest {max(abs(x) for x in d):.0f}%')
        logx(ax)
        ax.set_title(f'{NAMES[S]}: node ages against MrBayes 3.2.7', loc='left', fontsize=10, color=INK)
        ax.set_xlabel('Node age in MrBayes (Ma, log scale)')
        ax.set_ylabel('Difference from MrBayes (%)')
        ax.legend(loc='upper right', fontsize=8)
        ax = axes[1][j]
        grid(ax)
        sa = {}
        for x in r:
            if x['node'].startswith('num_samp_anc='):
                sa.setdefault(x['program'], {})[int(x['node'].split('=')[1])] = float(x['median'])
        ks = [k for k in sorted(set().union(*sa.values())) if max(sa[p].get(k, 0) for p in sa) >= 0.001]
        w = 0.27
        for i, (p, lab, c) in enumerate([('ffbd', 'ffbd', BLUE), ('mb', 'MrBayes', AQUA), ('rb', 'RevBayes', ORANGE)]):
            ax.bar([k + (i - 1) * w for k in ks], [sa[p].get(k, 0) for k in ks], width=w - 0.03, color=c, label=lab)
        ax.set_title(f'{NAMES[S]}: number of sampled ancestors', loc='left', fontsize=10, color=INK)
        ax.set_xlabel('Sampled ancestors in a posterior tree')
        ax.set_ylabel('Share of posterior samples')
        ax.set_xticks(ks)
        ax.legend(loc='upper left')
    fig.suptitle('Test 3a. Same fixed tree, fossils and fixed rates, no data: ffbd and MrBayes agree, RevBayes differs',
                 x=0.01, ha='left', fontsize=10.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(f'{O}/three_programs_no_data.png', dpi=160)
    print('wrote three_programs_no_data.png')

# ---- test 4e: controlled test, phylograms against the 23-gene full DNA
ARMS = [('fullDNA noise', 'Full DNA, chains 1-3 vs 4-6 (noise)', GREY, 'o'), ('g23', '23-gene phylogram', BLUE, 'o'),
        ('allgene', 'All-gene phylogram (~580k sites)', ORANGE, 's'), ('large_tree_a', 'Larger-tree phylogram, pruned', AQUA, '^')]
ctrl = {S: rows(f'ctrl_{S}.tsv') for S in NAMES}
if all(ctrl.values()):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.6))
    for j, S in enumerate(NAMES):
        ax = axes[j]
        grid(ax)
        ax.axhline(0, color=GREY, lw=1, ls='--', zorder=1)
        for arm, lab, c, m in ARMS:
            s = [x for x in ctrl[S] if x['arm'] == arm]
            xs = [float(x['dna_median']) for x in s]
            ys = [100 * (float(x['arm_median']) - float(x['dna_median'])) / float(x['dna_median']) for x in s]
            med = st.median(abs(y) for y in ys)
            ax.scatter(xs, ys, s=30, marker=m, color=c, edgecolor='white', lw=1, zorder=3 if arm != 'fullDNA noise' else 2,
                       label=f'{lab}: median |diff| {med:.1f}%')
        logx(ax)
        ax.set_title(NAMES[S], loc='left', fontsize=10, color=INK)
        ax.set_xlabel('Node age, 23-gene full-DNA FBD (Ma, log scale)')
        ax.set_ylabel('Difference from full DNA (%)')
        ax.legend(loc='upper left', bbox_to_anchor=(0, -0.16), fontsize=8)
    fig.suptitle('Test 4e. Same tree, fossils and priors; only the phylogram given to ffbd changes. '
                 'Reference: RevBayes full-DNA FBD, 23 genes', x=0.01, ha='left', fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(f'{O}/phylogram_vs_dna_controlled.png', dpi=160)
    print('wrote phylogram_vs_dna_controlled.png')

# ---- test 2b: SBC rank histograms
sbc = rows('sbc_ranks.tsv')
if sbc:
    LAB = {'speciation_rate': 'Speciation', 'extinction_rate': 'Extinction', 'psi': 'Fossil sampling', 'branch_rates_mean': 'Mean clock rate',
           'root': 'Root age', 'num_samp_anc': 'Sampled ancestors', 'sbc_side1': 'Crown age, root side 1',
           'sbc_side2': 'Crown age, root side 2', 'sbc_fossil': 'One fossil age'}
    fig, axes = plt.subplots(3, 3, figsize=(10, 8))
    for ax, (p, lab) in zip(axes.flat, LAB.items()):
        r = [int(x['rank']) for x in sbc if x['parameter'] == p]
        n, B = len(r), 20
        e = n / B
        sd = math.sqrt(n * (1 / B) * (1 - 1 / B))
        for y in (e - 2.58 * sd, e + 2.58 * sd):
            ax.axhline(y, color=INK, lw=0.9, ls=':', zorder=3)
        ax.axhline(e, color=INK, lw=0.9, ls='--', zorder=3)
        ax.hist(r, bins=[i * 5 for i in range(B + 1)], color=BLUE, edgecolor='white', lw=1, zorder=2)
        ax.set_title(f'{lab} (n = {n})', loc='left', fontsize=9.5, color=INK)
        ax.set_xlim(0, 100)
        ax.set_xlabel('Rank of the true value (0-99)', fontsize=8)
    fig.suptitle('Test 2b. Simulation-based calibration, 500 replicates (2 to 8,056 living species).\n'
                 'Dashed: expected count; dotted: 99% range if ranks are uniform', x=0.01, ha='left', fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(f'{O}/sbc_rank_histograms.png', dpi=160)
    print('wrote sbc_rank_histograms.png')

# ---- test 4c: simulation study, coverage and error by node age
sim = rows('sim_study.tsv')
if sim:
    BINS = [(0, 5), (5, 10), (10, 20), (20, 40), (40, 60), (60, 101)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for arm, lab, c, dx in [('exact', 'Exact phylogram (tests the sampler)', BLUE, -0.12), ('dna', 'Phylogram from simulated DNA', ORANGE, 0.12)]:
        s = [x for x in sim if x['arm'] == arm and x['group'] == 'node']
        cov, err, ns = [], [], []
        for lo, hi in BINS:
            b = [x for x in s if lo <= float(x['true']) < hi]
            ns.append(len(b))
            cov.append(st.fmean(int(x['in95']) for x in b))
            err.append(100 * st.median(float(x['relerr']) for x in b))
        xs = [i + dx for i in range(len(BINS))]
        axes[0].plot(xs, cov, marker='o', color=c, lw=2, ms=6, label=lab)
        axes[1].plot(xs, err, marker='o', color=c, lw=2, ms=6, label=lab)
    tick = [f'{lo}-{hi if hi < 100 else "100"}\n(n={n})' for (lo, hi), n in zip(BINS, ns)]
    for ax in axes:
        grid(ax)
        ax.set_xticks(range(len(BINS)))
        ax.set_xticklabels(tick, fontsize=7.5)
        ax.set_xlabel('True node age (Ma)')
    axes[0].axhline(0.95, color=GREY, ls='--', lw=1)
    axes[0].set_ylim(0.6, 1.0)
    axes[0].set_ylabel('Share of true ages inside the 95% interval')
    axes[0].set_title('Coverage of node ages', loc='left', fontsize=10, color=INK)
    axes[0].legend(loc='lower right', fontsize=8)
    axes[1].axhline(0, color=GREY, ls='--', lw=1)
    axes[1].set_ylabel('Median error of the posterior median (%)')
    axes[1].set_title('Bias of node ages (positive = too old)', loc='left', fontsize=10, color=INK)
    fig.suptitle('Test 4c. 100 simulated trees (30 to 120 living species, 6,911 nodes), known true ages',
                 x=0.01, ha='left', fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f'{O}/simulation_coverage_by_age.png', dpi=160)
    print('wrote simulation_coverage_by_age.png')
