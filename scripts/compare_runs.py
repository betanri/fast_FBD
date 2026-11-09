#!/usr/bin/env python3
"""Compare two groups of chains parameter by parameter (ffbd vs RevBayes, or RevBayes run vs run).

  compare_runs.py --set NAME --outdir DIR --burn 0.1 --a LABEL PFX1 PFX2 [PFXn..] --b LABEL PFX1 PFX2 [..]
                  [--clades clades.tsv --b-node-trees]   read the extant node ages of group B from .trees (DNA runs log no ages)
                  [--thin N] for the tree reading

A prefix P means P.log (and P.trees if --b-node-trees).  Output: DIR/NAME_LABELa_vs_LABELb.tsv and a printed summary.  Columns
compared: every column present in both logs among node ages (age_*), root/origin age, rates, plus branch_rates[i] if present.
The %diff of the medians is reported against group B; chain-to-chain differences within each group give the noise floor.
"""
import argparse, re, sys, pathlib
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--set", required=True); ap.add_argument("--outdir", required=True); ap.add_argument("--burn", type=float, default=0.1)
ap.add_argument("--a", nargs="+", required=True); ap.add_argument("--b", nargs="+", required=True)
ap.add_argument("--clades", default=None); ap.add_argument("--b-node-trees", action="store_true")
a_args = ap.parse_args()
la, pa = a_args.a[0], a_args.a[1:]; lb, pb = a_args.b[0], a_args.b[1:]; burn = a_args.burn
out = pathlib.Path(a_args.outdir); out.mkdir(parents=True, exist_ok=True)

def ess(x):
    x = np.asarray(x, float); n = len(x)
    if n < 10 or np.std(x) == 0: return float("nan")
    x = x - x.mean(); f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f), 2 * n)[:n]; ac = ac / ac[0]
    tau = 1.0; k = 1; prev = np.inf
    while k + 1 < n:
        g = ac[k] + ac[k + 1]
        if g <= 0: break
        g = min(g, prev); prev = g; tau += 2 * g; k += 2
    return n / tau

def read_table(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
    rows = []
    with open(p) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(hdr):
                continue
            try:
                rows.append([float(x) if x not in ("nan", "NA", "inf", "-inf") else np.nan for x in parts])
            except ValueError:
                continue
    d = np.array(rows, float)
    return {k: d[:, i] for i, k in enumerate(hdr)}

TOK = re.compile(r"\(|\)|,|;|\[[^\]]*\]|:[-0-9.eE+]+|[^\[\]():,;]+")
def tree_ages(nw, clades, extant):
    """Return {clade name: age} for one newick, with age = depth of an extant tip minus depth of the MRCA of the clade's extant tips"""
    nodes = []; stack = []; last = None
    for tok in TOK.findall(nw):
        if tok == "(":
            n = {"ch": [], "len": 0.0, "name": None}
            if stack: stack[-1]["ch"].append(n)
            stack.append(n); nodes.append(n); last = n
        elif tok == ")":
            last = stack.pop()
        elif tok in (",", ";") or tok.startswith("["):
            pass
        elif tok.startswith(":"):
            last["len"] = float(tok[1:])
        else:
            n = {"ch": [], "len": 0.0, "name": tok}
            stack[-1]["ch"].append(n); nodes.append(n); last = n
    root = nodes[0]
    def rec(n, d):
        n["depth"] = d + n["len"]
        n["tips"] = frozenset([n["name"]]) if not n["ch"] else frozenset().union(*[rec(c, n["depth"]) for c in n["ch"]])
        return n["tips"]
    rec(root, 0.0)
    ref = next(n["depth"] for n in nodes if n["name"] in extant)
    res = {}
    for nm, cl in clades:
        best = min((n for n in nodes if cl <= n["tips"]), key=lambda n: len(n["tips"]))
        res[nm] = ref - best["depth"]
    return res

def load(prefix, node_from_trees, clades):
    d = read_table(prefix + ".log")
    for a_, b_ in (("root_age", "cond_age"), ("ucln_mean", "branch_rates_mean"), ("mean_rate", "mean_rt"), ("origin_time", "cond_age")):
        if a_ in d and b_ not in d: d[b_] = d[a_]
    n = len(next(iter(d.values()))); keep = int(n * burn)
    d = {k: v[keep:] for k, v in d.items()}
    if node_from_trees:
        extant = clades[-1][1]
        ages = {nm: [] for nm, _ in clades}
        with open(prefix + ".trees") as f:
            hdr = f.readline().rstrip("\n").split("\t"); k = len(hdr) - 1
            trees = [ln.rstrip("\n").split("\t") for ln in f]
        trees = [t for t in trees if len(t) > k and t[k].endswith(";")]
        keep2 = int(len(trees) * burn)
        for t in trees[keep2:]:
            r = tree_ages(t[k], clades, extant)
            for nm in ages: ages[nm].append(r[nm])
        for nm, v in ages.items():
            d[("age_" + nm)] = np.array(v)
    return d

clades = None
if a_args.clades:
    clades = []
    for line in open(a_args.clades).read().splitlines()[1:]:
        nm, k, t = line.split("\t")
        clades.append((nm, frozenset(t.split(","))))

runs_a = [load(p, False, None) for p in pa]
runs_b = [load(p, a_args.b_node_trees, clades) for p in pb]
print(f"{a_args.set}: {la} chains {len(runs_a)} samples per chain {[len(next(iter(r.values()))) for r in runs_a]}; {lb} chains {len(runs_b)} samples {[len(next(iter(r.values()))) for r in runs_b]}")

cols_a = set.intersection(*[set(r) for r in runs_a]); cols_b = set.intersection(*[set(r) for r in runs_b])
common = sorted(cols_a & cols_b)
special = ["cond_age", "root_age_tree", "origin_time", "speciation_rate", "extinction_rate", "psi", "diversification", "turnover", "num_samp_anc", "branch_rates_mean", "mean_rt", "mean_rate", "Posterior", "Likelihood"]
groups = [("ages", sorted(c for c in common if c.startswith("age_"))), ("params", [c for c in special if c in common]), ("rates", sorted((c for c in common if c.startswith("branch_rates[")), key=lambda s: int(s[13:-1])))]

def summ(x):
    return float(np.median(x)), float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))

rows = ["group\tparam\tmed_%s\tlo_%s\thi_%s\tmed_%s\tlo_%s\thi_%s\tpct_diff_medians\twidth_ratio_%s_over_%s\t%s_chain_diff_pct\t%s_chain_diff_pct\tess_min_%s\tess_min_%s" % (la, la, la, lb, lb, lb, la, lb, la, lb, la, lb)]
report = {}
for gname, cols in groups:
    stats = []
    for c in cols:
        A = np.concatenate([r[c] for r in runs_a]); B = np.concatenate([r[c] for r in runs_b])
        ma, la2, ha = summ(A); mb, lb2, hb = summ(B)
        pct = 100 * (ma - mb) / mb if mb != 0 else np.nan
        wr = (ha - la2) / (hb - lb2) if hb > lb2 else np.nan
        mas = [float(np.median(r[c])) for r in runs_a]; mbs = [float(np.median(r[c])) for r in runs_b]
        da = 100 * (max(mas) - min(mas)) / abs(ma) if ma != 0 else np.nan
        db = 100 * (max(mbs) - min(mbs)) / abs(mb) if mb != 0 else np.nan
        ea = min(ess(r[c]) for r in runs_a); eb = min(ess(r[c]) for r in runs_b)
        rows.append(f"{gname}\t{c}\t{ma:.5g}\t{la2:.5g}\t{ha:.5g}\t{mb:.5g}\t{lb2:.5g}\t{hb:.5g}\t{pct:.3f}\t{wr:.3f}\t{da:.3f}\t{db:.3f}\t{ea:.0f}\t{eb:.0f}")
        stats.append((c, abs(pct), wr, da, db, ea, eb))
    report[gname] = stats

fn = out / f"{a_args.set}_{la}_vs_{lb}.tsv"
fn.write_text("\n".join(rows) + "\n")
print("wrote", fn)
# summary of ages
for gname in ("ages", "params", "rates"):
    st = report.get(gname, [])
    if not st: continue
    pct = np.array([s[1] for s in st]); wr = np.array([s[2] for s in st]); db = np.array([s[4] for s in st])
    print(f"--- {gname}: n={len(st)}  |%diff medians|: median {np.nanmedian(pct):.3f}, max {np.nanmax(pct):.3f} ({st[int(np.nanargmax(pct))][0]}); "
          f"#>1%: {int(np.sum(pct > 1.0))}; width ratio: median {np.nanmedian(wr):.3f}, range {np.nanmin(wr):.3f}-{np.nanmax(wr):.3f}; "
          f"{lb} chain-to-chain %diff median {np.nanmedian(db):.3f} max {np.nanmax(db):.3f}; ESS(min over chains) {la}: {int(np.nanmin([s[5] for s in st]))}, {lb}: {int(np.nanmin([s[6] for s in st]))}")
