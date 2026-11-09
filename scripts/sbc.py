#!/usr/bin/env python3
"""Simulation-based calibration for ffbd (validation step 2).

  sbc.py run  NREP OUTDIR FFBD FFBD_SIM SEED [--no-data] [--rb-sa-move] [--gens G] [--print P] [--burnin B]
  sbc.py analyze OUTDIR

For each replicate: ffbd_sim draws (rates, tree, data) exactly from the prior conditional on the tip configuration and
writes an input whose start is the true tree; ffbd then samples the posterior; we record the rank of the true value among
the posterior draws (ties broken at random).  Uniform ranks = correct moves and density.  Results: OUTDIR/ranks.tsv.
"""
import argparse, os, random, subprocess, sys, pathlib, math

PARAMS = ["speciation_rate", "extinction_rate", "psi", "cond_age", "branch_rates_mean", "age_n1", "age_n2", "mean_rt", "num_samp_anc",
          "age_FOSSIL_1", "age_FOSSIL_2", "age_FOSSIL_3", "diversification", "root_age_tree"]

def run_one(ffbd, sim, inpath, prefix, seed, rep, gens, printe, burnin, nodata, rbsa, origin=False):
    r = subprocess.run([sim, inpath, str(seed), str(rep)] + (["origin"] if origin else []), capture_output=True, text=True)
    if r.returncode or not r.stdout.startswith("TRUTH"):
        raise RuntimeError("sim failed: " + r.stderr[:300] + r.stdout[:300])
    t = r.stdout.split()
    # TRUTH rep tries lam mu psi root m n1 n2 n3 F1 F2 F3 meanR nsa
    truth = {"speciation_rate": float(t[3]), "extinction_rate": float(t[4]), "psi": float(t[5]), "cond_age": float(t[6]),
             "branch_rates_mean": float(t[7]), "age_n1": float(t[8]), "age_n2": float(t[9]), "age_FOSSIL_1": float(t[11]),
             "age_FOSSIL_2": float(t[12]), "age_FOSSIL_3": float(t[13]), "mean_rt": float(t[14]), "num_samp_anc": float(t[15]), "root_age_tree": float(t[16])}
    truth["diversification"] = truth["speciation_rate"] - truth["extinction_rate"]
    cmd = [ffbd, "run", inpath, "--prefix", prefix, "--seed", str(seed + 1), "--gens", str(gens), "--print", str(printe),
           "--tree-print", str(gens), "--burnin", str(burnin), "--log-fossils"]
    if nodata:
        cmd.append("--no-data")
    if rbsa:
        cmd.append("--rb-sa-move")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError("ffbd failed: " + r.stderr[:300])
    with open(prefix + ".log") as f:
        hdr = f.readline().strip().split("\t")
        rows = [line.strip().split("\t") for line in f if line.strip()]
    col = {h: i for i, h in enumerate(hdr)}
    ranks = {}
    for p in PARAMS:
        vals = [float(x[col[p]]) for x in rows[1:]]  # first log line is the pre-burn-in/start state; drop it
        tv = truth[p]
        less = sum(v < tv for v in vals)
        equal = sum(v == tv for v in vals)
        rank = less + random.randint(0, equal)
        ranks[p] = (rank, len(vals) + 1)
    ranks["tries"] = (int(t[2]), 1)
    return ranks

def chisq_sf(x, df):
    # survival function of chi-square via regularized upper incomplete gamma (series/continued fraction)
    a = df / 2.0; x2 = x / 2.0
    if x2 <= 0: return 1.0
    def gser():
        s = term = 1.0 / a; n = a
        for _ in range(10000):
            n += 1; term *= x2 / n; s += term
            if abs(term) < abs(s) * 1e-15: break
        return s * math.exp(-x2 + a * math.log(x2) - math.lgamma(a))
    def gcf():
        b = x2 + 1 - a; c = 1e300; d = 1 / b; h = d
        for i in range(1, 10000):
            an = -i * (i - a); b += 2; d = an * d + b
            if abs(d) < 1e-300: d = 1e-300
            c = b + an / c
            if abs(c) < 1e-300: c = 1e-300
            d = 1 / d; delta = d * c; h *= delta
            if abs(delta - 1) < 1e-15: break
        return math.exp(-x2 + a * math.log(x2) - math.lgamma(a)) * h
    return 1 - gser() if x2 < a + 1 else gcf()

def ks_uniform(u):
    u = sorted(u); n = len(u)
    d = max(max(abs((i + 1) / n - v), abs(i / n - v)) for i, v in enumerate(u))
    # asymptotic p-value
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    p = 0.0
    for j in range(1, 101):
        p += 2 * (-1) ** (j - 1) * math.exp(-2 * j * j * lam * lam)
    return d, max(0.0, min(1.0, p))

def analyze(outd):
    lines = open(os.path.join(outd, "ranks.tsv")).read().splitlines()
    hdr = lines[0].split("\t")
    data = {h: [] for h in hdr}
    for ln in lines[1:]:
        for h, v in zip(hdr, ln.split("\t")):
            data[h].append(float(v))
    res = []
    n = len(lines) - 1
    bins = max(5, min(20, n // 25))
    for p in PARAMS:
        L = data["L_" + p] if ("L_" + p) in data else None
        r = data[p]; L_val = data["L"]
        # continuous rank in [0,1]: (rank + uniform0-1)/L is already handled by tie-randomization; bin rank into `bins` bins
        counts = [0] * bins
        us = []
        for ri, Li in zip(r, L_val):
            u = ri / Li  # rank in 0..L-1 out of L possible values, so u in [0,1)
            us.append(u)
            counts[min(bins - 1, int(u * bins))] += 1
        exp = n / bins
        chi = sum((c - exp) ** 2 / exp for c in counts)
        pv = chisq_sf(chi, bins - 1)
        d, ks = ks_uniform(us)
        res.append((p, chi, pv, d, ks, counts))
    print(f"replicates {n}, bins {bins}")
    print("param\tchi2\tp_chi2\tKS_D\tp_KS\tbin_counts")
    for p, chi, pv, d, ks, counts in res:
        print(f"{p}\t{chi:.2f}\t{pv:.3f}\t{d:.3f}\t{ks:.3f}\t{','.join(map(str, counts))}")
    worst = min(r[2] for r in res)
    print(f"smallest chi2 p-value: {worst:.4f} (with {len(res)} parameters, Bonferroni threshold 0.05/{len(res)} = {0.05/len(res):.4f})")

if __name__ == "__main__":
    if sys.argv[1] == "analyze":
        analyze(sys.argv[2]); sys.exit(0)
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("nrep", type=int); ap.add_argument("outdir"); ap.add_argument("ffbd"); ap.add_argument("sim"); ap.add_argument("seed", type=int)
    ap.add_argument("--no-data", action="store_true"); ap.add_argument("--rb-sa-move", action="store_true"); ap.add_argument("--gens", type=int, default=1000000); ap.add_argument("--print", dest="printe", type=int, default=5000); ap.add_argument("--burnin", type=int, default=50000); ap.add_argument("--origin", action="store_true")
    a = ap.parse_args()
    o = pathlib.Path(a.outdir); o.mkdir(parents=True, exist_ok=True)
    random.seed(a.seed)
    ranksf = open(o / "ranks.tsv", "w")
    ranksf.write("rep\ttries\tL\t" + "\t".join(PARAMS) + "\n")
    for rep in range(a.nrep):
        inp = str(o / f"tmp_{a.seed}.txt"); pre = str(o / f"tmp_{a.seed}")
        rnk = run_one(a.ffbd, a.sim, inp, pre, a.seed + 2 * rep + 1, rep, a.gens, a.printe, a.burnin, a.no_data, a.rb_sa_move, a.origin)
        L = rnk["speciation_rate"][1]
        ranksf.write(f"{rep}\t{rnk['tries'][0]}\t{L}\t" + "\t".join(str(rnk[p][0]) for p in PARAMS) + "\n")
        ranksf.flush()
    for suf in [".txt", ".log", ".trees", ".ckp"]:
        try: os.remove(str(o / f"tmp_{a.seed}{suf}"))
        except OSError: pass
    ranksf.close()
    analyze(str(o))
