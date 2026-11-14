#!/usr/bin/env python3
"""Simulation-based calibration of ffbd with rate shifts (skyline FBD); same design as sbc.py.

  sbc_sky.py run NREP OUTDIR FFBD FFBD_SIM SEED --shifts T1,T2 [--vary lambda,mu,psi] [--origin] [--no-data] [--gens G] [--print P] [--burnin B]
  sbc_sky.py analyze OUTDIR

ffbd_sim draws lambda_i, mu_i, psi_i ~ Exp(10) per interval (shared rates: one draw), the tree and the phylogram exactly from the model,
conditional on the tip configuration, and writes an input whose start is the truth. ffbd samples the posterior; the rank of each true
value among the posterior draws must be uniform. Results: OUTDIR/ranks.tsv.
"""
import argparse, os, random, subprocess, sys, pathlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sbc import chisq_sf, ks_uniform

BASE = ["cond_age", "root_age_tree", "age_n1", "age_n2", "num_samp_anc", "age_FOSSIL_1", "age_FOSSIL_2", "age_FOSSIL_3", "branch_rates_mean", "mean_rt"]
RN = ["speciation_rate", "extinction_rate", "psi"]

def rate_params(nI, vary):
    out = []
    for w, nm in enumerate(RN):
        if vary[w]: out += [f"{nm}[{i + 1}]" for i in range(nI)]
        else: out.append(nm)
    return out

def run_one(a, inpath, prefix, seed, rep, vary, nI):
    cmd = [a.sim, inpath, str(seed), str(rep), "origin" if a.origin else "root", a.shifts]
    if a.vary != "lambda,mu,psi": cmd.append(a.vary)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode or not r.stdout.startswith("TRUTH"):
        raise RuntimeError("sim failed: " + r.stderr[:300] + r.stdout[:300])
    t = r.stdout.split()
    truth = {"cond_age": float(t[6]), "branch_rates_mean": float(t[7]), "age_n1": float(t[8]), "age_n2": float(t[9]),
             "age_FOSSIL_1": float(t[11]), "age_FOSSIL_2": float(t[12]), "age_FOSSIL_3": float(t[13]), "mean_rt": float(t[14]),
             "num_samp_anc": float(t[15]), "root_age_tree": float(t[16])}
    rv = [float(x) for x in t[17:17 + 3 * nI]]
    for w, nm in enumerate(RN):
        vals = rv[w * nI:(w + 1) * nI]
        if vary[w]:
            for i in range(nI): truth[f"{nm}[{i + 1}]"] = vals[i]
        else: truth[nm] = vals[0]
    cmd = [a.ffbd, "run", inpath, "--prefix", prefix, "--seed", str(seed + 1), "--gens", str(a.gens), "--print", str(a.printe),
           "--tree-print", str(a.gens), "--burnin", str(a.burnin), "--log-fossils"]
    if a.no_data: cmd.append("--no-data")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode: raise RuntimeError("ffbd failed: " + r.stderr[:300])
    with open(prefix + ".log") as f:
        hdr = f.readline().strip().split("\t")
        rows = [line.strip().split("\t") for line in f if line.strip()]
    col = {h: i for i, h in enumerate(hdr)}
    ranks = {}
    for p in PARAMS:
        vals = [float(x[col[p]]) for x in rows[1:]]  # drop the start state
        tv = truth[p]
        less = sum(v < tv for v in vals); equal = sum(v == tv for v in vals)
        ranks[p] = (less + random.randint(0, equal), len(vals) + 1)
    ranks["tries"] = (int(t[2]), 1)
    return ranks

def analyze(outd):
    lines = open(os.path.join(outd, "ranks.tsv")).read().splitlines()
    hdr = lines[0].split("\t"); params = hdr[3:]
    data = {h: [] for h in hdr}
    for ln in lines[1:]:
        for h, v in zip(hdr, ln.split("\t")): data[h].append(float(v))
    n = len(lines) - 1; bins = max(5, min(20, n // 25))
    print(f"replicates {n}, bins {bins}")
    print("param\tchi2\tp_chi2\tKS_D\tp_KS\tmean_rank\tbin_counts")
    worst = 1.0; worstks = 1.0
    for p in params:
        us = [ri / Li for ri, Li in zip(data[p], data["L"])]
        counts = [0] * bins
        for u in us: counts[min(bins - 1, int(u * bins))] += 1
        e = n / bins; chi = sum((c - e) ** 2 / e for c in counts); pv = chisq_sf(chi, bins - 1)
        d, ks = ks_uniform(us); worst = min(worst, pv); worstks = min(worstks, ks)
        print(f"{p}\t{chi:.2f}\t{pv:.3f}\t{d:.3f}\t{ks:.3f}\t{sum(us) / n:.3f}\t{','.join(map(str, counts))}")
    print(f"smallest p-value: chi2 {worst:.4f}, KS {worstks:.4f} ({len(params)} quantities; Bonferroni threshold 0.05/{len(params)} = {0.05 / len(params):.4f})")

if __name__ == "__main__":
    if sys.argv[1] == "analyze":
        analyze(sys.argv[2]); sys.exit(0)
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("nrep", type=int); ap.add_argument("outdir"); ap.add_argument("ffbd"); ap.add_argument("sim"); ap.add_argument("seed", type=int)
    ap.add_argument("--shifts", required=True); ap.add_argument("--vary", default="lambda,mu,psi"); ap.add_argument("--origin", action="store_true")
    ap.add_argument("--no-data", action="store_true"); ap.add_argument("--gens", type=int, default=1000000); ap.add_argument("--print", dest="printe", type=int, default=5000); ap.add_argument("--burnin", type=int, default=50000)
    a = ap.parse_args()
    nI = len(a.shifts.split(",")) + 1
    vary = [w in a.vary.split(",") for w in ["lambda", "mu", "psi"]]
    PARAMS = rate_params(nI, vary) + BASE
    o = pathlib.Path(a.outdir); o.mkdir(parents=True, exist_ok=True)
    random.seed(a.seed)
    ranksf = open(o / "ranks.tsv", "w")
    ranksf.write("rep\ttries\tL\t" + "\t".join(PARAMS) + "\n")
    for rep in range(a.nrep):
        inp = str(o / f"tmp_{a.seed}.txt"); pre = str(o / f"tmp_{a.seed}")
        rnk = run_one(a, inp, pre, a.seed + 2 * rep + 1, rep, vary, nI)
        L = rnk[PARAMS[0]][1]
        ranksf.write(f"{rep}\t{rnk['tries'][0]}\t{L}\t" + "\t".join(str(rnk[p][0]) for p in PARAMS) + "\n"); ranksf.flush()
    for suf in [".txt", ".log", ".trees", ".ckp"]:
        try: os.remove(str(o / f"tmp_{a.seed}{suf}"))
        except OSError: pass
    ranksf.close()
    analyze(str(o))
