#!/usr/bin/env python3
"""Density equality between ffbd and RevBayes (validation step 1).
  density_check.py gen DUMPFILE DATA_DIR OUTDIR RB_DATA_DIR  -> OUTDIR/state_K.Rev + tree_K.nwk (paths for the cluster)
  density_check.py compare DUMPFILE RBOUT                     -> table; exits nonzero if |diff| > tol
RBOUT: concatenated stdout of rb runs containing lines "DENS<TAB>k<TAB>fbd root rates_mean_rates hyper data win".
"""
import sys, os, pathlib, re

def read_dump(p):
    out = []
    for line in open(p):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 15:
            continue
        d = {"it": parts[0], "lam": parts[1], "mu": parts[2], "psi": parts[3], "cond": parts[4], "m": parts[5],
             "rates": parts[6].split(",") if parts[6] else [], "fbd": float(parts[7]), "root": float(parts[8]),
             "hyper": float(parts[9]), "rate_prior": float(parts[10]), "data": float(parts[11]), "win": float(parts[12]),
             "nsa": int(parts[13]), "nwk": parts[14], "sigma": parts[15] if len(parts) > 15 else "0.5"}
        out.append(d)
    return out

TEMPLATE = r'''
setOption("outputPrecision", "15", FALSE)
seed(1)
data_dir <- "@DATA@"
source(data_dir + "/settings.Rev")
taxa <- readTaxonData(data_dir + "/fbd_taxa.tsv")
n_taxa <- taxa.size()
initial_tree <- readTrees("@TREE@", treetype="clock")[1]
speciation_rate ~ dnExponential(10)
extinction_rate ~ dnExponential(10)
psi ~ dnExponential(10)
speciation_rate.setValue(@LAM@)
extinction_rate.setValue(@MU@)
psi.setValue(@PSI@)
root_time_real ~ dnNormal(root_mean, root_sd)
root_time_real.setValue(initial_tree.rootAge())
cond_age := abs(root_time_real)
fbd_dist = dnBDSTP(rootAge=cond_age, lambda=speciation_rate, mu=extinction_rate, psi=psi, rho=rho, r=0, taxa=taxa)
source(data_dir + "/constraints.Rev")
backbone_tree <- readTrees(data_dir + "/backbone.nwk", treetype="clock")[1]
fbd_tree ~ dnConstrainedTopology(fbd_dist, constraints=constraints, backbone=backbone_tree, initialTree=initial_tree)
nf <- 0
for (i in 1:n_taxa) {
    if (taxa[i].getMaxAge() > 0.0) {
        nf <- nf + 1
        fossil_age[nf] := tmrca(fbd_tree, clade(taxa[i]))
        win_lo[nf] <- taxa[i].getMinAge() - 0.0001
        win_hi[nf] <- taxa[i].getMaxAge() + 0.0001
        fossil_win[nf] ~ dnUniform(fossil_age[nf] - win_hi[nf], fossil_age[nf] - win_lo[nf])
        fossil_win[nf].clamp(0.0)
    }
}
source(data_dir + "/extant_nodes.Rev")
branch_rates_mean ~ dnExponential(10.0)
branch_rates_mean.setValue(@M@)
n_branches <- n_extant_branches
rates_val <- v(@RATES@)
for (i in 1:n_branches) {
    branch_rates[i] ~ dnExp(1/branch_rates_mean)
    branch_rates[i].setValue(rates_val[i])
}
source(data_dir + "/fast_terms.Rev")
ln_fbd <- fbd_tree.lnProbability()
ln_root <- root_time_real.lnProbability()
ln_hyper <- speciation_rate.lnProbability() + extinction_rate.lnProbability() + psi.lnProbability()
ln_rp <- branch_rates_mean.lnProbability()
for (i in 1:n_branches) { ln_rp <- ln_rp + branch_rates[i].lnProbability() }
ln_data <- 0.0
for (j in 1:b_obs.size()) { ln_data <- ln_data + b_obs[j].lnProbability() }
ln_win <- 0.0
for (j in 1:nf) { ln_win <- ln_win + fossil_win[j].lnProbability() }
nsa <- fbd_tree.numSampledAncestors()
print("DENS", @K@, ln_fbd, ln_root, ln_hyper, ln_rp, ln_data, ln_win, nsa, sep="\t")
q()
'''


BEARS_CONS = '''constraint_ursinae = clade("Melursus_ursinus", "Ursus_arctos", "Ursus_maritimus", "Helarctos_malayanus", "Ursus_americanus", "Ursus_thibetanus", "Ursus_abstrusus", "Ursus_spelaeus")
constraints = v(constraint_ursinae)'''
BEARS = TEMPLATE.replace('source(data_dir + "/settings.Rev")\ntaxa <- readTaxonData(data_dir + "/fbd_taxa.tsv")', 'rho <- 1.0\ntaxa <- readTaxonData(data_dir + "/bears_taxa_specimen.tsv")')
BEARS = BEARS.replace('root_time_real ~ dnNormal(root_mean, root_sd)\nroot_time_real.setValue(initial_tree.rootAge())\ncond_age := abs(root_time_real)\nfbd_dist = dnBDSTP(rootAge=cond_age,', 'root_time_real ~ dnUniform(37.0, 55.0)\nroot_time_real.setValue(@COND@)\ncond_age := abs(root_time_real)\nfbd_dist = dnBDSTP(originAge=cond_age,')
BEARS = BEARS.replace('source(data_dir + "/constraints.Rev")', BEARS_CONS)
CAT = TEMPLATE.replace('taxa <- readTaxonData(data_dir + "/fbd_taxa.tsv")', 'taxa <- readTaxonData(data_dir + "/taxa.tsv")')
CAT = CAT.replace('root_time_real ~ dnNormal(root_mean, root_sd)', 'root_time_real ~ dnUniform(root_min, root_max)')
CAT = CAT.replace('source(data_dir + "/constraints.Rev")', 'source(data_dir + "/extant_clades.Rev")\nsource(data_dir + "/constraints.Rev")')
CAT = CAT.replace('source(data_dir + "/extant_nodes.Rev")\n', 'source(data_dir + "/extant_nodes.Rev")\nsource(data_dir + "/calibrations.Rev")\n')
a = CAT.index('branch_rates_mean ~ dnExponential(10.0)'); b = CAT.index('ln_fbd <- fbd_tree.lnProbability()')
CAT = CAT[:a] + '''moves = VectorMoves()
ucln_mean ~ dnExponential(10.0)
ucln_sigma ~ dnExponential(3.0)
ucln_sigma.setValue(@SIGMA@)
ucln_mu := ln(ucln_mean) - ucln_sigma * ucln_sigma * 0.5
source(data_dir + "/fast_terms_rates.Rev")
ucln_mean.setValue(@M@)
n_branches <- n_extant_branches
@RATESLINES@
for (i in 1:n_branches) { branch_rates[i].setValue(rates_val[i]) }
''' + CAT[b:]
CAT = CAT.replace('ln_rp <- branch_rates_mean.lnProbability()', 'ln_rp <- ucln_mean.lnProbability() + ucln_sigma.lnProbability()')
CAT = CAT.replace('for (j in 1:nf) { ln_win <- ln_win + fossil_win[j].lnProbability() }', 'for (j in 1:nf) { ln_win <- ln_win + fossil_win[j].lnProbability() }\nln_win <- ln_win + cal_otophysi.lnProbability() + cal_stem_characiformes.lnProbability()')
TEMPLATES = {"fish": TEMPLATE, "bears": BEARS, "cat": CAT}

if sys.argv[1] == "gen":
    dump, data, outd = sys.argv[2:5]
    kind = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] in TEMPLATES else "fish"
    T = TEMPLATES[kind]
    rows = read_dump(dump)
    o = pathlib.Path(outd); o.mkdir(parents=True, exist_ok=True)
    for k, r in enumerate(rows):
        (o / f"tree_{k}.nwk").write_text(r["nwk"] + "\n")
        txt = T.replace("@DATA@", data).replace("@TREE@", str(o / f"tree_{k}.nwk")).replace("@LAM@", r["lam"]).replace("@MU@", r["mu"]).replace("@PSI@", r["psi"]).replace("@M@", r["m"]).replace("@K@", str(k)).replace("@RATES@", ",".join(r["rates"])).replace("@COND@", r["cond"]).replace("@SIGMA@", r["sigma"]).replace("@RATESLINES@", "\n".join(f"rates_val[{i+1}] <- {v}" for i, v in enumerate(r["rates"])))
        (o / f"state_{k}.Rev").write_text(txt)
    print("wrote", len(rows), "state scripts to", o)
else:
    dump, rbout = sys.argv[2:4]
    rows = read_dump(dump)
    res = {}
    for line in open(rbout, errors="replace"):
        if line.startswith("DENS"):
            parts = line.split()
            k = int(parts[1])
            res[k] = [float(x) for x in parts[2:9]]
    names = ["fbd", "root", "hyper", "rate_prior", "data", "win"]
    worst = {n: 0.0 for n in names}; worstt = 0.0
    print("state\tnsa_ffbd/nsa_rb\t" + "\t".join("d_" + n for n in names) + "\td_total")
    nr = 0
    for k, r in enumerate(rows):
        if k not in res:
            print(k, "MISSING RevBayes result"); continue
        nr += 1
        v = res[k]
        ds = [r[n] - v[i] for i, n in enumerate(names)]
        for n, d in zip(names, ds):
            if abs(d) > abs(worst[n]): worst[n] = d
        tot = sum(ds)
        if abs(tot) > abs(worstt): worstt = tot
        print(f"{k}\t{r['nsa']}/{int(v[6])}\t" + "\t".join(f"{d:.3g}" for d in ds) + f"\t{tot:.3g}")
    print(f"states compared: {nr}/{len(rows)}; worst per component: " + " ".join(f"{n}={worst[n]:.3g}" for n in names) + f"; worst total = {worstt:.3g}")
