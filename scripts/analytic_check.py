#!/usr/bin/env python3
"""Prior-only analytic check (validation step 2): small tree with fixed FBD parameters and fixed root, compare the MCMC
marginals (P(fossil is a sampled ancestor), mean fossil age, mean node age, attachment edge frequencies) against numerical
quadrature of the exact density.  Usage: analytic_check.py FFBD OUTDIR"""
import sys, re, subprocess, os, math
import numpy as np
ffbd, outd = sys.argv[1], sys.argv[2]
os.makedirs(outd, exist_ok=True)
lam, mu, psi, rho, T = 0.1, 0.05, 0.05, 0.5, 20.0
w1, w2 = 3.0, 7.0
inp = f"""FFBD 1
cond root
mode control
n 3
m 1
rho {rho}
root_prior normal 20 1000
origin_prior none
init_lambda {lam}
init_mu {mu}
init_psi {psi}
init_cond_age {T}
init_clock 0.05
tips 4
0 A E 0 0
1 B E 0 0
2 C E 0 0
3 FOSSIL_F F {w1} {w2}
initial_tree 7
0 6 -1 -1 0 0
1 5 -1 -1 0 0
2 4 -1 -1 0 0
3 6 -1 -1 5 0
4 -1 5 2 20 0
5 4 6 1 12 0
6 5 3 0 8 0
constraints 1
cl_F 3 0 1 3
extant_nodes 2
5 n1
4 n2
branches 0
terms 0
"""
p = os.path.join(outd, "input.txt")
open(p, "w").write(inp)
pre = os.path.join(outd, "run")
res = subprocess.run([ffbd, "run", p, "--prefix", pre, "--seed", "5", "--gens", "300000000", "--print", "2000", "--tree-print", "2000", "--burnin", "20000", "--mode", "control", "--fix-fbd", "--log-fossils"], capture_output=True, text=True)
print(res.stdout[-900:]); print(res.stderr[-300:])
d = np.genfromtxt(pre + ".log", names=True, delimiter="\t", deletechars="")
a = d["age_FOSSIL_F"]; s = d["age_n1"]; sa = d["sa_FOSSIL_F"]
# edge frequencies from trees (F on edge A, on edge B, or on the stem of (A,B))
stem = 0; ea = 0; eb = 0; n = 0
ind_a = []; ind_b = []; ind_stem = []
for line in open(pre + ".trees"):
    if line.startswith("Iter"): continue
    nw = line.rstrip("\n").split("\t")[-1]
    # find sibling of FOSSIL_F: compare what clade groups with it via the shape of the newick around F
    # parse newick into nested tuples
    def parse(s):
        s = s.rstrip(";"); i = 0
        def node():
            nonlocal i
            if s[i] == "(":
                i += 1; L = node(); assert s[i] == ","; i += 1; R = node(); assert s[i] == ")"; i += 1
                m = re.match(r":[-0-9.eE+]+", s[i:]); ln = 0.0
                if m: ln = float(m.group(0)[1:]); i += len(m.group(0))
                return ("I", L, R, ln)
            m = re.match(r"[A-Za-z0-9_]+", s[i:]); name = m.group(0); i += len(name)
            m = re.match(r":[-0-9.eE+]+", s[i:]); ln = 0.0
            if m: ln = float(m.group(0)[1:]); i += len(m.group(0))
            return ("T", name, ln)
        return node()
    t = parse(nw)
    def tips(n): return {n[1]} if n[0] == "T" else tips(n[1]) | tips(n[2])
    # find F's sibling tip set
    sib = None
    def rec(n):
        global sib
        if n[0] == "I":
            for x, y in ((n[1], n[2]), (n[2], n[1])):
                if x[0] == "T" and x[1] == "FOSSIL_F": sib = tips(y)
            rec(n[1]); rec(n[2])
    rec(t)
    n += 1
    ind_a.append(sib == {"A"}); ind_b.append(sib == {"B"}); ind_stem.append(sib == {"A", "B"})
    if sib == {"A"}: ea += 1
    elif sib == {"B"}: eb += 1
    elif sib == {"A", "B"}: stem += 1
    else: raise SystemExit("unexpected sibling " + str(sib))
# quadrature of the exact density
A = math.sqrt((lam - mu - psi) ** 2 + 4 * lam * psi)
B = ((2 * rho - 1) * lam + mu + psi) / A
def lnD(t):
    e = np.exp(-A * t)
    return 2 * np.log(2) - A * t - 2 * np.log(1 + B + e * (1 - B))
def E(t):
    e = np.exp(-A * t)
    r = (1 + B - e * (1 - B)) / (1 + B + e * (1 - B))
    return (lam + mu + psi - A * r) / (2 * lam)
Ns = 3000; Na = 600
sg = (np.arange(Ns) + 0.5) * T / Ns; ds = T / Ns
ag = w1 + (np.arange(Na) + 0.5) * (w2 - w1) / Na; da = (w2 - w1) / Na
# cumulative of D on a fine grid (for ∫_a^b D(p) dp) by trapezoid via a dense grid
fine = np.linspace(0, T, 400001); Dfine = np.exp(lnD(fine)); cum = np.concatenate([[0.0], np.cumsum(0.5 * (Dfine[1:] + Dfine[:-1]) * np.diff(fine))])
def intD(lo, hi):
    if np.ndim(lo) == 0:
        lo = np.array([lo]); hi = np.array([hi])
    ci = np.interp(hi, fine, cum) - np.interp(lo, fine, cum)
    return ci
Es, Ea = E(ag), lnD(ag)
Ds = np.exp(lnD(sg))
Z = mSA = ma = ms = eAB = eStem = 0.0
for si, s_val in enumerate(sg):
    sw = lam * Ds[si]            # the node n1 bifurcation: lambda D(s)
    # sampled ancestor weight: psi x number of edges spanning a  (A,B when a<s; stem when a>s)
    nsae = np.where(ag < s_val, 2.0, 1.0)
    wsa = psi * nsae
    # non-SA: psi E(a)/D(a) * lambda * 2 * integral of D(p) over the attachment range
    non = psi * Es * np.exp(-Ea) * lam * 2.0
    # edge A or B: a<s -> p in (a,s)
    edgeAB = np.where(ag < s_val, intD(ag, np.full(Na, s_val)), 0.0)
    # stem: p in (max(a,s), T)
    lo = np.maximum(ag, s_val)
    stemI = intD(lo, np.full(Na, T))
    wns_each_AB = non * edgeAB
    wns_stem = non * stemI
    total = wsa + 2 * wns_each_AB + wns_stem   # A and B edges both available (factor 2) ; stem one
    Z += sw * ds * total.sum() * da
    mSA += sw * ds * wsa.sum() * da
    ma += sw * ds * (total * ag).sum() * da
    ms += sw * ds * total.sum() * da * s_val
    # edge posterior probabilities (sibling tip of the fossil): non-SA attachments and SA on the edge share the same sibling identity?  In the SA case the
    # fossil is sampled ancestor on the edge above A (or B, or stem): sibling tip set is the same as the edge's lower node.
    eAB += sw * ds * ((np.where(ag < s_val, psi, 0.0) + wns_each_AB).sum() * da)   # either edge A or B counted once = one edge
    eStem += sw * ds * ((np.where(ag >= s_val, psi, 0.0) + wns_stem).sum() * da)
PSA = mSA / Z; Ea_q = ma / Z; Es_q = ms / Z
Pst = eStem / Z; PAB1 = eAB / Z
print(f"SAMPLES {len(a)}; trees {n}")
rows_old = [("P(SA)", sa.mean(), PSA), ("mean fossil age", a.mean(), Ea_q), ("mean age_n1", s.mean(), Es_q),
        ("P(on edge A)", ea / n, PAB1), ("P(on edge B)", eb / n, PAB1), ("P(on stem of AB)", stem / n, Pst)]
def bm_se(v):
    v = np.asarray(v, float); b = 25; k = len(v) // b
    if k < 2: return np.nan
    m = np.array([v[i*k:(i+1)*k].mean() for i in range(b)])
    return m.std(ddof=1) / math.sqrt(b)
ta = np.array(a); ts = np.array(s); tsa = np.array(sa)
# per-sample edge indicators from trees were accumulated above; re-derive batch SEs for the indicators
inds = {"P(on edge A)": ind_a, "P(on edge B)": ind_b, "P(on stem of AB)": ind_stem}
rows = [("P(SA)", tsa, PSA), ("mean fossil age", ta, Ea_q), ("mean age_n1", ts, Es_q), ("P(on edge A)", ind_a, PAB1), ("P(on edge B)", ind_b, PAB1), ("P(on stem of AB)", ind_stem, Pst)]
print("quantity\tmcmc\tquadrature\tdiff\tbatch_se\tz")
ok = True
for nm, arr, q in rows:
    v = float(np.mean(arr)); se = bm_se(arr); z = (v - q) / se
    print(f"{nm}\t{v:.5f}\t{q:.5f}\t{v-q:+.5f}\t{se:.5f}\t{z:+.2f}")
    if abs(z) > 4: ok = False
print("ANALYTIC_CHECK", "PASS" if ok else "CHECK (some |z|>4; see MCSE note in VALIDATION.md)")
