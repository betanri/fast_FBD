#!/usr/bin/env python3
"""Convert a folder of RevBayes-style inputs into one ffbd input file.

  ffbd_convert.py --set DIR --out FILE [--mode fast|control] [--root-prior normal MEAN SD | uniform MIN MAX | soft MIN MAX SD P]
                  [--bears REVFILE] [--seqrev DIR]

DIR holds: fbd_taxa.tsv (or bears_taxa_specimen.tsv), settings.Rev, backbone.nwk, initial.nwk, constraints_members.tsv or
constraints.Rev, extant_nodes.Rev, fast_terms.Rev.  --bears gives the Rev script of the bear runs (origin prior, constraint).
Without --root-prior the prior comes from settings.Rev (root_mean, root_sd: Normal).  Only reads files; never edits the data.
"""
import argparse, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("--set", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--mode", default="fast"); ap.add_argument("--root-prior", nargs="+", default=None)
ap.add_argument("--bears", default=None)
ap.add_argument("--terms-file", default="fast_terms.Rev")
a = ap.parse_args()
d = a.set


def rd(fn):
    return open(d + "/" + fn).read()


# ---- taxa ----
import os
taxfile = "bears_taxa_specimen.tsv" if a.bears else ("fbd_taxa.tsv" if os.path.exists(d + "/fbd_taxa.tsv") else "taxa.tsv")
lines = rd(taxfile).splitlines()
hdr = lines[0].split("\t")
ix = {h: i for i, h in enumerate(hdr)}
taxa = {}
for ln in lines[1:]:
    p = ln.split("\t")
    taxa[p[ix["taxon"]]] = (float(p[ix["min_age"]]), float(p[ix["max_age"]]))
extant = sorted(t for t, (mn, mx) in taxa.items() if mx <= 0.0)
fossils = sorted(t for t, (mn, mx) in taxa.items() if mx > 0.0)
name2id = {t: i for i, t in enumerate(extant)}
name2id.update({t: len(extant) + i for i, t in enumerate(fossils)})

# ---- settings ----
sett = {"rho": None, "cond": "root", "root_age_init": None, "clock_rate_init": None, "root_mean": None, "root_sd": None, "root_min": None, "root_max": None}
init_lam, init_mu, init_psi = 0.04, 0.02, 0.001
origin_prior = None
if a.bears:
    sett["rho"] = 1.0
    sett["cond"] = "origin"
    init_cond = 50.0
    init_lam, init_mu, init_psi = 0.2, 0.1, 0.1
    sett["clock_rate_init"] = None
    origin_prior = (37.0, 55.0)
    rtxt = open(a.bears).read()
    m = re.search(r"clock_rate_init\s*<-\s*([0-9.eE+-]+)", rtxt)
    # bear clock init: read if present (the Rev script may hardwire a value), else below
else:
    txt = rd("settings.Rev")
    for key in sett:
        m = re.search(r"^\s*%s\s*<-\s*(.+?)\s*(#.*)?$" % key, txt, re.M)
        if m and key != "cond":
            v = m.group(1)
            try:
                sett[key] = float(eval(v))
            except Exception:
                sys.exit("cannot parse setting " + key + ": " + v)
        elif m:
            sett[key] = m.group(1).strip().strip('"')
    init_cond = sett["root_age_init"]
if a.root_prior:
    rp = a.root_prior
else:
    if sett["root_mean"] is not None:
        rp = ["normal", str(sett["root_mean"]), str(sett["root_sd"])]
    elif sett["root_min"] is not None:
        rp = ["uniform", str(sett["root_min"]), str(sett["root_max"])]
    else:
        rp = ["none"]
if a.bears:
    rp = ["none"]


# ---- newick parsing ----
def parse(s):
    s = s.strip().rstrip(";")
    i = 0
    def node():
        nonlocal i
        n = {"ch": [], "name": None, "len": 0.0}
        if s[i] == "(":
            i += 1
            while True:
                n["ch"].append(node())
                if s[i] == ",":
                    i += 1
                    continue
                if s[i] == ")":
                    i += 1
                    break
                raise ValueError(s[i:i + 30])
        m = re.match(r"[^,():;\[]+", s[i:])
        if m:
            if not n["ch"]:
                n["name"] = m.group(0)
            i += len(m.group(0))
        if i < len(s) and s[i] == "[":
            i = s.index("]", i) + 1
        if i < len(s) and s[i] == ":":
            i += 1
            m = re.match(r"[-0-9.eE+]+", s[i:])
            n["len"] = float(m.group(0))
            i += len(m.group(0))
        if i < len(s) and s[i] == "[":
            i = s.index("]", i) + 1
        return n
    return node()


init = parse(rd("initial.nwk"))
if len(init["ch"]) != 2:
    sys.exit("initial tree root is not bifurcating")

# ids: tips as name2id; internals numbered after tips in pre-order
nodes = []
def walk(n, parent, depth):
    n["parent"] = parent
    n["depth"] = depth + n["len"]
    nodes.append(n)
    for c in n["ch"]:
        walk(c, n, n["depth"])
init["len"] = 0.0
walk(init, None, 0.0)
N = len(taxa)
if len([x for x in nodes if not x["ch"]]) != N:
    sys.exit("initial tree tip count != taxa count")
# depth of extant tips (take the max)
exd = [x["depth"] for x in nodes if not x["ch"] and x["name"] in extant]
D = max(exd)
for x in nodes:
    x["age"] = D - x["depth"]
    if not x["ch"]:
        if x["name"] in extant:
            x["age"] = 0.0
        else:
            mn, mx = taxa[x["name"]]
            x["age"] = min(mx, max(mn, x["age"]))
next_int = N
for x in nodes:
    if x["ch"]:
        x["id"] = next_int
        next_int += 1
    else:
        x["id"] = name2id[x["name"]]
if next_int != 2 * N - 1:
    sys.exit("binary tree expected")

# keep ages consistent after the fossil clamping: parent age >= child age; SA if pendant length is exactly 0
sa_count = 0
for x in nodes:
    if x["parent"] is not None and x["age"] > x["parent"]["age"]:
        # a fossil clamped up against its parent: pull the parent up (rounding only)
        cur = x
        while cur["parent"] is not None and cur["age"] > cur["parent"]["age"]:
            cur["parent"]["age"] = cur["age"]
            cur = cur["parent"]
for x in nodes:
    x["sa"] = 0
    if x["parent"] is not None and not x["ch"] and x["name"] in fossils:
        sib = [c for c in x["parent"]["ch"] if c is not x][0]
        if x["len"] == 0.0 and sib["len"] > 0:
            x["sa"] = 1
            sa_count += 1

def tipset(n):
    if not n["ch"]:
        return {n["name"]}
    return set().union(*[tipset(c) for c in n["ch"]])
for x in nodes:
    x["tips"] = tipset(x)

# ---- constraints ----
cons = []
if a.bears:
    rtxt = open(a.bears).read()
    m = re.search(r"constraint_ursinae\s*=\s*clade\(([^)]*)\)", rtxt)
    cons.append(("constraint_ursinae", {s.strip().strip('"') for s in m.group(1).split(",")}))
else:
    try:
        ct = rd("constraints_members.tsv").splitlines()
        h = ct[0].split("\t")
        for ln in ct[1:]:
            p = ln.split("\t")
            cons.append((p[0], set(p[h.index("members")].split(","))))
    except FileNotFoundError:
        txt = rd("constraints.Rev")
        defs = {}
        if os.path.exists(d + "/extant_clades.Rev"):
            for m in re.finditer(r"^(\w+)\s*=\s*clade\((.*?)\)\s*$", rd("extant_clades.Rev"), re.M):
                sset = set()
                for arg in re.split(r",\s*", m.group(2)):
                    arg = arg.strip()
                    if arg.startswith('"'):
                        sset.add(arg.strip('"'))
                    elif arg in defs:
                        sset |= defs[arg]
                    else:
                        sys.exit("unknown reference in extant_clades.Rev: " + arg)
                defs[m.group(1)] = sset
        for m in re.finditer(r"^(\w+)\s*=\s*clade\((.*?)\)\s*$", txt, re.M):
            s = set()
            for arg in re.split(r",\s*", m.group(2)):
                arg = arg.strip()
                if arg.startswith('"'):
                    s.add(arg.strip('"'))
                elif arg in defs:
                    s |= defs[arg]
                else:
                    sys.exit("unknown clade reference " + arg)
            defs[m.group(1)] = s
        vec = re.search(r"constraints\s*=\s*v\(([^)]*)\)", txt, re.S).group(1)
        for name in re.split(r"[,\s]+", vec.strip()):
            if name:
                cons.append((name, defs[name]))

allset = set(taxa)
uniq = {}
for name, s in cons:
    if not s <= allset:
        sys.exit("constraint %s has taxa outside the set: %s" % (name, sorted(s - allset)[:3]))
    fs = frozenset(s)
    if fs == frozenset(allset) or len(s) < 2:
        continue
    uniq.setdefault(fs, name)
cons = [(nm, s) for s, nm in uniq.items()]

# check the initial tree satisfies the constraints
tsets = [x["tips"] for x in nodes]
for name, s in cons:
    if set(s) not in tsets:
        sys.exit("initial tree violates constraint " + name)

# ---- extant nodes and branches ----
ent = rd("extant_nodes.Rev")
# fish sets may have extant_nodes.Rev with clade references; resolve definitions
defs = {}
if os.path.exists(d + "/extant_clades.Rev"):
    ent = rd("extant_clades.Rev") + "\n" + ent
for m in re.finditer(r"^(clade_\w+)\s*=\s*clade\((.*?)\)\s*$", ent, re.M):
    s = set()
    for arg in re.split(r",\s*", m.group(2)):
        arg = arg.strip()
        if arg.startswith('"'):
            s.add(arg.strip('"'))
        elif arg in defs:
            s |= defs[arg]
        else:
            sys.exit("unknown reference in extant_nodes: " + arg)
    defs[m.group(1)] = s
names = sorted({m.group(1) for m in re.finditer(r"^age_(\w+)\s*:=", ent, re.M)} | {k[6:] for k in defs})
# keep the file order (matches age_n1.. style)
ordered = []
for m in re.finditer(r"^age_(\w+)\s*:=", ent, re.M):
    if m.group(1) not in ordered:
        ordered.append(m.group(1))
for k in defs:
    k2 = k[6:]
    if k2 not in ordered:
        ordered.append(k2)
ext_ids = {}
ext_nodes = []
exset = set(extant)
for nm in ordered:
    if "clade_" + nm not in defs:
        continue
    s = defs["clade_" + nm]
    # MRCA of the clade's extant tips, restricted to extant taxa (fossil constraint nodes are not extant nodes)
    target = s & exset
    best = None
    for x in nodes:
        if x["ch"] and target <= (x["tips"] & exset) and (best is None or len(x["tips"]) < len(best["tips"])):
            best = x
    if best is None:
        sys.exit("extant node not found: " + nm)
    ext_nodes.append((best["id"], nm))
    ext_ids[nm] = best["id"]
    best["extant_internal"] = True

# branches
branches = {}
for m in re.finditer(r"^branch_time\[(\d+)\]\s*:=\s*abs\(\s*age_(\w+)(?:\s*-\s*age_(\w+))?\s*\)", ent, re.M):
    i = int(m.group(1)) - 1
    pa = m.group(2)
    ch = m.group(3)
    if pa not in ext_ids or (ch and ch not in ext_ids):
        sys.exit("branch endpoint unknown " + m.group(0))
    branches[i] = (ext_ids[pa], ext_ids[ch] if ch else -1)
K = 0
if branches:
    K = max(branches) + 1
    if len(branches) != K:
        sys.exit("branch indices are not contiguous")

# ---- data terms ----
terms = []
if a.mode == "fast":
    if not K:
        sys.exit("no branch_time lines found")
    ftxt = rd(a.terms_file)
    mus = {}
    for m in re.finditer(r"^b_mu\[(\d+)\]\s*:=\s*abs\((.*)\)", ftxt, re.M):
        j = int(m.group(1))
        comps = m.group(2).split("#")[0].split("+")
        br = []
        for c in comps:
            mm = re.search(r"branch_rates\[(\d+)\]\s*\*\s*branch_time\[(\d+)\]", c)
            if mm.group(1) != mm.group(2):
                sys.exit("rate and time indices differ: " + c)
            br.append(int(mm.group(1)) - 1)
        mus[j] = br
    obs = {}
    for m in re.finditer(r"^b_obs\[(\d+)\]\s*~\s*dnNormal\(\s*b_mu\[\d+\]\s*,\s*([0-9.eE+-]+)\s*\)", ftxt, re.M):
        obs[int(m.group(1))] = float(m.group(2))
    clamp = {}
    for m in re.finditer(r"^b_obs\[(\d+)\]\.clamp\(\s*([0-9.eE+-]+)\s*\)", ftxt, re.M):
        clamp[int(m.group(1))] = float(m.group(2))
    jj = sorted(mus)
    if not (set(jj) == set(obs) == set(clamp)):
        sys.exit("terms do not line up")
    for k, j in enumerate(jj):
        terms.append((k, clamp[j], obs[j], mus[j]))
    if sum(len(t[3]) for t in terms) != K:
        sys.exit("terms do not cover every rooted branch exactly once: %d of %d" % (sum(len(t[3]) for t in terms), K))

# ---- UCLN clock (catfish), starting rates and node calibrations ----
ucln = False; init_sigma = 0.5; init_rates = None
ucln_mean_init = None
if a.mode == "fast":
    if "dnLognormal(ucln_mu" in ftxt:
        ucln = True
        m = re.search(r"ucln_mean\.setValue\(\s*([0-9.eE+-]+)", ftxt)
        if m: ucln_mean_init = float(m.group(1))
        m = re.search(r"ucln_sigma\.setValue\(\s*([0-9.eE+-]+)", ftxt)
        if m: init_sigma = float(m.group(1))
        init_rates = [None] * K
        for m in re.finditer(r"^branch_rates\[(\d+)\]\.setValue\(\s*([0-9.eE+-]+)", ftxt, re.M):
            init_rates[int(m.group(1)) - 1] = float(m.group(2))
        if any(v is None for v in init_rates):
            init_rates = None
cals = []
if os.path.exists(d + "/calibrations.Rev"):
    ctxt = rd("calibrations.Rev")
    for m in re.finditer(r"dnUniform\(\s*age_(\w+)\s*-\s*([0-9.eE+-]+)\s*,\s*age_\1\s*-\s*([0-9.eE+-]+)\s*\)", ctxt):
        nm, hi, lo = m.group(1), float(m.group(2)), float(m.group(3))
        if nm not in ext_ids:
            sys.exit("calibration on unknown node " + nm)
        cals.append((ext_ids[nm], lo, hi))
        nodeage = [x for x in nodes if x["ch"] and x["id"] == ext_ids[nm]][0]["age"]
        if not (lo <= nodeage <= hi):
            sys.exit("initial tree violates calibration on " + nm)
# clock init
clock_init = sett["clock_rate_init"]
if clock_init is None and ucln_mean_init is not None:
    clock_init = ucln_mean_init
if clock_init is None:
    # bears: the observed rates are of the order of the summed tree length over path; fall back on 0.01
    clock_init = 0.01
    if a.bears:
        rt = open(a.bears).read()
        m = re.search(r"branch_rates_mean\.setValue\(\s*([0-9.eE+-]+)", rt)
        if m:
            clock_init = float(m.group(1))

# ---- write ----
with open(a.out, "w") as f:
    f.write("FFBD 1\n")
    f.write("cond %s\nmode %s\nn %d\nm %d\nrho %.17g\n" % (sett["cond"], a.mode, len(extant), len(fossils), float(sett["rho"])))
    f.write("root_prior " + " ".join(rp) + "\n")
    f.write("origin_prior " + ("uniform %.17g %.17g" % origin_prior if origin_prior else "none") + "\n")
    f.write("init_lambda %g\ninit_mu %g\ninit_psi %g\ninit_cond_age %.17g\ninit_clock %.17g\n" % (init_lam, init_mu, init_psi, float(init_cond), float(clock_init)))
    if a.mode == "fast" and ucln:
        f.write("rate_prior ucln\ninit_sigma %g\n" % init_sigma)
    f.write("tips %d\n" % N)
    for t in extant + fossils:
        mn, mx = taxa[t]
        f.write("%d %s %s %.17g %.17g\n" % (name2id[t], t, "E" if t in extant else "F", mn, mx))
    f.write("initial_tree %d\n" % len(nodes))
    for x in sorted(nodes, key=lambda z: z["id"]):
        p = x["parent"]["id"] if x["parent"] is not None else -1
        ch = [c["id"] for c in x["ch"]] if x["ch"] else [-1, -1]
        f.write("%d %d %d %d %.17g %d\n" % (x["id"], p, ch[0], ch[1], x["age"], x["sa"]))
    f.write("constraints %d\n" % len(cons))
    for nm, s in cons:
        ids = sorted(name2id[t] for t in s)
        f.write("%s %d %s\n" % (nm, len(ids), " ".join(map(str, ids))))
    f.write("extant_nodes %d\n" % len(ext_nodes))
    for nid, nm in ext_nodes:
        f.write("%d %s\n" % (nid, nm))
    if a.mode == "fast":
        f.write("branches %d\n" % K)
        for i in range(K):
            f.write("%d %d %d\n" % (i, branches[i][0], branches[i][1]))
        f.write("terms %d\n" % len(terms))
        for k, o, se, br in terms:
            f.write("%d %.17g %.17g %d %s\n" % (k, o, se, len(br), " ".join(map(str, br))))
    else:
        f.write("branches 0\nterms 0\n")
    if cals:
        f.write("calibrations %d\n" % len(cals))
        for v, lo, hi in cals:
            f.write("%d %.17g %.17g\n" % (v, lo, hi))
    if a.mode == "fast" and init_rates is not None:
        f.write("init_rates %d\n%s\n" % (K, " ".join("%.17g" % r for r in init_rates)))
print("wrote", a.out, ": extant", len(extant), "fossils", len(fossils), "constraints", len(cons), "extant nodes", len(ext_nodes),
      "branches", K, "terms", len(terms), "SA in start", sa_count)
