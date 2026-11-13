#!/usr/bin/env python3
"""Convert a folder of RevBayes-style inputs into one ffbd input file.

  ffbd_convert.py --set DIR --out FILE [--mode fast|control] [--root-prior normal MEAN SD | uniform MIN MAX | soft MIN MAX SD P]
                  [--auto-start] [--write-start FILE]
                  [--bears REVFILE] [--seqrev DIR]

DIR holds: fbd_taxa.tsv (or bears_taxa_specimen.tsv), settings.Rev, backbone.nwk, constraints_members.tsv or constraints.Rev,
extant_nodes.Rev, fast_terms.Rev; optional: initial.nwk (a starting tree with fossils; built automatically when absent or with
--auto-start, see auto_start), crown.tsv, stem.tsv, calibrations.Rev.  --bears gives the Rev script of the bear runs (origin prior, constraint).
Without --root-prior the prior comes from settings.Rev (root_mean, root_sd: Normal).  Only reads files; never edits the data.
"""
import argparse, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("--set", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--mode", default="fast"); ap.add_argument("--root-prior", nargs="+", default=None)
ap.add_argument("--bears", default=None)
ap.add_argument("--terms-file", default="fast_terms.Rev")
ap.add_argument("--auto-start", action="store_true", help="build the starting tree even if initial.nwk exists")
ap.add_argument("--write-start", default=None, help="also save the automatic starting tree (newick) to this file")
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

# ---- starting tree: initial.nwk if the folder has one (and --auto-start is not given), else built automatically ----
def _post(root):
    out, st = [], [root]
    while st:
        n = st.pop()
        out.append(n)
        st.extend(n["ch"])
    return out[::-1]


def extant_defs():
    """extant clade name (as in extant_nodes.Rev, without 'clade_') -> set of tip names"""
    ent_ = (rd("extant_clades.Rev") + "\n" if os.path.exists(d + "/extant_clades.Rev") else "") + rd("extant_nodes.Rev")
    dd = {}
    for m in re.finditer(r"^(clade_\w+)\s*=\s*clade\((.*?)\)\s*$", ent_, re.M):
        st_ = set()
        for arg in re.split(r",\s*", m.group(2)):
            arg = arg.strip()
            if arg.startswith('"'):
                st_.add(arg.strip('"'))
            elif arg in dd:
                st_ |= dd[arg]
            else:
                sys.exit("unknown reference in extant_nodes.Rev: " + arg)
        dd[m.group(1)] = st_
    return {k[6:]: v for k, v in dd.items()}


def auto_start():
    """A valid starting tree from the tree of living species (backbone.nwk), the fossil windows and the clade constraints.
    Start ages: mean path height of the backbone branch lengths (edge counts if it has none), scaled to the initial root age.
    Each fossil goes on the stem of the smallest constraint that contains it (the node of its living members), at its minimum
    age; fossils on the same stem are stacked from the smallest constraint up. A crown-rule fossil goes inside the crown, on
    the stem of the crown node's first child. Nodes are then pushed up wherever needed: above their children and fossils,
    and to the minimum of any node calibration. All checks below then run on this tree as on a supplied one."""
    EPS = 0.01
    bb = parse(rd("backbone.nwk"))
    bb["parent"] = None
    order = _post(bb)
    for x in order:
        for c in x["ch"]:
            c["parent"] = x
        x["tips"] = frozenset([x["name"]]) if not x["ch"] else frozenset().union(*[c["tips"] for c in x["ch"]])
        if x["ch"] and len(x["ch"]) != 2:
            sys.exit("auto start: backbone.nwk must be fully bifurcating")
    if bb["tips"] != frozenset(extant):
        sys.exit("auto start: backbone.nwk tips differ from the living species of the taxa file")
    exset_ = set(extant)
    byset = {x["tips"]: x for x in order}
    r0 = init_cond
    if r0 is None:
        r0 = {"normal": lambda: float(rp[1]), "uniform": lambda: 0.5 * (float(rp[1]) + float(rp[2])),
              "soft": lambda: 0.5 * (float(rp[1]) + float(rp[2]))}.get(rp[0], lambda: None)()
    if r0 is None:
        sys.exit("auto start: give root_age_init in settings.Rev or a root prior")
    uselen = any(c["len"] > 0 for x in order for c in x["ch"])
    h, w = {}, {}
    for x in order:
        if not x["ch"]:
            h[id(x)], w[id(x)] = 0.0, 1
        else:
            w[id(x)] = sum(w[id(c)] for c in x["ch"])
            h[id(x)] = sum(w[id(c)] * (h[id(c)] + (c["len"] if uselen else 1.0)) for c in x["ch"]) / w[id(x)]
    k_ = r0 / h[id(bb)]

    def node_of(E, what):
        x = byset.get(frozenset(E))
        if x is None:
            sys.exit("auto start: the living species of %s are not a clade of backbone.nwk" % what)
        return x

    def below(x, v):                                   # v is a proper ancestor of x
        u = x["parent"]
        while u is not None:
            if u is v:
                return True
            u = u["parent"]
        return False
    defs_ = None
    rules = {}
    for kind in ("crown", "stem"):
        if os.path.exists(d + "/%s.tsv" % kind):
            defs_ = defs_ or extant_defs()
            for ln in rd("%s.tsv" % kind).splitlines():
                if ln.strip() and not ln.startswith("#"):
                    fname, nm = ln.split("\t")[:2]
                    if nm not in defs_:
                        sys.exit("%s entry unknown: %s" % (kind, ln))
                    rules[fname] = (kind, node_of(defs_[nm] & exset_, "%s node %s" % (kind, nm)))
    cal_lo, cal_hi = {}, {}
    if os.path.exists(d + "/calibrations.Rev"):
        defs_ = defs_ or extant_defs()
        for m in re.finditer(r"dnUniform\(\s*age_(\w+)\s*-\s*([0-9.eE+-]+)\s*,\s*age_\1\s*-\s*([0-9.eE+-]+)\s*\)", rd("calibrations.Rev")):
            x = node_of(defs_[m.group(1)] & exset_, "calibrated node " + m.group(1))
            cal_lo[id(x)] = max(cal_lo.get(id(x), 0.0), float(m.group(3)))
            cal_hi[id(x)] = min(cal_hi.get(id(x), float("inf")), float(m.group(2)))
    lad = {id(x): [] for x in order}
    crown_on = {id(x): [] for x in order}
    for f in fossils:
        homes = [(len(st_), nm) for nm, st_ in cons if f in st_]
        if not homes:
            sys.exit("auto start: fossil %s is in no clade constraint (every fossil needs one below the root)" % f)
        size, nm = min(homes)
        E = dict(cons)[nm] & exset_
        if not E:
            sys.exit("auto start: constraint %s has no living species; give initial.nwk" % nm)
        n_ = node_of(E, "constraint " + nm)
        kind, v = rules.get(f, (None, None))
        if kind == "crown" and (n_ is v or below(v, n_)):
            if not v["ch"]:
                sys.exit("auto start: crown rule of %s is on a single species" % f)
            crown_on[id(v["ch"][0])].append((taxa[f][0], f))
            continue
        if kind == "stem" and below(n_, v):
            sys.exit("auto start: stem rule of %s conflicts with its constraint %s" % (f, nm))
        lad[id(n_)].append((size, taxa[f][0], f))
    if lad[id(bb)]:
        sys.exit("auto start: fossils would sit above the root: %s" % [f for _, _, f in lad[id(bb)]][:3])
    top, ladder, age = {}, {}, {}
    for x in order:
        a = 0.0 if not x["ch"] else max([h[id(x)] * k_, cal_lo.get(id(x), 0.0)] + [top[id(c)] + EPS for c in x["ch"]])
        if a > cal_hi.get(id(x), float("inf")):
            sys.exit("auto start: a calibrated node had to be older than its maximum (%g > %g)" % (a, cal_hi[id(x)]))
        age[id(x)] = a
        L = []
        for _, lo, f in sorted(lad[id(x)]):
            a = max(a, lo) + EPS
            L.append((a, f))
        for lo, f in sorted(crown_on[id(x)]):
            a = max(a, lo) + EPS
            L.append((a, f))
        top[id(x)], ladder[id(x)] = a, L
    if rp[0] == "uniform" and age[id(bb)] > float(rp[2]):
        sys.exit("auto start: the root had to be %.2f, above the root prior maximum %s" % (age[id(bb)], rp[2]))
    sys.setrecursionlimit(max(10000, 4 * len(taxa)))

    def nw(x):
        if not x["ch"]:
            s_, a = x["name"], 0.0
        else:
            parts = []
            for c in x["ch"]:
                cs, ca = nw(c)
                parts.append("%s:%r" % (cs, age[id(x)] - ca))
            s_, a = "(" + ",".join(parts) + ")", age[id(x)]
        for a2, f in ladder[id(x)]:
            s_ = "(%s:%r,%s:%r)" % (s_, a2 - a, f, a2 - taxa[f][0])
            a = a2
        return s_, a
    return nw(bb)[0] + ";", age[id(bb)]


auto = a.auto_start or not os.path.exists(d + "/initial.nwk")
if auto:
    start_nwk, init_cond = auto_start()
    print("starting tree built automatically (root %.2f Ma)" % init_cond)
    if a.write_start:
        open(a.write_start, "w").write(start_nwk + "\n")
init = parse(start_nwk if auto else rd("initial.nwk"))
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

# tip sets bottom-up once (the old recursive version recomputed every subtree for every node: hours on 25k tips)
for x in _post([n for n in nodes if n["parent"] is None][0]):
    x["tips"] = frozenset([x["name"]]) if not x["ch"] else frozenset().union(*[c["tips"] for c in x["ch"]])

# check the initial tree satisfies the constraints
tsets = {x["tips"] for x in nodes}
for name, s in cons:
    if frozenset(s) not in tsets:
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
# MRCA lookup: the MRCA of an extant clade is the smallest internal node whose extant tips are exactly that clade
by_ext = {}
for x in nodes:
    if x["ch"]:
        k = frozenset(x["tips"] & exset)
        if k not in by_ext or len(x["tips"]) < len(by_ext[k]["tips"]):
            by_ext[k] = x
for nm in ordered:
    if "clade_" + nm not in defs:
        continue
    s = defs["clade_" + nm]
    # MRCA of the clade's extant tips, restricted to extant taxa (fossil constraint nodes are not extant nodes)
    target = s & exset
    best = by_ext.get(frozenset(target))
    if best is None:  # not an exact extant clade: fall back to the full search
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

# ---- crown option: crown.tsv (fossil<TAB>extant node name, e.g. n123): the fossil may attach only below that extant node ----
crown, stem = [], []
byname = {x["name"]: x for x in nodes if not x["ch"]}
for kind, dst in (("crown", crown), ("stem", stem)):
    if not os.path.exists(d + "/%s.tsv" % kind):
        continue
    for ln in rd("%s.tsv" % kind).splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        fname, nm = ln.split("\t")[:2]
        if fname not in name2id or nm not in ext_ids:
            sys.exit("%s entry unknown: %s" % (kind, ln))
        v = ext_ids[nm]
        u = byname[fname]["parent"]
        while u is not None and u["id"] != v:
            u = u["parent"]
        if (u is None) == (kind == "crown"):
            sys.exit("initial tree breaks the %s rule of %s" % (kind, fname))
        dst.append((name2id[fname], v))

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
    for kind, lst in (("crown", crown), ("stem", stem)):
        if lst:
            f.write("%s %d\n" % (kind, len(lst)))
            for fid, v in lst:
                f.write("%d %d\n" % (fid, v))
print("wrote", a.out, ": extant", len(extant), "fossils", len(fossils), "constraints", len(cons), "extant nodes", len(ext_nodes),
      "branches", K, "terms", len(terms), "SA in start", sa_count)
