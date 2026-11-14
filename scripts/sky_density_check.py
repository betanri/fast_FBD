#!/usr/bin/env python3
"""FBD log density of ffbd with rate shifts (skyline) against RevBayes 1.3.2 dnBDSTP(timeline=...).

  sky_density_check.py gen INPUT DUMPFILE OUTDIR SHIFTS [--reverse] [--width W]
      INPUT: the ffbd input file (tip types, rho, cond); DUMPFILE: from ffbd --dump (lambda, mu, psi fields are comma lists, one
      value per interval, youngest first); SHIFTS: comma list of shift ages, or "none". Writes OUTDIR/state_K.Rev (+ tree, taxa).
      --reverse passes the rate vectors oldest first (only to test the RevBayes ordering convention).
  sky_density_check.py compare DUMPFILE RBOUT [SHIFTS]
      RBOUT: stdout of the rb runs, lines "DENS<TAB>k<TAB>lnP". Prints the largest |ffbd - RevBayes| of the FBD density.
The density compared is ffbd's lnFBDFull (dnBDSTP lnProbability: process density and the tree-shape term).
"""
import sys, pathlib, re, math

def read_dump(p):
    out = []
    for line in open(p):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 15:
            continue
        out.append({"lam": parts[1].split(","), "mu": parts[2].split(","), "psi": parts[3].split(","), "cond": parts[4],
                    "fbd": float(parts[7]), "nsa": int(parts[13]), "nwk": parts[14]})
    return out

def read_input(p):
    rho, cond, types = None, None, {}
    toks = open(p).read().split("\n")
    i = 0
    while i < len(toks):
        w = toks[i].split()
        if w and w[0] == "rho": rho = w[1]
        if w and w[0] == "cond": cond = w[1]
        if w and w[0] == "tips":
            n = int(w[1])
            for j in range(n):
                q = toks[i + 1 + j].split(); types[q[1]] = q[2]
            break
        i += 1
    return rho, cond, types

def tip_ages(nwk):
    # depths from the root; age = max extant depth - depth (extant tips at age 0)
    s = nwk.strip().rstrip(";")
    depth = {}
    pos = [0]
    def parse(d0):
        if s[pos[0]] == "(":
            pos[0] += 1
            kids = []
            while True:
                kids.append(parse_child(d0))
                if s[pos[0]] == ",": pos[0] += 1; continue
                if s[pos[0]] == ")": pos[0] += 1; break
            return kids
        m = re.match(r"[^:,()]+", s[pos[0]:]); nm = m.group(0); pos[0] += len(nm)
        return nm
    def parse_child(d0):
        node = parse(None)
        bl = 0.0
        if pos[0] < len(s) and s[pos[0]] == ":":
            m = re.match(r":([-0-9.eE+]+)", s[pos[0]:]); bl = float(m.group(1)); pos[0] += len(m.group(0))
        return (node, bl)
    root = parse(None)
    def walk(node, d):
        if isinstance(node, str):
            depth[node] = d; return
        for ch, bl in node: walk(ch, d + bl)
    walk(root, 0.0)
    return depth

def node_ages(nwk):
    """ages of all nodes and tips: deepest tip depth (an extant tip) minus depth"""
    seq = [m.group(0) for m in re.finditer(r"\(|\)|,|:[-0-9.eE+]+|[^:,()]+", nwk.strip().rstrip(";"))]
    sys.setrecursionlimit(100000)
    def rec(i):  # returns (next index, children or None for a tip, branch length)
        kids = None
        if seq[i] == "(":
            i += 1; kids = []
            while True:
                i, sub, bl = rec(i); kids.append((sub, bl))
                i += 1
                if seq[i - 1] == ")": break
        else:
            i += 1
        bl = 0.0
        if i < len(seq) and seq[i].startswith(":"): bl = float(seq[i][1:]); i += 1
        return i, kids, bl
    _, tree, _ = rec(0)
    depths = []
    def walk(kids, dd):
        depths.append(dd)
        if kids is not None:
            for sub, bl in kids: walk(sub, dd + bl)
    walk(tree, 0.0)
    H = max(depths)
    return [H - dd for dd in depths]

REV = r'''setOption("outputPrecision", "15", FALSE)
taxa <- readTaxonData("@TAXA@")
tr <- readTrees("@TREE@", treetype="clock")[1]
lam <- @LAM@
mu <- @MU@
psi <- @PSI@
fbd ~ dnBDSTP(@COND@=@CONDV@, lambda=lam, mu=mu, psi=psi, rho=@RHO@, r=0@TL@, taxa=taxa)
fbd.setValue(tr)
print("DENS", @K@, fbd.lnProbability(), fbd.numSampledAncestors(), sep="\t")
'''

def vec(xs):
    return xs[0] if len(xs) == 1 else "v(" + ",".join(xs) + ")"

if __name__ == "__main__":
    if sys.argv[1] == "gen":
        inp, dump, outd, shifts = sys.argv[2:6]
        rev = "--reverse" in sys.argv
        W = float(sys.argv[sys.argv.index("--width") + 1]) if "--width" in sys.argv else 1e-7  # half-width of the fossil taxon age ranges
        rho, cond, types = read_input(inp)
        rows = read_dump(dump)
        o = pathlib.Path(outd); o.mkdir(parents=True, exist_ok=True)
        T = [] if shifts == "none" else shifts.split(",")
        allrev = []
        for k, r in enumerate(rows):
            d = tip_ages(r["nwk"])
            ext = [d[t] for t in d if types[t] == "E"]
            h = max(ext)
            with open(o / f"taxa_{k}.tsv", "w") as f:
                f.write("taxon\tmin_age\tmax_age\n")
                for t in d:
                    a = max(0.0, h - d[t])
                    if types[t] == "E": f.write(f"{t}\t0.0\t0.0\n")
                    else: f.write(f"{t}\t{max(0.0, a - W):.15g}\t{a + W:.15g}\n")
            (o / f"tree_{k}.nwk").write_text(r["nwk"] + "\n")
            L, M, P = r["lam"], r["mu"], r["psi"]
            if not T: L, M, P = L[:1], M[:1], P[:1]   # no timeline: constant-rate dnBDSTP (states must have equal rates in all intervals)
            if not T and len(set(r["lam"])) + len(set(r["mu"])) + len(set(r["psi"])) > 3: sys.exit("SHIFTS none needs equal rates in all intervals")
            if rev: L, M, P = L[::-1], M[::-1], P[::-1]
            txt = (REV.replace("@TAXA@", f"taxa_{k}.tsv").replace("@TREE@", f"tree_{k}.nwk")  # rb is run from inside OUTDIR
                   .replace("@LAM@", vec(L)).replace("@MU@", vec(M)).replace("@PSI@", vec(P))
                   .replace("@COND@", "originAge" if cond == "origin" else "rootAge")
                   .replace("@CONDV@", r["cond"] if cond == "origin" else "tr.rootAge()")
                   .replace("@RHO@", rho).replace("@TL@", (", timeline=v(" + ",".join(T) + ")") if T else "").replace("@K@", str(k)))
            allrev.append(txt)
        # one script per state (an rb error in one state does not stop the others)
        for k, txt in enumerate(allrev):
            (o / f"state_{k}.Rev").write_text(txt + "q()\n")
        print("wrote", len(rows), "states to", o)
    else:
        dump, rbout = sys.argv[2:4]
        rows = read_dump(dump)
        res = {}
        for line in open(rbout, errors="replace"):
            if line.startswith("DENS"):
                p = line.split(); res[int(p[1])] = (float(p[2]), int(float(p[3])))
        # RevBayes 1.3.2 computes ln N! in the tree-shape term by Stirling's series up to the 1/(360 N^3) term (error 4.7e-8 at N = 7,
        # 1.5e-10 at N = 22, below 1e-12 from N = 150); "corrected" removes that known constant from the difference.
        shifts = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 and sys.argv[4] != "none" else []
        N = len(tip_ages(rows[0]["nwk"])) if rows else 0
        st = N * math.log(N) - N + 0.5 * math.log(2 * math.pi * N) + 1 / (12 * N) - 1 / (360 * N ** 3) if N else 0.0
        corr = st - math.lgamma(N + 1) if N else 0.0
        worst = {"regular": 0.0, "at_shift": 0.0, "degenerate": 0.0}; cnt = {k: 0 for k in worst}; nr = 0; bad = 0; maxabs = 0.0
        for k, r in enumerate(rows):
            if k not in res: print(k, "MISSING RevBayes result"); continue
            nr += 1; dd = r["fbd"] - res[k][0] - corr; maxabs = max(maxabs, abs(r["fbd"]))
            # degenerate: a non-SA fossil tip with branch length exactly 0 (its age equals its parent's, by floating-point rounding)
            if re.search(r"[(,][^:,()]+:0(?![.0-9eE])", r["nwk"]): cls = "degenerate"
            elif shifts and any(abs(a - T) < 1e-9 for a in node_ages(r["nwk"]) for T in shifts): cls = "at_shift"
            else: cls = "regular"
            cnt[cls] += 1
            if abs(dd) > abs(worst[cls]): worst[cls] = dd
            if r["nsa"] != res[k][1]: bad += 1
        nsas = sorted(set(r["nsa"] for r in rows))
        extra = "".join(f"; {c} states {cnt[c]}: worst {worst[c]:.3g}" for c in ("at_shift", "degenerate") if cnt[c])
        print(f"states {nr}/{len(rows)}  N {N}  SA {min(nsas)}-{max(nsas)} (counts differ in {bad})  max|lnP| {maxabs:.4g}  RevBayes ln N! error {corr:.3g}  |  corrected worst ffbd-RevBayes: regular states {cnt['regular']}: {worst['regular']:.3g}{extra}")
