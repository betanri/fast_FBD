#!/usr/bin/env python3
"""MrBayes 3.2.7 full-DNA FBD on A06 and A12: an independent DNA check of ffbd (user request, 12 Nov 2026).

Same as the RevBayes full-DNA runs (rb_e2) wherever MrBayes allows: the 23-gene alignment (16,107 sites), GTR+G4 unpartitioned, the
extant tree fixed (one partial constraint per extant node, fossils free), the clade constraints (hard), uniform fossil windows
(calibrate), rho, random sampling, and the same uniform root prior (treeagepr). Differences that MrBayes forces:
  - clock: IGR (independent gamma rates; MrBayes has no exponential-per-branch clock), clockratepr exponential(10), igrvarpr exponential(10)
  - FBD priors on MrBayes's parameters: net diversification exponential(10), turnover beta(1,1), sampling proportion beta(1,1)
So a prior-only run (data=no, same priors) is written too: it shows how much of each age comes from the DNA.

  make_mb_dna.py OUTDIR      (run from FBD/tree)  -> OUTDIR/{SET}_mbdna.nex, {SET}_mbprior.nex  (SEED and OUTPREFIX to fill in)
"""
import csv, os, re, sys


def newick_clades(s):
    """extant clades (frozensets of tip names) of a newick string"""
    s = re.sub(r':[-0-9.eE+]+', '', s.strip().rstrip(';'))
    stack, out = [[]], []
    tok = re.findall(r'\(|\)|,|[^(),]+', s)
    for t in tok:
        if t == '(':
            stack.append([])
        elif t == ')':
            c = stack.pop()
            fs = frozenset(x for y in c for x in (y if isinstance(y, frozenset) else [y]))
            out.append(fs)
            stack[-1].append(fs)
        elif t != ',':
            stack[-1].append(t.strip())
    return out

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
ROOTPR = {'A06': (69.5, 77.5), 'A12': (132.2, 146.2)}     # the uniform root priors of the RevBayes and ffbd runs


def fasta(p):
    d = {}
    for b in open(p).read().split('>')[1:]:
        n, s = b.split('\n', 1)
        d[n.split()[0]] = s.replace('\n', '').upper()
    return d


for S, (lo, hi) in ROOTPR.items():
    d = f'fast_fbd_test/fish_root/data/{S}_root'
    num, den = re.search(r'^rho <- *([0-9.]+) */ *([0-9.]+)', open(f'{d}/settings.Rev').read(), re.M).groups()
    rho = float(num) / float(den)
    aln = fasta(f'{d}/aln.fasta')
    L = {len(v) for v in aln.values()}
    assert len(L) == 1
    L = L.pop()
    taxa = list(csv.DictReader(open(f'{d}/fbd_taxa.tsv'), delimiter='\t'))
    ext = [r['taxon'] for r in taxa if float(r['max_age']) == 0]
    fos = [r for r in taxa if float(r['max_age']) > 0]
    assert set(ext) == set(aln), 'alignment and extant taxa differ'
    cons = [(r['constraint'], r['members'].split(',')) for r in csv.DictReader(open(f'{d}/constraints_members.tsv'), delimiter='\t')]
    bb = [c for c in newick_clades(open(f'{d}/backbone.nwk').read()) if 1 < len(c) < len(ext)]
    assert len(bb) == len(ext) - 2
    for kind in ('dna', 'prior'):
        M = ['#NEXUS', 'begin data;', f'  dimensions ntax={len(taxa)} nchar={L};', '  format datatype=dna missing=? gap=-;', '  matrix']
        M += [f'    {r["taxon"]} {aln[r["taxon"]] if r["taxon"] in aln else "?" * L}' for r in taxa]
        M += ['  ;', 'end;', '', 'begin mrbayes;', '  set autoclose=yes nowarn=yes seed=SEED swapseed=SEED;', '  lset nst=6 rates=gamma ngammacat=4;']
        M += [f'  calibrate {r["taxon"]} = uniform({float(r["min_age"]):g},{float(r["max_age"]):g});' for r in fos]
        names = []
        for nm, mem in cons:
            M.append(f'  constraint {nm} hard = {" ".join(mem)};'); names.append(nm)
        for k, c in enumerate(bb):
            M.append(f'  constraint bb{k + 1} partial = {" ".join(sorted(c))} : {" ".join(sorted(set(ext) - c))};'); names.append(f'bb{k + 1}')
        M += [f'  prset topologypr = constraints({",".join(names)});', '  prset brlenspr = clock:fossilization;',
              '  prset speciationpr = exponential(10);', '  prset extinctionpr = beta(1,1);', '  prset fossilizationpr = beta(1,1);',
              f'  prset sampleprob = {rho:.10g};', '  prset samplestrat = random;',
              f'  prset treeagepr = uniform({lo:g},{hi:g});', '  prset nodeagepr = calibrated;',
              '  prset clockratepr = exponential(10);', '  prset clockvarpr = igr;', '  prset igrvarpr = exponential(10);',
              f'  mcmcp data={"yes" if kind == "dna" else "no"} ngen=NGEN samplefreq=SFREQ printfreq=100000 nruns=1 nchains=1 '
              'diagnfreq=1000000 checkpoint=yes checkfreq=1000000 filename=OUTPREFIX;',
              '  mcmc;', 'end;', '']
        open(f'{OUT}/{S}_mb{kind}.nex', 'w').write('\n'.join(M))
    print(S, 'extant', len(ext), 'fossils', len(fos), 'sites', L, 'hard constraints', len(cons), 'backbone partial constraints', len(bb),
          'rho', round(rho, 6), 'root prior', lo, hi)
