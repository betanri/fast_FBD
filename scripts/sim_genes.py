#!/usr/bin/env python3
"""Few genes against many genes, with known true ages (test 5e of VALIDATION.md).

Reuses the true trees of the simulation study (sim_study.py simulate: sim1/rN/phylo_subs.nwk = true rate x time per branch).
Simulates a multi-gene data set on each true tree, gene by gene, so genes differ the way real genes do:
  - nuclear genes: own rate multiplier ~ Gamma(shape 2, mean 1), own GTR exchangeabilities and base frequencies, own gamma
    shape alpha ~ Uniform(0.3, 1.5); length Uniform(500, 1000) sites for the 21 "exon" genes, Uniform(300, 700) for the rest
  - mitochondrial-like genes: rate multiplier 10 (they saturate on deep branches), AT-rich composition, alpha 0.3
Gene sets (the few set is part of the many set, as in the real data):
  few  = 21 nuclear exons + 2 mitochondrial-like genes (1,500 and 1,100 sites)     ~ 18,000 sites (like the 23-gene set)
  many = few + 8 more mitochondrial-like genes + 969 more nuclear genes = 1,000 genes  ~ 500,000 sites (like the all-gene set)
Both are analyzed the same way as the real clades: concatenated, one GTR+G4 model, ML branch lengths on the true topology.

  sim_genes.py REPDIR SEED       writes REPDIR/genes.tsv and REPDIR/gene_trees/gNNNN.nwk (one scaled true phylogram per gene),
                                 then the sbatch script runs AliSim per gene and concatenates few.fa and many.fa
"""
import os, random, re, sys

D, seed = sys.argv[1], int(sys.argv[2])
rng = random.Random(seed)
os.makedirs(f'{D}/gene_trees', exist_ok=True)
tree = open(f'{D}/phylo_subs.nwk').read().strip()


def scaled(t, m):
    return re.sub(r':([0-9.eE+-]+)', lambda x: f':{float(x.group(1)) * m:.10f}', t)


def nuclear(length):
    rates = [round(rng.lognormvariate(0, 0.3) * r, 4) for r in (1.0, 3.0, 1.0, 1.0, 3.0, 1.0)]
    f = [rng.gammavariate(20, 1) for _ in range(4)]
    f = [round(x / sum(f), 4) for x in f]
    f[3] = round(1 - sum(f[:3]), 4)
    return dict(kind='nuclear', length=length, mult=rng.gammavariate(2.0, 0.5), alpha=round(rng.uniform(0.3, 1.5), 3),
                model=f'GTR{{{",".join(map(str, rates[:5]))}}}+F{{{",".join(map(str, f))}}}')


def mito(length):
    return dict(kind='mito', length=length, mult=10.0, alpha=0.3, model='GTR{1.0,8.0,1.0,1.0,8.0}+F{0.32,0.28,0.12,0.28}')


genes = [nuclear(rng.randint(500, 1000)) for _ in range(21)] + [mito(1500), mito(1100)]          # few set (23 genes)
genes += [mito(rng.randint(300, 1800)) for _ in range(8)] + [nuclear(rng.randint(300, 700)) for _ in range(969)]   # many set adds
with open(f'{D}/genes.tsv', 'w') as fh:
    fh.write('gene\tset\tkind\tlength\trate_multiplier\talpha\tmodel\n')
    for i, g in enumerate(genes, 1):
        fh.write(f'g{i:04d}\t{"few" if i <= 23 else "many"}\t{g["kind"]}\t{g["length"]}\t{g["mult"]:.4f}\t{g["alpha"]}\t'
                 f'{g["model"]}+G4{{{g["alpha"]}}}\n')
        open(f'{D}/gene_trees/g{i:04d}.nwk', 'w').write(scaled(tree, g['mult']) + '\n')
few = sum(g['length'] for g in genes[:23]); many = sum(g['length'] for g in genes)
print(f'{D}: 23 few genes {few} sites; 1,000 many genes {many} sites')
