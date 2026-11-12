# fast_FBD

A fast Bayesian program for dating phylogenies with fossils, using the **fossilized birth-death (FBD)** process.

> **Status: research software.** The math and sampling have passed strict tests, and Fast FBD matches full-DNA FBD on two test clades. See [Validation](#validation) before using results.

## What it does

FBD dating treats fossils as real tips on the tree, instead of turning them into minimum ages on nodes. That uses all the fossils, not just the oldest ones, but standard programs become extremely slow on large trees.

`ffbd` makes this practical with two ideas:

1. **Fast FBD.** The DNA alignment is summarized once, as a tree of branch lengths with a standard error for each branch. The program then asks how well each proposed set of dates and rates explains those branch lengths, instead of recalculating the full DNA likelihood.
2. **Local updates.** When one age changes, only the affected terms are recalculated. Cost per step does not grow with tree size.

On one CPU core it runs about **2.6 million MCMC proposals per second**, with no slowdown between 1,185 and 2,424 tips.

## What is fixed and what is estimated

| Fixed | Estimated |
|---|---|
| The relationships among living species | Ages of all nodes |
| Sampling fraction of living species (rho) | Speciation, extinction and fossil sampling rates |
| Fossil clade constraints and age windows | Where each fossil attaches, within its constraint |
| | Each fossil's age, within its window |
| | Whether a fossil is a sampled ancestor |
| | Rates of molecular evolution, per branch |

**Use it when** the tree of living species is settled, you have many fossils, and full-DNA FBD is too slow.

**Do not use it when** the topology of living species is in question, or for very short or poorly supported branches where a Normal approximation of branch length is doubtful.

## Install

You need a C++17 compiler and nothing else.

```bash
git clone https://github.com/betanri/fast_FBD.git
cd fast_FBD
make
```

## Quick start

A tiny simulated data set is included (4 living species, 3 fossils):

```bash
mkdir -p runs
./ffbd run examples/example_input.txt --prefix runs/example --seed 1 --gens 2000000 --print 2000 --tree-print 20000 --burnin 100000
```

This takes under a second and writes:

| File | Contents |
|---|---|
| `runs/example.log` | One row per sample: posterior, rates, and the age of each node (`age_n1`, ...), readable in Tracer |
| `runs/example.trees` | Sampled trees with fossils as tips |

Run at least two chains with different `--seed` values and check that they agree.

## What you need

| Input | Notes |
|---|---|
| Rooted tree of living species | Fixed topology |
| Branch lengths and their standard errors | For example ML branch lengths on the fixed tree, with bootstrap standard errors |
| Fossil list | Name, minimum and maximum age |
| Clade constraints | Which fossils belong in which clades |
| Crown or stem rules (optional) | Fossils that must sit inside a crown, or on its stem only (see below) |
| Sampling fraction (rho) | Share of living species included |
| Root age prior | `uniform MIN MAX`, `soft MIN MAX SD P` (soft bounds) or `normal MEAN SD` |
| A valid starting tree | With fossils attached, and consistent with the constraints |

These go in one text file. `examples/example_input.txt` shows every section: the settings, then `tips`, `initial_tree`, `constraints`, `extant_nodes`, `branches` and `terms`. If your inputs are already in RevBayes-style files, `scripts/ffbd_convert.py` builds this file for you.

### Crown and stem rules

A clade constraint says a fossil belongs to a group, and that includes the group's stem. Two optional sections narrow this:

- `crown K` followed by K lines `fossil node`: the fossil may attach only below that extant node, inside the crown. The crown node must then be at least as old as the fossil.
- `stem K` followed by K lines `fossil node`: the fossil may not attach below that extant node. Its clade constraint keeps it in the total group, so it sits on the stem.

`fossil` is the tip index and `node` the extant node index used in the input file. The converter reads them from `crown.tsv` and `stem.tsv` (fossil name, extant node name) in the input folder. The rules only restrict where a fossil can attach; they add no age prior. `examples/example_crown_stem.txt` shows both.

## Main options

| Option | Meaning |
|---|---|
| `--prefix P` | Start of the output file names |
| `--seed S` | Random seed |
| `--gens G` | Number of proposals |
| `--burnin N` | Proposals used to tune moves, not logged |
| `--print N`, `--tree-print N` | How often to log parameters and trees |
| `--checkpoint N`, `--resume` | Save progress and continue a run exactly where it stopped |
| `--mode control` | Run without molecular data, to see what the fossils and priors alone imply |
| `--log-rates` | Also log branch rates |
| `--check` | Check internal consistency while running (slow) |

## The model

- Constant-rate fossilized birth-death process with sampled ancestors, conditioned on the root age. The density matches RevBayes 1.3.2 `dnBDSTP` inside `dnConstrainedTopology`.
- Fossil ages are uniform within their windows.
- Clock: one rate per branch of the living-species tree, exponential by default or uncorrelated lognormal.
- Data: for each branch, observed length ~ Normal(rate x time, standard error). The two branches at the root share one observation.

## Validation

Full details and result files are in [VALIDATION.md](VALIDATION.md) and `validation/`. Status as of 11 November 2026.

| Test | Result |
|---|---|
| Probabilities equal RevBayes 1.3.2 on fixed trees | **Passed**: difference at most 4e-9 on 89 states across 4 data sets |
| Simulation-based calibration, data simulated from known truth (500 replicates, 3 settings) | **Passed** |
| Prior-only run against exact analytic results | **Passed** |
| Internal consistency and exact restart from a checkpoint | **Passed** |
| Agreement with RevBayes with sampled-ancestor changes switched off | **Passed**: node ages differ by 0.3 to 0.4% (median) |
| Agreement with RevBayes with sampled-ancestor changes switched on | **Does not match**, apparently because of a RevBayes move; see below |
| Fast FBD against full-DNA FBD (RevBayes, that move removed), two fish clades | **Preliminary, agrees**: node ages r = 0.998 and 0.999, mean difference +0.5 and -0.3 Myr; root ages 75.2 vs 75.1 Ma and 139.0 vs 139.1 Ma |
| Catfish data set (713 living tips, 94 fossils) against the original full analysis | **Close**: node ages r = 0.995, Fast 1.8 Myr younger on average. The original used a lognormal clock and the RevBayes move above, so some difference is expected |
| Convergence on large trees (1,185 and 2,424 tips) | **Passed**: after 20 billion proposals, main-parameter ESS 500 to 800; two chains agree on node ages to 0.9% (median) |
| Three programs on the same fixed-rate prior: `ffbd`, MrBayes 3.2.7, RevBayes 1.3.2 (new) | **`ffbd` = MrBayes**: node ages differ by 0.4 to 0.6% (median), sampled ancestors 7.01 vs 7.02 and 9.55 vs 9.58. RevBayes differs (3 to 7%, about one extra sampled ancestor) |
| Simulation study, 100 trees of 30 to 120 tips, branch lengths from simulated DNA (new) | **Ages close to nominal**: 95% coverage 0.91 for node ages, 0.95 for fossil ages; nodes under 5 Ma about 1 Myr too old |
| Simulation-based calibration at realistic size: 500 replicates, up to 8,056 tips (new) | **Passed**: all 9 quantities uniform (lowest p = 0.087) |
| Crown and stem rules (new) | **Passed**: an exact conditioning test (below); runs without these rules are byte-identical to the previous release |

![Node ages, full-DNA FBD against Fast FBD](validation/node_ages_full_vs_fast.png)

**Preliminary: full-DNA chains still too short; more chains in progress.** Convergence of each side (two chains each, first 10% discarded):

| Data set | Side | Run time per chain | Node-age ESS (median, lowest) | Root-age ESS | Two chains differ (median, largest) |
|---|---|---|---|---|---|
| Gobiaria | Full DNA (RevBayes) | ~25 h | 19, 10 | 51-52 | 1.6%, 6.6% |
| Gobiaria | Fast FBD (`ffbd`) | ~2.5 min | 901, 756 | 764-812 | 0.26%, 1.1% |
| Protacanthopterygii | Full DNA (RevBayes) | ~22 h | 35, 14 | 557-1,531 | 4.7%, 15.3% |
| Protacanthopterygii | Fast FBD (`ffbd`) | ~2.5 min | 901, 778 | 901 | 0.79%, 3.7% |

The full-DNA ESS values are low, so the agreement shown below is within the noise of the full-DNA runs. Four more full-DNA chains per clade are running; this test will be updated when they finish.

| Data set | Living species | Fossils | DNA |
|---|---|---|---|
| Gobiaria | 24 | 10 | 16,107 sites |
| Protacanthopterygii | 22 | 14 | 16,107 sites |
| Catfishes | 713 | 94 | 21 genes, 4,614 sites |

*Each point is one node of the extant time tree (fossils pruned; median ages), matched by clade. Gobiaria and Protacanthopterygii: full-DNA FBD in RevBayes 1.3.2 with `mvCollapseExpandFossilBranch` removed, against `ffbd` on the same data, priors and fossils. Catfishes: the original full analysis (21-gene alignment, partitioned GTR+Gamma, uncorrelated lognormal clock; run with the original RevBayes moves) against `ffbd`. "Full FBD" means the full DNA sequence likelihood; Fast FBD (`ffbd`) uses only branch lengths and their standard errors, estimated once from the same DNA. Dashed line: equal ages.*

**Crown and stem test.** On the small example (4 living species, 3 fossils), FOSSIL_1 was given a crown rule and FOSSIL_3 a stem rule on the same node. The run with these rules was compared with unconstrained runs kept only where both rules happened to hold (10.3% of samples). If the rules are implemented correctly, the two must give the same distribution. With 11 seeds per side, the rules held in all 1,100,011 samples. Every quantity agreed within 1.4 standard errors: speciation, extinction, fossil sampling, number of sampled ancestors and all three node ages (`validation/crown_stem_test_11seeds.txt`).

## Bottom line: ffbd, MrBayes, RevBayes and BEAST2

- **`ffbd` and MrBayes agree.** Two independently written programs, given the same tree, fossils and fixed rates, give the same node ages (within about 0.5%) and the same number of sampled ancestors (7.01 vs 7.02; 9.55 vs 9.58).
- **RevBayes 1.3.2 disagrees with both.** It finds about one extra sampled ancestor and makes young nodes older, by up to a third on some two-species nodes. With one RevBayes move switched off (`mvCollapseExpandFossilBranch`), RevBayes agrees with `ffbd` (0.3 to 0.4%). Our reading is that this move is not reversible when fossils have age ranges. The RevBayes developers have not yet confirmed this.
- **BEAST2 was not run.** `ffbd` keeps the tree of living species fixed while each fossil moves freely within its clade. BEAST2's sampled-ancestor package cannot express that, because its clade constraints count fossils as clade members. A BEAST2 run would therefore be a different model, not a check of this one, so MrBayes was used as the independent program.
- **`ffbd` recovers known truth.** On 500 simulated data sets with trees of up to 8,056 living species, the true values sit where they should in `ffbd`'s estimates. On 100 simulated DNA data sets, node and fossil ages have close to 95% coverage.
- **Practical meaning.** For FBD analyses with fossil age ranges and sampled ancestors, `ffbd` and MrBayes give the same answer. RevBayes 1.3.2 with its default moves can give more sampled ancestors and older young nodes.

### A difference from RevBayes

`ffbd` and RevBayes 1.3.2 give different posteriors when fossils have age ranges and sampled ancestors can be added or removed. The evidence points to the RevBayes move `mvCollapseExpandFossilBranch`. In the 1.3.2 source, it draws the new attachment age above the *maximum* of the fossil's age window, while the reverse move can start from attachments inside the window. Those states can be left but not re-entered in the same way, so the move is not reversible.

With that move removed from both programs, they agree. An independent program, MrBayes 3.2.7, agrees with `ffbd` and not with RevBayes on the same fixed-rate prior (VALIDATION.md, section 9). `ffbd` also recovers known truth in simulations. We have not yet confirmed this with the RevBayes developers. `--rb-sa-move` imitates that proposal, for testing only.

## Known limits

- The topology of living species cannot change.
- Fossils must fall inside a clade constraint below the root in root-conditioned runs.
- Rates are constant through time (no shifts or mass extinctions).
- The result depends on the quality of the branch lengths and their standard errors.

## Citation and license

Not yet published. No license has been chosen yet, so please contact the author before reusing the code.
