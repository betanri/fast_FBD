# fast_FBD

A fast Bayesian program for dating phylogenies with fossils, using the **fossilized birth-death (FBD)** process.

> **Status: research software.** The math and sampling have passed strict tests, and Fast FBD (given a phylogram) matches full-DNA FBD (given the DNA) on two test clades. See [Validation](#validation) before using results.

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

Status as of 12 November 2026. Details, numbers and result files: [VALIDATION.md](VALIDATION.md) and `validation/`.

### In plain language

`ffbd` puts ages on the nodes of a phylogeny using fossils. It uses the fossilized birth-death (FBD) model, which treats the origin of species, their extinction and the preservation of fossils as one process. Each fossil has an age range, and the program works out where on the tree each fossil fits. The tree of living species stays fixed. The program never reads the DNA itself. Instead it reads a phylogram: the branch lengths of that fixed tree, estimated once from the DNA, each with a standard error (SE). This shortcut is what lets it date trees of thousands of species in hours rather than months. Because the method is new, it has to earn trust. The validation follows the standard recipe for Bayesian software (Cook et al. 2006; Talts et al. 2018; Mendes et al. 2025) and for FBD simulation studies (Heath et al. 2014; Barido-Sottani et al. 2019). That recipe asks four questions. Each test is labeled by the data it uses: no data (fossils and priors only), a phylogram, or the full DNA.

**Does it compute the right probabilities?** For fixed trees and settings, every term was compared with RevBayes, the established FBD program. The totals matched to within 4e-9, which is rounding error (test 1).

**Does it explore those probabilities correctly?** Hundreds of data sets were simulated with known true answers, to check that the truth lands where it should among the program's estimates. This check is simulation-based calibration (SBC). It passed on small trees, on realistic trees of up to 8,056 living species, and on a 3-species case with an exact answer (tests 2a to 2c).

**Does it agree with other programs given the same model?** With no data and the same fixed rates, `ffbd` and MrBayes gave the same node ages (within about 0.5%) and the same number of sampled ancestors (test 3a). The two programs were written completely independently. Sampled ancestors are fossils that sit directly on a branch rather than on a side branch. RevBayes disagreed with both: it found about one extra sampled ancestor and made young nodes up to a third older. The likely cause is one RevBayes move, `mvCollapseExpandFossilBranch`, which is not reversible when fossils have age ranges. With that move switched off, RevBayes agrees with `ffbd` (test 3b). The RevBayes developers have not yet confirmed the problem. BEAST2 was not used, because it cannot hold the living-species tree fixed while fossils move, so it would be testing a different model.

**Does the phylogram shortcut give the same answer as analyzing the full DNA?** This is the question that matters most in practice.
- **Same DNA (test 4a).** On two fish clades, RevBayes analyzed the full 23-gene DNA in 6 long chains each, while `ffbd` got the phylogram from that same DNA. The ages agreed to 1.3 to 1.8% at the median. That is no larger than the differences between the DNA chains themselves.
- **MrBayes (test 4b).** A second, independent DNA check in MrBayes is running. It uses a slightly different clock, because MrBayes does not offer the one `ffbd` uses.
- **Simulated DNA (test 4c).** On 100 simulated trees, the true ages fell inside the 95% interval 91% of the time with the shortcut and 93% without it.
- **Catfishes (test 4d).** `ffbd` came close to an earlier full-DNA analysis, but that comparison is not like for like.
- **Controlled test (test 4e).** Everything stayed identical except the phylogram. A phylogram from the same 23 genes reproduced the 23-gene DNA result. Phylograms built from about 1,100 genes, or taken from a larger genome-scale tree, gave ages 2.5 to 3.8 million years younger on average, still inside the DNA uncertainty. How the SEs were computed made no difference.

Two further checks passed. Very large runs converged, with two independent chains agreeing (test 5). The optional rules that force a fossil inside or outside a crown group behave exactly as intended (test 6).

**What it all means.** The program does the math correctly, explores it correctly and agrees with an independent program. The phylogram shortcut costs essentially nothing compared with analyzing the same DNA. Which genes the branch lengths come from matters more than the shortcut. The known weak spot is shared with full-DNA FBD: very young nodes, under about 5 million years, tend to come out roughly 1 million years too old. This comes from the model and its priors, not from `ffbd` or its shortcut.

### What kind of data each test uses

`ffbd` never reads DNA. Its data are a **phylogram**: the branch lengths of the fixed tree, in substitutions per site, each with a standard error, estimated once from the DNA. Every test below says which of three kinds of data it uses:

- **No data**: fossils and priors only. These tests check the tree and fossil part of the model, which is where programs differ most.
- **Phylogram**: branch lengths and standard errors, either estimated from DNA or simulated.
- **DNA**: the full sequence alignment, used by RevBayes or MrBayes as the reference. `ffbd` is then given a phylogram estimated from that same DNA.

### How a Bayesian dating program is validated

The standard checks for a new Bayesian method (Cook et al. 2006; Talts et al. 2018; Mendes et al. 2025) ask four questions. Simulation studies of FBD dating follow the same logic (Heath et al. 2014; Barido-Sottani et al. 2019).

1. **Does it calculate the right probabilities?** For fixed trees and parameter values, compare every term with an established program.
2. **Does it sample those probabilities correctly?** Simulate many data sets from known true values and check that the truth falls inside the 95% interval about 95% of the time. In simulation-based calibration, the rank of the true value among the posterior samples must be uniform.
3. **Does it agree with other programs on the same model?**
4. **Does the shortcut work?** Fast FBD replaces the DNA with a phylogram. On the same data, it should give the same ages as an analysis of the full DNA.

### Results

| # | Question | Test | Data | Result |
|---|---|---|---|---|
| 1 | Right probabilities | Every probability term against RevBayes 1.3.2, on fixed trees (4 data sets, 89 states) | No data and phylogram | **Passed**: difference at most 4e-9 |
| 2a | Correct sampling | Simulation-based calibration, small trees (500 replicates in each of 3 settings) | No data, and simulated phylogram | **Passed** |
| 2b | Correct sampling | Simulation-based calibration, realistic trees (500 replicates, 2 to 8,056 living species) | Simulated phylogram | **Passed**: all 9 quantities uniform |
| 2c | Correct sampling | A 3-species prior with an exact answer | No data | **Passed** |
| 3a | Same answer as other programs | ffbd, MrBayes 3.2.7 and RevBayes 1.3.2 on the same fixed-rate prior (two fish clades) | No data | **ffbd = MrBayes** (node ages within 0.4 to 0.6% at the median, same number of sampled ancestors). **RevBayes differs** (3 to 7% at the median, up to 34%), because of one RevBayes move (below) |
| 3b | Same answer as other programs | ffbd and RevBayes with that RevBayes move switched off | No data | **Agree**: 0.3 to 0.4% |
| 4a | Shortcut works | ffbd against RevBayes full-DNA FBD, two fish clades (23 genes, 16,107 sites; 6 DNA chains per clade) | **DNA** (RevBayes) vs phylogram from the same DNA (ffbd) | **Agree**: node ages r = 0.9995 and 0.9997, median difference 1.3% and 1.8%, mean +0.5 and +0.7 Myr. That is within the noise between the DNA chains themselves (1.5% and 3.2%) |
| 4b | Shortcut works | ffbd against MrBayes full-DNA FBD, the same two clades and genes | **DNA** (MrBayes) vs phylogram (ffbd) | **Running**: Protacanthopterygii chains finished, Gobiaria chains about half done |
| 4c | Shortcut works | 100 simulated trees (30 to 120 tips); ffbd given phylograms estimated from simulated DNA (2,000 sites) | Phylogram from **simulated DNA** | **Close to nominal**: true node ages inside the 95% interval 91% of the time (93% when the exact phylogram is used); nodes under 5 Ma about 1 Myr too old |
| 4d | Shortcut works | Catfishes (713 living species, 94 fossils) against an earlier full-DNA RevBayes analysis | **DNA** (RevBayes, 21 genes) vs phylogram (ffbd) | **Close**: r = 0.994, ffbd 1.8 Myr younger on average. Not like for like: the earlier analysis used another clock and the RevBayes move below |
| 4e | Which phylogram to use | Controlled test on the two clades of 4a: same tree, fossils and priors, only the phylogram changes. Reference: the 23-gene full DNA of 4a | **DNA** (RevBayes, 23 genes) vs phylograms from 23 genes, from all genes (about 580,000 sites), and from a genome-scale phylogram of a larger tree pruned to the clade | **23-gene phylogram = 23-gene DNA** (median 1.0% and 1.5%, within DNA chain noise). **More genes give younger ages**: all-gene phylogram 5 to 6% at the median, 2.5 to 2.7 Myr younger on average; pruned larger-tree phylogram 8 to 9%, 3.7 to 3.8 Myr younger. Every node stays inside the 23-gene DNA 95% interval. The standard-error rule makes no difference |
| 5 | Large trees | Two chains each at 1,185 and 2,424 tips, 20 billion proposals | Phylogram | **Converged**: effective sample sizes 500 to 800 for the main parameters; chains agree on node ages to 0.9% (median) |
| 6 | Crown and stem rules | Constrained runs against unconstrained runs filtered to the same rules | Phylogram (small example) | **Passed**: all quantities within 1.4 standard errors |

### Figures

![Node ages, full-DNA FBD against Fast FBD](validation/node_ages_full_vs_fast.png)

*Tests 4a and 4d: DNA against phylogram. Each point is one node of the living-species time tree (median ages), matched by clade. Gobiaria and Protacanthopterygii: RevBayes 1.3.2 full-DNA FBD (23 genes, 16,107 sites; 6 chains; collapse-expand move removed) against `ffbd` given the phylogram from the same alignment, with the same fossils and priors. Catfishes: the earlier full-DNA RevBayes analysis (21 genes, lognormal clock, original moves) against `ffbd`. Dashed line: equal ages.*

![Phylograms against 23-gene DNA](validation/phylogram_vs_dna_controlled.png)

*Test 4e: DNA against several phylograms. Difference of each node's median age from the 23-gene full-DNA FBD. Grey: two halves of the DNA chains against each other (the noise floor). Blue: the phylogram from the same 23 genes. Orange: the phylogram from all genes. Green: a genome-scale phylogram of a larger tree, pruned to the clade.*

![ffbd, MrBayes and RevBayes on the same prior](validation/three_programs_no_data.png)

*Test 3a: no data. Top: difference of each node's median age from MrBayes. Bottom: how many fossils are sampled ancestors in each posterior tree. ffbd and MrBayes overlap; RevBayes finds more sampled ancestors and older young nodes.*

![Simulation-based calibration](validation/sbc_rank_histograms.png)

*Test 2b: simulated phylograms. Rank of the true value among 99 posterior draws, 500 replicates. A correct sampler gives flat histograms; dotted lines mark the range expected by chance (99%).*

![Simulation study, coverage by node age](validation/simulation_coverage_by_age.png)

*Test 4c: simulated phylograms, exact and estimated from simulated DNA. Left: share of true node ages inside the 95% interval, by true age. Right: median error of the node-age estimate. Young nodes are too old in both versions, so the cause is the model and prior, not the phylogram shortcut.*

### Why RevBayes differs, and what MrBayes and BEAST2 show

With no data, ffbd and RevBayes 1.3.2 disagree whenever fossils have age ranges and can become sampled ancestors, meaning fossils that lie directly on a branch rather than on a side branch. RevBayes finds more sampled ancestors and makes young nodes older, by up to a third on some two-species nodes. The probabilities are identical (test 1), and the programs agree once one RevBayes move, `mvCollapseExpandFossilBranch`, is switched off (test 3b). In the RevBayes 1.3.2 source, that move draws the new attachment age above the maximum of the fossil's age window. The reverse move can start from attachments inside the window, so the move is not reversible.

MrBayes 3.2.7 (Zhang et al. 2016) is independent code. Given the same fixed tree, fossils and fixed rates, with no data, it agrees with ffbd and not with RevBayes (test 3a). Test 3a uses no data because the disagreement comes from how trees and fossils are proposed, not from the DNA. A full-DNA comparison with MrBayes is test 4b.

BEAST2's sampled-ancestor package (Gavryushkina et al. 2014) was not used. Its clade constraints count fossils as clade members, so it cannot hold the living-species tree fixed while fossils move within their clades.

The RevBayes developers have not yet confirmed the move problem. `--rb-sa-move` imitates that proposal, for testing only.

### What this means

- **ffbd computes and samples the FBD model correctly** (tests 1, 2): its probabilities match RevBayes to 4e-9, and simulations with known answers pass on trees of up to 8,056 living species.
- **With no data, an independent program, MrBayes, gives the same answer as ffbd.** RevBayes 1.3.2 with its default moves does not.
- **With DNA, the phylogram shortcut gives the same ages as full-DNA FBD** on two clades (test 4a): the differences are as small as those between DNA chains. The MrBayes DNA check (test 4b) is running. On simulated DNA, interval coverage is close to the target 95%.
- **The genes behind the phylogram matter more than the shortcut** (test 4e). A phylogram from the same 23 genes reproduces the 23-gene DNA analysis. Phylograms from many more genes give ages a few million years younger, inside the 23-gene intervals.
- **The weak spot is very young nodes**: in simulations, nodes under 5 Ma come out about 1 Myr too old, with or without the shortcut.

### References

- Barido-Sottani J, Aguirre-Fernández G, Hopkins MJ, Stadler T, Warnock R (2019). Ignoring stratigraphic age uncertainty leads to erroneous estimates of species divergence times under the fossilized birth-death process. *Proceedings of the Royal Society B* 286: 20190685.
- Cook SR, Gelman A, Rubin DB (2006). Validation of software for Bayesian models using posterior quantiles. *Journal of Computational and Graphical Statistics* 15: 675-692.
- Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014). Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. *PLoS Computational Biology* 10: e1003919.
- Heath TA, Huelsenbeck JP, Stadler T (2014). The fossilized birth-death process for coherent calibration of divergence-time estimates. *PNAS* 111: E2957-E2966.
- Höhna S, Landis MJ, Heath TA, et al. (2016). RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language. *Systematic Biology* 65: 726-736.
- Mendes FK, Bouckaert R, Carvalho LM, Drummond AJ (2025). How to validate a Bayesian evolutionary model. *Systematic Biology*, doi:10.1093/sysbio/syae064.
- Ronquist F, Teslenko M, van der Mark P, et al. (2012). MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. *Systematic Biology* 61: 539-542.
- Talts S, Betancourt M, Simpson D, Vehtari A, Gelman A (2018). Validating Bayesian inference algorithms with simulation-based calibration. arXiv:1804.06788.
- Zhang C, Stadler T, Klopfstein S, Heath TA, Ronquist F (2016). Total-evidence dating under the fossilized birth-death process. *Systematic Biology* 65: 228-249.

## Known limits

- The topology of living species cannot change.
- Fossils must fall inside a clade constraint below the root in root-conditioned runs.
- Rates are constant through time (no shifts or mass extinctions).
- The result depends on the quality of the branch lengths and their standard errors.

## Citation and license

Not yet published. No license has been chosen yet, so please contact the author before reusing the code.
