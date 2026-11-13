# Validation of ffbd

Status as of 12 November 2026. "RevBayes" means RevBayes 1.3.2 and "MrBayes" means MrBayes 3.2.7. Result files and figures are in `validation/`. The numbering follows the Validation section of the [README](README.md), which also explains the approach and gives the references.

**Data.** `ffbd` never reads DNA. Its data are a phylogram: branch lengths of the fixed living-species tree with standard errors. Each test says which data it uses:

- **No data**: fossils and priors only.
- **Phylogram**: branch lengths and standard errors, estimated from DNA or simulated.
- **DNA**: the full alignment, analyzed by RevBayes or MrBayes as the reference, while `ffbd` gets a phylogram estimated from the same DNA.

**Data sets.** Two small fish clades: Gobiaria (A06: 24 living species, 10 fossils) and Protacanthopterygii (A12: 22 living species, 14 fossils), each with a 23-gene alignment of 16,107 sites; the RevBayes bear tutorial data; and catfishes (713 living species, 94 fossils).

## Summary

| # | Test | Data | Result |
|---|---|---|---|
| 1 | Probabilities equal RevBayes on fixed trees | No data and phylogram | Passed |
| 2a | Simulation-based calibration, small trees | No data and simulated phylogram | Passed |
| 2b | Simulation-based calibration, realistic trees (up to 8,056 tips) | Simulated phylogram | Passed |
| 2c | Exact 3-species prior | No data | Passed |
| 3a | ffbd, MrBayes and RevBayes on the same prior | No data | ffbd and MrBayes agree; RevBayes differs |
| 3b | ffbd and RevBayes with one RevBayes move off | No data | Agree |
| 4a | Fast FBD against RevBayes full-DNA FBD | DNA vs phylogram | Agree, within DNA chain noise |
| 4b | Fast FBD against MrBayes full-DNA FBD | DNA vs phylogram | Gobiaria agrees; on Protacanthopterygii MrBayes's different clock gives younger ages than both RevBayes DNA and ffbd |
| 4c | Simulation study, 100 trees of 30 to 120 tips | Phylogram from simulated DNA | Node and fossil ages near 95% coverage |
| 4d | Catfishes against an earlier full-DNA analysis | DNA vs phylogram | Close, not identical |
| 4e | Controlled test: which phylogram | DNA vs several phylograms | 23-gene phylogram = 23-gene DNA; more genes give younger ages |
| 5 | Convergence on large trees | Phylogram | Passed |
| 6 | Crown and stem rules | Phylogram | Passed |

## 1. Probabilities on fixed trees (passed)

*Data: no data and phylogram.* For fixed trees and parameter values, every term was compared with RevBayes: FBD density, root prior, rate priors, data terms and fossil windows.

| Data set | States compared | Largest difference in the total |
|---|---|---|
| A06 | 24 | 3.7e-10 |
| A12 | 24 | 2.8e-9 |
| Bears | 17 of 18 | 2.0e-9 |
| Catfish | 24 | 3.7e-9 |

RevBayes refused to build one bear state, because its origin was below a fossil's window maximum. That is a check at model construction, not a difference in density.

## 2. Correct sampling

### 2a. Simulation-based calibration, small trees (passed)

*Data: no data, and simulated phylograms.* Data sets were simulated exactly from the model, with known true values, using `ffbd_sim` (separate code from the sampler). If the sampler is correct, the rank of each true value among the posterior samples is uniform. There were 500 replicates per setting, and 13 parameters were tested by chi-square and Kolmogorov-Smirnov tests with a Bonferroni threshold of 0.0038.

| Setting | Smallest p-value | Result |
|---|---|---|
| Root-conditioned, with data | 0.014 | Passed |
| Root-conditioned, no data | 0.079 | Passed |
| Origin-conditioned, with data | 0.014 | Passed |

Testing found and fixed two bugs: a cached sum of rates that could become corrupted, and a checkpoint restart that was not exactly reproducible. All tests were rerun after the fixes. A restarted run now gives byte-identical output.

### 2b. Simulation-based calibration at realistic size (passed)

*Data: simulated phylograms.* `scripts/sim_study.py sbc`: each of 500 replicates draws speciation, extinction, fossil sampling and the mean clock rate from the `ffbd` priors (exponential, rate 10) and the root from the analysis prior Normal(30, 2), then simulates an FBD tree with no size filter. Only information that is exact for this test is used: fossil windows on a fixed 4-Myr grid, one constraint per side of the root, and a constant observation standard deviation (0.01). Trees had 2 to 8,056 living tips (median 20) and 0 to 15,838 fossils (median 25). 100 million proposals per replicate, 99 thinned draws, rank of the true value among them (`scripts/sbc_ranks.py`; `validation/sbc_realistic_500reps.txt`; figure `validation/sbc_rank_histograms.png`). With a correct sampler the ranks are uniform. Bonferroni threshold for 9 tests: p < 0.0056.

| Quantity | n | Uniformity p (20 bins) | Mean rank (0.5 expected) | Median ESS |
|---|---|---|---|---|
| Speciation | 500 | 0.66 | 0.490 | 361 |
| Extinction | 500 | 0.50 | 0.503 | 361 |
| Fossil sampling | 500 | 0.65 | 0.488 | 361 |
| Mean clock rate | 500 | 0.087 | 0.485 | 320 |
| Root age | 500 | 0.79 | 0.493 | 361 |
| Number of sampled ancestors | 500 | 0.91 | 0.504 | 361 |
| Crown age, root side 1 | 407 | 0.24 | 0.522 | 172 |
| Crown age, root side 2 | 407 | 0.17 | 0.513 | 180 |
| One random fossil's age | 482 | 0.87 | 0.494 | 361 |

All ranks are uniform: `ffbd` samples the FBD posterior correctly on trees of up to thousands of tips.

### 2c. Exact 3-species prior (passed)

*Data: no data.* An exact analytic check of a 3-species, 1-fossil prior agreed: sampled-ancestor probability 0.6871 against 0.6889, mean fossil age 4.913 against 4.914.

## 3. Agreement with other programs

### 3a. ffbd, MrBayes and RevBayes on the same prior (ffbd agrees with MrBayes)

*Data: no data.* `scripts/make_fixed3.py` writes the same model for all three programs: the A06 or A12 living-species tree fixed (in MrBayes, one partial constraint per living-species node, fossils free), the same clade constraints, uniform fossil age windows (MrBayes `calibrate`), speciation 0.10, extinction 0.04, fossil sampling 0.03 and the root age all **fixed**. With every rate fixed, the three programs must give the same node ages and the same number of sampled ancestors. Two chains per program and data set; comparison `scripts/compare_fixed3.py` (`validation/three_program_prior_A06_A12.txt`; figure `validation/three_programs_no_data.png`).

| Data set | Mean sampled ancestors, ffbd / MrBayes / RevBayes | Node ages, ffbd vs MrBayes (median, largest) | RevBayes vs MrBayes (median, largest) | MrBayes chain 1 vs chain 2 (median, largest) |
|---|---|---|---|---|
| A06 | 7.01 / 7.02 / 7.95 | 0.36%, 2.7% | 3.2%, 34% | 0.89%, 7.4% |
| A12 | 9.55 / 9.58 / 10.13 | 0.59%, 4.2% | 6.5%, 28.5% | 0.86%, 6.4% |

`ffbd` and MrBayes agree within chain noise. RevBayes gives about one extra sampled ancestor and older young nodes (largest on two-tip nodes, for example 10.6 Ma in MrBayes and 14.2 Ma in RevBayes).

The same pattern appears when RevBayes and `ffbd` are compared directly, with estimated rates:

| Comparison | Data | Median difference in node ages | Largest difference |
|---|---|---|---|
| A06 | No data | 1.64% | 12.8% |
| A12 | No data | 4.97% | 8.1% |
| A06 | Phylogram | 1.13% | 21.5% |
| A12 | Phylogram | 2.59% | 13.0% |
| Bears | Phylogram | 8.96% | 12.6% |

RevBayes node-age ESS was low on these runs (32 to 35 per chain on the phylogram runs), so part of each difference is RevBayes sampling noise.

**Explanation.** The RevBayes move `mvCollapseExpandFossilBranch` (`CollapseExpandFossilBranchProposal.cpp`, v1.3.2). When it adds a sampled ancestor's branch, it draws the new age from above the maximum of the fossil's age window. The reverse move can start from attachments inside the window, and those cannot be reached in the forward direction, so the move is not reversible for fossils with age ranges. Four results support this: the probabilities are identical (test 1), the agreement is good once the move is removed (test 3b), `ffbd` recovers truth in simulation (test 2), and MrBayes, written independently, agrees with `ffbd` (this test). This has not yet been confirmed with the RevBayes developers.

**BEAST2.** The sampled-ancestor package of BEAST2 was not used. Its clade constraints count fossils as clade members, so it cannot hold the living-species tree fixed while fossils move within their clades. A BEAST2 run would be a different model, so MrBayes was used as the independent program.

### 3b. ffbd and RevBayes with the collapse-expand move off (passed)

*Data: no data.*

| Data set | Node-age medians, median difference | Largest difference | Two RevBayes chains differ by (median) |
|---|---|---|---|
| A06 | 0.28% | 1.44% | 0.59% |
| A12 | 0.41% | 1.44% | 0.70% |

## 4. Does the phylogram shortcut work?

### 4a. Fast FBD against RevBayes full-DNA FBD (agree)

*Data: DNA (RevBayes) against a phylogram from the same DNA (`ffbd`).* RevBayes full-DNA FBD (GTR+G, exponential clock, fixed living-species tree, 23 genes, 16,107 sites) and `ffbd` on A06 and A12, with the same fossils, constraints, priors and uniform root prior. The RevBayes collapse-expand move is removed (test 3), so the only remaining difference is the branch-length shortcut. Each side was summarized as an extant time tree (median node ages), and nodes were matched by clade (`scripts/branch_scatter.py --ages-only`; `validation/node_ages_full_vs_fast.tsv`; figure `validation/node_ages_full_vs_fast.png`).

| Data set | Nodes | r | Median difference | Mean difference (Fast minus full) | Root age, full vs Fast | DNA chains 1-3 vs 4-6 (median) |
|---|---|---|---|---|---|---|
| A06 | 23 | 0.9995 | 1.3% | +0.5 Myr | 75.2 vs 75.2 Ma | 1.5% |
| A12 | 21 | 0.9997 | 1.8% | +0.7 Myr | 139.0 vs 139.0 Ma | 3.2% |

The differences are no larger than those between the DNA chains themselves.

| Set | Side | Chains | Length per chain | Run time per chain |
|---|---|---|---|---|
| A06, A12 | Full DNA, RevBayes | 6 | 250,000 generations | ~22 to 25 h |
| A06, A12 | Phylogram, `ffbd` | 2 | 400 million proposals | ~2.5 min |

Full-DNA node-age ESS per chain is low (median 19 for A06 and 35 for A12 on the first chains), which is why six chains were pooled.

### 4b. Fast FBD against MrBayes full-DNA FBD (Gobiaria agrees; Protacanthopterygii shows clock sensitivity)

*Data: DNA (MrBayes) against a phylogram (`ffbd`), with RevBayes DNA (test 4a) alongside.* `scripts/make_mb_dna.py` sets up MrBayes 3.2.7 on the same 23-gene alignment (16,107 sites, GTR+G4), the fixed living-species tree (partial constraints), the same clade constraints, fossil windows, sampling fraction and root prior. Two settings differ because MrBayes requires it:

- **Clock.** MrBayes has no exponential-per-branch clock, so the independent gamma rates (IGR) clock is used.
- **FBD priors.** These are set on MrBayes's own parameters (net diversification, turnover, sampling proportion).

A no-data run with the same priors shows how much of each age comes from the DNA. There are two DNA chains and two no-data chains per clade, 50 million generations each, sampled every 5,000, with the first 10% dropped. MrBayes trees are in substitutions, so each node age is the node height divided by that sample's clock rate. Scoring: `scripts/compare_mbdna.py` (`validation/mrbayes_dna_A06_A12.txt`; per node and program: `validation/plotdata_mbdna_A06.tsv`, `plotdata_mbdna_A12.tsv`).

| Comparison (node-age medians) | A06 median % (largest), mean shift | A12 median % (largest), mean shift |
|---|---|---|
| Noise: MrBayes DNA chain 1 vs chain 2 | 0.5% (1.7%), -0.1 Myr | 4.1% (17%), +1.3 Myr |
| Noise: RevBayes DNA chains 1-3 vs 4-6 | 1.4% (10.6%), -0.9 Myr | 3.3% (13.7%), +1.4 Myr |
| MrBayes DNA vs RevBayes DNA | 3.1% (21%), -0.2 Myr | 20% (42%), -10.3 Myr |
| `ffbd` 23-gene phylogram vs MrBayes DNA | 2.8% (25%), +0.6 Myr; 23 of 23 inside the MrBayes 95% interval | 28% (77%), +10.9 Myr; 14 of 21 inside |
| `ffbd` 23-gene phylogram vs RevBayes DNA | 1.0% (4.6%), +0.5 Myr | 1.5% (3.4%), +0.7 Myr |
| MrBayes no data vs MrBayes DNA | 8% | 54% |
| Root age (Ma): RevBayes DNA / MrBayes DNA / `ffbd` | 75.2 / 74.0 / 75.2 | 138.9 / 138.8 / 139.0 |

**Gobiaria (A06).** All three agree, within about twice the chain noise.

**Protacanthopterygii (A12).** The two DNA programs differ beyond their chain noise.

- MrBayes puts the inner nodes 5 to 24 Myr younger. For example, one 11-tip node is 91 Ma in RevBayes and 68 Ma in MrBayes (95% interval 53 to 86).
- MrBayes intervals are about twice as wide.
- The root, held by its prior, agrees.
- `ffbd` stays within 1.5% of RevBayes, whose clock and priors it shares.
- The DNA moves the MrBayes ages far from its no-data run (54% at the median), so the difference comes from how the DNA is read under each clock, not from the priors alone.

**Reading.** Test 4b does not single out the shortcut: given the same model, the shortcut reproduces full-DNA FBD (tests 4a and 4e). It shows that the clock model can change node ages more than the shortcut does.

### 4c. Simulation study at realistic size

*Data: simulated phylograms, exact and estimated from simulated DNA.* After Heath et al. (2014). `scripts/sim_study.py simulate` grows FBD trees forward from a known root (100 Ma; speciation 0.10, extinction 0.06, fossil sampling 0.01, extant sampling 0.5; trees kept at 30 to 120 living tips with at least 5 fossils), gives each fossil an age window of 2 to 8 Myr containing its true age and a clade constraint that is true on the true tree, and draws branch rates from an exponential with mean 0.004. Two versions of every tree:

- **exact**: branch lengths drawn from the `ffbd` likelihood itself (tests the sampler).
- **dna**: 2,000 sites simulated with AliSim (GTR+G4); ML branch lengths and 100 bootstrap standard errors on the true topology, exactly as for the real clades (`scripts/fish_fast.py`; tests the shortcut).

Root prior soft-bounded on 90 to 110 Ma; 200 million proposals per run; 100 trees. Scoring: `scripts/sim_study_summary.py` (`validation/sim_study_100reps.txt`; figure `validation/simulation_coverage_by_age.png`). Coverage of the 95% interval:

| Version | Node ages (6,911) | Fossil ages (3,631) | Root | Speciation | Extinction | Fossil sampling | Sampled ancestors | Mean clock rate |
|---|---|---|---|---|---|---|---|---|
| exact | 0.931 | 0.956 | 0.99 | 0.58 | 0.60 | 0.89 | 0.83 | 0.87 |
| dna | 0.914 | 0.954 | 0.99 | 0.56 | 0.54 | 0.92 | 0.82 | 0.38 |

- Node and fossil ages are close to 95% coverage, and the DNA shortcut costs only 1 to 2 points.
- Nodes younger than 5 Ma come out too old in both versions (median error +56% in exact and +59% in dna, about 1.2 Myr). The error falls with age: +10 to 12% at 15 to 30 Ma, +1 to 2% above 60 Ma. Because both versions show it, this comes from the model and prior, not from the shortcut.
- Speciation and extinction are under-covered, and the dna version overestimates the mean clock rate. Test 2b shows the sampler itself is calibrated, so this comes from the design of this study: one fixed set of true values, trees kept only at 30 to 120 tips, and exponential priors that pull the rates down.

### 4d. Catfishes (close, not identical)

*Data: DNA (an earlier RevBayes analysis, 21 genes) against a phylogram (`ffbd`).* 553 nodes matched by clade, r = 0.994, `ffbd` 1.8 Myr younger on average (median difference 10.5%, largest on young nodes). The earlier analysis used an uncorrelated lognormal clock in RevBayes with the collapse-expand move, so it is not a like-for-like comparison.

### 4e. Controlled test: which phylogram (23-gene phylogram = 23-gene DNA)

*Data: DNA (RevBayes, 23 genes; the reference of test 4a) against several phylograms given to `ffbd`.* For A06 and A12, every run has the same taxa, fixed living-species tree, fossils, constraints, priors and root prior; only the phylogram changes. Phylograms (IQ-TREE 2.3.6, GTR+G4, branch lengths on the fixed tree):

- **23 genes**: the 16,107 sites of test 4a, 500 bootstrap standard errors.
- **All genes**: every gene of the genome-scale data set (1,132 to 1,133 genes, about 580,000 sites), 500 bootstrap standard errors.
- **All genes, SE rule**: the same lengths with the standard errors used for large trees, sqrt(b / 20,000) with a floor of 1 / 20,000.
- **Larger-tree phylogram, pruned**: branch lengths of a genome-scale phylogram estimated on a much larger tree, pruned to the clade, with the SE rule. Two such phylograms were tested and have identical branch lengths within these clades, so they are shown once.

Each input was checked: the converter rebuilds the original 4a input exactly, and each phylogram input differs from it only in its branch-length lines. Two `ffbd` chains per phylogram, 400 million proposals. Comparison `scripts/ctrl_compare.py` (`validation/phylogram_vs_dna_controlled.txt`; figure `validation/phylogram_vs_dna_controlled.png`).

| Data set | Phylogram | Median difference from DNA | Largest | Mean shift | Nodes inside the DNA 95% interval | Root (Ma) |
|---|---|---|---|---|---|---|
| A06 | DNA chains 1-3 vs 4-6 (noise) | 1.5% | 11.9% | | 23 of 23 | 75.4 |
| A06 | 23 genes | 1.0% | 4.6% | +0.5 Myr | 23 of 23 | 75.2 |
| A06 | All genes | 6.3% | 26.8% | -2.7 Myr | 23 of 23 | 75.4 |
| A06 | All genes, SE rule | 6.0% | 26.8% | -2.7 Myr | 23 of 23 | 75.4 |
| A06 | Larger tree, pruned | 8.1% | 38.7% | -3.8 Myr | 23 of 23 | 75.4 |
| A12 | DNA chains 1-3 vs 4-6 (noise) | 3.2% | 12.0% | | 21 of 21 | 138.8 |
| A12 | 23 genes | 1.5% | 3.4% | +0.7 Myr | 21 of 21 | 139.0 |
| A12 | All genes | 5.0% | 46.2% | -2.5 Myr | 21 of 21 | 139.4 |
| A12 | All genes, SE rule | 4.7% | 45.6% | -2.5 Myr | 21 of 21 | 139.5 |
| A12 | Larger tree, pruned | 8.8% | 48.6% | -3.7 Myr | 21 of 21 | 139.4 |

- **Same genes, DNA or phylogram: no difference.** With the 23-gene phylogram, ages match the 23-gene DNA analysis as closely as the DNA chains match each other.
- **More genes give younger ages.** The all-gene and larger-tree phylograms make nodes 2.5 to 3.8 Myr younger on average, most on young nodes, but every node stays inside the 23-gene DNA 95% interval. This is a difference in the data, not in the method: the all-gene trees are 17 to 24% shorter than the 23-gene trees and have different relative branch lengths (log-log r 0.81 to 0.91 against the 23-gene lengths).
- **The standard errors hardly matter.** Bootstrap standard errors and the large-tree rule give the same ages (within 0.3% at the median).
- **Precision of the inputs.** Median SE/b is 0.116 (A06) and 0.121 (A12) for the 23-gene phylograms against 0.023 and 0.019 for the all-gene phylograms. Mitochondrial genes (COI, CYTB in the 23-gene set; 10 mitochondrial genes in the all-gene set) make up less than 2% of the all-gene sites. The true ages of these clades are unknown, so this test shows that the input changes the ages, not which input is more accurate (see the README, "Why ffbd may do better than full-DNA FBD on real data").

## 5. Convergence on large trees (passed)

*Data: phylogram.* Two chains each at 1,185 and 2,424 tips, 20 billion proposals each (about 2 hours on one core): effective sample sizes 500 to 800 for the main parameters, and at least 139 for node ages. Node-age medians agree between chains to 0.9% (median), and 95% of nodes agree within 5%.

## 6. Crown and stem rules (passed)

*Data: phylogram (small example).* The rules restrict where a fossil may attach: only below the crown node, or never below it. They add no density term. Proposals that would break a rule are rejected, which leaves the target distribution as the unconstrained one restricted to the allowed states. Two checks:

1. **Exact conditioning test.** On the small example (`examples/example_crown_stem.txt`), the run with a crown rule (FOSSIL_1) and a stem rule (FOSSIL_3) was compared with unconstrained runs keeping only samples where both rules held. With 11 seeds per side and 100,001 samples each, the rules held in every sample of the constrained runs. All quantities agreed within 1.4 standard errors (`validation/crown_stem_test_11seeds.txt`). A crown-only test with 2 seeds also agreed. `--check` runs, which validate the full state after every proposal, raised no error.
2. **Regression.** With no crown or stem rules, the new version gives byte-identical log and tree files to the previous release on A06, A12, bears and catfish, with and without DNA (`validation/regression_previous_release.txt`).

## Speed

0.39 seconds per million proposals on one core of a compute cluster, about 2.6 million proposals per second, and 0.13 to 0.16 seconds on a laptop. Speed was the same at 1,185 and 2,424 tips.
