# Validation of ffbd

Status as of 11 November 2026. "RevBayes" means RevBayes 1.3.2. Result files are in `validation/`.

Data sets: two small fish clades (A06: 24 living species and 10 fossils; A12: 22 and 14), the RevBayes bear tutorial data, and a catfish data set (713 living tips, 94 fossils).

## Summary

| # | Test | Result |
|---|---|---|
| 1 | Probabilities equal RevBayes on fixed trees | Passed |
| 2 | Simulation-based calibration and analytic prior | Passed |
| 3 | Agreement with RevBayes, sampled-ancestor changes off | Passed |
| 4 | Agreement with RevBayes, sampled-ancestor changes on | Does not match (section 4) |
| 5 | Fast FBD against full-DNA FBD | Preliminary: agrees, full-DNA chains not yet well converged |
| 6 | Catfish comparison | Close, not identical (section 6) |
| 7 | Convergence on large trees | Passed |
| 8 | Crown and stem rules | Passed |
| 9 | Three programs on the same prior (ffbd, MrBayes 3.2.7, RevBayes) | ffbd and MrBayes agree; RevBayes differs |
| 10 | Simulation study, 100 trees of 30 to 120 tips | Node and fossil ages near 95% coverage; speciation and extinction under-covered (explained by test 11) |
| 11 | Simulation-based calibration at realistic size (500 replicates, up to 8,056 tips) | Passed |

## Bottom line: ffbd, MrBayes, RevBayes and BEAST2

- **`ffbd` and MrBayes agree.** Two independently written programs, given the same tree, fossils and fixed rates, give the same node ages (within about 0.5%) and the same number of sampled ancestors (7.01 vs 7.02; 9.55 vs 9.58).
- **RevBayes 1.3.2 disagrees with both.** It finds about one extra sampled ancestor and makes young nodes older, by up to a third on some two-species nodes. With one RevBayes move switched off (`mvCollapseExpandFossilBranch`), RevBayes agrees with `ffbd` (0.3 to 0.4%). Our reading is that this move is not reversible when fossils have age ranges. The RevBayes developers have not yet confirmed this.
- **BEAST2 was not run.** `ffbd` keeps the tree of living species fixed while each fossil moves freely within its clade. BEAST2's sampled-ancestor package cannot express that, because its clade constraints count fossils as clade members. A BEAST2 run would therefore be a different model, not a check of this one, so MrBayes was used as the independent program.
- **`ffbd` recovers known truth.** On 500 simulated data sets with trees of up to 8,056 living species, the true values sit where they should in `ffbd`'s estimates. On 100 simulated DNA data sets, node and fossil ages have close to 95% coverage.
- **Practical meaning.** For FBD analyses with fossil age ranges and sampled ancestors, `ffbd` and MrBayes give the same answer. RevBayes 1.3.2 with its default moves can give more sampled ancestors and older young nodes.

## 1. Probabilities on fixed trees (passed)

For fixed trees and parameter values, every term was compared with RevBayes: FBD density, root prior, rate priors, data terms and fossil windows.

| Data set | States compared | Largest difference in the total |
|---|---|---|
| A06 | 24 | 3.7e-10 |
| A12 | 24 | 2.8e-9 |
| Bears | 17 of 18 | 2.0e-9 |
| Catfish | 24 | 3.7e-9 |

RevBayes refused to build one bear state, because its origin was below a fossil's window maximum. That is a check at model construction, not a difference in density.

## 2. Simulation-based calibration (passed)

Data sets were simulated exactly from the model, with known true values, using `ffbd_sim` (separate code from the sampler). If the sampler is correct, the rank of each true value among the posterior samples is uniform. There were 500 replicates per setting, and 13 parameters were tested by chi-square and Kolmogorov-Smirnov tests with a Bonferroni threshold of 0.0038.

| Setting | Smallest p-value | Result |
|---|---|---|
| Root-conditioned, with data | 0.014 | Passed |
| Root-conditioned, no data | 0.079 | Passed |
| Origin-conditioned, with data | 0.014 | Passed |

An exact analytic check of a 3-species, 1-fossil prior also agreed: sampled-ancestor probability 0.6871 against 0.6889, mean fossil age 4.913 against 4.914.

Testing found and fixed two bugs: a cached sum of rates that could become corrupted, and a checkpoint restart that was not exactly reproducible. All tests above were rerun after the fixes. A restarted run now gives byte-identical output.

## 3. Agreement with RevBayes, sampled-ancestor changes off (passed)

Both programs were run on the prior with the collapse-expand move switched off.

| Data set | Node-age medians, median difference | Largest difference | Two RevBayes chains differ by (median) |
|---|---|---|---|
| A06 | 0.28% | 1.44% | 0.59% |
| A12 | 0.41% | 1.44% | 0.70% |

## 4. Agreement with RevBayes, sampled-ancestor changes on (does not match)

| Comparison | Median difference in node ages | Largest difference |
|---|---|---|
| A06, no data | 1.64% | 12.8% |
| A12, no data | 4.97% | 8.1% |
| A06, Fast FBD | 1.13% | 21.5% |
| A12, Fast FBD | 2.59% | 13.0% |
| Bears, Fast FBD | 8.96% | 12.6% |

Compared with `ffbd`, RevBayes gives lower extinction rates, more sampled ancestors and older node ages.

Our explanation is the RevBayes move `mvCollapseExpandFossilBranch` (`CollapseExpandFossilBranchProposal.cpp`, v1.3.2). When it adds a sampled ancestor's branch, it draws the new age from above the maximum of the fossil's age window. The reverse move can start from attachments inside the window, and those cannot be reached in the forward direction, so the move is not reversible for fossils with age ranges.

Three results support this: the probabilities are identical (section 1), the agreement is good once the move is removed (section 3), and `ffbd` recovers truth in simulation (section 2). This has not yet been confirmed with the RevBayes developers.

RevBayes node-age ESS was low on these runs (32 to 35 per chain on the Fast FBD runs), so part of each difference is RevBayes sampling noise.

## 5. Fast FBD against full-DNA FBD (preliminary)

RevBayes 1.3.2 full-DNA FBD and `ffbd` Fast FBD on A06 and A12, two chains each, same fossils, priors and uniform root prior. The RevBayes collapse-expand move is removed (section 4), so the only remaining difference is the branch-length approximation. Each side was summarized as an extant time tree (median node ages), and nodes were matched by clade (`scripts/branch_scatter.py --ages-only`; `validation/node_ages_full_vs_fast.tsv`).

| Data set | Nodes | r | Median difference | Mean difference (Fast minus full) | Root age, full vs Fast |
|---|---|---|---|---|---|
| A06 | 23 | 0.998 | 1.5% | +0.5 Myr | 75.1 vs 75.2 Ma |
| A12 | 21 | 0.999 | 2.7% | -0.3 Myr | 139.1 vs 139.0 Ma |

**Convergence (preliminary; longer full-DNA chains in progress).** Two chains each, first 10% discarded. Data: A06 24 living species, 10 fossils; A12 22 living species, 14 fossils; 16,107 aligned sites each.

| Set | Side | Generations or proposals per chain | Run time | Node-age ESS (median, lowest) | Root-age ESS | Chains differ (median, largest) |
|---|---|---|---|---|---|---|
| A06 | Full DNA, RevBayes | 250,000 generations | ~25 h | 19, 10 | 51-52 | 1.6%, 6.6% |
| A06 | Fast FBD, `ffbd` | 400 million proposals | ~2.5 min | 901, 756 | 764-812 | 0.26%, 1.1% |
| A12 | Full DNA, RevBayes | 250,000 generations | ~22 h | 35, 14 | 557-1,531 | 4.7%, 15.3% |
| A12 | Fast FBD, `ffbd` | 400 million proposals | ~2.3 min | 901, 778 | 901 | 0.79%, 3.7% |

The full-DNA node-age ESS values are low, so the agreement in the table above is within the noise of the full-DNA runs. Four more full-DNA chains per clade (250,000 generations each) are running; the test will be updated when they finish.


## 6. Catfish (close, not identical)

`ffbd` Fast FBD against the original full analysis of the same data: 553 nodes matched by clade, r = 0.995, Fast 1.8 Myr younger on average (median difference 10.2%, largest on young nodes). Fast FBD ages come from all logged samples of both chains. The original analysis used an uncorrelated lognormal clock in RevBayes with the collapse-expand move, so it is not a like-for-like comparison.

## 7. Convergence on large trees (passed)

Two chains each at 1,185 and 2,424 tips, 20 billion proposals each (about 2 hours on one core): effective sample sizes 500 to 800 for the main parameters, and at least 139 for node ages. Node-age medians agree between chains to 0.9% (median), and 95% of nodes agree within 5%.

## 8. Crown and stem rules (passed)

The rules restrict where a fossil may attach: only below the crown node, or never below it. They add no density term. Proposals that would break a rule are rejected, which leaves the target distribution as the unconstrained one restricted to the allowed states. Two checks:

1. **Exact conditioning test.** On the small example (`examples/example_crown_stem.txt`), the run with a crown rule (FOSSIL_1) and a stem rule (FOSSIL_3) was compared with unconstrained runs keeping only samples where both rules held. With 11 seeds per side and 100,001 samples each, the rules held in every sample of the constrained runs. All quantities agreed within 1.4 standard errors (`validation/crown_stem_test_11seeds.txt`). A crown-only test with 2 seeds also agreed. `--check` runs, which validate the full state after every proposal, raised no error.
2. **Regression.** With no crown or stem rules, the new version gives byte-identical log and tree files to the previous release on A06, A12, bears and catfish, with and without DNA (`validation/regression_previous_release.txt`).

## 9. Three programs on the same prior (ffbd agrees with MrBayes)

An independent referee for section 4. `scripts/make_fixed3.py` writes the same model for `ffbd`, MrBayes 3.2.7 and RevBayes 1.3.2: the A06 or A12 living-species tree fixed (in MrBayes, one partial constraint per living-species node, fossils free), the same clade constraints, uniform fossil age windows (MrBayes `calibrate`), speciation 0.10, extinction 0.04, fossil sampling 0.03 and the root age all **fixed**, and no data. With every rate fixed, the three programs must give the same node ages and the same number of sampled ancestors. Two chains per program and data set; comparison `scripts/compare_fixed3.py` (`validation/three_program_prior_A06_A12.txt`).

| Data set | Mean sampled ancestors, ffbd / MrBayes / RevBayes | Node ages, ffbd vs MrBayes (median, largest) | RevBayes vs MrBayes (median, largest) | MrBayes chain 1 vs chain 2 (median, largest) |
|---|---|---|---|---|
| A06 | 7.01 / 7.02 / 7.95 | 0.36%, 2.7% | 3.2%, 34% | 0.89%, 7.4% |
| A12 | 9.55 / 9.58 / 10.13 | 0.59%, 4.2% | 6.5%, 28.5% | 0.86%, 6.4% |

`ffbd` and MrBayes agree within chain noise. RevBayes gives about one extra sampled ancestor and older young nodes (largest on two-tip nodes, for example 10.6 Ma in MrBayes and 14.2 Ma in RevBayes), the pattern expected from the move described in section 4. BEAST2 (SA package) was not used: its monophyly constraints count fossils as clade members, so it cannot hold a living-species tree fixed while fossils move.

## 10. Simulation study at realistic size

After Heath et al. (2014). `scripts/sim_study.py simulate` grows FBD trees forward from a known root (100 Ma; speciation 0.10, extinction 0.06, fossil sampling 0.01, extant sampling 0.5; trees kept at 30 to 120 living tips with at least 5 fossils), gives each fossil an age window of 2 to 8 Myr containing its true age and a clade constraint that is true on the true tree (the total group of the clade it attaches to, or of its parent or grandparent clade), and draws branch rates from an exponential with mean 0.004. Two versions of every tree:

- **exact**: branch lengths drawn from the `ffbd` likelihood itself (tests the sampler).
- **dna**: 2,000 sites simulated with AliSim (GTR+G4); ML branch lengths and 100 bootstrap standard errors on the true topology, exactly as for the real clades (`scripts/fish_fast.py`; tests the Fast FBD approximation).

Root prior soft-bounded on 90 to 110 Ma; 200 million proposals per run; 100 trees. Scoring: `scripts/sim_study_summary.py` (`validation/sim_study_100reps.txt`). Coverage of the 95% interval:

| Version | Node ages (6,911) | Fossil ages (3,631) | Root | Speciation | Extinction | Fossil sampling | Sampled ancestors | Mean clock rate |
|---|---|---|---|---|---|---|---|---|
| exact | 0.931 | 0.956 | 0.99 | 0.58 | 0.60 | 0.89 | 0.83 | 0.87 |
| dna | 0.914 | 0.954 | 0.99 | 0.56 | 0.54 | 0.92 | 0.82 | 0.38 |

- Node and fossil ages are close to 95% coverage, and the DNA approximation costs only 1 to 2 points.
- Nodes younger than 5 Ma come out too old in both versions (median error +56% in exact and +59% in dna, about 1.2 Myr). The error falls with age: +10 to 12% at 15 to 30 Ma, +1 to 2% above 60 Ma. Because both versions show it, this comes from the model and prior, not from the approximation.
- Speciation and extinction are under-covered and the dna version overestimates the mean clock rate. Section 11 shows the sampler itself is calibrated, so this comes from the design of this study: one fixed set of true values, trees kept only at 30 to 120 tips, and exponential priors that pull the rates down.

## 11. Simulation-based calibration at realistic size (passed)

`scripts/sim_study.py sbc`: each of 500 replicates draws speciation, extinction, fossil sampling and the mean clock rate from the `ffbd` priors (exponential, rate 10) and the root from the analysis prior Normal(30, 2), then simulates an FBD tree with no size filter. Only information that is exact for this test is used: fossil windows on a fixed 4-Myr grid, one constraint per side of the root, and a constant observation standard deviation (0.01). Trees had 2 to 8,056 living tips (median 20) and 0 to 15,838 fossils (median 25). 100 million proposals per replicate, 99 thinned draws, rank of the true value among them (`scripts/sbc_ranks.py`; `validation/sbc_realistic_500reps.txt`). With a correct sampler the ranks are uniform. Bonferroni threshold for 9 tests: p < 0.0056.

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

## Speed

0.39 seconds per million proposals on one core of a compute cluster, about 2.6 million proposals per second, and 0.13 to 0.16 seconds on a laptop. Speed was the same at 1,185 and 2,424 tips.
