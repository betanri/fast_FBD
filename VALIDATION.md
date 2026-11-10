# Validation of ffbd

Status as of 10 November 2026. "RevBayes" means RevBayes 1.3.2. Result files are in `validation/`.

Data sets: two small fish clades (A06: 24 living species and 10 fossils; A12: 22 and 14), the RevBayes bear tutorial data, and a catfish data set (713 living tips, 94 fossils).

## Summary

| # | Test | Result |
|---|---|---|
| 1 | Probabilities equal RevBayes on fixed trees | Passed |
| 2 | Simulation-based calibration and analytic prior | Passed |
| 3 | Agreement with RevBayes, sampled-ancestor changes off | Passed |
| 4 | Agreement with RevBayes, sampled-ancestor changes on | Does not match (section 4) |
| 5 | Fast FBD against full-DNA FBD | Passed |
| 6 | Catfish comparison | Close, not identical (section 6) |
| 7 | Convergence on large trees | Passed |
| 8 | Crown and stem rules | Passed |

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

## 5. Fast FBD against full-DNA FBD (passed)

RevBayes 1.3.2 full-DNA FBD and `ffbd` Fast FBD on A06 and A12, two chains each, same fossils, priors and uniform root prior. The RevBayes collapse-expand move is removed (section 4), so the only remaining difference is the branch-length approximation. Each side was summarized as an extant time tree (median node ages), and branches were matched by clade and parent clade (`scripts/branch_scatter.py`; `validation/branch_lengths_full_vs_fast.tsv`).

| Data set | Branches | r (log scale) | Median difference | Within 10% | Root age, full vs Fast |
|---|---|---|---|---|---|
| A06 | 46 | 0.998 | 2.9% | 91% | 75.1 vs 75.2 Ma |
| A12 | 42 | 0.996 | 2.9% | 91% | 139.1 vs 139.0 Ma |

## 6. Catfish (close, not identical)

`ffbd` Fast FBD against the original full analysis of the same data: 1,106 branches, r = 0.958 (log scale), median difference 10.6%, 46% within 10%; node ages r = 0.994, Fast 1.8 Myr younger on average. Fast FBD ages come from all logged samples of both chains. The original analysis used an uncorrelated lognormal clock in RevBayes with the collapse-expand move, so it is not a like-for-like comparison. The differences are largest on short branches.

## 7. Convergence on large trees (passed)

Two chains each at 1,185 and 2,424 tips, 20 billion proposals each (about 2 hours on one core): effective sample sizes 500 to 800 for the main parameters, and at least 139 for node ages. Node-age medians agree between chains to 0.9% (median), and 95% of nodes agree within 5%.

## 8. Crown and stem rules (passed)

The rules restrict where a fossil may attach: only below the crown node, or never below it. They add no density term. Proposals that would break a rule are rejected, which leaves the target distribution as the unconstrained one restricted to the allowed states. Two checks:

1. **Exact conditioning test.** On the small example (`examples/example_crown_stem.txt`), the run with a crown rule (FOSSIL_1) and a stem rule (FOSSIL_3) was compared with unconstrained runs keeping only samples where both rules held. With 11 seeds per side and 100,001 samples each, the rules held in every sample of the constrained runs. All quantities agreed within 1.4 standard errors (`validation/crown_stem_test_11seeds.txt`). A crown-only test with 2 seeds also agreed. `--check` runs, which validate the full state after every proposal, raised no error.
2. **Regression.** With no crown or stem rules, the new version gives byte-identical log and tree files to the previous release on A06, A12, bears and catfish, with and without DNA (`validation/regression_previous_release.txt`).

## Speed

0.39 seconds per million proposals on one core of a compute cluster, about 2.6 million proposals per second, and 0.13 to 0.16 seconds on a laptop. Speed was the same at 1,185 and 2,424 tips.
