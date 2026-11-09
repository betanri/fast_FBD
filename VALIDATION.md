# Validation of ffbd

Status as of 9 November 2026. "RevBayes" means RevBayes 1.3.2. Result files are in `validation/`.

Data sets: two small fish clades (A06: 24 living species and 10 fossils; A12: 22 and 14), the RevBayes bear tutorial data, and a catfish data set (713 living tips, 94 fossils).

## Summary

| # | Test | Result |
|---|---|---|
| 1 | Probabilities equal RevBayes on fixed trees | Passed |
| 2 | Simulation-based calibration and analytic prior | Passed |
| 3 | Agreement with RevBayes, sampled-ancestor changes off | Passed |
| 4 | Agreement with RevBayes, sampled-ancestor changes on | Does not match (section 4) |
| 5 | Fast FBD against full-DNA FBD | Running |
| 6 | Catfish comparison | Running |
| 7 | Convergence on large trees | Running |

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

## 5 to 7. Running

- **Fast FBD against full DNA.** RevBayes full-DNA and `ffbd` runs on A06 and A12, with the collapse-expand move off in both, so the only remaining difference is the branch-length approximation.
- **Catfish.** Comparison with RevBayes on the same data.
- **Large trees.** Convergence with two or more chains at 1,185 and 2,424 tips. Interim effective sample sizes after 4.7 billion proposals were 134 to 192 for the main parameters.

## Speed

0.39 seconds per million proposals on one core of a compute cluster, about 2.6 million proposals per second, and 0.13 to 0.16 seconds on a laptop. Speed was the same at 1,185 and 2,424 tips.
