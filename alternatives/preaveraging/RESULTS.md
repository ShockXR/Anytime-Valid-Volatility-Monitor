# Results: pre-averaging before MAD truncation

All numbers come from `run_experiment.py` (seeds fixed) and are stored in `results/`. Everything here is [SYNTHETIC]: one-second returns, 23,400 per day, daily volatility 1%, price-level noise with noise-to-signal ratio `gamma` (noise sd divided by the per-second return sd), 200 paths per cell unless stated. Glossary: "raw" means the MAD-clipped one-second returns with no pre-averaging (the normalised, kappa-corrected version of Bapat's estimator); "pre-averaged" means window `k = 60` with weight `min(x, 1 - x)`; "clip" and "drop" are the two ways of truncating a pre-averaged return at `3 x scaled MAD`; "blocks" are every `k`-th pre-averaged return, which use disjoint prices.

## Q1. When does the plain truncated estimator break, and does pre-averaging fix it? [SYNTHETIC]

Integrated quarticity, mean ratio to the true value (relative RMSE in brackets), log-OU stochastic volatility, no jumps:

| gamma | 0 | 0.1 | 0.25 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|---|
| theory for raw, `(1 + 2 gamma^2)^2` | 1 | 1.04 | 1.27 | 2.25 | 9 | 81 | 1089 |
| raw returns, MAD clip | 0.99 (0.02) | 1.03 (0.04) | 1.26 (0.26) | 2.23 (1.24) | 9.0 (8.1) | 81.4 (81.5) | 1138 (1155) |
| pre-averaged, no truncation | 1.00 (0.12) | 0.99 (0.12) | 1.01 (0.13) | 1.01 (0.13) | 1.00 (0.10) | 1.01 (0.13) | 1.00 (0.13) |
| pre-averaged, MAD clip | 0.99 (0.11) | 0.99 (0.11) | 1.00 (0.11) | 1.00 (0.12) | 0.99 (0.10) | 0.99 (0.11) | 0.99 (0.11) |
| pre-averaged, MAD clip, no noise correction | 0.99 (0.11) | 0.99 (0.11) | 1.00 (0.11) | 1.00 (0.12) | 1.00 (0.10) | 1.01 (0.12) | 1.10 (0.15) |
| pre-averaged, MAD drop | 0.98 (0.11) | 0.99 (0.11) | 0.98 (0.11) | 0.99 (0.12) | 0.98 (0.11) | 0.97 (0.12) | 0.97 (0.12) |

The same with 4 jumps per day (size 50 per-second sd):

| gamma | 0 | 0.1 | 0.25 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|---|
| RQ, no truncation | 1047 (2246) | 1147 (2648) | 1134 (2194) | 844 (1566) | 1154 (2139) | 1336 (2413) | 2143 (2612) |
| raw returns, MAD clip | 1.00 (0.02) | 1.04 (0.04) | 1.26 (0.26) | 2.25 (1.26) | 9.0 (8.1) | 80.4 (80.5) | 1083 (1097) |
| pre-averaged, no truncation | 33.4 (69) | 34.6 (74) | 36.5 (68) | 27.0 (47) | 35.1 (63) | 39.3 (71) | 32.5 (52) |
| pre-averaged, MAD clip | 1.14 (0.22) | 1.12 (0.19) | 1.15 (0.22) | 1.14 (0.20) | 1.17 (0.22) | 1.18 (0.25) | 1.15 (0.24) |
| pre-averaged, MAD drop | 0.99 (0.12) | 0.99 (0.11) | 1.00 (0.13) | 1.00 (0.11) | 1.01 (0.12) | 1.02 (0.12) | 0.99 (0.12) |

Integrated variance behaves the same way with the square root of the factors: the raw truncated IV has mean ratio `1 + 2 gamma^2` (1.13 at `gamma = 0.25`, 3.0 at 1, 33.9 at 4), the pre-averaged versions 0.99 to 1.00 without jumps (relative RMSE 0.05 against 0.01 for raw at `gamma = 0`), clip 1.04 to 1.06 and drop 0.99 to 1.01 with jumps.

Window sensitivity (`results/q1_window_sensitivity.csv`, `gamma = 1`, 4 jumps): mean ratio of the clipped version grows with `k` (1.06, 1.11, 1.15, 1.24, 1.32 for `k = 20, 40, 60, 100, 150`) because a jump stays in `k - 1` windows; the dropped version stays at 0.99 to 1.07 with relative RMSE 0.07, 0.09, 0.11, 0.16, 0.20. At `gamma = 4` with `k = 20` the noise term is 47% of `sigma^2` before the correction, and the corrected drop estimator still has mean ratio 1.00 (RMSE 0.09).

**Reading.**

* **The plain truncated estimator breaks between `gamma = 0.1` and `0.25`:** 3% bias at 0.1, 26% at 0.25, a factor 2.2 at 0.5, 9 at 1, 81 at 2, exactly `(1 + 2 gamma^2)^2`. The 10% line is crossed at `gamma = 0.16` and the factor 2 at `gamma = 0.45`. Truncation does nothing against noise, because noise is in every return, not in a few.
* **Pre-averaging then truncation fixes it at every noise level**, with mean ratio 0.99 to 1.00 and relative RMSE 0.10 to 0.12. The price is efficiency when there is no noise: the `n^{1/4}` rate gives a relative RMSE of 0.11 against 0.02 for the raw estimator at `gamma = 0`. Pre-averaging should be switched on when the estimated noise ratio exceeds about 0.15, not always.
* **The noise correction for the fourth moment matters from `gamma` of about 2** (1% at 2, 10% at 4 for `k = 60`; 47% at `k = 20` and `gamma = 4`). It is cheap and should always be applied.
* **Clipping is the wrong truncation after pre-averaging when there are jumps.** A clipped jump keeps `tau^4` in each of the `k - 1` windows it enters, about 4% of the day's quarticity per jump at `k = 60`: +12% to +18% with 4 jumps per day, growing with `k`. Dropping the window instead leaves 0.99 to 1.02. On raw returns clipping is fine (a jump touches one return), which is why the main report found no such bias.
* **Pre-averaging without truncation is not jump robust at all** (mean ratio 27 to 39 with 4 jumps per day): the two stages are complementary, as the entry in `alternatives/ALTERNATIVES.md` assumed.

## Q2. The local variance used by the monitor [SYNTHETIC]

Trailing 300-second window, evaluated every 600 seconds, 100 paths, against the true spot variance. Mean ratio and median absolute log error (a typical relative error):

| gamma | 0 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| raw returns, MAD clip, no jumps | 1.00 (0.06) | 1.51 (0.41) | 3.07 (1.11) | 9.23 (2.21) | 33.6 (3.50) |
| pre-averaged, overlapping, no jumps | 1.00 (0.33) | 1.01 (0.31) | 0.98 (0.32) | 1.00 (0.32) | 0.99 (0.31) |
| pre-averaged, blocks, no jumps | 1.01 (0.46) | 1.02 (0.43) | 0.98 (0.45) | 1.00 (0.45) | 1.00 (0.47) |
| raw returns, MAD clip, 4 jumps | 1.00 (0.06) | 1.51 (0.41) | 3.04 (1.10) | 9.24 (2.21) | 33.7 (3.50) |
| pre-averaged, overlapping, 4 jumps | 1.05 (0.32) | 1.05 (0.33) | 1.05 (0.33) | 1.06 (0.33) | 1.06 (0.32) |
| pre-averaged, blocks, 4 jumps | 1.04 (0.46) | 1.05 (0.44) | 1.07 (0.48) | 1.04 (0.48) | 1.07 (0.47) |

**Reading.** The raw local variance is inflated by `1 + 2 gamma^2` (1.5 at `gamma = 0.5`, 3 at 1), with the same breaking point as Q1 (10% at `gamma = 0.22`). The pre-averaged local variance is unbiased at every noise level but much noisier over a 300-second window: a typical error of 32% (overlapping) or 45% (blocks) against 6% for raw returns without noise, because a window of 300 seconds holds only about 5 independent pre-averaged blocks of 60 seconds. For the monitor this is the real trade-off: pre-averaging buys unbiasedness under noise with a loss of about 60 in effective sample size at `k = 60`.

## Q3. The monitor's false alarm rate under noise [SYNTHETIC]

Constant volatility, two calibration days, one monitored day, each monitor at `A = steps per day / 0.05` (bound: 5% false alarm per day). One-day false alarm rate, increase side / decrease side:

| gamma | 0 | 0.5 | 1 | 2 | 4 | Conditional null holds under the model |
|---|---|---|---|---|---|---|
| raw returns (23,400 steps, `A = 468,000`) | 5.0% / 7.5% | 8.5% / 9.0% | **18.0% / 13.5%** | **31.0% / 29.0%** | **41.5% / 27.5%** | no under price noise |
| overlapping pre-averaged, `k = 60` (23,341 steps) | **100% / 100%** | 100% / 100% | 100% / 100% | 100% / 100% | 100% / 100% | no |
| blocks, `k = 60` (390 steps, `A = 7800`) | 2.5% / 0.0% | 1.5% / 0.5% | 2.5% / 1.0% | 2.0% / 2.0% | 4.0% / 1.5% | yes |
| blocks, `k = 20` (1170 steps, `A = 23,400`) | 4.5% / 3.0% | 3.0% / 5.0% | 2.5% / 1.0% | 5.5% / 3.5% | 3.0% / 2.5% | yes |
| blocks, `k = 10` (2340 steps, `A = 46,800`) | 5.5% / 1.5% | 5.0% / 4.0% | 5.5% / 2.0% | 2.5% / 2.0% | 3.0% / 4.5% | yes |

**Reading.**

* **Overlapping pre-averaged returns destroy the guarantee completely** (100% false alarms per day at every noise level, including zero noise). Consecutive windows share `k - 2` raw returns, so after a large value the next value is expected to be large: the conditional-mean null `E[Y_t | past] <= 1` fails even though the unconditional mean is exactly 1. Smoothing a statistic before betting on it is not allowed.
* **The raw-return monitor also loses its guarantee under price noise, and this was not known before.** Its false alarm rate is 5.0% at `gamma = 0` (at the bound, with a small excess on the decrease side from the estimated baseline) and rises to 8.5%, 18%, 31% and 41.5% at `gamma = 0.5, 1, 2, 4`. The reason is the same as for overlap: with price noise `r_t = e_t + u_t - u_{t-1}`, so `E[r_t^2 | past] = sigma^2 delta + omega^2 + E[u_{t-1}^2 | past]`, which exceeds the baseline after a large return. Consecutive squared returns are positively correlated (correlation `rho^2` with `rho = -gamma^2 / (1 + 2 gamma^2)`, 0.2 at `gamma = 2`), and the betting process is exploited by it. E4 of the main report saw only 3.7% at `gamma = 0.5` over 390 steps; over 23,400 steps the excess is visible at the same noise level (8.5%). The statement in `docs/math.md` section 5.1 that the null "allows any dependence" is true only of dependence that keeps the *conditional* moment on the null side; MA(1) noise does not.
* **Non-overlapping blocks restore the guarantee at every `k` tested** (0% to 5.5%, within Monte Carlo error of the bound), because blocks use disjoint sets of prices and are independent under i.i.d. noise. Dropping one return between blocks is what makes this work; adjacent blocks that share a price would again be dependent.

## Q4. Detection with noise: raw monitor against block monitors [SYNTHETIC]

Sd change at the start of the monitored day. Median delay in seconds (every monitor detected on 100% of paths except where stated):

| sd change | monitor | gamma = 0 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|
| x 1.5 | raw returns | **72.5** | **106.5** | 282.5 | 1566.5 | 8077 (94% detected) |
| x 1.5 | blocks, `k = 60` | 2699 | 2759 | 2579 | 2879 | 2729 |
| x 1.5 | blocks, `k = 20` | 949 | 979 | 1069 | 1199 | 1589 |
| x 1.5 | blocks, `k = 10` | 574 | 554 | **599** | **899** | **2109** |
| x 2 | raw returns | **32** | **43.5** | 86.5 | 449.5 | 2907 |
| x 2 | blocks, `k = 60` | 1319 | 1379 | 1379 | 1319 | 1379 |
| x 2 | blocks, `k = 20` | 519 | 479 | 529 | 579 | 679 |
| x 2 | blocks, `k = 10` | 259 | 269 | **299** | **359** | **649** |
| x 0.5 | raw returns | **147** | **218** | 517 | 2824.5 | 10,795 (70% detected) |
| x 0.5 | blocks, `k = 60` | 6719 | 6359 | 6599 | 6599 | 6779 |
| x 0.5 | blocks, `k = 20` | 2319 | 2359 | 2379 | 2579 | 3389 |
| x 0.5 | blocks, `k = 10` | 1244 | 1259 | **1369** | **1769** | **3939** |

Dilution check: the mean of the raw monitor's normalised stream after a doubling of the sd was 3.13, 2.59, 1.89, 1.32 and 1.09 at `gamma = 0, 0.5, 1, 2, 4`, against the theory `(4 + 2 gamma^2) / (1 + 2 gamma^2)` = 4 (3.13 after clipping), 2.67, 2.0, 1.33, 1.09. The `k = 60` blocks stayed at 3.04 to 3.15 at every noise level; the `k = 10` blocks fell to 2.61 and 1.92 at `gamma = 2` and 4, because their own noise term (`12 gamma^2 / k^2` in units of the signal) is no longer small.

**Reading.**

* **Without noise the raw monitor is 8 to 45 times faster** than the block monitors (32 s against 259 s for `k = 10` and 1319 s for `k = 60` at a doubling), simply because it keeps every observation.
* **Noise dilutes the raw monitor's signal** exactly as the formula predicts, so its delay grows 90 times between `gamma = 0` and 4 for a doubling, and a halving at `gamma = 4` is missed on 30% of days. The block monitors are almost unaffected up to `gamma = 2`.
* **The crossover is at `gamma` of about 1 to 2.** With `k = 10` the block monitor is valid at every noise level (Q3) and from `gamma = 2` on it is also faster than the raw monitor (359 s against 449.5 s at a doubling), which by then has a 31% false alarm rate anyway. Between `gamma = 1` and 2 the raw monitor is still faster but no longer valid. The practical rule: estimate `gamma` from the first autocovariance; below about 0.5 use raw returns, above it use blocks with a small `k` (10 to 20) and the noise correction.
* **Smaller blocks are better for the monitor than the estimation-optimal `k = 60`** (`theta = 0.4`), since the monitor cares about the number of observations, not the `n^{1/4}` constant; `k = 10` is 5 times faster than `k = 60` and still valid, at the cost of a noise term that must be corrected for and that dilutes the signal at `gamma >= 2`.

## Answers to the questions

1. **Breaking point of the plain truncated estimator:** about `gamma = 0.16` for a 10% bias in integrated quarticity (`(1 + 2 gamma^2)^2`), 0.22 for the local variance (`1 + 2 gamma^2`); a factor 2 at `gamma = 0.45`; the estimator is useless by `gamma = 1`. **Pre-averaging then truncation fixes both** at every noise level (mean ratio 0.99 to 1.00), with the noise correction for the fourth moment needed from `gamma` of about 2, and with "drop" rather than "clip" when there are jumps.
2. **The monitor's false alarm rate:** on raw returns it rises from 5% to 41.5% per day as `gamma` goes from 0 to 4 (price noise breaks the conditional null); on overlapping pre-averaged returns it is 100%; on non-overlapping blocks it stays within the 5% bound at every `k` and `gamma`. The block version is justified because blocks use disjoint prices and are therefore independent under the model; it costs a factor `k` in observations, which `k = 10` to 20 keeps affordable.

## Negative results and failed attempts

* Clipping (the paper's truncation) after pre-averaging has a jump bias that grows with `k`; it had to be replaced by dropping.
* The hypothesis that the raw-return monitor stays valid under price noise, which the main report's E4 supported at one noise level over one short day, is false over a full day at `gamma >= 0.5`.
* A first unit test compared a calibration path and monitored paths simulated with different numbers of steps, so their per-step scales differed by `sqrt(10)` and the detector "false alarmed" on every path. The detector was right; the test was fixed.
* A background watcher for the first run never exited because `pgrep -f` matched its own command line; the run itself had finished. No result was affected.

## Ranking

* **Research novelty:** moderate to high. The composition (pre-average, drop, noise-correct, then bet on blocks) is new as a monitor, and the finding that MA(1) price noise breaks the e-detector's conditional null on raw returns is a correction to the main package's claims.
* **Mathematical validity:** high for the estimator (unbiasedness verified to 1% across 14 cells) and for the block monitor (bound holds in 30 cells); the local-constancy derivation of the fourth-moment correction is heuristic and should be replaced by the exact pre-averaged quarticity of Jacod et al. (2009) in a paper.
* **Practical relevance:** high at one-second sampling, where noise ratios of 0.2 to 2 are common; low at one minute for BTC (noise ratio well below the 0.16 breaking point).
* **Feasibility:** high. About 250 lines, no new dependencies, five minutes to run.
* **Demonstration value:** high. Both figures make the point in one glance (RMSE flat against a line that rises by a factor 1000; false alarms rising against flat).
* **Likely interest to a sequential-analysis and high-frequency econometrics researcher:** high. The overlap and MA(1) failures are exactly the kind of assumption that e-detector users must check, and the pre-averaged block monitor is a clean, provable fix.
