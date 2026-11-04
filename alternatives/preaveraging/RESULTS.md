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
* **The crossover is at `gamma` of about 1 to 2.** With `k = 10` the block monitor is valid at every noise level (Q3) and from `gamma = 2` on it is also faster than the raw monitor (359 s against 449.5 s at a doubling), which by then has a 31% false alarm rate anyway. Between `gamma = 1` and 2 the raw monitor is still faster but no longer valid. The practical rule (sharpened by Q6): estimate `gamma` from the first autocovariance; if it is clearly above zero the raw monitor is already invalid (8.5% at `gamma = 0.5`), so use blocks with `k = max(4, ceil(4 gamma))` and the noise correction. With `k = 4` the blocks detect a doubling in 127 s at `gamma = 0.5`, a factor 3 behind the invalid raw monitor rather than the factor 8 of `k = 10`.
* **Smaller blocks are better for the monitor than the estimation-optimal `k = 60`** (`theta = 0.4`), since the monitor cares about the number of observations, not the `n^{1/4}` constant; `k = 10` is 5 times faster than `k = 60` and still valid, at the cost of a noise term that must be corrected for and that dilutes the signal at `gamma >= 2`.

## Q5. Noise and jumps together: benchmark against Podolskij and Vetter (2009) [SYNTHETIC]

`results/q5_benchmark_noise_jumps.csv`, `run_benchmark.py`. 200 paths per cell, log-OU stochastic volatility, `gamma` in `{0, 0.5, 1, 2}`, jumps per day x size in `{0, 2 x 25 sd, 2 x 100 sd, 8 x 25 sd, 8 x 100 sd}` (25 per-second sd is 0.16% of the price, 100 sd is 0.65%). The modulated bipower (variance) and tripower (quarticity) variation uses `K = 76`-step increments averaged over blocks of 153 observations (`c1 = 0.5`, `c2 = 1`, 152 blocks per day) with the noise correction through the exact block weights. "Ours" is pre-averaging (`k = 60`) then MAD drop with the noise correction.

**Integrated quarticity, mean ratio to the truth (relative RMSE) [SYNTHETIC].**

| cell | raw returns, MAD clip | pre-averaged, MAD drop (ours) | pre-averaged, MAD clip | modulated tripower (PV 2009) | pre-averaged tripower | pre-averaged, no truncation |
|---|---|---|---|---|---|---|
| no jumps, gamma = 0 | 0.99 (0.02) | 0.98 (0.11) | 0.99 (0.10) | 1.00 (0.30) | 1.00 (0.13) | 1.00 (0.12) |
| no jumps, gamma = 0.5 | 2.22 (1.23) | 0.98 (0.12) | 0.99 (0.12) | 1.01 (0.33) | 1.00 (0.14) | 1.00 (0.12) |
| no jumps, gamma = 1 | 8.91 (8.00) | 0.98 (0.11) | 0.99 (0.11) | 0.99 (0.33) | 1.00 (0.13) | 1.00 (0.12) |
| no jumps, gamma = 2 | 80.9 (81.1) | 0.98 (0.11) | 0.99 (0.11) | 1.00 (0.30) | 1.00 (0.12) | 1.00 (0.11) |
| 2 x 25 sd, gamma = 0 | 0.99 (0.02) | 0.99 (0.10) | 1.03 (0.11) | 1.12 (0.37) | 1.06 (0.15) | 2.20 (2.8) |
| 2 x 25 sd, gamma = 1 | 8.98 (8.06) | 0.98 (0.11) | 1.03 (0.12) | 1.03 (0.34) | 1.04 (0.16) | 2.12 (2.8) |
| 2 x 100 sd, gamma = 0 | 0.99 (0.02) | 0.99 (0.12) | 1.10 (0.19) | 1.63 (1.27) | 1.45 (0.86) | 269 (589) |
| 2 x 100 sd, gamma = 1 | 9.09 (8.17) | 0.98 (0.11) | 1.09 (0.16) | 1.57 (1.73) | 1.56 (1.63) | 289 (690) |
| 8 x 25 sd, gamma = 0 | 1.00 (0.02) | 1.02 (0.11) | 1.18 (0.24) | 1.27 (0.57) | 1.21 (0.33) | 5.57 (7.1) |
| 8 x 25 sd, gamma = 1 | 9.13 (8.22) | 1.05 (0.12) | 1.22 (0.27) | 1.35 (0.65) | 1.25 (0.34) | 5.36 (6.5) |
| 8 x 100 sd, gamma = 0 | 1.00 (0.02) | 1.00 (0.11) | 1.43 (0.49) | 4.94 (10.3) | 3.96 (7.8) | 973 (1414) |
| 8 x 100 sd, gamma = 1 | 9.06 (8.14) | 1.01 (0.12) | 1.43 (0.49) | 5.76 (13.8) | 3.58 (4.8) | 1200 (1811) |
| 8 x 100 sd, gamma = 2 | 81.2 (81.3) | 1.01 (0.12) | 1.43 (0.49) | 4.53 (6.9) | 3.77 (5.2) | 1013 (1604) |

The cells not shown (`gamma = 0.5` and the remaining `gamma = 2` rows) are in the CSV and fall between their neighbours; the only outlier is modulated tripower at `8 x 100 sd, gamma = 0.5`: 13.7 (117), one path with a block product of extreme size.

**Integrated variance [SYNTHETIC].** Ours: mean ratio 0.99 to 1.02 in all 20 cells, relative RMSE 0.047 to 0.061. Modulated bipower: 0.99 to 1.01 (RMSE 0.13 to 0.14) without jumps, 1.05 (0.14) with `2 x 25 sd`, 1.10 (0.19) with `8 x 25 sd`, 1.16 to 1.20 (0.31 to 0.33) with `2 x 100 sd`, 1.72 to 1.83 (0.94 to 1.27) with `8 x 100 sd`, at every noise level. Raw MAD clip: 1.50, 3.0, 9.0 at `gamma = 0.5, 1, 2`, jumps or not.

**CLT interval for the modulated bipower variance (95%, standard error from the independent Gaussian block approximation) [SYNTHETIC].** Coverage without jumps 0.91, 0.935, 0.905, 0.955 at `gamma = 0, 0.5, 1, 2` (nominal 0.95, Monte Carlo sd 0.015); with `2 x 25 sd` 0.945 to 0.97; with `8 x 25 sd` 0.91 to 0.955; with `2 x 100 sd` 0.905 to 0.94; with `8 x 100 sd` 0.57 to 0.62 (the jump bias is outside the interval). Median relative half-width 0.26 without jumps, 0.61 to 0.67 with `8 x 100 sd`. Coverage is 1 to 4 points below nominal without jumps, which is the known small-sample effect of a Gaussian block approximation with 152 blocks and a `sqrt(n)`-dependent bias; the paper's own constants could not be checked (full text unreachable here), so the interval is ours.

**Reading.**

* **Without jumps all noise-robust estimators are unbiased, and the modulated multipower is the least precise:** relative RMSE 0.30 for quarticity and 0.13 to 0.14 for variance against 0.11 and 0.05 to 0.06 for pre-averaging then truncation, a factor 2.5 to 3. The reason is structural: `K = c1 sqrt(n)` block averaging reduces the day to 152 nearly independent numbers, while pre-averaging with `k = 60` keeps 23,341 overlapping ones, equivalent to about 390 independent blocks of `k` with a variance-reducing overlap; and powers below 2 (bipower, tripower) cost efficiency even without jumps.
* **With jumps the multipower estimators are robust in the limit but not at `n = 23,400`.** A jump of size `J` enters one block average with weight close to 1, and the bipower / tripower products then keep `|J|` or `|J|^{4/3}` times the neighbouring block's average, which vanishes only as `n^{-1/4}`. Eight jumps of 100 sd (0.65% each, 5% of the day's variance) inflate the modulated tripower quarticity by a factor 4.5 to 5.8 and the modulated bipower variance by a factor 1.7 to 1.8, and the pre-averaged tripower by 3.6 to 4.0 and 1.5. Two jumps of 25 sd leave 3% to 12%. Truncation by dropping after pre-averaging removes all of this: ours is 0.98 to 1.05 across the whole grid.
* **Clipping after pre-averaging is the worst of the truncation rules at many jumps** (1.43 at `8 x 100 sd`), confirming Q1, but still far better than any multipower rule.
* **The raw MAD-clipped estimator is excellent against jumps and useless against noise** (1.00 with 8 jumps of 100 sd at `gamma = 0`; 2.2, 9, 81 at `gamma = 0.5, 1, 2`), which is the gap Bapat (2026) leaves and pre-averaging fills.
* **Noise and jumps do not interact** for ours or for the multipower estimators: every column is flat in `gamma` within Monte Carlo error. The two problems are handled by two separate operations (pre-averaging for noise, truncation for jumps) and the operations compose.

## Q6. The block size: delay against `k`, validity at small `k`, and a rule of thumb [SYNTHETIC]

`results/q6_block_size.csv`, `results/fig3_block_size_rule.png`. Constant volatility, 200 paths per cell, two calibration days then one monitored day, blocks `k` in `{4, 6, 10, 14, 20, 30, 60}` (5850 to 390 blocks per day), `A = (n / k) / 0.05`, `gamma` in `{0.5, 1, 2, 4}`.

**Median detection delay in seconds, sd x 2 at the start of the day [SYNTHETIC].**

| gamma | k = 4 | 6 | 10 | 14 | 20 | 30 | 60 | best k observed | `k* = 6 gamma / v` |
|---|---|---|---|---|---|---|---|---|---|
| 0.5 | **127** | 179 | 269 | 370 | 479 | 749 | 1319 | 4 | 1.5 (floored to 4) |
| 1 | **179** | 206 | 289 | 377 | 539 | 719 | 1439 | 4 | 3 (floored to 4) |
| 2 | 405 | **332** | 354 | 419 | 539 | 764 | 1439 | 6 | 6 |
| 4 | 2131 | 1154 | 689 | **615** | 659 | 869 | 1439 | 14 | 12 |

**Median detection delay in seconds, sd x 1.5 [SYNTHETIC].**

| gamma | k = 4 | 6 | 10 | 14 | 20 | 30 | 60 | best k observed | `k* = 6 gamma / v` |
|---|---|---|---|---|---|---|---|---|---|
| 0.5 | **299** | 398 | 589 | 790 | 999 | 1514 | 2879 | 4 | 2 (floored to 4) |
| 1 | **459** | 485 | 614 | 783 | 1039 | 1499 | 2699 | 4 | 4 |
| 2 | 1399 | 965 | **839** | 944 | 1089 | 1559 | 2699 | 10 | 8 |
| 4 | 8751 | 4145 | 2329 | 1672 | **1539** | 1754 | 2669 | 20 | 16 |

Every change was detected within the day in every cell except `gamma = 4`, `k = 4`, sd x 1.5 (97.5%).

**One-day false alarm rate (bound 5%) [SYNTHETIC].** 28 cells (4 noise levels x 7 block sizes): 0.5% to 6.5%, mean 3.4%; four cells are above 5% by at most 1.5 points (6.5% at `gamma = 0.5, k = 10`; 5.5% and 6.0% elsewhere), as expected from a 5% rate measured on 200 paths (Monte Carlo sd 1.5 points). The guarantee holds down to `k = 4`, where the noise share of a block is 0.17 to 10.7 and the block statistic is far from Gaussian; nothing in the derivation needs Gaussian blocks.

**Observed mean of the monitored value after the change against the dilution formula [SYNTHETIC].** For a doubling the untruncated formula `(4 + rho_k) / (1 + rho_k)` gives 1.26 at `gamma = 4, k = 4` (`rho_k = 10.7`) and the observed mean is 1.25; at `k = 60` the observed mean is 3.05 to 3.14, which is the truncated value of 4 (`4 kappa_2(1.5) / kappa_2(3) = 3.12`) rather than 4 itself.

**Reading.**

* **The delay is U-shaped in `k` with a minimum at `k* = 6 gamma / v`**, as the delay formula `D(k) proportional to k (v^2 + rho_k)^2` predicts (README, section 2): the observed minima are at 4, 4, 6, 14 for a doubling and 4, 4, 10, 20 for `x 1.5`, against 4, 4, 6, 12 and 4, 4, 8, 16 from the formula, within one grid step everywhere. Right of the minimum the delay grows in proportion to `k` (a factor 10 from `k = 4` to 60 at `gamma = 0.5`); left of it the noise dilution dominates (a factor 3.5 from `k = 14` down to 4 at `gamma = 4`).
* **The rule of thumb `k = max(4, ceil(4 gamma))`** (target `v = 1.5`, `preavg.rule_of_thumb_k`) gives `k = 4, 4, 8, 16` at `gamma = 0.5, 1, 2, 4` and lands within 10% of the best observed delay in all 8 cells (interpolating between grid points: 343 s against 332 s at `gamma = 2` and 630 s against 615 s at `gamma = 4` for a doubling; 902 s against 839 s and 1628 s against 1539 s for `x 1.5`). The first rule tried, `k = 7 gamma` from capping the dilution at 0.25, was conservative by a factor 2 in `k` and cost 25% to 40% in delay; the derivation replaced it.
* **At small noise the smallest block is best, and it is fast:** `k = 4` detects a doubling in 127 s at `gamma = 0.5` and 179 s at `gamma = 1`, against 43.5 s and 86.5 s for the raw monitor, which is invalid there (8.5% and 18% false alarms). So the price of validity under noise is a factor 2 to 3 in delay at `gamma <= 1`, not the factor 8 of the `k = 10` blocks reported in Q4.
* **Delay cost of `k` on the right branch:** about `k / 4` times the `k = 4` delay, so `k = 60` costs a factor 10 over `k = 4` when the noise is small; this is why the estimation-optimal `theta = 0.4` is the wrong choice for monitoring.

## Answers to the questions

1. **Breaking point of the plain truncated estimator:** about `gamma = 0.16` for a 10% bias in integrated quarticity (`(1 + 2 gamma^2)^2`), 0.22 for the local variance (`1 + 2 gamma^2`); a factor 2 at `gamma = 0.45`; the estimator is useless by `gamma = 1`. **Pre-averaging then truncation fixes both** at every noise level (mean ratio 0.99 to 1.00), with the noise correction for the fourth moment needed from `gamma` of about 2, and with "drop" rather than "clip" when there are jumps.
2. **The monitor's false alarm rate:** on raw returns it rises from 5% to 41.5% per day as `gamma` goes from 0 to 4 (price noise breaks the conditional null); on overlapping pre-averaged returns it is 100%; on non-overlapping blocks it stays within the 5% bound at every `k` and `gamma`. The block version is justified because blocks use disjoint prices and are therefore independent under the model (README, section 2, formal statement: no boundary term, since one return between blocks is skipped); it costs a factor `k` in observations, which `k = max(4, ceil(4 gamma))` keeps affordable.
3. **Noise and jumps together (Q5):** pre-averaging then MAD drop is unbiased for integrated quarticity (0.98 to 1.05) and variance (0.99 to 1.02) over the whole grid of 4 noise levels x 5 jump settings, with relative RMSE 0.10 to 0.13 and 0.05 to 0.06. The modulated bipower / tripower variation of Podolskij and Vetter (2009) is unbiased without jumps but 2.5 to 3 times less precise (RMSE 0.30 and 0.13), and at `n = 23,400` it keeps a jump bias of `|J|^{4/3}` per affected product: a factor 4.5 to 5.8 in quarticity and 1.7 to 1.8 in variance with 8 jumps of 100 sd. Its CLT interval for variance covers 0.905 to 0.955 without jumps (nominal 0.95) and 0.57 to 0.62 with 8 large jumps.
4. **Block size (Q6):** the delay is U-shaped in `k` with its minimum at `k* = 6 gamma / v` (derived from the growth rate of the bet and confirmed within one grid step at every noise level); the rule of thumb `k = max(4, ceil(4 gamma))` is within 10% of the best observed delay in all 8 cells, and the guarantee holds down to `k = 4` (28 cells, 0.5% to 6.5% false alarms).

## Negative results and failed attempts

* Clipping (the paper's truncation) after pre-averaging has a jump bias that grows with `k`; it had to be replaced by dropping.
* The multipower (bipower, tripower) estimators, modulated or pre-averaged, are not jump robust at `n = 23,400` in the sense that matters for quarticity: a 100 sd jump leaves a factor 1.5 (two jumps) to 5 (eight jumps). The asymptotic robustness holds only as `n^{-1/4}`.
* The full text of Podolskij and Vetter (2009) could not be fetched (arxiv.org, projecteuclid.org, repec.econ.au.dk, eldorado.tu-dortmund.de unreachable from this environment), so the modulated bipower constants `nu1`, `nu2` and the CLT standard error are our own derivation, verified by simulation (`tests/test_preavg.py`), and no equation numbers are quoted.
* The first rule of thumb for `k` (cap the block noise share at 0.25, `k = 7 gamma`) was ad hoc and Q6 showed it conservative by a factor 2; it was replaced by the minimiser of the delay formula, which fits the data.
* The hypothesis that the raw-return monitor stays valid under price noise, which the main report's E4 supported at one noise level over one short day, is false over a full day at `gamma >= 0.5`.
* A first unit test compared a calibration path and monitored paths simulated with different numbers of steps, so their per-step scales differed by `sqrt(10)` and the detector "false alarmed" on every path. The detector was right; the test was fixed.
* A background watcher for the first run never exited because `pgrep -f` matched its own command line; the run itself had finished. No result was affected.

## Ranking

* **Research novelty:** moderate to high. The composition (pre-average, drop, noise-correct, then bet on blocks) is new as a monitor, and the finding that MA(1) price noise breaks the e-detector's conditional null on raw returns is a correction to the main package's claims.
* **Mathematical validity:** high for the estimator (unbiasedness verified to about 1% across 14 cells of Q1 and within 5% across the 20 noise x jump cells of Q5) and for the block monitor (the conditional null and the guarantee are stated and proved for every `k >= 2` in the README, with no boundary term; the bound holds in 30 cells of Q3 and 28 of Q6); the local-constancy derivation of the fourth-moment correction is heuristic and should be replaced by the exact pre-averaged quarticity of Jacod et al. (2009) in a paper, and the delay formula behind `k*` ignores truncation.
* **Practical relevance:** high at one-second sampling, where noise ratios of 0.2 to 2 are common; low at one minute for BTC (noise ratio well below the 0.16 breaking point). Against the published noise-and-jump estimator (Podolskij and Vetter 2009) the prototype is 2.5 to 3 times more precise and has no finite-sample jump bias.
* **Feasibility:** high. About 450 lines including the benchmark, no new dependencies, five minutes for Q1 to Q4 and two minutes for Q5 and Q6.
* **Demonstration value:** high. The three figures make the point in one glance (RMSE flat against a line that rises by a factor 1000; false alarms rising against flat; a U-shaped delay curve whose minimum a one-line formula predicts).
* **Likely interest to a sequential-analysis and high-frequency econometrics researcher:** high. The overlap and MA(1) failures are exactly the kind of assumption that e-detector users must check, and the pre-averaged block monitor is a clean, provable fix.
