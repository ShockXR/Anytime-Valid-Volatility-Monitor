# Pre-averaging before MAD truncation (alternative A1)

A minimal working prototype of the idea in `alternatives/ALTERNATIVES.md`, entry A1: average returns over short overlapping blocks (pre-averaging, Jacod, Li, Mykland, Podolskij and Vetter 2009) before clipping them at a multiple of the MAD and forming local variance and quarticity, so that market microstructure noise is removed as well as jumps. Results are in `RESULTS.md`; every number there comes from `run_experiment.py` and is labelled [SYNTHETIC].

## 1. Research problem

At high sampling frequencies (seconds), observed prices are the efficient price plus a small error from the bid-ask bounce, discrete ticks and stale quotes: "microstructure noise". The noise does not shrink with the sampling interval, while the true price movement does, so at one-second steps noise can dominate every squared return. Bapat (2026) claims noise robustness for the MAD-truncated kernel estimator, but its noise model (eq. 13) adds an independent error to each *return*, which is indistinguishable from extra volatility; `docs/estimator_issues.md`, point 4b, shows that no estimator can remove such noise and that the truncated estimator does not try. Under the standard model, where the error is added to each *price*, returns become negatively autocorrelated and pre-averaging can remove the noise.

Questions:

1. At what noise-to-signal ratio does the plain truncated estimator break, for integrated quarticity (bias and RMSE against the true value) and for the local variance used by the monitor, and does pre-averaging followed by truncation fix it?
2. What does the monitor's false alarm rate do under noise, with and without pre-averaging? Pre-averaged returns overlap: what does that do to the e-detector's assumptions, and is a non-overlapping block version a valid substitute?

## 2. Mathematical formulation

**Model.** Observed log prices `Y_j = X_j + u_j`, `j = 0, ..., n`, over a day of length `T` with step `delta = T / n`. The efficient price `X` has spot variance `sigma_t^2` (stochastic) and possibly jumps; the noise `u_j` is i.i.d. with variance `omega^2`, independent of `X`. Observed returns `r_j = Y_j - Y_{j-1}` have `E[r_j^2] = sigma^2 delta + 2 omega^2` and `E[r_j r_{j-1}] = -omega^2`. The noise-to-signal ratio per return is `gamma = omega / (sigma sqrt(delta))`. Clipping at `3 x MAD` of the raw returns keeps `(1 + 2 gamma^2)` of the true variance and `(1 + 2 gamma^2)^2` of the true quarticity: a 10% quarticity bias at `gamma = 0.16`, a factor 2 at `gamma = 0.45`, a factor 81 at `gamma = 2`.

**Pre-averaging.** With window `k` and weight `g(x) = min(x, 1 - x)`,

`rbar_i = sum_{j=1}^{k-1} g(j / k) r_{i+j}`, `i = 0, ..., n - k`.

Write `psi2_k = (1/k) sum_j g(j/k)^2` (about 1/12) and `psi1_k = k sum_j (g(j/k) - g((j-1)/k))^2` (exactly 1 for this `g`). Then, for locally constant `sigma`,

`E[rbar_i^2] = k psi2_k delta sigma^2 + psi1_k omega^2 / k`.

With `k = theta sqrt(n)` the two terms have the same order, which is why the method converges at rate `n^{1/4}` instead of `n^{1/2}` (Jacod et al. 2009). The local pre-averaged variance `V_i = rbar_i^2 / (k psi2_k delta)` has mean `sigma^2 + c` with the noise term `c = omega^2 psi1_k / (psi2_k k^2 delta)`, which in units of `sigma^2` is `gamma^2 psi1_k / (psi2_k k^2)`, about `gamma^2 / 300` for `k = 60`. The noise variance is estimated by `omega^2 = -mean(r_j r_{j-1})` (minus the first autocovariance), computed on MAD-clipped raw returns so that a jump cannot distort it.

**Local quarticity with the noise correction.** If `rbar_i` is conditionally Gaussian, `E[rbar_i^4] = 3 (E[rbar_i^2])^2`, so `rbar_i^4 / (3 (k psi2_k delta)^2)` has mean `(sigma^2 + c)^2`. The corrected local quarticity

`Q_i = rbar_i^4 / (3 (k psi2_k delta)^2) - 2 c V_i + c^2`

has mean `(s + c)^2 - 2c (s + c) + c^2 = s^2 = sigma^4`. Integrated quantities are `IV = T x mean_i(V_i) - c` and `IQ = T x mean_i(Q_i)`. This is a "plug in the local estimate" construction in the spirit of Jacod and Rosenbaum (2013), simplified by the local-constancy assumption; the full pre-averaged quarticity of Jacod et al. (2009) has additional small-sample terms that are not needed here.

**Truncation after pre-averaging.** `rbar_i` is clipped at `tau = c_mad x scaled MAD(rbar)` with `c_mad = 3` (eq. 11 of Bapat 2026 applied to the pre-averaged returns). The clipped second and fourth moments are divided by the Gaussian survival factors `kappa_2(3) = 0.995` and `kappa_4(3) = 0.964` (`docs/math.md`, section 3), after which the noise correction above goes through unchanged. Two modes are compared: `clip` (replace `|rbar| > tau` by `tau`, as in the paper) and `drop` (replace by 0, a Mancini-style threshold), with the matching factors. The difference matters for jumps: a jump enters `k - 1` overlapping windows, and clipping leaves `tau^4` in each of them, which is about `(k - 1) x 27 / n` of the day's quarticity per jump (7% for `k = 60`, `n = 23,400`).

**Local variance for the monitor.** Over a trailing window of `h` seconds the monitor's local variance is the mean of `min(r^2, tau^2) / (kappa_2 delta)` on raw returns (mean `sigma^2 (1 + 2 gamma^2)`), or the mean of the noise-corrected `V_i` over the pre-averaged windows that end inside the last `h` seconds (mean `sigma^2`).

**Monitoring and the overlap problem.** The e-detector of the main package needs `E[Y_t | past] <= 1` under the null, where the past includes the previous monitored values. For overlapping pre-averaged returns, `rbar_t` shares `k - 2` raw returns with `rbar_{t-1}`, so after a large `rbar_{t-1}` the next value is expected to be large too: the conditional mean is not constant even with constant volatility and i.i.d. noise. The supermartingale argument fails, and the false alarm bound with it. (The marginal mean is unaffected, which is why the *estimators* above are fine.) Taking every `k`-th window (`rbar_0, rbar_k, rbar_2k, ...`) gives blocks that use disjoint sets of prices (`Y_i, ..., Y_{i+k-1}` and `Y_{i+k}, ..., Y_{i+2k-1}`), so under the model they are independent, the conditional-mean null holds, and the guarantee holds with `n / k` observations per day. The block monitor therefore runs at alarm level `A = (n / k) / 0.05` for a 5% per-day budget, against `n / 0.05` for the raw monitor.

**Formal statement for the block monitor.** Prices `Y_j = X_j + u_j`, `j = 0, 1, ...`, with `u_j` i.i.d. with mean zero and variance `omega^2` (any distribution), independent of `X`; `X` an Ito semimartingale with spot variance `sigma_t^2` and a finite-activity jump part. Block `b` (`b = 0, 1, ...`) uses the prices `Y_{bk}, ..., Y_{bk+k-1}` and the statistic `Z_b = min(rbar_b^2, tau^2)` with `rbar_b = sum_{j=1}^{k-1} g(j/k) (Y_{bk+j} - Y_{bk+j-1})`. Let `G_b` be all information up to the last price time of block `b` (the paths of `X`, `sigma`, the jumps and the noise up to `(bk + k - 1) delta`). Then:

1. *Independence of the past.* `rbar_{b+1}` is a function of the `X` increments over `((bk+k-1) delta, (bk+2k-1) delta]` and of the noise values `u_{(b+1)k}, ..., u_{(b+1)k+k-1}`. The return `Y_{(b+1)k} - Y_{(b+1)k-1}` between the blocks is skipped, so no price of block `b` is reused, and the noise values are fresh and independent of `G_b`. Given `G_b`, `rbar_{b+1}` therefore has the law of a pre-averaged return of a process started at time `(bk+k-1) delta` with its own spot-volatility and jump dynamics and independent noise. No boundary term remains: the independence of the noise parts is exact, not up to `O(1/k)`. (If adjacent blocks shared the price `Y_{bk+k-1}`, that noise value would enter `rbar_b` with weight `g((k-1)/k) = 1/k` and `rbar_{b+1}` with weight `-g(1/k) = -1/k`, giving `Cov(rbar_b, rbar_{b+1}) = -omega^2 / k^2`, an `O(1/k)` relative dependence against the noise part `omega^2 / k` of the variance. Dropping one return between blocks removes it.)
2. *Conditional null.* Write `m0 = E[Z]` under the calibration regime (spot variance `sigma_0^2`, the same noise and jump laws). The null is `H0: E[Z_{b+1} | G_b] <= m0` for every `b`. By 1 and the monotonicity of `min(s^2 e^2, tau^2)` in `s` for every innovation `e`, it holds whenever the spot variance over block `b+1` is at most `sigma_0^2` and the noise and jump laws are unchanged; the value of `m0` includes the noise term `psi1 omega^2 / k`.
3. *Guarantee.* Under `H0` the betting factor `f_{b+1} = 1 + lambda_{b+1} (Z_{b+1} / m0 - 1)` with `G_b`-measurable `lambda_{b+1} in [0, 1)` satisfies `E[f_{b+1} | G_b] <= 1`, so the Shiryaev-Roberts statistic over blocks is an e-detector with respect to `(G_b)`, and Theorems A and B of `docs/math.md` give `E[number of blocks to a false alarm] >= A` and `P(false alarm within m blocks) <= m / A`. With `m = floor(n / k)` blocks per day and `A = m / 0.05`, the per-day false alarm probability is at most 5%. Nothing requires Gaussian blocks, so the guarantee holds for every `k >= 2`.

**Why the raw-return monitor fails under price noise.** With `r_t = Delta X_t + u_t - u_{t-1}` and `F_{t-1}` containing `Y_{t-1} = X_{t-1} + u_{t-1}`, `E[r_t^2 | F_{t-1}] = sigma^2 delta + omega^2 + E[u_{t-1}^2 | F_{t-1}]`. The last term has mean `omega^2` but is not constant: for Gaussian noise and a constant `sigma`, conditioning on the previous return alone gives `E[r_t^2 | r_{t-1}] = sigma^2 delta + 2 omega^2 + rho_u^2 (r_{t-1}^2 - sigma^2 delta - 2 omega^2)` with `rho_u = gamma^2 / (1 + 2 gamma^2)`. After a larger than average squared return the conditional mean exceeds `m0`, so `H0` fails with no change in volatility, and the false alarm bound with it. The size of the violation grows with `gamma` (`rho_u^2` = 0.11, 0.20, 0.24 at `gamma` = 1, 2, 4).

**Delay cost of the block size.** After a change of the sd by a factor `v`, the block statistic's mean moves to about `(v^2 + rho_k) / (1 + rho_k)` times its baseline, where `rho_k = psi1 gamma^2 / (psi2 k^2)` (about `12 gamma^2 / k^2`) is the noise share of a block (`preavg.block_noise_share`). Small `k` gives more blocks per day but a diluted change; large `k` gives an undiluted change but few blocks and a delay that grows roughly in proportion to `k`. The trade-off can be solved. For Gaussian blocks (truncation ignored) the monitored value `Y = Z / m0` has, after the change, mean `(v^2 + rho_k) / (1 + rho_k)` and standard deviation `sqrt(2)` times that, so the betting e-detector's growth per block is about `snr^2 / 2` with `snr = (v^2 - 1) / (sqrt(2) (v^2 + rho_k))` (second-order expansion of `E[log f]` at the best fixed bet `lambda = (E[Y] - 1) / Var(Y)`), and the delay in seconds is about

`D(k) = k x 2 log(A) / snr^2 = 4 k log(A) (v^2 + rho_k)^2 / (v^2 - 1)^2`, `rho_k = 12 gamma^2 / k^2`.

`dD/dk = 0` at `rho_k = v^2 / 3`, that is at `k* = 6 gamma / v`: the best block carries a noise share equal to one third of the post-change variance ratio. The curve is flat near `k*` (within 10% for `k` between `k* / 1.5` and `1.5 k*`), so one target `v = 1.5` gives the rule of thumb `k = max(4, ceil(4 gamma))` (`preavg.rule_of_thumb_k`), with `gamma` estimated from `omega^2 = -mean(r_t r_{t-1})` divided by the diffusive return variance `mean(r_t^2) - 2 omega^2`. (A first, ad hoc rule that capped the dilution `rho_k` at 0.25, `k = 7 gamma`, turned out in Q6 to be conservative by a factor 2; the derivation above replaced it and predicts the observed minima within one grid step.) Q6 checks the formula against simulated delays.

## 3. Why it matters

* It fixes the one claim of Bapat (2026) that the main report could not support: noise robustness. The truncation handles jumps; pre-averaging handles noise; the two compose, with a small and explicit correction.
* It makes the monitor usable at one-second sampling, where the information per day is 60 times larger than at one minute but noise would otherwise inflate every statistic and dilute every change. The dilution is one practical point: under noise the raw monitor sees a variance doubling as a change of `(4 + 2 gamma^2) / (1 + 2 gamma^2)`, which is 1.33 at `gamma = 2`. The other is validity: Q3 shows that the raw monitor also loses its false alarm guarantee under price noise, because MA(1) returns break the conditional null (section 2), so at one-second sampling the block monitor is not only faster from `gamma` of about 2 but the only valid one from `gamma` of about 0.5.
* The overlap discussion is a general warning for e-detectors on smoothed statistics: smoothing before betting silently breaks the conditional null.

## 4. Comparison with existing approaches

| Approach | Noise | Jumps | Online | Guarantee for the monitor | Relation |
|---|---|---|---|---|---|
| MAD-truncated kernel (Bapat 2026) on raw returns | no (bias `(1 + 2 gamma^2)^2`) | yes (clipping) | yes | yes, but diluted | Baseline. |
| Pre-averaging (Jacod et al. 2009), no truncation | yes | no (a jump of 50 sd dominates the fourth power) | yes | fails with overlap | First stage here. |
| Pre-averaging, then MAD clipping or dropping (this prototype) | yes | yes | yes | yes with blocks | Composition, with a noise correction for the fourth moment. |
| Modulated bipower variation (Podolskij and Vetter 2009) | yes | yes | offline | none | The established estimator for noise plus jumps; a natural benchmark for a later version. |
| Two-scale realized variance (Zhang, Mykland and Ait-Sahalia 2005) | yes | no | offline | none | Alternative noise treatment for variance only. |
| Sampling sparsely (one minute instead of one second) | partly (noise ratio falls as `sqrt(60)`) | no | yes | yes | What the main package does; it throws away 59 of every 60 observations. |

## 5. Implementation

`preavg.py` (imports `avmonitor`, nothing in `src` imports it):

* `preavg_weights`, `psi_constants`, `choose_k`: the weight function, the discrete constants `psi1_k`, `psi2_k`, and `k = theta sqrt(n)`.
* `preaverage(r, k, step)`: pre-averaged returns, overlapping (`step = 1`) or in blocks (`step = k`).
* `noise_variance`: minus the first autocovariance of MAD-clipped raw returns.
* `truncation_factors(a, mode)`, `truncate_preavg`: Gaussian survival factors for clipping and dropping, and the truncation itself.
* `local_moments`: noise-corrected local variance `V_i` and local quarticity `Q_i`.
* `preavg_estimates`: integrated variance and quarticity from pre-averaged, truncated, noise-corrected returns. `raw_truncated_estimates`: the plain route for comparison (normalised, kappa-corrected truncated RV and RQ, which equal the averaged kernel estimator of eq. 16 up to edge effects).
* `spot_variance_raw`, `spot_variance_preavg`: trailing-window local variance on raw and on pre-averaged returns.
* `monitor_stream`: the bounded stream `min(rbar^2, tau^2)` and its bound for the e-detector, overlapping or in blocks.

`run_experiment.py` runs Q1 to Q4 and draws the two figures; `python run_experiment.py q3 q4` reruns parts.

## 6. Experimental methodology [SYNTHETIC]

One-second returns, `n = 23,400` per day, `T = 1` day, daily volatility 1% (per-second return sd `6.5 x 10^-5`), log-OU spot variance (mean reversion 5, vol of vol 0.5) for Q1 and Q2 and constant volatility for Q3 and Q4, price-level noise with `omega = gamma x per-second sd`, jumps as a compound Poisson process with 4 jumps per day of size 50 per-second sd (0.33% of the price). `k = 60` (`theta = 0.39`), `c_mad = 3`, fixed seeds, two worker processes.

* **Q1 Integrated quarticity and variance.** `gamma` in `{0, 0.1, 0.25, 0.5, 1, 2, 4}`, 0 and 4 jumps per day, 200 paths per cell. Estimators: plain RQ, raw returns with MAD clip (Bapat, normalised and kappa-corrected), pre-averaged without truncation, pre-averaged with clip, pre-averaged with drop, pre-averaged with clip but no noise correction. Mean ratio to the truth and relative RMSE. Window sensitivity: `k` in `{20, 40, 60, 100, 150}` at `gamma = 1` and `4`. Outputs `results/q1_iq_iv_vs_noise.csv`, `results/q1_window_sensitivity.csv`.
* **Q2 Spot variance.** Same paths (100 per cell, `gamma` in `{0, 0.5, 1, 2, 4}`), local variance over the last 300 seconds at every 600th second, compared with the true spot variance: raw truncated, pre-averaged overlapping, pre-averaged blocks. Mean ratio, median and 95th percentile of the absolute log error. Output `results/q2_spot_variance.csv`.
* **Q3 Monitor false alarms.** Constant volatility, `gamma` in `{0, 0.5, 1, 2, 4}`, 200 paths; two days of calibration (46,800 returns) then one monitored day. Monitors, each at `A = steps per day / 0.05`: raw truncated returns (23,400 steps, `A = 468,000`), overlapping pre-averaged returns with `k = 60` (included to measure the violation), and non-overlapping blocks with `k = 60`, `20` and `10` (390, 1170 and 2340 blocks per day, `A = 7800`, `23,400` and `46,800`); increase and decrease sides. One-day false alarm rate against the 5% bound. Output `results/q3_monitor_false_alarms.csv`.
* **Q4 Detection under noise.** As Q3 but the sd changes by `x 1.5`, `x 2` or `x 0.5` at the start of the monitored day. Raw monitor against the three block monitors: share detected within the day and median delay in seconds (counted from the change, including the `k - 1` seconds a block needs to fill). Output `results/q4_monitor_detection.csv`.

### Benchmark and block-size study (`run_benchmark.py`)

* **Q5 Noise and jumps together [SYNTHETIC].** Grid: `gamma` in `{0, 0.5, 1, 2}` x (jumps per day, jump size in per-second sd) in `{(0, 0), (2, 25), (2, 100), (8, 25), (8, 100)}`, log-OU stochastic volatility, 200 paths per cell. Estimators of integrated quarticity and variance: raw MAD clip; pre-averaged MAD drop and MAD clip (`k = 60`); the modulated bipower (variance) and tripower (quarticity) variation of Podolskij and Vetter (2009), block averages of `K = 0.5 sqrt(n)`-step increments over blocks of `sqrt(n)` observations (`K = 76`, block 153, 152 blocks per day), with the noise correction through the exact block weights (`preavg.mbv_block_weights`, asymptotically `nu1 = c1^2 (3 c2 - 4 c1) / (3 (c2 - c1)^2)` and `nu2 = 2 c1 / (c2 - c1)^2`); the pre-averaged bipower and tripower variation (`k = 60`, factors `k` apart); and pre-averaging without truncation. Bias (mean ratio), relative RMSE, and for the modulated bipower variance estimator the coverage and relative half-width of the 95% interval whose standard error comes from the independent Gaussian block approximation (`Var(sum |Ybar_m Ybar_{m+1}|) = (M - 1) (mu_2^2 + 2 mu_2 mu_1^2 - 3 mu_1^4) mean(v_m^2)`, our derivation; the paper's theorem could not be fetched to confirm its constants). Output `results/q5_benchmark_noise_jumps.csv`.
* **Q6 Block size [SYNTHETIC].** Constant volatility, `gamma` in `{0.5, 1, 2, 4}`, blocks `k` in `{4, 6, 10, 14, 20, 30, 60}`, two calibration days and one monitored day, `A = (n / k) / 0.05`; sd `x 2`, `x 1.5` and no change (one-day false alarm rate). Median delay in seconds against `k`, the noise share `rho_k`, the formula optimum `k* = 6 gamma / v` and the rule-of-thumb `k` (target `v = 1.5`). Output `results/q6_block_size.csv`, figure `results/fig3_block_size_rule.png`.

**Reference check.** Podolskij, M. and Vetter, M. (2009), Estimation of volatility functionals in the simultaneous presence of microstructure noise and jumps, *Bernoulli* 15(3), 634 to 658, is listed on arXiv as 0909.0827 (title, authors, the modulated bipower construction, jump robustness for powers below 2 and the `n^{-1/4}` rate confirmed from the abstract and search listings in November 2026). The full text could not be fetched in this environment (arxiv.org, projecteuclid.org, repec.econ.au.dk and eldorado.tu-dortmund.de are not reachable), so equation numbers are not quoted and the constants were derived from the block weights and verified by simulation (`tests/test_preavg.py`).

## 9. Folder structure

```
alternatives/preaveraging/
├── README.md             this file (items 1 to 6 and 9)
├── RESULTS.md            results with every number labelled [SYNTHETIC] (item 8)
├── preavg.py             prototype module, benchmarks and the block-size rule (items 5, 7)
├── run_experiment.py     Q1 to Q4 and figures 1 and 2 (items 6, 7)
├── run_benchmark.py      Q5 (Podolskij and Vetter benchmark) and Q6 (block size), figure 3
├── tests/
│   └── test_preavg.py    8 pytest tests, about 3 seconds
└── results/
    ├── q1_iq_iv_vs_noise.csv
    ├── q1_window_sensitivity.csv
    ├── q2_spot_variance.csv
    ├── q3_monitor_false_alarms.csv
    ├── q4_monitor_detection.csv
    ├── q5_benchmark_noise_jumps.csv
    ├── q6_block_size.csv
    ├── fig1_iq_rmse_vs_noise.png
    ├── fig2_monitor_under_noise.png
    └── fig3_block_size_rule.png
```

Run from the repository root with the package installed: `python alternatives/preaveraging/run_experiment.py` (about 5 minutes on two cores, no data needed), `python alternatives/preaveraging/run_benchmark.py` (about 2 minutes) and `python -m pytest -q alternatives/preaveraging/tests`.

## References

* Jacod, J., Li, Y., Mykland, P. A., Podolskij, M. and Vetter, M. (2009). Microstructure noise in the continuous case: the pre-averaging approach. *Stochastic Processes and their Applications*, 119(7), 2249 to 2276. doi:10.1016/j.spa.2008.11.004 (reference checked online, November 2026).
* Jacod, J. and Rosenbaum, M. (2013). Quarticity and other functionals of volatility: efficient estimation. *Annals of Statistics*, 41, 1462 to 1484.
* Podolskij, M. and Vetter, M. (2009). Estimation of volatility functionals in the simultaneous presence of microstructure noise and jumps. *Bernoulli*, 15(3), 634 to 658.
* Zhang, L., Mykland, P. A. and Ait-Sahalia, Y. (2005). A tale of two time scales. *JASA*, 100(472), 1394 to 1411.
* Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482.
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science*, 2(2), 229 to 260.
