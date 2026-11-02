# Report part ALT3: a robust local plug-in estimator of integrated quarticity (A5)

Code, tests, README (mathematics, proofs, comparison with the literature) and full results are in `alternatives/robust_iq_estimator/`. Every number below is copied from its `results/` folder. [SYNTHETIC] marks simulated data or closed-form theory, [REAL] marks Binance spot one-minute bars (BTCUSDT and DOGEUSDT, 1 July to 31 December 2025, 184 UTC days each). Iteration notes: `docs/iterations/ALT3.md`.

**Question.** Integrated quarticity (`IQ = int sigma^4`, the quantity that sets the width of every confidence interval for daily variance) is estimated in Bapat (2026) by clipping returns at `c x MAD` and averaging kernel-weighted local fourth powers. `docs/estimator_issues.md` lists four problems with it. Can the paper's idea (robust-scale truncation plus local smoothing) be turned into an estimator that is consistent, close to the efficient asymptotic variance `8 int sigma^8`, robust to jumps, equipped with a confidence interval, and better behaved on real minute data?

**Method.** Two changes of principle and one of detail. (1) Aggregate by **squaring a local variance** (the plug-in principle of Jacod and Rosenbaum 2013) instead of averaging local fourth powers: over blocks of k returns, `IQ_hat = sum over blocks of k delta (block mean of r^2 / delta)^2`, whose asymptotic variance is `8 (k + 3) / (k + 2) int sigma^8` (RQ-type estimators, including the paper's averaged kernel, are at `32/3`). (2) Truncate against a **local** scale (median absolute deviation, then a truncated root mean square, over about three blocks) with either a fixed cut-off a = 3 or the growing cut-off `a_n = sqrt(2 log n)`, and **drop** truncated returns rather than clip them. (3) Remove the squaring and truncation biases with an **exact finite-sample correction**, computed from the cumulants of a clipped or dropped squared normal, so the estimator is unbiased under locally constant volatility for every block size. The same algebra gives `int sigma^8` and hence a 95% interval. Main proposal **PI-LU** (local scale, `a_n`, drop, `k = ceil(sqrt n)`). Second proposal **PI-NN** (added in iteration 3, after the real data): PI-LU with the scale raised to 1.4826 times the median of the six nearest absolute returns when that is larger, and blocks of `round(n^(1/4))`. Simulations: log-OU stochastic volatility from `avmonitor.simulate`, vol of vol 0.5 or 1 with a U-shaped intraday pattern, 0 or 4 jumps a day of three sizes, n = 390, 1440, 23,400, 500 days per cell, 15 estimators, true IQ known exactly.

**Q1 Efficiency [SYNTHETIC].** Constant volatility, n = 1440, 8000 days. AVAR in units of `int sigma^8` (closed form / simulated) and efficiency `8 / AVAR`:

| RQ | TPQ | MedRQ | MinRQ | RK-k (Bapat, normalised, corrected) | PI-LU (k = 38) | PI-NN (k = 6) |
|---|---|---|---|---|---|---|
| 10.67 / 10.81, 0.74 | 13.65 / 13.48, 0.59 | 14.24 (sim.), 0.56 | 18.50 (sim.), 0.43 | 9.26 (sim.), 0.86 | 8.20 / 8.46, 0.95 | 9.00 / 9.17, 0.87 |

A jump clipped at a = 3 leaves 17.55 (plug-in) or 27.01 (fourth powers) times `sigma^4 delta` (simulation 17.56 and 27.02); a dropped jump leaves 0.003. At n = 390 each clipped jump is 4.5% to 6.9% of the day's IQ.

**Q2 Consistency [SYNTHETIC].** Strong stochastic volatility plus intraday pattern, no jumps (S5): relative bias at n = 390 / 1440 / 23,400 and slope of log RMSE on log n (minus 0.5 is the `sqrt(n)` rate):

| | RK (Bapat) | RK-k | PI-G3 (daily scale) | TRQ-M (Mancini) | MedRQ | PI-LU | PI-NN |
|---|---|---|---|---|---|---|---|
| bias | -0.209 / -0.159 / -0.144 | -0.179 / -0.128 / -0.112 | -0.069 / -0.061 / -0.055 | -0.038 / -0.028 / -0.024 | -0.017 / 0.005 / 0.002 | -0.002 / -0.001 / 0.000 | 0.008 / 0.006 / 0.001 |
| slope | -0.11 | -0.14 | -0.22 | -0.36 | -0.50 | -0.51 | -0.51 |

**Q3 Bandwidth [SYNTHETIC].** The paper's cross-validation (eq. 17 as printed) picks the smallest bandwidth offered (h = 5) on 58% to 97% of days. Its expected value under constant volatility is 74.3, 77.2, 99.2, 142.8 at h = 5, 10, 20, 30 (minimum at h = 5 whatever the data). A corrected target picks the oracle spot bandwidth (median h = 50 to 150), but that oversmooths for the integrated quantity (kernel plug-in RMSE 0.250 at the corrected choice against 0.198 at h = 10, S5, n = 390). The block length of PI-LU matters little in smooth designs (default within 1% to 8% of the best k).

**Q4 Comparison and intervals [SYNTHETIC].** Relative RMSE at n = 390 (n = 1440 in brackets):

| | S1 no jumps | S3 medium jumps | S4 large jumps | S5 strong SV | S6 S5 + jumps |
|---|---|---|---|---|---|
| RQ | 0.167 (0.091) | 28.2 (102) | 3017 (12917) | 0.197 (0.111) | 31.2 (95.1) |
| TPQ | 0.189 (0.100) | 0.630 (0.379) | 3.54 (2.99) | 0.239 (0.130) | 0.782 (0.383) |
| MedRQ | 0.186 (0.106) | 0.471 (2.10) | 177 (93.4) | 0.236 (0.127) | 19.3 (1.19) |
| RK (Bapat) | 0.156 (0.091) | 0.249 (0.091) | 0.312 (0.089) | 0.256 (0.184) | 0.198 (0.154) |
| PI-LU | 0.154 (0.078) | **0.179 (0.082)** | **0.163 (0.078)** | 0.183 (0.098) | 0.217 (0.106) |
| PI-NN | 0.162 (0.085) | 0.320 (0.121) | 0.236 (0.088) | 0.192 (0.105) | 0.375 (0.129) |

PI-LU's 95% interval covers the true IQ 92.6% to 97.2% of the time in all 18 cells (simulation standard error about 1%).

**Q5 Real data [REAL].** Median over days of each estimator's ratio to MedRQ of the same day, and rank correlation with MedRQ:

| | RK (Bapat) | RK-k | PI-G3 | TRQ-M | PI-L3 | PI-LU | KPI-L3 | PI-NN | MinRQ | TPQ |
|---|---|---|---|---|---|---|---|---|---|---|
| BTC | 0.33 | 0.35 | 0.38 | 0.63 | 0.70 | 0.66 | 0.72 | **0.95** | 1.03 | 0.99 |
| DOGE | 0.42 | 0.44 | 0.46 | 0.66 | 0.71 | 0.69 | 0.73 | **0.96** | 1.05 | 1.03 |
| Spearman, BTC / DOGE | 0.865 / 0.814 | 0.865 / 0.814 | 0.879 / 0.837 | 0.927 / 0.890 | 0.939 / 0.908 | 0.931 / 0.888 | 0.941 / 0.911 | **0.975 / 0.965** | 0.995 / 0.994 | 0.987 / 0.987 |

**Q6 Bursts or jumps [REAL + SYNTHETIC].** Locally standardised one-minute returns: 65% (BTC) and 44% (DOGE) of returns beyond `sqrt(2 log n)` have a neighbour beyond 2 sd, and the lag-1 autocorrelation of capped squared returns is 0.114 and 0.099; isolated jumps on a slowly varying volatility give 0% and minus 0.008. A two-factor volatility model with a fast factor (sd 0.75 in log variance, half-life 3 minutes), picked from a 3 by 3 grid on these diagnostics, gives 41% and 0.120 and reproduces the real ratios (PI-LU / MedRQ 0.68 against 0.66 to 0.69; PI-NN 0.95 against 0.95 to 0.96; 1-minute over 5-minute MedRQ 1.41 against 1.33 to 1.45). In that model, median relative error and RMSE at n = 390 / 1440 / 23,400:

| | RK (Bapat) | PI-LU | MedRQ | TPQ | PI-NN |
|---|---|---|---|---|---|
| median error, no jumps | -0.58 / -0.59 / -0.57 | -0.49 / -0.47 / -0.26 | -0.40 / -0.20 / -0.02 | -0.40 / -0.20 / -0.02 | -0.36 / -0.23 / -0.04 |
| RMSE, no jumps | 0.59 / 0.59 / 0.59 | 0.49 / 0.46 / 0.27 | 0.42 / 0.29 / 0.08 | 0.42 / 0.28 / 0.08 | **0.39 / 0.28 / 0.07** |
| RMSE, 4 jumps a day | 0.54 / 0.56 / 0.58 | 0.48 / 0.44 / 0.27 | 0.89 / 1.63 / 1.28 | 0.43 / 0.51 / 0.17 | **0.36 / 0.26 / 0.07** |

**Reading.**
* Squaring a local variance instead of averaging fourth powers raises efficiency from 0.74 (RQ, and Bapat's averaged kernel, which is RQ in disguise) to 0.95, and the exact correction makes it unbiased at every n (without it PI-LU is 9.2% too high at n = 390).
* Bapat's rule with one daily MAD is inconsistent when volatility moves during the day: its bias does not shrink (minus 14% at n = 23,400) and its RMSE falls like `n^-0.1`. A local scale restores the `sqrt(n)` rate for every truncation rule; the growing cut-off with dropping (PI-LU) also removes jumps without a residual and is the best estimator in the jump scenarios at n = 390 and 1440.
* The paper's cross-validation choice is a scale artefact, and even a corrected CV answers the wrong question for an integrated quantity: the IQ estimator wants less smoothing than the spot estimator.
* On real data the local estimators agree with MedRQ and MinRQ far better than the MAD-truncated ones, but only PI-NN closes the gap. The diagnosis: large one-minute crypto returns come in bursts of a few minutes, which are part of volatility. In a model calibrated to that, every estimator that truncates against a scale measured over hours (Bapat's, PI-LU, PI-L3, the kernel plug-in) removes 20% to 59% of IQ, while the nearest-neighbour fourth-power estimators keep bursts but break with jumps. PI-NN combines both properties and has the lowest RMSE among robust estimators there, with and without jumps.
* Practical rule: compute the clustering diagnostic first. Isolated large returns: use PI-LU. Clustered ones: use PI-NN, and do not trust agreement with a 5-minute estimate (in the model the 5-minute MedRQ is 0.57 of the true IQ at the median).

**Caveats.** PI-NN was designed after seeing the real data, with settings chosen in exploratory runs (other seeds), so its numbers are not strictly out of sample; it is worse than PI-LU with isolated jumps at small n (S3, n = 390: 0.320 against 0.179) and its `int sigma^8` estimate is not robust. The burst reading is a model fit, not a proof: clustered small jumps or microstructure effects could leave similar signatures, and the daily-scale estimators fit the real data less well (RK 0.53 in the model, 0.33 to 0.42 real). No central limit theorem is proved for the feasible estimators; coverage is numerical evidence. One venue, two coins, six months, no microstructure noise in the simulations.
