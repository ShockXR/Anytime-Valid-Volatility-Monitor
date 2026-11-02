# Results: robust local plug-in estimator of integrated quarticity (A5)

Every number below comes from `run_experiment.py` (seeds fixed in the code) and is stored in `results/`. Labels: **[SYNTHETIC]** for simulations and closed-form theory, **[REAL]** for Binance one-minute data. Estimator labels are defined in `README.md`, section 5. RK and RK-k are this repository's version of the paper's kernel estimator (normalised weights, scaled MAD, c = 3, h = 30), so every comparison is with that version, not with the formula as printed. PI-LU and PI-NN are proposals for discussion: they are more efficient or less biased in our simulations, which is not the same as having improved the published estimator. "Relative bias" and "relative RMSE" are `estimate / true IQ - 1` averaged over simulated days (RMSE: root mean square). The true IQ is the exact integral of the simulated volatility path.

## Summary

1. **Efficiency [SYNTHETIC].** Under constant volatility, squaring a local variance gives an asymptotic variance of `8 (k + 3) / (k + 2)` times `int sigma^8`: 8.20 for k = 38 (8.46 simulated, efficiency 0.95), close to the efficient bound 8. The gain is large against untruncated RQ (10.67, efficiency 0.74), TPQ (13.6), MedRQ (14.2) and MinRQ (18.5). It is modest against the paper's kernel with the kappa correction: that estimator is a truncated RQ with edge weights, and truncation already lowers its variance, to 8.19 with an oracle cut-off (closed form, efficiency 0.98) and 9.26 with the estimated daily MAD (simulated, efficiency 0.86).
2. **Consistency [SYNTHETIC].** With a fixed c and one daily MAD, the kernel estimator converges to `kappa_4(c)` times IQ under constant volatility (minus 3.6% with the scaled MAD, minus 29.5% with the raw MAD), a known constant that the kappa correction removes. When volatility moves within the day, the corrected estimator keeps a bias that does not shrink with n: minus 0.6% to minus 1.1% in our moderate design S1, and minus 11% at n = 23,400 in the strong design S5 (minus 14% without the correction), where its RMSE falls like `n^-0.14` instead of `n^-0.5`. This setting lies outside the conditions of the paper's Theorem 1, which assumes a threshold that shrinks more slowly than the returns; a multiple of the MAD shrinks at the same rate. All estimators with a **local** scale converge at the `sqrt(n)` rate (slopes minus 0.48 to minus 0.68).
3. **Jumps [SYNTHETIC].** Dropping with the growing cut-off `sqrt(2 log n)` (PI-LU) leaves no residual from detected jumps; clipping (the paper's choice) leaves `O(1/n)` per jump, which is 18% to 28% of IQ at n = 390 with four jumps a day (closed form, a = 3). PI-LU has the lowest RMSE of the 15 estimators in 6 of the 8 jump cells at n = 390 and 1440. In S2 and S6 at n = 390 the paper's kernel and PI-G3 are slightly better (S2: RK 0.156, PI-G3 0.158, PI-LU 0.165; S6: PI-G3 0.191, RK-k 0.197, RK 0.198, PI-LU 0.217), because their downward truncation bias offsets the upward jump bias. PI-LU's 95% intervals cover 93% to 97% in every one of 18 cells.
4. **Cross-validation [SYNTHETIC].** Eq. (17) as printed picks the smallest bandwidth on offer in 58% to 97% of simulated days, for scale reasons alone (its expected value is minimised at h = 5 under constant volatility). The corrected target tracks the oracle spot-optimal bandwidth, but that bandwidth is too long for the *integrated* quantity, and no CV variant can see volatility bursts of a few minutes.
5. **Real data [REAL].** On Binance BTCUSDT and DOGEUSDT (July to December 2025) the estimators with one daily MAD threshold (RK, RK-k, PI-G3) are 0.33 to 0.46 times MedRQ at one minute (RK alone 0.33 to 0.42, that is 2.4 to 3.0 times lower); the local plug-ins (PI-LU, PI-L3, KPI-L3) are 0.66 to 0.73, and the neighbour-guarded PI-NN 0.95 to 0.96 (Spearman 0.97 with MedRQ). The true IQ is unknown on real data, so these are ratios between estimators, not errors.
6. **A possible explanation of the gap [REAL + SYNTHETIC].** Large one-minute crypto returns are mostly not isolated: 44% to 65% of them have a neighbour beyond 2 local standard deviations, against 8.9% expected by chance for independent Gaussian returns, and squared standardised returns stay autocorrelated for about 5 to 10 minutes. That many apparent jumps are short bursts of volatility is already documented for equities and currencies (Christensen, Oomen and Podolskij 2014; Bajgrowicz, Scaillet and Treccani 2016). What this prototype adds is the consequence for truncated IQ estimators. A two-factor volatility model with a fast factor (half-life 3 minutes, sd 0.75 in log variance) calibrated to these diagnostics reproduces the real estimator ratios (for example PI-LU / MedRQ 0.68 simulated against 0.66 to 0.69 real; PI-NN 0.95 against 0.95 to 0.96; MedRQ 1-minute / 5-minute 1.41 against 1.33 to 1.45). If that model is right, **every estimator that truncates against a scale estimated over hours is biased down by about 35% to 59% at one minute (31% to 56% with jumps)**; for our version of the paper's kernel the figure is 56% to 59% at every n. PI-NN has the lowest RMSE among the robust estimators in five of the six model cells (0.26 to 0.28 at n = 1440 against 0.28 to 1.63 for TPQ, MedRQ and MinRQ; at n = 1440 without jumps it leads TPQ only by 0.277 to 0.280, and at n = 390 with jumps TRQ-M is slightly lower, 0.352 against 0.360).
7. **No single winner.** With isolated jumps and smooth volatility, PI-LU does best; with volatility bursts, PI-NN does best in our calibrated model; PI-NN pays for its burst robustness with a larger error from undetected jumps at n = 390 (RMSE 0.32 against 0.18 for PI-LU in S3).

## Q1. Efficiency and constants [SYNTHETIC]

Constant volatility, n = 1440, 8000 simulated days (`results/q1_avar.csv`). AVAR is `n Var(estimate)` in units of `int sigma^8` (T = 1). Simulation standard errors are about 1.6% of each value.

| Estimator | AVAR closed form | AVAR simulated | Efficiency `8 / AVAR` (simulated AVAR where there is one) |
|---|---|---|---|
| Efficient bound | 8 | | 1 |
| RQ | 10.67 | 10.81 | 0.74 |
| TPQ (p = 4/3) | 13.65 | 13.48 | 0.59 |
| QPQ | 14.92 | | 0.54 |
| MedRQ | | 14.24 | 0.56 |
| MinRQ | | 18.50 | 0.43 |
| TRQ-M (Mancini threshold) | 10.67 | 10.61 | 0.75 |
| RK (paper's kernel, our normalised version, no kappa correction) | | 8.60 | biased, mean 0.964 of IQ |
| RK-k (same, kappa-corrected) | 8.19 (oracle cut-off: efficiency 0.98) | 9.26 | 0.86 |
| TRQ-U (fourth powers, local threshold) | 10.67 | 10.52 | 0.76 |
| PI-LU (k = 38) | 8.20 (efficiency 0.976) | 8.46 | 0.95 |
| PI-L3 (k = 38) | 8.20 | 8.43 | 0.95 |
| KPI-L3 (kernel plug-in, h = 30) | | 8.53 | 0.94 |
| PI-NN (k = 6) | 9.00 | 9.17 | 0.87 |
| PI-LU without the finite-sample correction | | 9.26 | 0.86 (and biased, mean 1.048) |
| Block plug-in, no truncation, k = 2, 5, 10, 20, 40, 160 | 10.00, 9.14, 8.67, 8.36, 8.19, 8.05 | 10.07, 9.32, 8.83, 8.49, 8.31, 8.17 | 0.79 to 0.98 |
| Block plug-in, clip at 3 times the TRUE sigma (oracle), k = 38 | 7.79 | 7.91 | above 1, only because it uses the unknown sigma |

The simulated values range from 1.2% below the closed form (TPQ, 13.48 against 13.65) to 3.2% above it (PI-LU, 8.46 against 8.20), all within two simulation standard errors; the feasible PI-LU and PI-L3 sit highest, as expected from their estimated scale (README section 2.5). RK-k is not comparable with its closed form, which is for the oracle cut-off. The averaged kernel of the paper is a truncated RQ with edge weights (`docs/math.md`, section 4). Without truncation its AVAR would be 32/3; clipping at three standard deviations with the kappa correction lowers it to 8.19 with an oracle cut-off and to 9.26 with the estimated daily MAD. So under constant volatility PI-LU (0.95 simulated) is somewhat more efficient than the corrected kernel (0.86) and about as efficient as the kernel's oracle version (0.976 against 0.977 in closed form).

**Truncation constants** (`results/q1_truncation_constants.csv`). Expected number of diffusive returns beyond the cut-off per day:

| Rule | cut-off a (sd) | per 390-return day | per 23,400-return day | share of `sigma^4` kept by clipping |
|---|---|---|---|---|
| raw MAD, c = 3 | 2.02 | 16.8 | 1007 | 0.705 |
| scaled MAD, c = 3 | 3.00 | 1.05 | 63.2 | 0.964 |
| `sqrt(2 log n)`, n = 390 / 1440 / 23,400 | 3.45 / 3.81 / 4.49 | 0.22 (n = 390) | 0.17 (n = 23,400) | 0.990 / 0.997 / 1.000 |
| Mancini `4 n^0.01` (TRQ-M) | 4.25 to 4.42 | 0.0085 | 0.23 | 0.999 to 1.000 |

The fixed rules truncate a fixed share of returns at every n, so their bias cannot vanish. The universal rule truncates fewer than one diffusive return per day at every n: 0.22, 0.20 and 0.17 per day at n = 390, 1440 and 23,400 (n times `p_beyond` in the CSV). The count tends to zero only like `1 / sqrt(pi log n)`, but the relative bias it causes is `O(sqrt(log n) / n)`.

**Jump residuals** (`results/q1_jump_residuals.csv`), in units of `sigma^4 delta` per jump:

| Case | Plug-in (k = 20) | Fourth powers (k = 1) |
|---|---|---|
| Jump clipped at a = 3, closed form | 17.55 | 27.01 |
| Same, simulation (40,000 paths, one 50 sd jump, n = 400) | 17.56 | 27.02 |
| Jump dropped (re-counted block), simulation | 0.003 | |
| Undetected jump of 2 sd, closed form / simulation | 8.73 / 8.74 | 13.33 / 13.32 |
| Undetected jump of 3 sd, closed form | 21.68 | 45.00 |

At n = 390 a clipped jump costs 4.5% of IQ for the plug-in and 6.9% for fourth powers (closed form, a = 3). With the universal cut-off clipping would cost more (6.5% and 12.0% per jump at n = 390), which is why PI-LU drops instead.

**Local scale ablation** (`results/q1_scale_ablation.csv`), constant volatility, mean of IQ estimate / IQ (simulation standard error 0.001):

| n | scale | clip a = 3 | drop a = 3 | clip `a_n` | drop `a_n` |
|---|---|---|---|---|---|
| 390 | true sigma (oracle) | 1.000 | 0.999 | 1.000 | 0.999 |
| 390 | local scaled MAD | 0.993 | 0.987 | 0.996 | 0.988 |
| 390 | two-step (used) | 1.004 | 1.020 | 1.000 | 1.006 |
| 1440 | local scaled MAD | 0.996 | 0.992 | 0.999 | 0.996 |
| 1440 | two-step (used) | 1.002 | 1.010 | 0.999 | 1.000 |

The two-step scale removes the downward bias of the MAD-only scale except for dropping at a fixed a = 3, which then overshoots by 1% to 2%. That combination is not used.

## Q2. Which estimators are consistent, and at what rate [SYNTHETIC]

Log-OU stochastic volatility (vol of vol 0.5, or 1 with an intraday U-shape of amplitude 0.6), 500 days per cell, jumps of fixed size in price units (`results/q2_q4_sim_summary.csv`, `results/q2_rates.csv`). The slope is the least-squares slope of log RMSE on log n over n = 390, 1440, 23,400; minus 0.5 is the `sqrt(n)` rate, 0 means no convergence. Figure: `results/fig1_convergence.png`.

**S5: strong stochastic volatility plus intraday pattern, no jumps.**

| Estimator | bias 390 | bias 1440 | bias 23,400 | RMSE 390 | RMSE 1440 | RMSE 23,400 | slope |
|---|---|---|---|---|---|---|---|
| RK (paper's kernel, our normalised version) | -0.209 | -0.159 | -0.144 | 0.256 | 0.184 | 0.156 | -0.11 |
| RK-k (kappa-corrected) | -0.179 | -0.128 | -0.112 | 0.236 | 0.160 | 0.128 | -0.14 |
| PI-G3 (plug-in, daily threshold) | -0.069 | -0.061 | -0.055 | 0.178 | 0.109 | 0.069 | -0.22 |
| TRQ-M (Mancini, daily BV scale) | -0.038 | -0.028 | -0.024 | 0.183 | 0.104 | 0.041 | -0.36 |
| MedRQ | -0.017 | 0.005 | 0.002 | 0.236 | 0.127 | 0.031 | -0.50 |
| TRQ-U | 0.003 | 0.002 | 0.001 | 0.193 | 0.111 | 0.027 | -0.48 |
| PI-L3 | -0.003 | 0.000 | 0.001 | 0.180 | 0.097 | 0.023 | -0.51 |
| PI-LU | -0.002 | -0.001 | 0.000 | 0.183 | 0.098 | 0.023 | -0.51 |
| KPI-L3 | -0.090 | -0.021 | 0.000 | 0.190 | 0.097 | 0.023 | -0.52 |
| PI-NN | 0.008 | 0.006 | 0.001 | 0.192 | 0.105 | 0.024 | -0.51 |

**S3: vol of vol 0.5, four medium jumps a day** (about 6 one-minute sd at n = 390, 46 at n = 23,400).

| Estimator | bias 390 | bias 1440 | bias 23,400 | RMSE 390 | RMSE 1440 | RMSE 23,400 | slope |
|---|---|---|---|---|---|---|---|
| RK | 0.130 | 0.011 | -0.043 | 0.249 | 0.091 | 0.047 | -0.38 |
| RK-k | 0.172 | 0.049 | -0.007 | 0.280 | 0.106 | 0.021 | -0.63 |
| TRQ-M | 0.119 | 0.018 | -0.005 | 0.273 | 0.103 | 0.021 | -0.62 |
| MedRQ | 0.140 | 0.143 | 0.196 | 0.471 | 2.095 | 4.377 | 0.50 |
| MinRQ | 0.117 | 0.037 | 0.000 | 0.355 | 0.271 | 0.029 | -0.64 |
| TPQ | 0.322 | 0.206 | 0.083 | 0.630 | 0.379 | 0.129 | -0.39 |
| PI-L3 (clip) | 0.136 | 0.041 | 0.001 | 0.240 | 0.094 | 0.018 | -0.62 |
| PI-LUc (clip, `a_n`) | 0.177 | 0.067 | 0.005 | 0.283 | 0.115 | 0.019 | -0.66 |
| **PI-LU (drop, `a_n`)** | **0.046** | **0.007** | -0.002 | **0.179** | **0.082** | **0.018** | -0.55 |
| PI-LU without correction | 0.129 | 0.053 | 0.011 | 0.227 | 0.100 | 0.021 | -0.57 |
| PI-NN | 0.140 | 0.035 | -0.001 | 0.320 | 0.121 | 0.019 | -0.68 |

**Reading.**
* The kappa correction removes the constant part of the bias of the daily-scale kernel (`kappa_4(3) = 0.964`; in S1, RK is minus 4.2% to minus 4.7% at every n and RK-k minus 0.6% to minus 1.1%). What remains depends on how much volatility moves within the day. In S5 the daily-scale rules (RK, RK-k, PI-G3, and to a lesser degree TRQ-M) keep a bias that does not shrink with n (RK-k minus 11% at n = 23,400), as README section 2.5 predicts for a single daily scale.
* Every local-scale estimator converges at the `sqrt(n)` rate, whatever the clip or drop choice.
* With jumps, the clipping estimators converge faster than `sqrt(n)` (slopes minus 0.6 to minus 0.7) because the residual of a clipped jump shrinks like `1/n`; at n = 390 and 1440 that residual dominates (PI-L3 plus 13.6% at n = 390). Dropping (PI-LU) removes it.
* MedRQ breaks down in S3 at n = 23,400 (RMSE 4.4): with four jumps among 23,400 returns, a few simulated days put two jumps next to each other and the median of three adjacent returns is then a jump. Its median error is 0.000 there (`median_rel_err` in the CSV), so this is a rare-event failure, not a bias. MinRQ has the same weakness at n = 390 to 1440 in S4 and S6.
* The finite-sample correction matters at practical n: without it PI-LU is 9.2% too high in S1 at n = 390 and still 1.3% at n = 23,400.
* RK in S3 has the smallest RMSE at n = 1440 among the daily-scale estimators only because two biases cancel (clipped jumps push up, clipping of diffusive returns pushes down); at n = 23,400 the downward one remains (minus 4.3%).

## Q3. Bandwidth and block length [SYNTHETIC]

**Cross-validation of eq. (17)** (`results/q3_cv_choice.csv`, `results/q3_cv_theory.csv`, `results/q3_iq_by_bandwidth.csv`). Returns clipped at 3 times the scaled daily MAD as in the paper, grid h in {5, 10, 20, 30, 40, 50, 75, 100, 150}, 300 days per cell. Share of days on which each criterion picks a given h, and the median choice:

| Scenario, n | printed eq. (17) | corrected (r^4 target) | corrected (r^2 target) | oracle spot MSE |
|---|---|---|---|---|
| S1, 390 | h = 5 on 87% (median 5) | median 150 (58% at 150) | median 150 | median 150 (63% at 150) |
| S1, 1440 | h = 5 on 97% | median 150 (84%) | median 150 (84%) | median 150 (97%) |
| S3, 390 | h = 5 on 89% | median 150 | median 150 | median 150 |
| S5, 390 | h = 5 on 58%, h = 10 on 41% | median 50 | median 50 | median 50 (r^4), 40 (r^2) |
| S5, 1440 | h = 5 on 84% | median 100 | median 100 | median 100 |

Closed form under constant volatility (Proposition 5): the expected printed CV is 74.3 at h = 5, 77.2 at h = 10, 99.2 at h = 20, 142.8 at h = 30 and 295.0 at h = 50 (units `sigma^8 delta^4`), so on the paper's grid {10, ..., 50} it always prefers h = 10, whatever the data. The corrected expected CV decreases from 75.5 (h = 10) to 69.8 (h = 50), as it should when volatility is constant.

At the oracle spot-optimal bandwidth the printed choice h = 5 has 6 to 22 times the spot MSE of the best h (`oracle_mse4_rel`). The corrected criteria pick close to the oracle.

**But the spot-optimal bandwidth is not IQ-optimal.** Relative RMSE of the IQ estimate as a function of h:

| h | KPI-L3, S5, 390 | RK-k, S5, 390 | KPI-L3, S5, 1440 | KPI-L3, S1, 1440 |
|---|---|---|---|---|
| 5 | 0.205 | 0.199 | 0.107 | 0.084 |
| 10 | 0.198 | 0.207 | 0.106 | 0.082 |
| 30 | 0.204 | 0.237 | 0.104 | 0.081 |
| 50 | 0.231 | 0.269 | 0.105 | 0.082 |
| 100 | 0.332 | 0.344 | 0.125 | 0.084 |
| 150 | 0.413 | 0.396 | 0.158 | 0.086 |

When volatility moves (S5), long windows bias the squared local variance down (Jensen) and the averaged kernel loses its first h returns, which are the busy opening of the day. So for IQ, the corrected CV (median h = 50 to 100 in S5) oversmooths: KPI-L3 at the corrected-CV choice has RMSE 0.250 at n = 390, against 0.201 at the printed choice and 0.198 at h = 10. The printed criterion's choice works well here, but because of its scale term rather than by design. For the *integrated* quantity the bandwidth should undersmooth the spot one.

**Block length of PI-LU** (`results/q3_block_length.csv`), relative RMSE (bias in brackets):

| k | S1, 390 | S5, 390 | S5, 1440 | S5, 23,400 |
|---|---|---|---|---|
| 2 | 0.164 | 0.193 (0.001) | 0.108 | 0.026 |
| 10 | 0.156 | 0.185 (0.001) | 0.100 | 0.024 |
| 20 (default at 390) | 0.154 | 0.183 (-0.002) | 0.100 | 0.023 |
| 40 | 0.150 | 0.175 (-0.011) | 0.098 | 0.023 |
| 80 | 0.148 | 0.170 (-0.040) | 0.096 | 0.023 |
| 160 | 0.145 | 0.194 (-0.127) | 0.094 | 0.023 |

The curve is flat between k = 10 and k = 80 in these smooth-volatility designs, and the default `ceil(sqrt n)` is within 1% to 8% of the best k in every cell. A quick exploratory check (not saved) showed that the squared-return CV also picks long windows on the burst model of Q6, so CV cannot detect minute-scale bursts; Q6 shows why that matters.

## Q4. Full comparison and confidence intervals [SYNTHETIC]

Relative RMSE, 500 days per cell (bias and standard deviation in `results/q2_q4_sim_summary.csv`). S1 no jumps, S2 four small jumps (about 2 sd at n = 390), S3 medium (6 sd), S4 large (20 sd), S5 strong SV plus pattern, S6 = S5 plus medium jumps.

**n = 390**

| Estimator | S1 | S2 | S3 | S4 | S5 | S6 |
|---|---|---|---|---|---|---|
| RQ | 0.167 | 0.631 | 28.2 | 3017 | 0.197 | 31.2 |
| TPQ | 0.189 | 0.213 | 0.630 | 3.54 | 0.239 | 0.782 |
| MedRQ | 0.186 | 0.208 | 0.471 | 177 | 0.236 | 19.3 |
| MinRQ | 0.217 | 0.243 | 0.355 | 184 | 0.280 | 23.0 |
| TRQ-M | 0.166 | 0.214 | 0.273 | 0.353 | 0.183 | 0.296 |
| RK | 0.156 | 0.156 | 0.249 | 0.312 | 0.256 | 0.198 |
| RK-k | 0.156 | 0.170 | 0.280 | 0.350 | 0.236 | 0.197 |
| TRQ-U | 0.167 | 0.196 | 0.218 | 0.186 | 0.193 | 0.268 |
| PI-G3 | 0.149 | 0.158 | 0.222 | 0.261 | 0.178 | 0.191 |
| PI-L3 | 0.151 | 0.167 | 0.240 | 0.269 | 0.180 | 0.283 |
| PI-LUc | 0.150 | 0.172 | 0.283 | 0.341 | 0.178 | 0.330 |
| PI-LU | 0.154 | 0.165 | 0.179 | 0.163 | 0.183 | 0.217 |
| KPI-L3 | 0.154 | 0.167 | 0.238 | 0.265 | 0.190 | 0.215 |
| PI-NN | 0.162 | 0.208 | 0.320 | 0.236 | 0.192 | 0.375 |

**n = 1440**

| Estimator | S1 | S2 | S3 | S4 | S5 | S6 |
|---|---|---|---|---|---|---|
| RQ | 0.091 | 1.35 | 102 | 12917 | 0.111 | 95.1 |
| TPQ | 0.100 | 0.125 | 0.379 | 2.99 | 0.130 | 0.383 |
| MedRQ | 0.106 | 0.111 | 2.10 | 93.4 | 0.127 | 1.19 |
| MinRQ | 0.119 | 0.127 | 0.271 | 33.5 | 0.144 | 1.42 |
| TRQ-M | 0.090 | 0.105 | 0.103 | 0.090 | 0.104 | 0.108 |
| RK | 0.091 | 0.084 | 0.091 | 0.089 | 0.184 | 0.154 |
| RK-k | 0.082 | 0.091 | 0.106 | 0.107 | 0.160 | 0.134 |
| TRQ-U | 0.090 | 0.098 | 0.098 | 0.086 | 0.111 | 0.115 |
| PI-G3 | 0.078 | 0.082 | 0.090 | 0.091 | 0.109 | 0.106 |
| PI-L3 | 0.078 | 0.086 | 0.094 | 0.096 | 0.097 | 0.120 |
| PI-LUc | 0.078 | 0.094 | 0.115 | 0.120 | 0.096 | 0.138 |
| PI-LU | 0.078 | 0.081 | 0.082 | 0.078 | 0.098 | 0.106 |
| KPI-L3 | 0.078 | 0.086 | 0.095 | 0.097 | 0.097 | 0.111 |
| PI-NN | 0.085 | 0.119 | 0.121 | 0.088 | 0.105 | 0.129 |

**n = 23,400**: all local estimators are at 0.018 to 0.026 in every scenario; RK 0.046 to 0.156, RK-k 0.021 to 0.128, PI-G3 up to 0.069 and TRQ-M up to 0.041 (the S5 and S6 columns), TPQ up to 0.58 (S4) and MedRQ 4.38 in S3 (adjacent jumps).

PI-LU is best or within 0.01 of the best in 11 of the 12 cells at n = 390 and 1440; the exception is S6 at n = 390 (0.217 against 0.191). There PI-G3, RK-k and RK do well (0.191, 0.197, 0.198) because the downward bias from the daily scale offsets the upward bias from the jumps; at n = 23,400 they are biased by minus 5%, minus 10% and minus 14%.

**95% confidence intervals** (`results/q4_coverage.csv`), coverage of the true IQ, 500 days per cell:

| Scenario | PI-LU 390 / 1440 / 23,400 | PI-LU log-scale | PI-L3 (clip) | PI-NN |
|---|---|---|---|---|
| S1 | 0.944 / 0.944 / 0.938 | 0.946 / 0.956 / 0.946 | 0.948 / 0.952 / 0.938 | 0.946 / 0.948 / 0.948 |
| S2 | 0.952 / 0.948 / 0.968 | 0.946 / 0.952 / 0.964 | 0.958 / 0.944 / 0.956 | 0.960 / 0.942 / 0.966 |
| S3 | 0.942 / 0.936 / 0.970 | 0.926 / 0.940 / 0.972 | 0.904 / 0.924 / 0.970 | 0.946 / 0.936 / 0.962 |
| S4 | 0.930 / 0.952 / 0.950 | 0.948 / 0.960 / 0.952 | 0.880 / 0.926 / 0.956 | 0.942 / 0.954 / 0.942 |
| S5 | 0.926 / 0.950 / 0.960 | 0.936 / 0.952 / 0.958 | 0.928 / 0.948 / 0.960 | 0.932 / 0.950 / 0.958 |
| S6 | 0.944 / 0.940 / 0.948 | 0.930 / 0.944 / 0.950 | 0.946 / 0.938 / 0.948 | 0.940 / 0.950 / 0.936 |

With 500 days the simulation standard error of a coverage near 0.95 is about 0.01. PI-LU covers within 0.926 to 0.972 everywhere. The clipping version undercovers with jumps at n = 390 (0.880 in S4) because of the clipped-jump residual. The mean half-width of the PI-LU interval is 29% to 36% of IQ at n = 390, 15% to 19% at 1440 and 4% to 5% at 23,400; the ratio of the estimated to the actual standard deviation is 0.91 to 1.05.

**Robust `int sigma^8` (the paper's "eighth power" future-work item)**, PI-LU, relative bias / RMSE: S1 0.024 / 0.397 at n = 390, minus 0.001 / 0.045 at 23,400; S5 minus 0.016 / 0.611 and minus 0.001 / 0.069. With medium jumps at n = 390 it is biased up (0.128 / 0.495), since a surviving jump enters to the eighth power. It is accurate enough for interval widths at n of 1440 and above. PI-NN's eighth-power estimate is much less robust (bias 1.65 in S3 at n = 390).

## Q5. Real data: Binance BTCUSDT and DOGEUSDT, 1 July to 31 December 2025 [REAL]

184 UTC days per coin, 1440 one-minute log returns per day (one missing value in the whole sample). Daily IQ for every estimator (`results/q5_real_daily.csv`). The true IQ is unknown, so only agreement can be measured. Median over days of the ratio to MedRQ and MinRQ of the same day, median absolute log ratio to MedRQ, rank correlation (Spearman) with MedRQ, and the same ratio to MedRQ at 5 minutes (`results/q5_real_summary.csv`). Figure: `results/fig2_real_ratio_to_medrq.png`.

| Estimator | / MedRQ, BTC | / MedRQ, DOGE | / MinRQ, BTC | / MinRQ, DOGE | abs log ratio, BTC | abs log ratio, DOGE | Spearman BTC | Spearman DOGE | / MedRQ at 5 min, BTC | DOGE |
|---|---|---|---|---|---|---|---|---|---|---|
| RQ | 1.26 | 1.32 | 1.32 | 1.33 | 0.27 | 0.36 | 0.959 | 0.921 | 1.25 | 1.40 |
| TPQ | 0.99 | 1.03 | 1.01 | 1.03 | 0.13 | 0.13 | 0.987 | 0.987 | 1.02 | 1.06 |
| MinRQ | 1.03 | 1.05 | 1.00 | 1.00 | 0.09 | 0.10 | 0.995 | 0.994 | 1.02 | 1.04 |
| TRQ-M | 0.63 | 0.66 | 0.65 | 0.68 | 0.46 | 0.41 | 0.927 | 0.890 | 0.82 | 0.82 |
| RK (paper's kernel, our normalised version) | 0.33 | 0.42 | 0.33 | 0.40 | 1.10 | 0.86 | 0.865 | 0.814 | 0.56 | 0.64 |
| RK-k | 0.35 | 0.44 | 0.34 | 0.42 | 1.06 | 0.82 | 0.865 | 0.814 | 0.58 | 0.66 |
| PI-G3 | 0.38 | 0.46 | 0.39 | 0.46 | 0.97 | 0.77 | 0.879 | 0.837 | 0.62 | 0.71 |
| TRQ-U | 0.76 | 0.80 | 0.79 | 0.77 | 0.28 | 0.24 | 0.928 | 0.887 | 0.83 | 0.86 |
| PI-L3 | 0.70 | 0.71 | 0.73 | 0.72 | 0.35 | 0.34 | 0.939 | 0.908 | 0.86 | 0.89 |
| PI-LUc | 0.75 | 0.78 | 0.79 | 0.77 | 0.29 | 0.25 | 0.946 | 0.922 | 0.89 | 0.92 |
| PI-LU | 0.66 | 0.69 | 0.69 | 0.67 | 0.41 | 0.37 | 0.931 | 0.888 | 0.80 | 0.80 |
| KPI-L3 | 0.72 | 0.73 | 0.75 | 0.74 | 0.33 | 0.31 | 0.941 | 0.911 | 0.80 | 0.83 |
| **PI-NN** | **0.95** | **0.96** | **0.95** | **0.97** | **0.16** | **0.17** | **0.975** | **0.965** | **0.98** | **1.02** |

**Answer to the question asked.** Yes: every local-scale estimator agrees with MedRQ and MinRQ day by day much better than the estimators with one daily MAD threshold (median ratio 0.66 to 0.80 against 0.33 to 0.46; absolute log ratio 0.24 to 0.41 against 0.77 to 1.10; higher rank correlation). The direction of the E5 finding on Bitstamp (MAD-truncated estimators 3 to 7 times lower than the nearest-neighbour ones) is reproduced on Binance, with a smaller factor: RK and RK-k are 2.3 to 3.0 times lower than MedRQ at the median. Which estimator is closer to the truth cannot be read off real data. The local plug-ins with a two-hour scale window still sit 20% to 34% below MedRQ, and only PI-NN closes the gap; Q6 asks which side of the gap a plausible model supports.

**Diagnostics** (`results/q5_real_diagnostics.csv`). Returns divided by their local scale have a median daily kurtosis of 6.3 (BTC) and 5.2 (DOGE), 3.7 and 3.5 even after trimming at `sqrt(2 log n)` (2.98 for a trimmed normal). 1.5% and 1.1% of returns lie beyond 3 local sd (0.27% for a normal) and 0.62% and 0.42% beyond 3.81 (0.014%). Zero returns: 2.4% and 2.6% of minutes.

**Sampling frequency** (`results/q5_real_frequency.csv`). Median ratio of the 1-minute estimate to the 5-minute estimate of the same day: MedRQ 1.33 (BTC) and 1.45 (DOGE), MinRQ 1.43 and 1.47, TPQ 1.30 and 1.42, RQ 1.38 and 1.39; PI-LU 1.14 and 1.18, PI-NN 1.20 and 1.31; RK 0.88 and 0.99. If returns were locally Gaussian at both frequencies, every consistent estimator would give the same answer at 1 and 5 minutes. The fourth-power estimators rise by a third or more.

**Sensitivity of the plug-in on the real data** (`results/q5_real_sensitivity.csv`, local scale, drop, BTC first then DOGE): with the cut-off fixed at 3.81, shortening the block raises the estimate steadily, from 0.58 / 0.61 of MedRQ at k = 160 to 0.66 / 0.69 at k = 38 and 0.77 / 0.81 at k = 2. With k = 38, raising the cut-off from 3 to 3.81, 5, 8 and no truncation gives 0.59, 0.66, 0.75, 0.82 and 0.87 of MedRQ (BTC). In the slow-volatility simulations neither knob mattered much (Q3); on real data both do. So real volatility varies within blocks of tens of minutes, and a sizeable share of the quarticity sits in returns beyond 3.8 local sd.

## Q6. Bursts or jumps? [REAL + SYNTHETIC]

**Clustering of large returns [REAL]** (`results/q6_burst_diagnostics.csv`; z = return / local scale):

| | acf of `min(z^2, 25)`, lag 1 | lag 5 | lag 10 | trimmed kurtosis | share beyond 3 | exceedances with a neighbour beyond 2 |
|---|---|---|---|---|---|---|
| BTCUSDT | 0.114 | 0.034 | 0.014 | 3.69 | 1.46% | 65% of 5024 |
| DOGEUSDT | 0.099 | 0.026 | 0.002 | 3.53 | 1.08% | 44% of 1253 |
| Model without a fast factor [SYNTHETIC] | -0.008 | -0.006 | -0.010 | 2.93 | 0.21% | 0% of 16 |
| Chance level, independent Gaussian returns (closed form) | 0 | 0 | 0 | 2.98 | 0.27% | 8.9% |
| Calibrated fast factor (sd 0.75, half-life 3 min) [SYNTHETIC] | 0.120 | 0.032 | -0.005 | 3.68 | 1.25% | 41% of 606 |

Isolated jumps would leave squared returns uncorrelated and give exceedances ordinary neighbours: the chance that at least one of the two neighbours of a return lies beyond 2 sd is `1 - 0.9545^2 = 8.9%` for independent Gaussian returns. The 0 of 16 of the model without a fast factor is consistent with that (1.4 expected); the real 44% to 65% of 1253 to 5024 exceedances is far above it. Large moves come in bursts that fade within about 5 to 10 minutes. This agrees with the high-frequency literature: Christensen, Oomen and Podolskij (2014) find with tick data on equities and currencies that lower-frequency jump measures attribute bursts of volatility to the jump component, and Bajgrowicz, Scaillet and Treccani (2016) find that most news produces bursts of volatility rather than jumps (both read from their abstracts only). The calibration grid (fast-factor sd 0.5, 0.75, 1; half-life 1, 3, 10 minutes; 100 days each) picked sd 0.75 and half-life 3 minutes as closest to the average of the two coins.

**The calibrated model reproduces the real estimator gaps [SYNTHETIC vs REAL]** (`results/q6_fastvol_frequency.csv`, n = 1440, 300 days):

| Estimator | 1-min / 5-min, model | 1-min / 5-min, real BTC / DOGE | / MedRQ at 1 min, model | / MedRQ at 1 min, real BTC / DOGE |
|---|---|---|---|---|
| MedRQ | 1.41 | 1.33 / 1.45 | 1 | 1 |
| TPQ | 1.41 | 1.30 / 1.42 | 1.00 | 0.99 / 1.03 |
| RK | 1.00 | 0.88 / 0.99 | 0.53 | 0.33 / 0.42 |
| TRQ-U | 1.09 | 1.22 / 1.27 | 0.80 | 0.76 / 0.80 |
| PI-L3 | 1.03 | 1.10 / 1.14 | 0.71 | 0.70 / 0.71 |
| PI-LU | 1.04 | 1.14 / 1.18 | 0.68 | 0.66 / 0.69 |
| KPI-L3 | 1.19 | 1.18 / 1.22 | 0.71 | 0.72 / 0.73 |
| PI-NN | 1.21 | 1.20 / 1.31 | 0.95 | 0.95 / 0.96 |

The ratios of the local estimators to MedRQ match to within 0.03, the frequency ratios of the fourth-power estimators to within 0.1. The daily-scale estimators are lower on real data than in the model (RK 0.33 to 0.42 against 0.53), which suggests a stronger intraday pattern or slow volatility in the real data than in the model. In the model, the 5-minute estimates are far below the truth (median 5-minute MedRQ / IQ = 0.57), so agreement with a 5-minute benchmark is not evidence of accuracy.

**Which estimators recover IQ in that model [SYNTHETIC]** (`results/q6_fastvol_estimators.csv`, 300 days per cell; median relative error, then RMSE):

| Estimator | median err 390 / 1440 / 23,400 | RMSE 390 / 1440 / 23,400 | with jumps: median err | with jumps: RMSE |
|---|---|---|---|---|
| RQ | -0.25 / -0.09 / -0.01 | 0.38 / 0.25 / 0.07 | 2.1 / 8.4 / 176 | 11.9 / 54.7 / 822 |
| TPQ | -0.40 / -0.20 / -0.02 | 0.42 / 0.28 / 0.08 | -0.27 / -0.12 / 0.02 | 0.43 / 0.51 / 0.17 |
| MedRQ | -0.40 / -0.20 / -0.02 | 0.42 / 0.29 / 0.08 | -0.31 / -0.14 / -0.02 | 0.89 / 1.63 / 1.28 |
| MinRQ | -0.37 / -0.20 / -0.03 | 0.41 / 0.34 / 0.09 | -0.30 / -0.12 / -0.02 | 1.03 / 0.80 / 0.10 |
| TRQ-M | -0.39 / -0.35 / -0.31 | 0.41 / 0.37 / 0.33 | -0.30 / -0.32 / -0.30 | 0.35 / 0.33 / 0.32 |
| RK | -0.58 / -0.59 / -0.57 | 0.59 / 0.59 / 0.59 | -0.54 / -0.56 / -0.57 | 0.54 / 0.56 / 0.58 |
| RK-k | -0.57 / -0.57 / -0.56 | 0.57 / 0.57 / 0.57 | -0.52 / -0.54 / -0.56 | 0.53 / 0.55 / 0.57 |
| PI-L3 | -0.45 / -0.44 / -0.29 | 0.46 / 0.44 / 0.30 | -0.38 / -0.41 / -0.29 | 0.41 / 0.41 / 0.30 |
| PI-LU | -0.49 / -0.47 / -0.26 | 0.49 / 0.46 / 0.27 | -0.47 / -0.44 / -0.27 | 0.48 / 0.44 / 0.27 |
| KPI-L3 | -0.51 / -0.44 / -0.20 | 0.52 / 0.44 / 0.21 | -0.45 / -0.41 / -0.20 | 0.46 / 0.41 / 0.21 |
| **PI-NN** | -0.36 / -0.23 / -0.04 | **0.39 / 0.28 / 0.07** | -0.30 / -0.20 / -0.05 | 0.36 / **0.26 / 0.07** |

Figure: `results/fig3_fastvol_bias.png`.

**Reading.**
* Volatility bursts of a few minutes are part of `sigma_t`, so they are part of IQ. Any estimator that truncates against a scale estimated over hours removes them as if they were jumps, and any long block averages them away. In that model our version of the paper's kernel (RK, RK-k) misses 56% to 59% of IQ at every n (median error), and the local plug-ins PI-L3, PI-LU and KPI-L3 20% to 51%. Their bias shrinks only slowly with n, because their windows are fixed multiples of `sqrt(n)` returns, which is still 28 minutes at n = 23,400.
* The fourth-power nearest-neighbour estimators (MedRQ, MinRQ, TPQ) keep bursts because a large return inside a burst has large neighbours. But they break down when jumps are added (RMSE 0.51 to 1.63 at n = 1440).
* PI-NN combines the two ideas. Its neighbour guard keeps bursts and its threshold still removes isolated jumps. It has the lowest RMSE of the robust estimators in five of the six cells; at n = 1440 without jumps it is nearly tied with TPQ (0.277 against 0.280), and at n = 390 with jumps TRQ-M is slightly lower (0.352 against 0.360). Its block-length table (`results/q6_fastvol_block_length.csv`) shows the guard is what matters: at n = 1440 and k = 5 the bias is minus 19% with the guard against minus 34% without it.
* At one-minute sampling every estimator is still biased down in this model, even those that truncate nothing or keep bursts (median error at n = 1440: RQ minus 9%, MedRQ, MinRQ and TPQ minus 20%, PI-NN minus 23%), because some burst variation lies within a single minute. Only finer sampling resolves it (n = 23,400: PI-NN minus 4%, MedRQ minus 2%).

## Negative results and caveats

* **The first proposal (PI-LU) is not the best on real-data-like data.** Q6 shows that its scale window and block length are too long for minute-scale volatility bursts. PI-NN was designed after seeing Q5 and Q6, and its settings (six neighbours, `k = n^(1/4)`) were chosen in an exploratory run on the burst model and on S3 and S4 with other seeds (iteration 3). Its numbers in Q4 to Q6 are therefore not strictly out of sample.
* **PI-NN is weaker than PI-LU with isolated jumps at small n** (S3, n = 390: RMSE 0.320 against 0.179; S6: 0.375 against 0.217). The noisy six-neighbour scale sometimes lets a jump through, and short blocks amplify it (`2 b^2 + b^4 / (k + 2)`). Its eighth-power estimate is not robust (bias 1.65 in S3 at n = 390), so its intervals rely on a fragile ingredient. Coverage was still 0.932 to 0.966 in the grid (plain intervals).
* **The burst interpretation is a model, not a proof.** The two-factor model matches seven diagnostics and the estimator ratios, but small jumps that cluster with trading activity, or microstructure effects, could produce similar signatures. Christensen, Oomen and Podolskij (2014) also report that some tick-level jumps are induced by microstructure. The global-scale estimators fit worse (RK 0.53 in the model against 0.33 to 0.42 real).
* **Cross-validation, corrected or not, targets the spot quantity.** It oversmooths for IQ when volatility moves, and it cannot see minute-scale bursts. The printed criterion's preference for the smallest h comes from its scale term, and its good IQ performance in S5 is a side effect of that, not of the criterion's design.
* **The fixed cut-off with a local scale (PI-L3) is consistent in the model but relies on local Gaussianity.** On real data, where local kurtosis is 5 to 6, its Gaussian correction does not undo the truncation.
* **No central limit theorem is proved** for the feasible estimators. The coverage results are numerical evidence only.
* One venue (Binance), two coins, six months. No microstructure noise in the simulations.

## Ranking (novelty, validity, relevance, feasibility, interest)

**Novelty: moderate.** The plug-in principle and its efficiency (Jacod and Rosenbaum 2013), locally Gaussian block corrections (Mykland and Zhang 2009), vanishing thresholds (Mancini 2009) and their optimal order (Figueroa-Lopez and Nisen 2013), nearest-neighbour truncation (Andersen, Dobrev and Schaumburg 2012) and kernel plug-ins (Figueroa-Lopez, Pang and Wu 2025) are all known. That apparent jumps at lower frequencies are often bursts of volatility is documented by Christensen, Oomen and Podolskij (2014) and Bajgrowicz, Scaillet and Treccani (2016). Our web searches (snippets only, not exhaustive) found none of the following:
* exact finite-sample corrections for clipped or dropped block estimators of any power of volatility, via cumulants;
* a statement of when the MAD rule recovers IQ: a known constant under constant volatility, removed by the kappa correction; a bias that does not shrink with one daily scale when volatility moves within the day; consistency with a local scale and a correction, or with a `sqrt(2 log n)` cut-off;
* the clipped and undetected jump residuals of plug-in against fourth-power estimators (`b^2` against `b^4`);
* the closed-form scale artefact in the cross-validation criterion as printed;
* the neighbour-guarded threshold;
* the consequence of volatility bursts for truncated IQ estimators on one-minute crypto data: scale-truncated estimators come out at 0.33 to 0.73 times the nearest-neighbour ones, and a model calibrated to the clustering reproduces those ratios.

Each piece alone is incremental; together they could make a coherent estimation paper.

**Validity: good for what is claimed.** The elementary propositions are proved and checked to three digits by simulation. Consistency and the `sqrt(n)` rate are shown numerically over three frequencies and six designs. Coverage is 93% to 97%. The CLT for the feasible estimators is not proved. PI-NN's settings were chosen after looking at the data.

**Practical relevance: high.** Every jump test and every confidence interval for daily variance needs IQ. On real minute data the choice of estimator moves IQ by a factor of 3, so this matters.

**Feasibility: high.** The code exists and runs in about 10 minutes. A paper needs a CLT for the local-threshold plug-in, which should follow from Jacod and Rosenbaum (2013) with a data-driven threshold, plus a careful treatment of the neighbour guard. That is a few months of work for a statistician who works on estimation.

**Likely interest to the paper's author: high.** It keeps the author's central idea (robust-scale truncation plus local smoothing), states conditions under which it recovers IQ, proposes for discussion an aggregation that is somewhat more efficient in our simulations (0.95 against 0.86 for the corrected kernel under constant volatility), touches both future-work items of the paper (data-driven bandwidth, higher powers of volatility), and adds a real-data observation on how truncated estimators behave on crypto minute data. Two points need care in a conversation: our kernel results are for our normalised, scaled-MAD version, not the formula as printed; and on real data we only measure ratios between estimators (RK is 0.33 to 0.42 times MedRQ), so "too low" holds only if the burst model is right. Framed as questions ("what does the MAD idea need when volatility moves within the day?"), it is a natural basis for a joint follow-up.

**Overall ranking among the estimation directions: first.** PI-LU is the clean, provable proposal for the classical model. PI-NN plus the burst consequence is the more interesting one for real data, but needs more theory and an out-of-sample check.
