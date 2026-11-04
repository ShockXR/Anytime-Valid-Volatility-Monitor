# Results: correlation monitor (alternative A7)

Every number here was produced by `run_experiment.py` (fixed seeds) and is stored in `results/`. Labels: **[SYNTHETIC]** for simulations, **[REAL]** for Binance spot 1-minute bars (BTC, ETH, SOL, DOGE against USDT, 1 January 2024 to 31 December 2025, 1,052,640 minutes per coin, no missing minutes, no duplicates). Settings are those of E4 of the main report unless stated: per-step units, calibration from the same process as the pre-change data, `A = 7800` for 390-step synthetic days, each side at `2A`, bets at fractions (0.05, 0.1, 0.2, 0.35, 0.5, 0.75) of the largest admissible bet.

Detectors compared (details in `README.md`):

* **sign**: betting e-detector on `sign(r1) sign(r2)`, bounded in [-1, 1], scale-free (recommended).
* **product c3** (and **c2**): betting e-detector on the product of standardised returns truncated at `c` robust standard deviations, bounded in [-c^2, c^2].
* **gaussian LR**: Shiryaev-Roberts mixture of bivariate Gaussian likelihood ratios for `rho0 +/- (0.15, 0.3, 0.45)`, valid only under exact normality.
* **naive CUSUM**: Page CUSUM on the raw, untruncated product of standardised returns, drift 0.15, threshold `h = 24.75` tuned on Gaussian `rho = 0.3` pairs to 2.5% per side per 390-step day (`q1_naive_threshold.csv`).

## Q0. Theory table [SYNTHETIC] (`q0_theory.csv`)

Gaussian pairs, 400,000 draws per state. `m0`, `m1`: mean of the statistic before and after the change. `lam_max`: largest admissible bet. "Opt bet": unconstrained growth-optimal bet `(m1 - m0) / (var + (m1 - m0)^2)`. "Grid bet": best bet on the grid actually used. "Cold" delay: `log(2A) / growth`. "SR" delay: `(log(2A / 500) + log 6) / growth`, which allows for the Shiryaev-Roberts head start (the statistic is already near 500 at a change at step 500) and for the best bet carrying weight 1/6 in the mixture. This last column is the one to compare with Q2.

| Change | Statistic | m0 | m1 | lam_max | Opt bet | Grid bet | Growth per step | Cold delay | SR delay |
|---|---|---|---|---|---|---|---|---|---|
| 0.3 to 0.6 | sign | 0.193 | 0.412 | 0.838 | 0.249 | 0.168 | 0.0239 | 403 | 219 |
| 0.3 to 0.6 | product c3 | 0.296 | 0.595 | 0.108 | 0.214 | 0.081 | 0.0202 | 477 | 259 |
| 0.3 to 0.6 | product c2 | 0.271 | 0.547 | 0.234 | 0.268 | 0.176 | 0.0353 | 274 | 148 |
| 0.6 to 0.3 | sign | 0.410 | 0.197 | 1.694 | 0.212 | 0.339 | 0.0234 | 413 | 224 |
| 0.6 to 0.3 | product c3 | 0.597 | 0.300 | 0.119 | 0.255 | 0.089 | 0.0215 | 450 | 244 |
| 0.6 to 0.3 | product c2 | 0.548 | 0.275 | 0.290 | 0.288 | 0.217 | 0.0335 | 288 | 156 |
| 0.0 to 0.3 | sign | 0.001 | 0.193 | 0.999 | 0.192 | 0.200 | 0.0185 | 522 | 283 |
| 0.0 to 0.3 | product c3 | 0.001 | 0.296 | 0.111 | 0.257 | 0.083 | 0.0210 | 461 | 250 |
| 0.0 to 0.3 | product c2 | 0.001 | 0.271 | 0.250 | 0.288 | 0.188 | 0.0360 | 268 | 145 |
| 0.3 to 0.0 | sign | 0.193 | -0.001 | 1.239 | 0.187 | 0.248 | 0.0181 | 535 | 290 |
| 0.3 to 0.0 | product c3 | 0.297 | 0.001 | 0.115 | 0.276 | 0.086 | 0.0216 | 448 | 243 |
| 0.3 to 0.0 | product c2 | 0.272 | 0.001 | 0.268 | 0.295 | 0.201 | 0.0358 | 270 | 146 |
| 0.5 to 0.9 | sign | 0.333 | 0.713 | 0.750 | 0.596 | 0.375 | 0.0914 | 106 | 57 |
| 0.5 to 0.9 | product c3 | 0.495 | 0.894 | 0.105 | 0.211 | 0.079 | 0.0266 | 363 | 197 |
| 0.5 to 0.9 | product c2 | 0.455 | 0.825 | 0.224 | 0.290 | 0.168 | 0.0479 | 202 | 109 |
| 0.6 to 0.0 | sign | 0.409 | 0.001 | 1.692 | 0.350 | 0.592 | 0.0873 | 111 | 60 |
| 0.6 to 0.0 | product c3 | 0.595 | 0.000 | 0.119 | 0.445 | 0.089 | 0.0480 | 201 | 109 |
| 0.6 to 0.0 | product c2 | 0.546 | 0.000 | 0.290 | 0.479 | 0.217 | 0.0942 | 102 | 56 |

**Reading.** With `c = 3` the product's bound of 9 caps the bet at about 0.11, half the growth-optimal value, so the product e-detector is throttled by the very truncation that makes it valid. The sign statistic has no such cap (its bound is 1) and is predicted to be slightly faster than product c3 for moderate changes and about 2 to 3 times faster for large ones. Tightening the truncation to `c = 2` makes the product the fastest of the three on Gaussian data, at the cost of a larger bias and the volatility sensitivity shown in Q1.

## Q1. Validity: false alarms with no correlation change [SYNTHETIC] (`q1_validity.csv`, `q1_true_m0.csv`)

2000 paths per row, 2000-step calibration window per path, 5-day horizon (1950 steps). Bounds for the e-detectors: 5% within 1 day, 25% within 5 days. 95% Wilson intervals are about plus or minus 1 point at 5% and plus or minus 2 points at 25%. In the three volatility-change rows the correlation is 0.3 throughout and only the volatility changes at step 500, so every alarm is a false alarm. "Share up" is the share of alarms that said "dependence increased".

### Q1a. Estimated baseline, band 0

| Scenario | sign: 1 day / 5 days (share up) | product c3: 1 day / 5 days (share up) | gaussian LR: 1 day / 5 days (share up) | naive CUSUM: 1 day / 5 days (share up) |
|---|---|---|---|---|
| Gaussian rho 0.3 | 2.9% / 23.0% (53%) | 0.0% / 11.3% (46%) | 5.0% / 32.0% (49%) | 3.2% / 18.4% (83%) |
| Student t4 rho 0.3 | 2.4% / 23.1% (48%) | 0.2% / 15.7% (46%) | **79.8% / 100%** (17%) | **96.2% / 100%** (56%) |
| Jumps rho 0.3 (0.5% per step, 8 sd, independent) | 3.1% / 24.1% (49%) | 0.0% / 12.4% (47%) | **57.0% / 98.3%** (3%) | **22.4% / 74.0%** (61%) |
| Co-jumps rho 0.3 (plus 0.2% simultaneous, independent signs) | 2.3% / 21.5% (51%) | 0.0% / 12.1% (45%) | **71.5% / 99.6%** (8%) | **56.5% / 97.7%** (50%) |
| Gaussian rho 0.6 | 2.3% / 22.0% (49%) | 0.0% / 13.6% (47%) | 6.6% / 27.8% (51%) | 11.0% / 44.3% (87%) |
| Volatility x2 in asset 1 at step 500, rho fixed | 2.6% / 23.6% (48%) | 0.0% / **98.3%** (100% up) | 5.1% / **99.9%** (3% up) | 3.8% / **100%** (98% up) |
| Volatility x0.5 in asset 1 at step 500, rho fixed | 3.1% / 23.8% (52%) | 0.0% / **100%** (0% up) | 4.6% / **79.6%** (95% up) | 2.7% / **48.2%** (7% up) |
| Volatility x2 in both at step 500, rho fixed | 3.0% / 23.7% (52%) | 0.0% / **100%** (100% up) | 4.9% / **100%** (5% up) | 3.6% / **100%** (94% up) |
| GARCH(1,1) clustering in each asset (alpha 0.05, beta 0.94, half-life about 70 steps), rho fixed | 2.9% / 23.9% (51%) | 2.1% / **63.6%** (45%) | **25.4% / 83.7%** (40%) | **33.3% / 87.2%** (82%) |

### Q1b. Estimated baseline, band 0.05 (e-detectors only)

| Scenario | sign: 1 day / 5 days | product c3: 1 day / 5 days |
|---|---|---|
| Gaussian rho 0.3 | 0.5% / 3.6% | 0.0% / 0.4% |
| Student t4 rho 0.3 | 0.4% / 3.3% | 0.0% / 2.5% |
| Jumps rho 0.3 | 0.7% / 4.2% | 0.0% / 0.5% |
| Co-jumps rho 0.3 | 0.5% / 3.0% | 0.0% / 0.6% |
| Gaussian rho 0.6 | 0.3% / 2.6% | 0.0% / 0.8% |
| Volatility x2 in asset 1 | 0.4% / 3.7% | 0.0% / 89.4% |
| Volatility x0.5 in asset 1 | 0.6% / 4.0% | 0.0% / 96.9% |
| Volatility x2 in both | 0.5% / 4.5% | 0.0% / 100% |
| GARCH(1,1) clustering in each asset | 0.6% / 4.4% | 0.9% / 26.9% |

### Q1c. True null mean (2 million draws), band 0: the guarantee without baseline error

| Scenario | true m0 (sign) | sign: 1 day / 5 days | true m0 (product) | product c3: 1 day / 5 days |
|---|---|---|---|---|
| Gaussian rho 0.3 | 0.194 | 2.0% / 15.4% | 0.298 | 0.0% / 6.2% |
| Student t4 rho 0.3 | 0.194 | 2.3% / 16.8% | 0.383 | 0.0% / 9.1% |
| Jumps rho 0.3 | 0.192 | 1.9% / 16.4% | 0.293 | 0.0% / 5.5% |
| Co-jumps rho 0.3 | 0.192 | 1.6% / 15.5% | 0.291 | 0.0% / 6.3% |
| Gaussian rho 0.6 | 0.410 | 1.8% / 15.3% | 0.596 | 0.0% / 11.3% |

**Reading.**

1. The sign monitor is inside its bound in every row (2.3% to 3.1% per day against 5%; 21.5% to 24.1% per 5 days against 25%), and the rate is the same for Gaussian, Student t4 and jumpy pairs, as a nonparametric guarantee should give. With the true null mean the 5-day rate drops to 15% to 17%: the gap to 23% is baseline estimation error from a 2000-step window, the same effect as E4d in the main report, and a band of 0.05 removes it (3% to 4.5% per 5 days).
2. The sign monitor **does not react to a volatility change** (23.6% to 23.8% per 5 days, identical to the no-change rows). The product monitor raises a "dependence up" alarm on 98% of paths when the volatility of one asset doubles and a "dependence down" alarm on 100% when it halves; the band does not help (89% to 100%). The Gaussian LR also alarms on almost every path, and in the **opposite direction** (97% of its alarms say "down" when asset 1 becomes more volatile): for a large idiosyncratic move the Gaussian quadratic form favours a smaller correlation. This is Proposition 2 of the README made visible.
3. The Gaussian LR and the naive CUSUM are unusable under fat tails or jumps (57% to 96% false alarms per day). Note the Student t4 pairs here are elliptical, so their sign mean is exactly `(2 / pi) arcsin(0.3)`, the same as for the Gaussian pairs (0.194 in Q1c), while the product's mean moves from 0.298 to 0.383: the product statistic also confuses a change of tail shape with a change of correlation.
4. The product e-detector is very conservative (0% per day, 11% to 16% per 5 days with band 0) because its bets are capped (Q0).
5. Persistent volatility clustering is the case in which the main package's own strict null fails (E9, Q9c: 77% to 100% of GARCH days with an alarm at band 0). The sign correlation monitor is untouched by it (2.9% per day, 23.9% per 5 days, the same as under constant volatility), because its null is invariant to any predictable volatility path. The product monitor's 5-day rate rises from 11% to 64% (27% with band 0.05), the Gaussian LR's to 84% and the naive CUSUM's to 87%: with clustering, products of returns have a conditional mean that moves with the two variances even though the correlation never moves.

## Q2. Detection delay [SYNTHETIC] (`q2_delay.csv`, `fig1_synthetic_delays.png`)

1000 paths, change at step 500, 3500 post-change steps, band 0. "FA before" is the share of paths that alarmed before the change (those paths do not count in the delay). The share of post-change alarms in the correct direction is 0.999 or more for the sign and product e-detectors in every row; it drops to 0.91 for the Gaussian LR (jumps) and 0.83 for the naive CUSUM (Student t4), where many "detections" are jump-driven.

### Q2a. Median delay in steps (FA before the change in brackets)

| Innovations | Change | sign | product c3 | product c2 | gaussian LR | naive CUSUM |
|---|---|---|---|---|---|---|
| Gaussian | 0.3 to 0.6 | 174 (4.0%) | 246 (0.1%) | 158 (0.8%) | 74 (7.1%) | 125 (5.2%) |
| Gaussian | 0.6 to 0.3 | 187 (3.1%) | 243 (0.0%) | 154 (1.2%) | 52 (7.5%) | 119 (12.9%) |
| Gaussian | 0.0 to 0.3 | 241 (3.0%) | 238 (0.0%) | 156 (0.8%) | 106 (7.0%) | 126 (2.1%) |
| Gaussian | 0.3 to 0.0 | 216 (3.4%) | 234 (0.1%) | 154 (1.4%) | 97 (7.3%) | 123 (4.5%) |
| Gaussian | 0.5 to 0.9 | 65 (4.3%) | 194 (0.0%) | 121 (1.4%) | 20 (7.1%) | 78 (9.5%) |
| Gaussian | 0.6 to 0.0 | 66 (3.3%) | 116 (0.2%) | 70 (1.3%) | 23 (8.1%) | 48 (11.4%) |
| Student t4 | 0.3 to 0.6 | 188 (4.0%) | 240 (1.0%) | 171 (2.3%) | 83 (86.9%) | 47 (97.9%) |
| Student t4 | 0.6 to 0.3 | 190 (3.2%) | 245 (0.8%) | 180 (1.8%) | 29 (96.6%) | 28 (99.7%) |
| Student t4 | 0.0 to 0.3 | 238 (4.1%) | 235 (0.3%) | 174 (2.0%) | 72 (79.8%) | 40 (96.5%) |
| Student t4 | 0.3 to 0.0 | 238 (3.3%) | 240 (0.3%) | 180 (1.5%) | 49 (85.6%) | 46 (97.6%) |
| Student t4 | 0.5 to 0.9 | 63 (2.7%) | 173 (0.2%) | 115 (1.4%) | 21 (93.0%) | 18 (99.8%) |
| Student t4 | 0.6 to 0.0 | 61 (3.8%) | 104 (0.1%) | 70 (1.2%) | 18 (95.2%) | 25 (99.9%) |
| Jumps | 0.3 to 0.6 | 191 (2.1%) | 253 (0.1%) | 163 (0.5%) | 71 (69.9%) | 112 (27.3%) |
| Jumps | 0.6 to 0.3 | 192 (4.1%) | 242 (0.1%) | 158 (1.2%) | 46 (91.7%) | 110 (37.0%) |
| Jumps | 0.0 to 0.3 | 246 (3.1%) | 253 (0.0%) | 166 (0.5%) | 134 (7.8%) | 128 (24.2%) |
| Jumps | 0.3 to 0.0 | 224 (4.0%) | 241 (0.0%) | 157 (1.3%) | 78 (68.4%) | 123 (28.0%) |
| Jumps | 0.5 to 0.9 | 65 (3.5%) | 197 (0.0%) | 122 (1.2%) | 20 (88.8%) | 78 (32.8%) |
| Jumps | 0.6 to 0.0 | 66 (2.5%) | 117 (0.2%) | 70 (1.1%) | 20 (92.7%) | 46 (37.8%) |

90th percentile delays for the sign monitor: 368 to 487 steps for the moderate changes and 110 to 124 for the two large ones (Gaussian and jumps); 376 to 505 and 109 to 123 under t4.

**Reading.**

1. On clean Gaussian pairs the Gaussian LR is the fastest, as it must be (it is the likelihood ratio of the true model): 74 steps for 0.3 to 0.6 against 174 for the sign monitor, a factor 2.4. For large changes the factor is about 3 (20 against 65). This is the price of a model-free statistic and it is in line with the classical Pitman efficiency of sign-type tests.
2. Under Student t4 or jumps the Gaussian LR and the naive CUSUM have already raised a false alarm on 70% to 100% of paths before the change, so their delay columns are not meaningful there. The sign monitor's delays are unchanged (188 against 174, 191 against 174): innovations do not matter to it.
3. Among the valid robust detectors the sign monitor is slightly faster than the product at `c = 3` for moderate changes (174 against 246 for 0.3 to 0.6) and 2 to 3 times faster for large ones (65 against 194 for 0.5 to 0.9), as Q0 predicted. Product c2 is the fastest robust option for moderate changes (158) but inherits the full volatility sensitivity of Q1. The observed sign delays agree with the Q0 "SR delay" column within 25% (174 against 219, 65 against 57, 66 against 60), and the product delays within 10% (246 against 259, 194 against 197).

## Q3. Real data [REAL]

### Q3a. The Epps effect (`q3a_epps.csv`)

Full sample, 2024 to 2025. "Sign rho" is `sin(pi m / 2)` for the mean `m` of the sign statistic. The share of exactly zero returns at 1 minute is 3.1% (BTC), 1.0% (ETH), 3.4% (SOL) and 3.6% (DOGE); at 5 minutes 0.3%, 0.2%, 1.4% and 1.5%. The median return at 5 minutes is below 1e-5 in absolute value for every coin, so centring the sign at zero is harmless.

| Pair | 1 min: Pearson / sign rho | 5 min | 15 min | 30 min | 60 min |
|---|---|---|---|---|---|
| BTC/ETH | 0.778 / 0.792 | 0.799 / 0.820 | 0.799 / 0.832 | 0.801 / 0.828 | 0.799 / 0.812 |
| BTC/SOL | 0.713 / 0.720 | 0.731 / 0.753 | 0.745 / 0.773 | 0.750 / 0.760 | 0.751 / 0.758 |
| BTC/DOGE | 0.660 / 0.730 | 0.683 / 0.768 | 0.710 / 0.785 | 0.720 / 0.783 | 0.722 / 0.774 |

By year (Pearson at 5 minutes): BTC/ETH 0.824 (2024) and 0.789 (2025); BTC/SOL 0.732 and 0.731; BTC/DOGE 0.673 and 0.698.

**Reading.** The Epps effect is visible but modest for these very liquid pairs: the 1-minute Pearson correlation is 0.02 (ETH) to 0.06 (DOGE) below its 30-minute value, and most of the gap is closed at 5 minutes (ETH) or 15 minutes (SOL, DOGE). The sign-implied correlation is above the Pearson correlation at every horizon, which is what tail dependence (co-movement in large moves) does to a Gaussian-calibrated sign measure; this does not matter for the monitor, which is calibrated on the sign mean itself. We monitor 5-minute bars: 288 per day, `A = 5760` (5% per day), calibration 1440 bars (5 days).

### Q3b. The streaming monitor on 5-minute bars (`q3b_summary.csv`, `q3b_events.csv`, `fig2_real_rho_alarms.png`)

731 days per pair. Returns were deseasonalised with the 20-day rolling intraday profile (irrelevant for the sign monitor, see README 2.4).

**Alarm counts (up / down).**

| Monitor | Band | BTC/ETH | BTC/SOL | BTC/DOGE | Per day |
|---|---|---|---|---|---|
| sign | 0 | 81 (35 / 46) | 82 (41 / 41) | 84 (38 / 46) | 0.11 |
| sign | 0.1 | 23 (10 / 13) | 24 (11 / 13) | 24 (10 / 14) | 0.03 |
| sign | 0.2 | 4 (2 / 2) | 1 (0 / 1) | 3 (1 / 2) | 0.004 |
| product c3 | 0 | 115 (60 / 55) | 114 (51 / 63) | 115 (50 / 65) | 0.16 |
| product c3 | 0.1 | 109 (59 / 50) | 106 (40 / 66) | 108 (43 / 65) | 0.15 |
| product c3 | 0.2 | 105 (42 / 63) | 102 (37 / 65) | 100 (42 / 58) | 0.14 |
| gaussian LR | 0 | 138 (27 / 111) | 136 (30 / 106) | 136 (30 / 106) | 0.19 |
| gaussian LR | 0.1 | 112 (13 / 99) | 116 (23 / 93) | 113 (18 / 95) | 0.15 |
| sign, **1-minute bars**, A = 28800 | 0.1 | 57 (14 / 43) | | | 0.08 |

Median monitored steps before an alarm (band 0.1): sign up 1900 to 3150 bars (7 to 11 days), sign down 5000 to 7600 bars (18 to 26 days); product 370 to 520 bars; Gaussian LR 180 to 470 bars. The sign monitor's up alarms start from a low calibration correlation (median implied rho 0.63 to 0.74) and its down alarms from a high one (0.77 to 0.84).

**What happens after an alarm.** For every alarm: the sign mean and the Pearson correlation over the 5 days after the alarm, against the same quantities in that segment's own calibration window ("moved as flagged" = higher after an up alarm, lower after a down alarm). Control: 2000 random times per pair, 5 days after against 5 days before.

| Monitor (band) | Pair | Up: n, share moved as flagged (sign / Pearson), median sign change | Down: n, share moved as flagged (sign / Pearson), median sign change |
|---|---|---|---|
| sign (0.1) | BTC/ETH | 10, **100% / 100%**, +0.11 | 12, 75% / 83%, -0.06 |
| sign (0.1) | BTC/SOL | 11, **91% / 100%**, +0.08 | 13, 85% / 69%, -0.04 |
| sign (0.1) | BTC/DOGE | 10, **100% / 100%**, +0.10 | 13, 92% / 77%, -0.04 |
| sign (0) | BTC/ETH | 35, 71% / 66%, +0.03 | 45, 73% / 67%, -0.02 |
| sign (0) | BTC/SOL | 41, 83% / 78%, +0.04 | 40, 93% / 80%, -0.06 |
| sign (0) | BTC/DOGE | 38, 68% / 76%, +0.05 | 45, 73% / 69%, -0.04 |
| product c3 (0.1) | BTC/ETH | 59, 61% / 68%, +0.02 | 50, 62% / 58%, -0.02 |
| product c3 (0.1) | BTC/SOL | 40, 75% / 63%, +0.02 | 66, 62% / 56%, -0.02 |
| product c3 (0.1) | BTC/DOGE | 43, 63% / 58%, +0.02 | 64, 59% / 56%, -0.01 |
| gaussian LR (0.1) | BTC/ETH | 13, 77% / 69%, +0.05 | 99, 56% / 52%, -0.01 |
| gaussian LR (0.1) | BTC/SOL | 23, 78% / 70%, +0.04 | 92, 58% / 59%, -0.02 |
| gaussian LR (0.1) | BTC/DOGE | 18, 61% / 44%, +0.01 | 94, 53% / 53%, -0.01 |
| control (random times) | BTC/ETH, SOL, DOGE | 49%, 48%, 47% higher after; median change 0.00, 0.00, 0.00 | |

**Do the alarms coincide with volatility changes?** BTC robust scale over the 5 days after the alarm divided by the 5 days before (control: median 0.99, share above 1 equal to 49%).

| Monitor (band 0.1) | Up alarms: median ratio (share above 1), ETH / SOL / DOGE pair | Down alarms |
|---|---|---|
| sign | 0.93 (30%) / 1.08 (55%) / 1.06 (70%) | 1.10 (67%) / 1.14 (62%) / 1.13 (69%) |
| product c3 | 1.05 (61%) / 1.09 (65%) / 1.08 (60%) | **0.83 (14%) / 0.83 (27%) / 0.81 (30%)** |
| gaussian LR | **0.85 (15%) / 0.86 (35%) / 0.83 (22%)** | 1.10 (60%) / 1.05 (61%) / 1.06 (57%) |

**Do dependence-up alarms cluster in sell-offs?** BTC log return over the 24 hours before the alarm (control at random times: mean +0.1%, 46% to 50% negative, 10th percentile -2.7% to -2.8% by construction).

| Monitor (band 0.1) | Up alarms: mean 24h BTC return before, share negative, share below the 10th percentile (ETH / SOL / DOGE pair) | Down alarms: same |
|---|---|---|
| sign | **-2.2% / -1.5% / -2.4%; 80% / 82% / 80%; 50% / 36% / 50%** | +1.4% / +1.8% / +0.5%; 33% / 15% / 38%; 0% / 0% / 0% |
| product c3 | -1.1% / -1.4% / -0.6%; 68% / 70% / 60%; 36% / 33% / 35% | +0.5% / +0.8% / +0.3%; 42% / 33% / 42%; 6% / 2% / 6% |
| gaussian LR | -0.1% / -0.2% / 0.0%; 38% / 48% / 44%; 0% / 4% / 6% | 0.0% / -0.1% / -0.1%; 50% / 57% / 50%; 13% / 20% / 14% |

The mean BTC return over the 24 hours *after* a sign up alarm is +0.2%, +0.5% and -0.3% (control +0.1%): the alarm says nothing about the direction of BTC afterwards, which is as it should be.

**Alarms shared across pairs** (`q3b_common_up_days.csv`). Of the 31 up alarms of the sign monitor (band 0.1), 12 are matched by an up alarm in both other pairs within 2 days; they fall in four clusters: 5 to 6 March 2024 (BTC daily log return -6.9% on 5 March), 14 to 15 March 2024 (-2.3% on 14 March), 13 June 2024 (-2.2%) and 25 to 27 November 2024 (-5.1% on 25 November). Unmatched up alarms also include 5 August 2024 (BTC/SOL; BTC -7.4% that day) and 2 April 2024 (BTC/DOGE; -6.2%). These dates are listed as computed from the data; we do not attach news events to them here.

**1-minute against 5-minute bars.** On 1-minute bars the same sign monitor (band 0.1, `A = 28800`) raises 57 alarms on BTC/ETH, 43 of them "down", against 23 (13 down) on 5-minute bars. The extra down alarms are the Epps effect at work: at 1 minute the concordance falls whenever trading becomes more asynchronous, which looks like a dependence decrease to the monitor. Aggregating to 5 minutes removes most of this.

**Reading.**

1. With band 0.1 the sign monitor raises about one alarm a month per pair, and its increase alarms are followed by a persistently higher dependence in 31 of 31 cases (Pearson yardstick) or 30 of 31 (sign yardstick), with a median rise of about +0.1 on the sign scale (about +0.1 in correlation at these levels), against 47% to 49% at random times. Decrease alarms are confirmed in 75% to 92% of cases. The product and Gaussian monitors raise 4 to 5 times more alarms and far weaker ones (53% to 78%).
2. The volatility-invariance is visible on real data. Product c3 "down" alarms coincide with a drop in BTC volatility (median ratio 0.83, only 14% to 30% with higher volatility after), exactly the confusion predicted in Q1. The Gaussian LR shows the inverted pattern also seen in Q1 (its "up" alarms coincide with volatility decreases, ratio 0.83 to 0.86; 78% to 88% of its alarms are "down" and barely better than chance). The sign monitor's alarms have volatility ratios on both sides of 1 and are not tied to a volatility direction.
3. Dependence-up alarms cluster after BTC sell-offs: 80% of them follow a negative 24-hour BTC return, and 36% to 50% follow a 24-hour return in the worst decile (10% by construction). Dependence-down alarms follow positive 24-hour returns. This is the Longin and Solnik (2001) asymmetry at the intraday scale, detected online with a false alarm guarantee, and it is the practical case for the monitor: diversification between these coins weakens exactly after the market falls.
4. The naive 5-day comparison can partly reflect mean reversion of the rolling correlation (up alarms start from low calibration correlations), as the main report's E5 also cautioned for volatility. The control at random times sets the base rate, and the Pearson yardstick (not used by the monitor) confirms the moves, but a stronger test would compare against a forecast of the correlation rather than against its calibration value.

## Q4. A basket of four assets (`q4_basket_validity.csv`, `q4_basket_delay.csv`, `q4_basket_real_*.csv`)

The basket concordance `B_t` (README 2.7) is the average of the six pairwise sign statistics, bounded in [-1/3, 1] for k = 4. One e-detector at level A on `B_t` is compared with "any of the six pairwise sign monitors, each two-sided at level 6A" (same overall guarantee by a union bound).

### Q4a. Validity [SYNTHETIC] (2000 paths, 2000-step calibration, equicorrelation 0.3)

| Scenario | basket, band 0: 1 day / 5 days | any pair at 6A, band 0 | basket, band 0.05 | any pair at 6A, band 0.05 |
|---|---|---|---|---|
| Gaussian | 1.7% / 19.3% | 2.5% / 27.1% | 0.0% / 0.2% | 0.4% / 2.5% |
| Student t4 | 2.2% / 21.9% | 3.3% / 27.5% | 0.0% / 0.5% | 0.7% / 3.0% |
| Volatility x2 in asset 1 at step 500, correlations fixed | 2.8% / 23.6% | 3.1% / 27.0% | 0.0% / 0.6% | 0.6% / 2.5% |
| Volatility x2 in all assets at step 500, correlations fixed | 2.3% / 22.9% | 2.7% / 28.2% | 0.1% / 0.6% | 0.7% / 2.6% |

The basket monitor is inside its bounds (5% and 25%) in every row and ignores the volatility changes. The union of six plug-in pairwise monitors sits at 27% to 28% per 5 days, 2 to 3 points above the 25% bound: the bound assumes the true baselines, and six estimated baselines (each from the same 2000 steps) push the union over it, as Q1c against Q1a showed for one pair. A band of 0.05 removes the excess in both designs.

### Q4b. Delay [SYNTHETIC] (1000 paths, change at step 500; median delay, 90th percentile in brackets, false alarms before the change 2.8% to 4.7% in every row)

| Change | basket | any of 6 pairs at 6A |
|---|---|---|
| All six correlations 0.3 to 0.6 | **72** (136) | 132 (235) |
| All six correlations 0.3 to 0.0 | **77** (131) | 150 (267) |
| All six correlations 0.3 to 0.9 | **21** (32) | 37 (54) |
| One pair 0.3 to 0.8, the other five fixed | 451 (1088) | **91** (148) |

**Reading.** The basket statistic averages six noisy pair statistics, so for a basket-wide change it is about twice as fast as the best pairwise monitor (72 against 132 steps) and needs no Bonferroni split. For a change confined to one pair it is five times slower (451 against 91), because the pair's shift of 0.33 on the sign scale becomes 0.33 / 6 = 0.055 in the basket mean. The two designs answer different questions: "did the basket start moving together?" (basket) and "did any pair change?" (any pair); a user who wants both can run both at level 2A each.

### Q4c. Real data [REAL]: BTC, ETH, SOL and DOGE, raw 5-minute returns, 2024 to 2025

No deseasonalisation, no scales, no truncation: the basket statistic needs none of them. Full-sample mean of `B_t` 0.573 (an average implied pairwise correlation of about 0.78). `A = 5760`, calibration 1440 bars, control at 2000 random times (47.6% higher after than before, median change -0.003; BTC 24-hour return before: mean +0.06%, 48.7% negative, 10th percentile -2.65%).

| Band | Alarms (up / down), per day | Up: share with higher basket concordance in the next 5 days than in the calibration window, median change | Up: mean BTC 24h return before, share negative, share below the 10th percentile | Down: share confirmed, median change | Down: mean BTC 24h return before, share negative |
|---|---|---|---|---|---|
| 0.05 | 72 (31 / 41), 0.10 | 87%, +0.070 | -1.1%, 77%, 13% | 88%, -0.045 | +1.6%, 15% |
| 0.1 | 31 (15 / 16), 0.04 | **93%, +0.078** | **-2.5%, 80%, 53%** | **94%, -0.060** | +1.7%, 19% |
| 0.2 | 1 (0 / 1), 0.001 | | | 100%, -0.038 | +4.7%, 0% |

Increase alarms of the band-0.1 basket monitor (UTC): 5 March 2024, 15 March 2024, 10 April 2024, 28 May 2024, 13 June 2024, 6 August 2024, 15 November 2024, 27 November 2024, 10 December 2024, 10 January 2025, 23 January 2025, 31 March 2025, 1 August 2025, 19 August 2025, 16 October 2025. Eight of the fifteen follow a 24-hour BTC return below the unconditional 10th percentile. The March, June, August and November 2024 dates coincide with the pairwise clusters of Q3b.

**Reading.** The basket version behaves like the pairwise sign monitor and is the natural "multivariate" deliverable: one guarantee for the whole basket, about one alarm a month at band 0.1, confirmed in 93% to 94% of cases against 48% at random, with increase alarms concentrated after BTC sell-offs (53% in the worst decile of 24-hour returns). The mean BTC return in the 24 hours after an increase alarm is +1.2% (control +0.06%, from `q4_basket_real_summary.csv`), with 15 alarms this is not evidence of anything and is reported only for completeness.

## Negative results and open issues

* The sign monitor is 2 to 3 times slower than the Gaussian-optimal detector on Gaussian data (Q2). This is the expected price of a sign statistic; a rank-based statistic such as a running Spearman or a winsorised product with a *predictable* scale could recover part of it while staying bounded, and was not tried.
* The product e-detector at `c = 3` is both the slowest robust option and sensitive to volatility; it has no advantage over the sign statistic in any experiment here. Its only merit is being the literal "truncated covariance" of the A7 proposal.
* With band 0, the plug-in baseline from 2000 steps leaves the 5-day false alarm rate near the bound (23% against 25%) and in the streaming test a one-standard-error baseline error produced a pre-change alarm (iteration notes). A band of 0.05 to 0.1 or the confidence-sequence baseline of the main package (`docs/math.md`, section 6) is needed in practice; the latter was not ported.
* On real data the alarms are rare with band 0.2 (1 to 4 in two years), so the band is a coarse dial on the sign scale; a finer choice (0.05, 0.15) was not explored to avoid tuning on the reported data.
* The basket statistic (Q4) covers k assets with one guarantee but dilutes a single-pair change by a factor k(k-1)/2 (five times slower than the pairwise union for one changed pair out of six). A weighted or index-based version (each asset against the equal-weight index of the others) and the k-asset analogue of the Gaussian LR benchmark were not built.
* The union of six plug-in pairwise monitors at level 6A exceeded its 5-day bound by 2 to 3 points with band 0 (27% to 28% against 25%): baseline estimation error adds up over pairs. A band of 0.05 fixes it; the confidence-sequence baseline would fix it with a proof.
* Event attribution for the alarm dates is deliberately left out (no verified news sources were consulted).

## Honest ranking of this direction

| Criterion | Score (1 to 5) | Why |
|---|---|---|
| Novelty | 4 | The e-detector framework and sign-based correlation measures are both classical, but their combination for online correlation-break detection with a finite-sample, model-free guarantee, and the explicit volatility-invariance argument for choosing the sign statistic, are not something we found in the literature. The observation that Bapat's truncation bound caps the bets of a product-based monitor (Q0) is small but new, and the bounded basket concordance (Q4) is a simple multivariate statement of the same idea. |
| Validity | 5 | The guarantee is proved in two lines from the main package's theorems and holds in every synthetic scenario, including the volatility changes and the persistent GARCH clustering that break every other detector tested (and that break the volatility monitor's own strict null, E9). |
| Practical relevance | 4 | Correlation breaks after sell-offs are the risk that diversification does not cover; the monitor finds them online, about once a month per pair, with confirmed persistence. Minus one point because the 5-day confirmation is a weak yardstick and no portfolio decision rule was tested. |
| Feasibility | 5 | About 500 lines of numpy, no scale estimate needed for the recommended statistic, 15 microseconds per step, the whole study runs in about 8 minutes. |
| Demonstration value | 4 | The "pure volatility change does not trigger it, but triggers everything else" table (Q1) and the volatility-ratio table on real data (Q3b) make the point in one glance. The sell-off clustering is a memorable real-data finding. |
| Interest to a sequential-analysis and high-frequency econometrics researcher | 4 | For sequential analysis: a clean example where the choice of a bounded, scale-free statistic buys invariance at a quantified efficiency cost, with delays matching a growth-rate calculation, and a basket-against-union comparison that makes the dilution trade-off explicit. For high-frequency econometrics: an online, jump-robust answer to the correlation-break question that the retrospective fluctuation tests address offline, plus a measured Epps effect on crypto pairs. The efficiency gap, weighted basket statistics and a forecast-relative null are the natural follow-ups. |

Overall: a strong second direction after the main volatility monitor. It reuses the same machinery, answers the paper's "multivariate" future-work item in its most useful form (dependence rather than covariation of volatility, for pairs and for a basket), and the key design choice (sign over truncated product) is backed by a proof, a theory table, simulations and real data that all point the same way.
