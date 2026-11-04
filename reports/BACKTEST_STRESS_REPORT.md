# Backtest and stress report (Agent 4): avmonitor

Scope: independent Monte Carlo of the false-alarm guarantee with fresh seeds, power and delay against the package's own baselines, one-at-a-time sensitivity, a failure-injection campaign, and an integrity check of real BTC data. The monitor is not a trading bot, so there is no Sharpe ratio: the metrics are false-alarm rate against the proved bound, detection rate and delay, and failure classes.

## Summary in plain language

1. The central claim holds under fresh seeds. For independent returns, including infinite-variance Pareto tails, jump contamination and Student t with 3 degrees of freedom, none of the proved settings ever exceeded its bound: the known-baseline e-detector reached at most 1.7% within a day (bound 5%) and 14.7% over five days (bound 25%), with 3,000 paths per law. [VERIFIED, SYNTHETIC]
2. The README setting (plug-in baseline with a 25% band) is extremely conservative: 0.0% false alarms within a day on every independent law. The price is slowness: a doubling of variance takes a median of 112 steps to detect, 2.6 times slower than the Gaussian likelihood-ratio rule, and a 1.5x increase takes 506 steps (more than a day). A fall of variance to two thirds is missed in 23% of runs within five days. [VERIFIED]
3. The plug-in baseline without the band is only approximately valid: within a day it stays under 5% (at most 3.1%) but over five days it reaches 26.7% to 28.0% against the nominal 25%. This confirms the README's own caveat. [VERIFIED]
4. When volatility clusters (GARCH) or the intraday pattern is not removed, the guarantee fails for every variant, including the README setting (7% and 47% within a day against 5%). The README says so (E10); I measured it. [VERIFIED]
5. Robustness to bad input is mixed: 93 injected failures gave 22 clean errors, 34 safe handling, 25 silent acceptances and 12 loud but uninformative crashes. The worst silent cases are a calibration length of 1 to 10, a threshold of 0, a NaN truncation constant (the monitor is blind forever with no warning) and loaders that keep an arbitrary duplicate row.

## Phase status

| Phase | What | Status |
| --- | --- | --- |
| B1 | Fresh-seed false-alarm Monte Carlo, 8 laws, 3,000 paths each, 7 methods | VERIFIED |
| B2 | Power and delay, 2 laws, 6 variance factors, 3 change times, 500 paths per cell | VERIFIED |
| B3 | One-at-a-time sensitivity, 42 configurations, 400 paths per scenario | VERIFIED |
| B4 | Failure injection, 93 cases | VERIFIED (cases run); defects found |
| B5 | Real BTC integrity, 3 months, hashes against data/MANIFEST.csv | VERIFIED for BTCUSDT 2024-01 to 2024-03; other files not downloaded |

Labels: SYNTHETIC = simulated here; REAL = Binance BTCUSDT 1-minute bars from the public mirror; OUR CHECK = computed here.

## B1. False-alarm guarantee with fresh seeds [SYNTHETIC]

Seeds: law i uses `default_rng(9_100_000 + 1000 i)` for the monitored paths and `9_150_000 + 1000 i` for the truth sample that sets the known baseline. No file in `experiments/` or `src/` uses a seed above 26,060 apart from the date-like constants 20,261,103 and 20,261,111, so these draws are new. Paths: 3,000 per law, each of length 1,000 calibration steps plus 1,950 monitored steps (five 390-step days). Threshold A = 7,800 (= 390 / 0.05), so the bound within m steps is m / A. "Within the first m steps" counts the first alarm after calibration, as in Theorem B. Laws: standard Gaussian; Student t with 3 and 5 degrees of freedom (unit variance); Gaussian plus jumps (probability 1% per step, jump sd 8); symmetric Pareto with tail index 1.5 (infinite variance, the adversarial-ish heavy tail); GARCH(1,1) with (alpha, beta) = (0.05, 0.90) and (0.08, 0.90); and iid Gaussian multiplied by an intraday U-shape that is not removed.

Methods: "Known baseline" is the package's two-sided betting e-detector with the true truncation level and true mean of the clipped square (the setting the theorem covers). The "Plug-in" rows are the streaming `VolatilityMonitor` with cal_len 1,000. The confidence-sequence row is `baseline="cs"`, delta 0.01 (bound 0.01 + m/A). Gaussian LR is `GaussianLRMonitor`; Page CUSUM and the repeated 60-step z-test come from `detectors.py` and are tuned for Gaussian data (Page threshold set for ARL = A on Gaussian simulated data, seed 9,100,099). Page and z have no guarantee; their cells are for contrast and carry no star.

Within the first 390 monitored steps (bound: 0.050 proved/nominal; confidence-sequence bound 0.060). Cell = rate [Wilson 95%]; * = interval entirely above the bound.

| Method | gaussian | student_t3 | student_t5 | jumps_p1pct_sd8 | pareto_alpha1.5 | garch_0.05_0.90 | garch_0.08_0.90 | unremoved_seasonality |
|---|---|---|---|---|---|---|---|---|
| Known baseline (theorem) | 0.008 [0.005,0.011] | 0.017 [0.013,0.022] | 0.013 [0.010,0.018] | 0.011 [0.008,0.015] | 0.016 [0.012,0.021] | 0.379 [0.361,0.396]* | 0.869 [0.856,0.881]* | 0.540 [0.522,0.558]* |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.001] | 0.000 [0.000,0.002] | 0.000 [0.000,0.001] | 0.000 [0.000,0.002] | 0.000 [0.000,0.001] | 0.073 [0.064,0.083]* | 0.469 [0.451,0.487]* | 0.024 [0.019,0.030] |
| Plug-in, band 0 | 0.018 [0.014,0.024] | 0.025 [0.020,0.031] | 0.027 [0.022,0.034] | 0.020 [0.016,0.026] | 0.031 [0.025,0.037] | 0.430 [0.412,0.447]* | 0.884 [0.872,0.895]* | 0.569 [0.551,0.587]* |
| Confidence-sequence baseline | 0.000 [0.000,0.001] | 0.000 [0.000,0.002] | 0.000 [0.000,0.001] | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] | 0.089 [0.079,0.099]* | 0.553 [0.535,0.571]* | 0.032 [0.026,0.039] |
| Gaussian LR (baseline) | 0.064 [0.056,0.073]* | 1.000 [0.999,1.000]* | 0.933 [0.924,0.942]* | 0.887 [0.876,0.898]* | 1.000 [0.999,1.000]* | 0.589 [0.572,0.607]* | 0.944 [0.935,0.952]* | 0.837 [0.823,0.849]* |
| Page CUSUM up (baseline) | 0.084 [0.074,0.094] | 1.000 [0.999,1.000] | 0.958 [0.951,0.965] | 0.887 [0.875,0.898] | 1.000 [0.999,1.000] | 0.482 [0.464,0.500] | 0.720 [0.704,0.736] | 0.746 [0.730,0.761] |
| Repeated z-test (baseline) | 0.770 [0.754,0.784] | 1.000 [0.999,1.000] | 0.999 [0.997,0.999] | 0.976 [0.970,0.981] | 1.000 [0.999,1.000] | 0.979 [0.974,0.984] | 0.999 [0.997,0.999] | 1.000 [0.998,1.000] |

Within the first 1950 monitored steps (bound: 0.250 proved/nominal; confidence-sequence bound 0.260). Cell = rate [Wilson 95%]; * = interval entirely above the bound.

| Method | gaussian | student_t3 | student_t5 | jumps_p1pct_sd8 | pareto_alpha1.5 | garch_0.05_0.90 | garch_0.08_0.90 | unremoved_seasonality |
|---|---|---|---|---|---|---|---|---|
| Known baseline (theorem) | 0.116 [0.105,0.128] | 0.139 [0.127,0.152] | 0.134 [0.122,0.147] | 0.129 [0.118,0.142] | 0.147 [0.135,0.160] | 0.972 [0.965,0.977]* | 1.000 [0.999,1.000]* | 0.993 [0.989,0.995]* |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.001] | 0.001 [0.000,0.002] | 0.001 [0.000,0.003] | 0.000 [0.000,0.002] | 0.000 [0.000,0.001] | 0.339 [0.322,0.356]* | 0.969 [0.963,0.975]* | 0.136 [0.124,0.148] |
| Plug-in, band 0 | 0.267 [0.251,0.283]* | 0.272 [0.257,0.289]* | 0.274 [0.258,0.290]* | 0.280 [0.265,0.297]* | 0.273 [0.258,0.290]* | 0.980 [0.974,0.984]* | 1.000 [0.999,1.000]* | 0.997 [0.994,0.998]* |
| Confidence-sequence baseline | 0.000 [0.000,0.002] | 0.002 [0.001,0.005] | 0.001 [0.001,0.003] | 0.001 [0.000,0.002] | 0.001 [0.000,0.003] | 0.425 [0.407,0.443]* | 0.988 [0.983,0.991]* | 0.188 [0.174,0.202] |
| Gaussian LR (baseline) | 0.301 [0.285,0.317]* | 1.000 [0.999,1.000]* | 1.000 [0.999,1.000]* | 1.000 [0.999,1.000]* | 1.000 [0.999,1.000]* | 0.992 [0.988,0.995]* | 1.000 [0.999,1.000]* | 1.000 [0.999,1.000]* |
| Page CUSUM up (baseline) | 0.300 [0.284,0.317] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] | 0.921 [0.911,0.930] | 0.987 [0.982,0.990] | 0.988 [0.983,0.991] |
| Repeated z-test (baseline) | 0.998 [0.996,0.999] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] | 1.000 [0.999,1.000] |

How to read this table:
* Independent laws (columns 2 to 6, plus Gaussian): the known-baseline e-detector, the confidence-sequence monitor and the README setting never exceeded their bound at 98, 390 or 1,950 steps. Largest rates: known baseline 1.7% within a day and 14.7% within five days (bounds 5% and 25%), so the proof is about 3 to 6 times conservative per day on these laws, in line with the README's "about four times". [VERIFIED]
* The plug-in baseline without a band is only approximate, and the rates show it: within a day it stays under 5% (largest 3.1%), but within five days it is 26.7% to 28.0% on all five independent laws against the nominal 25% (interval lower bounds 25.1% to 26.4%, so the excess is real but small). The README labels this setting "approximate"; this quantifies it. [VERIFIED]
* Gaussian likelihood ratio: 6.4% within a day even on Gaussian data (interval 5.6% to 7.3%, from estimating sigma in calibration), 88% to 100% on heavy tails and jumps. Page CUSUM: 8.4% on Gaussian, 72% to 100% elsewhere. Repeated z-test: 77% on Gaussian within a day, 98% to 100% elsewhere. [VERIFIED]
* GARCH clustering and an unremoved intraday pattern break the null as E10 says: the known-baseline rate is 38% and 87% within a day (GARCH), 54% (seasonality). The README setting also fails on GARCH (7.3% and 46.9%), the band and the confidence sequence rescue only the milder seasonality case (2.4% and 3.2% within a day, 13.6% and 18.8% over five days, inside their bounds). [VERIFIED, and it is the honest limit of the guarantee]
* Root cause class for the GARCH and seasonality failures: invalid assumption (the theorem needs a conditional mean equal to the baseline), not an implementation defect.

## B2. Power and delay [SYNTHETIC]

Seeds `9_200_000 + 100 cell`; 500 paths per cell; `simulate.change_returns` (for the jump law: probability 0.5% per step, jump sd 8; the factor applies to the diffusive variance). Variance factors 1.5, 2, 3 (increases) and 2/3, 1/2, 1/3 (decreases) at monitored steps 50, 400 and 1,500 after calibration (1,000 steps); detection counts the first alarm within 1,950 steps of the change; delay = alarm step - change step + 1. Full grid: `reports/stress/b2_power_delay.csv`. The table shows the change at step 400 (the other times give the same delays; at step 1,500 the non-guaranteed rules have already false-alarmed in 20% to 27% of runs).

Matching: A = 7,800 for every e-detector and for Gaussian LR; Page CUSUM threshold set for ARL = A on Gaussian data (seed 9,200,099); "z-test tuned" has its level set to give 5% false alarms per 390 steps on Gaussian data (seed 9,200,098, level 5.6e-4). Which are valid-level controlled: the cs monitor (proved), the known-baseline e-detector (proved, not run here as a stream); the plug-in rows are approximate; Gaussian LR is valid only for exactly normal returns; Page, the tuned z-test and the naive z-test have no guarantee.

| Law | Variance factor | Method | Detected within 5 days | Median delay (IQR) | Alarm before the change |
|---|---|---|---|---|---|
| gaussian | 1.5 | E-monitor plugin band0.25 (README) | 0.96 | 506 (290 to 820) | 0.00 |
| gaussian | 1.5 | E-monitor plugin band0 | 0.97 | 124 (75 to 190) | 0.03 |
| gaussian | 1.5 | E-monitor cs baseline | 0.95 | 483 (272 to 782) | 0.00 |
| gaussian | 1.5 | Gaussian LR monitor | 0.93 | 114 (66 to 193) | 0.07 |
| gaussian | 1.5 | Page CUSUM up | 0.91 | 83 (48 to 141) | 0.09 |
| gaussian | 1.5 | z-test tuned to 5% per day | 0.90 | 86 (52 to 157) | 0.10 |
| gaussian | 1.5 | z-test level 0.05 (naive) | 0.24 | 43 (26 to 66) | 0.76 |
| gaussian | 2 | E-monitor plugin band0.25 (README) | 1.00 | 112 (77 to 153) | 0.00 |
| gaussian | 2 | E-monitor plugin band0 | 0.99 | 53 (36 to 70) | 0.01 |
| gaussian | 2 | E-monitor cs baseline | 1.00 | 111 (75 to 155) | 0.00 |
| gaussian | 2 | Gaussian LR monitor | 0.94 | 43 (26 to 61) | 0.06 |
| gaussian | 2 | Page CUSUM up | 0.92 | 33 (22 to 49) | 0.08 |
| gaussian | 2 | z-test tuned to 5% per day | 0.91 | 39 (26 to 50) | 0.09 |
| gaussian | 2 | z-test level 0.05 (naive) | 0.23 | 23 (14 to 35) | 0.77 |
| gaussian | 3 | E-monitor plugin band0.25 (README) | 1.00 | 41 (30 to 55) | 0.00 |
| gaussian | 3 | E-monitor plugin band0 | 0.98 | 25 (19 to 34) | 0.02 |
| gaussian | 3 | E-monitor cs baseline | 1.00 | 41 (29 to 55) | 0.00 |
| gaussian | 3 | Gaussian LR monitor | 0.95 | 17 (10 to 25) | 0.05 |
| gaussian | 3 | Page CUSUM up | 0.92 | 16 (9 to 23) | 0.08 |
| gaussian | 3 | z-test tuned to 5% per day | 0.93 | 20 (13 to 27) | 0.07 |
| gaussian | 3 | z-test level 0.05 (naive) | 0.24 | 12 (8 to 17) | 0.76 |
| gaussian | 0.667 | E-monitor plugin band0.25 (README) | 0.77 | 883 (600 to 1217) | 0.00 |
| gaussian | 0.667 | E-monitor plugin band0 | 0.98 | 176 (137 to 229) | 0.02 |
| gaussian | 0.667 | E-monitor cs baseline | 0.98 | 459 (319 to 656) | 0.00 |
| gaussian | 0.667 | Gaussian LR monitor | 0.94 | 138 (90 to 211) | 0.06 |
| gaussian | 0.667 | z-test tuned to 5% per day | 0.40 | 915 (361 to 1375) | 0.10 |
| gaussian | 0.667 | z-test level 0.05 (naive) | 0.24 | 52 (38 to 75) | 0.76 |
| gaussian | 0.5 | E-monitor plugin band0.25 (README) | 1.00 | 267 (230 to 309) | 0.00 |
| gaussian | 0.5 | E-monitor plugin band0 | 0.98 | 114 (93 to 131) | 0.02 |
| gaussian | 0.5 | E-monitor cs baseline | 1.00 | 196 (165 to 231) | 0.00 |
| gaussian | 0.5 | Gaussian LR monitor | 0.92 | 66 (50 to 84) | 0.08 |
| gaussian | 0.5 | z-test tuned to 5% per day | 0.91 | 178 (96 to 362) | 0.09 |
| gaussian | 0.5 | z-test level 0.05 (naive) | 0.23 | 37 (26 to 48) | 0.77 |
| gaussian | 0.333 | E-monitor plugin band0.25 (README) | 1.00 | 154 (141 to 165) | 0.00 |
| gaussian | 0.333 | E-monitor plugin band0 | 0.97 | 82 (70 to 93) | 0.03 |
| gaussian | 0.333 | E-monitor cs baseline | 1.00 | 126 (111 to 143) | 0.00 |
| gaussian | 0.333 | Gaussian LR monitor | 0.93 | 35 (26 to 44) | 0.07 |
| gaussian | 0.333 | z-test tuned to 5% per day | 0.91 | 56 (50 to 60) | 0.09 |
| gaussian | 0.333 | z-test level 0.05 (naive) | 0.24 | 33 (24 to 40) | 0.76 |
| jumps_p0.5pct_sd8 | 2 | E-monitor plugin band0.25 (README) | 1.00 | 114 (76 to 169) | 0.00 |
| jumps_p0.5pct_sd8 | 2 | E-monitor plugin band0 | 0.98 | 51 (36 to 71) | 0.02 |
| jumps_p0.5pct_sd8 | 2 | E-monitor cs baseline | 1.00 | 115 (75 to 170) | 0.00 |
| jumps_p0.5pct_sd8 | 2 | Gaussian LR monitor | 0.31 | 35 (18 to 54) | 0.69 |
| jumps_p0.5pct_sd8 | 2 | Page CUSUM up | 0.30 | 32 (17 to 44) | 0.70 |
| jumps_p0.5pct_sd8 | 2 | z-test tuned to 5% per day | 0.32 | 35 (22 to 48) | 0.68 |
| jumps_p0.5pct_sd8 | 2 | z-test level 0.05 (naive) | 0.05 | 20 (14 to 32) | 0.95 |
| jumps_p0.5pct_sd8 | 0.5 | E-monitor plugin band0.25 (README) | 1.00 | 306 (243 to 414) | 0.00 |
| jumps_p0.5pct_sd8 | 0.5 | E-monitor plugin band0 | 0.98 | 122 (97 to 151) | 0.02 |
| jumps_p0.5pct_sd8 | 0.5 | E-monitor cs baseline | 1.00 | 222 (170 to 288) | 0.00 |
| jumps_p0.5pct_sd8 | 0.5 | Gaussian LR monitor | 0.33 | 61 (40 to 86) | 0.67 |
| jumps_p0.5pct_sd8 | 0.5 | z-test tuned to 5% per day | 0.34 | 152 (70 to 274) | 0.66 |
| jumps_p0.5pct_sd8 | 0.5 | z-test level 0.05 (naive) | 0.05 | 35 (22 to 52) | 0.95 |

Reading:
* On clean Gaussian data the guaranteed monitors are slower than the unguaranteed ones, as expected: for a variance doubling the README setting needs a median of 112 steps against 43 (Gaussian LR), 33 (Page CUSUM) and 39 (tuned z-test). The band-0 plug-in monitor is much closer (53 steps) but is only approximate.
* The advantage of the e-detector is not speed but level control under realistic returns. With jumps, Gaussian LR, Page and the tuned z-test raise a false alarm before the change in 66% to 70% of runs, so their detection rates collapse to 30% to 34%; the guaranteed monitors keep 100% detection and 0% early alarms.
* Small changes are hard: the README setting detects a 1.5x variance increase in 96% of runs but only after a median of 506 steps (IQR 290 to 820), more than a trading day; a fall to 2/3 is detected in only 77% of runs within five days (median 883). The confidence-sequence monitor is faster for decreases (459, 98%).
* Decreases are slower than increases at the same ratio for all guaranteed monitors (variance x0.5: 267 steps against x2: 112 for the README setting), consistent with the bounded-statistic argument in the README.
* Status: VERIFIED. Comparisons are at the nominal threshold, not at an equal measured false-alarm rate: the README setting is far below its allowed rate, so part of its slowness is over-conservatism that a smaller band would trade away (see B3).

## B3. Sensitivity, one parameter at a time [SYNTHETIC]

Seeds `9_300_000 + 10 scenario` (the same simulated data for every configuration, so differences are due to the setting); 400 paths per scenario. Base configuration: plug-in baseline, A = 7,800, cal_len 1,000, c = 3, null_band 0.25, default bet grids. Scenarios: null Gaussian, null Student t3 (false alarms within 390 and 1,950 steps), variance x2 and x0.5 at monitored step 400 (detection rate within 1,950 steps and median delay). Full numbers: `reports/stress/b3_sensitivity.csv` and `b3_summary.csv`. With 400 paths the false-alarm rates near zero cannot separate 0% from 1%; the cliffs below are far larger than that.

| Parameter | Values tried | What happens | Shape |
| --- | --- | --- | --- |
| Calibration length | 100, 200, 390, 500, 1000, 2000, 4000 | Delay of a x2 increase 101 to 117 steps for every length; false alarms stay at or below 0.8% (Gaussian, 390 steps) and 6.8% over five days (t3, 100 steps, still under the 25% bound); detection dips to 97.8% to 98.8% at 100 | gradual |
| Truncation constant c | 1, 1.5, 2, 3, 4, 6, 10, no truncation | c = 1: only 7% of x2 increases and 47% of x0.5 decreases detected (delay 964 and 841). c = 1.5: 99.5% and 100%. c = 3 is the best compromise: increase delay 108, decrease delay 258. Larger c speeds up increases a little (97.5 at c >= 6) but slows decreases badly (406 at c = 4, 770 at 6, 1,698 at 10 with 84% detected) and with no truncation the decrease side detects nothing (0%) | CLIFF below c = 1.5 (power), GRADUAL above 3 (decreases); no truncation is catastrophic for decreases |
| MAD scaling constant | 0.5, 1.0 (raw MAD as printed in the paper), 1.4826, 2.0 | Equivalent to changing c in sd units (0.5 gives c of about 1.0 sd, 1.0 gives 2.0 sd, 2.0 gives 4.0 sd). 0.5: the same collapse as c = 1 (9% and 51%); 1.0: works, increase delay 174, decrease 191; 2.0: increase 100, decrease 410 | CLIFF at 0.5, otherwise gradual |
| Null band | 0, 0.1, 0.25, 0.5 | 0: fastest (50 and 112 steps) but false alarms over five days 25.2% (Gaussian) and 26.0% (t3) against the 25% nominal; 0.1: 71 and 152 steps, 1.2% over five days; 0.25: 109 and 258 steps, 0%; 0.5: increases still detected (289 steps) but only 12% of x0.5 decreases (a halving is no longer outside the band) | gradual in delay; CLIFF for decreases once the band approaches the effect size |
| Delta (confidence-sequence baseline, band 0) | 0.001, 0.01, 0.05, 0.2 | Delay 134 to 76 (increase) and 217 to 154 (decrease); false alarms at most 1.2% over five days | gradual |
| Increase bets | default 6, single 0.05, single 0.5, single 0.9, wide 9 | Delay 104 to 114 for the grid or single 0.5; 212 for single 0.05; 258 for single 0.9 (with 98.8% detection and 1.5% false alarms over five days under t3) | gradual |
| Decrease fractions | default 5, single 0.1, 0.5, 0.95 | Decrease delay 258 (default), 1,122 (single 0.1), 331 (0.5), 210 (0.95); increases unaffected | gradual |
| Threshold A | 1,950 to 78,000 | Delay grows smoothly with log A: increase 83 to 153, decrease 196 to 357; no false alarms in 400 null paths at any A | gradual |

Reading: the guarantee side (false alarms) is insensitive to everything except dropping the band, where the plug-in sits on its approximate bound. The power side has two cliffs, both avoidable: truncating too hard (c at or below 1, or a MAD scaling of 0.5) and a band as wide as the effect to be found. Nothing here breaks the null in this iid setting. Status: VERIFIED.


## B4. Failure injection (93 cases, SYNTHETIC inputs)

Seeds 9,500,000 + case. Classes: DETECTED 22 (the package raises its own error), SAFELY HANDLED 34, SILENTLY IGNORED 25, CRASHED 12 (loud failure from library internals or delayed far from the bad input). Full table with observed behaviour: `reports/stress/b4_failure_cases.csv`.

| Input | Behaviour | Class |
| --- | --- | --- |
| NaN (one point, a block of 300, 50 inside calibration), inf | skipped as missing, state unchanged, later alarms identical to the clean run apart from the index shift | SAFELY HANDLED |
| One return of 1e6 sd, a +50 sd then -50 sd bounce, 1e300 | truncation bounds the influence; no alarm. Realistic fat-finger (a 5% bar that reverts on 0.05% minute returns, README setting, 300 paths): 0 alarms | SAFELY HANDLED |
| Empty input, one point | no-op | SAFELY HANDLED |
| `calibrate()` on empty, one point, all zeros | ValueError with a clear message | DETECTED |
| Bad parameters: cal_len < 0, c <= 0, null_band < 0, delta outside (0,1) for cs, tau_frac outside, bet grids with 1.0 or negatives, c_down <= 0, unknown baseline | ValueError (for cal_len < 0 and c <= 0 the message names a different parameter, cal_max and c_down) | DETECTED |
| All-NaN stream, all-zero stream (constant price) | no alarm ever, calibration restarts forever, no warning; only `cal_stats["restarted"]` shows it | SILENTLY IGNORED |
| 60% zero returns | fallback scale path used, counted in `cal_stats` | SAFELY HANDLED |
| cal_len 0, 1, 2, 5, 10 | accepted; 12, 12, 12, 1 and 8 alarms in 3,000 null steps (a valid monitor should give about 0.4) | SILENTLY IGNORED |
| threshold 0, negative, 1, NaN | accepted; thresholds <= 1 alarm at every first monitored step, NaN never alarms | SILENTLY IGNORED |
| c = NaN | accepted; calibration never completes, blind forever | SILENTLY IGNORED |
| null_band >= 1 | decrease side switched off silently | SILENTLY IGNORED |
| Outliers inside the calibration window (README setting, cal_len 500, 300 paths per row) | 0 outliers 0.0% alarms, 10 outliers of 15 sd 0.7%, 25 outliers 68.7% within 2,500 steps | invalid assumption, SILENTLY IGNORED beyond about 5% contamination |
| Strings, None, 2-D arrays passed to `update`/`run` | accepted during calibration, TypeError or ValueError only when calibration ends (or at once while monitoring) | CRASHED (delayed) |
| Interrupt and resume | pickle and deepcopy resume gives identical alarms. There is no save/load API; a cold restart on the remaining data recalibrates and cannot see a shift that already happened | SAFELY HANDLED |
| Loader: unsorted rows | sorted | SAFELY HANDLED |
| Loader: duplicate timestamps | dropped silently, count not reported by the Bitstamp loader; which duplicate survives is arbitrary because `sort_index()` is not stable (on 286 duplicated minutes about half kept the later, absurd row, whether the file was shuffled or not) | SILENTLY IGNORED (defect) |
| Loader: zero or negative close (Bitstamp and Nifty loaders) | +-inf or NaN returns with a RuntimeWarning only (the Binance loader masks them and reports the count) | SILENTLY IGNORED (defect) |
| Loader: three missing rows in the Bitstamp file | a 4-minute return treated as a 1-minute return (the source normally fills gaps) | SILENTLY IGNORED |
| Loader: timestamps in milliseconds read as seconds | accepted, dates around year 55,000 | SILENTLY IGNORED |
| Loader: empty file, header only, missing column, text in a price, space before a quote, empty Nifty body | pandas/numpy errors (EmptyDataError, TypeError, KeyError, IndexError) with no message about the file | CRASHED (loud) |
| Loader: truncated last row, BOM header, flat bars | NaN or handled | SAFELY HANDLED |

Root causes of the 40 non-benign cases (silent or crashed 37, plus 3 accidental messages): 28 implementation gaps (no input validation; 21 silent, 7 crashes), 5 errors raised by pandas or numpy, 2 invalid assumptions (contaminated calibration), 1 data problem (duplicates), 3 accidental error messages. None of the failures is a wrong statistic on valid input.

## B5. Real-data integrity [REAL]

Downloaded with the package's own `fetch` into a scratch folder: Binance BTCUSDT 1-minute bars for 2024-01, 2024-02 and 2024-03 (3 files, 20 MB). Script `reports/stress/b5_data_integrity.py`, output `b5_btc_integrity.csv`.

| Check | 2024-01 | 2024-02 | 2024-03 |
| --- | --- | --- | --- |
| SHA-256 equals data/MANIFEST.csv | yes | yes | yes |
| Row count equals manifest | 44,640 | 41,760 | 44,640 |
| Duplicate timestamps | 0 | 0 | 0 |
| Strictly increasing timestamps, whole minutes | yes | yes | yes |
| Missing minutes | 0 | 0 | 0 |
| Non-positive or NaN price | 0 | 0 | 0 |
| OHLC inconsistent bars | 0 | 0 | 0 |
| Zero-volume bars, flat bars | 0, 0 | 0, 0 | 0, 0 |
| Zero returns share | 1.9% | 3.1% | 1.3% |
| Returns beyond 20 robust sd | 20 | 10 | 6 |
| Largest return / robust sd | 57.8 | 70.5 | 29.6 |

Reading: the three files are exactly the committed ones and clean; the 36 returns beyond 20 robust sd are genuine heavy tails of BTC, not data errors (no bar has high below low, and the largest high-low range is 4.4% of the close). Other months, other coins and the Nifty and Bitstamp files were not downloaded, so the manifest comparison covers 3 of 140 files. Status: VERIFIED for these three files, UNVERIFIED for the rest.

## Defects and weaknesses found (for the coordinator, owner Agent 1 for src/)

| Id | Where | What | Severity |
| --- | --- | --- | --- |
| A4-B-001 | `data.load_bitstamp_csv` | `sort_index()` is not stable, so the surviving duplicate is arbitrary; use `kind="mergesort"` as `_finish` does, and report the duplicate count | Medium |
| A4-B-002 | `data.load_bitstamp_csv`, `data.load_nifty` | non-positive prices give +-inf / NaN returns with a warning only; mask them and report them as `load_binance` does | Medium |
| A4-B-003 | `monitor.VolatilityMonitor.__init__` | no validation of cal_len (>= 2 or so), threshold (> 1), c (finite, > 0, message names c), null_band (< 1) | Medium |
| A4-B-004 | `monitor.VolatilityMonitor.update` | non-numeric input accepted during calibration, fails later; coerce with `float(r)` | Low |
| A4-B-005 | `monitor.VolatilityMonitor` | no warning when the monitor is blind (all NaN, constant price); expose a `status()` or emit a warning when calibration restarts repeatedly | Medium |
| A4-B-006 | docs | the README setting is 2.6 times slower than a Gaussian LR rule and misses a fall to 2/3 of the variance in 23% of runs; say so next to the false-alarm claim | Low |

## Commands and seeds

```
cd reports/stress    (all with OMP_NUM_THREADS=1, python -I)
python b1_false_alarm.py 3000 b1_false_alarm.csv           # seeds 9_100_000 + 1000 i, 9_150_000 + 1000 i, Page 9_100_099
python b2_power_delay.py 500 b2_power_delay.csv            # seeds 9_200_000 + 100 cell, 9_200_098, 9_200_099
python b3_sensitivity.py 400 b3_sensitivity.csv            # seeds 9_300_000 + 10 scenario
python b4_failure_cases.py b4_raw_cases.csv && python b4_classify.py   # seeds 9_500_000 and up
python b5_data_integrity.py <scratch_dir> b5_btc_integrity.csv 2024-01 2024-02 2024-03
```
Seeds from 9,000,000 upward are not used in `experiments/` or `src/`. Wilson intervals use z = 1.96.

## What this report does not show

Everything in B1 to B4 is simulated; B5 covers three months of one coin. Paths are independent across rows but the same simulated path is fed to all methods in a law (a paired comparison). Detection comparisons use the nominal threshold, not an equal measured false-alarm rate. The pre-registered live session and the real-data event studies in the README were not rerun.
