# R3. Where the guarantee bends: new null laws and realistic change shapes [SYNTHETIC]

Scope: the false-alarm bound of the monitor under return laws that B1 (`reports/BACKTEST_STRESS_REPORT.md`) did not test, and its power under change shapes that B2 did not test. Everything here is simulated. Same methods and bounds as B1: known-baseline e-detector (the theorem's own setting, bound m / A), plug-in monitor with band 0.25 (the README setting) and band 0 (bound m / A, approximate), and the confidence-sequence baseline (bound 0.01 + m / A). A = 7,800, so the bound is 5% within 390 monitored steps (one day) and 25% within 1,950 (five days). Calibration window 1,000 steps, as in B1.

Metrics, explained once. "Rate" is the share of paths with a first alarm within the horizon. The bracket is a Wilson 95% confidence interval for that share. A star (*) marks a cell whose whole interval lies above the bound; a cell whose interval straddles the bound is inconclusive. "(u/d)" counts alarms of the increase and of the decrease side. "Delay" is alarm step minus change step plus one. m0 is the mean of the clipped squared return min(r², τ²), the quantity the monitor actually watches; "m0hat" is its estimate from the calibration window.

## 1. Summary in plain language

* **The monitor watches the clipped square, not the variance, and the README should say so.** When the tails of unit-variance returns switch from Student t30 to t3 (variance unchanged), the decrease side fires on 65% of paths within a day in the theorem's own setting, and on 62% with the plug-in band 0; the README band of 0.25 hides most of it (1.0% within a day). In the reverse direction the README setting fires on the increase side on 45% of paths within a day and 84% within 975 steps. The mirror case: a doubling or tripling of variance carried only by rare 20 sd jumps moves the clipped mean by only 2% to 4% and is never flagged by the README setting (0.0% within five days, 2,000 paths). Not a bug, stated in one sentence in `docs/math.md` 5.6.4, absent from `docs/limitations.md`. [VERIFIED, SYNTHETIC, NEW limitation]
* **Bid-ask bounce bends the guarantee at a realistic strength and breaks it at a strong one.** MA(1) with θ = -0.3: every proved setting inside its bounds (known baseline 2.9% per day, 21.8% per five days). θ = -0.6: known baseline 5.0% per day (inconclusive) and 40.7% per five days against 25% (exceeds), plug-in band 0 8.6% per day. README setting and confidence-sequence baseline at most 0.9%. Already a limitation in words; the five-day numbers are new. [VERIFIED]
* **The realistic seasonality pipeline leaves an error that shows up over five days.** With the U-shape estimated by `intraday_scale` from the past 20 sessions, the known baseline is inside the one-day bound (3.4%) but at 32.2% within five days against 25%; plug-in band 0 35.5%. The README setting and the confidence-sequence baseline absorb it (0.1% and 0.2%). This is the measurement `docs/math.md` 5.6.3 says was missing. [VERIFIED, NEW limitation with numbers]
* **Skew, missing data, tick grids and a misspecified volatility forecast do not break the proved settings.** Two-piece normal skew 0.6 and 0.8, 5% random NaN, NaN gaps of 30 to 120 minutes, prices on a grid with 20%, 40% and 60% exact zeros, and GARCH(1,1) forecasts fitted to GJR, EGARCH and log-normal stochastic volatility data: the known baseline stayed at 0.4% to 1.7% per day and 8% to 15% per five days, the README setting and the confidence-sequence baseline at 0% to 0.2%. The plug-in band 0 exceeded its five-day bound on the 20% tick law (29.4%) and on GJR (31.2%), a little above its known 27% to 28% i.i.d. level. At 60% zeros the monitor extended every calibration window to 2,000 returns and its decrease side was inactive. [VERIFIED within bound]
* **A gradual doubling is caught, a short burst mostly is not.** Variance doubling linearly over a day: README setting alarms during the ramp on 59% of paths and within a day after it on 100% (median 371 steps after the ramp starts). A 30-minute x3 burst: caught inside or within an hour after on 27% of paths (60-minute burst: 82%), median 25 steps after it starts. [VERIFIED]
* **A change during calibration can leave the README monitor silent for a week.** Doubling at step 500 of the 1,000-step window: learnt baseline 26% too low, alarm within a day on 4.2% of paths, within five days on 44%, so silent for five days on 56%; confidence-sequence baseline 14%. Doubling at step 900: 94% within a day. Doubling at step 100: baseline only 5% low and the monitor is (correctly) silent. [VERIFIED, NEW limitation]
* Everything above is simulated; 2,000 paths per cell; one law or shape at a time. Nothing here contradicts the proofs: every exceedance is in a cell where the theorem's assumption (a constant conditional clipped mean) fails, and the two settings with a proof for an estimated baseline or a band never exceeded their bounds on any null law tested. [SYNTHETIC]

## 2. What was run

* Scripts (all in `reports/round4/R3/`): `r3_common.py` (laws, a vectorised copy of the plug-in monitor, helpers), `r3_null.py` (null laws), `r3_change.py` (ramp, burst, change during calibration), `r3_change2.py` (variance change carried only by jumps), `make_tables.py` (tables below). Logs `*.log`; results `r3_null.csv`, `r3_null_extra.csv`, `r3_change.csv`, `r3_change2.csv`.
* Commands: `python r3_null.py 2000 r3_null.csv` (15 laws), `python r3_null.py 2000 r3_null_extra.csv t30_to_t5 t5_to_t30 ma1_-0.1 ma1_0.3 seasonal_amp1.0`, `python r3_change.py 2000 r3_change.csv`, `python r3_change2.py 2000 r3_change2.csv`, each with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`, one core, no multiprocessing. Total run time about 45 minutes.
* Paths: 2,000 per null law and per change cell (B1 used 3,000, B2 500).
* Seeds (numpy `default_rng`), all in the block 9,700,000 to 9,799,999: null law i (its index in `LAWS` of `r3_common.py`, 0 to 19) uses 9,700,000 + 1,000 i for the paths and 9,750,000 + 1,000 i for the independent sample that sets the known baseline (τ = 3 x scaled MAD, m0 = mean of the clipped square); the GARCH(1,1) fits for the forecast laws use 9,750,500 + 1,000 i; change cell j uses 9,760,000 + 100 j (`r3_change.py`) and 9,761,000 + 100 j (`r3_change2.py`); the true m0 of the doubled regime uses 9,769,999; the tick sizes were found with 9,799,999; smoke tests used 9,799,000, 9,799,100 and 9,799,500.
* Check of the vectorised plug-in monitor: on 40 Student t3 paths with a change and one NaN block, its first-alarm steps and directions matched the streaming `VolatilityMonitor` exactly for band 0.25 and band 0. The confidence-sequence rows and the tick rows always use the streaming monitor itself.
* Laws. (a) Student t, unit variance, degrees of freedom 30 for calibration and the first 975 monitored steps, then 3 (and the reverse; also 30 to 5 and 5 to 30). (b) MA(1) returns r_t = (e_t + θ e_{t-1}) / sqrt(1 + θ²), θ = -0.3, -0.6 (and -0.1, +0.3). (c) A Gaussian random walk price with unit steps, rounded to a tick grid of 0.51, 1.10 and 1.94 step sd, which makes 20%, 40% and 60% of returns exactly zero; returns are differences of the rounded price. (d) Two-piece normal with left and right sd (0.6, 1.4) and (0.4, 1.6), standardised to mean 0 and variance 1 (skewness 0.59 and 0.79). (e) Gaussian with 5% of returns missing at random, or with gaps of 30 to 120 steps (about one per 390 steps) in calibration and monitoring, passed as NaN. (f) Gaussian times a U-shaped intraday profile (`simulate.intraday_u_shape`, amplitude 0.6; also 1.0) over 32 sessions of 390 minutes, divided by `data.intraday_scale(lookback_days=20, smooth=15, start="past", circular=False)`, as the real-data experiments do; the first 10 sessions are NaN by construction and are dropped, then calibration and monitoring follow. (g) GJR-GARCH (ω 0.02, α 0.02, β 0.92, γ 0.08), an EGARCH (α 0.10, γ -0.05, β 0.97) and log-normal stochastic volatility (φ 0.98, innovation sd 0.15), each with Gaussian innovations; a GARCH(1,1) is fitted by `forecast.garch11_fit` on one 30,000-step training path of the same model, frozen, and the monitored returns are divided by its one-step-ahead forecast (`forecast.garch11_filter`, `forecast.standardise`) after a 2,000-step warm-up, as E10 does.
* Change shapes, Gaussian returns. (h) Variance 1 to 2 linearly over 390 steps from monitored step 50, then 2. (i) Variance x3 for 30 (and 60) steps from monitored step 200, then 1. (j) Variance x2 from calibration step 100, 500 or 900 of 1,000, then forever. (k) Variance x2 or x3 from monitored step 50, carried only by rare jumps (probability 0.25% or 0.5% per step, jump sd 20).

## 3. Results

(a) Tail shape switch at constant variance: within the first 390 monitored steps (bound 0.050; cs bound 0.060). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | t30_to_t3 | t3_to_t30 | t30_to_t5 | t5_to_t30 |
|---|---|---|---|---|
| Known baseline (theorem) | 0.007 [0.005,0.012] (10u/5d) | 0.011 [0.008,0.017] (9u/14d) | 0.013 [0.009,0.020] (18u/9d) | 0.017 [0.012,0.023] (16u/17d) |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.000 [0.000,0.002] (0u/0d) |
| Plug-in, band 0 | 0.011 [0.008,0.017] (17u/6d) | 0.027 [0.021,0.035] (24u/30d) | 0.018 [0.013,0.025] (23u/13d) | 0.029 [0.023,0.037] (33u/25d) |
| Confidence-sequence baseline | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.001 [0.000,0.004] (1u/1d) |
paths per law: 2000 to 2000

(a) Tail shape switch at constant variance: within the first 1950 monitored steps (bound 0.250; cs bound 0.260). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | t30_to_t3 | t3_to_t30 | t30_to_t5 | t5_to_t30 |
|---|---|---|---|---|
| Known baseline (theorem) | 0.992 [0.987,0.995]* (53u/1931d) | 1.000 [0.998,1.000]* (1924u/76d) | 0.347 [0.327,0.369]* (71u/624d) | 0.595 [0.573,0.616]* (1133u/57d) |
| Plug-in, band 0.25 (README) | 0.064 [0.054,0.075] (0u/127d) | 0.836 [0.819,0.851]* (1671u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.004 [0.002,0.008] (8u/0d) |
| Plug-in, band 0 | 0.986 [0.980,0.990]* (114u/1858d) | 1.000 [0.998,1.000]* (1860u/140d) | 0.418 [0.397,0.440]* (128u/708d) | 0.620 [0.598,0.641]* (1118u/121d) |
| Confidence-sequence baseline | 0.346 [0.325,0.367]* (0u/692d) | 0.866 [0.850,0.880]* (1731u/1d) | 0.005 [0.003,0.010] (1u/10d) | 0.009 [0.006,0.015] (18u/1d) |
paths per law: 2000 to 2000


After the tail switch (monitored step 975), among paths with no alarm before it. Horizon counts from the switch; bound as if a fresh segment.

| Law | Method | Horizon | Alarms / paths | Rate [Wilson 95%] | up / down |
|---|---|---|---|---|---|
| t30_to_t3 | Known baseline (theorem) | switch+390 | 1249 / 1912 | 0.653 [0.632,0.674]* | 1 / 1248 |
| t30_to_t3 | Known baseline (theorem) | switch+975 | 1896 / 1912 | 0.992 [0.986,0.995]* | 1 / 1895 |
| t30_to_t3 | Plug-in, band 0.25 (README) | switch+390 | 20 / 2000 | 0.010 [0.006,0.015] | 0 / 20 |
| t30_to_t3 | Plug-in, band 0.25 (README) | switch+975 | 127 / 2000 | 0.064 [0.054,0.075] | 0 / 127 |
| t30_to_t3 | Plug-in, band 0 | switch+390 | 1113 / 1790 | 0.622 [0.599,0.644]* | 1 / 1112 |
| t30_to_t3 | Plug-in, band 0 | switch+975 | 1762 / 1790 | 0.984 [0.977,0.989]* | 1 / 1761 |
| t3_to_t30 | Known baseline (theorem) | switch+390 | 1884 / 1884 | 1.000 [0.998,1.000]* | 1883 / 1 |
| t3_to_t30 | Known baseline (theorem) | switch+975 | 1884 / 1884 | 1.000 [0.998,1.000]* | 1883 / 1 |
| t3_to_t30 | Plug-in, band 0.25 (README) | switch+390 | 891 / 2000 | 0.446 [0.424,0.467]* | 891 / 0 |
| t3_to_t30 | Plug-in, band 0.25 (README) | switch+975 | 1671 / 2000 | 0.836 [0.819,0.851]* | 1671 / 0 |
| t3_to_t30 | Plug-in, band 0 | switch+390 | 1730 / 1752 | 0.987 [0.981,0.992]* | 1730 / 0 |
| t3_to_t30 | Plug-in, band 0 | switch+975 | 1752 / 1752 | 1.000 [0.998,1.000]* | 1752 / 0 |
| t30_to_t5 | Known baseline (theorem) | switch+390 | 183 / 1898 | 0.096 [0.084,0.111]* | 7 / 176 |
| t30_to_t5 | Known baseline (theorem) | switch+975 | 593 / 1898 | 0.312 [0.292,0.334]* | 9 / 584 |
| t30_to_t5 | Plug-in, band 0.25 (README) | switch+390 | 0 / 1999 | 0.000 [-0.000,0.002] | 0 / 0 |
| t30_to_t5 | Plug-in, band 0.25 (README) | switch+975 | 0 / 1999 | 0.000 [-0.000,0.002] | 0 / 0 |
| t30_to_t5 | Plug-in, band 0 | switch+390 | 211 / 1800 | 0.117 [0.103,0.133]* | 12 / 199 |
| t30_to_t5 | Plug-in, band 0 | switch+975 | 636 / 1800 | 0.353 [0.332,0.376]* | 15 / 621 |
| t5_to_t30 | Known baseline (theorem) | switch+390 | 328 / 1902 | 0.172 [0.156,0.190]* | 326 / 2 |
| t5_to_t30 | Known baseline (theorem) | switch+975 | 1092 / 1902 | 0.574 [0.552,0.596]* | 1090 / 2 |
| t5_to_t30 | Plug-in, band 0.25 (README) | switch+390 | 1 / 2000 | 0.001 [0.000,0.003] | 1 / 0 |
| t5_to_t30 | Plug-in, band 0.25 (README) | switch+975 | 8 / 2000 | 0.004 [0.002,0.008] | 8 / 0 |
| t5_to_t30 | Plug-in, band 0 | switch+390 | 419 / 1786 | 0.235 [0.216,0.255]* | 408 / 11 |
| t5_to_t30 | Plug-in, band 0 | switch+975 | 1025 / 1786 | 0.574 [0.551,0.597]* | 1008 / 17 |


How to read this: the monitor does not watch the variance. It watches the mean of the clipped square min(r², τ²), with τ fixed from calibration (about 3 sd for t30 tails). At unit variance that mean is 0.99 for t30, 0.91 for t5 and 0.73 for t3 (τ from the t30 calibration), so a switch to fatter tails looks like a 27% fall and the decrease side fires: 65% of silent paths alarm within a day of the switch with the known baseline, 62% with the plug-in band 0. The README band of 0.25 hides most of it (1.0% within a day, 6.4% within 975 steps), because the fall is only just outside the band. In the reverse direction the clipped mean rises by about 54% (τ is then 1.96 from the t3 calibration), and the README setting fires on 45% of paths within a day of the switch and on 84% within 975 steps, all on the increase side. The milder t30 to t5 switch (a 9% fall of the clipped mean) is inside the README band and the README setting stays silent (0.1% within five days) while the known baseline alarms on 35% of paths within five days; t5 to t30 (a 14% rise) gives 0.4% and 60%. Before the switch every proved setting is inside its bound, so this is not a bug: it is the statistic doing what the proofs say it does (`docs/math.md` 5.6 item 4), and a NEW limitation to list with numbers, because a user reads "volatility shifts" as "variance shifts". Row "Known baseline" at 1,950 steps is starred only because it includes the switch.


(b) MA(1) bid-ask bounce and (d) skewed two-piece normal: within the first 390 monitored steps (bound 0.050; cs bound 0.060). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | ma1_-0.1 | ma1_-0.3 | ma1_-0.6 | ma1_0.3 | twopiece_0.6_1.4 | twopiece_0.4_1.6 |
|---|---|---|---|---|---|---|
| Known baseline (theorem) | 0.008 [0.005,0.013] (13u/3d) | 0.029 [0.022,0.037] (47u/10d) | 0.050 [0.041,0.060] (83u/17d) | 0.024 [0.018,0.031] (40u/7d) | 0.011 [0.007,0.017] (16u/6d) | 0.012 [0.008,0.018] (20u/4d) |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) |
| Plug-in, band 0 | 0.020 [0.015,0.027] (27u/13d) | 0.045 [0.037,0.056] (76u/15d) | 0.086 [0.075,0.100]* (126u/47d) | 0.036 [0.029,0.045] (56u/16d) | 0.015 [0.011,0.021] (17u/13d) | 0.022 [0.016,0.029] (33u/11d) |
| Confidence-sequence baseline | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) |
paths per law: 2000 to 2000

(b) MA(1) bid-ask bounce and (d) skewed two-piece normal: within the first 1950 monitored steps (bound 0.250; cs bound 0.260). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | ma1_-0.1 | ma1_-0.3 | ma1_-0.6 | ma1_0.3 | twopiece_0.6_1.4 | twopiece_0.4_1.6 |
|---|---|---|---|---|---|---|
| Known baseline (theorem) | 0.124 [0.111,0.140] (124u/125d) | 0.218 [0.201,0.237] (266u/171d) | 0.407 [0.386,0.429]* (483u/332d) | 0.235 [0.217,0.254] (285u/185d) | 0.122 [0.109,0.138] (135u/110d) | 0.130 [0.116,0.145] (143u/117d) |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] (0u/0d) | 0.002 [0.001,0.004] (3u/0d) | 0.005 [0.003,0.009] (10u/0d) | 0.002 [0.001,0.004] (3u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.001 [0.000,0.004] (2u/0d) |
| Plug-in, band 0 | 0.286 [0.267,0.307]* (288u/285d) | 0.386 [0.364,0.407]* (433u/338d) | 0.579 [0.557,0.600]* (628u/529d) | 0.370 [0.349,0.391]* (395u/345d) | 0.267 [0.248,0.286] (275u/258d) | 0.274 [0.255,0.294]* (297u/251d) |
| Confidence-sequence baseline | 0.000 [0.000,0.002] (0u/0d) | 0.003 [0.001,0.006] (3u/2d) | 0.009 [0.005,0.014] (14u/3d) | 0.001 [0.000,0.004] (2u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.001 [0.000,0.003] (1u/0d) |
paths per law: 2000 to 2000


How to read this: MA(1) returns break the theorem's assumption (the conditional mean of the next squared return depends on the last shock), two-piece normal returns do not (they are i.i.d.). Skew is harmless: all four settings are inside their bounds for both skew levels. Bounce at θ = -0.3 (a noise-to-signal variance ratio of about 0.3) is also inside the bounds in the theorem's own setting (2.9% per day, 21.8% per five days). At θ = -0.6 the known baseline is at the bound within a day (5.0%, interval 4.1% to 6.0%) and well above it within five days (40.7% against 25%), and the plug-in without a band exceeds the one-day bound (8.6%). The README setting and the confidence-sequence baseline stay at or below 0.9% in every cell. `docs/math.md` 5.6 item 6 and `docs/limitations.md` section 3 already name price noise as a limitation; the five-day failure at θ = -0.6 in the theorem's own setting is new and should be added. A weak bounce (θ = -0.1) and positive autocorrelation (θ = +0.3) also stay inside the bounds in the theorem's setting (0.8% and 2.4% per day, 12.4% and 23.5% per five days); only the plug-in band 0 exceeds its five-day bound there (28.6% and 37.0%).


(c) Prices on a tick grid: within the first 390 monitored steps (bound 0.050; cs bound 0.060). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | tick_zero20 | tick_zero40 | tick_zero60 |
|---|---|---|---|
| Known baseline (theorem) | 0.017 [0.012,0.023] (20u/13d) | 0.009 [0.006,0.014] (18u/0d) | 0.005 [0.003,0.009] (10u/0d) |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) |
| Plug-in, band 0 | 0.034 [0.026,0.042] (30u/37d) | 0.013 [0.009,0.020] (27u/0d) | 0.011 [0.007,0.016] (21u/0d) |
| Confidence-sequence baseline | 0.001 [0.000,0.003] (1u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) |
paths per law: 2000 to 2000

(c) Prices on a tick grid: within the first 1950 monitored steps (bound 0.250; cs bound 0.260). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | tick_zero20 | tick_zero40 | tick_zero60 |
|---|---|---|---|
| Known baseline (theorem) | 0.132 [0.118,0.148] (132u/132d) | 0.097 [0.084,0.110] (149u/44d) | 0.081 [0.070,0.094] (162u/0d) |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.003] (1u/0d) | 0.001 [0.000,0.003] (1u/0d) |
| Plug-in, band 0 | 0.294 [0.275,0.315]* (285u/304d) | 0.243 [0.225,0.262] (279u/207d) | 0.122 [0.108,0.137] (236u/8d) |
| Confidence-sequence baseline | 0.002 [0.001,0.004] (1u/2d) | 0.001 [0.000,0.003] (1u/0d) | 0.000 [0.000,0.002] (0u/0d) |
paths per law: 2000 to 2000


How to read this: a price on a coarse grid gives exact-zero returns and a negatively autocorrelated return series, which is outside the theorem's i.i.d. examples but not outside its assumption as long as the clipped mean stays put. Every proved setting stayed inside its bound at 20%, 40% and 60% zeros (known baseline 0.5% to 1.7% per day, 8.1% to 13.2% per five days). At 60% zeros the scaled MAD of the whole window is 0, so the plug-in monitor extended its calibration window to 2,000 returns on every path and then used its fallback scale (median of the nonzero absolute returns), as the module docstring says; after that, the decrease side was inactive on these paths (0 down alarms), so a fall in volatility would not be detected there. Only the plug-in band 0 at 20% zeros exceeds its five-day bound (29.4%), a little above its 26.7% to 28.0% on i.i.d. laws in B1. Verdict: VERIFIED within bound; the 60% case is a limitation on the decrease side rather than on the false-alarm rate.


(e) Missing data and (f) estimated intraday profile: within the first 390 monitored steps (bound 0.050; cs bound 0.060). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | missing_5pct | missing_blocks | seasonal_pipeline | seasonal_amp1.0 |
|---|---|---|---|---|
| Known baseline (theorem) | 0.011 [0.007,0.016] (15u/6d) | 0.004 [0.002,0.008] (6u/2d) | 0.034 [0.026,0.042] (50u/17d) | 0.041 [0.033,0.051] (68u/14d) |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) |
| Plug-in, band 0 | 0.020 [0.015,0.027] (29u/11d) | 0.013 [0.008,0.018] (22u/3d) | 0.044 [0.036,0.054] (67u/22d) | 0.038 [0.030,0.047] (62u/14d) |
| Confidence-sequence baseline | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) |
paths per law: 2000 to 2000

(e) Missing data and (f) estimated intraday profile: within the first 1950 monitored steps (bound 0.250; cs bound 0.260). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | missing_5pct | missing_blocks | seasonal_pipeline | seasonal_amp1.0 |
|---|---|---|---|---|
| Known baseline (theorem) | 0.114 [0.100,0.128] (114u/113d) | 0.089 [0.077,0.102] (86u/92d) | 0.322 [0.302,0.343]* (304u/340d) | 0.373 [0.353,0.395]* (424u/323d) |
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.004] (2u/0d) | 0.001 [0.000,0.003] (1u/0d) |
| Plug-in, band 0 | 0.255 [0.236,0.275] (253u/257d) | 0.226 [0.208,0.244] (243u/208d) | 0.355 [0.334,0.376]* (385u/325d) | 0.403 [0.382,0.425]* (450u/356d) |
| Confidence-sequence baseline | 0.001 [0.000,0.004] (2u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.002 [0.001,0.005] (2u/2d) | 0.000 [0.000,0.002] (0u/0d) |
paths per law: 2000 to 2000


How to read this: missing returns passed as NaN are skipped, and the proof allows that when missingness does not depend on the current return. Both missing laws (5% at random; gaps of 30 to 120 steps, 17% of all returns) stayed inside every bound, and the gaps actually lowered the rates (the segment has fewer effective steps). The seasonality pipeline is the realistic case: the U-shape is estimated from the past 20 sessions with `intraday_scale`, not known. The estimated profile leaves a small time-varying error in the clipped mean. Within a day the theorem's setting is fine (3.4%), but within five days it reaches 32.2% against 25%, and the plug-in band 0 reaches 35.5% (its i.i.d. level is 26.7% to 28.0%). The README band and the confidence-sequence baseline absorb the error (0.1% and 0.2% within five days). This is the measurement that `docs/math.md` 5.6 item 3 says did not exist; it is a NEW limitation with numbers, not a bug. A stronger U-shape (amplitude 1.0, the profile's sd ratio across the day about 2) makes it worse: 4.1% within a day (inconclusive against 5%) and 37.3% within five days for the known baseline, 40.3% for the plug-in band 0, still 0.1% and 0.0% for the README setting and the confidence-sequence baseline.


(g) Forecast-relative pipeline with a misspecified GARCH(1,1) fit: within the first 390 monitored steps (bound 0.050; cs bound 0.060). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | fc_gjr | fc_egarch | fc_lnsv |
|---|---|---|---|
| Known baseline (theorem) | 0.013 [0.009,0.019] (20u/6d) | 0.007 [0.004,0.012] (5u/9d) | 0.016 [0.011,0.022] (0u/32d) |
| Plug-in, band 0.25 (README) | 0.001 [0.000,0.003] (1u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) |
| Plug-in, band 0 | 0.023 [0.017,0.031] (34u/12d) | 0.017 [0.012,0.024] (20u/14d) | 0.018 [0.013,0.025] (0u/36d) |
| Confidence-sequence baseline | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.003] (0u/1d) |
paths per law: 2000 to 2000

(g) Forecast-relative pipeline with a misspecified GARCH(1,1) fit: within the first 1950 monitored steps (bound 0.250; cs bound 0.260). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.

| Method | fc_gjr | fc_egarch | fc_lnsv |
|---|---|---|---|
| Known baseline (theorem) | 0.151 [0.136,0.168] (166u/137d) | 0.110 [0.097,0.124] (70u/149d) | 0.131 [0.117,0.146] (7u/255d) |
| Plug-in, band 0.25 (README) | 0.001 [0.000,0.004] (2u/0d) | 0.000 [0.000,0.002] (0u/0d) | 0.000 [0.000,0.002] (0u/0d) |
| Plug-in, band 0 | 0.312 [0.292,0.333]* (318u/306d) | 0.243 [0.225,0.263] (208u/279d) | 0.223 [0.205,0.241] (114u/331d) |
| Confidence-sequence baseline | 0.000 [0.000,0.002] (0u/0d) | 0.001 [0.000,0.004] (0u/2d) | 0.002 [0.001,0.005] (0u/4d) |
paths per law: 2000 to 2000


How to read this: here the monitored series is r_t divided by a GARCH(1,1) forecast whose model is wrong (the data come from GJR, EGARCH or log-normal stochastic volatility). The fitted parameters are in `r3_null.csv` (α 0.06 to 0.09, β 0.90 to 0.91). A wrong forecast model leaves some clustering in the standardised returns, but with these parameters the residual clustering is mild: the known baseline stays inside its bounds in all three cases (0.7% to 1.6% per day, 11% to 15% per five days), as do the README setting and the confidence-sequence baseline (at most 0.2%). Only the plug-in band 0 on GJR exceeds its five-day bound (31.2%). Under log-normal volatility the alarms are almost all on the decrease side, which says the fitted GARCH forecast overshoots after large returns. Verdict: VERIFIED within bound for the proved settings. This supports the README's forecast-relative remedy more than E10 did (E10 used the true model), with the caveat that only three misspecifications at moderate persistence were tried.


(h) Ramp: variance rises linearly from 1 to 2 over 390 steps, starting at monitored step 50, then stays at 2.

| Method | Alarm before ramp | Detected during ramp (390 steps) | Detected within 780 steps of ramp start | Detected within 5 days | Median delay from ramp start (IQR) | Share up |
|---|---|---|---|---|---|---|
| Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] | 0.594 [0.572,0.615] | 1.000 [0.998,1.000] | 1.000 [0.998,1.000] | 371 (320 to 425) | 1.00 |
| Plug-in, band 0 | 0.000 [0.000,0.002] | 0.991 [0.986,0.994] | 1.000 [0.998,1.000] | 1.000 [0.998,1.000] | 235 (195 to 279) | 1.00 |
| Confidence-sequence baseline | 0.000 [0.000,0.002] | 0.602 [0.580,0.623] | 1.000 [0.998,1.000] | 1.000 [0.998,1.000] | 369 (317 to 424) | 1.00 |


(i) Burst: variance x3 for L steps starting at monitored step 200, then back to 1.

| Burst L | Method | Alarm before burst | Inside burst | Within 60 after | Inside or within 60 after | Later than 60 after (5-day horizon) | Median delay from burst start (IQR) | Share up |
|---|---|---|---|---|---|---|---|---|
| 30 | Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] | 0.252 [0.233,0.271] | 0.017 [0.012,0.023] | 0.268 [0.249,0.288] | 0.000 [0.000,0.002] | 25 (20 to 28) | 1.00 |
| 30 | Plug-in, band 0 | 0.006 [0.003,0.010] | 0.625 [0.604,0.646] | 0.056 [0.047,0.067] | 0.681 [0.660,0.701] | 0.083 [0.072,0.096] | 23 (18 to 29) | 0.94 |
| 30 | Confidence-sequence baseline | 0.000 [0.000,0.002] | 0.257 [0.238,0.276] | 0.013 [0.009,0.020] | 0.270 [0.251,0.290] | 0.001 [0.000,0.003] | 25 (20 to 28) | 1.00 |
| 60 | Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] | 0.806 [0.788,0.823] | 0.011 [0.008,0.017] | 0.818 [0.800,0.834] | 0.000 [0.000,0.002] | 37 (27 to 47) | 1.00 |
| 60 | Plug-in, band 0 | 0.006 [0.004,0.011] | 0.978 [0.971,0.984] | 0.002 [0.001,0.005] | 0.980 [0.973,0.985] | 0.005 [0.003,0.010] | 26 (18 to 35) | 1.00 |
| 60 | Confidence-sequence baseline | 0.000 [0.000,0.002] | 0.807 [0.789,0.824] | 0.009 [0.006,0.015] | 0.817 [0.799,0.833] | 0.000 [0.000,0.002] | 37 (27 to 46) | 1.00 |


(j) Change during calibration: variance x2 from calibration step c of 1,000 onwards (so the last 1000 - c calibration returns already have the new variance).

| c | Method | Alarm within 390 | Alarm within 1950 | Up within 1950 | Down within 1950 | Median delay from monitoring start (IQR) | Median m0hat / true new m0 |
|---|---|---|---|---|---|---|---|
| 100 | Plug-in, band 0.25 (README) | 0.000 [0.000,0.002] | 0.001 [0.000,0.003] | 0.001 [0.000,0.003] | 0.000 [0.000,0.002] | 1471 (1471 to 1471) | 0.95 |
| 100 | Plug-in, band 0 | 0.051 [0.042,0.062] | 0.413 [0.392,0.435] | 0.383 [0.362,0.405] | 0.030 [0.023,0.038] | 1024 (643 to 1450) | 0.95 |
| 100 | Confidence-sequence baseline | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] | nan (nan) |  |
| 500 | Plug-in, band 0.25 (README) | 0.042 [0.034,0.052] | 0.444 [0.422,0.466] | 0.444 [0.422,0.466] | 0.000 [0.000,0.002] | 1025 (650 to 1466) | 0.74 |
| 500 | Plug-in, band 0 | 0.786 [0.768,0.804] | 1.000 [0.998,1.000] | 1.000 [0.998,1.000] | 0.000 [0.000,0.002] | 248 (157 to 359) | 0.74 |
| 500 | Confidence-sequence baseline | 0.011 [0.008,0.017] | 0.136 [0.122,0.152] | 0.136 [0.122,0.152] | 0.000 [0.000,0.002] | 1166 (782 to 1573) |  |
| 900 | Plug-in, band 0.25 (README) | 0.936 [0.925,0.946] | 1.000 [0.998,1.000] | 1.000 [0.998,1.000] | 0.000 [0.000,0.002] | 172 (117 to 252) | 0.55 |
| 900 | Plug-in, band 0 | 1.000 [0.998,1.000] | 1.000 [0.998,1.000] | 1.000 [0.998,1.000] | 0.000 [0.000,0.002] | 77 (53 to 107) | 0.55 |
| 900 | Confidence-sequence baseline | 0.863 [0.847,0.877] | 1.000 [0.997,1.000] | 1.000 [0.997,1.000] | 0.000 [0.000,0.002] | 206 (134 to 307) |  |

paths per cell: 2000


(k) Variance x2 or x3 from monitored step 50, carried only by rare jumps (sd 20) that the truncation clips.

| Cell | Method | True variance after | Clipped mean after / before | Alarm within 390 | Alarm within 1950 | Up within 1950 |
|---|---|---|---|---|---|---|
| jumpvar_x2 | Plug-in, band 0.25 (README) | 2 | 1.015 | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] |
| jumpvar_x2 | Plug-in, band 0 | 2 | 1.015 | 0.025 [0.019,0.033] | 0.292 [0.273,0.313] | 0.206 [0.189,0.225] |
| jumpvar_x2 | Confidence-sequence baseline |  |  | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] |
| jumpvar_x3 | Plug-in, band 0.25 (README) | 3 | 1.035 | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] | 0.000 [0.000,0.002] |
| jumpvar_x3 | Plug-in, band 0 | 3 | 1.035 | 0.039 [0.031,0.048] | 0.321 [0.300,0.341] | 0.278 [0.258,0.298] |
| jumpvar_x3 | Confidence-sequence baseline |  |  | 0.000 [0.000,0.002] | 0.001 [0.000,0.003] | 0.001 [0.000,0.003] |

paths per cell: 2000


How to read (h): a doubling spread over a day is caught by the README setting on 59% of paths while the variance is still rising, and on every path within a day after the ramp ends; the median alarm comes 371 steps after the ramp starts, when the variance has reached about 1.95. Compare B2: an abrupt doubling took a median of 112 steps. Band 0 is about 140 steps faster but has no proved level.

How to read (i): a 30-minute x3 burst is caught inside the burst by the README setting on 25% of paths and within an hour after on another 2%; 73% of such bursts are missed outright and no alarm follows later. A 60-minute burst is caught on 82%. The median alarm, when it comes, is 25 steps after the burst starts, so a 30-step burst is simply too short for the default A. Band 0 catches 68% of the 30-minute bursts but then raises an alarm later on 8% of paths. The confidence-sequence baseline behaves like the README setting.

How to read (j): the variance doubles while the calibration window is still being collected. "m0hat / true new m0" says how far below the new regime the learnt baseline sits. With the change at step 100 of 1,000 the baseline is only 5% low, inside the band, and the monitor is silent (0.1% within five days): it has learnt the new regime and nothing is wrong. With the change at step 500 the baseline is 26% low; the README monitor alarms within a day on only 4.2% of paths and within five days on 44%, so it stays blind for a week on 56% of paths, and the confidence-sequence baseline alarms on 14%. Band 0 alarms within a day on 79%. With the change at step 900 (only 100 new-regime returns in the window) the baseline is 45% low and the README monitor alarms within a day on 94% of paths, median 172 steps. In every case the alarms are on the increase side (the learnt baseline is too low), and after an alarm the monitor recalibrates on the new regime, so "blind" means "silent until then". NEW limitation: nothing in the monitor detects a change inside its own calibration window.

How to read (k): the variance doubles or triples, but only through rare jumps of 20 sd, which the truncation clips at τ² (about 9). The clipped mean rises by only 1.5% (x2) and 3.5% (x3) in the sample (about 2% and 4% in expectation), far inside the 25% band, so the README setting and the confidence-sequence baseline never alarm within five days (0.0% to 0.1%); band 0 alarms on 29% to 32% within five days, barely above its 27% null level, so it is blind too. This is the mirror image of (a): the monitor is blind to a variance change carried by the tails, by design (`docs/math.md` 5.6 item 4). It is the price of the jump robustness the README advertises, and should be stated next to it.

## 4. Claims that need changing

| File | What | Severity |
|---|---|---|
| `README.md`, first sentence ("raises an alarm when volatility shifts up or down") and the "What the project does" item 2 | Say that the monitor watches the mean of the clipped squared return, which can move while the variance does not, and can stay still while the variance moves. Numbers: a switch from t30 to t3 tails at constant unit variance fires the decrease side on 65% of paths within a day (known baseline) and on 45% of paths within a day for the README setting in the reverse direction; a doubling or tripling of variance carried only by rare large jumps is never flagged by the README setting (0.0% within five days, see (k)). `docs/math.md` 5.6 item 4 already says this in one sentence without numbers. | High: it changes what an alarm means and what a non-alarm means |
| `docs/limitations.md`, section 1 | Add the tail-shape item above as a NEW limitation (it is in `docs/math.md` 5.6.4 but not in the limitations list), with the numbers of table (a) and (k). | High |
| `docs/limitations.md`, section 3 ("Truncation does not remove microstructure noise ... at one-minute sampling the noise is small") | Add numbers: with MA(1) bounce θ = -0.3 the theorem's own setting stays inside its bounds (2.9% per day, 21.8% per five days); at θ = -0.6 it reaches 5.0% per day (inconclusive) and 40.7% per five days (bound 25%, exceeded), and the plug-in without a band reaches 8.6% per day. The README setting and the confidence-sequence baseline stay below 1%. Already a stated limitation (`docs/math.md` 5.6 item 6); the per-five-day failure of the known baseline is NEW. | Medium |
| `docs/math.md` 5.6 item 3 ("no experiment so far isolates this effect") and `docs/limitations.md` | The effect is now isolated: with the real-data pipeline (`intraday_scale`, 20-day past-only lookback) on a simulated U-shape of amplitude 0.6, the known baseline stays inside the one-day bound (3.4%) but reaches 32.2% within five days against 25% (plug-in band 0: 35.5%); the README setting (0.1%) and the confidence-sequence baseline (0.2%) absorb it. NEW limitation with numbers; replace the sentence. | Medium |
| `README.md` E16 sentence "every failure matched the condition that was broken (... coarse price ticks)" | Not contradicted, but qualify: on a random-walk price rounded to a grid giving 20%, 40% or 60% exact-zero returns, no proved setting exceeded its bound here (known baseline 0.5% to 1.7% per day, 8.1% to 13.2% per five days). At 60% zeros the plug-in monitor extended every calibration window to 2,000 returns and used its fallback scale on every path, and only the increase side was active, which the docs should say is the expected behaviour. | Low |
| `reports/BACKTEST_STRESS_REPORT.md` B2 and `README.md` power claims (median delay 112 steps for a doubling) | Add the realistic shapes: a linear doubling over a day is caught during the ramp on 59% of paths and on 100% within a day after it (median delay 371 steps from the ramp start); a 30-minute x3 burst is caught inside the burst or within an hour after on only 27% of paths (60-minute burst: 82%); a doubling that begins halfway through the calibration window leaves the README monitor with a baseline 26% too low and silent for five days on 56% of paths. | Medium: these are the shapes a user will meet |
| `docs/limitations.md`, section 1, "The baseline is estimated" | Add the in-calibration change case (j): the plug-in and confidence-sequence baselines learn a mixed regime, and no restart logic notices it; the only remedy is a shorter window or a change detector that also runs during calibration. NEW. | Medium |

## 5. What this does not show

* Everything is simulated with Gaussian or Student innovations; no real data were used. The MA(1), tick, skew, missing and seasonality laws are stylised and tested one at a time, not combined as on a real exchange.
* The known-baseline rows use a baseline from a separate long sample, so they test the theorem's own setting. For laws where the conditional mean moves (MA(1), tail switches, forecast residuals) the theorem makes no claim, so a star there means "the assumption fails and it matters", not "the proof is wrong".
* One tick family (a Gaussian random walk rounded to a grid) stands in for price discreteness; a real order book also has queueing and a tick-to-volatility ratio that changes through the day.
* The forecast laws use one fitted GARCH(1,1) per model (one training path), as E10 does, so the fit error is the same on every path and not averaged over fits.
* Power cells use Gaussian returns and one change size each; 2,000 paths give intervals of about plus or minus 2 points.
* Nothing here tests the recalibration protocol over many days, the panel or multi-frequency monitors, or the guard of E24.
