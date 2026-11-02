# Results: quarticity-stream e-detector

All numbers come from `run_experiment.py` (seeds fixed) and are stored in `results/`. [SYNTHETIC] marks simulated data, [REAL] marks Binance spot one-minute bars. Settings are in `README.md`, section 6. Short glossary: "truncated" means clipped at `tau = 3 x scaled MAD`; "stream" is the per-step statistic fed to the detector (`min(r^2, tau^2)` for variance, `min(r^4, tau^4)` for quarticity); "A" is the alarm level, chosen so that the chance of a false alarm within one day is at most 5%; "delay" is the number of steps from the change to the alarm.

## Q1. Is the quarticity-stream detector valid? [SYNTHETIC]

No change, 300 paths of 20,000 steps, `A = 7800`. Bounds: false alarm within 390 steps at most 5%, average run length at least 7800 (run lengths censored at 20,000, so the averages below are underestimates).

| Scenario | Stream | Side | False alarm within 390 steps | Mean run length (censored) |
|---|---|---|---|---|
| Gaussian | variance | up / down | 1.7% / 1.0% | 9269 / 8761 |
| Gaussian | quarticity | up / down | 1.7% / 0.0% | 9471 / 9046 |
| jumps (p = 0.005, 8 sd) | variance | up / down | 1.7% / 0.7% | 8867 / 8799 |
| jumps | quarticity | up / down | 1.7% / 0.0% | 9348 / 9334 |
| Student t5 | variance | up / down | 1.7% / 0.3% | 9350 / 8392 |
| Student t5 | quarticity | up / down | 2.7% / 0.0% | 8858 / 8558 |

**Reading.** Both bounds hold in every cell, with the model-free `q0` (mean of `min(r^4, tau^4)` in the calibration window) and the bound `tau^4`. The decrease side of the quarticity stream never false-alarmed within a day: its bets are capped at `1 / (B_q - 1) = 0.037`, so it is conservative. Validity is as expected from the theory (`README.md`, section 2): the proof only needs a nonnegative bounded statistic.

## Q2. Pure variance change: how much slower is the quarticity stream? [SYNTHETIC]

Change at step 500, `A = 7800`, 300 paths, the same eight-point bet grid for both streams. Median delay in steps (share of paths with a false alarm before the change was 0% to 4.3% in every cell, and 96% to 100% of paths detected).

| Change | Variance stream, Gaussian | Quarticity stream, Gaussian | Variance stream, jumps | Quarticity stream, jumps |
|---|---|---|---|---|
| sd x 1.25 | 96 | 87.5 | 92 | 91 |
| sd x 1.5 | 38 | 33 | 40 | 37 |
| sd x 2 | 18 | 13 | 18 | 14 |
| sd x 0.8 | 153 | 271 | 151 | 299 |
| sd x 0.67 | 95 | 201.5 | 99 | 216 |
| sd x 0.5 | 68.5 | 168 | 74 | 189 |

**Reading.** The expected "fourth powers are noisier, so slower" is only half right. For *increases* the quarticity stream is 1% to 28% faster: a sd ratio `v` moves its mean by about `v^4` against `v^2`, which more than pays for the extra noise (theory: snr 0.58 against 0.67 at `v = 2`, but the detector's growth with a grid of bets favours the larger mean shift). For *decreases* it is 1.8 to 2.6 times slower, as predicted by the bet cap (`0.037` against `0.124`): the quarticity stream drops much further (`E[Y] = 0.065` against `0.25` at sd x 0.5) but the bettor may only stake a third as much. Jumps in the null change little for either stream.

## Q3. Variance held fixed, shape changed [SYNTHETIC]

Change at step 500, 20,000 monitored steps, `A = 400,000` (null false alarm within the horizon at most 5%; the realised pre-change false alarm share was at most 1% in every cell), 200 paths. "Detected" is the share of paths that alarmed within 20,000 steps after the change; the median delay is over those paths.

| Alternative (variance fixed at 1) | Variance stream, up | Variance stream, down | Quarticity stream, up | Quarticity stream, down |
|---|---|---|---|---|
| Student t3 | 0% | **100%** (650) | 6% (13,847) | 32% (9592) |
| Student t4 | 0% | **99.5%** (2033) | 57.5% (8441) | 2% |
| Student t5 | 0% | **96.5%** (3937) | 71.5% (6461) | 1.5% |
| Student t8 | 3% | 57.5% (8448) | **75%** (7626) | 0% |
| mixture d = 0.5 (variance 0.5 or 1.5) | 8% | 28% (9669) | **82.5%** (4915) | 0.5% |
| mixture d = 0.7 | 5% | 46% (7712) | **99.5%** (2678) | 0% |
| mixture d = 0.9 | 2% | 58% (8153) | **99%** (1033) | 0% |
| slow SV, log-sd 0.3 | 100% (715) | 100% (524) | 100% (664) | 100% (2113) |
| slow SV, log-sd 0.5 | 100% (496) | 100% (271) | 100% (442) | 100% (1205) |

Theory (`results/q3_theory.csv`), mean of the normalised stream and per-step signal-to-noise ratio, clipping at `a` pre-change standard deviations:

| Alternative | a = 3: variance `E[Y]` (snr) | a = 3: quarticity `E[Y]` (snr) | a = 5: variance | a = 5: quarticity |
|---|---|---|---|---|
| Student t3 | 0.73 (-0.18) | 0.95 (-0.01) | 0.83 (-0.07) | 1.96 (+0.07) |
| Student t5 | 0.92 (-0.06) | 1.12 (+0.03) | 0.98 (-0.01) | 1.70 (+0.06) |
| Student t8 | 0.97 (-0.02) | 1.11 (+0.03) | 1.00 (-0.00) | 1.39 (+0.05) |
| mixture d = 0.5 | 0.99 (-0.01) | 1.14 (+0.04) | 1.00 (-0.00) | 1.25 (+0.05) |
| mixture d = 0.9 | 0.96 (-0.03) | 1.45 (+0.10) | 1.00 (-0.00) | 1.80 (+0.10) |

**Reading.**

* **Fatter tails at fixed variance (Student t) are caught by the variance stream, as a decrease, not by the quarticity stream.** This is a surprise worth stating plainly. With `tau` fixed at 3 pre-change standard deviations, a Student t has more mass near zero and its tail mass is clipped, so the truncated variance falls (by 27% for t3, 8% for t5), while the truncated quarticity rises only modestly (by 9% to 12% for t4 to t8, and not at all for t3) because clipping removes most of a fat tail's fourth moment. The quarticity stream's noise is also 2.5 times larger. Only at `df = 8`, where the variance effect fades, does the quarticity stream win (75% against 57.5%). Raising the clipping level to `a = 5` reverses the picture in theory (quarticity snr 0.06 to 0.07 against variance 0.01 to 0.07), at the cost of less protection against jumps. The clipping level is therefore a design choice that decides *which* shape changes are visible, and should be reported with any result.
* **Volatility of volatility at fixed mean variance is what the quarticity stream detects and the variance stream cannot.** For the two-point mixtures the quarticity stream alarms on 82.5% to 99% of paths, with median delays of 4915, 2678 and 1033 steps (about 13, 7 and 3 trading days of 390 steps), while the variance stream's increase side stays at 2% to 8% (the null rate) and its decrease side alarms slowly on 28% to 58% of paths, in the wrong direction. The variance stream reads "more flicker" as "slightly less variance" (the clipped share of the high state).
* **A change that moves the variance up and down over time (slow SV) is detected by both streams in both directions**, because the conditional null is violated by construction; the quarticity stream is 7% to 11% faster on the increase side and 2.5 to 4.5 times slower on the decrease side, as in Q2.
* **Shape detection is slow.** Signal-to-noise ratios of 0.03 to 0.10 per step give expected log growth of 0.0005 to 0.005 per step, so reaching `log(400,000) = 12.9` takes thousands of steps. At one-minute sampling that is hours to days; at one second it is minutes.

## Q4. Real data: Binance BTCUSDT and DOGEUSDT, 1 July to 31 December 2025 [REAL]

264,960 minutes per coin, no missing minutes, zero-return share 4.3% (BTC) and 3.3% (DOGE). Returns deseasonalised by the 20-day intraday profile. `A = 28,800` per side, 1440-minute calibration, recalibration after each alarm. "Post / baseline" compares the next 1440 minutes with the calibration window using the same `tau`.

| Coin | Monitor | Band | Alarms (per day) | Up / down | Median post / baseline truncated variance, up | Share of up alarms followed by higher variance | Median post / baseline truncated variance, down | Share of down alarms followed by lower variance | Share within 60 min (240 min) of a variance-monitor alarm |
|---|---|---|---|---|---|---|---|---|---|
| BTCUSDT | variance | 0.25 | 127 (0.69) | 74 / 53 | 1.33 | 77% | 0.51 | 83% | (itself) |
| BTCUSDT | quarticity | 0.25 | 131 (0.71) | 73 / 58 | 1.45 | 80% | 0.59 | 84% | 62% (73%) |
| BTCUSDT | quarticity | 0.5625 | 101 (0.55) | 68 / 33 | 1.46 | 79% | 0.29 | 100% | 32% (50%) |
| BTCUSDT | control, 2000 random times | | | | 1.00 | 50% | 1.00 | 50% | |
| DOGEUSDT | variance | 0.25 | 126 (0.69) | 73 / 53 | 1.32 | 72% | 0.54 | 83% | (itself) |
| DOGEUSDT | quarticity | 0.25 | 129 (0.70) | 67 / 62 | 1.38 | 76% | 0.66 | 82% | 60% (73%) |
| DOGEUSDT | quarticity | 0.5625 | 101 (0.55) | 67 / 34 | 1.45 | 87% | 0.43 | 97% | 35% (51%) |
| DOGEUSDT | control, 2000 random times | | | | 1.00 | 50% | 1.00 | 50% | |

Quarticity ratios (post / baseline truncated quarticity): median 1.67 (BTC) and 1.64 (DOGE) after up alarms, 0.43 and 0.54 after down alarms, with band 0.25; the control gives 1.00 and 1.01. Monthly alarm counts are flat (20 to 23 per month per monitor at band 0.25, `results/q4_real_alarms_per_month.csv`).

**What the quarticity alarms coincide with.**

* **Mostly the same events as the variance monitor.** With band 0.25, 87% (BTC) and 87% (DOGE) of quarticity alarms fall on a UTC day that also has a variance alarm, 60% to 62% within an hour of one. Every up alarm of every monitor had at least one clipped return in the preceding hour: alarms are triggered by bursts of large moves, on both streams.
* **The alarms far from any variance alarm (more than 4 hours away; 35 for BTC and 34 for DOGE at band 0.25) are still followed by moves in the flagged direction** on 80% (BTC) and 88% (DOGE) of occasions, with median post / baseline variance ratios of 1.42 and 1.30 after the up alarms and 0.85 and 0.66 after the down alarms. So they are not noise, but they are not "shape at fixed level" events either: the level moved too, the variance monitor simply had not reached its threshold or had just recalibrated.
* **A pure shape change at fixed variance was not observed on this data.** The ratio `q_ratio / var_ratio^2` (1 for a scale change with the shape fixed) has median 0.81 to 0.90 after up alarms and 1.17 to 1.38 after down alarms for *every* monitor, including the variance one. This is the mechanical effect of a fixed `tau`: when the variance rises, clipping removes a larger share of the fourth moment than of the second. It should not be read as evidence of lighter tails after an increase.
* **The one dated event.** All three monitors on both coins alarmed between 15:28 and 15:31 UTC on 10 October 2025, with truncated variance 2.6 to 3.5 times the baseline and truncated quarticity 4.0 to 7.4 times over the next 24 hours, followed by decrease alarms on 11 October. This is the crypto market crash of 10 to 11 October 2025 after the United States announced 100% tariffs on Chinese imports (checked online; sources in `results/q4_real_event_check.csv`).
* **The band matters more than the stream.** With the same tolerance on the variance (band 0.5625 on quarticity), the quarticity monitor alarms 20% less often than the variance monitor and its alarms are stronger (100% and 97% of down alarms followed by lower variance, median ratio 0.29 and 0.43).

**Caveats [REAL].** One venue, two coins, six months. Both monitors fire roughly once a day because crypto volatility genuinely drifts, so the recalibration cadence, not the statistic, sets the alarm rate. "Followed by" is a sanity check against the random-time control, not ground truth. The real-data run does not separate level from shape; Q3 shows that the separation is possible on synthetic data but needs thousands of steps.

## Answers to the four questions

1. **Valid:** yes, with the model-free `q0` and the bound `tau^4`; the bounds hold in all 12 null cells.
2. **Pure variance change:** 1% to 28% *faster* for increases, 1.8 to 2.6 times slower for decreases (bet cap).
3. **Shape at fixed variance:** detects vol-of-vol changes that the variance stream cannot (82% to 99% against 2% to 8%); does *not* beat the variance stream on fat tails at `c = 3`, where the variance stream sees a decrease instead. The clipping level decides what is visible.
4. **Real data:** alarm rates match the variance monitor (0.70 against 0.69 per day); 87% of its alarms share a day with a variance alarm; the rest are still level changes. No pure shape event was found in six months of two coins.

## Negative results and failed attempts

* A first version of Q3 used `A = 7800` with a 10,000-step horizon. Under the null the Shiryaev-Roberts run length is close to `A`, so about half of the paths alarmed within the horizon with no change at all and every "detection" share was uninterpretable. The experiment was rerun with `A = 400,000` (`docs/iterations/ALT1.md`).
* The hypothesis "the quarticity stream detects fat tails faster" is false at `c = 3`.
* The real-data "shape ratio" turned out to measure the truncation mechanics, not the shape.

## Ranking

* **Research novelty:** moderate. The detector is the main package's with a different statistic; the new content is the null-mean and bet-range derivation, the finding that the clipping level decides which shape changes are visible, and the observation that a variance monitor reads fat tails as a volatility decrease.
* **Mathematical validity:** high. Same proof, confirmed in 12 null cells.
* **Practical relevance:** moderate. Vol-of-vol changes are detectable but slow (thousands of steps); at one-minute sampling that is days. The decrease side is weak by construction.
* **Feasibility:** high. About 250 lines on top of the package, no new dependencies, runs in five minutes.
* **Demonstration value:** high for the synthetic story (Q2 and Q3 figures), low for the real-data story (the two monitors behave alike).
* **Likely interest to a sequential-analysis and high-frequency econometrics researcher:** the Q3 interaction between truncation level and detectable alternatives, and the side-by-side level/shape monitor, are the parts worth a conversation; the real-data result is a check, not a finding.
