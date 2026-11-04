# Results: conformal rank e-detector against the mean-based monitor

Every number below comes from `run_experiment.py` (fixed seeds) and is stored under `results/`. Labels: [SYNTHETIC] for simulated data, [REAL] for Binance one-minute data. Method details are in `README.md`, section 6. Intervals are 95% Wilson intervals for shares and distribution-free 95% intervals (binomial order statistics) for medians.

## Summary

* **Exact validity, any law, any reference length, any horizon [SYNTHETIC].** No-change false alarm rates of the rank detector: 1.9% to 3.5% per day (bound 5%) and 22.2% to 24.3% per week (bound 35%) for Gaussian, Student t4 and jump returns, with one-day or two-hour references (Q1). The plug-in mean-based monitor breaks the weekly bound (51% to 57%) and, with a two-hour window, the daily one (42% to 46%).
* **The trap of a fixed reference [SYNTHETIC].** Ranking against the calibration window only (quantile bins learned once) makes a false alarm near certain within a week (98% to 100%), for the same reason the plug-in baseline fails (Q1c).
* **Power [SYNTHETIC].** Against the plug-in monitor: 4% to 7% slower for variance changes with rank scale bets alone, 11% to 19% with the full mixture; 2.9 times faster for a Gaussian to t4 change at fixed variance; slower for a change in jump intensity at fixed variance (a negative result); only the signed score sees a pure skew change (Q2). Against the repo's proved alternative (confidence-sequence baseline): 1.7 to 2.8 times faster, and it detects the tail change that the latter almost never sees (Q5).
* **Truncation does not matter for ranks [SYNTHETIC]** (Q3), except to bound the Gaussianised decrease bet at a fixed level.
* **Real data [REAL].** On BTCUSDT and ETHUSDT (Jul to Dec 2025) with matched tolerance (band 0.25) the two designs raise similar numbers of alarms (0.62 to 0.69 per day) and 49% to 61% of them coincide within an hour (chance 5% to 6%). Most single-design alarms happen while the other monitor is calibrating. While the main monitor was watching, the rank detector alone raised 17 (BTC) and 5 (ETH) alarms; on BTC these are shape changes at constant scale: the largest moves became about three times rarer and the smallest and zero moves more common (Q4). The main monitor alone raised 32 per coin, all volatility decreases.
* **Clustering breaks both designs alike** (Q1b); a short reference cuts the rank detector's blind time from about two thirds to about one tenth at a cost of 17% to 52% in delay (Q4, Q6).

Names used in the tables:

* **MAIN plug-in, band b**: the repo's `VolatilityMonitor` (plug-in baseline from a 1440-step calibration window, `c = 3`, two-sided Shiryaev-Roberts, each side at `2A`), run on the same paths (the vectorised copy in `tsb_sim.py` matches the class alarm for alarm, see the tests).
* **MAIN confidence-sequence baseline**: `VolatilityMonitor(baseline="cs")`, the repo's rigorous version (guarantee `delta + m / A` with `delta = 0.01`).
* **RANK omnibus (abs or signed), b**: this prototype with the frozen default configuration (scale-up bets 0.25, scale-down bets 0.25, histogram learners 0.5), score `|u|` or `u`, null band b.
* **RANK scale only**: the rank detector with the scale bets only (the main monitor's statistic computed on ranks, with an exactly known baseline).

All monitors use `A = 28,800`: under their nulls, at most 5% chance of a false alarm within one 1440-step day and at most 35% within seven days.

## Q1. Validity with no change [SYNTHETIC]

1000 paths per law; share of paths with a false alarm within 1 day (bound 5%) and, for the 1440-step calibration, within 7 days (bound 35%). File `results/q1_validity.csv`.

**Calibration window 1440 steps (one day)**

| Detector | Gaussian | Student t4 (unit variance) | Gaussian + jumps (0.5% per step, 8 sd) |
|---|---|---|---|
| MAIN plug-in, band 0 | 5.4% (4.2 to 7.0); 7 d: 57.0% (53.9 to 60.0) | 6.0% (4.7 to 7.6); 7 d: 52.5% (49.4 to 55.6) | 4.6% (3.5 to 6.1); 7 d: 51.4% (48.3 to 54.5) |
| MAIN plug-in, band 0.25 | 0.0% (0.0 to 0.4); 7 d: 0.0% (0.0 to 0.4) | 0.0% (0.0 to 0.4); 7 d: 0.0% (0.0 to 0.4) | 0.0% (0.0 to 0.4); 7 d: 0.0% (0.0 to 0.4) |
| RANK scale only (abs), b 0 | 2.1% (1.4 to 3.2); 7 d: 23.6% (21.1 to 26.3) | 2.1% (1.4 to 3.2); 7 d: 22.5% (20.0 to 25.2) | 2.0% (1.3 to 3.1); 7 d: 22.2% (19.7 to 24.9) |
| RANK omnibus (abs), b 0 | 2.4% (1.6 to 3.5); 7 d: 22.6% (20.1 to 25.3) | 3.1% (2.2 to 4.4); 7 d: 24.2% (21.6 to 27.0) | 2.4% (1.6 to 3.5); 7 d: 23.0% (20.5 to 25.7) |
| RANK omnibus (abs), b 0.25 | 0.0% (0.0 to 0.4); 7 d: 0.0% (0.0 to 0.4) | 0.0% (0.0 to 0.4); 7 d: 0.0% (0.0 to 0.4) | 0.0% (0.0 to 0.4); 7 d: 0.0% (0.0 to 0.4) |
| RANK omnibus (signed), b 0 | 2.2% (1.5 to 3.3); 7 d: 22.7% (20.2 to 25.4) | 2.0% (1.3 to 3.1); 7 d: 23.4% (20.9 to 26.1) | 1.9% (1.2 to 2.9); 7 d: 24.3% (21.7 to 27.1) |

**Calibration window 120 steps (two hours), one day**

| Detector | Gaussian | Student t4 (unit variance) | Gaussian + jumps (0.5% per step, 8 sd) |
|---|---|---|---|
| MAIN plug-in, band 0 | 42.7% (39.7 to 45.8) | 45.7% (42.6 to 48.8) | 42.0% (39.0 to 45.1) |
| MAIN plug-in, band 0.25 | 1.1% (0.6 to 2.0) | 1.9% (1.2 to 2.9) | 1.2% (0.7 to 2.1) |
| RANK scale only (abs), b 0 | 2.3% (1.5 to 3.4) | 2.7% (1.9 to 3.9) | 1.9% (1.2 to 2.9) |
| RANK omnibus (abs), b 0 | 3.5% (2.5 to 4.8) | 2.9% (2.0 to 4.1) | 2.2% (1.5 to 3.3) |
| RANK omnibus (abs), b 0.25 | 0.0% (0.0 to 0.4) | 0.0% (0.0 to 0.4) | 0.0% (0.0 to 0.4) |
| RANK omnibus (signed), b 0 | 2.4% (1.6 to 3.5) | 2.5% (1.7 to 3.7) | 2.0% (1.3 to 3.1) |

**Reading.**

* The rank detectors are valid in every row, with almost the same numbers for the three laws and for both calibration lengths: 1.9% to 3.5% per day against the 5% bound, 22.2% to 24.3% within seven days against 35%. That is what an exact, distribution-free null law buys. The exact null makes `M_t - t` a martingale, so the seven-day rate is closer to its bound (about two thirds of it) than the one-day rate (roughly half).
* The plug-in mean-based monitor at band 0 sits at its one-day bound with a one-day calibration (4.6% to 6.0%) but **breaks the seven-day bound: 51% to 57% of paths alarm within a week against at most 35%** (all lower interval limits above 48%). The reason is the baseline: each path's estimate of `m0` is off by a fixed amount, which acts as a constant drift that the bets find given enough time. With a two-hour calibration the plug-in fails within one day (42% to 46%). The rank detector has no baseline, so there is nothing to drift.
* Band 0.25 makes both designs silent on i.i.d. data (0 alarms in 1000 paths, over 7 days).

## Q1b. Volatility clustering, where both nulls are false [SYNTHETIC]

GARCH(1,1), unit unconditional variance, 500 one-day paths. Share of days with a false alarm (bound 5% if the null held). File `results/q1b_garch.csv`.

| Detector | persistence 0.95 (half-life 14 steps) | persistence 0.99 (half-life 69 steps) |
|---|---|---|
| MAIN plug-in, band 0 | 80.6% (76.9 to 83.8) | 99.0% (97.7 to 99.6) |
| MAIN plug-in, band 0.25 | 15.4% (12.5 to 18.8) | 67.0% (62.8 to 71.0) |
| RANK scale only (abs), b 0 | 76.8% (72.9 to 80.3) | 98.2% (96.6 to 99.1) |
| RANK omnibus (abs), b 0 | 68.8% (64.6 to 72.7) | 96.4% (94.4 to 97.7) |
| RANK omnibus (abs), b 0.25 | 9.8% (7.5 to 12.7) | 53.0% (48.6 to 57.3) |
| RANK histogram only (abs), b 0 | 22.6% (19.2 to 26.5) | 83.8% (80.3 to 86.8) |

**Reading.** Clustering breaks exchangeability, and the rank detector alarms on it about as often as the mean-based monitor does (slightly less often at the same nominal band, 9.8% against 15.4% for short clusters, but the two bands do not mean the same thing, so this is not a like-for-like comparison). The scale bets are what clustering triggers: the histogram learners alone alarm on 22.6% of days for short clusters. Neither design separates "a persistent cluster" from "a change"; the fix is the same for both (monitor returns relative to a volatility forecast, as in E10 of the main report).

## Q2. Power: variance, tails, jumps and skew [SYNTHETIC]

Change at monitored step 500, 4000 steps after it, 400 paths. Median delay in steps with its 95% interval; in brackets the share detected within 4000 steps when below 95%. False alarms before the change were at most 1.0% for every detector and scenario. File `results/q2_power.csv`, figure `results/fig1_synthetic_delays.png`.

| Change | MAIN plug-in, band 0 | RANK scale only (abs), b 0 | RANK omnibus (abs), b 0 | RANK omnibus (signed), b 0 |
|---|---|---|---|---|
| sd x1.25 (variance x1.56) | 127 (119 to 140) | 136 (125 to 150) | 152 (139 to 161) | 152 (138 to 162) |
| sd x1.5 (variance x2.25) | 49 (46 to 51) | 51 (48 to 55) | 55 (51 to 58) | 55 (51 to 58) |
| sd x0.8 (variance x0.64) | 211 (202 to 220) | 222 (213 to 228) | 237 (233 to 248) | 239 (234 to 250) |
| sd x0.67 (variance x0.44) | 125 (121 to 129) | 130 (125 to 133) | 139 (134 to 141) | 140 (135 to 143) |
| tails: Gaussian to t4, same variance | 1203 (1099 to 1316) | 1404 (1299 to 1679) [86%] | 417 (396 to 444) | 442 (417 to 458) |
| tails: Gaussian to t3, same variance | 448 (423 to 486) | 508 (468 to 551) | 230 (221 to 238) | 240 (233 to 251) |
| tails: t4 to Gaussian, same variance | 428 (386 to 466) | 540 (489 to 625) | 357 (341 to 377) | 373 (357 to 397) |
| jumps: 0.2% to 0.5% per step, same variance | 1393 (1253 to 1573) [88%] | 1698 (1465 to 1837) [57%] | 1688 (1521 to 1842) [54%] | 1626 (1461 to 1855) [52%] |
| jumps: 0.2% to 1% per step, same variance | 361 (338 to 388) | 418 (380 to 455) | 399 (377 to 425) | 422 (399 to 461) |
| skew: large moves 80% negative, same law of \|u\| | 2587 (2207 to 2855) [29%] | 2791 (2187 to 3149) [8%] | 2265 (1543 to 2830) [10%] | 809 (758 to 862) |

**Reading.**

* **Variance changes.** The rank scale bets cost 4% to 7% in median delay against the plug-in monitor (136 against 127 steps for sd x1.25), the price of an exact guarantee; the omnibus mixture, which spends half its wealth on the histogram learners, costs 11% to 19%. The plug-in's advantage is partly borrowed: Q1 shows its baseline error makes it anti-conservative over days.
* **Tails at fixed variance.** The omnibus detector is 2.9 times faster for Gaussian to t4 (417 against 1203 steps) and 1.9 times faster for Gaussian to t3, and also faster for the reverse change t4 to Gaussian (357 against 428). The scale bets alone are slow here, so the gain comes from the histogram learners, which see that both very small and very large moves became more frequent.
* **Jump intensity at fixed total variance: a negative result.** The mean-based monitor is faster (1393 against 1688 steps, 88% against 54% detected, for 0.2% to 0.5%; 361 against 399 for 0.2% to 1%). Holding the total variance fixed forces the diffusive variance down (from 0.93 to 0.82, or to 0.64), and that decrease of the truncated variance is what the mean-based monitor detects; it is a cleaner signal than the extra jumps, which affect only the top 0.5% to 1% of ranks.
* **Skew.** Only the signed score sees it (99% detected, median 809 steps). For every `|u|`-based detector the law of the statistic does not change at all, so their "detections" here are false alarms after step 500: 8% and 10% for the rank detectors (bound 15.6% for 4500 steps) and 29% for the plug-in, consistent with its multi-day anti-conservativeness in Q1.
* The signed score changes the median delay by at most 6% against the absolute score on the other rows (from 4% faster to 6% slower), so it is the better default when skew matters.

## Q3. Does truncation matter when the statistic is a rank? [SYNTHETIC]

Rank omnibus detector (absolute score) on the same paths with score `|u|` and with `min(|u|, 3)` (fixed `tau = 3`, so every move beyond 3 sd becomes a tie). File `results/q3_truncation.csv`.

| Scenario | score \|u\| | score min(\|u\|, 3) |
|---|---|---|
| sd x1.25, median delay | 169 (158 to 182) | 168 (158 to 181) |
| Gaussian to t4, median delay | 412 (387 to 442) | 437 (406 to 464) |
| jumps 0.2% to 1%, median delay | 413 (390 to 437) | 422 (392 to 447) |
| no change, Gaussian, false alarm within 1 day | 1.8% (1.1 to 2.8) | 1.9% (1.2 to 2.9) |
| no change, jumps 0.5% per step, false alarm within 1 day | 1.5% (0.9 to 2.5) | 1.7% (1.1 to 2.7) |

**Reading.** As section 2.6 of the README predicts, truncating the score of a rank test changes nothing that matters: validity is unchanged (ties are handled exactly), and the delays move by at most about 6%, with overlapping intervals, slightly against truncation for the tail and jump changes, where the size of the largest moves carries some information. For ranks, truncation is neither needed nor useful. (The delays differ slightly from Q2 because Q3 uses its own seeds.)

## Q1c. Why the reference must grow [SYNTHETIC]

Gaussian returns, no change, 500 paths, 14 days. The same bettors on p-values computed against the growing reference (the conformal construction of this prototype) and against the fixed reference window only ("quantile bins learned from the calibration window"). File `results/q1c_fixed_reference.csv`.

| Reference | Reference length | Detector | 1 day (bound 5%) | 7 days (bound 35%) | 14 days (bound 70%) |
|---|---|---|---|---|---|
| growing (conformal, exact) | 240 | RANK omnibus (abs), b 0 | 3.0% (1.8 to 4.9) | 22.8% (19.3 to 26.7) | 43.2% (38.9 to 47.6) |
| growing (conformal, exact) | 240 | RANK histogram only (abs), b 0 | 2.6% (1.5 to 4.4) | 26.4% (22.7 to 30.4) | 48.6% (44.2 to 53.0) |
| fixed window | 240 | RANK omnibus (abs), b 0 | 90.4% (87.5 to 92.7) | 100.0% (99.2 to 100.0) | 100.0% (99.2 to 100.0) |
| fixed window | 240 | RANK histogram only (abs), b 0 | 93.4% (90.9 to 95.3) | 100.0% (99.2 to 100.0) | 100.0% (99.2 to 100.0) |
| growing (conformal, exact) | 1440 | RANK omnibus (abs), b 0 | 2.4% (1.4 to 4.1) | 24.4% (20.8 to 28.4) | 46.6% (42.3 to 51.0) |
| growing (conformal, exact) | 1440 | RANK histogram only (abs), b 0 | 2.8% (1.7 to 4.6) | 27.4% (23.7 to 31.5) | 50.2% (45.8 to 54.6) |
| fixed window | 1440 | RANK omnibus (abs), b 0 | 21.2% (17.8 to 25.0) | 98.4% (96.9 to 99.2) | 99.8% (98.9 to 100.0) |
| fixed window | 1440 | RANK histogram only (abs), b 0 | 29.8% (26.0 to 34.0) | 98.6% (97.1 to 99.3) | 99.8% (98.9 to 100.0) |

**Reading.** This is the cleanest result of the prototype. Ranking new returns against a fixed calibration window looks harmless (each p-value is exactly uniform on its own), but the p-values share the same reference and are therefore dependent; the learning bettors find the reference's sampling error and a false alarm becomes certain: 90% to 93% of days with a 240-step reference, 21% to 30% of first days and 98% of weeks with a one-day reference. With the growing reference, every rate stays below its bound at every horizon. The fixed-window failure is the nonparametric twin of the plug-in baseline's failure in Q1: both bet against a quantity estimated once and then treated as known.

## Q5. The mean-based monitor with a proved guarantee [SYNTHETIC]

`VolatilityMonitor(baseline="cs")` (confidence-sequence baseline, `delta = 0.01`), run on the same paths as Q1 (Gaussian) and as five Q2 scenarios. File `results/q5_cs_baseline.csv`.

| | MAIN confidence-sequence baseline | RANK omnibus (abs), b 0 (from Q1, Q2) |
|---|---|---|
| no change, false alarm within 1 day | 0.0% (0.0 to 0.4) | 2.4% (1.6 to 3.5) |
| no change, false alarm within 7 days | 0.1% (0.0 to 0.6) | 22.6% (20.1 to 25.3) |
| sd x1.25, median delay | 348 (322 to 377) | 152 (139 to 161) |
| sd x0.8, median delay | 400 (383 to 419) | 237 (233 to 248) |
| Gaussian to t4, median delay | 2251 (1755 to 2974) [10% detected] | 417 (396 to 444) |
| jumps 0.2% to 1%, median delay | 1129.5 (1057 to 1226) [92% detected] | 399 (377 to 425) |
| skew at fixed law of \|u\| | 0% detected (correctly silent) | 10% "detected" (false alarms, see Q2) |

**Reading.** Among detectors whose guarantee holds with an estimated or absent baseline, the rank detector is much faster: 2.3 times for a 25% sd increase, 1.7 times for a 20% sd decrease, 2.8 times for the jump change, and it detects the t4 change that the confidence-sequence monitor almost never sees within 4000 steps. The confidence-sequence monitor pays for estimating `m0` with a wide, conservative interval (0.1% false alarms in a week against a 35% budget); the rank detector spends the budget instead (22.6%), because its null law is known exactly.

## Q6. The price of a short reference window [SYNTHETIC]

Same monitored data (change at monitored step 500, 4000 steps after it, 400 paths) with a reference window of 120, 240 or 1440 steps before monitoring starts. Median delay (95% interval), [share detected when below 95%], and the share of paths with a false alarm in the 500 pre-change monitored steps (bound 500 / 28,800 = 1.7% under a correct null). File `results/q6_short_reference.csv`.

| Change | Reference 120 | Reference 240 | Reference 1440 |
|---|---|---|---|
| **RANK omnibus (abs), b 0** | | | |
| sd x1.25 | 193 (176 to 212); FA 0.8% | 188.5 (172 to 204); FA 0.5% | 153 (142 to 169); FA 0.5% |
| sd x0.8 | 270 (253 to 287); FA 0.5% | 259 (249 to 273); FA 0.8% | 230 (225 to 237); FA 0.2% |
| Gaussian to t4 | 623 (541 to 713) [88%]; FA 0.5% | 565 (513 to 619) [94%]; FA 0.5% | 409 (386 to 433); FA 0.2% |
| jumps 0.2% to 1% | 588 (516 to 654) [94%]; FA 0.5% | 554 (497 to 582); FA 0.2% | 423 (401 to 452); FA 0.2% |
| **MAIN plug-in, band 0** | | | |
| sd x1.25 | 123.5 (108 to 139) [88%]; FA 11.5% | 121 (110 to 135) [94%]; FA 5.5% | 130 (123 to 139); FA 0.2% |
| sd x0.8 | 205 (194 to 218) [86%]; FA 13.5% | 199 (188 to 213) [94%]; FA 6.2% | 204 (196 to 211); FA 0.8% |
| Gaussian to t4 | 926.5 (803 to 1053) [70%]; FA 12.8% | 1056 (973 to 1157) [75%]; FA 6.0% | 1243 (1155 to 1406) [94%]; FA 0.8% |
| jumps 0.2% to 1% | 357.5 (333 to 399) [86%]; FA 14.0% | 373 (335 to 407) [93%]; FA 7.2% | 379.5 (353 to 402); FA 0.5% |

**Reading.** With a two-hour reference the rank detector stays valid (at most 0.8% false alarms before the change) and is 17% to 52% slower than with a one-day reference (the early p-values are coarse, and post-change returns soon make up a large part of the growing reference, which dilutes the contrast). The plug-in monitor with the same short window alarms falsely on 11.5% to 14.0% of paths within 500 steps (bound 1.7%) and its apparent speed is partly false alarms. A short reference is therefore a real option only for the rank detector; it trades detection speed for much less blind time after each alarm (Q4, item 6).

## Q4. Real data: Binance BTCUSDT and ETHUSDT, July to December 2025 [REAL]

184 days, 264,960 minutes per coin, no missing minutes (loader checks). Exact-zero one-minute returns: 4.3% of BTC minutes and 1.0% of ETH minutes. Each alarm of a rank detector is compared with the alarms of the main monitor at band 0.25 (its "partner"), and each main-monitor alarm with those of the rank detector at `b = 0.25`. "Matched": the partner alarms within 60 minutes. "Chance": the matched share after shifting the alarm times by random offsets of 1 to 30 days (200 shifts). Unmatched alarms are split by whether the partner was monitoring or calibrating at that moment. Labels give the bettor family with the largest share of `M_t` at the alarm. Files `results/q4_real_summary.csv`, `q4_real_events.csv`, `q4_real_profiles.csv`, `q4_real_random_times.csv`, `q4_real_rank_only_alarms.csv`; figures `fig2_real_alarm_profiles.png` and `fig3_real_density_ratio.png`.

| Coin | Detector | Alarms (per day) | Share of time monitoring | Up / down / shape | Matched within 60 min (chance) | Unmatched: partner watching / calibrating |
|---|---|---|---|---|---|---|
| BTC | MAIN plug-in, band 0 | 153 (0.83) | 16% | 67 / 86 / n.a. | 20 = 13% (5.8%) | 46 / 87 |
| BTC | MAIN plug-in, band 0.25 | 127 (0.69) | 31% | 72 / 55 / n.a. | 62 = 49% (5.7%) | 32 / 33 |
| BTC | RANK omnibus (abs), b 0 | 155 (0.84) | 15% | 57 / 26 / 72 | 26 = 17% (5.9%) | 39 / 90 |
| BTC | RANK omnibus (abs), b 0.25 | 123 (0.67) | 33% | 69 / 12 / 42 | 62 = 50% (5.7%) | 17 / 44 |
| BTC | RANK omnibus (abs, nonzero only), b 0.25 | 116 (0.63) | 33% | 66 / 10 / 40 | 53 = 46% (5.5%) | 8 / 55 |
| BTC | RANK omnibus (signed), b 0.25 | 114 (0.62) | 38% | 70 / 17 / 27 | 64 = 56% (5.8%) | 3 / 47 |
| BTC | RANK omnibus (abs), b 0.25, 120-minute reference | 255 (1.39) | 88% | 142 / 40 / 73 | 58 = 23% (5.8%) | 31 / 166 |
| ETH | MAIN plug-in, band 0 | 152 (0.83) | 17% | 75 / 77 / n.a. | 16 = 11% (5.1%) | 50 / 86 |
| ETH | MAIN plug-in, band 0.25 | 122 (0.66) | 33% | 76 / 46 / n.a. | 69 = 57% (5.1%) | 32 / 21 |
| ETH | RANK omnibus (abs), b 0 | 148 (0.80) | 20% | 76 / 29 / 43 | 31 = 21% (5.4%) | 48 / 69 |
| ETH | RANK omnibus (abs), b 0.25 | 114 (0.62) | 38% | 75 / 14 / 25 | 69 = 61% (5.7%) | 5 / 40 |
| ETH | RANK omnibus (abs, nonzero only), b 0.25 | 114 (0.62) | 37% | 71 / 16 / 27 | 66 = 58% (5.7%) | 5 / 43 |
| ETH | RANK omnibus (signed), b 0.25 | 114 (0.62) | 38% | 77 / 13 / 24 | 63 = 55% (5.4%) | 4 / 47 |
| ETH | RANK omnibus (abs), b 0.25, 120-minute reference | 204 (1.11) | 91% | 131 / 35 / 38 | 63 = 31% (5.5%) | 15 / 126 |

(For the main-monitor rows the partner is the rank detector at `b = 0.25`; for the rank rows it is the main monitor at band 0.25.)

**What the alarms look like.** Medians over the 240 minutes before each alarm, against the segment's calibration window (from `results/q4_real_profiles.csv`; random times use the same reference delays):

| Coin | Group | n | Robust scale ratio | Truncated variance ratio | Shape distance (KS) | Zero returns: reference / window | Moves beyond 3 robust sd: reference / window |
|---|---|---|---|---|---|---|---|
| BTC | RANK b 0.25, matched | 62 | 1.116 | 1.214 | 0.065 | 2.8% / 2.5% | 3.6% / 3.8% |
| BTC | RANK b 0.25, rank only, main monitor watching | 17 | 0.994 | 1.037 | 0.131 | 3.8% / 5.0% | 3.5% / 3.8% |
| BTC | RANK b 0.25, rank only, main monitor calibrating | 44 | 0.860 | 1.017 | 0.073 | 1.5% / 2.1% | 3.3% / 4.2% |
| BTC | MAIN band 0.25, main only, rank detector watching | 32 | 0.725 | 0.500 | 0.062 | 1.2% / 1.0% | 3.2% / 1.7% |
| BTC | random times | 2000 | 0.974 | 0.922 | 0.076 | 2.3% / 1.7% | 3.3% / 2.5% |
| ETH | RANK b 0.25, matched | 69 | 1.081 | 1.278 | 0.063 | 0.6% / 0.4% | 3.0% / 4.2% |
| ETH | RANK b 0.25, rank only, main monitor watching | 5 | 0.636 | 0.574 | 0.064 | 0.7% / 1.7% | 2.6% / 4.1% |
| ETH | MAIN band 0.25, main only, rank detector watching | 32 | 0.817 | 0.521 | 0.057 | 0.3% / 0.4% | 2.9% / 1.9% |
| ETH | random times | 2000 | 0.990 | 0.927 | 0.058 | 0.5% / 0.4% | 3.0% / 2.1% |

The learned density ratio at the alarm (`fig3_real_density_ratio.png`) shows the shape directly. For the 17 BTC alarms raised only by the rank detector while the main monitor was watching, the median ratio of new to reference frequency is 0.32 for the largest 0.2% of move sizes and 0.57 for the next 0.8%, about 1 in the middle, and 1.49 and 1.52 for moves smaller than 97% to 99.8% of the earlier ones: **the very largest moves became about three times rarer and very small moves about one and a half times as common, at an unchanged median move size** (robust scale ratio 0.994). For alarms matched by the main monitor the pattern is the familiar volatility increase (more large moves, fewer small ones). The 17 alarms are listed in `results/q4_real_rank_only_alarms.csv`.

**Reading.**

1. **With the same nominal tolerance the two designs agree more than they differ.** At band 0.25 both raise 0.62 to 0.69 alarms per day, and 49% to 61% of alarms are matched by the other design within an hour, against 5% to 6% by chance.
2. **Most "only one design" alarms are an artefact of blind time.** Each monitor spends about two thirds of the time calibrating after its alarms (monitored share 31% to 38%). Of the rank detector's 61 unmatched BTC alarms, 44 fell while the main monitor was calibrating; for ETH 40 of 45.
3. **What only the two-sample test raises.** While the main monitor was watching, the rank detector raised 17 BTC alarms and 5 ETH alarms that it did not. The BTC ones are mostly shape changes at constant scale: fewer large moves, more tiny and zero moves, a scale-free shape distance of 0.131 against 0.076 at random times. In words, the market became "stickier" (more minutes without a price change, fewer jumps) without the typical move size changing, which a variance monitor with a 25% tolerance does not see. Skipping zero returns halves this group (8 alarms), so about half of it is a change in the share of zero returns, that is a microstructure or liquidity change rather than a volatility change. The weekend share of these 17 alarms (18%) is below that of random times (30%), so they are not simply weekend lulls.
4. **What only the mean-based monitor raises.** While the rank detector was watching, the main monitor raised 32 alarms on each coin that the rank detector did not; they are volatility decreases (truncated variance ratio 0.50 and 0.52). The rank detector's band penalty `b E_U|w|` is relatively heavier on its decrease bets, so at `b = 0.25` it is the less sensitive design for decreases. The two bands do not mean the same thing, so this is a statement about these two settings, not about the designs in general.
5. **Signed scores add little here.** The signed rank detector raised 3 (BTC) and 4 (ETH) alarms that the main monitor missed while watching, fewer than the absolute score: no sign of skew-only changes in these data at this alarm budget.
6. **A short reference is now allowed.** The rank detector with a 120-minute reference is still exact under exchangeability and monitors 88% to 91% of the time instead of about a third. It raises 1.11 to 1.39 alarms per day (1.22 to 1.57 per monitored day, below the 1.65 to 2.04 of the one-day reference). The plug-in monitor cannot do this (Q1: 42% false alarms per day with a two-hour window).
7. **The real series are far from the null.** Even the banded detectors raise about two alarms per monitored day where the null allows at most 0.05. As in the main report, this is volatility clustering and genuine regime change together, and neither design can separate them on raw deseasonalised returns (Q1b).

## Ranking and when to use this design

**Honest ranking.** *Novelty: medium.* Conformal test martingales and conformal Shiryaev-Roberts procedures exist (Vovk et al. 2003, 2021; Volkhonskiy et al. 2017), and the two-sample betting framework is Shekhar and Ramdas (2024). What is new here is the combination for this problem: a Gaussianised, truncated rank bet that is the main monitor's statistic with an exactly known baseline; a null band on the density of the conformal p-values; the side-by-side demonstration that a fixed reference and a plug-in baseline fail in the same way; and alarm-level density-ratio diagnostics on intraday returns. *Validity: the strongest in the repository under its null.* Exact for any return law, any reference length and any horizon, which neither the plug-in (fails over days) nor the confidence-sequence baseline (valid but conservative) achieves. Its weak point is the null itself: exchangeability excludes volatility clustering, which real one-minute returns have, so on raw deseasonalised returns it is no more trustworthy than the main monitor. *Practical relevance: medium.* On BTC and ETH it mostly agrees with the main monitor, adds a handful of shape alarms in six months (fewer large moves and more tiny or zero moves at constant scale), misses some volatility decreases at the chosen band, and its short-reference option removes most blind time. *Feasibility: high.* About 80 microseconds per step in pure Python against about 10 for the main monitor (`results/t_timing.csv`), no parameter estimated, defaults frozen from one pilot. *Demonstration value: high.* Q1, Q1c and Q5 make the cost of estimating a baseline visible in a few tables, and the t4 result shows what a variance monitor cannot see. *Interest to a sequential-analysis researcher: high.* The "reference sample" problem in nonparametric change detection, the Polya-urn view of why the reference must grow, the e-detector guarantees in the filtration of the p-values, the rank analogue of truncation, and an open efficiency question (the Gaussianised score is a normal-scores statistic, so we conjecture that its asymptotic efficiency against the mean-based detector under Gaussian data is 1 as the reference grows; whether the measured 4% to 7% gap is a finite-reference effect is not yet checked).

**When to prefer it over the mean-based monitor.** (1) When the baseline cannot be trusted: short or recent calibration windows, long monitoring without recalibration, or a need for a guarantee that holds for days, not hours. (2) When the change of interest is in shape rather than level: tails, the share of jumps or zero returns, and (with the signed score) skew. (3) When the monitored series is close to exchangeable, for example residuals of a good volatility forecast. Prefer the mean-based monitor when the question is the level of volatility, when the data cluster and a tolerance band on a moment is the natural null, when an alarm must come with a direction and a size, and for decreases of volatility at a given tolerance.
