# R2. Detection speed at equal false-alarm rates [SYNTHETIC]

Agent R2, round 4, 4 November 2026. This extends yesterday's B2 (`reports/BACKTEST_STRESS_REPORT.md`). B2 compared delays with every method at its nominal threshold. Here every method is first held to the same measured false-alarm rate. Status: **VERIFIED**. All parts ran in full, with more paths than planned. I checked my vectorised code against the package's own monitors: 0 mismatches. All data are simulated.

## 1. Summary in plain language

1. **When false-alarm rates are equal, the guarantee costs almost nothing in the statistic.** Each method got the threshold that makes 5% of no-change Gaussian paths alarm within the first day (390 steps). For a doubling of variance, the four e-detector variants then need a median of 37 to 44 steps. The Gaussian likelihood ratio, Page CUSUM and z-test need 38 to 43. Over the four change sizes, the two proved variants (known baseline, confidence-sequence baseline) are between 17% faster and 33% slower than the best classical rule. [VERIFIED, SYNTHETIC]
2. **The README setting is slow because of its threshold, not its statistic or its band.** At A = 7,800 it raised 0 false alarms in 20,000 no-change days. It needs 114, 515, 265 and 860 steps (variance x2, x1.5, x0.5, x2/3). The same statistic at its 5% threshold needs 37, 99, 73 and 134 steps: 3.1 to 6.4 times faster. When rates are equal, the band even helps: 37 against 43 steps without it for x2, and 73 against 93 for x0.5. The band hurts only together with the proved threshold. Together they push the real rate from the allowed 5% to 0.0%. [VERIFIED]
3. **What a user with a proof actually pays.** The confidence-sequence monitor is the only proved option when the baseline is estimated. At its proved threshold it is 3.0 to 5.5 times slower than the best classical rule tuned to 5% on the right law (Gaussian: 114 against 38 steps for x2). The known-baseline e-detector at A uses 1.0% of its 5% budget and is 1.2 to 1.7 times slower. Against classical rules tuned to the same 1%, it is as fast: 52 against 51 to 58 steps for x2. [VERIFIED]
4. **On heavy tails the classical rules are slow even with perfect tuning.** Tuned to 5% on Student t3 data itself, for a doubling of variance: e-detectors 74 to 91 steps, Page CUSUM 203, Gaussian likelihood ratio 249, z-test 750. For a halving, the likelihood ratio finds 7% of changes within five days and the z-test 8%. Every e-detector finds 95%. With jumps: likelihood ratio 554 and Page 354 steps, against 39 to 46. [VERIFIED]
5. **Tuning on the wrong law is a disaster for the classical rules.** Tuned on Gaussian data and run on Student t3, they raise a false alarm within a day on 99.9% of paths (with jumps: 63% to 67%). So they almost never see the real change. The e-detectors at A stay at 0.0% to 2.6%. [VERIFIED]
6. **The proved threshold is what protects the e-detectors.** If their thresholds are also tuned on Gaussian simulations instead of using A, they reach 6.3% to 9.8% on Student t3 against the 5% target (4.7% to 5.4% with jumps). That is far better than the classical rules, but no longer within the target. [VERIFIED]
7. **Caveat on the matching horizon.** Tuned to 5% within the first day, the e-detectors reach 35% to 49% false alarms over five days. The classical rules reach 21% to 25%, because the e-detector statistic builds up over time. Matched instead to 25% over five days on Gaussian data, the e-detectors need 1.1 to 1.4 times the best classical delay for x2 and 1.2 to 1.8 times for x0.5. On Student t3 and jumps they stay 5.6 to 11 times faster for x2. [VERIFIED]

Words used in this report:
* A **false alarm** is an alarm on a path with no change.
* The **false-alarm rate** is the share of no-change paths that alarm within the first 390 monitored steps (one trading day of minutes).
* An **oracle threshold** is tuned by simulation on the very law that is then tested. No real user can do this.
* **Detected** means the first alarm came after the change and within 1,950 steps (five days) of it.
* **Delay** is the alarm step minus the change step, plus 1, counted only on detected paths.
* The **median** is the middle delay. The **IQR** (interquartile range) is the range of the middle half of delays.
* A **Wilson 95% interval** is the range of rates that the data cannot rule out at the 95% level.

## 2. What was run

**Design.** Three laws:
* Gaussian.
* Student t with 3 degrees of freedom, scaled to variance 1.
* Gaussian plus jumps (probability 0.5% per step, jump sd 8).

All paths come from `avmonitor.simulate.change_returns`, as in B2. Each path has 1,000 calibration steps, then monitoring. A change multiplies the variance by 2, 1.5, 0.5 or 2/3 at monitored step 400. For the jump law it multiplies only the diffusive part.

**Methods and how each alarms.** These calls are copied from `reports/stress/b2_power_delay.py`:
* README: `VolatilityMonitor(cal_len=1000, null_band=0.25)`.
* Band 0: the same with `null_band=0`.
* Confidence-sequence: `baseline="cs", delta=0.01`.
* Known baseline: `detectors.betting_edetector`, up and down sides, with the true tau and m0 of the law, taken from a 5,000,000-draw sample as in B1.
* Gaussian LR: `GaussianLRMonitor`.
* Page CUSUM (up only): `detectors.page_cusum` on squared returns divided by the calibration MAD scale, drift 0.25.
* z-test: `detectors.repeated_window_ztest`, window 60.

The e-detectors and the likelihood ratio alarm when max(up statistic, down statistic) reaches 2A = 15,600. Page CUSUM has no decrease side, so its decrease cells are n/a.

**Tuning.** Every alarm rule here means "alarm at the first step where a running statistic reaches h". So I ran each statistic on each path without stopping and kept its running maximum. The oracle threshold is the value that exactly 5% of the tuning paths reach within 390 steps (an order statistic). The first alarm for any threshold comes from that same simulation. Thresholds were set on the tuning paths. Every rate below is measured on separate fresh paths.

**Paths.**
* 20,000 tuning no-change paths per law.
* 20,000 fresh no-change paths per law.
* 2,000 paths per change cell (12 cells).

This is more than the brief's 1,000 / 1,000 / 500. The vectorised code ran the whole grid in 291 seconds on one core.

**Validation** (`r2_validate.py`, log `validate.log`). I compared my code with the package's own streaming monitors and array detectors on 40 paths in each of 4 cells: Gaussian x2, t3 x0.5, jumps with no change, Gaussian with no change. Each method was run at the nominal threshold and at a low threshold. That gives 48 comparisons, with **0 mismatched first alarms**.

Cross-checks with yesterday's runs:
* README at A, Gaussian x2: 114 steps here, 112 in B2.
* Confidence-sequence: 114 here, 111 in B2.
* Known-baseline false-alarm rates at A: 1.0% per day and 12.6% over five days here, 0.8% and 11.6% in B1.

**Seeds.** All use numpy `default_rng` and lie inside 9,600,000 to 9,699,999. "law" is the law index (0 to 2) and "chunk" is a block of 500 paths.
* Truth samples: 9,600,000 to 9,600,002. The validation reused 9,600,000 to 9,600,003 with 1,000,000 draws.
* Tuning no-change paths: 9,610,000 + 100 x law + chunk, that is 9,610,000 to 9,610,039, 9,610,100 to 9,610,139 and 9,610,200 to 9,610,239.
* Fresh no-change paths: 9,620,000 + 100 x law + chunk, same pattern.
* Change cells: 9,630,000 + 1,000 x law + 100 x factor + chunk. The factor index is 0 to 3 for x2, x1.5, x0.5, x2/3, and chunks run 0 to 3.
* Validation paths: 9,690,000 to 9,690,003.
* Bootstrap of the median intervals: 9,699,999.

**Commands.** Run from `reports/round4/R2` with the venv active and `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`:
* `python r2_validate.py 40 > validate.log`
* `python r2_main.py 20000 20000 2000 <dir>/r2 > r2_main.log`
* `python r2_tables.py r2_rows.csv > r2_tables.md`
* `python r2_roc.py r2_rows.csv .`

**Files.**
* Scripts: `r2_lib.py` (the statistics), `r2_main.py`, `r2_validate.py`, `r2_tables.py`, `r2_roc.py`.
* `r2_rows.csv`: every threshold, law and cell, with Wilson intervals and bootstrap intervals for the medians.
* `r2_thresholds.json`.
* `r2_tables.md`: all tables, plus one more with every method at A.
* `roc_gaussian_var_x2.png`, `roc_student_t3_var_x2.png`, and the logs.

## 3. Results [SYNTHETIC]

### 3.1 Oracle thresholds and the false-alarm rates they really give

How to read this: each method's threshold was tuned on 20,000 no-change paths of the same law and checked on 20,000 fresh paths. All rates within 390 steps are near 5%, so the delay tables below compare methods at the same real rate. The last column shows that the e-detectors' rate keeps rising after the first day. [VERIFIED]

| Law | Method | Oracle threshold | Proved/nominal threshold | False alarms within 390 steps [95% CI] | Within 1,950 steps |
|---|---|---|---|---|---|
| Gaussian | README (plug-in, band 0.25) | 349.6 | 15,600 (= 2A) | 4.7% [4.4, 5.0] | 35.5% [34.9, 36.2] |
| Gaussian | Plug-in, band 0 | 6456 | 15,600 (= 2A) | 4.7% [4.4, 5.0] | 48.4% [47.7, 49.1] |
| Gaussian | Confidence-sequence baseline | 506.4 | 15,600 (= 2A) | 4.3% [4.1, 4.6] | 38.2% [37.6, 38.9] |
| Gaussian | Known baseline (theorem) | 4172 | 15,600 (= 2A) | 5.3% [5.0, 5.6] | 49.3% [48.6, 50.0] |
| Gaussian | Gaussian LR | 2.123e+04 | 15,600 (= 2A) | 4.9% [4.6, 5.2] | 25.3% [24.7, 25.9] |
| Gaussian | Page CUSUM (up only) | 32.91 | none | 5.4% [5.1, 5.7] | 20.9% [20.4, 21.5] |
| Gaussian | Repeated 60-step z-test | 3.773 | none | 5.2% [4.9, 5.5] | 20.8% [20.2, 21.4] |
| Student t3 | README (plug-in, band 0.25) | 471.1 | 15,600 (= 2A) | 4.9% [4.6, 5.2] | 30.0% [29.4, 30.7] |
| Student t3 | Plug-in, band 0 | 8212 | 15,600 (= 2A) | 4.9% [4.6, 5.2] | 42.4% [41.7, 43.0] |
| Student t3 | Confidence-sequence baseline | 662.1 | 15,600 (= 2A) | 4.8% [4.5, 5.1] | 29.7% [29.0, 30.3] |
| Student t3 | Known baseline (theorem) | 5781 | 15,600 (= 2A) | 4.8% [4.6, 5.1] | 38.0% [37.3, 38.7] |
| Student t3 | Gaussian LR | 2.83e+122 | 15,600 (= 2A) | 4.8% [4.5, 5.1] | 80.2% [79.6, 80.7] |
| Student t3 | Page CUSUM (up only) | 936.6 | none | 4.8% [4.5, 5.1] | 97.0% [96.7, 97.2] |
| Student t3 | Repeated 60-step z-test | 57.76 | none | 4.9% [4.6, 5.2] | 22.9% [22.3, 23.5] |
| Gaussian + jumps | README (plug-in, band 0.25) | 360.4 | 15,600 (= 2A) | 4.9% [4.6, 5.2] | 35.6% [34.9, 36.2] |
| Gaussian + jumps | Plug-in, band 0 | 6575 | 15,600 (= 2A) | 5.0% [4.7, 5.3] | 47.3% [46.6, 48.0] |
| Gaussian + jumps | Confidence-sequence baseline | 503.6 | 15,600 (= 2A) | 4.7% [4.4, 5.0] | 37.9% [37.2, 38.5] |
| Gaussian + jumps | Known baseline (theorem) | 4362 | 15,600 (= 2A) | 5.0% [4.7, 5.3] | 46.9% [46.2, 47.6] |
| Gaussian + jumps | Gaussian LR | 8.42e+58 | 15,600 (= 2A) | 5.1% [4.8, 5.4] | 26.8% [26.2, 27.5] |
| Gaussian + jumps | Page CUSUM (up only) | 412.8 | none | 5.1% [4.8, 5.4] | 38.6% [38.0, 39.3] |
| Gaussian + jumps | Repeated 60-step z-test | 33 | none | 5.0% [4.7, 5.3] | 23.4% [22.8, 24.0] |

### 3.2 Detection at an equal false-alarm rate (oracle tuning)

How to read this: each cell is the median delay in steps (IQR), then the share of 2,000 paths detected within five days. About 5% of paths alarm before the change in every row, so 95% is the most any row can reach. Bootstrap 95% intervals for the medians are in `r2_rows.csv`; on Gaussian data they are typically plus or minus 2 to 5 steps. [VERIFIED]

**Gaussian**

| Method | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|
| README (plug-in, band 0.25) | 37 (24-58), 95% | 99 (54-165), 96% | 73 (58-90), 95% | 134 (97-193), 95% |
| Plug-in, band 0 | 43 (29-63), 94% | 100 (59-151), 96% | 93 (75-112), 95% | 147 (112-188), 95% |
| Confidence-sequence baseline | 44 (28-68), 95% | 122 (66-212), 96% | 70 (54-86), 96% | 120 (84-169), 95% |
| Known baseline (theorem) | 41 (27-57), 95% | 86 (53-133), 95% | 84 (68-101), 95% | 127 (97-161), 95% |
| Gaussian LR | 43 (26-65), 94% | 116 (65-191), 96% | 66 (49-87), 95% | 144 (92-225), 95% |
| Page CUSUM (up only) | 38 (24-56), 94% | 92 (53-158), 95% | n/a | n/a |
| Repeated 60-step z-test | 42 (29-56), 94% | 92 (54-182), 95% | 523 (248-976), 75% | 990 (505-1428), 7% |

**Student t3**

| Method | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|
| README (plug-in, band 0.25) | 80 (48-130), 95% | 226 (116-416), 95% | 122 (80-182), 95% | 296 (168-513), 95% |
| Plug-in, band 0 | 79 (51-119), 96% | 189 (107-301), 95% | 122 (86-170), 95% | 240 (152-366), 95% |
| Confidence-sequence baseline | 91 (55-147), 96% | 270 (133-484), 94% | 116 (75-173), 95% | 262 (151-454), 95% |
| Known baseline (theorem) | 74 (48-106), 94% | 171 (102-259), 95% | 113 (81-155), 95% | 218 (140-321), 95% |
| Gaussian LR | 249 (167-343), 95% | 422 (274-583), 95% | 936 (464-1555), 7% | 1092 (673-1568), 17% |
| Page CUSUM (up only) | 203 (134-269), 94% | 305 (200-429), 95% | n/a | n/a |
| Repeated 60-step z-test | 750 (341-1228), 67% | 856 (391-1343), 42% | 1015 (580-1569), 8% | 1026 (532-1515), 10% |

**Gaussian plus jumps**

| Method | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|
| README (plug-in, band 0.25) | 39 (24-61), 95% | 102 (58-174), 95% | 77 (60-100), 95% | 150 (103-225), 94% |
| Plug-in, band 0 | 45 (30-65), 95% | 103 (63-157), 95% | 99 (77-124), 95% | 160 (116-219), 94% |
| Confidence-sequence baseline | 46 (28-72), 94% | 130 (71-224), 95% | 72 (54-94), 95% | 128 (88-191), 95% |
| Known baseline (theorem) | 42 (28-60), 95% | 92 (57-139), 94% | 89 (69-110), 94% | 140 (103-186), 94% |
| Gaussian LR | 554 (384-736), 95% | 1079 (662-1439), 78% | 1202 (668-1583), 32% | 936 (414-1384), 22% |
| Page CUSUM (up only) | 354 (234-469), 95% | 671 (415-977), 92% | n/a | n/a |
| Repeated 60-step z-test | 888 (444-1352), 35% | 886 (425-1402), 28% | 878 (373-1358), 19% | 939 (418-1397), 20% |

### 3.3 Where the README setting's slowness comes from

How to read this: cells are median delays in steps. The columns are:
1. The README setting at its proved threshold.
2. The same statistic at its oracle 5% threshold.
3. and 4. The same without the band, at each of the two thresholds.
5. The known baseline at its oracle threshold.
6. "Best non-e competitor": the fastest of the likelihood ratio, Page CUSUM and z-test at their oracle thresholds.

Going from column 1 to column 2 removes the whole gap to the best competitor on Gaussian data. Going from column 2 to column 4 shows that the band itself is not the problem. [VERIFIED]

| Law | Variance | README at 2A | README at oracle | Band 0 at 2A | Band 0 at oracle | Known baseline at oracle | Best non-e competitor at oracle |
|---|---|---|---|---|---|---|---|
| Gaussian | x2 | 114 | 37 | 51 | 43 | 41 | 38 (Page CUSUM) |
| Gaussian | x1.5 | 515 | 99 | 119 | 100 | 86 | 92 (Page CUSUM) |
| Gaussian | x0.5 | 265 | 73 | 113 | 93 | 84 | 66 (Gaussian LR) |
| Gaussian | x2/3 | 860 | 134 | 179 | 147 | 127 | 144 (Gaussian LR) |
| Student t3 | x2 | 273 | 80 | 90 | 79 | 74 | 203 (Page CUSUM) |
| Student t3 | x1.5 | 966 | 226 | 221 | 189 | 171 | 305 (Page CUSUM) |
| Student t3 | x0.5 | 518 | 122 | 138 | 122 | 113 | 936 (Gaussian LR, 7% detected) |
| Student t3 | x2/3 | 1114 | 296 | 284 | 240 | 218 | 1026 (z-test, 10% detected) |
| Gaussian + jumps | x2 | 118 | 39 | 52 | 45 | 42 | 354 (Page CUSUM) |
| Gaussian + jumps | x1.5 | 559 | 102 | 123 | 103 | 92 | 671 (Page CUSUM) |
| Gaussian + jumps | x0.5 | 308 | 77 | 119 | 99 | 89 | 878 (z-test, 19% detected) |
| Gaussian + jumps | x2/3 | 979 | 150 | 202 | 160 | 140 | 936 (Gaussian LR, 22% detected) |

### 3.4 Practical tuning: what a user without a guarantee would do

How to read this: the classical rules use thresholds tuned to 5% on Gaussian data. The e-detectors keep A = 7,800. A "detected" share near 0% means the rule had almost always alarmed falsely before the change. [VERIFIED]

| Law | Method | Threshold used | False alarms within 390 steps [95% CI] | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|---|---|---|
| Student t3 | README (plug-in, band 0.25) | 2A (proved) | 0.0% [0.0, 0.0] | 273 (170-421), 100% | 966 (572-1367), 45% | 518 (328-776), 98% | 1114 (710-1551), 19% |
| Student t3 | Plug-in, band 0 | 2A | 2.6% [2.4, 2.9] | 90 (59-133), 98% | 221 (132-353), 97% | 138 (101-193), 97% | 284 (182-435), 97% |
| Student t3 | Confidence-sequence baseline | 2A (proved) | 0.0% [0.0, 0.0] | 254 (155-388), 100% | 891 (537-1291), 53% | 347 (229-530), 100% | 927 (573-1302), 54% |
| Student t3 | Known baseline (theorem) | 2A (proved) | 1.6% [1.5, 1.8] | 91 (61-128), 99% | 221 (140-335), 98% | 141 (104-188), 98% | 286 (189-417), 98% |
| Student t3 | Gaussian LR | tuned on Gaussian | 99.9% [99.9, 100.0] | 0.1% detected | 0.1% detected | 0.1% detected | 0.1% detected |
| Student t3 | Page CUSUM (up only) | tuned on Gaussian | 99.9% [99.9, 100.0] | 0.0% detected | 0.1% detected | n/a | n/a |
| Student t3 | Repeated 60-step z-test | tuned on Gaussian | 99.9% [99.9, 100.0] | 0.1% detected | 0.1% detected | 0.1% detected | 0.1% detected |
| Gaussian + jumps | README (plug-in, band 0.25) | 2A (proved) | 0.0% [0.0, 0.0] | 118 (77-172), 100% | 559 (336-900), 93% | 308 (240-403), 100% | 979 (652-1374), 64% |
| Gaussian + jumps | Plug-in, band 0 | 2A | 2.0% [1.9, 2.3] | 52 (36-75), 98% | 123 (77-188), 98% | 119 (96-149), 98% | 202 (145-271), 98% |
| Gaussian + jumps | Confidence-sequence baseline | 2A (proved) | 0.0% [0.0, 0.0] | 119 (76-176), 100% | 556 (323-908), 91% | 220 (172-285), 100% | 562 (375-874), 95% |
| Gaussian + jumps | Known baseline (theorem) | 2A (proved) | 1.0% [0.9, 1.2] | 53 (37-75), 99% | 125 (80-183), 99% | 119 (97-147), 99% | 198 (150-257), 99% |
| Gaussian + jumps | Gaussian LR | tuned on Gaussian | 66.6% [66.0, 67.3] | 38 (22-60), 32% | 88 (47-144), 33% | 60 (41-81), 33% | 110 (63-178), 34% |
| Gaussian + jumps | Page CUSUM (up only) | tuned on Gaussian | 65.4% [64.7, 66.1] | 34 (21-51), 33% | 75 (42-126), 35% | n/a | n/a |
| Gaussian + jumps | Repeated 60-step z-test | tuned on Gaussian | 62.6% [62.0, 63.3] | 39 (25-52), 36% | 71 (43-135), 37% | 260 (111-517), 36% | 359 (139-719), 37% |

Extra: the e-detectors also tuned on Gaussian simulations, so no longer at the proved threshold. How to read this: the target was 5% within 390 steps. On Student t3 every variant overshoots it. [VERIFIED]

| Law | Method | False alarms within 390 steps [95% CI] | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|---|---|
| Student t3 | README (plug-in, band 0.25) | 9.8% [9.4, 10.2] | 68 (38-108), 90% | 170 (86-308), 89% | 97 (63-142), 90% | 216 (123-361), 90% |
| Student t3 | Plug-in, band 0 | 6.3% [6.0, 6.6] | 75 (48-112), 94% | 178 (101-283), 94% | 117 (81-160), 93% | 224 (139-341), 94% |
| Student t3 | Confidence-sequence baseline | 8.3% [7.9, 8.7] | 79 (46-127), 92% | 213 (112-397), 91% | 98 (64-146), 91% | 214 (122-356), 92% |
| Student t3 | Known baseline (theorem) | 7.1% [6.7, 7.4] | 69 (44-99), 92% | 153 (89-232), 93% | 102 (73-141), 93% | 194 (124-289), 92% |
| Gaussian + jumps | README (plug-in, band 0.25) | 5.4% [5.1, 5.7] | 39 (23-60), 95% | 100 (56-172), 95% | 75 (58-99), 94% | 144 (98-216), 94% |
| Gaussian + jumps | Plug-in, band 0 | 5.1% [4.8, 5.5] | 45 (29-65), 95% | 103 (62-156), 95% | 98 (77-124), 95% | 159 (116-218), 94% |
| Gaussian + jumps | Confidence-sequence baseline | 4.7% [4.4, 5.0] | 46 (29-72), 94% | 130 (71-225), 95% | 72 (54-95), 95% | 129 (89-194), 95% |
| Gaussian + jumps | Known baseline (theorem) | 5.3% [5.0, 5.6] | 42 (27-59), 94% | 91 (57-137), 94% | 88 (68-109), 94% | 138 (101-185), 94% |

### 3.5 Matched over five days instead of one

How to read this: each threshold now gives 25% false alarms within 1,950 steps on the same law. That is the five-day bound of A = 7,800. Cells are median delay (IQR). Under this matching the e-detectors use less of their first-day budget (1.6% to 4.0%) and become somewhat slower on Gaussian data. The ranking on heavy tails does not change. [VERIFIED]

| Law | Method | False alarms within 390 / 1,950 steps | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|---|---|
| Gaussian | README (plug-in, band 0.25) | 3.4% / 25.2% | 40 (25-61) | 108 (58-180) | 80 (65-99) | 153 (110-219) |
| Gaussian | Plug-in, band 0 | 1.6% / 24.9% | 52 (36-73) | 123 (77-188) | 116 (97-135) | 183 (144-233) |
| Gaussian | Confidence-sequence baseline | 2.6% / 25.0% | 49 (30-74) | 140 (77-245) | 78 (61-96) | 136 (97-195) |
| Gaussian | Known baseline (theorem) | 2.2% / 25.8% | 46 (32-64) | 105 (66-154) | 99 (82-116) | 151 (121-191) |
| Gaussian | Gaussian LR | 5.0% / 25.4% | 43 (26-65) | 115 (65-191) | 65 (49-86) | 144 (92-224) |
| Gaussian | Page CUSUM (up only) | 6.7% / 25.3% | 36 (23-54) | 85 (48-146) | n/a | n/a |
| Gaussian | Repeated 60-step z-test | 6.4% / 25.0% | 41 (28-55) | 84 (52-164) | 356 (160-694) | 885 (405-1358) |
| Student t3 | README (plug-in, band 0.25) | 4.0% / 25.2% | 85 (50-136) | 240 (123-446) | 130 (85-192) | 323 (179-554) |
| Student t3 | Plug-in, band 0 | 2.2% / 25.1% | 93 (61-137) | 230 (138-368) | 144 (105-201) | 294 (191-453) |
| Student t3 | Confidence-sequence baseline | 3.9% / 25.1% | 96 (58-153) | 288 (147-524) | 123 (79-185) | 278 (160-489) |
| Student t3 | Known baseline (theorem) | 3.0% / 25.2% | 81 (53-118) | 193 (118-296) | 125 (92-169) | 249 (163-364) |
| Student t3 | Gaussian LR | 1.2% / 25.0% | 535 (398-689) | 929 (653-1175) | 2.5% detected | 4.7% detected |
| Student t3 | Page CUSUM (up only) | 0.6% / 24.9% | 662 (530-805) | 1011 (788-1214) | n/a | n/a |
| Student t3 | Repeated 60-step z-test | 5.4% / 24.9% | 695 (322-1196) | 840 (396-1332) | 1000 (624-1546) | 1015 (504-1511) |
| Gaussian + jumps | README (plug-in, band 0.25) | 3.4% / 24.7% | 42 (26-66) | 112 (63-188) | 85 (66-111) | 176 (118-259) |
| Gaussian + jumps | Plug-in, band 0 | 1.8% / 25.6% | 53 (36-76) | 126 (78-191) | 121 (97-151) | 206 (150-276) |
| Gaussian + jumps | Confidence-sequence baseline | 2.9% / 24.8% | 50 (31-80) | 147 (84-262) | 81 (61-106) | 150 (104-225) |
| Gaussian + jumps | Known baseline (theorem) | 2.3% / 25.2% | 47 (32-68) | 108 (68-162) | 104 (83-129) | 170 (129-223) |
| Gaussian + jumps | Gaussian LR | 4.8% / 25.5% | 570 (397-754) | 1102 (679-1457) | 1199 (654-1564) | 932 (414-1390) |
| Gaussian + jumps | Page CUSUM (up only) | 2.3% / 25.2% | 470 (330-605) | 841 (560-1219) | n/a | n/a |
| Gaussian + jumps | Repeated 60-step z-test | 5.5% / 25.5% | 888 (440-1366) | 888 (422-1392) | 860 (369-1344) | 940 (423-1378) |

### 3.6 Delay against false-alarm rate (variance x2)

The two figures are `roc_gaussian_var_x2.png` and `roc_student_t3_var_x2.png`. Each line is one method at nine thresholds. Their target rates on the tuning paths run from 0.1% to 30%, and each point is placed at the rate measured on the fresh paths. Stars mark the methods at A = 7,800. A text box lists the stars that are off the scale.

How to read this: lower and further left is better.
* **Gaussian:** all seven lines lie in one band (at 1%: 51 to 58 steps; at 0.1%: 66 to 85 steps). The README and confidence-sequence stars sit far to the left of the 5% line, at 114 steps.
* **Student t3:** the e-detector lines lie below the classical ones at every rate. The gap is about 2.2 to 3.4 times at 5% (against Page and the likelihood ratio) and 5 to 13 times at 0.1%. [VERIFIED]

## 4. README claims that need changing

| File, line | What it says | What this run shows | Severity |
|---|---|---|---|
| README.md line 34 (E21 bullet) | "Ville's threshold costs little delay at a matched rate, while thresholds that use the full 5% are faster still" | True for the known baseline and band 0, whose rate at A is 1.0% to 1.9%. For the README's own setting and for the confidence-sequence baseline, A uses almost none of the 5% budget (0 and 1 false alarms in 20,000 days). It costs 3.1 to 6.4 times (README) and 2.6 to 4.1 times (confidence sequence) the delay of the same statistic at a 5% threshold. The README should say so, with these numbers. | Medium (incomplete; a reader would think the recommended setting loses little) |
| README.md line 104 (usage example) | `null_band=0.25` with `threshold=7800`, labelled only "plug-in baseline (approximate)" | This setting raised 0 false alarms in 20,000 independent no-change days. At the same A it is 2.2 times slower than band 0 for x2 (114 against 51 steps) and 4.8 times slower for x2/3 (860 against 179), and it misses 22% of x2/3 falls within five days. Add a comment that it is very conservative and slow, and why it is used (B1: band 0 exceeds 25% over five days). | Medium |
| README.md line 31 | The confidence-sequence baseline restores a proof "at a cost of 1.8 to 7.6 times the delay" | Consistent: here it is 1.7 to 4.2 times band 0 at A, with 1,000 calibration steps. But at an equal measured rate the confidence-sequence statistic needs 0.75 to 1.2 times band 0's delay (44 against 43, 122 against 100, 70 against 93, 120 against 147). The cost is in the threshold, not in the baseline. Worth one clause. | Low |
| docs/limitations.md line 15 | "This costs detection speed" (conservative bound) | Quantify: about 1.3 to 1.4 times for the known baseline (52 against 41 steps for x2) and 3 to 6 times for the README setting. | Low |
| README.md lines 26 and 34 (new evidence, no claim broken) | The advantage "comes from the statistic and the design" | Supported, and it can be strengthened. Even tuned on the true heavy-tailed law, the classical rules are 2.2 to 14 times slower for x2 (Student t3: Page 203 and likelihood ratio 249 against 74 to 91 steps; jumps: 354 and 554 against 39 to 46). They also miss most variance falls. | Low (addition) |

No README claim was broken. The most important correction: the README's recommended setting uses almost none of its false-alarm budget, and that is the whole reason it is slow.

## 5. What this does not show

* Only independent returns. Under GARCH clustering or unremoved seasonality the null fails for every method (B1), and equal-rate tuning there was not attempted.
* Oracle thresholds are not available in practice, because they need the true law. They measure the statistic, not a usable procedure. Section 3.4 is the usable comparison.
* One change time (step 400, just after the first day) and four variance factors. B2 found delays insensitive to the change time at the nominal thresholds. That was not rechecked at the oracle thresholds.
* Matching on the first day favours the e-detectors over five days, because their rate keeps rising (section 3.1). Section 3.5 gives the five-day matching. Other matchings, such as average run length or the rate on later days, could rank the methods differently.
* Page CUSUM was run only as an upward rule, as in B2. A two-sided CUSUM was not tested for decreases.
* The rates at A for the README and confidence-sequence settings (0 and 1 in 20,000) are too small to tune the competitors to. The smallest matched rate is 0.1%.
* Median delays count only detected paths. With about 95% detected in every Gaussian oracle row, this is a fair comparison. It is not fair for rows with low detection.
* No real data were used.

Note from the coordinator: the agent's file-writing tool refused to create this report, so the coordinator saved the agent's text here unchanged.
