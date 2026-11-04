# V1. Independent check of R2 and R4 (round 4)

Checker V1, 4 November 2026. I did not edit anything in R2, R4, src, tests, results or docs. My scripts and logs are in this folder.
Seeds for my own simulated paths: 9,960,000 to 9,960,013. All runs single core.

## R2 (detection speed at equal false-alarm rates, SYNTHETIC)

What I checked in the code:
* **Running maximum shortcut.** Every method in the package alarms the first time its statistic reaches the threshold, so "first alarm" equals "first time the running maximum passes the threshold". The README monitor and the Gaussian LR monitor alarm when either side's mean reaches 2A; r2_lib.py uses max(up, down) against h = 2A, which is the same rule. Restarts after an alarm do not matter, because R2 only ever uses the first alarm (a false alarm within 390 steps is the first alarm; a delay is counted only when the first alarm comes after the change). The band (m_up = 1.25 m0, m_dn = 0.75 m0), the decrease-bet sizes, the truncation and the 1,000-step calibration match monitor.py line by line. The cap at 1e300 cannot change any threshold used (the largest is 2.8e122).
* **Independent paths.** Tuning seeds 9,610,000 + 100 law + chunk, checking seeds 9,620,000 + ..., change cells 9,630,000 + ...; with 40 chunks of 500 these never overlap.
* **Delay.** The first changed return is monitored step 400 (sig[:, nu:] in simulate.py). Alarms before step 400 are excluded and counted as "early"; delay = alarm - 400 + 1, the same rule as detectors.summarize_alarms.

Cells recomputed by streaming the package's own classes, with my own random number generation:

| Cell | R2 | V1 |
|---|---|---|
| README, Gaussian, null rate at oracle 349.6 | 4.7% | 4.85% (2,000 paths, +-0.9) |
| README, Gaussian, variance x2 | 37 (24-58), 95% | 37 (22-59), 94% (1,000 paths) |
| Page CUSUM, t3, null rate at oracle 936.6 | 4.8% | 4.9% (8,000 paths, +-0.5) |
| Page CUSUM, t3, variance x2 | 203 (134-269), 94% | 204 (138-276), 94% (2,000 paths) |
| Gaussian LR tuned on Gaussian, run on t3: null rate | 99.9% | 100% (666 of 666) |
| same, variance x2 detected | 0.1% | 0% |
| Gaussian LR, Gaussian null at 21,233 | 4.9% | 4.4% (+-0.9) |
| Gaussian LR, t3, oracle 2.8e122: null rate / x2 | 4.8% / 249 (167-343), 95% | 5.75% (+-1.0) / 250 (158-349), 95% |
| README, t3, oracle 471.1: null rate / x2 | 4.9% / 80 (48-130), 95% | 5.25% (+-1.0) / 80 (48-132), 94% |

All agree within Monte Carlo error. The headline ratios in the R2 summary (3.1 to 6.4 times, 17% faster to 33% slower, 3.0 to 5.5 times, 63% to 67%) follow from the tables.

Minor points (none moves a headline):
* r2_validate.py compares with the package on only 40 paths per cell. My larger checks above cover this.
* Page CUSUM is scaled by the calibration MAD, not by the mean of r^2. I tried the mean instead (t3, variance x2, own 5% threshold): it is slower (339 steps against 204), so R2's choice does not handicap the CUSUM.
* For the jump law, "variance x2" multiplies only the diffusive part (total variance x1.76). Section 2 of the R2 report says so; the table headers do not.
* z-test thresholds above about 8 (t3: 57.8, jumps: 33) cannot be passed through the package's `level` argument (the normal quantile underflows); the shortcut is still the same rule.

**Verdict R2: CONFIRMED.**

## R4 (real data per asset and year)

What I checked:
* **Blind share and alarms come from the monitor's real state.** run_monitors.py stores, before each update, whether the monitor is still calibrating; the minute that completes a calibration window counts as blind and the alarm minute counts as watched, which is right. I rebuilt the inputs and reran VolatilityMonitor myself (plug-in, c = 3, band 0.25; crypto A = 28,800, 1,440-minute calibration; Nifty A = 7,500, 750 minutes), and derived the blind minutes a second way, from the monitor's segment starts and alarm indices. The returns, the blind flags and every alarm index are identical to R4's saved runs.
* **60-minute windows.** "Post" covers the 60 minutes after the alarm minute (alarm minute excluded), "pre" the 60 minutes ending with it, and the baseline is the segment's own calibration window. No look-ahead or off-by-one.
* **Placebo.** Same pool (watched minutes of the same asset and year with valid windows) and, for Nifty, the same minute of the session. Code is correct for its definition.

| Asset-year | R4 blind / up / down / up hit | V1 |
|---|---|---|
| BTC 2024 | 71.00% / 133 / 112 / 0.940 | identical |
| BTC 2025 | 69.42% / 138 / 116 / 0.928 | identical |
| BTC 2026 | 65.99% / 117 / 78 / 0.906 | identical |
| Nifty 2023 | 42.60% / 25 / 27 / 0.960 | identical |
| Nifty 2024 | 64.88% / 43 / 37 / 0.977 | identical |
| Nifty 2025 | 55.43% / 36 / 33 / 0.861 | identical |

The claim "64% to 72% on crypto (17 asset-years), 43% to 65% on Nifty (52% overall)" matches asset_year.csv (63.7% to 71.7%; 42.6% to 64.9%).

Minor points:
* The placebo does not match on recent volatility. Up alarms come when the last hour was about 2.8 times the baseline. Among watched minutes with similarly high recent volatility, the hit share is 79% for BTC 2025 (alarms 93%, R4 placebo 30%) and 98% for Nifty 2024 (alarms 98%). So "useful years" mostly says alarms come when volatility is already high; R4's own post / pre result (0 of 27) points the same way. (v1_r4_placebo_probe.py)
* summarise.py line 59 computes year weights `w` but never uses them: the 2026 reference values are plain means over years. The years have similar lengths, so this changes little.
* experiments/e07_nifty.py line 150 fills the first Nifty sessions of 2017 with the median of the first 20 sessions, a small look-ahead limited to early 2017.

**Verdict R4: CONFIRMED (minor interpretation caveat on the placebo).**

Files: v1_r2_check.py / .log, v1_r2_check2.py / .log, v1_r4_check.py / .log, v1_r4_placebo_probe.py / .log.
