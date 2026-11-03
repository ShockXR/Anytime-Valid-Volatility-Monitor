# Methodology validation: anytime-valid volatility monitor (prototype B)

Written 3 November 2026 by validation agent 2 (an AI agent), at commit d2e41b6 plus the working-tree edits listed in `reports/findings/A2_B.csv`. Labels: REAL (market data), SYNTHETIC (simulated), PUBLISHED (printed in a paper), OUR CHECK (recomputed by this validation).

## Summary for a beginner

1. The project builds an alarm that watches one-minute price returns and shouts when volatility goes up or down, with a mathematically proved cap on how often it shouts falsely, however often you look.
2. It borrows one idea from Prof. Bapat's 2026 paper (clip each return at a multiple of its median absolute deviation, MAD) and puts it inside a betting-based sequential test (an "e-detector").
3. I checked 16 load-bearing formulas and rules against the code (by reading and by 20 independent numeric tests): all 20 tests pass (OUR CHECK, `reports/checks/equation_vs_code_check.py`). The guarantee holds on a fresh simulation (0.9% false alarms in the first day against a 5% bound) and the "day 40" counterexample in the docs reproduces exactly.
4. The proof is for the first day after each calibration only, not for every later day, and the real-data event result (FOMC announcements) failed on fresh data. The pitch already says both; I made sure the README and pitch also say that every review so far was done by AI agents and fixed a stale test count.
5. Open items are in section 7 and in `reports/findings/A2_B.csv`.

## 1. Problem

Detect, online and with a provable limit on false alarms, a rise or fall in the variance of high-frequency returns that contain jumps and fat tails. Classical detectors fire on jumps; a test checked every minute fires by chance eventually. Neither prototype is a trading system; metrics are false-alarm rate against the proved bound, detection delay, precision against simple baselines, calibration and reproducibility.

## 2. Assumptions

| Assumption | Where | Why it matters |
| --- | --- | --- |
| Null: the conditional mean of X_t = min(r_t^2, tau^2) stays at or below its baseline (rise) or at or above (fall) while "nothing changes". | `docs/math.md` 5.1 | Everything proved needs it. Volatility clustering, price-level noise and estimated baselines break it (E9, E10, E14, ALT1). |
| Baseline m0 is known, or comes from a confidence sequence. A plug-in baseline is only approximate. | `docs/math.md` section 6 | Plug-in gave 41% to 44% false alarms over five days against a 25% bound (E10, E14, SYNTHETIC). |
| Truncation level tau, bets and intraday profile are predictable (use only past data). | `detectors.py`, `data.intraday_scale` | A bet that peeks ahead voids the supermartingale argument. |
| An alarm bound m/A is "5% per day" only for the first day after a (re)calibration. | `docs/math.md` Theorem B and remark; `docs/theory_per_period.md` | The README and pitch used to be read as 5% every day (referee comment M2). |

## 3. Methodology (pipeline)

1. Calibrate on `cal_len` returns (default 1,440): scale = median|r| / 0.6745, tau = c x scale with c = 3, baseline m0 = mean of min(r^2, tau^2) (`detectors.calibrate`, `monitor.py`).
2. Each new return gives X_t = min(r_t^2, tau^2) and Y_t = X_t / m0. Increase side: bet factor 1 + lambda (Y_t - 1), lambda in [0, 1). Decrease side: 1 + lambda (1 - Y_t), lambda < 1 / (B - 1) with B = tau^2 / m0.
3. Shiryaev-Roberts statistic M_t = f_t (M_{t-1} + 1), averaged over a grid of bets (increase 0.02 to 0.5; decrease 0.1 to 0.95 of the largest legal bet). Alarm when either side reaches 2A. A = 7,800 for 390-step days (5% per day), 28,800 for 1,440-minute crypto days.
4. On an alarm, recalibrate with fresh data (the default monitor is therefore blind for about one calibration window; roughly 70% of the time on real crypto data, E12).
5. Baseline options: plug-in; plug-in with a margin or null band; confidence-sequence (tau from the first 25% of the window, bounds from the rest at delta = 0.01; proved bound delta + m/A).

Hyperparameters and their source: c = 3 (paper's convention, with the MAD scaled by 1.4826; the paper's tables fit the raw MAD, `docs/estimator_issues.md`), grid and A as above, cal_len 1,440, null band 0 by default (0.25 in E4c and the README example), delta 0.01.

## 4. Implementation: equation versus code

| # | Claim (document) | Code | Verdict | Evidence |
| --- | --- | --- | --- | --- |
| 1 | SR recursion M_t = f_t (M_{t-1} + 1) equals the sum over start times of products of bets (`math.md` 5.4) | `detectors.py:272` | MATCH | Brute-force sum equals recursion to 1e-10 (OUR CHECK, 60 steps) |
| 2 | CUSUM recursion C_t = f_t max(C_{t-1}, 1) = max over start times; C <= M | `detectors.py:274` | MATCH | Brute force equal; C < M |
| 3 | Increase factor 1 + lambda (Y - 1), lambda in [0, 1) | `detectors.py:221-227, 265` | MATCH | `betting_edetector` path equals hand recursion |
| 4 | Decrease factor 1 + lambda (1 - Y), lambda < 1 / (B - 1), B = tau^2 / m0 | `detectors.py:229-240` | MATCH, with a doc gap | Equal to hand recursion. For direction "down" a user-supplied grid is read as fractions of the largest legal bet (the code multiplies it by lam_max), while the docstring says "sequence of bets" (finding A2-B-005) |
| 5 | X = min(r^2, tau^2); tau = c x scaled MAD about zero; scaled MAD = median|r| / Phi^-1(3/4) | `detectors.py:100-102` | MATCH | tau equals 3 x median|r| / 0.6745 exactly; m0 equals mean of X |
| 6 | Surviving share of 2nd and 4th moments after clipping: kappa2(a), kappa4(a); table in `math.md` 3 (raw MAD c = 3: a = 2.023, 0.925, 0.705; scaled: 3.000, 0.995, 0.964) | `robust.py:199-214` | MATCH | Closed form equals 4 million Monte Carlo draws to 0.001 and the printed table to 3 decimals |
| 7 | Effective cut-off in standard deviations: raw MAD a = 0.6745 c; scaled a = c | `robust.py:217` | MATCH | 2.0235 and 3 |
| 8 | Two-sided: each side at 2A keeps the budget | `monitor.py` (module docstring, update loop), `detectors.py:288` | MATCH | Fresh simulation, Gaussian returns, known baseline, 6,000 paths x 390 steps: 0.9% false alarms in the first day against the 5% bound (OUR CHECK, SYNTHETIC) |
| 9 | Theorem B: P(false alarm within m) <= m / A via S_t = M_t + (m - t) and Ville's inequality | proof in `math.md` 5.4 | MATCH (proof read, simulation consistent) | I read the argument and found no gap; I did not formally verify it. AI review round 1 reported no error in Theorems A and B |
| 10 | TH6 Proposition 1b: law Y_t = 1 makes every SR statistic equal t, so the two-sided rule alarms with certainty at step 2A = 15,600, the last step of day 40 | `detectors.py` | MATCH | First alarm at step 15,600 = 40 x 390 (OUR CHECK) |
| 11 | Confidence-sequence baseline: upper and lower bounds each at delta/2 cover the true clipped mean | `confseq.py:287-` | MATCH (conservative) | Coverage 100.00% over 4,000 simulated windows for both bounds against a 99.5% target (SYNTHETIC); the bounds are generous, which is the price of a time-uniform proof |
| 12 | Decrease speed limit: no valid SR detector can alarm in fewer than ceil(log(1 + A/B) / log(B/(B-1))) steps; B = 9.045 at c = 3 | `design.py:411` | MATCH | 58 steps at A = 7,800 and 64 at 2A; default bets 73 and 79 (OUR CHECK; pitch quotes 58, 64, 73, 79) |
| 13 | Adaptive bet uses local moments of data strictly before t | `detectors.py:140-165` | MATCH (by reading) | m1[:, L:] = acc[:, :n-L] indexes only earlier steps |
| 14 | Missing values leave the statistic unchanged (no bet, no new chain) | `detectors.py:261-275` | MATCH (by reading) | |
| 15 | Calibration window rules: cal_len 1,440; degenerate windows extended up to cal_max = 2 cal_len; cs split tau_frac 0.25; delta 0.01 | `monitor.py:110-130`, `detectors.py:56-120` | MATCH | Signature equals module docstring |
| 16 | Intraday profile = median of |r| at that clock minute over previous 20 days, current day excluded | `data.py:78-118` | MATCH, with a legacy option | Default start="past" is predictable. start="pooled" uses up to 20 days of future data and is kept to reproduce E5, E7, E8; E7's report part says so |
| 17 | Equation map and paper estimators (Part 1 of `docs/equation_map.md`) | `robust.py` | NOT RE-CHECKED | Covered by tests and AI review round 1; outside the 20 formulas I sampled |

## 5. Experimental evidence (all numbers are in the committed reports and traced by `scripts/trace_numbers.py`)

* **Guarantee in simulation (SYNTHETIC).** Robust e-detector 0.3% to 2.3% false alarms in the first day against 5% under jumps and fat tails, where Gaussian likelihood-ratio and CUSUM detectors alarm on 70% to 97% of paths; median delay 36 to 49 steps for a 50% volatility rise (E4). Measured 1.1% to 1.4% for four return laws (E9). The worst legal adversary found reached 2.62% in the first day and 18.1% over five days against 25% (E16).
* **Where it breaks (SYNTHETIC).** Estimated plug-in baseline: 41% to 44% over five days against 25%. Price-level noise: 8.5% to 41.5% per day. Clustering: 82% to 100% of quiet days. Each has a documented fix that costs speed.
* **Real data (REAL).** On permuted calm days of five coins (real ticks, zeros and fat tails, exchangeable): first-day rate 2.26% (known baseline), 0.01% (confidence sequence), 3.38% (rank detector) against 5%, plug-in 41% over five days (E19, pre-registered). Block-resampled calm days alarm on 86% to 100% of streams, so the minute-level null is rarely true. Alarms beat random times by 21 to 38 points after multiplicity adjustment but add little to a strong HAR forecast (E20).
* **Event result.** FOMC statements on BTC: 8 of 16 against 10.9% of controls, Holm p = 0.0003 on already-seen 2024 to 2025 data (E22). Pre-registered on 16 unseen 2026 releases it failed: 2 of 16 against 5.9%, p = 0.25, planned power 0.82 (E23). A guard against blind stretches did not rescue it on seen data (4 of 16, p 0.056, post hoc, E24). A forward test to 1 June 2027 and an anytime-valid version (E25) are registered. The FOMC result is exploratory.
* **Live run (REAL).** Nifty session of 3 November 2026: no alarm, because the monitor was recalibrating after an alarm the day before; foreseeable, not evidence.
* **Tests and tracing (OUR CHECK).** 581 tests pass (baseline log). Number trace rerun after my edits: README 282 numbers, pitch 199, paper 2,552, all traced. The trace can match coincidentally (see section 7).

## 6. Results in one table

| Question | Answer | Strength |
| --- | --- | --- |
| Does the proved bound hold inside its conditions? | Yes, with a wide margin (about 4 times conservative with default bets) | Strong (proof, simulation, adversarial search) |
| Does it hold for every day? | No, first day after each calibration and long-run average only; restart every day gives every day with a known baseline | Proved (TH6) |
| Does it hold with an estimated baseline? | Only with the confidence-sequence baseline (delta + m/A); plug-in is approximate | Proved for cs baseline |
| Are real-data alarms informative? | Beat random times; add little to a strong forecast | Moderate |
| Does the FOMC result stand? | No: failed on fresh data | Exploratory only |
| Is it better than classical CUSUM? | No for precision once thresholds are calibrated; the edge is a calibration-free proved threshold | Stated in E21 |

## 7. Limitations and failure scenarios

* No human has reviewed the proofs or code; the 11 review rounds and 3 referee reports were AI agents (now stated in the README, pitch, strategy file and each review header).
* The null is about conditional means; volatility clustering, microstructure noise, coarse ticks and drifting baselines violate it. The default monitor is blind after each alarm.
* Data are public mirrors, not exchange feeds; the paper's own data source is not named.
* The trace script cannot tell a right match from a coincidence: before my edit "About 570 tests" counted as traced because 570 appears in E25 as an alarm count (finding A2-B-003).
* The committed trace CSVs under `results/trace/` were stale (README trace dated from the round 2 text); rerunning the script regenerated them (finding A2-B-004).
* README states that E18 onward was not rerun on a fresh clone; Agent 1's reproducibility report decides whether that is still true.
* Failure scenarios worth rehearsing: a stream with clustering and a plug-in baseline (false alarms far above 5%); a long run past day 40 with a calm exchangeable law; a calibration window with more than half zero returns; the monitor blind through the event of interest (47 days on BTC with the cs baseline).
