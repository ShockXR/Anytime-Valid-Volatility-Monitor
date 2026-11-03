# TH6 iteration notes: false alarms after the first day (referee M2) and the hypothesis of Theorem D (referee M4)

Owner: TH6. Files: `src/avmonitor/periodic.py`, `tests/test_periodic.py`, `experiments/th6_per_period.py`, `results/th6_*.csv`, `results/figures/th6_*.png`, `docs/theory_per_period.md`, `docs/math.md` (remark after Theorem B, the bullet after it, and section 6.3; see iterations 6 and 8 for small consistency edits in 6.5 and 6.7), `docs/report_parts/TH6.md`, this file.

---

## Iteration 1. What Theorem B bounds after the first day

**OBJECTIVE.** Decide whether "at most 5% per day" follows from Theorem B or Corollary C, and if not, what does.

**PROBLEM IDENTIFIED.** The proof of Theorem B starts the supermartingale `M_t + (m - t)` at time 0 with value m. At the start of day d the Shiryaev-Roberts statistic has expected value up to `(d - 1) m`, so the same argument gives only `d m / A` for day d. The paper and the README read the theorem as a per-day bound.

**HYPOTHESIS.** No bound below 1 on the day-d probability holds for the plain detector over the null class, and an i.i.d. law should show it.

**CHANGES IMPLEMENTED.** Proposition 1 of `docs/theory_per_period.md`: the conditional bound `(N_s + m) / A` for any window, and the law `|r| = √m₀` (Y = 1), for which every Shiryaev-Roberts statistic equals t and the two-sided rule alarms at step 2A exactly (day 40 for both settings). A remark after Theorem B in `docs/math.md` and a corrected bullet ("5% chance of any false alarm per trading day" now reads "within the first trading day after each (re)calibration", with the long-run statement of Corollary C).

**EXPERIMENTS PERFORMED.** Unit test `test_constant_null_law_alarms_at_2A` (added in iteration 5).

**RESULTS.** Proved. Corollary C is unaffected: it bounds the long-run count (one alarm in 40 days is 0.025 per day).

**IMPROVEMENTS.** The reading of Theorem B is now stated where the theorem is.

**REMAINING ISSUES.** How large the effect is for natural streams (iteration 4).

**NEXT ITERATION.** Constructions that give a bound for every day.

---

## Iteration 2. Two constructions with a bound for every period

**OBJECTIVE.** A detector with a proved bound for every day, under as large a null class as possible.

**PROBLEM IDENTIFIED.** The Shiryaev-Roberts statistic carries everything it has accumulated. Two ways to stop that: restart it at day boundaries, or use the CUSUM-type statistic `C_t = f_t max(C_{t-1}, 1)`, which already exists in `detectors.betting_edetector(kind="CUSUM")` and is the e-CUSUM of the e-detector literature.

**HYPOTHESIS.** (a) Restarting makes each day look like the first, so Theorem B applies to every day given the past. (b) For the CUSUM statistic, reversing time turns the maximum over start times into the running maximum of a wealth process, so under an exchangeable null `P(C_s ≥ x) ≤ 1 / x` at every s, and a window bound of about m / A follows.

**CHANGES IMPLEMENTED.** Theorem 1 (restart, with a carry κ: `(m + κ) / A` per period, given the past; also the expected count per period). Lemma 1 (time reversal), Lemma 2 (conditional window bound under martingale differences), Lemma 3 (an integral bound for Pareto-dominated variables), Theorem 2: one bet `m / A`; an average over bets `(m + log A + H(w)) / A`; the default two-sided rule `(m + log 2A + (log 6 + log 5) / 2) / A` = 5.146% (390, 7800) and 5.044% (1440, 28,800); the same with restarts after alarms; the one-sided hazard by Harris's inequality. Code: `restarted_sr_bound`, `window_bound_single`, `window_bound_mixture`, `window_bound_two_sided`, `threshold_for_window` (A = 8027.69 and 29,053.41 for exactly 5%), `conditional_window_bound`, `reversed_running_max`; `PeriodicVolatilityMonitor(period=, offset=, carry=, statistic="sr" or "cusum")`, which with its defaults delegates to `VolatilityMonitor.update` unchanged.

**EXPERIMENTS PERFORMED.** Unit tests of the closed forms, of the time-reversal identity on random factor sequences, of Lemma 1's tail at a late time (40,000 paths), of Theorem 2(a) on a late window, of Theorem 1 on the calm stream, and of the streaming monitor.

**RESULTS.** All tests pass. **Failure recorded:** the first version of `window_bound_two_sided` divided the sum of the two one-sided bounds by 2 and returned 2.57% for the default rule. The union over the two sides, each at level 2A, gives the sum of the two one-sided bounds at 2A, not half of it; re-deriving it caught the error before any result used it, and the test now checks the closed form `(390 + log 15600 + (log 6 + log 5) / 2) / 7800`. **Second failure:** the test comparing the CUSUM and Shiryaev-Roberts streaming monitors failed because the plug-in monitor raised an alarm and recalibrated, after which the two monitors had different baselines; the test now compares them only before the first event and requires more than 100 checked steps.

**IMPROVEMENTS.** Two proved per-period guarantees; a monitor class that implements the restart without changing any default.

**REMAINING ISSUES.** Theorem 2 needs exchangeability; whether it survives martingale differences.

**NEXT ITERATION.** Martingale-difference nulls.

---

## Iteration 3. Martingale-difference nulls: the counterexample

**OBJECTIVE.** Decide whether Theorem 2 holds under the null class of Theorem B (conditional mean fixed, any dependence).

**PROBLEM IDENTIFIED.** The time reversal needs the law of the factor sequence to be invariant under reversal. An adaptive stream can steer the statistic.

**HYPOTHESIS.** Because Y = 1 leaves `max(C, 1)` unchanged, a stream can park the statistic just below the threshold for as long as it likes and then strike.

**CHANGES IMPLEMENTED.** Proposition 3 (hunt with `{0, B}` until one draw could reach the threshold, park with Y = 1, strike): no alarm before the strike, the parking time is finite almost surely, and the window probability is at least `(1 - o(1)) / B`, about 11%, for every m and A. Code: `park_strike_policy`, `bold_policy`, `two_point`.

**EXPERIMENTS PERFORMED.** Unit test at small scale (m = 10, A = 60); exploration at full scale with 4000 paths (not reported, superseded by the 30,000-path run of iteration 4).

**RESULTS.** At small scale the window probability is above the exchangeable bound by more than five standard errors, and below the conditional bound of Lemma 2. **Failure recorded:** a first version of the full-scale part computed B through the default arguments of `default_bets` instead of the law and failed in the smoke run; fixed to `law_B("gaussian")`.

**IMPROVEMENTS.** A sharp division: the restart covers martingale differences; the e-CUSUM's window bound needs exchangeability.

**REMAINING ISSUES.** How close a martingale-difference stream can come to the restarted detector's bound after day 1 (bold play reaches 2.75%).

**NEXT ITERATION.** The full simulations.

---

## Iteration 4. Simulations

**OBJECTIVE.** The day-by-day numbers the referee asked for, with standard errors below 0.2 points where they matter.

**CHANGES IMPLEMENTED.** `experiments/th6_per_period.py`, parts `laws`, `hazard390`, `hazard1440`, `family`, `adversary`, `delay`, `renewal`, `tolerance`, `figures`. A vectorised engine (`periodic.run_engine`) that runs several detectors on the same draws, checked to reproduce `betting_edetector` plus `two_sided_alarm` exactly for both recursions (about 3.8 ms per step for 3 detectors on 30,000 paths). The true τ and m0 of the four E9 laws are computed by numerical integration (`results/th6_laws.csv`).

**EXPERIMENTS PERFORMED.** All with `nice -n 10` and at most two processes.
* Run A (01:59 to 02:52 UTC, 3 November 2026): `laws hazard390 family adversary delay tolerance`. hazard390: 91 to 180 s per stream; family: about 2 minutes per law at 390 steps and 10 at 1440; bold play: 53 to 69 s (390) and 169 to 220 s (1440) per detector; park and strike: 103 s and 210 s; the delay part at 4000 paths per cell (its CSV was kept in the scratch directory and replaced by run C).
* Run B (01:59 to 02:49): `hazard1440`, 378 to 795 s per stream.
* Run C (02:52 to 03:02): `delay` at 20,000 paths per cell (the 4000-path run took 3 to 13 s per cell, so the precision could be raised cheaply), then `renewal`.
* Run D (02:53): `tolerance` again after iteration 6's correction.
* `figures` last.

**RESULTS.** See `docs/report_parts/TH6.md`. Plain detector, 390-step days: hazard 0.98% on day 1, 4.70% on average over days 11 to 20 (Gaussian); bold play 5.25% over days 11 to 20; the referee's family 28.0% on day 10 and 65.1% on day 20; Corollary C's count 0.037 to 0.048 per day. Restarted: at most 2.75% (390) and 3.52% (1440) on any day. e-CUSUM: at most 0.18% (i.i.d.) and 0.40% (bold play). Park and strike: above the i.i.d. bound from day 17, 6.16% on day 20. **Failure recorded:** for the e-CUSUM without restarts the engine's count is the number of steps at or above the threshold, not of alarms, so its mean count and long-run rate are meaningless; they are written as missing values, and only the probability of a day with the statistic above the threshold is reported. The seed list in the script's docstring was wrong in the first version and was corrected.

**IMPROVEMENTS.** Every number in the report part comes from these CSV files.

**REMAINING ISSUES.** The pooled bold-play hazard over days 11 to 20 was chosen after seeing the table (post hoc).

**NEXT ITERATION.** Delay.

---

## Iteration 5. Delay, and why the restart costs what it costs

**OBJECTIVE.** Compare delays at the same proved level for every day, and explain the cost of the restart for changes at random times within a day.

**PROBLEM IDENTIFIED.** In the first delay run the restarted detector's cost depended on where in the day the change fell, and the plain detector was faster than any construction with a per-day bound.

**HYPOTHESIS.** A renewal argument: after a restart the detector is a fresh start on i.i.d. post-change data, period after period.

**CHANGES IMPLEMENTED.** Proposition 2: `E[D_r] ≤ E[min(D₀, m - k)] + P(D₀ > m - k) E[min(D₀, m)] / P(D₀ ≤ m)` for a change k steps into a day, and `D ≤ D₀` for the plain detector. Part `renewal`, which simulates `D₀` over one day and averages the bound over the change position. This part was added after the delay results were seen (post hoc). Unit test for Proposition 1(b).

**RESULTS.** The restarted detector's mean delay was 1.06 to 1.44 times the plain detector's (390-step days) and 1.02 to 1.14 times (1440), below the bound of Proposition 2 in all 12 cells; the e-CUSUM at its certified A was 1.24 to 1.83 times slower. The plain detector is 14% to 31% faster than a fresh start because of the statistic it builds up under the null.

**IMPROVEMENTS.** The trade-off is explicit: the head start that makes the plain detector fast is the accumulation that raises its later-day false alarm probability.

**REMAINING ISSUES.** A rate-matched comparison (each detector at the threshold that gives the same measured per-day rate) was considered and not run: it would be law-specific and would drop the proved level that the comparison is about.

**NEXT ITERATION.** M4.

---

## Iteration 6. Theorem D with an algorithm-free hypothesis (M4)

**OBJECTIVE.** State Theorem D (Theorem 6.4 of `docs/math.md`) under a null that does not involve the bets, check the proof, and say what the hypothesis excludes.

**PROBLEM IDENTIFIED.** The previous hypotheses compared the monitored conditional mean with bet-weighted calibration averages, whose weights the algorithm chooses from the data.

**HYPOTHESIS.** Since a weighted average lies between the smallest and the largest of its terms, `sup μ_t ≤ (1 + b) min_calibration μ_i` implies the bet-weighted condition, and the proof goes through unchanged.

**CHANGES IMPLEMENTED.** In `docs/math.md` section 6.3: (H↑) and (H↓) as tolerances with a horizon h (h = m for the bound within m steps, h = ∞ for the run length); Step 2 of the proof now passes through the min and max; a part 4 with no assumption (`δ + P(V_m) + m / A`, by stopping the statistic at the predictable time at which the tolerance first fails); a paragraph on what the tolerance excludes; the bet-weighted version kept as a remark; a paragraph on Shekhar and Ramdas (2024). For consistency, Corollaries 6.5 and 6.6, the paragraph "Choosing A and δ" and three sentences of 6.5 that described the old hypothesis were adjusted.

**EXPERIMENTS PERFORMED.** Part `tolerance`: GARCH(1,1) in the four cells of E9 and a slow log-normal stochastic volatility, calibration window and one monitored day of 390 or 1440 steps, `μ_t` exact, 2000 paths per cell.

**RESULTS.** The band needed has a median of 0.90 to 5.01 under GARCH and 0.58 to 1.22 under the stochastic volatility; no path of any cell satisfied the tolerance with b = 0.25. **Failure recorded:** in the first run the GARCH returns used for τ were drawn afresh given the volatility path instead of being the returns that drive the recursion. This only changes τ slightly (the medians moved by at most 0.031), but the part now uses the driving returns, and the CSV is from the corrected run.

**IMPROVEMENTS.** The null of Theorem D is now a condition on the data alone (for a fixed calibration length).

**REMAINING ISSUES.** For the stochastic volatility, `μ_t` is computed given the volatility path; given the returns alone it is smoother and the band needed is smaller.

**NEXT ITERATION.** References.

---

## Iteration 7. References

**OBJECTIVE.** Relate the guarantees to Lai (1998), Mei (2008), Pergamenchtchikov and Tartakovsky (2018) and Shekhar and Ramdas (2024), and check every reference.

**PROBLEM IDENTIFIED.** arxiv.org and the publishers' pages could not be opened from this machine (every direct fetch failed), so only web searches were possible.

**CHANGES IMPLEMENTED.** Section 6 and the verification note of `docs/theory_per_period.md`; the Shekhar and Ramdas paragraph in `docs/math.md`.

**RESULTS.** Verified by web search: bibliographic data of all four papers and of Harris (1960); Lai's criterion `sup_k P_∞(k ≤ T < k + m)` through a secondary source; Mei's result on mixture pre-change models and the alternative criteria named in his keywords; the local probability of false alarm of Pergamenchtchikov and Tartakovsky from their abstract; the construction and the run-length guarantee of Shekhar and Ramdas from their abstract. Not verified: the exact conditions on Lai's window, Mei's exact definitions, whether the local probability of Pergamenchtchikov and Tartakovsky conditions on no earlier alarm, and the full text of Shekhar and Ramdas.

**REMAINING ISSUES.** A reading of the full texts when they can be opened.

**NEXT ITERATION.** The independent review of round 8 (iteration 8).

---

## Iteration 8. Independent review, round 8

**OBJECTIVE.** Apply `docs/review/independent_review_round8_th6.md` (no Critical, 2 Major, 12 Minor) and the consistency items raised for the paper.

**PROBLEM IDENTIFIED.** (M1) The proved per-day bound was credited to `PeriodicVolatilityMonitor`, which had only the plug-in and confidence-sequence baselines, and Theorem 1 assumes a known one; the extension to the confidence-sequence baseline was asserted, not proved, and stated as if it held given the past. (M2) The summary said Theorem 2 holds "whatever happened before"; its bounds are unconditional, and given the state at the start of a window the next-step alarm chance can be about `P(Y > 1)`, 32% for Gaussian returns. Minor: the sharpness of `d m / A` was overstated (m1); a remedy in 6.5 works against the restated hypothesis (m2); 6.7 listed parts 1 to 3 of Theorem 6.4 (m3); a wrong range of daily hazards, mixing the 390 and 1440 tables (m4); part 4 needs the missingness assumption of 6.7 and `T < ∞` (m5); Proposition 2 left three conditions unstated (m6); the GARCH bands are per path while parts 1 to 3 need the tolerance almost surely (m7); wording and small numbers (m8 to m12, including an off-by-one in the descriptive parked share). From the paper: 2.75% is a hazard, not a probability (the largest probability is 2.66%); the README-29 row of the appendix had misread the awake-B figure of E14, which is 0.035 false alarms per day as a long-run rate; paper theorem numbers drop out of these files.

**HYPOTHESIS.** A known-baseline mode of the class, built on the existing calibration hook, reproduces the simulation engine step by step, so Theorem 1 covers it as stated. For the confidence-sequence baseline the proof of Theorem 6.4, Steps 1 to 3, applies on each period of a segment conditionally on `F_b`, which gives `1{Gᶜ} + 1{G} (m + κ) / A` and, after taking expectations, `δ + (m + κ) / A`, marginally only.

**CHANGES IMPLEMENTED.**
* `src/avmonitor/periodic.py`: `PeriodicVolatilityMonitor(baseline="known", tau=..., m0=..., tau_down=None, m0_down=None)`. It validates `tau > 0` and `0 < m0 ≤ tau²` on both sides, refuses `margin_z` and `target_rel_width`, sets the null values directly in `_reset` (so it never calibrates and restarts at once after an alarm), and refuses `tau` or `m0` without `baseline="known"`. The class and module docstrings now say which baseline carries which bound, and that Theorem 2's bound is unconditional.
* `tests/test_periodic.py`: three tests. The known mode equals the engine event by event (period 50 and no period, with a change in the stream); on a calm-then-Gaussian stream it raises no alarm while calm and its per-day rate is within three standard errors of `m / A`, while the plain detector's is more than five above; the errors above are raised.
* `docs/theory_per_period.md`: Corollary 1 (confidence-sequence baseline) with its proof and three differences from Theorem 1 (marginal only; tolerance to the end of the period; per segment); Theorem 1 and Theorem 2 state "known baseline"; Theorem 2 says its bounds are unconditional and that given `F_s` only Lemma 2 holds; the conditions of Proposition 2 (m6); the non-aligned window bound `(2m + κ) / A`; m1, m8, m9, m10, m12.
* `docs/math.md`: remark after Theorem B (m1, m4, m11, the two constructions with their baselines); "What the tolerance excludes" (m7); part 4 (m5); "Choosing A and δ" (the per-day statement for the confidence-sequence baseline); 6.5 (m2); 6.7 (m3).
* `docs/report_parts/TH6.md`: the same corrections, the 2.66% probability against the 2.75% hazard, the day-16 crossing for 1440-step days, and the README-29 appendix row.
* `experiments/th6_per_period.py`: `parked_before_day` uses `<=`; the park-and-strike runs are now their own part, `parkstrike`, so they can be rerun without bold play.

**EXPERIMENTS PERFORMED.** `python experiments/th6_per_period.py parkstrike` (same seeds) to apply the `<=` fix; the full test suite. Numbers read off existing CSVs for m1 (`th6_calm_family.csv`), m4 (`th6_hazard_390.csv`), m12 (`th6_adversary.csv`, `th6_park_strike.csv`).

**RESULTS** [SYNTHETIC].
* The known mode matches the engine exactly in both tests. **Failure recorded:** the first sizing of the calm-stream test (19 calm periods) made the plain detector alarm with certainty at step 2A, and nine calm periods with 800 paths left it too close to `m / A` to separate; 15 calm periods are used.
* Rerun of park and strike: every alarm probability is identical; the parked shares moved by at most 0.007 points (2 paths in 30,000) at full scale and 0.07 points in the scaled run, so no rounded figure changed (12.1% on day 20 at 390, 57% by period 120 in the scaled run).
* Changed numbers in the text: Gaussian SR hazard on days 10 to 20 at 390 steps from "4.4% to 4.9%" to "4.70% on average over days 11 to 20, single days 4.34% to 5.18%"; the sharpness of `d m / A`: at most 35.7% against 50% (d = 10) and 73.2% against 100% (d = 20) over both settings and both laws of the last day; pooled bold-play hazard over days 11 to 20 5.23% with standard error about 0.06 points (was "about 0.05"); park and strike passes the i.i.d. bound on day 16 at 1440 (5.30%, 2.0 standard errors); Mei paragraph 4.7% (was 4.6%); a window that straddles one boundary `(2m + κ) / A` (was `2 (m + κ) / A`).

**IMPROVEMENTS.** Every per-day claim now names its baseline: known (given the past, Theorem 1, implemented and tested), confidence sequence (marginal, Corollary 1), plug-in (none).

**REMAINING ISSUES.** Corollary 1 is marginal: no bound is proved for a period given the observed past with an estimated baseline. It also inherits the tolerance of Theorem 6.4, which the GARCH cells of the `tolerance` part never satisfy for b < 1. The known mode needs `m0` and `tau` from outside the stream, which in practice is rarely available.

**NEXT ITERATION.** None planned. Open questions are listed in section 7 of the theory note.
