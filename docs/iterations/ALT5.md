# ALT5 iterations: two-sample betting monitor (conformal rank e-detector)

Prototype folder: `alternatives/two_sample_betting/`. Times are UTC, night of 2 November 2026.

---

ITERATION NUMBER: ALT5-1 (formulation and design pilot)

OBJECTIVE: Turn entry A4 of `alternatives/ALTERNATIVES.md` (a sequential two-sample test by betting between the recent window and the calibration window) into a working change detector with an exact false alarm guarantee and no estimated baseline.

PROBLEM IDENTIFIED: The paired two-sample betting test of Shekhar and Ramdas (2024) needs a fresh reference draw per step, but the monitor has one finite calibration window. Two obvious adaptations fail: drawing reference points without replacement runs out after `n0` steps and doubles the noise of every payoff; ranking new returns against the fixed reference (quantile bins learned on the calibration window) gives p-values that are uniform one at a time but dependent, so a learning bettor is not betting fairly.

HYPOTHESIS: Ranking each new return among ALL earlier returns of the segment (calibration window plus monitored returns, a growing reference) gives smoothed conformal p-values that are exactly i.i.d. uniform under exchangeability (Vovk, Nouretdinov and Gammerman 2003). Any bet that is fair against the uniform law and uses only past p-values then gives a Shiryaev-Roberts e-detector with exactly the guarantees of `docs/math.md` 5.4. Mixing a rank analogue of the main monitor's truncated-variance bet (Gaussianised score, null mean known to be 1) with a histogram learner (plug-in density bettor) should cost little on variance changes and gain a lot on shape changes.

CHANGES IMPLEMENTED: `rankbet.py` (vectorised exact p-values with a binary indexed tree, streaming `SequentialRanker`, Gaussianised score, three bettor families, null band `b`, Shiryaev-Roberts mixture, streaming `RankMonitor`), `tsb_sim.py` (unit-variance laws, GARCH, a vectorised copy of the plug-in `VolatilityMonitor`), `run_experiment.py` part `pilot`, `tests/test_rankbet.py`.

EXPERIMENTS PERFORMED: Pilot on development seeds 1000 to 1002 (never reused): seven configurations (memory of the histogram learners 60/240, 240/1440, 1440/5760 steps; prior 50 or 200; family weights; scale bets only; histogram only) on four changes (200 paths) and two one-day nulls (400 paths). Output `results/p_pilot_design.csv`.

RESULTS: [SYNTHETIC] Scale bets alone were fastest for pure variance changes (median delay 149.5 against 167.5 steps for the default mixture for sd x1.25) but slow for Gaussian to t4 at fixed variance (1348, 74% detected); the histogram learner alone was 1.5 to 2 times slower on variance changes but fast on t4 (387). The mixtures C0 (default), C1, C3 and C4 were within about 10% of each other on every row. Null false alarms within one day: 1.0% to 4.0% (bound 5%). Under GARCH(1,1) with persistence 0.95 every configuration with scale bets alarmed on 70% to 77% of days; the histogram learner alone on 25%.

IMPROVEMENTS: C0 (scale up 0.25, scale down 0.25, histogram 0.5; half-lives 240 and 1440; prior 50) was frozen as the default before any reported run: it was never the worst mixture and never more than 8% behind the best mixture on any change.

REMAINING ISSUES: The pilot cannot tell whether real data are close enough to exchangeable for the exact null to be useful.

NEXT ITERATION: Unit tests of the proofs' ingredients, then the reported synthetic study.

---

ITERATION NUMBER: ALT5-2 (tests of the validity ingredients; two failed first versions)

OBJECTIVE: Check numerically every step of the validity proof and the equivalence of the benchmark with the repo's monitor.

PROBLEM IDENTIFIED: Two first test versions failed. (1) A test asserted that a Gaussian to t3 change (same variance) would be carried by the histogram learner; the alarm was carried by the scale-down bets (92% of `M_t`). (2) A test of the band null checked `E_h[f] <= 1` on `h = 1 + b sign(w)`, which is not a probability density, and failed.

HYPOTHESIS: (1) is not a bug: a variance-matched t3 puts more mass near zero, so the Gaussianised truncated variance falls (the same effect the quarticity prototype reports for the truncated variance), and the decrease bets pick it up. (2) is a wrong test: the band proof only needs `E_h[w] <= b E_U|w|` for densities `h`.

CHANGES IMPLEMENTED: (1) The test now checks detection soon after the change and that the family shares add up to 1. (2) The test uses the most favourable valid density `h = 1 + b sign(w - median w)`, renormalised.

EXPERIMENTS PERFORMED: 9 tests, about 3 seconds: p-values against brute force with ties, streaming against vectorised, uniformity and zero lag-one correlation for i.i.d. Poisson data, `E_U[Y] = 1`, positivity and fairness of every bet for a learner state driven far from uniform (bands 0 and 0.3), `E[M_t] = t` under the null, the vectorised plug-in monitor against `VolatilityMonitor` alarm for alarm, detection of a tail change, NaN handling and configuration checks.

RESULTS: All pass. The p-values agree with brute force exactly, and for i.i.d. Poisson scores (heavy ties) the mean, variance, 5% quantile and lag-one correlation of the p-values are within the test tolerances (0.005, 0.002, 0.004 and 0.01) of the uniform values.

IMPROVEMENTS: Each line of the proof has a test.

REMAINING ISSUES: None for validity under exchangeability.

NEXT ITERATION: Q1 to Q3, the synthetic comparison with the mean-based monitor.

---

ITERATION NUMBER: ALT5-3 (synthetic validity and power against the mean-based monitor)

OBJECTIVE: Measure, on the same simulated paths, the false alarm rates and detection delays of the frozen rank detector and of the repo's plug-in `VolatilityMonitor` (Q1, Q1b, Q2, Q3 of `run_experiment.py`).

PROBLEM IDENTIFIED: The plug-in monitor's guarantee is proved for a known baseline; its validity tests in the main report (E8) cover short horizons only. A fair comparison needs both detectors on identical paths, at the same `A`, and validity checked over several days as well as one.

HYPOTHESIS: (1) The rank detector holds both the one-day (5%) and seven-day (35%) bounds for every return law and reference length. (2) It is a little slower for pure variance changes (the cost of an exact null), clearly faster for tail changes at fixed variance (histogram learners), and faster for a change of jump intensity at fixed variance. (3) Truncating a rank score changes nothing that matters.

CHANGES IMPLEMENTED: Parts `q1` (three laws, references of 1440 and 120 steps, 1000 paths, 7 days), `q1b` (GARCH(1,1), persistence 0.95 and 0.99, 500 days), `q2` (ten changes, 400 paths, change at monitored step 500), `q3` (score `|u|` against `min(|u|, 3)`). The plug-in monitor runs through the vectorised copy in `tsb_sim.py`, checked alarm for alarm against the class by a test.

EXPERIMENTS PERFORMED: As above, fixed seeds per part between 100 and 599 (listed in the script), never the pilot seeds. Outputs `results/q1_validity.csv`, `q1b_garch.csv`, `q2_power.csv`, `q3_truncation.csv`, `fig1_synthetic_delays.png`.

RESULTS: [SYNTHETIC] (1) Confirmed: rank detectors 1.9% to 3.5% per day and 22.2% to 24.3% per week. Unexpected: the plug-in monitor at band 0 holds the one-day bound with a one-day calibration (4.6% to 6.0%, at the bound) but alarms within seven days on 51.4% to 57.0% of paths (bound 35%), and on 42.0% to 45.7% of days with a two-hour calibration; its baseline error acts as a constant drift. Band 0.25 silences both designs on i.i.d. data. (2) Partly confirmed: rank scale bets are 4% to 7% slower on variance changes (136 against 127 steps for sd x1.25), the omnibus mixture 11% to 19%; Gaussian to t4 at fixed variance is 2.9 times faster (417 against 1203). Refuted for jumps (a recorded negative result): the plug-in monitor is faster for 0.2% to 0.5% jumps per step (1393 against 1688, 88% against 54% detected), because holding total variance fixed lowers the diffusive variance and the truncated variance falls. Only the signed score sees a pure skew change (809 steps, 99% detected); for `|u|`-based detectors the "detections" there are false alarms after the change point. (3) Confirmed: delays move by at most about 6%, false alarm rates by noise only. Clustering (Q1b) breaks both designs about alike (rank b 0.25: 9.8% and 53.0% of days; plug-in band 0.25: 15.4% and 67.0%).

IMPROVEMENTS: The seven-day horizon is now part of every validity table. The plug-in's multi-day failure is reported as a finding for the main monitor.

REMAINING ISSUES: Is the rank detector's validity edge simply a matter of how the reference is used? Q1c tests a fixed-reference rank design. The comparison with the repo's proved variant (confidence-sequence baseline) is missing.

NEXT ITERATION: Q1c and Q5.

---

ITERATION NUMBER: ALT5-4 (fixed reference ablation and the confidence-sequence baseline)

OBJECTIVE: Show why the reference must grow, and compare with the main monitor's rigorous variant instead of only its plug-in form.

PROBLEM IDENTIFIED: A reviewer would ask whether ranking against the calibration window only (quantile bins learned once, a simpler design) works as well, and whether the rank detector's speed advantage over a valid mean-based detector survives against `VolatilityMonitor(baseline="cs")`.

HYPOTHESIS: The fixed-reference design is not a martingale: the p-values share one reference and are dependent, a learning bettor gains about `(K - 1) / (2 n0)` in log wealth per step, so false alarms become certain. The confidence-sequence monitor is valid but conservative and therefore slower.

CHANGES IMPLEMENTED: Part `q1c` with `pvalues_fixed_reference` (same bettors, same paths, references of 240 and 1440 steps, 14 days, 500 paths), part `q5` (the class with `baseline="cs"`, `delta = 0.01`, on the Q1 Gaussian paths and five Q2 scenarios).

EXPERIMENTS PERFORMED: Seeds 150, 151 (Q1c) and the Q1 and Q2 seeds (Q5). Outputs `results/q1c_fixed_reference.csv`, `q5_cs_baseline.csv`.

RESULTS: [SYNTHETIC] Fixed reference: 90.4% to 93.4% of first days with a false alarm for a 240-step reference, 21.2% to 29.8% for a one-day reference, 98.4% to 100% within a week. Growing reference: below the bound at 1, 7 and 14 days (43.2% to 50.2% at 14 days against 70%). Confidence-sequence monitor: 0.0% per day and 0.1% per week; median delays 348 (sd x1.25), 400 (sd x0.8), 1129.5 with 92% detected (jumps 0.2% to 1%), and the t4 change detected on only 10% of paths; the rank omnibus detector is 1.7 to 2.8 times faster on every row.

IMPROVEMENTS: The fixed-reference failure and the plug-in failure are now shown side by side as one phenomenon (betting against a quantity estimated once).

REMAINING ISSUES: Real data.

NEXT ITERATION: Q4 on BTCUSDT and ETHUSDT.

---

ITERATION NUMBER: ALT5-5 (real data: BTCUSDT and ETHUSDT, July to December 2025; one confound found and removed)

OBJECTIVE: Compare alarm counts and overlap with `VolatilityMonitor` at band 0.25 and describe the alarms raised only by the two-sample test.

PROBLEM IDENTIFIED: (1) The streaming `RankMonitor` was several times slower per step than needed because the truncated moment factor of the Gaussianised score (scipy quadrature) was recomputed at every step; fixed with a cached `_kappa2`, numerics unchanged. (2) The first version of the comparison counted every rank alarm without a main-monitor alarm within 60 minutes as "rank only": 61 for BTC at `b = 0.25`. Both monitors are blind for one calibration day after each alarm and monitor only about a third of the time, so most of these could not have been matched.

HYPOTHESIS: Most single-design alarms fall while the partner is calibrating; the genuine "rank only" alarms (partner watching) should look like shape changes at constant scale.

CHANGES IMPLEMENTED: Monitored intervals of every monitor; each unmatched alarm split by whether the partner was monitoring; chance level of matching from 200 random shifts of 1 to 30 days; 2000 random control windows with the same reference delays; per-alarm learned density ratio of new to reference p-values; weekend flag; a nonzero-returns-only variant; a 120-minute-reference variant. All detector settings were fixed before the run.

EXPERIMENTS PERFORMED: Part `q4`, outputs `results/q4_real_*.csv`, `fig2_real_alarm_profiles.png`, `fig3_real_density_ratio.png`.

RESULTS: [REAL] At band 0.25 the main monitor raises 127 (BTC) and 122 (ETH) alarms, the rank detector 123 and 114; 50% and 61% of rank alarms are matched within an hour (chance 5.7%). Of the 61 unmatched BTC rank alarms, 44 fell while the main monitor was calibrating (ETH: 40 of 45). The remaining 17 BTC alarms have a robust scale ratio of 0.994 and a shape distance of 0.131 against 0.076 at random times; the learned density ratio is 0.32 for the largest 0.2% of moves and 1.49 to 1.52 for the smallest moves; 18% fall on weekends against 30% at random times; dropping zero returns halves the group to 8. The main monitor alone raised 32 alarms per coin while the rank detector watched, all decreases (truncated variance ratio 0.50 and 0.52). With a 120-minute reference the rank detector monitors 88% to 91% of the time.

IMPROVEMENTS: The headline "only the two-sample test" count is now 17 and 5 instead of 61 and 45, with the blind-time artefact reported separately.

REMAINING ISSUES: Real returns cluster, so neither design's null holds; the 17 alarms are descriptive evidence, not tested discoveries.

NEXT ITERATION: Q6, the cost of the short reference that removes most blind time.

---

ITERATION NUMBER: ALT5-6 (short references; one rerun after a design error)

OBJECTIVE: Measure how much slower the rank detector is with a 120 or 240-step reference, since its guarantee does not depend on the reference length.

PROBLEM IDENTIFIED: The first version of part `q6` drew the reference before the monitored data from one generator, so the monitored paths differed across reference lengths and the delays were not paired. It was stopped before finishing and its output discarded.

HYPOTHESIS: The rank detector stays valid with a short reference and pays a moderate delay; the plug-in monitor with the same window is not valid.

CHANGES IMPLEMENTED: Monitored data drawn first from the same seed for every reference length, then the reference; parts `q6` and `timing` (cost per step, `results/t_timing.csv`).

EXPERIMENTS PERFORMED: Four Q2 changes, 400 paths, references 120, 240 and 1440 steps, seeds 600 to 603. Output `results/q6_short_reference.csv`.

RESULTS: [SYNTHETIC] Rank omnibus detector: at most 0.8% false alarms in the 500 pre-change steps (bound 1.7%) at every reference length; 17% to 52% slower with 120 steps than with 1440 (for example 193 against 153 steps for sd x1.25; 623 against 409, 88% detected, for t4). Plug-in monitor with 120 steps: 11.5% to 14.0% false alarms before the change. Cost: about 80 microseconds per step for `RankMonitor` against about 10 for `VolatilityMonitor`.

IMPROVEMENTS: The README and RESULTS state when a short reference is worth it (less blind time against slower detection).

REMAINING ISSUES: The conjectured asymptotic efficiency of 1 for the Gaussianised score under Gaussian data is not tested; ranks relative to a volatility forecast (to handle clustering) are not built.

NEXT ITERATION: None in this session.
