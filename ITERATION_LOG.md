# Iteration log

This log records what was tried, what went wrong and what changed. It follows the ten-step iteration framework of the project brief, in its required format (objective, problem, hypothesis, change, experiments, results, improvements, remaining issues, next step). Each framework iteration below consolidates several concrete rounds of work. Failed attempts are part of the record.

* **Round 1** (2 November 2026, daytime): first prototype, experiments E1 to E5. Its original 16-row log is kept unchanged in the appendix at the end.
* **Round 2** (night of 2 to 3 November 2026): real data for all five coins of the paper and Nifty 50, an independent review, the confidence-sequence baseline, forecast-relative monitoring, new theory and alternative prototypes (experiments E6 to E17).

Detailed notes per work stream, with every sub-iteration and failure, are in `docs/iterations/`: `MAIN.md` (integration and data), `E6.md` (five-coin reproduction), `E7.md` (Nifty 50 and India VIX), `E8.md` (confidence-sequence baseline), `ALT1.md` to `ALT5.md` (alternative prototypes), `EV1.md` (crypto event study), and the theory and benchmark tracks as they finish. The independent review is `docs/review/independent_review_round1.md`.

| Framework iteration | Primary objective | Main evidence |
|---|---|---|
| 1 | Initial implementation and baseline reproduction | E1, E5a, E6a, E6b |
| 2 | Mathematical correctness and equation verification | `docs/estimator_issues.md`, `docs/math.md` 5 and 6, review C1, E10 Q10d |
| 3 | Numerical stability and edge cases | `docs/iterations/MAIN.md` R2-2, E8 iteration 1, review M6 and m17 |
| 4 | Statistical validation and bias analysis | E2, E4d, E4e, E8 |
| 5 | Benchmarking against alternative methods | E4, E6e, E13, `alternatives/` |
| 6 | Sensitivity analysis and stress testing | E4c, E9, E10, ALT1 (price noise) |
| 7 | Larger-scale simulations and profiling | E9, E8 (10,000 paths per cell) |
| 8 | Architecture and implementation | `docs/iterations/MAIN.md` R2-3, `forecast.py`, `confseq.py` |
| 9 | Independent verification, regression testing, reproducibility | review round 1 and its fixes, CI, manifest, reruns |
| 10 | Final research audit and presentation readiness | clean-clone audit, one-command reproduction, review round 2 |

---

ITERATION NUMBER: 1 (initial mathematical implementation and baseline reproduction)

OBJECTIVE: Implement Bapat (2026) exactly as printed, with every benchmark of its Tables 3 and 4, and find out whether the published tables come from the printed formulas.

PROBLEM IDENTIFIED: The paper does not state the volatility used in its simulations, the exchange behind its real-data tables, or whether its kernel weights and MAD are normalised. Without those, "reproduction" could mean anything.

HYPOTHESIS: If the code follows the printed formulas, a single volatility value should reproduce Table 7, and the ratio Kernel / RQ in Table 6 should equal the raw weight sum divided by the sample length, whatever the data.

CHANGES IMPLEMENTED: `robust.py` (printed and normalised estimators, MinRQ, MedRQ, TPQ, QPQ, BPQ, TRQ), `simulate.py` (the paper's return model), `experiments/e01_reproduce_tables.py`, `e05_real_btc.py`; in round 2 `experiments/e06a_five_coins_tables.py` and `paper_tables.py` with all published rows.

EXPERIMENTS PERFORMED: E1 (fit one σ on 15 rows of Table 7, predict the other 75 numbers); E5a (Table 6 BTC column on Bitstamp); E6a and E6b (Tables 5 and 6 for all five coins on Binance, with five falsifiable hypotheses for every mismatch).

RESULTS: σ = 0.0215 reproduces the plain-kernel rows of Table 7 with a typical error of 3.5%, but the robust rows at σ_J = 0.10 come out about 5 times smaller than published (not reproduced; reported as a question). On real data Kernel / RQ equals W(30) / (n - 30) for every coin, so the code uses unnormalised weights and RQ without the factor n. The published ratios imply n between 262,854 and 263,052, about 2,000 returns fewer than a full half year, which fits a venue with missing minutes. The raw (unscaled) MAD fits the robust kernel (ratio 0.81 to 0.96) and the scaled MAD does not (1.86 to 2.48). Binance matches the published BTC and XRP quartiles within 0.4% and 1.7% and four of five means within 2.5%. The ETH mean has the opposite sign with the same size (a sign or ordering slip is consistent). DOGE's Q1 = Q3 = 0 needs a price grid of at least 0.00025 USD, and then the printed robust kernel must be exactly 0, not the published 2.00 × 10⁻¹². The real-data RMSE column is reproduced in 29 of 35 entries by |Q - 2.2 × 10⁻⁶|.

IMPROVEMENTS: A precise, evidence-based list of questions for the author instead of vague doubts. The first DOGE hypothesis (prices rounded to 4 decimals) failed and is recorded as failed.

REMAINING ISSUES: The σ_J = 0.10 robust rows and the DOGE row cannot be explained without the author's code.

NEXT ITERATION: Check the mathematics behind each formula.

---

ITERATION NUMBER: 2 (mathematical correctness and equation verification)

OBJECTIVE: Verify every formula used, in the paper and in our own code and proofs.

PROBLEM IDENTIFIED: Round 1 found that the MinRQ constant in our first implementation was wrong (estimates 1.247 times the truth), and that the paper's eq. (16) is algebraically a truncated RQ. Round 2's review found an omitted cross term in one of our own checks and a proof step that used the wrong condition.

HYPOTHESIS: Each claim can be checked by algebra plus a Monte Carlo test, and the checks will find errors on both sides.

CHANGES IMPLEMENTED: Corrected the MinRQ constant to π / (3π - 8) with a unit test for every constant; `docs/estimator_issues.md` (scale of the kernel estimator, non-vanishing truncation bias, eq. 16 equals truncated RQ, weight order in eq. 15); `docs/math.md` section 5 (test supermartingales, Theorems A and B, the impossibility of detecting a variance decrease nonparametrically); section 6 (Theorem 6.4: P(false alarm within m steps) ≤ δ + m / A with an estimated baseline); `docs/equation_map.md` (paper notation to code). After the review: the RQ bias check now includes the jump cross term 6σ²σ_J² (1.17 × 10⁻⁴, not 1.07 × 10⁻⁴); the proof of Theorem A uses Fatou's lemma; a restart corollary (Corollary C); the impossibility result is stated precisely and attributed (Bahadur and Savage 1956). In E10, a closed form ρ(k) = k(1 - β) / (1 - α - β + αk) for the share of a k-fold break that survives standardisation by a frozen GARCH forecast.

EXPERIMENTS PERFORMED: Monte Carlo checks of every constant; the identity of section 4 to machine precision; the restart corollary (E4e); ρ(k) against simulation (E10, Q10d).

RESULTS: All constants check out after the MinRQ fix. The reviewer confirmed sections 1 to 5 correct apart from precision issues. Corollary C holds numerically (26.5 false alarms per 10,000 steps with the true baseline and 30.1 with plug-in recalibration, against a bound of 33.3). ρ(k) matches simulation to three decimals at persistence 0.95 and 0.99.

IMPROVEMENTS: Every theorem in the docs now has a proof that a referee can check, with its assumptions stated (filtration, independence of the innovation, missing data).

REMAINING ISSUES: The proofs assume a conditional null that real data violate in known ways (iteration 6).

NEXT ITERATION: Edge cases that real data expose.

---

ITERATION NUMBER: 3 (numerical stability and edge-case handling)

OBJECTIVE: Make the estimators and monitor safe on coarse-tick, gappy and frozen data.

PROBLEM IDENTIFIED: With more than half of the returns exactly zero (coarse ticks as in the paper's DOGE, frozen prices during an outage), the MAD is 0, every return is clipped to 0 and estimators silently return 0; the monitor raised `ZeroDivisionError`. The online scale counted NaN as zero and could underflow to exactly 0. The review found that two scripts turned missing returns into zeros (which look like evidence for a volatility fall: a 120-minute gap triggered a false "down" alarm) and that one NaN silently disabled the vectorised detector.

HYPOTHESIS: Warnings, a tick-size floor, window extension and treating missing values as "no bet" fix these without changing any result on clean data.

CHANGES IMPLEMENTED: `robust.DegenerateThresholdWarning`, `robust.tick_floor`, a safe `ew_quantile_scale`; window extension and fallback in the monitor (E8, iteration 1); NaN passed through in E5, E8 and E10 and skipped in `betting_edetector`; `calibrate` drops NaN; a check that decrease bets stay admissible; `data.mark_frozen` for long runs of exact zeros. Earlier (round 1) the raw-r² and Gaussian benchmarks overflowed with large jumps; exponents are now capped.

EXPERIMENTS PERFORMED: `tests/test_edge_cases.py` (9 tests), regression tests that the safe online scale is bit-identical to the old rule on clean data, reruns of E5 and E8 real-data parts.

RESULTS: All tests pass. Reruns changed no real-data result (each file has one missing value). One test first failed for an instructive reason: with 60% zeros the true median of |r| is 0, so the tracker legitimately drifts toward 0 and can climb back by at most a factor e every 562 steps; recovery within a day needs a tick floor.

IMPROVEMENTS: No silent zeros and no crash path; missing data handled by a rule that keeps the guarantee (math.md 5.6, item 5).

REMAINING ISSUES: The tick floor needs the price tick, estimated from the data.

NEXT ITERATION: Validate the statistics and measure bias.

---

ITERATION NUMBER: 4 (statistical validation and bias analysis)

OBJECTIVE: Measure the bias of truncation and the effect of an estimated baseline on the guarantee.

PROBLEM IDENTIFIED: Theorems A and B assume the true null mean m₀. In practice it is estimated, and in round 1 a unit test with 390 calibration steps gave 67% false alarms against a 50% bound.

HYPOTHESIS: The estimation error explains the excess; a confidence-sequence baseline removes it with a proof.

CHANGES IMPLEMENTED: E2 (truncation bias); `margin_z` (round 1); `confseq.py` and `VolatilityMonitor(baseline="cs")` (E8); `experiments/e04e_review_checks.py` (plug-in against the true baseline).

EXPERIMENTS PERFORMED: E2a to E2d; E4d (calibration length); E4e (plug-in against true baseline, 4000 paths); E8 (10,000 paths per null cell, three laws, three calibration lengths, real data); E10 Q10b2 (five-day horizon).

RESULTS: Clipping at 3 raw MADs keeps 70.5% of the fourth moment and at 3 scaled MADs 96.4%, for every n. On real data truncation removes 52% to 85% of the quarticity measured by jump-robust estimators, because one-minute returns are fat-tailed (E6c). With 390 calibration steps the plug-in baseline roughly triples the one-day false alarm rate (4.6% to 5.1% against 1.6% to 1.8% with the true baseline); over five days it gave 44% against a guaranteed 25%. The confidence-sequence baseline kept the one-day rate at or below 0.27% against the proved 5%, at a cost of 1.8 to 7.6 times the plug-in delay with 390 calibration steps (first written as 7.2; corrected in review round 2a) and 1.2 to 1.7 times with 5000. On real data it removed 16% to 20% of the alarms, almost all decrease alarms, and raised the share of decrease alarms followed by lower volatility from 90% to 97% (band 0.25).

IMPROVEMENTS: A proved guarantee with an estimated baseline, and a quantitative account of what the default plug-in costs. A first unit test failure (67% against 50%) was the starting point.

REMAINING ISSUES: The cs baseline is very conservative with short calibration windows; a running truncation level is not yet covered (E14).

NEXT ITERATION: Compare with other methods.

---

ITERATION NUMBER: 5 (benchmarking against alternative methods)

OBJECTIVE: Compare the robust e-detector with parametric, classical and naive detectors, and with published change-point methods.

PROBLEM IDENTIFIED: A detector is only interesting if it beats a fair competitor, and on real data "alarms are followed by higher volatility" does not separate a detector from volatility clustering (review M3).

HYPOTHESIS: The robust detector is slower than parametric detectors on clean Gaussian data but the only valid one under jumps and fat tails; on real data its precision is similar to a well-tuned naive rule, and its value is the guarantee.

CHANGES IMPLEMENTED: E4 benchmarks (Gaussian likelihood ratio, Page CUSUM, repeated z-test, raw r²); the quarticity-adaptive bet; E6e (EWMA and HAR forecasts, rate-matched ratio rules); E13 (Bibinger, Jirak and Vetter 2017; Bayesian online change-point detection); alternative prototypes in `alternatives/`.

EXPERIMENTS PERFORMED: E4a and E4b (300 paths per cell); Q4c; E6e (five coins, in and out of sample); E13 (see its report part).

RESULTS: With jumps or Student t returns, Gaussian methods false alarm on 70% to 97% of days; the robust e-detector stays at 0.3% to 2.3% and detects a 50% increase in a median of 36 to 49 steps, about 10% to 70% slower than parametric methods on clean data. The adaptive bet was slightly slower than a fixed grid in 14 of 18 cells (negative result). On real data, alarms beat EWMA and HAR forecasts 71% to 83% of the time against about 50% at random times, but a ratio rule tuned with hindsight to the same number of alarms is about as precise (negative result for precision; the e-detector adds a guarantee with fixed settings). E13: Bayesian online change-point detection (BOCPD) false alarms on 73% (jumps) and 53% (Student t5) of days at a threshold giving 1% on Gaussian data, against 1.2% to 1.3% for the e-detector, but detects decreases two to three times faster on clean data (negative result for our speed on decreases). A rolling adjacent-window test in the style of Bibinger, Jirak and Vetter misses most small changes. On five coins at equal up and down alarm counts, the e-detector's alarms validate more often than all four competitors with 2 to 4 hour windows (78.5% and 80.5% against 58.7% to 71.7%; held out 82.0% and 69.7% against 56.7% to 71.5%, except truncated BOCPD on held-out down alarms, p = 0.47); on Nifty no difference is significant. Review round 3 then showed that the outcome (the next day against a one-day forecast) suits the e-detector's time scale: with day-long windows truncated BOCPD and the BJV-type tests match its down alarms (plain BOCPD stays lower) and the best competitor trails its up alarms by a borderline 6 points, and the BJV-type test's low offline power was a window-length effect (99% with 500-return windows).

IMPROVEMENTS: An honest statement of what the method adds, now against published methods as well as naive rules.

REMAINING ISSUES: The BJV implementation is a reconstruction not yet checked against the paper's equations (the full text was not reachable from our environment); BOCPD with a heavy-tailed likelihood was not tried; p-values on real data treat alarms as independent.

NEXT ITERATION: Stress the assumptions.

---

ITERATION NUMBER: 6 (sensitivity analysis and parameter stress testing)

OBJECTIVE: Find where the guarantee stops applying and what the knobs (band, calibration length, forecast) do.

PROBLEM IDENTIFIED: Round 1: with slow stochastic volatility every detector alarmed on about 25% of days. Round 2: under persistent GARCH clustering the raw-return monitor alarms on most days; under price-level noise returns become MA(1).

HYPOTHESIS: These are violations of the conditional null, not bugs. A null band absorbs small drift; persistent clustering needs the null to be restated relative to a forecast; noise needs block pre-averaging.

CHANGES IMPLEMENTED: `null_band` (E4c); E9 GARCH stress; `forecast.py` and E10; ALT1's block monitor.

EXPERIMENTS PERFORMED: E4c; E9 Q9c (2000 GARCH paths per cell); E10 Q10a to Q10e; ALT1 noise grid.

RESULTS: A band of 0.5 brought drift-driven false alarms below 5% in E4c (0.25 reduced them from 25% to 8.3%). Under GARCH with a half-life of an hour or more, no band up to 1 restores 5% (12% to 20% of days at band 1). Monitoring returns divided by a fitted GARCH forecast restores 3% to 6% per day, but a frozen persistent forecast passes on only a capped share of any break, (1 - β)/α at most (a doubling looks like +9% at persistence 0.99), and the Kullback-Leibler rate puts any detector's delay near two days. Price noise breaks the raw-return guarantee (8.5% to 41.5% per day); non-overlapping pre-averaged blocks restore it.

IMPROVEMENTS: A precise map of what the guarantee covers (`docs/limitations.md`) and a principled way to restate the null.

REMAINING ISSUES: The band is a heuristic, and the forecast-relative guarantee is only as good as the forecast.

NEXT ITERATION: Scale up the simulations.

---

ITERATION NUMBER: 7 (larger-scale simulations and computational profiling)

OBJECTIVE: Measure the guarantee's tightness at a scale where Monte Carlo error is small, and the cost per step.

PROBLEM IDENTIFIED: E4 used 300 paths per cell and estimated baselines, so it could not separate "the bound is loose" from "the baseline is off". Round 1's strided kernel sums took 30 s for 20,000 paths.

HYPOTHESIS: With the true baseline the false alarm rate is below m / A for every distribution with the right truncated mean, and roughly equal across distributions.

CHANGES IMPLEMENTED: E9 (4000 paths per distribution, 2000 GARCH paths per cell, profiling); E8 (10,000 paths per null cell); a loop over lags for kernel sums (round 1).

EXPERIMENTS PERFORMED: E9 Q9a to Q9d; E8 Q8a; timing of the streaming monitor.

RESULTS: One-day false alarm rates of 1.1% to 1.4% against the 5% bound for Gaussian, Student t3, jump and tick-rounded returns: valid, distribution free and about 4 times conservative. The streaming monitor runs at about 10 microseconds per step, so no compiled code is needed. Kernel sums went from 30 s to 4 s. A first run of the GARCH part crashed on duplicate pivot keys; the grid was fixed.

IMPROVEMENTS: Tight Monte Carlo error bars on every validity claim.

REMAINING ISSUES: Why the bound is 4 times conservative is studied by the theory track (E11).

NEXT ITERATION: Make the prototype usable by someone else.

---

ITERATION NUMBER: 8 (architecture, modularity and implementation)

OBJECTIVE: Make the prototype executable from a clean environment by someone who has not read the code.

PROBLEM IDENTIFIED: No exact dependency list, no command-line entry point, no end-to-end test; new features (forecasts, confidence sequences) needed their own modules.

HYPOTHESIS: A small CLI and file-to-alarm integration tests catch interface mistakes that unit tests miss.

CHANGES IMPLEMENTED: `requirements.txt` with exact versions; `python -m avmonitor simulate | run | estimators`; `tests/test_integration.py`; `forecast.py`, `confseq.py` as separate modules with backward-compatible options in `monitor.py`; `scripts/download_data.py` with a hashed manifest; an executed demo notebook; alternatives kept in `alternatives/` with their own tests.

EXPERIMENTS PERFORMED: CLI runs on synthetic data, on Binance ETHUSDT (October 2025) and on Nifty (May and June 2024, alarm at 09:24 IST on the 2024 general election counting day).

RESULTS: The first integration test failed because, under pandas 3, `DatetimeIndex.asi8` is not in nanoseconds; the fixture was rewritten. 164 tests pass from the repository root.

IMPROVEMENTS: One-command demonstrations and a clean package boundary.

REMAINING ISSUES: None blocking.

NEXT ITERATION: Independent verification.

---

ITERATION NUMBER: 9 (independent verification, regression testing and reproducibility)

OBJECTIVE: Have every claim checked by someone who did not write it, fix what is wrong, and make results reproducible.

PROBLEM IDENTIFIED: Authors check their own work badly. Before the review, several docs overstated claims (the impossibility result, noise robustness, what the monitor takes from the paper, the strength of the real-data check).

HYPOTHESIS: An adversarial review against the CSV files and the paper text will find wording problems and latent bugs, but no wrong mathematics.

CHANGES IMPLEMENTED: Independent review round 1 (2 critical, 7 major, 18 minor findings) and fixes to README, pitch, estimator notes, experiment report and math notes; the strictly predictable intraday profile (`start="past"`); E4e and E5c reruns; GitHub Actions running the full test suite on every push; `data/MANIFEST.csv` with SHA-256 hashes; E6 rerun from scratch after the profile change.

EXPERIMENTS PERFORMED: The reviewer's checks (now `experiments/e04e_review_checks.py`); reruns of E5c, E6 and E8 after the fixes.

RESULTS: No wrong theorem was found; every E1 to E5 and E9 number matched its CSV. The fixes changed wording, one bias figure (C1) and code paths that real data had not yet exercised; E5c, E6 and E8 real-data results moved by at most 2 points.

IMPROVEMENTS: The docs now say exactly what is proved, what is approximate and what is new.

REMAINING ISSUES: A second review round covers the work added overnight.

NEXT ITERATION: Final audit.

---

ITERATION NUMBER: 10 (final research audit, documentation and presentation readiness)

OBJECTIVE: Check that a researcher can reproduce every table from a clean clone, and that the pitch, README and report agree with the results.

PROBLEM IDENTIFIED: Many agents wrote parts of the report in parallel.

HYPOTHESIS: A clean-clone run with hash checks and a second review will catch inconsistencies.

CHANGES IMPLEMENTED: In progress (updated at the end of round 2).

EXPERIMENTS PERFORMED: In progress.

RESULTS: In progress.

IMPROVEMENTS: In progress.

REMAINING ISSUES: In progress.

NEXT ITERATION: Meeting with Prof. Bapat.

---

## Additional iterations (round 2)

ITERATION NUMBER: 11 (real data for all five coins and Nifty 50)

OBJECTIVE: Get real one-minute data for ETH, XRP, SOL, DOGE and Nifty 50.

PROBLEM IDENTIFIED: Every exchange API, Yahoo, Kaggle, Hugging Face, Zenodo, Stooq and www.nseindia.com were refused by the cloud network policy.

HYPOTHESIS: Public GitHub repositories mirror these data, and raw.githubusercontent.com is allowed.

CHANGES IMPLEMENTED: `data.load_binance`, `data.load_nifty`, `regular_sessions`, `scripts/download_data.py` and the hashed manifest (details in `docs/iterations/MAIN.md`, R2-1).

EXPERIMENTS PERFORMED: Offline loader tests; first comparison with Table 5 (followed up in E6).

RESULTS: 140 files: Binance spot one-minute klines for all five coins (2024 to 2025, no missing minutes in July to December 2025) and Nifty 50 one-minute bars (April 2017 to October 2026, 2,369 sessions). India VIX daily data from archives.nseindia.com (E7).

IMPROVEMENTS: Every claim of the paper about its five coins became testable, and an Indian market application became possible.

REMAINING ISSUES: The data are third-party mirrors; the loaders report every check and the manifest pins the bytes.

NEXT ITERATION: E6 and E7.

---

ITERATION NUMBER: 12 (Nifty 50 application)

OBJECTIVE: Apply the monitor to Nifty 50 one-minute data and check alarms against India VIX and scheduled events.

PROBLEM IDENTIFIED: A 375-bar session market has an opening and closing pattern that a 24-hour profile mishandles, and event days need verified dates.

HYPOTHESIS: A session-aware profile and a development and held-out split give an honest test.

CHANGES IMPLEMENTED, EXPERIMENTS, RESULTS: See `docs/report_parts/E7.md` and `docs/iterations/E7.md` (six iterations, including a flawed first control that was replaced).

IMPROVEMENTS, REMAINING ISSUES: As in those files.

NEXT ITERATION: Event studies for crypto (E12).

---

ITERATION NUMBER: 13 (forecast-relative monitoring)

See `docs/iterations/MAIN.md`, R2-5, and `docs/report_parts/E10.md`. In short: dividing returns by a fitted GARCH forecast restores the guarantee under clustering (3% to 6% against 82% to 100%), at the cost of a capped signal for structural breaks, with a closed-form cap and an information bound, both checked by simulation; on real data, alarm rates barely change, and the plug-in baseline was found to break the bound over five days.

---

ITERATION NUMBER: 14 (detection delay theory, and a design change it suggested)

OBJECTIVE: Explain how fast the detector reacts, what robustness costs in speed, and whether the theory points to a better design.

PROBLEM IDENTIFIED: E4 measured that the robust detector is 14%, 35% and 73% slower than a Gaussian likelihood-ratio detector on clean data, with no explanation of where the gap comes from or whether it can be closed.

HYPOTHESIS: The delay is about ln A divided by the growth rate of the best bet, so a ratio of information numbers should predict the price of robustness, and the bound 0 ≤ X ≤ τ² should limit the evidence against a fall.

CHANGES IMPLEMENTED: `docs/theory_delay.md` (growth rate, its dual form, a proved delay bound, a second-order formula, the information ladder, the decrease speed limit log(B / (B - 1)) with B = τ²/m₀, asymptotic optimality); `src/avmonitor/design.py`; then a new monitor option `c_down`, a separate, tighter clip for the decrease side, with two tests.

EXPERIMENTS PERFORMED: E11 (Q11a to Q11e: 20 delay cases, the price of robustness, bet grids, the slack in the bound, optimality as A grows to 10^8) and E11f (c_down = 3, 2, 1.5; four return laws, 2000 null and 600 change paths per cell; five coins, July to December 2025).

RESULTS: The second-order formula is within 3% of the simulated delay in all 20 cases; the information ratio predicts the penalties (1.16, 1.37, 1.95 against 1.14, 1.35, 1.73 simulated). For decreases truncation caps the evidence per step at log(B / (B - 1)), about 1/c², for every bet that is valid over the whole null: at c = 3 and A = 7800 even a total collapse of volatility takes any Shiryaev-Roberts e-detector built from valid bets at least 58 steps (64 at 2A), and the default five-bet grid 73 (79). c_down = 2 cut median decrease delays by 13% to 51% for Gaussian, Student t4 and jump returns and by 32% to 52% with ticks of 0.25 and 0.75 sd, but was 12% slower for a 20% fall at a tick of 0.5 sd (paired difference +29 steps, standard error 1.7); one-day false alarms were 5.0% to 6.7% with the plug-in baseline (empirical only); on real data it was neutral to slightly positive. (These numbers are after review round 2b, which found that the first version's 73-step floor holds for the default grid only and that its rounded-price model scaled the tick with the volatility.)

IMPROVEMENTS: A theory-driven, tested design option, and a quantitative account of the price of robustness and of the slack in the false alarm bound.

REMAINING ISSUES: Review round 2b is fixed (`docs/review/independent_review_round2b.md`; C1, C2, M1, M2 and the minors). c_down now works with the confidence-sequence baseline (each side at δ/2), but E11f has not been rerun with it; c_down is an option, not the default.

NEXT ITERATION: Panel monitoring (15).

---

ITERATION NUMBER: 15 (many assets at once)

OBJECTIVE: Monitor a panel of assets with a guarantee that does not degrade with the number of assets.

PROBLEM IDENTIFIED: K single-asset monitors at level α each give up to about Kα chance of a false alarm per day, and a Shiryaev-Roberts statistic divided by time is not an e-value, so e-BH cannot be applied to it directly.

HYPOTHESIS: Per-asset mixture e-processes with a sticky e-BH rule control the false discovery rate at every stopping time, under any dependence between assets.

CHANGES IMPLEMENTED: `docs/theory_panel.md` (Theorem 1 with proof), `src/avmonitor/panel.py`, `experiments/e15_panel.py`.

EXPERIMENTS PERFORMED: 48 synthetic designs (independent assets and a common factor with correlation 0.7 and common jumps), horizons up to ten blocks, the confidence-sequence baseline, and the five coins plus Nifty 50.

RESULTS: Observed false discovery rate at most 0.96% against 5%. On real data discoveries are confirmed 78% to 97% of the time against about 50% at random times; the gain over a Bonferroni split is mostly earlier flags, not more of them.

IMPROVEMENTS: A proved multi-asset guarantee, with its limits stated (the chance of any false discovery once real changes are present is not controlled).

REMAINING ISSUES: As above; the panel theory awaits an independent review (round 2b).

NEXT ITERATION: Review round 2a (16).

---

ITERATION NUMBER: 16 (independent review round 2a and its fixes)

OBJECTIVE: Have the overnight work (math.md sections 5 and 6, E8 to E12, ALT3, README, strategy) checked by a reviewer who did not write it.

PROBLEM IDENTIFIED: Many parts were written in parallel, and top-level documents quoted numbers from report parts that later changed.

HYPOTHESIS: Tracing every number to its CSV will find copying errors and overclaims, but no wrong theorem.

CHANGES IMPLEMENTED: `docs/review/independent_review_round2a.md` (8 critical, 13 major, 21 minor findings, 7 derivations). Fixes: `forecast.ewma_variance` rejects a non-finite start value and Q10c starts its EWMA on the first day with data; Q10d reports a pooled ratio of means with a bootstrap standard error; corrected numbers and wording in README, pitch, strategy, math.md, E8, E9, E10, E12 and ALT3 (for example 1.8 to 7.6 times the delay, not 7.2; the forecast cap described as a cap on a ratio, not a fraction of the shift; "restores the guarantee" only with a correct forecast model).

EXPERIMENTS PERFORMED: The reviewer's reruns (Theorem 6.4 checked step by step; the clipped-RQ variance 8.19 by integration and Monte Carlo); our reruns of Q10c (the summary CSV is reproduced exactly) and Q10d (every earlier column unchanged).

RESULTS: No wrong theorem. The critical findings were copying errors and overclaims; one result (Q10c) could not be reproduced from the committed script until the start-value fix.

IMPROVEMENTS: Every top-level number now traces to a CSV; the main project and ALT3 no longer contradict each other about improving the estimator.

REMAINING ISSUES: E12 control matching (weekday-matched pools, band confound) is being redone.

NEXT ITERATION: Benchmarks against published methods (E13), running truncation (E14), sequential estimation (ALT6), final audit (10).

---

ITERATION NUMBER: 17 (adversarial search for null streams, E16)

OBJECTIVE: Find how close a stream that satisfies every condition of Theorem B can push the full default detector to its bound, and where exactly the bound breaks when a condition fails.

PROBLEM IDENTIFIED: E9 found the bound about four times conservative for four natural laws (1.1% to 1.4% per day against 5%). Nobody had searched for the worst law for the full two-sided grid, or allowed the law to react to the detector's own state, which the theorem permits. "Four times conservative" could have been an artefact of friendly test laws.

HYPOTHESIS: Both statistics are affine in the next observation, so a two-point law with mean 1 is the adversary's best choice at each step; "bold play" (the widest legal spread, landing exactly on 2A when one draw can reach it) should be close to optimal. Outside the null class, each failure should line up with the condition that was broken.

CHANGES IMPLEMENTED: `experiments/e16_adversarial_null.py` (a vectorised engine with a per-step callback, bold play, an 8-parameter state-dependent policy family, hidden-Markov shape switching, GARCH, plug-in and price-tick designs) and `tests/test_e16_adversarial.py` (every searched law keeps the conditional mean at 1, each return law uses its own true τ and m0, and the engine matches `VolatilityMonitor` alarm for alarm).

EXPERIMENTS PERFORMED: Grid search, then Nelder-Mead on common random numbers, over nine families; the best law of each family re-run on seeds the search never saw (200,000 paths at one day, 60,000 at five days). Outside the class: 64 GARCH cells, the plug-in baseline, and returns rounded to price ticks.

RESULTS: Failed first: the in-sample search maximum was biased upwards (a small probability estimated on a few thousand paths, maximised over many candidates), so every reported rate now comes from fresh seeds. Inside the class the worst stream found is bold play at 2.62% per day (95% interval 2.50% to 2.75%) against 5%, and 18.1% over five days against 25%; the Gaussian reference is 0.97% and 12.4%. Under bold play about half of the expected wealth is stranded on paths that never alarm, and that accounts for almost all of the remaining slack under the laws tried; a worse law is not ruled out. Outside the class: GARCH is above the bound from a drift of 3.9% in the conditional mean at one day (1.6% over five days; 43 of 64 cells at one day, ordered by that drift); the plug-in baseline reaches 41% to 42% over five days, and 32.5% in one day under a shape change; price ticks fail only from 4 sd (over five days; 6 sd at one day).

IMPROVEMENTS: "About four times conservative" now reads "about four times for natural laws, about two times at one day against the worst stream found" in the README, limitations and the delay theory note.

REMAINING ISSUES: No proof that bold play is optimal for the two-sided grid, so the remaining slack is not shown to be structural. Independent review round 5 (`docs/review/independent_review_round5_e16.md`) found no critical problem; it asked for these qualifiers, a proper in-class control for price ticks, and a correction to the delay theory note, where an earlier dynamic programme had been read as the worst case.

NEXT ITERATION: A bound for false alarms as a function of drift in the conditional mean, the condition GARCH breaks (TH4).

---

ITERATION NUMBER: 18 (one monitor for several sampling frequencies, E17)

OBJECTIVE: Combine e-detectors run on 1, 5 and 15 minute returns into one monitor with one proved false alarm bound, and find out whether it beats each frequency.

PROBLEM IDENTIFIED: A coarse detector changes only at the end of its block, so between block ends it is not a fair bet. The usual proof does not cover the plain average checked every minute.

HYPOTHESIS: Counting each coarse detector only at the wealth it is sure to keep (a "certified" average) restores a proof. The combination should never be much worse than the best single frequency and should catch changes that one frequency misses.

CHANGES IMPLEMENTED: `src/avmonitor/multifreq.py` (`MultiFrequencyMonitor`, `proved_bound`, `level_for_budget`), `docs/theory_multifreq.md` (Theorems 1 to 3), `experiments/e17_multifrequency.py`, `tests/test_multifreq.py`.

EXPERIMENTS PERFORMED: Levels and handicaps; false alarms under four laws and with the plug-in baseline; persistent changes; short bursts and a jump in lag-one autocorrelation; price noise (γ = 1, then γ = 2, added after seeing γ = 1); the five coins with settings committed before the first real run; after review round 4, a rate-matched comparison on the coins.

RESULTS: A union bound proves only three times the budget for the plain average; the certified average has P(false alarm within m) ≤ Σ w_k ⌈m/k⌉ / A when each frequency's null holds with the true baselines. With independent simulated returns it never beats the best single frequency on persistent changes, staying within 0% to 5% (equal weights) or 11% to 18% (balanced weights) of it, and it catches 95% to 97% of paths after an autocorrelation jump from 0 to 0.3 (84% to 89% for Student t3), against 3.8% for the one-minute detector. On the coins it behaves like the one-minute monitor. Failed or corrected: review round 4 found a calibration bug with missing minutes (no reported number changed; every result file was byte-identical after the fix) and overstated wording; an unplanned comparison with E6 suggested the 5-minute detector was more precise, but at matched alarm counts it is not (up alarms confirmed 90.0% against 90.5%, p = 0.88); the first rate-matched run capped the levels too low and never reached the target counts, and was rerun.

IMPROVEMENTS: Every statement now carries its conditions (true baselines, independent returns); calibration is safe with missing minutes; `MultiFrequencyMonitor` takes separate increase and decrease levels.

REMAINING ISSUES: No proof with the plug-in baselines used on real data. The gain under price noise holds only at γ = 2, a level chosen after seeing γ = 1.

NEXT ITERATION: Landmark-event case study (E18) and the final audit (10).

---

ITERATION NUMBER: 19 (bounds when the null holds only approximately, TH4)

OBJECTIVE: Say how much the false alarm probability can grow when the conditional mean of the clipped squared return drifts above its baseline, the condition that volatility clustering, price ticks and the plug-in baseline break (E16).

PROBLEM IDENTIFIED: Theorem B is all or nothing: it holds when the conditional mean stays at its baseline and says nothing when it drifts by 1%. E16 showed failures outside the class but no theory for how fast the bound degrades.

HYPOTHESIS: Each bet's factor has conditional mean at most 1 plus its bet size times the drift, so an almost-supermartingale argument (in the spirit of Robbins and Siegmund 1971) should give the bound m/A times an explicit inflation factor, and a deterministic stream at the maximal drift should attain it.

CHANGES IMPLEMENTED: `docs/theory_robustness.md` (Lemma 1, Theorem 1 in sup form, Propositions 1 to 5, Lemma 2 and Theorems 2 and 3 in pathwise form, Corollary 1 for the plug-in baseline), `src/avmonitor/robustness.py`, `tests/test_robustness.py`, `experiments/th4_robustness.py`.

EXPERIMENTS PERFORMED: A design table of drifts at which the one-day bound reaches 10% and 25%; constant and switching drifts; every GARCH, plug-in and price-tick cell of E16 against the proved bounds; exact counterexamples.

RESULTS: The one-day 5% bound reaches 10% at a drift of 1.63% and 25% at 2.66% (downward drift only: 6.36% and 10.51%), and 1 at 3.740%, where a deterministic stream alarms within the day from 3.756%. No proved bound is below a measured rate anywhere. For GARCH the proved bound is below 1 in only 22 of 64 one-day cells, all with a small drift. The plug-in baseline is bounded by averaging over calibration windows: 41.8% at one day with 390 calibration steps, 1.29 times the worst stream E16 found. Failed or corrected: an independent review (round 6) found Theorem 3 false as first written (its pathwise term must be the maximum over time, built from bounds that do not respond to the data; two counterexamples), the claim that the mean drift gives a lower estimate false, and "sharp" overclaimed (the bound is attained only where it reaches 1). One test was vacuous (every bound equal to 1) and was replaced. No reported number changed.

IMPROVEMENTS: The guarantee now has a quantitative robustness statement, and the plug-in failure of E10, E14 and E16 has a proved upper bound.

REMAINING ISSUES: The exact worst case below the deterministic alarm step; a bound for price ticks (needs the drift given the whole past); a GARCH bound that is informative where the nominal bound fails.

NEXT ITERATION: The referee-style review of the paper (major revision) and the work that answers it (TH5, TH6, E19 to E21).

---

ITERATION NUMBER: 20 (landmark shocks as pre-registered case studies, E18)

OBJECTIVE: See how fast the monitor and published methods flag four well-known market shocks: the Indian election results of 4 June 2024, the global sell-off of 5 August 2024, the tariff sell-off of 7 April 2025 and the crypto liquidation cascade of 10 October 2025.

PROBLEM IDENTIFIED: Event studies on scheduled releases (E7, E12) say nothing about large unscheduled shocks, and case studies invite cherry-picking of dates, reference times and competitor settings.

HYPOTHESIS: With dates and reference times verified from news sources and every setting committed before any run, the comparison is honest even with four events; the e-detector should be competitive in speed with BOCPD and slower than nothing it is calibrated against.

CHANGES IMPLEMENTED: `docs/preregistration_e18_events.md` (committed before the script existed), `experiments/e18_landmark_events.py`, thresholds matched to the e-detector's alarm rate outside the event windows.

EXPERIMENTS PERFORMED: The primary streaming design on Nifty 50 and the five coins; three pre-registered sensitivity designs (fresh calibration just before each reference time, other start dates, whole-sample matching); a descriptive look at the hours before each crypto alarm (added after the run).

RESULTS: On Nifty the e-detector flagged all three shocks, 11, 104 and 172 minutes after the open. On crypto it flagged none of the fifteen coin and event pairs within 24 hours: on each event it had raised an up alarm 4.5 to 6.9 hours earlier, as the turbulence began, and was still recalibrating. Competitors were mostly paused for the same reason (1 to 5 of 45 configurations watching), and 22 of their 35 crypto catches came in the first minute after a pause ended; the day-long two-window variance-ratio scan flagged 10 October 2025 on all five coins within 17 to 73 minutes. With a fresh calibration the e-detector flagged all 18 cases in 15 to 207 minutes, earlier than the short-window competitors on 5 August 2024 and later than plain BOCPD on 10 October 2025 and on two Nifty events. A container restart killed the crypto BOCPD runs midway; only the lost jobs were rerun, with no setting changed.

IMPROVEMENTS: A concrete, pre-registered demonstration of the cost of the blind period after an alarm, the problem E14's monitors address.

REMAINING ISSUES: Four events are case studies, not evidence. Rate matching uses data from after the events, which favours the competitors.

NEXT ITERATION: E14's monitors without a blind period on the same events, exploratory and decided after this result: at the e-detector's null band and a similar alarm rate (0.59 to 0.94 against 0.66 to 0.71 per crypto day) they flagged 15 and 17 of the 18 cases within 24 hours, fast on 5 August 2024 and 10 October 2025 and slow on 7 April 2025, against 3 for the e-detector; their alarm lists match E14's exactly.

---

ITERATION NUMBER: 21 (the real-data inference, re-analysed after the referee report, E20)

OBJECTIVE: Answer the referee's point that the real-data inference used unadjusted subgroups, pooled dependent alarms and a weak notion of "information beyond a forecast".

PROBLEM IDENTIFIED: Several tests per table without adjustment; p-values that treat alarms on the same calendar day across coins as independent; post hoc comparisons reported next to pre-specified ones; and hit rates against slow forecasts, which a trailing ratio rule also achieves.

HYPOTHESIS: The comparisons with random times are strong enough to survive any adjustment, the comparisons with the ratio rule and the longer-window competitors are not, and against a forecast with an intraday component the alarms add little.

CHANGES IMPLEMENTED: `docs/preregistration_e20.md` (families, Holm, day-block bootstrap, regressions and loss comparison fixed before any adjusted result), `experiments/e20_inference.py`, `tests/test_e20_inference.py`.

EXPERIMENTS PERFORMED: Classification of 36 real-data claims; Holm within families; a bootstrap resampling calendar days jointly across coins; regressions of next-day log realised variance on log forecasts (EWMA, HAR, HAR with the last hour's realised variance) and alarm indicators with day-clustered errors; a QLIKE comparison of forecasts scaled at alarms.

RESULTS: One claim rests on a pre-specified primary test (FOMC on BTC holds; CPI in the same family does not). All 16 comparisons of alarms with random times or controls hold (Holm p ≤ 0.0016); 1 of 48 comparisons with the ratio rule does (the post hoc subgroup the paper had highlighted); the 6-point margin over the best day-long competitor does not (Holm 0.44 and 0.53). Against HAR with the last hour's realised variance, increase alarms add nothing significant (Holm 0.32) and decrease alarms a little (Holm 0.046), no more than a 240-minute ratio rule; no QLIKE gain out of sample. The rank restart of E14 had 6.7 points fewer confirmed up alarms than the plug-in (Holm 0.006). Failed first: a pandas key-matching bug (NaN keys dropped), a family of 12 tests where 8 were pre-registered (corrected; the 12-test version is not reported), and an alarm-count check that missed alarms off the grid.

IMPROVEMENTS: Every real-data claim now has a class (primary, secondary, exploratory) and an adjusted status; README, pitch and limitations now say that alarms beat random times but not a strong forecast.

REMAINING ISSUES: Only one pre-specified primary test exists; new confirmatory evidence would need a new pre-registered study (the live Nifty run is a demonstration, not that study).

NEXT ITERATION: The paper's real-data section is rewritten with these classifications; E19 (real-data null streams) and E21 (stronger benchmarks).

---

ITERATION NUMBER: 22 (limits for every valid procedure and the plug-in regime, TH5)

OBJECTIVE: Answer the referee's points M3 and sections 2 and 3: is the decrease speed limit a property of e-detectors or of the problem, how far are the e-detectors from the best possible delay, when does the plug-in baseline keep its level as n, m and A grow together, and why does the growth-optimal decrease bet sit at its upper limit?

PROBLEM IDENTIFIED: The 58-step speed limit of `docs/theory_delay.md` bounds bets, not rules; Lorden-type bounds are asymptotic; TH4's plug-in bound averages over windows but gives no rate in n.

HYPOTHESIS: The null class contains independent two-point laws that put the largest possible mass on any observed value, so a valid rule cannot stop on a path of small null probability; and the plug-in's drift is controlled by the Kolmogorov distance between the window and the return law, of order 1 / √n.

CHANGES IMPLEMENTED: `docs/theory_limits.md` (Theorem 2 pathwise floor, Theorem 5 delay bound, Proposition 6 asymptotics, Theorem 8 and Corollary 10 plug-in regime, Proposition 11 margin, Proposition 12 decrease bet at the boundary), `src/avmonitor/limits.py`, `experiments/th5_limits.py`, `tests/test_limits.py`.

EXPERIMENTS PERFORMED: Closed-form floors for c = 1.5 to 5; delay bounds for 18 post-change laws (10,000 paths of 4,000 steps, levels A and 2A); asymptotics up to A = 10⁸; plug-in bounds and simulations for n = 390 to 92,160 and m = 98 to 7,800; the hold-then-bold stream; the decrease bet for nine laws.

RESULTS: No valid procedure alarms on a complete collapse at c = 3 before step 45 (e-detectors 58, default grid 73). The default grid takes 1.69 to 2.16 times the lower bound on the mean delay and the best single bet 1.34 to 1.76 times; no alarm on 180,000 paths came before the pathwise floor. The plug-in keeps m / A up to 1 + O(m / √n) when m is small compared with √n (C = 0.037 Gaussian, 2.52 for every law) and fails by an unbounded factor beyond; at one day after a 390-step window the worst case lies between 33.8% and 42.4%, the hold-then-bold stream gives 36.5% and i.i.d. returns at most 5.14%. A margin z = Φ⁻¹(1 - δ / 2) restores m / A + δ. The decrease optimum is at its limit only for light tails (ε* = 2.4 × 10⁻⁸ for a Gaussian halving; λ* / λ_max 0.93 to 0.55 for Student t3). Failed first: the least favourable law taken from a tolerance-based maximiser had mass 0.982 (replaced by a root solve; no number changed), overflow in the plug-in bounds (fixed by integration in log space), and a wrong grid factor in a first collapse evaluation (caught before any number was written).

IMPROVEMENTS: The speed limit is now a statement about every valid procedure, the e-detectors' efficiency has a number, and the plug-in's failure has a location (m ≍ √n) and a remedy.

REMAINING ISSUES: The delay bound is for a change at the first step and is a Monte Carlo estimate; the distribution-free plug-in bound is informative only for m below about √n / 20; the i.i.d. boundary is a heuristic; three 2025 to 2026 references were checked from search listings only, because arxiv.org is blocked from the build machine.

NEXT ITERATION: TH6 (what Theorem B says per period, and an algorithm-free form of Theorem D).

---

ITERATION NUMBER: 23 (what Theorem B says per day, and a bound for every day, TH6)

OBJECTIVE: Answer the referee's points M2 (the documents read Theorem B as "5% per day") and M4 (Theorem D's hypothesis depended on the bets the algorithm chooses).

PROBLEM IDENTIFIED: Theorem B bounds the first false alarm within m steps of a (re)start. Read as a bound for every day it is wrong, and README, pitch, limitations and the paper used that reading in several places.

HYPOTHESIS: A null stream exists that pushes later days above m / A, so no per-day bound follows from Theorem B; restarting the statistics at each day boundary restores one, at a delay cost that a renewal argument can bound.

CHANGES IMPLEMENTED: `docs/theory_per_period.md` (Proposition 1, Theorem 1 restart at day boundaries, Proposition 2 delay bound, Theorem 2 e-CUSUM window bound, Proposition 3 park and strike), a remark after Theorem B and Theorem 6.4 restated in `docs/math.md`, `src/avmonitor/periodic.py` (`PeriodicVolatilityMonitor`), `experiments/th6_per_period.py`, `tests/test_periodic.py`; per-day sentences corrected in README, pitch and limitations.

EXPERIMENTS PERFORMED: Hazards over 20 days for the four laws of E9 at 390 and 1440 steps per day; the referee's calm-then-active family; E16's bold play for 20 days; park and strike; delays of the plain, restarted and e-CUSUM detectors (20,000 paths) and the renewal bound; the tolerance band under GARCH(1,1).

RESULTS: The constant null law Y = 1 makes the default detector alarm with certainty on day 40. For i.i.d. Gaussian returns the hazard is 0.98% on day 1 and 4.70% on average over days 11 to 20; bold play averages 5.25% (pooling chosen after seeing the table); the referee's family gives 11.3%, 28.0% and 65.1% on days 5, 10 and 20; Corollary C held (0.037 to 0.048 per day). The restarted detector stayed at most 2.75% (390) and 3.52% (1440) on any day, at 1.06 to 1.44 and 1.02 to 1.14 times the plain detector's delay. The e-CUSUM bound is 5.146% (390, 7800), measured at most 0.40%, and fails under martingale differences (6.16% on day 20 for park and strike). Under GARCH the tolerance band needed has median 0.90 to 5.01. Failed first: a two-sided window bound halved by mistake (fixed before use), a wrong B in park and strike, CUSUM counts without restart that meant nothing (now missing), GARCH tolerance on fresh returns instead of the driving ones (medians moved by at most 0.031). The renewal bound was added after the delay results were seen.

IMPROVEMENTS: Every per-day statement now says which day it covers; a monitor with a proved bound for every day exists, with its delay cost measured.

REMAINING ISSUES: The confidence-sequence baseline's tolerance is rarely met under clustering; the e-CUSUM needs exchangeability; a rate-matched delay comparison was dropped.

NEXT ITERATION: E19 and E21, then the paper's per-day sentences and TH6 slots.

---

Further work streams (E13, E14, E16, E17, the alternative prototypes, logs ALT1 to ALT6) keep their own iteration notes in `docs/iterations/`. Their main results are summarised in `docs/experiment_report.md` and `alternatives/ALTERNATIVES.md`.

---

## Appendix: round 1 log (2 November 2026, unchanged)

| # | What we tried | What happened | What we changed |
|---|---|---|---|
| 1 | Read the paper in full and checked the earlier summary's two concerns against the printed formulas. | Both concerns follow directly from the algebra. Also found that eq. (16) is a rescaled truncated RQ, and that eq. (15) reverses the weight order. | Planned E1 to test whether the published tables come from the printed formulas, instead of assuming it. |
| 2 | Tried exchange APIs (Binance, Coinbase, Kraken, Bybit, OKX, Yahoo) for real data. | All blocked by the network policy of this environment. | Found a public GitHub-hosted Bitstamp BTC/USD 1-minute file covering 2025 and 2026 (raw.githubusercontent.com is reachable). Wrote `data.load_bitstamp_csv`. Other coins and Nifty are left for later. |
| 3 | First unbiasedness check of all classical estimators under constant volatility. | MinRQ came out 1.247 times the truth. The constant π / (π - 2) belongs to MinRV, not MinRQ. | Corrected it to π / (3π - 8) (Andersen, Dobrev and Schaumburg 2012). A unit test now guards every constant. |
| 4 | Computed kernel sums with a strided sliding-window view and matrix product. | 20,000 × 390 paths took about 30 s. | Replaced with a loop over the h + 1 lags. It now takes about 4 s. |
| 5 | Unit test: threshold RQ should ignore a jump of 0.05 at σ = 0.3 over 390 steps. | Failed. That jump is only 3.3 diffusive standard deviations, inside any sensible threshold. | Moved the test to a 13 sd jump. **Lesson:** moderate jumps cannot be told apart from diffusion by any threshold rule, so "jump robust" only ever means robust to *large* jumps. |
| 6 | First null check of the e-detector: A = 1000 over 5000 steps. | 99% of paths raised a false alarm. This looked like a bug. | It is not a bug. The Shiryaev-Roberts statistic has ARL very close to A, so a horizon five times A is expected to alarm. We measured ARL over 30,000 steps (about 1160 for A = 1000, as the theory says) and switched to choosing A from a per-day bound: P(false alarm within m steps) ≤ m / A (Theorem B, which we proved in `docs/math.md`). |
| 7 | First detection experiment: A = 1000 with the change at step 1000. | About half the paths alarmed before the change, so the delay figures were meaningless. | Used A = 7800 (5% per 390-step day) and moved the change to step 500. |
| 8 | Early detection test calibrated on clean data and monitored data with jumps. | The robust detector alarmed. Correct behaviour, wrong experiment: adding jumps *is* a change. | Calibrate on data from the same process as the monitored pre-change data. Jumps and fat tails are then part of "normal", which is what the robustness claim is about. |
| 9 | Ran E4 with slow stochastic volatility inside the null. | Every detector, robust or not, alarmed on about 25% of days. | Recognised that the null "variance never exceeds baseline" is false under drift. Added `null_band` and ran E4c: a band of 0.25 brings false alarms to 5% to 8% per day and keeps detection of large changes. |
| 10 | Streaming monitor unit test with a 500-step calibration window. | The false alarm rate was 67% against a bound of 50%. | Traced to error in the estimated baseline (the theorems assume the true one). Added `margin_z` to `calibrate` and ran E4d. One day of calibration gives 5.4% against a 5% bound, and a margin of one standard error fixes it. Documented in `docs/math.md`, section 5.6. |
| 11 | Hypothesis: a bet sized from the local truncated quarticity (aGRAPA) detects faster than a fixed grid of bets. | No meaningful difference in any scenario (for example 86 against 87 steps). | Kept it as an option. The default is the simpler grid. Reported as a negative result. |
| 12 | E1: fitted one volatility to Table 7 and predicted the rest. | The plain kernel and the small-jump robust rows reproduce well. The large-jump (σ_J = 0.10) robust rows are about 5 times too small in our simulation. | Reported the discrepancy instead of tuning extra parameters to hide it. It is listed as a question for the author. |
| 13 | Raw-r² and Gaussian benchmark detectors overflowed with large jumps. | Warnings and `inf` statistics. | Capped exponents and silenced overflow inside the recursions. An infinite statistic still correctly means "alarm". |
| 14 | First real-data monitor run reported the share of alarms "followed by higher volatility". | Without a reference, 79% is hard to interpret because volatility clusters. | Added a control: the same ratio at 2000 random times (50%, median 1.00). |
| 15 | Real-data daily IQ on a common scale. | MAD-truncated estimators came out 3 to 7 times below MedRQ, MinRQ and TPQ. | Interpreted with E2b (truncation bias under intraday seasonality and fat tails). Real data strengthens the case for local, online MAD and for using truncation mainly to bound the monitored statistic. |
| 16 | Checked all written text for em and en dashes. | Removed every one. | Standing rule for this project. |

Notes added in round 2 to this appendix: row 9's "5% to 8%" is 8.3% and 5.0% for the two sides (E4c); row 10's 5.4% is within Monte Carlo error of 5% and the informative comparison is with the true-baseline rate (E4e); row 11 was refined by the review (slightly slower in 14 of 18 cells).
