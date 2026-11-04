# EV1 iteration notes: crypto event study (US macro releases)

Owner: EV1. Files: `experiments/e12_crypto_events.py`, `data/events/us_macro_events.csv`, `results/e12*`, `results/figures/e12*`, `docs/report_parts/E12.md`, this file. Times are UTC, night of 2 to 3 November 2026.

---

## Iteration 1. Verified event calendar and the streaming event study

**OBJECTIVE.** Build a verified list of FOMC statements and US CPI releases for 2024 and 2025 with exact UTC times, run the frozen E5c monitor on Binance BTC, ETH and SOL, and test whether up alarms cluster in the 6 hours after a release.

**PROBLEM IDENTIFIED.** (1) The official sources (federalreserve.gov, bls.gov) could not be fetched directly: the session's network policy does not resolve those hosts. (2) Event-time conversion must handle the EST/EDT switch, which differs by year. (3) The 2025 federal shutdown changed the CPI calendar, so a mechanical "second Wednesday" rule would be wrong. (4) A naive before/after comparison would be confounded by time of day and weekday, and alarms of the same monitor in nearby windows are dependent.

**HYPOTHESIS.** Up alarms are more likely in the 6 hours after a release than in matched control windows; they come after the release, not before; down alarms follow within two days as volatility normalises; the Gaussian LR benchmark alarms around releases because of the first big minute.

**CHANGES IMPLEMENTED.** Web search used as the verification route: the search tool indexes the official calendar and archive pages (titles such as "Consumer Price Index News Release - 2024 M01 Results" at `bls.gov/news.release/archives/cpi_02132024.htm`, the Fed meeting calendar, statement PDFs headed "For release at 2:00 p.m. EDT"), and each row of `us_macro_events.csv` carries the URL and what the search showed. Rows for the cancelled October 2025 CPI (status `cancelled`, original date 13 November 2025) and the December 2025 CPI (13 January 2026, outside the sample) are kept as documentation. The script converts "date + 14:00 / 08:30 America/New_York" to UTC with the IANA database and stops if the CSV's UTC column disagrees. Control windows: same weekday, same Eastern clock time, at least 2 days from any event day. Test: exact stratified permutation test (Poisson-binomial tail), pooled one-sided Fisher as a cross-check, primary family fixed as BTC x {FOMC, CPI} x up-within-6-hours with Holm, everything else secondary. Per-alarm diagnostics (largest |return| in the 5 minutes to the alarm in scale units; robust volatility of the next 1440 minutes over the baseline) for the Gaussian comparison.

**EXPERIMENTS PERFORMED.** Streaming monitor (A = 28,800, c = 3, cal_len 1440, band 0.25) and Gaussian LR (band 0) on the three coins, January 2024 to December 2025, 39 events, one run (about 1 minute).

**RESULTS.** (Superseded by iteration 3; the numbers below are from this run.) BTC, robust, up alarm within 6 h: FOMC 9/16 against 14% of control windows (p = 0.0001, Holm 0.0002); CPI 5/23 against 11% (p = 0.13). ETH 8/16 and 7/23, SOL 9/16 and 9/23, all p <= 0.004 except BTC CPI. No up alarm before the release in 0 to 1 of 16 FOMC and 0 to 1 of 23 CPI events per coin (controls 6% to 12%). Down alarms in 6 to 48 h: 61% to 69% of events against 62% to 70% of controls. Gaussian LR: 1 to 7 events with an alarm per cell, p >= 0.07.

**IMPROVEMENTS.** A verified calendar with the shutdown changes recorded; a test design that matches on clock time and weekday and accounts for the multiplicity.

**REMAINING ISSUES.** The streaming monitor spends 70% of the time recalibrating (0.7 alarms per day, each followed by a 1440-minute calibration) and the Gaussian monitor 92%. It was awake at only 9 of 16 FOMC and 9 of 23 CPI releases on BTC. Among those, it alarmed at 9/9 FOMC releases. So the event rates measure "awake and alarmed", and the Gaussian comparison is uninformative (awake at 0 to 4 events). The figure's rug marks and legends overlapped.

**NEXT ITERATION.** A design in which both detectors are awake at every anchor.

---

## Iteration 2. Fresh-calibration design, control diagnostics and figure

**OBJECTIVE.** Measure what the detector can do at a release when it is awake, for both detectors, without changing the frozen settings.

**PROBLEM IDENTIFIED.** The endogenous blindness above. Restricting to windows with the monitor awake (the `monitoring at start` sample) helps but selects on the stream's own history (9 to 12 events per cell) and leaves the Gaussian benchmark with 0 to 4 events.

**HYPOTHESIS.** With a monitor calibrated on the 24 hours before each anchor (event or control) and run for 48 hours, the robust event rate rises towards the awake rate (near 100% for FOMC), the control rate rises to the base rate of a fresh 6-hour window, and the Gaussian monitor alarms in most windows of either kind because single large returns are common.

**CHANGES IMPLEMENTED.** `fresh_window`: for each anchor a new monitor on `u[i0 - 1440 : i0 + 2880]` (post windows) and one on `u[i0 - 1800 : i0]` (pre window); the script checks that calibration used exactly 1440 returns (it did for every anchor). The same stratified test is run for the new design; a `design` column and a Holm correction over the six robust tests within each design were added to `e12_tests.csv`; `e12_control_windows.csv` exports every control window's outcomes; the delay table gained the conditional "down alarm after an up alarm" rate for events and for control windows; the Gaussian-versus-robust table gained rows for the first up alarm in event windows and in control windows under the fresh design. Figure: rug heights now use the row's shared y-maximum, one shared legend below the panels, alarm counts at the right of each rug row, four rug rows (streaming up and down, fresh robust up, fresh Gaussian up).

**EXPERIMENTS PERFORMED.** One full run (3.5 minutes, one process).

**RESULTS.** (Superseded by iteration 3.) Robust, fresh calibration, up alarm within 6 h: FOMC 14/16, 14/16, 15/16 against 30% to 31% of controls (p <= 0.000007); CPI 14/23, 16/23, 14/23 against 26% to 28% (p <= 0.0002). Median delay 23 to 34 minutes, minimum 14. Gaussian LR: alarms at 94% to 100% of releases but also in 65% to 78% of control windows; median delay 0 to 2 minutes.

**IMPROVEMENTS.** The Gaussian comparison is now informative, the FOMC and CPI contrasts are estimated on all events, and the blindness of the streaming protocol is reported as a finding rather than hidden in low rates.

**REMAINING ISSUES.** Review finding M6 (coordinator): the script filled missing returns with zeros before the intraday profile (`intraday_scale(r.fillna(0.0))`), and the profile's old default used future days for the first 20 days.

**NEXT ITERATION.** Fix the data handling and rerun.

---

## Iteration 3. Review fix: no zero-filled returns, past-only intraday profile, 20-day warm-up

**OBJECTIVE.** Apply review finding M6 and `intraday_scale`'s new default (`start="past"`) without changing the study design.

**PROBLEM IDENTIFIED.** Zero-filled missing minutes would look like evidence for a volatility fall (there are none in this Binance sample apart from the first return, but the code must not rely on that). With `start="past"` the profile is NaN for the first 10 days and uses fewer than 20 days until day 21, and Binance data start on 1 January 2024, so no earlier warm-up can be loaded.

**HYPOTHESIS.** Treating 1 to 20 January 2024 as warm-up and starting the monitored stream on 21 January 2024 (as E6b does) removes both problems at the cost of one event (the CPI release of 11 January 2024). Fresh-calibration results for the remaining events should be unchanged (the profile from day 21 on is the same 20-day rolling median as before); streaming results change because the stream starts later.

**CHANGES IMPLEMENTED.** `load_coin` passes returns with their NaN to `intraday_scale` (default `start="past"`), checks that the profile is finite for every monitored minute, and keeps the stream from `MONITOR_START = 2024-01-21`; the monitor receives NaN as NaN (it skips them). `read_events` excludes events inside the warm-up and writes them to `e12_events_excluded.csv`; control windows start on 21 January 2024. Docstring and report updated; the monitored period is 711 days.

**EXPERIMENTS PERFORMED.** One full run (2 minutes), `python -m pytest -q` (164 tests pass).

**RESULTS.** 38 events (16 FOMC, 22 CPI), 134 FOMC and 299 CPI control windows. Fresh-calibration rows for the remaining events are identical to iteration 2 (per-event delays unchanged); the control rates moved by at most one point. Streaming, robust, up alarm within 6 h: BTC FOMC 9/16 against 14% (p = 0.0001, Holm 0.0002, primary), CPI 5/22 against 11% (p = 0.11, primary); ETH 8/16 (p = 0.0002) and 7/22 (p = 0.003); SOL 9/16 (p = 0.0005) and 9/22 (p = 0.00004). Fresh calibration, robust: FOMC 14/16, 14/16, 15/16 against 29% to 31% (p <= 0.00001); CPI 13/22, 16/22, 14/22 against 25% to 27% (p <= 0.0004); Holm over six at most 0.0004. Delay: median 23 to 34 minutes (fresh), minimum 10 to 14. Pre-release up alarms 6% to 19% of events against 18% to 26% of controls. Down alarms in 6 to 48 h: no excess (p >= 0.17). Gaussian LR (fresh): 94% to 100% of releases, 64% to 77% of control windows, median delay 0 to 2 minutes, alarm on a return of 8 to 9 scale units, post / baseline 0.97 to 1.02 on control anchors.

**IMPROVEMENTS.** No zero-filled returns anywhere in the pipeline; every value of the intraday profile uses past days only; the excluded event is listed.

**REMAINING ISSUES.** (1) 16 and 22 events; Wilson intervals of about plus or minus 20 points. (2) Control windows at 08:30 ET contain other scheduled releases (weekly jobless claims on Thursdays, producer prices, payrolls), not verified or excluded, so the CPI contrast is conservative. (3) The fresh-calibration design was chosen after seeing iteration 1 and is marked secondary. (4) Q12c cannot be answered well under a protocol whose recalibration absorbs the burst; it needs either a shorter recalibration or a decrease detector that keeps the pre-release baseline. (5) Direct fetches of the official hosts remain blocked; verification is through search results that index those pages.

**NEXT ITERATION.** Not started: (a) a shorter or warm recalibration after alarms (a change in `monitor.py`, not owned here) and a rerun of the streaming design; (b) extend the control definition with verified dates of the other 08:30 ET releases; (c) add XRP and DOGE (the loaders exist) once the protocol change is decided, to avoid repeating runs.

---

## Iteration 4. Review round 2a fixes: within-hour counts, matched control rates, band-matched benchmark

**OBJECTIVE.** Apply findings C1, M2, M3, m10 and m11 of `docs/review/independent_review_round2a.md` without changing the frozen settings or the primary family.

**PROBLEM IDENTIFIED.** (1) C1: the Reading quoted the 6-hour counts as "within the hour". (2) M2: the control rates and the pooled Fisher test used the union of all events' control windows, which weights weekdays by the number of control days rather than by the events (FOMC: 65 Wednesdays and 69 Thursdays for 15 Wednesday statements and one Thursday statement); a 6-to-48-hour window that starts on a Thursday or Friday runs into the weekend, so the down-alarm base rate was overstated. (3) M3: the Gaussian monitor runs at band 0 and the robust monitor at band 0.25, so the rate-ratio contrast mixed the band with truncation. (4) m10 and m11: two sentences that the data did not support.

**HYPOTHESIS.** The stratified permutation test already compares each event with its own matched controls, so its p-values do not move; the matched rates change little for the up outcome and a lot for the down outcome; at equal band the robust monitor keeps a lower control rate than the Gaussian monitor.

**CHANGES IMPLEMENTED.** `stratified_tests` now reports the matched control rate (per event, the share of its own controls with the outcome, averaged over events), the expected count from matched controls, the matched rate ratio, a lower-tail permutation p-value (`p_stratified_permutation_fewer`) and an exact stratified Fisher test (`stratified_fisher_tail`, exact Cochran-Mantel-Haenszel with one stratum per weekday and clock time). The union rate and the pooled Fisher test stay as `union_*_superseded` and `p_fisher_pooled_superseded`. New detector "robust e-detector, band 0" in both designs. New secondary outcome `up_post_1h` (first up alarm at most 60 minutes after the anchor) with matched control rates in `e12_delay.csv` and `e12_event_windows.csv`. `weekend_flags` checks that every event and its controls agree on whether the pre, post and follow-up windows touch a weekend (they do; flags in `e12_events_checked.csv`). E12.md: headline table, all control rates matched, Q12c rewritten, band-0 rows in Q12d, wording of m10 and m11, revision note.

**EXPERIMENTS PERFORMED.** One full run (2.3 minutes, one process). Every column and row of the previous CSVs is reproduced exactly (only new columns and the band-0 rows were added; the figure is byte-identical), and the stratified permutation p-values are unchanged to the last digit.

**RESULTS.** Up alarm within 6 hours, robust, band 0.25: matched control rates 11% to 17% (streaming) and 23% to 32% (fresh), at most 4.2 points from the union rates; conclusions of Q12a unchanged. Within 1 hour: streaming FOMC 9, 6, 8 of 16 (9/9, 6/8, 8/9 when awake), fresh FOMC 14, 13, 11 of 16, CPI 8, 11, 11 of 22, against matched control rates of 1% to 8%. Down alarm in 6 to 48 hours after FOMC: fresh 63%, 56%, 56% against matched 49%, 43%, 42% (union had 71%, 66%, 65%; p = 0.20, 0.21, 0.17), streaming 69%, 63%, 69% against 45%, 43%, 41% (p = 0.04, 0.09, 0.02, Holm over six at least 0.15); "because the base rate is high" withdrawn. Down-alarm rates of control windows: 70% to 94% for Thursday and Friday anchors against 38% to 57% for Tuesday and Wednesday anchors. Band 0, fresh calibration: robust control rates 39% to 51% against 62% to 76% for the Gaussian monitor; rate ratios 1.7 to 2.4 against 1.3 to 1.6 (2.5 to 3.0 at band 0.25); about a third of the control-rate gap is the band. Robust minimum delay 8 to 11 minutes at band 0. Pre-release: fewer alarms on event days in 12 of 12 cells, one-sided p for "fewer" 0.11 to 0.65.

**IMPROVEMENTS.** Every control rate in E12 is matched on weekday and clock time; the headline numbers for both designs and both windows are in one table at the top of E12.md; the Gaussian comparison is like for like.

**REMAINING ISSUES.** The band-0 variant, the 1-hour outcome and the stratified Fisher test were added with the results known and are secondary. The robust and Gaussian monitors still differ in their betting rule as well as in truncation. The figure's control curve and the "controls" rows of the per-alarm table pool weekdays (descriptive only). The items of iteration 3 remain.

**NEXT ITERATION.** As in iteration 3.
