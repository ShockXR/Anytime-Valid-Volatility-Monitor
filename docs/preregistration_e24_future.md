# Pre-registration E24-future: the release test on releases after E23's data, with the guarded P1 or the rank monitor

*Written and committed before any Binance return after 24 October 2026 23:59 UTC was downloaded or loaded by anyone working on this project (section 2). The commit time of this file is the record. Anything decided after this commit will be marked post hoc in the report. The test cannot run before its end date, **1 June 2027 00:00 UTC** (section 5): `experiments/e24_blind_guard.py future` refuses to run before that date, on data that end before it, and while any release in the period is unverified.*

## 1. Why

E23's pre-registered primary test (P1, BTC, FOMC and CPI releases of 2026 pooled) failed: 2 of 16 releases with an up alarm within 6 hours against a matched control rate of 5.9%, p = 0.25. P1 was blind (unable to raise an increase alarm) or recalibrating at 5 of the 16 BTC releases, and E24 traced the blind and near-blind BTC windows to quiet weekend stretches in the part of the calibration window that fixes the truncation level, three of them at the scale of one price tick (`docs/report_parts/E24.md`). The rank monitor S2, which is never blind, flagged 4 of 6 FOMC statements on BTC (secondary family F23e). E24 proved a guard that removes blind windows without losing the per-day bound (`docs/theory_blind_guard.md`), and its post hoc reruns on 2024 to 2026 are exploratory only, because those data were seen. A confirmatory answer needs releases that nobody has looked at. This document fixes, before those data exist, which monitor and which test will be used, how many releases are needed, and when the test will be run: **the rank monitor S2, one stratified test over BTC's FOMC statements and CPI releases pooled, 11 planned releases, end date 1 June 2027** (section 5).

## 2. Freshness check: no committed work used crypto returns after 24 October 2026

Done on 3 November 2026, before this commit, without loading any market data.

* E23 (`experiments/e23_fresh_releases.py`) is the only committed work that loaded Binance bars of 2026. Its stream ends at 25 October 2026 00:00 UTC, exclusive (`results/e23_data_checks.csv`, column `stream_end_exclusive`): the mirror's October 2026 files stopped at 24 October 23:59, and it had no November file. The October files on this machine (`data/raw/crypto/binance/*/..._2026-10.csv`) still end there.
* Every other committed crypto load ends at 2026-01-01 (exclusive), as listed in `docs/preregistration_e23.md`, section 2; the git history since then adds no crypto load past that date except E23's, and E24's own runs (parts `synthetic`, `real`, `e19blind`, `future_gate`) read only data that end by 25 October 2026 00:00 UTC; the parts `theory`, `selection` and `power` read no market data. The Bitstamp BTC/USD file on disk runs to 2 November 2026, but every committed reader cuts it at 2026-01-01.
* The dates after 24 October 2026 that appear in committed results are Nifty 50 dates (E7, E19, the live Nifty study), not crypto.

So the fresh sample is all Binance returns of BTCUSDT, ETHUSDT and SOLUSDT from **25 October 2026 00:00 UTC** on. The 20 days before (5 to 24 October 2026, known data, used by E23) only warm up the intraday profile.

## 3. Events

`data/events/us_macro_events_e24_future.csv`, the columns of E23's file. Built on 3 November 2026 only from dates that web search tied to the official schedules; the official hosts (`www.federalreserve.gov`, `www.bls.gov`) could not be fetched from this environment, so, as for E12 and E23, verification used search results that index the official pages. Each row says what was found.

* **Released:** the FOMC statement of 28 October 2026, 14:00 ET (18:00 UTC), the row of E23's file (E23 excluded it because its data ended before it).
* **Scheduled, from the official schedules:** CPI releases at 08:30 ET on 10 November and 10 December 2026, and 13 January, 11 February, 10 March, 13 April, 12 May, 10 June, 14 July, 11 August, 14 September, 14 October, 10 November and 10 December 2027 (BLS schedule of releases for the Consumer Price Index; the 10 November 2026 date also on the BLS CPI home page). FOMC statements at 14:00 ET on 9 December 2026 and on the second day of each 2027 meeting of the Federal Reserve's tentative schedule (press release of 5 September 2025, `monetary20250905a`): 27 January, 17 March, 28 April, 9 June, 28 July, 15 September, 27 October and 8 December 2027. Third-party calendars list projected 2027 CPI dates that differ; only the BLS listing is used.
* **A 2027 meeting without a statement.** On 28 October 2026 the Fed chair said that one of the eight 2027 meetings is expected to be a strategy meeting with no policy vote and no statement, without saying which (news reports cited in the file). E12's rule (the statement on the second day of each scheduled meeting) then gives seven 2027 statements. For planning only, section 5 drops the earliest (27 January 2027), which gives the fewest statements for every end date.
* **Rule for updating the file before the test runs (the only change allowed after this commit).** After each release, change its `status` from `scheduled` to `released` and record the verification (official page or search results indexing it, access date). If a release moves, set `rescheduled`, put the new date in `date` and the old one in `original_date`. If a CPI release is cancelled (as in November 2025), set `cancelled`. If an FOMC meeting has no statement, set `no_statement`. An unscheduled FOMC statement is not added (E12's rule). The future part refuses to run while any row in the period still says `scheduled`, or while any month of the period has no CPI row.

## 4. Copied from E12, E22 and E23 without change

By importing `experiments/e12_crypto_events.py`, `experiments/e22_guaranteed_events.py` and `experiments/e23_fresh_releases.py`. Every (coin, arm) run is `E22.job` itself, with E12's period and event file switched by E23's `_set_period` (E23's gate showed that the switch, set to E12's own period, reproduces E22 exactly).

* **Coins and data:** Binance spot one-minute bars, BTCUSDT, ETHUSDT, SOLUSDT, from the mirror read by `avmonitor.data.load_binance`. Loaded from **5 October 2026** (20 days of intraday-profile warm-up), **monitored from 25 October 2026 00:00 UTC to the end date of section 5, exclusive.** Every monthly file from October 2026 to the month before the end date is downloaded; a file cached before its month was complete (October 2026) is downloaded again. The stream is cut at the end date even if the data go further. Data checks as in E23 (rows, duplicates, missing minutes and their runs, non-positive prices, NaN returns, SHA-256 of every file); missing returns stay NaN, never zero-filled.
* **Inclusion:** a release counts if it is at or after the stream start and its 6-hour window ends at or before the end date (E23's rule). **Controls:** E12's `control_anchors`: same weekday and US Eastern clock time, at least 2 days from every event day, anchor at least 1440 + 360 minutes after the stream start, 48-hour follow-up inside the period.
* **Design and outcome:** E12's streaming design; outcome an up alarm within [0, 6 h) of the release (`up_post`); a release at which the monitor is calibrating or blind counts as a release without one.
* **Test:** E12's exact stratified permutation test, one-sided (each release its own stratum with its matched control windows; T = number of releases with the outcome, Poisson-binomial with p_i = (y_i + k_i) / (1 + n_i)), through E23's `stratified_test` (pooled design) or E12's per-type test (one type), which agree to the last digit (E23's gate). Fisher cross-check reported, not used for decisions.

**Arms** (c = 3, calibration 1440 minutes, band b = 0.25, as in E22 and E23):

| Arm | Role | Constructor |
|---|---|---|
| G1 | candidate primary; secondary (F24c) after section 5 | `GuardedPeriodicVolatilityMonitor(threshold=36000, c=3, cal_len=1440, null_band=0.25, baseline="cs", delta=0.01, period=1440, offset=o, guard_tries=3, guard_eta=0.25, guard_fallback="rank", rank_seed=numpy.random.default_rng([2407]))` (`src/avmonitor/guard.py`): P1 with the retry guard of `docs/theory_blind_guard.md`; budgets δ/2, δ/4, δ/4 (halving); a window is discarded when τ² ≤ 1.25 m↑; after three such windows, E14's awake rank monitor for the rest of the segment; restart at 00:00 UTC (o checked by `E22.utc_day_offset`) |
| **S2** | **primary (chosen in section 5)** | E22's S2: E14's awake rank monitor, `ConformalMonitor(threshold=28800, cal_len=120, family="lr", band=0.25, seed=14000, restart_len=0)` |
| P1 | exploratory | E22's and E23's primary, `PeriodicVolatilityMonitor(threshold=36000, c=3, cal_len=1440, null_band=0.25, baseline="cs", delta=0.01, period=1440, offset=o)` |

**Guarantees** (each under its own hypothesis, which, as in E22 and E23, is not expected to hold on most real days; alarm rates on these data are therefore not false alarm rates). G1: for every period (a segment within one UTC day), P(an alarm in it) ≤ δ + m/A = 0.01 + 1440/36,000 = 5%, marginally, under the ±25% tolerance relative to the kept window's sequence part to the end of that day, or, in a segment that falls back to the rank monitor, exchangeability of the returns after the fallback given the past (Theorems 2 and 3 of `docs/theory_blind_guard.md`); long-run 0.0404 expected alarms per monitored day. S2: P(an alarm within the first 1440 steps after a restart) ≤ 1440/28,800 = 5% under exchangeability (E14). P1: as G1, without the guard.

## 5. Power, the primary test and the end date

Planning before any data after 24 October 2026 exist, from the BTC rates of the two candidate monitors in 2024 to 2026 (`results/e24_future_planning_inputs.csv`). These are **planning inputs** taken from data that were seen; G1's are from E24's post hoc rerun. Using them to plan a test on new data is sample-size planning, not data snooping: no data of the test period are touched.

| Monitor | Release type | 2024 to 2025 (E22's design) | 2026 (E23's design) | Pooled release rate | Pooled matched control rate |
|---|---|---|---|---|---|
| G1 | FOMC | 9 / 16 | 3 / 6 | 12 / 22 = 54.5% | 10.7% |
| G1 | CPI | 4 / 22 | 1 / 10 | 5 / 32 = 15.6% | 6.6% |
| S2 | FOMC | 13 / 16 | 4 / 6 | 17 / 22 = 77.3% | 16.5% |
| S2 | CPI | 10 / 22 | 2 / 10 | 12 / 32 = 37.5% | 12.5% |

(The pooled control rate is the mean of the per-release matched control rates over the releases of both periods.)

**Model.** As in E23 (`E23.simulate_power`, unchanged): each release has the outcome independently with probability p1 of its type and monitor, each matched control window with p0; releases of one weekday-and-time group share their control windows (control counts from `E12.control_anchors` on a gap-free minute index from 25 October 2026 to the end date); the test of section 4 at one-sided level 0.05. The calendar is the event file's released and scheduled rows minus the statement of 27 January 2027 (section 3: the earliest 2027 meeting is taken as the one without a statement, which gives the fewest statements for every end date); 2028 dates, where the grid needs them, are planning placeholders (2027 date plus 364 days, two FOMC meetings dropped), never event dates. Seed 24,000; common random numbers across designs. Reproduced by `nice -n 10 python experiments/e24_blind_guard.py power --procs=2` (`results/e24_future_power_grid.csv`, `e24_future_candidates.csv`, `e24_future_power.csv`, `e24_future_power_calendar.csv`, `e24_future_power_sensitivity.csv`; figure `results/figures/e24_future_power.png`).

**Decision rule** (in the docstring of `part_power`, committed in `09866ec` before the part was first run). Four candidates: monitor G1 or S2, design BTC FOMC alone or BTC FOMC and CPI pooled. For each, E* is the first end date on the grid of first days of a month, 1 January 2027 to 1 January 2029 (00:00 UTC, exclusive), at which the power under that monitor's pooled 2024 to 2026 rates is at least 0.80 (5000 replications), confirmed with 20,000 (else the next month). A candidate is eligible if its simulated size under the null (release rate = control rate) at E* is at most 0.055. The primary test is the eligible candidate with the earliest E*; ties: higher power, then the pooled design, then G1. The end date is that E*.

**Result** [SYNTHETIC; planning inputs REAL]:

| Candidate | E* | Releases planned (FOMC + CPI) | Power at E* (s.e.) | Size under the null at E* |
|---|---|---|---|---|
| G1, FOMC alone | 1 December 2027 | 8 + 13 | 0.831 (0.003) | 0.046 |
| G1, pooled | 1 November 2027 | 8 + 12 | 0.815 (0.003) | 0.040 |
| S2, FOMC alone | 1 July 2027 | 5 + 8 | 0.828 (0.003) | 0.041 |
| **S2, pooled** | **1 June 2027** | **4 + 7** | **0.822 (0.003)** | **0.031** |

**Choice: monitor S2, design BTC FOMC and CPI pooled, end date 1 June 2027 00:00 UTC (exclusive).** The sample for power 0.8 at the planning rates is **11 releases**: the FOMC statements of 28 October and 9 December 2026 and 17 March and 28 April 2027, and the CPI releases of 10 November and 10 December 2026 and 13 January, 11 February, 10 March, 13 April and 12 May 2027 (12 if the meeting of 26 to 27 January 2027 issues a statement). The period is about seven months of data, less than a year. The end date is fixed here and in `experiments/e24_blind_guard.py` (`FUTURE_END`, `FUTURE_PRIMARY_ARM`, `FUTURE_PRIMARY_DESIGN`).

All candidates at the chosen end date (20,000 replications; `e24_future_power.csv`):

| Scenario | G1 FOMC alone | G1 pooled | S2 FOMC alone | **S2 pooled** |
|---|---|---|---|---|
| pooled 2024 to 2026 rates (planning) | 0.506 | 0.527 | 0.696 | **0.822** |
| 2026 rates only (E23's design) | 0.454 | 0.329 | 0.687 | **0.596** |
| half of the pooled excess over the control rate | 0.232 | 0.224 | 0.261 | **0.369** |
| null (release rate = control rate) | 0.048 | 0.034 | 0.032 | **0.031** |

**What this power does not cover, stated now.** The planning rates include 2024 and 2025, when every monitor's release rates were higher than in 2026. If the 2026 rates are the truth, the chosen test has power 0.60 at the end date, and 0.80 would need data to 1 October 2027 (0.84 there; `e24_future_power_sensitivity.csv`, descriptive, computed after the first run and not part of the rule); if the true excess is half of the pooled one, power is 0.37, and 0.80 would need data to about August 2028, more than a year and partly on placeholder dates. The rule fixed in advance takes the earliest date with power 0.8 at the planning rates; a failure of the primary test will therefore be weak evidence against a smaller effect (section 9). The rejection rate under the null model is 3.1% at nominal 5% (Monte Carlo s.e. 0.12 points): the pooled test is conservative here, as in E23.

## 6. The primary test, the decision rule and the secondary families

All tests: outcome `up_post`, all releases included, E12's stratified permutation test, one-sided; Holm's step-down procedure within each family at 5%; no adjustment across families. Implemented in `future_families` and `future_pipeline` of `experiments/e24_blind_guard.py`.

| Family | Arm | Cells | Tests | Status |
|---|---|---|---|---|
| **F24a** | **S2** | **BTC, FOMC and CPI releases pooled** | **1** | **primary** |
| F24b | S2 | ETH and SOL, pooled | 2 | secondary (replications; the coins move together, so they are not independent evidence) |
| F24c | G1 | BTC, ETH and SOL, pooled | 3 | secondary (the guarded P1 of E24) |
| F24d | S2 | three coins, FOMC and CPI separately | 6 | secondary |

**Primary test and decision rule.** S2, BTC, all included FOMC statements and CPI releases from 25 October 2026 to 1 June 2027 in one stratified permutation test. **It passes if p ≤ 0.05.** Reported effect: releases with an up alarm within 6 hours out of the number included, the mean matched control rate, their difference and ratio, the 95% Wilson interval of the release rate, p, and the Fisher cross-check; the FOMC and CPI parts are shown separately as description (and tested in F24d).

**Exploratory, reported raw, no decision attached:** P1 (unguarded) in every cell; G1 per release type; the 1-hour outcome; the tests restricted to releases (and controls) at which the monitor was watching; E12's other outcomes; a comparison of each rate with its 2024 to 2025 and 2026 rates (description only).

## 7. Diagnostics

* **Per-day alarm share** against the 5% bound and alarms per monitored day against the long-run bound, per arm and coin (`E22.period_counts`), with Wilson intervals that treat periods as independent (descriptive).
* **Releases watched** (`E22.watched_table`): monitoring at the release minute with an increase side able to alarm; for G1 a rank-fallback segment counts as watched.
* **Guard bookkeeping for G1:** windows tried, discarded and kept, rank fallbacks, extra calibration minutes, and the τ-part and sequence-part scales of every window with τ²/m↑ ≤ 1.25 (E24's mechanism check).
* **Median delay** from release to the first up alarm among releases with one.

## 8. Refusals and gates (in code; the future part stops if any fails)

1. `FUTURE_END` in `experiments/e24_blind_guard.py` must be the end date of section 5, and the run date must be on or after it.
2. Every row of the event file inside the period must have status `released`, `rescheduled`, `cancelled` or `no_statement`, and every calendar month of the period after October 2026 must have a CPI row.
3. Every monthly file up to the end date must exist in the mirror, and the last bar present for all three coins must be at or after the end date minus one minute.
4. **Pipeline gate** (part `future_gate`, run before this commit): the future pipeline, set to E23's period and events, reproduces E23's committed P1 and S2 test rows exactly (192 rows: every outcome, sample and cell; 14 columns; 0 mismatches; `results/e24_future_gate.csv`).
5. The number of included releases and the control counts of the primary design are reported against the plan of section 5; if fewer releases are included than planned (a release cancelled, or a statement-free meeting), the test is run on those included and the power at the actual count is reported, recomputed with the same code.

## 9. Reading fixed in advance

* **If the primary test passes (p ≤ 0.05):** "On FOMC statements and CPI releases after 24 October 2026, which nobody in the project had seen when the design was fixed, E14's rank monitor (a proved bound of 5% on the probability of an alarm within the first day after each restart, under exchangeability) raised up alarms within six hours more often on BTC than at matched control times." The release contrast of E22 (2024 to 2025) and of E23's secondary family F23e (2026 FOMC statements) is then confirmed on fresh data for the pooled hypothesis. It does not say which release type carries the effect (F24d and the description do), that the bound holds on these data (the exchangeability hypothesis fails on most real days, as E14 and E19 found), or anything about the guarded P1 (F24c).
* **If it fails (p > 0.05):** the release contrast is not confirmed on fresh releases. With planned power 0.82 at the pooled 2024 to 2026 rates, a failure would happen with probability about 0.18 if those rates were true, about 0.40 at the 2026 rates and about 0.63 at half of the pooled excess; a failure therefore says that the contrast after October 2026 is smaller than in 2024 to 2026 or absent, and does not show that there is none. The release results of the project then rest only on data that were seen when the designs were fixed (E12, E22) or on secondary families (E23), and the report and paper must say so.
* **Secondary families.** F24b says whether the pooled contrast also appears on ETH and SOL (not independent of BTC). F24c says whether the guarded P1, which carries the per-day bound under the tolerance, shows the contrast on fresh releases; with the planning rates its power at this end date is about 0.53 on BTC (`e24_future_power.csv`), so a failure there is weak evidence. F24d tests the FOMC and CPI cells separately, with low power. None of them can replace or rescue the primary test.
* **If fewer releases are included than planned** (a cancelled CPI release, a postponed meeting), the test is run on those included; the power at the actual count will be reported, recomputed with the same code and seed.
* **If more are included** (a statement on 27 January 2027), the test is run on all included releases, as E12's rule says.

## 10. Outputs

`experiments/e24_blind_guard.py` (parts `power` and `future`); `data/events/us_macro_events_e24_future.csv`; `results/e24_future_planning_inputs.csv`, `e24_future_power_grid.csv`, `e24_future_candidates.csv`, `e24_future_power.csv`, `e24_future_power_calendar.csv` and `e24_future_power_sensitivity.csv`, `e24_future_gate.csv` (now); `results/e24_future_*.csv` from the future part (data files and checks, events used and excluded, tests with families, windows, monitor summary, watched, alarms); `results/figures/e24_future_power.png`; `tests/test_guard.py` (the refusals and the event file, synthetic, no market data); `docs/report_parts/E24.md`. Runs use `nice -n 10` and at most 2 worker processes. Random numbers: G1's rank fallback seed 2407, S2's seed 14,000, the power simulation seed 24,000.
