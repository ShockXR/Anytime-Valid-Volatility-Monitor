# Pre-registration E23: E22's release test on fresh releases, the FOMC statements and CPI releases of 2026

*Written and committed, together with the event file `data/events/us_macro_events_2026.csv` and nothing else, before any Binance return after 31 December 2025 was downloaded or loaded by anyone working on this project (section 2). The commit time of this file is the record. Anything decided after this commit will be marked post hoc in the report.*

## 1. Why

Every real-data result on FOMC statements and CPI releases so far (E12, E14 Q14g, E22) used the releases of 2024 and 2025. E22 reran E12's pre-specified primary test with monitors that carry a proved per-day bound, but it was designed after E12's FOMC result was known, on the same data (referee round 2, F15; E22's own caveat). E23 asks E22's question, with E22's primary monitor, on releases that nobody in this project has looked at: the FOMC statements and CPI releases of 2026, from January to the end of the data (expected: October 2026). This turns the one pre-specified real-data test into a confirmatory one. A negative result will be reported as fully as a positive one.

## 2. Freshness check: no committed experiment or result used crypto returns after 2025-12-31

Done on 3 November 2026, before this commit, without loading any market data.

* **Experiments and scripts.** Every committed call that loads crypto bars (`load_binance`, `load_bitstamp_csv`, `download_binance`, and the period tables that feed them) ends at `2026-01-01`, exclusive: E5 and E8 (Bitstamp BTC, June or July to December 2025), E6a, E6b (`PERIODS`: July to December 2025 and 21 January 2024 to 30 June 2025), E7 (its Bitstamp BTC comparison, July to December 2025), E10, E11f, E12 (`START, END = "2024-01-01", "2026-01-01"`), E13 (to 2025-07-01), E14, E15, E17, E19 and E21 (through E6's `prepare`, July to December 2025), E18 (E12's stream; its crypto landmark windows end by 12 October 2025), E20, E22 (through E12), the six studies in `alternatives/` (crypto loads end at 2026-01-01; `risk_switch` freezes `crypto_load` to 2024-01-01 to 2026-01-01), `scripts/download_data.py` (crypto to 2026-01-01) and the demo notebook (October 2025). The git history of `experiments/`, `alternatives/`, `scripts/`, `src/`, `notebooks/` and `examples/` has no added line that loads crypto bars past 2026-01-01.
* **Results.** The only committed result files with 2026 dates are Nifty 50 files: `e07a_daily_iq_nifty.csv`, `e07a_frozen_runs.csv`, `e07b_alarms.csv`, `e07b_session_vol.csv`, `e07d_event_detail.csv`, `e07d_expiry_left_out.csv`, the Nifty rows of `e19_calm_days.csv` (column `dataset` = "Nifty 50" for all 205 rows dated 2026) and `results/live/` (a Nifty dry run). No integer column of any result CSV (in `results/` or `alternatives/*/results/`) holds an epoch timestamp in 2026 (seconds, milliseconds, microseconds or nanoseconds), and the history of the result files has no added crypto row dated 2026.
* **Raw data on this machine.** No Binance file after December 2025 is present (`data/raw/crypto/binance/*/` ends with the 2025-12 files; `data/MANIFEST.csv` lists none). The Bitstamp BTC/USD file in `data/raw/crypto/` was downloaded on 2 November 2026 and runs to that day, but every committed reader cuts it at 2026-01-01; E23 does not read it.
* **Event files.** `data/events/us_macro_events.csv` lists the CPI release of 13 January 2026 as released after E12's sample (documentation only, never used).
* **Nifty.** Nifty 50 bars of 2026 were used by E7, E7c, E7d, E19 and the live Nifty study; that is why E23 is crypto only.

So the fresh sample is all Binance returns of BTCUSDT, ETHUSDT and SOLUSDT from 1 January 2026 on. The December 2025 bars (known data) enter E23 only through the intraday profile's 20-day warm-up (section 4), never as monitored returns, release windows or control windows.

## 3. Events

`data/events/us_macro_events_2026.csv`, same columns as `data/events/us_macro_events.csv` plus `access_date`. Every date was checked by web search on 3 November 2026 against the official pages: the Federal Reserve's press release for each statement (`monetaryYYYYMMDDa.htm`) and meeting page, and the BLS news release archive page for each CPI release (`cpi_MMDDYYYY.htm`) or, where that page was not surfaced, BLS regional releases of the same date and the BLS current-release page. The official hosts could not be fetched from this environment (`www.federalreserve.gov` and `www.bls.gov` do not resolve here), so, as for E12, verification used search results that index the official pages. No date was left out for lack of verification.

* **FOMC statements (7), 14:00 US Eastern, second day of each scheduled meeting:** 28 January (EST, 19:00 UTC), 18 March, 29 April, 17 June, 29 July, 16 September, 28 October 2026 (EDT, 18:00 UTC). All Wednesdays. E12's rule is "the second day of each scheduled meeting" (notation votes left out). The search found no unscheduled FOMC meeting or intermeeting statement in 2026; had there been one, E12's rule would not include it.
* **CPI releases (10), 08:30 US Eastern:** 13 January (December 2025 CPI, EST), 13 February (January CPI, EST; **rescheduled** from Wednesday 11 February after the partial lapse in federal appropriations that began on 31 January 2026; BLS archive page `cpi_02132026.htm`, BLS revised release dates page, news reports; the original date rests on secondary sources), 11 March (EDT from 8 March), 10 April, 12 May, 10 June, 14 July, 12 August, 11 September, 14 October 2026 (12:30 UTC). Weekdays: Tuesday 3, Wednesday 4, Friday 3. The 11 September date rests on BLS regional releases dated 20260911 and the BLS current-release page (indexed as "2026 M08 Results", embargoed until 11 September 2026); the national archive page was not surfaced. No 2026 CPI release in the sample was cancelled or postponed except the one on 13 February. The October 2026 CPI is scheduled for 10 November 2026, after the data.
* UTC times come from the IANA time zone database (America/New_York); E12's reader recomputes them and stops if the file disagrees.

## 4. Copied from E12 and E22 without change

By importing `experiments/e12_crypto_events.py` and `experiments/e22_guaranteed_events.py`. Every (coin, arm) run is `E22.job` itself; the only change is the period and the event list, made by setting E12's module constants (`EVENTS`, `START`, `MONITOR_START`, `END`) and replacing `E12.read_events` by E12's reader followed by the inclusion rule below. A gate (section 8) checks that this switch, set to E12's own period, reproduces E22's committed P1 rows on BTC exactly.

* **Coins and data:** Binance spot one-minute bars, BTCUSDT, ETHUSDT, SOLUSDT, from the mirror read by `avmonitor.data.load_binance`. Loaded from **12 December 2025** (20 days that only warm up the intraday profile, as E12's 1 to 21 January 2024; the profile at a minute uses earlier days only), **monitored stream from 1 January 2026 00:00 UTC**, which leaves 12 days (calibration of 1440 minutes included) before the first release. The stream does not need to start in December: E12's design calibrates on the first 1440 minutes of the stream. **End:** the monthly files are downloaded from December 2025 on, month by month until the first month the mirror does not have; the stream ends at the last minute present for all three coins (exclusive end = that minute + 1 minute). Expected: 31 October 2026 23:59 UTC.
* **Data checks, reported:** rows read, duplicate timestamps (and whether duplicated rows are identical; `load_binance` keeps the first), missing minutes (count, and every run of them), non-positive prices, the last timestamp per coin, NaN returns in the monitored period. **Missing returns stay NaN, never zero-filled**; the monitors skip them (as in E12). SHA-256, size and row count of every file read are written to `results/e23_data_files.csv`.
* **Inclusion rule:** a release counts only if its 6-hour window [t, t + 6 h) ends at or before the data end. **Controls:** E12's `control_anchors` over the 2026 period, unchanged: every day of the monitored period with the release's weekday at the release's US Eastern clock time, at least 2 days from every event day (FOMC and CPI), the anchor at least 1440 + 360 minutes after the stream start and its 48-hour follow-up window inside the data, so that the 6-hour window of every control also ends before the data end.
* **Design and outcome:** E12's streaming design (one monitor over the whole stream: calibrate 1440 minutes, monitor, alarm, recalibrate). Outcome: an up alarm within [0, 6 h) of the release (`up_post`); a release at which the monitor is calibrating or blind counts as a release without one.
* **Arms (E22's constructors, `E22.make_arm`; c = 3, calibration 1440 minutes, null band 0.25):**

| Arm | Role | Constructor |
|---|---|---|
| **P1** | primary | `PeriodicVolatilityMonitor(threshold=36000, c=3, cal_len=1440, null_band=0.25, baseline="cs", delta=0.01, period=1440, offset=o)`, o such that every restart falls on 00:00 UTC (o = 0 here; `E22.utc_day_offset` checks it) |
| P0 | secondary | E12's plug-in, `VolatilityMonitor(threshold=28800, c=3, cal_len=1440, null_band=0.25)` |
| S2 | secondary | E14's awake rank monitor, `make_monitor("awake-B", 0.25, 1440, 28800, 36000, seed=14000)` |
| S1 | secondary | P1 without the restart, `VolatilityMonitor(threshold=36000, c=3, cal_len=1440, null_band=0.25, baseline="cs", delta=0.01)` |

The guarantee each arm carries, and its hypothesis, are those of E22's pre-registration (section 5 there): P1, δ + m / A = 5% for every day of a segment, marginally, under the ±25% tolerance of `docs/math.md` Theorem 6.4 to the end of that day (Corollary 1 of `docs/theory_per_period.md`); S1, 5% for the first day of a segment; S2, 5% for the first day after each restart under exchangeability; P0, nothing proved.

* **Test:** E12's exact stratified permutation test, one-sided (more releases with the outcome than expected). Each release is its own stratum: its window and its matched control windows, within which the release label is exchangeable under the null; T = number of releases with an up alarm is Poisson-binomial with p_i = (y_i + k_i) / (1 + n_i) (k_i of the release's n_i control windows with the outcome); p = P(T ≥ observed). The release type enters through the strata: FOMC and CPI releases have their own clock times and so their own control windows. Cross-check: E12's exact stratified Fisher test (one group per weekday and clock time), reported, not used for decisions.

## 5. Power, and the choice of the primary test

Planning before any 2026 return is seen, from E22's committed P1 rates on BTC in 2024 and 2025 (`results/e22_tests.csv`): FOMC 8 of 16 releases with an up alarm within 6 hours (50.0%) against a matched control rate of 10.9%; CPI 5 of 22 (22.7%) against 5.9%. Using the rates of a previous experiment to plan the next one is ordinary sample-size planning, not data snooping: it touches no 2026 data. The two candidate designs fixed beforehand (by the instructions for E23) were (a) BTC FOMC alone and (b) one stratified test over BTC FOMC and CPI releases together.

**Model.** The 2026 calendar of section 3 with an assumed data end of 1 November 2026 00:00 UTC: 7 FOMC releases (Wednesday 14:00 ET, 29 control windows each, shared), 10 CPI releases (Tuesday 08:30 ET, 29; Wednesday 08:30 ET, 29; Friday 08:30 ET, 40; shared within each group), from `E12.control_anchors` on a gap-free minute index (`results/e23_power_calendar.csv`). Each release has the outcome independently with probability p1 of its type, each control window with p0 of its type; releases of one weekday-and-time group share their control windows, as in the data. Test as in section 4 at one-sided level 0.05. 20,000 replications, seed 23,000, common random numbers across designs. Reproduced by `nice -n 10 python experiments/e23_fresh_releases.py power` (`results/e23_power.csv`).

| Scenario | (a) BTC FOMC alone, 7 releases | (b) BTC FOMC and CPI together, 17 releases |
|---|---|---|
| **E22's rates** (FOMC 50% vs 10.9%, CPI 22.7% vs 5.9%) | power **0.690** (Monte Carlo s.e. 0.003) | power **0.818** (0.003) |
| Half of E22's excess over the control rate (FOMC 30.5%, CPI 14.3%) | 0.330 | 0.405 |
| E22's FOMC rates, no CPI effect (CPI 5.9% vs 5.9%) | 0.690 | 0.534 |
| Null (release rate = control rate) | 0.058 (s.e. 0.0017) | 0.039 (s.e. 0.0014) |

**Choice: design (b), the pooled test, power 0.82 at E22's rates.** It has the better power at the planning rates, as the rule for choosing required. Its cost, stated now: if the CPI effect is absent and only the FOMC effect is as in E22, the pooled test has power 0.53 against 0.69 for FOMC alone. E22's rates are probably optimistic (the FOMC result is the reason this question is asked), so the half-excess line is the more sober planning value: about 0.40 for the pooled test.

**A property of E12's test found while planning (not a 2026 result).** Under the null model above, FOMC alone rejects in 5.8% of replications at nominal 5% (Monte Carlo s.e. 0.17 points), because the seven releases share their 29 control windows while the test treats their strata as independent; the pooled test rejects in 3.9%. The decision rule below keeps E12's test unchanged; the excess size of the FOMC-only cells is recorded here and will be repeated in the report.

**Reference on the known data** (2024 and 2025, E22's committed P1 BTC windows, same function): FOMC alone 8 of 16, p = 0.00017 (E22's number); pooled 13 of 38 against a mean matched control rate of 8.0%, p = 6.5 × 10⁻⁶.

## 6. The primary test, the decision rule and the secondary families

All tests: outcome `up_post`, all releases, E12's stratified permutation test, one-sided; Holm's step-down procedure within each family at 5%; no adjustment across families.

| Family | Arm | Cells | Tests | Status |
|---|---|---|---|---|
| **F23a** | **P1** | **BTC, FOMC and CPI releases together (pooled)** | **1** | **primary** |
| F23b | P1 | ETH and SOL, pooled | 2 | secondary (replications; the coins move together, so they are not independent evidence and are never pooled with BTC) |
| F23c | P1 | three coins, FOMC and CPI separately | 6 | secondary (E22's cells on fresh releases) |
| F23d | P0 | three coins, FOMC and CPI separately | 6 | secondary |
| F23e | S2 | three coins, FOMC and CPI separately | 6 | secondary |
| F23f | S1 | three coins, FOMC and CPI separately | 6 | secondary |

**Primary test and decision rule.** P1, BTC, all included FOMC statements and CPI releases of 2026 in one stratified permutation test. **It passes if p ≤ 0.05** (one test, no adjustment). Reported effect: releases with an up alarm within 6 hours out of the number included, the mean matched control rate, their difference and ratio, the 95% Wilson interval of the release rate, p, and the Fisher cross-check; the FOMC and CPI parts are shown separately as description (and tested in F23c).

**Exploratory, reported raw, no decision attached:** the pooled test for P0, S2 and S1; the 1-hour outcome (`up_post_1h`), per type and pooled; the per-type tests restricted to releases (and controls) at which the monitor was watching; E12's other outcomes (`up_pre`, `down_follow`, `any_post`), which `E12.stratified_tests` computes anyway; and a comparison of each 2026 rate with E22's 2024 to 2025 rate for the same cell (description only).

## 7. Diagnostics

* **Per-day alarm share against the 5% bound,** per arm and coin (`E22.period_counts`, `E22.first_day`): for P1 the share of monitored periods (a segment within one UTC day) with an alarm, the quantity Corollary 1 bounds by 5%; for P0, S1 and S2 the share of segments with an alarm within their first 1440 monitored steps; alarms per UTC day and per monitored day against the long-run bounds (0.0404 per monitored day for P1 and S1, 0.05 for S2, nominal 0.05 for P0). Wilson intervals treat periods and segments as independent and are descriptive only. As in E22, the tolerance is not expected to hold on real days, so these are alarm rates, not false alarm rates.
* **Releases watched,** with E19's blind-segment criterion (`E22.watched_table`): monitoring at the release minute and, for P0, P1 and S1, the segment's increase side able to alarm (off when τ² ≤ m_up = U (1 + b)); for S2 watched means monitoring. Per arm, coin and type, with the controls' watched shares and E19's per-run blind summary.
* **Median delay** from release to the first up alarm, among releases with one: descriptive only.

## 8. Gates (run first; E23 stops if either fails)

1. **Period switch.** With E12's own period and event file, E23's switch must reproduce E22's committed P1 BTC rows exactly: every row of `e22_tests.csv` for P1 and BTC (all outcomes and samples, 14 columns), the P1 BTC rows of `e22_windows.csv`, and the alarm and period counts of `e22_monitor_summary.csv` (equality after reading at round-trip precision).
2. **Pooled function.** The pooled test reads each release's stratum (y, n, k) from E22's window rows; restricted to one release type it must give `E12.stratified_tests`' counts, expected values, permutation and Fisher p-values exactly, for every arm, coin, type and the outcomes `up_post` and `up_post_1h`.

## 9. Reading fixed in advance

* **If the primary test passes (p ≤ 0.05):** "On 2026 releases that were not seen when the design was fixed, a monitor with a proved bound of 5% on the probability of a false alarm in each day of a segment (marginally, under a ±25% tolerance on the conditional truncated variance up to the end of the day) raised up alarms within six hours of FOMC statements and CPI releases on BTC more often than at matched control times." The E12 and E22 release contrast is then confirmed on fresh data for the pooled hypothesis. It does not by itself say which release type carries the effect (F23c and the description do), that the bound holds on these data, or that the alarms mark lasting changes.
* **If it fails (p > 0.05):** the release contrast is not confirmed on fresh releases. With planned power 0.82 under E22's rates, a failure would happen with probability about 0.18 if those rates were true, and about 0.60 if the true excess were half of E22's; a failure therefore says that the 2026 contrast is smaller than E22's or absent, and does not show that there is none. The E12 and E22 FOMC result then rests on data that were known when E22 was designed, and the report, README and paper must say that the only confirmatory test did not confirm it. The watched shares and the rates among watched releases will be reported to say descriptively whether blind or recalibrating time or fewer alarms at watched releases lie behind the failure; no test is attached to that split.
* **Secondary families.** F23b says whether the pooled contrast also appears on ETH and SOL (not independent of BTC). F23c says, on fresh releases, whether FOMC and CPI each carry the contrast under P1, cell by cell as in E22; with 7 FOMC statements and about 10 CPI releases per coin these cells have low power (section 5), so a failure there is weak evidence. F23d to F23f say whether the plug-in, the rank monitor and the confidence-sequence monitor without restart agree with P1. None of them can replace or rescue the primary test.
* **If fewer releases are included than planned** (an earlier data end), the test is run on those included; the power at the actual count will be reported, recomputed with the same code.

## 10. Outputs

`experiments/e23_fresh_releases.py` (parts `power`, `gate`, `run`); `results/e23_power.csv`, `e23_power_calendar.csv`, `e23_gate.csv`, `e23_data_files.csv`, `e23_data_checks.csv`, `e23_data_gaps.csv`, `e23_events_used.csv`, `e23_events_excluded.csv`, `e23_pooled_check.csv`, `e23_tests.csv`, `e23_windows.csv`, `e23_monitor_summary.csv`, `e23_watched.csv`, `e23_delay.csv`, `e23_alarms.csv`, `e23_e22_reference.csv`; `results/figures/e23_release_rates.png`; `tests/test_e23.py` (synthetic data only: event file schema, inclusion rule, controls exclude event days, the power function on a toy case, the pooled function against E12's per-type test); `docs/report_parts/E23.md`; `docs/iterations/E23.md`. All results on market data are **[REAL]**; tests and the power calculation are **[SYNTHETIC]** (the power calculation uses E22's real rates as inputs to a simulation). Runs use `nice -n 10` and at most 2 worker processes. No random numbers enter P0, P1 or S1; S2 uses E14's seed 14,000 for its tie-breaking; the power simulation uses seed 23,000.
