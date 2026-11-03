# Pre-registration E22: the FOMC and CPI test of E12 with monitors that carry a proved per-day bound

*Written and committed before any part of E22 was run on real data: before the plug-in arm was rerun for the reproduction check, and before any of the three configurations below was run on the release streams. The commit time of this file is the record. Anything decided after this commit will be marked post hoc in the report.*

## 1. Why

The referee's second report (`docs/review/referee_report_paper_round2.md`, finding J1) says that every real-data analysis of the paper uses the plug-in baseline, which has no finite-sample guarantee, including the only pre-specified primary test of the real-data work: E12's test of up alarms after FOMC statements and CPI releases on BTC (classified as primary, "stated", in E20). The abstract nevertheless says that the FOMC result shows "what the guarantee buys". E22 reruns E12's primary test with monitor configurations that carry a proved per-day false alarm bound. Only the monitor configuration changes; everything else is E12's, by importing E12's functions.

## 2. What has been seen before (disclosure)

* **E12** (plug-in, A = 28,800, band 0.25, streaming): FOMC on BTC 9 of 16 releases with an up alarm within 6 hours against a matched control rate of 14%, p = 0.00012 (Holm 0.00025); CPI on BTC 5 of 22 against 12%, p = 0.11. ETH and SOL: all four cells p at most 0.0032. The streaming plug-in was monitoring at about half of the releases.
* **E14 Q14g** ran E12's streaming design with the rank monitor without blind period (awake-B, configuration S2 below) and published its result (`results/e14g_e12_tests.csv`): BTC FOMC 13 of 16, p = 2.0e-7; BTC CPI 10 of 22, p = 0.00027; every awake-B cell p at most 0.00042. **S2 is therefore a rerun of a published analysis**, included so that all three configurations are reported on the same footing (same tables, families, per-day rates and watched shares); it adds no confirmatory evidence.
* **E19 Part 2** committed alarm lists of `VolatilityMonitor(baseline="cs")` with A = 36,000 and δ = 0.01 at bands 0, 0.25 and 0.5 (`results/e19_cs_events.csv`) for these coins over 21 January 2024 to 30 June 2025 and 1 July to 31 December 2025, on E6's input, which is built like E12's (same mirror, same 20-day intraday profile). Over the first of these periods S1 below may give the same alarm times. Those alarms were never related to the release calendar, and I have not done so before this commit.
* **P1** (the daily restart with the confidence-sequence baseline) has never been run on real data; TH6 ran it on synthetic streams only.

The question of E22 was chosen after E12's FOMC result was known. What is fixed here is the configuration, the test, the families and the reading; the data are not new (compare referee round 2, F15).

## 3. Copied from E12 without change

By importing `experiments/e12_crypto_events.py` (functions `read_events`, `load_coin`, `control_anchors`, `window_stats`, `alarm_frame`, `stratified_tests`, `poisson_binomial_tail`, `stratified_fisher_tail`, `holm`, `wilson`; constants `COINS`, `CAL`, `BAND`, `POST`, `HOUR`, `TYPES`):

* **Releases:** the 16 FOMC statements (14:00 US Eastern) and 22 CPI releases (08:30 US Eastern) of 2024 and 2025 used by E12, from `data/events/us_macro_events.csv` (the CPI release of 11 January 2024 excluded as in E12).
* **Coins and data:** Binance spot one-minute bars, BTCUSDT, ETHUSDT, SOLUSDT, 1 January 2024 to 31 December 2025; returns divided by the predictable intraday profile (`intraday_scale`, 20-day lookback, past days only); monitored stream from 21 January 2024 (711 days, 1,023,840 minutes per coin). The profile at minute t uses only earlier days, so it is known before t and the deseasonalised return u_t is a valid input for the theorems below; the hypotheses are statements about u.
* **Design:** E12's *streaming* design: one monitor runs over the whole stream (calibrate, monitor, alarm, recalibrate or restart). E12's fresh-calibration design is not a streaming monitor and is not rerun.
* **Outcome:** an up alarm within [0, 6 h) of the release (E12's `up_post`).
* **Controls:** for each release, every day of the monitored period with the same weekday at the same US Eastern clock time, at least 2 days from every event day (E12's `control_anchors`); matched control rate = mean over releases of the share of its own control windows with the outcome.
* **Test statistic and sidedness:** E12's exact stratified permutation test, one-sided (more event windows with the outcome than expected): within each release's stratum (its window and its control windows) the event label is exchangeable under the null, T = number of releases with an up alarm, p = P(T ≥ observed) under the Poisson-binomial law with p_i = stratum share. Cross-check: E12's exact stratified Fisher test (reported, not used for decisions).

## 4. Monitor configurations

Common to the confidence-sequence arms: c = 3, cal_len = 1440 minutes, null band b = 0.25 (E12's frozen band), δ = 0.01 split equally between the two bounds, τ from the first 25% of the calibration window (360 minutes), the empirical-Bernstein sequence on the remaining 1080 (`tau_frac=0.25`, `cs_c=0.75`, `cs_n_plan=1080`, the defaults), default bets (six increase bets 0.02 to 0.5, five decrease bets), Shiryaev-Roberts statistics, each side alarming at 2A. A = 1440 / (0.05 - δ) = 36,000, as in E8, E14 and E19, so that the per-day budget δ + m / A is E12's 5%. After every alarm the monitor recalibrates on the next 1440 minutes, as in E12.

| Arm | Role | Constructor |
|---|---|---|
| **P1** | primary | `PeriodicVolatilityMonitor(threshold=36000, c=3, cal_len=1440, null_band=0.25, baseline="cs", delta=0.01, period=1440, offset=o)` |
| **S1** | secondary | `VolatilityMonitor(threshold=36000, c=3, cal_len=1440, null_band=0.25, baseline="cs", delta=0.01)` |
| **S2** | secondary | E14's awake rank monitor: `make_monitor("awake-B", 0.25, 1440, 28800, 36000, seed=14000)` of `experiments/e14_running_tau.py`, that is `ConformalMonitor(threshold=28800, cal_len=120, family="lr", band=0.25, seed=14000, restart_len=0)` |
| **P0** | reference, reproduction | E12's plug-in: `VolatilityMonitor(threshold=28800, c=3, cal_len=1440, null_band=0.25)` |

**P1 runs E12's design as it stands.** `PeriodicVolatilityMonitor` with `baseline="cs"` calibrates and recalibrates exactly as `VolatilityMonitor(baseline="cs")` and, in addition, sets both statistics to 0 at every step t with (t - o) mod 1440 = 0 after the end of calibration. The offset o is chosen so that the restart falls on the first minute of each UTC day (o = 0 if the stream starts at 00:00 UTC; the script checks that every restart minute is 00:00 UTC). No substitute configuration is needed. One consequence, the same for releases and their controls up to the daylight saving offset: a 6-hour window that starts at 19:00 UTC (FOMC statements in Eastern standard time) contains a restart 5 hours in.

## 5. The guarantee each arm carries on these data

Notation: u_t the deseasonalised return, X_t = min(u_t², τ²), μ_t = E[X_t | F_{t-1}], m = 1440 (one UTC day). (H↑) and (H↓) are the tolerance of `docs/math.md` Theorem 6.4 with b = 0.25: during monitoring μ_t never exceeds 1.25 times the smallest μ_i of the confidence-sequence part of the segment's calibration window, and never falls below 0.75 times the largest.

**P1 (Corollary 1 of `docs/theory_per_period.md`, with carry κ = 0).** For every segment and every period (b, b'] of that segment, a period being the monitored part of one UTC day (b' - b ≤ 1440; the first period starts at the end of calibration T): P(the segment raises an alarm in (b, b']) ≤ δ + m / A = 0.01 + 1440 / 36,000 = **5%**, marginally and not given the observed past, provided (H↑) and (H↓) hold for T < t ≤ b', that is **to the end of that day**. Given the observed past (for example no alarm so far) the bound is only P(Gᶜ | past) + m / A, with G the event that the segment's confidence bounds hold. The statement is per segment; after an alarm the next segment has its own G. Because recalibration takes 1440 minutes, the monitored minutes of each UTC day belong to a single segment, so P1 raises at most one alarm per UTC day. **Long-run count (derived here):** the proof of Corollary 6.6 of `docs/math.md` applies unchanged to P1, because a restart only lowers the current statistic M and hence the supermartingale W of that proof; so if (H↑) and (H↓) hold with h = ∞ in every segment, the expected number of alarms in the first M monitored minutes is at most (δ + M / A) / (1 - δ), that is 0.0404 per monitored day plus 0.0101 in total.

**S1 (Theorem 6.4 and Corollary 6.6 of `docs/math.md`).** For each segment, P(an alarm within the first 1440 monitored minutes) ≤ δ + m / A = **5%** if (H↑) and (H↓) hold for those 1440 minutes (part 1); with no assumption, ≤ 5% + P(V_1440), V the event that the tolerance fails within them (part 4). For later days of a segment no per-day bound below 1 is proved: on day d of a segment the same argument gives only δ + d m / A. Long-run count: as for P1 (Corollary 6.6 itself), 0.0404 per monitored day plus 0.0101, if the tolerance holds with h = ∞ in every segment.

**S2 (Theorem B.2 and Remark B.3 of `docs/theory_running_tau.md`).** If the returns of each segment are exchangeable given the returns up to its start (for the first segment: from the first return), P(an alarm within the first 1440 monitored minutes of a segment) ≤ m / A = 1440 / 28,800 = **5%**, exactly, with no baseline, no δ and no truncation level; the band divides every payoff by 1.25 and only lowers the statistic. Long-run count: at most M / A, 0.05 per monitored day (Corollary A.5 with δ = 0). For later days of a segment no per-day bound below 1 is proved (d m / A on day d).

**P0 (plug-in).** Nothing is proved. Nominal: m / A = 5% for the first day of a segment if the baseline were the true one; E14 Q14h measured about twice the nominal rate per day on i.i.d. synthetic data (band 0, A = 28,800).

**What none of these says about real data.** Every hypothesis above excludes volatility clustering beyond the band (TH6: under GARCH(1,1) no path satisfied the tolerance at b = 0.25; E19: calm days resampled in blocks already alarm far above every bound). On real days the null is not known to hold, so the measured rates in section 7 are alarm rates, not false alarm rates. E22 tests whether the release contrast survives a monitor that *carries* a proved per-day bound; it does not test whether the bound holds on these data.

## 6. Families, tests and decision rule

All tests: E12's stratified permutation test, one-sided, outcome `up_post`, all releases (sample "all"; a release at which the monitor is calibrating or blind counts as a release without an up alarm, as in E12's primary). Holm's step-down procedure within each family at 5%, as in E20; no adjustment across families.

| Family | Cells | Tests | Status |
|---|---|---|---|
| **F22a** | P1, BTC, FOMC and CPI | 2 | **primary** (E12's primary family, configuration P1) |
| F22b | P1, ETH and SOL, FOMC and CPI | 4 | secondary (replications, as E20's F12b) |
| F22c | S1, three coins, FOMC and CPI | 6 | secondary |
| F22d | S2, three coins, FOMC and CPI | 6 | secondary (published in E14 Q14g) |
| F22e | P0, as E12 published it | 2 | reference only (reproduction) |

**Primary test: FOMC on BTC under P1.** It passes if its Holm-adjusted p-value in F22a is at most 0.05. CPI on BTC under P1 is in the same Holm family (as in E12) but is a secondary result. Reported effect: releases with an up alarm within 6 hours out of 16, the matched control rate, their difference in percentage points and their ratio, the 95% Wilson interval of the release rate, p and Holm p.

**Exploratory, reported raw, no decision attached:** the same tests restricted to releases at which the monitor was watching (section 7) with controls also watching; the 1-hour outcome (`up_post_1h`); E12's other outcomes (`up_pre`, `down_follow`, `any_post`), which `stratified_tests` computes anyway.

## 7. Diagnostics

**Reproduction check (gate, run first).** P0 is run through the E22 script and its streaming results are compared with E12's published CSVs: every row of `results/e12_tests.csv` with detector "robust e-detector" and design "streaming" (all outcomes, both samples: counts, matched control rates, permutation and Fisher p-values), the matching rows of `e12_event_windows.csv`, and the alarm count and monitored share in `e12_monitor_summary.csv`. The new arms are run only if all agree exactly (floating-point equality after reading the CSVs with round-trip precision); otherwise E22 stops and reports the difference. S2 is also compared with E14's published awake-B rows (`results/e14g_e12_tests.csv`); a difference there is reported, not a stop.

**Alarm rate per day against the bound,** per arm and coin: alarms per UTC day (711 days); alarms per monitored day (monitored minutes / 1440); for P1, the share of monitored periods (segment and UTC day) with an alarm, against 5%; for S1, S2 and P0, the share of segments with an alarm within their first 1440 monitored minutes (segments with at least 1440 minutes left in the stream), against 5%; and alarms per monitored day against the long-run bound (0.0404 for P1 and S1, 0.05 for S2, nominal 0.05 for P0). 95% Wilson intervals, which treat periods and segments as independent and are therefore descriptive only.

**Share of releases watched.** A release is *watched* if at its minute the monitor is monitoring (not calibrating) and, for the X-based arms (P0, P1, S1), the segment's increase side can alarm: E19's blind-segment criterion, the increase side is off when τ² ≤ m_up = U (1 + b) (`side_activity` of `experiments/e19_real_null.py`). For S2 watched means monitoring. Reported per arm, coin and release type: releases watched, the mean share of the 6-hour window watched, the same for the matched control windows, and E19's per-run blind summary (segments with the increase side off, their minutes and share of monitored time, the longest, and the median and 10% quantile of τ² / m_up). If the confidence-sequence arms are blind or recalibrating at most releases, that is reported as the finding.

## 8. Reading fixed in advance

**If the primary test passes** (Holm p ≤ 0.05 for FOMC on BTC under P1): "A monitor configuration with a proved bound of 5% on the probability of a false alarm in each day of a segment (marginally, under a ±25% tolerance on the conditional truncated variance up to the end of the day) raised up alarms within six hours of FOMC statements on BTC more often than at matched control times." E12's FOMC result then does not depend on the plug-in baseline. It still does not show that the bound holds on these data, that the alarms are detections of a lasting change, or anything about CPI.

**If it fails:** the FOMC contrast is not established for a configuration with a per-day guarantee. E12's FOMC result then rests on the plug-in only, and the abstract must not present it as showing what the guarantee buys. The watched shares and the event rate among watched releases will be reported to say descriptively whether blind or recalibrating time or fewer alarms at watched releases is behind the failure; no test is attached to that split.

**Secondary arms.** S1 answers whether the daily restart matters for the release contrast (same baseline, no restart; only the first day of a segment has a bound). S2's result is known (E14 Q14g). Neither can replace or rescue the primary test.

## 9. Outputs

`experiments/e22_guaranteed_events.py`; `results/e22_reproduction.csv`, `e22_monitor_summary.csv`, `e22_windows.csv`, `e22_tests.csv`, `e22_watched.csv`, `e22_alarms.csv`; one figure `results/figures/e22_release_rates.png` (release rate against matched control rate per arm, coin and release type, with the number of releases watched); `tests/test_e22.py` (synthetic data only); `docs/report_parts/E22.md`; `docs/iterations/E22.md`. All results on market data are **[REAL]**; tests are **[SYNTHETIC]**. Runs use `nice -n 10` and at most 2 worker processes. No random numbers enter P0, P1 or S1; S2 uses E14's seed 14000 for its tie-breaking.
