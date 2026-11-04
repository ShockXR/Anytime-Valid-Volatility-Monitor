# R4 plan: walk-forward real-data backtest by asset and year (declared before any run)

Written 2026-11-04 02:01 UTC, before any monitor was run by R4. Nothing below is changed after seeing results; any later change is listed in REPORT.md as post hoc.

## Data and preparation (reused from the experiment scripts)

* Binance BTC, ETH, XRP, SOL, DOGE (USDT), loaded with `avmonitor.data.load_binance(download=False)` from 2024-01-01. Returns divided by `intraday_scale(lookback_days=20)` (start="past", circular), exactly as `e12_crypto_events.load_coin`. The first 20 days only warm up the profile; the monitored stream runs from 2024-01-21 00:00 UTC.
* Bitstamp BTC/USD: `load_bitstamp_csv` with the loader default `filled="missing"` (no zero fill). The file starts on 2025-01-07, so the stream runs from 20 days after the first bar.
* Nifty 50: `e07_nifty.load_clean()` (regular sessions, frozen runs set to missing) and `e07_nifty.deseasonalise` (session-aware profile, 20 past sessions). Stream from the first 2017 session.
* End of every stream: 2026-10-25 00:00 UTC (Nifty: last session before 25 October 2026). Reason: the pre-registered forward tests cover FOMC and CPI releases from 25 October 2026; stopping before that date keeps R4 away from them. This drops one week of October 2026.
* Missing returns are passed as NaN (skipped by the monitor), never filled with zeros.
* One continuous run per asset and monitor (no restarts at year boundaries). Results are then split by calendar year of the alarm or minute (UTC for crypto, IST session date for Nifty).

## Monitors

* README setting, plug-in: `VolatilityMonitor(threshold=A, c=3, cal_len=L, null_band=0.25)`; crypto A = 28,800, L = 1,440; Nifty A = 7,500, L = 750 (E7 primary).
* cs baseline: `VolatilityMonitor(threshold=A_cs, c=3, cal_len=L, null_band=0.25, baseline="cs", delta=0.01)`; crypto A_cs = 36,000; Nifty A_cs = 375 / 0.04 = 9,375 (E14, E19).

## Metrics (per asset, per calendar year, per monitor)

1. **Watched minute**: a step at which the monitor is monitoring (calibrated) and the return is finite. **Alarms per 1,000 watched minutes**, up and down separately.
2. **Blind share**: share of stream steps (minutes for crypto, session minutes for Nifty) spent calibrating. **Blind stretch**: a maximal run of calibrating steps; median and longest length (in steps, also in days of 1,440 or sessions of 375). A stretch is assigned to the year in which it starts. For cs, also the **increase-side-off share**: monitored steps in segments with tau² <= m_up (E19's blind mode), reported separately and added to the calibrating share as "blind for increases".
3. **Post / pre ratio**: for an alarm at step t, RV_pre = mean of u² over steps t-59..t and RV_post = mean of u² over t+1..t+60 (u = deseasonalised return; at least 50 of 60 finite, else the alarm is dropped from this metric). Statistic: median of log(RV_post / RV_pre) over the alarms of the asset-year, reported as exp(median).
4. **Placebo**: for each alarm, a random watched minute from the same asset and year (for Nifty also the same minute of the session), same windows. 2,000 placebo sets, each with one draw per alarm. One-sided Monte Carlo p-value: up alarms, P(placebo median >= alarm median); down alarms, P(placebo median <= alarm median). p = (1 + count) / (1 + 2,000).
5. **Hit share vs baseline**: share of up alarms whose RV_post exceeds the calibration baseline RV (mean of u² over the finite returns of the calibration window of the segment in which the alarm fired); for down alarms, share with RV_post below it. Placebo: the same share for the placebo draws against the baseline in force at the placebo minute; same 2,000 sets and one-sided p-value.
6. **Restart split (question 4)**: monitored steps from the end of calibration to the alarm: early <= 60, middle 61 to one day (1,440 crypto, 375 Nifty), late > one day. Hit share and post / pre ratio per bin, pooled over crypto assets and years, and Nifty separately. Permutation test (2,000 label shuffles) of early against late on the hit share.

Seeds (numpy default_rng): 9,800,001 (crypto placebo sets, offset by asset index 0 to 5 and year), 9,800,101 (Nifty placebo), 9,800,201 (restart permutation tests). Exact seeds per cell are listed in the output CSV.

## Pass / fail lines

* **Q1 stability.** For each asset and monitor, a year's alarm rate (up plus down per 1,000 watched minutes) is stable if it is within ±30% of that asset's median yearly rate. Usefulness holds in an asset-year if the up hit share beats its placebo with p < 0.05 and the down hit share beats its placebo with p < 0.05. Q1 PASS if usefulness holds in at least 80% of asset-years for each direction and every asset's yearly rate is stable; PARTIAL if usefulness holds in 50% to 80%; FAILED below 50%.
* **Q2 2026 like 2024 to 2025.** For each asset with 2026 data: 2026 rate within ±30% of the mean of its earlier full years, blind share within ±10 points, up and down hit shares within ±10 points. PASS if all four hold for at least 5 of the 7 assets (5 coins, Bitstamp BTC, Nifty; Bitstamp compared with 2025 only).
* **Q3 "blind about 70%".** The README claim holds in an asset-year if the plug-in blind share is between 60% and 80%. VERIFIED for crypto if at least 80% of crypto asset-years are in range; checked separately for Nifty (README line 166 says "on real data", not only crypto). cs is described, not tested.
* **Q4 early alarms.** Early alarms are "less useful" if their up hit share is at least 10 points below the late alarms' and the permutation p < 0.05.
* Order of work if time runs short: Nifty all years, BTC, ETH, then XRP, SOL, DOGE, Bitstamp.
