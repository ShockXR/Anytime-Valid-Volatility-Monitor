# W1: does an alarm say more than "the last hour was volatile"? [REAL]

Agent W1, round 4, 4 November 2026. Plan declared at 02:39 UTC in `PLAN.md`, before any W1 result. All results are [REAL] data. Monitor: plug-in (README setting), R4's runs. Labels: VERIFIED, PARTIAL, FAILED, UNVERIFIED.

## 1. Summary in plain words

* **Up alarms carry information beyond the last hour. VERIFIED.** Over all 27 asset-years (2,434 up alarms), 91.4% of up alarms were followed by an hour above the calibration level. Control minutes whose last hour was in the same decile of volatility: 64.0%. Random minutes (R4's placebo): 29.2%. All 27 asset-years pass (p < 0.05, also after Holm). With a tighter match (100 bins instead of 10) the controls reach 74.6%, and 24 of 27 asset-years still pass (Nifty 2022, 2024 and 2025 do not).
* **Down alarms carry almost nothing beyond the last hour. PARTIAL (deciles), FAILED (tighter match).** 95.2% of down alarms were followed by a calmer hour, against 91.9% (deciles) or 92.9% (100 bins) at matched controls, and 70.9% at random minutes. Only 6 and 5 of 27 asset-years pass, none after Holm. So most of R4's down-alarm advantage over random times comes from the last hour already being calm.
* **Post hoc: the up-alarm edge is partly about the last few minutes.** Matching on the last hour, 6 hours and day together leaves the edge at +16.0 points (27 of 27 pass). Adding the last 10 minutes halves it to +7.8 points (20 of 27 pass).
* **On fresh 2026 data the alarm does not improve a strong forecast of the next hour. FAILED for the next hour.** Against E20's HAR-I forecast (daily, weekly, monthly and last-hour realised variance), the alarm dummies are not significant for BTC, ETH or Nifty, and the out-of-sample loss (QLIKE) changes by -0.54% to +0.05%, with every interval covering 0.
* **For the next day, up alarms point the wrong way on crypto. PARTIAL by the declared rule.** After an up alarm, next-day realised variance on the five coins was about 10% *below* the HAR-I forecast (coefficient -0.104, interval -0.160 to -0.042). HAR-I overreacts to the burst that set off the alarm. Only one of six primary out-of-sample cells passed (ETH next day, QLIKE -0.10%), a negligible gain.
* **Post hoc: against a HAR forecast without the last hour, alarms do add information** (up coefficients +0.48 to +0.72 for the next hour, all intervals above 0). The alarm is a good summary of the recent burst, but a forecast that already sees the last hour has no use for it.

## 2. What was run

* **Data.** R4's run arrays (`reports/round4/R4/run_monitors.py`): Binance BTC, ETH, XRP, SOL, DOGE (21 Jan 2024 to 24 Oct 2026), Bitstamp BTC/USD (28 Jan 2025 to 24 Oct 2026), Nifty 50 (Apr 2017 to 23 Oct 2026). Part 2 also reads raw 1-minute log returns: `load_binance` (1 Jan 2024 to 24 Oct 2026) and `e07_nifty.load_clean` (the loader reads to 1 Nov; everything from 25 Oct on is dropped before any computation, as in R4). No FOMC or CPI time was used.
* **Part 1** (`p1_matched_placebo.py`, `p1_summarise.py`). R4's definitions: pre = realised variance (RV) of the 60 minutes ending at the minute, post = RV of the next 60, base = RV of the minute's calibration window, at least 50 finite returns. The match variable is pre / base. The pool is the watched minutes of the same asset and year, alarm minutes removed. Strata are deciles (primary) or 100 bins (declared sensitivity) of the pool. For Nifty, strata also use the session hour (7 buckets from 09:15 IST). This was feasible: no stratum was empty. For each alarm, one control is drawn from its stratum; this is repeated 2,000 times. One-sided p = (1 + sets at least as good) / 2,001, so the floor is 0.0005. The alarm hit shares equal R4's (for example BTC 2025 up 92.8%, Nifty 2024 up 97.7%).
* **Part 1 post hoc** (`p1_posthoc_knn.py`). The control is drawn from the 200 nearest pool minutes in log(pre / base) over 60 minutes, then also over 6 hours and 1 day (Nifty: 180 and 375 minutes), then also over 10 minutes.
* **Part 2** (`p2_forecast_2026.py`). E20's `ols_dk`, `qlike` and `cbb_weights` are imported. E20's hourly grid, HAR-I design and daily refitted expanding-window fit are copied, with horizons 60 minutes and one day. Nifty grid: ends of session hours, and a day is one session of 375 minutes. Nifty windows allow up to 10% missing returns (frozen stretches), rescaled; crypto follows E20's rule (no missing returns). U (D) = a plug-in up (down) alarm in the hour ending at the grid time. Evaluation: grid times from 1 Jan 2026 whose outcome window ends by 24 Oct 2026 23:59. Metric A: OLS of log RV on log F(HAR-I) + U + D, asset intercepts, 95% interval from a circular day-block bootstrap (7 days, 1,999 draws, days shared across coins). Metric B: HAR-I + U + D fitted the same past-only way, mean QLIKE difference against HAR-I.
* **Part 2 post hoc** (`p2_posthoc_har.py`): the same regression against plain HAR (no last-hour term).
* **Seeds** (numpy default_rng): Part 1 9,980,000 + 1,000 x asset + 10 x (year - 2017) + direction + 2 for 100 bins (in `p1_matched_placebo.csv`). Post hoc matching 9,987,000 to 9,988,819. Part 2 9,989,000 to 9,989,071, post hoc 9,989,100 to 9,989,121. Run: `python p1_matched_placebo.py; python p1_summarise.py; python p1_posthoc_knn.py; python p2_forecast_2026.py; python p2_posthoc_har.py` (each under 1 minute, single core). First rebuild R4's run arrays with `python reports/round4/R4/run_monitors.py NIFTY BTC ETH XRP SOL DOGE BITSTAMP` (about 2 minutes; they are not committed). The coordinator reran `p1_matched_placebo.py`, `p1_summarise.py` and `p2_forecast_2026.py` on a separate copy: the output files were identical byte for byte (Part 1) and to every digit (Part 2).

## 3. Results

### Part 1, pooled over asset-years [REAL]

How to read this: "hit" = share of alarms whose next hour was above (up) or below (down) the calibration level. "Cells passing" = asset-years with p < 0.05 out of 27. The RV ratio is next-hour RV after the alarm over that after its control (geometric mean); above 1 means more volatile.

| Direction, match | Alarms | Hit, alarms | Hit, matched controls | Hit, random (R4) | Difference, points (range over cells) | Cells passing (Holm) | RV ratio, median (range) | Verdict |
|---|---|---|---|---|---|---|---|---|
| Up, deciles | 2,434 | 91.4% | 64.0% | 29.2% | +27.4 (+14.5 to +39.0) | 27 (27) | 1.87 (1.52 to 2.32) | VERIFIED |
| Up, 100 bins | 2,434 | 91.4% | 74.6% | 29.2% | +16.8 (+0.9 to +26.9) | 24 (21) | 1.50 (1.13 to 1.85) | VERIFIED |
| Down, deciles | 2,012 | 95.2% | 91.9% | 70.9% | +3.3 (-11.1 to +12.0) | 6 (0) | 0.84 (0.73 to 1.03) | PARTIAL |
| Down, 100 bins | 2,012 | 95.2% | 92.9% | 70.9% | +2.3 (-10.4 to +12.8) | 5 (0) | 0.95 (0.82 to 1.13) | FAILED |

* The decile match is loose for up alarms: the alarms' median pre / base is 2.1 to 3.2, their decile controls' 1.6 to 2.1. With 100 bins the control median is within 10% of the alarm median in every cell (0.91 to 1.03 times; down alarms 0.95 to 1.05 times). So the 100-bin row is the fairer one.
* For down alarms the next-hour RV ratio is below 1 in 18 of 27 cells with deciles (p < 0.05), but only 6 of 27 with 100 bins.
* V1's quick probe (BTC 2025: 79%; Nifty 2024: 98% at similar recent volatility) agrees in direction. Here, with 100 bins: 75% and 94%.

### Part 1, per asset-year [REAL]

How to read this: n = alarms; then hit share of alarms, of decile controls (p), of 100-bin controls (p); last: RV ratio alarm over control, deciles / 100 bins. Full table: `p1_matched_placebo.csv`.

| Asset-year | Up n | Up hit | Decile control (p) | 100-bin control (p) | RV ratio | Down n | Down hit | Decile control (p) | 100-bin control (p) | RV ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| NIFTY 2017 | 28 | 96% | 59% (<0.001) | 73% (<0.001) | 1.93 / 1.59 | 17 | 94% | 84% (0.212) | 85% (0.223) | 0.80 / 0.84 |
| NIFTY 2018 | 30 | 97% | 63% (<0.001) | 80% (0.002) | 1.76 / 1.27 | 26 | 100% | 88% (0.037) | 87% (0.027) | 0.79 / 0.82 |
| NIFTY 2019 | 32 | 94% | 71% (0.001) | 77% (0.009) | 1.67 / 1.48 | 25 | 100% | 93% (0.170) | 95% (0.272) | 1.03 / 1.13 |
| NIFTY 2020 | 38 | 97% | 73% (<0.001) | 75% (<0.001) | 1.77 / 1.55 | 30 | 100% | 92% (0.089) | 93% (0.109) | 0.74 / 0.86 |
| NIFTY 2021 | 36 | 92% | 72% (0.005) | 75% (0.005) | 1.52 / 1.22 | 37 | 97% | 95% (0.433) | 96% (0.515) | 0.94 / 0.97 |
| NIFTY 2022 | 31 | 94% | 71% (<0.001) | 84% (0.080) | 1.70 / 1.35 | 27 | 93% | 93% (0.704) | 96% (0.900) | 0.82 / 0.94 |
| NIFTY 2023 | 25 | 96% | 68% (0.001) | 81% (0.025) | 1.75 / 1.50 | 27 | 89% | 91% (0.778) | 92% (0.825) | 0.91 / 0.97 |
| NIFTY 2024 | 43 | 98% | 83% (0.002) | 94% (0.214) | 1.67 / 1.13 | 37 | 97% | 93% (0.272) | 94% (0.352) | 0.95 / 1.04 |
| NIFTY 2025 | 36 | 86% | 71% (0.016) | 85% (0.553) | 1.63 / 1.31 | 33 | 94% | 91% (0.396) | 93% (0.564) | 0.80 / 1.01 |
| NIFTY 2026 | 25 | 96% | 77% (0.014) | 78% (0.006) | 1.55 / 1.58 | 22 | 82% | 93% (0.980) | 92% (0.989) | 1.03 / 1.13 |
| BTC 2024 | 133 | 94% | 63% (<0.001) | 70% (<0.001) | 1.98 / 1.67 | 112 | 93% | 91% (0.316) | 92% (0.403) | 0.85 / 0.96 |
| BTC 2025 | 138 | 93% | 67% (<0.001) | 75% (<0.001) | 1.88 / 1.52 | 116 | 94% | 90% (0.101) | 91% (0.220) | 0.80 / 0.92 |
| BTC 2026 | 117 | 91% | 60% (<0.001) | 71% (<0.001) | 1.93 / 1.47 | 78 | 97% | 91% (0.019) | 91% (0.024) | 0.75 / 0.86 |
| ETH 2024 | 143 | 87% | 61% (<0.001) | 70% (<0.001) | 1.77 / 1.50 | 104 | 91% | 90% (0.381) | 91% (0.580) | 0.88 / 0.97 |
| ETH 2025 | 142 | 89% | 63% (<0.001) | 74% (<0.001) | 1.94 / 1.61 | 103 | 92% | 92% (0.486) | 94% (0.794) | 0.79 / 0.91 |
| ETH 2026 | 116 | 81% | 54% (<0.001) | 62% (<0.001) | 1.86 / 1.65 | 78 | 91% | 90% (0.479) | 92% (0.666) | 0.88 / 1.05 |
| XRP 2024 | 127 | 93% | 54% (<0.001) | 66% (<0.001) | 2.32 / 1.85 | 120 | 97% | 91% (0.013) | 92% (0.027) | 0.91 / 0.94 |
| XRP 2025 | 131 | 91% | 70% (<0.001) | 82% (0.005) | 1.70 / 1.36 | 126 | 96% | 96% (0.519) | 96% (0.661) | 0.89 / 0.95 |
| XRP 2026 | 104 | 96% | 62% (<0.001) | 74% (<0.001) | 2.00 / 1.57 | 87 | 94% | 91% (0.245) | 93% (0.369) | 0.82 / 0.95 |
| SOL 2024 | 128 | 95% | 65% (<0.001) | 81% (<0.001) | 1.76 / 1.37 | 101 | 97% | 93% (0.058) | 94% (0.123) | 0.89 / 0.97 |
| SOL 2025 | 136 | 93% | 68% (<0.001) | 76% (<0.001) | 1.87 / 1.46 | 108 | 98% | 93% (0.012) | 94% (0.026) | 0.79 / 0.85 |
| SOL 2026 | 105 | 89% | 66% (<0.001) | 76% (<0.001) | 1.87 / 1.54 | 84 | 99% | 94% (0.043) | 96% (0.111) | 0.84 / 0.94 |
| DOGE 2024 | 112 | 94% | 61% (<0.001) | 78% (<0.001) | 1.87 / 1.41 | 127 | 99% | 95% (0.006) | 95% (0.013) | 1.00 / 1.05 |
| DOGE 2025 | 134 | 95% | 67% (<0.001) | 75% (<0.001) | 2.02 / 1.62 | 112 | 96% | 93% (0.197) | 94% (0.282) | 0.84 / 0.99 |
| DOGE 2026 | 110 | 88% | 62% (<0.001) | 70% (<0.001) | 1.88 / 1.59 | 86 | 94% | 93% (0.403) | 94% (0.569) | 0.91 / 1.03 |
| BITSTAMP 2025 | 122 | 88% | 65% (<0.001) | 80% (0.018) | 1.87 / 1.33 | 107 | 94% | 90% (0.079) | 91% (0.140) | 0.77 / 0.87 |
| BITSTAMP 2026 | 112 | 88% | 59% (<0.001) | 71% (<0.001) | 2.00 / 1.54 | 82 | 94% | 90% (0.141) | 91% (0.200) | 0.73 / 0.84 |

### Part 1 post hoc, nearest-neighbour matching [REAL, post hoc]

How to read this: pooled over 27 asset-years; "+" adds a look-back to the match. File `p1_posthoc_knn_summary.csv`.

| Matched on | Up: alarms vs controls | Up cells passing | Down: alarms vs controls | Down cells passing |
|---|---|---|---|---|
| last hour | 91.4% vs 75.0% | 25 of 27 | 95.2% vs 92.9% | 5 of 27 |
| + last 6 hours + last day | 91.4% vs 75.4% | 27 of 27 | 95.2% vs 93.9% | 2 of 27 |
| all of the above + last 10 minutes | 91.4% vs 83.5% | 20 of 27 | 95.2% vs 93.6% | 1 of 27 |

### Part 2, 2026 only, against HAR-I [REAL]

How to read this: b = slope on the log forecast (1 means the forecast is right on average). cU, cD = shift in log next-period RV in the hour after an up or down alarm, beyond the forecast (-0.10 means about 10% lower), with 95% day-block intervals. Last column: change in mean QLIKE loss when the alarm dummies are added to the forecast (negative = better), with interval. Pass lines: cU interval above 0, cD interval below 0, QLIKE interval below 0. SOL, XRP, DOGE and BTC+ETH are in `p2_forecast_2026.csv`.

| Asset | Outcome | Grid points | Up / down alarm points | b | cU (interval) | cD (interval) | QLIKE, HAR-I | QLIKE change (interval) |
|---|---|---|---|---|---|---|---|---|
| BTC | next hour | 7127 | 117 / 79 | 1.00 | -0.024 (-0.160 to +0.117) | +0.010 (-0.109 to +0.125) | 0.2849 | +0.05% (-0.01% to +0.11%) |
| BTC | next day | 7104 | 116 / 79 | 1.03 | -0.121 (-0.190 to -0.048) | -0.062 (-0.185 to +0.057) | 0.1505 | +0.03% (-0.07% to +0.15%) |
| ETH | next hour | 7127 | 116 / 78 | 1.02 | +0.005 (-0.128 to +0.162) | +0.054 (-0.077 to +0.195) | 0.3245 | +0.02% (-0.03% to +0.09%) |
| ETH | next day | 7104 | 116 / 78 | 1.06 | -0.088 (-0.165 to -0.008) | -0.116 (-0.220 to -0.027) | 0.1738 | -0.10% (-0.18% to -0.03%) |
| NIFTY | next hour | 1286 | 25 / 23 | 0.98 | +0.205 (-0.031 to +0.506) | +0.160 (-0.079 to +0.417) | 0.1279 | -0.54% (-1.87% to +0.36%) |
| NIFTY | next day | 1337 | 26 / 25 | 1.02 | +0.140 (+0.000 to +0.298) | +0.105 (-0.080 to +0.284) | 0.1001 | -0.09% (-0.26% to +0.05%) |
| CRYPTO_POOLED | next hour | 35635 | 552 / 414 | 1.00 | -0.029 (-0.110 to +0.058) | -0.020 (-0.072 to +0.032) | 0.2769 | +0.03% (-0.02% to +0.10%) |
| CRYPTO_POOLED | next day | 35520 | 550 / 414 | 1.01 | -0.104 (-0.160 to -0.042) | -0.050 (-0.117 to +0.020) | 0.1702 | -0.03% (-0.06% to +0.00%) |

* **Next hour:** no pass in any primary cell. Pooled crypto: cU -0.029, cD -0.020, both intervals cover 0.
* **Next day:** crypto up coefficients are negative (BTC -0.121, ETH -0.088, pooled five coins -0.104); the declared line asks for positive, so they fail. Down: ETH -0.116 passes, BTC and Nifty do not. Nifty up +0.140 has its lower end at 0.000.
* **Declared verdict:** Metric B passes in 0 of 3 primary assets for the next hour and 1 of 3 (ETH) for the next day, so **PARTIAL** by the letter of the plan. The one gain is 0.10% of the loss.
* **Post hoc, plain HAR without the last hour** (`p2_posthoc_har.csv`): cU = +0.71 (BTC), +0.72 (ETH), +0.48 (Nifty) for the next hour, +0.14 to +0.30 for the next day; every interval is above 0. Crypto down alarms: -0.20 to -0.36, intervals below 0.

## 4. README or limitations claims that need changing

| Where | Claim now | What W1 found | Suggested change |
|---|---|---|---|
| README.md line 37 | Alarms "beat random times but not a strong forecast"; against HAR-I "decrease alarms a little (about 6% lower)" | On Jan to Oct 2026, against HAR-I, neither direction helps the next hour; the down effect for the next day holds on ETH only (BTC -0.062, interval -0.185 to +0.057); up alarms point 10% the wrong way | Add: "On fresh 2026 data (W1) the alarms did not improve HAR-I for the next hour, and after up alarms next-day variance was about 10% below HAR-I" |
| docs/limitations.md line 19 | Random-time and forecast controls "reduce but do not remove" the clustering effect | At matched last-hour volatility, the down-alarm edge falls from 24 points to 2 to 3 points (mostly not significant); the up-alarm edge falls from 62 to 17 points but stays | Add the numbers; say that decrease alarms add little beyond "the last hour was calm" |
| R4 REPORT, "useful years" (27 of 27 up, 26 of 27 down) | Alarms beat a random-time placebo | True, but at matched recent volatility only 5 to 6 of 27 down cells pass | Note the W1 matched result beside R4's placebo line |

## 5. What this does not show

* **No ground truth.** "Hit" compares with the calibration level, not with a true change.
* **Dependence.** Part 1 p-values treat alarms in one asset-year as independent. Alarms are at least one calibration apart, but days are not independent. Part 2 uses day blocks.
* **The match is on the minute, not the process.** Controls can be minutes just before an alarm. The post hoc result (the edge halves when the last 10 minutes are matched) suggests much of the up-alarm edge is "the burst is happening right now", not a slow regime change.
* **Part 1 and Part 2 measure different things.** Part 1 checks the exact alarm minute against the calibration level. Part 2 checks clock-hour grid points (0 to 59 minutes after the alarm) against a forecast that already uses the last hour.
* **Small 2026 samples.** Nifty has 25 up and 23 down alarm points; BTC and ETH about 117 and 78. Only the plug-in monitor at band 0.25 was used, not E6's e-detector runs of E20, so the numbers are not directly comparable with E20's.
* **Post hoc items** (not in PLAN.md): the nearest-neighbour matching and the plain-HAR regression.

**Files in `reports/round4/W1/`:** `PLAN.md`, `REPORT.md`; scripts `p1_matched_placebo.py`, `p1_summarise.py`, `p1_posthoc_knn.py`, `p2_forecast_2026.py`, `p2_posthoc_har.py` with their `.log` files; results `p1_matched_placebo.csv`, `p1_verdicts.csv`, `p1_by_asset.csv`, `p1_posthoc_knn.csv`, `p1_posthoc_knn_summary.csv`, `p2_forecast_2026.csv`, `p2_posthoc_har.csv`.
