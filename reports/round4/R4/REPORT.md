# R4: walk-forward backtest by asset and year, 2017 to October 2026 [REAL]

Agent R4, round 4, 4 November 2026. The plan was declared at 02:01 UTC in `PLAN.md`, before any run. Everything here is [REAL] unless marked [SYNTHETIC]. Status labels: VERIFIED, PARTIAL, FAILED, UNVERIFIED.

## 1. Summary in plain words

* **"Blind about 70% of the time" holds for crypto but not for the Nifty 50. PARTIAL.** We used the README setting: plug-in baseline, band 0.25. The monitor spent 64% to 72% of the minutes recalibrating in each of the 17 crypto asset-years. On the Nifty 50 it was 43% to 65% (52% over 2017 to 2026). Only 2024 reached 60%. README line 166 says "on real data". It should say "on crypto".
* **Alarms land on the right side of the calibration level in every year and on every asset. VERIFIED.** After an up alarm, the next 60 minutes had more variance than the calibration window in 81% to 98% of cases per asset-year. At random watched minutes this was 25% to 38%. After a down alarm, the next 60 minutes had less variance in 82% to 100% of cases (random: 60% to 76%). Alarms beat the random-time control with p of about 0.001 in 27 of 27 asset-years for up alarms and 26 of 27 for down alarms.
* **But the hour after an up alarm is calmer than the hour before it. FAILED on the declared post / pre test.** The post / pre ratio is the realised variance of the 60 minutes after an alarm divided by that of the 60 minutes before. For plug-in up alarms its median was 0.59 to 0.96, against about 1.0 at random times, in all 27 asset-years. Down alarms show the mirror image (1.07 to 1.58). A post hoc simulation [SYNTHETIC] shows most of this is built in. Even when volatility rises for good, the ratio is 0.79 to 0.91, because an alarm fires right after a run of large returns. So an up alarm means "volatility is now well above the calibrated level". It does not mean "volatility will keep rising".
* **Alarm behaviour is stable across years on crypto, less so on the Nifty 50. PARTIAL.** On every coin, plug-in alarms per 1,000 watched minutes stayed within ±30% of that coin's median in every year (1.22 to 1.75). On the Nifty 50 they ranged from 0.99 (2023) to 2.48 (2024). 2021 and 2024 fall outside the ±30% line.
* **January to October 2026 looks like 2024 and 2025 for the plug-in. VERIFIED (6 of 7 assets).** In 2026, crypto alarm rates were 9% to 26% lower, mostly because down alarms fell (11% to 31% fewer). Blind time was 2 to 6 points lower. All of this is within the declared limits. The one miss: Nifty down alarms in 2026 were confirmed 82% of the time against 96% before (22 alarms).
* **The cs baseline does not pass the 2026 check (4 of 7 assets). FAILED.** On BTC its 2026 alarm rate fell 42%. From 16 August to 4 October 2026 it could not raise an increase alarm for 49 days. **New: Bitstamp BTC/USD shows the same blind mode at the same time (29 days, 16 August to 14 September 2026).** So it is not only a Binance tick effect.
* **cs is blind less often than the plug-in, but not much less for increases.** cs spent 41% to 57% of crypto minutes calibrating, and 14% to 44% on the Nifty 50. Some cs segments cannot raise an increase alarm at all. Counting those, BTC was blind for increases 54% to 58% of the time.
* **Alarms right after a restart are not worse.** An early alarm is an up alarm within 60 minutes of the end of a calibration (4% to 8% of up alarms). Early alarms were confirmed about as often as alarms more than a day later: Binance plug-in 96% against 93%, Nifty 91% against 95%, p = 0.33 to 0.90. No plug-in down alarm ever fired within 60 minutes of a restart.

## 2. What was run

**Data.** All real data. One continuous run per asset and monitor; results split by calendar year.

| Asset | Monitored stream | Steps | Missing returns |
|---|---|---|---|
| BTC, ETH, XRP, SOL, DOGE (Binance USDT, 1 minute) | 21 Jan 2024 to 24 Oct 2026 (1 to 20 Jan 2024 only warm up the profile) | 1,451,520 each | 0 |
| Bitstamp BTC/USD | 28 Jan 2025 to 24 Oct 2026 (file starts 7 Jan 2025) | 914,400 | 10,866 (filled minutes made missing, loader default) |
| Nifty 50 | 3 Apr 2017 to 23 Oct 2026, 2,356 regular sessions | 883,486 | 3,287 |

* **Data preparation reused, not rewritten.** Crypto: as `e12_crypto_events.load_coin` (`intraday_scale`, 20 past days, start="past"). Bitstamp: `load_bitstamp_csv(filled="missing")`. Nifty: `e07_nifty.load_clean` and `e07_nifty.deseasonalise` (session profile). Missing returns are passed as missing, never as zeros.
* **Cut-off.** All streams stop before 25 October 2026, when the pre-registered forward tests start. No FOMC or CPI time was used anywhere.
* **Fresh 2026 files.** R4 fetched XRP and DOGE for January to October 2026 with `avmonitor.data.fetch` (20 files, hashes in `fetched_2026_files.csv`). These files are not in `data/MANIFEST.csv`, and neither are the BTC, ETH and SOL 2026 files.
* **Monitors.** Plug-in (README setting): `VolatilityMonitor(threshold=A, c=3, cal_len=L, null_band=0.25)`. On crypto A = 28,800 and L = 1,440. On Nifty A = 7,500 and L = 750 (E7 primary). cs: the same with `baseline="cs", delta=0.01`, and A = 36,000 (crypto) or 9,375 (Nifty).
* **Reproduction check.** The plug-in alarms up to 31 December 2025 are identical, minute and direction, to E12's committed lists for BTC (499), ETH (492) and SOL (473). The 610 Nifty alarms are identical to E7's committed primary list (E7 has one more, on 27 October 2026, after R4's cut-off).

**Commands.** Run from the repository root with the venv active and `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`:

```
python reports/round4/R4/run_monitors.py NIFTY BTC ETH XRP SOL DOGE BITSTAMP   # about 2 minutes, writes runs/*.npz
cd reports/round4/R4
python analyse.py             # asset_year.csv, alarms.csv, blind_stretches.csv, asset_summary.csv
python summarise.py           # verdicts.txt, q1_stability.csv, q2_2026.csv, q3_blind.csv, q4_restart_*.csv
python extra_stats.py         # blind_for_up_stretches.csv (post hoc)
python posthoc_selection.py   # posthoc_selection.csv [SYNTHETIC, post hoc]
python figures.py             # fig_blind_share.png, fig_post_pre_ratio.png
```

The run arrays (`runs/*.npz`, 186 MB) are not committed; `run_monitors.py` rebuilds them in about 2 minutes.

**Seeds** (numpy `default_rng`):
* Placebo sets: 9,800,000 + 1,000 × asset + 100 × monitor + 10 × (year − 2017) + direction. That spans 9,800,000 to 9,806,191. Each cell's seed is in `asset_year.csv`, columns `up_seed` and `down_seed`.
* Restart permutation tests: 9,800,201.
* Post hoc simulation: 9,800,301 to 9,800,304.

This seed formula replaces the shorter list in PLAN.md. It was fixed while writing the code, before any result was seen.

**Metrics, defined once.**
* *Watched minute*: a minute in which the monitor is watching (not calibrating) and the return exists.
* *Blind share*: the share of minutes spent calibrating.
* *Blind stretch*: an unbroken run of calibrating minutes.
* *Realised variance (RV)* of a window: the average squared deseasonalised return.
* *Post / pre ratio*: RV of the 60 minutes after the alarm divided by RV of the 60 minutes up to it.
* *Hit share*: for up alarms, the share whose next 60 minutes have RV above the RV of the calibration window the alarm was measured against. For down alarms, below.
* *Placebo*: for every alarm, one random watched minute from the same asset and year (for Nifty, also the same minute of the session). We drew 2,000 such sets. The p-value is the share of sets that do at least as well as the alarms. The smallest possible value is 0.0005.

## 3. Results

### Q1. Is alarm behaviour stable across years and assets? PARTIAL

How to read this: one row per asset, plug-in monitor. "Rate" is up plus down alarms per 1,000 watched minutes, lowest and highest year. Hit shares are the range over years. The placebo column is the highest random-time value in any year. "Useful years" counts years where the hit share beat the placebo with p < 0.05.

| Asset (years) | Rate, low to high | Years within ±30% of median | Up hit share | Placebo (max) | Useful years, up | Down hit share | Placebo (max) | Useful years, down | Post / pre, up | Post / pre, down |
|---|---|---|---|---|---|---|---|---|---|---|
| Nifty 50 (10) | 0.99 to 2.48 | 8 of 10 | 86% to 98% | 38% | 10 of 10 | 82% to 100% | 70% | 9 of 10 | 0.73 to 0.96 | 1.07 to 1.58 |
| BTC (3) | 1.34 to 1.70 | 3 of 3 | 91% to 94% | 30% | 3 of 3 | 93% to 97% | 73% | 3 of 3 | 0.66 to 0.78 | 1.26 to 1.36 |
| ETH (3) | 1.31 to 1.75 | 3 of 3 | 81% to 89% | 29% | 3 of 3 | 91% to 92% | 71% | 3 of 3 | 0.59 to 0.71 | 1.15 to 1.37 |
| XRP (3) | 1.26 to 1.75 | 3 of 3 | 91% to 96% | 28% | 3 of 3 | 94% to 97% | 74% | 3 of 3 | 0.73 to 0.84 | 1.19 to 1.27 |
| SOL (3) | 1.22 to 1.40 | 3 of 3 | 89% to 95% | 31% | 3 of 3 | 97% to 99% | 71% | 3 of 3 | 0.67 to 0.78 | 1.08 to 1.32 |
| DOGE (3) | 1.36 to 1.55 | 3 of 3 | 88% to 95% | 30% | 3 of 3 | 94% to 99% | 76% | 3 of 3 | 0.69 to 0.87 | 1.25 to 1.41 |
| Bitstamp BTC (2) | 1.33 to 1.63 | 2 of 2 | 88% | 29% | 2 of 2 | 94% | 72% | 2 of 2 | 0.69 to 0.74 | 1.28 to 1.39 |

* **Declared line.** PASS needs usefulness in at least 80% of asset-years in each direction, and every asset's yearly rate stable.
  * Plug-in: usefulness in 27 of 27 asset-years (up) and 26 of 27 (down). But Nifty 2021 (1.95) and 2024 (2.48) break the rate line, so **PARTIAL**.
  * cs: usefulness in 27 of 27 (up) and 22 of 27 (down). The misses are Nifty years with only 5 to 9 down alarms. Rates are unstable for Nifty in 4 years and for BTC in 2026. **PARTIAL**.
* **Post / pre.** The alarms beat the placebo in **0 of 27** asset-years, in either direction, for both monitors. For cs the up range is 0.62 to 1.12 (`q1_stability.csv`). Read the post hoc check below before taking this as "alarms are useless".
* **Alarms per calendar day (plug-in)** stay in a narrow band: 0.64 to 0.71 on crypto in every year, 0.21 to 0.33 on Nifty. Blind time caps this number, which is why it hardly moves.

**Post hoc check of the post / pre metric [SYNTHETIC, not declared].** Independent returns, plug-in monitor as above. Volatility changes for good 300 steps after calibration. 300 paths per row (`posthoc_selection.csv`).

| Returns | Change in sd | Median post / pre of the alarms | Share with post below pre |
|---|---|---|---|
| Gaussian | ×1.5 (up) | 0.91 | 65% |
| Student t, 3 degrees of freedom | ×1.5 (up) | 0.79 | 67% |
| Gaussian | ×0.6 (down) | 1.05 | 44% |
| Student t, 3 degrees of freedom | ×0.6 (down) | 1.31 | 29% |

How to read this: even when the change lasts and returns do not cluster, the hour after an up alarm looks calmer than the hour before. The reason is that the large returns in that hour are what set off the alarm.
* On real data, 61% to 68% of up alarms had post below pre.
* The real ratios (0.59 to 0.96) are a little lower than the Student t row, so some real calming after the alarm is likely. The declared metric cannot separate that from the built-in effect.
* The hour after an up alarm still had a median of 2.1 to 2.6 times the calibration RV. After a down alarm it was 0.25 to 0.40 times.

### Q2. Does January to October 2026 look like 2024 and 2025? Plug-in VERIFIED, cs FAILED

How to read this: plug-in monitor. Each cell reads "earlier years → 2026". Earlier years are 2024 and 2025 for the coins, 2025 for Bitstamp, and 2018 to 2025 for Nifty. Limits declared in advance: rate within ±30%, blind share within ±10 points, each hit share within ±10 points.

| Asset | Rate per 1,000 watched min | Blind share | Up hit share | Down hit share | Plug-in within all limits | cs within all limits |
|---|---|---|---|---|---|---|
| Nifty 50 | 1.53 → 1.41 (−8%) | 52% → 51% | 94% → 96% | 96% → **82%** | no | no (down 100% → 88%, 8 alarms) |
| BTC | 1.64 → 1.34 (−18%) | 70% → 66% | 93% → 91% | 93% → 97% | yes | **no (rate −42%, blind 54% → 41%)** |
| ETH | 1.58 → 1.31 (−17%) | 69% → 65% | 88% → 81% | 92% → 91% | yes | yes |
| XRP | 1.70 → 1.26 (−26%) | 71% → 65% | 92% → 96% | 96% → 94% | yes | yes |
| SOL | 1.38 → 1.22 (−12%) | 67% → 64% | 94% → 89% | 98% → 99% | yes | yes |
| DOGE | 1.49 → 1.36 (−9%) | 68% → 66% | 94% → 88% | 97% → 94% | yes | yes |
| Bitstamp BTC | 1.63 → 1.33 (−19%) | 70% → 66% | 88% → 88% | 94% → 94% | yes | **no (rate −31%)** |

* **Plug-in: 6 of 7 assets within all limits** (the declared PASS line is at least 5).
  * In 2026, down alarms fell 11% to 31% per watched minute. Up alarms changed by +1% to −21%.
  * Every 2026 crypto value moves the same way: fewer alarms and less blind time. So 2026 was a year with fewer volatility shifts, not a year where the monitor broke.
* **cs: 4 of 7 (FAILED).**
  * On BTC the drop comes from one segment. From 16 August to 4 October 2026 (49 days), the increase side could not alarm (tau² ≤ m_up, E19's blind mode). This matches E23's 47-day segment.
  * Bitstamp BTC/USD was in the same mode from 16 August to 14 September 2026 (29 days).
  * Before 2026, the longest such stretch was 7.2 days (BTC, September 2025). No other coin had one.

### Q3. Does "blind about 70% of the time" hold, and how different is cs? Crypto VERIFIED, Nifty FAILED

How to read this: share of minutes spent calibrating, lowest and highest year, and over the whole run. "cs blind for increases" also counts the watched minutes in which cs cannot raise an up alarm. "Longest stretch" is the longest unbroken run in which cs is blind for increases. Declared line: an asset-year agrees with the README if the plug-in share is 60% to 80%.

| Asset | Plug-in blind share by year (all years) | Years in 60% to 80% | cs calibrating (all years) | cs blind for increases (all years) | Longest cs stretch blind for increases |
|---|---|---|---|---|---|
| Nifty 50 | 43% to 65% (52%) | 1 of 10 | 14% to 44% (27%) | same (27%) | 2 sessions |
| BTC | 66% to 71% (69%) | 3 of 3 | 41% to 56% (50%) | 54% to 58% (56%) | **49.0 days** (Aug to Oct 2026) |
| ETH | 65% to 72% (68%) | 3 of 3 | 49% to 56% (53%) | same | 1 day |
| XRP | 65% to 72% (69%) | 3 of 3 | 50% to 57% (54%) | same | 1 day |
| SOL | 64% to 67% (66%) | 3 of 3 | 47% to 50% (49%) | same | 1 day |
| DOGE | 66% to 69% (68%) | 3 of 3 | 49% to 54% (51%) | same | 1 day |
| Bitstamp BTC | 66% to 70% (68%) | 2 of 2 | 44% to 53% (49%) | 53% (53%) | **28.9 days** (Aug to Sep 2026) |

* **Crypto:** all 17 asset-years are between 64% and 72%. The README number is right for crypto (VERIFIED).
* **Nifty 50:** only 1 of 10 years is in range, 52% overall. The general wording at README line 166 does not hold for Nifty (FAILED there). Each Nifty calibration still costs two sessions. Nifty is less blind because it raises fewer alarms per session.
* **Plug-in blind stretches tell us nothing extra.** Each stretch is exactly one calibration window: 1 day on crypto, 2 sessions on Nifty, for both the median and the longest (`blind_stretches.csv`). Calibration never needed a long extension. The stretches that matter are cs's "blind for increases" runs above.
* **cs vs plug-in:** cs is blind 15 to 19 points less than the plug-in on crypto over the whole run, and 25 points less on Nifty, because it alarms less often. It raised 22% to 28% fewer alarms on crypto and 48% fewer on Nifty.

### Q4. Alarms right after a restart against later alarms (up alarms) [REAL]

How to read this: "early" means within 60 watched minutes after a calibration ends. "Late" means more than one day after (1,440 minutes on crypto, 375 on Nifty). p comes from a one-sided permutation test (2,000 shuffles) of "the early hit share is lower". Declared line: early alarms are "less useful" if their hit share is at least 10 points lower with p < 0.05.

| Group | Monitor | Early up alarms (share of up alarms) | Hit share early / late | p | Post / pre early / late | Next hour over calibration RV, early / late |
|---|---|---|---|---|---|---|
| 5 Binance coins | plug-in | 129 (6.9%) | 96% / 93% | 0.90 | 0.75 / 0.75 | 2.7 / 2.0 |
| 5 Binance coins | cs | 84 (5.2%) | 96% / 91% | 0.98 | 0.78 / 0.64 | 3.0 / 2.5 |
| Bitstamp BTC | plug-in | 18 (7.7%) | 83% / 89% | 0.46 | 0.53 / 0.68 | 2.1 / 2.0 |
| Nifty 50 | plug-in | 23 (7.1%) | 91% / 95% | 0.33 | 0.66 / 0.87 | 2.2 / 2.1 |
| Nifty 50 | cs | 7 (3.5%) | 100% / 94% | 1.00 | 0.80 / 0.82 | 2.1 / 2.4 |

* **No group meets the "less useful" line.** Alarms right after a restart are confirmed as often as later ones. Bitstamp cs has only 7 early alarms; it is in `q4_restart_tests.csv`.
* **On Nifty and Bitstamp, early alarms have a lower post / pre ratio** (0.66 against 0.87, and 0.53 against 0.68). So an alarm just after a restart more often comes at the top of a burst that is already fading. This was not tested.
* **Down alarms almost never fire early.** With the plug-in, none of about 2,000 down alarms fired within 60 minutes of a restart. With cs, 1 of 1,081 on crypto did. A decrease alarm needs many small returns in a row, so it cannot fire this fast.

## 4. README and docs claims that need changing

| File, line | Claim now | What R4 found | Suggested change | Severity |
|---|---|---|---|---|
| README.md, 166 | "blind for about 70% of the time on real data" | Crypto 64% to 72% in every asset-year (holds). Nifty 50 43% to 65%, 52% over 2017 to 2026 | "about 70% of the time on crypto data and about half the time on the Nifty 50" | Medium (the wording covers Nifty and is wrong there) |
| docs/limitations.md, 11 | cs blind mode: E19's 5 cases "all on BTC", and E23's 47-day BTC segment. The E19 section of the experiment report says it "was not checked" whether these start on the quiet BTC tick stretches | Bitstamp BTC/USD had a 28.9-day blind-for-increases stretch (16 Aug to 14 Sep 2026), at the same time as Binance BTC's 49-day one. On BTC in 2026, cs was blind for increases 58% of all minutes | Add the Bitstamp case and say the mode is not specific to one venue | Medium |
| README.md, 37 and docs/limitations.md section 2 | Alarms are "followed by volatility on the flagged side" | True against the calibration level: 81% to 98% (up) and 82% to 100% (down) per asset-year. But in 61% to 68% of cases the hour after an up alarm is calmer than the hour before, partly by construction | Add one sentence: an up alarm means volatility is above the calibrated level, not that it will keep rising | Low (clarification) |
| README.md, 32 and PITCH_BAPAT.md, 58 | "on real crypto data it is blind about 70% of the time" | Holds in all 17 crypto asset-years, including fresh 2026 data (64% to 66%) | None needed. Could add "64% to 72% in every year 2024 to 2026" | None |

## 5. What this does not show

* **No ground truth.** Nobody knows when volatility "really" changed. All usefulness checks compare with the calibration level or with random times. Volatility clusters, so any detector that reacts to a burst will look good on the hit share. The random-time control removes part of this effect, not all of it. No strong forecast or naive rule was run here. E6, E19 and E20 did that for earlier data.
* **p-values are optimistic.** Alarms close in time are not independent, but the placebo draws treat them as if they were. The hit-share p-values sit at the floor (0.0005) and would stay small under any sensible correction. Still, the exact values should not be quoted.
* **The post / pre metric is biased by construction** (see the post hoc simulation). Its failure does not show that alarms are useless. It shows that they do not predict a further rise.
* **No false-alarm statement.** Real minute returns are not a null, so every rate here is an alarm rate, not a false-alarm rate.
* **Short years.** 2024 starts on 21 January, 2026 ends on 24 October, Bitstamp starts on 28 January 2025, and Nifty 2017 starts in April. Cells with fewer than 10 alarms (several Nifty and cs down cells) are noisy.
* **Narrow coverage.** One venue per coin, one public mirror per data set, and one setting per monitor (band 0.25). Other bands, the awake monitors of E14 and the guard of E24 were not run.
* **Post hoc items** (not in PLAN.md): the synthetic post / pre check, the cs blind-for-increases stretches, the alarms-per-day line, and the early / late post / pre comparison.

**Files in `reports/round4/R4/`:**
* Plan and report: `PLAN.md`, `REPORT.md`.
* Main tables: `asset_year.csv` (full per asset-year table, both monitors), `alarms.csv` (every alarm with its windows), `asset_summary.csv`.
* Blind stretches: `blind_stretches.csv`, `blind_for_up_stretches.csv`.
* Per-question results: `q1_stability.csv`, `q2_2026.csv`, `q3_blind.csv`, `q4_restart_bins.csv`, `q4_restart_tests.csv`, `posthoc_selection.csv`, `verdicts.txt`.
* Data hashes: `fetched_2026_files.csv`.
* Figures: `fig_blind_share.png`, `fig_post_pre_ratio.png`.
* The scripts and `logs/`.

Note from the coordinator: the agent's file-writing tool refused to create this report, so the coordinator saved the agent's text here unchanged apart from the sentence about the run arrays in section 2.
