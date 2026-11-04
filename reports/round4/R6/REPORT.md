# R6. Is there a faster setting that keeps its false-alarm promise on REAL calm days?

Agent R6, round 4, 4 November 2026. **This question is post hoc**: it was chosen after R2 found that the README setting raises almost no false alarms and is 3 to 6 times slower than it needs to be. Status: **VERIFIED**. Every planned run finished in full, with twice E19's stream count. All outputs are in `reports/round4/R6/`. Nothing outside this folder was written.

Words used here: a **false alarm** is an alarm on a stream where volatility never changes. The **band** (null_band) is the share of change the monitor ignores: with band 0.1 it only bets on a variance that moved by more than 10%. The **plug-in** baseline estimates the calm level from the calibration window and has no proof. The **cs** baseline (confidence sequence, delta 0.01) has a proof. A **Wilson 95% interval** is the range of rates the data cannot rule out. A rate is **above** its bound when the whole interval lies above it (E19's rule).

## 1. Summary in plain language

* **Band 0.1 and band 0.05 keep the promise on real calm days, with room to spare.** [REAL marginal, constructed streams] On 40,000 permuted calm-day streams pooled over five coins, band 0.1 alarmed within the first day on 0.27% of streams (interval 0.22% to 0.32%, bound 5%) and within five days on 1.41% (1.30% to 1.53%, bound 25%). Band 0.05: 1.14% and 8.73%. Nifty 50 (8,000 streams): 0.34% and 2.81% for band 0.1; 1.07% and 10.66% for band 0.05. VERIFIED.
* **Band 0 breaks the promise, as E19 found.** 5.61% within a day (5.39% to 5.84%, above 5%) and 41.94% within five days (above 25%) on the five coins. On Nifty 2.88% within a day (holds) and 31.35% over five days (above). VERIFIED.
* **The README setting (band 0.25) and both cs settings almost never alarm.** README: 1 alarm in 40,000 streams within a day (0.00%), 9 within five days (0.02%). cs band 0: 0.03% and 0.13%. cs band 0.1: 0.00% and 0.01%. VERIFIED.
* **The faster bands buy real speed.** [SYNTHETIC] Median steps to detect a doubling of variance (Gaussian, 2,000 paths): README 115, band 0.1 70, band 0.05 59, band 0 50. For a halving: 262, 154, 132, 113. On Student t3 data: doubling 281, 140, 112, 89; halving 513, 214, 171, 139. So band 0.1 is 1.6 to 2.4 times faster than the README and band 0.05 is 1.9 to 3.0 times faster. VERIFIED.
* **The proved cs baseline at band 0 is as fast as today's README setting and holds with a huge margin.** Median delays 114, 195, 262 and 341 steps on the four cells against 115, 262, 281 and 513 for the README (at the same A = 7,800). cs with band 0.1 is slow (177 to 601 steps) and misses 6% to 8% of Student t3 changes. VERIFIED.
* **What this means for the README.** Band 0.1 is a defensible new default for the approximate (plug-in) setting: on real calm days it uses about a twentieth of the one-day budget and is 1.6 to 2.4 times faster. For a setting with a proof, cs band 0 costs nothing in speed against the current default. See section 4 for the caveats.

## 2. What was run

**Pools and streams.** `r6_pools.py` calls E19's own `prepare_pools` (unchanged code from `experiments/e19_real_null.py`), with `E19_OUT` and `E19_SCRATCH` pointed at `reports/round4/R6/pools/`, so `results/` was not touched. The pools reproduce `results/e19_pools.csv` exactly (BTC 81 calm days and 116,640 returns, ETH 95 and 136,800, XRP 95 and 136,800, SOL 100 and 144,000, DOGE 94 and 135,360, Nifty 50 639 and 238,522; same tau, m0 and B). Streams are E19's `StreamGen("perm")`: each stream is the start of a uniform random shuffle of a pool (sampling without replacement), E19's primary stream type. E19's full pools were used, no subset.

**Settings.** E19's data-set settings: crypto day 1,440 minutes, calibration 1,440, A = 28,800 for the plug-in and 36,000 for cs; Nifty day 375, calibration 750, A = 7,500 and 9,375. Six settings on the same streams: plug-in with band 0.25 (README), 0.1, 0.05, 0; cs (delta 0.01) with band 0 and 0.1. Each stream runs until its first alarm or 5 monitored days (`horizon_only=True` in E19's vectorised `x_engine_multi`, which E19's tests check alarm for alarm against `VolatilityMonitor`). Bounds from E19's `bounds()`: 5% within a day for all; 25% within five days for the plug-in (nominal, no proof) and 21% for cs (0.01 + 5 x 1,440 / 36,000).

**Streams.** 8,000 per data set (8 batches of 1,000), twice E19's 4,000: 40,000 pooled over the five coins, 8,000 for Nifty. No stream ended early.

**Delay grid.** `r6_delay.py`: paths from `avmonitor.simulate.change_returns` as in `reports/stress/b2_power_delay.py` (calibration 1,000, change at monitored step 400, horizon 1,950 steps), variance x2 and x0.5, Gaussian and Student t3 (`heavy_tail_df=3`, unit variance), 2,000 paths per cell (the brief allowed 300). Thresholds: A = 7,800 for the plug-in (B2 and README); cs at 7,800 (as B2) and at 9,750 = 390 / (0.05 - 0.01) (E19's rule). Same vectorised engine.

**Check.** `r6_validate.py`: the engine against the `VolatilityMonitor` class at band 0.1, band 0.05 and cs band 0.1 (bands E19's tests did not cover), 10 Student t3 paths with a change, thresholds 50, 7,800 and 9,750: 0 mismatched first alarms in 50 (`validate.log`). Cross-checks with committed work: E19's pooled plug-in band 0 was 5.50% and 41.0% (here 5.61% and 41.9%, fresh seeds); E19's cs band 0 was 0.01% and 0.13% (here 0.03% and 0.13%); R2's README delay for Gaussian x2 was 114 and B2's 112 (here 115); B3's band 0.1 five-day rate on synthetic data was 1.2% (here 1.4% on real calm days).

**Seeds** (numpy `default_rng`, all inside 9,970,000 to 9,979,999): streams 9,970,000 + 100 x data set + batch (data sets 0 to 5 in E19's order, batches 0 to 7); delay cells 9,975,000 + 10 x law + factor (law 0 Gaussian, 1 t3; factor 0 for x2, 1 for x0.5); validation 9,979,999.

**Commands** (from `reports/round4/R6`, venv active, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`): `python r6_pools.py`, `./run_null.sh` (48 calls of `python r6_null.py <dataset> <batch> 1000`, 3 to 17 s each), `python r6_delay.py 2000 r6_delay.csv` (60 s), `python r6_validate.py 10`, `python r6_summarise.py`, `python r6_tables.py`. Total about 12 minutes on one core.

**Files.** `r6_null_rows.csv` (every data set, setting, interval, bound and seed list), `r6_delay.csv`, `r6_tables.md` (all tables, per coin too), `null_batches/*.pkl` (per-stream first alarms), logs `pools.log`, `null.log`, `delay.log`, `validate.log`, `summarise.log`.

## 3. Tables

### 3a. False alarms on permuted real calm days [REAL marginal, constructed streams]

How to read this: each row is one setting on streams of shuffled real calm-day returns, where any alarm is false. "Day 1" is the share of streams with an alarm within the first monitored day, "5 days" within the first five. A rate "holds" unless its whole Wilson interval lies above the bound. Per-coin tables are in `r6_tables.md`; the five coins agree with each other to within their intervals.

**Five coins pooled, 40,000 streams per setting.** Bounds: 5% (day 1); 25% plug-in and 21% cs (5 days).

| Setting | A | Day 1 [Wilson 95%] | 5 days [Wilson 95%] | Verdict |
|---|---|---|---|---|
| plug-in, band 0.25 (README) | 28,800 | 0.00% [0.00%, 0.01%] (1 alarm) | 0.02% [0.01%, 0.04%] (9) | holds |
| plug-in, band 0.1 | 28,800 | 0.27% [0.22%, 0.32%] (106) | 1.41% [1.30%, 1.53%] (563) | holds |
| plug-in, band 0.05 | 28,800 | 1.14% [1.04%, 1.25%] (457) | 8.73% [8.46%, 9.02%] (3,494) | holds |
| plug-in, band 0 | 28,800 | 5.61% [5.39%, 5.84%] (2,245) | 41.94% [41.46%, 42.43%] (16,777) | above on both |
| cs baseline, band 0.1 | 36,000 | 0.00% [0.00%, 0.01%] (0) | 0.01% [0.00%, 0.02%] (3) | holds |
| cs baseline, band 0 | 36,000 | 0.03% [0.01%, 0.05%] (10) | 0.13% [0.10%, 0.17%] (51) | holds |

**Nifty 50, 8,000 streams per setting.** Same bounds.

| Setting | A | Day 1 [Wilson 95%] | 5 days [Wilson 95%] | Verdict |
|---|---|---|---|---|
| plug-in, band 0.25 (README) | 7,500 | 0.01% [0.00%, 0.07%] (1) | 0.06% [0.03%, 0.15%] (5) | holds |
| plug-in, band 0.1 | 7,500 | 0.34% [0.23%, 0.49%] (27) | 2.81% [2.47%, 3.20%] (225) | holds |
| plug-in, band 0.05 | 7,500 | 1.07% [0.87%, 1.33%] (86) | 10.66% [10.00%, 11.36%] (853) | holds |
| plug-in, band 0 | 7,500 | 2.88% [2.53%, 3.26%] (230) | 31.35% [30.34%, 32.38%] (2,508) | day 1 holds, 5 days above |
| cs baseline, band 0.1 | 9,375 | 0.00% [0.00%, 0.05%] (0) | 0.00% [0.00%, 0.05%] (0) | holds |
| cs baseline, band 0 | 9,375 | 0.01% [0.00%, 0.07%] (1) | 0.05% [0.02%, 0.13%] (4) | holds |

### 3b. Detection delay of the settings that hold [SYNTHETIC]

How to read this: 2,000 simulated paths per cell, 1,000 calibration steps, variance change at monitored step 400. "Median delay" is the middle number of steps from the change to the alarm among paths detected within 1,950 steps; the share detected is in brackets when it is below 99%. "Early" is the share of paths that alarmed before the change (a false alarm; it is larger for smaller bands). The plug-in rows use A = 7,800; cs rows use A = 7,800 as in B2, with E19's A = 9,750 in brackets. Band 0 is shown for scale only: it fails 3a. Full rows with quartiles and means: `r6_tables.md`, `r6_delay.csv`.

| Setting | Gaussian x2 | Gaussian x0.5 | Student t3 x2 | Student t3 x0.5 | Early (range over cells) |
|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 115 | 262 | 281 | 513 (98.3%) | 0.0% |
| plug-in, band 0.1 | 70 | 154 | 140 | 214 | 0.1% to 0.4% |
| plug-in, band 0.05 | 59 | 132 | 112 | 171 (98.6%) | 0.7% to 1.4% |
| plug-in, band 0 (fails 3a) | 50 (97.9%) | 113 (97.9%) | 89 (96.7%) | 139 (96.9%) | 2.1% to 3.4% |
| cs baseline, band 0 | 114 (119) | 195 (203) | 262 (275) | 341 (364) | 0.0% |
| cs baseline, band 0.1 | 177 (189) | 270 (280) | 516 (93.9%) (540) | 601 (93.4%) (637) | 0.0% |

Speed against the README setting, from the medians: band 0.1 is 1.6, 1.7, 2.0 and 2.4 times faster on the four cells; band 0.05 is 1.9, 2.0, 2.5 and 3.0 times faster; cs band 0 is 1.0, 1.3, 1.1 and 1.5 times faster.

## 4. Should the README recommend a different default?

**Yes, with one of two choices, and the choice depends on what the README wants to promise.**

* **If the README keeps the plug-in baseline as the default, band 0.1 is the better default.** On 40,000 real calm-day streams it used 0.27% of its 5% one-day budget and reached 1.4% (crypto) to 2.8% (Nifty) over five days against the nominal 25%. That is still a wide margin, and the same margin as band 0.25 in kind, only smaller. In return it detects a change 1.6 to 2.4 times sooner than band 0.25, and 2.0 to 2.4 times sooner on heavy-tailed data, where the README setting is slowest. Band 0.05 is faster still and also holds, but with less room: its five-day rate is 6 times that of band 0.1 (8.7% against 1.4%), and the rate climbs steeply from band 0.1 to 0.05 to 0 (1.4%, 8.7%, 42%). A data set a little noisier than these calm days could push band 0.05 over the nominal five-day line while band 0.1 stays well under it. So band 0.1 is the safer of the two.
* **If the README wants its default to have a proof, cs band 0 should be the default, and it costs nothing in speed.** It is the proved option when the baseline is estimated, it stayed at 0.03% (day 1) and 0.13% (five days) on real calm days, and its delays match or beat the current README setting on every cell (114 against 115 steps for a Gaussian doubling, 341 against 513 on Student t3 halving). Adding a band to cs only slows it (band 0.1: 177 to 601 steps) and makes it miss 6% to 8% of Student t3 changes, so cs should stay at band 0.
* **Band 0 for the plug-in should not be the default**: it is above both bounds on the five coins and above the five-day bound on Nifty, which repeats E19's finding with fresh seeds.
* **What should stay in the README.** The plug-in has no proof at any band; this test cannot give it one. The band also has a second job that this test does not measure (section 5), so the README should keep saying that 0.25 was chosen for real, non-calm monitoring in E6 and that smaller bands raise more alarms on real days.

## 5. What this does not show

* **Post hoc.** The question and the bands 0.1 and 0.05 were chosen after R2's result and after E19's committed results were known. The seeds are fresh, but the choice of what to test was not blind.
* **Only calm days, and only their shuffled versions.** The streams keep the real marginal law (zeros, ticks, fat tails) of days whose robust volatility was within 20% to 25% of the trailing median. Shuffling removes all clustering within a day. E19 showed that block-resampled real days break every setting, including the README one, at about the real alarm rate. So "holds on real calm days" does not mean "5% false alarms on a real exchange day". On real days the band's main job is to absorb ordinary within-day volatility movement, and a smaller band will raise more alarms there (E6, E8). That cost was not measured here.
* **Thresholds differ between the two parts.** The null test used E19's A = 28,800 with 1,440-minute days and 1,440-step calibration (Nifty: 7,500, 375, 750). The delay grid used B2's A = 7,800 with 1,000-step calibration. Both are "5% per day", and B3 found the same 1.2% five-day rate for band 0.1 at A = 7,800 on synthetic data, but the delays at A = 28,800 were not measured and will be longer in steps.
* **Nominal bounds only for the plug-in.** The 5% and 25% lines for the plug-in are what the known-baseline proof would give; the plug-in never had that proof, and the five-day rate of band 0 shows the gap. The cs rows are the only ones with a bound that is proved.
* **One change time, two change sizes.** Delays were measured at monitored step 400 for x2 and x0.5 only. B2 and B3 show the README setting is far slower for x1.5 and x2/3; the faster bands were not measured there.
* **Synthetic delays, Gaussian and t3 only.** No jumps, no GARCH, no seasonality in the delay grid.
