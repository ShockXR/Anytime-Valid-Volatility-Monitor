# Pre-registration E19: the guarantee on real data

*Written and committed before any real-data part of `experiments/e19_real_null.py` was run: before the calm days were identified, before any pool or stream was built from real returns, and before any monitor was run on real or real-based streams for this study. The commit time is the record of when the design was fixed. Anything decided after this commit will be marked as post hoc in the report.*

**Why.** The referee report on the paper (`docs/review/referee_report_paper.md`, major comment M1) says that the real-data evidence never exercises the anytime-valid guarantee: every main real-data run uses the plug-in baseline, which has no proof, and on real days the null is false most of the time. M1 asks for (a) the confidence-sequence baseline (Theorem D, `docs/math.md` section 6) as the primary real-data configuration, (b) null streams with the real marginal law, built by permuting returns from calm periods, and block-bootstrapped streams to show what clustering does, and (c) the forecast-relative null on real data with its cost. Minor comment m32 asks for the Gaussian row of Table S2 (the split clip, E11f) with the confidence-sequence baseline.

## 1. Data and calm stretches

**Crypto.** The five Binance coins of E6 (BTCUSDT, ETHUSDT, XRPUSDT, SOLUSDT, DOGEUSDT), deseasonalised exactly as in `experiments/e06b_five_coins_monitor.py`: `prepare(coin, "Jul-Dec 2025")` (history from 1 January 2024, returns divided by the predictable intraday profile with a 20-day lookback). Period: the 184 UTC days from 1 July to 31 December 2025.

**Nifty 50.** `load_clean("2021-10-01", "2026-11-01")` and `deseasonalise` of `experiments/e07_nifty.py` (E7's session profile). Period: the sessions from 1 January 2022 to 31 October 2026 (E7's held-out period); the three months before only warm up the profile and the trailing medians.

**Calm days (sessions).** For each day d (UTC day for crypto, session for Nifty), let v_d be the robust volatility of the deseasonalised returns of that day: the scaled median absolute value (1.4826 × median |u|) of its finite returns, defined when the day has at least 1000 finite returns (crypto, as E6's `MIN_FINITE`) or 300 (Nifty). Let R_d be the median of v over the 20 previous days (sessions) on which v is defined (these may lie before the period). **Day d is calm if 0.8 ≤ v_d / R_d ≤ 1.25.** No alarm, no monitor output and no event date enters this definition. All calm days of the period are used; no day is dropped for any other reason.

**Pool.** For each data set, the finite deseasonalised returns of all calm days, concatenated in calendar order (so a day keeps its internal order, which only the block bootstrap uses). Reported per pool: number of calm days and returns, share of zero returns, the known baseline (τ = 3 × scaled MAD of the pool, m0 = mean of min(u², τ²) over the pool), B = τ² / m0, the share of returns clipped at τ.

## 2. Part 1: null streams with the real marginal law (M1(b))

**Stream types.** For each pool of size N:

* **perm (primary): exchangeable.** A stream is the start of a uniform random permutation of the pool (sampling without replacement). It keeps the real zeros, ticks and jumps and removes all time dependence. Given the returns drawn so far, the rest of the stream is again a uniform random permutation of the rest of the pool, so the rank detector's guarantee holds exactly, including after restarts (`docs/theory_running_tau.md`, Remark B.3). For the known baseline and the confidence-sequence baseline the null holds up to a finite-population effect: E[X_t | past] is the mean of the not-yet-drawn part of the pool, a martingale that starts at m0. Its largest relative deviation from m0 along each stream will be computed and reported.
* **iid: exactness check.** Draws with replacement from the pool: i.i.d. returns from the pool's empirical law. On these streams the null holds exactly for (i), (iii) and (iv) below (Theorem B, Corollary 6.5, Theorem B.2 of the theory note).
* **sb30, sb390, sb1440: real clustering.** Stationary bootstrap (Politis and Romano 1994) of the pool in its calendar order, wrapped circularly, with mean block lengths 30, 390 and 1440 minutes (the same lengths, in returns, for Nifty). Within a block the real dependence is kept, so the conditional mean of X moves, and the nulls of the theorems do not hold. These streams show what real clustering does to each baseline; their rates are compared with the same nominal bounds, but no bound is claimed.

**Monitors** (default detector everywhere: c = 3, six increase bets 0.02 to 0.5, five decrease bets 0.1 to 0.95 of 1 / (B - 1), Shiryaev-Roberts statistics, each side alarming at 2A; settings of E6, E8 and E14, nothing tuned here):

| | Monitor | Calibration | A (crypto / Nifty) | After an alarm |
|---|---|---|---|---|
| (i) | known baseline: τ and m0 of the pool | none; the first calibration-length returns of a stream are skipped so that all monitors start at the same step | 28,800 / 7,500 | restart the statistics at 0 at once (Corollary C) |
| (ii) | plug-in baseline (`VolatilityMonitor`, default) | 1440 / 750 returns | 28,800 / 7,500 | recalibrate on the next 1440 / 750 returns |
| (iii) | confidence-sequence baseline (`VolatilityMonitor(baseline="cs")`, δ = 0.01, τ from the first 25% of the window, defaults otherwise) | 1440 / 750 returns | 36,000 / 9,375 | recalibrate on the next 1440 / 750 returns |
| (iv) | rank detector of E14, route B, LR bet family (`ConformalMonitor`), restarted at once (awake-B) | the first 1440 / 750 returns form the reference set | 28,800 / 7,500 | empty the reference set and bet from the next return |

Day length: 1440 minutes (crypto), one session of 375 minutes (Nifty). Calibration 1440 (crypto, E6) and 750 (Nifty, E7's pre-registered primary). Null band b = 0 (primary) and b = 0.25 (secondary; E6's real-data setting; for (iv) the band divides each rank payoff by 1 + b, as in E14).

The monitors are vectorised re-implementations of the classes. Before any real-data run they are checked alarm for alarm against `VolatilityMonitor` (plug-in and cs, including the extension rule for degenerate calibration windows) and `ConformalMonitor(restart_len=0)` on synthetic streams (`tests/test_e19_real_null.py`).

**Endpoints.** Each stream runs until it has had m monitored steps, m = 20 days (crypto: 28,800 steps) or 20 sessions (Nifty: 7,500 steps), with the restart rule above. Permutation streams hold at most cal + m + 10 cal returns (bootstrap streams cal + m + 40 cal); a stream that runs out first stops early, and the number of such streams is reported.

1. **One day (primary):** P(first alarm within the first day of monitoring), that is within 1440 (375) monitored steps after the first calibration. Bound: m / A = 5% for (i), (ii, nominal only) and (iv); δ + m / A = 5% for (iii).
2. **Five days:** P(first alarm within 7200 (1875) monitored steps). Bound: 25% for (i), (ii, nominal) and (iv); 0.01 + 7200 / 36,000 = 21% for (iii).
3. **False alarms per monitored day with a restart after each alarm,** over the first m monitored steps. Bound (expected number): m / A, that is 0.05 per monitored day, for (i) (Corollary C) and (iv) (Remark B.3), nominal for (ii); (δ + m / A) / (1 - δ) over m steps for (iii) (Corollary 6.6), that is 0.0409 per monitored day for m = 20 days. Also reported per calendar day (alarms per day of stream used, blind time included).

**Numbers of streams.** perm: 4000 streams per data set; iid, sb30, sb390, sb1440: 2000 per data set. The same streams (the same seeds) for all four monitors and both bands. Pooled over the five coins: 20,000 perm streams (standard error at most 0.35 points for any rate) and 10,000 of each other type (at most 0.5 points). Seeds: `numpy.random.SeedSequence([19, data set index, stream type index, batch index])`; rank tie-breakers from a separate stream of the same sequence.

**Intervals and reading.** Probabilities with 95% Wilson intervals; per-day rates with a 95% normal interval from the per-stream counts (streams are independent). A rate is **above its bound** if the lower end of its interval is above the bound, and **at or below** otherwise. For (i), (iii) and (iv) on perm and iid streams the theory predicts at or below the bound (up to the finite-population effect for (i) and (iii) on perm streams). For (ii) the theory predicts nothing; E14 found 41% to 43% over five days after 390 calibration steps on i.i.d. synthetic data.

**Primary endpoint.** Endpoint 1 (one day), perm streams, band 0, pooled over the five coins, for each of the four monitors, with its Wilson interval and bound. Everything else in Part 1 (five days, per day, per coin, iid streams, bootstrap streams, band 0.25, Nifty) is secondary.

## 3. Part 2: the confidence-sequence baseline as the primary real-data configuration (M1(a))

**Runs.** Exactly E6's monitoring (`e06b_five_coins_monitor.py`, imported, not changed): the five Binance coins and Bitstamp BTC/USD, Jul to Dec 2025, and the five Binance coins out of sample, 21 January 2024 to 30 June 2025; `VolatilityMonitor(threshold=36000, c=3, cal_len=1440, null_band=b, baseline="cs", delta=0.01)` for b = 0, 0.25 and 0.5, recalibrating after each alarm. A = 1440 / (0.05 - 0.01) as in E8.

**Validation.** E6's own functions on the cs alarms: Q6d (`summarise`: share of up alarms followed by a higher and of down alarms followed by a lower robust volatility over the next 1440 minutes than the calibration baseline, Wilson intervals, matched-delay random-time control) and Q6e (`forecast_hits`: next-day robust volatility above / below EWMA and HAR forecasts made at the alarm, random times, the rate-matched ratio rules). Random seeds: 1906 (Q6d) and 1907 (Q6e). The plug-in comparison is E6's published output (`results/e06d_*.csv`, `results/e06e_*.csv`), not a rerun.

**Primary quantity of Part 2.** Pooled over the five Binance coins, Jul to Dec 2025, band 0.25: alarms per day and the validated shares of Q6d (up and down) for cs, next to E6's plug-in. Intervals: Wilson as in E6, and a day-block bootstrap (calendar days resampled with replacement, all coins jointly, 2000 resamples, percentile intervals) for the cs share, the plug-in share and their difference. Everything else is secondary. No p-value will be called significant without a multiplicity statement.

**Monitoring time.** For each coin, period and band: (a) the share of minutes monitored (not calibrating) for cs and plug-in; (b) the median minutes from the end of calibration to an alarm; (c) for plug-in alarms with a same-direction cs alarm within 6 hours, the median of (cs time - plug-in time), and the share of plug-in alarms with no such cs alarm.

## 4. Part 3: the forecast-relative null on real data (M1(c))

**Monitor input.** z_t = u_t / σ̂_t, with σ̂_t² a predictable forecast (`src/avmonitor/forecast.py`, unchanged):

* **GARCH:** `garch11_fit` (Gaussian QML with variance targeting) on each coin's deseasonalised returns of 1 January to 30 June 2025, then frozen; `garch11_filter` run from 1 January 2025 with the fitted unconditional variance as start value.
* **EWMA:** `ewma_variance` with half-life 360 minutes (E10c's variant), run from the first finite deseasonalised return, start value the mean square of its first 1440 finite values.

Monitors on z, Jul to Dec 2025, five Binance coins: plug-in (A = 28,800) and cs (A = 36,000, δ = 0.01), bands 0 and 0.25, calibration 1440, recalibration after each alarm.

**Reported.** Cost: alarms per day, alarms per monitored day against the bound that would hold if the forecast-relative null were true (0.05 plug-in nominal, 0.0409 cs), the share of minutes monitored, the median minutes from the end of calibration to an alarm. Validated shares with E6's functions: Q6d with the calibration baseline taken as the robust volatility of the deseasonalised returns u over the alarm's 1440-minute calibration window (so that it is on the same scale as the next-day robust volatility), and Q6e (EWMA and HAR forecasts of the robust volatility of u, random times, seeds 1908 and 1909). Comparison: the raw deseasonalised monitor (E6 plug-in, Part 2 cs). All of Part 3 is descriptive.

## 5. Part 4: Table S2, Gaussian row, with the confidence-sequence baseline (m32) [SYNTHETIC]

E11f part a, Gaussian law, reproduced with its seed and draw order (`e11f_split_c.py`, SEED = 20261103: the null paths and the three decrease cells are the first Gaussian draws), so the plug-in rows must reproduce Table S2 (median delays 247 / 176 / 185, 158 / 85 / 73, 118 / 58 / 40). The same paths with `baseline="cs"`, δ = 0.01, A = 36,000, c = 3 and c_down = 3, 2, 1.5 (two confidence sequences when c_down ≠ c, each at δ / 2, as in the monitor). Reported: one-day false alarm rate on the 2000 null paths (bound 5%), and the median decrease delay, share detected within 3000 steps, mean delay and paired difference against c_down = 3 on the 600 paths of each cell. The vectorised cs split detector is checked against the class.

## 6. Reporting

* Files: `results/e19_pools.csv`, `results/e19_null_streams.csv` (one row per data set, stream type, monitor and band, with every endpoint, interval and bound), `results/e19_null_pooled.csv`, `results/e19_cs_*.csv` (Part 2), `results/e19_fr_*.csv` (Part 3), `results/e19_split_c_cs.csv` (Part 4), figures `results/figures/e19_*.png`, report `docs/report_parts/E19.md`, iteration notes `docs/iterations/E19.md`.
* Labels: Part 1 streams are built from real returns but are constructed streams; they are labelled **[REAL marginal, constructed streams]**. Parts 2 and 3 are **[REAL]**, Part 4 **[SYNTHETIC]**.
* Every endpoint above is reported, whatever its value, for every data set and stream type that could be built. If a pool cannot be built or a run fails, that is reported instead of a replacement. Nothing in sections 1 to 5 will be changed after the runs; any analysis added later is marked post hoc.

## 7. What was already known before this document

* E6 (plug-in, band 0.25): 0.66 to 0.71 alarms per day on every coin and period; 75% to 84% of up and 83% to 93% of down alarms validated (Q6d); 71% to 81% hit rates against forecasts (Q6e).
* E8: on Bitstamp BTC, July to December 2025, the cs baseline cut alarms by 16% to 20% at every band, mostly decrease alarms, and doubled the share of time monitored at band 0 (15% to 29%).
* E14 and E10: the plug-in exceeds the known-baseline bound over five days on i.i.d. synthetic data (41% to 43% after 390 calibration steps) and over five days for GARCH-standardised returns (44% after 1440); E14 found 0.084 to 0.094 false alarms per day over 30-day i.i.d. streams after 390 calibration steps.
* E10c: monitors on deseasonalised returns divided by an EWMA forecast (half-lives 360 and 1440 minutes, band 0.25, plug-in) raised 0.57 to 0.71 alarms per day on the five coins, July to December 2025.
* E6 notes that on quiet weekends Binance BTCUSDT often sits between two adjacent prices one tick apart.
* While writing this document the vectorised monitors were timed on a synthetic pool (Student t4 returns with 4% zeros, 130,000 values; scratch files outside the repository). No real return was loaded for this study before this commit.
