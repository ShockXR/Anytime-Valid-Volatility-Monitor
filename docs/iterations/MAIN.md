# Iterations by the integrating session (round 2, night of 2 to 3 November 2026)

These are merged into `ITERATION_LOG.md`. Times are UTC.

---

ITERATION NUMBER: R2-1 (real data expansion)

OBJECTIVE: Get real one-minute data for the four coins of the paper that were missing (ETH, XRP, SOL, DOGE) and for Nifty 50.

PROBLEM IDENTIFIED: Round 1 had only Bitstamp BTC/USD. Every exchange API (Binance, Coinbase, Kraken, Bitstamp, Bybit, OKX, KuCoin, Gate, Gemini), Yahoo, Kaggle, Hugging Face, Zenodo, Stooq and www.nseindia.com were refused by the cloud network policy (403 on CONNECT).

HYPOTHESIS: Public GitHub repositories mirror these data, and raw.githubusercontent.com is allowed.

CHANGES IMPLEMENTED: Searched GitHub (repository search API). Found per-coin Binance spot 1-minute archives (github.com/parshant231234-eng/<SYMBOL>-Historical-Data, one CSV per month) and a cleaned Nifty 50 1-minute archive (github.com/technovusin/nifty50-historical-data, April 2017 to October 2026). Wrote `data.load_binance` (missing minutes become NaN rows, never zero returns), `data.load_nifty` (returns within sessions only, a missing minute drops the two-minute return), `regular_sessions`, and `scripts/download_data.py`, which writes `data/MANIFEST.csv` with SHA-256 hashes.

EXPERIMENTS PERFORMED: Offline unit tests on synthetic files (gaps, duplicates, overnight gaps, missing minutes). A first comparison of Binance quartiles and means with the paper's Table 5.

RESULTS: 140 files, 803 MB of crypto and 55 MB of Nifty. No missing minutes or duplicates in Binance Jul to Dec 2025. Nifty: 2,369 sessions, median 375 bars, 6 short sessions. Binance quartiles match the paper for BTC and XRP closely, the means match for four coins, the paper's ETH mean has the opposite sign, and the paper's DOGE (Q1 = Q3 = 0) cannot come from Binance DOGE (3.2% zero returns). Followed up in E6.

IMPROVEMENTS: All five coins of the paper and an Indian index are now testable. The data are third-party mirrors, so the loaders report every check and the manifest pins the bytes.

REMAINING ISSUES: Mirrors are not the exchanges themselves. India VIX daily data are tried from archives.nseindia.com in E7.

NEXT ITERATION: Edge cases that real data exposes.

---

ITERATION NUMBER: R2-2 (framework iteration 3: numerical stability and edge cases)

OBJECTIVE: Make the estimators and the online scale safe on zero-heavy and gappy real data.

PROBLEM IDENTIFIED: Probing with synthetic edge cases showed four failures. (1) With more than half of the returns exactly zero (coarse ticks, as in the paper's DOGE, or frozen prices such as the NSE outage of 24 February 2021), the MAD is 0, every return is clipped to 0, and every truncated estimator silently returns 0. (2) `VolatilityMonitor` raised `ZeroDivisionError` in that case. (3) `ew_quantile_scale` counted NaN returns as zeros, and could underflow to exactly 0 and never recover. (4) If its first 50 values were zeros it started at 1e-300 and needed hundreds of thousands of steps to climb back.

HYPOTHESIS: A warning plus a tick-size floor fixes (1); extending the calibration window fixes (2); NaN skipping and a positive floor fix (3) and (4), without changing any result on clean data.

CHANGES IMPLEMENTED: `robust.DegenerateThresholdWarning` and `robust.tick_floor(prices)`; `ew_quantile_scale` skips NaN, never goes below `max(floor, 1e-300)`, and starts from the median of non-zero values when needed. The monitor fix (2) was assigned to the agent that owns `monitor.py` (see `docs/iterations/E8.md`).

EXPERIMENTS PERFORMED: `tests/test_edge_cases.py` (7 tests), including a regression test that the rewritten online scale gives bit-identical output to the old update rule on clean data.

RESULTS: All tests pass. One test first failed for an instructive reason: with 60% zeros the true median of |r| is 0, so the tracker legitimately drifts toward 0, and the multiplicative update can climb back by at most a factor e every 562 steps (half-life 390). Recovery within a day therefore needs a floor. The test now checks no underflow without a floor and fast recovery with a tick floor.

IMPROVEMENTS: No silent zeros, no crash path in the estimators, documented recovery speed.

REMAINING ISSUES: Choosing the floor needs the price tick, which `tick_floor` estimates from the data.

NEXT ITERATION: A demonstration and end-to-end tests.

---

ITERATION NUMBER: R2-3 (framework iteration 8: architecture and usability)

OBJECTIVE: Make the prototype executable from a clean environment by someone who has not read the code.

PROBLEM IDENTIFIED: No exact dependency list, no command-line entry point, no test that runs the whole pipeline.

HYPOTHESIS: A small CLI and file-to-alarm integration tests catch interface mistakes that unit tests miss.

CHANGES IMPLEMENTED: `requirements.txt` with exact versions, `python -m avmonitor simulate | run | estimators` (also installed as `avmonitor`), `tests/test_integration.py` (loader, deseasonalisation, monitor, CLI, module entry point).

EXPERIMENTS PERFORMED: Ran the CLI on synthetic data, on Binance ETHUSDT for October 2025 (20 alarms in 31 days) and on Nifty for May and June 2024, where it alarmed at 09:24 IST on 4 June 2024, the 2024 general election counting day, with next-session volatility 5.8 times the baseline.

RESULTS: The first version of the integration test failed. The cause was in the test fixture: under pandas 3 a `DatetimeIndex` is not stored in nanoseconds, so `asi8 // 10**6` gave wrong millisecond timestamps, every minute became missing, and the monitor then crashed on an all-NaN calibration window. The fixture now builds timestamps explicitly; the crash path is the same one fixed in R2-2.

IMPROVEMENTS: One-command demonstrations on synthetic and real data; 48 tests in total after this round.

REMAINING ISSUES: A notebook demonstration is planned for the end of the round.

NEXT ITERATION: Test the false alarm guarantee at scale.

---

ITERATION NUMBER: R2-4 (framework iterations 6 and 7: stress testing, larger simulations, profiling)

OBJECTIVE: Check Theorem B at a scale where its tightness can be measured, under distributions it should cover and under one it does not cover.

PROBLEM IDENTIFIED: E4 used 300 paths per cell and estimated baselines, so it could not separate "the bound is loose" from "the baseline is off". It also never tried a volatility model with persistent clustering.

HYPOTHESIS: With the true baseline, the false alarm rate is below m / A for any distribution with the right truncated mean, and roughly the same across distributions; under GARCH clustering the null is false and the rate exceeds the budget.

CHANGES IMPLEMENTED: `experiments/e09_stress_guarantee.py` (4000 paths per distribution, 2000 GARCH paths per cell, profiling).

EXPERIMENTS PERFORMED: Q9a to Q9d, see `docs/report_parts/E9.md`.

RESULTS: One-day false alarm rates of 1.1% to 1.4% against the 5% bound for Gaussian, Student t3, jump and tick-rounded returns, so valid, distribution free and about 4 times conservative. Under GARCH with half-life 69 steps or more, no band up to b = 1 restores 5% (12% to 20% at b = 1). The streaming monitor runs at 10 microseconds per step; no compiled code is needed. A first run of the GARCH part crashed when two parameter pairs had the same persistence (a pivot on duplicate keys); the grid was fixed to four distinct persistences.

IMPROVEMENTS: A precise statement of what the guarantee covers and what it does not, and a reason to restate the null relative to a forecast.

REMAINING ISSUES: The conservativeness is examined theoretically by the theory track (E11).

NEXT ITERATION: Forecast-relative monitoring (E10).

---

ITERATION NUMBER: R2-5 (framework iterations 4 and 6: statistical validation and sensitivity, forecast-relative null)

OBJECTIVE: Restore a meaningful false alarm guarantee under persistent volatility clustering, and measure what that costs.

PROBLEM IDENTIFIED: E9 showed that under GARCH clustering with a half-life of an hour or more, the raw-return monitor alarms on 12% to 100% of days whatever the band, because the conditional mean of the statistic is not constant.

HYPOTHESIS: Dividing returns by a predictable volatility forecast makes the conditional variance constant under the model, so the guarantee applies again. A forecast that reacts to returns will also absorb part of any real break, so power should fall.

CHANGES IMPLEMENTED: `forecast.py` (predictable EWMA and GARCH(1,1) filters, a quasi maximum likelihood GARCH fit with variance targeting, `standardise`) with tests; `experiments/e10_forecast_relative.py` (Q10a to Q10d).

EXPERIMENTS PERFORMED: Q10a validity under three GARCH persistences (600 paths each); Q10b power against scale and intercept breaks that double the long-run variance (300 paths, 5 days); Q10b2 the same monitors without a break, and with a 4-fold break; Q10d the exact signal that survives standardisation and the Kullback-Leibler rate; Q10e the transient right after a break; Q10c five Binance coins, July to December 2025.

RESULTS: A fitted GARCH forecast gives 3% to 6% one-day false alarms against 82% to 100% for raw returns. EWMA forecasts work only when their half-life matches the clustering. Breaks are detected within five days by the forecast-relative monitor at band 0 (93% to 94%) but slowly (median about 1600 to 1800 steps, more than a day), and never at band 0.25. Q10d explains why: scaling GARCH returns by c is the same model as multiplying the intercept by c², and a frozen pre-break filter absorbs the break except for a constant offset, so the visible ratio is ρ(k) = k(1 - β)/(1 - α - β + αk), capped at (1 - β)/α (1.2 at α = 0.05, β = 0.94; a doubling shows as +9%). The KL rate then puts any detector's delay at about 1.9 days for a doubling at persistence 0.99 (Lai 1998, asymptotic). Without a break the same monitors alarm within five days on 38% to 44% of paths, so part of that detection rate is the plug-in baseline: for exactly i.i.d. standardised returns the five-day false alarm rate was 44% against the 25% that Theorem B guarantees with the true baseline. A short transient after the break carries extra evidence (a scale break is fully visible at first and fades with a half-life of about 11 steps), which matters mainly for large breaks. On real data, dividing by an EWMA forecast changed alarm rates by between +2% and -17%; alarms remained informative relative to the forecast (next-day realised variance on the flagged side of the forecast 60% to 87% of the time, against 37% to 63% at random times).

IMPROVEMENTS: A clear statement of two different monitoring questions (any change in volatility, or a departure from a model) with a guarantee for each, and a closed-form explanation, checked by simulation, of why breaks are hard to see relative to a persistent model.

REMAINING ISSUES: The forecast-relative guarantee is only as good as the model; a misspecified forecast (EWMA with the wrong half-life) is not valid. The information bound is asymptotic. A forecast with an explicitly modelled long-run level, re-estimated slowly, is a possible middle ground that we have not tested.

NEXT ITERATION: Integrate the agents' results (E6, E7, E8, E11, alternatives) and the independent review.
