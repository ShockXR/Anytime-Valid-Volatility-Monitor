# W1 plan (declared before any result was computed)

Written at 02:39 UTC on 4 November 2026, before any W1 analysis was run. All data are [REAL].
Monitor: plug-in (README setting), the alarms and run arrays of R4 (reports/round4/R4/run_monitors.py).
Seeds: numpy default_rng, only in 9,980,000 to 9,989,999.

## Part 1. Volatility-matched placebo (every asset and year, up and down alarms)

Definitions as in R4: pre = realised variance (RV, mean squared deseasonalised return) of the 60 minutes ending
with the minute; post = RV of the next 60 minutes; base = RV of the calibration window of the minute's segment;
windows need at least 50 finite returns. Hit for an up alarm: post > base. Hit for a down alarm: post < base.

* Matching variable: x = pre / base at the alarm minute.
* Pool for one asset-year: watched minutes of that asset and year with valid pre, post and base, alarm minutes removed.
* Primary match: the decile of x, with decile edges taken from the pool of that asset-year. For Nifty also the same
  session hour (7 buckets of 60 minutes from 09:15 IST).
* Declared sensitivity (also before results): the same with 100 percentile bins instead of deciles (a tighter match).
* Randomisation ("permutation") test: for each alarm draw one control minute from its stratum, at random; repeat
  2,000 times. One-sided p = (1 + number of control sets at least as good as the alarms) / 2,001.
* Metrics per asset-year and direction: hit share (alarms vs mean over control sets); mean log(post / base),
  reported as a ratio exp(alarm mean - control mean). For down alarms "better" means lower.
* Seeds: 9,980,000 + 1,000 x asset + 10 x (year - 2017) + direction (0 up, 1 down) + 2 for the percentile version.
* Pass line per cell: hit-share p < 0.05. Verdict per direction over the asset-years with at least 10 alarms:
  VERIFIED (alarms carry information beyond the last hour) if at least 80% of cells pass and the pooled
  difference is at least 5 points; FAILED if at most 20% pass or the pooled difference is under 2 points;
  PARTIAL otherwise. Holm over the cells is reported too.

## Part 2. Strong-forecast check on 2026 data (1 January to 24 October 2026)

* Assets: BTC, ETH, Nifty 50 (primary); SOL, XRP, DOGE if time allows (secondary).
* Returns: crypto raw 1-minute log returns (load_binance, as E20); Nifty raw log returns of regular sessions
  (e07_nifty.load_clean). Nothing after 24 October 2026 is read. No FOMC or CPI time is used.
* Grid: crypto hourly (UTC hours), as E20. Nifty: the ends of session hours (session minutes 60, 120, ..., 360,
  and 375), and "a day" is one session (375 minutes).
* Forecast: E20's HAR-I (log RV of last day, mean of last 5 and 22 days, last hour, hour and weekday dummies),
  daily refitted expanding-window OLS from the start of the data (past only, E20's code copied and adapted).
* Alarm dummies U, D: a plug-in up (down) alarm in the hour ending at the grid time.
* Outcomes: log RV of the next 60 minutes and of the next day (1,440 minutes; Nifty 375 session minutes).
* Metric A (E20 style): on the 2026 grid only, regress log RV on log F(HAR-I) + U + D (asset intercepts).
  95% interval from a circular day-block bootstrap (7-day blocks, days shared across assets, 1,999 draws,
  seed 9,989,000 + case). Also the Driscoll-Kraay standard error of E20.
* Metric B (out of sample): HAR-I + U + D refitted the same expanding way. Difference in mean QLIKE on the 2026 grid,
  (HAR-I + alarms) minus HAR-I, with the same bootstrap interval. Negative = alarms help.
* Pass lines, per asset and horizon: "adds information" if the up coefficient has its 95% interval above 0 (down:
  below 0). "Improves the forecast" if the QLIKE difference interval lies below 0. Overall: the alarm improves the
  strong forecast (VERIFIED) if Metric B passes for at least 2 of the 3 primary assets at the next-hour horizon;
  FAILED if it passes for none at either horizon; PARTIAL otherwise.
