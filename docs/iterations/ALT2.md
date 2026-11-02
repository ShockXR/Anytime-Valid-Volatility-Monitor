# ALT2 iteration notes: change-triggered risk switch (alternatives/risk_switch)

One entry per substantial iteration, in order. All dates are 2 November 2026 (UTC).

## Iteration 1: design and first pipeline

**OBJECTIVE.** Turn A6 of `alternatives/ALTERNATIVES.md` into a working prototype that answers one question: does sizing a long position by the monitor's regime estimate, changed only at alarms, control risk better than EWMA volatility targeting after costs?

**PROBLEM IDENTIFIED.** The monitor's regime estimate (the robust scale of the calibration window) is only known 1440 bars after an alarm, so a pure "set exposure at the alarm" rule has nothing to set it to. Also the monitor works on deseasonalised returns `u = r / s`, whose scale has to be mapped back to raw annualised volatility, and the median-based scale understates the standard deviation of fat-tailed minute returns.

**HYPOTHESIS.** Two switch variants cover the gap: "hold" (change exposure only when the recalibration completes) and "provisional" (also set a bridge exposure at the alarm from the last 120 deseasonalised returns, which are mostly post-change because detection delays are shorter). A single conversion factor `k` = median daily ratio of realised standard deviation to robust scale, estimated on 2024 crypto and frozen, maps robust scale to standard deviation.

**CHANGES IMPLEMENTED.** `strategies.py` (`prepare`, EWMA / robust rolling / switch rules, `monitor_events`, `robust_to_sd_factor`), `backtest.py` (accounting with cost at the trade bar, metrics, circular block bootstrap), synthetic tests.

**EXPERIMENTS PERFORMED.** Smoke run on BTC 2024 (not reported as a result).

**RESULTS.** The first run crashed: `ZeroDivisionError` inside `VolatilityMonitor._calibrate`. Traced to the weekend of 22 to 23 June 2024, when more than half of the minutes in a 1440-bar window had exactly zero return on Binance BTCUSDT (price pinned at 64277.47 for many minutes with non-zero trade counts), so the median of `|u|` was 0, the truncation level 0 and the decrease side's null value 0. The main package already has `scale_floor` for this; it had never been needed on Bitstamp data.

**IMPROVEMENTS.** A numerical floor on the robust scale (`SCALE_FLOOR = 1e-6` per minute, 0.07% a year for a 24/7 market) passed to the monitor, the switch and the robust rolling rule. My first attempt set the floor at 0.2 under the wrong belief that `u` has unit scale; `u` is in return units (the profile is a unitless multiplier with median 1), so 0.2 froze every exposure at 0.002. Caught by looking at the mean exposure column.

**REMAINING ISSUES.** Which EWMA half-life is the fair benchmark, and how to pick it without tuning on the out-of-sample years.

**NEXT ITERATION.** A pre-stated selection rule on 2024 crypto only.

## Iteration 2: selection rule, degenerate band and the cost of following every alarm

**OBJECTIVE.** Fix the grid and the selection rule before any out-of-sample number is produced.

**PROBLEM IDENTIFIED.** (a) On BTC 2024 the minute-refreshed EWMA with a 6-hour half-life tracks the 40% target almost perfectly (tracking RMS 0.025) but turns over 529 units a year, a 26% annual cost drag at 5 bp; choosing by tracking alone would pick an unusable rule. (b) Band 1.0 is degenerate: the decrease side's null becomes "variance at least 0 times baseline", the monitor raises only a handful of increase alarms a year and the switch freezes at whatever its first calibration said (mean exposure 0.33 on BTC, tracking RMS 1.04). (c) The robust switch at bands 0.25 and 0.5 raises 190 to 260 alarms a year on BTC, moves exposure by about 0.35 at each, and so turns over 90 to 160 units a year (4.4% to 8% cost drag at 5 bp), more than daily EWMA targeting (16 to 43). A selection rule with a cost budget then has no eligible switch variant and the fallback "cheapest variant" returned the degenerate band 1.0.

**HYPOTHESIS.** A rule "lowest tracking RMS among variants with cost drag at most 2% a year, else lowest tracking RMS flagged as over budget" is simple, stated in advance, and keeps the switch's representative meaningful. Wider bands below 1 and a 1% false-alarm budget may reduce the alarm count enough to bring a switch variant within budget. A 3-day robust window is the fair analogue of the EWMA's smoothing (the 1-day rolling median's day-to-day log change has standard deviation 0.43 against 0.18 for the 1-day EWMA, mostly genuine weekend swings).

**CHANGES IMPLEMENTED.** Grid: bands 0.25, 0.5, 0.75; false-alarm budgets 5% and 1% per session; robust windows 1 and 3 days. `select()` with the over-budget flag carried into every output table. `run_experiment.py` saves the whole in-sample grid, the frozen set in and out of sample, paired bootstrap differences, cost sensitivity, alarm counts, a Nifty by-year table and three figures.

**EXPERIMENTS PERFORMED.** Quick runs (one coin, short Nifty sample) to debug; the full run follows in iteration 3.

**RESULTS.** See iteration 3.

**IMPROVEMENTS.** Metrics rewritten with `np.bincount` on precomputed session and week codes (about 20 times faster than the pandas version), which keeps the full run under ten minutes.

**REMAINING ISSUES.** None blocking.

**NEXT ITERATION.** Full run and write-up.
