# Risk switch: sizing a position only when the volatility monitor says the regime changed

A working minimal prototype of alternative A6 in [`../ALTERNATIVES.md`](../ALTERNATIVES.md).
It asks one question, answers it on real data with every setting frozen in advance,
and reports what it finds, including the parts that do not favour the idea.
Numbers live in [`RESULTS.md`](RESULTS.md); every one of them comes from
`run_experiment.py` and is saved under `results/`.

## 1. Research problem

A risk manager who holds a long position in a volatile asset usually wants the
position's *risk*, not its size, to stay steady. The standard tool is
**volatility targeting**: estimate the asset's volatility, then hold
`target / estimate` units of it, so that when volatility doubles the position
halves. The estimate is refreshed on a schedule (every day in the academic
versions, every hour or minute in high-frequency practice), so the position is
nudged all the time, whether or not anything has changed. Every nudge costs money.

This project's monitor (`avmonitor.monitor.VolatilityMonitor`) does something
different: it watches one-minute returns and raises an alarm only when it has
collected enough evidence that the volatility regime has moved, with a proven
bound on false alarms. That suggests a **change-triggered risk switch**: keep
the position fixed until an alarm, then reset it to `target / (new regime's
robust volatility)` and hold again. Fewer, better-justified trades, if it works.

The research question is therefore:

> Does scaling a position by the monitor's regime estimate, changed only at
> alarms, give better **risk control** than standard volatility targeting,
> **after transaction costs**?

This is a risk-management study. It makes no claim about predicting returns,
and the Sharpe ratios are reported with wide bootstrap intervals for that reason.

## 2. Mathematical formulation

Returns are one-minute log returns `r_t`. An exposure rule is a sequence `w_t`
that is measurable with respect to information up to the close of minute `t`.
The position earns `w_t R_{t+1}` over minute `t + 1`, where `R` is the simple
close-to-close return (overnight gaps included for Nifty). Changing exposure
costs `kappa |w_t - w_{t-1}|` with `kappa` the one-way cost per unit of turnover.
Wealth compounds as

```
NAV_{t+1} = NAV_t (1 + w_t R_{t+1} - kappa |w_t - w_{t-1}|).
```

All targeting rules share one map from a per-minute variance forecast `v_t` to an
exposure:

```
sigma_t = sqrt( (B v_t + o_t) * S ),      w_t = min(w_max, target / sigma_t),
```

where `B` is the number of one-minute bars in a session (1440 for 24/7 crypto,
375 for Nifty), `S` the sessions per year (365 or 252) and `o_t` a forecast of
the squared overnight gap (zero for crypto; an EWMA over sessions for Nifty).
`target` is 40% a year for crypto and 15% for Nifty; `w_max = 2`.

**Classic volatility targeting (benchmark).** `v_t` is an exponentially weighted
moving average of squared minute returns with half-life `H` minutes,
`v_t = lambda v_{t-1} + (1 - lambda) r_t^2`, `lambda = 2^{-1/H}`. The exposure is
refreshed every minute, every hour or once a day (at the last bar of the session,
so it applies from the next session). This is the high-frequency version of the
rules in Moreira and Muir (2017) and Harvey et al. (2018).

**Robust rolling targeting (control).** `v_t = (k m_t)^2 s2_t`, where `m_t` is the
scaled median of `|u|` over a trailing window of deseasonalised returns
`u = r / s` (`s` is the predictable intraday profile from `avmonitor.data.intraday_scale`),
`s2_t` the trailing mean of `s^2`, and `k` a fixed conversion from robust scale to
standard deviation (minute returns are fat tailed, so the median-based scale
understates the standard deviation). Refreshed daily. It separates "robust
estimate" from "change-triggered timing".

**Change-triggered switch (the proposal).** Run the monitor on `u`. It
calibrates on 1440 bars, computes the robust scale `m` of that window, then
monitors with a Shiryaev-Roberts e-detector at level `A = B / alpha`, which
guarantees that the chance of a false alarm within a session is at most
`alpha` (5% or 1% here). When it alarms it recalibrates on the next 1440 bars.
The exposure changes only at these events:

* at the end of each calibration, `v = (k m)^2 s2`, so `w = min(w_max, target / sigma)`;
* optionally ("provisional" variant) at the alarm itself, with `m` replaced by the
  scaled median of the last 120 deseasonalised returns, which mostly belong to
  the new regime because detection delays are shorter than that; this bridges the
  calibration period, at the price of a second trade per alarm.

Between events `w_t` is held fixed. The same construction driven by the Gaussian
likelihood-ratio monitor (`GaussianLRMonitor`, no robustness to jumps) is the
third benchmark.

**Metrics.** With `RV_week` the annualised realised volatility of net minute
returns over a calendar week,

```
tracking RMS = sqrt( mean_weeks  log(RV_week / target)^2 ),
```

the delivered-risk error (0.10 means the weekly risk is typically 10% away from
the target). Also: realised annual volatility, maximum drawdown, worst session,
annualised Sharpe ratio (zero risk-free rate) with a 95% circular block-bootstrap
interval (blocks of 10 sessions, 2000 draws, seed 0), paired bootstrap intervals
for the Sharpe difference against daily EWMA targeting, turnover per year
(`sum |w_t - w_{t-1}|`), trades per year, and cost drag (`kappa` times turnover).

## 3. Why it matters

* Volatility targeting is one of the most widely used risk rules in asset
  management, and its benefit is documented across asset classes (Harvey et al.
  2018: higher Sharpe ratios for risk assets, fewer extreme losses everywhere).
  Its cost is constant rebalancing. A rule that trades only on statistically
  justified regime changes would be cheaper and easier to explain to a committee
  ("we cut the position because the monitor flagged a 30% rise in volatility,
  with at most a 5% chance per day of a false alarm").
* The monitor's alarms were shown to be followed by persistent volatility
  changes on Bitstamp BTC (E5 in `docs/experiment_report.md`). This study asks
  whether that statistical property translates into a usable risk rule, which is
  the natural next step and the question a practitioner will ask first.
* It is also a stress test of the monitor itself: how often it alarms on five
  coins and on an equity index over nine years, and whether jump-robustness
  matters for the trigger.

## 4. Comparison with existing approaches

| Approach | Estimate | When exposure changes | Robust to jumps | Guarantee on trades |
|---|---|---|---|---|
| Buy and hold | none | never | n/a | n/a |
| Volatility targeting, monthly realised variance (Moreira and Muir 2017) | realised variance of the previous month's daily returns | monthly | no | none |
| Volatility targeting, exponentially weighted (Harvey et al. 2018; RiskMetrics) | EWMA of squared daily returns | daily (here: minute, hour or day on minute data) | no | none |
| Robust rolling targeting (control here) | scaled median of absolute deseasonalised minute returns | daily | yes | none |
| Change-triggered switch (this study) | robust scale of the post-alarm calibration window | only at monitor events | yes (robust e-detector) or no (Gaussian LR) | false alarms at most 5% (or 1%) per session |

Moreira and Muir (2017) choose the scaling constant so that the managed
portfolio has the same unconditional volatility as buy and hold, which uses the
whole sample; later work (Cederburg et al. 2020, Liu et al. 2019) criticised that
as look-ahead. Here every constant is fixed in advance or estimated on the 2024
crypto window only, and the 2025 crypto and Nifty results are out of sample.

## 5. Implementation

```
alternatives/risk_switch/
├── README.md               this file
├── RESULTS.md              results and reading, every number from results/
├── strategies.py           exposure rules (EWMA, robust rolling, switch), data preparation
├── backtest.py             accounting, metrics, block bootstrap
├── run_experiment.py       the experiment (research questions in its docstring)
├── tests/test_risk_switch.py   synthetic tests, about 4 seconds, no data needed
└── results/                CSV tables, frozen_parameters.json, figures/
```

The code imports `avmonitor` (the monitor, the data loaders, the intraday
profile and the MAD constant) and nothing in `src/` imports it. `strategies.py`
builds predictable per-bar arrays once (`prepare`), then each rule is a pure
function of those arrays. `backtest.py` is vectorised numpy; the only Python
loop is the monitor's own per-minute update (about 1.5 s per coin-year).

Look-ahead is prevented by construction (exposure at `t` is applied to `t + 1`)
and checked by a test that perturbs the data after a bar and verifies that no
rule's exposure before that bar moves.

## 6. Experimental methodology

* **Data [REAL].** Binance spot one-minute bars for BTCUSDT, ETHUSDT, XRPUSDT,
  SOLUSDT and DOGEUSDT, 1 January 2024 to 31 December 2025 (no missing minutes,
  no duplicates, checks saved in `results/frozen_parameters.json`), and Nifty 50
  index one-minute bars, April 2017 to October 2026, regular sessions only.
* **Windows.** Tuning: crypto, 1 February to 31 December 2024 (January is
  warm-up for the 20-day profile and the first calibration). Out of sample:
  crypto, 1 January to 31 December 2025, with the frozen rules carried over
  from 2024 without a reset; Nifty, 1 June 2017 to 30 October 2026.
* **Grid on 2024 crypto.** EWMA half-lives 0.25, 1, 3, 10, 30 days at minute,
  hour and day refresh; robust windows 1 and 3 days; switch with bands 0.25,
  0.5, 0.75, false-alarm budgets 5% and 1% per session, with and without the
  provisional step; Gaussian LR switch with both budgets.
* **Selection rule (stated before the run).** For each family (EWMA minute,
  EWMA hour, EWMA day, robust rolling, robust switch): among variants whose mean
  cost drag at the base cost is at most 2% a year, the one with the lowest mean
  tracking RMS across the five coins. If no variant of a family is within
  budget, the lowest tracking RMS is taken anyway and the choice is flagged
  "over budget" in every table. The Gaussian switch inherits the robust
  switch's budget and provisional flag. The robust-to-standard-deviation
  conversion `k` is the median over coins of the median daily ratio in 2024.
* **Costs.** Base case 5 basis points per unit of turnover for crypto
  (Binance spot taker fees are about 10 bp before discounts, but perpetual
  futures and maker fills are cheaper; 5 bp is a middle value), 2.5 bp for Nifty
  futures (exchange and statutory charges plus half a tick). Sensitivity: 0, 2,
  5, 10, 20 bp for crypto and 0, 1, 2.5, 5, 10 bp for Nifty. No financing cost on
  leverage, no market impact.
* **Statistics.** Metrics per coin and for the equal-weight portfolio of the
  five coins (daily returns averaged). Sharpe intervals and paired differences
  by circular block bootstrap. Nifty is also broken down by calendar year.

Run everything with `python alternatives/risk_switch/run_experiment.py` (about
seven minutes on two cores; `--quick` runs one coin in two minutes). Tests:
`python -m pytest -q alternatives/risk_switch/tests`.

## 7. Limitations

* Two years of crypto and one equity index. The 2025 crypto year is a single
  out-of-sample year; Nifty gives nine years but a single asset.
* The cost model is one number per unit of turnover. Real costs rise with
  size and with volatility, which would hurt every active rule, most of all the
  ones that trade when volatility has just jumped.
* The switch's regime estimate is a one-day robust scale times a fixed constant
  `k`. On days with jumps the true standard deviation is larger than `k` times
  the robust scale, so the switch (and the robust rolling control) will sit
  above target on such days. That is a property of robust estimation, not a bug,
  and it is visible in the tracking numbers.
* The monitor is used exactly as shipped in `src/avmonitor` (plug-in baseline,
  no confidence-sequence baseline). Its alarm frequency, not its validity, is
  what drives the result.

## References

* Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482.
* Harvey, C. R., Hoyle, E., Korgaonkar, R., Rattray, S., Sargaison, M. and Van Hemert, O. (2018). The impact of volatility targeting. *Journal of Portfolio Management*, 45(1), 14 to 33. doi:10.3905/jpm.2018.45.1.014. (Checked by web search on 2 November 2026: https://scholars.duke.edu/publication/1370354.)
* Moreira, A. and Muir, T. (2017). Volatility-managed portfolios. *Journal of Finance*, 72(4), 1611 to 1644. doi:10.1111/jofi.12513. (Checked by web search on 2 November 2026: https://www.nber.org/papers/22208.)
* Cederburg, S., O'Doherty, M. S., Wang, F. and Yan, X. S. (2020). On the performance of volatility-managed portfolios. *Journal of Financial Economics*, 138(1), 95 to 117. (Cited for the look-ahead critique as summarised in the search results above; not independently verified here.)
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science*, 2(2), 229 to 260.
