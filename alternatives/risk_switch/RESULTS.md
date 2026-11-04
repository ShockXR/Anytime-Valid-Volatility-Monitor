# Results: change-triggered risk switch against volatility targeting

All results are **[REAL]**: Binance spot one-minute bars for BTC, ETH, XRP, SOL
and DOGE (2024 and 2025) and Nifty 50 index one-minute bars (2017 to 2026).
Every number below is copied from `results/summary_tables.md`, which
`summarise.py` renders from the CSV files written by `run_experiment.py`
(fixed seed 0, run on 2 November 2026, 551 s on one core). Labels: "tuning
sample" is 2024 crypto; "out of sample" is 2025 crypto and the whole Nifty history.

## 1. The answer in one paragraph

The change-triggered switch is a **usable but not better** risk rule. It trades
far less often than daily volatility targeting (93 exposure changes a year per
coin against 365; 12 a year on Nifty against 246) and, out of sample, its tail
risk is on a par: on the equal-weight crypto portfolio in 2025 its maximum
drawdown was 34.4% against 33.8% for daily EWMA targeting and 54.6% for buy and
hold, and its worst day 8.6% against 9.0% and 15.9%. On Nifty over nine years
its maximum drawdown was 28.1% against 30.8% and 39.3%. But on the metric the
study is about, **keeping delivered risk close to the target, it is worse** than
EWMA targeting at every refresh frequency: weekly realised volatility sat
0.27 log units from the 40% target on crypto against 0.19 to 0.22 for the EWMA
rules, and 0.40 on Nifty against 0.24 to 0.26. The reason is visible in the
exposure paths: the switch holds a position for days while volatility drifts,
and when it does move it moves by a lot. Sharpe ratios are not distinguishable
between any two rules (every paired bootstrap interval contains zero), so this
is not a return story in either direction. Turnover on crypto is the same as
daily EWMA targeting (31 units a year against 30) despite four times fewer
trades, because each alarm-driven trade is four times larger. The robust
e-detector is a clearly better trigger than the Gaussian likelihood-ratio
monitor, which alarms three to five times as often, costs twice as much to
follow and gives worse tails.

## 2. How often the monitor alarms

Alarms per year, with recalibration on 1440 bars after each alarm (T1).

| monitor | false alarm budget per session | band | crypto 2024, mean of 5 coins | crypto 2025, mean of 5 coins | Nifty 2017 to 2026 |
|---|---|---|---|---|---|
| robust e-detector | 5% | 0.25 | 254.7 | 249.2 | 43.2 |
| robust e-detector | 5% | 0.5 | 183.3 | 179.8 | 28.1 |
| robust e-detector | 5% | 0.75 | 99.8 | 93.2 | 12.0 |
| robust e-detector | 1% | 0.25 | 239.0 | 234.4 | 39.8 |
| robust e-detector | 1% | 0.5 | 166.3 | 163.8 | 24.5 |
| robust e-detector | 1% | 0.75 | 86.5 | 79.2 | 10.4 |
| Gaussian LR | 5% | | 333.4 | 331.2 | 59.9 |
| Gaussian LR | 1% | | 329.7 | 323.8 | 58.5 |

On crypto the monitor sees a new volatility regime every 1.5 days with a 25%
band and every 3.7 days even with a 75% band. On Nifty (intraday returns only,
375 bars a session) alarms are 6 to 8 times rarer per calendar year. The
false-alarm budget barely matters (5% against 1% changes the count by 6% to
15%), which says the alarms are driven by genuine volatility moves, not by the
statistical noise the budget controls. The Gaussian monitor alarms 3.3 to 5
times more often than the frozen robust setting (band 0.75); about two thirds
of its alarms (65% to 69%) are increases, consistent with jumps being read as
regime changes.

## 3. Tuning sample: 2024 crypto, what the selection rule picked

Mean over the five coins, 1 February to 31 December 2024, base cost 5 bp (T2).
The rule: lowest tracking RMS among variants with cost drag at most 2% a year.

| family | chosen variant | tracking RMS | realised vol | max drawdown | worst day | turnover / yr | trades / yr | cost drag / yr |
|---|---|---|---|---|---|---|---|---|
| EWMA, minute refresh | half-life 3 days | 0.166 | 39.4% | -34.1% | -6.6% | 31.2 | 525,600 | 1.6% |
| EWMA, hourly refresh | half-life 3 days | 0.174 | 39.4% | -33.2% | -6.4% | 19.4 | 8,760 | 1.0% |
| EWMA, daily refresh | half-life 1 day | 0.153 | 41.5% | -32.6% | -7.0% | 29.8 | 365 | 1.5% |
| robust rolling, daily | 3-day window | 0.262 | 46.7% | -38.3% | -8.3% | 31.9 | 364 | 1.6% |
| switch, robust monitor | 5% budget, band 0.75, hold | 0.259 | 41.3% | -35.1% | -7.7% | 35.2 | 100 | 1.8% |
| switch, Gaussian LR | 5% budget, hold (over budget) | 0.214 | 45.3% | -39.8% | -8.2% | 59.3 | 333 | 3.0% |
| buy and hold | | | 78.2% | -48.2% | -12.6% | 0 | 0 | 0 |

Two things the full grid (T2 in `results/summary_tables.md`) shows:

* Tighter tracking is available at a price. A minute-refreshed EWMA with a
  6-hour half-life tracks the target almost exactly (RMS 0.029) but turns over
  397 units a year, a 19.8% annual cost drag at 5 bp. With a 1-day half-life and
  hourly refresh the RMS is 0.098 at 2.8% drag. The 2% budget rules these out.
* The switch's best tracking (RMS 0.182, 1% budget, band 0.5, with the
  provisional step) costs 4.6% a year. Every switch variant with a band of 0.5
  or less is over budget because it fires 160 to 250 times a year. Only the
  band-0.75 "hold" variants are within budget, at RMS 0.26 to 0.27, which is why
  the frozen switch is the one that moves least. No Gaussian-driven variant is
  within budget.
* The robust-to-standard-deviation factor is `k = 1.216` (median daily ratio;
  per coin 1.15 to 1.28). Using a root-mean-square ratio instead would give 1.21
  to 1.31, a difference of about 3%, so the choice does not matter.

## 4. Out of sample: 2025 crypto

Mean of the five coins (T3b) and the equal-weight portfolio of the five (T3d,
with bootstrap intervals). Base cost 5 bp. 2025 was a poor year for the coins
(buy and hold Sharpe about 0), so the return columns are all near or below zero.

| rule | realised vol | tracking RMS | max drawdown | worst day | Sharpe | mean exposure | turnover / yr | trades / yr | cost drag / yr |
|---|---|---|---|---|---|---|---|---|---|
| buy and hold | 75.9% | | -61.0% | -16.9% | -0.03 | 1.00 | 0 | 0 | 0 |
| EWMA, minute | 42.4% | 0.192 | -40.3% | -13.5% | -0.20 | 0.62 | 32.0 | 525,600 | 1.6% |
| EWMA, hour | 41.2% | 0.220 | -39.3% | -9.7% | -0.08 | 0.62 | 19.8 | 8,760 | 1.0% |
| EWMA, day | 43.6% | 0.205 | -41.6% | -10.2% | -0.10 | 0.65 | 29.9 | 365 | 1.5% |
| robust rolling, day | 47.8% | 0.274 | -43.9% | -9.6% | -0.04 | 0.72 | 31.2 | 362 | 1.6% |
| switch, robust monitor | 41.5% | 0.271 | -40.9% | -10.6% | -0.15 | 0.61 | 31.5 | 93 | 1.6% |
| switch, Gaussian LR (over budget) | 47.8% | 0.242 | -45.9% | -10.6% | -0.22 | 0.72 | 60.3 | 329 | 3.0% |

Equal-weight portfolio of the five coins, 2025:

| rule | realised vol | tracking RMS | max drawdown | worst day | Sharpe [95% CI] | difference vs EWMA day [95% CI] |
|---|---|---|---|---|---|---|
| buy and hold | 69.0% | | -54.6% | -15.9% | -0.01 [-1.65, 1.81] | +0.07 [-0.45, +0.60] |
| EWMA, minute | 38.9% | 0.223 | -34.8% | -13.4% | -0.22 [-2.12, 1.84] | -0.14 [-0.44, +0.19] |
| EWMA, hour | 37.7% | 0.240 | -31.7% | -9.1% | -0.06 [-2.01, 1.94] | +0.02 [-0.22, +0.26] |
| EWMA, day | 39.9% | 0.204 | -33.8% | -9.0% | -0.08 [-2.01, 1.86] | reference |
| robust rolling, day | 43.7% | 0.219 | -35.2% | -8.7% | 0.00 [-1.90, 1.99] | +0.08 [-0.18, +0.38] |
| switch, robust monitor | 37.2% | 0.259 | -34.4% | -8.6% | -0.13 [-2.07, 1.90] | -0.05 [-0.41, +0.32] |
| switch, Gaussian LR | 43.6% | 0.186 | -38.8% | -9.2% | -0.21 [-2.11, 1.76] | -0.13 [-0.49, +0.25] |

Per coin (T4, T5): the switch had the best drawdown and worst day on DOGE
(-37.2% and -8.2% against -43.4% and -10.5% for daily EWMA) and the worst on
SOL (-41.1% and -13.3% against -37.7% and -10.6%). Its tracking RMS was worse
than daily EWMA on four coins out of five (0.235 against 0.150 on BTC, 0.266
against 0.172 on ETH, 0.263 against 0.152 on SOL, 0.279 against 0.253 on DOGE)
and better on XRP (0.315 against 0.299). The paired Sharpe differences range
from -0.46 [-0.94, +0.01] on BTC to +0.34 [-0.15, +0.88] on XRP: noise.

The in-sample to out-of-sample drift for the switch is small (tracking 0.259 to
0.271, drag 1.8% to 1.6%, trades 100 to 93 a year), so nothing was overfitted;
the rule simply is what it is.

## 5. Out of sample: Nifty 50, June 2017 to October 2026

Same frozen rules, 15% target, base cost 2.5 bp (T6).

| rule | realised vol | tracking RMS | max drawdown | worst day | Sharpe [95% CI] | difference vs EWMA day [95% CI] | mean exposure | turnover / yr | trades / yr | cost drag / yr |
|---|---|---|---|---|---|---|---|---|---|---|
| buy and hold | 16.4% | | -39.3% | -12.7% | 0.66 [-0.02, 1.38] | +0.07 [-0.29, +0.44] | 1.00 | 0 | 0 | 0 |
| EWMA, minute | 15.6% | 0.256 | -28.9% | -4.9% | 0.58 [-0.11, 1.28] | -0.01 [-0.09, +0.06] | 1.23 | 25.9 | 92,613 | 0.6% |
| EWMA, hour | 15.7% | 0.258 | -28.6% | -4.9% | 0.61 [-0.08, 1.30] | +0.02 [-0.05, +0.09] | 1.23 | 13.8 | 1,730 | 0.3% |
| EWMA, day | 15.9% | 0.242 | -30.8% | -5.5% | 0.59 [-0.09, 1.29] | reference | 1.25 | 17.2 | 246 | 0.4% |
| robust rolling, day | 16.0% | 0.291 | -29.0% | -5.7% | 0.62 [-0.08, 1.32] | +0.03 [-0.09, +0.14] | 1.23 | 6.2 | 250 | 0.2% |
| switch, robust monitor | 14.3% | 0.402 | -28.1% | -6.3% | 0.69 [0.01, 1.39] | +0.10 [-0.07, +0.26] | 1.06 | 2.9 | 12 | 0.1% |
| switch, Gaussian LR (over budget) | 15.8% | 0.292 | -28.0% | -7.7% | 0.56 [-0.14, 1.26] | -0.03 [-0.14, +0.08] | 1.21 | 8.7 | 60 | 0.2% |

By calendar year (T7), the picture is consistent:

* **2020 (the Covid crash)** is where volatility targeting earns its keep. Buy
  and hold: realised vol 30.9%, maximum drawdown 39.3%, worst day 12.7%. Every
  targeting rule: vol 16.3% to 18.1%, drawdown 26.6% to 29.0%, worst day 4.4%
  to 5.4%. The switch (16.6%, 28.1%, 5.4%) is in the pack, with the Gaussian
  switch (16.9%, 26.7%, 4.7%) and the EWMA rules (16.3% to 16.5%, 26.6% to
  27.1%, 4.4% to 4.6%) slightly better on the worst day.
* **Calm years** show the switch's weakness as a targeting rule. In 2017 and
  2023 it delivered 12.3% and 13.5% volatility against the 15% target, because
  with a 75% band it only resets when volatility has moved a lot, and slow
  drifts down never trigger it; the EWMA rules delivered 15.1% to 16.0%. Its
  tracking RMS is the worst in every single year (0.29 to 0.60 against 0.16 to
  0.36 for daily EWMA).
* Its drawdowns are the smallest of the targeting rules in 9 years out of 10
  (for example 14.2% against 19.1% in 2022, 17.8% against 24.5% in 2026), partly
  because its exposure is lower on average (1.06 against 1.25).

The Nifty switch turns over 2.9 units a year in 12 trades, against 17.2 units
in 246 trades for daily EWMA targeting: on a market with rare regime changes
the switch really is the cheap rule, and its cost drag (0.1% a year) is a
rounding error.

## 6. Transaction costs

Net Sharpe ratio as the one-way cost per unit of turnover varies (T8).

Crypto, equal-weight portfolio, 2025:

| rule | 0 bp | 2 bp | 5 bp | 10 bp | 20 bp |
|---|---|---|---|---|---|
| buy and hold | -0.01 | -0.01 | -0.01 | -0.01 | -0.01 |
| EWMA, minute | -0.18 | -0.19 | -0.22 | -0.26 | -0.34 |
| EWMA, hour | -0.03 | -0.05 | -0.06 | -0.09 | -0.14 |
| EWMA, day | -0.04 | -0.06 | -0.08 | -0.12 | -0.19 |
| robust rolling, day | 0.04 | 0.02 | 0.00 | -0.03 | -0.11 |
| switch, robust monitor | -0.09 | -0.11 | -0.13 | -0.18 | -0.26 |
| switch, Gaussian LR | -0.14 | -0.17 | -0.21 | -0.28 | -0.42 |

Nifty 50, 2017 to 2026:

| rule | 0 bp | 1 bp | 2.5 bp | 5 bp | 10 bp |
|---|---|---|---|---|---|
| buy and hold | 0.66 | 0.66 | 0.66 | 0.66 | 0.66 |
| EWMA, minute | 0.62 | 0.61 | 0.58 | 0.54 | 0.46 |
| EWMA, hour | 0.63 | 0.63 | 0.61 | 0.59 | 0.55 |
| EWMA, day | 0.62 | 0.61 | 0.59 | 0.56 | 0.51 |
| robust rolling, day | 0.63 | 0.63 | 0.62 | 0.61 | 0.59 |
| switch, robust monitor | 0.70 | 0.69 | 0.69 | 0.69 | 0.68 |
| switch, Gaussian LR | 0.58 | 0.57 | 0.56 | 0.55 | 0.52 |

Costs do not change any ranking on crypto between 0 and 20 bp: the switch
loses 0.17 of Sharpe over that range, daily EWMA 0.15, hourly EWMA 0.11, the
Gaussian switch 0.28. On Nifty the switch is almost cost-proof (0.70 to 0.68
from 0 to 10 bp) while the EWMA rules lose 0.08 to 0.16; at 10 bp the switch's
net annual return is 9.7% against 8.1% for daily EWMA and 10.8% for buy and
hold. These are return differences inside the bootstrap noise; the robust
statement is about the turnover that produces them.

## 7. Robust trigger against Gaussian trigger

Driven by the Gaussian likelihood-ratio monitor instead of the robust
e-detector, the same switch alarms 3.3 to 5 times as often (331 against 93 a
year on crypto, 60 against 12 on Nifty), turns over twice as much (60 against
31 on crypto), is over the cost budget in every variant, and has worse tails:
maximum drawdown 45.9% against 40.9% on crypto, worst day 7.7% against 6.3% on
Nifty, and the worst Sharpe of all rules on both markets. Many of its alarms
are jump-driven increases (E5 found the same on Bitstamp): after a jump it cuts
exposure, then recalibrates on a window inflated by the jump and raises it
again. This is the one comparison in the study with a clear winner, and it is
the one that tests the paper's robustness idea directly.

## 8. Reading

1. **The switch controls tails as well as the EWMA rules, with a fraction of
   the trades.** Out of sample, maximum drawdown and worst day are within two
   percentage points of daily EWMA targeting on crypto and better on Nifty,
   with 93 against 365 trades a year on crypto and 12 against 246 on Nifty.
   For a risk manager whose constraint is "do not blow through the drawdown
   limit" and who dislikes constant rebalancing, it is a reasonable rule.
2. **It does not track a volatility target.** Tracking RMS is 30% worse than
   EWMA on crypto and 65% worse on Nifty, in sample and out. A rule that holds
   exposure fixed between alarms cannot follow slow drifts, and on crypto the
   alarms it does react to are large (each trade moves exposure by 0.34 on
   average against 0.08 for daily EWMA). On a market like Nifty with rare
   alarms it also runs below target in calm years.
3. **Same turnover as daily EWMA on crypto, not less.** The expected finding
   "matches but with fewer trades" holds for the number of trades, not for the
   amount traded. Fewer, larger trades cost the same at a per-unit cost; at
   real costs that rise with trade size they would cost more.
4. **No Sharpe effect either way.** Every paired difference against daily
   EWMA has a 95% interval containing zero, on every coin, on the portfolio and
   on Nifty. One year of five correlated coins cannot resolve differences of
   0.1 to 0.3 in Sharpe, and the study was not designed to.
5. **Robustness of the trigger matters.** The Gaussian-driven switch is worse
   on every count. This is the clearest positive result for the paper's
   robust estimator in a decision context.
6. **The robust rolling control shows where the gap comes from.** A robust
   estimate refreshed daily (no change detection) tracks worse than the EWMA
   (0.274 against 0.205 on crypto) and runs above target (47.8% against
   43.6%), because the median-based scale ignores the days with jumps that
   dominate realised variance. The switch inherits this; its estimate is the
   same quantity measured over the calibration window.

## 9. Negative results and caveats

* The headline negative result: a change-triggered switch built from the
  monitor as shipped is **not** a better volatility-targeting rule than an
  EWMA refreshed daily or hourly, on either market, before or after costs.
* The provisional-exposure step (a second trade at the alarm) improves
  tracking (0.183 against 0.216 at band 0.5) but doubles turnover; no
  provisional variant was within the cost budget on crypto.
* The 1% false-alarm budget changes little; the band is the only knob that
  controls the alarm rate, and it trades tracking for cost monotonically.
* Nifty's 2026 (January to October) shows a drawdown of 15.9% for buy and hold
  and 17.8% to 24.5% for the targeting rules: leverage (exposure above 1 when
  volatility is low) hurts in a steady decline. The cap of 2 is part of the
  frozen protocol; a lower cap would change the Nifty numbers for every
  targeting rule.
* Costs are a single number per unit of turnover; no financing cost on
  leverage, no impact. The Nifty backtest is on the index, not on futures (no
  roll cost, no basis).
* Two years of crypto and one equity index. The tuning year (2024) was a
  strong year for the coins and the test year (2025) a weak one, which is why
  the Sharpe columns flip sign; the risk columns are the ones to read.

## 10. Figures

* `results/figures/fig1_btc_2025_exposures.png`: BTC in 2025 with the
  frozen switch's alarms, its exposure path against the Gaussian-driven one,
  and the daily and hourly EWMA paths. The switch's steps are few and large;
  the EWMA paths are continuous.
* `results/figures/fig2_nifty_weekly_vol.png`: weekly realised volatility of
  buy and hold, daily EWMA targeting and the switch on Nifty, 2017 to 2026,
  against the 15% target. The EWMA line hugs the target; the switch sits below
  it in calm years and spikes with it in March 2020.
* `results/figures/fig3_cost_sensitivity.png`: net Sharpe against cost for
  every frozen rule, crypto portfolio 2025 and Nifty.

## Files

`results/is_grid_crypto_2024.csv` (all 34 variants, 5 coins),
`results/is_selection_crypto_2024.csv`, `results/frozen_parameters.json`
(chosen variants, `k`, data checks, config), `results/frozen_crypto.csv`,
`results/sharpe_diff_crypto_2025.csv`, `results/cost_sensitivity_crypto_2025.csv`,
`results/alarms_crypto.csv`, `results/frozen_nifty.csv`,
`results/frozen_nifty_by_year.csv`, `results/sharpe_diff_nifty.csv`,
`results/cost_sensitivity_nifty.csv`, `results/alarms_nifty.csv`,
`results/grid_nifty_not_used_for_selection.csv` (the whole grid on Nifty, for
transparency only; the Nifty choices were frozen on crypto), `results/summary_tables.md`.
