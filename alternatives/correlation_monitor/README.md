# Correlation monitor: anytime-valid detection of dependence breaks between two assets

Alternative direction A7 of `alternatives/ALTERNATIVES.md`. This folder is a self-contained prototype. It imports `avmonitor` (the main package) but nothing in `src/` imports it. Every number quoted in this folder lives in `RESULTS.md` and comes from `run_experiment.py`, labelled [SYNTHETIC] or [REAL].

## 1. Research problem

The main package watches the volatility of one asset: it clips each return at a robust threshold (Bapat's `c x MAD` rule), which makes the squared return bounded, and runs a betting e-detector on the bounded statistic. The guarantee ("at most a 5% chance of a false alarm per day, however often you look") holds without any model of the returns. Bapat (2026) names a multivariate extension as future work. The question here is:

> Can the same recipe (robust scale, truncation, betting e-detector) detect a change in the **dependence between two assets** online, with the same finite-sample false alarm guarantee, without mistaking a change in the volatility of one asset for a change in their correlation, and robustly to jumps?

The last two clauses are the hard part. Correlation is estimated from products of returns, and products of returns move for three unrelated reasons: the correlation changed, the volatility of one asset changed, or one asset jumped. A monitor for *correlation* must respond to the first only.

## 2. Mathematical formulation

### 2.1 Set-up

Two return series `r1_t, r2_t` (one-minute or five-minute log returns, after the same intraday deseasonalisation as the main package). `F_t` is the information up to time t. A calibration window gives, for each asset, the robust scale `s_k = 1.4826 x median|r_k|` (scaled MAD about zero) and the truncation level `tau_k = c s_k` with `c = 3`, exactly as in `avmonitor.monitor`.

### 2.2 Two bounded statistics

* **Truncated product** (the natural "truncated covariance" candidate)

  `P_t = clip(r1_t / s1, -c, c) x clip(r2_t / s2, -c, c)`, with `-c^2 <= P_t <= c^2`.

  For Gaussian pairs its mean is close to the Pearson correlation `rho` (0.299 for `rho = 0.3` at `c = 3`).

* **Sign concordance** (the statistic we recommend)

  `S_t = sign(r1_t) x sign(r2_t)`, with `S_t in {-1, 0, +1}` (0 when either return is exactly zero).

  Its mean is `P(same sign) - P(opposite signs)`, which is Blomqvist's beta (Blomqvist 1950), also called the median correlation, when the returns have median zero. The sign statistic does not need the scales at all: it is the same for `r1` and for `a r1` with any `a > 0`.

**Proposition 1 (link to correlation).** If `(r1_t, r2_t)` given `F_{t-1}` is elliptical with centre zero and correlation parameter `rho` (for example bivariate normal or bivariate Student t, both with any scales), then

`E[S_t | F_{t-1}] = (2 / pi) arcsin(rho)`.

*Proof.* An elliptical vector with centre zero can be written as `R A U`, where `U` is uniform on the unit circle, `R > 0` is independent of `U`, and `A A' = Sigma`. The signs of the two coordinates do not depend on `R`, so the probability of each sign pattern is the same as in the Gaussian case, where Sheppard's formula gives `P(both positive) = 1/4 + arcsin(rho) / (2 pi)`. Then `E[S] = 2 [P(++) + P(--)] - 1 = (2 / pi) arcsin(rho)`. (For Kendall's tau the same identity for elliptical distributions is in Lindskog, McNeil and Nyfeler 2003.) The map `rho -> (2 / pi) arcsin(rho)` is strictly increasing, so "correlation went up" and "mean of `S` went up" are the same statement for elliptical pairs.

Outside the elliptical family `E[S]` is still a concordance measure, so "the dependence changed" is what an alarm means in general. The code reports `rho_hat = sin(pi m / 2)` as the implied correlation.

### 2.3 Null hypotheses and bets

Let `Y_t` be the chosen statistic with almost sure bounds `lo <= Y_t <= hi` (`[-1, 1]` for the sign, `[-c^2, c^2]` for the product), and `m0` the calibration value of its mean. The nulls are **conditional moment statements**, as in `docs/math.md` section 5.1:

* increase: `H0_up: E[Y_t | F_{t-1}] <= m0 + b` for all t,
* decrease: `H0_down: E[Y_t | F_{t-1}] >= m0 - b` for all t,

with an additive null band `b >= 0` (the statistic lives on a fixed scale, so the band is additive; `b = 0.1` on the sign scale is about 0.15 in correlation near zero). Nothing else is assumed: no distribution, no independence, no moments of the returns, and any jump process is allowed as long as it keeps the conditional mean on the null side.

For a predictable bet `lam_t in [0, 1 / (m_up - lo))`, with `m_up = m0 + b`, define `f_t = 1 + lam_t (Y_t - m_up)`. Then `f_t >= 1 - lam_t (m_up - lo) > 0` and `E[f_t | F_{t-1}] <= 1` under `H0_up`. For a decrease take `f_t = 1 + lam_t (m_dn - Y_t)` with `lam_t in [0, 1 / (hi - m_dn))`. Products of factors started at any time j are therefore test supermartingales, and the Shiryaev-Roberts statistic

`M_t = f_t (M_{t-1} + 1)`, `M_0 = 0`, alarm at `inf{t : M_t >= A}`,

is an e-detector (Shin, Ramdas and Rinaldo 2024).

**Theorem (validity).** Under the null, `E[alarm time] >= A` and `P(alarm within m steps) <= m / A`.

*Proof.* Theorems A and B of `docs/math.md` section 5.4 use only two facts: `f_t >= 0` and `E[f_t | F_{t-1}] <= 1`. Both were shown above. The two-sided monitor runs each side at `2A`, so `(M_up + M_down) / 2` is again an e-detector at level A and the first alarm of either side inherits both bounds. An equal-weight average over a grid of bets (fractions 0.05 to 0.75 of the largest admissible bet, as in `corrmon.DEFAULT_FRACTIONS`) is also an e-detector, by the same linearity argument as in the main package.

The sign statistic is bounded without any truncation, so the decrease side is valid by itself. For the product the bound comes from Bapat's truncation, as for the volatility monitor (`docs/math.md`, 5.3).

### 2.4 Invariance to volatility changes (what separates a correlation break from a volatility break)

**Proposition 2.** Suppose `r_k,t = sigma_k,t eps_k,t` for k = 1, 2, where `sigma_k,t > 0` is `F_{t-1}`-measurable (any predictable volatility process, including one that jumps to a new regime) and the conditional law of `(eps_1t, eps_2t)` given `F_{t-1}` does not change. Then the conditional law of `S_t` does not change either, so neither `H0_up` nor `H0_down` for the sign statistic is affected by any change in `sigma_1` or `sigma_2`.

*Proof.* `sign(sigma eps) = sign(eps)` for `sigma > 0`. The conditional law of `(sign eps_1t, sign eps_2t)` is fixed by assumption.

The product statistic has no such property: `clip(sigma_1 eps_1 / s1, -c, c)` grows with `sigma_1` for every `eps_1`, so a volatility increase in one asset raises `E[P_t]` (doubling the volatility of one asset roughly doubles the mean, see the tests) and the product monitor raises a "correlation up" alarm. The same holds for a Gaussian likelihood ratio and for any Pearson-type statistic. This is the central design point of the prototype, and the experiments in `RESULTS.md` Q1 measure it.

A one-line consequence for practice: the sign monitor needs no scale estimate, no truncation level and no intraday deseasonalisation (all of these are positive multiplicative factors). Its only calibration quantity is the mean `m0`.

### 2.5 Robustness to jumps

A jump in one asset changes one `S_t` by at most 2 and the log-wealth of a bet by at most `log(1 + lam (1 - m_up))` or `log(1 + lam (hi - m_dn))`: bounded influence, which is the betting version of Theorem 3 of Bapat (2026). Independent jumps, and simultaneous jumps with independent signs, leave `E[S]` unchanged (checked in the tests). Simultaneous jumps with the *same* sign are a genuine increase in dependence, and the monitor is meant to flag them.

### 2.6 Where the truncation level enters the speed

The betting speed is governed by the admissible bet. For the sign statistic the largest increase bet is `1 / (1 + m0)`, about 0.84 at `rho = 0.3`. For the product at `c = 3` it is `1 / (9 + m0)`, about 0.11, far below the growth-optimal bet `(m1 - m0) / (var + (m1 - m0)^2)` of about 0.25 for a shift from 0.3 to 0.6. So the bound that makes the product valid also throttles it. The theory table (Q0 in `RESULTS.md`) quantifies this: on clean Gaussian data the sign monitor is predicted to be about as fast as the product with `c = 3` for moderate changes and faster for large ones, and the product with `c = 2` is fastest. This is a concrete reason to prefer the scale-free statistic: it gives up little speed and gains invariance.

## 3. Why it matters

* **Portfolios.** Diversification relies on correlations staying moderate. In sell-offs they rise (Longin and Solnik 2001 for equity markets), which is when the diversification is needed. A monitor that flags "dependence has increased" early, with a controlled false alarm rate, is a risk control input: it can trigger a hedge review, a reduction of gross exposure or a switch of the risk model.
* **Separating the two breaks.** Risk systems typically react to volatility. A correlation break without a volatility break (or the reverse) needs a different response (position sizing against concentration, not against variance). A statistic that cannot confuse the two is useful for that reason alone.
* **Validity under crypto conditions.** Crypto pairs have fat tails, co-jumps and frequent volatility regime changes, so Gaussian sequential tests are not credible there. The sign monitor's guarantee is model-free and finite-sample.
* **For the paper's author.** This is the sequential, bounded-influence version of the multivariate extension named as future work, and the one place in this repository where a scale-free statistic makes Bapat's truncation unnecessary for validity (and shows what the truncation costs in speed when it is used).

## 4. Comparison with existing approaches

| Approach | Online? | Finite-sample false alarm control? | Robust to jumps and fat tails? | Confuses a volatility change with a correlation change? |
|---|---|---|---|---|
| Retrospective fluctuation test for a constant correlation (Wied, Kraemer and Dehling 2012) | No (fixed sample, asymptotic Brownian bridge critical values) | No (asymptotic) | Needs finite fourth moments | Partly (it standardises by full-sample variances) |
| Sequential monitoring with boundary functions (Chu, Stinchcombe and White 1996) | Yes | Asymptotic only | Model dependent | Depends on the statistic |
| Gaussian likelihood ratio e-detector for the correlation (benchmark here) | Yes | Yes, but only if pairs are exactly bivariate normal with the calibration scales | No | Yes |
| Page CUSUM on the raw product of standardised returns (naive benchmark here) | Yes | No (threshold tuned by simulation under an assumed model) | No | Yes |
| Betting e-detector on the truncated product (this folder) | Yes | Yes, model-free | Yes (bounded) | Yes |
| **Betting e-detector on the sign concordance (this folder)** | Yes | Yes, model-free | Yes (bounded, no scale needed) | **No (Proposition 2)** |

Related estimation literature that this prototype does not replace: realized covariance and the Hayashi and Yoshida (2005) estimator deal with asynchronous observation for *estimating* covariation; the Epps (1979) effect (correlations estimated from very short intervals are biased toward zero because the two assets do not trade at the same instants) is handled here by aggregating to 5-minute bars before monitoring, and its size is measured in `RESULTS.md` Q3a. Sign-based and rank-based correlation measures are classical (Blomqvist 1950); what is new here is their use as the bounded statistic of an anytime-valid change detector and the volatility-invariance argument for choosing them.

## 5. Implementation

`corrmon.py` (about 520 lines including docstrings, numpy only, imports `avmonitor.detectors.DetectorResult` and the MAD constant):

* `sign_concordance`, `truncated_product`, `gaussian_sign_mean`, `rho_from_sign_mean`, `gaussian_product_moments`: the statistics and their Gaussian means.
* `calibrate_pair(r1, r2, stat, c, band, margin_z) -> PairBaseline`: scales, `m0`, null values `m_up`, `m_dn`, bounds and the largest admissible bets.
* `bounded_edetector(y, m0, lo, hi, direction, threshold, fractions)`: vectorised Shiryaev-Roberts (or CUSUM) betting e-detector for a bounded statistic with an additive null; `m0` may be one value per path; NaNs are skipped.
* `gaussian_corr_lr_edetector` and `cusum_additive` + `calibrate_cusum_threshold`: the two benchmarks.
* `CorrelationMonitor(threshold, stat, c, cal_len, null_band, benchmark)`: the streaming calibrate, monitor, alarm, recalibrate loop, same protocol as `avmonitor.monitor.VolatilityMonitor`, each side at `2A`, events as `PairEvent` (index, direction, `m0`, implied correlation, scales).
* `aggregate_returns(logret, bars)`: non-overlapping 5-minute (or any) bars from 1-minute log returns, NaN unless every minute is present.

`simulate_pairs.py`: Gaussian or multivariate Student t pairs (common mixing variable, so tail dependence), a correlation change at a known step, independent volatility changes per asset, independent jumps and co-jumps.

Cost: one step of the streaming monitor is a handful of numpy operations on a 6-vector, about 15 microseconds; two years of 5-minute bars for one pair take about 3 seconds. The vectorised detector handles 2000 paths at once.

## 6. Experimental methodology

Synthetic settings mirror E4 of the main report: per-step units, a 2000-step calibration window drawn from the same process as the pre-change data (so jumps and fat tails are part of "normal"), `A = 7800` for 390-step days (at most 5% false alarms per day), both sides at `2A`, fixed seeds.

* **Q0 theory table [SYNTHETIC].** Null mean, variance, admissible bet, best grid bet and the predicted delay `log(2A) / growth rate` for six correlation changes, for the sign statistic and the product at `c = 3` and `c = 2`. Says which statistic should win before any detection run.
* **Q1 validity [SYNTHETIC].** 2000 paths of 5 days under eight scenarios: Gaussian (`rho` 0.3 and 0.6), Student t4, independent jumps, co-jumps, and three pure volatility changes at step 500 (asset 1 x2, asset 1 x0.5, both x2) with the correlation fixed. Reported: P(false alarm within 1 day) against the bound 5%, P(within 5 days) against 25%, and the share of alarms that are "up". Four detectors: sign, product (`c = 3`), Gaussian LR, naive CUSUM (threshold tuned on Gaussian pairs to 2.5% per side per day). The e-detectors are run with band 0 and 0.05 and, in a second table, with the true null mean (2 million draws) to check the guarantee free of baseline error.
* **Q2 delay [SYNTHETIC].** 1000 paths, change at step 500, 3500 post-change steps, six transitions (0.3 to 0.6, 0.6 to 0.3, 0.0 to 0.3, 0.3 to 0.0, 0.5 to 0.9, 0.6 to 0.0) and three innovation types (Gaussian, t4, jumps). Reported: false alarms before the change, detection share, share of alarms in the right direction, median, mean and 90th percentile delay. Five detectors (the product also at `c = 2`).
* **Q3 real data [REAL].** Binance spot 1-minute bars for BTC, ETH, SOL and DOGE against USDT, 1 January 2024 to 31 December 2025, from the mirror described in `avmonitor.data`. (a) Epps effect: Pearson and sign-implied correlation at 1, 5, 15, 30 and 60 minutes, by year. (b) The streaming monitor on 5-minute bars (sum of five 1-minute log returns, NaN if any minute is missing), deseasonalised by the 20-day rolling intraday profile of `avmonitor.data.intraday_scale` (this only matters for the product and Gaussian monitors, see 2.4), `A = 5760` (288 bars per day, 5% per day), calibration 1440 bars (5 days), bands 0, 0.1, 0.2. For every alarm: the sign mean and Pearson correlation over the 5 days after against the 5 days before, the BTC robust scale after / before (does the alarm coincide with a volatility change?), and the BTC return over the 24 hours before and after the alarm. The same quantities at 2000 random times are the control. Days on which all three pairs raise an increase alarm are listed with the BTC daily return. A 1-minute run of the sign monitor on BTC/ETH shows the effect of the bar length.

What is deliberately not done: no tuning on the real data (all settings were fixed from the synthetic runs and the main report), no claim about specific news events (dates are listed, not interpreted), and no trading backtest.

## 7. Folder structure

```
alternatives/correlation_monitor/
  README.md            this file: problem, formulation, proofs, comparison, implementation, methodology
  RESULTS.md           every number, labelled [SYNTHETIC] or [REAL], and an honest ranking
  corrmon.py           statistics, calibration, bounded betting e-detector, benchmarks, streaming monitor
  simulate_pairs.py    synthetic bivariate returns with correlation and volatility changes, jumps, t tails
  run_experiment.py    Q0 to Q3; writes results/ (about 5 minutes on one core; --only q0,q1,q2 or --only q3)
  tests/test_corrmon.py   9 synthetic tests, about 3 seconds, no data files or network needed (CI safe)
  results/             CSV outputs and two figures (fig1 synthetic delays, fig2 real alarms)
```

Run from the repository root: `python -m pytest -q alternatives/correlation_monitor/tests` and `python alternatives/correlation_monitor/run_experiment.py`.

## References

* Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. International Review of Economics and Finance, 109, 105482.
* Blomqvist, N. (1950). On a measure of dependence between two random variables. Annals of Mathematical Statistics, 21(4), 593 to 600. doi:10.1214/aoms/1177729754
* Chu, C.-S. J., Stinchcombe, M. and White, H. (1996). Monitoring structural change. Econometrica, 64(5), 1045 to 1065.
* Epps, T. W. (1979). Comovements in stock prices in the very short run. Journal of the American Statistical Association, 74(366), 291 to 298. doi:10.2307/2286325
* Hayashi, T. and Yoshida, N. (2005). On covariance estimation of non-synchronously observed diffusion processes. Bernoulli, 11(2), 359 to 379. doi:10.3150/bj/1116340299
* Lindskog, F., McNeil, A. and Nyfeler, U. (2003). Kendall's tau for elliptical distributions. In Credit Risk: Measurement, Evaluation and Management, Physica-Verlag, 149 to 156.
* Longin, F. and Solnik, B. (2001). Extreme correlation of international equity markets. Journal of Finance, 56(2), 649 to 676. doi:10.1111/0022-1082.00340
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. New England Journal of Statistics in Data Science, 2(2), 229 to 260.
* Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. Journal of the Royal Statistical Society Series B, 86(1), 1 to 27.
* Wied, D., Kraemer, W. and Dehling, H. (2012). Testing for a change in correlation at an unknown point in time using an extended functional delta method. Econometric Theory, 28(3), 570 to 589. doi:10.1017/S0266466611000661

The six references not already in `docs/literature_review.md` (Blomqvist, Chu et al., Epps, Hayashi and Yoshida, Longin and Solnik, Wied et al.) were checked by web search on 2 November 2026 (publisher or Project Euclid records).
