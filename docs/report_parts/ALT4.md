## ALT4. Can the recipe detect a change in the correlation between two assets, without mistaking a volatility change for it?

Prototype in `alternatives/correlation_monitor/` (README with the formulation and proofs, `RESULTS.md` with every number, Q0 to Q4). Script `run_experiment.py`, results in `alternatives/correlation_monitor/results/`.

**Question.** Bapat (2026) names a multivariate extension as future work. The product of two assets' returns moves when their correlation changes, but also when the volatility of either asset changes or when one of them jumps. Is there a bounded statistic whose conditional mean moves with the correlation only, so that the e-detector of section 5 of `docs/math.md` applies unchanged?

**Method.** Two candidates. (i) The truncated product of standardised returns, `clip(r1 / s1, -c, c) clip(r2 / s2, -c, c)`, bounded in [-c^2, c^2] (the literal "truncated covariance"). (ii) The sign concordance `sign(r1) sign(r2)`, bounded in [-1, 1], with mean `(2 / pi) arcsin(rho)` for any elliptical pair, and invariant to any positive predictable rescaling of either asset (a two-line proof in the README). Both run in the Shiryaev-Roberts betting e-detector with an additive null band, each side at `2A`, so Theorems A and B hold verbatim. Benchmarks: a bivariate Gaussian likelihood ratio e-detector and a Page CUSUM on the raw product tuned to 5% per day under Gaussian pairs.

**Results [SYNTHETIC]** (2000 paths per row, A = 7800, 2000-step calibration, band 0). Share of paths with a false alarm within one day (bound 5%) and within five days (bound 25%):

| Scenario (correlation 0.3 throughout) | sign | product c3 | Gaussian LR | naive CUSUM |
|---|---|---|---|---|
| Gaussian | 2.9% / 23.0% | 0.0% / 11.3% | 5.0% / 32.0% | 3.2% / 18.4% |
| Student t4 | 2.4% / 23.1% | 0.2% / 15.7% | 79.8% / 100% | 96.2% / 100% |
| Jumps (0.5% per step, 8 sd) | 3.1% / 24.1% | 0.0% / 12.4% | 57.0% / 98.3% | 22.4% / 74.0% |
| Volatility of asset 1 doubles at step 500 | 2.6% / 23.6% | 0.0% / 98.3% (all "up") | 5.1% / 99.9% (97% "down") | 3.8% / 100% |
| Volatility of asset 1 halves at step 500 | 3.1% / 23.8% | 0.0% / 100% (all "down") | 4.6% / 79.6% | 2.7% / 48.2% |
| GARCH(1,1) clustering in each asset (alpha 0.05, beta 0.94), correlation fixed | 2.9% / 23.9% | 2.1% / 63.6% | 25.4% / 83.7% | 33.3% / 87.2% |

Median delay for a correlation change from 0.3 to 0.6 at step 500 (1000 paths): sign 174 steps, product c3 246, product c2 158, Gaussian LR 74 (but 70% to 100% false alarms before the change under t4 or jumps), naive CUSUM 125. For 0.5 to 0.9: sign 65, product c3 194, Gaussian LR 20. A growth-rate calculation predicts these delays within 25% and shows why the product is slow: its bound of 9 caps the admissible bet at 0.11, less than half the growth-optimal bet.

**Results [REAL]** (Binance 5-minute bars, 2024 to 2025, A = 5760 for 5% per day, calibration 5 days, band 0.1). The sign monitor raises 23 to 24 alarms per pair in two years. All 31 increase alarms are followed by a higher Pearson correlation over the next 5 days than in the calibration window (47% to 49% at random times), and 75% to 92% of the decrease alarms by a lower one. Increase alarms follow a negative 24-hour BTC return in 80% of cases and a bottom-decile return in 36% to 50% (10% by construction). The product monitor's decrease alarms coincide with drops in BTC volatility (median after / before ratio 0.83, against 0.99 at random times), the confusion predicted by the simulation; the Gaussian LR raises 4 to 5 times more alarms, 78% to 88% of them "down" and close to chance. At one minute the Epps effect (Pearson 0.02 to 0.06 lower than at 30 minutes) triples the number of "down" alarms, so 5-minute bars are used.

**Basket of four assets [SYNTHETIC and REAL].** The average pairwise sign concordance of k assets is bounded in [-1/(k-1), 1], so one e-detector at level A covers the whole basket with the same guarantee. Against the union of six pairwise monitors at level 6A it is twice as fast for a basket-wide change (72 against 132 steps for all correlations 0.3 to 0.6) and five times slower for a change confined to one pair (451 against 91), and it keeps its bound under t4 and volatility changes (19% to 24% per 5 days against 25%). On BTC, ETH, SOL and DOGE (raw 5-minute returns, band 0.1) it raises 31 alarms in two years; 93% of increase and 94% of decrease alarms are confirmed over the next 5 days (48% at random), and 8 of the 15 increase alarms follow a bottom-decile 24-hour BTC return.

**Reading.** The sign statistic turns the volatility monitor into a correlation monitor with the same guarantee, with no scale estimate and no truncation needed, and it is the only detector tested that a pure volatility change cannot trigger. It pays for this with a delay 2 to 3 times that of the Gaussian-optimal detector on clean Gaussian data, and is about as fast as the truncated product at c = 3, whose bets are throttled by Bapat's truncation bound.

**Caveats.** Plug-in baseline (a band of 0.05 to 0.1 or the confidence-sequence baseline of section 6 is needed; with band 0 the 5-day false alarm rate is 23% against a bound of 25%). The 5-day confirmation on real data can partly reflect mean reversion of the rolling correlation. The basket statistic dilutes a change confined to one pair; weighted variants and a forecast-relative null were not built.
