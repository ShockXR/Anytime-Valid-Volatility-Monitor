## ALT5. Two-sample betting: can the monitor compare the recent window with the calibration window without estimating a baseline?

Code, tests, README (formulation, proof, comparison with the literature) and full results are in `alternatives/two_sample_betting/`. Every number below is copied from its `results/` folder. [SYNTHETIC] marks simulated data, [REAL] marks Binance spot one-minute bars (BTCUSDT and ETHUSDT, 1 July to 31 December 2025, 184 days each). Iteration notes: `docs/iterations/ALT5.md`.

**Question.** The main monitor bets that the mean of the truncated squared return equals a baseline `m0` estimated from the calibration window; its guarantee is proved for a known `m0`. Alternative A4 asks for a sequential two-sample test by betting that detects any change in the law of the returns, with nothing estimated. Is it exactly valid, how much speed does it give up on variance changes, what does it see that a variance monitor cannot, and what does it flag on real data?

**Method.** Each new return is ranked among all earlier returns of the segment (calibration window plus monitored returns, a growing reference), giving a smoothed conformal p-value. Under exchangeability these p-values are exactly i.i.d. uniform (Vovk, Nouretdinov and Gammerman 2003), so any bet `1 + lam (w(p) - b E_U|w|)` with a predictable witness `w` of mean 0 under the uniform law is fair, and the Shiryaev-Roberts mixture of such bets is an e-detector with the guarantees of `docs/math.md` 5.4 (average run length at least `A`, false alarm within `m` steps at most `m / A`), in the filtration of the p-values. The null band `b` allows the density of each p-value to be within `1 +- b` of uniform. Three bettor families: scale up and scale down bets on the Gaussianised truncated score `min(z^2, 9) / kappa2` of the p-value (the main monitor's statistic, but with null mean exactly 1), and histogram learners that estimate the density ratio of new to reference p-values (the rank-space analogue of a kernel two-sample witness, Shekhar and Ramdas 2024). Weights 0.25, 0.25 and 0.5, frozen from a pilot on development seeds. Truncation does not matter for ranks (monotone invariance, ties handled exactly) except to bound the Gaussianised decrease bet. Volatility clustering breaks exchangeability, as it breaks the main monitor's null. `A = 28,800` everywhere (5% per 1440-step day).

**Validity with no change [SYNTHETIC].** 1000 paths per law; share with a false alarm within 1 day (bound 5%) and 7 days (bound 35%), one-day calibration window, and within 1 day with a two-hour window:

| Detector | Gaussian, 1 d / 7 d | Student t4, 1 d / 7 d | Jumps, 1 d / 7 d | Two-hour window, 1 d (three laws) |
|---|---|---|---|---|
| main monitor, plug-in, band 0 | 5.4% / **57.0%** | 6.0% / **52.5%** | 4.6% / **51.4%** | **42.0% to 45.7%** |
| rank detector, omnibus, b 0 | 2.4% / 22.6% | 3.1% / 24.2% | 2.4% / 23.0% | 2.2% to 3.5% |
| rank detector, fixed reference (240 / 1440 steps) | 90.4% / 100% and 21.2% / 98.4% (Gaussian, 500 paths) | | | |

The plug-in baseline's error acts as a constant drift that the bets find within days; ranking against the calibration window only (quantile bins learned once) fails the same way, because the p-values then share one reference and are dependent. With the growing reference every rate is below its bound at 1, 7 and 14 days. The repo's proved variant (`baseline="cs"`) is valid but conservative: 0.0% per day and 0.1% per week.

**Power [SYNTHETIC].** Change at monitored step 500, 400 paths, median delay in steps [share detected within 4000 steps when below 95%]:

| Change | main, plug-in, band 0 | main, confidence-sequence baseline | rank, scale bets only | rank, omnibus |
|---|---|---|---|---|
| sd x1.25 | 127 | 348 | 136 | 152 |
| sd x0.8 | 211 | 400 | 222 | 237 |
| Gaussian to t4, same variance | 1203 | 2251 [10%] | 1404 [86%] | **417** |
| jumps 0.2% to 1% per step, same variance | **361** | 1129.5 [92%] | 418 | 399 |
| jumps 0.2% to 0.5% per step, same variance | **1393** [88%] | not run | 1698 [57%] | 1688 [54%] |

The signed score (rank of `u` rather than `|u|`) alone detects a pure skew change (809 steps, 99%); for every `|u|`-based detector the law of the statistic does not change there. Truncating the rank score at 3 sd moves delays by at most about 6%. With a 120-step reference the rank detector stays valid (at most 0.8% false alarms before the change, bound 1.7%) and is 17% to 52% slower; the plug-in monitor with the same window has 11.5% to 14.0%. Under GARCH(1,1) clustering both designs alarm on most days (band 0.25, persistence 0.99: main 67.0%, rank 53.0%).

**Real data [REAL].** Returns divided by `intraday_scale`, one-day calibration, recalibration after each alarm. At band 0.25 the main monitor raises 127 (BTC) and 122 (ETH) alarms, the rank detector 123 and 114; 50% and 61% of the rank alarms are matched by a main-monitor alarm within an hour (chance 5.7%). Each monitor is calibrating about two thirds of the time, and 44 of the 61 unmatched BTC rank alarms (40 of 45 for ETH) fell while the main monitor was calibrating. While the main monitor was watching, the rank detector alone raised 17 BTC and 5 ETH alarms. The BTC ones are shape changes at constant scale: robust scale ratio 0.994, scale-free shape distance 0.131 against 0.076 at random times, and a learned density ratio of 0.32 for the largest 0.2% of moves and about 1.5 for the smallest ones, with more zero returns (5.0% against 3.8% in the reference). Dropping zero returns halves the group to 8, so part of it is a change in price stickiness rather than volatility. The main monitor alone raised 32 alarms per coin while the rank detector watched, all volatility decreases (truncated variance ratio 0.50 and 0.52). With a 120-minute reference the rank detector monitors 88% to 91% of the time.

**Reading.**
* The rank detector has the strongest guarantee of the designs tested: exact for any return law, any reference length and any horizon. The plug-in main monitor holds its bound for a day but not for a week (51% to 57% against 35%), which matters for every deployment that does not recalibrate daily.
* On variance changes it costs 4% to 7% in delay with scale bets alone and 11% to 19% with the full mixture; against the proved confidence-sequence variant it is 1.7 to 2.8 times faster. It is 2.9 times faster than the plug-in monitor for a tail change at fixed variance, which the variance monitors barely see. A change of jump intensity at fixed total variance is a negative result: the variance monitor is faster, because the diffusive variance falls.
* On BTC and ETH the two designs mostly agree; the rank detector adds a few shape alarms and misses some decreases at this band. Its practical gain is the short reference, which the plug-in design cannot use.

**Caveats.** Exchangeability fails on real one-minute returns (clustering), so real alarms are descriptive, not tested discoveries; ranks relative to a volatility forecast are the natural next step and were not built. The bands of the two designs are not on the same scale, so the real-data comparison is between two settings, not two designs in general. The rank detector gives no direction or size for shape alarms beyond the density ratio. About 80 microseconds per step in pure Python against about 10 for the main monitor (`results/t_timing.csv`). The conjecture that the Gaussianised score is asymptotically as efficient as the main monitor under Gaussian data is untested. One venue, two coins, six months.
