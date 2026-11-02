# Alternative prototype directions

These are other ways to build on Bapat (2026). Each one was first written down as an idea (round 1) and then built as a small, separate prototype (round 2) with its own README, code, tests and RESULTS file. They are kept apart from the main package on purpose, except A2, which proved useful enough to move into `src/avmonitor/confseq.py`.

Every number below is copied from the RESULTS file of the prototype, which in turn is produced by its `run_experiment.py` with fixed seeds. Labels: [SYNTHETIC] simulation, [REAL] market data.

## Summary table

Scores run from 1 (weak) to 10 (strong). They are our judgement, explained in each entry, not measurements. "Feasibility" means how easy it would be to turn the prototype into a finished study. "Interest" means our guess of how interesting the direction is to a researcher in sequential analysis and high-frequency econometrics.

| | Direction | Folder | Novelty | Validity | Relevance | Feasibility | Demonstration | Interest | Total (of 60) |
|---|---|---|---|---|---|---|---|---|---|
| A7 | Correlation-break monitor | `correlation_monitor/` | 7 | 9 | 8 | 9 | 8 | 7 | **48** |
| A2 | Confidence-sequence baseline | `src/avmonitor/confseq.py` (E8) | 6 | 10 | 7 | 9 | 6 | 8 | **46** |
| A5 | Robust local plug-in estimator of quarticity | `robust_iq_estimator/` | 7 | 7 | 8 | 7 | 7 | 9 | **45** |
| A4 | Rank (conformal) two-sample e-detector | `two_sample_betting/` | 6 | 10 | 7 | 8 | 7 | 7 | **45** |
| A1 | Pre-averaging before truncation, block monitor | `preaveraging/` | 6 | 8 | 6 | 9 | 8 | 7 | **44** |
| A3 | Quarticity-stream monitor | `quarticity_stream/` | 5 | 9 | 5 | 9 | 6 | 8 | **42** |
| A6 | Change-triggered risk switch | `risk_switch/` | 3 | 7 | 6 | 9 | 6 | 3 | **34** |

For comparison, the main direction (the robust volatility monitor itself) scores 57 of 70 on the seven criteria used in `RESEARCH_STRATEGY_BAPAT.md`, section 8.

**Recommended follow-ups.** A5 is the natural estimation paper with Prof. Bapat, because it stays inside his estimator and answers his own future-work items. A7 is the natural second detection paper, because it reuses the same machinery for the multivariate question the paper lists as future work. A1 matters for any use of the monitor at second-level sampling, because it found a failure of the main monitor's guarantee and a provable fix.

**One lesson that appears in four directions.** The e-detector's guarantee needs the *conditional* mean of the statistic, given the past, to stay at the baseline. Anything that makes consecutive values of the statistic depend on each other breaks it, even when the *unconditional* mean is right: volatility clustering (E10), price-level noise and overlapping pre-averaged windows (A1), a reference window that is reused for every new observation (A4) and a plug-in baseline (E4e). Each prototype found its own fix: a forecast, non-overlapping blocks, a growing reference, a confidence sequence.

---

## A1. Noise-robust first stage: pre-averaging, then truncation

**Question.** Truncation removes jumps but not microstructure noise. Does averaging returns over small blocks first (Jacod et al. 2009) fix that, for the estimator and for the monitor?

**Results [SYNTHETIC].**
* The plain truncated quarticity estimator breaks between noise ratios 0.1 and 0.25: its bias is exactly `(1 + 2 gamma^2)^2`, 26% at 0.25 and a factor 9 at 1. Pre-averaging then truncation is unbiased at every noise level (mean ratio 0.99 to 1.00), at the price of a relative RMSE of 0.11 against 0.02 when there is no noise.
* After pre-averaging, the paper's *clipping* leaves a jump bias that grows with the window (+12% to +18% with four jumps a day); *dropping* the window removes it.
* **A correction to the main package:** price-level noise breaks the monitor's conditional null on raw returns. Its one-day false alarm rate rises from 5.0% to 8.5%, 18%, 31% and 41.5% at noise ratios 0.5, 1, 2 and 4. Overlapping pre-averaged returns are worse (100%). Non-overlapping blocks restore the guarantee at every block size tested (0% to 5.5%), because they use disjoint prices.

**Scores.** Novelty 6: the estimator is known, the composition as a monitor and the noise failure are new to us. Validity 8: the estimator's fourth-moment correction is derived under local constancy, a heuristic. Relevance 6: high at one-second sampling, low at one minute for BTC. Feasibility 9, Demonstration 8 (two figures make the point), Interest 7.

## A2. Confidence-sequence baseline (now in the main package)

**Question.** The default monitor plugs in an estimated baseline, which voids the proof. Can a running confidence sequence for the pre-change level restore it?

**Results.** Yes, with a proof: P(false alarm within m steps) ≤ δ + m / A (`docs/math.md`, Theorem 6.4). In simulation the false alarm rate is at most 0.27% against a 5% bound, at a cost of 1.8 to 7.2 times the delay with a 390-step calibration (E8). On real data the share of decrease alarms followed by lower volatility rose from 90% to 97%.

**Scores.** Novelty 6 (betting confidence sequences exist; their use as a detector baseline with this bound is ours), Validity 10, Relevance 7, Feasibility 9, Demonstration 6, Interest 8.

## A3. Monitor quarticity itself

**Question.** Does the same e-detector on `min(r^4, tau^4)`, the paper's own target, see changes the variance monitor cannot?

**Results.** Valid in all 12 null cells [SYNTHETIC]. For a pure variance change it is 1% to 28% faster for increases and 1.8 to 2.6 times slower for decreases (its bets are capped). It detects changes in the volatility of volatility at a fixed mean variance (82.5% to 99% of paths, against 2% to 8% for the variance monitor), but slowly (thousands of steps). A surprise: at c = 3 fatter tails at fixed variance are seen by the *variance* monitor, as a decrease, because clipping removes the tails. On BTC and DOGE [REAL] it alarms at the same rate as the variance monitor (0.70 against 0.69 per day) and mostly on the same days; no pure shape change was found in six months.

**Scores.** Novelty 5, Validity 9, Relevance 5, Feasibility 9, Demonstration 6 (strong synthetic story, weak real-data story), Interest 8 (closest to the paper's own quantity, and the clipping-level effect is worth a discussion).

## A4. Rank (conformal) two-sample e-detector

**Question.** Instead of betting on the mean of the truncated statistic against an estimated baseline, bet on the *ranks* of new returns among all past returns. Ranks need no baseline at all.

**Results [SYNTHETIC].** Exactly valid for any return law, reference length and horizon: 1.9% to 3.5% per day (bound 5%) and 22.2% to 24.3% per week (bound 35%), where the plug-in mean monitor reaches 51% to 57% per week. Ranking against a *fixed* reference window instead makes a false alarm near certain within a week (98% to 100%): the same trap as the plug-in baseline. Against the main monitor it is 4% to 19% slower for variance changes, 2.9 times faster for a change from Gaussian to fat tails, and 1.7 to 2.8 times faster than the confidence-sequence version (A2). Truncation does not matter for ranks. [REAL] On BTC and ETH both designs raise similar numbers of alarms and about half coincide within an hour; the rank detector alone found shape changes at constant scale on BTC.

**Scores.** Novelty 6 (conformal test martingales and two-sample betting exist; the combination and the fixed-reference demonstration are ours), Validity 10 (the strongest guarantee in the repository), Relevance 7, Feasibility 8, Demonstration 7, Interest 7.

## A5. Robust local plug-in estimator of integrated quarticity

**Question.** Averaging local fourth powers (which is what the paper's estimator does once averaged) is known to be inefficient. Does squaring a local variance instead (Jacod and Rosenbaum 2013), with robust truncation, give a better estimator?

**Results.**
* [SYNTHETIC] The plug-in with a growing cut-off `sqrt(2 log n)` and dropped returns (PI-LU) has asymptotic variance 8.2 against the efficient bound 8 (efficiency 0.95), against 0.74 for RQ-type estimators. With jumps it has the lowest RMSE of all 15 estimators at n = 390 and 1440, and its 95% intervals cover 93% to 97% in all 18 cells.
* [SYNTHETIC] Truncating against one daily MAD stops converging when volatility moves during the day (bias still minus 14% at n = 23,400 under stochastic volatility with an intraday pattern); a local scale restores the usual `sqrt(n)` rate. Under the paper's constant-volatility design this problem does not arise.
* [REAL] On BTC and DOGE, large one-minute returns come in bursts (44% to 65% have a large neighbour). Every estimator that truncates against a scale measured over hours, the paper's included, removes a third to a half of the quarticity as a result. A neighbour-guarded version (PI-NN) is designed for this; its settings were chosen after seeing the data, so it still needs an out-of-sample check.

**Scores.** Novelty 7 (each ingredient is known; exact finite-sample corrections, the analysis of the MAD rule and the burst finding are new to us), Validity 7 (no central limit theorem proved; coverage is numerical), Relevance 8, Feasibility 7 (a paper needs a CLT), Demonstration 7, Interest 9 (it is the paper's own estimation problem). This direction says the most about the published estimator, so it must be presented as "what the MAD idea needs in order to be consistent and efficient", not as a correction.

## A6. Change-triggered risk switch

**Question.** Does scaling a position down after an increase alarm (and up after a decrease alarm) control risk better than ordinary volatility targeting?

**Results [REAL].** Tails are on a par with daily EWMA targeting (2025 crypto portfolio: maximum drawdown 34.4% against 33.8%, and 54.6% for buy and hold; Nifty over nine years: 28.1% against 30.8% and 39.3%) with far fewer trades (93 against 365 changes a year per coin). But it keeps delivered risk further from the target (weekly volatility 0.27 log units off target against 0.19 to 0.22), and Sharpe ratios cannot be told apart.

**Scores.** Novelty 3, Validity 7 (an honest backtest with costs and paired bootstrap intervals), Relevance 6, Feasibility 9, Demonstration 6, Interest 3 (a trading rule is far from the paper's topic). Kept as a negative result.

## A7. Correlation-break monitor

**Question.** The paper names a multivariate extension as future work. Can the same machinery detect breaks in the dependence between two assets, online and with a guarantee, without being fooled by volatility changes?

**Results.** A betting e-detector on the sign concordance `sign(r1) sign(r2)` is bounded and scale-free, so pure volatility changes and persistent GARCH clustering do not trigger it (false alarms 2.3% to 3.1% against 5% [SYNTHETIC]), while every other detector tested is triggered. It is also faster than a product of standardised returns truncated at c = 3, whose bets are capped by the truncation bound (174 against 246 steps for a correlation change from 0.3 to 0.6). Truncating the product at c = 2 is a little faster still (158 steps) but inherits the product's sensitivity to volatility changes. A basket version monitors several pairs at once. On Binance coins [REAL] it finds correlation breaks about once a month per pair, clustered after sell-offs.

**Scores.** Novelty 7, Validity 9, Relevance 8 (correlation breaks after sell-offs are the risk diversification does not cover), Feasibility 9, Demonstration 8, Interest 7.
