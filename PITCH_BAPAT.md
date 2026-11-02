# Pitch for Prof. Sudeep Bapat

*A short brief for a first meeting, written as answers to ten questions. Every number has a label: [PUBLISHED] from your paper, [SYNTHETIC] from our simulations, [REAL] from our runs on market data. Details are in the [experiment report](docs/experiment_report.md); proofs are in [`docs/math.md`](docs/math.md) and [`docs/theory_delay.md`](docs/theory_delay.md).*

## 1. What problem are we solving?

We want an **online alarm for volatility changes**: a tool that reads one-minute returns as they arrive and says "volatility has gone up" or "volatility has gone down" soon after it happens, while **guaranteeing** how rarely it raises a false alarm, however long it runs and however often someone looks at it.

## 2. Why does it matter?

Risk limits, position sizes, margin requirements and model refits should all react when volatility changes regime. Two things make an online alarm hard. First, a few jumps or fat-tailed returns distort ordinary volatility measures, so classical detectors fire on jumps rather than on regime changes. Second, a test that is checked every minute will eventually fire by chance unless it is built for repeated looks. Existing tools solve one problem or the other, not both with a finite-sample guarantee.

## 3. What did the paper propose?

Bapat (2026) estimates integrated quarticity (the quantity behind confidence intervals for realized volatility) robustly: each return is clipped at c times its median absolute deviation (MAD), and the clipped fourth powers are smoothed with a kernel. Clipping bounds the influence of any single return, so jumps and outliers cannot dominate. The paper shows good performance in simulations with jumps and noise and on one-minute data for five cryptocurrencies.

## 4. What did we implement?

* Your estimators exactly as printed, normalised versions, and all the benchmark estimators of your Tables 3 and 4 (`src/avmonitor/robust.py`).
* A streaming monitor (`src/avmonitor/monitor.py`) that clips each squared return at c × MAD from past data only, X_t = min(r_t², τ²), and runs **betting e-detectors** on X_t, one for increases and one for decreases. An e-detector keeps a running "wealth" from bets against "nothing has changed"; it alarms when the wealth reaches A.
* Three ways to set the "normal" level: plug-in, plug-in with a margin, and a confidence sequence with a proof (`confseq.py`); forecast-relative monitoring for clustered volatility (`forecast.py`).
* Real-data pipelines for Binance BTC, ETH, XRP, SOL and DOGE (your five coins), Bitstamp BTC/USD and the Nifty 50 index, with a command line and an executed notebook. About 190 tests.

## 5. What did we independently validate?

* **Your tables.** With one fitted volatility, the printed formulas reproduce your simulation tables for the plain kernel, and the published crypto table is consistent with unnormalised kernel weights, the raw MAD and a sample with about 2,000 missing minutes (E1, E6) [PUBLISHED vs SYNTHETIC and REAL]. Once both kernels are put on the same scale, **your robust kernel beats the plain kernel in most simulation cells** (E1).
* **Our guarantee.** Proved: average time to a false alarm at least A, and P(false alarm within m steps) ≤ m / A (`docs/math.md`, Theorems A and B). Measured: 1.1% to 1.4% per day against a 5% bound for four very different return laws (E9) [SYNTHETIC]. The bound itself is not loose: an adversarial but legitimate return law reaches 95% of it for a single bet (E11, exact dynamic programming).
* **Robustness.** With jumps or fat tails, Gaussian likelihood-ratio and CUSUM detectors raise a false alarm within a day on 70% to 97% of paths; the robust e-detector stays at 0.3% to 2.3% and detects a 50% volatility rise in a median of 36 to 49 steps (E4) [SYNTHETIC].

## 6. What did our experiments demonstrate?

* **The clipping step of your estimator is what makes a model-free alarm possible.** No nonparametric test can detect a *fall* in the variance itself (a one-sided form of Bahadur and Savage 1956, `docs/math.md` 5.3), but it can detect a fall in a bounded quantity such as the clipped variance.
* **The price of robustness is small and explained.** On clean normal data the robust detector is 14%, 35% and 73% slower than a Gaussian likelihood-ratio detector for volatility rises of 25%, 50% and 100%. A ratio of information numbers predicts these penalties (1.16, 1.37, 1.95), and splits them into a truncation cost, a cost of making no distributional assumption (the largest part) and a cost of a finite set of bets (2% to 4%) (E11) [SYNTHETIC].
* **On real data the alarms mean something.** On your five coins, alarms are followed by volatility on the flagged side of EWMA and HAR forecasts 71% to 83% of the time against about 50% at random times, in and out of sample (E6). 14 to 15 of 16 FOMC statements raise an alarm within six hours, against about 30% of matched control times (E12). On Nifty 50, 9 of 11 Union Budget sessions raise an alarm against 25% of ordinary sessions, and 78% of increase alarms in a held-out period are followed by a more volatile session against 47% at random (E7) [REAL].

## 7. What did we discover?

* **Clipping sets a speed limit for detecting decreases**: no valid detector can gather more than about 1 / c² of evidence per step against a fall, so at c = 3 even a total collapse of volatility takes at least 73 steps. A tighter clip on the decrease side (c ≈ 1.5 to 2) cuts decrease delays by 26% to 49% (E11) [SYNTHETIC].
* **An estimated baseline quietly breaks the guarantee over long horizons**: 44% false alarms against a guaranteed 25% over five days (E10), 51% to 57% against 35% over seven days (ALT5). A confidence-sequence baseline restores a proof, P(false alarm within m) ≤ δ + m / A (E8).
* **Volatility clustering is itself a change in conditional variance.** Dividing by a fitted GARCH forecast restores the guarantee, but a forecast that keeps adapting lets through at most (1 − β)/α of any level shift, a fact we derived in closed form: a doubling of volatility looks like a 9% rise at persistence 0.99 (E10) [SYNTHETIC].
* **Price-level microstructure noise breaks the guarantee on raw returns** (false alarms rise from 5% to 41.5% per day as noise grows), and non-overlapping pre-averaged blocks restore it (ALT1) [SYNTHETIC].
* **On one-minute crypto data, large returns come in bursts**, not as isolated jumps: 44% to 65% have a large neighbour. Any estimator that clips against a scale measured over hours then removes a third to a half of the quarticity; a local scale fixes most of it (ALT3) [REAL and SYNTHETIC].
* **The honest bottom line on real data**: a simple ratio rule tuned with hindsight is about as precise as the e-detector (E6e). What the e-detector adds is a false-alarm guarantee with settings fixed in advance, not extra precision.

## 8. What are the limitations?

* The guarantee needs the clipped squared returns to keep their conditional mean while nothing changes. Clustering, price noise and estimated baselines violate this; each has a fix, and each fix costs speed.
* The bound is about four times conservative with the default set of bets (E9, E11).
* After an alarm the monitor recalibrates and is blind about 70% of the time on real crypto data (E12); a restart without a blind period is in progress (E14).
* Our crypto and Nifty data come from public mirrors, not from the exchanges, and your paper does not name its data source.
* Full list: [`docs/limitations.md`](docs/limitations.md).

## 9. What should be done next?

1. A short joint paper on the monitor (a draft is in `paper/`): the bounded-statistic argument, the delay theory and the decrease speed limit, the confidence-sequence baseline, and the crypto and Nifty applications.
2. An estimation follow-up on your estimator: a local-scale plug-in version that reaches 95% efficiency in simulation (ALT3), which needs a central limit theorem.
3. A sequential-estimation bridge to your own field: purely sequential fixed-accuracy estimation of volatility with your quarticity estimator as the variance plug-in, compared with an anytime-valid confidence sequence (ALT6, in progress).
4. Panel monitoring of many assets with false-discovery control (E15) and a benchmark against published change-point methods (E13), both in progress.

## 10. What feedback or collaboration are we seeking?

**Questions about the paper (asked to learn, not to correct).** Each is documented with evidence in [`docs/estimator_issues.md`](docs/estimator_issues.md).

1. Are the kernel weights normalised to sum to one? Your Table 6 fits the raw weight sum (Kernel / RQ ≈ 4.0 × 10⁻⁵ for all five coins).
2. Is the MAD multiplied by 1.4826? The tables fit the unscaled MAD, which at c = 3 clips at about 2 standard deviations.
3. Which per-step volatility did the simulations use (σ ≈ 0.0215 fits most of Table 7)? Are the Table 8 noise levels variances or standard deviations? What differs in the σ_J = 0.10 design, where our robust-kernel numbers come out lower than yours?
4. How was the RMSE in Table 6 computed? It matches |estimate − 2.2 × 10⁻⁶| for every coin.
5. Table 5 gives Q1 = Q3 = 0 for DOGE, which makes the MAD 0 and the clipped kernel exactly 0, yet Table 6 prints 2.00 × 10⁻¹² (the same value as XRP). Could the two rows come from different series? And the ETH mean in Table 5 has the size we find but the opposite sign: could that be a typo?
6. In the crypto application, was the MAD computed per day, over the whole sample, or locally?

**Collaboration.** Your view on what the monitor should watch (clipped variance, clipped quarticity, or both); data for Indian markets at one minute or finer; and whether the monitor paper (item 9.1) or the estimation follow-up (9.2) could become joint work.

## Five-minute demo

```bash
pip install -e ".[dev]"
python -m pytest -q                                   # about 190 tests, no data needed
python -m avmonitor simulate --ratio 2                # synthetic stream with a known change
python -m avmonitor run --source binance --symbol BTCUSDT --download \
    --start 2025-10-08 --end 2025-10-14 --plot btc.png  # real data around the 10 October 2025 crash
```
