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
* Real-data pipelines for Binance BTC, ETH, XRP, SOL and DOGE (your five coins), Bitstamp BTC/USD and the Nifty 50 index, with a command line and an executed notebook. About 240 tests, and one command (`scripts/reproduce.py`) that checks the data hashes, runs the tests and reruns the experiments.

## 5. What did we independently validate?

* **Your tables.** With one fitted volatility, the printed formulas reproduce your simulation tables for the plain kernel, and the published crypto table is consistent with unnormalised kernel weights, the raw MAD and a sample with about 2,000 missing minutes (E1, E6) [PUBLISHED vs SYNTHETIC and REAL]. Once both kernels are put on the same scale, **your robust kernel beats the plain kernel in most simulation cells** (E1).
* **Our guarantee.** Proved: average time to a false alarm at least A, and P(false alarm within m steps) ≤ m / A (`docs/math.md`, Theorems A and B). Measured: 1.1% to 1.4% per day against a 5% bound for four very different return laws (E9) [SYNTHETIC]. The bound itself is not loose: an adversarial but legitimate return law reaches 95% of it for a single bet (E11, exact dynamic programming).
* **Robustness.** With jumps or fat tails, Gaussian likelihood-ratio and CUSUM detectors raise a false alarm within a day on 70% to 97% of paths; the robust e-detector stays at 0.3% to 2.3% and detects a 50% volatility rise in a median of 36 to 49 steps (E4) [SYNTHETIC].

## 6. What did our experiments demonstrate?

* **The clipping step of your estimator is what makes a model-free alarm possible.** No nonparametric test can detect a *fall* in the variance itself (a one-sided form of Bahadur and Savage 1956, `docs/math.md` 5.3), but it can detect a fall in a bounded quantity such as the clipped variance.
* **The price of robustness is small and explained.** On clean normal data the robust detector is 14%, 35% and 73% slower than a Gaussian likelihood-ratio detector for volatility rises of 25%, 50% and 100%. A ratio of information numbers predicts these penalties (1.16, 1.37, 1.95), and splits them into a truncation cost, a cost of making no distributional assumption (the largest part) and a cost of a finite set of bets (2% to 4%) (E11) [SYNTHETIC].
* **On real data the alarms mean something.** On your five coins, alarms are followed by volatility on the flagged side of EWMA and HAR forecasts 71% to 83% of the time against about 50% at random times, in and out of sample (E6). Around FOMC statements, the streaming monitor (the design fixed in advance) flags 8 to 9 of 16 within six hours, against 11% to 17% of control windows on the same weekday and clock time: every statement at which it was monitoring rather than recalibrating. A monitor calibrated on the day before each statement flags 14 to 15 of 16, against 29% to 32% of controls (E12). On Nifty 50, 9 of 11 Union Budget sessions raise an alarm against 25% of ordinary sessions, and 78% of increase alarms in a held-out period are followed by a more volatile session against 47% at random (E7) [REAL].
* **Against published methods.** Bayesian online change-point detection (BOCPD, Adams and MacKay 2007) raises false alarms on 53% to 73% of days with fat tails or jumps at a threshold that gives 1% on normal data; our monitor stays at 1.2% to 1.3%. On clean data BOCPD detects decreases two to three times faster [SYNTHETIC]. On your five coins, with every method tuned to the same numbers of up and down alarms, our alarms are followed by a day beyond the forecast more often (78.5% and 80.5% against 58.7% to 71.7%; in a held-out half-year 82.0% and 69.7% against 56.7% to 71.5%, with one tie on down alarms); on Nifty no difference is significant (E13) [REAL]. The comparison with the volatility jump test of Bibinger, Jirak and Vetter (2017) uses our reconstruction of it, not yet checked against the paper's equations.

## 7. What did we discover?

* **Clipping sets a speed limit for detecting decreases**: clipping caps the evidence that any valid bet can gather against a fall at log(B / (B - 1)), about 1 / c² nats per step. At c = 3 and A = 7800, even a total collapse of volatility takes any Shiryaev-Roberts e-detector built from valid bets at least 58 steps (64 at the two-sided level 2A), and our default bets at least 73 (79). A separate clip for the decrease side (c_down = 2, keeping c = 3 for increases) cut simulated median decrease delays by 13% to 51% for Gaussian, Student t4 and jump returns, and by 32% to 52% for prices rounded to ticks of 0.25 and 0.75 standard deviations, but was 12% slower for a 20% fall at a tick of 0.5 standard deviations; on five coins, which have few zero returns, it was neutral to slightly positive. Its false alarm guarantee holds with a known baseline or with the confidence-sequence baseline, which now supports it; the simulations used the plug-in baseline, whose one-day false alarm rates (5.0% to 6.7%) are empirical only (E11, E11f) [SYNTHETIC and REAL].
* **An estimated baseline quietly breaks the guarantee over long horizons**: 41% to 44% false alarms against a guaranteed 25% over five days (E10, E14), 51% to 57% against 35% over seven days (ALT5). We proved that no bound of the form C m / A can hold for a plug-in baseline at a fixed calibration length, because a fixed share of calibration windows gives some bet a positive drift (E14, Proposition A.6). A confidence-sequence baseline restores a proof, P(false alarm within m) ≤ δ + m / A (E8), and a rank-based detector needs no baseline at all (E14).
* **The blind period after an alarm can be removed with a proof.** A rank (conformal) e-detector restarted at once after each alarm keeps an exact m / A bound and catches 98% to 99.6% of changes in synthetic streams with about one change a day, against 56% to 66% for the default protocol, which is blind 57% to 71% of the time (E14) [SYNTHETIC].
* **Volatility clustering is itself a change in conditional variance.** Dividing by a fitted GARCH forecast brings false alarms back to the level of the true model, which restores the guarantee when the forecast model is right. But a forecast that keeps adapting also absorbs most of a real level shift, and we derived by how much in closed form: a doubling of the long-run variance raises the ratio of squared returns to the forecast by only 9% at persistence 0.99, and no shift, however large, can push it above (1 − β)/α = 1.2 (E10) [SYNTHETIC].
* **Combining time scales needs care.** Averaging detectors run on 1, 5 and 15 minute returns has no proof if checked every minute, because a coarse detector is not a fair bet between its block ends. Counting it only at the wealth it is sure to keep gives a proved bound, Σ w_k ⌈m/k⌉ / A; the combination is never faster than the best single frequency on a persistent change, but catches 95% to 97% of a slow trend the one-minute detector almost never sees (E17) [SYNTHETIC].
* **Price-level microstructure noise breaks the guarantee on raw returns** (false alarms rise from 5% to 41.5% per day as noise grows), and non-overlapping pre-averaged blocks restore it (ALT1) [SYNTHETIC].
* **On one-minute crypto data, large returns come in bursts**, not as isolated jumps: 44% to 65% have a large neighbour, against 8.9% by chance. Our normalised version of your robust kernel then comes out at 0.33 to 0.42 times MedRQ. The true quarticity is unknown on real data, but if a calibrated burst model is right, any estimator that clips against a scale measured over hours is biased down by about 35% to 59%, and a local scale removes most of that (ALT3) [REAL and SYNTHETIC].
* **The honest bottom line on real data**: a simple ratio rule tuned with hindsight is about as precise as the e-detector (E6e). What the e-detector adds over that rule is a false-alarm guarantee with settings fixed in advance; the published change-point methods of E13 are less precise on crypto at the same alarm counts.

## 8. What are the limitations?

* The guarantee needs the clipped squared returns to keep their conditional mean while nothing changes. Clustering, price noise and estimated baselines violate this; each has a fix, and each fix costs speed.
* The bound is about four times conservative with the default set of bets (E9, E11).
* After an alarm the default monitor recalibrates and is blind about 70% of the time on real crypto data (E12). The awake monitors of E14 watch every release instead, but without a null band the rank version raises 3.6 to 4.3 alarms per crypto day.
* Our crypto and Nifty data come from public mirrors, not from the exchanges, and your paper does not name its data source.
* Full list: [`docs/limitations.md`](docs/limitations.md).

## 9. What should be done next?

1. A short joint paper on the monitor (a draft is in `paper/`): the bounded-statistic argument, the delay theory and the decrease speed limit, the confidence-sequence baseline, and the crypto and Nifty applications.
2. An estimation follow-up on your estimator: a local-scale plug-in version that reaches 95% efficiency in simulation (ALT3), which needs a central limit theorem.
3. A sequential-estimation bridge to your own field (ALT6, built): purely sequential fixed-accuracy estimation of the volatility level. Chow and Robbins on clipped squared returns, which runs on a running quarticity estimate, covers 94% to 96% at about the ideal sample size for every return law we tried, where raw squares fail under fat tails; an anytime-valid confidence sequence costs an exact factor 2 log(2/α) / z² = 1.92 more (1.55 with a running intersection); on BTC, about 50 minutes give plus or minus 20% on volatility. A natural next step is to use your robust kernel estimator as the variance plug-in and study its stopping time.
4. Finish the benchmark against published methods (E13): check our reconstruction of the Bibinger, Jirak and Vetter test against the paper's equations, and try BOCPD with a heavy-tailed likelihood, since BOCPD is faster than our monitor for decreases on clean data.

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
python -m pytest -q                                   # about 240 tests, no data needed
python -m avmonitor simulate --ratio 2                # synthetic stream with a known change
python -m avmonitor run --source binance --symbol BTCUSDT --download \
    --start 2025-10-08 --end 2025-10-14 --plot btc.png  # real data around the 10 October 2025 crash
```
