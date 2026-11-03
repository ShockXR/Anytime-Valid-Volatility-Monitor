# Defence summary for Prof. Bapat: ten questions, plain answers

Prepared 3 November 2026 by validation agent 2 (an AI agent). Every number is PUBLISHED (printed in your 2026 paper), SYNTHETIC (simulated), REAL (market data) or OUR CHECK (recomputed during validation), and every one agrees with the project's latest verified results. Files named are in this repository.

**Five-line version.** (1) We turned the clipping idea of your quarticity estimator (clip each return at c times its median absolute deviation) into an online alarm for volatility rises and falls. (2) The alarm carries a proved cap on false alarms: within the first m steps after calibration, at most m / A, with A = 7,800 giving 5% for a 390-step day. (3) The cap holds in simulation with large margin, but only for the first day after each calibration, and only if the clipped squared returns keep their conditional mean. (4) On real data the alarms beat random times but add little to a strong forecast, and the one confirmatory event test on fresh data failed. (5) No human has reviewed the proofs; all reviews were by AI agents.

## 1. How do you know this is correct?

Three layers. Proofs: Theorems A and B (average run length at least A; false alarm within m steps at most m / A) rest on a supermartingale argument and Ville's inequality (`docs/math.md` 5.4); I read it and found no gap, and an AI referee in review round 1 reported none, but no human has checked. Code: I checked 16 formulas and rules against the code with 20 independent numeric tests and all pass (OUR CHECK, `reports/checks/equation_vs_code_check.py`): the Shiryaev-Roberts and CUSUM recursions equal brute-force sums over start times; the clipped-moment factors equal 4 million Monte Carlo draws and your table to three decimals (raw MAD at c = 3 keeps 0.925 of the second moment and 0.705 of the fourth); the decrease speed limit of 58 steps reproduces. Numbers: a script traces 2,552 numbers in the short paper, 282 in the README and 199 in the pitch to committed sources; 581 tests pass. The tracer can match coincidentally, so I do not lean on it alone.

## 2. Why this evaluation?

Because there is no ground truth for "a real volatility change", the main yardstick is the quantity the theory controls: false alarms against the proved bound, plus detection delay against Gaussian likelihood-ratio and CUSUM detectors. For real data we cannot show when the null holds, so we permute the one-minute returns of calm days (keeping the real ticks, zeros and fat tails but making them exchangeable) and ask whether the proved monitors stay under their bounds. Alarm precision is measured as "followed by volatility on the flagged side of a forecast", against about 50% at random times. We do not use Sharpe ratios; this is not a trading system.

## 3. Could it be overfitting?

There is little to fit: the monitor has few knobs (c = 3, a grid of bets, A fixed by the bound). The worry is rather selection of reported results. Defences: the main real-data tests were pre-registered in git before the results were read (plans exist for E18 to E20, E22 to E25 and the live Nifty run; there is no external registry); the primary test was fixed in advance (BTC on FOMC) with Holm correction over the family; the best competitor was "picked after seeing the results" and the pitch says so; and the fresh-data test (E23, 16 releases of January to October 2026 that no analysis had seen) was exactly the guard against this, and it failed. Several analyses (E24 guard, E25 anytime-valid test) were added after seeing the failure; they are labelled post hoc, and a forward test to 1 June 2027 is registered instead of rescuing the claim.

## 4. What assumptions matter most?

The null is about the conditional mean of X_t = min(r_t^2, tau^2): it must stay at or below the baseline for a rise alarm (at or above for a fall alarm). Three things break it: volatility clustering (raw monitor alarms on 82% to 100% of quiet days, SYNTHETIC), price-level microstructure noise (false alarms rise from 5% to 41.5% per day), and an estimated baseline (plug-in baseline gives 41% to 44% over five days against a 25% bound). Also: the tau, bets and intraday profile must use only past data; one-minute bars are assumed; the proved baseline must be known or come from the confidence-sequence construction.

## 5. What if conditions are different?

Other return laws: four very different laws (Gaussian, Student t3 with infinite fourth moment, jumps, and coarse price ticks with 20% zero returns) gave 1.1% to 1.4% within the first day (E9, SYNTHETIC). Other days: the bound is for the first day after each calibration only; a stream of i.i.d. returns with |r| exactly sqrt(m0) is in the null and makes the default monitor alarm with certainty at step 2A, which is the last step of day 40 (reproduced, OUR CHECK; `docs/theory_per_period.md`). A restart at every day boundary with a known baseline gives a proved bound for every day at a cost of 2% to 14% in delay for 1,440-minute days. Other assets and frequencies: five coins, Bitstamp BTC and Nifty 50 at one minute were tried; tick data and other markets were not.

## 6. What is the strongest evidence?

The proof plus its stress tests. Inside its conditions the measured false-alarm rate is far below the bound: 0.9% in a fresh 6,000-path simulation against 5% (OUR CHECK), 0.3% to 2.3% under jumps and fat tails (E4), and 2.62% for the worst legal adversary we could construct against 5% (E16, SYNTHETIC). On real exchangeable calm days the proved monitors stayed under every bound (2.26%, 0.01%, 3.38% against 5%; E19, REAL). The delay price of robustness is explained: 14%, 35% and 73% slower than a Gaussian detector for volatility rises of 25%, 50% and 100% on clean normal data, predicted by a ratio of information numbers (E11).

## 7. What would invalidate it?

- A gap in a proof (no human check yet; the first thing a statistician should read is `docs/math.md` section 5.4 and section 6).
- A real data stream whose clipped squared returns have a drifting conditional mean while the monitor assumes none: block-resampled calm days alarm on 86% to 100% of streams, so the minute-level null is rarely true on real data. The proved guarantee is then about the monitor, not the data.
- The forward event test: if releases keep failing to separate from controls (fixed date 1 June 2027, power 0.82), the FOMC story is dead; it is already exploratory.
- A bug in the bets: factors must stay positive (increase lambda below 1; decrease lambda below 1 / (B - 1)); the code enforces both and I checked both extremes.

## 8. Is it reproducible?

Yes in design, partly confirmed. Seeds are fixed; raw data are not committed but `scripts/download_data.py` fetches pinned public mirrors and writes SHA-256 hashes; `scripts/reproduce.py` checks hashes, runs the tests and reruns experiments against `results/MANIFEST.csv`. The README says the quick tier and E1 to E17 reproduced on a fresh clone and E18 onward has not been rerun on a clean clone yet. I rebuilt the experiment report from its parts and traced the paper, README and pitch; another agent checks clean-clone reproduction (`reports/REPRODUCIBILITY.md`).

## 9. What is the most serious unresolved weakness?

Real data cannot show the guarantee is relevant: the null rarely holds at the minute scale, and the one confirmatory real-data test failed. The honest claim is "a monitor with a proved cap under stated conditions", not "a better detector of real regime changes". A classical CUSUM on the same clipped statistic is about as precise once its threshold is calibrated; the edge is a calibration-free threshold. Second: after each alarm the default monitor is blind while it recalibrates (about 70% of the time on real crypto data). Third: no human has reviewed the proofs.

## 10. Questions for you (tailored)

**"What is actually new, and do you claim to improve my estimator?"** No claim to improve it. The tools are known (Ville's inequality, e-detectors, confidence sequences, Bahadur and Savage). The new part is their combination with your clipping for jump-robust online volatility monitoring, plus specific results: the proved decrease speed limit, the daily-restart construction, the plug-in baseline failure and the confidence-sequence fix. Our literature search found no method with this combination, but it was not exhaustive, and some related papers were read only as abstracts. Once the two kernels are on the same scale your robust kernel beats the plain kernel in most simulation cells (E1).

**Questions we would ask you** (to learn, not to correct; each is documented in `docs/estimator_issues.md`):
1. Are the kernel weights normalised to sum to one? Your Table 6 fits the raw weight sum (kernel over RQ about 4.0 x 10^-5 for all five coins).
2. Is the MAD multiplied by 1.4826? The tables fit the unscaled MAD, which at c = 3 clips at about 2 standard deviations.
3. Which per-step volatility did the simulations use (0.0215 fits most of Table 7), and are the Table 8 noise levels variances or standard deviations?
4. How was the RMSE in Table 6 computed? It matches |estimate - 2.2 x 10^-6| for every coin.
5. DOGE: Table 5 has Q1 = Q3 = 0, which makes the MAD 0 and the clipped kernel exactly 0, yet Table 6 prints 2.00 x 10^-12, the same value as XRP. Could the two rows come from different series?
6. Was the MAD computed per day, over the whole sample, or locally? And what data source did the crypto application use?

## Things not to say

Do not say the guarantee is "5% per day" on every day. Do not say the plug-in baseline is proved. Do not present the FOMC result as established; it failed on fresh data (2 of 16 against 5.9%, p = 0.25) and the follow-ups are post hoc. Do not say the reviews were independent experts; they were AI agents. Do not call the live Nifty run of 3 November evidence: no alarm, because the monitor was still recalibrating after the previous day's alarm.
