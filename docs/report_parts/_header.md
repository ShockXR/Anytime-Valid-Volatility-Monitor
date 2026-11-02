# Experiment report

Every experiment answers one research question. Scripts are in `experiments/` (and in `alternatives/<name>/` for the alternative prototypes), raw outputs in `results/`, figures in `results/figures/`. Seeds are fixed in the scripts. Each result is labelled:

* **[SYNTHETIC]** for our simulations.
* **[REAL]** for market data: Bitstamp BTC/USD one-minute bars (E5), Binance spot one-minute bars for BTC, ETH, XRP, SOL and DOGE (E6, E10, E12 and later), Nifty 50 one-minute bars and India VIX daily closes (E7). Sources and SHA-256 hashes are in `data/MANIFEST.csv`; raw files are not committed (`data/README.md`).
* **[PUBLISHED]** for numbers copied from Bapat (2026).

To rerun: `pip install -r requirements.txt && pip install -e .`, then `python scripts/download_data.py`, then `python experiments/<script>.py`. This report is assembled from `docs/report_parts/` by `python scripts/build_report.py`. What the results do not establish is collected in `docs/limitations.md`.

## Index

| ID | Question | Main finding |
|---|---|---|
| E1 | Do the printed formulas reproduce the published simulation tables? | Yes for the plain kernel with one fitted volatility; the large-jump robust rows do not reproduce. After normalisation the robust kernel wins in most cells. |
| E2 | How large is the truncation bias? | 3 raw MADs keep 70.5% of the fourth moment, 3 scaled MADs 96.4%, at every n. |
| E3 | Can the truncation run online? | Yes: rolling or exponentially weighted MADs keep the jump protection. |
| E4 | False alarms and delay against parametric and classical detectors | The robust e-detector is the only valid one under jumps and fat tails, 10% to 70% slower on clean Gaussian data. The plug-in baseline roughly triples the false alarm rate with short calibration (E4e). |
| E5 | Real data, Bitstamp BTC | 81% and 90% of up and down alarms are followed by volatility in the flagged direction, but clustering explains much of the up figure. |
| E6 | The paper's five coins on Binance | Published tables imply unnormalised weights, the raw MAD and about 2,000 missing minutes; an ETH sign slip and a DOGE contradiction. Alarms beat EWMA and HAR forecasts, but a hindsight-tuned ratio rule is as precise. |
| E7 | Nifty 50 and India VIX | Out of sample, 78% of up alarms are followed by higher volatility (47% at random); every budget day watched from the open raised an alarm; the monitor answers "has volatility changed?", not "is it high?". India VIX moves with the alarms on the same day but the alarms do not lead it. |
| E8 | A guarantee with an estimated baseline | Proved: P(false alarm within m steps) ≤ δ + m / A with a confidence-sequence baseline; observed at most 0.27% against 5%. |
| E9 | How tight is the guarantee, and when does it fail? | About 4 times conservative and distribution free; persistent GARCH clustering breaks it. |
| E10 | Monitoring relative to a volatility forecast | Restores the guarantee under clustering; a frozen persistent forecast lets through at most (1 - β)/α of any break, and an information bound puts the delay near two days. |
| E11 | Detection delay and bet design (theory) | Proved delay bound and asymptotic optimality (the delay attains log A / KL_inf); the price of robustness is a ratio of information numbers that predicts E4's 14% to 73% penalties, mostly the cost of making no distributional assumption; truncation caps decrease detection at log(B / (B - 1)) per step (at least 73 steps at c = 3); the false alarm bound is nearly attained by an adversarial null for a single bet. |
| E12 | Crypto alarms around FOMC and CPI releases | Up alarms within 6 hours of 14 to 15 of 16 FOMC statements against about 30% of matched control windows; the streaming monitor is blind about 70% of the time while recalibrating. |
| E13 | Benchmark against Bibinger, Jirak and Vetter (2017) and Bayesian online change points | In progress. |
| E14 | A guarantee with a running truncation level | In progress. |
| E15 | Panel monitoring with false discovery control | Proved FDR control at every stopping time with e-BH on per-asset e-processes (a Shiryaev-Roberts statistic divided by t is not an e-value); FDR at most 0.96% in 48 designs; on five coins plus Nifty, discoveries are confirmed 78% to 97% of the time against about 50% at random, and e-BH's gain over Bonferroni is mostly earlier flags. |
| ALT1 | Quarticity stream; pre-averaging for noise | Price noise breaks the raw-return guarantee; non-overlapping pre-averaged blocks restore it. |
| ALT2 | Risk switch | Tails on par with daily volatility targeting, worse tracking, no Sharpe difference; the robust trigger beats the Gaussian one. |
| ALT3 | Robust local IQ estimator | Squaring a local variance is 95% efficient against 74% for RQ-type estimators; one daily MAD stops converging when volatility moves within the day; on BTC and DOGE large returns come in bursts and scale-truncated estimators remove a third to a half of IQ. |
| ALT4 | Correlation monitor | Sign concordance gives a valid correlation monitor that volatility changes cannot fool. |
| ALT5 | Two-sample (rank) betting test | Ranks against a growing reference give an exactly valid detector for any return law and horizon (1.9% to 3.5% per day against 5%); ranks against a fixed reference fail like a plug-in baseline; 2.9 times faster for a change in tail shape, 4% to 19% slower for a scale change. |
| ALT6 | Sequential fixed-accuracy estimation of volatility (classical purely sequential rules against anytime-valid confidence sequences) | In progress. |
