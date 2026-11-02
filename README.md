# Anytime-Valid Volatility Monitor

**An online monitor that flags volatility regime breaks in high-frequency returns, with false-alarm guarantees that hold no matter how often you check it.** It combines jump- and noise-robust realized measures with e-detectors, a nonparametric framework for sequential change detection.

> Status: working prototype (November 2026). Estimators, detectors, a streaming monitor, 26 tests, simulation studies and a first real-data run on Bitstamp BTC/USD are in place. Start with [`PITCH_BAPAT.md`](PITCH_BAPAT.md) for the short version.

## Key results so far

* **The detector keeps its promise when returns are messy [synthetic].** With jumps or fat-tailed returns, Gaussian likelihood-ratio and CUSUM detectors raise a false alarm within a day on 70% to 97% of paths. The robust e-detector stays at 0.3% to 2.3%, inside its 5% guarantee, and still detects a 50% volatility increase in a median of 36 to 49 steps ([experiment report](docs/experiment_report.md), E4).
* **Truncation is what makes "volatility went down" detectable.** Without an upper bound on the monitored statistic, no nonparametric e-test can detect a decrease ([math notes](docs/math.md), 5.3).
* **Alarms on real data mean something [real].** On Bitstamp BTC/USD one-minute bars (July to December 2025), 81% of increase alarms and 90% of decrease alarms are followed by a day of volatility in the flagged direction, against 50% at random times (E5).
* **Reproducing the paper.** The printed kernel formulas are not on the integrated-quarticity scale, and the c × MAD truncation leaves a bias that does not shrink with n. Both have one-line fixes. The published tables are consistent with the printed formulas ([estimator investigation](docs/estimator_issues.md)).

## Motivation

Volatility in high-frequency markets shifts between regimes. Risk limits, position sizing and model retraining should all react when it does. Two problems make this hard to do online:

1. **Jumps and microstructure noise.** Raw realized measures are distorted by a few extreme returns. Robust estimators fix this offline. For example, Bapat (2026) proposes a kernel-weighted estimator of quarticity that truncates returns at a multiple of the median absolute deviation (MAD).
2. **Repeated looks.** A classical test checked every minute will eventually raise a false alarm. Anytime-valid methods based on e-processes control false alarms at every stopping time (Shin, Ramdas and Rinaldo, 2024).

This project joins the two into one online monitor and adds a simple risk switch that acts when a change is detected. We are not aware of anytime-valid change detection applied to Indian markets (Nifty 50, NSE stocks, India VIX).

## Research questions

1. **Online robustness:** Does replacing the full-sample MAD with a rolling or exponentially weighted MAD keep the robustness of the kernel estimator, so that it can run live?
2. **Detection:** How do e-detectors built on robust local variance compare with CUSUM and with detectors on raw realized variance, in false alarm rate and detection delay?
3. **Application:** Does a change-triggered risk switch (for example, cutting exposure after a detected volatility increase) improve drawdowns and risk-adjusted returns out of sample?

## Method

| Step | What | Notes |
| --- | --- | --- |
| 1 | Robust local measures | Truncate returns at c times MAD, then kernel-weight squared returns (local variance) and fourth powers (local quarticity), following Bapat (2026). |
| 2 | Online truncation | Replace full-sample MAD with a rolling or exponentially weighted MAD, with a floor for assets whose prices move in discrete ticks. |
| 3 | Bet sizing | Local truncated quarticity is the variance of the bet, so it can size the bets (optional: it did not beat a fixed grid of bets in E4). |
| 4 | E-detector | Shiryaev-Roberts betting e-detectors on the truncated squared return, one for increases and one for decreases. Proven: average run length at least A, and P(false alarm within m steps) at most m / A. |
| 5 | Benchmarks | Gaussian likelihood-ratio e-detector, Page's CUSUM, a repeated window z-test, and the same e-detector on raw r². |
| 6 | Simulation | Jump diffusion with Bernoulli jumps and microstructure noise, with known change points, to measure false alarms and delay. |
| 7 | Risk switch | Scale exposure on detection and backtest on held-out data (not started). |

## Data

| Source | Use | Access |
| --- | --- | --- |
| Crypto 1-minute bars (BTC, ETH, XRP, SOL, DOGE) | Comparable with the empirical setting of Bapat (2026) | Free exchange APIs |
| Nifty 50 1-minute bars | Indian market application | Free or low cost vendors |
| India VIX daily | Check of detected regimes | Free from NSE |

Raw data is not committed. See [`data/README.md`](data/README.md).

## Repository layout

```
.
├── PITCH_BAPAT.md              short brief for Prof. Bapat
├── RESEARCH_STRATEGY_BAPAT.md  questions, decisions, roadmap, risks
├── ITERATION_LOG.md            what was tried and what changed
├── alternatives/               other directions, kept out of the main package
├── data/                       raw and processed data (git-ignored), download steps
├── docs/
│   ├── paper_review.md         the paper in plain words (17-point review)
│   ├── estimator_issues.md     checks of the estimator, with evidence
│   ├── math.md                 derivations and proofs (truncation, e-detector guarantees)
│   ├── literature_review.md    verified references and the gap
│   └── experiment_report.md    E1 to E5 results, labelled synthetic or real
├── experiments/                one script per research question
├── results/                    CSV outputs and figures
├── src/avmonitor/
│   ├── robust.py               all quarticity estimators (as printed and normalised), MAD truncation, online MAD
│   ├── detectors.py            betting e-detectors, Gaussian LR, Page CUSUM, repeated z-test
│   ├── monitor.py              streaming calibrate, monitor, alarm, recalibrate loop
│   ├── simulate.py             paper design, stochastic volatility with jumps and noise, change points
│   ├── data.py                 Bitstamp loader and intraday deseasonalisation
│   └── strategy.py             placeholder for the risk switch (not started)
└── tests/
```

## Getting started

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
python -m pytest -q              # 26 tests
python experiments/e04_detection.py
```

Minimal use of the monitor:

```python
import numpy as np
from avmonitor.monitor import VolatilityMonitor

r = np.concatenate([np.random.randn(3000), 2 * np.random.randn(1000)])  # sd doubles at step 3000
mon = VolatilityMonitor(threshold=7800, cal_len=1000, null_band=0.25)   # 5% false alarm budget per 390 steps
for e in mon.run(r):
    print(e.index, e.direction)
```

## Roadmap

- [x] Robust kernel estimator with full-sample MAD, checked against the simulation setup in Bapat (2026)
- [x] Rolling and exponentially weighted MAD versions
- [x] Simulator with known change points
- [x] E-detector and CUSUM benchmarks, compared on false alarms and detection delay
- [x] Crypto 1-minute application (Bitstamp BTC/USD)
- [ ] ETH, XRP, SOL, DOGE and Nifty 50
- [ ] Confidence-sequence baseline so the guarantee also covers an estimated baseline
- [ ] Change-triggered risk switch backtest
- [x] Short write-up of first results

## References

The full, verified list is in [`docs/literature_review.md`](docs/literature_review.md).


- Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482. doi:10.1016/j.iref.2026.105482
- Barndorff-Nielsen, O. E. and Shephard, N. (2004). Power and bipower variation with stochastic volatility and jumps. *Journal of Financial Econometrics*, 2(1), 1 to 37.
- Lee, S. S. and Mykland, P. A. (2008). Jumps in financial markets: a new nonparametric test and jump dynamics. *Review of Financial Studies*, 21(6), 2535 to 2563.
- Mancini, C. (2009). Non-parametric threshold estimation for models with stochastic diffusion coefficient and jumps. *Scandinavian Journal of Statistics*, 36, 270 to 296.
- Page, E. S. (1954). Continuous inspection schemes. *Biometrika*, 41(1/2), 100 to 115.
- Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *Journal of the Royal Statistical Society Series B*, 86(1), 1 to 27.
- Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science*, 2(2), 229 to 260. doi:10.51387/23-NEJSDS51

## Author

Shreyan Thanki ([@ShockXR](https://github.com/ShockXR))
