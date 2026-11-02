# Anytime-Valid Volatility Monitor

**An online monitor that flags volatility regime breaks in high-frequency returns, with false-alarm guarantees that hold no matter how often you check it.** It combines jump- and noise-robust realized measures with e-detectors, a nonparametric framework for sequential change detection.

> Status: early prototype. The research design is set, and the code is being built.

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
| 3 | Standardisation | Use local quarticity to scale local variance, since quarticity sets the asymptotic variance of realized variance. |
| 4 | E-detector | Build e-processes for an increase or decrease in variance and raise an alarm when the e-detector crosses its threshold. |
| 5 | Benchmarks | CUSUM, Page-Hinkley and detectors on raw realized variance. |
| 6 | Simulation | Jump diffusion with Bernoulli jumps and microstructure noise, with known change points, to measure false alarms and delay. |
| 7 | Risk switch | Scale exposure on detection and backtest on held-out data. |

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
├── data/            raw and processed data (git-ignored)
├── docs/            notes and write-ups
├── notebooks/       exploratory analysis
├── src/avmonitor/   package code
│   ├── data.py          loaders for 1-minute bars
│   ├── robust.py        MAD truncation and kernel-weighted local variance and quarticity
│   ├── detectors.py     e-detectors and CUSUM benchmarks
│   ├── simulate.py      jump diffusion with noise and known change points
│   └── strategy.py      change-triggered risk switch and backtest
└── tests/
```

## Getting started

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

## Roadmap

- [ ] Robust kernel estimator with full-sample MAD, checked against the simulation setup in Bapat (2026)
- [ ] Rolling and exponentially weighted MAD versions
- [ ] Simulator with known change points
- [ ] E-detector and CUSUM benchmarks, compared on false alarms and detection delay
- [ ] Crypto 1-minute application, then Nifty 50
- [ ] Change-triggered risk switch backtest
- [ ] Short write-up of first results

## References

- Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482. doi:10.1016/j.iref.2026.105482
- Barndorff-Nielsen, O. E. and Shephard, N. (2004). Power and bipower variation with stochastic volatility and jumps. *Journal of Financial Econometrics*, 2(1), 1 to 37.
- Lee, S. S. and Mykland, P. A. (2008). Jumps in financial markets: a new nonparametric test and jump dynamics. *Review of Financial Studies*, 21(6), 2535 to 2563.
- Page, E. S. (1954). Continuous inspection schemes. *Biometrika*, 41(1/2), 100 to 115.
- Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science*, 2(2), 229 to 260. doi:10.51387/23-NEJSDS51

## Author

Shreyan Thanki ([@ShockXR](https://github.com/ShockXR))
