# Anytime-Valid Volatility Monitor

**An online monitor that watches one-minute returns and raises an alarm when volatility shifts up or down, with a proven limit on false alarms that holds however often you look at it.** It puts the jump-robust truncation of Bapat (2026) inside e-detectors, a nonparametric framework for sequential change detection (Shin, Ramdas and Rinaldo 2024).

> Status: research prototype, round 2 (November 2026). About 190 tests, 15 experiments on simulated and real data, five alternative directions, proofs in [`docs/math.md`](docs/math.md). For the short version read [`PITCH_BAPAT.md`](PITCH_BAPAT.md); for every number read the [experiment report](docs/experiment_report.md).

Every result in this repository carries a label: **[PUBLISHED]** (taken from a paper), **[SYNTHETIC]** (our simulation) or **[REAL]** (our run on market data). Nothing is reported that was not run.

## What the project does

1. **Checks the estimator.** It re-implements the robust kernel quarticity estimator of Bapat (2026) exactly as printed, finds what it converges to, and reproduces the published simulation and crypto tables.
2. **Turns the estimator's key idea into a monitor.** The paper clips each return at c times its median absolute deviation (MAD). We monitor the clipped squared return X_t = min(r_t², τ²) with a betting e-detector that alarms when its running wealth passes a threshold A.
3. **Proves and measures the false-alarm guarantee.** Average time to a false alarm is at least A steps, and the chance of any false alarm within m steps is at most m / A, under jumps, fat tails and any return distribution, as long as the conditions in [`docs/math.md`](docs/math.md) section 5 hold. We also show where those conditions fail on real data and what fixes them.
4. **Applies it to real markets:** Binance BTC, ETH, XRP, SOL and DOGE (the paper's five coins), Bitstamp BTC/USD and the Nifty 50 index, with event studies on FOMC statements, US CPI releases and Indian Union Budget days.

## Research question

*Can a robust, jump-resistant volatility measure drive an online detector of volatility changes whose false alarm rate is guaranteed in finite samples, without a model for the returns, and what does that guarantee cost in detection speed?*

Sub-questions (status in [`RESEARCH_STRATEGY_BAPAT.md`](RESEARCH_STRATEGY_BAPAT.md), section 3): correctness of the printed estimators (Q1); online truncation from past data only (Q2); validity and delay against classical and published detectors (Q3); why a bounded statistic is needed (Q4); estimated baselines (Q5); when the guarantee stops applying (Q6); meaning of alarms on real data (Q7); use in a risk rule (Q8).

## Key results

* **The guarantee holds where classical detectors break [SYNTHETIC].** With jumps or fat tails, Gaussian likelihood-ratio and CUSUM detectors raise a false alarm within a day on 70% to 97% of paths. The robust e-detector stays at 0.3% to 2.3% against its 5% bound, and detects a 50% volatility increase in a median of 36 to 49 steps (E4). Across four very different return laws the measured rate is 1.1% to 1.4% (E9), so the bound is valid but about four times conservative.
* **Why truncation matters for sequential testing.** No nonparametric e-test can detect a *fall* in the variance (a one-sided form of Bahadur and Savage 1956, [`docs/math.md`](docs/math.md) 5.3). A bounded functional can, and the truncated variance is one with a direct volatility meaning.
* **Estimated baselines need care [SYNTHETIC].** Plugging in an estimated "normal" level roughly triples false alarms with a 390-step calibration window and gives 44% against a guaranteed 25% over five days (E4e, E10). A confidence-sequence baseline restores a proof, P(false alarm within m) ≤ δ + m / A, at a cost of 1.8 to 7.2 times the delay with short calibration (E8).
* **Volatility clustering breaks the null, and fixing it caps the signal [SYNTHETIC].** Under persistent GARCH clustering the raw monitor alarms on 82% to 100% of quiet days. Dividing returns by a fitted GARCH forecast brings this back to 3% to 6%, but a forecast that keeps adapting passes at most (1 − β)/α of any level shift: a doubling of volatility looks like a 9% rise at persistence 0.99. We derived this cap in closed form and checked it by simulation (E10).
* **Alarms on real data carry information, but the main contribution is the guarantee [REAL].** On the five coins, alarms are followed by volatility on the flagged side of EWMA and HAR forecasts 71% to 83% of the time, against about 50% at random times, in and out of sample (E6). A simple ratio rule tuned with hindsight is about as precise (E6e), so what the e-detector adds is a false-alarm guarantee with fixed settings, not extra precision.
* **Many assets at once [SYNTHETIC and REAL].** A Shiryaev-Roberts statistic divided by time is not an e-value, so the usual multiple-testing recipes do not apply directly. Per-asset e-processes with e-BH give a proved false discovery rate at every stopping time (at most 0.96% observed in 48 designs). On five coins plus Nifty its discoveries are confirmed 78% to 97% of the time against about 50% at random times; its gain over a Bonferroni split is mostly earlier flags, not more of them (E15, `docs/theory_panel.md`).
* **Events [REAL].** 14 to 15 of 16 FOMC statements are flagged within six hours, against about 30% of matched control times (E12). On Nifty 50, Union Budget days are flagged (E7).
* **Reproducing the paper [PUBLISHED vs SYNTHETIC].** The published tables are consistent with the printed formulas, unnormalised kernel weights, the raw (unscaled) MAD and a crypto sample with about 2,000 missing minutes. Two findings we would like to discuss: the printed kernel estimators are not on the integrated-quarticity scale, and truncation at c × MAD leaves a bias that does not shrink with n. Both have one-line fixes ([`docs/estimator_issues.md`](docs/estimator_issues.md), E1, E6).

We do **not** claim to have improved the estimator. Once both kernels are put on the same scale, the robust kernel beats the plain kernel in most simulation cells (E1).

## Method in brief

| Step | What | Where |
| --- | --- | --- |
| 1 | Truncation level τ = c × MAD × 1.4826 from past data only (rolling or exponentially weighted), c = 3 by default; optionally a tighter c_down = 2 for the decrease side, which detects falls faster (E11f) | `robust.py`, `monitor.py` |
| 2 | Bounded statistic X_t = min(r_t², τ²), after an intraday profile estimated from past days only | `robust.py`, `data.py` |
| 3 | Baseline μ₀ = mean of X on a calibration window: plug-in, plug-in with a margin (band), or a confidence sequence | `detectors.py`, `confseq.py` |
| 4 | Betting e-detector (Shiryaev-Roberts form): wealth M_t = Σ_j Π_{i=j..t} (1 + λ (X_i − μ₀)) over start times j, averaged over a grid of bets λ; separate detectors for increases and decreases | `detectors.py` |
| 5 | Alarm when M_t ≥ A; A = (steps per day) / (false alarm budget per day); then recalibrate | `monitor.py` |
| 6 | Optional forecast-relative mode: divide returns by a GARCH or EWMA variance forecast first | `forecast.py` |

The proofs (Theorem A: average run length ≥ A; Theorem B: P(false alarm within m) ≤ m / A; Corollary C for restarts; Theorem 6.4 for the confidence-sequence baseline) and the exact conditions are in [`docs/math.md`](docs/math.md). [`docs/equation_map.md`](docs/equation_map.md) maps every symbol to the line of code that computes it.

## Data

| Source | Use | Access |
| --- | --- | --- |
| Binance spot 1-minute bars: BTC, ETH, XRP, SOL, DOGE | The paper's five coins (E6) | Public GitHub mirror (exchange APIs are blocked from our cloud machine) |
| Bitstamp BTC/USD 1-minute bars | First real-data run (E5) | Public GitHub mirror |
| Nifty 50 1-minute bars, April 2017 to October 2026 | Indian market application (E7) | Public GitHub mirror |
| India VIX daily | Check of Nifty alarms (E7) | NSE archive server |

Raw data is not committed. `python scripts/download_data.py` fetches everything and writes `data/MANIFEST.csv` with a SHA-256 hash of each file. Details and caveats: [`data/README.md`](data/README.md).

## Install

Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
python -m pytest -q              # about 190 tests, no data needed
```

## Run

```bash
python -m avmonitor simulate --ratio 2                       # synthetic stream, sd doubles at step 3000
python -m avmonitor estimators                               # all quarticity estimators on one simulated day
python -m avmonitor run --source binance --symbol ETHUSDT --download \
    --start 2025-03-01 --end 2025-03-31 --plot eth.png       # real data, alarms printed and plotted
```

`notebooks/demo.ipynb` is an executed walk-through (rebuild it with `python scripts/build_demo_notebook.py`).

Minimal use from Python:

```python
import numpy as np
from avmonitor.monitor import VolatilityMonitor

rng = np.random.default_rng(0)
r = np.concatenate([rng.standard_normal(3000), 2 * rng.standard_normal(1000)])  # sd doubles at step 3000
mon = VolatilityMonitor(threshold=7800, cal_len=1000, null_band=0.25)          # 5% false-alarm budget per 390 steps
for e in mon.run(r):
    print(e.index, e.direction)
```

## Reproduce the experiments

Each experiment is one script in `experiments/` that writes CSV files to `results/` and figures to `results/figures/`. Synthetic experiments need no data; real-data experiments need `python scripts/download_data.py` first. Seeds are fixed.

| Experiment | Question | Script | Label |
| --- | --- | --- | --- |
| E1 | Reproduce the paper's simulation tables | `e01_reproduce_tables.py` | SYNTHETIC vs PUBLISHED |
| E2 | Truncation bias, fair comparison of estimators | `e02_truncation_and_fair_comparison.py` | SYNTHETIC |
| E3 | Online (past-only) truncation | `e03_online_truncation.py` | SYNTHETIC |
| E4 | False alarms and delay against benchmarks; band, calibration length, plug-in checks | `e04_detection.py`, `e04c_null_band.py`, `e04d_calibration_length.py`, `e04e_review_checks.py` | SYNTHETIC |
| E5 | Bitstamp BTC/USD run | `e05_real_btc.py` | REAL |
| E6 | Five coins: paper tables, monitor, forecasts, rate-matched competitor | `e06a_five_coins_tables.py`, `e06b_five_coins_monitor.py` | REAL vs PUBLISHED |
| E7 | Nifty 50, India VIX, Union Budget and policy days | `e07_nifty.py`, `e07c_india_vix.py`, `e07d_events.py` | REAL |
| E8 | Confidence-sequence baseline | `e08_cs_baseline.py` | SYNTHETIC and REAL |
| E9 | Stress tests of the guarantee, runtime | `e09_stress_guarantee.py` | SYNTHETIC |
| E10 | Clustering and forecast-relative monitoring | `e10_forecast_relative.py` | SYNTHETIC and REAL |
| E11 | Detection delay theory and bet design; a tighter clip for decreases | `e11_delay_theory.py`, `e11f_split_c.py` | SYNTHETIC and REAL |
| E12 | FOMC and CPI event study on crypto | `e12_crypto_events.py` | REAL |
| E13 | Benchmark against published change-point methods | `e13_benchmarks.py` | SYNTHETIC |
| E14 | Running truncation level and monitors without a blind period | `e14_running_tau.py` | SYNTHETIC |
| E15 | Panel of assets with false discovery control | `e15_panel.py synthetic`, `horizon`, `real`, `cs`, `figures` | SYNTHETIC and REAL |

Then rebuild the report from its parts:

```bash
python scripts/build_report.py   # writes docs/experiment_report.md from docs/report_parts/
python experiments/make_figures.py
```

Most scripts take part names as arguments (for example `python experiments/e10_forecast_relative.py a b`) so that one part can be rerun alone. The largest simulations take tens of minutes on one core.

## How to read the output

* **An alarm** means the running evidence against "volatility is at its calibrated level" has reached A. The direction says whether volatility moved up or down. The guarantee is about false alarms: with budget 0.05 per day, a monitor on a stream with no change raises a false alarm on at most about 5% of days in the long run.
* **The guarantee is conditional.** It needs the clipped squared returns to have the calibrated conditional mean while nothing changes. Strong volatility clustering and price-level noise violate this; the forecast-relative mode and non-overlapping blocks restore it (E10, ALT1).
* **Delay** is measured in steps (minutes). It grows like ln A divided by an information rate that we compute in closed form (E10, E11).
* **Real-data hit rates** compare the next day's realised volatility with a forecast. They are compared with random times and with a rate-matched rule, never read alone.

## Known limitations

The full list, with what each one affects, is in [`docs/limitations.md`](docs/limitations.md). The main ones:

* The guarantee is proved for a known baseline and for the confidence-sequence baseline. The default plug-in baseline is only approximately valid (E4e).
* Volatility clustering is itself a change in the conditional variance; without a forecast the monitor alarms on it. With a forecast it is valid but much slower to see level shifts (E10).
* The monitor is blind for about 70% of the time on real data while it recalibrates after alarms (E12).
* The Binance and Nifty data come from public mirrors, not from the exchanges. The paper does not name its crypto data source, so our five-coin data may differ from the paper's.
* The bound is about four times conservative (E9).

## Repository layout

```
.
├── PITCH_BAPAT.md              short brief: ten questions answered
├── RESEARCH_STRATEGY_BAPAT.md  questions, hypotheses, scores, risks, roadmap
├── ITERATION_LOG.md            iterations in a fixed format: problem, hypothesis, change, result
├── alternatives/               five other directions, each with its own README and RESULTS
├── data/                       raw data (git-ignored), MANIFEST.csv with hashes, README
├── docs/
│   ├── paper_review.md         the paper in plain words
│   ├── estimator_issues.md     checks of the estimator, with evidence
│   ├── math.md                 derivations and proofs
│   ├── theory_delay.md         detection delay and bet design
│   ├── theory_panel.md         many assets at once: e-BH and false discovery control
│   ├── equation_map.md         symbol to code map
│   ├── literature_review.md    verified references and the gap
│   ├── limitations.md          what the results do not show
│   ├── experiment_report.md    all results (generated from report_parts/)
│   ├── report_parts/           one file per experiment
│   ├── iterations/             per-workstream iteration notes
│   └── review/                 independent reviews and responses
├── experiments/                one script per experiment
├── notebooks/demo.ipynb        executed demonstration
├── results/                    CSV outputs and figures
├── scripts/                    data download, report build, notebook build
├── src/avmonitor/
│   ├── robust.py               quarticity estimators (as printed and normalised), MAD truncation, online MAD
│   ├── detectors.py            betting e-detectors, Gaussian LR, Page CUSUM, repeated z-test
│   ├── confseq.py              confidence-sequence baseline
│   ├── monitor.py              streaming calibrate, monitor, alarm, recalibrate loop
│   ├── forecast.py             EWMA and GARCH variance forecasts for forecast-relative monitoring
│   ├── competitors.py          published change-point methods used as benchmarks
│   ├── running.py              running truncation level
│   ├── design.py               bet grid design from the delay theory
│   ├── panel.py                panel monitor with e-BH across assets
│   ├── simulate.py             paper design, stochastic volatility with jumps and noise, GARCH, change points
│   ├── data.py                 loaders, predictable intraday profile, frozen-price detection
│   ├── strategy.py             interface for the risk switch (implemented in alternatives/risk_switch)
│   └── cli.py                  command line
└── tests/
```

## Roadmap

- [x] Robust kernel estimators, as printed and normalised, checked against the paper's simulation tables
- [x] Online (past-only) truncation
- [x] E-detector with proofs, compared with Gaussian, CUSUM and repeated-test benchmarks
- [x] Real data: Bitstamp BTC, the five Binance coins, Nifty 50, event studies
- [x] Confidence-sequence baseline with a proof
- [x] Forecast-relative monitoring under clustering, with the signal cap
- [x] Risk switch backtest (alternatives/risk_switch)
- [ ] Warm restart after alarms to remove the blind period
- [ ] Short paper draft
- [ ] Discussion with Prof. Bapat of the open questions in [`docs/estimator_issues.md`](docs/estimator_issues.md)

## References

The full, verified list is in [`docs/literature_review.md`](docs/literature_review.md).

- Bahadur, R. R. and Savage, L. J. (1956). The nonexistence of certain statistical procedures in nonparametric problems. *Annals of Mathematical Statistics*, 27(4), 1115 to 1122.
- Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482. doi:10.1016/j.iref.2026.105482
- Barndorff-Nielsen, O. E. and Shephard, N. (2004). Power and bipower variation with stochastic volatility and jumps. *Journal of Financial Econometrics*, 2(1), 1 to 37.
- Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. *Journal of Econometrics*, 31(3), 307 to 327. doi:10.1016/0304-4076(86)90063-1
- Lai, T. L. (1998). Information bounds and quick detection of parameter changes in stochastic systems. *IEEE Transactions on Information Theory*, 44(7), 2917 to 2929.
- Lamoureux, C. G. and Lastrapes, W. D. (1990). Persistence in variance, structural change, and the GARCH model. *Journal of Business and Economic Statistics*, 8(2), 225 to 234.
- Lee, S. S. and Mykland, P. A. (2008). Jumps in financial markets: a new nonparametric test and jump dynamics. *Review of Financial Studies*, 21(6), 2535 to 2563.
- Mancini, C. (2009). Non-parametric threshold estimation for models with stochastic diffusion coefficient and jumps. *Scandinavian Journal of Statistics*, 36, 270 to 296.
- Mikosch, T. and Starica, C. (2004). Nonstationarities in financial time series, the long-range dependence, and the IGARCH effects. *Review of Economics and Statistics*, 86(1), 378 to 390.
- Page, E. S. (1954). Continuous inspection schemes. *Biometrika*, 41(1/2), 100 to 115.
- Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science*, 2(2), 229 to 260. doi:10.51387/23-NEJSDS51
- Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *Journal of the Royal Statistical Society Series B*, 86(1), 1 to 27.

## Author

Shreyan Thanki ([@ShockXR](https://github.com/ShockXR))
