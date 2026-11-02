# Pitch for Prof. Sudeep Bapat

*A short brief for the first meeting. Plain language first, details on request.*

## The idea in three sentences

Your 2026 estimator clips each return at a multiple of the MAD and then smooths it locally with a kernel. The clipping step (eq. 11) makes every return's contribution **bounded**. We use that step to build a **sequential monitor** that raises an alarm when volatility shifts up or down, with a false alarm guarantee that holds however often it is checked (an anytime-valid e-detector). Your local kernel quarticity can size the bets, but it did not beat a fixed grid of bets in our tests, so the default monitor uses the truncation step only.

## Why it fits your work

* It extends your quarticity paper from *estimation over a fixed sample* to *online detection*. You list adaptive and data-driven extensions as future work, and the online version goes in that direction.
* The core tool is sequential analysis (test martingales and optional stopping), which is close to your work on purely sequential procedures.
* It gives truncation a clear role in sequential testing. A fall in the variance itself cannot be detected by any nonparametric e-test (`docs/math.md` 5.3, a one-sided form of the impossibility result of Bahadur and Savage, 1956). A detector must therefore monitor a bounded functional of the returns. Truncation at c × MAD gives one with a direct volatility meaning: if the shape of the innovations is fixed, the truncated variance falls exactly when σ falls.

## What the prototype does today

1. Implements your estimators exactly as printed, plus normalised versions and all the benchmarks in your Tables 3 and 4 (`src/avmonitor/robust.py`).
2. Runs a streaming monitor: a MAD and baseline estimated from a calibration window (re-estimated after each alarm), truncated variance, and Shiryaev-Roberts betting e-detectors for increases and decreases (`src/avmonitor/monitor.py`). The guarantee is proved for a known baseline and for a confidence-sequence baseline; with the default plug-in baseline it holds approximately.
3. **Simulation [synthetic]:** with jumps or Student-t returns, Gaussian likelihood-ratio and CUSUM detectors raise a false alarm in a day on 70% to 97% of paths. The robust e-detector stays at 0.3% to 2.3%, within its 5% guarantee, and still detects a 50% volatility increase in a median of 36 to 49 steps.
4. **Real data [Bitstamp BTC/USD, Jul to Dec 2025]:** 81% of its "volatility up" alarms are followed by a day of higher volatility, and 90% of its "down" alarms by lower volatility, against 50% at random times. Gaussian benchmark alarms are weaker (63%). Volatility clustering explains much of the increase figure (a naive last-hour rule scores about the same); the decrease alarms are more informative than such a rule.

## Questions about the paper's code (asked to learn, not to correct)

While reproducing the tables we noticed a few things that we would like to understand. Each one is documented with evidence in `docs/estimator_issues.md`.

1. **Kernel weights.** Are they normalised to sum to one? With raw weights, the estimate scales with the weight sum (about 10.6 at h = 30). In Table 6, Kernel / RQ is 4.0 × 10⁻⁵ for all five coins, which is within 1% of the raw weight sum divided by n. We think this explains why RMSE rises with h in Tables 7 and 8.
2. **MAD scaling.** Is the MAD multiplied by 1.4826? Without it, c = 3 clips at about 2 standard deviations and keeps about 70% of the fourth moment, at every sample size. Our reproductions of Tables 6 and 7 fit the unscaled version better.
3. **Simulation settings.** What per-step volatility was used? One value, σ ≈ 0.0215, reproduces most of Table 7. Were the noise levels in Table 8 variances or standard deviations? And what differs in the σ_J = 0.10 design, where our robust-kernel numbers come out lower than yours?
4. **RQ in the tables.** Was it computed without the factor n? Tables 3, 4 and 6 fit `(1/3) Σ r⁴`.

If our reading is right, the fixes are one line each, and the robustness idea stands. Once the estimator is normalised, its most useful form is the local estimate, which is what the monitor uses.

## What we would like from you

* A quick check of the questions above, or a look at your code, so the comparison is fair.
* Your view on the target: should the monitor watch the variance stream (as now), the quarticity stream, or both?
* Data suggestions for Indian markets (Nifty 50 or stock one-minute bars), and whether the IITB-Citadel lab has access.
* Whether this could become a joint short paper. The tools are known: Shiryaev-Roberts e-detectors and their average-run-length guarantee (Shin, Ramdas and Rinaldo 2024), Ville's inequality, and the Bahadur and Savage impossibility. What would be new: a truncated-variance null with a matching decrease detector, the audit of the estimator's scale, the guarantee with an estimated baseline, and the crypto and Nifty applications.

## One-minute demo

```bash
pip install -e ".[dev]"
python -m pytest -q                       # all tests (about 2 minutes)
python experiments/e04_detection.py       # detection table (synthetic, about 1 minute)
python experiments/e05_real_btc.py        # real BTC run (needs the data file, see data/README.md)
```
