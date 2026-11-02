# Pitch for Prof. Sudeep Bapat

*A short brief for the first meeting. Plain language first, details on request.*

## The idea in three sentences

Your 2026 estimator clips each return at a multiple of the MAD and then smooths it locally with a kernel. That makes it a **robust, bounded, local** measure of volatility at every minute. We use exactly that property to build a **sequential monitor** that raises an alarm when volatility shifts up or down, with a false alarm guarantee that holds however often it is checked (an anytime-valid e-detector).

## Why it fits your work

* It extends your quarticity paper from *estimation over a fixed sample* to *online detection*. You list adaptive and data-driven extensions as future work, and the online version goes in that direction.
* The core tool is sequential analysis (test martingales and optional stopping), which is close to your work on purely sequential procedures.
* It gives a clean theoretical role to truncation. Without an upper bound on the monitored statistic, **no** nonparametric e-test can ever detect a fall in volatility (a two-line proof in `docs/math.md`, 5.3). Your truncation provides the bound.

## What the prototype does today

1. Implements your estimators exactly as printed, plus normalised versions and all the benchmarks in your Tables 3 and 4 (`src/avmonitor/robust.py`).
2. Runs a streaming monitor: online MAD, truncated variance, Shiryaev-Roberts betting e-detectors for increases and decreases, and recalibration after each alarm (`src/avmonitor/monitor.py`).
3. **Simulation [synthetic]:** with jumps or Student-t returns, Gaussian likelihood-ratio and CUSUM detectors raise a false alarm in a day on 70% to 97% of paths. The robust e-detector stays at 0.3% to 2.3%, within its 5% guarantee, and still detects a 50% volatility increase in a median of 36 to 49 minutes.
4. **Real data [Bitstamp BTC/USD, Jul to Dec 2025]:** 81% of its "volatility up" alarms are followed by a day of higher volatility, and 90% of its "down" alarms by lower volatility, against 50% at random times. Gaussian benchmark alarms are weaker (63%).

## Questions about the paper's code (asked to learn, not to correct)

While reproducing the tables we noticed a few things that we would like to understand. Each one is documented with evidence in `docs/estimator_issues.md`.

1. **Kernel weights.** Are they normalised to sum to one? With raw weights, the estimate scales with the weight sum (about 10.6 at h = 30). In Table 6, Kernel / RQ is 4.0 × 10⁻⁵ for all five coins, which is exactly the raw weight sum divided by n. We think this explains why RMSE rises with h in Tables 7 and 8.
2. **MAD scaling.** Is the MAD multiplied by 1.4826? Without it, c = 3 clips at about 2 standard deviations and keeps about 70% of the fourth moment, at every sample size. Our reproductions of Tables 6 and 7 fit the unscaled version better.
3. **Simulation settings.** What per-step volatility was used? One value, σ ≈ 0.0215, reproduces most of Table 7. Were the noise levels in Table 8 variances or standard deviations? And what differs in the σ_J = 0.10 design, where our robust-kernel numbers come out lower than yours?
4. **RQ in the tables.** Was it computed without the factor n? Tables 3, 4 and 6 fit `(1/3) Σ r⁴`.

If our reading is right, the fixes are one line each, and the robustness idea stands. Once the estimator is normalised, its most useful form is the local estimate, which is what the monitor uses.

## What we would like from you

* A quick check of the questions above, or a look at your code, so the comparison is fair.
* Your view on the target: should the monitor watch the variance stream (as now), the quarticity stream, or both?
* Data suggestions for Indian markets (Nifty 50 or stock one-minute bars), and whether the IITB-Citadel lab has access.
* Whether this could become a joint short paper: theory (guarantees and the impossibility result), simulations, and a Nifty and crypto application.

## One-minute demo

```bash
pip install -e ".[dev]"
python -m pytest -q                       # 26 tests
python experiments/e04_detection.py       # detection table (synthetic, about 1 minute)
python experiments/e05_real_btc.py        # real BTC run (needs the data file, see data/README.md)
```
