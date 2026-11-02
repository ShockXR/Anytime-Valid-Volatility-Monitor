# Quarticity-stream e-detector (alternative A3)

A minimal working prototype of the idea in `alternatives/ALTERNATIVES.md`, entry A3: run the same anytime-valid betting e-detector as the main package, but on the truncated **fourth** power of the return, `min(r^4, tau^4)`, instead of the truncated square `min(r^2, tau^2)`. Results are in `RESULTS.md`; every number there comes from `run_experiment.py` and is labelled [SYNTHETIC] or [REAL].

## 1. Research problem

The main monitor watches the truncated variance of returns and alarms when it moves. Volatility is not the only thing that changes in a market. The *shape* of the return distribution can change while the variance stays put: tails can get fatter (more frequent large moves of either sign) or volatility can start to flicker from minute to minute (a higher "volatility of volatility"). Both show up in the fourth moment of returns, which is what quarticity measures. Bapat (2026) built a robust estimator of exactly this quantity. The question here is: can the same bounded, truncated statistic drive a sequential detector of changes in quarticity, with the same false alarm guarantee, and what does it see that the variance stream does not?

Concretely:

1. Is the quarticity-stream detector valid (same guarantee as the variance stream)?
2. Under a pure variance change, how much slower (or faster) is it than the variance stream?
3. Under a change that keeps the variance fixed but changes the tails or the vol of vol, does it detect what the variance stream cannot?
4. On real one-minute crypto data, how often does it alarm, and what do its alarms coincide with?

## 2. Mathematical formulation

Returns `r_t` arrive one at a time. A calibration window gives a robust scale (the scaled median absolute deviation, MAD, a standard-deviation estimate that ignores outliers) and a clipping level `tau = c x MAD` with `c = 3`. Define the two bounded streams

* variance stream: `X_t = min(r_t^2, tau^2)` in `[0, tau^2]`,
* quarticity stream: `Q_t = min(r_t^4, tau^4) = X_t^2` in `[0, tau^4]`.

The mean of `Q_t` is the truncated local quarticity of Bapat (2026), eq. (12), for a window of one return.

**Null hypotheses.** With `q0` the baseline mean of `Q_t`, the increase null is `H0+: E[Q_t | past] <= q0` for all t, and the decrease null is `H0-: E[Q_t | past] >= q0`. These are nonparametric: they restrict one conditional moment and allow any dependence, any shape and any jump process that keeps the moment on the null side.

**Null mean.** The model-free choice is `q0 = mean of min(r^4, tau^4)` over the calibration window. For normal returns with standard deviation `sigma` and `tau = a sigma` the closed form is `q0 = 3 kappa_4(a) sigma^4`, where `kappa_4(a)` is the share of the fourth moment that survives clipping (`docs/math.md`, section 3). Two wrong choices are worth naming. `q0 = 3 m0^2` (the Gaussian plug-in from the variance baseline `m0`) is 2.7% too high at `a = 3` and 21% too high with the raw MAD (`a = 2.02`), because clipping removes a larger share of the fourth moment than of the second. `q0 = m0^2` is 2.9 times too low. A `q0` that is too high makes the decrease detector anti-conservative and the increase detector blind; too low does the reverse.

**Bets and bet range.** Let `Y_t = Q_t / q0`, which lies in `[0, B_q]` with `B_q = tau^4 / q0`. For an increase the factor is `f_t = 1 + lambda_t (Y_t - 1)` with a predictable `lambda_t` in `[0, 1)`; it is positive because `Y_t >= 0`. For a decrease the factor is `f_t = 1 + lambda_t (1 - Y_t)`, positive only if `lambda_t < 1 / (B_q - 1)`. For normal returns and `a = 3`, `B_q = 81 / 2.89 = 28.0`, so the largest admissible decrease bet is `0.037`, against `B_var = 9 / 0.995 = 9.05` and `0.124` for the variance stream. The quarticity stream therefore has a 3.4 times smaller bet range on the decrease side. This is a structural cost of fourth powers: the bound `tau^4` is far above the typical value.

**Guarantee.** Under either null, `L_t^(j) = prod_{s=j}^t f_s` is a nonnegative supermartingale with initial value 1 for every start `j`, the Shiryaev-Roberts statistic `M_t = sum_j L_t^(j) = f_t (M_{t-1} + 1)` satisfies `E[M_t - t | past] <= M_{t-1} - (t - 1)`, and so (Theorems A and B in `docs/math.md`, section 5.4) the alarm time `inf{t : M_t >= A}` has `E[run length] >= A` and `P(alarm within m steps) <= m / A`. Nothing in the proof uses the form of the statistic beyond "nonnegative and bounded", so the quarticity stream inherits the guarantee unchanged. Two-sided monitoring runs each side at `2A`.

**What moves each stream.** For an alternative density `p` of the standardised return (pre-change standard normal, clipped at `a` standard deviations), the mean of the normalised stream is `E[Y] = E_p[min(|Z|^k, a^k)] / E_norm[min(|Z|^k, a^k)]` with `k = 2` or `4`, and its per-step signal-to-noise ratio is `snr = (E[Y] - 1) / sd(Y)`. The expected log growth of the best constant bet is about `snr^2 / 2` per step, so `snr` is the right single number to compare streams. `qstream.stream_shift` computes it by numerical integration for any density. Three facts that the experiments then confirm:

* A pure variance increase by a factor `v` moves `E[Y]` by about `v^2` on the variance stream and `v^4` on the quarticity stream (less, once clipping bites). The quarticity stream's noise is larger too, so the two snr values are close.
* A variance-matched Student t (fatter tails, same variance, `tau` fixed) *lowers* the truncated variance (more mass near zero, tail mass clipped) and raises the truncated quarticity only modestly, because clipping at `3 sigma` removes most of a fat tail. At `a = 3` the variance stream sees a larger shift (downwards) than the quarticity stream (upwards). Raising `a` to 4 or 5 changes the balance in favour of quarticity.
* A variance-matched scale mixture (two volatility states, equal chance, same mean variance) raises the truncated quarticity by 14% to 45% while the truncated variance moves by only 1% to 5%. This is where the quarticity stream is the right tool.

## 3. Why it matters

* It is the sequential counterpart of the quantity Bapat (2026) estimates. A detector of changes in truncated quarticity is a natural "part two" to a robust quarticity estimator, and sequential analysis is the author's core field.
* Changes in tail thickness and vol of vol are risk events in their own right. A desk that sizes positions by variance alone will be surprised by a regime where the variance is unchanged but large moves come in clusters. Options traders price vol of vol directly.
* It clarifies the role of truncation. Truncation makes decrease detection possible at all (`docs/math.md`, 5.3), but the bound it supplies grows with the power: `tau^2 / m0 = 9` against `tau^4 / q0 = 28`. The price of a higher power is paid on the decrease side.
* It shows what a variance monitor *misreads*: a fattening of the tails at fixed variance registers on the variance stream as a volatility *decrease*. Running both streams side by side separates "level" from "shape".

## 4. Comparison with existing approaches

| Approach | What it detects | Guarantee | Relation |
|---|---|---|---|
| Main package (variance stream) | Changes in truncated variance | Anytime-valid (Theorems A and B) | Same detector, second power. Faster for decreases. |
| This prototype (quarticity stream) | Changes in truncated quarticity: variance, tails, vol of vol | Same | Fourth power. Slightly faster for increases, slower for decreases, sees vol of vol. |
| Both streams together | Level and shape | Each side at `2A` per stream, or all four at `4A` | The combination `(M_var+ + M_var- + M_q+ + M_q-) / 4` is again an e-detector. |
| Two-sample betting test between windows (Shekhar and Ramdas 2024, A4) | Any change of distribution | Anytime-valid | More general, less interpretable, not direction-specific. |
| Kurtosis CUSUM or GARCH-type monitoring of tail parameters | Tail changes | Asymptotic or parametric | Needs a model for the pre-change law. The quarticity stream needs only a bounded moment. |
| Realized kurtosis tracking (RQ / RV^2 by day) | Shape changes | None (descriptive) | Offline and daily. |

Why not simply monitor the ratio `Q / X^2` (a truncated kurtosis)? A ratio of two streams is not bounded in a useful way and its conditional mean is not a clean function of the regime. Monitoring the two bounded moments separately keeps the guarantee and still lets the reader compare level with shape after the fact (`q_ratio / var_ratio^2` in the real-data output).

## 5. Implementation

`qstream.py` (imports `avmonitor`, nothing in `src` imports it):

* `variance_stream`, `quarticity_stream`: the two bounded statistics.
* `calibrate_streams`: scale, `tau`, `m0`, `q0`, the bounds `B_var`, `B_q` and the largest admissible decrease bets, from a calibration window, model free.
* `gaussian_null_moments`: the closed forms `kappa_2 sigma^2` and `3 kappa_4 sigma^4`, used to check the calibration.
* `run_stream_detector`: wraps `avmonitor.detectors.betting_edetector` with the right null mean, bound and grid for each stream. The increase grid `Q_UP_GRID = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5)` adds two small bets to the package default, because fourth powers are noisier; the experiments use this grid for both streams so that the comparison is fair.
* `stream_shift`, `student_pdf`, `mixture_pdf`, `scaled_normal_pdf`: the theory of section 2.
* `shape_change_returns`: simulator of unit-variance returns whose distribution switches at a known index (pure sd change, Student t, two-point scale mixture, slow stochastic volatility), with optional jumps in the whole path.
* `QuarticityMonitor`: streaming calibrate, monitor, alarm, recalibrate loop on `min(r^4, tau^4)`, same protocol as `avmonitor.monitor.VolatilityMonitor` (subclass), with `null_band` on `q0`.

`run_experiment.py` runs Q1 to Q4 and draws the two figures. `python run_experiment.py q3 q4` reruns parts.

## 6. Experimental methodology

Per-step units (`sigma0 = 1`), `tau = 3 x scaled MAD` from a 2000-step calibration window drawn from the same process as the pre-change data (jumps and fat tails are part of "normal", as in E4 of the main report), Shiryaev-Roberts statistics, fixed seeds, two worker processes.

* **Q1 Validity [SYNTHETIC].** No change, 20,000 monitored steps, `A = 7800` (so the bound is `P(false alarm within 390 steps) <= 5%` and `E[run length] >= 7800`). Gaussian, Gaussian with jumps (probability 0.005, size 8 sd), Student t5. 300 paths. Both streams, both sides. Output: `results/q1_validity.csv`.
* **Q2 Pure variance change [SYNTHETIC].** Change at step 500 of 3500, sd ratios 1.25, 1.5, 2 (increase side) and 0.8, 0.67, 0.5 (decrease side), Gaussian and jumps, 300 paths, `A = 7800`. Median delay, share detected, share of false alarms before the change. Output: `results/q2_variance_change.csv`.
* **Q3 Shape change with the variance fixed [SYNTHETIC].** Change at step 500, then 20,000 monitored steps, `A = 400,000` so that the null false alarm rate within the horizon is at most 5% (with `A = 7800` the null itself would alarm on about half the paths within 10,000 steps, which made a first version of this experiment unreadable; see `docs/iterations/ALT1.md`). Alternatives: Student t with 3, 4, 5, 8 degrees of freedom scaled to unit variance; two-point scale mixtures with `d = 0.5, 0.7, 0.9` (variance `1 - d` or `1 + d` with equal chance each step); slow log-volatility AR(1) with stationary sd 0.3 and 0.5 and half-life 60 steps, mean variance 1. 200 paths. Both streams, both sides. A theory table (`results/q3_theory.csv`) gives `E[Y]`, `sd(Y)` and `snr` for each stream at clipping levels `a = 3, 4, 5`.
* **Q4 Real data [REAL].** Binance spot BTCUSDT and DOGEUSDT one-minute bars, 1 July to 31 December 2025 (loader `avmonitor.data.load_binance`, 264,960 minutes each, no missing minutes), deseasonalised by the predictable 20-day intraday profile of the main package. Three monitors, each with `A = 28,800` per side (at most a 5% chance of any false alarm per 1440-minute day under the null) and a 1440-minute calibration: variance monitor with band 0.25, quarticity monitor with band 0.25, quarticity monitor with band 0.5625 (the same tolerance on the variance, since `1.25^2 = 1.5625`). For each alarm: the truncated variance and truncated quarticity over the next 1440 minutes relative to the calibration window (same `tau`), the distance in minutes to the nearest variance-monitor alarm, and whether a clipped return occurred in the hour before. Control: the same ratios at 2000 random times. The one event date named in `RESULTS.md` was checked online and is cited in `results/q4_real_event_check.csv`.

## 9. Folder structure

```
alternatives/quarticity_stream/
├── README.md             this file (items 1 to 6 and 9)
├── RESULTS.md            results with every number labelled [SYNTHETIC] or [REAL] (item 8)
├── qstream.py            prototype module (item 5, 7)
├── run_experiment.py     Q1 to Q4 and figures (items 6, 7)
├── tests/
│   └── test_qstream.py   7 pytest tests, about 1 second
└── results/
    ├── q1_validity.csv
    ├── q2_variance_change.csv
    ├── q3_shape_change.csv
    ├── q3_theory.csv
    ├── q4_real_summary.csv
    ├── q4_real_events.csv
    ├── q4_real_alarms_per_month.csv
    ├── q4_real_event_check.csv
    ├── fig1_synthetic_delays.png
    └── fig2_real_doge_alarms.png
```

Run from the repository root with the package installed (`pip install -e .`): `python alternatives/quarticity_stream/run_experiment.py` (about 5 minutes on two cores; the real-data part needs the Binance files under `data/raw`, see `data/README.md`) and `python -m pytest -q alternatives/quarticity_stream/tests`.

## References

* Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482.
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science*, 2(2), 229 to 260.
* Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *JRSS B*, 86(1), 1 to 27.
* Shekhar, S. and Ramdas, A. (2024). Nonparametric two-sample testing by betting. *IEEE Transactions on Information Theory*, 70.
