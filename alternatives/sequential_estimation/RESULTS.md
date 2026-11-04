# Results: sequential fixed-accuracy estimation of a volatility level

Every number below is copied from a CSV file in `results/`, produced by `run_experiment.py` with fixed seeds. [SYNTHETIC] marks simulated data; [REAL] marks Binance spot BTCUSDT one-minute bars (1 July to 31 December 2025). Formulas and definitions are in `README.md`. Iteration notes are in `docs/iterations/ALT6.md`.

## Summary

The question: keep reading one-minute returns until the volatility level is known to within plus or minus `d` (relative) with 95% (or 90%) confidence, then stop. The benchmark is the oracle size `n* = z² CV² / d²`, which needs the unknown shape of the returns. Main findings:

1. **Truncation makes the classical rules work on returns** [SYNTHETIC]. Chow and Robbins (1965) on `min(r², τ²)`, with `τ` three times the pilot's scaled MAD, covers at 0.936 to 0.964 (nominal 0.95) with 0.96 to 1.01 times the oracle sample, for normal, t5, t3 and jump returns, at every `d` from 0.4 to 0.05 (R3) and for any pilot from 10 to 100 returns (R6). On raw squares the same rule covers 0.72 to 0.78 for t3 at every `d`, 0.78 for jumps at `d = 0.2`, and its N grows like `d^{-3}` for t3 (R3, R4). No rule can guarantee coverage of the raw variance for every law (Proposition 1, README 2.2), so the truncated variance `θ_τ` is the right target; mapping it back to `σ²` costs about 1% under a normal shape but is badly biased under the wrong shape (2.5% coverage for t3, R9).
2. **The price of a guarantee that holds whenever one stops** [SYNTHETIC, closed form]. A betting confidence sequence, stopped when its width reaches the target, covers in every synthetic cell (at least 0.969), whatever the law. It needs `2 log(2/α) / z²` = 1.92 times the oracle as `d → 0` without the running intersection, and 1.548 with it (a Brownian first-passage limit, README 2.9), confirmed by simulation (1.545 at `d = 0.025`, R2); at practical accuracies 1.57 to 1.58 (`d = 0.05`) and 1.73 to 1.84 (`d = 0.2`) (R3). Time-uniform bets cost 5.6 to 14.6 times the oracle; bets matched to the target width are essential.
3. **Relative accuracy is the right way to state it** [SYNTHETIC]: in absolute accuracy the classical rules under-cover at `d = 0.2` even for normal returns (0.90, nearly all misses below), because of a skewness bias of the stopped mean that the relative rule removes (R5). Robbins' minimum risk rule on the truncated statistic is within 10% of the optimal risk at `n0 = 50` and within 3% from `n0 = 200`, with bounded regret (R8). Stein's two-stage rule needs a pilot of at least 30 (R6).
4. **Volatility clustering changes the question** [SYNTHETIC]. Under GARCH(1,1) every rule built for independent data, classical or anytime-valid, stops 7 to 47 times too early for the stationary level (coverage 0.28 to 0.82); batch means repairs it, but only at 2,000 to 19,000 returns (R7). What the classical rule does cover, at 0.948 to 0.966, is the average conditional truncated variance over the window it used (martingale central limit theorem). The running intersection of the confidence sequence is invalid for such a moving target (0.68 to 0.79); without it the sequences cover the window average at 0.986 to 1.000.
5. **BTCUSDT one-minute returns** [REAL, no coverage claim]. For about plus or minus 20% on volatility at 95%: about 50 minutes with truncated Chow and Robbins, 66 with raw squares, 140 with the betting sequence without intersection; for plus or minus 20% on the variance: 220, 411 and 505 minutes. The raw rule has not stopped by the end of its start hour in 58% of hours (42% to 69% by hour), against 27% for the truncated rule. The start hour changes the truncated rules' needs by less than 10% (R10).
6. **Bootstrap of quiet days** [REAL, resampled]: i.i.d. resampling of real returns reproduces the synthetic results (truncated Chow and Robbins 0.949, betting sequence 0.979, raw 0.919); keeping 30-minute blocks breaks every procedure (0.70 to 0.82), as clustering does in R7 (R11).

**Recommendation.** To quote a volatility level with a stated accuracy after an event, use truncated Chow and Robbins with relative accuracy and report the target as the window's average truncated variance; it is cheap (about the oracle sample) and well calibrated. When a guarantee at a data-dependent stopping time is required, use the betting confidence sequence without the running intersection, at 2.3 to 2.8 times the classical sample on real data. Never use raw squared returns, and never claim the stationary level from a few hours of clustered data.

## Names used in the tables

* **CR-raw, CR-trunc**: Chow and Robbins (1965) on `X = r²` or on `X = min(r², τ²)`, `τ = 3 ×` the pilot's scaled MAD (relative accuracy: Nadas' version on returns divided by the pilot scale). Interval `X̄_N ± d` or `X̄_N (1 ± d)`.
* **CR-BM-trunc**: the same with a batch-means long-run variance (GARCH part only).
* **Stein-raw, Stein-trunc**: Stein's (1945) two-stage rule with the pilot as first stage.
* **Oracle-raw, Oracle-trunc**: a fixed sample of `⌈n*⌉` with the true `ν` (a reference, not a procedure).
* **CS-wsr, CS-tuned, CS-opt**: the empirical Bernstein confidence sequence with the time-uniform bets of Waudby-Smith and Ramdas, with bets tuned to the target (`d_Z / σ̂_Z²`), and with constant bets that minimise the predicted stopping time; **CS-bet**: the betting confidence sequence. All four run on the truncated statistic after the pilot, use the running intersection unless marked "no intersection" (R2, R7, R10), and stop at the target width. Their `N` includes the pilot.
* **coverage**: share of paths with `|θ̂_N - θ| ≤ d θ` (relative) or `≤ d` (absolute), where θ is the procedure's own target (`E r²` for raw, `θ_τ` for truncated). **E[N]/n\***: mean over paths of `N / n*`, with `n*` the oracle size for the procedure's own target on that path. **regret**: mean of `N - n*`.
* Monte Carlo standard error of a coverage near 0.95 with 2000 paths: 0.005.

## R1. The price of anytime validity in closed form [SYNTHETIC, closed form]

`results/a_closed_form.csv`, `results/a_brownian_limit.csv`. Predicted `E[T]/n*` for the sequence alone (without the pilot), relative accuracy, truncation at `c = 3` times the true scaled MAD.

| law | quantity | d = 0.4 | d = 0.2 | d = 0.1 | d = 0.05 | d = 0.01 | limit, no intersection | limit, intersection (κ = 1) |
|---|---|---|---|---|---|---|---|---|
| normal, α = 0.05 | empirical Bernstein, optimised bets | 2.74 | 2.43 | 2.20 | 2.07 | 1.95 | 1.921 | 1.548 |
| normal, α = 0.05 | best constant-bet betting sequence for the known law | 1.80 | 1.89 | 1.91 | 1.92 | 1.92 | 1.921 | |
| t3, α = 0.05 | empirical Bernstein, optimised bets | 2.39 | 2.20 | 2.07 | 2.00 | 1.94 | 1.921 | 1.548 |
| t3, α = 0.05 | betting oracle | 1.68 | 1.86 | 1.91 | 1.92 | 1.92 | 1.921 | |
| normal, α = 0.10 | empirical Bernstein, optimised bets | 3.16 | 2.80 | 2.53 | 2.38 | 2.25 | 2.215 | 1.708 |
| normal, α = 0.10 | betting oracle | 2.07 | 2.18 | 2.21 | 2.21 | 2.21 | 2.215 | |

Brownian limit with the running intersection, by bet multiplier κ (40,000 paths, standard errors 0.0006 to 0.0018):

| α | κ = 0.8 | κ = 1 | κ = 1.2 | κ = 1.5 | coverage at the stop, κ = 1 |
|---|---|---|---|---|---|
| 0.05 | 1.721 | **1.548** | **1.475** | 1.512 | 0.976 |
| 0.10 | 1.919 | 1.708 | 1.603 | **1.601** | 0.952 |

**Reading.**
* Without intersection every sequence tends to `2 log(2/α) / z²`: 1.92 times the oracle at `α = 0.05`, 2.21 at `α = 0.10`. This factor is the ratio of the squared Chernoff quantile `2 log(2/α)` to the squared normal quantile `z²`, so it is the price of an exponential tail bound in place of the central limit theorem, not of time-uniformity (README 2.9).
* The running intersection lowers the limit to 1.548 and 1.708, and slightly larger bets (κ = 1.2 to 1.5) lower it further to about 1.48 and 1.60. The bets used by `CS-opt` and `CS-bet` are tuned for the rule without intersection (κ near 1), so about 5% is left on the table.
* At practical accuracies the empirical Bernstein bound pays an extra penalty (2.4 at `d = 0.2` for normal returns) because its bets hit the cap and `ψ` is far from quadratic; a betting sequence does not pay it (the oracle is even below the limit at `d = 0.4`, where it uses the skewness of Z). The penalty is smaller for t3 than for normal returns because truncation at the smaller `τ` of a t3 law gives Z a smaller range relative to its spread.
* The closed forms were checked by the tests (`2 log(2/α)/z² = 1.92056` and `2.21451`, coverage `0.99340`, the prediction tends to the limit, the Brownian limit without intersection equals its closed form) and by simulation in R2.

## R2. Simulated E[N]/n* against the closed forms [SYNTHETIC]

`results/j_verify_limits.csv`. Normal returns, relative accuracy, `α = 0.05`, pilot 100, 1000 paths. "Predicted" is the first-order prediction of README 2.9 with the exact `ψ` and caps, computed path by path from that path's `τ` and including the pilot. Cells: E[N]/n* (coverage at the stop).

| sequence | d = 0.1 | d = 0.05 | d = 0.025 | d → 0 |
|---|---|---|---|---|
| CS-tuned, no intersection: simulated | 2.740 (0.996) | 2.163 (0.995) | 2.018 (0.996) | |
| CS-tuned, no intersection: predicted | 2.626 | 2.150 | 2.015 | 1.921 (coverage 0.993) |
| CS-opt, no intersection: simulated | 2.362 (0.995) | 2.109 (0.994) | 2.007 (0.994) | |
| CS-opt, no intersection: predicted | 2.341 | 2.104 | 2.005 | 1.921 (coverage 0.993) |
| CS-bet, no intersection: simulated | 2.084 (0.993) | 1.963 (0.994) | 1.934 (0.992) | |
| CS-tuned, intersection: simulated | 2.091 (0.991) | 1.699 (0.983) | 1.597 (0.973) | 1.548 (coverage 0.976) |
| CS-opt, intersection: simulated | 2.023 (0.993) | 1.733 (0.988) | 1.629 (0.976) | 1.548 |
| CS-bet, intersection: simulated | 1.702 (0.978) | 1.576 (0.983) | **1.545** (0.968) | 1.548 |
| CR-trunc (for reference) | 0.998 (0.956) | 0.999 (0.953) | 1.001 (0.942) | 1 |

**Reading.** The closed forms are confirmed. Without intersection, the per-path prediction of the empirical Bernstein sequences is within 1% of the simulation for optimised bets at every `d` (4% for tuned bets at `d = 0.1`), and both tend to 1.92 with coverage at the stop of 0.99, as derived. With intersection the betting sequence reaches the Brownian limit already at `d = 0.025` (1.545 against 1.548) with coverage 0.968 (limit 0.976); the empirical Bernstein sequences approach it more slowly (1.60 and 1.63). So `E[N_cs]/n*` is, to a good approximation, **1.55 at α = 0.05 and 1.71 at α = 0.10 for a well-designed sequence at small d**, and the formulas of README 2.9 give it at practical `d`.

## R3. Relative accuracy, four laws [SYNTHETIC]

`results/b_main.csv`. 2000 paths per cell, pilot `m = 30`, `n_max = 600,000`. Each cell is "coverage / E[N]/n*" (for t3 raw rules `n*` is infinite, so the mean N is shown). Nominal coverage 0.95.

| law | procedure | d = 0.4 | d = 0.2 | d = 0.1 | d = 0.05 |
|---|---|---|---|---|---|
| normal | CR-raw | 0.929 / 0.93 | 0.948 / 0.97 | 0.951 / 0.99 | 0.942 / 1.00 |
| normal | CR-trunc | 0.937 / 0.98 | 0.952 / 0.99 | 0.953 / 1.00 | 0.942 / 1.00 |
| normal | Stein-raw | 0.934 / 1.04 | 0.937 / 1.03 | 0.927 / 1.01 | 0.931 / 1.02 |
| normal | Stein-trunc | 0.939 / 1.08 | 0.941 / 1.07 | 0.931 / 1.06 | 0.938 / 1.06 |
| normal | Oracle-trunc | 0.942 / 1.01 | 0.952 / 1.00 | 0.946 / 1.00 | 0.941 / 1.00 |
| normal | CS-wsr | 1.000 / 5.67 | 0.999 / 5.86 | 0.998 / 7.27 | 0.998 / 9.24 |
| normal | CS-tuned | 0.999 / 3.89 | 0.995 / 2.70 | 0.984 / 2.09 | 0.974 / 1.71 |
| normal | CS-opt | 1.000 / 3.56 | 0.997 / 2.45 | 0.987 / 1.97 | 0.976 / 1.73 |
| normal | **CS-bet** | 0.996 / 2.95 | 0.982 / 1.84 | 0.974 / 1.64 | 0.969 / **1.57** |
| t5 | CR-raw | 0.905 / 0.42 | 0.903 / 0.59 | 0.930 / 0.71 | 0.948 / 0.81 |
| t5 | CR-trunc | 0.936 / 0.97 | 0.949 / 1.00 | 0.956 / 1.00 | 0.944 / 1.00 |
| t5 | Stein-raw | 0.884 / 0.43 | 0.829 / 0.45 | 0.826 / 0.45 | 0.811 / 0.46 |
| t5 | Stein-trunc | 0.939 / 1.07 | 0.942 / 1.08 | 0.944 / 1.07 | 0.940 / 1.07 |
| t5 | CS-wsr | 0.998 / 6.23 | 0.998 / 7.47 | 0.994 / 9.72 | 0.997 / 12.59 |
| t5 | CS-opt | 0.998 / 3.21 | 0.992 / 2.23 | 0.984 / 1.83 | 0.975 / 1.67 |
| t5 | **CS-bet** | 0.992 / 2.34 | 0.980 / 1.75 | 0.976 / 1.63 | 0.970 / 1.57 |
| t3 | CR-raw | **0.782** / N 215 | **0.765** / N 2225 | **0.753** / N 18280 | **0.721** / N 126482 |
| t3 | CR-trunc | 0.938 / 1.01 | 0.944 / 1.01 | 0.954 / 1.00 | 0.949 / 1.00 |
| t3 | Stein-raw | 0.725 / N 132 | 0.589 / N 523 | 0.445 / N 2051 | **0.355** / N 8437 |
| t3 | Stein-trunc | 0.943 / 1.10 | 0.946 / 1.10 | 0.957 / 1.12 | 0.946 / 1.11 |
| t3 | CS-wsr | 0.998 / 6.80 | 0.996 / 8.67 | 0.994 / 11.29 | 0.999 / 14.60 |
| t3 | CS-opt | 0.995 / 3.04 | 0.990 / 2.15 | 0.986 / 1.79 | 0.978 / 1.66 |
| t3 | **CS-bet** | 0.990 / 2.23 | 0.973 / 1.73 | 0.976 / 1.62 | 0.973 / 1.57 |
| jumps | CR-raw | 0.918 / 0.25 | **0.783** / 0.47 | 0.832 / 0.81 | 0.929 / 0.96 |
| jumps | CR-trunc | 0.964 / 0.96 | 0.942 / 0.98 | 0.955 / 1.00 | 0.945 / 1.00 |
| jumps | Stein-raw | 0.877 / 0.23 | 0.656 / 0.24 | 0.590 / 0.24 | 0.595 / 0.24 |
| jumps | Stein-trunc | 0.968 / 1.07 | 0.929 / 1.05 | 0.932 / 1.05 | 0.933 / 1.05 |
| jumps | CS-wsr | 0.998 / 5.61 | 0.998 / 6.19 | 0.998 / 7.85 | 0.997 / 10.19 |
| jumps | CS-opt | 0.999 / 3.43 | 0.993 / 2.40 | 0.989 / 1.94 | 0.981 / 1.71 |
| jumps | **CS-bet** | 0.994 / 2.66 | 0.976 / 1.79 | 0.974 / 1.67 | 0.973 / 1.58 |

At `α = 0.10` (nominal 0.90) the picture is the same: CR-trunc 0.888 to 0.916 with ratios 0.96 to 1.07; CR-raw 0.655 to 0.703 for t3 and 0.719 to 0.861 for jumps; CS-bet coverage 0.939 to 0.995 with ratios 1.73 to 1.75 at `d = 0.05` (Brownian limit 1.708), 1.94 to 2.10 at `d = 0.2`; CS-tuned 1.78 to 1.89 at `d = 0.05`. No procedure reached `n_max` except CR-raw under t3 (2.9% of paths at `d = 0.05`, 0.2% at `d = 0.1`).

**Reading.**
* **Truncation is what makes the classical rules usable on returns.** CR-trunc has coverage between 0.936 and 0.964 in all 16 cells and `E[N]/n*` between 0.96 and 1.01: it behaves exactly as Chow and Robbins' asymptotics promise, already at plus or minus 40%. On raw squares the same rule fails for heavy tails and jumps, and the failure does not go away as `d` shrinks for t3 (0.78, 0.77, 0.75, 0.72: the self-normalised limit is not normal, README 2.3). For jumps it stops far too early (median N 342 against `n* = 1333` at `d = 0.2`) because a 1% jump rate is invisible in a few hundred returns. Stein's rule on raw squares is the worst procedure tested (0.36 for t3 at `d = 0.05`).
* **The confidence sequences hold their guarantee everywhere** (smallest coverage 0.969 at nominal 0.95 and 0.939 at nominal 0.90), as Theorem 1 says, for every law including t3.
* **Price.** CS-bet costs 1.57 to 1.58 times the oracle at `d = 0.05` for every law (Brownian limit 1.548), 1.73 to 1.84 at `d = 0.2`, and 2.2 to 3.0 at `d = 0.4`, where the 30-return pilot alone is 0.5 to 0.7 of `n*`. Empirical Bernstein bounds cost 5% to 10% more at `d = 0.05` and 20% to 40% more at `d = 0.2`. The time-uniform bets of Waudby-Smith and Ramdas are a poor choice for a known target: 5.6 to 14.6 times the oracle, and worse as `d` shrinks, because a bet schedule that spreads its power over all times wastes it when the stopping width is known in advance.
* The confidence sequences over-cover at their stopping time (0.97 to 1.00 at nominal 0.95), as predicted in README 2.9 (0.976 in the intersected limit).

## R4. Distribution of N and regret [SYNTHETIC]

`results/b_main.csv`, relative accuracy, `α = 0.05`.

| law | d | procedure | n* | E[N] | sd N | N 5% | N 50% | N 95% | regret |
|---|---|---|---|---|---|---|---|---|---|
| normal | 0.2 | CR-raw | 192 | 186 | 33 | 137 | 184 | 242 | -6 |
| normal | 0.2 | CR-trunc | 178 | 174 | 28 | 132 | 173 | 223 | -3 |
| normal | 0.2 | Stein-trunc | 178 | 186 | 55 | 106 | 180 | 288 | 8 |
| normal | 0.2 | CS-opt | 178 | 438 | 88 | 303 | 436 | 589 | 260 |
| normal | 0.2 | CS-bet | 178 | 328 | 62 | 233 | 327 | 434 | 150 |
| normal | 0.05 | CR-raw | 3073 | 3063 | 139 | 2843 | 3061 | 3297 | -10 |
| normal | 0.05 | CR-trunc | 2859 | 2853 | 271 | 2307 | 2908 | 3183 | -6 |
| normal | 0.05 | Stein-trunc | 2859 | 2978 | 900 | 1702 | 2870 | 4605 | 119 |
| normal | 0.05 | CS-bet | 2859 | 4482 | 659 | 3270 | 4558 | 5430 | 1623 |
| t3 | 0.1 | CR-raw | infinite | 18280 | 33311 | 3680 | 10002 | 53570 | |
| t3 | 0.1 | CR-trunc | 981 | 980 | 268 | 587 | 956 | 1439 | -1 |
| t3 | 0.1 | CS-bet | 981 | 1594 | 510 | 891 | 1526 | 2527 | 613 |
| jumps | 0.2 | CR-raw | 1333 | 624 | 595 | 150 | 342 | 1869 | -709 |
| jumps | 0.2 | CR-trunc | 203 | 198 | 42 | 141 | 191 | 277 | -5 |
| jumps | 0.2 | CS-bet | 203 | 363 | 90 | 239 | 352 | 530 | 160 |

(`n*` is averaged over paths; for truncated rules it varies from path to path with the pilot's `τ`, which is most of the spread of N for CR-trunc and the oracle.)

**Reading.** The classical rules have a small negative regret (they stop a few returns early on average, the usual behaviour of Chow and Robbins with the `1/n` term). Stein's rule has a spread of N three times that of Chow and Robbins, because its plan rests on 30 squared returns. Under t3 the raw rule's N has a standard deviation larger than its mean and grows like `d^{-3}` (215, 2225, 18280, 126482 for `d` = 0.4 to 0.05: factors of 10.3, 8.2 and 6.9 per halving of `d`, against the predicted 8), so in practice it has no planning value. The confidence sequence's regret grows like `n*` (a constant factor), the classical rules' regret stays bounded.

## R5. Absolute accuracy and the skewness bias [SYNTHETIC]

`results/c_absolute.csv`, `σ² = 1`, `α = 0.05`, 2000 paths. Each cell is "coverage (miss below, miss above) / E[N]/n*".

| law | procedure | d = 0.2 | d = 0.1 | d = 0.05 |
|---|---|---|---|---|
| normal | CR-raw | **0.899 (0.088, 0.012)** / 0.92 | 0.945 (0.037, 0.018) / 0.98 | 0.947 (0.031, 0.022) / 1.00 |
| normal | CR-trunc | 0.900 (0.086, 0.013) / 0.94 | 0.945 (0.035, 0.021) / 0.99 | 0.940 (0.033, 0.028) / 1.00 |
| normal | Stein-trunc | 0.905 (0.082, 0.013) / 1.11 | 0.917 (0.053, 0.029) / 1.10 | 0.906 (0.054, 0.040) / 1.10 |
| normal | CS-opt | 0.997 (0.001, 0.003) / 2.70 | 0.989 (0.003, 0.009) / 2.06 | 0.976 (0.011, 0.013) / 1.76 |
| normal | CS-bet | 0.983 (0.001, 0.016) / 1.96 | 0.974 (0.006, 0.019) / 1.73 | 0.967 (0.013, 0.020) / 1.61 |
| t5 | CR-raw | 0.777 (0.219, 0.004) / 0.50 | 0.901 (0.091, 0.007) / 0.69 | 0.920 (0.064, 0.017) / 0.80 |
| t5 | CR-trunc | 0.875 (0.117, 0.008) / 0.91 | 0.940 (0.044, 0.016) / 0.99 | 0.941 (0.040, 0.020) / 1.00 |
| t5 | CS-bet | 0.984 (0.006, 0.011) / 1.92 | 0.972 (0.011, 0.018) / 1.76 | 0.971 (0.013, 0.015) / 1.62 |
| t3 | CR-raw | **0.473 (0.527, 0.000)** | 0.641 (0.357, 0.001) | 0.692 (0.294, 0.013) |
| t3 | CR-trunc | 0.903 (0.087, 0.009) / 0.96 | 0.920 (0.070, 0.010) / 0.98 | 0.945 (0.034, 0.021) / 1.00 |
| t3 | CS-bet | 0.981 (0.004, 0.015) / 2.25 | 0.974 (0.008, 0.018) / 1.80 | 0.974 (0.012, 0.015) / 1.67 |
| jumps | CR-raw | 0.622 (0.378, 0.001) / 0.41 | 0.830 (0.161, 0.009) / 0.82 | 0.927 (0.054, 0.019) / 0.97 |
| jumps | CR-trunc | 0.901 (0.089, 0.009) / 0.92 | 0.944 (0.038, 0.018) / 0.99 | 0.954 (0.029, 0.017) / 0.99 |
| jumps | CS-bet | 0.979 (0.005, 0.015) / 1.88 | 0.971 (0.013, 0.015) / 1.73 | 0.972 (0.017, 0.011) / 1.61 |

**The skewness bias, measured [SYNTHETIC].** For normal returns, the first-order formula of README 2.3 predicts a relative bias of `-μ₃/(ν n*) = -4/n*` for the absolute rule and 0 for the relative rule. Measured mean of `X̄_N / θ - 1` for CR-raw (standard error about `d / (1.96 √2000)`):

| rule | d | n* | predicted | measured | standard error |
|---|---|---|---|---|---|
| absolute | 0.2 | 192 | -0.021 | -0.036 | 0.002 |
| absolute | 0.1 | 768 | -0.005 | -0.004 | 0.001 |
| absolute | 0.05 | 3073 | -0.001 | -0.001 | 0.001 |
| relative | 0.2 | 192 | 0 | -0.002 | 0.002 |
| relative | 0.1 | 768 | 0 | -0.001 | 0.001 |

**Reading.** In absolute accuracy the classical rules under-cover at `d = 0.2` even for normal returns (0.90), and almost every miss is an underestimate (0.088 below against 0.012 above): paths with a low running quarticity stop early and also have a low running variance. The first-order formula has the right sign and is accurate at `d = 0.1`; at `d = 0.2` (`n*` = 192) higher-order terms add about 70%. The relative rule removes the effect for normal returns, as predicted (misses 0.022 below and 0.030 above at `d = 0.2`, R3), which is one more reason to state volatility accuracy in relative terms. The confidence sequences do not have this problem: their misses are balanced and rare.

## R6. Pilot size [SYNTHETIC]

`results/d_pilot.csv`. Relative accuracy, `α = 0.05`, 2000 paths per cell, pilot `m` of 10, 30 or 100 returns. The pilot sets the truncation level `τ` and the scale of the confidence sequences, and is counted in `N`. Each cell is "coverage / E[N]/n*" (mean N for t3 raw, where `n*` is infinite).

| law | procedure | d = 0.2, m = 10 | m = 30 | m = 100 | d = 0.1, m = 10 | m = 30 | m = 100 |
|---|---|---|---|---|---|---|---|
| normal | CR-raw | 0.948 / 0.98 | 0.951 / 0.97 | 0.945 / 0.97 | 0.953 / 1.00 | 0.950 / 0.99 | 0.949 / 0.99 |
| normal | CR-trunc | 0.952 / 1.00 | 0.951 / 0.99 | 0.946 / 0.99 | 0.957 / 1.00 | 0.949 / 1.00 | 0.952 / 1.00 |
| normal | Stein-trunc | 0.939 / **1.33** | 0.944 / 1.05 | 0.948 / 1.02 | 0.940 / **1.32** | 0.934 / 1.08 | 0.950 / 1.01 |
| normal | CS-opt | 0.996 / 2.81 | 0.997 / 2.47 | 0.995 / 2.74 | 0.989 / 2.15 | 0.986 / 1.97 | 0.989 / 2.03 |
| normal | CS-bet | 0.979 / 1.88 | 0.981 / **1.84** | 0.983 / 2.16 | 0.976 / 1.75 | 0.976 / **1.65** | 0.976 / 1.71 |
| t3 | CR-raw | 0.753 / N 2325 | 0.746 / N 2059 | 0.749 / N 2187 | 0.735 / N 17840 | 0.746 / N 19432 | 0.750 / N 18673 |
| t3 | CR-trunc | 0.953 / 1.01 | 0.947 / 1.00 | 0.948 / 1.01 | 0.956 / 1.00 | 0.953 / 1.00 | 0.950 / 1.00 |
| t3 | Stein-trunc | 0.932 / 1.38 | 0.947 / 1.09 | 0.948 / 1.03 | 0.934 / 1.34 | 0.947 / 1.08 | 0.941 / 1.03 |
| t3 | CS-opt | 0.990 / 2.86 | 0.988 / 2.14 | 0.988 / 2.33 | 0.985 / 2.06 | 0.984 / 1.79 | 0.982 / 1.83 |
| t3 | CS-bet | 0.978 / 1.80 | 0.970 / 1.73 | 0.977 / 1.96 | 0.976 / 1.87 | 0.977 / 1.63 | 0.974 / 1.65 |
| jumps | CR-raw | 0.793 / 0.46 | 0.774 / 0.44 | 0.794 / 0.48 | 0.830 / 0.79 | 0.817 / 0.78 | 0.833 / 0.80 |
| jumps | CR-trunc | 0.952 / 0.98 | 0.940 / 0.98 | 0.952 / 0.99 | 0.949 / 0.99 | 0.955 / 0.99 | 0.945 / 1.00 |
| jumps | Stein-trunc | 0.919 / 1.23 | 0.929 / 1.05 | 0.945 / 1.02 | 0.917 / 1.24 | 0.928 / 1.04 | 0.942 / 1.01 |
| jumps | CS-opt | 0.997 / 2.80 | 0.994 / 2.39 | 0.993 / 2.63 | 0.992 / 2.15 | 0.985 / 1.93 | 0.985 / 1.98 |
| jumps | CS-bet | 0.985 / 1.77 | 0.972 / 1.78 | 0.978 / 2.09 | 0.974 / 1.80 | 0.969 / 1.66 | 0.976 / 1.70 |

**Reading.**
* **Chow and Robbins does not care about the pilot.** CR-trunc keeps coverage 0.940 to 0.957 and a ratio of 0.98 to 1.01 for every pilot from 10 to 100 returns: its stopping rule re-estimates the variance at every step, so the pilot only fixes `τ`, and a noisy `τ` changes the target (`θ_τ`) but not the quality of the rule for that target. CR-raw for jumps stays at 0.77 to 0.83 whatever the pilot: a bigger pilot does not reveal a 1% jump rate either.
* **Stein's rule pays for a small pilot**, as its theory says (its whole plan rests on the pilot variance): 23% to 38% more returns than the oracle with `m = 10`, 4% to 9% with 30, 1% to 3% with 100, with coverage 0.917 to 0.950.
* **The confidence sequences have an interior optimum near `m = 30`.** With `m = 10` the scale used to place the bets is rough (ratio 1.77 to 1.88 at `d = 0.2`); with `m = 100` the pilot is itself 0.4 to 0.6 of `n*` at `d = 0.2` and is spent without contributing to the sequence (sample splitting, README 2.8), so the ratio rises to 1.96 to 2.16. At `d = 0.1` the pilot is a smaller share and `m = 30` and `m = 100` are within 4%.
* Before the grid extension of the betting sequence (iteration ALT6-5), the run with `m = 10` had 1 to 6 paths in 2000 that never stopped and mean ratios of 8.4 (normal) and 15.9 (t3) at `d = 0.2` (`results/d_pilot_core_grid_only.csv`). The table above is after the fix: no path reached the cap, and the largest single-path ratio for normal returns at `m = 10`, `d = 0.2` was 4.0.

## R7. Volatility clustering [SYNTHETIC]

`results/e_garch.csv`. GARCH(1, 1) returns with normal innovations and unconditional variance 1 (Bollerslev 1986): `(a, b) = (0.10, 0.85)` (persistence 0.95) and `(0.05, 0.94)` (persistence 0.99), relative accuracy, `α = 0.05`, pilot 30, 1000 paths per cell. Two targets are scored:

* the **stationary level** `θ_τ`, the long-run truncated variance (what one wants to know when asking "what is the volatility of this asset"), whose oracle size uses the long-run variance of the truncated squares, `n*` = 1346 and 5370 (persistence 0.95, `d` = 0.2 and 0.1) and 5840 and 22979 (persistence 0.99), between 8 and 32 times the i.i.d. value;
* the **window average**: the mean of the conditional truncated variances `E[min(r_t², τ²) | past]` over the returns the procedure actually used (an analogue of integrated variance, Barndorff-Nielsen and Shephard 2002).

Each cell is "coverage of the stationary level / coverage of the window average / E[N]".

| procedure | persistence 0.95, d = 0.2 | 0.95, d = 0.1 | 0.99, d = 0.2 | 0.99, d = 0.1 |
|---|---|---|---|---|
| CR-raw | 0.521 / 0.966 / 212 | 0.491 / 0.958 / 952 | 0.382 / 0.955 / 197 | 0.293 / 0.958 / 951 |
| CR-trunc | 0.561 / 0.961 / 186 | 0.514 / 0.956 / 783 | 0.409 / 0.951 / 180 | 0.317 / 0.948 / 769 |
| CR-BM-trunc (batch means) | 0.971 / 1.000 / 2109 | 0.911 / 1.000 / 4992 | 0.993 / 1.000 / 10255 | 0.917 / 1.000 / 18916 |
| Oracle-trunc (long-run `n*`) | 0.946 / 1.000 / 1347 | 0.932 / 1.000 / 5371 | 0.951 / 1.000 / 5841 | 0.960 / 1.000 / 22980 |
| CS-opt, intersection | 0.626 / 0.949 / 383 | 0.510 / **0.793** / 996 | 0.440 / 0.967 / 368 | 0.295 / **0.696** / 824 |
| CS-bet, intersection | 0.577 / 0.924 / 277 | 0.455 / **0.713** / 774 | 0.421 / 0.969 / 281 | 0.279 / **0.682** / 666 |
| CS-opt, no intersection | 0.815 / 1.000 / 584 | 0.734 / 0.993 / 1900 | 0.549 / 0.999 / 580 | 0.449 / 0.996 / 2008 |
| CS-bet, no intersection | 0.756 / 0.997 / 448 | 0.695 / 0.992 / 1726 | 0.516 / 0.999 / 442 | 0.427 / 0.986 / 1815 |

The batch-means rule uses batches of 200 (persistence 0.95) and 1000 returns (0.99) and needs at least 10 batches, so it cannot stop before 2000 or 10000 returns; as a check, its estimate of the long-run variance of the raw squares is 20.5 against the closed form 22.7 (persistence 0.95) and 93.0 against 96.2 (0.99). No procedure reached the cap of 200,000.

**Reading.**
* **For the stationary level, every rule built for independent data stops far too early.** CR-trunc stops at 1/7 (persistence 0.95) to 1/32 (0.99) of the long-run `n*` and covers 0.32 to 0.56; CR-raw is worse (1/9 to 1/47, coverage 0.29 to 0.52). The confidence sequences fail in the same way (0.28 to 0.82): their guarantee is for a fixed target and i.i.d.-type martingale structure, and the stationary level is not the conditional mean of the next truncated square. Only a long-run variance (batch means, 0.91 to 0.99 at 0.82 to 1.76 times `n*`, partly because of its 10-batch minimum) or a known one (the oracle, 0.93 to 0.96) gets there, at 1,300 to 23,000 one-minute returns: about one to sixteen days for plus or minus 10% to 20% on a level that itself moves within days. Fixed-accuracy estimation of the *stationary* level is not a realistic intraday goal under clustering.
* **For the window average, the classical rule is right.** CR-trunc covers the average conditional truncated variance over the returns it used at 0.948 to 0.961, CR-raw at 0.955 to 0.966. This is the martingale central limit theorem: `Σ (X_t − E[X_t | past])` is a martingale, and the i.i.d. variance estimate `ν̂` is at least the average conditional variance, so the usual normal interval is valid, slightly conservative, for this moving target (README 2.10).
* **The running intersection is not valid for a moving target.** The intersected sequences cover the window average at only 0.68 to 0.79 at `d = 0.1`: an old bound, valid for the average over the returns seen up to that time, is kept although the target has moved since. Without intersection, Theorem 1 applies to the bet-weighted average of the conditional means (coverage of that weighted average 0.997 to 1.000 in `cs_theorem_cover`), and the sequences cover the plain window average at 0.986 to 1.000, at 1.5 to 2.7 times the stopping time of the intersected version. **Under volatility clustering, use the sequences without intersection, and state the target as the window average.**

## R8. Robbins' minimum risk rule [SYNTHETIC]

`results/f_minrisk.csv`. Loss `A (X̄_N - θ)² + N`, `A` chosen so that the optimal fixed sample is `n0` (for the truncated statistic `n0` uses the population `τ`; the benchmark `2 n0` is recomputed per path from that path's `τ`), pilot 10, `γ = 1`, 4000 paths. Each cell is "E[N]/n0 / risk efficiency `E[loss]/(2 n0)` / regret `E[loss] - 2 n0` (standard error)".

| law | statistic | n0 = 50 | n0 = 200 | n0 = 800 |
|---|---|---|---|---|
| normal | raw | 0.92 / 1.151 / 15.1 (1.7) | 0.98 / 1.036 / 14.6 (4.9) | 1.00 / 1.008 / 12.3 (18.6) |
| normal | trunc | 0.95 / 1.099 / 9.3 (1.4) | 0.99 / 1.006 / 2.1 (4.5) | 1.00 / 0.996 / -6.0 (16.8) |
| t5 | raw | 0.52 / 0.752 / -24.8 (0.8) | 0.66 / 0.895 / -42.0 (5.4) | 0.78 / 0.878 / -194.9 (17.6) |
| t5 | trunc | 0.93 / 1.101 / 9.9 (1.4) | 0.99 / 1.030 / 11.8 (5.2) | 1.00 / 1.005 / 7.1 (18.4) |
| jumps | raw | 0.36 / 0.454 / -54.6 (0.6) | 0.45 / 0.751 / -99.6 (4.3) | 0.65 / 1.015 / 23.5 (19.2) |
| jumps | trunc | 0.91 / 1.082 / 8.1 (1.4) | 0.98 / 1.018 / 7.2 (4.9) | 0.99 / 0.990 / -16.3 (18.9) |

**Reading.** On the truncated statistic Robbins' rule does what the theory says for every law: risk efficiency 1.08 to 1.10 at `n0 = 50`, within 3% of 1 at `n0 = 200` and within 1% at 800, with a regret of about 2 to 12 observations that does not grow with `n0` (bounded regret). On raw squares it does so only for normal returns. For t5 and jump returns the rule stops at 0.36 to 0.78 of `n0` and its realised risk is *below* the fixed-sample optimum `2 n0` (0.45 to 0.90). This is not a free lunch: the rule stops early on paths that have not yet seen a jump or a large t5 move, where it estimates the diffusive part and pays a squared bias that is smaller than the variance `ν/n` a fixed sample pays on average; `2 n0` is the best *fixed* sample, not a bound for sequential rules. For t5 the loss itself has an infinite variance (it involves `E r⁸`), so the Monte Carlo average and its standard error are not reliable. Truncation turns a rule whose behaviour depends on rare events into one that is close to its asymptotics at `n0 = 200`.

## R9. Back to σ² [SYNTHETIC]

`results/g_mapping.csv`. Inversion of `θ_τ = σ² g(τ/σ)` under a normal shape (README 2.4).

Closed form at the population `τ` (`c` times the true scaled MAD): the mapped value divided by the true variance (for jumps, by the diffusive variance and by the total variance).

| law | c = 2 | c = 3 | c = 4 |
|---|---|---|---|
| normal | 1.000 | 1.000 | 1.000 |
| t5 | 0.793 | 0.877 | 0.931 |
| t3 | 0.528 | 0.621 | 0.695 |
| jumps, diffusive variance | 1.031 | 1.055 | 1.088 |
| jumps, total variance | 0.825 | 0.844 | 0.871 |

Simulation, relative accuracy `d = 0.2`, `α = 0.05`, 2000 paths: the reported interval for `θ_τ` and its image for `σ²` (diffusive variance for jumps).

| law | procedure | covers θ_τ | image covers σ² | mapped estimate within 20% of σ² | median mapped / σ² |
|---|---|---|---|---|---|
| normal | CS-bet | 0.980 | 0.980 | 0.976 | 0.991 |
| normal | CR-trunc | 0.938 | 0.938 | 0.928 | 0.998 |
| t5 | CS-bet | 0.976 | 0.808 | 0.748 | 0.859 |
| t5 | CR-trunc | 0.929 | 0.680 | 0.723 | 0.868 |
| t3 | CS-bet | 0.977 | **0.025** | 0.018 | 0.616 |
| t3 | CR-trunc | 0.939 | 0.032 | 0.047 | 0.624 |
| jumps | CS-bet | 0.974 | 0.931 | 0.919 | 1.049 |
| jumps | CR-trunc | 0.940 | 0.926 | 0.871 | 1.047 |

**Reading.** With the right shape the map costs nothing in coverage (normal: 0.980 for `θ_τ` and for `σ²`) and about 1% in sample size at `c = 3` (elasticity 0.976, README 2.4). With the wrong shape the image is biased by exactly the closed-form factor: an exactly valid statement about the truncated variance becomes a 2.5% coverage statement about the variance of t3 returns. For jumps the truncated variance is close to the diffusive variance (5.5% above it), which is often the quantity of interest. The practical conclusion: report `θ_τ` (a robust scale with a guarantee) or a volatility in the units of a stated shape, never "σ²" without one.

## R10. Real data: minutes needed by hour of day [REAL]

`results/h_real_summary.csv`, `results/h_real_hourly.csv`, `results/fig3_real_minutes_by_hour.png`. Binance spot BTCUSDT one-minute log returns, 1 July to 31 December 2025 (loaded and checked by `avmonitor.data`: no duplicate or missing minutes), divided by the package's intraday volatility profile estimated from the previous 20 days only (`intraday_scale(..., start="past")`), so that a window that crosses a busy hour does not mix two volatility levels. One window per start hour: 4392 windows, each starting at the hour's first minute. The first 30 returns are the pilot; the rules then run for at most 1440 more returns (24 hours). Every count of minutes below includes the 30-minute pilot, and because no minute is missing, minutes equal returns. 4.3% of one-minute returns are exactly zero (the price did not move); in 1 window in 4392 the pilot's median absolute return was zero and the scale fell back to the median of the nonzero absolute returns.

**No coverage is claimed on real data**: the true level is unknown and, as R7 and R11 show, clustering moves it. The "n* proxy" is a descriptive yardstick: `z² CV² / d²` with `CV²` the squared coefficient of variation of `min(r², τ²)` over the 24 hours after the pilot (median over windows 2.44, against 1.92 for normal returns: BTC returns have heavier tails than normal returns even after truncation, and their volatility moves within a day).

`α = 0.05`. "d = 0.4" is plus or minus 40% on the truncated variance, which is about plus or minus 20% on volatility (a ratio between 0.77 and 1.18); "d = 0.2" is plus or minus 20% on the variance (about 10% on volatility). Cells: median over windows, with the range of the 24 hourly medians in brackets.

| d | n* proxy | procedure | minutes until stop | N / n* proxy | share not stopped within 60 minutes | not stopped within 24 hours |
|---|---|---|---|---|---|---|
| 0.5 | 38 | CR-raw | 38 (30 to 47) | 1.09 | **0.213** (0.12 to 0.38) | 0.002 |
| 0.5 | 38 | CR-trunc | 32 (30 to 33) | 0.94 | 0.017 (0.00 to 0.04) | 0 |
| 0.5 | 38 | Stein-trunc | 36 (32 to 37) | 1.00 | 0.011 (0.00 to 0.04) | 0 |
| 0.5 | 38 | CS-bet | 99 (95 to 106) | 2.86 | 0.933 (0.86 to 0.98) | 0 |
| 0.5 | 38 | CS-bet, no intersection | 106 (102 to 111) | 3.02 | 0.955 (0.88 to 0.99) | 0 |
| 0.4 | 59 | CR-raw | 66 (52 to 82) | 1.16 | **0.577** (0.42 to 0.69) | 0.003 |
| 0.4 | 59 | CR-trunc | **50** (46 to 52) | 0.91 | 0.272 (0.17 to 0.34) | 0 |
| 0.4 | 59 | Stein-trunc | 55 (49 to 58) | 0.97 | 0.375 (0.26 to 0.44) | 0 |
| 0.4 | 59 | CS-opt | 160 (147 to 175) | 2.92 | 0.995 (0.98 to 1.00) | 0.001 |
| 0.4 | 59 | CS-bet | 124 (114 to 132) | 2.27 | 0.976 (0.95 to 0.99) | 0 |
| 0.4 | 59 | CS-bet, no intersection | **140** (131 to 149) | 2.52 | 0.985 (0.96 to 0.99) | 0 |
| 0.3 | 104 | CR-raw | 138 (113 to 176) | 1.36 | 0.976 (0.95 to 0.99) | 0.012 |
| 0.3 | 104 | CR-trunc | 92 (86 to 100) | 0.94 | 0.914 (0.86 to 0.96) | 0 |
| 0.3 | 104 | CS-bet | 168 (149 to 180) | 1.75 | 0.995 (0.98 to 1.00) | 0.001 |
| 0.3 | 104 | CS-bet, no intersection | 221 (195 to 231) | 2.23 | 0.996 (0.98 to 1.00) | 0.001 |
| 0.2 | 234 | CR-raw | 411 (350 to 485) | **1.80** | 1.000 | **0.060** |
| 0.2 | 234 | CR-trunc | **220** (200 to 230) | 0.98 | 0.997 (0.98 to 1.00) | 0.001 |
| 0.2 | 234 | Stein-trunc | 220 (196 to 229) | 0.96 | 1.000 | 0 |
| 0.2 | 234 | CS-opt | 368 (336 to 401) | 1.71 | 1.000 | 0.003 |
| 0.2 | 234 | CS-bet | 284 (254 to 303) | **1.32** | 1.000 (0.99 to 1.00) | 0.001 |
| 0.2 | 234 | CS-opt, no intersection | 604 (558 to 662) | 2.65 | 1.000 | 0.046 |
| 0.2 | 234 | CS-bet, no intersection | **505** (447 to 551) | 2.20 | 1.000 (0.99 to 1.00) | 0.019 |

At `α = 0.10` everything is shorter by the ratio of squared quantiles (about 0.7): for d = 0.4, CR-trunc 35 minutes (32 to 37 by hour), CR-raw 42 (33 to 55), CS-bet without intersection 119 (115 to 126); for d = 0.2, CR-trunc 151, CR-raw 256 with 3.2% not stopped within 24 hours, CS-bet without intersection 420. At `α = 0.10`, `d = 0.5` the n* proxy (26) is below the pilot, so the classical rules stop as soon as the pilot ends.

**Reading.**
* **Minutes needed for about plus or minus 20% on volatility** (d = 0.4, 95%): about 50 minutes with truncated Chow and Robbins (46 to 52 by start hour), 55 with Stein, 66 with Chow and Robbins on raw squares (52 to 82), and 140 with the anytime-valid betting sequence without intersection (131 to 149). For plus or minus 20% on the variance (d = 0.2): 220 minutes truncated, 411 raw, 505 anytime-valid. The guarantee at the data-dependent stop costs a factor of 2.3 (d = 0.2) to 2.8 (d = 0.4) in minutes over truncated Chow and Robbins, of the same size as the synthetic factors of R2 and R3.
* **How often the raw rule fails to stop within the hour**: at d = 0.4 the raw rule has not stopped by the end of the start hour in 58% of hours (42% to 69% depending on the start hour), against 27% for the truncated rule; at d = 0.5, 21% (12% to 38%) against 1.7%. At d = 0.3 and below no rule stops within an hour (the n* proxy is above 60). At d = 0.2 the raw rule does not stop within 24 hours in 6.0% of windows, the signature of a variance estimate dominated by a few large squares (R3, R4). The raw rule needs 1.2 to 1.9 times the returns of the truncated rule, and the ratio grows as `d` shrinks, as it does for heavy-tailed synthetic laws.
* **The hour of day matters little** once returns are deseasonalised and accuracy is relative: the hourly medians of the truncated rules stay within 10% of their overall median, those of the raw rule within about 25%, with the raw rule slowest for windows starting at 21 to 22 UTC (and at 18 to 20 UTC for d = 0.2) and fastest around 13 to 16 UTC, when trading is heaviest. A relative accuracy depends on the shape of the truncated squares, not on the volatility level, so the time-of-day profile of volatility does not translate into a time-of-day profile of sample size.
* **The intersected confidence sequences stop suspiciously early on real data.** For normal i.i.d. returns CS-bet needs 1.84 times `n*` at d = 0.2 (R3); on BTC it needs 1.32 times the proxy with the running intersection and 2.20 without. The same gap appears under GARCH (R7: 277 against 448 returns) and in the block bootstrap (R11: 1.47 against 1.74), where the intersection is known to be invalid. The conservative reading is that on real data the intersected sequence is partly fitting the drift of the level, and the figures without intersection (505 minutes for d = 0.2, 140 for d = 0.4) are the ones to quote.
* **Windows that never stop** (1.9% for CS-bet without intersection at d = 0.2) are mostly those starting on Friday afternoons (50 of 85) whose 24 hours cover a quiet weekend: after the pilot the truncated mean is 0.30 in units of the pilot's squared scale (median 1.34 over all windows), 7.5% of returns are zero (median 2.4% over all windows) and the in-window `CV²` is 7.0 (2.4 over all windows), so the needed sample size is about three times larger. This was checked with a separate diagnostic run, not saved in `results/`.

## R11. Bootstrap of quiet days [REAL, resampled]

`results/i_bootstrap_summary.csv`, `results/i_bootstrap_days.csv`. The 30 quiet days of the same deseasonalised BTCUSDT series (complete days with the smallest dispersion of their 24 hourly robust scales: standard deviation of log scale 0.19 to 0.25) are resampled with the stationary bootstrap of Politis and Romano (1994), mean block length 1 (i.i.d. resampling of the day's returns) or 30 minutes (keeps volatility clustering within half an hour), 200 resamples per day, `d = 0.2`, `α = 0.05`, pilot 30. The target is the day's own truncated mean of `r²` at each path's `τ`, which is known exactly for the resampled law, so coverage can be measured. This is a check of the procedures on the empirical shape of real returns, not a coverage claim about real data. Cells: mean over the 30 days (smallest day in brackets).

| procedure | block 1: coverage | miss below | E[N]/n* | block 30: coverage | miss below | E[N]/n* |
|---|---|---|---|---|---|---|
| CR-raw | 0.919 (0.825) | 0.061 | 0.82 | 0.697 (0.520) | 0.200 | 0.79 |
| CR-trunc | **0.949** (0.920) | 0.021 | 1.00 | 0.748 (0.595) | 0.129 | 1.00 |
| Stein-trunc | 0.943 (0.895) | 0.027 | 1.08 | 0.742 (0.600) | 0.126 | 1.07 |
| CS-opt | 0.990 (0.970) | 0.008 | 2.22 | 0.823 (0.645) | 0.108 | 1.98 |
| CS-bet | **0.979** (0.950) | 0.010 | 1.74 | 0.748 (0.555) | 0.122 | **1.47** |

(`n*` here is the i.i.d. oracle size of the resampled day; with blocks the long-run variance is larger, so these ratios understate what a valid rule would need.)

**Reading.** With i.i.d. resampling the real shape of BTC returns behaves like the synthetic laws: truncated Chow and Robbins covers at 0.949 at the oracle size, the betting sequence at 0.979 with 1.74 times the oracle (synthetic: 1.73 to 1.84), and the raw rule under-covers (0.919, worst day 0.825) because it stops early on resamples without large moves (0.82 of `n*`), with three quarters of its misses below the target. Keeping 30-minute blocks breaks every procedure, the confidence sequences included (0.70 to 0.82): the dependence inflates the long-run variance, exactly as under GARCH (R7), and even on the quietest days of the sample. The intersected betting sequence stops earlier with blocks than without (1.47 against 1.74), which is the moving-target effect of R7 and R10.

## Negative results and failed attempts

Negative results (things that do not work, measured):

* **Classical rules on raw squared returns** under heavy tails or jumps (R3, R5): coverage 0.72 to 0.78 for t3 at every `d` (relative accuracy) and down to 0.47 (absolute), 0.78 for jumps at `d = 0.2`; N grows like `d^{-3}` for t3. Stein's rule on raw squares: down to 0.36. Robbins' minimum risk rule on raw squares stops at 36% to 78% of the optimal size for t5 and jumps (R8).
* **Absolute accuracy with the classical rules** under-covers at `d = 0.2` even for normal returns (0.90, nine misses in ten below the target), and the first-order bias formula is 70% too small there (R5).
* **The time-uniform bets of Waudby-Smith and Ramdas** cost 5.6 to 14.6 times the oracle for a known target width, and more as `d` shrinks (R3).
* **Bets tuned to the target** (`λ = d_Z / σ̂_Z²`, optimal as `d → 0`) are poor at practical `d`, where they hit the cap: 3.89 and 2.70 times the oracle at `d` = 0.4 and 0.2 (R3), against 2.95 and 1.84 for the betting sequence.
* **Mapping the truncated variance to σ² with the wrong shape** gives an interval for σ² that covers on 2.5% of paths for t3 returns (R9).
* **The bets used here are tuned for the rule without intersection**; with intersection a bet 1.2 to 1.5 times larger would save about 5% in the limit (R1). Not implemented.
* **Every i.i.d. rule, classical or anytime-valid, for the stationary level under GARCH**: stops 7 to 47 times too early with coverage 0.28 to 0.82 (R7). Batch means repairs it only at 2,000 to 19,000 returns.
* **The running intersection under a moving target**: the intersected confidence sequences cover the window average at 0.68 to 0.79 at `d = 0.1` under GARCH (R7), so the intersection, harmless for i.i.d. data, must be switched off under clustering.

Failed attempts during development (details in `docs/iterations/ALT6.md`):

* An optimiser grid for the constant bets that did not contain small bets gave a limit of 2.22 instead of 1.92 (found by a unit test, fixed).
* A first betting sequence chose its bets from a histogram with fixed bins of width 0.016 on [0, 1]; on a path whose pilot overestimated the scale 1.84 times (mean of Z 0.033) the upper bet collapsed to its minimum and the rule ran to the cap of 600,000 returns (ratio 195 on that path). A finer grid of candidate means did not help (wrong diagnosis); bins relative to the pilot mean did.
* Linear (chord) interpolation of the log capital between grid points was too loose; adding a second-order bound made the bounds exact to 0.05% with a 48-point grid.
* A grid of candidate means limited to `[ẑ_0/5, 5 ẑ_0]` stalls when the pilot is very small: with a 10-return pilot, 1 to 6 of 2000 paths per cell never stopped (ratios 8.4 for normal and 15.9 for t3 returns at `d = 0.2`, `results/d_pilot_core_grid_only.csv`). Extending the grid to the whole of (0, 1) with a per-point cap on the lower bet fixed it (R6) without changing any path whose bounds stay in the core.

## Ranking

Scores from 1 (low) to 10 (high).

| Criterion | Score | Why |
|---|---|---|
| Novelty | 6 | Every building block is classical: Chow and Robbins (1965), Stein (1945), Robbins (1959), Nadas (1967), betting confidence sequences (Waudby-Smith and Ramdas 2024), and truncation from the repository's own monitor. Our web searches (not exhaustive, README 4) did not find: the closed-form price of stopping a confidence sequence at a target width, `2 log(2/α) / z²`, with its reading as the Chernoff against central-limit price rather than a price of time-uniformity; the Brownian limit 1.548 with the running intersection; the sequential Bahadur and Savage argument for volatility (Proposition 1); or a sequential fixed-accuracy study of truncated volatility under heavy tails, jumps and clustering. Each piece is small. |
| Validity | 8 | Theorem 1 is proved for any law and holds in every i.i.d. cell; the closed forms are confirmed to within 1% by simulation and by the tests; the classical rules match their asymptotics. Minus two points: under clustering only the window average is covered, and only without the intersection; the classical rules are asymptotic, with a measured under-coverage of up to 0.10 in absolute accuracy at `d = 0.2`; the mapping to `σ²` needs a shape assumption. |
| Practical relevance | 6 | Answers a concrete question on real data (about 50 minutes for plus or minus 20% on BTC volatility, 140 with a guarantee), and gives a clear warning about raw squares and about clustering. Minus points because the level moves within hours, so a fixed-accuracy estimate of it is a snapshot; for detecting that it moves, the main monitor is the right tool. |
| Feasibility | 9 | About 1300 lines of numpy and scipy, no new dependency, every rule streams in constant time per return (the betting sequence in time proportional to its grid of about 136 points), and the whole study runs in about 17 minutes on two worker processes. |
| Demonstration value | 8 | Figure 1 (simulated cost against the two closed-form limits), Figure 2 (truncated against raw coverage by law) and Figure 3 (minutes by start hour on BTC) carry the three messages; the GARCH table shows the moving-target issue in one glance. |
| Interest to a sequential-analysis and high-frequency econometrics researcher | 8 | For sequential analysis: an exact, interpretable price of anytime validity at a target width, a Brownian first-passage limit for the intersected sequence, and the observation that the running intersection breaks for a moving target. For high-frequency econometrics: the link between a sequential rule under clustering and an integrated-variance-type target, the `d^{-3}` growth of raw rules under infinite quarticity, and real-data sample sizes by hour. Open: bets tuned for the intersected rule (about 5% left), a long-run variance inside the sequence, and the pilot's sample-split cost. |

Overall: a solid, mostly negative-to-clarifying direction. It shows that the classical theory is fine once returns are truncated and accuracy is relative, quantifies exactly what the anytime-valid guarantee costs, and makes clear that the real limitation is not the stopping rule but the moving level.
