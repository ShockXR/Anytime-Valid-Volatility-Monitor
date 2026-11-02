# From the paper's equations to the code

This table lets a reader go from each equation in Bapat (2026) to the function that implements it, and to the test or experiment that checks it. Equation numbers are the paper's. "As printed" means the formula is coded literally, including the features questioned in `docs/estimator_issues.md`. "Normalised" means the version on the integrated quarticity scale used by the prototype.

## Part 1. The paper

| Eq. | What it is | Code | Checked by |
|---|---|---|---|
| (1), (8), (9) | Price model: diffusion plus jumps, returns as increments | `simulate.sv_jump_returns`, `simulate.change_returns` | `tests/test_robust.py`, E2, E4 |
| (2) | Integrated quarticity `IQ = ∫ σ_t⁴ dt` | target computed from the simulated σ path in `simulate.SVPaths` | E2c |
| (3) | Realized quarticity `RQ = (n/3) Σ r⁴` | `robust.realized_quarticity` (`T` sets the scale) | `test_classical_estimators_unbiased_constant_vol` |
| (4) | Bipower quarticity, as printed | `robust.bipower_quarticity` (also the corrected constant, see `docs/estimator_issues.md` 4d) | E5a, E6 |
| (5) | Tripower quarticity | `robust.tripower_quarticity` (p = 4/3; other p in E5a, E6) | unit tests, E5a |
| (Table 3) | MinRQ and MedRQ (Andersen, Dobrev and Schaumburg 2012) | `robust.min_rq`, `robust.med_rq` | constants guarded by unit tests (iteration 2 fixed MinRQ) |
| (Table 4) | Threshold RQ (Mancini type) | `robust.threshold_rq` | `test_threshold_rq_ignores_large_jump` |
| (6) | Local kernel quarticity, as printed (unnormalised weights, no `1/Δ²`) | `robust.local_quarticity_printed` | `test_printed_kernel_scale_is_weight_sum_over_n_squared` |
| (6), normalised | Spot quarticity `Σ w_k r⁴ / (3Δ²)` with `Σ w_k = 1` | `robust.local_quarticity` | `docs/math.md` 2, unit tests |
| (7) | Average of eq. (6) over the sample | `robust.kernel_quarticity_printed`; normalised: `robust.kernel_iq(c=None)` | E1, E5a, E6 |
| (10), (13) | Simulation design: `r_t = σ ε_t + J_t (+ α_t)` | `simulate.paper_returns` | E1 |
| (11) | Truncation `r̃ = min(max(r, -τ), τ)`, `τ = c · MAD(r)` | `robust.truncate`, `robust.mad`, `robust.mad_threshold` (raw or scaled MAD; optional floor; warning when MAD = 0) | E2a, `tests/test_edge_cases.py` |
| (12), (20) | Local robust kernel quarticity | `robust.local_quarticity_printed(truncate(r, τ))`; normalised: `robust.local_quarticity` on truncated returns | E1, E3 |
| (14) | Gaussian kernel weights `K(i/h)` | `robust.kernel_weights`, `robust.kernel_weight_sum`; all kernels of Table 2 in `robust._kernel_values` | `test_kernel_weight_sum` |
| (15) | Weighted sum with index `w_{s-t+h}` (reverses the weight order) | `robust.local_quarticity_printed(..., reverse_weights=True)` | `docs/estimator_issues.md` 4a |
| (16), (21) | Robust kernel estimator, as printed | `robust.robust_kernel_quarticity_printed` (raw or scaled MAD) | E1, E5a, E6 |
| (16), normalised | Robust kernel IQ, optionally divided by the surviving fourth-moment share κ₄ | `robust.kernel_iq(c=3, mad_scaled=True, bias_correct=True)` | E2b, E2c |
| (16) without the kernel | Truncated realized quarticity (shown equal to eq. 16 up to edge effects) | `robust.truncated_rq_mad` | E2d, `docs/math.md` 4 |
| (17), (18) | Cross-validated bandwidth | not implemented (the averaged estimator does not depend on h beyond edge effects, `docs/math.md` 4) | |
| (19), Thm. 1 | Consistency | not a code item. `docs/math.md` 3.1 shows that with `τ = c · MAD` the truncation bias does not vanish, so the limit is `κ₄ · IQ` | E2a |
| (22), Thm. 2 | Asymptotic normality | not checked (the variance AVAR is not given in the paper) | |

The share of the fourth moment that survives truncation, `κ₄(a)`, is `robust.truncated_moment_factor`, and `robust.effective_cutoff(c, scaled)` converts `c` into the cut-off in standard deviations (`docs/math.md` 3).

## Part 2. The monitor (our extension)

| Item in `docs/math.md` | What it is | Code | Checked by |
|---|---|---|---|
| 5.1 | Statistic `X_t = min(r_t², τ²)` and the null `E[X_t ∣ past] ≤ m₀` (or `≥`) | `detectors.calibrate`, `monitor.VolatilityMonitor._calibrate` | `tests/test_detectors.py` |
| 5.2 | Bets `1 + λ (X_t / m₀ - 1)` (increase) and `1 + λ (1 - X_t / m₀)` (decrease) | `detectors.betting_edetector`, `VolatilityMonitor.update` | `test_bets_have_mean_at_most_one` |
| 5.3 | No decrease detection without a bound on X | (a proof; the decrease bet needs `upper`) | `betting_edetector` raises without `upper` |
| 5.4, Theorems A and B | Shiryaev-Roberts recursion `M_t = f_t (M_{t-1} + 1)`, ARL at least A, `P(false alarm within m) ≤ m / A` | `detectors.betting_edetector(kind="SR")`, two sides at `2A` in `detectors.two_sided_alarm` and in the monitor | E4a, E9 (bound tightness at scale) |
| 5.5 | Bet sizing from local truncated quarticity (optional) | `betting_edetector(adaptive_h=...)` | E4, Q4c |
| 5.6 | Estimated baseline: margin and null band | `calibrate(margin_z=..., null_band=...)`, `VolatilityMonitor(null_band=...)` | E4c, E4d, E9c |
| Intraday profile | Predictable per-minute scale from past days | `data.intraday_scale` | E5, E6, E7 |
| Online MAD | Rolling and exponentially weighted scale | `robust.rolling_mad_zero`, `robust.ew_quantile_scale` | E3, `tests/test_edge_cases.py` |

## Part 3. Keeping the guarantee with an estimated baseline (`docs/math.md`, section 6)

| Item | What it is | Code | Checked by |
|---|---|---|---|
| Lemma 6.1 | Exponential inequality `exp(λ(ν - Z) - ψ_a(λ) e) ≤ 1` in expectation for a bounded Z | `confseq.psi_coef` | `tests/test_confseq.py` |
| Theorem 6.2, Corollary 6.3 | Time-uniform empirical-Bernstein confidence sequence for the mean of a bounded variable | `confseq.eb_confseq`, `confseq.StreamingEBConfSeq` | `tests/test_confseq.py` (coverage on i.i.d. and time-varying means) |
| Section 6.3, Theorem 6.4 | Detector baselines from the confidence sequence; `P(false alarm within m) ≤ δ + m / A` | `confseq.cs_baseline`, `VolatilityMonitor(baseline="cs")` | `tests/test_monitor_cs.py`, E8 |
| Corollary 6.6 | Recalibration after every alarm | `VolatilityMonitor` restart loop with `baseline="cs"` | E8 (alarm counts) |
| Section 6.4 | τ from data not used by the confidence sequence | `confseq.split_tau` | `tests/test_confseq.py` |

## Part 4. Forecast-relative monitoring (E10)

| Item | What it is | Code | Checked by |
|---|---|---|---|
| Standardised returns `r_t / sqrt(h_t)` with a predictable variance forecast | EWMA and GARCH(1,1) forecasts | `forecast.ewma_variance`, `forecast.garch11_filter`, `forecast.garch11_fit`, `forecast.standardise` | `tests/test_forecast.py` |
| Signal cap `ρ(k) = k(1 - β) / (1 - α - β + αk) ≤ (1 - β) / α` | Share of a k-fold variance break that a frozen GARCH forecast passes to the standardised returns | `experiments/e10_forecast_relative.rho` | E10 Q10d (simulated against closed form) |
| Information rate `I = E[½(q - 1 - ln q)]` and delay `ln A / I` | Lai (1998) first-order delay for the standardised stream | `experiments/e10_forecast_relative.q10d` | E10 Q10d, Q10e |

## Part 5. Detection delay and bet design (`docs/theory_delay.md`)

| Item | What it is | Code | Checked by |
|---|---|---|---|
| Section 2, Lemma 0 | Growth rate `g(λ) = E_Q log(1 + λU)` of one bet | `design.growth_rate`, `design.log_factor_moments` | `tests/test_delay_theory.py` |
| Lemma 2 | Best nonparametric growth `KL_inf` (oracle bet) | `design.best_bet`, `design.approx_best_bet` | `tests/test_delay_theory.py`, E11 |
| Corollary 3 | Information ladder: KL, KL of the truncated statistic, KL_inf, grid growth | `design.kl_gaussian_scale`, `design.kl_truncated`, `design.growth_rate` | E11b |
| Corollary 4 | Speed limit for decreases: at most `log(B / (B - 1))` nats per step for every valid bet; floor of a given grid, and floor `log(1 + A / B) / log(B / (B - 1))` for any valid bets | `design.lam_max`, `design.speed_limit` (a given grid), `design.decrease_speed_limit` and `design.sr_steps_constant_factor` (any valid bets) | `tests/test_delay_theory.py::test_speed_limit_is_attained_by_extreme_data`, `tests/test_decrease_speed_limit.py` |
| Theorem 1 | Non-asymptotic delay bound `min_k (log(A / w_k) + overshoot_k) / g_k` | `design.wald_delay_bound`, `design.overshoot_bound` | E11 (1.3 to 2.7 times the simulated delay) |
| Section 4.3 | Second-order delay `(log A + κ - ζ) / g` | `design.delay_second_order`, `design.renewal_constants`, `design.zeta_brownian` | E11 (within 3% in 20 cases) |
| Section 7 | Optimised bet grids | `design.optimise_grid`, `design.design_ratio` | E11c |
| Section 8 | Worst-case false alarm probability over all null laws, by dynamic programming | `design.worst_case_false_alarm` | E11d |
