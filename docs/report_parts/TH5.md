## TH5. Limits for every valid procedure, and the plug-in baseline when n, m and A grow together

Theory `docs/theory_limits.md` (statements, proofs and the full tables); code `src/avmonitor/limits.py`; script `experiments/th5_limits.py`; outputs `results/th5_*.csv` and `results/figures/th5_floors.png`, `results/figures/th5_plugin_regime.png`; tests `tests/test_limits.py`; iteration notes `docs/iterations/TH5.md`. The theorems are proved in the theory note; all numerical results are **[SYNTHETIC]**, with seeds fixed in the script.

**Question.** Sections 2, 3 and M3 of the referee report. (1) The decrease speed limit of `docs/theory_delay.md` (58 steps on a complete collapse for any Shiryaev-Roberts e-detector with linear bets at c = 3, A = 7800) is a statement about e-detectors. How fast can *any* procedure with the guarantee of Theorem B detect a change, and how far are the e-detectors from that? (2) When does the plug-in baseline keep the nominal false alarm level if the calibration length n, the horizon m and the alarm level A grow together, and what margin restores it? (3) Why is the growth-optimal decrease bet so close to its upper limit `λ_max = 1 / (B - 1)`?

**Set-up.** Monitored statistic `Y = min(r², τ²) / m₀ ∈ [0, B]`, nonparametric null classes (conditional mean of Y at least 1, at most 1, or equal to 1; laws need not be independent or identical), validity `P(T ≤ n) ≤ n / A` for every n (Theorem B). Default detector: six increase bets, five decrease bets, each side at 2A, A = 7800, c = 3. Post-change laws for delays: Gaussian, Student t3 and the jump law of E4, sd ratios 0.5, 0.67, 0.8, 1.25, 1.5 and 2, change at the first step, true baseline. Plug-in: Gaussian and Student t3 returns, n = 390 to 92,160, m = 98 to 7800.

**Summary.**
* **A speed limit for every valid procedure (Theorem 2, proved).** On a path y, a rule valid at levels `n / A` can stop at step n only if `W_n(y) = Π_{s ≤ n} 1 / π(y_s) ≥ A / n`, where `π(y)` is the largest probability a null law can put on the value y: `(B - 1) / (B - y)` for `y < 1` (decreases), `1 / y` for `y > 1` (increases). The proof places an independent two-point null law at each step, under which the path has probability exactly `1 / W_n(y)`. The floor is exact path by path, holds for rules that see the raw returns, and a randomised rule can beat it only with probability at most `(n / A) W_n(y)`.
* **Complete collapse: 45 steps for every deterministic procedure (Corollary 4, proved).** With all returns zero at c = 3 (B = 9.045), no deterministic valid procedure alarms before step 45 at A = 7800 (50 at 2A), against 57.7 for every Shiryaev-Roberts e-detector with linear bets and 73 for the default grid (79 at 2A). A randomised valid procedure can alarm by step n on this path with probability at most `min(1, n W_n / A)` (0.98 at step 44), and its mean delay there is at least 43.4 (48.5 at 2A; a linear programme). The 13 steps between 45 and 58 are the factor `n / B` in the level that one path allows and an e-detector, which must spread `E M_n ≤ n` over all paths, cannot use. At c = 2 the numbers are 23, 28.7 and 36. Without truncation no deterministic valid procedure alarms before step A on any path (a randomised one by step n < A with probability at most n / A), the sequential form of `docs/math.md` 5.3.
* **The e-detectors are within a factor of about 2 of the best possible delay (Theorem 5, proved; numbers Monte Carlo).** Combining the pathwise floor with a Neyman-Pearson bound at the least favourable null law gives a lower bound on the mean delay of every valid deterministic procedure. Over the 18 laws, one detector side at level A takes 1.69 to 2.16 times the bound (default grid) and 1.34 to 1.76 times (growth-optimal single bet); against the bound for randomised procedures, 1.81 to 2.24 and 1.53 to 1.80. These ratios bound the inefficiency from above, since no rule is known to attain the bound. As a consistency check of the code, no alarm of either, on 180,000 paths at two levels, came before the pathwise floor.
* **The plug-in keeps its level exactly when m is small compared with √n (Theorem 8, Corollary 10, proved).** The drift of the plug-in baseline is at most `B̂ D_n`, with `D_n` the Kolmogorov distance between the window and the return law. So, for every stream in the class `C(W)` of streams whose conditional truncated mean at τ̂ equals that of the return law, `P(T ≤ m) ≤ (m / A)(1 + C (m + 1) / √n (1 + o(1)))` when `m / √n → 0`, with `C = 0.037` for Gaussian and 0.047 for Student t3 returns (first order) and `C = 1.79` for every law. Conversely, a stream in `C(W)` that holds the drift and then bets boldly makes `P(T ≤ m) / (m / A) → ∞` when `m / √n → ∞`, `m / log A → ∞` and `m / A → 0`. The boundary is `m ≍ √n` for any A with `log A = o(m)`. At one day after a 390-step window of Gaussian returns, the largest false alarm probability over `C(W)` is at least about 33% and at most about 42% (nominal 5%): these are window averages of proved per-window bounds, estimated by Monte Carlo (32.3% and 41.6% over 20,000 windows in the CSV; 32.3% to 33.2% and 41.6% to 42.2% over six runs of 20,000 windows; standard errors 0.33 and 0.30 points). The hold-then-bold stream, simulated, gives 36.5%. For i.i.d. monitored returns the excess is much smaller and falls quickly with n (simulated, Q-TH5c).
* **The margin for the plug-in is the normal quantile (Proposition 11, proved).** Baselines `m̂ ± z se` keep `P(T ≤ m) ≤ m / A + 2 Φ̄(z) + o(1)` uniformly in m and A, and the worst stream attains `2 Φ̄(z)` when `m / (√n log A) → ∞`; so `z = Φ⁻¹(1 - δ / 2)`. The confidence-sequence baseline (Theorem D) guarantees the same without asymptotics at an effective z of 3.8 to 5.7.
* **The decrease optimum is at the boundary only for light tails (Proposition 12, proved).** The growth of the bet `λ_max (1 - ε)` is `G_c + q_B log ε` up to small terms, with `q_B` the probability that a return reaches the truncation level, and `ε* = q_B / h(ε*)`. For Gaussian halving `q_B = 2 × 10⁻⁹` and `ε* = 2.4 × 10⁻⁸`, the referee's value. For Student t3 `q_B` is 0.6% to 2.4% and `λ* / λ_max` is 0.93 to 0.55; a bet at `λ_max (1 - 10⁻⁶)` has negative growth there and is up to 53% slower.

### Q-TH5a. Floors on the most extreme paths

Gaussian truncation, `results/th5_collapse.csv`. At A = 7800 (at 2A in brackets).

| c | B | complete collapse: every deterministic valid procedure | randomised valid procedure, smallest mean | SR e-detector, linear bets | default grid | complete explosion: every deterministic valid procedure | randomised, smallest mean | SR e-detector | default grid |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 4.345 | 23 (25) | 22.1 (24.4) | 28.7 (31.3) | 36 (39) | 6 (6) | 5.0 (5.5) | 6 (7) | 11 (12) |
| 2.5 | 6.393 | 33 (36) | 31.8 (35.3) | 41.8 (45.9) | 52 (57) | 5 (5) | 4.1 (4.6) | 5 (6) | 8 (9) |
| 3 | 9.045 | 45 (50) | 43.4 (48.5) | 57.7 (63.6) | 73 (79) | 4 (4) | 3.7 (3.9) | 5 (5) | 7 (7) |
| 4 | 16.00 | 73 (82) | 70.9 (80.0) | 95.9 (106.7) | 122 (134) | 3 (4) | 3.0 (3.2) | 4 (4) | 5 (6) |

A rule that promises only `P(T ≤ 390) ≤ 0.05` (one day) could alarm on the collapse after 26 steps at c = 3: controlling every n at once costs 19 steps on this path.

### Q-TH5b. Delay floor for every valid procedure against the e-detectors

10,000 paths of 4000 steps per law, A = 7800, mean delays in steps (`results/th5_floors.csv`). "Bound" is the lower bound for deterministic valid rules (Theorem 5: the larger of the pathwise floor and the Neyman-Pearson bound restricted to the floor).

| Law | sd ratio | pathwise floor | bound | default grid | best single bet | grid / bound | best bet / bound |
|---|---|---|---|---|---|---|---|
| Gaussian | 0.5 | 55.2 | 55.2 | 93.2 | 74.2 | 1.69 | 1.34 |
| Gaussian | 0.8 | 70.1 | 90.4 | 194.9 | 156.5 | 2.16 | 1.73 |
| Gaussian | 1.25 | 16.8 | 65.4 | 120.7 | 108.3 | 1.85 | 1.66 |
| Gaussian | 2 | 9.1 | 12.2 | 23.0 | 18.4 | 1.88 | 1.51 |
| Student t3 | 0.5 | 39.7 | 44.9 | 81.5 | 69.3 | 1.82 | 1.55 |
| Student t3 | 0.8 | 47.8 | 147.0 | 281.3 | 259.3 | 1.91 | 1.76 |
| Student t3 | 1.25 | 18.7 | 117.3 | 218.4 | 196.5 | 1.86 | 1.67 |
| jumps | 0.5 | 54.1 | 56.9 | 99.0 | 82.6 | 1.74 | 1.45 |
| jumps | 2 | 9.2 | 12.5 | 23.5 | 18.9 | 1.89 | 1.51 |

All 18 laws are in the theory note (section 3.3). Over all of them, one side at A against the bound for rules valid at `n / A`: grid / bound 1.69 to 2.16, best bet / bound 1.34 to 1.76 (at 2A against `n / (2A)`: 1.65 to 2.09 and 1.33 to 1.72). The two-sided monitor runs each side at 2A and is valid at `n / A`; its sides take 1.85 to 2.40 (grid) and 1.49 to 2.03 (best bet) times the bound at A, and the bound for two-sided validity, which was not computed, can only be smaller. Randomised rules gain at most 15% on the bound (grid / randomised bound 1.81 to 2.24, best bet 1.53 to 1.80). Batch-means Monte Carlo standard errors (10 batches of 1000 paths, `se_*` columns) are at most 0.7% of the bound and at most 0.012 for a ratio. As A grows the gap closes slowly: the best bet over the bound is 1.38 to 1.93 at `A = 10³` and 1.22 to 1.40 at `A = 10⁸` (`results/th5_asymptotic.csv`, six laws). The pathwise floor alone grows like `log A / γ_Q` with `γ_Q ≥ KL_inf`, with equality only when the post-change law has no mass beyond 1 (Proposition 6), so for increases and fat tails the Neyman-Pearson part of the bound does the work.

### Q-TH5c. The plug-in baseline: bounds and simulation

A = 7800, c = 3, each side at 2A. Bounds from calibration windows only (20,000 windows per cell for n ≤ 23,040, 5000 above; `results/th5_plugin_bound.csv`); simulation with 10,000 paths per cell, each with its own window and 7800 monitored i.i.d. returns of the same law (`results/th5_plugin_sim.csv`). All entries are ratios to the nominal `m / A`, and every ratio is capped at `A / m` (79.6 at m = 98, 20 at m = 390, 4 at m = 1950).

| Law | n | m | `m / √n` | upper bound, every stream in the class (Corollary 1 of TH4) | first-order formula (Theorem 8c) | proved lower bound, worst stream (Proposition 9) | i.i.d. returns, plug-in (simulated) | i.i.d. returns, true baseline (simulated) |
|---|---|---|---|---|---|---|---|---|
| Gaussian | 390 | 98 | 4.96 | 3.03 | 1.18 | 0.63 | 0.13 | 0.04 |
| Gaussian | 390 | 390 | 19.75 | 8.32 | 1.73 | 6.46 | 0.72 | 0.19 |
| Gaussian | 390 | 1,950 | 98.74 | 3.70 | 4.64 | 3.48 | 1.63 | 0.49 |
| Gaussian | 1,440 | 98 | 2.58 | 1.21 | 1.10 | 0.01 | 0.05 | 0.02 |
| Gaussian | 1,440 | 390 | 10.28 | 5.16 | 1.38 | 3.31 | 0.32 | 0.18 |
| Gaussian | 1,440 | 1,950 | 51.39 | 3.45 | 2.90 | 3.03 | 0.94 | 0.49 |
| Gaussian | 5,760 | 98 | 1.29 | 1.06 | 1.05 | 0.01 | 0.09 | 0.07 |
| Gaussian | 5,760 | 390 | 5.14 | 2.24 | 1.19 | 0.57 | 0.24 | 0.21 |
| Gaussian | 5,760 | 1,950 | 25.69 | 2.95 | 1.95 | 2.23 | 0.59 | 0.49 |
| Gaussian | 23,040 | 98 | 0.65 | 1.03 | 1.02 | 0.00 | 0.07 | 0.06 |
| Gaussian | 23,040 | 390 | 2.57 | 1.19 | 1.09 | 0.05 | 0.20 | 0.17 |
| Gaussian | 23,040 | 1,950 | 12.85 | 2.32 | 1.47 | 1.34 | 0.53 | 0.50 |
| Gaussian | 92,160 | 98 | 0.32 | 1.01 | 1.01 | 0.00 | not run | not run |
| Gaussian | 92,160 | 390 | 1.28 | 1.06 | 1.05 | 0.04 | not run | not run |
| Gaussian | 92,160 | 1,950 | 6.42 | 1.74 | 1.24 | 0.61 | not run | not run |
| Student t3 | 390 | 390 | 19.75 | 10.79 | 1.93 | 8.88 | 1.03 | 0.33 |
| Student t3 | 390 | 1,950 | 98.74 | 3.80 | 5.65 | 3.66 | 1.69 | 0.58 |
| Student t3 | 23,040 | 390 | 2.57 | 1.28 | 1.12 | 0.07 | 0.32 | 0.30 |
| Student t3 | 23,040 | 1,950 | 12.85 | 2.68 | 1.61 | 1.80 | 0.57 | 0.54 |

* Where `m / √n < 1` (m = 98, n ≥ 11,520), the upper bound is at most 4.1% (Gaussian) and 5.3% (t3) above nominal, and the first-order formula `1 + C_F (m + 1) / √n` is within 1% of it.
* Where `m / √n` is large, the lower bound for the worst stream (a window average of proved per-window bounds, Monte Carlo) can be a large multiple of nominal (6.5 at n = m = 390 for Gaussian returns, 8.9 for t3) and is 25% to 96% of the upper bound (24 cells with `m / √n > 5`, median 67%). It is above nominal in 19 of these 24 cells; the five others (0.57 to 0.95 of nominal) have `m / √n` between 5.1 and 9.1. Run through the detector on 4000 fresh windows, that stream alarmed within one day of a 390-step window in 36.5% of windows (95% interval 35.0% to 38.0%), between the per-window lower bound averaged over the same windows (33.8%) and Corollary 1's upper bound for them (42.4%) (`results/th5_holdbold.csv`; all four cells checked are sandwiched).
* The distribution-free bound (Theorem 8b) is informative only for horizons far below √n: a factor of 1.15 at m = 39 and 2.3 at m = 98 for n = 10⁶, and 1.021 at m = 98 and 1.15 at m = 390 for n = 10⁸ (`results/th5_distribution_free.csv`; a first version of the integration under-evaluated the bound for n ≥ 10⁶ and reported 1.013 and 1.07).
* `√n sd(d)` is within 1.3% of its limit `σ_Y` at every n from 390 (`results/th5_plugin_drift.csv`), so the normal approximation of Theorem 8c is accurate at these window lengths.

**Independent monitored returns, false alarm rates in percent** (`results/th5_plugin_sim.csv`; Wilson 95% half-widths at most 1.0 point).

| Law | n | plug-in, one day (nominal 5%) | plug-in, five days (25%) | z = 2.576, one day | z = 2.576, five days | cs, one day | cs, five days | true baseline, one day | true baseline, five days |
|---|---|---|---|---|---|---|---|---|---|
| Gaussian | 390 | 3.61 | 40.77 | 0.03 | 0.42 | 0.00 | 0.02 | 0.94 | 12.15 |
| Gaussian | 1,440 | 1.62 | 23.40 | 0.13 | 0.96 | 0.02 | 0.07 | 0.92 | 12.36 |
| Gaussian | 5,760 | 1.19 | 14.70 | 0.31 | 2.66 | 0.14 | 1.07 | 1.06 | 12.34 |
| Gaussian | 23,040 | 1.02 | 13.33 | 0.46 | 5.93 | 0.32 | 3.96 | 0.85 | 12.51 |
| Student t3 | 390 | 5.14 | 42.34 | 0.06 | 0.55 | 0.01 | 0.01 | 1.64 | 14.50 |
| Student t3 | 1,440 | 2.26 | 23.81 | 0.15 | 1.53 | 0.00 | 0.25 | 1.60 | 13.96 |
| Student t3 | 5,760 | 1.73 | 16.96 | 0.49 | 4.00 | 0.30 | 2.07 | 1.57 | 14.96 |
| Student t3 | 23,040 | 1.61 | 14.36 | 0.75 | 6.62 | 0.63 | 4.78 | 1.48 | 13.56 |

* For i.i.d. returns the true baseline is itself well below nominal (0.85% to 1.64% at one day, 12.2% to 15.0% at five days), so the plug-in's excess appears first as a ratio to the true baseline. At five days that ratio is 3.4, 1.9, 1.5, 1.19, 1.13 and 1.07 for Gaussian returns as n goes from 390 to 23,040 (2.9 to 1.06 for t3), although `m / √n` is still 13 at n = 23,040. The worst stream in the class is much worse than i.i.d. returns; the comparison with the heuristic scale `log A / √n` of the theory note (section 4.6) was made after seeing the results [post hoc].
* The plug-in is not significantly above nominal at one day in any cell (largest 5.14% for t3 at n = 390, interval 4.72% to 5.59%), but it is far above at five days after a one-day window (40.8% and 42.3% against 25%).
* No simulated rate is above its Corollary 1 bound (all 288 rows). With z = 2.576 the inflated plug-in's bound is at most 1.061 times nominal in every cell, and the share of windows outside the shifted baselines is 0.009 to 0.013 against `2 Φ̄(2.576) = 0.01` (Proposition 11b). The inflated plug-in and the cs baseline are conservative for i.i.d. returns (at most 0.75% and 0.63% at one day); their price is a later alarm after a change, which TH5 does not measure.

### Q-TH5d. The decrease bet at its boundary

`results/th5_boundary.csv`, 4000 paths of 6000 steps, single bet at threshold A = 7800.

| Law | sd ratio | `q_B` | `ε*` | `λ* / λ_max` | delay at `ε*` | at ε = 10⁻⁶ | growth at ε = 10⁻⁶ (`KL_inf`) |
|---|---|---|---|---|---|---|---|
| Gaussian | 0.5 | 2.0 × 10⁻⁹ | 2.4 × 10⁻⁸ | 1 - 2.4 × 10⁻⁸ | 74.1 | 74.1 | 0.0881 (0.0881) |
| Gaussian | 0.8 | 1.8 × 10⁻⁴ | 8.3 × 10⁻³ | 0.9917 | 156.4 | 155.8 | 0.0332 (0.0346) |
| Student t3 | 0.5 | 6.5 × 10⁻³ | 0.069 | 0.931 | 69.7 | 71.3 | 0.0249 (0.0913) |
| Student t3 | 0.67 | 1.5 × 10⁻² | 0.223 | 0.777 | 122.8 | 139.9 | -0.1260 (0.0450) |
| Student t3 | 0.8 | 2.4 × 10⁻² | 0.455 | 0.545 | 258.8 | 396.7 | -0.2817 (0.0172) |
| jumps | 0.8 | 3.7 × 10⁻³ | 0.123 | 0.877 | 183.5 | 184.7 | -0.0127 (0.0273) |

Where the growth is negative the delay stays finite: the Shiryaev-Roberts statistic restarts after each return at the truncation level and alarms on a long enough run without one (Proposition 12e gives an ε-free bound, which is loose).

### Reading

The guarantee of Theorem B has a price that no procedure can avoid, and the e-detectors pay at most about twice the minimum. On the hardest path for the decrease side, a complete collapse, every deterministic procedure with the guarantee needs 45 steps at c = 3 (a randomised one at least 43.4 on average), the best e-detectors 58 and the default grid 73. Across realistic post-change laws one side of the default grid is at most 1.7 to 2.2 times slower than the best valid deterministic procedure could be on average (at most 1.8 to 2.2 against randomised ones), and the best single bet at most 1.3 to 1.8 times. These are upper bounds on the inefficiency, because the lower bound need not be attained. So at least half of the delay is the problem's, through the truncation level B, for the best single bet in all 18 laws and for the default grid in 16 of them (at least 46% in every law: the grid's ratio is 2.16 at Gaussian sd ratio 0.8 and 2.03 for jumps at 0.8); the rest is an unknown split between the detector and the slack in the bound.

For the plug-in baseline the relevant comparison is the horizon against the square root of the calibration length. For every stream in `C(W)`, the bound of Theorem B is inflated by at most `1 + C (m + 1) / √n` to first order when m is small compared with √n. The scale √n is sharp: for `m ≫ √n` (and `m ≫ log A`) some stream in `C(W)` breaks the nominal level by an unbounded factor. The constant C is not shown to be sharp. With minute data, a one-day horizon equals √n at n = 152,100 returns (390 trading days), which is the boundary itself, and is short compared with √n only for calibration windows of order 10⁶ returns or more. A margin of `z = Φ⁻¹(1 - δ / 2)` standard errors restores `m / A + δ` asymptotically, and the confidence-sequence baseline restores it exactly at a wider margin.

### Caveats

* The floors are for the change at the first step (no pre-change data), which is the worst case in Lorden's sense; they do not bound the conditional delay after a pre-change history (E4's delays, which used one, are smaller than E11's).
* The Neyman-Pearson parts of the delay bound are Monte Carlo estimates (importance weights from 10,000 paths per law), not certified values; batch-means standard errors are in `results/th5_floors.csv`. The pathwise floor is exact on each path.
* The delay floors of Q-TH5b are for one detector side and one-sided validity; the bound for the two-sided guarantee over `H=` is smaller and was not computed.
* The bounds concern rules valid over the full nonparametric null. A procedure valid only over a smaller model (for example a fixed return shape) is not covered.
* The plug-in lower bound uses an adversarial stream that knows the drift of the window. For i.i.d. monitored returns the boundary is not proved; a heuristic suggests the scale `log A / √n` (theory note, section 4.6), and the simulation of Q-TH5c shows a much smaller excess than the worst case.
* The distribution-free constant (1.79) is 48 times the Gaussian one (0.037); its bound is informative only for `m` well below `√n / 20`.
* Everything is synthetic and for the default grid with c = 3 and A = 7800, except where a table varies c or A.
