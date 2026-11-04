## TH4. Robustness of the guarantee: how much can a small drift of the conditional mean inflate the false alarm probability?

Theory `docs/theory_robustness.md` (statements, proofs and the full tables); code `src/avmonitor/robustness.py`; script `experiments/th4_robustness.py`; outputs `results/th4_*.csv` and `results/figures/th4_sup_bound.png`, `results/figures/th4_e16_comparison.png`; tests `tests/test_robustness.py`; iteration notes `docs/iterations/TH4.md`. The theorems are proved in the theory note; all numerical results are **[SYNTHETIC]**, with seeds fixed in the script.

**Question.** Theorem B proves `P(false alarm within m steps) ≤ m / A` for the default two-sided detector when the conditional mean of the clipped squared return is exactly the baseline m0. E16 mapped where the bound fails when that condition fails: GARCH clustering, a plug-in baseline, price ticks. What can be proved when the condition holds only approximately, `E[X_t | past] ≤ m0 (1 + δ_t)` on the increase side and `≥ m0 (1 - δ_t)` on the decrease side? How fast does the guarantee degrade with the drift δ? Which forms of drift are harmless? What does the plug-in baseline cost? And how do the bounds compare with E16's measured rates?

**Set-up.** The default detector of Theorem B: six increase bets λ = 0.02 to 0.5, five decrease bets 0.1 to 0.95 times 1 / (B - 1), Shiryaev-Roberts statistics, each side alarming at 2A, A = 7800, so the bound is 5% for one 390-step day and 25% for five days. The drift `δ_t ≥ 0` may change over time but must be known one step ahead. Its smallest values are the positive and negative parts of `E[X_t | past] / m0 - 1`.

**Summary.**
* **A bounded drift inflates the bound by a closed-form factor (Theorem 1, proved).** If the drift never exceeds δ, then `P(T ≤ m) ≤ (m / A) I_m(δ)` with `I_m(δ) = ½ [Σ_k w_k G_m(λ_k δ↑) + Σ_l v_l G_m(φ_l δ↓)]` and `G_m(x) = (1 + x) ((1 + x)^m - 1) / (m x)`. The factor is about `1 + (m + 1)(λ̄ δ↑ + φ̄ δ↓) / 4` while `λ_max δ m` is small. Once `λ_max δ m` is large it grows like `(1 + λ_max δ)^m`. At one day a drift of 0.1% on both sides gives 5.14% instead of 5%, 0.5% gives 5.80% and 1% gives 7.06%. **Design guide:** the one-day bound reaches 10% at a drift of 1.63% and 25% at 2.66% (upward drift only, a mean between m0 and `m0 (1 + δ)`: 1.72% and 2.70%; downward drift only: 6.36% and 10.51%). Over five days a drift of 0.5% on both sides, or upward, already gives 1. The bound cannot be lowered where it reaches 1, because a deterministic stream at `m0 (1 + δ)` attains it at its alarm step. Below that, in one day, the strongest adaptive stream found reaches 0.33 to 0.63 of it, and whether the bound is attained there is open.
* **A drift that averages out is cheap, but only a pathwise bound can see it (Theorems 2 and 3, proved).** What matters is the largest product of `(1 + bet × drift)` over a window. A calm spell resets it. For any drift known one step ahead (Theorem 2), `P(T ≤ m) ≤ K m / A + P(max window product > K)`. When the drift does not respond to the data (Theorem 3, with the largest window product over the whole horizon) the bound is the detector run on the drift itself. For a ±10% switching drift that flips at every step, Theorem 3 gives 15.3% at one day, where the sup form gives 1 (measured 0.73%). **What cannot be said:** when the drift responds to the data, Theorem 3's formula is not a bound (a stream where it gives 10.2% at one day against an alarm probability of 11.1%). No bound of the form `C(mean drift) × m / A` exists (a seven-step burst decided before monitoring). A bound in the mean drift exists (Proposition 4), but its extra term is weighted by A and is larger than the whole one-day budget once the mean upward drift exceeds 0.026%.
* **Plug-in baseline (Corollary 1, proved).** The baseline error is a drift fixed by the calibration window, so `P(T ≤ m) ≤ E_W[min{1, (m / A) I_m(δ(W))}]`. For Gaussian returns and 390 calibration steps this is 41.8% at one day and 92.4% at five days. E16 measured 4.0% and 41.8% for i.i.d. returns, and 32.5% and 85.8% for its noise-free shape adversary, so the bound is within a factor of 1.3 of that stream. The one-day bound is 7.3% at 11,520 calibration steps and 6.0% at 23,040. This is consistent with Proposition A.6 (E14): the constant in front of m / A has to grow with m.
* **Against E16, every proved bound is above every measured rate where its assumption holds.** That covers 128 GARCH cells, 36 plug-in cells and 80 decile cells, plus 112 new constant-drift and switching cells. For GARCH the sup form gives 1 everywhere, because the conditional mean comes arbitrarily close to `B m0`. The proved pathwise bound is informative in 22 of 64 cells at one day (6 to 29 times the measured rate) and 6 at five days, all with a small drift. Theorem 3 used as a heuristic is 1.2 to 8.3 times the measured rate and was never below it.

### Q-TH4a. The sup form: design guide and tightness

Inflation factor `I_m` for a drift of δ on both sides (`results/th4_sup_inflation.csv`; the bound is `min(1, (m / A) I_m)`):

| drift δ | one day | five days | first horizon at which the bound reaches 1 |
|---|---|---|---|
| 0.10% | 1.027 | 1.160 | 4,880 |
| 0.50% | 1.160 | 4.239 | 1,919 |
| 1.00% | 1.412 | 158 | 1,140 |
| 2.00% | 2.638 | 1.2 × 10^6 | 654 |
| 3.00% | 7.452 | 1.2 × 10^10 | 468 |

Drift at which the one-day bound reaches 10% and 25% (`results/th4_design.csv`):

| | 10% in one day | 25% in one day |
|---|---|---|
| both sides drift by δ | 1.63% | 2.66% |
| upward drift only (mean between m0 and m0 (1 + δ)) | 1.72% | 2.70% |
| downward drift only (mean between m0 (1 - δ) and m0) | 6.36% | 10.51% |
| both sides, 1440-step day with A = 28,800 | 0.44% | 0.72% |

**Tightness** (`results/th4_constant_drift.csv`, `results/th4_sharpness.csv`, figure `th4_sup_bound.png`). The setting is 80 cells with a known constant drift, 10,000 paths per i.i.d. law and 20,000 for the adaptive stream. The adaptive stream holds Y at its mean for a searched number of steps, then plays bold as in E16. No rate is above the bound. The deterministic stream alarms within one day from an upward drift of 3.756%, and the bound reaches 1 at 3.740%, so Theorem 1 is attained where it reaches 1. Below that drift the deterministic stream does not alarm within the day (at the 10% and 25% design drifts it alarms at steps 515 to 795), the adaptive stream is the strongest found, and the worst case is open. At an upward drift of 2% it gives 6.6% in one day against a bound of 12.4%, and at 3% it gives 22.7% against 35.9%. Noisy i.i.d. streams barely notice a small drift within a day: the two-point law {0, B} gives 2.0% at a drift of 0.25% and 2.4% at 5%.

### Q-TH4b. A drift that averages out

The volatility switches between two levels with truncated means `m0 (1 ± ε)`, driven by a Markov chain independent of the returns. 6,000 paths per cell (`results/th4_switching.csv`; 16 cells in the CSV file).

| ε | stay | one day: measured | sup form | Theorem 3 | five days: measured | sup form | Theorem 3 |
|---|---|---|---|---|---|---|---|
| 2% | 0.5 | 0.98% | 13.2% | 6.3% | 12.0% | 1 | 42.3% |
| 5% | 0.5 | 1.20% | 1 | 8.9% | 13.0% | 1 | 73.5% |
| 10% | 0.5 | 0.73% | 1 | 15.3% | 12.2% | 1 | 96.0% |
| 10% | 0.99 | 3.28% | 1 | 80.6% | 32.9% | 1 | 1 |
| 10% | 0.999 | 5.25% | 1 | 66.0% | 67.0% | 1 | 1 |

The pathwise bound is below 1 in all 16 cells at one day, against 8 for the sup form. Persistence turns the same drift size into a constant drift: at ε = 10% the one-day rate rises from 0.73% to 5.25% as the half-life grows from 1 to 693 steps.

### Q-TH4c. GARCH clustering (E16 Q16d): bound against measured rate

These are E16's 64 cells, on 6,000 fresh paths each with the exact drift (`results/th4_garch.csv`, figure `th4_e16_comparison.png`). The drift responds to the data, so only Theorem 2 and Propositions 2b and 4 are proved; "best proved" is the smallest of the three.

| | one day | five days |
|---|---|---|
| cells above the nominal bound (E16) | 43 | 51 |
| sup form below 1 | 0 | 0 |
| best proved bound below 1 | 22 | 6 |
| best proved / measured, where below 1 | 6.4 to 29 | 2.7 to 6.3 |
| Theorem 3 (heuristic) below 1 | 44 | 19 |
| Theorem 3 (heuristic) / measured, where below 1 | 1.2 to 8.3 | 2.0 to 4.9 |
| any bound or heuristic below the measured rate | 0 | 0 |

Examples at one day (measured, then best proved, then Theorem 3): α + β = 0.5, α = 0.002 (mean drift 0.23%): 0.97%, 6.2%, 5.2%. α + β = 0.9, α = 0.005 (1.23%): 1.53%, 25.0%, 6.1%. α + β = 0.98, α = 0.005 (2.75%): 3.95%, 79.0%, 15.3%. α + β = 0.999, α = 0.002 (4.83%, above the nominal bound): 5.58%, 1, 36.9%.

### Q-TH4d. Plug-in baseline (E16 Q16e)

Corollary 1 is evaluated on 100,000 fresh calibration windows per cell (`results/th4_plugin.csv`, `results/th4_plugin_deciles.csv`, `results/th4_plugin_caln.csv`).

| calibration | monitored stream | one day: E16 measured | bound | five days: E16 measured | bound |
|---|---|---|---|---|---|
| 390 | Gaussian, i.i.d. | 4.0% | 41.8% | 41.8% | 92.4% |
| 390 | shape adversary, q = 0 (noise free) | 32.5% | 41.8% | 85.8% | 92.4% |
| 1440 | Gaussian, i.i.d. | 1.6% | 26.1% | 23.1% | 86.1% |
| 1440 | shape adversary, q = 0 (noise free) | 16.1% | 26.1% | 73.1% | 86.1% |

Calibration length for the Gaussian plug-in (20,000 windows each): the one-day bound is 42.1% (390 steps), 25.7% (1,440), 17.9% (2,880), 11.1% (5,760), 7.3% (11,520) and 6.0% (23,040). The five-day bound is 92.4% down to 57.7%. (The 390-step value, 42.1% on 20,000 windows, is the 41.8% of the table on other windows.)

### Reading

The guarantee of Theorem B degrades smoothly, not abruptly. A drift of the conditional mean of a few tenths of a percent costs a few percent of the bound over a day. The cost then compounds: the one-day bound doubles at a drift of about 1.6% and reaches 1 at about 3.7% (upward drift). Over five days a drift of 0.5% exhausts it. This gives a concrete scale for "the null holds approximately". Where the bound reaches 1 the scale is sharp: a deterministic stream at m0 (1 + δ) alarms within the day from an upward drift of 3.756%, where the bound reaches 1 at 3.740%. At the design drifts the deterministic stream does not alarm within the day, the strongest stream found reaches 0.33 to 0.63 of the bound, and whether the 10% and 25% values are attained is open. For noisy streams it is conservative: with a constant upward drift of 2%, i.i.d. Gaussian returns alarm on 1.5% of days and the two-point law {0, B} on 1.9%, against a bound of 12.4%. Drift that alternates is far cheaper than its size suggests, but only a pathwise bound can show it. When the drift responds to the data, as under GARCH, the proved pathwise bound has to cover the case where a high drift comes exactly when the wealth is high, which makes it loose. The plug-in baseline is a persistent drift fixed before monitoring starts, the case where nothing averages out, and there Corollary 1 is close to the worst stream E16 found. In practice the drift must be kept small and short-lived relative to the baseline. That is what E8 (a confidence-sequence baseline), E10 (a forecast-relative null) and E14 (a rank detector) are for.

### Caveats

* The proved bounds are statements about all streams with a given drift. Their numerical values for GARCH, switching and the plug-in are Monte Carlo estimates of expectations over the drift, and the threshold K of Theorem 2 is chosen on the same paths (optimistic by at most 1.75 points with 95% probability, by the Dvoretzky-Kiefer-Wolfowitz inequality; adding this changes no count).
* Theorem 3 for GARCH and Theorem 1 at E16's mean drift are heuristics. They held in every cell here, but neither is a bound in general: section 3.4 of the theory note gives, for each, a stream with a drift that responds to the data or is decided before monitoring where it is below the alarm probability at one day.
* Price ticks are not covered. The theorems need the drift given the whole past, which for rounded prices is a filtering problem, and E16 stored the drift given the previous return only. Theorem 1 at that drift (`results/th4_ticks.csv`) is a heuristic that can be too small.
* The worst case over streams with a given constant drift, at horizons shorter than the deterministic alarm step, is not known. The adaptive stream is a lower bound on it.
* All statements are for band 0. With the monitor's null band b, Theorem 1 holds with the drifts reduced to `((δ↑ - b) / (1 + b))⁺` and `((δ↓ - b) / (1 - b))⁺` (remark after Theorem 1 in the theory note).
* Everything is synthetic and for the default grid with A = 7800 (one design row uses A = 28,800).
