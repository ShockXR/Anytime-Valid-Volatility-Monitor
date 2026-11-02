# Mathematical notes

This document collects the definitions, derivations and proofs the prototype relies on. Each section ends with what is proved, what is checked only numerically, and what is assumed. A plain-language version is in `docs/paper_review.md` and `PITCH_BAPAT.md`.

## 1. Set-up and notation

* Log price on `[0, T]`: `dp_t = μ_t dt + σ_t dW_t + dJ_t` (Bapat 2026, eq. 8).
* Equally spaced returns `r_i = p_{iΔ} - p_{(i-1)Δ}`, i = 1, ..., n, with `Δ = T / n`.
* Integrated variance `IV = ∫₀ᵀ σ_t² dt` and integrated quarticity `IQ = ∫₀ᵀ σ_t⁴ dt`.
* Spot variance `σ_t²` and spot quarticity `σ_t⁴`.
* For a standard normal Z: `φ`, `Φ` are its density and distribution function, and `μ_p = E|Z|^p = 2^{p/2} Γ((p+1)/2) / √π`.

Why quarticity matters for variance: with locally constant volatility and no jumps, `r_i ≈ σ_i √Δ Z_i`, so `Var(r_i²) = 2 σ_i⁴ Δ²`. The variance of realized variance `Σ r_i²` is therefore `2 Δ Σ σ_i⁴ Δ ≈ 2 Δ IQ`. This is the classical result of Barndorff-Nielsen and Shephard (2002): `√n (RV - IV) → N(0, 2 T IQ)`. Any test or confidence interval for variance needs an estimate of IQ.

## 2. Normalised local estimators

Let `w_0, ..., w_h` be one-sided kernel weights with `w_k ∝ K(-k/h)` (lag k = 0 is the newest return) and `Σ_k w_k = 1`. Define

* local variance `V̂_t = (1/Δ) Σ_k w_k r_{t-k}²`,
* local quarticity `Q̂_t = (1 / (3 Δ²)) Σ_k w_k r_{t-k}⁴`.

If σ is constant over the window and there are no jumps, `E[V̂_t] = σ²` and `E[Q̂_t] = σ⁴`, for every h and every kernel. Also `Var(V̂_t) = 2 σ⁴ Σ_k w_k²`, which for the uniform kernel is `2 σ⁴ / (h + 1)`. So the local estimator has rate √h, not √n.

An integrated estimate is `IQ̂ = T · mean_t Q̂_t`. Section 4 shows this is essentially realized quarticity.

**The printed estimator.** Eq. (6) uses raw weights `K((i - t)/h)` that sum to `W(h) = Σ_{k=0}^h K(k/h)`, and has no `1/Δ²`. Its mean is `W(h) σ⁴ Δ²`, that is `W(h) / n²` times `σ⁴` when T = 1. For the Gaussian kernel `W(h) ≈ 0.341 h + 0.32`.

*Status:* proved (elementary). Checked numerically in `tests/test_robust.py`.

## 3. Truncated moments

Clip at `a` standard deviations: `Z_a = clip(Z, -a, a)`. Integration by parts gives

* `E[Z² 1{|Z| < a}] = (2Φ(a) - 1) - 2 a φ(a)`,
* `E[Z⁴ 1{|Z| < a}] = 3 E[Z² 1{|Z| < a}] - 2 a³ φ(a)`.

So the share of the second and fourth moments that survives clipping is

* `κ₂(a) = E[Z_a²] = (2Φ(a) - 1) - 2 a φ(a) + 2 a² (1 - Φ(a))`,
* `κ₄(a) = E[Z_a⁴] / 3 = [3 ((2Φ(a) - 1) - 2 a φ(a)) - 2 a³ φ(a) + 2 a⁴ (1 - Φ(a))] / 3`.

With `τ = c · MAD`, the cut-off in standard deviations is `a = c · Φ⁻¹(3/4) ≈ 0.6745 c` for the raw MAD and `a = c` for the scaled MAD. For c = 3:

| | a | `P(|Z| > a)` | κ₂ | κ₄ |
|---|---|---|---|---|
| raw MAD | 2.023 | 0.0430 | 0.925 | 0.705 |
| scaled MAD | 3.000 | 0.0027 | 0.995 | 0.964 |

*Status:* proved. Checked against 2 million draws in `test_truncated_moment_factor_matches_monte_carlo`.

### 3.1 Why the truncation bias does not vanish

Let `s_n = σ √Δ` be the scale of a diffusive return. `MAD(r)` is a consistent estimate of `0.6745 s_n` (for constant σ), so `τ_n / s_n → 0.6745 c`, a constant. The clipped fraction of returns therefore converges to `P(|Z| > a) > 0`, and the truncated realized quarticity converges in probability to `κ₄(a) · IQ`, not IQ.

For the truncation bias to vanish, the threshold must satisfy `τ_n / s_n → ∞`. Mancini (2009) uses `τ_n = c Δ^ϖ` with `0 < ϖ < 1/2`. To still remove jumps, the threshold must also go to 0. A practical alternative at fixed n is to divide by `κ₄(a)`. That correction is exact only if the diffusive returns are locally normal.

*Status:* proved for constant volatility. With stochastic volatility and one global MAD, the bias depends on the volatility path (numerical evidence in E2b).

## 4. The averaged kernel estimator equals truncated RQ up to edge effects

Let `x_i = r̃_i⁴` and let the raw weight at lag k be `v_k`. Then

`(1 / (n - h)) Σ_{t=h+1}^{n} Σ_{k=0}^{h} v_k x_{t-k} = (1 / (n - h)) Σ_{i=1}^{n} c_i x_i`, where `c_i = Σ_{k : h+1 ≤ i+k ≤ n} v_k`.

For `h < i ≤ n - h`, every k qualifies, so `c_i = W(h)`. Only the first h and last h returns get partial weights. So eq. (16) equals `W(h) / (3 (n - h)) Σ_i r̃_i⁴` plus an edge term over 2h returns. In other words, it is a rescaled truncated realized quarticity. The same is true for normalised weights, with W replaced by 1.

*Status:* exact identity, checked in `test_aggregated_kernel_is_weighted_realized_quarticity`.

## 5. Anytime-valid detection of variance changes

### 5.1 The statistic and the null

Fix a truncation level τ and a baseline `m₀ > 0` from a calibration window. For each new return define

`X_t = min(r_t², τ²) ∈ [0, τ²]`, and `Y_t = X_t / m₀ ∈ [0, B]` with `B = τ² / m₀`.

Let `F_t` be the information up to time t. The nulls are

* for an increase: `H₀↑: E[Y_t | F_{t-1}] ≤ 1` for all t,
* for a decrease: `H₀↓: E[Y_t | F_{t-1}] ≥ 1` for all t.

These are nonparametric. They restrict only one conditional moment of the truncated squared return and allow any shape, any dependence and any jump process that keeps the moment on the null side.

**Link to volatility.** If `r_t = σ_t ε_t` and ε has a fixed distribution, then `min(σ² ε², τ²)` is nondecreasing in σ for every ε. So `E[Y_t | F_{t-1}]` is nondecreasing in `σ_t`, and `H₀↑` holds whenever `σ_t` stays at or below its calibration level. This monotonicity holds for every innovation distribution. Only the *size* of the change in `E[Y]` for a given change in σ depends on that distribution.

### 5.2 Bets and test supermartingales

For an increase, take a predictable bet `λ_t ∈ [0, 1)` and the factor `f_t = 1 + λ_t (Y_t - 1)`. Since `Y_t ≥ 0`, `f_t ≥ 1 - λ_t > 0`. Under `H₀↑`, `E[f_t | F_{t-1}] = 1 + λ_t (E[Y_t | F_{t-1}] - 1) ≤ 1`.

For a decrease, take `f_t = 1 + λ_t (1 - Y_t)` with `λ_t ∈ [0, 1/(B - 1))`. Since `Y_t ≤ B`, `f_t ≥ 1 - λ_t (B - 1) > 0`. Under `H₀↓`, `E[f_t | F_{t-1}] ≤ 1`.

In both cases `L_t^{(j)} = Π_{s=j}^{t} f_s` (with `L_{j-1}^{(j)} = 1`) is a nonnegative supermartingale with initial value 1, a *test supermartingale*, for every start time j.

### 5.3 Why a decrease detector needs a bounded statistic

**Claim.** Let `H` be the set of all distributions on `[0, ∞)` with mean at least 1. If `E: [0, ∞) → [0, ∞)` satisfies `E_P[E(Y)] ≤ 1` for every `P ∈ H` (an e-variable for H), then `E(y) ≤ 1` for every `y < 1`.

*Proof.* Fix `y < 1` and `ε ∈ (0, 1)`. Pick `z` such that the two-point law `P_ε = (1 - ε) δ_y + ε δ_z` has mean exactly 1, that is `z = (1 - (1 - ε) y) / ε`. Then `P_ε ∈ H`, so `(1 - ε) E(y) + ε E(z) ≤ 1`. Since `E ≥ 0`, this gives `E(y) ≤ 1 / (1 - ε)`. Let ε go to 0. ∎

So without an upper bound on the statistic, no nonparametric evidence against "variance has not fallen" can ever build up from small observations. Truncation supplies the bound `Y ≤ B`, and with it the bets in 5.2. This is the most direct link between Bapat's truncation and the detector: **clipping is what makes a robust detector of volatility decreases possible at all.** For increases, clipping is not needed for validity, but it bounds the effect of a single jump on the wealth (`f_t ≤ 1 + λ (B - 1)`), which is the betting version of the bounded influence of Theorem 3.

### 5.4 E-detectors and the false alarm guarantees

Following Shin, Ramdas and Rinaldo (2024), define the Shiryaev-Roberts type statistic

`M_t = Σ_{j=1}^{t} L_t^{(j)}`, which satisfies `M_t = f_t (M_{t-1} + 1)` with `M_0 = 0`

(when all chains share the same predictable bet `λ_t`), and the CUSUM type statistic `C_t = max_j L_t^{(j)}`, which satisfies `C_t = f_t max(C_{t-1}, 1)`. The detector raises an alarm at `τ_A = inf{t : M_t ≥ A}`.

**Lemma 1.** Under the null, `M_t - t` is a supermartingale.

*Proof.* `E[M_t | F_{t-1}] = E[f_t | F_{t-1}] (M_{t-1} + 1) ≤ M_{t-1} + 1`. ∎

**Theorem A (average run length).** Under the null, `E[τ_A] ≥ A`.

*Proof.* If `E[τ_A] = ∞` there is nothing to show. Otherwise, the increments of `M_t - t` are integrable on `{t ≤ τ_A}` and optional stopping applies to Lemma 1, giving `E[M_{τ_A}] ≤ E[τ_A]`. Since `M_{τ_A} ≥ A`, `A ≤ E[τ_A]`. ∎

**Theorem B (false alarm within a horizon).** Under the null, for every integer m, `P(τ_A ≤ m) ≤ m / A`.

*Proof.* For `t ≤ m`, set `S_t = M_t + (m - t)`. By Lemma 1, S is a supermartingale on `{0, ..., m}`, it is nonnegative, and `S_0 = m`. Ville's inequality gives `P(max_{t ≤ m} S_t ≥ A) ≤ m / A`. Since `M_t ≥ A` implies `S_t ≥ A`, the claim follows. ∎

**Consequences used in the code.**

* Choosing `A = 7800` for 390-step days gives at most a 5% chance of any false alarm per trading day. Choosing `A = 28800` does the same for 1440-minute crypto days.
* `C_t ≤ M_t` at all times, so the CUSUM version alarms no earlier and inherits both guarantees. It is more conservative (E4: run length above 2.4 times A on average, with runs cut off at 20,000 steps).
* An average of several M statistics (several fixed bets, or a fixed grid plus an adaptive bet) has the same recursion property. It is again an e-detector and inherits Theorems A and B.
* Two-sided monitoring: `(M↑ + M↓) / 2` is an e-detector, and `max(M↑, M↓) ≥ 2A` implies `(M↑ + M↓) / 2 ≥ A`. Running each side at level `2A` therefore keeps both guarantees for the first alarm of either side.

*Status:* proved, assuming the null holds with the true `m₀`.

### 5.5 Choosing the bet: where quarticity enters

For a bet on `Y` with conditional mean `μ` and variance `v`, the expected log-growth `E log(1 + λ (Y - 1))` is maximised approximately at

`λ* ≈ (μ - 1) / (v + (μ - 1)²)`

(the aGRAPA approximation of Waudby-Smith and Ramdas 2024). Since `E[Y²] = E[r̃⁴] / m₀²`, the variance v is a *truncated local quarticity* (up to the factor 3 and the scale `m₀²`). The adaptive bettor in `detectors.betting_edetector(adaptive_h=...)` estimates μ and `E[Y²]` with the normalised one-sided Gaussian kernel of section 2, applied to data strictly before t, so `λ_t` is predictable and validity is untouched.

*Status:* validity proved (any predictable bet works). Whether the adaptive bet improves detection is an empirical question. In E4 it did not improve on a fixed grid of bets by a meaningful margin.

### 5.6 What is not covered by the proofs

1. **Estimated baseline.** The theorems assume the true `m₀`. In practice `m₀` comes from a calibration window. If the estimate is too low, the increase detector is anti-conservative. E4d measures the effect: with 390 calibration steps the one-day false alarm rate was 5.4% (bound 5%), with 1440 or more it was at most 2.4%. A margin of one standard error removes the excess. A fully rigorous treatment would use a confidence sequence for `m₀` (Shekhar and Ramdas 2024).
2. **Volatility that drifts inside a regime.** The null asks the conditional truncated variance never to exceed the baseline. Slow stochastic volatility violates this by design, so the detector rightly alarms (E4a: about 25% per day for a log-volatility standard deviation of 0.15). The "null band" `H₀↑: E[Y_t | F_{t-1}] ≤ 1 + b` restores control at the cost of detecting only changes larger than about b (E4c).
3. **Intraday seasonality.** Real returns are deseasonalised by a predictable intraday profile (`data.intraday_scale`). Errors in that profile act like a small time-varying violation of the null and are absorbed by the band.
4. **Interpretation.** An alarm means the truncated variance moved. Mapping that to a change in σ requires a fixed innovation shape (section 5.1). A change in the tail shape alone, with σ unchanged, can also trigger an alarm.
