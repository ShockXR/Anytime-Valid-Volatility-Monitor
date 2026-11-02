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

Let `F_t` be the information up to time t: the monitored returns up to t together with the calibration data. Everything the monitor uses (τ, m₀, the bets and the deseasonalising profile) must be known at the time it is used, that is adapted to this filtration (`data.intraday_scale` with its default `start="past"` satisfies this). The nulls are

* for an increase: `H₀↑: E[Y_t | F_{t-1}] ≤ 1` for all t,
* for a decrease: `H₀↓: E[Y_t | F_{t-1}] ≥ 1` for all t.

These are nonparametric. They restrict only one conditional moment of the truncated squared return and allow any shape, any dependence and any jump process that keeps the moment on the null side. "Any dependence" is a statement about the shape of the joint law, not a licence for any data: the *conditional* mean must stay on the null side at every step. Two common features of real data violate this even when the unconditional variance is constant. Persistent volatility clustering moves the conditional mean above and below m₀ (E9, E10). Price-level microstructure noise makes returns MA(1), so after a large return the next squared return has a higher conditional mean; the raw-return monitor then loses its guarantee (8.5% to 41.5% false alarms per day for noise-to-signal ratios 0.5 to 4 in the pre-averaging alternative, `docs/report_parts/ALT1.md`), and non-overlapping pre-averaged blocks restore it.

**Link to volatility.** Suppose `r_t = σ_t ε_t`, where `σ_t` is known at time t - 1 and the innovation `ε_t` is independent of `(σ_t, F_{t-1})` with a fixed distribution. Then `min(σ² ε², τ²)` is nondecreasing in σ for every ε, so `E[Y_t | F_{t-1}]` is a nondecreasing function of `σ_t`, and `H₀↑` holds whenever `σ_t ≤ σ₀` almost surely, where σ₀ is the calibration level (and `H₀↓` whenever `σ_t ≥ σ₀`). This monotonicity holds for every innovation distribution; only the *size* of the change in `E[Y]` for a given change in σ depends on it. Without the independence (for example with a contemporaneous leverage effect, where the sign or size of ε_t is linked to σ_t), the conditional mean need not be monotone in σ.

### 5.2 Bets and test supermartingales

For an increase, take a predictable bet `λ_t ∈ [0, 1)` and the factor `f_t = 1 + λ_t (Y_t - 1)`. Since `Y_t ≥ 0`, `f_t ≥ 1 - λ_t > 0`. Under `H₀↑`, `E[f_t | F_{t-1}] = 1 + λ_t (E[Y_t | F_{t-1}] - 1) ≤ 1`.

For a decrease, take `f_t = 1 + λ_t (1 - Y_t)` with `λ_t ∈ [0, 1/(B - 1))`. Since `Y_t ≤ B`, `f_t ≥ 1 - λ_t (B - 1) > 0`. Under `H₀↓`, `E[f_t | F_{t-1}] ≤ 1`.

In both cases `L_t^{(j)} = Π_{s=j}^{t} f_s` (with `L_{j-1}^{(j)} = 1`) is a nonnegative supermartingale with initial value 1, a *test supermartingale*, for every start time j.

### 5.3 Why a decrease detector needs a bounded statistic

**Claim.** Let `H` be the set of all distributions on `[0, ∞)` with mean at least 1. If `E: [0, ∞) → [0, ∞)` satisfies `E_P[E(Y)] ≤ 1` for every `P ∈ H` (an e-variable for H), then `E(y) ≤ 1` for every `y ≥ 0`. So no e-variable for H ever exceeds 1.

*Proof.* Fix `y < 1` and `ε ∈ (0, 1)`. Pick `z` such that the two-point law `P_ε = (1 - ε) δ_y + ε δ_z` has mean exactly 1, that is `z = (1 - (1 - ε) y) / ε`. Then `P_ε ∈ H`, so `(1 - ε) E(y) + ε E(z) ≤ 1`. Since `E ≥ 0`, this gives `E(y) ≤ 1 / (1 - ε)`. Let ε go to 0. For `y ≥ 1` the point mass `δ_y` is itself in H, so `E(y) ≤ 1` directly. ∎

This is a one-sided relative of the classical impossibility result of Bahadur and Savage (1956, Annals of Mathematical Statistics 27(4), 1115 to 1122): without restrictions on the tails, the mean cannot be tested. The asymmetry matters. The *increase* null "mean at most 1" on `[0, ∞)` does admit e-variables (`E(y) = y`), which is why the increase detector also works on raw r² (E4).

So a fall in the variance itself cannot be detected by any nonparametric e-test, whatever statistic it computes. A detector of decreases must therefore change the target to a *bounded* functional of the returns. Truncation does exactly that: the truncated variance `E[min(r², τ²)]` is bounded and testable, and under the fixed-shape assumption of 5.1 it falls exactly when σ falls. This is the most direct link between Bapat's truncation and the detector. Truncation is not the only route (exceedance indicators `1{|r| > q}`, ranks or any bounded increasing transform of r² would also work; E14 studies a rank-based detector), but it is the one with a direct volatility meaning and a known bias (section 3). For increases, clipping is not needed for validity, but it bounds the effect of a single jump on the wealth (`f_t ≤ 1 + λ (B - 1)`), which is the betting version of the bounded influence of Theorem 3.

### 5.4 E-detectors and the false alarm guarantees

Following Shin, Ramdas and Rinaldo (2024), define the Shiryaev-Roberts type statistic

`M_t = Σ_{j=1}^{t} L_t^{(j)}`, which satisfies `M_t = f_t (M_{t-1} + 1)` with `M_0 = 0`

(when all chains share the same predictable bet `λ_t`), and the CUSUM type statistic `C_t = max_j L_t^{(j)}`, which satisfies `C_t = f_t max(C_{t-1}, 1)`. The detector raises an alarm at `τ_A = inf{t : M_t ≥ A}`.

**Lemma 1.** Under the null, `M_t - t` is a supermartingale.

*Proof.* `E[M_t | F_{t-1}] = E[f_t | F_{t-1}] (M_{t-1} + 1) ≤ M_{t-1} + 1`. ∎

**Theorem A (average run length).** Under the null, `E[τ_A] ≥ A`.

*Proof.* If `E[τ_A] = ∞` there is nothing to show. Otherwise write τ = τ_A. For every fixed t, the stopped process `M_{τ∧t} - (τ∧t)` is a supermartingale started at 0 (Lemma 1 and optional stopping at the bounded time τ∧t), so `E[M_{τ∧t}] ≤ E[τ∧t] ≤ E[τ]`. Since `E[τ] < ∞`, τ is finite almost surely, so `M_{τ∧t} → M_τ ≥ A` almost surely as t grows. M is nonnegative, so Fatou's lemma gives `A ≤ E[M_τ] ≤ liminf E[M_{τ∧t}] ≤ E[τ]`. ∎ (Corrected after the independent review; the earlier version appealed to integrable increments, which is not the right condition.)

**Theorem B (false alarm within a horizon).** Under the null, for every integer m, `P(τ_A ≤ m) ≤ m / A`.

*Proof.* For `t ≤ m`, set `S_t = M_t + (m - t)`. By Lemma 1, S is a supermartingale on `{0, ..., m}`, it is nonnegative, and `S_0 = m`. Ville's inequality gives `P(max_{t ≤ m} S_t ≥ A) ≤ m / A`. Since `M_t ≥ A` implies `S_t ≥ A`, the claim follows. ∎

**Consequences used in the code.**

* Choosing `A = 7800` for 390-step days gives at most a 5% chance of any false alarm per trading day. Choosing `A = 28800` does the same for 1440-minute crypto days.
* `C_t ≤ M_t` at all times, so the CUSUM version alarms no earlier and inherits both guarantees. It is more conservative (E4a: average run length 2.31 to 2.40 times A, with runs cut off at 20,000 steps).
* An average `M̄_t` of several M statistics (several fixed bets, or a fixed grid plus an adaptive bet) does not obey the recursion itself, but it satisfies `E[M̄_t | F_{t-1}] ≤ M̄_{t-1} + 1`, which is all Lemma 1 needs. It is again an e-detector and inherits Theorems A and B.
* Two-sided monitoring: under the two-sided null `(1 - b) m₀ ≤ E[X_t | F_{t-1}] ≤ (1 + b) m₀`, `(M↑ + M↓) / 2` is an e-detector, and `max(M↑, M↓) ≥ 2A` implies `(M↑ + M↓) / 2 ≥ A`. Running each side at level `2A` therefore keeps both guarantees for the first alarm of either side.

**Corollary C (restarts after alarms).** The monitor restarts after every alarm. Suppose the null holds in every segment relative to that segment's baseline (true baselines). Let `G_t` be the sum of the two-sided statistic `(M↑ + M↓)/2` at all past alarms, plus its current value, minus the number of monitored steps so far. G is a supermartingale with `G₀ = 0`, and every alarm adds at least A to the banked sum. Hence the expected number of false alarms in m monitored steps is at most `m / A`. Numerical check (E4e) [SYNTHETIC]: with A = 300, the true baseline gave 26.5 alarms per 10,000 monitored steps and the streaming monitor with plug-in recalibration on 390 steps gave 30.1, against the bound of 33.3. Section 6 (Corollary 6.6) gives the version with an estimated baseline. (Corollary C was suggested by the independent review.)

*Status:* proved, assuming the null holds with the true `m₀`.

### 5.5 Choosing the bet: where quarticity enters

For a bet on `Y` with conditional mean `μ` and variance `v`, the expected log-growth `E log(1 + λ (Y - 1))` is maximised approximately at

`λ* ≈ (μ - 1) / (v + (μ - 1)²)`

(the aGRAPA approximation of Waudby-Smith and Ramdas 2024). Since `E[Y²] = E[r̃⁴] / m₀²`, the second moment of Y (which with μ gives v) is a *truncated local quarticity* (up to the factor 3 and the scale `m₀²`). The adaptive bettor in `detectors.betting_edetector(adaptive_h=...)` estimates μ and `E[Y²]` with the normalised one-sided Gaussian kernel of section 2, applied to data strictly before t, so `λ_t` is predictable and validity is untouched.

*Status:* validity proved (any predictable bet works). Whether the adaptive bet improves detection is an empirical question. In E4 it was slightly slower than a fixed grid of bets in 14 of 18 cells (by up to 9%).

### 5.6 What is not covered by the proofs

1. **Estimated baseline.** The theorems assume the true `m₀`. In practice `m₀` comes from a calibration window. If the estimate is too low, the increase detector is anti-conservative. E4d and E4e measure the effect. The informative comparison is with the rate under the true baseline, not with the bound: with 390 calibration steps the one-day false alarm rate of the increase side was 4.6% to 5.1% against 1.6% to 1.8% with the true baseline (E4e; E8 found 4.5% to 4.9%), so the plug-in error roughly triples the rate and only the slack in Ville's inequality keeps it near the bound. With 1440 steps it was 2.0% to 2.8%. Over longer horizons the excess grows: in E10 (Q10b2), exactly i.i.d. standardised returns monitored for five days with a one-day plug-in baseline gave 44% false alarms against the 25% the theorem guarantees with the true baseline. A margin of one standard error removes most of the one-day excess. A fully rigorous treatment would use a confidence sequence for `m₀` (Shekhar and Ramdas 2024). Section 6 does this: with a confidence-sequence baseline (`baseline="cs"`) the guarantee becomes P(false alarm within m steps) ≤ δ + m / A, proved with the baseline estimated from the data.
2. **Volatility that drifts inside a regime.** The null asks the conditional truncated variance never to exceed the baseline. Slow stochastic volatility violates this by design, so the detector rightly alarms (E4a: about 25% per day for a log-volatility standard deviation of 0.15). The "null band" `H₀↑: E[Y_t | F_{t-1}] ≤ 1 + b` reduces false alarms at the cost of detecting only changes larger than about b (E4c: from about 25% to 8.3% at b = 0.25 and 1.7% at b = 0.5). It restores the guarantee only if the conditional truncated variance stays inside the band, which slow log-normal drift does not do, and persistent GARCH clustering is not absorbed by any reasonable band (E9). Monitoring relative to a volatility forecast is the principled alternative (E10).
3. **Intraday seasonality.** Real returns are deseasonalised by a predictable intraday profile (`data.intraday_scale`). Errors in that profile act like a small time-varying violation of the null. We expect a band to absorb small errors, but no experiment so far isolates this effect (E7 compares deseasonalisation choices on Nifty).
4. **Interpretation.** An alarm means the truncated variance moved. Mapping that to a change in σ requires a fixed innovation shape (section 5.1). A change in the tail shape alone, with σ unchanged, can also trigger an alarm.
5. **Missing returns.** The monitor skips a missing return. That amounts to a zero bet chosen after seeing that r_t is missing, which is valid if the null holds conditionally on the return being observed: missingness completely at random, or decided by the past (for example scheduled breaks). Missing returns must never be replaced by zeros, which look like evidence for a fall in volatility (review M6; fixed in E5 and E8). In the plug-in calibration, missing values in the window reduce the effective sample below `cal_len`.

## 6. Keeping the guarantee with an estimated baseline

Section 5 assumes the true baseline `m₀`. This section proves a guarantee that holds when the baseline is *estimated* from a calibration window, at the price of a small failure probability δ that the user chooses:

> If there is no change, P(any false alarm within m monitored steps) ≤ δ + m / A.

The idea in one sentence: replace the point estimate of `m₀` by a *confidence bound* that is valid at every sample size at once (a *confidence sequence*), so that the increase detector bets against a value that is, with probability at least 1 - δ, not below the truth, and the decrease detector against a value that is not above it. Code: `src/avmonitor/confseq.py`, `VolatilityMonitor(baseline="cs")` and `detectors.calibrate(..., baseline="cs")`. Experiment: E8 (`docs/report_parts/E8.md`).

### 6.1 Set-up

* Returns `r₁, r₂, ...` and the information up to time t, `F_t = σ(r₁, ..., r_t)`.
* **Sample split.** The first `n_τ` returns of the calibration window fix the truncation level, `τ = c × scaled MAD(r₁, ..., r_{n_τ})`. So τ is known at time `n_τ` (it is `F_{n_τ}`-measurable). For `t > n_τ` put `X_t = min(r_t², τ²) ∈ [0, τ²]` and `Z_t = X_t / τ² ∈ [0, 1]`. Both are known at time t (they are *adapted*). Section 6.4 explains why the split is needed.
* `μ_t = E[X_t | F_{t-1}]` is the conditional truncated variance and `ν_t = μ_t / τ²`.
* **End of calibration.** T is a *stopping time* with `T > n_τ`: the decision to stop at t uses only data up to t. Examples: T = cal_len, or the adaptive rule of 6.6.
* **Monitoring** starts at T + 1 and uses the same τ.

### 6.2 A confidence sequence for the mean of a bounded variable

**Lemma 6.1.** Let `a > 0` and `λ ≥ 0` with `λ a < 1`, and put `ψ_a(λ) = (-log(1 - λa) - λa) / a²` (which is ≥ 0). Then for every `ξ ≥ -a`,

`log(1 + λξ) ≥ λξ - ψ_a(λ) ξ²`, equivalently `exp(λξ - ψ_a(λ) ξ²) ≤ 1 + λξ`.

*Proof.* For `x > -1`, Taylor's formula with integral remainder for `log(1 + x)` gives `log(1 + x) - x = x² φ(x)` with `φ(x) = -∫₀¹ (1 - s) / (1 + s x)² ds`. For each `s ∈ [0, 1]` the function `x ↦ 1 / (1 + s x)²` is nonincreasing on `(-1, ∞)`, so φ is nondecreasing. Put `x = λξ ≥ -λa > -1`. Then `log(1 + λξ) - λξ = (λξ)² φ(λξ) ≥ (λξ)² φ(-λa)`, and `φ(-λa) = (log(1 - λa) + λa) / (λa)²`, so the right side equals `ξ² (log(1 - λa) + λa) / a² = -ψ_a(λ) ξ²`. (For λ = 0 both sides are 0.) ∎

With a = 1 this is the inequality of Fan, Grama and Liu (2015) that underlies the empirical Bernstein confidence sequences of Howard et al. (2021) and Waudby-Smith and Ramdas (2024). The version with a general a follows by rescaling.

**Construction.** Index the observations of the sequence by `i = 1, 2, ...` (global time `n_τ + i`). In this subsection `F_{i-1}` means the information before the i-th observation (global time `n_τ + i - 1`), and `ν_i = E[Z_i | F_{i-1}]`. Fix error levels `δ_u, δ_l ∈ (0, 1)` and a cap `c ∈ (0, 1)`. Choose, using only data before step i (*predictably*):

* a centre `ẑ_{i-1} ∈ (0, 1)`,
* one-sided ranges `a↑_i = 1 - ẑ_{i-1}` and `a↓_i = ẑ_{i-1}`, so that `ẑ_{i-1} - Z_i ≥ -a↑_i` and `Z_i - ẑ_{i-1} ≥ -a↓_i` always hold,
* bets `λ↑_i, λ↓_i > 0` with `λ↑_i a↑_i ≤ c` and `λ↓_i a↓_i ≤ c`.

Put `e_i = (Z_i - ẑ_{i-1})²`, `ψ↑_i = ψ_{a↑_i}(λ↑_i)`, `ψ↓_i = ψ_{a↓_i}(λ↓_i)` and

* `U_t = [Σ_{i≤t} λ↑_i Z_i + log(1/δ_u) + Σ_{i≤t} ψ↑_i e_i] / Σ_{i≤t} λ↑_i`,
* `L_t = [Σ_{i≤t} λ↓_i Z_i - log(1/δ_l) - Σ_{i≤t} ψ↓_i e_i] / Σ_{i≤t} λ↓_i`.

They estimate bet-weighted averages of the conditional means:

`ν̄↑_t = Σ_{i≤t} λ↑_i ν_i / Σ_{i≤t} λ↑_i` and `ν̄↓_t = Σ_{i≤t} λ↓_i ν_i / Σ_{i≤t} λ↓_i`.

**Theorem 6.2 (time-uniform coverage).** Without any assumption on the dependence between observations or on how `ν_i` varies over time,

`P(ν̄↑_t ≥ U_t for some t ≥ 1) ≤ δ_u` and `P(ν̄↓_t ≤ L_t for some t ≥ 1) ≤ δ_l`.

*Proof.* Upper bound. Let `K_0 = 1` and `K_t = Π_{i≤t} exp(λ↑_i (ν_i - Z_i) - ψ↑_i e_i)`. Write `ν_i - Z_i = ξ_i + (ν_i - ẑ_{i-1})` with `ξ_i = ẑ_{i-1} - Z_i ≥ -a↑_i`, and note `e_i = ξ_i²`. Lemma 6.1 gives

`E[K_t / K_{t-1} | F_{t-1}] ≤ E[1 + λ↑_t ξ_t | F_{t-1}] · exp(λ↑_t (ν_t - ẑ_{t-1})) = (1 - y) e^y ≤ 1`, where `y = λ↑_t (ν_t - ẑ_{t-1})`,

because `E[ξ_t | F_{t-1}] = ẑ_{t-1} - ν_t` and `(1 - y) e^y ≤ 1` for every real y (the maximum is at y = 0). So K is a nonnegative supermartingale with `K_0 = 1`. Ville's inequality gives `P(K_t ≥ 1/δ_u for some t) ≤ δ_u`. Taking logarithms, `K_t ≥ 1/δ_u` is the event `Σ λ↑_i ν_i ≥ Σ λ↑_i Z_i + Σ ψ↑_i e_i + log(1/δ_u)`, that is `ν̄↑_t ≥ U_t`.
Lower bound. The same argument with `K'_t = Π_{i≤t} exp(λ↓_i (Z_i - ν_i) - ψ↓_i e_i)`, `ξ_i = Z_i - ẑ_{i-1} ≥ -a↓_i`, and `(1 + y) e^{-y} ≤ 1`. ∎

**Corollary 6.3 (a constant level).** If `ν_i ≥ ν*` for every i, for a constant ν* (or one known when the sequence starts), then `P(U_t ≤ ν* for some t) ≤ δ_u`. If `ν_i ≤ ν**` for every i, then `P(L_t ≥ ν** for some t) ≤ δ_l`. In particular a constant conditional mean ν satisfies `P(ν ∉ (L_t, U_t) for some t) ≤ δ_u + δ_l`.

*Proof.* `ν̄↑_t` is an average of the `ν_i`, so `ν̄↑_t ≥ ν*`, and `U_t ≤ ν*` implies `ν̄↑_t ≥ U_t`. The lower case is symmetric. ∎

Since `0 ≤ ν̄ ≤ 1`, the clipped bounds `min(U_t, 1)` and `max(L_t, 0)` still satisfy `ν̄↑_t ≤ min(U_t, 1)` and `ν̄↓_t ≥ max(L_t, 0)` on the same events. In the units of X: `Ū_t = τ² min(U_t, 1)`, `L̄_t = τ² max(L_t, 0)`, `μ̄↑_t = τ² ν̄↑_t` and `μ̄↓_t = τ² ν̄↓_t`.

**Remarks.**

1. *Relation to the literature.* With `a↑ = a↓ = 1` these are the predictable plug-in empirical Bernstein (PrPl-EB) bounds of Waudby-Smith and Ramdas (2024, Theorem 2), written with `e_i (-log(1 - λ) - λ)` in place of their `v_i ψ_E(λ)`. Using the one-sided ranges lets the lower bound use bets up to `c / ẑ` instead of c. For a truncated squared normal return `ẑ ≈ 0.11`, so the cap rises from 0.75 to about 7, and at the sample sizes used here it no longer binds for the lower bound. The upper bound gains a little (`c / (1 - ẑ)`). The weighted-average form for time-varying means comes from the same proof. We prove it here and make no claim that it is new.
2. *Any predictable choice is valid.* The code uses the choices of Waudby-Smith and Ramdas: `ẑ_t = (1/2 + Σ_{i≤t} Z_i) / (t + 1)`, `σ̂²_t = (1/4 + Σ_{i≤t} (Z_i - ẑ_i)²) / (t + 1)`, and `λ_t = min(√(2 log(1/δ) / (d_t σ̂²_{t-1})), c / a_t)` with c = 0.75 and either `d_t = n_plan` (the bound is tightest near a planned sample size) or `d_t = t log(1 + t)` (their time-uniform schedule). Validity does not depend on these choices; they only decide where the bound is tight.
3. *Width.* When the cap does not bind, the half-width near `n_plan` is about `√(2 σ_Z² log(1/δ) / n)`, the empirical Bernstein rate. The truncated squared return is a hard case: its standard deviation is about 1.4 times its mean, and its range is about 6.5 standard deviations (normal returns, c = 3). E8 (Q8d) compares the bound with the best that any betting bound can achieve for the true law. The comparison is in the report; the summary is that most of the width is the price of a finite-sample guarantee for this skewed variable, not a weakness of this particular bound.

### 6.3 The false alarm guarantee with an estimated baseline

**Monitoring rule.** At time T set `m↑ = (1 + b) Ū_T` and `m↓ = (1 - b) L̄_T`, where `b ∈ [0, 1)` is the null band, fixed in advance. For `t > T`:

* the increase statistic uses factors `f↑_t = 1 + λ↑ (X_t / m↑ - 1)` with bets λ↑ from a fixed grid in `[0, 1)`;
* if `L̄_T > 0` and `B = τ² / m↓ > 1`, the decrease statistic uses `f↓_t = 1 + λ↓ (1 - X_t / m↓)` with bets λ↓ from a grid of fractions of `1 / (B - 1)`, fixed at time T. Otherwise the decrease statistic is 0 for the whole segment. (A smaller lower bound means a larger B and therefore smaller admissible bets, so an uncertain baseline costs the decrease side power, and a lower bound of 0 switches it off.)
* Both statistics are the Shiryaev-Roberts recursions of 5.4, started at 0 at time T. A missing return leaves them unchanged.
* Two-sided alarm: the first `t > T` with `max(M↑_t, M↓_t) ≥ 2A`. One-sided use: one statistic at level A.
* Implementation details that keep the theorem applicable (`monitor.py`): when a window is degenerate, τ is set later, at a stopping time, which is still before the sequence starts; once the sequence has started it is never restarted, so each segment uses exactly one sequence; and every rule that ends calibration (fixed length, width target, cap, too few finite returns) uses only data seen so far, so T is a stopping time.

**Null hypothesis "no change".** In terms of the calibration averages `μ̄↑_T`, `μ̄↓_T` of 6.2:

* (H↑) `μ_t ≤ (1 + b) μ̄↑_T` for every `t > T`;
* (H↓) `μ_t ≥ (1 - b) μ̄↓_T` for every `t > T`.

Nothing is assumed about the calibration window beyond `0 ≤ X ≤ τ²`.

**Theorem 6.4.** Let `δ = δ_u + δ_l`.

1. Under (H↑) and (H↓), for every integer `m ≥ 1`, `P(an alarm at one of the times T + 1, ..., T + m) ≤ δ + m / A`. (More precisely `≤ P(Gᶜ) + P(G) m / A ≤ δ + (1 - δ) m / A` for m < A, with G defined in the proof.)
2. Under (H↑) and (H↓), the expected number of monitored steps before the first alarm is at least `(1 - δ) A`.
3. One side alone, at level A: under (H↑) the increase side alone has `P(alarm within m) ≤ δ_u + m / A`; under (H↓) the decrease side alone has `≤ δ_l + m / A`.

*Proof.* Let `G↑ = {ν̄↑_T < U_T}`, `G↓ = {ν̄↓_T > L_T}` and `G = G↑ ∩ G↓`. These events are known at time T.

*Step 1: the bounds hold at the random time T.* Apply Theorem 6.2 to the observations after `n_τ` (τ is fixed when they start). The theorem bounds the chance that the bound fails at *some* t, and `G↑ᶜ` is contained in that event whatever T is. So `P(G↑ᶜ) ≤ δ_u`, `P(G↓ᶜ) ≤ δ_l` and `P(Gᶜ) ≤ δ`.

*Step 2: on G both statistics are e-detectors.* Let `t > T`. On G, (H↑) and the clipping remark give `μ_t ≤ (1 + b) μ̄↑_T ≤ (1 + b) Ū_T = m↑`, so `E[f↑_t | F_{t-1}] = 1 + λ↑ (μ_t / m↑ - 1) ≤ 1`, and `f↑_t ≥ 1 - λ↑ > 0`. If the decrease side is on, (H↓) gives `μ_t ≥ (1 - b) μ̄↓_T > (1 - b) L̄_T = m↓`, so `E[f↓_t | F_{t-1}] ≤ 1`, and `f↓_t ≥ 1 - λ↓ (B - 1) > 0` because `X_t ≤ τ²`. The bets are fixed at time T, so they are predictable. As in Lemma 1 of 5.4, the average `M_t = (M↑_t + M↓_t) / 2` (each `M` itself an average over its grid) satisfies `E[M_t | F_{t-1}] ≤ M_{t-1} + 1` on `G ∩ {t > T}`. A missing return leaves M unchanged, so `M_t = M_{t-1} + 1{r_t observed} (f_t (M_{t-1} + 1) - M_{t-1})` factor by factor; under the missingness assumption of 6.7 the bound on `E[f_t | F_{t-1}]` also holds given that r_t is observed, which gives `E[M_t | F_{t-1}] ≤ M_{t-1} + P(r_t observed | F_{t-1}) ≤ M_{t-1} + 1`. If the decrease side is off, `M↓ = 0` and the same holds.

*Step 3: Ville's inequality after a random start.* Fix s and let `E_s = G ∩ {T = s}`, an event known at time s. For `s ≤ t ≤ s + m` put `S_t = 1_{E_s} (M_t + m - (t - s))`, with `M_s = 0`. Because `E_s` is known at time `t - 1`, Step 2 gives `E[S_t | F_{t-1}] ≤ S_{t-1}`. Also `S ≥ 0` and `S_s = m 1_{E_s}`. Ville's inequality (here Doob's maximal inequality for a nonnegative supermartingale on a finite horizon) gives `P(max_{s<t≤s+m} S_t ≥ A) ≤ E[S_s] / A = m P(E_s) / A`. An alarm needs `max(M↑, M↓) ≥ 2A`, which implies `M ≥ A` and hence `S ≥ A`. Summing over s: `P(G and an alarm within m steps) ≤ m P(G) / A`. (If `T = ∞` there is no monitoring and no alarm.)

*Step 4.* `P(alarm within m) ≤ P(Gᶜ) + m P(G) / A ≤ δ + m / A`. The function `p ↦ p + (1 - p) m / A` is increasing when `m < A`, which gives the sharper form.

*Part 2.* Assume `T < ∞` almost surely (true for every rule in the monitor as long as finite returns keep arriving; otherwise read the statement on `{T < ∞}`). Let `τ_A` be the first `t > T` with `M_t ≥ A`; the actual alarm comes no earlier. With `E_s` as above and `n ≥ 1`, `N = min(τ_A, s + n)` is a bounded stopping time, and `1_{E_s} (M_t - (t - s))`, `t ≥ s`, is a supermartingale by Step 2. Optional stopping gives `E[1_{E_s} (M_N - (N - s))] ≤ 0`, so `E[1_{E_s} (N - s)] ≥ E[1_{E_s} M_N] ≥ A P(E_s ∩ {τ_A ≤ s + n})`. Letting `n → ∞` (monotone convergence) gives `E[1_{E_s} (τ_A - s)] ≥ A P(E_s ∩ {τ_A < ∞})`; if `τ_A = ∞` with positive probability on `E_s`, the left side is infinite. Either way `E[1_{E_s} (τ_A - s)] ≥ A P(E_s)`. Summing over s: `E[τ_A - T] ≥ E[1_G (τ_A - T)] ≥ A P(G) ≥ (1 - δ) A`.

*Part 3.* The same proof with one statistic, level A, and only `G↑` (or `G↓`). ∎

**Corollary 6.5 (constant conditional law).** Suppose that, from step `n_τ + 1` until the change, the conditional law of `r_t` given `F_{t-1}` is one fixed law P₀ (for example i.i.d. returns, possibly with jumps and fat tails). Let `m(s) = E_{P₀}[min(r², s²)]`. Because τ is fixed at time `n_τ`, `μ_t = m(τ)` for every `t > n_τ` before the change, so `μ̄↑_T = μ̄↓_T = m(τ)` and (H↑), (H↓) hold with b = 0. Hence, with no change, **P(false alarm within m monitored steps) ≤ δ + m / A with the baseline estimated from the data.** More generally, if `r_t = σ_t ε_t` with i.i.d. innovations independent of the past and a predictable `σ_t` (5.1), then `μ_t = g_τ(σ_t)` with `g_τ` nondecreasing. If `σ_t = σ₀` during calibration, then `σ_t ≤ σ₀` afterwards implies (H↑) and `σ_t ≥ σ₀` afterwards implies (H↓). ∎

**Corollary 6.6 (recalibration after every alarm).** The monitor starts a new calibration after each alarm. Number the segments k = 1, 2, ...; segment k has its own τ, bounds, calibration averages and error budget δ = δ_u + δ_l. Suppose that (H↑) and (H↓) hold in every segment, relative to that segment's own calibration averages (no change anywhere). Let `N_m` be the number of alarms within the first m *monitored* steps (calibration steps are not counted). Then

`E[N_m] ≤ (δ + m / A) / (1 - δ)`, and therefore `P(N_m ≥ 1) ≤ (δ + m / A) / (1 - δ)`.

With δ = 0 (a known baseline in every segment) this reads `E[N_m] ≤ m / A`.

*Proof.* Let `G_k` be the event that both bounds of segment k hold at its end of calibration, and split `N_m = N_bad + N_good` according to whether the alarm's segment has `G_k`.
*Bad segments.* Segment k starts at a stopping time `S_k`. Its bounds use only returns after `S_k`, so Theorem 6.2 applied given the information at `S_k` gives `P(G_kᶜ | F_{S_k}) ≤ δ` on `{S_k < ∞}`. Let `C_k` be the event that fewer than m monitored steps have passed before `S_k`; it is known at `S_k`. A segment raises at most one alarm, and only segments with `C_k` can contribute to `N_m`, so `N_bad ≤ Σ_k 1{G_kᶜ} 1{C_k}` and `E[N_bad] ≤ Σ_k E[1{C_k} P(G_kᶜ | F_{S_k})] ≤ δ E[#{k : C_k}]`. Segment k + 1 starts right after alarm k, so `C_{k+1}` implies that alarm k is among the first m monitored steps. Hence `#{k : C_k} ≤ 1 + N_m` and `E[N_bad] ≤ δ (1 + E[N_m])`.
*Good segments.* Let `W_t` = (sum of the two-sided statistic `M = (M↑ + M↓) / 2` at all past alarms of good segments) + (the current M, if the current segment is good and monitoring) - (the number of monitored steps so far in good segments). Whether the current segment is good is known at its end of calibration, so by Step 2 of Theorem 6.4 every monitored step of a good segment has `E[ΔW | F_{t-1}] ≤ 0`. At an alarm the statistic moves into the sum unchanged, and in all other steps W does not move. So W is a supermartingale with `W_0 = 0`. Let σ be the time of the m-th monitored step. Each alarm of a good segment banks at least A (it needs `max(M↑, M↓) ≥ 2A`, so `M ≥ A`), and at most m monitored steps have passed by σ, so `W_{σ∧n} ≥ A N_good(σ∧n) - m` for every n. Optional stopping at the bounded time `σ ∧ n` gives `A E[N_good(σ∧n)] - m ≤ 0`, and monotone convergence as n → ∞ gives `E[N_good] ≤ m / A`.
Adding, `E[N_m] ≤ δ (1 + E[N_m]) + m / A`. Since `N_m ≤ m` is finite, this rearranges to the claim. ∎

**Choosing A and δ.** For a target α for one monitored "day" of m steps, choose `δ < α` and `A = m / (α - δ)`. E8 uses δ = 0.01 and α = 5%: `A = 390 / 0.04 = 9750` for a 390-step day (instead of 7800) and `A = 1440 / 0.04 = 36000` for a 1440-minute crypto day (instead of 28800). The larger A costs little, because the delay grows like log A (`log 9750 / log 7800 ≈ 1.025`). The main cost is the wider baseline (6.2, remark 3). In the two-sided monitor δ is split equally, `δ_u = δ_l = δ / 2`.

### 6.4 Subtle point (a): τ comes from the same calibration data

In the plug-in calibration, τ is c × MAD of the *whole* window. Then `X_i = min(r_i², τ²)` for a return i in that window depends, through τ, on returns that come after i. So `X_i` is not adapted, and the process K in Theorem 6.2 need not be a supermartingale. The plug-in mean is still a consistent estimate, but no finite-sample statement follows. Three valid fixes:

1. **Sample splitting (implemented).** τ from the first `n_τ` returns (default: 25% of the window), the confidence sequence on the rest. Needs nothing beyond 6.1. Cost: the sequence sees only 75% of the window, so the width grows by about `√(4/3) ≈ 1.15`. τ need not be precise: any τ known at `n_τ` keeps validity, it only changes power.
2. **A grid of truncation levels with a union bound.** Fix levels `τ_1 < ... < τ_K` before calibration (for example multiples of a rough scale known in advance). Run K sequences, each at level δ / K, on the whole window. At the end pick τ by any rule, for example the grid point nearest c × MAD of the whole window. All K bounds hold at once with probability at least 1 - δ, so the chosen one does too. Cost: `log(K/δ)` in place of `log(1/δ)` (K = 20, δ = 0.01: 7.6 against 4.6, so about 29% wider). Not implemented.
3. **A predictable running τ_t** (for example c × MAD of `r_1, ..., r_{t-1}`, as in E3, capped at a level fixed in advance so that one bound `B ≥ τ_t²` can normalise every `X_t`). Then `X_t = min(r_t², τ_t²)` is adapted, and Theorem 6.2 (applied to `X_t / B`) covers the average of `E[min(r_t², τ_t²) | F_{t-1}]`, a mean whose truncation level keeps changing, while monitoring uses one level. Under a constant conditional law (Corollary 6.5) `m(s)` is nondecreasing in s, so the comparison works in one direction only: the increase side must monitor with a level at most `min_i τ_i`, and the decrease side with a level at least `max_i τ_i`. This is valid under the constant-law assumption but not under (H↑), (H↓) alone. Not implemented.

### 6.5 Subtle point (b): what "no change" means when volatility clusters

* Theorem 6.2 needs neither stationarity nor independence. Under volatility clustering (GARCH, stochastic volatility) it still covers `μ̄_T`, the bet-weighted average of the conditional truncated variances *within the calibration window*. By contrast, `margin_z` uses the i.i.d. standard error, which is too small when squared returns are positively autocorrelated, so it has no guarantee there (E8, Q8f: under GARCH with persistence 0.99 the upper bound built from 2.33 i.i.d. standard errors missed on 1.6% of windows at a nominal 1%, while the upper confidence sequence at δ_u = 0.01 missed its target at some step on 0.05% of windows, and the lower one at δ_l = 0.01 on 0.175%; `results/e08f_garch_coverage.csv`).
* The monitoring null (H↑) then says: *the conditional truncated variance never rises above (1 + b) times its calibration-window average.* With clustering, `μ_t` moves all the time, so with b = 0, (H↑) fails whenever volatility rises above its calibration average, even if nothing structural happened. Alarms at those moments are not false alarms in the sense of Theorem 6.4. They are correct rejections of a null that is too strict for clustered data. The confidence sequence protects against *estimation error*; it cannot protect against this *model mismatch*.
* The remedy is the band b: "within a factor 1 ± b of the calibration average". Theorem 6.4 holds for any b fixed in advance, but whether (H↑) and (H↓) hold for a given b on real data is an empirical question (E4c; the GARCH experiment of E9). A longer calibration window that spans several clustering cycles also makes the calibration average closer to the long-run level.
* The target is weighted by the bets. With `d_t = n_plan`, `λ_i` is proportional to `1 / σ̂_{i-1}` until the cap binds, so the weights are close to uniform after the first few dozen observations. With the time-uniform schedule they decrease like `1 / √(i log i)`, which puts more weight on early calibration data. The weighting matters only when `μ_i` drifts within the window.

### 6.6 Subtle point (c): a calibration length chosen from the data

Theorem 6.2 bounds the chance that the bound fails at *any* t. The event "the bound holds at every t" is contained in "the bound holds at T" for every random T. So Theorem 6.4 holds for every stopping rule. The monitor offers one (`target_rel_width = w`): stop at the first `t ≥ cal_len` with relative width `(Ū_t - L̄_t) / X̄_t ≤ w` (`X̄_t` is the plain mean of the sequence part), or at `cal_max`.

This would not be justified with a fixed-n interval `I_n`, one with `P(μ ∈ I_n) ≥ 1 - δ` for each fixed n. Such intervals say nothing about `I_T` at a data-dependent T. For example, the normal-approximation interval checked at every n excludes the true mean at some n with probability one (law of the iterated logarithm), so a rule that waits for a favourable look will get one. The width rule is milder (it waits for a small width, which depends on the data mostly through the variance estimate), but a fixed-n interval still gives no proof for it. The time-uniform sequence pays for this freedom with a wider bound at any single n than a fixed-n interval with the same δ (6.2, remark 3; E8, Q8d).

The planned size `n_plan` may be any value fixed before the sequence starts, even one computed from the first `n_τ` returns, because the bets need only be predictable.

### 6.7 Status

* **Proved here:** Lemma 6.1; Theorem 6.2 and Corollary 6.3 (given Ville's inequality); Theorem 6.4, parts 1 to 3 (given Theorem 6.2, Ville's inequality and optional stopping at bounded stopping times); Corollaries 6.5 and 6.6.
* **Cited:** Ville's inequality (Ville 1939; a modern statement is in Howard et al. 2020); the predictable plug-in empirical Bernstein construction and its bet schedule (Waudby-Smith and Ramdas 2024, Theorem 2); the inequality of Lemma 6.1 for a = 1 (Fan, Grama and Liu 2015), which is also the basis of the empirical Bernstein confidence sequences in Howard et al. (2021); e-detectors (Shin, Ramdas and Rinaldo 2024).
* **Checked numerically:** coverage on bounded i.i.d. data and with time-varying means (`tests/test_confseq.py`); the alarm count bound of Corollary 6.6 on long i.i.d. streams (`tests/test_monitor_cs.py`); false alarm rates, bound failures, delays and adaptive lengths (E8). In E8, with δ = 0.01 and A = 9750, the one-day false alarm rate of the cs baseline was at most 0.27% in nine synthetic settings against a proved 5%, and the confidence bound failed on at most 0.08% of paths against an allowed 1% (`results/e08a_null.csv`).
* **Assumed, not proved:** (H↑) and (H↓) on real data (6.5); the deseasonalising profile is treated as known (5.6, item 3); missing returns are skipped, which is valid when, given the past, whether a return is missing says nothing about its size (the conditional means in (H↑), (H↓) are then the same with or without conditioning on "observed"). No bound on the detection delay is proved; delays are measured in E8.

**References for this section.** Fan, X., Grama, I. and Liu, Q. (2015). Exponential inequalities for martingales with applications. *Electronic Journal of Probability* 20, paper 1, 1 to 22. doi:10.1214/EJP.v20-3496. Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2020). Time-uniform Chernoff bounds via nonnegative supermartingales. *Probability Surveys* 17, 257 to 317. doi:10.1214/18-PS321. Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2021). Time-uniform, nonparametric, nonasymptotic confidence sequences. *Annals of Statistics* 49(2), 1055 to 1080. doi:10.1214/20-AOS1991. Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science* 2(2), 229 to 260. doi:10.51387/23-NEJSDS51. Ville, J. (1939). *Étude critique de la notion de collectif.* Thèse, Paris (numdam.org/item/THESE_1939__218__1_0). Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *JRSS B* 86(1), 1 to 27. doi:10.1093/jrsssb/qkad009.
