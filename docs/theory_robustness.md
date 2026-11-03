# The false alarm guarantee when the null holds only approximately (TH4)

Tracker B84. `docs/math.md` (Theorem B) proves that the default detector raises a false alarm within m steps with probability at most m / A, provided the conditional mean of the clipped squared return never leaves the baseline. E16 (`docs/report_parts/E16.md`) mapped what happens when that condition fails slightly: GARCH clustering, a plug-in baseline, price ticks. This document turns that map into proved statements: how much a small drift of the conditional mean can inflate the false alarm probability, when a drift that averages out is harmless, what cannot be said, and what the plug-in baseline costs. Code: `src/avmonitor/robustness.py`. Experiment: `experiments/th4_robustness.py` (outputs `results/th4_*.csv`, figures `results/figures/th4_*.png`, all **[SYNTHETIC]**). Tests: `tests/test_robustness.py`. Report part: `docs/report_parts/TH4.md`.

Each claim carries one of three labels: **[PROVED]** (a complete proof is given here, using only Ville's inequality and optional stopping at bounded stopping times, both already cited in `docs/math.md`), **[NUMERICAL]** (simulation or numerical evaluation) and **[HEURISTIC]** (a formula used outside the conditions under which it is proved, checked numerically).

## 0. Summary in plain words

The monitor's guarantee (Theorem B) assumes that, given everything seen so far, the expected clipped squared return is exactly the baseline m₀. Real data never satisfy this exactly. This document asks how much larger the false alarm probability can get when the expected value is off by a small relative amount, called the **drift** δ: the expected value is at most `m₀ (1 + δ)` (an upward drift, which pushes the increase side towards a false alarm) or at least `m₀ (1 - δ)` (a downward drift, which pushes the decrease side).

1. **A drift that never exceeds a known δ (Theorem 1, proved).** The probability of a false alarm within m steps is at most `m / A` times an inflation factor `I_m(δ)`, which has a closed form. The factor stays close to 1 as long as the drift times the horizon is small: for one day of 390 steps (A = 7800, bound 5%), a drift of 0.1% on both sides gives 5.14%, 0.5% gives 5.80% and 1% gives 7.06%. Beyond that it grows exponentially. The one-day bound reaches 10% at a drift of 1.63% and 25% at 2.66% (upward drift only: 1.72% and 2.70%; downward only: 6.36% and 10.51%, because the decrease bets are smaller). Over five days a drift of 0.5% already makes the bound useless. The bound cannot be lowered in general: a stream that sits exactly at `m₀ (1 + δ)` attains it at the step where it alarms. In one day an adaptive stream that watches the detector reaches a third to two thirds of it.
2. **A drift that varies (Theorems 2 and 3, Propositions 2b, 4 and 5, proved).** Under GARCH clustering the drift has no useful upper limit, so Theorem 1 says nothing. When the drift goes up and down, a calm spell undoes a burst. What matters is then the largest product of `(1 + bet × drift)` over a window of time, not the largest drift. This gives a bound for any drift that is known one step ahead, including a drift that responds to the data (Theorem 2). A sharper bound holds when the drift does not respond to the data (Theorem 3). Two things cannot be done. The expected product of these factors does not bound the expected wealth when the drift responds to the data (a two-step counterexample). And the average drift alone cannot give a bound of the form "constant × m / A" (Proposition 5: a short burst decided before monitoring starts). A bound in terms of the average drift does exist, but it carries an extra term that does not shrink as A grows (Proposition 4).
3. **Plug-in baseline (Corollary 1, proved).** A baseline estimated from a calibration window is wrong by a random factor that then stays fixed during monitoring. So Theorem 1 applies to each window, and the false alarm probability is at most the average of the capped bound over calibration windows. For Gaussian returns with 390 calibration steps this is 42% for one day and 92% for five days. E16 measured 4.0% and 41.8% for i.i.d. returns, and 32.5% and 85.8% for the worst stream it found with the same conditional mean, so the bound is within a factor of 1.3 of the worst case. Bringing the one-day bound down to 6% needs about 23,000 calibration steps. This agrees with Proposition A.6 of `docs/theory_running_tau.md`: the constant in front of `m / A` has to grow with the horizon.
4. **Numerical check against E16 (numerical).** Every proved bound lies above every measured rate where its assumption holds. That covers 128 GARCH cells, 36 plug-in cells, 80 decile cells, 32 switching cells and 80 constant-drift cells, with no exception at the 95% level. For GARCH the proved bounds are informative only where the drift is small: below 1 in 22 of the 64 cells at one day (6 to 29 times the measured rate) and in 6 at five days. The constant-drift bound gives 1 in every GARCH cell. Applying the drift-path bound of Theorem 3 to GARCH is not proved (point 2), but it is closer (1.2 to 8 times the measured rate) and was never below it.

In short, a drift of a few tenths of a percent costs little over a day. A persistent drift of 2% or more, or any drift sustained over several days, uses up the guarantee. The plug-in baseline's error is exactly such a persistent drift. The remedies are the ones the repository already has: a confidence-sequence baseline (E8), a forecast-relative null (E10) and the rank detector (E14).

## 1. Set-up and notation

As in `docs/math.md`, sections 5.1 to 5.4.

* `X_t = min(r_t², τ²)`, baseline `m₀ > 0`, `Y_t = X_t / m₀ ∈ [0, B]` with `B = τ² / m₀ > 1`. `F_t` is the information up to time t; `F_0` contains the calibration data. τ, m₀ and all bets are known at time 0 (`F_0`-measurable).
* **Increase bets** `λ_1, ..., λ_K ∈ [0, 1)` with weights `w_k > 0`, `Σ w_k = 1`; factors `f_{k,t} = 1 + λ_k (Y_t - 1) ≥ 1 - λ_k > 0`. **Decrease bets** `φ_1, ..., φ_L ∈ [0, 1 / (B - 1))` with weights `v_l`; factors `g_{l,t} = 1 + φ_l (1 - Y_t) ≥ 1 - φ_l (B - 1) > 0`. The default (`detectors.DEFAULT_UP_GRID`, `DEFAULT_DOWN_FRACTIONS`): `λ = 0.02, 0.05, 0.1, 0.2, 0.35, 0.5` and `φ_l = s_l / (B - 1)` with `s = 0.1, 0.25, 0.5, 0.75, 0.95`, equal weights. So `λ_max = 0.5` and `φ_max = 0.95 / (B - 1)`, which is 0.118 at the Gaussian value `B = 9.045` (c = 3). Mean bets: `λ̄ = Σ w_k λ_k = 0.2033` and `φ̄ = Σ v_l φ_l = 0.51 / (B - 1)` (0.0634 at B = 9.045).
* **Shiryaev-Roberts statistic of one bet** (called M in `docs/math.md`): `R_t = f_t (R_{t-1} + 1)`, `R_0 = 0`, so `R_t = Σ_{j=1}^{t} Π_{i=j}^{t} f_i`. Side averages `M↑_t = Σ_k w_k R↑_{k,t}` and `M↓_t = Σ_l v_l R↓_{l,t}`, and `N_t = (M↑_t + M↓_t) / 2`.
* **Alarm.** Two-sided default rule: `T = inf{t ≥ 1 : max(M↑_t, M↓_t) ≥ 2A}`. One side alone at level A: `T↑ = inf{t : M↑_t ≥ A}`, and `T↓` likewise. Since `M↑, M↓ ≥ 0`, an alarm of the two-sided rule at time T implies `N_T ≥ A`.
* **Signed drift.** `d_t = E[Y_t | F_{t-1}] - 1 ∈ [-1, B - 1]`. It is predictable (known at time t - 1). `E[f_{k,t} | F_{t-1}] = 1 + λ_k d_t` and `E[g_{l,t} | F_{t-1}] = 1 - φ_l d_t`. Theorem B is the case `d_t ≤ 0` (increase side) or `d_t ≥ 0` (decrease side) at every step.

**The approximate null.** Predictable `δ↑_t, δ↓_t ≥ 0` with

* (A↑) `E[X_t | F_{t-1}] ≤ m₀ (1 + δ↑_t)`, equivalently `d_t ≤ δ↑_t`, so `E[f_{k,t} | F_{t-1}] ≤ 1 + λ_k δ↑_t`;
* (A↓) `E[X_t | F_{t-1}] ≥ m₀ (1 - δ↓_t)`, equivalently `-d_t ≤ δ↓_t`, so `E[g_{l,t} | F_{t-1}] ≤ 1 + φ_l δ↓_t`.

The smallest choices are the positive and negative parts `δ↑_t = d_t⁺ = max(d_t, 0)` and `δ↓_t = d_t⁻ = max(-d_t, 0)`. At each step at most one of them is positive, but the bounds below allow both sides to drift over a horizon, as under GARCH, where calm spells push the conditional mean below m₀ and bursts push it above. E16 summarises the drift of a GARCH cell by the mean absolute drift `E|d_t|` (`drift_mean_abs` in `results/e16b_garch.csv`, averaged over time and paths, with the 1% and 99% quantiles of `d_t` in `drift_q01`, `drift_q99`).

**Two growth functions.** For `x > -1` and an integer m ≥ 1,

`S_m(x) = Σ_{j=1}^{m} (1 + x)^j = (1 + x) ((1 + x)^m - 1) / x` (with `S_m(0) = m`), and `G_m(x) = S_m(x) / m`.

## 2. The sup form

### 2.1 Expected wealth under a drift

**Lemma 1.** [PROVED] Let `f_t ≥ 0` be adapted, and let `ρ_1, ρ_2, ...` be numbers at least 1 with `E[f_t | F_{t-1}] ≤ ρ_t` for every t. Let `R_t = f_t (R_{t-1} + 1)`, `R_0 = 0`, and assume every `f_t` is bounded. Then for every stopping time T (possibly infinite) and every t ≥ 0,

`E[R_{T∧t}] ≤ s_t := Σ_{j=1}^{t} Π_{i=j}^{t} ρ_i`.

In particular (T = ∞) `E[R_t] ≤ s_t`, and with a constant `ρ_i = 1 + λδ`: `E[R_t] ≤ Σ_{j=1}^{t} (1 + λδ)^{t-j+1} = S_t(λδ)`.

*Proof.* Put `a_t = E[R_{T∧t}]` (finite, since the factors are bounded). Split `R_{T∧t} = R_T 1{T ≤ t-1} + R_t 1{T > t-1}`. The event `{T > t-1}` is known at time t - 1, so

`E[R_t 1{T > t-1}] = E[1{T > t-1} E[f_t | F_{t-1}] (R_{t-1} + 1)] ≤ ρ_t E[1{T > t-1} (R_{t-1} + 1)]`.

Because `ρ_t ≥ 1` and `R ≥ 0`, also `E[R_T 1{T ≤ t-1}] ≤ ρ_t E[R_T 1{T ≤ t-1}]`. Adding, and using `R_{T∧(t-1)} = R_T 1{T ≤ t-1} + R_{t-1} 1{T > t-1}`, gives `a_t ≤ ρ_t (a_{t-1} + P(T > t-1)) ≤ ρ_t (a_{t-1} + 1)`. Since `a_0 = 0 = s_0` and `s_t = ρ_t (s_{t-1} + 1)`, induction gives `a_t ≤ s_t`. ∎

The condition `ρ_t ≥ 1` is used only for the part of the mass that has already stopped. Section 3 shows what changes when the drift can be negative.

### 2.2 Theorem 1

**Theorem 1 (sup form).** [PROVED] Suppose (A↑) and (A↓) hold with `δ↑_t ≤ δ↑` and `δ↓_t ≤ δ↓` for every t, where `δ↑, δ↓ ≥ 0` are constants (or known at time 0). Then for every m ≥ 1:

1. Two-sided default rule: `P(T ≤ m) ≤ (1 / (2A)) [Σ_k w_k S_m(λ_k δ↑) + Σ_l v_l S_m(φ_l δ↓)] = (m / A) I_m(δ↑, δ↓)`, with the **inflation factor**

   `I_m(δ↑, δ↓) = ½ [Σ_k w_k G_m(λ_k δ↑) + Σ_l v_l G_m(φ_l δ↓)]`, `G_m(x) = (1 + x) ((1 + x)^m - 1) / (m x)`.

2. Increase side alone at level A: `P(T↑ ≤ m) ≤ (m / A) Σ_k w_k G_m(λ_k δ↑)`. Decrease side alone: `P(T↓ ≤ m) ≤ (m / A) Σ_l v_l G_m(φ_l δ↓)`.
3. With `δ↑ = δ↓ = 0` both reduce to Theorem B, since `G_m(0) = 1`. If δ↑ and δ↓ are known at time 0 but random, the statements hold conditionally on `F_0`: `P(T ≤ m | F_0) ≤ min{1, (m / A) I_m(δ↑, δ↓)}`.

*Proof.* On `{T ≤ m}`, `N_{T∧m} = N_T ≥ A`, and N ≥ 0, so `A P(T ≤ m) ≤ E[N_{T∧m}]`. By Lemma 1 with `ρ = 1 + λ_k δ↑`, `E[R↑_{k,T∧m}] ≤ S_m(λ_k δ↑)`, and with `ρ = 1 + φ_l δ↓`, `E[R↓_{l,T∧m}] ≤ S_m(φ_l δ↓)`. Average with the weights. Part 2 is the same argument with one side and level A. For part 3 apply Lemma 1 with conditional expectations given `F_0` (the δ's are then fixed numbers and `{T > t-1}` is still known at time t - 1). ∎

Two remarks. First, the bound is per bet: each bet is inflated by its own `G_m(λ_k δ)`. A bound that inflates the whole average by its largest bet (a single factor `1 + λ_max δ` in a recursion for `M↑`) is true but weaker. Second, the decrease side uses `φ_l = s_l / (B - 1)`, so the same drift costs it about `λ_max (B - 1) / 0.95 ≈ 4.2` times less at the Gaussian B.

### 2.3 When the factor is close to 1, and how it blows up

**Proposition 1.** [PROVED] For x ≥ 0 and m ≥ 1,

`1 + x (m + 1) / 2 ≤ (1 + x)^{(m+1)/2} ≤ G_m(x) ≤ exp(x (m + 1) / 2 + x² m² / 24)`.

Hence, with `x* = max(λ_max δ↑, φ_max δ↓)`,

`1 + (m + 1) (λ̄ δ↑ + φ̄ δ↓) / 4 ≤ I_m(δ↑, δ↓) ≤ exp(x* (m + 1) / 2 + x*² m² / 24)`,

and `I_m(δ↑, δ↓) ≥ (w_max / 2) ((1 + λ_max δ↑)^{m+1} - (1 + λ_max δ↑)) / (λ_max δ↑ m)`, where `w_max` is the weight of the largest increase bet (1/6 for the default).

*Proof.* Lower bound: `G_m(x)` is the average of `(1 + x)^j` over j = 1, ..., m, which is at least their geometric mean `(1 + x)^{(m+1)/2}` (arithmetic and geometric means); Bernoulli's inequality `(1 + x)^r ≥ 1 + r x` for r ≥ 1 gives the first inequality. Upper bound: `(1 + x)^j ≤ e^{x j}`, and for x > 0

`(1/m) Σ_{j=1}^{m} e^{x j} = e^{x (m+1)/2} sinh(m x / 2) / (m sinh(x / 2)) ≤ e^{x (m+1)/2} sinh(z) / z` with `z = m x / 2`,

using the geometric sum and `sinh(x/2) ≥ x/2`. Finally `sinh(z) / z = Σ_n z^{2n} / (2n+1)!` and `exp(z² / 6) = Σ_n z^{2n} / (6^n n!)`, and `(2n+1)! ≥ 6^n n!` for every n ≥ 0 (equality for n = 0, 1; going from n to n + 1 multiplies the left side by `(2n+2)(2n+3)` and the right side by `6 (n+1)`, and `2 (2n+3) ≥ 6`), so `sinh(z) / z ≤ exp(z² / 6) = exp(m² x² / 24)`. For I: the lower bound follows by averaging `1 + x(m+1)/2` over the bets (it is linear in x); the upper bound because `G_m` is increasing in x and each `λ_k δ↑`, `φ_l δ↓` is at most x*. The last claim keeps only the largest increase bet. ∎

**Reading.** The factor is close to 1 exactly when `x* m` is small: to first order `I_m ≈ 1 + (m + 1)(λ̄ δ↑ + φ̄ δ↓) / 4`, and `I_m ≤ e^{0.051} ≈ 1.05` whenever `x* (m + 1) ≤ 0.1`. When `λ_max δ m` is large, `I_m` grows exponentially in m, like `(1 + λ_max δ)^m / (λ_max δ m)`: the bound reaches 1 after a number of steps of order `log(A) / (λ_max δ)`. This is the same time scale as the drift of the plug-in in Proposition A.6 of `docs/theory_running_tau.md` (`log A / η` there).

### 2.4 Sharpness

**Proposition 2.** [PROVED]

1. *One side, one bet, level A.* The stream `Y_t = 1 + δ` at every step (allowed if `1 + δ ≤ B`) satisfies (A↑) with equality, and `R_t = S_t(λδ)` exactly, so `T↑ = m* := min{t : S_t(λδ) ≥ A}`. Hence `P(T↑ ≤ m) = 1` for m ≥ m* and 0 before. Since `S_{m*}(λδ) = (1 + λδ)(S_{m*-1}(λδ) + 1) < (1 + λδ)(A + 1)`, no bound of the form `c S_m(λδ) / A` valid for this stream can have `c < 1 / ((1 + λδ)(1 + 1/A))`: the inflation factor of Theorem 1 cannot be lowered at the horizon m*.
2. *Two-sided default.* For the same stream, at its alarm step m*, Theorem 1 (with `δ↓ = 0`) gives at most `(1 + λ_max δ)(1 + 1/(2A)) + m* / (2A)`, while the probability is 1. With A = 7800 and δ = 0.03 the stream alarms at step 472 and the bound there is 1.038.

*Proof.* Part 1 is a direct computation. Part 2: the increase average is `M↑_t = Σ_k w_k S_t(λ_k δ)`, which is below 2A at m* - 1 and at most `(1 + λ_max δ)(M↑_{m*-1} + 1)` at m*; the decrease side contributes `Σ_l v_l S_{m*}(0) = m*`. The numbers come from `robustness.deterministic_alarm_step` and `sup_bound`. ∎

So Theorem 1 is attained, up to these small factors, at the horizon where the deterministic stream alarms. At shorter horizons the worst case is not known; section 5.2 compares the bound with an adaptive adversary that has the same drift.

## 3. Pathwise forms: drift that varies, and what cannot be said

Theorem 1 needs a bound on the drift at every step. Under GARCH the conditional mean comes arbitrarily close to `B m₀` when σ_t is large (almost every return is then clipped), so the sup form only gives `δ↑ = B - 1` and a useless bound, although the drift averages to 0 over time. This section replaces `(1 + λδ)^k` by the product of `(1 + λ δ_i)` over the window.

### 3.1 A drift-corrected detector and a pathwise inequality

**Lemma 2.** [PROVED] Let `u_t` and `ℓ_t` be predictable with `u_t ≥ d_t` and `ℓ_t ≥ -d_t` (for example the exact drift, `u_t = d_t` and `ℓ_t = -d_t`; signs are allowed), and with `u_t ≥ -1` and `ℓ_t ≥ 1 - B`. Put `ρ↑_{k,t} = 1 + λ_k u_t > 0` and `ρ↓_{l,t} = 1 + φ_l ℓ_t > 0`. Define the **drift-corrected** factors `f_{k,t} / ρ↑_{k,t}` and `g_{l,t} / ρ↓_{l,t}`, their Shiryaev-Roberts statistics and the two-sided average `Ñ_t` built in the same way as `N_t`. Then:

1. `E[Ñ_t | F_{t-1}] ≤ Ñ_{t-1} + 1`, so for every a > 0 and m ≥ 1, `P(max_{t ≤ m} Ñ_t ≥ a) ≤ m / a`.
2. Path by path, `N_t ≤ D*_t Ñ_t` for every t, where `D*_t` is the largest of the window products `Π_{i=j}^{t} ρ↑_{k,i}` and `Π_{i=j}^{t} ρ↓_{l,i}` over all bets k, l and all starts `1 ≤ j ≤ t`. `D*_t` is predictable and obeys the CUSUM recursion `C_t = ρ_t max(C_{t-1}, 1)`, `C_0 = 0`, for each bet, `D*_t = max over bets of C_t`.

*Proof.* (1) `ρ` is positive and predictable, so `E[f_{k,t} / ρ↑_{k,t} | F_{t-1}] = (1 + λ_k d_t) / (1 + λ_k u_t) ≤ 1` (as `λ_k ≥ 0` and `d_t ≤ u_t`), and likewise for the decrease factors. The rest is Lemma 1 of `docs/math.md` (5.4), the remark there on averages of several bets, and the proof of Theorem B with the level a in place of A (Ville's inequality applied to `Ñ_t + m - t`). (2) The wealth of bet k started at j is `Π_{i=j}^{t} f_{k,i} = (Π_{i=j}^{t} f_{k,i} / ρ↑_{k,i}) × Π_{i=j}^{t} ρ↑_{k,i}`, and the second factor is at most `D*_t`. Sum over j and average over bets and sides. The CUSUM recursion is the usual one for `max_{j ≤ t} Π_{i=j}^{t} ρ_i`. ∎

The drift-corrected detector is itself a valid detector whenever the drift is known one step ahead; it is the multiplicative cousin of E16's control that bets against the true conditional mean, and of E10's forecast-relative monitor. The experiment reports its false alarm rate as a check (it must stay below m / A).

### 3.2 Theorem 2: any predictable drift

**Theorem 2 (threshold form).** [PROVED] Under the conditions of Lemma 2, for every m ≥ 1 and every K > 0,

`P(T ≤ m) ≤ K m / A + P(max_{t ≤ m} D*_t > K)`.

*Proof.* On `{T ≤ m}`, `A ≤ N_T ≤ D*_T Ñ_T`. Either `D*_T > K`, which lies in `{max_{t ≤ m} D*_t > K}`, or `Ñ_T ≥ A / K`, which has probability at most `K m / A` by Lemma 2 (1) with a = A / K. ∎

Because the drift is signed, a calm spell (`d < 0`) shrinks the increase products, and the CUSUM restarts: a drift that averages out over a window does not accumulate, and `D*` is governed by the largest excursion of the drift, not by its total. No independence is needed: the drift may respond to the data, as under GARCH. Under a constant drift Theorem 2 is weaker than Theorem 1 (with `K = (1 + x*)^m` it loses a factor of about `x* m` when that is large), so it is a tool for varying drift. The bound needs the law of `max_{t ≤ m} D*_t`, which depends on the drift path only, not on the detector.

**Proposition 2b (truncated sup form).** [PROVED] For any predictable drift and any `δ↑, δ↓ ≥ 0`,

`P(T ≤ m) ≤ (m / A) I_m(δ↑, δ↓) + P(max_{t ≤ m} d_t⁺ > δ↑ or max_{t ≤ m} d_t⁻ > δ↓)`.

*Proof.* Let `σ = inf{t : d_t⁺ > δ↑ or d_t⁻ > δ↓}`; since `d_t` is known at time t - 1, the event `{t < σ}` is known at time t - 1. Replace every factor by 1 from step σ on. The modified factors have conditional means at most `1 + λ_k δ↑` and `1 + φ_l δ↓` at every step, so by Theorem 1 the modified detector alarms within m with probability at most `(m / A) I_m`. If `T ≤ m` and `T < σ`, the modified and the original statistics agree up to T, so the modified detector alarms at T; otherwise `σ ≤ T ≤ m`. ∎

### 3.3 Theorem 3: exogenous drift

**Theorem 3 (pathwise form).** [PROVED] Suppose there is a σ-field G (the "drift driver") such that the drift bounds `δ↑_t, δ↓_t ≥ 0` are G-measurable and, for every t, `E[Y_t | F_{t-1} ∨ G] ≤ 1 + δ↑_t` and `E[Y_t | F_{t-1} ∨ G] ≥ 1 - δ↓_t`. Then

`P(T ≤ m | G) ≤ min{1, (1 / (2A)) [Σ_k w_k Σ_{j=1}^{m} Π_{i=j}^{m} (1 + λ_k δ↑_i) + Σ_l v_l Σ_{j=1}^{m} Π_{i=j}^{m} (1 + φ_l δ↓_i)], D*_m m / A}`,

where in the last term `D*_m` is the largest window product of Lemma 2 for the signed drift given G. Hence `P(T ≤ m)` is at most the expectation of the right side over the drift path.

*Proof.* Work with the larger filtration `H_t = F_t ∨ G`. The factors are adapted to it, T is a stopping time for it, and given `H_0` the numbers `ρ_t = 1 + λ_k δ↑_t ≥ 1` are fixed. Lemma 1 with conditional expectations given `H_0` gives `E[R_{T∧m} | H_0] ≤ Σ_j Π_{i=j}^{m} ρ_i`, and the proof of Theorem 1 gives the first term. For the second, Lemma 2 holds for H (the conditional means given `H_{t-1}` are bounded by the G-measurable drift), and Ville's inequality holds conditionally on `H_0` (apply it to `1_E (Ñ_t + m - t)` for each event E in `H_0`, as in Step 3 of the proof of Theorem 6.4 in `docs/math.md`), so Theorem 2 holds given `H_0` with the `H_0`-measurable choice `K = D*_m`. ∎

The sum `Σ_j Π_{i=j}^{m} (1 + λ δ_i)` obeys the Shiryaev-Roberts recursion `s_t = (1 + λ δ_t)(s_{t-1} + 1)`: the bound is the detector run on the drift itself (`robustness.drift_sr_sums`). Examples of an exogenous drift:

* a deterministic drift path, for instance a known error in the intraday profile used for deseasonalising (`docs/math.md` 5.6, item 3), with G trivial;
* a drift decided at time 0, for instance by the calibration window (the plug-in baseline, section 4), with `G = F_0`;
* stochastic volatility independent of the innovations: `r_t = σ_t ε_t` with i.i.d. `ε_t` independent of the whole volatility path and of `F_0`. Given G = the volatility path (and `F_0`), the returns are independent with `E[X_t | F_{t-1} ∨ G] = E[min(σ_t² ε², τ²)]`, a fixed function of σ_t. Section 5.3 uses a two-state version.

### 3.4 What cannot be said

**(a) The expected product does not bound the expected wealth when the drift responds to the data.** [PROVED] For an exogenous drift, Theorem 3 gives `E[R_t] ≤ Σ_j E[Π_{i=j}^{t} (1 + λ δ_i)]`. Without exogeneity this fails. Two steps: `Y_1` is 0 or 2 with probability ½ each (so `d_1 = 0`), and `d_2 = c 1{Y_1 = 2}` with `Y_2 = 1 + d_2` (any `c ≤ B - 1`). Then

`E[f_1 f_2] = ½ (1 - λ) + ½ (1 + λ)(1 + λ c) = 1 + λ c (1 + λ) / 2 > 1 + λ c / 2 = E[(1 + λ d_1)(1 + λ d_2)]`.

With λ = 0.5 and c = 2: 1.75 against 1.5, and `E[R_2] = 1.75 + 1.5 = 3.25` against `Σ_j E[Π (1 + λ δ_i)] = 3`. The drift is high exactly on the paths where the wealth is already high. This is the GARCH mechanism (a large return raises the next conditional variance), so for GARCH only Theorem 2 and Propositions 2b and 4 apply; Theorem 3 there is a heuristic.

**(b) The mean drift gives a lower estimate, not a bound (Jensen).** [PROVED] Even for an exogenous drift, the bound of Theorem 3 depends on the law of the window products, not on the mean drift. Since `log(1 + x) ≤ x`, `Π_{i=j}^{m} (1 + λ δ_i) ≤ exp(λ Σ_{i=j}^{m} δ_i)`, so the bound is controlled by the moment generating function of the cumulative drift over windows. By Jensen's inequality `E[exp(λ S)] ≥ exp(λ E[S])`: knowing E[S] bounds this quantity from below only. Two cases are exact. If the `δ_i` are independent across steps (a drift that averages out at every step), `E[Π_{i=j}^{m} (1 + λ δ_i)] = Π_{i=j}^{m} (1 + λ E[δ_i])`, so the expected Theorem 3 bound equals the sup form evaluated at the mean positive drift, and with the cap at 1 inside the expectation it can only be smaller. If instead the drift is persistent, for example `δ_i = ε` for the whole window with probability ½ and 0 otherwise (mean ε/2 at every step), the expected bound is about half of the sup form at ε, which is much larger than the sup form at ε/2 when `λ ε m` is large. Section 5.3 shows both regimes.

**(c) A bound in terms of the mean drift exists, but only with an additive term that does not shrink with A.**

**Proposition 4 (mean drift).** [PROVED] For any predictable drift,

`P(T ≤ m) ≤ m / A + (1 + 1/(2A)) E[Σ_{t ≤ T∧m} (λ_max δ↑_t + φ_max δ↓_t)] ≤ m / A + (1 + 1/(2A)) m (λ_max δ̄↑ + φ_max δ̄↓)`,

where `δ̄↑ = (1/m) Σ_{t ≤ m} E[δ↑_t]` is the mean drift over time and paths (and `δ̄↓` likewise). For one side alone at level A, the factor is `(1 + 1/A)` and only that side's term appears.

*Proof.* On `{T > t-1}` both side averages at t - 1 are below 2A. Using `E[f_{k,t} | F_{t-1}] ≤ 1 + λ_k δ↑_t ≤ 1 + λ_max δ↑_t` and the same for g,

`E[N_t | F_{t-1}] ≤ N_{t-1} + 1 + ½ [λ_max δ↑_t (M↑_{t-1} + 1) + φ_max δ↓_t (M↓_{t-1} + 1)] ≤ N_{t-1} + 1 + c_t`, with `c_t = (A + ½)(λ_max δ↑_t + φ_max δ↓_t)`.

So `V_t = N_{T∧t} - Σ_{s ≤ T∧t} (1 + c_s)` is a supermartingale with `V_0 = 0` (its increment is `1{T > t-1}` times a term with conditional mean at most 0), and `E[V_m] ≤ 0` gives `E[N_{T∧m}] ≤ E[T∧m] + E[Σ_{s ≤ T∧m} c_s]`. Finally `A P(T ≤ m) ≤ E[N_{T∧m}]` and `E[T∧m] ≤ m`. ∎

The drift enters weighted by A, because the proof must allow the drift to act when the statistic is already near the threshold. At one day (m = 390) the second term is `(1 + 1/(2A)) 195 δ̄↑` from the increase side alone: it is larger than the whole budget m / A = 5% as soon as the mean upward drift exceeds 0.026%, so the bound at least doubles. The next proposition shows that the additive form is not an artefact of the proof.

**Proposition 5 (no multiplicative bound in terms of the mean drift).** [PROVED] Let `s* = s*(A, B)` be the step at which the stream `Y_t = B` makes the two-sided default detector alarm (`s* = 7` for A = 7800 and B = 9.045). Let ξ be a coin with `P(ξ = 1) = p`, independent of everything else and known at time 0. If ξ = 1, put `Y_t = B` for `t ≤ s*` and `Y_t = 1` afterwards; if ξ = 0, `Y_t = 1` throughout. Then for every m with `s* ≤ m < 2A`:

`P(T ≤ m) = p`, while the mean absolute drift over the m steps is `δ̄ = p (B - 1) s* / m`, so `P(T ≤ m) = δ̄ m / ((B - 1) s*)`.

Consequently no inequality `P(T ≤ m) ≤ C(δ̄) m / A` with C depending only on the mean drift can hold for all A: take `m = s*(A)` and `p = δ̄ / (B - 1)`; then `P(T ≤ m) = δ̄ / (B - 1)` while `C(δ̄) s*(A) / A → 0`, because `s*(A) ≤ ⌈log(2A K) / log(1 + λ_max (B - 1))⌉` grows only like log A (K the number of increase bets with equal weights).

*Proof.* If ξ = 0, every factor equals 1, so `M↑_t = M↓_t = t < 2A` for t < 2A: no alarm. If ξ = 1, the stream is deterministic and alarms at s* by definition; the drift is `d_t = B - 1` for t ≤ s* and 0 after. For the bound on s*: `M↑_t ≥ (1/K)(1 + λ_max (B - 1))^t`. ∎

The coin is known at time 0, so this drift is exogenous, and Theorem 3 handles it exactly: given ξ = 1 the bound is 1, given ξ = 0 it is m / A, so `P(T ≤ m) ≤ p + (1 - p) m / A`. Proposition 4 gives `m / A + (1 + 1/(2A)) λ_max p (B - 1) s*`, about 28 times p at the default: the right order in `δ̄ m`, loose by the factor `λ_max (B - 1) s*`. The construction is the plug-in mechanism in miniature: a drift decided by something seen before monitoring starts. The mean drift is therefore not a sufficient summary across families of laws. Within the GARCH family E16 found that it orders the failures almost perfectly (Spearman 0.99), because the shape of the drift path changes little within the family.

## 4. The plug-in baseline

**Set-up.** A calibration window W (part of `F_0`) gives `τ̂`, `m̂` and `B̂ = τ̂² / m̂`; the decrease bets are `φ_l = s_l / (B̂ - 1)`. Monitoring uses `X_t = min(r_t², τ̂²)` and `Y_t = X_t / m̂`. Suppose that after calibration the conditional mean of `X_t` stays at a level `μ*(W)` that is fixed by the window, for example `μ* = m_F(τ̂) = E_F[min(r², τ̂²)]` for i.i.d. returns with law F (E16's τ-hat oracle bets against exactly this). More generally it is enough that `E[X_t | F_{t-1}] ≤ μ*↑` and `≥ μ*↓` with both known at time 0. Then the plug-in detector has the drift

`δ↑(W) = (μ*↑ / m̂ - 1)⁺`, `δ↓(W) = (1 - μ*↓ / m̂)⁺`,

which is decided by the window and constant afterwards. If the estimate is below the true baseline by a factor `1 / (1 + δ)`, that is `m̂ = μ* / (1 + δ)`, the increase side is exactly in the sup form with that δ.

**Corollary 1.** [PROVED] `P(T ≤ m | F_0) ≤ min{1, (m / A) I_m(δ↑(W), δ↓(W); B̂(W))}`, and therefore

`P(T ≤ m) ≤ β_m := E_W[min{1, (m / A) I_m(δ↑(W), δ↓(W); B̂(W))}]`,

where `I_m(·; B̂)` is the inflation factor of Theorem 1 with the decrease bets `s_l / (B̂ - 1)`. More generally, for any thresholds x, y ≥ 0, `P(T ≤ m) ≤ P(δ↑(W) > x or δ↓(W) > y) + (m / A) E_W[I_m(x, y; B̂(W))]`.

*Proof.* Theorem 1, part 3, with the drift and the bets known at time 0; the second form splits on the event `{δ↑ ≤ x, δ↓ ≤ y}` and uses that `I_m` increases in the drifts. ∎

**What the bound needs.** The law of `(δ↑(W), δ↓(W), B̂(W))` over calibration windows. For a given return law F this is a computation over windows only, with no detector (`robustness.plugin_bound`); section 5.5 does it for E16's laws. A distribution-free version needs, for every law in the class, a tail bound `P(δ↑(W) > x) ≤ ε_n(x)`. Two things stand in the way.

1. *Adaptedness.* In the monitor's plug-in, τ̂ and m̂ come from the same window, so the summands `min(r_i², τ̂²)` are not adapted and the confidence sequences of `docs/math.md` 6.2 do not apply to m̂ (6.4 there). Sample splitting (the cs baseline) or a band over all truncation levels (route A of `docs/theory_running_tau.md`) removes this obstacle.
2. *Precision.* Even with a valid upper confidence bound, the certified drift is far too large for a one-day guarantee: the one-day two-sided bound reaches 10% at a drift of 1.6% (section 5.1), while the best distribution-free upper bound E8 found for the truncated mean is 37% above the truth at 390 calibration steps, 17% at 1440 and 8% at 5000 (Gaussian returns, δ = 0.01, `docs/report_parts/E8.md`, Q8d). With the certified drift in place of δ the bound is 1.

So the distribution-free route is the confidence-sequence baseline of `docs/math.md` section 6, which moves the baseline up (down) by the half-width so that on an event of probability at least `1 - δ` the drift is at most 0: Corollary 1 with x = y = 0 on that event is Theorem 6.4, `P ≤ δ + m / A`. A plug-in with a margin, `m̂ (1 + x)`, sits in between: its drift is positive only when the estimation error exceeds x.

**Connection to Proposition A.6 of `docs/theory_running_tau.md`.** A.6 proves that at a fixed calibration length no inequality `P(T ≤ m) ≤ C m / A` with a constant C can hold for all A and m: `liminf_{A→∞} P(T_A ≤ ⌈(log A)²⌉) ≥ p_{η,κ} > 0`. Corollary 1 is consistent with this and says how the constant must depend on m.

* For every window with a positive median of `|r|`, `1.90 ≈ 2 a² / (1 + a²) ≤ B̂ ≤ 2 a² ≈ 39.6` with `a = c × 1.4826` (c = 3), because τ̂ = a × median, at least half of the `|r_i|` are at least the median (so `m̂ ≥ median² / 2`) and at least half are at most the median (so `m̂ ≤ (median² + τ̂²) / 2`). Also `δ↑(W) ≤ B̂ - 1` because `μ* ≤ τ̂²`, `δ↓(W) ≤ 1`, and `φ_l δ↓ < 1 / (B̂ - 1)`. So for each fixed m, `I_m` is bounded and Corollary 1 gives `P(T ≤ m) ≤ C_m m / A` with `C_m = E_W[I_m(δ↑(W), δ↓(W); B̂(W))]`, a constant for each m that grows exponentially in m (Proposition 1).
* Along `m_A = ⌈(log A)²⌉`, for every window with `δ↑(W) > 0`, `(m_A / A) I_{m_A}(δ↑(W), ·) ≥ (m_A / A) G_{m_A}(λ_max δ↑) / 12 → ∞`, because `(1 + λ_max δ↑)^{m_A} / A = exp(m_A log(1 + λ_max δ↑) - log A) → ∞`. The same holds on the decrease side for windows with `δ↓(W) > 0`. So, by Fatou's lemma, the capped bound `β_{m_A}` has lower limit at least `P(δ↑(W) > 0 or δ↓(W) > 0)`: it stays bounded away from 0, as A.6 requires.
* The two results fit together. In the set-up of A.6 (i.i.d. returns with law F, `μ* = m_F(τ̂)`), Jensen's inequality gives `g↑(W) ≤ max_λ log(1 + λ (m_F(τ̂) / m̂ - 1))`, which is positive only if `δ↑(W) > 0`, and likewise `g↓(W) > 0` only if `δ↓(W) > 0`. Hence `P(g(W) > 0) ≤ P(δ↑(W) > 0 or δ↓(W) > 0)`: the lower limit that A.6 proves for the probability is never above the lower limit of the upper bound of Corollary 1. A.6 needs a positive growth rate under F, which takes a drift larger than about `λ_min Var(Y) / 2` (the heuristic in `docs/theory_running_tau.md` 3.7). Corollary 1 covers every stream with the same conditional mean, including the noise-free stream of E16's shape adversary, for which any positive drift gives positive growth; section 5.5 shows that this stream comes within a factor of 1.3 of the bound at one day.

## 5. Numerical check

All numbers in this section are **[SYNTHETIC]** or exact evaluations of the closed forms, from `experiments/th4_robustness.py` (seeds in its docstring). A = 7800, one day = 390 steps (bound m / A = 5%), five days = 1950 steps (25%), Gaussian B = 9.045 unless stated. Rates come with 95% Wilson intervals in the CSV files. "Above the bound" would mean that the lower end of the interval is above it; this never happens.

### 5.1 The design guide: Theorem 1 in numbers

The inflation factor `I_m` of Theorem 1, for a drift of δ on both sides (`δ↑ = δ↓ = δ`), on the increase side only and on the decrease side only (`results/th4_sup_inflation.csv`; the bound is `min(1, (m / A) I_m)`; the last column is the first horizon at which `(m / A) I_m ≥ 1`):

| drift δ | one day: both sides | one day: increase only | one day: decrease only | five days: both sides | steps until the bound reaches 1 (both sides) |
|---|---|---|---|---|---|
| 0.05% | 1.013 | 1.010 | 1.003 | 1.072 | 6,072 |
| 0.10% | 1.027 | 1.021 | 1.006 | 1.160 | 4,880 |
| 0.20% | 1.056 | 1.044 | 1.013 | 1.412 | 3,483 |
| 0.50% | 1.160 | 1.127 | 1.033 | 4.239 | 1,919 |
| 1.00% | 1.412 | 1.342 | 1.070 | 158 | 1,140 |
| 2.00% | 2.638 | 2.480 | 1.158 | 1.2 × 10^6 | 654 |
| 3.00% | 7.452 | 7.180 | 1.273 | 1.2 × 10^10 | 468 |
| 5.00% | 147 | 146 | 1.619 | 1.4 × 10^18 | 305 |

The first-order formula of Proposition 1 is accurate while the factor is close to 1: at δ = 0.1%, `1 + 391 × (0.2033 + 0.0634) × 0.001 / 4 = 1.026` against the exact 1.027.

**Design drifts.** The drift at which the one-day bound of the default detector reaches 10% (twice the nominal 5%) and 25% (five times), found by root finding on the closed form (`results/th4_design.csv`):

| | 10% in one day | 25% in one day |
|---|---|---|
| both sides drift by δ | 1.63% | 2.66% |
| upward drift only | 1.72% | 2.70% |
| downward drift only | 6.36% | 10.51% |
| both sides, 1440-step day with A = 28,800 (same 5% per day) | 0.44% | 0.72% |

The downward drift costs about four times less because the decrease bets are smaller (section 2.2). A finer clock with the same budget per day is less robust: the 1440-step day at the same 5% tolerates only about a quarter of the drift (upward only 0.47% and 0.73%; downward only 1.72% and 2.84%), because the drift compounds over 1440 steps instead of 390.

### 5.2 Constant drift: is Theorem 1 tight?

Streams with a known constant drift, `E[Y_t | past] = 1 + δ` (increase side) or `1 - δ` (decrease side) at every step: a deterministic stream (`Y_t` equal to its mean), the i.i.d. two-point law {0, B} (the widest spread), Gaussian returns with the variance scaled so that the truncated mean is exactly `m₀ (1 ± δ)`, 10,000 paths each, and an adaptive adversary with the same drift: hold `Y_t` at its mean for h steps (every bet then gains without risk), then play bold as in E16 (the points 0 and B, or the exact point that lands a side on 2A), with h searched over 15 values on 4,000 paths and the best h rerun on 20,000 fresh paths (`results/th4_constant_drift.csv`, figure `results/figures/th4_sup_bound.png`).

| side | δ | horizon | Theorem 1 bound | deterministic | two-point {0, B} | Gaussian, variance shifted | hold, then bold play | bound / adversary |
|---|---|---|---|---|---|---|---|---|
| increase | 0.25% | 390 | 5.3% | 0 | 2.0% | 1.1% | 2.3% (hold 249) | 2.33 |
| increase | 0.25% | 1950 | 37.7% | 0 | 14.8% | 12.6% | 24.7% (hold 1663) | 1.53 |
| increase | 1.00% | 390 | 6.7% | 0 | 1.9% | 1.1% | 2.7% (hold 305) | 2.48 |
| increase | 1.00% | 1950 | 1 | 1 | 15.1% | 12.9% | 100.0% (hold 1247) | 1.00 |
| increase | 2.00% | 390 | 12.4% | 0 | 1.9% | 1.5% | 6.6% (hold 360) | 1.88 |
| increase | 3.00% | 390 | 35.9% | 0 | 2.0% | 1.7% | 22.7% (hold 360) | 1.58 |
| increase | 5.00% | 390 | 1 | 1 | 2.4% | 2.7% | 100.0% (hold 305) | 1.00 |
| decrease | 1.00% | 390 | 5.3% | 0 | 1.8% | 0.7% | 2.4% (hold 55) | 2.21 |
| decrease | 1.00% | 1950 | 40.5% | 0 | 15.8% | 13.4% | 22.3% (hold 1802) | 1.81 |
| decrease | 2.00% | 1950 | 1 | 0 | 16.0% | 15.9% | 94.8% (hold 1940) | 1.06 |
| decrease | 5.00% | 390 | 8.1% | 0 | 1.9% | 1.4% | 2.7% (hold 222) | 3.03 |
| decrease | 10.00% | 390 | 21.9% | 0 | 3.5% | 4.9% | 12.4% (hold 360) | 1.76 |

(The CSV file has all 80 cells: δ = 0.25%, 0.5%, 1%, 2%, 3%, 5% upward and 1%, 2%, 5%, 10% downward, at both horizons.)

**Reading.** No measured rate is above the bound, in any of the 80 cells. The deterministic stream attains the bound where Proposition 2 says it should, at the drift where the bound reaches 1 (`results/th4_sharpness.csv`): the one-day bound reaches 1 at an upward drift of 3.740% and the deterministic stream alarms within the day from 3.756% on (downward: 15.01% and 15.08%; over five days 0.499% and 0.520% upward, 1.925% and 2.016% downward; figure, top left). Below that drift the deterministic stream does not alarm at all within the horizon (for δ = 3% upward it alarms at step 472), so it says nothing about the bound there. At shorter horizons than the deterministic alarm step, the adaptive adversary is the strongest stream we found: at one day it reaches 0.33 to 0.63 of the bound wherever the bound is below 1, and over five days 0.55 to 0.65 where the bound is below 1 (it alarms on 85% to 100% of paths where the bound is 1). Noisy i.i.d. streams are far below the bound: with the widest spread {0, B} the one-day rate barely moves with the drift (2.0% at 0.25%, 2.4% at 5%), because the noise of `Y_t` drowns a small drift within a day; over five days the drift shows (Gaussian: 12.6% at 0.25%, 31.9% at 5%). So Theorem 1 is the right worst case over all streams with a given drift, and very conservative for noisy ones; the worst case at horizons shorter than the deterministic alarm step is not known exactly (section 6).

### 5.3 An exogenous drift that averages out: sup form against the pathwise forms

The volatility switches between two levels whose truncated means are `m₀ (1 + ε)` and `m₀ (1 - ε)`, by a symmetric Markov chain independent of the Gaussian innovations (probability `stay` of keeping the state), so the drift is `±ε` and averages to 0 over time. This is the third example of section 3.3, so Theorems 1, 2 and 3 and Proposition 4 all apply. 6,000 paths per cell (`results/th4_switching.csv`). The pathwise bounds are expectations over the drift paths, evaluated on the same paths.

| ε | stay | half-life (steps) | one day: measured | sup form | Theorem 3 | Theorem 2 | five days: measured | sup form | Theorem 3 | Theorem 2 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1% | 0.5 | 1 | 0.75% | 7.1% | 5.6% | 7.1% | 12.3% | 1 | 32.7% | 44.5% |
| 1% | 0.99 | 69 | 0.90% | 7.1% | 5.9% | 30.2% | 12.7% | 1 | 87.0% | 1 |
| 2% | 0.5 | 1 | 0.98% | 13.2% | 6.3% | 9.2% | 12.0% | 1 | 42.3% | 64.8% |
| 2% | 0.99 | 69 | 1.05% | 13.2% | 7.6% | 64.5% | 13.4% | 1 | 100.0% | 1 |
| 5% | 0.5 | 1 | 1.20% | 1 | 8.9% | 17.9% | 13.0% | 1 | 73.5% | 1 |
| 5% | 0.9 | 7 | 1.32% | 1 | 21.9% | 61.4% | 13.0% | 1 | 100.0% | 1 |
| 5% | 0.999 | 693 | 1.77% | 1 | 52.1% | 1 | 27.0% | 1 | 1 | 1 |
| 10% | 0.5 | 1 | 0.73% | 1 | 15.3% | 34.5% | 12.2% | 1 | 96.0% | 1 |
| 10% | 0.99 | 69 | 3.28% | 1 | 80.6% | 1 | 32.9% | 1 | 1 | 1 |
| 10% | 0.999 | 693 | 5.25% | 1 | 66.0% | 1 | 67.0% | 1 | 1 | 1 |

(The CSV file has all 16 cells: ε = 1%, 2%, 5%, 10% and stay = 0.5, 0.9, 0.99, 0.999. Proposition 4 gives 1 in every cell.)

**Reading.** The sup form uses the largest drift ε and ignores that it averages out; it is below 1 in 8 of the 16 cells at one day and in none at five days. The pathwise bound of Theorem 3 is below 1 in all 16 cells at one day and in 11 at five days, and it is always at least as good as the sup form. When the state flips at every step (stay = 0.5) a drift of 10% costs almost nothing: the measured one-day rate is 0.73%, and Theorem 3 gives 15.3% against 1 for the sup form. The more persistent the drift, the closer it acts to a constant drift, as section 3.4 (b) predicts: at ε = 10% the one-day rate rises from 0.73% (stay 0.5) to 5.25% (stay 0.999, above the nominal 5%), and the five-day rate from 12.2% to 67.0%. Theorem 2, which also covers a drift that responds to the data, pays for that generality: it is never below Theorem 3 here, and it is below 1 in 11 of the 16 cells at one day (every cell with ε of at most 2%, and the cells with stay 0.5 or 0.9 at 5% and stay 0.5 at 10%) and in 3 at five days (stay 0.5 or 0.9 at 1%, stay 0.5 at 2%). The drift-corrected detector of Lemma 2 stays within its bound m / A in every cell, and the pathwise inequality `N_t ≤ D*_t Ñ_t` held at every step of every path (largest ratio 0.9989).

### 5.4 GARCH clustering (E16 Q16d): bound against measured rate

The 64 GARCH(1,1) cells of E16, with E16's parameters, τ and m₀ (from `results/e16b_garch.csv`), on 6,000 fresh paths per cell. The exact drift `d_t = σ_t² h(τ / σ_t) / m₀ - 1` is known (Gaussian innovations), which gives the pathwise quantities on the same paths. Because the drift responds to the data, only Theorem 2, Proposition 2b and Proposition 4 are proved here; "best proved" is the smallest of the three. The threshold K of Theorem 2 is chosen on the same paths, which makes its value slightly optimistic; the bounds are Monte Carlo estimates of proved quantities. Two heuristics are also reported: Theorem 3's pathwise form, which section 3.4 (a) shows is not proved for GARCH, and the sup form at E16's stored mean absolute drift (`drift_mean_abs`), which Proposition 5 shows is not a bound in general. The measured rate on the new paths agrees with E16's rate within Monte Carlo error (`results/th4_garch.csv`, figure `results/figures/th4_e16_comparison.png`). The table shows 20 of the 64 cells.

| α + β | α | mean abs. drift | drift 99% quantile | one day: measured (E16) | best proved | Theorem 3 (heuristic) | sup form at mean drift (heuristic) | five days: measured (E16) | best proved | Theorem 3 (heuristic) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.002 | 0.23% | 1.2% | 0.97% (0.92%) | 6.2% | 5.2% | 5.3% | 13.3% (12.7%) | 35.8% | 28.6% |
| 0.5 | 0.005 | 0.58% | 2.9% | 1.03% (1.18%) | 7.9% | 5.4% | 6.0% | 14.4% (13.4%) | 54.0% | 36.1% |
| 0.5 | 0.01 | 1.15% | 5.8% | 1.25% (1.48%) | 11.5% | 6.0% | 7.6% | 14.6% (14.9%) | 82.7% | 52.1% |
| 0.5 | 0.02 | 2.29% | 11.7% | 1.62% (1.53%) | 20.2% | 7.5% | 17.0% | 18.1% (17.4%) | 1 | 80.4% |
| 0.5 | 0.05 | 5.64% | 29.3% | 2.90% (3.02%) | 51.4% | 23.3% | 1 | 26.7% (27.1%) | 1 | 99.8% |
| 0.9 | 0.002 | 0.49% | 1.9% | 1.53% (1.38%) | 11.3% | 5.4% | 5.8% | 14.8% (15.7%) | 81.9% | 39.7% |
| 0.9 | 0.005 | 1.23% | 4.7% | 1.53% (1.52%) | 25.0% | 6.1% | 7.9% | 19.2% (18.0%) | 1 | 87.9% |
| 0.9 | 0.01 | 2.45% | 9.5% | 2.42% (2.90%) | 51.1% | 8.5% | 20.0% | 24.2% (25.3%) | 1 | 99.6% |
| 0.9 | 0.02 | 4.87% | 19.4% | 5.30% (5.48%) | 92.2% | 30.5% | 1 | 39.6% (41.4%) | 1 | 1 |
| 0.9 | 0.05 | 11.87% | 51.2% | 19.40% (19.73%) | 1 | 87.5% | 1 | 81.4% (81.8%) | 1 | 1 |
| 0.98 | 0.002 | 1.10% | 3.7% | 2.02% (1.78%) | 40.2% | 6.1% | 7.4% | 21.6% (21.8%) | 1 | 87.4% |
| 0.98 | 0.005 | 2.75% | 9.4% | 3.95% (3.92%) | 79.0% | 15.3% | 27.6% | 35.9% (35.9%) | 1 | 100.0% |
| 0.98 | 0.01 | 5.47% | 19.3% | 9.43% (8.78%) | 1 | 41.9% | 1 | 60.8% (60.4%) | 1 | 1 |
| 0.98 | 0.02 | 10.77% | 40.8% | 22.40% (25.07%) | 1 | 72.8% | 1 | 90.2% (91.0%) | 1 | 1 |
| 0.98 | 0.05 | 25.57% | 116.3% | 66.70% (66.47%) | 1 | 99.8% | 1 | 99.9% (100.0%) | 1 | 1 |
| 0.999 | 0.002 | 4.83% | 15.7% | 5.58% (5.22%) | 1 | 36.9% | 1 | 51.0% (52.4%) | 1 | 100.0% |
| 0.999 | 0.005 | 12.30% | 45.6% | 25.72% (24.93%) | 1 | 65.2% | 1 | 86.8% (86.8%) | 1 | 1 |
| 0.999 | 0.01 | 23.40% | 93.5% | 52.15% (52.00%) | 1 | 87.0% | 1 | 98.2% (98.2%) | 1 | 1 |
| 0.999 | 0.02 | 41.84% | 204.7% | 79.22% (79.13%) | 1 | 98.5% | 1 | 100.0% (100.0%) | 1 | 1 |
| 0.999 | 0.05 | 64.93% | 331.1% | 95.63% (96.35%) | 1 | 1 | 1 | 100.0% (100.0%) | 1 | 1 |

**Over all 64 cells.**

| | one day | five days |
|---|---|---|
| cells above the nominal bound m / A (E16) | 43 | 51 |
| sup form (Theorem 1 with the true largest drift, B - 1) below 1 | 0 | 0 |
| best proved bound below 1 | 22 | 6 |
| best proved bound / measured rate, where below 1 | 6.4 to 29 | 2.7 to 6.3 |
| proved bound below the measured rate (lower end of the interval) | 0 | 0 |
| Theorem 3 (heuristic) below 1 | 44 | 19 |
| Theorem 3 (heuristic) / measured rate, where below 1 | 1.2 to 8.3 | 2.0 to 4.9 |
| sup form at the mean drift (heuristic) below 1 | 18 | 2 |
| either heuristic below the measured rate (lower end of the interval) | 0 | 0 |
| drift-corrected detector (Lemma 2), largest rate | 1.22% | 14.4% |

**Reading.** The sup form is useless for GARCH: the conditional mean comes arbitrarily close to `B m₀`, so the only valid constant is `δ↑ = B - 1` and Theorem 1 gives 1 in every cell. The pathwise form of Theorem 2 is what makes a proved statement possible at all: it is below 1 in 22 cells at one day, all with a mean absolute drift of at most 11% (and up to 1.15% at five days), because a calm spell resets the window products. It is loose there (6 to 29 times the measured rate at one day), for two reasons visible in section 5.2: a bound valid for every stream with the same drift path must cover streams with no noise, and GARCH returns are noisy; and the threshold form pays `K m / A` with K between 1.2 and 7.8. Only 2 of the 43 cells that E16 found above the nominal bound at one day get a proved bound below 1 ((0.5, 0.1) with 6.7% measured against 93.6%, and (0.99, 0.005) with 5.5% against 97.7%), and none of the 51 at five days: the proved bounds certify the cells where the drift is small, which are also the cells where the nominal bound already holds. Truncated sup (Proposition 2b) and mean drift (Proposition 4) were never the best proved bound: at one day Proposition 2b is below 1 in only 7 cells and Proposition 4 in 5, and at five days neither is below 1 in any cell. Theorem 3, used as a heuristic, is much closer (1.2 to 8.3 times the measured rate) and was above every measured rate; section 3.4 (a) shows why this cannot be promoted to a proof (the drift is high exactly on the paths where the wealth is already high). The sup form at E16's mean drift, also a heuristic, was above every measured rate as well, which agrees with E16's finding that within the GARCH family the mean drift orders the failures. The pathwise inequality of Lemma 2 held at every step of every path (largest ratio `N_t / (D*_t Ñ_t)` 0.99999994, so on some steps it is close to an equality). The drift-corrected detector, which is E16's control in multiplicative form, stayed below m / A everywhere.

### 5.5 The plug-in baseline (E16 Q16e): Corollary 1 against measured rates

For each law and calibration length of E16's plug-in cells, 100,000 fresh calibration windows give the drift `δ(W) = m(τ̂) / m̂ - 1` (with `m(τ̂)` the true mean of the law at the estimated truncation level) and `B̂(W)`; the bound `β_m` of Corollary 1 is the average of the capped per-window bound. Under i.i.d. returns and under E16's shape adversary (which keeps the conditional mean at `m(τ̂)`) the assumption of Corollary 1 holds exactly (`results/th4_plugin.csv`).

| calibration | monitored law | horizon | E16 measured (95% interval) | Corollary 1 bound | bound / measured |
|---|---|---|---|---|---|
| 390 | Gaussian, i.i.d. | 390 | 4.0% (3.8% to 4.2%) | 41.8% | 10.5 |
| 390 | Student t (2.5), i.i.d. | 390 | 5.0% (4.7% to 5.2%) | 55.9% | 11.2 |
| 390 | Gaussian, i.i.d. | 1950 | 41.8% (41.2% to 42.4%) | 92.4% | 2.21 |
| 1440 | Gaussian, i.i.d. | 390 | 1.6% (1.5% to 1.7%) | 26.1% | 16.4 |
| 1440 | Gaussian, i.i.d. | 1950 | 23.1% (22.6% to 23.5%) | 86.1% | 3.73 |
| 390 | Gaussian window, shape adversary q = 0.2 | 390 | 14.8% (14.4% to 15.2%) | 41.8% | 2.83 |
| 390 | Gaussian window, shape adversary q = 0 | 390 | 32.5% (32.0% to 33.1%) | 41.8% | 1.29 |
| 390 | Gaussian window, shape adversary q = 0 | 1950 | 85.8% (85.4% to 86.2%) | 92.4% | 1.08 |
| 1440 | Gaussian window, shape adversary q = 0 | 390 | 16.1% (15.7% to 16.6%) | 26.1% | 1.61 |
| 1440 | Gaussian window, shape adversary q = 0 | 1950 | 73.1% (72.6% to 73.6%) | 86.1% | 1.18 |

(The CSV file has all 36 cells: four i.i.d. families at each calibration length and horizon, and the shape adversary at q = 1, 0.5, 0.2, 0.05, 0). In E16's shape adversary, q = 1 is the Gaussian law itself and q = 0 the noise-free stream `Y_t = 1 + δ(W)`.)

**Reading.** The bound lies above every measured rate, in all 36 cells, and also within each decile of the baseline error that E16 reports (80 cells, `results/th4_plugin_deciles.csv`, with the bound recomputed on the fresh windows whose error falls in the same range). It is nearly attained by the noise-free shape adversary (1.08 to 1.61 times its rate) and 10 to 16 times above the i.i.d. rate at one day, for the same reason as in section 5.2. So for the plug-in, Corollary 1 is close to the worst case over all streams with the plug-in's conditional mean, and E16's shape adversary is close to the worst such stream. In the deciles (390 calibration steps, i.i.d. Gaussian, one day), the bound is 5.3% to 11.2% in deciles 5 to 9, where E16 measured 0.8% to 2.2%. It is 1 in the three deciles where m̂ is most too low, where E16 measured 2.8% to 16.0%, and 37% in the decile where m̂ is most too high, where E16 measured 7.5%.

**How long must the calibration be?** The bound for the Gaussian plug-in as the calibration grows (20,000 windows each, `results/th4_plugin_caln.csv`):

| calibration steps | sd of the drift | windows with abs. drift above 1.63% | one-day bound | five-day bound |
|---|---|---|---|---|
| 390 | 7.08% | 81.5% | 42.1% | 92.4% |
| 1,440 | 3.65% | 65.3% | 25.7% | 85.8% |
| 2,880 | 2.58% | 52.4% | 17.9% | 80.4% |
| 5,760 | 1.82% | 37.4% | 11.1% | 74.1% |
| 11,520 | 1.28% | 20.4% | 7.3% | 65.7% |
| 23,040 | 0.91% | 7.4% | 6.0% | 57.7% |

The one-day bound falls like the share of windows whose drift exceeds the one-day design drift of 1.63%, and comes close to the nominal 5% only after about 23,000 calibration steps (59 days of 390 steps). The five-day bound stays above 50% at every length tried. Over five days the bound stays within 1.42 times the nominal 25% only for a drift of at most 0.2% (section 5.1), which the plug-in does not reach. This is the quantitative form of Proposition A.6 (section 4).

### 5.6 Price-level ticks (E16 Q16f): what is missing

The theorems need the drift given the whole past. For prices rounded to a tick, that is a filtering problem: the conditional mean depends on where the unobserved efficient price sits inside its tick, which the whole history of rounded prices informs. E16 stored the conditional mean given the previous return only. That is an average of the drift given the whole past, so its range lies inside the true range, and Theorem 1 evaluated at it (`results/th4_ticks.csv`) is a heuristic that can be too small, not a bound. For what it is worth, it was above the measured rate at every tick and horizon. It is informative (below 1) only for ticks of at most 1 sd at one day (at most 5.9%, against measured rates of 0.7% to 1.5%) and at most 0.75 sd over five days (26% to 40%, against 12% to 14%). From a tick of 1.5 sd the stored upward drift after a move of two or more ticks is 4.3% or more, and the heuristic is 1. A proved bound for ticks would need the drift given the whole past, for example by filtering the position of the efficient price, and is left open.

### 5.7 What the tests check

`tests/test_robustness.py` (21 tests): the closed form and the elementary bounds of Proposition 1; reduction to Theorem B at δ = 0; Lemma 1 by Monte Carlo with an arbitrary stopping time; Theorem 1 by simulation on streams with a known constant drift (two laws, increase and decrease sides, at a small A where the bound is below 1 and the drift matters); the deterministic stream reaching the bound (Proposition 2); the vectorised path engine against the package's own detector and alarm rule; the pathwise inequality of Lemma 2; the CUSUM and Shiryaev-Roberts recursions against brute force; Theorems 2 and 3 and Propositions 2b and 4 on a switching drift, with the bounds below 1; the two-step counterexample, exactly and by Monte Carlo; the burst of Proposition 5 (s* = 7, the rate equals the share of bursts, the mean drift identity); Corollary 1 per window against the noise-free stream.

## 6. Status

| Claim | Status | Checked by |
|---|---|---|
| Lemma 1: under a drift, the expected stopped Shiryaev-Roberts statistic is at most the drift's Shiryaev-Roberts sum | **[PROVED]** | Monte Carlo test with a stopping time |
| Theorem 1 (sup form): `P(T ≤ m) ≤ (m / A) I_m(δ↑, δ↓)`, two-sided and one side, closed form of `I_m` | **[PROVED]** | simulation tests with a known constant drift (two laws, both sides); 80 cells of section 5.2 |
| Proposition 1: `I_m` is close to 1 when `x* m` is small; exponential growth when it is large; first-order formula | **[PROVED]** | tests of the elementary bounds |
| Proposition 2: Theorem 1 is attained by the deterministic stream at its alarm step | **[PROVED]** | test; section 5.2 and `th4_sharpness.csv` (the bound reaches 1 within 0.02 points of drift of where the stream alarms, at one day) |
| Lemma 2: the drift-corrected detector is valid, and `N_t ≤ D*_t Ñ_t` path by path | **[PROVED]** | test; held at every step of every simulated path (sections 5.3, 5.4) |
| Theorem 2 (threshold form, any predictable drift, including one that responds to the data) | **[PROVED]** | tests; above every measured rate (sections 5.3, 5.4) |
| Proposition 2b (truncated sup form) | **[PROVED]** | tests; section 5.4 |
| Theorem 3 (pathwise form, exogenous drift) | **[PROVED]** | tests; section 5.3 |
| 3.4 (a): the expected product does not bound the expected wealth when the drift responds to the data | **[PROVED]** (counterexample) | exact and Monte Carlo test (1.75 against 1.5) |
| 3.4 (b): the mean drift gives a lower estimate of the pathwise bound, not an upper bound | **[PROVED]** (Jensen) | section 5.3 (persistence raises the rate at a fixed drift size) |
| Proposition 4 (mean drift, with an additive term weighted by A) | **[PROVED]** | tests; section 5.4 |
| Proposition 5: no bound `C(δ̄) m / A` in terms of the mean drift alone | **[PROVED]** (counterexample) | test (s* = 7, rate equals the share of bursts) |
| Corollary 1 (plug-in baseline): average over calibration windows of the capped Theorem 1 bound | **[PROVED]** | test; 36 cells and 80 decile cells of section 5.5 |
| Consistency with Proposition A.6 (section 4) | **[PROVED]** | |
| Design drifts (10% and 25% in one day) and the inflation table | **[NUMERICAL]** (root finding and evaluation of a proved closed form) | `th4_design.csv`, `th4_sup_inflation.csv` |
| How close Theorem 1 is to the worst case at horizons shorter than the deterministic alarm step | **[NUMERICAL]**: the best adaptive stream found reaches 0.33 to 0.63 of it at one day; the exact worst case is **open** | section 5.2 |
| Theorem 3's pathwise form applied to GARCH | **[HEURISTIC]**: above every measured rate in all 128 cells, 1.2 to 8.3 times the rate where below 1 at one day; not a proof (3.4 (a)) | section 5.4 |
| Theorem 1 at E16's mean absolute drift, for GARCH | **[HEURISTIC]**: above every measured rate in all 128 cells; not a bound in general (Proposition 5) | section 5.4 |
| Theorem 1 at E16's one-step drift, for price ticks | **[HEURISTIC]**: above every measured rate; a proved bound needs the drift given the whole past (filtering), **open** | section 5.6 |
| A distribution-free bound for the plug-in baseline | **Open**: needs a tail bound on the baseline error for every law in the class and an adapted estimate (section 4); E8's confidence-sequence baseline is the proved alternative | |
| A proved bound for GARCH that is informative when E16 sees the nominal bound fail (mean drift of a few percent) | **Open**: Proposition 5 rules out a bound in the mean drift alone in general, but not within a parametric family such as GARCH; the proved pathwise bounds certify only the cells with a small drift (22 of 64 at one day) | section 5.4 |

## 7. References

* Ville, J. (1939). *Étude critique de la notion de collectif.* Thèse, Paris. Ville's inequality, as cited in `docs/math.md` (6.7), with the modern statement of Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2020), Time-uniform Chernoff bounds via nonnegative supermartingales, *Probability Surveys* 17, 257 to 317.
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science* 2(2), 229 to 260 (the Shiryaev-Roberts e-detector, as cited in `docs/math.md`).

Everything else (Lemmas 1 and 2, Theorems 1 to 3, Propositions 1 to 5, Corollary 1) is proved here from these and from elementary inequalities (arithmetic and geometric means, Bernoulli's inequality, Jensen's inequality, the series of sinh).
