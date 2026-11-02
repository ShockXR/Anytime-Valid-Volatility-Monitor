# A false alarm guarantee with a running truncation level

This note extends `docs/math.md` (sections 5 and 6) to monitors whose truncation level τ_t changes at every step. Code: `src/avmonitor/running.py`. Tests: `tests/test_running.py`. Experiments: `experiments/e14_running_tau.py`, report `docs/report_parts/E14.md`.

## 0. Summary in plain words

The detector of `docs/math.md` bets on the truncated squared return `X_t = min(r_t², τ²)` against its null mean `m(τ)`. Its guarantee, P(false alarm within m steps) ≤ m / A, needs the true `m(τ)` for one fixed τ. Section 6 of `docs/math.md` (E8) keeps the guarantee with an estimated mean, at the price of a small failure probability δ, but only for a τ fixed before the estimation starts. Real data want a τ that follows the current scale. A rule that sets τ_t from past data only is called *predictable*. This note gives two ways to let τ_t be any predictable rule and keep a proof.

* **Route A (a band for all truncation levels at once).** If the calibration returns are independent draws from one law F, estimate the whole function `s ↦ m_F(s)` with a band `L(s) ≤ m_F(s) ≤ U(s)` that holds for every s simultaneously with probability at least 1 - δ. Then bet against `U(τ_t)` (increase side) and `L(τ_t)` (decrease side), whatever τ_t is. Result: P(false alarm within m) ≤ δ + m / A (Theorem A.1). Two bands are implemented: the time-uniform Dvoretzky-Kiefer-Wolfowitz (DKW) band of Howard and Ramdas (2022), integrated (simple, but wide), and a band made of empirical Bernstein confidence sequences on a grid of truncation levels (much tighter).
* **Route B (ranks).** If the returns are *exchangeable* (their joint law does not change when their order is permuted; i.i.d. returns with any law, including ties and jumps, are exchangeable), the rank of `|r_t|` among all earlier values, randomised to break ties, is an exactly uniform p-value, and these p-values are independent (Proposition B.1, a classical result of Vovk and coauthors). Betting on them needs no baseline, no δ and no truncation level, and the guarantee P(false alarm within m) ≤ m / A is exact (Theorem B.2). The price is that post-change returns enter the reference set and dilute the evidence (section 4.3).
* **No blind period.** Both routes allow monitoring to resume immediately after an alarm instead of after a full calibration day: route A with a band that keeps learning and Shiryaev-Roberts chains that each bet against the band available when they start (Theorem A.4); route B simply restarts its reference set (Remark B.3). Corollary A.5 bounds the expected number of false alarms when the monitor restarts after every alarm.

* **Why the plug-in is not enough.** At a fixed calibration length, the plug-in detector has a positive drift under the null on a fixed share of calibration windows, so its run length grows only like log A there and no bound of the form C m / A can hold for all A (Proposition A.6). E14 sees it at A = 7800 already: within five days (1950 steps) the plug-in alarms on 41% to 43% of null paths at 390 calibration steps, above the bound 1950 / 7800 = 25% that the oracle (12% to 14%) and the rank detectors (14% to 18%) respect.

What is proved, what is cited and what is only checked numerically is listed in section 8.

## 1. Set-up and notation

* Returns `r_1, r_2, ...`. `F_t` is the information at time t: `r_1, ..., r_t` and any independent randomisation used up to t (tie-breakers, for example). A random quantity known at time t - 1 is *predictable*.
* For a law F of a return, `S_F(x) = P_F(r² > x)` (x ≥ 0) and the truncated second moment at level s > 0 is
  `m_F(s) = E_F[min(r², s²)]`.
* Truncation level τ_t > 0, predictable; `X_t = min(r_t², τ_t²) ∈ [0, τ_t²]`.
* Monitoring starts after the end of calibration T (a fixed time or a stopping time). Bets: increase side `f↑_t = 1 + λ (X_t / m↑_t - 1)` with λ from the fixed grid `[0.02, 0.05, 0.1, 0.2, 0.35, 0.5]` (`detectors.DEFAULT_UP_GRID`); decrease side `f↓_t = 1 + λ_t (1 - X_t / m↓_t)` with `λ_t = φ / (B_t - 1)`, `B_t = τ_t² / m↓_t`, φ from `[0.1, 0.25, 0.5, 0.75, 0.95]`, and `λ_t = 0` when `m↓_t ≤ 0` or `B_t ≤ 1`. Here `m↑_t`, `m↓_t` are predictable.
* Shiryaev-Roberts statistic `M_t = f_t (M_{t-1} + 1)`, `M_T = 0`, averaged over the bet grid; two-sided alarm when `max(M↑_t, M↓_t) ≥ 2A`. When all chains share the bet at time t, `M_t = Σ_{j ≤ t} Π_{i=j}^{t} f_i` is the sum of the wealths of chains started at every j (`docs/math.md`, 5.4).
* `docs/math.md` Lemma 1 and Theorems A and B: if `E[f_t | F_{t-1}] ≤ 1` for every t, then `M_t - t` is a supermartingale, P(τ_A ≤ m) ≤ m / A and E[τ_A] ≥ A. The two-sided rule inherits both because `(M↑ + M↓) / 2` is again such a statistic.

## 2. Facts about truncated means

**Lemma 2.1.** For any law of r and s > 0:

1. `m(s) = ∫_0^{s²} S(x) dx`.
2. `m` is nondecreasing in s, and `m(s) / s²` is nonincreasing in s.
3. If `S_G(x) ≤ S_F(x)` for every x ≥ 0 (|r| is stochastically smaller under G), then `m_G(s) ≤ m_F(s)` for every s.
4. Conditional form: if s is `F_{t-1}`-measurable, then `E[min(r_t², s²) | F_{t-1}] = ∫_0^{s²} P(r_t² > x | F_{t-1}) dx` almost surely.

*Proof.* (1) For V ≥ 0 and a ≥ 0, `min(V, a) = ∫_0^a 1{V > x} dx`; take expectations (Tonelli). (2) The integrand is nonnegative, so m is nondecreasing. `min(r², s²) / s² = min(r² / s², 1)` is nonincreasing in s for every r, so its expectation is too. (3) Integrate the inequality over `[0, s²]`. (4) The identity in (1) holds pointwise for `a = s²`; take conditional expectations, using a regular conditional distribution of r_t given `F_{t-1}` and the fact that s is fixed given `F_{t-1}`. ∎

Part 2 gives bounds between grid points. If `s_k < s < s_{k+1}`:

`max(m(s_k), m(s_{k+1}) s² / s_{k+1}²) ≤ m(s) ≤ min(m(s_{k+1}), m(s_k) s² / s_k²)`. (2.1)

## 3. Route A: a band for all truncation levels

### 3.1 The null hypotheses

Let F be a law of returns, unknown to the analyst.

* **(C) Calibration.** The calibration returns `r_1, ..., r_T` are i.i.d. with law F.
* **(H↑_b) Increase null.** For every `t > T` and every s > 0, `E[min(r_t², s²) | F_{t-1}] ≤ (1 + b) m_F(s)`.
* **(H↓_b) Decrease null.** For every `t > T` and every s > 0, `E[min(r_t², s²) | F_{t-1}] ≥ (1 - b) m_F(s)`.

Here `b ∈ [0, 1)` is the null band, fixed in advance. By Lemma 2.1(3, 4), (H↑_0) holds whenever the conditional law of `|r_t|` given the past is stochastically no larger than F, for example `r_t = σ_t ε_t` with i.i.d. innovations ε_t independent of the past, `σ_t ≤ σ_0` predictable, and calibration returns `σ_0 ε_i`. (H↑_0) and (H↓_0) together hold if all returns are i.i.d. with law F. With b = 0 they are the natural "no change" hypotheses; they restrict only truncated conditional means, not the dependence between returns.

### 3.2 Bands

**Definition 3.1.** A *band at time s for F* is a pair of functions `U_s, L_s : (0, ∞) → [0, ∞]`, computed from `r_1, ..., r_s` (and the randomisation up to s). Its *failure event* is `D_s = {there is a level u > 0 with m_F(u) > U_s(u) or m_F(u) < L_s(u)}`. A sequence of bands is *time-uniform at level δ* if `P(∪_{s ≥ s_0} D_s) ≤ δ`; a band at one fixed time T has level δ if `P(D_T) ≤ δ`.

If F is fixed, `D_s` is determined by `r_1, ..., r_s`, so `D_s ∈ F_s`, and the first failure time `N = inf{s ≥ s_0 : D_s occurs}` is a stopping time.

Only the levels actually used need to be covered. If τ_t is restricted to a finite grid `{s_1, ..., s_K}`, the failure event may be taken over the grid only. With (2.1) a grid band extends to every level in `[s_1, s_K]`.

### 3.3 The guarantee with a frozen band

**Monitoring rule (frozen band).** At the end of calibration T a band `(U, L) = (U_T, L_T)` is computed. For `t > T`, with any predictable τ_t: `m↑_t = (1 + b) U(τ_t)` and `m↓_t = (1 - b) L(τ_t)`, bets as in section 1.

**Theorem A.1.** Assume (C), (H↑_b) and (H↓_b), and let the band at T have level δ (with T a stopping time, let the band be time-uniform at level δ). Then for every m ≥ 1:

1. P(an alarm at one of the times T + 1, ..., T + m) ≤ δ + m / A;
2. the expected number of monitored steps before the first alarm is at least (1 - δ) A;
3. one side alone at level A: under (H↑_b), P(increase alarm within m) ≤ δ + m / A, and symmetrically for the decrease side.

*Proof.* Let `G = D_Tᶜ` (the band holds at T). G is known at time T, and `P(Gᶜ) ≤ δ`: for a fixed T by definition; for a stopping time T because `D_T ⊆ ∪_s D_s`.

*Step 1 (the bets are admissible on every path).* `X_t ≥ 0` and λ < 1, so `f↑_t ≥ 1 - λ > 0`. On the decrease side λ_t is computed from the same `m↓_t` and τ_t that enter the factor, and `X_t ≤ τ_t²`, so `f↓_t ≥ 1 - φ > 0`; a step with no admissible bet has factor 1. Admissibility therefore never depends on whether the band holds, which matters because the band failure is not observable.

*Step 2 (on G, every factor has conditional mean at most 1).* Let `t > T`. τ_t is `F_{t-1}`-measurable, so by (H↑_b) applied at the level `s = τ_t` (allowed because (H↑_b) holds for every s and the conditional expectation in Lemma 2.1(4) is taken with s fixed), `μ_t := E[X_t | F_{t-1}] ≤ (1 + b) m_F(τ_t)`. On G, `m_F(τ_t) ≤ U(τ_t)`, so `μ_t ≤ m↑_t` and `E[f↑_t | F_{t-1}] = 1 + λ (μ_t / m↑_t - 1) ≤ 1`. Likewise, on G with the decrease side on, `μ_t ≥ (1 - b) m_F(τ_t) ≥ m↓_t` and `E[f↓_t | F_{t-1}] = 1 + λ_t (1 - μ_t / m↓_t) ≤ 1` (λ_t is predictable). Both statements hold on `G ∩ {t > T}`, an event known at time t - 1.

*Step 3.* Steps 3 and 4 and parts 2 and 3 of the proof of Theorem 6.4 in `docs/math.md` use only "on an event known at T, every factor after T has conditional mean at most 1", with the bet grid fixed at T. They apply word for word with predictable bets: for `E_s = G ∩ {T = s}` and `s ≤ t ≤ s + m`, `S_t = 1_{E_s} (M_t + m - (t - s))` is a nonnegative supermartingale with `S_s = m 1_{E_s}`, Ville's inequality gives `P(E_s ∩ {alarm within m}) ≤ m P(E_s) / A`, and summing over s, `P(alarm within m) ≤ P(Gᶜ) + m P(G) / A ≤ δ + m / A`. ∎

**Remark (what the band buys).** Theorem 6.4 of `docs/math.md` needs one confidence sequence for one level fixed before the sequence starts, hence the sample split. Theorem A.1 needs a band for all levels that τ_t may take, and then τ_t may depend on anything in the past, including the monitored returns. With the grid band of 3.5 the levels must still be fixed before the band's data start (a short pilot fixes the grid), but the level used at each step is free.

### 3.4 Band 1: the time-uniform DKW band, integrated

**Lemma A.2 (from a distribution band to a band on m).** Let `V_i = r_i²`, `F_n` the empirical distribution function of `V_1, ..., V_n` and `S_n = 1 - F_n`. If `sup_x |F_n(x) - P(V ≤ x)| ≤ ε`, then for every s > 0, `L(s) ≤ m_F(s) ≤ U(s)` with

`U(s) = ∫_0^{s²} min(1, S_n(x) + ε) dx` and `L(s) = ∫_0^{s²} max(0, S_n(x) - ε) dx`.

*Proof.* `S_F(x) ∈ [S_n(x) - ε, S_n(x) + ε] ∩ [0, 1]` for every x; integrate over `[0, s²]` and use Lemma 2.1(1). ∎

With sorted values `v_(1) ≤ ... ≤ v_(n)`, `m_n(a) = ∫_0^a S_n = mean(min(V, a))`, `J = ⌊n ε⌋` and `J' = ⌈n (1 - ε)⌉ - 1`: `U(√a) = a` if `a ≤ v_(J+1)`, else `v_(J+1) + m_n(a) - m_n(v_(J+1)) + ε (a - v_(J+1))`; `L(√a) = m_n(a') - ε a'` with `a' = min(a, v_(J'+1))` (`truncated_means_from_band`; checked against numerical integration in the tests).

**Proposition A.3 (Howard and Ramdas 2022, Theorem 2).** Let `V_1, V_2, ...` be i.i.d. with any distribution function F, `t_min ≥ 1`, `A > 1/√2` and C > 0. With `γ² = (2 / η) (A - √(2 (η - 1) / C))²`,

`P(there is t ≥ t_min with sup_x |F_t(x) - F(x)| > A √((log(1 + log(t / t_min)) + C) / t)) ≤ inf_η 4 e^{-γ² C} (1 + 1 / ((γ² - 1) log η))`,

the infimum over η > 1 with γ² > 1. (`log(1 + log(t / t_min)) = log log(e t / t_min)`.) For a target α, `running.hr_dkw_constant` solves for C exactly as the authors' library does; their closed-form choice is A = 0.85 and C ≈ 0.8 log(1612 / α).

*How the statement was checked.* The paper's full text could not be fetched from this environment (arxiv.org and projecteuclid.org are blocked). The statement above was taken from the authors' own implementation (`EmpiricalProcessLILBound` in `src/confseq/uniform_boundaries.h` of github.com/gostevehoward/confseq, whose README says this function "is based on Theorem 2"), and the closed form `log log(et) + 0.8 log(1612/α)` was confirmed in search-engine extracts of the published paper. A reader with access should compare with Theorem 2 of the published version. Our test `test_hr_band_is_time_uniform_on_uniform_streams` checks the crossing probability by simulation (it is a check, not a proof).

**Atoms (price ticks).** If F has atoms, write `V_i = F⁻(W_i)` with `W_i` i.i.d. uniform and `F⁻` the generalised inverse. Then `V_i ≤ x` if and only if `W_i ≤ F(x)`, so `F_t(x) = G_t(F(x))` with `G_t` the empirical distribution function of the W's, and `sup_x |F_t(x) - F(x)| ≤ sup_u |G_t(u) - u|` for every t on the same probability space. The bound for continuous laws therefore implies the bound for every law.

**Width.** The DKW band has the same width ε at every x, so `U(s) - m_F(s)` can be close to `ε s²`. At s = 3 standard deviations of a normal return, `s² ≈ 9 m(s)`, so the relative width is about `9 ε`. This is why the integrated DKW band is much wider than a confidence sequence for one level (measured in E14, `results/e14w_bands.csv`).

### 3.5 Band 2: confidence sequences on a grid of levels

Fix a ratio ρ > 1 and integers `k_min < 0 < k_max`. After a pilot of `n_p` returns, put `s_c = c × scaled MAD(r_1, ..., r_{n_p})` and levels `s_k = s_c ρ^k`, `k = k_min, ..., k_max` (K levels). For each level run the upper and the lower empirical Bernstein confidence sequences of `docs/math.md` 6.2 on `Z_i = min(r_i², s_k²) / s_k²`, `i > n_p`, each at level `δ / (2K)`.

**Proposition A.4a (grid band).** Under (C), with `U_k(t)`, `L_k(t)` the bounds after the observations up to t (in units of r²), `P(there are k and t with m_F(s_k) ∉ [L_k(t), U_k(t)]) ≤ δ`.

*Proof.* The levels are fixed at time `n_p`, before the sequences start. Under (C) the conditional mean of `Z_i` given the past is the constant `m_F(s_k) / s_k²` for every `i > n_p`, so Corollary 6.3 of `docs/math.md` gives failure probability at most `δ / (2K)` for each of the 2K one-sided sequences; add them up. ∎

The monitor snaps τ_t to the grid (the level nearest in log scale to `c ×` the running scaled MAD), so the band is needed at the grid levels only. Remark: (2.1) extends the grid band to every s in `[s_{k_min}, s_{k_max}]` if a user prefers not to snap.

**Defaults (`running.py`).** ρ = 2^{1/4}, k from -6 to 6 (13 levels from 0.35 to 2.83 times the centre), `n_p = 50`, `c = 3`, bets of the sequences tuned for the window length (`n_plan`), cap 0.75. Running scale: scaled MAD of the last 120 finite returns, refreshed every 10 steps. These were fixed before the experiments (section 7).

**Why the grid band beats the DKW band.** Each grid sequence adapts to the variance of `min(r², s_k²)`, while the DKW band pays its full width at every x. Two costs remain: the union bound (`log(2K/δ)` instead of `log(2/δ)`, with K = 13 and δ = 0.01 about 7.9 instead of 5.3) and, at high levels, the range: an upper bound for a [0, 1] variable whose values are mostly tiny cannot bet more than about 0.75 against the possibility of a value near 1, so the upper bound at the top levels is loose (E14, `results/e14w_bands.csv`).

### 3.6 No blind period: a learning band and staggered chains

The standard protocol recalibrates for `cal_len` steps after each alarm, and the monitor is blind meanwhile. Here monitoring starts immediately after an alarm.

**Monitoring rule (awake).** A segment starts at time S (a stopping time, for example right after an alarm). The grid centre `s_c` is computed from returns before S (`c ×` the scaled MAD of the last 120 finite returns), so the levels are known at S. The grid confidence sequences run on the segment's returns `r_{S+1}, r_{S+2}, ...` (for the first `n_learn` of them, then they stop). Snapshot sizes `0 = e_0 < e_1 < ... < e_last = n_learn` are fixed in advance (`running.awake_epochs`: 0, 16, 23, 32, 45, ..., 1024, 1440). The band after e observations is `B(e) = (U(e), L(e))`, with `B(0)` the trivial band `(s_k², 0)`. The chain started at segment step j (it bets from the j-th segment return on) uses the snapshot `B(e(j))`, `e(j) = max{e_i ≤ j - 1}` for j ≥ 1 (all its band's data come before j). All chains use the same predictable τ_t. Chains with the same snapshot share one Shiryaev-Roberts recursion, so the statistic is `M_t = Σ_i M^{(i)}_t` with `M^{(i)}_t = f^{(i)}_t (M^{(i)}_{t-1} + 1{e(t - S) = e_i})`, where `f^{(i)}_t` uses `U(e_i)(τ_t)` or `L(e_i)(τ_t)`.

**Theorem A.4.** Suppose the segment's returns `r_{S+1}, r_{S+2}, ...` are i.i.d. with some law F given `F_S` (no change), and the grid band has total level δ. Then P(an alarm within m steps of S) ≤ δ + m / A, and the expected number of steps to the first alarm is at least (1 - δ) A.

*Proof.* Work given `F_S` and number the segment's steps 1, 2, .... Let `N` be the first segment step at which the band fails at some grid level (`N = ∞` if never). Proposition A.4a (applied to the segment's returns, which are i.i.d. F; the levels are fixed at S) gives `P(N < ∞) ≤ δ`, and N is a stopping time because each failure event is determined by the data so far.

The wealth of the chain started at j is `W^{(j)}_t = Π_{i=j}^{t} f^{(j)}_i`, with `f^{(j)}_i` using `B(e(j))`. Define the stopped chains `W'^{(j)}_t = 1{N > e(j)} W^{(j)}_t`. The event `{N > e(j)}` is known at segment step `e(j) ≤ j - 1`. On it, the snapshot `B(e(j))` covers `m_F` at every grid level, and for every `t ≥ j`, `E[X_t | F_{t-1}] = m_F(τ_t)` because `r_t` is independent of `F_{t-1}` with law F and τ_t is predictable; so, as in Step 2 of Theorem A.1, `E[f^{(j)}_t | F_{t-1}] ≤ 1` on `{N > e(j)}`. Hence `W'^{(j)}` is a nonnegative supermartingale from step j on with `W'^{(j)}_{j-1} = 1{N > e(j)} ≤ 1`, and `M'_t = Σ_{j ≤ t} W'^{(j)}_t` satisfies `E[M'_t | F_{t-1}] ≤ M'_{t-1} + 1`. Theorem B of `docs/math.md` uses only this inequality, so `P(max_{t ≤ m} M'_t ≥ A) ≤ m / A`. On `{N = ∞}`, `M' = M`. Hence `P(M reaches A within m) ≤ P(N < ∞) + m / A ≤ δ + m / A`, and the two-sided rule (each side at 2A) is handled as in section 1.

Run length. Let τ be the first step with `M_τ ≥ A`. For every n, `M'_t - t` stopped at the bounded time `τ ∧ n` gives `E[τ ∧ n] ≥ E[M'_{τ∧n}] ≥ A P(τ ≤ n, M'_τ = M_τ)`. Every chain started at `j ≤ τ` has `e(j) ≤ τ - 1`, so `M'_τ = M_τ` on `{N ≥ τ}`, and `P(τ ≤ n, N ≥ τ) ≥ P(τ ≤ n) - P(N < ∞)`. Letting n grow, `E[τ] ≥ A (P(τ < ∞) - δ)`; if `P(τ = ∞) > 0` then `E[τ] = ∞`. Either way `E[τ] ≥ (1 - δ) A`. ∎

Two points the proof makes precise. First, the good event "the band never fails" is not known at any finite time, so the proof does not condition on it; it multiplies each chain by the indicator of an event that IS known when the chain starts, `{N > e(j)}`, which keeps every stopped chain a supermartingale. Second, a chain never uses data from after its own start, so post-change returns never enter the band of a chain that starts before the change; chains that start after the change use bands that contain some post-change returns, but the chain that starts right after the change does not.

**Remark (what the learning band cannot do).** Validity is free, power is not. A chain that starts before a change at segment step ν bets against a band from fewer than ν returns, which is wide when ν is small; a chain that starts after the change bets against a band that already contains post-change returns, and once most of the band's data are post-change it tests the new level against itself. So a change early in a segment is hard to detect for route A, whatever the bets. E14 Q14f measures it: nothing is detected when ν = 30, and only large decreases when ν = 120. The rank detector of route B has the same dilution in its reference set (4.3) but its p-values are exact from the first return, and its early reference values still separate the two regimes; it detects most of these changes.

**Corollary A.5 (restart after every alarm).** Run segments back to back: segment k starts right after alarm k - 1 (standard protocol: after `cal_len` calibration steps; awake protocol: at once) and has its own band with level δ. Suppose that in every segment the null of Theorem A.1 (standard) or Theorem A.4 (awake) holds relative to that segment's own law. Let `N_m` be the number of alarms within the first m monitored steps. Then `E[N_m] ≤ (δ + m / A) / (1 - δ)`.

*Proof.* This is Corollary 6.6 of `docs/math.md` with "good segment" redefined. Call alarm k *good* if its segment's band has not failed before the alarm (for the frozen band: the band at T holds; for the awake band: `N_k ≥` the alarm's segment step). Bad alarms: a segment raises at most one alarm, and its band fails with probability at most δ given the information at its start, so `E[#bad] ≤ δ E[#segments started within the first m monitored steps] ≤ δ (1 + E[N_m])`. Good alarms: the process `W_t` = (sum of the stopped statistics M' at past good alarms) + (the current stopped statistic) - (number of monitored steps) is a supermartingale with `W_0 = 0`; each good alarm banks at least A because `M' = M` up to that alarm; optional stopping at the m-th monitored step gives `E[#good] ≤ m / A`. Add and rearrange. ∎

### 3.7 Why a plug-in baseline cannot carry a bound of this form

The default monitor bets against the plug-in mean. E8 found its one-day false alarm rate close to, but not above, the bound m / A. The next result shows that this is an accident of the horizon: at a fixed calibration length, no bound of the form `C m / A` can hold for the plug-in as A grows.

*Set-up.* Calibration window `W = (r_1, ..., r_n)`, `τ̂ = c ×` scaled MAD of W, `m̂ = n⁻¹ Σ min(r_i², τ̂²)`; the monitor uses `τ_t = τ̂` and `m↑_t = m↓_t = m̂` (band 0). Assume that all returns, calibration and monitored, are i.i.d. with a law F (the null holds exactly). Write `X = min(r², τ̂²)` for a monitored return, `B̂ = τ̂² / m̂`, and, given W,

`g↑(W) = max_{λ ∈ Λ↑} E_F[log(1 + λ (X / m̂ - 1)) | W]`, `g↓(W) = max_φ E_F[log(1 + λ_φ (1 - X / m̂)) | W]` with `λ_φ = φ / (B̂ - 1)` (`g↓ = -∞` if `B̂ ≤ 1`),

`g(W) = max(g↑(W), g↓(W))`, K the number of bets on the side that attains it (6 or 5) and `f_max(W)` the largest factor of the attaining bet (`1 + λ (B̂ - 1)` or `1 + λ_φ`).

**Proposition A.6 (the plug-in drifts).** Let `T_A` be the first alarm of the two-sided plug-in rule (each side at 2A). If `g(W) > 0`, then

`E[T_A | W] ≤ (log(2 A K) + log f_max(W)) / g(W)`.

Hence for every η > 0 and κ < ∞ with `p_{η,κ} = P(g(W) ≥ η, f_max(W) ≤ e^κ)` and every m ≥ 1, `P(T_A ≤ m) ≥ p_{η,κ} (1 - (log(12 A) + κ) / (η m))`. The probability `p_{η,κ}` depends on n and F but not on A, so with `m_A = (log A)²`, `P(T_A ≤ m_A) → p_{η,κ}` as A → ∞ while `m_A / A → 0`: if `p_{η,κ} > 0` for some η, κ, then no inequality `P(T_A ≤ m) ≤ C m / A` with a constant C holds for all A and m.

*Proof.* Fix W. Given W, the monitored returns are i.i.d. F and independent of W, so the factors `f_i` of the attaining bet are i.i.d. with `E log f_i = g(W) > 0`. They are bounded: `f_i ≤ f_max(W)` and `f_i ≥ 1 - λ ≥ 1/2` (increase side) or `f_i ≥ 1 - λ_φ (B̂ - 1) = 1 - φ ≥ 0.05` (decrease side). The Shiryaev-Roberts statistic of this bet is `Σ_{j ≤ t} Π_{i=j}^{t} f_i ≥ W_t := Π_{i=1}^{t} f_i`, and the side's statistic, an average over K bets, is at least `W_t / K`. So `T_A ≤ T' = inf{t : S_t ≥ b}` with `S_t = log W_t` and `b = log(2 A K)`. The process `S_t - g t` is a martingale with bounded increments; optional stopping at the bounded time `T' ∧ k` gives `g E[T' ∧ k | W] = E[S_{T' ∧ k} | W] ≤ b + log f_max(W)`, because S is below b before T' and the last step adds at most `log f_max`. Let k → ∞ (monotone convergence). The second claim follows from Markov's inequality `P(T' > m | W) ≤ E[T' | W] / m` on the event `{g ≥ η, f_max ≤ e^κ}` (K ≤ 6 gives `log(2AK) ≤ log(12A)`), integrated over W. The third is immediate. ∎

*When is `p_{η,κ} > 0`?* `g↑(W) > 0` as soon as the true mean `m_F(τ̂)` exceeds `m̂` by more than about `λ_min Var_F(X / m̂) / 2` in relative terms (λ_min = 0.02, the smallest bet), and symmetrically for `g↓` when `m̂` is too large. The relative error of `m̂` has a standard deviation of order `n^{-1/2}`, so at n = 390 or 1440 a large share of calibration windows gives one side a positive drift. E14 computes this share, the median run length bound and the Markov lower bound on the false alarm rate for A = 7800 and A = 10⁶ (`results/e14w_plugin_drift.csv`, expectations over F by quadrature on 4 million draws), and Q14a shows the consequence at A = 7800: within 1950 steps (bound 25%) the plug-in alarms on 41% to 43% of null paths at n = 390.

*Why the other monitors escape.* With a band, on the event that the band holds, every bet has `E[f_t | F_{t-1}] ≤ 1`, hence `E log f_t ≤ 0`: no chain drifts up, which is what Theorems A.1 and A.4 and `docs/math.md` Theorem 6.4 use. The rank detector of route B has `E[f(p_t) | G_{t-1}] = 1` exactly. The frozen rank reference (`B-freeze`, section 4.3) has the same defect as the plug-in, with the frozen reference in the role of W.

## 4. Route B: ranks

### 4.1 Exact p-values

Let `a_t = |r_t|` (any fixed nonconformity score works; ranks of |r| and of r² are the same) and let `θ_1, θ_2, ...` be i.i.d. uniform on [0, 1], independent of the returns. The *smoothed conformal p-value* of step t is

`p_t = ( #{i ≤ t : a_i > a_t} + θ_t #{i ≤ t : a_i = a_t} ) / t`

(the count of equal values includes i = t). A small p_t means `|r_t|` is large compared with all earlier values; with no ties `p_t = (#{i < t: a_i > a_t} + θ_t) / t`.

**Proposition B.1 (Vovk, Nouretdinov and Gammerman 2003; Vovk, Gammerman and Shafer 2005).** If `r_1, r_2, ...` is exchangeable, then `p_1, p_2, ...` are independent and uniform on [0, 1].

*Proof (included because the whole route rests on it).* Fix n. For t ≤ n let `B_t` be the bag (multiset) of `a_1, ..., a_t` and `H_t = σ(B_t, a_{t+1}, ..., a_n, θ_{t+1}, ..., θ_n)`. Then (i) `H_1 ⊇ H_2 ⊇ ... ⊇ H_n` (from `B_t` and `a_{t+1}` one gets `B_{t+1}`); (ii) `p_{t+1}, ..., p_n` are `H_t`-measurable (each `p_s`, s > t, is a function of `B_s`, `a_s` and `θ_s`); (iii) given `H_t`, `a_t` is a uniformly random element of the bag `B_t` (by exchangeability, the order of the first t values given their bag and the later values is uniform over the orderings) and `θ_t` is independent uniform. If the bag has distinct values `v_1 > v_2 > ...` with multiplicities `k_1, k_2, ...` and `a_t = v_j` (probability `k_j / t`), then `p_t = (k_1 + ... + k_{j-1} + θ_t k_j) / t` is uniform on the j-th of consecutive intervals of lengths `k_j / t` that tile [0, 1]. So `p_t` is uniform given `H_t`, hence independent of `H_t` and of `p_{t+1}, ..., p_n`. Then `P(p_1 ≤ u_1, ..., p_n ≤ u_n) = E[1{p_2 ≤ u_2, ...} P(p_1 ≤ u_1 | H_1)] = u_1 P(p_2 ≤ u_2, ..., p_n ≤ u_n)`, and induction gives `Π u_t`. ∎

Two consequences matter here. Ties (price ticks, zero returns) are handled exactly by the randomisation. And the p-values are independent of each other, **but not of the raw past**: given `r_1, ..., r_{t-1}`, the law of `p_t` is the law of the rank of a new draw against a fixed sample, which is not uniform. The martingale property therefore holds in the filtration `G_t = σ(p_1, ..., p_t)`, and bets must be chosen from past p-values only (here they are fixed functions).

### 4.2 The exact detector

A *betting function* is a measurable `f : [0, 1] → [0, ∞)` with `∫_0^1 f(u) du = 1`. Choose betting functions `f_t` predictably in G (fixed functions are the simplest case).

**Theorem B.2.** If the returns are exchangeable, then for the statistic `M_t = f_t(p_t) (M_{t-1} + 1)` started at any fixed time W (a warm-up whose returns only enter the reference set), `M_t - (t - W)` is a G-martingale; hence P(an alarm within m steps of W) ≤ m / A and the expected time to the first alarm is at least A, with no estimated quantity and no failure probability. Averages over several betting functions, and the two-sided rule (increase and decrease families, each at 2A), keep both bounds.

*Proof.* By Proposition B.1, `p_t` is uniform and independent of `G_{t-1}`, so `E[f_t(p_t) | G_{t-1}] = ∫ f_t = 1` and `E[M_t | G_{t-1}] = M_{t-1} + 1`. The alarm time is a G-stopping time, and Theorems A and B of `docs/math.md` apply in the filtration G. ∎

This is the conformal Shiryaev-Roberts procedure of Vovk (2021) and Vovk et al. (2021) applied to the score |r|; what is new here is only the choice of score and betting functions for volatility, the decrease side, and the comparison with route A.

**Betting functions used (`running.py`).** *Likelihood ratio in rank space* ("lr"): if the reference law is normal and the new return is ρ times a draw from it, `p_t` has density `g_ρ(u) = ρ⁻¹ exp(-q(u)² (ρ⁻² - 1) / 2)` with `q(u) = Φ⁻¹(1 - u/2)` (the |Z| quantile). We use `f(u) = 1 + λ (g_ρ(max(u, u_0)) / C_ρ - 1)` with `u_0 = P(|Z| > 3)` and `C_ρ = u_0 g_ρ(u_0) + 1 - 2 (1 - Φ(q(u_0)/ρ))` (so that `∫ f = 1`), ρ ∈ {1.25, 1.5, 2} with λ = 1 (increase) and ρ ∈ {0.8, 0.67, 0.5} with λ = 0.95 (decrease). Clipping at `u_0` is the rank version of truncating at 3 standard deviations: one extreme return multiplies the wealth by at most `g_ρ(u_0) / C_ρ` (`results/e14w_constants.csv`). Validity never uses normality; normality only decides which alternatives the bets are tuned for. *Indicators* ("ind"): `1 + λ (1{p ≤ α} / α - 1)` (|r_t| above the running (1 - α) quantile of the past) for increases and `1 + λ (1{p > 1 - β} / β - 1)` (below the running β quantile) for decreases, small grids of (α, λ) and (β, λ). *Null band* b: the payoff is divided by 1 + b, `f / (1 + b)`, which has mean `1 / (1 + b) < 1`; it is a tolerance for mild departures from exchangeability and never hurts validity.

**Remark B.3 (no blind period).** After an alarm the reference set restarts empty and betting starts at the next return (warm-up 0). The p-values are exact from the first one, so Theorem B.2 holds for the new segment if its returns are exchangeable. Early p-values are coarse (the first is just θ), which costs power, not validity. Corollary A.5 holds with δ = 0: `E[N_m] ≤ m / A`.

### 4.3 What is lost: dilution of the reference set

Let the reference set hold `n_0` pre-change values, and let k post-change values have arrived. For continuous laws and large `n_0`, the p-value of a new post-change return with `|r| = a` is close to

`π_k(a) = (n_0 S_F(a) + k S_G(a)) / (n_0 + k)`,

a mixture of the pre-change and the post-change tail probabilities. Under the post-change law G, the expected log-growth of a bet f at step k is `g_k = E_G[log f(π_k(a))]`. At k = 0 it is the full rate `E_G[log f(S_F(a))]` (for the likelihood-ratio bet tuned to the true ρ this is the Kullback-Leibler divergence, the best possible rate). As `k / n_0 → ∞`, `π_k(a) → S_G(a)`, which is uniform under G, and `g_k → E[log f(U)] ≤ log E[f(U)] = 0` (Jensen). The evidence per step therefore decays on the time scale `n_0`: changes that are detected within a small fraction of `n_0` steps lose little, slow detections lose a lot. `results/e14w_dilution.csv` gives `g_k` for a normal scale change by numerical integration.

**Freezing the reference set** (stop adding values once a statistic reaches a candidate level, here 20) removes the dilution but also the exactness: with a frozen reference R of size n, later p-values are, given R, i.i.d. with a law that depends on R; their conditional mean payoff `E[f(p) | R]` differs from 1 by a random amount of order `n^{-1/2}`, and `E_R[(E[f(p) | R])^m] ≥ 1` by Jensen, so the product of bets is not a supermartingale. This is the plug-in problem again, in rank space. For the indicator bet the conditional exceedance probability given R is `S_F` at an order statistic, whose law is known exactly (a Beta law, for any continuous F); betting against an upper confidence bound for it restores a proof, and that is route A with one level. `B-freeze` in E14 measures how much validity the heuristic loses.

## 5. Comparison

| | Route A, frozen band | Route A, awake | Route B | E8 (fixed τ) |
|---|---|---|---|---|
| Null | calibration i.i.d. F; monitored conditional truncated means ≤ (≥) those of F at every level | segment i.i.d. | segment exchangeable | conditional truncated mean ≤ (≥) calibration average at one level |
| Allows dependence in monitored data | yes (only conditional truncated means are restricted) | no | no | yes |
| τ during monitoring | any predictable rule (snapped to the grid) | same | not needed | fixed |
| Bound | δ + m / A | δ + m / A | m / A, exact | δ + m / A |
| Main cost | band width (union over levels; DKW much wider) | wide band early in a segment | dilution after a change; exchangeability fails under volatility clustering | sample split; fixed τ |
| Blind after an alarm | cal_len | 0 (120 at start-up only) | cal_len, or 0 with restart | cal_len |

## 6. Interpretation of an alarm

Route A alarms say that a truncated conditional second moment left the band of the calibration law at the current truncation level. Route B alarms say that the returns are not exchangeable: the most recent |r| values rank too high (increase family) or too low (decrease family) among the earlier ones. Under volatility clustering both nulls are false (`docs/math.md`, 6.5), so on real data an alarm is a statement about the data, not a "false alarm" in the sense of the theorems; E14 Q14c and Q14d measure what happens.

## 7. Choices fixed before the experiments

Grid ratio 2^{1/4} and 13 levels (centre ± 1.5 octaves) cover decreases to half and increases to about 2.8 times the calibration scale; pilot 50 returns; running scale from 120 returns refreshed every 10 steps (a predictable rule; for normal returns the relative standard deviation of a sample median of 120 absolute values is about 1.17 / √120, roughly 0.11, below the grid spacing 2^{1/4} - 1, roughly 0.19); δ = 0.01; bet grids of `detectors.py`; rank bets as in 4.2; freezing level 20; awake schedule `n_learn = 1440`, snapshots from 16 observations growing by √2, confidence sequence bets tuned for 300 observations (a compromise between early and late tightness chosen in a pilot on other seeds, `docs/iterations/TH2.md`).

## 8. Status

* **Proved here:** Lemma 2.1 and (2.1); Theorem A.1; Lemma A.2; Proposition A.4a; Theorem A.4; Corollary A.5 (given the arguments of Corollary 6.6 in `docs/math.md`); Proposition A.6 (the plug-in has no bound of the form C m / A); the reduction of Proposition A.3 from continuous laws to laws with atoms; Proposition B.1 (a known result; the proof is given for completeness); Theorem B.2 (given Proposition B.1 and Theorems A and B of `docs/math.md`); the dilution limit `g_∞ ≤ 0` of 4.3.
* **Cited:** Proposition A.3 (Howard and Ramdas 2022, Theorem 2; statement checked against the authors' code, see 3.4); the empirical Bernstein confidence sequences and Ville's inequality as in `docs/math.md` section 6; conformal test martingales and the conformal Shiryaev-Roberts and CUSUM procedures (Vovk 2021; Vovk et al. 2021); e-detectors (Shin, Ramdas and Rinaldo 2024).
* **Checked numerically:** coverage of both bands, uniformity and independence of the p-values with ties, the exact false alarm bound of route B, the δ + m / A bound of route A with running τ and of the awake monitor, streaming equals vectorised (`tests/test_running.py`); validity at scale, delays, seasonality and drift, real data and cost (E14); the share of calibration windows with positive plug-in drift and the run length bound of Proposition A.6 (E14 Q14w); false alarms per day of the restarted monitors against Corollary A.5 on 30-day streams (E14 Q14h).
* **Assumed, not proved:** the nulls on real data (volatility clustering violates all of them); the deseasonalising profile is treated as known; missing returns are skipped (valid if missingness says nothing about the size of the return); the dilution formula `π_k` is an approximation for large reference sets.

## 9. References (each checked by a web search or fetch on 2 November 2026)

* Howard, S. R. and Ramdas, A. (2022). Sequential estimation of quantiles with applications to A/B testing and best-arm identification. *Bernoulli* 28(3), 1704 to 1728. doi:10.3150/21-BEJ1388. https://projecteuclid.org/journals/bernoulli/volume-28/issue-3/Sequential-estimation-of-quantiles-with-applications-to-A-B-testing/10.3150/21-BEJ1388.full ; arXiv:1906.09712. Authors' code: https://github.com/gostevehoward/confseq (README and `src/confseq/uniform_boundaries.h`, fetched).
* Vovk, V. (2021). Testing randomness online. *Statistical Science* 36(4), 595 to 611. doi:10.1214/20-STS817. https://projecteuclid.org/journals/statistical-science/volume-36/issue-4/Testing-Randomness-Online/10.1214/20-STS817.full (conformal CUSUM and Shiryaev-Roberts procedures and their validity).
* Vovk, V., Petej, I., Nouretdinov, I., Ahlberg, E., Carlsson, L. and Gammerman, A. (2021). Retrain or not retrain: conformal test martingales for change-point detection. *Proceedings of Machine Learning Research* 152 (COPA 2021), 191 to 210. https://proceedings.mlr.press/v152/vovk21b.html ; arXiv:2102.10439.
* Vovk, V., Nouretdinov, I. and Gammerman, A. (2003). Testing exchangeability on-line. *Proceedings of the 20th International Conference on Machine Learning (ICML 2003)*, 768 to 775. https://mlanthology.org/icml/2003/vovk2003icml-testing (smoothed conformal p-values and exchangeability martingales).
* Vovk, V., Gammerman, A. and Shafer, G. (2005; 2nd edition 2022). *Algorithmic Learning in a Random World.* Springer. https://link.springer.com/book/10.1007/978-3-031-06649-8 (second edition).
* Fedorova, V., Gammerman, A., Nouretdinov, I. and Vovk, V. (2012). Plug-in martingales for testing exchangeability on-line. *ICML 2012*. arXiv:1204.3251. https://icml.cc/imls/conferences/2012/papers/808.pdf (describes the uniform and independent p-values under exchangeability).
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science* 2(2), 229 to 260. doi:10.51387/23-NEJSDS51. https://nejsds.nestat.org/journal/NEJSDS/article/59
* Waudby-Smith, I. and Ramdas, A. (2024), Ville (1939), Howard et al. (2020, 2021): as cited in `docs/math.md`, section 6.
