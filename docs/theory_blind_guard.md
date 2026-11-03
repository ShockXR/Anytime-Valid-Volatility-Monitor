# Blind segments of the confidence-sequence baseline, and a proved guard (E24)

E19 found that a segment of the monitor with the confidence-sequence baseline (`baseline="cs"`, `docs/math.md`, section 6) can be unable to raise an increase alarm for its whole life, and E23 found the same failure on fresh data: on BTC in 2026 the primary monitor P1 (`PeriodicVolatilityMonitor(baseline="cs", delta=0.01, A = 36,000, period=1440)`, band b = 0.25) could not alarm upward for 28% of its monitored time, one segment for 47 days, and three of the sixteen releases fell in such segments (`docs/report_parts/E23.md`). E19 proposed a guard (discard the window and calibrate again) and left it unproved, because the windows kept are then selected by their bounds. This note proves a version of that guard, shows by a construction that the selection does matter, and states what the alternatives (a time-out and a fallback to the rank monitor) guarantee. Code: `src/avmonitor/guard.py`. Tests: `tests/test_guard.py`. Experiment: `experiments/e24_blind_guard.py`, report part `docs/report_parts/E24.md`.

Labels: **[PROVED]** (a complete proof is given here, using only Theorem 6.2 of `docs/math.md`, Ville's inequality, optional stopping at bounded stopping times, and Proposition B.1 of `docs/theory_running_tau.md` for the rank fallback), **[NUMERICAL]**, **[SKETCH]** (an argument that is not written out in full).

## 0. Summary

1. **The blind condition is decided when calibration ends** [PROVED, Proposition 1]. A segment is blind when `τ² ≤ m↑ = (1 + b) Ū_T`; then every increase factor is at most 1, the increase side cannot alarm within 2A monitored steps of a (re)start, and with P1's daily restart (period m < 2A) it can never alarm. Both τ² and Ū_T are known at the end of calibration T, a stopping time, so the event is in `F_T`. More generally, with `ρ = τ² / m↑ ≤ 1 + η` ("near-blind") every increase factor is at most `1 + λ η`, the expected log growth of every increase bet is at most `log(1 + λ_max η)` under any law, and no increase alarm can come within `d_min(ρ)` steps of a restart, an exact floor attained by `X ≡ τ²` (Proposition 2). For P1 (A = 36,000, the default grid with λ_max = 0.5) the increase side cannot alarm within a day for any data when `ρ ≤ 1.0105`; `d_min` is 1503 steps at ρ = 1.01, 371 at 1.05, 91 at 1.25 and 30 at 2.
2. **A retry guard with error spending** [PROVED, Theorem 1]. If the window just calibrated is (near-)blind, discard it and calibrate again on the next returns, with error budget δ_k for the k-th window of the segment, `δ_1 + ... + δ_K ≤ δ`, at most K windows, then fall back. Then P(the bounds of the kept window fail) ≤ `Σ_k δ_k P(window k is tried) ≤ δ`, given everything before the segment. The proof uses only that the windows are disjoint and each confidence sequence starts at a stopping time.
3. **Choosing which window to keep does break per-window coverage** [PROVED, Proposition 3]. There is a law (adaptive, inside the tolerance relative to every window) under which a window is kept exactly when its upper bound fails: given that a window before the K-th is kept, its bound fails with probability 1, and with the same budget δ for every window the probability that the kept window's bound fails is `1 - Π_k (1 - p_k)`, about K times the single-window failure probability. So E19's worry is correct for the naive guard; spending δ over the windows is what fixes it, and the union bound cannot be improved in general.
4. **Theorem 6.4, Corollary 6.6 and Corollary 1 hold unchanged** [PROVED, Theorem 2], with δ the total budget, T the end of the kept window, and the tolerance (H↑), (H↓) relative to the kept window's confidence-sequence part. The tolerance is an almost-sure hypothesis about conditional means, so the selection does not enter the proofs; what changes is the null itself (no change relative to the kept window, not relative to the first window), and neither null implies the other (Remark 3). A sufficient condition that does not depend on the selection is a constant conditional law from the end of the first window's τ part on (Corollary 6.5 of `docs/math.md`), which covers a step in volatility inside that τ part.
5. **Alternatives** [PROVED, Propositions 4 and 5, Theorem 3]. A time-out after D periods without an alarm keeps every per-segment bound and costs `δ / D` alarms per monitored period in the long-run bound: `E[N_M] ≤ (δ (1 + ⌊M / (D m)⌋) + M / A) / (1 - δ)`. A fallback to E14's awake rank monitor after K blind windows gives the mixed guarantee `P(an alarm within m monitored steps of a segment) ≤ Σ_k δ_k + m / A`, where the cs branch needs the tolerance relative to the kept window and the rank branch needs the returns after the fallback time to be exchangeable given the past; the per-period and long-run bounds carry over, the latter with a per-segment proof because the two branches live in different filtrations.
6. **A calendar day is not a period of a segment** [PROVED upper bound, SKETCH of the lower bound; not specific to the guard, Remark 6]. Corollary 1 bounds a given period of a given segment. For a fixed calendar day d the proof gives `δ E[number of segments started by the end of day d] + m / A`, and an adaptive null can push the probability above `δ + m / A`.
7. **Recommended guard**: section 6.

## 1. Set-up

As in `docs/math.md`, 6.1 to 6.3. Returns `r_1, r_2, ...`, `F_t = σ(r_1, ..., r_t)` (with any independent randomisation the monitor uses added to `F_0`). A **segment** starts at a stopping time S (the start of the stream, or the step after an alarm). Its calibration consists of one or more **windows**. Window k starts at a stopping time `S_{k-1}` (`S_0 = S`). Its τ part is the next `n_τ` returns (longer if they are degenerate, `monitor.py`), which fix `τ_k` at a stopping time `σ_k`. Its confidence sequence (CS) part runs on `X_t = min(r_t², τ_k²)` for `t > σ_k` and ends at a stopping time `T_k > σ_k`, where the upper and lower bounds `Ū^{(k)}`, `L̄^{(k)}` are read; their error levels are `δ_k / 2` each. The next window, if any, starts at `S_k = T_k`. So the windows are disjoint, consecutive, and each starts at a stopping time.

The window that is used for monitoring is the **kept window** κ, with `T = T_κ`, `m↑ = (1 + b) Ū^{(κ)}_T` and `m↓ = (1 - b) L̄^{(κ)}_T`. The monitoring rule, the bets, the two-sided alarm and the optional restart at period boundaries are those of `docs/math.md` 6.3 and `docs/theory_per_period.md` section 3. `μ_t = E[X_t | F_{t-1}]` with the kept window's τ, and `ν̄↑`, `ν̄↓` are the bet-weighted calibration averages of Theorem 6.2. The **good event** of window k is `G_k = {ν̄↑_{T_k} < U^{(k)}_{T_k}} ∩ {ν̄↓_{T_k} > L^{(k)}_{T_k}}` (normalised units, as in Theorem 6.4), and `G = G_κ`.

The **tolerance relative to the kept window** (horizon h) is (H↑) and (H↓) of `docs/math.md` 6.3 with T = `T_κ`, the minimum and maximum taken over the CS part `σ_κ < i ≤ T_κ` of the kept window, and `μ` computed with `τ_κ`:

* (H↑) `sup_{T_κ < t ≤ T_κ + h} μ_t ≤ (1 + b) min_{σ_κ < i ≤ T_κ} μ_i`, and (H↓) `inf_{T_κ < t ≤ T_κ + h} μ_t ≥ (1 - b) max_{σ_κ < i ≤ T_κ} μ_i`, almost surely.

## 2. The blind condition

**Lemma 1 (a blind segment cannot alarm upward) [PROVED].** If `τ² ≤ m↑`, every increase factor `1 + λ (X_t / m↑ - 1)` is at most 1, and each increase bet's Shiryaev-Roberts statistic satisfies `R_t ≤ R_{t-1} + 1`, so `R_t ≤ t - t₀` for t monitored steps after a (re)start at `t₀`. The increase side alarms only when the average of its bets reaches 2A, so it cannot alarm within 2A - 1 monitored steps of a (re)start; with a restart every m < 2A steps (P1: m = 1440, 2A = 72,000) it can never alarm.

*Proof.* `X_t ≤ τ² ≤ m↑` gives `X_t / m↑ - 1 ≤ 0`, and λ ≥ 0. Then `R_t = f_t (R_{t-1} + 1) ≤ R_{t-1} + 1`, and `R_{t₀} = 0`. An average of numbers each at most `t - t₀` is at most `t - t₀`. ∎

The decrease side is not affected: in E23 both blind BTC segments ended with a down alarm.

**Proposition 1 (decided at a stopping time) [PROVED].** For a fixed `η ≥ 0`, the near-blind event `B = {τ² ≤ (1 + η) m↑}` of a window is in `F_T`, and `T_B = T` on B, `T_B = ∞` off B, is a stopping time. With η = 0 this is the blind event of Lemma 1, E19's criterion.

*Proof.* τ is `F_σ`-measurable with σ < T, and `Ū_T = τ² min(U_T, 1)` is a function of the CS observations up to T and of the predictable bets, hence `F_T`-measurable; b and η are constants. For every t, `{T_B = t} = {T = t} ∩ B ∈ F_t` because `B ∩ {T = t} ∈ F_t` for `B ∈ F_T`. ∎

Since τ, `m↑` and the bets stay fixed for the life of the segment, B is also the event that the segment is (near-)blind from start to end: blindness is never discovered later.

In normalised units the condition reads `min(U_T, 1) ≥ 1 / ((1 + η)(1 + b))`: with b = 0.25 and η = 0, the upper bound must reach 80% of τ². For i.i.d. normal returns `E[X] / τ² ≈ 0.11` (c = 3), so this needs most CS-part returns to exceed τ: the CS part must be far more volatile than the τ part. With normal returns a step in standard deviation by a factor k after the τ part gives `τ² / m↑` near 2.1 (k = 2), 1.2 (k = 4), 1.0 (k = 6) and 0.95 (k = 8) for P1's configuration [NUMERICAL, `results/e24_synthetic.csv`, column `ratio_first_window_median`].

**Proposition 2 (power of a near-blind increase side) [PROVED].** Let `ρ = τ² / m↑` and let λ ∈ [0, 1) be an increase bet.

(a) Every factor satisfies `1 + λ (X_t / m↑ - 1) ≤ 1 + λ (ρ - 1)`, with equality exactly when `X_t = τ²`.

(b) For bets `λ_1, ..., λ_K` averaged with equal weights, the largest possible value of the side average d monitored steps after a (re)start, over all data, is `a_d(ρ) = (1/K) Σ_k Σ_{i=1}^{d} (1 + λ_k (ρ - 1))^i`, attained by `X ≡ τ²`. So the increase side cannot alarm before `d_min(ρ) = min{d : a_d(ρ) ≥ 2A}` steps after a (re)start, and this floor is attained. With a restart every m steps, the increase side can never alarm if `d_min(ρ) > m`.

(c) For any law, `E[log f_t | F_{t-1}] ≤ log(1 + λ (μ_t / m↑ - 1)) ≤ log(1 + λ (ρ - 1))`, so each increase bet's expected log growth per step, the rate at which it gathers evidence after any change, is at most `log(1 + λ_max η)` when `ρ ≤ 1 + η`. An increase can be detected at all only if the post-change truncated mean exceeds `m↑`, and the truncated mean can never exceed τ²: the side sees at most a factor ρ above its threshold.

*Proof.* (a) `X_t ≤ τ² = ρ m↑`. (b) The recursion `R_t = f_t (R_{t-1} + 1)` is nondecreasing in `f_t` and in `R_{t-1}` (both nonnegative), so by induction every `R_t` is at most the value obtained with every factor at its maximum `g_k = 1 + λ_k (ρ - 1)`, which is `Σ_{i=1}^{d} g_k^i`; the data `X ≡ τ²` attain it for every bet at once. (c) Jensen's inequality for the concave logarithm, then (a) applied to `μ_t ≤ τ²`. ∎

For P1 (A = 36,000, the default grid 0.02 to 0.5): `d_min = 72,000` at ρ = 1, 9854 at 1.001, 2703 at 1.005, 1503 at 1.01, 827 at 1.02, 371 at 1.05, 202 at 1.1, 91 at 1.25, 51 at 1.5 and 30 at 2; P1's increase side is blind in effect (no alarm possible within a day, whatever the data) for `ρ ≤ 1.0105` [NUMERICAL, exact recursion, `results/e24_near_blind.csv`]. The floor of (b) is weak above ρ ≈ 1.05, because it is attained only when every return exceeds τ; the rate bound (c) is the useful one there: with ρ = 1.25 the largest bet gains at most `log 1.125 = 0.118` nats per step against `log 3 = 1.10` at a typical ρ = 5.

## 3. The retry guard

**The guard.** Fix in advance an integer `K ≥ 1`, budgets `δ_1, ..., δ_K > 0` with `Σ_k δ_k ≤ δ`, and `η ≥ 0`. For k = 1, 2, ...: calibrate window k with budget δ_k; if it is not near-blind (`τ_k² > (1 + η) m↑^{(k)}`), keep it; if it is near-blind and k < K, discard it and start window k + 1 at the next return; if k = K, fall back: either keep window K anyway ("keep"), or switch to the rank monitor for the segment ("rank", Theorem 3). After an alarm the next segment starts again with window 1. Default schedule (`guard.spending_schedule`, "halving"): `δ_k = δ 2^{-k}` for k < K and `δ_K = δ_{K-1}`, so the budgets sum to δ; K = 1 gives `δ_1 = δ`, the unguarded monitor.

**Lemma 2 (Theorem 6.2 at a stopping time) [PROVED].** Let σ be a stopping time, let τ be `F_σ`-measurable, and run the sequence of `docs/math.md` 6.2 on `Z_t = min(r_t², τ²) / τ²`, `t > σ`, with predictable bets. Then on `{σ < ∞}`, `P(ν̄↑_t ≥ U_t for some t > σ | F_σ) ≤ δ_u` and the same for the lower bound with `δ_l`.

*Proof.* Put `K_t = 1` for `t ≤ σ` and multiply, for `t > σ`, by the factor `exp(λ↑_t (ν_t - Z_t) - ψ↑_t e_t)` of the proof of Theorem 6.2. Since `{σ < t} ∈ F_{t-1}` and τ, the centre and the bet are `F_{t-1}`-measurable on that event, the computation of Theorem 6.2 gives `E[K_t / K_{t-1} | F_{t-1}] ≤ 1`, so K is a nonnegative supermartingale with `K_σ = 1`. For `A ∈ F_σ`, `(1_A K_{σ+j})_{j ≥ 0}` is a nonnegative supermartingale in `(F_{σ+j})` (optional sampling for nonnegative supermartingales), and Ville's inequality gives `P(A ∩ {sup_{t > σ} K_t ≥ 1/δ_u}) ≤ δ_u P(A ∩ {σ < ∞})`. The event in braces is the event that the bound fails at some t > σ (proof of Theorem 6.2). ∎

**Theorem 1 (the kept window's bounds) [PROVED].** For every segment with start S, on `{S < ∞}`,

`P(Gᶜ | F_S) ≤ Σ_{k=1}^{K} δ_k P(window k is tried | F_S) ≤ Σ_k δ_k ≤ δ`,

where `Gᶜ` is the event that a window is kept and its bounds fail (with the "keep" fallback a window is always kept; with the "rank" fallback none is kept when all K windows are blind). The discard rule may be any rule decided at `T_k` from the data up to `T_k`; near-blindness is one.

*Proof.* Let `W_k = {window k is tried} = {the windows 1, ..., k - 1 were all discarded}`, an event in `F_{T_{k-1}} ⊆ F_{σ_k}` (empty for k > K). If window k is kept and its bounds fail, then `W_k ∩ G_kᶜ` occurred, and `G_kᶜ` is contained in the event that window k's sequence fails at some t > σ_k (it is read at `T_k > σ_k`). By Lemma 2 applied at `σ_k` with `A = W_k ∩ A'`, `A' ∈ F_S`, `P(A' ∩ W_k ∩ G_kᶜ) ≤ (δ_k / 2 + δ_k / 2) P(A' ∩ W_k)`. A union over k gives `P(A' ∩ Gᶜ) ≤ Σ_k δ_k P(A' ∩ W_k)` for every `A' ∈ F_S`. ∎

What the proof uses: the windows do not overlap, so window k's sequence is a fresh sequence started at its own stopping time; whether window k is tried is known before its sequence starts; and each window has its own share of δ. It does not use, and does not give, coverage of the kept window given that it was kept.

**Proposition 3 (the selection is real; the union bound is sharp) [PROVED].** Fix c, b, η, K and the window length (τ part `n_τ`, CS part `n_cs`), and put `q = 1 / ((1 + η)(1 + b)) < 1` (any b > 0 or η > 0). Consider the adaptive law: in the τ part of every window the returns are `±s` with fair random signs (s > 0), so `τ = c × 1.4826 s` exactly; every other return (CS parts and monitoring) is `±τ` with probability q and 0 otherwise, independently of everything before. Then:

(a) the tolerance (H↑), (H↓) holds relative to every window, for every horizon and every b ≥ 0 (μ_t = q τ² in every CS part and during monitoring);
(b) a window is near-blind exactly when `min(U_T, 1) ≥ q`, and its upper bound fails exactly when `U_T ≤ q = ν̄↑_T`; so a window that is not near-blind has a failing upper bound, and `P(Gᶜ | window k kept) = 1` for every k < K with positive probability of being kept;
(c) with budgets `δ_1, ..., δ_K`, `P(Gᶜ) ≥ 1 - Π_{k=1}^{K} (1 - p(δ_k))`, where `p(δ) = P(U_T < q)` for one window's upper sequence at level δ / 2 on i.i.d. Bernoulli(q) data. With the same budget δ for every window (no spending) this is about `K p(δ)`, K times the single-window failure probability, and the bound `P(Gᶜ) ≤ Σ_k δ_k` of Theorem 1 cannot be improved beyond the gap between `p(δ)` and `δ / 2`.

*Proof.* (a) `μ_t = E[min(r_t², τ²) | F_{t-1}] = q τ²` outside the τ parts, which the tolerance excludes. (b) In normalised units `ν_i = q` for every CS observation, so `ν̄↑_T = q` whatever the bets; near-blind means `τ² ≤ (1 + η)(1 + b) τ² min(U_T, 1)`, that is `min(U_T, 1) ≥ q`; and `G↑ᶜ = {U_T ≤ ν̄↑_T} = {U_T ≤ q}`. A window k < K that is kept is not near-blind, so `U_T < q` and its upper bound fails. (c) The CS data of the windows are i.i.d. Bernoulli(q) sequences (times τ²), independent across windows, so the events `{U^{(k)}_{T_k} < q}` are independent with probabilities `p(δ_k)`. If one of them occurs, let k be the first: windows 1 to k - 1 have `U ≥ q` and are discarded, window k is not near-blind and is kept, and its upper bound fails. Hence `Gᶜ` contains the union of these events, whose probability is `1 - Π_k (1 - p(δ_k))`. ∎

So E19's worry is exactly right for the naive guard: conditioning on being kept destroys the per-window coverage, and only the unconditional union over the windows tried survives. The numbers for a small configuration are in `results/e24_selection.csv` [NUMERICAL]. The construction makes the coverage event fail; whether that produces false alarms depends on how far the bound undershoots (the increase side drifts upward only if `(1 + b) U_T < q`), so the false alarm rate itself need not rise by the same factor. Theorem 6.4 charges the coverage event, which is what the proof needs.

**Theorem 2 (the guarantees hold unchanged) [PROVED].** Run the guard with any K, budgets with `Σ_k δ_k ≤ δ` and the "keep" fallback, and let T be the end of the kept window.

1. (Theorem 6.4 of `docs/math.md`.) Parts 1 to 4 hold as stated, with δ the total budget and (H↑), (H↓) relative to the kept window. In particular, under the tolerance with h = m, `P(an alarm within the first m monitored steps of the segment) ≤ δ + m / A`.
2. (Corollary 6.6.) If the tolerance holds with h = ∞ in every segment relative to that segment's kept window, then `E[N_M] ≤ (δ + M / A) / (1 - δ)` for the number `N_M` of alarms in the first M monitored steps. Calibration steps, including those of discarded windows, are not monitored steps.
3. (Corollary 1 of `docs/theory_per_period.md`.) With a restart at every period boundary (carry at most κ), for every period `(b, b']` of the segment with `T ≤ b` and `b' - b ≤ m`, if the tolerance holds relative to the kept window for `T < t ≤ b'`, then `P(the segment raises an alarm in (b, b']) ≤ δ + (m + κ) / A`, marginally.

*Proof.* Each of the three proofs uses the calibration only through: (i) T is a stopping time; (ii) τ, `m↑`, `m↓` and the bets are `F_T`-measurable; (iii) G is in `F_T` and `P(Gᶜ) ≤ δ` (given `F_S` for Corollary 6.6); (iv) on G, the tolerance gives `μ_t ≤ (1 + b) min_i μ_i ≤ (1 + b) μ̄↑_T ≤ m↑` and the mirror inequality for the decrease side (Step 2 of Theorem 6.4). Under the guard, (i) `T = T_κ` with `{κ = k} = W_k ∩ B_kᶜ ∈ F_{T_k}` (B_k the near-blind event of window k; with the "keep" fallback `{κ = K} = W_K`), so `{T = t} = ∪_k {κ = k, T_k = t} ∈ F_t`; (ii) on `{κ = k}` these are window k's quantities, each `F_{T_k}`-measurable; (iii) `G = ∪_k ({κ = k} ∩ G_k) ∈ F_T`, and Theorem 1 gives the bound; (iv) on `{κ = k} ∩ G_k` the inequalities are those of Step 2 for window k, with the tolerance relative to window k. Steps 3 and 4, parts 2 to 4 of Theorem 6.4, the good-segment supermartingale of Corollary 6.6 and the supermartingale of Corollary 1 use nothing else. In Corollary 6.6 the bad-segment term uses `P(Gᶜ_j | F_{S_j}) ≤ δ` for segment j, which is Theorem 1 conditionally on the segment's start. ∎

**Remark 1 (cost of spending).** The first window gets `δ_1 < δ` (δ / 2 with halving), so its bounds are a little wider than the unguarded ones: the half-width scales roughly like `√log(2 / δ_1)`, a factor `√(log 400 / log 200) ≈ 1.06` for δ = 0.01. The total δ, and so the choice `A = m / (α - δ)`, is unchanged.

**Remark 2 (a retry without spending is a time-out with D = 0).** If each window kept its own δ, each window would be a segment of its own in the sense of Corollary 6.6; the long-run bound would charge δ per window tried, `E[N_bad] ≤ δ E[number of windows tried among the first M monitored steps]`, which no bound on the number of retries controls unless K is capped, and the per-segment bound would be `K δ + m / A` (Proposition 3 shows this order is attained for the coverage event).

**Remark 3 (the tolerance relative to the kept window).** (H↑), (H↓) are almost-sure conditions on conditional means, which are properties of the law given the past; the proofs need them on every path, for the window that is kept on that path, and no conditional statement given the selection is used. So the selection does not break the hypothesis as a hypothesis. It does change what the hypothesis says:

* A sufficient condition that does not depend on which window is kept is the tolerance relative to every window that may be tried. Under a constant conditional law from `σ_1` on (`docs/math.md`, Corollary 6.5; for example i.i.d. returns with any law), it holds for every window with b = 0, because `σ_k ≥ σ_1`. This covers the motivating case, a step in volatility inside the first window's τ part (design `step` of E24), for both the guarded and the unguarded monitor.
* The guarded and unguarded monitors test different nulls, and neither implies the other. If volatility rises in the middle of window 1's CS part and stays, the tolerance fails relative to window 1 but holds relative to window 2 (the guarded monitor is covered, the unguarded is not, though if window 1 is blind it cannot alarm upward anyway). If volatility falls during window 2 and comes back, the tolerance holds relative to window 1 but fails relative to window 2, and the guarded monitor may alarm, correctly relative to its own reference, where the unguarded would not.
* In part 4 of Theorem 6.4 the event `V_m` (the tolerance fails within m steps) is relative to the kept window, and the selection can change its probability; the bound `δ + P(V_m) + m / A` holds as stated.

**Remark 4 (decrease side).** The guard looks only at the increase side. Any other rule decided at `T_k`, for example "discard if the decrease side is off", is covered by the same Theorem 1.

**Remark 5 (missing returns, degenerate windows).** Both are handled inside a window as in `monitor.py` (a window is extended until usable; the τ part may restart before any sequence has started, which spends no δ). Each window still has exactly one sequence per side, started at a stopping time.

**Remark 6 (a fixed calendar day; not specific to the guard).** Corollary 1 (and part 3 of Theorem 2) bounds a given period of a given segment. For a fixed period `(b, b']` of the stream clock (for example a calendar day), with calibration at least m steps so that at most one segment monitors in it, the same argument gives

`P(an alarm in (b, b']) ≤ δ E[#{segments whose calibration starts before b'}] + (m + κ) / A` [PROVED].

*Proof.* Let J be the segment that monitors in the period. On `G_J`, at the first monitored step of the period its statistic is at most κ (restart at b, or a fresh start), and Ville's inequality from that stopping time gives `(m + κ) / A`. The event `G_Jᶜ` is contained in the union over segments j started before b' of `G_jᶜ`, and Theorem 1 at each start gives `δ` times the expected number of such segments. ∎

The first term grows with the number of earlier segments. [SKETCH] An adaptive null can make it bite: the adversary sees the data, so it knows whether the current segment's upper bound undershoots its calibration mean; in such a segment it holds the conditional mean at the lower edge allowed by the tolerance until the chosen day, then raises it to the upper edge, which exceeds `m↑`; in the other segments it plays E16's bold play to end them quickly with alarms. The probability that the segment active on that day is a bad one then grows with the number of segments before it. So "δ + m / A for every day" is a statement about the periods of each segment, not about calendar days; the long-run statement for calendar time is Corollary 6.6. We have not simulated this adversary.

## 4. Alternatives

### 4.1 A time-out

**Rule.** A segment that has monitored D full periods (D m monitored steps) without an alarm ends, and a new segment calibrates (each segment with budget δ).

**Proposition 4 [PROVED].** (a) Theorem 6.4 and Corollary 1 hold for each segment as stated (a time-out only shortens a segment). (b) If the tolerance holds with h = ∞ relative to every segment's window, `E[N_M] ≤ (δ (1 + ⌊M / (D m)⌋) + M / A) / (1 - δ)`. Per monitored period, the long-run bound becomes `(δ / D + m / A) / (1 - δ)` as M grows: the time-out costs δ / D alarms per period.

*Proof.* (a) The time-out is a stopping time, and Ville's inequality for the first m steps or for a period does not care whether the segment would have continued. (b) In the proof of Corollary 6.6 the bad-segment term is `δ E[#{k : C_k}]` with `C_k` = "segment k starts before M monitored steps have passed". A segment starts after an alarm or after a time-out, and a time-out among the first M monitored steps needs D m monitored steps of its own, so `#{k : C_k} ≤ 1 + N_M + ⌊M / (D m)⌋`. In the good-segment supermartingale W, a time-out drops the current statistic without banking it, which only lowers W. So `E[N_M] ≤ δ (1 + E[N_M] + ⌊M / (D m)⌋) + M / A`. ∎

With δ = 0.01 and D = 7 the extra term is 0.0014 per day against `m / A = 0.04`. The costs are elsewhere: a blind segment still lasts up to D periods, and each time-out spends one calibration window (for P1, one day in D + 1 if no alarm comes). Applied only to segments that are near-blind at T (a decision in `F_T`), the time-out is a delayed retry without spending, with the accounting of Remark 2.

### 4.2 A fallback to the awake rank monitor

**Rule.** If all K windows of a segment are near-blind, monitor the rest of the segment with E14's awake rank monitor (awake-B: smoothed conformal p-values of `|r_t|` among the segment's earlier returns, likelihood-ratio betting functions in rank space, payoffs divided by 1 + b, `docs/theory_running_tau.md` 4.1, 4.2 and Remark B.3), started with an empty reference set at `T_K`, with the same A and, with P1, restarted at the same period boundaries. After its alarm, the next segment calibrates again with window 1.

**Theorem 3 (the mixed guarantee) [PROVED].** Let `F = {all K windows are near-blind} ∈ F_{T_K}`. Suppose that on `Fᶜ` the tolerance holds relative to the kept window (horizon h), and that on F, conditionally on `F_{T_K}`, the returns observed after `T_K` are exchangeable (and the tie-breaking draws are independent uniforms). Then:

1. `P(an alarm within the first m monitored steps of the segment) ≤ P(Fᶜ ∩ Gᶜ) + m / A ≤ Σ_k δ_k + m / A ≤ δ + m / A` (h = m).
2. With period restarts: for every period `(b, b']` of the segment with `b' - b ≤ m`, `P(an alarm in it) ≤ δ + m / A` (h up to b'); on F the bound `m / A` holds given `F_{T_K}` and the p-values before b.
3. If the hypotheses hold in every segment (h = ∞; exchangeability of the whole remainder of a fallback segment), `E[N_M] ≤ (δ + M / A) / (1 - δ)`.

*Proof.* 1. Split on `Fᶜ ∩ Gᶜ`, `Fᶜ ∩ G` and F. The first has probability at most `Σ_k δ_k` by Theorem 1 (on `Fᶜ` some window k ≤ K is kept). On `Fᶜ ∩ G`, Step 3 of Theorem 6.4 with `E_s = Fᶜ ∩ G ∩ {T = s}` gives `m P(Fᶜ ∩ G) / A`. On F: given `F_{T_K}`, Proposition B.1 of `docs/theory_running_tau.md` applied to the exchangeable returns after `T_K` makes the p-values i.i.d. uniform and independent of `F_{T_K}`, so in `H_t = F_{T_K} ∨ σ(p-values up to t)` each payoff has conditional mean `1 / (1 + b) ≤ 1` and the two-sided Shiryaev-Roberts average N satisfies `E[N_t | H_{t-1}] ≤ N_{t-1} + 1`; Ville's inequality on `[T_K, T_K + m]` given `F_{T_K}` gives `m P(F) / A`. The three terms add to at most `Σ_k δ_k + m / A`. 2. The same split with the period's supermartingale (Corollary 1 on `Fᶜ`; on F a restart at b sets N to 0 and the argument runs in H given `H_b`). 3. Corollary 6.6's single supermartingale W mixes the filtrations F and H, so we use a per-segment argument instead. For segment j let `R_j` be the number of monitored steps still to come when it starts monitoring before the M-th one (known at its start of monitoring), and `ℓ_j ≤ R_j` its monitored steps among the first M; `Σ_j ℓ_j ≤ M`. On its good branch (`Fᶜ_j ∩ G_j`, or `F_j`), `1_E (N_t - (t - T_j))` is a supermartingale (in F, or in H given `F_{T_j}`; restarts only lower N), and optional stopping at `min(alarm, T_j + R_j)` gives `A P(E and segment j's alarm is among the first M monitored steps) ≤ E[1_E ℓ_j]`. On its bad branch, `P(Fᶜ_j ∩ Gᶜ_j | F_{S_j}) ≤ δ` by Theorem 1, and only segments that start before M monitored steps have passed (`C_j`, known at `S_j`) can contribute. So `E[N_M] ≤ δ E[#{j : C_j}] + M / A ≤ δ (1 + E[N_M]) + M / A`. ∎

The two hypotheses are about different things: a tolerance on conditional truncated variances for the cs branch, exchangeability of the returns after the fallback time for the rank branch. Neither needs the returns of the blind windows to be exchangeable, since the rank monitor's reference set starts empty at the fallback.

### 4.3 A grid of truncation levels (not implemented)

`docs/math.md` 6.4, fix 2: fix levels `τ_1 < ... < τ_J` when the CS part starts (for example multiples of the τ-part scale), run J sequences at level δ / J each, and choose τ at the end of calibration by any rule, for example the level with the largest headroom `τ_j² / m↑_j`. All J bounds hold together with probability at least 1 - δ, so the choice is free and no calibration time is lost. Proved there; the cost is `log(J / δ)` in place of `log(1 / δ)` in every segment (J = 5, δ = 0.01: 6.2 against 4.6), and it cannot help when the CS part is more volatile than the largest level allows.

## 5. What is proved, checked and assumed

* **Proved here:** Lemma 1, Propositions 1 and 2, Lemma 2, Theorem 1, Proposition 3, Theorem 2, Proposition 4, Theorem 3, the upper bound of Remark 6.
* **Sketch only:** the lower-bound side of Remark 6.
* **Checked numerically** (`tests/test_guard.py`, E24): the blind condition on constructed data; the retry and spending bookkeeping; identity with the unguarded monitor when K = 1, and when no window is blind with `δ_1` equal to the unguarded δ; Monte Carlo null checks of the first-m and per-period bounds on i.i.d. Gaussian and Student t3 data with a step in volatility; Proposition 3's coverage inflation (`results/e24_selection.csv`); false alarms, blind shares, detection and calibration cost in E24's synthetic designs.
* **Assumed, not proved:** the tolerance (or, for the rank fallback, exchangeability after the fallback) on real data. As in E22 and E23, it fails on most real days, so real-data alarm rates are not false alarm rates.

## 6. Recommended guard

See `docs/report_parts/E24.md` for the evidence. The recommendation, fixed after the synthetic experiments and before the real-data reruns, is in section 6.1 below.

### 6.1 Recommendation

(Filled in iteration 2 of `docs/iterations/E24.md`.)
