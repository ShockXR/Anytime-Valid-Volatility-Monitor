# False alarms after the first day, and two guarantees for every period (TH6)

Referee comment M2 (`docs/review/referee_report_paper.md`). Theorem B of `docs/math.md` bounds the probability of a false alarm within the first m steps after the statistics start (at the start of monitoring or after a recalibration). The paper and the README have read it as "5% per day". This note says what Theorem B and Corollary C do and do not bound for later days, proves two constructions that do give a bound for every period, shows by a counterexample where the second one stops, and relates the guarantees to the window-limited false alarm criteria of the change-detection literature. Code: `src/avmonitor/periodic.py`. Experiment: `experiments/th6_per_period.py` (outputs `results/th6_*.csv`, figures `results/figures/th6_*.png`, all **[SYNTHETIC]**). Tests: `tests/test_periodic.py`. Report part: `docs/report_parts/TH6.md`.

Labels: **[PROVED]** (a complete proof is given here, using only Ville's inequality, optional stopping at bounded stopping times and Harris's inequality for product measures), **[NUMERICAL]** (simulation), **[CITED]**.

## 0. Summary

1. **What Theorem B bounds.** Let T be the first false alarm and m the length of a day. Theorem B gives `P(T ≤ m) ≤ m / A` and, more generally, `P(T ≤ d m) ≤ d m / A`. It says nothing better about day d alone. This is not slack in the proof. The i.i.d. law `|r_t| = √m₀` (so `Y_t = 1`, inside the null class) makes both Shiryaev-Roberts statistics equal to t, so the default two-sided monitor alarms with certainty on day 40 (Proposition 1). For i.i.d. Gaussian returns the conditional probability of a first false alarm on day d, given none before, rises from 0.98% on day 1 to 4.70% on average over days 11 to 20 (390-step days). E16's bold play, a null stream that watches the detector, is at 5.25% on average over days 11 to 20, above 5%. For the referee's stream (calm for d - 1 days, then Gaussian) the probability of a false alarm on day d is 11.3% at d = 5, 28.0% at d = 10 and 65.1% at d = 20 [NUMERICAL, section 5]. The long-run statement is Corollary C: with a restart after every alarm, the expected number of false alarms over the first D days is at most `D m / A`, so 0.05 per day **on average**; on a single day it can be larger.
2. **Restarted Shiryaev-Roberts (Theorem 1, [PROVED]).** Restart both statistics at 0 at every period boundary. Then `P(false alarm in period i | everything before the period) ≤ m / A` for every period of at most m steps, under every null law of Theorem B, and the expected number of false alarms in each period is at most `m / A` too. Keeping at most κ of the statistic across a boundary gives `(m + κ) / A`. The cost is delay for changes late in a period (Proposition 2, section 5.4).
3. **The CUSUM-type e-detector (Theorem 2, [PROVED]).** `C_t = f_t max(C_{t-1}, 1)`, the maximum over start times of the product of bets (the e-CUSUM of Shin, Ramdas and Rinaldo 2024). If the returns are i.i.d. (or exchangeable and in the null class) and the bets are fixed, time reversal shows that `C_s` has the law of the running maximum of a test supermartingale, so `P(C_s ≥ x) ≤ 1 / x` at every s (Lemma 1). Hence one bet has `P(an alarm in any window of m steps) ≤ m / A`, whatever happened before. The monitor's average over bets has `(m + log A + H(w)) / A`, H the entropy of the weights, and the default two-sided rule (6 and 5 equal weights, each side at 2A) has `(m + log(2A) + (log 6 + log 5) / 2) / A`: **5.146%** for 390-step days at A = 7800 and **5.044%** for 1440-step days at A = 28,800. One side with i.i.d. returns also has the conditional hazard bounded by the same number.
4. **Martingale-difference nulls (Proposition 3, [PROVED]).** Without exchangeability the time reversal fails, and so does the window bound: an adaptive null stream (every draw has conditional mean exactly the baseline) can park the CUSUM statistic just below the threshold and then strike. For every m ≥ 1 and every A, the probability of a false alarm in the window `(s, s + m]` can be made at least `(1 - o(1)) / B` as s grows (B = 9.045 for Gaussian clipping at c = 3, so about 11%), while Theorem 2's bound goes to 0 as A grows. What remains under martingale differences is the conditional bound of Lemma 2.
5. **Delay (section 5.4, [NUMERICAL]).** At the same proved level of 5% per day, the restarted detector costs 6% to 44% (390-step days) and 2% to 14% (1440-step days) in mean delay against the plain Shiryaev-Roberts detector for Gaussian changes on day 3, mostly for changes late in the day, and stays below the renewal bound of Proposition 2; the e-CUSUM at its certified A costs 30% to 83% and 24% to 66%, because a maximum over start times is smaller than a sum and its proved bound is far above its actual false alarm rate. The plain detector's speed comes partly from the statistic it builds up under the null, which is what raises its later-day false alarm probability.

## 1. Set-up

As in `docs/math.md`, 5.1 to 5.4 and `docs/theory_robustness.md`, section 1.

* `Y_t = X_t / m₀ ∈ [0, B]` with the true baseline m₀ and `B = τ² / m₀ > 1`. `F_t` is the information up to t.
* **Null class of Theorem B** (two-sided, band 0): `E[Y_t | F_{t-1}] = 1` for every t (for one side, `≤ 1` or `≥ 1`). This allows any dependence that keeps the conditional mean fixed. We call it the **martingale-difference null**. The **exchangeable null** adds: for every s, the law of `(Y_1, ..., Y_s)` is invariant under permutations. I.i.d. returns are the main example; a return law whose shape is drawn once at the start, every shape having the same truncated mean, is another.
* **Bets.** Increase bets `λ_k ∈ [0, 1)` and decrease bets `φ_l ∈ [0, 1 / (B - 1))`, fixed in advance, with factors `f_{k,t} = 1 + λ_k (Y_t - 1)` and `g_{l,t} = 1 + φ_l (1 - Y_t)`, all nonnegative with conditional mean 1 under the null. Weights `w_k`, `v_l` (equal in the default: 6 increase bets 0.02 to 0.5, 5 decrease bets 0.1 to 0.95 times `1 / (B - 1)`).
* **Statistics of one bet.** Shiryaev-Roberts `R_t = f_t (R_{t-1} + 1) = Σ_{j ≤ t} Π_{i=j}^{t} f_i`; CUSUM type `C_t = f_t max(C_{t-1}, 1) = max_{j ≤ t} Π_{i=j}^{t} f_i` (with `R_0 = C_0 = 0`). Since a maximum of nonnegative terms is at most their sum, `C_t ≤ R_t`.
* **Two-sided rule.** Each side averages its bets (`M↑_t = Σ_k w_k R↑_{k,t}`, and the same with C); an alarm is raised when either side reaches 2A. Then `N_t = (M↑_t + M↓_t) / 2 ≥ A` at an alarm.
* **Days.** Day d is the steps `((d - 1) m, d m]`. For the first alarm T: `p_d = P(T ∈ day d)` and the hazard `h_d = P(T ∈ day d | T > (d - 1) m)`. With a restart after every alarm (a recalibration with the true baseline): `q_d = P(at least one alarm on day d)` and `e_d = E[number of alarms on day d]`.

What was proved before: Theorem B gives `p_1 ≤ m / A` and `p_1 + ... + p_D = P(T ≤ D m) ≤ D m / A`. Corollary C gives `e_1 + ... + e_D ≤ D m / A`.

## 2. The plain Shiryaev-Roberts detector has no per-day bound

**Proposition 1 [PROVED].** (a) For every s ≥ 0 and m ≥ 1, under the martingale-difference null, `P(an alarm in (s, s + m] | F_s) ≤ min{1, (N_s + m) / A}`, and therefore `P(an alarm on day d) ≤ min{1, d m / A}`.
(b) Over the i.i.d. null (even over i.i.d. laws), no bound on `p_d`, `h_d` or `q_d` below 1 holds for all d: the law `Y_t = 1` (that is `|r_t| = √m₀`; `X_t = m₀` needs `√m₀ ≤ τ`, which holds since `m₀ ≤ τ²`) gives `R↑_{k,t} = R↓_{l,t} = t` for every bet, so the two-sided rule alarms at exactly `t = ⌈2A⌉`. With A = 7800 and m = 390, or A = 28,800 and m = 1440, that is the last step of day 40: `p_40 = h_40 = q_40 = e_40 = 1`.

*Proof.* (a) For `s ≤ t ≤ s + m` put `S_t = N_t + (s + m - t)`. As in Lemma 1 of `docs/math.md`, `E[N_t | F_{t-1}] ≤ N_{t-1} + 1`, so S is a nonnegative supermartingale on `[s, s + m]` with `S_s = N_s + m`, and Ville's inequality applied given `F_s` bounds the chance that S reaches A. An alarm implies `N ≥ A`, hence `S ≥ A`. Take expectations and use `E[N_s] ≤ s`. (b) With `f = g = 1` the recursion gives `R_t = R_{t-1} + 1`. ∎

So Theorem B and Corollary C are the right statements for the Shiryaev-Roberts detector: a bound for the first m steps after a (re)start, and a bound on the long-run count. The example in (b) is consistent with Corollary C: one alarm in 40 days is 0.025 per day. Section 5 gives the size of the effect for natural laws.

## 3. Theorem 1: restart at every period boundary

Let `0 = b_0 < b_1 < b_2 < ...` be the period boundaries, stopping times with `b_{i+1} - b_i ≤ m` (fixed day starts in practice). At each boundary, set every bet's statistic to 0 (restart), or, more generally, scale each side's bets down by a common factor so that each side average is at most κ ≥ 0 (carry at most κ). Alarms are raised by the two-sided rule; the statistics may or may not restart after an alarm.

**Theorem 1 [PROVED].** Under the martingale-difference null of Theorem B, for every period i:

1. `P(an alarm in (b_i, b_{i+1}] | F_{b_i}) ≤ min{1, (m + κ) / A}` almost surely;
2. if the statistics restart after every alarm, `E[number of alarms in (b_i, b_{i+1}] | F_{b_i}] ≤ (m + κ) / A`.

Hence, with κ = 0 (restart) and `A = m / α`, every day has `q_d ≤ α`, `e_d ≤ α` and, because `{T > b_i} ∈ F_{b_i}`, the hazard `h_d ≤ α`, whatever happened on earlier days. A window of m consecutive steps that is not aligned with the boundaries meets at most two periods and has probability at most `2 (m + κ) / A`.

*Proof.* (1) At `b_i` each side average is at most κ, so `N_{b_i} ≤ κ`. Until the first alarm in the period, the statistic with and without restarts after alarms coincide, so it suffices to bound the chance that the statistic without restarts reaches `N ≥ A` in the period. For `b_i ≤ t ≤ b_i + m` put `S_t = N_t + (b_i + m - t)`; as in Proposition 1(a) it is a nonnegative supermartingale given `F_{b_i}`, with `S_{b_i} ≤ κ + m`. Ville's inequality given `F_{b_i}` gives the claim (the period ends at `b_{i+1} ≤ b_i + m`). (2) Let `W_t` be the sum of `N` at the alarms in `(b_i, t]` (each at least A), plus the current `N_t`, minus `t - b_i`. At a non-alarm step `E[W_t - W_{t-1} | F_{t-1}] = E[N_t | F_{t-1}] - N_{t-1} - 1 ≤ 0`; at an alarm the value of N moves into the sum and N restarts at 0, which leaves W unchanged at that step relative to the computation without restart. So W is a supermartingale with `W_{b_i} = N_{b_i} ≤ κ`. Optional stopping at the bounded time `b_{i+1}` gives `A E[count | F_{b_i}] - m ≤ E[W_{b_{i+1}} | F_{b_i}] ≤ κ`. ∎

The proof is that of Theorem B and Corollary C, started afresh at each boundary; the point is only that a restart makes the start of every period look like the start of monitoring. With a confidence-sequence baseline (`docs/math.md`, Theorem 6.4) the bound becomes `δ + m / A` for every period of a segment, since the event that the confidence bound fails is one event per segment. `PeriodicVolatilityMonitor(period=m)` implements the restart on the stream's clock (for example at every session open).

**What the restart costs in delay.** Let the returns after a change at step ν be i.i.d. with some post-change law and independent of the past. Let `D₀` be the delay (steps up to and including the alarm) of the two-sided detector started at 0 at the change, and `p = P(D₀ ≤ m)`.

**Proposition 2 [PROVED].** If the change falls k steps after the start of a period of length m (`0 ≤ k < m`) and `p > 0`, the delay `D_r` of the detector restarted at every period boundary satisfies

`E[D_r] ≤ E[min(D₀, m - k)] + P(D₀ > m - k) E[min(D₀, m)] / p`,

with equality for k = 0 in the form `E[D_r] = E[min(D₀, m)] / p`. The plain detector has `D ≤ D₀` pathwise. If `p = 0`, the restarted detector can detect the change only in the period in which it occurs.

*Proof.* The recursion `R_t = f_t (R_{t-1} + 1)` is nondecreasing in `R_{t-1}` because `f_t ≥ 0`, and every statistic is at least 0 at the change, so on the same data every statistic of the restarted (or the plain) detector is at least that of the detector started at 0 at the change, until the next restart. Hence on `{D₀ ≤ m - k}` the restarted detector alarms by step `D₀`, and the plain detector always does. Otherwise the restarted detector is set to 0 at the next boundary, `m - k` steps after the change. From then on, period after period, it is the detector started at 0 on i.i.d. post-change data, independent of the first partial period, and a period contains an alarm with probability p. The number of periods without an alarm is geometric, so the expected delay counted from that boundary is `(1 - p) m / p + E[D₀ | D₀ ≤ m] = E[min(D₀, m)] / p`. ∎

So the restart costs little when a change is usually detected within a day from scratch (p near 1) and does not fall in the last `D₀` or so steps of a day; it costs about one more fresh detection for changes late in the day; and it loses changes whose detection from scratch needs more than a day. Keeping a carry κ (scaling the statistic down instead of to 0) keeps part of the evidence across the boundary at the price of `(m + κ) / A` instead of `m / A`. The plain detector is faster than `D₀` because it enters the change with the statistic accumulated under the null, `E[R_ν] ≤ ν`; that head start is the same accumulation that raises its false alarm probability on later days (section 2). Section 5.4 evaluates the bound (`results/th6_delay_renewal.csv`).

## 4. The CUSUM-type e-detector

### 4.1 Time reversal

**Lemma 1 [PROVED].** Assume the exchangeable null and fixed bets. For every s ≥ 1 the vector `(C_{k,s})_k` of CUSUM statistics of the bets has the same law as `(max_{1 ≤ j ≤ s} W_{k,j})_k`, where `W_{k,j} = Π_{i ≤ j} f_{k,i}` is the wealth of bet k started at time 1. In particular `P(C_{k,s} ≥ x) ≤ 1 / x` for every x > 0, every bet and every s, and likewise for the decrease bets.

*Proof.* Write the factors at time i as a vector `F_i = φ(Y_i)` for a fixed function φ (the bets are fixed). Then `C_{k,s} = max_{j ≤ s} Π_{i=j}^{s} F_{k,i} = Ψ_k(F_1, ..., F_s)` with `Ψ_k(v_1, ..., v_s) = max_j Π_{i=j}^{s} v_{k,i}`. Reversing the order of the arguments gives `Ψ_k(v_s, ..., v_1) = max_{j' ≤ s} Π_{i=1}^{j'} v_{k,i}`. By exchangeability `(F_1, ..., F_s)` and `(F_s, ..., F_1)` have the same law, which gives the identity in law, jointly over k. Under the null, `E[F_{k,j} | F_{j-1}] ≤ 1`, so `W_k` is a nonnegative supermartingale with `W_{k,0} = 1`, and Ville's inequality gives `P(max_{j ≤ s} W_{k,j} ≥ x) ≤ 1 / x`. ∎

The same reversal applied to the Shiryaev-Roberts statistic gives `R_{k,s} =ᵈ Σ_{j ≤ s} W_{k,j}`, a sum instead of a maximum, whose tail is not `1 / x` uniformly in s; Proposition 1(b) is the extreme case.

### 4.2 Theorem 2

**Lemma 2 (conditional window bound) [PROVED].** Under the martingale-difference null, with fixed or predictable bets, for one side averaging K bets with weights w and alarming at level A,

`P(an alarm in (s, s + m] | F_s) ≤ min{1, (Ū_s + m - 1) / A}`, where `Ū_s = Σ_k w_k max(C_{k,s}, 1)`.

*Proof.* For each bet put `Z_{k,s+1} = f_{k,s+1} max(C_{k,s}, 1) = C_{k,s+1}` and `Z_{k,t} = f_{k,t} (Z_{k,t-1} + 1)` for t > s + 1. By induction `C_{k,t} = f_{k,t} max(C_{k,t-1}, 1) ≤ f_{k,t} (C_{k,t-1} + 1) ≤ Z_{k,t}`. With `Z̄_t = Σ_k w_k Z_{k,t}`, the process `S_s = Ū_s + m - 1`, `S_t = Z̄_t + (s + m - t)` for `s < t ≤ s + m` is a nonnegative supermartingale: `E[S_{s+1} | F_s] ≤ Ū_s + m - 1` and `E[Z̄_t | F_{t-1}] ≤ Z̄_{t-1} + 1` afterwards. Ville's inequality given `F_s`, and `Σ_k w_k C_{k,t} ≤ S_t`. ∎

**Lemma 3 [PROVED].** If `U ≥ 1` and `P(U ≥ x) ≤ 1 / x` for x ≥ 1, then for `0 < w ≤ A`, `E[min(1, w U / A)] ≤ (w / A)(1 + log(A / w))`.

*Proof.* `E[min(1, wU/A)] = ∫₀¹ P(U ≥ vA / w) dv ≤ w / A + ∫_{w/A}^{1} w / (vA) dv`. ∎

**Theorem 2 [PROVED].** Assume the exchangeable null and fixed bets. Let s ≥ 0 and m ≥ 1.

(a) *One bet, one side at level A.* `P(C_t ≥ A for some t ∈ (s, s + m]) ≤ m / A`.
(b) *One side, average of K bets with weights w, level A.* `P(Σ_k w_k C_{k,t} ≥ A for some t ∈ (s, s + m]) ≤ (m + log A + H(w)) / A`, with `H(w) = -Σ_k w_k log w_k ≤ log K`.
(c) *The default two-sided rule* (each side at 2A): `P(an alarm in (s, s + m]) ≤ (m + log(2A) + (H(w) + H(v)) / 2) / A`. For 6 and 5 equal weights: 5.146% at m = 390, A = 7800, and 5.044% at m = 1440, A = 28,800. The A that makes it exactly 5% is 8027.7 and 29,053.4 (`threshold_for_window`).
(d) The same bounds hold for "at least one alarm in the window" when the monitor restarts after every alarm, and for the first alarm.
(e) *Hazard, one side, i.i.d. returns.* `P(T ∈ (s, s + m] | T > s)` is bounded by the bound of (a) or (b). For the two-sided rule, `P(T ∈ (s, s + m] | T > s) ≤ (m - 1) / A + (1 + log(2A) + (H(w) + H(v)) / 2) / (A P(T > s))`.

*Proof.* (a) Union bound over the m steps of the window and Lemma 1 at each step. (b) By Lemma 2, `P(...) ≤ E[min{1, (Ū_s + m - 1) / A}] ≤ (m - 1) / A + E[min(1, Ū_s / A)]`, and `min(1, Σ_k a_k) ≤ Σ_k min(1, a_k)` for `a_k ≥ 0` gives `E[min(1, Ū_s / A)] ≤ Σ_k E[min(1, w_k U_k / A)]` with `U_k = max(C_{k,s}, 1)`, which satisfies Lemma 3 by Lemma 1. Summing `(w_k / A)(1 + log(A / w_k))` over k gives `(1 + log A + H(w)) / A`. (c) Union over the two sides, each at level 2A. (d) After a restart at time a, the restarted statistic is `max_{a < j ≤ t} Π_{i=j}^{t} f_i ≤ C_t`, bet by bet, so every alarm of the restarted monitor in the window is a time at which the statistic without restarts is at or above the threshold. The first alarm is such a time. (e) With i.i.d. returns `Y_1, Y_2, ...` are independent. On the increase side every `C_{k,t}` is a nondecreasing function of `(Y_1, ..., Y_t)` (each factor is nondecreasing in Y and nonnegative, and `max(·, 1)` and products of nonnegative nondecreasing functions are nondecreasing), and `{T > s}` is a nonincreasing event. Harris's inequality for product measures gives `P(C_t ≥ A, T > s) ≤ P(C_t ≥ A) P(T > s)` and `E[min(1, w_k U_k / A) 1{T > s}] ≤ E[min(1, w_k U_k / A)] P(T > s)`. Applying the proofs of (a) and (b) on `{T > s}` (Lemma 2 holds given `F_s`, and `{T > s} ∈ F_s`) gives the one-sided claim; the decrease side is the same with the order reversed. For two sides the statistics are monotone in opposite directions and Harris's inequality does not apply; dividing the Lemma 3 part by `P(T > s)` gives the stated form. ∎

**Remarks.**

1. *Sharpness of (a).* Over the class of all nonnegative factors with mean 1 (an unclipped statistic), the law `f ∈ {0, A}` with `P(f = A) = 1 / A` alarms at a step exactly when `f = A`, so the window probability is `1 - (1 - 1/A)^m`, which is m / A to first order. The extra `log A + H(w)` in (b) comes from bounding the average of dependent maxima; whether it can be removed for the average is open. It adds 2.9% (relative) to the one-day bound at A = 7800 and 0.9% at A = 28,800.
2. *Other ways to mix the bets.* (i) The maximum over start times of the bet average, `max_j Σ_k w_k Π_{i=j}^{t} f_{k,i}`, is a single CUSUM of the mixture wealth, so Lemma 1 applies to it directly and the window bound is exactly m / A (two sides at 2A: m / A). It has no finite recursion (a start j can be dropped once another start dominates it in every bet). (ii) The maximum over bets of the per-bet statistics needs a union bound: K m / A at level A.
3. *Predictable bets* (the adaptive bet of `detectors.betting_edetector`) break the reversal: the bet at time i depends on `Y_1, ..., Y_{i-1}`, and the reversed sequence is not a betting sequence. Lemma 2 still holds; Theorem 2 is not claimed.
4. *Estimated baselines.* With the confidence-sequence baseline of `docs/math.md` section 6 and i.i.d. returns, the monitored draws are, given the calibration window, i.i.d. with conditional mean on the null side on the event G that the bounds hold, so the window bound becomes `δ + (bound of (c))`. The plug-in baseline is not covered (Proposition A.6 of `docs/theory_running_tau.md`).
5. *Theorems A and B and Corollary C hold for the CUSUM statistic too*, since `C_t ≤ R_t` and `E[f_t max(C_{t-1}, 1) | F_{t-1}] ≤ C_{t-1} + 1`.

### 4.3 Martingale-difference nulls: the reversal fails

**Proposition 3 [PROVED].** Fix the bets (one bet, or the default two sides) and A. There is a null stream with `E[Y_t | F_{t-1}] = 1` and `Y_t ∈ [0, B]` at every step, and stopping times ρ with `P(ρ ≤ s) → 1` as s → ∞, such that `T > s` almost surely and `P(T = s + 1) ≥ P(ρ ≤ s) / B`. Hence, for every m ≥ 1 and every A,

`liminf_{s → ∞} sup_{null laws} P(an alarm in (s, s + m] | T > s) ≥ 1 / B`,

while the bound of Theorem 2 tends to 0 as A grows. No bound on the window probability that tends to 0 with A holds over the martingale-difference null.

*Proof.* *Hunt.* Draw `Y_t` from `{0, B}` with `P(Y_t = B) = 1 / B` (mean 1), independently, until the first time ρ at which one draw could take a side to 2A: on the increase side `Σ_k w_k (1 + λ_k (B - 1)) max(C↑_{k,ρ}, 1) ≥ 2A`, or on the decrease side `Σ_l v_l (1 + φ_l) max(C↓_{l,ρ}, 1) ≥ 2A`. *Park.* From ρ + 1 to s, play `Y = 1`. *Strike.* At s + 1, draw from `{0, B}` again. Every draw has conditional mean 1 and lies in `[0, B]` (`|r| = 0`, `|r| = √m₀` or `|r| ≥ τ`).
No alarm before s + 1: before ρ, no value of `Y ∈ [0, B]` takes a side to 2A, because the increase statistics are nondecreasing in Y (largest at B) and the decrease statistics nonincreasing (largest at 0). At ρ the statistic is below 2A on both sides by the same argument applied at ρ - 1, and since `f(B) ≥ 1` and `g(0) ≥ 1`, `Σ_k w_k max(C↑_{k,ρ}, 1) ≤ Σ_k w_k f_k(B) max(C↑_{k,ρ-1}, 1) < 2A`; playing `Y = 1` replaces each `C` by `max(C, 1)` and then leaves it unchanged, so no alarm while parked.
ρ is finite almost surely: from any state, n consecutive draws of B make the increase side at least `Σ_k w_k (1 + λ_k (B - 1))^n`, so ρ comes no later than the end of the first run of `n*` consecutive B's, with `n*` depending only on A and the bets, and such runs occur almost surely (disjoint blocks of `n*` steps are independent, each all-B with probability `B^{-n*} > 0`).
Strike: if the increase side parked, the draw B (probability 1 / B) raises an alarm at s + 1; if the decrease side parked, the draw 0 (probability `1 - 1/B ≥ 1/B` since B ≥ 2) does. ∎

What fails is Lemma 1: under the parking stream `P(C_s ≥ x)` tends to 1 for every x below the parking level. Lemma 2 is still true and predicts the strike: at a parked state `Ū_s` is close to the level, so its bound is large. Corollary C also still holds, as it does for every null stream: the expected number of false alarms over D days is at most `D m / A`. Section 5.3 shows the effect at full scale; at A = 7800 the hunt is slow, so the rate passes the i.i.d. bound only on day 17 (390-step days).

## 5. Numerical check [NUMERICAL]

Default detector, true baseline (τ and m₀ computed exactly for each law, `results/th6_laws.csv`), 20 days. Full tables: `docs/report_parts/TH6.md` and the CSV files named there.

### 5.1 The plain detector after the first day

The conditional probability `h_d` of a first false alarm on day d given none before, in percent, for the plain Shiryaev-Roberts detector; averages of `h_d` over the days shown; 40,000 paths per stream (bold play 30,000); standard errors of single `h_d` at most 0.17 points (i.i.d. laws), 0.21 (bold play) and 0.24 (calm stream). The last column is the expected number of false alarms per day over days 1 to 20 with a restart after every alarm, which Corollary C bounds by `m / A = 0.05`. Files: `results/th6_hazard_390.csv`, `results/th6_hazard_1440.csv`, `results/th6_adversary.csv`.

| day length, A | stream | day 1 | days 2 to 5 | days 6 to 10 | days 11 to 20 | largest `h_d` (day) | false alarms per day, days 1 to 20 |
|---|---|---|---|---|---|---|---|
| 390, 7800 | Gaussian | 0.98 | 3.02 | 4.13 | 4.70 | 5.18 (19) | 0.0376 |
| 390, 7800 | Student t3 | 1.49 | 3.39 | 4.11 | 4.68 | 5.06 (19) | 0.0392 |
| 390, 7800 | jumps | 0.99 | 3.03 | 4.05 | 4.60 | 5.03 (19) | 0.0371 |
| 390, 7800 | ticks | 1.52 | 3.29 | 4.16 | 4.66 | 4.88 (18) | 0.0388 |
| 390, 7800 | bold play (E16) | 2.53 | 4.18 | 4.91 | 5.25 | 5.44 (14) | 0.0463 |
| 390, 7800 | `Y = 1` for 9 days, then Gaussian | 0 | 0 | 28.16 on day 10 | 5.47 | 28.16 (10) | 0.0370 |
| 1440, 28,800 | Gaussian | 2.04 | 3.64 | 4.50 | 4.68 | 4.94 (18) | 0.0416 |
| 1440, 28,800 | Student t3 | 2.44 | 3.92 | 4.46 | 4.57 | 4.77 (13) | 0.0419 |
| 1440, 28,800 | jumps | 2.07 | 3.67 | 4.34 | 4.67 | 4.93 (17) | 0.0411 |
| 1440, 28,800 | ticks | 2.37 | 3.95 | 4.58 | 4.78 | 4.95 (13) | 0.0430 |
| 1440, 28,800 | bold play (E16) | 3.13 | 4.56 | 5.00 | 5.02 | 5.27 (14) | 0.0479 |
| 1440, 28,800 | `Y = 1` for 9 days, then Gaussian | 0 | 0 | 34.03 on day 10 | 5.28 | 34.03 (10) | 0.0407 |

For i.i.d. returns the hazard rises over the first week and then levels off just below `m / A`; no single i.i.d. day is significantly above 5%. Bold play, which is inside the null class, is above 5% on average over days 11 to 20 at 390-step days: 5.25%, with a pooled standard error of about 0.05 points (the pooling over days 11 to 20 was chosen after seeing the table). So "at most 5% on every day" is false for the plain detector on a null stream that is not artificial. Every long-run count is below 0.05, as Corollary C requires.

The referee's family (`results/th6_calm_family.csv`, 40,000 paths per day d, standard errors at most 0.25 points): `Y = 1` for d - 1 days, then one Gaussian day. After the calm days every Shiryaev-Roberts bet sits at `(d - 1) m` exactly.

| day d | 1 | 2 | 5 | 10 | 15 | 20 |
|---|---|---|---|---|---|---|
| Ville's bound `min(1, d m / A)` | 5 | 10 | 25 | 50 | 75 | 100 |
| Shiryaev-Roberts, 390-step days | 0.90 | 3.04 | 11.33 | 28.02 | 46.82 | 65.07 |
| Shiryaev-Roberts, 1440-step days | 1.82 | 5.06 | 15.13 | 34.30 | 52.84 | 70.82 |
| restarted daily, 390 (bound 5) | 0.90 | 0.95 | 0.91 | 0.85 | 0.97 | 0.93 |
| e-CUSUM, 390 (bound 5.146) | 0.02 | 0.05 | 0.04 | 0.03 | 0.08 | 0.07 |

With a two-point day `{0, B}` instead of a Gaussian one the plain detector reaches 70.09% (390) and 73.16% (1440) on day 20.

### 5.2 The per-period constructions on the same streams

Probability `q_d` of at least one false alarm on day d (restart after every alarm), smallest and largest over days 1 to 20, in percent. The e-CUSUM column "no restart" is the probability that the statistic is at or above the threshold at some step of day d without any restart, which is the event Theorem 2(c) bounds. Same files as 5.1.

| day length | stream | restarted daily (bound 5) | e-CUSUM (bound 5.146 or 5.044) | e-CUSUM, no restart |
|---|---|---|---|---|
| 390 | Gaussian | 0.87 to 1.03 | 0.04 to 0.10 | at most 0.12 |
| 390 | Student t3 | 1.47 to 1.65 | 0.08 to 0.17 | at most 0.18 |
| 390 | jumps | 0.93 to 1.16 | 0.05 to 0.10 | at most 0.12 |
| 390 | ticks | 1.38 to 1.63 | 0.06 to 0.14 | at most 0.15 |
| 390 | bold play | 2.23 to 2.66 | 0.28 to 0.38 | not run |
| 390 | calm 9 days, then Gaussian | 0 to 1.03 | 0 to 0.10 | at most 0.11 |
| 1440 | Gaussian | 1.78 to 2.08 | 0.06 to 0.11 | at most 0.11 |
| 1440 | Student t3 | 2.21 to 2.57 | 0.10 to 0.17 | at most 0.18 |
| 1440 | jumps | 1.90 to 2.10 | 0.05 to 0.10 | at most 0.10 |
| 1440 | ticks | 2.29 to 2.57 | 0.09 to 0.15 | at most 0.15 |
| 1440 | bold play | 3.11 to 3.52 | 0.26 to 0.40 | not run |
| 1440 | calm 9 days, then Gaussian | 0 to 2.06 | 0 to 0.11 | at most 0.11 |

Both constructions keep their bound on every day of every stream. The restarted detector repeats its first day: its `q_d` is flat in d at the plain detector's day-1 rate. The e-CUSUM is far below its bound; its proved bound is a union over the steps of the window and is loose for natural laws (it is attained only by laws with a large jump, Remark 1 after Theorem 2).

### 5.3 The counterexample

Park and strike against the default two-sided e-CUSUM, Gaussian clipping (`B = 9.045`): hunt with `{0, B}` draws until one draw could reach the threshold, park with `Y = 1`, and on day d play bold (E16) from the parked or unparked state. 30,000 paths; standard errors at most 0.14 points (`results/th6_park_strike.csv`).

| day d | 1 | 5 | 10 | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|---|
| share parked before day d, 390 | 0 | 2.7 | 5.9 | 9.1 | 9.7 | 10.3 | 10.9 | 11.5 | 12.1 |
| P(false alarm on day d), 390 (Theorem 2 bound 5.146 if i.i.d.) | 0.26 | 1.60 | 3.19 | 4.70 | 5.10 | 5.43 | 5.67 | 6.07 | 6.16 |
| Lemma 2 averaged over paths, 390 | 5.00 | 6.38 | 8.01 | 9.53 | 9.81 | 10.12 | 10.42 | 10.70 | 11.00 |
| P(false alarm on day d), 1440 (bound 5.044 if i.i.d.) | 0.28 | 1.63 | 3.42 | 4.96 | 5.30 | 5.52 | 5.90 | 6.12 | 6.44 |

At full scale the hunt is slow, so the probability passes the i.i.d. bound only on day 17 (390-step days; 5.43%, 2.2 standard errors above 5.146%) and is 7.3 standard errors above it on day 20 (6.16%). Scaled down to periods of m = 20 steps and A = 400 (`results/th6_park_scaled.csv`, 20,000 paths), where the i.i.d. bound is 7.1%, the share parked grows to 57% by period 120 and the probability of a false alarm in period d grows to 2.8% (d = 10), 11.3% (40), 20.3% (80) and 27.0% (120), well above the lower bound of Proposition 3 (the parked share times `1 / B` = 11.1%, so 6.3% at period 120, and at least 11.1% in the limit, when every path is parked). The mean of Lemma 2's conditional bound is above the measured probability in every row, as it must be.

### 5.4 Delay cost

A change of the return standard deviation by the factor shown, at a uniformly random step of day 3 (Gaussian returns, true baseline), 20,000 paths per cell, horizon 12 days (`results/th6_delay.csv`, `results/th6_delay_by_position.csv`, `results/th6_delay_renewal.csv`). The delay counts the steps from the change to the first alarm at or after it, including the alarm step; every detector found every change within the horizon. Standard errors of the mean delays are at most 0.7 steps. The plain detector uses A (a proved bound for the first day after a (re)start only); the restarted detector A; the restarted detector with carry m uses 2A, so `(m + m) / (2A) = 5%`; the e-CUSUM uses A = 8027.7 (390) or 29,053.4 (1440), where its window bound is 5% under i.i.d. returns. `D₀` is the fresh-start delay of Proposition 2 (20,000 further paths), and the Proposition 2 bound is averaged over the uniform position of the change.

| day length | sd factor | plain SR, mean delay | restarted daily / plain | restarted, carry m / plain | e-CUSUM certified / plain | fresh start `E[min(D₀, m)] / p` | restarted daily: measured / Proposition 2 bound | restarted daily by position: first half / 50% to 90% / last 10% |
|---|---|---|---|---|---|---|---|---|
| 390 | 0.5 | 70.8 | 1.26 | 1.23 | 1.72 | 101.9 | 89.2 / 115.0 | 82 / 91 / 122 |
| 390 | 0.67 | 97.5 | 1.32 | 1.26 | 1.75 | 136.8 | 128.4 / 160.4 | 112 / 142 / 156 |
| 390 | 0.8 | 161.9 | 1.44 | 1.34 | 1.83 | 217.8 | 232.6 / 274.9 | 196 / 275 / 239 |
| 390 | 1.25 | 113.4 | 1.21 | 1.22 | 1.61 | 135.0 | 137.1 / 151.8 | 124 / 149 / 155 |
| 390 | 1.5 | 43.2 | 1.09 | 1.14 | 1.39 | 51.9 | 47.3 / 54.7 | 45 / 46 / 63 |
| 390 | 2 | 20.1 | 1.06 | 1.11 | 1.30 | 24.7 | 21.2 / 25.4 | 21 / 20 / 27 |
| 1440 | 0.5 | 87.4 | 1.09 | 1.13 | 1.58 | 118.3 | 94.8 / 123.1 | 92 / 88 / 135 |
| 1440 | 0.67 | 120.7 | 1.10 | 1.14 | 1.60 | 159.7 | 132.3 / 168.4 | 127 / 123 / 198 |
| 1440 | 0.8 | 201.7 | 1.14 | 1.16 | 1.66 | 257.9 | 229.3 / 279.9 | 214 / 226 / 323 |
| 1440 | 1.25 | 139.7 | 1.07 | 1.13 | 1.50 | 163.0 | 149.8 / 170.0 | 145 / 146 / 193 |
| 1440 | 1.5 | 51.8 | 1.03 | 1.10 | 1.32 | 60.1 | 53.6 / 61.1 | 53 / 52 / 64 |
| 1440 | 2 | 23.7 | 1.02 | 1.08 | 1.24 | 27.9 | 24.1 / 28.1 | 24 / 24 / 26 |

Three things explain the table. (i) The plain detector is 14% to 31% faster than a fresh start, because it enters the change with the statistic built up under the null over two days; that is the same accumulation that raises its false alarm probability on later days (section 2). (ii) The restarted detector is close to a fresh start and below the bound of Proposition 2 in every cell; its cost sits in changes whose remaining part of the day is shorter than a typical detection (for the 0.8 factor at 390 steps, changes in the second half of the day, where the restart discards most of the evidence). Over a 1440-step day the cost is 2% to 14%. A carry of m, at 2A, helps a little for decreases with 390-step days and costs a little elsewhere. (iii) The e-CUSUM is 1.24 to 1.83 times slower: it is a maximum, not a sum, over start times, so it is smaller than a fresh Shiryaev-Roberts statistic, and its proved bound (a union over the m steps of a window) is far above its actual rate (0.04% to 0.17% per day under i.i.d. returns against 5.146%). At equal proved level the restarted Shiryaev-Roberts detector is the faster of the two per-period constructions, and it covers the martingale-difference null, where the e-CUSUM's window bound fails (Proposition 3).

The Proposition 2 part (`renewal`) was added after the first delay results were seen, to explain them; its bound is proved, but its choice as a summary is post hoc.

## 6. Relation to window-limited false alarm criteria

* **Lai (1998)** studies detection under a constraint on the probability of a false alarm in a window, in the unconditional form `sup_k P_∞(k ≤ T < k + m_α) ≤ α`, with a window `m_α` chosen according to α, and shows that a window-limited CUSUM rule is asymptotically optimal under it, with delays of order `|log α|` divided by an information number [CITED; see the verification note below]. Theorem 2 gives exactly this criterion for the e-CUSUM, with `α = (m + log(2A) + H̄) / A`, uniformly over the exchangeable null class instead of for one pre-change law. Theorem 1 gives it for the restarted Shiryaev-Roberts detector for windows aligned with the periods (and `2 (m + κ) / A` for the others), over the martingale-difference class. The plain Shiryaev-Roberts detector satisfies it only for the first window (k = 1, Theorem B); Proposition 1(b) shows that it fails for later k even over i.i.d. nulls. Lai's asymptotics let α go to 0 with the window growing; in our use m is fixed by the trading day and α = 5%, so they do not decide which rule is faster here, and section 5.4 compares the rules by simulation.
* **The local conditional version**, `sup_k P_∞(T < k + m | T ≥ k) ≤ β`, is the per-window hazard. Pergamenchtchikov and Tartakovsky (2018) build detection procedures around a maximal "local" probability of a false alarm within a period of fixed length and show that the Shiryaev-Roberts procedure is asymptotically optimal under it for general dependent data, with a known pre-change law [CITED]. Theorem 1 bounds the hazard for every period by m / A (with restarts), Theorem 2(e) for one side by the window bound.
* **Mei (2008)** asks whether the average run length to false alarm is always an informative criterion and shows that, with mixture pre-change models, detection schemes with finite detection delays can have infinite average run lengths to false alarm; the paper discusses other criteria, among them quantiles of the run length and an expected false alarm rate [CITED]. Our setting has the same feature: the null is composite and contains mixtures (an exchangeable law is a mixture of i.i.d. laws). Theorem A (`E[T] ≥ A`) holds for every null law, but it coexists with per-day hazards that rise to 4.6% under i.i.d. Gaussian returns and to 1 on day 40 under `Y = 1` (Proposition 1). Corollary C (at most one false alarm per A monitored steps in expectation, with restarts) is a statement about the rate of false alarms over time, which is the kind of criterion Mei argues for; we could not open the paper to compare his exact definition.

**How the references were checked (3 November 2026).** Web searches only; arxiv.org, IEEE Xplore, Taylor and Francis and the authors' pages could not be opened from this machine, so each statement above rests on the abstract, keywords and the parts of the text shown by the search engine, and on secondary sources that cite the paper. Lai (1998): title, journal, volume, issue, pages and the form of the criterion (`sup_k P_∞(k ≤ T < k + m)`, from a 2021 *Statistica Sinica* paper that cites it) verified this way; the exact conditions on the window `m_α` were not verified. Mei (2008): title, journal, volume, issue, pages, DOI, the main result on mixture pre-change models and the alternative criteria named in the keywords verified; the paper's examples and exact definitions were not. Pergamenchtchikov and Tartakovsky (2018): journal, volume, pages and the statement quoted from the abstract verified; whether their criterion conditions on no earlier alarm was not verified. Shin, Ramdas and Rinaldo (2024): that they define e-SR (sum over start times) and e-CUSUM (maximum) verified from the abstract and secondary sources; their mixture construction was not checked in the full text.

## 7. Status

* **Proved here:** Proposition 1, Theorem 1, Proposition 2, Lemmas 1 to 3, Theorem 2 (a) to (e), Proposition 3.
* **Checked numerically:** section 5 and `tests/test_periodic.py` (engine against the package's detector; Lemma 1 as an identity and as a tail bound; Theorem 1 on the calm stream; Theorem 2(a) on a late window; Lemma 2 and Proposition 3 on a small scale; the streaming monitor).
* **Open:** whether `log A + H(w)` in Theorem 2(b) can be removed for the average over bets; a window bound for predictable bets under an exchangeable null; the two-sided hazard without the factor `1 / P(T > s)`; how close any null stream comes to Theorem 1 after the first period (bold play reaches about half of it on day 1, E16).

## 8. References

* Harris, T. E. (1960). A lower bound for the critical probability in a certain percolation process. *Proceedings of the Cambridge Philosophical Society* 56, 13 to 20. (Verified by web search. The lemma there is for increasing events of independent binary coordinates; the form used in Theorem 2(e), for nondecreasing and nonincreasing functions of independent real coordinates, follows by induction on the number of coordinates from Chebyshev's association inequality for one coordinate.)
* Lai, T. L. (1998). Information bounds and quick detection of parameter changes in stochastic systems. *IEEE Transactions on Information Theory* 44(7), 2917 to 2929.
* Mei, Y. (2008). Is average run length to false alarm always an informative criterion? *Sequential Analysis* 27(4), 354 to 376. doi:10.1080/07474940802445790.
* Pergamenchtchikov, S. and Tartakovsky, A. G. (2018). Asymptotically optimal pointwise and minimax quickest change-point detection for dependent data. *Statistical Inference for Stochastic Processes* 21(1), 217 to 259. doi:10.1007/s11203-016-9149-x.
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science* 2(2), 229 to 260. doi:10.51387/23-NEJSDS51.
* Tartakovsky, A., Nikiforov, I. and Basseville, M. (2014). *Sequential Analysis: Hypothesis Testing and Changepoint Detection.* CRC Press. (Title, authors and year verified; the chapter on false alarm criteria was not consulted.)
* Ville, J. (1939), as in `docs/math.md`.
