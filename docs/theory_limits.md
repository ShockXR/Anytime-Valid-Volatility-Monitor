# Limits for every valid procedure, the plug-in baseline when n, m and A grow together, and the decrease bet at its boundary (TH5)

This note answers three questions left open by `docs/math.md`, `docs/theory_delay.md`, `docs/theory_running_tau.md` and `docs/theory_robustness.md` (TH4), and raised in sections 2, 3 and M3 of `docs/review/referee_report_paper.md`.

1. How fast can **any** procedure detect a change, if it must keep the false alarm guarantee of Theorem B over the nonparametric null? Sections 2 and 3 give a floor that holds path by path for every valid procedure, deterministic or randomised, and a delay floor built from it.
2. When does the **plug-in baseline** keep the nominal level, if the calibration length n, the horizon m and the alarm level A grow together? Section 4 proves a dichotomy at `m ≍ √n` for the worst stream in the null class, gives the margin z that the inflated plug-in needs, and relates it to the confidence-sequence baseline (Theorem 6.4 of `docs/math.md`, Theorem D of the paper).
3. Why is the growth-optimal decrease bet so close to its upper limit `λ_max = 1 / (B - 1)`, and what does a bet at `λ_max (1 - ε)` cost? Section 5.

Labels: **[PROVED]** with a complete proof here; **[PROVED, standard argument]** when the argument is known and is given for completeness; **[NUMERICAL]** for quantities computed by quadrature or Monte Carlo; **[HEURISTIC]** for arguments that are not proofs. All simulated numbers are **[SYNTHETIC]**: script `experiments/th5_limits.py` (seeds in its docstring), functions `src/avmonitor/limits.py`, tests `tests/test_limits.py`, outputs `results/th5_*.csv`.

## 0. Summary in plain words

* **A speed limit for every valid procedure (Theorem 2).** Every stopping rule with `P(T ≤ n) ≤ α_n` under every law of the null class (for Theorem B, `α_n = n / A`) can alarm at step n on a path y only if `W_n(y) = Π_{s ≤ n} 1 / π(y_s) ≥ 1 / α_n`, where `π(y)` is the largest probability that a null law can put on the value y: `π(y) = (B - 1) / (B - y)` for `y < 1` on the decrease side and `π(y) = 1 / y` for `y > 1` on the increase side. The proof uses one null law per path, an independent two-point law at each step, under which the path has probability exactly `1 / W_n(y)`. The floor is exact: for each path some valid rule alarms exactly at it. A randomised rule can alarm earlier, but only with probability at most `α_n W_n(y)`.
* **Complete collapse.** If every monitored return is zero (c = 3, Gaussian truncation, B = 9.045), no procedure valid at the levels `n / A` can alarm before step 45 at A = 7800 (50 at 2A). Shiryaev-Roberts e-detectors with linear bets need 58 (64) and the default grid 73 (79). The gap between 45 and 58 is the factor `n / B` in the level that a single path allows and an e-detector cannot use. Without truncation (`B = ∞`) no valid procedure alarms before step A on any path: the sequential form of the impossibility in `docs/math.md` section 5.3.
* **A delay floor (Theorem 5).** Combining the pathwise floor with a Neyman-Pearson bound at the least favourable null law gives a lower bound on the mean delay of every valid procedure. Over 18 post-change laws (Gaussian, Student t3, jumps; sd ratios 0.5 to 2), change at the first step, A = 7800: the default grid takes 1.69 to 2.16 times this bound, and the growth-optimal single bet 1.34 to 1.76 times. On all 180,000 simulated paths (10,000 per law), at both levels, neither the default grid nor the best single bet alarmed before the pathwise floor.
* **Asymptotics (Proposition 6).** The pathwise floor grows like `log A / γ_Q` with `γ_Q = E_Q log max(1, (B - Y) / (B - 1))` (decrease) and is weaker than Lorden's rate `log A / KL_inf` unless Q puts no mass beyond 1; the change-of-measure bound gives `log A / KL_inf` for every procedure valid at the levels `n / A`.
* **The plug-in baseline (Theorem 8, Proposition 9, Corollary 10).** With a calibration window of n returns, the drift of the plug-in baseline is at most `B̂ D_n`, where `D_n` is the Kolmogorov distance between the window and the return law. Hence, for every return law, `P(T ≤ m) ≤ (m / A)(1 + C (m + 1) / √n (1 + o(1)))` when `m / √n → 0`, with `C = 0.037` for Gaussian returns (first order) and `C = 2.52` for every law (distribution-free). Conversely, for the worst stream with the right conditional mean, `P(T ≤ m) / (m / A) → ∞` when `m / √n → ∞`, and `P(T ≤ m) → 1` when `m / (√n log A) → ∞`. The boundary is `m ≍ √n`, whatever A.
* **The inflated plug-in (Proposition 11).** Baselines `m̂ ± z se` keep `P(T ≤ m) ≤ m / A + 2 Φ̄(z) + o(1)` uniformly in m and A, and `2 Φ̄(z)` is attained by the worst stream when `m / (√n log A) → ∞`. So `z = Φ⁻¹(1 - δ / 2)` (2.576 for δ = 0.01) is the margin. The confidence-sequence baseline does the same without asymptotics, at an effective z of 3.8 to 5.7.
* **The decrease bet at its boundary (Proposition 12).** The growth of the bet `λ_max (1 - ε)` is `G_c + q_B log ε` up to small terms, where `q_B` is the probability that a return reaches the truncation level. For light tails `q_B` is tiny, so the optimum is extremely close to `λ_max` (`ε* = 2.4 × 10⁻⁸` for a Gaussian halving, as the referee found). For Student t3, `q_B` is 0.6% to 2.4% and the optimum is far inside (`λ* / λ_max = 0.55` at sd ratio 0.8); a bet at `λ_max (1 - 10⁻⁶)` then has negative growth and is 53% slower.

## 1. Set-up

As in `docs/math.md` sections 5.1 to 5.4 and `docs/theory_robustness.md` section 1.

* Monitored statistic `Y_t = min(r_t², τ²) / m₀ ∈ [0, B]`, `B = τ² / m₀ > 1`, with τ and m₀ known (sections 2, 3 and 5) or estimated from a calibration window (section 4). `F_t` is the information up to time t.
* **Null classes.** `H↓` (no decrease): all laws of `(Y_t)` on `[0, B]` with `E[Y_t | F_{t-1}] ≥ 1` for every t. `H↑` (no increase): `E[Y_t | F_{t-1}] ≤ 1`. `H=` (two-sided): `E[Y_t | F_{t-1}] = 1`. The laws in these classes need not be identically distributed or independent; independent, non-identical product laws are allowed, and the proofs below use them.
* **Validity.** A stopping rule T is *valid at levels* `α_1 ≤ α_2 ≤ ...` for a class H if `P(T ≤ n) ≤ α_n` for every `P ∈ H` and every n. Theorem B of `docs/math.md` gives `α_n = n / A` for one side at level A, and for the default two-sided rule (each side at 2A) over `H=`; each side alone at 2A is valid at `n / (2A)` over its own class.
* **Deterministic and randomised rules.** A deterministic rule is a stopping time of the filtration of `(Y_t)`: whether `T ≤ n` is a function of `y_{1:n}`. A randomised rule is `T(y, U)` with U uniform on (0, 1), independent of Y under every law considered, and `T(·, u)` a stopping time for each u (a sequence of independent uniforms can be generated from one).
* **Bets.** Increase factor `1 + λ (Y - 1)`, `0 ≤ λ ≤ 1`; decrease factor `1 + λ (1 - Y)`, `0 ≤ λ ≤ λ_max = 1 / (B - 1)`. Write `U = Y - 1` (increase) or `U = 1 - Y` (decrease). Default grid: `λ = 0.02, 0.05, 0.1, 0.2, 0.35, 0.5` and `λ = s / (B - 1)` with `s = 0.1, 0.25, 0.5, 0.75, 0.95`, equal weights.

## 2. A speed limit for every valid procedure

### 2.1 The largest probability of one value

**Lemma 1.** [PROVED] For `y ∈ [0, B]` let

`π↓(y) = (B - 1) / (B - y)` if `y < 1`, and `π↓(y) = 1` if `y ≥ 1`;  `π↑(y) = 1 / y` if `y > 1`, and `π↑(y) = 1` if `y ≤ 1`;  `π=(y) = π↓(y) π↑(y)`.

Then `π↓(y) = max {P(Y = y) : P a law on [0, B] with E_P Y ≥ 1}`, and likewise `π↑` for `E_P Y ≤ 1` and `π=` for `E_P Y = 1`. The maximum is attained by a two-point law with mean exactly 1: on `{y, B}` with mass `π↓(y)` at y (decrease side, y < 1), on `{0, y}` with mass `1 / y` at y (increase side, y > 1), and by the point mass at y when `π = 1` and the mean constraint allows it (for `π=`, `y = 1`).

*Proof.* Decrease side, y < 1: with `p = P(Y = y)`, `1 ≤ E Y ≤ p y + (1 - p) B` because `Y ≤ B`, so `p ≤ (B - 1) / (B - y)`. The law `p δ_y + (1 - p) δ_B` with `p = (B - 1) / (B - y)` has mean `B - p (B - y) = 1`. For y ≥ 1 the point mass `δ_y` has mean at least 1. Increase side, y > 1: `1 ≥ E Y ≥ p y` because `Y ≥ 0`, so `p ≤ 1 / y`, attained by `(1/y) δ_y + (1 - 1/y) δ_0`. For y ≤ 1, `δ_y` has mean at most 1. The two-sided case: for y < 1 the decrease-side law has mean exactly 1, and an upper bound on P(Y = y) over a smaller class is still `π↓(y)`; likewise for y > 1; for y = 1 the point mass. ∎

`tests/test_limits.py` checks Lemma 1 against a linear programme on a grid of `[0, B]`. For the default c = 3 and Gaussian truncation (B = 9.045), `π↓(0) = 0.8894`: a null law can put 89% of its mass on a zero return.

### 2.2 The floor lemma

For a path `y = (y_1, y_2, ...)` in `[0, B]` write `W_n(y) = Π_{s ≤ n} 1 / π(y_s)` with π the function of Lemma 1 for the class at hand. So `log W_n(y) = Σ_{s ≤ n} log max(1, (B - y_s) / (B - 1))` on the decrease side and `Σ_{s ≤ n} log max(1, y_s)` on the increase side; for `H=` it is the sum of both. Write `N_α(y) = min {n : W_n(y) ≥ 1 / α_n}` (`= ∞` if there is none), and `N_A = N_α` for `α_n = n / A`.

**Theorem 2 (floor lemma).** [PROVED] Let H be one of `H↓`, `H↑`, `H=` and π the matching function of Lemma 1.

(a) *Deterministic rules.* If T is valid at levels `α_n` for H, then for every path y and every n with `T(y) ≤ n`, `W_n(y) ≥ 1 / α_n`. In particular `T(y) ≥ N_α(y)` for every path.

(b) *Randomised rules.* If T is randomised and valid at levels `α_n`, then for every path y and every n, `P_U(T(y, U) ≤ n) ≤ α_n W_n(y)`.

(c) *Sharpness.* `sup_{P ∈ H} P(Y_{1:n} = y_{1:n}) = 1 / W_n(y)` for every y and n, even over adaptive (non-product) laws. If α is nondecreasing and `W_n(y) ≥ 1 / α_n`, the rule "stop at n if `Y_{1:n} = y_{1:n}`, otherwise never" is valid at levels α and stops at n on y. Hence `N_α(y) = min {T(y) : T deterministic and valid at levels α}` for every path y.

*Proof.* Fix a path y and n. For s ≤ n let `Q_s` be the two-point law of Lemma 1 for the value `y_s` (or `δ_{y_s}` when `π(y_s) = 1`), and for s > n let `Q_s = δ_1`. Let `P_y` be the product law under which `Y_1, Y_2, ...` are independent with laws `Q_1, Q_2, ...`. By independence `E_{P_y}[Y_s | F_{s-1}] = E_{Q_s} Y`, which is at least 1 (decrease), at most 1 (increase) or equal to 1 (two-sided) by Lemma 1. So `P_y ∈ H`. Under `P_y`, `P_y(Y_{1:n} = y_{1:n}) = Π_{s ≤ n} π(y_s) = 1 / W_n(y)`, which is positive: `π↓ ≥ (B - 1) / B` and `π↑ ≥ 1 / B`.

(a) Since T is a stopping time, the event `{T ≤ n}` is a set of values of `Y_{1:n}`. If `T(y) ≤ n`, it contains `y_{1:n}`, so `α_n ≥ P_y(T ≤ n) ≥ P_y(Y_{1:n} = y_{1:n}) = 1 / W_n(y)`. With `n = T(y)` this gives `T(y) ≥ N_α(y)`.

(b) By Fubini (U is independent of Y under `P_y`), `α_n ≥ P_y(T ≤ n) = ∫_0^1 P_y(T(·, u) ≤ n) du ≥ ∫_0^1 1{T(y, u) ≤ n} P_y(Y_{1:n} = y_{1:n}) du = P_U(T(y, U) ≤ n) / W_n(y)`, where the middle inequality uses (a)'s observation for each fixed u.

(c) "≥" is attained by `P_y`. For "≤", take any `P ∈ H` and write `E_s = {Y_{1:s} = y_{1:s}}`. If `P(E_{s-1}) > 0`, then `E[Y_s 1_{E_{s-1}}] = E[E[Y_s | F_{s-1}] 1_{E_{s-1}}]`, which is at least `P(E_{s-1})` on the decrease side; so the conditional law of `Y_s` given the event `E_{s-1}` (an elementary conditional probability) is a law on `[0, B]` with mean at least 1, and Lemma 1 gives `P(Y_s = y_s | E_{s-1}) ≤ π(y_s)`. The same holds on the other sides. Multiplying, `P(E_n) ≤ Π_{s ≤ n} π(y_s) = 1 / W_n(y)` (if some `P(E_{s-1}) = 0` then `P(E_n) = 0`). The singleton rule therefore has `P(T ≤ k) = 0` for k < n and `P(T ≤ k) = P(E_n) ≤ 1 / W_n(y) ≤ α_n ≤ α_k` for k ≥ n. ∎

**What the lemma needs, and what it does not.**

* *Positive probability.* The path has probability exactly `1 / W_n(y) ≥ B^{-n}` under the law `P_y` built for it. The lemma does not need the path to have positive probability under the post-change law Q; it is a statement about each path.
* *Continuous data and ties.* Nothing changes: each value `y_s`, repeated or not, gets its own two-point law. A path from a continuous law has probability 0 under that law but positive probability under `P_y`.
* *Non-identical laws.* The proof needs the class to contain independent, non-identically distributed laws. If validity is required only over i.i.d. null laws, the same argument gives `sup_P Π_s P({y_s})` over laws with the mean constraint, an empirical likelihood of the path, which is never larger than `1 / W_n(y)` and is equal to it on paths with a single repeated value, such as the complete collapse below.
* *Rules on raw returns.* A monitor sees `r_t`, not only `Y_t`. If the null is stated for the returns (`E[min(r_t², τ²) | F_{t-1}] ≥ m₀`), the two-point laws can be placed on `{r_s, τ}` (decrease) or `{r_s, 0}` (increase), which have the same masses. So the floor holds for every rule that uses the raw returns, and for every rule that also uses other information whose law the null leaves free.
* *Valid at one horizon only.* A rule that only promises `P(T ≤ m) ≤ α` for one m is valid at the constant levels `α_n = α` for n ≤ m. Its floor on a complete collapse is `log(1/α) / log(B / (B - 1))`: 26 steps for α = 0.05 at c = 3, against 45 for the levels `n / A` with A = 7800 (which equal 0.05 at one day). This is the price, on this path, of controlling every n at once.

### 2.3 Consequences

**Corollary 3 (no truncation, no detection of decreases).** [PROVED] Let `H↓_∞` be the class of all laws on `[0, ∞)` with `E[Y_t | F_{t-1}] ≥ 1`. A deterministic rule valid at levels `α_n` for `H↓_∞` never stops at a step n with `α_n < 1`, on any path. A randomised one stops by such an n with probability at most `α_n` on every path, which is what it could achieve by ignoring the data. With `α_n = n / A`: no deterministic valid rule alarms before step A on any path.

*Proof.* `H↓_∞` contains the class `H↓` built with any finite bound `B' > max_{s ≤ n} y_s`, so Theorem 2 applies with `π↓(y; B') = (B' - 1) / (B' - y)`, which tends to 1 as `B' → ∞`. Hence `W_n(y; B') → 1`, and (a) would need `1 ≥ 1 / α_n`; (b) gives `P_U(T ≤ n) ≤ α_n`. ∎

This is the sequential form, for every procedure, of the claim in `docs/math.md` section 5.3 (no e-variable for the untruncated decrease null exceeds 1), a one-sided relative of Bahadur and Savage (1956). The increase side has no such problem: `π↑(y) = 1 / y` does not depend on B.

**Corollary 4 (the most extreme paths).** [PROVED] (a) On a *complete collapse*, `y_s = 0` for every s, every rule valid at levels `n / A` over `H↓` (or over `H=`) satisfies `T ≥ n_↓ = min {n : n log(B / (B - 1)) ≥ log(A / n)}`. (b) On a *complete explosion*, `y_s = B` for every s, every rule valid at `n / A` over `H↑` (or `H=`) satisfies `T ≥ n_↑ = min {n : n log B ≥ log(A / n)}`. Both floors are attained by valid rules (Theorem 2c).

For comparison, `docs/theory_delay.md` (Corollary 4) shows that every Shiryaev-Roberts e-detector whose factors are valid e-variables satisfies `M_n ≤ B ((B / (B - 1))^n - 1)` on the collapse, so it cannot alarm before `log(1 + A / B) / log(B / (B - 1))` steps, and the default grid needs more. Gaussian truncation, `results/th5_collapse.csv`:

| c | B | collapse: every valid rule, A = 7800 (2A) | fractional root | SR e-detector, linear bets | default grid | `log A / log(B / (B - 1))` |
|---|---|---|---|---|---|---|
| 1.5 | 2.890 | 15 (17) | 14.77 (16.18) | 18.61 (20.24) | 23 (25) | 21.1 |
| 2 | 4.345 | 23 (25) | 22.38 (24.66) | 28.65 (31.30) | 36 (39) | 34.3 |
| 2.5 | 6.393 | 33 (36) | 32.26 (35.74) | 41.79 (45.86) | 52 (57) | 52.7 |
| 3 | 9.045 | 45 (50) | 44.16 (49.16) | 57.71 (63.62) | 73 (79) | 76.5 |
| 4 | 16.00 | 73 (82) | 72.50 (81.44) | 95.94 (106.67) | 122 (134) | 138.9 |
| 5 | 25.00 | 106 (120) | 105.43 (119.37) | 140.76 (157.70) | 180 (199) | 219.5 |

| c | explosion: every valid rule, A (2A) | SR e-detector, linear bets | default grid (largest bet 0.5) |
|---|---|---|---|
| 1.5 | 7 (8) | 9 (9) | 15 (16) |
| 2 | 6 (6) | 6 (7) | 11 (12) |
| 2.5 | 5 (5) | 5 (6) | 8 (9) |
| 3 | 4 (4) | 5 (5) | 7 (7) |
| 4 | 3 (4) | 4 (4) | 5 (6) |
| 5 | 3 (3) | 3 (3) | 5 (5) |

**How far is the λ_max e-detector from the floor, and why.** On the collapse at c = 3 the single bet at `λ_max` (the SR column) needs 57.7 steps against 44.2 for the floor. The Shiryaev-Roberts statistic sums the products from every start time, but on this path the sum is dominated by the longest product and equals `B (F^n - 1)` with `F = B / (B - 1)`; it reaches A when `F^n ≈ A / B`. The floor needs `F^n ≥ A / n`. The difference is `log(n / B) / log F ≈ log(45 / 9.045) / 0.1172 = 13.7` steps. In words: a rule valid at the levels `n / A` may, on one particular path, spend the whole budget `n / A` accumulated by step n, while an e-detector has to spread its expected value `E M_n ≤ n` over all paths. The default grid adds 15 more steps (73) because four of its five bets are smaller than the largest. Corollary 4 does not say that some practical rule detects collapses after 45 steps from every starting point: a rule that tests disjoint blocks of 45 steps for all-zero returns is valid (`(A / 45) × 0.8894^45 ≤ 1`, so the union bound gives `P(T ≤ n) ≤ n / A`), but a collapse that starts inside a block is found only at the end of the next block, up to 89 steps later.

**The two-sided class.** For the default monitor the two-sided class `H=` is the relevant one, since its guarantee `P(T ≤ m) ≤ m / A` holds when the conditional mean is exactly m₀. Every rule valid at `n / A` over `H=` needs 45 steps on the collapse at c = 3; the default two-sided monitor, whose decrease side alarms at 2A, needs 79 (a factor 1.76).

### 2.4 Is the lemma known?

We found no statement of Theorem 2 in this form, for sequential rules valid over a nonparametric mean-constrained class with a level that grows with n, and we make no claim of priority: the argument is elementary. Its ingredients are standard. Lemma 1 is the extreme-point argument that also gives the speed limit of `docs/theory_delay.md` (Corollary 4) and the claim of `docs/math.md` 5.3. In the game-theoretic framework of Shafer and Vovk (2019), `1 / W_n(y)` is the largest probability that the null class assigns to the path, and the bound `T(y) ≥ N_α(y)` says that a valid rule cannot reject a path whose largest null probability exceeds the level. For a single bounded mean, Clerico (2025) and Waudby-Smith, Sandoval and Jordan (2025) prove optimality properties of betting procedures among all valid ones (minimal confidence-set width; growth and expected rejection time); Theorem 2 is a pathwise statement in the same spirit for the change detection problem. The searches we ran (search-engine listings only, since arxiv.org could not be reached from this machine) did not show a pathwise floor of this kind for change detection.

## 3. A delay floor for every valid procedure

### 3.1 Theorem 5

**Theorem 5 (delay floor).** [PROVED] Let T be valid at levels `n / A` for H (one of `H↓`, `H↑`, `H=`), let Q be any law of the monitored stream (the post-change law, change at the first step), and let U be independent of the stream under Q.

(a) If T is deterministic, `T ≥ N_A` Q-almost surely, so `Q(T ≤ n) ≤ Q(N_A ≤ n)` for every n and `E_Q[T] ≥ E_Q[N_A]`.

(b) If T is randomised, `Q(T ≤ n) ≤ E_Q[min{1, (n / A) W_n}]`.

(c) Let `P° ∈ H` be a law with `Q|_{F_n} ≪ P°|_{F_n}` and `Λ_n = log dQ/dP°` on `F_n`. Then `Q(T ≤ n) ≤ β_n`, where

`β_n = sup {E_Q φ : φ a function of y_{1:n}, 0 ≤ φ ≤ c_n, E_{P°} φ ≤ n / A}`,

with `c_n = 1` for every rule, `c_n = 1{N_A ≤ n}` for deterministic rules and `c_n = min{1, (n / A) W_n}` for randomised rules. The supremum is attained by the test that fills `c_n` in decreasing order of `Λ_n` until the `P°`-budget `n / A` is spent. Consequently `E_Q[T] ≥ 1 + Σ_{n ≥ 1} (1 - β_n)` and every p-quantile of T is at least `min {n : β_n ≥ p}`.

(d) If Q is i.i.d. and the growth-optimal bet `λ*` of Q is interior, `0 < λ* < λ_max`, then `P* = Q / (1 + λ* U)`, taken i.i.d., is in `H↓ ∩ H↑ ∩ H=` (mean of Y exactly 1) and `Λ_n = Σ_{s ≤ n} log(1 + λ* U_s)`, with `E_Q Λ_1 = KL_inf(Q)`.

Since Lorden's worst-case delay is a supremum over change times that includes the change at the first step, (a) to (c) are lower bounds on `WADD(T)` for every valid rule.

*Proof.* (a) and (b) are Theorem 2 (a) and (b), integrated over Q. (c) For a deterministic rule, `φ = 1{T ≤ n}` is a function of `y_{1:n}` with `E_{P°} φ = P°(T ≤ n) ≤ n / A` and, by (a), `φ ≤ 1{N_A ≤ n}`. For a randomised rule take `φ(y) = P_U(T(y, U) ≤ n)`, which satisfies `φ ≤ min{1, (n / A) W_n}` by Theorem 2b. In both cases `Q(T ≤ n) = E_Q φ ≤ β_n`. Optimality of the greedy test `φ*` (the constrained Neyman-Pearson lemma): let k be its threshold, so `φ* = c_n` where `Λ_n > k` and `φ* = 0` where `Λ_n < k`. For any feasible φ, `(φ* - φ)(e^{Λ_n} - e^k) ≥ 0` pointwise; integrating against `P°` gives `E_Q φ* - E_Q φ ≥ e^k (E_{P°} φ* - E_{P°} φ) ≥ 0` when `φ*` spends the whole budget (if it does not, `φ* = c_n` everywhere and is trivially optimal). The delay follows from `E_Q T = Σ_{n ≥ 0} Q(T > n)`. (d) is Lemma 2b of `docs/theory_delay.md`: at an interior maximiser `E_Q[U / (1 + λ* U)] = 0`, so `P*` has total mass `E_Q[1 / (1 + λ* U)] = 1 - λ* E_Q[U / (1 + λ* U)] = 1` and `E_{P*} U = 0`. ∎

**Numerical evaluation.** [NUMERICAL] `β_n` is estimated from N = 10,000 paths drawn from Q: the `P°`-expectations are importance-weighted averages with weights `e^{-Λ_n} / N` (`limits.np_power`, `limits.np_cdf`), so the bounds below are Monte Carlo estimates of the bounds, not certified values. `λ*` is the root of `E_Q[U / (1 + λU)] = 0` found in `log ε` with `λ = λ_max (1 - ε)` (`limits.least_favourable_bet`); for the Gaussian halving `ε* = 2.4 × 10⁻⁸`, and a maximiser found to a tolerance in λ (as `design.best_bet` does) leaves `P*` with mass 0.982 instead of 1, which is not a law in the null class. (The missing mass sits at `Y = B`, which no simulated value of that law reached, so the estimates did not change when this was corrected; see `docs/iterations/TH5.md`.) The growth-optimal bets of all 18 laws are interior: every decrease law has an atom at B, and every increase law has `E_Q[1 / Y] = ∞`.

### 3.2 Asymptotics

Write `ξ(y) = log(1 / π(y))` and `γ_Q = E_Q ξ(Y)`: `γ_Q = E_Q log max(1, (B - Y) / (B - 1))` on the decrease side and `E_Q log max(1, Y)` on the increase side.

**Proposition 6.** [PROVED] Let Q be i.i.d. with `γ_Q > 0`.

(a) `N_A / log A → 1 / γ_Q` almost surely and in L¹ as `A → ∞`.

(b) `γ_Q ≥ KL_inf(Q)`, with equality exactly when `Q(Y > 1) = 0` (decrease side) or `Q(Y < 1) = 0` (increase side).

(c) If `λ*` is interior, then for every rule T valid at levels `n / A` over a class containing `P*` (deterministic or randomised) and every `ε ∈ (0, 1)`, `Q(T ≤ (1 - ε) log A / KL_inf) → 0`, so `E_Q[T] ≥ (log A / KL_inf(Q))(1 - o(1))`.

*Proof.* (a) `log W_n = Σ_{s ≤ n} ξ(Y_s)` with i.i.d. increments in `[0, ξ_max]`, `ξ_max = log(B / (B - 1))` or `log B`. Since `log W_n ≤ n ξ_max`, `N_A → ∞`. By the definition of `N_A`, `log W_{N_A} ≥ log A - log N_A` and `log W_{N_A - 1} < log A - log(N_A - 1)`. Divide by `N_A` and use the strong law (`log W_n / n → γ_Q`) and `log N_A / N_A → 0`: `log A / N_A → γ_Q`. For L¹: `N_A ≤ τ_A = min {n : log W_n ≥ log A}`, and `τ_A / log A → 1 / γ_Q` in L¹ by the elementary renewal theorem, so dominated convergence (with a dominating sequence that converges in L¹) applies.

(b) Every bet factor `f_λ` is an e-variable for every law in the class, in particular for the two-point law of Lemma 1 at y: `π(y) f_λ(y) ≤ E f_λ ≤ 1`, so `log f_λ(y) ≤ ξ(y)` for every y (for y with `π(y) = 1` the point mass gives `f_λ(y) ≤ 1`). Hence `g(λ) = E_Q log f_λ ≤ γ_Q` for every λ, and `KL_inf = max_λ g(λ) ≤ γ_Q` (Lemma 2 of `docs/theory_delay.md`). Equality forces `log f_{λ*}(Y) = ξ(Y)` Q-almost surely. On the decrease side, on `{Y > 1}` the right side is 0 and the left side is negative when `λ* > 0`, so `Q(Y > 1) = 0`. Conversely, if `Q(Y > 1) = 0` then `U ≥ 0`, g is nondecreasing, `λ* = λ_max` and `f_{λ_max}(y) = (B - y) / (B - 1) = 1 / π↓(y)` for `y ≤ 1`, so `g(λ_max) = γ_Q`. The increase side is the same with `f_1(y) = y`.

(c) Let `φ = P_U(T ≤ n | y_{1:n})`. For any a, `Q(T ≤ n) = E_Q[φ; Λ_n ≤ a] + E_Q[φ; Λ_n > a] ≤ e^a E_{P*}[φ] + Q(Λ_n > a) ≤ e^a n / A + Q(Λ_n > a)`. Take `n = n_A = ⌊(1 - ε) log A / KL_inf⌋` and `a = (1 - ε/2) log A`. The first term is `n_A A^{-ε/2} → 0`; the second tends to 0 by the weak law of large numbers, because `a / n_A → KL_inf (1 - ε/2) / (1 - ε) > KL_inf = E_Q Λ_1`. So `E_Q T ≥ n_A Q(T > n_A) = (1 - ε)(log A / KL_inf)(1 - o(1))` for every ε. ∎

This is the change-of-measure argument behind the lower bounds of Lorden (1971) and Lai (1998). Theorem 2b of `docs/theory_delay.md` states it for rules whose average run length is at least A under the i.i.d. null laws; (c) is the same rate for the validity notion of Theorem B, for randomised rules, and needs only one null law, `P*`. Lai (1998) introduced a false alarm constraint of this kind, on the probability of an alarm within a window (as described in secondary sources; we could not open the paper). By (b), the pathwise floor alone gives the slower rate `log A / γ_Q`, and gives Lorden's rate exactly when the post-change law has no mass on the far side of 1, as on the complete collapse.

The asymptotics have a second-order term that matters at A = 7800: because `α_n = n / A`, the floor is close to the root `n*` of `n γ_Q = log(A / n)`, not to `log A / γ_Q`. Part `asymptotic` of the script (`results/th5_asymptotic.csv`, 4000 paths per law) shows this. [NUMERICAL]

| Law | sd ratio | A | floor `E N_A` | floor / `n*(γ_Q)` | floor / `(log A / γ_Q)` | NP, deterministic / `n*(KL_inf)` | best single bet / `n*(KL_inf)` |
|---|---|---|---|---|---|---|---|
| gaussian | 0.5 | 1e+03 | 37.1 | 1.016 | 0.487 | 0.993 | 1.367 |
| gaussian | 0.5 | 1e+04 | 57.5 | 1.010 | 0.567 | 0.985 | 1.319 |
| gaussian | 0.5 | 1e+06 | 101.9 | 1.006 | 0.669 | 0.979 | 1.242 |
| gaussian | 0.5 | 1e+08 | 148.4 | 1.004 | 0.731 | 0.976 | 1.194 |
| gaussian | 0.8 | 1e+03 | 46.1 | 1.015 | 0.454 | 0.731 | 1.334 |
| gaussian | 0.8 | 1e+04 | 73.1 | 1.010 | 0.540 | 0.758 | 1.308 |
| gaussian | 0.8 | 1e+06 | 132.1 | 1.006 | 0.651 | 0.815 | 1.241 |
| gaussian | 0.8 | 1e+08 | 194.1 | 1.004 | 0.717 | 0.852 | 1.192 |
| gaussian | 1.5 | 1e+03 | 9.6 | 1.145 | 0.792 | 0.766 | 1.284 |
| gaussian | 1.5 | 1e+04 | 13.1 | 1.106 | 0.809 | 0.800 | 1.265 |
| gaussian | 1.5 | 1e+06 | 20.5 | 1.070 | 0.842 | 0.855 | 1.215 |
| gaussian | 1.5 | 1e+08 | 28.1 | 1.051 | 0.863 | 0.885 | 1.180 |
| student t3 | 0.5 | 1e+03 | 27.4 | 1.021 | 0.535 | 0.823 | 1.304 |
| student t3 | 0.5 | 1e+04 | 41.3 | 1.014 | 0.606 | 0.828 | 1.279 |
| student t3 | 0.5 | 1e+06 | 71.3 | 1.008 | 0.697 | 0.851 | 1.220 |
| student t3 | 0.5 | 1e+08 | 102.7 | 1.006 | 0.753 | 0.871 | 1.177 |
| student t3 | 0.8 | 1e+03 | 32.4 | 1.020 | 0.509 | 0.627 | 1.213 |
| student t3 | 0.8 | 1e+04 | 49.6 | 1.013 | 0.585 | 0.704 | 1.244 |
| student t3 | 0.8 | 1e+06 | 86.8 | 1.006 | 0.682 | 0.793 | 1.212 |
| student t3 | 0.8 | 1e+08 | 125.8 | 1.004 | 0.741 | 0.845 | 1.178 |
| student t3 | 1.5 | 1e+03 | 11.0 | 1.132 | 0.759 | 0.719 | 1.260 |
| student t3 | 1.5 | 1e+04 | 15.1 | 1.092 | 0.780 | 0.768 | 1.238 |
| student t3 | 1.5 | 1e+06 | 24.0 | 1.063 | 0.823 | 0.830 | 1.189 |
| student t3 | 1.5 | 1e+08 | 33.0 | 1.045 | 0.849 | 0.862 | 1.160 |

`n*(r)` is the root of `n r = log(A / n)`. The floor is within 0.4% to 5% of `n*(γ_Q)` at `A = 10⁸` but only 0.73 to 0.86 of `log A / γ_Q`: the `log n` correction is large at every practical A. The deterministic Neyman-Pearson bound approaches `n*(KL_inf)` from below (0.85 to 0.98 of it at `A = 10⁸`), and the best single bet from above (1.16 to 1.19). The ratio of the best single bet to the bound falls from 1.38 to 1.93 at `A = 10³` to 1.22 to 1.40 at `A = 10⁸`. Four of the 4000 paths of Student t3 at sd ratio 0.8 had no alarm of the best bet within 2112 steps at `A = 10⁸`; their delay is counted as 2112, so that mean is slightly too small.

### 3.3 The floor against the e-detectors

[NUMERICAL, SYNTHETIC] Part `floors` of the script, `results/th5_floors.csv` and `results/th5_floor_cdf.csv`: 10,000 paths of 4000 steps per law, change at the first step, true baseline, c = 3, A = 7800. "Floor" is `E_Q[N_A]` (Theorem 5a); "NP" are the bounds of Theorem 5c at `P*`, for every rule, for deterministic rules (combined with the floor) and for randomised rules; "grid" and "best bet" are the simulated mean delays of the default grid and of the single bet `λ*` at threshold A, on the same paths. The last two columns divide them by the deterministic bound `max(floor, NP deterministic)`. No path was left without an alarm or a floor crossing.

| Law | sd ratio | side | `KL_inf` | `γ_Q` | floor | NP, every rule | NP, deterministic | NP, randomised | default grid | best single bet | grid / bound | best bet / bound |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gaussian | 0.5 | down | 0.0881 | 0.0908 | 55.2 | 45.6 | 55.2 | 47.0 | 93.2 | 74.2 | 1.69 | 1.34 |
| gaussian | 0.67 | down | 0.0626 | 0.0769 | 63.3 | 57.2 | 64.7 | 58.8 | 124.2 | 98.7 | 1.92 | 1.52 |
| gaussian | 0.8 | down | 0.0346 | 0.0680 | 70.1 | 85.5 | 90.4 | 86.9 | 194.9 | 156.5 | 2.16 | 1.73 |
| gaussian | 1.25 | up | 0.0514 | 0.4077 | 16.8 | 65.3 | 65.4 | 65.4 | 120.7 | 108.3 | 1.85 | 1.66 |
| gaussian | 1.5 | up | 0.1645 | 0.5669 | 12.7 | 26.1 | 26.2 | 26.2 | 47.8 | 41.7 | 1.82 | 1.59 |
| gaussian | 2 | up | 0.4302 | 0.8422 | 9.1 | 12.0 | 12.2 | 12.0 | 23.0 | 18.4 | 1.88 | 1.51 |
| student t3 | 0.5 | down | 0.0913 | 0.1350 | 39.7 | 40.6 | 44.9 | 41.6 | 81.5 | 69.3 | 1.82 | 1.55 |
| student t3 | 0.67 | down | 0.0450 | 0.1193 | 44.1 | 70.4 | 72.5 | 70.9 | 138.4 | 123.3 | 1.91 | 1.70 |
| student t3 | 0.8 | down | 0.0172 | 0.1085 | 47.8 | 146.2 | 147.0 | 146.3 | 281.3 | 259.3 | 1.91 | 1.76 |
| student t3 | 1.25 | up | 0.0243 | 0.3623 | 18.7 | 117.3 | 117.3 | 117.3 | 218.4 | 196.5 | 1.86 | 1.67 |
| student t3 | 1.5 | up | 0.0850 | 0.4745 | 14.8 | 44.4 | 44.4 | 44.4 | 80.3 | 71.9 | 1.81 | 1.62 |
| student t3 | 2 | up | 0.2521 | 0.6742 | 10.9 | 18.5 | 18.6 | 18.5 | 34.3 | 28.9 | 1.84 | 1.55 |
| jumps | 0.5 | down | 0.0752 | 0.0931 | 54.1 | 49.1 | 56.9 | 50.2 | 99.0 | 82.6 | 1.74 | 1.45 |
| jumps | 0.67 | down | 0.0515 | 0.0792 | 61.8 | 63.8 | 69.9 | 65.1 | 134.3 | 112.6 | 1.92 | 1.61 |
| jumps | 0.8 | down | 0.0273 | 0.0701 | 68.4 | 102.6 | 106.2 | 103.5 | 215.5 | 185.3 | 2.03 | 1.74 |
| jumps | 1.25 | up | 0.0480 | 0.4016 | 17.0 | 68.8 | 68.9 | 68.8 | 127.3 | 114.2 | 1.85 | 1.66 |
| jumps | 1.5 | up | 0.1568 | 0.5578 | 12.9 | 27.0 | 27.1 | 27.1 | 49.2 | 43.1 | 1.81 | 1.59 |
| jumps | 2 | up | 0.4164 | 0.8291 | 9.2 | 12.2 | 12.5 | 12.3 | 23.5 | 18.9 | 1.89 | 1.51 |

* **No violation.** On all 180,000 paths (18 laws, 10,000 paths each), at both levels A and 2A, the alarms of the default grid and of the best single bet came at or after the pathwise floor, as Theorem 2 requires (`pathwise_violations_*` columns).
* **Efficiency.** At A the default grid takes 1.69 to 2.16 times the lower bound for deterministic valid rules, and the best single bet 1.34 to 1.76 times; at 2A (each side's level in the two-sided monitor) 1.65 to 2.09 and 1.33 to 1.72. The ratio is closest to 1 for the strongest decreases, where the pathwise floor is the binding bound, and largest for small decreases (Gaussian sd ratio 0.8: 2.16).
* **Which bound binds.** For strong decreases with light tails the pathwise floor is the larger bound (Gaussian 0.5: 55.2 against 45.6 for the Neyman-Pearson bound over every rule), and it sets the bound for deterministic rules. For increases and for fat tails the Neyman-Pearson bound is far above the floor (Student t3 at 0.8: floor 47.8, NP 146), because values of Y beyond 1 cost nothing in `W_n` (they cannot be ruled out by any single null law) but carry information about the change.
* **Randomisation buys little.** The randomised bounds are at most 15% below the deterministic ones (largest gap: Gaussian 0.5, 47.0 against 55.2) and less than 1% below for every increase.
* **Consistency with E11.** The default grid's delays agree with E11a's change-at-the-first-step delays (`results/e11a_delay.csv`, `sim_nu0_grid_mean`): for example 93.2 against 93.2 for the Gaussian halving and 120.7 against 122.1 for the Gaussian sd ratio 1.25. E4's delays (`e4_mean`) are smaller because E4 changed at step 501 after a pre-change history; Theorem 5 does not apply to them, since pre-change data can already have raised `W_n`.
* **What the gap is made of.** The ratio of the best single bet to the bound falls slowly with A (1.38 to 1.93 at `A = 10³`, 1.22 to 1.40 at `A = 10⁸`, table of section 3.2). The bound treats each n separately and may spend the whole budget `n / A` at each of them, while a detector must control all n at once; how much of the remaining gap is this price, and how much is the e-detector's own loss, is not settled here.

## 4. The plug-in baseline when n, m and A grow together

### 4.1 Set-up

The monitor's plug-in (`docs/theory_running_tau.md`, set-up of Proposition A.6; `docs/theory_robustness.md`, section 4): a calibration window `W = (r_1, ..., r_n)`, `τ̂ = a × median_i |r_i|` with `a = c × 1.4826`, `m̂ = n⁻¹ Σ_i min(r_i², τ̂²)`, `B̂ = τ̂² / m̂`, decrease bets `s_l / (B̂ - 1)`, both sides normalised by m̂, each side alarming at 2A.

**(S)** The calibration returns are i.i.d. with a law F with `P(r = 0) = 0`. Write `m_F(s) = E_F[min(r², s²)]`. The monitored stream belongs to the class `C(W)` of streams with `E[min(r_t², τ̂²) | F_{t-1}] = m_F(τ̂)` for every t, where `F_0` contains W. The stream may depend on W. Monitored returns that are i.i.d. F and independent of W are one member of `C(W)`.

The **drift** of the window is `d(W) = m_F(τ̂) / m̂ - 1`: under any stream in `C(W)`, `E[Y_t | F_{t-1}] = 1 + d(W)`. TH4's Corollary 1 gives, for every stream in `C(W)`,

`P(T ≤ m) ≤ E_W[min{1, (m / A) I_m(d⁺, d⁻; B̂)}]`, `I_m(δ↑, δ↓; B̂) = ½ [Σ_k w_k G_m(λ_k δ↑) + Σ_l v_l G_m(s_l δ↓ / (B̂ - 1))]`, `G_m(x) = ((1 + x)^{m+1} - (1 + x)) / (m x)`.

Write `λ̄ = 0.2033` and `s̄ = 0.51` for the mean bets, `λ_max = 0.5`, `s_max = 0.95`, and `D_n = sup_x |F̂_n(x) - F(x)|` for the Kolmogorov distance between the empirical law of the `|r_i|` and the law of `|r|` (the same as for `r²`).

### 4.2 Two lemmas

**Lemma 7.** [PROVED]

(a) For every window, `|d(W)| ≤ B̂(W) D_n`.

(b) If the median of the `|r_i|` is positive, `b_lo = 2a² / (1 + a²) ≤ B̂ ≤ b_hi = 2a²` (1.904 and 39.57 at c = 3).

(c) For `x ≥ 0` and `m ≥ 1`, `G_m(x) ≤ 1 + x (m + 1) e^{m x} / 2`. Hence `I_m(δ↑, δ↓; B̂) ≤ 1 + ((m + 1) / 4)(λ̄ δ↑ + φ̂ δ↓) e^{m x*}` with `φ̂ = s̄ / (B̂ - 1)` and `x* = max(λ_max δ↑, s_max δ↓ / (B̂ - 1))`, when at most one of `δ↑, δ↓` is positive.

*Proof.* (a) For a nonnegative variable Z and `t ≥ 0`, `E min(Z, t) = ∫_0^t P(Z > x) dx`. Applied to the empirical law and to F of `Z = r²`, for every `s > 0`, `m̂(s) - m_F(s) = ∫_0^{s²} (F^{(2)}(x) - F̂^{(2)}_n(x)) dx`, where `F^{(2)}` is the distribution function of `r²`. So `|m̂(s) - m_F(s)| ≤ s² D_n` for all s at once, and at `s = τ̂` (where `m̂(τ̂) = m̂`), `|m̂ - m_F(τ̂)| ≤ τ̂² D_n`. Divide by m̂. (b) is in `docs/theory_robustness.md` section 4: at least half of the `|r_i|` are at least the median, so `m̂ ≥ median² / 2`, and at least half are at most the median, so `m̂ ≤ (median² + τ̂²) / 2`. (c) `(1 + x)^j - 1 ≤ e^{jx} - 1 ≤ j x e^{jx}` (because `e^u - 1 ≤ u e^u` for `u ≥ 0`), so `G_m(x) - 1 = m⁻¹ Σ_{j=1}^{m} ((1 + x)^j - 1) ≤ (x / m) e^{mx} Σ_{j ≤ m} j = x (m + 1) e^{mx} / 2`. Average over the bets with their weights. ∎

`tests/test_limits.py` checks (a) on simulated Gaussian and t3 windows (the largest ratio `|d| / (B̂ D_n)` over 20,000 windows per cell is reported in `results/th5_plugin_drift.csv`), (b) on Cauchy, exponential, two-valued and t1.5 samples, and (c) on a grid.

### 4.3 Upper bounds

**Theorem 8.** Assume (S). For every m ≥ 1 and A ≥ m:

(a) [PROVED] `P(T ≤ m) ≤ (m / A) E_W[min{A / m, 1 + ((m + 1) / 4)(λ̄ d⁺ + φ̂ d⁻) e^{m x*(W)}}]`, with `x*(W) = max(λ_max d⁺, s_max d⁻ / (B̂ - 1))`.

(b) [PROVED, distribution-free] With `k1 = max(λ̄ b_hi, s̄ b_lo / (b_lo - 1))` and `k2 = max(λ_max b_hi, s_max b_lo / (b_lo - 1))` (8.045 and 19.78 at c = 3),

`P(T ≤ m) ≤ (m / A)(1 + ((m + 1) / 4) k1 M_n(k2 m))`, `M_n(t) = ∫_0^1 (1 + t x) e^{t x} min{1, 2 e^{-2 n x²}} dx`,

for every F with `P(r = 0) = 0`. In particular, uniformly over all such F, `P(T ≤ m) ≤ (m / A)(1 + C_DF (m + 1) / √n (1 + o(1)))` as `m / √n → 0`, with `C_DF = k1 √(π / 2) / 4 = 2.52`.

(c) [PROVED] If, in addition, F is fixed and the median of `|r|` is unique, then `√n d(W) ⇒ N(0, σ_Y²)` with `σ_Y² = Var_F(min(r², τ²)) / m₀²` (τ the population level, `m₀ = m_F(τ)`), and for every sequence `m = m_n` with `m / √n → 0`,

`P(T ≤ m) ≤ (m / A)(1 + C_F (m + 1) / √n (1 + o(1)))`, `C_F = (λ̄ + φ̄) σ_Y / (4 √(2π))`, `φ̄ = s̄ / (B - 1)`,

uniformly in `A ≥ m`. For Gaussian returns at c = 3, `σ_Y = 1.38585`, B = 9.045 and `C_F = 0.0369`; for Student t3, `σ_Y = 1.58267`, B = 6.362 and `C_F = 0.0471`.

*Proof.* (a) TH4's Corollary 1 and Lemma 7c.

(b) By Lemma 7a and 7b, `λ̄ d⁺ ≤ λ̄ B̂ D_n ≤ λ̄ b_hi D_n`, and `φ̂ d⁻ ≤ s̄ (B̂ / (B̂ - 1)) D_n ≤ s̄ (b_lo / (b_lo - 1)) D_n` because `B̂ / (B̂ - 1)` decreases in B̂. So `λ̄ d⁺ + φ̂ d⁻ ≤ k1 D_n` and likewise `x* ≤ k2 D_n`. Drop the cap in (a): `P(T ≤ m) ≤ (m / A)(1 + ((m + 1) / 4) k1 E[D_n e^{k2 m D_n}])`. For `g(x) = x e^{tx}`, `E g(D_n) = ∫_0^1 g'(x) P(D_n > x) dx` (as `D_n ≤ 1`), and the Dvoretzky-Kiefer-Wolfowitz inequality with Massart's constant, `P(D_n > x) ≤ 2 e^{-2 n x²}` for every x (Dvoretzky, Kiefer and Wolfowitz 1956; Massart 1990), gives `M_n(t)`. (Massart's inequality is stated for a continuous distribution function; for a law of `|r|` with atoms, the quantile transform writes `D_n` as a supremum of the uniform empirical process over a subset of [0, 1], which is at most the supremum over all of [0, 1], so the same tail bound holds.) For the rate, `M_n(t) = n^{-1/2} ∫_0^{√n} (1 + (t / √n) u) e^{(t / √n) u} min{1, 2 e^{-2u²}} du`, and with `t / √n = k2 m / √n → 0` the integral tends to `∫_0^∞ min{1, 2 e^{-2u²}} du ≤ √(π / 2)`.

(c) *The drift is asymptotically normal.* Write `m̂(s) - m_F(s) = Δ(s)`. Then `m̂ - m_F(τ̂) = Δ(τ) + (Δ(τ̂) - Δ(τ))` and, as in Lemma 7a, `|Δ(τ̂) - Δ(τ)| = |∫_{τ²}^{τ̂²} (F^{(2)} - F̂^{(2)}_n) dx| ≤ |τ̂² - τ²| D_n`. The sample median of `|r_i|` is consistent when the population median is unique, so `τ̂ → τ` in probability, and `√n D_n` is bounded in probability (DKW), so `√n |Δ(τ̂) - Δ(τ)| → 0` in probability. By the central limit theorem `√n Δ(τ) ⇒ N(0, Var_F min(r², τ²))`, and `m̂ → m₀` in probability. Hence `√n d = -√n (m̂ - m_F(τ̂)) / m̂ ⇒ N(0, σ_Y²)`. The estimation of τ does not enter at first order, because the null is the truncated mean at the estimated level.
*The bound.* By (a), `P(T ≤ m) ≤ (m / A)(1 + ((m + 1) / (4 √n)) E[V_n])` with `V_n = √n (λ̄ d⁺ + φ̂ d⁻) e^{m x*}`. Here `φ̂ → φ̄` (since `B̂ → B`), `m x* ≤ k2 (m / √n) √n D_n → 0` in probability, and `√n d ⇒ N(0, σ_Y²)`, so `V_n ⇒ (λ̄ + φ̄) σ_Y × (Z⁺ or Z⁻)`, a variable with mean `(λ̄ + φ̄) σ_Y E[Z⁺] = (λ̄ + φ̄) σ_Y / √(2π)` (Z standard normal; `d⁺` and `d⁻` are never both positive). The family `V_n` is uniformly integrable: `V_n ≤ k1 √n D_n e^{k2 (m / √n) √n D_n}`, and DKW bounds `E[(√n D_n)² e^{2 k2 (m / √n) √n D_n}]` uniformly in n once `m / √n ≤ 1`. So `E V_n → (λ̄ + φ̄) σ_Y / √(2π)`, which is the claim. A enters only through the cap, which can only lower the bound. ∎

**Reading.** The upper bound is close to the nominal `m / A` exactly when the horizon is short compared with √n. For Gaussian returns the first-order excess is 3.7% of `m / A` per unit of `(m + 1) / √n`; the distribution-free constant is 68 times larger, because it must allow for any law and uses `D_n` rather than the drift. (b) is informative only for small m: for n = 10⁶ it gives a factor 1.15 at m = 39 and 2.3 at m = 98, for n = 10⁸ a factor 1.013 at m = 98 (`results/th5_distribution_free.csv`).

### 4.4 Lower bound: the worst stream in the class

**Proposition 9 (hold, then bet boldly).** [PROVED] Assume (S), fix a window with `d = d(W) ≠ 0`, and take the side that d favours: the increase side if `d > 0` (largest bet `β = λ_max = 0.5` with weight `w = 1/6`), the decrease side if `d < 0` (`β = s_max / (B̂ - 1)`, `w = 1/5`). Let `K = 2A / w`, let F be the factor of β at its favourable extreme (`1 + β (B̂ - 1)` at `Y = B̂`, or `1 + β` at `Y = 0`), and let `ρ = E√f / √(E f) < 1` for the factor f of β under the two-point law on `{0, B̂}` with mean `1 + d`. For `0 ≤ h ≤ m` consider the stream

`Y_t = 1 + d` for `t ≤ h`, then `Y_t` i.i.d. on `{0, B̂}` with `P(Y_t = B̂) = (1 + d) / B̂` for `h < t ≤ m`

(returns `±√(m_F(τ̂))`, then `0` or `±τ̂`). It belongs to `C(W)`, and

`P(T ≤ m | W) ≥ (S_h(β |d|) - √(K S_h(β |d|)) ρ^{m - h}) / (F K)`,

where `S_h(x) = Σ_{j=1}^{h} (1 + x)^j`; the probability is 1 if `S_h(β|d|) ≥ K`. Also, the constant stream `Y_t = 1 + d` (h = m) alarms by m exactly when `Σ_k w_k S_m(λ_k d⁺) ≥ 2A` or `Σ_l v_l S_m(s_l d⁻ / (B̂ - 1)) ≥ 2A`.

*Proof.* Every value is in `[0, B̂]` and every conditional mean is `1 + d`, so the stream is in `C(W)`. Let R be the Shiryaev-Roberts statistic of the bet β; the side's statistic is at least `w R`, so `R ≥ K` forces an alarm. After the hold phase `R_h = S_h(β|d|)`, deterministic. In the bold phase `R_{h+j} = f_{h+j} (R_{h+j-1} + 1) ≥ f_{h+j} R_{h+j-1}`, so `R_{h+j} ≥ V_j = R_h Π_{i ≤ j} f_{h+i}`. Put `μ = E f = 1 + β|d| ≥ 1` and `V'_j = V_j / μ^j`, a positive martingale with `V'_0 = R_h`. Let `σ = min {j ≤ L : V_j ≥ K}` with `L = m - h`, and `σ = L` if there is none. Optional stopping at the bounded time σ gives `R_h = E V'_σ`. On `{V_σ ≥ K}`, `V'_σ ≤ V_σ ≤ F V_{σ-1} < F K` (if `σ ≥ 1`). Otherwise `V'_L ≤ V_L < K`, so `V'_L = min(V'_L, K) ≤ √(K V'_L)`. Hence `R_h ≤ F K P(V_σ ≥ K) + √K E√(V'_L)`, and `E √(V'_L) = √(R_h) (E √(f / μ))^L = √(R_h) ρ^L`. Rearranging gives the bound, since `{V_σ ≥ K}` implies an alarm by `h + σ ≤ m`. `ρ < 1` by Jensen's inequality, because f takes two values. For the constant stream the side statistics are deterministic and equal to the sums shown, and they increase in t. ∎

`limits.hold_bold_lower_bound` evaluates the larger of the two bounds over a grid of h for each window; averaging over windows gives a lower bound on the worst-case false alarm probability over the class (`results/th5_plugin_bound.csv`, `hold_bold_lower`). Part `holdbold` runs the stream itself through the detector (`results/th5_holdbold.csv`).

**Corollary 10 (the dichotomy at m ≍ √n).** [PROVED] Assume (S) with F fixed, the median of `|r|` unique and `Var_F min(r², τ²) > 0`. Let `n → ∞`, `A = A_n`, `m = m_n` with `m / A → 0`, and let `β*_m = sup P(T ≤ m)` over streams in `C(W)` chosen as functions of W.

(i) If `m / √n → 0`, then `limsup β*_m / (m / A) ≤ 1`.

(ii) If `m / √n → ∞` and `m / log A → ∞`, then `β*_m / (m / A) → ∞`.

(iii) If `m / (√n log A) → ∞`, then `β*_m → 1`.

*Proof.* (i) is Theorem 8c. (ii) Fix `z₀ > 0` and `0 < η` small, and let `E_n = {|d| ≥ z₀ / √n, |B̂ - B| ≤ η}`; `P(E_n) → P(|N(0, σ_Y²)| ≥ z₀) > 0` by Theorem 8c. On `E_n`, β, w and F are bounded away from 0 and ∞, and ρ is at most some `ρ̄ < 1` (it is a continuous function of `(d, B̂)` that is below 1 on the compact set `|d| ≤ 1/2`, `|B̂ - B| ≤ η`). Take `L_A = ⌈log(2 √K) / log(1 / ρ̄)⌉ = O(log A)` and `h = m - L_A ≥ m / 2` (true eventually, as `m / log A → ∞`). Then `√(K S_h) ρ̄^{L_A} ≤ √(S_h) / 2 ≤ S_h / 2` because `S_h ≥ 1`, so Proposition 9 gives `P(T ≤ m | W) ≥ S_h / (2 F K)`. Proposition 1 of TH4 gives `S_h(x) ≥ h (1 + x (h + 1) / 2)`, so on `E_n`

`P(T ≤ m | W) / (m / A) ≥ (w / (8 F)) (1 + β z₀ m / (4 √n))`,

whose right side tends to ∞; with the cap at 1 and `m / A → 0`, `β*_m / (m / A) ≥ P(E_n) min{A / m, (w / (8F))(1 + β z₀ m / (4√n))} → ∞`. (iii) For the constant stream on `{|d| ≥ z₀ / √n}`, the favoured side's statistic is at least `w S_m(β z₀ / √n) ≥ w (1 + β z₀ / √n)^m`, which exceeds 2A as soon as `m log(1 + β z₀ / √n) ≥ log(2A / w)`; this holds eventually because `m / (√n log A) → ∞`. So `β*_m ≥ P(|d| ≥ z₀ / √n) - o(1) → P(|N(0, σ_Y²)| ≥ z₀)`, which is arbitrarily close to 1 for small z₀. ∎

So, for the worst stream with the right conditional mean, the plug-in keeps its nominal level when the horizon is short compared with √n and loses it by an unbounded factor when the horizon is long compared with √n, whatever the alarm level. One day (m = 390) is short compared with √n only for n of order 10⁶ or more: √n = 390 at n = 152,100 (390 trading days of minute data).

### 4.5 The inflated plug-in and the confidence-sequence baseline

`detectors.calibrate(margin_z=z)` moves the two baselines apart: `m↑ = m̂ + z se` for the increase side and `m↓ = m̂ - z se` for the decrease side, `se = sd(min(r_i², τ̂²)) / √n`, with decrease bets `s_l / (τ̂² / m↓ - 1)` (the decrease side is off if `m↓ ≤ 0` or `τ̂² ≤ m↓`). Let `G_z = {m̂ - z se ≤ m_F(τ̂) ≤ m̂ + z se}`.

**Proposition 11.** Assume (S).

(a) [PROVED] For every stream in `C(W)` and all m, A: `P(T ≤ m) ≤ (m / A) P(G_z) + P(G_z^c) ≤ m / A + P(G_z^c)`.

(b) [PROVED] Under the conditions of Corollary 10, `P(G_z^c) → 2 Φ̄(z)`, so `limsup_n sup_{m, A} [P(T ≤ m) - m / A] ≤ 2 Φ̄(z)`, and `z = Φ⁻¹(1 - δ / 2)` gives `P(T ≤ m) ≤ m / A + δ + o(1)` uniformly in m and A.

(c) [PROVED] If `m / (√n log A) → ∞` and `m / A → 0`, then `sup_{C(W)} P(T ≤ m) ≥ 2 Φ̄(z) - o(1)`. So the margin costs, in the worst case, exactly its two-sided non-coverage, and `Φ⁻¹(1 - δ / 2)` is the smallest fixed z with an asymptotic excess of at most δ in this regime.

*Proof.* (a) On `G_z` both drifts relative to the shifted baselines are at most 0 (`E[X_t | F_{t-1}] = m_F(τ̂) ≤ m↑` and `≥ m↓`), so Theorem B, conditionally on `F_0` (Theorem 1 of TH4 with zero drift), gives `P(T ≤ m | W) ≤ m / A`; on `G_z^c` use 1. (b) From the proof of Theorem 8c, `√n (m̂ - m_F(τ̂)) ⇒ N(0, Var_F min(r², τ²))`, and `√n se → sd_F min(r², τ²)` in probability (the sample variance of `min(r_i², τ̂²)` is consistent because `τ̂ → τ`), so `(m̂ - m_F(τ̂)) / se ⇒ N(0, 1)` by Slutsky's lemma. (c) Fix `η > 0`. On `{m_F(τ̂) ≥ m↑ (1 + η / √n)}` the constant stream at `m_F(τ̂)` gives the increase side a drift of at least `η / √n` and alarms by m eventually, as in Corollary 10 (iii); this event has probability tending to `Φ̄(z + η / σ_Y)`. The decrease side is symmetric. Let `η → 0`. ∎

**Relation to the confidence-sequence baseline.** The cs baseline of `docs/math.md` section 6 uses the first quarter of the window for τ and puts the two baselines at the upper and lower ends of a confidence sequence for the truncated mean, at level δ split between the sides. Theorem 6.4 (Theorem D) proves `P(T ≤ m) ≤ m / A + δ` for every n, m and A and every law, with no asymptotics: it is Proposition 11a with `P(G^c) ≤ δ` guaranteed. The price is the width. Part `plugin_bound` measures it as an effective z, `(cs upper / m_F(τ_cs) - 1) √n / σ_Y` (median over windows), with δ = 0.01: 3.95 to 5.74 on the increase side and 3.78 to 3.99 on the decrease side for Gaussian returns, 3.90 to 5.56 and 3.76 to 3.97 for Student t3, falling with n from 390 to 23,040 (`results/th5_plugin_drift.csv`), against 2.576 for the asymptotic margin with the same δ. Since the cs baseline estimates the mean from three quarters of the window, part of this factor (√(4/3) = 1.15) is the sample split.

### 4.6 Independent returns, Proposition A.6 and Lemma A.7

The class `C(W)` contains streams that a market need not produce: the hold-then-bold stream knows `d(W)`. For monitored returns that are i.i.d. F, the picture is different, and only partly proved.

* *What is proved.* The upper bounds of Theorem 8 apply. Proposition A.6 of `docs/theory_running_tau.md` gives a lower bound through windows with a positive growth rate under F, `g(W) > 0`, and Lemma A.7 proves that such windows have positive probability at every fixed n, through very calm windows whose probability is of order `π_ε^n` (about 10⁻⁵⁴⁶ for Gaussian returns at n = 390). In the joint regime that lower bound disappears: a positive growth rate needs a drift of a fixed size (heuristically `d > λ_min σ_Y² / 2`, 0.019 on the increase side and 0.012 on the decrease side for Gaussian returns), while the drift is of order `n^{-1/2}`, so the share of windows with `g(W) > 0` should tend to 0 as n grows. Part `plugin_bound` computes `g(W)` by quadrature for every window: the share is 0.83 at n = 390, 0.41 at 5760, 0.11 at 23,040 and 0.005 at 92,160 for Gaussian returns, and 0.76, 0.25, 0.021 and 0 for Student t3 (`results/th5_plugin_drift.csv`, 20,000 windows per cell; at 390 and 1440 steps these are within 1.1 points of E14 Q14w, `results/e14w_plugin_drift.csv`, which used 2000 windows and a different quadrature). The calm windows of Lemma A.7 are the large-deviation version of the same mechanism; Corollary 10 replaces them by the typical windows of the central limit theorem, at the price of an adversarial stream.
* *What is heuristic.* [HEURISTIC] For one bet with a small drift `d` of order `n^{-1/2}` and negative growth, `E f = 1 + λ d > 1` but `E log f < 0`. The recursion `R = f (R + 1)` then has a stationary law with a power tail `P(R > x) ≈ C x^{-κ}`, where `E f^κ = 1` (Kesten 1973), and `κ ≈ 1 - 2 d / (λ σ_Y²)` for small d. The alarm rate per step is then about `(2A)^{-κ}`, an inflation of `A^{1-κ} = exp(2 d log A / (λ σ_Y²))` over the nominal rate. With `d ≈ σ_Y Z / √n` this is `exp(2 Z log A / (λ σ_Y √n))`: for i.i.d. returns the relevant scale appears to be `log A / √n` rather than `m / √n`, once m is long enough for the statistic to reach its stationary law. We have not turned this into a proof.

### 4.7 Numbers

[NUMERICAL, SYNTHETIC] A = 7800, c = 3, each side at 2A, ratios to the nominal `m / A`. Part `plugin_bound` draws calibration windows only (20,000 per cell for n ≤ 23,040 and 5000 above) and evaluates the bounds window by window (`results/th5_plugin_drift.csv`, `results/th5_plugin_bound.csv`, `results/th5_distribution_free.csv`). Part `holdbold` runs the stream of Proposition 9 through the detector (`results/th5_holdbold.csv`). Part `plugin_sim` monitors 7800 i.i.d. returns of the calibration law after each of 10,000 windows per cell (`results/th5_plugin_sim.csv`). `results/figures/th5_plugin_regime.png` plots the bounds and the simulated rates against `m / √n`.

**The drift (Lemma 7a, Theorem 8c, section 4.6).**

| Law | n | windows | `√n sd(d)` (limit `σ_Y`) | `√n E d⁺` (limit `σ_Y / √(2π)`) | max `abs(d) / (B̂ D_n)` | share with positive growth under F (up, down, either) | cs effective z (up, down) |
|---|---|---|---|---|---|---|---|
| gaussian | 390 | 20,000 | 1.398 (1.386) | 0.596 (0.553) | 0.544 | 0.401, 0.425, 0.826 | 5.74, 3.99 |
| gaussian | 1,440 | 20,000 | 1.392 (1.386) | 0.579 (0.553) | 0.494 | 0.308, 0.365, 0.673 | 5.18, 3.86 |
| gaussian | 2,880 | 20,000 | 1.382 (1.386) | 0.574 (0.553) | 0.498 | 0.236, 0.313, 0.549 | 4.54, 3.83 |
| gaussian | 5,760 | 20,000 | 1.382 (1.386) | 0.567 (0.553) | 0.509 | 0.155, 0.251, 0.406 | 4.26, 3.80 |
| gaussian | 11,520 | 20,000 | 1.391 (1.386) | 0.569 (0.553) | 0.503 | 0.075, 0.174, 0.249 | 4.08, 3.78 |
| gaussian | 23,040 | 20,000 | 1.386 (1.386) | 0.563 (0.553) | 0.499 | 0.020, 0.092, 0.113 | 3.95, 3.80 |
| gaussian | 46,080 | 5,000 | 1.378 (1.386) | 0.550 (0.553) | 0.497 | 0.002, 0.029, 0.031 | not run |
| gaussian | 92,160 | 5,000 | 1.386 (1.386) | 0.558 (0.553) | 0.473 | 0.000, 0.005, 0.005 | not run |
| t3 | 390 | 20,000 | 1.589 (1.583) | 0.687 (0.631) | 0.732 | 0.390, 0.369, 0.759 | 5.56, 3.97 |
| t3 | 1,440 | 20,000 | 1.578 (1.583) | 0.668 (0.631) | 0.747 | 0.288, 0.270, 0.558 | 4.40, 3.86 |
| t3 | 2,880 | 20,000 | 1.586 (1.583) | 0.661 (0.631) | 0.741 | 0.208, 0.202, 0.410 | 4.17, 3.82 |
| t3 | 5,760 | 20,000 | 1.588 (1.583) | 0.654 (0.631) | 0.715 | 0.126, 0.122, 0.248 | 4.01, 3.80 |
| t3 | 11,520 | 20,000 | 1.562 (1.583) | 0.644 (0.631) | 0.745 | 0.048, 0.047, 0.095 | 3.93, 3.79 |
| t3 | 23,040 | 20,000 | 1.587 (1.583) | 0.630 (0.631) | 0.752 | 0.011, 0.011, 0.021 | 3.90, 3.76 |
| t3 | 46,080 | 5,000 | 1.601 (1.583) | 0.647 (0.631) | 0.702 | 0.000, 0.001, 0.001 | not run |
| t3 | 92,160 | 5,000 | 1.592 (1.583) | 0.637 (0.631) | 0.755 | 0.000, 0.000, 0.000 | not run |

* `√n sd(d)` is within 1.3% of `σ_Y` at every n for both laws, already at n = 390, so the normal limit of Theorem 8c is accurate at these sizes. `√n E d⁺` is 8% (Gaussian) and 9% (t3) above its limit at n = 390 and within 1% of it at n = 92,160.
* Lemma 7a holds with room: the largest `|d| / (B̂ D_n)` in any cell is 0.55 (Gaussian) and 0.76 (t3).
* The share of windows whose plug-in has a positive growth rate under F falls from 0.83 to 0.005 (Gaussian) and from 0.76 to 0 (t3) as n goes from 390 to 92,160.

**Upper and lower bounds (Theorems 8a and 8c, Proposition 9).** The table shows all Gaussian cells and four of the eight values of n for t3; the CSV has all of them. Every ratio is capped at `A / m` (79.6 at m = 98, 20 at m = 390, 4 at m = 1950).

| Law | n | m | `m / √n` | upper: Corollary 1 | upper: first order `1 + C_F (m+1)/√n` | upper: Theorem 8a | lower: Proposition 9 | noise-free stream, P(alarm) | inflated z = 2.576, upper |
|---|---|---|---|---|---|---|---|---|---|
| gaussian | 390 | 98 | 4.96 | 3.034 | 1.185 | 11.41 | 0.635 | 0.007 | 1.004 |
| gaussian | 390 | 390 | 19.75 | 8.316 | 1.730 | 12.25 | 6.460 | 0.314 | 1.056 |
| gaussian | 390 | 1,950 | 98.74 | 3.704 | 4.642 | 3.84 | 3.477 | 0.856 | 1.030 |
| gaussian | 1,440 | 98 | 2.58 | 1.212 | 1.096 | 2.72 | 0.009 | 0.000 | 1.001 |
| gaussian | 1,440 | 390 | 10.28 | 5.156 | 1.380 | 8.36 | 3.314 | 0.159 | 1.021 |
| gaussian | 1,440 | 1,950 | 51.39 | 3.446 | 2.895 | 3.70 | 3.031 | 0.734 | 1.021 |
| gaussian | 2,880 | 98 | 1.83 | 1.106 | 1.068 | 1.46 | 0.006 | 0.000 | 1.000 |
| gaussian | 2,880 | 390 | 7.27 | 3.587 | 1.269 | 6.56 | 1.732 | 0.080 | 1.005 |
| gaussian | 2,880 | 1,950 | 36.34 | 3.225 | 2.340 | 3.58 | 2.664 | 0.635 | 1.017 |
| gaussian | 5,760 | 98 | 1.29 | 1.064 | 1.048 | 1.16 | 0.005 | 0.000 | 1.000 |
| gaussian | 5,760 | 390 | 5.14 | 2.242 | 1.190 | 4.94 | 0.567 | 0.023 | 1.002 |
| gaussian | 5,760 | 1,950 | 25.69 | 2.953 | 1.948 | 3.41 | 2.227 | 0.519 | 1.010 |
| gaussian | 11,520 | 98 | 0.91 | 1.041 | 1.034 | 1.08 | 0.005 | 0.000 | 1.000 |
| gaussian | 11,520 | 390 | 3.63 | 1.480 | 1.134 | 3.37 | 0.117 | 0.002 | 1.001 |
| gaussian | 11,520 | 1,950 | 18.17 | 2.646 | 1.670 | 3.20 | 1.774 | 0.403 | 1.009 |
| gaussian | 23,040 | 98 | 0.65 | 1.027 | 1.024 | 1.04 | 0.005 | 0.000 | 1.000 |
| gaussian | 23,040 | 390 | 2.57 | 1.190 | 1.095 | 2.08 | 0.048 | 0.000 | 1.000 |
| gaussian | 23,040 | 1,950 | 12.85 | 2.322 | 1.474 | 2.92 | 1.341 | 0.299 | 1.005 |
| gaussian | 46,080 | 98 | 0.46 | 1.018 | 1.017 | 1.02 | 0.004 | 0.000 | 1.000 |
| gaussian | 46,080 | 390 | 1.82 | 1.102 | 1.067 | 1.40 | 0.040 | 0.000 | 1.000 |
| gaussian | 46,080 | 1,950 | 9.08 | 1.997 | 1.335 | 2.56 | 0.954 | 0.209 | 1.004 |
| gaussian | 92,160 | 98 | 0.32 | 1.013 | 1.012 | 1.02 | 0.004 | 0.000 | 1.000 |
| gaussian | 92,160 | 390 | 1.28 | 1.062 | 1.047 | 1.15 | 0.037 | 0.000 | 1.000 |
| gaussian | 92,160 | 1,950 | 6.42 | 1.745 | 1.237 | 2.25 | 0.613 | 0.127 | 1.002 |
| t3 | 390 | 98 | 4.96 | 4.127 | 1.236 | 14.44 | 1.279 | 0.015 | 1.004 |
| t3 | 390 | 390 | 19.75 | 10.787 | 1.933 | 14.67 | 8.880 | 0.429 | 1.061 |
| t3 | 390 | 1,950 | 98.74 | 3.804 | 5.654 | 3.89 | 3.656 | 0.906 | 1.031 |
| t3 | 2,880 | 98 | 1.83 | 1.146 | 1.087 | 1.78 | 0.007 | 0.000 | 1.000 |
| t3 | 2,880 | 390 | 7.27 | 4.389 | 1.343 | 8.36 | 2.411 | 0.110 | 1.006 |
| t3 | 2,880 | 1,950 | 36.34 | 3.483 | 2.713 | 3.72 | 3.096 | 0.753 | 1.018 |
| t3 | 23,040 | 98 | 0.65 | 1.035 | 1.031 | 1.06 | 0.005 | 0.000 | 1.000 |
| t3 | 23,040 | 390 | 2.57 | 1.281 | 1.121 | 2.53 | 0.074 | 0.000 | 1.001 |
| t3 | 23,040 | 1,950 | 12.85 | 2.681 | 1.605 | 3.24 | 1.797 | 0.406 | 1.008 |
| t3 | 92,160 | 98 | 0.32 | 1.016 | 1.015 | 1.02 | 0.005 | 0.000 | 1.000 |
| t3 | 92,160 | 390 | 1.28 | 1.082 | 1.061 | 1.22 | 0.047 | 0.000 | 1.000 |
| t3 | 92,160 | 1,950 | 6.42 | 1.952 | 1.303 | 2.64 | 0.810 | 0.164 | 1.002 |

* *Short horizons.* Where `m / √n < 1` (m = 98 and n ≥ 11,520), Corollary 1 is at most 4.1% (Gaussian) and 5.3% (t3) above nominal, the first-order formula of Theorem 8c is within 1% of it (1.041 against 1.034 at n = 11,520, 1.013 against 1.012 at n = 92,160, Gaussian), and the explicit bound of Theorem 8a is within 5% of it.
* *Long horizons.* In 19 of the 24 cells with `m / √n > 5` the proved lower bound of Proposition 9 is above nominal, up to 6.5 times (Gaussian) and 8.9 times (t3) at n = m = 390, and it is 25% to 96% of Corollary 1's upper bound (median 67%). The five exceptions (lower bound 0.57 to 0.95 times nominal) are at `m / √n` between 5.1 and 9.1; Corollary 10 (ii) is a limit statement, and at a fixed A the hold phase needs a large enough `m |d|` to lift the statistic. The first-order formula is far below Corollary 1 in this range (1.19 against 2.24 at n = 5760, m = 390): outside `m / √n → 0` the constant `C_F` does not describe the excess, as Corollary 10 requires.
* The noise-free stream alone alarms by m with probability 0.31 at n = m = 390 and 0.86 at n = 390, m = 1950 (Gaussian).
* *Inflated plug-in.* With z = 2.576, Corollary 1 for the shifted baselines is at most 1.061 times nominal in every cell, and the share of windows with `m_F(τ̂)` outside `[m̂ - z se, m̂ + z se]` is 0.009 to 0.013 against the limit `2 Φ̄(2.576) = 0.01` of Proposition 11b (0.305 to 0.323 against 0.317 for z = 1, and 0.041 to 0.050 against 0.046 for z = 2).

**The distribution-free bound (Theorem 8b), ratio to nominal, capped at `A / m`.**

| n | m = 5 | m = 10 | m = 39 | m = 98 | m = 390 | m = 1950 |
|---|---|---|---|---|---|---|
| 390 | 255.120 | above 1000 | above 1000 | above 1000 | above 1000 | above 1000 |
| 1,440 | 5.597 | 339.105 | above 1000 | above 1000 | above 1000 | above 1000 |
| 5,760 | 1.524 | 5.214 | above 1000 | above 1000 | above 1000 | above 1000 |
| 23,040 | 1.134 | 1.480 | 249.441 | above 1000 | above 1000 | above 1000 |
| 92,160 | 1.048 | 1.122 | 4.540 | above 1000 | above 1000 | above 1000 |
| 1,000,000 | 1.008 | 1.015 | 1.153 | 2.316 | above 1000 | above 1000 |
| 10,000,000 | 1.002 | 1.004 | 1.017 | 1.053 | 1.670 | above 1000 |
| 100,000,000 | 1.001 | 1.001 | 1.005 | 1.013 | 1.073 | 10.210 |

It is informative only when m is far below √n: a factor below 1.2 needs n ≥ 10⁶ at m = 39 and n = 10⁸ at m = 390 (1.073).

**Proposition 9 checked by simulation** (Gaussian, 4000 windows per cell; in each window the stream with the hold length that maximises the proved bound; probabilities, not ratios).

| n | m | nominal | proved lower bound (mean over windows) | simulated rate (95% interval) | noise-free stream | Corollary 1 upper bound |
|---|---|---|---|---|---|---|
| 390 | 390 | 0.05 | 0.338 | 0.365 (0.350 to 0.380) | 0.330 | 0.424 |
| 390 | 1,950 | 0.25 | 0.864 | 0.893 (0.883 to 0.902) | 0.850 | 0.924 |
| 2,880 | 390 | 0.05 | 0.087 | 0.111 (0.101 to 0.121) | 0.080 | 0.180 |
| 2,880 | 1,950 | 0.25 | 0.663 | 0.723 (0.709 to 0.736) | 0.632 | 0.805 |

In all four cells the simulated rate lies between the proved lower bound and Corollary 1's upper bound, and the lower bound is 78% to 97% of the simulated rate.

**Independent monitored returns** (rates in percent; 10,000 paths per cell; the Wilson 95% intervals in the CSV have half-widths of at most 1.0 point).

| Law | n | m (nominal) | plug-in | inflated z = 1 | inflated z = 2 | inflated z = 2.576 | confidence sequence | true baseline |
|---|---|---|---|---|---|---|---|---|
| gaussian | 390 | 390 (5%) | 3.61 | 0.65 | 0.11 | 0.03 | 0.00 | 0.94 |
| gaussian | 390 | 1950 (25%) | 40.77 | 9.58 | 1.58 | 0.42 | 0.02 | 12.15 |
| gaussian | 390 | 7800 (100%) | 85.23 | 29.28 | 5.76 | 2.11 | 0.03 | 55.41 |
| gaussian | 1440 | 390 (5%) | 1.62 | 0.57 | 0.25 | 0.13 | 0.02 | 0.92 |
| gaussian | 1440 | 1950 (25%) | 23.40 | 7.39 | 2.07 | 0.96 | 0.07 | 12.36 |
| gaussian | 1440 | 7800 (100%) | 76.20 | 31.58 | 9.25 | 4.41 | 0.35 | 55.89 |
| gaussian | 2880 | 390 (5%) | 1.43 | 0.75 | 0.32 | 0.23 | 0.09 | 1.16 |
| gaussian | 2880 | 1950 (25%) | 18.54 | 7.84 | 3.13 | 1.91 | 0.42 | 12.63 |
| gaussian | 2880 | 7800 (100%) | 70.20 | 34.28 | 13.64 | 7.88 | 1.71 | 56.04 |
| gaussian | 5760 | 390 (5%) | 1.19 | 0.67 | 0.37 | 0.31 | 0.14 | 1.06 |
| gaussian | 5760 | 1950 (25%) | 14.70 | 7.70 | 3.94 | 2.66 | 1.07 | 12.34 |
| gaussian | 5760 | 7800 (100%) | 64.27 | 37.60 | 18.95 | 12.96 | 5.24 | 56.10 |
| gaussian | 11520 | 390 (5%) | 1.11 | 0.69 | 0.42 | 0.33 | 0.24 | 1.01 |
| gaussian | 11520 | 1950 (25%) | 14.06 | 8.86 | 5.38 | 3.92 | 2.03 | 12.47 |
| gaussian | 11520 | 7800 (100%) | 61.54 | 40.46 | 25.30 | 18.79 | 9.89 | 56.23 |
| gaussian | 23040 | 390 (5%) | 1.02 | 0.75 | 0.55 | 0.46 | 0.32 | 0.85 |
| gaussian | 23040 | 1950 (25%) | 13.33 | 9.92 | 7.10 | 5.93 | 3.96 | 12.51 |
| gaussian | 23040 | 7800 (100%) | 58.90 | 44.29 | 31.71 | 26.01 | 16.90 | 55.87 |
| t3 | 390 | 390 (5%) | 5.14 | 1.14 | 0.17 | 0.06 | 0.01 | 1.64 |
| t3 | 390 | 1950 (25%) | 42.34 | 10.45 | 1.64 | 0.55 | 0.01 | 14.50 |
| t3 | 390 | 7800 (100%) | 86.06 | 32.23 | 6.44 | 2.24 | 0.06 | 57.12 |
| t3 | 1440 | 390 (5%) | 2.26 | 0.84 | 0.28 | 0.15 | 0.00 | 1.60 |
| t3 | 1440 | 1950 (25%) | 23.81 | 8.42 | 2.72 | 1.53 | 0.25 | 13.96 |
| t3 | 1440 | 7800 (100%) | 75.82 | 34.67 | 11.70 | 6.32 | 1.23 | 56.67 |
| t3 | 2880 | 390 (5%) | 1.82 | 0.92 | 0.49 | 0.37 | 0.14 | 1.54 |
| t3 | 2880 | 1950 (25%) | 19.58 | 8.92 | 3.89 | 2.63 | 0.97 | 14.29 |
| t3 | 2880 | 7800 (100%) | 70.53 | 37.54 | 16.96 | 10.66 | 3.85 | 56.63 |
| t3 | 5760 | 390 (5%) | 1.73 | 1.10 | 0.63 | 0.49 | 0.30 | 1.57 |
| t3 | 5760 | 1950 (25%) | 16.96 | 9.96 | 5.77 | 4.00 | 2.07 | 14.96 |
| t3 | 5760 | 7800 (100%) | 64.86 | 40.27 | 22.93 | 16.36 | 8.35 | 57.18 |
| t3 | 11520 | 390 (5%) | 1.72 | 1.17 | 0.84 | 0.64 | 0.48 | 1.62 |
| t3 | 11520 | 1950 (25%) | 15.34 | 9.95 | 6.67 | 5.18 | 3.03 | 13.85 |
| t3 | 11520 | 7800 (100%) | 61.75 | 42.37 | 28.69 | 22.76 | 13.51 | 56.34 |
| t3 | 23040 | 390 (5%) | 1.61 | 1.20 | 0.89 | 0.75 | 0.63 | 1.48 |
| t3 | 23040 | 1950 (25%) | 14.36 | 10.87 | 8.04 | 6.62 | 4.78 | 13.56 |
| t3 | 23040 | 7800 (100%) | 58.83 | 44.99 | 33.98 | 28.57 | 20.56 | 55.80 |

* *The true baseline is conservative for i.i.d. returns.* At one day it alarms in 0.85% to 1.64% of paths against the nominal 5%, and at five days in 12.2% to 15.0% against 25%: for i.i.d. returns the bound `m / A` of Theorem B is loose by a factor of 3 to 6 at m = 390 and 1.7 to 2.1 at m = 1950.
* *The plug-in.* At one day it is at or below nominal in every cell (largest: 5.14% for t3 at n = 390, interval 4.72% to 5.59%). At five days with one day of calibration it is not: 40.8% (Gaussian) and 42.3% (t3) against 25% at n = 390, and 23.4% and 23.8% at n = 1440.
* *The excess over the true baseline falls with n.* At m = 1950 the ratio of the plug-in's rate to the true baseline's is 3.4, 1.9, 1.5, 1.19, 1.13 and 1.07 (Gaussian) and 2.9, 1.7, 1.4, 1.13, 1.11 and 1.06 (t3) for n = 390 to 23,040, while `m / √n` falls only from 99 to 13 and the worst stream in `C(W)` is still above nominal at n = 23,040 (Proposition 9: 1.34 times). For i.i.d. returns the joint regime is therefore much milder than for the worst stream. This is consistent with the heuristic scale `log A / √n` of section 4.6 (`log(2A) / √n` is 0.49 at n = 390 and 0.064 at n = 23,040), but the comparison was made after seeing the results and is not a test of it [post hoc].
* No simulated rate is above its Corollary 1 bound: in all 288 rows of the CSV the lower end of the Wilson interval is below the bound averaged over the same windows.
* The inflated plug-ins and the confidence-sequence baseline are conservative for i.i.d. returns: at one day z = 2.576 gives 0.03% to 0.75% and the cs baseline 0% to 0.63%. Their price is a later alarm after a change, which this note does not measure.

## 5. The decrease bet at its boundary

For a decrease bet write `λ = λ_max (1 - ε)`, `0 < ε ≤ 1`, so that its factor is `f_ε(y) = ε + (1 - ε) f₀(y)` with `f₀(y) = (B - y) / (B - 1)`, the factor at `λ_max`. Let `q_B = Q(Y = B)` (the probability that a post-change return reaches the truncation level), `G_c = E_Q[log f₀(Y); Y < B]` (the growth at `λ_max` without the atom) and `g(ε) = E_Q log f_ε(Y)`.

**Proposition 12.** [PROVED] Let Q be i.i.d. with `E_Q[log f₀(Y); Y < B] > -∞`.

(a) `(1 - ε) G_c + q_B log ε ≤ g(ε) ≤ G_c + q_B log ε + E_Q[log(1 + ε / f₀(Y)); Y < B]`.

(b) g is concave in ε, and the growth-optimal bet is strictly inside, `λ* < λ_max`, exactly when `E_Q[(B - 1) / (B - Y)] > 1`. This holds whenever `q_B > 0`, and also whenever Y has a positive density near B.

(c) At an interior optimum `ε* = q_B / h(ε*)` with `h(ε) = E_Q[(f₀(Y) - 1) / (ε + (1 - ε) f₀(Y)); Y < B]`. So `ε*` is proportional to the atom, with a factor `1 / h(ε*)` that depends on the bulk of the law.

(d) With `ε₀ = q_B / G_c` (if `G_c > q_B`), `g(ε₀) ≥ G_c - q_B (1 + log(G_c / q_B))`.

(e) For the Shiryaev-Roberts statistic of this single bet at level A, for every change time and pre-change history: `WADD ≤ (log A + log(B / (B - 1))) / g(ε)` if `g(ε) > 0` (Theorem 1 of `docs/theory_delay.md`), and `WADD ≤ min_n n / p_n(ε)` for every ε, where `p_n(ε) = Q(Σ_{i ≤ n} log f_ε(Y_i) ≥ log A)`. Since `f_ε ≥ (1 - ε) f₀` on `{Y < B}`, `p_n(ε) ≥ (1 - q_B)^n Q(Σ_{i ≤ n} log f₀(Y'_i) ≥ log A + n log(1 / (1 - ε)))` with `Y'` drawn from Q given `Y < B`, which does not tend to 0 as `ε → 0`. So the delay stays bounded as `ε → 0`, even where `g(ε) → -∞`.

*Proof.* (a) On `{Y = B}`, `f_ε = ε`. On `{Y < B}`, `f_ε` is the average of `f₀` and 1 with weights `1 - ε` and ε, so `log f_ε ≥ (1 - ε) log f₀` by concavity of the logarithm; and `f_ε ≤ f₀ + ε`, so `log f_ε ≤ log f₀ + log(1 + ε / f₀)`. (b) `log f_ε` is concave in ε (the log of an affine function), and `g'(0+) = E[(1 - f₀) / f₀] = E[(B - 1) / (B - Y)] - 1` (by monotone convergence; `+∞` if `q_B > 0`, since `1 / f₀ = ∞` at B, and also if the density of Y is positive near B, since `∫ dy / (B - y)` diverges). The optimum is at ε = 0 (`λ* = λ_max`) exactly when `g'(0+) ≤ 0`. (c) `g'(ε) = q_B / ε - h(ε) = 0`. (d) Put `ε₀ G_c = q_B` into the lower bound of (a). (e) The first bound is Theorem 1 of `docs/theory_delay.md` with one bet (weight 1, largest log factor `log(1 + λ) ≤ log(B / (B - 1))`). For the second, the statistic at the end of any block of n consecutive steps is at least the product of the factors in the block, so the detector has alarmed by the end of the first block whose log factors sum to at least log A. Disjoint blocks are independent with success probability `p_n`, so `P(T > k n) ≤ (1 - p_n)^k` and `E T ≤ n Σ_{k ≥ 0} (1 - p_n)^k = n / p_n`. The same holds after any pre-change history, since the statistic is nonnegative. ∎

**Numbers.** [NUMERICAL, SYNTHETIC] Part `boundary` (`results/th5_boundary.csv`, 4000 paths of 6000 steps per law, A = 7800):

| Law | sd ratio | `q_B` | `G_c` | `ε*` | `λ* / λ_max` | `KL_inf` | growth at ε = 10⁻⁶ | delay at `ε*` | at ε = 0.05 | at ε = 10⁻⁶ | Wald bound at `ε*` | block bound at `ε*` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gaussian | 0.5 | 1.97e-09 | 0.0881 | 2.37e-08 | 1.0000 | 0.0881 | 0.0881 | 74.1 | 77.2 | 74.1 | 103.1 | 115 |
| gaussian | 0.67 | 7.55e-06 | 0.0627 | 1.34e-04 | 0.9999 | 0.0626 | 0.0626 | 99.0 | 103.0 | 99.0 | 145.1 | 186 |
| gaussian | 0.8 | 1.77e-04 | 0.0357 | 8.35e-03 | 0.9917 | 0.0346 | 0.0332 | 156.4 | 160.1 | 155.8 | 262.0 | 396 |
| student t3 | 0.5 | 6.49e-03 | 0.1146 | 6.90e-02 | 0.9310 | 0.0913 | 0.0249 | 69.7 | 69.3 | 71.3 | 100.0 | 155 |
| student t3 | 0.67 | 1.48e-02 | 0.0778 | 2.23e-01 | 0.7767 | 0.0450 | -0.1260 | 122.8 | 126.3 | 139.9 | 203.2 | 338 |
| student t3 | 0.8 | 2.38e-02 | 0.0475 | 4.55e-01 | 0.5451 | 0.0172 | -0.2817 | 258.8 | 327.8 | 396.7 | 530.8 | 909 |
| jumps | 0.5 | 3.53e-03 | 0.0899 | 4.18e-02 | 0.9582 | 0.0752 | 0.0411 | 82.4 | 82.8 | 82.8 | 120.8 | 172 |
| jumps | 0.67 | 3.54e-03 | 0.0649 | 6.01e-02 | 0.9399 | 0.0515 | 0.0160 | 112.8 | 112.4 | 113.5 | 176.4 | 267 |
| jumps | 0.8 | 3.70e-03 | 0.0384 | 1.23e-01 | 0.8773 | 0.0273 | -0.0127 | 183.5 | 179.7 | 184.7 | 332.3 | 538 |

Delays are simulated means (change at the first step, single bet at threshold A); the bounds are those of (e).

* **Light tails: the boundary is harmless.** For Gaussian returns `q_B` is `2.0 × 10⁻⁹` (sd ratio 0.5) to `1.8 × 10⁻⁴` (0.8), and `ε*` is `2.4 × 10⁻⁸` to `8.3 × 10⁻³`, the values `λ* / λ_max = 1 - 2.4 × 10⁻⁸` to 0.9917 that the referee reported. Every ε of 10⁻³ or less gives a growth within 4.2% of `KL_inf` and a simulated delay within 0.4% of the delay at `ε*`; the choice `ε₀ = q_B / G_c` of (d) gives a delay within 0.2% of it (within 3% in all nine laws).
* **Fat tails: the atom decides.** For Student t3 `q_B` is 0.65% to 2.4% and `λ* / λ_max` is 0.93, 0.78 and 0.55 at sd ratios 0.5, 0.67 and 0.8. A bet at `λ_max (1 - 10⁻⁶)` pays `q_B log(10⁶)` nats per step for the atom: its growth is negative at sd ratios 0.67 and 0.8 (-0.126 and -0.282), and its delay is 14% and 53% longer than at `ε*` (140 against 123 steps, 397 against 259). For the jump law, whose atom is small (0.35% to 0.37%) but whose bulk is Gaussian, the boundary costs at most 1% in delay even where the growth is negative.
* **The growth rate is the wrong summary near the boundary.** Where `g(ε) < 0` the Wald bound is infinite but the delay is not: the Shiryaev-Roberts statistic restarts after each atom and alarms on a long enough run without one. The block bound of (e) captures this but is loose (for t3 at 0.8 it is about 12,000 steps against a simulated 397).

So the referee's observation that the decrease optimum sits at the boundary holds for light tails, where it follows from `q_B ≈ 0`; it does not hold for fat tails, where the default grid's largest decrease bet, `0.95 λ_max`, is already beyond the optimum.

## 6. Status

* **Proved here:** Lemma 1; Theorem 2 (deterministic rules, randomised rules, sharpness); Corollaries 3 and 4; Theorem 5 (given the Neyman-Pearson argument, proved in the text); Proposition 6 (given the strong law, the elementary renewal theorem and the weak law); Lemma 7; Theorem 8 (given TH4's Corollary 1, the DKW inequality with Massart's constant and the central limit theorem); Proposition 9; Corollary 10; Proposition 11; Proposition 12.
* **Standard arguments, given for completeness:** the change-of-measure step of Proposition 6c (Lorden 1971, Lai 1998).
* **Numerical:** the Neyman-Pearson bounds of section 3 are Monte Carlo estimates (importance weights from 10,000 paths per law); the floor means are Monte Carlo averages of an exact pathwise quantity; the growth-rate share of section 4.6 uses quadrature per window.
* **Heuristic:** the `log A / √n` scale for i.i.d. monitored returns (section 4.6); the expression of `ε*` in terms of the bulk is exact but `h(ε*)` is computed numerically.
* **Not done:** the delay floor for a change after a pre-change history (Theorem 5 covers the change at the first step, which is the worst case in Lorden's sense); a proof for i.i.d. monitored returns in the joint regime; a rule that attains the delay floor uniformly over change times; the delay cost of the margin z of Proposition 11 and of the cs baseline.
* **Failures recorded in `docs/iterations/TH5.md`:** the least favourable law from `design.best_bet` had mass 0.982 for the Gaussian halving (fixed by solving the first-order condition in `log ε`; the estimates did not change); the explicit and distribution-free bounds overflowed before they were computed in log space and capped per window; a wrong factor in a first evaluation of the default grid's collapse floor.

## 7. References

* Bahadur, R. R. and Savage, L. J. (1956). The nonexistence of certain statistical procedures in nonparametric problems. *Annals of Mathematical Statistics* 27(4), 1115 to 1122. (Cited as in `docs/math.md`.)
* Clerico, E. (2025). On the optimality of coin-betting for mean estimation. *International Journal of Approximate Reasoning* 187, 109550. arXiv:2412.02640. (Verified from search listings; the paper itself could not be opened from this machine.)
* Dvoretzky, A., Kiefer, J. and Wolfowitz, J. (1956). Asymptotic minimax character of the sample distribution function and of the classical multinomial estimator. *Annals of Mathematical Statistics* 27(3), 642 to 669. doi:10.1214/aoms/1177728174. (Verified from search listings.)
* Kesten, H. (1973). Random difference equations and renewal theory for products of random matrices. *Acta Mathematica* 131, 207 to 248. doi:10.1007/BF02392040. (Verified from search listings; used only for the heuristic of section 4.6.)
* Lai, T. L. (1998). Information bounds and quick detection of parameter changes in stochastic systems. *IEEE Transactions on Information Theory* 44(7), 2917 to 2929. (Verified from search listings; the description of its window constraint is from secondary sources.)
* Lorden, G. (1971). Procedures for reacting to a change in distribution. *Annals of Mathematical Statistics* 42(6), 1897 to 1908. (Cited as in `docs/theory_delay.md`.)
* Massart, P. (1990). The tight constant in the Dvoretzky-Kiefer-Wolfowitz inequality. *Annals of Probability* 18(3), 1269 to 1283. doi:10.1214/aop/1176990746. (Verified from search listings.)
* Ram, A. and Ramdas, A. (2026). Asymptotically optimal sequential change detection for bounded means. arXiv:2602.05272. (Cited as in `docs/theory_delay.md`.)
* Shafer, G. and Vovk, V. (2019). *Game-Theoretic Foundations for Probability and Finance.* Wiley. (Verified from search listings.)
* Waudby-Smith, I., Sandoval, R. and Jordan, M. I. (2025). Universal log-optimality for general classes of e-processes and sequential hypothesis tests. arXiv:2504.02818. (Verified from search listings; the paper itself could not be opened from this machine.)
