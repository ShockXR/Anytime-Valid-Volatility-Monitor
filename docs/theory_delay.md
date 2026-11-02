# Detection delay theory for the betting e-detector

`docs/math.md` (section 5) proves that the detector rarely raises false alarms. This document is the other half of a sequential detection theory: how fast the detector reacts to a real change, how the bets should be chosen, what robustness costs compared with the best parametric method, why the false alarm bound looks loose, and whether the detector is optimal as false alarms become rare. Every number comes from `experiments/e11_delay_theory.py` (outputs `results/e11*.csv`, all **[SYNTHETIC]**), except the test of a separate decrease clip in section 6.3, which comes from `experiments/e11f_split_c.py` (outputs `results/e11f_split_c_synthetic.csv` **[SYNTHETIC]** and `results/e11f_split_c_real.csv` **[REAL]**). The helper functions are in `src/avmonitor/design.py`. The experiment write-up is `docs/report_parts/E11.md`.

Each claim carries one of four labels:

* **[PROVED]**: a complete proof is given here.
* **[STANDARD]**: a known result from the cited literature, restated for this setting (we checked that the references exist and what they cover, but we do not repeat their proofs).
* **[HEURISTIC]**: an approximation derived by a standard but non-rigorous argument. Its accuracy is checked by simulation.
* **[NUMERICAL]**: established only by numerical integration or simulation.

## 0. Summary in plain words

* **Speed is governed by one number, the growth rate.** After a change, the logarithm of the detector's "wealth" grows by a fixed amount per step on average, `g`. The detector alarms when the log wealth reaches `log A` (A is the alarm level, 7800 in the experiments). So the delay is roughly `log(A) / g` steps. This first-order formula is right in the limit of very rare false alarms, but at practical alarm levels it overstates the delay of small changes by up to 112%; it looks accurate for large changes only because two errors cancel. The reason is that the Shiryaev-Roberts statistic adds up evidence from many possible change times at once. A second-order formula, `(log(A) + κ - ζ) / g`, with a "memory" term `ζ ≈ log(1/g) + 0.58` and an overshoot term `κ`, predicts the simulated delay of a single bet to within 3% in all 20 cases tried.
* **A proved, non-asymptotic delay bound.** For every alarm level, every change time and every pre-change history, the expected delay is at most `min_k (log(A / w_k) + overshoot_k) / g_k`, where `w_k` is the weight of bet k. In the simulations it is 1.3 to 2.7 times the true delay.
* **The price of robustness is a ratio of information numbers.** The Gaussian likelihood ratio detector (the best method when returns are exactly normal) grows at the Kullback-Leibler divergence KL; the robust detector grows at `g ≤ KL`. The ratio `KL / g` splits into three proved factors: information lost by truncating, information lost by betting only on the mean (the price of making no distributional assumption), and the cost of a finite grid of bets. On clean Gaussian data with c = 3 the total is 1.16, 1.37 and 1.95 for standard deviation ratios 1.25, 1.5 and 2. This explains the 10% to 70% delay penalty of E4: the simulated penalties against the Gaussian likelihood ratio mixture of E4 are 1.14, 1.35 and 1.73. The largest factor is the nonparametric one, not truncation.
* **Truncation sets a speed limit for decreases.** No bet that is valid for every return distribution with the same truncated variance can earn more than `log(B / (B - 1)) ≈ 1/B` nats of evidence in one step against a decrease, where `B = τ² / m₀ ≈ c²`. This cap applies to every valid bet, not only to our grid. With c = 3 and A = 7800, a Shiryaev-Roberts e-detector built from such bets needs at least 58 steps to alarm even when volatility collapses to zero (64 at the two-sided level 2A), and the default five-bet decrease grid needs 73 (79), because it spreads its weight over smaller bets (Corollary 4). This is why E4's decrease delays are slow, and why a smaller c (about 1.5 to 2) on the decrease side should cut them by a quarter to a half. E11f confirms this for continuous return laws; with prices rounded to a tick the effect depends on the tick size and can reverse for a small fall (section 6.3).
* **Better bets help a little, and only if the tails are known.** Bets tuned on Gaussian data (increase bets 0.25 and 0.5; decrease bets at 0.9 and 0.99 of the largest admissible bet) are up to 18% faster than the default grid on Gaussian or jump data, but under Student t3 data they are up to 55% slower than the default, and one decrease case has negative growth. A design pooled over four tail shapes (increase bets 0.16 and 0.39; decrease bets at 0.5 and 0.85 of the largest admissible bet; equal weights) is at least as fast as the default in all 48 validation cases inside the design range, by 0.4% to 10.5% (median 5%), but up to 17% slower for small increases (sd ratio 1.1) under fat tails. The default grid is already within 26% of the best single bet everywhere inside the design range, and within 36% at sd ratio 3, the largest change tried (section 7).
* **The false alarm bound is not loose, the default detector is.** The one-day false alarm rate of the default two-sided detector is about one fifth of the bound `m / A`. One fifth of the gap is overshoot (the statistic jumps past A); the rest is expected wealth that sits on paths that never alarm, mostly held by the decrease side and the smallest increase bets, which cannot build wealth within a day. An adversarial but legitimate null distribution reaches 95% of the bound at one day for a single bet of size 0.5 (a lower bound on the worst case, computed by dynamic programming), so a sharper threshold that is provable for every null distribution does not exist for single bets; for the default mixture it would need a new argument (section 8).
* **Asymptotic optimality (tracker B73).** As A grows, the worst-case delay of a fixed grid divided by `log A / KL_inf` tends to the constant `KL_inf / g_grid` (1.00 to 1.08 for the default grid in our cases), where `log A / KL_inf` is, to first order, the lower bound for any detector that is valid over the whole nonparametric null. So the default grid is first-order optimal only for changes whose best bet is one of its bets. A grid refined with about `√(log A)` bets, or a continuous mixture over the bet, is first-order optimal: we prove that the delay is at most `(log A + ½ log log A + O(1)) / KL_inf` (refined grid) or at most `(log A + log log A + O(1)) / KL_inf` (continuous mixture; `½ log log A` by a heuristic). These are upper bounds; the matching lower bound is proved to first order only. First-order optimality for bounded means is already known (Ram and Ramdas 2026); what is added here is the refinement rate and the constants for this detector, and we make no claim of priority for the first-order result. Numerically the asymptotic regime arrives very slowly: even the oracle bet is 5% to 24% *faster* than `log A / KL_inf` at A = 10^8, and the second-order formula explains this to within 0.024. The robust detector keeps `KL_inf / KL` = 53% to 88% of the asymptotic speed of the Gaussian likelihood ratio detector for increases (sd ratio 2 to 1.25) and 28% to 80% for decreases (sd ratio 0.5 to 0.8) (section 9).

## 1. Set-up and notation

* `X_t = min(r_t², τ²)`, `Y_t = X_t / m₀ ∈ [0, B]` with `B = τ² / m₀` (`B = ∞` without truncation). With `τ = c` times the scaled MAD and Gaussian returns, `τ = c σ₀`, `m₀ = κ₂(c) σ₀²` (`docs/math.md`, section 3) and `B = c² / κ₂(c)`, for example `B = 9.045` for c = 3.
* Write `U = Y - 1` for the increase side and `U = 1 - Y` for the decrease side. U ranges over `[u_lo, u_hi]`, which is `[-1, B - 1]` (increase) or `[1 - B, 1]` (decrease). A bet is `f_λ = 1 + λ U` with `λ ∈ [0, λ_max)`, where `λ_max = 1 / |u_lo|`, that is `λ_max = 1` (increase) and `λ_max = 1 / (B - 1)` (decrease). Then `f_λ > 0`.
* **Null.** `H₀`: the conditional law of `Y_t` given the past lies in `H = {P on [0, B] : E_P[U] ≤ 0}`.
* **Post-change model.** From time ν on, `Y_ν, Y_{ν+1}, ...` are i.i.d. with law Q and independent of the past `F_{ν-1}`, with `E_Q[U] > 0`.
* **Detector.** Bets `λ_1, ..., λ_K` with weights `w_k > 0`, `Σ w_k = 1`. Shiryaev-Roberts (SR) components `M^k_t = f_{λ_k}(Y_t) (M^k_{t-1} + 1)`, `M^k_0 = 0`; mixture `M_t = Σ_k w_k M^k_t`; alarm `T_A = inf{t : M_t ≥ A}`. The CUSUM-type version uses `C^k_t = f_{λ_k}(Y_t) max(C^k_{t-1}, 1)`. `L^{(j,k)}_t = Π_{s=j}^{t} f_{λ_k}(Y_s)` is the wealth of bet k started at time j, and `M^k_t = Σ_{j ≤ t} L^{(j,k)}_t`.
* **Delay.** Lorden's (1971) worst-case delay `WADD(T) = sup_ν ess sup E_ν[(T - ν + 1)⁺ | F_{ν-1}]`, and the conditional delay `E_ν[T - ν + 1 | T ≥ ν]` used in simulations. In the code the first post-change index is `nu` and delay is `alarm - nu + 1`. "Change at the first monitored step" (ν = 1) is the case with no pre-change history; E4 used ν = 501. Lemma 5 (section 9) shows that for the SR mixture the ν = 1 delay *is* the worst-case delay.

## 2. The growth rate

**Definition.** `g(λ) = E_Q[log f_λ(Y)] = E_Q[log(1 + λ U)]`.

**Lemma 0 (shape).** [PROVED, elementary] `g(0) = 0`, `g'(0) = E_Q[U]`, and `g''(λ) = -E_Q[U² / (1 + λU)²] < 0` unless U is degenerate. So g is strictly concave, has a unique maximiser `λ*(Q)` in `[0, λ_max]`, and `λ* > 0` exactly when `E_Q[U] > 0`.

**The oracle bet and its growth.** Write `KL_inf(Q) = max_λ g(λ)` (the name is explained by Lemma 2). `design.best_bet` computes it by numerical integration over the exact law of Y (`design.y_law`: Gauss-Legendre nodes on `|r| < τ` plus the atom at `τ`). The integration agrees with `scipy.integrate.quad` to 10 digits and with 40 million Monte Carlo draws within 0.7 standard errors; `tests/test_delay_theory.py` checks every model and both sides against Monte Carlo.

**The quadratic approximation is biased.** Expanding `log(1 + λU) = λU - λ²U²/2 + λ³U³/3 - ...` and stopping at the square gives `λ ≈ E[U] / E[U²] = d / (v + d²)`, the aGRAPA bet of Waudby-Smith and Ramdas (2024) that the adaptive bettor in `detectors.py` uses. For a truncated squared return, U is strongly skewed to the right on the increase side (`E[U³]` is large), so the cubic term matters. **[NUMERICAL]** Under Gaussian data with c = 3:

| sd ratio | E_Q[Y] | oracle bet λ* | KL_inf = g(λ*) | aGRAPA bet | its growth / KL_inf | best default-grid bet | its growth / KL_inf |
|---|---|---|---|---|---|---|---|
| 1.1 | 1.202 | 0.102 | 0.0093 | 0.074 | 0.94 | 0.10 | 1.00 |
| 1.25 | 1.524 | 0.241 | 0.0514 | 0.123 | 0.80 | 0.20 | 0.98 |
| 1.5 | 2.082 | 0.423 | 0.1645 | 0.145 | 0.64 | 0.35 | 0.98 |
| 2 | 3.129 | 0.643 | 0.4302 | 0.145 | 0.49 | 0.50 | 0.96 |

Over the ten increase cases of `results/e11a_delay.csv` (three models at sd ratios 1.25, 1.5 and 2, plus Gaussian at 1.1), the plug-in aGRAPA bet is too small by a factor of 1.4 to 4.4 and keeps only 49% to 94% of the growth (1.7 to 4.4 and 49% to 86% without the Gaussian case at 1.1), even when the mean and variance are known exactly. This is one reason why the quarticity-based adaptive bet did not help in E4 (Q4c), on top of the noise in estimating them and the cap at half the largest bet. For decreases the bias has the opposite sign: `U = 1 - Y` has a long left tail, the quadratic rule asks for a bet beyond `λ_max` and is clipped to `0.999 λ_max`. On Gaussian data at sd ratios 0.5 to 0.8 this costs at most 1%, because the oracle bet sits at the edge anyway (0.992 to 0.99999 of `λ_max`); but at sd ratio 0.9 it keeps only 39% of the growth, with the jump model at sd ratio 0.8 47%, and under Student t5 data at sd ratio 0.8 its growth is negative. A bet at the edge loses `log(1 - λ/λ_max)` at every capped return. The fixed default grid reaches 96% to 100% of the optimal growth for increases in this range.

**Growth of a mixture.** [STANDARD] With fixed weights, `(1/n) log M_n → g_grid := max_k g(λ_k)` almost surely on post-change data (strong law of large numbers for each `log L^{(j,k)}_n`, together with `max_k (log w_k + log L^{(1,k)}_n) ≤ log M_n ≤ log n + max_{j,k} log L^{(j,k)}_n`: the SR sum also contains products that start after step 1, so the upper bound must take the maximum over start times j too; the proof of Theorem 3a below gives the details).

## 3. What limits the growth

**Lemma 1 (growth is bounded by information).** [PROVED] Let P be a law of Y and `E ≥ 0` a function with `E_P[E(Y)] ≤ 1` (an e-variable for P). Then `E_Q[log E(Y)] ≤ KL(Q ‖ P)`.

*Proof.* If Q is not absolutely continuous with respect to P, `KL = ∞`. Otherwise let `q = dQ/dP`. On `{q > 0}`, `log E - log q = log(E / q)`, so by Jensen's inequality `E_Q[log E] - KL(Q ‖ P) = E_Q[log(E / q)] ≤ log E_Q[E / q] = log E_P[E · 1{q > 0}] ≤ log E_P[E] ≤ 0`. ∎

Every bet `f_λ` is an e-variable for every P in H, in particular for the true pre-change law `P₀` of Y. Since Y is a function of the return r, the data processing inequality gives the chain

`g(λ) ≤ KL(Q_Y ‖ P₀_Y) ≤ KL(Q_r ‖ P₀_r)`.

**Lemma 2 (the best nonparametric growth, KL_inf duality).** [PROVED] Let B be finite, `KL_inf(Q) := inf_{P ∈ H} KL(Q ‖ P)`.

(a) Every e-variable for the whole class H (`E ≥ 0` and `sup_{P ∈ H} E_P[E] ≤ 1`) is dominated by a bet: there is `λ ∈ [0, λ_max]` with `E(y) ≤ 1 + λ u(y)` for all y.

(b) `max_λ g(λ) = KL_inf(Q)` whenever the maximiser `λ*` is interior (`0 < λ* < λ_max`).

*Proof.* (a) Work with u = U(y), which ranges over `[u_lo, u_hi]`, with `u_lo < 0 < u_hi`. Let φ be the upper concave envelope of E as a function of u. The largest value of `E_P[E]` over laws with mean `E_P[U] = a` equals φ(a) (two-point laws suffice), so the assumption says `φ(a) ≤ 1` for every `a ≤ 0`. φ is concave and finite and 0 is an interior point, so it has a supergradient s at 0: `φ(u) ≤ φ(0) + s u ≤ 1 + s u`. If `s ≥ 0`, take `λ = s`. Since `0 ≤ E(u_lo) ≤ φ(u_lo) ≤ 1 + s u_lo`, we get `s ≤ 1 / |u_lo| = λ_max`. If `s < 0`, φ decreases on `[0, u_hi]`, so `φ(u) ≤ φ(0) ≤ 1` there, and `φ(u) ≤ 1` for `u ≤ 0` by assumption; take `λ = 0`.

(b) "≤": each `f_λ` is an e-variable for each P in H, so Lemma 1 gives `g(λ) ≤ KL(Q ‖ P)`; take the infimum over P and the maximum over λ. "≥": at an interior maximiser, `g'(λ*) = E_Q[U / (1 + λ* U)] = 0`. Define `P*(dy) = Q(dy) / (1 + λ* U(y))`. Its total mass is `E_Q[1 / (1 + λ*U)] = 1 - λ* E_Q[U / (1 + λ*U)] = 1`, its mean of U is `E_Q[U / (1 + λ*U)] = 0`, so `P* ∈ H`, and `KL(Q ‖ P*) = E_Q[log(1 + λ*U)] = g(λ*)`. ∎

Part (b) is the dual form of the "K_inf" information number of Honda and Takemura (2010) for bounded support, and P* is the reverse information projection of Q onto the null (Grünwald, de Heide and Koolen 2024; Larsson, Ramdas and Ruf 2025, who show in general that the growth-optimal e-variable is the likelihood ratio against this projection). `tests/test_delay_theory.py::test_kl_inf_duality_certificate` checks numerically that P* has mass 1, mean exactly on the null boundary and `KL(Q ‖ P*) = g(λ*)`. In the boundary case `λ* = λ_max` the same holds with the leftover mass of P* placed at `u_lo`; we do not need it.

**Corollary 3 (the information ladder).** [PROVED] With `g_grid = max_k g(λ_k)`:

`KL(Q_r ‖ P₀_r) ≥ KL(Q_X ‖ P₀_X) ≥ KL_inf(Q) ≥ g_grid`.

The three ratios are the **truncation loss** `KL_r / KL_X`, the **nonparametric loss** `KL_X / KL_inf` (the price of a null that constrains only the mean of Y, not its shape) and the **grid loss** `KL_inf / g_grid`. `test_information_ladder_and_gaussian_kl` checks the ladder in 36 cases.

**Corollary 4 (speed limits).** [PROVED] Since `f_λ ≤ 1 + λ_max u_hi`, every bet and therefore every nonparametric e-variable (Lemma 2a) satisfies `log E ≤ log B` (increase) and `log E ≤ log(B / (B - 1))` (decrease). Hence `KL_inf ≤ log B` for increases and `KL_inf ≤ log(B / (B - 1)) ≈ 1 / (B - 1/2)` for decreases, however large the change.

For a Shiryaev-Roberts detector the cap gives a floor on the delay that holds path by path.

* *A given grid.* If every step earns the largest factor `F_k = 1 + λ_k u_hi`, then `M_n ≤ Σ_k w_k Σ_{j=1}^{n} F_k^{n-j+1}`, so `T_A` can never be smaller than the first n where this sum reaches A (`design.speed_limit`).
* *Every valid bet.* On the decrease side every valid factor is at most `F = B / (B - 1)`. So every SR e-detector whose factors are valid e-variables for the whole null, with any bets and any weights (fixed, or chosen from past data), satisfies `M_n ≤ Σ_{j=1}^{n} F^{n-j+1} = F (F^n - 1) / (F - 1) = B (F^n - 1)`, because `(F - 1) / F = 1 / B`. It cannot alarm before `F^n ≥ 1 + A / B`, that is before `log(1 + A / B) / log(B / (B - 1))` steps (`design.decrease_speed_limit`). A single bet close to `λ_max` reaches this floor when every return is zero.

For Gaussian returns:

| c | B | cap per step, `log(B / (B - 1))` nats | floor for any valid bets, A = 7800 (2A = 15,600) | default five-bet grid, A = 7800 (2A) | `log A / log(B / (B - 1))` |
|---|---|---|---|---|---|
| 1.5 | 2.890 | 0.425 | 19 (21) | 23 (25) | 21.1 |
| 2 | 4.345 | 0.262 | 29 (32) | 36 (39) | 34.3 |
| 3 | 9.045 | 0.117 | 58 (64) | 73 (79) | 76.5 |
| 4 | 16.00 | 0.0645 | 96 (107) | 122 (134) | 138.9 |
| 5 | 25.00 | 0.0408 | 141 (158) | 180 (199) | 219.5 |

At c = 3, for example, `log(1 + 7800 / 9.045) / 0.11716 = 57.7`, so 58 steps. The default grid is slower because it gives weight 1/5 to each of five bets (0.1, 0.25, 0.5, 0.75 and 0.95 of `λ_max`), and the smaller ones earn less than the largest factor at every step. A volatility collapse to zero returns is detected after exactly 73 steps by the default grid at c = 3, and after 58 steps by a single bet at `0.999 λ_max` (`test_speed_limit_is_attained_by_extreme_data`, `tests/test_decrease_speed_limit.py`). The last column ignores the sum over start times: it is the first-order delay `log A / KL_inf` of a total collapse (section 4.2), a limit as A grows and not a floor at A = 7800. For detectors that are not e-detectors, only that asymptotic lower bound is proved here (Theorem 2b).

## 4. Detection delay

### 4.1 A non-asymptotic bound

**Theorem 1.** [PROVED] Under the post-change model, for the SR or the CUSUM-type mixture, every alarm level A, every change time ν and every pre-change history,

`E[(T_A - ν + 1)⁺ | F_{ν-1}] ≤ min over k with g_k > 0 of (log(A / w_k) + o_k) / g_k`,

where `g_k = g(λ_k)`, `Z_k = log f_{λ_k}(Y)` with Y from Q, and `o_k = min(sup Z_k, E[(Z_k⁺)²] / E[Z_k])` bounds the mean overshoot.

*Proof.* Fix k. For `t ≥ ν`, all terms of the SR sum are nonnegative, so `M_t ≥ w_k M^k_t ≥ w_k L^{(ν,k)}_t`; the CUSUM-type statistic also satisfies `C^k_t ≥ L^{(ν,k)}_t`. Let `S_n = Σ_{i=0}^{n-1} Z_{k, ν+i}` and `N = inf{n ≥ 1 : S_n ≥ log(A / w_k)}`. Then `w_k L^{(ν,k)}_{ν+N-1} = w_k e^{S_N} ≥ A`, so either `T_A < ν` or `T_A ≤ ν + N - 1`; in both cases `(T_A - ν + 1)⁺ ≤ N`. N depends only on post-change data, which are independent of `F_{ν-1}`. S is a random walk with drift `g_k > 0` and increments bounded below by `log(1 - λ_k |u_lo|) > -∞`, so `E[N] < ∞` and Wald's identity gives `g_k E[N] = E[S_N] = log(A / w_k) + E[S_N - log(A / w_k)]`. The overshoot is at most `sup Z_k`, and its mean is at most `E[(Z⁺)²] / E[Z]` uniformly in the level by Lorden's (1970) inequality. ∎

The bound is the elementary version, for linear bets on a bounded statistic, of the delay bounds that Shin, Ramdas and Rinaldo (2024, section 4) give for mixtures of exponential baseline e-detectors. It holds for every A, including the two-sided use with each side at 2A (replace A by 2A). It pays the full `log(1 / w_k) = log K` for not knowing the right bet and ignores the help of the other bets and of the other start times; section 5 shows it is 1.3 to 2.7 times the simulated delay.

### 4.2 First-order asymptotics and optimality

**Theorem 2.** (a) [PROVED, section 9] For a fixed grid with `g_grid > 0`, `WADD(T_A) = (log A / g_grid)(1 + o(1))` as `A → ∞` (Theorem 3a).

(b) [PROVED given Lorden 1971] Let T be any stopping rule (not necessarily an e-detector) whose average run length is at least A under every i.i.d. pre-change law in H. Then `WADD_Q(T) ≥ (log A / KL_inf(Q))(1 + o(1))`. *Proof:* the infimum in Lemma 2 is attained by `P* ∈ H` (section 3), the rule has average run length at least A under P*, and Lorden's (1971) lower bound for detecting a change from P* to Q gives `(log A / KL(Q ‖ P*))(1 + o(1))`. ∎ Ram and Ramdas (2026) prove a universal lower bound of this form with `KL_inf` for composite pre- and post-change classes, and a matching upper bound for bounded means, which covers our truncated statistic. So first-order optimality of a suitable mixture is already known for our setting; Theorems 2b and 3 are instances of it with explicit constants for this detector.

(c) Consequently the SR e-detector with the oracle bet is first-order optimal among all detectors that control the run length over the nonparametric null (its average run length is at least A under every null law, `docs/math.md` section 5), and a grid detector has asymptotic efficiency `g_grid / KL_inf(Q)`. Against the Gaussian likelihood ratio CUSUM or SR, which are first-order optimal when returns are exactly normal (Lorden 1971; Pollak 1985), the asymptotic delay ratio is `KL(Q_r ‖ P₀_r) / g_grid`. Section 9 states how the grid must be refined to recover full efficiency.

### 4.3 Second order: why `log(A) / g` is not enough at A = 7800

[HEURISTIC; accuracy checked in sections 5 and 9] Take one bet and the change at the first step (`M_0 = 0`). Because `M_n = Σ_{j ≤ n} L^{(j)}_n = e^{S_n} V_n` with `S_n = Σ_{i ≤ n} Z_i` and `V_n = Σ_{j=0}^{n-1} e^{-S_j}`,

`log M_n = S_n + log V_n`, and `V_n` increases to `V_∞ = Σ_{j ≥ 0} e^{-S_j} < ∞`.

So `log M_n` is a random walk plus a slowly changing term. Nonlinear renewal theory (Lai and Siegmund 1977; Woodroofe 1982; Pollak 1987 for the SR statistic) then gives

`E[T_A] = (log A + κ - ζ) / g + o(1)`, with `κ = E[H²] / (2 E[H])` and `ζ = E[log V_∞]`,

where H is the first ascending ladder height of S (κ is the limiting mean overshoot). We have not checked the regularity conditions of those theorems for our law of Y (continuous part plus an atom at the cap), so we label the formula heuristic. `design.renewal_constants` estimates κ and ζ by Monte Carlo.

**A closed form for the memory term.** Approximate S by a Brownian motion with drift g and variance `s²` per step. Dufresne's (1990) identity, `∫_0^∞ e^{2(W_u - a u)} du` has the law of `1 / (2 Γ_a)` with `Γ_a ~ Gamma(a, 1)`, gives after a time change

`V_∞ ≈ 2 / (s² Γ_{2g/s²})`, so `ζ ≈ log(2 / s²) - ψ(2g / s²)`

(ψ is the digamma function; `design.zeta_brownian`). For a growth-optimal bet and a small change, `s² ≈ 2g`, and then `ζ ≈ log(1 / g) + γ` with Euler's constant `γ = 0.5772`. Therefore

**`E[delay] ≈ (log(A g) - γ + κ) / g`.**

The memory term is the reason the first-order formula fails: for a 25% increase in standard deviation `g ≈ 0.05`, so `ζ ≈ 3.6` while `log A = 8.96`; the first-order formula overstates the delay by about 60%. The relative error of `log(A)/g` is roughly `(log(1/g) + γ - κ) / log(A g)`, which grows as the change gets smaller and only disappears as `A → ∞` (very slowly, since it decays like `1 / log A`). The Brownian ζ agreed with the Monte Carlo ζ to within 0.08 for every case with `g < 0.2`, and within 0.21 for the large changes (`g` from 0.31 to 0.43) (`results/e11a_delay.csv`, columns `zeta_oracle` and `zeta_oracle_bm`).

**Head start.** When the change happens after a pre-change period (ν = 501 in E4), the SR statistic enters it with `M_{ν-1} > 0` (under the null `E[M_{ν-1}] ≤ ν - 1`), and `log M_n = S_n + log(M_{ν-1} + V_n)`. This shortens the delay further: in the simulations by 17% to 34% compared with ν = 1. Under the null, small bets keep `M^k` close to its mean `ν - 1` on most paths, while for large bets the mean is carried by rare paths and `M^k` is usually small, so the head start helps grids with small bets more (section 7).

**Mixtures.** Pathwise, `T_A` is at most the first time any single component k reaches `A / w_k`, so `(log(A / w_k) + κ_k - ζ_k) / g_k` for the best k approximates an upper bound. The real mixture is faster, because on each path the bet that happens to fit best takes the lead (the components are driven by the same data), and the overhead is less than `log K`. We have no closed form for this effect; the simulations measure it.

## 5. How good are the approximations?

[NUMERICAL] `results/e11a_delay.csv`. True baseline (population τ and m₀), c = 3, A = 7800, default grid (six increase bets, five decrease bets), change at the first monitored step, 3000 paths per case (2000 for ν = 501). "Second order, oracle" uses Monte Carlo κ and ζ and is compared with a simulation of the single oracle bet.

| Scenario | sd ratio | g_grid | first order log(A)/g_grid | Wald bound (proved) | simulated grid, ν = 1 | second order, oracle | simulated oracle, ν = 1 | simulated grid, ν = 501 | E4 (estimated baseline, ν = 501) |
|---|---|---|---|---|---|---|---|---|---|
| Gaussian | 1.1 | 0.0093 | 964 | 1221 | 455 | 405 | 417 | 378 | |
| Gaussian | 1.25 | 0.0502 | 179 | 233 | 122 | 109 | 107 | 98 | 100 |
| Gaussian | 1.5 | 0.161 | 56 | 75 | 48 | 42 | 42 | 39 | 40 |
| Gaussian | 2 | 0.413 | 22 | 30 | 23 | 18.4 | 18.7 | 18 | 18 |
| Jumps | 1.25 | 0.0473 | 190 | 247 | 128 | 115 | 115 | 105 | 114 |
| Jumps | 2 | 0.401 | 22 | 31 | 23 | 18.9 | 19.2 | 18 | 19 |
| Student t5 | 1.25 | 0.0320 | 280 | 362 | 178 | 158 | 159 | 144 | 141 |
| Student t5 | 2 | 0.309 | 29 | 39 | 29 | 24.1 | 24.0 | 24 | 25 |
| Gaussian | 0.9 (down) | 0.0101 | 886 | 1054 | 434 | 389 | 392 | 342 | |
| Gaussian | 0.8 (down) | 0.0340 | 264 | 314 | 194 | 157 | 156 | 143 | 151 |
| Gaussian | 0.67 (down) | 0.0598 | 150 | 179 | 124 | 99 | 99 | 86 | 91 |
| Gaussian | 0.5 (down) | 0.0839 | 107 | 127 | 93 | 74 | 74 | 62 | 67 |

(Mean delays in steps. Standard errors are below 1.5% of the mean. The full table has all 20 cases.)

**Reading.**

* **Second order works.** Across all 20 cases the second-order prediction for a single bet is within 3% of its simulated mean delay (ratios 0.97 to 1.02).
* **First order fails for small changes.** Across the 20 cases `log(A) / g_grid` divided by the simulated grid delay ranges from 0.94 to 2.12, and it exceeds 1.35 in every small-change case (sd ratios 1.1, 1.25, 0.9 and 0.8). For large changes it is close (0.94 to 0.99 at sd ratio 2), but only because two errors cancel: it ignores the memory term (which shortens the delay) and the grid overhead `log K` and the overshoot (which lengthen it).
* **The proved bound holds everywhere** and is 1.3 (large changes) to 2.7 (small changes) times the simulated delay.
* **The grid costs 9% to 26% against the oracle bet** (ν = 1), less than its `log K` overhead suggests, and the dominant-bet second-order formula with weight 1/K overstates the grid delay by 3% to 32%: the other bets help.
* **E4 is reproduced.** With the true baseline and ν = 501 the simulated delays match E4's (which used estimated baselines and only 300 paths) to within about 10%. The pre-change head start explains why the ν = 501 delays, and E4's, are 17% to 34% shorter than the ν = 1 values.

## 6. The price of robustness

[PROVED decomposition, NUMERICAL values] `results/e11b_information.csv`, `e11b_price_simulated.csv`, `e11b_c_curve.csv`, `e11b_jump_burst.csv`.

### 6.1 Where the information goes (clean Gaussian data, c = 3)

| sd ratio | KL (Gaussian LR) | KL of truncated X | KL_inf (best robust bet) | g_grid | truncation loss | nonparametric loss | grid loss | total KL / g_grid |
|---|---|---|---|---|---|---|---|---|
| 1.1 | 0.0097 | 0.0096 | 0.0093 | 0.0093 | 1.01 | 1.03 | 1.00 | 1.04 |
| 1.25 | 0.0581 | 0.0565 | 0.0514 | 0.0502 | 1.03 | 1.10 | 1.02 | 1.16 |
| 1.5 | 0.220 | 0.203 | 0.165 | 0.161 | 1.08 | 1.23 | 1.02 | 1.37 |
| 2 | 0.807 | 0.638 | 0.430 | 0.413 | 1.27 | 1.48 | 1.04 | 1.95 |
| 3 | 2.90 | 1.56 | 0.861 | 0.747 | 1.86 | 1.81 | 1.15 | 3.89 |
| 0.9 (down) | 0.0104 | 0.0103 | 0.0101 | 0.0101 | 1.00 | 1.02 | 1.00 | 1.02 |
| 0.8 (down) | 0.0431 | 0.0431 | 0.0346 | 0.0340 | 1.00 | 1.24 | 1.02 | 1.27 |
| 0.67 (down) | 0.125 | 0.125 | 0.063 | 0.060 | 1.00 | 2.00 | 1.05 | 2.09 |
| 0.5 (down) | 0.318 | 0.318 | 0.088 | 0.084 | 1.00 | 3.61 | 1.05 | 3.79 |

**The largest cost is not truncation.** For increases up to sd ratio 2 the nonparametric loss (betting only on the mean of Y, which is what makes the detector valid for any return distribution) is larger than the truncation loss; the grid costs 2% to 4%. For decreases the nonparametric loss is the speed limit of Corollary 4: KL_inf cannot exceed `log(B / (B - 1)) = 0.117` per step, while the Gaussian KL keeps growing.

### 6.2 Does `KL / g` explain the E4 penalty?

[NUMERICAL] Mean delay, change at the first step, 3000 paths with common random numbers, c = 3:

| | sd 1.25 | sd 1.5 | sd 2 | sd 0.8 | sd 0.67 | sd 0.5 |
|---|---|---|---|---|---|---|
| Gaussian LR, true ratio known (oracle) | 99.8 | 34.3 | 12.4 | 126.1 | 53.0 | 24.1 |
| Gaussian LR mixture of E4 (three ratios) | 105.8 | 35.6 | 13.3 | 130.7 | 56.5 | 27.0 |
| Robust, oracle bet | 108.6 | 41.6 | 18.5 | 156.0 | 98.7 | 74.2 |
| Robust, default grid | 121.1 | 48.1 | 23.0 | 193.9 | 124.0 | 93.1 |
| **Penalty: grid / LR mixture** | **1.14** | **1.35** | **1.73** | 1.48 | 2.19 | 3.45 |
| Prediction `KL / g_grid` (first order) | 1.16 | 1.37 | 1.95 | 1.27 | 2.09 | 3.79 |
| Penalty: oracle bet / LR oracle | 1.09 | 1.21 | 1.48 | 1.24 | 1.86 | 3.08 |
| Prediction, second-order delays | 1.10 | 1.23 | 1.52 | 1.24 | 1.86 | 3.07 |

E4 reported median penalties of about 1.12, 1.35 and 1.7 for increases (robust grid against the Gaussian LR mixture). **Yes, the information ratio explains them:** `KL / g_grid` predicts 1.16, 1.37 and 1.95, and the second-order formula, applied to both detectors, predicts the oracle-against-oracle penalty to within 0.04 in all six cases (the largest difference, 0.035, is at sd ratio 2). The first-order ratio overstates the penalty for large changes because the memory and overshoot terms are of similar size for both detectors and dilute the ratio.

### 6.3 The truncation constant c

[NUMERICAL] Growth relative to the Gaussian KL (efficiency `g_grid / KL`) and simulated mean delays, default grid:

| c | clipped share of normal returns | increase: delay at sd 1.25 / 1.5 / 2 (clean) | growth at sd 1.25 with jumps (E4 jump model) | jump burst growth / growth of a real 25% increase | consecutive capped returns needed to alarm from zero | decrease: delay at sd 0.8 / 0.67 / 0.5 (clean) | decrease speed limit in steps, default grid (any valid bets) |
|---|---|---|---|---|---|---|---|
| 2 | 4.6% | 143 / 59 / 30 | 0.0400 | 0.013 | 11 | 143 / 70 / 47 | 36 (29) |
| 3 | 0.27% | 121 / 48 / 23 | 0.0473 | 0.035 | 7 | 194 / 124 / 93 | 73 (58) |
| 4 | 0.006% | 118 / 46 / 21 | 0.0458 | 0.074 | 5 | 299 / 203 / 156 | 122 (96) |
| 5 | 0.00006% | 117 / 45 / 21 | 0.0420 | 0.112 | 5 | 427 / 297 / 229 | 180 (141) |
| none | 0 | 117 / 45 / 21 | 0.0135 | 0.39 | 1 | not possible | |

"Jump burst": the standard deviation stays the same but the jump probability rises from 0.005 to 0.02 per step (jumps of 8 standard deviations); the column gives the growth this produces for the increase detector, relative to the growth produced by a genuine 25% rise in standard deviation. It measures how strongly jumps alone are mistaken for higher volatility.

**Reading and a recommended balance.**

* **Increases.** On clean data the delay barely changes for `c ≥ 4`, and c = 3 costs 3%, 6% and 12% at sd ratios 1.25, 1.5 and 2 against no truncation. With jumps in the data, growth at sd ratio 1.25 peaks at c between 3 and 3.5 (`e11b_c_curve.csv`: best c = 3.25 at sd 1.25, 3.5 at 1.5 and 4.25 at 2); without truncation it collapses to 29% of the c = 3 value (E4: delay 246 against 102). A jump burst moves the c = 3 detector only 3.5% as much as a real 25% volatility rise, against 7% at c = 4 and 39% without truncation. **c = 3 to 4 balances the two**; c = 3 if jumps are frequent, c = 4 if delay on clean data matters more.
* **Decreases.** Here truncation is what makes detection possible at all, and the speed limit dominates: smaller c is faster. The growth-maximising c is 1.25 for sd ratio 0.5, 1.5 for 0.67, 2 for 0.8 and about 3 for 0.9 (`e11b_c_curve.csv`); c = 1.5 keeps about 80% or more of the Gaussian information for all ratios from 0.5 to 0.9 (efficiency 0.84, 0.94, 0.88, 0.80), c = 2 at least 57%. With c = 2 the simulated decrease delays fall by 26% to 49% compared with c = 3.
* **The two sides need not share c.** Each side is an e-detector for its own truncated statistic, and the two-sided rule only needs each side's guarantee, so `c ≈ 3.5` for increases and `c ≈ 1.5 to 2` for decreases is a legitimate design. The guarantee needs each side's baseline to be known, or to come from the confidence-sequence baseline of `docs/math.md` section 6: `VolatilityMonitor(baseline="cs", c_down=...)` runs the lower confidence sequence on `min(r², τ_dn²)` and the upper one on `min(r², τ²)`, each at level δ/2, so Theorem 6.4 of `docs/math.md` applies to each side with the same total δ. With the plug-in baseline there is no proof, with or without the split. The null of the decrease side then reads "the mean of returns squared and clipped at `τ_dn` (1.5 to 2 MADs) has not fallen". For a fixed shape of the innovations this is implied by "volatility has not fallen" (`docs/math.md`, 5.1), but a smaller clip is more sensitive to changes of shape at constant volatility, for example a larger share of zero returns. E11f tested the split with the plug-in baseline (`docs/report_parts/E11.md`, E11f). For three continuous laws (Gaussian, Student t4, jumps) c_down = 2 cut the median decrease delay in all nine cells, by 13% to 51%. With prices rounded to a tick that is fixed in units of the pre-change standard deviation, the result depends on the tick size: c_down = 2 was 32% to 52% faster at ticks of 0.25 and 0.75 sd, but 12% slower for a 20% fall at a tick of 0.5 sd (median 211 against 188 steps), where the calibration scale is one tick times 1.4826 = 0.74 sd and c = 3 already clips at about 2.2 sd. On five coins of real data, which have few ties, c_down = 2 was neutral to slightly positive. The no-change false alarm rates of these plug-in runs were 5.0% to 6.7% per day, and 5.0% to 6.05% with c_down = c in the same runs, so the excess over the 5% bound (which holds with the true baseline) comes from the estimated baseline, not from the split.

### 6.4 Plain-language summary

On clean normal data the robust detector with the default bets takes 14% to 73% longer than the Gaussian likelihood ratio detector of E4 to detect an increase of 25% to 100% in standard deviation (section 6.2). The first-order prediction of that gap, `KL / g_grid` (section 6.1), splits on a log scale as follows: about 60% (59% to 68%) is the cost of making no assumption about the shape of the return distribution, 19% to 35% is the cost of clipping large returns, and the remaining 6% to 16% is the cost of not knowing the size of the change in advance (sd ratios 1.25 to 2). This splits the prediction, not the simulated penalty. In exchange, the robust detector keeps its false alarm guarantee when returns have jumps or fat tails, where the parametric detector raises a false alarm on almost every path (E4). For decreases the price is higher, because clipping caps how much evidence one step can provide. A tighter clip on the decrease side recovers part of it on continuous data (E11f: a quarter to a half of the delay for Gaussian data); with prices rounded to a coarse tick the gain depends on the tick size and can turn into a loss for small falls.

## 7. Choosing the bets

[NUMERICAL] `results/e11c_design.csv`, `e11c_design_sim.csv`, `e11c_null_check.csv`, `e11c2_variants.csv`, `e11c3_final_validation.csv`, `e11c4_pooled_design.csv`, `e11c4_robust_validation.csv`; part E of `e11e_asymptotic_optimality.csv`.

### 7.1 Objective and method

For a prior range of change sizes (sd ratio 1.2 to 2 for increases, 0.5 to 0.85 for decreases), choose bets and weights that minimise the **worst-case ratio of the mixture's delay to the oracle delay** (the delay of the single best bet for that change). The design uses formulas only: each bet's delay is approximated by the second-order SR formula with the Brownian memory term and weight cost `log(1 / w_k)` (`design.sr_delay_approx`), and the mixture by its best bet (`design.design_ratio`). A lattice search over geometric grids is followed by a Nelder-Mead polish (`design.optimise_grid`). Candidates are then selected by simulation on one seed and validated on fresh seeds. The decrease bets are expressed as fractions of `λ_max`, because `λ_max = 1/(B - 1)` depends on the tail shape through B.

This objective is exactly the finite-A trade-off of Theorem 3b: each extra bet costs `log K` nats of evidence (a second-order term), while gaps in the grid cost growth, which is multiplied by `log A`. At A = 7800 the trade-off favours two or three bets.

### 7.2 Design on Gaussian data (parts C, C2, C3)

The Gaussian design selected two increase bets (0.2565 and 0.4893 with weights 0.569 and 0.431; predicted worst ratio 1.18 against 1.37 for the default) and a single decrease bet at `0.99 λ_max` (predicted 1.01 against 1.33). The rounded increase grid (0.25, 0.5) with equal weights did as well in simulation. The single decrease bet at `0.99 λ_max` turned out to be fragile under Student t5 data (ratio to the oracle 1.12 at sd ratio 0.8, 1.26 at 0.85 and 1.51 at 0.9, `e11c2_variants.csv`), because a bet at the edge loses `log(1 - 0.99) = -4.6` nats at every capped return; two bets at 0.9 and 0.99 of `λ_max` cut the worst in-range ratio over the three models to 1.12.

Fresh-seed validation of this "Gaussian-tuned" pair (3000 paths, `e11c3_final_validation.csv`): with the change at the first step it is 7% to 10% faster than the default for increases and 17% faster for decreases on Gaussian data, and similar on the jump model. With the change at step 501 the gains shrink (3% to 8% for increases on Gaussian data), and under t5 data at sd ratios 1.2 and 0.85 the tuned grid is 12% and 15% *slower* than the default. Part E then showed a worse failure under Student t3 data: at sd ratio 0.8 both decrease bets have negative growth, 22% of paths had not alarmed at A = 10^8, and for increases at sd ratio 1.25 its asymptotic delay ratio is 1.43.

### 7.3 Design pooled over tail shapes (part C4) and its validation (part C4V)

Because the Gaussian-tuned grid had already been seen to fail, a second design was made before any further simulation, pooling 36 post-change laws for increases (Gaussian, jumps, Student t5, Student t3; 9 sd ratios each) and 32 for decreases (8 each). The selection rule was fixed in advance: the smallest predicted pooled worst ratio, ties within 1% going to fewer bets, with optimised bets rounded to two digits (increase) or to multiples of 0.05 (decrease fractions). Predicted pooled worst ratios:

| Design | Increase bets | Predicted worst ratio, pooled | Decrease bets (fraction of λ_max) | Predicted worst ratio, pooled |
|---|---|---|---|---|
| default | 0.02, 0.05, 0.1, 0.2, 0.35, 0.5 | 1.54 | 0.1, 0.25, 0.5, 0.75, 0.95 | 1.44 |
| Gaussian-tuned | 0.25, 0.5 | 4.13 | 0.9, 0.99 | no finite value (negative growth) |
| **robust (selected)** | **0.16, 0.39** | **1.27** | **0.5, 0.85** | **1.26** |
| best three bets | 0.13, 0.23, 0.44 | 1.27 | 0.45, 0.65, 0.9 | 1.29 |

Fresh-seed validation, 2000 paths per case, four tail shapes, change at step 1 and at step 501, ratio of mean delay to the oracle bet's mean delay (`e11c4_robust_validation.csv`; worst case over sd ratios inside the design range and over both change times):

| Model | Increase: default | Increase: Gaussian-tuned | Increase: robust | Decrease: default | Decrease: Gaussian-tuned | Decrease: robust |
|---|---|---|---|---|---|---|
| Gaussian | 1.25 | 1.13 | 1.24 | 1.26 | 1.04 | 1.25 |
| Jumps | 1.25 | 1.14 | 1.24 | 1.21 | 1.06 | 1.18 |
| Student t5 | 1.21 | 1.24 | 1.19 | 1.23 | 1.34 | 1.21 |
| Student t3 | 1.19 | 1.39 | 1.17 | 1.17 | 1.78 | 1.13 |

**Reading.**

* The robust grid is never slower than the default inside the design range: across all 48 in-range cases (4 models, 3 sd ratios per side, 2 change times) its delay is 0.895 to 0.996 times the default's (median 0.948). The Gaussian-tuned grid ranges from 0.82 to 1.55 times the default (median 0.91).
* **The formula overstated the gain.** It predicted a worst-case ratio of 1.54 for the default against 1.27 for the robust grid; in simulation both are about 1.25. The dominant-bet approximation is accurate for two bets but too pessimistic for the six-bet default, whose other bets help (section 5).
* **Outside the design range the robust grid can be slower.** For a small increase (sd ratio 1.1) under t3 data it is 13% (change at step 1) and 17% (change at step 501) slower than the default, and 14% slower under t5 data at step 501, because its smallest bet (0.16) is too large for that change. For sd ratio 3, sd ratio 0.4 and sd ratio 0.9 it is within 0.90 to 1.04 of the default.
* **Head start.** With the change at step 501 the default gains more than the two-bet grids (section 4.3): its small bets carry wealth close to `ν - 1` into the change. This is why the gains of both designed grids shrink at ν = 501.
* **False alarms before the change** (500 null steps, bound 500/7800 = 6.4%): default 1.0% to 3.2%, robust 1.6% to 4.6%, Gaussian-tuned 2.9% to 5.5%. All are within the bound; the designed grids use more of it because they waste less wealth (section 8).
* **Larger A favours more bets.** At A = 10^8 the six-bet default is faster than the robust grid in three of the five Gaussian cases of part E, as Theorem 3 predicts: the asymptotic delay ratio of the robust grid (`KL_inf / g_grid`) is 1.10 to 1.17 in those cases, against 1.02 to 1.05 for the default.

### 7.4 Recommendation

* For alarm levels around `10^4` and changes of sd ratio 1.2 to 2 (up) or 0.5 to 0.85 (down), **increase bets (0.16, 0.39) and decrease bets (0.5, 0.85) × `λ_max`, equal weights**, are at least as fast as the default for every tail shape tried, by up to 10%. The gain is modest; the default grid is already a sound robust choice.
* If smaller increases (sd ratio below 1.2) matter, or the alarm level is much larger than `10^4`, keep the default grid (or add a small bet to the robust grid; not tested).
* Do not use the Gaussian-tuned bets unless the tails of the monitored returns are known to be close to normal: their extra speed on Gaussian data (up to 18% against the default) is not worth losses of up to 55% under t3 tails.
* The tail shape matters through `λ_max` (decrease side) and through the size of the oracle bet: under t3 data the oracle increase bets are 20% to 37% smaller than under Gaussian data (0.15 against 0.24 at sd ratio 1.25, 0.52 against 0.64 at sd ratio 2; `e11e_asymptotic_optimality.csv`, column `lam_oracle`).

## 8. Why the false alarm bound looks loose

[PROVED identity, NUMERICAL values] `results/e11d_decomposition.csv`, `e11d_two_sided_split.csv`, `e11d_worst_case_dp.csv`, `e11d_two_point_adversary.csv`; figure `results/figures/e11_false_alarm_slack.png`.

### 8.1 An exact identity

Let `σ = T_A` and m a horizon. Under any null law in H, `M_t - t` is a supermartingale with `M_0 = 0` (`docs/math.md`, section 5), so optional stopping at `σ ∧ m` gives

`P(σ ≤ m) E[M_σ | σ ≤ m] + E[M_m ; σ > m] ≤ E[σ ∧ m]`,

with equality under an exact null (`E_P[U] = 0` at every step, as for Gaussian data with the true baseline). [PROVED] Under an exact null this can be written

`P(σ ≤ m) = (m / A) × F1 × F2 × F3`, with `F1 = E[σ ∧ m] / m`, `F2 = 1 - E[M_m ; σ > m] / E[σ ∧ m]`, `F3 = A / E[M_σ | σ ≤ m]`.

All three factors are at most 1. The bound `m / A` is attained only if the detector never stops early (F1), leaves no expected wealth on paths that have not alarmed (F2) and does not overshoot A (F3). For the two-sided rule the same holds for `N_t = (M^up_t + M^down_t) / 2` with each side at 2A.

### 8.2 Where the slack is (Gaussian null, true baseline, c = 3, one day m = 390)

| Detector | false alarm rate | rate / bound | F1 | F2 (wealth that crossed) | F3 (overshoot) |
|---|---|---|---|---|---|
| two-sided default grids (E9 setting), 10,000 paths | 0.88% | 0.18 | 1.00 | 0.24 | 0.80 |
| two-sided default grids, 40,000 paths (part D2) | 1.03% | 0.21 | | | |
| increase side only, default grid, at A | 1.58% | 0.32 | 0.99 | 0.45 | 0.78 |
| decrease side only, default grid, at A | 0.61% | 0.12 | 1.00 | 0.14 | 0.96 |
| increase side, single bet 0.5 | 3.52% | 0.70 | 0.98 | 0.92 | 0.67 |
| increase side, single bet 0.2 | 2.88% | 0.58 | 0.99 | 0.68 | 0.81 |
| increase side, single bet 0.05 | 0.03% | 0.006 | 1.00 | 0.004 | 0.96 |

(The identity holds up to Monte Carlo noise: its two sides differ by 0% to 16%, because `M_m` is heavy-tailed; the factors are therefore approximate.) The two-sided split does not waste budget: with each side at 2A, the increase side alarms on 0.79% and the decrease side on 0.24% of 40,000 null days, and both sides alarm on the same day in 0.00%, so the union bound is tight.

**Answer to E9's factor of about 4 (here 4.9 with the true baseline).** Stopping costs nothing at one day (F1 ≈ 1). The overshoot costs a factor 0.80: an alarm typically arrives with the statistic at 1.25 A. The main factor is F2 ≈ 0.24: three quarters of the expected wealth sits on paths that have not alarmed. On those paths the decrease side holds 60% of it and the three smallest increase bets (0.02, 0.05, 0.1) another 31%. These bets cannot build wealth quickly (the decrease bets are capped by `λ_max`, the small increase bets by their size), so they keep their expected wealth `w_k m` without ever reaching A within a day. A single bet of 0.5 has almost no uncrossed wealth (F2 = 0.92) and uses 70% of the bound; its only loss is overshoot.

### 8.3 Is the bound sharp?

* **At m = A, exactly, for each side alone or for a one-sided detector at level A.** [PROVED] The constant law `Y ≡ 1` is in H; under it every bet pays `f = 1`, `M_t = t` and the alarm comes at step A exactly. So `P(T_A ≤ A) = 1 = A / A` and the average run length is exactly A: neither guarantee can be improved for any grid. The two-sided rule alarms when one side reaches 2A, so under this law it alarms at step 2A, and this argument does not show that the two-sided bound is sharp.
* **At one day, for a single bet, almost.** [NUMERICAL, rigorous lower bound up to discretisation] Dynamic programming over the whole nonparametric null (at each step the adversary picks any law with `E[U] ≤ 0`; two-point laws suffice by the concave envelope argument of Lemma 2a) gives the worst-case false alarm probability within m steps as a fraction of `m / A`:

| Bet | m = 98 | m = 390 | m = 1950 | m = 3900 |
|---|---|---|---|---|
| increase 0.5 | 0.875 | 0.946 | 0.879 | 0.785 |
| increase 0.2 | 0.517 | 0.853 | 0.880 | 0.794 |
| increase 0.05 | 0.0002 | 0.181 | 0.694 | 0.716 |
| decrease 0.95 λ_max | 0.464 | 0.847 | 0.866 | 0.781 |

The values are achievable by an explicit adversary (they are lower bounds for the worst case, computed on a 700-point log grid of wealth states); an upper bound computed by rounding the other way was useless (rounding errors accumulate), so the only upper bound is `m / A`.

* **For the default mixture, at least 42%.** The best i.i.d. two-point null found for the increase default grid (values 0 and 5 with probabilities 0.8 and 0.2) gives 2.12% at one day (against 1.58% for Gaussian data and the bound 5%); for the decrease grid 1.98% (values 0 and B). A single bet of 0.5 reaches 3.79% under a two-point law.

### 8.4 Is there a sharper threshold that is still provable?

* **For single bets, hardly.** The adversary reaches 85% to 95% of `m / A` at one day for bets 0.2 and 0.5 and for the decrease bet at `0.95 λ_max`, so the bound leaves only 5% to 15% of slack against the worst null law and a sharper threshold that is valid for every null law can gain at most that much. [NUMERICAL]
* **For the default mixture, unknown.** The adversary we found reaches 42% of the bound; a certified upper bound below `m / A` would need an exact computation over an 11-dimensional state (one wealth per bet), which we did not attempt. The small-bet result (18% at one day for the bet 0.05) suggests that grids with small bets could carry a provably sharper short-horizon bound, but we have no proof.
* **Provable improvements that do exist.** (i) The budget can be split unevenly between the sides: thresholds `A_up` and `A_down` with `1/A_up + 1/A_down = 1/A` keep `P(false alarm within m) ≤ m / A` by the union bound. Giving the increase side more than half of the budget, for example `A_up = 4A/3`, saves `log(2A / (4A/3)) = log 1.5 = 0.41` nats on that side, that is `0.41 / g` steps. (ii) A time-varying threshold `A - (m - t)` for `t ≤ m` is valid with the same bound (Ville's inequality applied to the nonnegative supermartingale `M_t + m - t`), but at one day it lowers the threshold by at most 5% (0.05 nats): negligible. (iii) Bets that cannot use the budget can be removed: the designed grids of section 7.2 use 48% of the bound at one day (2.40% with 40,000 paths), against 21% for the default.

## 9. Asymptotic optimality as the false alarm rate goes to zero (tracker B73)

[PROVED theorem; NUMERICAL check] `results/e11e_asymptotic_optimality.csv`, `e11e2_second_order_check.csv`, `e11e2_e10_crosscheck.csv`; figure `results/figures/e11_asymptotic_optimality.png`.

### 9.1 The question and the lower bound

Lai (1998) showed that for a change from a simple pre-change law to a simple post-change law, any detector whose average run length to a false alarm is at least A has worst-case delay at least `(log A / I)(1 + o(1))`, where I is the Kullback-Leibler information of the change (for i.i.d. data this is Lorden's 1971 bound). For our composite null H, Theorem 2b applies Lorden's bound at the least favourable null law P* and gives the lower bound `(log A / KL_inf(Q))(1 + o(1))` for every detector whose run length is at least A under every i.i.d. null law in H. Ram and Ramdas (2026) prove this lower bound in general for composite pre- and post-change classes and show that it is attained for bounded means. Shin, Ramdas and Rinaldo (2024) give non-asymptotic delay upper bounds for mixtures of e-detectors that match the oracle rate `log A / KL` when the baseline e-detectors are built from the right one-parameter family; we could not confirm from the sources we could open (abstract and introduction) that they also state the `KL_inf` lower bound, so we cite Ram and Ramdas for it. The question here is whether our detector, a mixture of linear bets on the truncated statistic, attains `log A / KL_inf`, and at what cost for a finite grid.

**Lemma 5 (the worst case is the change at the first step).** [PROVED] For the SR mixture and for the CUSUM-type mixture, `WADD(T_A) = E_Q[T⁰_A]`, where `T⁰_A` is the alarm time of the detector started from `M_0 = 0` on i.i.d. data from Q.

*Proof.* Each component recursion `x ↦ f_λ(y)(x + 1)` (or `f_λ(y) max(x, 1)`) is nondecreasing in x. By induction, the statistic after ν, started from any state `M^k_{ν-1} ≥ 0`, is at least the statistic started from 0 and driven by the same post-change data, so on `{T_A ≥ ν}` the delay `T_A - ν + 1` is at most `T⁰_A` computed from `Y_ν, Y_{ν+1}, ...`. These are independent of `F_{ν-1}`, so `E_ν[(T_A - ν + 1)⁺ | F_{ν-1}] ≤ E_Q[T⁰_A]`, with equality for ν = 1. ∎

So the simulations with the change at the first step measure the worst-case delay exactly, and any pre-change history (the SR head start of section 4.3) can only help.

### 9.2 Theorem 3

Fix Q with `E_Q[U] > 0` and an interior optimum `0 < λ* < λ_max`, and write `g* = g(λ*) = KL_inf(Q)` (Lemma 2b). For `0 < b < λ_max` let `ū = max(|u_lo|, u_hi)`, `C_b = ū² / (1 - b |u_lo|)²` (a bound on `|g''|` and on `|∂²_λ log(1 + λU)|` over `[0, b]`, for every Q on `[0, B]`) and `o_b = log(1 + b u_hi)` (a bound on every log factor `log f_λ`, `λ ≤ b`, hence on every overshoot).

**Theorem 3.** [PROVED]

**(a) Fixed grid.** For a fixed grid with `g_grid = max_k g(λ_k) > 0`,

`lim_{A→∞} WADD(T_A) / (log A / KL_inf) = KL_inf / g_grid ≥ 1`.

So a fixed grid is first-order optimal at Q exactly when λ* is one of its bets. When it is, `WADD(T_A) ≤ (log A + log(1/w_k) + o_k) / KL_inf` for that bet k (Theorem 1): the penalty of the finite grid is then of second order, at most `(log(1/w_k) + o_k) / KL_inf` steps, bounded in A.

**(b) Grid refined at the right rate.** Fix `0 < a < b < λ_max`. For `K ≥ 2` let `λ_k = a + (k - 1)h`, `h = (b - a)/(K - 1)`, with weights `1/K`. For every Q with `λ*(Q) ∈ [a, b]` and every A,

`WADD(T_A) ≤ (log A + log K + o_b) / (g* - C_b h² / 8)`.

With `K_A = 1 + ⌈(b - a) √(log A)⌉` bets, and uniformly over the laws Q with `λ* ∈ [a, b]` and `g* ≥ g₀ > 0`,

`WADD(T_A) ≤ (log A + ½ log log A + O(1)) / KL_inf`.

More generally every `K_A → ∞` with `log K_A = o(log A)` is first-order optimal; `K_A ≍ √(log A)` minimises the second-order term of this bound.

**(c) Continuous mixture.** Let π be a probability measure on `[0, b]` with a density at least `p₀ > 0` on `[a, b]`, and `M^π_t = ∫ M^λ_t π(dλ)` (a valid SR e-detector by Fubini's theorem, so all false alarm guarantees of `docs/math.md` hold unchanged). If `λ* ∈ (a, b)`, then

`WADD(T_A) ≤ (log A + log log A + O(1)) / KL_inf`.

**What is new here and what is not.** Theorem 3b and 3c are upper bounds on the delay; the lower bound proved in this note (Theorem 2b) is first order only, `(log A / KL_inf)(1 + o(1))`, so the second-order terms `½ log log A` and `log log A` are not shown to be exact. First-order optimality for bounded means is known (Ram and Ramdas 2026). This section gives the refinement rate of the grid and the constants for our detector; we make no claim of priority for the first-order result. We did not check whether Ram and Ramdas also give the `√(log A)` refinement rate.

*Proof of (a).* The upper half is Theorem 1 with the best bet: `WADD ≤ (log A + log(1/w_k) + o_k) / g_grid`. For the lower half take ν = 1 (Lemma 5). For bet k let `S^k_n = Σ_{i ≤ n} log f_{λ_k}(Y_i)`, `S^k_0 = 0`, and `R^k_N = max_{0 ≤ i < n ≤ N} (S^k_n - S^k_i)`. Write `S^k_n = g_k n + D^k_n` with `D^k` a zero-mean random walk. Then `S^k_n - S^k_i ≤ max(g_k, 0)(n - i) + 2 max_{j ≤ N} |D^k_j|`, and `max_{j ≤ N} |D^k_j| / N → 0` almost surely by the strong law of large numbers, so `R^k_N ≤ max(g_k, 0) N + o(N)` almost surely. Since `M_n = Σ_k w_k Σ_{j ≤ n} L^{(j,k)}_n ≤ n max_{k, j} L^{(j,k)}_n`, for all `n ≤ N`

`log M_n ≤ log N + max_k R^k_N ≤ g_grid N + o(N)` almost surely.

Fix ε > 0 and put `N_A = (1 - ε) log A / g_grid`. Then `max_{n ≤ N_A} log M_n ≤ (1 - ε) log A + o(log A) < log A` for all large A, almost surely, so `P(T_A ≤ N_A) → 0` and `E[T_A] ≥ N_A P(T_A > N_A) = N_A (1 - o(1))`. As ε is arbitrary, `E[T_A] ≥ (log A / g_grid)(1 - o(1))`. Divide both halves by `log A / KL_inf`. ∎

*Proof of (b).* Some bet satisfies `|λ_k - λ*| ≤ h/2`. Since `g'(λ*) = 0` and `|g''| ≤ C_b` on `[0, b]`, Taylor's theorem gives `g(λ_k) ≥ g* - C_b (h/2)² / 2 = g* - C_b h² / 8`. Theorem 1 with this bet, `w_k = 1/K` and `o_k ≤ o_b` gives the first display. With `K = K_A` we have `h ≤ 1/√(log A)` and, for `log A ≥ 1`, `K_A ≤ (2 + b - a) √(log A)`, so `log K_A ≤ ½ log log A + log(2 + b - a)`. The denominator is at least `g* (1 - η_A)` with `η_A = C_b / (8 g₀ log A)`, and `(log A + ½ log log A + c) / (g* (1 - η_A)) = (log A + ½ log log A + c) / g* + C_b / (8 g* g₀) + o(1)`. The general statement follows in the same way: the bound is `(log A)(1 + O(h_A²)) / g* + O(log K_A)`. Balancing `log K` against `(log A) C_b (b - a)² / (8 g* K²)` gives `K ≍ √(log A)`. ∎

*Proof of (c).* Take ν = 1 (Lemma 5) and ε > 0 with `[λ* - ε, λ* + ε] ⊂ [a, b]`. Write `h_i(λ) = log(1 + λ U_i)`. Each `h_i` is concave in λ, so on the interval it is at least its value at one of the end points: `h_i(λ) ≥ Z_i := min(h_i(λ* - ε), h_i(λ* + ε))`. Keeping only the start time 1 and the interval,

`M^π_t ≥ ∫_{λ*-ε}^{λ*+ε} Π_{i ≤ t} (1 + λ U_i) π(dλ) ≥ 2 ε p₀ exp(Σ_{i ≤ t} Z_i)`.

The `Z_i` are i.i.d. and bounded. By Taylor's theorem `h_i(λ* ± ε) ≥ h_i(λ*) ± ε h_i'(λ*) - C_b ε² / 2`, so `Z_i ≥ h_i(λ*) - ε |h_i'(λ*)| - C_b ε² / 2` and `E[Z] ≥ g* - ε D - C_b ε² / 2` with `D = E_Q |U / (1 + λ* U)| ≤ ū / (1 - b |u_lo|)`. The Wald argument of Theorem 1, applied to the random walk `Σ Z_i` and the level `log(A / (2 ε p₀))`, gives

`E[T_A] ≤ (log A + log(1/ε) + log(1/(2 p₀)) + o_b) / (g* - ε D - C_b ε² / 2)`.

With `ε = 1 / log A` the numerator is `log A + log log A + O(1)` and the denominator is `g* (1 - O(1/log A))`, which gives the claim. ∎

**The continuous mixture should do better than (c).** [HEURISTIC] The proof gives away a factor in the window: it controls the worst λ in a window of width `1/log A` path by path. A Laplace approximation, `log ∫ e^{S_n(λ)} π(dλ) ≈ max_λ S_n(λ) - ½ log(n σ*²) + log p(λ*) + O(1)` with `σ*² = -g''(λ*)`, suggests `(log A + ½ log log A + O(1)) / KL_inf`, the same order as the refined grid of (b). The same `½ log log` term appears in the expected sample size of mixture-type tests of power one in exponential families (Robbins and Siegmund 1974; Pollak and Siegmund 1975; we checked their bibliographic records, not the exact form of their expansions).

**Corollary 6 (the cost of a gap in the grid).** [PROVED expansion; the small-change form is a Taylor approximation] If λ* is not a bet, the grid loses a first-order fraction of the speed, `1 - g_grid / KL_inf`, and `KL_inf - g(λ_k) = ½ σ*² (λ_k - λ*)² + O(|λ_k - λ*|³)` with `σ*² = E_Q[U² / (1 + λ* U)²]`, the variance under Q of `U / (1 + λ* U)`. For a small change (`E_Q U = d`, `Var_Q U = v`, d small) `λ* ≈ d/v`, `KL_inf ≈ d² / (2v)` and `g(λ) / KL_inf ≈ 1 - (λ / λ* - 1)²`: a bet that is wrong by a factor r keeps the fraction `1 - (r - 1)²` of the growth. A geometric grid whose neighbouring bets differ by the factor ρ therefore keeps at least `1 - ((ρ - 1)/(ρ + 1))²` of the growth for every λ* between its smallest and largest bets (the worst case is λ* at the arithmetic midpoint of two neighbours): 0.89 for ρ = 2 and 0.97 for ρ = 1.43. In delay, a gap costs `log A (1/g_grid - 1/KL_inf)` steps, which grows without bound in A, while an extra bet costs about `log(K + 1) - log K` nats, which does not. This is the trade-off behind the design of section 7 and the refinement rate of Theorem 3b.

### 9.3 Numerical check

[NUMERICAL] Mean delay divided by `log A / KL_inf`, change at the first step (equal to the worst-case delay by Lemma 5), c = 3, 1000 paths per case and the same paths for all detectors, A = 10^2 to 10^8. "Fine mixture" is 100 bets equally spaced on `(0, cap)` with equal weights (cap = 0.95 for increases, `0.99 λ_max` for decreases), a discretised continuous mixture; "robust" is the grid of section 7.3. In brackets: the fixed-grid limit `KL_inf / g_grid` of Theorem 3a. No path was left without an alarm, except for the Gaussian-tuned grid under t3 data (section 7.2).

| Case | oracle bet, A = 10^2 / 10^4 / 10^8 | second-order formula, 10^4 / 10^8 | fine mixture, 10^4 / 10^8 | default grid, 10^4 / 10^8 (limit) | robust grid, 10^8 (limit) |
|---|---|---|---|---|---|
| Gaussian, up 1.25 | 0.36 / 0.66 / 0.83 | 0.63 / 0.82 | 0.71 / 0.90 | 0.73 / 0.90 (1.02) | 0.88 (1.10) |
| Gaussian, up 1.5 | 0.58 / 0.79 / 0.89 | 0.77 / 0.89 | 0.83 / 0.93 | 0.90 / 0.95 (1.02) | 0.93 (1.01) |
| Gaussian, up 2 | 0.80 / 0.90 / 0.94 | 0.89 / 0.95 | 0.96 / 0.97 | 1.11 / 1.08 (1.04) | 1.12 (1.14) |
| Gaussian, down 0.8 | 0.31 / 0.62 / 0.82 | 0.62 / 0.81 | 0.80 / 0.94 | 0.77 / 0.91 (1.02) | 0.93 (1.10) |
| Gaussian, down 0.5 | 0.50 / 0.74 / 0.87 | 0.74 / 0.87 | 0.96 / 1.03 | 0.92 / 1.00 (1.05) | 1.05 (1.17) |
| Student t3, up 1.25 | 0.24 / 0.57 / 0.78 | 0.54 / 0.77 | 0.65 / 0.89 | 0.64 / 0.87 (1.08) | 0.81 (1.00) |
| Student t3, up 1.5 | 0.43 / 0.67 / 0.84 | 0.69 / 0.85 | 0.72 / 0.90 | 0.74 / 0.90 (1.02) | 0.88 (1.08) |
| Student t3, up 2 | 0.68 / 0.83 / 0.90 | 0.82 / 0.91 | 0.87 / 0.94 | 0.97 / 0.98 (1.00) | 0.99 (1.05) |
| Student t3, down 0.8 | 0.19 / 0.53 / 0.76 | 0.51 / 0.76 | 0.57 / 0.80 | 0.58 / 0.82 (1.01) | 0.78 (1.01) |
| Student t3, down 0.5 | 0.48 / 0.72 / 0.86 | 0.72 / 0.86 | 0.87 / 0.96 | 0.84 / 0.94 (1.00) | 0.93 (1.03) |
| Jumps, up 1.25 | 0.35 / 0.61 / 0.80 | 0.63 / 0.81 | 0.67 / 0.88 | 0.68 / 0.87 (1.02) | 0.86 (1.09) |
| Jumps, up 1.5 | 0.57 / 0.77 / 0.88 | 0.77 / 0.88 | 0.81 / 0.93 | 0.88 / 0.94 (1.02) | 0.92 (1.00) |
| Jumps, up 2 | 0.79 / 0.90 / 0.95 | 0.88 / 0.94 | 0.95 / 0.99 | 1.11 / 1.09 (1.04) | 1.12 (1.13) |
| Jumps, down 0.8 | 0.27 / 0.58 / 0.79 | 0.58 / 0.79 | 0.69 / 0.86 | 0.67 / 0.85 (1.05) | 0.83 (1.00) |
| Jumps, down 0.5 | 0.45 / 0.70 / 0.85 | 0.70 / 0.85 | 0.87 / 0.96 | 0.84 / 0.93 (1.00) | 0.95 (1.07) |

**Reading.**

* **The asymptotic regime is far away.** Even the oracle bet, which attains the lower bound to first order (Theorem 3a with `g_grid = KL_inf`), is 5% to 24% *faster* than `log A / KL_inf` at A = 10^8, and 10% to 47% faster at A = 10^4. This is not a contradiction: the lower bound is `(log A / KL_inf)(1 + o(1))`, and the `o(1)` is the SR memory term of section 4.3. The second-order formula `1 + (κ - ζ) / log A` (Monte Carlo κ and ζ) matches the simulated oracle ratio to within 0.024 in all 15 cases for every `A ≥ 10^4` (within 0.037 from `A = 10^3`), so the ratio approaches 1 like `1 - (ζ - κ)/log A`, with `ζ - κ` between 1.0 and 4.5 nats.
* **Fixed grids behave as Theorem 3a says.** The default grid's ratio to the oracle falls toward its limit `KL_inf / g_grid` (1.00 to 1.08): for example Gaussian up 2, where λ* = 0.64 is not in the grid, 1.23 at A = 10^4 and 1.15 at 10^8 against the limit 1.04. Its excess over the oracle, measured in nats as `(delay - oracle delay) × KL_inf`, grows between A = 10^4 and 10^8 by 0.12 to 0.95 nats.
* **The fine mixture is close to the oracle and slowly getting closer.** Its delay is 1.04 to 1.13 times the oracle's for increases and 1.07 to 1.19 for decreases at A = 10^8 (1.05 to 1.30 at 10^4). Its excess in nats grows between A = 10^4 and 10^8 by 0.18 to 1.13 (`e11e2_second_order_check.csv`), against 0.35 for `½ log log A` and 0.69 for `log log A`. So the growth is of the predicted order, but at these alarm levels it is not yet the asymptotic `½ log log A` (the excess grows fastest for small increases, where the posterior over λ concentrates slowly). For the Gaussian decreases λ* lies at or beyond `0.99 λ_max`, outside the mixture's range, so the interior condition of Theorem 3c fails and a small first-order loss remains (limit 1.001 and 1.014).
* **Numerically, first-order optimality cannot be confirmed or refuted at A ≤ 10^8.** What the numbers do confirm is the structure of the theory: the second-order formula for the oracle, the fixed-grid limits, and the slow decline of the mixture's relative overhead.

**What fraction of the likelihood ratio detector's speed is kept.** [NUMERICAL] As A → ∞ the delay ratio between the correctly specified likelihood ratio detector (for Gaussian data, the Gaussian LR detector of E4; for the other models, the LR detector for a scale change of the true model) and the robust oracle bet is `KL_inf / KL`, where KL is the information of the raw returns:

| Model | up 1.25 | up 1.5 | up 2 | down 0.8 | down 0.5 |
|---|---|---|---|---|---|
| Gaussian | 0.88 | 0.75 | 0.53 | 0.80 | 0.28 |
| Jumps | 0.89 | 0.81 | 0.69 | 0.65 | 0.24 |
| Student t3 | 0.94 | 0.98 | 0.96 | 0.72 | 0.44 |

For a grid, multiply by `g_grid / KL_inf` (the inverse of the limits in the table above). Heavy tails carry little information about scale (the t3 KL at sd ratio 1.25 is 0.026, against 0.058 for Gaussian data), so the robust detector loses little against the correctly specified parametric detector under t3 data; the Gaussian LR detector itself is not valid there (E4). Section 6.1 splits the Gaussian ratio into truncation and nonparametric losses.

### 9.4 Cross-reference with E10

E10 (`docs/report_parts/E10.md`, Q10d and Q10e) computed the per-step KL divergence of a volatility break seen through a frozen GARCH forecast and used `ln A / I` (Lai 1998) as a guide. Its observed median delays were about 0.6 times the step at which the cumulative KL reaches `ln A`, which E10 attributed mainly to the plug-in baseline. Section 4.3 gives a second explanation that needs no baseline error. In an i.i.d. Gaussian version of E10's steady state (variance ratio s1/s0 at t = 3000 from `results/e10e_transient.csv`, true baseline, increase side, c = 3, A = 28,800, 1000 paths; `e11e2_e10_crosscheck.csv`):

| Break (E10) | variance ratio | `ln A / KL` (steps) | memory term ζ (Brownian) | second-order prediction, oracle bet | simulated oracle bet: mean / median | simulated default grid: mean / median | mean / (`ln A / KL`), oracle / default |
|---|---|---|---|---|---|---|---|
| k = 2 | 1.117 | 3218 | 6.38 (6.39) | 1275 | 1323 / 1177 | 1493 / 1304 | 0.41 / 0.46 |
| k = 4 | 1.190 | 1285 | 5.51 (5.50) | 634 | 650 / 571 | 741 / 687 | 0.51 / 0.58 |

At these small growth rates (`KL_inf` = 0.0031 and 0.0077 per step, 98% and 96% of KL) the memory term is 54% to 62% of `ln A = 10.27`, so even with the true baseline the SR detector alarms at about 0.4 to 0.6 times `ln A / KL`. The memory term alone is therefore large enough to produce E10's factor of about 0.6; the plug-in baseline is not needed to explain it, although it may contribute. This is a hypothesis for E10, not a result: we did not rerun E10's GARCH setting, whose transient and time-varying KL differ from this i.i.d. version.

## 10. What is proved, what is heuristic, what is numerical

| Claim | Status | Where |
|---|---|---|
| Growth is concave; unique oracle bet | PROVED | Lemma 0 |
| Growth ≤ KL; information ladder; speed limits for decreases | PROVED | Lemma 1, Corollaries 3 and 4 |
| `max g = KL_inf`, every nonparametric e-variable is dominated by a bet | PROVED | Lemma 2 |
| Non-asymptotic delay bound for every A, ν and history | PROVED | Theorem 1 |
| Worst-case delay of an SR mixture = delay with the change at step 1 | PROVED | Lemma 5 |
| Fixed grid: delay ratio to `log A / KL_inf` tends to `KL_inf / g_grid` | PROVED | Theorem 3a |
| Refined grid with about `√(log A)` bets: delay at most `(log A + ½ log log A + O(1)) / KL_inf` (an upper bound; first-order optimality for bounded means is known from Ram and Ramdas 2026) | PROVED | Theorem 3b |
| Continuous mixture: delay at most `(log A + log log A + O(1)) / KL_inf` | PROVED | Theorem 3c |
| Lower bound `log A / KL_inf` for any detector valid over H | PROVED given Lorden (1971); STANDARD (Ram and Ramdas 2026) | Theorem 2b |
| Continuous mixture: `½ log log A` | HEURISTIC (Laplace) | section 9.2 |
| Second-order delay `(log A + κ - ζ) / g`, Brownian ζ | HEURISTIC, checked to 3% (A = 7800) and 0.024 in ratio (A ≥ 10^4) | sections 4.3, 5, 9.3 |
| False alarm identity with three factors (exact null) | PROVED | section 8.1 |
| Bound attained at m = A; average run length bound sharp (one side alone, or a one-sided detector at level A) | PROVED | section 8.3 |
| Worst null reaches 85% to 95% of `m / A` at one day for single bets | NUMERICAL (dynamic programming lower bound) | section 8.3 |
| Price-of-robustness values, c recommendations, bet designs, E10 cross-check | NUMERICAL | sections 6, 7, 9 |

## 11. References

(Bibliographic records checked by web search. We could not open the full texts of the arXiv papers from this environment, so statements about their content rest on their abstracts and introductions.)

* Dufresne, D. (1990). The distribution of a perpetuity, with applications to risk theory and pension funding. *Scandinavian Actuarial Journal* 1990(1), 39-79.
* Grünwald, P., de Heide, R. and Koolen, W. (2024). Safe testing. *Journal of the Royal Statistical Society Series B* 86(5).
* Honda, J. and Takemura, A. (2010). An asymptotically optimal bandit algorithm for bounded support models. *Proceedings of COLT 2010*.
* Lai, T. L. (1998). Information bounds and quick detection of parameter changes in stochastic systems. *IEEE Transactions on Information Theory* 44(7), 2917-2929.
* Lai, T. L. and Siegmund, D. (1977). A nonlinear renewal theory with applications to sequential analysis I. *Annals of Statistics* 5(5), 946-954.
* Larsson, M., Ramdas, A. and Ruf, J. (2025). The numeraire e-variable and reverse information projection. arXiv:2402.18810 (version 4, February 2025; reported elsewhere as published in the *Annals of Statistics*, journal version not checked).
* Lorden, G. (1970). On excess over the boundary. *Annals of Mathematical Statistics* 41(2), 520-527.
* Lorden, G. (1971). Procedures for reacting to a change in distribution. *Annals of Mathematical Statistics* 42(6), 1897-1908.
* Pollak, M. (1985). Optimal detection of a change in distribution. *Annals of Statistics* 13(1), 206-227.
* Pollak, M. (1987). Average run lengths of an optimal method of detecting a change in distribution. *Annals of Statistics* 15(2), 749-779.
* Pollak, M. and Siegmund, D. (1975). Approximations to the expected sample size of certain sequential tests. *Annals of Statistics* 3(6), 1267-1282.
* Ram, A. and Ramdas, A. (2026). Asymptotically optimal sequential change detection for bounded means. arXiv:2602.05272.
* Robbins, H. and Siegmund, D. (1974). The expected sample size of some tests of power one. *Annals of Statistics* 2(3).
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science* 2(2), 229-260.
* Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *Journal of the Royal Statistical Society Series B* 86(1).
* Woodroofe, M. (1982). *Nonlinear Renewal Theory in Sequential Analysis*. CBMS-NSF Regional Conference Series in Applied Mathematics 39, SIAM.
