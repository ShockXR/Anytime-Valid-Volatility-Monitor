# One guarantee for detectors at 1, 5 and 15 minutes

This note derives the false alarm guarantee of the multi-frequency monitor (`src/avmonitor/multifreq.py`, experiment E17). It builds on `docs/math.md` section 5 (Lemma 1, Theorems A and B) and uses the same notation.

**In plain words.** A volatility change can be easiest to see at one minute, or only at fifteen minutes (for example when one-minute prices are noisy). Running e-detectors at several sampling frequencies and averaging their wealth (Vovk and Wang 2021: an average of e-values is an e-value) is the obvious way to avoid choosing. The catch is that the detectors look at the data on different clocks. A fifteen-minute detector updates once per block, and one minute before the block ends most of its next bet is already decided by the returns seen so far. On the one-minute clock it is therefore not a fair bet, and the averaging argument breaks if the average is checked every minute (section 3). We fix this by crediting each coarse detector, between its block ends, only with the smallest value its wealth can take at the next block end. The resulting monitor is checked every minute and keeps a proved bound that is better than the one-minute monitor's `m / A` (Theorem 1). We also prove weaker bounds for the plain average and for a weighted maximum (Theorem 3), say exactly which null each frequency needs and when the one-minute null implies it (section 2.1), handle the two sides and recalibration (Corollary 4), and show that averaging cannot remove the price of not knowing the right frequency: in the worst case it is a factor J in the threshold for J frequencies, the same as a Bonferroni split (section 5).

---

## 1. Set-up

* One-minute returns `r_t`. `F_t` is the information up to minute t (the calibration data and the monitored returns up to t). Monitoring starts after a calibration window that ends at minute T; everything fixed at T (truncation levels, baselines, bets) is `F_T`-measurable.
* A finite set of frequencies `K` (block lengths in minutes, here {1, 5, 15}), `J = |K|`, and weights `w_k ≥ 0` with `Σ_k w_k = 1`.
* **Blocks.** For frequency k, block b (b = 1, 2, ...) covers minutes `T + (b - 1) k + 1, ..., T + b k`. Its *start* is `s = T + (b - 1) k` (the last minute before it) and its *end* is `e = s + k`. The block return `R_b^(k)` is the sum of the one-minute returns of the block (**plain** blocks), or of its last `k - 1` returns (**gap** blocks, k ≥ 2: consecutive blocks then share no price, as in the non-overlapping block monitor of `alternatives/preaveraging/README.md`, section 2). `L_k` is the number of returns summed (k or k - 1).
* **Statistic.** `X_b^(k) = min((R_b^(k))², τ_k²)` with a truncation level `τ_k` and baselines `m↑_k`, `m↓_k`, all fixed at T. In the code `τ_k = c` times the scaled MAD of k-minute returns in the calibration window and `m↑_k = (1 + b) m_k`, `m↓_k = (1 - b) m_k` with `m_k` the mean of the clipped squares (plug-in); the theorems assume the true values (section 8).
* **Bets.** For each frequency a grid of bets, fixed at T: increase factors `f = 1 + λ_j (X / m↑_k - 1)` with `λ_j ∈ [0, 1)`; decrease factors `f = 1 + λ_j (1 - X / m↓_k)` with `λ_j = φ_j / (B_k - 1)`, `φ_j ∈ [0, 1)` and `B_k = τ_k² / m↓_k > 1`. Every factor is at least `ℓ_j > 0`, where `ℓ_j = 1 - λ_j` (increase, because `X ≥ 0`) or `ℓ_j = 1 - φ_j` (decrease, because `X ≤ τ_k²`).
* **Statistics.** For each bet the Shiryaev-Roberts recursion over blocks, `M_b = f_b (M_{b-1} + 1)`, `M_0 = 0`; the frequency's statistic `M^(k)` is the mean over its grid. A block that contains a missing one-minute return is skipped (`M_b = M_{b-1}`), as the one-minute monitor skips a missing return (`docs/math.md` 5.6, item 5). In minute time, `M_t^(k)` is the value after the last completed block (the **frozen** value).

## 2. The null each frequency needs

For frequency k the increase and decrease nulls are

* `H↑_k`: `E[X_b^(k) | F_s] ≤ m↑_k` for every block b, where s is the block's start,
* `H↓_k`: `E[X_b^(k) | F_s] ≥ m↓_k` for every block b.

The conditioning is on **all one-minute information at the start of the block**, not only on earlier k-minute blocks. This is what makes the frequencies comparable on one clock, and it is what the code needs (nothing in the monitor uses information that the null does not condition on).

**Lemma 1′ (one frequency, block clock).** Under `H↑_k`, with the missingness condition of `docs/math.md` 6.7 (given `F_s`, whether the block is complete says nothing about the size of its return), each increase statistic satisfies `E[M_b | F_s] ≤ M_{b-1} + 1` for every block b with start s. The same holds for the decrease side under `H↓_k`.

*Proof.* `M_{b-1}` is `F_s`-measurable. If the block is complete, `M_b = f_b (M_{b-1} + 1)`, otherwise `M_b = M_{b-1}`. So `E[M_b | F_s] = (M_{b-1} + 1) E[f_b 1{complete} | F_s] + M_{b-1} P(incomplete | F_s)`. Under the null and the missingness condition, `E[f_b 1{complete} | F_s] ≤ P(complete | F_s)`, which gives `E[M_b | F_s] ≤ M_{b-1} + P(complete | F_s) ≤ M_{b-1} + 1`. ∎

On its own clock (the information at its block ends), a k-minute detector is therefore an e-detector, and Theorems A and B of `docs/math.md` hold for it with m counted in blocks: P(false alarm within ⌊m / k⌋ blocks) ≤ ⌊m / k⌋ / A. This is the block monitor of the pre-averaging alternative.

### 2.1 When does the one-minute null imply the k-minute null?

Not in general. Two things separate them.

1. **Truncation is not additive.** `min(R², τ_k²)` is not a sum of one-minute pieces, so a statement about the conditional means of `min(r_t², τ_1²)` says nothing exact about `E[min(R², τ_k²) | F_s]`.
2. **The variance of a sum depends on the dependence inside the block.** The one-minute null restricts each return's conditional size, not how returns inside a block line up. Example: let `e_t` be i.i.d. N(0, 1) and `r_t = S_b |e_t|`, where `S_b = ±1` is a fair random sign drawn at the start of each block and kept for the whole block. Then `r_t² = e_t²`, so the one-minute null holds with equality relative to an i.i.d. N(0, 1) calibration period. But `R_b = S_b Σ_t |e_t|` has `E[R_b²] = k + k (k - 1) (2 / π)`, which is 17.7 instead of 5 for k = 5. The five-minute variance has risen 3.5-fold with no change at one minute; only the coarse detector can see this, and relative to the calibration period its null is false.

It does hold in three useful cases.

* **(a) A constant conditional law.** If the conditional law of `r_t` given `F_{t-1}` is one fixed law `P₀` for every t (as in Corollary 6.5 of `docs/math.md`), then `r_t` is independent of `F_{t-1}` for every t, so the returns are i.i.d. with law `P₀`, the block `(r_{s+1}, ..., r_{s+k})` is independent of `F_s`, and `E[X_b^(k) | F_s] = E_{P₀}[min(R², τ_k²)]` for every b. With that value as the true baseline, `H↑_k` and `H↓_k` hold with equality, for every k at once. This covers jumps, fat tails and price ticks, as long as they are i.i.d.
* **(b) A scale model with the scale fixed inside each block.** If `r_t = σ_b ε_t` for the minutes of block b, with `ε_t` i.i.d., independent of the past, and `σ_b` known at the block's start, then `E[X_b^(k) | F_s] = g_k(σ_b)` with `g_k(σ) = E[min(σ² (Σ_t ε_t)², τ_k²)]`, which is nondecreasing in σ because `min(σ² y², τ_k²)` is nondecreasing in σ for every y. So `σ_b ≤ σ₀` for every block implies `H↑_k` (and `σ_b ≥ σ₀` implies `H↓_k`) with `m_k = g_k(σ₀)`. If the scales are constant over blocks of the coarsest frequency, this holds for all frequencies at once. (If σ changes inside a block, the same conclusion holds for Gaussian ε and scales that are fixed at the block's start, because R is then N(0, Σ σ_t²); we make no general claim beyond that.)
* **(c) Raw squares under martingale differences.** If `E[r_t | F_{t-1}] = 0` and `E[r_t² | F_{t-1}] ≤ μ` for every t, then for any block, `E[R² | F_s] = Σ_t E[r_t² | F_s] + 2 Σ_{t<u} E[r_t E[r_u | F_{u-1}] | F_s] = Σ_t E[E[r_t² | F_{t-1}] | F_s] ≤ L_k μ`. So the *untruncated* increase null at one minute implies the untruncated increase null at k minutes with baseline `L_k μ`. This does not extend to truncated statistics (point 1) or to the decrease side (`docs/math.md` 5.3).

**Price noise.** Observed log prices `Y_j = X_j + u_j` with i.i.d. noise of variance `ω²`, independent of the efficient price X, and `γ = ω / (σ √δ)` the noise-to-signal ratio of one-minute returns (the model of `docs/report_parts/ALT1.md`). The one-minute null fails (ALT1). For plain k-minute blocks, `R_b` and `R_{b-1}` share the price `Y_s`, so `Cov(R_b, R_{b-1}) = -ω²` and, in units of `σ² δ`, `Var(R_b) = k + 2 γ²` and their correlation is `ρ_k = -γ² / (k + 2 γ²)`. For Gaussian noise and conditioning on the previous block alone, `E[R_b² | R_{b-1}] = Var(R) + ρ_k² (R_{b-1}² - Var(R))`, so the conditional mean moves with the previous block, by a factor `ρ_k² = γ⁴ / (k + 2 γ²)²`: at γ = 1 this is 0.111 at one minute (the value of ALT1), 0.020 at five minutes and 0.0035 at fifteen. Plain coarse blocks are therefore much closer to the null but not exactly on it. **Gap blocks are exactly on it**: block b then uses only the prices `Y_{s+1}, ..., Y_{s+k}`, so `R_b = (X_{s+k} - X_{s+1}) + u_{s+k} - u_{s+1}`. If X has independent increments (constant σ, i.i.d. jumps), both parts are independent of `F_s`, which is generated by `Y_0, ..., Y_s`. Then `E[X_b^(k) | F_s]` is a constant and the null holds with equality with the true baseline, which includes the noise term `2 ω²`. (`tests/test_multifreq.py` checks the two correlations by simulation.)

**What the combination needs.** Every frequency with a positive weight needs its own two-sided null. Under i.i.d. returns (case a) all of them hold with the true baselines. Under price noise the one-minute null fails and only gap blocks are exactly valid, so a combination that keeps weight on one minute inherits a violation (E17 measures how much).

## 3. Why the plain average checked every minute has no proof from Lemma 1

The natural claim, "on the one-minute filtration `E[M_t^(k) | F_{t-1}] ≤ M_{t-1}^(k) + 1{t ends a k-block}`", is **false**, and so is the optional stopping property that Theorems A and B rest on.

**Proposition 2.** Take k = 2, i.i.d. N(0, 1) returns, no truncation (τ = ∞), the true baseline `m_2 = 2` and one increase bet λ ∈ (0, 1), so `f_b = 1 + λ (R_b² / 2 - 1)`.

1. At the end of a block (minute `e = s + 2`), `E[f_b | F_{e-1}] = 1 + λ (r_{s+1}² + 1 - 2) / 2 = 1 + λ (r_{s+1}² - 1) / 2`, which exceeds 1 whenever `|r_{s+1}| > 1`. So `E[M_e | F_{e-1}] > M_{e-1} + 1` on that event.
2. Let τ stop in the middle of block 2 (minute 3) if `|r_3| ≤ 1`, and at its end (minute 4) otherwise. τ is a stopping time on the one-minute filtration, and `E[M_τ] = 1 + p + λ E[(Z² - 1) 1{|Z| > 1}] = 1 + p + 0.484 λ` with `p = P(|Z| > 1) = 0.317` and Z standard normal, while the expected number of completed blocks is `1 + p`. So `M_t - (completed blocks)` is not a supermartingale on the one-minute filtration.
3. Counting blocks *started* instead does not rescue it. Stop in block B at its first minute if `|r_first| ≤ 1` and at its end otherwise. Then `E[M_τ] - B = 0.242 λ B - 0.683`, which is positive for `B > 2.82 / λ` (from the 6th block on at λ = 0.5).

*Proof.* 1 is direct, since `E[(r_{s+1} + r_{s+2})² | F_{s+1}] = r_{s+1}² + 1`. For 2, `M_1 = f_1` is independent of `(r_3, r_4)` with mean 1, so `E[M_τ] = E[f_1] P(|r_3| ≤ 1) + E[f_2 1{|r_3| > 1}] E[f_1 + 1]`, and `E[f_2 1{|r_3| > 1}] = p + (λ / 2) E[(Z² - 1) 1{|Z| > 1}]`, with `E[Z² 1{|Z| > 1}] = 2 (φ(1) + 1 - Φ(1)) = 0.8013`, so `E[(Z² - 1) 1{|Z| > 1}] = 0.4840`. For 3, the same computation with `E[M_{B-1}] = B - 1` gives `E[M_τ] = (B - 1)(1 - p) + B (p + 0.242 λ) = B - (1 - p) + 0.242 λ B`. ∎ (`tests/test_multifreq.py` checks part 2 by simulation.)

The coarse detector is an e-detector in a *smaller* filtration (the information at its block ends), and validity in a smaller filtration does not carry over to a larger one. This is the general problem studied by Choe and Ramdas (arXiv:2402.09698), who lift such e-processes to the larger filtration with "adjusters", functions of the running maximum, at a logarithmic cost in power. Here the structure is simpler (the coarse detector is behind by less than one block), and section 4 uses it directly.

Proposition 2 does not show that the plain every-minute average *raises more false alarms* than `Σ_k w_k ⌈m / k⌉ / A`; it shows that Lemma 1 cannot prove that it does not. Theorem 3 gives the bound we can prove for it, and E17 measures its actual rate.

## 4. The certified average and its guarantee

**Certified values.** For a bet j of frequency k and side `↑` or `↓`, let s be the last block end at or before minute t. Define

* `C_t = M_t` if t is a block end (or t = T), and
* `C_t = min(M_s, ℓ_j (M_s + 1))` if t lies strictly inside the next block.

`C` is computable at t. For k = 1 every minute is a block end, so `C = M`. The **certified average** of a side is `C̄_t = Σ_k w_k (grid mean of C_t^(k))`. The monitor alarms at the first minute `t > T` with `max(C̄↑_t, C̄↓_t) ≥ 2A`. At minutes that end a block of every frequency (multiples of 15 minutes for {1, 5, 15}), the certified average equals the plain average.

**Lemma 2 (a one-minute supermartingale above the certified value).** For a bet j of frequency k, let `Q_t = E[M_{e(t)} | F_t]`, where `e(t)` is the end of the block that contains minute t (`e(t) = t` at block ends, so `Q_t = M_t` there, and `Q_T = 0`). Then, under the frequency's null and the missingness condition,

1. `Q_t ≥ C_t` for every t, and
2. `E[Q_t | F_{t-1}] ≤ Q_{t-1} + 1{t is the first minute of a block}`.

*Proof.* 1. At block ends both equal `M_t`. Strictly inside a block with start s, the value at the block's end is either `f (M_s + 1)` with `f ≥ ℓ_j`, or `M_s` if the block is skipped, so `M_{e(t)} ≥ min(ℓ_j (M_s + 1), M_s) = C_t` surely. `C_t` is `F_t`-measurable, so taking `E[· | F_t]` gives `Q_t ≥ C_t`.
2. If t is not the first minute of its block, `e(t) = e(t - 1)` and `E[Q_t | F_{t-1}] = E[E[M_{e(t)} | F_t] | F_{t-1}] = Q_{t-1}` (tower property; this includes the block's end). If t is the first minute of a block with start `s = t - 1`, then `Q_{t-1} = M_s` and `E[Q_t | F_s] = E[M_{s+k} | F_s] ≤ M_s + 1` by Lemma 1′ (for k = 1 this is Lemma 1 of `docs/math.md` itself). ∎

**Theorem 1 (certified average).** Suppose that for every frequency k with `w_k > 0` the nulls `H↑_k` and `H↓_k` hold with the baselines used by the monitor, and the missingness condition holds. Put `ρ = Σ_k w_k / k` and `κ = Σ_k w_k (1 - 1 / k)`. Then for every integer m ≥ 1:

1. `P(an alarm at one of the minutes T + 1, ..., T + m) ≤ Σ_k w_k ⌈m / k⌉ / A`, which is at most `m ρ / A + κ / A` and equals `m ρ / A` when every k divides m;
2. the expected number of monitored minutes before the first alarm is at least `(A - κ) / ρ`.

*Proof.* Let `S_t = Σ_k w_k (grid mean of Q^(k↑)_t + grid mean of Q^(k↓)_t) / 2` and `c(t) = Σ_k w_k ⌈(t - T) / k⌉`, the weighted number of blocks started by minute t. By Lemma 2 and linearity, `E[S_t | F_{t-1}] ≤ S_{t-1} + c(t) - c(t - 1)` (a weighted average of the indicators of Lemma 2). `S ≥ 0` and `S_T = 0`.
1. For `T ≤ t ≤ T + m` the process `S_t + c(T + m) - c(t)` is a nonnegative supermartingale (c is not random) with initial value `c(T + m) = Σ_k w_k ⌈m / k⌉`. By Ville's inequality the probability that it ever reaches A before `T + m` is at most `c(T + m) / A`. An alarm at t means `max(C̄↑_t, C̄↓_t) ≥ 2A`, so `(C̄↑_t + C̄↓_t) / 2 ≥ A`, and by Lemma 2 (part 1, summed with the weights) `S_t ≥ A`, so the supermartingale is at least A too. Finally `⌈m / k⌉ ≤ m / k + (k - 1) / k`.
2. As in Theorem A of `docs/math.md`: if the alarm time τ has infinite mean there is nothing to prove. Otherwise optional stopping at the bounded times `τ ∧ n` gives `E[S_{τ∧n}] ≤ E[c(τ ∧ n)]`, and Fatou's lemma with `S_τ ≥ A` gives `A ≤ E[c(τ)] ≤ ρ E[τ - T] + κ`. ∎

**Remarks.**

* *Better than m / A.* With the one-minute monitor's level A, the bound `Σ_k w_k ⌈m / k⌉ / A` is below `m / A` as soon as a coarse frequency has positive weight. We use the gain the other way: E17 lowers A until the bound equals the one-minute monitor's 5% per day (section 5).
* *What the certification costs.* On the increase side `ℓ_j = 1 - λ_j ≥ 0.5` for the default grid, so between block ends a coarse frequency is credited with at least about half its wealth. On the decrease side the default fractions reach 0.95, so `ℓ_j` can be 0.05: coarse evidence for a decrease is credited mostly at block ends. Nothing is lost at block ends.
* *Special cases.* One frequency alone is the block monitor (bound `⌈m / k⌉ / A`; a lone coarse detector can only alarm at block ends, because its certified value never exceeds its frozen value). Checking the plain average only at common block ends is a monitor that alarms no earlier than the certified one, so it is covered too.
* *One side alone.* Under `H↑_k` for all k, the increase side alone at level A has `P(alarm within m) ≤ Σ_k w_k ⌈m / k⌉ / A` (the same proof with `S = Σ_k w_k Q^(k↑)`).

## 5. Choosing A and the weights; the price of not knowing the frequency

**Level.** For a budget α per horizon m (E17: α = 5% per 1440-minute day), `A = Σ_k w_k ⌈m / k⌉ / α`. For one frequency alone this is `A_k = ⌈m / k⌉ / α`: 28,800, 5,760 and 1,920 for 1, 5 and 15 minutes (`multifreq.level_for_budget`).

**Handicap.** Suppose only frequency j carries evidence. The certified average then alarms (at a block end of j) when `w_j M^(j) ≥ 2A`, that is when `M^(j) ≥ h_j × 2 A_j` with

`h_j = Σ_k w_k ⌈m / k⌉ / (w_j ⌈m / j⌉)`.

**Proposition 3.** For any weights, `max_j h_j ≥ J`, with equality exactly when `w_j ⌈m / j⌉` is the same for every j ("balanced" weights; `w_j ∝ j` when every k divides m).

*Proof.* With `v_j = w_j ⌈m / j⌉ / Σ_k w_k ⌈m / k⌉` (nonnegative, summing to 1), `Σ_j v_j h_j = Σ_j 1 = J`, so the `v`-weighted average of the handicaps is J and the largest is at least J, with equality only if all equal J. ∎

The e-detector's delay grows like the logarithm of its threshold divided by the information it collects per minute (`docs/theory_delay.md`), so a handicap h adds about `log h` to the logarithm of the threshold. For J = 3 and balanced weights this is `log 3 = 1.10` against `log(2 × 28,800) = 10.96` for one minute (about 10%) and `log(2 × 1,920) = 8.25` for fifteen minutes (about 13%). Equal weights (1/3 each) give handicaps `k (1 + 1/5 + 1/15)`: 1.27, 6.3 and 19, so they are almost free when one minute is the right frequency and costly when fifteen minutes is (`log 19 = 2.9`).

The weighted maximum of Theorem 3 has exactly the same single-frequency thresholds (`w_j M^(j) ≥ 2A`) at the same proved budget, so a factor J is also the price of a plain Bonferroni split over the J frequencies. The average can only do better than the maximum when several frequencies carry comparable evidence at the same moment. Because a detector's wealth grows exponentially at a rate that differs between frequencies, the fastest one usually dominates the sum, so we expect the average to track the best single frequency with a delay penalty of about `log h` rather than to beat it. E17 tests this.

## 6. Weaker bounds for the plain average and for the weighted maximum

**Theorem 3.** Under the nulls of Theorem 1, with `N_k = ⌊m / k⌋` (blocks completed within m minutes):

1. (Weighted maximum.) The monitor that alarms when `max_k w_k M^(k)_t ≥ 2A` on either side, checked every minute with frozen values, has `P(false alarm within m) ≤ Σ_k w_k N_k / A`.
2. (Plain average.) The monitor that alarms when `Σ_k w_k M^(k)_t ≥ 2A` on either side, checked every minute with frozen values, has `P(false alarm within m) ≤ (Σ_k √(w_k N_k))² / A`.

*Proof.* For each k let `D^(k) = (M^(k↑) + M^(k↓)) / 2`. On the block clock it is an e-detector (Lemma 1′ for each side), so by Theorem B of `docs/math.md` `P(max_{b ≤ N_k} D_b^(k) ≥ x) ≤ N_k / x` for every x > 0. The frozen value at any minute up to `T + m` is the value after one of the first `N_k` blocks, and either side reaching `2y` implies `D^(k) ≥ y`.
1. An alarm implies `M^(k, side) ≥ 2A / w_k` for some k and side, so `D^(k) ≥ A / w_k`. A union bound over k gives `Σ_k N_k w_k / A`.
2. Fix any `a_k > 0` with `Σ_k w_k a_k = A`. If `Σ_k w_k M^(k, side) ≥ 2A`, then `M^(k, side) ≥ 2 a_k` for some k (otherwise the weighted sum would be below `2 Σ_k w_k a_k = 2A`), so `D^(k) ≥ a_k`. A union bound gives `Σ_k N_k / a_k`. By the Cauchy-Schwarz inequality, `(Σ_k √(w_k N_k))² = (Σ_k √(w_k a_k) √(N_k / a_k))² ≤ (Σ_k w_k a_k)(Σ_k N_k / a_k) = A Σ_k N_k / a_k`, with equality for `a_k ∝ √(N_k / w_k)`; that choice gives the bound. ∎

For balanced weights and m a multiple of every k, `(Σ_k √(w_k N_k))² = J × Σ_k w_k N_k`: the plain average needs a J times larger A than the certified average for the same proved budget, which would make its single-frequency handicap `J²`. E17 runs the plain average at the certified level, where it has no proof, to see whether the proof is merely conservative.

## 7. Recalibration

After every alarm the monitor recalibrates all frequencies on the next calibration window, and all block grids restart at the new monitoring start.

**Corollary 4 (restarts).** Suppose that in every segment the nulls of Theorem 1 hold relative to that segment's own baselines (no change anywhere, true baselines). Let `N_m` be the number of alarms of the certified average within the first m *monitored* minutes (calibration minutes are not counted). Then `E[N_m] ≤ (m ρ + κ) / (A - κ)`.

*Proof.* Let `W_t` = (sum of `S` at all past alarms) + (the current `S`, while monitoring) - (the weighted number of blocks started in monitored minutes so far, summed over segments). A segment starts at a stopping time and its `S` starts at 0, so by the proof of Theorem 1 every monitored minute has `E[ΔW | F_{t-1}] ≤ 0`. At an alarm the current `S` moves into the sum unchanged (an abandoned partial block was already counted when it started), and during calibration nothing changes, so W is a supermartingale with `W_0 = 0`. Each alarm banks `S ≥ A` (Lemma 2). Let σ be the minute of the m-th monitored minute. Optional stopping at `σ ∧ n` and monotone convergence give `A E[N_m] ≤ E[blocks started by σ]`. A segment with ℓ monitored minutes starts `Σ_k w_k ⌈ℓ / k⌉ ≤ ℓ ρ + κ` weighted blocks, and at most `N_m + 1` segments have monitored minutes among the first m. So `A E[N_m] ≤ m ρ + κ (E[N_m] + 1)`, and `N_m ≤ m` is finite, which rearranges to the claim. ∎

We do not claim a restart bound for the plain average or the weighted maximum: when the monitor restarts in the middle of a coarse block because another frequency alarmed, the decision to abandon the block depends on that block's own returns, which is the problem of section 3 again.

## 8. Estimated baselines

The theorems assume that each frequency's τ and baseline are the true values for the null. The code estimates them from the same calibration window (plug-in), so, as for the one-minute monitor (`docs/math.md` 5.6, item 1, and Proposition A.6 of `docs/theory_running_tau.md`), the guarantee is then approximate. Coarse frequencies have fewer independent observations in the window (96 non-overlapping fifteen-minute blocks in a day), so their plug-in error is larger; the code uses every window of the block length in the calibration data (all phases of the block grid), which reduces but does not remove this. A confidence-sequence baseline for each frequency, with error `δ_k` and τ fixed from a first part of the window, should give `Σ_k δ_k + Σ_k w_k ⌈m / k⌉ / A` by the argument of Theorem 6.4 of `docs/math.md` applied on the event that all J bounds hold; this is not implemented or checked here.

## 9. Status

* **Proved here:** Lemma 1′, Proposition 2, Lemma 2, Theorem 1, Proposition 3, Theorem 3 and Corollary 4, given Ville's inequality and optional stopping at bounded stopping times; the implications of section 2.1 and the gap-block independence under i.i.d. price noise.
* **Cited:** averaging e-values (Vovk and Wang 2021); e-detectors (Shin, Ramdas and Rinaldo 2024); Ville's inequality (`docs/math.md` section 6); lifting e-processes across filtrations (Choe and Ramdas), mentioned for context only.
* **Checked numerically** (`tests/test_multifreq.py`): the streaming monitor equals the weighted sum of separately computed detectors at every minute; the vectorised code used in E17 equals the streaming monitor; k = 1 alone reproduces `VolatilityMonitor` event for event; first-alarm rates on i.i.d. streams are below all three bounds; alarm counts with restarts are below Corollary 4; Proposition 2, part 2; the block correlations under price noise.
* **Assumed, not proved:** the nulls on real data (as for the one-minute monitor, volatility clustering violates them); the plug-in baselines.

## References

* Choe, Y. J. and Ramdas, A. Combining evidence across filtrations. arXiv:2402.09698 (version 5, March 2026), listed there as accepted in the *Journal of the Royal Statistical Society, Series B*; volume and pages not found. Only the search listing and abstract-level descriptions could be read (arxiv.org is not reachable from this machine).
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science* 2(2), 229 to 260. doi:10.51387/23-NEJSDS51.
* Vovk, V. and Wang, R. (2021). E-values: calibration, combination and applications. *Annals of Statistics* 49(3), 1736 to 1754. doi:10.1214/20-AOS2020. Bibliographic details checked against the Project Euclid and Royal Holloway listings found by web search; the full text could not be opened here (arxiv.org and the authors' site are not reachable from this machine). We use only its basic points: the arithmetic mean of e-values is an e-value whatever their dependence, and among symmetric ways of merging e-values the arithmetic mean essentially dominates the others.
