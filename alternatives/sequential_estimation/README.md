# Sequential fixed-accuracy estimation of a volatility level (alternative ALT6)

A research prototype that asks a classical question of sequential analysis about one-minute returns: **how many minutes of data do we need before we can state the current volatility level to a prescribed accuracy, with a prescribed confidence, if we decide when to stop by looking at the data as they arrive?** It compares the classical stopping rules (Chow and Robbins 1965, its relative-accuracy version, Stein's two-stage rule, Robbins' minimum risk rule) on squared returns and on MAD-truncated squared returns (the truncation of Bapat 2026), with an anytime-valid confidence sequence that is simply stopped when it is narrow enough. Every number quoted from an experiment is in `RESULTS.md`, labelled [SYNTHETIC] or [REAL], and comes from `run_experiment.py`.

Contents: 1 problem, 2 mathematical formulation (with derivations), 3 why it matters, 4 comparison with existing approaches, 5 implementation, 6 experimental methodology, 7 folder structure and how to run, references.

## 1. Research problem

One-minute returns `r_1, r_2, ...` arrive one at a time. After dividing by a predictable intraday profile (`avmonitor.data.intraday_scale`, built from earlier days only), we treat the returns inside a window of a few hours as independent draws from one law. The volatility level is a number that summarises the scale of that law, for example the variance `θ = E r²`.

A **sequential procedure** looks at the returns one by one and, after each one, decides whether to stop. Its **stopping time** `N` is the random number of returns it has seen when it stops; the decision to stop at `n` may use only the first `n` returns. At `N` it reports an estimate `θ̂_N`. The goal is **fixed accuracy**:

* absolute accuracy: `P(|θ̂_N - θ| ≤ d) ≥ 1 - α`, or
* relative accuracy: `P(|θ̂_N - θ| ≤ d θ) ≥ 1 - α`, which says "within plus or minus 100d percent" and does not depend on the unit of the returns.

The probability `P(|θ̂_N - θ| ≤ d)` is called the **coverage** at the stopping time. A fixed sample size cannot do this when the spread of the data is unknown, because the right sample size depends on that spread. The sequential answer, sixty years old, is to estimate the spread on the fly and stop when the estimated required size is reached.

Three facts make the volatility version interesting:

1. The statistic is `X = r²`, whose own variance `ν = Var(r²) = E r⁴ - (E r²)²` involves the fourth moment, the **quarticity**. Any fixed-accuracy rule for volatility is therefore a rule driven by a running estimate of quarticity, which is the quantity Bapat (2026) estimates robustly.
2. Returns are heavy tailed and jumpy. For Student t returns with 3 degrees of freedom `ν` is infinite, and a single jump can multiply the sample quarticity by ten. Classical rules come with guarantees only as `d → 0` and with finite variance; what happens at the accuracies a trader would ask for (plus or minus 20% to 50%) is an empirical question.
3. A **confidence sequence** (Howard et al. 2021; Waudby-Smith and Ramdas 2024) is a sequence of intervals that contain the target at all times simultaneously with probability `1 - α`. Stopping it the first time its width is at most `2d` gives a fixed-accuracy procedure that is valid at every `d`, with no asymptotics. The price is extra data. How much extra, and why, is the central quantitative question of this prototype.

## 2. Mathematical formulation

### 2.1 Targets and the optimal fixed sample size

Two statistics are studied.

* **Raw:** `X_i = r_i²`, target `θ = E r²` (the per-minute variance), spread `ν = Var(r²)`.
* **Truncated:** `X_i = min(r_i², τ²)` with `τ = c × (scaled MAD of a pilot sample)`, `c = 3`. The **pilot** is the first `m` returns (default `m = 30`); the scaled MAD is `1.4826 × median |r|` over the pilot (the median absolute deviation about zero, scaled to equal σ for normal returns). The target is the **truncated variance** `θ_τ = E[min(r², τ²) | τ]` at the realised `τ`, and its spread `ν_τ = E[min(r⁴, τ⁴) | τ] - θ_τ²` is a truncated quarticity: up to a constant, the averaged robust kernel estimator of Bapat (2026) estimates exactly this (`docs/estimator_issues.md`, Issue 3).

If `ν` were known, a fixed sample of size

`n* = z² ν / d²` (absolute) or `n* = z² ν / (d θ)² = z² CV² / d²` (relative),

with `z = z_{1-α/2}` the normal quantile and `CV² = ν / θ²` the squared coefficient of variation, would have coverage close to `1 - α` by the central limit theorem. `n*` is the benchmark: every procedure is judged by `E[N] / n*` (the **efficiency ratio**) and by `E[N] - n*` (the **regret**, extra observations on average).

For normal returns `r = σ Z` and `τ = a σ`, with `Φ̄` the normal upper tail and `φ` the normal density,

* `κ₂(a) = E[min(Z², a²)] = 1 - 2Φ̄(a) - 2aφ(a) + 2a²Φ̄(a)`,
* `κ₄(a) = E[min(Z⁴, a⁴)] = 3(1 - 2Φ̄(a) - 2aφ(a)) - 2a³φ(a) + 2a⁴Φ̄(a)`,

so `θ_τ = σ² κ₂(a)` and `CV²_τ = κ₄(a) / κ₂(a)² - 1`. Raw squared normal returns have `CV² = 2`; at `a = 3`, `κ₂ = 0.99501` and `CV²_τ = 1.92058`, so the truncated statistic needs 4% fewer observations for its own target (`seqest.kappa2`, `seqest.kappa4`, both checked against the package's truncated-moment functions and Monte Carlo in the tests). For Student t5 returns raw `CV² = 8` (kurtosis 9) against truncated 2.41; for t3, raw `CV²` is infinite against truncated 2.51.

### 2.2 No rule can be exactly valid for the raw variance

**Proposition 1.** Let `P0` be any law on `[0, ∞)` for `X = r²` (for example the law of a squared normal return) and let `(N, θ̂_N)` be any sequential procedure (a stopping rule and an estimate computed from `X_1, ..., X_N`, possibly with extra independent randomisation) that stops with probability one under `P0`. Then for every `d > 0` and `η > 0` there is a law `P` with bounded support (so every moment is finite) such that

`P0(|θ̂_N - θ(P0)| ≤ d) + P(|θ̂_N - θ(P)| ≤ d) ≤ 1 + η`.

In particular, if the coverage under `P0` is at least `1 - α`, the coverage under `P` is at most `α + η`. The same holds for relative accuracy with `d < 1`.

*Proof.* Pick `n0` with `P0(N > n0) ≤ η`, a small `ε > 0`, and a point `M` that is not an atom of `P0`, with `εM = 3d + εθ(P0)` (relative accuracy: `εM` large enough that `(1 - d) θ(P) > (1 + d) θ(P0)`). Put `P = (1 - ε) P0 + ε δ_M`, so `θ(P) = θ(P0) + 3d`. Let `A` be the event that none of `X_1, ..., X_{n0}` equals `M`, and `E = {N ≤ n0, |θ̂_N - θ(P0)| ≤ d}`. Given `A`, the first `n0` observations under `P` are i.i.d. `P0`, and `E` depends only on them, so `P(E ∩ A) = (1 - ε)^{n0} P0(E) ≥ (1 - ε)^{n0} (P0(cover) - η)`. On `E`, `|θ̂_N - θ(P)| ≥ 3d - d > d`, a miss under `P`. Hence `P(cover) ≤ 1 - (1 - ε)^{n0} (P0(cover) - η)`, and letting `ε → 0` (with `n0` fixed) gives the claim. ∎

This is the sequential form of the impossibility result of Bahadur and Savage (1956) and the same contamination argument as `docs/math.md` 5.3. Its message: Chow and Robbins' guarantee for the raw variance is, and must be, a statement about `d → 0` for one fixed law, never a guarantee at a given `d` over a class of laws. A rare large return that the procedure has not yet seen moves the variance without moving the data seen so far. **Truncation removes the obstruction**: `min(r², τ²) ≤ τ²`, so contamination of mass ε moves `θ_τ` by at most `ετ²`, and section 2.8 gives a procedure that is exactly valid for `θ_τ` at every `d`, for every law of the returns.

### 2.3 The Chow and Robbins rule

With `X̄_n` and `ν̂_n` the running mean and sample variance of `X_1, ..., X_n`,

`N = inf{ n ≥ m : n ≥ (z / d)² (ν̂_n + 1/n) }`, report `X̄_N` (interval `X̄_N ± d`).

The rule says: stop as soon as the sample size passes the estimated `n*`. The `1/n` term stops the rule from firing early when `ν̂_n` happens to be tiny. Chow and Robbins (1965) proved, for i.i.d. data with `0 < ν < ∞`, as `d → 0`: `N / n* → 1` almost surely, `E[N] / n* → 1` (**asymptotic efficiency**) and coverage `→ 1 - α` (**asymptotic consistency**). In this prototype the pilot returns are part of the sums (the rule is first checked at `n = m`).

*Relative accuracy* (Nadas 1967): `N = inf{ n ≥ m : n ≥ (z/d)² (ν̂_n + 1/n) / X̄_n² }`, report `X̄_N ± d X̄_N`. The `1/n` term has the units of `ν`, so this rule is not unit-free. We make it unit-free by dividing every return by the pilot's scaled MAD first (a predictable normalisation); all relative-accuracy runs do this.

**Skewness and the stopping bias.** The rule stops early when `ν̂_n` is low, and for squared returns a low `ν̂_n` comes with a low `X̄_n` (their correlation is `μ₃ / √(ν (μ₄ - ν²))`, about 0.76 for squared normal returns, where `μ₃` and `μ₄` are the third and fourth central moments of X). A first-order calculation: Wald's identity gives `E[S_N] = 0` for `S_n = Σ_{i≤n} (X_i - θ)`, so `E[X̄_N - θ] = E[S_N / N] ≈ -E[S_N (N - n*)] / n*²`. In the absolute rule `N - n* ≈ n* (ν̂ - ν) / ν` and `E[S_{n*} n* (ν̂_{n*} - ν)] ≈ n* μ₃`, so

`E[X̄_N] - θ ≈ -μ₃ / (ν n*)` (absolute rule).

In the relative rule `N - n* ≈ n* [(ν̂ - ν)/ν - 2(X̄ - θ)/θ]`, which adds `+2ν/θ` and gives

`E[X̄_N] - θ ≈ -(μ₃ - 2ν²/θ) / (ν n*)` (relative rule).

For any gamma law, squared normal returns included, `μ₃ = 2ν²/θ` and the first-order bias of the relative rule is zero. These are our own first-order calculations (second-order theory is in Woodroofe 1977 and Aras and Woodroofe 1993, which we did not read in full); the simulations in `RESULTS.md` show the predicted asymmetry of misses.

**Infinite quarticity.** Under t3 returns `P(r² > x)` decays like `x^{-3/2}`, so `Σ r⁴` grows like `n^{4/3}` and `ν̂_n` like `n^{1/3}`. The rule `n ≥ C d^{-2} n^{1/3}` then stops at `N ≈ d^{-3}` instead of `d^{-2}`, and the studentised mean converges to a non-normal limit (self-normalised sums outside the domain of attraction of the normal law: Logan, Mallows, Rice and Shepp 1973; Giné, Götze and Mason 1997). The coverage of the raw rule then has no reason to tend to `1 - α`; in our simulations it falls from 0.78 to 0.72 as `d` goes from 0.4 to 0.05 (`RESULTS.md`, R3).

### 2.4 The truncated statistic: what it buys and what it costs

Truncation gives a bounded statistic with the same rules (`CR-trunc`, `Stein-trunc`) and makes the confidence sequence possible. Its target is `θ_τ`, not `σ²`.

**Mapping back to σ².** If the returns are `r = σ ε` with a *known* unit-variance shape `ε` (for example normal), then `θ_τ = φ_τ(σ²) = σ² g(τ/σ)` with `g(a) = E[min(ε², a²)]`. Since `dφ_τ/dσ² = g(a) - (a/2) g'(a) = E[ε² 1{|ε| ≤ a}] > 0`, the map is increasing and can be inverted (`seqest.sigma2_from_theta_tau`, Brent's method). Monotone maps carry intervals to intervals, so an interval for `θ_τ` with coverage `1 - α` becomes one for `σ²` with the same coverage, *if the shape is right*. Two costs:

* *Accuracy.* The elasticity `e(a) = d log θ_τ / d log σ² = E[ε² 1{|ε| ≤ a}] / E[min(ε², a²)]` is 0.976 at `a = 3` for normal returns (`seqest.elasticity_gaussian`), so a relative accuracy `d` on `θ_τ` is about `d / 0.976` on `σ²`. In sample size: the effective `CV²` for `σ²` is `1.9206 / 0.976² = 2.018` against 2 for the raw statistic, a loss of 0.9%. At `c = 2` the elasticity is 0.80 and the loss is larger (`CV²` 1.454 / 0.802² = 2.26).
* *Shape.* If the shape is wrong, the mapped value is biased. Mapping with a normal shape gives, at `c = 3`, 0.877 `σ²` for t5 returns, 0.621 `σ²` for t3 returns, and 1.055 times the diffusive variance (0.844 times the total variance) for normal returns with 1% Bernoulli jumps of 5 sd (`results/g_mapping.csv`, closed form). This is not a defect of the procedures: the truncated variance measures the core of the distribution, and the tail part of `σ²` is exactly what Proposition 1 says cannot be pinned down at a given accuracy.

### 2.5 Stein's two-stage rule

Take the pilot of `m` values, compute its sample variance `s_m²`, and sample up to

`N = max(m, ⌊(t_{m-1} / d)² s_m²⌋ + 1)` (relative: `d` replaced by `d X̄_m`),

then report the mean of all `N` values. For normal X, `√N (X̄_N - θ) / s_m` has exactly a Student t law with `m - 1` degrees of freedom, so the coverage is at least `1 - α` exactly (Stein 1945). For squared returns X is far from normal, so the guarantee does not apply. Also, for the absolute rule, `E[N] / n* ≥ (t_{m-1}/z)²` (1.09 at `m = 30`) and `N` is very variable, because `s_m²` from 30 squared returns is a poor estimate of `ν`. It is included as the cheapest classical baseline.

### 2.6 Robbins' minimum risk rule

Point estimation with loss `A (X̄_n - θ)² + n` (squared error weighted by A, plus a cost of one per observation). For known `ν` the risk `A ν / n + n` is smallest at `n0 = √(A ν)`, with risk `2 n0`. Robbins (1959) replaced `ν` by its estimate:

`N = inf{ n ≥ m : n ≥ √A (s_n + n^{-γ}) }`, with `s_n = √ν̂_n` and `γ = 1`.

Performance measures: **risk efficiency** `E[loss] / (2 n0)` (tends to 1) and **regret** `E[loss] - 2 n0` (bounded under conditions; Starr 1966; Woodroofe 1977 for normal data; Ghosh and Mukhopadhyay 1979 and Chow and Yu 1981 beyond normal data). Here `X = r²` or `min(r², τ²)`, so the rule is a minimum risk estimator of the variance or of the truncated variance.

### 2.7 Absolute or relative accuracy

Volatility has no natural unit, so relative accuracy is the natural goal ("know the variance to within plus or minus 20%"). Relative accuracy on the variance `d` corresponds to about `d / 2` on the volatility `σ = √θ`. All main experiments use relative accuracy; one part (`abs`) repeats them in absolute accuracy with `σ² = 1`, where the skewness bias of 2.3 is not cancelled.

### 2.8 Anytime-valid confidence sequences stopped at the target width

The pilot fixes `τ` (and, for relative accuracy, the unit). The sequence then runs on the later returns only (**sample splitting**, `docs/math.md` 6.4), on the bounded variable `Z_i = min(r_i², τ²) / τ² ∈ [0, 1]` with mean `μ_Z = θ_τ / τ²`. Given the pilot, the `Z_i` are i.i.d. under the working model.

A **confidence sequence** is a pair of processes `L_t ≤ U_t`, each computed from `Z_1, ..., Z_t`, with `P(L_t ≤ μ_Z ≤ U_t for all t ≥ 1) ≥ 1 - α`. Because the guarantee holds at all times at once, it also holds at any random time, in particular at a time chosen by looking at the data. Replacing `[L_t, U_t]` by the running intersection `[max_{s≤t} L_s, min_{s≤t} U_s]` keeps the guarantee when the target is constant (on the good event every interval contains `μ_Z`, so the intersection does too); section 2.10 explains why this fails when the conditional mean moves.

**Stopping rules.** Absolute: `T = inf{t : U_t - L_t ≤ 2 d / τ²}`, report the midpoint `θ̂ = τ² (L_T + U_T) / 2`. Relative: `T = inf{t : L_t > 0 and U_t - L_t ≤ d (U_t + L_t)}`, report the harmonic mean `θ̂ = τ² · 2 L_T U_T / (L_T + U_T)`. The total number of returns used is `N = m + T`.

**Theorem 1 (validity at the stop).** For every law of the returns, every `d > 0`, every pilot size and every `α`, given the pilot,

`P(T < ∞ and |θ̂ - θ_τ| > d) ≤ α` (absolute) and `P(T < ∞ and |θ̂ - θ_τ| > d θ_τ) ≤ α` (relative).

*Proof.* On the event that `μ_Z ∈ [L_t, U_t]` for all t (probability at least `1 - α`), at `T` we have `μ_Z ∈ [L_T, U_T]`. Absolute: an interval of width at most `2d/τ²` containing `μ_Z` puts its midpoint within `d/τ²` of `μ_Z`. Relative: for any `μ ∈ [L, U]`, `θ̂/μ ≤ 2LU / ((L+U) L) = 2U/(L+U) ≤ 1 + d` and `θ̂/μ ≥ 2L/(L+U) ≥ 1 - d`, and both inequalities are the stopping condition `U - L ≤ d (U + L)` rewritten. ∎

`T` can be infinite if the width never falls below the target (a constant bet gives a width that does not shrink to zero); the experiments cap the sample at `n_max` and report how often the cap is reached. No moment condition appears: a t3 law or a jump law is covered. What the theorem does not give is a statement about `σ²` (section 2.4) or about dependent data (section 2.10).

Two constructions are implemented.

**(a) Empirical Bernstein sequences** (`CS-wsr`, `CS-tuned`, `CS-opt`). These are the predictable plug-in empirical Bernstein bounds of Waudby-Smith and Ramdas (2024) with the one-sided ranges and the proof of `docs/math.md` 6.2: with a predictable centre `ẑ_{i-1}`, `e_i = (Z_i - ẑ_{i-1})²`, ranges `a↑ = 1 - ẑ`, `a↓ = ẑ`, bets `λ_i` with `λ_i a_i ≤ 0.75` and `ψ_a(λ) = (-log(1 - λa) - λa) / a²`,

`U_t = [Σ λ↑_i Z_i + log(2/α) + Σ ψ↑_i e_i] / Σ λ↑_i`, `L_t = [Σ λ↓_i Z_i - log(2/α) - Σ ψ↓_i e_i] / Σ λ↓_i`

(each side at level `α/2`). The running centre and variance start from the pilot's `Z` values with weight `m` (known before the sequence starts, so still predictable). The bet schedules differ:

* `wsr`: the time-uniform schedule `λ_t = √(2 log(2/α) / (σ̂²_{t-1} t log(1 + t)))` of Waudby-Smith and Ramdas, with their start (centre 1/2, variance 1/4). With the intersection switched off, the code reproduces `avmonitor.confseq.eb_confseq` exactly (a test).
* `tuned`: `λ = d_Z / σ̂²` (capped), with `d_Z` the target half-width in Z units; section 2.9 shows that this is optimal as `d → 0`.
* `opt`: the constant bets per side that minimise the predicted stopping time with the exact `ψ` and the caps (`seqest.opt_plan`), re-solved from the running estimates at the start of every block of 128 steps.

**(b) A betting sequence** (`CS-bet`). For each candidate mean `μ` the capital processes

`K↑_t(μ) = Π_{i≤t} (1 + λ↑_i (μ - Z_i))`, `K↓_t(μ) = Π_{i≤t} (1 + λ↓_i (Z_i - μ))`

are nonnegative martingales when `μ` is the true mean (each factor has conditional mean one), so by Ville's inequality `P(K_t(μ_Z) ≥ 2/α for some t) ≤ α/2` for each side. The confidence set is every `μ` that neither side has rejected. Implementation details that matter:

* One bet per side, shared by all candidate means, chosen predictably once per block: the growth-rate-optimal constant bet against `ẑ ± w` (with `w` the target half-width) under the empirical law of the Z values seen so far, binned relative to the pilot mean (bins of width `ẑ_0/8` up to `10 ẑ_0`, one overflow bin, one bin for `Z = 1`; per-bin sums keep the binned mean exact).
* Because the bet does not depend on `μ`, `F↑(μ) = log K↑_t(μ)` is nondecreasing and concave in `μ`, and `|F↑''(μ)|` is nonincreasing; the lower side is the mirror image. The candidate means are a fine log grid of 96 points on `[ẑ_0/5, 5 ẑ_0]` (the core) plus 8 log-spaced points down to `10⁻⁶` and 32 up to `1 - 10⁻⁶`. Between two grid points the code uses the larger of two valid lower bounds on `F`: the chord (concavity) and the second-order bound `F(m_0) + F'(m_0) y - |F''(m_0)| y²/2` with the curvature taken at the end where it is largest. Every `μ` where a lower bound on `F` reaches `log(2/α)` is rejected, so the reported bounds are conservative. A test checks them against exact root finding: never inside, and at most 0.05% wider.
* Bets are capped so that every factor stays positive: `λ↑ ≤ 0.75 / (1 - ẑ_0/5)`, which keeps `λ↑ (1 - μ) < 1` for every `μ ∈ [0, 1)`, and `λ↓ ≤ 0.75 / (5 ẑ_0)` on the core. Above the core the lower bet is capped per point, `min(λ↓, 0.75/μ)`. Each term `log(1 + min(λ, 0.75/μ)(Z - μ))` is still nonincreasing in `μ`, so the rejected set is still an interval and the largest rejected grid point is a valid lower bound; the terms with a capped bet are convex rather than concave, so above the core the code does not interpolate. The extensions matter only when the pilot is far off: the first version had a core-only grid, and with a 10-return pilot 3 of 2000 normal paths had their true mean above the grid, so the lower bound could never pass it and the rule never stopped (`RESULTS.md`, R6). With the extensions the same paths stop at most 4 times `n*`, and paths whose bounds stay in the core are unchanged.

### 2.9 How much data does the confidence sequence need?

**Without intersection, as `d → 0`.** Write `L = log(2/α)` and `σ_Z²` for the variance of Z. When the bets are small, `ψ_a(λ) ≈ λ²/2`, and with a constant bet `λ` the upper bound is about `Ẑ_t + L/(λt) + λσ_Z²/2`. The half-width `L/(λt) + λσ_Z²/2` reaches `d_Z` first, at `t = 2 L σ_Z² / d_Z²`, when `λ = d_Z / σ_Z²` (the `tuned` bet). Dividing by `n* = z² σ_Z² / d_Z²`:

`E[N] / n* → 2 log(2/α) / z²` = **1.9206** at `α = 0.05` and **2.2145** at `α = 0.10` (`seqest.cs_ratio_limit`).

(At `α = 0.05` this happens to equal, to four digits, the `CV²_τ = 1.9206` of section 2.1; the two numbers are unrelated.) The reading: a fixed-sample interval built from the central limit theorem has half-width `z σ/√n`; the sequence has `√(2L) σ/√n`, which is the half-width of a fixed-sample Chernoff (exponential) bound at `α/2` per side. At a target width known in advance, Ville's inequality applied to a linear boundary is tight at one time, so **time-uniformity itself costs nothing in the limit; the factor 1.92 is the price of an exponential tail bound instead of a normal quantile.** At the stop, the midpoint is approximately normal with standard deviation `σ_Z / √t`, so the coverage there tends to `1 - 2Φ̄(√(2L))` = 0.9934 (`α = 0.05`) or 0.9856 (`α = 0.10`): the sequence over-covers at its stopping time. With a bet `λ = κ d_Z / σ_Z²` the limit is `2L / (κ(2 - κ) z²)`, smallest at `κ = 1`.

**With the running intersection, as `d → 0`.** Measure time as `u = t / n*` and the bounds in units of `d_Z`. By Donsker's theorem the bounds converge to

`U(u) = W(u) / (z u) + L / (κ z² u) + κ/2`, `L(u) = W(u) / (z u) - L / (κ z² u) - κ/2`

with W a standard Brownian motion, and the intersected sequence stops at `u* = inf{u : min_{v≤u} U(v) - max_{v≤u} L(v) ≤ 2}`. `E[u*]` has no closed form we know of; `seqest.brownian_limit` computes it by simulation (step 0.001, 40,000 paths): **1.548** at `α = 0.05` and **1.708** at `α = 0.10` for `κ = 1` (`results/a_brownian_limit.csv`). The intersection helps because the upper bound keeps its luckiest past value and the lower bound its luckiest past value, at different times. Coverage at the stop in this limit is 0.976 and 0.952. The bet that is optimal without intersection is not optimal with it: `κ = 1.2` gives 1.475 at `α = 0.05`.

**At a practical d.** At `d = 0.2` the tuned bet hits the cap (`λ a ≤ 0.75`) and `ψ` is far from quadratic (`ψ_a(λ) / λ² = h(λa)` with `h(x) = (-log(1 - x) - x)/x²`, which is 1.13 at `x = 0.75` against 0.5 for small x). Replacing every running estimate by its limit and keeping the exact `ψ` and caps gives the first-order prediction (`seqest.cs_n_predict`)

`t ≈ L (1/λ↑ + 1/λ↓) / (2 d_Z - σ_Z² [λ↑ h(λ↑ a↑) + λ↓ h(λ↓ a↓)])`

(relative accuracy: `(1 - d)` and `(1 + d)` weights from linearising `(U - L) ≤ d(U + L)`), and `opt_plan` minimises it over the two bets. For normal returns at `c = 3`, `α = 0.05` it predicts `E[N]/n*` of 2.43 (`d = 0.2`), 2.20 (0.1), 2.07 (0.05), tending to 1.92 (`results/a_closed_form.csv`).

**The best a betting sequence can do for a known law.** For a known law of Z, the best constant bet against `μ ± w` grows the log-capital at rate `g↑(w) = max_λ E log(1 + λ(μ + w - Z))` per step (and `g↓` for the lower side), so the side rejects `μ ± w` after `L / g(w)` steps; `seqest.betting_oracle_n` finds the smallest t at which both sides together meet the stopping rule. This "betting oracle" is a benchmark, not a lower bound for every sequence (the intersected sequences can beat it), and it also tends to 1.92 as `d → 0`. For normal returns at `c = 3` it is 1.80 (`d = 0.4`), 1.89 (0.2), 1.91 (0.1).

### 2.10 Dependent returns: volatility clustering

Under GARCH(1,1) clustering, `r_t = √h_t ε_t`, `h_t = ω + a r_{t-1}² + b h_{t-1}` with `ω = 1 - a - b` (unit variance), squared returns are positively autocorrelated: kurtosis `3(1 - (a+b)²) / (1 - (a+b)² - 2a²)` (Bollerslev 1986), `Corr(r_t², r_{t-k}²) = ρ₁ (a+b)^{k-1}` with `ρ₁ = a (1 - ab - b²) / (1 - 2ab - b²)` (Bollerslev 1988). The variance of a sample mean of `n` squared returns is then about `ν_LR / n` with the **long-run variance** `ν_LR = ν (1 + 2 Σ_k ρ_k) = ν (1 + 2ρ₁ / (1 - a - b))`: a factor 8.16 for `(a, b) = (0.10, 0.85)` and 32.0 for `(0.05, 0.94)` (`seqest.Garch.lrv_factor_raw`, checked by simulation in the tests).

* The Chow and Robbins rule uses the i.i.d. variance `ν̂`, so it stops about 8 to 32 times too early for the stationary variance. The fix implemented here is **batch means** (`CR-BM`): split the data into consecutive batches of `B` returns and use `B` times the sample variance of the batch means, which estimates `ν_LR` when B and the number of batches both grow; the rule is checked at batch ends after at least 10 batches.
* The confidence sequence needs no independence (Theorem 6.2 of `docs/math.md`), but then it covers the **bet-weighted average of the conditional truncated variances** `E[min(r_i², τ²) | past]` over the returns it has seen, not the stationary level. Both targets are reported: coverage of the stationary `θ_τ` (computed for each `τ` from long simulated paths with batch means of length 5000, `run_experiment.garch_targets`) and coverage of the theorem's own target (computed exactly from the simulated conditional variances, using `E[min(h ε², τ²) | h] = h κ₂(τ / √h)`).
* **What the classical rule covers under clustering.** Write `X_t = min(r_t², τ²)` and `μ_t = E[X_t | past]`. The sum `Σ (X_t - μ_t)` is a martingale, so by the martingale central limit theorem `X̄_n - μ̄_n` (with `μ̄_n` the average of the `μ_t` over the sample) is approximately normal with variance `(1/n²) Σ Var(X_t | past)`. The i.i.d. variance estimate `ν̂` estimates the average of `E[(X_t - μ̄)²]`, which is the average conditional variance plus the spread of the `μ_t` around their mean, so it is at least as large as what is needed. Hence the classical interval, with `ν̂`, is asymptotically valid and slightly conservative for the **window average** `μ̄_N`, an analogue of integrated variance (Barndorff-Nielsen and Shephard 2002), even though it is far too short for the stationary level. Results R7: coverage of the window average 0.948 to 0.966.
* **The running intersection is not valid for a moving target.** Intersecting the bounds over time is justified only when every bound in the sequence is about the same number. With clustering the target of the bound at time `t` is an average over the first `t` returns, which changes with `t`; an old, narrow bound can exclude the current average. Under GARCH the intersected sequences cover the window average at only 0.68 to 0.79 at `d = 0.1`, while the plain sequences (`intersect=False`; for the empirical Bernstein sequence this is what `avmonitor.confseq.eb_confseq` computes) cover it at 0.986 to 1.000. Under clustering, use the sequences without intersection.

## 3. Why it matters

* **For a trader or risk system** the question "how long must I wait after a regime change before I can trust a volatility number to plus or minus 20%?" has a direct answer in minutes, by hour of day, and the answer depends on the procedure: classical rules answer it only asymptotically, and heavy tails make them stop too early.
* **For sequential analysis** the problem joins two literatures that rarely meet: fixed-accuracy stopping rules (an active field since the 1960s) and anytime-valid inference by betting (an active field since about 2020). The comparison is concrete: the same data, the same accuracy target, an asymptotic guarantee against an exact one, and a closed-form account of the price.
* **For high-frequency econometrics** the truncated variance and truncated quarticity are standard objects (truncated realised variance, Bapat 2026), but stopping rules for them are not; and the impossibility result (Proposition 1) says that for the raw variance no stopping rule can be exactly valid at a given accuracy.

## 4. Comparison with existing approaches

| Approach | Target | Guarantee at a given d | Needs | Expected sample size |
|---|---|---|---|---|
| Fixed sample of `n*` (oracle) | any mean | asymptotic (CLT), needs known `ν` | `ν` | `n*` by definition |
| Chow and Robbins (1965), Nadas (1967) | mean of X | none; coverage `→ 1 - α` as `d → 0` | finite variance | `E[N]/n* → 1` |
| Stein (1945) two-stage | mean of X | exact for normal X only | normal X | `≥ (t_{m-1}/z)² n*` |
| Robbins (1959) minimum risk | mean of X, point estimate | risk efficiency `→ 1` | finite moments (Chow and Yu 1981) | `≈ n0` |
| Empirical Bernstein confidence sequence stopped at width `2d` | bounded mean (`θ_τ`) | **exact**, any law | bounded X (truncation) | `→ 1.92 n*`, 1.55 `n*` with intersection (`α = 0.05`) |
| Betting confidence sequence stopped at width `2d` | bounded mean (`θ_τ`) | **exact**, any law | bounded X | close to the intersected limit already at moderate d |
| Asymptotic confidence sequences (Waudby-Smith, Arbour, Sinha, Kennedy and Ramdas 2024) | mean | asymptotic, time-uniform | finite moments | not studied here |
| Realised variance with a feasible CLT (Barndorff-Nielsen and Shephard 2002) | integrated variance over a fixed window | asymptotic in the sampling frequency | quarticity estimate | fixed window, no stopping |
| Renewal and duration-based sampling (Li, Nolte and Nolte 2019) | integrated variance | asymptotic | stopping times of the price, not of the estimate | not a fixed-accuracy rule |
| Sequential fixed-accuracy estimation for AR and threshold AR/ARCH models (Lai and Siegmund 1983; Konev and Nazarenko 2020; Vorobeychikov, Burkatovskaya and Sergeeva 2016) | regression or autoregression parameters | non-asymptotic for some designs | a parametric model | model specific |

**What we searched for and did not find.** We searched (web search only; projecteuclid.org and api.crossref.org could not be fetched from this machine) for fixed-width or fixed-precision stopping rules for volatility from high-frequency data, for sequential estimation of realised or integrated variance, and for confidence sequences applied to financial volatility. We found none. The nearest work is listed in the table: renewal-time sampling, which uses stopping times of the price path to build estimators (not stopping rules for accuracy), sequential fixed-accuracy estimation of autoregressive and threshold AR/ARCH parameters, and the general theory of fixed-accuracy intervals for a positive parameter (Banerjee and Mukhopadhyay 2016) and of distribution-free minimum risk estimation (Khariton and Mukhopadhyay 2024). This absence is a statement about our searches, not a proof that no such paper exists.

**What is new here** (to our knowledge, and modest): the sequential impossibility statement for fixed-accuracy estimation of a variance (Proposition 1, an easy extension of Bahadur and Savage 1956); the use of a MAD-truncated statistic to make an exactly valid fixed-accuracy procedure possible; the closed-form account of the price of the confidence sequence (`2 log(2/α) / z²`, the Chernoff reading, the Brownian limit with intersection, and the practical-d prediction with the exact `ψ`); the observation that under volatility clustering the classical rule is still valid for the window's average conditional truncated variance while the running intersection of a confidence sequence is not; and the measurement of all procedures on the same heavy-tailed, jumpy and clustered paths and on real one-minute data.

## 5. Implementation

`seqest.py` (numpy and scipy; imports `avmonitor.confseq.psi_coef` and `avmonitor.robust.MAD_NORMAL_SCALE`; nothing in `src` imports it):

* Closed forms: `zq`, `n_star`, `cs_ratio_limit`, `cs_limit_coverage`, `brownian_limit`, `kappa2`, `kappa4`, `opt_plan`, `cs_n_predict`, `betting_oracle_n`, `sigma2_from_theta_tau`, `elasticity_gaussian`, `garch_target_table`.
* Laws with exact truncated moments (`Law.trunc_moments` by 400-point Gauss-Legendre quadrature of `∫ P(r² > y) dy` and `∫ 2y P(r² > y) dy`, closed forms where they exist): `Gaussian`, `StudentT` (unit variance), `Jumps` (normal plus Bernoulli normal jumps), `Garch` (GARCH(1,1), returns the conditional variances too), `Bootstrap` (stationary bootstrap of Politis and Romano 1994 from one data vector; its truncated moments are the exact empirical ones) and `Windows` (real data, one row per window).
* `Setup` (accuracy `d`, `α`, pilot `m`, `c`, relative or absolute, cap `n_max`, block length, bet cap, batch length, minimum risk `A` and `γ`, grid size and span of the betting sequence, intersection on or off, procedures) and `run(law, setup, reps, seed)`, which runs every procedure on the same simulated paths (common random numbers) block by block, drawing new returns only for paths where some procedure is still running. Procedures: `CR-raw`, `CR-trunc`, `CR-BM-trunc`, `Stein-raw`, `Stein-trunc`, `Oracle-raw`, `Oracle-trunc`, `MinRisk-raw`, `MinRisk-trunc`, `CS-wsr`, `CS-tuned`, `CS-opt`, `CS-bet`.
* `summarise`: coverage at N, misses below and above, `E[N]`, its standard deviation and 5%, 50%, 95% quantiles, `n*`, `E[N/n*]`, regret `E[N - n*]`, the share that reached `n_max`, the coverage of the reported interval, the relative bias, and for the mean-variance rules the coverage of the studentised interval `X̄_N ± z √(ν̂_N/N)` (which separates the error of the stopping rule from the error of the normal approximation).

`run_experiment.py`: parts `limit`, `verify`, `map`, `main`, `abs`, `pilot`, `garch`, `minrisk`, `real`, `boot`, `figs` (section 6). Seeds are fixed per cell from a hash of the cell's name; at most two worker processes.

`tests/test_seqest.py`: exact checks of the formulas (normal truncated moments against the package and Monte Carlo, quadrature against closed forms for t3, t5 and jumps, GARCH autocorrelation and long-run variance against simulation, the limits 1.92056, 2.21451 and 0.99340, the Brownian limit, the inversion to `σ²`), brute-force replays of every stopping rule on fixed data (Chow and Robbins, Nadas, Stein, batch means, minimum risk), the empirical Bernstein sequence against `avmonitor.confseq.eb_confseq`, the betting bounds against exact root finding, time-uniform coverage of both sequences on Beta data, and smoke tests of the engine.

## 6. Experimental methodology

All synthetic laws have unit variance (jumps: unit diffusive variance plus jumps). Default: relative accuracy, pilot `m = 30`, `c = 3`, `α ∈ {0.05, 0.10}`, 2000 paths per cell, `n_max = 600,000`.

* **limit [SYNTHETIC, closed form].** For normal (`c` = 2, 3, 4), t5, t3 and jump returns: `n*`, the empirical Bernstein prediction with tuned and optimised bets, and the betting oracle, for `d` from 0.4 to 0.01 and both `α`; and the Brownian limit with intersection for `κ` = 0.8, 1, 1.2, 1.5.
* **verify [SYNTHETIC].** Normal returns, pilot 100, `d` = 0.1, 0.05, 0.025, 1000 paths: `E[N]/n*` of `CS-tuned`, `CS-opt`, `CS-bet` with and without the running intersection, against the per-path prediction of section 2.9 and the two limits.
* **main [SYNTHETIC].** Normal, t5, t3 and jump returns (1% chance per minute of a normal jump with standard deviation 5); `d` = 0.4, 0.2, 0.1, 0.05; both `α`; all ten fixed-accuracy procedures (the raw oracle is skipped for t3, where `n*` is infinite).
* **abs [SYNTHETIC].** The same in absolute accuracy (`σ² = 1`), `d` = 0.2, 0.1, 0.05, `α = 0.05`.
* **pilot [SYNTHETIC].** Pilot sizes 10, 30, 100 for normal, t3 and jump returns at `d` = 0.2 and 0.1.
* **garch [SYNTHETIC].** GARCH(1,1) with `(a, b)` = (0.10, 0.85) and (0.05, 0.94), 1000 paths, `n_max = 200,000`, `d` = 0.2 and 0.1; batch means with `B = 200` and `B = 1000`.
* **minrisk [SYNTHETIC].** Robbins' rule on raw and truncated squared returns, normal, t5 and jump returns, `n0` = 50, 200, 800, pilot 10, 4000 paths.
* **map [SYNTHETIC].** Bias of the normal-shape inversion `θ_τ → σ²` for every law and `c` (closed form), and the coverage of `σ²` by the mapped `CS-bet` and `CR-trunc` intervals at `d = 0.2`.
* **real [REAL].** Binance spot BTCUSDT one-minute bars, 1 July to 31 December 2025 (history from 1 June 2025 for the profile), log returns divided by `intraday_scale` (20 earlier days, smoothing over 15 minutes, past only). For every day and every hour, a window starts at the first valid minute of that hour and holds the next `30 + 1440` valid returns. Each procedure runs on each window (`d` = 0.5, 0.4, 0.3, 0.2; both `α`; the two confidence sequences also without the running intersection, see 2.10); we record how many returns and how many clock minutes it needs, by start hour, the share that stops within 60 minutes and the share that never stops within the window, and, as a yardstick only, an `n*` proxy `z² CV² / d²` from the window's own truncated squares. If the pilot's median absolute return is zero (price did not move in half the pilot minutes), the scale falls back to the median of the nonzero absolute returns, as `VolatilityMonitor` does. **No coverage is claimed on real data**: the true level is unknown and volatility moves within hours.
* **boot [REAL, resampled].** The 30 days with the smallest dispersion of hourly robust scale. For each, a stationary bootstrap of that day's deseasonalised returns (mean block length 1, that is i.i.d. resampling, and 30), 200 resampled paths, `d = 0.2`, `α = 0.05`; coverage of the day's own truncated variance (known exactly for the bootstrap law). This is a check of the procedures under the real shape of returns, not a claim about real coverage.

## 7. Folder structure and how to run

```
alternatives/sequential_estimation/
├── README.md            this file
├── RESULTS.md           all numbers, labelled [SYNTHETIC] or [REAL], negative results, ranking
├── seqest.py            laws, closed forms, stopping rules, confidence sequences, engine
├── run_experiment.py    parts limit, verify, map, main, abs, pilot, garch, minrisk, real, boot, figs
├── tests/
│   └── test_seqest.py   fast synthetic tests (about 10 seconds)
└── results/
    ├── a_closed_form.csv             predicted E[N]/n* of the sequences, every law, d, alpha
    ├── a_brownian_limit.csv          d -> 0 limit with intersection, by bet multiplier
    ├── b_main.csv                    relative accuracy, four laws, all procedures
    ├── c_absolute.csv                absolute accuracy
    ├── d_pilot.csv                   pilot size 10, 30, 100
    ├── d_pilot_core_grid_only.csv    the same before the grid extension of CS-bet (failed attempt, kept for the record)
    ├── e_garch.csv                   volatility clustering, batch means
    ├── f_minrisk.csv                 Robbins' minimum risk rule
    ├── g_mapping.csv                 truncated variance back to sigma^2
    ├── j_verify_limits.csv           simulated E[N]/n* against the closed forms
    ├── h_real_hourly.csv             BTCUSDT, by start hour
    ├── h_real_summary.csv            BTCUSDT, overall
    ├── i_bootstrap_days.csv          bootstrap of quiet days, per day
    ├── i_bootstrap_summary.csv       bootstrap of quiet days, summary
    ├── _garch_table_*.csv            stationary truncated variance and long-run variance on a tau grid (cache)
    ├── fig1_price_of_anytime_validity.png
    ├── fig2_coverage_by_law.png
    └── fig3_real_minutes_by_hour.png
```

From the repository root, with the package installed (`pip install -e .`):

```
python -m pytest -q alternatives/sequential_estimation/tests
python alternatives/sequential_estimation/run_experiment.py limit verify map main abs pilot garch minrisk
python alternatives/sequential_estimation/run_experiment.py real boot figs
```

Any subset of parts can be named. On two worker processes of a shared machine the synthetic parts took about 15 minutes in all (most of it in `main`, `abs` and `pilot`, where raw Chow and Robbins under t3 needs `d^{-3}` returns, and in `garch`). The real-data parts need the Binance files under `data/raw` (see `data/README.md`) and take about two minutes.

## References

All entries were checked by web search; Project Euclid and Crossref pages could not be opened from this machine, so page ranges come from catalogue and citing pages.

* Anscombe, F. J. (1952). Large-sample theory of sequential estimation. *Proceedings of the Cambridge Philosophical Society*, 48, 600 to 607.
* Aras, G. and Woodroofe, M. (1993). Asymptotic expansions for the moments of a randomly stopped average. *Annals of Statistics*, 21(1), 503 to 519. doi:10.1214/aos/1176349039.
* Bahadur, R. R. and Savage, L. J. (1956). The nonexistence of certain statistical procedures in nonparametric problems. *Annals of Mathematical Statistics*, 27(4), 1115 to 1122.
* Banerjee, S. and Mukhopadhyay, N. (2016). A general sequential fixed-accuracy confidence interval estimation methodology for a positive parameter: illustrations using health and safety data. *Annals of the Institute of Statistical Mathematics*, 68(3), 541 to 570. doi:10.1007/s10463-015-0504-2.
* Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482.
* Barndorff-Nielsen, O. E. and Shephard, N. (2002). Econometric analysis of realized volatility and its use in estimating stochastic volatility models. *Journal of the Royal Statistical Society, Series B*, 64(2), 253 to 280. doi:10.1111/1467-9868.00336.
* Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. *Journal of Econometrics*, 31(3), 307 to 327. doi:10.1016/0304-4076(86)90063-1.
* Bollerslev, T. (1988). On the correlation structure for the generalized autoregressive conditional heteroskedastic process. *Journal of Time Series Analysis*, 9(2), 121 to 131. doi:10.1111/j.1467-9892.1988.tb00459.x.
* Chow, Y. S. and Robbins, H. (1965). On the asymptotic theory of fixed-width sequential confidence intervals for the mean. *Annals of Mathematical Statistics*, 36, 457 to 462. doi:10.1214/aoms/1177700156.
* Chow, Y. S. and Yu, K. F. (1981). The performance of a sequential procedure for the estimation of the mean. *Annals of Statistics*, 9(1), 184 to 189. doi:10.1214/aos/1176345345.
* Ghosh, M. and Mukhopadhyay, N. (1979). Sequential point estimation of the mean when the distribution is unspecified. *Communications in Statistics, Theory and Methods (Series A)*, 8(7), 637 to 652.
* Ghosh, M., Mukhopadhyay, N. and Sen, P. K. (1997). *Sequential Estimation*. Wiley, New York.
* Giné, E., Götze, F. and Mason, D. M. (1997). When is the Student t-statistic asymptotically standard normal? *Annals of Probability*, 25(3), 1514 to 1531. doi:10.1214/aop/1024404523.
* Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2021). Time-uniform, nonparametric, nonasymptotic confidence sequences. *Annals of Statistics*, 49(2), 1055 to 1080. doi:10.1214/20-AOS1991.
* Khariton, Y. and Mukhopadhyay, N. (2024). Distribution-free sequential minimum risk point estimation of an unknown mean (exact title not verified). *Communications in Statistics, Theory and Methods*, 53(7), 2617 to 2644. doi:10.1080/03610926.2022.2147796.
* Konev, V. and Nazarenko, B. (2020). Sequential fixed accuracy estimation for nonstationary autoregressive processes. *Annals of the Institute of Statistical Mathematics*, 72(1), 235 to 264. doi:10.1007/s10463-018-0689-2.
* Lai, T. L. and Siegmund, D. (1983). Fixed accuracy estimation of an autoregressive parameter. *Annals of Statistics*, 11(2), 478 to 485.
* Li, Y., Nolte, I. and Nolte, S. (2019). Renewal based volatility estimation. SSRN working paper (no journal version verified).
* Logan, B. F., Mallows, C. L., Rice, S. O. and Shepp, L. A. (1973). Limit distributions of self-normalized sums. *Annals of Probability*, 1(5), 788 to 809. doi:10.1214/aop/1176996846.
* Mukhopadhyay, N. and de Silva, B. M. (2009). *Sequential Methods and Their Applications*. Chapman and Hall/CRC, Boca Raton.
* Nadas, A. (1967). *An extension of a theorem of Chow and Robbins on sequential confidence intervals for the mean*. PhD dissertation, Columbia University. (A journal version is cited in the literature as Annals of Mathematical Statistics, 1969; we could not verify its volume and pages.)
* Politis, D. N. and Romano, J. P. (1994). The stationary bootstrap. *Journal of the American Statistical Association*, 89(428), 1303 to 1313. doi:10.1080/01621459.1994.10476870.
* Robbins, H. (1959). Sequential estimation of the mean of a normal population. In U. Grenander (ed.), *Probability and Statistics (the Harald Cramér volume)*, Almqvist and Wiksell, Stockholm, 235 to 245.
* Starr, N. (1966). On the asymptotic efficiency of a sequential procedure for estimating the mean. *Annals of Mathematical Statistics*, 37(5), 1173 to 1185. doi:10.1214/aoms/1177699263.
* Stein, C. (1945). A two-sample test for a linear hypothesis whose power is independent of the variance. *Annals of Mathematical Statistics*, 16(3), 243 to 258.
* Ville, J. (1939). *Étude critique de la notion de collectif*. Thèse, Paris.
* Vorobeychikov, S. E., Burkatovskaya, Y. B. and Sergeeva, E. E. (2016). TAR(p)/ARCH(1) process with an arbitrary threshold: guaranteed parameter estimation and change-point detection. *IAENG International Journal of Applied Mathematics*, 46(3), 353 to 366.
* Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *Journal of the Royal Statistical Society, Series B*, 86(1), 1 to 27. doi:10.1093/jrsssb/qkad009.
* Waudby-Smith, I., Arbour, D., Sinha, R., Kennedy, E. H. and Ramdas, A. (2024). Time-uniform central limit theory and asymptotic confidence sequences. *Annals of Statistics*, 52(6), 2613 to 2640. doi:10.1214/24-AOS2408.
* Woodroofe, M. (1977). Second order approximations for sequential point and interval estimation. *Annals of Statistics*, 5(5), 984 to 995. doi:10.1214/aos/1176343953.
