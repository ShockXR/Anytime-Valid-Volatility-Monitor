# Robust local plug-in estimator of integrated quarticity (alternative A5)

A corrected and more efficient version of the robust quarticity estimator of Bapat (2026). It keeps the paper's two ideas, truncation of returns at a multiple of a robust scale and local smoothing, but changes three things: the scale is estimated **locally**, the local estimates are **aggregated by squaring a local variance** (the plug-in principle of Jacod and Rosenbaum 2013) instead of averaging local fourth powers, and every truncation and squaring bias is removed with an **exact finite-sample correction**. A growing cut-off `sqrt(2 log n)` makes the estimator consistent while jumps are still removed. The same machinery gives an estimate of `int sigma^8` and therefore a confidence interval, and it extends to any power of volatility, which the paper names as future work.

A second finding came from the real data. Large one-minute crypto returns arrive in short bursts, not as isolated jumps, so any estimator that truncates against a scale measured over hours (Bapat's included) removes real volatility along with the jumps. A neighbour-guarded variant (`PI-NN`, section 2.9) keeps the bursts and still removes isolated jumps; it is the best robust estimator on a burst model calibrated to the data, and the first proposal (`PI-LU`) remains the best one when volatility is smooth.

Results are in `RESULTS.md`. Every number there comes from `run_experiment.py` and is labelled [SYNTHETIC] or [REAL]. Iteration notes are in `docs/iterations/ALT3.md`.

## 1. Research problem

Integrated quarticity, `IQ = int_0^T sigma_t^4 dt`, is the quantity that sets the width of every confidence interval for daily variance: realized variance has asymptotic variance `2 T IQ / n` (Barndorff-Nielsen and Shephard 2002). It must be estimated from the same intraday returns, and it is much harder to estimate than variance because fourth powers blow up any jump.

Bapat (2026) clips returns at `tau = c x MAD` (MAD: median absolute deviation, a robust spread measure) and averages kernel-weighted local fourth powers. The repository found four problems (`docs/estimator_issues.md`):

* (a) the printed kernel weights are not normalised, so the estimate grows with the bandwidth;
* (b) with `tau = c x MAD` the truncation bias does not vanish as sampling gets finer: only 70.5% (raw MAD) or 96.4% (scaled MAD) of the fourth moment survives, for every n;
* (c) the averaged estimator equals truncated realized quarticity up to edge effects, so the kernel adds nothing for the integrated quantity;
* (d) a clipped jump still contributes `tau^4`.

The paper itself lists data-driven bandwidths and higher powers of volatility as future work. The research question here is:

> Can Bapat's robust-scale truncation be turned into an IQ estimator that is (i) consistent as n grows, (ii) close to the efficient asymptotic variance, (iii) robust to jumps, (iv) equipped with a confidence interval, and (v) better behaved on real one-minute crypto data than the MAD-truncated estimators?

## 2. Mathematical formulation

### 2.1 Set-up

Log price `dp_t = mu_t dt + sigma_t dW_t + dJ_t` on `[0, T]` (one day, T = 1), J a finite-activity jump process. Returns `r_i = p_{i delta} - p_{(i-1) delta}`, i = 1..n, `delta = T / n`. Targets: `IQ = int sigma_t^4 dt` and, for intervals, `IQ8 = int sigma_t^8 dt`. For a standard normal Z write `phi`, `Phi` for its density and distribution function and

`T_q(a) = E[Z^(2q) 1{|Z| <= a}]`, with `T_0 = 2 Phi(a) - 1` and `T_q = (2q - 1) T_{q-1} - 2 a^(2q-1) phi(a)` (integration by parts).

### 2.2 Why aggregation matters

Averaging local fourth powers, `mean_t sum_k w_k r_{t-k}^4 / (3 delta^2)`, is realized quarticity (RQ) whatever the weights (`docs/math.md`, section 4), because every interior return gets total weight one. Its asymptotic variance (AVAR, the limit of `n Var(estimate)` in units of `int sigma^8`) is `Var(Z^4) / 9 = (105 - 9) / 9 = 32/3`.

The plug-in principle estimates the local **variance** first and squares it: `IQ_hat = sum_j (duration of block j) x sigma_hat_j^4` with `sigma_hat_j^2` from the k returns of block j. Squared returns are far less noisy than fourth powers (`Var(Z^2) = 2` against `Var(Z^4) = 96`), and the efficient bound for estimating `int g(c_t) dt` with `c = sigma^2` is `int 2 g'(c_t)^2 c_t^2 dt`, which is `8 int sigma^8` for `g(c) = c^2` (Jacod and Rosenbaum 2013, who also show that plug-in estimators with a bias correction reach it).

### 2.3 The block plug-in with an exact local correction

Cut the day into blocks of about k returns (default `k = ceil(sqrt(n))`: 20, 38 and 153 for n = 390, 1440 and 23,400). In block j let `X_i = g(r_i)^2` be the truncated squared returns (clip or drop, section 2.5), `m_j` the number used, `S_j = sum X_i`, and

`sigma_hat_j^2 = S_j / (m_j delta)`, `sigma_hat_j^4 = (sigma_hat_j^2)^2 / e_2(m_j, a)`, `sigma_hat_j^8 = (sigma_hat_j^2)^4 / e_4(m_j, a)`,

where `e_2 = E[(S/m)^2]` and `e_4 = E[(S/m)^4]` for i.i.d. standard normal increments truncated at `a` standard deviations. They follow from the cumulants of `X`: with raw moments `mu_1..mu_4` of X, S has cumulants `m kappa_q(X)`, and `E[S^2] = K_2 + K_1^2`, `E[S^4] = K_4 + 4 K_3 K_1 + 3 K_2^2 + 6 K_2 K_1^2 + K_1^4`. Without truncation `e_2 = (m + 2) / m` (the classical squared chi-square correction) and `e_4 = (m + 2)(m + 4)(m + 6) / m^3`.

**Proposition 1 (exact unbiasedness, proved).** If volatility is constant within the block, increments are Gaussian and the cut-off `a` is fixed in units of the true local standard deviation, then `E[sigma_hat_j^4] = sigma^4` and `E[sigma_hat_j^8] = sigma^8` for every m, and `Var(sigma_hat_j^4) = sigma^8 (e_4 / e_2^2 - 1)`. Without truncation `e_4 / e_2^2 - 1 = 8 (m + 3) / (m (m + 2))`.
*Proof.* `sigma_hat_j^2 / sigma^2 = S_std / m` where `S_std` is the sum of m i.i.d. copies of X. Take expectations of its square and fourth power. For the variance, `(k + 4)(k + 6) / (k (k + 2)) - 1 = 8 (k + 3) / (k (k + 2))`. For the drop rule, condition on which returns are kept: given the kept set, the kept standardised returns are i.i.d. truncated normals, so the same argument applies with `X = Z^2` given `|Z| <= a`. QED.

**Proposition 2 (efficiency, proved for constant volatility).** With blocks of k returns, `n Var(IQ_hat) = (n / k) (k delta)^2 sigma^8 8 (k + 3) / (k (k + 2)) n = 8 (k + 3) / (k + 2) sigma^8` (T = 1). So:

| estimator | AVAR / `int sigma^8` | efficiency 8 / AVAR |
|---|---|---|
| block plug-in, k = 1 (this is RQ) | 32/3 = 10.67 | 0.75 |
| block plug-in, k = 10 | 8.67 | 0.92 |
| block plug-in, k = 20 | 8.36 | 0.96 |
| block plug-in, k -> infinity | 8 | 1 |
| tripower quarticity, p = 4/3 | 13.65 (closed form below) | 0.59 |
| quad-power quarticity | 14.92 | 0.54 |

For tripower quarticity the summands are 2-dependent, so with `a = mu_{8/3}`, `b = mu_{4/3}` (`mu_p = E|Z|^p`) the AVAR is `(a^3 + 2 a^2 b^2 + 2 a b^4 - 5 b^6) / b^6`; quad-power is analogous with 3-dependence. MedRQ and MinRQ are only simulated (`results/q1_avar.csv`). With stochastic volatility the same calculation holds block by block as long as `k_n sqrt(delta) -> 0`, which keeps the bias from volatility moving inside a block (`about -xi^2 k delta / 6` in relative terms for a log-volatility diffusion coefficient xi, from `E[within-block variance of sigma^2]`) below the standard error. The default `k = sqrt(n)` sits on that boundary: the bias is then a fixed fraction `xi^2 / (6 sqrt 8)` of the standard error (6% for xi = 1). Section 2.8 and the block-length experiment check that this is harmless; a formally centred limit theorem needs `k_n = o(sqrt(n))`, for example `n^0.4`.

### 2.4 Local scale

The threshold needs the local standard deviation of a return, `s_i = sigma_{t_i} sqrt(delta)`. Bapat uses one MAD for the whole day, which is wrong whenever volatility moves within the day (stochastic volatility, the intraday U-shape): in quiet hours almost nothing is clipped and in busy hours far too much. Here the day is cut into scale blocks of `ceil(sqrt(n))` returns and each block's scale uses the block and its two neighbours (about `3 sqrt(n)` returns, so the window grows but shrinks in time, as consistency needs). Two steps: the scaled MAD `s0 = 1.4826 median|r|`, then the truncated root mean square `s = sqrt(mean(r^2 1{|r| <= 4 s0}) / T_1(4))`. The second step matters: with the MAD alone (efficiency 37%) the realised cut-off `a s_hat / s` is noisy, more returns are truncated than the correction assumes, and the estimate was 0.4% to 1.3% too low at n = 390 under constant volatility; with the two-step scale the clipping variants and the universal drop rule are within 0.6% (`results/q1_scale_ablation.csv`, iteration 1 in `docs/iterations/ALT3.md`). Dropping at a fixed a = 3 then overshoots by 2%, because the scale window contains the block itself (a block with large returns gets a larger threshold and fewer drops than the correction assumes); that is one reason to pair dropping with the growing cut-off.

### 2.5 Truncation rules: which ones are consistent

Per-return threshold `tau_i = a_n s_hat_i`. Two maps: **clip** `g(r) = max(-tau, min(r, tau))` (Bapat's eq. 11) and **drop** `g(r) = r 1{|r| <= tau}` with the block re-counted (Mancini 2009 style).

* **Fixed cut-off, one daily scale (Bapat's rule, `a = c = 3` with the scaled MAD).** The realised cut-off at time t is `a_t = c x MAD_day / s_t`, which varies with the volatility path. Any correction uses one `kappa(a)`, so the limit is `int sigma^4 kappa_4(a_t) dt / kappa_4(c)`, not IQ. **Inconsistent** under stochastic volatility or an intraday pattern, with or without the kappa correction (the correction helps only if volatility is constant).
* **Fixed cut-off, local scale, exact correction (`PI-L3`).** Under the model, increments are locally Gaussian, so Proposition 1 applies block by block: **consistent**, at rate `sqrt(n)`. The price: the correction relies on the *shape* of the local return distribution. If one-minute returns are locally fat-tailed (they are, see RESULTS Q5), a fixed `a = 3` keeps biting into genuine variation and the Gaussian correction does not undo it.
* **Growing cut-off `a_n = sqrt(2 log n)` with a local scale (`PI-LU`, drop; `PI-LUc`, clip).** The expected number of diffusive returns beyond the cut-off in a day is `n P(|Z| > a_n) ~ 1 / sqrt(pi log n) -> 0` (0.22, 0.05 and 0.003 for n = 390, 1440, 23,400), and the relative truncation bias in the local variance is `2 (a_n phi(a_n) + 1 - Phi(a_n)) = O(sqrt(log n) / n) = o(n^(-1/2))`, so it is negligible even *without* a correction and whatever the local shape. Meanwhile `tau_n ~ s sqrt(2 log n) = sigma sqrt(2 delta log(1/delta)) -> 0`, so every jump of fixed size is eventually beyond the threshold. This is the slowest growth for which the expected number of wrongly truncated returns vanishes, and it is the order of the optimal threshold derived for finite-activity jump models (Figueroa-Lopez and Nisen 2013, as reported in later work by the same group; we could not read the full text). Any rule `a_n = c n^gamma` with `0 < gamma < 1/2` (Mancini's `Delta^varpi` with `varpi = 1/2 - gamma`) is also consistent; at practical n, Mancini's usual `varpi = 0.49, c = 4` gives `a_n` of 4.2 to 4.4, close to `sqrt(2 log n)`.
* **Drop or clip?** A clipped jump still contributes `tau^2` to its block. **Proposition 3 (proved, constant sigma, oracle cut-off):** the expected extra IQ from one clipped jump is, in units of `sigma^4 delta`,
  `(a^4 + 2 a^2 (k - 1) mu_1 - mu_2 - 2 (k - 1) mu_1^2) / (mu_2 + (k - 1) mu_1^2)`,
  with `mu_1, mu_2` the first two moments of the clipped squared normal. That is about `2 (a^2 - 1)` for the plug-in (17.5 at a = 3, k = 20) and `(a^4 - mu_2) / mu_2` for fourth powers (27.0 at a = 3). Per jump this is a relative bias of `O(1/n)`, so clipping is still consistent, but at n = 390 four jumps a day cost 18% (plug-in) to 28% (fourth powers). Dropping with a re-count leaves no residual (Proposition 1 conditional on the kept set). Hence `PI-LU` drops.
* **Undetected jumps.** **Proposition 4 (proved):** a jump of `b` local standard deviations that stays below the threshold adds `2 b^2 + b^4 / (k + 2)` (units `sigma^4 delta`) to the plug-in, against `2 b^2 + b^4 / 3` for fourth powers (k = 1). The plug-in damps the `b^4` term by `3 / (k + 2)`: for b = 3 and k = 20 it is 21.7 against 45.

**What the estimated scale does to the variance.** With the true sigma in the threshold (an "oracle" cut-off at a = 3) the clipped block estimator has an AVAR *below* the bound 8 (closed form `k (e_4 / e_2^2 - 1)`; see `results/q1_avar.csv`). That is not a contradiction: the oracle uses the unknown sigma. With an estimated scale the truncation indicator is correlated with the block's own returns (a large block gets a larger threshold), which adds variance, and the simulated AVAR of the feasible estimator is above 8 again (same file). The confidence interval therefore uses the untruncated Gaussian factor `8 (m + 3) / (m (m + 2))`, which is the right order for the feasible estimator.

### 2.6 Bapat's kernel, aggregated correctly (`KPI-L3`)

The kernel idea survives if it is used for the local variance: `V_t = sum_k w_k g(r_{t-k})^2 / delta` with the paper's one-sided Gaussian kernel (normalised, h = 30), then `IQ_hat = T mean_t V_t^2 / (mu_1^2 + (mu_2 - mu_1^2) sum_k w_k^2)`. The denominator is `E[V_t^2] / sigma^4` for clipped Gaussian increments (proved as in Proposition 1), and `1 / sum w_k^2` is the kernel's effective sample size. Overlapping windows behave like blocks of that effective size (simulated AVAR in `results/q1_avar.csv`).

### 2.7 Confidence interval and higher powers

`Var(IQ_hat)` is estimated by `sum_j (duration_j)^2 sigma_hat_j^8 8 (m_j + 3) / (m_j (m_j + 2))`, using the unbiased block estimate of `sigma^8` from Proposition 1. The 95% interval is `IQ_hat +/- 1.96 se` (plain) or `IQ_hat exp(+/- 1.96 se / IQ_hat)` (log scale, which respects positivity and the right skew of a variance-type statistic). The same construction estimates `int sigma^(2p) dt` for any p by dividing `(sigma_hat_j^2)^p` by `E[(S/m)^p]`, which answers the paper's "sixth or eighth power" future-work item; `IQ8_hat` is reported in RESULTS.

### 2.8 Bandwidth and cross-validation

Eq. (17) of the paper compares `r_t^4` with a leave-one-out kernel quarticity `Q_{t,-t}(h)` built from raw weights. **Proposition 5 (proved, constant volatility, i.i.d. clipped Gaussian returns).** With `q_1 = E[g(Z)^4]`, `q_2 = E[g(Z)^8]` and leave-one-out raw weights `v_k = K(-k/h)`, k = 1..h, in units of `sigma^8 delta^4`:

* printed: `E[CV(h)] = (q_2 - q_1^2)(1 + sum v^2 / 9) + q_1^2 (1 - sum v / 3)^2`;
* corrected (target `r_t^4`, weights renormalised to sum to one): `E[CV(h)] = (q_2 - q_1^2)(1 + sum w^2)`.

The printed criterion has a squared bias `q_1^2 (1 - sum v / 3)^2` and a variance term `sum v^2 / 9`, both growing linearly in h because the raw weights do not shrink. Its minimum is at h = 5, below the paper's grid {10, ..., 50}, and on that grid it always picks h = 10: the choice is a scale artefact. The corrected criterion decreases in h when volatility is constant, as it should. More generally, if `r_t` is independent of the past given the volatility path, `E[(r_t^4 - P_t)^2] = E[Var(r_t^4 | sigma)] + E[(E[r_t^4 | sigma] - P_t)^2]` for any predictor P built from the past, so the corrected CV is the mean squared error of the local estimate plus a constant. The constant (the noise of `r_t^4`, a 96 `sigma^8 delta^4` term) is large, so a squared-return target (`r_t^2` against the leave-one-out local variance) is a less noisy version for the plug-in estimators.

For the **integrated** quantity the relevant tuning parameter is the block length k, and the trade-off is different from the spot one: variance falls like `8 (1 + 1 / (k + 2))` while the within-block volatility bias grows like `xi^2 k delta / 6`. Minimising `8 delta / k + (xi^2 k delta / 6)^2` gives `k* = (144 n / xi^4)^(1/3)` (38 for n = 390 and xi = 1), and the curve is flat around it. The experiment measures that flatness directly.

### 2.9 Volatility bursts and the neighbour guard (`PI-NN`)

The real-data check (RESULTS, Q5 and Q6) showed that large one-minute crypto returns are not isolated: they come in bursts lasting a few minutes. Under the model these bursts are part of `sigma_t` and therefore of IQ. A threshold built from a scale over about two hours treats a burst like a run of jumps and removes it, and a block of 38 minutes averages it away (Jensen: the square of a block mean of `sigma^2` is below the block mean of `sigma^4`). Two changes address this while keeping the theory of sections 2.3 to 2.5:

* **Neighbour guard.** The scale for return i becomes `max(s_local_i, s_nbr_i)`, where `s_nbr_i = 1.4826 median(|r_{i-3}|, ..., |r_{i-1}|, |r_{i+1}|, ..., |r_{i+3}|)` is the nearest-neighbour scale of Andersen, Dobrev and Schaumburg (2012), without the return itself. Inside a burst the neighbours are large, so the return is kept; an isolated jump has ordinary neighbours and is removed. Because the guard only raises the threshold, it never truncates a diffusive return that the local rule would keep, so the truncation bias stays `O(sqrt(log n) / n)`. A fixed-size jump is still removed once `tau_n -> 0`, unless another jump falls among its six neighbours, an event whose probability vanishes for finite-activity jumps. Unlike MedRQ, one neighbouring jump does not disarm the guard, because the median of six neighbours moves little.
* **Shorter blocks**, `k_n = round(n^(1/4))` (4, 6, 12). The efficiency `(k + 2) / (k + 3)` is lower (0.86 to 0.93), but block durations shrink like `n^(-3/4)`, so `k_n sqrt(delta) -> 0` (the formal centring condition that `sqrt(n)` only meets at the boundary) and bursts of a few minutes are resolved at fine sampling.

The finite-sample correction is the one of `PI-LU` (cut-off `a_n`); since the guard truncates less, it slightly over-corrects, by at most the drop correction at `a_n` (about 1.5% of `sigma^4` at n = 390, less at larger n).

### 2.10 What is proved, checked and assumed

* **Proved (elementary):** Propositions 1 to 5, the AVAR constants of RQ, TPQ, QPQ and the block plug-in, the truncation constants.
* **Checked numerically:** all of the above against simulation (`results/q1_*.csv`, tests); MedRQ and MinRQ AVARs; the feasible estimators with estimated scales; behaviour under stochastic volatility, intraday pattern and jumps; coverage of the intervals.
* **Not proved:** the neighbour guard's effect on the variance and the over-correction of `PI-NN` are only checked by simulation.
* **Assumed, not proved here:** the full stable central limit theorem with stochastic volatility and jumps for the feasible estimator. It should follow from Jacod and Rosenbaum (2013) for the drop rule with `a_n = sqrt(2 log n)` (their threshold `alpha delta^varpi` is replaced by a data-driven one of the same order up to logs), but that extension has not been written out.

## 3. Why it matters

* Every confidence interval and jump test for daily variance needs IQ. A robust IQ that is too low by a factor 3 to 7 on real data (E5 of the main report) gives intervals that are far too narrow and jump tests that reject far too often.
* It answers the paper's own future-work items (data-driven tuning, higher powers) and fixes the four issues of `docs/estimator_issues.md` while keeping the paper's core idea, robust-scale truncation.
* It turns Bapat's local kernel quantity into a principled building block: the same local variance `V_t` drives the online monitor of the main package, and squaring it with the right correction gives spot quarticity.

## 4. Comparison with existing approaches

| Approach | Robust to jumps | Consistent | AVAR / `int sigma^8` | Notes |
|---|---|---|---|---|
| RQ (Barndorff-Nielsen and Shephard 2002) | no | yes, no jumps | 10.67 | baseline |
| Tripower / quad-power quarticity | yes | yes | 13.65 / 14.92 | needs no tuning |
| MinRQ, MedRQ (Andersen, Dobrev and Schaumburg 2012) | yes | yes | simulated | nearest-neighbour truncation |
| Truncated RQ (Mancini 2009) | yes | yes | 10.67 | drops returns above `c Delta^varpi` |
| Bapat (2026), normalised | partly (clipped jumps leave `tau^4`) | no (fixed c, daily MAD) | about 8.6 to 9.3 (biased) | this repository's E2 |
| Jacod and Rosenbaum (2013) | yes, with truncation | yes | 8 (efficient) | asymptotic bias correction, Mancini threshold |
| Mykland and Zhang (2009) | no | yes | block estimators, local Gaussian approximation | the exact chi-square block correction comes from this "locally constant" view |
| Figueroa-Lopez, Pang and Wu (2025, Econometric Theory) | not in the abstract | yes | kernel spot estimators plugged into functionals | closest to the kernel plug-in `KPI` |
| Guo (2026, arXiv 2604.27665) | not stated in the abstract | yes | 56 `int c^4` for twice IQ | squared returns times local spot volatility |
| **This prototype, PI-LU** | yes (drop + local scale) | yes (`a_n = sqrt(2 log n)`) | 8 (k + 3) / (k + 2) | exact finite-sample corrections, robust local scale, CI, `int sigma^8` |
| **This prototype, PI-NN** | yes for isolated jumps, keeps volatility bursts | yes | 8 (k + 3) / (k + 2), k = n^(1/4) | PI-LU plus a nearest-neighbour guard and short blocks |

Literature details were checked by web search (search-result snippets; the arXiv and Project Euclid full texts could not be fetched from this environment). See RESULTS.md, "Novelty", for what is and is not new.

## 5. Implementation

`riq.py` (imports only `avmonitor.robust` helpers; nothing in `src` imports it):

* `trunc_even_moments`, `square_moments`, `mean_power_moments`, `gaussian_block_relvar`: the moment algebra of Proposition 1.
* `local_scale` (two-step MAD then truncated RMS over a block and its neighbours), `neighbour_scale` (6 nearest neighbours), `nn_block`, `thresholds` (rules `fixed`, `universal`, `none`; scale `local`, `global` or `local+neighbour`).
* `block_plugin` (IQ, `int sigma^8`, standard error, share truncated; clip or drop; NaN-safe), `kernel_plugin` (KPI), `truncated_rq_local` (fourth powers with the same threshold).
* Theory helpers: `avar_block`, `avar_block_exact`, `avar_tpq`, `avar_qpq`, `avar_truncated_rq`, `avar_block_truncated`, `clipped_jump_residual`, `undetected_jump_residual`, `printed_cv_expectation`, `cv_curves`.

Estimator labels used in the results:

| label | definition |
|---|---|
| RQ, TPQ, MedRQ, MinRQ | `avmonitor.robust` (TPQ with p = 4/3) |
| TRQ-M | `avmonitor.robust.threshold_rq`: drop above `4 sqrt(BV) delta^0.49` (daily bipower scale) |
| RK, RK-k | `avmonitor.robust.kernel_iq`, h = 30, c = 3, scaled daily MAD, without / with the kappa_4 correction (Bapat's estimator, normalised) |
| TRQ-U | fourth powers, local scale, `a_n = sqrt(2 log n)`, drop, corrected |
| PI-G3 | block plug-in, Bapat's daily threshold (3 x scaled MAD), clip, corrected for a = 3 |
| PI-L3 | block plug-in, local scale, a = 3, clip, corrected |
| PI-LUc | block plug-in, local scale, `sqrt(2 log n)`, clip, corrected |
| **PI-LU** | block plug-in, local scale, `sqrt(2 log n)`, **drop**, corrected (main proposal) |
| PI-LU-raw | PI-LU without the finite-sample correction |
| KPI-L3 | kernel plug-in (Gaussian kernel h = 30 on the local variance), local scale, a = 3, clip, corrected |
| **PI-NN** | block plug-in, `k = round(n^(1/4))`, scale `max(local, 6-neighbour)`, `sqrt(2 log n)`, drop, corrected (second proposal, iteration 3) |

## 6. Methodology

* **Q1 theory** (`part_theory`): closed forms, and constant-volatility simulations (n = 1440, 8000 days) of every estimator's AVAR; truncation constants; jump residuals (closed form and a 40,000-path check); expected CV curves.
* **Q2 and Q4 simulation grid** (`part_sim`): log-OU stochastic volatility from `avmonitor.simulate.sv_jump_returns` (mean reversion 5 per day, daily volatility 0.02, vol of vol xi = 0.5, or xi = 1 plus a U-shaped intraday pattern of amplitude 0.6), with no jumps or 4 jumps a day of size `N(0, (f x 0.02)^2)`, f = 0.1, 0.3, 1 (about 2, 6 and 20 one-minute standard deviations at n = 390; jumps keep their size in price units as n grows, which is the asymptotic regime of the theory). n = 390, 1440 and 23,400; 500 days per cell; fixed seeds. Relative bias, standard deviation, RMSE, `sqrt(n) x RMSE`, slope of log RMSE on log n, coverage of the 95% intervals and accuracy of `IQ8_hat`. The true IQ is the exact integral of the simulated volatility path; `int sigma^8` is computed from the per-step averaged variance.
* **Q3 bandwidth** (`part_cv`, and the block-length table from `part_sim`): printed and corrected CV curves on the paper's robust (clipped) returns, against the oracle spot MSE computed from the known volatility, plus the IQ error of the kernel estimators as a function of h; IQ error of PI-LU for k from 2 to 320.
* **Q6 volatility bursts** (`part_fastvol`): diagnostics of clustering and tail thickness of locally standardised one-minute returns (autocorrelation of capped `z^2` at lags 1, 5, 10; trimmed kurtosis; share beyond 3; share of exceedances with a neighbour beyond 2), on the real data and on a grid of two-factor volatility models (the S5 slow factor plus a fast log-OU factor with sd 0.5 to 1 and half-life 1 to 10 minutes). The grid point closest to the real diagnostics is then used to score every estimator against the exact IQ at n = 390, 1440, 23,400, with and without jumps, and to compare 1-minute with 5-minute estimates.
* **Q5 real data** (`part_real`): Binance BTCUSDT and DOGEUSDT, 1 July to 31 December 2025, UTC days, 1-minute and 5-minute log returns (missing minutes dropped). Ratios to MedRQ and MinRQ day by day, and diagnostics of local tail thickness.

The defaults of PI-LU (`k = ceil(sqrt n)`, scale window of three blocks, b = 4 in the two-step scale, c = 3, h = 30) were fixed before the grid was run; the one change after a first check (two-step scale instead of MAD alone) was made on constant-volatility Gaussian data. PI-NN was designed *after* seeing Q5 and Q6: its settings (six neighbours, `k = n^(1/4)`) were chosen in an exploratory run on the burst model and on S3 and S4 with other seeds, so its numbers in the grid and on the real data are not out of sample in the strict sense. Both points are recorded in the iteration notes.

## 7. Folder structure

```
alternatives/robust_iq_estimator/
  README.md            this file
  RESULTS.md           every number, labelled [SYNTHETIC] or [REAL], and the ranking paragraph
  riq.py               estimators and theory helpers
  run_experiment.py    Q1 to Q6 (python run_experiment.py [theory|ablation|sim|cv|real|fastvol|figures|all])
  make_figures.py      figures from the CSV files
  tests/test_robust_iq_estimator.py   synthetic tests, about 3 seconds, no data files needed
  results/             CSV outputs and figures
```

Run from the repository root: `python alternatives/robust_iq_estimator/run_experiment.py all` (about 10 to 15 minutes on two cores; the real-data part needs `data/raw`). Tests: `python -m pytest -q alternatives/robust_iq_estimator/tests`.

## References

* Andersen, T. G., Dobrev, D. and Schaumburg, E. (2012). Jump-robust volatility estimation using nearest neighbor truncation. *Journal of Econometrics* 169, 75 to 93.
* Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance* 109, 105482.
* Barndorff-Nielsen, O. E. and Shephard, N. (2002). Econometric analysis of realized volatility and its use in estimating stochastic volatility models. *Journal of the Royal Statistical Society B* 64, 253 to 280.
* Corsi, F., Pirino, D. and Reno, R. (2010). Threshold bipower variation and the impact of jumps on volatility forecasting. *Journal of Econometrics* 159(2), 276 to 288.
* Figueroa-Lopez, J. E. and Nisen, J. (2013). Optimally thresholded realized power variations for Levy jump diffusion models. *Stochastic Processes and their Applications* 123(7), 2648 to 2677.
* Figueroa-Lopez, J. E., Pang, J. and Wu, B. (2025). Estimation of integrated volatility functionals with kernel spot volatility estimators. *Econometric Theory*, doi:10.1017/S0266466625100170 (arXiv 2407.09759).
* Guo, Y. (2026). A note on estimation of quarticity based on spot volatility. arXiv 2604.27665.
* Jacod, J. and Rosenbaum, M. (2013). Quarticity and other functionals of volatility: efficient estimation. *Annals of Statistics* 41(3), 1462 to 1484. doi:10.1214/13-AOS1115.
* Mancini, C. (2009). Non-parametric threshold estimation for models with stochastic diffusion coefficient and jumps. *Scandinavian Journal of Statistics* 36, 270 to 296.
* Mykland, P. A. and Zhang, L. (2009). Inference for continuous semimartingales observed at high frequency. *Econometrica* 77(5), 1403 to 1445. doi:10.3982/ECTA7417.
