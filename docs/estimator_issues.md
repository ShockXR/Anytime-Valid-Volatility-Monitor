# Investigation of the estimator in Bapat (2026)

This note checks two possible problems that an earlier summary raised about the robust kernel quarticity estimator, plus several smaller points found while reproducing the paper. The aim is to understand the estimator, not to find fault with it. Every claim is labelled with the kind of evidence behind it:

* **Algebra**: follows directly from the printed formulas.
* **Simulation**: our own synthetic experiments (scripts in `experiments/`, numbers in `results/`).
* **Published**: numbers printed in the paper.
* **Real data**: Bitstamp BTC/USD 1-minute bars, 1 July to 31 December 2025 (E5), and Binance spot 1-minute bars for the paper's five coins (E6, ALT3).
* **Inference**: a reading of the published numbers that fits the evidence but cannot be confirmed without the author's code.

Plain-language background is in `docs/paper_review.md`. The formal statements are in `docs/math.md`.

## Summary

| # | Concern | Verdict | Size | Easy fix? |
|---|---|---|---|---|
| 1 | Kernel estimator (eqs. 6, 7, 12, 16) is not on the scale of integrated quarticity | **Genuine, in the printed formulas.** The published tables behave closely as the printed formulas predict (except the σ_J = 0.10 robust rows, 4h), so the code very likely matches them. | The estimate is off by a factor that depends on the bandwidth h (about 10.6 for h = 30 on the per-step scale, and 10.6 / n² on the scale of eq. 3). | Yes: use weights that sum to one and divide by the step length squared. |
| 2 | Truncation at c times MAD leaves a bias that does not vanish | **Genuine limitation of Theorem 1 as stated.** Its size depends on whether the MAD is scaled, which the paper does not say. | Keeps 70.5% of the fourth moment with the raw MAD (bias of about minus 30%), or 96.4% with the scaled MAD (about minus 3.6%), for every n. | Yes: divide by a known constant, or let the threshold shrink more slowly than the returns. |
| 3 | The averaged estimator (eq. 16) is the same as truncated realized quarticity | **Genuine (algebra).** The kernel and bandwidth only matter for the local estimate. | Correlation 0.98 to 0.99 with plain truncated RQ in every scenario we ran. | Nothing to fix: it means the kernel's real value is local, which is what an online monitor needs. |
| 4 | Smaller points: weight order in eq. (15), noise model, TPQ and BPQ constants, one claim about Table 7, the real-data RMSE | Mixed: see below. | Mostly minor. | Yes. |
| 5 | One daily MAD when volatility moves within the day, and bursts of large returns on real data (round 2, ALT3) | **Genuine outside the paper's simulation design**, which keeps volatility constant within the day. | Bias still about minus 11% at n = 23,400 with the kappa correction (minus 14% without) under stochastic volatility with an intraday pattern [SYNTHETIC]; on BTC and DOGE one-minute data, estimators with one daily MAD threshold come out at 0.33 to 0.46 times MedRQ [REAL]. | Partly: a local scale restores convergence; bursts need a neighbour guard, still under study. |

None of this undermines the core idea, which is to truncate returns with a robust scale and then kernel-smooth. Points 1 and 2 change the numbers in the simulation tables. Point 3 changes where the method is most useful.

---

## Issue 1. The scale of the kernel estimator

### The concern

Integrated quarticity over a day is `IQ = ∫ σ_t^4 dt`. For one-minute returns with step length Δ = T / n, a normal return has `E[r_i^4] = 3 σ_i^4 Δ²`. Realized quarticity (eq. 3) therefore multiplies the sum of fourth powers by `n / 3`, which puts it on the IQ scale (with T = 1).

The kernel estimator as printed in eq. (6) is

`Q_t = (1/3) Σ_{i=t-h}^{t} K((i - t)/h) r_i^4`

with no normalisation of the weights and no factor of n.

### Algebra

Write `W(h) = Σ_{k=0}^{h} K(k/h)` for the sum of the raw weights. For the Gaussian kernel, `W(h)` is about `0.341 h + 0.32`. So `W(10) = 3.7`, `W(30) = 10.6` and `W(50) = 17.4`. Under locally constant volatility,

`E[Q_t] = W(h) σ_t^4 Δ²`.

Averaging over t (eq. 7) keeps the same factor. So the printed estimator targets `W(h) / n²` times IQ (with T = 1), or `W(h)` times the per-step quarticity `σ^4 Δ²`. Either way the result **grows with the bandwidth**, which a consistent estimator should not do.

### Simulation (E1, E2)

* With constant volatility and n = 390, the printed eq. (7) returns 6.94 × 10⁻⁵ times the true IQ at h = 30, matching `W(30) / 390² = 6.94 × 10⁻⁵` (test `test_printed_kernel_scale_is_weight_sum_over_n_squared`).
* The normalised version (weights summing to one, times `1 / (3 Δ²)`) is unbiased: its mean ratio to the truth is 0.9998.

### Do the published tables come from the printed formulas? (E1, Inference)

The paper does not report the diffusive volatility used in its simulations. If the code follows the printed formulas and measures error against the per-step quarticity `σ^4`, then the kernel RMSE in Table 7 should grow roughly in proportion to `W(h) - 1`. We tested this:

* We fitted one number, the per-step volatility σ, on the 15 rows of Table 7 with σ_J = 0.01, where jumps barely matter. The best fit is **σ = 0.0215** (σ⁴ = 2.14 × 10⁻⁷). With that single number, our simulated RMSEs for the printed formula match the published ones with a typical error of about 3.5% (mean squared log error 0.0012).
* Out of sample, the same σ gives a median simulated / published ratio of 1.11 for the σ_J = 0.05 rows and 1.40 for the σ_J = 0.10 rows, where 1000 paths with large jumps make RMSE noisy. See `results/figures/e01_table7_reproduction.png`.
* Table 3 and Table 4 report an RQ bias of about 1.0 × 10⁻⁴. With σ = 0.0215, `(1/3) Σ r^4` (eq. 3 without the factor n) has expected value 8.3 × 10⁻⁵ from the diffusion plus 3.3 × 10⁻⁵ from the jumps (including the cross term 6σ²σ_J²), giving 1.17 × 10⁻⁴. That is within 8% of Table 4 and 17% of Table 3, the right order of magnitude for RQ without the factor n (with it, the value would be 390 times larger). (Corrected after the independent review: an earlier version left out the cross term and reported 1.07 × 10⁻⁴.)

**The real-data tables give the cleanest check.** If the kernel in Table 6 is eq. (7) over the whole sample and RQ is `(1/3) Σ r^4`, then their ratio is about `W(30) / (n - 30)` regardless of the data. With n = 264,960 one-minute returns (the length of 1 July to 31 December 2025), this is 3.99 × 10⁻⁵. The published ratios Kernel / RQ in Table 6 are:

| Coin | Kernel | RQ | Ratio |
|---|---|---|---|
| BTC | 3.98 × 10⁻¹¹ | 9.92 × 10⁻⁷ | 4.01 × 10⁻⁵ |
| ETH | 2.01 × 10⁻¹⁰ | 5.02 × 10⁻⁶ | 4.00 × 10⁻⁵ |
| XRP | 9.42 × 10⁻⁹ | 2.34 × 10⁻⁴ | 4.03 × 10⁻⁵ |
| SOL | 3.10 × 10⁻¹⁰ | 7.73 × 10⁻⁶ | 4.01 × 10⁻⁵ |
| DOGE | 2.67 × 10⁻⁹ | 6.66 × 10⁻⁵ | 4.01 × 10⁻⁵ |

All five coins agree with the prediction to within 1%. On our Bitstamp BTC data the same two formulas give a ratio of 3.99 × 10⁻⁵ (E5a). This is strong evidence that the code uses unnormalised weights, and that RQ in the tables is computed without the factor n.

### Consequences

1. **The RMSE comparisons mix scales.** RQ as computed in the tables, `(1/3) Σ r^4`, is about n times the per-step target (with the factor n of eq. 3 it would be n² times), the kernel about W(h) times, and the robust kernel about 0.705 W(h) times (see Issue 2). Differences in RMSE therefore mostly reflect different scale factors, not different accuracy.
2. **"Robust beats plain kernel" in Tables 7 and 8 is mostly a scale effect.** The plain kernel overshoots by about W(h). Truncation multiplies by about 0.705, which pulls the overshoot towards the truth. That is why the gap widens with h. Once both estimators are normalised (E1 counterfactual, `results/e01_normalised_counterfactual.csv`), the comparison is fair, and the robust idea holds up well. With the raw MAD, the normalised robust kernel loses to the plain kernel only when jumps are small (σ_J = 0.01: relative RMSE 0.32 against 0.17). It wins at σ_J = 0.05 (medians 0.29 against 0.83, in all 15 cells) and at σ_J = 0.10 (0.28 against 10.6). With the scaled MAD it wins in all 45 cells, including σ_J = 0.01 (0.148 to 0.154 against 0.164 to 0.172).
3. **Cross-validation (eqs. 17, 18) is affected too.** CV compares `r_t^4`, whose mean is `3 σ^4 Δ²`, with `Q_t`, whose mean is `W(h) σ^4 Δ²`. The mismatch changes with h, so the chosen h partly reflects the scale factor. With normalised weights, the target in eq. (17) should be `r_t^4 / 3`, or equivalently `3 Q_t` should be compared with `r_t^4`.

### Verdict

Genuine in the printed formulas, and very likely in the code behind the tables. It is a scaling problem, not a conceptual one. The fix is one line, and `avmonitor.robust.local_quarticity` and `kernel_iq` implement it.

**Question for Prof. Bapat:** does the code normalise the kernel weights, and is the simulation target the per-step quarticity or the integrated one?

---

## Issue 2. Truncation bias that does not vanish

### The concern

Eq. (11) clips returns at `τ = c · MAD(r)` with c = 3. For normal returns, the raw MAD equals `0.6745 σ` (σ here is the per-step standard deviation), so `τ = 2.02 σ`. Clipping at 2.02 standard deviations changes the fourth moment a lot.

### Algebra

For a standard normal Z and cut-off a, the share of the fourth moment that survives clipping is

`κ₄(a) = E[clip(Z, -a, a)^4] / 3`,

which has a closed form (`docs/math.md`, section 3). The key values are:

| MAD convention | Cut-off a (in sd) | Returns clipped | κ₄ (fourth moment kept) | κ₂ (variance kept) |
|---|---|---|---|---|
| raw MAD (e.g. `scipy.stats.median_abs_deviation` default) | 2.02 | 4.3% | 0.705 | 0.925 |
| scaled MAD (× 1.4826, R's `mad()` default) | 3.00 | 0.27% | 0.964 | 0.995 |

**Why the bias does not vanish.** The proof of Theorem 1 says that `τ_n → 0` more slowly than the returns shrink, so that eventually almost no return is clipped. But `MAD(r)` shrinks at exactly the same rate as the returns (both are proportional to `σ √Δ`). So `τ / (σ √Δ)` stays at 2.02 (or 3) for every n, the clipped share stays at 4.3% (or 0.27%), and the estimator converges to `κ₄ × IQ`, not to IQ. The threshold estimators of Mancini (2009) avoid this by using `τ_n = c Δ^ϖ` with `ϖ < 1/2`, so that `τ_n / √Δ → ∞`.

### Simulation (E2a)

The simulated ratio of the MAD-truncated estimator to the true IQ matches κ₄ at n = 390, 1440 and 23,400 (`results/e02a_truncation_factor.csv`). For example, with the raw MAD and c = 3 the ratios are 0.705, 0.705 and 0.705. Dividing by κ₄ removes the bias (test `test_raw_mad_truncation_bias_does_not_vanish_with_n`).

With stochastic volatility the bias gets worse, because one full-sample MAD clips calm periods too little and volatile periods too much. After the κ correction, the raw-MAD estimator's mean ratio falls to 0.97, 0.90 and 0.78 as the volatility of volatility rises from 0.5 to 1.0 to 1.5 (E2b). A local (rolling) MAD avoids this.

### Which convention did the paper use? (Inference)

The paper does not say. Two independent checks both point to the **raw** MAD:

* E1: on the Table 7 rows used for fitting, the raw-MAD simulation matches the published robust RMSE with a median ratio of 0.90, against 1.29 for the scaled MAD.
* E5a (real data): the printed robust kernel on Bitstamp BTC gives 1.28 × 10⁻¹³ with the raw MAD and 3.41 × 10⁻¹³ with the scaled MAD. The published value is 1.54 × 10⁻¹³.

### Verdict

This is a genuine limitation of Theorem 1 as stated. The truncation bias does not vanish for any fixed c. Its size is about 30% with the raw MAD and about 3.6% with the scaled MAD. Both are easy to correct: divide by κ₄, or use a threshold that shrinks more slowly than √Δ. The robustness idea is not affected. A bounded influence function (Theorem 3) and asymptotic unbiasedness simply need different threshold rates, and a finite-sample correction can give both in practice.

**Question for Prof. Bapat:** is the MAD in eq. (11) scaled by 1.4826?

---

## Issue 3. The averaged estimator is truncated realized quarticity

### Algebra

Swap the two sums in eq. (16). Each return `r_i` enters every window that contains it, and an interior return (`h < i ≤ n - h`) collects every weight exactly once. So

`Q_robust = (1 / (3 (n - h))) Σ_i c_i r̃_i^4`, where `c_i = W(h)` for every interior i.

Apart from the first and last h returns, the averaged kernel estimator is a constant times the realized quarticity of the truncated returns. The kernel shape and the bandwidth only change that constant and the treatment of the edges. Test `test_aggregated_kernel_is_weighted_realized_quarticity` checks the identity exactly.

### Simulation (E2d)

The normalised robust kernel IQ and plain truncated RQ have correlation 0.982 to 0.990 and mean ratio 0.996 to 1.001 across all seven scenarios (`results/e02d_kernel_vs_truncated_rq.csv`).

### Consequence

For integrated quarticity, the sensitivity to bandwidth and kernel shown in Tables 7 and 8, Figs. 3 to 7 and Table 9 mostly comes from the constant `W(h)` (Issue 1), not from smoothing. The kernel's real contribution is the **local** estimate `Q_t` itself. That is the quantity an online monitor needs, and it is what this repository builds on.

---

## Smaller points

**4a. Weight order in eq. (15).** Eq. (14) defines `w_i = K(i/h)` for i = 0, ..., h, and eq. (15) uses `w_{s - t + h} r*_s` for s = t - h, ..., t. That gives the newest return (s = t) the weight `w_h = K(1)`, the smallest one, and the oldest return the largest. Eq. (6) does the opposite, and the text says recent returns should count more. This looks like a typo. It does not affect the averaged estimator (Issue 3) but it does affect the local one. `local_quarticity_printed(..., reverse_weights=True)` reproduces eq. (15) literally.

**4b. The noise model (eq. 13).** Adding i.i.d. normal noise to each *return* is the same as adding extra volatility. No estimator can separate the two from the returns alone, and truncation at 3 MAD does not try: in E2c, return-level noise with a noise to signal ratio of 0.5 multiplies every estimator's mean by about 1.55 (that is `(1 + 0.5²)²`), robust or not. For example RQ moves from 0.995 to 1.555 times the truth and the raw-MAD robust kernel from 0.680 to 1.070. The second number looks good only because the two biases happen to cancel. The standard microstructure model adds noise to *prices*, which makes returns negatively autocorrelated. Pre-averaging and two-scale estimators exploit that structure. Truncation alone does not address noise under either model, so "robust to noise" should be read as "robust to occasional large noisy prints".

**4c. Noise levels in Table 8.** Table 8 shows RMSE almost unchanged across σ_α² = 0, 0.0005 and 0.001. With σ = 0.0215, a noise *variance* of 0.001 would be larger than the diffusive variance (4.6 × 10⁻⁴), and our simulation then gives RMSEs about 6 times the published ones. Reading 0.0005 and 0.001 as *standard deviations* reproduces Table 8 (median simulated / published ratio 1.12 for the kernel and 0.92 for the robust kernel, E1). Either way the noise in Case 3 is small next to the diffusive returns, so Case 3 does not stress-test noise robustness.

**4d. TPQ and BPQ.** Eq. (5) allows `p ∈ (0, 1)`, but a product of three `|r|^p` terms is on the quarticity scale only when `3p = 4`, so p = 4/3 is the standard choice. Eq. (4) multiplies `|r_i| |r_{i-1}|³` by `μ₁⁻⁴`, but the constant that makes it unbiased is `1 / (μ₁ μ₃)`, and the factor n is missing. The magnitudes in Table 6 fit lower powers better (Inference): the cross-coin ratios of TPQ match volatility to the power of about 1.5 (p = 0.5), and BPQ is within a factor of 3 of the bipower *variation* (a σ² quantity) on our BTC data. Note that μ₁⁻⁴ is the constant of quad-power quarticity, μ₁⁻⁴ n Σ |r_i||r_{i-1}||r_{i-2}||r_{i-3}|; was eq. (4) perhaps meant to be QPQ? Either way, the very large RMSEs of BPQ and TPQ in Tables 3, 4 and 6 come from scale, not from a lack of robustness. With p = 4/3 and the correct constant, TPQ has a relative RMSE of 0.19 without jumps and 0.47 with moderate jumps (E2c). That is far better than the tables suggest, although large jumps still hurt it more than truncation does.

**4e. Table 7 text.** Section 4.2 says that "in all the cases the robust estimator proves to be better than the standard estimator". Table 7 has three exceptions, all at σ_J = 0.10 and h = 10 (for example 2.6 × 10⁻⁶ against 1.8 × 10⁻⁶ at p = 0.005).

**4f. Real-data RMSE (Table 6).** The true quarticity of real data is unknown, and the paper assumes 2.2 × 10⁻⁶. The robust kernel's estimates are about 10⁻¹² to 10⁻¹³, so a single estimate of about 1.5 × 10⁻¹³ cannot have an RMSE of 2.12 × 10⁻⁶ against an assumed 2.2 × 10⁻⁶, while the plain kernel's 3.98 × 10⁻¹¹ would give 2.20 × 10⁻⁶. How was RMSE computed for Table 6? Several estimators have a smaller published RMSE than the robust kernel for some coins: MedRQ for ETH and SOL, MinRQ for ETH (1.25 × 10⁻⁶) and SOL (8.46 × 10⁻⁷), RQ for BTC, and the plain kernel for SOL (2.19 against 2.20 × 10⁻⁶). With every estimator on a common daily-IQ scale (E5b), the jump-robust estimators MedRQ, MinRQ, TPQ and QPQ agree with each other to within about 15%. The MAD-truncated estimators come out 3 to 7 times lower. This is consistent with Issue 2: real one-minute returns are fat-tailed and seasonal within the day, so clipping at a few MADs removes genuine diffusive variation.

**4g. Rates in Theorem 2.** For the averaged estimator a √n rate is plausible, because it is a truncated RQ (Issue 3). The *local* estimator `Q_t` averages about h returns, so its rate is √h, not √n. Step 4 of the proof says "as τ → ∞", while the assumptions have τ → 0. A complete proof would need the threshold rate from Issue 2.

**4h. A result we could not reproduce.** For σ_J = 0.10, the published robust-kernel RMSEs are 1.3 to 22 times larger than ours (median 5.3, quartiles 1.9 to 7.2; median simulated / published ratio 0.19, E1), while the plain kernel matches. Some unreported detail of the large-jump design seems to differ, for example how MAD is computed or how jump sizes are drawn. We report this rather than tune parameters to hide it.

---

## Issue 5. One daily scale when volatility moves within the day (round 2)

This issue comes from the alternative prototype `alternatives/robust_iq_estimator` (ALT3); details, code and every number are in its `RESULTS.md`.

### The concern

The paper truncates every return of a day at 3 times one MAD computed from the whole day. Its simulation design keeps volatility constant within the day, and there this works: the bias is a known constant (Issue 2) and can be divided out. Real one-minute volatility moves during the day (an intraday pattern, news, bursts). Then one daily threshold is too tight in the busy hours and too loose in the quiet ones.

### Simulation (ALT3, Q2) [SYNTHETIC]

With stochastic volatility and an intraday pattern, the daily-MAD rule's bias does not shrink as the sampling gets finer: about minus 11% at n = 23,400 even after the Issue 2 correction (minus 14% without it), and its RMSE falls like n to the power minus 0.1 instead of minus 0.5. Estimators that truncate against a **local** scale (a few dozen neighbouring returns) converge at the usual square-root-of-n rate (fitted slopes minus 0.48 to minus 0.68). So the robustness idea survives; it needs a local scale when volatility is not constant.

### Real data (ALT3, Q5 and Q6) [REAL]

On Binance BTCUSDT and DOGEUSDT one-minute bars (July to December 2025), the estimators with one daily MAD threshold come out at 0.33 to 0.46 times MedRQ (our normalised kernel alone 0.33 to 0.42), local plug-ins at 0.66 to 0.73, and a version that does not truncate a large return when its neighbours are also large at 0.95 to 0.96. Large one-minute returns come in bursts: 44% to 65% of them have a large neighbour, against 8.9% by chance for independent Gaussian returns. (That apparent jumps are often short bursts of volatility is documented for other markets by Christensen, Oomen and Podolskij 2014 and Bajgrowicz, Scaillet and Treccani 2016.) The true IQ is unknown on real data, so these are ratios between estimators, not errors. A model with a fast volatility factor (half-life about 3 minutes), calibrated to these diagnostics, reproduces the real ratios. In that model, every estimator that truncates against a scale measured over hours is biased down by about 35% to 59% at one minute (31% to 56% with jumps), and even a 5-minute MedRQ is only 0.57 of it. So agreement between two estimators, or with a 5-minute estimate, is not evidence of accuracy on this data.

### Verdict

Outside the paper's design, a daily scale is not enough. This strengthens the case made in Issue 3 for the **local** estimate. The burst reading rests on a fitted model, not a proof, and the neighbour-guarded estimator's settings were chosen after seeing the data, so it still needs an out-of-sample check.

**Question for Prof. Bapat:** in the crypto application, was the MAD computed per day, per coin over the whole sample, or locally?

---

## What this means for the prototype

* Use **normalised** local estimators with an explicit time scale (`robust.local_variance`, `robust.local_quarticity`).
* Use the **scaled** MAD, estimated from past data only. E3 shows that a rolling or exponentially weighted MAD costs nothing measurable in accuracy and keeps the jump protection. (The streaming monitor itself estimates the MAD once per calibration window and re-estimates it after each alarm; a running threshold with a proof is studied in E14.)
* Treat truncation primarily as a way to make the monitored statistic **bounded**. A fall in the variance itself cannot be detected by any nonparametric e-test (`docs/math.md` 5.3, a one-sided form of Bahadur and Savage 1956), so a detector must monitor a bounded functional of the returns. Truncation at c × MAD gives one with a direct volatility meaning: if the shape of the innovations is fixed, the truncated variance falls exactly when σ falls. (Other bounded statistics, such as exceedance indicators or ranks, would also work; `alternatives/two_sample_betting` (ALT5) studies a rank-based one.) Truncation bias then matters less, because detection compares like with like (truncated variance before against truncated variance after).
