# Paper review: Bapat (2026)

**Paper.** Sudeep R. Bapat (2026). "A jump and noise robust kernel based estimator of quarticity for high frequency data." *International Review of Economics and Finance* 109, 105482. doi:10.1016/j.iref.2026.105482. Single author, Shailesh J. Mehta School of Management, IIT Bombay. Received 16 March 2026, accepted 3 June 2026, online 8 June 2026. Open access under CC BY-NC 4.0.

Checked against the PDF in the project files. The title, author, journal, volume, article number, DOI and dates above are printed on the paper's first page. Our literature search could confirm the title, author, journal and date online, but not the volume or DOI (the publisher's page could not be fetched), so those come from the PDF only.

This review follows the 17-point structure of the research brief. It is written for a reader new to the area. For the detailed checks of the estimator, see `docs/estimator_issues.md`.

## Key ideas in plain words

* **Volatility** is how much a price moves. **Variance** is volatility squared.
* **Quarticity** is volatility to the fourth power, added up over a day (`IQ = ∫ σ_t⁴ dt`). It sounds abstract, but it has one very practical job: it tells you **how uncertain your estimate of variance is**. If you measure today's variance from one-minute returns, the error bar on that number is set by quarticity.
* **Jumps** are sudden large moves (news, liquidations). **Microstructure noise** is the small jitter added by how trading works (bid-ask bounce, price ticks). Both distort estimates, and fourth powers make it much worse: one jump raised to the fourth power can swamp a whole day.
* **Truncation** (clipping) caps every return at plus or minus τ before computing anything, so a jump cannot dominate.
* **MAD** (median absolute deviation) is a robust measure of spread: the typical distance from the median. Because it uses medians, a few huge returns barely move it.
* **Kernel smoothing** is a weighted moving average where nearby observations get more weight.

## The 17 points

**1. Motivation.** Quarticity is needed for confidence intervals on realized variance, for jump tests (the paper cites Aït-Sahalia and Jacod 2009, Barndorff-Nielsen and Shephard 2005, and Lee and Mykland 2008), for vol-of-vol in option models, and for risk measures. Existing estimators each handle jumps or noise, but not both, and fourth powers make the problem acute.

**2. Problem statement.** Estimate integrated quarticity from high-frequency returns that contain jumps and microstructure noise.

**3. Context.** The log price is a continuous semimartingale (eq. 1), or a jump diffusion (eq. 8). Returns are modelled as `r_t = σ_t ε_t + J_t` (eq. 10), with Bernoulli(p) jumps of size `N(μ_J, σ_J²)`, and optionally `+ α_t` noise (eq. 13).

**4. Existing methods (Section 2).** Realized quarticity RQ (eq. 3), bipower BPQ (eq. 4), tripower TPQ (eq. 5), threshold variants (Mancini 2009, Corsi, Pirino and Renò 2010), MinRQ and MedRQ (Andersen, Dobrev and Schaumburg 2012; the paper cites the 2008 working paper version), pre-averaging (Jacod et al. 2009), plus efficient (Jacod and Rosenbaum 2013), noise-and-jump robust (Podolskij and Vetter 2009) and Fourier (Mancino and Sanfelici 2012) approaches. Table 1 compares them on noise robustness, jump robustness and simplicity.

**5. Gap.** According to the paper, no simple estimator is robust to both noise and jumps.

**6. Proposed method.**
1. Clip returns at `τ = c · MAD(r)` (eq. 11), with MAD over the whole sample and c = 3.
2. Local quarticity at time t: a kernel-weighted sum of the past h clipped fourth powers, divided by 3 (eq. 12).
3. Integrated quarticity: the average of the local values over all t (eq. 16).
4. Bandwidth h chosen by leave-one-out cross-validation over {10, 20, 30, 40, 50} (eqs. 17 and 18).

**7. Mathematical formulation.** Eqs. (6) and (7) give the plain kernel version, (11), (12) and (16) the robust version, and (14) and (15) the Gaussian-weight version with noise.

**8. Estimators used.** RQ, BPQ, TPQ, MinRQ, MedRQ, TRQ, kernel and robust kernel. Eleven kernels are listed in Table 2.

**9. Assumptions.** Theorem 1: σ bounded and continuous almost everywhere, `MAD(r) → 0`, `h → ∞` and `h/n → 0`. Theorem 2: σ càdlàg and bounded, the same bandwidth conditions, and uniformly bounded fourth moments of the clipped returns. Noise is i.i.d. normal, added to returns.

**10. Data requirements.** Equally spaced intraday returns. The simulations use 390 one-minute steps (one US trading day). The empirical part uses one-minute closes of five cryptocurrencies.

**11. Algorithm.** Compute MAD, clip, slide a kernel window, average. Cross-validate h on a grid.

**12. Complexity.** O(n h) for the kernel sums, which is linear in n for a fixed h. The full-sample MAD is O(n) on average.

**13. Experiments.**
* Case 1 (Table 3): RQ, BPQ, TPQ and kernel, with h = 30, p = 0.01, σ_J = 0.05 and 1000 paths.
* Case 2 (Table 4, Fig. 2, Table 7, Figs. 3 and 4): adds MinRQ, MedRQ, TRQ and the robust kernel (c = 3), plus a Model Confidence Set test that keeps only the robust kernel. A sensitivity grid over p, σ_J and h.
* Case 3 (Table 8, Figs. 5 and 6): noise variance 0, 0.0005 and 0.001.
* Case 4 (Fig. 7): kernel choice. Logistic is best and Gaussian second.
* Case 5 (Table 9): CV picks h = 20 for the kernel and h = 40 for the robust kernel.
* Real data (Tables 5 and 6, Fig. 8): BTC, ETH, XRP, SOL and DOGE, 1 July to 31 December 2025, with h = 30 and c = 3, and RMSE against an assumed true quarticity of 2.2 × 10⁻⁶.

**14. Results claimed.** The robust kernel has the lowest bias and RMSE in simulation and is the only member of the superior set in the MCS test. It is "consistently" better than the plain kernel. On real data it ranks best overall. For DOGE, MAD is close to zero because of price discreteness.

**15. Limitations (stated and found).** Stated: kernels are light-tailed, robustness comes only from truncation, and DOGE breaks the MAD rule. Found by us: see `docs/estimator_issues.md`. In brief, the kernel weights are not normalised, the truncation bias does not vanish, the averaged estimator equals truncated RQ, the noise model makes noise look like volatility, and the real-data RMSE depends mainly on the assumed true value.

**16. Conclusions.** The robust kernel estimator is simple and outperforms standard estimators across jump intensities, jump sizes, bandwidths and kernels.

**17. Extensions suggested by the paper.** Adaptive (data-driven) bandwidth and kernel, sixth and eighth moments, and a multivariate version.

## Treatment of time in the paper (what is local and what is global)

The paper's local estimator `Q_t` looks only backwards (from t - h to t), which is exactly what a real-time system needs. But two ingredients are global. The MAD uses the whole sample, including the future, and the reported quantity is the average over the whole day. The prototype keeps the backward-looking idea, estimates the MAD from past data only (a calibration window, re-estimated after each alarm; a rolling or exponentially weighted MAD in E3), and uses local, per-minute statistics directly instead of averaging them away.

## What we take from the paper

1. Clipping with a robust scale is simple, fast and works well against large jumps (E2c: relative RMSE 0.24 to 0.30 for the robust kernel variants, against 14 for RQ, with four jumps a day of about 5 sd).
2. The local kernel estimate is a natural streaming statistic.
3. Its bounded output (Theorem 3) is exactly the property that a nonparametric sequential test needs (`docs/math.md`, section 5.3).
