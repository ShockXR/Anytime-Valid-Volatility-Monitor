# Alternative prototype directions

These are other ways to build on Bapat (2026). They are kept apart from the main package on purpose. None of them is implemented in `src/avmonitor`. Each entry says what it is, why it might matter, what it would take, and how we rate it.

## A1. Noise-robust first stage: pre-averaging, then truncation

**What.** Average returns over small blocks (Jacod et al. 2009) before applying the MAD clip and kernel.
**Why.** Truncation does not remove microstructure noise under either noise model (`docs/estimator_issues.md`, 4b). Pre-averaging does, at the cost of a slower rate.
**Effort.** Medium. Pre-averaged returns overlap, so the bias corrections and the e-detector's independence-free null need care. The null remains a conditional moment statement, but the statistic is then autocorrelated by construction.
**Rating.** High value for liquid assets sampled at seconds. Low value at one minute for BTC, where the noise ratio is small.

## A2. Confidence-sequence baseline (Shekhar and Ramdas 2024)

**What.** Replace the plug-in `m₀` with a running confidence sequence for the pre-change truncated variance, and declare a change when the post-change estimate leaves it.
**Why.** Makes the false alarm guarantee hold *with* an estimated baseline (`docs/math.md`, 5.6).
**Effort.** Medium. Betting confidence sequences for bounded means already exist (Waudby-Smith and Ramdas 2024).
**Rating.** Highest priority theoretical upgrade. Planned for Phase 2.

## A3. Monitor quarticity itself

**What.** Run the same e-detector on `min(r⁴, τ⁴)`, the truncated local quarticity stream.
**Why.** Detects changes in vol-of-vol or tail behaviour. This is closer to the paper's own target and may interest Prof. Bapat.
**Effort.** Low, since the code is generic in the statistic.
**Rating.** Worth a quick experiment. Expect slower detection, because fourth powers are noisier.

## A4. Two-sample betting test between windows

**What.** Bet on whether recent truncated squared returns come from the same distribution as the calibration window (Shekhar and Ramdas, IEEE Trans. Inf. Theory 2024).
**Why.** Detects any distribution change, not just a shift in variance.
**Effort.** Medium to high.
**Rating.** Interesting but less interpretable for risk users.

## A5. Bias-corrected local plug-in for IQ (Jacod and Rosenbaum 2013)

**What.** Estimate IQ as `Σ Δ σ̂_t⁴` from local variance estimates, with the standard bias correction, using MAD-truncated returns.
**Why.** The principled way to aggregate local estimates (rather than averaging local fourth powers, which is truncated RQ). Could yield an efficient and robust IQ estimator.
**Effort.** Medium. The theory exists for threshold-truncated versions.
**Rating.** Good *estimation* follow-up paper with Prof. Bapat, separate from detection.

## A6. Risk switch backtest

**What.** Scale a long BTC (or Nifty futures) position by `baseline sd / current regime sd` after each alarm. Compare with buy and hold and with trailing volatility targeting.
**Why.** Shows practical value.
**Effort.** Low to code. A credible result needs several years and several assets, plus transaction costs.
**Rating.** Do it last, and only with enough data.

## A7. Multivariate monitor

**What.** Truncated covariance (or a robust correlation) e-detector across several assets.
**Why.** The paper names a multivariate extension as future work. Correlation breaks matter for portfolios.
**Effort.** High.
**Rating.** Longer term.
