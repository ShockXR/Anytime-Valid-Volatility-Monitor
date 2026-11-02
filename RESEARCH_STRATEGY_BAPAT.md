# Research strategy: Prof. Sudeep Bapat project

## One-line goal

Turn Bapat's jump-robust local volatility measures into an **online monitor** that flags volatility regime changes in real time, with a **proven false alarm guarantee that holds however often you look**.

## What the paper gives us, and what it leaves open

| The paper provides | It leaves open |
|---|---|
| A simple robust recipe: clip returns at c × MAD, then kernel-smooth fourth powers (Bapat 2026, eqs. 11 and 12). | The formulas as printed are not on the scale of integrated quarticity, and the truncation bias does not vanish (`docs/estimator_issues.md`, Issues 1 and 2). |
| A *local* (backward-looking) estimate at every time t. | The paper averages the local values into one number per sample, which turns out to be truncated RQ (Issue 3). The local signal itself is not used. |
| Bounded influence (Theorem 3). | No sequential use. The MAD uses the whole sample, including the future. |

The most valuable and least explored part of the paper is the **local, bounded, robust statistic**. Sequential analysis is Prof. Bapat's core field, so an anytime-valid extension sits naturally between his two interests.

## Research questions

1. **Correctness.** What exactly does the estimator converge to, and how should it be scaled and bias-corrected? *Status: answered (E1, E2, `docs/math.md`, sections 2 to 4).*
2. **Online robustness.** Can the full-sample MAD be replaced by a rolling or exponentially weighted MAD without losing robustness? *Status: yes in simulation (E3).*
3. **Detection.** Can the bounded truncated variance drive a nonparametric e-detector with guaranteed false alarm control under jumps, fat tails and noise, and how fast is it compared with parametric and classical detectors? *Status: yes, guarantee proved (`docs/math.md`, section 5) and confirmed (E4). It costs 10% to 70% extra delay on clean Gaussian data and is the only valid option among those tested once jumps or fat tails appear.*
4. **Why truncation is essential.** *Status: proved that no nonparametric e-variable can detect a decrease of an unbounded nonnegative statistic (`docs/math.md`, 5.3). Truncation is what makes decrease detection possible. This is the cleanest theoretical link to the paper.*
5. **Real data.** Do alarms on real minute data correspond to persistent volatility changes? *Status: on Bitstamp BTC (July to December 2025), 81% of increase alarms and 90% of decrease alarms (band 0.25) are followed by a move in the flagged direction, against 50% at random times (E5c).*
6. **Use.** Does a change-triggered risk rule improve drawdowns out of sample? *Status: not started (see below).*

## Decisions taken

| Decision | Why |
|---|---|
| Use normalised local estimators and the **scaled** MAD, computed **online**. | Removes the scale problem and most of the truncation bias. E3 shows no loss from going online. |
| Monitor the **truncated variance** `min(r², τ²)`, not quarticity. | It is the quantity that moves when volatility changes, and its boundedness gives validity for both directions. Quarticity enters as the variance of the bet. |
| Shiryaev-Roberts e-detector with a fixed grid of bets as the default. | Simple, exact guarantee (Theorems A and B), recursive in O(1) per step. The quarticity-adaptive bet did not help (E4, Q4c). |
| Pick the alarm level A from a per-day false alarm budget (`A = steps per day / 0.05`). | Easy to explain: "at most a 5% chance of any false alarm in a day". |
| Recalibrate after each alarm, and use a null band of 0.25 on real data. | Real volatility drifts. Without a band, the strict null is false and alarms come too often (E4c, E5c). |
| Keep alternative directions in `alternatives/`. | Keeps the main package focused and easy to review. |
| Defer the risk switch backtest. | Six months of one asset cannot support a trading claim, and a weak backtest would distract from the statistical contribution. It is better done after the method is agreed with Prof. Bapat and more data is available. |

## What is genuinely new (as far as our search found)

* Combining a jump-robust, high-frequency local volatility measure with an **anytime-valid** sequential change detector (`docs/literature_review.md`, section 9).
* The impossibility result for decrease detection without a bound, which explains *why* truncation is needed in the sequential setting. It is short and elementary, but we have not seen it stated in this context.
* A reproducibility analysis of the paper's tables. We present this as questions, not corrections.

What is **not** new: the e-detector framework (Shin, Ramdas and Rinaldo 2024), betting bets (Waudby-Smith and Ramdas 2024), threshold estimation (Mancini 2009) and kernel spot volatility (Kristensen 2010).

## Roadmap

**Phase 1, done in this prototype.** Estimators (printed and corrected), simulator, detectors, streaming monitor, 26 tests, experiments E1 to E5, documentation.

**Phase 2, after the first meeting (about 2 to 3 weeks).**
1. Confirm the code conventions with Prof. Bapat (weight normalisation, MAD scaling, simulation σ, σ_J = 0.10 design) and rerun E1 with his settings.
2. Add data: ETH, XRP, SOL, DOGE (as in the paper) and Nifty 50 one-minute bars. Handle discrete prices (DOGE) with a tick-size floor on τ, which is already supported (`scale_floor`).
3. Replace the plug-in baseline with a confidence-sequence baseline (Shekhar and Ramdas 2024), so the guarantee holds with an estimated baseline too.
4. Event study: do alarms line up with scheduled news (FOMC, CPI, RBI policy days for Nifty)?

**Phase 3 (stretch).** A short paper: "Anytime-valid detection of volatility regime changes with robust high-frequency measures". Theory: guarantees, the impossibility result, delay bounds under a change. Experiments: E4 plus Nifty and crypto applications. Then an optional risk switch backtest on several years of data.

## Risks and how we handle them

| Risk | Mitigation |
|---|---|
| The scale findings could come across as criticism of the professor's paper. | Frame them as questions about his code (`PITCH_BAPAT.md`). The prototype builds on his idea and does not replace it. |
| Too many alarms on real crypto data (0.5 to 0.9 per day). | Null band, longer calibration, or a minimum effect size. Alarms are a dial, not a defect. |
| The guarantee needs the true baseline. | Margin (E4d), then confidence sequences (Phase 2). |
| Intraday seasonality is estimated. | Predictable rolling medians. The band absorbs small errors. |
| Only one asset and venue so far. | Phase 2 data work. |
