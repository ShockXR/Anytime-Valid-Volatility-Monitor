# Research strategy: anytime-valid robust volatility monitoring

*Written for a reader who is new to finance and statistics: jargon is explained the first time it appears. Status lines are updated with results as of 3 November 2026.*

## 1. Research objective

Turn the robust, jump-resistant volatility measures of Bapat (2026, International Review of Economics and Finance 109, 105482) into an **online monitor**: a tool that watches a stream of one-minute returns and raises an alarm when volatility (how strongly prices move) shifts up or down, with a **proven limit on false alarms that holds however often the monitor is checked** (an "anytime-valid" guarantee).

## 2. Problem statement

Volatility in high-frequency markets moves between regimes (calm, stressed, calm again). Risk limits, position sizes and model refits should react when it does. Two things make this hard to do online.

1. **Jumps and fat tails.** A few extreme returns distort ordinary volatility estimates. Bapat (2026) proposes clipping each return at c times its median absolute deviation (MAD) and kernel-smoothing the result, which bounds the influence of any single return.
2. **Repeated looks.** A classical test checked every minute will eventually raise a false alarm by chance. E-processes and e-detectors (Shin, Ramdas and Rinaldo 2024) control false alarms at every possible stopping time.

The paper solves the first problem for a fixed sample. It does not address the second: its averaged estimator is a daily summary, and the local quantity it averages is what an online monitor needs.

## 3. What exactly are we validating? (research questions)

| # | Question | Status |
|---|---|---|
| Q1 | **Correctness.** What do the printed estimators converge to, how should they be scaled, and do the published tables follow the printed formulas? | Answered (E1, E2, E5a, E6). The printed kernel is not on the quarticity scale, truncation at c × MAD leaves a bias that does not shrink with n, and the published tables are consistent with the printed (unnormalised) formulas, the raw MAD and a sample with about 2,000 missing minutes. |
| Q2 | **Online robustness.** Can the full-sample MAD be replaced by one computed from past data only, without losing the jump protection? | Yes in simulation (E3). |
| Q3 | **Detection.** Can a bounded truncated statistic drive a nonparametric e-detector with a proven false alarm guarantee under jumps, fat tails and noise, and how fast is it compared with parametric, classical and published change-point methods? | Guarantee proved (`docs/math.md` 5) and confirmed at scale (E4, E9: 1.1% to 1.4% against a 5% bound for four very different return laws). 10% to 70% slower than parametric detectors on clean Gaussian data; the only one of E4's detectors that is valid in both directions under jumps and fat tails (an e-detector on raw r² also stays valid for increases, but truncation is needed for decreases and for speed under jumps: 102 against 246 steps for a 25% increase). Against published methods (E13): BOCPD is not valid under jumps or fat tails (false alarms on 53% to 73% of days at a 1% Gaussian threshold, against 1.2% to 1.3% measured for the plug-in e-detector) but detects large decreases two to three times faster on clean data; a two-window variance-ratio scan (our construction inspired by Bibinger, Jirak and Vetter) misses most small changes with short windows and matches the e-detector offline with windows of 500 returns, at 1.9 to 3.5 times its delay online; on crypto at equal alarm counts the e-detector's alarms validate more often than competitors with 2 to 4 hour windows (78.5% and 80.5% against 58.7% to 71.7%; held out 82.0% and 69.7% against 56.7% to 71.5%), but with day-long windows truncated BOCPD and the variance-ratio scans match it on down alarms (plain BOCPD stays lower) and the best competitor trails by a borderline 6 points on up alarms; on Nifty no difference is significant. |
| Q4 | **Why a bounded statistic.** Why does the paper's truncation matter for sequential testing? | A fall in the variance cannot be detected by any nonparametric e-test (math.md 5.3, related to Bahadur and Savage 1956); a bounded functional such as the truncated variance can. |
| Q5 | **Estimated baseline.** Does the guarantee survive when the "normal" level is estimated from data? | Not with a plug-in estimate: it roughly triples the false alarm rate with short calibration (E4e), gives 41% to 44% over five days where a known baseline would guarantee at most 25% (E10, E14), and admits no bound of the form C m / A at a fixed calibration length (Proposition A.6 and Lemma A.7, E14). Yes, with a proof, with a confidence-sequence baseline (E8, math.md 6), with a band that allows any running truncation level, and with a rank detector that needs no baseline (E14). |
| Q6 | **When does the guarantee stop applying?** | Under persistent volatility clustering (E9) and price-level noise (ALT1), the conditional null is false. Monitoring relative to a volatility forecast restores it when the forecast model is right, at a quantified cost (E10); with a fitted GARCH and a plug-in baseline the one-day rate is 3% to 6%, close to the 4.7% to 5.3% with the true model. Non-overlapping pre-averaged blocks restore it under price noise (ALT1). |
| Q7 | **Real data.** Do alarms on real one-minute data mean anything beyond volatility clustering? | Alarms are followed by volatility on the flagged side of EWMA and HAR forecasts 71% to 83% of the time against about 50% at random times, on five coins in and out of sample (E6), and on Nifty 50 (E7). But a ratio rule tuned with hindsight is about as precise (E6e), so the contribution on real data is a guarantee with fixed settings, not extra precision. Event studies: FOMC statements (E12), Union Budget days (E7). The FOMC and CPI contrast of 2024 and 2025 held with a guaranteed monitor (E22) but was not confirmed on the fresh releases of 2026 (E23, pre-registered: 2 of 16 on BTC, p = 0.25). |
| Q8 | **Use.** Does a change-triggered risk rule help? | Tails on par with daily volatility targeting, worse tracking, no Sharpe difference (ALT2). |

## 4. Hypotheses (stated before the experiments, with outcomes)

| Hypothesis | Outcome |
|---|---|
| H1. The published simulation tables follow the printed formulas with unnormalised weights. | Supported (E1, E5a, E6b), except the large-jump robust rows. |
| H2. A bounded truncated statistic gives a detector whose false alarm rate stays below m / A under jumps and fat tails, while Gaussian detectors fail. | Supported (E4, E9). |
| H3. Sizing bets with the local truncated quarticity speeds up detection. | Rejected: slightly slower in 14 of 18 cells (E4, Q4c). |
| H4. A plug-in baseline costs little. | Rejected for short calibration windows and long horizons (E4e, E10, E14), with a proof that no bound of the form C m / A can hold (E14). |
| H5. Alarms on real data carry information beyond a volatility forecast. | Supported against EWMA and HAR (E6e, E10c); not beyond a well-tuned naive rule. |
| H6. Dividing returns by a volatility forecast restores the guarantee under clustering without losing power. | Half supported: validity restored when the forecast model is right (in simulation the fitted GARCH matches the true model's false alarm rate; a proof also needs a true or confidence-sequence baseline), power against level shifts strongly reduced, with a closed-form explanation (E10). |

## 5. What exactly are we building?

| Item | Specification |
|---|---|
| Inputs | One-minute log returns (a file or a live stream), optional volatility forecast, settings: c, A (false alarm budget), calibration length, band, baseline type. |
| Outputs | Alarm events (time, direction up or down, baseline level, delay since calibration), segment summaries, per-step statistics on request. |
| Mathematical model | Truncated squared return X_t = min(r_t², τ²) with τ = c × scaled MAD from past data; nulls on the conditional mean of X; Shiryaev-Roberts betting e-detectors for increases and decreases; plug-in, margin or confidence-sequence baseline. |
| Algorithms | O(1) recursive updates per step (about 10 microseconds); calibrate, monitor, alarm, recalibrate loop; predictable intraday deseasonalisation. |
| Data | Synthetic jump diffusions, Student t, ticks, GARCH, price noise; Bitstamp BTC/USD, Binance spot BTC, ETH, XRP, SOL, DOGE (the paper's five coins), Nifty 50, India VIX (`data/MANIFEST.csv`). |
| Evaluation metrics | False alarm rate per day against the proved bound, detection probability and delay, share of alarms followed by volatility on the flagged side of a forecast against random-time and rate-matched controls, event alignment, runtime. |
| Visualisations | Alarm timelines on price and volatility charts, false alarm and delay tables, event profiles (`results/figures/`). |
| Demonstration | `python -m avmonitor simulate | run | estimators`, `notebooks/demo.ipynb`. |
| Success criteria | Section 9. |

## 6. Minimum scientifically credible implementation, and extensions

1. **Core (done).** Estimators as printed and normalised; the e-detector with a proof; simulation of false alarms and delay against benchmarks; one real-data run.
2. **Research extensions (done or in progress).** Five coins and Nifty; the confidence-sequence baseline (E8); stress tests (E9); forecast-relative monitoring (E10); event studies (E7, E12); delay theory (E11); running truncation (E14); panel monitoring (E15); an adversarial search for null streams (E16); several sampling frequencies in one monitor (E17); benchmark against published change-point methods (E13).
3. **Advanced computational improvements.** Not needed: the monitor runs at about 10 microseconds per step in plain NumPy (E9d).
4. **Optional demonstration features.** CLI, executed notebook; a live dashboard was judged unnecessary.

## 7. What would make this interesting to Prof. Bapat

* Sequential analysis is his main field; the prototype puts his estimator into that setting, with proofs.
* An independent reproduction of his tables on the five coins, with precise questions where they do not reproduce (asked to learn, not to correct).
* A clear theoretical role for truncation in sequential testing.
* Real-data applications on crypto and on an Indian index.
* Several natural joint-paper directions (section 11).

The main project does not claim to improve his estimator. Where the normalised robust kernel does well (E1: it beats the plain kernel in most cells once both are on the same scale), we say so. One alternative (ALT3) explores a different quarticity estimator that, in simulation, is somewhat more efficient than the κ-corrected robust kernel under constant volatility and stays consistent when volatility moves within the day; it is a proposal for discussion, not a correction.

## 8. Choosing the primary direction

Five candidate directions were scored from 1 (weak) to 10 (strong) on the seven criteria of the project brief. "Ease" scores implementation complexity the other way round, so that 10 means simplest. Scores are judgements, explained below the table; they are not measurements.

| Direction | Mathematical significance | Relevance to Prof. Bapat | Feasibility | Empirical value | Ease | Demonstration quality | Collaboration potential | Total (of 70) |
|---|---|---|---|---|---|---|---|---|
| **P1. Anytime-valid robust volatility monitor (chosen)** | 8 | 9 | 8 | 7 | 7 | 9 | 9 | **57** |
| P2. Reproduction and correction study of the estimator only | 5 | 8 | 10 | 4 | 9 | 5 | 6 | 47 |
| P3. Efficient, consistent robust IQ estimator (alternative ALT3) | 8 | 9 | 6 | 6 | 6 | 5 | 8 | 48 |
| P4. Change-triggered risk switch (alternative ALT2) | 3 | 4 | 8 | 6 | 8 | 7 | 3 | 39 |
| P5. Jump-robust correlation-break monitor (alternative ALT4) | 7 | 7 | 6 | 6 | 5 | 7 | 6 | 44 |

**Why these scores.**

* **P1, mathematical significance 8.** It needs real proofs (test supermartingales, Shiryaev-Roberts bounds, the impossibility result for decreases without a bound, the confidence-sequence baseline) rather than only algebra. Not 10, because the e-detector framework itself is known (Shin, Ramdas and Rinaldo 2024); the contribution is the combination and the role of truncation.
* **P1, relevance 9.** It uses exactly the part of the paper that is most original (a bounded, robust, local statistic) and it lives in sequential analysis, Prof. Bapat's main field.
* **P1, feasibility 8 and ease 7.** Everything runs on a laptop in minutes. The hard parts are mathematical, not computational. Real data had to come from public mirrors, which cost a point.
* **P1, empirical value 7.** Real-data validation is possible (five coins, Nifty, several years), but "was this alarm a true regime change?" has no ground truth, so validation is indirect (forecasts, events, rate-matched rules).
* **P1, demonstration 9.** A live monitor with alarms on real charts is easy to show and easy to understand.
* **P1, collaboration 9.** It naturally becomes a joint short paper: theory (guarantees with an estimated baseline, the role of a bounded statistic, the forecast-relative signal cap), simulation, and a crypto and Nifty application.
* **P2.** Very feasible and relevant, but on its own it reads as a list of corrections to the author's paper, which is a poor opening for a collaboration and has little new science. Its content is kept as the "questions about the paper" part of the pitch.
* **P3.** Strong estimation question in the author's own lane, but the efficient-estimation theory already exists (Jacod and Rosenbaum 2013), so the novelty is narrower and the demonstration is a table, not a tool. Kept as the leading alternative.
* **P4.** Easy to code and visually appealing, but it is an application, not research, and a credible trading result needs much more data and care with costs. Built as an alternative with an honest risk-management framing; it showed no Sharpe gain.
* **P5.** The paper names the multivariate extension as future work, and correlation breaks matter. It needs more new theory and data handling. Built as an alternative; the sign-concordance statistic turned out to be simple and valid.

**Decision.** P1 is the primary prototype. P3 is the recommended follow-up if Prof. Bapat prefers estimation over detection. P4 and P5 are working minimal prototypes in `alternatives/`.

## 9. Success criteria (set before round 2, with outcomes)

| Criterion | Met? |
|---|---|
| Every theorem has a complete proof with stated assumptions, checked by an independent reviewer. | Yes for sections 5 and 6 of `docs/math.md` (section 5: review round 1; section 6: review round 2a, no errors found, a few statements sharpened). The delay theory (`docs/theory_delay.md`) and the panel theory (`docs/theory_panel.md`) were checked in review round 2b: no proof step was wrong, but the 73-step decrease floor holds for the default bets only (58 steps for any valid bets), the rounded-price test of c_down was mis-modelled (c_down = 2 is slower for a 20% fall at a tick of 0.5 sd), and the optimality claim was narrowed and credited to Ram and Ramdas (2026); all fixed. The running-truncation theory (`docs/theory_running_tau.md`) is in review round 3. |
| The false alarm rate stays below the proved bound in every null scenario the theorem covers, with Monte Carlo error small enough to see the slack. | Yes (E9: 4000 paths; E8: 10,000 paths). |
| Every published number we discuss is either reproduced or has a stated, tested hypothesis for why not. | Yes for Tables 5 to 8, with open questions listed (E1, E6). |
| Real-data alarms beat a random-time control and a forecast, out of sample, on several assets. | Yes (E6, E7); not beyond a hindsight-tuned naive rule (E6e). |
| A researcher can rerun every table from a clean clone. | In progress (final audit, iteration 10 of the log). |

## 10. Risks and how we handled them

| Risk | Handling |
|---|---|
| Findings about the paper could read as criticism. | Framed as questions about the code, with evidence (`PITCH_BAPAT.md`). The prototype builds on the idea and does not replace it. |
| Overstated claims. | Independent review round 1 corrected the wording of the impossibility result, the role of the kernel, noise robustness and the strength of the real-data check. |
| The guarantee needs the true baseline. | Confidence-sequence baseline with a proof (E8); plug-in effect measured (E4e, E10). |
| Real volatility clusters, so the null is often false. | Stated precisely (`docs/limitations.md`); forecast-relative monitoring (E10). |
| Data come from third-party mirrors. | Loaders report gaps and duplicates; SHA-256 manifest pins the bytes; cross-venue check BTC Binance against Bitstamp (E6d: 72% of alarms matched within an hour against 3% by chance). |
| The monitor is blind while recalibrating. | Found by the event studies (E7, E12: about 70% of the time on crypto). E14 gives monitors with no blind period and a proof; on E12's releases they watch 114 of 114 instead of 55, at the price of more alarms per day without a band. |

## 11. Development roadmap

**Done in rounds 1 to 3.** Everything in sections 3 to 6, about 500 automated tests, experiments E1 to E23, theory work streams TH4 to TH6 (robustness to drift, limits for every valid procedure, per-day guarantees), alternatives ALT1 to ALT6, eight rounds of independent review and a referee report on the paper with a response, a one-command reproduction script (`scripts/reproduce.py`), a short paper draft with a supplement (`paper/`), a release-ready package (version 0.3.0, not released) and a replay page of two landmark days.

**Next, after a first meeting with Prof. Bapat (about 2 to 3 weeks).**
1. Confirm the code conventions (weight normalisation, MAD scaling, simulation σ, the σ_J = 0.10 design, the RMSE definition, the DOGE data) and rerun E1 and E6 with his settings.
2. Decide the target with him: variance stream, quarticity stream (ALT1), or both.
3. Finish the theory that is still open: the alarm rate of the awake monitors on real data without a band, and why they flag fewer of the releases at which they are awake (E14); the chance of any false discovery in the panel once real changes are present (E15); and a central limit theorem for the ALT3 estimator.

**Then (stretch).** A short paper, "Anytime-valid detection of volatility regime changes with robust high-frequency measures": guarantees with an estimated baseline, the bounded-statistic argument, the forecast-relative signal cap and information bound, and the crypto and Nifty applications. In parallel, the estimation paper of ALT3 if he prefers that direction.
