# Anytime-Valid Volatility Monitor

**An online monitor that watches one-minute returns and raises an alarm when volatility shifts up or down, with a proven limit on false alarms that holds however often you look at it.** It puts the MAD-based truncation used by Bapat (2026), an idea for jump robustness that goes back to Mancini (2009), inside e-detectors, a nonparametric framework for sequential change detection (Shin, Ramdas and Rinaldo 2024).

> Status: research prototype, round 3 (November 2026). 656 tests, 25 experiments (E1 to E25) on simulated and real data, three theory work streams (TH4 to TH6), eight alternative directions, eleven rounds of review by separate AI agents, not people (the ninth an audit of the package code, the tenth of the blind-segment guard, the eleventh of the anytime-valid release test) and three AI referee passes on the paper, proofs in [`docs/math.md`](docs/math.md) and the theory notes, and a short paper draft in [`paper/`](paper/). For the short version read [`PITCH_BAPAT.md`](PITCH_BAPAT.md); for every number read the [experiment report](docs/experiment_report.md).

Every result in this repository carries a label: **[PUBLISHED]** (taken from a paper), **[SYNTHETIC]** (our simulation) or **[REAL]** (our run on market data). Nothing is reported that was not run.

**How this work was done.** Shreyan Thanki directed the project; the code, experiments, proofs and write-ups were produced with extensive AI assistance (Claude, by Anthropic). Every review and referee report in `docs/review/` was written by a separate AI agent instructed to act as a sceptical referee, not by a human expert, so "independent" there means "not the agent that wrote the files". The proofs have therefore not been checked by a human statistician. "Pre-registered" means a plan committed to git before the data it tests were read; there is no external registry. Every number is produced by code in this repository, and `scripts/trace_numbers.py` checks that the numbers printed in this file, the pitch and the short paper appear in a committed result or report.

## What the project does

1. **Checks the estimator.** It re-implements the robust kernel quarticity estimator of Bapat (2026) exactly as printed, finds what it converges to, and reproduces the published simulation and crypto tables.
2. **Turns the estimator's key idea into a monitor.** The paper clips each return at c times its median absolute deviation (MAD). We monitor the clipped squared return X_t = min(r_t², τ²) with a betting e-detector that alarms when its running wealth passes a threshold A.
3. **Proves and measures the false-alarm guarantee.** Average time to a false alarm is at least A steps, and the chance of any false alarm within m steps is at most m / A, under jumps, fat tails and any return distribution, as long as the conditions in [`docs/math.md`](docs/math.md) section 5 hold. We also show where those conditions fail on real data and what fixes them.
4. **Applies it to real markets:** Binance BTC, ETH, XRP, SOL and DOGE (the paper's five coins), Bitstamp BTC/USD and the Nifty 50 index, with event studies on FOMC statements, US CPI releases and Indian Union Budget days.

## Research question

*Can a robust, jump-resistant volatility measure drive an online detector of volatility changes whose false alarm rate is guaranteed in finite samples, without a model for the returns, and what does that guarantee cost in detection speed?*

Sub-questions (status in [`RESEARCH_STRATEGY_BAPAT.md`](RESEARCH_STRATEGY_BAPAT.md), section 3): correctness of the printed estimators (Q1); online truncation from past data only (Q2); validity and delay against classical and published detectors (Q3); why a bounded statistic is needed (Q4); estimated baselines (Q5); when the guarantee stops applying (Q6); meaning of alarms on real data (Q7); use in a risk rule (Q8).

## Key results

* **The guarantee holds where classical detectors break [SYNTHETIC].** With jumps or fat tails, Gaussian likelihood-ratio and CUSUM detectors raise a false alarm within a day on 70% to 97% of paths. The robust e-detector stays at 0.3% to 2.3% against its 5% bound, and detects a 50% volatility increase in a median of 36 to 49 steps (E4). Across four very different return laws the measured rate is 1.1% to 1.4% (E9), so the bound is valid but about four times conservative for ordinary returns. An adversarial search over null streams that satisfy the theorem's conditions (E16) found nothing above the bound: the worst stream found, which watches the detector's state, always uses the widest legal spread and, when one draw can reach the threshold, lands exactly on it, reaches 2.62% within the first day after calibration (95% interval 2.50% to 2.75%) against 5%, and 18.1% over five days against 25%. Outside the conditions every failure matched the condition that was broken (clustering, an estimated baseline, coarse price ticks) and disappeared when it was restored (for ticks, the same values drawn independently).
* **How fast the guarantee degrades when the null holds only approximately [proof and SYNTHETIC].** If the expected clipped squared return can exceed its baseline by at most a relative drift δ, the probability of a false alarm within m steps is at most m/A times a closed-form inflation factor, attained by a deterministic stream where the bound reaches 1. The one-day 5% bound doubles at a drift of 1.63% and reaches 1 at about 3.7% upward drift. A drift that alternates costs far less, but only a pathwise bound shows it, and that bound is proved only when the drift does not respond to the data (or in a looser form). The plug-in baseline's estimation error acts as a persistent drift: its proved bound is 42% for one day with 390 calibration steps, within a factor of 1.3 of the worst stream found in E16 (TH4, [`docs/theory_robustness.md`](docs/theory_robustness.md)).
* **How fast any valid procedure can be, and when the plug-in keeps its level [proof and SYNTHETIC].** A deterministic rule with the guarantee of Theorem B, even one that sees the raw returns, can stop on a path only where a product of per-step likelihood bounds reaches A / n; a randomised rule can stop earlier only with a correspondingly small probability. So on a complete collapse of volatility at c = 3 and A = 7800 no valid deterministic rule alarms before step 45 (a randomised one has a mean delay of at least 43.4 steps on that path); the e-detectors need 58 steps and the default grid 73. Over 18 post-change laws the default grid's mean delay is at most 1.7 to 2.2 times that of the best valid deterministic rule, and the best single bet's at most 1.3 to 1.8 times, because these ratios are taken against a lower bound that may not be attained. For the plug-in baseline, every stream whose conditional mean matches the return law keeps the nominal level, up to a factor 1 + O(m / √n), when the horizon m is short compared with the square root of the calibration length n; when m is long compared with both √n and log A, some such stream breaks it by an unbounded factor. A one-day horizon after a one-day window is already in the second regime: the worst case lies between about 33% and 42% against a nominal 5% (averages over 20,000 simulated windows of bounds proved for each window), while i.i.d. returns give at most 5.1%. A margin of Φ⁻¹(1 - δ / 2) standard errors on the baseline restores m / A + δ. Finally, the growth-optimal decrease bet sits at its upper limit only for light tails (TH5, [`docs/theory_limits.md`](docs/theory_limits.md)).
* **What the guarantee says about each day, and a monitor with a bound for every day [proof and SYNTHETIC].** Theorem B bounds the chance of a false alarm within the first day after each (re)calibration; on later days only the long-run average of 0.05 false alarms per day is proved (Corollary C), and a constant null stream makes the default detector alarm with certainty on day 40. Measured on later days, the conditional chance of a false alarm was 4.7% on average over days 11 to 20 for i.i.d. Gaussian returns, 5.25% for a null stream that watches the detector, and 28% on a day after nine calm days. With a known baseline, restarting the statistics at every day boundary (`PeriodicVolatilityMonitor(baseline="known")`) gives a proved 5% for every day, given everything before it, at a cost of 6% to 44% in mean delay for 390-step days and 2% to 14% for 1440-step days; with an estimated baseline the restart (`PeriodicVolatilityMonitor`) carries a per-day proof only through the confidence-sequence baseline (δ + m / A for each day of a segment, not given the past, and only while Theorem D's tolerance holds to the end of that day), and none with the plug-in. The confidence-sequence baseline's theorem is now stated with a tolerance that does not depend on the bets, plus a version with no assumption at all (TH6, [`docs/theory_per_period.md`](docs/theory_per_period.md)).
* **Why truncation matters for sequential testing.** No nonparametric e-test can detect a *fall* in the variance (a one-sided form of Bahadur and Savage 1956, [`docs/math.md`](docs/math.md) 5.3). A bounded functional can, and the truncated variance is one with a direct volatility meaning.
* **Estimated baselines need care [SYNTHETIC].** Plugging in an estimated "normal" level roughly triples false alarms with a 390-step calibration window and gives 41% to 44% false alarms over five days, where a known baseline would guarantee at most 25% (E4e, E10, E14); over 30-day streams it raises 0.084 to 0.094 false alarms per day against a nominal 0.05. We proved that at a fixed calibration length no bound of the form C m / A can hold for it, for each of the four simulated return laws of E14 (Proposition A.6 and Lemma A.7 in [`docs/theory_running_tau.md`](docs/theory_running_tau.md)); how often a calibration window makes it drift (77% to 83% of windows at 390 steps) is a numerical estimate. A confidence-sequence baseline restores a proof, P(false alarm within m) ≤ δ + m / A, at a cost of 1.8 to 7.6 times the delay with short calibration (E8).
* **A monitor with no blind period [SYNTHETIC and REAL].** After each alarm the default monitor spends a day recalibrating, and on real crypto data it is blind about 70% of the time (E12). A rank-based e-detector, restarted with an empty reference set right after each alarm, keeps an exact m / A bound with no baseline at all (0.035 false alarms per day in simulation, against a proved 0.05 per day) and catches 98% to 99.6% of changes in streams where volatility moves about once a day, against 56% to 66% for monitors that recalibrate for a day after each alarm (65% to 66% for the default plug-in). Rerun on the FOMC and CPI design of E12, it watches all 114 releases instead of 55; this rerun was done after E12, so it is a check, not a pre-registered test (E14). On the real coins its up alarms were confirmed less often than the plug-in monitor's (6.7 points fewer, Holm p = 0.006; E20).
* **Against published methods [SYNTHETIC and REAL].** We compared the monitor with Bayesian online change-point detection (BOCPD, Adams and MacKay 2007) and a two-window variance-ratio scan inspired by Bibinger, Jirak and Vetter (2017), which is our construction and not their procedure (the paper could not be obtained from our environment). At thresholds that give 1% false alarms per day on Gaussian data, BOCPD raises a false alarm on 73% of days with jumps and 53% with Student t5 returns, while the e-detector stays at 1.2% to 1.3% (plug-in baseline, measured; 0.0% to 0.1% with the confidence-sequence baseline, whose 5% bound is proved). On clean Gaussian data BOCPD detects large decreases two to three times faster (26 against 81 steps for a halving of the sd, 51 against 109 for sd x 0.67); for sd x 0.8, and with ticks or price noise, it is 1.0 to 1.9 times faster. Run on clipped returns it keeps that speed under jumps but not under fat tails (7.5% false alarms per day). With windows of 50 to 125 returns the variance-ratio scan finds a 25% change in 6% to 30% of offline windows against 100%; with windows of 500 returns it finds 99%, so the gap is window length, but run online it then needs 1.9 to 3.5 times the e-detector's delay. On the five coins, with every method tuned to raise the same numbers of up and down alarms, the e-detector's alarms are followed by a day beyond an EWMA forecast more often than those of competitors with windows of 2 to 4 hours: 78.5% (up) and 80.5% (down) against 58.7% to 71.7% in Jul to Dec 2025, and 82.0% and 69.7% against 56.7% to 71.5% in a held-out half-year, where BOCPD on clipped returns is not significantly different on down alarms (67.0% against 69.7%, p = 0.47); the other differences hold with tests clustered by week (p ≤ 0.010). With windows of half a day to a day, fixed before the run, BOCPD on clipped returns and the variance-ratio scan with day-long windows match the e-detector's down alarms (plain BOCPD stays lower), and the best competitor, picked after seeing the results, trails its up alarms by about 6 points (p = 0.028 to 0.065 depending on the test, about 0.2 if the two methods are treated as independent samples), and on a four-hour outcome the short-window competitors do as well or better. On Nifty no difference at the first-run windows is significant, and longer windows do not help the competitors (E13).
* **Stronger benchmarks: the advantage comes from the statistic and the design, not from Ville's threshold [SYNTHETIC and REAL].** At the e-detector's own false alarm rate on Gaussian data, a threshold calibrated by simulation makes the same betting statistic no faster for increases and about 10% faster for decreases, and a CUSUM on the same truncated statistic is 15% to 40% faster for decreases; so Ville's threshold costs little delay at a matched rate, while thresholds that use the full 5% are faster still. Thresholds calibrated on Gaussian data do not carry over to fat tails even on the truncated statistic (CUSUM 5.0% to 8.2% against a 1.3% target under Student t), and a bootstrap of the calibration window repairs this for i.i.d. returns. A Student t BOCPD removes plain BOCPD's false alarms under jumps (0.7% against 73%) but finds only 56% of 25% increases (e-detector 98%). On the five coins, at matched alarm counts, CUSUM and Shiryaev-Roberts rules on the same truncated statistic are as precise as the e-detector (never significantly worse in 24 comparisons), while a Student t BOCPD and CUSUM monitoring with a training sample are less precise in most comparisons (unadjusted p-values). The e-detector's real advantage is that its threshold needs no calibration and carries a proof (E21).
* **Volatility clustering breaks the null, and fixing it caps the signal [SYNTHETIC].** Under persistent GARCH clustering the raw monitor alarms on 82% to 100% of quiet days. Dividing returns by a fitted GARCH forecast brings this back to 3% to 6% per day, close to the rate with the true forecast model, but a forecast that keeps adapting also absorbs most of a real level shift. A doubling of the long-run variance raises the ratio of squared returns to the forecast by only 9% at persistence 0.99, and no shift, however large, can push that ratio above (1 − β)/α = 1.2. We derived this cap in closed form and checked it by simulation (E10).
* **The guarantee on real data [REAL marginal and REAL].** On streams that permute the one-minute returns of calm real days, so that they keep the real ticks, zeros and heavy tails but are exchangeable, the proved monitors never exceeded a bound in 252 cells: within the first day 2.26% with the known baseline, 0.01% with the confidence-sequence baseline and 3.38% for the rank detector, against 5%; the plug-in baseline gave 5.50% within the first day and 41% over five days against a bound of 25% (pre-registered, E19). Real calm days are still not a null at the minute scale: resampled in blocks they make every monitor alarm within a day on 86% to 100% of streams. Used as the primary configuration on the five coins, the confidence-sequence baseline raises about a quarter fewer alarms, keeps the share of up alarms followed by higher volatility (78% against 80%) and raises it for down alarms (96% against 88%); its matched up alarms come 4 to 10 minutes later, and because it recalibrates less often it watches more of the time (46% to 49% of minutes against 30% to 33% at band 0.25); in 5 of 7031 segments, all on BTC, its increase side could not alarm at all, once for 136 days (found after the run; the guard of E24 removes such segments with a proof).
* **Alarms on real data beat random times but not a strong forecast; the main contribution is the guarantee [REAL].** On the five coins, alarms are followed by volatility on the flagged side of EWMA and HAR forecasts 71% to 83% of the time, against about 50% at random times, in and out of sample (E6). A simple ratio rule tuned with hindsight is about as precise (E6e), so what the e-detector adds over the simplest rule is a false-alarm guarantee with fixed settings; the published methods of E13, tuned to the same alarm counts, are less often followed by a day beyond the forecast on crypto when their windows are 2 to 4 hours long, but with day-long windows BOCPD on clipped returns and the variance-ratio scan match the e-detector on down alarms and the best of them trails it by a borderline 6 points on up alarms. A re-analysis with Holm adjustment within families and a bootstrap that resamples calendar days jointly across coins (E20, methods fixed before its results) confirms every comparison of alarms with random times or controls (Holm p ≤ 0.0016), but only 1 of 48 comparisons with the ratio rule, and the 6-point margin does not survive adjustment. Against a HAR forecast that includes the last hour's realised variance, increase alarms add no significant information about next-day realised variance (Holm p = 0.32) and decrease alarms a little (about 6% lower, Holm p = 0.046), no more than a 240-minute ratio rule, and scaling the forecast at alarms does not lower its out-of-sample QLIKE loss. Of 36 real-data claims, one rests on a pre-specified primary test (FOMC on BTC, which holds; CPI in the same family does not), and that one was not confirmed on fresh releases (E23, below).
* **Many assets at once [SYNTHETIC and REAL].** A Shiryaev-Roberts statistic divided by time is not an e-value, so the usual multiple-testing recipes do not apply directly. Per-asset e-processes with e-BH give a proved false discovery rate at every stopping time (at most 0.96% observed in 48 designs). On five coins plus Nifty its discoveries are confirmed 78% to 97% of the time against about 50% at random times; its gain over a Bonferroni split is mostly earlier flags, not more of them (E15, `docs/theory_panel.md`).
* **One monitor for several time scales [proof, SYNTHETIC and REAL].** Averaging e-detectors run on 1, 5 and 15 minute returns needs care: a coarse detector changes only at the end of its block, so between block ends it is not a fair bet, and the usual proof does not cover the plain average checked every minute (a union bound proves three times the budget for balanced weights). Counting each coarse detector only at the wealth it is sure to keep gives a proved bound, P(false alarm within m minutes) ≤ Σ w_k ⌈m/k⌉ / A, whenever each frequency's own null holds with the true baselines (for example for independent returns with any tails); with the plug-in baselines used on real data there is no proof, and in simulation the one-minute detector and the equal-weight average then exceed 5% per day. With independent simulated returns the combination never beats the best single frequency, but stays within 5% (equal weights) to 18% (balanced weights) of it on persistent changes, and it catches 95% to 97% of paths (Gaussian returns; 84% to 89% for Student t3) after a jump in the lag-one autocorrelation of returns from 0 to 0.3, which leaves the one-minute variance unchanged and which the one-minute detector almost never sees (3.8%). On five coins it behaves like the one-minute monitor, and at matched alarm counts the 5 and 15 minute detectors alone are no more precise than it (a check added after review; E17, [`docs/theory_multifreq.md`](docs/theory_multifreq.md)).
* **Events in 2024 and 2025 [REAL].** Around the 16 FOMC statements of 2024 and 2025, the streaming monitor (the design fixed in advance) raised an increase alarm within six hours at 9, 8 and 9 statements for BTC, ETH and SOL, against 11% to 17% of control windows on the same weekday at the same clock time. Those are exactly the statements at which it was monitoring rather than recalibrating. The pre-specified primary test (BTC, FOMC) passes after Holm correction (p = 0.0002); the CPI test on BTC does not (p = 0.11). That run used the plug-in baseline, which has no finite-sample guarantee. Rerun under a pre-registration with the confidence-sequence baseline restarted every UTC day, which is proved to keep the chance of a false alarm on each day at or below 5% (marginally, while volatility stays within the ±25% tolerance of its calibrated level), the primary test still passes: 8 of 16 statements against 10.9% of matched control windows, Holm p = 0.0003 (a rerun of E12's test on the same data, so designed after the E12 result was known); this monitor was watching at 80 of the 114 coin and release pairs against the plug-in's 55 (E22). The same monitor alarmed on about half of its monitored days, so real days break the tolerance: E22 shows a monitor that carries the guarantee finding the releases, not the guarantee being exercised. A monitor calibrated on the day before each statement flags 14 to 15 of 16, against 29% to 32% of controls (E12). On Nifty 50, 9 of 11 Union Budget sessions raise an alarm against 25% of ordinary sessions (E7).
* **On fresh releases the FOMC contrast was not confirmed, so it is exploratory [REAL].** In a pre-registered test on the 16 FOMC statements and CPI releases of January to October 2026, which no analysis in the project had seen, the same monitor raised an up alarm within six hours of 2 of 16 releases on BTC against 5.9% of matched control times (p = 0.25, planned power 0.82; E23). Two things stand out: the confidence-sequence monitor could not raise an increase alarm for 28% of its monitored time on BTC (one stretch of 47 days, the blind mode of E19), and the rank monitor, which is never blind, still flagged 4 of the 6 FOMC statements on each coin (a pre-registered secondary result, Holm p = 0.016 to 0.055). No CPI contrast held in 2026 under any monitor. A guard proved in E24 (discard a blind calibration window and calibrate again, splitting the error budget over the windows tried) removes every blind segment on BTC, but on these already-seen data it does not rescue the test (4 of 16, p = 0.056, post hoc). The question now has a pre-registered forward test with a fixed date: the rank monitor on the BTC FOMC statements and CPI releases from 25 October 2026 to 1 June 2027, with power 0.82 at the planning rates (`docs/preregistration_e24_future.md`). A second pre-registered test asks the same question anytime-validly (E25): each release, with its own control windows within four weeks that no other stratum's window overlaps, contributes an e-value, and their running product may be checked after every release and continued without an end date; if releases are like their controls, it ever reaches 20 with probability at most 5%. At the planning rates it is less powerful than the fixed-date test on 1 June 2027 (0.43 against 0.82) but reaches 0.82 by 2028 and 0.98 by 2029 [SYNTHETIC]; its first release enters on 26 November 2026. On data already seen it would have crossed 20 in mid 2024 on all three coins and ends at 3.8 on BTC in 2026 (post hoc; `docs/preregistration_e25_sequential.md`, `docs/theory_sequential_events.md`).
* **Landmark shocks: an honest case study [REAL].** In four pre-registered case studies (the Indian election results of 4 June 2024, the sell-offs of 5 August 2024 and 7 April 2025, and the crypto liquidation cascade of 10 October 2025), the streaming monitor flagged all three Nifty 50 shocks 11 to 172 minutes after the open, but none of the fifteen crypto coin and event pairs within 24 hours: on each crypto event it had raised an up alarm 4.5 to 6.9 hours earlier, as the turbulence began, and was still recalibrating. Its competitors, at matched alarm rates, were mostly paused for the same reason. With a fresh calibration just before each event it flagged all 18 cases in 15 to 207 minutes. The lesson is the blind period after an alarm: in an exploratory run decided after seeing this, E14's monitors with no blind period, at the same null band and a similar alarm rate, flagged 15 and 17 of the 18 cases within 24 hours, within about an hour on 5 August 2024 and 10 October 2025 but late or not at all on 7 April 2025 (E18; case studies, not evidence).
* **One live session, settings frozen before the open [REAL].** Every setting of a live run on the Nifty 50 session of 3 November 2026 was committed before the open (`docs/preregistration_nifty_2026-11-03.md`). The pipeline ran end to end on bars it had never seen, but the monitor raised no alarm because it could not: a down alarm at 11:42 IST the day before had started a recalibration of 750 returns, which lasted past the close. The dry run on 2 November had shown that alarm, so this was foreseeable, and we missed it. The run shows the pipeline working and, once more, the blind period after an alarm; it is not evidence (`results/live/2026-11-03/`).
* **Reproducing the paper [PUBLISHED vs SYNTHETIC].** The published tables are consistent with the printed formulas, unnormalised kernel weights, the raw (unscaled) MAD and a crypto sample with about 2,000 missing minutes. Two findings we would like to discuss: the printed kernel estimators are not on the integrated-quarticity scale, and truncation at c × MAD leaves a bias that does not shrink with n. Both have simple fixes under constant volatility (normalise the weights; divide by κ₄). When volatility moves within the day, the truncation also needs a local scale (E2b, ALT3; [`docs/estimator_issues.md`](docs/estimator_issues.md), E1, E6).

The main project does **not** claim to improve the estimator. Once both kernels are put on the same scale, the robust kernel beats the plain kernel in most simulation cells (E1). One alternative direction (ALT3, `alternatives/robust_iq_estimator/`) explores a different quarticity estimator: in simulation it is somewhat more efficient than the κ-corrected robust kernel under constant volatility (0.95 against 0.86; the kernel with an oracle cut-off reaches 0.98) and stays consistent when volatility moves within the day. It is a proposal for discussion, not a correction.

## Method in brief

| Step | What | Where |
| --- | --- | --- |
| 1 | Truncation level τ = c × MAD × 1.4826, estimated on each calibration window from past data only and held fixed until the next recalibration (a running level is studied in E3 and E14), c = 3 by default; optionally a tighter c_down = 2 for the decrease side, which detected simulated falls 13% to 51% faster for continuous return laws; with prices rounded to a coarse tick the gain depends on the tick size, and at a tick of half a standard deviation it was 12% slower for a 20% fall (E11f) | `robust.py`, `monitor.py` |
| 2 | Bounded statistic X_t = min(r_t², τ²), after an intraday profile estimated from past days only | `robust.py`, `data.py` |
| 3 | Baseline μ₀ = mean of X on a calibration window: plug-in, plug-in with a margin (band), or a confidence sequence | `detectors.py`, `confseq.py` |
| 4 | Betting e-detectors (Shiryaev-Roberts form): wealth M_t = Σ_j Π_{i=j..t} (1 + λ (X_i / μ₀ − 1)) over start times j for increases, with λ in [0, 1), and the same with factors 1 + λ (1 − X_i / μ₀), λ < 1 / (B − 1), B = τ² / μ₀, for decreases; each averaged over a grid of bets λ | `detectors.py` |
| 5 | Alarm when either side's wealth reaches 2A (so the two sides together keep the budget); A = m / α with m the steps per day, which gives at most α chance of a false alarm within the first day after each (re)calibration and at most α expected false alarms per day in the long run, not α on every later day (TH6); then recalibrate | `monitor.py` |
| 6 | Optional forecast-relative mode: divide returns by a GARCH or EWMA variance forecast first | `forecast.py` |

The proofs (Theorem A: average run length ≥ A; Theorem B: P(false alarm within m) ≤ m / A; Corollary C for restarts; Theorem 6.4 for the confidence-sequence baseline) and the exact conditions are in [`docs/math.md`](docs/math.md). [`docs/equation_map.md`](docs/equation_map.md) maps every symbol to the line of code that computes it.

## Data

| Source | Use | Access |
| --- | --- | --- |
| Binance spot 1-minute bars: BTC, ETH, XRP, SOL, DOGE | The paper's five coins (E6) | Public GitHub mirror (exchange APIs are blocked from our cloud machine) |
| Bitstamp BTC/USD 1-minute bars | First real-data run (E5) | Public GitHub mirror |
| Nifty 50 1-minute bars, April 2017 to October 2026 | Indian market application (E7) | Public GitHub mirror |
| India VIX daily | Check of Nifty alarms (E7) | NSE archive server |

Raw data is not committed. `python scripts/download_data.py` fetches everything and writes `data/MANIFEST.csv` with a SHA-256 hash of each file. Details and caveats: [`data/README.md`](data/README.md).

## Install

Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
python -m pytest -q              # the full test suite, no data needed
```

Or install the package alone from GitHub: `pip install "git+https://github.com/ShockXR/Anytime-Valid-Volatility-Monitor"`. Version 0.3.0 is release-ready but not released; see `CHANGELOG.md`. The public API is listed in [`docs/api.md`](docs/api.md) (rebuild it with `python scripts/build_api_docs.py`), and `examples/` holds four short runnable scripts on synthetic data: a quickstart, the confidence-sequence baseline, the monitor restarted every day, and a panel of assets with e-BH.

## Run

```bash
python -m avmonitor simulate --ratio 2                       # synthetic stream, sd doubles at step 3000
python -m avmonitor estimators --source binance --symbol BTCUSDT --day 2025-10-10 --download   # all quarticity estimators on one real day
python -m avmonitor run --source binance --symbol ETHUSDT --download \
    --start 2025-03-01 --end 2025-03-31 --plot eth.png       # real data, alarms printed and plotted
```

`notebooks/demo.ipynb` is an executed walk-through (rebuild it with `python scripts/build_demo_notebook.py`). [`docs/demo/replay.html`](docs/demo/replay.html) is a self-contained page that replays two landmark days of E18 minute by minute (Nifty 50 on 4 June 2024, when the monitor flagged the shock 11 minutes after the open, and BTC on 10 October 2025, when it was recalibrating and blind): the statistic, both e-processes against the threshold, the blind periods, and the alarms of every competitor, checked against E18's saved alarm lists when the page is built (`python scripts/build_replay.py`, about 20 minutes, or seconds with `--cache`). [REAL]

Minimal use from Python:

```python
import numpy as np
from avmonitor.monitor import VolatilityMonitor

rng = np.random.default_rng(0)
r = np.concatenate([rng.standard_normal(3000), 2 * rng.standard_normal(1000)])  # sd doubles at step 3000
mon = VolatilityMonitor(threshold=7800, cal_len=1000, null_band=0.25)          # A = 390 / 0.05; plug-in baseline (approximate)
for e in mon.run(r):
    print(e.index, e.direction)
```

## Reproduce the experiments

Each experiment is one script in `experiments/` that writes CSV files to `results/` and figures to `results/figures/`. Synthetic experiments need no data; real-data experiments need `python scripts/download_data.py` first. Seeds are fixed. `python scripts/reproduce.py` checks the data hashes, runs the tests and reruns the experiments against `results/MANIFEST.csv` (`--tier full` for every experiment, `--from E18` to resume). On a fresh clone the quick tier reproduced every file, and the full tier reproduced every file of E1 to E17 (identical, within a relative tolerance of 1e-9, or equal except run-time columns); E18 onward has not yet been rerun on a fresh clone.

| Experiment | Question | Script | Label |
| --- | --- | --- | --- |
| E1 | Reproduce the paper's simulation tables | `e01_reproduce_tables.py` | SYNTHETIC vs PUBLISHED |
| E2 | Truncation bias, fair comparison of estimators | `e02_truncation_and_fair_comparison.py` | SYNTHETIC |
| E3 | Online (past-only) truncation | `e03_online_truncation.py` | SYNTHETIC |
| E4 | False alarms and delay against benchmarks; band, calibration length, plug-in checks | `e04_detection.py`, `e04c_null_band.py`, `e04d_calibration_length.py`, `e04e_review_checks.py` | SYNTHETIC |
| E5 | Bitstamp BTC/USD run | `e05_real_btc.py` | REAL |
| E6 | Five coins: paper tables, monitor, forecasts, rate-matched competitor | `e06a_five_coins_tables.py`, `e06b_five_coins_monitor.py` | REAL vs PUBLISHED |
| E7 | Nifty 50, India VIX, Union Budget and policy days | `e07_nifty.py`, `e07c_india_vix.py`, `e07d_events.py` | REAL |
| E8 | Confidence-sequence baseline | `e08_cs_baseline.py` | SYNTHETIC and REAL |
| E9 | Stress tests of the guarantee, runtime | `e09_stress_guarantee.py` | SYNTHETIC |
| E10 | Clustering and forecast-relative monitoring | `e10_forecast_relative.py` | SYNTHETIC and REAL |
| E11 | Detection delay theory and bet design; a tighter clip for decreases | `e11_delay_theory.py`, `e11f_split_c.py` | SYNTHETIC and REAL |
| E12 | FOMC and CPI event study on crypto | `e12_crypto_events.py` | REAL |
| E13 | Benchmark against published change-point methods (BOCPD) and a two-window variance-ratio scan | `e13_benchmarks.py offline-null`, `offline-change`, `online-null`, `online-change 0`, `online-change 1`, `real-calib`, `real-crypto`, `real-crypto-h1-2025`, `real-nifty`, `cost`, `figures` | SYNTHETIC and REAL |
| E14 | Running truncation level and monitors without a blind period | `e14_running_tau.py` | SYNTHETIC |
| E15 | Panel of assets with false discovery control | `e15_panel.py synthetic`, `horizon`, `real`, `cs`, `figures` | SYNTHETIC and REAL |
| E16 | Adversarial search for null streams that push the false alarm rate up | `e16_adversarial_null.py a`, `b`, `figures` | SYNTHETIC |
| E17 | One monitor for 1, 5 and 15 minutes (certified average of e-detectors) | `e17_multifrequency.py baselines levels null change real ratematch figures` | SYNTHETIC and REAL |
| E18 | Landmark shocks: pre-registered case studies against published methods | `e18_landmark_events.py run`, `figure`, `context`, `awake` | REAL |
| E19 | The guarantee on real data: permuted and block-bootstrapped calm days, the confidence-sequence baseline on five coins, forecast-relative null, Table S2 with cs | `e19_real_null.py pools p1 p2 p3 p4 blind p1std p2boot figs` | REAL and SYNTHETIC |
| E20 | Re-analysis of the real-data inference: claim classification, Holm, day-block bootstrap, forecast regressions and QLIKE | `e20_inference.py checks bootstrap holm events forecast classify` | REAL |
| E21 | Benchmarks not calibrated on Gaussian data only: classical CUSUM and Shiryaev-Roberts rules on the same statistic with simulated or bootstrap thresholds, Student t BOCPD, CUSUM monitoring with a training sample | `e21_benchmarks.py calib`, `null`, `change 0`, `change 1`, `change-merge`, `e13-calib`, `e13-null`, `e13-change 0`, `e13-change 1`, `e13-change-merge`, `real-calib`, `real-bocpd 250` (also 720, 1440), `real h2`, `real h1`, `real-pooled`, `figures` | SYNTHETIC and REAL |
| E22 | E12's FOMC and CPI test rerun with monitors that carry a proved per-day bound (pre-registered) | `e22_guaranteed_events.py run` | REAL |
| E23 | The same test on fresh releases: FOMC statements and CPI releases of January to October 2026 (pre-registered) | `e23_fresh_releases.py power`, `run` | REAL |
| E24 | A proved guard against blind segments of the confidence-sequence baseline; post hoc reruns; forward pre-registration (data to 1 June 2027) | `e24_blind_guard.py theory`, `selection`, `synthetic --paths=800`, `real`, `e19blind`, `future_gate`, `power`, `figs` (`future` after 1 June 2027) | SYNTHETIC and REAL |
| E25 | An anytime-valid version of the forward release test: an e-process over releases with disjoint four-week strata (pre-registered, one amendment before any fresh data); gates, power, post hoc illustration | `e25_sequential_events.py theory`, `strata`, `gates`, `power`, `posthoc`, `figs` (`status` from 26 November 2026) | SYNTHETIC and REAL |
| TH4 | False alarm bounds when the conditional mean drifts (approximate null), checked against E16 | `th4_robustness.py sup constant switching garch plugin ticks counterexamples figures` | proof and SYNTHETIC |
| TH5 | Limits for every valid procedure (pathwise floor and mean delay bound), the plug-in baseline when n, m and A grow together, the decrease bet at its boundary | `th5_limits.py collapse floors asymptotic boundary plugin_bound plugin_sim holdbold figures` | proof and SYNTHETIC |
| TH6 | What Theorem B bounds after the first day; restart at every day boundary and the e-CUSUM; Theorem D with a tolerance hypothesis | `th6_per_period.py laws hazard390 hazard1440 family adversary delay renewal tolerance figures` | proof and SYNTHETIC |

Then rebuild the report from its parts:

```bash
python scripts/build_report.py   # writes docs/experiment_report.md from docs/report_parts/
python experiments/make_figures.py
```

Most scripts take part names as arguments (for example `python experiments/e10_forecast_relative.py a b`) so that one part can be rerun alone. The largest simulations take tens of minutes on one core.

## How to read the output

* **An alarm** means the running evidence against "volatility is at its calibrated level" has reached A. The direction says whether volatility moved up or down. The guarantee is about false alarms: with budget 0.05 per day, a monitor on a stream with no change has at most a 5% chance of a false alarm in the first day after each (re)calibration, and at most 0.05 expected false alarms per day in the long run. On a particular later day the chance can be higher (TH6: about 5% for i.i.d. returns, 28% on a day after nine calm days), unless the statistics restart every day with a known baseline (TH6; with an estimated baseline see the caveat above). This is proved for the true baseline and the confidence-sequence baseline (Corollary C and Theorem 6.4 in `docs/math.md`); for the default plug-in baseline it holds only approximately (E4e), and with short calibration windows it can be well above nominal (11% to 19% instead of 5% with 500 calibration steps in E13's offline design).
* **The guarantee is conditional.** It needs the clipped squared returns to have the calibrated conditional mean while nothing changes. Strong volatility clustering and price-level noise violate this; the forecast-relative mode (when the forecast model is right) and non-overlapping blocks restore it (E10, ALT1).
* **Delay** is measured in steps (minutes). It grows like ln A divided by an information rate, which we compute in closed form for the e-detector's own bets (E11) and by simulation for forecast-relative monitoring (E10).
* **Real-data hit rates** compare the next day's realised volatility with a forecast. They are compared with random times and with a rate-matched rule, never read alone.

## Known limitations

The full list, with what each one affects, is in [`docs/limitations.md`](docs/limitations.md). The main ones:

* The guarantee is proved for a known baseline, for the confidence-sequence baseline and for the rank monitor of E14. The default plug-in baseline is only approximately valid within a day and fails over several days (E4e, E10, E14).
* Volatility clustering is itself a change in the conditional variance; without a forecast the monitor alarms on it. With a forecast it is valid but much slower to see level shifts (E10).
* With the default protocol the monitor is blind for about 70% of the time on real data while it recalibrates after alarms (E12). The rank monitor of E14 (`ConformalMonitor` and `AwakeRunningTauMonitor` in `running.py`) removes the blind period with a proof, but without a null band it raises 3.6 to 4.3 alarms per crypto day.
* The two-window variance-ratio scan of E13 is our construction inspired by Bibinger, Jirak and Vetter (2017), not their procedure: the paper's full text could not be opened from our environment. Likewise the CUSUM monitoring rules of E21 are of the type of Horváth, Hušková, Kokoszka and Steinebach (2004), not their exact procedure, and a β-divergence BOCPD was not implemented. Their conclusions use critical values calibrated by simulation or matched alarm counts, which do not depend on the published constants.
* The Binance and Nifty data come from public mirrors, not from the exchanges. The paper does not name its crypto data source, so our five-coin data may differ from the paper's.
* The bound is about four times conservative for ordinary returns (E9) and about two times at one day (1.4 times over five days) against the worst null stream we found (E16); under that stream the rest is wealth that the default bets leave on paths that never alarm, and whether a worse stream exists is open.

## Repository layout

```
.
├── PITCH_BAPAT.md              short brief: ten questions answered
├── RESEARCH_STRATEGY_BAPAT.md  questions, hypotheses, scores, risks, roadmap
├── ITERATION_LOG.md            iterations in a fixed format: problem, hypothesis, change, result
├── CHANGELOG.md                versions 0.1.0 to 0.3.0 (0.3.0 unreleased)
├── LICENSE                     MIT
├── alternatives/               eight other directions (ranked in ALTERNATIVES.md), each with its own README and RESULTS
├── data/                       raw data (git-ignored), MANIFEST.csv with hashes, README
├── docs/
│   ├── paper_review.md         the paper in plain words
│   ├── estimator_issues.md     checks of the estimator, with evidence
│   ├── math.md                 derivations and proofs
│   ├── theory_delay.md         detection delay and bet design
│   ├── theory_panel.md         many assets at once: e-BH and false discovery control
│   ├── theory_running_tau.md   running truncation level, rank detector, no blind period
│   ├── theory_multifreq.md     one monitor for several sampling frequencies (certified average)
│   ├── theory_robustness.md    false alarm bounds when the null holds only approximately (TH4)
│   ├── theory_limits.md        limits for every valid procedure; plug-in regime in n, m and A (TH5)
│   ├── theory_per_period.md    false alarms after the first day; bounds for every day (TH6)
│   ├── theory_blind_guard.md   a proved guard against blind confidence-sequence segments (E24)
│   ├── theory_sequential_events.md  an anytime-valid release test: e-values over releases (E25)
│   ├── equation_map.md         symbol to code map
│   ├── api.md                  public API reference (generated by scripts/build_api_docs.py)
│   ├── literature_review.md    verified references and the gap
│   ├── limitations.md          what the results do not show
│   ├── experiment_report.md    all results (generated from report_parts/)
│   ├── demo/replay.html        minute-by-minute replay of two landmark days (E18)
│   ├── report_parts/           one file per experiment
│   ├── iterations/             per-workstream iteration notes
│   └── review/                 independent reviews and responses
├── examples/                   four short runnable scripts on synthetic data
├── experiments/                one script per experiment
├── paper/                      short paper draft (Markdown source, LaTeX, verified BibTeX)
├── notebooks/demo.ipynb        executed demonstration
├── results/                    CSV outputs and figures
├── scripts/                    data download, report build, notebook build
├── src/avmonitor/
│   ├── robust.py               quarticity estimators (as printed and normalised), MAD truncation, online MAD
│   ├── detectors.py            betting e-detectors, Gaussian LR, Page CUSUM, repeated z-test
│   ├── confseq.py              confidence-sequence baseline
│   ├── monitor.py              streaming calibrate, monitor, alarm, recalibrate loop
│   ├── forecast.py             EWMA and GARCH variance forecasts for forecast-relative monitoring
│   ├── competitors.py          published change-point methods used as benchmarks
│   ├── benchmarks_robust.py    classical rules with simulated or bootstrap thresholds, Student t BOCPD, training-sample CUSUM (E21)
│   ├── running.py              running truncation level, rank (conformal) monitor, monitors with no blind period
│   ├── multifreq.py            certified average of e-detectors at several sampling frequencies
│   ├── robustness.py           drift inflation factors and pathwise bounds for an approximate null (TH4)
│   ├── limits.py               pathwise floors, delay lower bounds, plug-in bounds (TH5)
│   ├── periodic.py             monitor restarted every day, e-CUSUM, per-day bounds (TH6)
│   ├── guard.py                confidence-sequence monitors with a guard against blind segments (E24)
│   ├── design.py               bet grid design from the delay theory
│   ├── panel.py                panel monitor with e-BH across assets
│   ├── simulate.py             paper design, stochastic volatility with jumps and noise, GARCH, change points
│   ├── data.py                 loaders, predictable intraday profile, frozen-price detection
│   ├── strategy.py             interface for the risk switch (implemented in alternatives/risk_switch)
│   └── cli.py                  command line
└── tests/
```

## Roadmap

- [x] Robust kernel estimators, as printed and normalised, checked against the paper's simulation tables
- [x] Online (past-only) truncation
- [x] E-detector with proofs, compared with Gaussian, CUSUM and repeated-test benchmarks
- [x] Real data: Bitstamp BTC, the five Binance coins, Nifty 50, event studies
- [x] Confidence-sequence baseline with a proof
- [x] Forecast-relative monitoring under clustering, with the signal cap
- [x] Risk switch backtest (alternatives/risk_switch)
- [x] Monitors with no blind period after alarms, with a proof (E14)
- [x] Short paper draft (`paper/`)
- [ ] Discussion with Prof. Bapat of the open questions in [`docs/estimator_issues.md`](docs/estimator_issues.md)

## References

The full, verified list is in [`docs/literature_review.md`](docs/literature_review.md).

- Bahadur, R. R. and Savage, L. J. (1956). The nonexistence of certain statistical procedures in nonparametric problems. *Annals of Mathematical Statistics*, 27(4), 1115 to 1122.
- Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482. doi:10.1016/j.iref.2026.105482
- Barndorff-Nielsen, O. E. and Shephard, N. (2004). Power and bipower variation with stochastic volatility and jumps. *Journal of Financial Econometrics*, 2(1), 1 to 37.
- Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. *Journal of Econometrics*, 31(3), 307 to 327. doi:10.1016/0304-4076(86)90063-1
- Lai, T. L. (1998). Information bounds and quick detection of parameter changes in stochastic systems. *IEEE Transactions on Information Theory*, 44(7), 2917 to 2929.
- Lamoureux, C. G. and Lastrapes, W. D. (1990). Persistence in variance, structural change, and the GARCH model. *Journal of Business and Economic Statistics*, 8(2), 225 to 234.
- Lee, S. S. and Mykland, P. A. (2008). Jumps in financial markets: a new nonparametric test and jump dynamics. *Review of Financial Studies*, 21(6), 2535 to 2563.
- Mancini, C. (2009). Non-parametric threshold estimation for models with stochastic diffusion coefficient and jumps. *Scandinavian Journal of Statistics*, 36, 270 to 296.
- Mikosch, T. and Starica, C. (2004). Nonstationarities in financial time series, the long-range dependence, and the IGARCH effects. *Review of Economics and Statistics*, 86(1), 378 to 390.
- Page, E. S. (1954). Continuous inspection schemes. *Biometrika*, 41(1/2), 100 to 115.
- Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science*, 2(2), 229 to 260. doi:10.51387/23-NEJSDS51
- Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *Journal of the Royal Statistical Society Series B*, 86(1), 1 to 27.

## Author

Shreyan Thanki ([@ShockXR](https://github.com/ShockXR))
