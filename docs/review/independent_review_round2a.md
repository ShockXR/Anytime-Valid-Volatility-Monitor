# Independent review, round 2a

> *Note added 3 November 2026: the "independent reviewer" or "referee" below was a separate AI agent (Claude, by Anthropic) instructed to act as a sceptical referee. It was not a human expert, and the review has not been confirmed by a human.*

**Reviewer role.** An adversarial reader: a statistics journal referee and a code auditor. The goal is that every claim shown to Prof. Bapat survives an expert reading.

**Scope.**

1. `docs/math.md` sections 5 and 6 (Theorems A and B, Corollary C, Lemma 6.1, Theorem 6.2, Theorem 6.4, Corollaries 6.5 and 6.6), and whether `detectors.py`, `confseq.py` and `monitor.py` meet the conditions the proofs need.
2. `docs/report_parts/E8.md`, `E9.md`, `E10.md`, `E12.md`, their scripts and the CSV files they cite.
3. `alternatives/robust_iq_estimator/` (README, RESULTS, `riq.py`, results) and `docs/report_parts/ALT3.md`.
4. `README.md` and `RESEARCH_STRATEGY_BAPAT.md`: every number traced to a CSV or a report part, and every claim of novelty or improvement checked.

**Read first.** `docs/review/independent_review_round1.md`. Points fixed since round 1 are not repeated. Where a round 1 point comes back in a file that round 1 did not cover, I say so.

**What I ran.** Every number in the four report parts and in ALT3 was recomputed from its CSV. Small reruns, all with output in `/tmp/claude-0/review2a/` and nothing written to `results/`:

* `check_e10c_init.py`, `check_e10c_rerun.py`: the E10 real-data part (Q10c) with the committed script and with a corrected start value.
* `check_e10d_pooled.py`: the E10 signal-cap simulation at persistence 0.999, same seed as the script.
* `check_riq_theory.py`, `check_prop1_se.py`: the ALT3 closed forms (truncation constants, AVARs, jump residuals, cross-validation) and a Monte Carlo check of ALT3 Proposition 1.
* `python -m pytest -q`: **220 passed in 46.7 s.**

No file in the repository was changed except this one. No git command that changes the tree was run.

**Severity.** *Critical*: wrong mathematics, a wrong number, or a conclusion the evidence does not support. *Major*: an overclaim, a missing condition or a misleading comparison. *Minor*: wording or clarity.

---

## 1. Summary verdict

**The mathematics holds.** I re-derived Lemma 1, Theorems A and B, Corollary C, Lemma 6.1, Theorem 6.2, Corollary 6.3, Theorem 6.4 (all three parts) and Corollaries 6.5 and 6.6, and found no error. The code meets the conditions the proofs need:

* τ is fixed before the confidence sequence starts;
* bets are fixed at the end of calibration;
* a sequence that has started is never restarted;
* missing returns are skipped;
* every rule that ends calibration uses only data seen so far.

The E10 signal cap ρ(k) = k(1 - β)/(1 - α - β + αk), its limit (1 - β)/α, the Kullback-Leibler (KL) rate and the Lai delay figures are also correct (section 5 shows the derivations).

**The problems are in the reporting.** Eight numbers or claims are wrong (Critical). Most are in the robust IQ alternative (ALT3) and in headline sentences of the README and E12. None of them overturns a main conclusion of the monitor project, but an expert reader would notice several at once. Three of them (C3, C8, and the wording covered in M5 to M7) make the author's estimator look worse than the project's own tables support. That matters most for a meeting with the author.

**Most important fix before the meeting.** Make ALT3 agree with its own Q1 table, and separate two statements that are currently merged:

* "biased by a known constant under constant volatility" (removed by the κ₄ correction);
* "inconsistent when volatility moves within the day".

**Reproducibility.** The E10 real-data part cannot be reproduced from the committed script (M1), although the CSV itself is correct.

---

## 2. Critical findings

### C1. E12: "alarms within the hour" is wrong for ETH, SOL and the fresh design

* **Where.** `docs/report_parts/E12.md:98`.
* **Says.** "Whenever the robust monitor is awake at an FOMC statement it alarms within the hour (9/9, 8/8, 9/9 across coins in the streaming design; 14/16, 14/16, 15/16 with a fresh calibration)".
* **Evidence.** These counts are the "up alarm within **6 hours**" outcome (Q12a and line 31; `results/e12_delay.csv`, column `with_up_in_6h`). The within-60-minute counts in the same CSV (`with_up_in_first_60min`) are:
  * streaming: 9 of 9 (BTC), 6 of 8 (ETH, latest alarm after 123 minutes), 8 of 9 (SOL, latest after 66 minutes);
  * fresh calibration: 14, 13 and 11 of 16.

  E12's own Q12b table (lines 48 to 55) shows the 14, 13 and 11.
* **Fix.** "... it alarms within six hours (9/9, 8/8, 9/9; 14/16, 14/16, 15/16 with a fresh calibration), usually within the first hour (streaming 9/9, 6/8, 8/9; fresh 14, 13 and 11 of 16)".

### C2. ALT3: "minus 14% even after the kappa correction" is minus 11%

* **Where.** `alternatives/robust_iq_estimator/RESULTS.md:8`.
* **Evidence.** `results/q2_rates.csv`, scenario S5 at n = 23,400: RK (no correction) has bias -0.1438; RK-k (kappa corrected) has -0.1117.
* **Fix.** "(minus 11% at n = 23,400 after the kappa correction, minus 14% without it)".
* **Related.** `docs/report_parts/ALT3.md:57` quotes minus 14% without naming the correction. It should quote the corrected figure, because the uncorrected one includes the constant minus 3.6% that the correction removes: 0.964 × 0.888 = 0.856 (see M5).

### C3. ALT3: the efficiency of the author's estimator is understated

* **Where.**
  * `RESULTS.md:7`: "lowers the asymptotic variance from 10.7 (RQ, and every RQ-type estimator, including Bapat's averaged kernel)".
  * `docs/report_parts/ALT3.md:56`: "raises efficiency from 0.74 (RQ, and Bapat's averaged kernel, which is RQ in disguise) to 0.95".
  * The same framing appears in `ALT3.md:7` and `RESULTS.md:38` ("RQ in disguise").
* **Evidence.**
  * The identity in `docs/math.md` section 4 says the averaged kernel equals a rescaled **truncated** RQ plus an edge term. At n = 390, the 2h = 60 edge returns (15% of the sample) get partial weight.
  * Truncation changes the variance. For a fourth power clipped at three standard deviations, Var(min(Z⁴, 81)) / E[min(Z⁴, 81)]² = 8.19. I got this by numerical integration; `riq.avar_truncated_rq(3, "clip")` gives 8.186, and a Monte Carlo with 2 × 10⁷ draws gives 8.188. Untruncated RQ gives 96/9 = 10.67.
  * ALT3's own Q1 table (`ALT3.md:13`, `results/q1_avar.csv`) lists RK-k with a simulated AVAR of 9.26 and efficiency 0.86. The ALT3 README comparison table (`README.md:129`) says "about 8.6 to 9.3".
  * So under constant volatility the fair comparison for the author's kappa-corrected kernel is:
    * 0.86 (feasible, estimated MAD, simulated);
    * 0.98 (oracle cut-off, closed form);
    * against 0.95 for PI-LU (simulated AVAR 8.46; closed form 8.20).
  * Squaring a local variance gains a lot over untruncated RQ (10.67 to 8.20) but little over the author's corrected estimator, and nothing over its oracle version.
* **Fix.**
  * `RESULTS.md:7`: "lowers the asymptotic variance from 10.7 (untruncated RQ) to 8(k + 3)/(k + 2) ... The author's averaged kernel is a truncated RQ (math.md section 4) and already gains from truncation: 8.19 with an oracle cut-off and 9.26 (efficiency 0.86) with the estimated MAD and the kappa correction, against 8.46 (0.95) for PI-LU."
  * Make the same change in `ALT3.md:56`.
  * Replace "RQ in disguise" with "a truncated RQ with edge weights".

### C4. ALT3 README: wrong count of wrongly truncated returns

* **Where.** `alternatives/robust_iq_estimator/README.md:72`.
* **Says.** n P(|Z| > a_n) ~ 1/sqrt(π log n) → 0 "(0.22, 0.05 and 0.003 for n = 390, 1440, 23,400)".
* **Evidence.** With a_n = sqrt(2 log n):
  * the exact normal tail gives n P(|Z| > a_n) = 0.215, 0.197 and 0.170;
  * the approximation 1/sqrt(π log n) gives 0.231, 0.209 and 0.178 (section 5, D5).

  The count falls only like 1/sqrt(log n). `RESULTS.md:49` states it correctly ("fewer than one diffusive return per day at every n").
* **Fix.** "(0.22, 0.20 and 0.17 for n = 390, 1440 and 23,400; the count falls only like 1/sqrt(log n), but the relative bias it causes is O(sqrt(log n)/n))".

### C5. README: "a doubling of volatility looks like a 9% rise"

* **Where.** `README.md:27`.
* **Evidence.**
  * The 9% is ρ(2) = 1.091 at α = 0.05, β = 0.94 (`E10.md:79`, `results/e10d_signal_cap.csv`). That is a doubling of the long-run **variance**.
  * A doubling of volatility (the standard deviation) is k = 4: ρ(4) = 1.143 as a ratio of means, and 1.187 as the mean ratio the monitor sees.
  * "Passes at most (1 - β)/α of any level shift" misdescribes the cap. (1 - β)/α caps the post-break ratio of the mean squared return to the mean forecast; it is not a fraction of the shift.
* **Fix.** "A doubling of the long-run variance looks like a 9% rise in the ratio of squared returns to the forecast at persistence 0.99, and no break, however large, can push that ratio above (1 - β)/α = 1.2."

### C6. ALT3: "underestimates IQ by 26% to 59% at one minute"

* **Where.** `RESULTS.md:12`.
* **Evidence.** `results/q6_fastvol_estimators.csv` at n = 1440 (one minute), for the estimators that truncate against a scale estimated over hours (TRQ-M, TRQ-U, RK, RK-k, PI-G3, PI-L3, PI-LU, PI-LUc, KPI-L3):
  * without jumps: relative bias -0.346 to -0.582, median error -0.353 to -0.587;
  * with medium jumps: relative bias -0.311 to -0.558.

  The 26% is PI-LU at n = 23,400 (-0.265), not at one minute.
* **Fix.** "by about 35% to 59% at one minute (31% to 56% with jumps)".

### C7. E8 and README: "1.8 to 7.2 times the delay"

* **Where.** `docs/report_parts/E8.md:70`; `README.md:26`.
* **Evidence.** In `results/e08b_change.csv`, the ratio of the cs median delay to the plug-in median delay at cal_len 390, over all 18 cells, runs from 1.82 to 7.56. The maximum is jumps at ×1.25 (695.5 / 92 steps). E8 says "All 54 cells are in e08b_change.csv", but its printed table leaves out that row, and the range was read off the printed rows.
* **Fix.** "1.8 to 7.6 times" in both places.

### C8. ALT3: "PI-LU has the lowest RMSE of all 15 estimators in the jump scenarios at n = 390 and 1440"

* **Where.** `RESULTS.md:9`; `docs/report_parts/ALT3.md:57` ("is the best estimator in the jump scenarios at n = 390 and 1440").
* **Evidence.** `results/q2_q4_sim_summary.csv` has 8 jump cells (S2, S3, S4 and S6 at n = 390 and 1440). PI-LU is lowest in 6 of them. The exceptions:
  * S2 at n = 390: PI-LU is **third** (RK 0.156, PI-G3 0.158, PI-LU 0.165);
  * S6 at n = 390: PI-LU is **fifth** (PI-G3 0.191, RK-k 0.197, RK 0.198, KPI-L3 0.215, PI-LU 0.217).

  RK and RK-k are the author's estimator. The summary's own range, "0.163 to 0.217 at n = 390", includes the S6 value where PI-LU is fifth.

  `RESULTS.md:204` says "best or within 0.01 of the best in 10 of the 12 cells". By my count it is 11 of 12: only S6 at n = 390 is further off, by 0.026.
* **Fix.** "PI-LU has the lowest RMSE in 6 of the 8 jump cells at n = 390 and 1440. In S2 and S6 at n = 390, the author's kernel and PI-G3 are slightly better, because their downward truncation bias offsets the upward jump bias." Change "10 of the 12" to "11 of the 12".

---

## 3. Major findings

### M1. E10 real-data part (Q10c) cannot be reproduced from the committed script

* **Where.** `experiments/e10_forecast_relative.py:211` and `:219`; `results/e10c_real_forecast_relative.csv`; `docs/report_parts/E10.md:104` to `121`.
* **Evidence.**
  * Both lines start the forecast with `init=float(np.nanmean(u_all[:DAY] ** 2))`.
  * `data.intraday_scale` now defaults to `start="past"` (changed after the CSV was written). So the first 10 days of `u_all` are NaN, and `init` is NaN.
  * `forecast.ewma_variance` uses `max(float(init), floor)`, which returns NaN, so every forecast value is NaN.
* **What the rerun gives.** Rerunning BTC with the committed script (`/tmp/claude-0/review2a/check_e10c_rerun.py`):
  * "deseasonalised only": 127 alarms, but 0 kept for validation;
  * "divided by EWMA, half-life 360": 0 alarms, because its input is divided by the NaN forecast and the monitor never calibrates.

  With `init` taken from the first day with finite values, the rerun reproduces the CSV row exactly: 127 alarms, 72 up with 66.67% above the forecast, 54 down with 72.22% below.
* **Why it matters.** The CSV is right, but a rerun from a clean clone silently gives empty results. That contradicts `RESEARCH_STRATEGY_BAPAT.md:108` ("A researcher can rerun every table from a clean clone").
* **Fix.**
  * Start from the first finite day: `np.nanmean(u_all[i0:i0 + DAY] ** 2)`, where i0 is the first index with a finite value.
  * Make `ewma_variance` raise an error on a non-finite `init`, and add a test.
  * Regenerate e10c and confirm that the numbers are unchanged.

### M2. E12: pooled control rates and the pooled Fisher test mix weekdays

* **Where.** `docs/report_parts/E12.md:16` (control definition), the Q12a table (lines 22 to 27), the Q12c table and text (lines 63 to 72), and line 98 ("on ordinary Wednesdays").
* **How the controls are built.** Controls are matched to each event by weekday. But the "control rate" columns and the pooled Fisher test use the **union** of all events' control windows:
  * FOMC: 65 Wednesdays and 69 Thursdays, although only one statement (7 November 2024) fell on a Thursday;
  * CPI: 71 Tuesday, 65 Wednesday, 69 Thursday and 94 Friday windows, and only one release (24 October 2025) fell on a Friday.
* **Why it biases the down outcome.** A window anchored on a Thursday afternoon looks 6 to 48 hours ahead into the weekend, when crypto volatility falls. For BTC with a fresh calibration, a down alarm followed in 94% of Thursday 14:00 controls but in 46% of Wednesday 14:00 controls.
* **The matched rates already exist.** `results/e12_delay.csv`, column `control_down_6h_to_48h_rate`, has the event-weighted matched rates. For FOMC (fresh calibration):

  | | BTC | ETH | SOL |
  |---|---|---|---|
  | matched control rate | 49% | 43% | 42% |
  | union rate printed in Q12c | 71% | 66% | 65% |
  | event rate | 63% | 56% | 56% |

* **Consequences.**
  1. Q12c's explanation "because the base rate is high" is largely a weekday artefact. Against matched controls, FOMC events are followed by down alarms more often than controls are. The stratified permutation test, which matches correctly, still gives p = 0.17 to 0.21, so "not significant" stands, but the stated reason is wrong.
  2. The pooled Fisher p-values for the down outcome do not test the matched design.
  3. For the "up within 6 hours" outcome the effect is small: matched FOMC fresh rates are 29%, 32% and 32%, against union rates of 31%, 31% and 29%. So Q12a and the README's "about 30%" are fine. Line 98 should still say "matched weekdays", not "ordinary Wednesdays".
* **Fix.**
  * Print the event-weighted matched control rate instead of the union rate.
  * Drop the pooled Fisher test, or stratify it (Cochran-Mantel-Haenszel).
  * Rewrite the explanation in Q12c.

### M3. E12: the Gaussian against robust comparison is confounded by the band

* **Where.** `E12.md:44` ("Its rate ratio is 1.2 to 1.6, against 2.2 to 3.2 for the robust detector"), `:102`; caveat at `:113`.
* **Evidence.**
  * The Gaussian LR monitor runs with band 0 and the robust monitor with band 0.25. `GaussianLRMonitor` ignores `null_band` (`monitor.py:327` to `367`).
  * The band alone moves the alarm rate a lot: 156 plug-in alarms at band 0 against 126 at band 0.25 (E8 Q8e).
  * So the Gaussian control rate of 64% to 77% is partly the narrower null and partly the missing truncation, and the rate-ratio contrast cannot be credited to truncation.
  * The Q12d "jump only" column does not separate the detectors at releases either: robust 19%, 40%, 35% against Gaussian 33%, 41%, 30%.
  * The caveat exists, but the Reading draws the conclusion anyway.
* **Fix.** Either run the robust monitor at band 0 (or give the Gaussian monitor a matching band) for Q12a and Q12d, or limit the Reading to what is supported. What is supported: the Gaussian alarm is triggered by one large minute (median largest recent return 8 to 9 scale units against about 4; delays of 0 to 2 minutes). Then drop the rate-ratio comparison.

### M4. E8 Q8e repeats round 1 M3

* **Where.** `docs/report_parts/E8.md:149`: "Against a 50% chance level (E5c control), all variants are informative."
* **Evidence.**
  * Round 1 M3 showed that a naive trailing-window rule scores 82% and 73% on the same check.
  * E6e (`E6.md:169` to `176`) shows that a rate-matched ratio rule is about as precise as the detector.
  * 50% is the random-time rate, not the relevant competitor.
* **Fix.** "All variants beat random times (about 50%). A rate-matched ratio rule does about as well (E6e), so this check shows that the alarms are not noise, not that they add precision."

### M5. ALT3: "inconsistent" must separate the two settings

* **Where.** `RESULTS.md:8`, `:111`, `:315`; `docs/report_parts/ALT3.md:57`; `alternatives/robust_iq_estimator/README.md:64` and `:70`.
* **Evidence.** There are three different situations:
  * **The paper's simulation design.** σ is constant per step (E1, `E01_E05.md:7`). With one daily MAD and c = 3, the estimator converges to κ₄(a) times IQ, where κ₄(a) = E[min(Z⁴, a⁴)]/3. That gives 0.964 with the scaled MAD (a = 3) and 0.705 with the raw MAD (a = 2.02) (my integration; `docs/estimator_issues.md`, table at line 100). This is inconsistent, but by a **known constant**, and dividing by κ₄ removes it exactly.
  * **The paper's Theorem 1.** The theorem allows σ to change over time, but its proof needs τ to shrink more slowly than the returns. A fixed multiple of the MAD does not (`docs/estimator_issues.md:105`). So the MAD rule falls **outside the theorem's conditions**. The theorem is not wrong; its conditions are not met.
  * **Volatility that moves within the day.** The kappa-corrected estimator keeps a bias that depends on how much σ moves. At n = 23,400 it is minus 0.6% to minus 1.1% in S1 (vol of vol 0.5, no intraday pattern; `RESULTS.md:111`) and minus 11% in S5 (vol of vol 1 with an intraday pattern). Only S5 shows a large bias.
* **What goes wrong in the text.** ALT3's "minus 14%" mixes the constant part and the moving-volatility part (C2). The README sentence "which is wrong whenever volatility moves within the day" (`README.md:64`) reads as a judgement on the author rather than a description.
* **Fix.** Use one wording everywhere:

  > "With a fixed c and one daily MAD, the estimator converges to κ₄(c) times IQ under constant volatility (minus 3.6% with the scaled MAD, minus 29.5% with the raw MAD), which the kappa correction removes. When volatility moves within the day, the corrected estimator keeps a bias that does not shrink with n: about minus 1% in our moderate design S1 and minus 11% in the strong design S5. A local scale removes it."

  Replace "which is wrong whenever volatility moves within the day" with "which does not adapt when volatility moves within the day".

### M6. ALT3: "the published estimator is ... about three times too low on real data"

* **Where.**
  * `RESULTS.md:329`;
  * `alternatives/robust_iq_estimator/README.md:117` ("too low by a factor 3 to 7 on real data");
  * `RESULTS.md:12` and `:296` ("Bapat's included", "costs Bapat's estimator 56% to 59% of IQ");
  * the same idea in `PITCH_BAPAT.md:42`.
* **Evidence.**
  * **The truth is unknown on real data.** What is measured is a ratio to MedRQ: RK / MedRQ = 0.33 to 0.42 (`RESULTS.md:275`), so RK is 2.4 to 3.0 times lower. In the burst model, MedRQ is itself 15% (bias) to 20% (median error) below the truth at n = 1440 (`q6_fastvol_estimators.csv`). So "too low" depends on the model.
  * **RK is not the published estimator.** RK is this repository's normalised kernel with the scaled MAD. The published estimator has unnormalised weights, so it is not on the IQ scale at all (`docs/estimator_issues.md`, issue 1), and E1 and E5a suggest it used the raw MAD.
  * **"3 to 7" is a Bitstamp figure.** It comes from E5. On Binance the figure is 2.3 to 3.0 (`RESULTS.md:241`).
* **Fix.** "Our normalised version of the robust kernel (scaled MAD, c = 3) is 0.33 to 0.42 times MedRQ on Binance one-minute data. If the burst model is right, it misses about 56% to 59% of IQ." Remove "the published estimator is inconsistent and about three times too low".

### M7. "We do not claim to have improved the estimator" conflicts with ALT3

* **Where.** `README.md:32` and `RESEARCH_STRATEGY_BAPAT.md:71`, against:
  * `alternatives/robust_iq_estimator/README.md:118` ("fixes the four issues of docs/estimator_issues.md");
  * `RESULTS.md:329` ("shows exactly when it is consistent and how to make it efficient");
  * `ALT3.md:56`.
* **Evidence.** ALT3 proposes estimators that it ranks above the author's, and the efficiency claim behind that ranking is overstated (C3, C8). A reader of `README.md:32` who then opens ALT3 sees a contradiction.
* **Fix.**
  * `README.md:32`: "The main project does not claim to improve the estimator. ALT3 explores an alternative IQ estimator. In simulation it is somewhat more efficient than the kappa-corrected robust kernel under constant volatility (0.95 against 0.86) and stays consistent when volatility moves within the day; it is a proposal for discussion."
  * ALT3 `README.md:118`: replace "fixes the four issues" with "addresses the four issues, in simulation".

### M8. E10 Q10d: the "ratio of means, simulated" column is biased and appears to break the cap

* **Where.** `docs/report_parts/E10.md:72` to `83`; `experiments/e10_forecast_relative.py:312`.
* **Evidence.**
  * The script averages per-path ratios mean(s1)/mean(s0) over 200 paths. A mean of ratios is a biased estimate of a ratio of means.
  * At persistence 0.999 (half-life 693 steps, 18,000 steps per path) the per-path means are noisy, so the bias is visible. The column shows 1.034 (k = 10) and 1.038 (k = 1000), above the proved cap of 1.0333 printed in the same table.
  * My rerun with the same seed (`/tmp/claude-0/review2a/check_e10d_pooled.py`) reproduces the per-path averages exactly. The pooled ratio, sum of mean(s1) divided by sum of mean(s0), agrees with the formula:

    | k | per-path average (as in the CSV) | pooled ratio | formula ρ(k) |
    |---|---|---|---|
    | 1.5 | 1.0122 | 1.0107 | 1.0109 |
    | 2 | 1.0187 | 1.0163 | 1.0164 |
    | 4 | 1.0291 | 1.0261 | 1.0248 |
    | 10 | 1.0336 | 1.0300 | 1.0299 |
    | 1000 | 1.0377 | 1.0338 | 1.0333 |

  * The "standard errors below 0.002" (line 72) belong to the mean-of-ratio column (`mean_of_ratio_se`), not to this one.
* **Fix.** Report the pooled ratio with its own standard error (bootstrap over paths), and say which column the standard errors refer to.

### M9. E10 and the top-level documents: "restores the guarantee"

* **Where.**
  * `E10.md:24`: "A well specified forecast restores the guarantee";
  * `E10.md:59`: "At one day the plug-in error costs little (about 5% against the 5% budget)";
  * `README.md:125`;
  * `RESEARCH_STRATEGY_BAPAT.md:27` (Q6) and `:40` (H6, "validity restored").
* **Evidence.**
  * Q10a uses the plug-in baseline, which the theorem does not cover. The fitted GARCH rates are 5.8%, 3.0% and 6.2% over 600 paths (6.2% has a Wilson interval of about 4.5% to 8.4%). The oracle rates are 5.3%, 4.7% and 5.3%.
  * These rates sit near the 5% bound only because of the slack in Ville's inequality. With the true baseline, E9 gives 1% to 1.4%, so the plug-in roughly quadruples the rate.
  * `math.md` 5.6(1) already says so ("the informative comparison is with the rate under the true baseline, not with the bound"), following round 1 M4. `E10.md:59` goes back to comparing with the bound.
  * The guarantee also comes back only if the forecast model is right. E10 itself shows that the EWMA forecast with half-life 360 fails (14% to 83% false alarms at band 0.25).
* **Fix.**
  * `E10.md:24`: "With the fitted GARCH forecast, the false alarm rate falls to the oracle's level (3% to 6% with a plug-in baseline). The proof needs the true or confidence-sequence baseline and a correct forecast model."
  * `E10.md:59`: "At one day the plug-in error already raises the rate from about 1.4% (true baseline, E9) to about 5%."
  * README and strategy: "restores it when the forecast model is right".

### M10. README: the FOMC headline quotes the secondary design

* **Where.** `README.md:29`.
* **Evidence.**
  * "14 to 15 of 16 FOMC statements are flagged within six hours" comes from the fresh-calibration design. E12 lists that design as secondary and says it was added after seeing the streaming results (`E12.md:14`, `:18`, `:111`).
  * The design fixed in advance as primary (streaming) flags 9, 8 and 9 of 16 (BTC, ETH, SOL), because the monitor was recalibrating at about half of the releases. At the releases where it was monitoring, it flagged 9 of 9, 8 of 8 and 9 of 9 (`E12.md:31`).
* **Fix.** "With a monitor calibrated on the day before each statement, 14 to 15 of 16 FOMC statements are flagged within six hours, against about 30% of matched control times. The streaming monitor flags 8 to 9 of 16: every statement at which it was monitoring (E12)."

### M11. Strategy: section 6 of math.md was not reviewed in round 1

* **Where.** `RESEARCH_STRATEGY_BAPAT.md:104`.
* **Evidence.** `docs/review/independent_review_round1.md:5`: "I did not review the confidence-sequence path, `confseq.py`, math.md section 6, or E6 to E8 and E10."
* **Fix.** "Section 5: review round 1. Section 6: review round 2a, no errors found."

### M12. Strategy: "the only valid one under jumps and fat tails"

* **Where.** `RESEARCH_STRATEGY_BAPAT.md:24`.
* **Evidence.**
  * In E4 (`E01_E05.md:84`), the e-detector on raw r² also stays within its bound for increases under jumps and Student t5 (2.3% and 1.3%). Truncation is needed for decreases, and for speed under jumps (102 against 246 steps).
  * The comparison with published change-point methods (E13) has no report part yet; `results/` holds only partial E13 CSVs (offline size, thresholds, online null).
* **Fix.** "the only one of E4's detectors that is valid in both directions under jumps and fat tails; the comparison with published methods (E13) is in progress."

### M13. README: "Both have one-line fixes"

* **Where.** `README.md:30`.
* **Evidence.**
  * The kappa correction is a one-line fix for the bias under constant volatility.
  * ALT3 (S5) shows that with one daily MAD, the corrected estimator still has a bias of minus 11% at n = 23,400 when volatility moves within the day.
  * The fix there is a local scale (`docs/estimator_issues.md:111`: "A local (rolling) MAD avoids this"), which takes more than one line.
* **Fix.** "Both have simple fixes under constant volatility (normalise the weights; divide by κ₄). When volatility moves within the day, the truncation also needs a local scale (E2b, ALT3)."

---

## 4. Minor findings

**E10.**

* **m1. Rounding.**
  * `E10.md:116`: DOGE "deseasonalised only" alarms per day is 0.6848, printed 0.69; it should be 0.68.
  * `E10.md:115`: SOL "EWMA 360" down share is 0.8148, printed 82%; it should be 81%.

  Source: `results/e10c_real_forecast_relative.csv`. All other Q10c numbers match, including "+2% to -17%" and the interval statements.
* **m2. Range and bound wording.**
  * `E10.md:121`: "60% to 71%" for up alarms covers the rows in the table. The half-life 1440 variant in the same CSV reaches down to 57.5% (XRP). Say which variants the range covers.
  * "far above the 0.05 per day that a correct forecast would give" should read "far above the **at most** 0.05 per day that a correct forecast would give".
* **m3. Lai bound wording.** `E10.md:87`: "no detector ... can be expected to notice it much faster than within about two days". Lai's bound is an asymptotic statement, as A grows, about the worst-case delay, and the next sentences report observed delays of about 0.6 times the figure. Write: "the information bound suggests delays of the order of two days; at this A the observed delays are about 0.6 times that (E11)".

**E8 and E9.**

* **m4. Q8d wording.**
  * `E8.md:100`: "at 390 steps it gives away 7.6 to 10.6 points" refers to the upper bound only; say so.
  * The Q8d "oracle" is the best constant-bet bound for the true law, as a median over 500 paths (`E8.md:89`). `math.md:216` calls it "the best that any betting bound can achieve". Use "best constant bet" in both places: it is a benchmark, not a strict limit.
* **m5.** `E8.md:148`: "the larger A alone changes nothing (126 alarms at both A values)". The up/down split moves (73/53 to 74/52) and the up hit rate moves from 81% to 77%. Write "changes almost nothing".
* **m9.** `E9.md:29`: "Short-lived clustering (half-life of a few steps)". The shortest case in the table has a half-life of 14 steps; say so.

**math.md.**

* **m6. A null that depends on the algorithm.** `math.md:228` to `231` state (H↑) and (H↓) relative to the bet-weighted calibration averages μ̄↑_T and μ̄↓_T, whose weights the algorithm chooses from the data. This is correct but unusual: whether "no change" holds depends on the bets. Add a sufficient condition that does not depend on the algorithm. If sup over t > T of μ_t ≤ (1 + b) × min over n_τ < i ≤ T of μ_i, then (H↑) holds for every bet sequence, because a weighted average is at least the minimum. The same holds for (H↓), with max in place of min.
* **m7. Missing returns in Part 2.** `math.md:238` and `:251`, Theorem 6.4 part 2: "monitored steps" also counts steps with a missing return. M is unchanged at those steps, and the proof uses E[M_t | F_{t-1}] ≤ M_{t-1} + P(observed) ≤ M_{t-1} + 1. Say this in the statement.
* **m8. Corollary C.** `math.md:131`:
  * With plug-in recalibration, each segment k has its own τ_k. "True baselines" should read: "m₀ equals the true conditional mean of min(r², τ_k²) for that segment's τ_k".
  * The 30.1 per 10,000 check is for the plug-in baseline, which the corollary does not cover. Label it as an empirical check.

**E12.**

* **m10.** `E12.md:59`: "The pre-release hours are quieter than usual, not noisier". This rests on 6% to 19% of 16 to 22 events against 18% to 26% of controls, and it is not tested. Write "there is no sign of more alarms before releases".
* **m11.** `E12.md:91`: "its alarms mark a short volatility regime rather than a jump". This is supported by the median largest return and by the delays, but not by the "jump only" column at releases (robust 19%, 40%, 35% against Gaussian 33%, 41%, 30%). Mention that.

**ALT3.**

* **m12. Small numerical points.**
  * `RESULTS.md:8`: slopes "minus 0.50 to minus 0.68", but TRQ-U in S5 has -0.48.
  * `RESULTS.md:11` says "0.33 to 0.44 times MedRQ", while `:241` says "0.33 to 0.46".
  * `RESULTS.md:38`: "simulated values sit about 1% above the closed forms throughout". They range from 1.2% **below** (TPQ, 13.48 against 13.65) to 3.2% above (PI-LU, 8.46 against 8.20).
  * The Q1 efficiency column uses the simulated AVAR for PI-LU (8/8.46 = 0.95), while the closed form gives 0.976. Say which one is used.
  * `RESULTS.md:298` and `ALT3.md:59`: "PI-NN has the lowest RMSE of all robust estimators at every n, with and without jumps". At n = 390 with medium jumps, TRQ-M is lower (0.352 against 0.360). At n = 1440 without jumps, PI-NN and TPQ are nearly tied (0.277 against 0.280).
* **m13. Chance level for the burst diagnostic.** `ALT3.md:47` and `RESULTS.md:257` ("isolated jumps on a slowly varying volatility give 0%") rest on 16 exceedances. Under i.i.d. Gaussian returns, the chance that at least one of a return's two neighbours lies beyond 2 sd is 1 - 0.9545² = 8.9%. The expected count among 16 is about 1.4, so 0 of 16 is not unusual. Report the 8.9% chance level. The real 44% to 65% is far above it, so the finding stands.
* **m14. Novelty of the burst finding.** `RESULTS.md:12` calls the burst diagnosis "the main new finding", and `RESULTS.md:319` lists it "above all". There is literature on apparent jumps at high frequency being bursts of volatility, for example Christensen, Oomen and Podolskij (2014, *Journal of Financial Economics*, "Fact or friction: jumps at ultra high frequency") and Bajgrowicz, Scaillet and Treccani (2016, *Management Science*, "Jumps in high-frequency data: spurious detections, dynamics, and news"). I have not read them for this review. Check them before claiming novelty, and cite them if they apply.

**README and strategy.**

* **m15. Betting factor and alarm level.**
  * `README.md:41` writes the betting factor as 1 + λ (X_i - μ₀). The code (`monitor.py:306`, `:309`) uses 1 + λ (X_i / μ₀ - 1) for increases, with λ in [0, 1), and 1 + λ (1 - X_i / μ₀) for decreases, with λ < 1/(B - 1). The README form is valid only with λ < 1/μ₀, and it hides the decrease factor.
  * `README.md:42`: in the two-sided monitor each side alarms at 2A (`math.md:129`, `monitor.py:312`), not at A.
* **m16.** `README.md:38`: "τ ... from past data only (rolling or exponentially weighted)". In the monitor, τ is fixed for each calibration window (`monitor.py:157` to `158`, `:241` to `242`). The rolling τ belongs to E3 and E14.
* **m17.** `README.md:126`: "an information rate that we compute in closed form". In E10 the KL rate is a simulation average (`E10.md:85`); only ρ(k) is in closed form.
* **m18. Counts.**
  * `README.md:5` and `:66` say "About 190 tests"; `RESEARCH_STRATEGY_BAPAT.md:123` says "164 automated tests". The suite now has 220 tests, and all pass.
  * `README.md:5` says "15 experiments", but when this review was written E13 and E15 had no report part, and E14's report part was still being written (E13 and E14 had partial CSVs, E15 none). Count only finished experiments, or say "13 reported, 2 in progress".
* **m19.** `README.md:3`: "the jump-robust truncation of Bapat (2026)". Truncation for jump robustness goes back to Mancini (2009), which is already in the reference list. "The MAD-based truncation used by Bapat (2026)" credits both.
* **m20.** `RESEARCH_STRATEGY_BAPAT.md:16`: "its averaged estimator throws away the local, minute-by-minute information" sounds dismissive of the author. Suggested: "its averaged estimator is a daily summary; the local quantity it averages is what an online monitor needs".
* **m21.** `README.md:124`: "raises a false alarm on at most about 5% of days in the long run". This holds for the true baseline and the confidence-sequence baseline (Corollaries C and 6.6), and only approximately for the plug-in. Say so.

---

## 5. Derivations

### D1. Theorem 6.4: the steps that carry the proof

* **G is known at time T.** ν_i = E[Z_i | F_{i-1}] is known one step ahead, U_t and L_t are adapted, and T is a stopping time. So ν̄↑_T, U_T, and therefore G↑ = {ν̄↑_T < U_T}, are F_T-measurable (known at time T).
* **Step 1.** The event "the bound fails at time T" is contained in "the bound fails at some time". Theorem 6.2 bounds the probability of the second event by δ_u. This works for every random T, which is also why the width-based stopping rule of 6.6 is allowed.
* **Step 2.** On G and for t > T, (H↑) gives μ_t ≤ (1 + b) μ̄↑_T ≤ (1 + b) Ū_T = m↑. So E[f↑_t | F_{t-1}] = 1 + λ (μ_t / m↑ - 1) ≤ 1, and f↑_t ≥ 1 - λ > 0. The decrease side is symmetric, with f↓_t ≥ 1 - λ (B - 1) > 0 because X_t ≤ τ². The bets depend only on data up to T.
* **Step 3.** On E_s = G ∩ {T = s}, the process S_t = 1_{E_s} (M_t + m - (t - s)) is a nonnegative supermartingale on s ≤ t ≤ s + m, with S_s = m 1_{E_s}. Ville's inequality gives P(E_s and an alarm) ≤ m P(E_s) / A. Summing over s gives m P(G) / A. An alarm needs max(M↑, M↓) ≥ 2A, which implies (M↑ + M↓)/2 ≥ A.
* **Part 2.** Optional stopping at the bounded time min(τ_A, s + n) on E_s, then monotone convergence, gives E[1_{E_s}(τ_A - s)] ≥ A P(E_s).
* **Corollary 6.6.** E[N] ≤ δ (1 + E[N]) + m/A. Since N ≤ m is finite, this rearranges to E[N] ≤ (δ + m/A)/(1 - δ).

### D2. The signal cap ρ(k)

* **Two recursions.** The pre-break filter is s_t = ω + α r²_{t-1} + β s_{t-1}. The true post-break variance is v_t = kω + α r²_{t-1} + β v_{t-1}.
* **Their difference.** d_t = v_t - s_t satisfies d_t = (k - 1)ω + β d_{t-1}, so d settles at (k - 1)ω/(1 - β).
* **Means.** With ω = 1 - α - β, E[v] = k and E[s] = k - (k - 1)(1 - α - β)/(1 - β). Hence

  ρ(k) = E[v]/E[s] = k(1 - β) / (k(1 - β) - (k - 1)(1 - α - β)) = k(1 - β) / (kα + 1 - α - β).

* **Shape.** The derivative with respect to k is (1 - β)(1 - α - β)/(kα + 1 - α - β)² > 0, and the limit as k grows is (1 - β)/α.
* **Values.** α = 0.05, β = 0.94: ρ(2) = 0.12/0.11 = 1.091, cap 1.2. α = 0.03, β = 0.969: cap 1.0333. All match `results/e10d_signal_cap.csv`.

### D3. KL rate and Lai delay

* **KL rate.** For two zero-mean normal laws with variances v (after the break) and s (before), KL = ½(q - 1 - ln q) with q = v/s. That is the script's formula, in the right direction.
* **Lai delay.** ln A = ln 28,800 = 10.268. At persistence 0.99 and k = 2: 10.268 / 0.003684 = 2787 steps = 1.94 days (CSV: 2787.4). The other rows agree to rounding.
* **Interval check.** The Wilson interval for 44% of 300 is 38.5% to 49.7%, matching "39% to 50%".

### D4. Truncation constants and the AVAR of clipped RQ

* **Limit under constant volatility.** With constant σ and a cut-off at a standard deviations, E[min(r⁴, τ⁴)] = σ⁴Δ² E[min(Z⁴, a⁴)]. So (n/3) Σ min(r⁴, τ⁴) converges to κ₄(a) times IQ, with κ₄(a) = E[min(Z⁴, a⁴)]/3:
  * scaled MAD (a = 3): κ₄ = 0.9638;
  * raw MAD (a = 3 × 0.6745 = 2.0235): κ₄ = 0.7045.
* **AVAR after the κ₄ correction.** Var(min(Z⁴, a⁴)) / E[min(Z⁴, a⁴)]² = 8.186 at a = 3 and 3.94 at a = 2.02. For untruncated RQ it is 96/9 = 10.67.
* **Caveat for a = 2.02.** The second number does not mean a better estimator: the a = 2.02 version is heavily biased before correction, and the correction is exact only under constant Gaussian volatility.

### D5. Expected number of diffusive returns beyond sqrt(2 log n)

* **Approximation.** n × 2(1 - Φ(a_n)) ≈ n × 2φ(a_n)/a_n = 2n e^{-log n} / (sqrt(2π) sqrt(2 log n)) = 1/sqrt(π log n).
* **Exact and approximate values.** Exact: 0.215, 0.197 and 0.170 at n = 390, 1440 and 23,400. Approximation: 0.231, 0.209 and 0.178.

### D6. Chance level of the neighbour share

Under i.i.d. N(0, 1), the probability that at least one of the two neighbours exceeds 2 in absolute value is 1 - (1 - 2(1 - Φ(2)))² = 1 - 0.9545² = 0.0889.

### D7. ALT3 Proposition 1, Monte Carlo

Blocks of m = 5, with 10⁷ blocks per case. The ratio of the corrected estimate to the target:

| Case | σ⁴ ratio | σ⁸ ratio |
|---|---|---|
| clipping at a = 3 | 0.99969 (s.e. 0.00040) | 0.99984 (s.e. 0.00106) |
| dropping at a = 3.81 | 0.99986 (s.e. 0.00042) | 0.99920 (s.e. 0.00118) |

All four are within one standard error of 1, so the exact unbiasedness claim checks out.

---

## 6. Checked and found correct

### Mathematics (math.md sections 5 and 6)

* **Section 5.4.** Lemma 1; Theorem A (the corrected proof using Fatou's lemma); Theorem B (Ville's inequality on S_t = M_t + m - t); CUSUM ≤ Shiryaev-Roberts; averages of e-detectors; the two-sided split at 2A.
* **Section 5.3.** The impossibility claim and its proof (the two-point law argument), and the hedged link to Bahadur and Savage.
* **Section 6.2.**
  * Lemma 6.1: φ is nondecreasing, and the bound is evaluated at x = -λa.
  * Theorem 6.2: the supermartingale step uses (1 - y)e^y ≤ 1, and the lower bound uses (1 + y)e^{-y} ≤ 1.
  * Corollary 6.3, and the remark on clipping the bounds to [0, 1].
* **Sections 6.3 to 6.4.** Theorem 6.4, parts 1 to 3 (D1); Corollaries 6.5 and 6.6.
* **Numbers in section 6.**
  * A = 9750 and 36,000; ln 9750 / ln 7800 = 1.025.
  * The union bound over K = 20 levels widens the bound by sqrt(7.6/4.6) = 1.29 ("29% wider"). The sample split widens it by sqrt(4/3) = 1.15.
  * ẑ ≈ 0.11, so the lower-bound cap is 0.75/0.11 ≈ 7. The truncated squared normal has sd/mean = 1.39 and a range of 6.5 sd (c = 3).
* **math.md numbers against the CSVs.**
  * 5.6(1): 4.6% to 5.1% against 1.6% to 1.8%, and 2.0% to 2.8% (`e04e_plugin_vs_true.csv`); E8's 4.5% to 4.9% (`e08a_null.csv`).
  * Corollary C: 26.5 and 30.1 per 10,000 against 33.3 (`e04e_restarts.csv`).
  * 6.7: "at most 0.27%" and "at most 0.08%" (`e08a_null.csv`).
  * 6.5: the GARCH coverage figures 1.6%, 0.05% and 0.175% (`e08f_garch_coverage.csv`).

### Code

* **Calibration and bounds.** `detectors.calibrate`, `detectors.betting_edetector`, `confseq.eb_confseq` and `confseq.StreamingEBConfSeq`:
  * the centre, variance estimate and bets use only past data;
  * the ranges are one-sided, a↑ = 1 - ẑ and a↓ = ẑ;
  * the bet cap is c/a;
  * δ/2 is spent on each side.
* **`VolatilityMonitor` with the confidence-sequence baseline.**
  * τ is computed from the first n_τ returns, including the current one, and that return is not fed to the sequence (`monitor.py:222` to `247`).
  * A started sequence is never restarted (`:256` to `261`).
  * Every rule that ends calibration uses only past data.
  * The decrease side switches off when L = 0 or B ≤ 1 (`:274` to `279`).
  * Monitoring starts on the step after calibration (`:280`).
  * Missing returns are skipped (`:303`), and the monitor restarts after each alarm (`:318`).
* **Forecasts and profile.** `forecast.ewma_variance` and `forecast.garch11_filter` use only past data. `data.intraday_scale` with the default `start="past"` has no look-ahead (round 1 m16 is fixed).
* **Wilson intervals.** Correct in the E10 and E12 scripts. Example: 44% of 300 gives 38.5% to 49.7%.
* **Tests.** 220 passed.

### Experiments

* **E8.** All table numbers in e08a to e08g, except C7, m4 and m5. Path counts: 10,000 for the null study, as the strategy says.
* **E9.** All numbers (`e09ab_bound_tightness.csv`, `e09c_garch_band.csv`, `e09d_profile.csv`).
* **E10.**
  * Q10a, Q10b, Q10b2 (including the 39% to 50% interval).
  * The Q10d formula, KL and Lai columns.
  * Q10e, including the scale-break half-life ln 2 / (-ln 0.94) = 11.2 steps.
  * The GARCH half-lives 14, 69 and 693 steps.
  * Q10c: all numbers except m1 and m2, including the alarm-rate changes of +2% to -17% and the random-time "below forecast" rates of 55% to 63%.
* **E12.**
  * Q12a, Q12b, Q12c and Q12d numbers against `e12_tests.csv`, `e12_delay.csv` and `e12_gaussian_vs_robust.csv`.
  * The logic of the stratified permutation test (a Poisson-binomial distribution within weekday-matched strata), and the Holm values.

### ALT3

* Proposition 1 (D7); the closed forms for e2 and e4; 8(k + 3)/(k + 2) and the optimal block length k* = (144n/ξ⁴)^{1/3}; Propositions 2 to 4 by hand.
* TPQ 13.65 and QPQ 14.92.
* The clipped and undetected jump residuals (17.55, 27.01, 8.73, 13.33, 21.68, 45.0), and 6.5% and 12.0% per jump at a_n.
* The oracle block AVAR 7.79, and `avar_block(6)` = 9.0.
* The printed cross-validation expectations, both the printed values (74.3 to 295.0) and the corrected ones.
* The κ₄ values; the definitions in the burst diagnostic; the burst model's match to the real estimator ratios; "56% to 59%" for RK and RK-k in the burst model; and "20% to 51%" for the local plug-ins.

### README and strategy numbers traced to their source and found correct

* E4: 70% to 97% (Gaussian LR and CUSUM), 0.3% to 2.3% (robust e-detector), 36 to 49 steps; robustness costs 10% to 70% extra delay.
* E9: 1.1% to 1.4%, and "about four times conservative".
* E4e: the plug-in roughly triples the false alarm rate.
* E10: 44% against 25%; 82% to 100% (raw returns) and 3% to 6% (fitted GARCH).
* E6: 71% to 83% (Q6e); E6d: 72% against 3%.
* E12: about 30% control rate for the up outcome (both union and matched); blind 70% of the time (`E12.md:31`).
* H3: 14 of 18 cells; E9d: 10 microseconds per step; E9: 4000 paths.

### References

The bibliographic details in the README and in math.md section 6 agree with what I know of these works: Bahadur and Savage 1956; Lai 1998; Lamoureux and Lastrapes 1990; Mikosch and Starica 2004; Mancini 2009; Lee and Mykland 2008; Page 1954; Shin, Ramdas and Rinaldo 2024; Waudby-Smith and Ramdas 2024; Fan, Grama and Liu 2015; Howard et al. 2020 and 2021. I did not re-fetch them.
