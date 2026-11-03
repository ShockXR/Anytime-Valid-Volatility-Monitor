# Independent review, round 1

> *Note added 3 November 2026: the "independent reviewer" or "referee" below was a separate AI agent (Claude, by Anthropic) instructed to act as a sceptical referee. It was not a human expert, and the review has not been confirmed by a human.*

Reviewer: REV1 (independent research review, not an author of the reviewed files). Date: 2 November 2026.

**Scope.** `docs/math.md` sections 1 to 5, `docs/estimator_issues.md`, `docs/paper_review.md`, `src/avmonitor/{robust,detectors,monitor,simulate,data,forecast}.py`, E1 to E5 in `docs/experiment_report.md`, E9 in `docs/report_parts/E9.md`, `README.md` and `PITCH_BAPAT.md`. The reviewed documents did not change between commit `b08b7a3` (start of review) and `98a1e83`. During the review, `monitor.py` and `detectors.py` gained a confidence-sequence baseline (`e90cfd0`, `d0c539e`) and `math.md` gained section 6. I reviewed the plug-in path of the current code. I did not review the confidence-sequence path, `confseq.py`, math.md section 6, or E6 to E8 and E10.

**Method.** I read every claim, checked it against the paper text and `results/*.csv`, and re-ran cheap checks. As the review brief asks, the scratch scripts are kept outside the repository, in `/tmp/claude-0/rev1/`: `check_constants.py`, `check_identity.py`, `check_code.py`, `check_stats.py`, `check_e5_data.py` and `check_e5c_control.py`. Each one runs in under 30 seconds. Labels: [SYNTHETIC] for my simulations (fixed seeds in the scripts), [REAL] for the Bitstamp file used in E5, [PUBLISHED] for numbers printed in Bapat (2026). If the authors want any of these numbers in the main docs, the scripts should first be moved into `experiments/` and their outputs saved under `results/`.

## Summary verdict

* **The mathematics of sections 1 to 5 is correct.** I found no wrong theorem, formula or identity. The truncated-moment formulas, the truncation-bias argument, the identity "averaged kernel = truncated RQ up to edge effects", the test supermartingales, the Shiryaev-Roberts and CUSUM recursions, Theorems A and B, the two-sided split at 2A and the impossibility proof all hold. Several statements need more precision (m1 to m4, m6, m7). One useful result is missing: a corollary that covers recalibration after alarms (m5). I state it below and checked it numerically.
* **The claims about the paper are well founded and mostly fairly worded.** One arithmetic slip needs fixing (C1). One reading of our own counterfactual undersells the paper's estimator (M1).
* **The code has no look-ahead in the detector or monitor path.** The decrease bets stay inside their admissible range everywhere, and the 2A logic is right. There are latent bugs: E5c turns missing returns into zeros, and that line has been copied into E8 (M6). The vectorised API lets a single NaN silently disable detection (m17). There are also two small start-up look-aheads (m16).
* **The numbers in E1 to E5 and E9 match the CSV files,** apart from rounding. Three "Reading" paragraphs claim more than their tables show: the E4c band (M2), the E5c real-data validation (M3: a naive trailing-window rule scores the same 82% on increases), and the adaptive bet in Q4c (m10).
* **README and the pitch overstate several things.** They overstate the impossibility result (C2), the role of the kernel estimator and of "online MAD" in the monitor (M5), noise robustness (M5), the false alarm guarantee under the default plug-in baseline (M4), and the novelty of the theory (M7).

---

## Critical (wrong math or wrong numbers)

### C1. The RQ bias check omits a cross term (`docs/estimator_issues.md:57`)

**Text.** "With σ = 0.0215, (1/3) Σ r⁴ (eq. 3 without the factor n) has expected value 8.3 × 10⁻⁵ from the diffusion plus 2.4 × 10⁻⁵ from the jumps, giving 1.07 × 10⁻⁴. This matches."

**Evidence.** Take r = σε + BZ with B ~ Bernoulli(p) and Z ~ N(0, σ_J²). Then E[r⁴] = 3σ⁴ + p(6σ²σ_J² + 3σ_J⁴). Over n = 390 returns, the jump part of (1/3) Σ r⁴ is n p (6σ²σ_J² + 3σ_J⁴) / 3 = 3.34 × 10⁻⁵. The quoted 2.44 × 10⁻⁵ is this value without the 6σ²σ_J² term. The correct total is 1.167 × 10⁻⁴. A Monte Carlo check with `simulate.paper_returns` (20,000 paths, σ = 0.0215, p = 0.01, σ_J = 0.05) gives a mean of 1.169 × 10⁻⁴ (standard error 4 × 10⁻⁷) [SYNTHETIC]. The published RQ biases are 9.94 × 10⁻⁵ (Table 3) and 1.07 × 10⁻⁴ (Table 4) [PUBLISHED]. So σ = 0.0215 overshoots them by 17% and 8%. A σ of 0.0204 or 0.0209 would match each one exactly (`check_constants.py`).

**Impact.** Low. The conclusion still holds: RQ in the tables lacks the factor n, because with that factor the value would be 390 times larger. But the author can check "This matches" in one line and will find the slip.

**Fix.** Replace the sentence with: "...plus 3.3 × 10⁻⁵ from the jumps (including the cross term 6σ²σ_J²), giving 1.17 × 10⁻⁴. That is within 8% of Table 4 and 17% of Table 3, the right order of magnitude for RQ without the factor n (with it, the value would be 390 times larger)."

### C2. README and the pitch misstate the impossibility result (`README.md:10`, `PITCH_BAPAT.md:13`)

**Text.** README: "Truncation is what makes 'volatility went down' detectable. Without an upper bound on the monitored statistic, no nonparametric e-test can detect a decrease." Pitch: "Without an upper bound on the monitored statistic, **no** nonparametric e-test can ever detect a fall in volatility ... Your truncation provides the bound."

**Problem.** `math.md` 5.3 correctly proves this: for the null "E[Y] ≥ 1 over all laws on [0, ∞)", every e-variable satisfies E(y) ≤ 1. The paraphrases say something different, and in two ways it is wrong.
1. Suppose "volatility" means the variance E[r²]. Then no nonparametric e-test can detect a fall, *whatever statistic it computes*, truncated or not. The null "E[r²] ≥ m₀" itself admits no nontrivial e-variable. Truncation does not make a fall in variance detectable. It changes the target to the truncated variance E[min(r², τ²)], which is bounded and therefore testable. Linking that target back to σ needs the fixed-shape assumption of `math.md` 5.1.
2. Truncation is not the only route to a bounded statistic. Exceedance indicators 1{|r| > q}, ranks and any bounded increasing transform of r² also allow nonparametric detection of a decrease. "Truncation is what makes ... detectable" and "Your truncation provides the bound" overstate the role of truncation.

A reader who works on sequential analysis will raise both points at once, so these sentences matter most to get right.

**Fix.** Pitch: "A fall in the variance itself cannot be detected by any nonparametric e-test (`math.md` 5.3, a one-sided form of the impossibility result of Bahadur and Savage, 1956). A detector must therefore monitor a bounded functional of the returns. Truncation at c × MAD gives one with a direct volatility meaning: if the shape of the innovations is fixed, the truncated variance falls exactly when σ falls." README: the same in one sentence.

---

## Major (unsupported claims, latent bugs that would change results)

### M1. The normalised counterfactual is misread, and the misreading undersells the paper (`docs/estimator_issues.md:74`, `docs/experiment_report.md:31`)

**Text.** "The plain kernel has a lower relative RMSE than the raw-MAD robust kernel when jumps are small (0.17 against 0.32 at σ_J = 0.01). The robust version wins clearly only when jumps are large (10.6 against 0.28 at σ_J = 0.10)."

**Evidence** (`results/e01_normalised_counterfactual.csv`) [SYNTHETIC]. At σ_J = 0.05, the raw-MAD robust kernel already wins in all 15 cells: relative RMSE 0.27 to 0.31 against 0.48 to 1.28, a factor of 1.6 to 4.7 (medians 0.29 against 0.83). The scaled-MAD robust kernel beats the plain kernel in **all 45 cells**, including σ_J = 0.01 (0.148 to 0.154 against 0.164 to 0.172).

**Why it matters.** This is the one place where the docs compare the paper's estimator with its plain counterpart on a fair scale, and the sentence makes the estimator look worse than our own table shows.

**Fix.** "With the raw MAD, the normalised robust kernel loses to the plain kernel only when jumps are small (σ_J = 0.01: 0.32 against 0.17). It wins at σ_J = 0.05 (0.29 against 0.83, medians) and at σ_J = 0.10 (0.28 against 10.6). With the scaled MAD it wins in all 45 cells."

### M2. The null band does not bring false alarms to target at b = 0.25, and is not a guarantee (`docs/experiment_report.md:137`, `docs/math.md:144`)

**Text.** "A band of 0.25 to 0.5 brings drift-driven false alarms back to the target range." `math.md`: the band "restores control".

**Evidence.** In `results/e04c_null_band.csv`, b = 0.25 gives 8.3% (increase side) and 5.0% (decrease side) per day, with each side run alone at A = 7800. The increase side is above the 5% target. Only b = 0.5 is inside it (1.7% and 0.0%). Nor does the band make the null true in this design. With a log-volatility standard deviation of 0.15 (`change_returns`, `sv_xi = 0.15`), the variance exceeds 1.25 times its mean on about 19% of steps and 1.5 times on about 7%. So a static band cannot restore the guarantee. Low false alarm rates here are an empirical finding about this design.

**Fix.** E4c: "A band of 0.5 brought drift-driven false alarms below 5% in this design (1.7% and 0.0%). A band of 0.25 reduced them from about 25% to 8.3% and 5.0%." `math.md` 5.6(2): "The band reduces false alarms (E4c). It restores the guarantee only if the conditional truncated variance stays inside the band, which slow log-normal drift does not do."

### M3. The real-data validation cannot tell the detector from volatility clustering (`docs/experiment_report.md:174`, `README.md:11`, `PITCH_BAPAT.md:20`)

**Text.** "Followed by a persistent shift in volatility in the direction it flagged, far more often than chance (81% and 90% with band 0.25, against 50%)". README: "Alarms on real data mean something."

**Problem.** The control compares the next day with the previous day at random times. An alarm, by construction, fires when recent volatility has already moved, so random times are not a matched control. The caveats paragraph says that volatility clustering makes any reacting detector look predictive, but the control does not address it.

**Evidence** (`check_e5c_control.py`) [REAL]. I used the same deseasonalised Bitstamp series and the same measure (robust sd of the next 1440 minutes divided by the baseline). At 20,000 random minutes I applied a naive rule: "flag *up* if the robust sd of the last hour is more than 1.5 times that of the preceding day". It is followed by higher volatility 82% of the time (median ratio 1.46, n = 3538). The mirror rule (last hour below baseline / 1.5) is followed by lower volatility 73% of the time (median 0.78). The robust e-detector with band 0.25 scores 81% (median 1.43) and 90% (median 0.65). The increase figure is therefore what volatility clustering alone delivers. The decrease alarms do look more informative than the naive rule (90% against 73%). One caveat: the naive rule fires far more often and has no false alarm control, so this does not make the detector useless. It shows that the share metric cannot separate the two.

**Fix.** Report a naive trailing-window rule (ideally matched on alarm rate) next to the detector, or compare forecasts out of sample (E10 may do this). Reword README:11 and pitch item 4. For example: "81% of increase alarms and 90% of decrease alarms were followed by a day of volatility in the flagged direction. A naive last-hour rule scores 82% and 73% on the same measure, so the increase figure mostly reflects volatility clustering."

### M4. With the default plug-in baseline, the guarantee is approximate and the effect is larger than the docs suggest (`README.md:3`, `PITCH_BAPAT.md:7`, `docs/math.md:143`, `docs/experiment_report.md:139`)

**Evidence** (`check_stats.py`, T4) [SYNTHETIC]. Increase side, A = 7800, one 390-step day, 4000 paths, Monte Carlo standard error 0.34 points:

| Law | true baseline | plug-in, 390 calibration steps | plug-in, 1440 steps |
|---|---|---|---|
| Gaussian | 1.75% | 4.57% | 2.75% |
| Student t3 | 1.73% | 4.67% | 2.55% |
| Jumps (0.5%, 8 sd) | 1.63% | 5.10% | 2.00% |

With 390 calibration steps, the plug-in baseline roughly triples the false alarm rate, and with jumps it reaches the bound. The 5% bound is "kept" only because Ville's inequality has about a factor of 3 of slack (E9). E4d's 5.4% is within one Monte Carlo standard error of 5%, so E4d alone does not show anti-conservativeness. The informative comparison is with the rate under the true baseline.

**Fix.** In `math.md` 5.6(1) and E4d, compare with the true-baseline rate as above. In the headlines (README:3 "guarantees that hold no matter how often you check it", pitch:7), add: "proved for a known baseline or with the confidence-sequence baseline (section 6). With the default plug-in baseline it holds approximately (E4d)." Consider making `baseline="cs"` the default once E8 is final, or else defaulting to `cal_len ≥ 1440` with `margin_z = 1`.

### M5. README and the pitch describe the monitor inaccurately (`README.md:3`, `PITCH_BAPAT.md:7`, `PITCH_BAPAT.md:18`, `docs/paper_review.md:66`)

1. README:3 says the monitor "combines jump- and noise-robust realized measures". The repository's own E2c reading says "No estimator here handles noise", and `estimator_issues.md` 4b explains why truncation does not address noise.
2. Pitch:7 says "We use exactly that property [a robust, bounded, local kernel measure] to build a sequential monitor". The default monitor uses no kernel smoothing. It bets on single truncated squared returns min(r², τ²) with a fixed grid. The kernel enters only through the optional adaptive bet (`betting_edetector(adaptive_h=...)`). `VolatilityMonitor` does not use that bet, and it was slower than the fixed grid in 14 of 18 E4 cells. What the monitor takes from the paper is eq. (11), the c × MAD truncation.
3. Pitch:18 says "online MAD", and paper_review:66 says the prototype "replaces the global MAD with an online one". In fact the monitor estimates the MAD once per calibration window and keeps it fixed for the segment (`monitor.py`, `_plugin_stats`). An online MAD appears only in E3.

**Fix.** Pitch: "We use the truncation step (eq. 11) to bound the monitored statistic. Your local kernel quarticity can size the bets, but it did not beat a fixed grid in our tests. The MAD comes from a calibration window and is re-estimated after each alarm." Remove "noise-robust" from README:3 or qualify it.

### M6. Missing returns are turned into zeros before the monitor sees them (`experiments/e05_real_btc.py:106`, copied into `experiments/e08_cs_baseline.py:354` "as in E5c")

`mon.run(np.nan_to_num(u, nan=0.0))` bypasses the monitor's own NaN skip. A missing minute becomes a zero return, which counts as evidence for a volatility decrease.

**Evidence** (`check_stats.py`, T3) [SYNTHETIC]. Gaussian stream, `cal_len` 1440, A = 28,800, then 120 missing minutes. With NaN passed through, there is no alarm. With NaN replaced by 0, a "down" alarm fires 69 minutes into the gap. On the Bitstamp file only 1 value of `u` is non-finite (`check_e5_data.py`) [REAL], so the E5 results are unaffected. On Binance (a NaN row for every missing minute) or Nifty (NaN at every intraday gap), E8 or any other copy of this line will produce false decrease alarms during outages. A related problem: `intraday_scale(r.fillna(0.0))` (e05:98, also in e08, e10 and `cli.py:81`) counts missing minutes as zero returns in the profile medians. E6b passes NaN, which is the right choice.

**Fix.** Pass `u` unchanged (the monitor skips NaN), and call `intraday_scale` on the returns with their NaN.

### M7. The novelty of the theory is overstated (`PITCH_BAPAT.md:38`, `PITCH_BAPAT.md:13`)

The pitch proposes "a joint short paper: theory (guarantees and the impossibility result)". Theorem A is the average-run-length guarantee of Shiryaev-Roberts e-detectors (Shin, Ramdas and Rinaldo 2024). Theorem B follows from the same supermartingale and Ville's inequality; check whether that paper already states it before presenting it as new. Section 5.3 is a one-sided relative of Bahadur and Savage (1956). `math.md` rightly says "Following Shin, Ramdas and Rinaldo (2024)", but the pitch presents these results as contributions.

**Fix.** List what is new: a truncated-variance null with a matching decrease detector, the estimator audit, the real-data application, and the confidence-sequence baseline if E8 holds up. Cite the known results as tools.

---

## Minor (precision, wording, small or latent code issues)

**m1. The proof of Theorem A uses the wrong condition (`math.md:116`).** "The increments of M_t - t are integrable on {t ≤ τ_A} and optional stopping applies" is not, in general, the right condition. The theorem is still true. One fix: on {τ_A ≥ t} we have M_{t-1} < A, so E[|ΔX_t| | F_{t-1}] ≤ 2(A + 1) + 1, and bounded conditional increments with E[τ_A] < ∞ allow optional stopping. A simpler proof: E[M_{τ∧t}] ≤ E[τ∧t] ≤ E[τ] for every t. Since M ≥ 0 and M_{τ∧t} → M_τ ≥ A almost surely when E[τ] < ∞, Fatou's lemma gives A ≤ E[τ_A].

**m2. State the assumptions behind the link to volatility (`math.md:84`, 5.1).** The innovation ε_t must be independent of (σ_t, F_{t-1}) with a fixed law, and σ_t ≤ σ₀ must hold almost surely. Without independence (for example, contemporaneous leverage), the conditional mean need not be monotone in σ. Also name the filtration: the returns plus the calibration data, with τ, m₀, the bets and the deseasonalising profile all adapted to it. Everything in the monitor is adapted, except the intraday profile in its first days (m16).

**m3. Three small corrections in `math.md:125` to `:127`.**
- An average of several M statistics does not obey the same recursion. It obeys E[M̄_t | F_{t-1}] ≤ M̄_{t-1} + 1, which is all Lemma 1 needs.
- The two-sided guarantee holds under the *two-sided* null, (1 - b) m₀ ≤ E[X_t | F_{t-1}] ≤ (1 + b) m₀. Say so.
- "Run length above 2.4 times A": E4a gives 2.31 to 2.40 times A, censored at 20,000 steps.

**m4. Section 5.3 can be stated more strongly (`math.md:96` to `:100`).** For y ≥ 1, the point mass δ_y is in H, so E(y) ≤ 1 for *every* y. No e-variable for H ever exceeds 1. Cite Bahadur and Savage (1956). It is also worth stating the asymmetry: the increase null E[Y] ≤ 1 on [0, ∞) does admit e-variables (E(y) = y). That is why the increase detector works on raw r² in E4.

**m5. Recalibration after alarms is missing from the theory.** The monitor restarts after each alarm, but `math.md` covers only the first alarm. One corollary covers the restarts. Let G_t be the sum of the two-sided statistic (M↑ + M↓)/2 at past alarms, plus its current value, minus the number of monitored steps. G is a supermartingale with G₀ = 0, and every alarm banks at least A. So E[number of false alarms in m monitored steps] ≤ m / A, provided the null holds in every segment relative to that segment's baseline. Numerical check (`check_stats.py`, T5) [SYNTHETIC]: with the true baseline, A = 300 and 400 paths of 60,000 steps, there were 159 alarms per path against a bound of 200. With the streaming monitor, plug-in recalibration and `cal_len` 390, there were 30.1 alarms per 10,000 monitored steps against a bound of 33.3. Add this to section 5.4.

**m6. State the condition for skipping NaN (`monitor.py`, `update`).** Skipping a step amounts to a zero bet, chosen after seeing that r_t is missing. This is valid if the null holds conditionally on the return being observed: missing completely at random, or missingness decided by the past. Add this to `math.md` 5.6. Also, calibration counts NaN toward `cal_len`, so `cal_len` overstates the effective sample.

**m7. "Variance" should be "second moment" (`math.md:137`, 5.5).** It is the second moment E[Y²] = E[r̃⁴] / m₀² that is a truncated local quarticity. The variance v also subtracts μ².

**m8. "Absorbed by the band" is unsupported (`math.md:145`, 5.6(3)).** No experiment in E1 to E5 shows that errors in the deseasonalising profile are absorbed by the band (E5c also runs with band 0). Cite E7's deseasonalisation check if it supports this. Otherwise soften the sentence.

**m9. Price noise makes the conditional null false (`experiment_report.md:99`).** The text says the e-detectors "stay inside their guarantee whenever the null is true, even with ... noise". Price-level noise makes returns MA(1). With noise sd 0.5, r_t = a_t + θ a_{t-1} with θ = -0.172 and σ_a² = 1.457, so E[r_t² | past] = 1.457 + (θ a_{t-1})². That exceeds the unconditional 1.5 whenever |θ a_{t-1}| > 0.21, which happens on about 32% of steps. The conditional null is false, so the 3.7% shows robustness to a mild violation (like the GARCH case in E9). It is not a case where the null holds.

**m10. The adaptive bet was slightly slower, not equal (`experiment_report.md:127`, Q4c).** "Within a few percent of the fixed grid in every cell (for example 86 against 87, 101 against 102, 74 against 71)." In fact the adaptive mixture was slower in 14 of 18 cells, by up to 9% (Student t5 at sd × 2: 24 against 22). Two of the three examples quoted (86 against 87, 101 against 102) are among the 4 cells where it was faster. Write "slightly slower in 14 of 18 cells (by up to 9%)".

**m11. Rounding and counting in `docs/experiment_report.md`.**
- `:107` to `:111`: the text says "≤ 4%" where the CSV has 4.3% (robust grid, Gaussian × 1.5; raw r², Gaussian × 1.25 and × 1.5); "(0%)" where it has 0.3% (CUSUM, jumps × 1.25 and Student t5 × 1.5); and "≤ 5%" where it has 5.3% (raw r², jumps × 1.25).
- `:65` (E2d): the correlations are 0.982 to 0.989, not up to 0.990. The mean ratios are 0.996 to 1.001, not from 0.997.
- `:164`: "A = 28,800 per side" is wrong. The monitor runs each side at 2A = 57,600, and 28,800 is the two-sided level.
- `:168` to `:170`: "Up / down" does not add up to "Alarms" for band 0.25 (73 + 52 = 125, not 126) or band 0.5 (59 + 32 = 91, not 92). The cause is `e05_real_btc.py:112`, which drops alarms with fewer than 1000 post-alarm minutes from the direction counts. (Line 118 aligns directions by truncating the list, which is correct only because the dropped alarms are the last ones.) Add a footnote or count all alarms.

**m12. Wording in `docs/estimator_issues.md`.**
- `:17`: "behave exactly as the printed formulas predict" should be "closely, except the σ_J = 0.10 robust rows (4h)".
- `:163` (4h): "2 to 10 times larger" should be "1.3 to 22 times (median 5.3, quartiles 1.9 to 7.2)".
- `:159` (4f): MinRQ also beats the robust kernel for ETH (1.25 × 10⁻⁶) and SOL (8.46 × 10⁻⁷), and the plain kernel does for SOL (2.19 against 2.20 × 10⁻⁶). The quotes around "RMSE" can read as sarcastic. Turn the point into a question instead: a single estimate of 1.5 × 10⁻¹³ cannot have an RMSE of 2.12 × 10⁻⁶ against an assumed 2.2 × 10⁻⁶, while the plain kernel's 3.98 × 10⁻¹¹ gives 2.20 × 10⁻⁶. So how was RMSE computed for Table 6?
- `:73`: "RQ is about n times the per-step target" holds for (1/3) Σ r⁴, as in the tables. With the factor n it is n² times. Say which one is meant.
- `:155` (4d): μ₁⁻⁴ is the constant of quad-power quarticity, μ₁⁻⁴ n Σ |r_i||r_{i-1}||r_{i-2}||r_{i-3}|. It is worth asking whether eq. (4) was meant to be QPQ.

**m13. Pitch wording (`PITCH_BAPAT.md:26`, `:19`).** Line 26 says "exactly the raw weight sum divided by n". The published ratios are 4.00 to 4.03 × 10⁻⁵, against W(30) / (n - 30) = 3.986 × 10⁻⁵ for n = 264,960. Say "within 1% of". On line 19, "36 to 49 minutes" are synthetic steps.

**m14. `docs/paper_review.md`.** Line 26 says "Andersen, Dobrev and Schaumburg 2012, extended by Schulz 2010": the chronology is odd, and the paper cites the 2008 working paper. Line 58 says "The real-data RMSE is uninformative". Say instead that it "depends mainly on the assumed true value".

**m15. `docs/report_parts/E9.md:17`.** The statement that a user "could lower A by a factor of up to about 3 and still see about 5% per day" was not run; it is an extrapolation from other horizons. Run it or write "we expect". Also, "Only the mean of the bounded statistic matters" is true of the bound, not of the observed rates.

**m16. Small start-up look-aheads in the code.**
- `robust.py:439` against `:455`: the docstring of `ew_quantile_scale` says that entry t uses only data before t. But the start value is the median of the first 50 values, so entries 0 to 49 use future data. Check: multiplying r[10:50] by 100 changes out[0] from 0.96 to 64.4 (`check_code.py`). It is used only in E3.
- `data.py:80`: `intraday_scale` fills the first days (until 10 earlier days exist) with the median profile of days 0 to 19, so days 1 to 10 use up to 19 days of future data. Check: changing only day 15 changes the scale on days 1 to 10 (and, correctly, on days 16 to 35). The docstring says so, but `experiment_report.md:164` calls the profile "predictable". Use NaN or a flat profile until enough history exists, or say so in E5.
- `data.py:84`: the circular smoothing wraps the end of the day into the start. For session markets (Nifty) it blends the close into the open. E7 seems aware of this. A parameter would make it explicit.

**m17. Silent NaN handling in the vectorised API (`detectors.py:97`, `betting_edetector`).** With one NaN in the window, `calibrate` (plug-in) returns NaN for both the scale and m₀. `betting_edetector` then never alarms. Check: one NaN at step 500 followed by a × 3 volatility increase gives alarm -1, while the same data with the NaN set to 0 alarm at step 506. For now only synthetic data use this API. Fix: drop NaN in `calibrate`, as the cs branch already does, and either skip NaN steps (factor 1) or raise an error. Also, `betting_edetector(direction="down")` does not check that x ≤ upper. A value above the bound makes a bet factor negative. Add a check.

**m18. Stale facts in README.** "26 tests" (`README.md:5`, `:84`; `PITCH_BAPAT.md:44`): there are now 62 test functions, or 117 collected test cases, and all pass. The data table says "Free exchange APIs", but the data come from public GitHub mirrors because exchange APIs are blocked. Both are in flux with E6 and E7.

**Fixed during the review.** At `b08b7a3`, the monitor raised `ZeroDivisionError` when a calibration window had more than 50% zero returns (checked in `check_code.py`). The current monitor (`e90cfd0`) extends the window instead. Keep the regression test.

---

## What I checked and found correct

* **`math.md` sections 1 to 2.** The CLT for realized variance; Var(V̂_t) = 2σ⁴ Σ w_k²; the mean of the printed estimator W(h) σ⁴ Δ²; W(10), W(30), W(50) = 3.73, 10.56, 17.39 (the linear approximation 0.341h + 0.32 is within 0.02).
* **Section 3.** Both truncated-moment closed forms (numerical integration gives κ₄(2.0235) = 0.704525, equal to the closed form), the κ table, and a = 0.6745c for the raw MAD.
* **Section 3.1.** The MAD shrinks at the rate of the returns, so the clipped share does not vanish. The rate condition of Mancini (2009) is stated correctly.
* **Section 4.** The identity is exact: on a test case the two sides agree to machine precision, and every interior return gets weight W(h) (`check_identity.py`). Edge returns are 2h / n = 15% of the sample at n = 390, which fits the correlation of 0.98 to 0.99 in E2d.
* **Sections 5.2 to 5.4.** Positivity of both bet types and the admissible decrease range; the test supermartingales for every start time; M_t = f_t(M_{t-1} + 1) and C_t = f_t max(C_{t-1}, 1); Lemma 1; Theorem A (with the m1 fix); Theorem B through Ville's inequality; C ≤ M; mixtures; the two-sided split at 2A; the 5.3 proof; and the aGRAPA form of the adaptive bet.
* **Paper claims, checked against the extracted text.** Eq. (6) gives K(0) to the newest return. Eq. (15) indexes the weights so that the newest return gets K(1) (the printed sum index is i while the summand uses s, so "looks like a typo" is a fair description). The proof of Theorem 1 does assume a slowly shrinking τ, and Step 4 of Theorem 2 does say "τ → ∞". All five Table 6 Kernel / RQ ratios (4.004 to 4.026 × 10⁻⁵) lie within 1% of W(30) / (n - 30) = 3.986 × 10⁻⁵. Table 7 has exactly three exceptions, all at σ_J = 0.10 and h = 10. The Table 8 noise reading checks out (medians on the noise > 0 rows: 1.12 and 0.92 if read as sd, 6.3 if read as variance). The TPQ cross-coin exponent from the Table 5 interquartile ranges is 1.45 to 1.49 (so p ≈ 0.5). The BPQ unbiasing constant 1 / (μ₁ μ₃) is correct. In `robust.py`, the MinRQ and MedRQ constants check out by Monte Carlo (1.0016 and 0.9997), and the TPQ and QPQ constants are standard. E6a (committed during the review) finds the raw-MAD robust kernel at 1.49 × 10⁻¹³ on Binance BTC against the published 1.54 × 10⁻¹³, which supports the raw-MAD inference of Issue 2.
* **Code.** Predictability holds in `_predictable_local_moments` (changing y[t] leaves entry t unchanged and changes entry t + 1), in `rolling_mad_zero` and in the forecast filters (the output for t is written before the update). In the monitor, τ and m₀ are fixed per segment, monitoring starts on the step after calibration, and recalibration uses fresh data. The decrease bets stay inside their admissible range everywhere: grid below 1/(B - 1), adaptive bet at most 0.5/(B - 1), monitor fractions at most 0.95. E4 uses the largest bound across paths, which is conservative, and E4c and E4d use the correct bounds. The 2A logic in both monitor classes is right. In `simulate.py`, the stationary mean of the log-OU process, the stochastic-volatility normalisation E[exp(2z - 2ξ²)] = 1 and the unit-variance Student t all check out.
* **Results against the CSV files.** E1 (every median recomputed from `e01_reproduce_tables.csv`), E2a to E2d, E3, E4a, E4b (up to rounding), E4c, E4d, E5a, E5b (all nine median ratios reproduced exactly), the E5c summary and both E9 tables. The E5 data description reproduces from the raw file: 264,960 rows, no time gaps, 0.30% zero-volume minutes, 5.73% zero returns, quartiles -2.11 and 2.13 × 10⁻⁴, mean -7.65 × 10⁻⁷, extremes -0.0171 and 0.0117 [REAL].
* **Tests.** `python -m pytest -q`: 117 passed in 46 s.

## Suggested order of work

1. Before anything goes to the author: C2, M5, M7 and C1 (one-sentence edits in the pitch, README and `estimator_issues.md`), then M1.
2. M3 (add the naive-rule benchmark to E5c, about a 30-line change) and M4 (one extra column in E4d, and the headline wording).
3. M6 (one line in E5c and E8; rerun E8 if it uses Binance or Nifty), M2, then the minor items.
