# Independent review, round 8: TH6 (false alarms after the first day; Theorem D restated)

> *Note added 3 November 2026: the "independent reviewer" or "referee" below was a separate AI agent (Claude, by Anthropic) instructed to act as a sceptical referee. It was not a human expert, and the review has not been confirmed by a human.*

**Reviewer role.** A sceptical referee for a statistics journal and a code auditor. I did not write any of the files under review and I did not take their claims on trust: every statement below was rederived, and every number I quote was either read from the committed CSV files with my own code or recomputed with an independent simulator.

**Scope.**

1. `docs/theory_per_period.md`: Propositions 1 to 3, Theorem 1, Lemmas 1 to 3, Theorem 2 (a) to (e), the remarks, sections 5 to 7.
2. `docs/math.md`: the remark after Theorem B, the corrected bullet under "Consequences used in the code", and the restated Theorem 6.4 (Theorem D) with its hypotheses, part 4, the bet-weighted remark, Corollaries 6.5 and 6.6, "Choosing A and δ", sections 6.5 and 6.7 (`git show 8e43ce8 -- docs/math.md`).
3. Code and results: `src/avmonitor/periodic.py`, `tests/test_periodic.py`, `experiments/th6_per_period.py`, `results/th6_*.csv`.
4. `docs/report_parts/TH6.md`.

**Read first.** README, `docs/math.md` sections 5 and 6, `docs/review/referee_report_paper.md` (M2 and M4), `src/avmonitor/monitor.py` and `src/avmonitor/detectors.py`.

**State of the tree.** TH6 was committed in 8e43ce8. Later commits (6dec15a, d3fd745, 9be152c and the paper commits up to 8ccff9b) did not change any file in scope: `git diff 8e43ce8 HEAD` on the files listed above is empty, so every line number below holds at both commits. Other agents had modified and untracked files in the working tree; none is in scope and none was touched.

**What I ran.** Scripts are in the session's scratchpad directory, folder `r8th6/` (not committed); nothing was written to `results/`. Everything ran under `nice -n 10`, one process at a time, in about two minutes in total.

* `python -m pytest -q tests/test_periodic.py`: **14 passed in 1.6 s.**
* Closed forms, by hand and with scipy: the window bound, the certified thresholds, τ, m₀, B and 1/B for Gaussian clipping (section 3, item 4).
* CSV checks with pandas: every number named in the brief, and every row of the tables in `docs/theory_per_period.md` section 5 and `docs/report_parts/TH6.md`.
* `r8th6/indep_check.py` and `r8th6/indep_check2.py`: an **independent simulator** that does not import `avmonitor.periodic` (my own two-sided Shiryaev-Roberts and CUSUM recursions with the default bets, my own bold play, my own GARCH tolerance computation), with new seeds. Results in section 4.

## 1. Verdict in brief

All theorems and propositions of TH6 are **correct**. I found no false statement among Proposition 1, Theorem 1, Proposition 2, Lemmas 1 to 3, Theorem 2 (a) to (e) and Proposition 3, and the restated Theorem D (parts 1 to 4) is correct with its proof. The constants 5.146% and 5.044% are right. Every headline number I checked matches its CSV, and the ones I recomputed independently agree within sampling error. There is **no Critical finding**.

There are **two Major findings**, both overclaims about *what is guaranteed where*, which is the subject of M2 itself, and a list of Minor ones. The most important Minor finding is a wrong number in the new remark after Theorem B (m4): the cited range "4.4% to 4.9% on days 10 to 20" belongs to the 1440-step table, and the 390-step file it cites gives 4.34% to 5.18%.

| Item | Verdict |
|---|---|
| 1. Proposition 1 | Correct (wording on sharpness overstated: m1) |
| 2. Theorem 1 | Correct; κ = 0 for the default restart; the conditional statement is justified. The extension to estimated baselines is overstated (M1) |
| 3. Proposition 2 | Correct upper bound; some conditions not stated (m6) |
| 4. Theorem 2 | Correct; constants and numbers right. The summary phrase "whatever happened before" claims a conditional bound that is false (M2) |
| 5. Proposition 3 | Correct; the stream is in the martingale-difference null; 1/B = 11.06% is a valid lower bound on the limit |
| 6. Theorem D | Correct; the hypothesis no longer involves the bets; part 4 is correct with the right events. Consistency gaps in 6.3, 6.5 and 6.7 (m2, m3, m5, m7) |
| 7. Numbers | All spot checks match, except the remark's 10-to-20-day range (m4) |
| 8. Overclaims | M1, M2, and Minor m1, m8, m9, m10 |

## 2. Findings ranked by severity

### Critical

None.

### Major

**M1. The per-day guarantee is attributed to `PeriodicVolatilityMonitor`, which has no mode in which Theorem 1 applies as stated; the extension to the confidence-sequence baseline is stated more strongly than it holds and is not proved.**

*Locations.* `docs/math.md` line 124 (last sentence of the remark after Theorem B: "restarting the statistics at every day boundary (`PeriodicVolatilityMonitor(period=m)`) gives at most `m / A` for every day under the null of Theorem B"); `docs/math.md` line 281 (end of "Choosing A and δ": "with `PeriodicVolatilityMonitor(period=m)` ... it holds for every day of a segment"); `docs/theory_per_period.md` line 50 ("With a confidence-sequence baseline ... the bound becomes `δ + m / A` for every period of a segment, since the event that the confidence bound fails is one event per segment"); `docs/report_parts/TH6.md` line 134 (Reading) and line 138 (Caveats, first bullet); `src/avmonitor/periodic.py` lines 348 to 351 (docstring: "With a known baseline, ... at most period / A for every period, whatever happened before (Theorem 1); a confidence-sequence baseline adds delta once per segment"). Outside the files under review, README line 27 says "Restarting the statistics at every day boundary (`PeriodicVolatilityMonitor`) gives a proved 5% for every day, given everything before it" with no qualification on the baseline.

*Problem.* Theorem 1 is proved for a known baseline m₀ (the null of Theorem B). `PeriodicVolatilityMonitor` inherits `VolatilityMonitor`, whose only baselines are `"plugin"` (the default) and `"cs"` (`monitor.py` lines 121 and 144); there is no known-baseline mode. So:

* with the default plug-in baseline nothing per day is proved, and Proposition A.6 of `docs/theory_running_tau.md` says no bound of the form C m / A can hold; a daily restart caps the accumulation of the statistic but not the persistent drift that the estimation error creates;
* with `baseline="cs"` the correct statement is weaker than "δ + m / A for every period, given everything before it". Let G be the event that the confidence bounds hold at the end of calibration T (Theorem 6.4, proof). G is in F_T, hence in F_{b_i} for every later boundary, so the argument of Theorem 1 gives, under (H↑) and (H↓) with a horizon h that reaches the end of period i,
  `P(alarm in period i | F_{b_i}) ≤ 1{Gᶜ} + 1{G} m / A`,
  hence the **marginal** bound `P(alarm in period i) ≤ δ + m / A`, and for the hazard only `P(alarm in period i | no alarm before b_i) ≤ P(Gᶜ | no alarm before b_i) + m / A`, which is not bounded by δ + m / A in general. Two hypotheses are missing from the one-line claim: the tolerance must hold with h ≥ (end of the period) - T, not with the h = m of Theorem 6.4 part 1, and the statement is marginal, not conditional on the past. None of this is written as a proof anywhere.

This matters because the point of M2 is exactly the difference between a conditional per-day bound and a marginal or first-day one. A reader of the remark in `docs/math.md` or of the README will believe that turning on `period=m` in the shipped monitor buys a proved 5% for every day given the past; it does not, for either baseline the class offers.

*Fix.* (i) In `docs/math.md` line 124, line 281 and README line 27, attach the guarantee to the construction, not to the class: "restarting the statistics at every day boundary gives at most m / A for every day, given everything before it, when the baseline is known (Theorem 1 of `docs/theory_per_period.md`); `PeriodicVolatilityMonitor(period=m, baseline="cs")` gives P(false alarm on a given day of a segment) ≤ δ + m / A when (H↑) and (H↓) hold up to the end of that day; with the default plug-in baseline no per-day bound is proved". (ii) Add the confidence-sequence version as a stated corollary of Theorem 1 with the three-line proof above. (iii) In `periodic.py` lines 348 to 351, replace "With a known baseline" by "For a known baseline (as in `run_engine`; this class itself estimates the baseline)" and spell out the marginal δ + m / A. Optionally (iv): add a known-baseline option (`m0=`, `tau=`) to the class so that Theorem 1 applies to it verbatim.

**M2. The summary of Theorem 2 claims a bound "whatever happened before"; Theorem 2 gives a marginal window bound, and the conditional probability is not bounded by it, even for i.i.d. returns.**

*Location.* `docs/theory_per_period.md` line 11 (summary item 3): "Hence one bet has `P(an alarm in any window of m steps) ≤ m / A`, whatever happened before."

*Problem.* Theorem 2 (a) to (d) bound the unconditional probability of an alarm in a window (s, s + m]. Conditionally on F_s they fail even for i.i.d. returns. The conditional probability is the quantity that Lemma 2 bounds, and it depends on the state at s: if the increase side's average `Σ_k w_k max(C_{k,s}, 1)` is just below 2A, any draw slightly above Y = 1 raises an alarm at s + 1, which for i.i.d. Gaussian returns has probability close to `P(Y > 1) = P(|Z| > 0.9975) = 32%` (for the two-point law {0, B} it is 1/B = 11%). Such states have positive probability at every s past the first few dozen steps (Y has a continuous part, so the side average has a density near 2A without having crossed it). The only conditional statements proved are Lemma 2 (which depends on the state at s) and the one-sided hazard of Theorem 2(e), which conditions on no earlier alarm, not on everything before. Theorem 1, by contrast, is genuinely conditional on F_{b_i}; the summary blurs that difference between the two constructions, which is the distinction M2 asks the authors to make.

*Fix.* Replace "whatever happened before" by "for every window, wherever it starts (an unconditional bound; for one side and i.i.d. returns also the hazard given no earlier alarm, Theorem 2(e))". The parallel phrase on line 46 ("whatever happened on earlier days") is correct, because Theorem 1 is conditional; keep it.

### Minor

**m1. Sharpness of "d m / A" is overstated.** `docs/math.md` line 124 ("and nothing better holds over the null class") and `docs/theory_per_period.md` line 9 ("It says nothing better about day d alone. This is not slack in the proof."). Proposition 1(b) shows that no bound below 1 holds *uniformly in d* (probability 1 on day 40, where d m / A = 2). It does not show that d m / A is attained, or approached, for a given d: the largest probabilities found are 28.0% (two-point day: 31.3%) at d = 10 against d m / A = 50%, and 65.1% (70.1%) at d = 20 against 100%. *Fix:* "and no bound below 1 holds for all days over the null class: ...; for a single day d the largest probability found is about half of d m / A (28% against 50% at d = 10)".

**m2. Section 6.5 recommends a remedy that works against the restated hypothesis.** `docs/math.md` line 295: "A longer calibration window that spans several clustering cycles also makes the calibration average closer to the long-run level." Under the new main hypothesis, (H↑) compares with the *minimum* and (H↓) with the *maximum* of μ_i over the calibration window. For a fixed monitoring window, enlarging the calibration window can only lower the minimum and raise the maximum, so it makes (H↑) and (H↓) harder to satisfy, not easier; TH6's own table agrees (median band needed 0.90 to 2.58 at 390 steps, 2.05 to 5.01 at 1440). The sentence is true only for the bet-weighted version of the remark. *Fix:* "In the bet-weighted version of the remark after Theorem 6.4, a longer calibration window ... long-run level; under the tolerance (H↑), (H↓) a longer window can only widen the band needed".

**m3. Section 6.7 was not updated.** `docs/math.md` line 308: "Proved here: ... Theorem 6.4, parts 1 to 3". Part 4 is new and proved. *Fix:* "parts 1 to 4".

**m4. Wrong range in the remark after Theorem B.** `docs/math.md` line 124: "rose from 0.98% on day 1 to 4.4% to 4.9% on days 10 to 20 (`results/th6_hazard_390.csv`)". In that file the Gaussian hazards on days 10 to 20 run from **4.34% (day 12) to 5.18% (day 19)**, with 5.10% on day 17 (standard errors 0.12 to 0.16 points). The quoted range 4.4% to 4.9% is the 1440-step one (4.45% to 4.94%, `results/th6_hazard_1440.csv`), whose day 1 is 2.04%, not 0.98%. So the sentence mixes the two tables, and in doing so hides that two i.i.d. days have point estimates above 5% (neither significantly, as `docs/theory_per_period.md` line 142 correctly says). My independent run (section 4) also has two days above 5% (5.09% and 5.07%). *Fix:* "rose from 0.98% on day 1 to 4.70% on average over days 11 to 20 (single days 4.3% to 5.2%, none significantly above 5%)".

**m5. Part 4 of Theorem D is "with no assumption on the law" except two.** `docs/math.md` line 244. The proof reuses Step 2, which handles missing returns only under the missingness assumption of 6.7, and the statement needs T < ∞ (or to be read on {T < ∞}). *Fix:* "With no assumption on the conditional means (missing returns as in 6.7): ...". Optionally note that the proof gives the slightly sharper `P(alarm within m) ≤ P(Gᶜ ∪ V_m) + m / A`, and even `P(ρ ≤ T + m) + m / A` without G, where ρ is the first t > T with μ_t outside [m↓, m↑]; that event involves the algorithm, so the stated form is the right headline.

**m6. Proposition 2: conditions that the proof uses but the statement does not state.** `docs/theory_per_period.md` lines 52 to 60. The renewal step needs (i) a restart after an alarm that resets the statistics to 0 with the same baseline and bets, with no blind period (so before the change the statistics are at least 0 and running; the streaming monitor's recalibration after an alarm breaks this); (ii) every later period of length exactly m (Theorem 1 allows periods of at most m, and with shorter periods `p` and `E[min(D₀, m)]` change from period to period); (iii) fixed bets. All three hold in the engine and in the experiment, so the measured comparison is valid. *Fix:* add them to the statement ("for the detector of Theorem 1 with a known baseline, fixed bets, a reset to 0 after each alarm and periods of exactly m steps").

**m7. "What the tolerance excludes" reports a per-path band; parts 1 to 3 need the hypothesis almost surely.** `docs/math.md` line 237 and `docs/report_parts/TH6.md` section Q-TH6f. The quantity computed (`experiments/th6_per_period.py`, `part_tolerance`) is, for each path, the smallest b for which that path satisfies the inequalities; its median is what the text reports as "the band needed". Parts 1 to 3 of Theorem 6.4 require (H↑) and (H↓) almost surely, so the band they need is the essential supremum over paths. The table itself settles this for every admissible band: b must be below 1, and the share of paths that need b > 1 is 44% to 100% in the GARCH cells (`p_holds_b1` = 0 to 0.556), so the inequalities fail with positive probability for every b < 1. (Without the restriction b < 1 the essential supremum is still infinite, because volatility near the GARCH floor during calibration and large during monitoring, with a large τ from the first quarter of the window, all have positive probability.) So parts 1 to 3 never apply to these models; the per-path numbers are what enters part 4, through `P(V_m) = 1 - P(path within band)`, which is 1 at b = 0.25 for every cell. The conclusion drawn ("effectively excludes GARCH clustering", TH6.md line 126) is right, and in fact stronger than stated. *Fix:* one sentence in both places: "Parts 1 to 3 need the tolerance almost surely, which no band achieves under GARCH; the per-path bands above are what part 4 charges through P(V_m), which is 1 at b = 0.25 in every cell."

**m8. "At the same proved level of 5% per day".** `docs/theory_per_period.md` line 13. The plain detector has no proved per-day level; that is the point of section 2. TH6.md line 14 words it correctly ("At a proved 5% for every day, the restarted SR detector's mean delay was ..."). *Fix:* "With the same A (a proved 5% for the first day for the plain detector, for every day for the restarted one), ...".

**m9. Bold play called "not artificial".** `docs/theory_per_period.md` line 142: "So 'at most 5% on every day' is false for the plain detector on a null stream that is not artificial." Bold play reads the detector's internal state at every step; it is a legitimate null stream but a deliberately adversarial one. The conclusion does not need the adjective, since Proposition 1(b) and the calm-then-Gaussian family already prove it. *Fix:* "on an adaptive null stream that watches the detector (and, with certainty, on the constant stream of Proposition 1(b))".

**m10. Mei paragraph quotes the referee's figure, not TH6's.** `docs/theory_per_period.md` line 217: "per-day hazards that rise to 4.6% under i.i.d. Gaussian returns". TH6's own table has 4.70% on average over days 11 to 20 and a largest single day of 5.18%. *Fix:* "that average 4.7% over days 11 to 20".

**m11. Theorem 2 needs fixed bets; the remark in `docs/math.md` does not say so.** Line 124: "the CUSUM statistic gives at most ... for every window of m steps when the returns are i.i.d. (or exchangeable)". The time reversal fails with the adaptive bet (`docs/theory_per_period.md`, remark 3 after Theorem 2). *Fix:* "... when the returns are i.i.d. (or exchangeable) and the bets are fixed (the default grid)".

**m12. Small items in numbers and wording.**
* `docs/report_parts/TH6.md` line 13: "At full scale it passes the i.i.d. bound on day 17" holds for 390-step days; with 1440-step days it passes on day 16 (5.30% against 5.044%, `results/th6_park_strike.csv`). *Fix:* "on day 17 (390) and day 16 (1440)".
* The pooled standard error of the bold-play hazard over days 11 to 20 is 0.058 points from the CSV (`docs/theory_per_period.md` line 142 and TH6.md line 11 say "about 0.05"). *Fix:* "about 0.06".
* `experiments/th6_per_period.py`, `park_strike`: `parked_before_day` uses `parked_at < (d - 1) * day`, but a path whose state at the start of day d already allows a landing gets `parked_at == (d - 1) * day` in the hunting loop, so it is counted as not parked. The descriptive share is low by about one step in 390; nothing else depends on it. *Fix:* `<=`.
* Theorem 1, last sentence (`docs/theory_per_period.md` line 46): for a window that straddles one boundary the argument of Proposition 1(a) gives `(2m + κ) / A` (the restart only lowers N), slightly better than `2 (m + κ) / A` when κ > 0. Optional.

## 3. Item by item

### Item 1. Proposition 1

**Verdict: correct** (with m1 on wording).

*(a)* For s ≤ t ≤ s + m, `S_t = N_t + (s + m - t)` with `N = (M↑ + M↓) / 2`. Under the two-sided martingale-difference null each side satisfies `E[M_t | F_{t-1}] ≤ M_{t-1} + 1`, so S is a nonnegative supermartingale with `S_s = N_s + m`; conditional Ville gives `(N_s + m) / A`; an alarm needs a side at 2A, hence N ≥ A. With restarts after alarms N only drops, so `E[N_s] ≤ s` still holds and the day-d bound `d m / A` follows. Theorem B indeed gives nothing better than `P(T ≤ d m) ≤ d m / A` for day d, by the same argument.

*(b)* With Y ≡ 1 every increase factor is `1 + λ (1 - 1) = 1` and every decrease factor `1 + φ (1 - 1) = 1`, so all 11 Shiryaev-Roberts statistics equal t, and both side averages equal t. The two-sided rule alarms at the first t ≥ 2A: t = 15,600 = 40 × 390 at A = 7800 and t = 57,600 = 40 × 1440 at A = 28,800, the last step of day 40 in both cases, on both sides at once. The law is in the null on both sides (E[Y | F] = 1), and `X = m₀` is attainable because `m₀ = E[min(r², τ²)] ≤ τ²`; it is also self-consistent (|r| = √m₀ has median |r| = √m₀, so a τ computed from these returns would be 3 × 1.4826 × √m₀ ≥ √m₀). With a reset after the alarm the next alarm is on day 80, so p₄₀ = h₄₀ = q₄₀ = e₄₀ = 1. I checked that the engine's dot product with weights 1/6 and 1/5 returns exactly 15,600 and 57,600 at t = 2A (no floating point shortfall), and `test_constant_null_law_alarms_at_2A` confirms the step at a small scale.

### Item 2. Theorem 1 (restart at every day boundary)

**Verdict: correct.** The extension to estimated baselines is the subject of M1.

*What κ is.* κ is an upper bound on each side's average at a boundary, so `N_{b_i} ≤ κ`. The default restart sets every bet's statistic to 0: in the engine `Engine.boundary` does `S[:] = 0.0` when `carry <= 0` (`periodic.py` line 192, default `carry: float = 0.0` on line 157), and in the monitor `_boundary` uses `np.zeros_like` when `carry <= 0` (line 385, default `carry=0.0` on line 363). So **κ = 0 and the bound is m / A** for the default. The variant "carry m at 2A" in the delay experiment has κ = 390 and threshold A' = 2A = 15,600, so `(m + κ) / A' = 780 / 15,600 = 5%`, as stated.

*Conditional on F_{b_i}.* The boundaries are stopping times. Because `E[N_t | F_{t-1}] ≤ N_{t-1} + 1` holds at every deterministic t, it holds on each event {b_i = s} ∈ F_s, so `(S_{b_i + k})_{0 ≤ k ≤ m}` with `S_t = N_t + (b_i + m - t)` is a nonnegative supermartingale for the filtration `(F_{b_i + k})`. Ville's inequality in its conditional form gives `P(max S ≥ A | F_{b_i}) ≤ S_{b_i} / A ≤ (κ + m) / A` almost surely. Until the first alarm in the period the statistics with and without restarts after alarms coincide, and an alarm implies N ≥ A, hence S ≥ A. The proof in the note compresses this to one sentence; it is the standard argument and it is valid. Since {T > b_i} ∈ F_{b_i}, the hazard bound h_d ≤ α follows, as stated.

*Expected number of alarms.* With W_t = (sum of N at the alarms in (b_i, t]) + N_t - (t - b_i): at an alarm step the pre-reset value of N is banked and N restarts at 0, so `W_t - W_{t-1} = N*_t - N_{t-1} - 1`, whose conditional mean is at most 0. Optional stopping at `b_{i+1} ≤ b_i + m` gives `A E[count | F_{b_i}] - m ≤ κ`. Correct.

*Implementation.* Engine: `boundary(t)` runs before step t (0-based) when `t % period == 0` and t > 0, so periods are the 0-based steps `[k m, (k + 1) m)`, exactly the days counted by `di = t // day` in `run_engine`. Monitor: `self.t` starts at -1 and is incremented before each return, so it is 0-based too, and the restart happens before returns 0, m, 2m, ... (the first minute of each day when the stream starts at a session open, as the docstring says); the restart is skipped at `t == start`, where the statistics are already 0. No off-by-one.

### Item 3. Proposition 2 (renewal bound on the delay)

**Verdict: correct upper bound, conditions partly unstated** (m6).

The recursion `R_t = f_t (R_{t-1} + 1)` with `f_t ≥ 0` is nondecreasing in `R_{t-1}`, for each bet on each side, and the alarm rule is monotone in the side averages, so on the same data the restarted (and the plain) detector alarms no later than the detector started at 0 at the change, until the next restart. On {D₀ ≤ m - k} this gives `D_r ≤ D₀`. Otherwise `D_r ≤ (m - k) + D'`, where D' is the delay from the next boundary; D' depends only on data after that boundary, so it is independent of the event {D₀ > m - k}, and the geometric count of alarm-free periods gives `E[D'] = (1 - p) m / p + E[D₀ | D₀ ≤ m] = E[min(D₀, m)] / p`. Adding, `E[D_r] ≤ E[min(D₀, m - k)] + P(D₀ > m - k) E[min(D₀, m)] / p`. At k = 0 the restart coincides with the change and `D_r = D'` exactly, so equality holds; the general formula reduces to the same value at k = 0. `part_renewal` implements this formula faithfully (with `d0 = first + 1`, `min(d0, day)` and the uniform average over k = 0, ..., m - 1, which matches the change positions `nu = 2 day + U{0, ..., day - 1}` of `part_delay`). Every measured mean delay is below its bound (for example 232.6 against 274.9 for sd × 0.8 at 390 steps). The labelling of the renewal comparison as post hoc is accurate.

### Item 4. Theorem 2 (e-CUSUM)

**Verdict: correct; constants and numbers right.** The summary phrase is M2.

*Lemma 1 (time reversal).* `C_s = max_{1 ≤ j ≤ s} Π_{i=j}^{s} f_i` (checked against the recursion `C_t = f_t max(C_{t-1}, 1)`, `C_0 = 0`). With `u_i = v_{s+1-i}`, `Π_{i=j}^{s} v_i = Π_{i'=1}^{s+1-j} u_{i'}`, so C_s is the running maximum of the forward partial products of the reversed sequence. The factors are a fixed function of Y_i (fixed bets), so exchangeability of (Y_1, ..., Y_s) gives the identity in law, jointly over the bets. The forward wealth is a nonnegative supermartingale with initial value 1 under the null (the conditional mean bound for F_{t-1} implies it for the Y-filtration by the tower property), and Ville gives `P(C_s ≥ x) ≤ 1 / x`. **Valid under the exchangeable null with fixed bets.** For an infinite exchangeable sequence the added requirement E[Y_t | F_{t-1}] = 1 means a mixture of i.i.d. laws each with mean 1, which matches the note's example.

*Lemma 2.* I rederived it. `C_{k,t} ≤ f_{k,t} (C_{k,t-1} + 1) ≤ Z_{k,t}` by induction, `E[Z̄_{s+1} | F_s] ≤ Ū_s` and `E[Z̄_t | F_{t-1}] ≤ Z̄_{t-1} + 1`, so `S` is a nonnegative supermartingale started at `Ū_s + m - 1`. It holds under martingale differences, also with predictable bets. Correct.

*Lemma 3.* `E[min(1, wU/A)] = ∫₀¹ P(U ≥ vA/w) dv ≤ w/A + (w/A) log(A/w)`. Correct.

*(b).* `min(1, a + c) ≤ c + min(1, a)` for c ≥ 0, then `min(1, Σ a_k) ≤ Σ min(1, a_k)`, then Lemma 3 with `U_k = max(C_{k,s}, 1)` (which has `P(U_k ≥ x) ≤ 1/x` for x ≥ 1 by Lemma 1, and `w_k ≤ 1 ≤ A`): `(m - 1)/A + Σ_k (w_k / A)(1 + log(A / w_k)) = (m + log A + H(w)) / A`. Correct.

*(c) and the constants.* One side at level 2A gives `(m + log 2A + H) / (2A)`; the union over the two sides gives `(m + log 2A + (H(w) + H(v)) / 2) / A`. The default has 6 equal increase weights and 5 equal decrease weights (`DEFAULT_UP_GRID`, `DEFAULT_DOWN_FRACTIONS`), so H(w) = log 6 and H(v) = log 5. **Right.** Recomputed: at (390, 7800), `(390 + log 15,600 + 1.700599) / 7800 = (390 + 9.655090 + 1.700599) / 7800 = 0.0514558`, so **5.146%**; at (1440, 28,800), `(1440 + 10.961278 + 1.700599) / 28,800 = 0.0504396`, so **5.044%**. The thresholds that make the bound exactly 5% are 8027.69 and 29,053.41 (brentq), as stated. The relative cost of `log 2A + H̄` is 2.91% at 390 and 0.88% at 1440, as in remark 1.

*(d)* After a restart at a, the statistic is `max_{a < j ≤ t} Π f ≤ C_t` bet by bet, so every alarm of the restarted monitor is a time at which the unrestarted side average is at or above 2A. Correct.

*(e)* With independent Y's, every `C_{k,t}` on the increase side is nondecreasing in (Y_1, ..., Y_t), and {T > s} is a decreasing event, so Harris's inequality gives the factorisations used; dividing by P(T > s) gives the one-sided hazard bound and, for two sides (where the monotonicity goes in opposite directions), the stated form with `1 / P(T > s)` on the Lemma 3 part. I rederived the two-sided expression: `(m - 1)/A + (1 + log 2A + H̄) / (A P(T > s))`. Correct.

*Remarks.* Remark 1 (the law f ∈ {0, A} needs an unclipped factor with λ = 1, as the note says), remark 2(i) (a single CUSUM of the mixture wealth reverses to the running maximum of one supermartingale, hence m / A) and remark 5 (C ≤ R) are correct. Remark 4 (confidence-sequence baseline with i.i.d. returns: given F_T the monitored draws are i.i.d. with mean on the null side on G) is correct.

### Item 5. Proposition 3 (park and strike)

**Verdict: correct.**

*Membership in the null class.* Hunt draws are on {0, B} with P(B) = 1/B (mean 1); parking plays Y = 1; the strike (and the experiment's bold play) uses two-point draws with mean exactly 1 on [0, B]. Which law is used at step t is decided by the state at t - 1 and by t, so it is predictable, and `E[Y_t | F_{t-1}] = 1` at every step. The three values correspond to |r| = 0, |r| = √m₀ and |r| ≥ τ. So the stream is **in the martingale-difference null class** of Theorem B (band 0), on both sides. It is not exchangeable, which is why Theorem 2 does not apply.

*No alarm before s + 1.* Before ρ no legal draw reaches 2A on either side (increase statistics are increasing in Y, decrease statistics decreasing). At ρ, because `f_k(B) ≥ 1`, `Σ_k w_k max(C_{k,ρ}, 1) ≤ Σ_k w_k f_k(B) max(C_{k,ρ-1}, 1) < 2A`, and likewise on the decrease side; parking replaces C by max(C, 1) and then holds it. I checked each inequality. The experiment enforces it (`park_strike` raises an error on any alarm during the hunt) and none occurred.

*ρ finite, strike probability.* A run of n* draws of B forces the trigger; such runs occur almost surely. If the increase side triggered, Y = B (probability 1/B) alarms; if the decrease side triggered, Y = 0 (probability 1 - 1/B ≥ 1/B, since B = 9.045 ≥ 2) alarms; if both, either draw does. So `P(T = s + 1) ≥ P(ρ ≤ s) / B → 1 / B`.

*The 11% limit.* `B = τ² / m₀ = 9 / 0.995007 = 9.04516` for Gaussian clipping at c = 3 (τ = 3 exactly, because the MAD scale is 1/Φ⁻¹(0.75)), and `1 / B = 0.11056`. **1/B ≈ 11% is right as a lower bound** on the limit; the actual limit of this stream is a mixture of 1/B and 1 - 1/B according to which side parks, so it is larger. The note and TH6.md say "at least", which is correct. In the scaled run the measured 27.0% at period 120 is above both the parked share times 1/B (6.3%) and 1/B itself, and below the mean of Lemma 2's bound (32.9%), as it must be.

### Item 6. Theorem D (Theorem 6.4 restated)

**Verdict: correct, with consistency gaps** (m2, m3, m5, m7).

*Does the hypothesis depend on the algorithm?* (H↑) `sup_{T<t≤T+h} μ_t ≤ (1 + b) min_{n_τ<i≤T} μ_i` and (H↓) involve the law of the returns, τ (a function of the first n_τ returns, through c and the split fraction), T and b. They do **not** involve the confidence-sequence bets, which was the referee's complaint, so the restatement does what M4 asked. "Algorithm-free" is accurate in that sense only: τ, n_τ and T are still algorithm settings, and τ is random, so the hypothesis must hold for the realised τ (almost surely). The note says this. With the width-target stopping rule, T depends on the sequence and hence on its bets; the sufficient form with `T_min < t ≤ T_max + h` and `n_τ < i ≤ T_max` removes this, and I checked that it implies (H↑) and (H↓) for every T in [T_min, T_max].

*Does the proof go through?* Yes. `μ̄↑_T = Σ λ↑_i μ_i / Σ λ↑_i` is a convex combination with positive (predictable) weights over the observed returns in (n_τ, T], so `min ≤ μ̄↑_T` and `μ̄↓_T ≤ max` pathwise; the minimum over all of (n_τ, T] (including steps with a missing return) is at most the minimum over the observed ones, so the inequality chain `μ_t ≤ (1 + b) min_i μ_i ≤ (1 + b) μ̄↑_T ≤ (1 + b) Ū_T = m↑` holds on G, and the decrease side likewise with the strict inequality `L̄_T < μ̄↓_T` on G. `{T < t ≤ T + h}` is in F_{t-1}, so Step 3 with h = m and Part 2 with h = ∞ go through unchanged.

*Part 4.* ρ (first t > T with μ_t > m↑, or μ_t < m↓ when the decrease side is on) is predictable because μ_t is F_{t-1}-measurable. M' (equal to M before ρ, then increasing by 1 per step) satisfies the one-step bound at every step; on G ∩ {ρ > T + m}, M' = M on the horizon; and on G, ρ ≤ T + m implies that (H↑) or (H↓) fails at ρ, because `μ_ρ > m↑ ≥ (1 + b) μ̄↑_T ≥ (1 + b) min_i μ_i` (and `μ_ρ < m↓ < (1 - b) μ̄↓_T ≤ (1 - b) max_i μ_i`). Hence `P(alarm within m) ≤ P(Gᶜ) + P(V_m) + m / A`. **Correct, and stated with the right events**: V_m is the failure of the tolerance at some t ∈ (T, T + m] relative to the calibration window, which is the event the proof needs. Including (H↓) failures when the decrease side is off only enlarges V_m, so the bound stays valid. Caveat m5 (missing returns, T < ∞).

*Edits elsewhere.* Corollary 6.5: the sup, min, inf and max all equal m(τ) under a constant conditional law, so both hypotheses hold with b = 0 for every horizon; the "more generally" part needs σ_t constant on (n_τ, T], as written. Consistent. Corollary 6.6: requiring h = ∞ in each segment is sufficient (only the monitored steps up to that segment's alarm matter) and its proof cites Step 2, which now carries h; consistent. Section 6.5, bullets 1, 2 and 4: consistent. Bullet 3: inconsistent (m2). Section 6.7: not updated for part 4 (m3). "What the tolerance excludes": numbers right, framing per path versus almost surely (m7).

### Item 7. Numbers

**Verdict: all spot checks match their CSV files, except m4.** In the table, "CSV" is my own pandas computation from the committed file; "independent" is my own simulator with new seeds (section 4).

| Claim (location) | Stated | From the CSV | Independent |
|---|---|---|---|
| Gaussian hazard, day 1, 390 (`th6_hazard_390.csv`) | 0.98% | 0.978% (se 0.049) | 1.085% (se 0.073, 20,000 paths) |
| Gaussian hazard, mean of days 11 to 20, 390 | 4.70% | 4.700% | 4.730% |
| Gaussian, days 2 to 5 and 6 to 10 | 3.02, 4.13 | 3.024, 4.127 | 3.074, 3.984 |
| Gaussian, single days 10 to 20 (math.md remark) | 4.4% to 4.9% | **4.34% to 5.18%** | 4.38% to 5.09% |
| Bold play, mean of days 11 to 20 (`th6_adversary.csv`) | 5.25% | 5.246% (pooled se 0.058) | 5.29% (15,000 paths) |
| Bold play, day 1 | 2.53% | 2.527% | 2.57% |
| Family, d = 5, 10, 20 (`th6_calm_family.csv`) | 11.3, 28.0, 65.1% | 11.33, 28.02, 65.07% | 11.36, 28.16, 64.86% (40,000 paths each) |
| Restarted, largest per-day probability or hazard, 390 / 1440 | 2.75%, 3.52% | hazard 2.751% (day 20), 3.516%; probability 2.663%, 3.517% (all bold play) | bold play against the restarted detector, days 1 to 6: 2.27% to 2.67% |
| Delay ratio restarted / plain, 390 (`th6_delay.csv`) | 1.06 to 1.44 | 1.056 to 1.437 | 1.056 (×2), 1.097 (×1.5), 1.436 (×0.8), 8000 paths |
| Delay ratio restarted / plain, 1440 | 1.02 to 1.14 | 1.018 to 1.137 | not rerun |
| e-CUSUM certified / plain | 1.30 to 1.83, 1.24 to 1.66 | 1.303 to 1.827, 1.243 to 1.656 | not rerun |
| Plain faster than a fresh start (`th6_delay_renewal.csv`) | 14% to 31% | 13.9% to 30.6% | |
| GARCH median band, both sides (`th6_tolerance.csv`) | 0.90 to 5.01 | 0.903 to 5.014 (390: 0.90 to 2.58; 1440: 2.05 to 5.01) | 390: 0.93 to 2.58; 1440: 2.05 to 4.87 (1000 paths) |
| No path within b = 0.25; at most 1.8% within b = 0.5 | as stated | 0 in every cell; 1.8% | 0; 1.3% |
| Decrease side alone, median band | 0.37 to 0.82 | 0.367 to 0.823 | |
| Park and strike, day 17 and day 20 (`th6_park_strike.csv`) | 5.43%, 6.16% (2.2 and 7.3 se) | 5.430%, 6.157% (2.17, 7.29 se) | |
| Corollary C, false alarms per day over 20 days | 0.037 to 0.048 | 0.0370 to 0.0479 | |
| e-CUSUM per day, i.i.d., at most | 0.18% | 0.182% | 0.14% (Gaussian, no restart) |
| Window bound, certified A | 5.146%, 5.044%; 8027.7, 29,053.4 | | 0.0514558, 0.0504396; 8027.69, 29,053.41 |

The differences between the CSV and my runs are within sampling error: where I computed standard errors the largest gap is 1.75 combined standard errors (the mean hazard over days 6 to 10, 3.98% against 4.13%), and the day-1 hazard differs by 1.2 (1.085% against 0.978%); the tolerance medians come from 1000 paths against 2000 and differ by at most 0.14. I also checked every row of the four tables in `docs/theory_per_period.md` section 5 and of TH6.md Q-TH6a to Q-TH6f against the CSVs; apart from the items in m4 and m12 they match to the printed precision.

*Simulation engine, read line by line.* No bug that affects a reported number.

* Recursions: SR adds 1 then multiplies by the factor; CUSUM takes max(·, 1) then multiplies. Both match `detectors.betting_edetector`, which `test_engine_matches_package_detector` checks path by path for both kinds.
* Alarm side threshold: every alarm and every bound uses 2A per side (`Engine.step`, `landing`, `PeriodicVolatilityMonitor.update`, `conditional_window_bound`, `window_bound_two_sided`). Correct.
* Restart timing and day boundaries: see item 2. In `part_family` the restarted detector is reset by `boundary(t0)` at t0 = (d - 1) m and the SR state `(d - 1) m` and CUSUM state 1 after the calm days are exact.
* Delays: `delay = first_after - nu + 1` counts the alarm step; changes at `2 day + U{0, ..., day - 1}` are a uniform step of day 3 (0-based); `detected_within_horizon` is 1 in every cell, so no censoring.
* Laws: τ, m₀ and B in `th6_laws.csv` agree with closed forms (Gaussian τ = 3, m₀ = 0.995007; jump mixture with jump sd 10; ticks of 0.5 with median |r| = 0.5 and τ = 2.2239).
* Bold play and park and strike: the landing draws solve the linear equation for the side average exactly; two-point draws have mean 1 (`test_two_point_draws_have_mean_one`). The only defect is the cosmetic off-by-one in `parked_before_day` (m12).

### Item 8. Overclaims

**Verdict:** M1 and M2 above; Minor m1, m8, m9 and m10. Things that are labelled correctly and need no change: every proposition labelled [PROVED] has a complete proof; the pooling of bold play over days 11 to 20 and the renewal comparison are labelled post hoc in both TH6.md and the theory note; the stochastic volatility tolerance is labelled as computed with the volatility path in the filtration; all simulated results carry [SYNTHETIC] or [NUMERICAL]; TH6.md's Caveats say that everything is for the true baseline (but see M1 for the sentence that follows). The README sentence in M1 is outside the files under review but repeats the same overclaim without the caveat.

## 4. Independent recomputation (details)

`r8th6/indep_check.py`: τ = 3, m₀ = 0.9950072780, B = 9.0451599689, 1/B = 0.1105564. Two-sided rule, 6 increase bets (0.02 to 0.5) and 5 decrease bets (0.1 to 0.95 of 1/(B - 1)), each side averaged and compared with 2A, reset to 0 after each alarm.

* `family 40000` (seeds 1005, 1010, 1020): P(alarm on day d) after d - 1 calm days with every SR bet at (d - 1) m: **11.36%, 28.16%, 64.86%** for d = 5, 10, 20.
* `hazard 20000` (seed 777), 20 Gaussian days of 390 steps, A = 7800: hazard by day (percent) 1.085, 2.411, 3.015, 3.418, 3.451, 3.643, 3.739, 4.063, 4.100, 4.376, 4.604, 4.378, 4.578, 4.902, 4.723, 5.090, 4.857, 4.543, 5.069, 4.557. The e-CUSUM without restarts crossed 2A on some step of a day with probability at most 0.14%; the detector restarted at every day boundary had 0.84% to 1.13% per day.
* `delay 8000`: mean delays plain / restarted for a change at a uniform step of day 3: sd × 2: 20.19 / 21.31; × 1.5: 43.43 / 47.63; × 0.8: 161.94 / 232.51.

`r8th6/indep_check2.py`:

* `bold 15000` (seed 4242): bold play against the plain detector, 20 days: day 1 2.57%, mean of days 11 to 20 5.29%, largest 5.60%. Against the detector restarted daily (seed 4343, 6 days): 2.27% to 2.67%.
* `tol`: GARCH(1,1) cells of E9 (1000 paths, burn-in 4000), τ from the first quarter of the window, exact μ_t: medians (increase, decrease, both) at 390: (1.44, 0.59, 1.44), (2.58, 0.72, 2.58), (1.64, 0.61, 1.64), (0.93, 0.49, 0.93); at 1440: (2.05, 0.67, 2.05), (4.87, 0.83, 4.87), (4.11, 0.82, 4.11), (2.28, 0.75, 2.28). No path within b = 0.25; at most 1.3% within b = 0.5.

## 5. References

I could not open full texts from this session; arxiv.org is blocked and I did not try to get around it. The bibliographic details of Lai (1998, IEEE Transactions on Information Theory 44(7), 2917 to 2929), Mei (2008, Sequential Analysis 27(4), 354 to 376), Pergamenchtchikov and Tartakovsky (2018, Statistical Inference for Stochastic Processes 21(1), 217 to 259), Harris (1960, Proceedings of the Cambridge Philosophical Society 56, 13 to 20) and Shin, Ramdas and Rinaldo (2024, New England Journal of Statistics in Data Science 2(2), 229 to 260) agree with my own knowledge of these papers, but **I could not verify the content claims** made about them (the exact form of Lai's window criterion and conditions on m_α, Mei's definitions, whether the criterion of Pergamenchtchikov and Tartakovsky conditions on no earlier alarm). I could not verify the page range given for Shekhar and Ramdas (2024, PMLR 235). The note's own verification paragraph (`docs/theory_per_period.md` line 219) says the same, so nothing is overclaimed there. The form of Harris's inequality used in Theorem 2(e) (nondecreasing and nonincreasing functions of independent real coordinates) is the standard product-measure form and needs no citation beyond what the note gives.

## 6. Summary of fixes

| # | Severity | Location | Fix |
|---|---|---|---|
| M1 | Major | `docs/math.md` 124 and 281; `docs/theory_per_period.md` 50; TH6.md 134, 138; `periodic.py` 348 to 351; README 27 | Attach the per-day bound to the known-baseline construction; state the cs version as marginal `δ + m / A` per day under the tolerance up to the end of that day, with a proof; say the plug-in default has no per-day bound |
| M2 | Major | `docs/theory_per_period.md` 11 | Replace "whatever happened before" by an unconditional statement (plus Theorem 2(e) for the one-sided hazard) |
| m1 | Minor | `docs/math.md` 124; `docs/theory_per_period.md` 9 | "No bound below 1 for all days", not "nothing better" |
| m2 | Minor | `docs/math.md` 295 | Restrict the longer-window sentence to the bet-weighted remark |
| m3 | Minor | `docs/math.md` 308 | "parts 1 to 4" |
| m4 | Minor | `docs/math.md` 124 | 390-step range is 4.3% to 5.2% (mean 4.70% over days 11 to 20) |
| m5 | Minor | `docs/math.md` 244 | Part 4 still needs the missingness assumption and T < ∞ |
| m6 | Minor | `docs/theory_per_period.md` 54 | State reset to 0 with the same baseline, periods of exactly m, fixed bets |
| m7 | Minor | `docs/math.md` 237; TH6.md Q-TH6f | Parts 1 to 3 need the tolerance almost surely; per-path bands enter part 4 only |
| m8 | Minor | `docs/theory_per_period.md` 13 | The plain detector has no proved per-day level |
| m9 | Minor | `docs/theory_per_period.md` 142 | Drop "not artificial" |
| m10 | Minor | `docs/theory_per_period.md` 217 | 4.7%, not 4.6% |
| m11 | Minor | `docs/math.md` 124 | Add "and the bets are fixed" |
| m12 | Minor | TH6.md 13; theory note 142; `th6_per_period.py` `park_strike` | Day 16 for 1440; se about 0.06; `<=` |
