# Independent review, round 5: E16 (adversarial search for null streams)

> *Note added 3 November 2026: the "independent reviewer" or "referee" below was a separate AI agent (Claude, by Anthropic) instructed to act as a sceptical referee. It was not a human expert, and the review has not been confirmed by a human.*

**Reviewer role.** An adversarial reader: a statistics journal referee and a code auditor. I did not write any of the files below. The goal is that every claim about E16 survives an expert reading.

**Scope.**

1. Code: `experiments/e16_adversarial_null.py`, `tests/test_e16_adversarial.py`.
2. Results: `docs/report_parts/E16.md`, `docs/iterations/E16.md`, `results/e16a_*.csv`, `results/e16b_*.csv`, `results/figures/e16_inside_class.png`, `results/figures/e16_outside_class.png`.
3. How E16 is stated elsewhere: `README.md:24` (first "Key results" bullet) and `README.md:146` (limitations list), `PITCH_BAPAT.md:43` (section 7) and `PITCH_BAPAT.md:52` (section 8), `docs/report_parts/_header.md:30` (row E16), `docs/limitations.md:12` ("The bound is conservative"), `docs/theory_delay.md:340` ("For the default mixture, unknown").
4. Statements in `docs/theory_delay.md` (section 8) and `docs/report_parts/E11.md` that E16 changes.

**Read first.** `docs/math.md` section 5 (the null class of Theorem B), `src/avmonitor/detectors.py` and `src/avmonitor/monitor.py`, `docs/theory_delay.md` sections 3 and 8, `docs/report_parts/E9.md` and `E11.md`, and the earlier reviews in `docs/review/`. Their points are not repeated.

**State of the tree.** E16 was committed in 18c12e5, 40561ab and eab03b9. The review is against HEAD 68a9cc5 on branch `claude/project-thread-ytr1da`. Between eab03b9 and 68a9cc5 only E17 lines of `README.md` and `_header.md` changed, so every line number below holds at both commits. Other agents had uncommitted changes in the working tree while I worked (E17 and paper files); none touches a file in scope.

**What I ran.** Scripts are in `/tmp/claude-0/-home-user/e30ce406-b0e8-5b23-b612-d49eb104d0ee/scratchpad/r5e16/`; nothing was written to `results/`. Everything ran under `nice -n 10` with at most 2 processes.

* `python -m pytest -q tests/test_e16_adversarial.py`: **20 passed in 4.3 s.**
* `rerun.py`, exact reruns with the script's own functions, seeds and path counts. **Every number equals its CSV cell:**
  * inside the class, the B scan of bold play at B = 2 (seed 16930, 60,000 paths of 1,950 steps, 23 s): 2.6200% in one day and 17.835% in five, as in `e16a_curves.csv`;
  * outside the class, the GARCH cell α + β = 0.99, α = 0.005 (seed 16242, 9 s): mean drift 0.038843, rates 5.5667% and 49.933%, oracle 0.7833% and 12.300%, as in `e16b_garch.csv`;
  * outside the class, price ticks of 4 sd (seed 16409, 12 s): τ = 17.7912, m0 = 3.1890, conditional means 0.8806 and 1.4797, price-rounded 5.28% and 28.85%, twin 0.84% and 6.58%, as in `e16b_price_ticks.csv`.
* `rerun.py`, fewer paths and new seeds (consistent with the CSVs within Monte Carlo error):
  * bold play at B = 2, 20,000 paths: 2.735% (2.52% to 2.97%) and 17.99% (17.46% to 18.53%);
  * the same GARCH cell, 2,000 paths, with the CSV's τ and m0: 5.85% (4.90% to 6.97%) and 52.1% (49.9% to 54.3%); oracle 0.90% and 13.5%.
* `rerun.py tick4_iid_marginal`: a control for price ticks with the same marginal law (M2).
* `dp_check.py`: E11's dynamic programme (`design.worst_case_false_alarm`) for the single bet 0.05 on finer wealth grids (M3).
* `highB.py`: bold play at B = 100, 150 and 1,000, 20,000 fresh paths each (m1).
* pandas checks of every table cell and every number in the text of `E16.md`, and of the E16 statements in the README, the pitch, the report index, the limitations, `theory_delay.md` and `ITERATION_LOG.md`.

**Severity.**

* *Critical*: an adversary that is not in the null class, a wrong count, or a wrong headline number.
* *Major*: an overclaim, a misleading comparison or control, or a statement elsewhere that E16 contradicts.
* *Minor*: a wrong number outside the headlines, a weak test, wording or clarity.

---

## 1. Summary verdict

**Every adversary called "inside the class" is inside the null class of Theorem B, the false alarms are counted correctly, every headline number traces to its CSV, and three cells rerun exactly. I found no critical problem.** Bold play draws Y_t from a two-point law {y_low, y_high} with 0 ≤ y_low ≤ 1 ≤ y_high ≤ B, both points fixed by the detector's state before the draw and the weight chosen so that the mean is exactly 1. This holds at the landing step too, and the path stops at its first alarm (section 5). The detector is the monitor's default (same grid, same recursion, each side at 2A, increase side checked first). The decomposition of the slack is the exact optional-stopping identity of `theory_delay.md` 8.1, computed correctly.

**The problems are in what the results are said to show.**

* **M1.** "The rest is structural", "why the rest of the slack cannot be removed" and "nearly sharp bet by bet but not for the mixture" claim more than a search can show. The decomposition says where the slack sits under the laws tried, not that no null law can move it. A one-line floor argument proves only 2% of slack at one day.
* **M2.** The price-tick control is not "the same paths with the condition restored", as `E16.md:182`, the iteration notes, the README, the pitch and the report index say. It is a different i.i.d. law with a very different marginal (95% zero returns against 80% at a tick of 4 sd). A control with the same marginal law, which I ran, confirms the conclusion: 0.92% in one day and 7.2% in five days.
* **M3.** `theory_delay.md` 8.3 and 8.4 still present the dynamic programme's 18% for the bet 0.05 as the worst case and use it as evidence for a sharper bound, and `theory_delay.md` 8.3 and `E11.md` Q11d still give 42% as the best adversary for the default mixture. E16 shows that the 18% is a loose lower bound: bold play reaches 26%, and refining E11's wealth grid raises the dynamic programme itself from 0.181 to 0.232. For the two-sided default the best adversary now reaches 52%.

**Most important fixes before the work is shown.**

1. Replace "structural" and "cannot be removed" by "under every law tried", and state the only slack that is proved (M1).
2. Replace the tick twin by i.i.d. draws from the marginal law of the price-rounded returns, or reword every "same paths" and "when it was restored" for ticks (M2).
3. Correct the small-bet and default-mixture statements in `theory_delay.md` 8.3 and 8.4 and in `E11.md` Q11d (M3).

---

## 2. Critical findings

None.

* **Null class.** Every adversary in part a keeps E[Y_t | past] = 1 exactly and Y_t in [0, B] (section 5), including bold play at and after the landing step.
* **Counting.** Two-sided default detector, each side at 2A, horizons in steps, Wilson intervals, fresh seeds: all correct (section 5).
* **Headline numbers.** 2.62% (2.50% to 2.75%), 18.1%, 2.40%, 0.97% and 12.4%, 0.52 and 0.73, 43 of 64 cells, 3.9%, 41% to 42%, 32.5%, 0.59, 4 sd and 0.88 all agree with their CSV cells.

---

## 3. Major findings

### M1. "The rest of the slack is structural" is not shown

* **Where.**
  * `docs/report_parts/E16.md:13`: "**The adversary removes about half of E9's slack, and the rest is structural.** ... 36% by the three smallest increase bets, which cannot build wealth fast enough to reach 2A within a day while the same draws also serve the other bets."
  * `E16.md:84`: "So the bound is nearly sharp bet by bet but not for the mixture of eleven bets on two sides within a day;"
  * `E16.md:180`: "The rest of the slack is structural for the default grid within short horizons:"
  * `experiments/e16_adversarial_null.py:28`: "a3  Why the rest of the slack cannot be removed:"
  * In the same spirit, without a qualifier: `README.md:146` ("the rest is wealth that the default bets leave on paths that never alarm") and `_header.md:30` ("the rest of the slack is wealth stranded in the decrease and small increase bets").
* **Evidence.**
  * The identity `P(σ ≤ m) = (m / A) F1 F2 F3` holds for every exact null law. Its value under bold play (F1 = 0.99, F2 = 0.49, F3 = 1.00, `e16a_slack.csv`) says where the slack sits *for that law*. To say that the rest "cannot be removed" one needs an upper bound on F2 over all null laws, and E16 has none. `E16.md:84` says so itself ("a provably sharper threshold ... would need a proof that F2 stays well below 1 for every null law, which we do not have"), and so does caveat 1 at `E16.md:186`. The summary bullet and the closing reading contradict these sentences.
  * What *can* be proved is small. On a path with no alarm by step m, each statistic is at least its value after m draws that all go against it: `(1 - λ)(1 - (1 - λ)^m) / λ` for an increase bet λ (Y = 0 at every step) and `(1 - φ)(1 - (1 - φ)^m) / φ` for a decrease fraction φ (Y = B at every step). For the default grid at m = 390 this gives `N_m ≥ c = 8.3` on every such path. Optional stopping then gives `A P + c (1 - P) ≤ m`, so `P(false alarm within m) ≤ (m - c) / (A - c) = 4.90%` at one day (0.98 of the bound, F2 ≤ 0.98). So 2% of the slack is provably structural; the remaining factor of about 1.9 between 2.62% and 4.90% is not explained by any proof.
  * "Cannot build wealth fast enough" is also not shown bet by bet. E16's own single-bet check finds that the same adversary takes the bet 0.05 *alone* to 0.26 of its bound in a day (`e16a_single_bet_check.csv`), and the bet 0.02 was not checked.
* **Fix.**
  * `E16.md:13`: replace "**The adversary removes about half of E9's slack, and the rest is structural.**" with "**The adversary removes about half of E9's slack; under every law tried, the rest is wealth left on paths that never alarm.**" Replace "which cannot build wealth fast enough to reach 2A within a day while the same draws also serve the other bets." with "which did not reach 2A within a day under any law tried while the same draws also serve the other bets; whether some null law can move this wealth onto alarmed paths is open."
  * `E16.md:84`: replace "So the bound is nearly sharp bet by bet but not for the mixture of eleven bets on two sides within a day;" with "So the bound is nearly sharp for the single bets 0.5 and 0.2 and the decrease bet at one day or more, while no law tried takes the mixture of eleven bets on two sides above about half of it within a day;". After "which we do not have" add "(the smallest value the two-sided statistic can take on a day without an alarm, c = 8.3, proves only P(false alarm within m) ≤ (m - c) / (A - c), 4.90% at one day)".
  * `E16.md:180`: replace "The rest of the slack is structural for the default grid within short horizons: the decrease bets and the small increase bets did not use their share of the budget within a day under any law we tried, and one draw has to serve all eleven bets on both sides at once." with "Under every law we tried, the rest of the slack is wealth left in the decrease bets and the small increase bets, which did not use their share of the budget within a day while one draw has to serve all eleven bets on both sides at once; we have not shown that no null law can do better."
  * `e16_adversarial_null.py:28`: replace "Why the rest of the slack cannot be removed:" with "Where the rest of the slack sits:".
  * `README.md:146`: replace "the rest is wealth that the default bets leave on paths that never alarm." with "under that stream the rest is wealth that the default bets leave on paths that never alarm, and whether a worse stream exists is open." `_header.md:30`: replace "the rest of the slack is wealth stranded in the decrease and small increase bets." with "under that adversary the rest of the slack is wealth stranded in the decrease and small increase bets."

### M2. The price-tick control is a different law, not the same paths with the condition restored

* **Where.**
  * `E16.md:182`: "each has a control inside the class on the same paths that stays below the bound".
  * `E16.md:176`: "it is the conditional mean condition that fails there, not the tails or the discreteness (the twin takes values on the same grid)."
  * `docs/iterations/E16.md:109`: "and a control inside the class on the same paths."
  * `README.md:24`: "and disappeared when it was restored." `PITCH_BAPAT.md:43`: "and vanished when it was restored". `_header.md:30`: "each failure matches a broken condition and disappears when it is restored".
* **Evidence.**
  * For GARCH and the plug-in the controls are what is claimed: the same X on the same paths, with the baseline replaced by the true conditional mean μ_t (`e16_adversarial_null.py:852` to `854`) or by the true mean at τ-hat (`:945`).
  * For price ticks the "twin" is new draws `z = h round(Z / h)` (`:1062`), a different random stream with a different law. Its marginal law is far from that of the price-rounded returns. At a tick of 4 sd my exact rerun gives 95.5% zero returns for the twin against 80.1% for rounded prices. The twin's B = τ² / m0 is 435 against 99 for rounded prices; at 6 sd it is 7,328 against 149 (`tick_law` and `e16b_price_ticks.csv`). So the twin does not hold "the tails or the discreteness" fixed; it changes the whole marginal law.
  * The twin is also not needed for the attribution. τ and m0 are the true unconditional values of the price-rounded stream, so if E[X_t | past] were m0 at every step, Theorem B would apply. The only condition that can fail is the conditional mean.
  * A control with the same marginal law confirms the conclusion. I drew the monitored returns i.i.d. from the pooled returns of 2,000 price-rounded paths of 10,000 steps at a tick of 4 sd (80.0% zeros; mean of X within 0.11% of the CSV's m0) and used the CSV's τ and m0 (10,000 paths, new seed): **0.92% (0.75% to 1.13%) in one day and 7.2% (6.7% to 7.7%) in five days**, against 5.28% and 28.85% for the price-rounded stream.
* **Fix.**
  * Preferred: in part b3, replace or complement the twin by i.i.d. draws from the marginal law of the price-rounded returns (pooled from the 1,000 long paths already drawn), with the same τ and m0; it costs a few seconds per tick. Then the README and pitch wording can stay, with "(for ticks, the same values drawn independently)" added.
  * In any case, `E16.md:176`: replace "and it is the conditional mean condition that fails there, not the tails or the discreteness (the twin takes values on the same grid)." with "and, since τ and m0 are the true unconditional values, the only condition of Theorem B that can fail is the conditional mean; the twin rounds each return independently and has a different marginal law (95% zero returns at a tick of 4 sd against 80% for rounded prices), so it checks the same grid, not the same values."
  * `E16.md:182`: replace "and each has a control inside the class on the same paths that stays below the bound:" with "and each has a control inside the class that stays below the bound (on the same paths for GARCH and the plug-in; for price ticks a separate i.i.d. law on the same tick grid):".
  * `iterations/E16.md:109`: replace "and a control inside the class on the same paths." with "and a control inside the class (on the same paths for GARCH and the plug-in; a separate i.i.d. law on the same tick grid for price ticks)."
  * If the control is not added: `README.md:24`, replace "and disappeared when it was restored." with "and an in-class control stayed below the bound (for clustering and the estimated baseline, the same paths with the condition restored)." `PITCH_BAPAT.md:43`: replace "and vanished when it was restored:" with "and an in-class control stayed below the bound:". `_header.md:30`: replace "and disappears when it is restored:" with "and an in-class control stays below the bound:".

### M3. `theory_delay.md` and `E11.md` still treat the dynamic programme's small-bet value as the worst case

* **Where.**
  * `docs/theory_delay.md:324`: "Dynamic programming over the whole nonparametric null ... gives the worst-case false alarm probability within m steps as a fraction of `m / A`:" and the table row at `:330` ("increase 0.05 | 0.0002 | 0.181 | 0.694 | 0.716").
  * `theory_delay.md:340`: "The small-bet result (18% at one day for the bet 0.05) suggests that grids with small bets could carry a provably sharper short-horizon bound, but we have no proof." This bullet was edited in the E16 integration commit (eab03b9) and still quotes 18%.
  * `theory_delay.md:335`: "**For the default mixture, at least 42%.**" and `docs/report_parts/E11.md:56`: "for the default mixture the best adversary found reaches 42% of the bound and a certified sharper bound is open."
* **Evidence.**
  * `E16.md:79` and `:84`: bold play, a legitimate adversary, reaches 0.26 of the bound for the bet 0.05 in one day against the dynamic programme's 0.18. On 41,020 fresh paths the interval for the bold-play ratio is 0.241 to 0.285 (`e16a_single_bet_check.csv`: 1.312%, Wilson 1.206% to 1.426%). Since every adversary's value is a lower bound on the worst case, the 0.181 cannot be the worst case.
  * **E16's reading that the dynamic programme's small-bet values are loose lower bounds is right, and the cause is the wealth grid.** The function rounds each transition down to a 700-point log grid from 0.95 to A, whose neighbouring states differ by a factor 1.013. Once the wealth exceeds about 77, the step M → M + 1 (a draw at Y = 1) rounds back down to M. A small bet builds wealth mainly through that step. Refining the grid with the repository's own function (`dp_check.py`), the one-day lower bound for the bet 0.05 rises from **0.181 (700 points) to 0.205, 0.220 and 0.232 (1,400, 2,800 and 5,600 points)**. At five days it rises from 0.694 to **0.773 (2,800 points)**, inside the 95% interval of bold play (0.766, 0.750 to 0.781).
  * So the 18% at `:340` is not evidence about the worst case, and "gives the worst-case false alarm probability" at `:324` is stated as a value, not a bound.
  * `:335` and `E11.md:56` give 42% as the best adversary "for the default mixture". That is the i.i.d. law for the increase grid alone at level A. For the two-sided default detector E16's bold play reaches 52% at one day.
* **Fix.**
  * `theory_delay.md:324`: replace "gives the worst-case false alarm probability within m steps as a fraction of `m / A`:" with "gives a lower bound on the worst-case false alarm probability within m steps, as a fraction of `m / A`:".
  * `theory_delay.md:333`: add at the end of the paragraph: "For the bet 0.05 the 700-point grid makes the lower bound loose (a step M → M + 1 rounds back down to M once M exceeds about 77): refining the grid raises the one-day value to 0.205, 0.220 and 0.232 with 1,400, 2,800 and 5,600 points, and E16's bold play reaches 0.26 at one day, 0.77 at five days and 0.75 at ten days (`results/e16a_single_bet_check.csv`)."
  * `theory_delay.md:335`: replace "**For the default mixture, at least 42%.**" with "**For the default mixture, at least 42% (increase grid alone, i.i.d. law) and 52% (two-sided detector, adaptive law, E16).**"
  * `theory_delay.md:340`: replace "The small-bet result (18% at one day for the bet 0.05) suggests that grids with small bets could carry a provably sharper short-horizon bound, but we have no proof." with "For the bet 0.05 alone the best adversary known reaches 26% of the bound at one day (E16's bold play; the 18% of the table in 8.3 is a loose lower bound from a coarse wealth grid), so grids with small bets may carry a provably sharper short-horizon bound, but we have no proof."
  * `E11.md:56`: replace "for the default mixture the best adversary found reaches 42% of the bound and a certified sharper bound is open." with "for the default mixture the best i.i.d. adversary found reaches 42% of the bound for the increase grid, an adaptive adversary reaches 52% for the two-sided detector (E16), and a certified sharper bound is open."

---

## 4. Minor findings

### m1. The range of B ignores the monitor's fallback for zero returns

* **Where.** `E16.md:7`: "With τ = 3 times the scaled MAD of the law itself, B can only lie between 1.90 and 39.6"; `e16_adversarial_null.py:793` to `794` ("B lies in [1.90, 39.6] (report, Q16b)"; the derivation is in the set-up, not in Q16b).
* **Evidence.** The derivation is right when the median absolute return is positive. With more than half of the returns at zero, the protocol rule that E16 itself uses (`tick_tau`, the median of the nonzero absolute returns) gives larger B. E16's own price-tick part has B = 47.6, 74.3, 99.3 and 148.6 at ticks of 2, 3, 4 and 6 sd (`e16b_price_ticks.csv`). The i.i.d. twin has B = 60.8, 435 and 7,328 at 2, 4 and 6 sd. The inside-class search stops at B = 39. Bold play beyond it does not change the conclusions: with 20,000 fresh paths each (`highB.py`), B = 100, 150 and 1,000 give 1.98%, 2.08% and 2.21% in one day and 15.1%, 14.7% and 11.5% in five, below the values at B = 2.
* **Fix.** `E16.md:7`: replace "With τ = 3 times the scaled MAD of the law itself, B can only lie between 1.90 and 39.6" with "With τ = 3 times the scaled MAD of the law itself and a positive median absolute return, B can only lie between 1.90 and 39.6". Add after "the Gaussian value is 9.045": "When more than half of the returns are zero, the monitor's fallback (the median of the nonzero absolute returns) allows larger B, up to 149 for the price ticks of Q16f." Consider adding B = 100 and 150 to part a4. `e16_adversarial_null.py:794`: replace "(report, Q16b)" with "(report, set-up; when the median absolute return is positive)".

### m2. Shape switching ran on 102,560 and 41,020 paths, not 100,000 and 40,000

* **Where.** `E16.md:9`: "(100,000 and 40,000 for shape switching)"; `e16_adversarial_null.py:609` to `615`.
* **Evidence.** The loop `for s in range(0, Pf, ch)` draws full chunks of `ch = 10,256` (one day) and `2,051` (five days) paths, so it runs 10 × 10,256 = 102,560 and 20 × 2,051 = 41,020 paths. `e16a_dependent.csv` says so (`paths`).
* **Fix.** `E16.md:9`: replace "(100,000 and 40,000 for shape switching)" with "(102,560 and 41,020 for shape switching)".

### m3. "The i.i.d. law {0, B}" is the searched two-point law

* **Where.** `E16.md:34`: "(bold play 2.40% and 17.6% against 1.92% and 14.5% for the i.i.d. law {0, B} at the same B)".
* **Evidence.** 1.92% and 14.5% are the best two-point laws of the search, {0.00005, 8.989} and {0.102, 9.041} (`e16a_iid.csv`). The law {0, B} itself gives 1.88% and 14.9% (`e16a_curves.csv`, 30,000 paths). The gain at five days is then 2.7 points, not three.
* **Fix.** Replace "against 1.92% and 14.5% for the i.i.d. law {0, B} at the same B" with "against 1.92% and 14.5% for the best i.i.d. two-point laws at the same B (the law {0, B} itself: 1.88% and 14.9%)", and "three points at five days" with "about three points at five days".

### m4. The middle deciles of the baseline error

* **Where.** `E16.md:156`: "but 0.8% to 1.3% in the four middle deciles, close to the oracle's 0.9%."
* **Evidence.** With a 390-step calibration and q = 1 (`e16b_plugin.csv`), deciles 4 to 7 give 1.70%, 1.13%, 0.83% and 0.97%. The range 0.8% to 1.3% is deciles 5 to 8 (1.13%, 0.83%, 0.97%, 1.27%).
* **Fix.** Replace "but 0.8% to 1.3% in the four middle deciles" with "but 0.8% to 1.7% in the four middle deciles (4 to 7)".

### m5. "Never above 0.59" and "at or below the bound in every case" are for the fresh runs only

* **Where.** `E16.md:14` and `:156` ("never exceeds 0.59 of the bound", "never above 0.59 of the bound"); `E16.md:156` ("a 1440-step calibration stays at or below the bound in every case found").
* **Evidence.** In the search runs (6,000 paths per law, `e16b_plugin_search.csv`) the τ-hat oracle reaches 0.63 of the bound (calibration 390, five days). Three laws with a 1440-step calibration have point estimates of 25.08% to 25.27% at five days. None of this is significant (standard error about 0.56 points), so the conclusions stand, but the statements are about the fresh runs.
* **Fix.** `E16.md:156`: replace "(never above 0.59 of the bound)" with "(never above 0.59 of the bound on the fresh paths; 0.63 at most in the 6,000-path search)". Replace "stays at or below the bound in every case found" with "is not significantly above the bound for any law (fresh paths: at most 24.5%; search: at most 25.3% on 6,000 paths)".

### m6. "Close to optimal" and "nearly sharp bet by bet" go beyond lower bounds

* **Where.** `E16.md:84`: "For a single bet the same policy is close to optimal:"; `E16.md:180`: "For single bets the same adversary comes within a few hundredths of the bound fraction reached by E11's dynamic programme, so the bound is nearly sharp bet by bet."
* **Evidence.** Both bold play and the dynamic programme give lower bounds. Matching a lower bound does not show optimality, and the optimum is unknown. For the bet 0.05 at five days a finer dynamic programme (0.773, M3) is at least as high as bold play (0.766). At 98 steps bold play is 0.08 below the dynamic programme for the decrease bet (0.38 against 0.46). For the bet 0.05 the best adversary known reaches only 0.26 of the bound in a day, so "nearly sharp bet by bet" does not hold for that bet within a day.
* **Fix.** `E16.md:84`: replace "For a single bet the same policy is close to optimal:" with "For a single bet the same policy is about as strong as the best adversary known:". `E16.md:180`: replace the sentence with "For the single bets 0.5 and 0.2 and the decrease bet the same adversary comes within a few hundredths of the bound fraction reached by E11's dynamic programme, so the bound is nearly sharp for these bets at one day or more; for the bet 0.05 the best adversary known reaches only 0.26 of the bound in a day."

### m7. The decomposition: Monte Carlo noise, and what "overshoot" contains

* **Where.** `E16.md:67` ("P(alarm within m) = (m / A) × F1 × F2 × F3"), the table at `:71` to `:73`, and `:84` ("alarms land on 2A instead of about 1.27 times 2A").
* **Evidence.**
  * In the sample the two sides of the identity differ by -1.0% (bold play: 382.6 against 386.3), +1.2% (Gaussian) and +0.6% (two-point), so F1 F2 F3 = 0.487 next to the measured 0.477 for bold play (`e16a_slack.csv`, `identity_lhs`, `identity_rhs`). `theory_delay.md` 8.2 says that the identity holds up to Monte Carlo noise; `E16.md` does not.
  * F3 = A / E[N_σ] with `N = (M↑ + M↓) / 2`, so the 1.27 includes the statistic of the side that did not alarm, not only the overshoot of the alarming side.
* **Fix.** After the table at `:73` add: "The identity holds up to Monte Carlo noise: its two sides differ by 0.6% to 1.2%, which is why F1 F2 F3 and the measured ratio differ in the second digit." At `:84` replace "alarms land on 2A instead of about 1.27 times 2A" with "at an alarm the two sides together hold 2A instead of about 1.27 times 2A".

### m8. How bold play and the factor of two are described in the README and the pitch

* **Where.** `README.md:24` ("which watches the detector's state and bets everything on landing exactly at the threshold") and `PITCH_BAPAT.md:43` ("which watches the detector and bets everything on landing exactly at the threshold"); `README.md:146`, `PITCH_BAPAT.md:52` and `docs/limitations.md:12` ("about two times").
* **Evidence.** The adversary chooses data, not bets. Bold play draws from the widest legal spread {0, B} at every step and lands exactly on 2A only when one draw can reach it (`e16_adversarial_null.py:258` to `287`). The factor of about two is for one day (5% / 2.62% = 1.9); over five days it is 1.4 (25% / 18.1%).
* **Fix.** `README.md:24`: replace "which watches the detector's state and bets everything on landing exactly at the threshold" with "which watches the detector's state, always uses the widest legal spread and, when one draw can reach the threshold, lands exactly on it". `PITCH_BAPAT.md:43`: the same replacement for "which watches the detector and bets everything on landing exactly at the threshold". `README.md:146`, `PITCH_BAPAT.md:52`: replace "about two times" with "about two times at one day (1.4 times over five days)". `limitations.md:12`: replace "about 2 times larger" with "about 2 times larger at one day (1.4 times over five days)".

### m9. The GARCH and tick thresholds are for one horizon

* **Where.** `PITCH_BAPAT.md:43` ("clustering from a mean drift of 3.9%"), `_header.md:30` ("GARCH from a mean conditional-mean drift of 3.9%, ... price ticks from 4 sd"), `E16.md:14` ("from a mean drift of 3.9%").
* **Evidence.** 3.9% is the smallest drift with a failure at one day. Cells with drifts of 4.8%, 4.9% and 5.6% stay below the bound at one day. At five days failures start from a drift of 1.6% (`e16b_garch.csv`). Price ticks fail from 4 sd at five days and from 6 sd at one day (`e16b_price_ticks.csv`).
* **Fix.** In the three places, replace "from a mean drift of 3.9%" (or "from a mean conditional-mean drift of 3.9%") with "from a mean drift of 3.9% at one day (1.6% over five days)". In `_header.md:30` replace "price ticks from 4 sd" with "price ticks from 4 sd (five days; 6 sd for one day)".

### m10. The landing test does not test landing

* **Where.** `tests/test_e16_adversarial.py:68` to `78`.
* **Evidence.** The docstring says "no alarm overshoots by more than one unit", but the test only checks `n_at_alarm ≥ thr / 2` and that the median of `M_alarm + M_other` is below `1.2 thr`. A policy without landing would pass whenever the other side is small. Bold play does land by construction (`e16_adversarial_null.py:277` to `284`), and F3 = 0.9986 in `e16a_slack.csv` agrees; the test just does not check it.
* **Fix.** In `engine(..., decompose=True)` also return the alarming side's statistic, and assert `thr ≤ M_alarm ≤ thr + 1` for every alarm of `run_bold(E.BOLD, ...)`.

### m11. Seed overlaps not covered by the note

* **Where.** `e16_adversarial_null.py:62` to `64`; `docs/iterations/E16.md:43`; caveat `E16.md:187`.
* **Evidence.** Part b3 uses seeds 16400 to 16410 (`:1043`), the same seeds as the b2 search with a 1440-step calibration for laws 0 to 10 (`:959`). Seed 16911 is used by a1 (ticks, five days, `:558`) and by a2 (shape switching, one day, `:610`). Seeds 16921 and 16925 are shared by the three adaptive families. No reported rate is selected on its own seeds, so nothing is biased. But the note "a few fresh seeds are shared by different families" does not mention that a reported part (b3) reuses search seeds.
* **Fix.** Move b3 to unused seeds (for example `SEED + 500 + i`) and rerun it (132 s), or add to the docstring: "b3 reuses the seeds of the b2 search with a 1440-step calibration; the two parts are unrelated, so this only makes them dependent."

### m12. A related misstatement in the iteration log (outside the scope)

* **Where.** `ITERATION_LOG.md:369`: "About half of the remaining slack is expected wealth stranded on paths that never alarm."
* **Evidence.** F2 = 0.49 means that about half of the expected *wealth* is stranded, and that accounts for almost all of the remaining slack (F1 = 0.99, F3 = 1.00).
* **Fix.** Replace with "Under bold play about half of the expected wealth is stranded on paths that never alarm, and that accounts for almost all of the remaining slack."

---

## 5. What I checked and found correct

**The null class (question 1).**

* **Bold play and the 8-parameter family** (`e16_adversarial_null.py:248` to `297`). At step t the two points `(y_low, y_high)` are functions of `(Su, Sd)` after step t - 1, of t and of fixed parameters, so they are known before the draw. They satisfy `0 ≤ y_low ≤ 1 ≤ y_high ≤ B`: the far points are `fl` in [0, 1] and `1 + (B - 1) fh` with `fh` in [0, 1]; the landing points are used only when `ysu ≤ B` (then `max(ysu, 1)`) or `ysd ≥ 0` (then `min(ysd, 1)`). `two_point_draw` puts weight `(1 - y_low) / (y_high - y_low)` on the high point with a fresh uniform, so `E[Y_t | everything before t, including all earlier uniforms] = 1` exactly. The degenerate cases (`y_high = 1` or `y_low = 1`) give Y = 1. So `X_t = m0 Y_t` lies in [0, τ²] (a return of 0 or one at least τ in size, or a landing value), and the conditional mean is m0 for both sides.
* **The landing.** Each side's mean statistic is affine in y: `mean((1 + λ(y - 1))(S + 1)) = a + (y - 1) b` with `a = mean(S) + 1`, `b = mean(λ (S + 1))`, so `ysu = 1 + (2A(1 + 10⁻⁹) - a) / b` takes the increase side exactly to the target, and likewise for the decrease side. If both landings are available, both points alarm. If `a` already exceeds the target, the policy plays Y = 1 and alarms with an overshoot below one unit. The path is removed at its first alarm (`engine`, `:227` to `:236`), and a false alarm is the first alarm, so nothing after the landing step enters any count.
* **Other families.** Two-point and three-point laws (the three-point weights give mean `y1 + p2 (y2 - y1) + p3 (y3 - y1) = 1` exactly); hold-then-bold (Y = 1, then bold play); the hidden Markov chain is independent of the draws and both shapes have mean 1 with the same B; the return laws use their own τ and m0: the scale-mixture median by root finding and `m0 = (1 - w) h(τ) + w k² h(τ / k)`; the Student t truncated mean in closed form, `E[T²; |T| < s] = ν / (ν - 2) I(s² / (ν + s²); 3/2, ν/2 - 1)`, which I re-derived from `T² / (ν + T²) ~ Beta(1/2, ν/2)`; ticks with the exact discrete law. `h(k) = 2Φ(k) - 1 - 2kφ(k) + 2k²(1 - Φ(k))` is right, and `TAU_G = 3`, `B_GAUSS = 9.045`.
* **Two points suffice** at each step for any value function, because the next state is affine in Y_t (the concave envelope argument of Lemma 2a in `theory_delay.md`).
* **The range 1.90 to 39.6** follows from `P(|r| ≥ med) ≥ 1/2` and `P(|r| ≤ med) ≥ 1/2` with `τ² / med² = (3 × 1.4826)² = 19.78`, when the median is positive (m1).
* **Floating point.** The search draws float32 uniforms (`crn`, `:439`), which moves the conditional mean by at most about `B × 2⁻²⁴`, 5e-7, in the search only. Every reported rate uses float64 uniforms.

**The detector and the counting (question 2).**

* `UPG` and `DNF` are imported from `detectors.py`. The recursions `S↑ = (1 + λ(Y - 1))(S↑ + 1)` and `S↓ = (1 + φ / (B - 1) (1 - Y))(S↓ + 1)`, the alarm when a side's mean reaches 2A, and the order (increase side first) are those of `VolatilityMonitor.update` (`monitor.py:329` to `351`, default `c_down = c`). The tests confirm alarm-for-alarm agreement with the monitor for the true baseline (three laws, including bold play) and the plug-in baseline.
* `rate_row` counts alarms at indices 0 to m - 1, so "within m steps". The Wilson interval `(p + z²/2n) / (1 + z²/n) ± z √(p(1 - p)/n + z²/4n²) / (1 + z²/n)` is the standard one (I recomputed 2.495% to 2.751% for 1,572 of 60,000). "Above the bound" means the lower end is above m / A, and for the one-sided rows the bound is changed to m / (2A).
* Seeds: no reported rate in part a or b2 uses a seed its search used. The search seeds are 16101 to 16135 and 16300 to 16436; the fresh seeds are 16900 to 16988 and 16991 to 17016 (overlaps in m11). The 8-parameter policy and the hold length were chosen on the training paths only.

**The decomposition (question 3).**

* `decomposition` computes `F1 = E[σ ∧ m] / m` with σ the 1-based alarm step, `F2 = 1 - E[N_m; σ > m] / E[σ ∧ m]` from the survivors' end states, and `F3 = A / E[N_σ | σ ≤ m]` with `N = (M↑ + M↓) / 2`. These are the factors of `theory_delay.md` 8.1. Under an exact null `N_t - t` is a martingale, so optional stopping at the bounded time `σ ∧ m` gives the identity with equality. All three factors are at most 1 (F3 because `N_σ ≥ A`), which is Theorem B.
* The stranded shares by bet use the weights 1/12 (increase) and 1/10 (decrease). 60% and 36% (bold play), 61% and 31% (Gaussian), 60% and 36% (two-point) agree with `e16a_slack.csv`. Bold play's gain over the i.i.d. law {0, B} is all in F3 (0.79 to 1.00), as the reading says.

**The comparison with E11 (question 4).**

* The single-bet runs use the same B, A, recursion and side as E11's programme, and the table in `E16.md:77` to `:82` agrees with both CSVs. Bold play beating the programme for the bet 0.05 is consistent, and E16's statement that the programme's small-bet values are loose lower bounds is right (M3 gives the reason and the refined values).

**Outside the class (question 5).**

* GARCH: `s2` at index t is the predictable variance of the return at index t, so `μ_t = s2 h(τ / √s2)` is the exact conditional mean for Gaussian innovations. The oracle bets against μ_t with decrease bets `φ / (τ² / μ_t - 1)`; these are predictable and below their cap, so the oracle is inside the class on the same paths. The drift statistics, 43 of 64 and 51 of 64 cells, Spearman 0.989 (one day) and 0.970 (five days), the cuts 3.9% and 5.6% (one day) and 1.6% and 2.5% (five days), the increase-side share 33% to 86%, and the oracle maxima 0.24 (1.2%) and 0.56 (14.1%) all agree with `e16b_garch.csv`.
* Plug-in: `plugin_baseline` is the monitor's calibration. The τ-hat oracle bets against `m(τ-hat)`, the true mean at the estimated level, on the same paths, so it is inside the class. The shape adversary keeps that mean, and at q = 0 the oracle's statistic is `M_t = t` and cannot reach 2A = 15,600 within 1,950 steps.
* Price ticks: starting the efficient price at a uniform position within a tick keeps that position uniform, so the returns are stationary and the pooled τ and m0 are the true unconditional values. The conditional means and the correlations agree with the CSV.

**Numbers (question 6).** Every cell of the tables in Q16a to Q16f agrees with its CSV after rounding, including every bold mark and every interval. The numbers in the text agree too, apart from m2 to m5:

* 0.48 and 0.71, ±0.07 points, 2.14% at B = 2.05, 43%, 2.82% against 2.63% in sample, 0.45 points, 2.39% and 17.8% for hold-then-bold;
* the horizon curve and its ratios, 2.0% to 2.6% and 16.5% to 18.1% over B, 2.67% (2.55% to 2.80%) and 17.7% in the independent run;
* 0.99, 0.49 and 1.00, 1.27, 0.01 to 0.03, 0.26 against 0.18 and 0.77 against 0.69;
* 40% to 44% (search) and 41% to 42% (fresh), 5.10% (4.86% to 5.36%), 24.5% (24.1% to 25.0%), 16.0% and 7.5% at the extreme deciles, 32.5%, 85.8%, 16.1% and 73.1%;
* 0.98 to 1.03, 0.93, 0.88 and 0.82, 1.18, 1.48 and 2.20, 2.33% against 1.03%, 28.8% (28.0% to 29.7%), 17.6% and 66.7%, 0.24 and 0.56;
* 2.7 and 1.5 times the Gaussian use of the bound, and the factor falling from 2.0 to 1.4 at five days.

The E16 numbers in `README.md:24` and `:146`, `PITCH_BAPAT.md:43` and `:52`, `_header.md:30`, `limitations.md:12` and `theory_delay.md:340` (2.62%, 2.50% to 2.75%, 18.1%, 0.97% and 12.4%, 52% and 73%, 3.9%, 42%, 32.5%, 4 sd, about two times) all agree with the CSVs. The report part and the index row are copied verbatim into `docs/experiment_report.md`. The two figures show the CSV values. The iteration notes agree with the CSVs (20 tests, 41,020 paths, run times not checked).

**Reproducibility (question 7).** `pytest`: 20 passed. Three cells (one inside the class, two outside) rerun exactly. Two cells rerun with fewer paths and new seeds agree within Monte Carlo error. Over three independent runs of bold play at B = 2 (60,000, 60,000 and 20,000 paths) the one-day rate is 2.62%, 2.67% and 2.735%, so the headline value is not inflated by being the largest of the B scan.
