# TH5 iteration notes: limits for every valid procedure, the plug-in baseline in the joint regime, and the decrease bet at its boundary

Owner: TH5. Files: `src/avmonitor/limits.py`, `tests/test_limits.py`, `experiments/th5_limits.py`, `results/th5_*.csv`, `results/figures/th5_*.png`, `docs/theory_limits.md`, `docs/report_parts/TH5.md`, this file. All numerical results are [SYNTHETIC].

---

## Iteration 1. A floor for every valid procedure

**OBJECTIVE.** Section 2 and M3 of the referee report ask whether the decrease speed limit of `docs/theory_delay.md` (Corollary 4: 58 steps at c = 3 for any Shiryaev-Roberts e-detector with linear bets, 73 for the default grid) is a property of e-detectors or of the problem. The aim was a bound that holds for every procedure with the guarantee of Theorem B, `P(T ≤ n) ≤ n / A` over the nonparametric null.

**PROBLEM IDENTIFIED.** Corollary 4 of `docs/theory_delay.md` bounds the factor of each valid bet, so it says nothing about rules that are not built from bets. Lorden-type lower bounds are asymptotic in A and concern the average delay, not a path.

**HYPOTHESIS.** The null class contains independent, non-identical product laws. For a fixed path one can place, at each step, a two-point law with mean 1 that puts the largest possible mass on the observed value. The path then has null probability `1 / W_n(y)`, and a valid rule cannot stop on it before `W_n(y) ≥ A / n`.

**CHANGES IMPLEMENTED.** `docs/theory_limits.md` section 2: Lemma 1 (largest null mass), Theorem 2 (deterministic rules, randomised rules, sharpness through the singleton rule and the chain rule over adaptive laws), Corollary 3 (no truncation: no valid rule alarms before step A), Corollary 4 (collapse and explosion floors), remarks on raw returns, i.i.d.-only nulls, one-horizon validity and the two-sided class. Code: `max_null_mass`, `log_floor_factor`, `floor_time`, `collapse_floor`, `randomised_floor_cdf`. Literature search for the lemma (search listings only; arxiv.org was not reachable).

**EXPERIMENTS PERFORMED.** Part `collapse` (closed form, no simulation), c = 1.5 to 5, levels A and 2A.

**RESULTS.** At c = 3: every deterministic valid rule needs 45 steps on a complete collapse (corrected in iteration 5, which also gives the randomised value 43.4; the first version said every valid rule) at A = 7800 (50 at 2A), against 57.7 (63.6) for every SR e-detector with linear bets and 73 (79) for the default grid; at c = 2: 23 (25), 28.7 (31.3), 36 (39). Explosion floor at c = 3: 4 steps, against 5 and 7. A rule valid only at one horizon with α = 0.05 needs 26 steps on the collapse. `results/th5_collapse.csv`.

**FAILURE.** The first evaluation of the default grid's collapse floor used the factors `1 + s_l` instead of `1 + s_l / (B - 1)` and gave nonsense; corrected before any number was written (the corrected values 36, 73, 122, 180 match Corollary 4 of `docs/theory_delay.md`). A test now pins the value 73.

**IMPROVEMENTS.** A floor that holds for every deterministic valid procedure, path by path, and is exact for each path; for randomised procedures, a bound on the alarm probability.

**REMAINING ISSUES.** A pathwise floor is not a delay bound for a law; it ignores values beyond 1 (which no single null law can rule out).

**NEXT ITERATION.** A delay floor that also uses the information in those values.

---

## Iteration 2. A delay floor and its asymptotics

**OBJECTIVE.** Turn the floor into a lower bound on the delay of every valid procedure for the post-change laws of E4 and E11, and measure the efficiency of the e-detectors against it.

**PROBLEM IDENTIFIED.** For increases and for fat tails the pathwise floor is far below the e-detectors' delays (Student t3 at sd ratio 0.8: about 48 steps against 281), because it uses only the values on one side of 1.

**HYPOTHESIS.** The least favourable null law `P* = Q / (1 + λ* U)` of `docs/theory_delay.md` (Lemma 2b) is in the null class with mean exactly 1. A valid rule's stopping event by n is a test of `P*` at level `n / A`, so the Neyman-Pearson lemma bounds `Q(T ≤ n)`. For deterministic rules the test can be restricted to the set where the floor allows stopping.

**CHANGES IMPLEMENTED.** Theorem 5 (pathwise floor, randomised version, constrained Neyman-Pearson bound), Proposition 6 (`N_A ~ log A / γ_Q`, `γ_Q ≥ KL_inf` with equality exactly when Q has no mass beyond 1, `log A / KL_inf` for every valid rule by change of measure). Code: `np_power`, `np_cdf`, `delay_lower_bound`, `quantile_lower_bound`, `sr_alarm_times`, `least_favourable_bet`.

**EXPERIMENTS PERFORMED.** Prototype (4000 paths) in the scratchpad, then part `floors`: 18 laws (Gaussian, Student t3, jumps; sd ratios 0.5, 0.67, 0.8, 1.25, 1.5, 2), 10,000 paths of 4000 steps each, levels A and 2A; part `asymptotic`: six laws at A = 10³ to 10⁸, 4000 paths.

**RESULTS.** Over 18 laws at A = 7800 (`results/th5_floors.csv`): the default grid takes 1.69 to 2.16 times the lower bound for deterministic valid rules and the growth-optimal single bet 1.34 to 1.76 times (1.65 to 2.09 and 1.33 to 1.72 at 2A). Randomised rules gain at most 15% on the bound. No alarm of either detector, on 180,000 paths at two levels, came before the pathwise floor. As A grows from 10³ to 10⁸ the floor approaches `n*(γ_Q)` (within 0.4% to 5% at 10⁸), the Neyman-Pearson bound approaches `n*(KL_inf)` from below (0.85 to 0.98) and the best bet from above (1.16 to 1.19) (`results/th5_asymptotic.csv`). The default grid's delays agree with E11a's change-at-the-first-step delays (93.2 against 93.2 for the Gaussian halving).

**FAILURE.** The first full run took `λ*` from `design.best_bet`, which maximises the growth to a tolerance in λ. For the Gaussian halving the optimum is at `λ_max (1 - 2.4 × 10⁻⁸)`, and the maximiser was at `1 - 3.0 × 10⁻⁸`: the growth agreed to ten digits, but `P* = Q / (1 + λU)` had total mass 0.982 instead of 1, because the atom at B contributes `q_B / ε`. A law of mass 0.982 is not in the null class, so the Neyman-Pearson bound for that cell was not valid. Found by checking the mass column (`p_star_mass`) after the run. Fixed by `limits.least_favourable_bet`, which solves the first-order condition in `log ε`; the mass is now 1 to 10⁻⁸ in every cell (`p_star_mass` in `results/th5_floors.csv`), a test checks it, and parts `floors`, `asymptotic` and `boundary` were run again. The Monte Carlo bounds did not change at the printed precision: the two measures agree on `[0, B)` to a relative 10⁻⁹ and differ only in the mass at `Y = B`, which none of the 4 × 10⁷ simulated values of that cell reached. The fix matters for the proof, not for these numbers. The other 17 cells had mass 1 to six digits before the fix.

**IMPROVEMENTS.** A lower bound on the delay of every valid procedure, with an efficiency ratio for the e-detectors.

**REMAINING ISSUES.** The Neyman-Pearson bounds are Monte Carlo estimates. The floor covers the change at the first step (the worst case in Lorden's sense), not a change after a pre-change history.

**NEXT ITERATION.** The plug-in baseline.

---

## Iteration 3. The plug-in baseline when n, m and A grow together

**OBJECTIVE.** Section 3 and M3 of the referee report: an upper bound of the form `(m / A)(1 + C m / √n)` with C explicit, a matching lower bound or counterexample for `m / √n → ∞`, the margin z for an inflated plug-in, and its relation to the confidence-sequence baseline.

**PROBLEM IDENTIFIED.** TH4's Corollary 1 bounds the false alarm probability by an average over windows, but gives no rate in n. Proposition A.6 and Lemma A.7 of `docs/theory_running_tau.md` work at fixed n through windows of probability about 10⁻⁵⁴⁶.

**HYPOTHESIS.** The drift of the plug-in is controlled by the Kolmogorov distance of the window, `|d| ≤ B̂ D_n` (because the truncated mean is the integral of the survival function up to τ̂²), so DKW gives a distribution-free rate and the central limit theorem a first-order constant. For a lower bound, a stream in the null class can first hold the drift deterministically (building up the Shiryaev-Roberts statistic like `S_h(β|d|)`) and then bet boldly on a two-point law, turning a drift of order `n^{-1/2}` into an alarm with probability of order `(m / A)(1 + m / √n)`.

**CHANGES IMPLEMENTED.** Lemma 7, Theorem 8 (a) explicit, (b) distribution-free, (c) first order with `C_F = (λ̄ + φ̄) σ_Y / (4 √(2π))`; Proposition 9 (hold, then bold) and Corollary 10 (dichotomy at `m ≍ √n`); Proposition 11 (inflated plug-in: `z = Φ⁻¹(1 - δ/2)` is sufficient for every stream and necessary for the worst stream when `m / (√n log A) → ∞`); the relation to the cs baseline; the connection with A.6 and A.7. Code: `plugin_windows`, `kolmogorov_distance`, `side_bound`, `first_order_constant`, `distribution_free_constants`, `dkw_log_moment`, `distribution_free_bound`, `noise_free_alarm`, `hold_bold_lower_bound`.

**EXPERIMENTS PERFORMED.** Part `plugin_bound` (windows only; n = 390 to 92,160; 20,000 windows up to n = 23,040 and 5000 beyond; m = 98, 390, 1950, 7800; Gaussian and t3), part `plugin_sim` (monitoring: 10,000 paths per cell, n = 390 to 23,040, 7800 monitored steps, six baselines on the same monitored returns), part `holdbold` (the stream of Proposition 9 run through the detector, 4000 windows per cell).

**RESULTS.** [SYNTHETIC] A = 7800. (1) `√n sd(d)` is within 1.3% of `σ_Y` at every n from 390 to 92,160 for both laws; `max |d| / (B̂ D_n)` is 0.55 (Gaussian) and 0.76 (t3); the share of windows with a positive growth rate under F falls from 0.83 to 0.005 (Gaussian) and 0.76 to 0 (t3) (`results/th5_plugin_drift.csv`). (2) Where `m / √n < 1`, Corollary 1 is at most 4.1% (Gaussian) and 5.3% (t3) above nominal and the first-order formula is within 1% of it; where `m / √n > 5`, the proved lower bound is above nominal in 19 of 24 cells, up to 6.5 (Gaussian) and 8.9 (t3) times at n = m = 390, and 25% to 96% of the upper bound (`results/th5_plugin_bound.csv`). (3) The hold-then-bold stream, simulated, lies between the proved lower bound and Corollary 1 in all four cells (n = m = 390: 0.338, 0.365, 0.424; `results/th5_holdbold.csv`). (4) Distribution-free bound: trivial for n ≤ 92,160 and m ≥ 98; 1.153 at n = 10⁶, m = 39; 1.073 at n = 10⁸, m = 390 (`results/th5_distribution_free.csv`; under-integrated, corrected to 1.149 in iteration 5). (5) i.i.d. monitored returns (`results/th5_plugin_sim.csv`): the plug-in is at or below nominal at one day in every cell (largest 5.14%, t3, n = 390, interval 4.72% to 5.59%), but at five days after a one-day window it alarms in 40.8% (Gaussian) and 42.3% (t3) of paths against 25%; its ratio to the true baseline at five days falls from 3.4 to 1.07 (Gaussian) and 2.9 to 1.06 (t3) as n goes from 390 to 23,040; the true baseline itself is at 0.85% to 1.64% at one day; z = 2.576 and the cs baseline give at most 0.75% and 0.63% at one day; no rate is above its Corollary 1 bound (288 rows). (6) Inflated plug-in, z = 2.576: non-coverage share 0.009 to 0.013 against the limit 0.01; bound at most 1.061 times nominal. The cs baseline's effective z is 3.95 to 5.74 (increase) and 3.78 to 3.99 (decrease) for Gaussian returns. Reading the i.i.d. excess as consistent with the `log A / √n` heuristic was done after seeing the results [post hoc].

**FAILURE.** In the smoke run the explicit bound of Theorem 8a averaged `e^{m x*}` without the cap and returned values up to 10³⁰³, and the distribution-free bound returned NaN because `e^{k2 m x}` overflowed inside the integral. Fixed by capping per window at `A / m` (which the theorem allows) and by integrating in log space (`dkw_log_moment`). The distribution-free bound is honest but weak: it is informative only for `m` well below `√n / 20` (`results/th5_distribution_free.csv`).

**IMPROVEMENTS.** The question "when does the plug-in keep its level" now has a proved answer for the worst stream in the null class, and the margin needed by the inflated plug-in is identified.

**REMAINING ISSUES.** For i.i.d. monitored returns the boundary is not proved; a heuristic based on the stationary tail of the Shiryaev-Roberts recursion suggests the scale `log A / √n` (section 4.6 of the theory note). The delay cost of the margin z, and of the cs baseline, after a change was not measured.

**NEXT ITERATION.** The decrease bet near `λ_max`.

---

## Iteration 4. The decrease bet at its boundary

**OBJECTIVE.** The referee reports `λ* / λ_max` between 0.9917 and `1 - 2.4 × 10⁻⁸` for the Gaussian decrease cases. Explain this, and give the delay of the bet `λ_max (1 - ε)`.

**PROBLEM IDENTIFIED.** At `λ_max` the factor is 0 when the return reaches the truncation level, so the growth is `-∞` whenever that event has positive probability, but the optimum can be extremely close to it.

**HYPOTHESIS.** The growth splits into a bulk term `G_c` and an atom term `q_B log ε`, so `ε*` is proportional to `q_B`, which is tiny for light tails and not for fat tails.

**CHANGES IMPLEMENTED.** Proposition 12: two-sided bounds on the growth, the condition `E_Q[(B - 1) / (B - Y)] > 1` for an interior optimum, the identity `ε* = q_B / h(ε*)`, the choice `ε₀ = q_B / G_c`, the Wald delay bound and an ε-free block bound for the Shiryaev-Roberts statistic. Code: `boundary_growth`, `boundary_terms`, `boundary_growth_bounds`, `block_delay_bound`.

**EXPERIMENTS PERFORMED.** Part `boundary`: nine decrease laws, ε from 0.5 to 10⁻⁶ plus `ε*` and `q_B / G_c`, 4000 paths of 6000 steps.

**RESULTS.** `results/th5_boundary.csv`. Gaussian: `q_B` from 2.0 × 10⁻⁹ to 1.8 × 10⁻⁴, `ε*` from 2.4 × 10⁻⁸ to 8.3 × 10⁻³ (the referee's range), and every ε of 10⁻³ or less gives a delay within 0.4% of the delay at `ε*`. Student t3: `q_B` from 0.65% to 2.4%, `λ* / λ_max` = 0.93, 0.78, 0.55; at `ε = 10⁻⁶` the growth is negative at sd ratios 0.67 and 0.8 and the delay is 14% and 53% longer than at `ε*`. Jumps: within 1% even where the growth is negative. The Wald bound is infinite where the growth is negative, but the simulated delays stay finite, which the block bound explains.

**IMPROVEMENTS.** The referee's observation is explained and shown to be specific to light tails.

**REMAINING ISSUES.** The block bound is loose (an order of magnitude for t3 at sd ratio 0.8).

---

## Iteration 5. Independent review, round 7

**OBJECTIVE.** Apply `docs/review/independent_review_round7_th5.md` (no Critical finding, five Major, 14 Minor) to the TH5 files.

**PROBLEM IDENTIFIED.** M1: Corollary 4 and its summaries claimed the 45-step collapse floor for every valid rule; Theorem 2b shows it holds for deterministic rules only. M2: `dkw_log_moment` integrated over all of x in [0, 1] with `quad`, and for n ≥ 10⁶ it under-resolved a tail of width a few `n^{-1/2}`, so 14 cells of `th5_distribution_free.csv` were too small (the bound looked tighter than it is). M3: "proved to be at least 34% and at most 42%" was a Monte Carlo average over 4000 windows of per-window proved bounds. M4: the summaries of Corollary 10 dropped `m / log A → ∞` and `m / A → 0`, said "whatever A", called `C(W)` the null class and called the factor `1 + C m / √n` sharp. M5: the Reading paragraph of the report part treated the delay lower bound as the best achievable delay. Minor: 0.73 for 0.72; "at or below nominal" for a point estimate of 5.14%; the comparison at 2A; the limit of `V_n`; `|d| ≤ 1/2` missing from `E_n`; the loose `C_DF`; a test that still used `design.best_bet`, a test threshold that sat on the corrected value; the wording of `ε*`; 150,000 against 10⁶; the five cells below nominal; no Monte Carlo error for the delay bounds; the role of `C(W)` in the cancellation of τ̂; "no violation" called a check of the theorem; the role of the citations in Proposition 6c.

**HYPOTHESIS.** None of the findings changes a proof. The randomised floor on the collapse is a linear programme: with `q_j = P_U(T = j | path)`, validity on the path is exactly `Σ_{j ≤ k} q_j F^{-j} ≤ k / A` for every k. The distribution-free moment is accurate if integrated in `u = √n x`, where the mass sits on a fixed scale.

**CHANGES IMPLEMENTED.**
* `limits.randomised_path_delay` (linear programme, HiGHS) and Corollary 4b and 4c in the theory note; new columns `rand_min_mean_delay` and `rand_max_prob_by_floor_minus_1` in `th5_collapse.csv`; "deterministic" added wherever the 45-step floor is quoted.
* `limits.dkw_log_moment` rewritten: exact integral `u₀ e^{s u₀}` on `[0, u₀]`, then `quad` in u on `[u₀, min(√n, max(u₀, s/4) + 30)]` with breakpoints around the maximum. `limits.distribution_free_rate_constant` gives `C_DF = k1 I₀ / 4 = 1.79` with `I₀ = u₀ + √(2π) Φ̄(2 u₀) = 0.888` (was 2.52 from the bound `I₀ ≤ √(π/2)`). New part `dfree` of the script writes `th5_distribution_free.csv` (also run by `plugin_bound`).
* Part `floors` now also writes the bound for randomised rules and its ratios (`efficiency_*_rand_*`), the two-sided monitor's sides at 2A against the bound at A (`efficiency_*_2A_vs_A`), and batch-means standard errors (`se_*`, 10 batches of 1000 of the same paths; no new random numbers, so every earlier column is unchanged, which was checked).
* Tests: `test_randomised_rules_on_the_collapse`; `test_dkw_moment_against_closed_form` (seven cells up to n = 10⁸, against the closed form of the u-integral, a Gaussian integral); `test_distribution_free_table_values`; the rate test now checks that the excess divided by `(m + 1) / √n` decreases to `C_DF` (at n = 10¹²) instead of a threshold; the Neyman-Pearson test uses `least_favourable_bet`.
* Text: every item of the review in `docs/theory_limits.md` and `docs/report_parts/TH5.md`.

**EXPERIMENTS PERFORMED.** Parts `collapse` and `dfree` (no randomness, under a second each). Part `floors` again (237 s, same seeds). A scratch check of the corrected moment against the closed form and a 4-million-point log-space trapezoid in all 48 cells (agreement to 10⁻⁹ in the log). Two further 20,000-window runs of the hold-then-bold lower bound and Corollary 1 at n = m = 390 (Gaussian, new seeds, scratch only).

**RESULTS.**
* Randomised rules on the collapse at c = 3, A = 7800: smallest mean delay 43.39 (48.46 at 2A), 74% of the mass at step 42; largest alarm probability by step 44 is 0.978. Explosion: 3.72 (3.91), probability 0.28 at step 3. This reproduces the review's 43.40 and 0.978.
* Distribution-free bound, ratio to nominal, old and new: n = 10⁸: 1.013 to 1.021 (m = 98), 1.073 to 1.149 (m = 390), 10.21 to 29.99 (m = 1950); n = 10⁷: 1.053 to 1.102 (m = 98), 1.670 to 3.941 (m = 390); n = 10⁶: 1.0075 to 1.0118 (m = 5), 1.0146 to 1.0238 (m = 10). Unchanged for n ≤ 92,160 and for n = 10⁶, m = 39 and 98; the `distribution_free_ratio` column of `th5_plugin_bound.csv` (n ≤ 92,160) is unchanged, so `plugin_bound` was not rerun. This reproduces the review's table.
* Hold-then-bold lower bound and Corollary 1 at n = m = 390, 20,000 windows: 33.0% and 33.2%, and 42.1% and 42.2%, in two new runs (standard errors 0.33 and 0.30 points), against 32.3% and 41.6% in the CSV and 32.4% to 32.9% and 41.7% to 42.0% in the review.
* Delay bounds: standard error at most 0.7% of the bound, at most 0.012 for a ratio. Against the randomised bound: 1.81 to 2.24 (grid) and 1.53 to 1.80 (best bet). Two-sided monitor's sides at 2A against the bound at A: 1.85 to 2.40 and 1.49 to 2.03.

**FAILURE.** The review's suggested sentence "at least half of the delay is the problem's" does not hold for the default grid at Gaussian sd ratio 0.8 (ratio 2.16, so only 46%) or for jumps at sd ratio 0.8 (2.03); the text gives "at least 46% (at least half in 16 of the 18 laws)". The first claim of the two-sided paragraph in the review (that the floor over `H=` is larger) is reversed: over `H=`, `W_n` is larger and the floor earlier, so the two-sided ratios are at least the one-sided ones; the text says so.

**IMPROVEMENTS.** The collapse floors now state their scope and give the randomised value; the distribution-free bound is computed correctly, with a constant 1.79 instead of 2.52; the plug-in summaries carry the hypotheses of Corollary 10; every efficiency ratio is labelled as an upper bound on inefficiency, with a Monte Carlo error.

**REMAINING ISSUES.** The bound for the two-sided guarantee over `H=` (with `π=`) was not computed. No rule is known to attain the delay floor. The block bound of Proposition 12e is loose.

**NEXT ITERATION.** None planned.

---

## Runs and housekeeping

* Smoke runs at `TH5_SCALE=0.02` into the scratchpad, then full runs with `nice -n 10`, at most two processes: run 1 (`collapse floors asymptotic boundary holdbold`), run 2 (`plugin_bound plugin_sim`), run 3 (`boundary` again after adding the block bound), run 4 (`floors asymptotic boundary` again after the fix of iteration 2), then `figures`. Wall times from the script's log lines: collapse 0 s, floors 290 s (run 4: 426 s), asymptotic 11 s (13 s), boundary 23 s (run 3: 81 s, run 4: 73 s), holdbold 8 s, plugin_bound 1452 s, plugin_sim 951 s, figures about 2 s. After the review (iteration 5): run 5 (`collapse dfree`, under a second each; `floors` again, 237 s). No run was killed..
* `git pull -q --rebase --autostash origin <branch>` failed with "Cannot rebase onto multiple branches"; `git fetch` followed by `git rebase` worked.
* A scratch script run from the scratchpad root imported another agent's `numbers.py` (a module name clash with the standard library) and failed with a circular import; scratch scripts were moved to a subdirectory.
* arxiv.org and Wikipedia could not be reached from this machine; references were checked from search listings only, and the theory note says so for each.
