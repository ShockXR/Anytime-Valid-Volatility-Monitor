# TH4 iteration notes: the false alarm bound when the null holds only approximately (tracker B84)

Owner: TH4. Files: `src/avmonitor/robustness.py`, `tests/test_robustness.py`, `experiments/th4_robustness.py`, `results/th4_*.csv`, `results/figures/th4_*.png`, `docs/theory_robustness.md`, `docs/report_parts/TH4.md`, this file.

---

## Iteration 1. The sup form

**OBJECTIVE.** Bound the false alarm probability of the default two-sided detector when `E[X_t | past] ≤ m₀ (1 + δ)` (increase side) and `≥ m₀ (1 - δ)` (decrease side), with the inflation factor over m / A in closed form.

**PROBLEM IDENTIFIED.** Theorem B uses that each bet's Shiryaev-Roberts statistic minus t is a supermartingale. Under a drift it is not. A first idea was to discount the statistic by `(1 + λ_max δ)^{-t}` and apply Ville's inequality. That puts the largest bet's factor on every bet, and the discount at the end of the horizon is the full `(1 + λ_max δ)^m`, which is far weaker than necessary.

**HYPOTHESIS.** Taking the expectation of each bet's statistic at the stopped time `T ∧ m`, with the drift factor entering the Shiryaev-Roberts recursion, gives `E[R_{T∧m}] ≤ S_m(λδ)` per bet. Averaging over bets then gives a per-bet inflation `G_m(λ_k δ)`, which is exact for a deterministic stream.

**CHANGES IMPLEMENTED.** Lemma 1 (the part of the mass that has already stopped needs `ρ_t ≥ 1`), Theorem 1 (two-sided and one side), Proposition 1 (geometric mean and sinh series bounds on `G_m`), Proposition 2 (sharpness). Code: `sr_growth_sum`, `growth_factor`, `growth_factor_bounds`, `sup_bound`, `sup_inflation`, `design_drift`, `deterministic_alarm_step`.

**EXPERIMENTS PERFORMED.** Closed-form evaluation only (part `sup`), and unit tests of the closed form, of the elementary bounds and of Lemma 1 by Monte Carlo (replaced by an exact test in iteration 6).

**RESULTS.** The one-day bound reaches 10% at a drift of 1.63% (both sides) and 25% at 2.66%. The decrease side tolerates about four times more drift, because its bets are smaller.

**IMPROVEMENTS.** A proved, closed-form, per-bet inflation factor that reduces to Theorem B at δ = 0.

**REMAINING ISSUES.** It needs a bound on the drift at every step. Under GARCH the conditional mean comes arbitrarily close to `B m₀`, so the only valid constant is `B - 1` and the bound is 1.

**NEXT ITERATION.** A pathwise form for a drift that varies.

---

## Iteration 2. Pathwise forms and what cannot be said

**OBJECTIVE.** Replace `(1 + λδ)^k` by the product of `(1 + λ δ_i)` over the window, and decide what can be said in terms of the mean drift.

**PROBLEM IDENTIFIED.** Lemma 1 needs `ρ_t ≥ 1`, and its expectation bound uses the drift as a number. When the drift is random and responds to the data, the expectation of the product of factors is not the product of the expectations.

**HYPOTHESIS.** (1) Dividing each factor by its predictable drift factor gives a valid detector `Ñ`. The wealth is then at most `Ñ` times the largest window product of the drift factors, path by path. A threshold on that product gives a bound for any predictable drift, signed, so that a calm spell resets the product. (2) For a drift that does not respond to the data, conditioning on the drift path gives Theorem 1 with the drift path in place of the constant. (3) No bound of the form `C(mean drift) × m / A` exists.

**CHANGES IMPLEMENTED.** Lemma 2, Theorem 2 (threshold form), Proposition 2b (truncated sup form), Theorem 3 (exogenous drift, conditional on a σ-field G), the two-step counterexample of 3.4 (a), the Jensen remark of 3.4 (b), Proposition 4 (mean drift with an additive term weighted by A) and Proposition 5 (a seven-step burst decided before monitoring). Code: `drift_sr_sums`, `drift_cusum`, `threshold_bound`, `exogenous_bound`, `truncated_sup_bound`, `mean_drift_bound`, `two_step_counterexample`, and the one-pass path engine `run_paths`. The engine computes the alarm, the drift-corrected alarm, the window products, the largest drifts and the summed drifts, and with `check=True` also the largest ratio `N_t / (D*_t Ñ_t)`.

**EXPERIMENTS PERFORMED.** Unit tests: the CUSUM and Shiryaev-Roberts recursions against brute force, the pathwise inequality on random drift streams, the counterexample exactly and by Monte Carlo, and the burst of Proposition 5 (s* = 7 at A = 7800 and the Gaussian B).

**RESULTS.** All pass. The first version of the test for Theorems 2 and 3 used A = 60 and m = 120. There `m / A = 2`, so every bound was trivially 1 and the test checked nothing. It now uses A = 500, m = 150 and a drift of ±8% that flips with probability ½. The measured rate is 9.9%, Theorem 2 gives 79% and Theorem 3 gives 50%, and the test asserts that both are below 1.

**IMPROVEMENTS.** A proved bound for GARCH-type drift (Theorem 2). A sharper bound for exogenous drift (Theorem 3). Proofs that the mean drift alone does not suffice.

**REMAINING ISSUES.** Theorem 2 has to cover a drift that is high exactly when the wealth is high, so it is expected to be loose for noisy data.

**NEXT ITERATION.** The plug-in baseline.

---

## Iteration 3. The plug-in baseline and Proposition A.6

**OBJECTIVE.** Express the plug-in baseline as a drift, combine it with the distribution of the estimate over calibration windows, and connect the result to Proposition A.6 of `docs/theory_running_tau.md` without contradicting it.

**PROBLEM IDENTIFIED.** With a plug-in baseline, the error of m̂ is a drift that stays constant during monitoring but is random across windows. The decrease bets also depend on B̂.

**HYPOTHESIS.** Theorem 1, applied conditionally on the window (part 3), gives `β_m = E_W[min{1, (m / A) I_m(δ(W); B̂(W))}]`. Since `B̂` lies between 1.90 and 39.6, the constant `C_m` in `C_m m / A` is finite for each m but grows exponentially with m, consistent with A.6.

**CHANGES IMPLEMENTED.** Corollary 1, the list of what a distribution-free version would need (an adapted estimate and a tail bound on its error; E8's confidence sequences certify drifts far too large for one day), and the connection to A.6. Along `m_A = ⌈(log A)²⌉` the capped bound tends to `P(δ(W) ≠ 0)` (dominated convergence), which is 1 for continuous laws, so there it is the trivial bound. By Jensen, A.6's event `g(W) > 0` lies inside `{δ(W) ≠ 0}`; this is consistent with A.6, as any proved upper bound must be. Code: `plugin_drift`, `plugin_bound`.

**EXPERIMENTS PERFORMED.** A unit test that the capped bound per window is at least the alarm indicator of the noise-free stream `Y = 1 + δ(W)`.

**RESULTS.** The test passes. The connection is consistent with the current statement of A.6 (`E[T_A | W] ≤ (log(2AK) + log f_max(W)) / g(W)` and its lower limit `P(g(W) > 0)`).

**IMPROVEMENTS.** A proved bound for the plug-in that holds for every stream with the plug-in's conditional mean.

**REMAINING ISSUES.** No distribution-free tail bound for the plug-in error.

**NEXT ITERATION.** The numerical check against E16.

---

## Iteration 4. Numerical check (smoke runs) and an adaptive adversary

**OBJECTIVE.** Check every bound against simulated streams with a known drift and against E16's cells, and measure how tight the bounds are.

**PROBLEM IDENTIFIED.** Smoke runs at a small scale (`TH4_SCALE`, results not reported) showed that Theorem 1 is far above the rate of i.i.d. streams, which the full run confirmed: Gaussian returns with a 5% upward drift of the truncated mean alarm on 2.7% of days against a bound of 1. The deterministic stream attains the bound only at its alarm step, so i.i.d. streams cannot tell whether the bound is tight at shorter horizons.

**HYPOTHESIS.** An adaptive stream with the same constant drift comes much closer. It holds Y at its mean for h steps, so that every bet gains without risk, then plays bold as in E16 (spread {0, B}, landing exactly on 2A when it can).

**CHANGES IMPLEMENTED.** `drift_bold_fn` in the experiment, using E16's engine. The hold h is searched over 15 values on 4,000 paths and the best is rerun on 20,000 fresh paths. Also the GARCH part, which re-simulates E16's 64 cells with the exact drift, and the plug-in part (100,000 fresh windows per cell, E16's deciles, a calibration length scan). The GARCH part reports two heuristics, labelled as such: Theorem 3 applied to GARCH, and Theorem 1 at E16's stored mean drift. One bug was fixed: `sup_bound` did not broadcast a vector of horizons against the bet grid.

**EXPERIMENTS PERFORMED.** Smoke runs of all parts, then the full runs (part constant 580 s, switching 140 s, garch 834 s, plugin 89 s, with at most 2 processes).

**RESULTS.** These are the full runs, reported in section 5 of the theory note. No measured rate is above any proved bound. The adaptive stream reaches 0.33 to 0.63 of Theorem 1 at one day. For GARCH, the proved bounds are below 1 in 22 of 64 cells at one day (6 to 29 times the measured rate), and the sup form is 1 everywhere. Corollary 1 is within a factor of 1.08 to 1.61 of E16's noise-free shape adversary.

**IMPROVEMENTS.** A measured answer to "how tight", with the gap split into the noise of the stream (i.i.d. against adaptive) and the price of covering a drift that responds to the data (Theorem 2 against Theorem 3).

**REMAINING ISSUES.** Price ticks: the drift given the whole past is not available.

**NEXT ITERATION.** Ticks as a labelled heuristic, a sharpness table, figures.

---

## Iteration 5. Ticks, sharpness, figures, write-up

**OBJECTIVE.** Close the remaining items of the brief and write the theory note, the report part and these notes.

**PROBLEM IDENTIFIED.** (1) The first constant-drift figure drew the deterministic stream as straight lines between grid points, which hid where it starts to alarm. Its legend also covered a label. (2) The name "Corollary 5" collided with "Proposition 5". (3) E16 stored the tick drift given the previous return only.

**HYPOTHESIS.** Not applicable (presentation and bookkeeping).

**CHANGES IMPLEMENTED.** The constant-drift figure is now a 2 × 2 grid by horizon and side. It shows the deterministic stream as an exact step function on a fine grid of drifts, adds the adaptive adversary and the nominal m / A, and puts the legend below the panels. The sup part also writes `th4_sharpness.csv`: the drift at which the bound reaches 1, against the smallest drift at which the deterministic stream alarms. The tick part writes `th4_ticks.csv`, with Theorem 1 at E16's one-step drift, labelled a heuristic. The corollary is renamed Corollary 1 (and the CSV column `bound_corollary1`).

**EXPERIMENTS PERFORMED.** Parts `sup`, `ticks` and `figures` rerun. The rerun of `sup` reproduced the earlier CSV files byte for byte.

**RESULTS.** Theorem 1 is attained where it reaches 1: at one day the bound reaches 1 at an upward drift of 3.740% and the deterministic stream alarms from 3.756%. The tick heuristic is informative only for ticks of at most 1 sd (one day) and was above every measured rate.

**IMPROVEMENTS.** The figure shows the sharpness directly. Every number in the theory note comes from a CSV file under `results/`.

**REMAINING ISSUES.** See section 6 of the theory note. Still open: the worst case at horizons shorter than the deterministic alarm step, a distribution-free plug-in bound, a proved tick bound, and a proved GARCH bound that is informative where the nominal bound fails.

**NEXT ITERATION.** None planned. The natural next step would be a family-specific bound for GARCH, using its Markov structure, that is informative at a mean drift of a few percent.

---

## Iteration 6. Independent review, round 6

**OBJECTIVE.** Apply the findings of `docs/review/independent_review_round6_th4.md` (no critical, 4 major, 9 minor).

**PROBLEM IDENTIFIED.**
* M1: Theorem 3 was false as written. Its `D*` term used `D*_m`, which covers only the windows ending at m, where the proof needs `max_{t ≤ m} D*_t`. It also built `D*` from the drift given G, which need not be G-measurable. The code already used the maximum, and every cell reported as proved had a drift that is a function of G.
* M2: section 3.4 (a) showed that a step of the proof of Theorem 3 fails, not its conclusion. Section 3.4 (b) ("the mean drift gives a lower estimate [PROVED]") was false and contradicted by section 5.3.
* M3: "sharp" and "the right worst case" went beyond Proposition 2. The bound is attained only where it reaches about 1.
* M4: the link to A.6 was called a sandwich. Along `m_A` Corollary 1 tends to `P(δ(W) ≠ 0)`, which is 1 for continuous laws.
* Minors:
  * m1: "upward drift only" wording.
  * m2: the null band was not covered.
  * m3: several test checks could not fail.
  * m4: the factor Theorem 2 loses was misstated.
  * m5: the threshold was chosen on the same paths without a quantified margin.
  * m6: the calibration-length wording.
  * m7: "above the nominal 5%" was not significant.
  * m8: two wording slips.
  * m9: no literature position for Lemma 1.

**HYPOTHESIS.** Not applicable (corrections).

**CHANGES IMPLEMENTED.**
* Theorem 3 is restated with G-measurable `u_t, ℓ_t` and `max_{t ≤ m} D*_t`, with the burst followed by a calm spell as the reason for the maximum. `drift_cusum(..., running=False)` returns `D*_m` itself, for that example.
* Section 3.4 (a) gains a stream with a drift that responds to the data, where Theorem 3's formula gives 10.2% against an alarm probability of 11.1% at one day. Section 3.4 (b) is rewritten: the mean drift does not summarise the pathwise bound in either direction (independent, alternating and persistent examples). A one-day violation of the sup form at the mean drift is added to Proposition 5 (8.9% against 10%).
* The sharpness wording now says the bound is attained where it reaches 1, and that tightness at the design drifts is open. The A.6 connection now states the limit (dominated convergence) instead of a sandwich.
* Minor fixes:
  * a remark on the null band and the margin after Theorem 1;
  * the factor `12 x* m / (1 + x*)`;
  * the Dvoretzky-Kiefer-Wolfowitz margin of 1.75 points;
  * the calibration-length and significance wording;
  * `w_K` in Proposition 1;
  * the Robbins and Siegmund (1971) device, with its volume details marked as not checked.
* New experiment part `counterexamples` (exact, no random draws), which writes `results/th4_counterexamples.csv`.
* Tests:
  * Lemma 1 is now checked exactly, by enumerating every path of a two-point law, with and without a stopping time.
  * The Theorem 1 simulation adds the adaptive hold-then-bold stream at A = 7800, whose rate lies strictly between m / A and the bound.
  * A switching case with a small drift makes Propositions 2b and 4 informative.
  * A test for the maximum over t in Theorem 3.
  * A test for the drift that responds to the data.
* Corollary 1's statement and number are unchanged.

**EXPERIMENTS PERFORMED.** The part `counterexamples` and the test suite. No other CSV file changed, and no reported number changed.

**RESULTS.** The reviewer's values were reproduced exactly:
* 0.876 and 80,628 for the burst;
* 0.523, 0.05, 10.2% and 11.1% (5.3% and 11.1% at m = 7) for the drift that responds to the data;
* 36.3%, 37.3% and 5.15% for the alternating drift;
* 1.44%, 8.9% and 10% for Proposition 5 at one day;
* 5.8%, 7.1%, 21.0% and 1 for the sup form at the mean parts of the switching cells with stay 0.5.

**IMPROVEMENTS.** Every statement labelled proved is now proved as written, and the claims of section 3.4 and of the summaries match what is proved.

**REMAINING ISSUES.** As in section 6 of the theory note. In addition, the Robbins and Siegmund reference lacks verified volume details.

**NEXT ITERATION.** None planned.
