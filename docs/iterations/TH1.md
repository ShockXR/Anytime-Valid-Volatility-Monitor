# Iteration notes: TH1 (detection delay theory, E11)

Theory work on the delay of the betting e-detector: `docs/theory_delay.md`, `experiments/e11_delay_theory.py`, `src/avmonitor/design.py`, `tests/test_delay_theory.py`, `docs/report_parts/E11.md`. One entry per substantial iteration, failed attempts included. Dates are UTC, 2 November 2026. All results **[SYNTHETIC]**.

## Iteration 1. Growth rates and delay predictions (part A)

* **OBJECTIVE.** Derive the delay of the mixture Shiryaev-Roberts e-detector from the growth rate `g(λ) = E log(1 + λU)` and check it against simulation in the E4 scenarios.
* **PROBLEM IDENTIFIED.** The textbook first-order delay `log(A) / g` overstated the simulated delays by up to 112% for small changes at A = 7800, while looking accurate for large changes.
* **HYPOTHESIS.** The Shiryaev-Roberts statistic sums over start times; `log M_n = S_n + log V_n` with a memory term `ζ = E log V_∞` that nonlinear renewal theory turns into `(log A + κ - ζ) / g`, and a Brownian approximation with Dufresne's identity gives `ζ ≈ log(1/g) + γ`.
* **CHANGES IMPLEMENTED.** `design.py`: exact law of Y by Gauss-Legendre quadrature plus the atom at the cap, growth rates, oracle bets, information numbers, Wald bound with Lorden's overshoot bound, Monte Carlo renewal constants, Brownian ζ, mixture simulator. Tests against Monte Carlo and exact identities.
* **EXPERIMENTS PERFORMED.** 20 cases (three models, both sides), 3000 paths with the change at step 1 and 2000 at step 501.
* **RESULTS.** Second order within 3% in all 20 cases; proved bound 1.3 to 2.7 times the truth; E4 reproduced within about 10% at ν = 501. A side finding: the plug-in quadratic (aGRAPA) bet is 1.4 to 4.4 times too small for increases.
* **IMPROVEMENTS.** A first version of `kl_truncated` returned infinity without truncation because the density underflowed in the tails; rewritten on the log scale (log-sum-exp).
* **REMAINING ISSUES.** The regularity conditions of the renewal theorems were not checked for our law of Y.
* **NEXT ITERATION.** Price of robustness.

## Iteration 2. Price of robustness and the truncation constant (part B)

* **OBJECTIVE.** Explain E4's 10% to 70% penalty against the Gaussian likelihood ratio detector and study c.
* **PROBLEM IDENTIFIED.** No decomposition of the penalty existed; truncation was the suspected cause.
* **HYPOTHESIS.** The penalty is `KL / g_grid`, which splits into truncation, nonparametric (KL_inf duality) and grid losses.
* **CHANGES IMPLEMENTED.** Proofs of Lemma 1 (growth ≤ KL), Lemma 2 (KL_inf duality and domination of every nonparametric e-variable by a bet) and the speed limits; numerical decomposition; simulated penalties with common random numbers.
* **EXPERIMENTS PERFORMED.** Information numbers for c in {2, 3, 4, 5, none}; simulated delays for both detectors; jump-burst sensitivity.
* **RESULTS.** The ratio explains E4 (predicted 1.16, 1.37, 1.95; simulated 1.14, 1.35, 1.73). The hypothesis about truncation was wrong: the largest loss is the nonparametric one. For decreases truncation sets a speed limit; c = 1.5 to 2 is best there, c = 3 to 4 for increases.
* **IMPROVEMENTS.** Recommendation that the two sides use different c (needs no new proof).
* **REMAINING ISSUES.** Not tested on real data.
* **NEXT ITERATION.** Bet design.

## Iteration 3. Bet design on formulas (part C), including two failed attempts

* **OBJECTIVE.** Grid and weights minimising the worst-case ratio of delay to the oracle delay over sd ratios 1.2 to 2 (up) and 0.5 to 0.85 (down).
* **PROBLEM IDENTIFIED.** (a) The first objective gave ratios near zero for absurd grids: the Brownian ζ is huge for mismatched or tiny bets, so `log(A/w) - ζ` went negative and hit the floor. (b) Nelder-Mead stopped in local optima (a three-bet result worse than two bets). (c) The decrease optimisation with free weights did not finish and was stopped.
* **HYPOTHESIS.** Cap ζ at `log(1/g) + γ` and use the form `log1p((A/w) e^{-ζ}) / g`; search a lattice of geometric grids before polishing; use equal weights for decreases.
* **CHANGES IMPLEMENTED.** `sr_delay_approx`, `design_ratio`, `optimise_grid` (lattice then polish).
* **EXPERIMENTS PERFORMED.** Design on Gaussian laws; selection by simulation (seed S1, 2000 paths); validation on a fresh seed (4000 paths) for three models; null check.
* **RESULTS.** Two increase bets (0.2565, 0.4893), one decrease bet at `0.99 λ_max`. Validation: 8% to 11% faster than the default for increases on Gaussian data (change at step 1).
* **IMPROVEMENTS.** The objective now matches the finite-A trade-off later formalised in Theorem 3b.
* **REMAINING ISSUES.** The design used Gaussian laws only.
* **NEXT ITERATION.** Robustness of the chosen grids.

## Iteration 4. Variants and validation of the Gaussian-tuned grid (parts C2, C3)

* **OBJECTIVE.** Check simpler variants and validate at the E4 change time.
* **PROBLEM IDENTIFIED.** The single decrease bet at `0.99 λ_max` loses `log(0.01)` at every capped return.
* **HYPOTHESIS.** Two decrease bets (0.9, 0.99) × `λ_max` are safer; rounded increase bets (0.25, 0.5) are as good.
* **CHANGES IMPLEMENTED.** Parts C2 and C3.
* **EXPERIMENTS PERFORMED.** Fresh seeds, 4000 and 3000 paths, three models, ν = 1 and 501.
* **RESULTS.** The single bet was fragile under t5 (ratio 1.51 at sd 0.9). The pair (0.9, 0.99) cut the worst in-range ratio to 1.12. At ν = 501 the gains shrank and under t5 the tuned grid was 12% to 15% slower than the default in two cases: the default's small bets carry a larger head start into the change.
* **IMPROVEMENTS.** Head-start mechanism added to the theory.
* **REMAINING ISSUES.** Fat tails beyond t5 not yet tested.
* **NEXT ITERATION.** False alarm slack.

## Iteration 5. Why the false alarm bound looks loose (parts D, D2)

* **OBJECTIVE.** Explain E9's factor of about 4 and decide whether a sharper provable threshold exists.
* **PROBLEM IDENTIFIED.** Without a decomposition, "conservative" could mean a loose bound or a wasteful detector.
* **HYPOTHESIS.** An exact optional-stopping identity separates stopping, stranded wealth and overshoot; dynamic programming over the nonparametric null gives the worst case for single bets.
* **CHANGES IMPLEMENTED.** `decomposition` (with a per-bet split of the stranded wealth), `worst_case_false_alarm` (dynamic programming on a log grid of wealth states), two-point i.i.d. adversaries, a 40,000-path two-sided split.
* **EXPERIMENTS PERFORMED.** Gaussian null, horizons 98 to 7800, 10,000 and 40,000 paths.
* **RESULTS.** Factor about 5 = overshoot (0.80) times stranded wealth (0.24), the latter held mostly by the decrease side and the smallest increase bets. Worst null reaches 95% of the bound at one day for a single bet of 0.5; the bound is exact at m = A.
* **IMPROVEMENTS.** Failed attempt recorded: the dynamic-programming upper bound (rounding the other way) was useless because rounding errors accumulate, so only the lower bound is reported. The identity holds only up to 0% to 16% Monte Carlo noise because the wealth is heavy-tailed; the factors are reported as approximate.
* **REMAINING ISSUES.** A certified worst case for the default mixture (11-dimensional state) is open.
* **NEXT ITERATION.** Tracker B73.

## Iteration 6. Asymptotic optimality (tracker B73, parts E and E2)

* **OBJECTIVE.** Show whether the delay attains `log A / KL_inf` as `A → ∞`, state the theorem and check it numerically for A = 10^2 to 10^8.
* **PROBLEM IDENTIFIED.** Theorem 2 covered only fixed grids; the rate at which a grid must be refined, and the second-order cost of a finite grid, were missing.
* **HYPOTHESIS.** Worst-case delay equals the delay from zero (monotonicity); a fixed grid has limit ratio `KL_inf / g_grid`; refining with about `√(log A)` bets gives `½ log log A`; a continuous mixture gives at most `log log A` by a concavity argument.
* **CHANGES IMPLEMENTED.** Lemma 5 and Theorem 3 with proofs; part E (one pass computes all thresholds), part E2 (renewal constants per case; E10 cross-check).
* **EXPERIMENTS PERFORMED.** 15 cases (Gaussian, t3, jumps; up 1.25, 1.5, 2; down 0.5, 0.8), 1000 paths, five detectors.
* **RESULTS.** The asymptotic regime is far away (best bet 5% to 24% faster than the bound at A = 10^8), but the second-order formula explains the gap within 0.024. Fixed grids approach their limits. A failure surfaced: the Gaussian-tuned grid has negative growth under t3 at sd 0.8 (22% of paths never alarm) and limit 1.43 at sd 1.25.
* **IMPROVEMENTS.** E10's factor 0.6 explained by the memory term (i.i.d. version, true baseline).
* **REMAINING ISSUES.** At A ≤ 10^8 first-order optimality cannot be confirmed numerically; the mixture's `½ log log A` is heuristic.
* **NEXT ITERATION.** Redesign the bets across tail shapes.

## Iteration 7. Bet design pooled over tail shapes (parts C4, C4V)

* **OBJECTIVE.** A grid that does not fail under fat tails.
* **PROBLEM IDENTIFIED.** The Gaussian-tuned grid was up to 55% slower than the default under t3 data.
* **HYPOTHESIS.** Pooling Gaussian, jump, t5 and t3 laws in the design objective (decrease bets as fractions of each model's `λ_max`) gives a grid that is never much worse than the default.
* **CHANGES IMPLEMENTED.** `part_c4` with a selection rule fixed before validation; `part_c4v` on a fresh seed.
* **EXPERIMENTS PERFORMED.** Design over 36 (up) and 32 (down) laws; validation with 2000 paths, four models, ν = 1 and 501, inside and outside the design range.
* **RESULTS.** Increase bets (0.16, 0.39), decrease bets (0.5, 0.85) × `λ_max`. At least as fast as the default in all 48 in-range cases (0.4% to 10.5% faster); up to 17% slower for sd ratio 1.1 under fat tails. Negative result: the formulas predicted a worst-case improvement from 1.54 to 1.27, the simulation shows 1.25 to 1.24, because the dominant-bet approximation overstates the cost of many-bet grids.
* **IMPROVEMENTS.** Recommendation made conditional on the alarm level and on the importance of small increases.
* **REMAINING ISSUES.** A third small bet (to cover sd 1.1) was not tested, to avoid tuning on the validation data. The tail shape of real deseasonalised returns, which decides between the designs, was not studied.
* **NEXT ITERATION.** For the project: apply the recommended c per side and grid only after a real-data check (files owned by others).
