# ALT6 iterations: sequential fixed-accuracy estimation of a volatility level

Prototype folder: `alternatives/sequential_estimation/`. Times are UTC, night of 2 November 2026. Every result is labelled [SYNTHETIC] (simulated) or [REAL] (Binance BTCUSDT one-minute bars).

---

ITERATION NUMBER: ALT6-1 (formulation, classical rules, first confidence sequence)

OBJECTIVE: Turn "how many minutes until the volatility level is known to plus or minus d, with confidence 1 - alpha?" into a fixed-accuracy sequential estimation problem; implement the classical stopping rules (Chow and Robbins, its relative version, Stein's two-stage rule, Robbins' minimum risk rule) on raw and on MAD-truncated squared returns, and an anytime-valid confidence sequence stopped at the target width; compare them on the same simulated paths.

PROBLEM IDENTIFIED: Classical rules for the variance are driven by a running estimate of quarticity, which is infinite or very unstable for heavy-tailed returns, and their guarantees are asymptotic (d to 0). A confidence sequence is exact at every d but only for a bounded statistic, and nobody had said how much more data it needs than the classical rules in this setting.

HYPOTHESIS: (1) Truncating at 3 times the pilot's scaled MAD makes the classical rules behave as their asymptotics promise even for t3 returns and jumps. (2) A confidence sequence stopped at width 2d is exactly valid for the truncated variance at every d; its price is a constant factor that can be computed in closed form.

CHANGES IMPLEMENTED: `seqest.py` (laws with exact truncated moments by quadrature, the stopping rules on running sums vectorised over paths, the predictable plug-in empirical Bernstein sequence with one-sided ranges and three bet schedules, a common-random-numbers engine), `tests/test_seqest.py`. Derivation of the limit `E[N]/n* -> 2 log(2/alpha) / z^2` (1.9206 at alpha = 0.05) for tuned constant bets without intersection, and of the first-order prediction with the exact `psi` and the bet caps (`cs_n_predict`, `opt_plan`). Commit ed57477.

EXPERIMENTS PERFORMED: Quick runs at d = 0.2 (normal, t5, t3; absolute and relative accuracy; 1000 paths); unit tests of every formula and brute-force replays of every rule.

RESULTS: [SYNTHETIC] Chow and Robbins on truncated squares had coverage 0.90 (absolute) and 0.95 (relative) with `E[N]/n*` 0.94 and 0.99; on raw squares under t3 the coverage was 0.44 (absolute). The first confidence sequence (tuned bets, no intersection) needed 4.40 (absolute) and 3.61 (relative) times `n*`, far above the limit 1.92. Optimised constant bets brought this to 3.02 and 2.78 against a prediction of 2.52 and 2.43. The unit test of the limit failed at first: the optimiser's bet grid (1e-3 to 1 times the cap) did not contain the small optimal bets and returned 2.22 instead of 1.92.

IMPROVEMENTS: The bet grid of `opt_plan` was centred on the quadratic optimum (1/30 to 30 times `d_Z / var_Z`, capped), and the limit test passed. A running intersection of the bounds (valid because the guarantee holds at all times) cut the relative-accuracy ratio from 2.78 to 2.47 for optimised bets.

REMAINING ISSUES: Even optimised empirical Bernstein bets cost 2.5 times `n*` at d = 0.2. With intersection the simulated ratios started to fall *below* the "limit" 1.92 at small d, which the derivation did not explain.

NEXT ITERATION: A betting confidence sequence (capital processes on a grid of candidate means) and an explanation of the intersected limit.

---

ITERATION NUMBER: ALT6-2 (a betting confidence sequence, and two bugs that made it stall)

OBJECTIVE: Add a betting confidence sequence (`CS-bet`) that uses the shape of the law of the truncated square, not just its variance, and check that it is valid and fast.

PROBLEM IDENTIFIED: The first version worked at d = 0.2 (`E[N]/n*` 1.97 absolute and 1.85 relative, coverage 0.98) but at d = 0.05 one normal path ran to the cap of 600,000 returns (ratio 195 on that path). Its pilot scale was 1.84 times the true sigma, so tau was large and the mean of Z = min(r^2, tau^2) / tau^2 was only 0.033.

HYPOTHESIS: First guess: the grid of candidate means was too coarse near such a small mean, so an adaptive finer grid would fix it. This was wrong: the finer grid did not help. Second guess, confirmed by printing the bets on that path: (1) the bets were chosen from a histogram of Z with fixed bins of width 0.016 on [0, 1], which lumps half the mass of Z into one bin when the mean is 0.033, so the growth-rate-optimal upper bet collapsed to the smallest bet allowed and the upper bound never moved; (2) linear interpolation of the log-capital between grid points (a valid bound by concavity) was too loose near the bound.

CHANGES IMPLEMENTED: (1) Bins relative to the pilot mean (width `zhat_0 / 8` up to `10 zhat_0`, an overflow bin and a bin for Z = 1), with per-bin sums so that the binned mean is exact. (2) A second valid lower bound on the log-capital between grid points: `F(m0) + F'(m0) y - |F''(m0)| y^2 / 2` with the curvature taken at the end of the cell where it is largest (it is monotone in m because one bet is shared by all candidate means); the code takes the larger of this and the chord. (3) New tests: bounds against exact root finding (never inside, at most 0.05% wider with a 48-point grid), time-uniform coverage on Beta data. (4) The engine was rewritten to draw returns only for paths where some rule is still running, because raw Chow and Robbins under t3 needs about `d^-3` returns. Commit e296f7a.

EXPERIMENTS PERFORMED: The stalled path and the d = 0.05 normal and t3 cells (2000 paths); the new tests.

RESULTS: [SYNTHETIC] Before the fixes, the d = 0.05 normal cell had a mean ratio of 1.94 for `CS-bet` (the stalled path alone added about 0.1; 1.63 on t3 paths). After them the stalled path stopped normally, and in the full study (ALT6-4) `CS-bet` needs 1.57 times `n*` at d = 0.05 on normal paths with coverage 0.969. Time-uniform misses on truncated t3 data were 2.2% (bet) and 2.1% (empirical Bernstein) against a bound of 5%.

IMPROVEMENTS: The betting sequence is now the fastest valid procedure at every d tested.

REMAINING ISSUES: Why intersected sequences beat 1.92; whether `CS-bet` at d = 0.05 (1.85 in this check) was converging to anything.

NEXT ITERATION: The d to 0 limit with intersection, then the full synthetic study.

---

ITERATION NUMBER: ALT6-3 (the limit with intersection is a Brownian first-passage functional)

OBJECTIVE: Explain and compute the d to 0 limit of `E[N]/n*` for an intersected confidence sequence.

PROBLEM IDENTIFIED: The derivation of 2 log(2/alpha) / z^2 assumes the bounds at one time only. With the running intersection, the upper bound keeps its best past value and the lower bound its best past value, at different times, so the stopping time is earlier.

HYPOTHESIS: In units of d and of time `u = t / n*`, the bounds converge (Donsker) to `W(u)/(z u) +- (L/(kappa z^2 u) + kappa/2)` with W a Brownian motion, `L = log(2/alpha)` and `kappa` the bet in units of the tuned bet; the intersected sequence stops when the running minimum of the upper curve minus the running maximum of the lower curve reaches 2. The limit is the mean of that first-passage time.

CHANGES IMPLEMENTED: `seqest.brownian_limit` (Brownian paths on a grid of step 0.001 up to u = 4, chunks of 2000 paths), its closed-form special case without intersection (`2L / (kappa (2 - kappa) z^2)`), a test, and the `limit` part of `run_experiment.py` (`results/a_brownian_limit.csv`, 40,000 paths per cell).

EXPERIMENTS PERFORMED: alpha 0.05 and 0.10; kappa 0.8, 1, 1.2, 1.5.

RESULTS: [SYNTHETIC] With intersection and kappa = 1 the limit is 1.548 (standard error 0.001) at alpha = 0.05 and 1.708 at alpha = 0.10, against 1.921 and 2.215 without intersection; the coverage at the stop in the limit is 0.976 and 0.952 (0.993 and 0.986 without intersection). The bet that is optimal without intersection is not optimal with it: kappa = 1.2 gives 1.475 at alpha = 0.05 and 1.603 at alpha = 0.10.

IMPROVEMENTS: The simulated ratios of the intersected sequences at d = 0.05 (about 1.57 for `CS-bet`) are now explained: they sit just above the Brownian limit, which is the right target.

REMAINING ISSUES: The bets of `CS-opt` and `CS-bet` are tuned for the non-intersected rule; a bet about 1.2 times larger should save about 5% in the limit. Not implemented.

NEXT ITERATION: The full synthetic study.

---

ITERATION NUMBER: ALT6-4 (the full synthetic study and a check of every closed form)

OBJECTIVE: Measure coverage at N, `E[N]/n*`, the distribution of N and the regret for every procedure on the same paths, for four return laws, four accuracies, two error levels and both kinds of accuracy; check the closed forms of README 2.9 by simulation.

PROBLEM IDENTIFIED: The quick runs of ALT6-1 to ALT6-3 used one accuracy at a time and few paths; the closed-form predictions had been checked only against themselves.

HYPOTHESIS: (1) Chow and Robbins on truncated squares is close to nominal for every law at every d tested, and on raw squares fails for t3 and jumps at every d. (2) The skewness of squared returns makes the absolute rule under-cover with one-sided misses, and the relative rule does not (first-order bias zero for gamma laws). (3) Without intersection the empirical Bernstein sequences follow the first-order prediction and tend to 1.92; with intersection the betting sequence approaches the Brownian limit 1.55.

CHANGES IMPLEMENTED: Parts `main`, `abs`, `verify` of `run_experiment.py` (fixed seeds per cell, two worker processes, common random numbers across procedures), the per-path prediction in `verify`, the studentised coverage and the misses below and above the target in `summarise`.

EXPERIMENTS PERFORMED: `main`: 32 cells of 2000 paths (normal, t5, t3, jumps; d = 0.4, 0.2, 0.1, 0.05; alpha = 0.05, 0.10), ten procedures each. `abs`: 12 cells. `verify`: 6 cells of 1000 paths (d = 0.1, 0.05, 0.025; with and without intersection).

RESULTS: [SYNTHETIC] (1) Confirmed: CR-trunc coverage 0.936 to 0.964 in all 16 cells at alpha = 0.05 with `E[N]/n*` 0.96 to 1.01; CR-raw 0.72 to 0.78 for t3 at every d and 0.78 for jumps at d = 0.2; Stein-raw down to 0.36. (2) Confirmed: absolute rule, normal returns, d = 0.2: coverage 0.90 with misses 0.088 below and 0.012 above; measured bias -0.004 at d = 0.1 against -0.005 predicted, but -0.036 against -0.021 at d = 0.2 (higher-order terms); relative rule: misses 0.022 below and 0.030 above, bias -0.002. (3) Confirmed: without intersection `CS-opt` gives 2.362, 2.109, 2.007 against predictions 2.341, 2.104, 2.005; with intersection `CS-bet` gives 1.702, 1.576, 1.545 against the limit 1.548. Coverage of every confidence sequence at least 0.969 (nominal 0.95) and 0.939 (nominal 0.90).

IMPROVEMENTS: The question "what does anytime validity cost?" now has a verified answer: about 1.55 times the oracle at small d (1.71 at alpha = 0.10) for the betting sequence, 1.6 to 1.8 at d = 0.1 to 0.2, and the formula of README 2.9 for the empirical Bernstein sequences at any d.

REMAINING ISSUES: The pilot-size study (ALT6-5) found a stall of the betting sequence with a 10-return pilot. Clustering (GARCH) and real data not yet run.

NEXT ITERATION: Pilot size and the grid of candidate means.

---

ITERATION NUMBER: ALT6-5 (a 10-return pilot can leave the truth off the grid)

OBJECTIVE: Measure the effect of the pilot size m = 10, 30, 100 (normal, t3, jumps; d = 0.2 and 0.1).

PROBLEM IDENTIFIED: With m = 10 the betting sequence had mean ratios of 8.4 (normal) and 15.9 (t3) at d = 0.2, because 1 to 6 of 2000 paths per cell never stopped (`results/d_pilot_core_grid_only.csv`). A 10-return pilot can underestimate the scale badly; then the true mean of Z lies above the top of the grid `5 zhat_0`, every grid point is rejected by the lower side, the lower bound stops at the top grid point, and the width can never shrink below the gap.

HYPOTHESIS: Extending the grid of candidate means to the whole of (0, 1) removes the stall. The cap on the lower bet (`0.75 / (5 zhat_0)`) makes factors negative for candidate means above the core, so the extension needs a per-point cap `min(lam, 0.75 / m)`, under which each log-capital term is still nonincreasing in m (so the confidence set is still an interval) but convex instead of concave (so no interpolation above the core). Widening the core instead (span 20) also removes the stall but lowers the lower bet cap and costs 14% at d = 0.2 (2.10 against 1.85 in a check on 1000 normal paths).

CHANGES IMPLEMENTED: 8 log-spaced grid points down to 1e-6 and 32 up to 1 - 1e-6; per-point cap of the lower bet above the core; interpolation only where no past bet was capped; a test with a pilot 20 times too small (the truth outside the core): coverage at all times and every path stops.

EXPERIMENTS PERFORMED: The check on 2000 normal paths with m = 10, d = 0.2; the full `pilot` part; reruns of every part that uses the betting sequence (same seeds).

RESULTS: [SYNTHETIC] With m = 10 at d = 0.2 the worst path now stops at 4.0 times `n*` (mean ratio 1.88), and paths whose bounds stay inside the core are unchanged (1.850 and 1.580 at m = 30, d = 0.2 and 0.05, identical before and after). Full pilot results in `RESULTS.md`, R6.

IMPROVEMENTS: The betting sequence no longer depends on the pilot being roughly right; it only gets slower when the pilot is far off.

REMAINING ISSUES: The pilot is used only for `tau`, the unit and the starting values; its returns are not used by the sequence (sample splitting), which costs `m` returns.

NEXT ITERATION: Volatility clustering.

---

ITERATION NUMBER: ALT6-6 (volatility clustering: which target, and the running intersection)

OBJECTIVE: Run every procedure on GARCH(1,1) returns with persistence 0.95 and 0.99 (d = 0.2 and 0.1, alpha = 0.05, pilot 30, 1000 paths) and find out what, if anything, each procedure covers.

PROBLEM IDENTIFIED: With clustering there are two different targets. The stationary level `theta_tau` is what one usually means by "the volatility of the asset", but the mean of `n` truncated squares has a long-run variance 8 to 32 times the i.i.d. one, so every rule that uses the i.i.d. variance stops too early for it. The confidence sequence's theorem covers a different quantity, a weighted average of conditional means, which moves with time.

HYPOTHESIS: (1) The i.i.d. rules, classical and anytime-valid, fail for the stationary level by a factor close to the long-run variance factor. (2) Batch means fixes the classical rule for the stationary level, at the cost of a minimum of 10 batches. (3) The classical rule is asymptotically valid for the average conditional truncated variance over its own window, by the martingale central limit theorem, because the i.i.d. variance estimate is at least the average conditional variance. (4) The running intersection, harmless when the target is fixed, can break the confidence sequence when the target moves.

CHANGES IMPLEMENTED: The engine records, for dependent laws, the average conditional truncated variance over the returns used (exact, from the simulated conditional variances) and the bet-weighted averages of Theorem 1; the stationary `theta_tau` and the long-run variance for every `tau` are computed once from long simulated paths and cached (`results/_garch_table_*.csv`); `run_experiment.py garch` runs `CS-opt` and `CS-bet` with and without the intersection; a `CR-BM-trunc` rule with batches of 200 and 1000 returns.

EXPERIMENTS PERFORMED: The `garch` part (32 cells), and a check of the batch-means long-run variance of raw squares against the closed form.

RESULTS: [SYNTHETIC] (1) Confirmed: for the stationary level CR-trunc stops at 1/7 to 1/32 of the long-run `n*` and covers 0.32 to 0.56; CR-raw 0.29 to 0.52; the confidence sequences 0.28 to 0.82. (2) Confirmed with a caveat: batch means covers 0.91 to 0.99 at 0.82 to 1.76 times the long-run `n*`, that is 2,000 to 19,000 returns, partly because of its 10-batch minimum; its long-run variance estimate of raw squares is 20.5 against 22.7 and 93.0 against 96.2. (3) Confirmed: CR-trunc covers the window average at 0.948 to 0.961 and CR-raw at 0.955 to 0.966. (4) Confirmed: with the intersection the sequences cover the window average at only 0.68 to 0.79 at d = 0.1 (0.92 to 0.97 at d = 0.2); without it 0.986 to 1.000, and Theorem 1's own weighted target 0.997 to 1.000, at 1.5 to 2.7 times the stopping time of the intersected version.

IMPROVEMENTS: The study now says precisely what can be promised under clustering: an accurate statement about the recent window's average truncated variance, classically or with a guarantee (without intersection), but not about the stationary level within hours. README 2.10 states the martingale argument and the moving-target failure of the intersection.

REMAINING ISSUES: No long-run variance inside the confidence sequence (a guarantee for the stationary level would need a mixing assumption, and the measured costs suggest it would not be useful intraday). The bets are not re-tuned for the non-intersected sequence.

NEXT ITERATION: Real data.

---

ITERATION NUMBER: ALT6-7 (real BTCUSDT: minutes by hour, zero returns, and an early stop that is too good)

OBJECTIVE: Measure, on Binance BTCUSDT one-minute returns (July to December 2025, divided by the past-only intraday profile), how many minutes each rule needs for d = 0.5, 0.4, 0.3, 0.2 by start hour, how often the raw rule fails to stop within the hour, and check the procedures on the real shape of returns with a stationary bootstrap of quiet days.

PROBLEM IDENTIFIED: The first real run crashed with "zero pilot scale": 4.3% of one-minute returns are exactly zero, and in 1 of 4392 windows more than half of the 30 pilot returns were zero, so the scaled MAD was zero. Second, the first complete run showed the intersected betting sequence needing only 1.3 times the minutes of truncated Chow and Robbins at d = 0.2 (median 284 against 220), where i.i.d. normal returns give about 1.85 (R3): on real data it looked cheaper than the theory allows.

HYPOTHESIS: (1) The package's own fallback (median of the nonzero absolute returns, as in `VolatilityMonitor._fallback_scale`) is the right repair and affects almost no window. (2) The low ratio of the intersected sequence is the moving-target effect of ALT6-6, not a gain: without the intersection the ratio should be above the i.i.d. value.

CHANGES IMPLEMENTED: `seqest.pilot_scale` with the zero-median fallback, used by every procedure; the real part also runs `CS-opt` and `CS-bet` without intersection, records an in-window `n*` proxy, the share of zero returns and the fallback share; Figure 3 shows median minutes by start hour at d = 0.4 (about plus or minus 20% on volatility) and d = 0.2 (plus or minus 20% on the variance), and the share not stopped within 60 minutes. A separate diagnostic run (not saved) looked at the windows that never stop.

EXPERIMENTS PERFORMED: The `real` part (4392 windows, both alpha), the `boot` part (30 quiet days, mean blocks 1 and 30, 200 resamples per day), Figure 3.

RESULTS: [REAL] (1) Confirmed: fallback in 0.02% of windows. Minutes needed at 95%, for about plus or minus 20% on volatility: CR-trunc 50 (46 to 52 by start hour), Stein 55, CR-raw 66 (52 to 82), CS-bet without intersection 140; for plus or minus 20% on the variance: 220, 220, 411 and 505. The raw rule has not stopped within its start hour in 58% of hours at d = 0.4 (42% to 69% by hour) and 21% at d = 0.5, against 27% and 1.7% for CR-trunc, and fails to stop within 24 hours in 6.0% of windows at d = 0.2. (2) Confirmed: without intersection CS-bet needs 2.20 times the proxy at d = 0.2 (1.32 with it); the block bootstrap shows the same drop (1.74 with i.i.d. resampling, 1.47 with 30-minute blocks), and under GARCH the intersected sequence stops at 277 against 448 returns. The windows that never stop (1.9% for CS-bet without intersection at d = 0.2) start mostly on Friday afternoons (50 of 85) and run into quiet weekend hours with 7.5% zero returns and an in-window `CV^2` of 7.0. [REAL, resampled] i.i.d. resampling of quiet days: CR-trunc 0.949, CS-bet 0.979, CR-raw 0.919; 30-minute blocks: 0.70 to 0.82 for every procedure.

IMPROVEMENTS: The real-data answer is now stated with the version of the confidence sequence that remains valid under clustering, and the apparent real-data advantage of the intersected sequence is explained rather than reported as a gain.

REMAINING ISSUES: No coverage on real data can be claimed. The deseasonalisation removes the intraday pattern but not the weekly one (quiet weekends). The pilot's returns are not used by the sequence. One coin, one venue, six months.

NEXT ITERATION: None planned within this prototype; the open items are listed in `RESULTS.md` (ranking) and in the report part.
