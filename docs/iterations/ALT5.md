# ALT5 iterations: two-sample betting monitor (conformal rank e-detector)

Prototype folder: `alternatives/two_sample_betting/`. Times are UTC, night of 2 November 2026.

---

ITERATION NUMBER: ALT5-1 (formulation and design pilot)

OBJECTIVE: Turn entry A4 of `alternatives/ALTERNATIVES.md` (a sequential two-sample test by betting between the recent window and the calibration window) into a working change detector with an exact false alarm guarantee and no estimated baseline.

PROBLEM IDENTIFIED: The paired two-sample betting test of Shekhar and Ramdas (2024) needs a fresh reference draw per step, but the monitor has one finite calibration window. Two obvious adaptations fail: drawing reference points without replacement runs out after `n0` steps and doubles the noise of every payoff; ranking new returns against the fixed reference (quantile bins learned on the calibration window) gives p-values that are uniform one at a time but dependent, so a learning bettor is not betting fairly.

HYPOTHESIS: Ranking each new return among ALL earlier returns of the segment (calibration window plus monitored returns, a growing reference) gives smoothed conformal p-values that are exactly i.i.d. uniform under exchangeability (Vovk, Nouretdinov and Gammerman 2003). Any bet that is fair against the uniform law and uses only past p-values then gives a Shiryaev-Roberts e-detector with exactly the guarantees of `docs/math.md` 5.4. Mixing a rank analogue of the main monitor's truncated-variance bet (Gaussianised score, null mean known to be 1) with a histogram learner (plug-in density bettor) should cost little on variance changes and gain a lot on shape changes.

CHANGES IMPLEMENTED: `rankbet.py` (vectorised exact p-values with a binary indexed tree, streaming `SequentialRanker`, Gaussianised score, three bettor families, null band `b`, Shiryaev-Roberts mixture, streaming `RankMonitor`), `tsb_sim.py` (unit-variance laws, GARCH, a vectorised copy of the plug-in `VolatilityMonitor`), `run_experiment.py` part `pilot`, `tests/test_rankbet.py`.

EXPERIMENTS PERFORMED: Pilot on development seeds 1000 to 1002 (never reused): seven configurations (memory of the histogram learners 60/240, 240/1440, 1440/5760 steps; prior 50 or 200; family weights; scale bets only; histogram only) on four changes (200 paths) and two one-day nulls (400 paths). Output `results/p_pilot_design.csv`.

RESULTS: [SYNTHETIC] Scale bets alone were fastest for pure variance changes (median delay 149.5 against 167.5 steps for the default mixture for sd x1.25) but slow for Gaussian to t4 at fixed variance (1348, 74% detected); the histogram learner alone was 1.5 to 2 times slower on variance changes but fast on t4 (387). The mixtures C0 (default), C1, C3 and C4 were within about 10% of each other on every row. Null false alarms within one day: 1.0% to 4.0% (bound 5%). Under GARCH(1,1) with persistence 0.95 every configuration with scale bets alarmed on 70% to 77% of days; the histogram learner alone on 25%.

IMPROVEMENTS: C0 (scale up 0.25, scale down 0.25, histogram 0.5; half-lives 240 and 1440; prior 50) was frozen as the default before any reported run: it was never the worst mixture and never more than 8% behind the best mixture on any change.

REMAINING ISSUES: The pilot cannot tell whether real data are close enough to exchangeable for the exact null to be useful.

NEXT ITERATION: Unit tests of the proofs' ingredients, then the reported synthetic study.

---

ITERATION NUMBER: ALT5-2 (tests of the validity ingredients; two failed first versions)

OBJECTIVE: Check numerically every step of the validity proof and the equivalence of the benchmark with the repo's monitor.

PROBLEM IDENTIFIED: Two first test versions failed. (1) A test asserted that a Gaussian to t3 change (same variance) would be carried by the histogram learner; the alarm was carried by the scale-down bets (92% of `M_t`). (2) A test of the band null checked `E_h[f] <= 1` on `h = 1 + b sign(w)`, which is not a probability density, and failed.

HYPOTHESIS: (1) is not a bug: a variance-matched t3 puts more mass near zero, so the Gaussianised truncated variance falls (the same effect the quarticity prototype reports for the truncated variance), and the decrease bets pick it up. (2) is a wrong test: the band proof only needs `E_h[w] <= b E_U|w|` for densities `h`.

CHANGES IMPLEMENTED: (1) The test now checks detection soon after the change and that the family shares add up to 1. (2) The test uses the most favourable valid density `h = 1 + b sign(w - median w)`, renormalised.

EXPERIMENTS PERFORMED: 9 tests, about 3 seconds: p-values against brute force with ties, streaming against vectorised, uniformity and zero lag-one correlation for i.i.d. Poisson data, `E_U[Y] = 1`, positivity and fairness of every bet for a learner state driven far from uniform (bands 0 and 0.3), `E[M_t] = t` under the null, the vectorised plug-in monitor against `VolatilityMonitor` alarm for alarm, detection of a tail change, NaN handling and configuration checks.

RESULTS: All pass. The p-values agree with brute force exactly, and for i.i.d. Poisson scores (heavy ties) the mean, variance, 5% quantile and lag-one correlation of the p-values are within the test tolerances (0.005, 0.002, 0.004 and 0.01) of the uniform values.

IMPROVEMENTS: Each line of the proof has a test.

REMAINING ISSUES: None for validity under exchangeability.

NEXT ITERATION: Q1 to Q3, the synthetic comparison with the mean-based monitor.
