# TH2 iteration notes: a guarantee with a predictable, running truncation level (tracker B68)

Owner: TH2. Files: `src/avmonitor/running.py`, `tests/test_running.py`, `experiments/e14_running_tau.py`, `results/e14_*.csv`, `results/figures/e14_*.png`, `docs/theory_running_tau.md`, `docs/report_parts/E14.md`, this file.

---

## Iteration 1. Route A: a band for all truncation levels

**OBJECTIVE.** Let the truncation level τ_t follow the current scale (any predictable rule) and keep a proved false alarm bound, by estimating the whole function `s ↦ m_F(s) = E_F[min(r², s²)]` with a band that holds for all s at once.

**PROBLEM IDENTIFIED.** `docs/math.md` section 6 (E8) covers one level fixed before its confidence sequence starts; section 6.4 lists a running τ as not covered. The plug-in with a running τ has no proof at all.

**HYPOTHESIS.** (1) A band on the distribution function of r² integrates to a band on m(s) for every s (Lemma A.2), so the time-uniform DKW band of Howard and Ramdas (2022) gives a guarantee for any predictable τ_t. (2) That band has the same width at every x, so after integration up to s² = 9 m(s) it will be much wider than a confidence sequence for one level. (3) A union of per-level empirical Bernstein confidence sequences on a geometric grid of levels (13 levels, δ / 26 each) will be close to E8's single-level width, because the union bound costs `log(2600)` against `log(200)` while no sample split is needed.

**CHANGES IMPLEMENTED.** `running.py`: `hr_dkw_constant` and `hr_dkw_radius` (the constant C solved numerically exactly as in the authors' library, github.com/gostevehoward/confseq), `massart_dkw_radius`, `truncated_means_from_band` (closed form of Lemma A.2), `calibrate_band` (methods "grid", "dkw", "plugin"), `running_level_index` (c × scaled MAD of the last 120 returns, refreshed every 10 steps, snapped to the grid), `two_sided_sr` (Shiryaev-Roberts with per-step truncation, null means and decrease bets), `route_a_inputs`, the streaming `RunningTauMonitor`.

**EXPERIMENTS PERFORMED.** Pilot on seeds not used later (scratch script `pilot1.py`, 400 Gaussian windows): median upper bound / true m at the centre level and at other levels, for 390 and 1440 calibration steps.

**RESULTS.** The pilot numbers are not reported as results; they only guided the design. The DKW band was about twice as wide as the grid band at the centre level, and the grid band was slightly wider than E8's single-level bound at the same length (the union bound costs a little more than the sample split saves). The grid band's lower bound at low levels was tighter than at the centre, which is what the decrease side needs. The final measured widths are in `results/e14w_bands.csv` (E14, Q14w).

**IMPROVEMENTS.** A proved bound for any predictable τ_t (Theorem A.1), under a null that restricts only conditional truncated means.

**REMAINING ISSUES.** The DKW route is too wide to be useful; it is kept because it is the cleanest "uniform in s" statement and is cited by the brief. The upper bound at high levels is loose (range problem, theory note 3.5).

**NEXT ITERATION.** Route B.

---

## Iteration 2. Route B: exact rank-based detector

**OBJECTIVE.** A statistic whose null law is known exactly, so that no baseline is estimated.

**PROBLEM IDENTIFIED.** Every X-based method pays for estimating m(τ): the plug-in with validity, E8 and route A with width.

**HYPOTHESIS.** Smoothed conformal p-values of |r_t| among all earlier values are i.i.d. uniform under exchangeability (Vovk, Nouretdinov and Gammerman 2003). Betting functions with mean 1 then give an exact Shiryaev-Roberts e-detector. A likelihood ratio in rank space (normal scale family), clipped at P(|Z| > 3), should be near-optimal for normal data and robust to single jumps.

**CHANGES IMPLEMENTED.** `dense_ranks`, a per-row Fenwick tree (`_RowFenwick`) for vectorised exact ranks with ties, `conformal_pvalues`, `lr_bet_factory`, `indicator_bet_factory`, `bet_family`, `family_factors` (one-pass evaluation; `scipy.special.ndtri` instead of `stats.norm.isf` made the streaming monitor several times faster), `conformal_detector` (with warm-up, null band f / (1 + b), optional heuristic freeze), streaming `ConformalMonitor`.

**EXPERIMENTS PERFORMED.** Unit tests: p-values equal a brute-force computation with many ties; uniformity (Kolmogorov-Smirnov) and lag-1 independence for Student t3 and for coarse ticks (38% ties at zero); betting functions integrate to 1; exact false alarm bound under ties and jumps; agreement between the detector and the p-value function.

**RESULTS.** All tests pass. The test bound of one test was first written as "the clipped bet is below 15"; the actual maximum factor for ρ = 2 is about 16 (the test was wrong, not the code, and now checks the exact value `g(u0) / C`).

**IMPROVEMENTS.** A detector with an exact guarantee and no tuning constant except the betting functions.

**REMAINING ISSUES.** Dilution after a change (theory note 4.3) and the fact that exchangeability is violated by volatility clustering.

**NEXT ITERATION.** The coordinator asked for a monitor with no blind period after an alarm (E12 found the streaming monitor blind at 13 of 22 CPI releases).

---

## Iteration 3. No blind period

**OBJECTIVE.** Resume monitoring immediately after an alarm, with a proof.

**PROBLEM IDENTIFIED.** After every alarm the standard protocol calibrates for 1440 minutes; at about 0.7 alarms per day the monitor is blind about 70% of the time (E12).

**HYPOTHESIS.** (A) A band that learns from the post-alarm returns, with each Shiryaev-Roberts chain betting against the snapshot of the band taken before the chain starts, is valid on the event that the band has not failed by the chain's start, an event known at that time; stopping each chain on that event keeps it a supermartingale without conditioning on the unobservable "the band never fails" (Theorem A.4). (B) Route B needs nothing new: restart the reference set and bet from the first return (Remark B.3).

**CHANGES IMPLEMENTED.** `awake_epochs` (snapshots at 0, 16, 23, 32, ..., 1024, 1440 observations), `_LevelsCS` (the E8 confidence sequence for 13 levels at once, checked to reproduce `eb_confseq` exactly), `awake_band_snapshots`, `awake_route_a` (vectorised), `AwakeRunningTauMonitor` (streaming), `ConformalMonitor(restart_len=0)`. Pilot (`pilot2.py`, other seeds): the bet schedule of the learning sequence. Tuned for 1440 observations the upper bound is loose early; the time-uniform schedule makes the lower bound loose late; `n_plan = 300` was chosen as a compromise before any reported run.

**EXPERIMENTS PERFORMED.** Smoke run of E14 Q14f with 200 to 400 paths (not reported).

**RESULTS.** The smoke run already showed the main finding, confirmed by the full run (E14 report): route A with a learning band is valid but nearly powerless when the change comes within the first few dozen steps of a segment, because a band from a few dozen returns is very wide; later chains use bands that already contain post-change returns. Route B restarted without warm-up detects such changes.

**IMPROVEMENTS.** Two proved monitors with no blind period; one of them (awake-B) is also fast.

**REMAINING ISSUES.** See iteration 4.

**NEXT ITERATION.** Full E14 runs.

---

## Iteration 4. Full E14 runs

(Filled in after the runs; see below.)
