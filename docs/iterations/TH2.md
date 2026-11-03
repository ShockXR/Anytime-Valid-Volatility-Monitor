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

**EXPERIMENTS PERFORMED.** Unit tests: p-values equal a brute-force computation with many ties; uniformity (Kolmogorov-Smirnov) and lag-1 independence for Student t3 and for coarse ticks (about 20% of values tied at zero); betting functions integrate to 1; exact false alarm bound under ties and jumps; agreement between the detector and the p-value function.

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

## Iteration 4. Full E14 runs, and a negative result for the plug-in

**OBJECTIVE.** Validity at scale, delay, robustness, real data and cost for every method (E14 Q14a to Q14g, Q14w), with the settings fixed in iterations 1 to 3.

**PROBLEM IDENTIFIED.** E8 judged the plug-in on one day of monitoring only, where it sits at its nominal 5%.

**HYPOTHESIS.** (1) Route B uses its budget and stays within it; route A and E8 are conservative. (2) Over longer horizons the plug-in's error does not average out, because one estimated mean is used for the whole segment: a window that underestimates the mean gives the increase bets a positive drift, a window that overestimates it does the same for the decrease bets.

**CHANGES IMPLEMENTED.** The validity runs report five-day as well as one-day rates. Part w was extended with `plugin_drift` (the share of calibration windows with positive null drift, by quadrature on 20,000 quantile points of 4 million draws, and the run length bound). The theory note gained Proposition A.6 with a proof: for a fixed calibration length, `E[T_A | W] ≤ (log(2AK) + log f_max) / g(W)` whenever the window W gives some bet a positive drift g, so no bound of the form C m / A can hold for the plug-in if such windows have positive probability (proved later, in iteration 6, as Lemma A.7; here it was only checked numerically). The labels of the rank-bet extremes in `e14w_constants.csv` were corrected (for ρ < 1 the largest factor is at u = 1, not at u0), and the run was repeated, with identical bands and drift files.

**EXPERIMENTS PERFORMED.** Q14a (5000 paths per cell), Q14b (1000), Q14c (1000), Q14f (4000 null, 1000 per delay cell), Q14w (2000 windows), Q14d (six real series, 15 monitor settings), Q14e.

**RESULTS.** (Numbers from `results/e14*.csv`; details in `docs/report_parts/E14.md`.) Hypothesis 2 confirmed. The plug-in alarms on 40.6% to 42.9% of null paths within five days after 390 calibration steps, against the 25% bound; 77% to 83% of windows give it a positive drift. Hypothesis 1 confirmed: B-lr 1.7% to 3.2% per day; route A and cs at most 0.26% over five days. Unexpected: B-lr is faster than the plug-in on every decrease (by a factor of two to three for x0.5). A bet-design confound is noted in the report. The awake route A monitor is valid but cannot detect changes in the first few dozen steps of a segment (iteration 3); awake-B can.

**IMPROVEMENTS.** A proved negative result that explains a numerical one, and a clear recommendation: the rank route, restarted with no warm-up, for a monitor with no blind period.

**REMAINING ISSUES.** The coordinator asked for blind time and false alarms per day on streams, and for the E12 releases.

**NEXT ITERATION.** Streaming measurements.

---

## Iteration 5. Streaming: blind time, false alarms per day, macro releases

**OBJECTIVE.** Measure the share of time blind, false alarms per day under a true null, and delays on streams; rerun E12's streaming design with the awake monitors.

**PROBLEM IDENTIFIED.** Q14f measures one segment after one restart; a stream with repeated alarms and restarts is what a user runs.

**HYPOTHESIS.** In a stream whose volatility moves about once a day, the standard protocol is blind most of the time and misses many changes. The awake monitors are almost never blind; their false alarm rate per day stays within the Corollary A.5 bound under a true null.

**CHANGES IMPLEMENTED.** E14 Q14h (`h_job`, `regime_path`: 30-day minute-scale streams, null and regime-switching, five streaming monitors); Q14g (`g_job`: E12's data, events, control windows and stratified test, with the monitor swapped); an option to set the pool size (`--procs`) so the runs respect the two-process limit.

**EXPERIMENTS PERFORMED.** Q14h, 40 paths of 30 days per cell; Q14g, three coins and four monitors.

**RESULTS.** Blind time: the standard protocol is blind 57% to 71% of the time, the awake monitors 0.28%. Under the null, awake-B raises 0.029 to 0.035 false alarms per day (bound 0.050) and awake-A none; the plug-in raises 0.084 to 0.094 per day, with 95% lower limits above its nominal 0.05, which confirms iteration 4 on streams. Under regime switching, awake-B detects 98.0% to 99.6% of the changes within a day (median delay 34 to 60 steps), against 56% to 66% for the standard protocol. E12 rerun: the plug-in column reproduces E12 (awake at 9, 8, 9 FOMC and 9, 8, 12 CPI releases). The awake monitors watch all 114 releases. Up alarms within 6 hours rise from 21 to 25 (awake-A) and 31 (awake-B) of 66 CPI releases, and from 26 to 30 and 35 of 48 FOMC statements. BTC CPI goes from p = 0.11 to 0.003 and 0.0003. (The awake-A numbers changed slightly in iteration 6: 24 instead of 25 CPI releases.) The awake monitors flag a smaller share of the releases at which they are awake than the plug-in does. A possible reason is a restart shortly before a release; this was not tested.

**IMPROVEMENTS.** The no-blind-period deliverable now has synthetic and real evidence.

**REMAINING ISSUES.** Without a null band, awake-B alarms 3.6 to 4.3 times per crypto day. The rank reference set costs O(n) per insertion when no alarm comes for a long time. The bet-design confound in the B-lr speed comparison remains.

**NEXT ITERATION.** None planned in this task; needed changes in files owned by others are listed in the handback.

---

## Iteration 6. Review round 3: Lemma A.7, an "off" decrease side, wording

**OBJECTIVE.** Fix the E14 findings of `docs/review/independent_review_round3.md`: M1, m1 to m8.

**PROBLEM IDENTIFIED.** (M1) Proposition A.6 stated a limit where its proof gives only a lower limit, and its hypothesis `p_{η,κ} > 0` was supported by simulation only, while the status table of the theory note called the result proved. (m4) A decrease side with no admissible bet (lower bound 0, or τ_t² at or below it) grew by 1 per step instead of staying at 0 as in `docs/math.md` 6.3; in the awake route A monitor, the first snapshot (no data, lower bound 0) added a constant 16 to the decrease statistic. (m1 to m3, m5 to m8) Wording, missing conditions and one missing test.

**HYPOTHESIS.** The reviewer's derivation D1 proves `p > 0` for the four E14 laws. The new decrease convention keeps every proof, because an unchanged statistic still satisfies `E[M_t | F_{t-1}] ≤ M_{t-1} + 1`, and changes numbers only where awake-A alarms near its threshold.

**CHANGES IMPLEMENTED.** Theory note: Proposition A.6 now states `liminf P(T_A ≤ m_A) ≥ p_{η,κ}` (and notes that the lower limit is at least `P(g(W) > 0)`); new Lemma A.7 with a full proof. Every step of D1 was checked and holds: the event that all calibration returns lie in (ε/2, ε], the range of B̂ (between a²/4 ≈ 4.946 and 4a² ≈ 79.13), the drift bound for the bet 0.5, the largest factor below 41, the general case (no atom at 0, 0 in the support) and the tick case (drift 0.472, largest factor 10.39). The status table now lists Proposition A.6 and Lemma A.7 as proved and the size of the share as numerical. Remark B.3 states exchangeability given the returns up to the alarm (m2); section 4.3 separates the unclipped likelihood-ratio rate (the Kullback-Leibler divergence) from the clipped one (m1; 0.573 against 0.807 at ρ = 2, both checked by numerical integration). `running.py`: an off decrease step leaves the decrease statistic unchanged in `two_sided_sr`, `RunningTauMonitor`, `awake_route_a` and `AwakeRunningTauMonitor`. Tests: an off decrease side stays at 0 (vectorised, streaming route A, awake route A), mixed on and off steps against a direct recursion, the streaming rank monitor against `conformal_detector` alarm by alarm across restarts with the same tie-breakers (m3), and the arithmetic of Lemma A.7 on simulated windows. Report: M1 wording, m5 to m8, and the m18 and m19 qualifiers (a known baseline would guarantee 25%; the default plug-in alone catches 65.4% to 65.9% of the changes).

**EXPERIMENTS PERFORMED.** With the new code, every E14 cell that uses the changed functions was rerun (scratch script, nothing written to `results/`) and compared with the CSVs: the X-based methods of Q14a, Q14b, Q14c and of Q14f's recalibrating monitors, awake-A in Q14f, Q14d (route A and awake-A), Q14g and Q14h. Where awake-A differed, the E14 parts were rerun with the script: Q14d, the delay part of Q14f, Q14g and Q14h.

**RESULTS.** The plug-in, cs and route A results are identical (no path of Q14a or Q14c had an off decrease step; 46 paths of one Q14b cell had one, without effect). The awake-A null rates of Q14f are identical. awake-A changed a little: one median delay in Q14f (966 to 967 steps); in Q14h, Student t3 regime switching, 77.8% instead of 77.9% of changes detected (x0.67 median 314 instead of 324, x1.5 348 instead of 343); in Q14d a few alarms (BTC band 0.25: 0.76 instead of 0.77 per day); in Q14g one CPI release fewer with an up alarm (24 instead of 25 of 66; SOL p 0.0059 instead of 0.0013, BTC p 0.0026 instead of 0.0032). No other monitor's numbers changed. The real-data costs per return of Q14d were remeasured on a busier machine and are higher. One point of the review was corrected: for m5 the reviewer gave the growth at k = 1940 as -0.001 to 0.013, but over all six ratios `results/e14w_dilution.csv` gives -0.028 to 0.013.

**IMPROVEMENTS.** The negative result for the plug-in is now proved for all four simulated laws, not only checked numerically; the decrease side matches `docs/math.md`.

**REMAINING ISSUES.** Lemma A.7 proves only that the share of drifting windows is positive; its size (77% to 83% at 390 calibration steps) remains a numerical estimate. The figures `results/figures/e14_*.png` were not redrawn after the reruns; the changed awake-A numbers are too small to be visible in them.

**NEXT ITERATION.** None planned.
