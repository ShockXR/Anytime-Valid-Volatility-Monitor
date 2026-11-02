# Independent review, round 3

**Reviewer role.** An adversarial reader: a statistics journal referee and a code auditor. The goal is that every claim in the files below survives an expert reading.

**Scope.**

1. E13, benchmarks against published methods:
   * `docs/report_parts/E13.md` and `docs/iterations/E13.md`;
   * `src/avmonitor/competitors.py` and `tests/test_competitors.py`;
   * `experiments/e13_benchmarks.py` and `results/e13*.csv`.
2. E14, a running truncation level, the rank (conformal) detector and monitoring with no blind period:
   * `docs/theory_running_tau.md`, including Proposition A.6;
   * `docs/report_parts/E14.md` and `docs/iterations/TH2.md`;
   * `src/avmonitor/running.py`, `experiments/e14_running_tau.py` and `results/e14*.csv`.
3. How E13 and E14 are stated elsewhere:
   * `README.md`, the bullets "Estimated baselines need care", "A monitor with no blind period" and "Against published methods" (lines 26 to 28), and line 30;
   * `PITCH_BAPAT.md` sections 6 to 9;
   * `docs/report_parts/_header.md`, rows E13 and E14.

**Read first.** `docs/math.md` sections 5 and 6, and the earlier reviews in `docs/review/`. Their points are not repeated.

**State of the tree.** The review is against commit 2387173. Commits up to 77662ba were made while it ran. None of them touched a file in scope (they changed `results/e11f_*` and files of other experiments), so the line numbers below hold at 77662ba.

**What I ran.** Every number in `E13.md` and `E14.md` was compared with its CSV, and every range, ratio and test was recomputed with pandas. The same was done for the E13 and E14 statements in the README, the pitch and the report index. Small checks, with scripts in `/tmp/claude-0/-home-user/e30ce406-b0e8-5b23-b612-d49eb104d0ee/scratchpad/r3/` and nothing written to `results/`:

* `check_e13_ticks.py`: reruns the E13 offline null job for the tick law with the script's own seed (66 s). Every method's matched rejection rate equals `results/e13a_offline_size.csv` exactly.
* `check_e14_drift.py`: reruns the Q14w plug-in drift cell (Gaussian, 390 calibration steps) with the script's seeds.
  * It reproduces the CSV: shares 0.406, 0.4255 and 0.8315.
  * It also compares the quadrature with exact numerical integration (scipy `quad`) on 60 random windows. The sign of the best increase and decrease growth agrees in 60 of 60.
* `check_restart.py`: the rank monitor restarted with an empty reference set after each alarm (`ConformalMonitor`, `cal_len = 0`, `restart_len = 0`, A = 100), 30 streams of 20,000 steps per law (D5).
* `check_bjv_k.py`: offline power of the overlapping BJV-style statistic at simulated 5% size for window lengths k = 50, 125, 250 and 500 (D3).
* `check_bocpd_cap.py` and `check_bocpd_cap2.py`: the capped BOCPD filter against the exact recursion without a cap (D4).
* Closed forms by hand (section 5).
* `python -m pytest -q tests/test_competitors.py tests/test_running.py`: **35 passed in 5.7 s.**

No file in the repository was changed except this one.

**Severity.**
* *Critical*: a wrong proof step or a wrong headline number.
* *Major*: an overclaim, a missing condition or a misleading comparison.
* *Minor*: a wrong number outside the headlines, wording or clarity.

---

## 1. Summary verdict

**The proofs hold, and no headline number is wrong. I found no critical problem.** I re-derived the following and found no error in any proof step:

* in `theory_running_tau.md`:
  * Lemma 2.1 and (2.1), Theorem A.1, Lemma A.2 (the closed forms) and the reduction from continuous laws to laws with atoms;
  * Proposition A.4a, Theorem A.4 (including the run-length part) and Corollary A.5;
  * the proof of Proposition A.6 (optional stopping and Markov's inequality);
  * Proposition B.1 with ties and the randomised p-value, Theorem B.2, and Remark B.3 for i.i.d. segments;
  * the dilution limit and the freeze argument of section 4.3;
* the BOCPD predictive and recursion, and the Gumbel and MOSUM constants used for BJV.

The code implements the stated formulas. The rank detector is exact with ties and after a restart, in theory (D5) and in simulation: 145 to 148 alarms per 20,000 steps against the bound 200.

**Proposition A.6 is correct in substance but misstated, and its hypothesis is only checked numerically (M1).** It says that a probability converges to p, where the proof gives a lower bound. Its conclusion needs p > 0, which the note supports only by simulation, yet the status table, the README, the pitch and the report index call the result proved. D1 gives a short proof of p > 0 for all four E14 laws, so the claim becomes true once that proof is added.

**The other problems are in how E13 is read.**

* **M2.** "Whose bound of 5% is a theorem" is said of the plug-in e-detector, which has no proved bound. The 1.2% to 1.3% quoted in the headlines is also the plug-in's.
* **M3.** BJV's low offline power comes from its short windows. With windows of 500 returns the same statistic finds 99.7% of 25% changes at matched size. So "accumulating evidence wins" is not shown.
* **M4.** The real-data score asks whether the next day stays beyond a one-day forecast. That is the time scale of the e-detector's 1440-minute baseline. The competitors' windows (120 and 250 minutes) were never varied on real data. "The published methods are less precise on crypto" is therefore too general, and "favours the competitors" is not supported.
* **M5.** "Two to three times faster for decreases on clean data" holds in 2 of the 9 clean decrease cells. This errs against the e-detector, but it is a headline in four files.

**Most important fixes before the work is shown.**

1. Restate Proposition A.6 with a lower limit, and add the proof that p > 0 (M1, D1).
2. Correct `E13.md:168` and say which baseline the 1.2% to 1.3% belongs to (M2).
3. Scope the BJV and real-data conclusions to the window lengths that were tried, and add one longer-window run (M3, M4).

---

## 2. Critical findings

None.

* **Proofs.** No proof step is wrong in a way that breaks a stated result. The misstated limit in Proposition A.6 (M1) is not used in its conclusion.
* **Headline numbers.** Every headline number in the README, the pitch and the report index agrees with its CSV. The problems are missing conditions (M2, M5) and readings the data do not support (M3, M4).

---

## 3. Major findings

### M1. Proposition A.6: the limit is only a lower bound, and p > 0 is not proved in the note

* **Where.**
  * `docs/theory_running_tau.md:154`: "`P(T_A ≤ m_A) → p_{η,κ}` as A → ∞".
  * `:158`, "When is p > 0?": a heuristic and a numerical share.
  * `:223`, the status table: "Proved here: ... Proposition A.6 (the plug-in has no bound of the form C m / A)".
  * `README.md:26`: "We proved that no bound of the form C m / A can hold for it".
  * `PITCH_BAPAT.md:40`: "We proved ... because a fixed share of calibration windows gives some bet a positive drift".
  * `docs/report_parts/_header.md:28`: "The plug-in has no bound of the form C m / A (Proposition A.6; ...)".
  * `docs/report_parts/E14.md:16` ("Proposition A.6 explains why"), `:165` and `theory_running_tau.md:158` ("quadrature on 4 million draws"), `E14.md:167` (table header) and `:178` ("proves a false alarm probability of at least 7.6% to 8.0%").
* **Evidence.**
  * **The limit.** The proof gives `P(T_A ≤ m_A) ≥ p_{η,κ} (1 - (log(12A) + κ) / (η (log A)²))`. The right side tends to p_{η,κ}, so the proof shows `liminf P(T_A ≤ m_A) ≥ p_{η,κ}`, not convergence to p_{η,κ} (D2). The conclusion needs only the lower limit, so it stands. This is why I rate the point major and not critical.
  * **The hypothesis.** The proposition is conditional: "if p_{η,κ} > 0 for some η, κ".
    * The note supports this only with a heuristic (line 158).
    * E14 estimates the share by Monte Carlo over 2000 windows, with expectations from quadrature.
    * So "proved" in the status table, the README, the pitch and the report index covers a statement whose hypothesis is only checked numerically.
  * **The gap closes.** D1 proves p > 0 for each of the four E14 laws and every calibration length n. An explicit event of positive probability gives the increase bet λ = 0.5 a drift of at least 0.19, with every factor below 41. The claim is true; the note just does not prove it.
  * **"Explains why" (E14.md:16).** A.6 gives the mechanism, but not the size of the effect. Its quantitative lower bound in the five-day setting of the 41% to 43% (A = 7800, m = 1950, 390 calibration steps) is 0.75% (`results/e14w_plugin_drift.csv`, `markov_lower_bound_false_alarm` = 0.0075 for gaussian, 390, A 7800, horizon 1950). `E14.md:378` says the bound is crude; line 16 should not suggest that it accounts for 41%.
  * **The quadrature.** The expectations are not computed on 4 million draws. `plugin_drift` uses Q = 20,000 quantile points of the 4 million draws (`experiments/e14_running_tau.py:836`). My check against exact integration supports the quadrature (60 of 60 signs agree), so this is a description error only.
  * **The median.** The "median run length bound" column is the median over the windows with positive drift only (`np.median(tb[pos])`). The table header does not say so.
  * **"Proves" (E14.md:178).** The 7.6% to 8.0% is a Monte Carlo estimate (2000 windows, quadrature) of a lower bound, not a proved value.
* **Fix.**
  * `theory_running_tau.md:154`: replace "so with `m_A = (log A)²`, `P(T_A ≤ m_A) → p_{η,κ}` as A → ∞ while `m_A / A → 0`" with "so with `m_A = (log A)²`, `liminf_{A→∞} P(T_A ≤ m_A) ≥ p_{η,κ}` while `m_A / A → 0`".
  * After the proof, add Lemma A.7 with the proof of D1:
    > **Lemma A.7 (positive drift has positive probability).** Let `q_ε = P(|r| < 4.45 ε)`. If `P(ε/2 < |r| ≤ ε) > 0` and `q_ε < 0.61` for some ε > 0, then for every n, `p_{η,κ} > 0` with `η = 1.089 (1 - q_ε) - 0.693 q_ε` and `κ = log 41`. This holds for every law with no atom at 0 whose support contains 0, including the Gaussian, Student t3 and jump laws of E14 (ε = 0.1). For the tick law, the event that every calibration return is ±0.5 gives drift 0.472.

    Then change the status line to "Proposition A.6 and Lemma A.7".
  * `theory_running_tau.md:158` and `E14.md:165`: replace "by quadrature on 4 million draws" with "by quadrature on 20,000 quantile points of 4 million draws". `E14.md:167`: write "median run length bound (windows with positive drift)".
  * `E14.md:16`: replace "Proposition A.6 explains why" with "Proposition A.6 gives the mechanism".
  * `E14.md:178`: replace "proves a false alarm probability of at least 7.6% to 8.0%" with "gives an estimated lower bound of 7.6% to 8.0% on the false alarm probability".

### M2. "Whose bound of 5% is a theorem" is said of the plug-in, which has no proved bound

* **Where.**
  * `docs/report_parts/E13.md:168`: "(7.5% per day under t5, against 1.2% for the e-detector, whose bound of 5% is a theorem)".
  * The headlines that quote the plug-in rate next to BOCPD's without naming the baseline:
    * `README.md:28`: "while the e-detector stays at 1.2% to 1.3%";
    * `PITCH_BAPAT.md:35`: "our monitor stays at 1.2% to 1.3%";
    * `docs/report_parts/_header.md:27`: "the e-detector stays at 1.2% to 1.3%".
* **Evidence.**
  * **Which detector gives 1.2%.** In Q13b (`E13.md:67`, `results/e13b_online_null.csv`) the 1.2% under Student t5 and the 1.3% under jumps are the plug-in's. The CS baseline gives 0.0% and 0.1%.
  * **The report contradicts itself.** `E13.md:63` says: "The e-detector's bound is 5% (proved for the CS baseline; for the plug-in only with the true baseline)". E14 (Proposition A.6) shows that the plug-in has no bound of the form C m / A at all. So line 168 attaches the theorem to the wrong detector.
* **Fix.**
  * `E13.md:168`: "(7.5% per day under t5, against 1.2% for the plug-in e-detector, which has no proved bound, and 0.0% for the CS-baseline e-detector, whose 5% bound is a theorem)".
  * `README.md:28`, `PITCH_BAPAT.md:35`, `_header.md:27`: "the e-detector stays at 1.2% to 1.3% (plug-in baseline, measured; 0.0% to 0.1% with the confidence-sequence baseline, whose 5% bound is proved)".

### M3. BJV's low offline power comes from its short windows, not from the lack of betting

* **Where.**
  * `docs/report_parts/E13.md:169`: "Its power is much lower in our designs ... For a lasting change after a known calm stretch, accumulating evidence (betting or CUSUM) wins."
  * `docs/report_parts/_header.md:27`: "BJV ... finds a 25% change in 6% to 30% of offline windows against 100%".
* **Evidence.**
  * **The check.** Same offline design as Q13a (D3): Gaussian, 2000 steps, change at step 1000. The overlapping statistic uses critical values at the simulated 95% quantile. Detection of a 25% increase / 20% decrease:
    * k = 50: 8.7% / 7.3%;
    * k = 125: 32.3% / 26.7%;
    * k = 250: 79.3% / 78.0% (truncated: 78.7% / 78.0%);
    * k = 500: 99.7% / 99.7%.
  * **Agreement with E13.** The k = 125 row agrees with Q13a's overlapping k = 125 row (30% and 33%) within Monte Carlo error (300 windows). Q13a has no overlapping k = 50 row; its block version at k = 50 gives 8% and 7%.
  * **What this shows.**
    * With k = 500 the test uses 500 returns on each side of the change, as much as the e-detector's calibration window, and it accumulates nothing. It matches the e-detector's 100% at matched size.
    * The gap in Q13a is therefore window length.
    * The e-detector's real advantages are elsewhere. It does not need to choose a window or know where the change is, and it runs online. A k = 500 test centred on a change in the middle of the window is optimally placed, which the e-detector is not.
  * `E13.md:169` partly says this ("its windows must be short"). The final sentence then credits accumulation, which the data do not separate from window length.
* **Fix.**
  * `E13.md:169`: replace "For a lasting change after a known calm stretch, accumulating evidence (betting or CUSUM) wins." with "The gap is mostly window length: with windows of 500 returns, as long as the e-detector's calibration window, the overlapping statistic finds 99.7% of 25% increases at matched size (Gaussian, plain version, which is within one point of the truncated one at k = 250; a check in review round 3, not in the CSVs). The e-detector does not have to choose the window, and it is online."
  * `_header.md:27`: "BJV with windows of 50 and 125 returns (a reconstruction ...) finds a 25% change in 6% to 30% of offline windows against 100%; longer windows close most of the gap".
  * Add k = 250 and k = 500 rows to the Q13a power table.

### M4. The real-data comparison is scored on the e-detector's time scale, so "less precise on crypto" is too general

* **Where.**
  * `docs/report_parts/E13.md:121`: "(this uses the real data, by design, and favours the competitors)".
  * `E13.md:170`: "clearly more often than the competitors'".
  * `README.md:28`; `README.md:30`: "the published methods of E13, tuned to the same alarm counts, are less precise on crypto".
  * `PITCH_BAPAT.md:35`; `PITCH_BAPAT.md:45`: "the published change-point methods of E13 are less precise on crypto at the same alarm counts".
  * `docs/report_parts/_header.md:27`.
* **Evidence.**
  * **What the outcome measures.** The outcome asks whether the robust volatility of the next 1440 minutes is above (or below) an EWMA forecast with a half-life of one day (`experiments/e13_benchmarks.py:465` to `494`).
    * The e-detector compares the clipped variance with a baseline calibrated on the previous 1440 minutes, so its alarms mark day-scale level shifts by construction.
    * Rolling BJV compares two 120-minute windows (`K_ONLINE = 120`, `e13_benchmarks.py:90`). BOCPD alarms on runs shorter than 250 minutes (`W_SHORT = 250`, `:89`). Their alarms flag changes over 2 to 4 hours, which can be real and still not last a day.
    * So the score measures persistence over the next day.
  * **The competitors' windows were not varied.** k = 120 and w = 250 were chosen on a Gaussian pilot (`E13.md:3`; `docs/iterations/E13.md`, iteration 2) and never changed on real data. What the comparison shows: at these window lengths, the competitors' alarms are less often followed by a day beyond a one-day forecast. It does not show that the methods are less precise on crypto.
  * **"Favours the competitors" is not supported.**
    * Matching uses the alarm counts only, not the outcome, so it gives the competitors no evident advantage.
    * The matching is per series and per direction (`e13_benchmarks.py:588`), and the forecast uses past data only. I found no leakage of the outcome into either.
  * **The tests.**
    * The two-proportion z-tests treat all alarms as independent.
    * Alarms of one method cluster in time, and the methods alarm on the same days, so the samples are dependent both within and between methods.
    * The caveat at `E13.md:178` says that alarms are not independent. But the reading says "clearly", and the tables print p < 0.001.
* **Fix.**
  * `README.md:30`: "the published methods of E13, at the window lengths we used (2 to 4 hours) and tuned to the same alarm counts, are less often followed by a day beyond the forecast on crypto". Same change at `PITCH_BAPAT.md:45`.
  * `E13.md:121`: replace "(this uses the real data, by design, and favours the competitors)" with "(this uses the real data, but not the outcome)".
  * `E13.md:170`:
    * replace "clearly more often" with "more often";
    * add: "The outcome asks whether the next day stays beyond a one-day forecast, the time scale of the e-detector's 1440-minute baseline; the competitors' windows (120 and 250 minutes) were not varied on real data, and the p-values treat clustered alarms as independent."
  * A check worth running: Q13d with BJV at k = 720 or 1440 and BOCPD with w = 1440. Alternatively, score every method also on a 2 to 4 hour horizon.

### M5. "Two to three times faster for decreases on clean data" holds in 2 of the 9 clean decrease cells

* **Where.**
  * `docs/report_parts/E13.md:167`: "On clean data (Gaussian, ticks, price noise) BOCPD detects *decreases* two to three times faster".
  * `README.md:28`: "two to three times for decreases".
  * `PITCH_BAPAT.md:35`: "On clean data BOCPD detects decreases two to three times faster".
  * `docs/report_parts/_header.md:27`.
* **Evidence.** Median delays, plug-in e-detector / matched BOCPD (`results/e13c_online_change.csv`):

  | Law | sd x 0.5 | sd x 0.67 | sd x 0.8 |
  |---|---|---|---|
  | Gaussian | 81 / 26 = 3.12 | 109 / 51 = 2.14 | 169 / 137 = 1.23 |
  | Ticks | 53 / 28 = 1.89 | 72 / 56 = 1.29 | 132 / 132 = 1.00 |
  | Price noise | 122 / 66 = 1.85 | 166.5 / 123 = 1.35 | 286 / 270 = 1.06 |

  * **The range fits only two cells.** Only the two largest Gaussian decreases fall in "two to three times". The other seven cells give 1.0 to 1.9.
  * **Detection shares differ.** BOCPD's medians are over fewer detected paths in the small cells: 86.5%, 85.5% and 64.5% at x 0.8, against 97.5% to 98.5% for the e-detector.
  * **Severity.** The sentence understates the e-detector. It is major because it appears as a headline in four files.
* **Fix.** In all four places: "On clean Gaussian data BOCPD detects large decreases two to three times faster (26 against 81 steps for a halving of the sd, 51 against 109 for sd x 0.67); for sd x 0.8, and with ticks or price noise, it is 1.0 to 1.9 times faster."

---

## 4. Minor findings

**theory_running_tau.md.**

* **m1. Section 4.3 (line 198): the k = 0 growth is not the Kullback-Leibler divergence.** The note says that "for the likelihood-ratio bet tuned to the true ρ this is the Kullback-Leibler divergence". The bet of E14 is the *clipped* likelihood ratio (clip at u₀ = P(|Z| > 3)). Its k = 0 growth is lower for increases (D6):
  * ρ = 2: 0.573 against 0.807;
  * ρ = 1.5: 0.194 against 0.220.

  Write: "for the unclipped likelihood-ratio bet tuned to the true ρ this would be the Kullback-Leibler divergence, the best possible rate; the clip lowers it for increases (0.573 against 0.807 at ρ = 2)".
* **m2. Remark B.3 (line 190) needs the condition given the past.** "Theorem B.2 holds for the new segment if its returns are exchangeable" should read "if its returns are exchangeable conditionally on the returns up to the alarm". The restart time is data dependent. The condition holds for i.i.d. returns and for any infinite exchangeable sequence. It is the form that the E14 caveat "exchangeable returns within a segment" (`E14.md:374`) needs.
* **m3. The status line 225 claims more than the tests check.** "Streaming equals vectorised (`tests/test_running.py`)" is tested for route A and the awake route A monitor (`test_streaming_route_a_matches_the_vectorised_version`, `test_awake_streaming_matches_vectorised`). For `ConformalMonitor` only one p-value is checked (`test_streaming_conformal_pvalue_matches_brute_force`). Add a test that `ConformalMonitor.run` and `conformal_detector` give the same alarms with the same sequence of θ.

**running.py.**

* **m4. An "off" decrease side grows by 1 per step.**
  * **The difference.** `src/avmonitor/running.py:336` and `:722` to `724` give an off step (m_dn ≤ 0 or τ² ≤ m_dn) the factor 1, so `S_dn` grows by 1 per step. `docs/math.md:223` sets the decrease statistic to 0 instead.
  * **Validity is unaffected.** E[M_t | F_{t-1}] ≤ M_{t-1} + 1 still holds.
  * **But a side that stays off alarms for certain.** If it stays off for 2A steps, it alarms deterministically. That is 72,000 steps for the awake route A monitor, far longer than any segment in E14, so no E14 number changes.
  * **Fix.** Use 0, as in `math.md`, or document the convention in the docstring.

**E14.md.**

* **m5. Line 192: "as the theory predicts".** The theory gives only the limit as k / n₀ → ∞ (g_k → a limit ≤ 0). That the growth is already near zero at k = n₀ = 1940 (-0.001 to 0.013, `results/e14w_dilution.csv`) is a numerical finding. Write: "The growth falls to about zero once the post-change values are as many as the pre-change values (k ≈ 1940); the theory gives only the limit ≤ 0 as k / n₀ → ∞."
* **m6. Line 116: "costs a little on small increases".** True for Gaussian returns (336 against 326 steps). Under Student t3 at x 1.25, the running level is 20% slower and detects less (937 against 783 steps, 76% against 88% detected; line 105). Write: "costs a little for Gaussian returns (336 against 326) and more under Student t3 tails (937 against 783 steps, 76% against 88% detected)".
* **m7. Line 347: compare like with like at the releases.** The paragraph compares the plug-in's 26 of 26 with the awake monitors' 30 and 35 of 48. At the releases where the plug-in was monitoring, the awake monitors flag fewer (from the table at lines 325 to 332 and the counts at line 345):
  * FOMC, plug-in awake at 26 statements: plug-in 26, awake-A 30 - 11 = 19, awake-B 35 - 12 = 23;
  * CPI, plug-in awake at 29 releases: plug-in 19, awake-A 25 - 15 = 10, awake-B 31 - 17 = 14.

  Add: "At the 26 FOMC statements at which the plug-in was monitoring, awake-A flagged 19 and awake-B 23 (plug-in 26); at the 29 such CPI releases, 10 and 14 (plug-in 19)."
* **m8. "Null band 0.25" means two different things.**
  * For the plug-in, the band 0.25 is a ±25% band on the baseline level.
  * For awake-B, it divides each rank payoff by 1.25 (`E14.md:319`). That tolerates a conditional mean payoff of up to 1.25 and has no direct volatility meaning.
  * Lines 286 and 313 treat it as one setting ("E12's frozen setting"). Say that the same number is used with a different meaning for the rank monitors.

**E13: competitors.py and E13.md.**

* **m9. The BOCPD cap is not checked by the test, and it does touch the short runs.**
  * **The claims.** `E13.md:13` says "a unit test checks the capped filter against the exact recursion". `competitors.py:313` to `314` says the merge "never [touches] the short runs the alarm rule looks at".
  * **The test.** `tests/test_competitors.py:125` runs 30 steps with r_max = 64, so the cap never binds. It checks the recursion, not the cap.
  * **My checks (D4).**
    * With r_max = 64 and w = 20, the capped short-run mass differs from the exact one by up to 0.41.
    * In the E13 setting with no change, the largest difference is 3.4 × 10⁻⁴.
    * With a change, the first crossing of the alarm level is at the same step in both filters. Later the short-run mass differs by up to 0.56 and the log odds by up to 5.6.
  * **Where it matters.** On real data the filter runs continuously through changes, so the approximation matters after the first alarm.
  * **Fix.**
    * `E13.md:13`: "We cap run lengths at 512 (an approximation: the merged bucket keeps the statistics of the heavier run; the unit test checks the recursion without the cap)".
    * Docstring: delete "and never the short runs the alarm rule looks at".
    * Add a test in which the cap binds.
* **m10. Line 3: "Every tuning constant was fixed before the reported runs."** True of the constants. But three methods were added after the first runs:
  * the studentised BJV and the CS-baseline e-detector after the first offline size run (`docs/iterations/E13.md:68`);
  * BOCPD on truncated returns after the first online runs (`:112`).

  The iteration log discloses this. Add to line 3: "Three methods marked as our modifications or as the CS baseline were added after the first runs (iterations 3 and 5); no constant was tuned on the reported runs."
* **m11. Line 121: the random-time shares cover only one side.** "At random times the share is 46.5% to 50.1% for each coin" is the share *above* the forecast (both half-years, `results/e13d_real_crypto_random.csv` and `..._h1_2025_random.csv`). The comparison for down alarms is the share *below*, 49.9% to 53.5%. Write: "At random times the next day is above the forecast 46.5% to 50.1% of the time and below it 49.9% to 53.5% (both half-years), and 50.3% and 49.7% for Nifty."
* **m12. Line 133: "p ≤ 0.006".** The largest of the eight p-values is 0.00615 (BOCPD, up alarms, band 0; recomputed from `results/e13d_real_pooled.csv`). Write "p ≤ 0.007". The held-out "p ≤ 0.007" (line 145, largest 0.00690) is right.
* **m13. Line 168: "81% to 90% of 25% increases against 98% to 99%".**
  * The 81% to 90% holds with Gaussian returns (BOCPD 90%, truncated 85%) and jumps (truncated 81%).
  * With ticks or price noise the two BOCPD versions detect 46.5% to 84.5% (`results/e13c_online_change.csv`).
  * The e-detector detects 97% to 99% over the five laws.

  Write: "(81% to 90% of 25% increases with Gaussian returns or jumps and 46.5% to 84.5% with ticks or price noise, against 97% to 99%)".
* **m14. Line 171: "roughly 1.5 to 2.5 times the delay".** The ratios of CS to plug-in median delays in Q13c run from 1.25 (Gaussian x 2: 25 / 20) to 2.57 (Student t5 x 1.25: 413 / 161). Write "1.25 to 2.6 times".
* **m15. Q13c "detected" counts first alarms in the wrong direction** (caption at line 80; `correct_direction` in the CSV). At sd x 1.25 and x 0.8, the share of detections whose first alarm is in the right direction is:
  * plug-in e-detector: 99.5% to 100%;
  * matched rolling BJV: 77% to 95%;
  * plain BOCPD: 69% to 100%.

  Line 169's "9% to 35% detected" for BJV therefore overstates correct detections a little. Add the correct-direction share to the tables, or add "(of which 77% to 95% in the right direction)".
* **m16. The nominal rolling BJV critical value uses the wrong length.** It is `C.mosum_crit(day, k, 0.05)` with day = 390 (`e13_benchmarks.py:355`). A day of 390 window positions spans 390 + 2k - 1 returns, which is the series length the MOSUM formula expects. The row is reported only as a side result and was conservative anyway (0.7% against 5%). Use 390 + 2k - 1, or say that the nominal level is approximate.

**README, pitch and report index.**

* **m17. "Ties" is not the right word.**
  * `README.md:28` says the truncated BOCPD "ties on down alarms (67.0%, p = 0.47)". 67.0% against 69.7% is not a tie.
  * `PITCH_BAPAT.md:35` ("with one tie on down alarms") and `_header.md:27` ("with one tie") mix this up with the only exact tie, which is per coin: DOGE down alarms, 41 of 58 each (`E13.md:145`).

  Write "where BOCPD on clipped returns is not significantly different on down alarms (67.0% against 69.7%, p = 0.47)", and drop "with one tie" in the pitch and the index or replace it the same way.
* **m18. The 25% is guaranteed only with a known baseline.** `README.md:26` and `PITCH_BAPAT.md:40` say "41% to 44% against a guaranteed 25%". The README also omits the condition that A.6 needs. Write: "41% to 44% over five days, where a known baseline would guarantee at most 25%". At `README.md:26` add "at a fixed calibration length" after "We proved that".
* **m19. `README.md:27` and `PITCH_BAPAT.md:41`: "56% to 66% for the default protocol".** The range covers three detectors under the standard recalibration protocol: plug-in, cs and B-lr (`results/e14h_streaming.csv`). The default monitor, the plug-in, catches 65.4% to 65.9%. Write: "against 56% to 66% for monitors that recalibrate for a day after each alarm (66% for the default plug-in)".

**Outside this scope: points from round 2b still in the pitch and the report index.** They are listed so that they are not lost; a fix may already be under way.

* `PITCH_BAPAT.md:39`: "at c = 3 even a total collapse of volatility takes at least 73 steps" (round 2b, C1) and "in all 12 simulated cells" (round 2b, C2).
* `docs/report_parts/_header.md:25`: "the delay attains log A / KL_inf" (round 2b, M2) and "(at least 73 steps at c = 3)" (round 2b, C1).

---

## 5. Derivations

### D1. Proposition A.6: p > 0 for every E14 law

* **The plug-in.**
  * `τ̂ = c × s × median |r_i|`, with c = 3 and s = 1.482602 (`running.scaled_mad`, which centres at zero).
  * `m̂` is the mean of `min(r_i², τ̂²)` over the same window W of n returns.
  * `B̂ = τ̂² / m̂`, and the bets are λ ∈ {0.02, ..., 0.5} (increase) and φ ∈ {0.1, ..., 0.95} (decrease).
* **The event.**
  * Let E = {ε/2 < |r_i| ≤ ε for all i ≤ n}, with P(E) = P(ε/2 < |r| ≤ ε)^n.
  * On E, the median of |r_i| lies in (ε/2, ε], so τ̂ lies in (2.224 ε, 4.448 ε].
  * Since τ̂ > ε, no calibration return is clipped, so m̂ = mean r_i² lies in (ε²/4, ε²].
  * Hence B̂ lies in (4.946, 79.13).
* **The drift of the bet λ = 0.5.** For a monitored return let Y = min(r², τ̂²) / m̂ ∈ [0, B̂].
  * Always, log(1 + 0.5 (Y - 1)) ≥ log 0.5 = -0.693.
  * On {|r| ≥ 4.45 ε} ⊂ {|r| ≥ τ̂}, Y = B̂ > 4.946, so the log factor is above log(1 + 0.5 × 3.946) = log 2.973 = 1.089.
  * With q_ε = P(|r| < 4.45 ε): g↑(W) ≥ 1.089 (1 - q_ε) - 0.693 q_ε. This is positive when q_ε < 0.611.
* **The largest factor.** Every increase factor is at most 1 + 0.5 (79.13 - 1) = 40.07, and every decrease factor at most 1 + 0.95 / (B̂ - 1) < 1.25. So f_max(W) ≤ 41 whichever bet attains g(W).
* **Conclusion.** p_{η,κ} ≥ P(E) > 0 with η = 1.089 (1 - q_ε) - 0.693 q_ε and κ = log 41, for every n.
* **The four E14 laws.**
  * **Gaussian, ε = 0.1.** q_ε = 0.344, so η = 0.477. P(0.05 < |r| ≤ 0.1) = 0.0398.
  * **Student t3 with unit variance, ε = 0.1.** q_ε = 0.503, so η = 0.193. P(0.05 < |r| ≤ 0.1) = 0.063.
  * **Jumps, ε = 0.1.** q_ε ≤ 0.995 × 0.344 + 0.005 < 0.35.
  * **Any law with no atom at 0 whose support contains 0.** As ε → 0, q_ε → P(r = 0) = 0, and some interval (ε/2, ε] with small ε has positive mass.
  * **Ticks of 0.5.** The lemma does not apply: there is an atom at 0, and ε = 0.5 gives q_ε = 0.976. Use the event that every calibration return is ±0.5, of probability 0.3493^n.
    * Then τ̂² = (1.5 s)² = 4.946, m̂ = 0.25 and B̂ = 19.78.
    * The exact drift of λ = 0.5 under the tick law, summed over the tick values, is 0.472.
    * The largest factor is 1 + 0.5 × 18.78 = 10.39.
* **What is left to numerics.** These events are very unlikely: 10^-546 for the Gaussian at n = 390, 10^-178 for the ticks. The proof only shows p > 0. The realistic size of the share, 56% to 83% of windows, remains the numerical result of Q14w.

### D2. The limit in Proposition A.6

* **The lower bound.** For fixed η and κ, the proof gives P(T_A ≤ m) ≥ p_{η,κ} (1 - (log(12A) + κ) / (η m)).
* **With m_A = (log A)².** The bracket tends to 1, so `liminf_A P(T_A ≤ m_A) ≥ p_{η,κ}`.
* **Why "→ p_{η,κ}" is wrong.** Since η and κ are arbitrary, any limit is at least the supremum of p_{η,κ} over η > 0 and κ < ∞. That supremum is P(g(W) > 0), which is in general larger than p_{η,κ} for a given η.
* **The conclusion survives.** A bound C m_A / A tends to 0, so the lower limit p_{η,κ} > 0 (D1) excludes it.

### D3. BJV offline power by window length

* **Design.** Gaussian, n = 2000, change at step 1000. The critical value is the 95% quantile over 2000 null windows (1000 for the truncated version). There are 300 change windows per cell, so the standard error is at most 2.9 points.
* **Results.** Overlapping statistic of `competitors.bjv_test`, rejection rates at sd x 1.25 / x 0.8:

  | k | plain | truncated at 4 bipower sd |
  |---|---|---|
  | 50 | 8.7% / 7.3% | |
  | 125 | 32.3% / 26.7% | |
  | 250 | 79.3% / 78.0% | 78.7% / 78.0% |
  | 500 | 99.7% / 99.7% | |

### D4. The BOCPD run-length cap

* **Setting.** Exact recursion without a cap, against `BOCPDVariance`. The predictive is the Student t with 2α degrees of freedom and scale √(β/α), the same log density as `competitors.py`.
* **Small cap (`check_bocpd_cap.py`).** r_max = 64, w = 20, H = 1/300, five paths of 400 steps with sd x 1.5 at step 250. The largest difference in short-run mass is 0.413.
* **E13 setting (`check_bocpd_cap2.py`).** r_max = 512, w = 250, H = 1/7800, α₀ = 2, warm start on 2000 calibration returns, 1500 monitored steps, one path per case:
  * **No change.** Short-mass difference at most 3.4 × 10⁻⁴; log odds 5.0 × 10⁻³.
  * **x 1.25 at step 500.** First crossing of the alarm level at step 650 in both filters. Later, short-mass difference up to 0.495 and log odds up to 5.60.
  * **x 0.8 at step 500.** First crossing at step 579 in both. Later, short-mass difference up to 0.560 and log odds up to 2.69.

### D5. The rank detector with ties and restarts

* **The p-value.** Let `size` earlier values be in the reference set, `le` of them ≤ a_t and `lt` of them < a_t. Then `p_t = (size - le + θ (le - lt + 1)) / (size + 1)`. This counts the strictly larger earlier values, plus θ times the ties including the current value. It is the smoothed conformal p-value of Proposition B.1, so it is uniform and independent under exchangeability, ties included. `conformal_detector` and `ConformalMonitor.pvalue` both implement this formula.
* **After a restart.** The first p-value is θ, which is uniform. Theorem B.2 then applies to the new segment (with the condition of m2).
* **Simulation (`check_restart.py`).** A = 100, restart with an empty reference set after each alarm, 30 streams of 20,000 steps. Alarms per stream:
  * Gaussian: 144.6;
  * ticks of 0.5: 146.6;
  * Student t3: 147.5.

  All are below the bound E[N_m] ≤ m / A = 200 of Remark B.3 and Corollary A.5 (δ = 0).

### D6. Growth of the clipped likelihood-ratio bet at k = 0

* **The divergence.** For a normal scale change by ρ, KL = (ρ² - 1 - log ρ²) / 2. This gives 0.0581, 0.2195 and 0.8069 for ρ = 1.25, 1.5 and 2, and 0.0431, 0.1275 and 0.3181 for ρ = 0.8, 0.667 and 0.5.
* **The clipped bet.** Its k = 0 growth in `results/e14w_dilution.csv` is 0.0553, 0.194 and 0.573 for the increases, and 0.043, 0.126 and 0.313 for the decreases.
* **Reading.** The clip at u₀ = P(|Z| > 3) costs little for decreases and up to 29% at ρ = 2.

---

## 6. Checked and found correct

### Mathematics (theory_running_tau.md)

* **Lemma 2.1 and (2.1).** The truncated mean as an integral of the tail function.
* **Theorem A.1.** A band for all levels gives a valid detector for any predictable τ_t.
* **Lemma A.2.** The closed forms of the integrated band, and the reduction from continuous laws to laws with atoms.
* **Proposition A.4a, Theorem A.4 and Corollary A.5.**
  * The snapshot band and chains stopped on 1{N > e(j)}, each a supermartingale.
  * The run-length part of Theorem A.4.
  * The restart count E[N_m] ≤ (δ + m / A) / (1 - δ), given Corollary 6.6 of `docs/math.md`.
* **Proposition A.6, the proof.**
  * The side statistic is at least W_t / K.
  * Optional stopping at T' ∧ k with bounded increments: factors at least 1/2 on the increase side and 0.05 on the decrease side.
  * The overshoot bound log f_max.
  * Markov's inequality on the event {g ≥ η, f_max ≤ e^κ}.
* **Proposition B.1 and Theorem B.2.** With ties and the randomised p-value. The betting functions have mean one, and the constant C_ρ of the clipped likelihood-ratio density is right.
* **Remark B.3.** For i.i.d. segments (m2 for the general case).
* **Section 4.3.** The dilution formula π_k, the limit g_∞ ≤ 0 by Jensen's inequality, and the freeze argument: E_R[(E[f(p) | R])^m] ≥ 1.

### Code

* **`competitors.py`, BOCPD.**
  * `BOCPDVariance` is Algorithm 1 of Adams and MacKay with the conjugate inverse gamma update and the exact Student t predictive.
  * The growth and change probabilities, the normalisation and the warm start are right.
  * `short_log_odds` is the log odds of a run shorter than w.
* **`competitors.py`, BJV.**
  * `gumbel_crit_max_abs_normal` uses a = √(2 log 2N) and b = a - (log log 2N + log 4π) / (2a), the Gumbel constants for the maximum of N absolute standard normals.
  * `mosum_crit` follows the cited formula.
  * The standardisation is 2σ⁴ for the Gaussian version and the sample variance of the squared returns for the studentised version.
  * The bipower scale and the truncation are as described.
* **`competitors.py`, rate matching.** `match_alarm_counts` matches by coordinate bisection, per series and per direction.
* **`e13_benchmarks.py`.**
  * Thresholds are matched per series, per band and per direction, with the pause equal to the calibration length (`:554`).
  * The `Validation` forecast uses past data only.
  * The pooled z-tests are as described (but see M4).
* **`running.py`.**
  * `two_sided_sr` matches the betting e-detector with fixed inputs.
  * `_LevelsCS` reproduces `eb_confseq`.
  * The running level is predictable: 3 × scaled MAD of the last 120 returns, every 10 steps, snapped to the grid.
  * The p-value formulas are as in D5.
  * `awake_route_a` and `AwakeRunningTauMonitor` agree with Theorem A.4.
* **`e14_running_tau.py`.** `plugin_drift` computes g↑, g↓, f_max, the run-length bound and the Markov bound as in Proposition A.6. My rerun and exact integration confirm one cell.

### Numbers that match their CSVs

All except those listed in sections 3 and 4.

* **E13.md.**
  * Q13a: both size tables (matched and nominal), the power table, the change-point accuracy, and plain BJV under jumps.
  * Q13b: the table and the two-day sentence.
  * Q13c: the three tables and the price-noise and tick sentence.
  * Q13d:
    * the three rate-matched tables;
    * "highest in all ten cells", and "nine of ten" with the DOGE tie at 41 of 58;
    * the band-0 sentences (except m12) and the alarm-rate paragraph;
    * the lag-one autocorrelations: -0.016 to 0.027, 0.013 and 0.132.
  * Q13e: the costs.
  * The reading, except M2, M5, m13 and m14. The offline nominal sizes, 11.4% to 18.6%, are right.
* **E14.md.**
  * Q14a:
    * the one-day and five-day tables;
    * route B at 1.7% to 3.2% and 14% to 18%; the oracle at 0.8% to 1.7% and 12% to 14%;
    * B-freeze at 6.4% to 7.0% and 31% to 34%; the plug-in at 40.6% to 42.9%; the band failure rates.
  * Q14b: the delay tables and the pre-change false alarm ranges.
  * Q14w: the band widths, the drift table (one cell rerun exactly), and the dilution table with "3% to 14%".
  * Q14h: the table and its reading:
    * blind 57% to 71%, awake at 31% to 49% of changes;
    * 0.029 and 0.035 per day; 99.6% and 98.0% detected, with medians 34 and 60.
  * Q14d: alarm rates, awake shares, "3.6 to 4.3 per day" and the confirmation shares.
  * Q14g: the release table and its counts: 114 against 55; 21 to 25 and 31; 26 to 30 and 35; 15 and 17; 11 and 12.
  * Q14e: the cost ratios.
* **README, pitch and index (E13 and E14 statements).** Every number agrees with the report parts; the points above concern only wording and conditions.

### Not checked

* **The cited papers.** I did not check them against their sources: the BJV paper is blocked here, and I took the Adams and MacKay recursion and Vovk's results from their standard statements.
* **Full reruns and data.** I did not rerun any experiment script in full, and I did not check the data hashes.
