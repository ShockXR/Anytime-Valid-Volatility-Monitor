# Independent review, round 2b

> *Note added 3 November 2026: the "independent reviewer" or "referee" below was a separate AI agent (Claude, by Anthropic) instructed to act as a sceptical referee. It was not a human expert, and the review has not been confirmed by a human.*

**Reviewer role.** An adversarial reader: a statistics journal referee and a code auditor. The goal is that every claim in the files below survives an expert reading.

**Scope.**

1. `docs/theory_delay.md` and `docs/report_parts/E11.md` (Q11a to Q11e and E11f).
2. The decrease-side truncation option `c_down`:
   * `src/avmonitor/monitor.py` and the `c_down` tests in `tests/test_monitor_cs.py`;
   * `experiments/e11f_split_c.py`, including whether its vectorised `first_alarms` matches the class;
   * `results/e11f_split_c_synthetic.csv` and `results/e11f_split_c_real.csv`.
3. Panel monitoring:
   * `docs/theory_panel.md` and `docs/report_parts/E15.md`;
   * `src/avmonitor/panel.py` and `experiments/e15_panel.py`;
   * e-BH on per-asset e-processes, sticky e-BH, FDR at stopping times, and the claim that the Shiryaev-Roberts statistic divided by time is not an e-value.
4. `src/avmonitor/design.py` and `docs/equation_map.md` Part 5.

**Read first.** `docs/review/independent_review_round2a.md`, which covered `docs/math.md` sections 5 and 6, E8 to E12, ALT3, the README and the strategy document. Those points are not repeated. Where a claim from this scope has spread into a file that round 2a covered (README, pitch, report index, iteration log), I list the places so that one fix can be applied everywhere, but I did not review those files again.

**State of the tree.** The review is against commit 527af52. While it was running, commit 4f9d870 corrected the share of zero returns in the E11f tick law from 38% to 20% (`E11.md` and the script comment). That text is now right (19.7%, D2). Finding C2 is about a different problem with the same law.

**What I ran.** Every number in `theory_delay.md`, `E11.md` (including E11f) and `E15.md` was recomputed from its CSV with pandas. Small checks, with scripts in `/tmp/claude-0/review2b/` and nothing written to `results/`:

* `check_speed_limit.py`: the decrease speed limit for the default grid and for single bets, c = 2 to 5, with `design.speed_limit` and in closed form.
* `check_klinf.py`: m₀, B, the best bet λ*, KL_inf and the Gaussian KL by numerical integration (scipy `quad`), written independently of the E11 code.
* `check_e11f.py`: the vectorised `first_alarms` of `e11f_split_c.py` against `VolatilityMonitor`, 216 runs.
* `check_e11f_ticks.py`, `check_e11f_ticks2.py`: the zero share of the tick law, and reruns of the tick cell with a fixed tick size (600 and 1000 paths).
* Closed forms recomputed by hand (section 5).
* `python -m pytest -q`: **240 passed in 58.5 s.** `tests/test_monitor_cs.py`, `tests/test_panel.py` and `tests/test_delay_theory.py` alone: 63 passed in 11.9 s.

No file in the repository was changed except this one.

**Severity.** *Critical*: wrong mathematics, a wrong number, or a conclusion the evidence does not support. *Major*: an overclaim, a missing condition or a misleading comparison. *Minor*: wording or clarity.

---

## 1. Summary verdict

**The proofs hold.** I re-derived the following and found no error in any proof:

* in `theory_delay.md`: Lemmas 0, 1, 2 and 5, Corollaries 3, 4 and 6, Theorems 1, 2 and 3, the time change behind the memory term ζ, and the false alarm identity of section 8;
* in `theory_panel.md`: Lemma 1, Propositions 1 and 2, Lemma 2, Theorem 1 and Remarks 3 to 6.

The code implements the stated formulas: `design.py`, the `c_down` option in `monitor.py`, `panel.py` and the metrics of `e15_panel.py`. The vectorised E11f code matches the class alarm for alarm.

**The problems are in four statements that have become headlines.**

* **C1.** "No valid detector can detect a total collapse of volatility in fewer than 73 steps" is false. 73 is the floor of the default five-bet decrease grid. A single bet, valid over the whole null, needs 58. The sentence is in the summary of `theory_delay.md`, in E11, in the report index, in the pitch and in the iteration log. `theory_delay.md:80` and the short paper state it correctly.
* **C2.** The E11f "ticks" law rounds returns *before* the volatility change, so the tick shrinks with volatility. With a fixed tick of half a standard deviation, `c_down = 2` is about 12% *slower* than c = 3 for a 20% fall. So "faster in every simulated cell, including rounded prices" is not supported.
* **M1.** `c_down` exists only with the plug-in baseline, which has no proof. "Valid by construction" and "the proofs are unchanged" need that condition.
* **M2.** Headlines say the delay "attains" log A / KL_inf, and they state the refined-grid upper bound as an equality. First-order optimality for bounded means is already known from a result the note itself cites.

The panel theory and code are correct as far as I checked, and every E15 number matches its CSV. The panel findings are minor wording points.

**Most important fixes before the work is shown.**

1. Replace the 73-step sentence everywhere (one line, C1).
2. Rerun the E11f tick cell with a fixed tick and soften the E11f conclusion (C2).
3. Scope "valid by construction" to a known baseline (M1).

---

## 2. Critical findings

### C1. "No valid detector can detect a total collapse in fewer than 73 steps" is false

* **Where.**
  * `docs/theory_delay.md:17`;
  * `docs/report_parts/E11.md:30`.

  The claim has spread to:
  * `docs/report_parts/_header.md:25` and `docs/experiment_report.md:27` ("at least 73 steps at c = 3");
  * `docs/experiment_report.md:923`, a copy of `E11.md:30`;
  * `PITCH_BAPAT.md:38` ("no valid detector can gather more than about 1 / c² of evidence per step against a fall, so at c = 3 even a total collapse of volatility takes at least 73 steps");
  * `ITERATION_LOG.md:303`.
* **Says.**
  * `theory_delay.md:17`: "With c = 3 this means at least 73 steps to detect even a total collapse of volatility."
  * `E11.md:30`: "with c = 3 no detector valid over the whole null can detect even a total collapse of volatility in fewer than 73 steps".
* **Evidence.**
  * The per-step cap log(B / (B - 1)) = 0.1172 nats at c = 3 (Corollary 4) is correct for every nonparametric e-variable. But 73 does not follow from the cap alone.
  * 73 is the floor of the **default decrease grid** (five bets at 0.1, 0.25, 0.5, 0.75 and 0.95 of λ_max, weight 1/5 each) at A = 7800, as `theory_delay.md:80` says correctly. The grid spends most of its weight on bets that earn less than the largest factor per step.
  * A Shiryaev-Roberts e-detector with **one** bet at 0.999 λ_max is valid over the whole null: its average run length is at least A under every null law (`docs/math.md`, section 5). When every return is zero, it reaches A = 7800 after **58** steps. This is deterministic (D1).

    `design.speed_limit` gives the floors below, at the one-sided level A = 7800, with the two-sided level 2A in brackets:

    | c | default grid | single bet at λ_max | log A / log(B / (B - 1)) |
    |---|---|---|---|
    | 2 | 36 (39) | 29 (32) | 34.3 |
    | 3 | 73 (79) | 58 (64) | 76.5 |
    | 4 | 122 (134) | 96 (107); 97 at 0.999 λ_max | 138.9 |
    | 5 | 180 (199) | 141 (158) | 219.5 |

  * For detectors that are not e-detectors, the only lower bound proved here is asymptotic. Theorem 2b gives (log A / KL_inf)(1 + o(1)), which at c = 3 is 76.5 × (1 + o(1)). That is not a floor at a finite A: the single-bet detector is already below 76.5.
* **Fix.**
  * `theory_delay.md:17` and `E11.md:30`: "With c = 3 the default decrease grid needs at least 73 steps (79 at the two-sided level 2A), even when volatility collapses to zero. No Shiryaev-Roberts e-detector whose bets are valid over the whole null can do it in fewer than 58 (64)."
  * `_header.md:25`, `experiment_report.md:27`: "(at least 73 steps at c = 3 for the default bets)".
  * Pitch: "no valid e-detector can gather more than about 1 / c² of evidence per step against a fall; with the default bets at c = 3 even a total collapse of volatility takes at least 73 steps".
  * Apply the same change to `ITERATION_LOG.md:303` and `experiment_report.md:923`.

### C2. E11f: in the tick law the tick shrinks with volatility; with a fixed tick, c_down = 2 is slower for a 20% fall

* **Where.**
  * `experiments/e11f_split_c.py:55` (the law) and `:117` (the change);
  * `docs/report_parts/E11.md:78` (method), `:89` (reading) and `:103` (conclusion);
  * `docs/theory_delay.md:214`.

  The claim has spread to:
  * `PITCH_BAPAT.md:38` and `ITERATION_LOG.md:303`;
  * `README.md:40` ("a tighter c_down = 2 for the decrease side, which detects falls faster");
  * `docs/experiment_report.md:982` and `:996`.
* **Says.**
  * `E11.md:89`: "With c_down = 2 the decrease delay falls in all 12 cells, by 4% to 51%".
  * `E11.md:103`: "it is faster in every simulated cell".
  * `theory_delay.md:214`: "in all 12 simulated cells, including fat tails, jumps and rounded prices".
* **Evidence.**
  * **The law.** The script draws rounded returns, `np.round(z / 0.5) * 0.5`, for the whole path. It then multiplies the part after the change by the ratio (`base[:, CAL:] *= ratio`). So after the change the returns are multiples of 0.4, 0.335 or 0.25, not of 0.5. The tick shrinks in proportion to volatility, and the share of zero returns stays at 19.7% before and after the change.
  * **Real prices have a fixed tick.** When volatility falls, more returns become zero. Rounding *after* the change of scale, with the tick fixed at 0.5, gives a zero share after the change of 24.5%, 29.1% and 38.3% for sd × 0.8, 0.67 and 0.5 (D2).
  * **Rerun with a fixed tick.** The calibration, band, A and `c_down` values are the script's. The code is my own, not the committed script.

    | Median decrease delay, c_down 3 / 2 / 1.5 | sd × 0.8 | sd × 0.67 | sd × 0.5 |
    |---|---|---|---|
    | the script's law (CSV, 600 paths) | 186 / 178 / 254 | 96 / 69 / 80 | 68 / 38 / 33 |
    | fixed tick of 0.5 sd (600 paths) | **192 / 219 / 256** | 99 / 83 / 87 | 69 / 40 / 35 |

    A second run with 1000 paths (seed 2026) at sd × 0.8 gave median / mean delays of 192 / 209, 215 / 233 and 254 / 275 for c_down 3, 2 and 1.5. With a fixed tick of half a standard deviation, c_down = 2 is about 12% **slower** than c = 3 for a 20% fall, not 4% faster.
  * **The answer depends on where the clip falls on the tick grid.** At sd × 0.8 the median / mean delays for c_down 3, 2 and 1.5 are:
    * fixed tick of 0.25 sd: 293 / 303, 185 / 202 and 186 / 205;
    * fixed tick of 0.75 sd: 301 / 311, 204 / 220 and 210 / 229.

    So c_down = 2 is faster at these two tick sizes and slower at 0.5 sd. "Faster in every cell, including rounded prices" holds only for a law that does not model a fixed tick.
  * **The "ticks" lesson is model-specific.** "c_down = 1.5 should not be used on data with many ties" (`E11.md:103`) also rests on that law. With a fixed tick, c_down = 1.5 is the slowest of the three at sd × 0.8 for a tick of 0.5 sd, but close to c_down = 2 at ticks of 0.25 and 0.75 sd.
* **Fix.**
  * Round after scaling: pass the ratio to `draw`, or round in the loop after `*= ratio`.
  * Rerun the tick cells at two or three tick sizes and report what they show.
  * Until then, replace "in all 12 cells" and "faster in every simulated cell" with: "faster for the three continuous laws; with rounded prices the result depends on the tick size and can reverse for small falls".
  * Make the same change in `theory_delay.md:214`, the pitch, the iteration log and `README.md:40`.

---

## 3. Major findings

### M1. c_down: "valid by construction" needs a known baseline, which the option does not allow

* **Where.**
  * `docs/report_parts/E11.md:89` ("The proofs are unchanged: each side is valid for its own clipped statistic") and `:103` ("valid by construction");
  * `docs/theory_delay.md:214`;
  * `src/avmonitor/monitor.py:99` to `100` (comment: "the two-sided guarantee is unchanged") and `:128` to `129` (the cs baseline is rejected).
* **Evidence.**
  * **The argument is correct with a known baseline.** Each side is an e-detector for its own truncated statistic, and the two-sided rule needs only each side's guarantee. This holds when each side's m₀ is the true null mean, or the confidence bound of `docs/math.md` section 6.
  * **The option rules that out.** `VolatilityMonitor(baseline="cs", c_down=...)` raises `ValueError`, so `c_down` exists only with the plug-in baseline, which has no guarantee (E4e; round 2a).
  * **The measured rates exceed the bound.** In E11f all 12 no-change cells have one-day false alarm rates of 5.4% to 6.7%, above the 5% bound. This includes `c_down = c`, so the excess comes from the baseline, not from the split.
  * **Wording of the bound.** `E11.md:89` says "above the 5% that the true baseline would guarantee". That reads as if the true baseline gives 5%; it guarantees *at most* 5%.
  * **The decrease-side null changes.** It now reads "the mean of min(r², τ_dn²) has not fallen". `theory_delay.md:214` says correctly that "volatility has not fallen" implies this for a fixed innovation shape. But the smaller clip is also more sensitive to changes of shape at constant volatility. A larger share of zero returns is one example (C2).
* **Fix.**
  * **Either implement c_down for the cs baseline.** The existing cs baseline already spends δ/2 on each side. Feeding the lower sequence min(r², τ_dn²) instead of min(r², τ²) keeps the total at δ and the proof of `math.md` section 6 for each side.
  * **Or reword.** Write: "valid when the baseline is known; with the plug-in baseline the false alarm rate was 5.4% to 6.7% per day, as for c_down = c". Change "the 5% that the true baseline would guarantee" to "the 5% bound that holds with the true baseline". Add one sentence on the narrower null of the decrease side.

### M2. Optimality headlines state an upper bound as an equality and do not credit a known result

* **Where.**
  * `docs/report_parts/_header.md:25` and `docs/experiment_report.md:27`: "Proved delay bound and asymptotic optimality (the delay attains log A / KL_inf)".
  * `docs/theory_delay.md:20`: "its delay is (log A + ½ log log A + O(1)) / KL_inf (refined grid, proved)".
  * `docs/theory_delay.md:448`, the status table.
  * `docs/report_parts/E11.md:59`: "a grid refined with about √(log A) bets has delay ...".
* **Evidence.**
  * **The default grid does not attain log A / KL_inf.** By Theorem 3a its limit ratio is KL_inf / g_grid = 1.00 to 1.08 in the E11 cases (`theory_delay.md:20`). At A ≤ 10^8 the simulated delays are 5% to 24% *below* log A / KL_inf, because of the Shiryaev-Roberts memory term (section 9.3). "Attains" is true, to first order, only for a grid that contains λ* or is refined as A grows.
  * **Theorem 3b proves an upper bound only.** It gives WADD ≤ (log A + ½ log log A + O(1)) / KL_inf. The lower bound proved (Theorem 2b) is first order only, (log A / KL_inf)(1 + o(1)). "Is" and "has delay" present the second-order term as exact.
  * **Known result.** `theory_delay.md:332` says that Ram and Ramdas (2026) "prove this lower bound in general for composite pre- and post-change classes and show that it is attained for bounded means". The truncated statistic has a bounded mean. So first-order optimality of a suitable mixture is already known, and Theorems 2b and 3 are instances of it with explicit constants for this detector. The status table lists them as PROVED with no such note. I did not check whether Ram and Ramdas also give the √(log A) refinement rate.
* **Fix.**
  * `_header.md:25` and `experiment_report.md:27`: "Proved delay bound; a refined grid is first-order asymptotically optimal (its delay divided by log A / KL_inf tends to 1); the default grid's limit is 1.00 to 1.08 times that".
  * `theory_delay.md:20`, `:448` and `E11.md:59`: "has delay at most".
  * After Theorem 3: "First-order optimality for bounded means is known (Ram and Ramdas 2026). This section gives the refinement rate and the constants for our grid; we make no claim of priority for the first-order result."

---

## 4. Minor findings

**theory_delay.md.**

* **m1. The mixture sandwich at line 50 is wrong as written.** It says log M_n lies between max_k (log w_k + log L^{(1,k)}_n) and that plus log(Kn).
  * The upper half fails: with n = 2, one bet, f₁ = 0.1 and f₂ = 2, M₂ = 2.2 but the bound gives 0.4 (D9). A Shiryaev-Roberts sum also contains products that start later.
  * The correct form is log M_n ≤ log n + max_{j,k} log L^{(j,k)}_n.
  * The proof of Theorem 3a is correct, so only this sentence needs the change.
* **m2. "Computed exactly by dynamic programming."** Line 19 says this, and `E11.md:47` says "exact worst case". The dynamic programme tracks a rounded-down state, so it gives a valid *lower* bound on the worst case, as `E11.md:70` itself says. Write "a lower bound on the worst case, computed by dynamic programming". The same applies to the Section 8 row of `docs/equation_map.md`.
* **m3. "At m = A, exactly"** (line 306; `E11.md:55`). The argument (Y ≡ 1 gives M_t = t) is right for a one-sided detector with threshold A. The two-sided rule alarms at 2A per side, so under Y ≡ 1 it alarms at step 2A. Say "for each side, or for a one-sided detector".
* **m4. The log-scale shares in section 6.4 (line 218) and `E11.md:30` were computed from rounded ratios.** From the CSV values (D7), for sd ratios 1.25, 1.5 and 2:
  * truncation 0.192, 0.250 and 0.352;
  * nonparametric 0.648, 0.676 and 0.588;
  * grid 0.160, 0.074 and 0.060.

  So the ranges are 19% to 35%, 59% to 68% and 6% to 16%, not 20% to 36%, 59% to 66% and 6% to 13%. Also say that this splits the first-order *prediction* KL / g_grid, not the simulated penalty.
* **m5. Section 6.2 (line 194): "within 0.03 in all six cases".** The differences between the simulated oracle penalty and the second-order prediction are 0.015, 0.012, 0.035, 0.008, 0.001 and 0.008, so the largest is 0.035, at sd 2. The simulated oracle penalty at sd 2 is 18.470 / 12.449 = 1.484, printed as 1.49; it should be 1.48 (D10). Write "within 0.04".
* **m6. Rounding** (`results/e11*.csv`):

  | Where | CSV value | Printed | Should be |
  |---|---|---|---|
  | line 140, jumps sd 1.25, Wald bound | 247.45 | 248 | 247 |
  | line 145, Gaussian sd 0.8 (down), grid at ν = 1 | 194.48 | 195 | 194 |
  | line 204, c = 4, decrease delay at sd 0.67 | 203.49 | 204 | 203 |
  | line 173, nonparametric loss at sd 0.8 | 1.2447 | 1.25 | 1.24 |
  | line 46, E_Q[Y] at sd 2 | 3.1295 | 3.130 | 3.129 |

* **m7. Line 48: "Over all nine increase cases".** The ranges quoted, "1.4 to 4.4" and "49% to 94%", include Gaussian sd 1.1, which makes ten cases. The nine cases alone give 1.69 to 4.44 and 49% to 86%. Either say ten, or quote the nine-case ranges.
* **m8. Line 18: "within 26% of the best single bet everywhere".** This holds inside the design range (largest ratio 1.258). At sd 3 the default grid reaches 1.36 (1.362 with jumps, 1.355 Gaussian). Write "everywhere inside the design range".
* **m9. Line 3: "Every number comes from `experiments/e11_delay_theory.py`".** Section 6.3 also cites E11f (`e11f_split_c.py`). Name both scripts.

**E11.md.**

* **m10. The Q11d table (lines 49 to 53) mixes two runs.** The row "two-sided default (E9), 40,000 paths" takes its false alarm rate and rate / bound from the 40,000-path run, but F2 and F3 from the 10,000-path run. Use one run, or label the columns.
* **m11. Line 5 quotes E4 medians (87, 35 and 17); the theory tables use E4 means (99.7, 39.6 and 18.0).** Label them "median delays" so that a reader comparing with section 5 is not confused.
* **m12. Line 89: "c_down = 1.5 is ... slower than c_down = 2 for a 20% fall".** False for Student t4, where c_down = 1.5 is faster (277 against 288.5). Write "slower for three of the four laws".
* **m13. Line 89 compares E11f medians with the theory's means.** "29%, 46% and 51% for Gaussian data, against the theory's 26% to 49%": the 29%, 46% and 51% are reductions in median delay, while the theory's 26% to 49% are in mean delay. The E11f mean delays give 25%, 44% and 50%. Compare like with like.
* **m14. Line 78: "The vectorised code is checked alarm for alarm against the class".** `check_against_class` covers one law (jumps), c_down 3 and 1.5, 4 paths and the first alarm only. My wider check found no mismatch:
  * runs: 216 (four laws, ratios 1, 0.8 and 0.5, c_down 3, 2 and 1.5, 6 paths each);
  * mismatches: 0.

  The statement is therefore true, but the committed check is narrower than the sentence. Move a check of this size into the test suite.
* **m15. The real-data percentages use a smaller denominator than the printed count.** They are shares among alarms that have a following day. For example, BTC at c = 3 shows "(55)", but the 85% is 46 of 54. State the denominator.

**c_down code and tests.**

* **m16. The c_down equality test is vacuous.** `tests/test_monitor_cs.py:400` to `410` compares the alarms of the default and of `c_down=3.0` on a no-change path. Both runs give 0 alarms, so the test would pass even if `c_down` were ignored. Use a path with alarms on both sides, or compare the internal bets.
* **m17. Tests to add:**
  * `c_down` with `margin_z` (the decrease-side standard error `se_dn`);
  * fallback and degenerate calibration windows, and missing (NaN) returns with `c_down ≠ c`;
  * that `lam_dn = down_fractions / (tau2_dn / m0_dn - 1)`, which I verified by hand;
  * the `segments` record.
* **m18. Document c_down in the code.** The module docstring of `monitor.py` does not mention `c_down`, and `segments` records `tau2_dn` but not `m0_dn`. Add both, so that a run can be audited from its record.

**Panel.**

* **m19. The window identity (`theory_panel.md:90`).** It says "the naive rule `E ≥ 1/α` with h = one block of m steps is the single-asset monitor with `A = m / α`".
  * Under the window identity, E ≥ 1/α means (M↑ + M↓)/2 ≥ m/α - (m - n).
  * The monitor alarms when max(M↑, M↓) ≥ 2A. That implies the mean is at least A, but not the other way round, and the panel threshold is lower by (m - n).

  Write "is close to the single-asset monitor".
* **m20. The pathwise ordering (`E15.md:16`).** "Path by path, naive flags an asset no later than e-BH, and e-BH no later than Bonferroni" is what `test_rule_ordering_is_pathwise` checks. That test uses the three panel rules on the same e-values (`run_block` with `rule="naive"`, `"ebh"` and `"bonferroni"`). The methods list just above names "naive SR (A)" and "Bonferroni SR (A × K)", which are the per-asset detectors. Say that the ordering is for the panel rules on the same e-values.
* **m21. `E15.md:136`: the two-proportion tests compare overlapping sets.** The discovery sets of e-BH and of the rate-matched threshold overlap, so the two samples are not independent and the p-values are not exact. The hit rates differ by at most 1.4 percentage points, so the conclusion stands. Add "(the samples overlap, so the tests are only indicative)".
* **m22. `E15.md:143`: "The naive rule's extra discoveries are similar".** The confirmation rates of up discoveries made only by the naive rule range from 36% to 77% across cells (for example, July to December 2025 at band 0.5: 10 of 13). That is wider than the 42% to 57% of e-BH's marginal discoveries. Give the range.
* **m23. Which result gives the cost of late changes.** `E15.md:110` speaks of "the log(1/w) cost of Proposition 2", and `theory_panel.md:144` says "the delay of Proposition 2, growing like 2 log(ν/h)". Proposition 2 is a lower bound for *any* e-process. The cost log(1/w_ν) of the Pareto weights, about log h + 2 log(1 + ν/h), belongs to the construction. Write "the cost of the Pareto weights, which by Proposition 2 cannot be avoided in order".

**equation_map.md Part 5.**

* **m24. Part 5 (lines 64 to 75) is incomplete.**
  * It lists `design.approx_best_bet` under Lemma 2 (best nonparametric growth). That function is the second-order aGRAPA approximation discussed in section 2.1, not the oracle bet. Give it its own row.
  * It has no rows for Lemma 5, Theorems 2 and 3, Corollary 6, or the `c_down` option.
  * There is no part for `docs/theory_panel.md` and `src/avmonitor/panel.py`. Add one: e-BH, `change_weights`, `PanelEProcess`, sticky e-BH (`run_block`), and the Lemma 2 and Proposition 1 tests.

---

## 5. Derivations

### D1. The decrease speed limit, single bet

* **Setting.** At c = 3 for Gaussian returns, B = τ² / m₀ = 9.04516 (m₀ = 0.995007). On the decrease side U = 1 - Y, λ_max = 1 / (B - 1), and with every return zero Y = 0 and U = 1.
* **Factor per step.** A single bet at λ_max earns F = 1 + 1/(B - 1) = B / (B - 1) = 1.12430 on every step, so log F = 0.11716. For comparison, 1 / (B - ½) = 0.11703.
* **Shiryaev-Roberts sum.** M_n = Σ_{j=1}^{n} F^{n-j+1} = F (F^n - 1) / (F - 1). Since (F - 1)/F = 1/B, M_n ≥ A is the same as F^n ≥ 1 + A/B.
* **At A = 7800.** n ≥ log(1 + 7800/9.04516) / 0.11716 = log(863.3) / 0.11716 = 6.761 / 0.11716 = 57.7, so 58 steps.
* **At 2A = 15,600.** log(1725.7) / 0.11716 = 63.6, so 64 steps.
* **At 0.999 λ_max.** The factor is 1.124176, which also gives 58.
* **Validity.** Any λ in [0, λ_max] gives a factor with null mean at most 1 under every null law, so the single-bet Shiryaev-Roberts detector has average run length at least A over the whole null.
* **The grid.** It spends weight 1/5 on each of five bets, the largest earning 1 + 0.95/(B - 1). That is why its floor is 73.
* **Without the sum over start times.** A single product would need log A / log F = 76.5 steps. This is also the first-order asymptotic bound of Theorem 2b for a total collapse.

### D2. The zero share of the tick law

* **Before the change.** round(Z / 0.5) = 0 exactly when |Z| < 0.25, so the zero share is 2Φ(0.25) - 1 = 0.1974. The simulation gives 0.1974.
* **After the change, with the script's law.** The script multiplies the rounded values by s, so a zero stays a zero and the share stays at 0.1974.
* **After the change, with a fixed tick.** Rounding sZ to the tick 0.5 gives zero when |Z| < 0.25/s. The share is 2Φ(0.25/s) - 1 = 0.2453, 0.2910 and 0.3829 for s = 0.8, 0.67 and 0.5.
* **Where 38% came from.** The earlier "38%" was 2Φ(0.5) - 1, the share of |Z| < 0.5.

### D3. The memory term ζ (section 4.3)

* **Time change.** In the Brownian approximation, log L of the bet is a Brownian motion with drift g and variance s² per step. Then the integral over lags t of exp(-log L_t) has the law of V = 2 / (s² Γ_a), with Γ_a a gamma variable of shape a = 2g / s² (Dufresne's identity).
* **The memory term.** ζ = E log V = log(2 / s²) - ψ(a), with ψ the digamma function.
* **Small changes.** When s² ≈ 2g, a ≈ 1 and ψ(1) = -γ, so ζ ≈ log(1/g) + γ. This matches the note.

### D4. Proposition 1 of the panel note

* **The two values.** M₁ = f₁ and M₂ = f₂ (f₁ + 1).
* **The stopping time.** σ = 2 only on {f₁ < 1}, which is known at step 1.
* **Taking expectations.** With E[f₂ | F₁] = 1:

  E[M_σ / σ] = E[f₁ 1{f₁ ≥ 1}] + E[1{f₁ < 1}(f₁ + 1)] / 2
  = E[f₁] - E[f₁ 1{f₁ < 1}] + E[1{f₁ < 1}(f₁ + 1)] / 2
  = 1 + E[1{f₁ < 1}(1 - f₁)] / 2
  = 1 + E[(1 - f₁)⁺] / 2.

* **Conclusion.** This exceeds 1 whenever P(f₁ < 1) > 0, so M_t / t is not an e-value at a stopping time even though E[M_t / t] ≤ 1 at every fixed t. Correct.

### D5. Theorem 1 (3) of the panel note and α'

* **Where the bound comes from.** By Lemma 1, FDP_t ≤ (α/K) Σ_{null} e_{k,t}. The average of the null e-processes, stopped at their change times, is a test supermartingale (Lemma 2).
* **Ville's inequality.** It gives P(sup_t FDP_t ≥ x) ≤ min(1, α K₀ / (K x)) ≤ min(1, α / x).
* **Integrating.** E[sup_t FDP_t] ≤ ∫₀¹ min(1, α/x) dx = α + α log(1/α).
* **Numbers.** At α = 0.05 this is 0.1998. The level α' with α'(1 + log(1/α')) = 0.05 is α' = 0.0087, since 0.0087 × (1 + log(1/0.0087)) = 0.04998. Both match Remark 6.

### D6. Lemma 2b and the KL_inf values

* **Mass 1.** At an interior optimum the first-order condition is E_Q[U / (1 + λ* U)] = 0. P* = Q / (1 + λ* U) then has mass E_Q[1 / (1 + λ* U)] = 1 - λ* E_Q[U / (1 + λ* U)] = 1.
* **Mean on the boundary.** Its mean of U is E_Q[U / (1 + λ* U)] = 0, so P* is on the boundary of the null.
* **The divergence.** KL(Q ‖ P*) = E_Q log(1 + λ* U) = g(λ*).
* **Numerical check.** Gaussian, c = 3, by `scipy.integrate.quad` (my code); every value matches `results/e11b_information.csv` and `results/e11a_delay.csv`:

  | case | λ* | KL_inf | Gaussian KL |
  |---|---|---|---|
  | sd 1.25 | 0.241213 | 0.051381 | 0.058106 |
  | sd 2 | 0.643042 | 0.430230 | 0.806853 |
  | sd 0.8 (down) | 0.123261 (0.99165 λ_max) | 0.034650 | 0.043144 |
  | sd 0.5 (down) | at λ_max | 0.088097 | 0.318147 |

### D7. Exact log-scale shares of the price of robustness

For each sd ratio, the share of a factor is the log of that factor divided by log(KL / g_grid). From the unrounded CSV values:

| sd | truncation | nonparametric | grid |
|---|---|---|---|
| 1.25 | 0.192 | 0.648 | 0.160 |
| 1.5 | 0.250 | 0.676 | 0.074 |
| 2 | 0.352 | 0.588 | 0.060 |

The rounded ratios printed in the table (for example 1.03, 1.10 and 1.02 at sd 1.25) give visibly different shares, because log(1.02) and log(1.03) are small numbers known to only one significant figure.

### D8. Theorem 3b constants

* **Grid size and spacing.** With K_A = 1 + ⌈(b - a)√(log A)⌉ equally spaced bets on [a, b], the bet nearest λ* is within 1/(2√(log A)).
* **Cost of the gap.** By Corollary 6, KL_inf - g(λ_k) ≤ C / log A.
* **Cost of the weight.** It is log K_A = ½ log log A + O(1).
* **Theorem 1 bound.** WADD ≤ (log A + log K_A + O(1)) / (KL_inf - C / log A).
* **Expanding.** This equals (log A + ½ log log A + O(1)) / KL_inf, since log A × (C / log A) / KL_inf² = O(1).
* **Conclusion.** This is an upper bound (M2).

### D9. Counterexample to the sandwich at line 50

* **Setting.** One bet (K = 1, w = 1), n = 2, f₁ = 0.1, f₂ = 2.
* **The statistic.** L^{(1)}_2 = 0.2 and L^{(2)}_2 = 2, so M₂ = f₂ (f₁ + 1) = 2.2.
* **The claimed upper bound.** log(0.2) + log(2) = log(0.4).
* **The correct bound.** M₂ ≤ n max_j L^{(j)}_2 = 4.

### D10. Section 6.2, oracle penalties

| | sd 1.25 | sd 1.5 | sd 2 | sd 0.8 | sd 0.67 | sd 0.5 |
|---|---|---|---|---|---|---|
| \|simulated - predicted\| | 0.015 | 0.012 | 0.035 | 0.008 | 0.001 | 0.008 |

The simulated penalty is the robust oracle-bet delay divided by the Gaussian LR oracle delay. At sd 2 it is 18.470 / 12.449 = 1.484.

---

## 6. Checked and found correct

### Mathematics (theory_delay.md)

* **Lemmas 0 and 1.** Shape of g and its unique maximiser.
* **Lemma 2.**
  * (a): the supergradient argument, and the boundedness of the e-variable for u > 0.
  * (b): P* has mass 1 and mean 0, and KL = g(λ*) (D6).
* **Corollaries 3 and 4.** log(B / (B - 1)) = 0.11716 ≈ 1 / (B - ½) = 0.11703.
* **Theorem 1.** Pathwise domination by the best single bet, Wald's identity, and Lorden's overshoot bound.
* **Theorem 2.** The proof of 2b, which applies Lorden (1971) at the least favourable law P*.
* **Theorem 3.** Parts (a), (b) (the bounds on K_A and the expansion of η_A, D8) and (c).
* **Lemma 5 and Corollary 6.** Lemma 5: the change at step 1 is the worst case. Corollary 6 includes the grid factors 0.89 and 0.97.
* **Section 4.3.** The time change behind ζ (D3).
* **"The first-order formula overstates by 60%".** 8.962 / 5.591 = 1.603.
* **Section 8.**
  * The false alarm identity and its factors F1, F2 and F3.
  * The uneven split (log 1.5 = 0.41).
  * The time-varying threshold (5%, 0.05 nats).
  * The validity of the dynamic-programming lower bound: the rounded-down tracked state is dominated by the actual process.

### Numbers in theory_delay.md and E11.md that match their CSVs

All except those listed in sections 2 to 4.

* **Section 5 and its ranges.**
  * Second-order prediction 0.972 to 1.018; first order 0.938 to 2.119; Wald bound 1.294 to 2.682.
  * Grid / oracle 1.093 to 1.261; second-order grid 1.029 to 1.320.
  * ν = 501: 0.666 to 0.831; E4 / ν = 501: 0.955 to 1.105.
  * ζ agreement 0.080 for g < 0.2 and 0.209 otherwise.
* **Section 6.**
  * 6.1 and the 6.2 table (apart from m5).
  * 6.3 and its optima: c = 3.25, 3.5 and 4.25 with jumps; 1.25, 1.5, 2 and 3 for decreases.
  * Efficiencies 0.84, 0.94, 0.88 and 0.80; "at least 57%"; c = 3 costs 3%, 6% and 12%; 26% to 49%.
* **Section 7.**
  * Design 0.2565 / 0.4893, 1.18 and 1.37; pooled design 1.27 / 1.26; the validation table.
  * 0.895 to 0.996 (median 0.948); 0.82 to 1.55 (median 0.91).
  * False alarms 1.0% to 3.2%, 1.6% to 4.6% and 2.9% to 5.5%; "three of five"; 1.10 to 1.17 and 1.02 to 1.05.
* **Section 8.** The tables, the 60% / 31% shares, 0.79 / 0.24 / 0.00, and the two-point laws 2.12%, 1.98% and 3.79%.
* **Section 9.3.**
  * The table; 5% to 24% and 10% to 47%; 0.024 and 0.037.
  * ζ - κ between 1.0 and 4.5; excess nats 0.12 to 0.95 and 0.18 to 1.13.
  * Fine mixture 1.04 to 1.13 and 1.07 to 1.19; the KL_inf / KL table; t3 KL 0.026.
* **E10 cross-check.** 3218 / 1285; 0.41, 0.46, 0.51 and 0.58; 54% to 62%.
* **E11f.**
  * The synthetic and real tables; 4% to 51%; "at most 0.7 points".
  * "Four of five coins"; the interval 73.8% to 92.4% for 85% of 55 (printed as "about 74% to 92%").
  * Calibration time 64% to 70%.

### Code

* **`design.py`.**
  * The formulas: `growth_rate`, `best_bet`, `lam_max`, `zeta_brownian`, `wald_delay_bound`, `overshoot_bound`, `speed_limit` and `sr_delay_approx`.
  * `worst_case_false_alarm` returns a probability reached by an admissible adversary, hence a valid lower bound on the worst case.
* **`monitor.py`, c_down.**
  * Validation of `c_down`, and the cs baseline rejected.
  * `_calibrate`: `tau2_dn`, `m0_dn` and `se_dn` for `margin_z`.
  * `_finish`: b = `tau2_dn / m_dn` and `lam_dn = down_fractions / (b - 1)`.
  * `update`: xd = min(r², `tau2_dn`). The increase side is untouched.
* **`e11f_split_c.py`.** `first_alarms` reproduces the class (m14: 216 runs, no mismatch).
* **`panel.py`.**
  * `ebh`: R = max{r : r e_(r) ≥ K/α}, and the rejected set is e ≥ K/(αR).
  * `change_weights`: the window and Pareto weights, with the tail h/(n + h).
  * `PanelEProcess.step`: own clock, frozen assets, the switched-off decrease side contributing a constant ½, and the cap at 1e250, which keeps the supermartingale property.
  * `run_block` (sticky freezing), `fdp_summary`, `calibrate_panel` (the cs baseline at δ/K) and `PanelMonitor`.
* **`e15_panel.py`.**
  * `metrics` counts the FDP as 0 when there is no discovery.
  * `null_mean` / `kappa2`, `simulate`, `run_horizon`, `run_panel` with its offline e-BH check, and `validate`.

### Mathematics (theory_panel.md)

* **The results.** Lemma 1 (as in Wang and Ramdas 2022), Proposition 1 (D4), Lemma 2 (the predictable change time makes ν - 1 a stopping time), Proposition 2 (data processing on disjoint windows), and Theorem 1 parts (1) to (4) (D5).
* **The remarks and constructions.**
  * Remarks 3 to 6, with α' = 0.0087.
  * The window identity E_n = M_n / h + (h - n)/h.
  * The Pareto telescoping sum, with 2 log 11 = 4.8.
  * The block statements of 6(b).

### Numbers in E15.md that match their CSVs

* **Q15a.**
  * The no-change table.
  * The worst cases 4.57 / 0.58 / 0.96, 7.83 / 0.67 / 0.67, 13.28 / 0.75 / 1.00, 33.3 / 1.17 / 8.0 and 0.44 / 0.012 / 0.115.
  * The standard error 0.29.
* **Q15b.**
  * The table; a mean of 5.0 points and a maximum of 21.3 points; 0.17; 99.88%.
  * 47 to 57, 51 to 66 and 41 to 42.5.
  * The cs table and the oracle values 93.3, 88.8, 94.1 and 88.4.
* **Q15c.** The tables.
* **Q15d and Q15e.**
  * The tables; 3.7% to 8.7%; 71% to 77%; 2 to 4 and 24 to 38 minutes; 72% to 78%.
  * 86% to 87% and 58% to 60%; 170 of 184; p-values 0.56 to 1.00; thresholds 33 to 48.
  * The Nifty table; 226 / 104; joint discoveries 329 / 144 / 4 and 35 / 212 / 28 (XRP alone 28).
  * 472 sessions of 374 returns each (checked against the data). The Nifty intraday profile uses past data only.

### Not checked

I did not check the cited papers against their sources in this round, and I did not rerun any experiment script in full.
