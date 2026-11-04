# Independent review, round 7: TH5 (limits for every valid procedure, the plug-in baseline when n, m and A grow together, the decrease bet at its boundary)

> *Note added 3 November 2026: the "independent reviewer" or "referee" below was a separate AI agent (Claude, by Anthropic) instructed to act as a sceptical referee. It was not a human expert, and the review has not been confirmed by a human.*

**Reviewer role.** A sceptical referee for a statistics journal and a code auditor. I did not write any of the files below and did not take their claims on trust: every statement in the scope was re-derived, and the headline numbers were recomputed with my own code where that was cheap.

**Scope.**

1. Theory: `docs/theory_limits.md` (Lemma 1, Theorem 2, Corollaries 3 and 4, Theorem 5, Proposition 6, Lemma 7, Theorem 8, Proposition 9, Corollary 10, Proposition 11, Proposition 12 and the numerical sections).
2. Code and tests: `src/avmonitor/limits.py`, `tests/test_limits.py`, `experiments/th5_limits.py`.
3. Results: `results/th5_*.csv` (all ten files).
4. Summary: `docs/report_parts/TH5.md`.

**Read first.** `docs/math.md` sections 5 and 6 (the null classes, Theorems A and B, Theorem 6.4), `docs/theory_delay.md` (Lemma 2 and 2b, Corollary 4, Theorems 1 and 2), and `docs/theory_robustness.md` (TH4: Lemma 1, Theorem 1, Proposition 1, Corollary 1), on which TH5 builds.

**State of the tree.** TH5 was committed in ff6efe4 and integrated in d88fcef and bbc013b. While I worked, the branch moved to 9be152c (TH6); `git diff bbc013b origin/claude/project-thread-ytr1da` touches none of the files in scope, `results/` or `README.md`, so every line number below holds at both commits. Other agents had modified and untracked files in the working tree; none is in scope.

**What I ran.** Scripts are in `/tmp/claude-0/-home-user/e30ce406-b0e8-5b23-b612-d49eb104d0ee/scratchpad/rv7th5/` (outside the repository); nothing was written to `results/`. Everything ran under `nice -n 10`, one process at a time, a few minutes in total.

* `python -m pytest -q tests/test_limits.py`: **29 passed in 3.8 s.**
* `indep1.py`: B, the collapse and explosion floors, the Shiryaev-Roberts floors and the default grid's speed limits for c = 1.5 to 5 at A and 2A, from closed forms, without the repository's code. **Every entry of both tables of section 2.3 and of `th5_collapse.csv` agrees.**
* `indep2.py`: `q_B`, `ε*`, `λ* / λ_max`, `KL_inf`, `G_c`, the growth at ε = 10⁻⁶ and the total mass of `P*` for the nine decrease laws, by adaptive quadrature over the exact law of `|r|` (my own densities, not `design.y_law`). **All nine rows of the table of section 5 agree to the digits printed.**
* `indep3.py`, `indep4.py`, `indep5.py`, `indep9.py`: σ_Y, B and C_F for Gaussian and Student t3 returns, the constants of Theorem 8(b), and every cell of `th5_distribution_free.csv` by two independent integrations. **σ_Y, B, C_F, k1, k2 and C_DF agree; 13 cells of the distribution-free table do not (M2).**
* `indep6.py`: my own simulation of the pathwise floor, the deterministic Neyman-Pearson bound at `P*` (my own greedy test), the default grid and the best single bet, for the Gaussian halving (10,000 paths) and Student t3 at sd ratio 0.8 (6000 paths), new seeds.
* `indep7.py`: three fresh 20,000-window estimates of the hold-then-bold lower bound and of Corollary 1 at n = m = 390, with standard errors; the mass of `P*` from `design.best_bet` and from `limits.least_favourable_bet`.
* `indep8.py`: a linear programme for the fastest randomised valid rule on the complete collapse (M1).
* pandas checks of every range and count quoted in sections 3 and 4 of the note against the CSV files.
* Literature: arxiv.org is blocked from this machine and I did not try to get around it, so I could not open Clerico (2025), Waudby-Smith, Sandoval and Jordan (2025) or Ram and Ramdas (2026). I could not open Lorden (1971) or Lai (1998) either; my remarks on them rest on the statements in the note and on standard accounts of those papers.

**Severity.**

* *Critical*: a false theorem on which a reported number depends, or a wrong headline number.
* *Major*: a statement that is false as written, a gap in a proof, a claim labelled proved that is not, a wrong reported number outside the headlines, or an overclaim.
* *Minor*: a small wrong number, a weak test, a missing condition that changes no number, wording.

---

## 1. Summary verdict

The mathematics of TH5 is sound. Lemma 1, Theorem 2 (deterministic, randomised, sharpness), Corollary 3, Theorem 5, Proposition 6, Lemma 7, Theorem 8 (a) to (c), Proposition 9, Corollary 10, Proposition 11 and Proposition 12 are correct as stated and proved, apart from two small gaps in the proofs of Theorem 8c and Corollary 10 (m4, m5). The least favourable law `P*` is a valid law in all three null classes, and the Neyman-Pearson step is valid for a sequential rule under the composite, non-i.i.d. null, because it uses one null law and one fixed n at a time. Every headline number I was asked to check reproduces: 45, 50, 57.7, 73 and 79 at c = 3; the ranges 1.69 to 2.16 and 1.34 to 1.76; 33.8%, 42.4% and 36.5%; 5.14%; `ε* = 2.4 × 10⁻⁸`; `λ* / λ_max` from 0.93 to 0.55.

**No Critical finding. Five Major findings**, all fixable without changing a theorem's proof:

* **M1.** Corollary 4 and every summary of it say that *every* valid procedure needs 45 steps on the complete collapse. That is true only for deterministic rules. A randomised valid rule can alarm by step 44 with probability 0.978, and its smallest possible mean delay on that path is 43.4.
* **M2.** The code under-evaluates the distribution-free bound of Theorem 8(b) for n ≥ 10⁶. Thirteen cells of the table in section 4.7 are too small, which makes the reported bound look tighter than it is: for example 1.07 instead of 1.15 at n = 10⁸, m = 390, and 1.67 instead of 3.94 at n = 10⁷, m = 390.
* **M3.** "Proved to be at least 34% and at most 42%" is a Monte Carlo average over 4000 windows. The same quantity over 20,000 windows, already in the note, is 32.3%.
* **M4.** The summaries of the plug-in dichotomy drop the hypotheses of Corollary 10: "whatever A" is false, and "this is sharp" refers to the scale √n, not to the factor `1 + C m / √n`.
* **M5.** The "Reading" paragraph of `TH5.md` treats a lower bound as the best achievable delay: "1.7 to 2.2 times slower than the best any valid procedure could do".

---

## 2. Findings ranked by severity

### Major

**M1. The 45-step floor holds for deterministic rules only; Corollary 4 and the summaries say "every rule" and "every procedure".**

* *Where.* `docs/theory_limits.md:81` (Corollary 4: "every rule valid at levels `n / A` over `H↓` (or over `H=`) satisfies `T ≥ n_↓`", and the same for `n_↑`); `:14` ("no procedure valid at the levels `n / A` can alarm before step 45", and "no valid procedure alarms before step A" for `B = ∞`); `:105` ("Every rule valid at `n / A` over `H=` needs 45 steps"); `docs/report_parts/TH5.md:11` ("45 steps for every procedure", "no valid procedure alarms before step 45", "no valid procedure alarms before step A"); `TH5.md:113` ("every procedure with the guarantee needs 45 steps"). Outside the review scope, `README.md:26` makes it explicit: "A rule with the guarantee of Theorem B, randomised or not, ... So ... no valid procedure alarms before step 45".
* *What is wrong.* Section 1 of the note defines "rule" to include randomised rules, and Corollary 4's proof uses only Theorem 2(a), which is for deterministic rules. Theorem 2(b) itself shows that the corollary fails for randomised rules. On the collapse at c = 3, A = 7800, `W_n = (B / (B - 1))^n`, and the singleton randomised rule "alarm at step n if `y_{1:n} = 0` and `U ≤ p`" is valid whenever `p ≤ n W_n / A`. With n = 44, `p = 0.978` is allowed, so a valid rule alarms at step 44 with probability 0.978. The largest probability of an alarm by step n is `min(1, n W_n / A)`: 0.129 at n = 30, 0.556 at n = 40, 0.850 at n = 43, 0.978 at n = 44. I minimised `Σ_j j q_j` subject to `Σ_{j ≤ k} q_j / W_j ≤ k / A` for every k (the product law of Theorem 2c attains every `P(Y_{1:j} = 0) = 1 / W_j` at once, so these constraints are exactly validity on this path). **The smallest mean delay of any randomised valid rule on the complete collapse is 43.40**, with 74% of the mass at step 42. The explosion side behaves the same way: probability 0.28 at step 3 against the deterministic floor of 4. The deterministic numbers 45 (50) and 4 (4) are correct, and the practical difference is under two steps. But Corollary 4 is labelled [PROVED] and, as written, is false for part of the class it names. `README.md:26` says "randomised or not" in the same sentence.
* *Fix.* Write "every deterministic rule" in Corollary 4, at `:14`, `:105`, `TH5.md:11` and `TH5.md:113`. Add one sentence: "a randomised valid rule alarms by step n on this path with probability at most `min(1, n W_n / A)` (0.98 at step 44), and its mean delay on the path is at least 43.4". Apply the same change to Corollary 3's summary ("no deterministic valid rule ..."), and remove "randomised or not" from `README.md:26`.

**M2. The distribution-free bound (Theorem 8b) is under-evaluated for n ≥ 10⁶, which makes it look tighter than it is.**

* *Where.* `src/avmonitor/limits.py:320` to `335` (`dkw_log_moment`); `results/th5_distribution_free.csv`; the table at `docs/theory_limits.md:405` to `414` and the sentences at `:271` and `:416`; `docs/report_parts/TH5.md:76`.
* *What is wrong.* `dkw_log_moment` integrates `exp(h(x) - hmax)` over all of [0, 1] with `scipy.integrate.quad`, with breakpoints only at `x0 = √(log 2 / (2n))` and `t / (4n)`. For large n, the mass of the integrand sits in an interval of width a few `n^{-1/2}` (10⁻⁴ at n = 10⁸) at the left end of a unit interval. The adaptive rule then under-resolves the sub-Gaussian tail beyond `x0`. I computed the bound two ways: by substituting `u = √n x` and integrating piecewise on `[0, u0]`, `[u0, c/4]` and `[c/4, c/4 + 30]`, and by a 4-million-point trapezoid rule in log space. The trapezoid rule agrees with the CSV to four digits in every cell with n ≤ 92,160, and at n ≥ 10⁶ the two methods agree with each other to four digits. Both disagree with the CSV here (ratio to nominal):

  | n | m | CSV and note | independent |
  |---|---|---|---|
  | 10⁶ | 5 | 1.0075 | 1.0118 |
  | 10⁶ | 10 | 1.0146 | 1.0238 |
  | 10⁷ | 5 | 1.0023 | 1.0035 |
  | 10⁷ | 10 | 1.0043 | 1.0066 |
  | 10⁷ | 39 | 1.0173 | 1.0286 |
  | 10⁷ | 98 | 1.0532 | 1.1018 |
  | 10⁷ | 390 | 1.6704 | 3.9409 |
  | 10⁸ | 5 | 1.0007 | 1.0011 |
  | 10⁸ | 10 | 1.0013 | 1.0020 |
  | 10⁸ | 39 | 1.0050 | 1.0077 |
  | 10⁸ | 98 | 1.0131 | 1.0213 |
  | 10⁸ | 390 | 1.0729 | 1.1491 |
  | 10⁸ | 1950 | 10.2102 | 29.9886 |

  Restricting the same `quad` call to `[0, min(1, max(x0, t / (4n)) + 12 / √n)]` reproduces the independent values, which confirms the cause. A sanity check that needs no code: at n = 10⁶, m = 5, the leading term alone is `1 + (6 / 4) × 8.045 × 0.888 / 1000 = 1.0107`, already above the CSV's 1.0075. The theorem is unaffected; the bound is valid, just larger than reported. The qualitative reading ("informative only for m far below √n") survives and is if anything strengthened. But the quoted values "1.013 at m = 98" (`:271`), "(1.073)" (`:416`) and "1.07 at m = 390 for n = 10⁸" (`TH5.md:76`) are wrong, in the direction that flatters the bound.
* *Fix.* Integrate in `u = √n x` on `[0, min(√n, max(u0, c/4) + 30)]` with breakpoints at `u0` and `c/4`, or cap the upper limit as above. Regenerate the table, and change `:271`, `:416` and `TH5.md:76` to 1.021 and 1.15. Add a test of `dkw_log_moment` against the u-integral at n = 10⁷ (see m7).

**M3. "Proved to be at least 34% and at most 42%" is a Monte Carlo estimate from the smaller of two samples.**

* *Where.* `docs/report_parts/TH5.md:13` ("the largest false alarm probability over the class is proved to be at least 34% and at most 42%"). The same figure appears in `README.md:26` ("between 34% and 42%").
* *What is wrong.* Both ends are averages over calibration windows of per-window bounds. Each per-window bound is proved, but the average over windows is estimated. The 33.8% and 42.4% come from part `holdbold`, with 4000 windows (`th5_holdbold.csv`). The same two quantities are computed in part `plugin_bound` with 20,000 windows and are already in the note: 6.460 × 5% = 32.3% and 8.316 × 5% = 41.6% (`docs/theory_limits.md:362`). Three fresh 20,000-window runs with the repository's own functions gave 32.9%, 32.7% and 32.4% for the lower bound (standard error 0.33 points) and 42.0%, 41.9% and 41.7% for the upper bound (0.30 points). So the population lower bound is about 32.7%. "At least 34%" is neither proved nor right in its second digit. The sandwich of the simulated 36.5% between 33.8% and 42.4% on the *same* 4000 windows (`TH5.md:75`, `theory_limits.md:420` to `427`) is a correct check and is not affected.
* *Fix.* Write: "at least about 33% and at most about 42% (window averages of proved per-window bounds, Monte Carlo over 20,000 windows, standard errors 0.3 points)". Keep the 4000-window figures only for the sandwich.

**M4. The summaries of the plug-in dichotomy drop the hypotheses of Corollary 10.**

* *Where.* `docs/theory_limits.md:17` and `:303`; `docs/report_parts/TH5.md:13` and `:115`.
* *What is wrong.*
  1. Corollary 10(ii) needs `m / log A → ∞` and `m / A → 0`, not only `m / √n → ∞`. The growth condition is necessary, not technical. Every per-step factor of the default detector is at most `1 + 0.5 (B̂ - 1) ≤ 20.3` (increase side, `B̂ ≤ 39.57`) or `1 + 0.95 / (B̂ - 1) ≤ 2.05` (decrease side). So no stream in `C(W)` can alarm before about `log(2A) / log 20.3` steps, and if log A grows faster than a fixed multiple of m, `P(T ≤ m) = 0` for every stream. "The boundary is `m ≍ √n`, whatever A" (`:17`, `TH5.md:13`) is therefore false as stated. The corollary itself (`:289` to `:295`) states the hypotheses correctly.
  2. "The bound of Theorem B is then inflated by a factor `1 + C m / √n` for every stream with the right conditional mean, and this is sharp" (`TH5.md:115`). The factor is a first-order upper bound valid only when `m / √n → 0`. At `m / √n = 5.1` Corollary 1 is 2.24 against the first-order 1.19 (`:399`). Nothing proves that the excess `C_F (m + 1) / √n` is attained: in the cells with `m / √n < 1`, Proposition 9's lower bound is 0.004 to 0.005 of nominal (`th5_plugin_bound.csv`, m = 98, n ≥ 11,520). What is sharp is the scale `m ≍ √n` at which the nominal level breaks.
  3. `TH5.md:13` calls `C(W)` "the null class". It is the class of streams whose conditional truncated mean at τ̂ equals `m_F(τ̂)`, not the null class of Theorem B.
* *Fix.* At `:17` and `TH5.md:13`, write "when `m / √n → ∞`, `m / log A → ∞` and `m / A → 0`", and replace "whatever A" by "for any A with `log A = o(m)`". At `TH5.md:115`, write "inflated by at most `1 + C (m + 1) / √n` to first order when m is small compared with √n; the scale √n is sharp: for `m ≫ √n` (and `m ≫ log A`) some such stream breaks the nominal level by an unbounded factor".

**M5. The "Reading" paragraph turns a lower bound into the best achievable delay.**

* *Where.* `docs/report_parts/TH5.md:113` ("the e-detectors pay about twice the minimum"; "the default grid is 1.7 to 2.2 times slower than the best any valid procedure could do on average"; "the slowness of decrease detection is mostly the problem"). Outside the scope, `README.md:26` has the same claim ("1.7 to 2.2 times the smallest mean delay that any valid procedure could achieve").
* *What is wrong.* Theorem 5 gives a lower bound. Whether any valid rule attains it is open, as the note says itself (`:216`: "how much of the remaining gap is this price, and how much is the e-detector's own loss, is not settled here"; `:522`: "a rule that attains the delay floor uniformly over change times" is not done). The correct statement is "at most 1.7 to 2.2 times the best valid deterministic procedure". The bound is also for deterministic rules. Against the randomised Neyman-Pearson bound (`np_rand_A`), the ratios are 1.81 to 2.24 for the grid and 1.53 to 1.80 for the best bet, so "any valid procedure" needs "deterministic" or these wider ranges. "Mostly the problem" is not supported where the ratio exceeds 2: at Gaussian sd ratio 0.8 the data show only that at least 46% (1 / 2.16) of the grid's mean delay is unavoidable. The heading at `TH5.md:12` ("within a factor of about 2 of the best possible delay") is correct, because it bounds the inefficiency from above.
* *Fix.* Write "pay at most about twice the minimum", "at most 1.7 to 2.2 times slower than the best valid deterministic procedure (at most 1.8 to 2.2 against randomised ones)", and "at least half of the delay is the problem's, the rest is an unknown split between the detector and the slack in the bound". Apply the same to `README.md:26`.

### Minor

* **m1.** `docs/theory_limits.md:184`: "only 0.73 to 0.86 of `log A / γ_Q`" at A = 10⁸. The CSV gives 0.717 (Gaussian, sd ratio 0.8) to 0.863, so write 0.72 to 0.86.
* **m2.** `docs/theory_limits.md:471` and `TH5.md:93`: "at or below nominal in every cell (largest: 5.14% ...)". The largest point estimate is above 5%. Write "not significantly above nominal (largest 5.14%, interval 4.72% to 5.59%)".
* **m3.** `docs/theory_limits.md:212` and `TH5.md:46`, "at 2A (each side's level in the two-sided monitor) 1.65 to 2.09". These ratios compare one side at 2A with the bound for rules valid at `n / (2A)`. The two-sided monitor's guarantee is validity at `n / A` over `H=`. Against the bound at A, its decrease or increase side at 2A takes **1.85 to 2.40** times the bound (my computation from `th5_floors.csv`), which is the comparison made for the collapse at `:105` (79 against 45). Moreover, over `H=` the pathwise floor uses `π= = π↓ π↑`, which also charges values above 1, so it is larger than the one-sided floor on these laws; it was not computed. State which comparison is meant.
* **m4.** `docs/theory_limits.md:269`: "`V_n ⇒ (λ̄ + φ̄) σ_Y × (Z⁺ or Z⁻)`" is not a well-defined limit. The limit is `σ_Y (λ̄ Z⁺ + φ̄ Z⁻)`. Its mean, `(λ̄ + φ̄) σ_Y / √(2π)`, and therefore C_F, are unchanged.
* **m5.** `docs/theory_limits.md:297` (proof of Corollary 10(ii)): ρ̄ is taken on the compact set `|d| ≤ 1/2`, but `E_n` does not restrict d. Intersect `E_n` with `{|d| ≤ 1/2}`; its probability still tends to `P(|N(0, σ_Y²)| ≥ z₀)`, because `d → 0` in probability.
* **m6.** `docs/theory_limits.md:256` and `:266`: C_DF = 2.52 is valid but loose. The proof bounds `∫_0^∞ min{1, 2 e^{-2u²}} du` by `√(π/2) = 1.253`, while the integral equals `u0 + √(2π) Φ̄(2 u0) = 0.888` with `u0 = √(log 2 / 2)`. The same argument then gives **C_DF = 8.045 × 0.888 / 4 = 1.79**, 48 times C_F instead of 68.
* **m7.** Tests. (i) `tests/test_limits.py:145` to `159` builds `S` from `D.best_bet` for the Gaussian halving. That is the case where, as the note says, `P*` has mass 0.982 (I measured 0.98244; `least_favourable_bet` gives 1 to 3 × 10⁻¹¹). Use `L.least_favourable_bet`. (ii) `:201` asserts `r[-1] < 1.002` at n = 10⁸, m = 10. The correct value is 1.00200, so a corrected `dkw_log_moment` (M2) sits on the threshold, and the threshold appears to have been tuned to the faulty values. (iii) The moment test at `:202` to `207` checks only n = 390, where the function is accurate; add n = 10⁷.
* **m8.** `docs/theory_limits.md:486`: "`ε*` is proportional to the atom" is loose, because `ε* = q_B / h(ε*)` is an implicit equation. When ε* is small, `h(ε*)` is close to `h(0) = 1 - q_B - E[(B - 1)/(B - Y); Y < B]` (0.083 for the Gaussian halving, against `G_c = 0.088`). Say that.
* **m9.** `TH5.md:115`: "a one-day horizon is short compared with √n only for calibration windows of more than about 150,000 returns". At 152,100 returns `m = √n`, which is the boundary, not "short compared with". The note says "of order 10⁶ or more" (`:303`). Align the two.
* **m10.** `TH5.md:75`: "the proved lower bound for the worst stream is a large multiple of nominal ... (24 cells with `m / √n > 5` ...)". This omits that in 5 of these 24 cells the lower bound is *below* nominal (0.57 to 0.95), which the note reports (`:399`). Add the qualifier.
* **m11.** No Monte Carlo error is given for the Neyman-Pearson bounds or the efficiency ratios. My independent rerun of Student t3 at sd ratio 0.8 (6000 paths, new seed) gave a deterministic bound of 148.6 against 147.0 and ratios 1.918 and 1.764 against 1.91 and 1.76, so the second decimal is uncertain. A bootstrap over paths would settle the reported ranges.
* **m12.** `docs/theory_limits.md:224` and `:268`. That the estimation of τ does not enter at first order follows from the definition of `C(W)`, which pins the conditional truncated mean *at τ̂* to `m_F(τ̂)`. That is right for i.i.d. F streams. A null stated at the population level, `E[min(r_t², τ²) | F_{t-1}] = m₀`, would not give the cancellation. One sentence would prevent a misreading.
* **m13.** `docs/theory_limits.md:211` and `TH5.md:12`: the absence of pathwise violations on 180,000 paths is a necessary consequence of Theorem B and Theorem 2(a) for these valid detectors. It checks the code, not the theorem; call it a consistency check.
* **m14.** `docs/theory_limits.md:153` (citations). The use of Lorden (1971) and Lai (1998) is not load-bearing: Proposition 6(c) is proved in full in the note. A remark that would help the reader: validity at the levels `n / A` already implies `E_P T ≥ Σ_{n < A} (1 - n / A) = (A + 1) / 2` under every null law. So for deterministic rules the ARL-based bound quoted at Theorem 2b of `docs/theory_delay.md` gives the same first-order rate. The new content of (c) is randomised rules and the need for only one null law. I could not verify Lai's window-limited constraint against the paper (the note says the same), nor Clerico (2025), Waudby-Smith, Sandoval and Jordan (2025) or Ram and Ramdas (2026) (arxiv.org blocked).

---

## 3. Item by item

### Item 1. Lemma 1 and Theorem 2 (the pathwise floor): **correct**; "sharp path by path" is justified in the sense stated

* *Two-point laws in each class* (`:35` to `:41`). Decrease side, y < 1: the law `p δ_y + (1 - p) δ_B` with `p = (B - 1) / (B - y)` has support `{y, B} ⊂ [0, B]` and mean `py + (1 - p) B = B - (B - 1) = 1`, so it lies in `H↓`, `H↑` and `H=` at once. Increase side, y > 1: `(1/y) δ_y + (1 - 1/y) δ_0` has mean 1 and support in `[0, B]` because y ≤ B. For y ≥ 1 (decrease) or y ≤ 1 (increase) the point mass is in the one-sided class, and at y = 1 in `H=`. **Correct.**
* *π is the true maximum.* The upper bounds `1 ≤ E Y ≤ p y + (1 - p) B` and `1 ≥ E Y ≥ p y` are tight. For `H=` the class is smaller and the maximiser has mean exactly 1, so `π= = π↓ π↑`. The linear-programme test passes. **Correct.**
* *Theorem 2(a)* (`:51`, `:59`). `{T ≤ n}` is a set of values of `Y_{1:n}` and contains `y_{1:n}` when `T(y) ≤ n`, so `α_n ≥ P_y(T ≤ n) ≥ 1 / W_n(y)`. The product law is in the class because independence makes `E[Y_s | F_{s-1}] = E_{Q_s} Y`. **Correct.**
* *Randomised rules* (`:53`, `:61`). `φ(y_{1:n}) = P_U(T(y, U) ≤ n)` is a function of `y_{1:n}` because each `T(·, u)` is a stopping time, and `P_y(T ≤ n) = E_{P_y} φ ≥ φ(y_{1:n}) / W_n(y)`. This needs U independent of Y under `P_y` (assumed) and joint measurability of T (implicit). **Correct.** The consequence for Corollary 4 is M1.
* *Rules that see raw returns* (`:70`). With the null stated for returns, the laws on `{r_s, ±τ}` (decrease) and `{r_s, 0}` (increase) have the same masses and truncated mean exactly `m₀`. Any law of returns gives `P(r_s = x | E_{s-1}) ≤ P(Y_s = y(x) | E_{s-1}) ≤ π(y(x))`. **Correct.** "Other information whose law the null leaves free" is vague. It holds if the null puts no constraint on side information Z, so that `P_y` can make Z deterministic and equal to the observed z.
* *Sharpness (c)* (`:55`, `:63`). The chain-rule argument with elementary conditional probabilities given `E_{s-1}` is right: `E_{s-1} ∈ F_{s-1}` and `E[Y_s 1_{E_{s-1}}] = E[E[Y_s | F_{s-1}] 1_{E_{s-1}}] ≥ P(E_{s-1})`. The singleton rule is valid when α is nondecreasing. So `N_α(y)` is the pointwise minimum of `T(y)` over deterministic valid rules. The attaining rule depends on the path, and no single rule attains the floor on every path. The note says so (`:103`, the block rule), so "sharp path by path" is accurate. `E_Q[N_A]` is a bound that no single rule need attain.

### Item 2. Corollaries 3 and 4, and the numbers: **correct for deterministic rules**; scope wrong as stated (M1)

* *Corollary 3* (`:75` to `:77`). `H↓_∞` contains `H↓(B')` for every finite `B' > max_{s ≤ n} y_s`, and `W_n(y; B') → 1`. Correct for deterministic rules (no stop while `α_n < 1`) and for randomised rules (probability at most `α_n`). The summaries drop "deterministic" (M1).
* *B = 9.045.* Here `τ = c × 1.4826 × Φ⁻¹(0.75) = 2.9999955` at c = 3, because the scaled MAD of N(0, 1) is 0.9999985. Then `m₀ = E min(Z², τ²) = 0.995007` and `B = τ² / m₀ = 9.04513`. With τ = c exactly, `c² / κ₂(3) = 9.04516` (the value the test uses). The difference is in the fifth digit and changes no floor.
* *The numbers at c = 3, A = 7800* (my closed forms, `indep1.py`). `π↓(0) = 0.88944`, `log(B / (B - 1)) = 0.117160`.
  * Collapse floor: n = 44 gives `5.155 < log(7800 / 44) = 5.178`, n = 45 gives `5.272 ≥ 5.155`, so **45**. At 2A, **50**. Fractional roots 44.16 and 49.16.
  * Shiryaev-Roberts floor: `log(1 + A / B) / log(B / (B - 1)) = 6.7608 / 0.11716 =` **57.71** (63.62 at 2A).
  * Default grid: factors `1 + s / (B - 1)` with s = 0.1, 0.25, 0.5, 0.75, 0.95, weights 1/5, recursion until the mean reaches A: **73** (79 at 2A).
  * Explosion: 4 (4), Shiryaev-Roberts 5 (5), grid 7 (7).

  All rows for c = 1.5 to 5 in both tables agree.
* *Side checks.* The block rule (`:103`) gives `(A / 45) × 0.88944^45 = 0.890 ≤ 1`, and a collapse starting just after a block start is found up to 89 steps later. The single-horizon floor is `log 20 / 0.11716 = 25.6`, so 26 (`:71`). The gap is `log(45 / B) / log F = 13.7` steps (`:103`).

### Item 3. Theorem 5: **correct**

* *P* is a null law* (`:127`, `:131`). At an interior maximiser, `E_Q[U / (1 + λ* U)] = 0`, so `E_Q[1 / (1 + λ* U)] = 1` and `E_{P*} U = 0`. `P*` is a probability law with mean of Y exactly 1, so it lies in `H=`, hence in `H↓` and `H↑`. `1 + λ* U` is bounded and positive, so `P* ≡ Q` and the condition `Q ≪ P*` holds. Interiority is right for all 18 laws: the decrease laws have an atom at B, so `g'(λ) → -∞` as `λ → λ_max`, and the increase laws have `E_Q[1/Y] = ∞`. I checked the mass of `P*` under the *continuous* law, with my own quadrature and root, for all nine decrease laws: it is 1 to within 2 × 10⁻¹⁰. **Correct.**
* *The Neyman-Pearson step* (`:121` to `:125`, `:131`). It is valid for a sequential rule under a composite, non-i.i.d. null, because it is a relaxation. For each fixed n, `φ = 1{T ≤ n}` (or `P_U(T ≤ n | y)`) is `F_n`-measurable, satisfies `E_{P*} φ ≤ n / A` because `P* ∈ H`, and satisfies the cap by Theorem 2. Dropping all other null laws and the coupling between different n can only enlarge the feasible set. The greedy test is optimal by the pointwise inequality `(φ* - φ)(e^{Λ_n} - e^k) ≥ 0`, which uses `Q ≪ P*`. **Correct**, and loose for the reason the note gives (`:216`).
* *`E T ≥ 1 + Σ (1 - β_n)`.* It follows from `E_Q T = Σ_{n ≥ 0} Q(T > n)` and `Q(T > 0) = 1`. `β_n ≤ 1` by construction, since `0 ≤ φ ≤ 1`. In the code, `np_cdf` caps at 1 and sets every later `β_n` to 1 once one exceeds 0.99999, and `delay_lower_bound` sums only n = 1 to H - 1. Both changes lower the bound, so both are conservative. The deterministic cap `1{W_n ≥ A / n}` equals `1{N_A ≤ n}` because `W_n` is nondecreasing and `A / n` decreasing. As a consequence `NP det ≥ E[N_A]`, which the CSV shows. **Correct.**
* *Numbers.* My own implementation of the floor, the greedy test at `P*` and the two detectors, on new paths: Gaussian halving, floor 55.16, deterministic bound 55.16, grid 93.11, best bet 74.18, ratios **1.688 and 1.345**; Student t3 at 0.8, floor 47.74, bound 148.6, grid 285.1, best bet 262.1, ratios **1.918 and 1.764**. The note reports 1.69, 1.34, 1.91 and 1.76. The ranges 1.69 to 2.16 and 1.34 to 1.76 (A) and 1.65 to 2.09 and 1.33 to 1.72 (2A) match `th5_floors.csv`. The randomised gap of at most 15% (Gaussian halving, 14.9%) matches too. The E11 consistency check holds (93.19 and 122.07 in `e11a_delay.csv`).

### Item 4. Proposition 6: **correct**

* *(a)* Both inequalities at `N_A` and `N_A - 1`, together with the strong law and `log N_A / N_A → 0`, give `log A / N_A → γ_Q`. For L¹, `N_A ≤ τ_A`, the renewal theorem holds for i.i.d. nonnegative increments with positive mean, and the generalised dominated convergence theorem (Pratt) applies. **Correct.**
* *(b)* `π(y) f_λ(y) ≤ E f_λ ≤ 1` under the two-point law gives `log f_λ ≤ ξ`, hence `KL_inf ≤ γ_Q`. Equality forces `log f_{λ*}(Y) = ξ(Y)` Q-almost surely, which fails on `{Y > 1}` when `λ* > 0`, and `λ* > 0` holds because `γ_Q > 0`. The converse with `λ* = λ_max` holds on the closed interval. **Correct**, with equality exactly as stated.
* *(c)* The change of measure with `a = (1 - ε/2) log A` and `n_A = ⌊(1 - ε) log A / KL_inf⌋` is right. `Λ_n` has bounded increments, so the weak law applies. **Correct.**
* *Citations.* "The change-of-measure argument behind the lower bounds of Lorden (1971) and Lai (1998)" is a fair description: Lorden's bound under an ARL constraint, Lai's general information bound with a window-limited false alarm constraint. I could not open either paper (see m14). The bibliographic details (Annals of Mathematical Statistics 42(6), 1897 to 1908; IEEE Transactions on Information Theory 44(7), 2917 to 2929) match what I know of them.
* *Table of section 3.2.* It matches `th5_asymptotic.csv` except "0.73" (m1).

### Item 5. Theorem 8: **correct with gaps** (m4, m5, m6, m12); the code under-evaluates (b) (M2)

* *Drift bound* (`:236`, `:242`). For a nonnegative Z, `E min(Z, t) = ∫_0^t P(Z > x) dx`. So `|m̂(s) - m_F(s)| ≤ s² D_n` holds simultaneously for all s, and taking `s = τ̂` gives `|d| ≤ B̂ D_n`. The Kolmogorov distance of `|r|` equals that of `r²`. **Correct.** In the CSV the largest ratio is 0.55 for Gaussian and 0.76 for t3.
* *Range of B̂* (`:238`). With the sample median (the average of the two middle values for even n), at least half of the `|r_i|` lie on each side of it, so the bounds hold. Then `b_lo = 1.9038`, `b_hi = 39.566`, `k1 = 8.0451`, `k2 = 19.783`. **Correct.**
* *Lemma 7c* uses `(1 + x)^j - 1 ≤ j x e^{jx}`; summed over the bets it gives the factor `(m + 1) / 4` when only one side drifts. **Correct.**
* *C_F and the cancellation of τ̂* (`:258` to `:268`). Here `Δ(τ̂) - Δ(τ) = ∫_{τ²}^{τ̂²} (F^{(2)} - F̂^{(2)}) dx`, which is bounded by `|τ̂² - τ²| D_n = o_p(1) O_p(n^{-1/2})`. So `√n (m̂ - m_F(τ̂))` has the same limit as `√n Δ(τ)`, and the estimation of τ cancels at first order. Only consistency of τ̂ is needed, not a rate, and the reason is the definition of `C(W)` (m12). The limit of `V_n` is `σ_Y (λ̄ Z⁺ + φ̄ Z⁻)`, whose mean is `(λ̄ + φ̄) σ_Y / √(2π)` (m4). Uniform integrability holds through DKW once `m / √n ≤ 1`. I recomputed: Gaussian, **σ_Y = 1.38585, B = 9.0451, C_F = 0.03687**; Student t3, **σ_Y = 1.58267, B = 6.3625, C_F = 0.04711**. The Gaussian moments are from closed forms; the t3 moments are by quadrature. The simulated `√n sd(d)` agrees with σ_Y within 1.3% at every n.
* *C_DF and DKW-Massart* (`:256`, `:266`). The integration by parts `E g(D_n) = ∫_0^1 g'(x) P(D_n > x) dx` is right for `g(x) = x e^{tx}` and `D_n ≤ 1`. Massart's constant 2 is used for all x, which is fine because the bound exceeds 1 where Massart's condition fails. The reduction for laws with atoms, through the quantile transform, is standard. The rate `C_DF = k1 √(π/2) / 4 = 2.52` is valid but loose: 1.79 with the same argument (m6). The *values* in the table are wrong (M2).
* *The (1 + o(1)) regime* (`:260`). It needs `m x* → 0` in probability, which follows from `m / √n → 0`, and `(m + 1) / √n → 0`. Both hold. No condition on A is needed, because A enters only through the cap. **Correct.**
* *Uniformity over streams.* Corollary 1 of TH4 bounds `P(T ≤ m | F_0)` for every stream whose conditional truncated mean at τ̂ is `m_F(τ̂)`, including streams that depend on W, since `F_0` contains W. The bound depends only on the law of the window. So (a), (b) and (c) hold uniformly over `C(W)`. (b) is also uniform over all F with `P(r = 0) = 0`; (c) is not uniform in F, because its o(1) depends on F through the central limit theorem and the uniqueness of the median. **Correct as stated.**

### Item 6. Proposition 9, Corollary 10, Proposition 11: **correct**; the summaries overstate them (M4)

* *Proposition 9* (`:275` to `:285`). The stream is in `C(W)`: the hold value `m_F(τ̂) ≤ τ̂²`, and the bold law on `{0, τ̂}` has truncated mean `m_F(τ̂)`. The optional-stopping argument holds as written:
  * `V'_σ ≤ V_σ ≤ F V_{σ-1} < F K` on `{V_σ ≥ K}`, because `μ ≥ 1`;
  * `V'_L < K` implies `V'_L ≤ √(K V'_L)`;
  * `E √(V'_L) = √(R_h) ρ^L` by independence of the bold phase.

  `ρ < 1` by Jensen. I checked the code's ρ for both sides (`limits._bold_side`, `hold_bold_lower_bound`) against the formula: the win probability is `(1 + d) / B̂` for the increase side and `1 - (1 - |d|) / B̂` for the decrease side. **Correct.** The constant-stream claim is exact. Its decrease-side sum for d > 0 is the conservative value m, which cannot reach 2A > m.
* *Corollary 10* (`:289` to `:301`). (i) is Theorem 8c. (ii) holds with `L_A = O(log A)`, `h ≥ m/2`, `S_h ≥ h(1 + x (h + 1)/2)` (TH4 Proposition 1), and the final ratio `(w / (8F))(1 + β z₀ m / (4 √n))`. It needs `|d| ≤ 1/2` inside `E_n` (m5). (iii) holds with the deterministic lower bound `β ≥ 0.95 / (b_hi - 1)` on the decrease side. **Correct.** "Ratio → ∞" and "probability → 1" hold exactly under the stated hypotheses. The growth condition on A cannot be dropped (M4).
* *Proposition 11* (`:309` to `:317`). (a): on `G_z` both sides have non-positive drift relative to their own normalisations, which is TH4's band remark with separate `m↑` and `m↓`. (b): Slutsky, with `√n se → sd_F min(r², τ²)`; the bound `P(T ≤ m) - m/A ≤ P(G_z^c)` is uniform in m and A. (c): `P(m_F(τ̂) ≥ m↑ (1 + η / √n)) → Φ̄(z + η / σ_Y)`, then let η → 0. So `2 Φ̄(z)` is attained in the limit, and `z = Φ⁻¹(1 - δ/2)` (2.5758 for δ = 0.01) is the smallest fixed z with an asymptotic excess of at most δ in that regime. **Correct.** The shares 0.009 to 0.013 against `2 Φ̄(2.576) = 0.009995`, and the cs effective z from 3.78 to 5.74, match `th5_plugin_bound.csv` and `th5_plugin_drift.csv`.

### Item 7. Proposition 12: **correct**

* (a): `f_ε = ε` on the atom; on `{Y < B}`, `f_ε` is a convex combination of `f₀` and 1, which gives the lower bound, and `f_ε ≤ f₀ + ε` gives the upper bound. Every row of `th5_boundary.csv` satisfies both. (b): `g'(0+) = E[(B - 1) / (B - Y)] - 1`, which is `+∞` with an atom at B or a positive density near B. (c): `g'(ε) = q_B / ε - h(ε)`, so `ε* = q_B / h(ε*)` (m8 on wording). (d): direct substitution.
* (e): the Wald bound is Theorem 1 of `docs/theory_delay.md` with one bet and overshoot at most `log(1 + λ) ≤ log(B / (B - 1))`. The block bound holds because the Shiryaev-Roberts statistic is at least the product over any block, whatever the history, and disjoint blocks are independent under i.i.d. Q: `E T ≤ n Σ_k (1 - p_n)^k = n / p_n`. On the no-atom event `f_ε ≥ (1 - ε) f₀`, which gives a lower bound on `p_n` that does not depend on ε. **Correct.** The block bounds in the table are Monte Carlo estimates of the bound, as the code's docstring says. The table should say so too.

### Item 8. Numbers: **all requested numbers reproduce**; the distribution-free table does not (M2)

| Number (location) | Note / CSV | Independent | How |
|---|---|---|---|
| B at c = 3 (`:43`) | 9.045 | 9.04513 | closed form |
| `π↓(0)` (`:43`) | 0.8894 | 0.88944 | closed form |
| Collapse floor, A and 2A (`:90`) | 45 (50) | 45 (50) | closed form |
| Fractional root (`:90`) | 44.16 (49.16) | 44.16 (49.16) | root finding |
| SR floor, linear bets (`:90`) | 57.71 (63.62) | 57.71 (63.62) | closed form |
| Default grid on the collapse (`:90`) | 73 (79) | 73 (79) | own recursion |
| Explosion: floor, SR, grid (`:99`) | 4, 5, 7 | 4, 5, 7 | closed form |
| Grid / bound, best bet / bound at A (`:212`) | 1.69 to 2.16, 1.34 to 1.76 | CSV ranges equal; own sims 1.688, 1.345 (Gaussian 0.5) and 1.918, 1.764 (t3 0.8) | own simulation |
| Hold-then-bold lower bound, Corollary 1, simulated rate, n = m = 390 (`TH5.md:75`) | 33.8%, 42.4%, 36.5% | same in `th5_holdbold.csv`; 20,000 windows give 32.4% to 32.9% and 41.7% to 42.0% | repository functions, new seeds (M3) |
| Plug-in, t3, n = 390, one day (`:471`) | 5.14% (4.72% to 5.59%) | same in `th5_plugin_sim.csv` | pandas |
| `ε*`, Gaussian halving (`:498`) | 2.37 × 10⁻⁸ | 2.365 × 10⁻⁸ | own quadrature and root |
| `λ* / λ_max`, t3 at 0.5, 0.67, 0.8 (`:501` to `:503`) | 0.931, 0.777, 0.545 | 0.9310, 0.7767, 0.5451 | own quadrature |
| `q_B`, t3 (`:511`) | 0.65% to 2.4% | 0.649% to 2.383% | own quadrature |
| Growth at ε = 10⁻⁶, t3 at 0.67 and 0.8 (`:511`) | -0.126, -0.282 | -0.1260, -0.2817 | own quadrature |
| C_F, Gaussian and t3 (`:262`) | 0.0369, 0.0471 | 0.03687, 0.04711 | closed form and quadrature |
| k1, k2, C_DF (`:252`, `:256`) | 8.045, 19.78, 2.52 | 8.0451, 19.783, 2.5207 | closed form |
| Distribution-free, n = 10⁸, m = 390 (`:416`) | 1.073 | **1.149** | two integrations (M2) |
| Distribution-free, n = 10⁶, m = 39 and 98 (`:271`) | 1.153, 2.316 | 1.1525, 2.3159 | two integrations |
| `√n` = 390 (`:303`) | n = 152,100 | 390² = 152,100 | arithmetic |

I also checked these ranges and counts in the CSV files:

* section 4.7: "19 of the 24 cells" (cells with `m / √n > 5` and m < 7800), the ratio of lower to upper bound from 0.253 to 0.961 with median 0.674, the five exceptions at `m / √n` from 5.1 to 9.1, the short-horizon claims (4.1% and 5.3%, within 1% and within 5%), and the ratios to the true baseline of 3.4 down to 1.07 and 2.9 down to 1.06;
* section 3.2: within 0.4% to 5% of the refined root, 0.85 to 0.98, 1.16 to 1.19, 1.38 to 1.93 and 1.22 to 1.40;
* section 5: within 4.2%, 0.4%, 0.2% and 3%, 14% and 53%, and at most 1% for the jump law.

All agree, except "0.73" (m1).

### Item 9. Overclaims in `docs/report_parts/TH5.md`

| Line | Statement | Status |
|---|---|---|
| 11 | "45 steps for every procedure", "no valid procedure alarms before step 45 / before step A" | False for randomised rules (M1) |
| 13 | "The boundary is `m ≍ √n`, whatever A"; "`→ ∞` when `m / √n → ∞`" without `m / log A → ∞` | Hypotheses dropped (M4) |
| 13 | "proved to be at least 34% and at most 42%" | Monte Carlo average; 20,000 windows give 32.3% to 32.9% (M3) |
| 13 | "a stream in the null class" | It is in `C(W)`, not the null class of Theorem B (M4) |
| 75 | "large multiple of nominal ... 24 cells" | 5 of the 24 are below nominal (m10) |
| 76 | "1.07 at m = 390 for n = 10⁸" | 1.15 (M2) |
| 93 | "at or below nominal in every cell (largest 5.14%)" | Wording (m2) |
| 113 | "pay about twice the minimum", "1.7 to 2.2 times slower than the best any valid procedure could do", "mostly the problem" | Lower bound read as the optimum (M5) |
| 115 | "inflated by a factor `1 + C m / √n` ... and this is sharp" | The factor is first order for `m / √n → 0` only; the scale is sharp, not the factor (M4) |
| 115 | "more than about 150,000 returns" | That is where `m = √n` (m9) |

The labels are otherwise accurate. The Neyman-Pearson bounds are called Monte Carlo (`:12`, `:120`). The `log A / √n` scale for i.i.d. returns is called heuristic and post hoc (`:92`, `:122`). The i.i.d. simulation is called simulated. Nothing else is presented as proved when it is only simulated.

---

## 4. What I did not check

* The figures `results/figures/th5_floors.png` and `th5_plugin_regime.png`, and `docs/iterations/TH5.md`.
* Section 4.6's numbers that refer to `docs/theory_running_tau.md` (Lemma A.7's 10⁻⁵⁴⁶) and to E14 (Q14w).
* A full rerun of `experiments/th5_limits.py`. I reran pieces with new seeds instead, as described above.
* The content of the references that I could not open (m14).
