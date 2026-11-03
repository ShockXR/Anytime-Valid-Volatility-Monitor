# Second referee report on "Anytime-valid detection of volatility changes in high-frequency returns with robust truncation" (revised version)

Manuscript read: `paper/short_paper.md` at commit 25806aa (1121 lines; main text about 9,980 words by my count, then Supplement S1 to S9) and the companion note `paper/extensions.md`. First report: `docs/review/referee_report_paper.md` (points M1 to M8 and m1 to m35).

What I did. I read the main text and the supplement in full, and checked every numbered theorem, proposition and corollary against its proof (Supplement S1 to S6, and `docs/theory_per_period.md` and `docs/theory_limits.md` where the paper defers to them). I recomputed the closed-form constants, solved the linear programme behind the randomised collapse delay, and checked 36 reported numbers against the committed CSV files in `results/` or against closed forms. I ran `scripts/trace_numbers.py` (2379 numbers, 0 untraced; the script rewrites the tracked file `results/trace/paper_trace.csv`, which I restored). I could not open arxiv.org or CRAN from my environment and did not try to get around that; claims that depend on those sources are marked as unverified below.

Severity scale used for findings: **Critical** (a false statement or a wrong headline number), **Major** (must be fixed before acceptance), **Minor** (should be fixed; does not change a conclusion).

---

## 1. Summary and recommendation

The revision is a substantial improvement. Every theorem I checked is correct as proved, and I found no false theorem and no wrong headline number. The new later-day results (Proposition 7, Theorem 8, Proposition S6, Theorem S7, Proposition S8), the pathwise floor and delay floor (Theorems 12 and 23), the plug-in theorem with growing $n$, $m$ and $A$ (Theorem 15) and the confidence-sequence theorem with a stated tolerance (Theorem 16, with its assumption-free part 4) answer my main theoretical concerns. The real-data sections now separate primary from exploratory analyses, use Holm adjustment and a day-block bootstrap, add null streams built from real returns, and report negative results honestly.

Three matters remain, none of them a false statement:

1. **(J1) The real-data evidence still rests on the plug-in baseline**, which has no finite-sample guarantee, while the abstract sentence on FOMC statements reads as if it showed "what the guarantee buys".
2. **(J2) Novelty and positioning rest partly on papers the authors could not read** (and that I could not read either), but the text states their contents with the same confidence as checked facts.
3. **(J3) The main text uses objects defined or proved only outside the paper**: the rank detector, which is one of the three "proved monitors" of the abstract and of the primary null-stream endpoint, and the running-truncation restarts of Section 6.5.

**Recommendation: minor revision**, conditional on J1 to J3. For a theory-and-methods venue, fit is the editor's call: the theory is now correct and complete for what it claims, but several results are short consequences of Ville's inequality and the Neyman-Pearson lemma; the contribution is the combination, the explicit constants and the honest negative results rather than any single deep theorem.

---

## 2. Response to the major points of the first report

| Point | Verdict | Evidence in the revision | What remains |
|---|---|---|---|
| M1. Real-data claims rest on an estimator without a guarantee | **Partly** | Section 7.1 separates primary from exploratory analyses and adds null streams from permuted calm days (S7.6); Section 7.2 adds a pre-registered confidence-sequence configuration; S7.6 adds a forecast-relative null. | The forecasts, the Driscoll-Kraay regressions, the matched-competitor comparisons, the FOMC and CPI test, the Nifty analysis and the landmark shocks all still use the plug-in. No real-data run uses a configuration with a per-day guarantee. The band $b$ is not motivated. See J1. |
| M2. Over-many-days reading of the first-$m$-steps bound | **Answered** | Proposition 7 ($(\bar M_s + m)/A$ given $\mathcal F_s$; $dm/A$ for day $d$), Theorem 8 (restart, $(m+\kappa)/A$ given $\mathcal F_{b_i}$), S1.7 (Proposition S6, Theorem S7, Proposition S8, hazards to day 20), the Lai and Mei discussion. | Residual "per day" wording for first-day rates (Minor F7). |
| M3. Limits of any valid procedure | **Answered** for items (i) to (iii) | (i) Theorem 12, Corollary 13, Theorem 23, Proposition S11; (ii) Theorem 15; (iii) Proposition S12. | Items (iv) and (v) of M3 are still open and now listed as open problems, which is acceptable. Depth is modest (Section 1). |
| M4. Confidence-sequence baseline: the null was algorithm-dependent | **Answered** | Theorem 16 states a tolerance on the conditional truncated variance fixed in advance, with horizon $h$; part 4 needs no assumption; "What the tolerance excludes" paragraph; abstract and contribution 2 updated; relation to Shekhar and Ramdas stated. | Nothing material; one scope remark in Section 4. |
| M5. Competitors | **Partly** | Classical rules on the e-detector's statistic $X$ (Page CUSUM, likelihood-ratio CUSUM), simulated and per-path bootstrap thresholds, a Student $t$ BOCPD, HKZ-type and CSW-type training-sample rules (Section 6.2, S5). | No published Bibinger-Jirak-Vetter procedure, no beta-divergence robust BOCPD, and HKZ and CSW only "of the type of" since the papers could not be read. Acceptable if said in Section 6.2, which it is in S5 but not in the main text. |
| M6. Multiple comparisons and dependence | **Answered** | S7.4: Holm adjustment, day-block bootstrap, Driscoll-Kraay standard errors, QLIKE loss, HAR with the indicator. | Section 7.2 mixes Holm-adjusted and unadjusted counts (Minor F12). |
| M7. Related work | **Answered** | The related-work paragraph (line 28) covers retrospective tests, monitoring with a training sample, Bayesian online detection, e-detectors and conformal martingales. | Several statements about those papers are unverified (J2). |
| M8. Length and self-containment | **Partly** | Main text about 10,000 words (my count 9,977); experiment labels replaced by supplement references; provenance tags hidden as comments. | Abstract 216 words by the authors' count (231 by mine) against about 200. The rank detector and the running-truncation restarts are used in the main text but defined or proved only in the companion note or repository notes (J3). |

### Minor points m1 to m35

| Verdict | Points |
|---|---|
| Answered | m3 to m9, m11 to m21, m25, m27 to m29, m31 to m35 |
| Answered in substance, but the text now contradicts itself | m22: S1.7 checks the restarted monitor at $A = 7800$ and $28{,}800$ (0.037 to 0.048 false alarms per day) and S7.6 at $A = 36{,}000$, yet S1.1 (line 556) still says "Corollary 6 itself was not checked at an operating level". Delete that sentence or point to S1.7. |
| Partly | m1 (abstract still over 200 words, see above); m2 (notation still overloaded: see Section 8); m10 (`paper/short_paper.pdf` is older than the Markdown source and does not contain E19 to E21); m24 (the renewal conditions are discussed but the single-bet result is not stated as a proposition); m30 (DOIs and issue numbers still inconsistent) |
| Not answered | m23 (the share $n_\tau$ of the window used for $\tau$ is not varied); m26 (no experiment isolates errors in the intraday profile) |

---

## 3. Correctness of the stated results

### 3.1 Result by result

| Result | Verdict | Remarks |
|---|---|---|
| Lemma 1, Lemma 2, Theorems 3 and 4 | Correct | Known baseline; unconditional; first $m$ steps. |
| Proposition 5, Corollary 6 | Correct | Corollary 6 is marginal over days. |
| Proposition 7 | Correct | Given $\mathcal F_s$; $dm/A$ for day $d$, and no bound below 1 holds for all days (day-40 example). Assumes $b = 0$ and $Y = X/m_0$ with known $m_0$; say so in the statement. |
| Theorem 8 (restart) | Correct | Known baseline only; given $\mathcal F_{b_i}$. The window remark is checked in Section 6(b) below. |
| Theorem 9 (drift), Lemma S1, Proposition S2, Lemma S3, Theorems S4 and S5 | Correct | Table S1 recomputed exactly. |
| Corollary 14 (plug-in through the drift bound) | Correct | Bounds 42%, 7.3%, 6.0% recomputed. |
| Proposition 10, Corollary 11 (no e-process for "has not fallen") | Correct | |
| Theorem 12 (pathwise floor) | Correct | (a) and (c) for deterministic rules, (b) for rules randomised independently of the data. Known baseline. |
| Corollary 13 | Correct, one gap in the text | See Finding F3: the mean delay 43.4 does not follow from the cap stated in the corollary. |
| Theorem 15 (plug-in, growing $n$, $m$, $A$) | Correct | Hypotheses incomplete in the statement: see F4. |
| Theorem 16 (confidence-sequence baseline), Lemmas S9 and S10 | Correct | Parts 1 to 3 need the tolerance almost surely over the horizon $h$; part 4 gives $\delta + \mathbb P(V_m) + m/A$ with no assumption. Proof in S2 is complete. |
| Corollaries 17 and 18 | Correct | |
| Lemmas 19 and 20, Corollaries 21 and 22, Theorems 24 and 25 (growth rate, speed limits, delay) | Correct | As in round 1; constants recomputed. |
| Theorem 23 (delay floor) | Correct | (a) deterministic; (b) both deterministic and randomised. Proposition S11 correct. |
| Proposition S6 (renewal bound) | Correct | |
| Theorem S7 (CUSUM-type window bound) | Correct | Exchangeable null, fixed bets, known baseline, unconditional. "Known baseline" is missing from the statement (F5). |
| Proposition S8 (no window bound under martingale differences) | Correct in substance | Definition of $\rho$ ambiguous: see F6. |
| Propositions S12 and S13 | Correct | $\rho$ values of S13 recomputed. |
| Corollary 1 of `docs/theory_per_period.md` (restart with confidence-sequence baseline, $\delta + m/A$ per day) | Correct | Marginal only, not given the past; cited in the paper (Section 3.1, line 114) without a number. Make it a numbered result in the supplement. |

### 3.2 Scope of each guarantee

The paper now states most scopes, but they are spread over many paragraphs. A table of this kind in Section 3 or S1 would prevent the remaining misreadings (F7, F8):

| Result | Baseline | Rule class | Conditioning | Horizon |
|---|---|---|---|---|
| Theorems 3, 4 | known | the e-detector | unconditional | first $m$ steps |
| Proposition 7 | known | the e-detector | given $\mathcal F_s$ | any $m$ steps after $s$ |
| Theorem 8 | known | restarted e-detector | given $\mathcal F_{b_i}$ | each day |
| Corollary 1 (theory note) | confidence sequence | restarted e-detector | marginal | each day |
| Theorem S7 | known | CUSUM-type, fixed bets | unconditional, exchangeable null | every window of $m$ steps |
| Proposition S8 | known | any SR or CUSUM rule | martingale-difference null | counterexample |
| Theorem 12 | known | (a), (c) deterministic; (b) randomised | pathwise | all $n$ |
| Theorem 23 | known | (a) deterministic; (b) both | post-change law $Q$ | delay |
| Theorem 15 | plug-in | the e-detector | unconditional | asymptotic in $n$, $m$, $A$ |
| Theorem 16 | confidence sequence | the e-detector | unconditional | parts 1 to 3: tolerance a.s. over $h$; part 4: none |

### 3.3 Findings on the theory

**F1. Minor. Abstract (line 9), contributions 1 and 3 (lines 24, 26), Section 5.1 heading (line 240).** "Limits the speed of any valid rule" and "a delay floor for every valid rule" omit "deterministic", and "valid" is never defined outside Theorem 12. The hard floor (Theorem 12(a), Corollary 13(b), Theorem 23(a)) is for deterministic rules; randomised rules obey only the softer bound of Theorem 12(b) and Theorem 23(b). All floors also assume a known baseline. *Fix:* define "valid" once in the introduction as "$\mathbb P(T \le n) \le n/A$ for every $n$ and every law in the null, with a known baseline", and write "every deterministic valid rule (randomised rules obey a slightly weaker bound)".

**F2. Minor. Proposition 7 and Theorem 12 (lines 105, 143).** Both use $Y = X/m_0$ and the three nulls with $b = 0$, and Theorem 12 assumes the baseline known. This is stated in Section 2 but not in the statements. *Fix:* add "with a known baseline and $b = 0$" to both statements.

**F3. Minor. Corollary 13 (line 147) and Section 5.1 (line 237).** The sentence "A randomised valid rule alarms by step $n$ on that path with probability at most $\min(1, nW_n/A)$, which is 0.98 at step 44, and its mean delay there is at least 43.4" suggests that 43.4 follows from that cap. It does not: summing $1 - \min(1, nW_n/A)$ over $n$ gives only 37.8. The value 43.4 is correct, but it is the optimum of a linear programme that imposes validity at every $n$ under the i.i.d. two-point law on $\lbrace 0, B\rbrace$, namely $\sum_{k \le n} p_k (B/(B-1))^{-k} \le n/A$; my solution gives 43.39 (the repository's `docs/theory_limits.md` section 3 gives the same). The cap itself, 0.978 at step 44, is correct. *Fix:* write "its mean delay there is at least 43.4 (a linear programme over the stopping distribution, Supplement S1.6)" and give the programme in S1.6.

**F4. Minor. Theorem 15 (line 170).** The proof uses a unique median of the calibration law (for the consistency of the scaled MAD), but the statement does not. The class $\mathcal C(W)$ depends on the calibration window $W$, so the guarantee is relative to a null fixed after calibration, not to a population-level null; this is correct but should be said where the theorem is stated, since a reader will compare it with Theorem 16. The symbol $T$ is used in the statement without being defined there. *Fix:* add "the calibration law has a unique median", "conditionally on $W$", and define $T$.

**F5. Minor. Theorem S7 (line 663).** The statement omits "with a known baseline"; the proof uses it. *Fix:* add it.

**F6. Minor. Proposition S8 (line 669).** The park-and-strike construction parks at $Y = 1$, which sets the CUSUM-type statistic to $\max(C, 1)$. The definition of $\rho$ (the ratio that the counterexample drives to 1) does not say whether a side's running average may exceed $2A$ when parking starts within 1 of the level. The conclusion is unaffected (the construction can start parking below $2A - 1$), but the statement should be unambiguous. *Fix:* define $\rho$ explicitly and start the parking phase at $C \le 2A - 1$.

**What is correct (theory).** Lemma 2 and Ville's inequality give Theorems 3 and 4 exactly as stated. Proposition 7 is sharp in the sense now stated (no bound below 1 for all days), and the day-40 path is a valid null path. Theorem 8's conditioning on $\mathcal F_{b_i}$ is right, since $\lbrace T_A > b_i\rbrace \in \mathcal F_{b_i}$. The time-reversal proof of Theorem S7 is correct under exchangeability with fixed bets, and Proposition S8 shows that the exchangeability assumption cannot be dropped. Theorem 12(b) follows by integrating the deterministic bound over the independent randomisation. Theorem 15's dichotomy at $m \asymp \sqrt n$ and its constants (0.037, 0.047, 1.79 via the Dvoretzky-Kiefer-Wolfowitz-Massart inequality) are right. Theorem 16 part 4 is an honest assumption-free bound and is the right thing to have added. Corollary 1 of the per-period note is correct but marginal only.

---

## 4. Numbers checked

All 36 numbers below agree with their source to the stated precision. "Closed form" means I recomputed from the formula in the paper.

| # | Paper location | Number(s) | Source | Result |
|---|---|---|---|---|
| 1 | Section 2 | $\kappa_2 = 0.995$, $B = 9.045$ | closed form | agrees |
| 2 | Proposition 7 | certain alarm on day 40 | `th6_calm_family.csv`; closed form | agrees |
| 3 | S1.7 | hazards 0.98% rising to 4.70% | `th6_hazard_390.csv`, `th6_hazard_1440.csv` | agrees |
| 4 | S1.7 | 3.0%, 4.1%, 4.7% | same | agrees |
| 5 | S1.7 | 28.02%, 65.07% (day 10, 20; 390 steps) | `th6_calm_family.csv` | agrees |
| 6 | Section 3.1, Proposition S6 | delay ratios 1.06 to 1.44 and 1.02 to 1.14 | `th6_delay.csv` | agrees |
| 7 | Theorem S7 | CUSUM delay ratios 1.30 to 1.83 and 1.24 to 1.66 | `th6_delay.csv` | agrees |
| 8 | Theorem S7 | 5.146%, 5.044%; 8027.7 and 29,053.4 | closed form | agrees |
| 9 | S1.7 | bold play 5.25% and 5.02% | `th6_adversary.csv` | agrees |
| 10 | S1.7 | restarted 2.66% and 3.52% (hazard 2.75%) | `th6_adversary.csv` | agrees |
| 11 | S1.7 | CUSUM-type 0.18% and 0.40% | `th6_adversary.csv` | agrees |
| 12 | S1.7 | 0.037 to 0.048 false alarms per day | `th6_calm_family.csv` | agrees |
| 13 | Table S1 | drift thresholds 1.63%/2.66%, 1.72%/2.70%, 6.36%/10.51%, 0.44%/0.72% | closed form | agrees |
| 14 | S1.6 | collapse floors 23, 33, 45, 73, 106 (25, 36, 50, 82, 120 at $2A$) | `th5_floors.csv` | agrees |
| 15 | S1.6 | e-detector floors 28.7 to 140.8 | `th5_floors.csv` | agrees |
| 16 | S1.6 | explosion 4; single horizon 26; $\pi(0) = 0.889$ | `th5_floors.csv`; closed form | agrees |
| 17 | Corollary 13 | caps 0.98, 0.13, 0.56, 0.85 | closed form | agrees |
| 18 | Corollary 13 | randomised mean delay 43.4 | my linear programme: 43.39 | agrees (see F3) |
| 19 | Corollary 14 | bounds 42%, 7.3%, 6.0% | closed form | agrees |
| 20 | Theorem 15 | $C$ = 0.0369, 0.0471, 1.79 | closed form | agrees |
| 21 | Section 5 | $I_0 = 0.888$; $b_{lo} = 1.904$, $b_{hi} = 39.57$; $k_1 = 8.045$, $k_2 = 19.78$ | `e11b_information.csv`; closed form | agrees |
| 22 | S2 | non-asymptotic factors 1.15, 2.3, 1.021, 1.15 | `th4_plugin_caln.csv` | agrees |
| 23 | Theorem 15(b) | hold-then-bold 36.5% (35.0 to 38.0); bounds 33.8% and 42.4% | `th5_holdbold.csv` | agrees |
| 24 | Theorem 16 | tolerance medians 0.90 to 2.58 and 2.05 to 5.01 | `th6_tolerance.csv` | agrees |
| 25 | Section 5.3 | $A$ = 9750 and 36,000; ratio 1.025 | closed form | agrees |
| 26 | Theorem 23 | ratios 1.687 to 2.156 and 1.344 to 1.764 (46% and 45%); floor 55.2 against 45.6 | `e11b_information.csv` | agrees |
| 27 | Table 1 | all entries | `e04a_null.csv` | agrees |
| 28 | Table 2 | all entries | `e08b_change.csv` | agrees |
| 29 | Section 6.3 | confidence-sequence delay 1.8 to 7.6, 1.4 to 3.0, 1.2 to 1.7 times the plug-in | E16 results | agrees |
| 30 | Section 6.2 | calibrated rule 106, 41, 20 (e-detector 100, 39, 19); decrease 144 against 162 | `e21_change.csv` | agrees (see F11) |
| 31 | Proposition S13 | $\rho$ values | closed form | agrees |
| 32 | Section 7.1 | null streams 2.26% (2.06 to 2.48), 0.01%, rank 3.38%, plug-in 5.50% and 41% | `e19_null_pooled.csv` | agrees |
| 33 | Section 7.2 | $-1.3$ ($-5.8$ to 3.2) and $+8.0$ (2.7 to 13.7) | `e19_cs_primary_bootstrap.csv` | agrees |
| 34 | S7.4 | F1 differences 27 (22 to 33) and 38 (34 to 44) | `e20_bootstrap.csv` | agrees |
| 35 | S7.4 | ratio rule Holm 0.0096, 1 of 48 | `e20_holm_existing.csv` | agrees |
| 36 | Section 7.3, S7.2 | FOMC Holm 0.00025; CPI 0.11; ETH and SOL 0.0032 and 0.00015; Budget 0.0006 | `e20_holm_existing.csv` | agrees |

I also checked the Nifty precisions (78% and 76% against 47% and 53%) in S7.1 and the time-monitored shares (Section 6(d) below); both agree.

**A surprising number, explained.** In Section 6.2 the calibrated betting rule is slower than the e-detector for increases (106 against 100 steps) although its pre-change false alarm rate is higher (3.2% against 2.6%). I read the calibration code path: the calibrated rule splits its target first-day rate equally between the two sides, whereas the e-detector's realised first-day rate falls mostly on the increase side, so the calibrated rule's increase side is in effect held to a lower rate. The number is right; the explanation should be one sentence in S5.

---

## 5. Overclaims and wording

**F7. Minor. "Per day" for first-day rates.** Lines 309, 328, 904, 906, 908, and Table S11 (lines 1074, 1081) call first-day (first-$m$-steps) rates "false alarms per day". After Proposition 7 the paper itself shows that later days can be much worse for the plain detector. *Fix:* write "on the first day" or "in the first 390 steps" wherever the rate is a first-day rate, and keep "per day" for the restarted monitor.

**F8. Minor. Abstract, last sentence (line 9).** "Simulations, these data and FOMC statements show what the guarantee buys and what it does not." The FOMC and CPI test (Section 7.3) used the plug-in, which has no finite-sample guarantee. *Fix:* "show what the monitor does and does not deliver", or run the primary test with the confidence-sequence baseline as well (part of J1).

**F9. Minor. Contribution 1 (line 24).** The guarantees are stated without "with a known baseline". *Fix:* add it; the estimated-baseline story is contribution 2.

**F10. Minor. Contribution 2 (line 25).** "The plug-in baseline keeps its level for horizons short compared with the square root of the calibration length" is asymptotic. With one 390-step day this needs $n$ well above $390^2 = 152{,}100$ calibration steps, far more than the 2000 used. *Fix:* "asymptotically, for horizons $m = o(\sqrt n)$", and add one sentence on what this means at the operating point.

**F11. Minor. Section 6.2 (line 328) and S5 (line 904).** "Per-path bootstrap thresholds held 1.4% to 2.2% (target 1.3%) and 4.9% to 6.3% (target 5%)": in `e21_null.csv`, the matched bootstrap CUSUM on $X$ was significantly above target (lower 95% limit above 1.305%) under Student $t_3$ (1.85%), jumps (1.85%) and price noise (2.2%), and at the 5% target the price-noise cells (6.05% and 6.3%) were also above. Five cells in all, among the laws without clustering. *Fix:* "were 1.4% to 2.2% against 1.3% and 4.9% to 6.3% against 5%, significantly above target in 5 of 24 cells".

**F12. Minor. Section 7.2 (line 376).** The paragraph mixes counts with and without Holm adjustment ("Holm $p \le 0.021$" for 15 of 16, but "(6 significant)" and "(14 significant)" without saying which), and "the advantage over BOCPD and the scan thus comes from the clipped statistic and the restart design" is a causal reading of an observational comparison. *Fix:* label every count "Holm-adjusted" or "unadjusted", and write "is consistent with coming from".

**F13. Minor. Section 4.1 (line 172).** "At most 5.14% for i.i.d." is the largest simulated rate, with Monte Carlo interval 4.72% to 5.59%, not a bound. *Fix:* "the largest simulated rate was 5.14% (4.72% to 5.59%)".

**F14. Minor. Class labels missing.** The Nifty analysis (Section 7.4) is called pre-registered, but it is secondary to the two primary analyses of line 356; the Budget event analysis, the forecast regressions of line 373 and the real-data comparisons of E21 carry no confirmatory or exploratory label. *Fix:* label each as in line 356.

**F15. Minor. Section 7.2 (line 367).** The "pre-registered primary configuration" is run on July to December 2025, data the authors had already analysed with the plug-in. The pre-registration fixes the configuration, not the data's novelty. *Fix:* say so, and give the out-of-sample period's numbers next to it (they are in S7.6).

**F16. Minor. Section 7.2 (line 367) and S7.6.** "Matched up alarms came a median 4 to 10 minutes later" pools both periods; for July to December 2025 alone it is 4 to 9. In S7.6 the range "40.5 to 142" is the July to December 2025 figure, placed after a sentence about the out-of-sample period. *Fix:* name the period each time.

**F17. Minor. Section 3.1 (line 114) and S1.7.** The delay ratios of the restarted monitor (1.06 to 1.44, 1.02 to 1.14) are given without the plain detector's false alarms before the change in those runs (4.6% to 4.9% and 6.6% to 7.0%), so a reader cannot tell the comparison is not rate-matched. *Fix:* report them side by side.

**F18. Minor. Section 7.4 ("Two lessons").** "Alarm rates on real data, 0.66 to 0.71 per day on crypto with the plug-in, are far above the budget of 0.05 false alarms per day of Corollary 6" compares alarms per calendar day with a budget that Corollary 6 gives per monitored day. Since the plug-in monitors only about 30% of the minutes, the comparable figure is 2.0 to 2.4 alarms per monitored day. The conclusion is strengthened, not weakened. *Fix:* compare like with like.

**F19. Minor. Section 8, open problem (1).** The constant 0.98 (the floor $c_{\min} = 8.3$ of the companion note) is not defined in the paper. *Fix:* define it or drop the number.

**F20. Minor. Section 5.1 (line 237).** "The gap between 45 and 58 is the factor $n/B$ in the level" is not derived there or in S1.6. *Fix:* one line of derivation in S1.6.

---

## 6. The four open inconsistencies

**(a) 35.7% (day 10) and 50% in `docs/theory_per_period.md` line 9, against 28.0% (31.3%) in point m1 of the round-8 review.** Both are right; they refer to different day lengths and day laws. The quantity is the probability of a false alarm on day $d$ for the plain two-sided Shiryaev-Roberts detector with a known baseline, after $d - 1$ calm days with $Y_t = 1$; since no alarm can occur before day 40 on calm days, it equals the probability of a first alarm on day $d$. From `th6_calm_family.csv`, on day 10: with 390-step days, 28.02% (Gaussian day) and 31.29% (two-point day on $\lbrace 0, B\rbrace$); with 1440-step days, 34.30% (Gaussian) and 35.735% (two-point). On day 20 the two-point 1440-step value is 73.16%. The bound $dm/A$ is 50% on day 10 for both day lengths ($10 \times 390/7800 = 10 \times 1440/28{,}800 = 0.5$). The round-8 figure (28.0%) is the 390-step Gaussian day; 35.7% is the 1440-step two-point day. **Fix (Minor):** the S1.7 sentence "The largest single-day probabilities found, 35.7% ... and 73.2% ..." (line 675) and the per-period note should say "two-point day, 1440-step days".

**(b) Is $(2m + \kappa)/A$ correct for a window across one boundary?** Yes. Let the window of $m$ steps start at $s$ in day $i$ and end in day $i + 1$, with $b_i \le s < b_{i+1} < s + m$. Given $\mathcal F_{b_i}$, Theorem 8 bounds the probability of an alarm in day $i$ by $(m + \kappa)/A$, and the restart at $b_{i+1}$ bounds an alarm in the remaining part of the window by $m/A$ given $\mathcal F_{b_{i+1}}$, hence also given $\mathcal F_{b_i}$ by the tower property. The union bound gives $(2m + \kappa)/A$ given $\mathcal F_{b_i}$, and so unconditionally. It does not hold given $\mathcal F_s$ in general, since the statistic at $s$ may already be near the level; the paper does not claim that. A slightly sharper bound is $(\kappa + m + (s - b_i))/A \le (\kappa + 2m - 1)/A$, because the part of the window in day $i+1$ has length $s + m - b_{i+1} \le m - (b_{i+1} - s)$. **What is correct:** the bound as stated. **Fix (Minor):** state the conditioning ($\mathcal F_{b_i}$ or unconditional) next to the bound in line 114.

**(c) strucchange version: 1.6-0 in `docs/report_parts/E21.md` against 1.5-4 on CRAN.** I could not open CRAN or its mirrors (hosts not reachable from my environment), so I cannot verify this from the primary page. Search-engine listings show 1.6-0, dated 4 August 2026, as the current release and 1.5-4, dated 2 September 2024, as the previous one; listings that show 1.5-4 appear to be stale mirrors. On that evidence 1.6-0 is the current version and both figures were true at different times, but this is **unverified against CRAN itself**. The paper cites no version (line 902 and reference line 544). **Fix (Minor):** give the version and the date read in S5, since the boundary was taken from the package source.

**(d) E19 shares of minutes monitored.** The paper's "46% to 49% of the minutes, against 30% to 33%" (line 367) matches `e19_cs_cost.csv` for July to December 2025, band 0.25, five Binance coins: confidence-sequence baseline 0.462 to 0.489, plug-in 0.298 to 0.332. The E19 report's summary range "43% to 50% against 29% to 34%" is the out-of-sample period (0.429 to 0.505 against 0.287 to 0.336), which is also the union over both periods. Neither figure is wrong. **Fix (Minor):** the E19 summary should name its period; the paper's figure is correct for the period it names.

---

## 7. Presentation

* **Length and self-containment.** The main text is about 10,000 words and reads as a self-contained paper, with three exceptions: the rank detector (J3), the running-truncation restarts of Section 6.5 (J3), and the undefined 0.98 of open problem (1) (F19). The landmark shocks of Section 7.3 (four case studies) add little to the argument and could move to S7.5, which would leave room for the definitions of J3.
* **Nothing essential is buried**, except that the scope of each guarantee is spread across many paragraphs; the table of Section 3.2 above, placed at the end of Section 3, would help.
* **Cross-references.** Tables 1, 2 and S1 to S12 and Figures 1, S1 and S2 all exist, and every reference I followed points to the right place. Supplement sections S1 to S9 and S7.1 to S7.6 are referenced correctly.
* **Notation.** Still overloaded (m2 partly answered): $W$ (calibration window and pathwise product $W_n$), $\rho$ (Proposition S8 and Proposition S13), $\kappa$ (truncation moments $\kappa_2$ and restart slack), $C$ (Theorem 15 constants and the CUSUM-type statistic), $F$ (calibration distribution and filtrations in places), $h$ (tolerance horizon and other uses), $V$ (event in Theorem 16 part 4 and variance). Renaming four of these would remove most of the friction.
* **Abstract.** 216 words by the authors' count, 231 by mine; cut to 200 (m1).
* **Files.** Regenerate `paper/short_paper.pdf` (m10).

---

## 8. Major findings

**J1. Major. Abstract (line 9), Sections 7.1 to 7.4. The real-data evidence rests on the plug-in baseline.** The revision adds a pre-registered confidence-sequence configuration (Section 7.2) and permuted null streams (Section 7.1, S7.6), and these are well done. But the forecast regressions, the matched-competitor comparisons, the FOMC and CPI primary test, the Nifty analysis and the landmark shocks all use the plug-in, which has no finite-sample guarantee, and no real-data run uses a configuration with a per-day guarantee (the restarted monitor with the confidence-sequence baseline). The permuted-stream check is close to tautological for monitors with a known baseline, since permutation makes the null hold by construction (the paper says this, in fairness). The abstract's FOMC sentence then reads as confirmatory evidence for the guarantee. *Fix:* (i) rerun the FOMC and CPI primary test with the confidence-sequence baseline and the day restart, reporting it alongside the plug-in result, or say plainly in the abstract and Section 7.3 that the primary test used a monitor without a guarantee; (ii) give one sentence motivating the band $b$; (iii) reword the abstract's last sentence (F8).

**J2. Major. Related work (line 28), Section 3.3, Section 4.2, Section 6.2, Discussion (line 397). Novelty and positioning rest on papers that were not read.** By the authors' own notes, several papers whose contents the text states could not be opened: Ram and Ramdas (2026), compared with Theorem 25; Clerico (2025), compared with Theorem 12 and Lemma 20(a); Shekhar and Ramdas (2024), where the paper says "we add one-sided ranges, the sample split and the explicit tolerance"; Lai (1998), read only through secondary sources; Mei (2008); Horváth, Kokoszka and Zhang (2006); Chu, Stinchcombe and White (1996); Knoblauch, Jewson and Damoulas (2018); Bibinger, Jirak and Vetter (2017). I could not open them either, so I cannot confirm or refute these statements. The closing statement that no method was found with the same combination of properties, and "What is new and what is not" in the Discussion, depend on them. *Fix:* for each paper not read in full, either read it before resubmission or mark the comparison as based on the abstract or a secondary source; soften "we add" claims about unread papers to "as far as we can tell from its abstract". The editor may wish to ask a referee with access to check the Clerico and Shekhar-Ramdas comparisons in particular, since they bear directly on contributions 1 and 2.

**J3. Major. Abstract (line 9), Section 7.1 (line 359), Section 6.5 (line 346). Objects used in the main text are defined or proved only outside the paper.** The rank detector is one of the three "proved monitors" in the abstract and in the primary null-stream endpoint (3.38% against 5%), but it is defined only in the companion note, which defers its validity proof to the repository notes. The running-truncation restarts of Section 6.5 "keep a proof", again given only in the companion note or the repository notes, and the companion note itself defers its proofs to the repository. A journal reader cannot check either claim. *Fix:* add to the supplement a definition of the rank detector and a short validity proof (conformal p-values from ranks, which are exactly uniform under exchangeability, fed to the betting martingale, then Theorem 4), and either the proof for the running-truncation restarts or a statement that it is outside the paper.

**Critical findings: none.**

**Minor findings:** F1 to F20 above, the m22 contradiction in S1.1 (line 556), the strucchange version (6c), the S1.7 day-length label (6a), the E19 period label (6d), and the presentation items of Section 7.

---

## 9. The three most important changes

1. **J1:** run the primary FOMC and CPI test also with the confidence-sequence baseline and the day restart, or state plainly in the abstract and Section 7.3 that it used a monitor without a guarantee; reword the abstract's last sentence.
2. **J2:** read, or mark as unread, every paper whose contents are used for novelty or positioning, and soften the corresponding claims.
3. **J3:** define the rank detector and prove its validity in the supplement, and do the same for the running-truncation restarts or drop the claim that they "keep a proof".
