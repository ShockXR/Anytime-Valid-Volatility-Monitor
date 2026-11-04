# Referee report on "Anytime-valid detection of volatility changes in high-frequency returns with robust truncation"

> *Note added 3 November 2026: the "independent reviewer" or "referee" below was a separate AI agent (Claude, by Anthropic) instructed to act as a sceptical referee. It was not a human expert, and the review has not been confirmed by a human.*

**Manuscript.** `paper/short_paper.md`, draft of November 2026, as of commit `8f71a2e` on branch `claude/project-thread-ytr1da` (879 lines). Line numbers below refer to that file.

**Reviewer.** Anonymous referee, written as for *Sequential Analysis* or the Theory and Methods section of *JASA*. I did not write any part of the paper or of the repository.

**Scope.** The whole paper, including Appendices A to D and Supplement S1 to S6. As instructed, a robustness theorem for drift in the conditional mean (TH4) and a landmark-event case study (E18) are not yet in the paper and are not reviewed. I read the earlier internal reviews (`docs/review/`, rounds 1 to 5) first; points that they raised and that the current text has fixed are not repeated. None of them raised M2, M3, M5 or M7 below, or the per-day reading of Theorem B.

**What I ran.** Short scripts in `/tmp/claude-0/-home-user/e30ce406-b0e8-5b23-b612-d49eb104d0ee/scratchpad/referee/` (outside the repository): `constants.py` (truncated moments, $B$, speed limits and floors), `delay.py` (growth rates, $\mathrm{KL}_{\inf}$, the bound of Theorem 8(a), the losses of Table 1), `lamstar.py` (where the growth-optimal decrease bet lies), `later_days.py` and `later_days_cusum.py` (false alarm probability on later days under a null stream), `wordcount.py` (length by part). I also ran the repository's own `scripts/trace_numbers.py` on the paper: all 2224 numbers trace to a source file. Web searches for references were made on 3 November 2026; arxiv.org, ideas.repec.org and some publisher pages could not be opened from this machine, so for several references I rely on search-engine listings, and I say so.

---

## 1. Summary and recommendation

**Summary.** The paper proposes an online monitor for changes in the volatility of one-minute returns. Each return is clipped at $\tau = c\,\hat s$, where $\hat s$ is a scaled median absolute deviation estimated from a calibration window, and the clipped squared return $X_t = \min(r_t^2, \tau^2)$ is monitored by two Shiryaev-Roberts e-detectors built from linear bets, one against "the conditional mean of $X_t$ has risen above $(1+b)m_0$" and one against "it has fallen below $(1-b)m_0$". With a known baseline the authors prove $\mathbb{E}[T_A] \ge A$ and $\mathbb{P}(T_A \le m) \le m/A$ over a nonparametric class defined only by the conditional mean of $X_t$. A one-sided Bahadur-Savage argument shows that no nontrivial e-process exists for "the variance has not fallen", which motivates the bounded statistic, and the bound $B = \tau^2/m_0$ then caps the evidence per step for decreases. A confidence-sequence baseline gives $\delta + m/A$. The delay theory comprises a Wald-type bound, first-order optimality of refined bet grids (credited to Ram and Ramdas 2026), an "information ladder" that splits the delay cost of robustness into truncation, nonparametric and grid losses, and a closed-form bound on the signal left after GARCH standardisation. Simulations compare the monitor with Gaussian likelihood-ratio, CUSUM, repeated-test, BOCPD and a reconstructed Bibinger-Jirak-Vetter (BJV) test, and include an adversarial search over null streams. Real-data studies cover five cryptocurrencies, the Nifty 50 and macro announcements. A final section raises questions about the estimator of Bapat (2026).

**Recommendation: major revision.** For *JASA* Theory and Methods, I would not support acceptance without substantially more theory (M3). The paper is careful, unusually honest about its limitations, and reproducible, and its proofs are correct as stated: I found no false theorem. It has four problems that a strong referee would not let pass:

1. The anytime-valid guarantee, which is the paper's reason to exist, is not exercised by any of the real-data evidence. All main real-data runs use the plug-in baseline, which the paper itself shows has no proof and exceeds the bound over multi-day horizons. Moreover, on real data the null is false on most days (M1).
2. The paper repeatedly reads Theorem B as "at most 5% per day". The theorem bounds the first $m$ steps after (re)calibration only. Under a null stream in the paper's own class, I find a 28% chance of a false alarm on day 10 (M2).
3. Most of the theory is a known result or a direct instance of one, which the paper largely concedes. The most interesting new statements, such as the speed limit and the impossibility, are proved only for e-detectors, while an equally short argument covers every procedure (M3).
4. The confidence-sequence guarantee holds under a null that depends on the algorithm's own bets and is much narrower than the abstract suggests (M4).

The benchmarks are partly strawmen (M5), the real-data inference needs multiplicity control and a clear primary/secondary split (M6), important literature is missing (M7), and at about 28,000 words by my count the manuscript is not a short paper (M8). I recommend major revision for *Sequential Analysis*. For *JASA* T&M, a substantially strengthened theoretical core would be needed (see M3 for concrete options). Alternatively, the work could be retargeted, in a much shorter form, to an applied or financial-econometrics venue.

**Strengths worth keeping.** Every number is traceable. Negative results are reported (the plug-in failure, BOCPD's speed on clean data, ratio rules matching the precision on real data). The adversarial search of E16 and the slack decomposition are a model of how to probe a conservative bound. The information ladder (Corollary 6, Table 1) is a clear and useful way to explain the price of robustness.

---

## 2. Novelty

### 2.1 Contribution by contribution (introduction, lines 27 to 32)

| Claimed contribution | Verdict | Reasoning |
|---|---|---|
| Robust statistic and conditional-mean nulls, Lemma 1 (line 27) | Straightforward instance of known ideas | Clipping at a robust scale is Huber winsorisation; truncated realised variance is Mancini (2009); bounded influence in sequential tests goes back to Huber's (1965) censored probability ratio test (not cited). A conditional-mean null for a bounded statistic is the standard setting of betting e-processes. Lemma 1 is a one-line monotonicity argument. What is new is the specific combination for one-minute returns with a predictable MAD threshold. That is an application contribution, not a methodological one. |
| Validity: Theorems A, B, Proposition 1, Corollary C (line 28) | Known | Theorem A is the e-detector guarantee of Shin, Ramdas and Rinaldo (2024), as the paper says. Theorem B is Ville's inequality applied to $M_t + (m - t)$, an elementary argument that the paper also concedes is not new. Proposition 1 is a Bonferroni split. Corollary C is a supermartingale bookkeeping argument. None is a contribution on its own. |
| Bounded statistic: Proposition 2, Corollary 3, Corollary 7 (line 29) | Known in substance; the interpretation is new | The impossibility of testing "mean at least $m$" for nonnegative unbounded variables is the one-sided Bahadur-Savage phenomenon. In the betting literature the asymmetry is folklore: lower confidence sequences exist for nonnegative heavy-tailed means (Mineiro 2022) but upper ones do not. The paper says it does not claim novelty. Lemma 5(a) (every e-variable for a bounded mean is dominated by a bet) is in the spirit of the optimality-of-coin-betting results of Clerico (2025), not cited. The speed limit $\log(B/(B-1))$ follows directly from the bounded range. The application-specific reading (a speed limit for decreases set by $c$) is useful, and could become a genuine contribution if extended to all procedures (M3). |
| Estimated baseline: Theorem D (line 30) | Modest | The one-sided-range empirical Bernstein bound is a direct use of the inequality of Fan, Grama and Liu (2015) within the framework of Howard, Ramdas, McAuliffe and Sekhon (2020, 2021). Combining a confidence sequence for the pre-change parameter with a detector via a union bound is close in spirit to Shekhar and Ramdas (2024), who reduce change detection to sequential estimation. That paper is in the authors' own `docs/literature_review.md` but is not cited in the paper. The negative result for the plug-in (Proposition A.6 of the theory note) is new, but it is elementary and set in a regime (fixed $n$, $A \to \infty$) of limited practical relevance (M3). |
| Detection delay: Theorems 8, 9, Corollary 6, Proposition 10 (line 31) | Mostly known; ladder and cap modestly new | Theorem 8(a) is the Wald and Lorden (1970) argument applied to one bet; Theorem 8(b) is a monotone coupling. Theorem 9 is, as the paper says, an instance of Ram and Ramdas (2026). Corollary 6 combines data processing with the $K_{\inf}$ duality of Honda and Takemura (2010): standard tools, though the "ladder" framing is helpful. Proposition 10 is two lines of algebra on the GARCH recursion. Its message, that a GARCH filter absorbs level shifts, is qualitatively known (Diebold 1986; Lamoureux and Lastrapes 1990, both cited). The closed form is new as far as I know, but it bounds the ratio of means, not the statistic the monitor uses. |
| Evidence (line 32) | The most original part | Extensive, honest and reproducible. Its main real-data message is that the e-detector is not more precise than simple or published rules (lines 9, 419, 441). |

**Overall.** The paper's own summary (line 475) is accurate: the tools are known and the contribution is their combination for this application. That is a legitimate contribution, but for a theory-and-methods venue the paper needs either a sharper theoretical result or a stronger demonstration that the guarantee matters in practice (M1, M3).

### 2.2 Cited references I checked

| Reference as cited | Status | How checked |
|---|---|---|
| Shin, Ramdas and Rinaldo (2024), NEJSDS 2(2), 229-260 | Verified | Web search listing of the NEJSDS article page (vol. 2, issue 2, pp. 229-260, doi 10.51387/23-NEJSDS51) |
| Ram and Ramdas (2026), arXiv:2602.05272 | Exists; content not opened | Web search listing of arxiv.org/abs/2602.05272v1 and a secondary summary (posted 5 February 2026; universal lower bound attained for bounded means, i.i.d. data). I could not open arXiv to check that the paper's Theorem 9 is an instance of it, as claimed on line 233. |
| Horváth, Trapani and Wang (2024), arXiv:2404.17885 | Apparently published under another title | A web search returns a *Journal of Time Series Analysis* 46(5), 981-996 (2025) article, "Sequential monitoring for changes in GARCH(1,1) models without assuming stationarity", doi 10.1111/jtsa.12824, by the same authors, as the published version. I could not open the RePEc page to confirm. Please check and cite the published version. |
| Andersen, Todorov and Zhou (2025), *Quantitative Economics* 16(2), 459-495 | Verified | Econometric Society publication page and RePEc listing in search results |
| Larsson, Ramdas and Ruf (2024), arXiv:2402.18810 | Exists; venue not verified | arXiv listing (v4, February 2025) and an LSE working-paper listing; I found no journal version |
| Howard, Ramdas, McAuliffe and Sekhon (2021), *Annals of Statistics* 49(2) | Verified | Project Euclid link in search results |
| Waudby-Smith and Ramdas (2024), JRSSB 86(1) | Verified (paper); theorem number not verified | Search listing of the JRSSB article. I could not open the text to confirm that the predictable plug-in empirical Bernstein bound is "Theorem 2" (line 148). |
| Bapat (2026), *IREF* 109, 105482 | Not verifiable online by me | No search hit. `docs/paper_review.md` states that the details were checked against the PDF and that the volume and DOI could not be confirmed online. Before submission, confirm them against the publisher. |

### 2.3 Important literature the paper misses

Each entry below was found by a web search on 3 November 2026. I list what I could and could not confirm. All are relevant enough that a referee in the area would expect them to be discussed.

**Sequential monitoring of variance and GARCH parameters in econometrics.** This is a whole literature on online monitoring with a training sample, with asymptotic control of the probability of any false alarm over an open or closed horizon. It is the natural comparator for the "estimated baseline" problem of Section 3.4 and for Section 4.5.

* Chu, C.-S. J., Stinchcombe, M. and White, H. (1996). Monitoring structural change. *Econometrica* 64(5), 1045-1065. Verified via the Econometric Society publication page in search results. Its abstract makes the paper's "repeated looks" point (line 21) for structural change.
* Berkes, I., Gombay, E., Horváth, L. and Kokoszka, P. (2004). Sequential change-point detection in GARCH(p,q) models. *Econometric Theory* 20(6), 1140-1167. Verified via RePEc and Cambridge Core listings in search results.
* Horváth, L., Kokoszka, P. and Zhang, A. (2006). Monitoring constancy of variance in conditionally heteroskedastic time series. *Econometric Theory* 22(3), 373-402. Verified via Cambridge Core and RePEc listings in search results. This is the closest prior work on online variance monitoring and must be discussed.
* Dette, H. and Gösmann, J. (2020). A likelihood ratio approach to sequential change point detection for a general class of parameters. *JASA* 115(531), 1361-1377, doi 10.1080/01621459.2019.1630562. Verified via RePEc and publisher figshare listings in search results. A *JASA* T&M paper on monitoring with a training sample, including variances.

**Robust change detection.**

* Huber, P. J. (1965). A robust version of the probability ratio test. *Annals of Mathematical Statistics* 36(6), 1753-1758, doi 10.1214/aoms/1177699803. Verified via Project Euclid listing. Censoring likelihood-ratio factors is the direct ancestor of clipping the statistic so that one observation has bounded effect on the wealth (line 131).
* Unnikrishnan, J., Veeravalli, V. V. and Meyn, S. Minimax robust quickest change detection. arXiv:0911.2551. The arXiv version is verified via search; I could not confirm the journal version (often cited as *IEEE Transactions on Information Theory*, 2011), so the authors should check it. It is the robust (least favourable distribution) answer to the same problem.
* Fearnhead, P. and Rigaill, G. (2019). Changepoint detection in the presence of outliers. *JASA* 114(525), 169-183, doi 10.1080/01621459.2017.1385466. Verified via RePEc and arXiv listings. Its central argument is that only bounded losses are robust to arbitrary outliers, a close cousin of the paper's argument for a bounded statistic.
* Knoblauch, J., Jewson, J. and Damoulas, T. (2018). Doubly robust Bayesian inference for non-stationary streaming data with β-divergences. *NeurIPS 31*. Verified via the NeurIPS proceedings page. A robust BOCPD, the obvious fairer BOCPD competitor (M5).

**E-values, confidence sequences and change detection.**

* Shekhar, S. and Ramdas, A. (2024). Reducing sequential change detection to sequential estimation. *ICML*, PMLR 235, 44628-44642. Verified via PMLR listing. Directly relevant to Theorem D.
* Vovk, V., Petej, I., Nouretdinov, I., Ahlberg, E., Carlsson, L. and Gammerman, A. (2021). Retrain or not retrain: conformal test martingales for change-point detection. PMLR 152, 191-210. Verified via PMLR and Royal Holloway listings. Conformal CUSUM and Shiryaev-Roberts procedures, directly relevant to the rank detector of E14 (Section 5.5).
* Dandapanthula, S. and Ramdas, A. (2025). Multiple testing in multi-stream sequential change detection. arXiv:2501.04130. Verified via arXiv listing in search results. Relevant to the panel monitor (Section 5.5, Supplement S4.2), which it does not cite.
* Wang, H. and Ramdas, A. (2023). Huber-robust confidence sequences. AISTATS, PMLR 206, 9662-9679. Verified via PMLR listing. Anytime-valid inference under heavy tails and contamination.
* Mineiro, P. (2022). A lower confidence sequence for the changing mean of non-negative right heavy-tailed observations with bounded mean. arXiv:2210.11133. Verified via arXiv listing. The constructive side of the asymmetry in Section 3.3.
* Clerico, E. On the optimality of coin-betting for mean estimation. arXiv:2412.02640. Verified via arXiv listing; a search summary gives *International Journal of Approximate Reasoning* 187 (2025), 109550, which I could not open. Relevant to Lemma 5(a).

**E-values in finance.** Both are in the authors' `docs/literature_review.md` but not in the paper.

* Wang, Q., Wang, R. and Ziegel, J. E-backtesting. *Management Science*, doi 10.1287/mnsc.2023.01659. Verified via the INFORMS DOI listing in search results; volume and pages not confirmed.
* Choe, Y. J. and Ramdas, A. Combining evidence across filtrations. arXiv:2402.09698. Verified via arXiv listing. Its Section 4.2 tests the randomness of high-volatility days in financial returns, and it is the general setting of the multi-frequency monitor (E17).

**False alarm criteria beyond the ARL.** These matter for M2.

* Lai, T. L. (1998), already cited, introduced the window-limited false alarm probability $\sup_k \mathbb{P}(k \le T < k + m)$ (verified via the IEEE listing and secondary descriptions in search results). The paper cites Lai only for the delay lower bound.
* Mei, Y. (2008). Is average run length to false alarm always an informative criterion? *Sequential Analysis* 27(4). Verified via the Georgia Tech prize announcement and discussion articles in search results; the page range and DOI were not confirmed.
* Tartakovsky, A., Nikiforov, I. and Basseville, M. (2014). *Sequential Analysis: Hypothesis Testing and Changepoint Detection*. Chapman and Hall/CRC. Verified via publisher and library listings. The standard reference for Shiryaev-Roberts procedures and their operating characteristics.

**The novelty statement on line 36.** "Our literature search ... did not find a method that combines jump-robust high-frequency volatility measures with anytime-valid sequential change detection." My own search found none either. The statement should nevertheless be followed by a discussion of the nearest neighbours above: econometric variance monitoring, robust quickest detection, and e-processes applied to volatility and to risk backtesting.

---

## 3. Correctness

### 3.1 Every stated result against its proof

| Result | Where proved | Verdict | Remarks |
|---|---|---|---|
| Lemma 1 | line 77 | Correct | "Falls exactly when $\sigma$ falls" (line 131) needs strict monotonicity of $g_\tau$, which holds iff $\mathbb{P}(0 < \lvert\varepsilon\rvert < \tau/s) > 0$. Rounded prices with a fixed tick are not of the form $\sigma_t\varepsilon_t$, so Lemma 1 does not apply to the tick laws used in Sections 4.4 and 5 (m15). |
| Lemma 2 | line 95 | Correct | For missing returns (line 63) the filtration must contain the missingness indicator, and the condition is $\mathbb{E}[X_t \mid \mathcal{F}_{t-1}, \text{observed}] \le m$; state it as an assumption (m14). |
| Theorem A | line 99 | Correct | |
| Theorem B | line 105 | Correct as stated | Its reading as "per day" on lines 9, 107 and 396 is not (M2). |
| Proposition 1 | line 111 | Correct | The rule $\max(M^\uparrow, M^\downarrow) \ge 2A$ is dominated by $\bar M \ge A$ with the same guarantee (m6). |
| Corollary C | line 115 | Correct (sketch adequate) | With $b = 0$ the hypothesis requires each segment's baseline to equal the true conditional mean, which no estimated baseline satisfies. The text says so (line 117). |
| Proposition 2 | line 123 | Correct | |
| Corollary 3 | line 129 | Correct | Can be strengthened to every stopping rule with a finite-horizon guarantee by the same construction (M3). |
| Lemma 4 | Appendix B | Correct (standard) | |
| Lemma 5 | Appendix B | Correct | (b) needs an interior maximiser. For decreases the maximiser is interior only barely: $\lambda^\ast/\lambda_{\max}$ = 0.9917, 0.99987 and $1 - 2.4\times 10^{-8}$ for sd ratios 0.8, 0.67 and 0.5 (`lamstar.py`), because the atom of $Y$ at $B$ has probability $1.8\times10^{-4}$, $7.5\times10^{-6}$ and $2.0\times10^{-9}$ (m19). |
| Corollary 6 | Appendix B | Correct | |
| Corollary 7 | Appendix B | Correct; numbers verified | The floor of the default grid for the two-sided monitor, at $2A$ per side, is 79 steps, not the 73 quoted at level $A$ (m20). |
| Theorem 8 | line 215 | Correct | The bound is for one side at level $A$; the monitor runs each side at $2A$ (m18). The overshoot term uses Lorden's (1970) bound correctly. |
| Theorem 9 | line 231 and Appendix B | Correct | (a) The worst-case delay is evaluated with pre-change law $P^\ast$, not the true $P_0$; say so (m16). (c) The constant $C = \bar v^2/(1 - b'\lvert v_{\mathrm{lo}}\rvert)^2$ explodes for decreases because $\lambda^\ast$ is almost at $\lambda_{\max}$ (m19). |
| Proposition 10 | line 267 | Correct; all six values of $\rho(k)$ in Table 2 verified | The "cap" bounds the ratio of means. The simulated mean of the ratio exceeds it at persistence 0.95, $k = 1000$: 2.204 against 2.0 (`results/e10d_signal_cap.csv`, visible in Figure 1). The paper says it has no proof for the mean of the ratio, but the abstract and line 31 still call it a cap on the visible signal (Section 4, items 3 and 6). |
| Lemma A.1, Theorem D.1 | Appendix A | Correct | |
| Theorem D | Appendix A | Correct as a conditional statement | Hypotheses $(H^\uparrow)$, $(H^\downarrow)$ depend on the algorithm's bets and on $T_0$ (M4). The proof's Step 2 uses them only on $G$, so the theorem could be stated with the algorithm-free sufficient condition as its hypothesis. |
| Corollaries E, F | Appendix A | Correct | |
| Proposition A.6, Lemma A.7 (cited on line 139) | only in `docs/theory_running_tau.md` | Correct as stated in the note | Results used in the main text must be stated and proved in the paper or its supplement (M8). |
| Supplement S4.3 constants | | Verified | 17.7, the ratios 1.27, 6.33 and 19, and $A$ = 12,160 and 4114 |
| Supplement S6 floor | line 877 | Verified | $c = 8.325$ and $P \le 4.8985\%$, 0.9797 of the bound. This $c$ clashes with the truncation constant $c$ (m2). |
| Second-order formula (S1) | heuristic | The digamma form of $\zeta$ follows from Dufresne's identity under the Brownian approximation | For a single bet, the conditions of nonlinear renewal theory (non-arithmetic bounded increments, a slowly changing $\log V_n$) look routine to check (m24). |

### 3.2 Issues of substance found in checking

1. **Theorem B is a statement about the first $m$ steps only (M2).** $\mathbb{E}[M_t] \le t$ under the null, so at the start of day $d$ the statistic carries up to $(d-1)m$ in expectation. Ville's inequality then gives only $\mathbb{P}(\text{alarm on day } d) \le dm/A$. I simulated a null stream inside the class of Theorem B: $Y_t \equiv 1$ (that is, $\lvert r_t\rvert = \sqrt{m_0}$, a legal law) for $d - 1$ days, then i.i.d. clipped Gaussian returns with the same conditional mean, with the default two-sided monitor, $c = 3$, $A = 7800$ and 40,000 paths (`later_days.py`). The probability of a false alarm on day $d$ was 0.97% ($d = 1$), 3.2% (2), 11.5% (5), 28.0% (10), 46.2% (15) and 61.5% (19), against a nominal 5% per day. With the i.i.d. law $\lbrace 0, B\rbrace$ on day $d$ it was 1.8%, 4.8%, 14.1%, 31.4%, 49.6% and 66.8%. Even under natural i.i.d. Gaussian returns, Table 4 implies a conditional per-day hazard that rises from 1.1% on day 1 to about 3.1% (days 2 to 5), 4.2% (days 6 to 10) and 4.6% (days 11 to 20). The CUSUM-type statistic does not accumulate in this way: 0.04% on day 1 and 0.07% on day 10 on the same stream (`later_days_cusum.py`). The paper's per-period statements must be restated (M2).
2. **The impossibility and the speed limit hold for every procedure, not only for e-detectors.** Take any stopping rule $T$ with $\mathbb{P}(T \le n) \le n/A$ for all $n$ and every law in the decrease null class on $[0, B]$, and a path $y$ with $T(y) = n$. Consider the law that at each step with $y_s < 1$ draws $y_s$ with probability $(B-1)/(B-y_s)$ and $B$ otherwise, and draws $y_s$ with probability 1 when $y_s \ge 1$. It is in the class, so $\prod_{s \le n} \max\lbrace 1, (B - y_s)/(B-1)\rbrace \ge A/n$. The left side is the wealth of the bet at $\lambda_{\max}$ with losing factors replaced by 1. For $B = \infty$ (unclipped), the factor is 1 and $T \ge A$ on every path, which is Corollary 3 for all procedures. For complete collapse with $c = 3$ and $A = 7800$ the bound gives $n \ge 45$ (`constants.py`; 23, 73 and 106 for $c$ = 2, 4 and 5), against 57.7 for any Shiryaev-Roberts e-detector with linear bets and 73 for the default grid. The increase side is symmetric: $\prod_s \max(1, y_s) \ge A/n$. This short argument fills the gap stated on line 201 ("For rules that are not e-detectors the only lower bound we have is the asymptotic one") and would make Section 3.3 a genuine contribution (M3).
3. **The decrease side lives at the boundary.** For the decreases of Table 1 the growth-optimal bet is within 0.9% (sd 0.8), 0.013% (sd 0.67) and $2.4\times 10^{-6}$% (sd 0.5) of $\lambda_{\max}$. The asymptotic theory of Section 4.2, built on an interior $\lambda^\ast$ and a grid on $[a, b']$ with $b' < \lambda_{\max}$, is therefore uninformative for decreases: its constants blow up. The relevant object is the bet at (or just below) $\lambda_{\max}$, whose growth $\mathbb{E}_Q \log((B - Y)/(B-1))$ is finite iff $Q$ has no atom at $B$. This also explains the "grid loss" of 1.02 to 1.05 for decreases: the top decrease bet is $0.95\lambda_{\max}$.
4. **The null of Theorem D (M4).** With $b = 0$, the algorithm-free sufficient condition (line 162) requires the conditional truncated variance during monitoring never to exceed its minimum over the calibration window. Any volatility clustering, however mild, then violates it. In effect the confidence-sequence guarantee covers Corollary E (a constant conditional law) plus whatever slack $b$ buys. That is far narrower than the "every form of dependence" of the abstract, which applies only to a known baseline.

### 3.3 Constants recomputed

| Quantity | Paper | Recomputed | Script |
|---|---|---|---|
| $\kappa_2(3)$, $\kappa_4(3)$ | 0.995, 0.964 | 0.9950, 0.9638 (Monte Carlo with $4\times10^6$ draws: 0.9944, 0.9638) | `constants.py` |
| Share clipped, $a = 3$ and $a = 2.023$ | 0.27%, 4.3% | 0.270%, 4.302% | `constants.py` |
| $\kappa_2$, $\kappa_4$ at $a = 2.023$ | 0.925, 0.705 | 0.9247, 0.7045 | `constants.py` |
| $B = c^2/\kappa_2(c)$, $c = 3$ | 9.045 | 9.0452 | `constants.py` |
| $\log(B/(B-1))$ | 0.117 | 0.1172 | `constants.py` |
| Floor for any linear-bet e-detector, $c = 3$, at $A$ and $2A$ | 57.7, 63.6 | 57.71, 63.62 | `constants.py` |
| Same floor, $c$ = 2, 4, 5 | 29, 96, 141 | 28.65, 95.94, 140.76 | `constants.py` |
| Default decrease grid floor, $c$ = 2 to 5, at $A$ | 36, 73, 122, 180 | 36, 73, 122, 180; at $2A$: 39, 79, 134, 199 | `constants.py` |
| Floor for any procedure with $\mathbb{P}(T \le n) \le n/A$ (new) | not in paper | 23, 45, 73, 106 | `constants.py` |
| S6 floor and bound | 8.3, 4.90%, 0.98 | 8.325, 4.8985%, 0.9797 | `constants.py` |
| Table S1, seven Gaussian rows: $g_{\mathrm{grid}}$, $\log A/g$, Theorem 8(a) | e.g. 0.0502, 179, 233 | all match; e.g. 0.0502, 178.6, 233.4; 0.0839, 106.8, 127.2 | `delay.py` |
| Table 1, all 18 loss entries | e.g. 1.03, 1.10, 1.02, 1.16 | all match to rounding (sd 1.5 total 1.365, printed 1.37) | `delay.py` |
| Table 2, $\rho(k)$ | six values | exact | by hand |
| Table 2, $\ln A/I$ at persistence 0.999 | $I$ = 0.0003, 21 days | $I = 0.0003$ gives 23.8 days; the CSV value $I = 0.000333$ gives 21.4 | `results/e10d_signal_cap.csv` |
| $A$ for $\delta = 0.01$, $\alpha = 5\%$ | 9750, 36,000 | 9750, 36,000 | by hand |
| Table 3 Monte Carlo standard errors | 0.8 and 1.3 points | 0.81, 1.26 | by hand |

---

## 4. Claims versus evidence

Each entry gives the line, the sentence, the problem and replacement text.

1. **Line 9 (abstract):** "(1.1% to 1.4% false alarms per day against a bound of 5%)". *Problem:* Table 4 measures the first 390 steps from $M_0 = 0$; later days are not bounded by 5% (M2). *Replace:* "(1.1% to 1.4% of paths with a false alarm within the first 390-step day, against a bound of 5%)".
2. **Line 9:** "while Gaussian likelihood-ratio, CUSUM and Bayesian online change-point detectors raise false alarms on half or more of the days under jumps or fat tails". *Problem:* the thresholds of these competitors were calibrated on Gaussian data (lines 294, 385). *Replace:* "... when their thresholds are calibrated on Gaussian data, raise false alarms ...".
3. **Line 9:** "and a closed-form cap on the signal that remains after dividing returns by a GARCH forecast". *Problem:* Proposition 10 bounds the ratio of means; the mean of the ratio, which the monitor sees, exceeds the bound at persistence 0.95 (2.20 against 2.0). *Replace:* "and a closed-form bound on the post-break ratio of the mean squared return to the mean GARCH forecast".
4. **Line 9:** "On five cryptocurrencies the alarms carry information beyond EWMA and HAR volatility forecasts". *Problem:* the evidence is a hit rate against two slow forecasts (EWMA with a one-day half-life; HAR on 1-, 5- and 22-day averages), which a trailing ratio rule matches (Table 7). *Replace:* "On five cryptocurrencies, 71% to 81% of alarms were followed by a next day on the flagged side of an EWMA or HAR forecast, against about 50% at random times, as often as for a simple ratio rule tuned with hindsight".
5. **Line 9:** "they follow scheduled events far more often than control times do". *Problem:* Reserve Bank of India policy days did not differ ($p = 0.36$), and the primary CPI test on BTC was not significant ($p = 0.11$). *Replace:* "they follow Union Budget speeches and FOMC statements far more often than matched control times; Reserve Bank of India policy days and, in the pre-specified design, CPI releases on BTC showed no significant difference".
6. **Line 31:** "monitoring returns divided by a correct volatility forecast makes the null hold again, but caps the visible signal of a variance break at $(1-\beta)/\alpha$". *Replace:* "... but bounds the post-break ratio of the mean squared return to the mean forecast by $(1-\beta)/\alpha$ (Proposition 10); the mean of the ratio, which the monitor sees, is larger and is not covered by the proof".
7. **Line 107:** "bounds the probability of any false alarm per day by 5%". *Replace:* "bounds by 5% the probability of a false alarm within the first day after (re)calibration. Over longer runs Corollary C bounds the expected number of false alarms by 0.05 per day, but the probability of a false alarm on a given later day is not bounded by 5%".
8. **Line 396:** "so that with a known baseline the probability of any false alarm is at most 5% per day or session". *Replace:* "so that with a known baseline the expected number of false alarms is at most 0.05 per day or session (Corollary C), and the probability of a false alarm on the first monitored day is at most 5% (Theorem B)".
9. **Line 396:** "In all runs the settings were fixed before the data were examined". *Problem:* several analyses were added after first results, as flagged elsewhere: the fresh-calibration FOMC design (line 449), the long-window competitors (line 439), the restarted-monitor rerun of the event study (line 379), $\gamma = 2$ in E17 (line 381) and the rate-matched E17 comparison. *Replace:* "The e-detector's settings were fixed before these data were examined ...; analyses added after first results are marked as secondary where they appear and carry no confirmatory weight."
10. **Line 201:** "For rules that are not e-detectors the only lower bound we have is the asymptotic one of Theorem 9(a)." *Problem:* a non-asymptotic bound for all procedures follows in a few lines (Section 3.2, item 2). *Replace* with that bound: "Any rule with $\mathbb{P}(T \le n) \le n/A$ over the null cannot alarm on a path before $\prod_{s \le n}\max\lbrace 1, (B - y_s)/(B-1)\rbrace \ge A/n$; for complete collapse with $c = 3$ and $A = 7800$, not before 45 steps."
11. **Line 253:** "The information ratio explains the 10% to 70% delay penalty measured in E4". *Problem:* the prediction is against the oracle likelihood ratio and the simulation against a Gaussian likelihood-ratio mixture; they differ by up to 17% (1.27 against 1.48 at sd 0.8). *Replace:* "The information ratio is consistent with the 10% to 70% delay penalty measured in E4 (to within about 15%)".
12. **Line 286:** "At persistence 0.99 no detector with this false alarm budget can be expected to see a doubling of the long-run variance much faster than within about two days, and at persistence 0.999 it is practically invisible for weeks." *Problem:* the paper's own monitors had median delays of 1523 to 1791 steps (1.06 to 1.24 days) at persistence 0.99 (line 371), and Lai's bound is asymptotic. *Replace:* "Lai's asymptotic bound puts the delay at about $\ln A/I$: 1.9 days at persistence 0.99 and about three weeks at 0.999. At $A = 28{,}800$ this is a guide to orders of magnitude, not a bound; the monitors of Section 5.4 detected a doubling at persistence 0.99 with median delays of 1.06 to 1.24 days."
13. **Line 286:** "This is a property of the question, not of the monitor". *Replace:* "Under a Gaussian GARCH model with known parameters, this is a property of the detection problem, not of the monitor".
14. **Line 347:** "For decreases the Gaussian likelihood-ratio detector was faster and close to its nominal false alarm rate in these designs". *Problem:* over the 500 pre-change steps, 6% to 9% of Gaussian paths had a false alarm, against $500/7800 = 6.4\%$ for a known baseline; 9% is about two standard errors above. *Replace:* "... was faster, with false alarm rates before the change of 3% to 9% over 500 steps, at or somewhat above the 6.4% that Theorem B gives with a known baseline".
15. **Lines 379 and 449:** "the BTC CPI test, not significant there, gave $p = 0.003$ and $0.0003$". *Problem:* a different monitor, run after the primary test failed, cannot rescue it. *Replace:* "(exploratory, outside the pre-specified family) the BTC CPI comparison gave nominal $p = 0.003$ and $0.0003$".
16. **Line 389:** "the e-detector's advantages are that it needs no window and no change location, and that it runs online". *Problem:* it needs a change-free calibration window of 390 to 5000 steps (Table 6). *Replace:* "the e-detector needs no comparison window and no change location, only a change-free calibration stretch, and it runs online".
17. **Line 419:** "Only at the strictest band (0.5), out of sample and against the EWMA forecast, was the e-detector better, by about 4 points for increases ($p = 0.03$ to 0.06) and 10 points for decreases against the 60-minute rule ($p < 0.001$)". *Problem:* one subgroup out of 3 bands × 2 periods × 2 forecasts × 2 ratio rules × 2 directions (48 pooled cells; `results/e06e_competitor.csv` has bands 0, 0.25 and 0.5 and windows of 60 and 240 minutes), with unadjusted p-values that treat correlated alarms as independent. *Replace:* "In one of the 48 cells examined (band 0.5, out of sample, EWMA forecast) the e-detector was better by about 4 points for increases and 10 points for decreases; we report this as exploratory, since the p-values are unadjusted for the 48 comparisons and assume independent alarms".
18. **Line 441:** "No competitor reached the e-detector's up alarms: the best longer-window competitor, picked after seeing the results, was about 6 points lower". *Problem:* the difference is borderline ($p$ from 0.028 to 0.21 depending on the test), and the comparator was chosen after seeing the results. *Replace:* "No competitor's up alarms were confirmed significantly more often. The best longer-window competitor, chosen after seeing the results, was about 6 points lower; depending on how the dependence between alarms is modelled, $p$ ranges from 0.03 to 0.21, so we regard the difference as inconclusive."
19. **Line 27:** "so the nulls have a volatility meaning without a model for the return distribution". *Replace:* "so the nulls have a volatility meaning without a parametric model for the innovation law (Lemma 1 assumes $r_t = \sigma_t\varepsilon_t$ with i.i.d. innovations independent of the past)".
20. **Line 457:** "Its value on real data is the proof, not better alarms." *Problem:* the main real-data runs (Tables 7 and 8, Sections 6.2 to 6.4) do not use the confidence-sequence baseline, and the null is false on most days. *Replace:* "Its value would be the proof, not better alarms; but the main real-data analyses above use the plug-in baseline, and on most real days the null is false (Section 6.4)."
21. **Line 9, last sentence:** "on real data the contribution is settings fixed in advance and, with the confidence-sequence baseline, a proved false alarm bound, not extra precision." *Problem:* as in item 20 (M1). *Replace* after the revision of M1, or qualify: "... a proved false alarm bound, which our main real-data runs, made with a plug-in baseline, do not carry".

### 4.1 Multiple comparisons, post hoc analyses and p-values that assume independence

* **Table 7 and line 419.** Wilson intervals and two-proportion tests on pooled alarms. Alarms on the same calendar day across coins, and up and down alarms within a coin, are dependent (the paper notes this on line 417). There are 48 cells, and the only significant ones form a post hoc subgroup.
* **Table 8 and line 441.** Twelve competitor rows × 2 directions × 2 periods, with "two-proportion tests, $p \le 0.02$" and no multiplicity adjustment. The long-window rows were added after the first run, and the "best competitor" was picked after seeing the results. The same comparison is reported with three p-values (0.065, 0.047 and 0.21), and the text does not say which test is primary. That the week-clustered p-value (0.047) is smaller than the naive one (0.065) suggests a paired test; please define each test.
* **Line 445 and Supplement S5.1.** Fisher tests for five event types (Budget, RBI, expiry, election, COVID), with no correction and no statement of which were pre-specified. The Budget result ($p = 0.0001$) survives any correction; the expiry result ($p = 0.002$) is, by the authors' own reading, possibly a profile artefact.
* **Line 449 and Supplement S5.2.** The primary family is well specified (Holm). The fresh-calibration design, the restarted monitors (E14, "every awake cell had a stratified $p \le 0.006$") and the replications on ETH and SOL are secondary and should be reported without confirmatory language.
* **Supplement S4.3.** $p = 0.88$, 0.43 and 0.0014 treat alarms as independent; the text says so. Fine, but drop the p-values or use a block bootstrap.
* **Supplement S5.3.** "No significant difference in any of 12 cells" is a power statement; give intervals instead.

**Fix.** For each real data set, pre-specify one primary endpoint, or state that none was pre-specified. Adjust within the family (Holm). Use a block bootstrap over calendar days across all coins jointly for intervals. Report everything else as exploratory, with intervals rather than p-values.

---

## 5. Presentation

**Length.** By my count, with display math removed and inline math counted as one word (`wordcount.py`): abstract 434 words, main text 15,650, references 860, appendices 3,101, supplement 8,087, about 28,000 in total. Plain `wc -w` gives 32,239. The committed PDF has 53 pages, and it predates the last change to the Markdown, which added E16. This is not a short paper; the main text alone is longer than a typical full-length article. The abstract is more than twice the 200-word limit of the ASA style guide (a copy of the guide found by web search; I did not check the current *JASA* instructions), and it cites references, which the same guide discourages.

**What to cut or move.**

* *Section 7 (Bapat's estimator)* is unrelated to the paper's question. It reads as a commentary on another paper, and no result of Sections 3 to 6 uses it (line 461 says so). Remove it; if the questions matter, send them to the author or write a separate note.
* *Appendix D (repository map), the "Status of this draft" paragraph (line 13) and the [SYNTHETIC]/[REAL]/[PUBLISHED] tags, and labels such as "(E14, Q14f)" and "ALT1".* Replace them with a reproducibility statement and with cross-references to supplement sections. ALT1, ALT3 and ALT5 are cited for substantive claims (lines 52, 79, 138, 373, 469), but only ALT1 has settings in Appendix C (m29).
* *Sections 5.5 and 5.6 and Supplements S4 and S6* (running truncation, restarts, panel, multiple frequencies, adversarial search). Each is a paper-sized topic. Keep one paragraph each in the main text, or move them to a companion paper; at present they dilute the main message.
* *Proposition 10 and Table 2* can move to the supplement, with one paragraph in the main text.
* *Results proved only in repository notes.* Proposition A.6 and Lemma A.7 (line 139), the multi-frequency theorems (S4.3) and the panel theorem (S4.2) must be stated and proved in the supplement, or not used.

**Notation (m2, m3).** The overloading is heavy enough to slow a careful reader:

* $\kappa$ means truncated moments ($\kappa_2$, $\kappa_4$), the mean overshoot (Section 4.3), the bound on the largest factor in $p_{\eta,\kappa}$ (line 139) and $\sum_k w_k(1 - 1/k)$ (S4.3).
* $\rho$ means the signal ratio $\rho(k)$, stopping times (proofs of Corollary C and Theorem D), $\sum_k w_k/k$ (S4.3) and a grid ratio (Appendix B).
* $c$ means the truncation constant and the floor $c = 8.3$ (S6); $c_0$ is a third constant.
* $\mathcal{H}$ is two different classes (Proposition 2 and Section 4.1), alongside $H_0$, $(H^\uparrow)$ and the ladder height $H$.
* $G$ is the innovation law (Lemma 1) and the good event (Theorem D).
* $Z$ is a standard normal (Section 2), the scaled clipped return (Section 3.4) and log factors (Theorems 8 and 9).
* $B$ is $\tau^2/m_0$ (Section 2) and $\tau^2/m^\downarrow$ (Section 3.2).
* $k$ is the bet index, the GARCH break factor, the block length, the frequency and the segment index.
* $\sigma$ is the volatility and a stopping time (proof of Corollary F).
* $E$ is an e-variable and the events $E_s$.
* The numbering mixes letters (Theorems A, B, D; Corollaries C, E, F) with a numeric counter that Lemmas 1 and 2 and Proposition 1 do not share (there are both a Lemma 2 and a Proposition 2). "Theorem D.1" clashes with "Appendix D".

**Figures and tables.** The main text has eight dense tables and two figures.

* Figure 1 shows the triangles (mean of the ratio) above the dashed line labelled "cap" at persistence 0.95, which invites the reader to think the proposition is wrong. Relabel the line as the bound on the ratio of means, and say in the caption that the monitored quantity can exceed it.
* Figure 2 (left) plots a zero rate (Gaussian, $m = 98$) on a log axis, which draws a vertical artefact down to $10^{-4}$ near $m = 390$. Plot rate$/(m/A)$ on a linear axis instead: it shows the conservatism factor directly.
* Both figures carry internal labels in their titles ("Q9a, Q9b [synthetic]", "(E10, synthetic)").
* Figures are read from `../results/figures/`; a submission needs self-contained files.
* A figure of the monitor on one real day (the two statistics, the threshold, the alarms and the realised volatility) would help a reader more than Table 8, which would itself read better as a dot plot with intervals.
* Table 5 packs false alarm shares into brackets; split it into two tables or a figure.

**Would a reader in sequential analysis follow it?** Sections 3 and 4 yes, with effort. The obstacles are the notation, the forward references to supplement results, the repository labels, and prose dense with numbers: many sentences in Sections 5 and 6 carry five or more. The introduction's six-bullet contribution list mixes theory with an inventory of experiments.

**Concrete restructuring proposal** (main text of about 9,000 to 10,000 words, abstract under 200):

1. *Introduction* (about 900 words): the problem, the two obstacles, three contributions stated honestly (an anytime-valid robust monitor with a characterisation of what is and is not detectable for every procedure; an estimated-baseline guarantee with a clean null; evidence on simulated and real data), and related work including Section 2.3 above.
2. *Model, statistic and nulls*: current Sections 2 and 3.1.
3. *Validity*: Theorems A and B with the correct finite-horizon and expected-count statements (M2), Proposition 1, Corollary C in two lines, and the strengthened impossibility and speed-limit result for all procedures (M3).
4. *Estimated baseline*: Theorem D with the algorithm-free null as the main statement (M4); the plug-in negative result stated precisely, with its proof in the supplement.
5. *Delay*: Lemmas 4 and 5, Corollaries 6 and 7, Theorem 8, Theorem 9 in one paragraph (credited), Table 1; the boundary behaviour of decreases.
6. *Simulations*: validity (Tables 3 and 4, Figure 2), delay against strengthened benchmarks (M5), and the cost of the confidence sequence (Table 6).
7. *Application*: one coherent crypto study with the confidence-sequence baseline as primary (M1), the ratio rule and robust BOCPD as comparators, and the FOMC study as a pre-specified validation. Nifty 50 and the panel go to the supplement.
8. *Discussion.*

The supplement then holds Supplements S1 to S3, the extensions (current Sections 5.5 and 5.6, S4 and S6), the Nifty study and the experiment settings.

---

## 6. Major comments (ranked by importance)

**M1. The real-data evidence does not exercise the guarantee that the paper is about.** *Location:* Section 6 (lines 396 to 457), abstract (line 9), Section 8 (lines 481 to 483). *Problem:* Every main real-data analysis (Tables 7 and 8, Sections 6.2 to 6.4) uses the plug-in baseline. Section 3.4 shows that this baseline has no false alarm bound of the form $Cm/A$, and that it exceeds the known-baseline bound over five days (41% to 44% against 25%). Moreover, the null is false on most real days: crypto alarms come every 1.4 to 1.5 days, and since the monitor is blind for a day after each alarm, alarms come on average about half a day after calibration ends; e-BH discovers 86% to 87% of coin-days (line 453). So the anytime-valid property plays no role in the real-data evidence, and the precision comparisons show no advantage over a ratio rule. The paper says much of this, but it does not draw the conclusion: as written, the empirical case for an anytime-valid monitor is not made. *Fix:*

* (a) Make the confidence-sequence baseline (Theorem D) the primary configuration in Section 6, with the plug-in as a secondary analysis; E8 shows it is feasible on Bitstamp BTC.
* (b) Show the guarantee on realistic null streams built from real data. For example, take returns from calm periods and permute them at random within the stream: this gives an exchangeable null with the real marginal law, including zeros, ticks and jumps. Compare the false alarm rate with the bound for both baselines. Block-bootstrapped streams would then show what clustering does.
* (c) State what null is scientifically meaningful on real data. "Volatility equals yesterday's truncated variance up to 25%" is almost always false at one-minute resolution. Either motivate $b$ from an economic magnitude, or make the forecast-relative null of Section 4.5 the default for real data and report its cost.

**M2. Theorem B bounds the first $m$ steps after (re)calibration, not "each day".** *Location:* lines 9, 107 and 396. (Line 290 and the captions of Tables 3 and 4 correctly speak of the first day.) *Problem:* Under the null, $\mathbb{E}[M_t] \le t$, so the statistic accumulates. Ville's inequality bounds the probability of an alarm on day $d$ only by $dm/A$. In the paper's own null class (Section 3.2, item 1), I find 11.5% on day 5 and 28.0% on day 10 against the nominal 5%. Even for i.i.d. Gaussian returns, the conditional per-day hazard implied by Table 4 rises to about 4.6%. Corollary C gives the right long-run statement, $\mathbb{E}[N_m] \le m/A$ with restarts, which is an expected count. *Fix:*

* Restate every "per day" claim as "within the first day after (re)calibration" or as "expected number of false alarms per day", following items 1, 7 and 8 of Section 4.
* Add a remark relating the guarantee to the window-limited criterion of Lai (1998) and to Mei (2008).
* If a per-day probability guarantee is wanted, give a construction that has one: for example, restart the statistic at each period boundary, which gives $m/A$ per period at the cost of carrying no evidence across the boundary, or use the CUSUM-type statistic, which did not accumulate in my simulation, with a proved window bound.
* Report the conditional per-day hazard alongside the cumulative rates of Table 4.

**M3. The theoretical contribution is thin for a theory-and-methods venue; an easy strengthening is available.** *Location:* Sections 3 and 4, and line 475. *Problem:* As Section 2.1 shows, every theorem is either known (Theorems A and B, Proposition 2, Theorem 9), an elementary instance (Theorem 8, Corollary 7), a union-bound combination (Theorem D) or algebra (Proposition 10). The paper candidly says so. That suits an applications paper, but not *JASA* T&M or a theory paper in *Sequential Analysis*. *Fix (any two of these would change my view):*

* (i) *Impossibility and speed limits for every procedure* (Section 3.2, item 2). For any stopping rule with $\mathbb{P}(T \le n) \le n/A$ over the null, an alarm at $n$ on path $y$ requires $\prod_{s \le n}\max\lbrace 1, (B - y_s)/(B-1)\rbrace \ge A/n$ (decreases) or $\prod_s \max(1, y_s) \ge A/n$ (increases). This contains Corollary 3 ($B = \infty$) and gives a non-asymptotic floor of 45 steps for complete collapse at $c = 3$, $A = 7800$, against 57.7 for linear-bet e-detectors. It identifies the bet at $\lambda_{\max}$ as the extremal bettor, and it answers the open question on line 201.
* (ii) *The plug-in in the joint regime.* The relevant asymptotics for practice let the calibration length $n$ grow with $A$, as in the training-sample asymptotics of Chu, Stinchcombe and White (1996) and Horváth, Kokoszka and Zhang (2006). Characterise $\mathbb{P}(T_A \le m)$ for the plug-in, or for a plug-in inflated by $z$ standard errors, as a function of $(n, m, A)$. The present Proposition A.6 ($n$ fixed, $A \to \infty$) is the least informative regime.
* (iii) *Decreases at the boundary.* Since $\lambda^\ast$ is within 1% of $\lambda_{\max}$ in every decrease case (Section 3.2, item 3), give the delay theory for the bet at $\lambda_{\max}(1 - \epsilon)$ and the trade-off between growth and ruin when $Q$ puts mass at the cap.
* (iv) *A proof of the cap for the mean of the ratio*, or a counterexample (Figure 1 suggests that at persistence 0.95 the mean of the ratio exceeds $(1-\beta)/\alpha$), and the same for the clipped standardised statistic that the monitor actually uses.
* (v) *A rigorous second-order expansion for a single bet* (Section 4.3), whose conditions look routine to verify.

**M4. Theorem D's null depends on the algorithm and is much narrower than presented.** *Location:* lines 9, 30, 152 to 162, 170. *Problem:* $(H^\uparrow)$ and $(H^\downarrow)$ compare the monitored conditional mean with bet-weighted averages whose weights the algorithm picks from the data. The null hypothesis therefore depends on the test statistic, which is unusual and hard to interpret. The algorithm-free sufficient condition (line 162) requires, with $b = 0$, that no conditional mean during monitoring exceed the smallest one during calibration, which rules out any clustering. The abstract presents $\delta + m/A$ next to "every form of dependence", which applies only to a known baseline. *Fix:* State Theorem D with the algorithm-free condition as its hypothesis, written as a tolerance: "$\sup_{t > T_0}\mu_t \le (1+b)\min_{n_\tau < i \le T_0}\mu_i$". Move the bet-weighted version to a remark. Say plainly in the abstract and introduction that with an estimated baseline the guarantee holds when the conditional truncated variance varies by less than a factor $1+b$ across calibration and monitoring. Discuss Shekhar and Ramdas (2024), the closest prior construction.

**M5. The benchmarks favour the proposed method.** *Location:* Sections 5.1, 5.3 and 5.6; Tables 3, 5 and 8. *Problem:*

* The Gaussian likelihood-ratio, CUSUM and BOCPD competitors are calibrated on Gaussian data and then run on jumps and fat tails. Their failure there is expected and says little about what anytime-valid betting buys.
* BOCPD uses a Gaussian observation model, although robust versions exist (Knoblauch, Jewson and Damoulas 2018).
* The BJV-type test is "our reconstruction from the published idea and was not checked against the paper's equations, whose full text we could not reach" (line 385). A published comparison with a method whose paper the authors have not read is not acceptable. The paper is in the *Annals of Statistics* (2017), and the authors' own literature review lists an arXiv version (1502.00043).

*Fix:*

* Add a classical robust detector on the same statistic: a CUSUM or Shiryaev-Roberts rule on $X_t$, with a threshold calibrated by simulation, or by a bootstrap from the calibration window, to the same false alarm rate. This isolates the price and the value of a distribution-free threshold.
* Add a Student-$t$ or β-divergence BOCPD.
* Obtain BJV and implement it as published, or drop the comparison.
* Add an econometric monitoring procedure for variance (Horváth, Kokoszka and Zhang 2006) with its training-sample calibration.

**M6. Real-data inference: multiplicity, post hoc choices, dependence and the outcome measure.** *Location:* Section 6, Supplements S4.3 and S5; see Section 4.1 above. *Problem:* Several tests per table without adjustment; p-values from pooled alarms that are dependent across coins and days; secondary designs and reruns reported with p-values next to the primary ones; one comparison reported with three p-values. The outcome (next-day robust volatility on the flagged side of the baseline or of a slow forecast) partly measures volatility persistence, which any trailing rule captures, and it uses the same median-based scale as the monitor's truncation. *Fix:*

* Pre-specify one primary endpoint per data set, or say that none was.
* Use Holm within families, and block-bootstrap intervals by calendar day across coins.
* Mark everything else as exploratory.
* Evaluate "information beyond a forecast" with a forecast regression or a loss comparison (for example QLIKE), with standard errors robust to dependence, against a forecast that contains a short-horizon intraday component.
* Add an outcome based on realised variance, so that the result does not depend on the same robust scale.

**M7. Missing literature.** *Location:* lines 21, 25, 36 and Section 8. *Problem:* Section 2.3 lists verified references that a referee in the area would expect to see: econometric variance and GARCH monitoring (Chu, Stinchcombe and White 1996; Berkes et al. 2004; Horváth, Kokoszka and Zhang 2006; Dette and Gösmann 2020); robust change detection (Huber 1965; Unnikrishnan, Veeravalli and Meyn; Fearnhead and Rigaill 2019; Knoblauch et al. 2018); e-value change detection and robust confidence sequences (Shekhar and Ramdas 2024; Vovk et al. 2021; Dandapanthula and Ramdas 2025; Wang and Ramdas 2023; Mineiro 2022; Clerico 2025); e-values in finance (Wang, Wang and Ziegel; Choe and Ramdas); and false alarm criteria (Lai 1998 window criterion; Mei 2008; Tartakovsky, Nikiforov and Basseville 2014). Shekhar and Ramdas, Dandapanthula and Ramdas, Choe and Ramdas, and Wang, Wang and Ziegel are already in `docs/literature_review.md` but were dropped from the paper. *Fix:* Add a related-work paragraph that positions the paper against each group, and revise the novelty statement on line 36 accordingly. Update Horváth, Trapani and Wang to the published version if confirmed.

**M8. Length and structure.** *Location:* whole manuscript. *Problem:* About 28,000 words with appendices and supplement (main text about 15,700), a 434-word abstract, an off-topic Section 7, repository labels throughout, and results used in the main text but proved only in repository notes (Proposition A.6, Lemma A.7, the multi-frequency and panel theorems). *Fix:* Follow the restructuring proposal of Section 5: a main text of about 9,000 to 10,000 words and a 200-word abstract; remove Section 7 and Appendix D; put every result used in the main text, with its proof, in the supplement; consider splitting the extensions (running truncation and restarts, panel, multiple frequencies, adversarial search) into a companion paper.

---

## 7. Minor comments

* **m1.** *Lines 9 and 11.* The abstract cites three references and runs to 434 words; cut to about 200 words, with no citations.
* **m2.** *Notation.* Rename the overloaded symbols listed in Section 5; for example, use $\gamma_{\mathrm{ov}}$ for the overshoot, $c_{\min}$ for the S6 floor, $\mathcal{N}_\infty$ and $\mathcal{N}_B$ for the two null classes, and $\Omega_G$ for the good event.
* **m3.** *Numbering.* Use one counter for all results, or consistent letters; rename Theorem D.1 to Lemma A.2 to avoid the clash with Appendix D.
* **m4.** *Line 13.* Remove "Status of this draft" and the bracketed provenance tags from the submission version; keep them in the repository.
* **m5.** *Experiment labels* (E4, Q14f, ALT1). Replace them with references to sections or tables of the supplement.
* **m6.** *Proposition 1, line 109.* Alarming when $\bar M_t \ge A$ is valid by Lemma 2 and is never later than $\max(M^\uparrow, M^\downarrow) \ge 2A$; the direction is the side with the larger statistic. Use that rule, or explain why not. The gain is small, since the inactive side carries little wealth, but the current rule is dominated.
* **m7.** *Table 2.* Give $I$ to two significant figures: $I = 0.0003$ gives 23.8 days, not 21; the CSV value 0.000333 gives 21.4.
* **m8.** *Figure 1.* The "cap" line lies below the simulated mean of the ratio at persistence 0.95 and $k = 1000$ (2.20 against 2.0). Relabel it and say so in the caption; remove "(E10, synthetic)" from the title; rename the y-axis.
* **m9.** *Figure 2, left.* There is a log-scale artefact from a zero rate. Plot rate$/(m/A)$ on a linear scale; remove "Q9a, Q9b [synthetic]" from the titles.
* **m10.** *Files.* Figures are read from `../results/figures/`, and `paper/short_paper.pdf` (53 pages, committed before E16 was added) is out of date. Regenerate, and make the figure files part of the submission.
* **m11.** *Table 5.* The false alarm shares in brackets are measured over the 500 steps before the change, not over a day, so the reference value is $500/7800 = 6.4\%$ per side, not 5%; say so in the caption. See also item 14 of Section 4.
* **m12.** *Table 6.* When fewer than 90% of paths detect within 2500 steps, say how the median treats undetected paths: excluded, or censored at 2500. If excluded, the medians of the confidence-sequence column are biased downwards.
* **m13.** *Line 148.* Check that the predictable plug-in empirical Bernstein bound is Theorem 2 of Waudby-Smith and Ramdas (2024). I could not open the paper.
* **m14.** *Line 63, missing returns.* State the assumption that missingness given the past is independent of the size of the return, and include the missingness indicator in $\mathcal{F}_t$. Note that exchange outages tend to coincide with volatile periods, so the assumption is not innocuous. Say in the main text how the first return after a gap is handled; E7 sets it to missing, which is the right choice.
* **m15.** *Line 131 and Lemma 1.* "Falls exactly when $\sigma$ falls" needs strict monotonicity of $g_\tau$; give the condition. Also note that returns rounded to a fixed price tick are not of the form $\sigma_t\varepsilon_t$, so Lemma 1 does not give the tick laws a volatility meaning; this matters for the reading of Table S2.
* **m16.** *Theorem 9(a).* State that $\mathrm{WADD}_Q$ is evaluated with pre-change law $P^\ast$, and that $\tau$ and $m^\uparrow$ are treated as fixed (conditioning on $\mathcal{F}_0$).
* **m17.** *Section 4.1, post-change model.* "Independent of $\mathcal{F}_{\nu-1}$" should read "independent of $\mathcal{F}_{\nu-1}$ given $\tau$ and the baseline", since $Y$ is defined through them.
* **m18.** *Theorem 8 and Table S1.* The bound is for one side at level $A$, while the monitor runs each side at $2A$. Give the two-sided form, with $\log(2A/w_k)$ in the numerator, or state that E11 simulated one side at $A$. The repository confirms the latter: `results/e11a_delay.csv` matches the one-sided bound.
* **m19.** *Section 4.2 and Table 1.* Report $\lambda^\ast/\lambda_{\max}$ for the decrease cases (0.9917, 0.99987 and $1 - 2.4\times10^{-8}$ by my computation), and say that the interior condition of Theorem 9 holds only barely there.
* **m20.** *Line 201.* Give the default-grid floor for the two-sided monitor as well: 79 steps at $2A$ (my computation; the repository's E11 report has it), next to 73 at $A$.
* **m21.** *Line 139.* The share of calibration windows with positive drift (77% to 83%) is called a "numerical estimate". Give its Monte Carlo uncertainty and the number of windows.
* **m22.** *Corollary C and line 117.* The numerical check uses $A = 300$, far from the operating range; add a check at $A = 7800$ or $28{,}800$.
* **m23.** *Section 3.4, construction.* The split $n_\tau$ = 25% of the window is a design choice whose effect on delay is not reported; give one sensitivity line.
* **m24.** *Section 4.3 and S1.* For a single bet, check the conditions of the nonlinear renewal theorem (non-arithmetic bounded increments; $\log V_n$ slowly changing) and state the result as a proposition. Keep the Brownian approximation of $\zeta$ labelled as heuristic.
* **m25.** *Line 237.* "Within 3% of the simulated mean delay of a single bet in all 20 cases" covers cases not shown in Table S1 (jumps and Student $t_3$ at several $A$). Show them in the supplement, or restrict the claim to the cases shown.
* **m26.** *Section 2, intraday profile (line 63).* Errors in the estimated profile violate the null directly, and the expiry-day result (S5.1) suggests they matter on real data. At least one experiment should isolate the effect; the paper lists this as untested.
* **m27.** *Section 5.2, price ticks.* "$\tau$ floored at three ticks" (line 311) is a design choice not mentioned in Section 2. State the rule for $\hat s = 0$ (many zero returns; DOGE in Section 7), which the theory implicitly excludes ($\tau > 0$).
* **m28.** *Section 5.4, line 373.* "5.0% without noise" at one-second sampling: say which baseline was used, since with the true baseline the first-day rate should be about 1%.
* **m29.** *Appendix C.* ALT3 and ALT5 are cited for substantive claims (lines 79, 138 and 469) but have no settings in the appendix; add them, or cite the supplement section that describes them.
* **m30.** *References.* Add issue numbers and DOIs consistently (for example Mancini 2009, Lorden 1970 and 1971, Bahadur and Savage 1956, Lai 1998, Inclán and Tiao 1994). Update Horváth, Trapani and Wang to the *Journal of Time Series Analysis* version if confirmed. Check whether Larsson, Ramdas and Ruf, and Choe and Ramdas, have appeared.
* **m31.** *Line 365.* "Most of the width is intrinsic rather than a weakness of the bound": the confidence-sequence width (37%) is about twice the normal-approximation width (19%), and 8 points above the oracle constant-bet benchmark. "About half of the excess over the normal approximation is the price of finite-sample, time-uniform validity" would be more accurate.
* **m32.** *Line 255 and Supplement S2.* The split clip $c_\downarrow$ was evaluated only with the plug-in baseline. Since the paper argues that the split is legitimate with the confidence-sequence baseline, run at least the Gaussian row of Table S2 with it.
* **m33.** *Section 6.2, line 445.* State whether the event-day analysis was part of the Nifty pre-registration (`docs/report_parts/E7.md` suggests it was not), and adjust the language accordingly.
* **m34.** *Section 6.1.* The cross-venue agreement of BTC alarms (72%) checks the data rather than the method, since arbitrage links the venues; say so, or drop it.
* **m35.** *Section 6.5, line 457.* "Bapat (2026) does not name the source of its crypto data": this belongs with Section 7, which I recommend removing.

---

## 8. Errors versus matters of taste

**Errors or overstatements that must be fixed:**

* the per-day reading of Theorem B (M2; Section 4, items 1, 7 and 8);
* the "cap on the visible signal" (Section 4, items 3 and 6);
* the "about two days" bound contradicted by the paper's own delays (Section 4, item 12);
* "needs no window" (Section 4, item 16);
* "In all runs the settings were fixed" (Section 4, item 9);
* the BTC CPI p-values presented as rescuing a failed primary test (Section 4, item 15);
* the unadjusted post hoc subgroups (Section 4.1);
* the BJV comparison without access to the paper (M5);
* the main-text reliance on results proved only in repository notes (M8);
* the stale PDF (m10).

None of these is an error in a proof.

**Matters of substance where I disagree with the authors' framing:** the value of the guarantee on real data (M1), the narrowness of Theorem D's null (M4), the depth of the theory for the target venue (M3), and the choice of benchmarks (M5).

**Matters of taste:** the length and the extensions (M8), Section 7, the provenance tags, notation and numbering, and the figure design. I have offered concrete alternatives, but the authors may reasonably prefer others.
