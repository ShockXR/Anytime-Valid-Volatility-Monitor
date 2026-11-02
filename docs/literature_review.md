# Literature review

Every reference below was checked online in November 2026 (publisher pages, arXiv, RePEc/IDEAS, journal issue pages or author repositories). Where a detail could not be confirmed it is marked *unverified* instead of guessed. Summaries describe what each paper does and how it relates to this project. They are not a substitute for reading the paper.

## 1. The paper this project builds on

| Reference | Why it matters |
|---|---|
| Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance* 109, 105482. doi:10.1016/j.iref.2026.105482 | Source of the MAD-truncated, kernel-weighted local quarticity estimator. Title, author, journal and date confirmed online (ScienceDirect listing, pii S1059056026005952). Volume, article number and DOI are taken from the PDF. |

Prof. Bapat's earlier work is mainly on **sequential estimation** (stopping rules that decide how much data to collect). Confirmed example: Bapat, S. R. (2018), On purely sequential estimation of an inverse Gaussian mean, *Metrika* 81(8), 1005 to 1024, doi:10.1007/s00184-018-0665-0. Other listed works on purely sequential and two-stage procedures (in *Sequential Analysis*, *Metrika* and *Annals of the Institute of Statistical Mathematics*) were seen only on ResearchGate. Their details are unverified here. A sequential, anytime-valid extension of his estimator is close to his core expertise.

## 2. Realized variance and quarticity

| Reference | What it does | Relation to this project |
|---|---|---|
| Barndorff-Nielsen, O. E. and Shephard, N. (2002). Econometric analysis of realized volatility and its use in estimating stochastic volatility models. *JRSS B* 64(2), 253 to 280. doi:10.1111/1467-9868.00336 | Asymptotic theory for the error of realized variance. Its variance involves integrated quarticity. | Why quarticity is the "error bar" of variance (`docs/math.md`, section 1). |
| Barndorff-Nielsen, O. E. and Shephard, N. (2004). Power and bipower variation with stochastic volatility and jumps. *Journal of Financial Econometrics* 2(1), 1 to 37. doi:10.1093/jjfinec/nbh001 | Bipower variation (robust to jumps). RV minus bipower isolates jumps. Multipower quarticity. | Benchmarks `bipower_variation`, `tripower_quarticity`, `quadpower_quarticity`. |
| Andersen, T. G., Dobrev, D. and Schaumburg, E. (2012). Jump-robust volatility estimation using nearest neighbor truncation. *Journal of Econometrics* 169(1), 75 to 93. doi:10.1016/j.jeconom.2012.01.011 | MinRV and MedRV, plus their quarticity versions. | Benchmarks `min_rq` and `med_rq`. The constants π / (3π - 8) and 3π / (9π + 72 - 52√3) are checked by simulation in `tests/test_robust.py`. |
| Mancini, C. (2009). Non-parametric threshold estimation for models with stochastic diffusion coefficient and jumps. *Scandinavian Journal of Statistics* 36, 270 to 296 (issue and DOI unverified) | Threshold (truncated) realized variance that drops returns above a vanishing threshold `c Δ^ϖ`, `ϖ < 1/2`. | The threshold rate that makes truncation bias vanish (Issue 2 in `docs/estimator_issues.md`). Benchmark `threshold_rq`. |
| Corsi, F., Pirino, D. and Renò, R. (2010). Threshold bipower variation and the impact of jumps on volatility forecasting. *Journal of Econometrics* 159(2), 276 to 288. doi:10.1016/j.jeconom.2010.07.008 | Threshold bipower variation with smaller finite-sample bias. | Alternative robust estimator, not yet implemented. |
| Jacod, J. and Rosenbaum, M. (2013). Quarticity and other functionals of volatility: efficient estimation. *Annals of Statistics* 41, 1462 to 1484. doi:10.1214/13-AOS1115 (issue unverified) | Sums of local spot-volatility estimates, with a bias correction, give efficient estimators of IQ. | The rigorous version of "average local estimates" (eq. 16). Shows how a local estimator *should* be aggregated: plug in `σ̂_t⁴` and correct the bias of the power. |
| Mancino, M. E. and Sanfelici, S. (2012). Estimation of quarticity with high-frequency data. *Quantitative Finance* 12(4), 607 to 622. doi:10.1080/14697688.2012.664936 | Fourier estimator of spot and integrated quarticity. Studies its efficiency under noise. | Alternative quarticity estimator cited by Bapat. |
| Kristensen, D. (2010). Nonparametric filtering of the realized spot volatility: a kernel-based approach. *Econometric Theory* 26(1), 60 to 93. doi:10.1017/S0266466609090616 | Kernel-weighted realized variance as a spot volatility estimator, with consistency, a CLT and bandwidth choice. | The established theory for kernel spot estimators like the local `Q_t`. Gives the √h rate (Issue 4g). |
| Foster, D. P. and Nelson, D. B. (1996). Continuous record asymptotics for rolling sample variance estimators. *Econometrica* 64(1), 139 to 174 (DOI unverified) | Optimal window lengths and weights for rolling variance estimators. | Guidance for the bandwidth of the local estimator and of the adaptive bet. |

## 3. Jumps and microstructure noise

| Reference | What it does | Relation |
|---|---|---|
| Lee, S. S. and Mykland, P. A. (2008). Jumps in financial markets: a new nonparametric test and jump dynamics. *Review of Financial Studies* 21(6), 2535 to 2563. doi:10.1093/rfs/hhm056 | Tests each return for a jump by scaling it with local bipower volatility. | The same "return divided by local robust scale" idea that our monitor uses. |
| Zhang, L., Mykland, P. A. and Aït-Sahalia, Y. (2005). A tale of two time scales. *JASA* 100(472), 1394 to 1411. doi:10.1198/016214505000000169 | Two-scales realized variance, consistent under noise. | A noise-robust alternative (truncation alone does not handle noise, Issue 4b). |
| Hansen, P. R. and Lunde, A. (2006). Realized variance and market microstructure noise. *Journal of Business and Economic Statistics* 24, 127 to 161 (DOI unverified) | Noise is time-dependent and correlated with efficient returns. | Real noise is not the i.i.d. return noise of eq. (13). |
| Jacod, J., Li, Y., Mykland, P. A., Podolskij, M. and Vetter, M. (2009). Microstructure noise in the continuous case: the pre-averaging approach. *Stochastic Processes and their Applications* 119, 2249 to 2276. doi:10.1016/j.spa.2008.11.004 | Pre-averaging estimator under noise. | Candidate first stage before truncation (`alternatives/`). |
| Podolskij, M. and Vetter, M. (2009). Estimation of volatility functionals in the simultaneous presence of microstructure noise and jumps. *Bernoulli* 15(3), 634 to 658. doi:10.3150/08-BEJ167 | Modulated bipower variation, robust to both noise and jumps. | The main existing answer to the paper's stated gap. A natural benchmark to add. |
| Christensen, K., Oomen, R. C. A. and Podolskij, M. (2014). Fact or friction: jumps at ultra high frequency. *Journal of Financial Economics* 114(3), 576 to 599. doi:10.1016/j.jfineco.2014.07.007 (read from the abstract only) | At tick frequency, many apparent jumps in equities and currencies turn out to be short bursts of volatility. | The same pattern ALT3 finds in one-minute crypto returns; ALT3 adds what it means for truncated quarticity estimators. |
| Bajgrowicz, P., Scaillet, O. and Treccani, A. (2016). Jumps in high-frequency data: spurious detections, dynamics, and news. *Management Science* 62(8), 2198 to 2217. doi:10.1287/mnsc.2015.2234 (read from the abstract only) | With a correction for testing many returns at once, many jumps flagged by standard tests turn out to be spurious; the remaining ones are rare and often linked to news. | Supports caution in reading large one-minute returns as jumps (ALT3, Issue 5 of `estimator_issues.md`). |

## 4. Intraday seasonality

| Reference | What it does | Relation |
|---|---|---|
| Andersen, T. G. and Bollerslev, T. (1997). Intraday periodicity and volatility persistence in financial markets. *Journal of Empirical Finance* 4(2 and 3), 115 to 158. doi:10.1016/S0927-5398(97)00004-2 | Documents strong intraday volatility patterns and models them. | Why returns must be deseasonalised before monitoring. |
| Boudt, K., Croux, C. and Laurent, S. (2011). Robust estimation of intraweek periodicity in volatility and jump detection. *Journal of Empirical Finance* 18(2), 353 to 367. doi:10.1016/j.jempfin.2010.11.005 | Jumps bias periodicity estimates. Robust estimators fix this. | Our `data.intraday_scale` uses medians for the same reason. |

## 5. Robust statistics

| Reference | What it does | Relation |
|---|---|---|
| Huber, P. J. (1964). Robust estimation of a location parameter. *Annals of Mathematical Statistics* 35, 73 to 101 (issue and DOI unverified) | M-estimators. Clipping (Huber's ψ) gives bounded influence. | The robustness idea behind eq. (11) and Theorem 3. |
| Rousseeuw, P. J. and Croux, C. (1993). Alternatives to the median absolute deviation. *JASA* 88(424), 1273 to 1283. doi:10.1080/01621459.1993.10476408 | The Sn and Qn scale estimators, which are more efficient than MAD. | Possible replacement for MAD in the truncation rule. |

## 6. Change detection: classical

| Reference | What it does | Relation |
|---|---|---|
| Page, E. S. (1954). Continuous inspection schemes. *Biometrika* 41(1 and 2), 100 to 115. doi:10.1093/biomet/41.1-2.100 | The CUSUM procedure. | Benchmark `page_cusum`. |
| Lai, T. L. (1998). Information bounds and quick detection of parameter changes in stochastic systems. *IEEE Transactions on Information Theory* 44(7), 2917 to 2929 | Asymptotic lower bound on detection delay, about log A divided by an information rate, for rules whose average run length to a false alarm is at least A, including dependent data. | The yardstick for our delay theory (E11) and the information bound for forecast-relative monitoring (E10). |
| Lorden, G. (1971). Procedures for reacting to a change in distribution. *Annals of Mathematical Statistics* 42(6), 1897 to 1908 (DOI unverified) | Minimax optimality of CUSUM-type rules (Lorden's criterion). | The standard yardstick for detection delay. |
| Pollak, M. (1985). Optimal detection of a change in distribution. *Annals of Statistics* 13(1), 206 to 227. doi:10.1214/aos/1176346587 | Shiryaev-Roberts-Pollak procedure, asymptotically almost minimax. | Our main statistic is of Shiryaev-Roberts type. |
| Inclán, C. and Tiao, G. C. (1994). Use of cumulative sums of squares for retrospective detection of changes of variance. *JASA* 89(427), 913 to 923 (DOI unverified) | Retrospective CUSUM of squares (ICSS) for variance breaks. | Offline counterpart of what we do online. |
| Bibinger, M., Jirak, M. and Vetter, M. (2017). Nonparametric change-point analysis of volatility. *Annals of Statistics* 45(4), 1542 to 1578. doi:10.1214/16-AOS1499; arXiv:1502.00043 | Minimax-optimal **retrospective** test for volatility jumps in high-frequency semimartingales. | The strongest offline method for the same question. Ours is the online, anytime-valid counterpart. E13 benchmarks a version built from the published idea; the paper's full text could not be opened from our environment, so its exact standardisation and equation numbers are not yet checked against our code. |

## 7. Change detection and testing with e-values (anytime-valid inference)

| Reference | What it does | Relation |
|---|---|---|
| Ville, J. (1939). *Étude critique de la notion de collectif.* Thèse, Paris (numdam.org/item/THESE_1939__218__1_0) | Ville's inequality: a nonnegative supermartingale starting at 1 exceeds 1/α with probability at most α. | The inequality behind every false alarm bound in `docs/math.md`. |
| Bahadur, R. R. and Savage, L. J. (1956). The nonexistence of certain statistical procedures in nonparametric problems. *Annals of Mathematical Statistics* 27(4), 1115 to 1122 | No nontrivial test or confidence interval for a mean exists over all distributions with finite mean. | A one-sided sequential form explains why a fall in the variance cannot be detected without a bounded statistic (`docs/math.md` 5.3). |
| Fan, X., Grama, I. and Liu, Q. (2015). Exponential inequalities for martingales with applications. *Electronic Journal of Probability* 20, paper 1, 1 to 22. doi:10.1214/EJP.v20-3496 | Exponential inequalities for martingales with bounded increments. | The inequality of Lemma 6.1 behind the empirical Bernstein confidence sequence (`docs/math.md` 6). |
| Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2020). Time-uniform Chernoff bounds via nonnegative supermartingales. *Probability Surveys* 17, 257 to 317. doi:10.1214/18-PS321 | A unified treatment of time-uniform concentration, including a modern statement of Ville's inequality. | Background for the confidence-sequence baseline. |
| Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2021). Time-uniform, nonparametric, nonasymptotic confidence sequences. *Annals of Statistics* 49(2), 1055 to 1080. doi:10.1214/20-AOS1991 | Confidence sequences for means that hold at all times, including empirical Bernstein versions. | The family our confidence-sequence baseline belongs to (E8). |
| Shafer, G. (2021). Testing by betting: a strategy for statistical and scientific communication. *JRSS A* 184(2), 407 to 431. doi:10.1111/rssa.12647 | Evidence as the outcome of a bet against the null. | The intuition of the detector: wealth grows only if the null is false. |
| Ramdas, A., Grünwald, P., Vovk, V. and Shafer, G. (2023). Game-theoretic statistics and safe anytime-valid inference. *Statistical Science* 38(4), 576 to 601. doi:10.1214/23-STS894 | Review of e-values, test martingales and anytime-valid inference. | Background for `docs/math.md`, section 5. |
| Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *JRSS B* 86(1), 1 to 27. doi:10.1093/jrsssb/qkad009 | Betting confidence sequences for bounded means, including the aGRAPA bet. | Our bets on the bounded statistic `Y_t`, and the adaptive bet formula. |
| Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science* 2(2), 229 to 260. doi:10.51387/23-NEJSDS51 | E-detectors built from e-processes started at each time, with ARL control. | The framework our detector instantiates. |
| Shekhar, S. and Ramdas, A. (2024). Reducing sequential change detection to sequential estimation. *Proceedings of the 41st ICML*, PMLR 235, 44628 to 44642. arXiv:2309.09111 | Change detection with unknown pre-change law via overlapping confidence sequences. | The rigorous way to handle an estimated baseline (`docs/math.md`, 5.6). |
| Shekhar, S. and Ramdas, A. (2024). Nonparametric two-sample testing by betting. *IEEE Transactions on Information Theory* 70, from p. 1178. doi:10.1109/TIT.2023.3305867 (issue and last page unverified) | Sequential two-sample tests as betting games. | Alternative design: compare the recent window with the calibration window directly. |
| Choe, Y. J. and Ramdas, A. Combining evidence across filtrations. arXiv:2402.09698, listed as accepted in *JRSS B* (volume and DOI unverified) | Combines e-processes built on different information sets. Its example tests whether high-volatility days in stock returns occur at random. | The closest finance application of e-processes to volatility we found. |
| Wang, Q., Wang, R. and Ziegel, J. E-backtesting. *Management Science*, doi:10.1287/mnsc.2023.01659 (year, volume and pages unverified) | Anytime-valid, model-free backtests of Expected Shortfall with e-processes. | Shows that e-processes are already used in financial risk management. |

| Wang, R. and Ramdas, A. (2022). False discovery rate control with e-values. *JRSS B* 84(3), 822 to 852. doi:10.1111/rssb.12489 | The e-BH procedure: FDR control from e-values under any dependence. | The panel monitor (E15, `docs/theory_panel.md`) runs e-BH on per-asset e-processes. |
| Dandapanthula, S. and Ramdas, A. (2025). Multiple testing in multi-stream sequential change detection. arXiv:2501.04130 | Shows that a finite average run length is incompatible with nontrivial worst-case FDR in multi-stream change detection, and proposes an error-over-patience criterion. | E15 positions its block-wise guarantee against this result. Statements checked against the arXiv listing only. |
| Wang, H., Dandapanthula, S. and Ramdas, A. (2025). Anytime-valid FDR control with the stopped e-BH procedure. *Statistics and Probability Letters*, arXiv:2502.08539 (volume and pages unverified) | FDR control of e-BH at stopping times. | Basis of the "sticky" e-BH variant in E15. |
| Tavyrikov, Y., Goeman, J. J. and de Heide, R. (2026). Carefree multiple testing with e-processes. *Electronic Journal of Statistics* 20, 3178 to 3189. doi:10.1214/26-EJS2546 | Multiple testing with e-processes that may be stopped and resumed. | Related to the time-uniform statements of E15. |
| Chen, J., Zhang, W. and Poor, H. V. (2020). A false discovery rate oriented approach to parallel sequential change detection problems. *IEEE Transactions on Signal Processing* 68, from p. 1823 (DOI unverified) | FDR-oriented detection in parallel streams, parametric. | Earlier, non-anytime-valid approach to the panel problem. |
| Chen, Y. and Li, X. (2023). Compound sequential change-point detection in parallel data streams. *Statistica Sinica* 33, 453 to 474. doi:10.5705/ss.202020.0508 | Detects changes that affect a subset of many streams. | Same problem, likelihood-based. |

The full list behind E15, with how each reference was checked, is in `docs/theory_panel.md`, section 9.

## 8. Online monitoring of volatility

| Reference | What it does | Relation |
|---|---|---|
| Andersen, T. G., Todorov, V. and Zhou, B. (2025). Real-time detection of local no-arbitrage violations. *Quantitative Economics* 16(2), 459 to 495. doi:10.3982/QE2585 | Sequential stopping rules on high-frequency prices that flag local departures from the Itô semimartingale model in real time. | Closest real-time, high-frequency monitoring work found. Different target (no-arbitrage violations, not volatility regimes). |
| Horváth, L., Trapani, L. and Wang, S. Sequential monitoring for explosive volatility regimes. arXiv:2404.17885 (published venue unverified) | Weighted CUSUM monitoring for a switch to explosive GARCH volatility, on daily data. | Online volatility monitoring, but parametric and daily. |

## 9. Sequential fixed-accuracy estimation (alternative A8)

Background for `alternatives/sequential_estimation/`, which asks how many minutes of returns are needed to state the volatility level to a given relative accuracy. Its README has the full list (about 25 entries, each checked by web search).

| Reference | What it does | Relation |
|---|---|---|
| Stein, C. (1945). A two-sample test for a linear hypothesis whose power is independent of the variance. *Annals of Mathematical Statistics* 16(3), 243 to 258 | Two-stage procedure: a pilot sample fixes the size of the second stage, giving a fixed-width interval with exact coverage under normality. | One of the classical rules compared in A8. Needs a pilot of at least 30 returns there. |
| Robbins, H. (1959). Sequential estimation of the mean of a normal population. In U. Grenander (ed.), *Probability and Statistics (the Harald Cramér volume)*, Almqvist and Wiksell, Stockholm, 235 to 245 | Minimum risk point estimation: stop when the estimated loss from sampling more exceeds its cost. | On the truncated statistic, within 3% of the optimal risk from a pilot of 200 (A8). |
| Chow, Y. S. and Robbins, H. (1965). On the asymptotic theory of fixed-width sequential confidence intervals for the mean. *Annals of Mathematical Statistics* 36, 457 to 462. doi:10.1214/aoms/1177700156 | Purely sequential fixed-width interval for any law with finite variance, valid as the width shrinks. | The main classical rule of A8. On min(r², τ²) it covers 0.94 to 0.96 at about the oracle sample size; on raw r² it fails for heavy tails. |
| Nadas, A. (1967). *An extension of a theorem of Chow and Robbins on sequential confidence intervals for the mean*. PhD dissertation, Columbia University (a 1969 journal version is cited in the literature; volume and pages not verified) | Relative (proportional) accuracy version of Chow and Robbins. | The form used in A8, since volatility is quoted in relative terms. |
| Banerjee, S. and Mukhopadhyay, N. (2016). A general sequential fixed-accuracy confidence interval estimation methodology for a positive parameter. *Annals of the Institute of Statistical Mathematics* 68(3), 541 to 570. doi:10.1007/s10463-015-0504-2 | Fixed-accuracy intervals for positive parameters (an interval whose ratio of end points is fixed). | The modern statement of the problem A8 solves for a volatility level. |
| Konev, V. and Nazarenko, B. (2020). Sequential fixed accuracy estimation for nonstationary autoregressive processes. *Annals of the Institute of Statistical Mathematics* 72(1), 235 to 264. doi:10.1007/s10463-018-0689-2 | Fixed-accuracy estimation for dependent data. | Relevant to the clustering case, where A8 finds that rules built for independent data stop 7 to 47 times too early for the stationary level. |
| Politis, D. N. and Romano, J. P. (1994). The stationary bootstrap. *Journal of the American Statistical Association* 89(428), 1303 to 1313. doi:10.1080/01621459.1994.10476870 | Resampling of dependent series in random-length blocks. | Used in A8 to resample quiet BTC days with and without their dependence. |

Waudby-Smith and Ramdas (2024), Howard et al. (2021) and Bahadur and Savage (1956), all in section 7, are the other main references of A8: the first two give the confidence sequences it stops at a target width, the third the reason no rule can be exactly valid for the raw variance.

## 10. The gap

Our search (November 2026, not exhaustive) found:

* robust *offline* estimators of quarticity and volatility (sections 2 and 3),
* *retrospective* nonparametric tests for volatility breaks (Bibinger, Jirak and Vetter 2017),
* *online* monitoring that is parametric or uses daily data (Horváth, Trapani and Wang),
* *anytime-valid* e-detectors in general (Shin, Ramdas and Rinaldo 2024), and e-processes in risk backtesting (Wang, Wang and Ziegel).

We did not find a method that combines **jump-robust high-frequency volatility measures** with **anytime-valid sequential change detection**. That combination is the contribution this prototype explores. Since the search was not exhaustive, this claim should be stated with care when it is pitched.
