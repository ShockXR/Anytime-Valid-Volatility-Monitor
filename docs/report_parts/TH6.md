## TH6. False alarms after the first day: what Theorem B bounds, and two detectors with a bound for every day (referee M2 and M4)

Theory `docs/theory_per_period.md` (statements, proofs and the full tables); `docs/math.md` (remark after Theorem B; Theorem 6.4 restated, which the referee report calls Theorem D); code `src/avmonitor/periodic.py` (closed-form bounds, a vectorised simulation engine, `PeriodicVolatilityMonitor`); script `experiments/th6_per_period.py`; outputs `results/th6_*.csv` and `results/figures/th6_hazard.png`, `results/figures/th6_delay.png`; tests `tests/test_periodic.py`; iteration notes `docs/iterations/TH6.md`. The propositions and theorems are proved in the theory note; all numerical results are **[SYNTHETIC]**, with seeds fixed in the script.

**Question.** (M2) Theorem B proves `P(T ≤ m) ≤ m / A` for the first false alarm T, and the paper and README read it as "at most 5% per day". What does it bound after the first day? How large is the probability of a false alarm on day d in practice? What detector gives a bound for every day, and at what cost in delay? How do these guarantees relate to the window-limited criteria of Lai (1998) and to Mei (2008)? (M4) Restate Theorem D with a null hypothesis that does not depend on the bets the algorithm chooses, and say what it excludes.

**Set-up.** The default two-sided detector with the true baseline: X = min(r², τ²), Y = X / m0 ∈ [0, B]; six increase bets 0.02 to 0.5 and five decrease bets 0.1 to 0.95 times 1 / (B - 1); each side the average of its bets, alarming at 2A. A = 7800 for 390-step days and 28,800 for 1440-step days (m / A = 5%). Twenty days. Streams: i.i.d. Gaussian, Student t3, jump and tick returns (the four laws of E9, with τ and m0 computed exactly, `results/th6_laws.csv`); the referee's stream (`|r| = √m0`, so Y = 1, for nine days, then Gaussian), and the referee's family (Y = 1 for d - 1 days, then one Gaussian or two-point day); E16's bold play run for 20 days. Detectors: the plain Shiryaev-Roberts (SR) detector; SR restarted at 0 at every day boundary; the CUSUM-type e-detector `C_t = f_t max(C_{t-1}, 1)` (e-CUSUM). With a restart after every alarm (a recalibration with the true baseline) we record the probability of at least one alarm on day d and the expected number of alarms per day; for the first alarm, `p_d = P(T on day d)` and the hazard `h_d = P(T on day d | T > (d - 1) m)`. 40,000 paths per stream (bold play and park and strike 30,000; delays 20,000 per cell).

**Summary.**
* **Theorem B is a bound for the first day after a (re)start, and no bound below 1 holds for all later days (Proposition 1, proved).** For day d it gives only `d m / A`, which no example reached (at most 35.7% against 50% at d = 10). The i.i.d. law Y = 1 is in the null and makes the default detector alarm with certainty on day 40. A remark after Theorem B in `docs/math.md` now says this, and the bullet "5% per trading day" there is corrected.
* **Measured.** For i.i.d. returns the hazard of the plain detector rises from 0.98% to 1.52% on day 1 to 4.60% to 4.70% on average over days 11 to 20 (390-step days); no single i.i.d. day is significantly above 5%. Bold play, which is inside the null class, averages 5.25% over days 11 to 20 (5.23% pooled over path-days, standard error about 0.06 points; the pooling over days 11 to 20 was chosen after seeing the results). For the referee's family the probability of a false alarm on day d is 11.3% at d = 5, 28.0% at d = 10 and 65.1% at d = 20. Corollary C held everywhere: 0.037 to 0.048 false alarms per day over 20 days, against 0.05.
* **Two detectors with a bound for every day (Theorems 1 and 2, proved).** (1) SR restarted at every day boundary, with a known baseline: at most m / A = 5% on every day, given everything before the day, under the martingale-difference null of Theorem B; `PeriodicVolatilityMonitor(period=m, baseline="known", tau=..., m0=...)` implements it exactly. With the confidence-sequence baseline the restart gives δ + m / A for each day of a segment, marginally and not given the past, while the tolerance of Theorem D holds to the end of that day (Corollary 1 of the theory note); with the plug-in baseline nothing per day is proved. Measured on every stream and day: the probability of a false alarm on a day was at most 2.66% (390) and 3.52% (1440), and the hazard at most 2.75% (390) and 3.52% (1440), all for bold play. (2) The e-CUSUM: by time reversal, `P(C_s ≥ x) ≤ 1 / x` at every s under an i.i.d. (or exchangeable) null with a known baseline and fixed bets, which gives at most `(m + log 2A + (log 6 + log 5) / 2) / A` for the probability of an alarm in any window of m steps: 5.146% (390, 7800) and 5.044% (1440, 28,800). This bound is unconditional; given the state at the start of the window it can fail badly (about 32% for the next step from a state just below the threshold). Measured: at most 0.18% (i.i.d.) and 0.40% (bold play).
* **The e-CUSUM bound fails without exchangeability (Proposition 3, proved).** A martingale-difference null stream that parks the statistic just below the threshold and then strikes has a window probability of at least about 1 / B = 11% in the limit, for every m and A. At full scale it passes the i.i.d. bound on day 17 (390) and day 16 (1440) and reaches 6.16% (390) and 6.44% (1440) on day 20.
* **Delay.** At a proved 5% for every day (known baseline), the restarted SR detector's mean delay was 1.06 to 1.44 times the plain detector's (390-step days) and 1.02 to 1.14 times (1440), within the renewal bound of Proposition 2 in every cell; the e-CUSUM at its certified A was 1.30 to 1.83 and 1.24 to 1.66 times slower. The plain detector is 14% to 31% faster than a fresh start because it enters the change with the statistic built up under the null, the same accumulation that raises its later-day false alarm probability.
* **Theorem D (M4).** Restated in `docs/math.md` with (H↑) `sup_{T<t≤T+h} μ_t ≤ (1 + b) min_{n_τ<i≤T} μ_i` and (H↓) `inf_{T<t≤T+h} μ_t ≥ (1 - b) max_{n_τ<i≤T} μ_i` (h = m for the bound within m steps), plus a part 4 with no assumption on the conditional means: `P(alarm within m) ≤ δ + P(tolerance fails within m) + m / A`. The proof goes through unchanged except that a weighted average lies between the smallest and largest of its terms; the bet-weighted version is kept as a remark. Under GARCH(1,1) clustering the band needed per path has a median of 0.90 to 5.01; no path of any cell satisfied the tolerance with b = 0.25, and parts 1 to 3, which need it almost surely, never apply.

### Q-TH6a. The plain detector after the first day

Hazard `h_d` of the plain SR detector in percent (averages of `h_d` over the days shown); standard errors of a single `h_d` at most 0.17 points for the i.i.d. laws, 0.21 for bold play and 0.24 for the calm stream. Last column: false alarms per day over days 1 to 20 with a restart after every alarm, which Corollary C bounds by 0.05 (`results/th6_hazard_390.csv`, `results/th6_hazard_1440.csv`, `results/th6_adversary.csv`; figure `th6_hazard.png`, left).

| day length, A | stream | day 1 | days 2 to 5 | days 6 to 10 | days 11 to 20 | largest `h_d` (day) | false alarms per day |
|---|---|---|---|---|---|---|---|
| 390, 7800 | Gaussian | 0.98 | 3.02 | 4.13 | 4.70 | 5.18 (19) | 0.0376 |
| 390, 7800 | Student t3 | 1.49 | 3.39 | 4.11 | 4.68 | 5.06 (19) | 0.0392 |
| 390, 7800 | jumps | 0.99 | 3.03 | 4.05 | 4.60 | 5.03 (19) | 0.0371 |
| 390, 7800 | ticks | 1.52 | 3.29 | 4.16 | 4.66 | 4.88 (18) | 0.0388 |
| 390, 7800 | bold play (E16) | 2.53 | 4.18 | 4.91 | 5.25 | 5.44 (14) | 0.0463 |
| 390, 7800 | Y = 1 for 9 days, then Gaussian | 0 | 0 | 28.16 on day 10 | 5.47 | 28.16 (10) | 0.0370 |
| 1440, 28,800 | Gaussian | 2.04 | 3.64 | 4.50 | 4.68 | 4.94 (18) | 0.0416 |
| 1440, 28,800 | Student t3 | 2.44 | 3.92 | 4.46 | 4.57 | 4.77 (13) | 0.0419 |
| 1440, 28,800 | jumps | 2.07 | 3.67 | 4.34 | 4.67 | 4.93 (17) | 0.0411 |
| 1440, 28,800 | ticks | 2.37 | 3.95 | 4.58 | 4.78 | 4.95 (13) | 0.0430 |
| 1440, 28,800 | bold play (E16) | 3.13 | 4.56 | 5.00 | 5.02 | 5.27 (14) | 0.0479 |
| 1440, 28,800 | Y = 1 for 9 days, then Gaussian | 0 | 0 | 34.03 on day 10 | 5.28 | 34.03 (10) | 0.0407 |

The cumulative probabilities agree with E9: `P(T ≤ 5 m)` = 12.4% (Gaussian), 14.2%, 12.5% and 13.8% here, against 12.8%, 14.3%, 13.2% and 13.7% in E9 (4000 paths). The first-alarm probability over 20 days is 55.5% to 56.9% for the i.i.d. laws (Ville's bound: 1).

### Q-TH6b. The referee's family

Y = 1 for d - 1 days (every SR bet then sits at `(d - 1) m` exactly), then one Gaussian day; probability of a false alarm on day d in percent, 40,000 paths per d, standard errors at most 0.25 points (`results/th6_calm_family.csv`; figure `th6_hazard.png`, middle).

| day d | 1 | 2 | 5 | 10 | 15 | 20 |
|---|---|---|---|---|---|---|
| Ville's bound `min(1, d m / A)` | 5 | 10 | 25 | 50 | 75 | 100 |
| plain SR, 390-step days | 0.90 | 3.04 | 11.33 | 28.02 | 46.82 | 65.07 |
| plain SR, 1440-step days | 1.82 | 5.06 | 15.13 | 34.30 | 52.84 | 70.82 |
| SR restarted daily, 390 (bound 5) | 0.90 | 0.95 | 0.91 | 0.85 | 0.97 | 0.93 |
| SR restarted daily, 1440 (bound 5) | 1.82 | 2.01 | 1.93 | 2.08 | 1.94 | 1.97 |
| e-CUSUM, 390 (bound 5.146) | 0.02 | 0.05 | 0.04 | 0.03 | 0.08 | 0.07 |
| e-CUSUM, 1440 (bound 5.044) | 0.06 | 0.06 | 0.04 | 0.06 | 0.05 | 0.06 |

With a two-point day `{0, B}` instead of a Gaussian one the plain detector reaches 70.09% (390) and 73.16% (1440) on day 20, the restarted detector stays at 1.73% to 2.00% (390) and 2.59% to 2.84% (1440), and the e-CUSUM at 0.14% to 0.27%.

### Q-TH6c. Per-day constructions on the same streams

Probability of at least one false alarm on day d (restart after every alarm), smallest and largest over the 20 days, in percent. "No restart" for the e-CUSUM is the probability that its statistic is at or above the threshold at some step of day d without any restart, the event Theorem 2(c) bounds (same files as Q-TH6a; figure `th6_hazard.png`, right).

| day length | stream | SR restarted daily (bound 5) | e-CUSUM (bound 5.146 or 5.044) | e-CUSUM, no restart |
|---|---|---|---|---|
| 390 | Gaussian | 0.87 to 1.03 | 0.04 to 0.10 | at most 0.12 |
| 390 | Student t3 | 1.47 to 1.65 | 0.08 to 0.17 | at most 0.18 |
| 390 | jumps | 0.93 to 1.16 | 0.05 to 0.10 | at most 0.12 |
| 390 | ticks | 1.38 to 1.63 | 0.06 to 0.14 | at most 0.15 |
| 390 | bold play | 2.23 to 2.66 | 0.28 to 0.38 | not run |
| 390 | Y = 1 for 9 days, then Gaussian | 0 to 1.03 | 0 to 0.10 | at most 0.11 |
| 1440 | Gaussian | 1.78 to 2.08 | 0.06 to 0.11 | at most 0.11 |
| 1440 | Student t3 | 2.21 to 2.57 | 0.10 to 0.17 | at most 0.18 |
| 1440 | jumps | 1.90 to 2.10 | 0.05 to 0.10 | at most 0.10 |
| 1440 | ticks | 2.29 to 2.57 | 0.09 to 0.15 | at most 0.15 |
| 1440 | bold play | 3.11 to 3.52 | 0.26 to 0.40 | not run |
| 1440 | Y = 1 for 9 days, then Gaussian | 0 to 2.06 | 0 to 0.11 | at most 0.11 |

The restarted detector repeats its first day on every day. For i.i.d. returns the e-CUSUM is about 30 to 130 times below its bound (13 times for bold play): the bound is a union over the m steps of the window, which only laws with a large jump approach.

### Q-TH6d. Martingale-difference nulls: park and strike

Proposition 3 against the default e-CUSUM: hunt with `{0, B}` draws until one draw could reach the threshold, park with Y = 1 (which leaves the statistic unchanged), then play bold on day d. Every draw has conditional mean exactly 1. 30,000 paths; standard errors at most 0.14 points (`results/th6_park_strike.csv`).

| day d | 1 | 5 | 10 | 15 | 17 | 20 |
|---|---|---|---|---|---|---|
| share parked before day d, 390 (%) | 0 | 2.7 | 5.9 | 9.1 | 10.3 | 12.1 |
| P(false alarm on day d), 390 (i.i.d. bound 5.146) | 0.26 | 1.60 | 3.19 | 4.70 | 5.43 | 6.16 |
| Lemma 2 (conditional bound) averaged over paths, 390 | 5.00 | 6.38 | 8.01 | 9.53 | 10.12 | 11.00 |
| P(false alarm on day d), 1440 (i.i.d. bound 5.044) | 0.28 | 1.63 | 3.42 | 4.96 | 5.52 | 6.44 |

On day 20 the probability is 7.3 (390) and 9.8 (1440) standard errors above the i.i.d. bound. Scaled down to periods of 20 steps with A = 400 (`results/th6_park_scaled.csv`, 20,000 paths; i.i.d. bound 7.1%), the parked share reaches 57% and the probability 27.0% by period 120.

### Q-TH6e. Delay cost

A change of the return standard deviation by the factor shown at a uniformly random step of day 3 (Gaussian returns), horizon 12 days; delay counted from the change to the first alarm at or after it, including the alarm step; every detector found every change. Ratios of mean delays to the plain SR detector (standard errors of the means at most 0.7 steps). The restarted detector uses A; the variant that keeps a carry of m at each boundary uses 2A, so its bound `(m + m) / 2A` is also 5%; the e-CUSUM uses A = 8027.7 (390) and 29,053.4 (1440), where its window bound is 5% (`results/th6_delay.csv`, `results/th6_delay_by_position.csv`, `results/th6_delay_renewal.csv`; figure `th6_delay.png`).

| day length | sd factor | plain SR mean delay | restarted daily | restarted, carry m | e-CUSUM certified | restarted daily: measured / Proposition 2 bound | restarted daily by position: first half / 50% to 90% / last 10% |
|---|---|---|---|---|---|---|---|
| 390 | 0.5 | 70.8 | 1.26 | 1.23 | 1.72 | 89.2 / 115.0 | 82 / 91 / 122 |
| 390 | 0.67 | 97.5 | 1.32 | 1.26 | 1.75 | 128.4 / 160.4 | 112 / 142 / 156 |
| 390 | 0.8 | 161.9 | 1.44 | 1.34 | 1.83 | 232.6 / 274.9 | 196 / 275 / 239 |
| 390 | 1.25 | 113.4 | 1.21 | 1.22 | 1.61 | 137.1 / 151.8 | 124 / 149 / 155 |
| 390 | 1.5 | 43.2 | 1.09 | 1.14 | 1.39 | 47.3 / 54.7 | 45 / 46 / 63 |
| 390 | 2 | 20.1 | 1.06 | 1.11 | 1.30 | 21.2 / 25.4 | 21 / 20 / 27 |
| 1440 | 0.5 | 87.4 | 1.09 | 1.13 | 1.58 | 94.8 / 123.1 | 92 / 88 / 135 |
| 1440 | 0.67 | 120.7 | 1.10 | 1.14 | 1.60 | 132.3 / 168.4 | 127 / 123 / 198 |
| 1440 | 0.8 | 201.7 | 1.14 | 1.16 | 1.66 | 229.3 / 279.9 | 214 / 226 / 323 |
| 1440 | 1.25 | 139.7 | 1.07 | 1.13 | 1.50 | 149.8 / 170.0 | 145 / 146 / 193 |
| 1440 | 1.5 | 51.8 | 1.03 | 1.10 | 1.32 | 53.6 / 61.1 | 53 / 52 / 64 |
| 1440 | 2 | 23.7 | 1.02 | 1.08 | 1.24 | 24.1 / 28.1 | 24 / 24 / 26 |

Proposition 2 bounds the restarted detector's expected delay by `E[min(D₀, m - k)] + P(D₀ > m - k) E[min(D₀, m)] / P(D₀ ≤ m)` for a change k steps into a day, where `D₀` is the delay from a fresh start at the change. The cost sits where the rest of the day is shorter than a typical detection. A carry of m at 2A helps a little for decreases with 390-step days and costs a little elsewhere. The plain detector also had 4.6% to 4.9% (390) and 6.6% to 7.0% (1440) false alarms in the two and a half days before the change, against 2.1% to 2.4% and 4.4% to 4.7% for the restarted detector and at most 0.2% for the e-CUSUM. The renewal part was added after the delay results had been seen, to explain them (post hoc).

### Q-TH6f. Theorem D: what the tolerance excludes (M4)

The smallest band b for which (H↑) and (H↓) hold over a calibration window and one monitored day of the same length, with τ from the first quarter of the window and `μ_t` exact (for the stochastic volatility, given the volatility path), 2000 paths per cell (`results/th6_tolerance.csv`).

| model | calibration and day | median b, increase side | median b, decrease side | median b, both | 90% quantile, both | share of paths with b ≤ 0.5 | share with b ≤ 1 |
|---|---|---|---|---|---|---|---|
| GARCH(1,1), persistence 0.95 | 390 | 1.46 | 0.58 | 1.46 | 2.08 | 0 | 8.2% |
| GARCH(1,1), persistence 0.99 | 390 | 2.58 | 0.71 | 2.58 | 4.82 | 0 | 4.2% |
| GARCH(1,1), persistence 0.999 | 390 | 1.62 | 0.62 | 1.62 | 3.54 | 0 | 28.4% |
| GARCH(1,1), persistence 0.9999 | 390 | 0.90 | 0.51 | 0.90 | 2.15 | 1.8% | 55.6% |
| log-normal SV, log-sd 0.15, half-life 390 | 390 | 0.54 | 0.37 | 0.58 | 1.07 | 32.7% | 87.3% |
| GARCH(1,1), persistence 0.95 | 1440 | 2.05 | 0.67 | 2.05 | 2.72 | 0 | 0 |
| GARCH(1,1), persistence 0.99 | 1440 | 5.01 | 0.82 | 5.01 | 7.82 | 0 | 0 |
| GARCH(1,1), persistence 0.999 | 1440 | 4.17 | 0.81 | 4.17 | 9.45 | 0 | 3.6% |
| GARCH(1,1), persistence 0.9999 | 1440 | 2.38 | 0.73 | 2.38 | 5.89 | 0 | 16.6% |
| log-normal SV, log-sd 0.15, half-life 390 | 1440 | 1.22 | 0.55 | 1.22 | 2.14 | 0.1% | 35.0% |

No path of any cell satisfied the tolerance with b = 0.25. These are per-path bands. Parts 1 to 3 of Theorem D need the tolerance almost surely, which no band below 1 achieves under GARCH (in every GARCH cell 44% to 100% of the paths need b > 1), so they never apply to these models; the per-path bands are what part 4 charges through the probability that the tolerance fails, which is 1 at b = 0.25 in every cell. Since b must be below 1 on the decrease side, the two-sided theorem excludes GARCH clustering at these persistences; the theorem's main case is a constant conditional law (Corollary 6.5). With a known baseline (Theorem B) only the conditional mean against m0 matters, and the calibration window plays no role.

### Relation to the literature

Lai (1998) bounds `sup_k P_∞(k ≤ T < k + m)`; Theorem 2 gives this criterion for the e-CUSUM uniformly over exchangeable nulls, Theorem 1 gives it for the restarted detector with a known baseline over the martingale-difference null for windows aligned with days (and `(2m + κ) / A` otherwise), and the plain detector satisfies it only for k = 1. Pergamenchtchikov and Tartakovsky (2018) use a maximal local probability of a false alarm in a window of fixed length, with a known pre-change law; Theorem 1 bounds the corresponding hazard by m / A for every day. Mei (2008) shows that with composite (mixture) pre-change models the average run length can mislead; Theorem A holds here for every null law while the day-40 example and the referee's family show day-level probabilities far above 5%, and Corollary C is a rate statement of the kind Mei recommends. Shekhar and Ramdas (2024) start a confidence sequence at every time and declare a change when their intersection is empty; Theorem D uses one confidence sequence on a calibration window as the baseline of an e-detector, needs a change-free calibration window, and gives a finite-horizon bound and a long-run count. All four were checked by web search only (arxiv.org and the publishers' pages could not be opened from this machine); the theory note says which statements rest on abstracts or secondary sources.

### Reading

"5% per day" is true of the first day after each (re)calibration, and of the long-run average number of false alarms per day, but not of every day. For ordinary i.i.d. returns the later-day hazard of the plain detector levels off just below 5%, so the misreading is mild there; for a null stream that is calm and then ordinary, or that watches the detector, it is wrong by a factor up to 13 (65% on day 20) or by a few tenths of a point. If a per-day statement is wanted, restarting the statistics at every day boundary gives it, with a known baseline, for every null law of Theorem B at a delay cost of 2% to 14% for 1440-step days and 6% to 44% for 390-step days; the e-CUSUM gives a window bound that also holds for windows that are not aligned with days, but only for exchangeable nulls, and it is slower. `PeriodicVolatilityMonitor(period=m, baseline="known", tau=..., m0=...)` implements this in the streaming monitor without changing `VolatilityMonitor`. With an estimated baseline only the confidence-sequence version carries a per-day proof (δ + m / A for each day, not given the past, while the tolerance of Theorem D holds), and the default plug-in version carries none.

### Caveats

* The simulations and Theorems 1 and 2 are for a known baseline. With the confidence-sequence baseline, Corollary 1 of the theory note gives `δ + m / A` for each day of a segment, marginally (given the observed past the probability is bounded only by `P(bounds failed | past) + m / A`) and only while the tolerance of Theorem D holds to the end of that day; the plug-in baseline is not covered by any per-day statement (Proposition A.6 of `docs/theory_running_tau.md`).
* Theorem 2 needs fixed bets; the adaptive bet of `detectors.betting_edetector` breaks the time reversal. Its bounds are unconditional; given the past only Lemma 2 holds. Whether `log A + H(w)` in Theorem 2(b) can be removed for the average over bets is open.
* The pooled bold-play hazard over days 11 to 20 and the renewal comparison were chosen after seeing the results. The tolerance numbers for stochastic volatility use the filtration that contains the volatility path; with returns only the conditional mean is smoother and the band needed is smaller.
* Park and strike is a lower-bound construction; how close a martingale-difference stream can come to the restarted detector's 5% after day 1 is not known (bold play reaches 2.75%).

### Appendix. Sentences that read Theorem B as a per-day or per-session probability

Line numbers are those of commit 126588d; each entry quotes the text so that it can be found after later edits. Only `docs/math.md` was edited here (both entries below are done). Measured one-day rates described as "per day" for simulations in which every path is a fresh first day (for example "BOCPD raises a false alarm on 73% of days") are correct and are not listed.

| file, line | current text | proposed replacement |
|---|---|---|
| `docs/math.md`, 126 (done) | "gives at most a 5% chance of any false alarm per trading day" | "gives at most a 5% chance of a false alarm within the first trading day after each (re)calibration, and at most 0.05 expected false alarms per trading day in the long run (Corollary C). It does not bound the chance of a false alarm on each later day (remark above)." |
| `docs/math.md`, after Theorem B (done) | (no remark) | "Remark (what Theorem B bounds)", with the day-40 example and the TH6 numbers |
| `README.md`, 24 | "reaches 2.62% per day (95% interval 2.50% to 2.75%) against 5%" | "reaches 2.62% within the first day after calibration (95% interval 2.50% to 2.75%) against 5%" |
| `README.md`, 29 | "(at most 3.5% false alarms per day in simulation, against a 5% bound)" | "(at most 0.035 false alarms per day over 30-day simulated streams, against a long-run bound of 0.05 per day)" (E14 Q14h: a long-run rate, not a first-day probability) |
| `README.md`, 49 | "A = (steps per day) / (false alarm budget per day)" | "A = m / α with m the steps per day: at most α chance of a false alarm within the first day after each (re)calibration and at most α expected false alarms per day in the long run, not α on every later day (TH6)" |
| `README.md`, 139 | "with budget 0.05 per day, a monitor on a stream with no change raises a false alarm on at most about 5% of days in the long run" | "with budget 0.05 per day, a monitor on a stream with no change has at most a 5% chance of a false alarm in the first day after each (re)calibration, and at most 0.05 expected false alarms per day in the long run (so at most 5% of days with one, on average); on a particular later day the chance can be higher (TH6: about 5% for i.i.d. returns, 28% on a day after nine calm days), unless the statistics restart every day with a known baseline (`PeriodicVolatilityMonitor(baseline="known")`)" |
| `PITCH_BAPAT.md`, 27 | "P(false alarm within m steps) ≤ m / A ... Measured: 1.1% to 1.4% per day against a 5% bound" | "P(false alarm within the first m steps after calibration) ≤ m / A ... Measured: 1.1% to 1.4% within the first day against a 5% bound; on later days the bound is on the long-run average, 0.05 false alarms per day, unless the monitor restarts every day (TH6)" |
| `PITCH_BAPAT.md`, 43 | "reached 2.62% per day against 5%" | "reached 2.62% within the first day against 5%" |
| `docs/limitations.md`, 9 | "E10a found about 5% per day with a one-day calibration window" | "E10a found about 5% within the first day after a one-day calibration window" |
| `docs/limitations.md`, 10 | (add after the first sentence) | "Theorem B itself bounds only the first day after each (re)calibration: on later days of a segment the chance of a false alarm can exceed 5% (TH6: about 5% after a week for i.i.d. returns, 5.25% for a null stream that watches the detector, 28% on a day after nine calm days); only the long-run average of 0.05 per day (Corollary C) holds. With a known baseline, restarting the statistics every day gives 5% for every day; with the confidence-sequence baseline it gives δ + m / A per day, not given the past; with the plug-in, nothing (`docs/theory_per_period.md`)." |
| `docs/limitations.md`, 11 | "gives up to about 25% chance of at least one false alarm per day across the group (union bound)" | "gives, by the union bound, up to about 25% chance of at least one false alarm across the group in the first day after calibration, and up to 0.25 expected false alarms per day across the group in the long run" |
| `docs/limitations.md`, 13 | "(E16: 2.62% per day; 1.4 times over five days)" | "(E16: 2.62% within the first day; 1.4 times over five days)" |
| `paper/short_paper.md`, 9 | "the probability of a false alarm within $m$ steps is at most $m/A$" | "the probability of a false alarm within the first $m$ steps after (re)calibration is at most $m/A$" |
| `paper/short_paper.md`, 91 | "the probability of a false alarm on day 5 is 11.5%, and on day 10 it is 28.0%. Even for i.i.d. Gaussian returns, Table 3 implies ... of 1.1% on day 1, about 3.1% on days 2 to 5, 4.2% on days 6 to 10 and 4.6% on days 11 to 20." | "is 11.3% on day 5, 28.0% on day 10 and 65.1% on day 20 (40,000 paths). Even for i.i.d. Gaussian returns the probability of a first false alarm on day d, given none before, is 0.98% on day 1 and 3.0%, 4.1% and 4.7% on average over days 2 to 5, 6 to 10 and 11 to 20, and for a null stream that watches the detector (bold play) 5.25% over days 11 to 20. With a known baseline, restarting the statistics at every day boundary bounds every day by $m/A$, given the past, and the CUSUM statistic bounds any window of $m$ steps by $(m + \log 2A + (\log 6 + \log 5)/2)/A$ (5.146%) when the returns are i.i.d., but not under martingale differences (TH6)." |

The paper's two TH6 slots: the M2 slot after its restatement of Theorem B can take Theorem 1 (restart at every day boundary with a known baseline: at most m / A for every day, given the past), Corollary 1 (confidence-sequence baseline: δ + m / A per day, marginal), Proposition 2 and the delay ratios of Q-TH6e, with Theorem 2 and Proposition 3 for the e-CUSUM; the slot after the paper's restated Theorem D can take `docs/math.md` Theorem 6.4 as it now stands: the tolerance with horizon h = m for the bound within m steps (the paper states it with h = ∞, which is stronger than needed), part 4 (no assumption: `δ + P(V_m) + m / A`), and the GARCH numbers of Q-TH6f.
