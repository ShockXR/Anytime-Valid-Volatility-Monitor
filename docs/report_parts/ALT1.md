# Report part ALT1: two alternative prototypes

Scripts, tests, READMEs and full results are in `alternatives/quarticity_stream/` (A3) and `alternatives/preaveraging/` (A1). Every number below is copied from the `results/` folders of those prototypes. [SYNTHETIC] marks simulated data, [REAL] marks Binance spot one-minute bars (1 July to 31 December 2025).

---

## A3. Monitoring quarticity itself: the e-detector on `min(r^4, tau^4)`

**Question.** The main monitor bets on the truncated squared return. The same bets on the truncated fourth power `Q_t = min(r^4, tau^4)` target the truncated local quarticity of Bapat (2026). Is that valid, how does it compare on pure variance changes, does it see changes of shape (tails, vol of vol) that the variance stream cannot, and what do its alarms mean on real data?

**Method.** Two things had to be derived: the null mean `q0` (the model-free sample mean of `min(r^4, tau^4)` over the calibration window; the Gaussian shortcut `3 m0^2` is 2.7% too high at `c = 3` and 21% too high with the raw MAD) and the decrease bet range `lambda < 1 / (tau^4 / q0 - 1) = 0.037` for normal returns, against 0.124 for the variance stream. With those, the proofs of `docs/math.md` section 5 apply unchanged. Settings as in E4: `tau = 3 x scaled MAD` from a 2000-step calibration window, `A = 7800`, 300 paths. A theory function computes the mean shift and signal-to-noise ratio of each stream under any alternative density.

**Q1 Validity [SYNTHETIC].** No change, 20,000 steps. False alarm within 390 steps (bound 5%) and censored mean run length (bound 7800): Gaussian 1.7% / 9471 (up) and 0.0% / 9046 (down); jumps 1.7% / 9348 and 0.0% / 9334; Student t5 2.7% / 8858 and 0.0% / 8558. Both bounds hold in every cell, for both streams.

**Q2 Pure variance change [SYNTHETIC].** Median delay in steps, Gaussian returns (jumps in brackets):

| sd ratio | x 1.25 | x 1.5 | x 2 | x 0.8 | x 0.67 | x 0.5 |
|---|---|---|---|---|---|---|
| variance stream | 96 (92) | 38 (40) | 18 (18) | 153 (151) | 95 (99) | 68.5 (74) |
| quarticity stream | 87.5 (91) | 33 (37) | 13 (14) | 271 (299) | 201.5 (216) | 168 (189) |

**Q3 Shape changed, variance fixed [SYNTHETIC].** 20,000 steps after the change, `A = 400,000` (null false alarm within the horizon at most 5%), 200 paths. Share detected (median delay):

| Alternative | variance stream, up | variance stream, down | quarticity stream, up |
|---|---|---|---|
| Student t3 / t5 / t8 | 0% / 0% / 3% | 100% (650) / 96.5% (3937) / 57.5% (8448) | 6% / 71.5% (6461) / 75% (7626) |
| two-point mixture, d = 0.5 / 0.7 / 0.9 | 8% / 5% / 2% | 28% / 46% / 58% | 82.5% (4915) / 99.5% (2678) / 99% (1033) |
| slow SV, log-sd 0.3 / 0.5 | 100% (715) / 100% (496) | 100% (524) / 100% (271) | 100% (664) / 100% (442) |

**Q4 Real data [REAL].** `A = 28,800` per side, 1440-minute calibration, recalibration after each alarm, band 0.25. Alarms per day: variance monitor 0.69 (BTC) and 0.69 (DOGE); quarticity monitor 0.71 and 0.70. Share of quarticity alarms on the same UTC day as a variance alarm: 87% and 87%; within one hour: 62% and 60%. Median truncated variance in the next 24 hours over the baseline: 1.45 and 1.38 after up alarms, 0.59 and 0.66 after down alarms (control at random times: 1.00). The 34 to 35 quarticity alarms per coin that are more than four hours from any variance alarm are still followed by a move in the flagged direction 80% and 88% of the time. All monitors on both coins alarmed at 15:28 to 15:31 UTC on 10 October 2025, the crypto crash after the US tariff announcement (checked online, `results/q4_real_event_check.csv`). With the same variance tolerance (band 0.5625 on quarticity), the quarticity monitor alarms 20% less often and its decrease alarms are followed by lower variance 100% and 97% of the time.

**Reading.**
* Valid, as the theory says.
* On pure variance changes the quarticity stream is 1% to 28% *faster* for increases (the fourth power amplifies the shift more than it adds noise) and 1.8 to 2.6 times slower for decreases (the bet cap).
* It detects vol-of-vol changes at fixed variance that the variance stream cannot (82% to 99% against 2% to 8%), but it does *not* beat the variance stream on fat tails at `c = 3`: clipping at 3 sd removes most of a fat tail's fourth moment, while the extra mass near zero lowers the truncated variance, so a Student t registers on the variance stream as a volatility *decrease*. At `c = 5` the theory reverses the ranking. The clipping level decides which shape changes are visible.
* On real data the two monitors behave alike; no pure shape event was found in six months of two coins. A planned shape diagnostic (`q_ratio / var_ratio^2`) measured only the mechanics of a fixed `tau` and is reported as such.

**Caveats.** Shape detection needs thousands of steps (signal-to-noise 0.03 to 0.10 per step). A first version of Q3 with `A = 7800` was unreadable because the null itself alarms within 10,000 steps on half the paths; the fix is recorded in `docs/iterations/ALT1.md`.

---

## A1. Pre-averaging before MAD truncation: microstructure noise at one-second sampling

**Question.** Under the standard price-noise model, at what noise-to-signal ratio `gamma = omega / (sigma sqrt(delta))` does the plain truncated estimator break, for integrated quarticity and for the monitor's local variance, and does pre-averaging (Jacod, Li, Mykland, Podolskij and Vetter 2009) followed by truncation fix it? What happens to the monitor's false alarm rate under noise, with and without pre-averaging, given that pre-averaged returns overlap?

**Method [SYNTHETIC].** One-second returns, 23,400 per day, daily volatility 1%, log-OU stochastic volatility, optional 4 jumps per day of 50 per-second sd, noise sd `gamma` times the per-second return sd. Pre-averaging window `k = 60`, weight `min(x, 1 - x)`; the local pre-averaged variance `V_i = rbar_i^2 / (k psi2 delta)` has mean `sigma^2 + c` with `c = omega^2 psi1 / (psi2 k^2 delta)`, and the local quarticity is corrected as `Q_i = rbar_i^4 / (3 (k psi2 delta)^2) - 2 c V_i + c^2`, whose mean is `sigma^4`; `omega^2` is minus the first autocovariance of clipped raw returns. Truncation of `rbar` at `3 x scaled MAD` by clipping (the paper's rule) or dropping, each with its Gaussian survival factors. Monitors run at `A = steps per day / 0.05`.

**Q1 Integrated quarticity, mean ratio to the truth (relative RMSE) [SYNTHETIC].**

| gamma | 0 | 0.25 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|
| raw returns, MAD clip, no jumps | 0.99 (0.02) | 1.26 (0.26) | 2.23 (1.24) | 9.0 (8.1) | 81.4 (81.5) | 1138 (1155) |
| pre-averaged, MAD clip, no jumps | 0.99 (0.11) | 1.00 (0.11) | 1.00 (0.12) | 0.99 (0.10) | 0.99 (0.11) | 0.99 (0.11) |
| pre-averaged, MAD clip, no noise correction | 0.99 (0.11) | 1.00 (0.11) | 1.00 (0.12) | 1.00 (0.10) | 1.01 (0.12) | 1.10 (0.15) |
| pre-averaged, MAD clip, 4 jumps | 1.14 (0.22) | 1.15 (0.22) | 1.14 (0.20) | 1.17 (0.22) | 1.18 (0.25) | 1.15 (0.24) |
| pre-averaged, MAD drop, 4 jumps | 0.99 (0.12) | 1.00 (0.13) | 1.00 (0.11) | 1.01 (0.12) | 1.02 (0.12) | 0.99 (0.12) |
| pre-averaged, no truncation, 4 jumps | 33.4 (69) | 36.5 (68) | 27.0 (47) | 35.1 (63) | 39.3 (71) | 32.5 (52) |

The raw estimator follows `(1 + 2 gamma^2)^2` exactly: 10% bias at `gamma = 0.16`, a factor 2 at 0.45. Integrated variance follows `1 + 2 gamma^2` (10% at 0.22). The monitor's 300-second local variance on raw returns is inflated by the same `1 + 2 gamma^2` (1.51 at `gamma = 0.5`, 3.07 at 1), while the pre-averaged local variance has mean ratio 0.98 to 1.07 at every noise level but a typical error of 32% (overlapping) to 45% (blocks) against 6% for raw returns without noise.

**Q3 One-day false alarm rate of the monitor, increase side (bound 5%) [SYNTHETIC].**

| gamma | 0 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| raw returns (23,400 steps) | 5.0% | 8.5% | **18.0%** | **31.0%** | **41.5%** |
| overlapping pre-averaged, k = 60 | **100%** | **100%** | **100%** | **100%** | **100%** |
| non-overlapping blocks, k = 60 / 20 / 10 | 2.5% / 4.5% / 5.5% | 1.5% / 3.0% / 5.0% | 2.5% / 2.5% / 5.5% | 2.0% / 5.5% / 2.5% | 4.0% / 3.0% / 3.0% |

**Q4 Median detection delay in seconds, sd doubles at the start of the day [SYNTHETIC].** Raw returns 32, 43.5, 86.5, 449.5, 2907 at `gamma = 0, 0.5, 1, 2, 4`; blocks with `k = 10`: 259, 269, 299, 359, 649; `k = 60`: 1319 to 1379 throughout. For a halving at `gamma = 4` the raw monitor detects on 70% of days (median 10,795 s); the `k = 10` blocks on 100% (3939 s).

**Reading.**
* Truncation alone does nothing against noise; the plain estimator is off by 26% at `gamma = 0.25` and useless by `gamma = 1`. Pre-averaging then truncation is unbiased at every noise level, at the cost of a relative RMSE of 0.11 instead of 0.02 when there is no noise. Switch it on when the estimated `gamma` exceeds about 0.15.
* After pre-averaging, clip must become drop: a clipped jump keeps `tau^4` in each of the `k - 1` windows it enters (+12% to +18% with 4 jumps per day, growing with `k`).
* Overlapping pre-averaged returns break the e-detector completely: the conditional-mean null fails because consecutive windows share `k - 2` returns. Blocks on disjoint prices restore it at every `k`.
* **The raw-return monitor also loses its guarantee under price noise** (8.5% at `gamma = 0.5`, 41.5% at 4). With price noise `E[r_t^2 | past] = sigma^2 delta + omega^2 + E[u_{t-1}^2 | past]`, so squared returns are positively autocorrelated and the conditional null fails. E4 of the main report saw 3.7% at `gamma = 0.5` only because it looked at 390 steps. The remark in `docs/math.md` 5.1 that the null allows any dependence should be qualified: only dependence that keeps the *conditional* moment on the null side.
* The practical monitor at one second is the `k = 10` block monitor: valid at every noise level, 8 times slower than raw returns without noise, faster than them (and valid, which they are not) from `gamma` of about 2.

**Caveats.** Constant volatility in Q3 and Q4 (so the noise effect is isolated); i.i.d. Gaussian noise, which real noise is not (Hansen and Lunde 2006); the fourth-moment correction assumes local Gaussianity of the pre-averaged return; no real one-second data was available in this environment.

## A1 follow-up (B69). Noise and jumps together, the published benchmark, and the block monitor's guarantee

**Question.** When microstructure noise and jumps are present at the same time, how does pre-averaging then MAD truncation compare with the estimator built for exactly that case, the modulated bipower and multipower variation of Podolskij and Vetter (2009, *Bernoulli* 15(3), 634 to 658), for integrated quarticity and variance? And for the monitor: what exactly is the guarantee of the non-overlapping block version under price noise, and how should the block size `k` be chosen?

**Method [SYNTHETIC].** `alternatives/preaveraging/run_benchmark.py` (two processes, about 2 minutes). Q5: one-second returns, 23,400 per day, log-OU stochastic volatility, 200 paths per cell, noise ratio `gamma` in `{0, 0.5, 1, 2}` x jumps per day and size in `{none, 2 x 25 sd, 2 x 100 sd, 8 x 25 sd, 8 x 100 sd}` (25 per-second sd is 0.16% of the price). Estimators: raw MAD clip (Bapat 2026); pre-averaged (`k = 60`) MAD drop (ours) and MAD clip; modulated bipower (variance) and tripower (quarticity) variation with `K = 76`-step increments averaged over 152 blocks of 153 observations (`c1 = 0.5`, `c2 = 1`), noise-corrected through the exact block weights; pre-averaged bipower and tripower; pre-averaged without truncation. Q6: constant volatility, `gamma` in `{0.5, 1, 2, 4}`, blocks `k` in `{4, 6, 10, 14, 20, 30, 60}`, `A = (n / k) / 0.05`, sd `x 2`, `x 1.5` and no change, 200 paths per cell.

**Reference check.** Title, authors, journal, volume and pages of Podolskij and Vetter (2009) and the arXiv listing 0909.0827 were confirmed from search listings; the full text could not be fetched (arxiv.org, projecteuclid.org, repec.econ.au.dk and eldorado.tu-dortmund.de are blocked here). The block-average variance constants `nu1 = c1^2 (3 c2 - 4 c1) / (3 (c2 - c1)^2)` and `nu2 = 2 c1 / (c2 - c1)^2`, and the CLT standard error from the independent Gaussian block approximation, are therefore our own derivation, verified by simulation in the tests; no equation numbers are quoted.

**The negative finding, in one table. One-day false alarm rate of the monitor under price noise, bound 5%, 200 days per cell [SYNTHETIC].**

| monitor (steps per day, `A`) | side | gamma = 0 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|
| raw truncated returns (23,400, 468,000) | up | 5.0% | **8.5%** | **18.0%** | **31.0%** | **41.5%** |
| raw truncated returns | down | 7.5% | **9.0%** | **13.5%** | **29.0%** | **27.5%** |
| overlapping pre-averaged, k = 60 (23,341, 466,820) | up and down | **100%** | **100%** | **100%** | **100%** | **100%** |
| non-overlapping blocks, k = 60 (390, 7800) | up | 2.5% | 1.5% | 2.5% | 2.0% | 4.0% |
| non-overlapping blocks, k = 60 | down | 0.0% | 0.5% | 1.0% | 2.0% | 1.5% |
| non-overlapping blocks, k = 20 (1170, 23,400) | up | 4.5% | 3.0% | 2.5% | 5.5% | 3.0% |
| non-overlapping blocks, k = 20 | down | 3.0% | 5.0% | 1.0% | 3.5% | 2.5% |
| non-overlapping blocks, k = 10 (2340, 46,800) | up | 5.5% | 5.0% | 5.5% | 2.5% | 3.0% |
| non-overlapping blocks, k = 10 | down | 1.5% | 4.0% | 2.0% | 2.0% | 4.5% |
| non-overlapping blocks, k = 4 (5850, 117,000), Q6 | up | n/a | 4.5% | 4.5% | 3.5% | 4.5% |

Monte Carlo sd of each entry is about 1.5 points at a true 5%. The raw-return monitor's rate rises monotonically with `gamma` on both sides; the block monitors stay within the bound in all 58 cells (28 of Q6, 30 of Q3, with the four entries above 5% within 1.5 points of it).

**Why: the conditional null under MA(1) noise, and what blocks restore.** With observed prices `Y_j = X_j + u_j` (noise i.i.d., mean zero, variance `omega^2`, any distribution) the raw return is `r_t = Delta X_t + u_t - u_{t-1}` and, with `F_{t-1}` containing `Y_{t-1}`, `E[r_t^2 | F_{t-1}] = sigma^2 delta + omega^2 + E[u_{t-1}^2 | F_{t-1}]`. The last term has mean `omega^2` but is not constant (for Gaussian noise and the previous return alone, `E[r_t^2 | r_{t-1}] = sigma^2 delta + 2 omega^2 + rho_u^2 (r_{t-1}^2 - sigma^2 delta - 2 omega^2)` with `rho_u = gamma^2 / (1 + 2 gamma^2)`), so after a large squared return the conditional mean exceeds the baseline with no change in volatility: the null `E[Y_t | F_{t-1}] <= 1` of the e-detector fails, and so does the bound. The block monitor uses `Z_b = min(rbar_b^2, tau^2)` with `rbar_b` pre-averaged over the prices `Y_{bk}, ..., Y_{bk+k-1}`; the return between `Y_{bk+k-1}` and `Y_{(b+1)k}` is skipped, so block `b+1` shares no price with block `b`, its noise values are fresh, and given everything up to the end of block `b` it has the law of a pre-averaged return of a process started there. The independence is exact, not up to a boundary term (if adjacent blocks shared one price the shared noise value would give `Cov(rbar_b, rbar_{b+1}) = -omega^2 / k^2`, an `O(1/k)` relative dependence against the noise part `omega^2 / k`; skipping one return removes it). The null `H0: E[Z_{b+1} | G_b] <= m0`, with `m0` the calibration mean that includes the noise term `psi1 omega^2 / k`, then holds whenever the spot variance over block `b+1` is at most its calibration value and the noise and jump laws are unchanged, by monotonicity of `min(s^2 e^2, tau^2)` in `s`. The betting factor `1 + lambda (Z_{b+1} / m0 - 1)` with a predictable `lambda in [0, 1)` has conditional mean at most 1, so the Shiryaev-Roberts statistic over blocks is an e-detector with respect to the block filtration and Theorems A and B of `docs/math.md` apply with `m = floor(n / k)` blocks per day and `A = m / 0.05`. Nothing requires Gaussian blocks, so the guarantee holds for every `k >= 2`, which Q6 confirms down to `k = 4`. Full statement in `alternatives/preaveraging/README.md`, section 2.

**Q5 Integrated quarticity with noise and jumps, mean ratio to the truth (relative RMSE) [SYNTHETIC].**

| cell | raw, MAD clip | pre-averaged, MAD drop (ours) | modulated tripower (PV 2009) | pre-averaged tripower |
|---|---|---|---|---|
| no jumps, gamma = 0 / 1 / 2 | 0.99 (0.02) / 8.9 (8.0) / 81 (81) | 0.98 (0.11) / 0.98 (0.11) / 0.98 (0.11) | 1.00 (0.30) / 0.99 (0.33) / 1.00 (0.30) | 1.00 (0.13) / 1.00 (0.13) / 1.00 (0.12) |
| 2 x 25 sd, gamma = 0 / 1 | 0.99 (0.02) / 9.0 (8.1) | 0.99 (0.10) / 0.98 (0.11) | 1.12 (0.37) / 1.03 (0.34) | 1.06 (0.15) / 1.04 (0.16) |
| 2 x 100 sd, gamma = 0 / 1 | 0.99 (0.02) / 9.1 (8.2) | 0.99 (0.12) / 0.98 (0.11) | 1.63 (1.27) / 1.57 (1.73) | 1.45 (0.86) / 1.56 (1.63) |
| 8 x 25 sd, gamma = 0 / 1 | 1.00 (0.02) / 9.1 (8.2) | 1.02 (0.11) / 1.05 (0.12) | 1.27 (0.57) / 1.35 (0.65) | 1.21 (0.33) / 1.25 (0.34) |
| 8 x 100 sd, gamma = 0 / 1 / 2 | 1.00 (0.02) / 9.1 (8.1) / 81 (81) | 1.00 (0.11) / 1.01 (0.12) / 1.01 (0.12) | 4.9 (10) / 5.8 (14) / 4.5 (6.9) | 4.0 (7.8) / 3.6 (4.8) / 3.8 (5.2) |

Integrated variance: ours 0.99 to 1.02 (RMSE 0.05 to 0.06) in all 20 cells; modulated bipower 0.99 to 1.01 (0.13 to 0.14) without jumps, 1.05 with `2 x 25 sd`, 1.10 with `8 x 25 sd`, 1.16 to 1.20 with `2 x 100 sd`, 1.72 to 1.83 with `8 x 100 sd`. Coverage of the modulated bipower 95% CLT interval: 0.905 to 0.955 without jumps (half-width 26% of the truth), 0.905 to 0.97 with 25 sd or two 100 sd jumps, 0.57 to 0.62 with `8 x 100 sd`.

**Q6 Median detection delay in seconds for a doubling of the sd, blocks of size `k` [SYNTHETIC].**

| gamma | k = 4 | 6 | 10 | 14 | 20 | 30 | 60 | best observed | `k* = 6 gamma / 2` | rule `max(4, 4 gamma)` |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | **127** | 179 | 269 | 370 | 479 | 749 | 1319 | 4 | 1.5 | 4 |
| 1 | **179** | 206 | 289 | 377 | 539 | 719 | 1439 | 4 | 3 | 4 |
| 2 | 405 | **332** | 354 | 419 | 539 | 764 | 1439 | 6 | 6 | 8 |
| 4 | 2131 | 1154 | 689 | **615** | 659 | 869 | 1439 | 14 | 12 | 16 |

For `x 1.5` the best observed `k` is 4, 4, 10, 20 against `k* = 6 gamma / 1.5` = 2, 4, 8, 16 (delays 299, 459, 839, 1539 s). The delay formula: after a sd change by `v` the block value has mean `(v^2 + rho_k) / (1 + rho_k)` with noise share `rho_k = 12 gamma^2 / k^2`, the bet grows by about `snr^2 / 2` per block with `snr = (v^2 - 1) / (sqrt(2) (v^2 + rho_k))`, so the delay in seconds is proportional to `k (v^2 + rho_k)^2`, minimised at `rho_k = v^2 / 3`, that is `k* = 6 gamma / v`. The curve is flat within 10% for `k` between `k* / 1.5` and `1.5 k*`, so the single rule `k = max(4, ceil(4 gamma))` (target `v = 1.5`, `gamma` estimated from `omega^2 = -mean(r_t r_{t-1})` over `mean(r_t^2) - 2 omega^2`) is within 10% of the best observed delay in all 8 cells.

**Reading.**
* The raw-return monitor's loss of validity under price noise is real, monotone in `gamma` and two-sided; it is not a Monte Carlo artefact (41.5% against 5%). It is the one claim of the main package that one-second data breaks. The fix is not pre-averaging as such (overlapping pre-averaged returns give 100% false alarms) but non-overlapping blocks, whose guarantee is exact for every `k >= 2`.
* Against the published noise-and-jump estimator, pre-averaging then MAD drop is unbiased over the whole grid and 2.5 to 3 times more precise without jumps; the modulated multipower estimators keep `|J|^{4/3}` per affected product at `n = 23,400`, a factor 4.5 to 5.8 in quarticity with eight 0.65% jumps, and its CLT interval then covers 60%. Noise and jumps do not interact for any estimator: pre-averaging handles one, truncation the other, and the two compose.
* The block size is a solved trade-off: `k* = 6 gamma / v`, in practice `k = max(4, 4 gamma)`. At `gamma <= 1` the smallest block (`k = 4`) is best and detects a doubling in 127 to 179 s, a factor 3 to 4 behind the raw monitor's 43.5 to 86.5 s, which is invalid there; the previously reported factor 8 came from `k = 10`.

**Needed changes outside ALT1's files.** `docs/math.md` 5.1 ("the null allows any dependence") and 5.6 should say that MA(1) price noise breaks the conditional null on raw returns and that the block monitor restores it; `ALTERNATIVES.md` A1 can record the benchmark outcome.

**Caveats.** All numbers are synthetic; i.i.d. Gaussian price noise (real noise is autocorrelated and heteroskedastic, Hansen and Lunde 2006); constant volatility in Q6; the delay formula ignores truncation; the modulated bipower constants and interval are our derivation, not the paper's printed ones; no one-second real data was available to estimate actual `gamma` values.
