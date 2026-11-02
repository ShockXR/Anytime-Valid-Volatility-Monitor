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
