# Two-sample betting monitor: conformal rank e-detector (alternative A4)

A minimal working prototype of entry A4 of `alternatives/ALTERNATIVES.md`: instead of betting that the mean of the truncated squared return has moved away from a baseline estimated on a calibration window, bet that **the newest return no longer looks like the returns seen before it**. The test compares the recent stream with the calibration window directly, has no baseline parameter to estimate, and reacts to any change in the distribution of returns (variance, tails, jump intensity and, with a signed score, skew). Every number quoted from experiments is in `RESULTS.md`, labelled [SYNTHETIC] or [REAL], and comes from `run_experiment.py`.

## 1. Research problem

The main monitor (`src/avmonitor/monitor.py`, `docs/math.md` section 5) watches one number: the mean of `X_t = min(r_t^2, tau^2)`, the truncated squared return. It has two limits that this prototype attacks.

1. **It needs a baseline.** Its guarantee holds when it bets against the TRUE pre-change mean `m0`. In practice `m0` is estimated from a calibration window, so the guarantee is only approximate (plug-in), or it needs a margin, or a confidence sequence for `m0` (`src/avmonitor/confseq.py`), which costs some power.
2. **It sees only one moment.** A change in the *shape* of the return distribution at fixed variance (fatter tails, more frequent jumps, a skew towards large down moves) moves the truncated variance only a little, sometimes in a misleading direction (a variance-matched fat-tailed law has a *lower* truncated variance).

The research question: can a sequential two-sample test by betting (Shekhar and Ramdas 2024) be turned into a change detector for intraday returns that (a) keeps an exact, finite-sample false alarm guarantee with **no estimated baseline at all**, (b) detects any change in the distribution of returns, and (c) costs little against the mean-based monitor when the change *is* a variance change? And what does it alarm on in real one-minute crypto data that the mean-based monitor does not?

## 2. Mathematical formulation

### 2.1 Set-up

`u_1, u_2, ...` are standardised returns (on real data: one-minute log returns divided by the predictable intraday volatility profile `avmonitor.data.intraday_scale`). A *segment* starts after each alarm. Its first `n0` returns form the reference (calibration) window; every later return is monitored. A *score* `s_t = s(u_t)` is fixed in advance:

* `kind="abs"`: `s_t = |u_t|`. Sees changes in the law of the size of moves (variance, tails, jumps, share of zero returns).
* `kind="signed"`: `s_t = u_t`. Also sees changes in the sign pattern (skew).

No scale, threshold, mean or quantile is estimated.

### 2.2 The null hypothesis

**H0 (no change):** the returns of the segment, reference window included, are *exchangeable*: their joint law does not change when their order is permuted. I.i.d. returns with any law (Gaussian, Student t, jumps, discrete prices with many exact zeros) are exchangeable. The null says nothing about the shape of that law.

This is a different trade from the main monitor. Its null `E[X_t | past] <= m0` restricts one conditional moment and allows any dependence, but needs `m0`. H0 here needs no parameter but forbids dependence: the order of the returns must carry no information. Volatility clustering (GARCH-type persistence) violates H0 even when the volatility level is stable (section 2.7).

### 2.3 Conformal p-values

At monitored time t, let `n = t - 1` be the number of earlier returns of the segment, `G_t` the number of them with a larger score, `E_t - 1` the number with an equal score, and `U_t` an independent uniform draw on (0, 1). The *smoothed conformal p-value* is

`p_t = (G_t + U_t E_t) / (n + 1)`.

A small `p_t` means "the new score is larger than almost everything seen in this segment". The reference sample *grows*: return t is ranked against the calibration window and against all monitored returns before it.

**Lemma 1 (Vovk, Nouretdinov and Gammerman 2003; Vovk, Gammerman and Shafer 2022, chapter 2).** Under H0, `p_1, p_2, ...` (computed over the whole segment) are independent and uniform on (0, 1), whatever the law of the returns, ties included.

*Proof sketch.* Let `B_t` be the unordered bag of the first t scores and `F^-_t = sigma(B_t, s_{t+1}, s_{t+2}, ..., U_{t+1}, U_{t+2}, ...)` (a decreasing filtration). Given `F^-_t`, exchangeability makes `s_t` a uniformly random element of the bag `B_t`, and `U_t` is independent of everything else, so `p_t`, which is the randomised position of `s_t` in the sorted bag, is uniform on (0, 1) given `F^-_t` (the randomisation spreads each tied block evenly over its share of (0, 1)). Also `p_{t+1}, p_{t+2}, ...` are `F^-_t`-measurable (each is a function of a later bag, a later score and a later `U`). Hence `p_t` is uniform and independent of `(p_{t+1}, p_{t+2}, ...)` for every t, which gives independence of the whole sequence by induction. ∎

The p-values are therefore an exactly known null law for *any* return distribution: this is what replaces the baseline `m0`. The reference window is not needed for validity (Lemma 1 holds from the first return of the segment); it only makes the early p-values fine enough to bet on.

**Why the reference must grow.** The tempting simpler design ranks each new return against the fixed reference window only (equivalently, "quantile bins learned from the reference"). Given the reference, those ranks are i.i.d. with bin probabilities equal to the random spacings of the reference quantiles (Dirichlet with all parameters 1 for continuous laws), not the uniform ones. Unconditionally each p-value is uniform, but they are dependent, and a bettor that learns the bin frequencies grows its wealth at the rate `KL(spacings, uniform)`, about `(K - 1) / (2 n0)` per step for K bins, under the null. A false alarm is then certain in the long run. This is the nonparametric form of the plug-in baseline problem of the main monitor, and Q1c measures it. Ranking against the growing sample (a Polya-urn update of the reference) is what makes the bets exactly fair.

### 2.4 Bets and the e-detector

Let `G_{t-1} = sigma(p_1, ..., p_{t-1})`. A *bet* at time t is a function `f_t(p) = 1 + lam_t (w_t(p) - b E_U|w_t|)`, where

* `w_t` is `G_{t-1}`-measurable (it may depend on past p-values only), has `E_U[w_t] = int_0^1 w_t(p) dp = 0`, and is bounded below, `w_t >= -L_t`;
* `0 <= lam_t < 1 / (L_t + b E_U|w_t|)`, so that `f_t > 0`;
* `b >= 0` is the null band of section 2.5 (`b = 0` for the exact null).

The prototype mixes three families of witnesses `w`:

1. **Scale up and scale down (the rank analogue of the main monitor).** With `z_t` the standard normal quantile that matches `p_t` (`z = Phi^{-1}(1 - p/2)` for the absolute score, `Phi^{-1}(1 - p)` for the signed one), define the Gaussianised score `Y_t = min(z_t^2, a^2) / kappa_2(a)`, with `a = 3` and `kappa_2(a) = E[min(Z^2, a^2)] = 0.9946` (`docs/math.md`, section 3). Under H0, `z_t` is exactly standard normal, so **`E[Y_t] = 1` exactly**. The bets are `w = Y - 1` (increase, `L = 1`) and `w = 1 - Y` (decrease, `L = B - 1` with `B = a^2 / kappa_2(a) = 9.05`), with the same grids as the main package (`DEFAULT_UP_GRID`, `DEFAULT_DOWN_FRACTIONS`). This is the main monitor's statistic computed on ranks: a truncated-variance bet with a baseline that is *known* rather than estimated.
2. **Histogram learner (the omnibus witness).** Fix bins of [0, 1] with widths `q_k` (finer near 0 and 1: edges 0, 0.002, 0.01, 0.03, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97, 0.99, 0.998, 1). Keep exponentially weighted counts `c_k` of past p-values in each bin (half-lives 240 and 1440 steps, one learner each), and smooth them towards the uniform law with 50 pseudo-observations: `pi_k = (c_k + 50 q_k) / (sum_j c_j + 50)`. The witness is `w(p) = pi_k / q_k - 1` for p in bin k (`L = 1`), the bet uses `lam = (1 - eps) / (1 + b E_U|w|)` with `eps = 0.05`. For `b = 0` this is the plug-in density bettor `f = (1 - eps) pi_k / q_k + eps` of Fedorova, Gammerman, Nouretdinov and Vovk (2012): it learns where the new p-values pile up and bets on that region, whatever its shape.

Each bettor `i` runs a Shiryaev-Roberts recursion `S^i_t = f^i_t(p_t) (S^i_{t-1} + 1)`, `S^i_0 = 0`, and the detector statistic is the weighted average `M_t = sum_i v_i S^i_t` (weights: 0.25 on the scale-up grid, 0.25 on the scale-down grid, 0.5 on the two histogram learners, equal within a family). The alarm is `tau_A = inf{t : M_t >= A}`.

**Theorem 1.** Under H0, for every A and m: `E[tau_A] >= A` and `P(tau_A <= m) <= m / A`. With `b > 0` the same holds under H0_b (section 2.5).

*Proof.* By Lemma 1, `p_t` is uniform and independent of `G_{t-1}`, and each `w^i_t` is `G_{t-1}`-measurable with mean 0 under the uniform law, so `E[f^i_t(p_t) | G_{t-1}] = 1 + lam^i_t (0 - b E_U|w^i_t|) <= 1`, with equality when `b = 0`. Then `E[S^i_t | G_{t-1}] <= S^i_{t-1} + 1`, hence `E[M_t - t | G_{t-1}] <= M_{t-1} - (t - 1)` since the weights sum to 1. This is Lemma 1 of `docs/math.md`, section 5.4, in the filtration of the p-values, and `tau_A` is a stopping time of that filtration. Theorems A (optional stopping) and B (Ville's inequality applied to `M_t + (m - t)`) of that section apply word for word. ∎

With `A = 1440 / 0.05 = 28,800` and 1440-minute days, the chance of any false alarm within a day of monitoring is at most 5%, the same budget as the main monitor (which runs each side at `2A`). No split between an up and a down side is needed: every direction lives inside one mixture. Under exact H0 the statistic is a martingale (`E[M_t] = t`), not just a supermartingale, so the bound is not loosened by the bets themselves.

### 2.5 A null band for mild dependence

**H0_b:** for every t, the conditional law of `p_t` given `G_{t-1}` has a density `h_t` with `1 - b <= h_t(p) <= 1 + b`. Exact H0 is `b = 0`. For a witness with `E_U[w] = 0`, `E_h[w] = int w (h - 1) <= b int |w| = b E_U|w|`, so the penalised bet `1 + lam (w - b E_U|w|)` again has conditional mean at most 1, and Theorem 1 holds. In words: "the rank of a new return may be up to `100 b`% more or less likely to fall in any region than under exchangeability". The band also bounds what the scale bets can be pushed by, because `|h - 1| <= b` implies `|E[Y] - 1| <= b E_U|Y - 1| = 0.966 b`. The prototype uses `b = 0` (exact) and `b = 0.25` (the same number as the main monitor's `null_band`, fixed before any real-data run).

### 2.6 What truncation does when the statistic is a rank

Three separate statements.

1. **The raw rank needs no truncation.** For any strictly increasing map g, `g(|u|)` has the same ranks as `|u|`. Clipping at a FIXED level tau is not strictly increasing above tau: all returns beyond tau become ties, and the randomised p-value spreads them uniformly over their block. Lemma 1 still holds (ties are allowed), but the test can no longer tell a 4 sd move from a 40 sd move; it only counts them. If tau were estimated from the ranked sample itself (for example `3 x MAD` of the calibration window), the clipped scores would no longer be exchangeable and Lemma 1 would hold only approximately. So for ranks, truncation is never needed for validity and can only discard information. Ranks are already bounded, and one return moves one rank, whatever its size: the bounded influence that clipping buys the main monitor (Bapat's Theorem 3 in betting form) comes for free. Q3 checks that clipping at a fixed `tau = 3` changes delays and false alarms only within noise.
2. **The Gaussianised score does need truncation, for decreases, exactly as in the main monitor.** `z^2` is unbounded above, and the argument of `docs/math.md` 5.3 shows that no nonparametric evidence that "the scale has not fallen" can build up from an unbounded statistic. Clipping `z^2` at `a^2` gives the bound `B = 9.05` and the decrease bets. The difference from the main monitor is that `a` is a fixed constant and the null mean `kappa_2(a)` is known exactly, so truncation no longer brings an estimation problem with it.
3. **Bounded influence.** One return changes each scale-up factor by at most `1 + 0.5 (B - 1) = 5.0`, each scale-down factor by at most `1 + 0.95 = 1.95`, and each histogram factor by at most `(1 - eps) max_k pi_k / q_k + eps`.

### 2.7 Where volatility clustering breaks the null

Exchangeability forbids any serial dependence in the scores. With GARCH-type clustering, a large move makes the next large move more likely, so `p_t` given the past is *not* uniform: after a burst of small p-values, small p-values are more likely. Every bettor that learns from recent p-values, and the scale bets inside a Shiryaev-Roberts recursion that restarts every step, can exploit that and alarm. This is the same mechanism that makes the main monitor alarm on persistent clusters (E9, Q9c of the main report): neither design can tell "a persistent volatility cluster" from "a change". The band `b` absorbs short-lived dependence only. A cleaner fix would rank returns divided by a predictable volatility forecast (the forecast-relative idea of E10), so that the null becomes "exchangeable relative to the forecast".

### 2.8 Relation to two-sample testing by betting

In Shekhar and Ramdas (2024), a bettor observes pairs `(X_t, Y_t)` from two streams and bets with payoff `lam_t (g_t(Y_t) - g_t(X_t))` for a predictable witness `g_t` (for example a plug-in kernel MMD witness); under `P_X = P_Y` the payoff has mean 0. Podkopaev and Ramdas (2023) use a sequentially trained classifier as the witness, and Podkopaev, Bloebaum, Kasiviswanathan and Ramdas (2023) apply the same principle to independence testing. Our problem has one fixed reference sample and one stream. A literal paired design would draw one unused reference return per step: it is valid, but it consumes the reference (it stops after `n0` steps) and it doubles the noise of each payoff, because both `g(Y_t)` and `g(X_t)` are random. The conformal construction uses the whole (growing) reference at every step and turns the comparison into a single p-value with an exactly uniform null law. The histogram learner then plays the role of the plug-in witness: `pi / q - 1` estimates the density ratio of the new data to the reference in rank space, the one-dimensional analogue of the MMD witness. Conformal versions of the Shiryaev-Roberts and CUSUM procedures were introduced by Vovk (2021) and Vovk et al. (2021); the e-detector bounds used here (average run length at least A and the per-horizon bound) are those of Shin, Ramdas and Rinaldo (2024), applied in the filtration of the p-values.

## 3. Why it matters

* **An exact guarantee without a baseline.** The main package needs a plug-in, a margin or a confidence sequence for `m0`. Here the null law of the bet input is exactly uniform by construction, so the guarantee is exact for any calibration length, even 2 hours of data (Q1).
* **One detector for every kind of change.** Variance, tails, jump intensity and skew enter one mixture. The family that carries the evidence at the alarm (scale up, scale down or histogram) says what kind of change it was.
* **A sharper view of what truncation does.** For ranks, truncation is unnecessary for validity and for bounded influence; it remains necessary only to bound the Gaussianised decrease bet, now at a fixed level (section 2.6).
* **A trap made visible.** The obvious "learn quantile bins on the reference, then bet" design is not a valid anytime test (section 2.3, Q1c), for the same reason a plug-in `m0` is not.
* **Clear limits.** Exchangeability is a stronger null than a conditional moment bound. On real one-minute returns the test reacts to things the mean-based monitor ignores (Q4), including changes that are not volatility changes at all.

## 4. Comparison with existing approaches

| Approach | Null | Baseline to estimate | Detects | Guarantee |
|---|---|---|---|---|
| Main monitor, plug-in (`VolatilityMonitor`) | `E[X_t | past] = m0` (within a band) | `m0`, `tau` | Changes in truncated variance | Exact only with the true `m0`; plug-in approximately |
| Main monitor, confidence-sequence baseline (`baseline="cs"`) | Same | Confidence bounds for `m0` | Same | `delta + m / A` with an estimated baseline |
| Quarticity stream (A3, `alternatives/quarticity_stream`) | `E[min(r^4, tau^4) | past] = q0` | `q0`, `tau` | Truncated fourth moment (tails, vol of vol) | As the main monitor |
| **This prototype (conformal rank e-detector)** | Exchangeability of the segment | **None** | Any change in the law of `|u|` (or of `u`) | **Exact**, `m / A` per horizon, ARL at least A |
| Paired two-sample betting (Shekhar and Ramdas 2024) with a fixed reference | `P_X = P_Y`, pairs i.i.d. | None | Any change the witness sees | Exact while unused reference points remain |
| Ranks against a fixed reference (quantile bins of the calibration window) | Not a martingale under exchangeability | None | Any change | **None** in the long run (Q1c) |
| Classical CUSUM or GLR on squared returns | Parametric pre-change law | All parameters | Changes of the modelled kind | Asymptotic or simulated thresholds |
| Kernel or classifier betting tests (Shekhar and Ramdas 2024; Podkopaev and Ramdas 2023) | Two-sample or independence null | None | Changes in high-dimensional data | Exact; heavier computation, not needed in one dimension |

When to prefer which: the mean-based monitor is the right tool when the question is "has volatility moved?" and returns cluster (its null tolerates dependence through the band, and an alarm has a direction and a size). The rank detector is the right tool when the question is "has anything about the return distribution changed?", when the calibration window is short or the baseline is unreliable, and when the stream is close to i.i.d. after standardisation (for example after dividing by a good volatility forecast).

## 5. Implementation

`rankbet.py` (imports `avmonitor`; nothing in `src` imports it):

* `score`, `pvalues_against_reference` (vectorised over paths with a binary indexed tree, exact with ties, `O(n log n)` per path), `SequentialRanker` (streaming, `bisect` on a sorted list).
* `gauss_score`, `gauss_score_abs_dev`: the Gaussianised rank score and its band penalty.
* `RankBetConfig` (score kind, `a`, grids, bins, half-lives, prior, `eps`, band, family weights), `RankBetCore` (the three bettor families and the Shiryaev-Roberts mixture, vectorised over paths), `run_rank_detector` (first alarm per path and the share of `M_t` from each family at the alarm).
* `RankMonitor`: streaming calibrate, monitor, alarm, restart loop with the same protocol as `VolatilityMonitor` (NaN returns skipped), about 80 microseconds per step in pure Python against about 10 for `VolatilityMonitor` (`results/t_timing.csv`, one run on a shared machine; the cost grows slowly with the segment length because each insertion into the sorted reference list is linear). Each `RankEvent` records the family shares, the mean of the last 240 p-values and the learned density ratio of new to reference p-values per bin (`RankBetCore.density_ratio`), which says what changed.

`tsb_sim.py`: unit-variance simulators (Gaussian, variance-matched Student t, jumps at fixed total variance, a sign-skew law with the law of `|u|` held fixed, GARCH(1,1)), and `main_monitor_first_alarm`, a vectorised copy of the plug-in `VolatilityMonitor` (a test checks it alarm-for-alarm against the class).

`run_experiment.py`: parts `pilot`, `q1`, `q1b`, `q1c`, `q2`, `q3`, `q4`, `q5`, `q6`, `timing`, `figs`.

## 6. Experimental methodology

All synthetic runs use per-step units (unit pre-change variance), a 1440-step calibration window (the main monitor's `cal_len`), `A = 28,800` (5% per 1440-step day), fixed seeds and at most two worker processes. The main monitor is the plug-in `VolatilityMonitor` (`c = 3`, two-sided, each side at `2A`), run on the same paths.

* **P (design pilot, development seeds 1000 to 1002).** Seven configurations (memory of the histogram learners, prior strength, family weights, scale only, histogram only) on four changes and two nulls, 200 to 400 paths. The default configuration C0 was frozen from this table before any reported run (`results/p_pilot_design.csv`). It was not the fastest on every row; it was the most balanced.
* **Q1 [SYNTHETIC] Validity.** No change; Gaussian, Student t4 and Gaussian with jumps (0.5% per step, 8 sd); 1000 paths; false alarm shares within 1 day and within 7 days (bounds 5% and 35%) with Wilson 95% intervals; and the same with a 120-step (2 hour) calibration window.
* **Q1b [SYNTHETIC] Clustering.** GARCH(1,1), unit variance, persistence 0.95 and 0.99 (half-lives 14 and 69 steps); 500 one-day paths; bands 0 and 0.25 for both designs.
* **Q1c [SYNTHETIC] Fixed versus growing reference.** Gaussian null, reference windows of 240 and 1440 steps, 14 days, 500 paths.
* **Q2 [SYNTHETIC] Power.** Change at monitored step 500, 4000 steps after it, 400 paths. Variance up (sd x1.25, x1.5), variance down (sd x0.8, x0.67), tails at fixed variance (Gaussian to variance-matched t4 and t3, and t4 to Gaussian), jump intensity at fixed total variance (0.2% to 0.5% and to 1% per step, jump sd 6, diffusion variance adjusted), skew at fixed `|u|` law (moves beyond 1.5 sd become negative with probability 0.8). Median delay with a distribution-free 95% interval (binomial order statistics), detection share, false alarms before the change.
* **Q3 [SYNTHETIC] Truncation of a rank score.** Same paths with score `|u|` and `min(|u|, 3)`.
* **Q5 [SYNTHETIC] The proved mean-based alternative.** `VolatilityMonitor(baseline="cs")` (`delta = 0.01`) on the Q1 Gaussian paths and five Q2 scenarios.
* **Q6 [SYNTHETIC] Short references.** The same monitored data with references of 120, 240 and 1440 steps, four changes, 400 paths.
* **T. Cost.** Wall-clock time per step of `RankMonitor` and `VolatilityMonitor` on the same 20,000 Gaussian returns, median of 3 repeats.
* **Q4 [REAL].** Binance spot BTCUSDT and ETHUSDT one-minute bars, 1 July to 31 December 2025 (`load_binance(download=False)`, history from 1 April 2025 for the intraday profile), divided by `intraday_scale` (20-day lookback), exactly as in E6. Monitors: main monitor with bands 0 and 0.25; rank detector (absolute score) with `b = 0` and `b = 0.25`; the same on nonzero returns only (`b = 0.25`); signed score (`b = 0.25`); absolute score with a 120-minute reference (`b = 0.25`). All were fixed before the run. For every alarm: the share of `M_t` from each family; the learned density ratio of new to reference p-values per bin; the distance to the nearest alarm of the partner design (main monitor at band 0.25 for rank detectors, rank detector at `b = 0.25` for main monitors), matched within 60 minutes, against a chance level from randomly shifted alarm times; whether the partner was monitoring or calibrating at that moment; and, over the 240 minutes before the alarm against the segment's calibration window, the robust scale ratio, the truncated variance ratio (what the main monitor sees), a scale-free shape distance (Kolmogorov-Smirnov distance between `|u|` divided by its median in the two windows), the shares of zero returns and of moves beyond 3 robust sd, and the share of up moves; the same over the next 1440 minutes. Control: 2000 random times with the same reference delays.

## 7. Folder structure

```
alternatives/two_sample_betting/
├── README.md              this file
├── RESULTS.md             all numbers, labelled [SYNTHETIC] or [REAL], and the ranking
├── rankbet.py             conformal p-values, bettors, e-detector, streaming RankMonitor
├── tsb_sim.py             simulators and the vectorised plug-in main monitor
├── run_experiment.py      parts pilot, q1, q1b, q1c, q2, q3, q4, q5, q6, timing, figs
├── tests/
│   └── test_rankbet.py    9 tests, about 2 seconds, no data needed
└── results/
    ├── p_pilot_design.csv              design pilot (development seeds)
    ├── q1_validity.csv                 no change: 1-day and 7-day false alarms
    ├── q1b_garch.csv                   volatility clustering
    ├── q1c_fixed_reference.csv         growing versus fixed reference
    ├── q2_power.csv                    delays for ten changes
    ├── q3_truncation.csv               truncated versus raw rank score
    ├── q5_cs_baseline.csv              confidence-sequence main monitor
    ├── q6_short_reference.csv          reference length 120, 240, 1440
    ├── q4_real_summary.csv             alarm counts, overlap, blind time
    ├── q4_real_events.csv              every alarm with its window metrics
    ├── q4_real_profiles.csv            medians by group of alarms
    ├── q4_real_random_times.csv        control windows
    ├── q4_real_rank_only_alarms.csv    alarms only the rank detector raised while the main monitor watched
    ├── t_timing.csv                    microseconds per step, rank against main monitor
    ├── fig1_synthetic_delays.png
    ├── fig2_real_alarm_profiles.png
    └── fig3_real_density_ratio.png
```

Run from the repository root with the package installed: `python alternatives/two_sample_betting/run_experiment.py` (all parts except the pilot; about 10 minutes on two cores of a shared machine) and `python -m pytest -q alternatives/two_sample_betting/tests`. The real-data part needs the Binance files under `data/raw` (see `data/README.md`).

## References

* Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance*, 109, 105482.
* Fedorova, V., Gammerman, A., Nouretdinov, I. and Vovk, V. (2012). Plug-in martingales for testing exchangeability on-line. *Proceedings of the 29th International Conference on Machine Learning (ICML)*, 1639 to 1646. arXiv:1204.3251.
* Podkopaev, A., Bloebaum, P., Kasiviswanathan, S. P. and Ramdas, A. (2023). Sequential kernelized independence testing. *Proceedings of the 40th International Conference on Machine Learning*, PMLR 202, 27957 to 27993.
* Podkopaev, A. and Ramdas, A. (2023). Sequential predictive two-sample and independence testing. *Advances in Neural Information Processing Systems 36 (NeurIPS 2023)*. arXiv:2305.00143.
* Shekhar, S. and Ramdas, A. (2024). Nonparametric two-sample testing by betting. *IEEE Transactions on Information Theory*, 70(2), 1178 to 1203. arXiv:2112.09162.
* Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science*, 2(2), 229 to 260.
* Volkhonskiy, D., Burnaev, E., Nouretdinov, I., Gammerman, A. and Vovk, V. (2017). Inductive conformal martingales for change-point detection. *Proceedings of the Sixth Workshop on Conformal and Probabilistic Prediction and Applications*, PMLR 60, 132 to 153.
* Vovk, V., Gammerman, A. and Shafer, G. (2022). *Algorithmic Learning in a Random World*, second edition. Springer. doi:10.1007/978-3-031-06649-8.
* Vovk, V. (2021). Testing randomness online. *Statistical Science*, 36(4), 595 to 611.
* Vovk, V., Nouretdinov, I. and Gammerman, A. (2003). Testing exchangeability on-line. *Proceedings of the 20th International Conference on Machine Learning (ICML)*, 768 to 775.
* Vovk, V., Petej, I., Nouretdinov, I., Ahlberg, E., Carlsson, L. and Gammerman, A. (2021). Retrain or not retrain: conformal test martingales for change-point detection. *Proceedings of the Tenth Symposium on Conformal and Probabilistic Prediction and Applications*, PMLR 152, 191 to 210.
