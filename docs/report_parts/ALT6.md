## ALT6. Sequential fixed-accuracy estimation: how many minutes until the volatility level is known to plus or minus d?

Code, tests, README (formulation, proofs, comparison with the literature) and full results are in `alternatives/sequential_estimation/`. Every number below is copied from its `results/` folder. [SYNTHETIC] marks simulated data, [REAL] marks Binance spot BTCUSDT one-minute bars (1 July to 31 December 2025). Iteration notes: `docs/iterations/ALT6.md`.

**Question.** The main monitor asks whether volatility has changed. A simpler question comes first: after a restart, how long must one watch before the current level is known to within plus or minus `d` with confidence `1 - α`, and can the stopping rule keep its promise for heavy-tailed returns, jumps and clustering? The benchmark is the oracle sample size `n* = z² CV² / d²` (relative accuracy; `CV²` is the squared coefficient of variation of the statistic, unknown in practice).

**Method.** All rules work on squared returns after a 30-return pilot, either raw or truncated at `τ = 3 ×` the pilot's scaled MAD (the main monitor's statistic, Bapat 2026). Classical rules: Chow and Robbins (1965) with Nadas' (1967) relative version (stop when the estimated `n*` is reached), Stein's (1945) two-stage rule, Robbins' (1959) minimum risk rule. Anytime-valid rules: confidence sequences for the mean of the bounded truncated statistic (empirical Bernstein and betting, after Waudby-Smith and Ramdas 2024), stopped as soon as their relative width reaches `d`; their coverage at the stop holds for any law and any stopping time (Theorem 1 of the README). Proposition 1 (a sequential Bahadur and Savage argument) shows that no rule can do this for the raw variance, so the target is the truncated variance `θ_τ`, with a map back to `σ²` under a stated shape. Laws: normal, t5, t3, normal with 1% jumps per minute, GARCH(1,1) with persistence 0.95 and 0.99; 1000 to 4000 paths per cell, fixed seeds. Real data: one window per start hour (4392 windows) of returns divided by the past-only intraday profile, and a stationary bootstrap (Politis and Romano 1994) of the 30 quietest days.

**Classical rules on raw and truncated squares [SYNTHETIC].** Relative accuracy, `α = 0.05`, 2000 paths; "coverage / E[N]/n*":

| law | CR-raw, d = 0.2 | CR-raw, d = 0.05 | CR-trunc, d = 0.2 | CR-trunc, d = 0.05 | CS-bet, d = 0.2 | CS-bet, d = 0.05 |
|---|---|---|---|---|---|---|
| normal | 0.948 / 0.97 | 0.942 / 1.00 | 0.952 / 0.99 | 0.942 / 1.00 | 0.982 / 1.84 | 0.969 / 1.57 |
| t3 | **0.765** / N 2225 | **0.721** / N 126482 | 0.944 / 1.01 | 0.949 / 1.00 | 0.973 / 1.73 | 0.973 / 1.57 |
| jumps | **0.783** / 0.47 | 0.929 / 0.96 | 0.942 / 0.98 | 0.945 / 1.00 | 0.976 / 1.79 | 0.973 / 1.58 |

**The price of anytime validity [SYNTHETIC, closed form].** As `d → 0`, a confidence sequence stopped at width `2d` needs `2 log(2/α) / z²` times the oracle without the running intersection of its bounds (1.921 at `α = 0.05`, 2.215 at 0.10; coverage at the stop 0.993 and 0.986) and a Brownian first-passage constant with it (1.548 and 1.708). Simulated with intersection (normal returns, pilot 100): 1.702, 1.576, **1.545** at `d` = 0.1, 0.05, 0.025. The factor is the ratio of the squared Chernoff quantile to the squared normal quantile, so it is the price of an exponential tail bound in place of the central limit theorem, not of time-uniformity. Time-uniform bets cost 5.6 to 14.6 times the oracle.

**Volatility clustering [SYNTHETIC].** GARCH(1,1), `d = 0.1`, 1000 paths; "coverage of the stationary level / coverage of the average conditional truncated variance over the window used / E[N]" (long-run `n*` 5370 and 22979):

| procedure | persistence 0.95 | persistence 0.99 |
|---|---|---|
| CR-trunc | 0.514 / 0.956 / 783 | 0.317 / 0.948 / 769 |
| CR-trunc with batch means | 0.911 / 1.000 / 4992 | 0.917 / 1.000 / 18916 |
| CS-bet, running intersection | 0.455 / **0.713** / 774 | 0.279 / **0.682** / 666 |
| CS-bet, no intersection | 0.695 / 0.992 / 1726 | 0.427 / 0.986 / 1815 |

**Real data [REAL].** BTCUSDT, `α = 0.05`, median minutes until stop over 4392 hourly windows (range of the 24 hourly medians), including the 30-minute pilot; no coverage is claimed:

| target accuracy | CR-raw | CR-trunc | Stein-trunc | CS-bet, no intersection | CR-raw not stopped within the start hour |
|---|---|---|---|---|---|
| d = 0.4 (about plus or minus 20% on volatility) | 66 (52 to 82) | **50** (46 to 52) | 55 (49 to 58) | 140 (131 to 149) | **58%** (42% to 69%), CR-trunc 27% |
| d = 0.2 (plus or minus 20% on variance) | 411 (350 to 485) | **220** (200 to 230) | 220 (196 to 229) | 505 (447 to 551) | 100%; 6.0% not stopped within 24 hours |

Stationary bootstrap of quiet days, `d = 0.2`: with i.i.d. resampling CR-trunc covers 0.949, CS-bet 0.979, CR-raw 0.919; with 30-minute blocks every procedure covers 0.70 to 0.82.

**Reading.**
* Truncation is what makes the classical rules usable: CR-trunc covers 0.936 to 0.964 at 0.96 to 1.01 times the oracle in all 16 law and accuracy cells, for any pilot from 10 to 100 returns. On raw squares the rule fails for t3 at every `d` (0.72 to 0.78, with N growing like `d^{-3}`) and for jumps (0.78), and Stein's rule is worse (down to 0.36). In absolute rather than relative accuracy even normal returns under-cover at `d = 0.2` (0.90, a skewness bias of the stopped mean), so accuracy should be stated relatively.
* The anytime-valid betting sequence is valid everywhere and costs a known factor: about 1.55 times the oracle at small `d`, 1.8 at `d = 0.2`; on BTC, 2.3 to 2.8 times the minutes of truncated Chow and Robbins.
* Under clustering the stationary level is out of reach within hours (every i.i.d. rule stops 7 to 47 times too early; batch means needs 2,000 to 19,000 returns). What the classical rule does cover is the window's average conditional truncated variance (0.948 to 0.966, martingale central limit theorem). The running intersection of the confidence sequence is invalid for that moving target; without it the sequence covers it at 0.986 to 1.000. On real data the intersected sequence stops suspiciously early (1.32 against 2.20 times the in-window `n*` proxy at `d = 0.2`), the same signature, so the figures without intersection are the ones to quote.
* The hour of day matters little: the truncated rules' hourly medians stay within 10% of the overall median, because relative accuracy depends on the shape of the truncated squares, not on the level.

**Caveats.** The guarantee is for the truncated variance; the map to `σ²` is exact only under the assumed shape (2.5% coverage for t3 under a normal shape). The classical rules are asymptotic. Under clustering only the window average is covered, never the stationary level. The bets of the betting sequence are tuned for the rule without intersection (about 5% left in the intersected limit) and were not re-tuned; the pilot's returns are not used by the sequence. Deseasonalisation removes the intraday pattern but not quiet weekends, which cause most of the windows that never stop. One coin, one venue, six months.
