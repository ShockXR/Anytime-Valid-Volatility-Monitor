# Independent review, round 4: E17 (one monitor for 1, 5 and 15 minutes)

**Reviewer role.** An adversarial reader: a statistics journal referee and a code auditor. I did not write any of the files below. The goal is that every claim about E17 survives an expert reading.

**Scope.**

1. Theory: `docs/theory_multifreq.md` (Lemma 1′, section 2.1, Proposition 2, Lemma 2, Theorem 1, Proposition 3, Theorem 3, Corollary 4, the gap blocks under price noise).
2. Code: `src/avmonitor/multifreq.py`, `tests/test_multifreq.py`, `experiments/e17_multifrequency.py`.
3. Results: `docs/report_parts/E17.md`, `results/e17_*.csv`, `results/figures/e17_delays.png`, `docs/iterations/E17.md`.
4. How E17 is stated elsewhere: `README.md:32` (bullet "One monitor for several time scales"), `PITCH_BAPAT.md:43` (section 7), `docs/report_parts/_header.md:31` (row E17), `docs/limitations.md:11`.

**Read first.** `docs/math.md` sections 5 and 6 (Lemma 1, Theorems A and B, Corollary C, Theorem 6.4) and the earlier reviews in `docs/review/`. Their points are not repeated.

**State of the tree.** The review is against commit 46a21da on branch `claude/project-thread-ytr1da`. Other agents had uncommitted changes in the working tree while I worked (E13, E14 and E16 files); none of them touches a file in scope, so the line numbers below hold at 46a21da.

**What I ran.** Every number in `E17.md` was compared with its CSV, and every range, ratio and share was recomputed with pandas. The same was done for the E17 statements in the README, the pitch, the report index and the limitations. Scripts are in `/tmp/claude-0/-home-user/e30ce406-b0e8-5b23-b612-d49eb104d0ee/scratchpad/r4e17/`; nothing was written to `results/`.

* `prop2.py`: the constants of Proposition 2 in closed form and by `scipy.integrate.quad`; part 2 by Monte Carlo (4 million paths at λ = 0.25, 0.5 and 0.9) and part 3 (400,000 paths at B = 4 and 8); the section 2.1 counterexample; the trend variance factors; the log handicaps.
* `rerun.py`: reruns four cells of the experiment with the script's own functions and seeds, on 2 processes under `nice -n 10` (29 s): null, Gaussian, true baseline; null, γ = 2, plug-in baseline; change, γ = 2, sd ×1.25; change, Gaussian, trend. **Every count, detection share and delay equals `results/e17_null.csv` and `results/e17_change.csv` exactly (0 mismatches over 50 monitor rows).**
* `check_misc.py`: the vectorised and streaming monitors with gap blocks (not covered by the tests); the vectorised and streaming plug-in calibrations; one missing minute in a calibration window (M6); the one-minute truncated mean under the trend design (m5).
* `nanfix.py`: a NaN-safe version of `calibration_block_returns` (the fix proposed in M6).
* A count of missing one-minute returns in the five coins with E6's `prepare` (0 for every coin).
* `python -m pytest -q tests/test_multifreq.py`: **17 passed in 8.9 s.**

No file in the repository was changed except this one.

**Severity.**

* *Critical*: a wrong proof step or a wrong headline number.
* *Major*: an overclaim, a missing condition, a misleading comparison, or a code defect in shipped functionality.
* *Minor*: a wrong number outside the headlines, wording or clarity.

---

## 1. Summary verdict

**The proofs hold, and every number I traced agrees with its CSV. I found no critical problem.** I re-derived Lemma 1′, Lemma 2, Theorem 1 (both parts), Corollary 4, Proposition 2 (all three parts), Theorem 3 (both parts), Proposition 3 and the three sufficient conditions of section 2.1, and checked small cases by hand (section 5). The certified value `min(M, ℓ (M + 1))` is a sure lower bound on the wealth at the next block end, including the Shiryaev-Roberts restart term and a skipped block, on both sides. The code computes what the note says, and four rerun cells reproduce the CSVs exactly.

**The problems are in how the results are stated.**

* **M1.** Proposition 2 shows that the usual *proof* fails for the plain every-minute average, not that its *guarantee* is false, and Theorem 3 of the same note proves a bound for it (three times the budget). Yet `E17.md` says "the obvious guarantee is false", and the README, the pitch, the report index and the limitations say the plain average "has no proof".
* **M2.** The README, the pitch and the report index state the bound for "independent returns with any tails" without saying that it needs the true baselines. With the plug-in baselines that E17 itself uses, the one-minute detector and the equal-weight average exceed 5% per day significantly on independent returns (Q17b).
* **M3.** "In simulation the combination never beats the best single frequency" (README) and "never faster ... on a persistent change" (pitch) drop the qualifier "independent returns" that `E17.md` and the report index keep. Under price noise the averages beat the best single frequency in several cells.
* **M4.** The only positive result for a proved combination under noise ("it helps for increases") comes from γ = 2, a level added after the γ = 1 results were seen, and holds there only for the two smaller persistent increases. The disclosure is in the caveats, not next to the claim.
* **M5.** The reason given for not reading the coarse detectors' higher confirmation rates on real data ("E6 found the same trade-off when the band is widened") holds for down alarms only. An approximately rate-matched comparison already in the repository points the other way for up alarms.
* **M6.** Code: one missing minute in a calibration window silently removes every coarse-frequency calibration window after it (a NaN in a cumulative sum), or delays monitoring by two days. No reported number is affected (the coins have no missing minutes), but the module and the note say missing returns are handled.
* **M7.** The README, the pitch and the report index call the trend design "a slow trend". It is an abrupt rise in the lag-one autocorrelation from 0 to 0.3 with the one-minute variance unchanged, and the 95% to 97% is for Gaussian returns only (84% to 89% for Student t3).

**Most important fixes before the work is shown.**

1. Replace "has no proof" and "the guarantee is false" by "the usual proof does not cover it; Theorem 3 proves three times the budget" (M1).
2. Add "with the true baselines" wherever the bound is stated for independent returns (M2), and "with independent returns" wherever the combination is said never to beat the best single frequency (M3).
3. Put the post-hoc status of γ = 2 and the burst exception next to every "helps for increases" (M4).
4. Fix `calibration_block_returns` for missing values and add a test (M6).

---

## 2. Critical findings

None.

* **Proofs.** No proof step is wrong. The details are in section 5.
* **Headline numbers.** Every number in the E17 statements of the README, the pitch, the report index and the limitations agrees with its CSV. The problems are missing conditions (M2, M3, M7) and readings the theory or the data do not support (M1, M4, M5).

---

## 3. Major findings

### M1. Proposition 2 refutes a proof step, not the guarantee, and Theorem 3 proves a bound for the plain average

* **Where.**
  * `README.md:32`: "Averaging e-detectors run on 1, 5 and 15 minute returns is not automatically valid: ... and the plain average checked every minute has no proof."
  * `PITCH_BAPAT.md:43`: "Averaging detectors run on 1, 5 and 15 minute returns has no proof if checked every minute".
  * `docs/report_parts/_header.md:31`: "The plain every-minute average has no proof".
  * `docs/limitations.md:11`: "The plain average checked every minute behaved like the certified one in every simulation, but has no proof."
  * `docs/report_parts/E17.md:11`: "*The obvious guarantee is false.*"; `E17.md:22`: "**The obvious guarantee is false, and a repaired one is proved.**"
  * `docs/iterations/E17.md:7`: "the suggested bound is false"; `:11`: "So the suggested bound cannot be proved for the plain average checked every minute."
* **Evidence.**
  * **What Proposition 2 shows.** It refutes the inequality `E[M_t | F_{t-1}] ≤ M_{t-1} + 1{t ends a block}` and the optional stopping property built on it (`theory_multifreq.md:52` to `60`). I confirmed its constants (section 5). It says nothing about the probability bound itself, and the note says so: "Proposition 2 does not show that the plain every-minute average *raises more false alarms* than `Σ_k w_k ⌈m / k⌉ / A`; it shows that Lemma 1 cannot prove that it does not" (`theory_multifreq.md:64`).
  * **A proof exists.** Theorem 3, part 2 (`theory_multifreq.md:120`, `:124`) proves `P(false alarm within m) ≤ (Σ_k √(w_k N_k))² / A` for exactly this monitor. At the E17 level this is 15% per day (`results/e17_levels.csv`, row "average, balanced, plain every minute", `proved_bound_per_day` = 0.15), as `E17.md:24` itself says ("which has only a 15% proved bound at this level").
  * **The data.** The plain average is within 0.07 points of the certified one in all 12 null cells (`results/e17_null.csv`; for example 0.78% against 0.77% for Gaussian returns).
  * So "the guarantee is false" is not shown, and "has no proof" contradicts Theorem 3. What is true: the bound `Σ_k w_k m / k / A` is not proved for the plain average, and the bound that is proved is J times larger for balanced weights.
* **Fix.**
  * `README.md:32`: replace "Averaging e-detectors run on 1, 5 and 15 minute returns is not automatically valid: a coarse detector changes only at the end of its block, so between block ends it is not a fair bet, and the plain average checked every minute has no proof." with "Averaging e-detectors run on 1, 5 and 15 minute returns needs care: a coarse detector changes only at the end of its block, so between block ends it is not a fair bet, and the usual proof does not cover the plain average checked every minute (a union bound proves three times the budget for balanced weights)."
  * `PITCH_BAPAT.md:43`: replace "Averaging detectors run on 1, 5 and 15 minute returns has no proof if checked every minute, because a coarse detector is not a fair bet between its block ends." with "The usual proof does not cover an average of detectors run on 1, 5 and 15 minute returns checked every minute, because a coarse detector is not a fair bet between its block ends; a union bound proves only three times the budget."
  * `_header.md:31`: replace "The plain every-minute average has no proof (a coarse detector is not a fair bet between its block ends)" with "The usual proof does not cover the plain every-minute average (a coarse detector is not a fair bet between its block ends; a union bound proves three times the budget for balanced weights)".
  * `limitations.md:11`: replace "but has no proof." with "but its proved bound is three times larger (Theorem 3, balanced weights)."
  * `E17.md:11`: replace "*The obvious guarantee is false.*" with "*The obvious proof fails.*". `E17.md:22`: replace "**The obvious guarantee is false, and a repaired one is proved.** The plain every-minute average has no proof from the standard argument (Proposition 2)." with "**The obvious proof fails, and a repaired guarantee is proved.** The standard argument does not cover the plain every-minute average (Proposition 2); a union bound proves three times the budget for it (Theorem 3)."
  * `iterations/E17.md:7`: replace "the suggested bound is false" with "the suggested proof fails". `:11`: replace "So the suggested bound cannot be proved for the plain average checked every minute." with "So the usual argument does not prove the suggested bound for the plain average checked every minute; Proposition 2 does not show that the bound itself fails."

### M2. The bound is stated for independent returns without the true-baseline condition

* **Where.**
  * `README.md:32`: "gives a proved bound, P(false alarm within m minutes) ≤ Σ w_k ⌈m/k⌉ / A, whenever each frequency's own null holds (for example for independent returns with any tails)."
  * `docs/report_parts/_header.md:31`: "a certified average has P(false alarm within m) ≤ Σ w_k ⌈m/k⌉ / A when each frequency's null holds (for example independent returns with any tails)."
  * `PITCH_BAPAT.md:43`: "Counting it only at the wealth it is sure to keep gives a proved bound, Σ w_k ⌈m/k⌉ / A;" (no condition at all).
* **Evidence.**
  * Theorem 1 assumes the nulls "with the baselines used by the monitor" (`theory_multifreq.md:83`), and section 8 says that with the plug-in baselines the guarantee is only approximate. `E17.md:240` and `iterations/E17.md:23` say so. The three summaries do not, and they give "independent returns" as the example where the bound holds.
  * With the plug-in baseline on independent returns (`results/e17_null.csv`):
    * one minute alone: 5.47% (Gaussian, Wilson 5.04% to 5.93%), 5.54% (Student t3) and 5.68% (jumps), all significantly above 5%;
    * the equal-weight average: 5.60% for Student t3 and for jumps (Wilson 5.17% to 6.07%), both significantly above 5%.
  * The real-data runs that the README cites in the next sentence use the plug-in baseline.
* **Fix.**
  * `README.md:32`: replace "(for example for independent returns with any tails)" with "(for example for independent returns with any tails and the true baselines; with the plug-in baselines used on real data there is no proof, and in simulation the one-minute detector and the equal-weight average then exceed 5% per day)".
  * `_header.md:31`: replace "(for example independent returns with any tails)" with "(for example independent returns with any tails and the true baselines; no proof with the plug-in baselines used on real data)".
  * `PITCH_BAPAT.md:43`: replace "gives a proved bound, Σ w_k ⌈m/k⌉ / A;" with "gives a proved bound, Σ w_k ⌈m/k⌉ / A, when each frequency's null holds with known baselines;".

### M3. "Never beats the best single frequency" is true for independent returns only

* **Where.**
  * `README.md:32`: "In simulation the combination never beats the best single frequency, but stays within 5% (equal weights) to 18% (balanced weights) of it on persistent changes".
  * `PITCH_BAPAT.md:43`: "the combination is never faster than the best single frequency on a persistent change".
  * The report keeps the qualifier: `E17.md:23` ("in any i.i.d. design"), `E17.md:236` and `_header.md:31` ("in independent-return designs").
* **Evidence.** Under price noise the averages beat the best single frequency in several persistent-change cells (`results/e17_change.csv`; detection share and median delay among detected paths):
  * γ = 1, sd ×1.25: equal-weight average 95.9% (577 minutes) against the best single frequency, one minute, 92.1% (600.5);
  * γ = 1, sd ×1.5: equal-weight average 99.5% (188) against 99.5% (193);
  * γ = 2, sd ×1.5: balanced average 97.3% (563) against the best single frequency, 5 minutes on plain blocks, 93.8% (598);
  * γ = 2, sd ×1.25: balanced average 44.1% against 37.4%.

  None of these monitors has a proof under noise, but neither has any of the single frequencies they are compared with. `E17.md:210` also reports the proved 5 and 15 average beating both proved singles at γ = 2.
* **Fix.**
  * `README.md:32`: replace "In simulation the combination never beats the best single frequency" with "With independent simulated returns the combination never beats the best single frequency".
  * `PITCH_BAPAT.md:43`: replace "the combination is never faster than the best single frequency on a persistent change" with "with independent simulated returns the combination is never faster than the best single frequency on a persistent change".

### M4. The noise result rests on a post-hoc cell and is stated too generally

* **Where.**
  * `docs/report_parts/E17.md:25`: "detects increases better than either valid single frequency at γ = 2 ... It helps for increases, not for decreases."
  * `E17.md:210` (the γ = 2 reading) and `E17.md:236`, answer (3): "The proved 5 and 15 gap average helps for increases at γ = 2 and costs for decreases."
  * `docs/report_parts/_header.md:31`: "under price noise only the 5 and 15 minute gap average is valid, and it helps for increases, not decreases."
  * `docs/iterations/E17.md:63`.
* **Evidence.**
  * **Post hoc.** γ = 2 was added after the first run, because "at γ = 1 the one-minute detector ... was still the fastest single frequency ... so the design had no noisy case in which the coarse frequencies win" (`iterations/E17.md:55` to `59`). This is disclosed in `E17.md:244` and `iterations/E17.md:67`, but not in the summary bullet, the answer or the report index, which all present the γ = 2 result as the answer to "does it help under price noise?". The seeds are unaffected (the law was appended at the end of the list; my rerun of a γ = 2 cell reproduces it exactly), so the disclosure is the only issue with the procedure.
  * **At the planned γ = 1 it does not help for increases** (proved 5 and 15 gap average against the 5-minute gap detector, `results/e17_change.csv`): sd ×1.25: 78.7% (803) against 79.2% (763); ×1.5: 100% at 328 against 303 minutes; ×2: 100% at 138 against 128; burst ×3 for 60 minutes: 38.2% against 47.9%.
  * **At γ = 2 it helps for two of the four increase cells only.** ×1.25: 28.1% against 21.4% and 25.3%; ×1.5: 93.3% against 85.5% and 82.8%. But ×2: 100% at 258 minutes against 100% at 253 for the 5-minute gap detector; and the burst ×3 for 60 minutes, also an increase: 10.9% against 14.8%.
* **Fix.**
  * `E17.md:25`: replace "It helps for increases, not for decreases." with "It helps for persistent increases of 25% and 50% at γ = 2, a noise level added after the γ = 1 results were seen; it does not help at the planned γ = 1 (×1.25: 78.7% against 79.2% for the 5-minute gap detector), for the ×3 burst of 60 minutes (γ = 2: 10.9% against 14.8%; γ = 1: 38.2% against 47.9%) or for decreases."
  * `E17.md:236`: replace "The proved 5 and 15 gap average helps for increases at γ = 2 and costs for decreases." with "The proved 5 and 15 gap average helps for persistent increases of 25% and 50% at γ = 2 (a level added after seeing γ = 1), not at γ = 1, not for bursts and not for decreases."
  * `E17.md:210`: after "The proved 5 and 15 average beats both proved singles for ×1.25 ... and ×1.5 ...", add "(γ = 2 was added after the γ = 1 results; for the ×3 burst of 60 minutes it is worse, 10.9% against 14.8%)".
  * `_header.md:31`: replace "and it helps for increases, not decreases" with "and it helps only for persistent increases at the strongest noise level, which was added after seeing the first results".

### M5. Real data: the reason for not reading the coarse detectors' higher confirmation rates holds for down alarms only

* **Where.**
  * `docs/report_parts/E17.md:234`: "E6 found the same trade-off when the band is widened, so these shares cannot be compared with the one-minute monitor's at different alarm rates."
  * `E17.md:27` and `:245` ("not a comparison at equal alarm rates").
* **Evidence.**
  * The averages are compared with the one-minute monitor at almost equal alarm rates (0.65 to 0.68 against 0.68 per day), which is fair. The coarse detectors are not (0.43 and 0.26 per day), and the report says so. The question is whether the lower rate explains their higher confirmation rates.
  * E6 already ran the one-minute monitor at band 0.5 on the same coins and period (`results/e06d_pooled.csv`, row "Jul-Dec 2025, robust e-detector, 0.5"): 0.47 to 0.51 alarms per day per coin, close to the 5-minute detector's 0.43.
    * Up alarms: 79.4% confirmed (291 alarms, Wilson 74.4% to 83.6%), unchanged from 79.6% at band 0.25. The 5-minute detector: 90.0% (Wilson 85.4% to 93.3%, `results/e17_real_pooled.csv`). The intervals do not overlap.
    * Down alarms: 96.7% (153 validated, Wilson 92.6% to 98.6%), up from 88.3% at band 0.25. The 5-minute detector: 94.6% (90.1% to 97.1%).
  * So widening the band reproduces the down-alarm gain but not the up-alarm gain. "E6 found the same trade-off" is true for down alarms only. The comparison has limits (different bands, alarms correlated across coins, and the frozen real settings, `experiments/e17_multifrequency.py:52` to `53`, allow no other comparison), so it is a lead, not a result. But the report should not dismiss the up-alarm difference with a reason the data contradict.
* **Fix.**
  * `E17.md:234`: replace "E6 found the same trade-off when the band is widened, so these shares cannot be compared with the one-minute monitor's at different alarm rates." with "An unplanned check against E6 suggests that for up alarms this is not only a rate effect: the one-minute monitor at band 0.5 alarms at a similar rate (0.47 to 0.51 per day per coin) and its up alarms are confirmed 79% of the time (74% to 84%), against 90% (85% to 93%) for the 5-minute detector at 0.43 per day; for down alarms the wider band closes the gap (97% against 95%). The bands differ and alarms are correlated, so this calls for a planned rate-matched comparison; it is not a result."
  * Suggested check: the one-minute monitor at a band chosen in advance to give 0.43 alarms per day on a different half-year, then compared with the 5-minute detector on July to December 2025.

### M6. Code: one missing minute in a calibration window breaks the coarse-frequency calibration

* **Where.**
  * `src/avmonitor/multifreq.py:188` to `189`, `calibration_block_returns`: `c = np.cumsum(...)` followed by `c[..., L:] - c[..., :-L]`.
  * Used by `MultiFrequencyMonitor._block_cal` (`:383`, followed by `R[np.isfinite(R)]` at `:384`, which hides the problem) and by `plugin_calibration` (`:199`).
  * The claims: module docstring `multifreq.py:55` to `57` ("Missing returns: a block that contains a missing one-minute return is skipped ... the certified value covers this case"); `theory_multifreq.md:16`.
* **Evidence.**
  * A NaN in a cumulative sum makes every later partial sum NaN, so every overlapping window that *ends* after the missing minute is NaN, not only the windows that contain it. Example: `calibration_block_returns(r, 5)` with `r = 1, ..., 20` and `r[3]` missing returns 16 NaN values; only windows 0 to 3 should be NaN.
  * In the streaming monitor (`check_misc.py`, i.i.d. N(0, 1), `cal_len = 1440`, balanced weights):
    * missing minute at position 700: monitoring starts on time, but the 5 and 15 minute calibrations use only the windows before minute 700. The 15-minute scale drops from 3.80 to 3.27 and its baseline from 15.27 to 12.34 (19% lower), silently;
    * missing minute at position 0 or 10: the 15-minute frequency has no usable window, the window is extended to `cal_max = 2880`, still unusable, and the monitor restarts calibration. Monitoring starts at minute 4320 instead of 1440, two days late, with `cal_stats` reporting one extension and one restart.
  * **No reported number is affected.** The five coins have no missing one-minute return in July to December 2025 (counted with E6's `prepare`), and the synthetic streams have none. But E17 is offered as a monitor (README, pitch), and the one-minute data of other experiments do have gaps.
  * **The tests do not catch it.** `tests/test_multifreq.py:125` runs a 5-minute monitor on streams with 1% missing returns and `cal_len = 300`, so nearly every calibration there is truncated, but the test only compares the certified with the stale check. The other tests with missing values use `freqs=(1,)` or `known_baseline`.
  * **A related wording point.** With `gap=True`, a missing return in the left-out first minute of a block does not skip the block (`multifreq.py:489`). That is correct (the block does not use that return), but the docstring and `theory_multifreq.md:16` say "a block that contains a missing one-minute return is skipped".
* **Fix.**
  * Replace `multifreq.py:188` to `189` with (checked in `nanfix.py`: identical output without missing values; with one missing minute at 700, 1431 of 1436 five-minute windows survive instead of 696):

    ```python
    fin = np.isfinite(r)
    z = np.zeros(r.shape[:-1] + (1,))
    c = np.cumsum(np.concatenate([z, np.where(fin, r, 0.0)], axis=-1), axis=-1)
    nbad = np.cumsum(np.concatenate([z, (~fin).astype(float)], axis=-1), axis=-1)
    out = c[..., L:] - c[..., :-L]
    return np.where(nbad[..., L:] - nbad[..., :-L] > 0, np.nan, out)
    ```
  * Add a test: with one missing minute at position 0 of a 1440-minute calibration window, the monitor calibrates at minute 1440, and the 15-minute scale is within 5% of the value without the missing minute.
  * `multifreq.py:55` and `theory_multifreq.md:16`: replace "a block that contains a missing one-minute return is skipped" with "a block in which one of the returns it sums is missing is skipped".

### M7. "A slow trend" misdescribes the trend design, and the 95% to 97% is for Gaussian returns only

* **Where.**
  * `README.md:32`: "catches 95% to 97% of a slow trend that the one-minute detector almost never sees (3.8%)".
  * `PITCH_BAPAT.md:43`: "catches 95% to 97% of a slow trend the one-minute detector almost never sees".
  * `docs/report_parts/_header.md:31`: "catches 95% to 97% of a slow trend the one-minute detector almost never sees (3.8%)".
* **Evidence.**
  * The design (`experiments/e17_multifrequency.py:147` to `158`, `:367`; `E17.md:17`) is a step change at minute 247 in the lag-one autocorrelation of the efficient returns, from 0 to 0.3, with the one-minute variance unchanged. Nothing in it is slow, and it is not a trend in prices or in volatility. A reader of the README will think of a gradual rise in volatility, which the one-minute detector would see.
  * The 95% to 97% (95.5% to 97.3%) and the 3.8% are for Gaussian returns. With Student t3 returns the averages catch 83.9% to 89.0% and the one-minute detector 17.5% (`results/e17_change.csv`; `E17.md:169`).
* **Fix.** In all three places, replace the phrase "catches 95% to 97% of a slow trend (that) the one-minute detector almost never sees" (the README and the report index add "(3.8%)") with "catches 95% to 97% of paths (Gaussian returns; 84% to 89% for Student t3) after a jump in the lag-one autocorrelation of returns from 0 to 0.3, which leaves the one-minute variance unchanged and which the one-minute detector almost never sees (3.8%)".

---

## 4. Minor findings

### m1. Proposition 2: a rounded constant and a loose test

* **Where.** `docs/theory_multifreq.md:60` ("`E[(Z² - 1) 1{|Z| > 1}] = 0.4840`"); `tests/test_multifreq.py:252` ("= 0.48396 lambda") and `:263` (`abs=0.03`).
* **Evidence.** `E[(Z² - 1) 1{|Z| > 1}] = 2 (φ(1) + 1 - Φ(1)) - P(|Z| > 1) = 0.80125 - 0.31731 = 0.48394` (closed form and `quad` agree). "0.484 λ" is right to three digits; "0.4840" and "0.48396" are not. The test tolerance 0.03 at λ = 0.5 is about 13 Monte Carlo standard errors (my run with 4 million paths: 0.24187 ± 0.00072 against the exact 0.24197; with the test's 400,000 paths the standard error is about 0.0023).
* **Fix.** `theory_multifreq.md:60`: "= 0.4839". `tests/test_multifreq.py:252`: "= 0.48394 lambda". `:263`: `abs=0.01`.

### m2. Section 2.1 (c): the decrease side does not fail for the stated reason

* **Where.** `docs/theory_multifreq.md:44`: "This does not extend to truncated statistics (point 1) or to the decrease side (`docs/math.md` 5.3)."
* **Evidence.** Under martingale differences `E[R² | F_s] = Σ_t E[r_t² | F_s]` exactly, so `E[r_t² | F_{t-1}] ≥ μ` for every t also gives `E[R² | F_s] ≥ L_k μ`. The implication extends; what fails is that no e-detector can test an untruncated decrease null (`math.md` 5.3).
* **Fix.** Replace the sentence with "This does not extend to truncated statistics (point 1). The same identity gives the untruncated decrease null at k minutes, but that null cannot be tested by any e-detector (`docs/math.md` 5.3)."

### m3. Price noise: two conditions are left out

* **Where.** `docs/theory_multifreq.md:46` ("For Gaussian noise and conditioning on the previous block alone, `E[R_b² | R_{b-1}] = ...`"); `docs/report_parts/E17.md:14` ("under independent noise, gap blocks are independent and Theorem 1 applies exactly") and `:241` ("only for independent noise").
* **Evidence.** The conditional-mean formula uses joint normality of `(R_b, R_{b-1})`, which needs Gaussian efficient returns as well as Gaussian noise. The exact validity of gap blocks needs the efficient price to have independent increments (the note says so: "constant σ, i.i.d. jumps"); `E17.md` drops this. Also, at `E17.md:14` "Theorem 1 applies exactly" is next to "gap versions", which include the three-frequency averages that still contain the one-minute detector.
* **Fix.** `theory_multifreq.md:46`: "For Gaussian noise and Gaussian efficient returns, conditioning on the previous block alone". `E17.md:14`: replace "under independent noise, gap blocks are independent and Theorem 1 applies exactly" with "under independent noise and an efficient price with independent increments, gap blocks are independent of the past, so Theorem 1 applies exactly to monitors built only from gap blocks". `E17.md:241`: "and only for independent noise on an efficient price with independent increments".

### m4. "Price ticks" are not covered by independent returns

* **Where.** `docs/report_parts/E17.md:10`: "If one-minute returns are independent with one fixed distribution (this allows fat tails, jumps and price ticks)".
* **Evidence.** Rounding prices to a tick grid makes returns dependent (the rounding error acts like price noise), so it is not covered by case (a). The note's own wording is correct ("as long as they are i.i.d.", `theory_multifreq.md:42`).
* **Fix.** Replace "(this allows fat tails, jumps and price ticks)" with "(this allows fat tails and jumps; prices rounded to a tick give dependent returns and are not covered)".

### m5. The trend design moves the one-minute truncated mean under Student t3

* **Where.** `docs/report_parts/E17.md:17` ("while the one-minute variance stays the same, so only block variances move") and `:169`.
* **Evidence.** The AR(1) filter changes the marginal law of a Student t3 return while keeping its variance. With the true τ² and baseline of `results/e17_true_baselines.csv`, the one-minute truncated mean rises to 1.045 times its baseline (Gaussian: 0.999, unchanged; `check_misc.py`, 2 million draws). This explains why the one-minute detector catches 17.5% of trends under t3 against 3.8% under Gaussian returns.
* **Fix.** `E17.md:17`: replace "so only block variances move" with "so only block variances move (for Gaussian returns; for Student t3 the one-minute truncated mean also rises by about 4%)". `E17.md:169`: add "because the change also raises the one-minute truncated mean by about 4% for this law".

### m6. "Certification costs nothing visible" is slightly too strong

* **Where.** `docs/report_parts/E17.md:24` and `:94`.
* **Evidence.** The plain average is never below the certified one, so it alarms no later. Where coarse evidence for a decrease matters (decrease bets have ℓ as low as 0.05), the cost is visible but small: γ = 2, sd ×0.5: 67.8% against 68.8%; ×2/3: 34.1% against 35.0% (`results/e17_change.csv`, balanced, certified against plain). Over all 44 change cells the detection shares differ by at most 0.95 points and the median delays by at most 5 minutes.
* **Fix.** `E17.md:24`: replace "**Certification costs nothing visible.**" with "**Certification costs little.**", and add after "the same delays as the certified one": "(at most about one point of detection, for decreases under strong noise)".

### m7. Theorem 3 is called "weaker" for the weighted maximum, whose bound is the same

* **Where.** `docs/report_parts/E17.md:12` ("Theorem 3 proves weaker bounds without certification: `Σ_k w_k ⌊m / k⌋ / A` for the weighted maximum"); section title `docs/theory_multifreq.md:115`.
* **Evidence.** `⌊m / k⌋ ≤ ⌈m / k⌉`, so the weighted maximum's bound is equal to or smaller than Theorem 1's. What it lacks is a restart result (`theory_multifreq.md:136`). Only the plain average's bound is weaker.
* **Fix.** `E17.md:12`: replace "Theorem 3 proves weaker bounds without certification:" with "Theorem 3 proves bounds for two monitors without certification:", and add "(the weighted maximum has the same budget but no restart result)".

### m8. The figure's title claims a proof that the noise panels do not have

* **Where.** `experiments/e17_multifrequency.py:570` (suptitle "... true baseline, 5% per day proved"), `:565` (panel titles "(1 min not valid)"); `results/figures/e17_delays.png`.
* **Evidence.** Under noise the figure plots plain coarse blocks and averages that contain the one-minute detector; none of them has a proof, not only the one-minute detector. In the γ = 2 panel the medians are among as few as 10% to 16% of paths (sd ×0.8) and close to the 1440-minute horizon, so the one-minute line looks fastest for decreases although it detects few of them.
* **Fix.** Panel titles: "(no plotted monitor has a proof)"; suptitle: "5% per day proved for the i.i.d. panels". Add the detection share as a second row of panels, or mark points with detection below 50%.

### m9. Two paths used in E17 are not covered by the tests

* **Where.** `tests/test_multifreq.py:190` to `205`; `docs/theory_multifreq.md:146` ("the vectorised code used in E17 equals the streaming monitor").
* **Evidence.** The test of vectorised against streaming code uses plain blocks and known baselines only. E17 also uses gap blocks and the vectorised plug-in calibration. I checked both (`check_misc.py`): with gap blocks the first alarms agree in 60 of 60 streams over the three combine and check modes, and the plug-in τ² and baselines agree for plain and gap blocks. So nothing is wrong, but the claim at `:146` is broader than the tests.
* **Fix.** Add `gap=True` to the parametrisation of `test_vectorised_matches_streaming`, and a test that `plugin_calibration` equals the streaming calibration.

### m10. The share of the statistic for the weighted maximum

* **Where.** `docs/report_parts/E17.md:217` (column "1-min share of the statistic"), `:225`, `:232`.
* **Evidence.** For the weighted maximum the share is computed in the weighted *sum* of the frequencies (`multifreq.py:520` to `523`), not in the statistic that alarms (the maximum).
* **Fix.** Add a footnote to the table: "for the weighted maximum, share in the weighted sum of the frequencies, not in the maximum".

### m11. The price of not knowing the frequency: the scope of "averaging cannot remove"

* **Where.** `docs/theory_multifreq.md:5`: "show that averaging cannot remove the price of not knowing the right frequency: in the worst case it is a factor J in the threshold".
* **Evidence.** Proposition 3 is about the certified average at the level set by Theorem 1's bound. It does not cover other ways of averaging or a sharper bound.
* **Fix.** Replace with "show that, at the level set by Theorem 1, averaging cannot remove the price of not knowing the right frequency: for any weights the worst case is a factor J in the threshold".

### m12. Small wording points

* `docs/theory_multifreq.md:89`: "the probability that it ever reaches A before `T + m`" should read "at some minute up to `T + m`" (the alarm at `T + m` itself is covered).
* `experiments/e17_multifrequency.py:35` to `36` says the balanced average of 5 and 15 minutes is run "under price noise also", but `:113` runs it for every law, as `E17.md:14` says. Move it to the list of monitors run for every law.
* `docs/report_parts/E17.md:18`: "the share of alarms beaten by (up) or beating (down) an EWMA volatility forecast" is hard to parse; `E17.md:215` has the clear version ("the share of up (down) alarms after which the realised volatility was above (below) the EWMA forecast"). Use that.

---

## 5. What I checked and found correct

**Proofs (`docs/theory_multifreq.md`).**

* **Lemma 1′**, including the missingness step: `E[f_b 1{complete} | F_s] ≤ P(complete | F_s)`, giving `E[M_b | F_s] ≤ M_{b-1} + 1`.
* **Lemma 2, part 1.** Inside a block with start s, the end value is `f (M_s + 1)` with `f ≥ ℓ_j` surely, or `M_s` if the block is skipped, so `M_{e(t)} ≥ min(ℓ_j (M_s + 1), M_s)` surely. This includes the Shiryaev-Roberts restart term (the `+ 1`). The certified value is `F_s`-measurable, so `Q_t ≥ C_t`.
* **The two sides.** Increase: `f = 1 + λ (X / m↑ - 1) ≥ 1 - λ` because `X ≥ 0`. Decrease: `f = 1 + λ (1 - X / m↓) ≥ 1 - λ (B - 1) = 1 - φ` because `X ≤ τ²`. The code uses `l_up = 1 - grid` and `l_dn = 1 - fractions` (`multifreq.py:228` to `229`, `:340` to `341`), and switches the decrease side off when `B ≤ 1`.
* **Lemma 2, part 2.** Inside a block (including its last minute) Q is a martingale by the tower property; at the first minute of a block `Q_{t-1} = M_s` and Lemma 1′ gives the `+ 1`. For k = 1 this is Lemma 1 of `math.md`.
* **Theorem 1, part 1.** `c(t) = Σ_k w_k ⌈(t - T) / k⌉` counts blocks whose first minute is at most t, so `c(t) - c(t - 1)` is the weighted indicator of Lemma 2; the factor 1/2 of the two sides gives exactly 1 per side pair. `S_t + c(T + m) - c(t)` is a nonnegative supermartingale started at `Σ_k w_k ⌈m / k⌉`. An alarm gives `(C̄↑ + C̄↓) / 2 ≥ A`, hence `S ≥ A`. `⌈m / k⌉ ≤ (m + k - 1) / k` gives `m ρ + κ`.
* **Small cases by hand.** m = 1, any weights: bound 1/A, and indeed `E[S_{T+1}] ≤ 1`. A lone 15-minute detector, m = 14: bound 1/A, true probability 0. A lone 15-minute detector, m = 15: 1/A, as Theorem B with one block. Frequencies {1, 2} with weights 1/2: m = 2 gives 1.5/A. Balanced {1, 5, 15}, m = 1440: `3 × 1440 / 21 / A`, so A = 4114.29 for 5%.
* **Theorem 1, part 2.** Optional stopping at `τ ∧ n`, Fatou with `S_τ ≥ A`, and monotonicity of c give `A ≤ E[c(τ)] ≤ ρ E[τ - T] + κ`. For k = 1 alone it is Theorem A. For a lone k-minute detector it gives `k A - (k - 1)`, slightly below the block monitor's `k A`: valid, a little loose.
* **Corollary 4.** At an alarm the current S is banked unchanged (it is `F_t`-measurable), an abandoned partial block was counted when it started, calibration minutes change nothing, and at most `N_m + 1` segments have monitored minutes among the first m. With k = 1 alone it reduces to Corollary C (`E[N_m] ≤ m / A`). The explanation at `:136` of why no restart result is claimed for the plain average or the maximum is right.
* **Proposition 2.** Part 1 directly. Part 2: `E[M_τ] = 1 + p + λ E[(Z² - 1) 1{|Z| > 1}]` with p = 0.31731 and the constant 0.48394 (closed form and quadrature); Monte Carlo with 4 million paths: 0.12105, 0.24187 and 0.43496 at λ = 0.25, 0.5 and 0.9, against 0.12099, 0.24197 and 0.43555. Part 3: `E[M_τ] - B = 0.24197 λ B - 0.68269`, positive for `B > 2.821 / λ` (from the 6th block at λ = 0.5); Monte Carlo at λ = 0.5: B = 8 gives 0.285 against 0.2852, B = 4 gives -0.186 ± 0.007 against -0.199.
* **Theorem 3.** Part 1: the union bound over k with thresholds `A / w_k` gives `Σ_k w_k N_k / A`. Part 2: the pigeonhole step and Cauchy-Schwarz in the right direction (the minimum of `Σ_k N_k / a_k` subject to `Σ_k w_k a_k = A` is `(Σ_k √(w_k N_k))² / A`). For balanced weights and k dividing m the ratio to Theorem 1 is exactly J (3), so 15% at the E17 level, and the single-frequency handicap at the level that proves 5% would be J² = 9.
* **Proposition 3.** `Σ_j v_j h_j = J` with v a probability vector, so the maximum is at least J, with equality iff all `w_j ⌈m / j⌉` are equal, which is `w_j ∝ j` when every k divides m. Handicaps at m = 1440: equal 1.267, 6.333, 19.0; balanced 3, 3, 3; sqrt 1.705, 3.813, 6.605; 5 and 15 balanced 2, 2 (all as in `results/e17_levels.csv`). Logs: log 3 = 1.099, log 57,600 = 10.961 (10.0%), log 3,840 = 8.253 (13.3%), log 19 = 2.944.
* **Section 2.1.** (a) A fixed conditional law makes the returns i.i.d. and independent of `F_T`, so every frequency's null holds with equality at the true baseline (τ_k is `F_T`-measurable). (b) Monotonicity of `min(σ² y², τ²)` in σ, and nesting of the block grids when the scale is constant over the coarsest blocks; the Gaussian extension for scales known at the block start. (c) The cross terms vanish under martingale differences (see m2 for the decrease side). The counterexample: `r_t² = e_t²` is independent of `F_{t-1}`, so the one-minute null holds with equality; `E[R²] = k + k (k - 1) (2 / π)` = 17.73 for k = 5 (3.55 times 5).
* **Price noise.** `Cov(R_b, R_{b-1}) = -ω²`, `Var(R_b) = k + 2 γ²`, `ρ_k = -γ² / (k + 2 γ²)`; `ρ_k²` = 0.111, 0.020, 0.0035 at γ = 1 and 0.095, 0.030 at γ = 2 (5 and 15 minutes); one-minute correlations -0.167, -0.333, -0.444; ALT1's 0.111. **Gap blocks are exactly valid under the stated model**: `R_b = (X_{s+k} - X_{s+1}) + u_{s+k} - u_{s+1}` is a function of variables independent of `(X_0, ..., X_s, u_0, ..., u_s)` when X has independent increments and the noise is i.i.d. and independent of X; no Gaussianity is needed. The calibration windows of length `k - 1` have the same law as the monitored gap blocks.
* **Trend factors.** For AR(1) with ρ = 0.3 and unit variance, `Var(Σ) / k` = 1.613 (k = 5) and 1.776 (k = 15), as stated.

**Code against the math (`src/avmonitor/multifreq.py`).**

* Blocks: minute i of a segment ends a block when `i % k == k - 1` (`:494`), the theory's `T + b k`; the gap leaves out the first return of each block in the streaming (`:489`) and vectorised (`:173` to `174`) code, which agree (60 of 60 streams).
* The certified value `mean_j min(S_j, ℓ_j (S_j + 1))` is used strictly between block ends and the full value at block ends (`:456` to `461`, `:494` to `503`; vectorised `:243` to `244`, `:248` to `266`), 0 before the first block ends.
* Calibration: τ_k = 3 × scaled median |R| of all windows of `L_k` returns in the calibration window, and `m_k` the mean of the clipped squares; vectorised and streaming agree. Restart after an alarm: `_reset` clears the statistics and the buffer, the next `cal_len` minutes recalibrate every frequency, and `pos = 0` restarts every block grid at the new monitoring start, as Corollary 4 assumes.
* `proved_bound`, `level_for_budget` and `handicap` implement Theorem 1, Theorem 3 and section 5. `results/e17_levels.csv` agrees with my hand computation (A = 28,800, 5,760, 1,920, 12,160, 4,114.29, 6,908.92 and 2,880).
* The experiment: the noise model (`u ~ N(0, γ²)` per price, so γ = ω / (σ √δ)), the block-law draws for the true baselines (two noise ends per block), the change at index 247 (inside a 5 and a 15 minute block), the burst window (burst end plus 15 minutes) and the delay count (`idx - 247 + 1`) match the report.

**Numbers.**

* Every cell of the Q17a to Q17f tables of `E17.md` equals its CSV value after rounding, including every bold mark (lower Wilson end above 5%), every bracketed detection share and every median delay.
* Every range and ratio in the readings and the summary: 0.32% to 2.62%; balanced 0.77% to 1.24%; within 0.07 points; 3.66% to 4.54%; 2.7 to 11.9 times slower (2.67 to 11.88); delay ratios 1.00 to 1.04 and 1.00 to 1.05 (equal), 1.04 to 1.08 and 1.06 to 1.10 (sqrt), 1.13 to 1.16 and 1.11 to 1.18 (balanced); 1.17 and 1.32 times the 5-minute trend delay; "at most 1.7%" (1.65%) of paths with an alarm before the change; Monte Carlo standard error at most 0.40 points; 0.65 to 0.68 alarms per day; 98.4% to 99.7% one-minute share; 72% to 94% with a one-minute alarm within 60 minutes; per-coin rates 0.66 to 0.70, 0.41 to 0.45, 0.22 to 0.33, 0.66 to 0.70 and 0.64 to 0.66 (184 days per coin); median lags 6 to 10 and 1 to 2 minutes; "equal or 1 to 4 points lower" (true for the rounded table values; the largest exact gap is 4.6 points, down alarms, 83.7% against 88.3%).
* E6 counts 127, 122, 129, 123 and 126 (`results/e06d_monitor_summary.csv`, band 0.25); E13 lag-one autocorrelations -0.016 to +0.027 and +0.132 (`results/e13d_real_crypto_random.csv`).
* README, pitch, report index and limitations: 5% and 18%, 0% to 5% and 11% to 18%, 95% to 97%, 3.8% all agree with the CSVs (the problems are the conditions, M2, M3 and M7). The report index row and `E17.md` are copied verbatim into `docs/experiment_report.md`.

**Reproducibility and procedure.**

* `pytest`: 17 passed. Four rerun cells reproduce the CSVs exactly.
* The real-data settings in the script's docstring were committed in 34c6155, before the real results (5da9ac2); between the two commits the script changed only in its runtime note and the figure code, so the real-data code was frozen as claimed.
* The γ = 2 addition is disclosed in `E17.md:244` and `iterations/E17.md:55` to `67`, and appending the law kept every earlier seed. The issue is where the disclosure sits (M4), not the procedure.
* The comparison of the averages with the one-minute monitor on real data is at almost equal alarm rates, and the report says the coarse singles are not (M5 is about the reason it gives).
