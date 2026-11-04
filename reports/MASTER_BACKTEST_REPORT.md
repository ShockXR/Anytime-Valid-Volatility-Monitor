# Master validation report: anytime-valid volatility monitor (avmonitor)

Date: 3 November 2026. Code state: PR #2 head (d2e41b6) plus the validation changes on this branch. All checks were run by AI agents (Sonnet), not by a human reviewer, and every claim below says where its evidence is. The monitor is not a trading system, so it is judged on false-alarm rate against a proved bound, detection delay and power, honesty about failure cases, and reproducibility.

## 1. Summary in plain language

1. **The main guarantee holds on fresh random numbers.** We ran 3,000 simulated paths for each of 8 return laws using seeds that no earlier experiment used. For independent returns, including infinite-variance Pareto tails, Student t with 3 degrees of freedom and jump contamination, no proved setting exceeded its bound: at most 1.7% false alarms within a day (bound 5%) and 14.7% over five days (bound 25%). [VERIFIED, SYNTHETIC]
2. **The guarantee has a price: speed.** The README setting needs a median of 112 steps to flag a doubling of variance, 2.6 times slower than a Gaussian likelihood-ratio rule (43 steps); a 1.5x increase takes 506 steps, more than a trading day; a fall of variance to two thirds is missed in 23% of runs within five days. [VERIFIED, SYNTHETIC]
3. **The guarantee breaks when the assumptions break, as the documents already say.** With GARCH volatility clustering or an uncorrected intraday pattern, every variant exceeds its bound, including the README setting (7% and 47% within a day against 5%). The plug-in baseline without a safety band is only approximately valid: 26.7% to 28.0% over five days against a nominal 25%. [VERIFIED]
4. **The code was hardened.** Six places where bad input was silently accepted and turned the detector into a dead switch now raise clear errors (threshold NaN or not positive, calibration length below 2, truncation constant NaN or infinite, nonpositive baseline, negative statistics, a crash in the repeated-window z-test). Each fix has a test that failed before the fix. No printed number changed. [VERIFIED]
5. **Robustness to messy data is still mixed.** Of 93 injected failures, 22 gave clean errors, 34 were handled safely, 25 were silently accepted and 12 crashed unhelpfully (counts measured before the six fixes above; section 6 lists what is still open). [VERIFIED]

**Presentation readiness: ready for a professor meeting WITH the listed limitations, not ready for release.** Gate checks are in section 7.

## 2. What was examined

| Item | Detail |
|---|---|
| Prototype | `avmonitor`: robust truncated squared returns, Shiryaev-Roberts and CUSUM betting e-detectors, confidence-sequence and plug-in baselines, rank monitor, guard, 25 experiments (E1 to E25), theory TH4 to TH6 |
| Data | REAL: Binance 1-minute bars (BTC, ETH, XRP, SOL, DOGE), Bitstamp BTC, Nifty 50 1-minute bars, India VIX (raw files not in the repository, hashes in `data/MANIFEST.csv`). SYNTHETIC: `avmonitor.simulate` |
| Environment | Python 3.13, 4 cores, exact versions in `requirements.txt` |
| Baselines compared | package versions of Gaussian likelihood-ratio e-detector, Page CUSUM, repeated z-test, plug-in baseline, confidence sequence |
| Not done | Full-tier rerun of E18 to E25 on a clean clone; real-data integrity of 137 of 140 data files (only BTC 2024-01 to 2024-03 was fetched and checked) |

## 3. Baseline against revised results

| Measure | Before this validation | After | Explained by |
|---|---|---|---|
| Tests | 581 passed (104 s) | 637 passed (98 s) | 56 new audit tests; 12 of them failed before the fixes |
| Quick reproduction tier | byte-identical to `results/MANIFEST.csv` | byte-identical; manifest rebuilt so 465 of 465 files are identical (three derived trace files had drifted by rounding only) | trace files rewritten by the number tracer |
| Numbers traced | 2,552 in the paper, 282 in the README, 199 in the pitch | same, 0 untraced | no number changed |
| Documents | test count "about 240" and "570"; no statement that the "independent" reviews were AI agents; the E22 rerun presented without saying it was designed after the E12 result | counts 637; AI authorship stated in the README, pitch, strategy file, limitations and all 15 review files; E22 caveat added | see `METHODOLOGY_VALIDATION.md` |
| Figures | E1 and E4 figures had typed constants in labels | regenerated from the code constants; pixel difference only in the title rows | `VISUAL_QA.md` |
| README | `avmonitor estimators` failed as written | now `estimators --source binance --symbol BTCUSDT --day 2025-10-10 --download` (run and checked) | |

## 4. Evidence by dimension

| Dimension | Finding | State |
|---|---|---|
| Functional correctness | Statistic nonnegative, mean about t under the null, alarm level exactly the threshold, no use of future steps; 16 equation rows checked against code, 20 of 20 numeric checks pass | VERIFIED |
| Guarantee | See section 1, points 1 to 3 | VERIFIED, SYNTHETIC |
| Power and delay | See point 2; against the package's own baselines at matched settings | VERIFIED, SYNTHETIC |
| Sensitivity | 42 one-at-a-time configurations: gradual change except cliffs when the truncation constant c is 1 or below and when the safety band is as wide as the effect | VERIFIED |
| Real-data value | Alarms beat random times; little gain against a strong HAR forecast (E20); the FOMC contrast failed its confirmation test on 2026 data (E23: 2 of 16 against 5.9%, p = 0.25), so it is exploratory | as reported in the repo, not rerun here |
| Reproducibility | Quick tier identical; full tier E1 to E17 identical on clean clones (authors); E18 to E25 not rerun on a clean clone | PARTIALLY VERIFIED |
| Interpretability and design | Core modules read line by line, no wrong formula found | VERIFIED (no human check of the proofs) |

## 5. Data checks

BTC 1-minute bars for 2024-01 to 2024-03 match the manifest hashes, with no duplicates, gaps or timestamp disorder. That is 3 of 140 files. Loader weaknesses found: the Bitstamp loader keeps an arbitrary one of duplicated minutes (unstable sort, 286 duplicated minutes in the test case), and the Bitstamp and Nifty loaders do not mask nonpositive prices (the Binance loader does). [PARTIALLY VERIFIED]

## 6. Open items, most serious first

1. **Loaders** (A4-B-001, A4-B-002): make the duplicate rule deterministic and mask nonpositive prices. Not changed here because it could move published numbers; run the reproduction manifest check after any change.
2. **Monitor input checks** (A4-B-003 to A4-B-005): calibration lengths 2 to 10 are still accepted (the null alarm rate is inflated), non-numeric input fails late, and all-NaN or constant streams leave the monitor blind with no warning.
3. **Benchmark detectors die after one NaN** (A1-B-007): the Gaussian likelihood-ratio and Page CUSUM comparators; the main detector skips missing values.
4. **CLI errors**: a mistyped date, unknown symbol or empty file prints a raw traceback.
5. **E18 to E25 not rerun on a clean clone**; the tracer accepts coincidental matches (A1-B-016).
6. Forward tests E24 (after 1 June 2027) and E25 (from 26 November 2026) are pre-registered and must be run exactly as written.

## 7. Presentation-readiness gates

| Gate | State |
|---|---|
| Primary functionality demonstrated | VERIFIED (CLI and replay page run, `DEMO_WALKTHROUGH.md`) |
| Critical calculations checked | VERIFIED (20 of 20), proofs read but not checked by a human |
| Important tests executed | VERIFIED (637 pass) |
| Major risks documented | VERIFIED (`RISK_REGISTER.md`) |
| Charts correspond to results | VERIFIED for the replay page and paper Figure 1 (alarm lists equal the result files); other figures sampled |
| Documentation matches implementation | VERIFIED for traced numbers, AI authorship now disclosed |
| Demo path validated | VERIFIED, replay page works offline with 0 console errors and no overflow at 390 pixels |
| Performance claims defensible | VERIFIED for the guarantee; the speed cost and the failed E23 confirmation are stated |

Reports in this folder: `CODE_AUDIT.md`, `REPRODUCIBILITY.md`, `METHODOLOGY_VALIDATION.md`, `PROFESSOR_DEFENSE.md`, `VISUAL_QA.md`, `DEMO_WALKTHROUGH.md`, `BACKTEST_STRESS_REPORT.md`, `RISK_REGISTER.md`, `EXECUTION_SUMMARY.md`, and `findings/` (50 findings in CSV).
