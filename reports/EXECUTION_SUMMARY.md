# Execution summary and restart checkpoint (avmonitor)

3 November 2026. Four AI agents (Sonnet) and a coordinator worked on this repository. A container restart stopped two agents partway; their files were on disk and both agents were resumed and finished. Nothing was lost.

## What was completed

| Workstream | Output | State |
|---|---|---|
| Code audit | `CODE_AUDIT.md`, `REPRODUCIBILITY.md`, `tests/test_audit_agent1.py` (56 tests), `findings/A1_B.csv` (17) | complete (full reproduction tier not rerun) |
| Documents and methodology | `METHODOLOGY_VALIDATION.md`, `PROFESSOR_DEFENSE.md`, `checks/equation_vs_code_check.py`, `findings/A2_B.csv` (11); AI-authorship disclosure, test counts, E22 caveat | complete |
| Figures and demo | `VISUAL_QA.md`, `DEMO_WALKTHROUGH.md`, `screenshots/`, two figures regenerated, `findings/A3_B.csv` (10) | complete |
| Backtest and stress | `BACKTEST_STRESS_REPORT.md`, `RISK_REGISTER.md`, `stress/`, `findings/A4_B.csv` | B1 to B4 done, B5 partial |
| Coordinator | `MASTER_BACKTEST_REPORT.md`, this file, manifest rebuild, README command fix, regression runs | complete |

## Tests

| Run | Result |
|---|---|
| Before any edit | 581 passed, 0 failed, 7 warnings, 104 s |
| After all agents | 637 passed (581 + 56 new), 0 failed |
| Quick reproduction tier | experiment outputs byte-identical before and after the fixes; manifest rebuilt, then 465 of 465 files identical |
| Number tracer | 2,552 numbers traced in the paper, 0 untraced |
| Not run | full tier from E18 onward, real-data experiments, B5 for 137 of 140 data files |

## Changes made to code, documents and figures

Code (validation only, no change to any computed result): detectors reject negative statistics, nonpositive or NaN baselines and unknown kind on empty input; the repeated-window z-test no longer crashes on short input; the monitor rejects a threshold that is NaN or not positive, a calibration length below 2 and a truncation constant that is not finite. Docstring for the "down" grid corrected. Figures E1 and E4 regenerated from code constants. Documents: test counts, AI-authorship disclosure, E22 and E12 order of design, README command.

## Verified, unverified, failed

VERIFIED (OUR CHECK): false-alarm guarantee on fresh seeds (8 laws, 3,000 paths each), power and delay, sensitivity, 20 of 20 equation checks, replay page, quick reproduction tier.
PARTIALLY VERIFIED: real-data integrity (BTC 2024-01 to 2024-03 only), full reproduction (E1 to E17 by the authors).
FAILED: loaders keep an arbitrary duplicate and do not mask nonpositive prices (Bitstamp, Nifty); monitor accepts calibration lengths 2 to 10 and silently goes blind on constant or all-NaN streams; the two benchmark detectors die after one NaN.
UNVERIFIED: the mathematical proofs have not been checked by a human.

## Follow-up on 4 November 2026 (after the merge of PRs 1 to 3)

Items 1 to 4 of the old "Next steps" list are done, each with tests in `tests/test_input_hardening.py` (19 tests, all written for this follow-up):

1. **Loaders.** The Bitstamp loader now sorts with a stable sort, so a duplicated minute always keeps the first row read, and it reports `rows_read`, `duplicates` and `nonpositive_prices` in `df.attrs["checks"]`. The Bitstamp and Nifty loaders now mask bars with a close of 0 or below (open, high, low and close become NaN, so the returns on both sides are missing, never minus infinity), as the Binance loader already did. Fixes A4-B-001 and A4-B-002 (R-B9, R-B10).
2. **Monitor.** `VolatilityMonitor` warns when `cal_len` is below 100 (`MIN_ADVISED_CAL_LEN`; the scaled MAD of n returns has a relative standard error of about 1.17/sqrt(n), 12 percent at n = 100, and every experiment uses 120 or more), warns at the end of `run` when it never finished calibrating (an all-NaN or frozen feed, so "no alarm" would otherwise look like "no change"), warns when `null_band` is 1 or more (the decrease side is then off), and `update` rejects a non-number at once. Fixes A4-B-003 to A4-B-005 (R-B7). Calibration lengths 2 to 99 still run, with the warning.
3. **Benchmark detectors.** `gaussian_lr_edetector` and `page_cusum` now skip a missing value (the statistic stays where it was), as `betting_edetector` does. Their only callers are E4 and E11, which use simulated data without NaN; E4 was rerun with the new code and its outputs are byte-identical. Fixes A1-B-007.
4. **Command line.** `python -m avmonitor` prints one line starting with `error:` and exits with status 2 for a missing data file (with a hint to add `--download`), `--end` not after `--start`, a `--budget` outside (0, 1], and other bad input, instead of a traceback.

Checks after the follow-up (OUR CHECK, 4 November 2026): full suite 656 passed; `scripts/reproduce.py --tier quick --steps run,compare` and the E4 rerun give 465 of 465 result files identical; `scripts/trace_numbers.py` 0 untraced; a fuzz run of 400 hostile streams (NaN, infinities, 1e6 outliers, frozen feeds, tiny scales) through the monitor and both benchmarks raised no exception, gave no non-finite statistic, and every monitor that never calibrated warned. The suite also passes with the pinned versions in `requirements.txt`.

## Next steps, in order

1. Run the full tier from E18 in a copy (`--tier full --steps run,compare --from E18`); not run in this follow-up.
2. Run the pre-registered forward tests exactly as written: E25 from 26 November 2026 (`python experiments/e25_sequential_events.py status`), E24 after 1 June 2027.

## How to resume

```bash
git checkout main
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e . --no-deps
python -m pytest -q -p no:cacheprovider                         # expect 656 passed
python scripts/reproduce.py --steps compare                     # expect 465 identical
python scripts/trace_numbers.py                                 # expect 0 untraced
```

Hazard: `scripts/reproduce.py --steps run` overwrites committed result files. Run it in a copy.
