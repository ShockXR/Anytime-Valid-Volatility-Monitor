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

## Next steps, in order

1. Make the duplicate rule deterministic and mask nonpositive prices in all loaders; run `python scripts/reproduce.py --steps compare` and `--tier full --steps run,compare --from E18` in a copy.
2. Add a warning when the monitor stays blind, and a minimum calibration length with a justified default.
3. Make the benchmark detectors skip NaN or document that they do not.
4. Wrap CLI errors in short messages.
5. Run the pre-registered forward tests exactly as written: E25 from 26 November 2026 (`python experiments/e25_sequential_events.py status`), E24 after 1 June 2027.

## How to resume

```bash
git checkout claude/project-thread-tddcw4
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e . --no-deps
python -m pytest -q -p no:cacheprovider                         # expect 637 passed
python scripts/reproduce.py --steps compare                     # expect 465 identical
python scripts/trace_numbers.py                                 # expect 0 untraced
```

Hazard: `scripts/reproduce.py --steps run` overwrites committed result files. Run it in a copy.
