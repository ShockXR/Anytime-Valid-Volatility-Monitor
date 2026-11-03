# Reproducibility check: anytime-valid volatility monitor

**Plain-language summary (for a beginner)**

1. "Reproducible" here means: a stranger runs one command and gets the same files as the ones stored in the repository.
2. We ran the repository's own one-command tool (`scripts/reproduce.py --tier quick`) in a throw-away copy, so the stored results were never touched. Both quick experiments (E3 and E10d) came out byte-for-byte identical to the stored files.
3. The tool still ends with "FAILED" because three small bookkeeping files (`results/trace/*.csv`, which record where each number in the documents comes from) do not match the hashes in `results/MANIFEST.csv`. This is a stale-manifest problem, not a numerical one (finding A1-B-010). Fix: regenerate the traces and the manifest together after the last document edit.
4. The full tier (all 25 experiments, several hours) and the real-data experiments were not run here: the raw data are not in the repository and downloads belong to another workstream. Those parts are BLOCKED, not failed.
5. Randomness is seeded everywhere we looked; replaying a simulation with the same seed gives bit-identical alarms.

## Results

| Check | Result | Status |
|---|---|---|
| `python scripts/reproduce.py --tier quick --steps run,compare` in an isolated copy (Python 3.13, numpy 2.5.4, pandas 3.1.0), before the code fixes | 2 experiments run (E3 online truncation, E10d signal cap), 462 result files identical to `results/MANIFEST.csv`, 3 different, 0 missing, exit status 1 | VERIFIED |
| Same command repeated after the input-validation fixes in `src/avmonitor` | identical outcome: 462 identical, the same 3 different | VERIFIED |
| The 3 different files | `results/trace/paper_trace.csv`, `pitch_bapat_trace.csv`, `readme_trace.csv`. At the committed head their SHA-256 already differs from the manifest (the working tree has additionally been regenerated since). Not caused by any experiment | VERIFIED |
| `data` step (SHA-256 of raw files in `data/MANIFEST.csv`) and the real-data experiments (E5 to E8, E10 to E13 real parts, E19 and others) | `data/raw` absent; not run | BLOCKED |
| Full tier (`--tier full`) | not run (several hours) | UNVERIFIED |
| Test suite | 581 passed at the baseline; 637 passed after adding 56 audit tests and the fixes (see CODE_AUDIT.md) | VERIFIED |
| Pinned versions | `requirements.txt` pins numpy 2.5.3 and pandas 3.0.5; the check ran on numpy 2.5.4 and pandas 3.1.0 and still matched byte for byte for the two quick experiments | VERIFIED (two experiments only) |
| Python 3.11 (the version CI uses) | all `.py` files parse under 3.11; the suite was executed on 3.13 only | PARTIALLY VERIFIED |
| Seeds | library defaults call `default_rng()` unseeded in `simulate.py` and `design.py`, but every experiment and test passes an explicit seed or generator; the two unseeded calls in `experiments/e20_inference.py` are immediately overwritten with a saved state. The guard's rank fallback is nondeterministic only if `guard_fallback="rank"` and no `rank_seed` is given | VERIFIED by grep and replay tests |
| Undocumented manual steps | none found in `reproduce.py`: data download is an explicit flag (`--download`), the order of steps is encoded in the script | VERIFIED |

## Hazards for the person presenting

* Running `reproduce.py` in the repository itself overwrites the committed result files ("run" step). Always run it in a copy, as done here.
* Tracked files under `results/trace/` change whenever `scripts/trace_numbers.py` is run. The six test files that might write to the tree (trace, reproduce, replay, public API, examples, integration) were run one by one and wrote nothing, so the changes seen during this session come from running the tracer, not from pytest. The whole suite was not checked file by file.
* Benchmark detectors `gaussian_lr_edetector` and `page_cusum` stop alarming after a single NaN (A1-B-007): relevant only if real-data benchmark rows contain missing returns.

## How to repeat

```
cp -r <repo> /tmp/avm_copy && cd /tmp/avm_copy
python scripts/reproduce.py --tier quick --steps run,compare     # about 2 minutes
python -m pytest -q                                              # about 2 minutes
```
