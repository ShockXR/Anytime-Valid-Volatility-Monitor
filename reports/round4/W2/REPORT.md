# W2: clean rerun of E18 (landmark case studies)

Agent W2, 4 November 2026, 02:38 to 03:25 UTC. Everything ran in the clean worktree `/home/user/wt-repro`
(commit 4de8a32, `data/raw` symlinked to the main checkout's data). The main checkout's `results/` was not touched;
this folder holds the only files written there.

## 1. What was run

Preparation: `git status --short` in the worktree showed only the `seconds` column changes left by R1's E22 and E23
reruns; `git checkout -- results/` restored the committed copies, so every result file on disk was the committed one.

Command (exactly the E18 entry of `scripts/reproduce.py`), started 02:38:02 UTC, log `e18_run.log` in this folder:

```
cd /home/user/wt-repro && source /home/user/.venv-bapat/bin/activate
PYTHONPATH=/home/user/wt-repro/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python experiments/e18_landmark_events.py run figure context awake
```

The script's own `Pool(2)` was left as it is (two worker processes; `E18_CACHE` unset, so every BOCPD array was
recomputed from the raw data). Afterwards: `PYTHONPATH=$PWD/src python scripts/reproduce.py --steps compare` in the
worktree, plus a per-file byte comparison against `git show HEAD:<file>` that also records each file's write time.

## 2. Wall time

| Part | Started | Finished | Time |
|---|---|---|---|
| `run` (Nifty 223 s; crypto BOCPD, 6 jobs on 15 series, 1,564 s; 15 series and S2, 377 s) | 02:38:02 | 03:14:06 | 2,164 s |
| `figure` and `context` | 03:14:06 | 03:14:12 | 6 s |
| `awake` (E14's monitors, exploratory) | 03:14:12 | 03:21:14 | 422 s |
| whole command (exit code 0) | 02:38:02 | 03:21:14 | 43 min 12 s |

The original run took 1,940 s (with `nice -n 10` and a surviving Nifty BOCPD file) and `awake` 1,191 s
(docs/iterations/E18.md). R1's attempt was stopped at 27 minutes, before the crypto BOCPD jobs ended.

## 3. Compare output

`scripts/reproduce.py --steps compare`, run at 03:14:13 UTC after `run`, `figure` and `context`, and again at 03:21:19 UTC after `awake` (same output both times):

```
== compare
compare: 465 identical, 0 equal within rtol 1e-09, 0 equal except timing columns, 0 different, 0 missing; 0 files not in the manifest
summary: compare ok
```

(465 = every file in `results/MANIFEST.csv`.) `git status --short` in the worktree was empty afterwards: not one byte
of any tracked result file changed.

## 4. Per-file table

"Rewritten" is read from the file's modification time (a file the rerun never reached would also count as "identical"
in the compare step, so the time matters).

| File | Rewritten by this run? | Result vs committed | Size |
|---|---|---|---|
| `results/e18_events.csv` | yes, 02:38:02 (start of `run`) | identical bytes | 5,706 B |
| `results/e18_outcomes.csv` | yes, 03:14:06 | identical bytes | 108,821 B |
| `results/e18_alarms.csv` | yes, 03:14:06 | identical bytes | 216,098 B |
| `results/e18_thresholds.csv` | yes, 03:14:06 | identical bytes | 38,785 B |
| `results/e18_rates.csv` | yes, 03:14:06 | identical bytes | 13,185 B |
| `results/figures/e18_events.png` (not in the manifest) | yes, 03:14:11 | identical bytes (SHA-256 1b916236...) | 627,697 B |
| `results/e18_context.csv` | yes, 03:14:12 | identical bytes | 1,354 B |
| `results/e18_awake.csv` | yes, 03:21:14 | identical bytes | 19,065 B |

The script's own checks also passed: the e-detector's alarms on the Nifty stream and the ten H1/H2 2025 coin series
are identical to E13's saved alarms (all 11 `True`), and the 96 whole-sample (S3) thresholds equal E13's exactly.

Numbers in README.md (the "Landmark shocks" bullet, line 42) and in the E18 section of docs/experiment_report.md
(lines 2555 to 2998: 11, 104 and 172 minutes on Nifty; 0 of 15 crypto pairs within 24 hours; up alarms 4.5 to 6.9
hours before T; 15 to 207 minutes with a fresh calibration; 17 to 73 minutes for the k = 1440 scan; 15 and 17 of 18
for the awake monitors): none would change, because the files they are read from are byte for byte the committed ones.

## 5. Replay page check (`scripts/build_replay.py`)

The full build (N1 and C3, about 21 minutes) did not fit before the deadline, so the Nifty-only build was run instead:
`python scripts/build_replay.py --days N1` (03:21:29 to 03:24:52 UTC, 203 s, exit 0; log `replay_n1.log`). Its check
against E18 passed: on 4 June 2024 all seven methods' alarms from 7 days before to 24 hours after T are identical to
`e18_alarms.csv` (e-detector 2, BOCPD w = 250 and 375 2 each, truncated BOCPD 2 each, scan k = 120 3, k = 375 1) and the
outcomes equal `e18_outcomes.csv`. The page it wrote holds only N1, so it was not compared with the committed
two-day page; the worktree's `docs/demo/replay.html` was restored with `git checkout`. The C3 (BTC, 10 October 2025)
half of the replay check is PARTIAL: not rerun.

## 6. Verdict

**VERIFIED.** The full E18 command (`run figure context awake`) finished with exit code 0 in 43 minutes on a clean
copy, and all seven E18 result files and the figure are byte for byte the committed ones. The script's internal
checks against E13 (alarm lists, 96 thresholds) and E14 (24 alarm lists, counts 112 to 785) all report identical.
No number in README.md or the experiment report's E18 section changes.

Limits: two worker processes, as the script insists (no `--procs` flag, R1's low-severity note); the data were the
files already on disk (checked by R1 against `data/MANIFEST.csv`); the C3 half of the replay build, E19 p1 and the
forward tests were not run.
