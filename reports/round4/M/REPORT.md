# Round 4, agent M: repairs to the package

Branch `round4-mend` in the worktree `/home/user/wt-mend`. It starts from PR #4 (`origin/claude/project-thread-2lfy7y`, commit `cd3dd4f`). Another thread opened PR #4 with fixes for the same audit items. On the coordinator's instruction I rebased onto it. Where the two overlapped, I kept PR #4's behaviour and wording. I kept only what PR #4 does not cover.

## Summary in plain words

The validation pass found places where the package accepted bad input without complaint, or crashed instead of giving a clear message. PR #4 had already fixed most of them. On top of it, this branch:

- refuses an alarm level (`threshold`) that is not a finite number above 1. With an infinite level, the monitor used to raise an alarm anyway once its statistic overflowed.
- refuses a truncation constant `c` or `c_down` that is not a finite positive number, with a message that names the argument.
- adds `status()`, a small report of where the monitor stands. It answers the question "is my monitor blind?". It shows how many steps the monitor spent calibrating (it cannot alarm then). It shows how many calibration windows in a row it had to throw away. It shows whether the "volatility up" side can still win any bet.
- makes the third benchmark detector (`repeated_window_ztest`) skip missing values, as PR #4 did for the other two.
- makes the data loaders also blank out a price that is infinite or missing, and count such bars.
- removes a crash (`OverflowError`) when returns are absurdly large (above about 1e153).
- makes `update()` of the monitor subclasses reject a non-number at once, as PR #4 did for the base class.

No committed result changed. I checked this three ways (see "Reproduction check").

## One row per fix

| Defect | File | What changed | Test (in `tests/test_round4_mending.py`) |
| --- | --- | --- | --- |
| A4-B-003, R5-B-002 | `monitor.py` | `threshold` must be finite and above 1. Before, any number above 0 was accepted, and so was inf. The check is in `VolatilityMonitor.__init__`, so it covers `GaussianLRMonitor` and every subclass. | `test_threshold_must_be_finite_and_above_one`, `test_threshold_just_above_one_is_accepted`, `test_infinite_threshold_no_longer_alarms_after_overflow` |
| A4-B-003 | `monitor.py` | `c` must be a finite positive number. A non-number such as `None` now gives a `ValueError` that starts "c must be a finite"; before it gave a `TypeError` from `math.isfinite`. `c_down` must also be finite (inf was accepted). | `test_c_errors_name_c`, `test_c_down_must_be_finite_and_positive` |
| A4-B-004 | `monitor.py`, `periodic.py`, `guard.py` | `update()` rejects a non-number at once, with PR #4's message. This now holds for `GaussianLRMonitor`, for `PeriodicVolatilityMonitor` with a period or the CUSUM statistic, and for the guard's rank fallback. A float or `numpy.float64` passes unchanged, so the arithmetic is exactly as before. | `test_subclass_update_rejects_non_numbers_at_once`, `test_gaussian_lr_monitor_update_with_python_floats_gives_exactly_the_same_alarms` |
| A4-B-005 | `monitor.py`, `guard.py` | New `status()` method. New attributes `cal_steps` (steps fed to calibration) and `restarts_in_a_row` (windows thrown away since the last finished calibration). `status()` also reports `tau2_over_m_up`, `up_can_gain`, `down_active` and `blind`. It only reads: no alarm changes. | `test_status_shows_a_blind_stream_and_its_end`, `test_status_on_ordinary_data_and_before_any_data`, `test_status_shows_an_increase_side_that_cannot_gain`, `test_status_works_for_the_subclasses` |
| A1-B-007 | `detectors.py` | `repeated_window_ztest` skips missing values (NaN or infinite). The window holds the last observed returns, and the alarm index points at the return that completed the window. Without missing values, the old code runs unchanged. | `test_ztest_without_missing_values_matches_the_old_rule`, `test_ztest_skips_a_missing_value_as_if_it_were_dropped`, `test_ztest_with_fewer_observed_returns_than_a_window_never_alarms` |
| A4-B-002 | `data.py` | New `_mask_nonfinite`. All three loaders set the prices of a bar to NaN when its close is infinite or missing, and report the count as `nonfinite_prices` in `df.attrs["checks"]`. It runs before PR #4's `_mask_nonpositive`, so no bar is counted twice. In `load_binance` it runs before the NaN rows for missing minutes are added, so those are not counted. | `test_bitstamp_masks_infinite_and_missing_closes_and_counts_each_bar_once`, `test_nifty_masks_an_infinite_close`, `test_binance_masks_an_infinite_close_but_does_not_count_missing_minutes` |
| R5-B-001 | `monitor.py` (`_plugin_stats`, `_calibrate`, `_feed_cs`, `GaussianLRMonitor.update`), `detectors.py` (`calibrate`) | `(c * scale) ** 2` on a Python float raised `OverflowError` for huge returns. A helper `_sq` returns `x ** 2` exactly as before, and `inf` on overflow. The rule for unusable calibration windows then rejects `inf`. | `test_calibration_on_huge_returns_takes_the_degenerate_path_not_an_overflow_error`, `test_squares_of_ordinary_values_are_exactly_as_before` |
| A3-B-002 (tests only) | none | PR #4's CLI handler already catches these errors. I added tests for the paths `test_input_hardening.py` does not reach: an empty file, HTTP 404, no network, a missing column, a bad date, and a parameter the monitor rejects. | `test_cli_data_errors_give_one_line_and_exit_code_two`, `test_cli_bad_date_gives_one_line_and_exit_code_two`, `test_cli_monitor_parameter_error_gives_one_line` |

Also changed:
- `CHANGELOG.md`: one new block under 0.3.0.
- `docs/api.md`: rebuilt with `scripts/build_api_docs.py`, so it lists `status()` and the new `update` text.
- the docstrings of the changed functions.

`competitors.py` and `benchmarks_robust.py` needed no change. I read every detector there, and each one already skips or carries NaN through `np.isfinite` masks (for example `page_cusum_paths`, `sr_log_paths`, `first_crossing`, the BOCPD filters and `hkz_paths`).

## Dropped because PR #4 already covers it

These were in my first version (local branch `round4-mend-pre-pr4`). They are not on `round4-mend`:

- The Bitstamp stable sort, the first row kept for a duplicated minute, and the duplicate count.
- Masking closes at or below zero in the Bitstamp and Nifty loaders.
- An error for `cal_len` below 20. PR #4 warns below 100 instead.
- An error for `null_band` of 1 or more. PR #4 warns instead.
- My warning after 3 calibration restarts in a row (`warn_after_restarts`, `BlindMonitorWarning`). PR #4 warns at the end of `run` when the monitor never calibrated. Keeping both would give two warnings for the same all-NaN stream, so I kept only `status()` and the counters.
- Turning the input of `VolatilityMonitor.update` into a float. PR #4's code for this is kept as it is.
- NaN skipping in `page_cusum` and `gaussian_lr_edetector`, and my `skip_nan=False` option that gave the old rule back. PR #4 always skips, with no option, so `repeated_window_ztest` now does the same.
- My rewrite of the CLI (date checks and a separate message for each error type). PR #4's handler covers the same cases.

## Reproduction check

| When | Command (in the worktree, `PYTHONPATH=/home/user/wt-mend/src`) | Files compared | Result |
| --- | --- | --- | --- |
| Before, at 4de8a32 | `python scripts/reproduce.py --steps compare` | 465 | 465 identical |
| After, first version | `python scripts/reproduce.py --tier quick --steps run,compare` | 465 (E3 and E10d rerun) | 465 identical, `results/` clean |
| After, rebased branch | same | 465 (E3 and E10d rerun) | 465 identical, `results/` clean |

The quick tier reruns only E3 and E10d. So I added two direct checks of the changed code. Each runs the old package (main checkout, 4de8a32) and the new one (this branch) on the same input and compares every output byte for byte. The scripts are in `reports/round4/M/checks/`.

1. **Loaders** (`loader_hashes.py`). I compared every column of the frames the loaders return for:
   - the whole Bitstamp file, in both `filled` modes, for all rows and for July to December 2025;
   - Nifty from April 2017 to October 2026, all bars and the May to June 2024 regular sessions;
   - each of the five Binance coins from January 2024 to December 2025.

   Result: 11 of 11 frames identical. Only the new keys in `checks` differ.
2. **Monitors and benchmarks** (`equiv.py`). I ran seven monitor settings: plug-in with band 0 and with band 0.25, Gaussian LR, the confidence-sequence baseline, the same with `c_down = 1.5`, periodic, and guarded with the rank fallback. I ran each on three real series:
   - BTCUSDT, June to December 2025;
   - Nifty, 2024;
   - Bitstamp, July to December 2025, with `filled="zero"`.

   The returns were first divided by the intraday scale, as the CLI does. I compared every alarm (index, side, statistic as `repr`) and a hash of every calibration record. I also ran E4's benchmark detectors (Gaussian LR on both sides, Page CUSUM and the 60-step z-test) on E4's own simulated paths: seeds 100 to 104, the first 30 of 300 paths, the full 20,000-step horizon.

   Result: 26 of 26 identical, both before and after the rebase.

PR #4 reports a byte-identical E4 rerun. I did not rerun E4 myself, because it uses four cores (`Pool(4)`). I did not rerun any other full-tier experiment either.

## Real data: duplicates and bad prices

I counted on the raw files in `data/raw/` (download of 4 November 2026), with `checks/raw_counts.py`. The output is in `checks/raw_counts.csv`.

| Source | Files | Rows | Duplicate timestamps | Out of order | Close at or below 0 | Close not finite | Open, high or low at or below 0, or not finite |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Bitstamp BTC/USD | 1 | 957,844 | 0 | 0 | 0 | 0 | 0 |
| Nifty 50 (2017 to October 2026) | 19 | 886,548 | 0 | 0 | 0 | 0 | 0 |
| Binance, five coins | 170 | 7,401,600 | 0 | 0 | 0 | 0 | 0 |

So neither the duplicate rule nor any of the price rules can fire on these files, and no returns series that an experiment uses can change. For that reason, no old behaviour had to stay reachable through an argument.

One caution: the Bitstamp mirror file is a rolling "latest" file. Today's copy has 957,844 rows, runs to 3 November 2026, and has a different SHA-256 from the copy in `data/MANIFEST.csv` (956,327 rows). Agent R1 checked that the first 956,327 rows of today's copy hash exactly to the committed SHA-256 (`reports/round4/R1/REPORT.md`).

## Property tests of agent R5

With this branch, `reports/round4/R5/test_round4_properties.py` (moved to `tests/test_round4_properties.py` by the coordinator) needed three edits, which the coordinator made:

- `test_calibration_window_of_huge_returns_does_not_crash` (R5-B-001): now passes, so its xfail mark was removed.
- `test_eprocess_stays_finite_and_nonnegative_on_extreme_inputs` (4 cases): its last line expected the old crash. Nothing crashes now, so it asserts that nothing crashed.
- `test_infinite_threshold_can_still_alarm_after_overflow` (R5-B-002): `VolatilityMonitor(threshold=math.inf)` now raises `ValueError`, so the test checks that.

## Test suite

Command: `PYTHONPATH=$PWD/src OMP_NUM_THREADS=1 python -m pytest -q -p no:cacheprovider`. Result in the worktree: **706 passed** in 140 s. That is PR #4's 656 plus the 50 new tests in `tests/test_round4_mending.py`.

A note for anyone who reruns in a worktree. If `data/raw` is a symbolic link to another checkout's data, `tests/test_public_api.py` (which rebuilds `docs/api.md` and compares it) sees one default path printed differently. So I ran the suite and `build_api_docs.py` with the link removed. In a normal checkout, where `data/raw` is a real folder, this does not arise.

## Not done

- No full-tier experiment was rerun. E4 needs four cores, the real-data experiments run long, and the Bitstamp file on disk differs from the manifest copy anyway. The direct old-against-new checks above cover the changed code on real data instead.
- The old behaviour of the three benchmark detectors, where one NaN stopped them for good, can no longer be switched back on. That was PR #4's choice. No committed caller passes NaN, so no result depends on it.
- `cal_len` below 100 and `null_band` of 1 or more give warnings, not errors (PR #4's choice). A calibration window of 2 returns is still accepted.
- The loaders report bad prices but do not repair them: they leave a NaN for the monitor to skip.
- The restart warning I first wrote (`warn_after_restarts`) is not on the branch. `status()` covers that need, together with PR #4's warning at the end of `run`.

Note from the coordinator: the agent's file-writing tool refused to create this report, so the coordinator saved the agent's text here, with the R5 section updated to say the three test edits were made.
