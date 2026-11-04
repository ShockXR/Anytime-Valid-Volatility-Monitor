# Round 4, agent R5: metamorphic and property tests

Files: `reports/round4/R5/test_round4_properties.py` (the tests), `reports/round4/R5/pytest_final.log` (the last run),
`reports/round4/R5/pytest_run4_with_audit.log` (the same file run together with `tests/test_audit_agent1.py`),
`reports/round4/R5/pytest_run1.log` to `pytest_run3.log` (earlier runs, kept for the record).

## 1. Summary in plain words

An example test checks one input and one expected answer. A property test checks a rule that must hold for
every input, for example "if you multiply all returns by 2, the alarms must not move". Such rules catch mistakes
that a single worked example cannot.

I read the monitor code (`monitor.py`, `detectors.py`, `confseq.py`, `cli.py`) and `docs/math.md` section 5
to confirm which rules are really part of the design, then wrote 28 test functions (56 test cases after
parameter expansion). Each case runs on random paths from five laws (Gaussian, Student t3, uniform, a two-point
law, Laplace), random lengths (600 to 2500 returns), random calibration lengths and alarm levels, with several
volatility regimes so that both the increase and the decrease side alarm. Seeds are 9,900,000 and up.

Result of the final run: **54 passed, 2 xfailed (expected failures that document real bugs), 0 failed, 28 s**
on one core. Run together with the audit file: 110 passed, 2 xfailed.

Two real bugs were found, both of low severity (they need absurd inputs or an infinite threshold). Three more
observations are listed below; none affects the false alarm guarantee. All the rules that the proof rests on
(expected betting factor exactly 1 under the exact null, statistic mean at most t, false alarm rate at most
m / A) held for every law tried, Gaussian or not.

The two bugs are marked `xfail(strict=True)` in the test file, so the suite turns red the moment they are fixed
and the marks can then be removed. No source file was changed.

## 2. Relations checked

| # | Relation | Test (prefix `test_`) | Result |
|---|---|---|---|
| 1 | Returns times k > 0: same alarms (plug-in, cs, c_down, band, margin, target width) | `scale_invariance_of_the_monitors[k]` (5 values of k) | pass (bit-identical for k a power of two) |
| 1 | Same for the Gaussian LR benchmark monitor | `scale_invariance_of_the_gaussian_lr_monitor[k]` | pass |
| 1 | Same for `betting_edetector` with m0 and tau^2 scaled by k^2, SR and CUSUM, both sides | `scale_invariance_of_the_known_baseline_detector[k]` | pass |
| 2 | r -> -r and random per-return sign flips: identical events and segments | `sign_invariance_global_and_per_return` | pass |
| 3 | update() one at a time = run() = random chunks (with empty chunks), 5 monitor kinds, final state too | `streaming_equals_batch_equals_random_chunks` | pass |
| 4 | NaN at a monitored step only shifts later indices | `nan_at_monitored_steps_only_shifts_later_events[3 kinds]` | pass |
| 4 | NaN inside a calibration window | (not tested) | not part of the design: the window closes after `cal_len` values including the NaN, so the baseline uses fewer returns (`docs/math.md` 5.6 item 5) |
| 5 | pickle / deepcopy mid-stream, continue, equals uninterrupted; original untouched | `pickle_or_deepcopy_mid_stream_then_continue_equals_uninterrupted[2]` | pass |
| 6 | Larger threshold never alarms earlier (first alarm), 4 monitor kinds | `larger_threshold_never_alarms_earlier` | pass |
| 6 | Same for the array detector, SR and CUSUM, both sides | `larger_threshold_never_alarms_earlier_for_the_array_detector` | pass |
| 6 | Wider null band or larger margin_z never alarms earlier (plug-in and cs) | `wider_null_band_or_larger_margin_never_alarms_earlier` | pass (proof of the monotone factor in the docstring) |
| 7 | E[one-step factor] = 1 under the exact null for every bet, both sides, 6 laws; below 1 with a band | `expected_one_step_factor_is_one_under_the_exact_null[6 laws]` | pass (true m0 by numerical integration; 400,000 draws; within 4 standard errors) |
| 7 | Gaussian LR factor: mean 1 for Gaussian data, far above 1 for t3 | `expected_gaussian_lr_factor_is_one_only_for_gaussian_data` | pass (documents why the benchmark is not valid for heavy tails) |
| 7b | Two-sided statistic (M_up + M_dn) / 2: mean at most t, and P(alarm within m) <= m / A, 4 laws, 20,000 paths | `two_sided_statistic_has_mean_at_most_t_and_false_alarm_rate_within_m_over_a[4]` | pass |
| 7c | Streaming monitor and `betting_edetector` are the same recursion (paths equal to 1e-10, same first alarm) | `streaming_monitor_equals_the_array_detector_on_one_segment` | pass |
| 8 | 1e-300, 1e300, subnormals, zeros, alternating huge and tiny, max float: statistics finite and >= 0 | `eprocess_stays_finite_and_nonnegative_on_extreme_inputs[4 kinds]` | pass for every stream that does not crash; the crash is bug R5-B-001 |
| 8 | A calibration window of huge returns must not raise | `calibration_window_of_huge_returns_does_not_crash` | **real bug R5-B-001** (xfail strict) |
| 8 | Gaussian LR monitor on the same streams: no NaN | `gaussian_lr_monitor_on_extreme_inputs_stays_finite_or_alarms` | pass (same crash skipped) |
| 8 | Array detector on x = 1e300, long zeros, subnormals | `array_detector_on_extreme_statistics` | pass |
| 8 | Array detector with x = +inf | `array_detector_treats_an_infinite_statistic_as_missing` | pass, but see observation R5-O-001 |
| 8 | threshold = inf means "never alarm" | `infinite_threshold_can_still_alarm_after_overflow` | **real bug R5-B-002** (xfail strict) |
| 9 | Restart after an alarm equals a fresh monitor started after the alarm (events, segments, state), 4 kinds | `restart_after_alarm_equals_a_fresh_monitor[4]` | pass |
| 9 | Restart after a degenerate (all zero) window equals a fresh monitor started at the restart point | `restart_after_a_degenerate_window_equals_a_fresh_monitor` | pass; the cs restart point is data dependent, see note; bookkeeping observation R5-O-002 |
| 10 | `python -m avmonitor simulate` (subprocess) prints the API's events, 3 parameter sets | `cli_simulate_subprocess_matches_the_python_api[3]` | pass |
| 10 | `avmonitor run --source bitstamp` on a small CSV prints the API's alarm times and sides, robust and gaussian | `cli_run_on_a_csv_matches_the_python_api[2]` | pass |
| + | c_down = c reproduces the default bit for bit | `c_down_equal_to_c_reproduces_the_default_exactly` | pass |
| + | Prefix property: run on r[:n] gives exactly the events with index < n | `events_before_a_cut_do_not_change_when_the_stream_is_cut` | pass |
| + | Permuting the first calibration window leaves the plug-in baseline unchanged | `permuting_a_calibration_window_does_not_change_the_plugin_baseline` | pass |
| + | list, tuple, Series, int array, int updates give the same events as a float array | `input_container_does_not_matter` | pass |

Notes on relations that needed care (the test comments say the same):

* NaN shift rule: an alarm index moves by the number of NaNs at or before it, but a segment's `start` moves only
  by the NaNs strictly before it, because `start` is set when the window closes, before the NaN at that position
  is seen. My first version of the test had this wrong; the code is right.
* Degenerate restart with `baseline="cs"`: the tau split window is extended up to `cal_max - (cal_len - n_tau)`
  values and then restarted, and the second attempt already holds some of the good returns, so the restart point
  is not `cal_max`. The test records where `_restart()` fired and starts the fresh monitor there. The plug-in path
  restarts exactly at `cal_max`, as documented.
* `scale_floor > 0` breaks scale invariance by design (a floor in return units) and is left out.

## 3. Real bugs

### R5-B-001 (low): OverflowError when a calibration window has a huge median absolute return

Minimal reproduction (from the installed package, numpy 2.5.3):

```python
from avmonitor import VolatilityMonitor, calibrate, GaussianLRMonitor
VolatilityMonitor(cal_len=2).run([1e300, 1e300, 1.0])        # OverflowError: (34, 'Numerical result out of range')
VolatilityMonitor(cal_len=4, baseline="cs").run([1e300] * 4)  # same
calibrate(np.full(10, 1e300))                                 # same
g = GaussianLRMonitor(cal_len=2); g.run([1.0, -1.0]); g.update(1e200)  # same (Python float through update only)
```

Where: `src/avmonitor/monitor.py` line 201 (`tau2 = (self.c * scale) ** 2` in `_plugin_stats`), line 239
(fallback in `_calibrate`), line 251 (`tau2_dn`), lines 280 and 287 (`_feed_cs`), line 407
(`GaussianLRMonitor.update`, only when `update` gets a Python float; `run` passes `np.float64` and gets `inf`);
`src/avmonitor/detectors.py` line 102 (`tau**2` in `calibrate`). `scale` is a Python float, and Python's `**`
raises `OverflowError` instead of returning `inf` (numpy would return `inf` with a warning). The threshold is
`c * scale` above about 1.34e154, so a window median above about 4.5e153 with c = 3.

Why it matters: the module docstring promises that a degenerate window is extended, falls back or restarts, and
`_usable()` already treats `inf` as unusable. The crash bypasses that path. The per-step `update()` of
`VolatilityMonitor` is safe (`r * r` gives `inf`, then `min(inf, tau2) = tau2`). Real returns never reach 1e154,
so this is a robustness defect, not a validity defect.

Suggested fix: compute the square with numpy, for example `tau2 = float(np.square(np.float64(self.c * scale)))`
(gives `inf`, which `_usable` rejects), in all five places of `monitor.py`; in `detectors.calibrate` use
`np.float64(tau) ** 2` and raise the existing `ValueError("degenerate calibration window ...")` when `tau**2` is
not finite; in `GaussianLRMonitor.update` use `z2 = float(np.float64(r / self.scale) ** 2)` or `r * r / scale**2`.

### R5-B-002 (low): threshold = inf is accepted but can still alarm

Reproduction:

```python
rng = np.random.default_rng(9900002)
m = VolatilityMonitor(threshold=float("inf"), cal_len=100)
m.run(np.r_[rng.standard_normal(100), np.full(2000, 1e6)])
m.events   # [Event(index=658, direction='up', ..., stat=inf)]
```

Where: `src/avmonitor/monitor.py` line 366 (`if mu >= 2 * self.A`) and line 415 (Gaussian LR). Audit A1-B-001
made `threshold` positive but kept `inf` allowed. Once the statistic overflows to `inf` (about 450 steps of
clipped returns at the largest up bet), `inf >= inf` is `True` and an alarm with `stat = inf` is recorded,
although the bound m / A = 0 promises no false alarm ever.

Suggested fix: reject a non-finite threshold in `__init__` (`if not (math.isfinite(self.A) and self.A > 0)`), or
test `math.isfinite(mu) and mu >= 2 * self.A`. The first is simpler and matches the audit's intent.

## 4. Observations (not bugs, worth a line in the docs)

* R5-O-001: `betting_edetector` treats `x = +inf` like a missing value (`np.isfinite` on the gain), so an
  untruncated squared return that overflowed carries no evidence. Validity is untouched (a zero bet is always
  allowed); detection is lost silently. The monitor never passes `inf` because it clips at tau^2 first. A
  `ValueError` for non-NaN non-finite `x` would make the choice explicit. `detectors.py`, loop in
  `betting_edetector`.
* R5-O-002: `_restart()` (`monitor.py` line 211) resets the buffer and the confidence sequences but not
  `_extended`, so the segment calibrated right after a restart is recorded with `extended=True` in `segments` and
  counted in `cal_stats["extended"]` although its window was not extended. Bookkeeping only; events are
  identical to a fresh monitor (test 9). Fix: `self._extended = False` in `_restart`.
* R5-O-003: with `baseline="cs"`, a window whose tau part is normal and whose second part is all zeros gives
  `m0 = m_dn = 0`, switches the decrease side off and does not extend the window; the plug-in path would extend
  and fall back. The state is valid and conservative (the increase side bets against the upper bound U), but the
  asymmetry is undocumented. Seen in the "normal then tiny / subnormal / zeros" streams of test 8.
* The Gaussian LR factor has an infinite expectation under t3 returns (the sample mean of 300,000 draws was
  above 1e187 with the code's cap of 700 on the log factor). This is the documented reason the benchmark is not
  valid for heavy tails; the test pins it.

## 5. What this does not show

* Random inputs with 8 to 10 paths per relation and five laws can miss rare paths; a relation that passed here
  is "not falsified", not proved. The proof-level checks (7, 7b) use the true baseline, so they say nothing about
  the plug-in baseline error, which `docs/math.md` 5.6 discusses and section 6 fixes.
* Scale invariance for k that is not a power of two holds up to floating point rounding; an alarm that sits
  within about 1e-12 of 2A could in principle move by one step. None did.
* The CLI `run` test uses a synthetic Bitstamp-format file and the in-process entry point with `DATA_ROOT`
  patched; `python -m avmonitor run` on the real data files, the Binance and Nifty loaders, and `--plot` were
  not exercised here (`tests/test_integration.py` covers the Binance loader).
* Not covered: `PeriodicVolatilityMonitor`, the guarded monitors, `MultiFrequencyMonitor`, `PanelMonitor`, the
  adaptive bet of `betting_edetector`, and the data cleaning functions.
* Extreme inputs were tested with finite thresholds 20 and inf only; no test feeds `inf` or `-inf` returns to the
  monitor (the audit covers `inf` in `calibrate`).
