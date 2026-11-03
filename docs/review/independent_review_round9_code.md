# Independent review, round 9: the package code

**Reviewer role.** A code auditor who did not write the package. I did not take the docstrings on trust: every claim below was tested by running the code, and every number comes from a script in `docs/review/round9_scripts/` that can be rerun.

**Scope.** In priority order:

1. `src/avmonitor/monitor.py` (`VolatilityMonitor`: calibration, τ and m₀, restart after an alarm, the plug-in and confidence-sequence paths).
2. `src/avmonitor/periodic.py` (`PeriodicVolatilityMonitor`, `restarted_sr_bound`, the simulation engine): period boundaries, restarts, the known, cs and plug-in modes.
3. `src/avmonitor/detectors.py` (`betting_edetector`, `calibrate`, `two_sided_alarm`, `Baseline`): predictability of every bet, the range [0, B] and the bet caps, the threshold used against the one documented, the Shiryaev-Roberts recursion, mixture weights, both sides.
4. `src/avmonitor/confseq.py` (`eb_confseq`, `StreamingEBConfSeq`, `cs_baseline`): predictable bets, one-sided ranges, the split of δ.
5. `src/avmonitor/multifreq.py`, `src/avmonitor/panel.py` (`ebh`), `src/avmonitor/cli.py`, and the public functions of `running.py`, `limits.py` and `design.py` as time allowed; `data.load_bitstamp_csv` for the handling of missing returns.

Cross-cutting checks: batch functions against streaming classes on the same data, missing values (the documents say never zero-fill), Y exactly 0 or exactly B, integer types, "first m steps" off by one, state that leaks across a restart, and Monte Carlo checks of the stated bounds through the streaming classes as a user would call them.

**Read first.** README (installation, use and key results), `docs/api.md`, `docs/math.md` sections 1 to 6, `docs/theory_per_period.md`, and the earlier reviews in `docs/review/` (round 1 audited the same modules at an early commit; its m17 is the origin of finding M1 below).

**State of the tree.** Reviewed at commit `949c1d5` on branch `claude/project-thread-ytr1da` and rechecked at `9c5786d`, the head when this report was written. The last commit that changed `src/avmonitor/` is `b99bd3f`; the commits after it touch experiments, results, documents and the paper only, so every line number below holds at all three commits, and `proposed_fixes.patch` applies to `9c5786d` (`git apply --check`). Other agents had uncommitted files in the working tree while I worked (E22, E23 and the paper); I did not touch them and no file in scope was modified.

**What I ran.** All scripts are in `docs/review/round9_scripts/` (copies of the scratchpad originals). Every run used `nice -n 10` and at most two processes.

* `r9_nan_batch.py`: missing values in `betting_edetector` against the streaming monitor (M1).
* `r9_batch_vs_stream.py`: batch functions against streaming classes, alarm by alarm, on 150 seeded paths per law (Gaussian and Student t3), with and without 10% missing returns; it also runs the proposed M1 fix.
* `r9_small_checks.py`: e-BH against a brute-force implementation, predictability of every bet, `eb_confseq` against `StreamingEBConfSeq`, Y at 0 and at B, the calm law Y = 1 (off by one), input validation, integer types, `c_down` with a known baseline.
* `r9_mc_null.py`: Monte Carlo under the null through `PeriodicVolatilityMonitor(baseline="known")` and `VolatilityMonitor(baseline="cs")`, i.i.d. Gaussian and t3 returns, bold play (the adaptive two-point null of E16), per-day rates with `period=390`, and the calm-then-Gaussian stream.
* `r9_multifreq.py`, `r9_panel.py`, `r9_running.py`: the same kind of checks for `MultiFrequencyMonitor`, the panel monitor and the three monitors of `running.py`.
* `r9_limits.py`: spot checks of `limits.py` (`sr_alarm_times`, `collapse_floor`, `floor_time`, `max_null_mass`).
* `r9_bitstamp.py`: the zero-volume minutes of the Bitstamp file (m5); it needs the raw file, which is not committed.
* `proposed_fixes.patch` and `round9_regression_tests.py`: every diff of this report as one patch, and regression tests for it (to be added as `tests/test_review_round9.py`; the file name in `round9_scripts/` keeps pytest from collecting it, since every test fails on the unpatched package). I applied the patch to a copy of the package outside the repository (never to `src/`), ran the full test suite on that copy with the new tests, and reran the reproductions against it (section 4).

## 1. Verdict in brief

**One Major finding, no Critical one.** No bug I found makes a stated guarantee of the streaming monitors fail. Under the null, every Monte Carlo rate that a document bounds (first day after a start, every day with `period=390` and a known baseline, the confidence-sequence baseline, the multi-frequency, panel and `running.py` monitors) is below its bound, most of them far below; bold play, the sharpest adversary known for this detector, comes close with one bet per side: 4.34% on the first day (0.87 of the bound) and, with daily restarts, between 4.28% and 4.93% on each of days 1 to 5 against 5% (section 5). This adversary nearly attains the bound of Theorem 1; no rate is above it. The batch functions and the streaming classes give identical alarms on complete data, alarm by alarm, for the known, plug-in and confidence-sequence baselines.

**M1** is in the batch API: `betting_edetector` (`detectors.py:241` to `252`) does not skip a missing value, as its docstring and `docs/math.md` 5.6 say it does. It sets the bet factor to 1 but still adds 1 to the Shiryaev-Roberts statistic, which is the update for an observation with Y = 1 exactly. So a stream with no data at all raises an alarm at step 2A, a degenerate calibration window turns into a deterministic alarm at step 2A, and on a 24-hour grid with session gaps the batch detector alarms on 89.5% of null paths against 32.3% for the streaming monitor on the same data, above the 60% that Theorem B allows for the observed returns. The fix changes four lines.

**Minor findings** (m1 to m7): `running.two_sided_sr` is silently disabled by one missing value (the failure mode of round 1, m17); `calibrate` returns τ = 0 and m₀ = 0 for a degenerate window without an error; bets that make a factor negative are accepted by the streaming monitors and by `betting_edetector`'s `lam_cap` and `adaptive_weight`; the known-baseline mode ignores `c_down` silently; the Bitstamp loader passes the provider's filled minutes on as zero returns (0.30% of the minutes of July to December 2025) [REAL]; `two_sided_alarm` returns 2,147,483,647 instead of -1 for alarm arrays of different integer types; `design.sr_steps_constant_factor` is one step too large at exact integer boundaries (no reported number is affected).

| Area | Verdict |
|---|---|
| `monitor.py` | Correct. Calibration, τ split, δ split, restart and the 2A rule checked; no state leaks across a restart |
| `periodic.py` | Correct. Known mode equals the engine and the batch detector; per-day rates below m / A; silently ignores `c_down` (m4) |
| `detectors.py` | Bets predictable, factors positive for valid inputs, weights sum to 1, threshold 2A per side as documented. Missing values not skipped (M1); degenerate `calibrate` (m2); unchecked `lam_cap`, `adaptive_weight` (m3); dtype (m6) |
| `confseq.py` | Correct. Predictable bets, one-sided ranges, running sums identical to the vectorised version, δ split as documented |
| `multifreq.py` | Correct. Certified value never above the next full value; vectorised equals streaming with missing returns |
| `panel.py` | Correct. e-BH right, streaming equals `run_block`, null rate below α |
| `cli.py` | Correct threshold; passes Bitstamp's filled minutes as zeros (m5, in `data.py`) |
| `running.py` | Monitors correct and below their bounds; `two_sided_sr` disabled by one NaN (m1) |
| `limits.py`, `design.py` | Spot checks correct (`sr_alarm_times`, `collapse_floor`, `floor_time`, `max_null_mass`, `decrease_speed_limit`); one off-by-one at exact boundaries in `design.py` (m7) |

## 2. Major

### M1. `betting_edetector` counts a missing value as an observation with Y = 1, so missing data alone raise alarms

*Location.* `src/avmonitor/detectors.py:241` to `252` (the loop), against its docstring at `:197` to `:200` ("Missing values (NaN) in ``x`` are skipped: the bet factor is 1 at that step"). The same claim is repeated in `src/avmonitor/benchmarks_robust.py:152` to `153` ("Missing values (NaN) leave the statistic unchanged, as in ``detectors.betting_edetector`` and ``monitor.VolatilityMonitor``"). The documented rule is `docs/math.md` 5.6 item 5 ("The monitor skips a missing return") and 6.3 ("A missing return leaves them unchanged").

*What is wrong.* At a NaN step the code sets the gain to 0, so the factor is 1, and then applies the Shiryaev-Roberts update `S = f * (S + 1.0)`. With f = 1 this is `S + 1`: a new chain is started and every statistic grows by 1, exactly as for an observation with Y = 1. The streaming monitors return before the update at a NaN step (`monitor.py:348` to `349`, `periodic.py:444` to `445`), and so do `benchmarks_robust.sr_log_paths`, `multifreq.component_paths` and the panel e-process. The fix of round 1 (m17) replaced "one NaN disables the detector" by "a NaN is imputed as Y = 1", which is not what the docstring says.

What this breaks:

1. **Alarms without data.** A stream with no observation at all alarms at step 2A - 1 (0-based), on both sides.
2. **A degenerate calibration gives a deterministic alarm.** `calibrate` on a window with more than half zero returns returns τ = 0 and m₀ = 0 without an error (m2). Then every Y = 0/0 is NaN, every step adds 1, and `betting_edetector(..., "up")` alarms at step 2A - 1 whatever the data.
3. **The guarantee in observed returns fails.** Theorem B bounds the probability of a false alarm within the first m steps; with the monitor's skip rule, m counts observed returns. For the batch function a step is a position in the array, missing or not. A series that starts with 2A - 1 missing values (an asset that starts trading later on a common clock) alarms at its first observation, with probability 1, even when that observation has Y = 1. On a realistic layout, a session market on a 24-hour minute grid (375 observed minutes and 1065 NaN per day, i.i.d. Gaussian returns, true baseline, A = 375 / 0.05 = 7500 as the CLI uses for Nifty), over 12 days (4500 observed returns) the batch detector alarmed on 89.5% of 400 null paths, after a median of 1703 observed returns, while the streaming monitor on the same data alarmed on 32.3%; Theorem B with m = 4500 observed returns allows 60% [SYNTHETIC].
4. **Batch and streaming disagree** on every path with missing values (150 of 150 in the comparison below), while they agree exactly on complete data.

What it does not break: with m counted as array positions (missing ones included) the bound m / A still holds, because a factor of 1 is a valid bet. So no result computed on complete data is affected. I found no committed experiment that passes NaN to `betting_edetector` (E4, E4c, E9, E13 and E15 use complete simulated data or calibrated arrays without gaps); the risk is in user code that applies the documented array API to real returns with gaps, as `examples/01_quickstart.py` shows how to do.

*Reproduction.* `r9_nan_batch.py`:

```
all-NaN stream of 15610 steps: batch alarm (up, down, two-sided) = 15599 15599 15599
all-NaN stream: streaming monitor events = 0
2A-1 NaN then one observation Y = 1: batch up alarm at 15599 statistic 15600.0
same stream, streaming monitor: events = 0 S_up mean = 1.0
5 NaN then Y = 1: sr_log_paths M_6 = 1.0  betting_edetector M_6 = 6.0
[SYNTHETIC] session market on a 24-hour grid, 12 days (4500 observed returns), A = 7500, 400 paths:
  batch betting_edetector: P(alarm) = 0.895; median observed returns before the alarm 1703
  streaming monitor:       P(alarm) = 0.323
  Theorem B bound with m = observed returns: 0.600
```

```
degenerate window: calibrate returns Baseline(scale=0.0, tau=0.0, m0=0.0, m0_up=np.float64(0.0), m0_down=1e-300, var_x=0.0, n_cal=500)
  betting_edetector(up, A = 50) alarm at 99 (= 2A - 1, whatever the data)
  down side raises: upper bound must exceed m0
```

`r9_batch_vs_stream.py`, comparison D (known baseline, restart after each alarm, 10% of returns missing, A = 60): `150 paths differ with betting_edetector, 0 with the M1 fix` for both Gaussian and t3 returns.

*Proposed fix.* Keep the statistic unchanged at a missing step (no bet and no new chain), for both recursions and for the adaptive bettor. The hunks below are copied from `proposed_fixes.patch` (line numbers are those of the whole patch); the two `lam_cap` lines in the first hunk belong to m3:

```diff
--- a/src/avmonitor/detectors.py
+++ b/src/avmonitor/detectors.py
@@ -194,11 +205,15 @@ def betting_edetector(
     start : int
         Monitoring starts at this index (earlier entries are ignored).
 
-    Missing values (NaN) in ``x`` are skipped: the bet factor is 1 at that step
-    (valid when missingness does not depend on the missing value, see
-    ``docs/math.md`` 5.6). Before November 2026 a single NaN silently disabled
-    the detector (independent review, m17).
+    Missing values (NaN) in ``x`` are skipped: the statistic is left unchanged at
+    that step (no bet and no new chain), as in ``VolatilityMonitor`` (valid when
+    missingness does not depend on the missing value, see ``docs/math.md`` 5.6).
+    Before November 2026 a single NaN silently disabled the detector (independent
+    review, m17), and until review round 9 (M1) a NaN step still added 1 to the
+    Shiryaev-Roberts statistic.
     """
+    if not 0 <= lam_cap <= 1 or not 0 <= adaptive_weight <= 1:
+        raise ValueError("lam_cap and adaptive_weight must lie in [0, 1]")
     x = np.atleast_2d(np.asarray(x, float))[:, start:]
     reps, n = x.shape
     y = x / m0
@@ -239,7 +254,8 @@ def betting_edetector(
     lam_path = np.empty(n) if (return_path and lam_adapt is not None) else None
     for t in range(n):
         g = sign * (y[:, t] - 1.0)  # gain per unit bet
-        g = np.where(np.isfinite(g), g, 0.0)  # missing value: no bet
+        ok = np.isfinite(g)  # missing value: no bet and no new chain
+        g = np.where(ok, g, 0.0)
         f = 1.0 + grid[None, :] * g[:, None]
         if lam_adapt is not None:
             f = np.concatenate([f, (1.0 + lam_adapt[:, t] * g)[:, None]], axis=1)
@@ -247,9 +263,9 @@ def betting_edetector(
                 lam_path[t] = lam_adapt[0, t]
         with np.errstate(over="ignore"):
             if kind == "SR":
-                S = f * (S + 1.0)
+                S = np.where(ok[:, None], f * (S + 1.0), S)
             elif kind == "CUSUM":
-                S = f * np.maximum(S, 1.0)
+                S = np.where(ok[:, None], f * np.maximum(S, 1.0), S)
             else:
                 raise ValueError("kind must be 'SR' or 'CUSUM'")
         M = S @ weights
```

With the fix, `M_t - (number of observed steps up to t)` is a supermartingale, so Theorem B holds with m counted in observed returns, as for the streaming monitor. On complete data the output is unchanged (the `np.where` selects the old value everywhere). With the patch applied to a copy of the package the full test suite passes (section 4), the comparison D above gives 0 differing paths, and the new test `test_m1_batch_and_streaming_agree_with_missing_returns` checks it.

## 3. Minor

### m1. `running.two_sided_sr` is disabled by a single missing value

*Location.* `src/avmonitor/running.py:354` to `369`. *What is wrong.* A NaN in `x` makes every later statistic NaN, so the function never alarms again: the failure mode of round 1, m17, for this vectorised route A detector. `RunningTauMonitor`, its streaming counterpart, skips missing returns (`running.py:711` to `712`). E14 calls it on complete simulated data, so no reported number is affected. *Reproduction* (`r9_small_checks.py`, item 9): `two_sided_sr first (up, down) alarm: complete data (8, -1); with x[0] = NaN (-1, -1)`. *Fix* (in `proposed_fixes.patch`):

```diff
--- a/src/avmonitor/running.py
+++ b/src/avmonitor/running.py
@@ -352,7 +352,8 @@ def two_sided_sr(x, tau2, m_up, m_dn, levels, up_grid=DEFAULT_UP_GRID, down_frac
     out = {L: [np.full(reps, -1), np.full(reps, -1)] for L in levels}
     top = levels[-1]
     for t in range(n):
-        xt = x[:, t]
+        ok = np.isfinite(x[:, t])  # a missing return leaves both statistics unchanged, as in RunningTauMonitor
+        xt = np.where(ok, x[:, t], 0.0)
         yu = xt / m_up[:, t]
         on = (m_dn[:, t] > 0) & (tau2[:, t] > m_dn[:, t])
         md = np.where(on, m_dn[:, t], 1.0)
@@ -362,11 +363,11 @@ def two_sided_sr(x, tau2, m_up, m_dn, levels, up_grid=DEFAULT_UP_GRID, down_frac
         fd = 1.0 + lam * (1.0 - (xt / md)[:, None])  # exactly 1 at an off step
         with np.errstate(over="ignore"):
             if kind == "SR":
-                Su = fu * (Su + 1.0)
-                Sd = fd * (Sd + on[:, None])  # an off step starts no chain
+                Su = np.where(ok[:, None], fu * (Su + 1.0), Su)
+                Sd = np.where(ok[:, None], fd * (Sd + on[:, None]), Sd)  # an off step starts no chain
             else:
-                Su = fu * np.maximum(Su, 1.0)
-                Sd = np.where(on[:, None], fd * np.maximum(Sd, 1.0), Sd)
+                Su = np.where(ok[:, None], fu * np.maximum(Su, 1.0), Su)
+                Sd = np.where((on & ok)[:, None], fd * np.maximum(Sd, 1.0), Sd)
         Mu = Su.mean(axis=1)
         Md = Sd.mean(axis=1)
         for L in levels:
```

### m2. `calibrate` (plug-in) returns τ = 0 and m₀ = 0 for a degenerate window without an error

*Location.* `src/avmonitor/detectors.py:100` to `108`. *What is wrong.* When more than half of the calibration returns are zero (a frozen feed, a coarse tick), the scaled MAD is 0, so τ = 0 and m₀ = 0, and `calibrate` returns them. The streaming monitor extends such a window (`monitor.py:224` to `234`) and `cs_baseline` raises ("tau is zero on some rows"), but the plug-in `calibrate` passes the problem on: `betting_edetector(x, 0.0, "up")` then sees Y = 0/0 = NaN at every step and, through M1, alarms at step 2A - 1 on any data; the decrease side raises "upper bound must exceed m0". *Reproduction* (`r9_nan_batch.py`, part 5, output quoted under M1): 300 zeros then 200 Gaussian returns of scale 10⁻³ give `Baseline(scale=0.0, tau=0.0, m0=0.0, ...)`, and the up detector at A = 50 alarms at index 99. `calibrate_panel` already catches `ValueError` and marks such an asset inactive, so raising is the consistent choice. *Fix:*

```diff
--- a/src/avmonitor/detectors.py
+++ b/src/avmonitor/detectors.py
@@ -101,6 +101,9 @@ def calibrate(
     tau = max(c * scale, floor)
     x = np.minimum(r_cal**2, tau**2)
     m0 = float(x.mean())
+    if not (tau > 0 and m0 > 0):
+        raise ValueError("degenerate calibration window (scale or truncated mean 0, for example more than half "
+                         "zero returns); pass floor > 0 or use a longer window")
     v = float(x.var(ddof=1))
     se = np.sqrt(v / len(x))
     up = (m0 + margin_z * se) * (1 + null_band)
```

### m3. Bets that make a factor negative are accepted

*Location.* `src/avmonitor/monitor.py:141` to `142` (`up_grid`, `down_fractions`, and `null_band` at `:140`), `multifreq.py:357` to `358`, `running.py:655` to `656` and `:1056` to `1057` (no check at all); `detectors.py:164` to `165` and `:229` to `234` (`adaptive_weight` and `lam_cap` unchecked; the fixed grid is checked at `:222`). `PanelEProcess` already refuses bets outside [0, 1) (`panel.py:165` to `170`). *What is wrong.* Validity needs every factor to be nonnegative (`docs/math.md` 5.2: λ ∈ [0, 1) up, λ ∈ [0, 1/(B - 1)) down). `VolatilityMonitor(up_grid=(1.5,))` or `down_fractions=(1.5,)` runs and produces negative statistics; `lam_cap > 1` lets the adaptive bet exceed the largest admissible bet; `adaptive_weight > 1` gives the grid a negative weight, so the average is no longer an e-detector; a negative `null_band` bets against a value below the baseline. These are input errors, not defaults, but the guarantee depends on them and nothing stops them. *Reproduction* (`r9_small_checks.py`, item 6): `S_dn [-1.09322375]` after two steps with `down_fractions=(1.5,)`; with `lam_cap=3` the adaptive bet reaches 2 and the statistic reaches -6.71e+07. *Fix:* a helper `check_bets` in `detectors.py` called by the four constructors, a range check on `null_band` (below), and the check of `lam_cap` and `adaptive_weight` shown in the first hunk of M1:

```diff
--- a/src/avmonitor/detectors.py
+++ b/src/avmonitor/detectors.py
@@ -116,6 +119,14 @@ DEFAULT_UP_GRID = (0.02, 0.05, 0.1, 0.2, 0.35, 0.5)
 DEFAULT_DOWN_FRACTIONS = (0.1, 0.25, 0.5, 0.75, 0.95)
 
 
+def check_bets(up_grid, down_fractions) -> None:
+    """Raise unless every increase bet and every decrease fraction lies in [0, 1), so that every factor is positive."""
+    for name, v in (("increase bets", up_grid), ("decrease fractions", down_fractions)):
+        v = np.asarray(v, float)
+        if v.size == 0 or not np.all(np.isfinite(v)) or np.any(v < 0) or np.any(v >= 1):
+            raise ValueError(f"{name} must be finite and lie in [0, 1)")
+
+
 @dataclass
 class DetectorResult:
     """Output of the array detectors: first alarm index per path (-1 if none) and optional paths."""
--- a/src/avmonitor/monitor.py
+++ b/src/avmonitor/monitor.py
@@ -76,7 +76,7 @@ from dataclasses import dataclass
 import numpy as np
 
 from .confseq import StreamingEBConfSeq
-from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
+from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID, check_bets
 from .robust import MAD_NORMAL_SCALE
 
 
@@ -140,6 +140,9 @@ class VolatilityMonitor:
         self.band = null_band
         self.up_grid = np.asarray(up_grid, float)
         self.down_fractions = np.asarray(down_fractions, float)
+        check_bets(self.up_grid, self.down_fractions)
+        if not null_band >= 0:
+            raise ValueError("null_band must be nonnegative")
         self.scale_floor = scale_floor
         if baseline not in ("plugin", "cs"):
             raise ValueError("baseline must be 'plugin' or 'cs'")
```

and the same `check_bets(self.up_grid, self.down_fractions)` line after the two assignments in `MultiFrequencyMonitor.__init__`, `RunningTauMonitor.__init__` and `AwakeRunningTauMonitor.__init__` (with the import added in each module; the exact hunks are in `proposed_fixes.patch`).

### m4. The known-baseline mode ignores `c_down` without an error

*Location.* `src/avmonitor/periodic.py:385` to `386` and `:413` to `421`. *What is wrong.* `PeriodicVolatilityMonitor(baseline="known", tau=3, m0=..., c_down=1.5)` is accepted; the decrease side uses `tau_down` (default τ), and `segments` records `c_down=1.5` next to `tau2_dn=9.0`, so the record says something that did not happen. *Reproduction* (`r9_small_checks.py`, item 8): `segment says c_down=1.5, but tau2_dn=9.0`. *Fix:*

```diff
--- a/src/avmonitor/periodic.py
+++ b/src/avmonitor/periodic.py
@@ -382,8 +382,9 @@ class PeriodicVolatilityMonitor(VolatilityMonitor):
         if baseline == "known":
             if tau is None or m0 is None:
                 raise ValueError("baseline='known' needs tau and m0")
-            if kw.get("margin_z") or kw.get("target_rel_width") is not None:
-                raise ValueError("margin_z and target_rel_width do not apply to a known baseline")
+            if kw.get("margin_z") or kw.get("target_rel_width") is not None or kw.get("c_down") is not None:
+                raise ValueError("margin_z, target_rel_width and c_down do not apply to a known baseline "
+                                 "(give tau_down and m0_down for the decrease side)")
             td = float(tau) if tau_down is None else float(tau_down)
             md = float(m0) if m0_down is None else float(m0_down)
             for name, (tt, mm) in (("increase", (float(tau), float(m0))), ("decrease", (td, md))):
```

### m5. The Bitstamp loader passes the provider's filled minutes on as zero returns

*Location.* `src/avmonitor/data.py:58` to `59` (`load_bitstamp_csv`), used by `cli.py:68` to `69` (`run --source bitstamp`) and by E5c, E6b and E8. *What is wrong.* The module docstring says that the Bitstamp file fills missing minutes with flat zero-volume candles (`data.py:7`). The loader flags them (`zero_volume`) but computes `logret` from the filled closes, so each filled minute becomes an exact zero return, and the return after it spans two minutes. `docs/math.md` 5.6 item 5 says missing returns must never be replaced by zeros; `load_binance` sets the returns on both sides of a gap to NaN. *Size* [REAL] (`r9_bitstamp.py`, which needs the raw file): from 2025-07-01 to 2026-01-01 the file has 264,960 minutes, of which 782 (0.30%) are zero-volume flat candles, all with a zero return, in 702 runs (longest 43 minutes); there are 15,180 zero returns in all, so most zeros are genuine unchanged prices. The effect on E5c, E6b and E8 is therefore small, but the CLI and those experiments do not follow the stated rule. Whether to rerun them is for the owner of those experiments. *Fix:*

```diff
--- a/src/avmonitor/data.py
+++ b/src/avmonitor/data.py
@@ -55,8 +55,10 @@ def load_bitstamp_csv(path: str | Path, start: str | None = None, end: str | Non
         df = df[df.index >= pd.Timestamp(start, tz="UTC")]
     if end is not None:
         df = df[df.index < pd.Timestamp(end, tz="UTC")]
-    df["logret"] = np.log(df["close"]).diff()
     df["zero_volume"] = df["volume"] <= 0
+    # A filled (tradeless) minute repeats the previous close: no price was observed, so the returns on both
+    # sides of it are missing, not zero (docs/math.md 5.6, item 5), as for the gaps of load_binance.
+    df["logret"] = np.log(df["close"].where(~df["zero_volume"])).diff()
     return df
 
 
```

### m6. `two_sided_alarm` with alarm arrays of different integer types

*Location.* `src/avmonitor/detectors.py:266` to `273`. *What is wrong.* "No alarm" is mapped to the largest value of each array's own type and mapped back with the first array's type only, so an int64 and an int32 array give 2,147,483,647 instead of -1. All package detectors create int64 arrays, so this needs user-built `DetectorResult` objects. *Reproduction* (`r9_small_checks.py`, item 7): `two_sided_alarm(int64 [-1, 5], int32 [-1, 7]) = [2147483647 5]`. *Fix:*

```diff
--- a/src/avmonitor/detectors.py
+++ b/src/avmonitor/detectors.py
@@ -265,12 +281,10 @@ def betting_edetector(
 
 def two_sided_alarm(up: DetectorResult, down: DetectorResult) -> np.ndarray:
     """First alarm of either side. Run each side at threshold ``2 A`` to keep ARL >= A overall."""
-    a, b = up.alarm.copy(), down.alarm.copy()
-    a[a < 0] = np.iinfo(a.dtype).max
-    b[b < 0] = np.iinfo(b.dtype).max
-    out = np.minimum(a, b)
-    out[out == np.iinfo(a.dtype).max] = -1
-    return out
+    a, b = np.asarray(up.alarm, np.int64), np.asarray(down.alarm, np.int64)
+    big = np.iinfo(np.int64).max
+    out = np.minimum(np.where(a < 0, big, a), np.where(b < 0, big, b))
+    return np.where(out == big, -1, out)
 
 
 # ----------------------------------------------------------------------------
```

### m7. `design.sr_steps_constant_factor` is one step too large at exact boundaries

*Location.* `src/avmonitor/design.py:403` to `406`. *What is wrong.* When `A` equals `M_n = F (F^n - 1) / (F - 1)` exactly, the logarithm can round up and `ceil` returns n + 1; the guard only corrects downward rounding ("guard against rounding at an exact integer"). In a grid of 195 exact boundaries (F in {1.01, 1.125, 1.5, 10/9, 2}, n = 1 to 39) the function was one step too large in 19 (`r9_small_checks.py`, item 10). `decrease_speed_limit` uses it, but at the settings reported (A = 300, 7800, 15600, 28800, 57600 with the Gaussian B = 9.045) it agrees with a brute-force recursion (31, 58, 64, 69, 75 steps), so no reported number changes. *Fix:*

```diff
--- a/src/avmonitor/design.py
+++ b/src/avmonitor/design.py
@@ -401,8 +401,10 @@ def sr_steps_constant_factor(A: float, F: float) -> int:
     if not F > 1:
         raise ValueError("F must exceed 1")
     n = int(np.ceil(np.log1p(A * (F - 1.0) / F) / np.log(F)))
-    while F * (F**n - 1.0) / (F - 1.0) < A:  # guard against rounding at an exact integer
+    while F * (F**n - 1.0) / (F - 1.0) < A:  # guard against rounding at an exact integer, both ways
         n += 1
+    while n > 1 and F * (F ** (n - 1) - 1.0) / (F - 1.0) >= A:
+        n -= 1
     return max(n, 1)
 
 
```

## 4. The proposed fixes, applied to a copy

`proposed_fixes.patch` holds every diff above (7 files, 18 hunks) and applies to `949c1d5` with `git apply` (checked with `git apply --check`). `round9_regression_tests.py` holds 10 regression tests (one or more per finding). I applied the patch to a copy of the package in the scratchpad, never to `src/`, and ran:

* the whole committed test suite at `9c5786d` (`python -m pytest -q` from the repository root, with the patched `src/` first on the path): **505 passed**, the same 505 that pass on the unpatched package; the 10 new tests against the patched copy: **10 passed**;
* the new tests against the unpatched package: **all 10 fail**, so each one detects its finding;
* `r9_nan_batch.py` against the patched copy: the all-NaN stream and the stream of 2A - 1 NaN then Y = 1 give no alarm, `betting_edetector` and `sr_log_paths` agree (M₆ = 1.0 for both), and on the session-market grid the batch detector alarms on 32.3% of the paths, the same 32.3% as the streaming monitor; the degenerate window now raises in `calibrate` (m2).

The fixes change no output on complete data with valid inputs: the existing suite, which includes recorded reference outputs of the monitor (`tests/test_monitor_cs.py::test_default_output_matches_recorded_reference`), passes unchanged.

## 5. What I checked and found correct

Each item names the code, what I checked and the result. All Monte Carlo results are [SYNTHETIC], with standard errors sqrt(p (1 - p) / paths).

1. **Predictability of every bet.** The adaptive bet of `betting_edetector` at step t uses Y over steps t - h - 1 to t - 1 (`detectors.py:149` to `151`); the confidence-sequence bets at observation i use the centre and variance after i - 1 observations, in both the vectorised (`confseq.py:158` to `166`) and the streaming version (`:214` to `226`). Changing Y_t (or Z_t) changed no bet at or before t: largest change 0 for both (`r9_small_checks.py`, item 2). The fixed grids of the streaming monitors are set when calibration ends; the running truncation level of `RunningTauMonitor` and `AwakeRunningTauMonitor` is computed from returns before the step (the window is updated after the bet), and the awake band's snapshot for a chain started at s uses the first `e(s) <= s` returns.
2. **Range [0, B] and the caps.** At r = 0 (Y = 0), r = ±τ (Y = B), r = 10⁶ and r = τ (1 + 10⁻¹⁵), the smallest increase factor was 0.5 (= 1 - 0.5) and the smallest decrease factor 0.050000000000000 (= 1 - 0.95), so every factor is positive at both ends (item 4). The decrease side is switched off when there is no admissible bet (B ≤ 1, or a lower bound of 0) and its statistic then stays at 0 (`monitor.py:319` to `324`). `betting_edetector` refuses x above the stated upper bound for the decrease side (`detectors.py:215` to `216`).
3. **Threshold used against threshold documented.** Every streaming class raises an alarm when a side's average reaches 2A (`monitor.py:357` to `360` and `:406` to `409`, `periodic.py:461` to `464`, `multifreq.py:531` to `534`, `running.py:733` to `736`, `:826` to `830`, `:1141` to `1144`), as the docstrings and `docs/math.md` 5.4 say; the CLI sets A = steps per day / budget (`cli.py:28` to `30`); `restarted_sr_bound` is min(1, (m + carry) / A), which matches Theorem 1 of `docs/theory_per_period.md` with the carry applied to each side's average (`periodic.py:427` to `434`).
4. **Recursions and weights.** Shiryaev-Roberts `R_t = f_t (R_{t-1} + 1)` and CUSUM `C_t = f_t max(C_{t-1}, 1)` in every class; the mixture weights are 1/K, and with the adaptive bettor (1 - w)/K and w, which sum to 1 (`detectors.py:226` to `229`); multi-frequency weights are normalised in `resolve_weights`.
5. **Off by one in "the first m steps".** Under the calm law Y = 1 (null on both sides) the known-baseline monitor with A = 50 alarmed at the 0-based indices 99, 199, 299, 399, that is at the 2A-th monitored step after each restart, as Proposition 1(b) of `docs/theory_per_period.md` says, and never with a period of 99 < 2A (item 5). The known mode starts at the first return (`segments[0]["start"] == 0`) and restarts at the step after an alarm.
6. **No state leaks across a restart.** `_reset` clears the buffer, both statistics, the extension and fallback flags, both confidence sequences and the split counter (`monitor.py:172` to `182`, `:202` to `207`); every quantity used for betting (τ², τ²_dn, m_up, m_dn, the decrease bets, `down_active`) is reassigned in `_finish` before the first monitored step. Comparison A below (about 10,000 restarts per law, every one compared alarm by alarm with a batch detector started afresh) would show any leak.
7. **Batch functions against streaming classes** (`r9_batch_vs_stream.py`, 150 paths per law with a rise and a fall of volatility, A = 60): (A) known baseline with a restart after every alarm, `PeriodicVolatilityMonitor` against repeated `betting_edetector` + `two_sided_alarm`: 0 paths differ, with 11,318 alarms (Gaussian) and 9,966 (t3); (B) plug-in, `VolatilityMonitor` against `calibrate` + `betting_edetector`: 0 of 150 differ, both laws; (C) confidence-sequence baseline, against `calibrate(baseline="cs")`: 0 of 150 differ, both laws; (D) 10% missing returns: 150 of 150 differ (M1), 0 with the fix. Also: `PanelMonitor` against `panel.run_block`, 20 panels with changes and 5% missing returns: 0 differ (45 discoveries) (`r9_panel.py`); `MultiFrequencyMonitor(freqs=(1,))` against `VolatilityMonitor`: 0 of 30 differ (160 alarms); the multi-frequency vectorised building blocks against the streaming class with 5% missing returns: 0 of 30 differ (`r9_multifreq.py`).
8. **Confidence sequences.** `eb_confseq` and `StreamingEBConfSeq` agree to 1.8e-15 (upper) and 4.4e-16 (lower) over 2000 observations; the monitor's confidence-sequence segment gives exactly the τ², m_up and m_dn of `calibrate(baseline="cs")` (differences 0) (item 3). The one-sided ranges are 1 - ẑ and ẑ with ẑ strictly inside (0, 1), so Lemma 6.1 applies; δ is split as δ/2 per bound (`monitor.py:279` to `284`, `detectors.py:89`), with `c_down != c` the upper bound comes from one sequence and the lower from the other, each at δ/2 (`:310` to `312`); the panel uses δ/K per asset; a started sequence is never restarted (`monitor.py:299` to `304`).
9. **e-BH** (`panel.ebh`): against a brute-force implementation of Wang and Ramdas (2022) on 20,000 random cases, 0 mismatches among 6,667 cases with continuous values and 6,666 with ties, zeros and infinite values; the 1,241 mismatches among 6,667 cases arise only when every value is computed as K / (α r), exactly on an e-BH boundary, where the rounding of e r against K / α decides (item 1). The rule rejects the R largest values with R = max{r : e_(r) ≥ K / (α r)}, as documented.
10. **Multi-frequency certified value.** Between block ends the certified value of a coarse frequency was never above its full statistic at the next block end: largest difference 0 over 40,704 minutes of 40 paths with 80 alarms and 5% missing returns (`r9_multifreq.py`, item 1). This is the property Theorem 1 of `docs/theory_multifreq.md` rests on.
11. **Rank bets of route B.** Each of the 16 betting functions of `running.bet_family` (`"lr"`: 3 per side, `"ind"`: 5 per side) integrates to 1 over a uniform p-value, within 6.3e-10 by numerical quadrature, and is at least 0.05 on [0, 1], so every factor of `ConformalMonitor` is a valid e-value under exchangeability (`r9_running.py`, item 1).
12. **`limits.py` spot checks** (`r9_limits.py`). `sr_alarm_times` gives the same first alarm as `betting_edetector` (same bets, equal weights, same threshold) on 200 paths per side and two thresholds: 0 differing paths, 800 alarms compared. `collapse_floor` equals `floor_time` on the constant extreme path and a direct search for the smallest n with n log F ≥ log(A / n) in 32 of 32 cases (A in {50, 300, 7800, 28800}, B in {2, 5, 9.045, 20}, both sides). The two-point law that `max_null_mass` names has mean 1 to within 8.9e-16.

**Monte Carlo under the null, through the streaming classes as a user calls them** (fresh monitor per path, one `update` per return; `r9_mc_null.py` unless stated). Known baseline: τ = 3 for Gaussian returns (m₀ = 0.9950) and τ = 3 times the population scaled MAD for t3 (τ = 3.4021, m₀ = 1.8191 by numerical integration; a Monte Carlo check with 4 million draws gave 1.8165).

| Stream (returns i.i.d. unless stated) | Monitor | m | A | Paths | False alarm rate (se) | Bound | Rate / bound |
|---|---|---|---|---|---|---|---|
| Gaussian, first m steps | `PeriodicVolatilityMonitor(baseline="known")`, no period | 390 | 7800 | 40,000 | 0.88% (0.05%) | m / A = 5% | 0.18 |
| Student t3, first m steps | `PeriodicVolatilityMonitor(baseline="known")`, no period | 390 | 7800 | 40,000 | 1.42% (0.06%) | m / A = 5% | 0.28 |
| Gaussian, first m steps | `PeriodicVolatilityMonitor(baseline="known")`, no period | 1440 | 28800 | 20,000 | 1.92% (0.10%) | m / A = 5% | 0.38 |
| Student t3, first m steps | `PeriodicVolatilityMonitor(baseline="known")`, no period | 1440 | 28800 | 20,000 | 2.58% (0.11%) | m / A = 5% | 0.52 |
| Bold play (E16), default grid, first m steps | `PeriodicVolatilityMonitor(baseline="known")`, no period | 390 | 7800 | 20,000 | 2.28% (0.11%) | m / A = 5% | 0.46 |
| Bold play (E16), one bet per side (0.5 up, 0.95 down), first m steps | `PeriodicVolatilityMonitor(baseline="known")`, no period | 390 | 7800 | 20,000 | 4.34% (0.14%) | m / A = 5% | 0.87 |
| Gaussian, alarm on day d = 1 to 5 | `PeriodicVolatilityMonitor(baseline="known", period=390)` | 390 | 7800 | 10,000 | 0.96%, 0.91%, 0.76%, 1.10%, 0.98% (se at most 0.10%) | m / A = 5% each day | 0.22 (largest) |
| Bold play, one bet per side, alarm on day d = 1 to 5 | `PeriodicVolatilityMonitor(baseline="known", period=390)` | 390 | 7800 | 10,000 | 4.73%, 4.28%, 4.93%, 4.67%, 4.42% (se at most 0.22%) | m / A = 5% each day | 0.99 (largest) |
| Calm (Y = 1) for 9 days, Gaussian day 10: alarm on day 10 | `PeriodicVolatilityMonitor(baseline="known")`, no period | 390 | 7800 | 4,000 | 27.75% (0.71%) | none; `docs/theory_per_period.md` reports 28.0% for this stream | n/a |
| The same stream: alarm on day 10 | `PeriodicVolatilityMonitor(baseline="known", period=390)` | 390 | 7800 | 4,000 | 0.68% (0.13%) | m / A = 5% | 0.13 |
| Gaussian, first m steps after calibration | `VolatilityMonitor(baseline="cs", cal_len=1440, delta=0.01)` | 1440 | 28800 | 4,000 | 0 of 4,000 | δ + m / A = 6% | 0 |
| Student t3, first m steps after calibration | `VolatilityMonitor(baseline="cs", cal_len=1440, delta=0.01)` | 1440 | 28800 | 4,000 | 0.05% (0.04%) | δ + m / A = 6% | 0.01 |
| Gaussian, first m steps | `MultiFrequencyMonitor` (1, 5, 15), known baseline (`r9_multifreq.py`) | 390 | 1114 | 4,000 | 0.30% (0.09%) | `proved_bound(390)` = 5% (A = `level_for_budget(0.05, 390)`) | 0.06 |
| The same, vectorised building blocks | `multifreq` building blocks | 390 | 1114 | 40,000 | 0.39% (0.03%) | 5% | 0.08 |
| Student t3, first m steps after 200 calibration returns | `ConformalMonitor(cal_len=200)` (`r9_running.py`) | 390 | 7800 | 2,000 | 2.75% (0.37%) | m / A = 5% | 0.55 |
| Gaussian, first m steps after calibration | `RunningTauMonitor(cal_len=1440)` (`r9_running.py`) | 1440 | 28800 | 1,000 | 0 of 1,000 | δ + m / A = 6% | 0 |
| Gaussian, first m steps after 120 pre-sample returns | `AwakeRunningTauMonitor` (`r9_running.py`) | 1440 | 28800 | 1,000 | 0 of 1,000 | δ + m / A = 6% | 0 |
| Gaussian, K = 5 assets, any discovery in one block | `panel.run_block`, e-BH, known baselines (`r9_panel.py`) | 390 | n/a (α = 5%) | 20,000 | 0.71% (0.06%) | α = 5% | 0.14 |
| Student t3, K = 5 assets, any discovery in one block | `panel.run_block`, e-BH, known baselines (`r9_panel.py`) | 390 | n/a (α = 5%) | 20,000 | 1.42% (0.08%) | α = 5% | 0.28 |
| Gaussian, 10% of returns missing, K = 5 assets, any discovery in one block | `panel.run_block`, e-BH, known baselines (`r9_panel.py`) | 390 | n/a (α = 5%) | 20,000 | 0.61% (0.06%) | α = 5% | 0.12 |

Reading the table:

* Every rate with a bound is below it. The i.i.d. rates are far below. Two rates can be set against E16 (same detector, m = 390, A = 7800): Gaussian returns, 0.88% (0.05%) here against 0.97% (0.93% to 1.02%) there; bold play with the default grid, 2.28% (0.11%) here against 2.40% (2.34% to 2.47%) there.
* Bold play with one bet per side is the only case near its bound: 4.34% on the first day, and with daily restarts 4.28% to 4.93% per day over five days (mean of the five 4.61%), against 5% each day. Theorem 1 of `docs/theory_per_period.md` bounds each day by m / A given the past, and no day is above it; the largest, 4.93% (se 0.22%), is within one standard error of the bound, so this adversary nearly attains it. These are single bets, not the default grid, so E16's proved limit of 0.98 of the bound at one day (which is for the default grid) does not apply to them.
* The calm-then-Gaussian row without a period is not a bound that fails: `docs/theory_per_period.md` (section 1, item 1) gives this stream as the example where a single day exceeds m / A, with 28.0%; this run gives 27.75% (0.71%). With `period=390` the same stream gives 0.68%, below m / A, as Theorem 1 says.
* The panel rates use the vectorised `run_block`, which `PanelMonitor` matches alarm by alarm (item 7).
* The confidence-sequence, `RunningTauMonitor` and `AwakeRunningTauMonitor` rates are near 0. No alarm in 4,000 paths puts the rate below about 0.075%, and none in 1,000 below about 0.3% (one-sided 95%, rule of three).
