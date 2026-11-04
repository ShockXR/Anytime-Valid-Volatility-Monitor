# Pre-registered live run: Nifty 50, session of Tuesday 3 November 2026

*Written and committed before the session opened (09:15 IST, 03:45 UTC, 3 November 2026). The commit time on GitHub is the record of when the settings were fixed.*

**What this is.** A demonstration that the monitor runs on a session it has never seen, with every setting and every reported quantity fixed in advance. It is **not evidence**: one session cannot confirm or refute anything about the method. The evidence is in E7 (2,361 sessions, with a development period and a held-out period).

## Frozen settings

All settings are E7's pre-registered primary configuration (`docs/report_parts/E7.md`), unchanged:

| Item | Setting |
| --- | --- |
| Data | Nifty 50 one-minute bars from the public mirror used in E7 (`avmonitor.data.load_nifty`), cleaned as in `experiments/e07_nifty.py` (`load_clean`: regular sessions only; runs of 5 or more exact zero returns set to missing) |
| Deseasonalising | E7's predictable session profile: median of \|r\| at each minute over the previous 20 sessions, 5-minute centred smoothing, current session excluded; profile data from 1 August 2026 |
| Monitor | `VolatilityMonitor`, A = 7500 (375 / 0.05; each side alarms at 2A), c = 3, calibration 750 returns, null band 0.25, plug-in baseline, recalibration after each alarm |
| Stream | Every regular session from 1 October 2026 up to and including 3 November 2026, in order, so the state of the monitor at the open is set by earlier sessions only |
| Code | `experiments/live_nifty_preregistered.py` as committed with this file |
| Settings hash | `python experiments/live_nifty_preregistered.py config` prints the settings and their SHA-256: `a25d76147b4807a7c66e2b93caaefb1be64485dd05217ae928ff08f4ccb6dfa4` |

## What will be reported, and nothing else

For the session of 3 November 2026, once its bars are in the mirror (usually the same evening), the command

```bash
python experiments/live_nifty_preregistered.py run --day 2026-11-03
```

writes `results/live/2026-11-03/alarms.csv` and `summary.json` with:

1. every alarm in the session: time (IST), direction (up or down) and the baseline scale it was measured against;
2. whether the monitor was monitoring or calibrating at 09:15, and the share of the session's minutes it spent monitoring;
3. the session's robust volatility (scaled median |u| of the deseasonalised returns) divided by the baseline scale in force at the open, or at the first minute the monitor was monitoring;
4. for each alarm, the robust volatility of the rest of the session after the alarm divided by the alarm's baseline; a ratio above 1 after an up alarm, or below 1 after a down alarm, agrees with the alarm;
5. the SHA-256 of every data file used.

No setting will be changed after the session, and no other quantity will be added to the report of this run. If the session has no alarm, that is the result. If the bars for 3 November never appear in the mirror, that will be reported instead.

**How to read it.** Under the null (volatility at its baseline) the chance of any false alarm in a session is at most 5% with a known baseline; with the plug-in baseline used here it is approximately that (E4e, E14). In E7's held-out period the monitor raised 0.26 alarms per session; 78% of increase alarms were followed by a more volatile session against 47% at random times. A single session tells us nothing beyond "the pipeline works on new data".

## Pipeline check on a past session

The same command with `--dry-run` was run on the session of 2 November 2026 (already in the mirror) before this file was committed, to check the pipeline end to end. Result [REAL]: the monitor was monitoring at the open; one down alarm at 11:42 IST against a baseline scale of 0.00030; the rest of the session had 0.84 times that robust volatility (agreeing with the alarm); the whole session had 0.77 times the baseline at the open; the monitor then recalibrated (750 returns, two sessions) and monitored 39% of the session's minutes. Files: `results/live/dryrun_2026-11-02/`. This check is not part of the pre-registered run.

## Result (added after the session, 3 November 2026; run at 11:54 UTC)

*Everything above this heading is unchanged since it was committed before the open.*

The session's bars were in the mirror at about 11:50 UTC. The pre-registered command was run once, with the frozen settings (the settings hash printed by `config` is `a25d76147b4807a7c66e2b93caaefb1be64485dd05217ae928ff08f4ccb6dfa4`, as above). Files: `results/live/2026-11-03/` (`alarms.csv`, `summary.json`). Result [REAL], item by item:

1. **No alarm** in the session.
2. The monitor was **calibrating at 09:15 and for the whole session**: it monitored 0% of the session's minutes.
3. Not defined: no baseline was in force at any minute of the session.
4. Not defined: no alarm.
5. Data files and their SHA-256 are in `summary.json` (the November file, which now holds the sessions of 2 and 3 November: `23c23a628441af156174906545fc3979f284bea4c77ab7e3777a334372ba3f9a`). Its rows for 2 November are identical to those the dry run used.

**Why.** The monitor's down alarm at 11:42 IST on 2 November (the alarm the dry run reported) started a recalibration of 750 returns. The rest of 2 November gave 213 usable returns and 3 November 360, 573 in all, so the monitor was still calibrating at the close of 3 November.

**What we missed.** This was predictable before the settings were frozen: the dry run had already shown the 11:42 alarm and the two-session recalibration that follows every alarm. We did not notice that it would cover the whole target session, and the settings were not changed after the fact. So this run shows only that the pipeline works on a session it had never seen. It also shows, live, a cost reported in E12 and handled in E14: the plug-in monitor is blind while it recalibrates after each alarm (750 returns, two Nifty sessions). E14's awake rank monitor has no such blind period; it was not part of this pre-registration and was not run on this session.

A small inconsistency in the frozen code: the docstring of `experiments/live_nifty_preregistered.py` gives the pre-registration's path as `live/PREREGISTRATION_2026-11-03.md`; the pre-registration is this file.
