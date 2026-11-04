# W3. The look-ahead in E7 (Nifty) at `experiments/e07_nifty.py` line 150: what it is and what it changes

Script: `reports/round4/W3/w3_e07_lookahead.py` (log: `w3_e07_lookahead.log`; alarm and summary tables per variant: `w3_alarms_*.csv`, `w3_summary_*.csv`; side by side: `w3_primary_comparison.csv`). The three `w3_alarms_*.csv` tables (about 1.2 MB each) are not committed; the script rebuilds them in a few minutes. Nothing under `results/`, `experiments/`, `src/` or `docs/` was changed. The rerun of the committed pipeline reproduces `results/e07b_alarms.csv` (every alarm time of all 8 configurations) and `results/e07b_monitor_summary.csv` exactly, so the comparison below is against the real committed numbers.

## Summary

* The checker is right. Line 150 gives the first 10 Nifty sessions (3 to 18 April 2017, 3,740 returns, 0.4% of the sample) an intraday profile computed from the first 20 sessions, so those 10 sessions are scaled with data from their own future (up to 19 sessions ahead). Every later session uses only past sessions.
* Removing the look-ahead (leave those 10 sessions unscaled, so the monitor skips them, as it does with any missing return) changes the primary monitor's alarms only in April 2017: 3 alarms disappear (6, 12 and 18 April, all "up") and 2 move by 1 and 60 minutes (21 and 26 April, both "down"). From 26 April 2017 10:08 onwards the 606 remaining primary alarms are identical in time and direction. Across all six robust configurations the alarm lists are identical from 2 June 2017 at the latest.
* Of the headline numbers, only the development-period alarm rate of the primary row flips at printed precision: 0.26 becomes 0.25 per session (300 to 297 alarms over 1,170 sessions). That flip is mostly a denominator effect: with the 10 unmonitored sessions left out it is 297 / 1,160 = 0.26 again. Every held-out number is unchanged to the last digit. The README sentence "9 of 11 Union Budget sessions raise an alarm against 25% of ordinary sessions" is unchanged (9 of 11; ordinary 540 to 538 of 2,178, both 25%).
* Seven other cells of the 16-row table in `docs/experiment_report.md` (all in the development period) move by one unit of the last printed digit, mostly Wilson interval endpoints and the random-time controls, whose random draws change once the pool of monitored minutes changes. No reading or conclusion changes.
* Recommended fix: delete line 150 (one line), which is exactly the package rule `start="past"`, and state in the E7 docs that the monitor starts at the 11th session. Expected cost: 3 fewer alarms in April 2017 and nothing else of note.

## What the look-ahead is

`session_profile` (`experiments/e07_nifty.py`, lines 137 to 156) builds a table `wide` with one row per session and one column per clock minute, holding |r|. Line 149 takes, for each minute, the median over the previous 20 sessions, requiring at least 10 of them (`min_periods=10`) and shifting by one session so the current session is excluded. For sessions 1 to 10 fewer than 10 earlier sessions exist, so the whole row is NaN. Line 150, `prof = prof.fillna(wide.iloc[:lookback].median())`, replaces every NaN cell by the median of |r| over sessions 1 to 20. For session k (k = 1 to 10) that median includes session k itself and sessions k+1 to 20, which lie in the future. This is the behaviour that `avmonitor.data.intraday_scale` calls `start="pooled"` and no longer uses by default (its docstring: "uses up to lookback_days days of future data"); the E7 docstring (lines 142 to 144) states it openly, so it is a documented shortcut rather than a hidden one, but it breaks the project rule that the deseasonalisation uses only past data.

Affected returns: the 10 sessions 3, 5, 6, 7, 10, 11, 12, 13, 17 and 18 April 2017 (there is no earlier Nifty minute data, so no warm-up is available). Nothing after 18 April 2017 uses the future. Line 150 also fills a second, harmless group of cells: from 17 August 2026 the last 13 or 14 minutes of each session are mostly frozen (set to NaN by `load_clean`), so their rolling median is NaN too, and line 150 fills them with the April 2017 value. That is stale, but it is past data, and it touches 18 real returns.

How wrong is the look-ahead profile? Compared with the first past-only profile (session 11), the pooled profile is within a factor 0.77 to 1.19 at 90% of the minutes (median ratio 0.97). The primary monitor calibrates on the first 750 returns (two sessions) and raised its first alarm on 6 April 2017 at 15:06, inside the affected window, so the alarms of April 2017 are the ones at risk.

## How big the effect is

Primary configuration (robust e-detector, calibration 750, band 0.25, A = 7,500, c = 3), full run, same seed 7 for the random-time control. "Committed" is line 150 as is. "Fixed" leaves the first 10 sessions unscaled (NaN, skipped by the monitor) and keeps line 150 for the later, past-only cells; deleting line 150 outright gives the same numbers in every cell below (it only shifts 7 alarms of August to October 2026 by 1 to 3 minutes, because 23 closing-minute returns become unmonitored).

| Quantity (where it is quoted) | Committed | Fixed | Printed number changes? |
|---|---|---|---|
| Primary alarms, all | 611 | 608 | n/a |
| Primary alarms in 2017 (April 2017) | 45 (5) | 42 (2) | n/a |
| Dev alarms per session (docs table, README "0.26 per session") | 300 / 1,170 = 0.256 | 297 / 1,170 = 0.254 | yes, 0.26 to 0.25 (0.26 if divided by the 1,160 monitored sessions) |
| Held-out alarms per session | 311 / 1,191 = 0.261 | 0.261 | no |
| Dev up alarms: next session above baseline | 73.3% (165 alarms) | 72.8% (162) | no (73%) |
| Dev down alarms: next session below baseline | 85.9% (135) | 85.9% (135) | no |
| Held-out up / down hit shares | 78.3% / 76.0% | 78.3% / 76.0% | no |
| Dev random-time controls up / down | 45.5% / 54.5% | 45.0% / 55.0% | yes, 46 to 45 and 54 to 55 (random draws change) |
| Held-out controls up / down | 47.0% / 53.0% | same | no |
| Dev / held-out median delay (steps) | 462.5 / 404 | 464 / 404 | no (463 / 404 rounded) |
| Dev next vs previous session, up / down | 72.7% / 66.7% | 72.2% / 66.7% | no (73% / 67%) |
| Budget sessions with an alarm (README, E7 Q7d) | 9 of 11 | 9 of 11 | no |
| Ordinary sessions with an alarm (README "25%") | 540 / 2,178 = 24.8% | 538 / 2,178 = 24.7% | no |

Alarm times of the primary configuration: committed only: 2017-04-06 15:06 up, 04-12 09:55 up, 04-18 14:52 up, 04-21 11:42 down, 04-26 10:08 down; fixed only: 2017-04-21 10:42 down, 04-26 10:07 down. All other 606 alarms are identical. The two late-April alarms move because the fixed monitor starts calibrating on 19 April instead of recalibrating after the 18 April alarm, so its baseline is slightly different for one cycle. Other robust configurations: 2 to 7 alarms lost and 2 to 4 gained, all before 2 June 2017; Gaussian LR: dev alarms per session 0.70 to 0.69 for calibration 375, otherwise unchanged at printed precision.

Of the 16 rows of the docs table, the 8 held-out rows are unchanged in every cell. In the 8 development rows the point estimates move by at most 2 percentage points (robust, calibration 750, band 0.5: up hit 77% to 79% on 70 alarms, one alarm fewer), and 8 cells, mostly interval endpoints, controls and median delays, move by one unit of the last digit (listed at the end of the log).

## Recommended minimal fix (not applied)

1. Delete line 150 of `experiments/e07_nifty.py`. Line 149 then leaves sessions 1 to 10 without a profile, `deseasonalise` turns them into NaN, and the monitor skips them as missing (the same path it uses for frozen minutes). This is the package's `start="past"` rule without a warm-up, since no Nifty data before 3 April 2017 exist.
2. Change the last sentence of the `session_profile` docstring (line 144) to "The first 10 sessions get no profile and are not monitored."
3. In `docs/experiment_report.md` Q7b and `README.md`: say the monitor starts at the 11th session (19 April 2017), replace the primary dev row by 297 alarms, 0.25 per session (or 0.26 over 1,160 monitored sessions), controls 45% and 55%, delay 464, and refresh the other 7 development cells from `w3_summary_no_lookahead.csv`. Then rerun e07 q7b, e07c and e07d to refresh `results/`.

If the late-2026 closing minutes should keep a scale (they lose it once line 150 goes), an alternative is to fill NaN cells with the profile of the previous session only where a past profile exists; this changes no number in the table above.

## What this does not show

* Only Q7b and the Budget-day share of Q7d were recomputed. The India VIX checks of Q7c use the same alarm lists; the dev up-alarm session count would drop from 165 to 162, so its dev percentages may move by about one point. The VIX files are not in this container (`data/raw/india_vix/` is git-ignored and absent), so this was not measured.
* E13, E15e and E18 reuse the Nifty input; their 2017 to 2021 development numbers could shift by the same April 2017 alarms. Nothing in 2022 to 2026 is affected, because by then both variants have the same input and the same monitor state.
* The 18 to 23 stale-filled closing minutes from August 2026 are a separate small untidiness, not a look-ahead; the table shows they do not change any printed number.
