# Round 4: backtests, stress tests and repairs (4 November 2026)

Start here for round 4. It ran from 01:48 to about 03:45 UTC on 4 November 2026, on top of `main` (commit 4de8a32, rounds 1 to 3 and the validation audit) and PR #4 (commit cd3dd4f, the audit fixes made by another thread). All the work was done by AI agents (Claude Opus and Claude Fable models), not by people: eight agents and a coordinator. Every number below comes from a report in this folder; open that report for the tables, seeds and commands.

Labels as in the rest of the repository: **[SYNTHETIC]** = our simulation, **[REAL]** = our run on market data. VERIFIED / PARTIAL / FAILED / UNVERIFIED say how far a claim was checked.

## The short version, in plain words

1. **The code is in better shape and no result moved.** PR #4 fixed most of the audit's open items. This round adds what it lacked: the alarm level must be a finite number above 1, the truncation constants must be finite and positive, a new `status()` shows when the monitor is blind, a third benchmark detector skips missing values, loaders also blank out infinite prices, and a crash on absurdly large returns is gone. The quick reproduction tier and direct old-against-new checks on real data gave identical outputs (465 of 465 result files; 11 of 11 loader frames; 26 of 26 monitor runs). [VERIFIED] (`M/REPORT.md`)
2. **The data are clean and the pre-registered tests reproduce.** All 139 fixed data files match their committed hashes byte for byte, and the daily Bitstamp file still contains the committed rows exactly. No duplicate timestamps, bad prices or bad prints in 170 files. E22 and E23, the two pre-registered real-data tests of round 3, reproduce exactly on a clean copy, and so do E25's gates and E19's pools (463 files identical, 2 equal except a timing column, 0 different). E18 was stopped at its time budget after writing one file, which is identical (partial). [REAL, VERIFIED] (`R1/REPORT.md`)
3. **At an equal false-alarm rate, the e-detector's statistic is as fast as the classical rules on normal returns and much faster on heavy tails.** When every method is tuned to the same 5% of false alarms per day, the e-detectors need 37 to 44 steps to catch a doubling of variance and the classical rules 38 to 43. On Student t returns the classical rules need 203 to 750 steps against 74 to 91. Tuned on normal data and then used on heavy tails, the classical rules raise a false alarm on 99.9% of days. [SYNTHETIC, VERIFIED, and checked independently by V1] (`R2/REPORT.md`)
4. **The README setting is slow because its threshold is far too safe, not because of its statistic.** With the null band of 0.25 and A = 7,800 it raised 0 false alarms in 20,000 simulated calm days, so it uses almost none of its 5% budget. The same statistic at a threshold that uses the budget is 3.1 to 6.4 times faster. [SYNTHETIC, VERIFIED] (`R2/REPORT.md`) A post hoc follow-up on permuted real calm days (E19's method, 40,000 streams on five coins, 8,000 on the Nifty 50) found that a band of 0.1 keeps the promise with room to spare (0.27% false alarms within a day against 5%, 1.41% over five days against 25%) and detects 1.6 to 2.4 times faster in simulation; the proved confidence-sequence baseline with band 0 is about as fast as today's README setting. Band 0 breaks the five-day promise (41.9%), as E19 found. The default was not changed. [REAL marginal and SYNTHETIC, VERIFIED, post hoc] (`R6/REPORT.md`)
5. **On real data from 2017 to October 2026, the alarms behave the same way in every asset and year, including fresh 2026 data, but they mostly say "volatility is already high".** After an up alarm the next hour is above the calibrated level in 81% to 98% of cases per asset and year, against 25% to 38% at random times. But compared with random minutes that were just as volatile, the gain mostly disappears, and the hour after an up alarm is calmer than the hour before it (partly by construction). So an alarm means "volatility is now well above its calibrated level", not "volatility will keep rising". W1 then matched the random minutes on last-hour volatility, which answers V1's caveat: up alarms keep a real edge (91.4% against 74.6% at equally volatile minutes, 24 of 27 asset-years pass), but down alarms add almost nothing (95.2% against 92.9%). [REAL, PARTIAL, checked independently by V1; matched placebo VERIFIED for up alarms, FAILED for down alarms], checked independently by V1] (`R4/REPORT.md`, `V1/REPORT.md`, `W1/REPORT.md`)
6. **"Blind about 70% of the time" is right for crypto, not for the Nifty 50.** Crypto: 64% to 72% in all 17 asset-years. Nifty 50: 43% to 65% by year, about half overall. The README now says so. [REAL, VERIFIED] (`R4/REPORT.md`)
7. **The confidence-sequence baseline's blind mode is not a Binance quirk.** From mid August 2026 it could not raise an increase alarm for 49 days on Binance BTC, and for 29 days on Bitstamp BTC at the same time. [REAL, VERIFIED] (`R4/REPORT.md`)
8. **The monitor watches the clipped squared return, not the variance, and three realistic situations bend the guarantee.** A switch from thin to fat tails at constant variance sets the known-baseline detector off on 65% of paths within a day; a variance increase carried only by rare large jumps can go unseen. Removing a U-shaped intraday pattern with the past-only estimator the real-data experiments use keeps the one-day bound (3.4%) but breaks the five-day one for the known baseline (32.2% against 25%); strong bid-ask bounce breaks it too (40.7%). The README setting and the confidence-sequence baseline stay far inside in all of these. Skewed returns, missing values, tick grids with up to 60% zero returns and misspecified GARCH forecasts all stayed within the bounds. A doubling during calibration leaves the README setting silent for five days on 56% of paths, and a 30-minute burst is caught on only 27%. These are new limitations, now in `docs/limitations.md`. [SYNTHETIC, VERIFIED] (`R3/REPORT.md`)
9. **The property tests found two small real bugs, both now fixed.** 56 new test cases check rules that must hold for every input (scaling, sign, streaming equals batch, missing values, pickling, threshold order, and the supermartingale property the proof rests on). All the rules the proof needs held for every law tried. The two bugs (a crash on returns above about 1e153, and an infinite threshold that could still alarm) are fixed by M. [SYNTHETIC, VERIFIED] (`R5/REPORT.md`)
10. **The look-ahead in E7 is real but tiny.** `experiments/e07_nifty.py` line 150 scales the first 10 Nifty sessions (3 to 18 April 2017) with a profile that partly uses their future. Removing it removes 3 up alarms in April 2017; the other 606 alarms are identical, every held-out number and the Budget result (9 of 11 sessions) are unchanged, and the development alarm rate goes from 0.26 to 0.25 per session. The one-line fix is not applied, so no committed number changed. [REAL, VERIFIED] (`W3/REPORT.md`)
11. **On fresh 2026 data the alarms do not beat a strong forecast.** From 1 January to 24 October 2026, adding the alarms to E20's HAR forecast (which already uses the last hour) did not improve the next-hour forecast for BTC, ETH or the Nifty 50 (loss change -0.54% to +0.05%, all intervals covering 0), and after an up alarm the next day on the five coins was about 10% calmer than that forecast. Against a forecast that does not see the last hour, the alarms do add information. So the alarm is a good summary of a burst that is happening, not a forecast. [REAL, PARTIAL by the rule declared in advance] (`W1/REPORT.md`)

## Who did what

| ID | Model | Task | Status | Report |
| --- | --- | --- | --- | --- |
| M | Opus | Repairs to the package on top of PR #4, with proof that no result changed | done | `M/REPORT.md` |
| R1 | Fable | Integrity of every data file; clean-copy rerun of the round 3 real-data experiments | done: E22, E23, E25 gates and E19 pools reproduce exactly; E18 stopped at its time budget (partial: one file written, identical) | `R1/REPORT.md` |
| R2 | Opus | Detection speed of every method at an equal measured false-alarm rate | done | `R2/REPORT.md` |
| R3 | Fable | New kinds of simulated "no change" data, and realistic change shapes | done | `R3/REPORT.md` |
| R4 | Opus | Walk-forward backtest on every asset and year, 2017 to October 2026 | done | `R4/REPORT.md` |
| R5 | Fable | Metamorphic and property tests of the code | done; moved to `tests/test_round4_properties.py` | `R5/REPORT.md` |
| V1 | Opus | Independent check of R2's and R4's code and numbers | done: both CONFIRMED | `V1/REPORT.md` |
| R6 | Fable | Is there a faster setting that still keeps its promise on real calm days? (post hoc, after R2) | done (its pools and batch files are not committed; `r6_pools.py` and `run_null.sh` rebuild them in about 12 minutes) | `R6/REPORT.md` |
| W1 | Opus | Placebo matched on recent volatility, and a strong-forecast check on 2026 data (second wave, plan written first) | done; rerun by the coordinator, identical | `W1/REPORT.md` |
| W3 | Fable | Measure the start-up look-ahead in `experiments/e07_nifty.py` (second wave) | done; fix proposed, not applied | `W3/REPORT.md` |

## What changed in the repository

* Code (`src/avmonitor/`): see `M/REPORT.md`, one row per fix. No committed result changed.
* Tests: 762 in total, up from 656 after PR #4 (50 in `tests/test_round4_mending.py`, 56 in `tests/test_round4_properties.py`).
* README: status line; the blind-share sentence now separates crypto and the Nifty 50; the E21 bullet says that the README setting and the confidence-sequence baseline leave almost all of the 5% budget unused; the usage example says the README setting is very conservative and slow; a new bullet points to this summary.
* `docs/limitations.md`: the confidence-sequence blind mode on Bitstamp; how conservative the README setting is; new limitations from R3 (what the monitor watches, the realistic seasonality pipeline, bid-ask bounce, a change during calibration, short bursts); what a real alarm shows (R4 and V1); two data vendor changes (R1); the test count.
* `CHANGELOG.md`: a round 4 entry.
* In `results/` only the README and pitch number traces were regenerated (0 untraced) and their two lines in `results/MANIFEST.csv` rebuilt; nothing in `experiments/` or the pre-registrations was changed, and neither forward test (E24 after 1 June 2027, E25 from 26 November 2026) was run.

## Open items, in order

1. Decide with Prof. Bapat whether the README should recommend a faster default: the plug-in with band 0.1 (approximate, about twice as fast) or the confidence-sequence baseline with band 0 (proved, as fast as today). Both were found post hoc (R2, R6), so test the chosen one on data not yet seen, for example as a secondary arm of the forward tests, without changing their pre-registered primary settings.
2. Done in this round (W1): the matched placebo and the 2026 forecast check. What is left is wording: the pitch and the paper should describe an alarm as "a burst is happening now", not as a forecast, and should not claim that decrease alarms add information beyond the last hour.
3. The late-2026 Nifty vendor change is now in `docs/limitations.md` (the last 14 minutes of each session are frozen in the 2026-08 to 2026-10 files, so about 360 returns per session are monitored, not 375; `R1/REPORT.md`); any future live run should state 360 in its configuration.
4. Rerun E18 and E19 p1 on a clean copy (each needs more than half an hour on two cores; E18 was stopped at its time budget in this round).
5. The small look-ahead in `experiments/e07_nifty.py` is now measured (W3): delete line 150 and rerun e07 q7b, e07c and e07d when a change of 0.26 to 0.25 in the development alarm rate is acceptable; E13, E15e and E18 reuse the same Nifty input and would need the same rerun.
6. Run the forward tests exactly as pre-registered: E25 from 26 November 2026 (`python experiments/e25_sequential_events.py status`), E24 after 1 June 2027.

## How to check this round

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e . --no-deps
python -m pytest -q -p no:cacheprovider          # expect 762 passed
python scripts/reproduce.py --steps compare      # expect 465 identical
```

Each report's section 2 lists its own commands and seeds. Seeds from 9,600,000 to 9,979,999 were used only in this round.
