# Demo walkthrough (Bapat prototype): about 6 minutes, every command run

Written by an AI agent (Claude), 3 November 2026, validation agent 3. Every command and output excerpt below was run in this session from the repository root with the project's Python environment (`python` below means that interpreter), `OMP_NUM_THREADS=1`. Each excerpt is labelled SYNTHETIC or REAL. The source files `src/avmonitor/detectors.py` and `monitor.py` were being edited by another agent while I ran these, so rerun this page's commands once after the final commit and compare the alarm steps.

## Summary for a beginner (5 lines)

1. The monitor watches a stream of price changes and shouts when the size of the moves (volatility) has changed, for example when a calm market turns wild.
2. "Anytime-valid" means the alarm threshold comes with a proof: if nothing changed, the chance of a false alarm is at most 5% per day, however often you look.
3. It is not a trading bot and does not predict crashes. It says "the volatility level is now different from the baseline it calibrated on".
4. The demo runs on a made-up stream with a known change (so you can see the alarm is right), then on real Bitcoin and Nifty minutes, then replays a famous day minute by minute in a web page.
5. If the internet is down, everything except the two real-data commands still works, and the replay page and notebook already contain the real-data results.

## Before the meeting (3 minutes, once)

| Check | Command | What I saw |
|---|---|---|
| Synthetic demo | `python -m avmonitor simulate --ratio 2` | 0.8 s, output below |
| Real-data demo, needs internet once | `python -m avmonitor run --source binance --symbol BTCUSDT --download --start 2025-10-08 --end 2025-10-14 --plot btc.png` | 2.8 s including the download from the public mirror; data lands in `data/raw/` (git-ignored) so the second run needs no internet |
| Replay page offline | open `docs/demo/replay.html` | works with the network off, no console errors, no external requests, no horizontal scroll at 390 px (checked) |
| Do not run live | `python -m pytest -q` | 581 passed in 104 s on the baseline; too slow for the room |

## The flow

### Step 1 (45 s): what and why

Say: "A volatility monitor that tells you when the market's jumpiness has really changed, with a mathematical guarantee on false alarms, using returns that are clipped so one crazy tick cannot fire it."

### Step 2 (1 min): the synthetic demo with a known change (SYNTHETIC)

```bash
python -m avmonitor simulate --ratio 2
```

Expected output (run today, default seed):

```
[SYNTHETIC] sd x2.0 at step 3000, back to 1 at step 6000, jumps with prob 0.003 per step
alarm level A = 7800: at most 5% chance of a false alarm per 390 steps under the null
robust e-detector: [(3010, 'up'), (6098, 'down'), (8287, 'down')]
Gaussian LR      : [(1000, 'up'), (2564, 'up'), (3581, 'up'), (4779, 'up'), (6036, 'down'), (7570, 'up'), (8823, 'up')]
```

Say: "At step 3000 I doubled the volatility. The robust monitor says 'up' at step 3010, ten steps later. At step 6000 I put it back and it says 'down' at 6098. The classical Gaussian test, which is not robust to jumps, fires seven times, including at step 1000 when nothing had happened."

Be ready for the third robust alarm (8287, 'down'), which looks like a false alarm. It is: the stream has been calm since step 6000. The guarantee is a budget, not zero: 5% per 390-step day. I measured this with `--ratio 1` (no change at all) over seeds 0 to 19, 9000 steps each (23.1 days): 26 robust alarms in 461.5 days, 0.056 per day (OUR CHECK, SYNTHETIC), about one false alarm per 18 quiet days, which is what a 5% per day budget gives. Say: "One false alarm every 18 or so quiet days is the price of the guarantee; the Gaussian test above fired 7 times in 9,000 steps (about one per 1,300 steps)."

Do not say "no false alarms". Do not say the guarantee holds forever after calibration: the proof covers the first day after each calibration (study TH6, limitations page).

### Step 3 (2 min): real Bitcoin around 10 October 2025 (REAL)

```bash
python -m avmonitor run --source binance --symbol BTCUSDT --download --start 2025-10-08 --end 2025-10-14 --plot btc.png
```

Expected output (run today):

```
[REAL] binance BTCUSDT, 2025-10-08 00:00:00+00:00 to 2025-10-13 23:59:00+00:00, 8640 bars
detector robust, A = 28800 (false alarm budget 5% per day), band 0.25
4 alarms in 6.0 days (0.67 per day)
2025-10-09 09:23:00+00:00  up    563 min since calibration  next_day_vol_over_baseline 0.85
2025-10-10 15:30:00+00:00  up    366                        2.07
2025-10-11 16:57:00+00:00  down   86                        0.75
2025-10-12 21:36:00+00:00  up    278                        1.08
```

Note the end date is exclusive (the data stop at 13 October 23:59), although the pitch text says `--end 2025-10-14`.

Open `btc.png`. Top panel: price. Bottom: rolling volatility. Red lines: alarms that volatility is up, green: down. Say: "The 10 October liquidation crash is the spike. The monitor said 'up' at 15:30 UTC that day, about four and a half hours before the 20:00 UTC reference time of the study; the last column says volatility over the following day was 2.07 times the baseline, but that day includes the crash itself. It did not predict the crash and I do not claim it did."

Two honest caveats to say out loud: (1) four alarms in six days, at least two (9 Oct, 12 Oct) have next-day volatility near the baseline, so they are the kind of alarm the 5% per day budget allows; (2) real data alarms add little over a good forecast (study E20), and the one fresh confirmatory real-data test did not succeed (E23, 2 of 16 releases, p = 0.25).

What the chart means for a beginner: price on top, how wild the minute-to-minute moves are underneath; a vertical line means "the wildness is now different from what this monitor was calibrated on".

Optional, Nifty (REAL), 2.2 s:

```bash
python -m avmonitor run --source nifty --start 2024-06-03 --end 2024-06-07 --download --plot nifty.png
```

Expected: `1500 bars`, `2 alarms in 4.0 days`, `2024-06-04 09:43 up` (5.8 times baseline), `2024-06-05 11:09 down`. The replay page shows the alarm at 09:26 for the same day because it starts from a week of history; say so if asked.

### Step 4 (1.5 min): the replay page (REAL)

Open `docs/demo/replay.html` in a browser (Chrome, Safari, Firefox; no server needed).

1. Tab "N1: Nifty 50, 4 June 2024" (election-result crash) is selected. Click Play, or use the "Jump to" menu: "e-detector up alarm, 4 Jun 09:26 IST".
2. Panel 1 shows the clipped squared return (grey bars; dark bars are the ones the monitor is using). Panel 2 shows the monitor's "wealth" on a log scale; an alarm fires when it crosses the dashed line (2A = 15,000).
3. Click "C3: BTC, 10 October 2025" to see the same on Bitcoin. The page says in its own words that this is one pre-registered case, not evidence that the method is better.

Say: "Each bar is one minute. The monitor places a bet that volatility is higher than baseline; if it is right the wealth explodes and crosses the line. A hatched stretch means it was recalibrating after its own alarm and could not alarm: a blind spot I show on purpose."

I checked the page against the committed results: for N1 every method's alarm count equals `results/e18_alarms.csv` (primary design), and for C3 on BTC the e-detector's 5 alarms equal the CSV's 5. Dark mode and the step, play and jump controls all worked, in desktop and phone widths (screenshots in `reports/screenshots/replay_*`).

### Step 5 (optional, 1 min): examples

```bash
python examples/01_quickstart.py
python examples/03_restart_every_day.py
```

Expected (SYNTHETIC): `streaming alarm: step 3023, direction up, 24 steps after the change`; and for example 3 the line `daily restart, known baseline: 0 alarm(s) before the change; first alarm after it on day 7, minute 216 (up), 17 minutes after the change; 39 alarm(s) after it`, explained by the note printed beneath it (a known baseline is never re-estimated, so it keeps alarming after a lasting change).

Do not run `examples/04_panel_ebh.py` live unless you explain it: it prints a deliberate false discovery ("asset_0 ... false discovery (no change)"), which is allowed by the false-discovery-rate rule but needs a sentence.

### Step 6: the closing sentence

"A proved false-alarm budget, a robust statistic, and honest limits: it is calibrated per day, it is blind while recalibrating, and the real-data evidence is mixed."

## Error states I tried

| Input | What happened | Say / do |
|---|---|---|
| `--start 2025-10-14 --end 2025-10-08` (reversed) | `no data in the requested period`, exit 1 | "Start must be before end." |
| `--start 2025-10-08 --end 2025-10-08` (same day) | same message | "End date is exclusive; use the next day." |
| `--start 2024-06-08 --end 2024-06-09 --source nifty` (Saturday) | `no data in the requested period` | "Market closed." |
| Binance on a weekend day | works: crypto trades 24/7, `1440 bars`, `0 alarms in 1.0 days` | n/a |
| `--start not-a-date` | a long Python traceback ending `DateParseError: Unknown datetime string format` | Avoid typing dates by hand; copy the commands. Reported as a finding. |
| unknown symbol `NOPEUSDT` or a future date with `--download` | traceback ending `urllib.error.HTTPError: HTTP Error 404: Not Found` | Same. |
| an empty data file (I emptied the 2025-10 file and restored it byte for byte) | traceback ending `pandas.errors.EmptyDataError: No columns to parse from file` | Delete the empty file and rerun with `--download`. |
| `python -m avmonitor estimators` (as written in the README) | `error: the following arguments are required: --day` | The README line is wrong; use `estimators --source binance --symbol BTCUSDT --day 2025-10-10 --download`, which prints ten estimators. The last, "as printed (raw MAD, eq. 16)", is 5.6e-13 against 3.5e-07 for the scaled version: the printed formula is on a different scale (see `docs/estimator_issues.md`); skip this command in a short demo. |

## Fallback if live execution or the internet fails

1. `docs/demo/replay.html` (self-contained, offline) covers the real-data story.
2. `notebooks/demo.ipynb` is executed (6 code cells with outputs, 0 error outputs, 3 figures); scroll through it.
3. `results/figures/demo_monitor.png` is the synthetic demo figure (`python experiments/demo.py` rewrites it; I ran it and `git status` showed the file unchanged, so the committed picture is what the code produces). The numbers in the notebook and that script (`3020 up, 6145 down, 8308 down`) use a different seed from the command-line demo (`3010, 6098, 8287`); both are right, do not mix them in one breath.
4. Committed result tables for the numbers: `results/e18_alarms.csv`, `results/e18_outcomes.csv`, `results/e09ab_bound_tightness.csv`.
