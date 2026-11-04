# R1. Real-data integrity and clean reproduction (round 4)

Date: 2026-11-04. Checkout: main at 4de8a32. All results below are [REAL]: they use the downloaded one-minute bars, not simulations. Scripts, logs and CSVs are in `reports/round4/R1/`.

## 1. Summary

* **All 139 fixed data files are exactly the committed ones.** 120 Binance months (BTC, ETH, XRP, SOL, DOGE, 2024-01 to 2025-12) and 19 Nifty 50 files match `data/MANIFEST.csv` in SHA-256 and row count. The 30 fresh 2026 Binance months that E23 used (BTC, ETH, SOL, 2026-01 to 2026-10) match `results/e23_data_files.csv`. VERIFIED.
* **The daily-updated Bitstamp file still contains the committed data byte for byte.** Its first 956,327 rows (header included, 53,449,859 bytes) hash to the committed SHA-256; the committed last bar is 2026-11-02 02:47 UTC. The file has no duplicate timestamps and no missing minutes from 2025-01-07 to 2026-11-02 (upstream fills missing minutes with flat zero-volume candles; there are 6,812 zero-volume bars, 0.7% of the rows). VERIFIED.
* **No bad prints found.** Across all 170 files: 0 duplicate timestamps, 0 out-of-order bars, 0 non-positive or NaN prices, 0 bars with high below low or open or close outside the range. Of the 1,700 largest one-minute returns (10 per file), none is larger than the high-low span of its two bars, and the biggest ones (XRP minus 24.8% and DOGE minus 26.9% on 2025-10-10 21:20 UTC, ETH minus 10% on 2025-02-03, BTC and Bitstamp minus 2% at 2026-01-31 18:43) appear on several coins or venues at the same minute with large volume. They are real market moves.
* **One data problem the pipeline must keep handling: the 2026-08 to 2026-10 Nifty files freeze the last 14 minutes of every session.** 97% of the bars from 15:15 to 15:28 IST are flat (open = high = low = close), and the 15:29 bar then carries the whole closing move. No earlier file has this. E7's cleaning rule (drop runs of 5 or more zero returns and the return after them) removes these minutes, so the live Nifty pre-registration monitors about 360 returns per session, not the 375 its config states. Severity: medium (documentation and vendor dependence).
* **Nifty files hold three exchange outage days with stale prices** (2017-07-10, 2021-02-24 and the 2020-03-13 circuit halt), one bogus mini-session (2021-11-04, a holiday, 43 identical bars), three Muhurat sessions outside normal hours and two Saturday test sessions. E7 already flags the frozen runs; the rest are dropped by `regular_sessions`.
* **The 2026-10 Binance files end on 2026-10-24 23:59 UTC** (34,560 rows, 10,080 minutes short). The mirror still has no later data today. Because of this E23 excluded the 28 October 2026 FOMC statement (recorded in `results/e23_events_excluded.csv`), so the README's "16 releases" is 16 of the 17 that occurred.
* **Clean reproduction: E22, E23 and the E25 gates reproduce the committed results exactly** (every file identical except the wall-clock "seconds" column; E25's 34 gate checks pass with 0 mismatches). E18 was stopped at the 25 minute budget after writing one identical file; E19 p1 was skipped as hours long. See section 3d.

## 2. What was run

Environment: `/home/user/.venv-bapat`, one core (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`).

Part A (main checkout, writes only to `reports/round4/R1/`):

```
python reports/round4/R1/download_2026.py            # E23.download_months for BTC, ETH, SOL: 2025-12 to 2026-10 (30 new files)
python reports/round4/R1/integrity_all.py reports/round4/R1   # all checks; writes integrity_files.csv, integrity_top10.csv, integrity_nifty_sessions.csv
python reports/round4/R1/make_tables.py               # tables.md
head -n 956328 data/raw/crypto/btcusd_bitstamp_1min_latest.csv | sha256sum   # Bitstamp prefix against the committed hash
```

Part B (separate worktree `/home/user/wt-repro` at 4de8a32, `data/raw` symlinked, `PYTHONPATH=/home/user/wt-repro/src`; the main checkout's `results/` was never touched):

```
python experiments/e22_guaranteed_events.py run --procs=1      # 2 min 44 s
python experiments/e23_fresh_releases.py run --procs=1         # 2 min 12 s
python experiments/e18_landmark_events.py run figure context awake   # the script has no --procs flag and uses its built-in Pool(2)
python experiments/e25_sequential_events.py gates --procs=1
python experiments/e19_real_null.py pools p1
python scripts/reproduce.py --steps compare                    # after each experiment
```

The runner is `reports/round4/R1/run_chain.sh`; logs are in `reports/round4/R1/logs/`.

## 3. Results

### 3a. Data integrity, one table per source

How to read these tables: "SHA-256 match" and "Rows match" compare the file on disk with the committed `data/MANIFEST.csv` (for the 2026 months, with `results/e23_data_files.csv`). "Dup. ts" is repeated timestamps. "Missing minutes" counts minutes with no bar inside the file's period (for Nifty, inside the regular 09:15 to 15:29 IST session of each trading day present). "Bad prices" are non-positive or NaN open, high, low or close. "Flat bars" have open = high = low = close. "Zero-return share" is the share of one-minute log returns that are exactly zero. "Robust sd" is 1.4826 times the median absolute return, a spread measure that heavy tails do not inflate; "Returns > 20 robust sd" counts one-minute moves larger than 20 of these units. "Top-10 flagged" is how many of the file's 10 largest moves were at least half reversed in the next minute (a rough warning sign, not proof of an error; see 3b).

**Binance, committed manifest months (2024-01 to 2025-12)**

| Coin | Files | SHA-256 match | Rows match | Rows | Dup. ts | Missing minutes | Bad prices | OHLC inconsistent | Zero-vol bars | Flat bars | Zero-return share (min to max) | Returns > 20 robust sd | Largest return / sd | Largest abs. return | Top-10 flagged |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BTCUSDT | 24 | 24/24 | 24/24 | 1,052,640 | 0 | 0 | 0 | 0 | 0 | 0 | 0.9% to 10.5% | 155 | 115 | 4.4% | 29/240 |
| DOGEUSDT | 24 | 24/24 | 24/24 | 1,052,640 | 0 | 0 | 0 | 0 | 1 | 169 | 1.0% to 11.5% | 216 | 329 | 26.9% | 44/240 |
| ETHUSDT | 24 | 24/24 | 24/24 | 1,052,640 | 0 | 0 | 0 | 0 | 0 | 0 | 0.4% to 2.9% | 178 | 133 | 10.1% | 46/240 |
| SOLUSDT | 24 | 24/24 | 24/24 | 1,052,640 | 0 | 0 | 0 | 0 | 0 | 0 | 2.0% to 6.9% | 122 | 61 | 7.4% | 46/240 |
| XRPUSDT | 24 | 24/24 | 24/24 | 1,052,640 | 0 | 0 | 0 | 0 | 0 | 118 | 1.3% to 24.3% | 296 | 363 | 24.8% | 38/240 |

**Binance, E23 fresh months (2026-01 to 2026-10), checked against results/e23_data_files.csv**

| Coin | Files | SHA-256 match | Rows match | Rows | Dup. ts | Missing minutes | Bad prices | OHLC inconsistent | Zero-vol bars | Flat bars | Zero-return share (min to max) | Returns > 20 robust sd | Largest return / sd | Largest abs. return | Top-10 flagged |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BTCUSDT | 10 | 10/10 | 10/10 | 427,680 | 0 | 10,080 | 0 | 0 | 0 | 0 | 0.8% to 7.4% | 65 | 88 | 2.3% | 10/100 |
| ETHUSDT | 10 | 10/10 | 10/10 | 427,680 | 0 | 10,080 | 0 | 0 | 0 | 0 | 0.5% to 4.2% | 95 | 105 | 4.5% | 9/100 |
| SOLUSDT | 10 | 10/10 | 10/10 | 427,680 | 0 | 10,080 | 0 | 0 | 0 | 1 | 6.0% to 13.7% | 60 | 139 | 5.7% | 7/100 |

**Nifty 50 (data/raw/nifty, one file per year, then per month in 2026)**

| File | SHA-256 | Rows (= manifest) | Sessions | Short sessions | Missing minutes in sessions | Bars outside 09:15 to 15:29 | Flat bars | Zero-return share | Returns > 20 sd | Largest return / sd | Largest abs. return | Top-10 flagged |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NIFTY50_1min_2017.csv | yes | 69,375 (yes) | 185 | 0 | 0 | 0 | 166 | 1.9% | 3 | 24 | 0.4% | 1/10 |
| NIFTY50_1min_2018.csv | yes | 91,874 (yes) | 245 | 1 | 1 | 0 | 0 | 1.7% | 4 | 41 | 1.0% | 0/10 |
| NIFTY50_1min_2019.csv | yes | 91,500 (yes) | 244 | 0 | 0 | 0 | 0 | 1.5% | 4 | 33 | 0.8% | 1/10 |
| NIFTY50_1min_2020.csv | yes | 94,125 (yes) | 251 | 0 | 0 | 0 | 101 | 1.0% | 26 | 160 | 6.0% | 3/10 |
| NIFTY50_1min_2021.csv | yes | 92,668 (yes) | 248 | 1 | 332 | 0 | 346 | 1.3% | 0 | 18 | 0.5% | 1/10 |
| NIFTY50_1min_2022.csv | yes | 92,617 (yes) | 247 | 6 | 8 | 0 | 38 | 0.5% | 0 | 18 | 0.5% | 1/10 |
| NIFTY50_1min_2023.csv | yes | 91,931 (yes) | 245 | 5 | 5 | 61 | 1 | 0.6% | 5 | 33 | 0.6% | 0/10 |
| NIFTY50_1min_2024.csv | yes | 92,523 (yes) | 248 | 2 | 537 | 60 | 3 | 0.3% | 16 | 38 | 0.9% | 2/10 |
| NIFTY50_1min_2025.csv | yes | 93,060 (yes) | 249 | 1 | 315 | 0 | 0 | 0.4% | 3 | 44 | 0.9% | 1/10 |
| NIFTY50_1min_2026-01.csv | yes | 7,500 (yes) | 20 | 0 | 0 | 0 | 0 | 0.3% | 0 | 14 | 0.3% | 1/10 |
| NIFTY50_1min_2026-02.csv | yes | 7,875 (yes) | 21 | 0 | 0 | 0 | 0 | 0.4% | 2 | 28 | 0.6% | 2/10 |
| NIFTY50_1min_2026-03.csv | yes | 7,125 (yes) | 19 | 0 | 0 | 0 | 0 | 0.2% | 0 | 7 | 0.3% | 1/10 |
| NIFTY50_1min_2026-04.csv | yes | 7,500 (yes) | 20 | 0 | 0 | 0 | 0 | 0.3% | 0 | 9 | 0.2% | 0/10 |
| NIFTY50_1min_2026-05.csv | yes | 7,125 (yes) | 19 | 0 | 0 | 0 | 0 | 0.3% | 0 | 14 | 0.4% | 1/10 |
| NIFTY50_1min_2026-06.csv | yes | 7,875 (yes) | 21 | 0 | 0 | 0 | 0 | 0.4% | 0 | 13 | 0.3% | 0/10 |
| NIFTY50_1min_2026-07.csv | yes | 8,625 (yes) | 23 | 0 | 0 | 0 | 0 | 0.5% | 0 | 14 | 0.3% | 2/10 |
| NIFTY50_1min_2026-08.csv | yes | 7,875 (yes) | 21 | 0 | 0 | 0 | 294 | 4.0% | 5 | 59 | 0.8% | 0/10 |
| NIFTY50_1min_2026-09.csv | yes | 7,875 (yes) | 21 | 0 | 0 | 0 | 294 | 3.9% | 1 | 22 | 0.4% | 0/10 |
| NIFTY50_1min_2026-10.csv | yes | 7,500 (yes) | 20 | 0 | 0 | 0 | 280 | 3.9% | 0 | 13 | 0.4% | 1/10 |

**Bitstamp BTC/USD (daily-updated file, rows up to the committed last bar 2026-11-02 02:47 UTC)**

| Check | Value |
|---|---|
| Rows on disk today | 957,844 |
| Rows up to 2026-11-02 23:59 UTC | 957,599 |
| Committed manifest rows | 956,327 |
| Last bar of the first 956,327 rows | 2026-11-02 02:47 UTC |
| SHA-256 of header + first 956,327 rows | 810a84b5... = committed manifest (byte-exact prefix, 53,449,859 bytes) |
| First bar | 2025-01-07 00:01 |
| Duplicate timestamps | 0 |
| Missing minutes (gaps) | 0 (0 gaps) |
| Strictly increasing, whole minutes | True, True |
| Bad prices, OHLC inconsistent | 0, 0 |
| Zero-volume bars (filled minutes) | 6,812 |
| Flat bars | 44,194 |
| Zero-return share | 5.7% |
| Returns > 20 robust sd | 223 |
| Largest return / sd, largest abs. return | 73, 2.4% |
| Top-10 flagged | 1/10 |

### 3b. The 10 largest moves of each file

How to read: a bad print (a wrong price that is corrected a minute later) would show a move larger than the two bars' high-low range, or a move that fully reverses in the next minute with little volume. A real move stays inside the range, has volume, and shows up at the same minute on other coins or venues.

**The 10 largest absolute one-minute returns of every file (1,700 returns in all)**

| Source | Returns looked at | Range covers the move and no instant reversal | At least half reversed in the next minute | Move larger than the two bars' high-low span | Largest return / robust sd |
|---|---|---|---|---|---|
| BTCUSDT | 340 | 301 | 39 | 0 | 115 |
| Bitstamp | 10 | 9 | 1 | 0 | 73 |
| DOGEUSDT | 240 | 196 | 44 | 0 | 329 |
| ETHUSDT | 340 | 285 | 55 | 0 | 133 |
| Nifty 50 | 190 | 172 | 18 | 0 | 160 |
| SOLUSDT | 340 | 287 | 53 | 0 | 139 |
| XRPUSDT | 240 | 202 | 38 | 0 | 363 |

Flagged reversals are concentrated in known turbulent minutes (10 October 2025 liquidation cascade, 5 August 2024, 3 February 2025, 31 January 2026 18:43 UTC, the Nifty circuit halt of 13 March 2020) and appear on several coins and on both BTC venues (Binance and Bitstamp) at the same minute. The largest, XRP minus 24.8% then plus 16.6% within two minutes on 2025-10-10 21:20 UTC, has a 30% high-low range and 5.6 million XRP of volume in that bar. I found no case that looks like a bad print. The high zero-return shares of XRP in 2024 (up to 24% of minutes) come from the tick size: 0.0001 USDT is 0.02% of the price, about two thirds of the robust sd, so many minutes end on the same tick. Status: VERIFIED (no data errors among the 1,700 largest moves).

### 3c. Nifty 50 oddities, one by one

How to read: a regular NSE session has 375 one-minute bars, 09:15 to 15:29 IST. Every row below is a day where the file departs from that. Full list of sessions with missing minutes: `integrity_nifty_sessions.csv`.

| Day | What the file shows | Likely cause | Handled by the code? |
|---|---|---|---|
| 2017-07-10 | 166 flat bars, 09:25 to 12:29 | NSE technical outage (trading halted 09:55 to 12:30) | yes, E7 frozen-run rule (`results/e07a_frozen_runs.csv`) |
| 2020-03-13 | 51 flat bars 09:22 to 10:20, then +6.0% in the 10:21 bar | lower circuit, 45-minute halt, rally on reopening | yes, frozen run plus the return after it is dropped |
| 2021-02-24 | 303 flat bars, 10:07 to 15:29 (70 distinct closes all day) | NSE outage (index feed frozen from about 10:08, halt from 11:40) | yes, frozen-run rule |
| 2021-11-04 | 43 bars 14:47 to 15:29, all with one price 17829.2 | Diwali holiday; vendor artifact (Muhurat session was 18:15 to 19:15) | yes, `regular_sessions` drops sessions under 370 bars |
| 2023-11-12, 2024-11-01, 2025-10-21 | 61, 60 and 60 bars outside or short of normal hours | Muhurat trading sessions | yes, short sessions dropped |
| 2024-03-02, 2024-05-18 | 107 and 106 bars, 09:15 to 12:30 | Saturday disaster-recovery test sessions | yes, short sessions dropped |
| 2022 and 2023 (11 days), 2018-07-09 | one missing minute (usually 15:29); 2022-03-07 three | vendor | the next return spans two minutes; E7 counts it as a gap |
| 2026-08, 2026-09, 2026-10 (every day) | 14 flat bars 15:15 to 15:28, closing move in the 15:29 bar; a `Volume` column appears | new vendor format | yes, frozen-run rule drops 15 returns per session (`results/e07a_frozen_runs.csv` lists them, for example 2026-10-29 15:15, 14 minutes); the pre-registered config says 375 steps per session and the dry run of 2026-11-02 reports 360 returns |

Also: the 2017 file starts on 2017-04-03, so "2017" is nine months. The `Timestamp` text column agrees with the `Epoch` column (IST) in every row of every file.

### 3d. Clean reproduction (Part B)

How to read: "identical" means the SHA-256 of every result file the experiment writes equals the committed one. "Timing only" means the only differing column is the measured wall-clock time. A headline number in README.md could only change if a file were "different".

| Entry | Command (in the worktree) | Wall time | Result vs committed | README number affected |
|---|---|---|---|---|
| E22 | `e22_guaranteed_events.py run --procs=1` | 2 min 44 s | all 6 result files and the figure rewritten: 5 identical, `e22_monitor_summary.csv` equal except `seconds` | none |
| E23 | `e23_fresh_releases.py run --procs=1` | 2 min 12 s | 17 result files and the figure rewritten: 16 identical, `e23_monitor_summary.csv` equal except `seconds` (`e23_power*.csv` belong to the `power` entry, not rerun) | none |
| E18 | `e18_landmark_events.py run figure context awake` | stopped after 27 min (started 02:04 UTC, killed 02:31 UTC by `run_chain2.sh`: over the 25 minute budget; the Nifty stream alone took 5 min 12 s and three crypto streams of five coins follow) | PARTIAL: only `e18_events.csv` was written before the stop, identical to the committed file; the alarm, outcome and rate files were not rerun | README line 42 not checked |
| E25 gates | `e25_sequential_events.py gates --procs=1` | 2 min 50 s | `e25_gates.csv` rewritten and identical to the committed file; the log reports gates 2 and 3: 34 checks, 0 mismatches (own runs reproduce the committed E22 and E23 alarm files) | none; VERIFIED |
| E19 pools (p1 skipped) | `e19_real_null.py pools` | 18 s | `e19_pools.csv` (the calm-day pools of the pre-registered null) rewritten and identical to the committed file; VERIFIED. p1 skipped: 2,000 streams per cell with 2 workers is hours, far over the 25 minute budget | README line 36 rests on p1, not rerun |

`scripts/reproduce.py --steps compare` after E22, E23 and E25 gates: "463 identical, 0 equal within rtol, 2 equal except timing columns, 0 different, 0 missing". E22's gate ("P0 reproduces E12 exactly", 60 rows and all columns, max difference 0) and E23's gates against E12 and E22 passed in the rerun. Status: VERIFIED for E22, E23, the E25 gates and the E19 pools; PARTIAL for E18; E19 p1 not run.

## 4. Defects and claims that need changing

| Where | What | Severity |
|---|---|---|
| Live Nifty pre-registration (`experiments/live_nifty_preregistered.py` config, `results/live/dryrun_2026-11-02/summary.json`: `steps_per_session: 375`) | Since August 2026 the Nifty source freezes 15:15 to 15:28 IST and books the closing move at 15:29. After E7 cleaning the session has about 360 returns, and the last 15 minutes are never monitored. The pre-registration text should say so, and the forward run should check each day that the frozen run is exactly the closing window and not an outage. | Medium |
| README.md line 41, E23 | "16 FOMC statements and CPI releases of January to October 2026": 17 occurred; the 28 October FOMC was excluded because the Binance mirror's October file stops at 2026-10-24 23:59 UTC (still true today). Say "16 of the 17" and give the reason. | Low |
| `data/README.md` | Says the committed Bitstamp copy has "last bar 2026-11-02"; add the exact last bar (02:47 UTC) and that the committed data is the first 956,327 rows, so a reader can verify the hash with `head -n 956328 ... | sha256sum`. | Low |
| `data/README.md`, Nifty section | Does not mention that the 2026 monthly files carry a `Volume` column and a different closing-minute convention from the yearly files. | Low |
| `experiments/e18_landmark_events.py`, `e19_real_null.py` | No `--procs` option; they always start 2 worker processes, so a one-core reproduction is not possible as `scripts/reproduce.py` suggests. | Low |

No README claim that rests on E22 or E23 changes. No data file differs from the committed manifest.

## 5. What this does not show

* The checks test internal consistency and agreement with the committed hashes. They cannot tell whether the Binance mirror or the Nifty vendor altered prices before publishing (for example whether the 2026 closing-minute freeze is the vendor's or the exchange's).
* "Real move" verdicts rest on range, volume and cross-venue agreement at the same minute; they are not a tick-by-tick audit against the exchange.
* Only the pre-registered real-data entries listed in 3d were rerun; E20, E21 real parts, E24 real parts and the E19 post hoc parts were not, and the forward tests (E24 future, E25 status) were deliberately not run.
* The Bitstamp check verifies the committed prefix; experiments cut this file to dates ending 2026-01-01 at the latest, so the rows added since 2 November do not enter any result.
