# Visual QA (Bapat prototype)

Written by an AI agent (Claude), 3 November 2026, validation agent 3. Labels: REAL = real market data, SYNTHETIC = simulated, OUR CHECK = recomputed by me in this session.

## Summary for a beginner (5 lines)

1. I checked what a professor will see: the replay web page, the live demo chart, the paper's three figures and a sample of the report figures (about 12 of the 51 report figures, plus the page and the notebook).
2. The replay page is in good shape: works offline, no errors, fits a phone screen, and its alarm lists match the committed result files.
3. The figures are labelled REAL or SYNTHETIC in their titles almost everywhere, and the paper figures state "Simulated data" in their captions.
4. I made one small fix in `experiments/make_figures.py`: two numbers typed into figure text (the bound 500/7800 and the factor 0.705) are now taken from the code that produced the data. Only the title pixels changed.
5. The bigger demo risks are in the command line, not the charts: a README command that fails as written, raw Python tracebacks for a mistyped date, and the budgeted false alarms that look wrong if not explained.

## What I checked and how

| Check | Method | Result |
|---|---|---|
| Inventory | figures by script: `make_figures.py` writes 3 figures (E1 two, E4 one); the other 48 report figures come from the individual `experiments/eNN_*.py`; paper figures from `paper/make_paper_figures.py`; `alternatives/*/results` hold 19 more | which doc embeds which figure listed from `grep`; README and PITCH embed none; `docs/experiment_report.md` names 41 distinct figure files |
| Labels | AST scan of every figure function for REAL or synthetic text | all report figures carry one except `e24_near_blind.png` (exact formula, no data) and the paper's `validity.png` (caption says "Simulated data") |
| Typed numbers | compared title and line constants with the code that made the data | E4 figure had `500 / 7800` typed (equal to `nu = 500` and `A = 7800.0` in `e04_detection.py`) and E1 had `0.705` typed (equals `truncated_moment_factor` = 0.7045): now computed |
| Plotted values vs files | paper Figure 1 vs `results/e09ab_bound_tightness.csv` and `e09c_garch_band.csv` | Gaussian, 390 steps: rate over bound 0.225 (plot about 0.22); GARCH persistence 0.95 share of alarm days 0.7815 at band 0 and 0.1625 at band 0.25 (plot 0.78 and 0.16, text says 78%): OUR CHECK, SYNTHETIC, all equal |
| Replay page | Chromium 1280x800 and 390x844, offline, clicked Step, Play, Jump, tab switch, theme toggle | 0 console errors, 0 external requests, 0 px horizontal overflow; embedded alarm lists vs `results/e18_alarms.csv` (primary design): all ten methods for N1 equal; BTC e-detector in C3: 5 equal 5 |
| Live demo chart | ran the BTC command and viewed `btc.png` | title says [REAL], red and green lines explained in the title, price and volatility panels aligned |

## Findings

| Figure or page | Issue | Severity | Status | Evidence |
|---|---|---|---|---|
| `docs/demo/replay.html` | none found: [REAL] tags, a plain "one pre-registered case study, not evidence that the method is better" statement, recalibration blind spots drawn on purpose, labelled controls, dark mode works | none | VERIFIED (alarm counts for N1 all methods, BTC e-detector in C3; the other C3 methods and the other assets were not compared) | `reports/screenshots/replay_desktop.png`, `replay_phone.png`, `replay_desktop_after_jump.png`, `replay_phone_c3_dark.png` |
| `experiments/make_figures.py` E4 and E1 | typed constants `500 / 7800` and `0.705` in figure text | Low | FIXED: read from `e04_detection.py` and `robust.py`; regenerated `e04_false_alarms.png` and `e01_scale_factor.png`; pixel diff limited to the title rows (E4 title now "E4 (synthetic): dashed line = guaranteed bound 500 / A = 0.064 for 500 pre-change steps, A = 7800", E1 title adds "(formula check, no data)"), `e01_table7_reproduction.png` identical | `git diff --stat results/figures` shows only those two PNGs |
| `paper/figures/validity.png` | The picture itself does not say SYNTHETIC; the caption in the paper does ("Simulated data."). Category x axis (98, 390, 1950, 3900, 7800) is equally spaced though the values are not | Low | OPEN, owner of `paper/make_paper_figures.py` (not in my ownership) | `paper/figures/validity.png` |
| `results/figures/e23_release_rates.png`, `e18_events.png` | Honest (the failed confirmatory test is shown with p values, controls and the 2024 to 2025 reference), but panel text is about 6 px at report page width | Low | OPEN | viewed; E23 primary line reads "2/16, watched 11, p 0.25 (primary)" |
| `python -m avmonitor run ... --plot btc.png` | Rolling volatility axis has no unit | Low | OPEN | `reports/screenshots/live_btc_2025-10-08_to_14.png` |
| `README.md` line 87 | `python -m avmonitor estimators` is documented as "all quarticity estimators on one simulated day" but fails: `error: the following arguments are required: --day` | Medium | OPEN, owner Agent 2 (README) and Agent 1 (CLI) | `reports/DEMO_WALKTHROUGH.md` error table |
| `python -m avmonitor run` | A mistyped date, unknown symbol, future date or an empty data file print raw Python tracebacks (`DateParseError`, `HTTPError 404`, `EmptyDataError`); a reversed or empty range prints only `no data in the requested period` | Medium | OPEN, owner Agent 1 | commands in `DEMO_WALKTHROUGH.md` error table |
| `python -m avmonitor simulate --ratio 2` and `experiments/demo.py` | The robust detector prints a third alarm (step 8287 down) when the stream has been calm since step 6000. It is the budgeted false alarm (5% per day): 26 alarms in 20 seeds of 9000 steps with no change at all, 0.056 per day (OUR CHECK, SYNTHETIC). The CLI never says so, so it looks like a failure live | Medium | OPEN: documented in the walkthrough; suggest one printed line "about one false alarm per 18 quiet days is expected" | `reports/DEMO_WALKTHROUGH.md` step 2 |
| `PITCH_BAPAT.md` lines 22 and 86 | "About 240 tests" but the baseline run has 581 passed | Medium | OPEN, owner Agent 2 | `baseline_bapat_pytest.log` tail: 581 passed |
| `python -m avmonitor run` summary line | "1 alarms in 2.0 days" (plural), mixed number formats (0.00205 vs 9.542e-06) in `estimators` | Low | OPEN | observed output |
| `notebooks/demo.ipynb` | Executed: 6 code cells with outputs, 0 error outputs, 3 figures; synthetic output tagged `[synthetic]`; the real-data cell says "real data files found: True", so it was built with data present and would need the download on a fresh machine | Low | VERIFIED | notebook JSON inspected |
| `experiments/demo.py` | Rewrites the committed `results/figures/demo_monitor.png` when run. Running it left `git status` clean for that file, so the committed picture is exactly what the code produces | Low | VERIFIED | `git status --short results/figures` empty after the run |
| Paper PDF `paper/short_paper.pdf` | 72 pages, 3 embedded raster images found on the figure pages; text figure captions say "Simulated data" where relevant. I did not render every page | Low | PARTIALLY VERIFIED | `pdfinfo`, `pdftotext` |

## Not done

I did not view `e05` to `e17` and `e19` to `e25` individually, nor any of the 21 figures under `alternatives/`. The scan above covers their labels only. The replay's C3 alarm lists for ETH, SOL, XRP and DOGE were not compared with the CSV (the page shows BTC only).
