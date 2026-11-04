# Changelog

All notable changes to the `avmonitor` package and its research record are listed here, newest first.
The format follows the Keep a Changelog convention, and version numbers follow semantic versioning.

No version has been tagged or published yet. Versions 0.1.0 and 0.2.0 group the work of the first two
rounds as recorded in `git log`; the date of each is the UTC date of its last commit. The rounds partly
overlapped in time, so some 0.3.0 work (the first paper drafts, E17) has commits dated before the end of
0.2.0.

## [0.3.0] - Unreleased

Round 3: theory TH4 to TH6, experiments E16 to E25, the paper, and the preparation of a release.

### Fixed (4 November 2026, after the validation audit)

- Loaders: the Bitstamp loader keeps the first row read for a duplicated minute (stable sort) and reports
  its checks; the Bitstamp and Nifty loaders mask nonpositive prices.
- `VolatilityMonitor` warns for `cal_len` below 100, when it never finished calibrating, and when
  `null_band` is 1 or more; `update` rejects a non-number at once.
- `gaussian_lr_edetector` and `page_cusum` skip missing values instead of freezing at NaN.
- The command line prints a short `error:` line and exits with status 2 instead of a traceback.
- No stored result changed (465 of 465 result files identical after the quick tier and an E4 rerun).

### Fixed and added (4 November 2026, round 4, on top of the entries above)

- `VolatilityMonitor` and its subclasses: `threshold` must be a finite number above 1 (was: positive, and
  infinity was accepted, which alarmed once a statistic overflowed to inf, R5-B-002); `c` and `c_down` must
  be finite and positive, with an error that names the argument (a non-number `c` gave a `TypeError`, and
  `c_down = inf` was accepted).
- New `VolatilityMonitor.status()` and the attributes `cal_steps` and `restarts_in_a_row`: they show at any
  moment whether the monitor is blind (steps spent calibrating, calibration restarts in a row, whether the
  increase side can gain at all). Nothing about the alarms changes.
- `update` of `GaussianLRMonitor`, of `PeriodicVolatilityMonitor` with a period or the CUSUM statistic, and
  of the guard's rank fallback rejects a non-number at once, as `VolatilityMonitor.update` already does.
- `repeated_window_ztest` skips missing values (NaN or infinite) like the other two benchmarks: the window
  holds the last observed returns. Without missing values the output is unchanged.
- The Bitstamp, Nifty and Binance loaders also mask bars whose close is infinite or missing and report them
  as `nonfinite_prices` in `df.attrs["checks"]`.
- A calibration window of huge returns (median above about 1e153) gave an `OverflowError` from
  `(c * scale) ** 2`; it now gives `tau^2 = inf`, which the degenerate-window rule rejects (R5-B-001). The
  same holds for `detectors.calibrate` and `GaussianLRMonitor.update`. Finite squares are computed exactly
  as before.
- No stored result changed: quick tier 465 of 465 identical, and the loaders and seven monitor
  configurations give byte-identical output on the real files (`reports/round4/M/REPORT.md`).

### Added

- Release preparation: complete package metadata in `pyproject.toml` (license, author, keywords,
  classifiers, project URLs, dependency lower bounds); the version is read from `avmonitor.__version__`
  only; an MIT `LICENSE`; this changelog.
- Top-level imports with `__all__`: the monitors (`VolatilityMonitor`, `PeriodicVolatilityMonitor`,
  `MultiFrequencyMonitor`, `PanelMonitor`), the array detectors, the confidence-sequence baseline and
  the robust estimators can be imported from `avmonitor` directly. Every existing import path still works.
- `docs/api.md`, an API reference generated from the docstrings by `scripts/build_api_docs.py` (standard
  library only), and short docstrings for the public classes and methods that had none.
- `examples/`: four runnable scripts on synthetic data (quickstart, confidence-sequence baseline,
  restart at every day boundary, panel with e-BH), each run by `tests/test_examples.py`.
- E16: adversarial search for null streams, inside the null class (i.i.d. laws, shape switching, adaptive
  bold play, slack decomposition, a low bound B) and outside it (GARCH, plug-in baseline, price-level
  ticks), with an in-class price-tick control added after review round 5.
- E17: `avmonitor.multifreq`, one monitor for 1, 5 and 15 minute detectors combined by a certified
  average of e-detectors, with its proved bound (`proved_bound`, `level_for_budget`, `handicap`);
  synthetic and five-coin results and a post hoc rate-matched comparison.
- E18: pre-registered case studies of four landmark market shocks (4 June 2024, 5 August 2024, 7 April
  2025, 10 October 2025), and an exploratory, not pre-registered run of the no-blind-period monitors on
  the same events.
- E19: pre-registered checks of the guarantee on real data: streams that permute calm real days (the
  proved monitors stay below every bound; the plug-in does not over five days), block-bootstrapped calm
  days (no monitor holds its bound), the confidence-sequence baseline as the primary configuration on the
  five coins, a forecast-relative null on real data, and the Gaussian row of Table S2 with that baseline;
  a post hoc diagnostic of confidence-sequence segments whose increase side cannot alarm.
- E20: pre-registered re-analysis of the real-data inference (claim classification, Holm correction
  within families, day-block bootstrap, forecast regressions on realised variance, QLIKE).
- E21: `avmonitor.benchmarks_robust`, benchmarks not calibrated on Gaussian data only: CUSUM and
  Shiryaev-Roberts rules on the truncated statistic with simulated or per-path bootstrap thresholds, a
  Student t BOCPD with an exact run-length recursion, and CUSUM monitoring with a training sample; on
  synthetic paths, E13's paths and the five coins at matched alarm counts.
- E22: E12's pre-registered FOMC and CPI test rerun with monitors that carry a proved per-day bound
  (`PeriodicVolatilityMonitor(baseline="cs")` restarted every UTC day, the same without restart, the
  awake rank monitor), after a gate that reproduces E12 exactly; pre-registered in
  `docs/preregistration_e22.md`.
- E23: the same test on the fresh FOMC and CPI releases of January to October 2026, pre-registered
  before any 2026 crypto return was loaded (`docs/preregistration_e23.md`), with the event dates in
  `data/events/us_macro_events_2026.csv`. The primary test fails.
- E24: `avmonitor.guard`, confidence-sequence monitors (plain and restarted every day) with a proved guard
  against blind segments: a near-blind calibration window is discarded and calibration repeated with the
  error budget split over the windows tried, with a fallback to the awake rank monitor; theory in
  `docs/theory_blind_guard.md`; a forward pre-registration (`docs/preregistration_e24_future.md`) for
  releases up to 1 June 2027.
- E25: an anytime-valid version of that forward test (`experiments/e25_sequential_events.py`): an
  e-process over releases with disjoint four-week strata of control windows, pre-registered before any
  fresh data (`docs/preregistration_e25_sequential.md`, with two amendments, both made before any fresh
  data), proof in `docs/theory_sequential_events.md`; the fresh test can run from 26 November 2026.
- `docs/demo/replay.html` and `scripts/build_replay.py`: a minute-by-minute replay of two landmark days
  of E18, rebuilt from the data and checked against E18's saved alarm lists.
- The pre-registered live run on the Nifty 50 session of 3 November 2026 (`results/live/2026-11-03/`,
  result added to `docs/preregistration_nifty_2026-11-03.md`): no alarm, because the monitor was still
  recalibrating after an alarm on 2 November for the whole session.
- TH4: `avmonitor.robustness`, false alarm bounds when the null holds only approximately (mean drift,
  plug-in baseline, GARCH), checked against E16.
- TH5: `avmonitor.limits`, limits for every valid procedure (pathwise and mean delay floors after a
  collapse of volatility), the plug-in baseline when n, m and A grow together, and the decrease bet at
  its boundary.
- TH6: `avmonitor.periodic`, what Theorem B bounds after the first day; `PeriodicVolatilityMonitor`
  restarts the statistics at every period boundary and offers an e-CUSUM statistic; Theorem D restated
  with a tolerance hypothesis. After review round 8: a known-baseline mode (`baseline="known"` with `tau`
  and `m0`) that implements Theorem 1 exactly (P(false alarm in a period | past) <= m / A), and
  Corollary 1 for the confidence-sequence baseline (delta + m / A for each period of a segment,
  marginally).
- The paper (`paper/`): short paper drafted in parts, then restructured after a referee report into a
  main text of about 10,000 words and a supplement (S1 to S9); `scripts/build_paper_pdf.py` builds the
  PDF without LaTeX.
- `scripts/trace_numbers.py`, which traces every number of a Markdown file to `docs/` or `results/`.
- Independent reviews, rounds 4 (E17), 5 (E16), 6 (TH4), 7 (TH5), 8 (TH6), 9 (the package code, with
  Monte Carlo checks of every documented guarantee and reproduction scripts in `docs/review/round9_scripts/`)
  10 (E24) and 11 (E25), and three referee passes on the paper draft.
- `avmonitor.detectors.check_bets`, which rejects bets that could make a factor negative.

### Changed

- `import avmonitor` now loads the top-level API, so it imports numpy and scipy.
- E13's "BJV-type" test is renamed the two-window variance-ratio scan (inspired by Bibinger, Jirak and
  Vetter 2017, not their procedure) everywhere; its real-data tables show the scored count of every share.
- Building the package needs setuptools 77 or later (SPDX license expression in `pyproject.toml`).
- `scripts/build_report.py --committed` builds the experiment report from the report parts as committed
  at HEAD, so unfinished parts of running work never enter it.
- README, limitations and the experiment report state per-day guarantees as TH6 proves them: Theorem B
  bounds the first day after each (re)calibration; the daily restart has a per-day bound with a known
  baseline, delta + m / A per day (not given the past) with the confidence-sequence baseline, and none
  with the plug-in baseline.

### Fixed

- `avmonitor.multifreq`: calibration with missing minutes (review round 4); every result file was
  byte-identical after the fix.
- TH4: Theorem 3 restated with the maximum over time of bounds that do not respond to the data, after
  review round 6 found two counterexamples to the first statement; a vacuous test replaced.
- TH5, after review round 7: the collapse floor holds for deterministic rules (a randomised rule has mean
  delay at least 43.4 on that path), the distribution-free plug-in bound is corrected for n of at least
  10^6, and the hypotheses of Corollary 10 are restored.
- TH6, after review round 8: the e-CUSUM window bounds of Theorem 2 are stated as unconditional (given
  the past only Lemma 2 holds).
- Paper Table S1: two cells corrected to their source files, found by `scripts/trace_numbers.py`.
- After review round 9: `betting_edetector` now leaves the statistic unchanged at a missing value (it
  used to add 1, as for an observation with Y = 1, so a stream of missing values alarmed at step 2A; no
  committed experiment passed missing values to it); `running.two_sided_sr` skips missing values instead
  of stopping; `calibrate` raises on a degenerate window; the monitors reject bets outside [0, 1); the
  known-baseline `PeriodicVolatilityMonitor` uses `c_down`; `two_sided_alarm` works with mixed integer
  types; `sr_steps_constant_factor` is exact at integer boundaries; `load_bitstamp_csv` marks the
  provider's filled minutes as missing returns instead of zero returns (0.30% of minutes in July to
  December 2025).
- Reproduction, after a full clean-clone rerun of E1 to E14 (no committed number changed):
  `scripts/reproduce.py` runs E13's review round 3 parts (`offline-long`, `online-long`, `fix-nominal`,
  `cap-check`, `real-long nifty`, `real-long crypto`, `long-pooled`) before E13's figures, and E14's
  parts a and b one calibration length or one law at a time, which writes the committed
  `e14a_validity_cal*.csv` and `e14b_delay_*.csv` (the unsuffixed files made E14's figures count every
  row twice and stop with an error); `compare` reports a file as "equal except timing columns" when only
  the listed run-time columns differ, treats missing text as equal to missing text, and the run step has
  a `--to` option. `e07_nifty.py` computes the package row of `results/e07b_deseason_check.csv` with
  `intraday_scale` as it was when the check was made (pooled start, the 09:15 bar left out); today's
  function leaves the first and last seven minutes of each session without a scale, so the row did not
  reproduce. `results/e13a_offline_size.csv` reproduces within a relative 1e-14 only: its first-run rows
  were read back once with pandas' default float parser when the label "null" became "no change".

## [0.2.0] - 2026-11-03

Round 2: real data for five coins and Nifty 50, the confidence-sequence baseline, forecast-relative
monitoring, delay theory, panel monitoring, benchmarks, alternative prototypes and independent reviews
(experiments E6 to E15, theory TH1 and TH2).

### Added

- Data: one-minute loaders for five Binance spot pairs (BTC, ETH, XRP, SOL and DOGE against USDT) and
  for Nifty 50, `scripts/download_data.py` with a SHA-256 manifest, and `data.mark_frozen`, which treats
  long runs of exact-zero returns (frozen feeds) as missing.
- Command-line interface: `python -m avmonitor simulate | run | estimators`, also installed as the
  `avmonitor` command, with end-to-end integration tests.
- `avmonitor.confseq`: one-sided empirical Bernstein confidence sequences (`eb_confseq`,
  `StreamingEBConfSeq`) and `cs_baseline`. `VolatilityMonitor(baseline="cs")` and
  `calibrate(baseline="cs")` keep the false alarm guarantee with an estimated baseline
  (delta + m / A), with an optional data-dependent calibration length (`target_rel_width`).
- `VolatilityMonitor`: calibration on degenerate windows (window extension, fallback scale, `cal_stats`,
  per-segment `segments`), `margin_z`, and `c_down`, a separate truncation constant for the decrease side
  (E11f), also with the confidence-sequence baseline.
- `robust.DegenerateThresholdWarning` and `robust.tick_floor` for zero-heavy returns; a `strict` option
  for `ew_quantile_scale`.
- `avmonitor.forecast`: predictable EWMA and GARCH(1,1) volatility forecasts, used for forecast-relative
  monitoring (E10); an executed demo notebook.
- `avmonitor.design` and the delay theory (E11, TH1): growth rates, delay predictions, bet design, false
  alarm slack, the decrease-side speed limit, and first-order optimality stated as proved upper bounds.
- `avmonitor.running` (TH2, E14): a running truncation level with a guarantee (uniform bands and
  conformal ranks) and monitors with no blind period.
- `avmonitor.panel` (E15): mixture e-processes per asset, e-BH across assets with sticky discoveries,
  Bonferroni and naive rules, and the streaming `PanelMonitor`.
- `avmonitor.competitors` (E13): an adjacent-window volatility jump test (BJV-type) and Bayesian online
  change-point detection for the variance, as benchmarks.
- Experiments E6 to E15: the published tables for five coins, streaming runs in and out of sample, Nifty
  50 and India VIX, the confidence-sequence baseline, a stress test of the false alarm bound at scale,
  forecast-relative monitoring, delay theory, an FOMC and CPI event study, benchmarks, the running
  truncation level and the panel monitor.
- Seven alternative prototypes under `alternatives/` (ALT1 to ALT6 in the report), with a ranking.
- Reproducibility: `scripts/reproduce.py` (one command: data hash check, tests, experiment tiers,
  comparison against a results manifest, report rebuild), `scripts/build_report.py` (experiment report
  from report parts), `requirements.txt` with exact versions, and an equation-to-code map.
- Independent reviews, rounds 1, 2a, 2b and 3, and the pre-registration of a live Nifty 50 run.

### Changed

- `pyproject.toml` version set to 0.2.0 (`avmonitor.__version__` still read 0.1.0 until 0.3.0).
- `data.intraday_scale` is strictly predictable by default (`start="past"`: past days only).
- A confidence sequence that has started is never restarted, so each segment uses one sequence per
  side.
- The CLI `run` command loads a warm-up period before `--start` and falls back to the requested window
  when warm-up files are missing.

### Fixed

- MAD truncation and the online scale on zero-heavy and gappy returns: no silent zeros, and
  `VolatilityMonitor` no longer fails when more than half of a calibration window is exactly zero.
- Missing values (NaN) are skipped by the detectors and by the monitor's confidence-sequence calibration
  instead of disabling them, and `calibrate` drops them.
- Running-truncation monitors: a decrease side that is switched off stays at 0 (review round 3).
- Competitors: the BOCPD run-length cap documented correctly and tested where it binds; an exact
  cumulative-sum path for the truncated BJV-type statistic with long windows.

## [0.1.0] - 2026-11-02

Round 1: the first prototype (experiments E1 to E5).

### Added

- Package `avmonitor` with a `src/` layout and tests.
- `avmonitor.robust`: realized quarticity estimators (RQ, MinRQ, MedRQ, tripower, quad-power, bipower,
  threshold RQ), MAD truncation, and kernel-weighted local variance and quarticity in the printed and in a
  normalised form.
- `avmonitor.simulate`: the return model of the reproduced paper and change-point simulators.
- `avmonitor.detectors`: the nonparametric betting e-detector on truncated squared returns
  (Shiryaev-Roberts and CUSUM forms), calibration of the null regime with `margin_z` and `null_band`, and
  benchmarks (Gaussian likelihood ratio e-detector, Page's CUSUM, a window z-test checked at every step).
- `avmonitor.monitor`: `VolatilityMonitor`, the streaming calibrate, monitor, alarm, recalibrate loop
  with a null band, and `GaussianLRMonitor`.
- `avmonitor.data`: a Bitstamp BTC/USD one-minute loader and intraday deseasonalisation.
- Experiments E1 to E5: reproduction of the published tables, truncation bias and a fair comparison of
  estimators, online truncation, detection on simulated data (with null band and calibration length
  studies), and a first real-data run on BTC/USD.
- Research notes (estimator issues, mathematics, literature, experiment report, strategy, iteration log,
  alternatives), a CI workflow that runs the tests on every push, and a synthetic demo script.
