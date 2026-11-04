# E20 pre-registration: a re-analysis of the real-data inference

*Written and committed before any adjusted result of E20 was computed: before any Holm adjustment, block bootstrap, forecast regression, loss comparison or event-calendar shift test of this study was run, and before the alarm lists were reloaded or rerun for E20. The commit time of this file is the record.*

## What this is, and what it is not

The data and the alarms analysed here have been seen before. Every alarm list, hit rate, p-value and Wilson interval that the paper reports was computed in earlier experiments (E6, E7, E12, E13, E14, E15, E17) and read by the authors, and by the referee whose report (`docs/review/referee_report_paper.md`, section 4 items 4, 5, 9, 15, 17 and 18, section 4.1 and major comment M6) motivates this study. **E20 is therefore a re-analysis, not a confirmatory study.** What is fixed in advance is only the method of the re-analysis: the families, the adjustment, the resampling scheme, the regression and loss specifications, and the rules that classify each claim. These are fixed before their own results are known. They cannot undo the fact that the questions were chosen by people who had seen the data; the classification of Part 1 records where that happened.

All results of E20 will be **[REAL]** (computed by us on market data). Nothing in this study changes a monitor setting.

## Part 1. Classification of every reported real-data test

**Scope.** Every real-data test or comparison with an inferential statement (a p-value, an interval, "significant", "not significant", "as often", "more often than") in `paper/short_paper.md` Section 6 and Supplement S5, and the real-data sentences of Sections 5.5 and S4.3 that Section 6 relies on (E14 Q14d and Q14g, E17 Q17f and Q17g), and the abstract (line 9). The paper is read as of commit `8f71a2e` (the version the referee read; it has not changed since).

**Evidence used.** The committed pre-registrations (`docs/preregistration*.md`), the docstrings of the experiment scripts, the iteration notes (`docs/iterations/`), and the commit history (`git log`): for each test, the commit that first describes it and the commit that first contains its results.

**Classes.**

* **Pre-specified primary (P).** A document names the test as the primary test, or as a member of the primary family, of its data set, and that document was committed either (i) in a commit strictly earlier than the first commit with the test's results ("verifiable"), or (ii) in the same commit as the results, with an explicit statement that it was fixed before the run ("stated, not verifiable from the commit order").
* **Pre-specified secondary (S).** The test (data, outcome, comparison) is described in a docstring or plan with the same timing evidence as above, but not designated primary; this includes replications named in advance.
* **Exploratory / post hoc (X).** The test was not described before its results; or a document says that it was added or changed after first results (including after a review); or the specific comparison the paper highlights was selected after seeing results (a subgroup, a "best" competitor); or the analysis was designed after the results of the same data set's primary test were known.

**Primary endpoint per data set.** For each data set (five coins in E6; five coins in E13; Nifty 50 monitor validation; Nifty 50 events; Nifty 50 against India VIX; FOMC and CPI on crypto; the panel of E15) the table gives the primary endpoint where one was pre-specified, and says "none was pre-specified" otherwise.

**Output.** `results/e20_classification.csv` and a table in `docs/report_parts/E20.md` with columns: claim, line in the paper, source (experiment and question), class, evidence for the class (commit hashes and the sentence of the docstring or note).

## Part 2. Multiplicity and dependence

### 2.1 Adjustment

Holm's step-down procedure (Holm 1979) within each family below, at a family-wise level of 5%. Adjusted p-values are reported next to the unadjusted bootstrap p-values and next to the original p-values of the experiment's CSV. No adjustment is made across families.

### 2.2 Families

All E6, E13, E14, E15 and E17 families are tests of a difference of two shares of alarms (or of an alarm share and a control share) whose p-value comes from the block bootstrap of 2.3. Band 0.25 unless stated.

| ID | Data | Tests in the family | Size |
|---|---|---|---|
| F1 | E6 Q6d, five coins | e-detector share in the flagged direction (next-day robust volatility against the calibration baseline) minus the matched-delay control share; up, down; Jul to Dec 2025, Jan 2024 to Jun 2025 | 4 |
| F2 | E6 Q6e, five coins | e-detector hit rate against a forecast minus the random-time hit rate; up, down; EWMA, HAR; two periods | 8 |
| F3 | E6 Q6e, five coins | e-detector hit rate minus ratio-rule hit rate; bands 0, 0.25, 0.5; two periods; EWMA, HAR; windows 60 and 240 minutes; up, down | 48 |
| F4 | E6 Q6e, five coins | robust e-detector (band 0) minus Gaussian likelihood-ratio monitor (band 0), EWMA forecast; up, down; two periods | 4 |
| F5 | E6 Q6d cross-venue | the existing circular-shift p-values (same-direction alarm within 60 minutes), four detectors, both directions of matching; Holm only (the shift test already keeps the serial dependence of each venue's alarms) | 8 |
| F6 | E13 crypto, first-run windows | e-detector minus competitor; four competitors (BOCPD, truncated BOCPD, BJV, studentised BJV); up, down; Jul to Dec 2025 and Jan to Jun 2025 | 16 |
| F7 | E13 crypto, longer windows | the same for the eight longer-window competitors | 32 |
| F8a | E13 Nifty, first-run windows | e-detector minus competitor, four competitors, up, down | 8 |
| F8b | E13 Nifty, longer windows | the same for the four session-long competitors | 8 |
| F9a | E7 Q7b, Nifty, primary configuration (calibration 750, band 0.25) | alarm share in the flagged direction minus the random-time control share; up, down; development, held-out | 4 |
| F9b | E7 Q7b | robust (primary configuration) minus Gaussian likelihood-ratio monitor (calibration 750); up, down; two periods | 4 |
| F10a | E7 Q7d events | Fisher p-values (as published) for "any alarm" on Union Budget, RBI policy, election result, March 2020 window and expiry sessions against ordinary sessions | 5 |
| F10b | E7 Q7d events | the same five event types with the circular-shift test of 2.4 | 5 |
| F10c | E7 Q7d events | Fisher p-values restricted to sessions monitoring at the open (as published) | 5 |
| F11 | E7 Q7c India VIX | Fisher p-values (as published), primary configuration: same day and next day; up and down; development and held-out | 8 |
| F12a | E12 | the pre-registered primary family, reported as published: BTC, streaming, up alarm within 6 hours, FOMC and CPI | 2 |
| F12b | E12 | replications, streaming: ETH and SOL, FOMC and CPI | 4 |
| F12c | E12 | fresh calibration, robust monitor: three coins, FOMC and CPI | 6 |
| F12d | E12 | streaming, releases at which the monitor was monitoring at the anchor: three coins, FOMC and CPI | 6 |
| F13 | E14 Q14g | stratified p-values (as published) of the three restarted or rank monitors (B-lr, awake-A, awake-B); three coins; FOMC and CPI | 18 |
| F14 | E14 Q14d, five coins | awake-A minus plug-in and awake-B minus plug-in, share of up (down) alarms followed by higher (lower) next-day robust volatility than the baseline; up, down | 4 |
| F15 | E15, S5.3 | e-BH minus the rate-matched per-asset threshold, share of discoveries confirmed; the 12 cells of `results/e15d_rate_matched_tests.csv` | 12 |
| F16 | E17 Q17g, five coins | matched one-minute monitor minus coarse detector (5 or 15 minutes); confirmed share, EWMA hit rate, HAR hit rate; up, down | 12 |
| F17 | E20 Part 3 | primary regression family (3.6) | 6 |
| F18 | E20 Part 3 | primary loss family (3.7) | 3 |

For F10a, F10c, F11, F12 and F13 the original test (Fisher exact test or E12's exact stratified permutation test) is kept and only the Holm adjustment is added; F10b adds a dependence-robust version for the Nifty events.

### 2.3 Block bootstrap over calendar days

* **Unit.** A calendar day (UTC) for crypto, a trading session for the Nifty 50. Every alarm, discovery, random time and control draw belongs to the day (session) of its own time stamp. The days of a period are all calendar days of that period, including days without alarms: Jul to Dec 2025, 184 days; 21 Jan 2024 to 30 Jun 2025, 527 days; Jan to Jun 2025, 181 days; the Nifty sessions of each period in E7's cleaned calendar (2024 and 2025: 494 sessions for E13). For E15, the days from the first to the last day with any discovery in that panel and period.
* **Scheme.** Circular block bootstrap (Politis and Romano 1992) with blocks of 7 consecutive days (crypto) or 5 consecutive sessions (Nifty), B = 9,999 replicates, seed 20. A replicate draws ceil(D / L) block starts uniformly from the D days, concatenates the blocks circularly and keeps the first D days; the weight of a day is the number of times it was drawn. **Days are resampled jointly**: one set of weights per replicate is applied to every coin, every method, both directions and the control draws of that period, so the dependence across coins and between methods on the same day is kept.
* **Statistics.** A share is (weighted hits) / (weighted scored alarms). A difference is the difference of two such shares computed with the same weights.
* **Intervals.** Percentile 95% intervals of the replicate statistics.
* **p-values.** For H0: difference = 0, p = min(1, 2 min(1 + #{Δ* ≤ 0}, 1 + #{Δ* ≥ 0}) / (B + 1)), the inversion of the percentile interval. The bootstrap standard error is also reported.
* **Controls that are not stored per draw.** E7's random-time control shares (`results/e07b_monitor_summary.csv`) are treated as fixed numbers (their own sampling error is about a third of the alarms', which this understates); the test is then share* minus the fixed control. E6's random times and matched-delay control draws are regenerated per draw (2.5) and resampled with the alarms.
* **Sensitivity (secondary).** Blocks of 14 days and 10 sessions.
* **Scored alarms.** Alarms without a scored outcome are dropped, as in the original experiments.

### 2.4 Circular-shift test for the Nifty events (F10b)

The 2,361 cleaned sessions in order; the indicator "the session has at least one alarm" of the primary configuration (from `results/e07b_alarms.csv`); for each event type the set of its sessions. Statistic T = number of event sessions with an alarm. Under the null that the alarm series is unrelated to the event calendar, shifting the calendar by k sessions (circularly, k = 0, ..., N - 1, all shifts) gives the reference distribution while keeping the serial dependence of the alarms. Two-sided p = #{k : |T_k - mean T| >= |T_0 - mean T|} / N.

### 2.5 Equivalence statements

Where the paper claims that two rules are "about as" precise or "as informative", E20 calls equivalence **shown** if the 95% bootstrap interval of the difference lies within plus or minus 5 percentage points, and **not shown** otherwise. The margin of 5 points is E20's choice, made here.

## Part 3. Information beyond a forecast (referee items 4 and M6)

### 3.1 Data and alarms

Five Binance coins (BTC, ETH, XRP, SOL, DOGE against USDT), one-minute closes from 1 January 2024 to 31 December 2025, raw log returns r_t (not deseasonalised, not truncated), loaded with `avmonitor.data.load_binance(..., download=False)`. Two periods, as in E6: **development** = 1 July to 31 December 2025, **held-out** = 21 January 2024 to 30 June 2025. Alarm lists (3.2): the E6 e-detector at band 0.25 (primary) and the two E6 ratio rules (windows S = 60 and 240 minutes) matched to it; band 0.5 and its matched ratio rules as a secondary block. Each period uses the alarms of that period's run.

### 3.2 Alarm lists and checks (task item 4)

* e-detector: `results/e06d_events.csv`. Check: counts per coin, period, band and direction equal `results/e06d_monitor_summary.csv`.
* Ratio rules: not stored per alarm, so rerun with the functions of `experiments/e06b_five_coins_monitor.py` (imported, not modified): `prepare`, `Series`, `match_rule`, `ratio_rule`. Check: thresholds k_up and k_down equal `results/e06e_competitor.csv` to a relative 1e-9 and the up and down counts equal it exactly. In addition, `summarise` and `forecast_hits` are called in E6's job order with E6's seed (6), and their output must reproduce `results/e06d_monitor_summary.csv` and `results/e06e_forecast_hits.csv` (counts exactly, shares to 1e-12); the random times and matched-delay control draws of E6 are regenerated from the same generator state and checked against those files.
* E13: `results/e13d_real_long_windows_crypto_alarms.csv` and `..._nifty_alarms.csv` (per-alarm outcomes stored). Check: pooled shares equal `results/e13d_real_long_windows_pooled.csv`.
* E7: `results/e07b_alarms.csv`. Check: counts and shares equal `results/e07b_monitor_summary.csv`; event-session counts equal `results/e07d_event_alarm_rates.csv`.
* E14: `results/e14d_real_alarms.csv`. Check: counts equal `results/e14d_real_summary.csv`.
* E15: `results/e15d_discoveries.csv.gz`. Check: the 12 cells equal `results/e15d_rate_matched_tests.csv`.
* E17: coarse detectors from `results/e17_real_events.csv`; the matched one-minute monitor is rerun once per coin and target with the stored levels of `results/e17_ratematch_levels.csv` and E17's `_one_minute_events` (imported, not modified). Check: counts equal `matched_up` and `matched_down`, and pooled shares equal `results/e17_ratematch_pooled.csv`.

If a check fails, that list is not used, and the failure is reported.

### 3.3 Outcomes

At an evaluation time t (information up to and including the return of minute t), with H = 1440 minutes (primary) or 360 minutes (secondary):

* realised variance RV_H(t) = sum of r_s^2 over s = t + 1, ..., t + H;
* bipower variation BV_H(t) = (pi / 2) (H / (H - 1)) sum of |r_s| |r_{s-1}| over s = t + 2, ..., t + H (Barndorff-Nielsen and Shephard 2004), the jump-robust alternative (secondary).

The outcome is log RV_H(t) (or log BV_H(t)). Times whose window has a missing return or a zero RV or BV are dropped.

### 3.4 Evaluation grid and alarm indicators

Grid: every full UTC hour of the period (t at minute 0 of each hour) such that the outcome window ends inside the period and every forecast in 3.5 is defined at t (one common sample for all forecasts). For a rule k, U_k(t) = 1 if rule k raised an up alarm at a minute in (t - 60, t], and D_k(t) = 1 if it raised a down alarm there; each alarm therefore enters exactly once, at most 59 minutes before its grid point.

### 3.5 Forecasts

Regressors at t, from returns up to t only:

* E(t): EWMA of r^2 with a half-life of 1440 minutes (E6's half-life), `pandas` `ewm(halflife=1440, ignore_na=True)`;
* RVd(t) = sum of r^2 over the 1440 minutes ending at t; RVw(t) = mean of RVd(t - 1440 j), j = 0, ..., 4; RVm(t) = mean of RVd(t - 1440 j), j = 0, ..., 21;
* RVh(t) = sum of r^2 over the 60 minutes ending at t (the intraday component), floored at 1% of RVd(t) / 24;
* calendar: dummies for the UTC hour of t (23) and the day of the week of t (6).

Three forecast models of y = log RV_H(t) (or log BV_H(t)), each an OLS regression on its regressors and the calendar dummies:

* **EWMA**: 1, log E(t);
* **HAR** (Corsi 2009, on rolling 1, 5 and 22 day sums anchored at t): 1, log RVd, log RVw, log RVm;
* **HAR-I**: HAR plus log RVh.

Each model is refitted for each coin, horizon and outcome measure once per UTC day d, on all hourly grid times of that coin from 1 January 2024 whose outcome window ended before d starts (expanding window), with at least 720 training points; the forecast for times in day d uses day d's coefficients. log F(t) is the fitted value; the level forecast is F(t) = exp(log F(t) + s^2 / 2), with s^2 the residual variance of that fit.

**The strongest forecast** is the model with the smallest mean QLIKE (3.7) of next-day RV on the development period, pooled over the five coins.

### 3.6 Forecast regression

For each period, horizon, outcome measure, forecast model f and rule k:

y(t) = a_coin + b log F_f(t) + c_up U_k(t) + c_dn D_k(t) + e(t),

pooled over coins and grid times, by OLS with coin intercepts. Standard errors: Driscoll and Kraay (1998), scores summed over all coins within each UTC calendar day, Bartlett weights with a bandwidth of 5 days (Newey and West 1987), multiplied by D / (D - 1) with D the number of days; normal reference distribution; two-sided p-values for c_up = 0 and c_dn = 0. Secondary standard errors: clustered by calendar day only (bandwidth 0). A rule adds information beyond the forecast in a direction if the coefficient has the sign of the alarm (c_up > 0, c_dn < 0) and its Holm-adjusted p-value is below 0.05.

* **Primary family F17**: held-out period, next-day log RV, the strongest forecast, rules {e-detector, ratio rule 60, ratio rule 240} at band 0.25, coefficients {c_up, c_dn}: 6 tests.
* **Secondary**: every other combination (development period, 6-hour horizon, bipower outcome, the other two forecasts, band 0.5), reported with nominal p-values and Holm within each block of 6 coefficients; and a joint model with the indicators of the e-detector and of one ratio rule together (does each add information given the other).

### 3.7 Loss comparison (QLIKE)

QLIKE loss L(y, F) = y / F - log(y / F) - 1 (Patton 2011), with y = RV_H(t) in levels. For rule k and forecast f, the adjusted forecast is F'(t) = F(t) kappa_up^{U_k(t)} kappa_dn^{D_k(t)}, where kappa_up is the mean of y / F over development-period grid times with U_k = 1 (pooled over coins; this is the QLIKE-minimising constant), and kappa_dn likewise. On the held-out period: the mean over all grid times of L(y, F') - L(y, F) (negative means the adjustment helps), its Driscoll and Kraay standard error as in 3.6 (a Diebold and Mariano 1995 type test of equal predictive accuracy), the two-sided p-value, the same mean over the grid times with an alarm only, and the relative change in mean QLIKE.

* **Primary family F18**: next-day RV, strongest forecast, rules {e-detector, ratio rule 60, ratio rule 240} at band 0.25: 3 tests.
* **Secondary**: the other forecasts, the 6-hour horizon, the reverse direction (kappas fitted on the held-out period, evaluated on the development period).

## Part 4. How each claim will be classified after adjustment

For each real-data claim of the paper, `docs/report_parts/E20.md` gives its adjusted status:

* **holds**: the claimed difference has an adjusted p-value below 0.05 within its family, and the bootstrap interval excludes 0 on the claimed side (for Part 3: the conditions of 3.6 or 3.7);
* **holds only unadjusted**: nominal p below 0.05, adjusted p at or above 0.05;
* **does not hold**: nominal p at or above 0.05, or the interval contains 0;
* **equivalence shown / not shown**: for "as often" claims, by the rule of 2.5;
* **exploratory**: claims of class X keep whatever status the numbers give, but are labelled exploratory and carry no confirmatory weight, whatever their p-value;

and a recommended replacement wording, which the authors may apply. E20 does not edit the paper.

## Files

`experiments/e20_inference.py` (parts `checks`, `bootstrap`, `holm`, `events`, `forecast`, `classify`), `tests/test_e20_inference.py` (helpers: Holm, circular block weights, bootstrap p-value, Driscoll and Kraay standard errors, QLIKE constant), `results/e20_*.csv`, `docs/report_parts/E20.md`, `docs/iterations/E20.md`. Runs use at most 2 worker processes under `nice -n 10`.

## References added by E20

Checked by web search on 3 November 2026 (publisher pages and arxiv.org cannot be opened from this machine; the details below come from index and library listings).

* Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics* 6(2), 65 to 70.
* Politis, D. N. and Romano, J. P. (1992). The circular block bootstrap, in R. LePage and L. Billard (eds.), *Exploring the Limits of Bootstrap*, Wiley, New York. The volume and the attribution were confirmed; the exact chapter title and pages could not be verified.
* Driscoll, J. C. and Kraay, A. C. (1998). Consistent covariance matrix estimation with spatially dependent panel data. *Review of Economics and Statistics* 80(4), 549 to 560.
* Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica* 55(3), 703 to 708.
* Patton, A. J. (2011). Volatility forecast comparison using imperfect volatility proxies. *Journal of Econometrics* 160(1), 246 to 256.
* Barndorff-Nielsen, O. E. and Shephard, N. (2004). Power and bipower variation with stochastic volatility and jumps. *Journal of Financial Econometrics* 2(1), 1 to 37.
* Diebold, F. X. and Mariano, R. S. (1995). Comparing predictive accuracy. *Journal of Business and Economic Statistics* 13(3), 253 to 263.
* Corsi, F. (2009). A simple approximate long-memory model of realized volatility. *Journal of Financial Econometrics* 7(2), 174 to 196 (already cited by the paper).
