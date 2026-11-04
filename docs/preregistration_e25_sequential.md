# Pre-registration E25: an anytime-valid version of the forward release test

*Written and committed on 3 November 2026, before any Binance return after 24 October 2026 23:59 UTC was downloaded or loaded by anyone working on this project (section 2), and before any code for this test was written. The commit time of this file is the record. Anything decided after this commit will be marked post hoc in the report.*

## 1. Why

E24's forward test (F24a, `docs/preregistration_e24_future.md`) is a fixed-sample test: it may be run once, on 1 June 2027, and its p-value is valid only if nobody acts on the data before then. That has two costs. Nobody may look at the evidence as releases arrive, and if the test fails because the effect is smaller than planned (power 0.37 at half the planned effect), the sample cannot be extended without spending more than its 5% error.

The monitor at the heart of this project solves the same problem for volatility: an e-detector may be checked after every minute and run forever, because its false alarm bound holds at every stopping time (Ville's inequality). E25 applies the same idea one level up, to the release test itself. Each release, once its matched control windows are complete, contributes one e-value; their running product is a test martingale under the null, so the evidence can be checked after every release, the test can stop the first time the evidence reaches 20, and it can continue past 1 June 2027 for as long as anyone wants, with the chance of ever reaching 20 under the null at most 5%.

E25 does not change F24a in any way. F24a runs on 1 June 2027 exactly as pre-registered, and its result is reported whatever E25 shows (section 7).

## 2. Freshness check

Done on 3 November 2026, before this commit, without loading any market data. Since E24's forward pre-registration (commit `22f1988`), no commit loaded a Binance return: the monthly files on this machine (`data/raw/crypto/binance/*/..._2026-10.csv`, downloaded by E23 at 06:37 UTC on 3 November) still end at 24 October 2026 23:59 UTC, and no November 2026 file has been downloaded. So every release in E25's sample, and every one of its control windows, lies in data that nobody in the project has seen.

## 3. Data, monitor and events (shared with F24a)

* **Coins and data:** as in F24a (`docs/preregistration_e24_future.md`, section 4): Binance spot one-minute bars, BTCUSDT (primary), ETHUSDT and SOLUSDT (secondary), loaded from 5 October 2026 (20 days of intraday-profile warm-up, E12's `load_coin` with E23's period switch), **monitored from 25 October 2026 00:00 UTC, with no end date.** Missing returns stay NaN. Data checks as in E23.
* **Monitor:** F24a's primary arm S2 without change: E22's `make_arm("S2")`, the awake rank monitor `awake-B` of E14 (A = 28,800, band 0.25, calibration 1440, seed 14,000), one continuous run per coin from the stream start. The run used at a look covers all complete data up to the look's cut-off. S2 is causal, so its alarms before any time t do not depend on data after t; this is checked by a gate (section 8), and it makes the S2 alarms before 1 June 2027 the same as those of F24a's S2 run.
* **Events:** `data/events/us_macro_events_e24_future.csv`, updated only under F24a's rule (section 3 of that file): FOMC statements and CPI releases with status `released` or `rescheduled` (at their actual time). `cancelled`, `no_statement` and `scheduled` rows are not releases. Rows after 2027 are added under the same rule, from the official schedules, before their month begins.
* **Outcome:** E12's `up_post`: an up alarm in [t, t + 6 h) after the release time t (a release at which the monitor is calibrating counts as one without), computed by E12's `window_stats`.

## 4. Strata: each release with its own control windows

Releases are taken in order of release time (UTC), i = 1, 2, .... For release i at Eastern clock time h on day d:

1. **Candidate controls:** the anchors at clock time h (US Eastern) on the days d + 7k, k in {-4, -3, -2, -1, 1, 2, 3, 4}: the same weekday and the same Eastern clock time, one to four weeks before and after.
2. **Exclusions, as in E12:** a candidate is dropped if its day is less than 2 days from the day of any release in the event file (any FOMC or CPI row with status `released` or `rescheduled`), if its anchor is earlier than 1440 + 360 minutes after the stream start (E12's `CAL + PRE`), or if its window [anchor, anchor + 6 h) is not inside the data.
3. **Disjoint strata:** a candidate already used as a control of an earlier release (smaller i) is dropped. No window is used twice, and no control window is a release window (rule 2).
4. **Completion:** stratum i is complete when the data reach t + 28 days + 6 hours (Eastern day arithmetic, so 4 weeks after the release at the same clock time, plus the window) and every event-file row within 30 days of d has a status other than `scheduled`. Only complete strata enter the test, and they enter in the order i.

Let n_i be the number of controls kept (0 to 8), Z_i the release's outcome, K_i the number of its controls with the outcome, S_i = Z_i + K_i and p_i = S_i / (n_i + 1).

## 5. The e-value of each stratum and the e-process

For an odds ratio ψ > 0, the e-value of stratum i is

    e_i(ψ) = ψ^{Z_i} / (1 - p_i + ψ p_i).

It is the likelihood ratio, for the release's outcome given S_i, of "the release window has odds ψ times those of its controls" against "the release window is one of the n_i + 1 windows chosen at random". If S_i = 0 or S_i = n_i + 1 (no window or every window with the outcome), e_i(ψ) = 1: such a stratum carries no evidence either way. If n_i = 0, e_i(ψ) = 1.

**E-process:** with the grid Ψ = {2, 4, 8, 16} and equal weights,

    M_j = (1/4) Σ_{ψ in Ψ} Π_{i ≤ j} e_i(ψ),   M_0 = 1.

**Decision rule (primary, E25a):** the null hypothesis of section 6 is rejected at the first complete stratum j with **M_j ≥ 20** (= 1/α, α = 0.05). M is reported after every complete stratum, whether or not it crosses.

The grid was fixed from the planning rates of F24a, section 5, before any data: S2's pooled 2024 to 2026 odds ratios are about 17 for FOMC statements (77.3% against 16.5%) and 4.2 for CPI releases (37.5% against 12.5%), and the 2026 rates alone are lower; the grid spans 2 to 16 so that a smaller effect is not ruled out.

## 6. The null hypothesis and the guarantee

**H0 (sequential exchangeability of each stratum):** for every i, given the outcomes of strata 1 to i - 1, the n_i + 1 outcomes of stratum i (the release window and its controls) are exchangeable.

Under H0, P(Z_i = 1 | strata 1 to i - 1, S_i) = p_i, so E[e_i(ψ) | strata 1 to i - 1, S_i] = p_i ψ / (1 - p_i + ψ p_i) + (1 - p_i) / (1 - p_i + ψ p_i) = 1 for every ψ. Each product Π_{i ≤ j} e_i(ψ) is then a nonnegative martingale with mean 1 in j, and so is their average M_j. **Ville's inequality gives P_H0(M_j ≥ 20 for some j ≥ 1) ≤ 1/20 = 5%,** over an unbounded number of releases, whenever anyone looks, and whatever rule they use to stop (the proof and a discussion of the null are in `docs/theory_sequential_events.md`).

Unlike F24a's Poisson-binomial p-value, this guarantee needs no independence between strata, only the conditional statement within each stratum given the earlier ones; the disjoint strata of section 4 are what make that statement meaningful (a window shared by two strata would be known when the later stratum is assessed). Like F24a's null, H0 is a statement about the monitor's outcomes, not about returns, and it is not implied by any theorem about S2: volatility that drifts over a month, or alarms that cluster at certain hours, can break it in either direction. Controls on both sides of the release (rule 1) remove a linear drift to first order.

## 7. Families, F24a, and how results will be read

| Family | Monitor | Cells | Threshold | Role |
|---|---|---|---|---|
| **E25a** | **S2** | **BTC, FOMC and CPI pooled** | **M ≥ 20** | **primary** |
| E25b | S2 | ETH and SOL, pooled, one e-process each | M ≥ 40 each (Bonferroni over the two) | secondary (the coins move together, so they are not independent evidence) |

**Exploratory, reported raw, no decision attached:** separate e-processes for FOMC and for CPI on each coin; the e-process of each single ψ in Ψ; E24's guarded P1 (G1) with the same strata.

**Relation to F24a.** F24a and E25a test related nulls with different controls (F24a: matched windows from the whole period; E25a: up to four weeks either side). Both are pre-registered at 5%. If the release contrast is ever called confirmed on the strength of either test, the report will state that two pre-registered tests were run on the same releases, so the chance that at least one of them rejects under the null is at most 10% (5% + 5%, with no assumption about their dependence). E25 cannot change F24a's verdict on 1 June 2027, and F24a's verdict does not stop E25.

**Reading fixed in advance.**
* **If M reaches 20 (at any complete stratum):** "On FOMC statements and CPI releases after 24 October 2026, which nobody in the project had seen when the design was fixed, the awake rank monitor's up alarms within six hours of a release were more frequent than in matched control windows within four weeks of it, by a sequential test valid at every look (level 5%)." It does not say which release type carries the effect, how large it is (the e-value is not an effect estimate; the release and control rates are reported beside it), or that the monitor's own bound holds on these data.
* **If M has not reached 20 when anyone looks:** M is reported with the number of complete strata, the release and control rates and the largest M so far, and the test continues. A small M is not evidence for H0. Under H0, M stays small; under the planning rates, its expected log growth per release is the quantity reported in the power section (section 9).
* **Nothing about the monitor's guarantee.** As for F24a: a crossing says release windows differ from control windows under this monitor, not that the monitor's 5% per-day bound holds on real days (it does not; E22, E23).

## 8. Gates and refusals (in code)

1. **Identity check:** for every n from 0 to 8, every S from 0 to n + 1 and every ψ in Ψ, the exact conditional mean of e(ψ) under H0 equals 1 to 1e-12.
2. **Causality gate:** on E23's known 2026 data, the S2 alarms of a run cut at an earlier time equal the alarms of the full run before that time (three cut points).
3. **Pipeline gate:** on E23's period and events, the release outcomes Z_i that E25's code computes for S2 equal the `up_post` values in E23's committed S2 windows (`results/e23_windows.csv`), for every coin and release.
4. **Refusals:** a stratum is not computed until it is complete (section 4); the status part refuses to report any stratum whose release or controls lie after the data's last complete minute; the code records the git HEAD and the changes to the shared pipeline since this pre-registration's commit.
5. **Simulated validity:** under H0 with the strata of the planned calendar and control rates from 5% to 40%, the share of 20,000 simulated release sequences on which M ever reaches 20 within 200 releases is at most 5% plus three Monte Carlo standard errors.

## 9. Power (descriptive; computed after this commit; changes nothing above)

Under the planning model of F24a, section 5 (`E23.simulate_power`'s assumptions: independent outcomes with S2's release and control rates; planning calendar from the event file without the statement of 27 January 2027, and 2028 placeholders), with n_i from rule 4's assignment on the planned calendar: the probability that M reaches 20 by 1 June 2027, 1 January 2028 and 1 January 2029, and the median number of releases to the crossing, for four scenarios (pooled 2024 to 2026 rates, 2026 rates, half of the pooled excess, null), seed 25,000, 20,000 replications. For comparison, F24a's power at 1 June 2027 is 0.822, 0.596, 0.369 and 0.031 in the same scenarios. E25 is expected to be less powerful than F24a at the same date (Ville's threshold costs power at a fixed sample, and its controls are fewer), and to gain where F24a cannot follow: after 1 June 2027.

## 10. Post hoc illustration on seen data (exploratory)

The same e-process, computed on releases the project has already seen, illustrates what the procedure does and is not evidence: S2 on BTC, ETH and SOL for the FOMC statements and CPI releases of 2024 to 2025 (E22's S2 run and period) and of January to October 2026 (E23's S2 run and period), with strata built by section 4 within each period. Reported as a path of M with its crossings, labelled post hoc.

## 11. Outputs

`experiments/e25_sequential_events.py` with parts `theory` (gate 1), `gates` (gates 2 and 3, known data), `power` (section 9, synthetic), `posthoc` (section 10, known data), `status` (the fresh test at a look: refuses incomplete strata; writes `results/e25_status_<cut-off date>_*.csv`) and `figs`; `docs/theory_sequential_events.md`; `tests/test_e25.py` (synthetic, no market data); `docs/report_parts/E25.md`. Random numbers: S2's seed 14,000 (as F24a), the power and validity simulations seed 25,000. Runs use `nice -n 10` and at most 2 worker processes.

## Amendment 1 (3 November 2026, about 10:50 UTC, before any fresh data were loaded)

The implementation (commits `d0c2d8f` to `80ee004`) found two gaps in section 4 and one wrong reference. No Binance return after 24 October 2026 23:59 UTC had been downloaded or loaded when this amendment was written, so it is still made before any data of the test. It was written after the power simulation (synthetic) and the post hoc illustration on seen data (section 10) had been run with the original rules; it is motivated by the logic of section 6, not by those results, and both are rerun with the amended rules, the original-rule runs being kept in the report for comparison. The text above is left as committed; where it differs, this amendment governs.

1. **Overlapping windows (rule 3).** Rule 3 forbids a window being used twice but not two windows overlapping in time: a CPI control at 08:30 ET and an FOMC control at 14:00 ET on the same day share the half hour from 14:00 to 14:30 ET (30 such pairs on the planned calendar). An alarm in the shared minutes would then count in two strata, and the outcome of the later stratum's window would be partly known from the earlier one, which is the problem section 6 says disjointness avoids. **Rule 3 now reads:** a candidate is dropped if it was already used as a control of an earlier release, or if its window [anchor, anchor + 6 h) overlaps the window of any release or the window of any control already kept for an earlier release. Release windows are never dropped; two releases whose windows overlap (an FOMC statement and a CPI release on the same day, as on 12 June 2024) both enter, and the overlap is reported. The planned forward calendar has no such pair.
2. **Drift (section 6).** The sentence "Controls on both sides of the release (rule 1) remove a linear drift to first order" holds only for strata whose kept controls are balanced around the release. After rules 2 and 3 most planned strata are unbalanced, and some have controls on one side only (the statement of 28 October 2026 has only later controls). **It now reads:** candidates on both sides of the release (rule 1) remove a linear drift to first order in strata whose kept controls are balanced around the release; rules 2 and 3 leave most strata unbalanced, and the number of balanced strata is reported. A linear drift can therefore move a stratum's e-value in either direction, which is part of what H0 assumes away.
3. **Reference (section 9).** "with n_i from rule 4's assignment" should read "with n_i from the assignment of rules 1 to 3 (section 4)".

## Amendment 2 (3 November 2026, after review round 11, before any fresh data)

Review round 11 (`docs/review/independent_review_round11_e25.md`) found no critical problem. This amendment adds refusals and records corrections; it changes no rule that decides which windows enter a stratum, no e-value and no threshold. No Binance return after 24 October 2026 23:59 UTC had been downloaded or loaded when it was written.

1. **Refusals (finding M1).** Rule 4 checked only that nearby event rows are not `scheduled`, not that they exist, and nothing froze a stratum once reported, so a missing row added later could change a stratum already reported. The status part now also refuses, for the event rows within 30 days of each stratum's windows, unless every calendar month in that range has a CPI row, every scheduled FOMC meeting in that range has a row (as `date` or `original_date`), and no `released` row has an access date before its date (F24a's refusals); (the list of scheduled FOMC meetings is F24a's, which ends in December 2027, so strata whose range reaches into 2028 wait until the 2028 schedule is added to `FOMC_SCHEDULED_LATER` in the script); and it keeps a ledger of every stratum it has reported (`results/e25_status_ledger.csv`) and refuses, writing nothing, if a later look computes any of them differently.
2. **Correction of Amendment 1's time (m6).** Amendment 1 was committed at 10:40 UTC (commit `e00c2b0`), not "about 10:50 UTC" as its heading says.
3. **Clarifications of section 6 and Amendment 1 item 2 (m2, m3).** A drift does not act "in either direction" across strata: every one-sided stratum has its controls after the release (strata are filled in release order and take their later candidates first), so a falling alarm rate favours rejection; in simulation, with no release effect and every window's rate halving every 91 days, M ever reached 20 over the 41 planned strata with probability 0.026 against 0.023 without drift. H0 needs no independence between strata, but it does need the earlier strata to carry no information that separates a release window from its controls; releases on consecutive days (4 pairs to December 2027) are where this could fail. Holiday controls (6 of 95 planned, around Thanksgiving and the year end) are also part of what H0 assumes away.

