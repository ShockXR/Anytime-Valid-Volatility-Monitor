"""E25. An anytime-valid version of the forward release test: an e-process over releases.

Pre-registration: docs/preregistration_e25_sequential.md (commit 30804bb, 3 November 2026, before any Binance return
after 24 October 2026 23:59 UTC was loaded and before this code was written). Proof: docs/theory_sequential_events.md.
Report part: docs/report_parts/E25.md. Iteration notes: docs/iterations/E25.md.

Each release i (FOMC statement or CPI release, in order of release time) gets a stratum: its own 6-hour window and up
to eight control windows at the same US Eastern weekday and clock time, one to four weeks before and after it, chosen
from the calendar only (pre-registration, section 4). With Z_i the release window's outcome (an up alarm of the rank
monitor S2 within 6 hours), n_i the number of controls kept, S_i the number of the n_i + 1 windows with the outcome and
p_i = S_i / (n_i + 1), the stratum's e-value for the odds ratio psi is e_i(psi) = psi^Z_i / (1 - p_i + psi p_i), and

    M_j = (1/4) sum over psi in {2, 4, 8, 16} of prod_{i <= j} e_i(psi),   M_0 = 1.

Under H0 (the n_i + 1 outcomes of each stratum are exchangeable given the earlier strata) M is a nonnegative martingale
with mean 1, so P(M_j >= 20 for some j) <= 5% (Ville's inequality), at every look and with no end date.

Amendment 1 (end of the pre-registration, commit e00c2b0, 3 November 2026, before any fresh data; made after the
power simulation and the post hoc run on seen data, motivated by section 6): rule 3 also drops a candidate whose
window [anchor, anchor + 6 h) overlaps the window of any release or of any control already kept for an earlier
release (release windows are never dropped; overlapping release windows are reported), and the number of balanced
strata is reported. These amended rules are the default of every part, and the only rules of `status`;
--rules=original reruns strata, power, posthoc and theory with the rules as first committed and writes
results/e25_*_original_rules.csv (no figures), kept for comparison.

Parts (each writes its own results/e25_*.csv):
  theory   [PROVED identity checked exactly; SYNTHETIC] gate 1: the exact conditional mean of e(psi) under H0 is 1 for
           n = 0..8, S = 0..n+1 and psi in the grid (rational arithmetic and floating point, results/e25_theory_identity.csv);
           gate 5: simulated validity under H0 with the strata sizes of the planned calendar, 20,000 sequences of 200
           releases, seed 25,000, exchangeable Bernoulli strata at four rates and four other laws of S
           (results/e25_validity.csv).
  strata   [calendar only, no market data] the deterministic assignment of section 4 (rules 1 to 3 with Amendment 1;
           completion, rule 4) on the planned forward calendar, with the balance of each stratum's controls
           (results/e25_strata_planned.csv).
  gates    [REAL, known data only, no download] gate 2 (causality: S2 on E23's 2026 BTC stream, full run against runs
           whose data end at three earlier cut-offs) and gate 3 (pipeline: E25's release outcomes on E23's period equal
           E23's committed up_post for S2, every coin and release; E25's own S2 runs reproduce the committed S2 alarm
           files of E22 and E23) (results/e25_gates.csv).
  power    [SYNTHETIC; planning inputs REAL rates of 2024 to 2026] section 9: P(M reaches 20 by 1 June 2027, 1 January
           2028, 1 January 2029), releases to the crossing, expected log growth, and an oracle that knows the control
           rate (results/e25_power.csv, e25_power_curve.csv, e25_power_calendar.csv, e25_power_growth.csv; figure
           results/figures/e25_power.png).
  posthoc  [REAL, post hoc, exploratory] section 10: the e-process on E22's period (2024 to 2025) and E23's period
           (2026), S2 on BTC, ETH and SOL, outcomes from the committed S2 alarm files (gate 3 shows they are E25's own);
           also the control rates before and after the release, M over the balanced strata alone, and the missing
           minutes inside each window (read from the files on disk, no download)
           (results/e25_posthoc_strata.csv, e25_posthoc_summary.csv; figure results/figures/e25_posthoc.png).
  status   the fresh test at a look (sections 3, 4, 7, 8.4 and Amendment 2). Refuses (SystemExit with the reasons)
           before any download when no stratum is complete at the cut-off (including, by Amendment 2, while the event
           file is incomplete within 30 days of a stratum's windows); downloads the Binance monthly files as E24's
           future part does; refuses again with the data's actual end; runs S2 (and G1, exploratory) from 25 October
           2026 (loaded from 5 October); refuses, writing nothing, when a stratum already in results/e25_status_ledger.csv
           is recomputed differently; then writes only complete strata (results/e25_status_<YYYY-MM-DD_HHMM of the
           cut-off>_*.csv, with the missing minutes of each window), the provenance file and the new ledger rows. NOT
           RUN in this commit: the first stratum (the FOMC statement of 28 October 2026) completes on 26 November 2026
           at 01:00 UTC.
  dryrun   [REAL, known data] gate 4's dry run: the status part on E23's known period (data to 24 October 2026), no
           download, two looks, outputs in a temporary directory; summaries in results/e25_dryrun_known_data_summary.csv.
  figs     results/figures/e25_*.png from the CSVs.

Usage: nice -n 10 python experiments/e25_sequential_events.py <part> [--procs=2] [--cutoff=YYYY-MM-DD]
       [--rules=amended|original]
"""

from __future__ import annotations

import itertools
import subprocess
import sys
import time
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

import e12_crypto_events as E12  # noqa: E402
import e22_guaranteed_events as E22  # noqa: E402
import e23_fresh_releases as E23  # noqa: E402
import e24_blind_guard as E24  # noqa: E402

OUT = ROOT / "results"
FIG = OUT / "figures"

# ----------------------------------------------------------------------------
# Fixed by the pre-registration (docs/preregistration_e25_sequential.md); nothing here may change after 30804bb
# ----------------------------------------------------------------------------
PREREG_COMMIT = "30804bb"
PSI = (2.0, 4.0, 8.0, 16.0)  # section 5, equal weights
ALPHA = 0.05
THRESHOLD = 1.0 / ALPHA  # 20: E25a, BTC
THRESHOLD_E25B = 2.0 / ALPHA  # 40: E25b, ETH and SOL each (Bonferroni over the two)
WEEKS = (-4, -3, -2, -1, 1, 2, 3, 4)  # rule 1
EXCLUDE_DAYS = E12.EXCLUDE_DAYS  # rule 2: a candidate day less than 2 days from a release day is dropped
START_GAP_MIN = E12.CAL + E12.PRE  # rule 2: anchor at least 1440 + 360 minutes after the stream start
WINDOW_MIN = E12.POST  # the outcome window [anchor, anchor + 6 h)
COMPLETE_DAYS = 28  # rule 4: data to t + 28 days + 6 hours (Eastern day arithmetic)
STATUS_DAYS = 30  # rule 4: no event-file row within 30 days of the release day may still be 'scheduled'
RELEASED = ("released", "rescheduled")
# Amendment 1 (docs/preregistration_e25_sequential.md, commit e00c2b0, before any fresh data): rule 3 also drops a
# candidate whose window overlaps the window of any release or of any control kept for an earlier release. The
# amended rules are the default everywhere (the status part uses nothing else); "original" reproduces the rules as
# first committed (30804bb), for the comparison in the report only, and its outputs carry the suffix below.
RULES = "amended"
ORIGINAL_SUFFIX = "_original_rules"


def _suffix() -> str:
    return "" if RULES == "amended" else ORIGINAL_SUFFIX
FUTURE_EVENTS = E24.FUTURE_EVENTS
STREAM_START = E24.FUTURE_MONITOR_START  # 2026-10-25 00:00 UTC
LOAD_START = E24.FUTURE_LOAD_START  # 2026-10-05: 20 days of intraday-profile warm-up (known data)
PRIMARY_COIN = "BTC"
PRIMARY_ARM = "S2"
STATUS_ARMS = ("S2", "G1")  # G1: exploratory (section 7)
PIPELINE_FILES = ("experiments/e12_crypto_events.py", "experiments/e14_running_tau.py",
                  "experiments/e22_guaranteed_events.py", "experiments/e23_fresh_releases.py",
                  "experiments/e24_blind_guard.py", "src/avmonitor", "data/events/us_macro_events_e24_future.csv")
SEED = 25_000
REPS = 20_000
VALIDITY_RELEASES = 200
VALIDITY_RATES = (0.05, 0.10, 0.20, 0.40)
POWER_DATES = ("2027-06-01", "2028-01-01", "2029-01-01")
PLANNING_END = "2029-01-01"
F24A_POWER_AT = "2027-06-01"
CAUSALITY_CUTS = ("2026-04-01", "2026-07-01", "2026-08-12 14:37", "2026-09-15")  # mid-day cut: review round 11, m8
# Amendment 2 (docs/preregistration_e25_sequential.md, after review round 11, before any fresh data), additional
# refusals of the status part only (they change no rule, e-value or threshold): M1 (i) the event file must be complete around each stratum, as F24a's refusals require
# for its own period: FOMC_SCHEDULED is F24a's list of scheduled meetings, known to the end of the last listed year;
# M1 (ii) strata already reported are frozen in LEDGER and a later look that recomputes any of them differently refuses.
FOMC_SCHEDULED = E24.FUTURE_FOMC_SCHEDULED
INTEGRITY_DAYS = COMPLETE_DAYS + STATUS_DAYS  # 30 days of the windows, which lie within 28 days of the release day
LEDGER = "e25_status_ledger.csv"
LEDGER_EXACT = ("date", "event_type", "release_utc", "controls", "n", "Z", "K", "S")
LEDGER_FLOAT = tuple(f"e_psi{int(p)}" for p in PSI) + ("M",)
SCENARIOS = ("pooled 2024 to 2026 (planning)", "2026 only (E23 design)", "half of the pooled excess",
             "null (release rate = control rate)")


# ----------------------------------------------------------------------------
# The e-value of a stratum and the e-process (section 5)
# ----------------------------------------------------------------------------


def e_value(z, s, n, psi: float):
    """e_i(psi) = psi^Z / (1 - p + psi p), p = S / (n + 1); vectorised. Equals 1 when S = 0, S = n + 1 or n = 0."""
    z, s, n = np.asarray(z, float), np.asarray(s, float), np.asarray(n, float)
    p = s / (n + 1.0)
    return psi**z / (1.0 - p + psi * p)


def e_value_oracle(z, p0, psi: float):
    """The oracle's e-value psi^Z / (1 - p0 + psi p0), with p0 the true control rate (descriptive only)."""
    z, p0 = np.asarray(z, float), np.asarray(p0, float)
    return psi**z / (1.0 - p0 + psi * p0)


def e_process(e_by_psi: np.ndarray, weights=None) -> tuple[np.ndarray, np.ndarray]:
    """From e-values of shape (..., J, len(PSI)): the running products (same shape) and M_j = sum_psi w_psi prod (...,
    J). Equal weights by default."""
    e_by_psi = np.asarray(e_by_psi, float)
    g = e_by_psi.shape[-1]
    w = np.full(g, 1.0 / g) if weights is None else np.asarray(weights, float)
    prod = np.cumprod(e_by_psi, axis=-2)
    return prod, prod @ w


def stratum_e(z, s, n, psi=PSI) -> np.ndarray:
    return np.stack([e_value(z, s, n, ps) for ps in psi], axis=-1)


# ----------------------------------------------------------------------------
# Strata (section 4): a pure function of the release list, the release days, the stream start and the data end
# ----------------------------------------------------------------------------


def et_anchor(day, time_et: str) -> pd.Timestamp:
    """The UTC time of the US Eastern clock time ``time_et`` on the Eastern calendar day ``day`` (DST by the IANA
    database, as E12)."""
    d = pd.Timestamp(day)
    if d.tzinfo is not None:
        d = d.tz_localize(None)
    return pd.Timestamp(f"{d:%Y-%m-%d} {time_et}").tz_localize(E12.TZ, ambiguous="raise",
                                                                 nonexistent="raise").tz_convert("UTC")


def _utc(t) -> pd.Timestamp:
    t = pd.Timestamp(t)
    return t.tz_localize("UTC") if t.tzinfo is None else t.tz_convert("UTC")


def _day(d) -> pd.Timestamp:
    d = pd.Timestamp(d)
    return (d.tz_localize(None) if d.tzinfo is not None else d).normalize()


def releases_frame(releases: pd.DataFrame) -> pd.DataFrame:
    """Releases in order of release time (UTC), numbered 1, 2, ... (the strata). Needs the columns date (US Eastern
    calendar day), event_type and release_time_et."""
    r = releases[["date", "event_type", "release_time_et"]].copy()
    if not len(r):
        return pd.DataFrame(columns=["stratum", "date", "event_type", "release_time_et", "release_utc", "day"])
    r["day"] = pd.to_datetime([_day(d) for d in r["date"]])
    r["release_utc"] = [et_anchor(d, h) for d, h in zip(r["day"], r["release_time_et"])]
    r = r.sort_values("release_utc", kind="stable").reset_index(drop=True)
    if r["release_utc"].duplicated().any():
        raise ValueError("two releases at the same time")
    r["stratum"] = np.arange(1, len(r) + 1)
    r["date"] = r["day"].dt.strftime("%Y-%m-%d")
    return r[["stratum", "date", "event_type", "release_time_et", "release_utc", "day"]]


def assign_strata(releases: pd.DataFrame, event_days, stream_start, data_end=None, rules: str | None = None) -> pd.DataFrame:
    """Rules 1 to 3 of the pre-registration, section 4, with Amendment 1. One row per (release, candidate control).

    Rule 1: candidates at the release's US Eastern clock time on the Eastern days d + 7k, k in +-1..4.
    Rule 2: a candidate is dropped if its day is less than 2 days from any day in ``event_days`` (the release days of
    the event file), if its anchor is earlier than stream start + 1440 + 360 minutes, or if its window
    [anchor, anchor + 6 h) does not end by ``data_end`` (exclusive; None: no data end, for planning).
    Rule 3: a candidate kept by an earlier release (smaller release time) is dropped; with the amended rules
    (default) also a candidate whose window overlaps the window of any release in ``releases`` (earlier or later)
    or of any control already kept for an earlier release. Release windows are never dropped. When every release
    day is in ``event_days``, rule 2 already removes every candidate that overlaps a release window (such a candidate
    is on the release's day or the day next to it); the release part of rule 3 then never binds.
    ``rules="original"`` gives rule 3 as first committed. Uses no market data and no outcome."""
    rules = RULES if rules is None else rules
    if rules not in ("amended", "original"):
        raise ValueError(rules)
    rel = releases_frame(releases)
    rel_ns = _ns(rel["release_utc"]) if len(rel) else np.zeros(0, np.int64)
    span_ns = WINDOW_MIN * 60 * 10**9
    days = sorted({_day(d) for d in event_days})
    days_ns = np.array([d.value for d in days], dtype=np.int64)
    min_anchor = _utc(stream_start) + pd.Timedelta(minutes=START_GAP_MIN)
    end = None if data_end is None else _utc(data_end)
    win = pd.Timedelta(minutes=WINDOW_MIN)
    used: dict[pd.Timestamp, int] = {}
    rows = []
    for r in rel.itertuples(index=False):
        for k in WEEKS:
            cday = r.day + pd.Timedelta(days=7 * k)
            a = et_anchor(cday, r.release_time_et)
            gap = np.abs(days_ns - cday.value) // (86_400 * 10**9) if days_ns.size else np.array([])
            near = [days[j] for j in np.flatnonzero(gap < EXCLUDE_DAYS)]
            if near:
                reason = f"rule 2: within {EXCLUDE_DAYS} days of the release day {near[0]:%Y-%m-%d}"
            elif a < min_anchor:
                reason = "rule 2: anchor earlier than the stream start + 1800 minutes"
            elif end is not None and a + win > end:
                reason = "rule 2: window not inside the data"
            elif a in used:
                reason = f"rule 3: used by stratum {used[a]}"
            elif rules == "amended" and (hit := np.flatnonzero(np.abs(rel_ns - a.value) < span_ns)).size:
                reason = f"rule 3 (Amendment 1): overlaps the release window of stratum {int(rel.stratum.iloc[hit[0]])}"
            elif rules == "amended" and (hit := [j for b, j in used.items() if abs(b.value - a.value) < span_ns]):
                reason = f"rule 3 (Amendment 1): overlaps a control of stratum {hit[0]}"
            else:
                reason = "kept"
                used[a] = int(r.stratum)
            rows.append(dict(stratum=int(r.stratum), date=r.date, event_type=r.event_type,
                             release_time_et=r.release_time_et, release_utc=r.release_utc, k=k,
                             candidate_day=f"{cday:%Y-%m-%d}", anchor_utc=a, anchor_eastern_zone=a.tz_convert(E12.TZ).tzname(),
                             kept=reason == "kept", reason=reason))
    cols = ["stratum", "date", "event_type", "release_time_et", "release_utc", "k", "candidate_day", "anchor_utc",
            "anchor_eastern_zone", "kept", "reason"]
    return pd.DataFrame(rows, columns=cols)


def strata_table(cands: pd.DataFrame) -> pd.DataFrame:
    """One row per stratum: the release, its kept controls (n of them) and the time the data must reach for the
    stratum to be complete, t + 28 days + 6 hours in Eastern day arithmetic (rule 4). ``balanced`` (Amendment 1,
    item 2): as many kept controls before the release as after it, and at least one on each side; a stratum without
    controls (e = 1 whatever happens) is neither balanced nor one-sided."""
    rows = []
    span = pd.Timedelta(minutes=WINDOW_MIN)
    rel_t = cands.groupby("stratum", sort=True).release_utc.first()
    for s, g in cands.groupby("stratum", sort=True):
        r = g.iloc[0]
        kept = g[g.kept].sort_values("anchor_utc")
        day = _day(r.date)
        before, after = int((kept.k < 0).sum()), int((kept.k > 0).sum())
        other = [f"stratum {int(j)} ({t:%Y-%m-%d %H:%M})" for j, t in rel_t.items()
                 if j != s and abs(t - r.release_utc) < span]
        rows.append(dict(stratum=int(s), date=r.date, event_type=r.event_type, release_time_et=r.release_time_et,
                         release_utc=r.release_utc, n=int(len(kept)),
                         controls=";".join(f"{a:%Y-%m-%d %H:%M}" for a in kept.anchor_utc),
                         controls_before=before, controls_after=after, balanced=before == after > 0,
                         release_window_overlaps=";".join(other),
                         data_needed_until=et_anchor(day + pd.Timedelta(days=COMPLETE_DAYS), r.release_time_et)
                         + pd.Timedelta(minutes=WINDOW_MIN)))
    cols = ["stratum", "date", "event_type", "release_time_et", "release_utc", "n", "controls", "controls_before",
            "controls_after", "balanced", "release_window_overlaps", "data_needed_until"]
    return pd.DataFrame(rows, columns=cols)


def _one_sided(tab: pd.DataFrame) -> pd.Series:
    """Strata with kept controls on one side of the release only (n > 0)."""
    return (tab.controls_before == 0) != (tab.controls_after == 0)


def calendar_integrity(table: pd.DataFrame, event_rows: pd.DataFrame, stream_start, fomc_scheduled=None) -> list[str]:
    """Amendment 2, item 1 (review round 11, M1; F24a's refusals, applied around each stratum). For each stratum of
    ``table``, the reasons the event file cannot yet be trusted within 30 days of the stratum's windows (an empty
    string when it can). The candidate windows of release day d lie from d - 28 to d + 28 days, so the range is the
    Eastern days [d - 58, d + 58]: every calendar month that meets it, from the month after the stream start's, must
    have a CPI row (as date or original_date, any status); every scheduled FOMC meeting of ``fomc_scheduled`` in it
    must have a row (as date or original_date; rows are never deleted), and the schedule must be known to its end (to
    the end of the last listed year); no 'released' row in it may have an access date before its date. A missing row
    would otherwise neither exclude nearby candidates (rule 2) nor block completion (rule 4), and a missing release
    within 56 days could later take candidates through Amendment 1's overlap rule."""
    sched = tuple(FOMC_SCHEDULED if fomc_scheduled is None else fomc_scheduled)

    def col(c):
        return event_rows[c].astype(str) if c in event_rows else pd.Series([""] * len(event_rows), index=event_rows.index,
                                                                           dtype=str)

    date, orig, acc = col("date"), col("original_date"), col("access_date")
    typ, status = col("event_type"), col("status")
    cpi_months = {pd.Period(x, "M") for x in list(date[typ == "cpi"]) + list(orig[typ == "cpi"]) if x}
    fomc_rows = set(date[typ == "fomc"]) | set(orig[typ == "fomc"])
    known_to = _day(f"{max(sched)[:4]}-12-31") if sched else None
    first_month = pd.Period(_day(stream_start), "M") + 1
    ev_days = np.array([_day(x).value for x in date], dtype=np.int64)
    span = INTEGRITY_DAYS * 86_400 * 10**9
    bad_access = (status == "released").to_numpy() & (acc < date).to_numpy()
    out = []
    for r in table.itertuples(index=False):
        d = _day(r.date)
        lo, hi = d - pd.Timedelta(days=INTEGRITY_DAYS), d + pd.Timedelta(days=INTEGRITY_DAYS)
        why = [f"no CPI row for {m}" for m in pd.period_range(max(pd.Period(lo, "M"), first_month), pd.Period(hi, "M"))
               if m not in cpi_months]
        if known_to is None or hi > known_to:
            why.append(f"the FOMC schedule (FOMC_SCHEDULED) is known only to "
                       f"{'nothing' if known_to is None else f'{known_to:%Y-%m-%d}'}; extend it to {hi:%Y-%m-%d}")
        why += [f"no FOMC row for the scheduled meeting of {d0}" for d0 in sched
                if lo <= _day(d0) <= hi and d0 not in fomc_rows]
        near = np.abs(ev_days - d.value) <= span
        why += [f"{date.iloc[i]} {typ.iloc[i]}: 'released' with an access date before the release"
                for i in np.flatnonzero(near & bad_access)]
        out.append("; ".join(why))
    return out


def completion(table: pd.DataFrame, event_rows: pd.DataFrame, data_end=None, integrity: list[str] | None = None
               ) -> pd.DataFrame:
    """Rule 4. A stratum is complete when the data reach t + 28 days + 6 hours (``data_end``, exclusive, is at or after
    it; None: no data end) and no row of the event file within 30 days of the release day still says 'scheduled'.
    With ``integrity`` (the reasons of ``calendar_integrity``, one per stratum; the status part only), a stratum is
    also not complete while the event file around it is incomplete. Strata enter the test in the order i, so
    ``in_test`` is the complete prefix."""
    t = table.copy()
    ev_days = np.array([_day(d).value for d in event_rows["date"]], dtype=np.int64)
    sched = (event_rows["status"].astype(str) == "scheduled").to_numpy()
    end = None if data_end is None else _utc(data_end)
    data_ok, status_ok, why = [], [], []
    for r in t.itertuples(index=False):
        d = _day(r.date).value
        near = np.abs(ev_days - d) <= STATUS_DAYS * 86_400 * 10**9
        pending = sorted(str(x) for x in event_rows.loc[near & sched, "date"])
        dok = end is None or end >= r.data_needed_until
        reasons = []
        if not dok:
            reasons.append(f"data end {end:%Y-%m-%d %H:%M} before {r.data_needed_until:%Y-%m-%d %H:%M} UTC")
        if pending:
            reasons.append("scheduled rows within 30 days: " + ", ".join(pending))
        data_ok.append(dok)
        status_ok.append(not pending)
        why.append("; ".join(reasons))
    t["data_complete"] = data_ok
    t["status_complete"] = status_ok
    t["complete"] = t["data_complete"] & t["status_complete"]
    if integrity is not None:
        t["calendar_complete"] = [not x for x in integrity]
        t["complete"] = t["complete"] & t["calendar_complete"]
        why = ["; ".join(x for x in (w, f"event file: {i}" if i else "") if x) for w, i in zip(why, integrity)]
    t["in_test"] = np.cumprod(t["complete"].to_numpy(bool)).astype(bool) if len(t) else np.zeros(0, bool)
    t["not_complete_because"] = why
    return t


def window_overlaps(cands: pd.DataFrame) -> pd.DataFrame:
    """Pairs of windows of different strata (release windows and kept controls) that overlap in time. Rule 3 makes
    every window belong to one stratum, but windows at different clock times on one day can overlap (a CPI control at
    08:30 ET and an FOMC control at 14:00 ET share 30 minutes)."""
    w = []
    for s, g in cands.groupby("stratum"):
        w.append((g.release_utc.iloc[0], int(s), "release"))
        w += [(a, int(s), "control") for a in g.loc[g.kept, "anchor_utc"]]
    w.sort(key=lambda x: x[0])
    span = pd.Timedelta(minutes=WINDOW_MIN)
    rows = []
    for i, (a, s, kind) in enumerate(w):
        for b, s2, kind2 in w[i + 1:]:
            if b >= a + span:
                break
            if s2 != s:
                rows.append(dict(window_1=f"{a:%Y-%m-%d %H:%M}", stratum_1=s, kind_1=kind, window_2=f"{b:%Y-%m-%d %H:%M}",
                                 stratum_2=s2, kind_2=kind2, overlap_minutes=int((a + span - b) / pd.Timedelta(minutes=1))))
    return pd.DataFrame(rows, columns=["window_1", "stratum_1", "kind_1", "window_2", "stratum_2", "kind_2",
                                       "overlap_minutes"])


# ----------------------------------------------------------------------------
# Outcomes of windows and strata
# ----------------------------------------------------------------------------


def _ns(x) -> np.ndarray:
    idx = pd.DatetimeIndex(pd.Series(x)) if not isinstance(x, pd.DatetimeIndex) else x
    if idx.tz is None:
        idx = idx.tz_localize("UTC")
    return idx.tz_convert("UTC").as_unit("ns").asi8


def window_outcomes(up_times, anchors) -> np.ndarray:
    """1 if an up alarm falls in [anchor, anchor + 6 h). On a gap-free minute index this is E12.window_stats' up_post
    (an up alarm with index in [i0, i0 + POST)); gate 3 checks the equivalence on E22's and E23's committed files."""
    up = np.sort(_ns(up_times)) if len(up_times) else np.zeros(0, np.int64)
    a = _ns(anchors) if len(anchors) else np.zeros(0, np.int64)
    span = WINDOW_MIN * 60 * 10**9
    return np.searchsorted(up, a + span, side="left") > np.searchsorted(up, a, side="left")


def stratum_outcomes(cands: pd.DataFrame, table: pd.DataFrame, up_times) -> pd.DataFrame:
    """Z (release window), K (kept controls with the outcome), n, S, p and the e-values of each stratum in ``table``."""
    t = table.copy()
    z = window_outcomes(up_times, t["release_utc"]).astype(int)
    kept = cands[cands.kept & cands.stratum.isin(t.stratum)]
    c_out = window_outcomes(up_times, kept["anchor_utc"]).astype(int)
    k_by = pd.Series(c_out, index=kept.stratum.to_numpy()).groupby(level=0).sum()
    t["Z"] = z
    t["K"] = [int(k_by.get(s, 0)) for s in t.stratum]
    t["S"] = t["Z"] + t["K"]
    t["p"] = t["S"] / (t["n"] + 1)
    t["control_outcomes"] = [";".join(str(int(v)) for v in c_out[(kept.stratum == s).to_numpy()]) for s in t.stratum]
    return t


def add_e_process(t: pd.DataFrame) -> pd.DataFrame:
    t = t.copy()
    e = stratum_e(t["Z"].to_numpy(), t["S"].to_numpy(), t["n"].to_numpy())
    prod, m = e_process(e)
    for j, ps in enumerate(PSI):
        t[f"e_psi{int(ps)}"] = e[:, j]
    for j, ps in enumerate(PSI):
        t[f"prod_psi{int(ps)}"] = prod[:, j]
    t["M"] = m
    t["running_max_M"] = np.maximum.accumulate(m) if len(m) else m
    return t


def crossing(t: pd.DataFrame, threshold: float) -> tuple[int | None, str]:
    hit = np.flatnonzero(t["M"].to_numpy() >= threshold) if len(t) else np.array([])
    if not hit.size:
        return None, ""
    r = t.iloc[int(hit[0])]
    return int(hit[0]) + 1, str(r.date)


# ----------------------------------------------------------------------------
# Part theory: gate 1 (identity) and gate 5 (simulated validity)
# ----------------------------------------------------------------------------


def identity_table() -> pd.DataFrame:
    """Gate 1. For n = 0..8, S = 0..n+1 and psi in the grid: under H0 every arrangement of the S ones among the n + 1
    windows is equally likely given S; the release window is window 0. The exact mean of e(psi) over the
    arrangements, in rational arithmetic, and the floating-point mean of the code's e_value."""
    rows = []
    for n in range(0, 9):
        for s in range(0, n + 2):
            arr = list(itertools.combinations(range(n + 1), s))
            zs = np.array([1 if 0 in a else 0 for a in arr])
            for psi in PSI:
                q = Fraction(int(psi))
                p = Fraction(s, n + 1)
                exact = sum((q ** int(z)) / (1 - p + q * p) for z in zs) / len(arr)
                flt = float(np.mean(e_value(zs, s, n, psi)))
                rows.append(dict(n=n, S=s, psi=psi, arrangements=len(arr), p=s / (n + 1),
                                 e_if_release_has_outcome=float(e_value(1, s, n, psi)) if s >= 1 else np.nan,
                                 e_if_release_has_not=float(e_value(0, s, n, psi)) if s <= n else np.nan,
                                 mean_exact=str(exact), mean_exact_is_1=exact == 1, mean_float=flt,
                                 abs_deviation_float=abs(flt - 1.0)))
    return pd.DataFrame(rows)


def _cycle(ns, length: int) -> np.ndarray:
    ns = np.asarray(ns, int)
    return np.resize(ns, length)


def simulate_h0(ns: np.ndarray, law: str, rate: float, reps: int, rng: np.random.Generator,
                batch: int = 2000) -> dict:
    """M over len(ns) strata under H0 (exchangeable windows within each stratum given the earlier strata), for one
    law of the number S of windows with the outcome. Returns the shares of sequences on which M ever reaches 10, 20 and
    100 and the mean of the final M."""
    J = len(ns)
    n = np.asarray(ns, int)
    nmax = int(n.max())
    hits = {10: 0, 20: 0, 100: 0}
    at = {1: 0, 10: 1, J: 2}  # M after 1 and 10 strata and at the end (E[M_j] = 1 under H0)
    msum, msq = np.zeros(3), np.zeros(3)
    done = 0
    while done < reps:
        b = min(batch, reps - done)
        if law == "exchangeable Bernoulli":
            w = rng.random((b, J, nmax + 1)) < rate  # the n + 1 windows of each stratum, i.i.d.
            valid = np.arange(nmax + 1)[None, None, :] < (n + 1)[None, :, None]
            w &= valid
            slot = np.floor(rng.random((b, J)) * (n + 1)).astype(int)  # the release is one of the n + 1 windows
            z = np.take_along_axis(w, slot[..., None], axis=2)[..., 0].astype(int)
            s = w.sum(axis=2)
        else:
            if law == "all or nothing":
                s = np.where(rng.random((b, J)) < rate, n + 1, 0)
            elif law == "exactly one window with the outcome":
                s = np.ones((b, J), int)
            elif law == "S uniform on 0..n+1":
                s = np.floor(rng.random((b, J)) * (n + 2)).astype(int)
            elif law == "S chosen from the earlier strata":
                s = np.zeros((b, J), int)
            else:
                raise ValueError(law)
            u = rng.random((b, J))
            if law == "S chosen from the earlier strata":
                # S_i = 1 after a release without the outcome, n_i after one with it (1 for the first stratum)
                z = np.zeros((b, J), int)
                prev = np.zeros(b, int)
                for j in range(J):
                    s[:, j] = np.where(prev == 1, n[j], 1) if n[j] > 0 else np.minimum(1, n[j] + 1)
                    z[:, j] = (u[:, j] * (n[j] + 1) < s[:, j]).astype(int)
                    prev = z[:, j]
            else:
                z = (u * (n + 1)[None, :] < s).astype(int)  # the release position is uniform given S
        _, m = e_process(stratum_e(z, s, n[None, :]))
        mx = m.max(axis=1)
        for th in hits:
            hits[th] += int((mx >= th).sum())
        for j, slot in at.items():
            msum[slot] += float(m[:, j - 1].sum())
            msq[slot] += float((m[:, j - 1] ** 2).sum())
        done += b
    mean = msum / reps
    se = np.sqrt(np.maximum(msq / reps - mean**2, 0.0) / reps)
    return dict(share_M_ever_ge_10=hits[10] / reps, share_M_ever_ge_20=hits[20] / reps,
                share_M_ever_ge_100=hits[100] / reps, mean_M_1=mean[0], mean_M_1_se=se[0], mean_M_10=mean[1],
                mean_M_10_se=se[1], mean_M_final=mean[2], mean_M_final_se=se[2])


def part_theory():
    idt = identity_table()
    idt.to_csv(OUT / "e25_theory_identity.csv", index=False)
    dev = float(idt.abs_deviation_float.max())
    ok1 = bool(idt.mean_exact_is_1.all()) and dev <= 1e-12
    print(f"gate 1: {len(idt)} cases (n = 0..8, S = 0..n+1, psi in {PSI}); exact mean 1 in every case: "
          f"{bool(idt.mean_exact_is_1.all())}; largest floating-point deviation {dev:.2e}; {'PASS' if ok1 else 'FAIL'}")
    cal = planning_table()
    ns = _cycle(cal["n"], VALIDITY_RELEASES)
    se_nom = np.sqrt(ALPHA * (1 - ALPHA) / REPS)
    rows = []
    laws = [("exchangeable Bernoulli", r) for r in VALIDITY_RATES] + [
        ("all or nothing", 0.2), ("exactly one window with the outcome", np.nan), ("S uniform on 0..n+1", np.nan),
        ("S chosen from the earlier strata", np.nan)]
    for i, (law, rate) in enumerate(laws):
        rng = np.random.default_rng([SEED, i])
        res = simulate_h0(ns, law, rate, REPS, rng)
        sh = res["share_M_ever_ge_20"]
        rows.append(dict(law=law, rate=rate, releases=VALIDITY_RELEASES, reps=REPS, seed=f"[{SEED}, {i}]",
                         n_source=f"planned calendar to 2029-01-01 (results/e25_power_calendar{_suffix()}.csv), cycled",
                         **res, mc_se_20=float(np.sqrt(sh * (1 - sh) / REPS)), bound=ALPHA,
                         limit_bound_plus_3se=ALPHA + 3 * se_nom, passes=sh <= ALPHA + 3 * se_nom))
        print(f"gate 5: {law}, rate {rate}: P(M ever >= 20 within {VALIDITY_RELEASES}) = {sh:.4f} "
              f"(>= 10: {res['share_M_ever_ge_10']:.4f}, >= 100: {res['share_M_ever_ge_100']:.4f}); mean of M after 1, 10 "
              f"and {VALIDITY_RELEASES} strata {res['mean_M_1']:.3f} ({res['mean_M_1_se']:.3f}), {res['mean_M_10']:.3f} "
              f"({res['mean_M_10_se']:.3f}), {res['mean_M_final']:.3f} ({res['mean_M_final_se']:.3f})", flush=True)
    val = pd.DataFrame(rows)
    val["rules"] = RULES
    val.to_csv(OUT / f"e25_validity{_suffix()}.csv", index=False)
    print(f"gate 5: {'PASS' if val.passes.all() else 'FAIL'} (limit {ALPHA + 3 * se_nom:.4f})")


# ----------------------------------------------------------------------------
# Part strata: the planned forward calendar
# ----------------------------------------------------------------------------


def future_event_rows() -> pd.DataFrame:
    return pd.read_csv(FUTURE_EVENTS, dtype=str, keep_default_na=False)


def planned_releases() -> pd.DataFrame:
    """The forward calendar for planning: the event file's released, rescheduled and scheduled rows (scheduled treated
    as released, for planning only), minus the statement of 27 January 2027 (F24a's planning convention)."""
    ev = future_event_rows()
    ev = ev[ev.status.isin(RELEASED + ("scheduled",)) & ~ev.date.isin(E24.PLANNING_NO_STATEMENT)]
    return ev[["date", "event_type", "release_time_et", "status"]].reset_index(drop=True)


def part_strata():
    rel = planned_releases()
    cands = assign_strata(rel, rel["date"], STREAM_START, None)
    tab = strata_table(cands)
    ev = future_event_rows()
    comp = completion(tab, ev.assign(status="released"), None)  # planning: every row treated as released
    ov = window_overlaps(cands)
    out = cands.merge(comp[["stratum", "n", "data_needed_until"]], on="stratum")
    sched = ev[ev.status == "scheduled"]
    out["scheduled_rows_within_30_days_now"] = [
        ";".join(sorted(d for d in sched.date if abs((_day(d) - _day(r)).days) <= STATUS_DAYS)) for r in out.date]
    ov_key = {}
    for r in ov.itertuples():
        ov_key.setdefault(r.window_1, []).append(f"{r.window_2} (stratum {r.stratum_2}, {r.overlap_minutes} min)")
        ov_key.setdefault(r.window_2, []).append(f"{r.window_1} (stratum {r.stratum_1}, {r.overlap_minutes} min)")
    out["overlaps_window_of_other_stratum"] = [";".join(ov_key.get(f"{a:%Y-%m-%d %H:%M}", [])) if k else ""
                                               for a, k in zip(out.anchor_utc, out.kept)]
    out = out.merge(comp[["stratum", "controls_before", "controls_after", "balanced", "release_window_overlaps"]],
                    on="stratum")
    out["rules"] = RULES
    out.to_csv(OUT / f"e25_strata_planned{_suffix()}.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_rows", 400, "display.max_colwidth", 120):
        print(comp[["stratum", "date", "event_type", "release_utc", "n", "controls_before", "controls_after",
                    "data_needed_until", "controls"]].to_string(index=False))
        print(out.groupby("reason").size().to_string())
        print(f"rules {RULES}: {len(comp)} strata, n from {comp.n.min()} to {comp.n.max()} "
              f"({';'.join(f'{v}:{int((comp.n == v).sum())}' for v in range(9) if (comp.n == v).any())}), "
              f"{int(comp.balanced.sum())} balanced, {int(_one_sided(comp).sum())} "
              f"with controls on one side only (kept controls before the release {int(comp.controls_before.sum())}, after "
              f"{int(comp.controls_after.sum())}); release windows that overlap: "
              f"{int((comp.release_window_overlaps != '').sum())} strata")
        print(f"windows of different strata that overlap: {len(ov)} pairs")
        if len(ov):
            print(ov.to_string(index=False))
    return out, comp


# ----------------------------------------------------------------------------
# Part power (section 9)
# ----------------------------------------------------------------------------


def planning_table(end: str = PLANNING_END) -> pd.DataFrame:
    """F24a's planning calendar (E24.planning_events: the event file's released and scheduled rows minus 27 January
    2027, plus the 2028 placeholders) with E25's strata (rules 1 to 3, no data end) and the time each stratum enters
    the test under the planning assumption that every event-file row is updated at its release time (rule 4)."""
    ev = E24.planning_events(end)
    cands = assign_strata(ev, ev["date"], STREAM_START, None)
    tab = strata_table(cands)
    tab["kind"] = ev.sort_values("release_utc")["kind"].to_numpy()
    # rows whose status must be settled: every event-file row (any status) and the placeholders, settled at release
    rows = future_event_rows()[["date", "event_type", "release_time_et"]]
    ph = ev[ev.kind != "event file"][["date", "event_type", "release_time_et"]].assign(
        date=lambda d: d["date"].dt.strftime("%Y-%m-%d"))
    rows = pd.concat([rows, ph], ignore_index=True)
    settle = np.array([et_anchor(d, h).value for d, h in zip(rows.date, rows.release_time_et)], dtype=np.int64)
    rdays = np.array([_day(d).value for d in rows.date], dtype=np.int64)
    entry = []
    for r in tab.itertuples(index=False):
        near = np.abs(rdays - _day(r.date).value) <= STATUS_DAYS * 86_400 * 10**9
        entry.append(max(r.data_needed_until.value, int(settle[near].max())))
    tab["enters_test"] = pd.to_datetime(np.maximum.accumulate(np.array(entry, dtype=np.int64)), utc=True)
    return tab


def _rates_matrix(tab: pd.DataFrame, rates: dict) -> tuple[np.ndarray, np.ndarray]:
    p1 = np.array([rates[t][0] for t in tab.event_type])
    p0 = np.array([rates[t][1] for t in tab.event_type])
    return p1, p0


def simulate_power(tab: pd.DataFrame, rates: dict, reps: int = REPS, seed: int = SEED, batch: int = 5000) -> dict:
    """Planning model of F24a, section 5: the release window has the outcome with probability p1 of its type, each
    control window with p0 of its type, all independently. Common random numbers across scenarios (uniforms compared
    with the rates). Returns, for E25's e-process and for the oracle that knows p0, the index (1-based) of the first
    stratum at which M >= 20 (0 if none) for every replication, and the final M."""
    rng = np.random.default_rng(seed)
    n = tab["n"].to_numpy(int)
    J, nmax = len(n), int(n.max())
    p1, p0 = _rates_matrix(tab, rates)
    first, first_or, fin, fin_or = [], [], [], []
    done = 0
    valid = np.arange(nmax)[None, None, :] < n[None, :, None]
    while done < reps:
        b = min(batch, reps - done)
        ur = rng.random((b, J))
        uc = rng.random((b, J, nmax))
        z = (ur < p1[None, :]).astype(int)
        k = ((uc < p0[None, :, None]) & valid).sum(axis=2)
        _, m = e_process(stratum_e(z, z + k, n[None, :]))
        eo = np.stack([e_value_oracle(z, p0[None, :], ps) for ps in PSI], axis=-1)
        _, mo = e_process(eo)
        for mm, out, f in ((m, first, fin), (mo, first_or, fin_or)):
            hit = mm >= THRESHOLD
            out.append(np.where(hit.any(axis=1), hit.argmax(axis=1) + 1, 0))
            f.append(mm[:, -1])
        done += b
    return dict(first=np.concatenate(first), first_oracle=np.concatenate(first_or), final=np.concatenate(fin),
                final_oracle=np.concatenate(fin_or))


def expected_log_growth(tab: pd.DataFrame, rates: dict, oracle: bool = False) -> np.ndarray:
    """Exact E[log e_i(psi)] per stratum and psi under the planning model (shape J x len(PSI))."""
    from scipy.stats import binom

    p1, p0 = _rates_matrix(tab, rates)
    out = np.zeros((len(tab), len(PSI)))
    for i, n in enumerate(tab["n"].to_numpy(int)):
        for j, ps in enumerate(PSI):
            tot = 0.0
            for z in (0, 1):
                pz = p1[i] if z else 1 - p1[i]
                if oracle:
                    tot += pz * np.log(e_value_oracle(z, p0[i], ps))
                    continue
                for k in range(n + 1):
                    tot += pz * binom.pmf(k, n, p0[i]) * np.log(e_value(z, z + k, n, ps))
            out[i, j] = tot
    return out


def part_power():
    tab = planning_table()
    inp = pd.read_csv(OUT / "e24_future_planning_inputs.csv")
    scen = E24.planning_scenarios(inp, PRIMARY_ARM)
    f24 = pd.read_csv(OUT / "e24_future_power.csv")
    f24 = f24[(f24.arm == PRIMARY_ARM) & (f24.design == "fomc+cpi")].set_index("scenario")
    check = {s: round(float(f24.loc[s, "power"]), 3) for s in SCENARIOS}
    print(f"F24a's power at {F24A_POWER_AT} (results/e24_future_power.csv, S2 pooled): {check}")
    if list(check.values()) != [0.822, 0.596, 0.369, 0.031]:
        raise SystemExit("F24a's committed power values are not 0.822, 0.596, 0.369, 0.031")
    entry = tab["enters_test"]
    entry_ns = _ns(entry)
    # the crossing probability changes only when a stratum enters: evaluate it at every entry time and month start
    grid = pd.DatetimeIndex(sorted(set(pd.date_range("2026-11-01", PLANNING_END, freq="MS", tz="UTC"))
                                   | set(entry[entry <= pd.Timestamp(PLANNING_END, tz="UTC")])))
    never = np.iinfo(np.int64).max
    rows, curve, growth = [], [], []
    for sname in SCENARIOS:
        rates = scen[sname]
        res = simulate_power(tab, rates)
        for kind, fk in (("E25 (p estimated from the stratum)", "first"), ("oracle (true control rate)", "first_oracle")):
            first = res[fk]
            ok = first > 0
            cdate = np.full(first.size, never, dtype=np.int64)  # time the crossing stratum enters the test
            cdate[ok] = entry_ns[first[ok] - 1]
            row = dict(scenario=sname, procedure=kind, p1_fomc=rates["fomc"][0], p0_fomc=rates["fomc"][1],
                       p1_cpi=rates["cpi"][0], p0_cpi=rates["cpi"][1], reps=REPS, seed=SEED)
            for d in POWER_DATES:
                hit = cdate <= pd.Timestamp(d, tz="UTC").value
                row[f"strata_entered_by_{d}"] = int((entry_ns <= pd.Timestamp(d, tz="UTC").value).sum())
                row[f"P_cross_by_{d}"] = float(hit.mean())
                row[f"se_{d}"] = float(np.sqrt(hit.mean() * (1 - hit.mean()) / REPS))
            okl = cdate <= pd.Timestamp(PLANNING_END, tz="UTC").value  # crossed by 1 January 2029
            srt = np.sort(np.where(okl, first, never))
            row["median_releases_to_crossing_all"] = float(srt[REPS // 2]) if srt[REPS // 2] < never else np.nan
            row["median_releases_to_crossing_given_crossing"] = float(np.median(first[okl])) if okl.any() else np.nan
            row["median_crossing_date_given_crossing"] = (
                str(pd.Timestamp(int(np.median(cdate[okl])), unit="ns", tz="UTC").date()) if okl.any() else "")
            row["share_never_crossing_by_2029-01-01"] = float(1 - okl.mean())
            if kind.startswith("E25"):
                row["F24a_power_at_2027-06-01"] = float(f24.loc[sname, "power"])
            rows.append(row)
            for g in grid:
                curve.append(dict(scenario=sname, procedure=kind, date=f"{g:%Y-%m-%d %H:%M}",
                                  P_cross_by_date=float((cdate <= g.value).mean()),
                                  strata_entered=int((entry_ns <= g.value).sum())))
        sub = tab[entry_ns <= pd.Timestamp(PLANNING_END, tz="UTC").value]
        for kind, orc in (("E25 (p estimated from the stratum)", False), ("oracle (true control rate)", True)):
            g = expected_log_growth(sub, rates, oracle=orc)
            for j, ps in enumerate(PSI):
                growth.append(dict(scenario=sname, procedure=kind, psi=ps, strata=len(g),
                                   mean_expected_log_e_per_release=float(g[:, j].mean()),
                                   mean_fomc=float(g[(sub.event_type == "fomc").to_numpy(), j].mean()),
                                   mean_cpi=float(g[(sub.event_type == "cpi").to_numpy(), j].mean())))
        print(f"{sname}: done", flush=True)
    pw = pd.DataFrame(rows)
    sfx = _suffix()
    pw.assign(rules=RULES).to_csv(OUT / f"e25_power{sfx}.csv", index=False)
    pd.DataFrame(curve).assign(rules=RULES).to_csv(OUT / f"e25_power_curve{sfx}.csv", index=False)
    pd.DataFrame(growth).assign(rules=RULES).to_csv(OUT / f"e25_power_growth{sfx}.csv", index=False)
    tab.assign(release_utc=tab.release_utc.astype(str), data_needed_until=tab.data_needed_until.astype(str),
               enters_test=tab.enters_test.astype(str), rules=RULES).to_csv(OUT / f"e25_power_calendar{sfx}.csv",
                                                                              index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 200):
        print(tab[["stratum", "date", "event_type", "kind", "n", "enters_test"]].to_string(index=False))
        cols = ["scenario", "procedure"] + [f"P_cross_by_{d}" for d in POWER_DATES] + [
            "median_releases_to_crossing_all", "median_releases_to_crossing_given_crossing",
            "median_crossing_date_given_crossing", "share_never_crossing_by_2029-01-01", "F24a_power_at_2027-06-01"]
        print(pw[cols].round(4).to_string(index=False))
        print(pd.DataFrame(growth).round(4).to_string(index=False))
    if RULES == "amended":
        fig_power(_plt())


# ----------------------------------------------------------------------------
# Market data: S2 runs on known periods (gates, post hoc) and at a look (status)
# ----------------------------------------------------------------------------


def period_e22() -> dict:
    return dict(E23.E12_PERIOD)


def period_e23(end: str | None = None) -> dict:
    p = E24._e23_period()
    if end is not None:
        p["end"] = end
    return p


def run_arm(args):
    """One continuous run of an arm on one coin and period, exactly as E22.job runs it (E12.load_coin through E23's
    period switch, E22.make_arm, E22.run_monitor, E12.alarm_frame). Returns the alarms (with their UTC times), the
    monitoring flag and the deseasonalised returns."""
    coin, arm, period = args
    t0 = time.time()
    if arm in E24.REAL_ARMS:
        E24._patch_e22()
    with E23.e12_period(**period):
        df, u = E12.load_coin(E12.COINS[coin])
    index = df.index
    mon = E22.make_arm(arm, E22.utc_day_offset(index))
    monitoring = E22.run_monitor(mon, u)
    ev = E12.alarm_frame(mon.events, 0, u)
    ev["time"] = index[ev["index"].to_numpy()] if len(ev) else pd.DatetimeIndex([], tz="UTC")
    return dict(coin=coin, arm=arm, period=period, alarms=ev, monitoring=monitoring, u=u, first=index[0],
                n=len(index), seconds=time.time() - t0)


def _run_many(jobs, procs: int):
    if procs <= 1:
        return [run_arm(j) for j in jobs]
    with Pool(min(procs, 2)) as pool:
        return pool.map(run_arm, jobs, chunksize=1)


def _up_times(alarms: pd.DataFrame):
    return pd.DatetimeIndex(alarms.loc[alarms.direction == "up", "time"])


def _committed_alarms(name: str, arm: str = PRIMARY_ARM) -> pd.DataFrame:
    a = pd.read_csv(OUT / name)
    a = a[a.arm == arm].copy()
    a["time"] = pd.to_datetime(a["time"], utc=True)
    return a


# ----------------------------------------------------------------------------
# Part gates: causality (gate 2) and pipeline (gate 3), known data only
# ----------------------------------------------------------------------------


def part_gates(procs: int):
    p23, p22 = period_e23(), period_e22()
    if pd.Timestamp(p23["end"]) > pd.Timestamp("2026-10-25"):
        raise SystemExit("E23's period must end by 2026-10-25 00:00 UTC (known data only)")
    jobs = [(c, PRIMARY_ARM, p23) for c in E12.COINS] + [(PRIMARY_COIN, PRIMARY_ARM, period_e23(cut))
                                                         for cut in CAUSALITY_CUTS] + [(c, PRIMARY_ARM, p22) for c in E12.COINS]
    res = _run_many(jobs, procs)
    for r in res:
        print(f"{r['coin']} {r['arm']} {r['period']['monitor_start']} to {r['period']['end']}: {len(r['alarms'])} alarms, "
              f"{r['seconds']:.0f} s", flush=True)
    rows = []
    full = {(r["coin"], r["period"]["end"]): r for r in res}
    # gate 2: causality
    f = full[(PRIMARY_COIN, p23["end"])]
    key = ["index", "direction", "delay", "baseline_scale", "max_abs_recent_over_scale"]
    for cut in CAUSALITY_CUTS:
        c = full[(PRIMARY_COIN, cut)]
        cut_t = pd.Timestamp(cut, tz="UTC")
        a_full = f["alarms"][f["alarms"].time < cut_t].reset_index(drop=True)
        a_cut = c["alarms"].reset_index(drop=True)
        same = len(a_full) == len(a_cut) and all(np.array_equal(a_full[k].to_numpy(), a_cut[k].to_numpy()) for k in key)
        m = c["n"]
        u_same = np.array_equal(f["u"][:m], c["u"], equal_nan=True)
        mon_same = np.array_equal(f["monitoring"][:m], c["monitoring"])
        rows.append(dict(gate="2 causality", period="2026 (E23 design)", coin=PRIMARY_COIN, check=f"cut at {cut}",
                         compared=len(a_full), mismatches=int(not same) + int(not u_same) + int(not mon_same),
                         detail=f"alarms before the cut: full run {len(a_full)}, cut run {len(a_cut)}; columns {key} "
                                f"identical: {same}; deseasonalised returns before the cut identical: {u_same}; "
                                f"monitoring flags identical: {mon_same}"))
    # gate 3: own runs reproduce the committed alarm files; release outcomes equal E23's (and E22's) up_post
    for name, per, wname, plabel in (("e23_alarms.csv", p23, "e23_windows.csv", "2026 (E23 design)"),
                                     ("e22_alarms.csv", p22, "e22_windows.csv", "2024 to 2025 (E22 design)")):
        com = _committed_alarms(name)
        win = pd.read_csv(OUT / wname)
        win = win[win.arm == PRIMARY_ARM]
        for coin in E12.COINS:
            own = full[(coin, per["end"])]["alarms"].reset_index(drop=True)
            cm = com[com.coin == coin].reset_index(drop=True)
            same = (len(own) == len(cm) and np.array_equal(own["index"].to_numpy(), cm["index"].to_numpy())
                    and np.array_equal(own["direction"].to_numpy(), cm["direction"].to_numpy())
                    and np.array_equal(_ns(own["time"]), _ns(cm["time"])))
            rows.append(dict(gate="3 pipeline", period=plabel, coin=coin,
                             check=f"E25's S2 run reproduces the committed S2 alarms ({name})", compared=len(cm),
                             mismatches=int(not same), detail=f"own {len(own)} alarms, committed {len(cm)}"))
            # the time column of the committed file is the minute of the index position (gap-free index)
            t_from_index = pd.Timestamp(per["monitor_start"], tz="UTC") + pd.to_timedelta(cm["index"], unit="min")
            rows.append(dict(gate="3 pipeline", period=plabel, coin=coin,
                             check=f"committed {name}: time = stream start + index minutes", compared=len(cm),
                             mismatches=int((_ns(t_from_index) != _ns(cm["time"])).sum()), detail=""))
            w = win[win.coin == coin]
            rel = pd.to_datetime(w["release_utc"], utc=True)
            want = w["up_post"].astype(bool).to_numpy()
            for src, up in (("E25's own S2 run", _up_times(own)), (f"committed {name}", _up_times(cm))):
                z = window_outcomes(up, rel)
                bad = int((z != want).sum())
                rows.append(dict(gate="3 pipeline", period=plabel, coin=coin,
                                 check=f"release outcomes Z from {src} equal up_post of {wname} (S2)",
                                 compared=len(want), mismatches=bad,
                                 detail=f"releases with the outcome: {int(z.sum())} (E25), {int(want.sum())} ({wname})"))
            # controls: E25's window outcome on E12's own matched control anchors reproduces control_up_post_rate
            with E23.e12_period(**per):
                events, _, _ = E12.read_events()
                index = pd.date_range(pd.Timestamp(E12.MONITOR_START, tz="UTC"), pd.Timestamp(E12.END, tz="UTC"),
                                      freq="1min", inclusive="left")
                days = events["date"].to_list()
                bad, cnt = 0, 0
                for _, r in events.iterrows():
                    anchors = E12.control_anchors(index, days, r["date"].weekday(), r["release_time_et"], index.size)
                    rate = float(np.mean(window_outcomes(_up_times(cm), index[anchors])))
                    ref = w[pd.to_datetime(w.release_utc, utc=True) == r["release_utc"]]["control_up_post_rate"]
                    cnt += 1
                    bad += int(len(ref) != 1 or not np.isclose(rate, float(ref.iloc[0]), rtol=0, atol=1e-12))
            rows.append(dict(gate="3 pipeline", period=plabel, coin=coin,
                             check=f"E25's window outcome on E12's control anchors equals control_up_post_rate of {wname}",
                             compared=cnt, mismatches=bad, detail="E12.control_anchors on a gap-free minute index"))
    g = pd.DataFrame(rows)
    g.to_csv(OUT / "e25_gates.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 200):
        print(g.to_string(index=False))
    print(f"gates 2 and 3: {len(g)} checks, {int((g.mismatches > 0).sum())} with mismatches")
    return g


# ----------------------------------------------------------------------------
# Part posthoc (section 10): seen data, exploratory
# ----------------------------------------------------------------------------


def period_strata(events_file, monitor_start, end) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Section 4 within one known period: the releases of the event file (released or rescheduled) at or after the
    stream start and before the data end, every release day of the file for rule 2, the period's end as the data end."""
    ev = pd.read_csv(events_file, dtype=str, keep_default_na=False)
    relrows = ev[ev.status.isin(RELEASED)]
    utc = [et_anchor(d, h) for d, h in zip(relrows.date, relrows.release_time_et)]
    inside = np.array([(t >= _utc(monitor_start)) and (t < _utc(end)) for t in utc], dtype=bool)
    rel = relrows[inside]
    cands = assign_strata(rel, relrows["date"], monitor_start, end)
    comp = completion(strata_table(cands), ev, end)
    return cands, comp


def period_gaps(per: dict) -> pd.DataFrame:
    """Runs of missing minutes of the three coins over a known period's monitored minutes (files already on disk,
    no download), for the missing minutes inside each window (review round 11, m10)."""
    from avmonitor import data as Dm

    out = []
    for coin, sym in E12.COINS.items():
        df = Dm.load_binance(sym, per["monitor_start"], per["end"], download=False)
        out.append(gap_runs(df.index, df["close"].isna().to_numpy(), coin))
    return pd.concat(out, ignore_index=True)


def side_rates(cands: pd.DataFrame, strata, up) -> dict:
    """Exploratory (review round 11, m2): kept controls before and after their release, and those with the outcome."""
    kept = cands[cands.kept & cands.stratum.isin(strata)]
    o = window_outcomes(up, kept["anchor_utc"]).astype(int)
    b = (kept.k < 0).to_numpy()
    nb, na = int(b.sum()), int((~b).sum())
    return dict(control_windows_before=nb, controls_before_with_outcome=int(o[b].sum()),
                control_rate_before=float(o[b].mean()) if nb else np.nan,
                control_windows_after=na, controls_after_with_outcome=int(o[~b].sum()),
                control_rate_after=float(o[~b].mean()) if na else np.nan)


def part_posthoc(procs: int):
    g = pd.read_csv(OUT / "e25_gates.csv")
    g3 = g[g.gate == "3 pipeline"]
    use_committed = len(g3) > 0 and int(g3.mismatches.sum()) == 0
    print(f"gate 3 holds: {use_committed}; outcomes from the committed S2 alarm files" if use_committed
          else "gate 3 failed or missing: rerunning S2")
    periods = (("2024 to 2025 (E22 design)", period_e22(), "e22_alarms.csv"),
               ("2026 (E23 design)", period_e23(), "e23_alarms.csv"))
    srows, summ = [], []
    for plabel, per, aname in periods:
        cands, comp = period_strata(per["events"], per["monitor_start"], per["end"])
        gaps = period_gaps(per)
        if use_committed:
            com = _committed_alarms(aname)
            ups = {c: _up_times(com[com.coin == c]) for c in E12.COINS}
        else:
            res = _run_many([(c, PRIMARY_ARM, per) for c in E12.COINS], procs)
            ups = {r["coin"]: _up_times(r["alarms"]) for r in res}
        for coin in E12.COINS:
            out = add_missing(stratum_outcomes(cands, comp, ups[coin]), cands, gaps, coin)
            t = add_e_process(out[out.in_test])
            for cell in ("fomc+cpi", "fomc", "cpi", "fomc+cpi, balanced strata only"):
                sub = (out[out.in_test] if cell == "fomc+cpi" else
                       out[out.in_test & out.balanced.astype(bool)] if cell.endswith("balanced strata only") else
                       out[out.in_test & (out.event_type == cell)])
                tt = add_e_process(sub)
                j20, d20 = crossing(tt, THRESHOLD)
                j40, d40 = crossing(tt, THRESHOLD_E25B)
                nk = tt["n"].to_numpy(int)
                summ.append(dict(
                    period=plabel, coin=coin, cell=cell, label="post hoc, exploratory (data seen before the design)",
                    releases_in_period=int(len(out) if cell.startswith("fomc+cpi") else (out.event_type == cell).sum()),
                    strata_in_test=len(tt), strata_not_complete=int((~out.in_test).sum()) if cell.startswith("fomc+cpi")
                    else int((~out.in_test & (out.event_type == cell)).sum()),
                    releases_with_outcome=int(tt.Z.sum()), release_rate=float(tt.Z.mean()) if len(tt) else np.nan,
                    control_windows=int(nk.sum()), controls_with_outcome=int(tt.K.sum()),
                    control_rate_pooled=float(tt.K.sum() / nk.sum()) if nk.sum() else np.nan,
                    control_rate_mean_over_strata=float(np.mean(tt.K[nk > 0] / nk[nk > 0])) if (nk > 0).any() else np.nan,
                    strata_uninformative=int(((tt.S == 0) | (tt.S == tt.n + 1)).sum()),
                    n_min=int(nk.min()) if len(nk) else np.nan, n_median=float(np.median(nk)) if len(nk) else np.nan,
                    n_max=int(nk.max()) if len(nk) else np.nan,
                    n_counts=";".join(f"{v}:{int((nk == v).sum())}" for v in range(9) if (nk == v).any()),
                    strata_balanced=int(tt.balanced.sum()),
                    strata_one_sided=int(_one_sided(tt).sum()),
                    strata_with_overlapping_release_window=int((tt.release_window_overlaps != "").sum()),
                    **side_rates(cands, tt.stratum, ups[coin]),
                    windows_with_missing_minutes=int(tt.windows_with_missing_minutes.sum()),
                    missing_minutes_in_windows=int(tt.missing_minutes_total.sum()),
                    final_M=float(tt.M.iloc[-1]) if len(tt) else 1.0, max_M=float(tt.M.max()) if len(tt) else 1.0,
                    first_crossing_20_stratum=j20, first_crossing_20_date=d20,
                    first_crossing_40_stratum=j40, first_crossing_40_date=d40,
                    **{f"final_prod_psi{int(ps)}": float(tt[f"prod_psi{int(ps)}"].iloc[-1]) if len(tt) else 1.0
                       for ps in PSI}))
            full = out.merge(t[["stratum"] + [c for c in t.columns if c.startswith(("e_psi", "prod_psi"))]
                               + ["M", "running_max_M"]], on="stratum", how="left")
            full.insert(0, "coin", coin)
            full.insert(0, "period", plabel)
            srows.append(full)
    st = pd.concat(srows, ignore_index=True)
    st["label"] = "post hoc, exploratory"
    for c in ("release_utc", "data_needed_until"):
        st[c] = st[c].astype(str)
    st["rules"] = RULES
    st.to_csv(OUT / f"e25_posthoc_strata{_suffix()}.csv", index=False)
    sm = pd.DataFrame(summ)
    sm["rules"] = RULES
    sm.to_csv(OUT / f"e25_posthoc_summary{_suffix()}.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 200):
        cols = ["period", "coin", "cell", "strata_in_test", "strata_not_complete", "releases_with_outcome",
                "release_rate", "control_rate_pooled", "n_counts", "strata_balanced", "strata_one_sided",
                "strata_with_overlapping_release_window", "control_rate_before", "control_rate_after",
                "missing_minutes_in_windows", "final_M", "max_M", "first_crossing_20_date", "first_crossing_40_date"]
        print(sm[cols].round(4).to_string(index=False))
        print(st[st.coin == PRIMARY_COIN][["period", "stratum", "date", "event_type", "n", "Z", "K", "S", "M",
                                           "in_test", "not_complete_because"]].round(3).to_string(index=False))
    if RULES == "amended":
        fig_posthoc(_plt())


# ----------------------------------------------------------------------------
# Part status: the fresh test at a look (never run before the first stratum is complete)
# ----------------------------------------------------------------------------


def fresh_plan(events: pd.DataFrame, data_end, stream_start=None, fomc_scheduled=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The strata of the fresh test at a look with data to ``data_end`` (exclusive): the releases of the event file
    (released or rescheduled, at or after the stream start), every release day of the file for rule 2, and the
    completion of rule 4 against the event file's statuses and (review round 11, M1) its completeness around each
    stratum (``calendar_integrity``). Uses no market data."""
    stream_start = STREAM_START if stream_start is None else stream_start
    relrows = events[events.status.isin(RELEASED)]
    utc = [et_anchor(d, h) for d, h in zip(relrows.date, relrows.release_time_et)]
    rel = relrows[np.array([t >= _utc(stream_start) for t in utc], dtype=bool)]
    cands = assign_strata(rel, relrows["date"], stream_start, data_end, rules="amended")  # the fresh test: Amendment 1
    tab = strata_table(cands)
    comp = completion(tab, events, data_end, integrity=calendar_integrity(tab, events, stream_start, fomc_scheduled))
    return cands, comp


def status_refusals(events: pd.DataFrame, cutoff, now=None, stream_start=None, fomc_scheduled=None) -> list[str]:
    """Reasons the status part must not compute anything (empty when it may). Checked before any market data are
    downloaded or read, and again with the data's actual end."""
    stream_start = STREAM_START if stream_start is None else stream_start
    now = _utc(pd.Timestamp.now(tz="UTC") if now is None else now)
    cut = _utc(cutoff)
    why = []
    if cut > now:
        why.append(f"the cut-off {cut:%Y-%m-%d %H:%M} UTC is after the time of the run {now:%Y-%m-%d %H:%M} UTC")
    if cut <= _utc(stream_start):
        why.append("the cut-off is not after the stream start")
        return why
    cands, comp = fresh_plan(events, cut, stream_start, fomc_scheduled)
    if not len(comp):
        why.append("no release in the event file at or after the stream start")
    elif not comp.in_test.any():
        first = comp.iloc[0]
        why.append(f"no stratum is complete at {cut:%Y-%m-%d %H:%M} UTC: stratum 1 ({first.date} {first.event_type}) "
                   f"is not complete because {first.not_complete_because}")
    return why


def missing_in_windows(gaps: pd.DataFrame | None, coin: str, anchors) -> np.ndarray:
    """Minutes of missing data inside each window [anchor, anchor + 6 h) (review round 11, m10; descriptive). ``gaps``
    lists runs of missing minutes as E23.data_checks writes them (coin, first_missing, minutes). S2 cannot alarm on a
    missing minute, so such a window counts as one without an alarm; the rules are not changed."""
    a = _ns(anchors) if len(anchors) else np.zeros(0, np.int64)
    if gaps is None or not len(gaps) or not len(a):
        return np.zeros(len(a), int)
    g = gaps[gaps.coin == coin]
    if not len(g):
        return np.zeros(len(a), int)
    g0 = _ns(pd.to_datetime(g["first_missing"], utc=True))
    g1 = g0 + g["minutes"].to_numpy(np.int64) * 60 * 10**9
    a1 = a + WINDOW_MIN * 60 * 10**9
    ov = np.minimum(a1[:, None], g1[None, :]) - np.maximum(a[:, None], g0[None, :])
    return (np.clip(ov, 0, None).sum(axis=1) // (60 * 10**9)).astype(int)


def gap_runs(index: pd.DatetimeIndex, missing: np.ndarray, coin: str) -> pd.DataFrame:
    """Runs of missing minutes in E23.data_checks' layout (coin, first_missing, minutes)."""
    m = np.asarray(missing, bool)
    if not m.any():
        return pd.DataFrame(columns=["coin", "first_missing", "minutes"])
    d = np.diff(np.concatenate(([0], m.astype(np.int8), [0])))
    starts, stops = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    return pd.DataFrame(dict(coin=coin, first_missing=[str(index[i]) for i in starts], minutes=stops - starts))


def add_missing(t: pd.DataFrame, cands: pd.DataFrame, gaps: pd.DataFrame | None, coin: str) -> pd.DataFrame:
    """Per stratum of ``t``: missing minutes in the release window, in each kept control window (in the order of
    ``controls``) and in all of them, and the number of the stratum's windows with any missing minute."""
    t = t.copy()
    rel = missing_in_windows(gaps, coin, t["release_utc"]) if len(t) else np.zeros(0, int)
    kept = cands[cands.kept & cands.stratum.isin(t.stratum)].sort_values(["stratum", "anchor_utc"])
    ctl = missing_in_windows(gaps, coin, kept["anchor_utc"])
    by = {s: ctl[(kept.stratum == s).to_numpy()] for s in t.stratum}
    t["missing_minutes_release_window"] = rel
    t["missing_minutes_controls"] = [";".join(str(int(v)) for v in by[s]) for s in t.stratum]
    t["missing_minutes_total"] = [int(r + by[s].sum()) for r, s in zip(rel, t.stratum)]
    t["windows_with_missing_minutes"] = [int((r > 0) + (by[s] > 0).sum()) for r, s in zip(rel, t.stratum)]
    return t


def status_tables(cands: pd.DataFrame, comp: pd.DataFrame, up_times_by_key: dict, data_end, gaps=None) -> tuple:
    """Outcomes and the e-process of the complete strata only (rule 4; refusal 8.4). ``up_times_by_key`` maps
    (arm, coin) to that run's up-alarm times; ``gaps`` (E23.data_checks' runs of missing minutes) gives each window's
    missing minutes. Returns (strata rows, path rows, summary rows, pending rows); pending rows carry no outcome.
    Raises if any window of a stratum in the test ends after ``data_end``."""
    end = _utc(data_end)
    span = pd.Timedelta(minutes=WINDOW_MIN)
    t = comp[comp.in_test]
    kept = cands[cands.kept & cands.stratum.isin(t.stratum)]
    ends = list(t.release_utc + span) + list(kept.anchor_utc + span)
    if any(e > end for e in ends):
        raise RuntimeError("a window of a stratum in the test ends after the data end; refused")
    srows, prows, summ = [], [], []
    for (arm, coin), up in up_times_by_key.items():
        # only alarms inside the windows of complete strata are ever looked at
        up = pd.DatetimeIndex(up)
        out = add_e_process(stratum_outcomes(cands, t, up)) if len(t) else t.assign(Z=[], K=[], S=[], p=[], M=[])
        out = add_missing(out, cands, gaps, coin)
        out.insert(0, "coin", coin)
        out.insert(0, "arm", arm)
        srows.append(out)
        for cell in ("fomc+cpi", "fomc", "cpi"):
            sub = out if cell == "fomc+cpi" else add_e_process(out[out.event_type == cell].drop(columns=["M"]))
            th = THRESHOLD if coin == PRIMARY_COIN else THRESHOLD_E25B
            role = ("E25a primary" if (arm, coin, cell) == (PRIMARY_ARM, PRIMARY_COIN, "fomc+cpi") else
                    "E25b secondary" if (arm == PRIMARY_ARM and cell == "fomc+cpi") else "exploratory, no decision")
            j, d = crossing(sub, th)
            nk = sub["n"].to_numpy(int)
            summ.append(dict(arm=arm, coin=coin, cell=cell, role=role, threshold=th, strata_in_test=len(sub),
                             last_stratum_date=str(sub.date.iloc[-1]) if len(sub) else "",
                             releases_with_outcome=int(sub.Z.sum()) if len(sub) else 0,
                             release_rate=float(sub.Z.mean()) if len(sub) else np.nan,
                             control_windows=int(nk.sum()), control_rate=float(sub.K.sum() / nk.sum()) if nk.sum() else np.nan,
                             M=float(sub.M.iloc[-1]) if len(sub) else 1.0, max_M=float(sub.M.max()) if len(sub) else 1.0,
                             windows_with_missing_minutes=int(sub.windows_with_missing_minutes.sum()) if len(sub) else 0,
                             missing_minutes_in_windows=int(sub.missing_minutes_total.sum()) if len(sub) else 0,
                             crossed=j is not None, first_crossing_stratum=j, first_crossing_date=d,
                             rejects_H0=bool(j is not None and role != "exploratory, no decision"),
                             **{f"prod_psi{int(ps)}": float(sub[f"prod_psi{int(ps)}"].iloc[-1]) if len(sub) else 1.0
                                for ps in PSI}))
        prows.append(out[["arm", "coin", "stratum", "date", "event_type", "M", "running_max_M"]])
    pend = comp[~comp.in_test][["stratum", "date", "event_type", "release_utc", "data_needed_until", "complete",
                                "not_complete_because"]]
    return (pd.concat(srows, ignore_index=True) if srows else pd.DataFrame(),
            pd.concat(prows, ignore_index=True) if prows else pd.DataFrame(), pd.DataFrame(summ), pend)


def _git(*args) -> str:
    try:
        res = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True)
        return res.stdout.strip()
    except (OSError, subprocess.CalledProcessError) as err:
        return f"git failed: {err}"


def provenance(cutoff=None, data_end=None) -> str:
    """The git HEAD, the changes to the shared pipeline since the pre-registration's commit, the commits to this
    script, and (review round 11, m11) uncommitted changes to this script and the shared pipeline."""
    git = _git
    own = ("experiments/e25_sequential_events.py",) + PIPELINE_FILES
    lines = [f"run at {pd.Timestamp.now(tz='UTC'):%Y-%m-%d %H:%M} UTC", f"HEAD {git('rev-parse', 'HEAD')}",
             f"cut-off {'' if cutoff is None else f'{_utc(cutoff):%Y-%m-%d %H:%M} UTC'}",
             f"data end (exclusive) {'' if data_end is None else f'{_utc(data_end):%Y-%m-%d %H:%M} UTC'}",
             f"pre-registration commit {PREREG_COMMIT}",
             f"git diff --stat {PREREG_COMMIT} -- {' '.join(PIPELINE_FILES)}:",
             git("diff", "--stat", PREREG_COMMIT, "--", *PIPELINE_FILES) or "(no change)",
             f"git log --oneline {PREREG_COMMIT}..HEAD -- experiments/e25_sequential_events.py:",
             git("log", "--oneline", f"{PREREG_COMMIT}..HEAD", "--", "experiments/e25_sequential_events.py") or "(none)",
             f"git status --porcelain -- {' '.join(own)} (uncommitted changes):",
             git("status", "--porcelain", "--", *own) or "(none)",
             f"git diff --stat HEAD -- {' '.join(own)}:",
             git("diff", "--stat", "HEAD", "--", *own) or "(none)",
             "written after every refusal check of the look passed (pre-download, data end, ledger)"]
    return "\n".join(lines) + "\n"


def ledger_rows(strata: pd.DataFrame, cutoff, data_end, head: str) -> pd.DataFrame:
    """Review round 11, M1 (ii): one row per reported stratum and (arm, coin), as first reported: the release, the kept
    control anchors, n, Z, K, S, the e-values and M, with the look's cut-off, data end and git HEAD. Exact fields as
    strings, so a later look compares them exactly."""
    if not len(strata):
        return pd.DataFrame(columns=["arm", "coin", "stratum", *LEDGER_EXACT, *LEDGER_FLOAT, "first_reported_cutoff",
                                     "first_reported_data_end", "first_reported_head"])
    r = strata[["arm", "coin", "stratum"]].copy()
    r["stratum"] = r["stratum"].astype(int).astype(str)
    r["date"] = strata["date"].astype(str)
    r["event_type"] = strata["event_type"].astype(str)
    r["release_utc"] = [f"{_utc(x):%Y-%m-%d %H:%M}" for x in strata["release_utc"]]
    r["controls"] = strata["controls"].fillna("").astype(str)
    for c in ("n", "Z", "K", "S"):
        r[c] = strata[c].astype(int).astype(str)
    for c in LEDGER_FLOAT:
        r[c] = strata[c].astype(float).map(repr)
    r["first_reported_cutoff"] = f"{_utc(cutoff):%Y-%m-%d %H:%M}"
    r["first_reported_data_end"] = f"{_utc(data_end):%Y-%m-%d %H:%M}"
    r["first_reported_head"] = head
    return r.reset_index(drop=True)


def read_ledger(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path, dtype=str, keep_default_na=False) if path.exists() else None


def ledger_differences(ledger: pd.DataFrame | None, new: pd.DataFrame) -> list[str]:
    """Every stratum already in the ledger must be recomputed identically at a later look: the same release, the same
    kept controls, n, Z, K and S exactly, the e-values and M to 1e-12 (relative). Returns the differences (empty when
    none). A data revision by the mirror, a reschedule into a reported stratum's 30 days, or an edit to the code would
    otherwise change a stratum that has been reported."""
    if ledger is None or not len(ledger):
        return []
    key = ["arm", "coin", "stratum"]
    cur = new.set_index(key)
    diffs = []
    for _, row in ledger.iterrows():
        k = tuple(row[c] for c in key)
        if k not in cur.index:
            diffs.append(f"{k[0]} {k[1]} stratum {k[2]} ({row['date']} {row['event_type']}): reported at the look of "
                         f"{row['first_reported_cutoff']}, not in the test now")
            continue
        now = cur.loc[k]
        bad = [c for c in LEDGER_EXACT if str(now[c]) != str(row[c])]
        bad += [c for c in LEDGER_FLOAT if not np.isclose(float(now[c]), float(row[c]), rtol=1e-12, atol=0.0)]
        if bad:
            diffs.append(f"{k[0]} {k[1]} stratum {k[2]} ({row['date']} {row['event_type']}): "
                         + ", ".join(f"{c} {row[c]} then, {now[c]} now" for c in bad))
    return diffs


def part_status(procs: int, cutoff: str | None):
    """The fresh test at a look. Refuses before any download when no stratum is complete at the cut-off (default: today
    00:00 UTC); refuses again with the data's actual end; refuses, writing nothing, when a stratum already in the ledger
    is recomputed differently (review round 11, M1). Writes only the complete strata, never an outcome of an incomplete
    one, tagged by the cut-off to the minute; the provenance file is written only after every refusal check passed."""
    import urllib.error

    from avmonitor import data as Dm

    if RULES != "amended":
        raise SystemExit("status refused: the fresh test uses the amended rules only (Amendment 1)")
    now = pd.Timestamp.now(tz="UTC")
    cut = _utc(cutoff) if cutoff else now.floor("D")
    ev = future_event_rows()
    why = status_refusals(ev, cut, now)
    if why:
        raise SystemExit("status refused:\n  " + "\n  ".join(why))
    tag = f"{cut:%Y-%m-%d_%H%M}"
    months = pd.period_range(pd.Period(LOAD_START, "M"), pd.Period(cut - pd.Timedelta(minutes=1), "M"))
    paths = {}
    for coin, sym in E12.COINS.items():
        paths[coin] = []
        for m in months:
            dest = Dm.binance_month_path(sym, m)
            try:  # a month cached before it was complete is downloaded again (E24's rule)
                paths[coin].append(Dm.fetch(Dm.binance_month_url(sym, m), dest,
                                            overwrite=not E24._month_complete(dest, m)))
            except urllib.error.HTTPError as err:
                if err.code == 404:
                    break
                raise
    d_end = min(_utc(E23.data_end(paths)), cut)
    why = status_refusals(ev, d_end, now)
    if why:
        raise SystemExit("status refused (data end " + f"{d_end:%Y-%m-%d %H:%M} UTC):\n  " + "\n  ".join(why))
    end_s = f"{d_end:%Y-%m-%d %H:%M}"
    saved = (E23.LOAD_START, E23.MONITOR_START)
    E23.LOAD_START, E23.MONITOR_START = LOAD_START, STREAM_START
    try:
        files, checks, gaps = E23.data_checks(paths, end_s)
    finally:
        E23.LOAD_START, E23.MONITOR_START = saved
    period = dict(events=str(FUTURE_EVENTS), start=LOAD_START, monitor_start=STREAM_START, end=end_s)
    res = _run_many([(c, a, period) for a in STATUS_ARMS for c in E12.COINS], procs)
    cands, comp = fresh_plan(ev, d_end)
    ups = {(r["arm"], r["coin"]): _up_times(r["alarms"]) for r in res}
    del res  # the alarm lists are never written or printed
    strata, path, summ, pend = status_tables(cands, comp, ups, d_end, gaps)
    head = _git("rev-parse", "HEAD")
    new_rows = ledger_rows(strata, cut, d_end, head)
    ledger = read_ledger(OUT / LEDGER)
    diffs = ledger_differences(ledger, new_rows)
    if diffs:
        raise SystemExit("status refused, nothing written: strata already reported are recomputed differently "
                         f"(ledger {OUT / LEDGER}):\n  " + "\n  ".join(diffs))
    prov = provenance(cut, d_end)
    (OUT / f"e25_status_{tag}_provenance.txt").write_text(prov)
    print(prov, flush=True)
    files.to_csv(OUT / f"e25_status_{tag}_data_files.csv", index=False)
    checks.to_csv(OUT / f"e25_status_{tag}_data_checks.csv", index=False)
    gaps.to_csv(OUT / f"e25_status_{tag}_data_gaps.csv", index=False)
    kept = cands[cands.stratum.isin(comp.stratum[comp.in_test])]
    kept.to_csv(OUT / f"e25_status_{tag}_candidates.csv", index=False)
    strata.to_csv(OUT / f"e25_status_{tag}_strata.csv", index=False)
    path.to_csv(OUT / f"e25_status_{tag}_path.csv", index=False)
    summ.assign(data_end=end_s, cutoff=f"{cut:%Y-%m-%d %H:%M}").to_csv(OUT / f"e25_status_{tag}_summary.csv", index=False)
    pend.to_csv(OUT / f"e25_status_{tag}_pending.csv", index=False)
    known = set() if ledger is None else set(map(tuple, ledger[["arm", "coin", "stratum"]].to_numpy()))
    add = new_rows[[tuple(k) not in known for k in new_rows[["arm", "coin", "stratum"]].to_numpy()]]
    out = add if ledger is None else pd.concat([ledger, add], ignore_index=True)
    out.to_csv(OUT / LEDGER, index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(summ.round(4).to_string(index=False))
        print(pend.to_string(index=False))
        print(f"ledger: {len(ledger) if ledger is not None else 0} rows verified, {len(add)} added ({OUT / LEDGER})")


DRYRUN_LOOKS = ("2026-10-10", "2026-10-25")


def part_dryrun(procs: int):
    """Gate 4's dry run on KNOWN data (review round 11, m9): the status part pointed at E23's period (event file
    us_macro_events_2026.csv followed by the forward file's later rows, which Amendment 2's check needs for the strata
    of September 2026; E23's FOMC meetings and F24a's as the schedule; stream from 1 January 2026, data to 24 October
    2026 23:59 UTC already on disk), with the downloader replaced by a guard that only returns files already on disk and
    none later than October 2026, and its outputs in a temporary directory. Two looks (10 and 25 October 2026), so
    the second verifies the ledger of the first. Writes the two looks' summaries to results/e25_dryrun_known_data_
    summary.csv. Not the fresh test: every number is on data seen before the design."""
    import tempfile

    from avmonitor import data as Dm

    def guard(url, dest, overwrite=False):
        dest = Path(dest)
        if not dest.exists() or dest.stem.split("_")[-1] > "2026-10":
            raise SystemExit(f"dry run refused: {dest.name} is not a known month on disk")
        return dest

    g = globals()
    names = ("OUT", "FUTURE_EVENTS", "FOMC_SCHEDULED", "STREAM_START", "LOAD_START")
    saved, fetch = {k: g[k] for k in names}, Dm.fetch
    p23 = period_e23()
    ev = pd.read_csv(p23["events"], dtype=str, keep_default_na=False)
    fwd = future_event_rows()
    ev = pd.concat([ev, fwd[fwd.date > ev.date.max()]], ignore_index=True)[ev.columns]
    sched = tuple(sorted(set(ev.loc[ev.event_type == "fomc", "date"]) | set(FOMC_SCHEDULED)))
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        try:
            ev.to_csv(Path(tmp) / "events.csv", index=False)
            g.update(OUT=Path(tmp), FUTURE_EVENTS=Path(tmp) / "events.csv", STREAM_START=p23["monitor_start"],
                     LOAD_START=p23["start"], FOMC_SCHEDULED=sched)
            Dm.fetch = guard
            for cut in DRYRUN_LOOKS:
                part_status(procs, cut)
                sm = pd.read_csv(Path(tmp) / f"e25_status_{_utc(cut):%Y-%m-%d_%H%M}_summary.csv")
                rows.append(sm.assign(look=cut))
            led = pd.read_csv(Path(tmp) / LEDGER, dtype=str)
        finally:
            g.update(saved)
            Dm.fetch = fetch
    out = pd.concat(rows, ignore_index=True)
    out.insert(0, "label", "dry run of the status part on known data (E23's period), not the fresh test")
    out["ledger_rows_after_look"] = [int((led.first_reported_cutoff <= f"{_utc(c):%Y-%m-%d %H:%M}").sum())
                                     for c in out.look]
    out.to_csv(OUT / "e25_dryrun_known_data_summary.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(out[["look", "arm", "coin", "cell", "role", "strata_in_test", "M", "missing_minutes_in_windows",
                   "ledger_rows_after_look"]].round(4).to_string(index=False))


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
SERIES = E24.SERIES
MARKERS = E24.MARKERS
INK, MUTED, GRID_C, SURFACE = E24.INK, E24.MUTED, E24.GRID_C, E24.SURFACE


def _plt():
    return E24._plt()


def fig_power(plt):
    c = pd.read_csv(OUT / "e25_power_curve.csv")
    pw = pd.read_csv(OUT / "e25_power.csv")
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    for i, s in enumerate(SCENARIOS):
        for proc, ls, lw in (("E25 (p estimated from the stratum)", "-", 2.0), ("oracle (true control rate)", "--", 1.1)):
            h = c[(c.scenario == s) & (c.procedure == proc)]
            ax.plot(pd.to_datetime(h.date), h.P_cross_by_date, color=SERIES[i], ls=ls, lw=lw, drawstyle="steps-post",
                    label=f"{s}" if proc.startswith("E25") else None)
        f = pw[(pw.scenario == s) & pw.procedure.str.startswith("E25")]["F24a_power_at_2027-06-01"].iloc[0]
        ax.scatter([pd.Timestamp(F24A_POWER_AT)], [f], s=46, facecolor=SURFACE, edgecolor=SERIES[i], lw=1.8, zorder=4,
                   marker=MARKERS[i])
    ax.axvline(pd.Timestamp(F24A_POWER_AT), color=MUTED, ls=":", lw=1)
    ax.axhline(ALPHA, color=MUTED, ls=(0, (1, 3)), lw=0.9)
    ax.annotate("F24a's fixed-sample power on 1 June 2027 (open markers)", xy=(pd.Timestamp(F24A_POWER_AT), 0.97),
                xytext=(4, 0), textcoords="offset points", fontsize=7.5, color=MUTED, va="top")
    ax.annotate("0.05", xy=(pd.Timestamp("2026-11-01"), ALPHA), xytext=(0, 3), textcoords="offset points", fontsize=7.5,
                color=MUTED)
    ax.plot([], [], color=MUTED, ls="--", lw=1.1, label="oracle that knows the control rate (descriptive)")
    ax.set_ylim(0, 1)
    ax.grid(color=GRID_C, lw=0.6)
    ax.set_ylabel("P(M has reached 20 by the date)")
    ax.set_xlabel("date (a stratum enters when complete: data to release + 4 weeks + 6 hours, and every event row "
                  "within 30 days settled)")
    ax.legend(frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(1.0, 0.06))
    ax.set_title("E25 (Amendment 1 rules): chance that the e-process has crossed 20, by date [SYNTHETIC; planning inputs REAL S2 BTC rates]",
                 loc="left", fontsize=9.5, color=INK)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "e25_power.png", dpi=130)
    plt.close(fig)


def fig_posthoc(plt):
    import matplotlib.dates as mdates

    st = pd.read_csv(OUT / "e25_posthoc_strata.csv")
    st = st[st.in_test.astype(bool)]
    periods = list(dict.fromkeys(st.period))
    fig, axes = plt.subplots(1, len(periods), figsize=(13.5, 4.8), gridspec_kw=dict(width_ratios=[2.3, 1]), sharey=True)
    for ax, per in zip(np.atleast_1d(axes), periods):
        for i, coin in enumerate(E12.COINS):
            h = st[(st.period == per) & (st.coin == coin)]
            x = pd.to_datetime(h.date)
            ax.step(x, h.M, where="post", color=SERIES[i], lw=1.8, label=coin)
            ax.scatter(x, h.M, color=SERIES[i], s=10, zorder=3, marker=MARKERS[i])
            hit = h[h.M >= THRESHOLD]
            if len(hit):
                ax.scatter([pd.to_datetime(hit.date.iloc[0])], [hit.M.iloc[0]], s=80, facecolor="none",
                           edgecolor=SERIES[i], lw=1.6, zorder=4)
        ax.axhline(THRESHOLD, color=INK, ls="--", lw=1)
        ax.axhline(THRESHOLD_E25B, color=MUTED, ls=":", lw=1)
        ax.axhline(1, color=GRID_C, lw=0.8)
        ax.set_yscale("log")
        ax.grid(color=GRID_C, lw=0.6)
        ax.set_title(per, loc="left", fontsize=9, color=INK)
        ax.set_xlabel("release date (strata in the test, in order)")
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 4, 7, 10)))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    axes[0].set_ylabel("M after each complete stratum (log scale)")
    axes[0].annotate("20 (E25a)", xy=(0.0, THRESHOLD), xycoords=("axes fraction", "data"), xytext=(2, 3),
                     textcoords="offset points", fontsize=7.5, color=INK)
    axes[0].annotate("40 (E25b, ETH and SOL)", xy=(0.0, THRESHOLD_E25B), xycoords=("axes fraction", "data"),
                     xytext=(2, 3), textcoords="offset points", fontsize=7.5, color=MUTED)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.0, 0.9))
    fig.suptitle("E25 post hoc, exploratory [REAL, data seen before the design], Amendment 1 rules: the e-process over FOMC and CPI "
                 "releases, rank monitor S2 (circles: first crossing of 20)", x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "e25_posthoc.png", dpi=130)
    plt.close(fig)


def part_figs():
    plt = _plt()
    for f, need in ((fig_power, "e25_power_curve.csv"), (fig_posthoc, "e25_posthoc_strata.csv")):
        if (OUT / need).exists():
            f(plt)
            print(f"{f.__name__}: done")
        else:
            print(f"{f.__name__}: {need} missing, skipped")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def main():
    part = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else ""
    opts = dict(a.lstrip("-").split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    procs = min(int(opts.get("procs", 2)), 2)
    global RULES
    RULES = opts.get("rules", "amended")
    if RULES not in ("amended", "original"):
        raise SystemExit("--rules must be amended (default) or original")
    t0 = time.time()
    if part == "theory":
        part_theory()
    elif part == "strata":
        part_strata()
    elif part == "gates":
        part_gates(procs)
    elif part == "power":
        part_power()
    elif part == "posthoc":
        part_posthoc(procs)
    elif part == "status":
        part_status(procs, opts.get("cutoff"))
    elif part == "dryrun":
        part_dryrun(procs)
    elif part == "figs":
        part_figs()
    else:
        raise SystemExit(__doc__)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
