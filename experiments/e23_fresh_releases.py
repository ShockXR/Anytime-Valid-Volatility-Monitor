"""E23. E22's release test on fresh releases: the FOMC statements and CPI releases of 2026.

Pre-registration: docs/preregistration_e23.md, committed on its own (with data/events/us_macro_events_2026.csv) before
any Binance return of 2026 was downloaded or loaded. Every real-data result before E23 on these releases (E12, E14
Q14g, E22) used the releases of 2024 and 2025, whose answer was known when E22 was designed; E23 asks E22's question
on releases nobody in this project had looked at.

Nothing is reimplemented. The monitors, the per-release and per-control window facts, the matched controls and the
exact stratified permutation test are E22's and E12's functions, imported: every (coin, arm) job is E22.job itself.
Only the period and the event list change: ``e12_period`` sets E12's module constants (events file, load start,
monitored start, end) and replaces E12.read_events by ``read_events``, which is E12's reader followed by E23's
inclusion rule (a release counts only if its 6-hour window ends before the data end). The matched controls are
E12.control_anchors over the 2026 period: same weekday and US Eastern clock time, at least 2 days from every event day,
the anchor at least 1440 + 360 minutes after the stream start and its 48-hour follow-up inside the data (so the 6-hour
window of every control also ends before the data end).

REAL MARKET DATA (Binance spot one-minute bars, BTCUSDT, ETHUSDT, SOLUSDT, from the mirror read by
avmonitor.data.load_binance): loaded from 12 December 2025 (20 days that only warm up the intraday profile, as E12's
1 to 21 January 2024), monitored from 1 January 2026 00:00 UTC to the last minute present in the mirror for all three
coins. Missing minutes stay NaN (never zero-filled); the monitors skip them.

Arms (E22's constructors, ``E22.make_arm``): P1 (primary) PeriodicVolatilityMonitor(baseline="cs", delta=0.01,
A = 36,000, period=1440), restart at 00:00 UTC; P0 E12's plug-in; S2 E14's awake rank monitor; S1 P1 without the
restart.

Parts (keep them separate):
  power  power of the two candidate primary designs, BTC FOMC alone and one stratified test over BTC FOMC and CPI
         together, from E22's committed P1 BTC rates and the 2026 release calendar; no market data of 2026
         (writes results/e23_power.csv). Run before the pre-registration; the numbers in it come from here.
  gate   E23's period switch, set to E12's own period and events, reproduces E22's committed P1 BTC rows exactly
         (2024 and 2025 data; writes results/e23_gate.csv). Run by ``run`` first.
  run    gate; downloads and checks the 2026 files; the four arms on the 2026 stream; tests, families, diagnostics,
         figure (results/e23_*.csv, results/figures/e23_release_rates.png).
Usage: nice -n 10 python experiments/e23_fresh_releases.py power|gate|run [--procs=2]
"""

from __future__ import annotations

import sys
import time
import urllib.error
from contextlib import contextmanager
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import e12_crypto_events as E12  # noqa: E402
import e22_guaranteed_events as E22  # noqa: E402
from avmonitor import data as D  # noqa: E402

OUT = ROOT / "results"
FIG = OUT / "figures"
EVENTS_2026 = ROOT / "data/events/us_macro_events_2026.csv"
LOAD_START = "2025-12-12"  # 20 days of intraday-profile warm-up before the stream, as E12's 2024-01-01 to 2024-01-21
MONITOR_START = "2026-01-01"  # first minute of the monitored stream (00:00 UTC, so P1's restart offset is 0)
PLANNED_END = "2026-11-01"  # planning value for the power part only (October 2026 complete); the run uses the data's end
FIRST_MONTH, LAST_MONTH_TRIED = "2025-12", "2026-12"  # monthly files looked for; the search stops at the first missing
ARMS = ("P1", "P0", "S2", "S1")
COINS = tuple(E12.COINS)
PRIMARY_COIN = "BTC"
POOLED = "fomc+cpi"
DESIGNS = ("fomc", POOLED)  # the two candidate primary designs of the pre-registration
PRIMARY_DESIGN = POOLED  # fixed by the power calculation (docs/preregistration_e23.md, section 5)
ALPHA = 0.05
POWER_REPS, POWER_SEED = 20_000, 23_000
# E12's own period, captured at import before anything is changed (used by the gate)
E12_PERIOD = dict(events=str(E12.EVENTS), start=E12.START, monitor_start=E12.MONITOR_START, end=E12.END)
_E12_READ_EVENTS = E12.read_events
CELLS = ("fomc", "cpi")
# Holm families of the pre-registration (section 6): (family, arm, coins, event types or designs)
FAMILIES = (("F23a", "P1", ("BTC",), (PRIMARY_DESIGN,)), ("F23b", "P1", ("ETH", "SOL"), (PRIMARY_DESIGN,)),
            ("F23c", "P1", COINS, CELLS), ("F23d", "P0", COINS, CELLS), ("F23e", "S2", COINS, CELLS),
            ("F23f", "S1", COINS, CELLS))
FAMILY_ROLE = {"F23a": "primary", "F23b": "secondary (P1, primary design on ETH and SOL)",
               "F23c": "secondary (P1, FOMC and CPI separately)", "F23d": "secondary (P0 plug-in)",
               "F23e": "secondary (S2 awake rank)", "F23f": "secondary (S1 cs, no restart)"}


# ----------------------------------------------------------------------------
# The period switch and the inclusion rule
# ----------------------------------------------------------------------------


def read_events():
    """E12's reader (status released or rescheduled, inside [MONITOR_START, END), UTC recomputed and checked), followed
    by E23's inclusion rule: a release counts only if its 6-hour window [t, t + 6 h) ends at or before the data end."""
    used, raw, excluded = _E12_READ_EVENTS()
    end = pd.Timestamp(E12.END, tz="UTC")
    # E12's reader lists a release at or after END with its warm-up reason; say what it is
    after = (pd.to_datetime(excluded["release_time_utc"]).dt.tz_localize("UTC") >= end).to_numpy()
    excluded.loc[after, "reason"] = f"at or after the data end ({E12.END} UTC)"
    late = (used["release_utc"] + pd.Timedelta(minutes=E12.POST) > end).to_numpy()
    if late.any():
        add = used.loc[late, ["date", "event_type", "release_time_utc", "description"]].copy()
        add["reason"] = f"its 6-hour window does not end before the data end ({E12.END} UTC)"
        excluded = pd.concat([excluded, add], ignore_index=True)
    return used[~late].reset_index(drop=True), raw, excluded


def _set_period(events, start, monitor_start, end):
    E12.EVENTS, E12.START, E12.MONITOR_START, E12.END = Path(events), start, monitor_start, end
    E12.read_events = read_events


@contextmanager
def e12_period(events, start, monitor_start, end):
    """Run E12's (and so E22's) functions on another period and event list; restores E12 afterwards."""
    saved = (E12.EVENTS, E12.START, E12.MONITOR_START, E12.END, E12.read_events)
    _set_period(events, start, monitor_start, end)
    try:
        yield
    finally:
        E12.EVENTS, E12.START, E12.MONITOR_START, E12.END, E12.read_events = saved


def run_jobs(jobs, procs: int, period: dict):
    """E22.job for every (coin, arm), with E12 switched to ``period`` (in each worker as well)."""
    if procs <= 1:
        with e12_period(**period):
            return [E22.job(j) for j in jobs]
    with Pool(min(procs, 2), initializer=_set_period, initargs=tuple(period.values())) as pool:
        return pool.map(E22.job, jobs, chunksize=1)


# ----------------------------------------------------------------------------
# The pooled test: one stratified test over FOMC and CPI releases together
# ----------------------------------------------------------------------------


def strata_from_windows(w: pd.DataFrame, outcome: str = "up_post") -> pd.DataFrame:
    """Each release's stratum from E22's window rows: its outcome y, its number of matched control windows n and the
    number k of them with the outcome (k = n times the matched control rate, an integer), and its weekday-and-time
    group (the release type fixes the clock time: FOMC 14:00, CPI 08:30 US Eastern)."""
    n = w["n_controls"].to_numpy(int)
    rate = w[f"control_{outcome}_rate"].to_numpy(float)
    k = np.rint(rate * n).astype(int)
    if not np.allclose(k / n, rate, rtol=0, atol=1e-12):
        raise ValueError("control rates are not k / n")
    return pd.DataFrame(dict(y=w[outcome].astype(bool).to_numpy(), n=n, k=k,
                             group=list(zip(w["weekday"], w["event_type"])), event_type=w["event_type"].to_numpy()))


def stratified_test(st: pd.DataFrame) -> dict:
    """E12's exact stratified permutation test (one-sided, more releases with the outcome) over any set of strata:
    within each release's stratum (its window and its matched control windows) the release label is exchangeable
    under the null, T = number of releases with the outcome is Poisson-binomial with p_i = (y_i + k_i) / (1 + n_i),
    p = P(T >= observed) (E12.poisson_binomial_tail). Cross-check: E12's exact stratified Fisher test, one group per
    weekday and clock time. For one release type this is E12.stratified_tests' row, to the last digit."""
    y, n, k = st["y"].to_numpy(bool), st["n"].to_numpy(int), st["k"].to_numpy(int)
    ps = (y + k) / (1 + n)
    k1, n1 = int(y.sum()), len(y)
    groups = {}
    for g, yi, ni, ki in zip(st["group"], y, n, k):
        rec = groups.setdefault(g, [0, 0, int(ki), int(ni)])
        if rec[2:] != [int(ki), int(ni)]:
            raise ValueError("releases of one weekday-and-time group must share their control windows")
        rec[0] += int(yi)
        rec[1] += 1
    cr = k / n
    lo, hi = E12.wilson(k1, n1)
    rate_matched = float(np.mean(cr)) if n1 else np.nan
    return dict(events=n1, events_with=k1, rate_event=k1 / n1 if n1 else np.nan, rate_event_lo=lo, rate_event_hi=hi,
                weekday_time_groups=len(groups), rate_control_matched=rate_matched,
                expected_from_matched_controls=float(cr.sum()),
                rate_ratio_matched=(k1 / n1) / rate_matched if n1 and rate_matched > 0 else np.nan,
                expected_under_null=float(ps.sum()),
                p_stratified_permutation=E12.poisson_binomial_tail(ps, k1) if n1 else np.nan,
                p_stratified_permutation_fewer=E12.poisson_binomial_tail(ps, k1, lower=True) if n1 else np.nan,
                p_stratified_fisher=E12.stratified_fisher_tail(groups.values()) if n1 else np.nan)


def pooled_rows(windows: pd.DataFrame, outcomes=("up_post", "up_post_1h")) -> pd.DataFrame:
    rows = []
    for (arm, coin), w in windows.groupby(["arm", "coin"], sort=False):
        for outcome in outcomes:
            st = strata_from_windows(w, outcome)
            rows.append(dict(arm=arm, coin=coin, detector=w["detector"].iloc[0], design="streaming", event_type=POOLED,
                             outcome=outcome, sample="all", control_windows=np.nan, **stratified_test(st)))
    return pd.DataFrame(rows)


def check_pooled_against_e12(tests: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    """The pooled function restricted to one release type must give E12.stratified_tests' numbers exactly."""
    rows = []
    cols = ["events", "events_with", "rate_control_matched", "expected_under_null", "p_stratified_permutation",
            "p_stratified_permutation_fewer", "p_stratified_fisher"]
    for (arm, coin, etype), w in windows.groupby(["arm", "coin", "event_type"], sort=False):
        for outcome in ("up_post", "up_post_1h"):
            mine = stratified_test(strata_from_windows(w, outcome))
            ref = tests[(tests.arm == arm) & (tests.coin == coin) & (tests.event_type == etype) & (tests.outcome == outcome)
                        & (tests["sample"] == "all")].iloc[0]
            bad = [c for c in cols if not (mine[c] == ref[c] or (np.isnan(mine[c]) and np.isnan(ref[c])))]
            rows.append(dict(check="pooled function on one type vs E12.stratified_tests", arm=arm, coin=coin,
                             event_type=etype, outcome=outcome, compared=len(cols), mismatches=len(bad),
                             columns_differing=";".join(bad)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Power (planning, no 2026 market data)
# ----------------------------------------------------------------------------


def planning_rates(path: Path = OUT / "e22_tests.csv") -> dict:
    """E22's P1 rates on BTC, 2024 and 2025: release rate and matched control rate per release type."""
    t = pd.read_csv(path)
    t = t[(t.arm == "P1") & (t.coin == PRIMARY_COIN) & (t.outcome == "up_post") & (t["sample"] == "all")]
    return {r.event_type: (float(r.rate_event), float(r.rate_control_matched)) for r in t.itertuples()}


def calendar_strata(events: pd.DataFrame, index: pd.DatetimeIndex) -> pd.DataFrame:
    """Each release's type, weekday-and-time group and number of matched control windows (E12.control_anchors on a
    gap-free minute index of the monitored period; E12.MONITOR_START and E12.END must be set to that period)."""
    event_days = events["date"].to_list()
    rows, cache = [], {}
    for _, r in events.iterrows():
        key = (r["date"].weekday(), r["release_time_et"])
        if key not in cache:
            cache[key] = E12.control_anchors(index, event_days, key[0], key[1], index.size)
        rows.append(dict(date=r["date"].date(), event_type=r["event_type"], group=key, n=len(cache[key])))
    return pd.DataFrame(rows)


def simulate_power(cal: pd.DataFrame, rates: dict, reps: int = POWER_REPS, seed: int = POWER_SEED,
                   alpha: float = ALPHA, designs=DESIGNS) -> dict:
    """Monte Carlo power of E12's stratified permutation test for each design, with common random numbers.

    Model: each release has the outcome independently with probability p1 of its type; each matched control window
    independently with p0 of its type; releases of one weekday-and-time group share their control windows (as in the
    data), so each group's number of control windows with the outcome is Binomial(n, p0), drawn once per replication.
    ``rates`` maps a release type to (p1, p0). Returns, per design, the rejection rate at one-sided level alpha, its
    Monte Carlo standard error, the mean of T and the median p-value."""
    rng = np.random.default_rng(seed)
    groups = sorted(set(cal["group"]))
    gtype = {g: cal.loc[cal["group"] == g, "event_type"].iloc[0] for g in groups}
    gn = {g: int(cal.loc[cal["group"] == g, "n"].iloc[0]) for g in groups}
    gi = np.array([groups.index(g) for g in cal["group"]])
    p1 = np.array([rates[t][0] for t in cal["event_type"]])
    n = cal["n"].to_numpy(int)
    sel = {d: (cal["event_type"].to_numpy() == d) if d != POOLED else np.ones(len(cal), bool) for d in designs}
    rej = {d: 0 for d in designs}
    tsum = {d: 0 for d in designs}
    pv = {d: [] for d in designs}
    for _ in range(reps):
        K = np.array([rng.binomial(gn[g], rates[gtype[g]][1]) for g in groups])
        y = rng.random(len(cal)) < p1
        k = K[gi]
        ps = (y + k) / (1 + n)
        for d in designs:
            s = sel[d]
            T = int(y[s].sum())
            p = E12.poisson_binomial_tail(ps[s], T)
            rej[d] += p <= alpha
            tsum[d] += T
            pv[d].append(p)
    out = {}
    for d in designs:
        pw = rej[d] / reps
        out[d] = dict(power=pw, mc_se=float(np.sqrt(pw * (1 - pw) / reps)), mean_T=tsum[d] / reps,
                      median_p=float(np.median(pv[d])), releases=int(sel[d].sum()))
    return out


def power_part(end: str = PLANNED_END, tag: str = "") -> tuple[pd.DataFrame, pd.DataFrame]:
    """The power table of the pre-registration (section 5), for the release calendar up to ``end`` (exclusive).
    With the default end it writes results/e23_power.csv and e23_power_calendar.csv (the pre-registered numbers); the
    run part calls it again with the data's end and tag "_actual" (pre-registration, section 9: the power at the actual
    count, recomputed with the same code)."""
    rates = planning_rates()
    with e12_period(EVENTS_2026, LOAD_START, MONITOR_START, end):
        events, _, _ = E12.read_events()
        index = pd.date_range(pd.Timestamp(MONITOR_START, tz="UTC"), pd.Timestamp(end, tz="UTC"), freq="1min",
                              inclusive="left")
        cal = calendar_strata(events, index)
    half = {t: (p0 + (p1 - p0) / 2, p0) for t, (p1, p0) in rates.items()}
    null = {t: (p0, p0) for t, (_, p0) in rates.items()}
    fomc_only = {t: (rates[t] if t == "fomc" else null[t]) for t in rates}
    rows = []
    for scen, r in (("E22 rates (P1, BTC, 2024 and 2025)", rates), ("half of E22's excess over the control rate", half),
                    ("E22's FOMC rates, no CPI effect", fomc_only), ("null (release rate = control rate)", null)):
        res = simulate_power(cal, r)
        for d in DESIGNS:
            rows.append(dict(scenario=scen, design=d, **res[d], reps=POWER_REPS, seed=POWER_SEED, alpha=ALPHA,
                             p1_fomc=r["fomc"][0], p0_fomc=r["fomc"][1], p1_cpi=r["cpi"][0], p0_cpi=r["cpi"][1],
                             calendar_end=end))
    pw = pd.DataFrame(rows)
    # reference, observed: the two designs on E22's committed P1 BTC windows (2024 and 2025; known data)
    w = pd.read_csv(OUT / "e22_windows.csv")
    w = w[(w.arm == "P1") & (w.coin == PRIMARY_COIN)]
    ref = []
    for d in DESIGNS:
        sub = w if d == POOLED else w[w.event_type == d]
        t = stratified_test(strata_from_windows(sub))
        ref.append(dict(scenario="observed on E22's 2024 and 2025 releases (known data, reference)", design=d,
                        releases=t["events"], observed_events_with=t["events_with"],
                        observed_rate_control_matched=t["rate_control_matched"],
                        observed_p=t["p_stratified_permutation"]))
    cal_rows = cal.assign(group=cal["group"].map(lambda g: f"{['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'][g[0]]} {g[1]} ET"))
    out = pd.concat([pw, pd.DataFrame(ref)], ignore_index=True)
    out.to_csv(OUT / f"e23_power{tag}.csv", index=False)
    cal_rows.to_csv(OUT / f"e23_power_calendar{tag}.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(cal_rows.to_string(index=False))
        print(out.round(4).to_string(index=False))
    return out, cal


# ----------------------------------------------------------------------------
# Gate: the period switch reproduces E22 on E12's own period
# ----------------------------------------------------------------------------


def gate_part(procs: int = 1) -> pd.DataFrame:
    res = run_jobs([(PRIMARY_COIN, "P1")], 1, E12_PERIOD)
    tests, windows, summ, _ = E22.collect(res)
    rd = dict(float_precision="round_trip")
    pub_t = pd.read_csv(OUT / "e22_tests.csv", **rd)
    pub_t = pub_t[(pub_t.arm == "P1") & (pub_t.coin == PRIMARY_COIN)]
    rows = E22._compare(tests, pub_t, ["coin", "event_type", "outcome", "sample"], E22.TEST_COLS,
                        "P1 BTC through the E23 period switch vs e22_tests.csv")
    pub_w = pd.read_csv(OUT / "e22_windows.csv", **rd)
    pub_w = pub_w[(pub_w.arm == "P1") & (pub_w.coin == PRIMARY_COIN)].copy()
    mine_w = windows.copy()
    for d in (pub_w, mine_w):
        d["date"] = d["date"].astype(str)
    wcols = E22.WINDOW_COLS + ["watched_at_start", "watched_share_post", "control_watched_at_start_rate"]
    rows += E22._compare(mine_w, pub_w, ["coin", "date", "event_type"], wcols, "P1 BTC windows vs e22_windows.csv")
    pub_s = pd.read_csv(OUT / "e22_monitor_summary.csv", **rd)
    pub_s = pub_s[(pub_s.arm == "P1") & (pub_s.coin == PRIMARY_COIN)]
    rows += E22._compare(summ, pub_s, ["coin"], ["alarms", "up", "down", "monitored_share", "periods", "periods_with_alarm"],
                         "P1 BTC summary vs e22_monitor_summary.csv")
    chk = check_pooled_against_e12(tests, windows)
    rows += [dict(check=r["check"] + f" ({r['event_type']}, {r['outcome']})", column="(all)", compared=r["compared"],
                  mismatches=r["mismatches"], max_abs_diff=np.nan) for _, r in chk.iterrows()]
    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "e23_gate.csv", index=False)
    print(rep.to_string(index=False))
    return rep


# ----------------------------------------------------------------------------
# Data: download, end of the stream, checks
# ----------------------------------------------------------------------------


def download_months(symbol: str) -> list[Path]:
    """Monthly files from FIRST_MONTH on, stopping at the first month the mirror does not have."""
    paths = []
    for m in pd.period_range(FIRST_MONTH, LAST_MONTH_TRIED, freq="M"):
        try:
            paths.append(D.fetch(D.binance_month_url(symbol, m), D.binance_month_path(symbol, m)))
        except urllib.error.HTTPError as err:
            if err.code == 404:
                break
            raise
    return paths


def data_end(paths: dict) -> str:
    """One minute after the last bar present for all three coins (the end of the stream, exclusive)."""
    last = []
    for coin, ps in paths.items():
        d = pd.read_csv(ps[-1], usecols=["timestamp"])
        last.append(pd.to_datetime(d["timestamp"].max(), unit="ms", utc=True))
    return (min(last) + pd.Timedelta(minutes=1)).strftime("%Y-%m-%d %H:%M")


def data_checks(paths: dict, end: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    files, checks, gaps = [], [], []
    for coin, ps in paths.items():
        sym = E12.COINS[coin]
        for p in ps:
            with open(p, "rb") as fh:
                nrows = sum(1 for _ in fh) - 1
            files.append(dict(file=str(Path(p).relative_to(ROOT)), bytes=Path(p).stat().st_size, rows=nrows, sha256=D.sha256(p)))
        raw = pd.concat([pd.read_csv(p) for p in ps], ignore_index=True)
        dup_ts = raw["timestamp"].duplicated(keep=False)
        dup_rows_identical = bool(raw[dup_ts].duplicated(keep=False).all()) if dup_ts.any() else True
        df = D.load_binance(sym, LOAD_START, end, download=False)
        c = df.attrs["checks"]
        mon = df.index >= pd.Timestamp(MONITOR_START, tz="UTC")
        miss = df["close"].isna().to_numpy()
        r_nan = df["logret"].isna().to_numpy()
        checks.append(dict(coin=coin, symbol=sym, files=len(ps), first_file=Path(ps[0]).name, last_file=Path(ps[-1]).name,
                           rows_read=c["rows_read"], duplicate_timestamps=c["duplicates"],
                           duplicate_rows_identical=dup_rows_identical,
                           last_bar_in_files=str(pd.to_datetime(raw["timestamp"].max(), unit="ms", utc=True)),
                           stream_end_exclusive=end, missing_minutes_loaded=c["missing_minutes"],
                           missing_minutes_monitored=int(miss[mon].sum()), nonpositive_prices=c["nonpositive_prices"],
                           nan_returns_monitored=int(r_nan[mon].sum()), monitored_minutes=int(mon.sum()),
                           first=c["first"], last=c["last"]))
        # runs of missing minutes
        i = 0
        idx = df.index
        while i < miss.size:
            if miss[i]:
                j = i
                while j < miss.size and miss[j]:
                    j += 1
                gaps.append(dict(coin=coin, first_missing=str(idx[i]), minutes=j - i, in_monitored_period=bool(mon[i])))
                i = j
            else:
                i += 1
    return pd.DataFrame(files), pd.DataFrame(checks), pd.DataFrame(gaps, columns=["coin", "first_missing", "minutes",
                                                                                   "in_monitored_period"])


# ----------------------------------------------------------------------------
# Tests, families, diagnostics, figure
# ----------------------------------------------------------------------------


def add_families(tests: pd.DataFrame) -> pd.DataFrame:
    t = tests.copy()
    t["family"] = ""
    t["family_role"] = ""
    t["p_holm_family"] = np.nan
    base = (t.outcome == "up_post") & (t["sample"] == "all")
    for fam, arm, coins, types in FAMILIES:
        sel = base & (t.arm == arm) & t.coin.isin(coins) & t.event_type.isin(types)
        t.loc[sel, "family"] = fam
        t.loc[sel, "family_role"] = FAMILY_ROLE[fam]
        t.loc[sel, "p_holm_family"] = E12.holm(t.loc[sel, "p_stratified_permutation"].to_numpy())
    t["primary_test"] = t.family == "F23a"
    t["primary_passes"] = t["primary_test"] & (t["p_holm_family"] <= ALPHA)
    t["status"] = np.where(t.family != "", "pre-registered",
                           np.where(t["sample"] == "all", "exploratory (raw p)", "exploratory, subset of releases (raw p)"))
    return t


def per_day_table(summ: pd.DataFrame) -> pd.DataFrame:
    s = summ.copy()
    s["bound_text"] = s["arm"].map(E22.BOUND_TEXT)
    for k, n_col, d_col in (("period", "periods_with_alarm", "periods"), ("first_day", "first_day_alarms", "segments_with_a_day")):
        ci = [E12.wilson(int(a), int(b)) for a, b in zip(s[n_col], s[d_col])]
        s[f"{k}_alarm_share"] = s[n_col] / s[d_col]
        s[f"{k}_alarm_share_lo"] = [x[0] for x in ci]
        s[f"{k}_alarm_share_hi"] = [x[1] for x in ci]
    s["per_day_bound"] = 0.05
    return s


def delay_table(windows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (arm, coin, etype), g in windows.groupby(["arm", "coin", "event_type"], sort=False):
        dl = g["first_up_delay_min"].dropna()
        rows.append(dict(arm=arm, coin=coin, event_type=etype, releases=len(g), with_up_in_6h=int(g["up_post"].sum()),
                         median_delay_min=float(dl.median()) if len(dl) else np.nan,
                         q25_delay_min=float(dl.quantile(0.25)) if len(dl) else np.nan,
                         q75_delay_min=float(dl.quantile(0.75)) if len(dl) else np.nan,
                         min_delay_min=float(dl.min()) if len(dl) else np.nan,
                         max_delay_min=float(dl.max()) if len(dl) else np.nan))
    return pd.DataFrame(rows)


def figure(tests: pd.DataFrame, watched: pd.DataFrame, e22_ref: pd.DataFrame, end: str):
    """Dot plot: share of release windows and of matched control windows with an up alarm within 6 hours, per arm, coin
    and design (pooled, FOMC, CPI), with E22's 2024 to 2025 release rate for the same cell as an open grey diamond."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    ink, muted, grid_c, surface, ref_c = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb", "#8a8984"
    ev_c, ct_c = "#2a78d6", "#eb6834"  # E22's pair (categorical slots 1 and 2); filled vs open as second cue
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
                         "axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted, "ytick.color": ink,
                         "font.size": 9, "figure.facecolor": surface, "axes.facecolor": surface})
    t = tests[(tests.outcome == "up_post") & (tests["sample"] == "all")]
    rows = (POOLED, "fomc", "cpi")
    fig, axes = plt.subplots(3, 3, figsize=(15.5, 9.4), sharex=True, sharey=True)
    ys = np.arange(len(ARMS))[::-1].astype(float)
    for r, etype in enumerate(rows):
        for c, coin in enumerate(COINS):
            ax = axes[r, c]
            for y, arm in zip(ys, ARMS):
                row = t[(t.arm == arm) & (t.coin == coin) & (t.event_type == etype)].iloc[0]
                ev, ct = float(row["rate_event"]), float(row["rate_control_matched"])
                ref = e22_ref[(e22_ref.arm == arm) & (e22_ref.coin == coin) & (e22_ref.event_type == etype)]
                ax.plot([ct, ev], [y, y], color=grid_c, lw=2, zorder=1)
                if len(ref):
                    ax.scatter([float(ref["rate_event"].iloc[0])], [y + 0.22], s=34, marker="D", facecolor=surface,
                               edgecolor=ref_c, lw=1.3, zorder=2)
                ax.scatter([ct], [y], s=46, facecolor=surface, edgecolor=ct_c, lw=2, zorder=2)
                ax.scatter([ev], [y], s=50, color=ev_c, edgecolor=surface, lw=1.5, zorder=3)
                p = float(row["p_stratified_permutation"])
                if etype == POOLED:
                    nw = int(watched[(watched.arm == arm) & (watched.coin == coin)]["watched_at_release"].sum())
                else:
                    nw = int(watched[(watched.arm == arm) & (watched.coin == coin) & (watched.event_type == etype)]["watched_at_release"].iloc[0])
                lab = f"{int(row['events_with'])}/{int(row['events'])}, watched {nw}\np {p:.2g}"
                if row["family"] == "F23a":
                    lab += " (primary)"
                ax.annotate(lab, xy=(1.02, y), xycoords=("axes fraction", "data"), va="center", ha="left", fontsize=7.5,
                            linespacing=1.1, color=ink if arm == "P1" else muted, annotation_clip=False,
                            fontweight="bold" if row["family"] == "F23a" else "normal")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(-0.6, len(ARMS) - 0.25)
            ax.grid(axis="x", color=grid_c, lw=0.6)
            ax.tick_params(axis="y", length=0)
            nrel = int(t[(t.arm == "P1") & (t.coin == coin) & (t.event_type == etype)]["events"].iloc[0])
            name = {POOLED: f"FOMC and CPI together ({nrel})", "fomc": f"FOMC statements ({nrel})", "cpi": f"CPI releases ({nrel})"}[etype]
            ax.set_title(f"{coin}: {name}", loc="left", fontsize=9.5, color=ink)
            if r == len(rows) - 1:
                ax.set_xlabel("share of windows with an up alarm within 6 hours")
    axes[0, 0].set_yticks(ys, [E22.LABEL[a] for a in ARMS])
    handles = [Line2D([], [], color=ev_c, marker="o", ls="", markersize=7, label="2026 release windows"),
               Line2D([], [], color=ct_c, marker="o", markerfacecolor=surface, markeredgewidth=2, ls="", markersize=7,
                      label="2026 matched control windows (same weekday and US Eastern clock time; mean over releases)"),
               Line2D([], [], color=ref_c, marker="D", markerfacecolor=surface, markeredgewidth=1.3, ls="", markersize=6,
                      label="E22: release windows of 2024 and 2025 (same arm and cell)")]
    fig.legend(handles=handles, loc="lower left", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.06, 0.0))
    fig.suptitle(f"E23: E22's release test on fresh releases, up alarms within 6 hours, streaming design [REAL] Binance, "
                 f"1 January 2026 to {pd.Timestamp(end) - pd.Timedelta(minutes=1):%d %B %Y}", x=0.01, ha="left",
                 fontsize=10.5, color=ink)
    fig.text(0.01, 0.948, "Beside each row: releases with an up alarm / releases, releases watched (monitoring and the "
             "increase side able to alarm), one-sided stratified permutation p (E12's test). The primary test (pre-registered) "
             "is P1 on BTC, FOMC and CPI together.", ha="left", fontsize=8, color=muted)
    fig.subplots_adjust(left=0.105, right=0.915, top=0.9, bottom=0.1, wspace=0.55, hspace=0.32)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "e23_release_rates.png", dpi=130)
    plt.close(fig)


def e22_reference() -> pd.DataFrame:
    """E22's committed release rates (2024 and 2025) per arm, coin and cell, with the pooled cell computed from
    E22's committed windows by the same pooled function; for the figure and the report only."""
    t = pd.read_csv(OUT / "e22_tests.csv")
    t = t[(t.outcome == "up_post") & (t["sample"] == "all") & t.arm.isin(ARMS)][["arm", "coin", "event_type", "events",
                                                                                     "events_with", "rate_event",
                                                                                     "rate_control_matched",
                                                                                     "p_stratified_permutation"]]
    w = pd.read_csv(OUT / "e22_windows.csv")
    pr = pooled_rows(w[w.arm.isin(ARMS)], outcomes=("up_post",))[t.columns]
    return pd.concat([t, pr], ignore_index=True)


def blind_segments_job(args) -> list[dict]:
    """Post hoc listing (the pre-registration asks for E19's per-run blind summary, not for the segments): every
    segment of an X-based arm whose increase side cannot alarm (tau^2 <= m_up), with its dates and the releases in it.
    Reruns the arm exactly as E22.job does (same loader, constructor and loop)."""
    coin, arm = args
    events, _, _ = E12.read_events()
    df, u = E12.load_coin(E12.COINS[coin])
    index = df.index
    mon = E22.make_arm(arm, E22.utc_day_offset(index))
    E22.run_monitor(mon, u)
    segs, _ = E22.segment_table(mon, u.size, True)
    rows = []
    for _, sg in segs[~segs["up_active"]].iterrows():
        t0, t1 = index[int(sg["start"])], index[int(sg["end"])]
        inside = events[(events["release_utc"] >= t0) & (events["release_utc"] <= t1)]
        rows.append(dict(arm=arm, coin=coin, segment=int(sg["segment"]), monitoring_from_utc=str(t0), last_minute_utc=str(t1),
                         days=float(sg["minutes"]) / E12.DAY, tau2_over_m_up=float(sg["tau2_over_m_up"]),
                         ended_by=f"{sg['direction']} alarm" if sg["alarm"] else "end of data",
                         releases_inside=";".join(f"{d:%Y-%m-%d} {t}" for d, t in zip(inside["date"], inside["event_type"]))))
    return rows


def run_part(procs: int):
    print("Gate: the period switch on E12's own period must reproduce E22's P1 BTC rows.", flush=True)
    rep = gate_part()
    if (rep.mismatches > 0).any():
        raise SystemExit("the period switch does not reproduce E22 exactly; E23 not run")
    print("REAL DATA. E23: downloading and checking the 2026 files.", flush=True)
    paths = {coin: download_months(sym) for coin, sym in E12.COINS.items()}
    end = data_end(paths)
    files, checks, gaps = data_checks(paths, end)
    files.to_csv(OUT / "e23_data_files.csv", index=False)
    checks.to_csv(OUT / "e23_data_checks.csv", index=False)
    gaps.to_csv(OUT / "e23_data_gaps.csv", index=False)
    print(checks.T.to_string(), flush=True)
    period = dict(events=str(EVENTS_2026), start=LOAD_START, monitor_start=MONITOR_START, end=end)
    with e12_period(**period):
        events, raw, excluded = E12.read_events()
    ev = events[["date", "event_type", "release_time_et", "release_time_utc", "eastern_zone", "weekday", "status",
                 "original_date", "source_url"]].copy()
    ev["date"] = ev["date"].dt.date
    ev.to_csv(OUT / "e23_events_used.csv", index=False)
    excluded.to_csv(OUT / "e23_events_excluded.csv", index=False)
    print(f"{len(events)} releases used ({int((events.event_type == 'fomc').sum())} FOMC, "
          f"{int((events.event_type == 'cpi').sum())} CPI); excluded: {len(excluded)}; stream end {end} UTC", flush=True)
    res = run_jobs([(c, a) for a in ARMS for c in COINS], procs, period)
    tests, windows, summ, alarms = E22.collect(res)
    chk = check_pooled_against_e12(tests, windows)
    chk.to_csv(OUT / "e23_pooled_check.csv", index=False)
    if (chk.mismatches > 0).any():
        raise SystemExit("the pooled function does not reproduce E12's per-type test; stop")
    tests = pd.concat([tests, pooled_rows(windows)], ignore_index=True)
    tests = add_families(tests)
    front = ["arm", "family", "family_role", "primary_test", "primary_passes", "status"]
    tests = tests[front + [c for c in tests.columns if c not in front]]
    tests.to_csv(OUT / "e23_tests.csv", index=False)
    windows.to_csv(OUT / "e23_windows.csv", index=False)
    summ = per_day_table(summ)
    summ.to_csv(OUT / "e23_monitor_summary.csv", index=False)
    watched = E22.watched_table(windows)
    watched.to_csv(OUT / "e23_watched.csv", index=False)
    alarms.to_csv(OUT / "e23_alarms.csv", index=False)
    delay = delay_table(windows)
    delay.to_csv(OUT / "e23_delay.csv", index=False)
    ref = e22_reference()
    ref.to_csv(OUT / "e23_e22_reference.csv", index=False)
    figure(tests, watched, ref, end)
    # pre-registration, section 9: the power at the actual count, same code and seed; its calendar must match the data's
    pw_actual, cal_actual = power_part(end, "_actual")
    w = windows[(windows.arm == "P1") & (windows.coin == PRIMARY_COIN)]
    if list(cal_actual["n"]) != list(w["n_controls"]):
        raise SystemExit("the power calendar's control counts differ from the data's")
    if procs <= 1:
        with e12_period(**period):
            blind = [r for j in [(c, a) for a in E22.X_ARMS for c in COINS] for r in blind_segments_job(j)]
    else:
        with Pool(min(procs, 2), initializer=_set_period, initargs=tuple(period.values())) as pool:
            blind = [r for rows in pool.map(blind_segments_job, [(c, a) for a in E22.X_ARMS for c in COINS], chunksize=1)
                     for r in rows]
    blind = pd.DataFrame(blind, columns=["arm", "coin", "segment", "monitoring_from_utc", "last_minute_utc", "days",
                                         "tau2_over_m_up", "ended_by", "releases_inside"])
    blind.to_csv(OUT / "e23_blind_segments.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 300):
        cols = ["arm", "coin", "alarms", "alarms_per_day", "alarms_per_monitored_day", "long_run_bound_per_monitored_day",
                "monitored_share", "watched_share", "period_alarm_share", "first_day_alarm_share", "segments_up_off",
                "up_off_share_of_monitored", "offset", "period_resets"]
        print(summ[cols].round(4).to_string(index=False))
        show = ["arm", "family", "coin", "event_type", "events", "events_with", "rate_event", "rate_control_matched",
                "p_stratified_permutation", "p_holm_family", "p_stratified_fisher"]
        print(tests[tests.family != ""][show].round(5).to_string(index=False))
        print(tests[(tests.family == "") & (tests.outcome == "up_post") & (tests["sample"] == "all")][show].round(5).to_string(index=False))
        print(watched.round(3).to_string(index=False))
        print(delay.to_string(index=False))
        print(blind.to_string(index=False))
        print(pw_actual.round(4).to_string(index=False))


def main():
    part = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "run"
    procs = int(dict(a.lstrip("-").split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a).get("procs", 2))
    procs = min(procs, 2)
    t0 = time.time()
    if part == "power":
        power_part()
    elif part == "gate":
        gate_part()
    elif part == "run":
        run_part(procs)
    else:
        raise SystemExit("usage: e23_fresh_releases.py power|gate|run [--procs=2]")
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
