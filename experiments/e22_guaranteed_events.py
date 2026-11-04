"""E22. E12's FOMC and CPI test rerun with monitors that carry a proved per-day false alarm bound.

Answers finding J1 of the referee's second report (docs/review/referee_report_paper_round2.md): every real-data
analysis, including E12's pre-specified primary test (FOMC and CPI on BTC, streaming plug-in monitor), used the
plug-in baseline, which has no finite-sample guarantee. The design was fixed in docs/preregistration_e22.md,
committed before this script was run on real data. Only the monitor configuration differs from E12: the
releases, coins, data preparation, outcome window, matched controls and the exact stratified permutation test are
E12's own functions, imported from experiments/e12_crypto_events.py.

REAL MARKET DATA (Binance spot one-minute bars, BTCUSDT, ETHUSDT, SOLUSDT, 2024 and 2025, deseasonalised by the
predictable intraday profile, loaded with download=False through E12.load_coin).

Arms (all streaming over the 711 monitored days; c = 3, calibration 1440 minutes, null band 0.25):
  P0  E12's plug-in, VolatilityMonitor(A = 28,800): reference; must reproduce E12's published numbers exactly
      (the gate: if it does not, the new arms are not run).
  P1  primary: PeriodicVolatilityMonitor(baseline="cs", delta=0.01, A = 36,000, period=1440, offset at 00:00 UTC):
      the confidence-sequence baseline with both statistics restarted at every UTC day boundary. Corollary 1 of
      docs/theory_per_period.md: P(alarm in a given day of a segment) <= delta + 1440 / A = 5%, marginally, if the
      tolerance (H-up), (H-down) of docs/math.md Theorem 6.4 holds to the end of that day.
  S1  VolatilityMonitor(baseline="cs", delta=0.01, A = 36,000), no daily restart: the bound 5% holds for the
      first monitored day of each segment (Theorem 6.4, part 1).
  S2  E14's awake rank monitor (make_monitor("awake-B", ...) of e14_running_tau.py, A = 28,800, payoffs divided
      by 1.25, seed 14000): m / A = 5% for the first day after each restart under exchangeability (Theorem B.2,
      Remark B.3 of docs/theory_running_tau.md). Its release result was published in E14 Q14g.

Outputs: results/e22_reproduction.csv, e22_monitor_summary.csv, e22_windows.csv, e22_tests.csv,
e22_watched.csv, e22_alarms.csv and results/figures/e22_release_rates.png.
Usage: nice -n 10 python experiments/e22_guaranteed_events.py reproduce|run [--procs=2]
  reproduce  P0 only, compared with E12's CSVs (writes e22_reproduction.csv)
  run        P0 (gate), then P1, S1 and S2; tests, families, diagnostics, figure
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import e12_crypto_events as E12  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402
from avmonitor.periodic import PeriodicVolatilityMonitor  # noqa: E402

OUT = ROOT / "results"
FIG = OUT / "figures"
DAY = E12.DAY  # 1440 minutes, one UTC day
DELTA = 0.01
A_CS = DAY / (0.05 - DELTA)  # 36,000, as in E8, E14 and E19
A_RANK = E12.A  # 28,800 for the rank monitor, as in E14 Q14g
RANK_SEED = 14_000  # E14's SEED
ARMS = ("P0", "P1", "S1", "S2")
NEW_ARMS = ("P1", "S1", "S2")
X_ARMS = ("P0", "P1", "S1")  # arms that bet on X = min(u^2, tau^2) and can have a blind increase side
DETECTOR = {"P0": E12.PRIMARY_DETECTOR, "P1": "cs, daily restart", "S1": "cs, no restart", "S2": "awake rank (awake-B)"}
LABEL = {"P0": "P0 plug-in (E12)", "P1": "P1 cs, daily restart", "S1": "S1 cs, no restart", "S2": "S2 awake rank"}
BOUND_TEXT = {"P0": "none proved (nominal m/A = 5% first day of a segment)",
              "P1": "delta + m/A = 5% for every day of a segment, marginal, tolerance to the end of the day",
              "S1": "delta + m/A = 5% for the first monitored day of a segment, tolerance for that day",
              "S2": "m/A = 5% for the first day after each restart, exchangeability"}
# Holm families of the pre-registration (section 6): (family, arm, coins)
FAMILIES = (("F22a", "P1", ("BTC",)), ("F22b", "P1", ("ETH", "SOL")), ("F22c", "S1", ("BTC", "ETH", "SOL")),
            ("F22d", "S2", ("BTC", "ETH", "SOL")), ("F22e", "P0", ("BTC",)))
FAMILY_ROLE = {"F22a": "primary", "F22b": "secondary (P1 replications)", "F22c": "secondary (S1)",
               "F22d": "secondary (S2, published in E14 Q14g)", "F22e": "reference (E12 as published)"}


# ----------------------------------------------------------------------------
# Monitors and per-minute flags
# ----------------------------------------------------------------------------


def utc_day_offset(index: pd.DatetimeIndex) -> int:
    """Offset o such that (t - o) % 1440 == 0 exactly at the 00:00 UTC minutes of a gap-free minute index."""
    idx = index.tz_convert("UTC")
    o = int((-(idx[0].hour * 60 + idx[0].minute)) % DAY)
    pos = np.arange(len(idx))
    midnight = (idx.hour == 0) & (idx.minute == 0)
    if not np.array_equal(midnight, (pos - o) % DAY == 0):
        raise ValueError("the minute index has gaps: day boundaries do not fall every 1440 steps")
    return o


def make_arm(arm: str, offset: int = 0):
    if arm == "P0":
        return VolatilityMonitor(threshold=E12.A, c=E12.C, cal_len=E12.CAL, null_band=E12.BAND)
    if arm == "P1":
        return PeriodicVolatilityMonitor(threshold=A_CS, c=E12.C, cal_len=E12.CAL, null_band=E12.BAND, baseline="cs",
                                         delta=DELTA, period=DAY, offset=offset)
    if arm == "S1":
        return VolatilityMonitor(threshold=A_CS, c=E12.C, cal_len=E12.CAL, null_band=E12.BAND, baseline="cs", delta=DELTA)
    if arm == "S2":
        import e14_running_tau as E14

        return E14.make_monitor("awake-B", E12.BAND, E12.CAL, A_RANK, A_CS, seed=RANK_SEED)
    raise ValueError(arm)


def run_monitor(mon, u: np.ndarray) -> np.ndarray:
    """Feed the stream; return the flag 'the return at t is monitored' (E12.run_streaming's loop)."""
    monitoring = np.zeros(u.size, bool)
    for t in range(u.size):
        monitoring[t] = mon.calibrated
        mon.update(u[t])
    return monitoring


def segment_table(mon, n: int, x_based: bool) -> tuple[pd.DataFrame, np.ndarray]:
    """One row per segment (first and last monitored step, alarm, whether the increase side can alarm) and the
    segment id of every monitored step (-1 while calibrating). The increase side is off when tau^2 <= m_up
    (E19's blind-segment criterion; every increase factor is then at most 1)."""
    alarm_of = {e.monitored_since: e for e in mon.events}
    rows = []
    seg_id = np.full(n, -1, np.int32)
    for k, seg in enumerate(mon.segments):
        s = int(seg["start"])
        e = alarm_of.get(s)
        end = int(e.index) if e is not None else n - 1
        ratio = float(seg["tau2"] / seg["m_up"]) if x_based and seg["m_up"] > 0 else (np.inf if x_based else np.nan)
        rows.append(dict(segment=k, start=s, end=end, minutes=end - s + 1, alarm=e is not None,
                         direction=e.direction if e is not None else "", tau2_over_m_up=ratio,
                         up_active=bool(ratio > 1.0) if x_based else True,
                         down_active=bool(seg.get("down_active", True))))
        seg_id[s : end + 1] = k
    return pd.DataFrame(rows), seg_id


def watched_flags(monitoring: np.ndarray, seg_id: np.ndarray, segs: pd.DataFrame) -> np.ndarray:
    """Monitoring, and the segment's increase side can alarm."""
    up = segs["up_active"].to_numpy(bool) if len(segs) else np.zeros(0, bool)
    w = monitoring.copy()
    m = seg_id >= 0
    w[m] &= up[seg_id[m]]
    return w


def period_counts(monitoring: np.ndarray, seg_id: np.ndarray, offset: int, alarm_index) -> dict:
    """Monitored periods (one segment within one UTC day) and the number of them with an alarm."""
    pos = np.flatnonzero(monitoring)
    day = (pos - offset) // DAY
    keys = seg_id[pos].astype(np.int64) * 100_000 + day
    periods = np.unique(keys).size
    a = np.asarray(alarm_index, np.int64)
    with_alarm = np.unique(seg_id[a].astype(np.int64) * 100_000 + (a - offset) // DAY).size if a.size else 0
    return dict(periods=int(periods), periods_with_alarm=int(with_alarm))


def first_day(segs: pd.DataFrame, events, n: int) -> dict:
    """Segments with at least one day of stream left, and those with an alarm within their first 1440 monitored
    steps (no return is missing in these series, so steps and monitored minutes coincide)."""
    alarm_of = {e.monitored_since: e.index for e in events}
    elig = hit = 0
    for s in segs["start"]:
        if s + DAY > n:
            continue
        elig += 1
        a = alarm_of.get(int(s))
        hit += int(a is not None and a - s < DAY)
    return dict(segments_with_a_day=elig, first_day_alarms=hit)


def long_run_bound_per_monitored_day(arm: str, monitored_minutes: int) -> float:
    """Bound on expected alarms per monitored day under each arm's null: (delta + M / A) / (1 - delta) over M
    monitored minutes for the cs arms (Corollary 6.6; for P1 by the same proof, a restart only lowers W), M / A for
    the rank monitor (Corollary A.5 with delta = 0), nominal M / A for the plug-in."""
    days = monitored_minutes / DAY
    if arm in ("P1", "S1"):
        return (DELTA + monitored_minutes / A_CS) / (1 - DELTA) / days
    return DAY / (A_RANK if arm == "S2" else E12.A)


# ----------------------------------------------------------------------------
# One coin and one arm
# ----------------------------------------------------------------------------


def job(args):
    coin, arm = args
    t_start = time.time()
    events, _, _ = E12.read_events()
    df, u = E12.load_coin(E12.COINS[coin])
    index, n = df.index, u.size
    events["i0"] = index.get_indexer(events["release_utc"])
    if (events["i0"] < 0).any():
        raise ValueError("a release time is not in the price index")
    offset = utc_day_offset(index)
    mon = make_arm(arm, offset)
    monitoring = run_monitor(mon, u)
    ev = E12.alarm_frame(mon.events, 0, u)
    segs, seg_id = segment_table(mon, n, arm in X_ARMS)
    if not np.array_equal(monitoring, seg_id >= 0):
        raise AssertionError("segment table and monitoring flag disagree")
    watched = watched_flags(monitoring, seg_id, segs)
    event_days = events["date"].to_list()
    cache = {}

    def controls(row):
        key = (row["date"].weekday(), row["release_time_et"])
        if key not in cache:
            cache[key] = E12.control_anchors(index, event_days, key[0], key[1], n)
        return cache[key]

    anchors = sorted({i for _, row in events.iterrows() for i in controls(row)})
    det = DETECTOR[arm]
    cstats = {i: E12.window_stats(ev, monitoring, u, i) for i in anchors}
    estats = {int(r["i0"]): E12.window_stats(ev, monitoring, u, int(r["i0"])) for _, r in events.iterrows()}
    tests = E12.stratified_tests(coin, det, "streaming", events, estats, cstats, controls)
    # exploratory sample: releases (and controls) at which the monitor was watching
    w_e = {i: {**s, "monitoring_at_start": bool(watched[i])} for i, s in estats.items()}
    w_c = {i: {**s, "monitoring_at_start": bool(watched[i])} for i, s in cstats.items()}
    wt = E12.stratified_tests(coin, det, "streaming", events, w_e, w_c, controls)
    tests += [dict(r, sample="watched at start") for r in wt if r["sample"] == "monitoring at start"]
    for r in tests:
        r["arm"] = arm

    def ratio_at(i):
        k = seg_id[i]
        return float(segs.loc[k, "tau2_over_m_up"]) if k >= 0 else np.nan

    windows = []
    for _, row in events.iterrows():
        i0 = int(row["i0"])
        st = {k: estats[i0][k] for k in E12.WINDOW_KEYS}
        cs = controls(row)
        cst = [cstats[i] for i in cs]
        cs_up = [c for c in cst if c["up_post"]]
        windows.append(dict(
            arm=arm, coin=coin, detector=det, design="streaming", date=row["date"].date(), event_type=row["event_type"],
            release_utc=row["release_utc"], weekday=row["weekday"], n_controls=len(cs), **st,
            watched_at_start=bool(watched[i0]), watched_share_post=float(watched[i0 : i0 + E12.POST].mean()),
            tau2_over_m_up_at_release=ratio_at(i0),
            control_up_post_rate=float(np.mean([c["up_post"] for c in cst])),
            control_up_post_1h_rate=float(np.mean([c["up_post_1h"] for c in cst])),
            control_up_pre_rate=float(np.mean([c["up_pre"] for c in cst])),
            control_down_follow_rate=float(np.mean([c["down_follow"] for c in cst])),
            control_down_after_up_given_up_rate=float(np.mean([c["down_after_up"] for c in cs_up])) if cs_up else np.nan,
            control_monitoring_at_start_rate=float(np.mean([c["monitoring_at_start"] for c in cst])),
            control_watched_at_start_rate=float(np.mean([watched[i] for i in cs])),
            control_watched_share_post=float(np.mean([watched[i : i + E12.POST].mean() for i in cs]))))
    mon_minutes = int(monitoring.sum())
    pc = period_counts(monitoring, seg_id, offset, ev["index"].to_numpy())
    fd = first_day(segs, mon.events, n)
    summ = dict(arm=arm, coin=coin, detector=det, A=float(mon.A), band=E12.BAND, days=n / DAY, alarms=len(ev),
                up=int((ev.direction == "up").sum()), down=int((ev.direction == "down").sum()),
                alarms_per_day=len(ev) / (n / DAY), monitored_share=float(monitoring.mean()),
                watched_share=float(watched.mean()), monitored_days=mon_minutes / DAY,
                alarms_per_monitored_day=len(ev) / (mon_minutes / DAY),
                long_run_bound_per_monitored_day=long_run_bound_per_monitored_day(arm, mon_minutes),
                **pc, **fd, segments=len(segs), offset=offset,
                period_resets=int(getattr(mon, "period_resets", 0)), seconds=time.time() - t_start)
    if arm in X_ARMS:
        r = segs["tau2_over_m_up"].to_numpy(float)
        off = segs[~segs["up_active"]]
        summ.update(segments_up_off=int(len(off)), up_off_minutes=int(off["minutes"].sum()),
                    up_off_share_of_monitored=float(off["minutes"].sum() / mon_minutes),
                    longest_up_off_days=float(off["minutes"].max() / DAY) if len(off) else 0.0,
                    down_off_minutes=int(segs.loc[~segs["down_active"], "minutes"].sum()),
                    median_tau2_over_m_up=float(np.median(r)), q10_tau2_over_m_up=float(np.quantile(r, 0.1)))
    alarms = ev.assign(arm=arm, coin=coin, time=index[ev["index"].to_numpy()])
    print(f"{coin} {arm}: {len(ev)} alarms, monitored {monitoring.mean():.3f}, watched {watched.mean():.3f}, "
          f"{summ['seconds']:.0f} s", flush=True)
    return tests, windows, summ, alarms


def run_jobs(jobs, procs: int):
    if procs <= 1:
        return [job(j) for j in jobs]
    with Pool(procs) as pool:
        return pool.map(job, jobs, chunksize=1)


# ----------------------------------------------------------------------------
# Reproduction checks
# ----------------------------------------------------------------------------

TEST_COLS = ["events", "events_with", "rate_event", "rate_event_lo", "rate_event_hi", "control_windows",
             "weekday_time_groups", "rate_control_matched", "expected_from_matched_controls", "rate_ratio_matched",
             "expected_under_null", "p_stratified_permutation", "p_stratified_permutation_fewer", "p_stratified_fisher"]
WINDOW_COLS = E12.WINDOW_KEYS + ["n_controls", "control_up_post_rate", "control_up_post_1h_rate", "control_up_pre_rate",
                                 "control_down_follow_rate", "control_down_after_up_given_up_rate",
                                 "control_monitoring_at_start_rate"]


def _compare(mine: pd.DataFrame, pub: pd.DataFrame, keys: list, cols: list, check: str) -> list[dict]:
    m = mine.merge(pub, on=keys, how="outer", suffixes=("_e22", "_pub"), indicator=True)
    rows = [dict(check=check, column="(row keys)", compared=len(m), mismatches=int((m["_merge"] != "both").sum()),
                 max_abs_diff=np.nan)]
    both = m[m["_merge"] == "both"]
    for c in cols:
        a, b = both[f"{c}_e22"], both[f"{c}_pub"]
        if a.dtype == bool or b.dtype == bool or a.dtype == object or b.dtype == object:
            a, b = a.astype(float), b.astype(float)
        a, b = a.to_numpy(float), b.to_numpy(float)
        same = (a == b) | (np.isnan(a) & np.isnan(b))
        diff = np.abs(a - b)
        rows.append(dict(check=check, column=c, compared=len(a), mismatches=int((~same).sum()),
                         max_abs_diff=float(np.nanmax(diff)) if np.isfinite(diff).any() else 0.0))
    return rows


def reproduction(tests: pd.DataFrame, windows: pd.DataFrame, summ: pd.DataFrame) -> pd.DataFrame:
    rd = dict(float_precision="round_trip")
    rows = []
    pub_t = pd.read_csv(OUT / "e12_tests.csv", **rd)
    pub_t = pub_t[(pub_t.detector == E12.PRIMARY_DETECTOR) & (pub_t.design == "streaming")]
    mine_t = tests[(tests.arm == "P0") & tests["sample"].isin(["all", "monitoring at start"])]
    rows += _compare(mine_t, pub_t, ["coin", "event_type", "outcome", "sample"], TEST_COLS, "P0 vs e12_tests.csv")
    prim = mine_t[(mine_t.coin == E12.PRIMARY_COIN) & (mine_t.outcome == "up_post") & (mine_t["sample"] == "all")]
    holm_mine = prim.assign(p_holm_primary=E12.holm(prim["p_stratified_permutation"]))
    pub_prim = pub_t[(pub_t.coin == E12.PRIMARY_COIN) & (pub_t.outcome == "up_post") & (pub_t["sample"] == "all")]
    rows += _compare(holm_mine, pub_prim, ["coin", "event_type", "outcome", "sample"], ["p_holm_primary"],
                     "P0 Holm (primary family) vs e12_tests.csv")
    pub_w = pd.read_csv(OUT / "e12_event_windows.csv", **rd)
    pub_w = pub_w[(pub_w.detector == E12.PRIMARY_DETECTOR) & (pub_w.design == "streaming")].copy()
    mine_w = windows[windows.arm == "P0"].copy()
    for d in (pub_w, mine_w):
        d["date"] = d["date"].astype(str)
    rows += _compare(mine_w, pub_w, ["coin", "date", "event_type"], WINDOW_COLS, "P0 vs e12_event_windows.csv")
    pub_s = pd.read_csv(OUT / "e12_monitor_summary.csv", **rd)
    pub_s = pub_s[pub_s.detector == E12.PRIMARY_DETECTOR]
    rows += _compare(summ[summ.arm == "P0"], pub_s, ["coin"], ["alarms", "up", "down", "monitored_share"],
                     "P0 vs e12_monitor_summary.csv")
    return pd.DataFrame(rows)


def s2_against_e14(tests: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    """S2 against E14 Q14g's published awake-B rows (a difference is reported, not a stop)."""
    pub = pd.read_csv(OUT / "e14g_e12_tests.csv", float_precision="round_trip")
    pub = pub[pub.method == "awake-B"]
    t = tests[(tests.arm == "S2") & (tests.outcome == "up_post") & (tests["sample"] == "all")]
    w = windows[windows.arm == "S2"].groupby(["coin", "event_type"])["monitoring_at_start"].sum().rename("awake_at_release")
    mine = t.set_index(["coin", "event_type"]).join(w).reset_index().rename(
        columns={"events": "releases", "events_with": "up_within_6h", "p_stratified_permutation": "stratified_p"})
    return pd.DataFrame(_compare(mine, pub, ["coin", "event_type"], ["releases", "awake_at_release", "up_within_6h",
                                                                      "stratified_p"], "S2 vs e14g_e12_tests.csv (awake-B)"))


# ----------------------------------------------------------------------------
# Families, watched shares, figure
# ----------------------------------------------------------------------------


def add_families(tests: pd.DataFrame) -> pd.DataFrame:
    t = tests.copy()
    t["family"] = ""
    t["family_role"] = ""
    t["p_holm_family"] = np.nan
    base = (t.outcome == "up_post") & (t["sample"] == "all")
    for fam, arm, coins in FAMILIES:
        sel = base & (t.arm == arm) & t.coin.isin(coins)
        t.loc[sel, "family"] = fam
        t.loc[sel, "family_role"] = FAMILY_ROLE[fam]
        t.loc[sel, "p_holm_family"] = E12.holm(t.loc[sel, "p_stratified_permutation"].to_numpy())
    t["primary_test"] = (t.family == "F22a") & (t.event_type == "fomc")
    t["status"] = np.where(t.family != "", "pre-registered",
                           np.where(t["sample"] == "all", "exploratory (raw p)", "exploratory, subset of releases (raw p)"))
    return t


def watched_table(windows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (arm, coin, etype), g in windows.groupby(["arm", "coin", "event_type"], sort=False):
        w = g["watched_at_start"].astype(bool)
        m = g["monitoring_at_start"].astype(bool)
        up = g["up_post"].astype(bool)
        rows.append(dict(arm=arm, coin=coin, event_type=etype, releases=len(g), monitoring_at_release=int(m.sum()),
                         blind_at_release=int((m & ~w).sum()), watched_at_release=int(w.sum()),
                         watched_share_of_releases=float(w.mean()), mean_watched_share_of_6h=float(g["watched_share_post"].mean()),
                         control_watched_at_start_rate=float(g["control_watched_at_start_rate"].mean()),
                         control_mean_watched_share_of_6h=float(g["control_watched_share_post"].mean()),
                         up_within_6h=int(up.sum()), up_within_6h_when_watched=int((up & w).sum()),
                         median_tau2_over_m_up_at_release=float(g["tau2_over_m_up_at_release"].median())))
    return pd.DataFrame(rows)


def figure(tests: pd.DataFrame, watched: pd.DataFrame):
    """Dot plot: share of release windows and of matched control windows with an up alarm within 6 hours, per arm,
    coin and release type, with the counts, releases watched and p written beside each row."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    ink, muted, grid_c, surface = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
    ev_c, ct_c = "#2a78d6", "#eb6834"  # categorical slots 1 and 2 (validated pair); filled vs open as second cue
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
                         "axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted, "ytick.color": ink,
                         "font.size": 9, "figure.facecolor": surface, "axes.facecolor": surface})
    names = {"fomc": "FOMC statements (16)", "cpi": "CPI releases (22)"}
    t = tests[(tests.outcome == "up_post") & (tests["sample"] == "all")]
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 6.8), sharex=True, sharey=True)
    ys = np.arange(len(ARMS))[::-1].astype(float)
    for r, etype in enumerate(E12.TYPES):
        for c, coin in enumerate(E12.COINS):
            ax = axes[r, c]
            for y, arm in zip(ys, ARMS):
                row = t[(t.arm == arm) & (t.coin == coin) & (t.event_type == etype)].iloc[0]
                wrow = watched[(watched.arm == arm) & (watched.coin == coin) & (watched.event_type == etype)].iloc[0]
                ev, ct = float(row["rate_event"]), float(row["rate_control_matched"])
                ax.plot([ct, ev], [y, y], color=grid_c, lw=2, zorder=1)
                ax.scatter([ct], [y], s=46, facecolor=surface, edgecolor=ct_c, lw=2, zorder=2)
                ax.scatter([ev], [y], s=50, color=ev_c, edgecolor=surface, lw=1.5, zorder=3)
                p = float(row["p_stratified_permutation"])
                lab = f"{int(row['events_with'])}/{int(row['events'])}, watched {int(wrow['watched_at_release'])}\np {p:.2g}"
                if row["family"] == "F22a":
                    lab += f", Holm {float(row['p_holm_family']):.2g}"
                ax.annotate(lab, xy=(1.02, y), xycoords=("axes fraction", "data"), va="center", ha="left", fontsize=7.5, linespacing=1.1,
                            color=ink if arm == "P1" else muted, annotation_clip=False,
                            fontweight="bold" if (arm == "P1" and coin == "BTC") else "normal")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(-0.6, len(ARMS) - 0.4)
            ax.grid(axis="x", color=grid_c, lw=0.6)
            ax.tick_params(axis="y", length=0)
            ax.set_title(f"{coin}: {names[etype]}", loc="left", fontsize=9.5, color=ink)
            if r == 1:
                ax.set_xlabel("share of windows with an up alarm within 6 hours")
    axes[0, 0].set_yticks(ys, [LABEL[a] for a in ARMS])
    handles = [Line2D([], [], color=ev_c, marker="o", ls="", markersize=7, label="release windows"),
               Line2D([], [], color=ct_c, marker="o", markerfacecolor=surface, markeredgewidth=2, ls="", markersize=7,
                      label="matched control windows (same weekday and US Eastern clock time; mean over releases)")]
    fig.legend(handles=handles, loc="lower left", ncol=2, frameon=False, fontsize=8.5, bbox_to_anchor=(0.08, 0.0))
    fig.suptitle("E12's release test rerun with monitors that carry a per-day bound: up alarms within 6 hours, "
                 "streaming design [REAL] Binance, 2024 to 2025", x=0.01, ha="left", fontsize=10.5, color=ink)
    fig.text(0.01, 0.935, "Beside each row: releases with an up alarm / releases, releases watched (monitoring and the "
             "increase side able to alarm), one-sided stratified permutation p (E12's test). P1 is the pre-registered "
             "primary arm; BTC FOMC under P1 is the primary test.", ha="left", fontsize=8, color=muted)
    fig.subplots_adjust(left=0.105, right=0.915, top=0.87, bottom=0.14, wspace=0.55, hspace=0.32)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "e22_release_rates.png", dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def collect(res):
    tests = pd.DataFrame([r for t, _, _, _ in res for r in t])
    windows = pd.DataFrame([w for _, ws, _, _ in res for w in ws])
    summ = pd.DataFrame([s for _, _, s, _ in res])
    alarms = pd.concat([a for _, _, _, a in res], ignore_index=True)
    return tests, windows, summ, alarms


def gate(procs: int):
    res = run_jobs([(c, "P0") for c in E12.COINS], procs)
    tests, windows, summ, alarms = collect(res)
    rep = reproduction(tests, windows, summ)
    bad = rep[rep.mismatches > 0]
    print(rep.to_string(index=False))
    return res, rep, bad


def main():
    part = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "run"
    procs = int(dict(a.lstrip("-").split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a).get("procs", 2))
    procs = min(procs, 2)
    print("REAL DATA. E22: E12's release test with P0 (plug-in, reproduction), P1, S1, S2.", flush=True)
    res0, rep, bad = gate(procs)
    rep.to_csv(OUT / "e22_reproduction.csv", index=False)
    if len(bad):
        raise SystemExit(f"P0 does not reproduce E12 exactly ({len(bad)} columns differ); new arms not run")
    print("P0 reproduces E12 exactly.", flush=True)
    if part == "reproduce":
        return
    res = res0 + run_jobs([(c, a) for a in NEW_ARMS for c in E12.COINS], procs)
    tests, windows, summ, alarms = collect(res)
    rep = pd.concat([rep, s2_against_e14(tests, windows)], ignore_index=True)
    rep.to_csv(OUT / "e22_reproduction.csv", index=False)
    tests = add_families(tests)
    front = ["arm", "family", "family_role", "primary_test", "status"]
    tests = tests[front + [c for c in tests.columns if c not in front]]
    tests.to_csv(OUT / "e22_tests.csv", index=False)
    windows.to_csv(OUT / "e22_windows.csv", index=False)
    summ["bound_text"] = summ["arm"].map(BOUND_TEXT)
    for k, n_col, d_col in (("period", "periods_with_alarm", "periods"), ("first_day", "first_day_alarms", "segments_with_a_day")):
        ci = [E12.wilson(int(a), int(b)) for a, b in zip(summ[n_col], summ[d_col])]
        summ[f"{k}_alarm_share"] = summ[n_col] / summ[d_col]
        summ[f"{k}_alarm_share_lo"] = [x[0] for x in ci]
        summ[f"{k}_alarm_share_hi"] = [x[1] for x in ci]
    summ["per_day_bound"] = 0.05
    summ.to_csv(OUT / "e22_monitor_summary.csv", index=False)
    watched = watched_table(windows)
    watched.to_csv(OUT / "e22_watched.csv", index=False)
    alarms.to_csv(OUT / "e22_alarms.csv", index=False)
    figure(tests, watched)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 200):
        print(rep.to_string(index=False))
        cols = ["arm", "coin", "alarms", "alarms_per_day", "alarms_per_monitored_day", "long_run_bound_per_monitored_day",
                "monitored_share", "watched_share", "period_alarm_share", "first_day_alarm_share", "segments_up_off",
                "up_off_share_of_monitored"]
        print(summ[cols].round(4).to_string(index=False))
        show = ["arm", "family", "coin", "event_type", "events", "events_with", "rate_event", "rate_control_matched",
                "p_stratified_permutation", "p_holm_family", "p_stratified_fisher"]
        print(tests[tests.family != ""][show].round(5).to_string(index=False))
        print(watched.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
