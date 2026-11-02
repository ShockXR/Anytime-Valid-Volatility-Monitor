"""E12. Event study: do crypto volatility alarms line up with scheduled US macro news?

REAL MARKET DATA: Binance spot 1-minute klines (BTCUSDT, ETHUSDT, SOLUSDT) via
``avmonitor.data.load_binance`` (download=False), 1 January 2024 to 31 December
2025 (no missing minutes in the mirror). The first 20 days are a warm-up for the
predictable intraday profile (``intraday_scale`` with its default ``start="past"``
uses past days only, so it is NaN at first); the monitored stream starts on
21 January 2024 (711 days), and returns are passed with their NaN (never filled
with zeros, the monitor skips them). Events inside the warm-up (the CPI release of
11 January 2024) are excluded and listed. Event dates and times are in
data/events/us_macro_events.csv: FOMC rate decisions (statement at 14:00 US
Eastern) and US CPI releases (08:30 US Eastern), each with a source URL found by
web search on 2026-11-02. US Eastern times are converted to UTC with the IANA
time zone database (America/New_York), so daylight saving time is handled by the
library; this script recomputes the UTC column from the Eastern time and stops if
the CSV disagrees. Rows marked cancelled (the October 2025 CPI, a casualty of the
federal government shutdown) or outside the sample are documentation only.

Settings are FROZEN from E5c (chosen on Bitstamp BTC, Jul to Dec 2025) and not
tuned here: A = 28,800 per side (at most a 5% chance of any false alarm per day
when the null holds), c = 3, cal_len = 1440 minutes, null band 0.25, returns
divided by the predictable intraday profile ``intraday_scale`` (20-day lookback).
Benchmark: the Gaussian likelihood ratio monitor (band 0), which reacts to single
large returns because it does not truncate.

Two designs
  "streaming": the monitor runs over the whole two years exactly as in E5c and E6
      (calibrate 1440 minutes, monitor, alarm, recalibrate). Its state at a release
      is whatever the stream left it in; after an alarm it is blind for a day.
  "fresh calibration" (added in iteration 2 after the streaming run showed that the
      monitor was in calibration mode at 40% to 60% of the releases): for every event
      and every control anchor a new monitor is calibrated on the 1440 minutes
      before the anchor and then monitors the 48 hours after it, so both detectors
      are awake at every anchor. For the pre-release window a monitor is calibrated
      on the 1440 minutes ending 6 hours before the anchor and monitors those 6 hours.

Research questions
  Q12a. Is an "up" alarm (volatility increase) more likely in the 6 hours after a
        release than in matched control windows (same US Eastern clock time, same
        weekday, on days at least 2 days away from any event)? Test: a stratified
        permutation test, exact. Within each event's stratum (the event window plus
        its matched control windows) the event label is exchangeable under the null,
        so the number of event windows with an alarm is a sum of independent
        Bernoulli draws with p_i = (alarms in stratum i) / (windows in stratum i);
        the one-sided p-value is P(T >= observed). A pooled one-sided Fisher exact
        test (event windows against the union of control windows) is reported as a
        cruder comparison. Multiple comparisons: the PRIMARY family, fixed before the
        first run, is {FOMC, CPI} x {up alarm within 6 hours} on BTC with the robust
        detector in the streaming design, Holm-corrected at 5%. ETH and SOL are
        replications (raw p-values; the three coins move together, so they are not
        independent evidence). Every other outcome, detector, design and sample is
        secondary and reported raw; a Holm correction over the six robust "up within
        6 hours" tests (3 coins x 2 types) is also given within each design.
  Q12b. How many minutes after the release does the first up alarm come, and how
        often does an up alarm fire in the 6 hours BEFORE the release (which would
        be anticipation or a false alarm, since the number is not yet public)?
  Q12c. Do "down" alarms (volatility normalising) follow in the next 1 to 2 days,
        that is in the window from 6 to 48 hours after the release?
  Q12d. Does the Gaussian LR monitor alarm around releases for the wrong reason
        (one jump rather than a regime change)? Per alarm: the largest absolute
        deseasonalised return in the 5 minutes up to the alarm, in units of the
        baseline scale (above 3 it would have been truncated by the robust monitor),
        and the robust volatility of the next 1440 minutes relative to the baseline
        (near 1 means no persistent change).

Outputs: results/e12_events_checked.csv, e12_events_excluded.csv, e12_event_windows.csv, e12_control_windows.csv,
e12_tests.csv, e12_delay.csv, e12_gaussian_vs_robust.csv, e12_monitor_summary.csv,
e12_all_alarms.csv, e12_profiles.csv and results/figures/e12_event_profiles.png.
No random numbers are used (the permutation distribution is computed exactly).
Usage: python experiments/e12_crypto_events.py   (about 5 minutes, one process)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from avmonitor.data import intraday_scale, load_binance
from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
EVENTS = ROOT / "data/events/us_macro_events.csv"
COINS = {"BTC": "BTCUSDT", "ETH": "ETHUSDT", "SOL": "SOLUSDT"}
START, END = "2024-01-01", "2026-01-01"  # data loaded (the first 20 days warm up the intraday profile)
MONITOR_START = "2024-01-21"  # first minute of the monitored stream; events before it are excluded
TZ = "America/New_York"
A, C, CAL, BAND = 28800.0, 3.0, 1440, 0.25
DETECTORS = {"robust e-detector": (VolatilityMonitor, BAND), "Gaussian LR": (GaussianLRMonitor, 0.0)}
PRIMARY_DETECTOR, PRIMARY_COIN, PRIMARY_DESIGN = "robust e-detector", "BTC", "streaming"
DESIGNS = ("streaming", "fresh calibration")
TYPES = ("fomc", "cpi")
OUTCOMES = ("up_post", "up_pre", "down_follow", "any_post")
POST = 6 * 60  # Q12a window: [release, release + 6 h)
PRE = 6 * 60  # anticipation window: [release - 6 h, release)
DOWN_TO = 48 * 60  # Q12c window: [release + 6 h, release + 48 h)
EXCLUDE_DAYS = 2  # a control day is at least this many days from every event day
PROFILE = (-6 * 60, 12 * 60)  # figure window in minutes
BIN = 5  # minutes per bin in the figure
JUMP_LOOKBACK = 5  # minutes up to the alarm scanned for the largest |return|
DAY = 1440
MIN_FINITE = 1000
WINDOW_KEYS = ["monitoring_at_start", "monitored_share_post", "any_post", "up_post", "n_up_post", "down_post",
               "first_up_delay_min", "first_up_max_abs_recent_over_scale", "first_up_post_over_baseline", "up_pre",
               "last_up_lead_min", "down_follow", "first_down_follow_min", "down_after_up", "abs_ret_first_hour_over_pre",
               "abs_ret_post_over_pre"]


# ----------------------------------------------------------------------------
# Events
# ----------------------------------------------------------------------------


def read_events() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Events used in the study (released or rescheduled, inside the sample), with a UTC check."""
    raw = pd.read_csv(EVENTS, dtype=str, keep_default_na=False)
    raw["date"] = pd.to_datetime(raw["date"])
    used = raw[raw["status"].isin(["released", "rescheduled"])].copy()
    et = pd.to_datetime(used["date"].dt.strftime("%Y-%m-%d") + " " + used["release_time_et"]).dt.tz_localize(TZ)
    used["release_utc"] = et.dt.tz_convert("UTC")
    csv_utc = pd.to_datetime(used["release_time_utc"]).dt.tz_localize("UTC")
    bad = used[(used["release_utc"] != csv_utc).to_numpy()]
    if len(bad):
        raise ValueError(f"release_time_utc disagrees with the time zone conversion:\n{bad[['date', 'release_time_utc']]}")
    used["eastern_zone"] = [t.tzname() for t in et]
    used["weekday"] = used["date"].dt.day_name()
    inside = (used["release_utc"] >= pd.Timestamp(MONITOR_START, tz="UTC")) & (used["release_utc"] < pd.Timestamp(END, tz="UTC"))
    warmup = (used["release_utc"] >= pd.Timestamp(START, tz="UTC")) & ~inside
    excluded = used[warmup][["date", "event_type", "release_time_utc", "description"]].copy()
    excluded["reason"] = f"inside the intraday-profile warm-up ({START} to {MONITOR_START})"
    return used[inside].reset_index(drop=True), raw, excluded


def control_anchors(index: pd.DatetimeIndex, event_days: list, weekday: int, time_et: str, n: int) -> list[int]:
    """Positions of matched control windows: same Eastern clock time and weekday, away from every event day."""
    out = []
    for d in pd.date_range(MONITOR_START, END, freq="D", inclusive="left"):
        if d.weekday() != weekday:
            continue
        if min(abs((d - e).days) for e in event_days) < EXCLUDE_DAYS:
            continue
        t = pd.Timestamp(f"{d:%Y-%m-%d} {time_et}", tz=TZ).tz_convert("UTC")
        i = int(index.get_indexer([t])[0])
        if i < CAL + PRE or i + DOWN_TO > n:
            continue
        out.append(i)
    return out


# ----------------------------------------------------------------------------
# Data and monitors
# ----------------------------------------------------------------------------


def load_coin(symbol: str):
    """Deseasonalised returns from MONITOR_START on; the days before it only warm up the intraday profile."""
    df = load_binance(symbol, START, END, download=False)
    r = df["logret"]
    s = intraday_scale(r, lookback_days=20)  # returns with their NaN; start="past" (default): no future days
    u = (r / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    keep = np.asarray(df.index >= pd.Timestamp(MONITOR_START, tz="UTC"))
    if not np.isfinite(s.to_numpy()[keep]).all():
        raise ValueError("the intraday profile is not available for every monitored minute")
    return df[keep], u[keep]


def alarm_frame(events, offset: int, u: np.ndarray) -> pd.DataFrame:
    """Alarms of a monitor as a frame with diagnostics; ``offset`` maps monitor positions to positions in ``u``."""
    a = np.abs(u)
    rows = []
    for e in events:
        i = e.index + offset
        seg = u[i + 1 : i + 1 + DAY]
        seg = seg[np.isfinite(seg)]
        post = float(np.median(np.abs(seg)) * MAD_NORMAL_SCALE / e.baseline_scale) if seg.size >= MIN_FINITE else np.nan
        recent = a[max(0, i - JUMP_LOOKBACK + 1) : i + 1]
        rows.append(dict(index=i, direction=e.direction, delay=e.index - e.monitored_since, baseline_scale=e.baseline_scale,
                         max_abs_recent_over_scale=float(np.nanmax(recent) / e.baseline_scale), post_over_baseline=post))
    cols = ["index", "direction", "delay", "baseline_scale", "max_abs_recent_over_scale", "post_over_baseline"]
    return pd.DataFrame(rows, columns=cols)


def run_streaming(name: str, u: np.ndarray) -> tuple[pd.DataFrame, np.ndarray]:
    """Run one monitor over the whole stream; return its alarms and a per-minute monitoring flag."""
    cls, band = DETECTORS[name]
    mon = cls(threshold=A, c=C, cal_len=CAL, null_band=band)
    monitoring = np.zeros(u.size, bool)
    for t in range(u.size):
        monitoring[t] = mon.calibrated  # True if the return at t is monitored (not used for calibration)
        mon.update(u[t])
    return alarm_frame(mon.events, 0, u), monitoring


def window_stats(ev: pd.DataFrame, monitoring: np.ndarray, u: np.ndarray, i0: int) -> dict:
    """Alarm facts for the windows anchored at position i0 (a release or a matched control time)."""
    idx = ev["index"].to_numpy()
    d = ev["direction"].to_numpy()
    rel = idx - i0
    post = (rel >= 0) & (rel < POST)
    pre = (rel >= -PRE) & (rel < 0)
    follow = (rel >= POST) & (rel < DOWN_TO)
    up_post, up_pre, dn_follow = post & (d == "up"), pre & (d == "up"), follow & (d == "down")
    first = int(np.argmax(up_post)) if up_post.any() else None
    after_up = (d == "down") & (rel < DOWN_TO) & (rel > (rel[first] if first is not None else np.inf))
    a = np.abs(u)
    pre_mean = float(np.nanmean(a[i0 - PRE : i0]))
    return dict(
        monitoring_at_start=bool(monitoring[i0]),
        monitored_share_post=float(monitoring[i0 : i0 + POST].mean()),
        any_post=bool(post.any()),
        up_post=bool(up_post.any()),
        n_up_post=int(up_post.sum()),
        down_post=bool((post & (d == "down")).any()),
        first_up_delay_min=int(rel[first]) if first is not None else np.nan,
        first_up_max_abs_recent_over_scale=float(ev["max_abs_recent_over_scale"].iloc[first]) if first is not None else np.nan,
        first_up_post_over_baseline=float(ev["post_over_baseline"].iloc[first]) if first is not None else np.nan,
        up_pre=bool(up_pre.any()),
        last_up_lead_min=int(-rel[up_pre].max()) if up_pre.any() else np.nan,
        down_follow=bool(dn_follow.any()),
        first_down_follow_min=int(rel[dn_follow].min()) if dn_follow.any() else np.nan,
        down_after_up=bool(after_up.any()),
        abs_ret_first_hour_over_pre=float(np.nanmean(a[i0 : i0 + 60]) / pre_mean) if pre_mean > 0 else np.nan,
        abs_ret_post_over_pre=float(np.nanmean(a[i0 : i0 + POST]) / pre_mean) if pre_mean > 0 else np.nan,
    )


def fresh_window(name: str, u: np.ndarray, i0: int) -> tuple[dict, pd.DataFrame]:
    """Fresh monitors around anchor i0: one calibrated on [i0 - 1440, i0) monitoring [i0, i0 + 48 h), one calibrated
    on [i0 - 1440 - 6 h, i0 - 6 h) monitoring the 6 hours before i0. Returns window facts and the alarms."""
    cls, band = DETECTORS[name]
    mon = cls(threshold=A, c=C, cal_len=CAL, null_band=band)
    mon.run(u[i0 - CAL : i0 + DOWN_TO])
    ok_post = bool(mon.segments) and mon.segments[0]["start"] == CAL
    post_ev = alarm_frame(mon.events, i0 - CAL, u)
    pre = cls(threshold=A, c=C, cal_len=CAL, null_band=band)
    pre.run(u[i0 - CAL - PRE : i0])
    ok_pre = bool(pre.segments) and pre.segments[0]["start"] == CAL
    pre_ev = alarm_frame(pre.events, i0 - CAL - PRE, u)
    ev = pd.concat([pre_ev, post_ev], ignore_index=True)
    monitoring = np.zeros(u.size, bool)
    monitoring[i0] = ok_post
    st = window_stats(ev, monitoring, u, i0)
    st["calibration_exact"] = ok_post and ok_pre
    return st, ev


# ----------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------


def poisson_binomial_tail(ps, k: int) -> float:
    """P(T >= k) for T a sum of independent Bernoulli(p_i), computed exactly by convolution."""
    dist = np.array([1.0])
    for p in ps:
        dist = np.convolve(dist, [1.0 - p, p])
    return float(dist[k:].sum()) if k < dist.size else 0.0


def holm(pvals) -> np.ndarray:
    p = np.asarray(pvals, float)
    m = p.size
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(np.argsort(p)):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (c - h, c + h)


def stratified_tests(coin, det, design, events, estats, cstats, controls) -> list[dict]:
    rows = []
    samples = ("all", "monitoring at start") if design == "streaming" else ("all",)
    for etype in TYPES:
        sub = events[events["event_type"] == etype]
        for outcome in OUTCOMES:
            for sample in samples:
                ys, ps, union = [], [], {}
                for _, row in sub.iterrows():
                    st = estats[int(row["i0"])]
                    cs = controls(row)
                    if sample == "monitoring at start":
                        if not st["monitoring_at_start"]:
                            continue
                        cs = [i for i in cs if cstats[i]["monitoring_at_start"]]
                    if not cs:
                        continue
                    pool = [st[outcome]] + [cstats[i][outcome] for i in cs]
                    ys.append(bool(st[outcome]))
                    ps.append(float(np.mean(pool)))
                    for i in cs:
                        union[i] = bool(cstats[i][outcome])
                k1, n1 = int(sum(ys)), len(ys)
                k0, n0 = int(sum(union.values())), len(union)
                p_strat = poisson_binomial_tail(ps, k1) if n1 else np.nan
                p_fisher = stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]], alternative="greater")[1] if n1 and n0 else np.nan
                lo, hi = wilson(k1, n1)
                lo0, hi0 = wilson(k0, n0)
                rows.append(dict(coin=coin, detector=det, design=design, event_type=etype, outcome=outcome, sample=sample,
                                 events=n1, events_with=k1, rate_event=k1 / n1 if n1 else np.nan, rate_event_lo=lo, rate_event_hi=hi,
                                 expected_under_null=float(np.sum(ps)), control_windows=n0, control_windows_with=k0,
                                 rate_control=k0 / n0 if n0 else np.nan, rate_control_lo=lo0, rate_control_hi=hi0,
                                 rate_ratio=(k1 / n1) / (k0 / n0) if n1 and n0 and k0 else np.nan,
                                 p_stratified_permutation=p_strat, p_fisher_pooled=p_fisher))
    return rows


# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------


def figure(profiles: dict, alarms: dict, events: pd.DataFrame, n_controls: dict):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    up_c, dn_c, ink, muted, ctrl_c, grid_c = "#eb6834", "#2a78d6", "#0b0b0b", "#52514e", "#6f6e6a", "#e6e5e0"
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": muted,
                         "axes.labelcolor": ink, "xtick.color": muted, "ytick.color": muted, "font.size": 9})
    names = {"fomc": "FOMC statement (14:00 ET)", "cpi": "CPI release (08:30 ET)"}
    lo, hi = PROFILE
    centres = (np.arange(lo, hi, BIN) + BIN / 2) / 60.0
    # rug rows: (design, detector, direction, colour, marker, height factor, legend label)
    rug = [("streaming", "robust e-detector", "up", up_c, "^", 1.46, "robust up alarm, streaming monitor"),
           ("streaming", "robust e-detector", "down", dn_c, "v", 1.36, "robust down alarm, streaming monitor"),
           ("fresh calibration", "robust e-detector", "up", up_c, "o", 1.26, "robust up alarm, fresh calibration at -24 h"),
           ("fresh calibration", "Gaussian LR", "up", muted, "|", 1.16, "Gaussian LR up alarm, fresh calibration at -24 h")]
    fig, axes = plt.subplots(len(TYPES), len(COINS), figsize=(13, 7.4), sharex=True, sharey="row")
    for r, etype in enumerate(TYPES):
        n_ev = int((events["event_type"] == etype).sum())
        ymax = float(max(np.nanmax(profiles[(coin, etype, lab)]) for coin in COINS for lab in ("event", "control")))
        for c, coin in enumerate(COINS):
            ax = axes[r, c]
            ev_prof, ct_prof = profiles[(coin, etype, "event")], profiles[(coin, etype, "control")]
            ax.plot(centres, ct_prof, color=ctrl_c, lw=1.4)
            ax.plot(centres, ev_prof, color=ink, lw=1.6)
            ax.axvline(0, color=muted, lw=0.8, ls=":")
            for design, det, direction, col, mk, h, _ in rug:
                rel = alarms[(coin, det, design, etype)]
                rel = rel[(rel["direction"] == direction)]["rel_min"].to_numpy() / 60.0
                rel = rel[(rel >= lo / 60) & (rel < hi / 60)]
                kw = dict(color=col, linewidths=0.9) if mk == "|" else dict(color=col, linewidths=0.8, edgecolor="white")
                ax.scatter(rel, np.full(rel.size, h * ymax), marker=mk, s=40 if mk == "|" else 28, zorder=3, **kw)
                ax.text(hi / 60 + 0.15, h * ymax, str(rel.size), va="center", ha="left", fontsize=7.5, color=muted)
            ax.set_ylim(0, 1.56 * ymax)
            ax.set_xlim(lo / 60 - 0.2, hi / 60 + 1.0)
            ax.grid(axis="y", color=grid_c, lw=0.6)
            ax.set_title(f"{coin}: {names[etype]}, {n_ev} events", loc="left", fontsize=9.5, color=ink)
            if r == len(TYPES) - 1:
                ax.set_xlabel("hours from the release")
            if c == 0:
                ax.set_ylabel("mean |deseasonalised return|")
    n_fomc, n_cpi = n_controls["fomc"], n_controls["cpi"]
    handles = [Line2D([], [], color=ink, lw=1.6, label="event windows, mean over events"),
               Line2D([], [], color=ctrl_c, lw=1.4, label=f"matched control windows, same weekday and Eastern clock time "
                                                          f"({n_fomc} for FOMC, {n_cpi} for CPI)")]
    for _, _, _, col, mk, _, lab in rug:
        handles.append(Line2D([], [], color=col, marker=mk, ls="", markersize=6 if mk != "|" else 9,
                              markeredgecolor="white" if mk != "|" else col, label=lab))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Average absolute deseasonalised one-minute return around US macro releases, with alarm times "
                 "(counts at the right of each row) [REAL] Binance, 2024 to 2025", x=0.01, ha="left", fontsize=10.5, color=ink)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "e12_event_profiles.png", dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def relative_alarms(ev: pd.DataFrame, events: pd.DataFrame, etype: str) -> pd.DataFrame:
    rels = []
    for _, row in events[events.event_type == etype].iterrows():
        rel = ev["index"].to_numpy() - int(row["i0"])
        keep = (rel >= -DAY) & (rel < DOWN_TO)
        rels.append(pd.DataFrame(dict(date=row["date"].date(), rel_min=rel[keep], direction=ev["direction"].to_numpy()[keep])))
    return pd.concat(rels, ignore_index=True) if rels else pd.DataFrame(columns=["date", "rel_min", "direction"])


def main():
    events, raw, excluded = read_events()
    print(f"REAL DATA. {len(events)} events used ({int((events.event_type == 'fomc').sum())} FOMC, "
          f"{int((events.event_type == 'cpi').sum())} CPI); {len(excluded)} excluded (inside the warm-up): "
          f"{[str(d.date()) for d in excluded['date']]}; {len(raw) - len(events) - len(excluded)} CSV rows are "
          f"documentation only (status: {sorted(set(raw['status']) - {'released', 'rescheduled'})})")
    OUT.mkdir(parents=True, exist_ok=True)
    excluded.to_csv(OUT / "e12_events_excluded.csv", index=False)
    event_days = events["date"].to_list()
    windows, tests, summary, profiles, alarm_rel, n_controls, all_alarms, fresh_ctrl_alarms = [], [], [], {}, {}, {}, [], []
    control_rows = []
    for coin, symbol in COINS.items():
        df, u = load_coin(symbol)
        n = u.size
        index = df.index
        events["i0"] = index.get_indexer(events["release_utc"])
        if (events["i0"] < 0).any():
            raise ValueError("a release time is not in the price index")
        cache = {}

        def controls(row):
            key = (row["date"].weekday(), row["release_time_et"])
            if key not in cache:
                cache[key] = control_anchors(index, event_days, key[0], key[1], n)
            return cache[key]

        events["n_controls"] = [len(controls(r)) for _, r in events.iterrows()]
        anchors = sorted({i for _, row in events.iterrows() for i in controls(row)})
        a = np.abs(u)
        for det in DETECTORS:
            # ---- streaming design
            ev, monitoring = run_streaming(det, u)
            days = n / DAY
            summary.append(dict(coin=coin, detector=det, band=DETECTORS[det][1], days=days, alarms=len(ev),
                                alarms_per_day=len(ev) / days, up=int((ev.direction == "up").sum()),
                                down=int((ev.direction == "down").sum()), monitored_share=float(monitoring.mean())))
            all_alarms.append(ev.assign(coin=coin, detector=det, time=index[ev["index"]]))
            cstats = {i: window_stats(ev, monitoring, u, i) for i in anchors}
            estats = {int(row["i0"]): window_stats(ev, monitoring, u, int(row["i0"])) for _, row in events.iterrows()}
            # ---- fresh calibration design
            fstats = {}
            for i in anchors:
                fstats[i], fev = fresh_window(det, u, i)
                fresh_ctrl_alarms.append(fev.assign(coin=coin, detector=det, anchor=i))
            festats, fresh_ev = {}, []
            for _, row in events.iterrows():
                i0 = int(row["i0"])
                festats[i0], fev = fresh_window(det, u, i0)
                fresh_ev.append(fev)
            fresh_ev = pd.concat(fresh_ev, ignore_index=True)
            n_inexact = sum(not s["calibration_exact"] for s in list(fstats.values()) + list(festats.values()))
            print(f"{coin} {det}: streaming alarms {len(ev)}, monitored share {monitoring.mean():.3f}; fresh windows with an "
                  f"extended calibration: {n_inexact}")
            for design, es, cs_ in (("streaming", estats, cstats), ("fresh calibration", festats, fstats)):
                for _, row in events.iterrows():
                    st = {k: es[int(row["i0"])][k] for k in WINDOW_KEYS}
                    cs = [cs_[i] for i in controls(row)]
                    cs_up = [c for c in cs if c["up_post"]]
                    windows.append(dict(coin=coin, detector=det, design=design, date=row["date"].date(), event_type=row["event_type"],
                                        release_utc=row["release_utc"], weekday=row["weekday"], n_controls=len(cs), **st,
                                        control_up_post_rate=float(np.mean([c["up_post"] for c in cs])),
                                        control_up_pre_rate=float(np.mean([c["up_pre"] for c in cs])),
                                        control_down_follow_rate=float(np.mean([c["down_follow"] for c in cs])),
                                        control_down_after_up_given_up_rate=float(np.mean([c["down_after_up"] for c in cs_up])) if cs_up else np.nan,
                                        control_monitoring_at_start_rate=float(np.mean([c["monitoring_at_start"] for c in cs]))))
                for i, c in cs_.items():
                    control_rows.append(dict(coin=coin, detector=det, design=design, anchor_utc=index[i], **{k: c[k] for k in WINDOW_KEYS}))
                tests += stratified_tests(coin, det, design, events, es, cs_, controls)
            for etype in TYPES:
                alarm_rel[(coin, det, "streaming", etype)] = relative_alarms(ev, events, etype)
                alarm_rel[(coin, det, "fresh calibration", etype)] = relative_alarms(fresh_ev, events, etype)
        # ---- profiles of |u| around releases and around the matched controls (unique anchors per type)
        lo, hi = PROFILE
        for etype in TYPES:
            sub = events[events.event_type == etype]
            ev_mat = np.array([a[i + lo : i + hi] for i in sub["i0"]])
            ctrl_idx = sorted({i for _, row in sub.iterrows() for i in controls(row)})
            ct_mat = np.array([a[i + lo : i + hi] for i in ctrl_idx])
            n_controls[etype] = len(ctrl_idx)
            for label, mat in (("event", ev_mat), ("control", ct_mat)):
                profiles[(coin, etype, label)] = np.nanmean(mat, axis=0).reshape(-1, BIN).mean(axis=1)
    # ---- save
    OUT.mkdir(parents=True, exist_ok=True)
    chk = events[["date", "event_type", "release_time_et", "release_time_utc", "eastern_zone", "weekday", "status",
                  "original_date", "i0", "n_controls", "source_url"]].copy()
    chk["date"] = chk["date"].dt.date
    chk.to_csv(OUT / "e12_events_checked.csv", index=False)
    win = pd.DataFrame(windows)
    win.to_csv(OUT / "e12_event_windows.csv", index=False)
    pd.DataFrame(control_rows).to_csv(OUT / "e12_control_windows.csv", index=False)
    tst = pd.DataFrame(tests)
    tst["family"] = "secondary"
    prim = ((tst.coin == PRIMARY_COIN) & (tst.detector == PRIMARY_DETECTOR) & (tst.design == PRIMARY_DESIGN)
            & (tst.outcome == "up_post") & (tst["sample"] == "all"))
    tst.loc[prim, "family"] = "primary"
    tst["p_holm_primary"] = np.nan
    tst.loc[prim, "p_holm_primary"] = holm(tst.loc[prim, "p_stratified_permutation"])
    tst["p_holm_six_coins_types"] = np.nan
    for design in DESIGNS:
        six = (tst.detector == PRIMARY_DETECTOR) & (tst.design == design) & (tst.outcome == "up_post") & (tst["sample"] == "all")
        tst.loc[six, "p_holm_six_coins_types"] = holm(tst.loc[six, "p_stratified_permutation"])
    tst.to_csv(OUT / "e12_tests.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "e12_monitor_summary.csv", index=False)
    alarms_df = pd.concat(all_alarms, ignore_index=True)
    alarms_df.to_csv(OUT / "e12_all_alarms.csv", index=False)
    # ---- Q12b and Q12c: delay table
    drows = []
    for (coin, det, design, etype), g in win.groupby(["coin", "detector", "design", "event_type"], sort=False):
        dl = g["first_up_delay_min"].dropna()
        drows.append(dict(coin=coin, detector=det, design=design, event_type=etype, events=len(g),
                          monitoring_at_release=int(g["monitoring_at_start"].sum()),
                          control_monitoring_at_start_rate=float(g["control_monitoring_at_start_rate"].mean()),
                          with_up_in_6h=int(g["up_post"].sum()), control_up_in_6h_rate=float(g["control_up_post_rate"].mean()),
                          with_up_in_first_30min=int((dl <= 30).sum()), with_up_in_first_60min=int((dl <= 60).sum()),
                          median_delay_min=float(dl.median()) if len(dl) else np.nan,
                          q25_delay_min=float(dl.quantile(0.25)) if len(dl) else np.nan,
                          q75_delay_min=float(dl.quantile(0.75)) if len(dl) else np.nan,
                          min_delay_min=float(dl.min()) if len(dl) else np.nan, max_delay_min=float(dl.max()) if len(dl) else np.nan,
                          with_up_in_6h_before=int(g["up_pre"].sum()), control_up_in_6h_before_rate=float(g["control_up_pre_rate"].mean()),
                          with_down_6h_to_48h=int(g["down_follow"].sum()), control_down_6h_to_48h_rate=float(g["control_down_follow_rate"].mean()),
                          with_up_then_down=int(g["down_after_up"].sum()),
                          down_after_up_given_up_rate=float(g.loc[g["up_post"], "down_after_up"].mean()) if g["up_post"].any() else np.nan,
                          control_down_after_up_given_up_rate=float(g["control_down_after_up_given_up_rate"].mean()),
                          median_first_down_follow_min=float(g["first_down_follow_min"].median()) if g["down_follow"].any() else np.nan,
                          median_abs_ret_first_hour_over_pre=float(g["abs_ret_first_hour_over_pre"].median()),
                          median_abs_ret_6h_over_pre=float(g["abs_ret_post_over_pre"].median()),
                          median_first_up_max_abs_recent_over_scale=float(g["first_up_max_abs_recent_over_scale"].median()) if g["up_post"].any() else np.nan,
                          median_first_up_post_over_baseline=float(g["first_up_post_over_baseline"].median()) if g["up_post"].any() else np.nan))
    delay = pd.DataFrame(drows)
    delay.to_csv(OUT / "e12_delay.csv", index=False)
    # ---- Q12d: Gaussian LR against the robust detector
    grows = []
    fca = pd.concat(fresh_ctrl_alarms, ignore_index=True)
    for coin in COINS:
        for det in DETECTORS:
            sets = []
            ups = alarms_df[(alarms_df.coin == coin) & (alarms_df.detector == det) & (alarms_df.direction == "up")]
            sets.append(("streaming: all up alarms", ups["max_abs_recent_over_scale"], ups["post_over_baseline"]))
            for design in DESIGNS:
                sub = win[(win.coin == coin) & (win.detector == det) & (win.design == design) & win["up_post"]]
                sets.append((f"{design}: first up alarm within 6 h of a release",
                             sub["first_up_max_abs_recent_over_scale"], sub["first_up_post_over_baseline"]))
            f = fca[(fca.coin == coin) & (fca.detector == det) & (fca.direction == "up")]
            f = f[(f["index"] - f["anchor"] >= 0) & (f["index"] - f["anchor"] < POST)]
            f = f.sort_values("index").groupby("anchor").first()
            sets.append(("fresh calibration: first up alarm within 6 h of a control anchor",
                         f["max_abs_recent_over_scale"], f["post_over_baseline"]))
            for where, m, p in sets:
                m, p = m.astype(float), p.astype(float)
                grows.append(dict(coin=coin, detector=det, alarms=where, n=len(m),
                                  median_max_abs_recent_over_scale=float(m.median()) if len(m) else np.nan,
                                  share_recent_return_above_3_scale=float((m > 3).mean()) if len(m) else np.nan,
                                  share_recent_return_above_5_scale=float((m > 5).mean()) if len(m) else np.nan,
                                  median_post_over_baseline=float(p.median()) if len(p) else np.nan,
                                  share_post_above_1=float((p > 1).mean()) if len(p) else np.nan,
                                  share_post_above_1_25=float((p > 1.25).mean()) if len(p) else np.nan,
                                  share_jump_only=float(((m > 3) & (p <= 1)).mean()) if len(m) else np.nan))
    gvr = pd.DataFrame(grows)
    gvr.to_csv(OUT / "e12_gaussian_vs_robust.csv", index=False)
    prof_rows = []
    lo, hi = PROFILE
    centres = np.arange(lo, hi, BIN) + BIN / 2
    for (coin, etype, label), prof in profiles.items():
        for cmin, v in zip(centres, prof):
            prof_rows.append(dict(coin=coin, event_type=etype, windows=label, minutes_from_release=cmin, mean_abs_u=v))
    pd.DataFrame(prof_rows).to_csv(OUT / "e12_profiles.csv", index=False)
    figure(profiles, alarm_rel, events, n_controls)
    # ---- print
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(pd.DataFrame(summary).round(3).to_string(index=False))
        show = ["coin", "detector", "design", "event_type", "outcome", "sample", "events", "events_with", "rate_event",
                "expected_under_null", "control_windows", "rate_control", "p_stratified_permutation", "p_fisher_pooled", "family",
                "p_holm_primary", "p_holm_six_coins_types"]
        print(tst[tst.outcome == "up_post"][show].round(4).to_string(index=False))
        print(tst[tst.outcome != "up_post"][show].round(4).to_string(index=False))
        print(delay.round(3).to_string(index=False))
        print(gvr.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
