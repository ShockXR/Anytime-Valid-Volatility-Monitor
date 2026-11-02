"""E7d. Event study for the Nifty monitor with verified event dates.

REAL MARKET DATA (Nifty 50 1-minute bars, see e07_nifty.py). Event dates are in
data/events/nifty_events.csv, each with a source URL found by web search on
2026-11-02. Dates that could not be verified were left out (see the report).

Research question
  Q7d. Do alarms of the streaming monitor (pre-registered primary configuration:
       calibration 750 steps, null band 0.25, A = 7500) occur more often on event
       sessions (Union Budget days, RBI monetary policy days, general election result
       days, the March 2020 COVID crash window, monthly index derivative expiry days)
       than on ordinary sessions (Fisher exact test)? At what time of day do they fire?

Event definitions
  * budget, rbi_policy, election_result: dates from the CSV.
  * covid_window: every session from 2 to 31 March 2020 (our own calendar definition of
    "the March 2020 crash"; it contains the two lower-circuit days in the CSV).
  * expiry: last Thursday of each month up to August 2025 and last Tuesday from
    September 2025 (rule and sources in data/events/nifty_expiry_rule.csv). If that
    weekday is not a session in the data (a holiday), the month is left out rather
    than moved to another day, because the holiday rule was not verified.
  Ordinary sessions are sessions that belong to none of these types.

Usage: python experiments/e07d_events.py   (about 1 minute)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from e07_nifty import A, C, OUT, PRIMARY, ROOT, deseasonalise, load_clean, robust_sd, wilson  # noqa: E402

from avmonitor.monitor import VolatilityMonitor  # noqa: E402

FIG = OUT / "figures"
EVENTS = ROOT / "data/events/nifty_events.csv"
TYPES = ["budget", "rbi_policy", "election_result", "covid_window", "expiry"]
BINS = [555, 600, 660, 720, 780, 840, 930]
BIN_LABELS = ["09:15-10:00", "10:00-11:00", "11:00-12:00", "12:00-13:00", "13:00-14:00", "14:00-15:30"]


def session_table():
    """Run the primary monitor once and summarise it per session."""
    df, info, _ = load_clean()
    u = deseasonalise(df).to_numpy()
    codes, sessions = pd.factorize(df["session"])
    cal, band = PRIMARY
    mon = VolatilityMonitor(threshold=A, c=C, cal_len=cal, null_band=band)
    monitoring = np.zeros(u.size, bool)
    for t, x in enumerate(u):
        mon.update(x)
        monitoring[t] = mon.calibrated
    ev = pd.DataFrame([dict(index=e.index, direction=e.direction) for e in mon.events])
    ev["session"] = sessions[codes[ev["index"]]]
    ev["minute"] = df["minute"].to_numpy()[ev["index"]]
    ev["time"] = df.index[ev["index"]]
    s = pd.DataFrame({"session": sessions})
    s["date"] = s["session"].dt.date
    s["monitored_share"] = pd.Series(monitoring).groupby(codes).mean().to_numpy()
    # monitoring mode at the first bar: otherwise the monitor is still recalibrating after an
    # earlier alarm and is blind for part of the session (the share above is endogenous,
    # because an alarm inside the session also starts a recalibration)
    s["monitoring_at_open"] = pd.Series(monitoring).groupby(codes).first().to_numpy()
    s["robust_vol"] = pd.Series(u).groupby(codes).apply(lambda x: robust_sd(x.to_numpy())).to_numpy()
    # level relative to the median of the previous 20 sessions (is the event day itself volatile?)
    s["vol_vs_prev20"] = s["robust_vol"] / s["robust_vol"].rolling(20, min_periods=10).median().shift(1)
    s["close"] = df.groupby("session")["close"].last().to_numpy()
    s["prev_close"] = s["close"].shift(1)
    for d in ("up", "down"):
        cnt = ev[ev.direction == d].groupby("session").size()
        s[f"n_{d}"] = s["session"].map(cnt).fillna(0).astype(int)
    s["n_alarms"] = s["n_up"] + s["n_down"]
    return df, s, ev


def tag_events(s: pd.DataFrame):
    e = pd.read_csv(EVENTS, parse_dates=["date"])
    dates = set(s["date"])
    rows = {t: set() for t in TYPES}
    missing = []
    for _, r in e.iterrows():
        t = "covid_window" if r.event_type == "covid_halt" else r.event_type
        if r.date.date() in dates:
            rows[t].add(r.date.date())
        else:
            missing.append((str(r.date.date()), r.event_type))
    covid = s[(s["session"] >= pd.Timestamp("2020-03-02", tz=s["session"].dt.tz))
              & (s["session"] <= pd.Timestamp("2020-03-31", tz=s["session"].dt.tz))]["date"]
    rows["covid_window"] |= set(covid)
    # monthly expiry from the verified weekday rule
    left_out = []
    for p in pd.period_range("2017-04", "2026-10", freq="M"):
        days = pd.date_range(p.start_time, p.end_time.normalize(), freq="D")
        wd = 3 if p <= pd.Period("2025-08", "M") else 1  # Thursday = 3, Tuesday = 1
        last = max(d for d in days if d.weekday() == wd).date()
        if last in dates:
            rows["expiry"].add(last)
        else:
            left_out.append(str(last))
    for t in TYPES:
        s[t] = s["date"].isin(rows[t])
    s["ordinary"] = ~s[TYPES].any(axis=1)
    return s, missing, left_out


def fisher_rows(s: pd.DataFrame, open_only: bool = False):
    out = []
    base = s[s["monitoring_at_open"]] if open_only else s
    o = base[base["ordinary"]]
    for t in TYPES:
        g = base[base[t]]
        for what, col in (("any alarm", "n_alarms"), ("up alarm", "n_up"), ("down alarm", "n_down")):
            k1, n1 = int((g[col] > 0).sum()), len(g)
            k0, n0 = int((o[col] > 0).sum()), len(o)
            odds, p = stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]]) if n1 else (np.nan, np.nan)
            lo, hi = wilson(k1, n1)
            lo0, hi0 = wilson(k0, n0)
            out.append(dict(event_type=t, outcome=what, sessions_used="monitoring at open" if open_only else "all", event_sessions=n1,
                            event_sessions_with=k1, rate_event=k1 / n1 if n1 else np.nan, rate_event_lo=lo,
                            rate_event_hi=hi, ordinary_sessions=n0, ordinary_sessions_with=k0, rate_ordinary=k0 / n0,
                            rate_ordinary_lo=lo0, rate_ordinary_hi=hi0, odds_ratio=odds, fisher_p=p,
                            median_vol_vs_prev20_event=float(g["vol_vs_prev20"].median()),
                            median_vol_vs_prev20_ordinary=float(o["vol_vs_prev20"].median())))
    return out


def time_of_day(s: pd.DataFrame, ev: pd.DataFrame) -> pd.DataFrame:
    ev = ev.copy()
    ev["date"] = ev["session"].dt.date
    lab = {}
    for t in TYPES + ["ordinary"]:
        for d in s.loc[s[t], "date"]:
            lab.setdefault(d, t)  # first matching type, TYPES order
    ev["group"] = ev["date"].map(lab).fillna("ordinary")
    ev["bin"] = pd.cut(ev["minute"], BINS, labels=BIN_LABELS, right=False)
    tab = pd.crosstab(ev["group"], ev["bin"], normalize="index").round(3)
    tab["alarms"] = ev.groupby("group").size()
    tab["median_minute"] = ev.groupby("group")["minute"].median()
    tab["median_time"] = tab["median_minute"].map(lambda m: f"{int(m // 60):02d}:{int(m % 60):02d}")
    return tab.reset_index(), ev


def figures(df: pd.DataFrame, s: pd.DataFrame, ev: pd.DataFrame, fis: pd.DataFrame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    up_c, dn_c, ink, muted = "#eb6834", "#2a78d6", "#0b0b0b", "#52514e"
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": muted,
                         "axes.labelcolor": ink, "xtick.color": muted, "ytick.color": muted, "font.size": 9})
    # Figure 1: price with alarms in two crisis windows
    wins = [("2020-02-17", "2020-04-30", "Feb to Apr 2020 (COVID crash; lower circuits 13 and 23 Mar)"),
            ("2024-05-27", "2024-06-14", "27 May to 14 Jun 2024 (election result day 4 Jun)")]
    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    for ax, (a, b, title) in zip(axes, wins):
        w = df.loc[a:b]
        x = np.arange(len(w))  # trading-minute axis, no overnight or weekend gaps
        ax.plot(x, w["close"].to_numpy(), lw=0.8, color=muted)
        pos = pd.Series(x, index=w.index)
        e = ev[(ev["time"] >= w.index[0]) & (ev["time"] <= w.index[-1])]
        for d, col, mk in (("up", up_c, "^"), ("down", dn_c, "v")):
            ed = e[e.direction == d]
            ax.scatter(pos[ed["time"]].to_numpy(), w.loc[ed["time"], "close"].to_numpy(), s=46, marker=mk, color=col,
                       edgecolor="white", linewidth=0.8, zorder=3, label=f"{d} alarm ({len(ed)})")
        sess = w["session"].drop_duplicates()
        step = max(1, len(sess) // 10)
        ticks = [pos[w.index[w["session"] == t][0]] for t in sess.iloc[::step]]
        ax.set_xticks(ticks, [t.strftime("%d %b") for t in sess.iloc[::step]])
        ax.set_ylabel("Nifty 50 index")
        ax.set_title(title, loc="left", fontsize=10, color=ink)
        ax.grid(axis="y", color="#e6e5e0", lw=0.6)
        ax.legend(frameon=False, loc="lower left")
    fig.suptitle("Monitor alarms on Nifty 50 one-minute data [REAL] (calibration 750, band 0.25, A = 7500)",
                 x=0.01, ha="left", fontsize=11, color=ink)
    fig.tight_layout()
    fig.savefig(FIG / "e07_nifty_alarms_crises.png", dpi=130)
    plt.close(fig)

    # Figure 2: alarm rate by event type with Wilson intervals
    f = fis[(fis["outcome"] == "any alarm") & (fis["sessions_used"] == "all")].copy()
    names = {"budget": "Union Budget", "rbi_policy": "RBI policy", "election_result": "Election result",
             "covid_window": "COVID window (Mar 2020)", "expiry": "Monthly expiry"}
    labels = [f"{names[t]} (n = {n})" for t, n in zip(f["event_type"], f["event_sessions"])]
    o = f.iloc[0]
    labels.append(f"Ordinary sessions (n = {int(o.ordinary_sessions)})")
    rate = list(f["rate_event"]) + [o.rate_ordinary]
    lo = list(f["rate_event_lo"]) + [o.rate_ordinary_lo]
    hi = list(f["rate_event_hi"]) + [o.rate_ordinary_hi]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    y = np.arange(len(labels))[::-1]
    cols = [up_c] * len(f) + [muted]
    ax.hlines(y, lo, hi, color=cols, lw=2)
    ax.scatter(rate, y, color=cols, s=50, zorder=3, edgecolor="white", linewidth=0.8)
    ax.axvline(o.rate_ordinary, color=muted, lw=0.8, ls=":")
    for yi, r, p in zip(y, rate, list(f["fisher_p"]) + [np.nan]):
        txt = f"{r:.0%}" + (f"  (Fisher p = {p:.2g})" if np.isfinite(p) else "")
        ax.text(1.03, yi, txt, va="center", fontsize=8.5, color=ink, transform=ax.get_yaxis_transform())
    ax.set_yticks(y, labels)
    ax.set_xlim(-0.01, 1.02)
    ax.set_xlabel("share of sessions with at least one alarm (dot) and Wilson 95% interval")
    ax.set_title("Alarm rate on event sessions vs ordinary sessions [REAL]", loc="left", fontsize=10, color=ink)
    ax.grid(axis="x", color="#e6e5e0", lw=0.6)
    fig.tight_layout()
    fig.savefig(FIG / "e07_nifty_event_alarm_rates.png", dpi=130)
    plt.close(fig)


def main():
    df, s, ev = session_table()
    s, missing, left_out = tag_events(s)
    print(f"REAL DATA: {len(s)} sessions, {len(ev)} alarms (primary config {PRIMARY})")
    print("event dates not in the session list:", missing)
    print("expiry months left out (rule date not a session):", len(left_out), left_out)
    fis = pd.DataFrame(fisher_rows(s, False) + fisher_rows(s, True))
    fis.to_csv(OUT / "e07d_event_alarm_rates.csv", index=False)
    show = ["event_type", "outcome", "sessions_used", "event_sessions", "event_sessions_with", "rate_event",
            "rate_ordinary", "odds_ratio", "fisher_p", "median_vol_vs_prev20_event", "median_vol_vs_prev20_ordinary"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(fis[show].round(3).to_string(index=False))
    tod, ev2 = time_of_day(s, ev)
    tod.to_csv(OUT / "e07d_time_of_day.csv", index=False)
    print(tod.to_string(index=False))
    # per-event detail for the dated events
    det = []
    for t in ("budget", "rbi_policy", "election_result", "covid_window"):
        for _, r in s[s[t]].iterrows():
            e = ev2[ev2["date"] == r["date"]]
            det.append(dict(date=r["date"], event_type=t, alarms=int(r["n_alarms"]), up=int(r["n_up"]),
                            down=int(r["n_down"]), alarm_times=" ".join(f"{x:%H:%M}{d[0]}" for x, d in zip(e["time"], e["direction"])),
                            monitored_share=round(r["monitored_share"], 3), monitoring_at_open=bool(r["monitoring_at_open"]), vol_vs_prev20=round(r["vol_vs_prev20"], 3),
                            close_to_close_return=round(np.log(r["close"] / r["prev_close"]), 4)))
    det = pd.DataFrame(det)
    det.to_csv(OUT / "e07d_event_detail.csv", index=False)
    print(det[det.event_type != "covid_window"].to_string(index=False))
    pd.DataFrame({"left_out_expiry_rule_dates": left_out}).to_csv(OUT / "e07d_expiry_left_out.csv", index=False)
    figures(df, s, ev, fis)


if __name__ == "__main__":
    main()
