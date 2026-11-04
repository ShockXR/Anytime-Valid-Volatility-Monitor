"""W3: measure the look-ahead at experiments/e07_nifty.py line 150.

Builds the E7 monitor input three ways and reruns every E7 Q7b configuration on each:
  committed      session_profile as committed (line 150 fills the first 10 sessions with the
                 median of the first 20 sessions, which uses the future).
  no_lookahead   the first 10 sessions get no profile (NaN, so the monitor skips them); any
                 later cell that line 150 fills (closing minutes from August 2026, which use
                 past data only) is filled as committed.  This isolates the look-ahead.
  drop_line150   line 150 removed entirely (the minimal fix): every cell without a past-only
                 median stays NaN.

Nothing under results/ is touched.  Outputs go to reports/round4/W3/.
Run:  python reports/round4/W3/w3_e07_lookahead.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments"))
import e07_nifty as E  # noqa: E402
import e07d_events as D  # noqa: E402  (tag_events, fisher_rows; importing writes nothing)

HERE = Path(__file__).resolve().parent
LB = E.LOOKBACK
MINP = max(5, LB // 2)
VARIANTS = ("committed", "no_lookahead", "drop_line150")


def profile_variant(df: pd.DataFrame, mode: str) -> pd.Series:
    """Copy of e07_nifty.session_profile with line 150 switchable."""
    a = df["logret"].abs()
    wide = pd.DataFrame({"a": a.values, "s": df["session"].values, "m": df["minute"].values}).pivot_table(
        index="s", columns="m", values="a", aggfunc="first", dropna=False)
    prof = wide.rolling(LB, min_periods=MINP).median().shift(1)
    fill = wide.iloc[:LB].median()
    if mode == "committed":
        prof = prof.fillna(fill)
    elif mode == "no_lookahead":
        later = prof.iloc[LB // 2:].fillna(fill)
        prof = pd.concat([prof.iloc[:LB // 2], later])
    elif mode != "drop_line150":
        raise ValueError(mode)
    if E.SMOOTH > 1:
        prof = prof.T.rolling(E.SMOOTH, center=True, min_periods=1).mean().T
    prof = prof.div(prof.median(axis=1), axis=0)
    long = prof.stack(future_stack=True)
    idx = pd.MultiIndex.from_arrays([df["session"].values, df["minute"].values])
    return pd.Series(long.reindex(idx).to_numpy(), index=df.index, name="profile")


def deseason(df: pd.DataFrame, mode: str) -> pd.Series:
    u = df["logret"] / profile_variant(df, mode)
    return u.where(np.isfinite(u))


def q7b_tables(df: pd.DataFrame, u: np.ndarray):
    """Lines 357 to 443 of e07_nifty.q7b, returning the tables instead of writing them."""
    sess_codes, sess_index = pd.factorize(df["session"])
    nsess = len(sess_index)
    svol = pd.Series(u).groupby(sess_codes).apply(lambda x: E.robust_sd(x.to_numpy())).reindex(range(nsess)).to_numpy()
    sess_dev = np.asarray(sess_index < E.DEV_END)
    minute = df["minute"].to_numpy()

    def nxt(i):
        return svol[i + 1] if i + 1 < nsess else np.nan

    def prv(i):
        return svol[i - 1] if i >= 1 else np.nan

    configs = [("robust e-detector", E.VolatilityMonitor, cal, band) for cal in E.CAL_LENS for band in E.BANDS]
    configs += [("Gaussian LR", E.GaussianLRMonitor, cal, 0.0) for cal in E.CAL_LENS]
    rng = np.random.default_rng(E.SEED)
    alarms, summary = [], []
    for name, cls, cal, band in configs:
        ev, base = E.run_stream(u, cls, cal, band)
        for e in ev:
            s = sess_codes[e.index]
            alarms.append(dict(detector=name, cal_len=cal, band=band, time=df.index[e.index],
                               session=sess_index[s].date(), period="dev" if sess_dev[s] else "held-out",
                               direction=e.direction, minute_of_day=int(minute[e.index]),
                               delay_steps=e.index - e.monitored_since, baseline=e.baseline_scale,
                               session_vol=svol[s], prev_vol=prv(s), next_vol=nxt(s)))
        a = pd.DataFrame([x for x in alarms if x["detector"] == name and x["cal_len"] == cal and x["band"] == band])
        for period, pmask in (("dev", sess_dev), ("held-out", ~sess_dev)):
            n_s = int(pmask.sum())
            ap = a[a.period == period] if len(a) else a
            row = dict(detector=name, cal_len=cal, band=band, period=period, sessions=n_s, alarms=len(ap),
                       alarms_per_session=len(ap) / n_s)
            for d, cmp in (("up", np.greater), ("down", np.less)):
                ad = ap[ap.direction == d] if len(ap) else ap
                row[f"{d}_alarms"] = len(ad)
                ok = ad.dropna(subset=["next_vol"]) if len(ad) else ad
                k = int(cmp(ok["next_vol"], ok["baseline"]).sum()) if len(ok) else 0
                lo, hi = E.wilson(k, len(ok))
                row[f"{d}_next_vs_baseline"] = k / len(ok) if len(ok) else np.nan
                row[f"{d}_next_vs_baseline_lo"], row[f"{d}_next_vs_baseline_hi"] = lo, hi
                ok2 = ok.dropna(subset=["prev_vol"]) if len(ok) else ok
                k2 = int(cmp(ok2["next_vol"], ok2["prev_vol"]).sum()) if len(ok2) else 0
                row[f"{d}_next_vs_prev"] = k2 / len(ok2) if len(ok2) else np.nan
            row["median_delay_steps"] = float(ap["delay_steps"].median()) if len(ap) else np.nan
            row["share_alarms_first_30min"] = float(np.mean(ap["minute_of_day"] < 9 * 60 + 45)) if len(ap) else np.nan
            pool = np.flatnonzero(np.isfinite(base) & pmask[sess_codes])
            draw = rng.choice(pool, size=min(5000, pool.size), replace=False)
            nv = np.array([nxt(sess_codes[t]) for t in draw])
            bv = base[draw]
            okc = np.isfinite(nv)
            for d, cmp in (("up", np.greater), ("down", np.less)):
                k = int(cmp(nv[okc], bv[okc]).sum())
                row[f"control_{d}_next_vs_baseline"] = k / okc.sum()
            summary.append(row)
    return pd.DataFrame(alarms), pd.DataFrame(summary), svol


def budget_numbers(df: pd.DataFrame, u: np.ndarray) -> dict:
    """e07d.session_table with the given u, then tag_events and fisher_rows (primary configuration)."""
    codes, sessions = pd.factorize(df["session"])
    cal, band = E.PRIMARY
    mon = E.VolatilityMonitor(threshold=E.A, c=E.C, cal_len=cal, null_band=band)
    monitoring = np.zeros(u.size, bool)
    for t, x in enumerate(u):
        mon.update(x)
        monitoring[t] = mon.calibrated
    ev = pd.DataFrame([dict(index=e.index, direction=e.direction) for e in mon.events])
    ev["session"] = sessions[codes[ev["index"]]]
    s = pd.DataFrame({"session": sessions})
    s["date"] = s["session"].dt.date
    s["monitoring_at_open"] = pd.Series(monitoring).groupby(codes).first().to_numpy()
    s["robust_vol"] = pd.Series(u).groupby(codes).apply(lambda x: E.robust_sd(x.to_numpy())).to_numpy()
    s["vol_vs_prev20"] = s["robust_vol"] / s["robust_vol"].rolling(20, min_periods=10).median().shift(1)
    for d in ("up", "down"):
        cnt = ev[ev.direction == d].groupby("session").size()
        s[f"n_{d}"] = s["session"].map(cnt).fillna(0).astype(int)
    s["n_alarms"] = s["n_up"] + s["n_down"]
    s, _, _ = D.tag_events(s)
    fis = pd.DataFrame(D.fisher_rows(s))
    b = fis[(fis.event_type == "budget") & (fis.outcome == "any alarm")].iloc[0]
    return dict(budget_sessions=int(b.event_sessions), budget_with_alarm=int(b.event_sessions_with),
                ordinary_sessions=int(b.ordinary_sessions), ordinary_with_alarm=int(b.ordinary_sessions_with),
                rate_ordinary=float(b.rate_ordinary), fisher_p=float(b.fisher_p))


def main():
    df, info, _ = E.load_clean()
    sessions = pd.DatetimeIndex(df["session"].unique())
    print("REAL DATA", info)
    affected = sessions[:LB // 2]
    print("sessions whose profile uses the future (line 150):", [str(s.date()) for s in affected])

    # the committed variant must equal e07_nifty.session_profile bit for bit
    p_c = profile_variant(df, "committed")
    p_e = E.session_profile(df)
    assert np.array_equal(np.isnan(p_c), np.isnan(p_e)) and np.allclose(p_c.dropna(), p_e.dropna(), rtol=0, atol=0)
    print("committed variant reproduces session_profile exactly: yes")

    # how different the look-ahead profile is from the first past-only profile
    p_n = profile_variant(df, "no_lookahead")
    first_past = p_c[df["session"] == sessions[LB // 2]]
    pooled = p_c[df["session"] == sessions[0]]
    ratio = (pooled.to_numpy() / first_past.to_numpy())
    print("pooled profile / first past-only profile (session 11), quantiles 5/50/95:",
          np.nanquantile(ratio, [0.05, 0.5, 0.95]).round(3))

    us = {m: deseason(df, m).to_numpy() for m in VARIANTS}
    for m in VARIANTS:
        print(f"{m}: monitored returns {int(np.isfinite(us[m]).sum())}")
    print("returns that differ committed vs no_lookahead:",
          int((np.isfinite(us["committed"]) != np.isfinite(us["no_lookahead"])).sum()
              + np.sum(np.isfinite(us["committed"]) & np.isfinite(us["no_lookahead"])
                       & (us["committed"] != us["no_lookahead"]))))
    print("returns that differ no_lookahead vs drop_line150:",
          int((np.isfinite(us["no_lookahead"]) != np.isfinite(us["drop_line150"])).sum()))

    out = {}
    for m in VARIANTS:
        al, sm, svol = q7b_tables(df, us[m])
        al.to_csv(HERE / f"w3_alarms_{m}.csv", index=False)
        sm.to_csv(HERE / f"w3_summary_{m}.csv", index=False)
        out[m] = (al, sm)
        print(f"== {m}: {len(al)} alarms over all 8 configurations")

    # reproduction check against the committed results
    al0 = pd.read_csv(ROOT / "results/e07b_alarms.csv")
    sm0 = pd.read_csv(ROOT / "results/e07b_monitor_summary.csv")
    alc = out["committed"][0]
    same_times = (len(al0) == len(alc)) and (pd.to_datetime(al0["time"]).astype(str).to_numpy()
                                             == alc["time"].astype(str).to_numpy()).all()
    print("committed rerun reproduces results/e07b_alarms.csv (all configs, same times):", bool(same_times))
    key = ["detector", "cal_len", "band", "period"]
    chk = sm0.merge(out["committed"][1], on=key, suffixes=("_res", "_w3"))
    cols = ["alarms", "up_next_vs_baseline", "down_next_vs_baseline", "control_up_next_vs_baseline",
            "control_down_next_vs_baseline", "median_delay_steps"]
    print("committed rerun reproduces results/e07b_monitor_summary.csv:",
          all(np.allclose(chk[c + "_res"], chk[c + "_w3"], equal_nan=True) for c in cols))

    # comparison of the primary configuration
    rows = []
    for m in VARIANTS:
        al, sm = out[m]
        a = al[(al.detector == "robust e-detector") & (al.cal_len == 750) & (al.band == 0.25)].copy()
        a["year"] = pd.to_datetime(a["session"]).dt.year
        p = sm[(sm.detector == "robust e-detector") & (sm.cal_len == 750) & (sm.band == 0.25)].set_index("period")
        r = dict(variant=m, alarms_total=len(a), alarms_2017=int((a.year == 2017).sum()),
                 alarms_apr2017=int(((a.year == 2017) & (pd.to_datetime(a["session"]).dt.month == 4)).sum()))
        for per in ("dev", "held-out"):
            q = p.loc[per]
            r[f"{per}_alarms"] = int(q.alarms)
            r[f"{per}_per_session"] = round(float(q.alarms_per_session), 4)
            r[f"{per}_up_hit"] = round(100 * float(q.up_next_vs_baseline), 2)
            r[f"{per}_up_control"] = round(100 * float(q.control_up_next_vs_baseline), 2)
            r[f"{per}_down_hit"] = round(100 * float(q.down_next_vs_baseline), 2)
            r[f"{per}_down_control"] = round(100 * float(q.control_down_next_vs_baseline), 2)
            r[f"{per}_up_next_vs_prev"] = round(100 * float(q.up_next_vs_prev), 2)
            r[f"{per}_down_next_vs_prev"] = round(100 * float(q.down_next_vs_prev), 2)
            r[f"{per}_median_delay"] = float(q.median_delay_steps)
            r[f"{per}_share_first30"] = round(100 * float(q.share_alarms_first_30min), 2)
        r.update(budget_numbers(df, us[m]))
        rows.append(r)
    cmp = pd.DataFrame(rows).set_index("variant").T
    cmp.to_csv(HERE / "w3_primary_comparison.csv")
    with pd.option_context("display.width", 200, "display.max_rows", 100):
        print(cmp.to_string())

    # alarm times: which ones differ?
    ac = out["committed"][0]
    an = out["no_lookahead"][0]
    for cal in E.CAL_LENS:
        for band in E.BANDS:
            tc = set(ac[(ac.detector == "robust e-detector") & (ac.cal_len == cal) & (ac.band == band)]["time"].astype(str))
            tn = set(an[(an.detector == "robust e-detector") & (an.cal_len == cal) & (an.band == band)]["time"].astype(str))
            only_c, only_n = sorted(tc - tn), sorted(tn - tc)
            print(f"robust cal {cal} band {band}: committed {len(tc)}, no_lookahead {len(tn)}, common {len(tc & tn)}, "
                  f"only committed {len(only_c)}, only no_lookahead {len(only_n)}, "
                  f"last differing time {max(only_c + only_n) if (only_c or only_n) else 'none'}")
            if (cal, band) == E.PRIMARY:
                print("   only committed:", only_c)
                print("   only no_lookahead:", only_n)

    # every printed cell of the docs table: does it change?
    sn = out["no_lookahead"][1]
    sc = out["committed"][1]
    mer = sc.merge(sn, on=key, suffixes=("_c", "_n"))
    chg = []
    for _, r in mer.iterrows():
        d = {}
        if round(r.alarms_per_session_c, 2) != round(r.alarms_per_session_n, 2):
            d["alarms_per_session"] = (round(r.alarms_per_session_c, 2), round(r.alarms_per_session_n, 2))
        for c in ("up_next_vs_baseline", "down_next_vs_baseline", "control_up_next_vs_baseline", "control_down_next_vs_baseline",
                  "up_next_vs_baseline_lo", "up_next_vs_baseline_hi", "down_next_vs_baseline_lo", "down_next_vs_baseline_hi"):
            a_, b_ = round(100 * r[c + "_c"]), round(100 * r[c + "_n"])
            if a_ != b_:
                d[c] = (a_, b_)
        if round(r.median_delay_steps_c) != round(r.median_delay_steps_n):
            d["median_delay"] = (r.median_delay_steps_c, r.median_delay_steps_n)
        if d:
            chg.append(dict(detector=r.detector, cal_len=r.cal_len, band=r.band, period=r.period, changes=d))
    print("docs-table cells that change at printed precision (committed -> no_lookahead):", len(chg))
    for c in chg:
        print("  ", c)


if __name__ == "__main__":
    main()
