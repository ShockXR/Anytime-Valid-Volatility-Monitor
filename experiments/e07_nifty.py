"""E7. Nifty 50 index, 1-minute bars, April 2017 to October 2026 (an Indian market application).

REAL MARKET DATA. Source: github.com/technovusin/nifty50-historical-data, loaded with
``avmonitor.data.load_nifty`` (IST, 09:15 to 15:29, returns within sessions only, so
overnight gaps never enter). Cross-checked against NSE daily closes in e07c.

Research questions
  Q7a. Describe Nifty minute returns the way Table 5 of Bapat (2026) describes crypto
       (quantiles, extremes, zero-return share) and compute daily integrated
       quarticity (IQ) on a common scale with the estimators of ``avmonitor.robust``
       (as in E5, Q5b). Is the MAD truncation bias on Nifty smaller or larger than on
       crypto, and why (intraday seasonality, fat tails)?
  Q7b. Run the streaming monitor on deseasonalised Nifty returns, with settings fixed
       from first principles: one session = 375 steps, A = 375 / 0.05 = 7500 (at most a
       5% chance of any false alarm per session under the null), calibration length
       375 or 750, null band 0, 0.25 or 0.5. Development period 2017 to 2021, held-out
       period 2022 to Oct 2026. Report alarms per session, the share of up alarms
       followed by a higher-volatility next session and of down alarms followed by a
       lower one (Wilson 95% intervals), against random times and against the Gaussian
       likelihood ratio benchmark.

The package function ``intraday_scale`` was written for 24-hour crypto and smooths the
profile circularly across the clock. On a 375-bar session that wraps the closing minutes
into the opening ones. This script therefore uses a session-aware profile
(``session_profile``) and reports the check that led to it (results/e07b_deseason_check.csv).

Usage
  python experiments/e07_nifty.py q7a
  python experiments/e07_nifty.py q7b-dev     # development period only (2017 to 2021)
  python experiments/e07_nifty.py q7b         # full run, dev and held-out reported separately
  python experiments/e07_nifty.py q7a --bitstamp-filled missing   # BTC rows of Q7a with the loader's default

The Bitstamp BTC/USD rows of Q7a (review round 9, m5). The committed rows were computed
(i) with the Bitstamp loader of before review round 9, which turned the provider's filled
(zero-volume) minutes into zero returns, and (ii) with the pooled start of the intraday
profile (``intraday_scale(..., start="pooled")``), which was the default until the
profile was made strictly predictable later the same day (commit 0e7d160). Q7a asks for
both explicitly (BITSTAMP_FILLED = "zero"), so that the script reproduces the committed
results; ``--bitstamp-filled missing`` makes the returns on both sides of a filled minute
missing. The effect on every committed cell is in results/r9_bitstamp_impact.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from avmonitor import robust as R
from avmonitor.data import intraday_scale, load_bitstamp_csv, load_nifty, regular_sessions
from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
TZ = "Asia/Kolkata"

# ----------------------------------------------------------------------------
# Settings fixed before any monitor run (first principles, see docstring)
# ----------------------------------------------------------------------------
STEPS_PER_SESSION = 375
A = STEPS_PER_SESSION / 0.05  # 7500
C = 3.0
CAL_LENS = (375, 750)
BANDS = (0.0, 0.25, 0.5)
PRIMARY = (750, 0.25)  # pre-registered primary configuration
DEV_END = pd.Timestamp("2022-01-01", tz=TZ)
# "zero": filled Bitstamp minutes give zero returns, as in the committed Q7a results; "missing": the loader's
# default (review round 9, m5). Set from the command line with --bitstamp-filled.
BITSTAMP_FILLED = "zero"
FROZEN_RUN = 5  # a 50-stock index that does not move for 5+ minutes is halted or frozen
MAX_FLAGGED_SHARE = 0.25  # drop sessions where more than this share of returns is flagged
LOOKBACK = 20  # sessions in the intraday profile
SMOOTH = 5  # minutes in the centred moving average of the profile (non-circular)
SEED = 7


# ----------------------------------------------------------------------------
# Data preparation (shared with e07c and e07d)
# ----------------------------------------------------------------------------


def load_clean(start: str = "2017-01-01", end: str = "2026-11-01"):
    """Regular Nifty sessions with frozen-price stretches set to NaN.

    A run of ``FROZEN_RUN`` or more exact zero returns in a row means the index did not
    update (trading halt, exchange outage, or the frozen last minutes seen from August
    2026). Those returns and the first return after the run (which spans the whole
    frozen stretch) are set to NaN. Sessions with more than ``MAX_FLAGGED_SHARE`` of
    their returns flagged are dropped. Returns (frame, info dict, flagged-runs table).
    """
    raw = load_nifty(start, end, download=False)
    df = regular_sessions(raw).copy()
    short = sorted(set(raw["session"].unique()) - set(df["session"].unique()))
    r = df["logret"]
    z = (r == 0).to_numpy()
    sess = df["session"].to_numpy()
    flag = np.zeros(len(df), bool)
    runs = []
    i, n = 0, len(df)
    while i < n:
        if z[i]:
            j = i
            while j + 1 < n and z[j + 1] and sess[j + 1] == sess[i]:
                j += 1
            L = j - i + 1
            if L >= FROZEN_RUN:
                flag[i : j + 1] = True
                if j + 1 < n and sess[j + 1] == sess[i]:
                    flag[j + 1] = True
                runs.append(dict(session=df.index[i].date(), start=df.index[i].strftime("%H:%M"), minutes=L))
            i = j + 1
        else:
            i += 1
    df["flagged"] = flag
    df.loc[flag, "logret"] = np.nan
    share = df.groupby("session")["flagged"].mean()
    bad = share[share > MAX_FLAGGED_SHARE].index
    df = df[~df["session"].isin(bad)].copy()
    df["minute"] = df.index.hour * 60 + df.index.minute
    info = dict(sessions_raw=int(raw["session"].nunique()), short_sessions_dropped=[str(s.date()) for s in short],
                frozen_sessions_dropped=[str(s.date()) for s in bad], sessions=int(df["session"].nunique()),
                returns=int(df["logret"].notna().sum()), flagged_returns=int(df["flagged"].sum()))
    return df, info, pd.DataFrame(runs)


def session_profile(df: pd.DataFrame, lookback: int = LOOKBACK, smooth: int = SMOOTH) -> pd.Series:
    """Predictable intraday volatility profile for sessions with a fixed clock (Nifty).

    For each minute of the session the profile is the median of ``|r|`` at that minute
    over the previous ``lookback`` sessions (current session excluded), smoothed with a
    centred moving average of ``smooth`` minutes that is truncated at the open and the
    close (no wrap-around) and skips missing values, then normalised to median 1 across
    the session. The first sessions get the profile of the first ``lookback`` sessions.
    """
    a = df["logret"].abs()
    wide = pd.DataFrame({"a": a.values, "s": df["session"].values, "m": df["minute"].values}).pivot_table(
        index="s", columns="m", values="a", aggfunc="first", dropna=False)
    prof = wide.rolling(lookback, min_periods=max(5, lookback // 2)).median().shift(1)
    prof = prof.fillna(wide.iloc[:lookback].median())
    if smooth > 1:
        prof = prof.T.rolling(smooth, center=True, min_periods=1).mean().T
    prof = prof.div(prof.median(axis=1), axis=0)
    long = prof.stack(future_stack=True)
    idx = pd.MultiIndex.from_arrays([df["session"].values, df["minute"].values])
    return pd.Series(long.reindex(idx).to_numpy(), index=df.index, name="profile")


def deseasonalise(df: pd.DataFrame, **kw) -> pd.Series:
    s = session_profile(df, **kw)
    u = df["logret"] / s
    return u.where(np.isfinite(u))


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (c - h, c + h)


def robust_sd(x) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.median(np.abs(x)) * R.MAD_NORMAL_SCALE) if x.size >= 50 else np.nan


# ----------------------------------------------------------------------------
# Q7a: descriptive statistics and daily IQ
# ----------------------------------------------------------------------------


def describe(r: np.ndarray, name: str) -> dict:
    r = r[np.isfinite(r)]
    q = np.quantile(r, [0, 0.0001, 0.25, 0.5, 0.75, 0.9999, 1])
    s = np.median(np.abs(r)) * R.MAD_NORMAL_SCALE
    return dict(series=name, n=r.size, minimum=q[0], q0001=q[1], Q1=q[2], median=q[3], mean=r.mean(), Q3=q[4],
                q9999=q[5], maximum=q[6], zero_share=float(np.mean(r == 0)), sd=r.std(),
                robust_sd=s, max_abs_over_robust_sd=np.abs(r).max() / s,
                excess_kurtosis=float(stats.kurtosis(r)), share_beyond_3_robust_sd=float(np.mean(np.abs(r) > 3 * s)),
                lag1_autocorr=float(np.corrcoef(r[1:], r[:-1])[0, 1]))


def iq_row(r: np.ndarray) -> dict:
    """Daily IQ estimators on a common scale (T = 1 session), as in E5 Q5b."""
    r = r[np.isfinite(r)]
    s = np.median(np.abs(r)) * R.MAD_NORMAL_SCALE
    return dict(n=r.size, RQ=R.realized_quarticity(r), TPQ=R.tripower_quarticity(r), QPQ=R.quadpower_quarticity(r),
                MinRQ=R.min_rq(r), MedRQ=R.med_rq(r), TRQ=R.threshold_rq(r), Kernel=R.kernel_iq(r, 30),
                RobustKernel_scaled_bc=R.kernel_iq(r, 30, c=3, mad_scaled=True, bias_correct=True),
                RobustKernel_raw=R.kernel_iq(r, 30, c=3, mad_scaled=False),
                TruncRQ_scaled_bc=R.truncated_rq_mad(r, c=3, mad_scaled=True, bias_correct=True),
                TruncRQ_raw=R.truncated_rq_mad(r, c=3, mad_scaled=False),
                RV=np.sum(r**2), clipped_share_scaled=float(np.mean(np.abs(r) > 3 * s)),
                excess_kurtosis=float(stats.kurtosis(r)), zero_share=float(np.mean(r == 0)))


IQ_COLS = ["RQ", "TPQ", "QPQ", "MinRQ", "MedRQ", "TRQ", "Kernel", "RobustKernel_scaled_bc", "RobustKernel_raw",
           "TruncRQ_scaled_bc", "TruncRQ_raw"]


def iq_table(groups) -> pd.DataFrame:
    rows = []
    for day, r in groups:
        if np.isfinite(r).sum() < 300:
            continue
        rows.append(dict(day=day, **iq_row(r)))
    return pd.DataFrame(rows).set_index("day")


def iq_summary(d: pd.DataFrame, label: str) -> pd.DataFrame:
    ratio = d[IQ_COLS].div(d["MedRQ"], axis=0)
    out = pd.DataFrame({"median_daily_IQ": d[IQ_COLS].median(), "median_ratio_to_MedRQ": ratio.median()})
    out["series"] = label
    out["days"] = len(d)
    out["median_clipped_share"] = d["clipped_share_scaled"].median()
    out["median_excess_kurtosis"] = d["excess_kurtosis"].median()
    return out.reset_index().rename(columns={"index": "estimator"})


def q7a():
    df, info, runs = load_clean()
    print("REAL DATA", info)
    runs.to_csv(OUT / "e07a_frozen_runs.csv", index=False)
    r = df["logret"].to_numpy()
    dev = (df.index < DEV_END)
    rows = [describe(r, "Nifty 50, Apr 2017 to Oct 2026 [REAL]"), describe(r[dev], "Nifty 50, 2017 to 2021 [REAL]"),
            describe(r[~dev], "Nifty 50, 2022 to Oct 2026 [REAL]")]
    # Nifty deseasonalised (scale-free: in units of the session's median-minute volatility)
    u = deseasonalise(df).to_numpy()
    rows.append(describe(u, "Nifty 50 deseasonalised, full sample [REAL]"))
    # Bitstamp BTC/USD, same period as E5 and Table 5 of the paper, for a like-for-like comparison
    btc = load_bitstamp_csv(ROOT / "data/raw/crypto/btcusd_bitstamp_1min_latest.csv", "2025-07-01", "2026-01-01",
                            filled=BITSTAMP_FILLED)
    rb = btc["logret"].to_numpy()
    rows.append(describe(rb, "BTC/USD Bitstamp, Jul to Dec 2025 [REAL]"))
    # circular profile is right for 24h crypto; pooled start as in the committed results (see the docstring)
    sb = intraday_scale(btc["logret"].fillna(0.0), lookback_days=20, start="pooled")
    ub = (btc["logret"] / sb).to_numpy(copy=True)
    ub[~np.isfinite(ub)] = np.nan
    rows.append(describe(ub, "BTC/USD deseasonalised [REAL]"))
    desc = pd.DataFrame(rows)
    desc.to_csv(OUT / "e07a_table5_style.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(desc.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    lr = df["logret"]
    print("Nifty extremes: min", lr.idxmin(), lr.min(), "| max", lr.idxmax(), lr.max())
    print("five largest |r|:", lr.abs().nlargest(5).to_dict())
    yr = pd.DataFrame(dict(r=r, y=df.index.year)).groupby("y")["r"].agg(
        zero_share=lambda x: float(np.mean(x[np.isfinite(x)] == 0)),
        robust_sd=lambda x: robust_sd(x), n="count")
    yr.to_csv(OUT / "e07a_by_year.csv")
    print(yr)

    # Daily IQ on raw and deseasonalised returns, Nifty and BTC
    sess = df["session"].to_numpy()
    groups_raw = [(k, g["logret"].to_numpy()) for k, g in df.groupby("session")]
    ud = pd.Series(u, index=df.index)
    groups_u = [(k, ud[g.index].to_numpy()) for k, g in df.groupby("session")]
    d_raw = iq_table(groups_raw)
    d_u = iq_table(groups_u)
    d_raw.to_csv(OUT / "e07a_daily_iq_nifty.csv", float_format="%.6g")
    bday = btc.index.floor("D")
    b_raw = iq_table([(k, g.to_numpy()) for k, g in btc["logret"].groupby(bday)])
    b_u = iq_table([(k, g.to_numpy()) for k, g in pd.Series(ub, index=btc.index).groupby(bday)])
    summ = pd.concat([iq_summary(d_raw, "Nifty raw returns"), iq_summary(d_u, "Nifty deseasonalised"),
                      iq_summary(b_raw, "BTC raw returns"), iq_summary(b_u, "BTC deseasonalised")])
    summ.to_csv(OUT / "e07a_iq_summary.csv", index=False)
    piv = summ.pivot(index="estimator", columns="series", values="median_ratio_to_MedRQ").loc[IQ_COLS]
    print(piv.round(3).to_string())
    print(summ.groupby("series")[["days", "median_clipped_share", "median_excess_kurtosis"]].first())
    # Back-of-envelope for exact zeros: with a share p0 of zero returns mixed into Gaussian returns,
    # median |r| falls, so tau = 3 scaled MADs sits at fewer than 3 sd of the non-zero returns.
    rows = []
    for name, p0 in (("Gaussian, no zeros", 0.0), ("Nifty zero share", float(np.mean(r[np.isfinite(r)] == 0))),
                     ("BTC zero share", float(np.mean(rb[np.isfinite(rb)] == 0)))):
        q = (0.5 - p0) / (1 - p0)
        a = 3.0 * stats.norm.ppf((1 + q) / 2) / stats.norm.ppf(0.75)
        rows.append(dict(case=name, zero_share=p0, cutoff_in_sd=a, clipped_share=2 * stats.norm.sf(a),
                         kappa4=R.truncated_moment_factor(a, 4)))
    zz = pd.DataFrame(rows)
    zz.to_csv(OUT / "e07a_zero_share_effect.csv", index=False)
    print(zz.round(4).to_string(index=False))
    del sess


# ----------------------------------------------------------------------------
# Q7b: streaming monitor
# ----------------------------------------------------------------------------


def deseason_check(df: pd.DataFrame, period_mask: np.ndarray) -> pd.DataFrame:
    """How flat is median |u| across the session after deseasonalising? (1 = perfectly flat)."""
    rows = []
    r = df["logret"]
    variants = {
        "package intraday_scale (circular, smooth 15)": intraday_scale(r, lookback_days=LOOKBACK),
        "session_profile smooth 1": session_profile(df, smooth=1),
        "session_profile smooth 5": session_profile(df, smooth=5),
        "session_profile smooth 15": session_profile(df, smooth=15),
        "no deseasonalisation": pd.Series(1.0, index=df.index),
    }
    m = df["minute"].to_numpy()
    for name, s in variants.items():
        u = (r / s).where(lambda x: np.isfinite(x)).to_numpy()
        sel = period_mask & np.isfinite(u)
        med = pd.Series(np.abs(u[sel])).groupby(m[sel]).median()
        med = med / med.median()
        rows.append(dict(variant=name, first_minute_0916=med.iloc[0], mean_0916_0920=med.iloc[:5].mean(),
                         mean_0921_0945=med.iloc[5:30].mean(), mean_1525_1529=med.iloc[-5:].mean(),
                         sd_log_across_minutes=float(np.log(med).std()),
                         max_over_min=float(med.max() / med.min())))
    return pd.DataFrame(rows)


def run_stream(u: np.ndarray, cls, cal_len: int, band: float):
    mon = cls(threshold=A, c=C, cal_len=cal_len, null_band=band)
    base = np.full(u.size, np.nan)
    for t, x in enumerate(u):
        mon.update(x)
        if mon.calibrated:
            base[t] = mon.scale
    return mon.events, base


def q7b(full: bool):
    df, info, _ = load_clean()
    if not full:
        df = df[df.index < DEV_END].copy()
    u_s = deseasonalise(df)
    u = u_s.to_numpy()
    sess_codes, sess_index = pd.factorize(df["session"])
    nsess = len(sess_index)
    dev_mask = (df.index < DEV_END)
    tag = "" if full else "_devonly"
    chk = deseason_check(df, dev_mask).assign(period="dev")
    if full:
        chk = pd.concat([chk, deseason_check(df, ~dev_mask).assign(period="held-out")])
    chk.to_csv(OUT / f"e07b_deseason_check{tag}.csv", index=False)
    print(chk.round(3).to_string(index=False))

    # per-session robust volatility of deseasonalised returns
    svol = pd.Series(u).groupby(sess_codes).apply(lambda x: robust_sd(x.to_numpy())).reindex(range(nsess)).to_numpy()
    sess_dev = np.asarray(sess_index < DEV_END)
    minute = df["minute"].to_numpy()

    def nxt(i):
        return svol[i + 1] if i + 1 < nsess else np.nan

    def prv(i):
        return svol[i - 1] if i >= 1 else np.nan

    configs = [("robust e-detector", VolatilityMonitor, cal, band) for cal in CAL_LENS for band in BANDS]
    configs += [("Gaussian LR", GaussianLRMonitor, cal, 0.0) for cal in CAL_LENS]
    rng = np.random.default_rng(SEED)
    alarms, summary = [], []
    for name, cls, cal, band in configs:
        ev, base = run_stream(u, cls, cal, band)
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
            if n_s == 0:
                continue
            ap = a[a.period == period] if len(a) else a
            row = dict(detector=name, cal_len=cal, band=band, period=period, sessions=n_s, alarms=len(ap),
                       alarms_per_session=len(ap) / n_s)
            for d, cmp in (("up", np.greater), ("down", np.less)):
                ad = ap[ap.direction == d] if len(ap) else ap
                row[f"{d}_alarms"] = len(ad)
                row[f"sessions_with_{d}"] = ad["session"].nunique() if len(ad) else 0
                ok = ad.dropna(subset=["next_vol"]) if len(ad) else ad
                k = int(cmp(ok["next_vol"], ok["baseline"]).sum()) if len(ok) else 0
                lo, hi = wilson(k, len(ok))
                row[f"{d}_next_vs_baseline"] = k / len(ok) if len(ok) else np.nan
                row[f"{d}_next_vs_baseline_lo"], row[f"{d}_next_vs_baseline_hi"] = lo, hi
                ok2 = ok.dropna(subset=["prev_vol"]) if len(ok) else ok
                k2 = int(cmp(ok2["next_vol"], ok2["prev_vol"]).sum()) if len(ok2) else 0
                lo2, hi2 = wilson(k2, len(ok2))
                row[f"{d}_next_vs_prev"] = k2 / len(ok2) if len(ok2) else np.nan
                row[f"{d}_next_vs_prev_lo"], row[f"{d}_next_vs_prev_hi"] = lo2, hi2
                row[f"{d}_median_next_over_baseline"] = float(np.median(ok["next_vol"] / ok["baseline"])) if len(ok) else np.nan
            row["median_delay_steps"] = float(ap["delay_steps"].median()) if len(ap) else np.nan
            row["share_alarms_first_30min"] = float(np.mean(ap["minute_of_day"] < 9 * 60 + 45)) if len(ap) else np.nan
            # random-time control: same run, random monitoring-mode minutes in this period, current baseline
            pool = np.flatnonzero(np.isfinite(base) & pmask[sess_codes])
            draw = rng.choice(pool, size=min(5000, pool.size), replace=False)
            nv = np.array([nxt(sess_codes[t]) for t in draw])
            bv = base[draw]
            okc = np.isfinite(nv)
            for d, cmp in (("up", np.greater), ("down", np.less)):
                k = int(cmp(nv[okc], bv[okc]).sum())
                row[f"control_{d}_next_vs_baseline"] = k / okc.sum()
                row[f"control_{d}_lo"], row[f"control_{d}_hi"] = wilson(k, int(okc.sum()))
            summary.append(row)
    # session-level control for next vs previous session (no detector involved)
    ctrl = []
    for period, pmask in (("dev", sess_dev), ("held-out", ~sess_dev)):
        idx = np.flatnonzero(pmask)
        idx = idx[(idx >= 1) & (idx + 1 < nsess)]
        nv, pv = svol[idx + 1], svol[idx - 1]
        ok = np.isfinite(nv) & np.isfinite(pv)
        for d, cmp in (("up", np.greater), ("down", np.less)):
            k = int(cmp(nv[ok], pv[ok]).sum())
            lo, hi = wilson(k, int(ok.sum()))
            ctrl.append(dict(period=period, direction=d, sessions=int(ok.sum()), share_next_vs_prev=k / ok.sum(), lo=lo, hi=hi))
    pd.DataFrame(ctrl).to_csv(OUT / f"e07b_control_sessions{tag}.csv", index=False)
    print(pd.DataFrame(ctrl).round(3).to_string(index=False))
    al = pd.DataFrame(alarms)
    al.to_csv(OUT / f"e07b_alarms{tag}.csv", index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(OUT / f"e07b_monitor_summary{tag}.csv", index=False)
    show = ["detector", "cal_len", "band", "period", "alarms", "alarms_per_session", "up_alarms", "down_alarms",
            "up_next_vs_baseline", "control_up_next_vs_baseline", "down_next_vs_baseline",
            "control_down_next_vs_baseline", "up_next_vs_prev", "down_next_vs_prev", "median_delay_steps",
            "share_alarms_first_30min"]
    with pd.option_context("display.width", 300, "display.max_columns", 40):
        print(sm[show].round(3).to_string(index=False))
    pd.DataFrame({"session": sess_index.date, "robust_vol": svol}).to_csv(OUT / f"e07b_session_vol{tag}.csv", index=False)


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--bitstamp-filled" in args:
        i = args.index("--bitstamp-filled")
        BITSTAMP_FILLED = args[i + 1]
        del args[i : i + 2]
    cmd = args[0] if args else "all"
    if cmd in ("q7a", "all"):
        q7a()
    if cmd == "q7b-dev":
        q7b(full=False)
    if cmd in ("q7b", "all"):
        q7b(full=True)
