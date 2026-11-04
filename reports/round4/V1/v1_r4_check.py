"""V1: recompute R4 blind share, up/down alarm counts and up hit share for BTC and Nifty asset-years,
running VolatilityMonitor directly (README setting used by R4: plug-in, c=3, band 0.25;
crypto A=28800, cal_len=1440; Nifty A=7500, cal_len=750). Blind time is derived from the monitor's own
segments/events (not from the 'calibrated' flag), as a second route. Usage: python v1_r4_check.py"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path("/home/user/anytime-valid-volatility-monitor")
sys.path.insert(0, str(ROOT / "experiments"))
from avmonitor.data import intraday_scale, load_binance
from avmonitor.monitor import VolatilityMonitor

SCR = Path("/tmp/claude-0/-home-user/413dde72-c15d-5c36-9682-3b5d711ff95d/scratchpad/R4_runs")
T = pd.read_csv(ROOT / "reports/round4/R4/asset_year.csv")


def btc():
    df = load_binance("BTCUSDT", "2024-01-01", "2026-10-25", download=False)
    r = df["logret"]
    s = intraday_scale(r, lookback_days=20)
    u = (r / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    keep = np.asarray(df.index >= pd.Timestamp("2024-01-21", tz="UTC"))
    return u[keep], df.index[keep].year.to_numpy(), 28800.0, 1440


def nifty():
    import e07_nifty as E7
    df, _, _ = E7.load_clean("2017-01-01", "2026-11-01")
    df = df[df.index < pd.Timestamp("2026-10-25", tz="UTC")]
    u = E7.deseasonalise(df).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    return u, np.asarray(pd.DatetimeIndex(df["session"]).year), 7500.0, 750


for name, prep, asset, years in (("BTC", btc, "BTC", [2024, 2025, 2026]), ("NIFTY", nifty, "NIFTY", [2023, 2024, 2025])):
    t0 = time.time()
    u, year, A, L = prep()
    mon = VolatilityMonitor(threshold=A, c=3.0, cal_len=L, null_band=0.25)
    for x in u:
        mon.update(x)
    n = u.size
    # monitored = from each segment's start to its alarm (inclusive) or to the end of data
    watched = np.zeros(n, bool)
    ev = mon.events
    ends = {e.monitored_since: e.index for e in ev}
    for sgm in mon.segments:
        st = sgm["start"]
        watched[st:ends.get(st, n - 1) + 1] = True
    blind = ~watched
    ev_idx = np.array([e.index for e in ev]); ev_up = np.array([e.direction == "up" for e in ev])
    # hit share of up alarms: mean u^2 over next 60 minutes (>= 50 finite) vs calibration-window mean u^2
    fin = np.isfinite(u); u2 = np.where(fin, u * u, 0.0)
    seg = {s_["start"]: s_ for s_ in mon.segments}
    # compare with R4 saved arrays
    zf = SCR / f"{asset}_plugin.npz"
    if zf.exists():
        z = np.load(zf)
        print(name, "u equal to R4 saved run:", bool(np.array_equal(np.nan_to_num(z["u"], nan=-9), np.nan_to_num(u, nan=-9))),
              "| calib equal to my blind:", bool(np.array_equal(z["calib"], blind)),
              "| alarm idx equal:", bool(np.array_equal(z["ev_idx"], ev_idx)), flush=True)
    for y in years:
        my = year == y
        ey = year[ev_idx] == y
        hits = []
        for e, isup in zip(ev_idx[ey], ev_up[ey]):
            if not isup or e + 60 >= n or e - 59 < 0:
                continue
            npost = fin[e + 1:e + 61].sum(); npre = fin[e - 59:e + 1].sum()
            if npost < 50 or npre < 50:
                continue
            post = u2[e + 1:e + 61].sum() / npost
            pre = u2[e - 59:e + 1].sum() / npre
            if not (post > 0 and pre > 0):
                continue
            since = [ev_k.monitored_since for ev_k in ev if ev_k.index == e][0]
            cl = seg[since]["cal_len"]
            a = max(since - cl, 0)
            base = u2[a:since].sum() / max(fin[a:since].sum(), 1)
            hits.append(post > base)
        r = T[(T.asset == asset) & (T.monitor == "plugin") & (T.year == y)].iloc[0]
        print(f"{name} {y}: steps {my.sum()} (R4 {r.steps}) | blind {blind[my].mean():.4f} (R4 {r.blind_share:.4f}) | "
              f"up {int((ey & ev_up).sum())} (R4 {r.n_up}) | down {int((ey & ~ev_up).sum())} (R4 {r.n_down}) | "
              f"up hit {np.mean(hits):.3f} over {len(hits)} (R4 {r.up_hit:.3f} over {r.up_n_valid})", flush=True)
    print(name, "time", round(time.time() - t0), "s", flush=True)
