"""R4: run the plug-in (README setting) and cs monitors continuously on one asset; save compact run files.

Usage: python run_monitors.py ASSET [ASSET ...]   (BTC ETH XRP SOL DOGE BITSTAMP NIFTY)
Data preparation is reused from experiments/e12_crypto_events.py (crypto) and experiments/e07_nifty.py (Nifty).
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments"))
from avmonitor.data import intraday_scale, load_binance, load_bitstamp_csv  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402

OUT = Path(__file__).resolve().parent / "runs"
OUT.mkdir(exist_ok=True)
COINS = {"BTC": "BTCUSDT", "ETH": "ETHUSDT", "XRP": "XRPUSDT", "SOL": "SOLUSDT", "DOGE": "DOGEUSDT"}
END_UTC = "2026-10-25"  # stop before the pre-registered forward-test window (see PLAN.md)


def prep_crypto(asset):
    if asset == "BITSTAMP":
        df = load_bitstamp_csv(ROOT / "data/raw/crypto/btcusd_bitstamp_1min_latest.csv", None, END_UTC)  # filled="missing"
        mon_start = df.index[0].floor("D") + pd.Timedelta(days=21)
    else:
        df = load_binance(COINS[asset], "2024-01-01", END_UTC, download=False)
        mon_start = pd.Timestamp("2024-01-21", tz="UTC")
    r = df["logret"]
    s = intraday_scale(r, lookback_days=20)  # start="past" (default): past days only; NaN returns kept
    u = (r / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    keep = np.asarray(df.index >= mon_start)
    if not np.isfinite(s.to_numpy()[keep]).all():
        raise ValueError("profile missing for some monitored minute")
    idx = df.index[keep]
    info = dict(first=str(idx[0]), last=str(idx[-1]), steps=int(keep.sum()), nan_returns=int(np.isnan(u[keep]).sum()),
                missing_minutes=int(df["close"].isna().to_numpy()[keep].sum()) if "close" in df else -1)
    year = idx.year.to_numpy()
    tod = (idx.hour * 60 + idx.minute).to_numpy()
    return u[keep], idx.asi8, year, tod, info


def prep_nifty():
    import e07_nifty as E7
    df, info0, _ = E7.load_clean("2017-01-01", "2026-11-01")
    df = df[df.index < pd.Timestamp(END_UTC, tz="UTC")]
    u = E7.deseasonalise(df).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    idx = df.index
    year = np.asarray(pd.DatetimeIndex(df["session"]).year)
    tod = df["minute"].to_numpy()
    info = dict(first=str(idx[0]), last=str(idx[-1]), steps=len(u), nan_returns=int(np.isnan(u).sum()),
                sessions=int(df["session"].nunique()), sessions_dropped_frozen=len(info0["frozen_sessions_dropped"]),
                sessions_dropped_short=len(info0["short_sessions_dropped"]))
    return u, idx.asi8, year, tod, info


def run(u, monitor):
    calib = np.zeros(u.size, bool)  # True while calibrating (blind)
    upd = monitor.update
    for t in range(u.size):
        calib[t] = not monitor.calibrated
        upd(u[t])
    ev = monitor.events
    seg = monitor.segments
    return dict(calib=calib,
                ev_idx=np.array([e.index for e in ev], np.int64),
                ev_up=np.array([e.direction == "up" for e in ev], bool),
                ev_since=np.array([e.monitored_since for e in ev], np.int64),
                seg_start=np.array([s["start"] for s in seg], np.int64),
                seg_cal_len=np.array([s["cal_len"] for s in seg], np.int64),
                seg_tau2=np.array([s["tau2"] for s in seg], float),
                seg_m_up=np.array([s["m_up"] for s in seg], float),
                seg_m0=np.array([s["m0"] for s in seg], float),
                seg_down_active=np.array([s["down_active"] for s in seg], bool))


def main(assets):
    for a in assets:
        t0 = time.time()
        if a == "NIFTY":
            u, ts, year, tod, info = prep_nifty()
            A, Acs, L = 7500.0, 375 / 0.04, 750
        else:
            u, ts, year, tod, info = prep_crypto(a)
            A, Acs, L = 28800.0, 1440 / 0.04, 1440
        print(a, "prepared", info, f"{time.time() - t0:.0f}s", flush=True)
        for name, mon in (("plugin", VolatilityMonitor(threshold=A, c=3.0, cal_len=L, null_band=0.25)),
                          ("cs", VolatilityMonitor(threshold=Acs, c=3.0, cal_len=L, null_band=0.25, baseline="cs", delta=0.01))):
            t1 = time.time()
            res = run(u, mon)
            np.savez_compressed(OUT / f"{a}_{name}.npz", u=u, ts=ts, year=year, tod=tod, **res)
            print(a, name, "alarms", res["ev_idx"].size, "up", int(res["ev_up"].sum()), "blind",
                  round(float(res["calib"].mean()), 3), f"{time.time() - t1:.0f}s", flush=True)
        pd.Series(info).to_json(OUT / f"{a}_info.json")


if __name__ == "__main__":
    main(sys.argv[1:])
