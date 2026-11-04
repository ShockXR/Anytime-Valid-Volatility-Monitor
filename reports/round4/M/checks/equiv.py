"""Run the changed code paths on real and simulated inputs and save every output, for an old-versus-new comparison.

usage: python -I equiv.py <src dir> <repo root for experiments> <data root> <out.json>
"""
import sys
src, repo, data_root, out_path = sys.argv[1:5]
sys.path.insert(0, src)
sys.path.insert(0, repo + "/experiments")
import json, hashlib, importlib
import numpy as np, pandas as pd
import avmonitor
from pathlib import Path
from avmonitor import data as D, detectors as Det
from avmonitor.monitor import VolatilityMonitor, GaussianLRMonitor
from avmonitor.periodic import PeriodicVolatilityMonitor
from avmonitor.guard import GuardedVolatilityMonitor, GuardedPeriodicVolatilityMonitor

res = {"avmonitor_file": avmonitor.__file__}

def ev(m):
    return [[int(e.index), e.direction, int(e.monitored_since), repr(float(e.baseline_scale)), repr(float(e.stat))] for e in m.events]

def seg_hash(m):
    return hashlib.sha256(repr([sorted((k, repr(v)) for k, v in s.items()) for s in m.segments]).encode()).hexdigest()

# 1. monitors on real data: BTCUSDT and Nifty, deseasonalised as in the CLI
root = Path(data_root)
btc = D.load_binance("BTCUSDT", "2025-06-01", "2026-01-01", root=root, download=False)
s = D.intraday_scale(btc["logret"], lookback_days=20)
u_btc = (btc["logret"] / s).to_numpy(copy=True)
u_btc[~np.isfinite(u_btc)] = np.nan
nif = D.regular_sessions(D.load_nifty("2024-01-01", "2025-01-01", root=root, download=False))
s = D.intraday_scale(nif["logret"], lookback_days=20, circular=False)
u_nif = (nif["logret"] / s).to_numpy(copy=True)
u_nif[~np.isfinite(u_nif)] = np.nan
bs = D.load_bitstamp_csv(root / "crypto" / "btcusd_bitstamp_1min_latest.csv", "2025-07-01", "2026-01-01", filled="zero")
s = D.intraday_scale(bs["logret"], lookback_days=20)
u_bs = (bs["logret"] / s).to_numpy(copy=True)
u_bs[~np.isfinite(u_bs)] = np.nan
configs = {
    "plugin_band25": lambda spd: VolatilityMonitor(threshold=spd / 0.05, c=3.0, cal_len=spd, null_band=0.25),
    "plugin_band0": lambda spd: VolatilityMonitor(threshold=spd / 0.05, c=3.0, cal_len=spd),
    "gaussian_lr": lambda spd: GaussianLRMonitor(threshold=spd / 0.05, c=3.0, cal_len=spd),
    "cs": lambda spd: VolatilityMonitor(threshold=36000.0, c=3.0, cal_len=spd, null_band=0.25, baseline="cs", delta=0.01),
    "cs_cdown": lambda spd: VolatilityMonitor(threshold=36000.0, c=3.0, cal_len=spd, baseline="cs", delta=0.01, c_down=1.5),
    "periodic_cs": lambda spd: PeriodicVolatilityMonitor(threshold=36000.0, c=3.0, cal_len=spd, null_band=0.25, baseline="cs", delta=0.01, period=spd),
    "guarded_rank": lambda spd: GuardedPeriodicVolatilityMonitor(threshold=36000.0, c=3.0, cal_len=spd, null_band=0.25, baseline="cs", delta=0.01, period=spd, guard_tries=3, guard_eta=0.25, guard_fallback="rank", rank_seed=7),
}
for name, u, spd in (("btc", u_btc, 1440), ("nifty", u_nif, 375), ("bitstamp_zero", u_bs, 1440)):
    for cname, mk in configs.items():
        m = mk(spd)
        m.run(u)
        res[f"{name}/{cname}"] = dict(events=ev(m), n_events=len(m.events), segments=seg_hash(m), cal_stats=m.cal_stats)

# 2. E4's benchmark detectors on E4's own simulated paths (first 30 of 300 rows, full 20,000-step horizon)
E4 = importlib.import_module("e04_detection")
from avmonitor.simulate import change_returns
for i, (nm, kw) in enumerate(E4.SCENARIOS.items()):
    rng = np.random.default_rng(100 + i)
    p = change_returns(300, E4.CAL + 20000, rng=rng, **kw)
    r, cal = p.r[:30, E4.CAL:], p.r[:30, :E4.CAL]
    bases = [Det.calibrate(cal[j], c=3.0) for j in range(30)]
    scale = np.array([b.scale for b in bases])[:, None]
    out = {
        "glr_up": Det.gaussian_lr_edetector(r / scale, 1.0, threshold=E4.A).alarm,
        "glr_down": Det.gaussian_lr_edetector(r / scale, 1.0, ratios=(0.5, 0.7, 0.8), threshold=E4.A).alarm,
        "page": Det.page_cusum(r**2 / scale**2, 1.0, 0.25, 12.0).alarm,
        "z60": Det.repeated_window_ztest(r / scale, 1.0, 2.0, 60).alarm,
    }
    res[f"e04/{nm}"] = {k: v.tolist() for k, v in out.items()}
json.dump(res, open(out_path, "w"), indent=0)
print("done", out_path)
