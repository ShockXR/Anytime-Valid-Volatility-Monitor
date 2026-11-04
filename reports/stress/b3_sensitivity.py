"""B3. One-at-a-time sensitivity of the streaming monitor (SYNTHETIC). Seeds: 9_300_000 + 10*scenario (same data stream for every config).
Usage: python b3_sensitivity.py <paths> <out.csv>.  Base config: plug-in baseline, A=7800, cal_len=1000, c=3, null_band=0.25, default bet grids."""
import sys, time
import numpy as np, pandas as pd
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/stress")
from common import wilson
import avmonitor.monitor as MON
from avmonitor.monitor import VolatilityMonitor
from avmonitor.simulate import change_returns
from avmonitor.detectors import DEFAULT_UP_GRID, DEFAULT_DOWN_FRACTIONS

A0 = 7800.0; H = 1950; NU = 400
BASE = dict(threshold=A0, cal_len=1000, c=3.0, null_band=0.25)
MAD0 = MON.MAD_NORMAL_SCALE
SWEEPS = []
for v in [100, 200, 390, 500, 1000, 2000, 4000]: SWEEPS.append(("cal_len", v, dict(cal_len=v)))
for v in [1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 1e6]: SWEEPS.append(("c", v, dict(c=v)))
for v in [0.0, 0.1, 0.25, 0.5]: SWEEPS.append(("null_band", v, dict(null_band=v)))
for v in [0.001, 0.01, 0.05, 0.2]: SWEEPS.append(("delta (cs baseline, band 0)", v, dict(baseline="cs", delta=v, null_band=0.0)))
for v in [0.5, 1.0, MAD0, 2.0]: SWEEPS.append(("MAD scaling constant", round(v, 4), dict(_mad=v)))
for nm, g in [("default 6 bets", DEFAULT_UP_GRID), ("single 0.05", (0.05,)), ("single 0.5", (0.5,)), ("single 0.9", (0.9,)), ("wide 9 bets", (0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9))]:
    SWEEPS.append(("up bet grid", nm, dict(up_grid=g)))
for nm, g in [("default 5 fractions", DEFAULT_DOWN_FRACTIONS), ("single 0.1", (0.1,)), ("single 0.5", (0.5,)), ("single 0.95", (0.95,))]:
    SWEEPS.append(("down bet fractions", nm, dict(down_fractions=g)))
for v in [1950.0, 3900.0, 7800.0, 15600.0, 31200.0, 78000.0]: SWEEPS.append(("threshold A", v, dict(threshold=v)))

def make(cfg):
    cfg = dict(BASE, **cfg); mad = cfg.pop("_mad", MAD0); MON.MAD_NORMAL_SCALE = mad
    return VolatilityMonitor(**cfg)

def scen_data(sc, n, N):
    rng = np.random.default_rng(9_300_000 + 10 * sc)
    if sc == 0: return rng.standard_normal((N, n)), None
    if sc == 1: return rng.standard_t(3, (N, n)) / np.sqrt(3.0), None
    f = 2.0 if sc == 2 else 0.5
    return None, f

def run(args):
    kind, val, cfg, N = args
    full = dict(BASE, **cfg); cal = full["cal_len"]; rows = []
    for sc, name in enumerate(["null gaussian", "null student_t3", "variance x2 at nu=400", "variance x0.5 at nu=400"]):
        n = cal + (H if sc < 2 else NU + H)
        if sc < 2:
            data, _ = scen_data(sc, n, N)
        else:
            rng = np.random.default_rng(9_300_000 + 10 * sc)
            data = change_returns(N, n, nu=cal + NU, ratio=float(np.sqrt(2.0 if sc == 2 else 0.5)), rng=rng).r
        al = np.full(N, -1)
        for i in range(N):
            m = make(cfg)
            for r in data[i]:
                ev = m.update(r)
                if ev is not None:
                    al[i] = ev.index - ev.monitored_since; break
        MON.MAD_NORMAL_SCALE = MAD0
        if sc < 2:
            for hz in (390, 1950):
                k = int(((al >= 0) & (al < hz)).sum()); lo, hi = wilson(k, N)
                rows.append(dict(sweep=kind, value=val, scenario=name, metric=f"false alarm within {hz}", estimate=k / N, ci_low=lo, ci_high=hi, bound=(full["threshold"] and (0.01 if full.get("baseline") == "cs" else 0.0) + hz / full["threshold"])))
        else:
            det = (al >= NU) & (al < NU + H); early = ((al >= 0) & (al < NU)).mean(); d = al[det] - NU + 1
            lo, hi = wilson(int(det.sum()), N)
            rows.append(dict(sweep=kind, value=val, scenario=name, metric="detection rate within 1950", estimate=det.mean(), ci_low=lo, ci_high=hi, bound=np.nan))
            rows.append(dict(sweep=kind, value=val, scenario=name, metric="median delay", estimate=float(np.median(d)) if d.size else np.nan, ci_low=float(np.quantile(d, .25)) if d.size else np.nan, ci_high=float(np.quantile(d, .75)) if d.size else np.nan, bound=np.nan))
            rows.append(dict(sweep=kind, value=val, scenario=name, metric="early false alarm before nu", estimate=early, ci_low=np.nan, ci_high=np.nan, bound=np.nan))
    print(kind, val, "done", flush=True)
    return rows

if __name__ == "__main__":
    N = int(sys.argv[1]); out = sys.argv[2]
    jobs = [(k, v, c, N) for k, v, c in SWEEPS]
    from multiprocessing import Pool
    with Pool(2) as p:
        res = p.map(run, jobs, chunksize=1)
    pd.DataFrame([r for rr in res for r in rr]).to_csv(out, index=False)
    print("EXIT", flush=True)
