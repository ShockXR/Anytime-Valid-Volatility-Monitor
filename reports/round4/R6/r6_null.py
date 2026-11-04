"""R6 step 1: false alarms of six settings on E19's permuted calm-day streams [REAL marginal, constructed streams].

Reuses E19's pool (built by r6_pools.py with E19.prepare_pools), E19's StreamGen("perm") and E19's vectorised
VolatilityMonitor (x_engine_multi). Only the seeds differ from E19: numpy default_rng(9_970_000 + 100 * d + b)
for data set index d and batch b (all inside 9,970,000 to 9,979,999).
Settings (mode, band, A): plug-in at A = day / 0.05 with bands 0.25 (README), 0.1, 0.05, 0; cs (delta 0.01)
at A_cs = day / (0.05 - 0.01) with bands 0 and 0.1. Each stream runs until its first alarm or 5 monitored days.
Usage: python r6_null.py <dataset> <batch> <rows>
"""
import os, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
os.environ.setdefault("E19_OUT", str(HERE / "pools"))
os.environ.setdefault("E19_SCRATCH", str(HERE / "pools"))
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/experiments")
import e19_real_null as E19  # noqa: E402

SEED0 = 9_970_000
CONFIGS = [("plugin", 0.25), ("plugin", 0.1), ("plugin", 0.05), ("plugin", 0.0), ("cs", 0.0), ("cs", 0.1)]


def run(name: str, b: int, rows: int):
    d = list(E19.DATASETS).index(name)
    day, cal, A, A_cs, _ = E19.DATASETS[name]
    m = 5 * day
    z = np.load(E19.pool_path(name))
    pool, tau2, m0 = z["pool"], float(z["tau2"]), float(z["m0"])
    length = cal + m + E19.EXTRA_CALS * cal
    seed = SEED0 + 100 * d + b
    gen = E19.StreamGen("perm", pool.size, rows, length, np.random.default_rng(seed))
    configs = [(mode, band, A_cs if mode == "cs" else A) for mode, band in CONFIGS]
    t0 = time.time()
    res = E19.x_engine_multi(pool, gen, configs, cal, m, known=(tau2, m0), horizon_only=True)
    out = []
    for (mode, band, Ath), R in zip(configs, res):
        out.append(dict(dataset=name, monitor=mode, band=band, A=Ath, batch=b, seed=seed, rows=rows,
                        first_t=R["first_t"], first_dir=R["first_dir"], mon_steps=R["mon_steps"], ended=R["ended"]))
    ck = HERE / "null_batches" / f"{name.replace(' ', '_')}_{b}.pkl"
    ck.parent.mkdir(exist_ok=True)
    pd.to_pickle(out, ck)
    print(f"{name} batch {b}: {rows} streams, seed {seed}, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
