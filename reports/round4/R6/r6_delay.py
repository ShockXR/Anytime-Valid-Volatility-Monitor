"""R6 step 2: detection delay on B2's synthetic grid [SYNTHETIC], with E19's vectorised VolatilityMonitor.

Grid (subset of reports/stress/b2_power_delay.py): paths from avmonitor.simulate.change_returns, 1000
calibration steps, change at monitored step nu = 400, variance factors 2 and 0.5 (ratio = sqrt(factor)),
laws Gaussian and Student t3 (heavy_tail_df=3, unit variance), horizon H = 1950 steps (five days of 390).
Thresholds: A = 7800 (README, B2) for the plug-in; cs at 7800 (as B2) and at 9750 = 390 / (0.05 - 0.01)
(E19's rule). Seeds: default_rng(9_975_000 + 10 * law + factor), law 0 = gaussian, 1 = t3; factor 0 = x2, 1 = x0.5.
Usage: python r6_delay.py <paths> <out.csv>
"""
import os, sys, time, itertools
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
os.environ.setdefault("E19_OUT", str(HERE / "pools"))
os.environ.setdefault("E19_SCRATCH", str(HERE / "pools"))
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/experiments")
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/stress")
import e19_real_null as E19  # noqa: E402
from common import wilson  # noqa: E402
from avmonitor.simulate import change_returns  # noqa: E402

A, CAL, H, NU = 7800.0, 1000, 1950, 400
LAWS = {"gaussian": dict(), "student_t3": dict(heavy_tail_df=3)}
FACTORS = [2.0, 0.5]
CONFIGS = [("plugin", 0.25, A), ("plugin", 0.1, A), ("plugin", 0.05, A), ("plugin", 0.0, A),
           ("cs", 0.0, A), ("cs", 0.1, A), ("cs", 0.0, 9750.0), ("cs", 0.1, 9750.0)]


def main(N: int, out: str):
    rows = []
    for (li, (law, kw)), (fi, f) in itertools.product(enumerate(LAWS.items()), enumerate(FACTORS)):
        seed = 9_975_000 + 10 * li + fi
        rng = np.random.default_rng(seed)
        p = change_returns(N, CAL + NU + H, nu=CAL + NU, ratio=float(np.sqrt(f)), rng=rng, **kw)
        t0 = time.time()
        res = E19.x_engine_multi(None, E19.ArrayGen(p.r), CONFIGS, CAL, NU + H, horizon_only=True)
        for (mode, band, Ath), R in zip(CONFIGS, res):
            alarm = np.where(R["first_t"] > 0, R["first_t"] - 1, -1)  # 0-based monitored step, as B2
            early = (alarm >= 0) & (alarm < NU)
            det = (alarm >= NU) & (alarm < NU + H)
            delays = alarm[det] - NU + 1
            lo, hi = wilson(int(det.sum()), N)
            rows.append(dict(law=law, var_factor=f, nu=NU, monitor=mode, band=band, A=Ath, paths=N, seed=seed,
                             detect_rate=det.mean(), det_ci_low=lo, det_ci_high=hi, early_false_alarm=early.mean(),
                             median_delay=float(np.median(delays)) if delays.size else np.nan,
                             q25_delay=float(np.quantile(delays, .25)) if delays.size else np.nan,
                             q75_delay=float(np.quantile(delays, .75)) if delays.size else np.nan,
                             mean_delay=float(delays.mean()) if delays.size else np.nan))
        print(f"{law} x{f}: {N} paths, seed {seed}, {time.time() - t0:.0f} s", flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("EXIT", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2])
