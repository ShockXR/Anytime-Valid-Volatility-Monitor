"""R6 check: E19's vectorised engine against the VolatilityMonitor class at the bands R6 adds (0.1, 0.05) and
cs band 0.1, on synthetic Student t3 paths with a change, at a low threshold so that alarms happen [SYNTHETIC].
Seed 9_979_999. Usage: python r6_validate.py <paths>"""
import os, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
os.environ.setdefault("E19_OUT", str(HERE / "pools"))
os.environ.setdefault("E19_SCRATCH", str(HERE / "pools"))
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/experiments")
import e19_real_null as E19  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402
from avmonitor.simulate import change_returns  # noqa: E402

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
CAL, NU, H = 1000, 400, 1950
rng = np.random.default_rng(9_979_999)
p = change_returns(N, CAL + NU + H, nu=CAL + NU, ratio=np.sqrt(2.0), rng=rng, heavy_tail_df=3)
CONFIGS = [("plugin", 0.1, 50.0), ("plugin", 0.05, 50.0), ("cs", 0.1, 50.0), ("plugin", 0.1, 7800.0), ("cs", 0.1, 9750.0)]
res = E19.x_engine_multi(None, E19.ArrayGen(p.r), CONFIGS, CAL, NU + H, horizon_only=True)
bad = 0
for (mode, band, A), R in zip(CONFIGS, res):
    for i in range(N):
        m = VolatilityMonitor(threshold=A, cal_len=CAL, null_band=band, baseline=mode, delta=0.01)
        cls = 0
        for r in p.r[i]:
            ev = m.update(r)
            if ev is not None:
                cls = ev.index - ev.monitored_since + 1
                break
        if cls != int(R["first_t"][i]):
            bad += 1
            print("MISMATCH", mode, band, A, i, "class", cls, "engine", int(R["first_t"][i]))
    print(mode, band, A, "first alarms (engine):", R["first_t"].tolist(), flush=True)
print("mismatches:", bad, "of", N * len(CONFIGS))
