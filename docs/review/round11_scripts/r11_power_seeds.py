"""Review round 11 (E25): the independent power simulation of r11_power_check.py with 100,000 replications per
scenario and a separate seed per scenario, to see whether the small positive differences of that check are noise.
Usage: nice -n 10 python docs/review/round11_scripts/r11_power_seeds.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT / "src"))
import r11_power_check as P  # noqa: E402
import e24_blind_guard as E24  # noqa: E402

st = P.own_planning_table()
scen = E24.planning_scenarios(pd.read_csv(ROOT / "results/e24_future_planning_inputs.csv"), "S2")
pw = pd.read_csv(ROOT / "results/e25_power.csv")
ent = np.array([pd.Timestamp(s["enters"]).value for s in st])
R = 100_000
for j, (sname, rates) in enumerate(scen.items()):
    fo = P.own_sim(st, rates, R, 7000 + j)
    cd = np.where(fo > 0, ent[np.maximum(fo, 1) - 1], np.iinfo(np.int64).max)
    r = pw[(pw.scenario == sname) & pw.procedure.str.startswith("E25")].iloc[0]
    msg = []
    for d in ("2027-06-01", "2028-01-01", "2029-01-01"):
        a = float((cd <= pd.Timestamp(d, tz="UTC").value).mean())
        b = float(r[f"P_cross_by_{d}"])
        se = np.sqrt(a * (1 - a) / R + b * (1 - b) / 20000)
        msg.append(f"{d}: own {a:.4f} CSV {b:.4f} z {(a - b) / se:+.2f}")
    print(f"{sname}: " + "; ".join(msg))
