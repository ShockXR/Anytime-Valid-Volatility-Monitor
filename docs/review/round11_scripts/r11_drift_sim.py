"""Review round 11 (E25), check 1d [SYNTHETIC]: Amendment 1, item 2 says a linear drift "can move a stratum's e-value in
either direction". On the amended planned calendar (41 strata to January 2029) every one-sided stratum has its
controls AFTER the release (19 of 19; 57 kept controls before releases against 115 after), because strata are
filled in release order and earlier strata claim their later candidates first. So a drift common to all windows
moves the unbalanced strata in the same direction. This simulates no release effect at all, with every window's
up-alarm probability following a common exponential trend p(t) = p0(type) * 2^(-(t - t0)/H) (H the half-life of
the alarm rate; negative H: a rising rate), and reports P(M ever >= 20 by 1 January 2029) and E[log e] per stratum.
Strata, entry order and p0 (pooled 2024 to 2026 control rates) as in E25's power part; 20,000 replications.
Usage: nice -n 10 python docs/review/round11_scripts/r11_drift_sim.py"""
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
import r11_strata_lib as L  # noqa: E402
import e24_blind_guard as E24  # noqa: E402

st = P.own_planning_table()
scen = E24.planning_scenarios(pd.read_csv(ROOT / "results/e24_future_planning_inputs.csv"), "S2")
p0 = {t: v[1] for t, v in scen["null (release rate = control rate)"].items()}
t0 = pd.Timestamp("2027-06-01", tz="UTC")
J = len(st)
nmax = max(s["n"] for s in st)
rel_t = np.array([(pd.Timestamp(s["release_utc"]) - t0).days for s in st], float)
ctl_t = np.full((J, nmax), np.nan)
for j, s in enumerate(st):
    for i, (_, a) in enumerate(sorted(s["kept"], key=lambda x: x[1])):
        ctl_t[j, i] = (pd.Timestamp(a) - t0).days
n = np.array([s["n"] for s in st])
base = np.array([p0[s["type"]] for s in st])
R = 20_000
for H in (None, 365.0, 182.0, 91.0, -182.0, -91.0):
    f = (lambda x: np.ones_like(x)) if H is None else (lambda x, H=H: 2.0 ** (-x / H))
    pr = np.clip(base * f(rel_t), 0, 1)
    pc = np.where(np.isnan(ctl_t), np.nan, np.clip(base[:, None] * f(ctl_t), 0, 1))
    rng = np.random.default_rng(11)
    hits = 0
    for b0 in range(0, R, 5000):
        z = (rng.random((5000, J)) < pr).astype(float)
        c = (rng.random((5000, J, nmax)) < np.nan_to_num(pc, nan=-1.0)).sum(2).astype(float)
        p = (z + c) / (n + 1)
        m = np.zeros((5000, J))
        for ps in L.PSI:
            m += np.cumprod(ps**z / (1 - p + ps * p), axis=1) / 4
        hits += int((m >= 20).any(1).sum())
    sh = hits / R
    lab = "no drift" if H is None else (f"falling rate, half-life {H:.0f} days" if H > 0 else
                                        f"rising rate, doubling time {-H:.0f} days")
    print(f"{lab}: P(M ever >= 20 over the {J} strata) = {sh:.4f} (s.e. {np.sqrt(sh*(1-sh)/R):.4f})")
