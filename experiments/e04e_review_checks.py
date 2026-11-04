"""E4e. Two checks requested by the independent review (round 1), moved into the repository.

All results [SYNTHETIC], fixed seeds. Ported from the reviewer's scratch script
so that every number quoted in the docs comes from a committed script.

Research questions
  Q4e-1 (review M4). How much does the default plug-in baseline raise the
        one-day false alarm rate of the increase side compared with the TRUE
        baseline? Gaussian, Student t3 (unit variance) and jump (0.5% per step,
        8 sd) returns; A = 7800, one 390-step day, 4000 paths, calibration of
        390 or 1440 steps.
  Q4e-2 (review m5). Recalibration after alarms. With restarts, is the expected
        number of false alarms in m monitored steps at most m / A? (a) Two-sided
        Shiryaev-Roberts statistics with the true baseline, A = 300, 400 paths
        of 60,000 steps; (b) the streaming monitor with plug-in recalibration
        (cal_len 390), Gaussian null, 40 paths of 60,000 steps.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import detectors as D
from avmonitor.monitor import VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE

OUT = Path(__file__).resolve().parents[1] / "results"


def law(name, rng, size):
    if name == "gaussian":
        return rng.standard_normal(size)
    if name == "t3":
        return rng.standard_t(3, size) / np.sqrt(3.0)
    if name == "jumps":
        z = rng.standard_normal(size)
        return z + (rng.random(size) < 0.005) * rng.normal(0, 8.0, size)
    raise ValueError(name)


def q1(A=7800.0, day=390, reps=4000):
    rows = []
    for name in ("gaussian", "t3", "jumps"):
        rng = np.random.default_rng(5)
        big = law(name, rng, 4_000_000)
        tau_true = 3 * np.median(np.abs(big)) * MAD_NORMAL_SCALE
        m0_true = np.mean(np.minimum(big**2, tau_true**2))
        for cal in (390, 1440):
            r = law(name, rng, (reps, cal + day))
            bases = [D.calibrate(r[i, :cal], 3.0) for i in range(reps)]
            tau2 = np.array([b.tau**2 for b in bases])[:, None]
            m0 = np.array([b.m0 for b in bases])[:, None]
            x = np.minimum(r[:, cal:] ** 2, tau2)
            fa_plug = float((D.betting_edetector(x / m0, 1.0, "up", threshold=A).alarm >= 0).mean())
            xt = np.minimum(r[:, cal:] ** 2, tau_true**2)
            fa_true = float((D.betting_edetector(xt / m0_true, 1.0, "up", threshold=A).alarm >= 0).mean())
            rows.append(dict(law=name, cal_len=cal, one_day_fa_plugin=fa_plug, one_day_fa_true_baseline=fa_true,
                             bound=day / A, mc_se=float(np.sqrt(0.05 * 0.95 / reps)), paths=reps))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e04e_plugin_vs_true.csv", index=False)
    print(df.round(4).to_string(index=False))


def q2():
    rows = []
    rng = np.random.default_rng(3)
    P, N, A = 400, 60000, 300.0
    up = np.asarray(D.DEFAULT_UP_GRID)
    big = rng.standard_normal(4_000_000)
    tau2 = (3 * np.median(np.abs(big)) * MAD_NORMAL_SCALE) ** 2
    m0 = np.mean(np.minimum(big**2, tau2))
    b = tau2 / m0
    dn = np.asarray(D.DEFAULT_DOWN_FRACTIONS) / (b - 1)
    Su = np.zeros((P, up.size))
    Sd = np.zeros((P, dn.size))
    count = np.zeros(P)
    for _ in range(N):
        y = np.minimum(rng.standard_normal(P) ** 2, tau2) / m0
        Su = (1 + up * (y[:, None] - 1)) * (Su + 1)
        Sd = (1 + dn * (1 - y[:, None])) * (Sd + 1)
        hit = (Su.mean(1) >= 2 * A) | (Sd.mean(1) >= 2 * A)
        count += hit
        Su[hit] = 0
        Sd[hit] = 0
    rows.append(dict(setting="true baseline, two-sided SR, restarts", A=A, alarms_per_10000_steps=1e4 * count.mean() / N,
                     bound_per_10000_steps=1e4 / A))
    rng = np.random.default_rng(4)
    tot, steps = 0, 0
    for _ in range(40):
        r = rng.standard_normal(60000)
        mon = VolatilityMonitor(threshold=A, cal_len=390)
        ev = mon.run(r)
        tot += len(ev)
        steps += 60000 - 390 * (len(ev) + 1)  # monitored steps (calibration after each alarm is not monitored)
    rows.append(dict(setting="streaming monitor, plug-in recalibration, cal_len 390", A=A,
                     alarms_per_10000_steps=1e4 * tot / steps, bound_per_10000_steps=1e4 / A))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e04e_restarts.csv", index=False)
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    q1()
    q2()
