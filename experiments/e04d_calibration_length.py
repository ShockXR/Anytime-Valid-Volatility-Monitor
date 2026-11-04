"""E4d. How long must the calibration window be?

Research question: the guarantee assumes the true null mean m0. In practice m0
is estimated from a calibration window, and an estimate that is too low makes
the increase detector anti-conservative. How does the one-day false alarm
rate (bound: 390 / A = 5% with the true m0) depend on the calibration length
and on a safety margin of z standard errors? Synthetic Gaussian data.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import detectors as D
from avmonitor.simulate import change_returns

OUT = Path(__file__).resolve().parents[1] / "results"
A, DAY = 7800.0, 390


def main(reps=2000):
    rows = []
    rng = np.random.default_rng(31)
    for cal in (390, 1440, 5000, 20000):
        for z in (0.0, 1.0, 2.0):
            p = change_returns(reps, cal + DAY, rng=rng)
            bases = [D.calibrate(p.r[i, :cal], 3.0, margin_z=z) for i in range(reps)]
            tau2 = np.array([b.tau**2 for b in bases])[:, None]
            up = np.array([b.m0_up for b in bases])[:, None]
            dn = np.array([b.m0_down for b in bases])[:, None]
            x = np.minimum(p.r[:, cal:] ** 2, tau2)
            a_up = D.betting_edetector(x / up, 1.0, "up", threshold=A).alarm
            a_dn = D.betting_edetector(x / dn, 1.0, "down", upper=float(np.max(tau2 / dn)), threshold=A).alarm
            rows.append(dict(cal_len=cal, margin_z=z, up_FA_one_day=float((a_up >= 0).mean()),
                             down_FA_one_day=float((a_dn >= 0).mean()), bound=DAY / A))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e04d_calibration_length.csv", index=False)
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
