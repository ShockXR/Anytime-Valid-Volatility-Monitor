"""E4c. Living with stochastic volatility: a null band.

Research question: when volatility drifts slowly inside a regime (slow SV),
the strict null "variance never exceeds the baseline" is false, so the
e-detector alarms on ordinary drift (E4a). Does widening the null to
"truncated variance within +/- b of baseline" restore a low false alarm rate
while keeping detection of large regime changes (sd x2 and x0.5)?
Synthetic data, slow SV with stationary log-sd 0.15 and half-life 390 steps.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import detectors as D
from avmonitor.simulate import change_returns

OUT = Path(__file__).resolve().parents[1] / "results"
A, CAL = 7800.0, 2000


def run(band, ratio, seed):
    rng = np.random.default_rng(seed)
    nu, n = (500, 3500) if ratio != 1.0 else (None, 20000)
    p = change_returns(300, CAL + n, nu=None if nu is None else CAL + nu, ratio=ratio, sv_xi=0.15, sv_halflife=390, rng=rng)
    alarms_up, alarms_dn = [], []
    for i in range(p.r.shape[0]):
        b = D.calibrate(p.r[i, :CAL], c=3.0, null_band=band)
        x = np.minimum(p.r[i, CAL:] ** 2, b.tau**2)
        alarms_up.append(D.betting_edetector(x, b.m0_up, "up", threshold=A).alarm[0])
        alarms_dn.append(D.betting_edetector(x, b.m0_down, "down", upper=b.tau**2, threshold=A).alarm[0])
    return np.array(alarms_up), np.array(alarms_dn), nu, n


def main():
    rows = []
    for band in (0.0, 0.25, 0.5):
        up, dn, _, n = run(band, 1.0, 5)
        rows.append(dict(band=band, case="no change", up_FA_within_390=float(((up >= 0) & (up < 390)).mean()),
                         down_FA_within_390=float(((dn >= 0) & (dn < 390)).mean()),
                         up_FA_within_20000=float((up >= 0).mean()), down_FA_within_20000=float((dn >= 0).mean())))
        for ratio in (2.0, 0.5):
            up, dn, nu, n = run(band, ratio, 6)
            a = up if ratio > 1 else dn
            s = D.summarize_alarms(a, nu, n)
            rows.append(dict(band=band, case=f"sd x{ratio}", **s))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e04c_null_band.csv", index=False)
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
