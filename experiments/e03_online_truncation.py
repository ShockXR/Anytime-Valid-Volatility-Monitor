"""E3. Can the truncation be made online without losing robustness?

Research question: eq. (11) uses the MAD of the whole sample, which is not
available in real time and is wrong for both regimes when volatility changes.
We compare local (spot) variance estimates from the one-sided Gaussian kernel
(h = 30) after truncation with (a) no truncation, (b) full-sample MAD,
(c) rolling 390-step MAD and (d) exponentially weighted MAD (half-life 390),
all with c = 3 and scaled MAD, on paths whose volatility doubles (sd x2) at
step 2000 and which contain jumps. Synthetic data, per-step units.

Metric: median absolute log error of the spot variance estimate against the
true per-step variance, before the change, in the 300 steps after it, and
long after it; and the share of jump steps whose window estimate is inflated
by more than 50 percent.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import robust as R
from avmonitor.simulate import change_returns

OUT = Path(__file__).resolve().parents[1] / "results"


def main(reps=200, n=4000, nu=2000, h=30, c=3.0):
    rng = np.random.default_rng(7)
    rows = []
    for ratio in (2.0, 0.5):
        p = change_returns(reps, n, nu=nu, ratio=ratio, p_jump=0.005, jump_scale=10.0, rng=rng)
        truth = p.sigma[:, h:] ** 2
        variants = {}
        variants["none"] = p.r
        tau_full = R.mad_threshold(p.r, c, scaled=True)
        variants["full-sample MAD"] = R.truncate(p.r, tau_full)
        roll = np.vstack([R.rolling_mad_zero(row, 390, min_periods=100) for row in p.r])
        roll = np.where(np.isnan(roll), np.nanmedian(roll, axis=1, keepdims=True), roll)
        variants["rolling MAD (390)"] = np.clip(p.r, -c * roll, c * roll)
        ew = np.vstack([R.ew_quantile_scale(row, halflife=390) for row in p.r])
        variants["EW MAD (half-life 390)"] = np.clip(p.r, -c * ew, c * ew)
        for name, rt in variants.items():
            est = R.local_variance(rt, h, T=n) * 1.0  # T = n gives per-step units
            if name != "none":
                est = est / R.truncated_moment_factor(c, 2)
            err = np.abs(np.log(est / truth))
            t = np.arange(h, n)
            rows.append(dict(ratio=ratio, truncation=name,
                             err_before=float(np.median(err[:, (t > 500) & (t < nu)])),
                             err_after_300=float(np.median(err[:, (t >= nu) & (t < nu + 300)])),
                             err_late=float(np.median(err[:, t > nu + 1000])),
                             p95_err_before=float(np.quantile(err[:, (t > 500) & (t < nu)], 0.95))))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e03_online_truncation.csv", index=False)
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
