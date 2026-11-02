"""E2. Truncation bias of tau = c * MAD, and a fair comparison of correctly scaled IQ estimators.

Research questions
  Q2a. How much of the fourth moment does tau = c * MAD remove, for raw and
       scaled MAD, and does the loss vanish as n grows (as Theorem 1 needs)?
  Q2b. With time-varying volatility, does a single full-sample MAD truncate
       high-volatility periods more?
  Q2c. Once every estimator is put on the IQ scale, how do they rank under
       jumps, return-level noise (the paper's eq. 13) and price-level noise?
  Q2d. Is the aggregated kernel estimator (eq. 16) distinguishable from plain
       truncated RQ?
All data are synthetic.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import robust as R
from avmonitor.simulate import sv_jump_returns

OUT = Path(__file__).resolve().parents[1] / "results"
rng = np.random.default_rng(2026)


def q2a():
    rows = []
    for c in (2.0, 2.5, 3.0, 4.0, 5.0):
        for scaled in (False, True):
            a = R.effective_cutoff(c, scaled)
            row = dict(c=c, mad="scaled" if scaled else "raw", cutoff_sd=a,
                       share_clipped=float(2 * (1 - __import__("scipy").stats.norm.cdf(a))),
                       kappa4_theory=R.truncated_moment_factor(a, 4), kappa2_theory=R.truncated_moment_factor(a, 2))
            for n in (390, 1440, 23400):
                r = rng.standard_normal((300, n)) / np.sqrt(n)
                row[f"kappa4_sim_n{n}"] = float(R.truncated_rq_mad(r, c, mad_scaled=scaled).mean())
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e02a_truncation_factor.csv", index=False)
    print(df.round(4).to_string(index=False))


def q2b():
    rows = []
    for xi in (0.0, 0.5, 1.0, 1.5):
        p = sv_jump_returns(500, 390, sigma_bar=0.2, kappa=5.0, xi=xi, rng=rng) if xi > 0 else None
        if p is None:
            r = rng.standard_normal((500, 390)) * 0.2 / np.sqrt(390)
            iq = np.full(500, 0.2**4)
        else:
            r, iq = p.r, p.iq
        for scaled in (False, True):
            est = R.truncated_rq_mad(r, 3.0, mad_scaled=scaled, bias_correct=True)
            rows.append(dict(vol_of_vol=xi, mad="scaled" if scaled else "raw",
                             mean_ratio_bias_corrected=float(np.mean(est / iq)),
                             mean_ratio_RQ=float(np.mean(R.realized_quarticity(r) / iq))))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e02b_sv_global_mad.csv", index=False)
    print(df.round(3).to_string(index=False))


ESTIMATORS = {
    "RQ": lambda r: R.realized_quarticity(r),
    "TPQ": lambda r: R.tripower_quarticity(r),
    "QPQ": lambda r: R.quadpower_quarticity(r),
    "MinRQ": lambda r: R.min_rq(r),
    "MedRQ": lambda r: R.med_rq(r),
    "TRQ": lambda r: R.threshold_rq(r),
    "Kernel": lambda r: R.kernel_iq(r, 30),
    "RobustKernel_rawMAD": lambda r: R.kernel_iq(r, 30, c=3, mad_scaled=False),
    "RobustKernel_scaledMAD": lambda r: R.kernel_iq(r, 30, c=3, mad_scaled=True),
    "RobustKernel_scaledMAD_bc": lambda r: R.kernel_iq(r, 30, c=3, mad_scaled=True, bias_correct=True),
    "TruncRQ_scaledMAD_bc": lambda r: R.truncated_rq_mad(r, 3, mad_scaled=True, bias_correct=True),
}


def q2c():
    scenarios = {
        "no jumps, no noise": dict(),
        "jumps 4/day, 5 sd": dict(jump_intensity=4, jump_size_sd=5 * 0.2 / np.sqrt(390)),
        "jumps 4/day, 15 sd": dict(jump_intensity=4, jump_size_sd=15 * 0.2 / np.sqrt(390)),
        "price noise, ratio 0.25": dict(noise_sd=0.25 * 0.2 / np.sqrt(390), noise_level="price"),
        "price noise, ratio 0.5": dict(noise_sd=0.5 * 0.2 / np.sqrt(390), noise_level="price"),
        "return noise (eq. 13), ratio 0.5": dict(noise_sd=0.5 * 0.2 / np.sqrt(390), noise_level="return"),
        "jumps 15 sd + price noise 0.25": dict(jump_intensity=4, jump_size_sd=15 * 0.2 / np.sqrt(390),
                                               noise_sd=0.25 * 0.2 / np.sqrt(390)),
    }
    rows, corr_rows = [], []
    for name, kw in scenarios.items():
        p = sv_jump_returns(1000, 390, sigma_bar=0.2, kappa=5.0, xi=0.5, rng=rng, **kw)
        for est_name, f in ESTIMATORS.items():
            e = f(p.r) / p.iq
            rows.append(dict(scenario=name, estimator=est_name, mean_ratio=float(e.mean()),
                             rel_rmse=float(np.sqrt(np.mean((e - 1) ** 2))), median_ratio=float(np.median(e))))
        a = R.kernel_iq(p.r, 30, c=3, mad_scaled=True)
        b = R.truncated_rq_mad(p.r, 3, mad_scaled=True)
        corr_rows.append(dict(scenario=name, corr=float(np.corrcoef(a, b)[0, 1]), mean_ratio=float(np.mean(a / b))))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e02c_fair_comparison.csv", index=False)
    piv = df.pivot(index="estimator", columns="scenario", values="rel_rmse")
    print("Relative RMSE (RMSE / IQ):")
    print(piv.round(3).to_string())
    piv2 = df.pivot(index="estimator", columns="scenario", values="mean_ratio")
    print("Mean estimate / IQ:")
    print(piv2.round(3).to_string())
    cr = pd.DataFrame(corr_rows)
    cr.to_csv(OUT / "e02d_kernel_vs_truncated_rq.csv", index=False)
    print(cr.round(4).to_string(index=False))


if __name__ == "__main__":
    q2a()
    q2b()
    q2c()
