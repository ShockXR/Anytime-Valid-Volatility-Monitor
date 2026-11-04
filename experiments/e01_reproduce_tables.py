"""E1. Do the formulas as printed reproduce the published Tables 7 and 8?

Research question: Bapat (2026) does not report the diffusive volatility used
in the simulations or whether the MAD is scaled. If the published RMSEs come
from eqs. (6), (7), (12) and (16) exactly as printed (unnormalised kernel
weights, target = per-step sigma^4), then one value of sigma should reproduce
the whole of Table 7, including the growth of RMSE with h.

Method: simulate eq. (10) with constant per-step sigma s, n = 390, 1000 paths,
Gaussian kernel. Fit s on the 15 rows with sigma_J = 0.01 (kernel estimator
only), then predict the other 75 published numbers out of sample, for raw and
for scaled MAD. All data here are synthetic.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from paper_tables import TABLE7, TABLE8  # noqa: E402

from avmonitor import robust as R  # noqa: E402
from avmonitor.simulate import paper_returns  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
HS = (10, 20, 30, 40, 50)


def rmse_rows(s, p, sj, noise_var=0.0, reps=1000, seed=0, normalised=False):
    rng = np.random.default_rng(seed)
    r, target = paper_returns(reps, 390, s, p, 0.0, sj, noise_var, rng)
    out = {}
    for h in HS:
        if normalised:
            k = R.local_quarticity(r, h, T=390).mean(axis=1)
            rows = {"kernel": k, "robust_raw": R.kernel_iq(r, h, T=390, c=3, mad_scaled=False) / 390,
                    "robust_scaled": R.kernel_iq(r, h, T=390, c=3, mad_scaled=True) / 390}
        else:
            rows = {
                "kernel": R.kernel_quarticity_printed(r, h),
                "robust_raw": R.robust_kernel_quarticity_printed(r, h, 3.0, mad_scaled=False),
                "robust_scaled": R.robust_kernel_quarticity_printed(r, h, 3.0, mad_scaled=True),
            }
        for name, est in rows.items():
            out[(h, name)] = float(np.sqrt(np.mean((est - target) ** 2)))
    return out


def fit_sigma():
    grid = np.round(np.arange(0.012, 0.032, 0.0005), 5)
    best = None
    for s in grid:
        err = []
        for p in (0.005, 0.01, 0.02):
            sim = rmse_rows(s, p, 0.01, reps=400, seed=1)
            for h in HS:
                pub = TABLE7[(p, 0.01, h)][0]
                err.append(np.log(sim[(h, "kernel")] / pub) ** 2)
        e = float(np.mean(err))
        if best is None or e < best[1]:
            best = (s, e)
    return best


def main():
    s_hat, fit_err = fit_sigma()
    print(f"fitted per-step sigma = {s_hat}  (sigma^4 = {s_hat**4:.3e}), mean sq log error {fit_err:.4f}")
    rows = []
    for (p, sj, h), (pub_k, pub_r) in TABLE7.items():
        sim = rmse_rows(s_hat, p, sj, reps=1000, seed=2)
        rows.append(dict(table=7, p_jump=p, sigma_J=sj, noise_var=0.0, h=h, pub_kernel=pub_k, pub_robust=pub_r,
                         sim_kernel=sim[(h, "kernel")], sim_robust_rawMAD=sim[(h, "robust_raw")],
                         sim_robust_scaledMAD=sim[(h, "robust_scaled")], in_fit=(sj == 0.01)))
    for (nv, h), (pub_k, pub_r) in TABLE8.items():
        for label, var in (("as_variance", nv), ("as_sd", nv**2)):
            sim = rmse_rows(s_hat, 0.01, 0.05, noise_var=var, reps=1000, seed=3)
            rows.append(dict(table=8, p_jump=0.01, sigma_J=0.05, noise_var=nv, noise_reading=label, h=h,
                             pub_kernel=pub_k, pub_robust=pub_r, sim_kernel=sim[(h, "kernel")],
                             sim_robust_rawMAD=sim[(h, "robust_raw")], sim_robust_scaledMAD=sim[(h, "robust_scaled")],
                             in_fit=False))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e01_reproduce_tables.csv", index=False)

    t7 = df[(df.table == 7) & (~df.in_fit)]
    def med_ratio(a, b):
        return float(np.median(a / b)), float(np.mean(np.abs(np.log(a / b)) < np.log(1.5)))
    print("Out-of-sample Table 7 rows (sigma_J = 0.05 and 0.10):")
    print("  kernel            median sim/pub, share within x1.5:", med_ratio(t7.sim_kernel, t7.pub_kernel))
    print("  robust raw MAD    median sim/pub, share within x1.5:", med_ratio(t7.sim_robust_rawMAD, t7.pub_robust))
    print("  robust scaled MAD median sim/pub, share within x1.5:", med_ratio(t7.sim_robust_scaledMAD, t7.pub_robust))
    t7in = df[(df.table == 7) & (df.in_fit)]
    print("In-fit rows robust raw/scaled:", med_ratio(t7in.sim_robust_rawMAD, t7in.pub_robust), med_ratio(t7in.sim_robust_scaledMAD, t7in.pub_robust))
    for lab in ("as_variance", "as_sd"):
        t8 = df[(df.table == 8) & (df.noise_reading == lab) & (df.noise_var > 0)]
        print(f"Table 8 noise read {lab}: kernel", med_ratio(t8.sim_kernel, t8.pub_kernel), "robust raw", med_ratio(t8.sim_robust_rawMAD, t8.pub_robust))

    # Counterfactual: normalised estimators on the same design
    cf = []
    for p in (0.005, 0.01, 0.02):
        for sj in (0.01, 0.05, 0.10):
            sim = rmse_rows(s_hat, p, sj, reps=1000, seed=2, normalised=True)
            for h in HS:
                cf.append(dict(p_jump=p, sigma_J=sj, h=h, rel_rmse_kernel=sim[(h, "kernel")] / s_hat**4,
                               rel_rmse_robust_rawMAD=sim[(h, "robust_raw")] / s_hat**4,
                               rel_rmse_robust_scaledMAD=sim[(h, "robust_scaled")] / s_hat**4))
    cf = pd.DataFrame(cf)
    cf.to_csv(OUT / "e01_normalised_counterfactual.csv", index=False)
    print(cf.groupby("sigma_J")[["rel_rmse_kernel", "rel_rmse_robust_rawMAD", "rel_rmse_robust_scaledMAD"]].median())
    with open(OUT / "e01_summary.txt", "w") as f:
        f.write(f"fitted_sigma={s_hat}\nsigma4={s_hat**4:.4e}\nfit_err={fit_err}\n")


if __name__ == "__main__":
    main()
