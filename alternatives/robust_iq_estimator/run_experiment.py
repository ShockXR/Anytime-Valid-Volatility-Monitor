"""Experiments for the robust local plug-in IQ estimator (alternative A5, agent ALT3).

Research questions
  Q1 [SYNTHETIC + theory] Efficiency. What are the asymptotic variances (AVAR, in
     units of int sigma^8) of RQ, TPQ, QPQ, MedRQ, MinRQ, the MAD-truncated
     kernel estimator and the new block / kernel plug-in estimators, and how
     close does each come to the efficient bound 8 int sigma^8? Closed forms
     are checked against constant-volatility simulations.
  Q2 [SYNTHETIC] Consistency of truncation rules. With log-OU stochastic
     volatility, intraday seasonality and jumps of fixed size, how do relative
     bias and RMSE change from n = 390 to 1440 to 23400 for: Bapat's fixed
     c = 3 x global MAD (clipping, with and without the kappa_4 correction),
     a Mancini threshold, and local thresholds with a fixed (c = 3) or growing
     (sqrt(2 log n)) cut-off, clipping or dropping? Which estimators converge,
     and at what rate (slope of log RMSE on log n)?
  Q3 [SYNTHETIC] Bandwidth. Does the paper's cross-validation (eqs. 17, 18)
     pick its bandwidth from a scale artefact, does the corrected target pick a
     sensible one, and how sensitive is the plug-in estimator to block length?
  Q4 [SYNTHETIC] Full comparison: relative bias, relative RMSE and coverage of
     95% intervals (plain and log scale) built from a robust estimate of
     int sigma^8, plus the accuracy of that int sigma^8 estimate.
  Q5 [REAL] Binance BTCUSDT and DOGEUSDT, 1 July to 31 December 2025, one day
     = 1440 one-minute returns. Do the new estimators agree with MedRQ and
     MinRQ day by day better than the MAD-truncated ones? Same at 5 minutes.

Usage: python run_experiment.py [theory|sim|cv|real|figures|all]
All seeds are fixed. Uses at most two processes. Outputs go to results/.
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import riq  # noqa: E402
from avmonitor import robust as R  # noqa: E402
from avmonitor.simulate import sv_jump_returns  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

SIGMA_BAR = 0.02  # daily volatility (T = 1 day); jump sizes are multiples of it
NS = (390, 1440, 23400)
SUBSTEPS = {390: 10, 1440: 4, 23400: 1}
REPS = 500
SCENARIOS = {
    "S1 no jumps": dict(xi=0.5),
    "S2 small jumps": dict(xi=0.5, lam=4.0, js=0.1),
    "S3 medium jumps": dict(xi=0.5, lam=4.0, js=0.3),
    "S4 large jumps": dict(xi=0.5, lam=4.0, js=1.0),
    "S5 strong SV + intraday pattern": dict(xi=1.0, season=0.6),
    "S6 strong SV + pattern + medium jumps": dict(xi=1.0, season=0.6, lam=4.0, js=0.3),
}
EST_ORDER = ["RQ", "TPQ", "MedRQ", "MinRQ", "TRQ-M", "RK", "RK-k", "TRQ-U",
             "PI-G3", "PI-L3", "PI-LUc", "PI-LU", "PI-LU-raw", "KPI-L3"]
K_LIST = (2, 5, 10, 20, 40, 80, 160, 320)
CV_H = (5, 10, 20, 30, 40, 50, 75, 100, 150)


# ----------------------------------------------------------------------------
# Estimator set
# ----------------------------------------------------------------------------


def estimate_all(r: np.ndarray, T: float = 1.0):
    """All estimators on the IQ scale for paths in the rows of ``r`` (no NaN)."""
    r = np.atleast_2d(np.asarray(r, float))
    n = r.shape[-1]
    est = {
        "RQ": R.realized_quarticity(r, T),
        "TPQ": R.tripower_quarticity(r, T),
        "MedRQ": R.med_rq(r, T),
        "MinRQ": R.min_rq(r, T),
        "TRQ-M": R.threshold_rq(r, T),
        "RK": R.kernel_iq(r, h=30, T=T, c=3.0, mad_scaled=True),
        "RK-k": R.kernel_iq(r, h=30, T=T, c=3.0, mad_scaled=True, bias_correct=True),
    }
    s = riq.local_scale(r)
    aU = riq.universal_cutoff(n)
    tauG, _ = riq.thresholds(r, rule="fixed", c=3.0, scale="global")
    est["TRQ-U"] = riq.truncated_rq_local(r, aU * s, aU, T, "drop")
    est["PI-G3"] = riq.block_plugin(r, tauG, 3.0, T=T, mode="clip")["iq"]
    pl3 = riq.block_plugin(r, 3.0 * s, 3.0, T=T, mode="clip")
    est["PI-L3"] = pl3["iq"]
    est["PI-LUc"] = riq.block_plugin(r, aU * s, aU, T=T, mode="clip")["iq"]
    plu = riq.block_plugin(r, aU * s, aU, T=T, mode="drop")
    est["PI-LU"] = plu["iq"]
    est["PI-LU-raw"] = riq.block_plugin(r, aU * s, aU, T=T, mode="drop", correct=False)["iq"]
    est["KPI-L3"] = riq.kernel_plugin(r, 3.0 * s, 3.0, h=30, T=T)
    extras = dict(se_LU=plu["se"], iq8_LU=plu["iq8"], drop_LU=plu["dropped"],
                  se_L3=pl3["se"], iq8_L3=pl3["iq8"], drop_L3=pl3["dropped"], scale=s, aU=aU)
    return est, extras


def simulate(scen: dict, n: int, reps: int, rng: np.random.Generator):
    return sv_jump_returns(
        reps, n, T=1.0, sigma_bar=SIGMA_BAR, kappa=5.0, xi=scen.get("xi", 0.5),
        jump_intensity=scen.get("lam", 0.0), jump_size_sd=scen.get("js", 0.0) * SIGMA_BAR,
        seasonal_amplitude=scen.get("season", 0.0), substeps=SUBSTEPS[n], rng=rng,
    )


# ----------------------------------------------------------------------------
# Q1: theory tables
# ----------------------------------------------------------------------------


def part_theory():
    """Closed-form constants and constant-volatility checks (Q1, plus jump residuals and CV theory)."""
    t0 = time.time()
    rows = []
    # Constant volatility simulation: sigma = 1, T = 1, n = 1440.
    n, reps = 1440, 8000
    rng = np.random.default_rng(11)
    sims = {}
    blocks = {}
    for b in range(8):
        r = rng.standard_normal((reps // 8, n)) / np.sqrt(n)
        est, _ = estimate_all(r)
        for key, v in est.items():
            sims.setdefault(key, []).append(v)
        for k in (2, 5, 10, 20, 40, 160):
            blocks.setdefault(k, []).append(riq.block_plugin(r, k=k)["iq"])
        s = riq.local_scale(r)
        blocks.setdefault("oracle-clip3", []).append(riq.block_plugin(r, np.full(r.shape, 3 / np.sqrt(n)), 3.0, mode="clip")["iq"])
        blocks.setdefault("feasible-clip3-k38", []).append(riq.block_plugin(r, 3.0 * s, 3.0, mode="clip")["iq"])
    sims = {k: np.concatenate(v) for k, v in sims.items()}
    blocks = {k: np.concatenate(v) for k, v in blocks.items()}

    def sim_avar(x):
        return n * np.var(x, ddof=1), n * np.var(x, ddof=1) * np.sqrt(2.0 / (len(x) - 1)), np.mean(x)

    closed = {
        "RQ": 32 / 3, "TPQ": riq.avar_tpq(), "MedRQ": np.nan, "MinRQ": np.nan, "TRQ-M": 32 / 3,
        "RK": np.nan, "RK-k": riq.avar_truncated_rq(3.0, "clip"), "TRQ-U": 32 / 3,
        "PI-G3": riq.avar_block(38), "PI-L3": riq.avar_block(38), "PI-LUc": riq.avar_block(38),
        "PI-LU": riq.avar_block(38), "PI-LU-raw": np.nan, "KPI-L3": np.nan,
    }
    notes = {
        "RQ": "closed form 32/3 (k = 1 block)", "TPQ": "closed form, 2-dependent summands",
        "MedRQ": "simulation only", "MinRQ": "simulation only",
        "TRQ-M": "no truncation limit 32/3", "RK": "biased (no kappa correction), AVAR of the biased estimator",
        "RK-k": "closed form with oracle cut-off a = 3: (E g^8 / E g^4^2 - 1)",
        "TRQ-U": "no truncation limit 32/3", "PI-G3": "Gaussian block formula 8(k+3)/(k+2), k = 38",
        "PI-L3": "Gaussian block formula, k = 38", "PI-LUc": "Gaussian block formula, k = 38",
        "PI-LU": "Gaussian block formula, k = 38", "PI-LU-raw": "biased by about 2/k",
        "KPI-L3": "overlapping Gaussian kernel h = 30, simulation only",
    }
    for key in EST_ORDER:
        av, se, mean = sim_avar(sims[key])
        rows.append(dict(estimator=key, avar_closed_form=closed[key], avar_simulated=av, avar_sim_se=se,
                         mean_over_truth=mean, efficiency_vs_bound=8.0 / av, note=notes[key]))
    rows.append(dict(estimator="QPQ", avar_closed_form=riq.avar_qpq(), avar_simulated=np.nan, avar_sim_se=np.nan,
                     mean_over_truth=np.nan, efficiency_vs_bound=8.0 / riq.avar_qpq(), note="closed form, 3-dependent summands"))
    for k in (2, 5, 10, 20, 40, 160):
        av, se, mean = sim_avar(blocks[k])
        rows.append(dict(estimator=f"block plug-in, no truncation, k = {k}", avar_closed_form=riq.avar_block(k),
                         avar_simulated=av, avar_sim_se=se, mean_over_truth=mean, efficiency_vs_bound=8.0 / av,
                         note="closed form 8(k+3)/(k+2)"))
    av, se, mean = sim_avar(blocks["oracle-clip3"])
    rows.append(dict(estimator="block plug-in, clip at a = 3 with the TRUE scale (oracle), k = 38",
                     avar_closed_form=riq.avar_block_exact(38, 3.0), avar_simulated=av, avar_sim_se=se,
                     mean_over_truth=mean, efficiency_vs_bound=8.0 / av,
                     note="closed form k (E[(S/k)^4]/E[(S/k)^2]^2 - 1) for clipped normals; below 8 only because it uses the true sigma"))
    av, se, mean = sim_avar(blocks["feasible-clip3-k38"])
    rows.append(dict(estimator="block plug-in, clip at 3 x local scale (feasible = PI-L3)", avar_closed_form=np.nan,
                     avar_simulated=av, avar_sim_se=se, mean_over_truth=mean, efficiency_vs_bound=8.0 / av,
                     note="same as PI-L3"))
    rows.append(dict(estimator="efficient bound", avar_closed_form=8.0, avar_simulated=np.nan, avar_sim_se=np.nan,
                     mean_over_truth=np.nan, efficiency_vs_bound=1.0, note="2 g'(c)^2 c^2 with g(c) = c^2"))
    pd.DataFrame(rows).to_csv(OUT / "q1_avar.csv", index=False)

    # Truncation constants for the cut-offs in use.
    trows = []
    cut = [("raw MAD, c = 3", 3 * stats.norm.ppf(0.75)), ("scaled MAD, c = 3", 3.0)]
    cut += [(f"universal sqrt(2 log n), n = {m}", riq.universal_cutoff(m)) for m in NS]
    cut += [(f"Mancini 4 n^0.01 (TRQ-M), n = {m}", 4 * m**0.01) for m in NS]
    for name, a in cut:
        p = 2 * stats.norm.sf(a)
        T = riq.trunc_even_moments(a, 2)
        clip = riq.square_moments(a, "clip")
        trows.append(dict(rule=name, a=a, p_beyond=p, expected_beyond_per_390=390 * p, expected_beyond_per_23400=23400 * p,
                          kappa2_clip=clip[0], kappa4_clip=clip[1] / 3, var_kept_drop=T[1], quart_kept_drop=T[2] / 3))
    pd.DataFrame(trows).to_csv(OUT / "q1_truncation_constants.csv", index=False)

    # Jump residuals: closed form and a simulation check (constant sigma, oracle threshold).
    jrows = []
    for n_, a in [(400, 3.0), (390, 3.0), (390, riq.universal_cutoff(390)), (1440, riq.universal_cutoff(1440)), (23400, riq.universal_cutoff(23400))]:
        k = 20 if n_ == 400 else riq.default_block(n_)
        jrows.append(dict(n=n_, a=a, k=k, plugin_clip_residual=riq.clipped_jump_residual(a, k),
                          fourth_power_clip_residual=riq.clipped_jump_residual(a, 1),
                          plugin_rel_bias_per_jump=riq.clipped_jump_residual(a, k) / n_,
                          fourth_power_rel_bias_per_jump=riq.clipped_jump_residual(a, 1) / n_, drop_residual=0.0, source="closed form"))
    rng = np.random.default_rng(12)
    n_, a, k, reps = 400, 3.0, 20, 40000
    z = rng.standard_normal((reps, n_)) / np.sqrt(n_)
    tau = np.full(z.shape, a / np.sqrt(n_))
    zj = z.copy()
    zj[:, 195] += 50.0 / np.sqrt(n_)
    d_pi = riq.block_plugin(zj, tau, a, k=k, mode="clip")["iq"] - riq.block_plugin(z, tau, a, k=k, mode="clip")["iq"]
    d_rq = riq.truncated_rq_local(zj, tau, a, mode="clip") - riq.truncated_rq_local(z, tau, a, mode="clip")
    d_dr = riq.block_plugin(zj, tau, a, k=k, mode="drop")["iq"] - riq.block_plugin(z, tau, a, k=k, mode="drop")["iq"]
    jrows.append(dict(n=n_, a=a, k=k, plugin_clip_residual=d_pi.mean() * n_, fourth_power_clip_residual=d_rq.mean() * n_,
                      plugin_rel_bias_per_jump=d_pi.mean(), fourth_power_rel_bias_per_jump=d_rq.mean(),
                      drop_residual=d_dr.mean() * n_, source=f"simulation, one 50 sd jump, {reps} paths, oracle threshold"))
    for b in (1.0, 2.0, 3.0):
        jrows.append(dict(n=np.nan, a=np.inf, k=20, plugin_clip_residual=riq.undetected_jump_residual(b, 20),
                          fourth_power_clip_residual=riq.undetected_jump_residual(b, 1), plugin_rel_bias_per_jump=np.nan,
                          fourth_power_rel_bias_per_jump=np.nan, drop_residual=np.nan,
                          source=f"closed form, undetected jump of b = {b:g} sd (not truncated)"))
    zj = z.copy()
    zj[:, 195] += 2.0 / np.sqrt(n_)
    d_pi = riq.block_plugin(zj, k=k)["iq"] - riq.block_plugin(z, k=k)["iq"]
    d_rq = R.realized_quarticity(zj) - R.realized_quarticity(z)
    jrows.append(dict(n=n_, a=np.inf, k=k, plugin_clip_residual=d_pi.mean() * n_, fourth_power_clip_residual=d_rq.mean() * n_,
                      plugin_rel_bias_per_jump=d_pi.mean(), fourth_power_rel_bias_per_jump=d_rq.mean(), drop_residual=np.nan,
                      source=f"simulation, undetected 2 sd jump, {reps} paths"))
    pd.DataFrame(jrows).to_csv(OUT / "q1_jump_residuals.csv", index=False)

    # Expected CV under constant volatility (eq. 17 printed versus corrected), scaled MAD c = 3.
    crows = []
    for h in range(2, 151):
        p, c = riq.printed_cv_expectation(h, 3.0)
        crows.append(dict(h=h, printed=p, corrected=c, raw_weight_sum_lags_1_to_h=float(R.kernel_weights(h, normalise=False)[:-1].sum())))
    pd.DataFrame(crows).to_csv(OUT / "q3_cv_theory.csv", index=False)
    print(f"theory done in {time.time() - t0:.0f} s")


# ----------------------------------------------------------------------------
# Q2 and Q4: simulation grid
# ----------------------------------------------------------------------------


def _summ(err: np.ndarray):
    err = err[np.isfinite(err)]
    m = len(err)
    return dict(rel_bias=err.mean(), rel_bias_se=err.std(ddof=1) / np.sqrt(m), rel_sd=err.std(ddof=1),
                rel_rmse=np.sqrt(np.mean(err**2)), median_rel_err=np.median(err), reps=m)


def run_cell(args):
    si, scen_name, n = args
    scen = SCENARIOS[scen_name]
    seed = 1000 + 10 * si + NS.index(n)
    rng = np.random.default_rng(seed)
    t0 = time.time()
    batch = 100 if n == 23400 else REPS
    acc = {k: [] for k in EST_ORDER}
    acc_k = {k: [] for k in K_LIST if k <= n // 2} if si in (0, 2, 4) else {}
    truth_iq, truth_iq8, ex = [], [], {}
    jumps = []
    for _ in range(REPS // batch):
        paths = simulate(scen, n, batch, rng)
        est, extras = estimate_all(paths.r)
        for k, v in est.items():
            acc[k].append(v)
        for key in ("se_LU", "iq8_LU", "drop_LU", "se_L3", "iq8_L3", "drop_L3"):
            ex.setdefault(key, []).append(extras[key])
        truth_iq.append(paths.iq)
        truth_iq8.append(np.sum(paths.spot_var**4, axis=-1) / n)
        jumps.append(paths.jump_count)
        if si in (0, 2, 4):  # block length study
            aU, s = extras["aU"], extras["scale"]
            for k in acc_k:
                acc_k[k].append(riq.block_plugin(paths.r, aU * s, aU, k=k, mode="drop")["iq"])
    iq = np.concatenate(truth_iq)
    iq8 = np.concatenate(truth_iq8)
    rows, cov, blk = [], [], []
    for k in EST_ORDER:
        e = np.concatenate(acc[k]) / iq - 1
        rows.append(dict(scenario=scen_name, n=n, estimator=k, **_summ(e)))
    ex = {k: np.concatenate(v) for k, v in ex.items()}
    z = stats.norm.ppf(0.975)
    for lab, est_key, se_key, i8_key in (("PI-LU", "PI-LU", "se_LU", "iq8_LU"), ("PI-L3", "PI-L3", "se_L3", "iq8_L3")):
        x = np.concatenate(acc[est_key])
        se = ex[se_key]
        lo, hi = x - z * se, x + z * se
        llo, lhi = x * np.exp(-z * se / x), x * np.exp(z * se / x)
        e8 = ex[i8_key] / iq8 - 1
        cov.append(dict(scenario=scen_name, n=n, estimator=lab,
                        coverage_plain=np.mean((iq >= lo) & (iq <= hi)), coverage_log=np.mean((iq >= llo) & (iq <= lhi)),
                        below_plain=np.mean(iq < lo), above_plain=np.mean(iq > hi),
                        mean_rel_halfwidth=np.mean(z * se / x), se_over_sd=np.mean(se / iq) / np.std(x / iq - 1),
                        iq8_rel_bias=e8.mean(), iq8_rel_rmse=np.sqrt(np.mean(e8**2)),
                        share_truncated=ex["drop_LU" if lab == "PI-LU" else "drop_L3"].mean(), reps=len(x)))
    for k, v in acc_k.items():
        e = np.concatenate(v) / iq - 1
        blk.append(dict(scenario=scen_name, n=n, k=k, **_summ(e)))
    mean_jumps = float(np.concatenate(jumps).mean())
    print(f"{scen_name:42s} n={n:6d} done in {time.time() - t0:5.1f} s (mean jumps {mean_jumps:.2f})", flush=True)
    return rows, cov, blk


def part_sim():
    cells = [(si, name, n) for si, name in enumerate(SCENARIOS) for n in NS]
    cells.sort(key=lambda c: -c[2])  # biggest first for load balance
    with Pool(2) as pool:
        res = pool.map(run_cell, cells, chunksize=1)
    rows = [r for a, _, _ in res for r in a]
    cov = [r for _, b, _ in res for r in b]
    blk = [r for _, _, c in res for r in c]
    df = pd.DataFrame(rows)
    df["sqrt_n_rel_rmse"] = np.sqrt(df["n"]) * df["rel_rmse"]
    order = {k: i for i, k in enumerate(SCENARIOS)}
    df = df.sort_values(["scenario", "estimator", "n"], key=lambda s: s.map(order) if s.name == "scenario" else s)
    df.to_csv(OUT / "q2_q4_sim_summary.csv", index=False)
    pd.DataFrame(cov).sort_values(["estimator", "scenario", "n"]).to_csv(OUT / "q4_coverage.csv", index=False)
    pd.DataFrame(blk).sort_values(["scenario", "n", "k"]).to_csv(OUT / "q3_block_length.csv", index=False)
    # Rates: slope of log relative RMSE and log |bias| on log n.
    rr = []
    for (sc, est), g in df.groupby(["scenario", "estimator"]):
        g = g.sort_values("n")
        x = np.log(g["n"].to_numpy(float))
        s_rmse = np.polyfit(x, np.log(g["rel_rmse"].to_numpy()), 1)[0]
        rr.append(dict(scenario=sc, estimator=est, rmse_390=g["rel_rmse"].iloc[0], rmse_1440=g["rel_rmse"].iloc[1],
                       rmse_23400=g["rel_rmse"].iloc[2], bias_390=g["rel_bias"].iloc[0], bias_1440=g["rel_bias"].iloc[1],
                       bias_23400=g["rel_bias"].iloc[2], slope_log_rmse=s_rmse))
    pd.DataFrame(rr).to_csv(OUT / "q2_rates.csv", index=False)


# ----------------------------------------------------------------------------
# Q3: cross-validation for the bandwidth
# ----------------------------------------------------------------------------


def _clip_moments(a):
    """Vectorised E[min(Z^2, a^2)] and E[min(Z^2, a^2)^2] for an array of cut-offs."""
    phi = stats.norm.pdf(a)
    p = 2 * stats.norm.sf(a)
    t1 = (1 - p) - 2 * a * phi
    t2 = 3 * t1 - 2 * a**3 * phi
    return t1 + a**2 * p, t2 + a**4 * p


def run_cv_cell(args):
    si, scen_name, n, reps = args
    rng = np.random.default_rng(5000 + 10 * si + NS.index(n))
    t0 = time.time()
    paths = simulate(SCENARIOS[scen_name], n, reps, rng)
    r = paths.r
    delta = 1.0 / n
    tau = riq.global_threshold(r, 3.0)[:, None]
    rt = np.clip(r, -tau, tau)
    cv = riq.cv_curves(rt, CV_H)
    v = paths.spot_var * delta  # per-step variance
    m1, m2 = _clip_moments(tau / np.sqrt(v))
    tgt2, tgt4 = v * m1, v**2 * m2  # E[r~^2 | sigma], E[r~^4 | sigma] (no jumps)
    hmax = max(CV_H)
    x2, x4 = rt**2, rt**4
    orc4 = np.empty((reps, len(CV_H)))
    orc2 = np.empty((reps, len(CV_H)))
    iq_err = {"RK-k": np.empty((reps, len(CV_H))), "KPI-L3": np.empty((reps, len(CV_H)))}
    s = riq.local_scale(r)
    for j, h in enumerate(CV_H):
        w = R.kernel_weights(h, normalise=False)[:-1]
        w = w / w.sum()
        M = n - h
        P4 = np.zeros((reps, M))
        P2 = np.zeros((reps, M))
        for i in range(h):
            P4 += w[i] * x4[:, i : i + M]
            P2 += w[i] * x2[:, i : i + M]
        sl = slice(hmax - h, None)
        orc4[:, j] = np.mean((P4[:, sl] - tgt4[:, h:][:, sl]) ** 2, axis=-1)
        orc2[:, j] = np.mean((P2[:, sl] - tgt2[:, h:][:, sl]) ** 2, axis=-1)
        iq_err["RK-k"][:, j] = R.kernel_iq(r, h=h, c=3.0, mad_scaled=True, bias_correct=True) / paths.iq - 1
        iq_err["KPI-L3"][:, j] = riq.kernel_plugin(r, 3.0 * s, 3.0, h=h) / paths.iq - 1
    rows = []
    best4 = CV_H[int(np.argmin(orc4.mean(0)))]
    best2 = CV_H[int(np.argmin(orc2.mean(0)))]
    for crit, arr in (("printed eq. 17 (r^4, raw weights)", cv["printed"]), ("corrected (r^4, normalised weights)", cv["corrected"]),
                      ("corrected, squared returns", cv["variance"]), ("oracle spot MSE, r^4 target", orc4), ("oracle spot MSE, r^2 target", orc2)):
        ch = np.asarray(CV_H)[np.argmin(arr, axis=1)]
        row = dict(scenario=scen_name, n=n, criterion=crit, median_h=float(np.median(ch)),
                   share_smallest_h=float(np.mean(ch == CV_H[0])), share_largest_h=float(np.mean(ch == CV_H[-1])),
                   h_minimising_mean_curve=CV_H[int(np.argmin(arr.mean(0)))])
        for h in CV_H:
            row[f"share_h{h}"] = float(np.mean(ch == h))
        # IQ error of the kernel plug-in at the chosen h, per path
        e = iq_err["KPI-L3"][np.arange(reps), np.argmin(arr, axis=1)]
        row["KPI_rel_rmse_at_choice"] = float(np.sqrt(np.mean(e**2)))
        rows.append(row)
    hrows = []
    for j, h in enumerate(CV_H):
        for key, arr in iq_err.items():
            e = arr[:, j]
            hrows.append(dict(scenario=scen_name, n=n, h=h, estimator=key, rel_bias=e.mean(), rel_rmse=np.sqrt(np.mean(e**2)),
                              oracle_mse4_rel=orc4[:, j].mean() / orc4.mean(0).min(), oracle_mse2_rel=orc2[:, j].mean() / orc2.mean(0).min(),
                              cv_printed_mean=cv["printed"][:, j].mean(), cv_corrected_mean=cv["corrected"][:, j].mean()))
    print(f"CV {scen_name:40s} n={n} done in {time.time() - t0:.1f} s (oracle h: r^4 {best4}, r^2 {best2})", flush=True)
    return rows, hrows


def part_cv():
    cells = [(si, name, n, 300) for si, name in enumerate(SCENARIOS) if si in (0, 2, 4) for n in (390, 1440)]
    with Pool(2) as pool:
        res = pool.map(run_cv_cell, cells, chunksize=1)
    pd.DataFrame([r for a, _ in res for r in a]).to_csv(OUT / "q3_cv_choice.csv", index=False)
    pd.DataFrame([r for _, b in res for r in b]).to_csv(OUT / "q3_iq_by_bandwidth.csv", index=False)


# ----------------------------------------------------------------------------
# Q5: real data
# ----------------------------------------------------------------------------


def _daily_returns(df: pd.DataFrame, minutes: int):
    """Dict day -> returns at the given bar length, NaN dropped (log returns summed within bars)."""
    out = {}
    lr = df["logret"]
    for day, g in lr.groupby(lr.index.floor("D")):
        x = g.to_numpy()
        if minutes > 1:
            m = len(x) // minutes * minutes
            x = x[:m].reshape(-1, minutes).sum(axis=1)  # NaN if any minute is missing
        x = x[np.isfinite(x)]
        if len(x) >= 0.9 * (1440 // minutes):
            out[day] = x
    return out


def part_real():
    from avmonitor.data import load_binance

    rows, diag = [], []
    for sym in ("BTCUSDT", "DOGEUSDT"):
        df = load_binance(sym, "2025-07-01", "2026-01-01", download=False)
        for minutes in (1, 5):
            for day, x in _daily_returns(df, minutes).items():
                est, ex = estimate_all(x[None, :])
                row = dict(symbol=sym, bar_minutes=minutes, day=day.date().isoformat(), n=len(x))
                row.update({k: float(v[0]) for k, v in est.items()})
                row["se_PI-LU"] = float(ex["se_LU"][0])
                row["iq8_PI-LU"] = float(ex["iq8_LU"][0])
                rows.append(row)
                if minutes == 1:
                    zz = x / ex["scale"][0]
                    aU = ex["aU"]
                    diag.append(dict(symbol=sym, day=row["day"], share_beyond_3=np.mean(np.abs(zz) > 3),
                                     share_beyond_aU=np.mean(np.abs(zz) > aU), kurtosis_local=float(np.mean(zz**4) / np.mean(zz**2) ** 2),
                                     kurtosis_local_trimmed=float(stats.kurtosis(zz[np.abs(zz) <= aU], fisher=False)),
                                     share_zero=np.mean(x == 0), share_dropped_PI_LU=float(ex["drop_LU"][0])))
    daily = pd.DataFrame(rows)
    daily.to_csv(OUT / "q5_real_daily.csv", index=False)
    dg = pd.DataFrame(diag)
    dg.to_csv(OUT / "q5_real_diagnostics_daily.csv", index=False)
    srows = []
    for (sym, mins), g in daily.groupby(["symbol", "bar_minutes"]):
        for k in EST_ORDER:
            for ref in ("MedRQ", "MinRQ"):
                lr = np.log(g[k] / g[ref])
                srows.append(dict(symbol=sym, bar_minutes=mins, estimator=k, reference=ref, days=len(g),
                                  median_ratio=float(np.exp(np.median(lr))), iqr_log_ratio=float(np.subtract(*np.percentile(lr, [75, 25]))),
                                  median_abs_log_ratio=float(np.median(np.abs(lr))), share_within_25pct=float(np.mean(np.abs(lr) <= np.log(1.25))),
                                  spearman=float(stats.spearmanr(g[k], g[ref])[0])))
    summ = pd.DataFrame(srows)
    summ.to_csv(OUT / "q5_real_summary.csv", index=False)
    # Frequency check: 1-minute estimate over 5-minute estimate, same day.
    one = daily[daily.bar_minutes == 1].set_index(["symbol", "day"])
    five = daily[daily.bar_minutes == 5].set_index(["symbol", "day"])
    common = one.index.intersection(five.index)
    frows = []
    for sym in ("BTCUSDT", "DOGEUSDT"):
        idx = [i for i in common if i[0] == sym]
        for k in EST_ORDER:
            lr = np.log(one.loc[idx, k].to_numpy() / five.loc[idx, k].to_numpy())
            frows.append(dict(symbol=sym, estimator=k, days=len(idx), median_ratio_1min_over_5min=float(np.exp(np.median(lr))),
                              iqr_log_ratio=float(np.subtract(*np.percentile(lr, [75, 25])))))
    pd.DataFrame(frows).to_csv(OUT / "q5_real_frequency.csv", index=False)
    # Cross-frequency agreement: each 1-minute estimator against the 5-minute MedRQ, MinRQ and TPQ of the same day.
    xrows = []
    for sym in ("BTCUSDT", "DOGEUSDT"):
        idx = [i for i in common if i[0] == sym]
        for k in EST_ORDER:
            for ref in ("MedRQ", "MinRQ", "TPQ"):
                lr = np.log(one.loc[idx, k].to_numpy() / five.loc[idx, ref].to_numpy())
                xrows.append(dict(symbol=sym, estimator_1min=k, reference_5min=ref, days=len(idx),
                                  median_ratio=float(np.exp(np.median(lr))), median_abs_log_ratio=float(np.median(np.abs(lr))),
                                  share_within_25pct=float(np.mean(np.abs(lr) <= np.log(1.25)))))
    pd.DataFrame(xrows).to_csv(OUT / "q5_real_cross_frequency.csv", index=False)
    # Sensitivity of the local plug-in on 1-minute data: block length k and cut-off a (local scale, drop).
    srows2 = []
    for sym in ("BTCUSDT", "DOGEUSDT"):
        df = load_binance(sym, "2025-07-01", "2026-01-01", download=False)
        days = _daily_returns(df, 1)
        med1 = one.loc[sym, "MedRQ"]
        med5 = five.loc[sym, "MedRQ"]
        grid = [("k", k, riq.universal_cutoff(1440)) for k in (2, 5, 10, 20, 38, 80, 160)]
        grid += [("a", 38, a) for a in (3.0, riq.universal_cutoff(1440), 5.0, 6.0, 8.0, np.inf)]
        for kind, k, a in grid:
            vals_pi, vals_rq = [], []
            for day, x in days.items():
                s = riq.local_scale(x)
                tau = a * s if np.isfinite(a) else None
                vals_pi.append(riq.block_plugin(x, tau, a, k=k, mode="drop")["iq"])
                vals_rq.append(riq.truncated_rq_local(x, a * s, a, mode="drop") if np.isfinite(a) else R.realized_quarticity(x))
            dkeys = [d.date().isoformat() for d in days]
            pi, rq4 = np.array(vals_pi), np.array(vals_rq)
            m1 = med1.loc[dkeys].to_numpy()
            m5 = med5.reindex(dkeys).to_numpy()  # NaN where the 5-minute day was incomplete
            srows2.append(dict(symbol=sym, vary=kind, k=k, a=a,
                               PI_over_MedRQ_1min=float(np.exp(np.median(np.log(pi / m1)))),
                               PI_over_MedRQ_5min=float(np.exp(np.nanmedian(np.log(pi / m5)))),
                               fourth_power_over_MedRQ_1min=float(np.exp(np.median(np.log(rq4 / m1)))),
                               fourth_power_over_MedRQ_5min=float(np.exp(np.nanmedian(np.log(rq4 / m5))))))
    pd.DataFrame(srows2).to_csv(OUT / "q5_real_sensitivity.csv", index=False)
    drows = []
    for sym, g in dg.groupby("symbol"):
        drows.append(dict(symbol=sym, days=len(g), median_share_beyond_3=g.share_beyond_3.median(), gaussian_share_beyond_3=2 * stats.norm.sf(3),
                          median_share_beyond_aU=g.share_beyond_aU.median(), gaussian_share_beyond_aU=2 * stats.norm.sf(riq.universal_cutoff(1440)),
                          median_kurtosis_local=g.kurtosis_local.median(), median_kurtosis_trimmed=g.kurtosis_local_trimmed.median(),
                          gaussian_kurtosis_trimmed=riq.trunc_even_moments(riq.universal_cutoff(1440), 2)[2] * riq.trunc_even_moments(riq.universal_cutoff(1440), 2)[0] / riq.trunc_even_moments(riq.universal_cutoff(1440), 2)[1] ** 2,
                          median_share_zero=g.share_zero.median(), median_share_dropped_PI_LU=g.share_dropped_PI_LU.median()))
    pd.DataFrame(drows).to_csv(OUT / "q5_real_diagnostics.csv", index=False)
    print(summ[(summ.reference == "MedRQ")].pivot_table(index="estimator", columns=["symbol", "bar_minutes"], values="median_ratio").round(3))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("theory", "all"):
        part_theory()
    if which in ("sim", "all"):
        part_sim()
    if which in ("cv", "all"):
        part_cv()
    if which in ("real", "all"):
        part_real()
    if which in ("figures", "all"):
        import make_figures

        make_figures.main()
