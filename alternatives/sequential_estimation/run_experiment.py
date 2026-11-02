"""Experiments for the sequential (stopping-rule) estimation prototype (ALT6).

Parts (run any subset: ``python alternatives/sequential_estimation/run_experiment.py main garch``)

  limit  [SYNTHETIC, closed form] Predicted E[N] / n* of the confidence sequences:
         the limit 2 log(2/alpha) / z^2, the empirical Bernstein prediction with tuned
         and optimised bets, and the best constant-bet betting sequence for the known
         law, for several d, alpha, truncation constants c and laws.
  main   [SYNTHETIC] Relative accuracy (half-width d theta). Gaussian, Student t5,
         Student t3 and Gaussian with Bernoulli jumps; d in {0.4, 0.2, 0.1, 0.05};
         alpha in {0.05, 0.10}; pilot m = 30. Coverage at N, E[N] / n*, distribution
         of N and regret E[N] - n* for every procedure on the same paths.
  abs    [SYNTHETIC] Absolute accuracy (half-width d, sigma^2 = 1), alpha = 0.05.
  pilot  [SYNTHETIC] Pilot size m in {10, 30, 100}.
  garch  [SYNTHETIC] GARCH(1, 1) clustering, persistence 0.95 and 0.99: coverage of the
         stationary truncated variance and of the average conditional truncated variance
         over the sample; batch-means (long-run variance) Chow-Robbins.
  minrisk [SYNTHETIC] Robbins (1959) minimum risk point estimation, raw and truncated:
         risk efficiency E[loss] / (2 n0) and regret E[loss] - 2 n0.
  verify [SYNTHETIC] E[N] / n* of the confidence sequences (Gaussian, relative accuracy, pilot
         100) with and without the running intersection, against the per-path first-order
         prediction, the limit 2 log(2/alpha) / z^2 and the Brownian limit with intersection.
  map    [SYNTHETIC, closed form and simulation] Mapping the truncated variance back to
         sigma^2 under a Gaussian shape: bias when the shape is wrong, and coverage of
         sigma^2 by the mapped betting interval.
  real   [REAL] Binance BTCUSDT one-minute returns, 1 Jul to 31 Dec 2025, divided by the
         predictable intraday profile: minutes needed per procedure for each start hour.
  boot   [REAL, resampled] Stationary bootstrap of the 30 most homogeneous days: coverage
         of the day's own truncated variance (i.i.d. and 30-minute blocks).
  figs   Figures from the CSV files.

Seeds are fixed per cell. At most two worker processes.
"""

from __future__ import annotations

import math
import sys
import time
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))

import seqest as S  # noqa: E402

RES = HERE / "results"
RES.mkdir(exist_ok=True)
WORKERS = 2


def seed_of(key: str) -> int:
    return zlib.crc32(key.encode()) % (2**31)


def law_of(name: str):
    return {
        "gauss": S.Gaussian(),
        "t5": S.StudentT(5),
        "t3": S.StudentT(3),
        "jumps": S.Jumps(0.01, 5.0),
        "garch95": S.Garch(0.10, 0.85),
        "garch99": S.Garch(0.05, 0.94),
    }[name]


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def _cell(args):
    """One simulation cell: every procedure on the same paths; returns summary rows."""
    part, lawname, d, alpha, m, relative, procs, reps, n_max, extra = args
    law = law_of(lawname)
    key = f"{part}|{lawname}|{d}|{alpha}|{m}|{relative}" + ("" if extra.get("intersect", True) else "|nointersect")
    st = S.Setup(d=d, alpha=alpha, m=m, relative=relative, procs=tuple(procs), n_max=n_max, **extra)
    t0 = time.time()
    gt = None
    if isinstance(law, S.Garch):
        gt = garch_targets(lawname)
    out = S.run(law, st, reps, seed_of(key), garch_targets=gt)
    rows = []
    for p in procs:
        res = out[p]
        sm = S.summarise(res, d, relative, alpha)
        row = dict(part=part, law=lawname, d=d, alpha=alpha, m=m, relative=relative, proc=p, reps=reps,
                   intersect=st.intersect, **sm)
        if res.cond_target is not None:
            row["cover_cond_avg"] = float(S.covered(res, d, relative, target=res.cond_target).mean())
            row["cond_over_uncond"] = float(np.mean(res.cond_target / res.target))
        if res.cond_target_w is not None:
            # the theorem's targets: lambda-weighted averages of the conditional means, one per side
            row["cs_theorem_cover"] = float(((res.lo <= res.cond_target_w[1]) & (res.cond_target_w[0] <= res.hi)).mean())
        if res.loss is not None:
            row["loss_mean"] = float(np.mean(res.loss))
        rows.append(row)
    log(f"{key} done in {time.time() - t0:.0f}s")
    return rows


_GT_CACHE = {}


def garch_targets(lawname):
    """Stationary truncated variance and long-run variance of the truncated statistic on a tau grid (cached)."""
    if lawname in _GT_CACHE:
        return _GT_CACHE[lawname]
    path = RES / f"_garch_table_{lawname}.csv"
    if path.exists():
        tab = pd.read_csv(path)
    else:
        law = law_of(lawname)
        rng = np.random.default_rng(seed_of("garch-table" + lawname))
        reps, n = 8, 400_000
        st = law.start(rng, reps)
        r, _, _ = law.draw(rng, st, np.arange(reps), n)
        taus = np.linspace(0.5, 8.0, 76)
        B = 5000
        m2, lrv = [], []
        for t in taus:
            x = np.minimum(r * r, t * t)
            m2.append(x.mean())
            bm = x[:, : n // B * B].reshape(reps, -1, B).mean(axis=2)
            lrv.append(B * bm.var(ddof=1))
        x = r * r
        bm = x[:, : n // B * B].reshape(reps, -1, B).mean(axis=2)
        lrv_raw = B * bm.var(ddof=1)
        log(f"{lawname}: long-run variance of r^2 by batch means {lrv_raw:.2f}, closed form "
            f"{law.nu_raw() * law.lrv_factor_raw():.2f}")
        tab = pd.DataFrame(dict(tau=taus, m2=m2, lrv=lrv))
        tab.to_csv(path, index=False)
    f = lambda tau: (np.interp(tau, tab.tau, tab.m2), np.interp(tau, tab.tau, tab.lrv))  # noqa: E731
    _GT_CACHE[lawname] = f
    return f


def run_cells(cells, out_name):
    with Pool(WORKERS) as pool:
        rows = [r for rs in pool.imap_unordered(_cell, cells) for r in rs]
    df = pd.DataFrame(rows).sort_values(["law", "d", "alpha", "m", "proc", "intersect"])
    df.to_csv(RES / out_name, index=False)
    log(f"wrote {out_name} ({len(df)} rows)")
    return df


PROCS_ALL = ["CR-raw", "CR-trunc", "Stein-raw", "Stein-trunc", "Oracle-raw", "Oracle-trunc",
             "CS-wsr", "CS-tuned", "CS-opt", "CS-bet"]

# ----------------------------------------------------------------------------- parts


def part_limit():
    rows = []
    for alpha in (0.05, 0.10):
        for lawname in ("gauss", "t5", "t3", "jumps"):
            law = law_of(lawname)
            for c in ((2.0, 3.0, 4.0) if lawname == "gauss" else (3.0,)):
                tau = law.tau_limit(c)
                m2, m4 = law.trunc_moments(np.array([tau]))
                mu_z, var_z = m2[0] / tau**2, (m4[0] - m2[0] ** 2) / tau**4
                zl, pl = law.z_law(tau, k=2000)
                for d in (0.4, 0.2, 0.1, 0.05, 0.02, 0.01):
                    ns = S.zq(alpha) ** 2 * var_z / (d * mu_z) ** 2
                    row = dict(alpha=alpha, law=lawname, c=c, d=d, tau=tau, theta_tau=m2[0],
                               cv2=var_z / mu_z**2, rho=d * mu_z / var_z, nstar=ns, limit=S.cs_ratio_limit(alpha))
                    row["eb_tuned"] = S.cs_n_predict(mu_z, var_z, d * mu_z, alpha, bets="tuned", relative=True,
                                                     d_rel=d) / ns
                    row["eb_opt"] = S.cs_n_predict(mu_z, var_z, d * mu_z, alpha, bets="opt", relative=True,
                                                   d_rel=d) / ns
                    row["bet_oracle"] = S.betting_oracle_n(zl, pl, d, alpha, relative=True) / ns
                    rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "a_closed_form.csv", index=False)
    log("wrote a_closed_form.csv")
    # d -> 0 limit with the bounds intersected over time: a Brownian first-passage functional
    rows = []
    for alpha in (0.05, 0.10):
        for kappa in (0.8, 1.0, 1.2, 1.5):
            res = S.brownian_limit(alpha, kappa, reps=40_000, seed=seed_of(f"brown{alpha}{kappa}"))
            no = S.brownian_limit(alpha, kappa, intersect=False)["ratio"]
            rows.append(dict(alpha=alpha, kappa=kappa, ratio_intersected=res["ratio"], se=res["se"],
                             coverage_at_stop=res["coverage"], ratio_not_intersected=no))
    pd.DataFrame(rows).to_csv(RES / "a_brownian_limit.csv", index=False)
    log("wrote a_brownian_limit.csv")


def part_main():
    cells = []
    for lawname in ("gauss", "t5", "t3", "jumps"):
        for d in (0.4, 0.2, 0.1, 0.05):
            for alpha in (0.05, 0.10):
                procs = [p for p in PROCS_ALL if not (lawname == "t3" and p == "Oracle-raw")]
                cells.append(("main", lawname, d, alpha, 30, True, procs, 2000, 600_000, {}))
    # biggest cells first so that the two workers finish together
    cells.sort(key=lambda c: (c[2], c[1] != "t3"))
    run_cells(cells, "b_main.csv")


def part_abs():
    cells = []
    for lawname in ("gauss", "t5", "t3", "jumps"):
        for d in (0.2, 0.1, 0.05):
            procs = [p for p in PROCS_ALL if not (lawname == "t3" and p == "Oracle-raw") and p != "CS-wsr"]
            cells.append(("abs", lawname, d, 0.05, 30, False, procs, 2000, 600_000, {}))
    cells.sort(key=lambda c: (c[2], c[1] != "t3"))
    run_cells(cells, "c_absolute.csv")


def part_pilot():
    cells = []
    for lawname in ("gauss", "t3", "jumps"):
        for d in (0.2, 0.1):
            for m in (10, 30, 100):
                procs = ["CR-raw", "CR-trunc", "Stein-trunc", "Oracle-trunc", "CS-opt", "CS-bet"]
                cells.append(("pilot", lawname, d, 0.05, m, True, procs, 2000, 600_000, {}))
    run_cells(cells, "d_pilot.csv")


def part_garch():
    for name in ("garch95", "garch99"):
        garch_targets(name)
    cells = []
    for lawname, batch in (("garch95", 200), ("garch99", 1000)):
        for d in (0.2, 0.1):
            procs = ["CR-raw", "CR-trunc", "CR-BM-trunc", "Oracle-trunc", "CS-opt", "CS-bet"]
            cells.append(("garch", lawname, d, 0.05, 30, True, procs, 1000, 200_000, {"batch": batch}))
            # without the running intersection the empirical Bernstein sequence covers the bet-weighted
            # average of the conditional truncated variances (docs/math.md Theorem 6.2); with it, not exactly
            cells.append(("garch", lawname, d, 0.05, 30, True, ["CS-opt", "CS-bet"], 1000, 200_000,
                          {"batch": batch, "intersect": False}))
    run_cells(cells, "e_garch.csv")


def _minrisk_cell(args):
    lawname, n0, kind, reps = args
    law = law_of(lawname)
    m = 10
    key = f"minrisk|{lawname}|{n0}|{kind}"
    if kind == "raw":
        nu = law.nu_raw()
        A = n0**2 / nu
    else:
        tau = law.tau_limit(3.0)
        m2, m4 = law.trunc_moments(np.array([tau]))
        A = n0**2 / (m4[0] - m2[0] ** 2)
    st = S.Setup(d=1.0, m=m, n_max=400_000, minrisk_A=A, procs=(f"MinRisk-{kind}",))
    out = S.run(law, st, reps, seed_of(key))
    res = out[f"MinRisk-{kind}"]
    if kind == "raw":
        n0_path = np.full(reps, math.sqrt(A * law.nu_raw()))
    else:
        tau = out["_meta"]["tau"]
        m2, m4 = law.trunc_moments(tau)
        n0_path = np.sqrt(A * (m4 - m2**2))
    risk_opt = 2 * n0_path
    loss = res.loss
    return dict(law=lawname, kind=kind, n0=n0, reps=reps, A=A, EN=float(res.N.mean()),
                EN_over_n0=float(np.mean(res.N / n0_path)), risk=float(loss.mean()),
                risk_eff=float(np.mean(loss) / np.mean(risk_opt)), regret=float(np.mean(loss - risk_opt)),
                regret_se=float(np.std(loss - risk_opt) / math.sqrt(reps)), capped=float(res.capped.mean()))


def part_minrisk():
    cells = [(lawname, n0, kind, 4000) for lawname in ("gauss", "t5", "jumps") for n0 in (50, 200, 800)
             for kind in ("raw", "trunc")]
    with Pool(WORKERS) as pool:
        rows = pool.map(_minrisk_cell, cells)
    pd.DataFrame(rows).to_csv(RES / "f_minrisk.csv", index=False)
    log("wrote f_minrisk.csv")


def part_map():
    rows = []
    for lawname in ("gauss", "t5", "t3", "jumps"):
        law = law_of(lawname)
        for c in (2.0, 3.0, 4.0):
            tau = law.tau_limit(c)
            m2, _ = law.trunc_moments(np.array([tau]))
            sig2_g = S.sigma2_from_theta_tau(m2, np.array([tau]))[0]
            diff = 1.0 if lawname != "jumps" else law.sigma**2  # diffusive variance
            total = law.theta_raw()
            rows.append(dict(kind="closed_form", law=lawname, c=c, tau=tau, theta_tau=m2[0],
                             theta_over_total_var=m2[0] / total, gauss_sigma2_over_diffusive=sig2_g / diff,
                             gauss_sigma2_over_total=sig2_g / total,
                             elasticity_gauss=S.elasticity_gaussian(tau / math.sqrt(diff))))
    # simulation: map the CS-bet interval endpoints to sigma^2 under the Gaussian shape (monotone map)
    for lawname in ("gauss", "t5", "t3", "jumps"):
        law = law_of(lawname)
        d, alpha = 0.2, 0.05
        st = S.Setup(d=d, alpha=alpha, m=30, relative=True, procs=("CS-bet", "CR-trunc"))
        out = S.run(law, st, 2000, seed_of("map" + lawname))
        tau = out["_meta"]["tau"]
        diff = 1.0 if lawname != "jumps" else law.sigma**2
        for p in ("CS-bet", "CR-trunc"):
            res = out[p]
            lo = S.sigma2_from_theta_tau(np.maximum(res.lo, 1e-12), tau)
            hi = S.sigma2_from_theta_tau(np.minimum(res.hi, tau**2 * (1 - 1e-12)), tau)
            est = S.sigma2_from_theta_tau(res.est, tau)
            rows.append(dict(kind="simulation", law=lawname, c=3.0, proc=p, d=d, alpha=alpha,
                             interval_cover_sigma2=float(((lo <= diff) & (diff <= hi)).mean()),
                             within_d_sigma2=float((np.abs(est / diff - 1) <= d).mean()),
                             median_ratio_sigma2=float(np.median(est / diff)),
                             interval_cover_theta_tau=float(((res.lo <= res.target) & (res.target <= res.hi)).mean())))
    pd.DataFrame(rows).to_csv(RES / "g_mapping.csv", index=False)
    log("wrote g_mapping.csv")


def _verify_cell(args):
    d, intersect, reps = args
    law = law_of("gauss")
    alpha, m = 0.05, 100
    procs = ("CR-trunc", "CS-tuned", "CS-opt", "CS-bet")
    key = f"verify|{d}|{intersect}"
    st = S.Setup(d=d, alpha=alpha, m=m, relative=True, procs=procs, n_max=400_000, intersect=intersect)
    t0 = time.time()
    out = S.run(law, st, reps, seed_of(key))
    tau = out["_meta"]["tau"]
    m2, m4 = law.trunc_moments(tau)
    mz, vz = m2 / tau**2, (m4 - m2**2) / tau**4
    lim_b = S.brownian_limit(alpha, 1.0, reps=20_000, seed=1)["ratio"]
    rows = []
    for p in procs:
        res = out[p]
        sm = S.summarise(res, d, True, alpha)
        row = dict(d=d, intersect=intersect, proc=p, reps=reps, m=m, **sm)
        if p in ("CS-tuned", "CS-opt"):
            # first-order prediction per path (this path's tau), N = m + predicted sequence length
            pred = np.array([S.cs_n_predict(a, b, d * a, alpha, bets=p.split("-")[1], relative=True, d_rel=d)
                             for a, b in zip(mz, vz)])
            row["pred_EN_over_nstar"] = float(np.mean((m + pred) / res.nstar))
        row["limit_not_intersected"] = S.cs_ratio_limit(alpha)
        row["limit_intersected_kappa1"] = lim_b
        rows.append(row)
    log(f"{key} done in {time.time() - t0:.0f}s")
    return rows


def part_verify():
    """[SYNTHETIC] E[N] / n* of the confidence sequences against the closed forms, with and without intersection."""
    cells = [(d, inter, 1000) for d in (0.025, 0.05, 0.1) for inter in (False, True)]
    with Pool(WORKERS) as pool:
        rows = [r for rs in pool.imap_unordered(_verify_cell, cells) for r in rs]
    df = pd.DataFrame(rows).sort_values(["intersect", "d", "proc"])
    df.to_csv(RES / "j_verify_limits.csv", index=False)
    log("wrote j_verify_limits.csv")


# ----------------------------------------------------------------------------- real data


def load_btc_deseasonalised():
    from avmonitor.data import intraday_scale, load_binance

    df = load_binance("BTCUSDT", "2025-06-01", "2026-01-01", download=False)
    r = df["logret"]
    prof = intraday_scale(r, lookback_days=20, smooth=15, start="past")
    u = (r / prof).loc["2025-07-01":"2025-12-31"]
    checks = df.attrs.get("checks", {})
    return u, checks


def real_windows(u: pd.Series, m: int, cap: int):
    """One window per (day, hour): the next m + cap valid deseasonalised returns from the hour's first minute."""
    valid = u.dropna()
    vals = valid.to_numpy()
    times = valid.index
    starts = pd.date_range(u.index[0].floor("D"), u.index[-1].floor("D"), freq="1h", tz="UTC")
    L = m + cap
    pos = np.searchsorted(times.values, starts.values)
    keep = pos + L <= vals.size
    pos, starts = pos[keep], starts[keep]
    idx = pos[:, None] + np.arange(L)[None, :]
    R = vals[idx]
    elapsed = (times.values[idx[:, -1]] - starts.values) / np.timedelta64(1, "m")
    return R, starts, elapsed, times.values[idx]


def part_real():
    u, checks = load_btc_deseasonalised()
    log(f"BTCUSDT checks: {checks}; deseasonalised returns {u.notna().sum()} of {u.size}")
    m, cap = 30, 1440
    R, starts, elapsed, tidx = real_windows(u, m, cap)
    log(f"{R.shape[0]} windows (start hours), median span of {m + cap} returns: {np.median(elapsed):.0f} minutes")
    rows, summ = [], []
    for alpha in (0.05, 0.10):
        for d in (0.5, 0.4, 0.3, 0.2):
            procs = ("CR-raw", "CR-trunc", "Stein-trunc", "CS-opt", "CS-bet")
            st = S.Setup(d=d, alpha=alpha, m=m, relative=True, procs=procs, n_max=m + cap)
            out = S.run(S.Windows(R), st, R.shape[0], 0)
            for p in procs:
                res = out[p]
                N = res.N
                # minutes on the clock from the hour's start to the stopping return (gaps make it longer than N)
                minutes = (tidx[np.arange(N.size), np.minimum(N, m + cap) - 1] - starts.values) / np.timedelta64(1, "m") + 1
                df = pd.DataFrame(dict(hour=starts.hour, N=N, minutes=minutes, capped=res.capped,
                                       est=res.est, lo=res.lo, hi=res.hi))
                for h, g in df.groupby("hour"):
                    rows.append(dict(alpha=alpha, d=d, proc=p, hour=h, windows=len(g), N_median=g.N.median(),
                                     N_mean=g.N.mean(), N_q90=g.N.quantile(0.9),
                                     minutes_median=g.minutes.median(),
                                     stop_within_60=float((g.minutes <= 60).mean()),
                                     stop_within_240=float((g.minutes <= 240).mean()),
                                     no_stop_within_cap=float(g.capped.mean())))
                summ.append(dict(alpha=alpha, d=d, proc=p, windows=len(df), N_median=df.N.median(),
                                 N_mean=df.N.mean(), N_q90=df.N.quantile(0.9),
                                 stop_within_60=float((df.minutes <= 60).mean()),
                                 no_stop_within_cap=float(df.capped.mean()),
                                 hour_with_largest_median=int(df.groupby("hour").N.median().idxmax()),
                                 largest_hour_median=float(df.groupby("hour").N.median().max()),
                                 smallest_hour_median=float(df.groupby("hour").N.median().min())))
            log(f"real alpha={alpha} d={d} done")
    pd.DataFrame(rows).to_csv(RES / "h_real_hourly.csv", index=False)
    pd.DataFrame(summ).to_csv(RES / "h_real_summary.csv", index=False)
    log("wrote h_real_hourly.csv and h_real_summary.csv")


def quiet_days(u: pd.Series, k: int = 30):
    """Days with the most homogeneous volatility: smallest dispersion of hourly robust scales (and complete)."""
    g = u.groupby(u.index.floor("D"))
    rows = []
    for day, x in g:
        if x.notna().sum() < 1430:
            continue
        hs = x.groupby(x.index.hour).apply(lambda v: np.nanmedian(np.abs(v)))
        rows.append(dict(day=day, disp=float(np.std(np.log(hs))), scale=float(np.nanmedian(np.abs(x)))))
    df = pd.DataFrame(rows).sort_values("disp")
    return df.head(k)


def _boot_cell(args):
    day, vals, block, d, alpha, reps = args
    law = S.Bootstrap(vals, block=block, name=f"boot{block}")
    procs = ("CR-raw", "CR-trunc", "Stein-trunc", "CS-opt", "CS-bet")
    st = S.Setup(d=d, alpha=alpha, m=30, relative=True, procs=procs, n_max=50_000)
    out = S.run(law, st, reps, seed_of(f"boot|{day}|{block}|{d}|{alpha}"))
    rows = []
    for p in procs:
        sm = S.summarise(out[p], d, True, alpha)
        rows.append(dict(day=str(day)[:10], block=block, d=d, alpha=alpha, proc=p, reps=reps, **sm))
    return rows


def part_boot():
    u, _ = load_btc_deseasonalised()
    q = quiet_days(u, 30)
    log(f"quiet days: dispersion {q.disp.min():.3f} to {q.disp.max():.3f}")
    cells = []
    for _, row in q.iterrows():
        vals = u.loc[row.day: row.day + pd.Timedelta("1D") - pd.Timedelta("1min")].to_numpy()
        for block in (1, 30):
            cells.append((row.day, vals, block, 0.2, 0.05, 200))
    with Pool(WORKERS) as pool:
        rows = [r for rs in pool.imap_unordered(_boot_cell, cells) for r in rs]
    df = pd.DataFrame(rows)
    df.to_csv(RES / "i_bootstrap_days.csv", index=False)
    agg = df.groupby(["block", "proc"]).agg(days=("day", "nunique"), coverage=("coverage", "mean"),
                                             cov_min=("coverage", "min"), miss_under=("miss_under", "mean"),
                                             EN=("EN", "mean"), EN_over_nstar=("EN_over_nstar", "mean"),
                                             capped=("capped", "mean")).reset_index()
    agg.to_csv(RES / "i_bootstrap_summary.csv", index=False)
    log("wrote i_bootstrap_days.csv and i_bootstrap_summary.csv")


# ----------------------------------------------------------------------------- figures


def part_figs():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"CR-raw": "#c0392b", "CR-trunc": "#2471a3", "Stein-trunc": "#7d3c98", "CS-wsr": "#bdbdbd",
              "CS-tuned": "#9e9e9e", "CS-opt": "#616161", "CS-bet": "#117a65", "Oracle-trunc": "#000000"}
    b = pd.read_csv(RES / "b_main.csv")
    a = pd.read_csv(RES / "a_closed_form.csv")
    # Figure 1: E[N] / n* of the confidence sequences against d, with the closed forms (Gaussian, alpha 0.05)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, alpha in zip(axes, (0.05, 0.10)):
        sub = b[(b.law == "gauss") & (b.alpha == alpha)]
        for p in ("CS-wsr", "CS-tuned", "CS-opt", "CS-bet", "CR-trunc"):
            s = sub[sub.proc == p].sort_values("d")
            ax.plot(s.d, s.EN_over_nstar, "o-", color=colors[p], label=f"{p} (simulated)")
        ca = a[(a.law == "gauss") & (a.c == 3.0) & (a.alpha == alpha)].sort_values("d")
        ax.plot(ca.d, ca.eb_opt, "--", color=colors["CS-opt"], label="CS-opt (predicted)")
        ax.plot(ca.d, ca.bet_oracle, "--", color=colors["CS-bet"], label="best betting bound (predicted)")
        ax.axhline(S.cs_ratio_limit(alpha), color="k", lw=0.8, ls=":", label="limit 2 log(2/alpha) / z^2")
        bl = pd.read_csv(RES / "a_brownian_limit.csv")
        bl = bl[(bl.alpha == alpha) & (bl.kappa == 1.0)].ratio_intersected.iloc[0]
        ax.axhline(bl, color="k", lw=0.8, ls="--", label="limit with intersection (Brownian)")
        ax.axhline(1.0, color="k", lw=0.5)
        ax.set_xscale("log")
        ax.set_xlabel("relative half-width d")
        ax.set_title(f"Gaussian returns, alpha = {alpha} [SYNTHETIC]")
    axes[0].set_ylabel("E[N] / n*")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(RES / "fig1_price_of_anytime_validity.png", dpi=130)
    plt.close(fig)
    # Figure 2: coverage at N against d, by law
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.6), sharey=True)
    for ax, lawname in zip(axes, ("gauss", "t5", "t3", "jumps")):
        sub = b[(b.law == lawname) & (b.alpha == 0.05)]
        for p in ("CR-raw", "CR-trunc", "Stein-trunc", "CS-bet"):
            s = sub[sub.proc == p].sort_values("d")
            ax.plot(s.d, s.coverage, "o-", color=colors[p], label=p)
        ax.axhline(0.95, color="k", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_ylim(0.4, 1.01)
        ax.set_title(f"{lawname} [SYNTHETIC]")
        ax.set_xlabel("relative half-width d")
    axes[0].set_ylabel("coverage at N (alpha = 0.05)")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(RES / "fig2_coverage_by_law.png", dpi=130)
    plt.close(fig)
    # Figure 3: real data, median minutes by start hour
    h = RES / "h_real_hourly.csv"
    if h.exists():
        r = pd.read_csv(h)
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
        sub = r[(r.alpha == 0.05) & (r.d == 0.3)]
        for p in ("CR-raw", "CR-trunc", "CS-opt", "CS-bet"):
            s = sub[sub.proc == p].sort_values("hour")
            axes[0].plot(s.hour, s.N_median, "o-", color=colors[p], label=p)
        axes[0].set_xlabel("start hour (UTC)")
        axes[0].set_ylabel("median returns until stop")
        axes[0].set_title("BTCUSDT, +/- 30%, alpha 0.05 [REAL]")
        axes[0].legend(fontsize=7)
        sub = r[(r.alpha == 0.05) & (r.d == 0.5)]
        for p in ("CR-raw", "CR-trunc", "CS-opt", "CS-bet"):
            s = sub[sub.proc == p].sort_values("hour")
            axes[1].plot(s.hour, 1 - s.stop_within_60, "o-", color=colors[p], label=p)
        axes[1].set_xlabel("start hour (UTC)")
        axes[1].set_ylabel("share not stopped within 60 minutes")
        axes[1].set_title("BTCUSDT, +/- 50%, alpha 0.05 [REAL]")
        fig.tight_layout()
        fig.savefig(RES / "fig3_real_minutes_by_hour.png", dpi=130)
        plt.close(fig)
    log("figures written")


PARTS = dict(limit=part_limit, main=part_main, abs=part_abs, pilot=part_pilot, garch=part_garch,
             minrisk=part_minrisk, map=part_map, verify=part_verify, real=part_real, boot=part_boot, figs=part_figs)

if __name__ == "__main__":
    names = sys.argv[1:] or [p for p in PARTS if p != "figs"] + ["figs"]
    for nm in names:
        t0 = time.time()
        log(f"part {nm}")
        PARTS[nm]()
        log(f"part {nm} finished in {time.time() - t0:.0f}s")
