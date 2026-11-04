"""E21. Benchmarks that are not calibrated on Gaussian data only (referee report, M5).

Competitors (``avmonitor.benchmarks_robust``):
  * Classical detectors on the e-detector's own statistic X_t = min(r_t^2, tau^2) (tau = 3 x scaled MAD,
    plug-in m0 from the calibration window, Y = X / m0): Page's two-sided CUSUM on Y with reference values
    from the Gaussian likelihood ratio of an sd change by 1.25 / 0.8 ("small") or 1.5 / 0.67 ("medium"),
    and the e-detector's own Shiryaev-Roberts betting statistic ("SR betting") with a calibrated threshold
    instead of the level 2A given by Ville's inequality. Thresholds (each side at half the target rate) are
    calibrated (a) by simulation on Gaussian data and (b) per path, by a bootstrap of the path's own
    calibration window (resampled X values, m0 recomputed on a resampled window).
  * BOCPD with a Student t3 observation model on a grid of variances (exact run-length recursion with
    runs of at least w = 250 lumped), calibrated as E13 calibrates BOCPD.
  * A CUSUM monitoring procedure with a training sample, of the type of Horvath, Kokoszka and Zhang (2006):
    Q(m, k) = sum of the k monitored Z minus k times the training mean, Z = r^2 (the published setting) or
    Z = X (our robust variant), standardised by a Bartlett long-run sd times sqrt(m), boundary
    c (1 + k/m) (k / (m + k))^gamma with gamma = 0.25 (0 and 0.45 as sensitivity), and the boundary of the
    OLS-CUSUM monitoring of the R package strucchange as a second sensitivity. Critical values: the
    simulated asymptotic one (open-ended monitoring, 5%) and one calibrated on Gaussian data.
Reference: the robust two-sided e-detector (plug-in baseline, each side at 2A, A = 7800), computed with
per-path decrease bets exactly as ``monitor.VolatilityMonitor`` does.

Research questions
  Q21a [SYNTHETIC] False alarms in the first 390-step day (and in two days) under Gaussian, Student t3 and t5,
        jumps, price noise, ticks, slow SV and GARCH 0.99, 4000 paths per law (bootstrap rules: 2000).
  Q21b [SYNTHETIC] Detection in the E4 design (calibration 2000, change at step 500 of 3500 monitored steps),
        sd ratios 1.25, 1.5, 2, 0.8, 0.67, 0.5, laws Gaussian, t3, t5, jumps, price noise, ticks, 1000 paths per cell
        (bootstrap rules: the first 500).
  Q21c [SYNTHETIC] The robust BOCPD in the E13 design, on the very paths of E13 (Q13b null: 1000 per law, seeds
        13310 to 13316; Q13c change: 200 per cell, seeds 13400 to 13429), plus Student t3 (1000 fresh paths).
  Q21d [REAL] Five Binance coins, Jul to Dec 2025 and Jan to Jun 2025, E13's real-data design: rate-matched to the
        e-detector's numbers of up and down alarms per series, validated on the next day (and next four hours)
        against an EWMA forecast, z-test and week-clustered test against the e-detector.

Fixed before any reported run (docs/iterations/E21.md, iteration 1): every setting above, the targets (the
e-detector's own first-day rate on clean Gaussian data, and 5%), nu = 3, the grids, gamma = 0.25 as the primary
HKZ-type boundary, the Bartlett bandwidth floor(m^(1/3)), the seeds and the numbers of paths.

Usage (each step reads what the previous ones wrote):
  python experiments/e21_benchmarks.py calib           # Gaussian thresholds and asymptotic critical values
  python experiments/e21_benchmarks.py null            # Q21a
  python experiments/e21_benchmarks.py change 0|1      # Q21b, two halves; then: python ... change-merge
  python experiments/e21_benchmarks.py e13-calib       # Q21c, BOCPD-t threshold on E13's calibration paths
  python experiments/e21_benchmarks.py e13-null        # Q21c
  python experiments/e21_benchmarks.py e13-change 0|1  # Q21c, two halves; then: python ... e13-change-merge
  python experiments/e21_benchmarks.py real-calib      # Gaussian thresholds for 1440-minute days (as E13)
  python experiments/e21_benchmarks.py real-bocpd W    # BOCPD-t log-odds on the ten coin series, W = 250, 720, 1440
  python experiments/e21_benchmarks.py real HALF       # Q21d, HALF = h2 (Jul to Dec 2025) or h1 (Jan to Jun 2025)
  python experiments/e21_benchmarks.py real-pooled     # pooled real-data table with tests
  python experiments/e21_benchmarks.py figures
Run heavy steps with nice -n 10; every step uses at most two worker processes.
"""

from __future__ import annotations

import os
import sys
import tempfile
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import e13_benchmarks as E13  # noqa: E402

from avmonitor import benchmarks_robust as BR  # noqa: E402
from avmonitor import detectors as D  # noqa: E402
from avmonitor.simulate import change_returns  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
CACHE = Path(os.environ.get("E21_CACHE", str(Path(tempfile.gettempdir()) / "e21_cache")))  # BOCPD arrays, not results

CAL, DAY, A = 2000, 390, 7800.0
LOG2A = float(np.log(2 * A))
NU_T, W_SHORT = 3.0, 250
GRID_SIM = BR.log_variance_grid(0.04, 25.0, 0.1)  # variances 0.04 to 25 (sd 0.2 to 5), 65 points
LAWS_NULL = ("gaussian", "student t3", "student t5", "jumps", "price noise", "ticks", "slow SV", "GARCH 0.99")
LAWS_CHANGE = ("gaussian", "student t3", "student t5", "jumps", "price noise", "ticks")
RATIOS = (1.25, 1.5, 2.0, 0.8, 0.67, 0.5)
N_CALIB, N_NULL, N_NULL_BOOT, N_CHANGE, N_CHANGE_BOOT, N_BOOT = 20000, 4000, 2000, 1000, 500, 1000
NU_CH, N_CH = 500, 3500
NOMINAL = 0.05

K = BR.gaussian_reference_value
CUSUM = {"CUSUM on X, small": (K(1.25), K(0.8)), "CUSUM on X, medium": (K(1.5), K(0.67))}
HKZ = {
    "HKZ-type on r2, gamma 0.25": dict(on="r2", gamma=0.25, boundary="hhks"),
    "HKZ-type on X, gamma 0.25": dict(on="x", gamma=0.25, boundary="hhks"),
    "HKZ-type on r2, gamma 0": dict(on="r2", gamma=0.0, boundary="hhks"),
    "HKZ-type on X, gamma 0": dict(on="x", gamma=0.0, boundary="hhks"),
    "HKZ-type on r2, gamma 0.45": dict(on="r2", gamma=0.45, boundary="hhks"),
    "HKZ-type on X, gamma 0.45": dict(on="x", gamma=0.45, boundary="hhks"),
    "CSW-type boundary on r2": dict(on="r2", gamma=0.0, boundary="csw"),
    "CSW-type boundary on X": dict(on="x", gamma=0.0, boundary="csw"),
}
SIDED = ["SR betting"] + list(CUSUM)  # thresholds per side at half the target rate
BOOT_RULES = {"SR betting": ("sr", {}), "CUSUM on X, small": ("cusum", dict(k_up=K(1.25), k_dn=K(0.8)))}


def gen(name, reps, n, rng, nu=None, ratio=1.0):
    """E13's laws (``e13_benchmarks.gen``) plus unit-variance Student t3."""
    if name == "student t3":
        return change_returns(reps, n, nu=nu, ratio=ratio, rng=rng, heavy_tail_df=3.0).r
    return E13.gen(name, reps, n, rng, nu=nu, ratio=ratio)


def wilson(k, n):
    return E13.wilson(int(k), int(n))


# ----------------------------------------------------------------------------
# Side statistics of every cheap rule on paths r_full = [calibration | monitored]
# ----------------------------------------------------------------------------


def side_stats(r_full, cal_len=CAL):
    """``{rule: (up, down)}`` statistics over the monitored part (larger = stronger evidence), and the inputs."""
    cal, mon = r_full[:, :cal_len], r_full[:, cal_len:]
    ti = BR.truncated_inputs(mon, cal)
    out = {"SR betting": (BR.sr_log_paths(ti.y, D.DEFAULT_UP_GRID, "up"), BR.sr_log_paths(ti.y, BR.down_bets(ti.B), "down"))}
    for nm, (ku, kd) in CUSUM.items():
        out[nm] = BR.page_cusum_paths(ti.y, ku, kd)
    z = {"r2": (cal**2, mon**2), "x": (np.minimum(cal**2, ti.tau2[:, None]), np.minimum(mon**2, ti.tau2[:, None]))}
    for nm, p in HKZ.items():
        Dk = BR.hkz_paths(*z[p["on"]], gamma=p["gamma"], boundary=p["boundary"])
        out[nm] = (Dk, -Dk)
    return out, ti


def threshold_sets(th):
    """``{(rule, calibration label): (h_up, h_dn)}`` from the thresholds table (Gaussian simulation and asymptotic)."""
    return {(r.rule, r.calibration): (r.h_up, r.h_dn) for r in th.itertuples()}


def edetector_alarm(stats_):
    """The e-detector: SR betting statistic at log(2A) on each side (per-path bets, as VolatilityMonitor)."""
    up, dn = stats_["SR betting"]
    return BR.first_alarm(up, dn, LOG2A, LOG2A)


# ----------------------------------------------------------------------------
# Calibration on Gaussian data
# ----------------------------------------------------------------------------


def calib_job(seed):
    r = gen("gaussian", 2000, CAL + DAY, np.random.default_rng(seed))
    st, _ = side_stats(r)
    a, _ = edetector_alarm(st)
    mx = {nm: BR.side_maxima(u, d, DAY) for nm, (u, d) in st.items()}
    return a, mx


def hkz_open_end_job(seed):
    """Finite-m check of the asymptotic HKZ-type critical value: 1000 Gaussian paths, m = 2000, 20,000 monitored steps."""
    rng = np.random.default_rng(seed)
    out = {nm: [] for nm, p in HKZ.items() if p["on"] == "r2"}
    for _ in range(4):
        r = gen("gaussian", 250, CAL + 20000, rng)
        cal, mon = r[:, :CAL], r[:, CAL:]
        for nm in out:
            p = HKZ[nm]
            Dk = BR.hkz_paths(cal**2, mon**2, gamma=p["gamma"], boundary=p["boundary"])
            out[nm].append(np.nanmax(np.abs(Dk), axis=1))
    return {nm: np.concatenate(v) for nm, v in out.items()}


def calib():
    seeds = [21000 + i for i in range(N_CALIB // 2000)]
    with Pool(2) as pool:
        res = pool.map(calib_job, seeds)
    a = np.concatenate([x[0] for x in res])
    fa = float(np.mean((a >= 0) & (a < DAY)))
    print(f"e-detector first-day false alarm rate on {a.size} Gaussian paths: {fa:.5f}", flush=True)
    rows = []
    for nm in res[0][1]:
        mu = np.concatenate([x[1][nm][0] for x in res])
        md = np.concatenate([x[1][nm][1] for x in res])
        for lab, alpha in (("Gaussian simulation, matched", fa), ("Gaussian simulation, 5%", NOMINAL)):
            if nm in SIDED:
                hu, hd, rate = BR.calibrate_two_sided(mu, md, alpha)
            else:  # the published HKZ rule is symmetric: one critical value for |D|
                hu = hd = float(np.quantile(np.maximum(mu, md), 1 - alpha))
                rate = float(np.mean(np.maximum(mu, md) >= hu))
            rows.append(dict(rule=nm, calibration=lab, target=alpha, h_up=hu, h_dn=hd, rate_on_calibration_set=rate,
                             paths=mu.size))
    rows.append(dict(rule="e-detector", calibration="Ville (2A per side)", target=np.nan, h_up=LOG2A, h_dn=LOG2A,
                     rate_on_calibration_set=fa, paths=a.size))
    # asymptotic critical values of the HKZ-type rules (open-ended monitoring, 5%)
    crit = {g: BR.wiener_sup_crit(g, NOMINAL, n_paths=20000, rng=np.random.default_rng(21050 + i))
            for i, g in enumerate((0.0, 0.25, 0.45))}
    cw = BR.csw_crit(NOMINAL)
    for nm, p in HKZ.items():
        c = cw if p["boundary"] == "csw" else crit[p["gamma"]]
        rows.append(dict(rule=nm, calibration="asymptotic, 5% open-ended", target=NOMINAL, h_up=c, h_dn=c,
                         rate_on_calibration_set=np.nan, paths=20000))
    th = pd.DataFrame(rows)
    # finite-m check of the asymptotic values (Gaussian, r^2, 20,000 monitored steps = 10 m), against the limit
    # restricted to t = k / (m + k) <= 10 / 11
    with Pool(2) as pool:
        oe = pool.map(hkz_open_end_job, [21060, 21061])
    chk = []
    for nm in oe[0]:
        p = HKZ[nm]
        c = cw if p["boundary"] == "csw" else crit[p["gamma"]]
        mx = np.concatenate([o[nm] for o in oe])
        chk.append(dict(rule=nm, critical_value=c, paths=mx.size, monitored_steps=20000, m=CAL,
                        rate_within_10m=float(np.mean(mx >= c))))
    for g in (0.0, 0.25, 0.45):
        s = BR.wiener_sup_samples(g, 20000, 10 / 11, np.random.default_rng(21070))
        chk.append(dict(rule=f"limit sup over t <= 10/11, gamma {g}", critical_value=crit[g], paths=20000,
                        monitored_steps=np.nan, m=np.nan, rate_within_10m=float(np.mean(s >= crit[g]))))
    pd.DataFrame(chk).to_csv(OUT / "e21_hkz_open_end_check.csv", index=False)
    th.to_csv(OUT / "e21_thresholds.csv", index=False)
    pd.set_option("display.width", 250)
    print(th.round(4).to_string(), flush=True)
    print(pd.DataFrame(chk).round(4).to_string(), flush=True)


# ----------------------------------------------------------------------------
# Q21a: false alarms
# ----------------------------------------------------------------------------


def _fa_rows(law, kind, method, calib_label, target, a, reps):
    k1 = int(np.sum((a >= 0) & (a < DAY)))
    k2 = int(np.sum((a >= 0) & (a < 2 * DAY)))
    lo, hi = wilson(k1, reps)
    return dict(law=law, method=method, calibration=calib_label, target=target, paths=reps, fa_one_day=k1 / reps,
                fa_one_day_lo=lo, fa_one_day_hi=hi, fa_two_days=k2 / reps, design=kind)


def null_job(args):
    i, law = args
    th = E13._read(OUT / "e21_thresholds.csv")
    ts = threshold_sets(th)
    targets = th.set_index(["rule", "calibration"]).target
    t0 = time.time()
    r = gen(law, N_NULL, CAL + 2 * DAY, np.random.default_rng(21100 + i))
    st, ti = side_stats(r)
    rows = []
    a, _ = edetector_alarm(st)
    rows.append(_fa_rows(law, "E21 null", "e-detector (plug-in, 2A per side)", "Ville", np.nan, a, N_NULL))
    for (nm, lab), (hu, hd) in ts.items():
        if nm == "e-detector":
            continue
        u, d = st[nm]
        a, _ = BR.first_alarm(u, d, hu, hd)
        rows.append(_fa_rows(law, "E21 null", nm, lab, targets[(nm, lab)], a, N_NULL))
    # bootstrap thresholds on the first N_NULL_BOOT paths
    nb = N_NULL_BOOT
    boot = BR.bootstrap_max_stats(ti.ycal[:nb] * ti.m0[:nb, None], ti.tau2[:nb], BOOT_RULES, DAY, N_BOOT,
                                  np.random.default_rng(21150 + i))
    alpha_m = float(th[(th.rule == "SR betting") & (th.calibration == "Gaussian simulation, matched")].target.iloc[0])
    for nm in BOOT_RULES:
        u, d = st[nm]
        for lab, alpha in (("bootstrap, matched", alpha_m), ("bootstrap, 5%", NOMINAL)):
            hu, hd = BR.bootstrap_thresholds(*boot[nm], alpha)
            a, _ = BR.first_alarm(u[:nb], d[:nb], hu, hd)
            rows.append(_fa_rows(law, "E21 null", nm, lab, alpha, a, nb))
    print(f"null {law}: {time.time() - t0:.0f} s", flush=True)
    return rows


def null():
    with Pool(2) as pool:
        rows = sum(pool.map(null_job, list(enumerate(LAWS_NULL))), [])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e21_null.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.pivot_table(index=["method", "calibration"], columns="law", values="fa_one_day").round(4).to_string())


# ----------------------------------------------------------------------------
# Q21b: detection in the E4 design
# ----------------------------------------------------------------------------


def _change_rows(law, ratio, method, lab, target, a, d, reps, nu, design):
    s = D.summarize_alarms(a, nu, 10**9)
    det = a >= nu
    want = 1.0 if ratio > 1 else -1.0
    return dict(design=design, law=law, ratio=ratio, method=method, calibration=lab, target=target, paths=reps, **s,
                correct_direction=float(np.mean(d[det] == want)) if det.any() else np.nan,
                detected_within_390=float(np.mean(det & (a < nu + DAY))))


def change_cell(law, ratio, seed, reps, n_mon, nu, n_boot_paths, design, boot_seed):
    th = E13._read(OUT / "e21_thresholds.csv")
    ts = threshold_sets(th)
    targets = th.set_index(["rule", "calibration"]).target
    r = gen(law, reps, CAL + n_mon, np.random.default_rng(seed), nu=CAL + nu, ratio=ratio)
    st, ti = side_stats(r)
    rows = []
    a, d = edetector_alarm(st)
    rows.append(_change_rows(law, ratio, "e-detector (plug-in, 2A per side)", "Ville", np.nan, a, d, reps, nu, design))
    for (nm, lab), (hu, hd) in ts.items():
        if nm == "e-detector":
            continue
        u, dn = st[nm]
        a, d = BR.first_alarm(u, dn, hu, hd)
        rows.append(_change_rows(law, ratio, nm, lab, targets[(nm, lab)], a, d, reps, nu, design))
    if n_boot_paths:
        nb = n_boot_paths
        boot = BR.bootstrap_max_stats(ti.ycal[:nb] * ti.m0[:nb, None], ti.tau2[:nb], BOOT_RULES, DAY, N_BOOT,
                                      np.random.default_rng(boot_seed))
        alpha_m = float(th[(th.rule == "SR betting") & (th.calibration == "Gaussian simulation, matched")].target.iloc[0])
        for nm in BOOT_RULES:
            u, dn = st[nm]
            for lab, alpha in (("bootstrap, matched", alpha_m), ("bootstrap, 5%", NOMINAL)):
                hu, hd = BR.bootstrap_thresholds(*boot[nm], alpha)
                a, d = BR.first_alarm(u[:nb], dn[:nb], hu, hd)
                rows.append(_change_rows(law, ratio, nm, lab, alpha, a, d, nb, nu, design))
    return rows


def change(part):
    cells = [(law, ratio) for law in LAWS_CHANGE for ratio in RATIOS]
    rows = []
    for i, (law, ratio) in enumerate(cells):
        if i % 2 != part:
            continue
        t0 = time.time()
        rows += change_cell(law, ratio, 21200 + i, N_CHANGE, N_CH, NU_CH, N_CHANGE_BOOT, "E4 design", 21250 + i)
        print(f"change {law} {ratio}: {time.time() - t0:.0f} s", flush=True)
        pd.DataFrame(rows).to_csv(OUT / f"e21_change_part{part}.csv", index=False)


def change_merge():
    df = pd.concat([E13._read(OUT / f"e21_change_part{p}.csv") for p in (0, 1)], ignore_index=True)
    df = df.sort_values(["law", "ratio", "method", "calibration"], kind="stable")
    df.to_csv(OUT / "e21_change.csv", index=False)
    for p in (0, 1):
        (OUT / f"e21_change_part{p}.csv").unlink()
    print(f"{len(df)} rows", flush=True)


# ----------------------------------------------------------------------------
# Q21c: robust BOCPD in the E13 design (the very paths of E13)
# ----------------------------------------------------------------------------


def bocpd_t(r_full, cal_len=CAL, hazard=1 / A, w=W_SHORT):
    """BOCPD-t (nu = 3) on the monitored part, warm-started and with the prior mean from the calibration window."""
    return BR.grid_bocpd_paths(r_full[:, cal_len:], hazard, GRID_SIM, w, nu=NU_T, cal=r_full[:, :cal_len])


def e13_calib_job(seed):
    r = E13.gen("gaussian", 2000, CAL + DAY, np.random.default_rng(seed))
    a, _ = E13.edet_alarm(r[:, CAL:], r[:, :CAL], A)
    L, _ = bocpd_t(r)
    return float(np.mean(a >= 0)), L.max(1)


def e13_calib():
    """As E13.matched_thresholds: the e-detector's rate on E13's 4000 Gaussian calibration paths, BOCPD-t matched to it."""
    with Pool(2) as pool:
        res = pool.map(e13_calib_job, [13300, 13301])
    fa = (res[0][0] + res[1][0]) / 2
    old = E13._read(OUT / "e13b_thresholds.csv").iloc[0]
    assert abs(fa - old.edetector_fa_per_day) < 1e-12, (fa, old.edetector_fa_per_day)
    lm = np.concatenate([x[1] for x in res])
    h = float(np.quantile(lm, 1 - fa))
    df = pd.DataFrame([dict(design="E13 online", day=DAY, A=A, cal_len=CAL, hazard=1 / A, w=W_SHORT, nu=NU_T,
                            grid_points=GRID_SIM.size, paths=lm.size, edetector_fa_per_day=fa,
                            bocpd_t_threshold_log_odds=h, rate_on_calibration_set=float(np.mean(lm >= h)))])
    df.to_csv(OUT / "e21_bocpd_thresholds.csv", index=False)
    print(df.to_string(), flush=True)


def _bocpd_rows(L, S, th):
    rows = np.arange(L.shape[0])
    out = {}
    for nm, h in (("BOCPD-t, matched", th), ("BOCPD-t, p >= 0.5", 0.0)):
        a = E13.C.first_crossing(L, h)
        out[nm] = (a, np.where(a >= 0, S[rows, np.maximum(a, 0)], 0.0))
    return out


def e13_null_job(args):
    law, seed = args
    th = E13._read(OUT / "e21_bocpd_thresholds.csv").bocpd_t_threshold_log_odds.iloc[0]
    t0 = time.time()
    r = gen(law, 1000, CAL + 2 * DAY, np.random.default_rng(seed))
    rows = []
    if law != "student t3":  # reproduction check of the E13 paths
        a, _ = E13.edet_alarm(r[:, CAL:], r[:, :CAL], A)
        rows.append(_fa_rows(law, "E13 null paths", "e-detector, E13 version (common decrease bets)", "Ville", np.nan, a,
                             1000))
    L, S = bocpd_t(r)
    for nm, (a, _) in _bocpd_rows(L, S, th).items():
        rows.append(_fa_rows(law, "E13 null paths", nm, "Gaussian simulation, matched (E13)", np.nan, a, 1000))
    print(f"e13 null {law}: {time.time() - t0:.0f} s", flush=True)
    return rows


def e13_null():
    jobs = [(law, 13310 + i) for i, law in enumerate(E13.NULL_SCEN)] + [("student t3", 21400)]
    with Pool(2) as pool:
        rows = sum(pool.map(e13_null_job, jobs), [])
    df = pd.DataFrame(rows)
    old = E13._read(OUT / "e13b_online_null.csv")
    ref = old[old.method == "e-detector, plug-in"].set_index("scenario").fa_one_day
    chk = df[df.method.str.startswith("e-detector")].set_index("law").fa_one_day
    assert np.allclose(chk.to_numpy(), ref.loc[chk.index].to_numpy(), atol=1e-12), "E13 null paths not reproduced"
    df.to_csv(OUT / "e21_e13_null.csv", index=False)
    print(df.pivot_table(index="method", columns="law", values="fa_one_day").round(4).to_string(), flush=True)


def e13_change(part):
    th = E13._read(OUT / "e21_bocpd_thresholds.csv").bocpd_t_threshold_log_odds.iloc[0]
    cells = [(nm, ratio) for nm in E13.CHANGE_SCEN for ratio in E13.RATIOS]
    rows = []
    for i, (law, ratio) in enumerate(cells):
        if i % 2 != part:
            continue
        t0 = time.time()
        r = E13.gen(law, 200, CAL + E13.N_ON, np.random.default_rng(13400 + i), nu=CAL + E13.NU_ON, ratio=ratio)
        a, d = E13.edet_alarm(r[:, CAL:], r[:, :CAL], A)
        rows.append(_change_rows(law, ratio, "e-detector, E13 version (common decrease bets)", "Ville", np.nan, a, d,
                                 200, E13.NU_ON, "E13 change paths"))
        L, S = bocpd_t(r)
        for nm, (a, d) in _bocpd_rows(L, S, th).items():
            rows.append(_change_rows(law, ratio, nm, "Gaussian simulation, matched (E13)", np.nan, a, d, 200, E13.NU_ON,
                                     "E13 change paths"))
        # the cheap rules on the same paths, at their own Gaussian-matched thresholds
        rows += [x for x in change_cell_on(r, law, ratio)]
        print(f"e13 change {law} {ratio}: {time.time() - t0:.0f} s", flush=True)
        pd.DataFrame(rows).to_csv(OUT / f"e21_e13_change_part{part}.csv", index=False)


def change_cell_on(r, law, ratio):
    th = E13._read(OUT / "e21_thresholds.csv")
    th = th[th.calibration.isin(["Gaussian simulation, matched", "asymptotic, 5% open-ended"])]
    st, _ = side_stats(r)
    rows = []
    a, d = edetector_alarm(st)
    rows.append(_change_rows(law, ratio, "e-detector (plug-in, 2A per side)", "Ville", np.nan, a, d, r.shape[0], E13.NU_ON,
                             "E13 change paths"))
    for row in th.itertuples():
        if row.rule == "e-detector" or (row.rule.startswith(("HKZ", "CSW")) and "0.25" not in row.rule
                                         and "CSW" not in row.rule):
            continue
        u, dn = st[row.rule]
        a, d = BR.first_alarm(u, dn, row.h_up, row.h_dn)
        rows.append(_change_rows(law, ratio, row.rule, row.calibration, row.target, a, d, r.shape[0], E13.NU_ON,
                                 "E13 change paths"))
    return rows


def e13_change_merge():
    df = pd.concat([E13._read(OUT / f"e21_e13_change_part{p}.csv") for p in (0, 1)], ignore_index=True)
    old = pd.concat([E13._read(OUT / f"e13c_online_change_part{p}.csv") for p in (0, 1)], ignore_index=True)
    ref = old[old.method == "e-detector, plug-in"].set_index(["scenario", "ratio"])
    chk = df[df.method == "e-detector, E13 version (common decrease bets)"].set_index(["law", "ratio"])
    for col in ("median_delay", "detected", "false_alarm_before_change"):
        assert np.allclose(chk[col].to_numpy(), ref.loc[chk.index, col].to_numpy(), equal_nan=True), col
    df.to_csv(OUT / "e21_e13_change.csv", index=False)
    for p in (0, 1):
        (OUT / f"e21_e13_change_part{p}.csv").unlink()
    print(f"{len(df)} rows; E13 e-detector reproduced", flush=True)


# ----------------------------------------------------------------------------
# Q21d: real data (E13's design)
# ----------------------------------------------------------------------------

RDAY, RA, RCAL = 1440, 28800.0, 1440
REAL_RULES = {  # restart rules on real data: name -> (kind, params, bisection range)
    "CUSUM on X, small": ("cusum", dict(k_up=K(1.25), k_dn=K(0.8)), (0.0, 5000.0)),
    "CUSUM on X, medium": ("cusum", dict(k_up=K(1.5), k_dn=K(0.67)), (0.0, 5000.0)),
    "SR betting": ("sr", dict(band=0.0), (0.0, 60.0)),
    "HKZ-type on r2, gamma 0.25": ("hkz", dict(on="r2", gamma=0.25, boundary="hhks"), (0.0, 500.0)),
    "HKZ-type on X, gamma 0.25": ("hkz", dict(on="x", gamma=0.25, boundary="hhks"), (0.0, 500.0)),
}
BOCPD_WS = (250, 720, 1440)
HALVES = {"h2": ("five Binance coins, Jul to Dec 2025", "crypto"), "h1": ("five Binance coins, Jan to Jun 2025",
                                                                         "crypto_h1_2025")}


def real_calib_job(seed):
    r = E13.gen("gaussian", 1000, RCAL + RDAY, np.random.default_rng(seed))
    a, _ = E13.edet_alarm(r[:, RCAL:], r[:, :RCAL], RA)
    st, _ = side_stats(r, RCAL)
    mx = {nm: BR.side_maxima(u, d, RDAY) for nm, (u, d) in st.items() if nm in REAL_RULES}
    L, _ = BR.grid_bocpd_paths(r[:, RCAL:], 1 / RA, GRID_REAL, W_SHORT, nu=NU_T, cal=r[:, :RCAL])  # the real-data grid
    return float(np.mean(a >= 0)), mx, L.max(1)


def real_calib():
    """Gaussian-matched thresholds for 1440-minute days on E13's calibration paths (seeds 13500, 13501)."""
    with Pool(2) as pool:
        res = pool.map(real_calib_job, [13500, 13501])
    fa = (res[0][0] + res[1][0]) / 2
    old = E13._read(OUT / "e13d_real_thresholds.csv").set_index("day").loc[1440]
    assert abs(fa - old.edetector_fa_per_day) < 1e-12, (fa, old.edetector_fa_per_day)
    rows = []
    for nm in res[0][1]:
        mu = np.concatenate([x[1][nm][0] for x in res])
        md = np.concatenate([x[1][nm][1] for x in res])
        if nm in SIDED:
            hu, hd, rate = BR.calibrate_two_sided(mu, md, fa)
        else:
            hu = hd = float(np.quantile(np.maximum(mu, md), 1 - fa))
            rate = float(np.mean(np.maximum(mu, md) >= hu))
        rows.append(dict(rule=nm, h_up=hu, h_dn=hd, target=fa, rate_on_calibration_set=rate, paths=mu.size))
    lm = np.concatenate([x[2] for x in res])
    h = float(np.quantile(lm, 1 - fa))
    rows.append(dict(rule=f"BOCPD-t, w={W_SHORT}", h_up=h, h_dn=h, target=fa, rate_on_calibration_set=float(np.mean(lm >= h)),
                     paths=lm.size))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e21_real_thresholds.csv", index=False)
    print(df.to_string(), flush=True)


def real_series(half):
    """The deseasonalised series of E13 for one half-year: {coin: (idx, u, i0)} (``h2``: Jul to Dec 2025, ``h1``: Jan to
    Jun 2025)."""
    from e06b_five_coins_monitor import COINS, prepare

    if half == "h2":
        return {c: prepare(c, "Jul-Dec 2025") for c in COINS}
    return E13.prepare_h1_2025()


def real_scale(u, i0):
    """Median trailing robust variance over the 30 days before i0 (predictable): the centre of the variance grid."""
    pv = E13.trailing_robust_var(u, RDAY)[max(0, i0 - 30 * RDAY) : i0]
    return float(np.nanmedian(pv))


GRID_REAL = BR.log_variance_grid(np.exp(-6.0), np.exp(4.5), 0.15)  # sd from 0.05 to 9.5 times the centre, 71 points


def real_bocpd(w):
    """BOCPD-t with short-run window w on the ten series (E13's real_bocpd conventions); arrays cached in CACHE.

    Each series is divided by the square root of its grid centre ``real_scale`` (a constant known before the
    monitored period), so that one relative grid serves every coin. The filter is scale-equivariant (likelihood,
    prior and predictable prior mean all scale with the series), so this equals running it on the raw series with
    the grid multiplied by the centre.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    out = {}
    t0 = time.time()
    for half in ("h2", "h1"):
        ser = real_series(half)
        coins = list(ser)
        lens = {ser[c][1].size - ser[c][2] for c in coins}
        assert len(lens) == 1, lens
        v0 = np.array([real_scale(ser[c][1], ser[c][2]) for c in coins])
        um = np.stack([ser[c][1][ser[c][2]:] for c in coins]) / np.sqrt(v0)[:, None]
        prior = np.stack([E13.trailing_robust_var(ser[c][1], RDAY)[ser[c][2]:] for c in coins]) / v0[:, None]
        L, S = BR.grid_bocpd_paths(um, 1 / RA, GRID_REAL, w, nu=NU_T, cal=um[:, :RCAL], prior_var_path=prior)
        for j, c in enumerate(coins):
            out[f"L_{half}_{c}"], out[f"S_{half}_{c}"] = L[j], S[j]
        print(f"BOCPD-t w={w} {half}: {time.time() - t0:.0f} s", flush=True)
    np.savez_compressed(CACHE / f"bocpd_t_w{w}.npz", **out)


def real_job(args):
    half, coin = args
    idx, u, i0 = real_series(half)[coin]
    market, tag = HALVES[half]
    um = u[i0:]
    n = um.size
    start = dead = RCAL
    t0 = time.time()
    events = {}
    from avmonitor.monitor import VolatilityMonitor

    old = E13._events_from_csv(OUT / f"e13d_real_{tag}_events.csv", idx, coin)
    for b in (0.0, 0.25):
        ev = [(e.index, e.direction) for e in VolatilityMonitor(threshold=RA, c=3.0, cal_len=RCAL, null_band=b).run(um)]
        assert [(t + i0, d) for t, d in ev] == old[f"e-detector, band {b}"], f"{coin}: e-detector not reproduced"
        events[f"e-detector, band {b}"] = ev
    gth = E13._read(OUT / "e21_real_thresholds.csv").set_index("rule")
    matched = []
    for nm, (kind, prm, lohi) in REAL_RULES.items():
        if nm in gth.index:
            events[f"{nm}, Gaussian-matched"] = BR.restart_events(um, RCAL, kind, gth.loc[nm, "h_up"], gth.loc[nm, "h_dn"],
                                                                  prm)
        for b in (0.0, 0.25):
            ref = events[f"e-detector, band {b}"]
            nu_ = sum(1 for e in ref if e[1] == "up")
            nd_ = len(ref) - nu_
            hu, hd, ev = BR.match_restart_counts(um, RCAL, kind, prm, nu_, nd_, lohi, iters=22, rounds=2, refine=6,
                                                 grid_points=15 if kind == "sr" else 41)
            events[f"{nm}, matched to band {b}"] = ev
            matched.append(dict(market=market, series=coin, band=b, method=nm, h_up=hu, h_down=hd, target_up=nu_,
                                target_down=nd_, got_up=sum(1 for e in ev if e[1] == "up"),
                                got_down=sum(1 for e in ev if e[1] == "down")))
        print(f"{half} {coin} {nm}: {time.time() - t0:.0f} s", flush=True)
    hb = float(gth.loc[f"BOCPD-t, w={W_SHORT}", "h_up"])
    for w in BOCPD_WS:
        z = np.load(CACHE / f"bocpd_t_w{w}.npz")
        L, S = z[f"L_{half}_{coin}"], z[f"S_{half}_{coin}"]
        Lu, Ld = np.where(S > 0, L, -np.inf), np.where(S < 0, L, -np.inf)
        if w == W_SHORT:
            events[f"BOCPD-t, w={w}, Gaussian-matched"] = E13.C.scan_alarms(Lu, Ld, hb, hb, start, dead)
        for b in (0.0, 0.25):
            ref = events[f"e-detector, band {b}"]
            nu_ = sum(1 for e in ref if e[1] == "up")
            nd_ = len(ref) - nu_
            hu, hd = E13.C.match_alarm_counts(Lu, Ld, nu_, nd_, start, dead, lo_hi=(-50.0, 1000.0), iters=50)
            ev = E13.C.scan_alarms(Lu, Ld, hu, hd, start, dead)
            events[f"BOCPD-t, w={w}, matched to band {b}"] = ev
            matched.append(dict(market=market, series=coin, band=b, method=f"BOCPD-t, w={w}", h_up=hu, h_down=hd,
                                target_up=nu_, target_down=nd_, got_up=sum(1 for e in ev if e[1] == "up"),
                                got_down=sum(1 for e in ev if e[1] == "down")))
    V, V4 = E13.Validation(u, RDAY), E13.Validation(u, E13.HORIZON_SHORT)
    alarms, summ = [], []
    for m, ev in events.items():
        evu = [(t + i0, d) for t, d in ev]
        h1, h4 = E13.score_events(V, evu), E13.score_events(V4, evu)
        alarms += [dict(market=market, series=coin, method=m, time=idx[t], direction=d, hit_next_day=x, hit_next_4h=y)
                   for (t, d), x, y in zip(evu, h1, h4)]
        rec = dict(market=market, series=coin, method=m, days=n / RDAY, alarms=len(ev), alarms_per_day=len(ev) / (n / RDAY))
        rec.update(V.score(evu))
        summ.append(rec)
    print(f"{half} {coin} done: {time.time() - t0:.0f} s", flush=True)
    return alarms, summ, matched


def real(half):
    from e06b_five_coins_monitor import COINS

    with Pool(2) as pool:
        res = pool.map(real_job, [(half, c) for c in COINS])
    pd.DataFrame(sum([r[0] for r in res], [])).to_csv(OUT / f"e21_real_{half}_alarms.csv", index=False)
    pd.DataFrame(sum([r[1] for r in res], [])).to_csv(OUT / f"e21_real_{half}_summary.csv", index=False)
    pd.DataFrame(sum([r[2] for r in res], [])).to_csv(OUT / f"e21_real_{half}_matched_thresholds.csv", index=False)


def real_pooled():
    """Pooled shares per market, method and direction, with the z-test and the week-clustered test against the e-detector."""
    al = pd.concat([E13._read(OUT / f"e21_real_{h}_alarms.csv") for h in ("h2", "h1")], ignore_index=True)
    t = pd.to_datetime(al.time.str[:-6])
    al["cl_week"] = t.dt.strftime("%G-W%V")
    rows = []
    for (mk, m), g in al.groupby(["market", "method"], sort=False):
        if "matched to band" in m:
            band = "0.0" if m.endswith("band 0.0") else "0.25"
        elif m.startswith("e-detector"):
            band = m.split("band ")[1]
        else:
            band = "0.0"  # Gaussian-matched rows are compared with band 0
        ref = al[(al.market == mk) & (al.method == f"e-detector, band {band}")]
        rec = dict(market=mk, method=m, band=float(band), alarms=len(g))
        for dn in ("up", "down"):
            for hz, col in (("day", "hit_next_day"), ("4h", "hit_next_4h")):
                a = g[(g.direction == dn) & g[col].notna()]
                b = ref[(ref.direction == dn) & ref[col].notna()]
                k, nn = int(a[col].sum()), len(a)
                rec[f"n_{dn}_{hz}"] = nn
                rec[f"{dn}_rate_{hz}"] = k / nn if nn else np.nan
                lo, hi = wilson(k, nn)
                rec[f"{dn}_{hz}_lo"], rec[f"{dn}_{hz}_hi"] = lo, hi
                if m.startswith("e-detector") or nn == 0 or len(b) == 0:
                    continue
                rec[f"{dn}_{hz}_p_ztest"] = E13.two_prop_p(int(b[col].sum()), len(b), k, nn)
                rec[f"{dn}_{hz}_p_week_clustered"] = E13.cluster_z(b[col].to_numpy(), b.cl_week.to_numpy(),
                                                                   a[col].to_numpy(), a.cl_week.to_numpy())
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "e21_real_pooled.csv", index=False)
    pd.set_option("display.width", 300)
    cols = ["market", "method", "alarms", "n_up_day", "up_rate_day", "up_day_p_ztest", "up_day_p_week_clustered",
            "n_down_day", "down_rate_day", "down_day_p_ztest", "down_day_p_week_clustered", "up_rate_4h", "down_rate_4h"]
    print(out[[c for c in cols if c in out]].round(4).to_string(index=False), flush=True)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------


def figures():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(exist_ok=True)
    nd = E13._read(OUT / "e21_null.csv")
    show = [("e-detector (plug-in, 2A per side)", "Ville", "e-detector (Ville)", "#1b4f72", "o"),
            ("SR betting", "Gaussian simulation, matched", "SR betting, Gaussian-calibrated", "#5dade2", "s"),
            ("SR betting", "bootstrap, matched", "SR betting, bootstrap", "#2e86c1", "D"),
            ("CUSUM on X, small", "Gaussian simulation, matched", "CUSUM on X, Gaussian-calibrated", "#e67e22", "s"),
            ("CUSUM on X, small", "bootstrap, matched", "CUSUM on X, bootstrap", "#a04000", "D"),
            ("HKZ-type on r2, gamma 0.25", "Gaussian simulation, matched", "HKZ-type on r², Gaussian-calibrated", "#7d3c98", "^"),
            ("HKZ-type on X, gamma 0.25", "Gaussian simulation, matched", "HKZ-type on X, Gaussian-calibrated", "#c39bd3", "v")]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    laws = list(LAWS_NULL)
    for j, (m, lab, name, col, mk) in enumerate(show):
        g = nd[(nd.method == m) & (nd.calibration == lab)].set_index("law").loc[laws]
        x = np.arange(len(laws)) + (j - len(show) / 2) * 0.09
        ax.errorbar(x, 100 * g.fa_one_day, yerr=[100 * (g.fa_one_day - g.fa_one_day_lo), 100 * (g.fa_one_day_hi - g.fa_one_day)],
                    fmt=mk, color=col, ms=5, lw=1, label=name)
    ax.axhline(5, color="0.4", ls="--", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(0.2, 100)
    tk = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]
    ax.set_yticks(tk, [f"{t:g}%" for t in tk])
    ax.minorticks_off()
    ax.set_xticks(np.arange(len(laws)), laws, rotation=20)
    ax.set_ylabel("false alarm in the first day")
    ax.set_title("E21 [SYNTHETIC]: first-day false alarms; thresholds for the e-detector's Gaussian rate (simulation or bootstrap)",
                 fontsize=9)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "e21_null_false_alarms.png", dpi=150)
    plt.close(fig)

    cd = E13._read(OUT / "e21_change.csv")
    ref = cd[cd.method == "e-detector (plug-in, 2A per side)"].set_index(["law", "ratio"]).median_delay
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, ups in zip(axes, (True, False)):
        rs = [r for r in RATIOS if (r > 1) == ups]
        for m, lab, name, col, mk in show[1:]:
            g = cd[(cd.method == m) & (cd.calibration == lab) & (cd.law == "gaussian") & cd.ratio.isin(rs)]
            ax.plot(g.ratio, g.median_delay.to_numpy() / ref.loc[[("gaussian", r) for r in g.ratio]].to_numpy(), mk + "-",
                    color=col, label=name, ms=5)
        ax.axhline(1, color="0.4", ls="--", lw=0.8)
        ax.set_xlabel("sd ratio")
        ax.set_title("increases" if ups else "decreases")
    axes[0].set_ylabel("median delay / e-detector's (Gaussian)")
    axes[1].legend(fontsize=7, frameon=False)
    fig.suptitle("E21 [SYNTHETIC]: delay relative to the e-detector, E4 design, Gaussian returns", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "e21_delay_ratio.png", dpi=150)
    plt.close(fig)

    if (OUT / "e21_real_pooled.csv").exists():
        pdf = E13._read(OUT / "e21_real_pooled.csv")
        pdf = pdf[pdf.method.str.contains("matched to band 0.25") | (pdf.method == "e-detector, band 0.25")]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
        for ax, mk in zip(axes, ("five Binance coins, Jul to Dec 2025", "five Binance coins, Jan to Jun 2025")):
            g = pdf[pdf.market == mk]
            e = g[g.method == "e-detector, band 0.25"].iloc[0]
            g = g[g.method != "e-detector, band 0.25"]
            y = np.arange(len(g))
            ax.errorbar(100 * g.up_rate_day, y + 0.15, xerr=[100 * (g.up_rate_day - g.up_day_lo), 100 * (g.up_day_hi - g.up_rate_day)],
                         fmt="o", color="#c0392b", ms=4, lw=1, label="up alarms")
            ax.errorbar(100 * g.down_rate_day, y - 0.15, xerr=[100 * (g.down_rate_day - g.down_day_lo),
                                                                 100 * (g.down_day_hi - g.down_rate_day)],
                        fmt="s", color="#2471a3", ms=4, lw=1, label="down alarms")
            ax.axvline(100 * e.up_rate_day, color="#c0392b", ls="--", lw=0.8)
            ax.axvline(100 * e.down_rate_day, color="#2471a3", ls=":", lw=0.8)
            ax.set_yticks(y, [m.replace(", matched to band 0.25", "") for m in g.method], fontsize=7)
            ax.set_title(mk.replace("five Binance coins, ", ""), fontsize=9)
            ax.set_xlabel("next day beyond the EWMA forecast (%)")
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle("E21 [REAL]: rate-matched to the e-detector (band 0.25); lines: e-detector", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG / "e21_real_validation.png", dpi=150)
        plt.close(fig)


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "calib"
    arg = sys.argv[2] if len(sys.argv) > 2 else None
    {"calib": calib, "null": null, "change": lambda: change(int(arg)), "change-merge": change_merge,
     "e13-calib": e13_calib, "e13-null": e13_null, "e13-change": lambda: e13_change(int(arg)),
     "e13-change-merge": e13_change_merge, "real-calib": real_calib, "real-bocpd": lambda: real_bocpd(int(arg)),
     "real": lambda: real(arg), "real-pooled": real_pooled, "figures": figures}[cmd]()


if __name__ == "__main__":
    main()
