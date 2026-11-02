"""E13. Benchmarks against strong published competitors (tracker B67).

Competitors (``avmonitor.competitors``):
  * BJV: the adjacent-window volatility jump test in the spirit of Bibinger, Jirak and Vetter (2017,
    Annals of Statistics 45(4)). Local realised variances over two adjacent windows of k returns,
    standardised difference, maximum over the sample, Gumbel-type critical value; jump-robust variant
    with truncation; change point at the maximum. Offline (a posteriori) by design.
  * BOCPD: Bayesian online change-point detection (Adams and MacKay 2007, arXiv:0710.3742) with a
    zero-mean Gaussian likelihood, inverse gamma prior on the variance, constant hazard and a capped
    run-length posterior. Alarm when the posterior mass of runs shorter than w steps exceeds a threshold.
Our method: the robust two-sided betting e-detector (``detectors.betting_edetector``, each side at 2A),
with the default plug-in baseline and with the confidence-sequence baseline (``calibrate(baseline="cs")``).

Research questions
  Q13a. Offline power at matched size. On a fixed window of 2000 steps with at most one change at step
        1000, what is the power of BJV (non-overlapping blocks and overlapping windows, plain and
        truncated), BOCPD (maximum posterior short-run mass over the window) and the e-detector (maximum
        of its statistic over the window) when every critical value gives exactly 5% size on clean
        Gaussian data? What size do these critical values give under jumps, Student t5, price noise,
        ticks, slow stochastic volatility and GARCH? Do the asymptotic (Gumbel) critical values of BJV
        hold their nominal 5%? How accurate is the BJV change-point estimate?
  Q13b. Online false alarms. Used online (BJV as a rolling test, BOCPD with its posterior-mass rule),
        with thresholds matched to the e-detector's false alarm rate per 390-step day on clean Gaussian
        data, what are the competitors' false alarm rates per day under the same contaminations?
  Q13c. Online detection. Detection probability, median delay and direction for sd ratios 1.25, 1.5, 2,
        0.8, 0.67 and 0.5 under Gaussian, jumps, Student t5, price noise and ticks.
  Q13d. Real data. Five Binance coins (Jul to Dec 2025) and Nifty 50 (2024 to 2025, regular sessions):
        alarms per day with the Gaussian-matched thresholds, and a validation matched on alarm counts
        (reviewer point M3): after each alarm, is the robust volatility of the next day above (up alarm)
        or below (down alarm) an EWMA forecast and the previous day, compared with the e-detector at the
        same numbers of up and down alarms?
  Q13e. Cost. Time per step of each method in a one-series streaming loop.

Fixed in advance (not tuned on the reported data):
  * e-detector: c = 3, plug-in baseline from the calibration window, A = 7800 for 390-step days
    (A = 28,800 for 1440-minute crypto days, A = 7500 for 375-minute Nifty sessions), band 0, and on real
    data the frozen E6 / E7 settings (cal 1440 and band 0 / 0.25 for crypto, cal 750 and band 0 / 0.25
    for Nifty).
  * BOCPD: alpha0 = 2, prior mean of the variance = calibration mean square (simulation) or the robust
    variance of the previous day (real data, predictable), hazard 1 / A (one change per 20 days a
    priori), short-run window w = 250 and run-length cap 512, chosen on a Gaussian pilot with seed 999
    (docs/iterations/E13.md, iteration 1) to give the competitor its best detection.
  * BJV online: k = 120, truncation at 4 times the larger bipower scale of the two windows (pilot as
    above). BJV offline: k = 50 (about sqrt(n)) and k = 125 (the value used in the paper's
    simulations, as reported by later work), truncation constant 4.
  * Thresholds of the competitors come from simulated Gaussian data only. On real data the only
    adjustment is the explicit matching of alarm counts to the e-detector, as in E6b.

All simulated results are [SYNTHETIC]; Binance and Nifty results are [REAL]. Seeds are fixed below.

Usage
  python experiments/e13_benchmarks.py offline-null    # Q13a: critical values and size
  python experiments/e13_benchmarks.py offline-change  # Q13a: power and change-point error
  python experiments/e13_benchmarks.py online-null  # Q13b (and the per-day matched thresholds)
  python experiments/e13_benchmarks.py online-change [part]   # Q13c, part 0 or 1 (two halves)
  python experiments/e13_benchmarks.py real-calib   # Gaussian thresholds for 1440- and 375-step days
  python experiments/e13_benchmarks.py real-crypto  # Q13d, five Binance coins
  python experiments/e13_benchmarks.py real-crypto-h1-2025  # Q13d, the same coins, 1 Jan to 30 Jun 2025
  python experiments/e13_benchmarks.py real-nifty   # Q13d, Nifty 50
  python experiments/e13_benchmarks.py cost         # Q13e
  python experiments/e13_benchmarks.py figures      # tables and figures from the saved CSV files
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from avmonitor import competitors as C
from avmonitor import detectors as D
from avmonitor.robust import MAD_NORMAL_SCALE
from avmonitor.simulate import change_returns

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
sys.path.insert(0, str(Path(__file__).resolve().parent))

CAL = 2000
DAY = 390
A = 7800.0  # 390 / 0.05
H_SIM = 1.0 / A
W_SHORT, R_MAX, ALPHA0 = 250, 512, 2.0
K_ONLINE, TRUNC = 120, 4.0
RATIOS = (1.25, 1.5, 2.0, 0.8, 0.67, 0.5)
NULL_SCEN = ("gaussian", "jumps", "student t5", "price noise", "ticks", "slow SV", "GARCH 0.99")
CHANGE_SCEN = ("gaussian", "jumps", "student t5", "price noise", "ticks")
KW = {
    "gaussian": dict(),
    "jumps": dict(p_jump=0.005, jump_scale=8.0),
    "student t5": dict(heavy_tail_df=5.0),
    "price noise": dict(noise_sd=0.5),
    "ticks": dict(),
    "slow SV": dict(sv_xi=0.15, sv_halflife=390.0),
}
TICK = 0.5  # in units of the pre-change sd: about 20% zero returns (as in E9)


def garch(reps, n, rng, alpha=0.05, beta=0.94, burn=3000):
    """GARCH(1,1) returns with unit unconditional variance (persistence 0.99, as in E9 Q9c)."""
    omega = 1 - alpha - beta
    s2 = np.ones(reps)
    r_prev = rng.standard_normal(reps)
    out = np.empty((reps, n))
    for t in range(burn + n):
        s2 = omega + alpha * r_prev**2 + beta * s2
        r_prev = np.sqrt(s2) * rng.standard_normal(reps)
        if t >= burn:
            out[:, t - burn] = r_prev
    return out


def gen(name, reps, n, rng, nu=None, ratio=1.0):
    """Paths of length n with an sd change by ``ratio`` at index ``nu`` (None: no change)."""
    if name == "GARCH 0.99":
        r = garch(reps, n, rng)
        if nu is not None:
            r[:, nu:] *= ratio
        return r
    r = change_returns(reps, n, nu=nu, ratio=ratio, rng=rng, **KW[name]).r
    if name == "ticks":
        r = np.round(r / TICK) * TICK
    return r


# ----------------------------------------------------------------------------
# The e-detector, as in E4 (each side at 2A; per-path plug-in or confidence-sequence baseline)
# ----------------------------------------------------------------------------


def edet_inputs(r, cal, baseline="plugin", c=3.0):
    """Normalised truncated squares for the up and down sides, and the upper bound for the down side."""
    reps = r.shape[0]
    bs = [D.calibrate(cal[i], c=c, baseline=baseline, delta=0.01) for i in range(reps)]
    tau2 = np.array([b.tau**2 for b in bs])[:, None]
    mu = np.array([b.m0_up for b in bs])[:, None]
    md = np.array([b.m0_down for b in bs])[:, None]
    x = np.minimum(r**2, tau2)
    return x / mu, x / md, float(np.max(tau2 / md))


def edet_alarm(r, cal, threshold, baseline="plugin"):
    """First alarm (index into r, -1 if none) and its direction (+1 up, -1 down) of the two-sided e-detector."""
    yu, yd, ub = edet_inputs(r, cal, baseline)
    up = D.betting_edetector(yu, 1.0, "up", threshold=2 * threshold)
    dn = D.betting_edetector(yd, 1.0, "down", upper=ub, threshold=2 * threshold)
    a = D.two_sided_alarm(up, dn)
    direction = np.where((up.alarm >= 0) & (up.alarm == a), 1.0, -1.0)
    return a, direction


def edet_max_stat(r, cal, baseline="plugin"):
    """max over t of max(M_up, M_down) / 2 for the window: the two-sided test rejects at level A iff this >= A.

    Same recursion as ``detectors.betting_edetector`` (kind "SR", default grids); checked against it in
    ``check_edet_max_stat``.
    """
    yu, yd, ub = edet_inputs(r, cal, baseline)
    gu = np.asarray(D.DEFAULT_UP_GRID)
    gd = np.asarray(D.DEFAULT_DOWN_FRACTIONS) / (ub - 1.0)
    reps, n = yu.shape
    Su, Sd = np.zeros((reps, gu.size)), np.zeros((reps, gd.size))
    best = np.zeros(reps)
    with np.errstate(over="ignore", invalid="ignore"):
        for t in range(n):
            Su = (1.0 + gu[None, :] * (yu[:, t, None] - 1.0)) * (Su + 1.0)
            Sd = (1.0 + gd[None, :] * (1.0 - yd[:, t, None])) * (Sd + 1.0)
            best = np.maximum(best, np.maximum(Su.mean(1), Sd.mean(1)) / 2.0)
    return best


def check_edet_max_stat(seed=5):
    rng = np.random.default_rng(seed)
    cal = rng.standard_normal((200, 500))
    r = change_returns(200, 800, nu=400, ratio=1.4, rng=rng).r
    m = edet_max_stat(r, cal)
    a, _ = edet_alarm(r, cal, 300.0)
    assert np.array_equal(m >= 300.0, a >= 0), "edet_max_stat disagrees with detectors.betting_edetector"


# ----------------------------------------------------------------------------
# Q13a: offline power at matched size
# ----------------------------------------------------------------------------

N_OFF, NU_OFF, CAL_OFF = 2000, 1000, 500
OFF_METHODS = [
    ("BJV blocks, k=50, truncated", dict(variant="blocks", k=50, trunc_c=TRUNC)),
    ("BJV overlap, k=50, truncated", dict(variant="overlap", k=50, trunc_c=TRUNC)),
    ("BJV overlap, k=50, plain", dict(variant="overlap", k=50, trunc_c=None)),
    ("BJV overlap, k=125, truncated", dict(variant="overlap", k=125, trunc_c=TRUNC)),
    # Our modification, not in the paper: studentise with the sample variance of the squared returns
    ("BJV overlap, k=50, truncated, studentised", dict(variant="overlap", k=50, trunc_c=TRUNC, norm="empirical")),
]


def offline_stats(r):
    """Test statistic of every offline method on windows r (reps, N_OFF); larger means more evidence."""
    out, cp, nominal = {}, {}, {}
    for name, kw in OFF_METHODS:
        res = C.bjv_test(r, kw["k"], 0.05, kw["variant"], kw["trunc_c"], kw.get("norm", "gaussian"))
        out[name], cp[name], nominal[name] = res.stat, res.cp_index, res.crit
    cal, mon = r[:, :CAL_OFF], r[:, CAL_OFF:]
    _, _, L = C.bocpd_paths(mon, H_SIM, cal=cal, alpha0=ALPHA0, r_max=R_MAX, short=W_SHORT, return_log_odds=True)
    out["BOCPD, w=250"] = L.max(1)  # log-odds of a run shorter than w (no saturation at p = 1)
    out["e-detector, plug-in"] = edet_max_stat(mon, cal)
    nominal["e-detector, plug-in"] = (N_OFF - CAL_OFF) / 0.05  # P(alarm within m steps) <= m / A with the true baseline
    out["e-detector, CS baseline"] = edet_max_stat(mon, cal, baseline="cs")
    nominal["e-detector, CS baseline"] = (N_OFF - CAL_OFF) / (0.05 - 0.01)  # delta + m / A = 5% with delta = 0.01
    return out, cp, nominal


def offline_job(args):
    kind, name, ratio, reps, seed = args
    rng = np.random.default_rng(seed)
    nu = None if kind == "null" else NU_OFF
    r = gen(name, reps, N_OFF, rng, nu=nu, ratio=ratio)
    t0 = time.time()
    s, cp, nominal = offline_stats(r)
    print(f"offline {kind} {name} {ratio}: {time.time() - t0:.0f} s", flush=True)
    return kind, name, ratio, s, cp, nominal


def offline_null():
    """Gaussian calibration set (critical values for 5% size) and size under the null scenarios."""
    check_edet_max_stat()
    jobs = [("calib", "gaussian", 1.0, 1000, 13100), ("calib", "gaussian", 1.0, 1000, 13101)]
    jobs += [("null", nm, 1.0, 1000, 13110 + i) for i, nm in enumerate(NULL_SCEN)]
    with Pool(2) as pool:
        res = pool.map(offline_job, jobs)
    cal = [x for x in res if x[0] == "calib"]
    s_cal = {m: np.concatenate([c[3][m] for c in cal]) for m in cal[0][3]}
    nominal = cal[0][5]
    crit = pd.DataFrame([dict(method=m, crit_matched=float(np.quantile(v, 0.95)), crit_nominal=nominal.get(m, np.nan),
                              calibration_paths=v.size) for m, v in s_cal.items()])
    crit.to_csv(OUT / "e13a_offline_crit.csv", index=False)
    rows = [dict(kind="calibration set", scenario="gaussian", ratio=1.0, method=m, paths=v.size) for m, v in s_cal.items()]
    stats_ = [s_cal[m] for m in s_cal]
    for kind, name, ratio, s, cp, _ in res:
        if kind == "null":
            for m, v in s.items():
                rows.append(dict(kind=kind, scenario=name, ratio=ratio, method=m, paths=v.size))
                stats_.append(v)
    _offline_table(rows, stats_, crit, "e13a_offline_size.csv")


def offline_change():
    crit = pd.read_csv(OUT / "e13a_offline_crit.csv")
    jobs = [("change", nm, ratio, 300, 13200 + 10 * i + j) for i, nm in enumerate(("gaussian", "jumps", "student t5"))
            for j, ratio in enumerate(RATIOS)]
    with Pool(2) as pool:
        res = pool.map(offline_job, jobs)
    rows, stats_ = [], []
    for kind, name, ratio, s, cp, _ in res:
        for m, v in s.items():
            row = dict(kind=kind, scenario=name, ratio=ratio, method=m, paths=v.size)
            if m in cp:
                hit = v > crit.set_index("method").loc[m, "crit_matched"]
                err = np.abs(cp[m][hit] - NU_OFF)
                row["cp_median_abs_error"] = float(np.median(err)) if err.size else np.nan
                row["cp_share_within_25"] = float(np.mean(err <= 25)) if err.size else np.nan
            rows.append(row)
            stats_.append(v)
    _offline_table(rows, stats_, crit, "e13a_offline_power.csv")


def _offline_table(rows, stats_, crit, fname):
    c = crit.set_index("method")
    for row, v in zip(rows, stats_):
        cm, cn = c.loc[row["method"], "crit_matched"], c.loc[row["method"], "crit_nominal"]
        row["crit_matched"] = cm
        row["reject_matched"] = float(np.mean(v > cm))
        row["crit_nominal"] = cn
        row["reject_nominal"] = float(np.mean(v > cn)) if np.isfinite(cn) else np.nan
        lo, hi = wilson(int(np.sum(v > cm)), v.size)
        row["reject_matched_lo"], row["reject_matched_hi"] = lo, hi
    df = pd.DataFrame(rows)
    df["kind"] = df.kind.replace({"null": "no change"})  # "null" would be read back as a missing value by pandas
    df.to_csv(OUT / fname, index=False)
    pd.set_option("display.width", 250)
    if (df.kind == "change").any():
        print(df.pivot_table(index="method", columns=["scenario", "ratio"], values="reject_matched").round(2).T.to_string())
        print(df.pivot_table(index="method", columns=["scenario", "ratio"], values="cp_median_abs_error").round(0).T.to_string())
    else:
        print(df.pivot_table(index="method", columns="scenario", values="reject_matched").round(3).to_string())
        print(df.pivot_table(index="method", columns="scenario", values="reject_nominal").round(3).to_string())


# ----------------------------------------------------------------------------
# Q13b and Q13c: online use
# ----------------------------------------------------------------------------

ONLINE_METHODS = ("e-detector, plug-in", "e-detector, CS baseline", "BOCPD, matched", "BOCPD, p >= 0.5",
                  "BJV rolling, matched", "BJV rolling, nominal 5% per day", "BJV rolling, studentised, matched",
                  "BOCPD on truncated returns, matched")


def clip_by_calibration(r_full, cal_len, c=3.0):
    """Clip every return at tau = c x scaled MAD of the calibration window (the e-detector's truncation)."""
    cal = r_full[:, :cal_len]
    tau = c * np.median(np.abs(cal), axis=1, keepdims=True) * MAD_NORMAL_SCALE
    return np.clip(r_full, -tau, tau)


def online_paths(r_full, cal_len):
    """Streams for the online competitors on paths r_full = [calibration | monitored]."""
    cal, mon = r_full[:, :cal_len], r_full[:, cal_len:]
    _, S, L = C.bocpd_paths(mon, H_SIM, cal=cal, alpha0=ALPHA0, r_max=R_MAX, short=W_SHORT, return_log_odds=True)
    # Our modification, not in Adams and MacKay: BOCPD on returns clipped at the e-detector's truncation level
    rc = clip_by_calibration(r_full, cal_len)
    _, St, Lt = C.bocpd_paths(rc[:, cal_len:], H_SIM, cal=rc[:, :cal_len], alpha0=ALPHA0, r_max=R_MAX, short=W_SHORT,
                              return_log_odds=True)
    hist = r_full[:, cal_len - 2 * K_ONLINE :]
    Dz = C.bjv_online_statistic(hist, K_ONLINE, TRUNC)[:, 2 * K_ONLINE :]
    De = C.bjv_online_statistic(hist, K_ONLINE, TRUNC, norm="empirical")[:, 2 * K_ONLINE :]
    return L, S, Dz, De, Lt, St


def online_calib_job(args):
    reps, seed, day, a_level, cal_len, k, hazard = args
    rng = np.random.default_rng(seed)
    r = gen("gaussian", reps, cal_len + day, rng)
    a, _ = edet_alarm(r[:, cal_len:], r[:, :cal_len], a_level)
    fa = float(np.mean(a >= 0))
    _, _, L = C.bocpd_paths(r[:, cal_len:], hazard, cal=r[:, :cal_len], alpha0=ALPHA0, r_max=R_MAX, short=W_SHORT,
                            return_log_odds=True)
    rc = clip_by_calibration(r, cal_len)
    _, _, Lt = C.bocpd_paths(rc[:, cal_len:], hazard, cal=rc[:, :cal_len], alpha0=ALPHA0, r_max=R_MAX, short=W_SHORT,
                             return_log_odds=True)
    Dz = C.bjv_online_statistic(r[:, cal_len - 2 * k :], k, TRUNC)[:, 2 * k :]
    De = C.bjv_online_statistic(r[:, cal_len - 2 * k :], k, TRUNC, norm="empirical")[:, 2 * k :]
    return fa, L.max(1), np.nanmax(np.abs(Dz), axis=1), np.nanmax(np.abs(De), axis=1), Lt.max(1)


def matched_thresholds(day, a_level, cal_len, k, hazard, reps, seed):
    """Thresholds giving the e-detector's per-day false alarm rate on clean Gaussian data (split over two processes)."""
    half = reps // 2
    with Pool(2) as pool:
        parts = pool.map(online_calib_job, [(half, seed, day, a_level, cal_len, k, hazard),
                                            (reps - half, seed + 1, day, a_level, cal_len, k, hazard)])
    fa = (parts[0][0] * half + parts[1][0] * (reps - half)) / reps
    pm = np.concatenate([p[1] for p in parts])
    dm = np.concatenate([p[2] for p in parts])
    em = np.concatenate([p[3] for p in parts])
    tm = np.concatenate([p[4] for p in parts])
    return dict(day=day, A=a_level, cal_len=cal_len, k=k, hazard=hazard, paths=reps, edetector_fa_per_day=fa,
                bocpd_threshold_log_odds=float(np.quantile(pm, 1 - fa)), bjv_threshold=float(np.quantile(dm, 1 - fa)),
                bjv_studentised_threshold=float(np.quantile(em, 1 - fa)),
                bocpd_truncated_threshold_log_odds=float(np.quantile(tm, 1 - fa)),
                bjv_nominal=C.mosum_crit(day, k, 0.05))


def online_alarms(r_full, th):
    """First alarm (index into the monitored part) and direction for every online method."""
    cal, mon = r_full[:, :CAL], r_full[:, CAL:]
    out = {}
    out["e-detector, plug-in"] = edet_alarm(mon, cal, A)
    out["e-detector, CS baseline"] = edet_alarm(mon, cal, A, baseline="cs")
    L, S, Dz, De, Lt, St = online_paths(r_full, CAL)
    rows = np.arange(r_full.shape[0])
    for name, h, LL, SS in (("BOCPD, matched", th["bocpd_threshold_log_odds"], L, S), ("BOCPD, p >= 0.5", 0.0, L, S),
                            ("BOCPD on truncated returns, matched", th["bocpd_truncated_threshold_log_odds"], Lt, St)):
        a = C.first_crossing(LL, h)
        out[name] = (a, np.where(a >= 0, SS[rows, np.maximum(a, 0)], 0.0))
    for name, h, Z in (("BJV rolling, matched", th["bjv_threshold"], Dz),
                       ("BJV rolling, nominal 5% per day", th["bjv_nominal"], Dz),
                       ("BJV rolling, studentised, matched", th["bjv_studentised_threshold"], De)):
        a = C.first_crossing(np.abs(Z), h)
        out[name] = (a, np.where(a >= 0, np.sign(Z[rows, np.maximum(a, 0)]), 0.0))
    return out


def null_job(args):
    name, reps, seed, th = args
    rng = np.random.default_rng(seed)
    r = gen(name, reps, CAL + 2 * DAY, rng)
    t0 = time.time()
    al = online_alarms(r, th)
    print(f"online null {name}: {time.time() - t0:.0f} s", flush=True)
    rows = []
    for m, (a, _) in al.items():
        k1 = int(np.sum((a >= 0) & (a < DAY)))
        k2 = int(np.sum((a >= 0) & (a < 2 * DAY)))
        lo, hi = wilson(k1, reps)
        rows.append(dict(scenario=name, method=m, paths=reps, fa_one_day=k1 / reps, fa_one_day_lo=lo, fa_one_day_hi=hi,
                         fa_two_days=k2 / reps))
    return rows


def online_null():
    th = matched_thresholds(DAY, A, CAL, K_ONLINE, H_SIM, 4000, 13300)
    pd.DataFrame([th]).to_csv(OUT / "e13b_thresholds.csv", index=False)
    print(th, flush=True)
    jobs = [(nm, 1000, 13310 + i, th) for i, nm in enumerate(NULL_SCEN)]
    with Pool(2) as pool:
        rows = sum(pool.map(null_job, jobs), [])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e13b_online_null.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.pivot(index="method", columns="scenario", values="fa_one_day").round(3).to_string())


NU_ON, N_ON = 500, 2000


def change_job(args):
    name, ratio, reps, seed, th = args
    rng = np.random.default_rng(seed)
    r = gen(name, reps, CAL + N_ON, rng, nu=CAL + NU_ON, ratio=ratio)
    t0 = time.time()
    al = online_alarms(r, th)
    print(f"online change {name} {ratio}: {time.time() - t0:.0f} s", flush=True)
    rows = []
    want = 1.0 if ratio > 1 else -1.0
    for m, (a, d) in al.items():
        s = D.summarize_alarms(a, NU_ON, N_ON)
        det = a >= NU_ON
        rows.append(dict(scenario=name, ratio=ratio, method=m, paths=reps, **s,
                         correct_direction=float(np.mean(d[det] == want)) if det.any() else np.nan,
                         detected_within_390=float(np.mean(det & (a < NU_ON + DAY)))))
    return rows


def online_change(part: int):
    th = pd.read_csv(OUT / "e13b_thresholds.csv").iloc[0].to_dict()
    cells = [(nm, ratio) for nm in CHANGE_SCEN for ratio in RATIOS]
    seeds = {c: 13400 + i for i, c in enumerate(cells)}
    mine = cells[part::2]
    rows = []
    for nm, ratio in mine:
        rows += change_job((nm, ratio, 200, seeds[(nm, ratio)], th))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"e13c_online_change_part{part}.csv", index=False)
    print(df.round(3).to_string())


# ----------------------------------------------------------------------------
# Q13d: real data
# ----------------------------------------------------------------------------


def wilson(k: int, n: int, z: float = 1.959964):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return centre - half, centre + half


def real_calib():
    rows = [matched_thresholds(1440, 28800.0, 1440, K_ONLINE, 1 / 28800.0, 2000, 13500),
            matched_thresholds(375, 7500.0, 750, K_ONLINE, 1 / 7500.0, 2000, 13510)]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e13d_real_thresholds.csv", index=False)
    print(df.to_string())


class Validation:
    """Robust volatility before and after each time, and an EWMA forecast, for one deseasonalised series.

    Same definitions as E6b (``e06b_five_coins_monitor.Series``) with a general day length: entry t of
    ``y`` is the scaled median |u| over (t, t + day], of ``prev`` over (t - day, t - 1], and the forecast is
    an EWMA of |u| (half-life one day) times the median ratio y / EWMA over the previous 30 days, using
    only data up to t.
    """

    def __init__(self, u, day, min_frac=0.7):
        self.n, self.day = u.size, day
        v1 = pd.Series(np.abs(u)).rolling(day, min_periods=int(min_frac * day)).median().to_numpy() * MAD_NORMAL_SCALE
        self.v1 = v1
        self.y = np.full(self.n, np.nan)
        self.y[: self.n - day] = v1[day:]
        self.prev = np.full(self.n, np.nan)
        self.prev[1:] = v1[:-1]
        self.ewma = pd.Series(np.abs(u)).ewm(halflife=day, ignore_na=True).mean().to_numpy()

    def forecast(self, t):
        t = np.atleast_1d(np.asarray(t, int))
        js = np.arange(1, 31)
        s = t[:, None] - self.day * js[None, :]
        valid = s >= 0
        s = np.clip(s, 0, None)
        ratio = np.where(valid, self.y[s] / self.ewma[s], np.nan)
        with np.errstate(all="ignore"):
            k = np.nanmedian(np.where(np.isfinite(ratio), ratio, np.nan), axis=1)
        enough = np.sum(np.isfinite(ratio), axis=1) >= 10
        return np.where(enough, k * self.ewma[t], np.nan)

    def score(self, events):
        t = np.array([e[0] for e in events], int)
        d = np.array([e[1] for e in events])
        out = {}
        F = self.forecast(t) if t.size else np.array([])
        y = self.y[t] if t.size else np.array([])
        pv = self.prev[t] if t.size else np.array([])
        for dn, cmp in (("up", np.greater), ("down", np.less)):
            m = (d == dn) if t.size else np.array([], bool)
            ok_f = m & np.isfinite(F) & np.isfinite(y)
            ok_p = m & np.isfinite(pv) & np.isfinite(y)
            out[f"n_{dn}"] = int(m.sum())
            out[f"{dn}_n_forecast"] = int(ok_f.sum())
            out[f"{dn}_hits_forecast"] = int(np.sum(cmp(y[ok_f], F[ok_f])))
            out[f"{dn}_n_prev"] = int(ok_p.sum())
            out[f"{dn}_hits_prev"] = int(np.sum(cmp(y[ok_p], pv[ok_p])))
            out[f"{dn}_median_log_y_over_forecast"] = float(np.median(np.log(y[ok_f] / F[ok_f]))) if ok_f.any() else np.nan
        return out

    def random_times(self, start, rng, n=3000):
        t = rng.integers(start + self.day, self.n - self.day - 1, n)
        F, y, pv = self.forecast(t), self.y[t], self.prev[t]
        okf = np.isfinite(F) & np.isfinite(y)
        okp = np.isfinite(pv) & np.isfinite(y)
        return dict(n_forecast=int(okf.sum()), share_above_forecast=float(np.mean(y[okf] > F[okf])),
                    share_below_forecast=float(np.mean(y[okf] < F[okf])), n_prev=int(okp.sum()),
                    share_above_prev=float(np.mean(y[okp] > pv[okp])), share_below_prev=float(np.mean(y[okp] < pv[okp])))


def trailing_robust_var(u, day):
    """Predictable prior variance for BOCPD: squared scaled median |u| over the previous ``day`` returns."""
    v = pd.Series(np.abs(u)).rolling(day, min_periods=day // 2).median().to_numpy() * MAD_NORMAL_SCALE
    out = np.full(u.size, np.nan)
    out[1:] = v[:-1] ** 2
    return out


def real_bocpd(um, prior, hazard, cal_len, truncated):
    """BOCPD on (rows of) monitored real returns um; prior: predictable prior variance per step (same shape).

    ``truncated=True`` clips each return at 3 times the square root of that predictable robust variance
    (the e-detector's truncation level, computed from the previous day instead of a calibration window).
    """
    x = np.clip(um, -3 * np.sqrt(prior), 3 * np.sqrt(prior)) if truncated else um
    x = np.where(np.isfinite(um), x, np.nan)
    if truncated:  # the first day has no trailing prior yet: use the calibration window's own robust scale
        sc = np.nanmedian(np.abs(um[:, :cal_len]), axis=1, keepdims=True) * MAD_NORMAL_SCALE
        x[:, :cal_len] = np.clip(um[:, :cal_len], -3 * sc, 3 * sc)
    return C.bocpd_paths(x, hazard, cal=x[:, :cal_len], alpha0=ALPHA0, r_max=R_MAX, short=W_SHORT, prior_var_path=prior,
                         return_log_odds=True)


def run_real_series(label, u, i0, day, a_level, cal_len, hazard, th, bands=(0.0, 0.25), batch_P=None, batch_T=None):
    """All methods on one deseasonalised series u; monitoring starts at i0 (calibration first, as in E6b)."""
    from avmonitor.monitor import VolatilityMonitor

    um = u[i0:]
    n = um.size
    start, dead = cal_len, cal_len  # first alarm possible after the calibration window; same pause after alarms
    events = {}
    for b in bands:
        mon = VolatilityMonitor(threshold=a_level, c=3.0, cal_len=cal_len, null_band=b)
        ev = mon.run(um)
        events[f"e-detector, band {b}"] = [(e.index, e.direction) for e in ev]
    Dz = C.bjv_online_statistic(um[None, :], K_ONLINE, TRUNC)[0]
    De = C.bjv_online_statistic(um[None, :], K_ONLINE, TRUNC, norm="empirical")[0]
    prior = trailing_robust_var(u, day)[i0:]
    if batch_P is None:
        P, S, L = (a[0] for a in real_bocpd(um[None, :], prior[None, :], hazard, cal_len, False))
    else:
        P, S, L = batch_P
    if batch_T is None:
        _, St, Lt = (a[0] for a in real_bocpd(um[None, :], prior[None, :], hazard, cal_len, True))
    else:
        _, St, Lt = batch_T
    hb = th["bocpd_threshold_log_odds"]
    events["BOCPD, Gaussian-matched"] = C.scan_alarms(np.where(S > 0, L, -np.inf), np.where(S < 0, L, -np.inf), hb, hb,
                                                      start, dead)
    events["BJV rolling, Gaussian-matched"] = C.scan_alarms(Dz, -Dz, th["bjv_threshold"], th["bjv_threshold"], start, dead)
    events["BJV rolling, nominal 5% per day"] = C.scan_alarms(Dz, -Dz, th["bjv_nominal"], th["bjv_nominal"], start, dead)
    he = th["bjv_studentised_threshold"]
    events["BJV rolling, studentised, Gaussian-matched"] = C.scan_alarms(De, -De, he, he, start, dead)
    ht = th["bocpd_truncated_threshold_log_odds"]
    events["BOCPD on truncated returns, Gaussian-matched"] = C.scan_alarms(np.where(St > 0, Lt, -np.inf),
                                                                           np.where(St < 0, Lt, -np.inf), ht, ht, start, dead)
    matched = []
    for b in bands:
        ref = events[f"e-detector, band {b}"]
        nu_ = sum(1 for e in ref if e[1] == "up")
        nd_ = len(ref) - nu_
        for nm, LL, SS in (("BOCPD", L, S), ("BOCPD on truncated returns", Lt, St)):
            Lu, Ld = np.where(SS > 0, LL, -np.inf), np.where(SS < 0, LL, -np.inf)
            hu, hd = C.match_alarm_counts(Lu, Ld, nu_, nd_, start, dead, lo_hi=(-50.0, 1000.0), iters=50)
            events[f"{nm}, matched to band {b}"] = C.scan_alarms(Lu, Ld, hu, hd, start, dead)
            matched.append(dict(series=label, band=b, method=f"{nm} (log-odds)", h_up=hu, h_down=hd, target_up=nu_,
                                target_down=nd_))
        for nm, Z in (("BJV rolling", Dz), ("BJV rolling, studentised", De)):
            hu, hd = C.match_alarm_counts(Z, -Z, nu_, nd_, start, dead, lo_hi=(0.0, 1000.0), iters=50)
            events[f"{nm}, matched to band {b}"] = C.scan_alarms(Z, -Z, hu, hd, start, dead)
            matched.append(dict(series=label, band=b, method=nm, h_up=hu, h_down=hd, target_up=nu_, target_down=nd_))
    # shift to positions in u
    events = {m: [(t + i0, d) for t, d in ev] for m, ev in events.items()}
    return events, matched, n


def summarise_real(label, u, i0, day, events, n_mon, seed):
    V = Validation(u, day)
    rng = np.random.default_rng(seed)
    rows = []
    for m, ev in events.items():
        rec = dict(series=label, method=m, days=n_mon / day, alarms=len(ev), alarms_per_day=len(ev) / (n_mon / day))
        rec.update(V.score(ev))
        rows.append(rec)
    rnd = dict(series=label, method="random times", days=n_mon / day, **V.random_times(i0, rng))
    return rows, rnd


def real_crypto():
    from e06b_five_coins_monitor import COINS, prepare

    th = pd.read_csv(OUT / "e13d_real_thresholds.csv").set_index("day").loc[1440].to_dict()
    prepared = {c: prepare(c, "Jul-Dec 2025") for c in COINS}
    # BOCPD for the five coins in one batch (same length, monitoring from 1 Jul 2025)
    lens = {c: prepared[c][1].size - prepared[c][2] for c in COINS}
    assert len(set(lens.values())) == 1, lens
    um = np.stack([prepared[c][1][prepared[c][2]:] for c in COINS])
    prior = np.stack([trailing_robust_var(prepared[c][1], 1440)[prepared[c][2]:] for c in COINS])
    t0 = time.time()
    P, S, L = real_bocpd(um, prior, 1 / 28800.0, 1440, False)
    print(f"BOCPD batch: {time.time() - t0:.0f} s for {um.shape}", flush=True)
    _, St, Lt = real_bocpd(um, prior, 1 / 28800.0, 1440, True)
    print(f"BOCPD on truncated returns, batch: {time.time() - t0:.0f} s", flush=True)
    rows, rnd, matched, ev_rows = [], [], [], []
    for j, c in enumerate(COINS):
        idx, u, i0 = prepared[c]
        events, mt, n_mon = run_real_series(c, u, i0, 1440, 28800.0, 1440, 1 / 28800.0, th, batch_P=(P[j], S[j], L[j]),
                                            batch_T=(None, St[j], Lt[j]))
        r1, r2 = summarise_real(c, u, i0, 1440, events, n_mon, 1300 + j)
        rows += r1
        rnd.append(r2)
        matched += mt
        for m, ev in events.items():
            ev_rows += [dict(series=c, method=m, time=idx[t], direction=d) for t, d in ev]
        x = u[i0:]
        x = x[np.isfinite(x)]
        rnd[-1]["lag1_autocorr"] = float(np.corrcoef(x[1:], x[:-1])[0, 1])
        print(f"done {c}", flush=True)
    save_real("crypto", rows, rnd, matched, ev_rows)


def real_crypto_h1_2025():
    """Same as real_crypto on a second, earlier half-year (1 Jan to 30 Jun 2025), with the same frozen settings."""
    from e06b_five_coins_monitor import COINS

    from avmonitor.data import intraday_scale, load_binance

    th = pd.read_csv(OUT / "e13d_real_thresholds.csv").set_index("day").loc[1440].to_dict()
    prepared = {}
    for c, sym in COINS.items():
        df = load_binance(sym, "2024-01-01", "2025-07-01", download=False)
        s = intraday_scale(df.logret, lookback_days=20)
        u = (df.logret / s).to_numpy(copy=True)
        u[~np.isfinite(u)] = np.nan
        i0 = int(np.searchsorted(df.index, pd.Timestamp("2025-01-01", tz="UTC")))
        prepared[c] = (df.index, u, i0)
    lens = {c: prepared[c][1].size - prepared[c][2] for c in COINS}
    assert len(set(lens.values())) == 1, lens
    um = np.stack([prepared[c][1][prepared[c][2]:] for c in COINS])
    prior = np.stack([trailing_robust_var(prepared[c][1], 1440)[prepared[c][2]:] for c in COINS])
    t0 = time.time()
    P, S, L = real_bocpd(um, prior, 1 / 28800.0, 1440, False)
    _, St, Lt = real_bocpd(um, prior, 1 / 28800.0, 1440, True)
    print(f"BOCPD batches: {time.time() - t0:.0f} s for {um.shape}", flush=True)
    rows, rnd, matched, ev_rows = [], [], [], []
    for j, c in enumerate(COINS):
        idx, u, i0 = prepared[c]
        events, mt, n_mon = run_real_series(c, u, i0, 1440, 28800.0, 1440, 1 / 28800.0, th, batch_P=(P[j], S[j], L[j]),
                                            batch_T=(None, St[j], Lt[j]))
        r1, r2 = summarise_real(c, u, i0, 1440, events, n_mon, 1500 + j)
        rows += r1
        rnd.append(r2)
        matched += mt
        for m, ev in events.items():
            ev_rows += [dict(series=c, method=m, time=idx[t], direction=d) for t, d in ev]
        print(f"done {c}", flush=True)
    save_real("crypto_h1_2025", rows, rnd, matched, ev_rows)


def real_nifty():
    from e07_nifty import deseasonalise, load_clean

    th = pd.read_csv(OUT / "e13d_real_thresholds.csv").set_index("day").loc[375].to_dict()
    df, info, _ = load_clean("2023-10-01", "2026-01-01")
    u = deseasonalise(df).to_numpy()
    i0 = int(np.searchsorted(df.index, pd.Timestamp("2024-01-01", tz="Asia/Kolkata")))
    t0 = time.time()
    events, mt, n_mon = run_real_series("Nifty 50", u, i0, 375, 7500.0, 750, 1 / 7500.0, th)
    print(f"Nifty: {time.time() - t0:.0f} s, sessions {df.session[i0:].nunique()}", flush=True)
    rows, rnd = summarise_real("Nifty 50", u, i0, 375, events, n_mon, 1400)
    x = u[i0:]
    ok = np.isfinite(x[1:]) & np.isfinite(x[:-1])
    rnd["lag1_autocorr"] = float(np.corrcoef(x[1:][ok], x[:-1][ok])[0, 1])
    rnd["sessions"] = int(df.session[i0:].nunique())
    ev_rows = []
    for m, ev in events.items():
        ev_rows += [dict(series="Nifty 50", method=m, time=df.index[t], direction=d) for t, d in ev]
    save_real("nifty", rows, [rnd], mt, ev_rows)


def save_real(tag, rows, rnd, matched, ev_rows):
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"e13d_real_{tag}_summary.csv", index=False)
    pd.DataFrame(rnd).to_csv(OUT / f"e13d_real_{tag}_random.csv", index=False)
    pd.DataFrame(matched).to_csv(OUT / f"e13d_real_{tag}_matched_thresholds.csv", index=False)
    pd.DataFrame(ev_rows).to_csv(OUT / f"e13d_real_{tag}_events.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.round(3).to_string())
    print(pd.DataFrame(rnd).round(3).to_string())


# ----------------------------------------------------------------------------
# Q13e: cost per step
# ----------------------------------------------------------------------------


def cost():
    import platform

    from avmonitor.monitor import VolatilityMonitor

    rng = np.random.default_rng(13600)
    n = 20000
    r = rng.standard_normal(n)
    rows = []
    t = time.perf_counter()
    VolatilityMonitor(threshold=28800, cal_len=1440).run(r)
    rows.append(dict(method="e-detector, VolatilityMonitor.update (two-sided, 11 bets)", steps=n, seconds=time.perf_counter() - t))
    for rmax in (256, 512):
        f = C.BOCPDVariance(1, 1 / 28800.0, ALPHA0, 1.0, rmax, min(W_SHORT, rmax - 2))
        t = time.perf_counter()
        for x in r:
            f.update(x)
        rows.append(dict(method=f"BOCPD, BOCPDVariance.update (run-length cap {rmax})", steps=n, seconds=time.perf_counter() - t))
    k = K_ONLINE
    t = time.perf_counter()
    for i in range(2 * k, n):
        C.bjv_overlap_statistics(r[None, i - 2 * k : i + 1], k, TRUNC)
    rows.append(dict(method=f"BJV rolling, one window pair per step (k={k}, truncated)", steps=n - 2 * k, seconds=time.perf_counter() - t))
    t = time.perf_counter()
    C.bjv_online_statistic(r[None, :], k, TRUNC)
    rows.append(dict(method=f"BJV rolling, vectorised over the whole series (k={k}, truncated)", steps=n, seconds=time.perf_counter() - t))
    t = time.perf_counter()
    C.bjv_online_statistic(r[None, :], k, None)
    rows.append(dict(method=f"BJV rolling, vectorised, plain (cumulative sums, k={k})", steps=n, seconds=time.perf_counter() - t))
    df = pd.DataFrame(rows)
    df["microseconds_per_step"] = df.seconds / df.steps * 1e6
    df["machine"] = f"{platform.machine()}, Python {platform.python_version()}, 4 cores shared with other jobs"
    df.to_csv(OUT / "e13e_cost.csv", index=False)
    print(df.round(2).to_string())


# ----------------------------------------------------------------------------
# Tables and figures
# ----------------------------------------------------------------------------


def pooled_real():
    parts = []
    labels = {"crypto": "five Binance coins, Jul to Dec 2025", "crypto_h1_2025": "five Binance coins, Jan to Jun 2025", "nifty": "Nifty 50"}
    for tag, lab in labels.items():
        p = OUT / f"e13d_real_{tag}_summary.csv"
        if p.exists():
            d = pd.read_csv(p)
            d["market"] = lab
            parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    g = df.groupby(["market", "method"], sort=False).agg(
        days=("days", "sum"), alarms=("alarms", "sum"), n_up=("n_up", "sum"), n_down=("n_down", "sum"),
        up_n_forecast=("up_n_forecast", "sum"), up_hits_forecast=("up_hits_forecast", "sum"),
        down_n_forecast=("down_n_forecast", "sum"), down_hits_forecast=("down_hits_forecast", "sum"),
        up_n_prev=("up_n_prev", "sum"), up_hits_prev=("up_hits_prev", "sum"),
        down_n_prev=("down_n_prev", "sum"), down_hits_prev=("down_hits_prev", "sum")).reset_index()
    g["alarms_per_day"] = g.alarms / g.days
    for dn in ("up", "down"):
        for v in ("forecast", "prev"):
            k, nn = g[f"{dn}_hits_{v}"], g[f"{dn}_n_{v}"]
            g[f"{dn}_rate_{v}"] = k / nn.replace(0, np.nan)
            ci = [wilson(int(a), int(b)) for a, b in zip(k, nn)]
            g[f"{dn}_{v}_lo"] = [c[0] for c in ci]
            g[f"{dn}_{v}_hi"] = [c[1] for c in ci]
    # two-proportion z-test of each matched competitor against the e-detector band it was matched to
    pv = {f"{dn}_p_vs_edetector": [] for dn in ("up", "down")}
    for _, r in g.iterrows():
        ref_name = None
        for b in ("0.0", "0.25"):
            if r.method.endswith(f"matched to band {b}"):
                ref_name = f"e-detector, band {b}"
        for dn in ("up", "down"):
            if ref_name is None:
                pv[f"{dn}_p_vs_edetector"].append(np.nan)
                continue
            ref = g[(g.market == r.market) & (g.method == ref_name)].iloc[0]
            k1, n1, k2, n2 = ref[f"{dn}_hits_forecast"], ref[f"{dn}_n_forecast"], r[f"{dn}_hits_forecast"], r[f"{dn}_n_forecast"]
            if n1 == 0 or n2 == 0:
                pv[f"{dn}_p_vs_edetector"].append(np.nan)
                continue
            pp = (k1 + k2) / (n1 + n2)
            se = np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
            pv[f"{dn}_p_vs_edetector"].append(float(2 * stats.norm.sf(abs(k1 / n1 - k2 / n2) / se)) if se > 0 else np.nan)
    for c, v in pv.items():
        g[c] = v
    g.to_csv(OUT / "e13d_real_pooled.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["market", "method", "alarms_per_day", "n_up", "up_rate_forecast", "up_forecast_lo", "up_forecast_hi",
            "up_p_vs_edetector", "n_down", "down_rate_forecast", "down_forecast_lo", "down_forecast_hi",
            "down_p_vs_edetector", "up_rate_prev", "down_rate_prev"]
    print(g[cols].round(3).to_string(index=False))
    return g


COLORS = {"e-detector, plug-in": "#2a78d6", "e-detector, CS baseline": "#1baf7a", "BOCPD, matched": "#eb6834",
          "BOCPD, p >= 0.5": "#f2a477", "BJV rolling, matched": "#7a52c7", "BJV rolling, nominal 5% per day": "#b9a3e3",
          "BJV rolling, studentised, matched": "#c43c8a", "BOCPD on truncated returns, matched": "#a83c0c"}


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#e6e5e1", lw=0.5)


def figures():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(exist_ok=True)
    nd = pd.read_csv(OUT / "e13b_online_null.csv")
    fig, ax = plt.subplots(figsize=(11, 4))
    scen = list(NULL_SCEN)
    meths = list(ONLINE_METHODS)
    w = 0.8 / len(meths)
    for i, m in enumerate(meths):
        sub = nd[nd.method == m].set_index("scenario").reindex(scen)
        x = np.arange(len(scen)) + (i - (len(meths) - 1) / 2) * w
        ax.bar(x, sub.fa_one_day, width=w, color=COLORS[m], label=m)
        ax.errorbar(x, sub.fa_one_day, yerr=np.clip([sub.fa_one_day - sub.fa_one_day_lo, sub.fa_one_day_hi - sub.fa_one_day], 0, None),
                    fmt="none", ecolor="#0b0b0b", elinewidth=0.6, capsize=1.5)
    ax.axhline(0.05, color="#0b0b0b", ls="--", lw=0.8)
    ax.text(len(scen) - 0.5, 0.06, "e-detector bound (5% per day)", ha="right", fontsize=8)
    ax.set_xticks(np.arange(len(scen)))
    ax.set_xticklabels(scen)
    ax.set_ylabel("false alarm within one 390-step day")
    ax.set_title("No change: false alarms per day [SYNTHETIC]. Competitor thresholds matched to the e-detector on clean Gaussian data",
                 fontsize=9, loc="left")
    ax.legend(fontsize=7, ncol=3, frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "e13_online_false_alarms.png", dpi=130)
    plt.close(fig)

    parts = [pd.read_csv(p) for p in sorted(OUT.glob("e13c_online_change_part*.csv"))]
    if parts:
        cd = pd.concat(parts, ignore_index=True)
        cd.to_csv(OUT / "e13c_online_change.csv", index=False)
        fig, axes = plt.subplots(2, len(CHANGE_SCEN), figsize=(14, 6), sharey="row")
        for j, sc in enumerate(CHANGE_SCEN):
            for i, (lab, sel) in enumerate((("increase", lambda d: d.ratio > 1), ("decrease", lambda d: d.ratio < 1))):
                ax = axes[i, j]
                for m in meths:
                    s = cd[(cd.scenario == sc) & (cd.method == m) & sel(cd)].sort_values("ratio")
                    x = s.ratio if lab == "increase" else 1 / s.ratio
                    ax.plot(x, s.median_delay, marker="o", ms=3, lw=1, color=COLORS[m], label=m)
                    low = s.detected < 0.8
                    ax.scatter(x[low], s.median_delay[low], s=30, marker="x", color=COLORS[m], lw=1.2)
                    bad = s.false_alarm_before_change > 0.2
                    ax.scatter(x[bad], s.median_delay[bad], s=40, facecolors="none", edgecolors="#d03b3b", lw=1)
                ax.set_yscale("log")
                ax.set_title(f"{sc}, {lab}", fontsize=8, loc="left")
                ax.set_xlabel("sd ratio" if lab == "increase" else "1 / sd ratio", fontsize=8)
                _style(ax)
            axes[0, 0].set_ylabel("median delay (steps)")
            axes[1, 0].set_ylabel("median delay (steps)")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.955))
        fig.suptitle("Online detection delay among detected paths [SYNTHETIC]. Crosses: fewer than 80% of paths detected the change. "
                     "Red rings: more than 20% raised a false alarm before it", fontsize=9, x=0.01, ha="left")
        fig.tight_layout(rect=(0, 0, 1, 0.9))
        fig.savefig(FIG / "e13_online_delay.png", dpi=130)
        plt.close(fig)

    off = pd.read_csv(OUT / "e13a_offline_power.csv")
    ch = off[off.kind == "change"]
    om = list(ch.method.unique())
    pal = ["#7a52c7", "#9b7fd8", "#c9b8ec", "#4b2a9a", "#eb6834", "#2a78d6"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for ax, sc in zip(axes, ("gaussian", "jumps", "student t5")):
        for c_, m in zip(pal, om):
            s = ch[(ch.scenario == sc) & (ch.method == m)].sort_values("ratio")
            ax.plot(np.log(s.ratio), s.reject_matched, marker="o", ms=3, lw=1, color=c_, label=m)
        ax.set_xticks(np.log(RATIOS))
        ax.set_xticklabels([str(x) for x in RATIOS], fontsize=7)
        ax.set_title(f"{sc}", fontsize=9, loc="left")
        ax.set_xlabel("sd ratio (log scale)", fontsize=8)
        _style(ax)
    axes[0].set_ylabel("power at 5% size (matched on Gaussian)")
    axes[0].legend(fontsize=6, frameon=False)
    fig.suptitle("Offline test on a 2000-step window, change at step 1000 [SYNTHETIC]", fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "e13_offline_power.png", dpi=130)
    plt.close(fig)

    g = pooled_real()
    markets = [m for m in ("five Binance coins, Jul to Dec 2025", "five Binance coins, Jan to Jun 2025", "Nifty 50") if m in set(g.market)]
    fig, axes = plt.subplots(1, len(markets), figsize=(6 * len(markets), 4.2), sharey=True)
    show = ["e-detector, band 0.25", "BOCPD, matched to band 0.25", "BOCPD on truncated returns, matched to band 0.25",
            "BJV rolling, matched to band 0.25", "BJV rolling, studentised, matched to band 0.25"]
    cols = ["#2a78d6", "#eb6834", "#a83c0c", "#7a52c7", "#c43c8a"]
    for ax, mk in zip(np.atleast_1d(axes), markets):
        gm = g[g.market == mk].set_index("method")
        for i, dn in enumerate(("up", "down")):
            for j, (m, c_) in enumerate(zip(show, cols)):
                if m not in gm.index:
                    continue
                r = gm.loc[m]
                x = i + (j - 2) * 0.16
                v, lo, hi = r[f"{dn}_rate_forecast"], r[f"{dn}_forecast_lo"], r[f"{dn}_forecast_hi"]
                ax.bar(x, v, width=0.15, color=c_, label=m if i == 0 else None)
                ax.errorbar(x, v, yerr=np.clip([[v - lo], [hi - v]], 0, None), fmt="none", ecolor="#0b0b0b", elinewidth=0.6,
                            capsize=1.5)
        ax.axhline(0.5, color="#9a9893", ls="--", lw=0.8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["after up alarms:\nnext day above EWMA forecast", "after down alarms:\nnext day below EWMA forecast"],
                           fontsize=8)
        ax.set_title(f"{mk} [REAL]", fontsize=9, loc="left")
        ax.set_ylim(0, 1)
        _style(ax)
    axes[0].set_ylabel("share of alarms (Wilson 95%)")
    axes[0].legend(fontsize=6, frameon=False, loc="lower left")
    fig.suptitle("Rate-matched validation: competitors tuned to the e-detector's numbers of up and down alarms", fontsize=9,
                 x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "e13_real_validation.png", dpi=130)
    plt.close(fig)


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "offline-null":
        offline_null()
    elif cmd == "offline-change":
        offline_change()
    elif cmd == "online-null":
        online_null()
    elif cmd == "online-change":
        online_change(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    elif cmd == "real-calib":
        real_calib()
    elif cmd == "real-crypto":
        real_crypto()
    elif cmd == "real-crypto-h1-2025":
        real_crypto_h1_2025()
    elif cmd == "real-nifty":
        real_nifty()
    elif cmd == "cost":
        cost()
    elif cmd == "figures":
        figures()
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
