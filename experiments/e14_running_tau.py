"""E14. A false alarm guarantee with a predictable, running truncation level (tracker B68).

Theory: docs/theory_running_tau.md. Code: src/avmonitor/running.py. Report: docs/report_parts/E14.md.

Methods compared (all two-sided: each side alarms at 2A; "proved" means a theorem covers the
false alarm bound under the stated null; "approx." means no proof):
  plugin          fixed tau = 3 x scaled MAD of the calibration window, plug-in mean (approx.), A = 7800
  plugin-run      running tau, plug-in mean at the running level (approx.), A = 7800
  cs              E8: tau from the first 25% of the window, confidence-sequence bounds (proved), A = 9750
  A-grid-fixed    route A, grid band, tau fixed at the grid centre (proved), A = 9750
  A-grid-run      route A, grid band, running tau (proved), A = 9750
  A-dkw-run       route A, Howard-Ramdas DKW band, running tau (proved), A = 9750
  B-lr            route B, rank bets (clipped Gaussian likelihood ratio in rank space) (proved, exact), A = 7800
  B-ind           route B, indicator bets on running quantiles (proved, exact), A = 7800
  B-freeze        B-lr with the reference set frozen once a statistic reaches 20 (heuristic, approx.), A = 7800
  oracle          true tau and true m(tau) (E9 reference; proved), A = 7800
  awake-A         route A with no blind period (learning band, staggered chains) (proved), A = 9750
  awake-B         route B restarted with no warm-up (proved, exact), A = 7800
With delta = 0.01 every "proved" bound is 5% for one 390-step day (390 / 7800 or 0.01 + 390 / 9750).

Research questions
  Q14a [SYNTHETIC] Validity at scale. With i.i.d. returns (Gaussian, Student t3, jumps, price
       ticks) and calibration windows of 390 and 1440 steps, what is the false alarm rate within
       one day (390 steps) and five days (1950) for every method, against the proved bound, and
       how often does each route A band fail?
  Q14b [SYNTHETIC] Delay and power. Standard deviation changes x1.25, x1.5, x2, x0.8, x0.67, x0.5
       500 steps into monitoring (Gaussian, t3, jumps; calibration 390 and 1440). Median delay,
       detection within 2500 steps, false alarms before the change.
  Q14c [SYNTHETIC] Seasonality and drift (the null is false): U-shaped intraday profile on raw
       returns, slow stochastic volatility, GARCH(1,1). Alarm rates per day with null band 0 and 0.25.
  Q14d [REAL] Five Binance coins (July to December 2025, deseasonalised as in E6) and Nifty 50
       (2024 to 2025, deseasonalised as in E7): alarms per day and the E5c check (is volatility in the
       next day above / below the baseline after an up / down alarm).
  Q14e [SYNTHETIC] Cost per step of every monitor.
  Q14f [SYNTHETIC] No blind period. (i) false alarms within one and five days of a restart for
       awake-A and awake-B; (ii) delay when the change comes 30, 120 or 360 steps after a restart,
       against the standard protocol (1440 steps of calibration after every alarm).
  Q14g [REAL] E12 rerun: BTC, ETH, SOL, 21 January 2024 to 31 December 2025, streaming design, band
       0.25 (E12's frozen setting). Share of time awake, awake at the 16 FOMC and 22 CPI releases,
       up alarm within 6 hours of a release, against matched control windows.

All settings (c = 3, delta = 0.01, the grid, the running window 120 refreshed every 10 steps, the
awake schedule) were fixed before the reported runs from first principles and two pilots on other
seeds (docs/iterations/TH2.md); nothing was tuned on the reported data. Seeds are fixed below.
  Q14w [SYNTHETIC] Supporting numbers: band widths at every grid level (grid, DKW time-uniform, DKW
       fixed-n, E8), constants, and the dilution of the rank detector after a change (numerical
       integration of the expected log growth).
Usage: python experiments/e14_running_tau.py a|b|c|d|e|f|g|w|figures [--smoke] [--cal=390] [--law=jumps]
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import running as RT
from avmonitor.confseq import cs_baseline
from avmonitor.monitor import VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE
from avmonitor.simulate import change_returns

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
sys.path.insert(0, str(ROOT / "experiments"))

SMOKE = "--smoke" in sys.argv
OPTS = dict(a.lstrip("-").split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)  # --cal=390, --law=jumps
DAY = 390
C = 3.0
DELTA = 0.01
A_PLUG = 7800.0
A_DELTA = DAY / (0.05 - DELTA)  # 9750
A_RANK = 7800.0
FREEZE_AT = 20.0
CAL_LENS = (390, 1440)
LAWS = ("gaussian", "student_t3", "jumps", "ticks")
SEED = 14_000

A_OF = {"oracle": A_PLUG, "plugin": A_PLUG, "plugin-run": A_PLUG, "cs": A_DELTA, "A-grid-fixed": A_DELTA,
        "A-grid-run": A_DELTA, "A-dkw-run": A_DELTA, "B-lr": A_RANK, "B-ind": A_RANK, "B-freeze": A_RANK,
        "awake-A": A_DELTA, "awake-B": A_RANK}
PROVED = {"oracle": True, "plugin": False, "plugin-run": False, "cs": True, "A-grid-fixed": True, "A-grid-run": True,
          "A-dkw-run": True, "B-lr": True, "B-ind": True, "B-freeze": False, "awake-A": True, "awake-B": True}


# ----------------------------------------------------------------------------
# Laws and the truth
# ----------------------------------------------------------------------------


def draw(law: str, size, rng):
    """The E9 laws: unit-scale Gaussian, Student t3 with unit variance, jumps (0.5% per step, sd 10), ticks of 0.5."""
    if law == "gaussian":
        return rng.standard_normal(size)
    if law == "student_t3":
        return rng.standard_t(3, size) / np.sqrt(3.0)
    if law == "jumps":
        return rng.standard_normal(size) + (rng.random(size) < 0.005) * rng.normal(0, 10, size)
    if law == "ticks":
        return np.round(rng.standard_normal(size) / 0.5) * 0.5
    raise ValueError(law)


class Truth:
    """m(s) = E[min(r^2, s^2)] from 4 million draws, and the E9 oracle level."""

    def __init__(self, law: str, n: int = 4_000_000, seed: int = 999):
        big = draw(law, n, np.random.default_rng(seed))
        self.v = np.sort(big**2)
        self.cs = np.concatenate([[0.0], np.cumsum(self.v)])
        tau = 3 * np.median(np.abs(big)) * MAD_NORMAL_SCALE
        self.tau = max(tau, 1.5) if law == "ticks" else tau

    def m(self, s):
        s2 = np.asarray(s, float) ** 2
        k = np.searchsorted(self.v, s2, side="right")
        return (self.cs[k] + s2 * (self.v.size - k)) / self.v.size


# ----------------------------------------------------------------------------
# Per-step inputs of the X-based methods
# ----------------------------------------------------------------------------


def x_inputs(method: str, r_cal, r_hist, n_mon: int, truth: Truth | None = None, band: float = 0.0):
    """(tau2, m_up, m_dn, band_failed) for the monitored steps of each path (arrays (reps, n_mon) or (reps, 1))."""
    reps, cal = r_cal.shape
    failed = np.zeros(reps, bool)
    if method == "oracle":
        t2 = np.full((reps, 1), truth.tau**2)
        m = np.full((reps, 1), float(truth.m(truth.tau)))
        return t2, m * (1 + band), m * (1 - band), failed
    if method == "plugin":
        tau = C * RT.scaled_mad(r_cal)
        m = np.mean(np.minimum(r_cal**2, (tau**2)[:, None]), axis=1)
        return (tau**2)[:, None], (m * (1 + band))[:, None], (m * (1 - band))[:, None], failed
    if method == "cs":
        b = cs_baseline(r_cal, c=C, delta_upper=DELTA / 2, delta_lower=DELTA / 2)
        if truth is not None:
            mt = truth.m(b.tau)
            failed = (mt > b.upper) | (mt < b.lower)
        return (b.tau**2)[:, None], (b.upper * (1 + band))[:, None], (b.lower * (1 - band))[:, None], failed
    kind = {"plugin-run": "plugin", "A-grid-fixed": "grid", "A-grid-run": "grid", "A-dkw-run": "dkw"}[method]
    tb = RT.calibrate_band(r_cal, method=kind, c=C, delta=DELTA)
    if truth is not None and kind != "plugin":
        mt = truth.m(tb.levels)
        failed = np.any((mt > tb.upper) | (mt < tb.lower), axis=1)
    tau2, mu, md, _ = RT.route_a_inputs(tb, r_hist, cal, n_mon, c=C, null_band=band, fixed=(method == "A-grid-fixed"))
    return tau2, mu, md, failed


def run_x(method, r_cal, r_mon, truth=None, band=0.0):
    hist = np.concatenate([r_cal, r_mon], axis=1)
    tau2, mu, md, failed = x_inputs(method, r_cal, hist, r_mon.shape[1], truth, band)
    x = np.minimum(r_mon**2, tau2)
    L = 2 * A_OF[method]
    au, ad = RT.two_sided_sr(x, tau2, mu, md, [L])[L]
    return au, ad, failed


def run_b(method, r_cal, r_mon, rng, band=0.0):
    a = np.abs(np.concatenate([r_cal, r_mon], axis=1))
    cal = r_cal.shape[1]
    L = 2 * A_OF[method]
    fam = "ind" if method == "B-ind" else "lr"
    out, _ = RT.conformal_detector(a, rng, cal, [L], family=fam, band=band,
                                   freeze_at=FREEZE_AT if method == "B-freeze" else None)
    au, ad = out[L]
    return np.where(au >= 0, au - cal, -1), np.where(ad >= 0, ad - cal, -1), np.zeros(a.shape[0], bool)


def run_method(method, r_cal, r_mon, rng, truth=None, band=0.0):
    if method.startswith("B-"):
        return run_b(method, r_cal, r_mon, rng, band)
    return run_x(method, r_cal, r_mon, truth, band)


def first_alarm(au, ad):
    big = np.iinfo(np.int64).max
    a = np.minimum(np.where(au >= 0, au, big), np.where(ad >= 0, ad, big))
    return np.where(a == big, -1, a)


def wilson(k, n, z=1.959964):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


# ----------------------------------------------------------------------------
# Q14a: validity at scale
# ----------------------------------------------------------------------------

METHODS_A = ("oracle", "plugin", "plugin-run", "cs", "A-grid-fixed", "A-grid-run", "A-dkw-run", "B-lr", "B-ind",
             "B-freeze")


def q14a():
    reps, chunk = (400, 200) if SMOKE else (5000, 1000)
    horizons = (DAY, 5 * DAY)
    n_mon = max(horizons)
    rows = []
    for li, law in enumerate(LAWS):
        truth = Truth(law)
        for ci, cal in enumerate(CAL_LENS):
            if "cal" in OPTS and int(OPTS["cal"]) != cal:
                continue
            t0 = time.time()
            acc = {m: [] for m in METHODS_A}
            fails = {m: [] for m in METHODS_A}
            for j in range(reps // chunk):
                rng = np.random.default_rng([SEED, 1, li, ci, j])
                r = draw(law, (chunk, cal + n_mon), rng)
                r_cal, r_mon = r[:, :cal], r[:, cal:]
                for m in METHODS_A:
                    au, ad, f = run_method(m, r_cal, r_mon, np.random.default_rng([SEED, 2, li, ci, j]), truth)
                    acc[m].append(first_alarm(au, ad))
                    fails[m].append(f)
            for m in METHODS_A:
                a = np.concatenate(acc[m])
                f = np.concatenate(fails[m])
                for h in horizons:
                    k = int(np.sum((a >= 0) & (a < h)))
                    lo, hi = wilson(k, a.size)
                    bound = (DELTA if PROVED[m] and m not in ("oracle",) and not m.startswith("B") else 0.0) + h / A_OF[m]
                    rows.append(dict(law=law, cal_len=cal, method=m, proved=PROVED[m], A=A_OF[m], horizon=h,
                                     bound=bound, false_alarm_rate=k / a.size, ci_low=lo, ci_high=hi, paths=a.size,
                                     band_failed=float(f.mean())))
            print(f"  14a {law:10s} cal {cal:5d}: {time.time() - t0:6.1f} s", flush=True)
    df = pd.DataFrame(rows)
    tag = f"_cal{OPTS['cal']}" if "cal" in OPTS else ""
    df.to_csv(OUT / f"e14a_validity{tag}{'_smoke' if SMOKE else ''}.csv", index=False)
    print(df[df.horizon == DAY].pivot_table(index=["law", "cal_len"], columns="method", values="false_alarm_rate")
          .round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Q14b: delay and power
# ----------------------------------------------------------------------------

METHODS_B = ("plugin", "cs", "A-grid-fixed", "A-grid-run", "A-dkw-run", "B-lr", "B-ind", "B-freeze")
RATIOS = (1.25, 1.5, 2.0, 0.8, 1 / 1.5, 0.5)
NU, N_MON_B = 500, 3000


def q14b():
    reps = 200 if SMOKE else 1000
    rows = []
    for li, law in enumerate(("gaussian", "student_t3", "jumps")):
        if "law" in OPTS and OPTS["law"] != law:
            continue
        for ci, cal in enumerate(CAL_LENS):
            t0 = time.time()
            for ri, ratio in enumerate(RATIOS):
                rng = np.random.default_rng([SEED, 3, li, ci, ri])
                r = draw(law, (reps, cal + N_MON_B), rng)
                r[:, cal + NU :] *= ratio
                r_cal, r_mon = r[:, :cal], r[:, cal:]
                side = 0 if ratio > 1 else 1
                for m in METHODS_B:
                    out = run_method(m, r_cal, r_mon, np.random.default_rng([SEED, 4, li, ci, ri]))
                    a_side, a_any = out[side], first_alarm(out[0], out[1])
                    early = (a_any >= 0) & (a_any < NU)
                    det = (a_side >= NU) & ~early
                    d = a_side[det] - NU + 1
                    rows.append(dict(law=law, cal_len=cal, sd_ratio=round(ratio, 3), method=m, A=A_OF[m],
                                     false_alarm_before_change=float(early.mean()), detected_within_2500=float(det.mean()),
                                     median_delay=float(np.median(d)) if d.size else np.nan,
                                     mean_delay=float(d.mean()) if d.size else np.nan, paths=reps))
            print(f"  14b {law:10s} cal {cal:5d}: {time.time() - t0:6.1f} s", flush=True)
    df = pd.DataFrame(rows)
    tag = f"_{OPTS['law']}" if "law" in OPTS else ""
    df.to_csv(OUT / f"e14b_delay{tag}{'_smoke' if SMOKE else ''}.csv", index=False)
    print(df.pivot_table(index=["law", "cal_len", "sd_ratio"], columns="method", values="median_delay").round(0).to_string())
    return df


# ----------------------------------------------------------------------------
# Q14c: seasonality and drift (null false)
# ----------------------------------------------------------------------------

METHODS_C = ("plugin", "cs", "A-grid-fixed", "A-grid-run", "B-lr", "B-ind")


def garch(reps, n, alpha, beta, rng, burn=5000):
    omega = 1 - alpha - beta
    s2 = np.ones(reps)
    prev = rng.standard_normal(reps)
    out = np.empty((reps, n))
    for t in range(burn + n):
        s2 = omega + alpha * prev**2 + beta * s2
        prev = np.sqrt(s2) * rng.standard_normal(reps)
        if t >= burn:
            out[:, t - burn] = prev
    return out


def q14c():
    reps = 200 if SMOKE else 1000
    cal, n_mon = 1440, 5 * DAY
    rows = []
    scen = {
        "intraday U-shape (amplitude 0.6), raw": lambda rng: change_returns(reps, cal + n_mon, seasonal_amplitude=0.6,
                                                                           period=DAY, rng=rng).r,
        "slow SV (log-vol sd 0.15, half-life 390)": lambda rng: change_returns(reps, cal + n_mon, sv_xi=0.15,
                                                                               sv_halflife=390.0, rng=rng).r,
        "GARCH(1,1) 0.05 / 0.94": lambda rng: garch(reps, cal + n_mon, 0.05, 0.94, rng),
    }
    for si, (name, gen) in enumerate(scen.items()):
        t0 = time.time()
        r = gen(np.random.default_rng([SEED, 5, si]))
        r_cal, r_mon = r[:, :cal], r[:, cal:]
        for band in (0.0, 0.25):
            for m in METHODS_C:
                au, ad, _ = run_method(m, r_cal, r_mon, np.random.default_rng([SEED, 6, si]), band=band)
                a = first_alarm(au, ad)
                for h in (DAY, 5 * DAY):
                    rows.append(dict(scenario=name, band=band, method=m, horizon=h, alarm_rate=float(np.mean((a >= 0) & (a < h))),
                                     up_share=float(np.mean((au >= 0) & (au < h))), down_share=float(np.mean((ad >= 0) & (ad < h))),
                                     paths=reps))
        print(f"  14c {name}: {time.time() - t0:6.1f} s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"e14c_seasonality_drift{'_smoke' if SMOKE else ''}.csv", index=False)
    print(df[df.horizon == DAY].pivot_table(index=["scenario", "band"], columns="method", values="alarm_rate").round(3).to_string())
    return df


# ----------------------------------------------------------------------------
# Q14f: no blind period (synthetic)
# ----------------------------------------------------------------------------


def awake_b(r_seg, rng, L, band=0.0):
    out, _ = RT.conformal_detector(np.abs(r_seg), rng, 0, [L], family="lr", band=band)
    return out[L]


def q14f():
    rows_null, rows_delay = [], []
    reps, chunk = (400, 200) if SMOKE else (4000, 1000)
    pre, n = RT.AWAKE_PRE, 5 * DAY
    sub = OPTS.get("sub", "all")
    for li, law in enumerate(LAWS if sub in ("all", "null") else ()):
        t0 = time.time()
        acc = {"awake-A": [], "awake-B": []}
        for j in range(reps // chunk):
            rng = np.random.default_rng([SEED, 7, li, j])
            r = draw(law, (chunk, pre + n), rng)
            out, _ = RT.awake_route_a(r[:, :pre], r[:, pre:], [2 * A_DELTA], delta=DELTA, c=C)
            acc["awake-A"].append(first_alarm(*out[2 * A_DELTA]))
            acc["awake-B"].append(first_alarm(*awake_b(r[:, pre:], np.random.default_rng([SEED, 8, li, j]), 2 * A_RANK)))
        for m, v in acc.items():
            a = np.concatenate(v)
            for h in (DAY, 5 * DAY):
                k = int(np.sum((a >= 0) & (a < h)))
                lo, hi = wilson(k, a.size)
                rows_null.append(dict(law=law, method=m, A=A_OF[m], horizon=h,
                                      bound=(DELTA if m == "awake-A" else 0.0) + h / A_OF[m], false_alarm_rate=k / a.size,
                                      ci_low=lo, ci_high=hi, paths=a.size))
        print(f"  14f null {law}: {time.time() - t0:6.1f} s", flush=True)
    # delay after a restart: standard protocol (calibrate 1440 steps after every alarm) against awake monitors
    reps = 200 if SMOKE else 1000
    cal, n_tot = 1440, 1440 + 2500
    for li, law in enumerate(("gaussian", "student_t3") if sub in ("all", "delay") else ()):
        for ni, nu in enumerate((30, 120, 360)):
            t0 = time.time()
            for ri, ratio in enumerate((2.0, 1.5, 0.5, 1 / 1.5)):
                rng = np.random.default_rng([SEED, 9, li, ni, ri])
                r = draw(law, (reps, pre + n_tot), rng)
                r[:, pre + nu :] *= ratio
                rp, seg = r[:, :pre], r[:, pre:]
                side = 0 if ratio > 1 else 1
                res = {}
                out, _ = RT.awake_route_a(rp, seg, [2 * A_DELTA], delta=DELTA, c=C)
                res["awake-A"] = out[2 * A_DELTA]
                res["awake-B"] = awake_b(seg, np.random.default_rng([SEED, 10, li, ni, ri]), 2 * A_RANK)
                for m in ("plugin", "cs", "A-grid-run", "B-lr"):
                    au, ad, _ = run_method(m, seg[:, :cal], seg[:, cal:], np.random.default_rng([SEED, 11, li, ni, ri]))
                    res[m + " (1440-step recalibration)"] = (np.where(au >= 0, au + cal, -1), np.where(ad >= 0, ad + cal, -1))
                for m, (au, ad) in res.items():
                    a_side, a_any = (au, ad)[side], first_alarm(au, ad)
                    early = (a_any >= 0) & (a_any < nu)
                    det = (a_side >= nu) & (a_side < nu + 2000) & ~early
                    d = a_side[det] - nu + 1
                    rows_delay.append(dict(law=law, change_after_restart=nu, sd_ratio=round(ratio, 3), method=m,
                                           alarm_before_change=float(early.mean()), detected_within_2000=float(det.mean()),
                                           median_delay=float(np.median(d)) if d.size else np.nan,
                                           q90_delay=float(np.quantile(d, 0.9)) if d.size else np.nan, paths=reps))
            print(f"  14f delay {law} nu {nu}: {time.time() - t0:6.1f} s", flush=True)
    a = pd.DataFrame(rows_null)
    b = pd.DataFrame(rows_delay)
    sfx = "_smoke" if SMOKE else ""
    if len(a):
        a.to_csv(OUT / f"e14f_awake_null{sfx}.csv", index=False)
        print(a.pivot_table(index=["law", "horizon"], columns="method", values="false_alarm_rate").round(4).to_string())
    if len(b):
        b.to_csv(OUT / f"e14f_awake_delay{sfx}.csv", index=False)
        print(b.pivot_table(index=["law", "change_after_restart", "sd_ratio"], columns="method", values="median_delay").round(0).to_string())
    return a, b


# ----------------------------------------------------------------------------
# Streaming monitors for real data
# ----------------------------------------------------------------------------


def make_monitor(method: str, band: float, cal: int, a_plug: float, a_delta: float, seed: int = 0):
    if method == "plugin":
        return VolatilityMonitor(threshold=a_plug, c=C, cal_len=cal, null_band=band)
    if method == "cs":
        return VolatilityMonitor(threshold=a_delta, c=C, cal_len=cal, null_band=band, baseline="cs", delta=DELTA)
    if method in ("A-grid-run", "A-grid-fixed", "A-dkw-run"):
        return RT.RunningTauMonitor(threshold=a_delta, c=C, cal_len=cal, method="dkw" if "dkw" in method else "grid",
                                    delta=DELTA, null_band=band, running=method != "A-grid-fixed")
    if method in ("B-lr", "B-ind"):
        return RT.ConformalMonitor(threshold=a_plug, cal_len=cal, family=method[2:], band=band, seed=seed)
    if method == "awake-A":
        return RT.AwakeRunningTauMonitor(threshold=a_delta, c=C, delta=DELTA, null_band=band)
    if method == "awake-B":
        return RT.ConformalMonitor(threshold=a_plug, cal_len=RT.AWAKE_PRE, family="lr", band=band, seed=seed, restart_len=0)
    raise ValueError(method)


def stream(mon, u):
    """Run a monitor; return its events and the per-step flag 'the return at t is monitored'."""
    awake = np.zeros(u.size, bool)
    for t in range(u.size):
        awake[t] = mon.calibrated
        mon.update(u[t])
    return mon.events, awake


def robust_sd(x, min_n):
    x = x[np.isfinite(x)]
    return float(np.median(np.abs(x)) * MAD_NORMAL_SCALE) if x.size >= min_n else np.nan


REAL_RUNS = [(m, b) for m in ("plugin", "cs", "A-grid-run", "B-lr", "awake-A", "awake-B") for b in (0.0, 0.25)] + \
            [(m, 0.25) for m in ("A-grid-fixed", "A-dkw-run", "B-ind")]


def real_job(job):
    source, method, band = job
    if source == "Nifty 50":
        from e07_nifty import deseasonalise, load_clean

        df, _, _ = load_clean("2023-11-01", "2026-01-01")
        u_all = deseasonalise(df).to_numpy()
        keep = np.asarray(df.index >= pd.Timestamp("2024-01-01", tz="Asia/Kolkata"))
        u = u_all[keep]
        cal, day, a_plug = 750, 375, 7500.0
    else:
        from avmonitor.data import intraday_scale, load_binance

        df = load_binance(source, "2025-06-01", "2026-01-01", download=False)
        s = intraday_scale(df["logret"], lookback_days=20)
        uu = (df["logret"] / s).to_numpy(copy=True)
        uu[~np.isfinite(uu)] = np.nan
        keep = np.asarray(df.index >= pd.Timestamp("2025-07-01", tz="UTC"))
        u = uu[keep]
        cal, day, a_plug = 1440, 1440, 28800.0
    a_delta = day / (0.05 - DELTA)
    mon = make_monitor(method, band, cal, a_plug, a_delta, seed=SEED)
    t0 = time.time()
    events, awake = stream(mon, u)
    secs = time.time() - t0
    n_days = u.size / day
    min_n = int(0.7 * day)
    rec = dict(source=source, method=method, band=band, A=a_plug if method in ("plugin", "B-lr", "B-ind", "awake-B") else a_delta,
               steps=u.size, days=n_days, alarms=len(events), alarms_per_day=len(events) / n_days,
               up=sum(e.direction == "up" for e in events), down=sum(e.direction == "down" for e in events),
               awake_share=float(awake.mean()), seconds=secs, us_per_step=secs / u.size * 1e6)
    ev_rows = []
    for e in events:
        post = robust_sd(u[e.index + 1 : e.index + 1 + day], min_n) / e.baseline_scale if e.baseline_scale > 0 else np.nan
        ev_rows.append(dict(source=source, method=method, band=band, index=e.index, direction=e.direction,
                            delay=e.index - e.monitored_since, post_over_baseline=post))
    ev = pd.DataFrame(ev_rows, columns=["source", "method", "band", "index", "direction", "delay", "post_over_baseline"])
    for d, cmp in (("up", np.greater), ("down", np.less)):
        p = ev.loc[ev.direction == d, "post_over_baseline"].dropna().to_numpy()
        k = int(np.sum(cmp(p, 1.0)))
        lo, hi = wilson(k, p.size)
        rec.update({f"{d}_checked": p.size, f"{d}_hit_rate": k / p.size if p.size else np.nan, f"{d}_hit_lo": lo,
                    f"{d}_hit_hi": hi, f"{d}_median_post_over_baseline": float(np.median(p)) if p.size else np.nan})
    return rec, ev


def q14d():
    from multiprocessing import Pool

    sources = ["BTCUSDT", "ETHUSDT", "XRPUSDT", "SOLUSDT", "DOGEUSDT", "Nifty 50"]
    runs = REAL_RUNS if not SMOKE else [("plugin", 0.25), ("awake-A", 0.25), ("B-lr", 0.25)]
    if SMOKE:
        sources = ["BTCUSDT"]
    jobs = [(s, m, b) for s in sources for m, b in runs]
    with Pool(2) as pool:
        res = pool.map(real_job, jobs, chunksize=1)
    summ = pd.DataFrame([r for r, _ in res])
    ev = pd.concat([e for _, e in res], ignore_index=True)
    sfx = "_smoke" if SMOKE else ""
    summ.to_csv(OUT / f"e14d_real_summary{sfx}.csv", index=False)
    ev.to_csv(OUT / f"e14d_real_alarms{sfx}.csv", index=False)
    cols = ["source", "method", "band", "alarms_per_day", "up", "down", "up_hit_rate", "down_hit_rate", "awake_share",
            "us_per_step"]
    print(summ[cols].round(3).to_string(index=False))
    return summ


# ----------------------------------------------------------------------------
# Q14g: E12 rerun with monitors that are awake after alarms
# ----------------------------------------------------------------------------

G_METHODS = ("plugin", "awake-A", "awake-B", "B-lr")


def g_job(job):
    import e12_crypto_events as E12

    coin, method = job
    events, _, _ = E12.read_events()
    df, u = E12.load_coin(E12.COINS[coin])
    index = df.index
    events["i0"] = index.get_indexer(events["release_utc"])
    mon = make_monitor(method, E12.BAND, E12.CAL, E12.A, E12.DAY / (0.05 - DELTA), seed=SEED)
    evs, awake = stream(mon, u)
    ev = E12.alarm_frame(evs, 0, u)
    event_days = events["date"].to_list()
    cache = {}
    rows = []
    for _, row in events.iterrows():
        key = (row["date"].weekday(), row["release_time_et"])
        if key not in cache:
            cache[key] = E12.control_anchors(index, event_days, key[0], key[1], u.size)
        ctrl = [E12.window_stats(ev, awake, u, i) for i in cache[key]]
        st = E12.window_stats(ev, awake, u, int(row["i0"]))
        rows.append(dict(coin=coin, method=method, date=row["date"].date(), event_type=row["event_type"],
                         awake_at_release=st["monitoring_at_start"], up_within_6h=st["up_post"],
                         first_up_delay_min=st["first_up_delay_min"], n_controls=len(ctrl),
                         control_up_rate=float(np.mean([c["up_post"] for c in ctrl])),
                         control_awake_rate=float(np.mean([c["monitoring_at_start"] for c in ctrl]))))
    summ = dict(coin=coin, method=method, days=u.size / 1440, alarms=len(ev), alarms_per_day=len(ev) / (u.size / 1440),
                up=int((ev.direction == "up").sum()), down=int((ev.direction == "down").sum()),
                awake_share=float(awake.mean()))
    return pd.DataFrame(rows), summ


def q14g():
    from multiprocessing import Pool

    import e12_crypto_events as E12

    coins = ["BTC"] if SMOKE else list(E12.COINS)
    jobs = [(c, m) for c in coins for m in G_METHODS]
    with Pool(2) as pool:
        res = pool.map(g_job, jobs, chunksize=1)
    win = pd.concat([w for w, _ in res], ignore_index=True)
    summ = pd.DataFrame([s for _, s in res])
    tests = []
    for (coin, method, etype), g in win.groupby(["coin", "method", "event_type"]):
        ps = g["control_up_rate"].to_numpy()  # stratum probability (event label exchangeable within its stratum)
        n_c = g["n_controls"].to_numpy()
        hits = g["up_within_6h"].to_numpy().astype(int)
        p_strat = (ps * n_c + hits) / (n_c + 1)
        T = int(hits.sum())
        tests.append(dict(coin=coin, method=method, event_type=etype, releases=len(g),
                          awake_at_release=int(g["awake_at_release"].sum()), up_within_6h=T,
                          control_up_rate=float(np.average(ps, weights=n_c)),
                          control_awake_rate=float(np.average(g["control_awake_rate"], weights=n_c)),
                          stratified_p=E12.poisson_binomial_tail(p_strat, T),
                          median_first_up_delay_min=float(g["first_up_delay_min"].median())))
    tests = pd.DataFrame(tests)
    sfx = "_smoke" if SMOKE else ""
    win.to_csv(OUT / f"e14g_e12_windows{sfx}.csv", index=False)
    summ.to_csv(OUT / f"e14g_e12_summary{sfx}.csv", index=False)
    tests.to_csv(OUT / f"e14g_e12_tests{sfx}.csv", index=False)
    print(summ.round(3).to_string(index=False))
    print(tests.round(4).to_string(index=False))
    return tests


# ----------------------------------------------------------------------------
# Q14e: cost per step
# ----------------------------------------------------------------------------


def q14e():
    import platform

    rows = []
    n = 20_000 if SMOKE else 200_000
    rng = np.random.default_rng([SEED, 12])
    r = rng.standard_normal(n)
    for name, mk in [
        ("VolatilityMonitor plugin", lambda: VolatilityMonitor(threshold=1e12, cal_len=1440)),
        ("VolatilityMonitor cs", lambda: VolatilityMonitor(threshold=1e12, cal_len=1440, baseline="cs")),
        ("RunningTauMonitor grid (route A)", lambda: RT.RunningTauMonitor(threshold=1e12, cal_len=1440)),
        ("RunningTauMonitor dkw (route A)", lambda: RT.RunningTauMonitor(threshold=1e12, cal_len=1440, method="dkw")),
        ("ConformalMonitor lr (route B)", lambda: RT.ConformalMonitor(threshold=1e12, cal_len=1440)),
        ("AwakeRunningTauMonitor (route A, no blind period)", lambda: RT.AwakeRunningTauMonitor(threshold=1e12)),
    ]:
        mon = mk()
        t = time.perf_counter()
        mon.run(r)
        sec = time.perf_counter() - t
        rows.append(dict(task=name + ", one series, no alarm", steps=n, seconds=sec, us_per_step=sec / n * 1e6))
    # learning-phase cost: restart every 1440 steps by forcing alarms with a tiny threshold is not possible
    # without changing the statistic, so time the awake monitor on its first 1440 steps only
    mon = RT.AwakeRunningTauMonitor(threshold=1e12)
    t = time.perf_counter()
    mon.run(r[: RT.AWAKE_PRE + RT.AWAKE_LEARN])
    sec = time.perf_counter() - t
    rows.append(dict(task="AwakeRunningTauMonitor, learning phase (first 1560 steps)", steps=RT.AWAKE_PRE + RT.AWAKE_LEARN,
                     seconds=sec, us_per_step=sec / (RT.AWAKE_PRE + RT.AWAKE_LEARN) * 1e6))
    reps, cal, m = 1000, 1440, 1950
    x = rng.standard_normal((reps, cal + m))
    for meth in ("plugin", "cs", "A-grid-run", "A-dkw-run", "B-lr"):
        t = time.perf_counter()
        run_method(meth, x[:, :cal], x[:, cal:], np.random.default_rng(0))
        sec = time.perf_counter() - t
        rows.append(dict(task=f"vectorised {meth} (calibration + monitoring), {reps} paths x {cal + m} steps",
                         steps=reps * (cal + m), seconds=sec, us_per_step=sec / (reps * (cal + m)) * 1e6))
    t = time.perf_counter()
    RT.awake_route_a(x[:, :120], x[:, 120:], [2 * A_DELTA])
    sec = time.perf_counter() - t
    rows.append(dict(task=f"vectorised awake-A, {reps} paths x {cal + m - 120} steps", steps=reps * (cal + m - 120),
                     seconds=sec, us_per_step=sec / (reps * (cal + m - 120)) * 1e6))
    df = pd.DataFrame(rows)
    df["machine"] = f"{platform.machine()}, Python {platform.python_version()}, 4 cores shared with other jobs"
    df.to_csv(OUT / f"e14e_cost{'_smoke' if SMOKE else ''}.csv", index=False)
    print(df[["task", "steps", "seconds", "us_per_step"]].round(3).to_string(index=False))
    return df


# ----------------------------------------------------------------------------
# Q14w: band widths, constants and the dilution of route B (supporting numbers for the theory note)
# ----------------------------------------------------------------------------


def q14w():
    from scipy import integrate, stats

    rows = []
    reps = 200 if SMOKE else 2000
    for li, law in enumerate(LAWS):
        truth = Truth(law)
        for ci, cal in enumerate(CAL_LENS):
            r = draw(law, (reps, cal), np.random.default_rng([SEED, 13, li, ci]))
            for kind in ("grid", "dkw", "dkw-fixed-n"):
                tb = RT.calibrate_band(r, method="grid" if kind == "grid" else "dkw", c=C, delta=DELTA,
                                       time_uniform=kind != "dkw-fixed-n")
                mt = truth.m(tb.levels)
                for k in range(tb.levels.shape[1]):
                    rows.append(dict(law=law, cal_len=cal, band=kind, level_index=k + RT.GRID_KMIN,
                                     level_over_centre=RT.GRID_RATIO ** (k + RT.GRID_KMIN),
                                     median_upper_over_truth=float(np.median(tb.upper[:, k] / mt[:, k])),
                                     median_lower_over_truth=float(np.median(tb.lower[:, k] / mt[:, k])),
                                     failed_any_level=float(np.mean(np.any((mt > tb.upper) | (mt < tb.lower), axis=1))),
                                     eps=tb.info.get("eps", np.nan), paths=reps))
            b = cs_baseline(r, c=C, delta_upper=DELTA / 2, delta_lower=DELTA / 2)
            mt = truth.m(b.tau)
            rows.append(dict(law=law, cal_len=cal, band="E8 cs (one level, sample split)", level_index=0, level_over_centre=1.0,
                             median_upper_over_truth=float(np.median(b.upper / mt)),
                             median_lower_over_truth=float(np.median(b.lower / mt)),
                             failed_any_level=float(np.mean((mt > b.upper) | (mt < b.lower))), eps=np.nan, paths=reps))
    pd.DataFrame(rows).to_csv(OUT / f"e14w_bands{'_smoke' if SMOKE else ''}.csv", index=False)
    K = RT.GRID_KMAX - RT.GRID_KMIN + 1
    const = [dict(name="Howard-Ramdas C for alpha = 0.01 (A = 0.85)", value=RT.hr_dkw_constant(DELTA)),
             dict(name="closed form 0.8 log(1612 / 0.01)", value=0.8 * np.log(1612 / DELTA))]
    for n in CAL_LENS:
        const += [dict(name=f"Howard-Ramdas radius at t = t_min = {n}, alpha = 0.01", value=float(RT.hr_dkw_radius(n, DELTA, t_min=n))),
                  dict(name=f"Howard-Ramdas radius at t = {n}, t_min = 1, alpha = 0.01", value=float(RT.hr_dkw_radius(n, DELTA))),
                  dict(name=f"Massart fixed-n DKW radius, n = {n}, alpha = 0.01", value=RT.massart_dkw_radius(n, DELTA))]
    const += [dict(name=f"log(2K / delta), K = {K}, delta = 0.01", value=np.log(2 * K / DELTA)),
              dict(name="log(2 / delta), delta = 0.01", value=np.log(2 / DELTA))]
    for f in RT.bet_family("lr")[0] + RT.bet_family("lr")[1]:
        const.append(dict(name=f"rank LR bet rho = {f.rho:.3f}: largest factor f(u0)", value=float(f(f.u0))))
        const.append(dict(name=f"rank LR bet rho = {f.rho:.3f}: smallest factor", value=float(min(f(1.0), f(f.u0)))))
    pd.DataFrame(const).to_csv(OUT / "e14w_constants.csv", index=False)
    # dilution: expected log growth of the matching rank LR bet at post-change step k, reference of n0 pre-change values
    drows = []
    n0 = 1940
    lr = {round(f.rho, 3): f for f in RT.bet_family("lr")[0] + RT.bet_family("lr")[1]}
    for rho, f in lr.items():
        def g(k):
            dens = lambda a: 2 * stats.norm.pdf(a / rho) / rho
            pik = lambda a: (n0 * 2 * stats.norm.sf(a) + k * 2 * stats.norm.sf(a / rho)) / (n0 + k)
            val, _ = integrate.quad(lambda a: dens(a) * np.log(f(pik(a))[()] if np.ndim(f(pik(a))) else f(pik(a))), 0, 12 * rho,
                                    limit=200)
            return val
        ks = np.arange(0, 6001)
        gk = np.array([g(k) for k in ks[::50]])
        gfull = np.interp(ks, ks[::50], gk)
        cum = np.cumsum(gfull)
        target = np.log(2 * A_RANK)
        hit = int(np.argmax(cum >= target)) + 1 if np.any(cum >= target) else np.nan
        hit0 = target / gk[0] if gk[0] > 0 else np.nan
        for k in (0, 50, 200, 500, 1000, 1940, 5000):
            drows.append(dict(rho=rho, n0=n0, k=k, growth_per_step=float(np.interp(k, ks[::50], gk)),
                              steps_to_log_2A_with_dilution=hit, steps_to_log_2A_without_dilution=hit0))
    pd.DataFrame(drows).to_csv(OUT / "e14w_dilution.csv", index=False)
    print(pd.DataFrame(drows).round(4).to_string(index=False))
    print(pd.DataFrame(const).round(4).to_string(index=False))


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------


def figures():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = pd.concat([pd.read_csv(f) for f in sorted(OUT.glob("e14a_validity*.csv")) if "smoke" not in f.name])
    b = pd.concat([pd.read_csv(f) for f in sorted(OUT.glob("e14b_delay*.csv")) if "smoke" not in f.name])
    order = ["oracle", "plugin", "plugin-run", "cs", "A-grid-fixed", "A-grid-run", "A-dkw-run", "B-lr", "B-ind", "B-freeze"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, cal in zip(axes, CAL_LENS):
        g = a[(a.cal_len == cal) & (a.horizon == DAY)]
        xs = np.arange(len(order))
        for k, law in enumerate(LAWS):
            gg = g[g.law == law].set_index("method").loc[order]
            ax.errorbar(xs + (k - 1.5) * 0.15, 100 * gg.false_alarm_rate, yerr=[100 * (gg.false_alarm_rate - gg.ci_low),
                        100 * (gg.ci_high - gg.false_alarm_rate)], fmt="o", ms=3, capsize=2, label=law)
        ax.axhline(5, color="k", ls="--", lw=1)
        ax.text(len(order) - 0.6, 5.1, "proved bound 5%", ha="right", va="bottom", fontsize=8)
        ax.set_xticks(xs)
        ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
        ax.set_title(f"Q14a [synthetic]: one-day false alarms, calibration {cal}")
        ax.grid(axis="y", lw=0.4, color="#dddddd")
    axes[0].set_ylabel("false alarm rate within 390 steps (%)")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "e14_validity.png", dpi=130)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    meths = ["plugin", "cs", "A-grid-fixed", "A-grid-run", "A-dkw-run", "B-lr", "B-ind", "B-freeze"]
    for ax, cal in zip(axes, CAL_LENS):
        g = b[(b.cal_len == cal) & (b.law == "gaussian")]
        for m in meths:
            gg = g[g.method == m].sort_values("sd_ratio")
            ax.plot(gg.sd_ratio, gg.median_delay, marker="o", ms=3, lw=1, label=m)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks([0.5, 0.67, 0.8, 1.25, 1.5, 2.0])
        ax.set_xticklabels(["0.5", "0.67", "0.8", "1.25", "1.5", "2"])
        ax.set_xlabel("sd ratio after / before")
        ax.set_title(f"Q14b [synthetic, Gaussian]: median delay, calibration {cal}")
        ax.grid(lw=0.4, color="#dddddd")
    axes[0].set_ylabel("median delay (steps)")
    axes[0].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "e14_delay.png", dpi=130)
    plt.close(fig)


def main():
    parts = [p for p in sys.argv[1:] if not p.startswith("--")] or ["a", "b", "c", "f", "e"]
    FIG.mkdir(parents=True, exist_ok=True)
    for p in parts:
        t0 = time.time()
        {"a": q14a, "b": q14b, "c": q14c, "d": q14d, "e": q14e, "f": q14f, "g": q14g, "w": q14w, "figures": figures}[p]()
        print(f"part {p}: {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
