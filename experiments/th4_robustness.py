"""TH4 (tracker B84). The false alarm bound when the null holds only approximately, checked against E16.

All results are SYNTHETIC. Seeds are fixed (listed below). Theory: docs/theory_robustness.md. Bound functions:
src/avmonitor/robustness.py. The detector is the monitor's default: X_t = min(r_t^2, tau^2), increase bets 0.02
to 0.5, decrease bets 0.1 to 0.95 times 1 / (B - 1) with B = tau^2 / m0, Shiryaev-Roberts statistics, each side
alarming at 2A, A = 7800 (one day = 390 steps, five days = 1950 steps).

Parts
  sup        Theorem 1 in closed form: the inflation factor over m / A against the drift delta, the horizon at
             which the bound reaches 1, and the design drifts at which the one-day bound reaches 10% and 25%
             (th4_sup_inflation.csv, th4_design.csv); the drift at which the bound reaches 1 against the
             smallest drift at which the deterministic stream alarms (th4_sharpness.csv). No simulation.
  constant   Simulation check of Theorem 1 on streams with a known constant drift: E[Y_t | past] = 1 + delta
             (or 1 - delta) at every step, for a deterministic stream, the two-point law {0, B} and Gaussian
             returns with a shifted variance; and an adaptive adversary with the same drift (hold Y at its mean,
             then bold play as in E16; the hold is searched on 4,000 paths and the best is rerun on 20,000 fresh
             paths) (th4_constant_drift.csv).
  switching  Exogenous drift that averages out: the volatility switches between two levels whose truncated
             means are m0 (1 + eps) and m0 (1 - eps), by a Markov chain independent of the innovations.
             Sup form (Theorem 1) against the pathwise forms (Theorems 2 and 3) (th4_switching.csv).
  garch      The 64 GARCH(1,1) cells of E16 (same parameters, tau and m0 taken from results/e16b_garch.csv,
             fresh seeds): the measured rate on the same paths, E16's measured rate, the proved bounds for a
             drift that responds to the data (Theorem 2, Propositions 2b and 4) and two heuristics that are not
             proved for GARCH (the sup form at E16's stored drift statistics, and Theorem 3's pathwise form)
             (th4_garch.csv).
  plugin     The plug-in cells of E16: Corollary 1 evaluated on 100,000 fresh calibration windows per law and
             calibration length, against E16's measured plug-in rates (i.i.d. laws, the shape adversary and the
             rate by decile of the baseline error); and the calibration length the Gaussian plug-in would need
             (th4_plugin.csv, th4_plugin_deciles.csv, th4_plugin_caln.csv).
  ticks      E16's price-level ticks: Theorem 1 at the drift range E16 stored, which conditions on the previous
             return only, so the result is a heuristic, not a proved bound (th4_ticks.csv). No simulation.
  figures    results/figures/th4_sup_bound.png and th4_e16_comparison.png.

Usage (each part writes its own CSV files; run with at most 2 processes on a shared machine):
  python experiments/th4_robustness.py sup constant switching garch plugin ticks figures
  TH4_SCALE=0.05 TH4_OUT=/tmp/th4 python experiments/th4_robustness.py ...   # smoke run

Seeds (SEED = 18400): constant SEED + 1..70 (30 laws, then a search seed and a fresh seed for each of the 20
adaptive cells); switching SEED + 100..115; garch SEED + 200..263; plugin windows
SEED + 400..431; calibration length scan SEED + 500..505.
"""

from __future__ import annotations

import ast
import importlib.util
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize

from avmonitor import robustness as R

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("TH4_OUT", ROOT / "results"))
FIG = OUT / "figures"
E16_RESULTS = ROOT / "results"  # E16's stored results are always read from the repository
SCALE = float(os.environ.get("TH4_SCALE", "1"))
SEED = 18400
A = 7800.0
DAY, WEEK = 390, 1950
HORIZONS = (DAY, WEEK)
B_G = R.gaussian_B()
TAU_G = 3.0
M0_G = float(R.trunc_mean_normal(TAU_G))


def _n(P):
    return int(max(200, round(P * SCALE))) if SCALE != 1 else int(P)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(c - h, 0.0), min(c + h, 1.0)


def _e16():
    spec = importlib.util.spec_from_file_location("e16_adversarial_null", ROOT / "experiments" / "e16_adversarial_null.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rate(alarm, m):
    k = int(np.sum((alarm >= 0) & (alarm < m)))
    lo, hi = wilson(k, alarm.size)
    return k / alarm.size, lo, hi


def merge_runs(parts):
    """Concatenate the per-chunk outputs of R.run_paths."""
    out = dict(alarm=np.concatenate([p["alarm"] for p in parts]),
               alarm_corrected=np.concatenate([p["alarm_corrected"] for p in parts]), horizons={})
    for m in parts[0]["horizons"]:
        out["horizons"][m] = {k: np.concatenate([p["horizons"][m][k] for p in parts]) for k in parts[0]["horizons"][m]}
    if "ratio_max" in parts[0]:
        out["ratio_max"] = np.concatenate([p["ratio_max"] for p in parts])
    return out


def proved_bounds(res, m, B, side="two-sided"):
    """The bounds that hold for any predictable drift, from one run of R.run_paths at horizon m."""
    h = res["horizons"][m]
    thr, K = R.threshold_bound(h["dstar"], m, A)
    trs, a, b = R.truncated_sup_bound(h["max_up"], h["max_dn"], m, A, B)
    md = R.mean_drift_bound(h["sum_up"], h["sum_dn"], m, A, B=B)
    exo = float(np.mean(np.minimum(1.0, np.minimum(h["exo"], h["dstar"] * m / A))))
    return dict(bound_threshold=thr, threshold_K=K, bound_truncated_sup=trs, trunc_a=a, trunc_b=b, bound_mean_drift=md,
                bound_proved=min(thr, trs, md), exo_pathwise=exo)


# ----------------------------------------------------------------------------
# Part sup: Theorem 1 in closed form
# ----------------------------------------------------------------------------


def part_sup():
    rows = []
    deltas = [0.0, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1, 0.2]
    lam, phi = R.default_bets(B_G)
    for d in deltas:
        for direction in ("both", "up", "down"):
            du = d if direction in ("both", "up") else 0.0
            dd = d if direction in ("both", "down") else 0.0
            row = dict(delta=d, direction=direction)
            for m in (98, DAY, WEEK, 7800):
                row[f"inflation_m{m}"] = float(R.sup_inflation(du, dd, m, B_G))
                row[f"bound_m{m}"] = float(min(1.0, R.sup_bound(du, dd, m, A, B_G)))
            # horizon at which the bound first reaches 1
            mm = np.arange(1, 20001)
            b = R.sup_bound(du, dd, mm, A, B_G)
            row["m_bound_reaches_1"] = int(mm[np.argmax(b >= 1)]) if np.any(b >= 1) else np.nan
            # elementary bounds on G_m for the largest bet of each side (one day)
            x = max(lam.max() * du, phi.max() * dd)
            lo, hi = R.growth_factor_bounds(x, DAY)
            row.update(G_day_largest_bet=float(R.growth_factor(x, DAY)), G_day_lower=float(lo), G_day_upper=float(hi))
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th4_sup_inflation.csv", index=False)
    print(df[["delta", "direction", "inflation_m390", "inflation_m1950", "bound_m390", "bound_m1950",
              "m_bound_reaches_1"]].round(4).to_string())
    des = []
    for label, m, a_ in (("one day of 390 steps, A = 7800", DAY, A), ("one day of 1440 steps, A = 28800", 1440, 28800.0)):
        for target in (0.10, 0.25):
            for direction in ("both", "up", "down"):
                d = R.design_drift(target, m, a_, B_G, direction)
                des.append(dict(day=label, m=m, A=a_, B=B_G, target=target, direction=direction, delta=d,
                                deterministic_alarm_step=R.deterministic_alarm_step(1 + d if direction != "down" else 1 - d,
                                                                                    a_, B_G) if np.isfinite(d) else np.nan))
    des = pd.DataFrame(des)
    des.to_csv(OUT / "th4_design.csv", index=False)
    print(des.round(5).to_string())
    # sharpness (Proposition 2): the drift at which the bound reaches 1 against the smallest drift at which the
    # deterministic stream Y_t = 1 +/- delta alarms within m (bisection on the alarm step, which falls with delta)
    sh = []
    for m in HORIZONS:
        for direction in ("up", "down"):
            d_one = R.design_drift(1.0, m, A, B_G, direction)
            lo, hi = 0.0, 0.5
            for _ in range(50):
                mid = (lo + hi) / 2
                st = R.deterministic_alarm_step(1 + mid if direction == "up" else 1 - mid, A, B_G)
                lo, hi = (lo, mid) if 0 < st <= m else (mid, hi)
            sh.append(dict(horizon=m, direction=direction, delta_bound_reaches_1=d_one, delta_deterministic_alarms=hi))
    sh = pd.DataFrame(sh)
    sh.to_csv(OUT / "th4_sharpness.csv", index=False)
    print(sh.round(5).to_string())
    return df, des, sh


# ----------------------------------------------------------------------------
# Part constant: simulation check of Theorem 1
# ----------------------------------------------------------------------------


def gaussian_scale_for(ratio):
    """Standard deviation s with E[min((s Z)^2, 9)] = ratio x m0 (tau = 3, m0 = E[min(Z^2, 9)])."""
    return optimize.brentq(lambda s: s * s * float(R.trunc_mean_normal(TAU_G / s)) - ratio * M0_G, 0.05, 20.0)


def constant_law_y(law, mean, P, n, rng):
    if law == "deterministic":
        return np.full((P, n), mean)
    if law == "two-point {0, B}":
        return np.where(rng.random((P, n)) < mean / B_G, B_G, 0.0)
    if law == "Gaussian, variance shifted":
        s = gaussian_scale_for(mean)
        return np.minimum((s * rng.standard_normal((P, n))) ** 2, TAU_G**2) / M0_G
    raise ValueError(law)


def part_constant():
    rows = []
    k = 0
    cells = [("up", d) for d in (0.0025, 0.005, 0.01, 0.02, 0.03, 0.05)] + [("down", d) for d in (0.01, 0.02, 0.05, 0.1)]
    for law in ("deterministic", "two-point {0, B}", "Gaussian, variance shifted"):
        for direction, d in cells:
            k += 1
            rng = np.random.default_rng(SEED + k)
            mean = 1 + d if direction == "up" else 1 - d
            P = 1 if law == "deterministic" else _n(10000)
            parts = []
            for s0 in range(0, P, 2500):
                p = min(2500, P - s0)
                Y = constant_law_y(law, mean, p, WEEK, rng)
                parts.append(R.run_paths(Y, np.full(Y.shape, mean - 1), B_G, A, HORIZONS))
            res = merge_runs(parts)
            for m in HORIZONS:
                r, lo, hi = rate(res["alarm"], m)
                rc, _, _ = rate(res["alarm_corrected"], m)
                du, dd = (d, 0.0) if direction == "up" else (0.0, d)
                bound = float(min(1.0, R.sup_bound(du, dd, m, A, B_G)))
                rows.append(dict(law=law, direction=direction, delta=d, horizon=m, paths=P, rate=r, ci_low=lo, ci_high=hi,
                                 bound_sup=bound, bound_over_rate=bound / r if r > 0 else np.inf,
                                 above_bound=bool(lo > bound), corrected_rate=rc, m_over_A=m / A,
                                 deterministic_alarm_step=R.deterministic_alarm_step(mean, A, B_G)))
            log("constant", law, direction, d, [(round(rw["rate"], 4), round(rw["bound_sup"], 4)) for rw in rows[-2:]])
    # an adaptive adversary with the same drift: hold Y = mean for h steps (every bet gains without risk), then
    # bold play (spread {0, B} with mean 1 +/- delta, landing exactly on 2A when one draw can do it), as in E16
    E = _e16()
    for direction, d in cells:
        mean = 1 + d if direction == "up" else 1 - d
        for m in HORIZONS:
            k += 1
            holds = np.unique(np.round(np.linspace(0, 0.995, 15) * m).astype(int))
            rng = np.random.default_rng(SEED + k)
            U = rng.random((_n(4000), m))
            search = [float(np.mean(E.engine(drift_bold_fn(mean, U, h), U.shape[0], m, B_G)["alarm"] >= 0)) for h in holds]
            best = int(holds[int(np.argmax(search))])
            k += 1
            rng = np.random.default_rng(SEED + k)
            alarms = []
            for s0 in range(0, _n(20000), 5000):
                U = rng.random((min(5000, _n(20000) - s0), m))
                alarms.append(E.engine(drift_bold_fn(mean, U, best), U.shape[0], m, B_G)["alarm"])
            alarm = np.concatenate(alarms)
            r, lo, hi = rate(alarm, m)
            du, dd = (d, 0.0) if direction == "up" else (0.0, d)
            bound = float(min(1.0, R.sup_bound(du, dd, m, A, B_G)))
            rows.append(dict(law="hold, then bold play (best hold)", direction=direction, delta=d, horizon=m,
                             paths=alarm.size, rate=r, ci_low=lo, ci_high=hi, bound_sup=bound,
                             bound_over_rate=bound / r if r > 0 else np.inf, above_bound=bool(lo > bound),
                             corrected_rate=np.nan, m_over_A=m / A,
                             deterministic_alarm_step=R.deterministic_alarm_step(mean, A, B_G), hold=best,
                             search_rate=max(search)))
            log("constant hold-then-bold", direction, d, m, "hold", best, round(r, 4), "bound", round(bound, 4))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th4_constant_drift.csv", index=False)
    print(df[["law", "direction", "delta", "horizon", "rate", "ci_high", "bound_sup", "corrected_rate", "hold"]]
          .round(4).to_string())
    return df


def drift_bold_fn(mean, U, hold, B=B_G, thr=2 * A):
    """Adaptive two-point law with conditional mean ``mean``: Y = mean for ``hold`` steps, then bold play.

    Bold play: high point B, low point 0, except that a high point that takes the increase side exactly to 2A
    (or a low point that takes the decrease side there) is used when it is a legal value. For E16's engine.
    """
    from avmonitor.detectors import DEFAULT_UP_GRID

    up = np.asarray(DEFAULT_UP_GRID, float)
    Ut = np.ascontiguousarray(U.T)
    tgt = thr * (1 + 1e-9)

    def y_fn(t, idx, Su, Sd, lamd):
        if t < hold:
            return np.full(idx.size, mean)
        au, bu = Su.mean(axis=1) + 1.0, (Su + 1.0) @ up / up.size
        ad, bd = Sd.mean(axis=1) + 1.0, (Sd + 1.0) @ lamd / lamd.size
        ysu = 1.0 + (tgt - au) / bu
        yh = np.where(ysu <= B, np.maximum(ysu, mean), B)
        ysd = 1.0 - (tgt - ad) / bd
        yl = np.where(ysd >= 0.0, np.minimum(ysd, mean), 0.0)
        span = yh - yl
        p = np.where(span > 0, (mean - yl) / np.where(span > 0, span, 1.0), 0.0)
        return np.where(span > 0, np.where(Ut[t, idx] < p, yh, yl), mean)

    return y_fn


# ----------------------------------------------------------------------------
# Part switching: exogenous drift that averages out
# ----------------------------------------------------------------------------


def part_switching():
    rows = []
    k = 0
    for eps in (0.01, 0.02, 0.05, 0.1):
        s_hi, s_lo = gaussian_scale_for(1 + eps), gaussian_scale_for(1 - eps)
        for stay in (0.5, 0.9, 0.99, 0.999):
            k += 1
            rng = np.random.default_rng(SEED + 100 + k)
            P = _n(6000)
            parts = []
            for s0 in range(0, P, 2000):
                p = min(2000, P - s0)
                state = np.empty((p, WEEK), bool)
                state[:, 0] = rng.random(p) < 0.5
                flip = rng.random((p, WEEK)) > stay
                for t in range(1, WEEK):
                    state[:, t] = state[:, t - 1] ^ flip[:, t]
                sig = np.where(state, s_hi, s_lo)
                Y = np.minimum((sig * rng.standard_normal((p, WEEK))) ** 2, TAU_G**2) / M0_G
                d = np.where(state, eps, -eps)
                parts.append(R.run_paths(Y, d, B_G, A, HORIZONS, check=True))
            res = merge_runs(parts)
            for m in HORIZONS:
                r, lo, hi = rate(res["alarm"], m)
                rc, _, _ = rate(res["alarm_corrected"], m)
                pb = proved_bounds(res, m, B_G)
                sup = float(min(1.0, R.sup_bound(eps, eps, m, A, B_G)))
                rows.append(dict(eps=eps, stay=stay, half_life=np.log(0.5) / np.log(stay) if stay > 0.5 else 1.0,
                                 horizon=m, paths=P, rate=r, ci_low=lo, ci_high=hi, bound_sup=sup,
                                 bound_exogenous=pb["exo_pathwise"], bound_threshold=pb["bound_threshold"],
                                 bound_mean_drift=pb["bound_mean_drift"],
                                 bound_best=min(sup, pb["exo_pathwise"], pb["bound_threshold"], pb["bound_mean_drift"]),
                                 corrected_rate=rc, ratio_max=float(res["ratio_max"].max()), m_over_A=m / A))
            log("switching", eps, stay, [(round(rw["rate"], 4), round(rw["bound_sup"], 3), round(rw["bound_exogenous"], 3),
                                          round(rw["bound_threshold"], 3)) for rw in rows[-2:]])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th4_switching.csv", index=False)
    print(df[["eps", "stay", "horizon", "rate", "ci_high", "bound_sup", "bound_exogenous", "bound_threshold",
              "bound_mean_drift", "corrected_rate"]].round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Part garch: the 64 GARCH cells of E16
# ----------------------------------------------------------------------------


def part_garch():
    E = _e16()
    g16 = pd.read_csv(E16_RESULTS / "e16b_garch.csv")
    cells = g16.drop_duplicates(["persistence", "alpha"])
    rows = []
    for k, (_, c) in enumerate(cells.iterrows(), start=1):
        pers, alpha, beta, tau, m0 = c.persistence, c.alpha, c.beta, c.tau, c.m0
        B = tau**2 / m0
        burn = int(max(1000, 7 * np.log(2) / -np.log(pers)))
        rng = np.random.default_rng(SEED + 200 + k)
        P = _n(6000)
        parts = []
        for s0 in range(0, P, 2000):
            p = min(2000, P - s0)
            Rr, S = E.garch(p, WEEK, alpha, beta, rng, burn)
            Y = np.minimum(Rr * Rr, tau**2) / m0
            d = S * R.trunc_mean_normal(tau / np.sqrt(S)) / m0 - 1  # exact E[Y_t | past] - 1 (Gaussian innovations)
            del Rr, S
            parts.append(R.run_paths(Y, d, B, A, HORIZONS, check=True))
        res = merge_runs(parts)
        for m in HORIZONS:
            e16 = g16[(g16.persistence == pers) & (g16.alpha == alpha) & (g16.horizon == m)].iloc[0]
            r, lo, hi = rate(res["alarm"], m)
            rc, _, _ = rate(res["alarm_corrected"], m)
            pb = proved_bounds(res, m, B)
            heur_mean = float(min(1.0, R.sup_bound(e16.drift_mean_abs, e16.drift_mean_abs, m, A, B)))
            heur_q99 = float(min(1.0, R.sup_bound(max(e16.drift_q99, 0), max(-e16.drift_q01, 0), m, A, B)))
            rows.append(dict(persistence=pers, alpha=alpha, half_life=c.half_life, tau=tau, m0=m0, B=B, horizon=m,
                             bound=m / A, drift_mean_abs=e16.drift_mean_abs, drift_q01=e16.drift_q01,
                             drift_q99=e16.drift_q99, paths=P, rate=r, ci_low=lo, ci_high=hi, e16_rate=e16.rate,
                             e16_ci_low=e16.ci_low, e16_ci_high=e16.ci_high, **pb,
                             heuristic_sup_at_mean_drift=heur_mean, heuristic_sup_at_q99=heur_q99,
                             corrected_rate=rc, e16_oracle_rate=e16.oracle_rate,
                             ratio_max=float(res["ratio_max"].max()),
                             max_up_median=float(np.median(res["horizons"][m]["max_up"])),
                             dstar_median=float(np.median(res["horizons"][m]["dstar"]))))
        log("garch", pers, alpha, [(round(rw["rate"], 4), round(rw["e16_rate"], 4), round(rw["bound_proved"], 4),
                                    round(rw["exo_pathwise"], 4), round(rw["heuristic_sup_at_mean_drift"], 4))
                                   for rw in rows[-2:]])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th4_garch.csv", index=False)
    print(df[["persistence", "alpha", "horizon", "drift_mean_abs", "rate", "e16_rate", "bound_threshold",
              "bound_truncated_sup", "bound_mean_drift", "exo_pathwise", "heuristic_sup_at_mean_drift",
              "corrected_rate", "ratio_max"]].round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Part plugin: Corollary 1 against E16's plug-in cells
# ----------------------------------------------------------------------------


def plugin_windows(E, kind, par, cal, W, seed, chunk=10000):
    """Plug-in drift d = m(tau-hat) / m-hat - 1, B-hat and E16's error m-hat / m(tau-hat) - 1 for W windows."""
    rng = np.random.default_rng(seed)
    ds, bs, errs = [], [], []
    for s0 in range(0, W, chunk):
        p = min(chunk, W - s0)
        r, mfun = E.law_returns(kind, par, rng, (p, cal))
        tau2_hat, m_hat = E.plugin_baseline(r)
        m_true = np.asarray(mfun(np.sqrt(tau2_hat)), float)
        ds.append(R.plugin_drift(m_true, m_hat))
        bs.append(tau2_hat / m_hat)
        errs.append(m_hat / m_true - 1)
    return np.concatenate(ds), np.concatenate(bs), np.concatenate(errs)


def part_plugin():
    E = _e16()
    p16 = pd.read_csv(E16_RESULTS / "e16b_plugin.csv")
    main = p16[~p16.law.str.startswith("gaussian calibration")]
    W = _n(100000)
    cache = {}
    k = 0
    rows = []

    def windows(kind, par, cal):
        nonlocal k
        key = (kind, tuple(par), cal)
        if key not in cache:
            k += 1
            cache[key] = plugin_windows(E, kind, par, cal, W, SEED + 400 + k)
            d = cache[key][0]
            log("plugin windows", key, "drift sd", round(float(d.std()), 4), "share d > 0", round(float(np.mean(d > 0)), 3))
        return cache[key]

    for _, row in main.iterrows():
        par = tuple(ast.literal_eval(row.params))
        d, Bh, err = windows(row.law, par, int(row.cal))
        m = int(row.horizon)
        b = R.plugin_bound(d, Bh, m, A)
        rows.append(dict(cal=int(row.cal), law=row.law, params=row.params, monitored="i.i.d. (same law)", horizon=m,
                         bound_nominal=m / A, plugin_rate=row.plugin_rate, plugin_ci_low=row.plugin_ci_low,
                         plugin_ci_high=row.plugin_ci_high, bound_corollary1=b, windows=d.size,
                         drift_sd=float(d.std()), err_sd=float(err.std()), e16_err_sd=row.err_sd,
                         share_drift_up=float(np.mean(d > 0))))
    sh = p16[p16.law == "gaussian calibration, shape adversary"]
    for _, row in sh.iterrows():
        d, Bh, err = windows("gaussian", (), int(row.cal))
        m = int(row.horizon)
        b = R.plugin_bound(d, Bh, m, A)
        rows.append(dict(cal=int(row.cal), law="gaussian", params="[]", monitored=f"shape adversary, {row.params}",
                         horizon=m, bound_nominal=m / A, plugin_rate=row.plugin_rate, plugin_ci_low=row.plugin_ci_low,
                         plugin_ci_high=row.plugin_ci_high, bound_corollary1=b, windows=d.size, drift_sd=float(d.std()),
                         err_sd=float(err.std()), e16_err_sd=row.err_sd, share_drift_up=float(np.mean(d > 0))))
    df = pd.DataFrame(rows)
    df["bound_over_rate"] = df.bound_corollary1 / df.plugin_rate
    df["rate_above_bound"] = df.plugin_ci_low > df.bound_corollary1
    df.to_csv(OUT / "th4_plugin.csv", index=False)
    print(df[["cal", "law", "params", "monitored", "horizon", "plugin_rate", "plugin_ci_high", "bound_corollary1",
              "drift_sd", "e16_err_sd"]].round(4).to_string())
    # rate by decile of the baseline error: the bound restricted to windows with the error in the same range
    dec = p16[p16.law == "gaussian calibration, rate by baseline error decile"]
    drows = []
    for _, row in dec.iterrows():
        d, Bh, err = windows("gaussian", (), int(row.cal))
        sel = (err >= row.err_low) & (err <= row.err_high)
        m = int(row.horizon)
        b = R.plugin_bound(d[sel], Bh[sel], m, A) if sel.any() else np.nan
        drows.append(dict(cal=int(row.cal), params=row.params, horizon=m, err_low=row.err_low, err_high=row.err_high,
                          drift_low=1 / (1 + row.err_high) - 1, drift_high=1 / (1 + row.err_low) - 1,
                          plugin_rate=row.plugin_rate, e16_paths=int(row.paths), bound_corollary1=b,
                          windows=int(sel.sum())))
    ddf = pd.DataFrame(drows)
    ddf.to_csv(OUT / "th4_plugin_deciles.csv", index=False)
    print(ddf.round(4).to_string())
    # the calibration length the Gaussian plug-in would need
    crows = []
    for j, n in enumerate((390, 1440, 2880, 5760, 11520, 23040), start=1):
        Wn = _n(20000)
        d, Bh, err = plugin_windows(E, "gaussian", (), n, Wn, SEED + 500 + j, chunk=max(500, int(4e6 // n)))
        for m in HORIZONS:
            crows.append(dict(cal=n, horizon=m, bound_nominal=m / A, bound_corollary1=R.plugin_bound(d, Bh, m, A),
                              drift_sd=float(d.std()), share_drift_above_design10=float(np.mean(np.abs(d) > 0.0163)),
                              windows=Wn))
        log("plugin cal length", n, [round(c["bound_corollary1"], 4) for c in crows[-2:]])
    cdf = pd.DataFrame(crows)
    cdf.to_csv(OUT / "th4_plugin_caln.csv", index=False)
    print(cdf.round(4).to_string())
    return df, ddf, cdf


# ----------------------------------------------------------------------------
# Part ticks: price-level ticks of E16 (heuristic only)
# ----------------------------------------------------------------------------


def part_ticks():
    """Theorem 1 at the drift E16 stored for price ticks, which conditions on the previous return only.

    The theorems need the drift given the whole past. E[X_t | previous return] averages the drift given the
    whole past over the rest of the past, so its range is inside the range of the true drift: the bound below
    is a heuristic (it can be too small), not a proved bound.
    """
    t16 = pd.read_csv(E16_RESULTS / "e16b_price_ticks.csv")
    rows = []
    for _, r in t16.iterrows():
        cm = np.array([r.cond_mean_after_zero, r.cond_mean_after_one_tick, r.cond_mean_after_two_or_more], float)
        cm = cm[np.isfinite(cm)]
        du, dd = float(max(cm.max() - 1, 0)), float(max(1 - cm.min(), 0))
        m = int(r.horizon)
        rows.append(dict(tick=r.tick, B=r.B, horizon=m, bound_nominal=m / A, drift_up_previous_return=du,
                         drift_down_previous_return=dd, heuristic_sup=float(min(1.0, R.sup_bound(du, dd, m, A, r.B))),
                         price_rate=r.price_rounding_rate, price_ci_low=r.price_ci_low, price_ci_high=r.price_ci_high,
                         control_rate=r.marginal_iid_rate))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th4_ticks.csv", index=False)
    print(df.round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

INK, MUTED, GRIDC = "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def _style(ax):
    ax.grid(True, color=GRIDC, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def figures():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    plain = FuncFormatter(lambda v, _: f"{v:g}")
    FIG.mkdir(exist_ok=True)
    floor = 1e-3

    # Figure 1: Theorem 1 against constant-drift streams (rows: horizon, columns: direction)
    con = pd.read_csv(OUT / "th4_constant_drift.csv")
    des = pd.read_csv(OUT / "th4_design.csv")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.2), sharey=True)
    laws = (("hold, then bold play (best hold)", SERIES[3], "^", "adaptive adversary: hold, then bold play"),
            ("two-point {0, B}", SERIES[0], "o", "two-point {0, B}, i.i.d."),
            ("Gaussian, variance shifted", SERIES[2], "s", "Gaussian, variance shifted, i.i.d."))
    for j, direction in enumerate(("up", "down")):
        dmax = 0.06 if direction == "up" else 0.12
        dgrid = np.linspace(0, dmax, 400)
        du, dd = (dgrid, 0 * dgrid) if direction == "up" else (0 * dgrid, dgrid)
        dfine = np.linspace(0.0005, dmax, 120)
        steps = np.array([R.deterministic_alarm_step(1 + v if direction == "up" else 1 - v, A, B_G) for v in dfine])
        for i, m in enumerate((DAY, WEEK)):
            ax = axes[i, j]
            ax.plot(dgrid, np.minimum(1, R.sup_bound(du, dd, m, A, B_G)), color=INK, lw=1.6,
                    label="Theorem 1 bound (proved)")
            det = np.where((steps > 0) & (steps <= m), 1.0, floor)
            ax.step(dfine, det, where="post", color=SERIES[1], lw=1.4, label="deterministic stream Y = 1 + delta")
            for law, col, mk, lab in laws:
                s_ = con[(con.law == law) & (con.direction == direction) & (con.horizon == m)].sort_values("delta")
                ax.plot(s_.delta, np.maximum(s_.rate, floor), color=col, lw=1.2, marker=mk, ms=4.5, label=lab)
            ax.axhline(m / A, color=MUTED, lw=0.8, ls="--", label="m / A (exact null, Theorem B)")
            if m == DAY:
                for tgt in (0.10, 0.25):
                    v = des[(des.m == DAY) & (des.target == tgt) & (des.direction == direction)].delta.iloc[0]
                    ax.axvline(v, color=MUTED, lw=0.8, ls=":")
                    ax.text(v, 0.0013, f" bound {tgt:.0%} at {v:.4f}", fontsize=7, color=MUTED, rotation=90,
                            va="bottom", ha="right")
            ax.set_yscale("log")
            ax.set_ylim(floor * 0.9, 1.3)
            ax.set_xlim(0, dmax)
            sign = "+" if direction == "up" else "-"
            ax.set_xlabel(f"drift delta (E[X_t | past] = m0 (1 {sign} delta))", color=INK, fontsize=9)
            if j == 0:
                ax.set_ylabel(f"P(false alarm within {m} steps), zeros at {floor}", color=INK, fontsize=9)
            ax.set_title(f"{'Upward' if direction == 'up' else 'Downward'} drift, {m} steps "
                         f"({'one day' if m == DAY else 'five days'})", fontsize=9.5, color=INK, loc="left")
            ax.yaxis.set_major_formatter(plain)
            _style(ax)
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle("TH4: the sup-form bound (Theorem 1) against streams with a known constant drift, default two-sided "
                 "detector (synthetic, A = 7800, B = 9.045)", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(FIG / "th4_sup_bound.png", dpi=140)
    plt.close(fig)

    # Figure 2: bound against measured rate, E16 cells and the switching cells
    g = pd.read_csv(OUT / "th4_garch.csv")
    pl = pd.read_csv(OUT / "th4_plugin.csv")
    sw = pd.read_csv(OUT / "th4_switching.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.0))
    lim = (floor * 0.8, 1.3)
    ax = axes[0]
    for m, mk in ((DAY, "o"), (WEEK, "s")):
        s_ = g[g.horizon == m]
        ax.scatter(np.maximum(s_.e16_rate, floor), s_.bound_proved, s=20, marker=mk, color=SERIES[0], edgecolor="white",
                   lw=0.5, label=f"best proved bound, {m} steps", zorder=3)
        ax.scatter(np.maximum(s_.e16_rate, floor), s_.heuristic_sup_at_mean_drift, s=20, marker=mk, color=SERIES[3],
                   edgecolor="white", lw=0.5, label=f"sup form at the mean drift (not proved), {m} steps", zorder=2)
    ax.set_title("GARCH cells of E16 (drift responds to the data)", fontsize=9.5, color=INK, loc="left")
    ax = axes[1]
    for cal, col in ((390, SERIES[0]), (1440, SERIES[2])):
        for m, mk in ((DAY, "o"), (WEEK, "s")):
            s_ = pl[(pl.cal == cal) & (pl.horizon == m)]
            iid = s_[s_.monitored.str.startswith("i.i.d.")]
            shp = s_[s_.monitored.str.startswith("shape")]
            ax.scatter(np.maximum(iid.plugin_rate, floor), iid.bound_corollary1, s=22, marker=mk, color=col,
                       edgecolor="white", lw=0.5, label=f"cal. {cal}, i.i.d., {m} steps", zorder=3)
            ax.scatter(np.maximum(shp.plugin_rate, floor), shp.bound_corollary1, s=22, marker=mk, facecolor="none",
                       edgecolor=col, lw=1.0, label=f"cal. {cal}, shape adversary, {m} steps", zorder=3)
    ax.set_title("Plug-in cells of E16 (Corollary 1)", fontsize=9.5, color=INK, loc="left")
    ax = axes[2]
    for m, mk in ((DAY, "o"), (WEEK, "s")):
        s_ = sw[sw.horizon == m]
        ax.scatter(np.maximum(s_.rate, floor), s_.bound_sup, s=20, marker=mk, color=SERIES[1], edgecolor="white", lw=0.5,
                   label=f"sup form (Theorem 1), {m} steps", zorder=3)
        ax.scatter(np.maximum(s_.rate, floor), np.minimum(s_.bound_exogenous, s_.bound_threshold), s=20, marker=mk,
                   color=SERIES[2], edgecolor="white", lw=0.5, label=f"pathwise (Theorems 2, 3), {m} steps", zorder=3)
    ax.set_title("Exogenous switching drift (new cells)", fontsize=9.5, color=INK, loc="left")
    for ax in axes:
        ax.plot(lim, lim, color=INK, lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        ax.set_xlabel(f"measured false alarm rate (zeros shown at {floor})", color=INK, fontsize=9)
        ax.set_ylabel("bound (capped at 1)", color=INK, fontsize=9)
        ax.xaxis.set_major_formatter(plain)
        ax.yaxis.set_major_formatter(plain)
        _style(ax)
        ax.legend(fontsize=6.5, frameon=False, loc="lower right")
    fig.suptitle("TH4: bound against measured rate (above the diagonal: the bound holds; synthetic, A = 7800)",
                 fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "th4_e16_comparison.png", dpi=140)
    plt.close(fig)
    log("figures written")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    parts = sys.argv[1:] or ["sup", "constant", "switching", "garch", "plugin", "ticks", "figures"]
    fns = dict(sup=part_sup, constant=part_constant, switching=part_switching, garch=part_garch, plugin=part_plugin,
               ticks=part_ticks, figures=figures)
    for p in parts:
        t0 = time.time()
        log("start", p)
        fns[p]()
        log("done", p, f"({time.time() - t0:.0f} s)")
