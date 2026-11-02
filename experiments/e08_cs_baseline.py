"""E8. A confidence-sequence baseline: does it keep the false alarm guarantee, and what does it cost?

Background: the e-detector guarantee (P(false alarm within m steps) <= m / A)
assumes the TRUE pre-change mean of X = min(r^2, tau^2). The monitor plugs in
the calibration-window mean, which breaks the guarantee slightly (E4d: 5.4%
against 5% with a 390-step window). docs/math.md, section 6, proves that
betting against one-sided time-uniform confidence bounds U (increase side)
and L (decrease side) gives P(any false alarm within m steps) <= delta + m / A,
with tau taken from the first part of the window (sample splitting).

Research questions
  Q8a [SYNTHETIC]. Without a change, what is the one-day (390 monitored steps)
       false alarm rate of the plug-in baseline, the plug-in plus one standard
       error (margin_z = 1) and the confidence-sequence (cs) baseline, for
       calibration lengths 390, 1440 and 5000, under Gaussian, jump and
       Student-t5 returns? A is chosen so that each method's bound is 5%
       (plug-in and margin: A = 7800; cs: delta = 0.01, A = 390 / 0.04 = 9750).
       How often does the cs bound itself fail (U below the true mean)?
  Q8b [SYNTHETIC]. What does validity cost in detection delay, for standard
       deviation changes x1.25, x1.5, x2 (increase) and x0.8, x0.67, x0.5
       (decrease) 500 steps after calibration?
  Q8c [SYNTHETIC]. Adaptive calibration length: stop calibrating once the
       relative width (U - L) / mean is at most w (w = 0.3 and 0.5; minimum 390,
       maximum 5000 steps). What is the distribution of the length, and are the
       false alarm rate and the delays as expected?
  Q8d [SYNTHETIC]. How far is the empirical Bernstein bound from the best
       possible betting bound for these distributions (oracle that knows the
       law and bets optimally) and from a normal-approximation interval?
  Q8e [REAL]. Bitstamp BTC/USD, July to December 2025, the configuration of
       E5c (1440-minute calibration, c = 3, bands 0, 0.25, 0.5): how do the
       alarm counts and hit rates change with the cs baseline
       (delta = 0.01, A = 1440 / 0.04 = 36000)?

  Q8g [SYNTHETIC]. How should delta be chosen? For a fixed total budget of 5%
       per day (A = 390 / (0.05 - delta)), how do the delays for sd x1.25 and
       x1.5 change with delta (Gaussian, 1440-step calibration)?
  Q8f [SYNTHETIC]. Volatility clustering (GARCH(1,1), persistence 0.99): does
       the cs bound still cover its target, the bet-weighted average of the
       conditional truncated variances in the calibration window (docs/math.md,
       Theorem 6.2), while an interval built from the i.i.d. standard error
       does not cover the plain average as often as it claims?

All settings were fixed before the runs; nothing was tuned on the reported data.
Seeds are fixed below. Run time: a few minutes on one core.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar

from avmonitor import robust as R
from avmonitor.confseq import eb_confseq, width_stopping_time
from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from avmonitor.monitor import VolatilityMonitor
from avmonitor.simulate import change_returns

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
DAY = 390
C = 3.0
TAU_FRAC = 0.25
A_PLUGIN = 7800.0  # 390 / 7800 = 5%
DELTA = 0.01  # cs failure probability per side in the one-sided analysis
A_CS = DAY / (0.05 - DELTA)  # 9750, so delta + m / A = 5%
SCEN = {
    "gaussian": dict(),
    "jumps": dict(p_jump=0.005, jump_scale=8.0),
    "student t5": dict(heavy_tail_df=5.0),
}
CAL_LENS = (390, 1440, 5000)
UP_RATIOS = (1.25, 1.5, 2.0)
DOWN_RATIOS = (0.8, 1 / 1.5, 0.5)
NU = 500  # change point, in monitored steps
N_MON = 3000  # monitored steps in the change experiment
REPS_NULL = 10000
REPS_CHANGE = 2000
CHUNK = 1000
NULL_CHUNK = 2500
UP_GRID = np.asarray(DEFAULT_UP_GRID)
DN_FRAC = np.asarray(DEFAULT_DOWN_FRACTIONS)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def true_mean_fn(scen: str, n: int = 4_000_000, seed: int = 999):
    """m(tau) = E[min(r^2, tau^2)] under the pre-change law, from a large sample (vectorised in tau)."""
    r = change_returns(1, n, rng=np.random.default_rng(seed), **SCEN[scen]).r[0]
    s = np.sort(r**2)
    cs = np.concatenate([[0.0], np.cumsum(s)])

    def m(tau):
        t2 = np.asarray(tau, float) ** 2
        k = np.searchsorted(s, t2)
        return (cs[k] + t2 * (n - k)) / n

    return m


def sr_first_alarms(x, m_up, m_dn, tau2, levels):
    """Shiryaev-Roberts betting statistics for both sides with per-path baselines.

    x: (reps, n) truncated squared returns; m_up, m_dn, tau2: (reps,). m_dn <= 0
    switches the decrease side off. Returns dict level -> (first up alarm,
    first down alarm), indices in 0..n-1 or -1.
    """
    reps, n = x.shape
    yu = x / m_up[:, None]
    down_on = (m_dn > 0) & (tau2 / np.where(m_dn > 0, m_dn, 1.0) > 1.0)
    b = np.where(down_on, tau2 / np.where(m_dn > 0, m_dn, 1.0), 2.0)
    lam_dn = DN_FRAC[None, :] / (b[:, None] - 1.0)
    yd = np.where(down_on[:, None], x / np.where(m_dn > 0, m_dn, 1.0)[:, None], 1.0)
    Su = np.zeros((reps, UP_GRID.size))
    Sd = np.zeros((reps, DN_FRAC.size))
    out = {L: [np.full(reps, -1), np.full(reps, -1)] for L in levels}
    top = max(levels)
    for t in range(n):
        with np.errstate(over="ignore"):
            Su = (1.0 + UP_GRID[None, :] * (yu[:, t : t + 1] - 1.0)) * (Su + 1.0)
            Sd = (1.0 + lam_dn * (1.0 - yd[:, t : t + 1])) * (Sd + 1.0)
        Mu = Su.mean(axis=1)
        Md = np.where(down_on, Sd.mean(axis=1), 0.0)
        for L in levels:
            au, ad = out[L]
            au[(au < 0) & (Mu >= L)] = t
            ad[(ad < 0) & (Md >= L)] = t
        if t % 50 == 0 and (out[top][0] >= 0).all() and (out[top][1] >= 0).all():
            break
    return out


def baselines(rcal, delta_up, delta_dn, stop=None, n_plan=None, n_tau=None):
    """Plug-in, margin_z = 1 and cs baselines for calibration windows rcal (reps, n).

    cs: tau from the first n_tau returns (default round(TAU_FRAC * n)), bounds from the rest,
    read at the end or at ``stop(cs)`` (number of sequence observations per row).
    """
    reps, n = rcal.shape
    scale = np.median(np.abs(rcal), axis=1) * R.MAD_NORMAL_SCALE
    tau = C * scale
    x = np.minimum(rcal**2, (tau**2)[:, None])
    m0 = x.mean(axis=1)
    se = np.sqrt(x.var(axis=1, ddof=1) / n)
    out = {
        "plugin": dict(tau=tau, up=m0, dn=m0),
        "margin_z=1": dict(tau=tau, up=m0 + se, dn=m0 - se),
    }
    n_tau = int(round(TAU_FRAC * n)) if n_tau is None else n_tau
    tau_s = C * np.median(np.abs(rcal[:, :n_tau]), axis=1) * R.MAD_NORMAL_SCALE
    xs = np.minimum(rcal[:, n_tau:] ** 2, (tau_s**2)[:, None])
    plan = (n - n_tau) if n_plan is None else n_plan
    U, L, M, K = (np.empty(reps) for _ in range(4))
    for a in range(0, reps, CHUNK):
        sl = slice(a, a + CHUNK)
        cs = eb_confseq(xs[sl], tau_s[sl] ** 2, delta_up, delta_dn, n_plan=plan)
        k = np.full(cs.upper.shape[0], cs.upper.shape[1]) if stop is None else stop(cs)
        i = np.arange(k.size)
        U[sl], L[sl], M[sl], K[sl] = cs.upper[i, k - 1], cs.lower[i, k - 1], cs.mean[i, k - 1], k
    out["cs"] = dict(tau=tau_s, up=U, dn=L, mean=M, n_cal=n_tau + K)
    return out


# ----------------------------------------------------------------------------
# Q8a and Q8b: fixed calibration lengths
# ----------------------------------------------------------------------------


def q8ab():
    null_rows, change_rows = [], []
    for si, scen in enumerate(SCEN):
        mtrue = true_mean_fn(scen)
        for ci, cal in enumerate(CAL_LENS):
            t0 = time.time()
            rng = np.random.default_rng(1000 + 10 * si + ci)
            # ---- null, in chunks to bound memory
            acc = {}
            for _ in range(REPS_NULL // NULL_CHUNK):
                rcal = change_returns(NULL_CHUNK, cal, rng=rng, **SCEN[scen]).r
                rmon = change_returns(NULL_CHUNK, DAY, rng=rng, **SCEN[scen]).r
                one = baselines(rcal, DELTA, DELTA)
                two = baselines(rcal, DELTA / 2, DELTA / 2)
                del rcal
                for meth in ("plugin", "margin_z=1", "cs"):
                    A = A_CS if meth == "cs" else A_PLUGIN
                    for label, b, lev in (("one", one[meth], A), ("two", two[meth] if meth == "cs" else one[meth], 2 * A)):
                        x = np.minimum(rmon**2, (b["tau"] ** 2)[:, None])
                        au, ad = sr_first_alarms(x, b["up"], b["dn"], b["tau"] ** 2, [lev])[lev]
                        mt = mtrue(b["tau"])
                        part = dict(au=au, ad=ad, fail_up=b["up"] < mt, fail_dn=b["dn"] > mt,
                                    ru=b["up"] / mt, rd=b["dn"] / mt, off=b["dn"] <= 0)
                        for k, v in part.items():
                            acc.setdefault((meth, label, k), []).append(v)
            cat = {k: np.concatenate(v) for k, v in acc.items()}
            for meth in ("plugin", "margin_z=1", "cs"):
                A = A_CS if meth == "cs" else A_PLUGIN
                g = {k[2]: v for k, v in cat.items() if k[0] == meth and k[1] == "one"}
                for side, a, fail, ratio in (("up", g["au"], g["fail_up"], g["ru"]), ("down", g["ad"], g["fail_dn"], g["rd"])):
                    fa = float((a >= 0).mean())
                    null_rows.append(dict(scenario=scen, cal_len=cal, method=meth, side=side, A=A,
                                          delta=DELTA if meth == "cs" else 0.0, bound=0.05,
                                          false_alarm_one_day=fa, mc_se=np.sqrt(fa * (1 - fa) / a.size),
                                          baseline_wrong_side=float(fail.mean()),
                                          median_null_over_true=float(np.median(ratio)),
                                          down_side_off=float(g["off"].mean()) if side == "down" else np.nan))
                # two-sided monitor semantics: each side at 2A; cs uses delta / 2 per side
                g = {k[2]: v for k, v in cat.items() if k[0] == meth and k[1] == "two"}
                fa = float(((g["au"] >= 0) | (g["ad"] >= 0)).mean())
                null_rows.append(dict(scenario=scen, cal_len=cal, method=meth, side="either (2A per side)", A=A,
                                      delta=DELTA if meth == "cs" else 0.0, bound=0.05, false_alarm_one_day=fa,
                                      mc_se=np.sqrt(fa * (1 - fa) / g["au"].size),
                                      baseline_wrong_side=float((g["fail_up"] | g["fail_dn"]).mean()),
                                      median_null_over_true=np.nan, down_side_off=np.nan))
            # ---- change
            rcal = change_returns(REPS_CHANGE, cal, rng=rng, **SCEN[scen]).r
            one = baselines(rcal, DELTA, DELTA)
            for ratio in UP_RATIOS + DOWN_RATIOS:
                rmon = change_returns(REPS_CHANGE, N_MON, nu=NU, ratio=ratio, rng=rng, **SCEN[scen]).r
                side = 0 if ratio > 1 else 1
                for meth in ("plugin", "margin_z=1", "cs"):
                    b = one[meth]
                    A = A_CS if meth == "cs" else A_PLUGIN
                    x = np.minimum(rmon**2, (b["tau"] ** 2)[:, None])
                    a = sr_first_alarms(x, b["up"], b["dn"], b["tau"] ** 2, [A])[A][side]
                    early = (a >= 0) & (a < NU)
                    det = a >= NU
                    d = a[det] - NU + 1
                    change_rows.append(dict(scenario=scen, cal_len=cal, method=meth, sd_ratio=round(ratio, 3),
                                            side="up" if side == 0 else "down", A=A,
                                            false_alarm_before_change=float(early.mean()),
                                            detected_within_2500=float(det.mean()),
                                            median_delay=float(np.median(d)) if d.size else np.nan,
                                            mean_delay=float(d.mean()) if d.size else np.nan))
            print(f"  {scen:11s} cal {cal:5d}: {time.time() - t0:5.1f} s", flush=True)
    a = pd.DataFrame(null_rows)
    b = pd.DataFrame(change_rows)
    a.to_csv(OUT / "e08a_null.csv", index=False)
    b.to_csv(OUT / "e08b_change.csv", index=False)
    return a, b


# ----------------------------------------------------------------------------
# Q8c: adaptive calibration length
# ----------------------------------------------------------------------------

N_MIN, N_MAX, N_TAU_ADAPT, N_PLAN_ADAPT = 390, 5000, 195, 1440 - 195
TARGETS = (0.3, 0.5)
REPS_ADAPT = 5000


def q8c():
    rows, lens = [], []
    for si, scen in enumerate(SCEN):
        mtrue = true_mean_fn(scen)
        for wi, w in enumerate(TARGETS):
            rng = np.random.default_rng(3000 + 10 * si + wi)

            def stop(cs, w=w):
                return width_stopping_time(cs, w, n_min=N_MIN - N_TAU_ADAPT, n_max=N_MAX - N_TAU_ADAPT)

            parts = []
            for _ in range(REPS_ADAPT // CHUNK):
                rcal = change_returns(CHUNK, N_MAX, rng=rng, **SCEN[scen]).r
                parts.append(baselines(rcal, DELTA, DELTA, stop=stop, n_plan=N_PLAN_ADAPT, n_tau=N_TAU_ADAPT)["cs"])
            del rcal
            b = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
            n_cal = b["n_cal"].astype(int)
            lens.append(pd.DataFrame(dict(scenario=scen, target=w, n_cal=n_cal)))
            mt = mtrue(b["tau"])
            rmon = change_returns(REPS_ADAPT, DAY, rng=rng, **SCEN[scen]).r
            x = np.minimum(rmon**2, (b["tau"] ** 2)[:, None])
            al = sr_first_alarms(x, b["up"], b["dn"], b["tau"] ** 2, [A_CS])[A_CS]
            row = dict(scenario=scen, target_rel_width=w, n_min=N_MIN, n_max=N_MAX,
                       len_q10=np.quantile(n_cal, 0.1), len_median=np.median(n_cal), len_q90=np.quantile(n_cal, 0.9),
                       share_capped=float(np.mean(n_cal >= N_MAX)), share_at_min=float(np.mean(n_cal <= N_MIN)),
                       up_FA_one_day=float((al[0] >= 0).mean()), down_FA_one_day=float((al[1] >= 0).mean()),
                       up_bound_failed=float(np.mean(b["up"] < mt)), down_bound_failed=float(np.mean(b["dn"] > mt)),
                       median_U_over_true=float(np.median(b["up"] / mt)), median_L_over_true=float(np.median(b["dn"] / mt)))
            for ratio in (1.5, 1 / 1.5):
                rch = change_returns(REPS_ADAPT, N_MON, nu=NU, ratio=ratio, rng=rng, **SCEN[scen]).r
                x = np.minimum(rch**2, (b["tau"] ** 2)[:, None])
                side = 0 if ratio > 1 else 1
                a = sr_first_alarms(x, b["up"], b["dn"], b["tau"] ** 2, [A_CS])[A_CS][side]
                det = a >= NU
                key = "up_x1.5" if ratio > 1 else "down_x0.67"
                row[f"{key}_median_delay"] = float(np.median(a[det] - NU + 1)) if det.any() else np.nan
                row[f"{key}_detected"] = float(det.mean())
            rows.append(row)
            print(f"  adaptive {scen:11s} w={w}: median length {np.median(n_cal):.0f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e08c_adaptive.csv", index=False)
    L = pd.concat(lens)
    L.to_csv(OUT / "e08c_adaptive_lengths.csv", index=False)
    return df, L


# ----------------------------------------------------------------------------
# Q8d: how wide is the bound compared with the best possible?
# ----------------------------------------------------------------------------


def q8d(null_df=None):
    rows = []
    for scen in SCEN:
        r = change_returns(1, 400_000, rng=np.random.default_rng(77), **SCEN[scen]).r[0]
        s = np.median(np.abs(r)) * R.MAD_NORMAL_SCALE
        tau2 = (C * s) ** 2
        Z = np.minimum(r**2, tau2) / tau2
        mu, sd = Z.mean(), Z.std()
        for cal in CAL_LENS:
            n = cal - int(round(TAU_FRAC * cal))
            target = np.log(1 / DELTA) / n

            def growth(m, sgn):
                f = lambda lam: -np.mean(np.log1p(lam * sgn * (m - Z)))  # noqa: E731
                hi = 1 / (1 - m) if sgn > 0 else 1 / m
                return -minimize_scalar(f, bounds=(0, 0.999 * hi), method="bounded").fun - target

            up = brentq(lambda m: growth(m, 1), mu + 1e-9, 0.999)
            lo = brentq(lambda m: growth(m, -1), 1e-6, mu - 1e-9)
            # the empirical Bernstein sequence on fresh samples of the same law (median over 500 paths)
            rng = np.random.default_rng(78)
            zz = rng.choice(Z, (500, n))
            cs = eb_confseq(zz, 1.0, DELTA, DELTA, n_plan=n)
            rows.append(dict(scenario=scen, cal_len=cal, n_cs=n, mean_Z=mu, sd_over_mean=sd / mu, range_over_sd=1 / sd,
                             oracle_betting_U_over_m=up / mu, eb_median_U_over_m=float(np.median(cs.upper[:, -1])) / mu,
                             normal_approx_U_over_m=1 + 2.3263 * sd / np.sqrt(n) / mu,
                             oracle_betting_L_over_m=lo / mu, eb_median_L_over_m=float(np.median(cs.lower[:, -1])) / mu,
                             normal_approx_L_over_m=1 - 2.3263 * sd / np.sqrt(n) / mu))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e08d_width_oracle.csv", index=False)
    return df


# ----------------------------------------------------------------------------
# Q8e: real data, Bitstamp BTC/USD (configuration of E5c)
# ----------------------------------------------------------------------------


def q8e():
    from avmonitor.data import intraday_scale, load_bitstamp_csv

    raw = ROOT / "data/raw/crypto/btcusd_bitstamp_1min_latest.csv"
    df = load_bitstamp_csv(raw, "2025-07-01", "2026-01-01")
    r = df.logret
    # Missing returns stay missing (NaN): the profile medians and the monitor skip them. E5c filled them
    # with 0, which counts a missing minute as a zero return (independent review, M6); on this file only
    # the first return is missing, so the two choices give the same results here.
    s = intraday_scale(r, lookback_days=20)
    u = (r / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    stream = u
    A_e5, A_cs = 28800.0, 1440 / (0.05 - DELTA)
    configs = []
    for band in (0.0, 0.25, 0.5):
        configs.append(("plugin", band, dict(threshold=A_e5, null_band=band)))
        configs.append(("cs", band, dict(threshold=A_cs, null_band=band, baseline="cs", delta=DELTA)))
    configs.append(("margin_z=1", 0.25, dict(threshold=A_e5, null_band=0.25, margin_z=1.0)))
    # control: the plug-in with the larger A of the cs runs, to separate the effect of A from that of the baseline
    configs.append(("plugin, A = 36000", 0.25, dict(threshold=A_cs, null_band=0.25)))
    configs.append(("cs adaptive (w=0.3, max 5760)", 0.25, dict(threshold=A_cs, null_band=0.25, baseline="cs",
                                                                 delta=DELTA, target_rel_width=0.3, cal_max=5760)))
    results, events = [], []
    days = len(u) / 1440
    for name, band, kw in configs:
        mon = VolatilityMonitor(c=C, cal_len=1440, **kw)
        ev = mon.run(stream)
        post, dirs = [], []
        for e in ev:
            seg = u[e.index + 1 : e.index + 1441]
            seg = seg[np.isfinite(seg)]
            if seg.size < 1000:
                continue
            ratio = np.median(np.abs(seg)) * R.MAD_NORMAL_SCALE / e.baseline_scale
            post.append(ratio)
            dirs.append(e.direction)
            events.append(dict(baseline=name, band=band, time=df.index[e.index], direction=e.direction,
                               delay_minutes=e.index - e.monitored_since, post_over_baseline=ratio))
        post, dirs = np.array(post), np.array(dirs)
        segs = pd.DataFrame(mon.segments)
        up_m, dn_m = dirs == "up", dirs == "down"
        res = dict(baseline=name, band=band, A=kw["threshold"], alarms=len(ev), alarms_per_day=len(ev) / days,
                   up=int(sum(e.direction == "up" for e in ev)), down=int(sum(e.direction == "down" for e in ev)),
                   median_post_ratio_up=float(np.median(post[up_m])) if up_m.any() else np.nan,
                   median_post_ratio_down=float(np.median(post[dn_m])) if dn_m.any() else np.nan,
                   share_up_followed_by_higher_vol=float(np.mean(post[up_m] > 1)) if up_m.any() else np.nan,
                   share_down_followed_by_lower_vol=float(np.mean(post[dn_m] < 1)) if dn_m.any() else np.nan,
                   median_cal_len=float(segs.cal_len.median()),
                   median_m_up_over_m0=float((segs.m_up / segs.m0).median()),
                   median_m_dn_over_m0=float((segs.m_dn / segs.m0).median()),
                   share_down_side_off=float(1 - segs.down_active.mean()),
                   monitored_share=float(1 - segs.cal_len.sum() / len(u)),
                   **{f"cal_{k}": v for k, v in mon.cal_stats.items()})
        results.append(res)
        print(f"  real {name:30s} band {band}: {len(ev)} alarms", flush=True)
    res = pd.DataFrame(results)
    res.to_csv(OUT / "e08e_real_summary.csv", index=False)
    pd.DataFrame(events).to_csv(OUT / "e08e_real_events.csv", index=False)
    return res


# ----------------------------------------------------------------------------
# Q8f: volatility clustering
# ----------------------------------------------------------------------------


def kappa2(a):
    """E[min(Z^2, a^2)] for standard normal Z (docs/math.md, section 3)."""
    from scipy import stats

    return (2 * stats.norm.cdf(a) - 1) - 2 * a * stats.norm.pdf(a) + 2 * a**2 * stats.norm.sf(a)


def q8f(reps=4000, cal=1440, alpha=0.05, beta=0.94, seed=4000):
    rng = np.random.default_rng(seed)
    omega = 1 - alpha - beta  # unconditional variance 1
    rows = []
    n = cal
    s2 = np.empty((reps, n))
    r = np.empty((reps, n))
    v = omega / (1 - alpha - beta) * np.ones(reps)
    # burn in 5000 steps so that the calibration window starts from the stationary law
    for _ in range(5000):
        e = rng.standard_normal(reps) * np.sqrt(v)
        v = omega + alpha * e**2 + beta * v
    for t in range(n):
        s2[:, t] = v
        r[:, t] = rng.standard_normal(reps) * np.sqrt(v)
        v = omega + alpha * r[:, t] ** 2 + beta * v
    n_tau = int(round(TAU_FRAC * cal))
    # cs: tau from the first part, target = bet-weighted average of mu_t over the second part
    tau_s = C * np.median(np.abs(r[:, :n_tau]), axis=1) * R.MAD_NORMAL_SCALE
    sig = np.sqrt(s2[:, n_tau:])
    mu = s2[:, n_tau:] * kappa2(tau_s[:, None] / sig)  # E[min(r_t^2, tau^2) | past]
    xs = np.minimum(r[:, n_tau:] ** 2, (tau_s**2)[:, None])
    for delta in (0.01, 0.05):
        cs = eb_confseq(xs, tau_s**2, delta, delta, n_plan=cal - n_tau)
        tu = (cs.lam_upper * mu).sum(axis=1) / cs.lam_upper.sum(axis=1)
        tl = (cs.lam_lower * mu).sum(axis=1) / cs.lam_lower.sum(axis=1)
        # time-uniform check over the whole sequence, and at the end of the window
        nu_u = np.cumsum(cs.lam_upper * mu, axis=1) / np.cumsum(cs.lam_upper, axis=1)
        nu_l = np.cumsum(cs.lam_lower * mu, axis=1) / np.cumsum(cs.lam_lower, axis=1)
        rows.append(dict(method="cs (empirical Bernstein)", nominal_miss=delta, target="bet-weighted average of mu_t",
                         miss_upper_at_end=float(np.mean(cs.upper[:, -1] <= tu)),
                         miss_lower_at_end=float(np.mean(cs.lower[:, -1] >= tl)),
                         miss_upper_any_t=float(np.mean((cs.upper <= nu_u).any(axis=1))),
                         miss_lower_any_t=float(np.mean((cs.lower >= nu_l).any(axis=1))),
                         median_upper_over_target=float(np.median(cs.upper[:, -1] / tu))))
    # i.i.d. standard error: tau from the whole window, target = plain average of mu_t over the window
    tau = C * np.median(np.abs(r), axis=1) * R.MAD_NORMAL_SCALE
    x = np.minimum(r**2, (tau**2)[:, None])
    target = (s2 * kappa2(tau[:, None] / np.sqrt(s2))).mean(axis=1)
    m0, se = x.mean(axis=1), np.sqrt(x.var(axis=1, ddof=1) / n)
    from scipy import stats

    for z, nominal in ((1.0, stats.norm.sf(1.0)), (stats.norm.isf(0.05), 0.05), (stats.norm.isf(0.01), 0.01)):
        rows.append(dict(method=f"plug-in +/- {z:.2f} iid SE", nominal_miss=float(nominal),
                         target="plain average of mu_t", miss_upper_at_end=float(np.mean(m0 + z * se <= target)),
                         miss_lower_at_end=float(np.mean(m0 - z * se >= target)), miss_upper_any_t=np.nan,
                         miss_lower_any_t=np.nan, median_upper_over_target=float(np.median((m0 + z * se) / target))))
    rows.append(dict(method="plug-in (no margin)", nominal_miss=np.nan, target="plain average of mu_t",
                     miss_upper_at_end=float(np.mean(m0 <= target)), miss_lower_at_end=float(np.mean(m0 >= target)),
                     miss_upper_any_t=np.nan, miss_lower_any_t=np.nan,
                     median_upper_over_target=float(np.median(m0 / target))))
    # how much the conditional truncated variance moves inside one window
    spread = np.quantile(mu, 0.9, axis=1) / np.quantile(mu, 0.1, axis=1)
    df = pd.DataFrame(rows)
    df["garch_alpha"], df["garch_beta"], df["cal_len"], df["reps"] = alpha, beta, cal, reps
    df["median_q90_over_q10_of_mu_in_window"] = float(np.median(spread))
    df.to_csv(OUT / "e08f_garch_coverage.csv", index=False)
    return df


# ----------------------------------------------------------------------------
# Q8g: choice of delta for a fixed budget
# ----------------------------------------------------------------------------


def q8g(reps=2000, cal=1440, seed=5000):
    rng = np.random.default_rng(seed)
    rcal = change_returns(reps, cal, rng=rng).r
    rnull = change_returns(reps, DAY, rng=rng).r
    rch = {ratio: change_returns(reps, N_MON, nu=NU, ratio=ratio, rng=rng).r for ratio in (1.25, 1.5)}
    rows = []
    for delta in (0.0025, 0.005, 0.01, 0.02, 0.03, 0.04):
        A = DAY / (0.05 - delta)
        b = baselines(rcal, delta, delta)["cs"]
        x = np.minimum(rnull**2, (b["tau"] ** 2)[:, None])
        au = sr_first_alarms(x, b["up"], b["dn"], b["tau"] ** 2, [A])[A][0]
        row = dict(delta=delta, A=A, cal_len=cal, median_U_over_mean=float(np.median(b["up"] / b["mean"])),
                   up_FA_one_day=float((au >= 0).mean()))
        for ratio, r in rch.items():
            x = np.minimum(r**2, (b["tau"] ** 2)[:, None])
            a = sr_first_alarms(x, b["up"], b["dn"], b["tau"] ** 2, [A])[A][0]
            det = a >= NU
            row[f"up_x{ratio}_median_delay"] = float(np.median(a[det] - NU + 1)) if det.any() else np.nan
            row[f"up_x{ratio}_detected"] = float(det.mean())
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e08g_delta_choice.csv", index=False)
    return df


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------


def figures(a, b, L):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cols = {"plugin": "#888888", "margin_z=1": "#4a6fa5", "cs": "#c0392b"}
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for k, scen in enumerate(SCEN):
        d = a[(a.scenario == scen) & (a.side == "up")]
        for meth, g in d.groupby("method"):
            ax[k].errorbar(g.cal_len, 100 * g.false_alarm_one_day, yerr=200 * g.mc_se, marker="o", label=meth,
                           color=cols[meth], capsize=3)
        ax[k].axhline(5, color="k", lw=0.8, ls="--")
        ax[k].set_xscale("log")
        ax[k].minorticks_off()
        ax[k].set_xticks(CAL_LENS, [str(c) for c in CAL_LENS])
        ax[k].set_title(f"{scen}: increase side, no change")
        ax[k].set_xlabel("calibration length (steps)")
    ax[0].set_ylabel("one-day false alarm rate (%)")
    ax[0].legend(frameon=False)
    fig.suptitle("[SYNTHETIC] False alarms with an estimated baseline (dashed: 5% bound; bars: 2 Monte Carlo SE)")
    fig.tight_layout()
    fig.savefig(FIG / "e08_false_alarms.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    d = b[(b.scenario == "gaussian")]
    for k, side in enumerate(("up", "down")):
        for meth in ("plugin", "margin_z=1", "cs"):
            for cal, ls in zip(CAL_LENS, (":", "--", "-")):
                g = d[(d.side == side) & (d.method == meth) & (d.cal_len == cal)].sort_values("sd_ratio")
                ax[k].plot(g.sd_ratio, g.median_delay, ls=ls, marker="o", color=cols[meth],
                           label=f"{meth}, cal {cal}" if k == 0 else None)
        ax[k].set_yscale("log")
        ax[k].minorticks_off()
        ax[k].set_yticks([20, 50, 100, 200, 500], ["20", "50", "100", "200", "500"])
        ax[k].set_xlabel("sd ratio after the change")
        ax[k].set_title(f"Gaussian, {'increase' if side == 'up' else 'decrease'}: median delay (steps)")
    ax[0].legend(frameon=False, fontsize=7, ncol=3)
    fig.suptitle("[SYNTHETIC] Detection delay by baseline and calibration length")
    fig.tight_layout()
    fig.savefig(FIG / "e08_delays.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
    bins = np.arange(N_MIN, L.n_cal.max() + 40, 20)
    for k, w in enumerate(TARGETS):
        for scen, col in zip(SCEN, ("#4a6fa5", "#c0392b", "#27ae60")):
            ax[k].hist(L[(L.scenario == scen) & (L.target == w)].n_cal, bins=bins, histtype="step", lw=1.5,
                       color=col, label=scen, density=True)
        ax[k].set_title(f"target relative width {w}")
        ax[k].set_xlabel("calibration length (steps)")
    ax[0].legend(frameon=False)
    fig.suptitle("[SYNTHETIC] Adaptive calibration length (allowed range 390 to 5000)")
    fig.tight_layout()
    fig.savefig(FIG / "e08_adaptive_length.png", dpi=130)
    plt.close(fig)


def main(parts="abcde"):
    t0 = time.time()
    if "a" in parts or "b" in parts:
        print("Q8a, Q8b [SYNTHETIC]")
        a, b = q8ab()
    else:
        a, b = pd.read_csv(OUT / "e08a_null.csv"), pd.read_csv(OUT / "e08b_change.csv")
    if "c" in parts:
        print("Q8c [SYNTHETIC]")
        c, L = q8c()
    else:
        c, L = pd.read_csv(OUT / "e08c_adaptive.csv"), pd.read_csv(OUT / "e08c_adaptive_lengths.csv")
    if "d" in parts:
        print("Q8d [SYNTHETIC]")
        d = q8d()
        print(d.round(3).to_string(index=False))
    if "e" in parts:
        print("Q8e [REAL]")
        e = q8e()
        print(e.round(3).to_string(index=False))
    if "g" in parts:
        print("Q8g [SYNTHETIC]")
        g = q8g()
        print(g.round(4).to_string(index=False))
    if "f" in parts:
        print("Q8f [SYNTHETIC]")
        f = q8f()
        print(f.round(4).to_string(index=False))
    figures(a, b, L)
    pd.set_option("display.width", 200)
    print(a.round(4).to_string(index=False))
    print(b.round(3).to_string(index=False))
    print(c.round(3).to_string(index=False))
    print(f"total {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "abcdefg")
