"""Experiments for the correlation monitor (alternative A7).

Research questions
  Q0 [SYNTHETIC] Theory table. For the sign statistic and the truncated product (c = 3 and
     c = 2), what are the null mean, the variance, the largest admissible bet, the best bet
     on the grid and the predicted detection delay log(2A) / growth rate for correlation
     changes such as 0.3 to 0.6, under Gaussian pairs? This says which statistic should win
     on clean data before any simulation is run.
  Q1 [SYNTHETIC] Validity. With a 2000-step calibration window from the same process, do the
     sign and product e-detectors keep P(false alarm within 390 steps) <= 390 / A (A = 7800,
     5% per day) and P(within 1950) <= 25% under Gaussian, Student t4, jumpy and co-jumpy
     pairs, and under a pure volatility change in one or both assets with the correlation
     fixed? The volatility change must NOT trigger the sign monitor. Benchmarks: a Gaussian
     likelihood ratio e-detector and a naive Page CUSUM on the raw product of standardised
     returns tuned to 5% per day under Gaussian pairs. A second table uses the true null
     mean (from 2 million draws) to check the guarantee without baseline error.
  Q2 [SYNTHETIC] Delay. Median detection delay and detection share for correlation changes
     (0.3 to 0.6, 0.6 to 0.3, 0.0 to 0.3, 0.3 to 0.0, 0.5 to 0.9, 0.6 to 0.0) at step 500,
     for Gaussian, Student t4 and jumpy innovations, all four detectors.
  Q3 [REAL] Binance BTC/ETH, BTC/SOL, BTC/DOGE, 1 January 2024 to 31 December 2025.
     (a) The Epps effect: Pearson correlation and sign-implied correlation at 1, 5, 15, 30 and
     60 minutes. (b) The streaming monitor on 5-minute returns (A = 5760, 5% per day,
     calibration 1440 bars = 5 days, bands 0, 0.1, 0.2): alarm counts, what the dependence does
     after each alarm against random times, whether alarms coincide with volatility changes
     (BTC robust scale after / before), and whether increase alarms cluster in BTC sell-offs
     (BTC return over the 24 hours before and after the alarm against random times).
  Q4 [SYNTHETIC and REAL] A basket of k = 4 assets. One e-detector on the average pairwise sign
     concordance (bounded in [-1/(k-1), 1]) against six pairwise monitors at level 6A: false
     alarms under Gaussian, t4 and one- or all-asset volatility changes with the correlations
     fixed; delay for basket-wide and single-pair correlation changes. Then the basket monitor on
     BTC, ETH, SOL and DOGE 5-minute returns, 2024 to 2025.

Settings follow E4 of the main report: per-step units, calibration from the same process as
the pre-change data, tau = c x scaled MAD per asset, A = 7800 for 390-step synthetic days and
A = 5760 for 288-bar real days (5-minute bars). All seeds are fixed. Runs in about 8 minutes
on one core (use ``--only q0,q1,q2``, ``--only q3`` or ``--only q4`` to run parts).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import corrmon as C  # noqa: E402
from simulate_pairs import pair_returns  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

A_SYN = 7800.0  # 5% false alarm budget per 390-step day
DAY = 390
CAL = 2000
NAIVE_DRIFT = 0.15

NULLS = {
    "gaussian rho 0.3": dict(rho=0.3),
    "student t4 rho 0.3": dict(rho=0.3, df=4.0),
    "jumps rho 0.3": dict(rho=0.3, p_jump=0.005, jump_scale=8.0),
    "cojumps rho 0.3": dict(rho=0.3, p_jump=0.005, p_cojump=0.002, jump_scale=8.0),
    "gaussian rho 0.6": dict(rho=0.6),
    "vol x2 asset 1 (rho 0.3 fixed)": dict(rho=0.3, nu=500, vol1_post=2.0),
    "vol x0.5 asset 1 (rho 0.3 fixed)": dict(rho=0.3, nu=500, vol1_post=0.5),
    "vol x2 both (rho 0.3 fixed)": dict(rho=0.3, nu=500, vol1_post=2.0, vol2_post=2.0),
    "garch(0.05, 0.94) each asset (rho 0.3 fixed)": dict(rho=0.3, garch=(0.05, 0.94)),
}
TRANSITIONS = [(0.3, 0.6), (0.6, 0.3), (0.0, 0.3), (0.3, 0.0), (0.5, 0.9), (0.6, 0.0)]
INNOVATIONS = {"gaussian": dict(), "student t4": dict(df=4.0), "jumps": dict(p_jump=0.005, jump_scale=8.0)}


# ----------------------------------------------------------------------------
# Shared machinery
# ----------------------------------------------------------------------------


def _split(reps, n, seed, nu=None, **kw):
    """Calibration window (CAL steps) and monitored paths from one process."""
    rng = np.random.default_rng(seed)
    p = pair_returns(reps, CAL + n, nu=None if nu is None else CAL + nu, rng=rng, **kw)
    return (p.r1[:, :CAL], p.r2[:, :CAL]), (p.r1[:, CAL:], p.r2[:, CAL:])


def _scales(cal):
    s1 = np.median(np.abs(cal[0]), axis=1) * C.MAD_NORMAL_SCALE
    s2 = np.median(np.abs(cal[1]), axis=1) * C.MAD_NORMAL_SCALE
    return s1[:, None], s2[:, None]


def _stat(stat, r, s1, s2, c):
    if stat == "sign":
        return C.sign_concordance(r[0], r[1])
    return np.clip(r[0] / s1, -c, c) * np.clip(r[1] / s2, -c, c)


def run_detector(name, cal, mon, A, band=0.0, true_m0=None, naive_h=None):
    """Two-sided first alarm per path for one detector. Each e-detector side runs at 2A."""
    s1, s2 = _scales(cal)
    if name in ("sign", "product c3", "product c2"):
        stat = "sign" if name == "sign" else "product"
        c = 2.0 if name == "product c2" else 3.0
        lo, hi = (-1.0, 1.0) if stat == "sign" else (-(c**2), c**2)
        y = _stat(stat, mon, s1, s2, c)
        m0 = np.mean(_stat(stat, cal, s1, s2, c), axis=1) if true_m0 is None else np.full(y.shape[0], true_m0)
        up = C.bounded_edetector(y, m0 + band, lo, hi, "up", threshold=2 * A)
        dn = C.bounded_edetector(y, m0 - band, lo, hi, "down", threshold=2 * A)
    elif name == "gaussian LR":
        m_sign = np.mean(C.sign_concordance(cal[0], cal[1]), axis=1)
        rho0 = C.rho_from_sign_mean(m_sign)
        z1, z2 = mon[0] / s1, mon[1] / s2
        up = C.gaussian_corr_lr_edetector(z1, z2, rho0, "up", threshold=2 * A)
        dn = C.gaussian_corr_lr_edetector(z1, z2, rho0, "down", threshold=2 * A)
    elif name == "naive CUSUM":
        y = (mon[0] / s1) * (mon[1] / s2)
        m0 = np.mean((cal[0] / s1) * (cal[1] / s2), axis=1)
        up = C.cusum_additive(y, m0, NAIVE_DRIFT, naive_h, "up")
        dn = C.cusum_additive(y, m0, NAIVE_DRIFT, naive_h, "down")
    else:
        raise ValueError(name)
    both = C.two_sided_alarm(up, dn)
    side = np.where(both < 0, "", np.where((up.alarm >= 0) & ((dn.alarm < 0) | (up.alarm <= dn.alarm)), "up", "down"))
    return both, side


def calibrate_naive(seed=11, reps=4000):
    """Threshold of the naive CUSUM: 2.5% per side within one 390-step day under Gaussian rho 0.3 pairs."""
    rng = np.random.default_rng(seed)
    p = pair_returns(reps, DAY, rho=0.3, rng=rng)
    return C.calibrate_cusum_threshold(p.r1 * p.r2, 0.3, NAIVE_DRIFT, DAY, 0.025)


# ----------------------------------------------------------------------------
# Q0: theory table
# ----------------------------------------------------------------------------


def q0_theory(A=A_SYN, n=400_000, seed=0):
    rows = []
    rng = np.random.default_rng(seed)
    fr = np.asarray(C.DEFAULT_FRACTIONS)
    for rho0, rho1 in TRANSITIONS:
        pre = pair_returns(1, n, rho=rho0, rng=rng)
        post = pair_returns(1, n, rho=rho1, rng=rng)
        for stat, c in (("sign", None), ("product", 3.0), ("product", 2.0)):
            if stat == "sign":
                y0 = C.sign_concordance(pre.r1[0], pre.r2[0])
                y1 = C.sign_concordance(post.r1[0], post.r2[0])
                lo, hi = -1.0, 1.0
            else:
                y0 = C.truncated_product(pre.r1[0], pre.r2[0], 1.0, 1.0, c)
                y1 = C.truncated_product(post.r1[0], post.r2[0], 1.0, 1.0, c)
                lo, hi = -(c**2), c**2
            m0, m1 = y0.mean(), y1.mean()
            direction = "up" if m1 > m0 else "down"
            lam_max = 1.0 / (m0 - lo) if direction == "up" else 1.0 / (hi - m0)
            g = (y1 - m0) if direction == "up" else (m0 - y1)
            growth_grid = [np.mean(np.log1p(f * lam_max * g)) for f in fr]
            k = int(np.argmax(growth_grid))
            lam_fine = np.linspace(0.001, 0.999, 999) * lam_max
            growth_fine = np.array([np.mean(np.log1p(l * g)) for l in lam_fine])
            kf = int(np.argmax(growth_fine))
            rows.append(dict(rho0=rho0, rho1=rho1, statistic=stat if c is None else f"{stat} c{c:g}",
                             direction=direction, m0=m0, m1=m1, shift=m1 - m0, sd_pre=y0.std(), sd_post=y1.std(),
                             snr=(m1 - m0) / y1.std(), lam_max=lam_max, unconstrained_opt_bet=(m1 - m0) / (y1.var() + (m1 - m0) ** 2) * (1 if direction == "up" else -1),
                             best_grid_fraction=fr[k], best_grid_bet=fr[k] * lam_max, growth_best_grid=growth_grid[k],
                             best_fine_bet=lam_fine[kf], growth_best_fine=growth_fine[kf],
                             predicted_delay_grid=np.log(2 * A) / max(growth_grid[k], 1e-12),
                             predicted_delay_fine=np.log(2 * A) / max(growth_fine[kf], 1e-12),
                             # Q2 changes at step 500, where the SR statistic is already about 500 under the
                             # null, and the best bet has weight 1/6 in the mixture: log(2A / 500) + log 6
                             predicted_delay_sr_from_500=(np.log(2 * A / 500.0) + np.log(len(fr))) / max(growth_grid[k], 1e-12)))
    df = pd.DataFrame(rows).drop(columns=[], errors="ignore")
    df.to_csv(OUT / "q0_theory.csv", index=False)
    return df


# ----------------------------------------------------------------------------
# Q1: validity
# ----------------------------------------------------------------------------


def true_null_means(seed=5, n=2_000_000):
    out = {}
    for name, kw in NULLS.items():
        if "nu" in kw or "garch" in kw:  # (a GARCH path of 2 million steps would be slow; its sign mean is the Gaussian one)
            continue
        rng = np.random.default_rng(seed)
        p = pair_returns(1, n, rng=rng, **kw)
        s1 = np.median(np.abs(p.r1[0])) * C.MAD_NORMAL_SCALE
        s2 = np.median(np.abs(p.r2[0])) * C.MAD_NORMAL_SCALE
        out[name] = dict(sign=float(C.sign_concordance(p.r1[0], p.r2[0]).mean()),
                         product=float(C.truncated_product(p.r1[0], p.r2[0], s1, s2, 3.0).mean()), scale1=s1, scale2=s2)
    return out


def q1_validity(reps=2000, horizon=5 * DAY, seed=100, naive_h=None):
    naive_h = calibrate_naive() if naive_h is None else naive_h
    rows = []
    for i, (name, kw) in enumerate(NULLS.items()):
        cal, mon = _split(reps, horizon, seed + i, **kw)
        for det in ("sign", "product c3", "gaussian LR", "naive CUSUM"):
            bands = (0.0, 0.05) if det in ("sign", "product c3") else (0.0,)
            for band in bands:
                alarm, side = run_detector(det, cal, mon, A_SYN, band=band, naive_h=naive_h)
                fa = alarm >= 0
                rows.append(dict(scenario=name, detector=det, band=band, baseline="estimated (2000 steps)",
                                 fa_1day=float((fa & (alarm < DAY)).mean()), bound_1day=DAY / A_SYN,
                                 fa_5days=float(fa.mean()), bound_5days=horizon / A_SYN,
                                 share_up=float((side == "up").sum() / max(fa.sum(), 1)), paths=reps))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q1_validity.csv", index=False)
    # the guarantee with the true null mean (no baseline error), e-detectors only
    tm = true_null_means()
    rows = []
    for i, (name, kw) in enumerate(NULLS.items()):
        if name not in tm:
            continue
        cal, mon = _split(reps, horizon, seed + 50 + i, **kw)
        for det, key in (("sign", "sign"), ("product c3", "product")):
            alarm, side = run_detector(det, cal, mon, A_SYN, true_m0=tm[name][key])
            fa = alarm >= 0
            rows.append(dict(scenario=name, detector=det, true_m0=tm[name][key], fa_1day=float((fa & (alarm < DAY)).mean()),
                             bound_1day=DAY / A_SYN, fa_5days=float(fa.mean()), bound_5days=horizon / A_SYN, paths=reps))
    df2 = pd.DataFrame(rows)
    df2.to_csv(OUT / "q1_true_m0.csv", index=False)
    pd.DataFrame([dict(naive_h=naive_h, drift=NAIVE_DRIFT, target_per_side_per_day=0.025)]).to_csv(OUT / "q1_naive_threshold.csv", index=False)
    return df, df2, naive_h


# ----------------------------------------------------------------------------
# Q2: delay
# ----------------------------------------------------------------------------


def q2_delay(reps=1000, nu=500, post=3500, seed=200, naive_h=None):
    naive_h = calibrate_naive() if naive_h is None else naive_h
    rows = []
    k = 0
    for inn, ikw in INNOVATIONS.items():
        for rho0, rho1 in TRANSITIONS:
            k += 1
            cal, mon = _split(reps, nu + post, seed + k, nu=nu, rho=rho0, rho1=rho1, **ikw)
            for det in ("sign", "product c3", "product c2", "gaussian LR", "naive CUSUM"):
                alarm, side = run_detector(det, cal, mon, A_SYN, naive_h=naive_h)
                early = (alarm >= 0) & (alarm < nu)
                det_ok = alarm >= nu
                delays = alarm[det_ok] - nu + 1
                want = "up" if rho1 > rho0 else "down"
                rows.append(dict(innovations=inn, rho0=rho0, rho1=rho1, detector=det,
                                 false_alarm_before_change=float(early.mean()), detected=float(det_ok.mean()),
                                 right_direction=float((side[det_ok] == want).mean()) if det_ok.any() else np.nan,
                                 median_delay=float(np.median(delays)) if delays.size else np.nan,
                                 mean_delay=float(delays.mean()) if delays.size else np.nan,
                                 q90_delay=float(np.quantile(delays, 0.9)) if delays.size else np.nan, paths=reps))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q2_delay.csv", index=False)
    return df


def fig_delays(df):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    dets = ["sign", "product c3", "product c2", "gaussian LR", "naive CUSUM"]
    labels = [f"{a:g} to {b:g}" for a, b in TRANSITIONS]
    for ax, inn in zip(axes, INNOVATIONS):
        sub = df[df.innovations == inn]
        x = np.arange(len(TRANSITIONS))
        for j, det in enumerate(dets):
            d = sub[sub.detector == det].set_index(["rho0", "rho1"]).reindex(TRANSITIONS)
            ax.bar(x + (j - 2) * 0.16, d.median_delay.to_numpy(), 0.16, label=det)
        ax.set_xticks(x, labels, rotation=30)
        ax.set_title(f"{inn} innovations")
        ax.set_yscale("log")
    axes[0].set_ylabel("median delay (steps)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Median delay of a correlation change at step 500 [SYNTHETIC], among paths with no alarm before the change\n"
                 "(under t4 and jumps the Gaussian LR and naive CUSUM alarm BEFORE the change on most paths, see q2_delay.csv)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_synthetic_delays.png", dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Q3: real data
# ----------------------------------------------------------------------------

PAIRS = (("BTCUSDT", "ETHUSDT"), ("BTCUSDT", "SOLUSDT"), ("BTCUSDT", "DOGEUSDT"))
START, END = "2024-01-01", "2026-01-01"
BAR = 5
A_REAL = (1440 // BAR) / 0.05  # 288 bars per day, 5% budget per day
CAL_REAL = 1440  # bars (5 days)


def _load_minutes():
    from avmonitor.data import load_binance

    out = {}
    for sym in ("BTCUSDT", "ETHUSDT", "SOLUSDT", "DOGEUSDT"):
        df = load_binance(sym, START, END, download=False)
        out[sym] = df["logret"]
        print("  loaded", sym, df.attrs["checks"], flush=True)
    return out


def q3a_epps(minutes):
    rows = []
    for a, b in PAIRS:
        for bars in (1, 5, 15, 30, 60):
            ra = minutes[a] if bars == 1 else C.aggregate_returns(minutes[a], bars)
            rb = minutes[b] if bars == 1 else C.aggregate_returns(minutes[b], bars)
            ok = ra.notna() & rb.notna()
            for period, mask in (("2024-2025", ok), ("2024", ok & (ra.index.year == 2024)), ("2025", ok & (ra.index.year == 2025))):
                x, y = ra[mask].to_numpy(), rb[mask].to_numpy()
                m_sign = float(np.mean(np.sign(x) * np.sign(y)))
                rows.append(dict(pair=f"{a[:-4]}/{b[:-4]}", bar_minutes=bars, period=period, n=int(mask.sum()),
                                 pearson=float(np.corrcoef(x, y)[0, 1]), sign_mean=m_sign,
                                 sign_implied_rho=float(C.rho_from_sign_mean(m_sign)),
                                 zero_share_a=float((x == 0).mean()), zero_share_b=float((y == 0).mean()),
                                 median_a=float(np.median(x)), median_b=float(np.median(y))))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q3a_epps.csv", index=False)
    return df


def _deseason(r5: pd.Series) -> pd.Series:
    from avmonitor.data import intraday_scale

    prof = intraday_scale(r5, lookback_days=20, smooth=3)
    return r5 / prof


def _window_stats(idx, r1, r2, btc, n_post=CAL_REAL, n_pre=CAL_REAL, n_day=288):
    """Dependence and volatility before and after index ``idx`` (NaN if the windows are incomplete)."""
    n = len(r1)
    if idx - n_pre < 0 or idx + 1 + n_post > n:
        return None
    pre = slice(idx - n_pre, idx)
    post = slice(idx + 1, idx + 1 + n_post)
    s_pre = np.nanmean(C.sign_concordance(r1[pre], r2[pre]))
    s_post = np.nanmean(C.sign_concordance(r1[post], r2[post]))
    ok_pre = np.isfinite(r1[pre]) & np.isfinite(r2[pre])
    ok_post = np.isfinite(r1[post]) & np.isfinite(r2[post])
    p_pre = np.corrcoef(r1[pre][ok_pre], r2[pre][ok_pre])[0, 1]
    p_post = np.corrcoef(r1[post][ok_post], r2[post][ok_post])[0, 1]
    v_pre, v_post = C.robust_scale(btc[pre]), C.robust_scale(btc[post])
    before = np.nansum(btc[max(0, idx - n_day) : idx])
    after = np.nansum(btc[idx + 1 : idx + 1 + n_day])
    return dict(sign_pre=s_pre, sign_post=s_post, pearson_pre=p_pre, pearson_post=p_post,
                btc_vol_ratio=v_post / v_pre if v_pre > 0 else np.nan, btc_ret_24h_before=before, btc_ret_24h_after=after)


def _cal_window_stats(since, r1, r2, n_cal=CAL_REAL):
    """Sign mean and Pearson correlation of the ``n_cal`` finite pairs before index ``since`` (the calibration window)."""
    ok = np.isfinite(r1[:since]) & np.isfinite(r2[:since])
    pos = np.nonzero(ok)[0][-n_cal:]
    if pos.size < n_cal // 2:
        return np.nan, np.nan
    return float(np.mean(np.sign(r1[pos]) * np.sign(r2[pos]))), float(np.corrcoef(r1[pos], r2[pos])[0, 1])


def q3b_monitor(minutes, seed=300):
    rng = np.random.default_rng(seed)
    five = {s: C.aggregate_returns(minutes[s], BAR) for s in minutes}
    des = {s: _deseason(five[s]) for s in five}
    btc_raw = five["BTCUSDT"]  # raw 5-minute log returns of BTC for the 24-hour return around alarms
    ev_rows, sum_rows, ctrl_rows = [], [], []
    configs = [("sign", b, None) for b in (0.0, 0.1, 0.2)] + [("product", b, None) for b in (0.0, 0.1, 0.2)] + [("sign", b, "gaussian") for b in (0.0, 0.1)]
    for a, b in PAIRS:
        pair = f"{a[:-4]}/{b[:-4]}"
        r1, r2 = des[a].to_numpy(), des[b].to_numpy()
        btc_des = des["BTCUSDT"].to_numpy()
        btc = btc_raw.to_numpy()
        idx = des[a].index
        days = (idx[-1] - idx[0]).total_seconds() / 86400
        # control: random times with complete windows
        ctrl = []
        for i in rng.integers(CAL_REAL, len(r1) - CAL_REAL - 1, 2000):
            w = _window_stats(int(i), r1, r2, btc_des)
            if w is not None:
                w["btc_ret_24h_before"] = float(np.nansum(btc[max(0, i - 288) : i]))
                w["btc_ret_24h_after"] = float(np.nansum(btc[i + 1 : i + 289]))
                ctrl.append(w)
        cdf = pd.DataFrame(ctrl)
        ctrl_rows.append(dict(pair=pair, n=len(cdf), share_sign_up=float((cdf.sign_post > cdf.sign_pre).mean()),
                              median_sign_change=float((cdf.sign_post - cdf.sign_pre).median()),
                              median_pearson_change=float((cdf.pearson_post - cdf.pearson_pre).median()),
                              median_btc_vol_ratio=float(cdf.btc_vol_ratio.median()),
                              share_btc_vol_up=float((cdf.btc_vol_ratio > 1).mean()),
                              mean_btc_ret_24h_before=float(cdf.btc_ret_24h_before.mean()),
                              share_btc_down_24h_before=float((cdf.btc_ret_24h_before < 0).mean()),
                              p10_btc_ret_24h_before=float(cdf.btc_ret_24h_before.quantile(0.1)),
                              mean_btc_ret_24h_after=float(cdf.btc_ret_24h_after.mean())))
        p10 = float(cdf.btc_ret_24h_before.quantile(0.1))
        for stat, band, bench in configs:
            t0 = time.time()
            mon = C.CorrelationMonitor(threshold=A_REAL, stat=stat, cal_len=CAL_REAL, null_band=band, benchmark=bench)
            events = mon.run(r1, r2)
            name = "gaussian LR" if bench else stat
            recs = []
            for e in events:
                w = _window_stats(e.index, r1, r2, btc_des)
                s_cal, p_cal = _cal_window_stats(e.monitored_since, r1, r2)
                rec = dict(pair=pair, monitor=name, band=band, time=idx[e.index], index=e.index, direction=e.direction,
                           rho_hat=e.rho_hat, m0=e.m0, monitored_steps=e.index - e.monitored_since + 1,
                           sign_cal=s_cal, pearson_cal=p_cal)
                if w is not None:
                    w["btc_ret_24h_before"] = float(np.nansum(btc[max(0, e.index - 288) : e.index]))
                    w["btc_ret_24h_after"] = float(np.nansum(btc[e.index + 1 : e.index + 289]))
                    rec.update(w)
                recs.append(rec)
                ev_rows.append(rec)
            edf = pd.DataFrame(recs)
            row = dict(pair=pair, monitor=name, band=band, alarms=len(events), per_day=len(events) / days,
                       up=int(sum(e.direction == "up" for e in events)), down=int(sum(e.direction == "down" for e in events)),
                       segments=len(mon.segments), seconds=round(time.time() - t0, 1))
            if len(edf) and "sign_post" in edf:
                for d in ("up", "down"):
                    s = edf[(edf.direction == d) & edf.sign_post.notna()]
                    if len(s) == 0:
                        continue
                    # yardstick for every monitor: the sign mean (and Pearson correlation) of the 5 days after
                    # the alarm against the same quantity in the segment's own calibration window
                    better = (s.sign_post > s.sign_cal) if d == "up" else (s.sign_post < s.sign_cal)
                    better_p = (s.pearson_post > s.pearson_cal) if d == "up" else (s.pearson_post < s.pearson_cal)
                    row[f"{d}_n"] = len(s)
                    row[f"{d}_share_sign_moved_as_flagged"] = float(better.mean())
                    row[f"{d}_share_pearson_moved_as_flagged"] = float(better_p.mean())
                    row[f"{d}_median_sign_change_post_minus_cal"] = float((s.sign_post - s.sign_cal).median())
                    row[f"{d}_median_pearson_change_post_minus_cal"] = float((s.pearson_post - s.pearson_cal).median())
                    row[f"{d}_median_sign_change_post_minus_pre"] = float((s.sign_post - s.sign_pre).median())
                    row[f"{d}_median_btc_vol_ratio"] = float(s.btc_vol_ratio.median())
                    row[f"{d}_share_btc_vol_up"] = float((s.btc_vol_ratio > 1).mean())
                    row[f"{d}_mean_btc_ret_24h_before"] = float(s.btc_ret_24h_before.mean())
                    row[f"{d}_share_btc_down_24h_before"] = float((s.btc_ret_24h_before < 0).mean())
                    row[f"{d}_share_btc_below_p10_24h_before"] = float((s.btc_ret_24h_before < p10).mean())
                    row[f"{d}_mean_btc_ret_24h_after"] = float(s.btc_ret_24h_after.mean())
            sum_rows.append(row)
            print(f"  {pair} {name} band {band}: {len(events)} alarms ({row['up']} up, {row['down']} down) in {row['seconds']} s", flush=True)
    # 1-minute sign monitor on BTC/ETH for the Epps comparison (A = 28800, 5% per 1440-minute day)
    r1, r2 = minutes["BTCUSDT"].to_numpy(), minutes["ETHUSDT"].to_numpy()
    mon = C.CorrelationMonitor(threshold=1440 / 0.05, stat="sign", cal_len=1440 * 5, null_band=0.1)
    ev = mon.run(r1, r2)
    days = (minutes["BTCUSDT"].index[-1] - minutes["BTCUSDT"].index[0]).total_seconds() / 86400
    sum_rows.append(dict(pair="BTC/ETH (1-minute bars)", monitor="sign", band=0.1, alarms=len(ev), per_day=len(ev) / days,
                         up=int(sum(e.direction == "up" for e in ev)), down=int(sum(e.direction == "down" for e in ev)),
                         segments=len(mon.segments)))
    events = pd.DataFrame(ev_rows)
    events.to_csv(OUT / "q3b_events.csv", index=False)
    summary = pd.DataFrame(sum_rows)
    summary.to_csv(OUT / "q3b_summary.csv", index=False)
    pd.DataFrame(ctrl_rows).to_csv(OUT / "q3b_control.csv", index=False)
    # increase alarms of the sign monitor (band 0.1) that are matched by an increase alarm in another pair
    # within 2 days (each pair recalibrates for 5 days after an alarm, so same-day matches are rare by design)
    s = events[(events.monitor == "sign") & (events.band == 0.1) & (events.direction == "up")].copy()
    s["t"] = pd.to_datetime(s.time)
    btc_daily = btc_raw.groupby(btc_raw.index.floor("D")).sum()
    rows = []
    for _, e in s.iterrows():
        near = s[(s.pair != e.pair) & ((s.t - e.t).abs() <= pd.Timedelta(days=2))]
        d0 = e.t.floor("D")
        rows.append(dict(pair=e.pair, time=e.t, matched_pairs=",".join(sorted(set(near.pair))), n_matched=near.pair.nunique(),
                         btc_log_return_that_day=float(btc_daily.get(d0, np.nan)),
                         btc_log_return_prev_day=float(btc_daily.get(d0 - pd.Timedelta(days=1), np.nan)),
                         btc_ret_24h_before=e.btc_ret_24h_before, btc_ret_24h_after=e.btc_ret_24h_after))
    cd = pd.DataFrame(rows)
    cd.to_csv(OUT / "q3b_common_up_days.csv", index=False)
    day = s.t.dt.floor("D")
    # all up-alarm days of the band-0.1 sign monitor per pair, with the BTC daily return (for the record)
    up_days = s.assign(day=day).groupby("day").agg(pairs=("pair", lambda x: ",".join(sorted(set(x)))), n=("pair", "size"))
    up_days["btc_log_return_that_day"] = [float(btc_daily.get(d, np.nan)) for d in up_days.index]
    up_days.to_csv(OUT / "q3b_up_alarm_days.csv")
    fig_real(des, events)
    return summary, events, cd


def fig_real(des, events):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
    for ax, (a, b) in zip(axes, PAIRS):
        pair = f"{a[:-4]}/{b[:-4]}"
        rho = C.rolling_sign_rho(des[a], des[b], CAL_REAL)
        ax.plot(rho.index, rho.to_numpy(), lw=0.7, color="0.3", label="5-day rolling sign-implied correlation")
        ev = events[(events.monitor == "sign") & (events.band == 0.1) & (events.pair == pair)]
        for d, col, mk in (("up", "tab:red", "^"), ("down", "tab:blue", "v")):
            e = ev[ev.direction == d]
            t = pd.to_datetime(e.time)
            ax.plot(t, rho.reindex(t, method="nearest").to_numpy(), mk, color=col, ms=6, label=f"{d} alarm")
        ax.set_title(f"{pair}, 5-minute bars, sign monitor, band 0.1, A = {A_REAL:.0f} [REAL]")
        ax.set_ylim(-0.1, 1.0)
        ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_real_rho_alarms.png", dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Q4: a basket of k assets
# ----------------------------------------------------------------------------

K = 4
PAIRS_K = [(i, j) for i in range(K) for j in range(i + 1, K)]


def _basket_split(reps, n, seed, corr0, nu=None, corr1=None, vol_post=None, df=None):
    from simulate_pairs import basket_returns

    rng = np.random.default_rng(seed)
    b = basket_returns(reps, CAL + n, corr0, None if nu is None else CAL + nu, corr1, vol_post, df, rng)
    return b.r[:, :CAL, :], b.r[:, CAL:, :]


def run_basket(cal, mon, A, band=0.0):
    """Two-sided basket e-detector at level A (each side at 2A); returns first alarm and side per path."""
    lo, hi = C.basket_bounds(K)
    y = C.basket_concordance(mon)
    m0 = np.mean(C.basket_concordance(cal), axis=1)
    up = C.bounded_edetector(y, m0 + band, lo, hi, "up", threshold=2 * A)
    dn = C.bounded_edetector(y, m0 - band, lo, hi, "down", threshold=2 * A)
    both = C.two_sided_alarm(up, dn)
    side = np.where(both < 0, "", np.where((up.alarm >= 0) & ((dn.alarm < 0) | (up.alarm <= dn.alarm)), "up", "down"))
    return both, side


def run_any_pair(cal, mon, A, band=0.0):
    """Six pairwise sign e-detectors, each two-sided at level 6A (sides at 12A), first alarm of any."""
    m = len(PAIRS_K)
    best = np.full(mon.shape[0], -1)
    for i, j in PAIRS_K:
        y = C.sign_concordance(mon[:, :, i], mon[:, :, j])
        m0 = np.mean(C.sign_concordance(cal[:, :, i], cal[:, :, j]), axis=1)
        up = C.bounded_edetector(y, m0 + band, -1.0, 1.0, "up", threshold=2 * m * A)
        dn = C.bounded_edetector(y, m0 - band, -1.0, 1.0, "down", threshold=2 * m * A)
        a = C.two_sided_alarm(up, dn)
        take = (a >= 0) & ((best < 0) | (a < best))
        best[take] = a[take]
    return best


def q4_basket_synthetic(reps=2000, horizon=5 * DAY, seed=400):
    from simulate_pairs import equicorrelation

    E = equicorrelation
    rows = []
    nulls = {
        "gaussian equicorrelation 0.3": dict(corr0=E(K, 0.3)),
        "student t4 equicorrelation 0.3": dict(corr0=E(K, 0.3), df=4.0),
        "vol x2 asset 1 at step 500 (correlations fixed)": dict(corr0=E(K, 0.3), nu=500, vol_post=np.array([2.0, 1, 1, 1])),
        "vol x2 all assets at step 500 (correlations fixed)": dict(corr0=E(K, 0.3), nu=500, vol_post=np.full(K, 2.0)),
    }
    for i, (name, kw) in enumerate(nulls.items()):
        cal, mon = _basket_split(reps, horizon, seed + i, **kw)
        for band in (0.0, 0.05):
            a, side = run_basket(cal, mon, A_SYN, band)
            fa = a >= 0
            rows.append(dict(scenario=name, detector="basket", band=band, fa_1day=float((fa & (a < DAY)).mean()),
                             fa_5days=float(fa.mean()), share_up=float((side == "up").sum() / max(fa.sum(), 1)), paths=reps))
            a = run_any_pair(cal, mon, A_SYN, band)
            fa = a >= 0
            rows.append(dict(scenario=name, detector="any of 6 pairs at 6A", band=band, fa_1day=float((fa & (a < DAY)).mean()),
                             fa_5days=float(fa.mean()), share_up=np.nan, paths=reps))
    val = pd.DataFrame(rows)
    val.to_csv(OUT / "q4_basket_validity.csv", index=False)
    # detection: basket-wide change, one-pair change, basket-wide decrease
    one = E(K, 0.3)
    one[0, 1] = one[1, 0] = 0.8
    changes = {
        "all pairs 0.3 to 0.6": dict(corr0=E(K, 0.3), corr1=E(K, 0.6)),
        "all pairs 0.3 to 0.0": dict(corr0=E(K, 0.3), corr1=E(K, 0.0)),
        "one pair 0.3 to 0.8, others fixed": dict(corr0=E(K, 0.3), corr1=one),
        "all pairs 0.3 to 0.9": dict(corr0=E(K, 0.3), corr1=E(K, 0.9)),
    }
    rows = []
    nu, post = 500, 3500
    for i, (name, kw) in enumerate(changes.items()):
        cal, mon = _basket_split(1000, nu + post, seed + 50 + i, nu=nu, **kw)
        for det in ("basket", "any of 6 pairs at 6A"):
            a = run_basket(cal, mon, A_SYN)[0] if det == "basket" else run_any_pair(cal, mon, A_SYN)
            early = (a >= 0) & (a < nu)
            ok = a >= nu
            d = a[ok] - nu + 1
            rows.append(dict(change=name, detector=det, false_alarm_before_change=float(early.mean()), detected=float(ok.mean()),
                             median_delay=float(np.median(d)) if d.size else np.nan, q90_delay=float(np.quantile(d, 0.9)) if d.size else np.nan,
                             paths=1000))
    delay = pd.DataFrame(rows)
    delay.to_csv(OUT / "q4_basket_delay.csv", index=False)
    return val, delay


def q4_basket_real(minutes, seed=500):
    syms = ("BTCUSDT", "ETHUSDT", "SOLUSDT", "DOGEUSDT")
    five = {s: C.aggregate_returns(minutes[s], BAR) for s in syms}
    R = np.column_stack([five[s].to_numpy() for s in syms])  # raw 5-minute returns: the statistic is scale-free
    idx = five["BTCUSDT"].index
    btc = five["BTCUSDT"].to_numpy()
    y = C.basket_concordance(R)
    days = (idx[-1] - idx[0]).total_seconds() / 86400
    rng = np.random.default_rng(seed)
    ctrl = rng.integers(CAL_REAL, len(y) - CAL_REAL - 1, 2000)
    c_before = np.array([np.nansum(btc[max(0, i - 288) : i]) for i in ctrl])
    p10 = float(np.quantile(c_before, 0.1))
    c_change = np.array([np.nanmean(y[i + 1 : i + 1 + CAL_REAL]) - np.nanmean(y[i - CAL_REAL : i]) for i in ctrl])
    rows, ev_rows = [], []
    for band in (0.05, 0.1, 0.2):
        mon = C.BasketMonitor(K, threshold=A_REAL, cal_len=CAL_REAL, null_band=band)
        events = mon.run(R)
        recs = []
        for e in events:
            ok = np.isfinite(y[: e.monitored_since])
            cal_pos = np.nonzero(ok)[0][-CAL_REAL:]
            y_cal = float(np.mean(y[cal_pos]))
            y_post = float(np.nanmean(y[e.index + 1 : e.index + 1 + CAL_REAL])) if e.index + 1 + CAL_REAL <= len(y) else np.nan
            recs.append(dict(band=band, time=idx[e.index], index=e.index, direction=e.direction, m0=e.m0,
                             basket_cal=y_cal, basket_post=y_post, btc_ret_24h_before=float(np.nansum(btc[max(0, e.index - 288) : e.index])),
                             btc_ret_24h_after=float(np.nansum(btc[e.index + 1 : e.index + 289]))))
        edf = pd.DataFrame(recs)
        ev_rows.append(edf)
        row = dict(band=band, alarms=len(events), per_day=len(events) / days, up=int((edf.direction == "up").sum()) if len(edf) else 0,
                   down=int((edf.direction == "down").sum()) if len(edf) else 0)
        for d in ("up", "down"):
            s = edf[(edf.direction == d) & edf.basket_post.notna()] if len(edf) else edf
            if len(s) == 0:
                continue
            moved = (s.basket_post > s.basket_cal) if d == "up" else (s.basket_post < s.basket_cal)
            row[f"{d}_n"] = len(s)
            row[f"{d}_share_moved_as_flagged"] = float(moved.mean())
            row[f"{d}_median_change_post_minus_cal"] = float((s.basket_post - s.basket_cal).median())
            row[f"{d}_mean_btc_ret_24h_before"] = float(s.btc_ret_24h_before.mean())
            row[f"{d}_share_btc_down_24h_before"] = float((s.btc_ret_24h_before < 0).mean())
            row[f"{d}_share_btc_below_p10_24h_before"] = float((s.btc_ret_24h_before < p10).mean())
            row[f"{d}_mean_btc_ret_24h_after"] = float(s.btc_ret_24h_after.mean())
        rows.append(row)
        print(f"  basket band {band}: {len(events)} alarms ({row['up']} up, {row['down']} down)", flush=True)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "q4_basket_real_summary.csv", index=False)
    pd.concat(ev_rows).to_csv(OUT / "q4_basket_real_events.csv", index=False)
    pd.DataFrame([dict(n=len(ctrl), share_up=float((c_change > 0).mean()), median_change=float(np.median(c_change)),
                       mean_btc_ret_24h_before=float(c_before.mean()), share_btc_down_24h_before=float((c_before < 0).mean()),
                       p10_btc_ret_24h_before=p10, basket_mean_full_sample=float(np.nanmean(y)))]).to_csv(OUT / "q4_basket_real_control.csv", index=False)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="q0,q1,q2,q3,q4")
    args = ap.parse_args()
    parts = args.only.split(",")
    t0 = time.time()
    naive_h = None
    if "q0" in parts:
        df = q0_theory()
        print(df[["rho0", "rho1", "statistic", "shift", "snr", "lam_max", "best_grid_bet", "predicted_delay_grid"]].to_string(), flush=True)
    if "q1" in parts:
        df, df2, naive_h = q1_validity()
        print("naive CUSUM threshold", naive_h)
        print(df.to_string(), flush=True)
        print(df2.to_string(), flush=True)
    if "q2" in parts:
        df = q2_delay(naive_h=naive_h)
        print(df.to_string(), flush=True)
        fig_delays(df)
    print(f"synthetic parts done in {time.time() - t0:.0f} s", flush=True)
    if "q3" in parts:
        minutes = _load_minutes()
        epps = q3a_epps(minutes)
        print(epps.to_string(), flush=True)
        summary, events, cd = q3b_monitor(minutes)
        print(summary.to_string(), flush=True)
        print(cd.to_string(), flush=True)
    if "q4" in parts:
        val, delay = q4_basket_synthetic()
        print(val.to_string(), flush=True)
        print(delay.to_string(), flush=True)
        minutes = _load_minutes() if "q3" not in parts else minutes
        print(q4_basket_real(minutes).to_string(), flush=True)
    print(f"all done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
