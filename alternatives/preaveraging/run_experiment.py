"""Experiments for pre-averaging before MAD truncation (alternative A1). All synthetic.

Research questions
  Q1 [SYNTHETIC] At what noise-to-signal ratio gamma = omega / (sigma sqrt(delta))
     does the plain MAD-truncated estimator of integrated quarticity break, and
     does pre-averaging followed by truncation fix it (bias and RMSE against the
     true IQ, with and without jumps)? Also: does the noise correction matter,
     does clipping or dropping handle a jump better, and how sensitive is the
     result to the window k?
  Q2 [SYNTHETIC] The same question for the local (spot) variance that the
     monitor uses: raw truncated returns against pre-averaged returns, over a
     300-second window, with and without jumps.
  Q3 [SYNTHETIC] What does the monitor's one-day false alarm rate do under noise
     (a) on raw truncated returns, (b) on overlapping pre-averaged returns, and
     (c) on non-overlapping pre-averaged blocks? Only (a) and (c) satisfy the
     conditional-mean null; (b) is included to show the size of the violation.
  Q4 [SYNTHETIC] With noise, how fast do the two valid monitors (a) and (c)
     detect a volatility change of sd x 1.5, x 2 and x 0.5 within one day?

Design: 1-second returns, n = 23,400 per day (T = 1 day), daily volatility 1%
(log-OU spot variance with vol-of-vol 0.5 for Q1 and Q2, constant for Q3 and
Q4), price-level noise N(0, omega^2) with omega = gamma x per-second return sd,
optional jumps (4 per day, size 50 per-second sd). Pre-averaging window k = 60
(theta = k / sqrt(n) = 0.39), weight g(x) = min(x, 1 - x), truncation at
3 x scaled MAD. Seeds fixed. Runs in about 4 minutes on two cores.
"""

from __future__ import annotations

import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import preavg as P  # noqa: E402
from avmonitor.detectors import betting_edetector  # noqa: E402
from avmonitor.robust import MAD_NORMAL_SCALE, realized_quarticity  # noqa: E402
from avmonitor.simulate import sv_jump_returns  # noqa: E402

OUT = HERE / "results"
N = 23400
T = 1.0
SIGMA = 0.01  # daily volatility
S_STEP = SIGMA / np.sqrt(N)  # per-second return sd
K = 60
JUMP_SIZE = 50 * S_STEP
GAMMAS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0)


def _sv_paths(reps, gamma, jumps, seed, xi=0.5):
    rng = np.random.default_rng(seed)
    return sv_jump_returns(reps, N, T, sigma_bar=SIGMA, kappa=5.0, xi=xi, jump_intensity=jumps,
                           jump_size_sd=JUMP_SIZE, noise_sd=gamma * S_STEP, noise_level="price", substeps=1, rng=rng)


# ----------------------------------------------------------------------------
# Q1 integrated quarticity and variance
# ----------------------------------------------------------------------------


def _summ(est, truth):
    ratio = est / truth
    return dict(mean_ratio=float(ratio.mean()), rel_rmse=float(np.sqrt(np.mean((ratio - 1) ** 2))),
                median_ratio=float(np.median(ratio)))


def q1_job(args):
    gamma, jumps, seed = args
    p = _sv_paths(200, gamma, jumps, seed)
    r = p.r
    rows = []
    ests = {
        "RQ, no truncation": dict(iq=realized_quarticity(r, T), iv=np.sum(r**2, axis=-1)),
        "raw returns, MAD clip (Bapat, normalised)": P.raw_truncated_estimates(r, T, 3.0),
        "pre-averaged, no truncation": P.preavg_estimates(r, K, T, c=None),
        "pre-averaged, MAD clip": P.preavg_estimates(r, K, T, c=3.0, mode="clip"),
        "pre-averaged, MAD drop": P.preavg_estimates(r, K, T, c=3.0, mode="drop"),
        "pre-averaged, MAD clip, no noise correction": P.preavg_estimates(r, K, T, c=3.0, mode="clip", noise_correct=False),
    }
    for name, e in ests.items():
        rows.append(dict(gamma=gamma, jumps_per_day=jumps, estimator=name, k=K, paths=200,
                         theory_raw_iq_factor=(1 + 2 * gamma**2) ** 2,
                         **{f"iq_{k}": v for k, v in _summ(e["iq"], p.iq).items()},
                         **{f"iv_{k}": v for k, v in _summ(e["iv"], p.iv).items()}))
    return rows


def q1_k_job(args):
    gamma, jumps, seed = args
    p = _sv_paths(200, gamma, jumps, seed)
    rows = []
    for k in (20, 40, 60, 100, 150):
        for mode in ("clip", "drop"):
            e = P.preavg_estimates(p.r, k, T, c=3.0, mode=mode)
            rows.append(dict(gamma=gamma, jumps_per_day=jumps, k=k, theta=k / np.sqrt(N), mode=mode,
                             noise_term_over_sigma2=float(np.mean(e["c_noise"]) / SIGMA**2),
                             **{f"iq_{kk}": v for kk, v in _summ(e["iq"], p.iq).items()}))
    return rows


# ----------------------------------------------------------------------------
# Q2 spot variance
# ----------------------------------------------------------------------------


def q2_job(args):
    gamma, jumps, seed = args
    p = _sv_paths(100, gamma, jumps, seed)
    r, delta, h = p.r, T / N, 300
    grid = np.arange(2 * h, N, 600)
    truth = p.spot_var[:, grid]
    tau = 3.0 * np.median(np.abs(r), axis=-1, keepdims=True) * MAD_NORMAL_SCALE
    x = np.minimum(r**2, tau**2)
    f2, _ = P.truncation_factors(3.0, "clip")
    raw = np.stack([x[:, t - h + 1 : t + 1].mean(axis=-1) for t in grid], axis=1) / (f2 * delta)
    rbar = P.preaverage(r, K, 1)
    V, _ = P.local_moments(rbar, K, delta, P.noise_variance(r, 3.0), 3.0, "clip", True)
    ends = np.arange(V.shape[-1]) + K - 1
    over = np.stack([V[:, (ends > t - h) & (ends <= t)].mean(axis=-1) for t in grid], axis=1)
    blk = (np.arange(V.shape[-1]) % K) == 0
    block = np.stack([V[:, blk & (ends > t - h) & (ends <= t)].mean(axis=-1) for t in grid], axis=1)
    rows = []
    for name, est in (("raw returns, MAD clip", raw), ("pre-averaged, overlapping", over), ("pre-averaged, blocks", block)):
        ratio = est / truth
        rows.append(dict(gamma=gamma, jumps_per_day=jumps, estimator=name, window_seconds=h, k=K,
                         mean_ratio=float(ratio.mean()), median_abs_log_error=float(np.median(np.abs(np.log(np.maximum(ratio, 1e-12))))),
                         p95_abs_log_error=float(np.quantile(np.abs(np.log(np.maximum(ratio, 1e-12))), 0.95))))
    return rows


# ----------------------------------------------------------------------------
# Q3 and Q4: the monitor under noise
# ----------------------------------------------------------------------------


def _const_paths(reps, days, gamma, rng, ratio_last_day=1.0):
    """Constant volatility, price noise, optional sd change on the last day."""
    n = days * N
    sd = np.full(n, S_STEP)
    if ratio_last_day != 1.0:
        sd[(days - 1) * N :] *= ratio_last_day
    x = rng.standard_normal((reps, n)) * sd
    u = rng.standard_normal((reps, n + 1)) * gamma * S_STEP
    return x + np.diff(u, axis=1)


def _streams(r_cal, r_mon, k, step):
    """Normalised monitored stream ``y`` (null mean 1) and conservative bound.

    ``k=None`` is the raw-return monitor; otherwise returns are pre-averaged
    with window ``k`` and taken every ``step`` windows (``step=k`` gives blocks).
    """
    if k is None:
        xc, xm = r_cal, r_mon
    else:
        xc, xm = P.preaverage(r_cal, k, step), P.preaverage(r_mon, k, step)
    tau = 3.0 * np.median(np.abs(xc), axis=-1, keepdims=True) * MAD_NORMAL_SCALE
    m0 = np.mean(np.minimum(xc**2, tau**2), axis=-1, keepdims=True)
    y = np.minimum(xm**2, tau**2) / m0
    return y, float(np.max(tau**2 / m0))


# (label, k, step, conditional-mean null holds under the model)
MONITORS = (("raw returns", None, 1, "no under price noise (MA(1) returns)"),
            ("overlapping pre-averaged, k = 60", 60, 1, "no"),
            ("pre-averaged blocks, k = 60", 60, 60, "yes"),
            ("pre-averaged blocks, k = 20", 20, 20, "yes"),
            ("pre-averaged blocks, k = 10", 10, 10, "yes"))


def q3_job(args):
    gamma, seed = args
    rng = np.random.default_rng(seed)
    r = _const_paths(200, 3, gamma, rng)
    cal, mon = r[:, : 2 * N], r[:, 2 * N :]
    rows = []
    for kind, k, step, valid in MONITORS:
        y, ub = _streams(cal, mon, k, step)
        m = y.shape[1]
        A = m / 0.05
        for direction in ("up", "down"):
            a = betting_edetector(y, 1.0, direction, upper=ub, threshold=A).alarm
            rows.append(dict(gamma=gamma, monitor=kind, direction=direction, steps_per_day=m, A=A, bound=0.05,
                             one_day_false_alarm_rate=float((a >= 0).mean()), conditional_null_holds=valid))
    return rows


def q4_job(args):
    gamma, ratio, seed = args
    rng = np.random.default_rng(seed)
    r = _const_paths(200, 3, gamma, rng, ratio_last_day=ratio)
    cal, mon = r[:, : 2 * N], r[:, 2 * N :]
    direction = "up" if ratio > 1 else "down"
    rows = []
    for kind, k, step, valid in MONITORS:
        if step == 1 and k is not None:
            continue  # the overlapping version has no false alarm control, so its delay is meaningless
        y, ub = _streams(cal, mon, k, step)
        m = y.shape[1]
        A = m / 0.05
        a = betting_edetector(y, 1.0, direction, upper=ub, threshold=A).alarm
        det = a >= 0
        delay_s = np.where(det, (a + 1) * step + (k - 1 if k else 0), np.nan)
        rows.append(dict(gamma=gamma, sd_ratio=ratio, direction=direction, monitor=kind, k=k if k else 1, A=A,
                         detected_within_day=float(det.mean()),
                         median_delay_seconds=float(np.nanmedian(delay_s)) if det.any() else np.nan,
                         diluted_variance_ratio_raw=(ratio**2 + 2 * gamma**2) / (1 + 2 * gamma**2),
                         observed_mean_y=float(y.mean())))
    return rows


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

COL = {"raw returns, MAD clip (Bapat, normalised)": "#2a78d6", "pre-averaged, MAD clip": "#eb6834",
       "pre-averaged, MAD drop": "#1baf7a", "pre-averaged, no truncation": "#eda100",
       "raw returns": "#2a78d6", "overlapping pre-averaged, k = 60": "#e87ba4",
       "pre-averaged blocks, k = 60": "#1baf7a", "pre-averaged blocks, k = 20": "#eda100",
       "pre-averaged blocks, k = 10": "#6250d6"}


def figures(q1, q3, q4):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for j, jumps in enumerate((0, 4)):
        d = q1[q1.jumps_per_day == jumps]
        for name in ("raw returns, MAD clip (Bapat, normalised)", "pre-averaged, no truncation",
                     "pre-averaged, MAD clip", "pre-averaged, MAD drop"):
            s = d[d.estimator == name].sort_values("gamma")
            ax[j].plot(s.gamma + 0.05, s.iq_rel_rmse, "-o", ms=4, lw=1.5, color=COL[name], label=name)
        ax[j].set_xscale("log")
        ax[j].set_yscale("log")
        ax[j].set_xlabel("noise-to-signal ratio gamma (0 plotted at 0.05)")
        ax[j].set_title(f"{jumps} jumps per day")
        for sp in ("top", "right"):
            ax[j].spines[sp].set_visible(False)
    ax[0].set_ylabel("relative RMSE of integrated quarticity")
    ax[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Synthetic: 1-second returns, 23,400 per day, 200 paths, k = 60", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_iq_rmse_vs_noise.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    d = q3[q3.direction == "up"]
    for name, _, _, _ in MONITORS:
        s = d[d.monitor == name].sort_values("gamma")
        ax[0].plot(s.gamma + 0.05, s.one_day_false_alarm_rate, "-o", ms=4, lw=1.5, color=COL[name], label=name)
    ax[0].axhline(0.05, color="#898781", lw=1, ls="--")
    ax[0].text(0.06, 0.055, "bound 5%", color="#898781", fontsize=8)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("noise-to-signal ratio gamma (0 plotted at 0.05)")
    ax[0].set_ylabel("one-day false alarm rate (increase side)")
    ax[0].set_title("Q3: false alarms under a constant volatility")
    ax[0].legend(frameon=False, fontsize=8)
    d = q4[q4.sd_ratio == 2.0]
    for name, k, step, _ in MONITORS:
        if step == 1 and k is not None:
            continue
        s = d[d.monitor == name].sort_values("gamma")
        ax[1].plot(s.gamma + 0.05, s.median_delay_seconds, "-o", ms=4, lw=1.5, color=COL[name], label=name)
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("noise-to-signal ratio gamma (0 plotted at 0.05)")
    ax[1].set_ylabel("median detection delay (seconds)")
    ax[1].set_title("Q4: sd doubles at the start of the day")
    ax[1].legend(frameon=False, fontsize=8)
    for a in ax:
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.suptitle("Synthetic: e-detector on 1-second data, 200 paths, 5% false alarm budget per day", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_monitor_under_noise.png", dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------


def main(parts=("q1", "q2", "q3", "q4")):
    OUT.mkdir(exist_ok=True)
    with Pool(2) as pool:
        if "q1" in parts:
            jobs = [(g, j, 1000 + 10 * i + jj) for i, g in enumerate(GAMMAS) for jj, j in enumerate((0, 4))]
            q1 = pd.DataFrame(sum(pool.map(q1_job, jobs), []))
            q1.to_csv(OUT / "q1_iq_iv_vs_noise.csv", index=False)
            print(q1.round(3).to_string(index=False))
            qk = pd.DataFrame(sum(pool.map(q1_k_job, [(1.0, 0, 1500), (1.0, 4, 1501), (4.0, 4, 1502)]), []))
            qk.to_csv(OUT / "q1_window_sensitivity.csv", index=False)
            print(qk.round(3).to_string(index=False))
        if "q2" in parts:
            jobs = [(g, j, 2000 + 10 * i + jj) for i, g in enumerate((0.0, 0.5, 1.0, 2.0, 4.0)) for jj, j in enumerate((0, 4))]
            q2 = pd.DataFrame(sum(pool.map(q2_job, jobs), []))
            q2.to_csv(OUT / "q2_spot_variance.csv", index=False)
            print(q2.round(3).to_string(index=False))
        if "q3" in parts:
            q3 = pd.DataFrame(sum(pool.map(q3_job, [(g, 3000 + i) for i, g in enumerate((0.0, 0.5, 1.0, 2.0, 4.0))]), []))
            q3.to_csv(OUT / "q3_monitor_false_alarms.csv", index=False)
            print(q3.round(3).to_string(index=False))
        if "q4" in parts:
            jobs = [(g, rt, 4000 + 10 * i + j) for i, g in enumerate((0.0, 0.5, 1.0, 2.0, 4.0)) for j, rt in enumerate((1.5, 2.0, 0.5))]
            q4 = pd.DataFrame(sum(pool.map(q4_job, jobs), []))
            q4.to_csv(OUT / "q4_monitor_detection.csv", index=False)
            print(q4.round(3).to_string(index=False))
    figures(pd.read_csv(OUT / "q1_iq_iv_vs_noise.csv"), pd.read_csv(OUT / "q3_monitor_false_alarms.csv"),
            pd.read_csv(OUT / "q4_monitor_detection.csv"))


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("q1", "q2", "q3", "q4"))
