"""Benchmark against Podolskij and Vetter (2009) and the block-size rule for the monitor. All synthetic.

Research questions
  Q5 [SYNTHETIC] Noise and jumps together. On a grid of noise ratios x jump
     intensities x jump sizes with stochastic volatility, compare for integrated
     quarticity (and variance) the raw MAD-clipped estimator, pre-averaging then
     MAD drop or clip (ours), the modulated bipower / tripower variation of
     Podolskij and Vetter (2009) with noise correction, and the pre-averaged
     bipower / tripower variation: bias (mean ratio), relative RMSE, and the
     coverage of the 95% interval for integrated variance from the modulated
     bipower CLT (variance from the independent Gaussian block approximation).
  Q6 [SYNTHETIC] The block monitor's delay as a function of the block size k
     and the noise ratio gamma, its validity at small k, and a rule of thumb
     for k from an estimated gamma.

Design as in run_experiment.py: 1-second returns, 23,400 per day, daily
volatility 1%, log-OU spot variance (vol of vol 0.5) for Q5, constant
volatility for Q6, price noise with ratio gamma, compound Poisson jumps.
Seeds fixed, two worker processes, about 2 minutes.
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
from avmonitor.robust import MAD_NORMAL_SCALE  # noqa: E402
from avmonitor.simulate import sv_jump_returns  # noqa: E402

OUT = HERE / "results"
N = 23400
T = 1.0
SIGMA = 0.01
S_STEP = SIGMA / np.sqrt(N)
K = 60
GRID_GAMMA = (0.0, 0.5, 1.0, 2.0)
GRID_JUMPS = ((0, 0.0), (2, 25.0), (2, 100.0), (8, 25.0), (8, 100.0))  # (jumps per day, size in per-second sd)


def _summ(est, truth):
    ratio = est / truth
    return dict(mean_ratio=float(ratio.mean()), rel_rmse=float(np.sqrt(np.mean((ratio - 1) ** 2))))


def q5_job(args):
    gamma, jumps, size, seed = args
    rng = np.random.default_rng(seed)
    p = sv_jump_returns(200, N, T, sigma_bar=SIGMA, kappa=5.0, xi=0.5, jump_intensity=jumps,
                        jump_size_sd=size * S_STEP, noise_sd=gamma * S_STEP, noise_level="price", substeps=1, rng=rng)
    r = p.r
    mbv = P.mbv_estimates(r, c1=0.5, c2=1.0)
    ests = {
        "raw returns, MAD clip (Bapat, normalised)": P.raw_truncated_estimates(r, T, 3.0),
        "pre-averaged, MAD drop (ours, k = 60)": P.preavg_estimates(r, K, T, c=3.0, mode="drop"),
        "pre-averaged, MAD clip (k = 60)": P.preavg_estimates(r, K, T, c=3.0, mode="clip"),
        "modulated bipower / tripower (Podolskij and Vetter 2009, c1 = 0.5, c2 = 1)": mbv,
        "pre-averaged bipower / tripower (k = 60)": P.preavg_multipower_estimates(r, K, T),
        "pre-averaged, no truncation (k = 60)": P.preavg_estimates(r, K, T, c=None),
    }
    rows = []
    for name, e in ests.items():
        row = dict(gamma=gamma, jumps_per_day=jumps, jump_size_sd=size, estimator=name, paths=200,
                   **{f"iq_{a}": b for a, b in _summ(e["iq"], p.iq).items()},
                   **{f"iv_{a}": b for a, b in _summ(e["iv"], p.iv).items()})
        if "iv_se" in e:
            row["iv_ci95_coverage"] = float(np.mean(np.abs(e["iv"] - p.iv) <= 1.96 * e["iv_se"]))
            row["iv_ci95_median_relative_halfwidth"] = float(np.median(1.96 * e["iv_se"] / p.iv))
        rows.append(row)
    return rows


# ----------------------------------------------------------------------------
# Q6 block size
# ----------------------------------------------------------------------------


def _const_paths(reps, days, gamma, rng, ratio_last_day=1.0):
    n = days * N
    sd = np.full(n, S_STEP)
    if ratio_last_day != 1.0:
        sd[(days - 1) * N :] *= ratio_last_day
    x = rng.standard_normal((reps, n)) * sd
    u = rng.standard_normal((reps, n + 1)) * gamma * S_STEP
    return x + np.diff(u, axis=1)


def _block_stream(r_cal, r_mon, k):
    xc, xm = P.preaverage(r_cal, k, k), P.preaverage(r_mon, k, k)
    tau = 3.0 * np.median(np.abs(xc), axis=-1, keepdims=True) * MAD_NORMAL_SCALE
    m0 = np.mean(np.minimum(xc**2, tau**2), axis=-1, keepdims=True)
    return np.minimum(xm**2, tau**2) / m0, float(np.max(tau**2 / m0))


KS = (4, 6, 10, 14, 20, 30, 60)


def q6_job(args):
    gamma, ratio, seed = args
    rng = np.random.default_rng(seed)
    r = _const_paths(200, 3, gamma, rng, ratio_last_day=ratio)
    cal, mon = r[:, : 2 * N], r[:, 2 * N :]
    direction = "up" if ratio > 1 else "down"
    rows = []
    for k in KS:
        y, ub = _block_stream(cal, mon, k)
        m = y.shape[1]
        A = m / 0.05
        a = betting_edetector(y, 1.0, direction, upper=ub, threshold=A).alarm
        det = a >= 0
        delay_s = np.where(det, (a + 1) * k + k - 1, np.nan)
        rows.append(dict(gamma=gamma, sd_ratio=ratio, k=k, blocks_per_day=m, A=A,
                         noise_share_rho=P.block_noise_share(k, gamma),
                         rule_of_thumb_k=P.rule_of_thumb_k(gamma),
                         k_star_formula=(6.0 * gamma / ratio) if ratio != 1.0 else np.nan,
                         detected_within_day=float(det.mean()) if ratio != 1.0 else np.nan,
                         one_day_false_alarm_rate=float(det.mean()) if ratio == 1.0 else np.nan,
                         median_delay_seconds=float(np.nanmedian(delay_s)) if (ratio != 1.0 and det.any()) else np.nan,
                         observed_mean_y=float(y.mean())))
    return rows


# ----------------------------------------------------------------------------


def figure(q6):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cols = {0.5: "#2a78d6", 1.0: "#eb6834", 2.0: "#1baf7a", 4.0: "#eda100"}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for j, ratio in enumerate((2.0, 1.5)):
        d = q6[q6.sd_ratio == ratio]
        for g, col in cols.items():
            s = d[d.gamma == g].sort_values("k")
            ax[j].plot(s.k, s.median_delay_seconds, "-o", ms=4, lw=1.5, color=col, label=f"gamma = {g:g}")
            kr = max(4.0, float(s.k_star_formula.iloc[0]))
            ax[j].axvline(kr, color=col, lw=0.8, ls=":", alpha=0.8)
        ax[j].set_xscale("log")
        ax[j].set_yscale("log")
        ax[j].set_xlabel("block size k (dotted: k* = max(4, 6 gamma / v) from the delay formula)")
        ax[j].set_title(f"sd x {ratio:g} at the start of the day")
        for sp in ("top", "right"):
            ax[j].spines[sp].set_visible(False)
    ax[0].set_ylabel("median detection delay (seconds)")
    ax[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Synthetic: block monitor delay against block size, 5% false alarm budget per day", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_block_size_rule.png", dpi=130)
    plt.close(fig)


def main(parts=("q5", "q6")):
    OUT.mkdir(exist_ok=True)
    with Pool(2) as pool:
        if "q5" in parts:
            jobs = [(g, j, sz, 5000 + 10 * i + jj) for i, g in enumerate(GRID_GAMMA) for jj, (j, sz) in enumerate(GRID_JUMPS)]
            q5 = pd.DataFrame(sum(pool.map(q5_job, jobs), []))
            q5.to_csv(OUT / "q5_benchmark_noise_jumps.csv", index=False)
            print(q5.round(3).to_string(index=False))
        if "q6" in parts:
            jobs = [(g, rt, 6000 + 10 * i + j) for i, g in enumerate((0.5, 1.0, 2.0, 4.0)) for j, rt in enumerate((2.0, 1.5, 1.0))]
            q6 = pd.DataFrame(sum(pool.map(q6_job, jobs), []))
            q6.to_csv(OUT / "q6_block_size.csv", index=False)
            print(q6.round(3).to_string(index=False))
    figure(pd.read_csv(OUT / "q6_block_size.csv"))


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("q5", "q6"))
