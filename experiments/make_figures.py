"""Figures for the experiment report, built from the CSV files in results/."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from avmonitor import robust as R  # noqa: E402

RES = Path(__file__).resolve().parents[1] / "results"
FIG = RES / "figures"


def fig_e01():
    df = pd.read_csv(RES / "e01_reproduce_tables.csv")
    t7 = df[df.table == 7]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for j, (col, pub, title) in enumerate([("sim_kernel", "pub_kernel", "Kernel estimator, eq. (7)"),
                                           ("sim_robust_rawMAD", "pub_robust", "Robust kernel, eq. (16), raw MAD")]):
        for sj, mk in zip((0.01, 0.05, 0.10), ("o", "s", "^")):
            s = t7[t7.sigma_J == sj]
            ax[j].scatter(s[pub], s[col], marker=mk, label=f"sigma_J = {sj}", alpha=0.8)
        lim = [1e-7, 1e-4]
        ax[j].plot(lim, lim, "k--", lw=0.8)
        ax[j].set_xscale("log"); ax[j].set_yscale("log")
        ax[j].set_xlabel("published RMSE (Table 7)"); ax[j].set_ylabel("our simulated RMSE, formula as printed")
        ax[j].set_title(title)
        ax[j].legend(fontsize=8)
    fig.suptitle("E1: one fitted volatility (sigma = 0.0215 per step) and the printed formulas reproduce Table 7 (synthetic)")
    fig.tight_layout(); fig.savefig(FIG / "e01_table7_reproduction.png", dpi=130)

    hs = np.arange(5, 61)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(hs, [R.kernel_weight_sum(h) for h in hs], label="sum of raw Gaussian kernel weights W(h)")
    ax.plot(hs, [R.kernel_weight_sum(h) * R.truncated_moment_factor(R.effective_cutoff(3, False), 4) for h in hs],
            label="W(h) x 0.705 (after c = 3 raw-MAD truncation)")
    ax.axhline(1, color="k", lw=0.8, ls="--"); ax.text(6, 1.3, "correct scale (weights sum to 1)")
    ax.set_xlabel("bandwidth h"); ax.set_ylabel("estimate / per-step quarticity")
    ax.set_title("Why RMSE grows with h in Tables 7 and 8"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "e01_scale_factor.png", dpi=130)


def fig_e04():
    cd = pd.read_csv(RES / "e04b_change.csv")
    dets_up = ["e-det up (grid)", "e-det up on raw r^2 (no truncation)", "Gaussian LR e-det up", "Page CUSUM up (calibrated)"]
    dets_dn = ["e-det down (grid)", "Gaussian LR e-det down"]
    fig, ax = plt.subplots(1, 3, figsize=(13, 4), sharey=False)
    for k, sc in enumerate(("gaussian", "jumps", "student t5")):
        s = cd[cd.scenario == sc]
        for d in dets_up + dets_dn:
            q = s[(s.detector == d) & ((s.ratio > 1) if d in dets_up else (s.ratio < 1))].sort_values("ratio")
            ax[k].plot(q.ratio, q.false_alarm_before_change, marker="o", label=d)
        ax[k].axhline(500 / 7800, color="k", lw=0.8, ls="--")
        ax[k].set_title(f"{sc}: false alarm before change"); ax[k].set_xlabel("sd ratio after change")
        ax[k].set_ylim(0, 1)
    ax[0].set_ylabel("share of paths"); ax[2].legend(fontsize=7)
    fig.suptitle("E4: dashed line = guaranteed bound 500 / A for 500 pre-change steps (synthetic)")
    fig.tight_layout(); fig.savefig(FIG / "e04_false_alarms.png", dpi=130)


if __name__ == "__main__":
    fig_e01()
    fig_e04()
