"""Figures for the robust IQ estimator study (reads the CSV files in results/)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"

# Reference categorical palette (light mode), fixed slot order.
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

SHOW = [("RK", "Robust kernel (c = 3, daily MAD)"), ("RK-k", "Same, kappa-corrected"),
        ("TRQ-M", "Mancini threshold RQ"), ("PI-LU", "New: local plug-in, sqrt(2 log n), drop")]


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def fig_convergence():
    df = pd.read_csv(OUT / "q2_q4_sim_summary.csv")
    scen = ["S3 medium jumps", "S5 strong SV + intraday pattern"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), sharey=True, facecolor=SURFACE)
    for ax, sc in zip(axes, scen):
        _style(ax)
        g = df[df.scenario == sc]
        ends = []
        for (key, label), col in zip(SHOW, SLOT):
            x = g[g.estimator == key].sort_values("n")
            ax.plot(x["n"], x["rel_rmse"], color=col, linewidth=2, marker="o", markersize=5, label=label)
            ends.append([np.log10(x["rel_rmse"].iloc[-1]), key])
        # direct labels at the right end, nudged apart so they never overlap
        ends.sort()
        for i in range(1, len(ends)):
            ends[i][0] = max(ends[i][0], ends[i - 1][0] + 0.07)
        for y, key in ends:
            ax.annotate(key, (23400, 10**y), xytext=(8, 0), textcoords="offset points", fontsize=8, color=INK, va="center")
        ns = np.array([390, 23400])
        ref = g[g.estimator == "PI-LU"].sort_values("n")["rel_rmse"].iloc[0] * np.sqrt(390 / ns)
        ax.plot(ns, ref, color=INK2, linestyle=":", linewidth=1.2, label="slope -1/2 (rate sqrt(n))")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(300, 60000)
        ax.set_title(sc, fontsize=10, color=INK)
        ax.set_xlabel("returns per day n", fontsize=9, color=INK2)
    axes[0].set_ylabel("relative RMSE of IQ estimate", fontsize=9, color=INK2)
    axes[0].legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.suptitle("[SYNTHETIC] Which estimators converge as sampling gets finer (500 days per point)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_convergence.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_real():
    path = OUT / "q5_real_daily.csv"
    if not path.exists():
        return
    d = pd.read_csv(path)
    d = d[d.bar_minutes == 1]
    keys = ["RQ", "TPQ", "MinRQ", "TRQ-M", "RK", "RK-k", "TRQ-U", "PI-G3", "PI-L3", "PI-LU", "KPI-L3"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), sharey=True, facecolor=SURFACE)
    for ax, sym in zip(axes, ("BTCUSDT", "DOGEUSDT")):
        _style(ax)
        g = d[d.symbol == sym]
        data = [np.log2(g[k] / g["MedRQ"]).to_numpy() for k in keys]
        bp = ax.boxplot(data, orientation="horizontal", widths=0.55, showfliers=False, patch_artist=True,
                        medianprops=dict(color=INK, linewidth=1.5), whiskerprops=dict(color=INK2),
                        capprops=dict(color=INK2), boxprops=dict(edgecolor=INK2))
        for k, box in zip(keys, bp["boxes"]):
            box.set_facecolor(SLOT[0] if k in ("TRQ-U", "PI-L3", "PI-LU", "KPI-L3") else "#c9c8c2")
        ax.axvline(0, color=INK2, linewidth=1)
        ax.set_yticks(range(1, len(keys) + 1), keys, fontsize=8)
        ax.set_title(f"{sym}, 184 days, 1-minute returns", fontsize=10, color=INK)
        ax.set_xlabel("log2(estimate / MedRQ), same day", fontsize=9, color=INK2)
    fig.suptitle("[REAL] Daily IQ relative to MedRQ (blue: estimators with a local threshold, this study)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_real_ratio_to_medrq.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main():
    fig_convergence()
    fig_real()


if __name__ == "__main__":
    main()
