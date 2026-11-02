"""Figures for paper/short_paper.md that have no counterpart in results/figures/.

Reads only CSV files written by the experiment scripts (no new simulation) and
writes PNG files to paper/figures/. Run from the repository root:

    python paper/make_paper_figures.py

Figure 1, signal_cap.png: the signal cap of forecast-relative monitoring (E10, Q10d).
The curve is the closed form rho(k) = k (1 - beta) / (1 - alpha - beta + alpha k);
the markers are the simulated ratio of means and mean of ratios from
results/e10d_signal_cap.csv; the dashed lines are the caps (1 - beta) / alpha.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "paper" / "figures"


def rho(k, a, b):
    """Ratio of the mean squared return to the mean pre-break forecast after a variance break of factor k."""
    return k * (1.0 - b) / (1.0 - a - b + a * k)


def fig_signal_cap():
    d = pd.read_csv(RES / "e10d_signal_cap.csv")
    colours = {0.95: "#1f77b4", 0.99: "#d95f02", 0.999: "#1b9e77"}
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    k = np.logspace(0, 3, 400)
    d["persistence"] = (d["alpha"] + d["beta"]).round(3)
    for (_, a, b), g in d.groupby(["persistence", "alpha", "beta"]):
        p = round(a + b, 3)
        col = colours.get(p, "k")
        ax.plot(k, rho(k, a, b), color=col, lw=1.8, label=f"closed form, alpha + beta = {p}")
        ax.axhline((1 - b) / a, color=col, ls="--", lw=0.9)
        ax.scatter(g["k"], g["ratio_of_means_sim"], color=col, marker="o", s=22, zorder=3)
        ax.scatter(g["k"], g["mean_of_ratio_sim"], color=col, marker="^", s=22, facecolors="none", zorder=3)
    ax.scatter([], [], color="grey", marker="o", s=22, label="simulated ratio of means")
    ax.scatter([], [], color="grey", marker="^", s=22, facecolors="none", label="simulated mean of ratio (monitored)")
    ax.plot([], [], color="grey", ls="--", lw=0.9, label="cap (1 - beta) / alpha")
    ax.set_xscale("log")
    ax.set_xlabel("factor k of the break in the long-run variance")
    ax.set_ylabel("post-break mean of r^2 / forecast")
    ax.set_ylim(0.95, 2.35)
    ax.legend(fontsize=7.5, loc="lower right", bbox_to_anchor=(1.0, 0.2), frameon=False)
    ax.set_title("Signal left after a frozen GARCH(1,1) forecast (E10, synthetic)", fontsize=10)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "signal_cap.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    fig_signal_cap()
    print("wrote", OUT / "signal_cap.png")
