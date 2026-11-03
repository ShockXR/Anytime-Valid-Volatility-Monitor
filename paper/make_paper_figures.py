"""Figures for paper/short_paper.md.

Reads only CSV files written by the experiment scripts (no new simulation) and
writes PNG files to paper/figures/. Run from the repository root:

    python paper/make_paper_figures.py

Every figure is drawn here, so that paper/figures/ is self-contained and no figure
carries an internal label in its title.

Figure 1, validity.png (simulated data).
    Left: false alarm rate within m steps divided by the bound m/A (column
    rate_over_bound of results/e09ab_bound_tightness.csv), with 95% intervals
    (ci_low and ci_high divided by the bound), on a linear axis.
    Right: one-day false alarm rate under GARCH(1,1) clustering as a function of
    the band b (column one_day_false_alarm of results/e09c_garch_band.csv), with
    95% intervals, one line per persistence.

Figure S1, delay_predictions.png (simulated data).
    Ratio of predicted to simulated mean delay for the 20 cases of
    results/e11a_delay.csv: first-order formula and proved bound of Theorem 18(a)
    for the default grid (divided by sim_nu0_grid_mean), and the second-order
    formula for the best single bet (divided by sim_nu0_oracle_mean).

Figure S2, signal_cap.png (simulated data).
    The ratio of means rho(theta) = theta (1 - beta) / (1 - alpha - beta + alpha theta)
    after a break of factor theta in the long-run variance; the filled circles are
    the simulated ratio of means, pooled over paths (column ratio_of_means_pooled,
    with two bootstrap standard errors as bars), and the open triangles the
    simulated mean of the ratio (column mean_of_ratio_sim), both from
    results/e10d_signal_cap.csv; the dashed lines are the bounds (1 - beta) / alpha
    on the ratio of means. The column ratio_of_means_sim (an average of per-path
    ratios, biased upward) is not plotted.
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

# Fixed categorical order (checked for colour-vision deficiency); every series
# also has its own marker, so identity never rests on colour alone.
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK = "#3a3a38"
GRID = "#d9d8d3"

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.6,
        "legend.frameon": False,
    }
)


def _grid(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def rho(theta, a, b):
    """Post-break ratio of the mean squared return to the mean pre-break forecast."""
    return theta * (1.0 - b) / (1.0 - a - b + a * theta)


def fig_validity():
    t = pd.read_csv(RES / "e09ab_bound_tightness.csv")
    g = pd.read_csv(RES / "e09c_garch_band.csv")
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(7.6, 3.4))

    laws = [
        ("gaussian", "Gaussian", C1, "o"),
        ("student_t3", "Student t3", C2, "s"),
        ("jumps", "jumps", C3, "^"),
        ("ticks", "price ticks", C4, "D"),
    ]
    offsets = [0.96, 0.987, 1.013, 1.04]  # small horizontal offsets so that the intervals do not overlap
    for (key, label, col, mk), off in zip(laws, offsets):
        d = t[t["distribution"] == key].sort_values("horizon")
        x = d["horizon"].to_numpy() * off
        y = d["rate_over_bound"].to_numpy()
        lo = y - d["ci_low"].to_numpy() / d["bound"].to_numpy()
        hi = d["ci_high"].to_numpy() / d["bound"].to_numpy() - y
        axl.errorbar(x, y, yerr=[lo, hi], color=col, marker=mk, ms=4.5, lw=1.4, capsize=2, label=label)
    axl.axhline(1.0, color=INK, ls="--", lw=0.9)
    axl.text(110, 1.02, "bound m/A", fontsize=8, color=INK, va="bottom")
    axl.set_xscale("log")
    axl.set_xticks([98, 390, 1950, 3900, 7800])
    axl.set_xticklabels(["98", "390", "1950", "3900", "7800"])
    axl.minorticks_off()
    axl.set_ylim(0, 1.1)
    axl.set_xlabel("horizon m (steps; 390 steps = one day)")
    axl.set_ylabel("false alarm rate / (m/A)")
    axl.set_title("Null holds: true baseline, A = 7800", fontsize=9, color=INK)
    axl.legend(fontsize=8, loc="lower right")
    _grid(axl)

    pers = [(0.95, C1, "o"), (0.99, C2, "s"), (0.999, C3, "^"), (0.9999, C4, "D")]
    for p, col, mk in pers:
        d = g[np.isclose(g["persistence"], p)].sort_values("band")
        y = d["one_day_false_alarm"].to_numpy()
        lo = y - d["ci_low"].to_numpy()
        hi = d["ci_high"].to_numpy() - y
        axr.errorbar(d["band"], y, yerr=[lo, hi], color=col, marker=mk, ms=4.5, lw=1.4, capsize=2,
                     label=f"persistence {p}")
    axr.axhline(0.05, color=INK, ls="--", lw=0.9)
    axr.text(0.62, 0.065, "5% (bound if the null held)", fontsize=8, color=INK, va="bottom")
    axr.set_ylim(0, 1.05)
    axr.set_xlabel("band b")
    axr.set_ylabel("share of days with an alarm")
    axr.set_title("Null false: GARCH(1,1), A = 28,800", fontsize=9, color=INK)
    axr.legend(fontsize=8, loc="upper right")
    _grid(axr)

    fig.tight_layout()
    fig.savefig(OUT / "validity.png", dpi=200)
    plt.close(fig)


def fig_delay_predictions():
    d = pd.read_csv(RES / "e11a_delay.csv")
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    series = [
        ("pred_first_order_grid", "sim_nu0_grid_mean", "first order, log A / g (default grid)", C1, "o"),
        ("pred_wald_bound_grid", "sim_nu0_grid_mean", "proved bound, Theorem 18(a) (default grid)", C2, "s"),
        ("pred_second_order_oracle", "sim_nu0_oracle_mean", "second order (best single bet)", C3, "^"),
    ]
    for pred, sim, label, col, mk in series:
        ax.scatter(d[sim], d[pred] / d[sim], color=col, marker=mk, s=28, label=label, zorder=3,
                   edgecolors="white", linewidths=0.6)
    ax.axhline(1.0, color=INK, ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200, 400])
    ax.set_xticklabels(["20", "50", "100", "200", "400"])
    ax.minorticks_off()
    ax.set_xlabel("simulated mean delay (steps)")
    ax.set_ylabel("predicted / simulated mean delay")
    ax.set_ylim(0.8, 2.9)
    ax.legend(fontsize=8, loc="upper left")
    _grid(ax)
    fig.tight_layout()
    fig.savefig(OUT / "delay_predictions.png", dpi=200)
    plt.close(fig)


def fig_signal_cap():
    d = pd.read_csv(RES / "e10d_signal_cap.csv")
    colours = {0.95: (C1, "o"), 0.99: (C2, "s"), 0.999: (C3, "D")}
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    theta = np.logspace(0, 3, 400)
    d["persistence"] = (d["alpha"] + d["beta"]).round(3)
    for (_, a, b), g in d.groupby(["persistence", "alpha", "beta"]):
        p = round(a + b, 3)
        col, mk = colours.get(p, ("k", "o"))
        ax.plot(theta, rho(theta, a, b), color=col, lw=1.8, label=f"closed form, persistence {p}")
        ax.axhline((1 - b) / a, color=col, ls="--", lw=0.9)
        ax.errorbar(
            g["k"],
            g["ratio_of_means_pooled"],
            yerr=2.0 * g["ratio_of_means_pooled_se"],
            color=col,
            fmt=mk,
            ms=4.5,
            capsize=2,
            zorder=3,
        )
        ax.scatter(g["k"], g["mean_of_ratio_sim"], color=col, marker="^", s=26, facecolors="none", zorder=3)
    ax.scatter([], [], color="grey", marker="o", s=22, label="filled markers: simulated ratio of means (pooled)")
    ax.scatter([], [], color="grey", marker="^", s=22, facecolors="none", label="open triangles: simulated mean of the ratio (monitored)")
    ax.plot([], [], color="grey", ls="--", lw=0.9, label="bound (1 - beta) / alpha on the ratio of means")
    ax.set_xscale("log")
    ax.set_xlabel("factor theta of the break in the long-run variance")
    ax.set_ylabel("post-break squared return relative to the forecast")
    ax.set_ylim(0.95, 2.35)
    ax.legend(fontsize=7.5, loc="lower right", bbox_to_anchor=(1.0, 0.25))
    _grid(ax)
    fig.tight_layout()
    fig.savefig(OUT / "signal_cap.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_validity()
    fig_delay_predictions()
    fig_signal_cap()
    for name in ("validity.png", "delay_predictions.png", "signal_cap.png"):
        print("wrote", OUT / name)
