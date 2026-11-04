"""R4 figures: blind share by asset and year; post / pre variance ratio of alarms by asset and year."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
T = pd.read_csv(HERE / "asset_year.csv")
ORDER = ["NIFTY", "BTC", "ETH", "XRP", "SOL", "DOGE", "BITSTAMP"]
NAMES = {"NIFTY": "Nifty 50", "BITSTAMP": "Bitstamp BTC"}
COL = dict(zip(ORDER, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]))
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURF, "axes.facecolor": SURF})


def legend_handles(ax):
    h = [plt.Line2D([], [], color=COL[a], marker="o", lw=2, ms=6, label=NAMES.get(a, a)) for a in ORDER]
    return h


fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for ax, mon, title in zip(axes, ["plugin", "cs"], ["Plug-in (README setting)", "Confidence-sequence (cs) baseline"]):
    ax.axhspan(0.6, 0.8, color="#e8e7e2", zorder=0)
    ax.text(2016.9, 0.79, "60% to 80%: 'blind about 70%'", color=INK2, va="top", fontsize=8)
    for a in ORDER:
        g = T[(T.asset == a) & (T.monitor == mon)]
        ax.plot(g.year, g.blind_share, color=COL[a], lw=2, marker="o", ms=6, mec=SURF, mew=1.5)
    ax.set_title(title, color=INK, loc="left", fontsize=10)
    ax.set_xticks(range(2017, 2027))
    ax.tick_params(axis="x", labelrotation=45)
    ax.grid(axis="y", color="#e8e7e2", lw=0.8)
axes[0].set_ylabel("Share of minutes spent calibrating (blind)")
axes[0].set_ylim(0, 0.85)
fig.legend(handles=legend_handles(axes[0]), loc="lower center", ncol=7, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig(HERE / "fig_blind_share.png", dpi=150)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for ax, k, title in zip(axes, ["up", "down"], ["Up alarms: next hour / previous hour", "Down alarms: next hour / previous hour"]):
    ax.axhline(1.0, color=INK2, lw=1, ls="--")
    ax.text(2016.9, 1.02, "random watched minutes: about 1.0", color=INK2, fontsize=8, va="bottom")
    for a in ORDER:
        g = T[(T.asset == a) & (T.monitor == "plugin")]
        ax.plot(g.year, g[f"{k}_ratio"], color=COL[a], lw=2, marker="o", ms=6, mec=SURF, mew=1.5)
    ax.set_title(title, color=INK, loc="left", fontsize=10)
    ax.set_xticks(range(2017, 2027))
    ax.tick_params(axis="x", labelrotation=45)
    ax.grid(axis="y", color="#e8e7e2", lw=0.8)
axes[0].set_ylabel("Median ratio of realised variance\n(60 min after / 60 min before), plug-in")
fig.legend(handles=legend_handles(axes[0]), loc="lower center", ncol=7, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig(HERE / "fig_post_pre_ratio.png", dpi=150)
print("ok")
