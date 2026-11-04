"""R2 ROC-style figures: median delay against measured false-alarm rate, variance x2 at step 400.
Usage: python r2_roc.py <rows.csv> <out_dir>"""
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = pd.read_csv(sys.argv[1], keep_default_na=False, na_values=[""])
out = sys.argv[2]
METHODS = ["README: plug-in, band 0.25", "plug-in, band 0", "confidence-sequence baseline", "known baseline (theorem)",
           "Gaussian LR monitor", "Page CUSUM up", "repeated 60-step z-test"]
LABEL = ["README setting (plug-in, band 0.25)", "Plug-in, band 0", "Confidence-sequence baseline", "Known baseline (theorem)",
         "Gaussian likelihood ratio", "Page CUSUM (up)", "Repeated 60-step z-test"]
COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]  # validated categorical order
MARK = ["o", "s", "^", "D", "v", "P", "X"]
LS = ["-", "-", "-", "-", "--", "--", "--"]  # solid: e-detectors, dashed: competitors without a guarantee
LAWN = {"gaussian": "Gaussian returns", "student_t3": "Student t, 3 degrees of freedom (unit variance)"}
FLOOR = 5e-5
nul = d[(d.kind == "nullcheck") & (d.horizon == 390)]
chg = d[(d.kind == "change") & np.isclose(d.factor, 2.0)]
n_null = int(nul.paths.iloc[0]); n_cell = int(chg.paths.iloc[0])
for law in ["gaussian", "student_t3"]:
    fig, ax = plt.subplots(figsize=(9, 7.4), dpi=130)
    notes = []
    for m, lab, c, mk, ls in zip(METHODS, LABEL, COL, MARK, LS):
        pts = []
        for _, r in nul[(nul.law == law) & (nul.method == m) & nul.label.str.startswith("roc_")].iterrows():
            q = chg[(chg.law == law) & (chg.method == m) & (chg.label == r.label)]
            if not q.empty and np.isfinite(q.iloc[0].median_delay):
                pts.append((max(r.rate, FLOOR), q.iloc[0].median_delay))
        pts.sort()
        pos = [p for p in pts if p[0] > FLOOR]
        zero = [p for p in pts if p[0] <= FLOOR]
        if pos:
            x, y = zip(*pos)
            ax.plot(x, y, ls=ls, color=c, lw=2, marker=mk, ms=8, mec="white", mew=1.0, label=lab)
        else:
            ax.plot([], [], ls=ls, color=c, lw=2, marker=mk, ms=8, mec="white", mew=1.0, label=lab)
        for x0, y0 in zero:  # no false alarm at all: shown at the left edge, not joined to the curve
            ax.plot(x0, y0, marker=mk, ms=8, mfc="white", mec=c, mew=1.5, ls="none")
        nom = nul[(nul.law == law) & (nul.method == m) & (nul.label == "nominal_2A")]
        if not nom.empty:
            q = chg[(chg.law == law) & (chg.method == m) & (chg.label == "nominal_2A")]
            if not q.empty and np.isfinite(q.iloc[0].median_delay) and q.iloc[0].rate >= 0.5:
                k = int(nom.iloc[0].k)
                xs, ys = (nom.iloc[0].rate if k else 3.5e-5), q.iloc[0].median_delay
                ax.plot(xs, ys, marker="*", ms=17 if k else 22, color=c, mec="black", mew=0.8, ls="none", zorder=5 - (k == 0))
                if xs <= 2e-4:
                    notes.append(f"{lab}: {k} false alarm{'s' if k != 1 else ''}, median delay {ys:.0f}")
            elif not q.empty:
                notes.append(f"{lab} at A: {100 * nom.iloc[0].rate:.2f}% false alarms, detects {100 * q.iloc[0].rate:.0f}% (off the chart)")
    ax.plot([], [], marker="*", ms=14, color="white", mec="black", ls="none", label="Star: the method's proved/nominal threshold A = 7,800")
    ax.axvline(0.05, color="#888888", lw=1, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(2.5e-5, 1.2)
    from matplotlib.ticker import FuncFormatter, FixedLocator, NullFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    yt = [t for t in (10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000) if ax.get_ylim()[0] <= t <= ax.get_ylim()[1]]
    ax.yaxis.set_major_locator(FixedLocator(yt)); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    if notes:
        ax.text(0.02, 0.02, "At A = 7,800 on 20,000 no-change paths:\n" + "\n".join(notes), transform=ax.transAxes, fontsize=8.5,
                va="bottom", ha="left", color="#333333", bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))
    ax.set_xlabel(f"False alarms within the first 390 monitored steps\n(share of {n_null:,} fresh no-change paths; stars at the far left: none or one)", fontsize=11)
    ax.set_ylabel(f"Median detection delay, steps (variance doubles at step 400; {n_cell:,} paths)", fontsize=11)
    ax.set_title(f"[SYNTHETIC] Delay against measured false-alarm rate: {LAWN[law]}", fontsize=12)
    ax.grid(True, which="major", color="#dddddd", lw=0.6)
    ax.grid(True, which="minor", color="#f0f0f0", lw=0.4)
    ax.annotate("5% per day", xy=(0.05, 1), xycoords=("data", "axes fraction"), xytext=(4, -14), textcoords="offset points", fontsize=9, color="#555555")
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{out}/roc_{law}_var_x2.png")
    plt.close(fig)
    print("wrote", f"{out}/roc_{law}_var_x2.png")
