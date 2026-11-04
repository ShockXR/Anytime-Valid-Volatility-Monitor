"""Draw the README's opening figure from the data embedded in docs/demo/replay.html.

The replay page is built by scripts/build_replay.py from E18's primary design, and its
embedded JSON holds, for every one-minute bar of the Nifty 50 case (N1, 4 June 2024):
the clipped and scaled squared return Y_t, the wealth of the increase and decrease sides,
the alarm threshold 2A and the recalibration intervals. This script only reads that JSON;
it runs no experiment. Run from the repository root:

    python docs/figures/make_readme_figure.py

It writes docs/figures/nifty_2024-06-04_alarm.png (the evidence panel alone, sized so its type
stays readable from a 375 px phone screen to a desktop README column) and
docs/figures/nifty_2024-06-04_alarm_detail.png (two panels, adding each minute's clipped return).
"""

import datetime as dt
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "docs" / "demo" / "replay.html"
OUT = ROOT / "docs" / "figures" / "nifty_2024-06-04_alarm_detail.png"
OUT_PHONE = ROOT / "docs" / "figures" / "nifty_2024-06-04_alarm.png"

INK = "#1d2433"
INK2 = "#5b6475"
GRID = "#e6e8ee"
EVIDENCE = "#1f4e8c"
THRESH = "#b23b2e"
BAR = "#8fa3c2"
BAR_RECAL = "#c9ccd3"
RECAL = "#f1f2f5"


def load_case(key: str = "N1") -> dict:
    html = SRC.read_text(encoding="utf-8")
    blob = re.search(r'<script[^>]*id="replay-data"[^>]*>(.*?)</script>', html, re.S).group(1)
    day = next(d for d in json.loads(blob)["days"] if d["key"] == key)
    t = [day["t0"]]
    for step, count in day["dt"]:
        for _ in range(count):
            t.append(t[-1] + step)
    day["clock"] = [dt.datetime(1970, 1, 1) + dt.timedelta(minutes=m + day["tzOffset"]) for m in t]
    return day


def main() -> None:
    d = load_case()
    clock = d["clock"]

    def idx(stamp: str) -> int:
        target = dt.datetime.strptime(stamp, "%Y-%m-%d %H:%M")
        return clock.index(target)

    lo, hi = idx("2024-06-03 14:50"), idx("2024-06-04 10:05")
    open_i, alarm_i = idx("2024-06-04 09:15"), idx("2024-06-04 09:26")
    alarms = [a for a, side in d["methods"][0]["alarms"] if lo <= a <= hi]
    assert alarms == [alarm_i], alarms  # the recorded e-detector up alarm in this window
    assert d["wu"][alarm_i] >= d["thr"] > d["wu"][alarm_i - 1]

    # x positions: one unit per bar, with a gap of 6 units for the overnight break
    gap = 6
    xpos = {i: (i - lo) + (gap if i >= open_i else 0) for i in range(lo, hi + 1)}
    recal = [(max(a, lo), min(b, hi + 1)) for a, b in d["recal"] if b > lo and a <= hi]

    plt.rcParams.update({"font.family": "Inter", "font.size": 15, "axes.edgecolor": GRID,
                         "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.4), dpi=150, sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1.6], "hspace": 0.12})
    fig.patch.set_facecolor("white")

    def shade(ax, label=False):
        for a, b in recal:
            x0, x1 = xpos[a] - 0.5, xpos[b - 1] + 0.5
            ax.axvspan(x0, x1, facecolor=RECAL, edgecolor="#b8bcc6", hatch="///", lw=0, zorder=0)
        ax.axvspan(xpos[open_i - 1] + 0.5, xpos[open_i] - 0.5, facecolor="white", lw=0, zorder=1)

    for ax in (ax1, ax2):
        shade(ax)
        ax.set_facecolor("white")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.axvline(xpos[open_i] - 0.5, color=INK2, lw=1, ls=(0, (2, 2)), zorder=2)

    # panel 1: the clipped, scaled squared return
    for i in range(lo, hi + 1):
        y = d["Y"][i]
        if y is None:
            continue
        monitored = not any(a <= i < b for a, b in d["recal"])
        ax1.bar(xpos[i], y, width=0.8, color=BAR if monitored else BAR_RECAL, zorder=3)
    cap = [s[1] for s in d["seg"] if s[0] <= open_i][-1]
    ax1.axhline(cap, color=INK2, lw=1, ls=":", zorder=2)
    ax1.text(xpos[hi] + 0.5, cap, f"cap B = {cap:.2f}", ha="right", va="bottom", color=INK2, fontsize=12.5)
    ax1.set_ylim(0, cap * 1.25)
    ax1.set_ylabel("clipped squared\nreturn / baseline", fontsize=13)
    ax1.yaxis.grid(True, color=GRID, lw=0.8)
    ax1.set_axisbelow(True)

    # panel 2: evidence (wealth) of the increase side, log scale, line broken while recalibrating
    seg_x, seg_y = [], []
    for i in range(lo, hi + 1):
        w = d["wu"][i]
        if i == open_i and seg_x:  # no line across the overnight break
            ax2.plot(seg_x, seg_y, color=EVIDENCE, lw=2.4, zorder=4, solid_capstyle="round")
            seg_x, seg_y = [], []
        if w is None:
            if seg_x:
                ax2.plot(seg_x, seg_y, color=EVIDENCE, lw=2.4, zorder=4, solid_capstyle="round")
            seg_x, seg_y = [], []
            continue
        seg_x.append(xpos[i])
        seg_y.append(w)
    if seg_x:
        ax2.plot(seg_x, seg_y, color=EVIDENCE, lw=2.4, zorder=4)
    ax2.set_yscale("log")
    ax2.set_ylim(0.3, 2e5)
    ax2.axhline(d["thr"], color=THRESH, lw=1.6, ls=(0, (6, 3)), zorder=3)
    ax2.text(xpos[lo] + 1, d["thr"] * 1.35, f"alarm threshold 2A = {d['thr']:,.0f}",
             color=THRESH, fontsize=14, va="bottom")
    ax2.plot([xpos[alarm_i]], [d["wu"][alarm_i]], "o", ms=9, mfc=THRESH, mec="white", mew=1.5, zorder=5)
    ax2.annotate("up alarm 09:26 IST\n11 minutes after the open", xy=(xpos[alarm_i], d["wu"][alarm_i]),
                 xytext=(xpos[alarm_i] + 3, 1.2e4 / 60), color=INK, fontsize=14, fontweight="semibold",
                 arrowprops={"arrowstyle": "-", "color": INK2, "lw": 1})
    ax2.text((xpos[alarm_i + 1] + xpos[hi]) / 2, 1.1, "recalibrating:\nnot monitoring",
             ha="center", va="center", color=INK2, fontsize=13,
             bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#b8bcc6", "lw": 0.8})
    first_resume = next(b for a, b in d["recal"] if a < lo < b)
    ax2.text((xpos[lo] + xpos[first_resume - 1]) / 2, 1.1, "recalibrating", ha="center", va="center",
             color=INK2, fontsize=13,
             bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#b8bcc6", "lw": 0.8})
    ax2.set_ylabel("evidence of a rise\n(wealth, log scale)", fontsize=13)
    ax2.yaxis.set_major_locator(FixedLocator([1, 10, 100, 1000, 10000, 100000]))
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax2.yaxis.grid(True, color=GRID, lw=0.8)
    ax2.set_axisbelow(True)

    ticks = ["2024-06-03 15:00", "2024-06-03 15:19", "2024-06-04 09:15", "2024-06-04 09:26",
             "2024-06-04 09:45", "2024-06-04 10:00"]
    labels = ["3 Jun\n15:00", "15:19", "4 Jun\n09:15 open", "09:26", "09:45", "10:00"]
    ax2.set_xticks([xpos[idx(s)] for s in ticks], labels)
    ax2.set_xlim(xpos[lo] - 1, xpos[hi] + 1)
    for ax in (ax1, ax2):
        ax.text((xpos[open_i - 1] + xpos[open_i]) / 2, 0.5, "overnight", transform=ax.get_xaxis_transform(),
                rotation=90, ha="center", va="center", color=INK2, fontsize=11)
    ax2.set_xlabel("India time (IST); one bar per minute; the overnight gap is shortened", fontsize=12.5)

    ax1.set_title("Nifty 50 on 4 June 2024, election results day",
                  loc="left", color=INK, fontsize=16, fontweight="semibold", pad=10)
    fig.text(0.01, 0.005, "REAL DATA, single historical case (E18). Source: docs/demo/replay.html, "
             "built by scripts/build_replay.py.", color=INK2, fontsize=11.5)
    fig.subplots_adjust(left=0.13, right=0.985, top=0.9, bottom=0.18)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, facecolor="white")
    print("wrote", OUT.relative_to(ROOT))
    plt.close(fig)
    phone(d, idx, alarm_i, open_i)


def phone(d: dict, idx, alarm_i: int, open_i: int) -> None:
    """The evidence panel alone, on a narrower window, with type sized for a 375 px screen."""
    lo, hi = idx("2024-06-03 15:19"), idx("2024-06-04 09:48")  # 15:19: monitoring resumed
    gap = 4
    xpos = {i: (i - lo) + (gap if i >= open_i else 0) for i in range(lo, hi + 1)}
    recal = [(max(a, lo), min(b, hi + 1)) for a, b in d["recal"] if b > lo and a <= hi]

    plt.rcParams.update({"font.family": "Inter", "font.size": 13, "axes.edgecolor": GRID,
                         "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
    fig, ax = plt.subplots(figsize=(5.4, 4.6), dpi=200)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for a, b in recal:
        ax.axvspan(xpos[a] - 0.5, xpos[b - 1] + 0.5, facecolor=RECAL, edgecolor="#b8bcc6", hatch="///",
                   lw=0, zorder=0)
    ax.axvline(xpos[open_i] - 0.5, color=INK2, lw=1, ls=(0, (2, 2)), zorder=2)

    seg_x, seg_y = [], []
    for i in range(lo, hi + 1):
        w = d["wu"][i]
        if (i == open_i or w is None) and seg_x:
            ax.plot(seg_x, seg_y, color=EVIDENCE, lw=2.6, zorder=4)
            seg_x, seg_y = [], []
        if w is None:
            continue
        seg_x.append(xpos[i])
        seg_y.append(w)
    if seg_x:
        ax.plot(seg_x, seg_y, color=EVIDENCE, lw=2.6, zorder=4)

    ax.set_yscale("log")
    ax.set_ylim(0.3, 3e5)
    ax.axhline(d["thr"], color=THRESH, lw=1.6, ls=(0, (6, 3)), zorder=3)
    ax.text(xpos[lo] + 0.3, d["thr"] * 1.5, f"alarm threshold\n2A = {d['thr']:,.0f}", color=THRESH,
            fontsize=13, va="bottom")
    ax.plot([xpos[alarm_i]], [d["wu"][alarm_i]], "o", ms=9, mfc=THRESH, mec="white", mew=1.5, zorder=5)
    ax.annotate("up alarm\n09:26 IST,\n11 min after\nthe open", xy=(xpos[alarm_i], d["wu"][alarm_i]),
                xytext=(xpos[alarm_i] + 3, 60), color=INK, fontsize=13, fontweight="semibold",
                arrowprops={"arrowstyle": "-", "color": INK2, "lw": 1})
    ax.text((xpos[alarm_i + 1] + xpos[hi]) / 2, 1.0, "recalibrating", ha="center", va="center",
            color=INK2, fontsize=12, bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#b8bcc6", "lw": 0.8})
    ax.text((xpos[open_i - 1] + xpos[open_i]) / 2 - 0.5, 0.5, "overnight", transform=ax.get_xaxis_transform(),
            rotation=90, ha="center", va="center", color=INK2, fontsize=11)
    ax.set_ylabel("evidence of a rise (log scale)", fontsize=12)
    ax.yaxis.set_major_locator(FixedLocator([1, 10, 100, 1000, 10000, 100000]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ticks = ["2024-06-03 15:19", "2024-06-04 09:15", "2024-06-04 09:45"]
    ax.set_xticks([xpos[idx(s)] for s in ticks], ["3 Jun\n15:19", "4 Jun\n09:15 open", "09:45"])
    ax.set_xlim(xpos[lo] - 1, xpos[hi] + 1)
    ax.set_xlabel("India time (IST), one point per minute", fontsize=12)
    ax.set_title("Nifty 50, 4 June 2024", loc="left", color=INK, fontsize=15, fontweight="semibold", pad=8)
    fig.text(0.02, 0.008, "Real data, single historical case (E18)", color=INK2, fontsize=11)
    fig.subplots_adjust(left=0.2, right=0.97, top=0.92, bottom=0.225)
    fig.savefig(OUT_PHONE, facecolor="white")
    print("wrote", OUT_PHONE.relative_to(ROOT))


if __name__ == "__main__":
    main()
