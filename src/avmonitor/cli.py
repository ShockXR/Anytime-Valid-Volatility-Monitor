"""Command-line demonstration of the volatility monitor.

Examples::

    python -m avmonitor simulate                       # synthetic path, change up then down
    python -m avmonitor run --source binance --symbol ETHUSDT --start 2025-10-01 --end 2025-11-01
    python -m avmonitor run --source nifty --start 2024-05-01 --end 2024-07-01 --plot nifty.png
    python -m avmonitor estimators --source binance --symbol BTCUSDT --day 2025-10-10

``run`` and ``estimators`` need the data files (``python scripts/download_data.py``).
Every output line says whether the data are synthetic or real.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from . import robust as R
from .monitor import GaussianLRMonitor, VolatilityMonitor

STEPS_PER_DAY = {"binance": 1440, "bitstamp": 1440, "nifty": 375}


def _threshold(steps_per_day: int, budget: float) -> float:
    """Alarm level A such that P(any false alarm within one day) <= budget (docs/math.md, Theorem B)."""
    return steps_per_day / budget


def _events_table(events, index, u, steps_per_day) -> pd.DataFrame:
    rows = []
    for e in events:
        seg = u[e.index + 1 : e.index + 1 + steps_per_day]
        seg = seg[np.isfinite(seg)]
        post = np.median(np.abs(seg)) * R.MAD_NORMAL_SCALE / e.baseline_scale if seg.size > steps_per_day // 2 else np.nan
        rows.append(dict(time=index[e.index] if index is not None else e.index, direction=e.direction,
                         minutes_since_calibration=e.index - e.monitored_since, next_day_vol_over_baseline=post))
    return pd.DataFrame(rows)


def cmd_simulate(a) -> int:
    rng = np.random.default_rng(a.seed)
    n = 9000
    sd = np.ones(n)
    sd[3000:6000] = a.ratio
    r = sd * rng.standard_normal(n)
    r += (rng.random(n) < a.jump_prob) * rng.normal(0, 10, n)
    A = _threshold(390, a.budget)
    rob = VolatilityMonitor(threshold=A, cal_len=1000).run(r)
    gau = GaussianLRMonitor(threshold=A, cal_len=1000).run(r)
    print(f"[SYNTHETIC] sd x{a.ratio} at step 3000, back to 1 at step 6000, jumps with prob {a.jump_prob} per step")
    print(f"alarm level A = {A:.0f}: at most {a.budget:.0%} chance of a false alarm per 390 steps under the null")
    print("robust e-detector:", [(e.index, e.direction) for e in rob])
    print("Gaussian LR      :", [(e.index, e.direction) for e in gau])
    return 0


def _load(a):
    from . import data as D

    if a.source == "binance":
        df = D.load_binance(a.symbol, a.start, a.end, download=a.download)
    elif a.source == "nifty":
        df = D.regular_sessions(D.load_nifty(a.start, a.end, download=a.download))
    elif a.source == "bitstamp":
        df = D.load_bitstamp_csv(D.DATA_ROOT / "crypto" / "btcusd_bitstamp_1min_latest.csv", a.start, a.end)
    else:
        raise SystemExit(f"unknown source {a.source}")
    return df


def cmd_run(a) -> int:
    from .data import intraday_scale

    df = _load(a)
    spd = STEPS_PER_DAY[a.source]
    r = df["logret"]
    s = intraday_scale(r.fillna(0.0), lookback_days=a.lookback_days)
    u = (r / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    A = _threshold(spd, a.budget)
    cls = GaussianLRMonitor if a.detector == "gaussian" else VolatilityMonitor
    mon = cls(threshold=A, c=a.c, cal_len=a.cal_len or spd, null_band=a.band)
    # gaps: a NaN return is skipped by the monitor (no information), which is the same as dropping it
    ev = mon.run(np.where(np.isfinite(u), u, np.nan))
    label = f"{a.source} {a.symbol if a.source == 'binance' else ''}".strip()
    print(f"[REAL] {label}, {df.index[0]} to {df.index[-1]}, {len(df)} bars")
    print(f"detector {a.detector}, A = {A:.0f} (false alarm budget {a.budget:.0%} per day), band {a.band}")
    tab = _events_table(ev, df.index, u, spd)
    days = len(df) / spd
    print(f"{len(ev)} alarms in {days:.1f} days ({len(ev) / max(days, 1e-9):.2f} per day)")
    if len(tab):
        with pd.option_context("display.width", 140, "display.max_rows", a.max_rows):
            print(tab.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    if a.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(2, 1, figsize=(11, 5.5), sharex=True)
        ax[0].plot(df.index, df["close"], lw=0.6, color="#333")
        ax[0].set_ylabel("price")
        hv = r.pow(2).rolling(spd // 24 if spd > 400 else 15, min_periods=5).sum().pow(0.5)
        ax[1].plot(df.index, hv, lw=0.5, color="#4a6fa5")
        ax[1].set_ylabel("rolling realized vol")
        for e in ev:
            col = "#c0392b" if e.direction == "up" else "#27ae60"
            for x in ax:
                x.axvline(df.index[e.index], color=col, lw=0.7, alpha=0.7)
        ax[0].set_title(f"[REAL] {label}: alarms (red = volatility up, green = down)")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=120)
        print("saved", a.plot)
    return 0


def cmd_estimators(a) -> int:
    a.start = a.day
    a.end = str((pd.Timestamp(a.day) + pd.Timedelta(days=1)).date())
    df = _load(a)
    r = df["logret"].dropna().to_numpy()
    with np.errstate(all="ignore"):
        rows = {
            "RQ": R.realized_quarticity(r), "TPQ": R.tripower_quarticity(r), "QPQ": R.quadpower_quarticity(r),
            "MinRQ": R.min_rq(r), "MedRQ": R.med_rq(r), "TRQ": R.threshold_rq(r),
            "Kernel (normalised, h=30)": R.kernel_iq(r, 30),
            "Robust kernel (scaled MAD, c=3)": R.kernel_iq(r, 30, c=3, mad_scaled=True),
            "Robust kernel (scaled MAD, c=3, bias corrected)": R.kernel_iq(r, 30, c=3, mad_scaled=True, bias_correct=True),
            "Robust kernel as printed (raw MAD, eq. 16)": R.robust_kernel_quarticity_printed(r, 30, 3.0),
        }
    print(f"[REAL] {a.source} {a.symbol if a.source == 'binance' else ''} {a.day}: {r.size} returns, "
          f"integrated quarticity with T = 1 day (the printed eq. 16 is on its own scale, see docs/estimator_issues.md)")
    for k, v in rows.items():
        print(f"  {k:<50s} {v:.4g}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m avmonitor", description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("simulate", help="synthetic demo with a known change")
    s.add_argument("--ratio", type=float, default=1.5, help="sd multiplier between steps 3000 and 6000")
    s.add_argument("--jump-prob", type=float, default=0.003)
    s.add_argument("--budget", type=float, default=0.05, help="false alarm budget per 390 steps")
    s.add_argument("--seed", type=int, default=42)
    s.set_defaults(func=cmd_simulate)

    for name, fn, hlp in (("run", cmd_run, "run the monitor on real data"),
                          ("estimators", cmd_estimators, "all quarticity estimators for one day")):
        q = sub.add_parser(name, help=hlp)
        q.add_argument("--source", choices=sorted(STEPS_PER_DAY), default="binance")
        q.add_argument("--symbol", default="BTCUSDT", help="Binance symbol, for example ETHUSDT")
        q.add_argument("--download", action="store_true", help="fetch missing files from the public mirrors")
        if name == "run":
            q.add_argument("--start", required=True)
            q.add_argument("--end", required=True)
            q.add_argument("--band", type=float, default=0.25, help="null band (0.25 is the real-data default)")
            q.add_argument("--budget", type=float, default=0.05, help="false alarm budget per day")
            q.add_argument("--c", type=float, default=3.0, help="truncation at c times the scaled MAD")
            q.add_argument("--cal-len", type=int, default=0, help="calibration length (default: one day)")
            q.add_argument("--lookback-days", type=int, default=20, help="days in the intraday profile")
            q.add_argument("--detector", choices=["robust", "gaussian"], default="robust")
            q.add_argument("--plot", default="", help="save a figure to this path")
            q.add_argument("--max-rows", type=int, default=40)
        else:
            q.add_argument("--day", required=True)
        q.set_defaults(func=fn)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
