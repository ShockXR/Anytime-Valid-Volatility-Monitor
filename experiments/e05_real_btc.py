"""E5. Real data: Bitstamp BTC/USD 1-minute bars, 1 July to 31 December 2025.

REAL MARKET DATA. Source: https://github.com/ff137/bitstamp-btcusd-minute-data
(Bitstamp only; Bapat 2026 does not name the venue of its crypto data).

Research questions
  Q5a. Which formulas reproduce the magnitudes in Table 6 (BTC column)? The
       printed RQ, BPQ, TPQ and kernel formulas imply different scales, so the
       published numbers can be used to infer how each was computed.
  Q5b. On a common, documented scale (daily IQ with T = 1 day), how do the
       estimators compare on real data?
  Q5c. Run the streaming monitor on deseasonalised returns. How often does it
       alarm, and are its alarms followed by a persistent change in volatility
       (unlike alarms of the Gaussian benchmark)?
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from paper_tables import TABLE5_BTC, TABLE6_BTC  # noqa: E402

from avmonitor import robust as R  # noqa: E402
from avmonitor.data import intraday_scale, load_bitstamp_csv  # noqa: E402
from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
RAW = ROOT / "data/raw/crypto/btcusd_bitstamp_1min_latest.csv"


def q5a(r: np.ndarray):
    n = r.size
    a = np.abs(r)
    mu = R._mu
    cands = {
        "RQ": {
            "(n/3) sum r^4, n = whole sample": n / 3 * np.sum(r**4),
            "mean over days of (n_d/3) sum_day r^4": None,
            "(1/3) sum r^4 (no n)": np.sum(r**4) / 3,
        },
        "BPQ": {
            "eq. (4) as printed: mu1^-4 sum |r_i||r_{i-1}|^3": mu(1) ** -4 * np.sum(a[1:] * a[:-1] ** 3),
            "mu1^-2 sum |r_i||r_{i-1}| (bipower variation)": mu(1) ** -2 * np.sum(a[1:] * a[:-1]),
        },
        "TPQ": {
            f"eq. (5) as printed, p = {p:.3g}": mu(p) ** -3 * np.sum(a[2:] ** p * a[1:-1] ** p * a[:-2] ** p)
            for p in (0.5, 2 / 3, 4 / 3)
        },
        "MinRQ": {"ADS constant, (n) scaling": R.min_rq(r), "ADS constant, no n": R.min_rq(r) / n},
        "MedRQ": {"ADS constant, (n) scaling": R.med_rq(r), "ADS constant, no n": R.med_rq(r) / n},
        "Kernel": {"eq. (7) as printed, h = 30": R.kernel_quarticity_printed(r, 30)},
        "Robust kernel": {
            "eqs. (11), (16) as printed, raw MAD": R.robust_kernel_quarticity_printed(r, 30, 3.0, mad_scaled=False),
            "eqs. (11), (16) as printed, scaled MAD": R.robust_kernel_quarticity_printed(r, 30, 3.0, mad_scaled=True),
        },
    }
    rows = []
    for est, d in cands.items():
        for formula, val in d.items():
            if val is None:
                continue
            rows.append(dict(estimator=est, formula=formula, ours=float(val), published=TABLE6_BTC[est],
                             ratio=float(val) / TABLE6_BTC[est]))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e05a_table6_btc.csv", index=False)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    q = np.quantile(r, [0, 0.25, 0.5, 0.75, 1])
    print("Table 5 check (ours vs published):", dict(zip(["min", "q1", "median", "q3", "max"], np.round(q, 6))),
          "mean", r.mean(), "| published", TABLE5_BTC)


def q5b(df: pd.DataFrame):
    rows = []
    for day, g in df.groupby(df.index.floor("D")):
        r = g.logret.dropna().to_numpy()
        if r.size < 1400:
            continue
        rows.append(dict(day=day, RQ=R.realized_quarticity(r), TPQ=R.tripower_quarticity(r), QPQ=R.quadpower_quarticity(r),
                         MinRQ=R.min_rq(r), MedRQ=R.med_rq(r), TRQ=R.threshold_rq(r), Kernel=R.kernel_iq(r, 30),
                         RobustKernel_scaled_bc=R.kernel_iq(r, 30, c=3, mad_scaled=True, bias_correct=True),
                         RobustKernel_raw=R.kernel_iq(r, 30, c=3, mad_scaled=False),
                         RV=np.sum(r**2), zero_share=float(np.mean(r == 0))))
    d = pd.DataFrame(rows).set_index("day")
    d.to_csv(OUT / "e05b_daily_iq.csv")
    cols = ["RQ", "TPQ", "QPQ", "MinRQ", "MedRQ", "TRQ", "Kernel", "RobustKernel_scaled_bc", "RobustKernel_raw"]
    summ = pd.DataFrame({"median_daily_IQ": d[cols].median(), "median_ratio_to_MedRQ": (d[cols].div(d.MedRQ, axis=0)).median(),
                         "share_days_ratio_RQ_over_TRQ_gt_2": [(d.RQ / d.TRQ > 2).mean()] * len(cols)})
    print(summ.to_string(float_format=lambda v: f"{v:.3g}"))
    print("Median IQ / RV^2 (equals 1 for constant volatility):", float((d.MedRQ / d.RV**2).median()))


def q5c(df: pd.DataFrame):
    r = df.logret
    s = intraday_scale(r.fillna(0.0), lookback_days=20)
    u = (r / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    results, events_all = [], []
    A = 28800.0  # 1440 / A = 5% chance of any false alarm per day under the null
    for name, cls, band in [("robust e-detector", VolatilityMonitor, 0.0), ("robust e-detector", VolatilityMonitor, 0.25),
                            ("robust e-detector", VolatilityMonitor, 0.5), ("Gaussian LR", GaussianLRMonitor, 0.0)]:
        mon = cls(threshold=A, c=3.0, cal_len=1440, null_band=band)
        ev = mon.run(np.nan_to_num(u, nan=0.0))
        # validation: robust sd in the 1440 minutes after each alarm relative to its baseline
        post = []
        for e in ev:
            seg = u[e.index + 1 : e.index + 1441]
            seg = seg[np.isfinite(seg)]
            if seg.size < 1000:
                continue
            post.append(np.median(np.abs(seg)) * R.MAD_NORMAL_SCALE / e.baseline_scale)
            events_all.append(dict(detector=name, band=band, time=df.index[e.index], direction=e.direction,
                                   delay_minutes=e.index - e.monitored_since, post_over_baseline=post[-1]))
        post = np.array(post)
        dirs = np.array([e.direction for e in ev][: len(post)])
        days = len(u) / 1440
        res = dict(detector=name, band=band, alarms=len(ev), alarms_per_day=len(ev) / days,
                   up=int(np.sum(dirs == "up")), down=int(np.sum(dirs == "down")),
                   median_post_ratio_up=float(np.median(post[dirs == "up"])) if np.any(dirs == "up") else np.nan,
                   median_post_ratio_down=float(np.median(post[dirs == "down"])) if np.any(dirs == "down") else np.nan,
                   share_up_followed_by_higher_vol=float(np.mean(post[dirs == "up"] > 1)) if np.any(dirs == "up") else np.nan,
                   share_down_followed_by_lower_vol=float(np.mean(post[dirs == "down"] < 1)) if np.any(dirs == "down") else np.nan)
        results.append(res)
    # Control: random times, next-day robust sd over previous-day robust sd
    rng = np.random.default_rng(0)
    ctrl = []
    for t in rng.integers(1440, len(u) - 1441, 2000):
        a, b = u[t - 1440 : t], u[t + 1 : t + 1441]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        ctrl.append(np.median(np.abs(b)) / np.median(np.abs(a)))
    ctrl = np.array(ctrl)
    results.append(dict(detector="control: 2000 random times", band=np.nan, alarms=np.nan, alarms_per_day=np.nan,
                        up=np.nan, down=np.nan, median_post_ratio_up=float(np.median(ctrl)),
                        median_post_ratio_down=float(np.median(ctrl)),
                        share_up_followed_by_higher_vol=float(np.mean(ctrl > 1)),
                        share_down_followed_by_lower_vol=float(np.mean(ctrl < 1))))
    res = pd.DataFrame(results)
    res.to_csv(OUT / "e05c_monitor_summary.csv", index=False)
    pd.DataFrame(events_all).to_csv(OUT / "e05c_events.csv", index=False)
    print(res.round(3).to_string(index=False))
    return u, s


def figure(df, u):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ev = pd.read_csv(OUT / "e05c_events.csv", parse_dates=["time"])
    ev = ev[(ev.detector == "robust e-detector") & (ev.band == 0.25)]
    r = df.logret
    hourly_vol = r.pow(2).resample("1h").sum().pow(0.5)
    fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    ax[0].plot(df.close.resample("1h").last(), lw=0.8, color="#333")
    ax[0].set_ylabel("BTC/USD (Bitstamp)")
    ax[1].plot(hourly_vol, lw=0.6, color="#4a6fa5")
    ax[1].set_ylabel("hourly realized vol")
    for _, e in ev.iterrows():
        col = "#c0392b" if e.direction == "up" else "#27ae60"
        for a in ax:
            a.axvline(e.time, color=col, lw=0.6, alpha=0.6)
    ax[0].set_title("Robust e-detector alarms (band 0.25): red = increase, green = decrease. Real data, Jul to Dec 2025")
    fig.tight_layout()
    fig.savefig(OUT / "figures/e05_btc_alarms.png", dpi=130)


def main():
    df = load_bitstamp_csv(RAW, "2025-07-01", "2026-01-01")
    r = df.logret.dropna().to_numpy()
    print(f"REAL DATA: {len(df)} minutes, zero-return share {np.mean(r == 0):.3f}")
    q5a(r)
    q5b(df)
    u, _ = q5c(df)
    figure(df, u)


if __name__ == "__main__":
    main()
