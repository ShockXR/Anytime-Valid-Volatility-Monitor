"""Experiments for the quarticity-stream e-detector (alternative A3).

Research questions
  Q1 [SYNTHETIC] Validity. With the model-free null mean q0 = mean of
     min(r^4, tau^4) and the bound tau^4, does the quarticity-stream detector keep
     P(false alarm within 390 steps) <= 390 / A and average run length >= A under
     Gaussian returns, jumps and Student t5 returns, on both sides?
  Q2 [SYNTHETIC] Pure variance change. How much slower is the quarticity stream
     than the variance stream for increases (sd x 1.25, 1.5, 2) and decreases
     (sd x 0.8, 0.67, 0.5)?
  Q3 [SYNTHETIC] Variance held fixed, shape changed. Gaussian to variance-matched
     Student t (df 3, 4, 5, 8), constant volatility to a variance-matched two-point
     mixture (vol of vol) and to slow stochastic volatility: which stream detects,
     in which direction, and how fast? A theory table gives the mean shift and
     signal-to-noise ratio of each stream for truncation at a = 3, 4, 5 sd.
  Q4 [REAL] Binance BTCUSDT and DOGEUSDT, 1-minute bars, 1 July to 31 December
     2025, deseasonalised. Alarm rates of the quarticity monitor against the
     variance monitor, and what the quarticity alarms coincide with: a variance
     alarm nearby, a level change, or a change of shape with the level fixed.

Settings follow E4 of the main report: per-step units, 2000-step calibration
window from the same process as the pre-change data, tau = 3 x scaled MAD,
A = 7800 for synthetic runs (5% false alarm budget per 390-step day) and
A = 28800 on real 1440-minute days (each side alarms at 2A). Seeds are fixed. Runs in about
4 minutes on two cores.
"""

from __future__ import annotations

import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import qstream as Q  # noqa: E402
from avmonitor.data import intraday_scale, load_binance  # noqa: E402
from avmonitor.detectors import summarize_alarms  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402
from avmonitor.robust import MAD_NORMAL_SCALE  # noqa: E402

OUT = HERE / "results"
A = 7800.0
CAL = 2000
REPS = 300
NULL_SCENARIOS = {
    "gaussian": dict(),
    "jumps": dict(p_jump=0.005, jump_scale=8.0),
    "student t5": dict(pre_df=5.0),
}
STREAMS = ("variance", "quarticity")


def _paths(reps, n, nu, kind, param, seed, **kw):
    """Calibration window and monitored returns from one process (jumps and tails in both)."""
    rng = np.random.default_rng(seed)
    r = Q.shape_change_returns(reps, CAL + n, None if nu is None else CAL + nu, kind, param, rng, **kw)
    return r[:, :CAL], r[:, CAL:]


def _run_all(cal, r, stream, direction, threshold=A):
    """Per-path baselines, then one vectorised detector run with the most conservative bound."""
    bases = [Q.calibrate_streams(cal[i], c=3.0) for i in range(r.shape[0])]
    tau = np.array([b.tau for b in bases])[:, None]
    if stream == "variance":
        y = np.minimum(r**2, tau**2) / np.array([b.m0 for b in bases])[:, None]
        ub = float(max(b.b_var for b in bases))
    else:
        y = np.minimum(r**4, tau**4) / np.array([b.q0 for b in bases])[:, None]
        ub = float(max(b.b_q for b in bases))
    from avmonitor.detectors import betting_edetector

    grid = Q.Q_UP_GRID if direction == "up" else None
    return betting_edetector(y, 1.0, direction, upper=ub, threshold=threshold, grid=grid).alarm


# ----------------------------------------------------------------------------
# Q1 validity
# ----------------------------------------------------------------------------


def q1_job(args):
    name, kw, seed = args
    n = 20000
    cal, r = _paths(REPS, n, None, "sd", 1.0, seed, **kw)
    rows = []
    for stream in STREAMS:
        for direction in ("up", "down"):
            a = _run_all(cal, r, stream, direction)
            rl = np.where(a >= 0, a + 1, n)
            rows.append(dict(scenario=name, stream=stream, direction=direction,
                             false_alarm_within_390=float(((a >= 0) & (a < 390)).mean()),
                             false_alarm_within_20000=float((a >= 0).mean()),
                             mean_run_length_censored=float(rl.mean()), bound_390=390 / A, A=A, paths=REPS))
    return rows


# ----------------------------------------------------------------------------
# Q2 pure variance change
# ----------------------------------------------------------------------------


def q2_job(args):
    name, kw, ratio, seed = args
    nu, n = 500, 3500
    cal, r = _paths(REPS, n, nu, "sd", ratio, seed, **kw)
    direction = "up" if ratio > 1 else "down"
    rows = []
    for stream in STREAMS:
        a = _run_all(cal, r, stream, direction)
        s = summarize_alarms(a, nu, n)
        rows.append(dict(scenario=name, sd_ratio=ratio, direction=direction, stream=stream, **s))
    return rows


# ----------------------------------------------------------------------------
# Q3 shape change with the variance fixed
# ----------------------------------------------------------------------------

SHAPE_ALTS = [("student", 3.0), ("student", 4.0), ("student", 5.0), ("student", 8.0),
              ("mixture", 0.5), ("mixture", 0.7), ("mixture", 0.9), ("sv", 0.3), ("sv", 0.5)]


A_Q3 = 400_000.0  # P(false alarm within the 20,000-step horizon) <= 5% under the null
HORIZON_Q3 = 20_000


def q3_job(args):
    kind, param, seed = args
    nu, reps = 500, 200
    n = nu + HORIZON_Q3
    cal, r = _paths(reps, n, nu, kind, param, seed)
    rows = []
    for stream in STREAMS:
        for direction in ("up", "down"):
            a = _run_all(cal, r, stream, direction, threshold=A_Q3)
            s = summarize_alarms(a, nu, n)
            rows.append(dict(kind=kind, param=param, stream=stream, direction=direction, A=A_Q3,
                             horizon_after_change=HORIZON_Q3, null_bound_within_horizon=HORIZON_Q3 / A_Q3, **s))
    return rows


def q3_theory():
    rows = []
    alts = [("sd", 1.5, Q.scaled_normal_pdf(1.5)), ("sd", 0.67, Q.scaled_normal_pdf(0.67))]
    alts += [("student", df, Q.student_pdf(df)) for df in (3.0, 4.0, 5.0, 8.0)]
    alts += [("mixture", d, Q.mixture_pdf(d)) for d in (0.5, 0.7, 0.9)]
    for a in (3.0, 4.0, 5.0):
        for kind, param, pdf in alts:
            rows.append(dict(cutoff_sd=a, kind=kind, param=param, **Q.stream_shift(pdf, a)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Q4 real data
# ----------------------------------------------------------------------------


def _deseasonalise(df):
    r = df.logret
    s = intraday_scale(r.fillna(0.0), lookback_days=20)
    u = (r / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    return u


def _post_ratios(u, e, c, cal_len, horizon=1440):
    """Truncated variance and quarticity in the next ``horizon`` minutes over the calibration window."""
    tau = c * e.baseline_scale
    cal = u[max(e.monitored_since - cal_len, 0): e.monitored_since]
    post = u[e.index + 1: e.index + 1 + horizon]
    cal, post = cal[np.isfinite(cal)], post[np.isfinite(post)]
    if post.size < 1000 or cal.size < 1000:
        return None
    m_c, m_p = np.mean(np.minimum(cal**2, tau**2)), np.mean(np.minimum(post**2, tau**2))
    q_c, q_p = np.mean(np.minimum(cal**4, tau**4)), np.mean(np.minimum(post**4, tau**4))
    pre = u[max(e.index - 60, 0): e.index + 1]
    pre = pre[np.isfinite(pre)]
    return dict(var_ratio=m_p / m_c, q_ratio=q_p / q_c, shape_ratio=(q_p / q_c) / (m_p / m_c) ** 2,
                clipped_in_last_hour=int(np.sum(np.abs(pre) >= tau)))


def q4_real():
    A_real, cal_len, c = 28800.0, 1440, 3.0
    summary, events, top_days = [], [], []
    runs = [("variance monitor", VolatilityMonitor, 0.25), ("quarticity monitor", Q.QuarticityMonitor, 0.25),
            ("quarticity monitor", Q.QuarticityMonitor, 0.5625)]
    for symbol in ("BTCUSDT", "DOGEUSDT"):
        df = load_binance(symbol, "2025-07-01", "2026-01-01", download=False)
        u = _deseasonalise(df)
        days = np.sum(np.isfinite(u)) / 1440
        stream = np.nan_to_num(u, nan=0.0)
        var_times = None
        for name, cls, band in runs:
            mon = cls(threshold=A_real, c=c, cal_len=cal_len, null_band=band)
            ev = mon.run(stream)
            idx = np.array([e.index for e in ev])
            if name == "variance monitor":
                var_times = idx
            rows = []
            for e in ev:
                pr = _post_ratios(u, e, c, cal_len)
                if pr is None:
                    continue
                near = int(np.min(np.abs(var_times - e.index))) if var_times is not None and var_times.size else -1
                rows.append(dict(symbol=symbol, monitor=name, band=band, time=df.index[e.index], direction=e.direction,
                                 delay_minutes=e.index - e.monitored_since, minutes_to_nearest_variance_alarm=near, **pr))
            ed = pd.DataFrame(rows)
            events.append(ed)
            if name == "variance monitor":
                var_days = set(ed.time.dt.floor("D"))
            up, dn = ed[ed.direction == "up"], ed[ed.direction == "down"]
            far = ed[ed.minutes_to_nearest_variance_alarm > 240]
            summary.append(dict(
                symbol=symbol, monitor=name, band=band, alarms=len(ev), alarms_per_day=len(ev) / days,
                share_same_utc_day_as_a_variance_alarm=float(ed.time.dt.floor("D").isin(var_days).mean()),
                alarms_more_than_240min_from_variance_alarm=len(far),
                far_alarms_share_up=float((far.direction == "up").mean()) if len(far) else np.nan,
                far_alarms_median_var_ratio=float(far.var_ratio.median()) if len(far) else np.nan,
                far_alarms_median_q_ratio=float(far.q_ratio.median()) if len(far) else np.nan,
                far_alarms_share_followed_in_direction=float(np.mean(np.where(far.direction == "up", far.var_ratio > 1, far.var_ratio < 1))) if len(far) else np.nan,
                up=int((np.array([e.direction for e in ev]) == "up").sum()),
                down=int((np.array([e.direction for e in ev]) == "down").sum()),
                median_var_ratio_up=up.var_ratio.median(), median_q_ratio_up=up.q_ratio.median(),
                share_up_followed_by_higher_var=(up.var_ratio > 1).mean(),
                share_up_followed_by_higher_q=(up.q_ratio > 1).mean(),
                median_var_ratio_down=dn.var_ratio.median(), median_q_ratio_down=dn.q_ratio.median(),
                share_down_followed_by_lower_var=(dn.var_ratio < 1).mean(),
                share_down_followed_by_lower_q=(dn.q_ratio < 1).mean(),
                median_shape_ratio_up=up.shape_ratio.median(), median_shape_ratio_down=dn.shape_ratio.median(),
                share_within_60min_of_variance_alarm=(ed.minutes_to_nearest_variance_alarm <= 60).mean(),
                share_within_240min_of_variance_alarm=(ed.minutes_to_nearest_variance_alarm <= 240).mean(),
                share_up_with_clipped_return_in_last_hour=(up.clipped_in_last_hour > 0).mean()))
            per_month = ed.groupby(ed.time.dt.to_period("M")).size()
            for month, cnt in per_month.items():
                top_days.append(dict(symbol=symbol, monitor=name, band=band, month=str(month), alarms=int(cnt)))
        # control: ratios at random times
        rng = np.random.default_rng(0)
        ctrl = []
        for t in rng.integers(cal_len, len(u) - 1441, 2000):
            a, b = u[t - cal_len: t], u[t + 1: t + 1441]
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            tau = c * np.median(np.abs(a)) * MAD_NORMAL_SCALE
            ctrl.append((np.mean(np.minimum(b**2, tau**2)) / np.mean(np.minimum(a**2, tau**2)),
                         np.mean(np.minimum(b**4, tau**4)) / np.mean(np.minimum(a**4, tau**4))))
        ctrl = np.array(ctrl)
        summary.append(dict(symbol=symbol, monitor="control: 2000 random times", band=np.nan, alarms=np.nan,
                            alarms_per_day=np.nan, share_same_utc_day_as_a_variance_alarm=np.nan,
                            alarms_more_than_240min_from_variance_alarm=np.nan, far_alarms_share_up=np.nan,
                            far_alarms_median_var_ratio=np.nan, far_alarms_median_q_ratio=np.nan,
                            far_alarms_share_followed_in_direction=np.nan, up=np.nan, down=np.nan,
                            median_var_ratio_up=np.median(ctrl[:, 0]), median_q_ratio_up=np.median(ctrl[:, 1]),
                            share_up_followed_by_higher_var=np.mean(ctrl[:, 0] > 1),
                            share_up_followed_by_higher_q=np.mean(ctrl[:, 1] > 1),
                            median_var_ratio_down=np.median(ctrl[:, 0]), median_q_ratio_down=np.median(ctrl[:, 1]),
                            share_down_followed_by_lower_var=np.mean(ctrl[:, 0] < 1),
                            share_down_followed_by_lower_q=np.mean(ctrl[:, 1] < 1),
                            median_shape_ratio_up=np.median(ctrl[:, 1] / ctrl[:, 0] ** 2),
                            median_shape_ratio_down=np.median(ctrl[:, 1] / ctrl[:, 0] ** 2),
                            share_within_60min_of_variance_alarm=np.nan, share_within_240min_of_variance_alarm=np.nan,
                            share_up_with_clipped_return_in_last_hour=np.nan))
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "q4_real_summary.csv", index=False)
    pd.concat(events).to_csv(OUT / "q4_real_events.csv", index=False)
    pd.DataFrame(top_days).to_csv(OUT / "q4_real_alarms_per_month.csv", index=False)
    print(s.round(3).to_string(index=False))
    return s


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

COL = {"variance": "#2a78d6", "quarticity": "#eb6834"}


def figures(q2: pd.DataFrame, q3: pd.DataFrame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    g = q2[q2.scenario == "gaussian"]
    xs = np.arange(6)
    labels = [f"x{r:g}" for r in g.sd_ratio.unique()]
    for j, stream in enumerate(STREAMS):
        d = g[g.stream == stream].sort_values("sd_ratio", ascending=False)
        ax[0].bar(xs + (j - 0.5) * 0.38, d.median_delay, width=0.36, color=COL[stream], label=stream + " stream")
    ax[0].set_xticks(xs, [f"x{r:g}" for r in g[g.stream == "variance"].sort_values("sd_ratio", ascending=False).sd_ratio])
    ax[0].set_xlabel("sd ratio after the change (Gaussian returns)")
    ax[0].set_ylabel("median detection delay (steps)")
    ax[0].set_title("Q2: pure variance change")
    ax[0].legend(frameon=False)
    for s in ("top", "right"):
        ax[0].spines[s].set_visible(False)
    d3 = q3[q3.direction == "up"].copy()
    d3["label"] = d3.kind + " " + d3.param.map(lambda v: f"{v:g}")
    order = [f"{k} {p:g}" for k, p in SHAPE_ALTS]
    xs = np.arange(len(order))
    for j, stream in enumerate(STREAMS):
        d = d3[d3.stream == stream].set_index("label").reindex(order)
        ax[1].bar(xs + (j - 0.5) * 0.38, d.detected, width=0.36, color=COL[stream], label=stream + " stream, up side")
    d3d = q3[(q3.direction == "down") & (q3.stream == "variance")].copy()
    d3d["label"] = d3d.kind + " " + d3d.param.map(lambda v: f"{v:g}")
    d = d3d.set_index("label").reindex(order)
    ax[1].plot(xs, d.detected, "v", color=COL["variance"], ms=7, label="variance stream, down side")
    ax[1].set_xticks(xs, order, rotation=45, ha="right")
    ax[1].set_ylabel(f"share detected within {HORIZON_Q3:,} steps")
    ax[1].set_title(f"Q3: variance fixed, shape changed (A = {A_Q3:,.0f})")
    ax[1].set_ylim(0, 1.45)
    ax[1].legend(frameon=False, fontsize=8, loc="upper center", ncol=3)
    for s in ("top", "right"):
        ax[1].spines[s].set_visible(False)
    fig.suptitle("Synthetic: variance stream against quarticity stream (A = 7800, tau = 3 sd)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_synthetic_delays.png", dpi=130)
    plt.close(fig)


def figure_real():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ev = pd.read_csv(OUT / "q4_real_events.csv", parse_dates=["time"])
    df = load_binance("DOGEUSDT", "2025-07-01", "2026-01-01", download=False)
    hourly = df.logret.pow(2).resample("1h").sum().pow(0.5)
    fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    ax[0].plot(hourly, lw=0.6, color="#898781")
    ax[0].set_ylabel("DOGEUSDT hourly realized vol")
    for name, a in (("variance monitor", ax[0]), ("quarticity monitor", ax[1])):
        e = ev[(ev.symbol == "DOGEUSDT") & (ev.monitor == name) & (ev.band == 0.25)]
        for _, row in e.iterrows():
            a.axvline(row.time, color=COL["variance" if name.startswith("variance") else "quarticity"],
                      lw=0.7, alpha=0.7, ymin=0.0 if row.direction == "down" else 0.5,
                      ymax=0.5 if row.direction == "down" else 1.0)
        a.set_title(f"{name} alarms (band 0.25): upper half = increase, lower half = decrease", fontsize=9)
    ax[1].plot(hourly, lw=0.6, color="#898781")
    ax[1].set_ylabel("DOGEUSDT hourly realized vol")
    fig.suptitle("Real data: Binance DOGEUSDT, 1 July to 31 December 2025", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_real_doge_alarms.png", dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------


def main(parts=("q1", "q2", "q3", "q4")):
    """``python run_experiment.py [q1 q2 q3 q4]`` runs the named parts (default all) and redraws the figures."""
    OUT.mkdir(exist_ok=True)
    with Pool(2) as pool:
        if "q1" in parts:
            q1 = pd.DataFrame(sum(pool.map(q1_job, [(nm, kw, 10 + i) for i, (nm, kw) in enumerate(NULL_SCENARIOS.items())]), []))
            q1.to_csv(OUT / "q1_validity.csv", index=False)
            print(q1.round(4).to_string(index=False))
        if "q2" in parts:
            jobs = []
            for i, (nm, kw) in enumerate([("gaussian", {}), ("jumps", NULL_SCENARIOS["jumps"])]):
                for j, ratio in enumerate((1.25, 1.5, 2.0, 0.8, 0.67, 0.5)):
                    jobs.append((nm, kw, ratio, 100 + 10 * i + j))
            q2 = pd.DataFrame(sum(pool.map(q2_job, jobs), []))
            q2.to_csv(OUT / "q2_variance_change.csv", index=False)
            print(q2.round(3).to_string(index=False))
        if "q3" in parts:
            q3 = pd.DataFrame(sum(pool.map(q3_job, [(k, p, 200 + i) for i, (k, p) in enumerate(SHAPE_ALTS)]), []))
            q3.to_csv(OUT / "q3_shape_change.csv", index=False)
            print(q3.round(3).to_string(index=False))
            th = q3_theory()
            th.to_csv(OUT / "q3_theory.csv", index=False)
            print(th.round(4).to_string(index=False))
    figures(pd.read_csv(OUT / "q2_variance_change.csv"), pd.read_csv(OUT / "q3_shape_change.csv"))
    if "q4" in parts:
        q4_real()
    figure_real()


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("q1", "q2", "q3", "q4"))
