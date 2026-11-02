"""E10. Forecast-relative monitoring: run the monitor on returns divided by a volatility forecast.

Motivation (E9, Q9c): under persistent volatility clustering the conditional
mean of the truncated squared return wanders, so the raw-return monitor alarms
on most days even though the unconditional variance never changes. Dividing
each return by a predictable volatility forecast restates the null as "the
forecast is right up to a constant". With the true conditional volatility this
null is exactly true for GARCH returns, so the guarantee applies again.

Research questions
  Q10a [SYNTHETIC]. Under a GARCH(1,1) null with no structural change, what is
        the one-day false alarm rate of the monitor on (i) raw returns,
        (ii) returns divided by the true conditional volatility (oracle),
        (iii) an EWMA forecast (half-lives 60 and 360 steps), and (iv) a GARCH
        forecast with parameters fitted on a separate training sample?
  Q10b [SYNTHETIC]. Power. A structural break doubles the variance level 200
        steps into the first monitored day, either as a scale shift (all
        returns multiplied by sqrt 2) or as a doubling of the GARCH intercept
        omega (the long-run variance doubles, reached gradually). Monitoring
        continues for five days. How often, and how fast, does each version
        detect it? A forecast that reacts to recent returns absorbs part of a
        break, so a trade-off between validity and power is expected.
  Q10b2 [SYNTHETIC]. Control for Q10b. The same monitors on paths with no
        break, over the same five days, so that "detected within 5 days" can be
        compared with how often the monitor alarms anyway. Also a larger break
        (long-run variance times 4), to test the signal cap of Q10d.
  Q10d [THEORY + SYNTHETIC]. How much of a break survives standardisation by a
        GARCH forecast with frozen parameters? Multiplying GARCH(omega, a, b)
        returns by c gives exactly GARCH(c^2 omega, a, b), so after the
        transient both kinds of break are the same model. In its stationary
        regime, the mean of r_t^2 / s2_t (s2_t the pre-break filter) is
        rho(k) = k (1 - b) / (1 - a - b + a k) for a variance factor k (unit
        pre-break variance), which never exceeds (1 - b) / a. We check the
        formula by simulation and compute the Kullback-Leibler rate between the
        post-break and pre-break conditional laws, which bounds how fast ANY
        detector with the same false alarm budget can react (Lai 1998).
  Q10e [SYNTHETIC]. Where in time is the evidence? After a break at time 0,
        the mean of s2_post / s2_pre along the path and the cumulative
        Kullback-Leibler divergence, and the step at which the average
        cumulative divergence reaches ln A (a rough predicted delay).
  Q10c [REAL]. On the five Binance coins (Jul to Dec 2025, deseasonalised
        one-minute returns), how do alarm rates change when the monitor runs on
        returns divided by an EWMA forecast, and after an "up" ("down") alarm
        does realised variance over the next day exceed (fall short of) the
        forecast made at the alarm, more often than at random times?
"""

from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import forecast as F
from avmonitor.data import intraday_scale, load_binance
from avmonitor.monitor import VolatilityMonitor

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
SEED = 1003
DAY = 1440
A = DAY / 0.05


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return float(c - h), float(c + h)


def garch_paths(reps, n, omega, a, b, rng, burn=3000, break_at=None, break_kind=None, factor=2.0):
    """GARCH(1,1) returns and true conditional variances, shape (reps, n), after a burn-in.

    A break at ``break_at`` multiplies the long-run variance by ``factor``:
    ``"omega"`` multiplies the intercept (the variance moves there gradually),
    ``"scale"`` multiplies every later return by ``sqrt(factor)`` (at once).
    """
    s2 = np.full(reps, omega / (1 - a - b))
    r = np.empty((reps, n))
    v = np.empty((reps, n))
    for t in range(-burn, n):
        om = omega
        if break_at is not None and t >= break_at and break_kind == "omega":
            om = factor * omega
        x = np.sqrt(s2) * rng.standard_normal(reps)
        scale = np.sqrt(factor) if (break_at is not None and t >= break_at and break_kind == "scale") else 1.0
        if t >= 0:
            r[:, t] = scale * x
            v[:, t] = s2 * scale**2
        s2 = om + a * x * x + b * s2
    return r, v


def first_event(stream, band, cal_len=DAY):
    mon = VolatilityMonitor(threshold=A, cal_len=cal_len, null_band=band)
    for i, x in enumerate(stream):
        ev = mon.update(x)
        if ev is not None:
            return ev
    return None


def standardised_versions(r, v_true, params_fit, pre):
    """Return dict name -> standardised series for one path; ``pre`` is the warm-up length to drop."""
    out = {"raw": r[pre:], "oracle": (r / np.sqrt(v_true))[pre:]}
    for hl in (60, 360):
        out[f"ewma{hl}"] = F.standardise(r, F.ewma_variance(r, hl, init=float(np.mean(r[:pre] ** 2))))[pre:]
    om, a, b = params_fit
    out["garch_fit"] = F.standardise(r, F.garch11_filter(r, om, a, b))[pre:]
    return out


SETTINGS = [(0.05, 0.90), (0.05, 0.94), (0.03, 0.969)]


def q10a(reps=600):
    rng = np.random.default_rng(SEED)
    rows = []
    pre = 2000  # warm-up for the forecasts, not monitored
    for a, b in SETTINGS:
        omega = 1 - a - b
        train, _ = garch_paths(1, 30000, omega, a, b, rng)
        fit = F.garch11_fit(train[0])
        r, v = garch_paths(reps, pre + 2 * DAY, omega, a, b, rng)
        counts = {}
        for i in range(reps):
            vers = standardised_versions(r[i], v[i], fit, pre)
            for name, z in vers.items():
                for band in (0.0, 0.25):
                    ev = first_event(z, band)
                    counts.setdefault((name, band), 0)
                    counts[(name, band)] += ev is not None
        for (name, band), k in counts.items():
            lo, hi = wilson(k, reps)
            rows.append(dict(alpha=a, beta=b, persistence=round(a + b, 6), fitted=str(tuple(round(x, 4) for x in fit)),
                             method=name, band=band, one_day_false_alarm=k / reps, ci_low=lo, ci_high=hi, paths=reps))
        print("Q10a", a + b, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e10a_null_false_alarms.csv", index=False)
    print(df.pivot_table(index=["method", "band"], columns="persistence", values="one_day_false_alarm").round(3).to_string())
    return df


def q10b(reps=300, post_days=5):
    """Breaks 200 steps into the first monitored day; monitoring continues for ``post_days`` days."""
    rng = np.random.default_rng(SEED + 1)
    rows = []
    pre = 2000
    brk = 200
    a, b = 0.05, 0.94
    omega = 1 - a - b
    train, _ = garch_paths(1, 30000, omega, a, b, rng)
    fit = F.garch11_fit(train[0])
    n_mon = DAY + post_days * DAY
    for kind in ("scale", "omega"):
        r, v = garch_paths(reps, pre + n_mon, omega, a, b, rng, break_at=pre + DAY + brk, break_kind=kind)
        res = {}
        for i in range(reps):
            # the oracle here knows the PRE-break model only, as any real forecaster would
            r_i = r[i]
            v_pre = F.garch11_filter(r_i, omega, a, b)
            vers = standardised_versions(r_i, v_pre, fit, pre)
            vers["oracle"] = F.standardise(r_i, v_pre)[pre:]
            for name, z in vers.items():
                for band in (0.0, 0.25):
                    ev = first_event(z, band)
                    res.setdefault((name, band), []).append(np.nan if ev is None else ev.index - (DAY + brk))
        for (name, band), d in res.items():
            d = np.asarray(d, float)
            hits = d[d >= 0]
            rows.append(dict(break_kind=kind, method=name, band=band, paths=reps,
                             false_alarm_before_break=float(np.mean(d < 0)),
                             detected_within_1_day=float(np.mean((d >= 0) & (d < DAY))),
                             detected_within_5_days=float(len(hits) / reps),
                             no_alarm=float(np.mean(np.isnan(d))),
                             median_delay=float(np.median(hits)) if len(hits) else np.nan))
        print("Q10b", kind, "done", flush=True)
    df = pd.DataFrame(rows).sort_values(["break_kind", "band", "method"])
    df.to_csv(OUT / "e10b_power.csv", index=False)
    print(df.round(3).to_string(index=False))
    return df


def _alarm_validation(z_alarm_stream, u, fc_var, index, band, spd=DAY, label=""):
    """Run the monitor on ``z_alarm_stream``; for each alarm compare next-day realised variance of ``u``
    with the EWMA forecast of ``u`` variance at the alarm time."""
    mon = VolatilityMonitor(threshold=A, cal_len=spd, null_band=band)
    ev = mon.run(z_alarm_stream)
    out = []
    for e in ev:
        seg = u[e.index + 1: e.index + 1 + spd]
        seg = seg[np.isfinite(seg)]
        if seg.size < spd // 2 or not np.isfinite(fc_var[e.index + 1]):
            continue
        out.append(dict(time=index[e.index], direction=e.direction, realised_over_forecast=float(np.mean(seg**2) / fc_var[e.index + 1])))
    return pd.DataFrame(out), len(ev)


def q10c():
    rows, events = [], []
    for sym in ("BTCUSDT", "ETHUSDT", "XRPUSDT", "SOLUSDT", "DOGEUSDT"):
        df = load_binance(sym, "2025-06-01", "2026-01-01", download=False)
        r = df.logret
        s = intraday_scale(r, lookback_days=20)
        u_all = (r / s).to_numpy(copy=True)
        u_all[~np.isfinite(u_all)] = np.nan
        keep = df.index >= pd.Timestamp("2025-07-01", tz="UTC")  # June is warm-up for the profile and the EWMA
        fc_all = F.ewma_variance(u_all, 1440, init=float(np.nanmean(u_all[:DAY] ** 2)))
        u, fc, idx = u_all[keep], fc_all[keep], df.index[keep]
        days = keep.sum() / DAY
        rng = np.random.default_rng(SEED + 7)
        rand = rng.integers(DAY, len(u) - DAY - 1, 2000)
        ctrl = np.array([np.nanmean(u[t + 1: t + 1 + DAY] ** 2) / fc[t + 1] for t in rand])
        variants = {"deseasonalised only (as E5)": u}
        for hl in (360, 1440):
            fch = F.ewma_variance(u_all, hl, init=float(np.nanmean(u_all[:DAY] ** 2)))[keep]
            variants[f"divided by EWMA forecast, half-life {hl} min"] = F.standardise(u, fch)
        for name, z in variants.items():
            ev, n_alarms = _alarm_validation(np.where(np.isfinite(z), z, np.nan), u, fc, idx, 0.25)
            up = ev[ev.direction == "up"] if len(ev) else ev
            dn = ev[ev.direction == "down"] if len(ev) else ev
            k_up = int((up.realised_over_forecast > 1).sum()) if len(up) else 0
            k_dn = int((dn.realised_over_forecast < 1).sum()) if len(dn) else 0
            rows.append(dict(symbol=sym, monitor_input=name, alarms=n_alarms, alarms_per_day=n_alarms / days,
                             up=len(up), down=len(dn),
                             share_up_realised_above_forecast=k_up / len(up) if len(up) else np.nan,
                             up_ci=str(tuple(round(x, 3) for x in wilson(k_up, len(up)))) if len(up) else "",
                             share_down_realised_below_forecast=k_dn / len(dn) if len(dn) else np.nan,
                             down_ci=str(tuple(round(x, 3) for x in wilson(k_dn, len(dn)))) if len(dn) else "",
                             control_share_above_forecast=float(np.mean(ctrl > 1))))
            if len(ev):
                ev = ev.assign(symbol=sym, monitor_input=name)
                events.append(ev)
        print("Q10c", sym, "done", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "e10c_real_forecast_relative.csv", index=False)
    if events:
        pd.concat(events).to_csv(OUT / "e10c_events.csv", index=False)
    with pd.option_context("display.width", 200):
        print(res.round(3).to_string(index=False))
    return res


def rho(k, a, b):
    """Stationary mean of r^2 / s2 after a k-fold break, s2 the pre-break GARCH filter (unit pre-break variance)."""
    return k * (1 - b) / (1 - a - b + a * k)


def q10b2(reps=300, post_days=5):
    """Q10b2: the Q10b monitors on paths with no break (control) and with a 4-fold omega break."""
    rng = np.random.default_rng(SEED + 1)
    pre, brk = 2000, 200
    a, b = 0.05, 0.94
    omega = 1 - a - b
    train, _ = garch_paths(1, 30000, omega, a, b, rng)
    fit = F.garch11_fit(train[0])  # identical to the Q10b fit (same seed, same first draw)
    rng = np.random.default_rng(SEED + 11)
    n_mon = DAY + post_days * DAY
    rows = []
    for kind, factor in (("none", 1.0), ("omega", 4.0)):
        r, v = garch_paths(reps, pre + n_mon, omega, a, b, rng,
                           break_at=None if kind == "none" else pre + DAY + brk, break_kind=kind, factor=factor)
        res = {}
        for i in range(reps):
            v_pre = F.garch11_filter(r[i], omega, a, b)
            vers = standardised_versions(r[i], v_pre, fit, pre)
            vers["oracle"] = F.standardise(r[i], v_pre)[pre:]
            for name, z in vers.items():
                for band in (0.0, 0.25):
                    ev = first_event(z, band)
                    res.setdefault((name, band), []).append(np.nan if ev is None else ev.index - (DAY + brk))
        for (name, band), d in res.items():
            d = np.asarray(d, float)
            hits = d[d >= 0]
            rows.append(dict(break_kind=kind, factor=factor, method=name, band=band, paths=reps,
                             alarm_before_break_time=float(np.mean(d < 0)),
                             alarm_within_1_day=float(np.mean((d >= 0) & (d < DAY))),
                             alarm_within_5_days=float(len(hits) / reps),
                             no_alarm=float(np.mean(np.isnan(d))),
                             median_delay=float(np.median(hits)) if len(hits) else np.nan))
        print("Q10b2", kind, factor, "done", flush=True)
    df = pd.DataFrame(rows).sort_values(["break_kind", "band", "method"])
    df.to_csv(OUT / "e10b2_control_and_large_break.csv", index=False)
    print(df.round(3).to_string(index=False))
    return df


def q10d(reps=200, n=20000):
    """Q10d: check rho(k) and compute the per-step Kullback-Leibler rate after a k-fold break.

    ``rho(k)`` is a ratio of means, E[r^2] / E[s2_pre]. The monitor sees the mean
    of the ratio, E[r^2 / s2_pre] = E[s2_post / s2_pre], which is larger (Jensen).
    After the transient, s2_post - s2_pre = (k - 1) omega / (1 - b) exactly, so
    as k grows the mean of the ratio tends to E[s / (s - omega / (1 - b))] for
    the pre-break GARCH variance s: it stays bounded however large the break.
    """
    rng = np.random.default_rng(SEED + 21)
    rows = []
    for a, b in SETTINGS:
        omega = 1 - a - b
        for k in (1.5, 2.0, 4.0, 10.0, 1000.0):
            r, v1 = garch_paths(reps, n, k * omega, a, b, rng, burn=5000)  # post-break steady state
            rom, mor, kl = [], [], []
            for i in range(reps):
                s0 = F.garch11_filter(r[i], omega, a, b, init=float(v1[i, 0]) / k)
                s1 = v1[i]
                s0, s1 = s0[2000:], s1[2000:]  # the pre-break filter forgets its start value
                q = s1 / s0
                rom.append(np.mean(s1) / np.mean(s0))
                mor.append(np.mean(q))
                kl.append(np.mean(0.5 * (q - 1 - np.log(q))))
            I = float(np.mean(kl))
            rows.append(dict(alpha=a, beta=b, persistence=round(a + b, 6), k=k,
                             rho_formula=rho(k, a, b), ratio_of_means_sim=float(np.mean(rom)),
                             mean_of_ratio_sim=float(np.mean(mor)), mean_of_ratio_se=float(np.std(mor) / np.sqrt(reps)),
                             ratio_of_means_cap=(1 - b) / a,
                             kl_rate_per_step=I, lai_delay_steps=float(np.log(A) / I) if I > 0 else np.inf,
                             lai_delay_days=float(np.log(A) / I / DAY) if I > 0 else np.inf))
        print("Q10d", a + b, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e10d_signal_cap.csv", index=False)
    with pd.option_context("display.width", 200):
        print(df.round(4).to_string(index=False))
    return df


def q10e(reps=400, n=4000, brk=500):
    """Q10e: the transient after a break, persistence 0.99."""
    a, b = 0.05, 0.94
    omega = 1 - a - b
    rows = []
    for kind in ("omega", "scale"):
        for k in (2.0, 4.0):
            rng = np.random.default_rng(SEED + 31)
            r, v = garch_paths(reps, n, omega, a, b, rng, break_at=brk, break_kind=kind, factor=k)
            q = np.empty_like(r)
            for i in range(reps):
                q[i] = v[i] / F.garch11_filter(r[i], omega, a, b)
            post = q[:, brk:]
            cum = np.cumsum(np.mean(0.5 * (post - 1 - np.log(post)), axis=0))
            hit = int(np.argmax(cum >= np.log(A))) if cum[-1] >= np.log(A) else np.nan
            row = dict(break_kind=kind, k=k, steps_until_mean_cum_kl_reaches_lnA=hit)
            for t in (0, 10, 30, 100, 300, 1000, 3000):
                row[f"mean_ratio_t{t}"] = float(post[:, t].mean())
            for t in (100, 300, 1000):
                row[f"cum_kl_t{t}"] = float(cum[t])
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e10e_transient.csv", index=False)
    with pd.option_context("display.width", 200):
        print(df.round(3).to_string(index=False))
    return df


def main(parts=("a", "b", "c")):
    if "a" in parts:
        q10a()
    if "b" in parts:
        q10b()
    if "b2" in parts:
        q10b2()
    if "c" in parts:
        q10c()
    if "d" in parts:
        q10d()
    if "e" in parts:
        q10e()


if __name__ == "__main__":
    import sys

    main(tuple(sys.argv[1:]) or ("a", "b", "b2", "c", "d", "e"))
