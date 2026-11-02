"""E17. A multi-frequency monitor: e-detectors at 1, 5 and 15 minutes combined into one.

Question: does a weighted combination of e-detectors run on 1-, 5- and 15-minute block returns, with a proved
false alarm bound (docs/theory_multifreq.md, Theorem 1), beat each single frequency? Is it never much worse
than the best one (robustness to not knowing the right frequency)? Does it help under price noise?

Parts, each runnable alone (``python experiments/e17_multifrequency.py baselines null change real``):

  baselines [SYNTHETIC] True truncation levels and baselines of each law, frequency and block type
            (Monte Carlo, 2 million blocks each, seed 1700). -> results/e17_true_baselines.csv
  levels    The level A of every monitor for a proved 5% per day and the handicaps of Proposition 3.
            -> results/e17_levels.csv
  null      [SYNTHETIC] One-day false alarm rate (share of 1440-minute monitored days with at least one alarm)
            under Gaussian, Student t3, jump (0.5% per minute, 8 sd) and price-noise (gamma = 0.5, 1 and 2, the
            model of ALT1) returns, with the true baseline and with the plug-in baseline estimated on the
            preceding 1440 minutes. 10,000 days per cell, seeds (17, law, baseline, chunk).
            -> results/e17_null.csv
  change    [SYNTHETIC] Detection with the true baseline: persistent sd changes x1.25, x1.5, x2, x0.8, x2/3, x0.5;
            bursts x2 and x3 lasting 15 or 60 minutes; a trend (lag-one autocorrelation 0 -> 0.3 with the
            one-minute variance unchanged, so only block variances move). Change at monitored minute 247 (not a
            multiple of 5 or 15, so coarse blocks straddle it), 1440 minutes after it; Gaussian, Student t3 and
            price noise gamma = 1 and 2. 2000 paths per cell, seeds (171, law, cell, chunk).
            -> results/e17_change.csv, results/e17_change_summary.csv
  real      [REAL] Five Binance coins, 1 July to 31 December 2025, E6 validation. -> results/e17_real*.csv
  figures   -> results/figures/e17_*.png

Monitors (fixed before any run). All have a proved bound of 5% per 1440-minute day under their nulls with the
true baseline, except the row marked "no proof at this level":
  * 1 min, 5 min, 15 min alone (A = 28,800, 5,760, 1,920; each side alarms at 2A);
  * certified averages of 1, 5 and 15 minutes with weights "equal" (1/3 each), "balanced" (proportional to k:
    1/21, 5/21, 15/21) and "sqrt" (proportional to sqrt(k)), A = 28,800 x sum_k w_k / k;
  * the plain average with balanced weights checked every minute at the certified level (Theorem 3 proves only
    3 times the budget for it);
  * the weighted maximum with balanced weights (Theorem 3, same 5% bound);
  * under price noise also the gap versions (first return of each coarse block left out) and the balanced
    average of 5 and 15 minutes only.
All with c = 3, the default bet grids, band 0 in the synthetic parts.

REAL settings, FROZEN before the real part was run for the first time (committed with this script):
  * data: Binance spot BTCUSDT, ETHUSDT, XRPUSDT, SOLUSDT, DOGEUSDT 1-minute bars via ``load_binance(...,
    download=False)``, monitored from 1 July 2025 00:00 to 31 December 2025 23:59 UTC; returns divided by the
    predictable intraday profile (``intraday_scale``, 20-day lookback), exactly E6's ``prepare``;
  * monitors: 1 min, 5 min, 15 min, the three certified averages, the balanced weighted maximum; plain blocks
    (no gap); band 0.25 (E6's main band); calibration 1440 minutes, recalibration after every alarm; A as above;
  * validation exactly as E6 (``e06b_five_coins_monitor.summarise`` and ``forecast_hits``): share of up alarms
    followed by higher robust volatility over the next 1440 minutes than the calibration baseline (one-minute
    scaled MAD of the calibration window) and of down alarms followed by lower, with Wilson intervals and the
    matched-delay random-time control (seed 17 + coin index); share above (up) or below (down) an EWMA forecast,
    against 3000 random times; alarms per day;
  * descriptive only: share of each frequency in the combined statistic at alarms; share of alarms with a
    same-direction 1-minute alarm within 60 minutes.
  * the reading compares each combination with the best single frequency per metric; no other comparison and no
    tuning on the real data.

Runtime on 2 cores of a shared machine: baselines about 10 seconds, null and change about 1 minute each, real about
1 minute (data already in data/raw), figures 1 second.
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from avmonitor.multifreq import (  # noqa: E402
    block_length,
    block_returns,
    component_paths,
    first_alarm,
    handicap,
    level_for_budget,
    minute_values,
    plugin_calibration,
    proved_bound,
    resolve_weights,
)
from avmonitor.robust import MAD_NORMAL_SCALE  # noqa: E402

OUT = ROOT / "results"
FIG = OUT / "figures"
DAY = 1440
ALPHA = 0.05
C = 3.0
FREQS = (1, 5, 15)
LAWS = ("gaussian", "student_t3", "jumps", "noise_0.5", "noise_1", "noise_2")
NOISE = {"noise_0.5": 0.5, "noise_1": 1.0, "noise_2": 2.0}
N_BASE = 2_000_000
CHUNK = 1000


# ----------------------------------------------------------------------------
# Monitors
# ----------------------------------------------------------------------------


def monitor_specs(noise: bool):
    """(name, freqs, weights, combine, certified, gap, level_check) for every monitor of the synthetic parts."""
    specs = [("1 min", (1,), [1.0], "average", True, False, "certified"),
             ("5 min", (5,), [1.0], "average", True, False, "certified"),
             ("15 min", (15,), [1.0], "average", True, False, "certified")]
    for wn in ("equal", "balanced", "sqrt"):
        specs.append((f"average, {wn}", FREQS, wn, "average", True, False, "certified"))
    specs.append(("average, balanced, plain every minute (no proof at this level)", FREQS, "balanced", "average",
                  False, False, "certified"))
    specs.append(("max, balanced", FREQS, "balanced", "max", False, False, "certified"))
    specs.append(("average 5 and 15, balanced", (5, 15), "balanced", "average", True, False, "certified"))
    if noise:
        specs += [("5 min, gap", (5,), [1.0], "average", True, True, "certified"),
                  ("15 min, gap", (15,), [1.0], "average", True, True, "certified")]
        for wn in ("equal", "balanced", "sqrt"):
            specs.append((f"average, {wn}, gap", FREQS, wn, "average", True, True, "certified"))
        specs.append(("max, balanced, gap", FREQS, "balanced", "max", False, True, "certified"))
        specs.append(("average 5 and 15, balanced, gap", (5, 15), "balanced", "average", True, True, "certified"))
    return specs


def spec_level(spec):
    name, freqs, w, combine, certified, gap, level_check = spec
    return level_for_budget(ALPHA, DAY, freqs, w, combine, level_check)


def spec_bound(spec):
    """Proved one-day bound at the level used (Theorem 1 or 3)."""
    name, freqs, w, combine, certified, gap, _ = spec
    check = "certified" if certified else "stale"
    return proved_bound(DAY, spec_level(spec), freqs, w, combine, check)


# ----------------------------------------------------------------------------
# Laws
# ----------------------------------------------------------------------------


def innovations(law, rng, shape):
    if law == "student_t3":
        return rng.standard_t(3, shape) / np.sqrt(3.0)
    return rng.standard_normal(shape)


def simulate(law, rng, reps, n, sd=None, ar=None):
    """One-minute returns (reps, n). ``sd``: multiplier of the efficient part per minute (n,) or (reps, n);
    ``ar``: lag-one autocorrelation per minute (n,) of the efficient part (variance kept at sd^2)."""
    z = innovations(law, rng, (reps, n))
    if ar is not None:
        x = np.empty_like(z)
        prev = np.zeros(reps)
        for t in range(n):
            a = ar[t]
            prev = a * prev + np.sqrt(1 - a * a) * z[:, t]
            x[:, t] = prev
        z = x
    if sd is not None:
        z = z * sd
    if law == "jumps":
        z = z + (rng.random((reps, n)) < 0.005) * rng.normal(0.0, 8.0, (reps, n))
    if law in NOISE:
        u = rng.normal(0.0, NOISE[law], (reps, n + 1))
        z = z + np.diff(u, axis=1)
    return z


def block_law_draws(law, rng, nb, L):
    """nb draws of a block return of L efficient minutes (with the two noise ends for noise laws)."""
    z = innovations(law, rng, (nb, L)).sum(axis=1)
    if law == "jumps":
        z = z + ((rng.random((nb, L)) < 0.005) * rng.normal(0.0, 8.0, (nb, L))).sum(axis=1)
    if law in NOISE:
        z = z + rng.normal(0.0, NOISE[law], nb) - rng.normal(0.0, NOISE[law], nb)
    return z


def components(law):
    """(k, gap) pairs needed for a law."""
    base = [(1, False), (5, False), (15, False)]
    return base + ([(5, True), (15, True)] if law in NOISE else [])


# ----------------------------------------------------------------------------
# Part: true baselines
# ----------------------------------------------------------------------------


def part_baselines():
    rows = []
    for li, law in enumerate(LAWS):
        for k, gap in components(law):
            rng = np.random.default_rng([1700, li, k, int(gap)])
            L = block_length(k, gap)
            draws = np.concatenate([block_law_draws(law, rng, 200_000, L) for _ in range(N_BASE // 200_000)])
            scale = float(np.median(np.abs(draws)) * MAD_NORMAL_SCALE)
            tau2 = (C * scale) ** 2
            x = np.minimum(draws**2, tau2)
            rows.append(dict(law=law, k=k, gap=gap, block_minutes=L, n_blocks=draws.size, scale=scale, tau2=tau2,
                             m0=float(x.mean()), sd_x=float(x.std()), clipped_share=float(np.mean(draws**2 > tau2))))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUT / "e17_true_baselines.csv", index=False)


def load_baselines():
    df = pd.read_csv(OUT / "e17_true_baselines.csv")
    return {(r.law, int(r.k), bool(r.gap)): (float(r.tau2), float(r.m0)) for r in df.itertuples()}


# ----------------------------------------------------------------------------
# Part: levels and handicaps
# ----------------------------------------------------------------------------


def part_levels():
    rows = []
    for spec in monitor_specs(True):
        name, freqs, w, combine, certified, gap, _ = spec
        ww = resolve_weights(w, freqs)
        h = handicap(freqs, ww, m=DAY)
        rows.append(dict(monitor=name, freqs=" ".join(map(str, freqs)), weights=" ".join(f"{x:.4f}" for x in ww),
                         combine=combine, certified=certified, gap=gap, A=spec_level(spec),
                         proved_bound_per_day=spec_bound(spec),
                         **{f"handicap_{k}": float(hk) for k, hk in zip(freqs, h)}))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e17_levels.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Shared simulation core
# ----------------------------------------------------------------------------


def alarms_for_chunk(r_mon, comp_baselines, specs, n):
    """First alarm (index, direction) of every monitor on the monitored returns r_mon (reps, n).

    comp_baselines: {(k, gap): (tau2, m0)} with scalars or (reps,) arrays.
    """
    vals = {}
    for (k, gap), (tau2, m0) in comp_baselines.items():
        tau2a = np.asarray(tau2, float)
        x = np.minimum(block_returns(r_mon, k, gap) ** 2, tau2a[..., None] if tau2a.ndim else tau2a)
        p = component_paths(x, m0, m0, tau2)
        for cert in (True, False):
            vals[(k, gap, cert)] = (minute_values(p["up_full"], p["up_cert"], k, n, cert),
                                    minute_values(p["dn_full"], p["dn_cert"], k, n, cert))
    out = {}
    for spec in specs:
        name, freqs, w, combine, certified, gap, _ = spec
        ww = resolve_weights(w, freqs)
        A = spec_level(spec)
        up = dn = None
        for wk, k in zip(ww, freqs):
            vu, vd = vals[(k, gap and k > 1, certified if combine == "average" else False)]
            if up is None:
                up, dn = wk * vu, wk * vd
            elif combine == "average":
                up, dn = up + wk * vu, dn + wk * vd
            else:
                up, dn = np.maximum(up, wk * vu), np.maximum(dn, wk * vd)
        out[name] = first_alarm(up, dn, A)
    return out


def wilson(k, n, z=1.959964):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return centre - half, centre + half


# ----------------------------------------------------------------------------
# Part: null
# ----------------------------------------------------------------------------

N_NULL = 10_000


def null_job(args):
    li, law, bi, baseline = args
    true = load_baselines()
    specs = monitor_specs(law in NOISE)
    comps = components(law)
    counts = {s[0]: [0, 0] for s in specs}  # up, down first alarms within the day
    t0 = time.time()
    for ch in range(N_NULL // CHUNK):
        rng = np.random.default_rng([17, li, bi, ch])
        n_cal = DAY if baseline == "plugin" else 0
        r = simulate(law, rng, CHUNK, n_cal + DAY)
        cb = {}
        for k, gap in comps:
            if baseline == "true":
                cb[(k, gap)] = true[(law, k, gap)]
            else:
                _, tau2, m0, _, _ = plugin_calibration(r[:, :DAY], k, C, gap, overlap=True)
                cb[(k, gap)] = (tau2, m0)
        res = alarms_for_chunk(r[:, n_cal:], cb, specs, DAY)
        for name, (idx, d) in res.items():
            counts[name][0] += int(np.sum(d == 1))
            counts[name][1] += int(np.sum(d == -1))
    rows = []
    for spec in specs:
        up, dn = counts[spec[0]]
        lo, hi = wilson(up + dn, N_NULL)
        rows.append(dict(law=law, baseline=baseline, monitor=spec[0], days=N_NULL, alarm_days=up + dn,
                         false_alarm_rate=(up + dn) / N_NULL, wilson_lo=lo, wilson_hi=hi, up=up, down=dn,
                         A=spec_level(spec), proved_bound=spec_bound(spec)))
    print(f"null {law} {baseline}: {time.time() - t0:.0f} s", flush=True)
    return rows


def part_null():
    jobs = [(li, law, bi, b) for li, law in enumerate(LAWS) for bi, b in enumerate(("true", "plugin"))]
    with Pool(2) as pool:
        rows = sum(pool.map(null_job, jobs), [])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e17_null.csv", index=False)
    piv = df.pivot_table(index="monitor", columns=["law", "baseline"], values="false_alarm_rate", sort=False)
    pd.set_option("display.width", 250)
    print((100 * piv).round(2).to_string())


# ----------------------------------------------------------------------------
# Part: change
# ----------------------------------------------------------------------------

NU = 247
H = DAY
N_CHANGE = 2000
CHANGE_LAWS = ("gaussian", "student_t3", "noise_1", "noise_2")
CELLS = [("shift", 1.25, None), ("shift", 1.5, None), ("shift", 2.0, None), ("shift", 0.8, None),
         ("shift", 2 / 3, None), ("shift", 0.5, None), ("burst", 2.0, 15), ("burst", 2.0, 60), ("burst", 3.0, 15),
         ("burst", 3.0, 60), ("trend", 0.3, None)]


def cell_label(cell):
    kind, v, D = cell
    if kind == "shift":
        return f"sd x{v:.3g}"
    if kind == "burst":
        return f"burst x{v:.3g} for {D} min"
    return f"trend, lag-one autocorrelation {v}"


def change_job(args):
    li, law, ci, cell = args
    true = load_baselines()
    specs = monitor_specs(law in NOISE)
    comps = components(law)
    kind, v, D = cell
    n = NU + H
    idx_all = {s[0]: [] for s in specs}
    dir_all = {s[0]: [] for s in specs}
    for ch in range(N_CHANGE // CHUNK):
        rng = np.random.default_rng([171, li, ci, ch])
        t = np.arange(n)
        if kind == "shift":
            r = simulate(law, rng, CHUNK, n, sd=np.where(t >= NU, v, 1.0))
        elif kind == "burst":
            r = simulate(law, rng, CHUNK, n, sd=np.where((t >= NU) & (t < NU + D), v, 1.0))
        else:
            r = simulate(law, rng, CHUNK, n, ar=np.where(t >= NU, v, 0.0))
        res = alarms_for_chunk(r, {c: true[(law,) + c] for c in comps}, specs, n)
        for name, (idx, d) in res.items():
            idx_all[name].append(idx)
            dir_all[name].append(d)
    rows = []
    want = -1 if (kind == "shift" and v < 1) else 1
    for spec in specs:
        idx = np.concatenate(idx_all[spec[0]])
        d = np.concatenate(dir_all[spec[0]])
        before = (idx >= 0) & (idx < NU)
        window_end = NU + D + 15 if kind == "burst" else n
        hit = (idx >= NU) & (idx < window_end) & (d == want)
        wrong = (idx >= NU) & (idx < window_end) & (d == -want)
        delay = idx[hit] - NU + 1
        rows.append(dict(law=law, change=cell_label(cell), kind=kind, size=v, duration=D, monitor=spec[0],
                         paths=idx.size, false_alarm_before=float(before.mean()), detected=float(hit.mean()),
                         wrong_direction=float(wrong.mean()),
                         median_delay=float(np.median(delay)) if delay.size else np.nan,
                         q25_delay=float(np.quantile(delay, 0.25)) if delay.size else np.nan,
                         q75_delay=float(np.quantile(delay, 0.75)) if delay.size else np.nan,
                         mean_delay=float(delay.mean()) if delay.size else np.nan))
    print(f"change {law} {cell_label(cell)} done", flush=True)
    return rows


def part_change():
    jobs = [(li, law, ci, cell) for li, law in enumerate(CHANGE_LAWS) for ci, cell in enumerate(CELLS)]
    with Pool(2) as pool:
        rows = sum(pool.map(change_job, jobs), [])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e17_change.csv", index=False)
    # summary: each monitor against the best single frequency of the same cell
    summ = []
    singles = ("1 min", "5 min", "15 min")
    for (law, change), g in df.groupby(["law", "change"], sort=False):
        g = g.set_index("monitor")
        kind = g["kind"].iloc[0]
        key = "detected" if kind == "burst" else "median_delay"
        s = g.loc[[m for m in singles if m in g.index]]
        if kind == "burst":
            best = s["detected"].idxmax()
        else:  # best = shortest median delay among singles that detect at least 90% of paths
            ok = s[s["detected"] >= 0.9]
            best = ok["median_delay"].idxmin() if len(ok) else s["detected"].idxmax()
        for m, row in g.iterrows():
            summ.append(dict(law=law, change=change, kind=kind, monitor=m, best_single=best,
                             detected=row["detected"], best_single_detected=g.loc[best, "detected"],
                             median_delay=row["median_delay"], best_single_median_delay=g.loc[best, "median_delay"],
                             delay_ratio_to_best_single=row["median_delay"] / g.loc[best, "median_delay"],
                             metric=key))
    pd.DataFrame(summ).to_csv(OUT / "e17_change_summary.csv", index=False)
    pd.set_option("display.width", 250)
    piv = df.pivot_table(index="monitor", columns=["law", "change"], values="median_delay", sort=False)
    print(piv.round(1).to_string())


# ----------------------------------------------------------------------------
# Part: real data
# ----------------------------------------------------------------------------

REAL_BAND = 0.25
REAL_PERIOD = "Jul-Dec 2025"
REAL_MONITORS = [("1 min", (1,), [1.0], "average"), ("5 min", (5,), [1.0], "average"),
                 ("15 min", (15,), [1.0], "average"), ("average, equal", FREQS, "equal", "average"),
                 ("average, balanced", FREQS, "balanced", "average"), ("average, sqrt", FREQS, "sqrt", "average"),
                 ("max, balanced", FREQS, "balanced", "max")]


def real_job(ci_coin):
    ci, coin = ci_coin
    from e06b_five_coins_monitor import Series, forecast_hits, prepare, summarise

    from avmonitor.multifreq import MultiFrequencyMonitor

    t0 = time.time()
    idx, u, i0 = prepare(coin, REAL_PERIOD)
    S = Series(idx, u, i0)
    events, extra = {}, []
    for name, freqs, w, combine in REAL_MONITORS:
        A = level_for_budget(ALPHA, DAY, freqs, w, combine, "certified")
        mon = MultiFrequencyMonitor(threshold=A, freqs=freqs, weights=w, c=C, cal_len=DAY, null_band=REAL_BAND,
                                    combine=combine, check="certified")
        ev = mon.run(u[i0:])
        events[(name, REAL_BAND)] = [(e.index + i0, e.direction, e.monitored_since + i0, e.baseline_scale) for e in ev]
        rec = dict(source=coin, monitor=name, A=A, alarms=len(ev), segments=len(mon.segments),
                   extended=mon.cal_stats["extended"], fallback=mon.cal_stats["fallback"],
                   restarted=mon.cal_stats["restarted"])
        for k in freqs:
            sh = [e.shares[k] for e in ev if np.isfinite(e.shares[k])]
            rec[f"mean_share_{k}min"] = float(np.mean(sh)) if sh else np.nan
        extra.append(rec)
    # same-direction 1-minute alarm within 60 minutes (descriptive)
    ref = events[("1 min", REAL_BAND)]
    rt = np.array([e[0] for e in ref])
    rd = np.array([e[1] for e in ref])
    for rec in extra:
        ev = events[(rec["monitor"], REAL_BAND)]
        hits, diffs = [], []
        for t, d, *_ in ev:
            m = (rd == d) & (np.abs(rt - t) <= 60) if rt.size else np.array([], bool)
            hits.append(bool(np.any(m)))
            if np.any(m):
                j = np.argmin(np.where(m, np.abs(rt - t), np.inf))
                diffs.append(t - rt[j])
        rec["share_with_1min_alarm_within_60"] = float(np.mean(hits)) if hits else np.nan
        rec["median_minutes_after_1min_alarm"] = float(np.median(diffs)) if diffs else np.nan
    rows, ev_rows = summarise(coin, REAL_PERIOD, S, events, np.random.default_rng(17 + ci))
    hits, _ = forecast_hits(coin, REAL_PERIOD, S, events, np.random.default_rng(170 + ci))
    print(f"real {coin}: {time.time() - t0:.0f} s", flush=True)
    return rows, ev_rows, hits, extra


def part_real():
    from e06b_five_coins_monitor import COINS

    with Pool(2) as pool:
        res = pool.map(real_job, list(enumerate(COINS)))
    summ = pd.DataFrame(sum([r[0] for r in res], []))
    evs = pd.DataFrame(sum([r[1] for r in res], []))
    hits = pd.concat([r[2] for r in res], ignore_index=True)
    extra = pd.DataFrame(sum([r[3] for r in res], []))
    summ = summ.rename(columns={"detector": "monitor"})
    summ.to_csv(OUT / "e17_real_summary.csv", index=False)
    evs.rename(columns={"detector": "monitor"}).to_csv(OUT / "e17_real_events.csv", index=False)
    hits = hits.rename(columns={"rule": "monitor"})
    hits.to_csv(OUT / "e17_real_forecast_hits.csv", index=False)
    extra.to_csv(OUT / "e17_real_alarm_details.csv", index=False)
    # pooled over coins
    pooled = []
    days = summ[summ.monitor == REAL_MONITORS[0][0]]["days"].sum()
    for name, *_ in REAL_MONITORS:
        s = summ[summ.monitor == name]
        rec = dict(monitor=name, alarms=int(s.alarms.sum()), alarms_per_day=float(s.alarms.sum() / days),
                   up=int(s.up.sum()), down=int(s.down.sum()))
        for dn in ("up", "down"):
            n = s[f"{dn}_validated_n"].sum()
            k = np.round(s[f"{dn}_share_in_flagged_direction"] * s[f"{dn}_validated_n"]).sum()
            lo, hi = wilson(k, n)
            rec[f"{dn}_share_in_flagged_direction"] = k / n if n else np.nan
            rec[f"{dn}_wilson_lo"], rec[f"{dn}_wilson_hi"] = lo, hi
            w = s[f"{dn}_validated_n"]
            rec[f"{dn}_matched_control_share"] = float(np.average(s[f"{dn}_matched_delay_control_share"], weights=w))
            h = hits[(hits.monitor == name) & (hits.forecast == "EWMA")]
            kh, nh = int(h[f"{dn}_hits"].sum()), int(h[f"n_{dn}"].sum())
            lo, hi = wilson(kh, nh)
            rec[f"{dn}_ewma_hit_rate"] = kh / nh if nh else np.nan
            rec[f"{dn}_ewma_wilson_lo"], rec[f"{dn}_ewma_wilson_hi"] = lo, hi
        e = extra[extra.monitor == name]
        rec["share_with_1min_alarm_within_60"] = float(np.average(e.share_with_1min_alarm_within_60,
                                                                  weights=e.alarms.clip(lower=1)))
        for k in FREQS:
            col = f"mean_share_{k}min"
            if col in e and e[col].notna().any():
                rec[col] = float(np.average(e[col].fillna(0), weights=e.alarms.clip(lower=1)))
        pooled.append(rec)
    rnd = hits[(hits.monitor == "random times") & (hits.forecast == "EWMA")]
    pooled.append(dict(monitor="random times (EWMA reference)",
                       up_ewma_hit_rate=rnd.up_hits.sum() / rnd.n_up.sum(),
                       down_ewma_hit_rate=rnd.down_hits.sum() / rnd.n_down.sum()))
    ctrl = summ[summ.monitor.str.startswith("control")]
    pooled.append(dict(monitor="random times (next day over previous day)",
                       up_share_in_flagged_direction=float(np.average(ctrl.up_share_in_flagged_direction, weights=ctrl.alarms)),
                       down_share_in_flagged_direction=float(np.average(ctrl.down_share_in_flagged_direction,
                                                                        weights=ctrl.alarms))))
    df = pd.DataFrame(pooled)
    df.to_csv(OUT / "e17_real_pooled.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.round(3).to_string())


# ----------------------------------------------------------------------------
# Part: figures
# ----------------------------------------------------------------------------


def part_figures():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(OUT / "e17_change_summary.csv")
    shifts = df[(df.kind == "shift")]
    mons = ["1 min", "5 min", "15 min", "average, equal", "average, balanced", "average, sqrt", "max, balanced"]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), sharey=True)
    titles = {"gaussian": "Gaussian", "student_t3": "Student t3", "noise_1": "price noise, gamma = 1",
              "noise_2": "price noise, gamma = 2"}
    for ax, law in zip(axes, ("gaussian", "student_t3", "noise_1", "noise_2")):
        g = shifts[shifts.law == law]
        changes = list(dict.fromkeys(g.change))
        for m in mons:
            y = [g[(g.change == c) & (g.monitor == m)].median_delay.iloc[0] for c in changes]
            single = not ("average" in m or "max" in m)
            ax.plot(range(len(changes)), y, marker="o", label=m, lw=2.5 if single else 1.5, ls="--" if single else "-")
        ax.set_xticks(range(len(changes)))
        ax.set_xticklabels(changes, rotation=30, fontsize=8)
        ax.set_yscale("log")
        ax.set_title(titles[law] + (" (1 min not valid)" if law.startswith("noise") else ""), fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("median delay among detections (minutes, log scale)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=7, fontsize=8, frameon=False)
    fig.suptitle("E17 [SYNTHETIC]: median detection delay, persistent sd changes, true baseline, 5% per day proved")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "e17_delays.png", dpi=130)
    print("saved", FIG / "e17_delays.png")


PARTS = {"baselines": part_baselines, "levels": part_levels, "null": part_null, "change": part_change,
         "real": part_real, "figures": part_figures}

if __name__ == "__main__":
    for p in sys.argv[1:] or list(PARTS):
        t0 = time.time()
        PARTS[p]()
        print(f"part {p}: {time.time() - t0:.0f} s", flush=True)
