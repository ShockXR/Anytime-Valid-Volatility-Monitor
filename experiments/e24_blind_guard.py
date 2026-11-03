"""E24. A proved guard against blind segments of the confidence-sequence baseline, and a forward pre-registration.

Theory: docs/theory_blind_guard.md. Code: src/avmonitor/guard.py. Report part: docs/report_parts/E24.md.

Why: in E23 the primary monitor P1 (PeriodicVolatilityMonitor(baseline="cs", delta=0.01, A = 36,000, period=1440),
band 0.25) could not raise an increase alarm for 28% of its monitored time on BTC, and three of the sixteen 2026
releases fell in such blind segments. A segment is blind when tau^2 <= m_up = (1 + b) U (E19).

Parts (keep them separate; each writes its own results/e24_*.csv):
  theory     [SYNTHETIC, exact recursions] near-blind delay floor d_min(rho) and rate bound (Proposition 2);
             results/e24_near_blind.csv.
  selection  [SYNTHETIC] Proposition 3: under the construction, the probability that the kept window's upper bound
             fails, with and without spending; results/e24_selection.csv.
  synthetic  [SYNTHETIC; the calm-day design uses REAL BTC returns, permuted by whole days] blind share, false alarms
             against the proved bounds, detection of a later increase and the calibration cost, for the unguarded
             P1, three guarded versions and E22's rank monitor S2; results/e24_synthetic.csv (summary) and
             results/e24_synthetic_paths.csv.gz (one row per path and monitor).
  real       [REAL, post hoc, exploratory] E22's design (2024 and 2025) and E23's design (2026) rerun with the guarded
             P1 through E22.job; results/e24_real_*.csv.
  e19blind   [REAL, post hoc] E19's blind diagnostic (Part 2 configurations, cs) with the guard, no daily restart;
             results/e24_e19_blind.csv.
  power      [SYNTHETIC; planning inputs are REAL rates of 2024 to 2026] the forward pre-registration's power
             calculation; results/e24_future_power.csv.
  future     the forward confirmatory test of docs/preregistration_e24_future.md; refuses to run before the data
             reach the fixed end date and before every release in the period is verified.
  figs       results/figures/e24_*.png from the CSVs.

Usage: nice -n 10 python experiments/e24_blind_guard.py <part> [--procs=2] [--paths=1000]
"""

from __future__ import annotations

import math
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from avmonitor.confseq import eb_confseq  # noqa: E402
from avmonitor.guard import (  # noqa: E402
    GuardedPeriodicVolatilityMonitor,
    GuardedVolatilityMonitor,
    near_blind_delay_floor,
    spending_schedule,
)
from avmonitor.periodic import PeriodicVolatilityMonitor  # noqa: E402
from avmonitor.running import ConformalMonitor  # noqa: E402

OUT = ROOT / "results"
FIG = OUT / "figures"
DAY = 1440
DELTA = 0.01
A_CS = DAY / (0.05 - DELTA)  # 36,000 (E8, E19, E22, E23)
A_RANK = 28_800.0  # E22's S2
C, CAL, BAND = 3.0, 1440, 0.25
NTAU = int(round(0.25 * CAL))  # 360
SEED = 24
GUARD_TRIES = 3  # K, fixed before any run: delta / 2, delta / 4, delta / 4 (halving)

# ----------------------------------------------------------------------------
# Monitors (fixed before the synthetic run)
# ----------------------------------------------------------------------------
MONITORS = ("U", "G", "GR", "GE", "R")
MONITOR_LABEL = {
    "U": "P1 unguarded (cs, daily restart)",
    "G": "guarded P1, K = 3, eta = 0, keep",
    "GR": "guarded P1, K = 3, eta = 0, rank fallback",
    "GE": "guarded P1, K = 3, eta = 0.25, rank fallback",
    "R": "S2 awake rank monitor (E22)",
}


def make_monitor(name: str, offset: int = 0, seed=None):
    base = dict(threshold=A_CS, c=C, cal_len=CAL, null_band=BAND, baseline="cs", delta=DELTA, period=DAY, offset=offset)
    if name == "U":
        return PeriodicVolatilityMonitor(**base)
    if name == "G":
        return GuardedPeriodicVolatilityMonitor(**base, guard_tries=GUARD_TRIES)
    if name == "GR":
        return GuardedPeriodicVolatilityMonitor(**base, guard_tries=GUARD_TRIES, guard_fallback="rank", rank_seed=seed)
    if name == "GE":
        return GuardedPeriodicVolatilityMonitor(**base, guard_tries=GUARD_TRIES, guard_eta=0.25, guard_fallback="rank",
                                                rank_seed=seed)
    if name == "R":  # E14's awake-B, exactly as E22's S2 (make_monitor("awake-B", ...) of e14_running_tau.py)
        return ConformalMonitor(threshold=A_RANK, cal_len=120, family="lr", band=BAND, seed=seed, restart_len=0)
    raise ValueError(name)


def bound_first_m(name: str) -> float:
    return DAY / A_RANK if name == "R" else DELTA + DAY / A_CS


def long_run_bound(name: str, minutes: float) -> float:
    """Bound on expected alarms over ``minutes`` monitored minutes (Corollary 6.6 / Theorem 2; rank: M / A)."""
    if name == "R":
        return minutes / A_RANK
    return (DELTA + minutes / A_CS) / (1 - DELTA)


# ----------------------------------------------------------------------------
# Part theory: Proposition 2
# ----------------------------------------------------------------------------


def part_theory():
    rows = []
    for A, m in ((A_CS, DAY), (390 / 0.04, 390)):
        lo, hi = 1.0, 2.0
        for _ in range(60):
            mid = (lo + hi) / 2
            if near_blind_delay_floor(mid, A, cap=m + 2) > m:
                lo = mid
            else:
                hi = mid
        for rho in (1.0, 1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0):
            rows.append(dict(A=A, period=m, rho=rho, d_min=near_blind_delay_floor(rho, A, cap=10**6),
                             max_log_growth_per_step=math.log(1 + 0.5 * (rho - 1)), eta_star=lo - 1))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "e24_near_blind.csv", index=False)
    print(d.to_string(index=False))


# ----------------------------------------------------------------------------
# Part selection: Proposition 3
# ----------------------------------------------------------------------------


def p_fail(n_cs: int, q: float, delta_window: float, reps: int, rng, batch: int = 4000) -> tuple[float, int]:
    """P(U_T < q) for one window's upper sequence (level delta_window / 2) on i.i.d. Bernoulli(q), n_plan = n_cs."""
    hits = 0
    done = 0
    while done < reps:
        b = min(batch, reps - done)
        z = (rng.random((b, n_cs)) < q).astype(float)
        U = eb_confseq(z, 1.0, delta_window / 2, delta_window / 2, n_plan=n_cs).upper[:, -1]
        hits += int((U < q).sum())
        done += b
    return hits / reps, hits


def part_selection():
    rng = np.random.default_rng([SEED, 3])
    rows = []
    q = 1 / (1 + BAND)  # eta = 0
    configs = (("P1 window (n_cs = 1080), delta = 0.01", CAL - NTAU, DELTA, 400_000),
               ("illustration (n_cs = 150), delta = 0.3", 150, 0.3, 200_000))
    for label, n_cs, delta, reps in configs:
        cache = {}
        for K in (1, 2, 3, 4):
            for scheme, deltas in (("same delta for every window (no spending)", (delta,) * K),
                                   ("halving (Theorem 1)", spending_schedule(delta, K))):
                ps = []
                for dk in deltas:
                    key = round(dk, 12)
                    if key not in cache:
                        cache[key] = p_fail(n_cs, q, dk, reps, rng)
                    ps.append(cache[key][0])
                fail = 1 - float(np.prod([1 - p for p in ps]))
                rows.append(dict(config=label, n_cs=n_cs, q=q, delta_total=delta, K=K, scheme=scheme,
                                 deltas=";".join(f"{d:.6g}" for d in deltas),
                                 p_single_window=";".join(f"{p:.6g}" for p in ps),
                                 p_kept_window_fails=fail, upper_side_budget=sum(deltas) / 2,
                                 reps_per_level=reps,
                                 mc_se_single=";".join(f"{math.sqrt(p * (1 - p) / reps):.2g}" for p in ps)))
    # direct check of the formula through the guarded monitor (illustration, K = 3, no spending, delta 0.3 per window)
    from avmonitor.robust import MAD_NORMAL_SCALE

    tau = C * MAD_NORMAL_SCALE
    small = dict(c=C, cal_len=200, null_band=BAND, baseline="cs")
    P, K, d = 4000, 3, 0.3
    kept = 0
    for _ in range(P):
        parts = []
        for _k in range(K):
            parts += [rng.choice([-1.0, 1.0], 50), np.where(rng.random(150) < q, tau, 0.0)]
        g = GuardedVolatilityMonitor(threshold=1e6, delta=0.95, guard_tries=K, guard_deltas=(d,) * K, **small)
        g.run(np.concatenate(parts))
        kept += g.segments[0]["guard"] == "kept"
    p1, _ = p_fail(150, q, d, 200_000, rng)
    rows.append(dict(config="direct check through GuardedVolatilityMonitor (cal 200, n_cs = 150)", n_cs=150, q=q,
                     delta_total=0.95, K=K, scheme="same delta 0.3 for every window", deltas="0.3;0.3;0.3",
                     p_single_window=f"{p1:.6g}", p_kept_window_fails=kept / P, upper_side_budget=0.45,
                     reps_per_level=P, mc_se_single=f"{math.sqrt(kept / P * (1 - kept / P) / P):.2g}",
                     formula_1_minus_prod=1 - (1 - p1) ** K))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "e24_selection.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 60):
        print(d.to_string(index=False))


# ----------------------------------------------------------------------------
# Part synthetic
# ----------------------------------------------------------------------------
# Designs (fixed before the run): a stream of 11 days; volatility as below until the change at t_c, a uniform minute
# of day 9 (t_c in [8 * 1440, 9 * 1440)), after which the standard deviation doubles until the end.
DESIGNS = {
    "step2_gauss": ("step", 2.0, "gauss"),
    "step4_gauss": ("step", 4.0, "gauss"),
    "step6_gauss": ("step", 6.0, "gauss"),
    "step8_gauss": ("step", 8.0, "gauss"),
    "step8_t3": ("step", 8.0, "t3"),
    "drift8_gauss": ("drift", 8.0, "gauss"),
    "calm_btc": ("calm", 1.0, "real"),
}
DESIGN_TEXT = {
    "step": "sd 1 for the first 360 minutes (the tau part of the first window), then sd k (constant law)",
    "drift": "sd rises log-linearly from 1 to k over the first calibration window (1440 minutes), then sd k",
    "calm": "whole calm days of BTC (E19's pool, July to December 2025, deseasonalised) in random order",
}
N_DAYS, CHANGE_DAY = 11, 8
N_STEPS = N_DAYS * DAY
BATCH = 50
_POOL = None  # calm-day pool (days x 1440), set in the main process before the workers fork


def calm_day_pool() -> np.ndarray:
    import e19_real_null as E19

    keys, u, in_period = E19.load_crypto("BTC")
    _, mask = E19.calm_flags(keys, u, in_period, E19.MIN_FINITE["crypto"])
    days = []
    k, start = np.unique(keys, return_index=True)
    order = np.argsort(start)
    start = start[order]
    end = np.r_[start[1:], len(u)]
    for a, b in zip(start, end):
        if b - a == DAY and mask[a:b].all() and np.isfinite(u[a:b]).all():
            days.append(u[a:b])
    return np.array(days)


def gen_path(design: str, rng: np.random.Generator):
    kind, k, law = DESIGNS[design]
    t_c = CHANGE_DAY * DAY + int(rng.integers(0, DAY))
    t = np.arange(N_STEPS)
    if kind == "calm":
        idx = rng.choice(_POOL.shape[0], N_DAYS, replace=False)
        u = _POOL[idx].reshape(-1).copy()
    else:
        z = rng.standard_normal(N_STEPS) if law == "gauss" else rng.standard_t(3, N_STEPS)
        if kind == "step":
            s = np.where(t < NTAU, 1.0, k)
        else:
            s = np.where(t < CAL, k ** (t / CAL), k)
        u = z * s
    u[t_c:] *= 2.0
    return u, t_c


def segment_rows(mon, n: int):
    """(start, end, alarm index or -1, direction, tau2 / m_up or inf) per segment, in order."""
    alarm_of = {e.monitored_since: e for e in mon.events}
    out = []
    for seg in mon.segments:
        s = int(seg["start"])
        e = alarm_of.get(s)
        m_up = seg.get("m_up", float("nan"))
        ratio = float(seg["tau2"] / m_up) if (m_up is not None and np.isfinite(m_up) and m_up > 0) else math.inf
        out.append((s, int(e.index) if e is not None else n - 1, int(e.index) if e is not None else -1,
                    e.direction if e is not None else "", ratio, seg.get("guard", ""), seg.get("try_index", 1)))
    return out


def analyse(mon, name: str, n: int, t_c: int, n_periods: int = 5) -> dict:
    segs = segment_rows(mon, n)
    rec = dict(monitor=name)
    if not segs:
        return rec
    s1, e1, a1, d1, r1, g1, k1 = segs[0]
    rec.update(seg1_start=s1, seg1_alarm=a1, seg1_dir=d1, seg1_ratio=r1, seg1_guard=g1, seg1_try=k1,
               first_m=bool(0 <= a1 < t_c and a1 - s1 < DAY))
    b0 = s1 + (-s1) % DAY
    edges = [s1] + [b0 + i * DAY for i in range(n_periods)]
    for i in range(n_periods):
        lo, hi = edges[i], edges[i + 1]
        rec[f"period{i + 1}_valid"] = hi <= t_c
        rec[f"period{i + 1}_alarm"] = bool(0 <= a1 < t_c and lo <= a1 < hi)
    disc = getattr(mon, "discarded", [])
    w1 = [d for d in disc if d["window_start"] == 0]
    rec["ratio_window1"] = w1[0]["tau2_over_m_up"] if w1 else r1
    rec["extra_cal_seg1"] = int(sum(d["window_len"] for d in disc if d["window_end"] < s1))
    rec["extra_cal_before_tc"] = int(sum(d["window_len"] for d in disc if d["window_end"] < t_c))
    mon_min = blind_min = alarms = 0
    awake = False
    for s, e, a, d, r, g, k in segs:
        lo, hi = s, min(e, t_c - 1)
        if hi >= lo:
            mon_min += hi - lo + 1
            if r <= 1.0:
                blind_min += hi - lo + 1
        if 0 <= a < t_c:
            alarms += 1
        if s <= t_c <= e and (a < 0 or a >= t_c):
            awake = r > 1.0
    rec.update(monitored_before_tc=mon_min, blind_before_tc=blind_min, alarms_before_tc=alarms, awake_at_tc=awake,
               seg1_blind=bool(r1 <= 1.0))
    ups = [e.index for e in mon.events if e.direction == "up" and e.index >= t_c]
    rec["up_delay"] = (ups[0] - t_c + 1) if ups else -1
    rec["detected_1d"] = bool(ups and ups[0] - t_c < DAY)
    rec["detected_2d"] = bool(ups and ups[0] - t_c < 2 * DAY)
    return rec


def synth_job(args):
    design, batch = args
    rows = []
    for p in range(batch * BATCH, (batch + 1) * BATCH):
        rng = np.random.default_rng([SEED, list(DESIGNS).index(design), p])
        u, t_c = gen_path(design, rng)
        for name in MONITORS:
            mon = make_monitor(name, 0, seed=np.random.default_rng([SEED, 99, list(DESIGNS).index(design), p]))
            mon.run(u)
            rows.append(dict(design=design, path=p, t_c=t_c, **analyse(mon, name, u.size, t_c)))
    print(f"synthetic {design} batch {batch} done", flush=True)
    return rows


def wilson_se(k: float, n: int) -> float:
    p = k / n if n else float("nan")
    return math.sqrt(p * (1 - p) / n) if n else float("nan")


def summarise_synthetic(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (design, name), g in d.groupby(["design", "monitor"], sort=False):
        P = len(g)
        rec = dict(design=design, design_text=DESIGN_TEXT[DESIGNS[design][0]], k=DESIGNS[design][1],
                   law=DESIGNS[design][2], monitor=name, monitor_label=MONITOR_LABEL[name], paths=P)
        fm = g["first_m"].mean()
        rec.update(first_m_rate=fm, first_m_se=wilson_se(g["first_m"].sum(), P), first_m_bound=bound_first_m(name))
        per = []
        for i in range(1, 6):
            v = g[g[f"period{i}_valid"]]
            per.append((v[f"period{i}_alarm"].mean() if len(v) else np.nan, len(v)))
        j = int(np.nanargmax([p for p, _ in per]))
        rec.update(per_period_max_rate=per[j][0], per_period_max_which=j + 1,
                   per_period_max_se=wilson_se(per[j][0] * per[j][1], per[j][1]),
                   per_period_rates=";".join(f"{p:.4f}" for p, _ in per), per_period_n=";".join(str(n) for _, n in per),
                   per_period_bound=bound_first_m(name) if name != "R" else np.nan)
        M = g["monitored_before_tc"].to_numpy(float)
        rec.update(alarms_per_path=g["alarms_before_tc"].mean(), alarms_per_path_se=g["alarms_before_tc"].std(ddof=1) / math.sqrt(P),
                   monitored_days_per_path=M.mean() / DAY,
                   long_run_bound_per_path=float(np.mean([long_run_bound(name, m) for m in M])))
        rec.update(blind_share=g["blind_before_tc"].sum() / M.sum(),
                   blind_share_se=float((g["blind_before_tc"] / np.maximum(M, 1)).std(ddof=1) / math.sqrt(P)),
                   seg1_blind_share=g["seg1_blind"].mean(), awake_at_change=g["awake_at_tc"].mean(),
                   awake_at_change_se=wilson_se(g["awake_at_tc"].sum(), P),
                   ratio_first_window_median=float(np.median(g["ratio_window1"])),
                   ratio_seg1_median=float(np.median(g["seg1_ratio"].replace(np.inf, np.nan).dropna()))
                   if g["seg1_ratio"].replace(np.inf, np.nan).notna().any() else np.nan,
                   seg1_try_mean=g["seg1_try"].mean(), seg1_rank_fallback=(g["seg1_guard"] == "rank fallback").mean(),
                   seg1_blind_kept=(g["seg1_guard"] == "blind kept").mean(),
                   extra_cal_seg1_minutes=g["extra_cal_seg1"].mean(), seg1_start_mean=g["seg1_start"].mean(),
                   extra_cal_before_change_minutes=g["extra_cal_before_tc"].mean())
        dl = g["up_delay"].to_numpy(float)
        dl = np.where(dl < 0, np.inf, dl)
        rec.update(detected_1d=g["detected_1d"].mean(), detected_1d_se=wilson_se(g["detected_1d"].sum(), P),
                   detected_2d=g["detected_2d"].mean(), median_delay=float(np.median(dl)),
                   median_delay_detected_1d=float(np.median(dl[dl <= DAY])) if (dl <= DAY).any() else np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


def part_synthetic(procs: int, paths: int):
    global _POOL
    if "calm_btc" in DESIGNS:
        _POOL = calm_day_pool()
        print(f"calm-day pool: {_POOL.shape[0]} days", flush=True)
    jobs = [(d, b) for d in DESIGNS for b in range(int(math.ceil(paths / BATCH)))]
    t0 = time.time()
    if procs <= 1:
        res = [synth_job(j) for j in jobs]
    else:
        with Pool(procs) as pool:
            res = pool.map(synth_job, jobs, chunksize=1)
    d = pd.DataFrame([r for rows in res for r in rows])
    d.to_csv(OUT / "e24_synthetic_paths.csv.gz", index=False)
    s = summarise_synthetic(d)
    s.to_csv(OUT / "e24_synthetic.csv", index=False)
    print(f"synthetic: {len(d)} rows in {time.time() - t0:.0f} s", flush=True)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        cols = ["design", "monitor", "first_m_rate", "per_period_max_rate", "alarms_per_path", "long_run_bound_per_path",
                "blind_share", "seg1_blind_share", "awake_at_change", "detected_1d", "median_delay", "seg1_try_mean",
                "seg1_rank_fallback", "extra_cal_seg1_minutes", "ratio_first_window_median"]
        print(s[cols].round(4).to_string(index=False))


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def main():
    part = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "theory"
    opts = dict(a.lstrip("-").split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    procs = min(int(opts.get("procs", 2)), 2)
    t0 = time.time()
    if part == "theory":
        part_theory()
    elif part == "selection":
        part_selection()
    elif part == "synthetic":
        part_synthetic(procs, int(opts.get("paths", 1000)))
    else:
        raise SystemExit(__doc__)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
