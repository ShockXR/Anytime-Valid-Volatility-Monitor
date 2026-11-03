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
    "frozen": ("frozen", 100.0, "ticks"),
    "drift8_gauss": ("drift", 8.0, "gauss"),
    "calm_btc_days": ("calm_days", 1.0, "real"),
    "calm_btc_perm": ("calm_perm", 1.0, "real"),
}
DESIGN_TEXT = {
    "step": "sd 1 for the first 360 minutes (the tau part of the first window), then sd k (constant law)",
    "frozen": "returns in price ticks: a frozen minute (0 or +-1 tick, equally likely) with probability 0.6 in the "
              "first 360 minutes and 0.1 afterwards, otherwise normal with sd k ticks (constant law after 360 minutes; "
              "the mechanism found in E23's blind BTC windows)",
    "drift": "sd rises log-linearly from 1 to k over the first calibration window (1440 minutes), then sd k",
    "calm_days": "whole calm days of BTC (E19's pool, July to December 2025, deseasonalised) in random order "
                 "(intraday structure kept; levels change at day boundaries, so no null holds exactly)",
    "calm_perm": "the returns of E19's BTC calm-day pool in a uniform random order (E19's primary stream: exchangeable)",
}
# The pre-change law is constant from T_NULL on: the null endpoints use each monitor's first segment that starts there
# or later (for the cs monitors this is always the first segment; the rank monitor starts after 120 returns and its
# first segment correctly alarms at the step).
T_NULL = {"step": 360, "frozen": 360, "drift": 1440, "calm_days": 0, "calm_perm": 0}
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
    if kind == "calm_days":
        idx = rng.choice(_POOL.shape[0], N_DAYS, replace=False)
        u = _POOL[idx].reshape(-1).copy()
    elif kind == "calm_perm":
        flat = _POOL.reshape(-1)
        u = flat[rng.choice(flat.size, N_STEPS, replace=False)].copy()
    elif kind == "frozen":
        p_frozen = np.where(t < NTAU, 0.6, 0.1)
        frozen = rng.random(N_STEPS) < p_frozen
        u = np.where(frozen, rng.integers(-1, 2, N_STEPS).astype(float), k * rng.standard_normal(N_STEPS))
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


def analyse(mon, name: str, n: int, t_c: int, t_null: int = 0, n_periods: int = 5) -> dict:
    """Per-path record. Null endpoints use the first segment starting at or after t_null ("segment 1")."""
    segs = segment_rows(mon, n)
    rec = dict(monitor=name)
    null_segs = [sg for sg in segs if sg[0] >= t_null and sg[0] < t_c]
    if not null_segs:  # no segment starts between the end of the pre-change law's settling and the change
        ups = [e.index for e in mon.events if e.direction == "up" and e.index >= t_c]
        rec.update(up_delay=(ups[0] - t_c + 1) if ups else -1, detected_1d=bool(ups and ups[0] - t_c < DAY),
                   detected_2d=bool(ups and ups[0] - t_c < 2 * DAY),
                   awake_at_tc=any(sg[0] <= t_c <= sg[1] and (sg[2] < 0 or sg[2] >= t_c) and sg[4] > 1.0 for sg in segs))
        return rec
    s1, e1, a1, d1, r1, g1, k1 = null_segs[0]
    rec["segments_before_null_start"] = len(segs) - len(null_segs)
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
    for s, e, a, d, r, g, k in null_segs:
        lo, hi = s, min(e, t_c - 1)
        if hi >= lo:
            mon_min += hi - lo + 1
            if r <= 1.0:
                blind_min += hi - lo + 1
        if 0 <= a < t_c:
            alarms += 1
    for s, e, a, d, r, g, k in segs:  # awake at the change: monitoring, with an increase side able to alarm
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
            rows.append(dict(design=design, path=p, t_c=t_c, **analyse(mon, name, u.size, t_c, T_NULL[DESIGNS[design][0]])))
    print(f"synthetic {design} batch {batch} done", flush=True)
    return rows


def wilson_se(k: float, n: int) -> float:
    p = k / n if n else float("nan")
    return math.sqrt(p * (1 - p) / n) if n else float("nan")


def cp_upper(k: int, n: int, level: float = 0.95) -> float:
    """One-sided Clopper-Pearson upper confidence limit for a binomial proportion (k of n)."""
    from scipy.stats import beta

    return float(beta.ppf(level, k + 1, n - k)) if k < n else 1.0


def summarise_synthetic(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (design, name), g0 in d.groupby(["design", "monitor"], sort=False):
        g = g0[g0["seg1_start"].notna()].copy()  # paths with a segment that starts after the pre-change law settles
        for c in g.columns:
            if c.endswith("_valid") or c.endswith("_alarm") or c in ("first_m", "seg1_blind", "awake_at_tc",
                                                                       "detected_1d", "detected_2d"):
                g[c] = g[c].astype(bool)
        P = len(g)
        rec = dict(design=design, design_text=DESIGN_TEXT[DESIGNS[design][0]], k=DESIGNS[design][1],
                   law=DESIGNS[design][2], monitor=name, monitor_label=MONITOR_LABEL[name], paths=len(g0),
                   paths_with_null_segment=P)
        fm = g["first_m"].mean()
        rec.update(first_m_rate=fm if P else np.nan, first_m_se=wilson_se(g["first_m"].sum(), P),
                   first_m_upper95=cp_upper(int(g["first_m"].sum()), P) if P else np.nan, first_m_bound=bound_first_m(name))
        # analyse() stores period 1 as [segment start, first day boundary), which is empty when the segment starts
        # on a boundary (every cs segment here starts at 1440 or 2880); per-period rates use the first four
        # non-empty periods of each path
        shift = (g["seg1_start"].astype(int) % DAY == 0).astype(int).to_numpy()
        per = []
        for i in range(1, 5):
            valid = np.array([bool(g[f"period{i + sh}_valid"].iloc[j]) for j, sh in enumerate(shift)], bool)
            alarm = np.array([bool(g[f"period{i + sh}_alarm"].iloc[j]) for j, sh in enumerate(shift)], bool)
            per.append((alarm[valid].mean() if valid.any() else np.nan, int(valid.sum())))
        j = int(np.nanargmax([p for p, _ in per])) if np.isfinite([p for p, _ in per]).any() else 0
        rec.update(per_period_max_rate=per[j][0], per_period_max_which=j + 1,
                   per_period_max_upper95=cp_upper(int(round(per[j][0] * per[j][1])), per[j][1]) if per[j][1] else np.nan,
                   per_period_max_se=wilson_se(per[j][0] * per[j][1], per[j][1]) if per[j][1] else np.nan,
                   per_period_rates=";".join(f"{p:.4f}" for p, _ in per), per_period_n=";".join(str(n) for _, n in per),
                   per_period_bound=bound_first_m(name) if name != "R" else np.nan)
        M = g["monitored_before_tc"].to_numpy(float)
        rec.update(alarms_per_path=g["alarms_before_tc"].mean(), alarms_per_path_se=g["alarms_before_tc"].std(ddof=1) / math.sqrt(P) if P > 1 else np.nan,
                   monitored_days_per_path=M.mean() / DAY if P else np.nan,
                   long_run_bound_per_path=float(np.mean([long_run_bound(name, m) for m in M])) if P else np.nan)
        rec.update(blind_share=g["blind_before_tc"].sum() / M.sum() if M.sum() > 0 else np.nan,
                   blind_share_se=float((g["blind_before_tc"] / np.maximum(M, 1)).std(ddof=1) / math.sqrt(P)) if P > 1 else np.nan,
                   seg1_blind_share=g["seg1_blind"].mean(),
                   ratio_first_window_median=float(np.median(g["ratio_window1"])) if P else np.nan,
                   ratio_seg1_median=float(np.median(g["seg1_ratio"].replace(np.inf, np.nan).dropna()))
                   if g["seg1_ratio"].replace(np.inf, np.nan).notna().any() else np.nan,
                   seg1_try_mean=g["seg1_try"].mean(), seg1_rank_fallback=(g["seg1_guard"] == "rank fallback").mean(),
                   seg1_blind_kept=(g["seg1_guard"] == "blind kept").mean(),
                   extra_cal_seg1_minutes=g["extra_cal_seg1"].mean(), seg1_start_mean=g["seg1_start"].mean(),
                   extra_cal_before_change_minutes=g["extra_cal_before_tc"].mean())
        dl = g0["up_delay"].to_numpy(float)
        dl = np.where(dl < 0, np.inf, dl)
        det1, det2 = g0["detected_1d"].astype(bool), g0["detected_2d"].astype(bool)
        aw = g0["awake_at_tc"].astype(bool)
        rec.update(awake_at_change=aw.mean(), awake_at_change_se=wilson_se(aw.sum(), len(g0)),
                   detected_1d=det1.mean(), detected_1d_se=wilson_se(det1.sum(), len(g0)),
                   detected_2d=det2.mean(), median_delay=float(np.median(dl)),
                   median_delay_detected_1d=float(np.median(dl[dl <= DAY])) if (dl <= DAY).any() else np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


def part_synthetic(procs: int, paths: int):
    global _POOL
    if any(v[0].startswith("calm") for v in DESIGNS.values()):
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
# Part real [REAL, post hoc, exploratory]: E22's and E23's designs with the guarded P1
# ----------------------------------------------------------------------------
# E24 arm name -> synthetic monitor configuration. Fixed after the synthetic part and before any real-data rerun
# (docs/iterations/E24.md, iteration 2): the recommended guard and the guard closest to E19's proposal.
REAL_ARMS = {"G1": "GE", "G0": "G"}
REAL_LABEL = {"G1": "G1 guarded P1 (recommended: K = 3, eta = 0.25, rank fallback)",
              "G0": "G0 guarded P1 (K = 3, eta = 0, keep)"}
REAL_SEED = 2407  # rank fallback tie-breakers
PERIOD_E22 = "2024-2025 (E22 design)"
PERIOD_E23 = "2026 (E23 design)"


def _patch_e22():
    """Teach E22.job the guarded arms (E22's code is imported, not changed)."""
    import e22_guaranteed_events as E22

    if getattr(E22, "_e24_patched", False):
        return E22
    orig = E22.make_arm

    def make_arm(arm, offset=0):
        if arm in REAL_ARMS:
            return make_monitor(REAL_ARMS[arm], offset, seed=np.random.default_rng([REAL_SEED]))
        return orig(arm, offset)

    E22.make_arm = make_arm
    E22.X_ARMS = tuple(E22.X_ARMS) + tuple(REAL_ARMS)
    for a in REAL_ARMS:
        E22.DETECTOR[a] = REAL_LABEL[a]
        E22.LABEL[a] = REAL_LABEL[a]
        E22.BOUND_TEXT[a] = "delta + m/A = 5% for every period of a segment, marginal, tolerance relative to the kept window"
    E22._e24_patched = True
    return E22


def _e23_period():
    import e23_fresh_releases as E23

    end = str(pd.read_csv(OUT / "e23_data_checks.csv")["stream_end_exclusive"].iloc[0])
    return dict(events=str(E23.EVENTS_2026), start=E23.LOAD_START, monitor_start=E23.MONITOR_START, end=end)


def _init_real(period):
    _patch_e22()
    if period is not None:
        import e23_fresh_releases as E23

        E23._set_period(*period.values())


def _real_job(args):
    import e22_guaranteed_events as E22

    return E22.job(args)


SCALE_COLS = ("tau_part_scale", "cs_part_scale", "tau_part_scale_logret", "cs_part_scale_logret", "tau_part_zero_share",
              "tick_logret", "cs_part_share_above_tau")


def guard_record_job(args):
    """Rerun one guarded arm (same loader, constructor and loop as E22.job) and return its guard bookkeeping: every
    window tried (kept or discarded), the rank-fallback segments, and the robust scale of the tau part and of the
    sequence part of every window with tau^2 / m_up <= 1.25 (the mechanism)."""
    import e12_crypto_events as E12
    from avmonitor.robust import MAD_NORMAL_SCALE

    E22 = _patch_e22()
    coin, arm = args
    df, u = E12.load_coin(E12.COINS[coin])
    index = df.index
    mon = E22.make_arm(arm, E22.utc_day_offset(index))
    E22.run_monitor(mon, u)
    rows = []

    r_raw = df["logret"].to_numpy(float)
    close = df["close"].to_numpy(float)

    def scales(a, b, nt):
        """Robust scales (1.4826 x median |.|) of the tau part and of the sequence part, in deseasonalised units (as
        the monitor sees them) and in raw log returns; the share of exact zeros in the tau part; one price tick
        (0.01 USDT) as a log return at the window's median price; the share of sequence-part returns above tau."""
        out = []
        for v in (u, r_raw):
            w = v[a:b + 1]
            x, y = w[:nt], w[nt:]
            x, y = x[np.isfinite(x)], y[np.isfinite(y)]
            out += [float(np.median(np.abs(x)) * MAD_NORMAL_SCALE), float(np.median(np.abs(y)) * MAD_NORMAL_SCALE)]
            if v is u:
                zero = float(np.mean(x == 0))
                above = float(np.mean(np.abs(y) > C * out[0]))
        tick = 0.01 / float(np.nanmedian(close[a:b + 1]))
        return out + [zero, tick, above]

    for d in mon.discarded:
        st, en = d["window_start"], d["window_end"]
        rows.append(dict(arm=arm, coin=coin, outcome="discarded", try_index=d["try_index"],
                         window_start_utc=str(index[st]), window_end_utc=str(index[en]), tau2_over_m_up=d["tau2_over_m_up"],
                         delta_window=d["delta_window"], **dict(zip(SCALE_COLS, scales(st, en, NTAU)))))
    for sg in mon.segments:
        st = int(sg["start"])
        L = int(sg["cal_len"])
        sc = dict(zip(SCALE_COLS, scales(st - L, st - 1, int(sg["n_tau"]) if sg.get("n_tau") else NTAU)))
        ratio = sg["tau2"] / sg["m_up"] if (sg.get("m_up") is not None and np.isfinite(sg["m_up"]) and sg["m_up"] > 0) \
            else sg.get("blind_tau2_over_m_up", np.nan)
        rows.append(dict(arm=arm, coin=coin, outcome=sg.get("guard", "kept"), try_index=sg.get("try_index", 1),
                         window_start_utc=str(index[st - L]), window_end_utc=str(index[st - 1]), tau2_over_m_up=ratio,
                         delta_window=sg.get("delta_window", np.nan), **sc,
                         monitoring_from_utc=str(index[st]) if st < len(index) else ""))
    stats = dict(arm=arm, coin=coin, **mon.guard_stats, extra_calibration_minutes=mon.extra_calibration_steps)
    return rows, stats


def run_real_period(period_label: str, period, procs: int):
    jobs = [(c, a) for a in REAL_ARMS for c in ("BTC", "ETH", "SOL")]
    with Pool(min(procs, 2), initializer=_init_real, initargs=(period,)) as pool:
        res = pool.map(_real_job, jobs, chunksize=1)
        rec = pool.map(guard_record_job, jobs, chunksize=1)
    E22 = _patch_e22()
    tests, windows, summ, alarms = E22.collect(res)
    for d in (tests, windows, summ, alarms):
        d.insert(0, "period", period_label)
    windows_rows = [r for rows, _ in rec for r in rows]
    stats = pd.DataFrame([st for _, st in rec])
    stats.insert(0, "period", period_label)
    wr = pd.DataFrame(windows_rows)
    wr.insert(0, "period", period_label)
    return tests, windows, summ, alarms, wr, stats


def pooled_tests(windows: pd.DataFrame) -> pd.DataFrame:
    import e23_fresh_releases as E23

    out = []
    for per, w in windows.groupby("period", sort=False):
        t = E23.pooled_rows(w.drop(columns=["period"]), outcomes=("up_post",))
        t.insert(0, "period", per)
        out.append(t)
    return pd.concat(out, ignore_index=True)


def committed_reference() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """P1 (and S2) as committed by E22 and E23, with the pooled test from their windows, for the comparison."""
    import e23_fresh_releases as E23

    tests, watched, summ = [], [], []
    for per, tag in ((PERIOD_E22, "e22"), (PERIOD_E23, "e23")):
        t = pd.read_csv(OUT / f"{tag}_tests.csv")
        t = t[t.arm.isin(["P1", "S2"]) & (t.outcome == "up_post") & (t["sample"] == "all")]
        w = pd.read_csv(OUT / f"{tag}_windows.csv")
        w = w[w.arm.isin(["P1", "S2"])]
        if tag == "e22":
            t = pd.concat([t, E23.pooled_rows(w, outcomes=("up_post",))], ignore_index=True)
        t = t.drop_duplicates(subset=["arm", "coin", "event_type"])
        t.insert(0, "period", per)
        tests.append(t)
        import e22_guaranteed_events as E22

        wt = E22.watched_table(w)
        wt.insert(0, "period", per)
        watched.append(wt)
        sm = pd.read_csv(OUT / f"{tag}_monitor_summary.csv")
        sm = sm[sm.arm.isin(["P1", "S2"])]
        sm.insert(0, "period", per)
        summ.append(sm)
    return pd.concat(tests, ignore_index=True), pd.concat(watched, ignore_index=True), pd.concat(summ, ignore_index=True)


def part_real(procs: int):
    print("REAL DATA, post hoc and exploratory: E22's and E23's designs with the guarded P1.", flush=True)
    E22 = _patch_e22()
    out = [run_real_period(PERIOD_E22, None, procs), run_real_period(PERIOD_E23, _e23_period(), procs)]
    tests = pd.concat([o[0] for o in out], ignore_index=True)
    windows = pd.concat([o[1] for o in out], ignore_index=True)
    summ = pd.concat([o[2] for o in out], ignore_index=True)
    alarms = pd.concat([o[3] for o in out], ignore_index=True)
    wrows = pd.concat([o[4] for o in out], ignore_index=True)
    stats = pd.concat([o[5] for o in out], ignore_index=True)
    pooled = pooled_tests(windows)
    tests = tests[(tests["sample"] == "all") & (tests.outcome.isin(["up_post", "up_post_1h"]))]
    tests = pd.concat([tests, pooled], ignore_index=True)
    tests["status"] = "post hoc, exploratory (raw p; the data were seen)"
    import e12_crypto_events as E12

    for k, n_col, d_col in (("period", "periods_with_alarm", "periods"),):
        ci = [E12.wilson(int(a), int(b)) for a, b in zip(summ[n_col], summ[d_col])]
        summ[f"{k}_alarm_share"] = summ[n_col] / summ[d_col]
        summ[f"{k}_alarm_share_lo"] = [x[0] for x in ci]
    summ["long_run_bound_per_monitored_day"] = [(DELTA + m / A_CS) / (1 - DELTA) / (m / DAY)
                                                 for m in summ["monitored_days"] * DAY]
    summ = summ.merge(stats, on=["period", "arm", "coin"], how="left")
    watched = []
    for per, w in windows.groupby("period", sort=False):
        wt = E22.watched_table(w.drop(columns=["period"]))
        wt.insert(0, "period", per)
        watched.append(wt)
    watched = pd.concat(watched, ignore_index=True)
    ref_t, ref_w, ref_s = committed_reference()
    tests.to_csv(OUT / "e24_real_tests.csv", index=False)
    windows.to_csv(OUT / "e24_real_windows.csv", index=False)
    summ.to_csv(OUT / "e24_real_summary.csv", index=False)
    watched.to_csv(OUT / "e24_real_watched.csv", index=False)
    alarms.to_csv(OUT / "e24_real_alarms.csv", index=False)
    wrows.to_csv(OUT / "e24_real_guard_windows.csv", index=False)
    # side-by-side contrast with the committed P1 and S2 rows
    cols = ["period", "arm", "coin", "event_type", "events", "events_with", "rate_event", "rate_control_matched",
            "p_stratified_permutation"]
    mine = tests[tests.outcome == "up_post"][cols]
    contrast = pd.concat([ref_t[cols], mine], ignore_index=True)
    wcols = ["period", "arm", "coin", "event_type", "releases", "monitoring_at_release", "blind_at_release",
             "watched_at_release", "up_within_6h", "up_within_6h_when_watched"]
    wsum = pd.concat([ref_w[wcols], watched[wcols]], ignore_index=True)
    contrast = contrast.merge(wsum, on=["period", "arm", "coin", "event_type"], how="left")
    contrast.to_csv(OUT / "e24_real_contrast.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 200):
        print(contrast.round(4).to_string(index=False))
        sc = ["period", "arm", "coin", "alarms", "monitored_share", "watched_share", "segments_up_off",
              "up_off_share_of_monitored", "longest_up_off_days", "period_alarm_share", "windows_tried",
              "windows_discarded", "fallback_rank", "fallback_keep", "extra_calibration_minutes"]
        print(summ[sc].round(4).to_string(index=False))
        print(ref_s[["period", "arm", "coin", "alarms", "monitored_share", "watched_share", "segments_up_off",
                     "up_off_share_of_monitored", "longest_up_off_days"]].round(4).to_string(index=False))
        print(wrows[wrows.tau2_over_m_up <= 1.25].round(6).to_string(index=False))


# ----------------------------------------------------------------------------
# Part e19blind [REAL, post hoc]: E19's blind diagnostic with the guard (no daily restart)
# ----------------------------------------------------------------------------


def e19_blind_job(job):
    import e06b_five_coins_monitor as E6
    import e19_real_null as E19

    source, period = job
    idx, u, i0 = E6.prepare(source, period)
    rows = []
    for band in E19.P2_BANDS:
        for name, kw in (("unguarded", None), ("guarded (K = 3, eta = 0, keep)", dict(guard_eta=0.0, guard_fallback="keep")),
                         ("guarded (K = 3, eta = 0.25, rank fallback)", dict(guard_eta=0.25, guard_fallback="rank"))):
            base = dict(threshold=E19.A_CS_CRYPTO, c=E19.C, cal_len=1440, null_band=band, baseline="cs", delta=E19.DELTA)
            if kw is None:
                from avmonitor.monitor import VolatilityMonitor

                mon = VolatilityMonitor(**base)
            else:
                mon = GuardedVolatilityMonitor(**base, guard_tries=GUARD_TRIES, rank_seed=np.random.default_rng([REAL_SEED]),
                                               **kw)
            ev = mon.run(u[i0:])
            n = len(u) - i0
            calib = n - sum(((e.index if e is not None else n - 1) - sg["start"] + 1)
                            for sg, e in ((sg, {x.monitored_since: x for x in mon.events}.get(sg["start"])) for sg in mon.segments))
            rows.append(dict(source=source, period=period, baseline="cs", band=band, monitor=name, alarms=len(ev),
                             **E19.side_activity(mon, n), calibrating_minutes=int(calib),
                             **({f"guard_{k}": v for k, v in mon.guard_stats.items()} if kw else {}),
                             extra_calibration_minutes=getattr(mon, "extra_calibration_steps", 0)))
    print(f"e19blind {source} {period} done", flush=True)
    return rows


def part_e19blind(procs: int):
    import e06b_five_coins_monitor as E6

    jobs = [(c, p) for p in E6.PERIODS for c in E6.COINS] + [("Bitstamp BTC", "Jul-Dec 2025")]
    with Pool(min(procs, 2)) as pool:
        res = pool.map(e19_blind_job, jobs, chunksize=1)
    d = pd.DataFrame([r for rows in res for r in rows])
    d.to_csv(OUT / "e24_e19_blind.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 300):
        print(d.drop(columns=["baseline"]).round(3).to_string(index=False))


# ----------------------------------------------------------------------------
# Forward pre-registration (docs/preregistration_e24_future.md): power and the future test
# ----------------------------------------------------------------------------
FUTURE_EVENTS = ROOT / "data/events/us_macro_events_e24_future.csv"
FUTURE_MONITOR_START = "2026-10-25"  # E23's stream ended at 2026-10-25 00:00 UTC (exclusive): the first unseen minute
FUTURE_LOAD_START = "2026-10-05"  # 20 days of intraday-profile warm-up (known data), as E12 and E23
FUTURE_END = None  # fixed by the power part's decision rule; written here and in the pre-registration, section 5
FUTURE_PRIMARY_ARM = None  # fixed by the power part's decision rule (pre-registration, section 6)
FUTURE_PRIMARY_DESIGN = None
FUTURE_ARMS = ("G1", "S2")  # the two candidate primary monitors
FUTURE_DESIGNS = ("fomc", "fomc+cpi")  # the two candidate primary designs (E23's)
FUTURE_COIN = "BTC"
POWER_TARGET = 0.80
POWER_REPS_GRID, POWER_REPS_CONFIRM, POWER_SEED_E24 = 5_000, 20_000, 24_000
NULL_SIZE_MAX = 0.055  # a candidate is eligible only if its simulated size under the planning null is at most this
# Planning only (never used by the future test): which 2027 meeting has no statement is not announced, so the plan
# drops the earliest one in the period (the fewest statements for every end date); 2028 rows are placeholders (the
# 2027 date plus 364 days, so the same weekday; FOMC: two 2028 meetings are expected to be strategy-only, so the first
# two placeholders are dropped). Placeholders are planning inputs, not event dates, and are never written to the file.
PLANNING_NO_STATEMENT = ("2027-01-27",)
PLANNING_FOMC_2028_DROPPED = 2
GRID_ENDS = [str(d.date()) for d in pd.date_range("2027-01-01", "2029-01-01", freq="MS")]


def _cell_rates(t: pd.DataFrame, arm: str, etype: str) -> tuple[int, int, float]:
    r = t[(t.arm == arm) & (t.coin == FUTURE_COIN) & (t.event_type == etype) & (t.outcome == "up_post")
          & (t["sample"] == "all")]
    if len(r) != 1:
        raise ValueError(f"expected one row for {arm} {etype}, got {len(r)}")
    r = r.iloc[0]
    return int(r.events), int(r.events_with), float(r.rate_control_matched)


def planning_inputs() -> pd.DataFrame:
    """[REAL rates used as planning inputs] per candidate monitor and release type on BTC: 2024 to 2025 (E22's design),
    2026 (E23's design) and the two together (releases summed; the matched control rate averaged over releases).
    S2: E22's and E23's committed tests. G1: this experiment's post hoc rerun (part real)."""
    e22, e23 = pd.read_csv(OUT / "e22_tests.csv"), pd.read_csv(OUT / "e23_tests.csv")
    e24 = pd.read_csv(OUT / "e24_real_tests.csv")
    rows = []
    for arm in FUTURE_ARMS:
        for etype in ("fomc", "cpi"):
            if arm == "S2":
                a, b = _cell_rates(e22, arm, etype), _cell_rates(e23, arm, etype)
                src = "results/e22_tests.csv; results/e23_tests.csv"
            else:
                a = _cell_rates(e24[e24.period == PERIOD_E22], arm, etype)
                b = _cell_rates(e24[e24.period == PERIOD_E23], arm, etype)
                src = "results/e24_real_tests.csv (post hoc rerun)"
            n = a[0] + b[0]
            rows.append(dict(arm=arm, event_type=etype, releases_2024_2025=a[0], with_up_2024_2025=a[1],
                             control_rate_2024_2025=a[2], releases_2026=b[0], with_up_2026=b[1], control_rate_2026=b[2],
                             p1_pooled=(a[1] + b[1]) / n, p0_pooled=(a[0] * a[2] + b[0] * b[2]) / n,
                             p1_2026=b[1] / b[0], p0_2026=b[2], source=src))
    return pd.DataFrame(rows)


def planning_scenarios(inp: pd.DataFrame, arm: str) -> dict:
    g = inp[inp.arm == arm].set_index("event_type")
    pooled = {t: (float(g.loc[t, "p1_pooled"]), float(g.loc[t, "p0_pooled"])) for t in ("fomc", "cpi")}
    y26 = {t: (float(g.loc[t, "p1_2026"]), float(g.loc[t, "p0_2026"])) for t in ("fomc", "cpi")}
    half = {t: (p0 + (p1 - p0) / 2, p0) for t, (p1, p0) in pooled.items()}
    null = {t: (p0, p0) for t, (_, p0) in pooled.items()}
    return {"pooled 2024 to 2026 (planning)": pooled, "2026 only (E23 design)": y26,
            "half of the pooled excess": half, "null (release rate = control rate)": null}


def planning_events(end: str) -> pd.DataFrame:
    """The release calendar of the plan: the event file's released and scheduled rows minus the planning
    no-statement meeting, plus 2028 placeholders; releases whose 6-hour window ends by ``end``."""
    import e12_crypto_events as E12

    ev = pd.read_csv(FUTURE_EVENTS, dtype=str, keep_default_na=False)
    ev = ev[ev.status.isin(["released", "rescheduled", "scheduled"]) & ~ev.date.isin(PLANNING_NO_STATEMENT)]
    ev = ev[["date", "event_type", "release_time_et"]].assign(kind="event file")
    last = pd.Timestamp(ev.date.max())
    if pd.Timestamp(end) > last:
        ph = []
        for etype in ("fomc", "cpi"):
            base = ev[(ev.event_type == etype) & ev.date.str.startswith("2027")].sort_values("date")
            if etype == "fomc":
                base = pd.concat([pd.DataFrame(dict(date=list(PLANNING_NO_STATEMENT), event_type="fomc",
                                                    release_time_et="14:00")), base]).sort_values("date")
                base = base.iloc[PLANNING_FOMC_2028_DROPPED:]
            for _, r in base.iterrows():
                ph.append(dict(date=str((pd.Timestamp(r.date) + pd.Timedelta(days=364)).date()), event_type=etype,
                               release_time_et=r.release_time_et, kind="planning placeholder (2027 date + 364 days)"))
        ev = pd.concat([ev, pd.DataFrame(ph)], ignore_index=True)
    ev["date"] = pd.to_datetime(ev["date"])
    et = pd.to_datetime(ev["date"].dt.strftime("%Y-%m-%d") + " " + ev["release_time_et"]).dt.tz_localize(E12.TZ)
    ev["release_utc"] = et.dt.tz_convert("UTC")
    lo, hi = pd.Timestamp(FUTURE_MONITOR_START, tz="UTC"), pd.Timestamp(end, tz="UTC")
    keep = (ev.release_utc >= lo) & (ev.release_utc + pd.Timedelta(minutes=E12.POST) <= hi)
    return ev[keep].sort_values("release_utc").reset_index(drop=True)


def planning_calendar(end: str) -> pd.DataFrame:
    import e23_fresh_releases as E23

    events = planning_events(end)
    with E23.e12_period(FUTURE_EVENTS, FUTURE_LOAD_START, FUTURE_MONITOR_START, end):
        index = pd.date_range(pd.Timestamp(FUTURE_MONITOR_START, tz="UTC"), pd.Timestamp(end, tz="UTC"), freq="1min",
                              inclusive="left")
        cal = E23.calendar_strata(events, index)
    cal["kind"] = events["kind"].to_numpy()
    return cal


def _power_job(args):
    import e23_fresh_releases as E23

    arm, scen, rates, end, reps = args
    cal = planning_calendar(end)
    res = E23.simulate_power(cal, rates, reps=reps, seed=POWER_SEED_E24, alpha=0.05, designs=FUTURE_DESIGNS)
    return [dict(arm=arm, scenario=scen, design=d, end=end, reps=reps, fomc=int((cal.event_type == "fomc").sum()),
                 cpi=int((cal.event_type == "cpi").sum()), placeholders=int((cal.kind != "event file").sum()), **res[d])
            for d in FUTURE_DESIGNS]


def part_power(procs: int):
    """[SYNTHETIC; planning inputs are REAL rates of 2024 to 2026] Decision rule of the pre-registration, section 5,
    fixed before this part was first run: for each candidate (monitor G1 or S2, design FOMC alone or pooled), E* is the
    first grid end (first day of a month, 00:00 UTC) at which the power under the candidate's pooled 2024 to 2026 rates
    is at least 0.80 with 5,000 replications, confirmed with 20,000 (if the confirmation falls below 0.80 the next month
    is tried). A candidate is eligible if its simulated size under the null scenario at E* is at most 0.055. The
    primary test is the eligible candidate with the earliest E*; ties: higher confirmed power, then the pooled design,
    then G1. The end date is that E*. If no candidate reaches 0.80 by 2029-01-01, the one with the highest power there
    is chosen and the end date is 2029-01-01."""
    inp = planning_inputs()
    scen = {arm: planning_scenarios(inp, arm) for arm in FUTURE_ARMS}
    primary = "pooled 2024 to 2026 (planning)"
    jobs = [(arm, primary, scen[arm][primary], end, POWER_REPS_GRID) for arm in FUTURE_ARMS for end in GRID_ENDS]
    with Pool(min(procs, 2)) as pool:
        grid = pd.DataFrame([r for rows in pool.map(_power_job, jobs, chunksize=1) for r in rows])
    cand = []
    for arm in FUTURE_ARMS:
        for d in FUTURE_DESIGNS:
            g = grid[(grid.arm == arm) & (grid.design == d)].sort_values("end")
            hits = [e for e in g.end if float(g.loc[g.end == e, "power"].iloc[0]) >= POWER_TARGET]
            e_star, conf = None, None
            for e in [x for x in GRID_ENDS if hits and x >= hits[0]]:
                conf = _power_job((arm, primary, scen[arm][primary], e, POWER_REPS_CONFIRM))
                conf = [r for r in conf if r["design"] == d][0]
                if conf["power"] >= POWER_TARGET:
                    e_star = e
                    break
            e_eval = e_star or GRID_ENDS[-1]
            if conf is None or conf["end"] != e_eval:
                conf = [r for r in _power_job((arm, primary, scen[arm][primary], e_eval, POWER_REPS_CONFIRM))
                        if r["design"] == d][0]
            nul = [r for r in _power_job((arm, "null (release rate = control rate)",
                                          scen[arm]["null (release rate = control rate)"], e_eval, POWER_REPS_CONFIRM))
                   if r["design"] == d][0]
            cand.append(dict(arm=arm, design=d, end_star=e_star, end_evaluated=e_eval, power=conf["power"],
                             power_se=conf["mc_se"], fomc=conf["fomc"], cpi=conf["cpi"], placeholders=conf["placeholders"],
                             null_size=nul["power"], null_size_se=nul["mc_se"],
                             eligible=nul["power"] <= NULL_SIZE_MAX))
    cand = pd.DataFrame(cand)
    el = cand[cand.eligible & cand.end_star.notna()].copy()
    if len(el):
        el["_d"] = (el.design != "fomc+cpi").astype(int)
        el["_a"] = (el.arm != "G1").astype(int)
        choice = el.sort_values(["end_star", "power", "_d", "_a"], ascending=[True, False, True, True]).iloc[0]
        end = choice.end_star
    else:
        ok = cand[cand.eligible] if cand.eligible.any() else cand
        choice = ok.sort_values("power", ascending=False).iloc[0]
        end = GRID_ENDS[-1]
    cand["chosen"] = (cand.arm == choice.arm) & (cand.design == choice.design)
    # every scenario at the chosen end date, for every candidate (the pre-registration's table)
    jobs = [(arm, s, scen[arm][s], end, POWER_REPS_CONFIRM) for arm in FUTURE_ARMS for s in scen[arm]]
    with Pool(min(procs, 2)) as pool:
        table = pd.DataFrame([r for rows in pool.map(_power_job, jobs, chunksize=1) for r in rows])
    for arm in FUTURE_ARMS:
        for s, r in scen[arm].items():
            sel = (table.arm == arm) & (table.scenario == s)
            table.loc[sel, "p1_fomc"], table.loc[sel, "p0_fomc"] = r["fomc"]
            table.loc[sel, "p1_cpi"], table.loc[sel, "p0_cpi"] = r["cpi"]
    cal = planning_calendar(end)
    cal_rows = cal.assign(group=cal["group"].map(
        lambda g: f"{['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'][g[0]]} {g[1]} ET"))
    inp.to_csv(OUT / "e24_future_planning_inputs.csv", index=False)
    grid.to_csv(OUT / "e24_future_power_grid.csv", index=False)
    cand.to_csv(OUT / "e24_future_candidates.csv", index=False)
    table.assign(chosen_end=end).to_csv(OUT / "e24_future_power.csv", index=False)
    cal_rows.to_csv(OUT / "e24_future_power_calendar.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 300):
        print(inp.round(4).to_string(index=False))
        print(grid.pivot_table(index="end", columns=["arm", "design"], values="power").round(3).to_string())
        print(cand.round(4).to_string(index=False))
        print(table.round(4).to_string(index=False))
        print(cal_rows.to_string(index=False))
    print(f"CHOICE: arm {choice.arm}, design {choice.design}, end {end}")


def _month_complete(path: Path, month: pd.Period) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    last = pd.to_datetime(pd.read_csv(path, usecols=["timestamp"])["timestamp"].max(), unit="ms", utc=True)
    return last >= (month + 1).start_time.tz_localize("UTC") - pd.Timedelta(minutes=1)


def future_refusals(end: str | None, data_end: str | None = None) -> list[str]:
    """Reasons the future test must not run (empty when it may). Checked before any market data are read."""
    import e12_crypto_events as E12

    why = []
    if end is None:
        return ["FUTURE_END is not set (the power part has not fixed it)"]
    if pd.Timestamp.now(tz="UTC") < pd.Timestamp(end, tz="UTC"):
        why.append(f"today is before the end date {end}")
    ev = pd.read_csv(FUTURE_EVENTS, dtype=str, keep_default_na=False)
    et = pd.to_datetime(ev["date"] + " " + ev["release_time_et"].replace("", "00:00")).dt.tz_localize(E12.TZ)
    utc = et.dt.tz_convert("UTC")
    inside = (utc >= pd.Timestamp(FUTURE_MONITOR_START, tz="UTC")) & (utc < pd.Timestamp(end, tz="UTC"))
    allowed = {"released", "rescheduled", "cancelled", "no_statement"}
    for _, r in ev[inside & ~ev.status.isin(allowed)].iterrows():
        why.append(f"{r.date} {r.event_type}: status {r.status!r}; verify the release and update the event file")
    months = pd.period_range(pd.Period(FUTURE_MONITOR_START, "M"), pd.Period(pd.Timestamp(end) - pd.Timedelta(minutes=1), "M"))
    cpi_months = set(pd.to_datetime(ev.loc[ev.event_type == "cpi", "date"]).dt.to_period("M"))
    for m in months:
        if m > pd.Period(FUTURE_MONITOR_START, "M") and m not in cpi_months:
            why.append(f"no CPI row for {m} (add the release, or a cancelled row, with its verification)")
    if data_end is not None and pd.Timestamp(data_end) < pd.Timestamp(end):
        why.append(f"the data end at {data_end} UTC, before the end date {end}")
    return why


def part_future(procs: int):
    """The confirmatory test of docs/preregistration_e24_future.md. Refuses to run before the end date, before every
    release in the period is verified, and when the data end before the end date. Nothing here may be changed after
    the pre-registration's commit except FUTURE_END's value as set by the power part in that commit."""
    import urllib.error

    import e12_crypto_events as E12
    import e23_fresh_releases as E23
    from avmonitor import data as Dm

    why = future_refusals(FUTURE_END)
    if why:
        raise SystemExit("future part refused:\n  " + "\n  ".join(why))
    months = pd.period_range(pd.Period(FUTURE_MONITOR_START, "M"), pd.Period(pd.Timestamp(FUTURE_END) - pd.Timedelta(minutes=1), "M"))
    paths = {}
    for coin, sym in E12.COINS.items():
        paths[coin] = []
        for m in months:
            dest = Dm.binance_month_path(sym, m)
            try:  # a month cached before it was complete (October 2026, cut by E23 on 25 October) is downloaded again
                paths[coin].append(Dm.fetch(Dm.binance_month_url(sym, m), dest, overwrite=not _month_complete(dest, m)))
            except urllib.error.HTTPError as err:
                raise SystemExit(f"future part refused: {sym} {m} is not in the mirror ({err.code})")
    E23.LOAD_START, E23.MONITOR_START = FUTURE_LOAD_START, FUTURE_MONITOR_START
    d_end = E23.data_end(paths)
    why = future_refusals(FUTURE_END, d_end)
    if why:
        raise SystemExit("future part refused:\n  " + "\n  ".join(why))
    files, checks, gaps = E23.data_checks(paths, FUTURE_END)
    files.to_csv(OUT / "e24_future_data_files.csv", index=False)
    checks.to_csv(OUT / "e24_future_data_checks.csv", index=False)
    gaps.to_csv(OUT / "e24_future_data_gaps.csv", index=False)
    period = dict(events=str(FUTURE_EVENTS), start=FUTURE_LOAD_START, monitor_start=FUTURE_MONITOR_START, end=FUTURE_END)
    tests = future_pipeline(period, procs, "e24_future")
    print(tests[tests.family != ""].round(5).to_string(index=False))


def future_pipeline(period: dict, procs: int, prefix: str | None) -> pd.DataFrame:
    """The arms G1, S2 and P1 on one period through E22.job, the pooled rows, the families of section 6 and the
    diagnostics; writes results/<prefix>_*.csv unless prefix is None. Returns the tests."""
    import e12_crypto_events as E12
    import e23_fresh_releases as E23

    with E23.e12_period(**period):
        events, _, excluded = E12.read_events()
    arms = ("G1", "S2", "P1")
    jobs = [(c, a) for a in arms for c in E12.COINS]
    with Pool(min(procs, 2), initializer=_init_real, initargs=(period,)) as pool:
        res = pool.map(_real_job, jobs, chunksize=1)
        rec = pool.map(guard_record_job, [(c, "G1") for c in E12.COINS], chunksize=1)
    E22p = _patch_e22()
    tests, windows, summ, alarms = E22p.collect(res)
    tests = pd.concat([tests, E23.pooled_rows(windows)], ignore_index=True)
    tests["family"], tests["family_role"], tests["p_holm_family"] = "", "", np.nan
    base = (tests.outcome == "up_post") & (tests["sample"] == "all")
    for name, role, arm, coins, cells in future_families():
        sel = base & (tests.arm == arm) & tests.coin.isin(coins) & tests.event_type.isin(cells)
        tests.loc[sel, "family"], tests.loc[sel, "family_role"] = name, role
        tests.loc[sel, "p_holm_family"] = E12.holm(tests.loc[sel, "p_stratified_permutation"].to_numpy())
    tests["primary_test"] = tests.family == "F24a"
    tests["primary_passes"] = tests.primary_test & (tests.p_holm_family <= 0.05)
    g = summ.arm.isin(["G1", "P1"])
    summ.loc[g, "long_run_bound_per_monitored_day"] = [(DELTA + m / A_CS) / (1 - DELTA) / (m / DAY)
                                                       for m in summ.loc[g, "monitored_days"] * DAY]
    if prefix is not None:
        events.to_csv(OUT / f"{prefix}_events_used.csv", index=False)
        excluded.to_csv(OUT / f"{prefix}_events_excluded.csv", index=False)
        tests.to_csv(OUT / f"{prefix}_tests.csv", index=False)
        windows.to_csv(OUT / f"{prefix}_windows.csv", index=False)
        summ.to_csv(OUT / f"{prefix}_monitor_summary.csv", index=False)
        E22p.watched_table(windows).to_csv(OUT / f"{prefix}_watched.csv", index=False)
        alarms.to_csv(OUT / f"{prefix}_alarms.csv", index=False)
        pd.DataFrame([r for rows, _ in rec for r in rows]).to_csv(OUT / f"{prefix}_guard_windows.csv", index=False)
        pd.DataFrame([st for _, st in rec]).to_csv(OUT / f"{prefix}_guard_stats.csv", index=False)
    return tests


def part_future_gate(procs: int):
    """Gate of the pre-registration, section 8 [REAL, known data]: the future pipeline, run on E23's period and
    events (no download), must reproduce E23's committed P1 and S2 test rows exactly (every outcome and sample,
    the 14 columns of E23's gate). Writes results/e24_future_gate.csv."""
    period = _e23_period()
    tests = future_pipeline(period, procs, None)
    ref = pd.read_csv(OUT / "e23_tests.csv")
    cols = ["events", "events_with", "rate_event", "rate_event_lo", "rate_event_hi", "rate_control_matched",
            "expected_from_matched_controls", "expected_under_null", "p_stratified_permutation",
            "p_stratified_permutation_fewer", "p_stratified_fisher", "weekday_time_groups", "rate_ratio_matched",
            "control_windows"]
    cols = [c for c in cols if c in ref.columns and c in tests.columns]
    rows = []
    for arm in ("P1", "S2"):
        r = ref[ref.arm == arm]
        for _, x in r.iterrows():
            m = tests[(tests.arm == arm) & (tests.coin == x.coin) & (tests.event_type == x.event_type)
                      & (tests.outcome == x.outcome) & (tests["sample"] == x["sample"])]
            if len(m) != 1:
                rows.append(dict(arm=arm, coin=x.coin, event_type=x.event_type, outcome=x.outcome, sample=x["sample"],
                                 compared=0, mismatches=1, columns_differing="row missing"))
                continue
            m = m.iloc[0]
            bad = [c for c in cols if not (pd.isna(x[c]) and pd.isna(m[c])) and not np.isclose(float(x[c]), float(m[c]),
                                                                                                rtol=1e-12, atol=0)]
            rows.append(dict(arm=arm, coin=x.coin, event_type=x.event_type, outcome=x.outcome, sample=x["sample"],
                             compared=len(cols), mismatches=len(bad), columns_differing=";".join(bad)))
    g = pd.DataFrame(rows)
    g.to_csv(OUT / "e24_future_gate.csv", index=False)
    print(f"future gate: {len(g)} rows of E23's P1 and S2 tests compared, {int(g.mismatches.sum())} mismatches")
    if g.mismatches.sum():
        print(g[g.mismatches > 0].to_string(index=False))


def future_families():
    """(family, role, arm, coins, cells) of the pre-registration, section 6, from the power part's choice."""
    if FUTURE_PRIMARY_ARM is None or FUTURE_PRIMARY_DESIGN is None:
        raise SystemExit("the primary arm and design are not set (the power part has not fixed them)")
    other = [a for a in FUTURE_ARMS if a != FUTURE_PRIMARY_ARM][0]
    cells = tuple(c for c in ("fomc", "cpi", "fomc+cpi") if c != FUTURE_PRIMARY_DESIGN)
    return (("F24a", "primary", FUTURE_PRIMARY_ARM, ("BTC",), (FUTURE_PRIMARY_DESIGN,)),
            ("F24b", "secondary (primary arm and design on ETH and SOL)", FUTURE_PRIMARY_ARM, ("ETH", "SOL"),
             (FUTURE_PRIMARY_DESIGN,)),
            ("F24c", "secondary (the other candidate monitor, primary design)", other, ("BTC", "ETH", "SOL"),
             (FUTURE_PRIMARY_DESIGN,)),
            ("F24d", "secondary (primary arm, the other cells)", FUTURE_PRIMARY_ARM, ("BTC", "ETH", "SOL"), cells))


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a5cc2"]  # categorical slots 1 to 5 (as in TH5 and TH6)
MARKERS = ["o", "s", "D", "^", "v"]
INK, MUTED, GRID_C, SURFACE = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK, "font.size": 9,
                         "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})
    return plt


def fig_near_blind(plt):
    d = pd.read_csv(OUT / "e24_near_blind.csv")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for i, ((A, m), g) in enumerate(d.groupby(["A", "period"], sort=False)):
        g = g[g.rho > 1.0]
        ax.plot(g.rho - 1, g.d_min, color=SERIES[i], marker=MARKERS[i], ms=4, lw=1.6,
                label=f"A = {A:,.0f}, period m = {m} (eta* = {g.eta_star.iloc[0]:.4f})")
        ax.axhline(m, color=SERIES[i], lw=0.8, ls=":")
        ax.annotate(f"m = {m}", xy=(g.rho.max() - 1, m), xytext=(2, 3), textcoords="offset points", fontsize=7.5,
                    color=SERIES[i], ha="right")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(color=GRID_C, lw=0.6)
    ax.set_xlabel("rho - 1, where rho = tau^2 / m_up of the window")
    ax.set_ylabel("fewest steps to an increase alarm, d_min(rho)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_title("E24: a near-blind window cannot alarm within a period (Proposition 2) [exact recursion]",
                 loc="left", fontsize=9.5, color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "e24_near_blind.png", dpi=130)
    plt.close(fig)


def fig_synthetic(plt):
    s = pd.read_csv(OUT / "e24_synthetic.csv")
    designs = list(dict.fromkeys(s.design))[::-1]
    panels = (("blind_share", None, "share of monitored time before the change\nin blind segments (tau^2 <= m_up)"),
              ("first_m_rate", "first_m_se", "alarm within the first day of the first\nnull segment (bound 0.05)"),
              ("detected_1d", "detected_1d_se", "up alarm within one day after the\nvolatility doubles"))
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.6), sharey=True)
    ys = np.arange(len(designs), dtype=float)
    off = np.linspace(-0.28, 0.28, len(MONITORS))
    for ax, (col, se, title) in zip(axes, panels):
        for j, name in enumerate(MONITORS):
            g = s[s.monitor == name].set_index("design").reindex(designs)
            x = g[col].to_numpy(float)
            ax.scatter(x, ys + off[j], color=SERIES[j], marker=MARKERS[j], s=26, zorder=3,
                       label=MONITOR_LABEL[name] if ax is axes[0] else None)
            if se is not None:
                e = 2 * g[se].to_numpy(float)
                ax.hlines(ys + off[j], x - e, x + e, color=SERIES[j], lw=1.0, zorder=2)
        if col == "first_m_rate":
            ax.axvline(0.05, color=MUTED, lw=1, ls="--")
            ax.annotate("delta + m/A = 0.05", xy=(0.05, len(designs) - 0.4), fontsize=7.5, color=MUTED, ha="left",
                        xytext=(3, 0), textcoords="offset points")
        ax.set_title(title, loc="left", fontsize=9, color=INK)
        ax.grid(axis="x", color=GRID_C, lw=0.6)
        ax.set_xlim(-0.02, 1.02)
        for y in ys[:-1]:
            ax.axhline(y + 0.5, color=GRID_C, lw=0.5)
    axes[0].set_yticks(ys, designs)
    fig.legend(loc="lower left", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.08, 0.0))
    fig.suptitle("E24 [SYNTHETIC]: blind share, first-day alarms and detection, per design (dots with +-2 Monte Carlo "
                 "s.e.); calm_btc_* designs use real BTC calm days", x=0.01, ha="left", fontsize=10, color=INK)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.86, bottom=0.17, wspace=0.12)
    fig.savefig(FIG / "e24_synthetic.png", dpi=130)
    plt.close(fig)


def fig_real(plt):
    c = pd.read_csv(OUT / "e24_real_contrast.csv")
    c = c[c.coin == "BTC"]
    arms = ["P1", "G0", "G1", "S2"]
    periods = [PERIOD_E22, PERIOD_E23]
    cells = ("fomc+cpi", "fomc", "cpi")
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 6.2), sharex=True, sharey=True)
    ys = np.arange(len(arms))[::-1].astype(float)
    for r, per in enumerate(periods):
        for k, cell in enumerate(cells):
            ax = axes[r, k]
            for y, arm in zip(ys, arms):
                row = c[(c.period == per) & (c.arm == arm) & (c.event_type == cell)]
                if not len(row):
                    continue
                row = row.iloc[0]
                ev, ct = float(row.rate_event), float(row.rate_control_matched)
                ax.plot([ct, ev], [y, y], color=GRID_C, lw=2, zorder=1)
                ax.scatter([ct], [y], s=44, facecolor=SURFACE, edgecolor=SERIES[1], lw=2, zorder=2)
                ax.scatter([ev], [y], s=48, color=SERIES[0], edgecolor=SURFACE, lw=1.5, zorder=3)
                lab = f"{int(row.events_with)}/{int(row.events)}, p {float(row.p_stratified_permutation):.2g}"
                if pd.notna(row.get("watched_at_release", np.nan)):
                    lab += f", watched {int(row.watched_at_release)}"
                ax.annotate(lab, xy=(1.02, y), xycoords=("axes fraction", "data"), va="center", fontsize=7.3,
                            color=INK if arm in ("G1", "G0") else MUTED, annotation_clip=False)
            ax.set_xlim(0, 1)
            ax.set_ylim(-0.6, len(arms) - 0.4)
            ax.grid(axis="x", color=GRID_C, lw=0.6)
            ax.set_title(f"BTC, {per}: {cell}", loc="left", fontsize=9, color=INK)
            if r == 1:
                ax.set_xlabel("share of windows with an up alarm within 6 hours")
    axes[0, 0].set_yticks(ys, ["P1 unguarded", "G0 guarded, eta 0, keep", "G1 guarded, eta 0.25, rank", "S2 rank"])
    fig.suptitle("E24 [REAL, post hoc, exploratory]: E22's and E23's release designs rerun with the guarded P1 "
                 "(filled: release windows; open: matched controls; raw one-sided p)", x=0.01, ha="left", fontsize=10,
                 color=INK)
    fig.subplots_adjust(left=0.12, right=0.9, top=0.89, bottom=0.1, wspace=0.55, hspace=0.35)
    fig.savefig(FIG / "e24_real_release_rates.png", dpi=130)
    plt.close(fig)


def fig_power(plt):
    g = pd.read_csv(OUT / "e24_future_power_grid.csv")
    cand = pd.read_csv(OUT / "e24_future_candidates.csv")
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    i = 0
    for arm in FUTURE_ARMS:
        for d in FUTURE_DESIGNS:
            h = g[(g.arm == arm) & (g.design == d)].sort_values("end")
            x = pd.to_datetime(h.end)
            chosen = bool(cand[(cand.arm == arm) & (cand.design == d)].chosen.iloc[0])
            ax.plot(x, h.power, color=SERIES[i], marker=MARKERS[i], ms=3.5, lw=2.2 if chosen else 1.3,
                    label=f"{arm}, {'FOMC and CPI pooled' if d == 'fomc+cpi' else 'FOMC alone'}" + (" (chosen)" if chosen else ""))
            i += 1
    ax.axhline(POWER_TARGET, color=MUTED, ls="--", lw=1)
    ch = cand[cand.chosen].iloc[0]
    if isinstance(ch.end_star, str):
        ax.axvline(pd.Timestamp(ch.end_star), color=MUTED, ls=":", lw=1)
        ax.annotate(f"end date {ch.end_star}", xy=(pd.Timestamp(ch.end_star), 0.05), xytext=(4, 0),
                    textcoords="offset points", fontsize=7.5, color=MUTED)
    ax.set_ylim(0, 1)
    ax.grid(color=GRID_C, lw=0.6)
    ax.set_ylabel("power, one-sided level 0.05")
    ax.set_xlabel("end of the forward period (releases from 28 October 2026)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("E24 forward pre-registration [SYNTHETIC; planning inputs: REAL BTC rates of 2024 to 2026]",
                 loc="left", fontsize=9.5, color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "e24_future_power.png", dpi=130)
    plt.close(fig)


def part_figs():
    plt = _plt()
    FIG.mkdir(parents=True, exist_ok=True)
    for f, need in ((fig_near_blind, "e24_near_blind.csv"), (fig_synthetic, "e24_synthetic.csv"),
                    (fig_real, "e24_real_contrast.csv"), (fig_power, "e24_future_power_grid.csv")):
        if (OUT / need).exists():
            f(plt)
            print(f"{f.__name__}: done")
        else:
            print(f"{f.__name__}: {need} missing, skipped")


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
    elif part == "synthetic_summary":  # re-summarise the saved per-path rows (no simulation)
        s = summarise_synthetic(pd.read_csv(OUT / "e24_synthetic_paths.csv.gz", low_memory=False))
        s.to_csv(OUT / "e24_synthetic.csv", index=False)
    elif part == "real":
        part_real(procs)
    elif part == "e19blind":
        part_e19blind(procs)
    elif part == "power":
        part_power(procs)
    elif part == "future":
        part_future(procs)
    elif part == "future_gate":
        part_future_gate(procs)
    elif part == "figs":
        part_figs()
    else:
        raise SystemExit(__doc__)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
