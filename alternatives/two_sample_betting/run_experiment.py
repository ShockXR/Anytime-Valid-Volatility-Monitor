"""Experiments for the conformal rank betting detector (alternative A4, two-sample betting).

Research questions
  P   (design, development seeds only) Which histogram memory, prior and mixture
      weights to use? Chosen on synthetic development paths (seeds 1000+), then
      frozen for every reported run. Nothing is tuned on real data.
  Q1  [SYNTHETIC] Validity. With no change, is the share of paths with a false
      alarm within one day (1440 steps) at most 1440 / A = 5%, and within seven
      days at most 35%, for Gaussian, Student t4 and jump returns? How does the
      plug-in mean-based monitor (VolatilityMonitor, band 0 and 0.25) compare?
  Q1b [SYNTHETIC] Volatility clustering (null broken for both designs): GARCH(1,1)
      with persistence 0.95 and 0.99. Share of days with a false alarm, with and
      without a null band.
  Q1c [SYNTHETIC] Why the reference must grow. Ranking each new return against a
      FIXED reference window (the "quantile bins learned from the reference"
      design) gives p-values that are uniform one at a time but dependent; does
      a learning bettor then break the false alarm guarantee over two weeks?
  Q2  [SYNTHETIC] Power. Median delay (with a distribution-free 95% interval)
      and detection share for variance increases and decreases, a change of
      tails at fixed variance (Gaussian to variance-matched Student t4 and t3),
      a change of jump intensity at fixed total variance, and a change of skew
      with the law of |u| held fixed (only the signed score can see it).
  Q3  [SYNTHETIC] Does truncation matter when the statistic is a rank? Same
      detector on the score min(|u|, tau) with a fixed tau = 3 (ties above tau).
  Q5  [SYNTHETIC] The repo's rigorous mean-based variant (VolatilityMonitor with
      baseline="cs", a confidence sequence for the baseline, guarantee delta + m / A)
      on the same paths as Q1 (Gaussian null) and five Q2 scenarios.
  Q4  [REAL] Binance BTCUSDT and ETHUSDT one-minute returns, 1 Jul to 31 Dec
      2025, divided by the predictable intraday profile of the main package.
      Alarm counts of the rank detector against the main monitor (band 0.25),
      overlap within 60 minutes, and what the alarms raised only by the rank
      detector look like (scale ratio, scale-free shape distance, zero-return
      share, tail share, sign balance), against random times.

Run from the repository root: ``python alternatives/two_sample_betting/run_experiment.py [parts]``
with parts among ``pilot q1 q1b q2 q3 q4 figs`` (default: all but pilot).
Seeds are fixed per part. At most two worker processes.
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import rankbet as R  # noqa: E402
import tsb_sim as S  # noqa: E402

ROOT = HERE.parents[1]
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

DAY = 1440
A = DAY / 0.05  # 28,800: at most 5% chance of a false alarm per day under the null
CAL = 1440

FROZEN = R.RankBetConfig()  # C0 of the pilot (results/p_pilot_design.csv), frozen before any reported run

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def wilson(k: int, n: int, z: float = 1.959964):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, max(0.0, centre - half), min(1.0, centre + half)


def median_ci(x: np.ndarray, level: float = 0.95):
    """Median with a distribution-free interval from binomial order statistics."""
    x = np.sort(np.asarray(x, float))
    n = x.size
    if n == 0:
        return np.nan, np.nan, np.nan
    lo_k = int(stats.binom.ppf((1 - level) / 2, n, 0.5))
    hi_k = int(stats.binom.isf((1 - level) / 2, n, 0.5))
    lo = x[max(lo_k - 1, 0)]
    hi = x[min(hi_k, n - 1)]
    return float(np.median(x)), float(lo), float(hi)


def rank_alarms(r: np.ndarray, cfg: R.RankBetConfig, rng, n_ref: int = CAL, tau: float | None = None):
    s = R.score(r, cfg.kind)
    if tau is not None:
        s = np.minimum(s, tau)
    p = R.pvalues_against_reference(s, n_ref, rng)
    return R.run_rank_detector(p, A, cfg)


def summarize_null(alarm: np.ndarray, horizons=(DAY, 7 * DAY)) -> dict:
    out = {}
    n = alarm.size
    for h in horizons:
        k = int(((alarm >= 0) & (alarm < h)).sum())
        p, lo, hi = wilson(k, n)
        d = h // DAY
        out[f"fa_{d}d"] = p
        out[f"fa_{d}d_lo"] = lo
        out[f"fa_{d}d_hi"] = hi
        out[f"bound_{d}d"] = h / A
    return out


def summarize_change(alarm: np.ndarray, nu: int) -> dict:
    n = alarm.size
    early = (alarm >= 0) & (alarm < nu)
    det = alarm >= nu
    delays = alarm[det] - nu + 1
    med, lo, hi = median_ci(delays)
    fa, fa_lo, fa_hi = wilson(int(early.sum()), n)
    return dict(false_alarm_before_change=fa, fa_lo=fa_lo, fa_hi=fa_hi, detected=float(det.mean()),
                median_delay=med, median_delay_lo=lo, median_delay_hi=hi,
                mean_delay=float(delays.mean()) if delays.size else np.nan)


# ----------------------------------------------------------------------------
# P. Design pilot on development seeds
# ----------------------------------------------------------------------------

PILOT_CONFIGS = {
    "C0 default hl(240,1440) prior50 w(.25,.25,.5)": R.RankBetConfig(),
    "C1 short memory hl(60,240)": R.RankBetConfig(half_lives=(60.0, 240.0)),
    "C2 long memory hl(1440,5760)": R.RankBetConfig(half_lives=(1440.0, 5760.0)),
    "C3 prior 200": R.RankBetConfig(prior=200.0),
    "C4 weights (1/3,1/3,1/3)": R.RankBetConfig(weights=(1 / 3, 1 / 3, 1 / 3)),
    "C5 scale only (.5,.5,0)": R.RankBetConfig(weights=(0.5, 0.5, 0.0)),
    "C6 histogram only (0,0,1)": R.RankBetConfig(weights=(0.0, 0.0, 1.0)),
}
PILOT_CHANGES = {
    "sd x1.25 (Gaussian)": ({"law": "gauss"}, {"law": "gauss", "sd": 1.25}),
    "sd x0.8 (Gaussian)": ({"law": "gauss"}, {"law": "gauss", "sd": 0.8}),
    "Gaussian to t4": ({"law": "gauss"}, {"law": "t", "df": 4}),
    "jumps p .002 to .01 (J=6)": ({"law": "jumps", "p": 0.002, "J": 6.0}, {"law": "jumps", "p": 0.01, "J": 6.0}),
}


def part_pilot():
    rows = []
    rng = np.random.default_rng(1000)
    nu, post, reps = 500, 3000, 200
    for name, (pre, postspec) in PILOT_CHANGES.items():
        r = S.change_paths(reps, CAL, nu, post, pre, postspec, rng)
        p = R.pvalues_against_reference(np.abs(r), CAL, rng)
        for cname, cfg in PILOT_CONFIGS.items():
            res = R.run_rank_detector(p, A, cfg)
            rows.append(dict(scenario=name, config=cname, **summarize_change(res.alarm, nu)))
        print(name, flush=True)
    for name, r in {
        "null Gaussian (1 day)": S.draw({"law": "gauss"}, (400, CAL + DAY), np.random.default_rng(1001)),
        "GARCH persistence 0.95 (1 day)": S.garch_paths(400, CAL + DAY, 0.95, 0.05, np.random.default_rng(1002)),
    }.items():
        p = R.pvalues_against_reference(np.abs(r), CAL, rng)
        for cname, cfg in PILOT_CONFIGS.items():
            res = R.run_rank_detector(p, A, cfg)
            k = int((res.alarm >= 0).sum())
            fa, lo, hi = wilson(k, res.alarm.size)
            rows.append(dict(scenario=name, config=cname, false_alarm_1d=fa, fa_lo=lo, fa_hi=hi))
        print(name, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "p_pilot_design.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Detectors compared in the synthetic parts
# ----------------------------------------------------------------------------

RANK_OMNI = FROZEN  # C0 of the pilot: scale up 0.25, scale down 0.25, histogram 0.5, half-lives 240 and 1440
RANK_SCALE = R.RankBetConfig(weights=(0.5, 0.5, 0.0))
RANK_SIGNED = R.RankBetConfig(kind="signed")


def detectors_on(r: np.ndarray, rng, cal: int = CAL, which=("main0", "main25", "scale", "omni", "signed"),
                 tau: float | None = None) -> dict:
    """First alarms (monitored step index) of every detector on the same paths."""
    out = {}
    if "main0" in which:
        out["MAIN plug-in, band 0"] = S.main_monitor_first_alarm(r, A, cal_len=cal, null_band=0.0)[0]
    if "main25" in which:
        out["MAIN plug-in, band 0.25"] = S.main_monitor_first_alarm(r, A, cal_len=cal, null_band=0.25)[0]
    need_abs = [k for k in ("scale", "omni", "omni25", "hist") if k in which]
    if need_abs:
        s = np.abs(r) if tau is None else np.minimum(np.abs(r), tau)
        p = R.pvalues_against_reference(s, cal, rng)
        cfgs = {"scale": ("RANK scale only (abs), b 0", RANK_SCALE),
                "omni": ("RANK omnibus (abs), b 0", RANK_OMNI),
                "omni25": ("RANK omnibus (abs), b 0.25", R.RankBetConfig(band=0.25)),
                "hist": ("RANK histogram only (abs), b 0", R.RankBetConfig(weights=(0.0, 0.0, 1.0)))}
        for k in need_abs:
            name, cfg = cfgs[k]
            out[name] = R.run_rank_detector(p, A, cfg).alarm
    if "signed" in which:
        p = R.pvalues_against_reference(r, cal, rng)
        out["RANK omnibus (signed), b 0"] = R.run_rank_detector(p, A, RANK_SIGNED).alarm
    return out


# ----------------------------------------------------------------------------
# Q1. Validity with no change
# ----------------------------------------------------------------------------

Q1_LAWS = {
    "Gaussian": {"law": "gauss"},
    "Student t4 (unit variance)": {"law": "t", "df": 4},
    "Gaussian + jumps (0.5% per step, 8 sd)": {"law": "jumps", "p": 0.005, "J": 8.0},
}


def _q1_job(args):
    law, seed, cal, days, reps = args
    rng = np.random.default_rng(seed)
    r = S.draw(Q1_LAWS[law], (reps, cal + days * DAY), rng)
    which = ("main0", "main25", "scale", "omni", "signed", "omni25")
    res = detectors_on(r, rng, cal=cal, which=which)
    rows = []
    for det, alarm in res.items():
        rows.append(dict(law=law, calibration=cal, paths=reps, detector=det,
                         **summarize_null(alarm, horizons=tuple(h * DAY for h in sorted({1, days})))))
    return rows


def part_q1():
    jobs = [(law, 100 + i, CAL, 7, 1000) for i, law in enumerate(Q1_LAWS)]
    jobs += [(law, 200 + i, 120, 1, 1000) for i, law in enumerate(Q1_LAWS)]  # short calibration: 2 hours
    with Pool(2) as pool:
        rows = [r for rr in pool.map(_q1_job, jobs) for r in rr]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q1_validity.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Q1b. Volatility clustering (GARCH), where the exchangeability null is false
# ----------------------------------------------------------------------------


def _q1b_job(args):
    persistence, alpha, seed, reps = args
    rng = np.random.default_rng(seed)
    r = S.garch_paths(reps, CAL + DAY, persistence, alpha, rng)
    res = detectors_on(r, rng, which=("main0", "main25", "scale", "omni", "omni25", "hist"))
    rows = []
    for det, alarm in res.items():
        k = int((alarm >= 0).sum())
        p, lo, hi = wilson(k, reps)
        rows.append(dict(persistence=persistence, alpha=alpha, half_life_steps=np.log(0.5) / np.log(persistence),
                         paths=reps, detector=det, fa_1d=p, fa_1d_lo=lo, fa_1d_hi=hi, bound_1d=DAY / A))
    return rows


def part_q1b():
    jobs = [(0.95, 0.05, 300, 500), (0.99, 0.04, 301, 500)]
    with Pool(2) as pool:
        rows = [r for rr in pool.map(_q1b_job, jobs) for r in rr]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q1b_garch.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Q1c. Why the reference must grow: ranks against a FIXED reference window
# ----------------------------------------------------------------------------


def pvalues_fixed_reference(s: np.ndarray, n_ref: int, rng) -> np.ndarray:
    """p-values of s[:, n_ref:] ranked against the first n_ref values only (the reference never grows)."""
    ref = np.sort(s[:, :n_ref], axis=1)
    mon = s[:, n_ref:]
    out = np.empty(mon.shape)
    for i in range(s.shape[0]):
        hi = np.searchsorted(ref[i], mon[i], side="right")
        lo = np.searchsorted(ref[i], mon[i], side="left")
        out[i] = (n_ref - hi + rng.random(mon.shape[1]) * (hi - lo + 1)) / (n_ref + 1)
    return out


def _q1c_job(args):
    n_ref, seed, reps, days = args
    rng = np.random.default_rng(seed)
    r = S.draw(G, (reps, n_ref + days * DAY), rng)
    rows = []
    for label, fn in (("growing reference (conformal, exact)", R.pvalues_against_reference),
                      ("fixed reference window (not a martingale)", pvalues_fixed_reference)):
        p = fn(np.abs(r), n_ref, rng)
        for det, cfg in (("RANK omnibus (abs), b 0", RANK_OMNI), ("RANK histogram only (abs), b 0", R.RankBetConfig(weights=(0.0, 0.0, 1.0)))):
            alarm = R.run_rank_detector(p, A, cfg).alarm
            rows.append(dict(reference=label, n_ref=n_ref, paths=reps, detector=det,
                             **summarize_null(alarm, horizons=(DAY, 7 * DAY, days * DAY))))
    return rows


def part_q1c():
    jobs = [(240, 150, 500, 14), (1440, 151, 500, 14)]
    with Pool(2) as pool:
        rows = [r for rr in pool.map(_q1c_job, jobs) for r in rr]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q1c_fixed_reference.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Q2. Power: variance, tails, jump intensity, skew
# ----------------------------------------------------------------------------

G = {"law": "gauss"}
Q2_SCENARIOS = {
    "sd x1.25 (variance x1.56)": (G, {"law": "gauss", "sd": 1.25}),
    "sd x1.5 (variance x2.25)": (G, {"law": "gauss", "sd": 1.5}),
    "sd x0.8 (variance x0.64)": (G, {"law": "gauss", "sd": 0.8}),
    "sd x0.67 (variance x0.44)": (G, {"law": "gauss", "sd": 0.67}),
    "tails: Gaussian to t4, same variance": (G, {"law": "t", "df": 4}),
    "tails: Gaussian to t3, same variance": (G, {"law": "t", "df": 3}),
    "tails: t4 to Gaussian, same variance": ({"law": "t", "df": 4}, G),
    "jumps: 0.2% to 0.5% per step, same variance": ({"law": "jumps", "p": 0.002, "J": 6.0},
                                                     {"law": "jumps", "p": 0.005, "J": 6.0}),
    "jumps: 0.2% to 1% per step, same variance": ({"law": "jumps", "p": 0.002, "J": 6.0},
                                                   {"law": "jumps", "p": 0.01, "J": 6.0}),
    "skew: large moves 80% negative, same |u| law": (G, {"law": "signskew", "thr": 1.5, "p_neg": 0.8}),
}
NU, POST, Q2_REPS = 500, 4000, 400


def _q2_job(args):
    name, seed = args
    pre, post = Q2_SCENARIOS[name]
    rng = np.random.default_rng(seed)
    r = S.change_paths(Q2_REPS, CAL, NU, POST, pre, post, rng)
    res = detectors_on(r, rng, which=("main0", "scale", "omni", "signed"))
    return [dict(scenario=name, paths=Q2_REPS, change_at=NU, horizon_after_change=POST, detector=det,
                 **summarize_change(alarm, NU)) for det, alarm in res.items()]


def part_q2():
    jobs = [(name, 400 + i) for i, name in enumerate(Q2_SCENARIOS)]
    with Pool(2) as pool:
        rows = [r for rr in pool.map(_q2_job, jobs) for r in rr]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q2_power.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Q3. Does truncation matter for a rank statistic?
# ----------------------------------------------------------------------------

Q3_SCENARIOS = ["sd x1.25 (variance x1.56)", "tails: Gaussian to t4, same variance",
                "jumps: 0.2% to 1% per step, same variance"]


def _q3_job(args):
    name, seed = args
    rows = []
    if name.startswith("null"):
        law = {"null Gaussian": G, "null jumps (0.5%, 8 sd)": {"law": "jumps", "p": 0.005, "J": 8.0}}[name]
        r = S.draw(law, (1000, CAL + DAY), np.random.default_rng(seed))
        for label, tau in (("score |u|", None), ("score min(|u|, 3)", 3.0)):
            alarm = detectors_on(r, np.random.default_rng(seed + 1), which=("omni",), tau=tau)["RANK omnibus (abs), b 0"]
            rows.append(dict(scenario=name, score=label, paths=1000, **summarize_null(alarm, horizons=(DAY,))))
        return rows
    pre, post = Q2_SCENARIOS[name]
    r = S.change_paths(Q2_REPS, CAL, NU, POST, pre, post, np.random.default_rng(seed))
    for label, tau in (("score |u|", None), ("score min(|u|, 3)", 3.0)):
        alarm = detectors_on(r, np.random.default_rng(seed + 1), which=("omni",), tau=tau)["RANK omnibus (abs), b 0"]
        rows.append(dict(scenario=name, score=label, paths=Q2_REPS, **summarize_change(alarm, NU)))
    return rows


def part_q3():
    jobs = [(name, 500 + i) for i, name in enumerate(Q3_SCENARIOS)] + [("null Gaussian", 510), ("null jumps (0.5%, 8 sd)", 512)]
    with Pool(2) as pool:
        rows = [r for rr in pool.map(_q3_job, jobs) for r in rr]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q3_truncation.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Q5. The mean-based monitor made rigorous: confidence-sequence baseline on the same paths
# ----------------------------------------------------------------------------


def cs_first_alarm(r: np.ndarray) -> np.ndarray:
    """First alarm (monitored step index) of VolatilityMonitor(baseline="cs") per path, -1 if none."""
    from avmonitor.monitor import VolatilityMonitor

    out = np.full(r.shape[0], -1)
    for i, row in enumerate(r):
        mon = VolatilityMonitor(threshold=A, cal_len=CAL, baseline="cs")
        for t, x in enumerate(row):
            ev = mon.update(x)
            if ev is not None:
                out[i] = ev.index - ev.monitored_since
                break
    return out


def _q5_job(args):
    kind, name, seed = args
    rng = np.random.default_rng(seed)
    if kind == "null":  # same paths as Q1 (same seed, drawn first)
        r = S.draw(Q1_LAWS[name], (1000, CAL + 7 * DAY), rng)
        return [dict(part="Q1 paths", scenario=name, detector="MAIN confidence-sequence baseline (delta 0.01), band 0",
                     **summarize_null(cs_first_alarm(r), horizons=(DAY, 7 * DAY)))]
    pre, post = Q2_SCENARIOS[name]  # same paths as Q2
    r = S.change_paths(Q2_REPS, CAL, NU, POST, pre, post, rng)
    return [dict(part="Q2 paths", scenario=name, detector="MAIN confidence-sequence baseline (delta 0.01), band 0",
                 **summarize_change(cs_first_alarm(r), NU))]


def part_q5():
    q2names = list(Q2_SCENARIOS)
    jobs = [("null", "Gaussian", 100)]
    jobs += [("change", nm, 400 + q2names.index(nm)) for nm in (
        "sd x1.25 (variance x1.56)", "sd x0.8 (variance x0.64)", "tails: Gaussian to t4, same variance",
        "jumps: 0.2% to 1% per step, same variance", "skew: large moves 80% negative, same |u| law")]
    with Pool(2) as pool:
        rows = [r for rr in pool.map(_q5_job, jobs) for r in rr]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q5_cs_baseline.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Q4. Real data: Binance BTCUSDT and ETHUSDT, Jul to Dec 2025
# ----------------------------------------------------------------------------

COINS = {"BTC": "BTCUSDT", "ETH": "ETHUSDT"}
HIST, START, END = "2025-04-01", "2025-07-01", "2026-01-01"
MATCH = 60  # minutes
WIN = 240  # evidence window before an alarm, minutes
N_RANDOM = 2000
# Detectors fixed before looking at any real alarm. "nonzero only" skips exact-zero returns (frequent in
# BTCUSDT one-minute data, a data property checked before the run, not an alarm outcome; the share is in
# q4_real_summary.csv): the subsequence of nonzero returns of an exchangeable sequence is exchangeable,
# so this variant is just as valid.
Q4_RANK = {
    "RANK omnibus (abs), b 0": (R.RankBetConfig(), False),
    "RANK omnibus (abs), b 0.25": (R.RankBetConfig(band=0.25), False),
    "RANK omnibus (abs, nonzero only), b 0.25": (R.RankBetConfig(band=0.25), True),
    "RANK omnibus (signed), b 0.25": (R.RankBetConfig(kind="signed", band=0.25), False),
}
SHORT_REF = "RANK omnibus (abs), b 0.25, 120-minute reference"  # validity does not need a long reference


def _intervals(segment_starts, events, n) -> list:
    """Monitored intervals [start, end] (inclusive) of a monitor; the rest of the time it was calibrating."""
    ends = {e.monitored_since: e.index for e in events}
    return [(s0, ends.get(s0, n - 1)) for s0 in segment_starts]


def _monitored_share(intervals, n) -> float:
    return sum(b - a + 1 for a, b in intervals) / n
REFERENCE_MAIN = "MAIN plug-in, band 0.25"


def q4_prepare(coin: str):
    from avmonitor.data import intraday_scale, load_binance

    df = load_binance(COINS[coin], HIST, END, download=False)
    s = intraday_scale(df.logret, lookback_days=20)
    u = (df.logret / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    i0 = int(np.searchsorted(df.index, pd.Timestamp(START, tz="UTC")))
    return df.index, u, i0, df.attrs["checks"]


def window_metrics(u, ref_a, ref_b, a, b) -> dict:
    ref = u[max(ref_a, 0):ref_b]
    ref = ref[np.isfinite(ref)]
    win = u[max(a, 0):b]
    win = win[np.isfinite(win)]
    nan = dict(scale_ratio=np.nan, trunc_var_ratio=np.nan, shape_ks=np.nan, zero_ref=np.nan, zero_win=np.nan,
               tail_ref=np.nan, tail_win=np.nan, pos_share_win=np.nan)
    if ref.size < 100 or win.size < 30:
        return nan
    mr, mw = np.median(np.abs(ref)), np.median(np.abs(win))
    if mr <= 0 or mw <= 0:
        return nan
    tau2 = (3.0 * mr * 1.482602218505602) ** 2
    nz = win[win != 0]
    return dict(
        scale_ratio=mw / mr,
        trunc_var_ratio=float(np.mean(np.minimum(win**2, tau2)) / np.mean(np.minimum(ref**2, tau2))),
        shape_ks=float(stats.ks_2samp(np.abs(win) / mw, np.abs(ref) / mr).statistic),
        zero_ref=float(np.mean(ref == 0)), zero_win=float(np.mean(win == 0)),
        tail_ref=float(np.mean(np.abs(ref) / mr > 3 / 0.6744897501960817)),
        tail_win=float(np.mean(np.abs(win) / mw > 3 / 0.6744897501960817)),
        pos_share_win=float(np.mean(nz > 0)) if nz.size else np.nan,
    )


def _q4_job(coin: str):
    from avmonitor.monitor import VolatilityMonitor

    idx, u, i0, checks = q4_prepare(coin)
    x = u[i0:]
    events, intervals = {}, {}
    for band in (0.0, 0.25):
        mon = VolatilityMonitor(threshold=A, c=3.0, cal_len=CAL, null_band=band)
        ev = mon.run(x)
        intervals[f"MAIN plug-in, band {band:g}"] = [(a + i0, b + i0) for a, b in
                                                     _intervals([sg["start"] for sg in mon.segments], ev, x.size)]
        events[f"MAIN plug-in, band {band:g}"] = [dict(index=e.index + i0, monitored_since=e.monitored_since + i0,
                                                       cal_start=e.monitored_since + i0 - CAL, direction=e.direction,
                                                       share_up=np.nan, share_down=np.nan, share_hist=np.nan,
                                                       mean_p_recent=np.nan) for e in ev]
    for name, (cfg, nonzero) in list(Q4_RANK.items()) + [(SHORT_REF, (R.RankBetConfig(band=0.25), False))]:
        mon = R.RankMonitor(threshold=A, cal_len=120 if name == SHORT_REF else CAL, config=cfg, seed=4)
        xx = np.where(x == 0, np.nan, x) if nonzero else x
        ev = mon.run(xx)
        intervals[name] = [(a + i0, b + i0) for a, b in _intervals([sg["start"] for sg in mon.segments], ev, x.size)]
        rows = []
        for e in ev:
            sh = np.array([e.share_up, e.share_down, e.share_hist])
            rows.append(dict(index=e.index + i0, monitored_since=e.monitored_since + i0, cal_start=e.cal_start + i0,
                             direction=("up", "down", "shape")[int(np.argmax(sh))], share_up=e.share_up,
                             share_down=e.share_down, share_hist=e.share_hist, mean_p_recent=e.mean_p_recent,
                             **{f"dr{k}": v for k, v in enumerate(e.density_ratio)}))
        events[name] = rows
    return coin, idx, u, i0, checks, events, intervals


PARTNER = {"MAIN": "RANK omnibus (abs), b 0.25", "RANK": REFERENCE_MAIN}


def _inside(t: np.ndarray, iv: list) -> np.ndarray:
    """True where time t lies in one of the (sorted, disjoint) monitored intervals."""
    if not iv:
        return np.zeros(t.size, bool)
    a = np.array([x[0] for x in iv])
    b = np.array([x[1] for x in iv])
    k = np.searchsorted(a, t, side="right") - 1
    return (k >= 0) & (t <= b[np.clip(k, 0, None)])


def _agg(col: pd.Series) -> float:
    """Median of a numeric column; share of True for a boolean one (used for the weekend share)."""
    return float(col.mean()) if col.dtype == bool else float(col.median())


def part_q4():
    with Pool(2) as pool:
        results = pool.map(_q4_job, list(COINS))
    rng = np.random.default_rng(44)
    ev_rows, sum_rows, ctrl_rows = [], [], []
    idx_of = {r[0]: r[1] for r in results}
    for coin, idx, u, i0, checks, events, intervals in results:
        n = u.size
        days = (n - i0) / DAY
        lo_t, hi_t = i0 + 2 * DAY + WIN, n - DAY - 1
        seg = u[i0:]
        seg = seg[np.isfinite(seg)]
        for det, evs in events.items():
            partner = PARTNER[det.split()[0]]
            pt = np.array([e["index"] for e in events[partner]], int)
            t = np.array([e["index"] for e in evs], int)
            dist = np.array([np.min(np.abs(pt - tt)) if pt.size else np.inf for tt in t])
            p_mon = _inside(t, intervals[partner])
            for e, d, pm in zip(evs, dist, p_mon):
                tt = e["index"]
                a = max(e["monitored_since"], tt - WIN + 1)
                m = window_metrics(u, e["cal_start"], e["monitored_since"], a, tt + 1)
                post = window_metrics(u, e["cal_start"], e["monitored_since"], tt + 1, tt + 1 + DAY)
                ev_rows.append(dict(coin=coin, detector=det, time=idx[tt], index=tt, delay_minutes=tt - e["monitored_since"],
                                    label=e["direction"], share_up=e["share_up"], share_down=e["share_down"],
                                    share_hist=e["share_hist"], mean_p_recent=e["mean_p_recent"],
                                    **{f"dr{k}": e.get(f"dr{k}", np.nan) for k in range(len(R.DEFAULT_EDGES) - 1)},
                                    partner=partner, minutes_to_nearest_partner=float(d), matched=bool(d <= MATCH),
                                    partner_monitoring=bool(pm),
                                    **m, **{f"post_{k}": v for k, v in post.items()}))
            # chance overlap: shift this detector's alarms by random offsets of 1 to 30 days (circular)
            span = n - i0
            chance = []
            for _ in range(200):
                sh = (t - i0 + rng.integers(DAY, 30 * DAY)) % span + i0
                dd = np.array([np.min(np.abs(pt - y)) if pt.size else np.inf for y in sh])
                chance.append(np.mean(dd <= MATCH) if dd.size else np.nan)
            labels = [e["direction"] for e in evs]
            unmatched = dist > MATCH
            mshare = _monitored_share(intervals[det], n - i0)
            sum_rows.append(dict(coin=coin, detector=det, days=days, zero_return_share=float(np.mean(seg == 0)),
                                 monitored_share=mshare, alarms=t.size,
                                 alarms_per_day=t.size / days, alarms_per_monitored_day=t.size / (days * mshare),
                                 up=labels.count("up"), down=labels.count("down"), shape=labels.count("shape"),
                                 partner=partner, matched_within_60min=int(np.sum(~unmatched)),
                                 matched_share=float(np.mean(~unmatched)) if t.size else np.nan,
                                 chance_matched_share=float(np.nanmean(chance)) if t.size else np.nan,
                                 unmatched=int(unmatched.sum()),
                                 unmatched_partner_monitoring=int(np.sum(unmatched & p_mon)),
                                 unmatched_partner_calibrating=int(np.sum(unmatched & ~p_mon))))
        # random times with a reference that ended as long before as a typical alarm's calibration did
        delays = [e["index"] - e["monitored_since"] for e in events["RANK omnibus (abs), b 0.25"]]
        for _ in range(N_RANDOM):
            d = int(rng.choice(delays)) if delays else 0
            tt = int(rng.integers(max(lo_t, i0 + d + DAY + WIN), hi_t))
            since = tt - d
            a = max(since, tt - WIN + 1)
            m = window_metrics(u, since - CAL, since, a, tt + 1)
            post = window_metrics(u, since - CAL, since, tt + 1, tt + 1 + DAY)
            ctrl_rows.append(dict(coin=coin, index=tt, delay_minutes=d, **m, **{f"post_{k}": v for k, v in post.items()}))
        print(coin, checks, flush=True)
    ev = pd.DataFrame(ev_rows)
    ev["weekend"] = pd.to_datetime(ev.time, utc=True).dt.dayofweek >= 5  # Saturday or Sunday, UTC
    ev.to_csv(OUT / "q4_real_events.csv", index=False)
    only = ev[(ev.detector == "RANK omnibus (abs), b 0.25") & ~ev.matched & ev.partner_monitoring]
    only[["coin", "time", "label", "share_up", "share_down", "share_hist", "scale_ratio", "trunc_var_ratio", "shape_ks",
          "zero_ref", "zero_win", "tail_ref", "tail_win", "dr0", "dr1", "dr9", "dr10", "dr11", "weekend"]].to_csv(
        OUT / "q4_real_rank_only_alarms.csv", index=False)
    ctrl = pd.DataFrame(ctrl_rows)
    ctrl.to_csv(OUT / "q4_real_random_times.csv", index=False)
    summ = pd.DataFrame(sum_rows)
    summ.to_csv(OUT / "q4_real_summary.csv", index=False)
    # profile: medians of the window metrics by group
    cols = ["scale_ratio", "trunc_var_ratio", "shape_ks", "zero_ref", "zero_win", "tail_ref", "tail_win", "pos_share_win",
            "post_scale_ratio", "post_shape_ks", "post_zero_win", "weekend"] + [f"dr{k}" for k in range(len(R.DEFAULT_EDGES) - 1)]
    prof = []
    for (coin, det), g in ev.groupby(["coin", "detector"]):
        groups = {"all alarms": g, "matched (partner alarm within 60 min)": g[g.matched],
                  "unmatched": g[~g.matched],
                  "unmatched, partner was monitoring": g[(~g.matched) & g.partner_monitoring],
                  "unmatched, partner was calibrating": g[(~g.matched) & ~g.partner_monitoring]}
        if det.startswith("RANK"):
            for lab in ("up", "down", "shape"):
                groups[f"unmatched, label {lab}"] = g[(~g.matched) & (g.label == lab)]
                groups[f"unmatched, partner monitoring, label {lab}"] = g[(~g.matched) & g.partner_monitoring & (g.label == lab)]
        for gname, gg in groups.items():
            prof.append(dict(coin=coin, detector=det, group=gname, n=len(gg), **{c: _agg(gg[c]) if len(gg) else np.nan for c in cols}))
    for coin, g in ctrl.groupby("coin"):
        g = g.assign(weekend=pd.to_datetime(idx_of[coin][g["index"].to_numpy()]).dayofweek >= 5)
        prof.append(dict(coin=coin, detector="control: random times (matched delay)", group="random", n=len(g),
                         **{c: _agg(g[c]) if c in g else np.nan for c in cols}))
    prof = pd.DataFrame(prof)
    prof.to_csv(OUT / "q4_real_profiles.csv", index=False)
    print(summ.to_string())
    print(prof.to_string())


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

# Categorical slots 1 to 4 of the reference palette (fixed order), neutral grey for controls.
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def part_figs():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Figure 1: median delays in Q2
    d = pd.read_csv(OUT / "q2_power.csv")
    dets = ["MAIN plug-in, band 0", "RANK scale only (abs), b 0", "RANK omnibus (abs), b 0", "RANK omnibus (signed), b 0"]
    scen = list(dict.fromkeys(d.scenario))
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    h = 0.19
    for j, det in enumerate(dets):
        g = d[d.detector == det].set_index("scenario").loc[scen]
        y = np.arange(len(scen)) + (j - 1.5) * h
        med = g.median_delay.to_numpy()
        lo = med - g.median_delay_lo.to_numpy()
        hi = g.median_delay_hi.to_numpy() - med
        ax.barh(y, med, height=h * 0.85, color=COLORS[j], label=det, edgecolor="white", linewidth=0.5)
        ax.errorbar(med, y, xerr=[lo, hi], fmt="none", ecolor=INK, elinewidth=0.7, capsize=1.5)
        for yy, top, det_share in zip(y, g.median_delay_hi.to_numpy(), g.detected.to_numpy()):
            if det_share < 0.9:
                ax.text(top * 1.06, yy, f"{det_share:.0%} detected", va="center", fontsize=6.5, color=MUTED)
    ax.set_yticks(np.arange(len(scen)))
    ax.set_yticklabels(scen, fontsize=8, color=INK)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("median detection delay in steps (log scale), 95% interval; shorter is better", fontsize=8, color=MUTED)
    ax.set_title("[SYNTHETIC] Q2: delay after a change, same alarm budget (A = 28,800)", fontsize=10, color=INK, loc="left")
    ax.set_xlim(right=d.median_delay_hi.max() * 2.2)
    _style(ax)
    fig.legend(*ax.get_legend_handles_labels(), fontsize=7.5, frameon=False, loc="lower center", ncol=2)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(OUT / "fig1_synthetic_delays.png", dpi=150)
    plt.close(fig)

    # Figure 2: what the real alarms look like
    ev = pd.read_csv(OUT / "q4_real_events.csv")
    ctrl = pd.read_csv(OUT / "q4_real_random_times.csv")
    det = "RANK omnibus (abs), b 0.25"
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
    for ax, coin in zip(axes, COINS):
        c = ctrl[ctrl.coin == coin]
        ax.scatter(np.log2(c.scale_ratio), c.shape_ks, s=5, color="#b5b4ad", alpha=0.5, label="random times (control)",
                   linewidths=0)
        e = ev[(ev.coin == coin) & (ev.detector == det)]
        for k, (lab, m) in enumerate((("rank alarm, main monitor also alarms within 60 min", e.matched),
                                      ("rank alarm only", ~e.matched))):
            ax.scatter(np.log2(e.scale_ratio[m]), e.shape_ks[m], s=16, color=COLORS[k], label=lab,
                       edgecolors="white", linewidths=0.5)
        out = int(np.sum(np.abs(np.log2(c.scale_ratio)) > 3))
        ax.set_xlim(-3, 3)
        ax.text(0.99, 0.98, f"{out} of {len(c)} random times lie outside the x range", transform=ax.transAxes,
                ha="right", va="top", fontsize=7, color=MUTED)
        ax.set_title(f"[REAL] {COINS[coin]}, Jul to Dec 2025", fontsize=10, color=INK, loc="left")
        ax.set_xlabel("log2 robust scale ratio (240 min before alarm / calibration window)", fontsize=8, color=MUTED)
        ax.axvline(0, color=MUTED, linewidth=0.6)
        _style(ax)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
    axes[0].set_ylabel("scale-free shape distance (KS of |u| / median |u|)", fontsize=8, color=MUTED)
    axes[0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_real_alarm_profiles.png", dpi=150)
    plt.close(fig)

    # Figure 3: learned density ratio of new versus reference move sizes at the alarm
    edges = np.asarray(R.DEFAULT_EDGES)
    labels = [f"{100 * a:g} to {100 * b:g}" for a, b in zip(edges[:-1], edges[1:])]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, coin in zip(axes, COINS):
        e = ev[(ev.coin == coin) & (ev.detector == det)]
        groups = (("alarm matched by the main monitor within 60 min", e[e.matched]),
                  ("rank only, main monitor watching", e[~e.matched & e.partner_monitoring]))
        for k, (lab, g) in enumerate(groups):
            med = g[[f"dr{i}" for i in range(len(labels))]].median().to_numpy()
            ax.plot(np.arange(len(labels)), med, marker="o", markersize=4, linewidth=2, color=COLORS[k],
                    label=f"{lab} (n = {len(g)})")
        ax.axhline(1, color=MUTED, linewidth=0.8)
        ax.set_yscale("log")
        ax.set_yticks([0.25, 0.5, 1, 2, 4])
        ax.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
        ax.minorticks_off()
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=60, fontsize=7)
        ax.set_xlabel("rank bin: % of earlier moves that were larger (0 = largest moves, 100 = smallest moves)", fontsize=8,
                      color=MUTED)
        ax.set_title(f"[REAL] {COINS[coin]}: what changed at the alarm (median)", fontsize=10, color=INK, loc="left")
        _style(ax)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
    axes[0].set_ylabel("new frequency / reference frequency (log scale)", fontsize=8, color=MUTED)
    axes[0].legend(fontsize=7, frameon=False, loc="upper center")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_real_density_ratio.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

PARTS = {"pilot": part_pilot, "q1": part_q1, "q1b": part_q1b, "q1c": part_q1c, "q2": part_q2, "q3": part_q3, "q4": part_q4, "q5": part_q5, "figs": part_figs}

if __name__ == "__main__":
    todo = sys.argv[1:] or [k for k in PARTS if k != "pilot"]
    for k in todo:
        t0 = time.time()
        PARTS[k]()
        print(f"[{k}] done in {time.time() - t0:.0f} s", flush=True)
