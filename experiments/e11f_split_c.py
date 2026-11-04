"""E11f: a separate, tighter truncation constant for the decrease side.

docs/theory_delay.md (section 6.3) shows that truncation caps the evidence a decrease
detector can gather per step at log(B / (B - 1)) with B = tau^2 / m0, about 1 / c^2, and
predicts that c of about 1.5 to 2 on the decrease side (keeping c = 3 for increases) cuts
decrease delays by a quarter to a half. Each side is an e-detector for its own truncated
statistic, so with a known baseline, or with the confidence-sequence baseline of
docs/math.md section 6, the two-sided guarantee is unchanged. This script uses the plug-in
baseline, which has no proof (E4e), so its false alarm rates are empirical only. It checks
the prediction and the price on data with ties and zeros, which a small c makes more visible.

Part a [SYNTHETIC]: false alarm rate within one day (no change) and decrease delay
(sd x 0.8, 0.67, 0.5 at the first monitored step) for c_down = 3, 2 and 1.5, with the
plug-in baseline from a 1440-step calibration window, A = 28,800 (each side at 2A).
Return laws:
  * continuous: Gaussian, Student t4, Gaussian with jumps (600 paths per change);
  * rounded prices ("ticks"): a Gaussian return is first multiplied by its standard deviation
    (1 before the change, the ratio after it) and then rounded to a tick that is FIXED in units
    of the pre-change standard deviation (0.25, 0.5 and 0.75), as on a real price grid, so
    more returns become zero when volatility falls (2000 paths per change, the same normal
    draws for every tick size and every c_down).
  * For comparison only, the law used before review 2b (round first, then multiply by the
    ratio, so the tick shrinks with volatility and the share of zeros does not change),
    at a tick of 0.5, on the same normal draws ("ticks_0.5_scaled").
The continuous laws use the same random numbers as the first version of this script, so
their rows are unchanged.
Part b [REAL]: the five Binance coins, July to December 2025, band 0.25: alarm counts and
the share of alarms followed by a next day on the flagged side of the baseline (with the
number of alarms that have a following day, which is the denominator of that share).

Usage: python experiments/e11f_split_c.py [a] [b]
Outputs: results/e11f_split_c_synthetic.csv, results/e11f_split_c_real.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from avmonitor import data as D  # noqa: E402
from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402
from avmonitor.robust import MAD_NORMAL_SCALE  # noqa: E402

RES = ROOT / "results"
SEED = 20261103
CAL = 1440
A = 28_800.0
C_UP = 3.0
C_DOWNS = (3.0, 2.0, 1.5)
LAWS = ("gaussian", "t4", "jumps")  # continuous laws
TICKS = (0.25, 0.5, 0.75)  # tick sizes in units of the pre-change standard deviation
OLD_TICK = 0.5  # tick of the law used before review 2b (kept for comparison only)
RATIOS = (0.8, 0.67, 0.5)
HORIZON = 3000
N_NULL = 2000
N_CHANGE = 600  # continuous laws
N_CHANGE_TICKS = 2000  # tick laws (the result depends on the tick size, so more paths)


def draw(law: str, size, rng: np.random.Generator) -> np.ndarray:
    """Continuous return laws with unit variance."""
    if law == "gaussian":
        return rng.standard_normal(size)
    if law == "t4":
        return rng.standard_t(4, size) / np.sqrt(2.0)  # unit variance
    if law == "jumps":
        z = rng.standard_normal(size)
        return z + (rng.random(size) < 0.005) * rng.normal(0.0, 8.0, size)
    raise ValueError(law)


def round_to_tick(r: np.ndarray, tick: float) -> np.ndarray:
    return np.round(r / tick) * tick


def tick_paths(z: np.ndarray, ratio: float, tick: float, law: str = "fixed", cal: int = CAL) -> np.ndarray:
    """Rounded returns from standard normal draws z (rows = paths); sd 1 before step ``cal``, ``ratio`` after.

    ``law="fixed"``: multiply by the sd first, then round to ``tick`` (a fixed price grid: when
    volatility falls, more returns are zero; the zero share is 2 Phi(tick / (2 sd)) - 1).
    ``law="scaled"``: round first, then multiply by the sd (the law used before review 2b: the tick
    shrinks with volatility and the zero share stays at 2 Phi(tick / 2) - 1). Comparison only.
    """
    prof = np.ones(z.shape[1])
    prof[cal:] = ratio
    if law == "fixed":
        return round_to_tick(z * prof[None, :], tick)
    if law == "scaled":
        return round_to_tick(z, tick) * prof[None, :]
    raise ValueError(law)


def first_alarms(r: np.ndarray, c_down: float, band: float = 0.0, cal: int = CAL, threshold: float = A,
                 c_up: float = C_UP):
    """Vectorised first alarm of each side for paths in the rows of r (calibration = first ``cal`` steps).

    Same arithmetic as VolatilityMonitor with the plug-in baseline up to the first alarm (checked
    against the class in check_against_class and in tests/test_e11f_split_c.py). Returns the
    0-based monitored index of the first up and the first down alarm (-1 if none; only the side that
    alarms first is recorded)."""
    cal_r, mon = r[:, :cal], r[:, cal:]
    scale = np.median(np.abs(cal_r), axis=1) * MAD_NORMAL_SCALE
    tau2 = (c_up * scale) ** 2
    tau2_dn = (c_down * scale) ** 2
    m_up = np.mean(np.minimum(cal_r**2, tau2[:, None]), axis=1) * (1 + band)
    m_dn = np.mean(np.minimum(cal_r**2, tau2_dn[:, None]), axis=1) * (1 - band)
    b = tau2_dn / m_dn
    lam_dn = np.asarray(DEFAULT_DOWN_FRACTIONS, float)[None, :] / (b[:, None] - 1.0)
    up_grid = np.asarray(DEFAULT_UP_GRID, float)[None, :]
    n, T = mon.shape
    S_up = np.zeros((n, up_grid.shape[1]))
    S_dn = np.zeros((n, lam_dn.shape[1]))
    t_up = np.full(n, -1)
    t_dn = np.full(n, -1)
    for t in range(T):
        r2 = mon[:, t] ** 2
        x = np.minimum(r2, tau2)
        xd = np.minimum(r2, tau2_dn)
        S_up = (1.0 + up_grid * (x / m_up - 1.0)[:, None]) * (S_up + 1.0)
        S_dn = (1.0 + lam_dn * (1.0 - xd / m_dn)[:, None]) * (S_dn + 1.0)
        hit_up = (S_up.mean(axis=1) >= 2 * threshold) & (t_up < 0) & (t_dn < 0)
        hit_dn = (S_dn.mean(axis=1) >= 2 * threshold) & (t_dn < 0) & (t_up < 0) & ~hit_up
        t_up[hit_up] = t
        t_dn[hit_dn] = t
    return t_up, t_dn


def class_first_alarm(r: np.ndarray, c_down: float, cal: int = CAL, threshold: float = A):
    """First alarm of VolatilityMonitor on one path, as (monitored index, direction), or None."""
    ev = VolatilityMonitor(threshold=threshold, cal_len=cal, c_down=c_down).run(r)
    return (ev[0].index - cal, ev[0].direction) if ev else None


def vectorised_first_alarm(t_up: np.ndarray, t_dn: np.ndarray, i: int):
    if t_up[i] >= 0:
        return int(t_up[i]), "up"
    if t_dn[i] >= 0:
        return int(t_dn[i]), "down"
    return None


def check_against_class(rng: np.random.Generator) -> None:
    r = draw("jumps", (4, CAL + 3000), rng) * np.r_[np.ones(CAL), 0.6 * np.ones(3000)][None, :]
    for cd in (3.0, 1.5):
        t_up, t_dn = first_alarms(r, cd)
        for i in range(r.shape[0]):
            assert class_first_alarm(r[i], cd) == vectorised_first_alarm(t_up, t_dn, i), (cd, i)


def _null_rows(r: np.ndarray, label: dict) -> list[dict]:
    rows = []
    for cd in C_DOWNS:
        t_up, t_dn = first_alarms(r, cd)
        rows.append(dict(**label, c_down=cd, scenario="no change", paths=r.shape[0],
                         false_alarm_1day=float(np.mean((t_up >= 0) | (t_dn >= 0))),
                         false_alarm_down_1day=float(np.mean(t_dn >= 0)),
                         zero_share_before=float(np.mean(r[:, :CAL] == 0)),
                         zero_share_after=float(np.mean(r[:, CAL:] == 0)), **_scale_cols(r)))
    return rows


def _scale_cols(r: np.ndarray) -> dict:
    """Calibration scale (scaled MAD) in units of the pre-change sd: with rounded returns it is a multiple
    of the tick times 1.4826, so the clip c * scale in units of the true sd depends on the tick size."""
    s = np.median(np.abs(r[:, :CAL]), axis=1) * MAD_NORMAL_SCALE
    return dict(cal_scale_median=float(np.median(s)), cal_scale_mean=float(np.mean(s)))


def _change_rows(r: np.ndarray, ratio: float, label: dict) -> list[dict]:
    rows, delays = [], {}
    horizon = r.shape[1] - CAL
    for cd in C_DOWNS:
        t_up, t_dn = first_alarms(r, cd)
        det = t_dn >= 0
        delays[cd] = np.where(det, t_dn + 1, horizon)
        d = delays[cd]
        diff = d - delays[C_DOWNS[0]]  # paired difference against c_down = 3 (same paths)
        rows.append(dict(**label, c_down=cd, scenario=f"sd x{ratio}", paths=r.shape[0],
                         detected=float(det.mean()), median_delay=float(np.median(d)),
                         mean_delay=float(d.mean()), se_mean_delay=float(d.std(ddof=1) / np.sqrt(d.size)),
                         mean_diff_vs_c3=float(diff.mean()),
                         se_diff_vs_c3=float(diff.std(ddof=1) / np.sqrt(d.size)),
                         share_faster_than_c3=float(np.mean(diff < 0)),
                         share_slower_than_c3=float(np.mean(diff > 0)),
                         zero_share_before=float(np.mean(r[:, :CAL] == 0)),
                         zero_share_after=float(np.mean(r[:, CAL:] == 0)), **_scale_cols(r)))
    return rows


def part_a() -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    check_against_class(np.random.default_rng(SEED + 1))
    rows = []
    for law in LAWS:  # same random numbers, in the same order, as the first version of the script
        label = dict(law=law, tick=np.nan, tick_law="")
        rows += _null_rows(draw(law, (N_NULL, CAL + CAL), rng), label)
        for ratio in RATIOS:
            rows += _change_rows(draw(law, (N_CHANGE, CAL + HORIZON), rng) * np.r_[np.ones(CAL), np.full(HORIZON, ratio)],
                                 ratio, label)
        print(law, "done", flush=True)
    # Rounded prices: one set of normal draws, rounded to each tick size (common random numbers).
    z = rng.standard_normal((N_NULL, CAL + CAL))
    for tick in TICKS:
        rows += _null_rows(tick_paths(z, 1.0, tick), dict(law=f"ticks_{tick:g}", tick=tick, tick_law="fixed"))
    for ratio in RATIOS:
        z = rng.standard_normal((N_CHANGE_TICKS, CAL + HORIZON))
        for tick in TICKS:
            rows += _change_rows(tick_paths(z, ratio, tick), ratio,
                                 dict(law=f"ticks_{tick:g}", tick=tick, tick_law="fixed"))
        rows += _change_rows(tick_paths(z, ratio, OLD_TICK, "scaled"), ratio,
                             dict(law=f"ticks_{OLD_TICK:g}_scaled", tick=OLD_TICK, tick_law="scales with sd"))
        print("ticks, sd x", ratio, "done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "e11f_split_c_synthetic.csv", index=False)
    return out


def part_b() -> pd.DataFrame:
    rows = []
    spd = 1440
    for sym in D.CRYPTO_SYMBOLS:
        df = D.load_binance(sym, "2025-06-01", "2026-01-01", download=False)
        r_all = df["logret"]
        s_all = D.intraday_scale(r_all, lookback_days=20)
        keep = df.index >= pd.Timestamp("2025-07-01", tz=df.index.tz)
        u = (r_all / s_all).to_numpy()[keep]
        u[~np.isfinite(u)] = np.nan
        zero_share = float(np.mean(u[np.isfinite(u)] == 0))
        for cd in C_DOWNS:
            mon = VolatilityMonitor(threshold=A, cal_len=CAL, null_band=0.25, c_down=cd)
            ev = mon.run(u)
            res = {"up": [], "down": []}
            for e in ev:
                seg = u[e.index + 1: e.index + 1 + spd]
                seg = seg[np.isfinite(seg)]
                if seg.size <= spd // 2:  # no following day to compare with (end of the sample)
                    continue
                post = np.median(np.abs(seg)) * MAD_NORMAL_SCALE / e.baseline_scale
                res[e.direction].append(post)
            up, dn = np.array(res["up"]), np.array(res["down"])
            n_cal = sum(s["cal_len"] for s in mon.segments)
            rows.append(dict(symbol=sym, c_down=cd, zero_share=zero_share, days=len(u) / spd,
                             alarms_up=int(sum(e.direction == "up" for e in ev)),
                             alarms_down=int(sum(e.direction == "down" for e in ev)),
                             up_followed_by_higher=float(np.mean(up > 1)) if up.size else np.nan,
                             down_followed_by_lower=float(np.mean(dn < 1)) if dn.size else np.nan,
                             median_post_ratio_down=float(np.median(dn)) if dn.size else np.nan,
                             share_time_calibrating=n_cal / len(u),
                             down_side_inactive_segments=int(sum(not s["down_active"] for s in mon.segments)),
                             up_checked=int(up.size), down_checked=int(dn.size)))
        print(sym, "done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "e11f_split_c_real.csv", index=False)
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["a", "b"]
    if "a" in parts:
        print(part_a().to_string(index=False))
    if "b" in parts:
        print(part_b().to_string(index=False))
