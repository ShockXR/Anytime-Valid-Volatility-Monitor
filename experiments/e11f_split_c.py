"""E11f: a separate, tighter truncation constant for the decrease side.

docs/theory_delay.md (section 6.3) shows that truncation caps the evidence a decrease
detector can gather per step at log(B / (B - 1)) with B = tau^2 / m0, about 1 / c^2, and
predicts that c of about 1.5 to 2 on the decrease side (keeping c = 3 for increases) cuts
decrease delays by a quarter to a half. Each side is a valid e-detector for its own
truncated statistic, so the two-sided guarantee is unchanged. This script checks the
prediction and the price on data with ties and zeros, which a small c makes more visible.

Part a [SYNTHETIC]: false alarm rate within one day (no change) and median decrease delay
(sd x 0.8, 0.67, 0.5 at the first monitored step) for c_down = 3, 2 and 1.5, with the
plug-in baseline from a 1440-step calibration window, A = 28,800 (each side at 2A),
four return laws: Gaussian, Student t4, Gaussian with jumps, and rounded ("tick") returns.
Part b [REAL]: the five Binance coins, July to December 2025, band 0.25: alarm counts and
the share of alarms followed by a next day on the flagged side of the baseline.

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
LAWS = ("gaussian", "t4", "jumps", "ticks")


def draw(law: str, size, rng: np.random.Generator) -> np.ndarray:
    if law == "gaussian":
        return rng.standard_normal(size)
    if law == "t4":
        return rng.standard_t(4, size) / np.sqrt(2.0)  # unit variance
    if law == "jumps":
        z = rng.standard_normal(size)
        return z + (rng.random(size) < 0.005) * rng.normal(0.0, 8.0, size)
    if law == "ticks":
        return np.round(rng.standard_normal(size) / 0.5) * 0.5  # about 20% zero returns (|z| < 0.25)
    raise ValueError(law)


def first_alarms(r: np.ndarray, c_down: float, band: float = 0.0):
    """Vectorised first alarm of each side for paths in the rows of r (calibration = first CAL steps).

    Same arithmetic as VolatilityMonitor with the plug-in baseline up to the first alarm (checked
    against the class in check_against_class)."""
    cal, mon = r[:, :CAL], r[:, CAL:]
    scale = np.median(np.abs(cal), axis=1) * MAD_NORMAL_SCALE
    tau2 = (C_UP * scale) ** 2
    tau2_dn = (c_down * scale) ** 2
    m_up = np.mean(np.minimum(cal**2, tau2[:, None]), axis=1) * (1 + band)
    m_dn = np.mean(np.minimum(cal**2, tau2_dn[:, None]), axis=1) * (1 - band)
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
        hit_up = (S_up.mean(axis=1) >= 2 * A) & (t_up < 0) & (t_dn < 0)
        hit_dn = (S_dn.mean(axis=1) >= 2 * A) & (t_dn < 0) & (t_up < 0) & ~hit_up
        t_up[hit_up] = t
        t_dn[hit_dn] = t
    return t_up, t_dn


def check_against_class(rng: np.random.Generator) -> None:
    r = draw("jumps", (4, CAL + 3000), rng) * np.r_[np.ones(CAL), 0.6 * np.ones(3000)][None, :]
    for cd in (3.0, 1.5):
        t_up, t_dn = first_alarms(r, cd)
        for i in range(r.shape[0]):
            ev = VolatilityMonitor(threshold=A, cal_len=CAL, c_down=cd).run(r[i])
            first = (ev[0].index - CAL, ev[0].direction) if ev else None
            mine = (int(t_up[i]), "up") if t_up[i] >= 0 else ((int(t_dn[i]), "down") if t_dn[i] >= 0 else None)
            assert first == mine, (cd, i, first, mine)


def part_a(n_null: int = 2000, n_change: int = 600) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    check_against_class(np.random.default_rng(SEED + 1))
    rows = []
    for law in LAWS:
        null = draw(law, (n_null, CAL + CAL), rng)
        for cd in C_DOWNS:
            t_up, t_dn = first_alarms(null, cd)
            rows.append(dict(law=law, c_down=cd, scenario="no change", paths=n_null,
                             false_alarm_1day=float(np.mean((t_up >= 0) | (t_dn >= 0))),
                             false_alarm_down_1day=float(np.mean(t_dn >= 0)),
                             detected=np.nan, median_delay=np.nan, mean_delay=np.nan))
        for ratio in (0.8, 0.67, 0.5):
            horizon = 3000
            base = draw(law, (n_change, CAL + horizon), rng)
            base[:, CAL:] *= ratio
            for cd in C_DOWNS:
                t_up, t_dn = first_alarms(base, cd)
                det = t_dn >= 0
                d = np.where(det, t_dn + 1, horizon)
                rows.append(dict(law=law, c_down=cd, scenario=f"sd x{ratio}", paths=n_change,
                                 false_alarm_1day=np.nan, false_alarm_down_1day=np.nan,
                                 detected=float(det.mean()), median_delay=float(np.median(d)),
                                 mean_delay=float(d.mean())))
        print(law, "done", flush=True)
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
                if seg.size <= spd // 2:
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
                             down_side_inactive_segments=int(sum(not s["down_active"] for s in mon.segments))))
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
