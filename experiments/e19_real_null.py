"""E19. The guarantee on real data: null streams built from real returns, the confidence-sequence baseline
as the primary real-data configuration, the forecast-relative null on real data, and the split clip with the
confidence-sequence baseline.

Answers referee comment M1 (docs/review/referee_report_paper.md) and minor comment m32. Every design choice
was fixed in docs/preregistration_e19.md, committed before any real-data part of this script was run.

Part 1 [REAL marginal law, constructed streams]. Null streams with the real marginal law. Deseasonalised
one-minute returns of the calm days of each coin (Binance, July to December 2025; calm = daily robust
volatility within [0.8, 1.25] times the median of the previous 20 days) are pooled. A stream is the start of a
uniform random permutation of the pool (exchangeable, primary), a draw with replacement from the pool (i.i.d.),
or a stationary bootstrap of the calm days in calendar order (mean block lengths 30, 390 and 1440 minutes).
Monitors: (i) the default detector with the known baseline of the pool (tau and m0 from the pool), (ii) the
plug-in baseline, (iii) the confidence-sequence baseline (docs/math.md, Theorem 6.4), (iv) the rank detector
of E14 restarted at once after each alarm (awake-B). Endpoints: the probability of a false alarm within the
first monitored day and within five monitored days, and false alarms per monitored day over 20 monitored days
with a restart after each alarm. The same for Nifty 50 sessions (375 minutes).
Part 2 [REAL]. E6's monitoring and validation (Q6d, Q6e; E6 code imported, not changed) with the
confidence-sequence baseline, next to E6's plug-in results, and the monitoring time it costs.
Part 3 [REAL]. E10's forecast-relative monitor on the five coins: cost and validated shares.
Part 4 [SYNTHETIC]. The Gaussian row of Table S2 (E11f, split clip c_down) with the confidence-sequence baseline.
Post hoc parts (added after seeing results, marked as such in the report): ``blind`` (segments in which the cs
increase side cannot alarm because tau^2 <= m_up), ``p1std`` (bootstrap of day-standardised calm days) and
``p2boot`` (the primary bootstrap of Part 2 recomputed from the saved events, with Bonferroni columns).

Usage: nice -n 10 python experiments/e19_real_null.py [pools] [p1] [p2] [p3] [p4] [blind] [p1std] [p2boot] [figs]
Environment: E19_SCALE (default 1) multiplies the numbers of streams (smoke runs only); E19_OUT redirects
outputs (default results/).
"""

from __future__ import annotations

import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from avmonitor.confseq import eb_confseq  # noqa: E402
from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID  # noqa: E402
from avmonitor.robust import MAD_NORMAL_SCALE  # noqa: E402
from avmonitor.running import family_factors  # noqa: E402

OUT = Path(os.environ.get("E19_OUT", ROOT / "results"))
FIG = OUT / "figures"
SCALE = float(os.environ.get("E19_SCALE", "1"))

# ----------------------------------------------------------------------------
# Settings fixed in docs/preregistration_e19.md
# ----------------------------------------------------------------------------

C = 3.0
DELTA = 0.01
SEED = 19
CALM_LO, CALM_HI = 0.8, 1.25  # calm day: robust vol / trailing median in [0.8, 1.25]
TRAIL = 20  # days (sessions) in the trailing median
BANDS = (0.0, 0.25)
STREAM_TYPES = ("perm", "iid", "sb30", "sb390", "sb1440")
BLOCK = {"sb30": 30, "sb390": 390, "sb1440": 1440}
MONITORS = ("known", "plugin", "cs", "rank")
EXTRA_CALS = 10  # permutation streams hold cal + m + 10 cal returns (more alarms: the stream is cut)
BOOT_EXTRA_CALS = 40  # bootstrap streams: cal + m + 40 cal returns at most
# dataset -> (day length, calibration length, A for known / plug-in / rank, A for cs, monitored days m)
DATASETS = {
    "BTC": (1440, 1440, 28800.0, 36000.0, 20),
    "ETH": (1440, 1440, 28800.0, 36000.0, 20),
    "XRP": (1440, 1440, 28800.0, 36000.0, 20),
    "SOL": (1440, 1440, 28800.0, 36000.0, 20),
    "DOGE": (1440, 1440, 28800.0, 36000.0, 20),
    "Nifty 50": (375, 750, 7500.0, 9375.0, 20),
}
N_STREAMS = 2000  # per dataset, stream type and monitor (both bands use the same streams)


def wilson(k: int, n: int, z: float = 1.959964):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, centre - half, centre + half


# ----------------------------------------------------------------------------
# Stream generators: indices into a pool, one row per stream, one column per step
# ----------------------------------------------------------------------------


class StreamGen:
    """Indices into a pool of size N for ``n`` independent streams, produced step by step.

    kind "perm": the first ``length`` entries of a uniform random permutation of the pool (exchangeable);
    kind "iid": independent uniform draws with replacement; kind "sb<l>": the stationary bootstrap of
    Politis and Romano (1994) on the pool in its original order, wrapped circularly, mean block length l
    (each step starts a new block at a uniform position with probability 1 / l, else moves to the next
    position). Index -1 means the stream is exhausted (only "perm" streams end, after ``length`` steps).
    The same seed gives the same streams, so every monitor sees the same streams.
    """

    def __init__(self, kind: str, N: int, n: int, length: int, seed):
        self.kind, self.N, self.n, self.length = kind, int(N), int(n), int(length)
        self.rng = np.random.default_rng(seed)
        self.t = 0
        if kind == "perm":
            L = min(self.length, self.N)
            self.prefix = np.empty((n, L), dtype=np.int32)
            for i in range(n):
                self.prefix[i] = self.rng.permutation(self.N)[:L]
            self.length = L
        elif kind.startswith("sb"):
            self.p_new = 1.0 / BLOCK[kind]
            self.pos = self.rng.integers(0, self.N, n)
        elif kind != "iid":
            raise ValueError(kind)

    def next(self) -> np.ndarray:
        t = self.t
        self.t += 1
        if t >= self.length:
            return np.full(self.n, -1, dtype=np.int64)
        if self.kind == "perm":
            return self.prefix[:, t].astype(np.int64)
        if self.kind == "iid":
            return self.rng.integers(0, self.N, self.n)
        if t > 0:
            new = self.rng.random(self.n) < self.p_new
            jump = self.rng.integers(0, self.N, self.n)
            self.pos = np.where(new, jump, (self.pos + 1) % self.N)
        return self.pos.copy()


class ArrayGen:
    """Streams given as an array of values (rows = streams); used by the tests and the synthetic checks."""

    def __init__(self, values: np.ndarray):
        self.v = np.asarray(values, float)
        self.n, self.length = self.v.shape
        self.t = 0

    def next(self) -> np.ndarray:
        t = self.t
        self.t += 1
        if t >= self.length:
            return np.full(self.n, -1, dtype=np.int64)
        return np.arange(self.n) * self.length + t  # index into v.ravel()


# ----------------------------------------------------------------------------
# Engine 1: the default detector on X = min(r^2, tau^2) with a known, plug-in or cs baseline
# ----------------------------------------------------------------------------


def _usable(*vals):
    return np.all([np.isfinite(v) & (v > 0) for v in vals], axis=0)


def x_engine(pool: np.ndarray, gen, mode: str, band: float, A: float, cal: int, m_target: int, **kw):
    """One configuration of ``x_engine_multi`` (see there)."""
    return x_engine_multi(pool, gen, [(mode, band, A)], cal, m_target, **kw)[0]


def x_engine_multi(pool: np.ndarray, gen, configs, cal: int, m_target: int, *, c: float = C, delta: float = DELTA,
                   known: tuple[float, float] | None = None, tau_frac: float = 0.25, cs_c: float = 0.75,
                   horizon_only: bool = False):
    """Vectorised ``VolatilityMonitor`` (c_down = c, no missing values) for several configurations at once.

    ``configs``: list of (mode, band, A). Every configuration runs on every stream of ``gen`` (rows are stacked
    configuration by configuration, so one pass over the streams serves all of them).
    mode "plugin" or "cs": exactly the class (baseline="plugin" or "cs", default settings), including the
    extension rule for degenerate calibration windows; after an alarm the row recalibrates on the next
    ``cal`` returns. mode "known": tau^2 and m0 given (``known``), the first ``cal`` returns of a stream are
    skipped (so that every monitor starts at the same step), and after an alarm the statistics restart at
    zero at once (Corollary C with delta = 0). Each row stops once it has had ``m_target`` monitored steps
    (``horizon_only``: after its first alarm, or after ``m_target`` monitored steps).

    Returns, per configuration, a dict of per-stream arrays: first alarm (monitored step count at the alarm,
    1-based; 0 if none) and its side (+1 up, -1 down), alarms within the first ``m_target`` monitored steps
    (up, down), monitored steps, stream steps used, calibration windows, and whether the stream ended early.
    """
    up_grid = np.asarray(DEFAULT_UP_GRID, float)
    fracs = np.asarray(DEFAULT_DOWN_FRACTIONS, float)
    ns = gen.n
    k = len(configs)
    n = ns * k
    src = np.tile(np.arange(ns), k)
    code = {"known": 0, "plugin": 1, "cs": 2}
    mode_r = np.repeat([code[m] for m, _, _ in configs], ns)
    band_r = np.repeat([float(b) for _, b, _ in configs], ns)
    A_r = np.repeat([float(a) for _, _, a in configs], ns)
    is_known, is_plugin, is_cs = mode_r == 0, mode_r == 1, mode_r == 2
    vals = pool if not isinstance(gen, ArrayGen) else gen.v.ravel()
    cal_max = 2 * cal
    n_tau = int(round(tau_frac * cal))
    n_cs = cal - n_tau
    mon = np.zeros(n, bool)
    done = np.zeros(n, bool)
    S_up = np.zeros((n, up_grid.size))
    S_dn = np.zeros((n, fracs.size))
    tau2 = np.ones(n)
    m_up = np.ones(n)
    m_dn = np.ones(n)
    lam = np.zeros((n, fracs.size))
    dn_on = np.zeros(n, bool)
    buf = np.zeros((n, cal_max + 1))
    wlen = np.zeros(n, dtype=np.int64)  # returns in the current calibration window
    ext = np.zeros(n, bool)  # window being extended
    tau_set = np.zeros(n, bool)  # cs: tau fixed
    n_tau_used = np.zeros(n, dtype=np.int64)
    first_t = np.zeros(n, dtype=np.int64)
    first_dir = np.zeros(n, dtype=np.int64)
    n_up = np.zeros(n, dtype=np.int64)
    n_dn = np.zeros(n, dtype=np.int64)
    mcount = np.zeros(n, dtype=np.int64)
    used = np.zeros(n, dtype=np.int64)
    ncal = np.zeros(n, dtype=np.int64)
    ended = np.zeros(n, bool)
    if is_known.any():
        k_tau2, k_m0 = known
    rows_all = np.arange(n)

    def set_baseline(rows, t2, m0u, m0d):
        """Start monitoring on ``rows`` against m_up = m0u (1 + b), m_dn = m0d (1 - b)."""
        mu = m0u * (1 + band_r[rows])
        md = m0d * (1 - band_r[rows])
        tau2[rows] = t2
        m_up[rows] = mu
        m_dn[rows] = md
        B = np.where(md > 0, t2 / np.where(md > 0, md, 1.0), np.inf)
        on = np.isfinite(B) & (B > 1.0)
        dn_on[rows] = on
        lam[rows] = np.where(on[:, None], fracs[None, :] / np.where(on, B - 1.0, 1.0)[:, None], 0.0)
        S_up[rows] = 0.0
        S_dn[rows] = 0.0
        mon[rows] = True
        ncal[rows] += 1
        wlen[rows] = 0
        ext[rows] = False
        tau_set[rows] = False

    t = -1
    while not done.all():
        t += 1
        idx = gen.next()[src]
        gone = (idx < 0) & ~done
        if gone.any():
            ended[gone] = True
            done[gone] = True
            if done.all():
                break
        live = ~done
        v = vals[np.where(idx >= 0, idx, 0)]
        used[live] += 1
        # ---- monitoring rows
        mrows = mon & live
        if mrows.any():
            x = np.minimum(v * v, tau2)
            fu = 1.0 + up_grid[None, :] * (x / m_up - 1.0)[:, None]
            S_up = np.where(mrows[:, None], fu * (S_up + 1.0), S_up)
            dmask = (mrows & dn_on)[:, None]
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                fd = 1.0 + lam * (1.0 - x / m_dn)[:, None]
                S_dn = np.where(dmask, fd * (S_dn + 1.0), S_dn)
            mcount[mrows] += 1
            mu = S_up.mean(axis=1)
            md = S_dn.mean(axis=1)
            a_up = mrows & (mu >= 2 * A_r)
            a_dn = mrows & ~a_up & (md >= 2 * A_r)
            alarm = a_up | a_dn
            if alarm.any():
                new_first = alarm & (first_t == 0)
                first_t[new_first] = mcount[new_first]
                first_dir[new_first] = np.where(a_up[new_first], 1, -1)
                inside = alarm & (mcount <= m_target)
                n_up[inside & a_up] += 1
                n_dn[inside & a_dn] += 1
                S_up[alarm] = 0.0
                S_dn[alarm] = 0.0
                recal = alarm & ~is_known
                mon[recal] = False
                wlen[recal] = 0
                ext[recal] = False
                tau_set[recal] = False
            fin = mrows & (mcount >= m_target)
            if horizon_only:
                fin = fin | (mrows & (first_t > 0))
            done[fin] = True
            crow = live & ~mrows & ~done
        else:
            crow = live
        # ---- calibrating rows (they did not monitor at this step)
        if not crow.any():
            continue
        rows = rows_all[crow]
        buf[rows, wlen[rows]] = v[rows]
        wlen[rows] += 1
        rk = rows[is_known[rows] & (wlen[rows] >= cal)]
        if rk.size:
            set_baseline(rk, np.full(rk.size, k_tau2), np.full(rk.size, k_m0), np.full(rk.size, k_m0))
        for i in rows[is_plugin[rows] & (wlen[rows] >= cal)]:
            r = buf[i, : wlen[i]]
            scale = float(np.median(np.abs(r)) * MAD_NORMAL_SCALE)
            t2 = (c * scale) ** 2
            m0 = float(np.mean(np.minimum(r * r, t2)))
            ok = bool(_usable(scale, t2, m0))
            if not (ok and (not ext[i] or r.size >= cal // 2)):
                if wlen[i] < cal_max:
                    ext[i] = True
                    continue
                if not ok:
                    nz = np.abs(r[r != 0])
                    scale = float(np.median(nz) * MAD_NORMAL_SCALE) if nz.size else 0.0
                    t2 = (c * scale) ** 2
                    m0 = float(np.mean(np.minimum(r * r, t2)))
                    if not bool(_usable(scale, t2, m0)):
                        wlen[i] = 0
                        ext[i] = False
                        continue
            set_baseline(np.array([i]), np.array([t2]), np.array([m0]), np.array([m0]))
        # cs: tau from the first part, then a confidence sequence on the next n_cs returns. Only rows at a
        # decision point (tau part complete, or sequence part complete) need work at this step.
        rc = rows[is_cs[rows]]
        due = rc[(~tau_set[rc] & (wlen[rc] >= n_tau)) | (tau_set[rc] & (wlen[rc] - n_tau_used[rc] >= n_cs))]
        for i in due:
            if not tau_set[i]:
                r = buf[i, : wlen[i]]
                scale = float(np.median(np.abs(r)) * MAD_NORMAL_SCALE)
                t2 = (c * scale) ** 2
                ok = bool(_usable(scale, t2))
                if not (ok and r.size >= max(1, n_tau // 2)):
                    if wlen[i] < cal_max - n_cs:
                        ext[i] = True
                        continue
                    if not ok:
                        nz = np.abs(r[r != 0])
                        scale = float(np.median(nz) * MAD_NORMAL_SCALE) if nz.size else 0.0
                        t2 = (c * scale) ** 2
                        if not bool(_usable(scale, t2)):
                            wlen[i] = 0
                            ext[i] = False
                            continue
                tau_set[i] = True
                tau2[i] = t2
                n_tau_used[i] = wlen[i]
                continue
            x = np.minimum(buf[i, n_tau_used[i] : wlen[i]] ** 2, tau2[i])
            cs = eb_confseq(x, tau2[i], delta / 2, delta / 2, n_plan=n_cs, c=cs_c)
            set_baseline(np.array([i]), np.array([tau2[i]]), np.array([float(cs.upper[0, -1])]),
                         np.array([float(cs.lower[0, -1])]))
    out = []
    for j in range(k):
        sl = slice(j * ns, (j + 1) * ns)
        out.append(dict(first_t=first_t[sl], first_dir=first_dir[sl], n_up=n_up[sl], n_dn=n_dn[sl],
                        mon_steps=mcount[sl], used=used[sl], calibrations=ncal[sl], ended=ended[sl]))
    return out


# ----------------------------------------------------------------------------
# Engine 2: the rank (conformal) detector of E14, route B, restarted at once after each alarm (awake-B)
# ----------------------------------------------------------------------------


class _Fenwick:
    """One Fenwick tree per row (counts of reference values with rank <= q), int32 storage."""

    def __init__(self, n: int, size: int):
        self.size = size
        self.tree = np.zeros((n, size + 1), dtype=np.int32)
        self.rows = np.arange(n)

    def add(self, q: np.ndarray, rows: np.ndarray):
        q = q[rows].copy()
        while rows.size:
            self.tree[rows, q] += 1
            q = q + (q & -q)
            keep = q <= self.size
            rows, q = rows[keep], q[keep]

    def count_le(self, q: np.ndarray, rows: np.ndarray) -> np.ndarray:
        s = np.zeros(self.rows.size, dtype=np.int64)
        q = q[rows].copy()
        while rows.size:
            s[rows] += self.tree[rows, q]
            q = q - (q & -q)
            keep = q > 0
            rows, q = rows[keep], q[keep]
        return s


def rank_engine(ranks: np.ndarray, gen, band: float, A: float, cal: int, m_target: int, *, theta_rng=None,
                theta: np.ndarray | None = None, family: str = "lr", horizon_only: bool = False):
    """Vectorised ``ConformalMonitor(cal_len=cal, restart_len=0, family=family, band=band)`` (awake-B, E14).

    ``ranks``: dense rank (1-based) of |value| for every pool entry, so that comparisons of |r| are
    comparisons of ranks. The first ``cal`` returns of a stream form the reference set; from then on every
    return is monitored: smoothed conformal p-value among the reference set, Shiryaev-Roberts update of the
    increase and decrease families (each at 2A), then the return joins the reference set. After an alarm the
    reference set is emptied and betting continues with the next return. Tie-breakers: ``theta`` (rows x
    steps, for the tests) or ``theta_rng`` (one uniform per row and step). Returns a dict as ``x_engine``.
    """
    fam = family_factors(family)
    n = gen.n
    fw = _Fenwick(n, int(ranks.max()) + 1)
    ref = np.zeros(n, dtype=np.int64)
    S_up = np.zeros((n, 3))
    S_dn = np.zeros((n, 3))
    R = dict(first_t=np.zeros(n, np.int64), first_dir=np.zeros(n, np.int64), n_up=np.zeros(n, np.int64),
             n_dn=np.zeros(n, np.int64), mon_steps=np.zeros(n, np.int64), used=np.zeros(n, np.int64),
             calibrations=np.ones(n, np.int64), ended=np.zeros(n, bool))
    done = np.zeros(n, bool)
    rows_all = np.arange(n)
    rk_vals = ranks if not isinstance(gen, ArrayGen) else ranks.ravel()
    t = -1
    while not done.all():
        t += 1
        idx = gen.next()
        gone = (idx < 0) & ~done
        if gone.any():
            R["ended"][gone] = True
            done |= gone
            if done.all():
                break
        q = rk_vals[np.where(idx >= 0, idx, 0)].astype(np.int64)
        live = ~done
        live_rows = rows_all[live]
        R["used"][live] += 1
        if t >= cal:
            le = fw.count_le(q, live_rows)
            lt = fw.count_le(q - 1, live_rows)
            th = theta[:, t] if theta is not None else theta_rng.random(n)
            p = (ref - le + th * (le - lt + 1)) / (ref + 1.0)
            fu, fd = fam(p)
            if band:
                fu = 1.0 + (fu - 1.0 - band) / (1.0 + band)
                fd = 1.0 + (fd - 1.0 - band) / (1.0 + band)
            with np.errstate(over="ignore"):
                S_up = np.where(live[:, None], fu * (S_up + 1.0), S_up)
                S_dn = np.where(live[:, None], fd * (S_dn + 1.0), S_dn)
            R["mon_steps"][live] += 1
            mu, md = S_up.mean(axis=1), S_dn.mean(axis=1)
            a_up = live & (mu >= 2 * A)
            a_dn = live & ~a_up & (md >= 2 * A)
            alarm = a_up | a_dn
            fw.add(q, live_rows)
            ref[live_rows] += 1
            if alarm.any():
                nf = alarm & (R["first_t"] == 0)
                R["first_t"][nf] = R["mon_steps"][nf]
                R["first_dir"][nf] = np.where(a_up[nf], 1, -1)
                inside = alarm & (R["mon_steps"] <= m_target)
                R["n_up"][inside & a_up] += 1
                R["n_dn"][inside & a_dn] += 1
                S_up[alarm] = 0.0
                S_dn[alarm] = 0.0
                R["calibrations"][alarm] += 1
                rr = rows_all[alarm]
                fw.tree[rr] = 0
                ref[rr] = 0
            fin = live & (R["mon_steps"] >= m_target)
            if horizon_only:
                fin |= live & (R["first_t"] > 0)
            done |= fin
        else:
            fw.add(q, live_rows)
            ref[live_rows] += 1
    return R


# ----------------------------------------------------------------------------
# Data: deseasonalised returns, calm days and pools
# ----------------------------------------------------------------------------

CRYPTO = ("BTC", "ETH", "XRP", "SOL", "DOGE")
NIFTY_LOAD = ("2021-10-01", "2026-11-01")
NIFTY_PERIOD = ("2022-01-01", "2026-11-01")
MIN_FINITE = {"crypto": 1000, "nifty": 300}
SCRATCH = Path(os.environ.get("E19_SCRATCH", ROOT / "data" / "processed" / "e19_pools"))  # derived data, git-ignored


def robust_vol(x: np.ndarray, min_finite: int) -> float:
    x = x[np.isfinite(x)]
    return float(np.median(np.abs(x)) * MAD_NORMAL_SCALE) if x.size >= min_finite else np.nan


def calm_flags(day_keys: np.ndarray, u: np.ndarray, in_period: np.ndarray, min_finite: int):
    """Calm days: robust vol of the day within [CALM_LO, CALM_HI] times the median of the previous TRAIL days.

    ``day_keys``: one key per return (UTC day or session, nondecreasing); ``in_period``: per return, whether its
    day belongs to the period. Returns a table with one row per day of the period (key, v, trailing median,
    ratio, calm) and the boolean mask of the returns that belong to calm days.
    """
    keys, start = np.unique(day_keys, return_index=True)
    order = np.argsort(start)
    keys, start = keys[order], start[order]
    end = np.r_[start[1:], len(u)]
    v = np.array([robust_vol(u[a:b], min_finite) for a, b in zip(start, end)])
    rows, mask = [], np.zeros(len(u), bool)
    for k in range(len(keys)):
        if not in_period[start[k]]:
            continue
        prev = v[:k][np.isfinite(v[:k])][-TRAIL:]
        R = float(np.median(prev)) if prev.size == TRAIL else np.nan
        ratio = v[k] / R if np.isfinite(v[k]) and np.isfinite(R) and R > 0 else np.nan
        calm = bool(np.isfinite(ratio) and CALM_LO <= ratio <= CALM_HI)
        rows.append(dict(day=keys[k], robust_vol=v[k], trailing_median=R, ratio=ratio, calm=calm,
                         finite_returns=int(np.isfinite(u[start[k]:end[k]]).sum())))
        if calm:
            mask[start[k]:end[k]] = True
    return pd.DataFrame(rows), mask


def load_crypto(coin: str):
    import e06b_five_coins_monitor as E6

    idx, u, i0 = E6.prepare(coin, "Jul-Dec 2025")
    keys = idx.floor("D").to_numpy()
    in_period = np.zeros(len(u), bool)
    in_period[i0:] = True
    return keys, u, in_period


def load_nifty_u():
    import e07_nifty as E7

    df, _, _ = E7.load_clean(*NIFTY_LOAD)
    u = E7.deseasonalise(df).to_numpy(float)
    sess = pd.DatetimeIndex(df["session"]).tz_localize(None).to_numpy()  # IST trading date at midnight, naive
    lo, hi = (np.datetime64(x, "D") for x in NIFTY_PERIOD)
    sd = sess.astype("datetime64[D]")
    return sess, u, (sd >= lo) & (sd < hi)


def build_pool(name: str):
    """Pool of the calm days of one data set: values, dense ranks of |value|, known baseline, descriptives."""
    if name == "Nifty 50":
        keys, u, in_period = load_nifty_u()
        mf = MIN_FINITE["nifty"]
    else:
        keys, u, in_period = load_crypto(name)
        mf = MIN_FINITE["crypto"]
    days, mask = calm_flags(keys, u, in_period, mf)
    pool = u[mask]
    pool = pool[np.isfinite(pool)]
    scale = float(np.median(np.abs(pool)) * MAD_NORMAL_SCALE)
    tau2 = (C * scale) ** 2
    X = np.minimum(pool**2, tau2)
    m0 = float(X.mean())
    _, inv = np.unique(np.abs(pool), return_inverse=True)
    info = dict(dataset=name, period_days=int(len(days)), calm_days=int(days.calm.sum()),
                days_undefined=int(days.ratio.isna().sum()), pool_returns=int(pool.size),
                zero_share=float(np.mean(pool == 0)), scale=scale, tau2=tau2, m0=m0, B=tau2 / m0,
                clipped_share=float(np.mean(pool**2 > tau2)), sd_X_over_m0=float(X.std() / m0),
                distinct_abs_values=int(inv.max() + 1))
    days.insert(0, "dataset", name)
    return pool, (inv + 1).astype(np.int32), info, days


def pool_path(name: str) -> Path:
    return SCRATCH / f"pool_{name.replace(' ', '_')}.npz"


def prepare_pools(names) -> pd.DataFrame:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    infos, days_all = [], []
    for name in names:
        pool, ranks, info, days = build_pool(name)
        np.savez(pool_path(name), pool=pool, ranks=ranks, tau2=info["tau2"], m0=info["m0"])
        infos.append(info)
        days_all.append(days)
        print(f"pool {name}: {info['calm_days']} calm of {info['period_days']} days, {info['pool_returns']} returns",
              flush=True)
    out = pd.DataFrame(infos)
    out.to_csv(OUT / "e19_pools.csv", index=False)
    pd.concat(days_all).to_csv(OUT / "e19_calm_days.csv", index=False)
    return out


# ----------------------------------------------------------------------------
# Part 1: null streams
# ----------------------------------------------------------------------------

BATCH = 1000


def stream_plan(names):
    """Jobs (data set, stream type, batch index, rows); perm 4000 streams, the other types 2000."""
    jobs = []
    for d, name in enumerate(names):
        for s, kind in enumerate(STREAM_TYPES):
            total = int(round((4000 if kind == "perm" else 2000) * SCALE))
            nb = max(1, int(math.ceil(total / BATCH)))
            for b in range(nb):
                rows = min(BATCH, total - b * BATCH)
                jobs.append((name, d, kind, s, b, rows))
    return jobs


def _drift(pool, tau2, m0, prefix, horizon):
    """Largest |E[X_t | past] / m0 - 1| over the first ``horizon`` steps of permutation streams (per row)."""
    X = np.minimum(pool**2, tau2)
    N, S = X.size, X.sum()
    h = min(horizon, prefix.shape[1])
    cs = np.cumsum(X[prefix[:, :h]], axis=1)
    drawn = np.concatenate([np.zeros((prefix.shape[0], 1)), cs[:, :-1]], axis=1)
    rem = (S - drawn) / (N - np.arange(h))[None, :]
    return np.max(np.abs(rem / m0 - 1.0), axis=1)


def p1_job(job):
    name, d, kind, s, b, rows = job
    ck = SCRATCH / "p1_jobs" / f"{name.replace(' ', '_')}_{kind}_{b}.pkl"
    if ck.exists():
        return pd.read_pickle(ck)
    day, cal, A, A_cs, mdays = DATASETS[name]
    m = mdays * day
    z = np.load(pool_path(name))
    pool, ranks, tau2, m0 = z["pool"], z["ranks"], float(z["tau2"]), float(z["m0"])
    N = pool.size
    length = cal + m + (EXTRA_CALS if kind == "perm" else BOOT_EXTRA_CALS) * cal
    seq = np.random.SeedSequence([SEED, d, s, b])
    seed_streams, seed_theta = seq.spawn(2)
    out = []
    t0 = time.time()
    drift = None
    configs = [(mon, band, A_cs if mon == "cs" else A) for mon in MONITORS if mon != "rank" for band in BANDS]
    gen = StreamGen(kind, N, rows, length, np.random.default_rng(seed_streams))
    if kind == "perm":
        drift = _drift(pool, tau2, m0, gen.prefix, cal + m)
    res = x_engine_multi(pool, gen, configs, cal, m, known=(tau2, m0))
    for (mon, band, _), R in zip(configs, res):
        out.append(dict(dataset=name, stream=kind, monitor=mon, band=band, batch=b, **R))
    t1 = time.time()
    for band in BANDS:
        gen = StreamGen(kind, N, rows, length, np.random.default_rng(seed_streams))
        R = rank_engine(ranks, gen, band, A, cal, m, theta_rng=np.random.default_rng(seed_theta))
        out.append(dict(dataset=name, stream=kind, monitor="rank", band=band, batch=b, **R))
    print(f"p1 {name} {kind} batch {b}: {rows} streams, X {t1 - t0:.0f}s, rank {time.time() - t1:.0f}s", flush=True)
    res = (out, (name, kind, b, drift))
    ck.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(res, ck)
    return res


def bounds(name: str, mon: str):
    day, cal, A, A_cs, mdays = DATASETS[name]
    m = mdays * day
    if mon == "cs":
        return DELTA + day / A_cs, DELTA + 5 * day / A_cs, (DELTA + m / A_cs) / (1 - DELTA) / mdays
    return day / A, 5 * day / A, day / A


def summarise_p1(results, drifts) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One row per data set, stream type, monitor and band (and pooled over the five coins)."""
    cat = {}
    for r in results:
        key = (r["dataset"], r["stream"], r["monitor"], r["band"])
        cat.setdefault(key, []).append(r)
    rows = []

    def row_for(label, parts, name_for_bounds):
        day = DATASETS[name_for_bounds][0]
        first = np.concatenate([p["first_t"] for p in parts])
        fdir = np.concatenate([p["first_dir"] for p in parts])
        alarms = np.concatenate([p["n_up"] + p["n_dn"] for p in parts])
        n_up = np.concatenate([p["n_up"] for p in parts])
        mon = np.concatenate([p["mon_steps"] for p in parts])
        used = np.concatenate([p["used"] for p in parts])
        ended = np.concatenate([p["ended"] for p in parts])
        n = first.size
        k1 = int(np.sum((first > 0) & (first <= day)))
        k5 = int(np.sum((first > 0) & (first <= 5 * day)))
        b1, b5, bd = bounds(name_for_bounds, label[2])
        p1, l1, h1 = wilson(k1, n)
        p5, l5, h5 = wilson(k5, n)
        days_mon = mon / day
        rate = alarms.sum() / days_mon.sum()
        # ratio estimator: linearised standard error over independent streams
        se = math.sqrt(np.sum((alarms - rate * days_mon) ** 2)) / days_mon.sum()
        rate_cal = alarms.sum() / (used / day).sum()
        return dict(dataset=label[0], stream=label[1], monitor=label[2], band=label[3], streams=n,
                    one_day=p1, one_day_lo=l1, one_day_hi=h1, one_day_bound=b1, one_day_alarms=k1,
                    five_days=p5, five_days_lo=l5, five_days_hi=h5, five_days_bound=b5, five_days_alarms=k5,
                    per_monitored_day=rate, per_monitored_day_lo=rate - 1.959964 * se,
                    per_monitored_day_hi=rate + 1.959964 * se, per_monitored_day_bound=bd,
                    per_calendar_day=rate_cal, monitored_share=float(mon.sum() / used.sum()),
                    first_alarm_up_share=float(np.mean(fdir[first > 0] == 1)) if np.any(first > 0) else np.nan,
                    alarms_up_share=float(n_up.sum() / alarms.sum()) if alarms.sum() else np.nan,
                    streams_ended_early=int(ended.sum()), mean_monitored_days=float(days_mon.mean()))

    for key, parts in sorted(cat.items(), key=lambda kv: (list(DATASETS).index(kv[0][0]), STREAM_TYPES.index(kv[0][1]),
                                                           MONITORS.index(kv[0][2]), kv[0][3])):
        rows.append(row_for(key, parts, key[0]))
    for kind in STREAM_TYPES:
        for mon in MONITORS:
            for band in BANDS:
                parts = [p for c in CRYPTO for p in cat.get((c, kind, mon, band), [])]
                if parts:
                    rows.append(row_for(("five coins", kind, mon, band), parts, "BTC"))
    df = pd.DataFrame(rows)
    for col in ("one_day", "five_days", "per_monitored_day"):
        df[f"{col}_vs_bound"] = np.where(df[f"{col}_lo"] > df[f"{col}_bound"], "above",
                                         "at or below")
    dr = []
    for name in DATASETS:
        vals = [x for (nm, kind, b, x) in drifts if nm == name and kind == "perm" and x is not None]
        if vals:
            v = np.concatenate(vals)
            dr.append(dict(dataset=name, streams=v.size, drift_median=float(np.median(v)),
                           drift_q99=float(np.quantile(v, 0.99)), drift_max=float(v.max())))
    return df, pd.DataFrame(dr)


# ----------------------------------------------------------------------------
# Post hoc (added after seeing Part 1): day-standardised pools
# ----------------------------------------------------------------------------
# The calm days of a pool differ in level (their robust volatility varies by a factor of about 1.5 to 2.1 between
# the 10% and 90% quantiles), so a bootstrap block from one day next to a block from another is a change of level.
# Dividing each calm day's returns by that day's robust volatility removes the between-day level differences and
# keeps the intraday dependence, which separates the two sources of alarms on the bootstrap streams.

STD_TYPES = ("perm", "sb30", "sb390", "sb1440")
STD_ROWS = 500


def build_std_pool(name: str):
    if name == "Nifty 50":
        keys, u, in_period = load_nifty_u()
        mf = MIN_FINITE["nifty"]
    else:
        keys, u, in_period = load_crypto(name)
        mf = MIN_FINITE["crypto"]
    days, mask = calm_flags(keys, u, in_period, mf)
    k, x = keys[mask], u[mask]
    fin = np.isfinite(x)
    k, x = k[fin], x[fin]
    calm = days[days.calm]
    v = pd.Series(calm.robust_vol.to_numpy(), index=calm.day.to_numpy()).reindex(k).to_numpy()
    pool = x / v
    scale = float(np.median(np.abs(pool)) * MAD_NORMAL_SCALE)
    tau2 = (C * scale) ** 2
    X = np.minimum(pool**2, tau2)
    m0 = float(X.mean())
    _, inv = np.unique(np.abs(pool), return_inverse=True)
    info = dict(dataset=name, pool="day-standardised", pool_returns=int(pool.size), zero_share=float(np.mean(pool == 0)),
                scale=scale, tau2=tau2, m0=m0, B=tau2 / m0, clipped_share=float(np.mean(pool**2 > tau2)),
                sd_X_over_m0=float(X.std() / m0))
    return pool, (inv + 1).astype(np.int32), info


def part_std(names=tuple(DATASETS), workers: int = 2):
    infos, jobs = [], []
    for d, name in enumerate(names):
        pool, ranks, info = build_std_pool(name)
        sname = f"{name} std"
        np.savez(pool_path(sname), pool=pool, ranks=ranks, tau2=info["tau2"], m0=info["m0"])
        infos.append(info)
        DATASETS[sname] = DATASETS[name]  # same day length, calibration and A; the workers are forked after this
        for kind in STD_TYPES:
            jobs.append((sname, 100 + d, kind, STREAM_TYPES.index(kind), 0, STD_ROWS))
    pd.DataFrame(infos).to_csv(OUT / "e19_pools_std.csv", index=False)
    with Pool(workers) as pool:
        res = pool.map(p1_job, jobs, chunksize=1)
    strip = lambda s: s[:-4] if s.endswith(" std") else s
    results = [dict(r, dataset=strip(r["dataset"])) for out, _ in res for r in out]
    drifts = [(strip(nm), kind, b, x) for _, (nm, kind, b, x) in res]
    df, drift = summarise_p1(results, drifts)
    df.insert(1, "pool", "day-standardised")
    df.to_csv(OUT / "e19_null_std.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(df[["dataset", "stream", "monitor", "band", "streams", "one_day", "five_days", "per_monitored_day",
                  "streams_ended_early"]].round(4).to_string(index=False))
        print(drift.to_string(index=False))
    return df


def p1_pooled() -> pd.DataFrame:
    """The pooled five-coin rows of e19_null_streams.csv, in their own file (as pre-registered)."""
    df = pd.read_csv(OUT / "e19_null_streams.csv")
    out = df[df.dataset == "five coins"]
    out.to_csv(OUT / "e19_null_pooled.csv", index=False)
    return out


def part1(names=tuple(DATASETS), workers: int = 2):
    pools = prepare_pools(names)
    print(pools.to_string(index=False), flush=True)
    jobs = stream_plan(names)
    print(f"part 1: {len(jobs)} jobs", flush=True)
    with Pool(workers) as pool:
        res = pool.map(p1_job, jobs, chunksize=1)
    results = [r for out, _ in res for r in out]
    drifts = [dr for _, dr in res]
    df, drift = summarise_p1(results, drifts)
    df.to_csv(OUT / "e19_null_streams.csv", index=False)
    drift.to_csv(OUT / "e19_perm_drift.csv", index=False)
    p1_pooled()
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(df[["dataset", "stream", "monitor", "band", "streams", "one_day", "one_day_hi", "five_days",
                  "five_days_hi", "per_monitored_day", "per_monitored_day_hi", "streams_ended_early"]].round(4).to_string(index=False))
        print(drift.to_string(index=False))
    return df


# ----------------------------------------------------------------------------
# Part 4 [SYNTHETIC]: Table S2, Gaussian row, with the confidence-sequence baseline (m32)
# ----------------------------------------------------------------------------


def first_alarms_cs(r: np.ndarray, c_down: float, band: float = 0.0, cal: int = 1440, threshold: float = 36000.0,
                    c_up: float = C, delta: float = DELTA, tau_frac: float = 0.25, cs_c: float = 0.75):
    """Vectorised first alarm of each side, ``VolatilityMonitor(baseline="cs", c_down=c_down)`` up to its first alarm.

    tau from the first round(tau_frac * cal) returns; the upper bound from a confidence sequence on min(r^2, tau^2)
    and the lower bound from one on min(r^2, tau_dn^2), each at delta / 2 per side (one sequence when
    c_down = c). No missing values and no degenerate windows (checked). Returns monitored indices (0-based) of
    the first up and the first down alarm (-1 if none; only the side that alarms first is recorded)."""
    n_tau = int(round(tau_frac * cal))
    cal_r, mon = r[:, :cal], r[:, cal:]
    scale = np.median(np.abs(cal_r[:, :n_tau]), axis=1) * MAD_NORMAL_SCALE
    if np.any(scale <= 0):
        raise ValueError("degenerate calibration window")
    tau2 = (c_up * scale) ** 2
    tau2_dn = (c_down * scale) ** 2
    half = delta / 2
    rest = cal_r[:, n_tau:] ** 2
    cs_u = eb_confseq(np.minimum(rest, tau2[:, None]), tau2, half, half, n_plan=cal - n_tau, c=cs_c)
    if c_down == c_up:
        cs_d = cs_u
    else:
        cs_d = eb_confseq(np.minimum(rest, tau2_dn[:, None]), tau2_dn, half, half, n_plan=cal - n_tau, c=cs_c)
    m_up = cs_u.upper[:, -1] * (1 + band)
    m_dn = cs_d.lower[:, -1] * (1 - band)
    on = m_dn > 0
    B = np.where(on, tau2_dn / np.where(on, m_dn, 1.0), np.inf)
    on &= B > 1
    lam_dn = np.where(on[:, None], np.asarray(DEFAULT_DOWN_FRACTIONS)[None, :] / np.where(on, B - 1.0, 1.0)[:, None], 0.0)
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
        with np.errstate(divide="ignore", invalid="ignore"):
            S_dn = np.where(on[:, None], (1.0 + lam_dn * (1.0 - xd / m_dn)[:, None]) * (S_dn + 1.0), S_dn)
        hit_up = (S_up.mean(axis=1) >= 2 * threshold) & (t_up < 0) & (t_dn < 0)
        hit_dn = (S_dn.mean(axis=1) >= 2 * threshold) & (t_dn < 0) & (t_up < 0) & ~hit_up
        t_up[hit_up] = t
        t_dn[hit_dn] = t
        if t % 64 == 0 and np.all((t_up >= 0) | (t_dn >= 0)):
            break
    return t_up, t_dn


def part4() -> pd.DataFrame:
    import e11f_split_c as E11F

    rng = np.random.default_rng(E11F.SEED)  # E11f's draw order: Gaussian null paths, then the three decrease cells
    CAL, HOR = E11F.CAL, E11F.HORIZON
    null = E11F.draw("gaussian", (E11F.N_NULL, CAL + CAL), rng)
    cells = [(ratio, E11F.draw("gaussian", (E11F.N_CHANGE, CAL + HOR), rng) * np.r_[np.ones(CAL), np.full(HOR, ratio)])
             for ratio in E11F.RATIOS]
    rows = []
    for baseline in ("plugin", "cs"):
        for cd in E11F.C_DOWNS:
            if baseline == "plugin":
                tu, td = E11F.first_alarms(null, cd)
            else:
                tu, td = first_alarms_cs(null, cd)
            k = int(np.sum((tu >= 0) | (td >= 0)))
            p, lo, hi = wilson(k, null.shape[0])
            rows.append(dict(baseline=baseline, A=E11F.A if baseline == "plugin" else 36000.0, c_down=cd,
                             scenario="no change", paths=null.shape[0], false_alarm_1day=p, wilson_lo=lo, wilson_hi=hi,
                             false_alarm_down_1day=float(np.mean(td >= 0)), bound=0.05))
        delays = {}
        for ratio, r in cells:
            for cd in E11F.C_DOWNS:
                if baseline == "plugin":
                    tu, td = E11F.first_alarms(r, cd)
                else:
                    tu, td = first_alarms_cs(r, cd)
                det = td >= 0
                d = np.where(det, td + 1, HOR)
                delays[(ratio, cd)] = d
                diff = d - delays[(ratio, E11F.C_DOWNS[0])]
                rows.append(dict(baseline=baseline, A=E11F.A if baseline == "plugin" else 36000.0, c_down=cd,
                                 scenario=f"sd x{ratio}", paths=r.shape[0], detected=float(det.mean()),
                                 up_alarm_first=float(np.mean(tu >= 0)), median_delay=float(np.median(d)),
                                 mean_delay=float(d.mean()), se_mean_delay=float(d.std(ddof=1) / np.sqrt(d.size)),
                                 mean_diff_vs_c3=float(diff.mean()),
                                 se_diff_vs_c3=float(diff.std(ddof=1) / np.sqrt(d.size))))
        print("part 4", baseline, "done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "e19_split_c_cs.csv", index=False)
    with pd.option_context("display.width", 250):
        print(out.round(4).to_string(index=False))
    return out


# ----------------------------------------------------------------------------
# Part 2 [REAL]: E6 with the confidence-sequence baseline (M1(a))
# ----------------------------------------------------------------------------

P2_BANDS = (0.0, 0.25, 0.5)
A_CS_CRYPTO = 36000.0  # 1440 / (0.05 - delta), as in E8
CS_NAME = "robust e-detector, cs"


def _cost_rows(mon, n: int, offset: int) -> dict:
    """Monitored minutes and delays from the end of calibration for one monitor run over n returns."""
    ev_by_start = {e.monitored_since: e for e in mon.events}
    monitored = 0
    for seg in mon.segments:
        e = ev_by_start.get(seg["start"])
        monitored += (e.index - seg["start"] + 1) if e is not None else (n - seg["start"])
    delays = np.array([e.index - e.monitored_since for e in mon.events], float)
    return dict(minutes=n, monitored_minutes=int(monitored), monitored_share=monitored / n,
                median_minutes_calibration_end_to_alarm=float(np.median(delays)) if delays.size else np.nan,
                calibrations=len(mon.segments))


def p2_job(job):
    import e06b_five_coins_monitor as E6
    from avmonitor.monitor import VolatilityMonitor

    source, period = job
    idx, u, i0 = E6.prepare(source, period)
    out, cost = {}, []
    for baseline in ("cs", "plugin"):
        for band in P2_BANDS:
            if baseline == "cs":
                mon = VolatilityMonitor(threshold=A_CS_CRYPTO, c=C, cal_len=1440, null_band=band,
                                        baseline="cs", delta=DELTA)
            else:
                mon = VolatilityMonitor(threshold=E6.A, c=E6.C, cal_len=E6.CAL, null_band=band)
            ev = mon.run(u[i0:])
            out[(baseline, band)] = [(e.index + i0, e.direction, e.monitored_since + i0, e.baseline_scale) for e in ev]
            cost.append(dict(source=source, period=period, baseline=baseline, band=band, alarms=len(ev),
                             up=sum(e.direction == "up" for e in ev), down=sum(e.direction == "down" for e in ev),
                             **_cost_rows(mon, len(u) - i0, i0)))
    print(f"p2 {source} {period} done", flush=True)
    return source, period, out, cost


def _pool_monitor(summ: pd.DataFrame, coins) -> pd.DataFrame:
    """E6's pooling of the Q6d validation (``pooled_monitor_summary``), applied to a given summary table."""
    m = summ[summ.source.isin(list(coins)) & ~summ.detector.str.startswith("control")]
    rows = []
    for (period, det, band), g in m.groupby(["period", "detector", "band"], sort=False):
        rec = dict(period=period, detector=det, band=band, coins=len(g), alarms_per_day_min=g.alarms_per_day.min(),
                   alarms_per_day_max=g.alarms_per_day.max(), up=int(g.up.sum()), down=int(g.down.sum()))
        for dn in ("up", "down"):
            n = int(g[f"{dn}_validated_n"].sum())
            k = int(np.round(g[f"{dn}_validated_n"] * g[f"{dn}_share_in_flagged_direction"].fillna(0)).sum())
            p_, lo, hi = wilson(k, n)
            rec.update({f"{dn}_n": n, f"{dn}_hits": k, f"{dn}_share_in_flagged_direction": p_, f"{dn}_wilson_lo": lo,
                        f"{dn}_wilson_hi": hi,
                        f"{dn}_matched_delay_control_mean": g[f"{dn}_matched_delay_control_share"].mean()})
        rows.append(rec)
    return pd.DataFrame(rows)


def _pool_hits(h: pd.DataFrame, coins, edet_rules) -> pd.DataFrame:
    """E6's pooling of the forecast hits (``pooled_summary``), with the e-detector rows named in ``edet_rules``."""
    h = h[h.source.isin(list(coins))].copy()
    h["short_window"] = h.short_window.fillna(0).astype(int)
    h["band"] = h.band.fillna(-1)
    keys = ["period", "forecast", "rule", "band", "short_window"] + (["config"] if "config" in h else [])
    g = h.groupby(keys, sort=False)[["up_hits", "n_up", "down_hits", "n_down"]].sum().reset_index()
    for dn in ("up", "down"):
        w = [wilson(k, n) for k, n in zip(g[f"{dn}_hits"], g[f"n_{dn}"])]
        g[f"{dn}_hit_rate"] = [x[0] for x in w]
        g[f"{dn}_wilson_lo"] = [x[1] for x in w]
        g[f"{dn}_wilson_hi"] = [x[2] for x in w]
    g["band"] = g.band.replace(-1, np.nan)
    return g


def _day_block_ci(ev: pd.DataFrame, days: np.ndarray, rng, B: int = 2000):
    """Day-block bootstrap of validated shares: calendar days resampled with replacement, all coins jointly.

    ``ev``: rows with columns day, set (label), direction, hit (0/1). Returns {(set, direction): (share, lo, hi,
    draws)}."""
    out = {}
    pos = {d: i for i, d in enumerate(days)}
    W = rng.multinomial(len(days), np.full(len(days), 1.0 / len(days)), size=B)  # counts per day, B x D
    for (lab, dn), g in ev.groupby(["set", "direction"]):
        hits = np.zeros(len(days))
        n = np.zeros(len(days))
        for d, h in zip(g.day, g.hit):
            hits[pos[d]] += h
            n[pos[d]] += 1
        draws = (W @ hits) / np.maximum(W @ n, 1)
        out[(lab, dn)] = (hits.sum() / n.sum(), *np.quantile(draws, [0.025, 0.975]), draws)
    return out


def primary_bootstrap(evs: pd.DataFrame) -> pd.DataFrame:
    """Primary quantity of Part 2: day-block bootstrap, Jul to Dec 2025, band 0.25, five Binance coins.

    95% percentile intervals (pre-registered) and, added post hoc as a multiplicity statement for the two
    directions, Bonferroni 97.5% intervals (columns lo_975, hi_975) from the same draws."""
    import e06b_five_coins_monitor as E6

    e6ev = pd.read_csv(ROOT / "results" / "e06d_events.csv")
    sel = lambda d, det: d[(d.period == "Jul-Dec 2025") & (d.detector == det) & (d.band == 0.25)
                           & d.source.isin(list(E6.COINS)) & np.isfinite(d.post_over_baseline)]
    rows = []
    for lab, d in (("cs", sel(evs, CS_NAME)), ("plugin", sel(e6ev, "robust e-detector"))):
        day = pd.to_datetime(d.time, utc=True).dt.floor("D").dt.strftime("%Y-%m-%d")
        hit = np.where(d.direction == "up", d.post_over_baseline > 1, d.post_over_baseline < 1).astype(int)
        rows.append(pd.DataFrame(dict(day=day.values, set=lab, direction=d.direction.values, hit=hit)))
    bev = pd.concat(rows, ignore_index=True)
    days = pd.date_range("2025-07-01", "2025-12-31", freq="D").strftime("%Y-%m-%d").to_numpy()
    ci = _day_block_ci(bev, days, np.random.default_rng(1910))
    brow = []
    q975 = lambda x: (float(np.quantile(x, 0.0125)), float(np.quantile(x, 0.9875)))
    for dn in ("up", "down"):
        for lab in ("cs", "plugin"):
            sh, lo, hi, dr = ci[(lab, dn)]
            brow.append(dict(quantity=f"{dn} share validated, {lab}", estimate=sh, lo=lo, hi=hi,
                             **dict(zip(("lo_975", "hi_975"), q975(dr)))))
        dd = ci[("cs", dn)][3] - ci[("plugin", dn)][3]
        brow.append(dict(quantity=f"{dn} share validated, cs minus plugin",
                         estimate=ci[("cs", dn)][0] - ci[("plugin", dn)][0], lo=float(np.quantile(dd, 0.025)),
                         hi=float(np.quantile(dd, 0.975)), **dict(zip(("lo_975", "hi_975"), q975(dd)))))
    boot = pd.DataFrame(brow)
    boot.to_csv(OUT / "e19_cs_primary_bootstrap.csv", index=False)
    return boot


def part2(workers: int = 2):
    import e06b_five_coins_monitor as E6

    jobs = [(c, p) for p in E6.PERIODS for c in E6.COINS] + [("Bitstamp BTC", "Jul-Dec 2025")]
    with Pool(workers) as pool:
        res = pool.map(p2_job, jobs, chunksize=1)
    events = {(s, p): o for s, p, o, _ in res}
    cost = pd.DataFrame([c for *_, cs in res for c in cs])
    rng_d, rng_e = np.random.default_rng(1906), np.random.default_rng(1907)
    summ, evs, hits, comps = [], [], [], []
    for src, period in jobs:
        idx, u, i0 = E6.prepare(src, period)
        S = E6.Series(idx, u, i0)
        ev_cs = {(CS_NAME, b): events[(src, period)][("cs", b)] for b in P2_BANDS}
        r, e = E6.summarise(src, period, S, ev_cs, rng_d)
        summ += r
        evs += e
        h, cp = E6.forecast_hits(src, period, S, ev_cs, rng_e)
        hits.append(h)
        comps.append(cp)
        print(f"p2 validation {src} {period} done", flush=True)
        del S
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "e19_cs_monitor_summary.csv", index=False)
    evs = pd.DataFrame(evs)
    evs.to_csv(OUT / "e19_cs_events.csv", index=False)
    hits = pd.concat(hits, ignore_index=True)
    hits.to_csv(OUT / "e19_cs_forecast_hits.csv", index=False)
    pd.concat(comps, ignore_index=True).to_csv(OUT / "e19_cs_competitor.csv", index=False)

    # plug-in rerun: alarm counts must equal E6's (deterministic); cost of the confidence sequence
    e6 = pd.read_csv(ROOT / "results" / "e06d_monitor_summary.csv")
    e6 = e6[e6.detector == "robust e-detector"]
    chk = cost[cost.baseline == "plugin"].merge(e6[["source", "period", "band", "alarms", "up", "down"]],
                                                on=["source", "period", "band"], suffixes=("", "_e6"))
    chk["reproduces_e6"] = (chk.alarms == chk.alarms_e6) & (chk.up == chk.up_e6) & (chk.down == chk.down_e6)
    print("plug-in rerun reproduces E6 alarm counts:", bool(chk.reproduces_e6.all()), flush=True)
    lag_rows = []
    for (src, period), o in events.items():
        for band in P2_BANDS:
            pl, cs = o[("plugin", band)], o[("cs", band)]
            for dn in ("up", "down"):
                tp = np.array([e[0] for e in pl if e[1] == dn])
                tc = np.array([e[0] for e in cs if e[1] == dn])
                lags = []
                for t in tp:
                    d = tc - t
                    d = d[np.abs(d) <= 360]
                    if d.size:
                        lags.append(d[np.argmin(np.abs(d))])
                lag_rows.append(dict(source=src, period=period, band=band, direction=dn, plugin_alarms=tp.size,
                                     cs_alarms=tc.size, matched_within_6h=len(lags),
                                     share_plugin_without_cs_within_6h=1 - len(lags) / tp.size if tp.size else np.nan,
                                     median_cs_minus_plugin_minutes=float(np.median(lags)) if lags else np.nan,
                                     mean_cs_minus_plugin_minutes=float(np.mean(lags)) if lags else np.nan))
    lag = pd.DataFrame(lag_rows)
    cost = cost.merge(chk[["source", "period", "band", "baseline", "reproduces_e6"]],
                      on=["source", "period", "band", "baseline"], how="left")
    cost.to_csv(OUT / "e19_cs_cost.csv", index=False)
    lag.to_csv(OUT / "e19_cs_lag.csv", index=False)

    # pooled tables, cs next to E6's plug-in
    pooled_cs = _pool_monitor(summ, E6.COINS)
    pooled_pl = pd.read_csv(ROOT / "results" / "e06d_pooled.csv")
    pooled_pl = pooled_pl[pooled_pl.detector == "robust e-detector"]
    pooled = pd.concat([pooled_cs.assign(baseline="cs", source_file="e19"),
                        pooled_pl.assign(baseline="plugin", source_file="e06d_pooled.csv")], ignore_index=True)
    pooled.to_csv(OUT / "e19_cs_pooled.csv", index=False)
    ph = _pool_hits(hits, E6.COINS, [CS_NAME])
    e6h = pd.read_csv(ROOT / "results" / "e06e_pooled.csv")
    ph = pd.concat([ph.assign(baseline="cs"), e6h.assign(baseline="plugin (E6)")], ignore_index=True)
    ph.to_csv(OUT / "e19_cs_forecast_pooled.csv", index=False)

    boot = primary_bootstrap(evs)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(pooled.round(3).to_string(index=False))
        print(boot.round(4).to_string(index=False))
        print(cost.round(3).to_string(index=False))
        print(lag.round(1).to_string(index=False))
    return pooled


# ----------------------------------------------------------------------------
# Part 3 [REAL]: the forecast-relative monitor of E10 on the five coins (M1(c))
# ----------------------------------------------------------------------------

FR_BANDS = (0.0, 0.25)
FR_TRAIN = ("2025-01-01", "2025-07-01")


def forecast_inputs(idx, u):
    """z = u / sigma-hat for the GARCH and the EWMA (half-life 360) forecasts (predictable, forecast.py)."""
    from avmonitor import forecast as F

    j0 = int(np.searchsorted(idx, pd.Timestamp(FR_TRAIN[0], tz="UTC")))
    j1 = int(np.searchsorted(idx, pd.Timestamp(FR_TRAIN[1], tz="UTC")))
    omega, a, b = F.garch11_fit(u[j0:j1])
    s2 = np.full(u.size, np.nan)
    s2[j0:] = F.garch11_filter(u[j0:], omega, a, b)
    f = int(np.argmax(np.isfinite(u)))
    head = u[f:][np.isfinite(u[f:])][:1440]
    v = np.full(u.size, np.nan)
    v[f:] = F.ewma_variance(u[f:], 360, init=float(np.mean(head**2)))
    return {"GARCH": F.standardise(u, s2), "EWMA 360": F.standardise(u, v)}, dict(omega=omega, alpha=a, beta=b)


def p3_job(coin):
    import e06b_five_coins_monitor as E6
    from avmonitor.monitor import VolatilityMonitor

    idx, u, i0 = E6.prepare(coin, "Jul-Dec 2025")
    Z, garch = forecast_inputs(idx, u)
    out, cost = {}, []
    for inp, z in Z.items():
        for baseline in ("plugin", "cs"):
            for band in FR_BANDS:
                if baseline == "cs":
                    mon = VolatilityMonitor(threshold=A_CS_CRYPTO, c=C, cal_len=1440, null_band=band,
                                            baseline="cs", delta=DELTA)
                else:
                    mon = VolatilityMonitor(threshold=E6.A, c=C, cal_len=1440, null_band=band)
                ev = mon.run(z[i0:])
                seg_len = {s["start"]: s["cal_len"] for s in mon.segments}
                name = f"robust e-detector, {inp}-relative, {baseline}"
                # baseline on the scale of u: robust vol of u over the alarm's calibration window
                out[(name, band)] = [(e.index + i0, e.direction, e.monitored_since + i0,
                                      E6.window_sd(u, e.monitored_since + i0 - seg_len[e.monitored_since],
                                                   e.monitored_since + i0)) for e in ev]
                c = _cost_rows(mon, len(u) - i0, i0)
                days_mon = c["monitored_minutes"] / 1440
                cost.append(dict(source=coin, input=inp, baseline=baseline, band=band, alarms=len(ev),
                                 up=sum(e.direction == "up" for e in ev), down=sum(e.direction == "down" for e in ev),
                                 alarms_per_day=len(ev) / ((len(u) - i0) / 1440), alarms_per_monitored_day=len(ev) / days_mon,
                                 bound_per_monitored_day=(0.05 if baseline == "plugin" else
                                                          (DELTA + days_mon * 1440 / 36000.0) / (1 - DELTA) / days_mon),
                                 **c, **garch))
    print(f"p3 {coin} done: GARCH {garch}", flush=True)
    return coin, out, cost


def part3(workers: int = 2):
    import e06b_five_coins_monitor as E6

    with Pool(workers) as pool:
        res = pool.map(p3_job, list(E6.COINS), chunksize=1)
    cost = pd.DataFrame([c for *_, cs in res for c in cs])
    cost.to_csv(OUT / "e19_fr_cost.csv", index=False)
    rng_d, rng_e = np.random.default_rng(1908), np.random.default_rng(1909)
    summ, evs, hits = [], [], []
    for coin, out, _ in res:
        idx, u, i0 = E6.prepare(coin, "Jul-Dec 2025")
        S = E6.Series(idx, u, i0)
        r, e = E6.summarise(coin, "Jul-Dec 2025", S, out, rng_d)
        summ += r
        evs += e
        for key, ev in out.items():
            h, _ = E6.forecast_hits(coin, "Jul-Dec 2025", S, {key: ev}, rng_e)
            hits.append(h.assign(config=key[0]))
        print(f"p3 validation {coin} done", flush=True)
        del S
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "e19_fr_monitor_summary.csv", index=False)
    pd.DataFrame(evs).to_csv(OUT / "e19_fr_events.csv", index=False)
    hits = pd.concat(hits, ignore_index=True)
    hits.to_csv(OUT / "e19_fr_forecast_hits.csv", index=False)
    pooled = _pool_monitor(summ, E6.COINS)
    pooled.to_csv(OUT / "e19_fr_pooled.csv", index=False)
    ph = _pool_hits(hits[hits.rule.str.startswith("robust") | (hits.rule == "random times")], E6.COINS, [])
    ph.to_csv(OUT / "e19_fr_forecast_pooled.csv", index=False)
    agg = cost.groupby(["input", "baseline", "band"]).agg(
        alarms_per_day_min=("alarms_per_day", "min"), alarms_per_day_max=("alarms_per_day", "max"),
        per_monitored_day_min=("alarms_per_monitored_day", "min"), per_monitored_day_max=("alarms_per_monitored_day", "max"),
        monitored_share_min=("monitored_share", "min"), monitored_share_max=("monitored_share", "max"),
        median_delay_min=("median_minutes_calibration_end_to_alarm", "min"),
        median_delay_max=("median_minutes_calibration_end_to_alarm", "max")).reset_index()
    agg.to_csv(OUT / "e19_fr_cost_summary.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(agg.round(3).to_string(index=False))
        print(pooled.round(3).to_string(index=False))
        print(ph[ph.forecast == "EWMA"].round(3).to_string(index=False))
    return pooled


# ----------------------------------------------------------------------------
# Post hoc diagnostic (added after seeing Part 2): segments whose increase side cannot alarm
# ----------------------------------------------------------------------------
# With the cs baseline tau is fixed from the first quarter of the calibration window and m_up = U (1 + band) from
# the rest. If tau^2 <= m_up, every increase factor 1 + lam (X / m_up - 1) is at most 1 because X <= tau^2, so
# the increase statistic cannot reach 2A in any realistic time: the segment is blind to increases until the next
# alarm, which can then only be a decrease alarm. The same can happen to the decrease side (down_active False).

def side_activity(mon, n: int) -> dict:
    """Monitored minutes in segments with the increase side (tau^2 <= m_up) or decrease side switched off."""
    ev_by_start = {e.monitored_since: e for e in mon.events}
    mon_total = up_off = dn_off = 0
    longest_up_off = 0
    b_up = []
    for seg in mon.segments:
        e = ev_by_start.get(seg["start"])
        length = (e.index - seg["start"] + 1) if e is not None else (n - seg["start"])
        mon_total += length
        ratio = seg["tau2"] / seg["m_up"] if seg["m_up"] > 0 else math.inf
        b_up.append(ratio)
        if ratio <= 1.0:
            up_off += length
            longest_up_off = max(longest_up_off, length)
        if not seg["down_active"]:
            dn_off += length
    b_up = np.array(b_up)
    return dict(segments=len(mon.segments), segments_up_off=int((b_up <= 1).sum()),
                monitored_minutes=int(mon_total), up_off_minutes=int(up_off),
                up_off_share_of_monitored=up_off / mon_total if mon_total else np.nan,
                longest_up_off_days=longest_up_off / 1440, down_off_minutes=int(dn_off),
                median_tau2_over_m_up=float(np.median(b_up)) if b_up.size else np.nan,
                q10_tau2_over_m_up=float(np.quantile(b_up, 0.1)) if b_up.size else np.nan)


def blind_job(job):
    import e06b_five_coins_monitor as E6
    from avmonitor.monitor import VolatilityMonitor

    part, source, period = job
    idx, u, i0 = E6.prepare(source, period)
    inputs = {"returns": u} if part == "p2" else forecast_inputs(idx, u)[0]
    rows = []
    for inp, z in inputs.items():
        for baseline in ("cs", "plugin"):
            for band in (P2_BANDS if part == "p2" else FR_BANDS):
                if baseline == "cs":
                    mon = VolatilityMonitor(threshold=A_CS_CRYPTO, c=C, cal_len=1440, null_band=band,
                                            baseline="cs", delta=DELTA)
                else:
                    mon = VolatilityMonitor(threshold=E6.A, c=C, cal_len=1440, null_band=band)
                ev = mon.run(z[i0:])
                rows.append(dict(part=part, source=source, period=period, input=inp, baseline=baseline, band=band,
                                 alarms=len(ev), **side_activity(mon, len(u) - i0)))
    print(f"blind {part} {source} {period} done", flush=True)
    return rows


def part_blind(workers: int = 1):
    import e06b_five_coins_monitor as E6

    jobs = ([("p2", c, p) for p in E6.PERIODS for c in E6.COINS] + [("p2", "Bitstamp BTC", "Jul-Dec 2025")]
            + [("p3", c, "Jul-Dec 2025") for c in E6.COINS])
    with Pool(workers) as pool:
        res = pool.map(blind_job, jobs, chunksize=1)
    d = pd.DataFrame([r for rows in res for r in rows])
    d.to_csv(OUT / "e19_blind.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 300):
        print(d.round(3).to_string(index=False))
    return d


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

MON_STYLE = {  # known baseline: neutral reference; the others: validated categorical slots 1 to 3 + markers
    "known": ("#52514e", "o", "(i) known baseline"),
    "plugin": ("#eb6834", "s", "(ii) plug-in"),
    "cs": ("#2a78d6", "D", "(iii) confidence sequence"),
    "rank": ("#1baf7a", "^", "(iv) rank, restarted at once"),
}
STREAM_LABEL = {"perm": "permuted", "iid": "i.i.d.", "sb30": "block 30", "sb390": "block 390", "sb1440": "block 1440"}


def _axes_style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color="#e6e5e1", lw=0.5)
    ax.set_axisbelow(True)


def fig_null(df: pd.DataFrame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sets = [("five coins", "Five coins pooled"), ("Nifty 50", "Nifty 50")]
    cols = [("one_day", "P(false alarm within the first day)"), ("five_days", "P(false alarm within five days)"),
            ("per_monitored_day", "False alarms per monitored day, restarts")]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4), sharex=True)
    x = np.arange(len(STREAM_TYPES))
    for i, (ds, title) in enumerate(sets):
        for j, (col, lab) in enumerate(cols):
            ax = axes[i, j]
            _axes_style(ax)
            d = df[(df.dataset == ds) & (df.band == 0.0)]
            for k, mon in enumerate(MONITORS):
                g = d[d.monitor == mon].set_index("stream").reindex(list(STREAM_TYPES))
                col_, mk, name = MON_STYLE[mon]
                xs = x + (k - 1.5) * 0.17
                v, lo, hi = g[col].to_numpy(), g[f"{col}_lo"].to_numpy(), g[f"{col}_hi"].to_numpy()
                ax.errorbar(xs, v, yerr=[np.maximum(v - lo, 0), np.maximum(hi - v, 0)], fmt=mk, color=col_, ms=6,
                            mec="white", mew=0.8, elinewidth=1, capsize=0, label=name)
                bnd = g[f"{col}_bound"].to_numpy()
                ax.hlines(bnd, xs - 0.08, xs + 0.08, color=col_, lw=1.6)
            ax.set_yscale("symlog", linthresh=1e-3, linscale=0.4)  # the cs rates are near 0
            ax.set_ylim(0, 1.5 if col != "per_monitored_day" else 10)
            ax.axvspan(-0.5, 1.5, color="#f3f2ee", zorder=0)
            ax.set_xlim(-0.5, len(STREAM_TYPES) - 0.5)
            ax.set_xticks(x)
            ax.set_xticklabels([STREAM_LABEL[s] for s in STREAM_TYPES], fontsize=8)
            ax.set_title(f"{title}: {lab}", fontsize=9, loc="left")
            if j == 0:
                ax.set_ylabel("rate (symmetric log scale)", fontsize=8)
            ax.tick_params(labelsize=8)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Null streams with the real marginal law, band 0 [REAL marginal, constructed streams]. Points: rate with "
                 "95% interval; short bars: each monitor's bound (nominal only for the plug-in).\nShaded: the null "
                 "holds (permuted pool, i.i.d. draws from it). Unshaded: stationary bootstrap of the calm days in "
                 "calendar order (mean block 30, 390, 1440 returns).", fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    fig.savefig(FIG / "e19_null_rates.png", dpi=130)
    plt.close(fig)


def fig_real(pooled: pd.DataFrame, cost: pd.DataFrame, fr: pd.DataFrame):
    """Alarms per day and validated shares: plug-in and cs on raw deseasonalised returns, forecast-relative inputs."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
    per = "Jul-Dec 2025"
    rows = []
    for base, det in (("plugin", "robust e-detector"), ("cs", CS_NAME)):
        for band in P2_BANDS:
            r = pooled[(pooled.period == per) & (pooled.detector == det) & (pooled.band == band)]
            if len(r):
                rows.append((f"raw, {base}, b={band:g}", base, r.iloc[0]))
    for inp in ("GARCH", "EWMA 360"):
        for base in ("plugin", "cs"):
            for band in FR_BANDS:
                det = f"robust e-detector, {inp}-relative, {base}"
                r = fr[(fr.detector == det) & (fr.band == band)]
                if len(r):
                    rows.append((f"{inp}, {base}, b={band:g}", base, r.iloc[0]))
    y = np.arange(len(rows))[::-1]
    for ax, (col, lab) in zip(axes[1:], (("up", "Up alarms followed by higher vol"),
                                          ("down", "Down alarms followed by lower vol"))):
        _axes_style(ax)
        for yy, (name, base, r) in zip(y, rows):
            c_, mk, _ = MON_STYLE[base]
            v = r[f"{col}_share_in_flagged_direction"]
            ax.errorbar([v], [yy], xerr=[[v - r[f"{col}_wilson_lo"]], [r[f"{col}_wilson_hi"] - v]], fmt=mk, color=c_,
                        ms=6, mec="white", mew=0.8, elinewidth=1)
        ax.axvline(0.5, color="#9a9893", lw=0.8, ls="--")
        ax.set_xlim(0.3, 1.03)
        ax.set_yticks(y)
        ax.set_yticklabels([])
        ax.set_title(f"{lab}: pooled share, 95% Wilson", fontsize=9, loc="left")
        ax.tick_params(labelsize=8)
    ax = axes[0]
    _axes_style(ax)
    for yy, (name, base, r) in zip(y, rows):
        c_, mk, _ = MON_STYLE[base]
        ax.plot([r["alarms_per_day_min"], r["alarms_per_day_max"]], [yy, yy], color=c_, lw=3, solid_capstyle="round")
    ax.axvline(0.05, color="#0b0b0b", lw=0.8, ls=":")
    ax.text(0.06, y[-1] - 0.45, "0.05", fontsize=7, color="#52514e")
    ax.set_xlim(0, 0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([n for n, _, _ in rows], fontsize=8)
    ax.set_title("Alarms per calendar day, range over the five coins", fontsize=9, loc="left")
    ax.tick_params(labelsize=8)
    fig.suptitle("Five Binance coins, July to December 2025 [REAL]. Orange: plug-in baseline; blue: confidence-sequence "
                 "baseline. raw = deseasonalised returns; GARCH, EWMA 360 = divided by a forecast.",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG / "e19_real_cs_forecast.png", dpi=130)
    plt.close(fig)


def figures():
    FIG.mkdir(parents=True, exist_ok=True)
    if (OUT / "e19_null_streams.csv").exists():
        fig_null(pd.read_csv(OUT / "e19_null_streams.csv"))
    if (OUT / "e19_cs_pooled.csv").exists() and (OUT / "e19_fr_pooled.csv").exists():
        fig_real(pd.read_csv(OUT / "e19_cs_pooled.csv"), pd.read_csv(OUT / "e19_cs_cost.csv"),
                 pd.read_csv(OUT / "e19_fr_pooled.csv"))


# ----------------------------------------------------------------------------


def main(argv):
    parts = argv or ["p1", "p2", "p3", "p4", "blind", "p1std", "figs"]
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if "pools" in parts:
        print(prepare_pools(tuple(DATASETS)).to_string(index=False))
    if "p1" in parts:
        part1()
    if "p1pooled" in parts:
        p1_pooled()
    if "p4" in parts:
        part4()
    if "p2" in parts:
        part2()
    if "p2boot" in parts:  # recompute the primary bootstrap from the saved events (same seed, same draws)
        print(primary_bootstrap(pd.read_csv(OUT / "e19_cs_events.csv")).round(4).to_string(index=False))
    if "p3" in parts:
        part3()
    if "blind" in parts:
        part_blind()
    if "p1std" in parts:
        part_std()
    if "figs" in parts:
        figures()
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main(sys.argv[1:])
