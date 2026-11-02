"""Conformal rank betting: a sequential "changed versus reference" test for the law of returns.

The main package bets on the truncated squared return ``X_t = min(r_t^2, tau^2)``
against a baseline mean ``m0`` estimated from a calibration window. This module
replaces "compare a moment with an estimated baseline" by "compare the new
return with every return seen so far in the segment, calibration window
included", which is a sequential two-sample test with a reference sample that
grows by one point per step.

Construction (proofs in README.md, section 2)
--------------------------------------------
1. Score. ``s_t = |u_t|`` (``kind="abs"``, sees variance, tails and jumps) or
   ``s_t = u_t`` (``kind="signed"``, also sees skew), where ``u_t`` is a
   (deseasonalised) return. No scale, threshold or mean is estimated.
2. Smoothed conformal p-value. With ``n`` earlier scores in the segment,
   ``G`` of them larger than ``s_t`` and ``E - 1`` of them equal to it,
   ``p_t = (G + U_t E) / (n + 1)`` with ``U_t`` uniform on (0, 1). Small ``p``
   means "larger than almost everything seen so far". If the segment is
   exchangeable (for example i.i.d. with any law, ties allowed), the p-values
   are independent and exactly uniform (Vovk, Nouretdinov and Gammerman 2003).
3. Bets. Each bettor pays ``f_t = 1 + lam (w_t(p_t) - b E_U|w_t|)`` where
   ``w_t`` is a function of PAST p-values only, has mean 0 under the uniform
   law and is bounded below. Three families:

   * scale up / scale down: ``w = +-(Y - 1)`` with the Gaussianised rank score
     ``Y = min(z^2, a^2) / kappa_2(a)``, ``z`` the standard normal quantile of
     ``p``. ``E_U[Y] = 1`` holds EXACTLY, so this is the main monitor's
     truncated-variance bet with a baseline that is known instead of estimated.
     The bets reuse the main package's grids.
   * histogram: ``w = pi_t / q - 1`` where ``q`` are the widths of fixed bins
     of [0, 1] and ``pi_t`` is an exponentially weighted, prior-smoothed
     frequency of past p-values in each bin (a plug-in density bettor in the
     spirit of Fedorova et al. 2012). It learns the shape of any departure.

4. Detector. Shiryaev-Roberts recursion ``S_t = f_t (S_{t-1} + 1)`` per bettor,
   ``M_t`` = weighted average over bettors, alarm when ``M_t >= A``. Under the
   null ``M_t - t`` is a martingale in the filtration of the p-values, so the
   average run length is at least ``A`` and P(alarm within m steps) <= m / A
   (Theorems A and B of docs/math.md, section 5.4, apply verbatim).

``band = b > 0`` widens the null to "the conditional density of each p-value
stays within a factor 1 +/- b of uniform", which tolerates mild dependence such
as short volatility clusters (README, section 2.5).
"""

from __future__ import annotations

import bisect
import functools
import math
from dataclasses import dataclass, field

import numpy as np
from scipy import integrate
from scipy.special import ndtri

from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from avmonitor.robust import truncated_moment_factor

# Bin edges for the histogram bettor: finer at both ends (for |u| scores, small p
# = large moves and large p = near-zero moves; for signed scores, the two tails).
DEFAULT_EDGES = (0.0, 0.002, 0.01, 0.03, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97, 0.99, 0.998, 1.0)
FAMILIES = ("scale_up", "scale_down", "histogram")


# ----------------------------------------------------------------------------
# Scores and conformal p-values
# ----------------------------------------------------------------------------


def score(u, kind: str = "abs"):
    """Conformity score of a return: ``|u|`` for ``kind="abs"``, ``u`` for ``"signed"``."""
    if kind == "abs":
        return np.abs(u)
    if kind == "signed":
        return np.asarray(u, float)
    raise ValueError("kind must be 'abs' or 'signed'")


def pvalues_against_reference(s: np.ndarray, n_ref: int, rng: np.random.Generator) -> np.ndarray:
    """Smoothed conformal p-values of ``s[:, n_ref:]``, each ranked among ALL earlier values of its row.

    Entry ``(i, t)`` is ``(G + U E) / (n_ref + t + 1)`` where ``G`` counts earlier
    values of row i (reference block and monitored values before t) that are
    larger than ``s[i, n_ref + t]`` and ``E`` counts the equal ones plus itself.
    Vectorised over rows with a Fenwick (binary indexed) tree, O(n log n) per row.
    Rows must not contain NaN.
    """
    s = np.atleast_2d(np.asarray(s, float))
    reps, n = s.shape
    m = n - n_ref
    if m <= 0:
        return np.empty((reps, 0))
    if np.isnan(s).any():
        raise ValueError("scores contain NaN")
    ref = np.sort(s[:, :n_ref], axis=1)
    mon = s[:, n_ref:]
    g_ref = np.empty((reps, m))
    e_ref = np.empty((reps, m))
    dense = np.empty((reps, m), np.int64)
    n_unique = 1
    for i in range(reps):
        hi = np.searchsorted(ref[i], mon[i], side="right")
        lo = np.searchsorted(ref[i], mon[i], side="left")
        g_ref[i] = n_ref - hi
        e_ref[i] = hi - lo
        uniq, inv = np.unique(mon[i], return_inverse=True)
        dense[i] = inv + 1  # 1-based ranks for the Fenwick tree
        n_unique = max(n_unique, uniq.size)
    levels = int(n_unique).bit_length()
    top = 1 << levels
    tree = np.zeros((reps, top + 2), np.int64)  # column top + 1 is a dump for out-of-range updates
    rows = np.arange(reps)
    g_mon = np.empty((reps, m))
    e_mon = np.empty((reps, m))

    def prefix(idx):
        acc = np.zeros(reps, np.int64)
        idx = idx.copy()
        for _ in range(levels + 1):
            acc += tree[rows, idx]  # tree[:, 0] stays 0
            idx -= idx & -idx
        return acc

    for t in range(m):
        v = dense[:, t]
        le = prefix(v)
        lt = prefix(v - 1)
        g_mon[:, t] = t - le
        e_mon[:, t] = le - lt
        idx = v.copy()
        for _ in range(levels + 1):
            tree[rows, np.where(idx <= top, idx, top + 1)] += 1
            idx = np.where(idx <= top, idx + (idx & -idx), idx)
    u = rng.random((reps, m))
    denom = n_ref + np.arange(m) + 1.0
    return (g_ref + g_mon + u * (e_ref + e_mon + 1.0)) / denom[None, :]


class SequentialRanker:
    """Streaming version of the smoothed conformal p-value (exact, with ties)."""

    def __init__(self, rng: np.random.Generator):
        self.rng = rng
        self.values: list[float] = []

    def add(self, s: float) -> None:
        bisect.insort(self.values, s)

    def pvalue(self, s: float, insert: bool = True, u: float | None = None) -> float:
        n = len(self.values)
        hi = bisect.bisect_right(self.values, s)
        lo = bisect.bisect_left(self.values, s)
        u = self.rng.random() if u is None else u
        p = (n - hi + u * (hi - lo + 1)) / (n + 1)
        if insert:
            self.values.insert(hi, s)
        return p


# ----------------------------------------------------------------------------
# Gaussianised rank score (the rank analogue of min(r^2, tau^2) / m0)
# ----------------------------------------------------------------------------


@functools.lru_cache(maxsize=32)
def _kappa2(a: float) -> float:
    return truncated_moment_factor(a, power=2)


def gauss_score(p, kind: str = "abs", a: float = 3.0):
    """``Y = min(z^2, a^2) / kappa_2(a)``, ``z`` the normal quantile matching ``p``. ``E_U[Y] = 1`` exactly."""
    p = np.asarray(p, float)
    z = -ndtri(p / 2.0) if kind == "abs" else -ndtri(p)
    return np.minimum(z * z, a * a) / _kappa2(float(a))


def gauss_score_abs_dev(a: float = 3.0) -> float:
    """``E_U |Y - 1|`` for the Gaussianised score (the band penalty of the scale bettors)."""
    k2 = truncated_moment_factor(a, power=2)

    def f(z):
        return abs(min(z * z, a * a) / k2 - 1.0) * math.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)

    z1 = math.sqrt(k2)  # where Y crosses 1
    pts = [0.0, z1, a, 12.0]
    val = sum(integrate.quad(f, pts[i], pts[i + 1], limit=200)[0] for i in range(3))
    return 2.0 * val


# ----------------------------------------------------------------------------
# Betting core (vectorised over paths)
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class RankBetConfig:
    kind: str = "abs"  # score: "abs" or "signed"
    a: float = 3.0  # truncation of the Gaussianised score, in standard normal units
    up_grid: tuple = tuple(DEFAULT_UP_GRID)  # scale-up bets (as in the main monitor)
    down_fractions: tuple = tuple(DEFAULT_DOWN_FRACTIONS)  # scale-down bets as fractions of the largest
    edges: tuple = DEFAULT_EDGES  # histogram bins on [0, 1]
    half_lives: tuple = (240.0, 1440.0)  # memory of the histogram learners, in steps
    prior: float = 50.0  # pseudo-observations of the uniform law in each learner
    eps: float = 0.05  # each histogram bet keeps this share of wealth out of the game
    band: float = 0.0  # null band b (0 = exact exchangeability null)
    weights: tuple = (0.25, 0.25, 0.5)  # mixture weights of scale_up, scale_down, histogram

    def __post_init__(self):
        if self.kind not in ("abs", "signed"):
            raise ValueError("kind must be 'abs' or 'signed'")
        if len(self.weights) != 3 or min(self.weights) < 0 or not math.isclose(sum(self.weights), 1.0):
            raise ValueError("weights must be three nonnegative numbers summing to 1")
        e = np.asarray(self.edges, float)
        if e[0] != 0.0 or e[-1] != 1.0 or np.any(np.diff(e) <= 0):
            raise ValueError("edges must increase from 0 to 1")
        if self.band < 0 or not 0 < self.eps < 1:
            raise ValueError("band must be >= 0 and eps in (0, 1)")


class RankBetCore:
    """Shiryaev-Roberts mixture of rank bettors for ``reps`` independent streams of p-values."""

    def __init__(self, cfg: RankBetConfig, reps: int = 1):
        self.cfg = cfg
        self.reps = reps
        b = cfg.band
        self.k2 = truncated_moment_factor(cfg.a, power=2)
        self.B = cfg.a**2 / self.k2  # upper bound of Y
        self.mY = gauss_score_abs_dev(cfg.a)  # E_U|Y - 1|
        self.lam_up = np.asarray(cfg.up_grid, float) / (1.0 + b * self.mY)
        self.lam_dn = np.asarray(cfg.down_fractions, float) / ((self.B - 1.0) + b * self.mY)
        if np.any(self.lam_up * (1.0 + b * self.mY) >= 1.0):
            raise ValueError("scale-up bets must be below 1")
        if np.any(np.asarray(cfg.down_fractions) >= 1.0):
            raise ValueError("down fractions must be below 1")
        self.edges = np.asarray(cfg.edges, float)
        self.q = np.diff(self.edges)
        self.K = self.q.size
        self.gam = 0.5 ** (1.0 / np.asarray(cfg.half_lives, float))
        self.H = self.gam.size
        self.w = np.asarray(cfg.weights, float)
        self.reset()

    def reset(self, rows=None):
        if rows is None:
            self.S_up = np.zeros((self.reps, self.lam_up.size))
            self.S_dn = np.zeros((self.reps, self.lam_dn.size))
            self.S_h = np.zeros((self.reps, self.H))
            self.counts = np.zeros((self.reps, self.H, self.K))
        else:
            self.S_up[rows] = 0.0
            self.S_dn[rows] = 0.0
            self.S_h[rows] = 0.0
            self.counts[rows] = 0.0

    def density_ratio(self, learner: int = 0) -> np.ndarray:
        """Current estimate of (frequency of new p-values in each bin) / (bin width), shape ``(reps, K)``."""
        c = self.counts[:, learner, :]
        pi = (c + self.cfg.prior * self.q) / (c.sum(axis=1, keepdims=True) + self.cfg.prior)
        return pi / self.q

    def family_stats(self) -> np.ndarray:
        """Weighted contribution of each family to ``M`` (shape ``(reps, 3)``)."""
        return np.column_stack([
            self.w[0] * self.S_up.mean(axis=1),
            self.w[1] * self.S_dn.mean(axis=1),
            self.w[2] * self.S_h.mean(axis=1),
        ])

    def step(self, p: np.ndarray) -> np.ndarray:
        """Bet on one p-value per stream; returns the mixture statistic ``M_t`` per stream."""
        cfg, b = self.cfg, self.cfg.band
        p = np.asarray(p, float).reshape(self.reps)
        y = gauss_score(p, cfg.kind, cfg.a) - 1.0
        with np.errstate(over="ignore"):
            if self.w[0] > 0:
                self.S_up = (1.0 + self.lam_up[None, :] * (y[:, None] - b * self.mY)) * (self.S_up + 1.0)
            if self.w[1] > 0:
                self.S_dn = (1.0 + self.lam_dn[None, :] * (-y[:, None] - b * self.mY)) * (self.S_dn + 1.0)
            if self.w[2] > 0:
                k = np.clip(np.searchsorted(self.edges, p, side="right") - 1, 0, self.K - 1)
                tot = self.counts.sum(axis=2, keepdims=True)
                pi = (self.counts + cfg.prior * self.q) / (tot + cfg.prior)  # (reps, H, K)
                wk = pi / self.q - 1.0
                mdev = np.abs(pi - self.q).sum(axis=2)  # E_U |w|
                lam = (1.0 - cfg.eps) / (1.0 + b * mdev)
                rows = np.arange(self.reps)
                f = 1.0 + lam * (wk[rows, :, k] - b * mdev)
                self.S_h = f * (self.S_h + 1.0)
                self.counts *= self.gam[None, :, None]
                self.counts[rows, :, k] += 1.0
        return (self.w[0] * self.S_up.mean(axis=1) + self.w[1] * self.S_dn.mean(axis=1)
                + self.w[2] * self.S_h.mean(axis=1))


@dataclass
class RankResult:
    alarm: np.ndarray  # first alarm (monitored step index) per path, -1 if none
    share: np.ndarray  # share of M from (scale_up, scale_down, histogram) at the alarm
    path: np.ndarray | None = None  # M_t for path 0 if requested
    info: dict = field(default_factory=dict)


def run_rank_detector(p: np.ndarray, threshold: float, cfg: RankBetConfig = RankBetConfig(),
                      return_path: bool = False) -> RankResult:
    """First alarm of the rank e-detector on rows of p-values (monitoring starts at column 0)."""
    p = np.atleast_2d(np.asarray(p, float))
    reps, m = p.shape
    core = RankBetCore(cfg, reps)
    alarm = np.full(reps, -1)
    share = np.full((reps, 3), np.nan)
    active = np.ones(reps, bool)
    path = np.empty(m) if return_path else None
    for t in range(m):
        M = core.step(p[:, t])
        if path is not None:
            path[t] = M[0]
        hit = active & (M >= threshold)
        if hit.any():
            alarm[hit] = t
            fam = core.family_stats()[hit]
            share[hit] = fam / fam.sum(axis=1, keepdims=True)
            active &= ~hit
            if not active.any() and path is None:
                break
    return RankResult(alarm, share, path)


# ----------------------------------------------------------------------------
# Streaming monitor (calibrate, monitor, alarm, restart), same protocol as VolatilityMonitor
# ----------------------------------------------------------------------------


@dataclass
class RankEvent:
    index: int  # position of the alarm in the stream
    monitored_since: int  # first monitored index of the segment
    cal_start: int  # first index of the segment's reference (calibration) window
    stat: float  # M_t at the alarm
    share_up: float  # share of M_t from the scale-up bettors
    share_down: float
    share_hist: float
    mean_p_recent: float  # mean p-value over the last 240 monitored steps (0.5 under the null)
    density_ratio: tuple = ()  # learned p-value frequency / bin width per bin (shortest-memory learner); 1 = no change


class RankMonitor:
    """Streaming conformal rank e-detector.

    The first ``cal_len`` finite returns of a segment only enter the reference
    (their p-values are not bet on); every later return is ranked among all
    earlier returns of the segment and bet on. After an alarm a new segment
    starts. If the returns of a segment are exchangeable, P(alarm within m
    monitored steps) <= m / threshold and the average run length is at least
    ``threshold``, whatever their distribution and for any ``cal_len``.
    """

    def __init__(self, threshold: float = 28800.0, cal_len: int = 1440, config: RankBetConfig = RankBetConfig(),
                 seed: int = 0, recent: int = 240):
        self.A = float(threshold)
        self.cal_len = int(cal_len)
        self.cfg = config
        self.rng = np.random.default_rng(seed)
        self.recent = recent
        self.t = -1
        self.events: list[RankEvent] = []
        self.segments: list[dict] = []
        self.core = RankBetCore(config, 1)
        self._reset()

    def _reset(self):
        self.ranker = SequentialRanker(self.rng)
        self.n_ref = 0
        self.calibrated = False
        self.cal_start = self.t + 1
        self.core.reset()
        self._recent_p: list[float] = []

    def update(self, u: float) -> RankEvent | None:
        self.t += 1
        if not np.isfinite(u):
            return None
        s = float(score(u, self.cfg.kind))
        if not self.calibrated:
            self.ranker.add(s)
            self.n_ref += 1
            if self.n_ref >= self.cal_len:
                self.calibrated = True
                self.start = self.t + 1
                self.segments.append(dict(cal_start=self.cal_start, start=self.start))
            return None
        p = self.ranker.pvalue(s)
        M = float(self.core.step(np.array([p]))[0])
        self._recent_p.append(p)
        if len(self._recent_p) > self.recent:
            self._recent_p.pop(0)
        if M >= self.A:
            fam = self.core.family_stats()[0]
            sh = fam / fam.sum()
            ev = RankEvent(self.t, self.start, self.cal_start, M, float(sh[0]), float(sh[1]), float(sh[2]),
                           float(np.mean(self._recent_p)), tuple(float(v) for v in self.core.density_ratio(0)[0]))
            self.events.append(ev)
            self._reset()
            return ev
        return None

    def run(self, returns) -> list[RankEvent]:
        for r in np.asarray(returns, float):
            self.update(r)
        return self.events
