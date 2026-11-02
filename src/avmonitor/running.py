"""Anytime-valid volatility monitoring with a predictable, running truncation level.

Why this module exists
----------------------
The e-detectors of ``detectors.py`` bet on ``X_t = min(r_t^2, tau^2)`` against
the null mean ``m(tau) = E[min(r^2, tau^2)]``. Their guarantee
(P(false alarm within m steps) <= m / A) needs the TRUE ``m(tau)`` for ONE
fixed ``tau``. ``confseq.py`` (E8) keeps the guarantee with an estimated mean,
but only for a ``tau`` that is fixed before the confidence sequence starts.
Real data want a truncation level that follows the current scale. This module
gives two ways to let ``tau_t`` change at every step (any rule that uses only
past data, a *predictable* rule) and keep a proved guarantee. Full statements
and proofs are in ``docs/theory_running_tau.md``.

Route A (uniform band). If the calibration returns are i.i.d. with law F,
a band that holds for ALL truncation levels at once,
``L(s) < m_F(s) < U(s)`` for every s with probability at least 1 - delta,
lets the increase side bet against ``U(tau_t)`` and the decrease side against
``L(tau_t)`` for ANY predictable ``tau_t``. If each monitored return is, given
the past, stochastically no larger than F (increase side) or no smaller
(decrease side), P(false alarm within m steps) <= delta + m / A. Two bands:

* ``"dkw"``: the time-uniform Dvoretzky-Kiefer-Wolfowitz band of Howard and
  Ramdas (2022, Bernoulli 28(3), Theorem 2) for the distribution function of
  ``r^2``, integrated: ``m_F(s) = int_0^{s^2} P(r^2 > x) dx``.
* ``"grid"``: one-sided empirical Bernstein confidence sequences
  (``confseq.eb_confseq``, docs/math.md Theorem 6.2) for ``m_F(s_k)`` on a
  geometric grid of K levels fixed by a short pilot, each at ``delta / (2K)``
  (union bound). ``tau_t`` is snapped to the grid. Much tighter than ``"dkw"``.

Route B (ranks, conformal). If the returns are exchangeable (for example
i.i.d. with ANY law, ties allowed), the smoothed conformal p-values of
``|r_t|`` among all past values are i.i.d. uniform on [0, 1] (Vovk,
Nouretdinov and Gammerman 2003; Vovk 2021). Any bet ``f(p_t)`` with
``int_0^1 f = 1`` then has conditional mean exactly 1, and the
Shiryaev-Roberts recursion gives P(false alarm within m steps) <= m / A
EXACTLY: no baseline, no delta, no truncation level. The price: post-change
returns enter the reference set and dilute the evidence.

Everything here is new code; nothing in ``detectors.py``, ``monitor.py`` or
``confseq.py`` is changed.
"""

from __future__ import annotations

import bisect
import math
from collections import deque
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize, special, stats

from .confseq import eb_confseq, psi_coef
from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from .monitor import Event
from .robust import MAD_NORMAL_SCALE

# ----------------------------------------------------------------------------
# Defaults (fixed before any experiment; see docs/theory_running_tau.md, section 7)
# ----------------------------------------------------------------------------

GRID_RATIO = 2.0 ** 0.25  # neighbouring truncation levels differ by about 19%
GRID_KMIN, GRID_KMAX = -6, 6  # levels from 2^-1.5 to 2^1.5 times the centre (13 levels)
N_PILOT = 50  # returns used only to centre the grid of route A ("grid")
RUN_WINDOW = 120  # running scale: scaled MAD of the last 120 finite returns
RUN_EVERY = 10  # ... recomputed every 10 steps (piecewise constant, predictable)
LR_UP_RATIOS = (1.25, 1.5, 2.0)
LR_DOWN_RATIOS = (0.8, 1 / 1.5, 0.5)
LR_DOWN_LAMBDA = 0.95  # mix the decrease bets with 1 so that one jump cannot zero them
IND_UP = ((0.02, 0.1), (0.02, 0.25), (0.1, 0.1), (0.1, 0.25), (0.1, 0.5))  # (alpha, lambda)
IND_DOWN = ((0.25, 0.1), (0.25, 0.25), (0.5, 0.1), (0.5, 0.25), (0.5, 0.5))  # (beta, lambda)


# ----------------------------------------------------------------------------
# 1. The time-uniform DKW band of Howard and Ramdas (2022)
# ----------------------------------------------------------------------------


def _hr_error(C: float, A: float) -> float:
    """Smallest crossing probability that Theorem 2 of Howard and Ramdas certifies for constants (A, C).

    Mirrors ``EmpiricalProcessLILBound::find_optimal_C`` of the authors' library
    (github.com/gostevehoward/confseq, ``src/confseq/uniform_boundaries.h``):
    ``alpha(C) = min_eta 4 exp(-gamma^2 C) (1 + 1 / ((gamma^2 - 1) log eta))`` with
    ``gamma^2 = (2 / eta) (A - sqrt(2 (eta - 1) / C))^2``, over ``eta > 1`` with ``gamma^2 > 1``.
    """

    def err(eta):
        g2 = 2.0 / eta * (A - math.sqrt(2.0 * (eta - 1.0) / C)) ** 2
        if g2 <= 1.0:
            return math.inf
        return 4.0 * math.exp(-g2 * C) * (1.0 + 1.0 / ((g2 - 1.0) * math.log(eta)))

    # eta_upper solves sqrt(eta / 2) + sqrt(2 (eta - 1) / C) = A (beyond it gamma^2 <= 1)
    f = lambda eta: math.sqrt(eta / 2.0) + math.sqrt(2.0 * (eta - 1.0) / C) - A
    eta_up = optimize.brentq(f, 1.0, 2.0 * A * A)
    res = optimize.minimize_scalar(err, bounds=(1.0 + 1e-12, eta_up), method="bounded",
                                   options={"xatol": 1e-10})
    return float(res.fun)


def hr_dkw_constant(alpha: float, A: float = 0.85) -> float:
    """The constant C with ``_hr_error(C, A) = alpha`` (solved numerically, as in the authors' code)."""
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0, 1)")
    if not A > 1 / math.sqrt(2):
        raise ValueError("A must exceed 1 / sqrt(2)")
    lo, hi = 0.5, 5.0
    while _hr_error(hi, A) > alpha:
        hi *= 2.0
    while _hr_error(lo, A) < alpha and lo > 1e-6:
        lo /= 2.0
    return float(optimize.brentq(lambda C: _hr_error(C, A) - alpha, lo, hi, xtol=1e-10))


def hr_dkw_radius(t, alpha: float, t_min: float = 1.0, A: float = 0.85):
    """Radius ``eps_t = A sqrt((log(1 + log(t / t_min)) + C) / t)`` for ``t >= t_min`` (infinite before).

    Howard and Ramdas (2022), Theorem 2: for i.i.d. observations with ANY distribution
    function F, ``P(exists t >= t_min: sup_x |F_t(x) - F(x)| > eps_t) <= alpha``
    (``F_t`` = empirical distribution function after t observations). Their closed-form
    choice is ``A = 0.85`` with ``C`` about ``0.8 log(1612 / alpha)``.
    """
    C = hr_dkw_constant(alpha, A)
    t = np.asarray(t, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        eps = A * np.sqrt((np.log1p(np.log(np.maximum(t, t_min) / t_min)) + C) / t)
    return np.where(t >= t_min, eps, np.inf)


def massart_dkw_radius(n: int, alpha: float) -> float:
    """Fixed-n two-sided DKW radius with Massart's constant: ``sqrt(log(2 / alpha) / (2 n))``."""
    return math.sqrt(math.log(2.0 / alpha) / (2.0 * n))


def truncated_means_from_band(v_sorted: np.ndarray, levels2: np.ndarray, eps: float):
    """Bounds on ``m(s) = int_0^{s^2} P(V > x) dx`` from a band ``|F_n - F| <= eps`` (one path).

    ``v_sorted``: sorted sample of ``V = r^2`` (length n); ``levels2``: values ``s^2``.
    Returns ``(upper, lower, plug_in)`` at each level:

    * upper ``U(a) = int_0^a min(1, S_n(x) + eps) dx``,
    * lower ``L(a) = int_0^a max(0, S_n(x) - eps) dx``,
    * plug-in ``m_n(a) = int_0^a S_n(x) dx = mean(min(V, a))``,

    where ``S_n(x) = #{V_i > x} / n`` is the empirical survival function.
    """
    v = np.asarray(v_sorted, float)
    a = np.asarray(levels2, float)
    n = v.size
    csum = np.concatenate([[0.0], np.cumsum(v)])

    def mhat(b):
        k = np.searchsorted(v, b, side="right")  # #{V <= b}
        return (csum[k] + b * (n - k)) / n

    m_a = mhat(a)
    # upper: S_n(x) + eps >= 1 while #{V <= x} <= n eps, i.e. for x < v_(J + 1), J = floor(n eps)
    J = int(math.floor(n * eps + 1e-12))
    if J >= n:
        upper = a.copy()
    else:
        b = v[J]  # v_(J + 1) in 1-based order statistics
        upper = np.where(a <= b, a, b + m_a - mhat(b) + eps * (a - b))
    # lower: S_n(x) - eps > 0 while #{V <= x} < n (1 - eps), i.e. for x < v_(J' + 1)
    Jp = int(math.ceil(n * (1.0 - eps) - 1e-12)) - 1
    if Jp < 0:
        lower = np.zeros_like(a)
    else:
        b = v[min(Jp, n - 1)]
        ap = np.minimum(a, b)
        lower = np.maximum(mhat(ap) - eps * ap, 0.0)
    return upper, lower, m_a


# ----------------------------------------------------------------------------
# 2. Bands on m(s) at a grid of truncation levels, for many calibration windows
# ----------------------------------------------------------------------------


@dataclass
class TruncBand:
    """Bounds on ``m_F(s_k)`` at levels ``s_k`` for each path. Arrays have shape (reps, K)."""

    levels: np.ndarray  # truncation levels s_k (units of r)
    upper: np.ndarray  # U_k: with probability >= 1 - delta, m_F(s_k) < U_k for all k (route A)
    lower: np.ndarray  # L_k
    plug_in: np.ndarray  # empirical mean of min(r^2, s_k^2) over the whole window (no guarantee)
    center: np.ndarray  # c x scale used to centre the grid, shape (reps,)
    method: str
    delta: float
    info: dict = field(default_factory=dict)


def scaled_mad(r, axis=-1):
    """Scaled median absolute value about zero, ignoring NaN (the robust sd used throughout the package)."""
    return np.nanmedian(np.abs(r), axis=axis) * MAD_NORMAL_SCALE


def grid_levels(center, ratio: float = GRID_RATIO, kmin: int = GRID_KMIN, kmax: int = GRID_KMAX):
    """Geometric grid ``center * ratio^k`` for k = kmin..kmax, shape (reps, K)."""
    k = np.arange(kmin, kmax + 1)
    return np.asarray(center, float)[..., None] * ratio ** k


def calibrate_band(
    r_cal,
    method: str = "grid",
    c: float = 3.0,
    delta: float = 0.01,
    ratio: float = GRID_RATIO,
    kmin: int = GRID_KMIN,
    kmax: int = GRID_KMAX,
    n_pilot: int = N_PILOT,
    cs_c: float = 0.75,
    time_uniform: bool = True,
    center=None,
    chunk: int = 100,
) -> TruncBand:
    """Route A bands from calibration windows ``r_cal`` (rows = paths, no NaN).

    ``delta`` is the TOTAL failure probability of the band (both sides, all levels).

    * ``method="grid"``: the first ``n_pilot`` returns give the centre ``c * scaled MAD``;
      on the remaining returns, for each of the K grid levels, an upper and a lower
      empirical Bernstein confidence sequence (``confseq.eb_confseq``, tuned for the
      window length) at level ``delta / (2 K)`` each, read at the end of the window.
    * ``method="dkw"``: the centre is ``c * scaled MAD`` of the whole window; the band is
      ``|F_n - F| <= eps`` for the distribution function of ``r^2`` with the
      Howard-Ramdas radius at ``t = n`` (``time_uniform=True``; ``t_min = n``) or the
      fixed-n Massart radius (``time_uniform=False``, valid only for a window length
      fixed in advance), and the bounds are its integrals (``truncated_means_from_band``).
    * ``method="plugin"``: no bound; ``upper = lower = plug_in`` (for comparison only).

    ``center`` (shape (reps,)) overrides the grid centre. For ``"grid"`` it must be computed
    from the pilot returns only, so that the levels are fixed before the sequences start.
    """
    r_cal = np.atleast_2d(np.asarray(r_cal, float))
    if np.any(~np.isfinite(r_cal)):
        raise ValueError("drop missing values before calibrating")
    reps, n = r_cal.shape
    x2 = r_cal**2
    info = {}
    if method == "grid":
        if not 1 <= n_pilot < n:
            raise ValueError("n_pilot must leave returns for the confidence sequences")
        center = c * scaled_mad(r_cal[:, :n_pilot]) if center is None else np.asarray(center, float).reshape(-1)
        if np.any(center <= 0):
            raise ValueError("pilot scale is zero on some rows; use a longer pilot")
        levels = grid_levels(center, ratio, kmin, kmax)
        K = levels.shape[1]
        d_each = delta / (2 * K)
        upper = np.empty((reps, K))
        lower = np.empty((reps, K))
        rest = x2[:, n_pilot:]
        m = rest.shape[1]
        for a in range(0, reps, chunk):
            sl = slice(a, a + chunk)
            lv2 = levels[sl] ** 2  # (b, K)
            b = lv2.shape[0]
            xx = np.minimum(rest[sl][:, None, :], lv2[:, :, None]).reshape(b * K, m)
            cs = eb_confseq(xx, lv2.reshape(-1), d_each, d_each, n_plan=m, c=cs_c)
            upper[sl] = cs.upper[:, -1].reshape(b, K)
            lower[sl] = cs.lower[:, -1].reshape(b, K)
        info.update(n_pilot=n_pilot, delta_each=d_each, n_cs=m)
    elif method in ("dkw", "plugin"):
        center = c * scaled_mad(r_cal) if center is None else np.asarray(center, float).reshape(-1)
        levels = grid_levels(center, ratio, kmin, kmax)
        if method == "dkw":
            eps = float(hr_dkw_radius(n, delta, t_min=n)) if time_uniform else massart_dkw_radius(n, delta)
            info.update(eps=eps, time_uniform=time_uniform)
        upper = np.empty_like(levels)
        lower = np.empty_like(levels)
    else:
        raise ValueError("method must be 'grid', 'dkw' or 'plugin'")
    plug = np.empty_like(levels)
    srt = np.sort(x2, axis=1)
    for i in range(reps):
        if method == "dkw":
            upper[i], lower[i], plug[i] = truncated_means_from_band(srt[i], levels[i] ** 2, eps)
        else:
            plug[i] = truncated_means_from_band(srt[i], levels[i] ** 2, 0.0)[2]
    if method == "plugin":
        upper, lower = plug.copy(), plug.copy()
    return TruncBand(levels=levels, upper=upper, lower=lower, plug_in=plug, center=center, method=method,
                     delta=delta if method != "plugin" else 0.0, info=info)


# ----------------------------------------------------------------------------
# 3. A predictable running truncation level, snapped to the grid
# ----------------------------------------------------------------------------


def running_level_index(r_hist, start: int, n_mon: int, center, c: float = 3.0, window: int = RUN_WINDOW,
                        every: int = RUN_EVERY, ratio: float = GRID_RATIO, kmin: int = GRID_KMIN,
                        kmax: int = GRID_KMAX, fixed: bool = False):
    """Grid index of ``tau_t`` for monitored steps ``t = start, ..., start + n_mon - 1`` (shape (reps, n_mon)).

    ``tau_t`` = the grid level nearest (in log scale) to ``c * scaled MAD`` of the
    ``window`` returns before the last refresh time, refreshed every ``every`` steps:
    the value used at step t depends only on returns with index < t (predictable).
    ``r_hist`` holds calibration and monitored returns of each path, (reps, start + n_mon).
    ``fixed=True`` returns the index of the centre (level 0) at every step.
    """
    r_hist = np.atleast_2d(np.asarray(r_hist, float))
    reps = r_hist.shape[0]
    K0 = -kmin
    if fixed:
        return np.full((reps, n_mon), K0, dtype=np.int64)
    out = np.empty((reps, n_mon), dtype=np.int64)
    center = np.asarray(center, float)
    for j in range(0, n_mon, every):
        t = start + j
        seg = r_hist[:, max(0, t - window) : t]
        s = c * scaled_mad(seg)
        with np.errstate(divide="ignore", invalid="ignore"):
            k = np.rint(np.log(s / center) / np.log(ratio))
        k = np.where(np.isfinite(k), k, 0)
        out[:, j : j + every] = (np.clip(k, kmin, kmax) + K0).astype(np.int64)[:, None]
    return out


# ----------------------------------------------------------------------------
# 4. Shiryaev-Roberts recursions with time-varying truncation (vectorised)
# ----------------------------------------------------------------------------


def two_sided_sr(x, tau2, m_up, m_dn, levels, up_grid=DEFAULT_UP_GRID, down_fracs=DEFAULT_DOWN_FRACTIONS,
                 kind: str = "SR"):
    """First alarms of the increase and the decrease statistic on ``x_t = min(r_t^2, tau_t^2)``.

    ``x, tau2, m_up, m_dn``: arrays (reps, n) (or broadcastable), all predictable except x.
    Increase factors ``1 + lam (x / m_up - 1)`` with lam in ``up_grid``; decrease factors
    ``1 + lam_t (1 - x / m_dn)`` with ``lam_t = frac / (tau2 / m_dn - 1)`` recomputed every step
    (a step with ``m_dn <= 0`` or ``tau2 <= m_dn`` gets factor 1). Statistics are averages
    over the grid. Returns dict level -> (first up alarm, first down alarm), -1 if none.
    """
    x = np.atleast_2d(np.asarray(x, float))
    reps, n = x.shape
    tau2 = np.broadcast_to(np.asarray(tau2, float), (reps, n))
    m_up = np.broadcast_to(np.asarray(m_up, float), (reps, n))
    m_dn = np.broadcast_to(np.asarray(m_dn, float), (reps, n))
    ug = np.asarray(up_grid, float)[None, :]
    df = np.asarray(down_fracs, float)[None, :]
    Su = np.zeros((reps, ug.shape[1]))
    Sd = np.zeros((reps, df.shape[1]))
    levels = sorted(levels)
    out = {L: [np.full(reps, -1), np.full(reps, -1)] for L in levels}
    top = levels[-1]
    for t in range(n):
        xt = x[:, t]
        yu = xt / m_up[:, t]
        on = (m_dn[:, t] > 0) & (tau2[:, t] > m_dn[:, t])
        md = np.where(on, m_dn[:, t], 1.0)
        B = np.where(on, tau2[:, t] / md, 2.0)
        lam = np.where(on[:, None], df / (B[:, None] - 1.0), 0.0)
        fu = 1.0 + ug * (yu[:, None] - 1.0)
        fd = 1.0 + lam * (1.0 - (xt / md)[:, None])
        with np.errstate(over="ignore"):
            if kind == "SR":
                Su = fu * (Su + 1.0)
                Sd = fd * (Sd + 1.0)
            else:
                Su = fu * np.maximum(Su, 1.0)
                Sd = fd * np.maximum(Sd, 1.0)
        Mu = Su.mean(axis=1)
        Md = Sd.mean(axis=1)
        for L in levels:
            au, ad = out[L]
            au[(au < 0) & (Mu >= L)] = t
            ad[(ad < 0) & (Md >= L)] = t
        if t % 64 == 0 and (out[top][0] >= 0).all() and (out[top][1] >= 0).all():
            break
    return out


def route_a_inputs(band: TruncBand, r_hist, start: int, n_mon: int, c: float = 3.0, null_band: float = 0.0,
                   fixed: bool = False, **kw):
    """Per-step ``tau_t^2``, ``m_up_t``, ``m_dn_t`` (each (reps, n_mon)) for route A (or the plug-in band)."""
    k = running_level_index(r_hist, start, n_mon, band.center, c=c, fixed=fixed, **kw)
    rows = np.arange(k.shape[0])[:, None]
    tau2 = band.levels[rows, k] ** 2
    m_up = band.upper[rows, k] * (1.0 + null_band)
    m_dn = band.lower[rows, k] * (1.0 - null_band)
    return tau2, m_up, m_dn, k


# ----------------------------------------------------------------------------
# 5. Route B: smoothed conformal p-values and rank bets
# ----------------------------------------------------------------------------


def dense_ranks(a: np.ndarray) -> np.ndarray:
    """1-based dense ranks of each row (equal values share a rank)."""
    a = np.atleast_2d(a)
    order = np.argsort(a, axis=1, kind="stable")
    srt = np.take_along_axis(a, order, axis=1)
    new = np.concatenate([np.ones((a.shape[0], 1), bool), srt[:, 1:] != srt[:, :-1]], axis=1)
    dr = np.cumsum(new, axis=1)
    out = np.empty_like(dr)
    np.put_along_axis(out, order, dr, axis=1)
    return out


class _RowFenwick:
    """One Fenwick (binary indexed) tree per row, for counts of values with rank <= q."""

    def __init__(self, reps: int, size: int):
        self.n = size
        self.tree = np.zeros((reps, size + 1), dtype=np.int64)
        self.rows = np.arange(reps)

    def add(self, q: np.ndarray, mask: np.ndarray | None = None):
        q = q.copy()
        rows = self.rows
        if mask is not None:
            rows, q = rows[mask], q[mask]
        while rows.size:
            self.tree[rows, q] += 1
            q = q + (q & -q)
            keep = q <= self.n
            rows, q = rows[keep], q[keep]

    def count_le(self, q: np.ndarray) -> np.ndarray:
        s = np.zeros(self.rows.size, dtype=np.int64)
        q = q.copy()
        rows = self.rows
        while rows.size:
            s[rows] += self.tree[rows, q]
            q = q - (q & -q)
            keep = q > 0
            rows, q = rows[keep], q[keep]
        return s


def conformal_pvalues(a, rng: np.random.Generator) -> np.ndarray:
    """Smoothed conformal p-values ``p_t = (#{i <= t: a_i > a_t} + theta_t #{i <= t: a_i = a_t}) / t``.

    Rows are independent paths; ``a`` is the nonconformity score (here ``|r|``). Under
    exchangeability of each row the p-values of a row are i.i.d. uniform on [0, 1]
    (docs/theory_running_tau.md, Proposition B.1). The count of equal values includes
    ``a_t`` itself, and ``theta_t`` are i.i.d. uniform tie-breakers.
    """
    a = np.atleast_2d(np.asarray(a, float))
    reps, n = a.shape
    rk = dense_ranks(a)
    fw = _RowFenwick(reps, int(rk.max()))
    theta = rng.random((reps, n))
    p = np.empty((reps, n))
    for t in range(n):
        q = rk[:, t]
        le = fw.count_le(q)
        lt = fw.count_le(q - 1)
        greater = t - le
        equal = le - lt + 1
        p[:, t] = (greater + theta[:, t] * equal) / (t + 1.0)
        fw.add(q)
    return p


def _abs_normal_quantile(u):
    """``q(u)``: the value with ``P(|Z| > q) = u`` for a standard normal Z (``-Phi^-1(u / 2)``)."""
    return -special.ndtri(np.asarray(u, float) / 2.0)


def lr_bet_factory(rho: float, u0: float = 2 * stats.norm.sf(3.0), lam: float = 1.0):
    """Betting function on a p-value: the Gaussian scale likelihood ratio in rank space, clipped.

    If ``|r|`` has the reference law and the new return is ``rho`` times a draw from it,
    and the reference law is normal, the p-value has density
    ``g(u) = (1 / rho) exp(-q(u)^2 (rho^-2 - 1) / 2)`` with ``q(u) = Phi^-1(1 - u / 2)``.
    The bet is ``f(u) = 1 + lam (g(max(u, u0)) / C - 1)`` with C chosen so that
    ``int_0^1 f = 1``. ``u0 = P(|Z| > 3)`` is the rank analogue of truncating at
    3 standard deviations: one extreme return multiplies the wealth by at most
    ``1 + lam (g(u0) / C - 1)``. Validity never depends on normality.
    """
    q0 = float(_abs_normal_quantile(u0))
    g0 = math.exp(-q0 * q0 * (rho**-2 - 1.0) / 2.0) / rho
    G0 = 2.0 * stats.norm.sf(q0 / rho)  # P_alt(p <= u0)
    C = u0 * g0 + (1.0 - G0)

    def f(u):
        u = np.maximum(np.asarray(u, float), u0)
        q = _abs_normal_quantile(u)
        g = np.exp(-q * q * (rho**-2 - 1.0) / 2.0) / rho
        return 1.0 + lam * (g / C - 1.0)

    f.rho, f.u0, f.lam, f.norm_const = rho, u0, lam, C
    return f


def indicator_bet_factory(level: float, lam: float, side: str):
    """``1 + lam (1{p <= level} / level - 1)`` (side "up") or ``1 + lam (1{p > 1 - level} / level - 1)`` ("down").

    "up" bets that ``|r_t|`` exceeds the running (1 - level) quantile of the past values,
    "down" that it falls below the running ``level`` quantile. Mean exactly 1 under a uniform p.
    """
    if not (0 < level < 1 and 0 <= lam <= 1):
        raise ValueError("need 0 < level < 1 and 0 <= lam <= 1")

    def f(u):
        u = np.asarray(u, float)
        hit = u <= level if side == "up" else u > 1.0 - level
        return 1.0 + lam * (hit / level - 1.0)

    return f


def bet_family(name: str):
    """Lists of betting functions (increase side, decrease side) for ``name`` in {"lr", "ind"}."""
    if name == "lr":
        up = [lr_bet_factory(r) for r in LR_UP_RATIOS]
        dn = [lr_bet_factory(r, lam=LR_DOWN_LAMBDA) for r in LR_DOWN_RATIOS]
    elif name == "ind":
        up = [indicator_bet_factory(a, l, "up") for a, l in IND_UP]
        dn = [indicator_bet_factory(b, l, "down") for b, l in IND_DOWN]
    else:
        raise ValueError("family must be 'lr' or 'ind'")
    return up, dn


def family_factors(name: str):
    """Vectorised evaluator ``p -> (factors of the increase bets, factors of the decrease bets)``.

    Same functions as ``bet_family(name)``, computed in one pass (for speed); p may be a
    scalar or a 1-d array, the output arrays have shape (len(p), K).
    """
    up, dn = bet_family(name)
    if name == "ind":
        lu = np.array([a for a, _ in IND_UP]); wu = np.array([l for _, l in IND_UP])
        ld = np.array([b for b, _ in IND_DOWN]); wd = np.array([l for _, l in IND_DOWN])

        def ev(p):
            p = np.atleast_1d(np.asarray(p, float))[:, None]
            return 1.0 + wu * ((p <= lu) / lu - 1.0), 1.0 + wd * ((p > 1.0 - ld) / ld - 1.0)

        return ev
    u0 = up[0].u0
    ku = np.array([f.rho**-2 - 1.0 for f in up]); cu = np.array([1.0 / (f.rho * f.norm_const) for f in up])
    lu = np.array([f.lam for f in up])
    kd = np.array([f.rho**-2 - 1.0 for f in dn]); cd = np.array([1.0 / (f.rho * f.norm_const) for f in dn])
    ld = np.array([f.lam for f in dn])

    def ev(p):
        p = np.atleast_1d(np.asarray(p, float))
        q2 = (-special.ndtri(np.maximum(p, u0) / 2.0))[:, None] ** 2
        return 1.0 + lu * (cu * np.exp(-q2 * ku / 2.0) - 1.0), 1.0 + ld * (cd * np.exp(-q2 * kd / 2.0) - 1.0)

    return ev


def conformal_detector(a, rng: np.random.Generator, warmup: int, levels, family: str = "lr",
                       freeze_at: float | None = None, kind: str = "SR", band: float = 0.0):
    """Route B on rows of ``a`` (scores ``|r|``): first alarms of the increase and the decrease statistic.

    Steps ``0 .. warmup - 1`` only enter the reference set (no bets). From ``warmup`` on,
    p-values are computed against all earlier scores of the row and both statistics are
    updated (Shiryaev-Roberts or CUSUM, averaged over the family). ``band = b > 0``
    divides each bet's payoff ``f(p)`` by ``1 + b`` before betting (``1 + (f - 1 - b) / (1 + b)``),
    a tolerance for mild departures from exchangeability; it never hurts validity.

    ``freeze_at = A0`` (heuristic, NOT covered by the exact guarantee): once either
    statistic of a row reaches A0, that row's reference set stops growing for good
    (later p-values use the frozen set). Returns (dict level -> (up alarms, down alarms),
    freeze times or None).
    """
    a = np.atleast_2d(np.asarray(a, float))
    reps, n = a.shape
    up, dn = bet_family(family)
    fam = family_factors(family)
    rk = dense_ranks(a)
    fw = _RowFenwick(reps, int(rk.max()))
    theta = rng.random((reps, n))
    Su = np.zeros((reps, len(up)))
    Sd = np.zeros((reps, len(dn)))
    levels = sorted(levels)
    out = {L: [np.full(reps, -1), np.full(reps, -1)] for L in levels}
    top = levels[-1]
    size = np.zeros(reps, dtype=np.int64)  # size of each row's reference set
    frozen = np.zeros(reps, bool)
    freeze_time = np.full(reps, -1) if freeze_at is not None else None
    for t in range(n):
        q = rk[:, t]
        if t >= warmup:
            le = fw.count_le(q)
            lt = fw.count_le(q - 1)
            p = (size - le + theta[:, t] * (le - lt + 1)) / (size + 1.0)
            fu, fd = fam(p)
            if band:
                fu = 1.0 + (fu - 1.0 - band) / (1.0 + band)
                fd = 1.0 + (fd - 1.0 - band) / (1.0 + band)
            with np.errstate(over="ignore"):
                if kind == "SR":
                    Su = fu * (Su + 1.0)
                    Sd = fd * (Sd + 1.0)
                else:
                    Su = fu * np.maximum(Su, 1.0)
                    Sd = fd * np.maximum(Sd, 1.0)
            Mu, Md = Su.mean(axis=1), Sd.mean(axis=1)
            for L in levels:
                au, ad = out[L]
                au[(au < 0) & (Mu >= L)] = t
                ad[(ad < 0) & (Md >= L)] = t
            if freeze_at is not None:
                new = ~frozen & (np.maximum(Mu, Md) >= freeze_at)
                freeze_time[new] = t
                frozen |= new
        grow = ~frozen
        fw.add(q, grow if frozen.any() else None)
        size += grow
        if t >= warmup and t % 64 == 0 and (out[top][0] >= 0).all() and (out[top][1] >= 0).all():
            break
    return out, freeze_time


# ----------------------------------------------------------------------------
# 6. Streaming monitors for real data (same protocol as monitor.VolatilityMonitor)
# ----------------------------------------------------------------------------


class RunningTauMonitor:
    """Route A in streaming form: calibrate, monitor with a running ``tau_t``, alarm, recalibrate.

    Calibration collects ``cal_len`` finite returns and builds a band (``calibrate_band``)
    with total failure probability ``delta``. Monitoring keeps the band fixed and sets
    ``tau_t`` to the grid level nearest ``c * scaled MAD`` of the last ``window`` finite
    returns (refreshed every ``every`` monitored steps; the first value comes from the end of
    the calibration window). The increase side bets against ``(1 + b) U(tau_t)``, the
    decrease side against ``(1 - b) L(tau_t)``; each side alarms at ``2 A``. Missing returns
    are skipped. ``running=False`` keeps ``tau`` at the grid centre (route A with a fixed level).
    If the calibration window has scale 0 (more than half zero returns), it is extended
    until usable (the extension rule uses only past data, so the stopping-time argument of
    the theory applies to the time-uniform bands).
    """

    def __init__(self, threshold: float = 36000.0, c: float = 3.0, cal_len: int = 1440, method: str = "grid",
                 delta: float = 0.01, null_band: float = 0.0, running: bool = True, window: int = RUN_WINDOW,
                 every: int = RUN_EVERY, up_grid=DEFAULT_UP_GRID, down_fractions=DEFAULT_DOWN_FRACTIONS,
                 n_pilot: int = N_PILOT):
        self.A = float(threshold)
        self.c = c
        self.cal_len = int(cal_len)
        self.method = method
        self.delta = delta
        self.band = null_band
        self.running = running
        self.window = window
        self.every = every
        self.up_grid = np.asarray(up_grid, float)
        self.down_fractions = np.asarray(down_fractions, float)
        self.n_pilot = n_pilot
        self.t = -1
        self.events: list[Event] = []
        self.segments: list[dict] = []
        self._reset()

    def _reset(self):
        self.buffer: list[float] = []
        self.calibrated = False
        self.S_up = np.zeros(len(self.up_grid))
        self.S_dn = np.zeros(len(self.down_fractions))
        self.recent: deque = deque(maxlen=self.window)
        self.k = None
        self.n_mon = 0

    def _calibrate(self):
        """Build the band. The centre is c x scaled MAD of the pilot ("grid") or of the window; if
        that is 0 (more than half zero returns), c x the scaled median of the nonzero |r| of the
        same returns; if all are zero, calibration restarts with fresh data. Every rule uses
        only returns already seen, so the end of calibration is a stopping time."""
        r = np.asarray(self.buffer)
        piece = r[: self.n_pilot] if self.method == "grid" else r
        s = scaled_mad(piece)
        if s <= 0:
            nz = np.abs(piece[piece != 0])
            s = float(np.median(nz) * MAD_NORMAL_SCALE) if nz.size else 0.0
        if s <= 0:
            self.buffer = []
            return
        b = calibrate_band(r[None, :], method=self.method, c=self.c, delta=self.delta, n_pilot=self.n_pilot,
                           center=np.array([self.c * s]))
        self.scale = float(scaled_mad(r))
        self.levels = b.levels[0]
        self.U = b.upper[0] * (1.0 + self.band)
        self.L = b.lower[0] * (1.0 - self.band)
        self.center = float(b.center[0])
        self.K0 = -GRID_KMIN
        self.recent.extend(r[-self.window :])
        self._set_level()
        self.start = self.t + 1
        self.calibrated = True
        self.segments.append(dict(start=self.start, cal_len=len(self.buffer), center=self.center,
                                  U0=float(self.U[self.K0]), L0=float(self.L[self.K0])))

    def _set_level(self):
        if not self.running:
            self.k = self.K0
            return
        s = self.c * scaled_mad(np.asarray(self.recent))
        k = round(math.log(s / self.center) / math.log(GRID_RATIO)) if s > 0 else 0
        self.k = int(min(max(k, GRID_KMIN), GRID_KMAX)) + self.K0

    def update(self, r: float) -> Event | None:
        self.t += 1
        if not math.isfinite(r):
            return None
        if not self.calibrated:
            self.buffer.append(r)
            if len(self.buffer) >= self.cal_len:
                self._calibrate()
            return None
        if self.n_mon > 0 and self.n_mon % self.every == 0:
            self._set_level()
        k = self.k
        tau2 = self.levels[k] ** 2
        x = min(r * r, tau2)
        self.S_up = (1.0 + self.up_grid * (x / self.U[k] - 1.0)) * (self.S_up + 1.0)
        md = self.L[k]
        if md > 0 and tau2 > md:
            lam = self.down_fractions / (tau2 / md - 1.0)
            self.S_dn = (1.0 + lam * (1.0 - x / md)) * (self.S_dn + 1.0)
        else:
            self.S_dn = self.S_dn + 1.0
        self.recent.append(r)
        self.n_mon += 1
        mu, mdn = self.S_up.mean(), self.S_dn.mean()
        ev = None
        if mu >= 2 * self.A:
            ev = Event(self.t, "up", self.start, self.scale, float(mu))
        elif mdn >= 2 * self.A:
            ev = Event(self.t, "down", self.start, self.scale, float(mdn))
        if ev is not None:
            self.events.append(ev)
            self._reset()
        return ev

    def run(self, returns) -> list[Event]:
        for r in np.asarray(returns, float):
            self.update(r)
        return self.events


class ConformalMonitor:
    """Route B in streaming form: exact rank-based e-detector, restarted after each alarm.

    Each segment collects ``cal_len`` finite returns into the reference set without
    betting, then for every new finite return computes the smoothed conformal p-value of
    ``|r_t|`` among all finite returns of the segment so far (itself included), updates
    the Shiryaev-Roberts statistics of the increase and the decrease family (each at
    ``2 A``) and adds ``|r_t|`` to the reference set. Under exchangeability within each
    segment, P(false alarm within m monitored steps) <= m / A exactly. ``baseline_scale``
    in an event is the scaled MAD of the reference set at the alarm (the returns the alarm
    compares against; the warm-up returns if fewer than 30), for the E5c validation.

    ``restart_len`` (default ``cal_len``) is the warm-up after an alarm. ``restart_len=0``
    gives a monitor with no blind period: after an alarm it bets from the very next
    return, against a reference set that starts empty and grows; the guarantee is
    unchanged because the p-values are exact from the first one.
    """

    def __init__(self, threshold: float = 28800.0, cal_len: int = 1440, family: str = "lr", band: float = 0.0,
                 seed: int = 0, kind: str = "SR", restart_len: int | None = None):
        self.A = float(threshold)
        self.cal_len = int(cal_len)
        self.restart_len = self.cal_len if restart_len is None else int(restart_len)
        self.up, self.dn = bet_family(family)
        self.fam = family_factors(family)
        self.band = band
        self.kind = kind
        self.rng = np.random.default_rng(seed)
        self.t = -1
        self.events: list[Event] = []
        self.segments: list[dict] = []
        self._reset()

    def _reset(self, warm: int | None = None):
        self.ref: list[float] = []
        self.calibrated = False
        self.warm = self.cal_len if warm is None else warm
        self.S_up = np.zeros(len(self.up))
        self.S_dn = np.zeros(len(self.dn))
        if self.warm == 0:
            self._start()

    def _start(self):
        self.calibrated = True
        self.start = self.t + 1
        self.scale = float(np.median(self.ref) * MAD_NORMAL_SCALE) if self.ref else float("nan")
        self.segments.append(dict(start=self.start, scale=self.scale))

    def pvalue(self, a: float) -> float:
        n = len(self.ref)
        lo = bisect.bisect_left(self.ref, a)
        hi = bisect.bisect_right(self.ref, a)
        return ((n - hi) + self.rng.random() * (hi - lo + 1)) / (n + 1.0)

    def update(self, r: float) -> Event | None:
        self.t += 1
        if not math.isfinite(r):
            return None
        a = abs(r)
        if not self.calibrated:
            bisect.insort(self.ref, a)
            if len(self.ref) >= self.warm:
                self._start()
            return None
        p = self.pvalue(a)
        fu, fd = self.fam(p)
        fu, fd = fu[0], fd[0]
        if self.band:
            fu = 1.0 + (fu - 1.0 - self.band) / (1.0 + self.band)
            fd = 1.0 + (fd - 1.0 - self.band) / (1.0 + self.band)
        if self.kind == "SR":
            self.S_up = fu * (self.S_up + 1.0)
            self.S_dn = fd * (self.S_dn + 1.0)
        else:
            self.S_up = fu * np.maximum(self.S_up, 1.0)
            self.S_dn = fd * np.maximum(self.S_dn, 1.0)
        mu, md = self.S_up.mean(), self.S_dn.mean()
        ev = None
        if mu >= 2 * self.A or md >= 2 * self.A:
            scale = float(np.median(self.ref) * MAD_NORMAL_SCALE) if len(self.ref) >= 30 else self.scale
            if mu >= 2 * self.A:
                ev = Event(self.t, "up", self.start, scale, float(mu))
            else:
                ev = Event(self.t, "down", self.start, scale, float(md))
        bisect.insort(self.ref, a)
        if ev is not None:
            self.events.append(ev)
            self._reset(self.restart_len)
        return ev

    def run(self, returns) -> list[Event]:
        for r in np.asarray(returns, float):
            self.update(r)
        return self.events


# ----------------------------------------------------------------------------
# 7. Monitoring with no blind period (route A "awake", staggered bands)
# ----------------------------------------------------------------------------

AWAKE_PRE = 120  # finite returns before a (re)start that centre the grid of truncation levels
AWAKE_LEARN = 1440  # the band learns from the first 1440 monitored returns of a segment, then is frozen
AWAKE_N_PLAN = 300  # bet schedule of the learning confidence sequences (tightest near 300 observations)
AWAKE_FIRST_EPOCH = 16
AWAKE_EPOCH_RATIO = 2.0 ** 0.5


def awake_epochs(n_learn: int = AWAKE_LEARN, first: int = AWAKE_FIRST_EPOCH, ratio: float = AWAKE_EPOCH_RATIO):
    """Snapshot sizes ``0 = e_0 < e_1 < ... < e_last = n_learn`` (observations in the band at each snapshot).

    The chain started at segment step s (s = 0, 1, ...) bets against the band built from the
    first ``e(s)`` returns of the segment, ``e(s)`` = the largest snapshot size ``<= s``, so the
    band uses only returns before s. Geometric sizes keep the number of snapshots at
    O(log n_learn) while each chain uses at least ``s / ratio`` returns.
    """
    sizes = {0, int(n_learn)}
    e = float(first)
    while e < n_learn:
        sizes.add(int(round(e)))
        e *= ratio
    return np.array(sorted(sizes), dtype=np.int64)


class _LevelsCS:
    """``confseq.StreamingEBConfSeq`` for K truncation levels at once (vectorised over levels).

    Observation ``z_k = min(r^2, s_k^2) / s_k^2``; upper and lower sequences at level
    ``delta_each`` each. Reproduces ``eb_confseq`` exactly (tested).
    """

    def __init__(self, levels, delta_each: float, n_plan: float | None = AWAKE_N_PLAN, c: float = 0.75):
        self.s2 = np.asarray(levels, float) ** 2
        K = self.s2.size
        self.log_d = math.log(1.0 / delta_each)
        self.n_plan = n_plan
        self.c = c
        self.t = 0
        self.sum_z = np.zeros(K)
        self.sum_dev2 = np.zeros(K)
        self.zhat = np.full(K, 0.5)
        self.var = np.full(K, 0.25)
        self.acc = np.zeros((6, K))  # sl_u, slz_u, spe_u, sl_l, slz_l, spe_l

    def update(self, r: float):
        z = np.minimum(r * r, self.s2) / self.s2
        self.t += 1
        t = self.t
        zprev, vprev = self.zhat, self.var
        d = float(self.n_plan) if self.n_plan is not None else t * math.log1p(t)
        e = (z - zprev) ** 2
        base = np.sqrt(2.0 * self.log_d / (d * vprev))
        a_u = 1.0 - zprev
        lu = np.minimum(base, self.c / a_u)
        a_l = zprev
        ll = np.minimum(base, self.c / a_l)
        acc = self.acc
        acc[0] += lu
        acc[1] += lu * z
        acc[2] += lu * lu * psi_coef(lu * a_u) * e
        acc[3] += ll
        acc[4] += ll * z
        acc[5] += ll * ll * psi_coef(ll * a_l) * e
        self.sum_z += z
        self.zhat = (0.5 + self.sum_z) / (t + 1.0)
        self.sum_dev2 += (z - self.zhat) ** 2
        self.var = (0.25 + self.sum_dev2) / (t + 1.0)

    def bounds(self):
        """(upper, lower) for ``m(s_k)`` in units of r^2; ``(s_k^2, 0)`` before any observation."""
        if self.t == 0:
            return self.s2.copy(), np.zeros_like(self.s2)
        a = self.acc
        up = np.minimum((a[1] + self.log_d + a[2]) / a[0], 1.0) * self.s2
        lo = np.maximum((a[4] - self.log_d - a[5]) / a[3], 0.0) * self.s2
        return up, lo


def awake_band_snapshots(r_seg, center, epochs=None, delta: float = 0.01, n_plan=AWAKE_N_PLAN, cs_c: float = 0.75,
                         ratio: float = GRID_RATIO, kmin: int = GRID_KMIN, kmax: int = GRID_KMAX, chunk: int = 100):
    """Vectorised snapshots of the learning band for many segments.

    ``r_seg`` (reps, n): the segment's returns (no NaN); ``center`` (reps,): grid centre known
    before the segment starts. Returns ``levels`` (reps, K) and ``upper``, ``lower`` of shape
    (reps, E, K): the band after ``epochs[e]`` observations (``epochs[0] = 0`` gives the trivial
    band ``(s_k^2, 0)``). Each one-sided sequence has level ``delta / (2 K)``.
    """
    r_seg = np.atleast_2d(np.asarray(r_seg, float))
    reps, n = r_seg.shape
    epochs = awake_epochs() if epochs is None else np.asarray(epochs)
    epochs = epochs[epochs <= n]
    levels = grid_levels(center, ratio, kmin, kmax)
    K = levels.shape[1]
    d_each = delta / (2 * K)
    E = epochs.size
    upper = np.empty((reps, E, K))
    lower = np.empty((reps, E, K))
    upper[:, 0, :] = levels**2
    lower[:, 0, :] = 0.0
    m = int(epochs[-1])
    if m > 0:
        x2 = r_seg[:, :m] ** 2
        for a in range(0, reps, chunk):
            sl = slice(a, a + chunk)
            lv2 = levels[sl] ** 2
            b = lv2.shape[0]
            xx = np.minimum(x2[sl][:, None, :], lv2[:, :, None]).reshape(b * K, m)
            cs = eb_confseq(xx, lv2.reshape(-1), d_each, d_each, n_plan=n_plan, c=cs_c)
            idx = epochs[1:] - 1
            upper[sl, 1:, :] = cs.upper[:, idx].reshape(b, K, E - 1).transpose(0, 2, 1)
            lower[sl, 1:, :] = cs.lower[:, idx].reshape(b, K, E - 1).transpose(0, 2, 1)
    return levels, upper, lower, epochs


def awake_route_a(r_pre, r_seg, levels_list, delta: float = 0.01, c: float = 3.0, null_band: float = 0.0,
                  running: bool = True, up_grid=DEFAULT_UP_GRID, down_fracs=DEFAULT_DOWN_FRACTIONS,
                  epochs=None, n_plan=AWAKE_N_PLAN, window: int = RUN_WINDOW, every: int = RUN_EVERY):
    """Route A with no blind period, vectorised: first alarms of both sides from segment step 0 on.

    ``r_pre`` (reps, n_pre): returns before the segment start (they centre the grid and seed
    the running window; they are NOT part of the band). ``r_seg`` (reps, n): the segment.
    The chain started at step s bets against the band snapshot of size ``e(s)``
    (``awake_epochs``); chains with the same snapshot share one Shiryaev-Roberts statistic,
    and the alarm statistic is the sum over snapshots (docs/theory_running_tau.md, Theorem A.4).
    Returns (dict level -> (up alarms, down alarms), centre).
    """
    r_pre = np.atleast_2d(np.asarray(r_pre, float))
    r_seg = np.atleast_2d(np.asarray(r_seg, float))
    reps, n = r_seg.shape
    center = c * scaled_mad(r_pre[:, -AWAKE_PRE:])
    levels, U, L, ep = awake_band_snapshots(r_seg, center, epochs, delta, n_plan)
    U = U * (1.0 + null_band)
    L = L * (1.0 - null_band)
    hist = np.concatenate([r_pre, r_seg], axis=1)
    k = running_level_index(hist, r_pre.shape[1], n, center, c=c, fixed=not running, window=window, every=every)
    E = ep.size
    which = np.searchsorted(ep, np.arange(n), side="right") - 1  # snapshot index of a chain started at s
    ug = np.asarray(up_grid, float)
    df = np.asarray(down_fracs, float)
    Su = np.zeros((reps, E, ug.size))
    Sd = np.zeros((reps, E, df.size))
    rows = np.arange(reps)
    levels_list = sorted(levels_list)
    out = {Lv: [np.full(reps, -1), np.full(reps, -1)] for Lv in levels_list}
    top = levels_list[-1]
    for t in range(n):
        e_now = which[t]
        kt = k[:, t]
        tau2 = levels[rows, kt] ** 2  # (reps,)
        x = np.minimum(r_seg[:, t] ** 2, tau2)
        cols = np.arange(e_now + 1)[None, :]
        u = U[rows[:, None], cols, kt[:, None]]  # (reps, e_now + 1): each snapshot's bound at tau_t
        l_ = L[rows[:, None], cols, kt[:, None]]
        fu = 1.0 + ug[None, None, :] * ((x[:, None] / u)[:, :, None] - 1.0)
        on = (l_ > 0) & (tau2[:, None] > l_)
        ls = np.where(on, l_, 1.0)
        B = np.where(on, tau2[:, None] / ls, 2.0)
        lam = np.where(on[:, :, None], df[None, None, :] / (B[:, :, None] - 1.0), 0.0)
        fd = 1.0 + lam * (1.0 - (x[:, None] / ls)[:, :, None])
        add = np.zeros((1, e_now + 1, 1))
        add[0, e_now, 0] = 1.0
        with np.errstate(over="ignore"):
            Su[:, : e_now + 1] = fu * (Su[:, : e_now + 1] + add)
            Sd[:, : e_now + 1] = fd * (Sd[:, : e_now + 1] + add)
        Mu = Su.mean(axis=2).sum(axis=1)
        Md = Sd.mean(axis=2).sum(axis=1)
        for Lv in levels_list:
            au, ad = out[Lv]
            au[(au < 0) & (Mu >= Lv)] = t
            ad[(ad < 0) & (Md >= Lv)] = t
        if t % 64 == 0 and (out[top][0] >= 0).all() and (out[top][1] >= 0).all():
            break
    return out, center


class AwakeRunningTauMonitor:
    """Route A with no blind period after an alarm (streaming).

    Start-up: the first ``n_pre`` finite returns of the stream centre the grid (the only
    blind stretch). A segment then monitors from its first return: a K-level learning band
    (``_LevelsCS``, total failure probability ``delta``, bet schedule ``n_plan``) is updated
    with each of the first ``n_learn`` returns and then frozen; snapshots are taken at the
    sizes ``awake_epochs(n_learn)``; the chain started at step s bets against the latest
    snapshot taken before s, with the running ``tau_t`` of ``RunningTauMonitor``; one
    Shiryaev-Roberts statistic per snapshot, summed; each side alarms at ``2 A``. After an
    alarm a new segment starts at the next return, its grid centred on ``c`` x the scaled MAD
    of the last ``n_pre`` finite returns (all before the new segment). Guarantee (per segment,
    i.i.d. returns within the segment): P(false alarm within m) <= delta + m / A.
    """

    def __init__(self, threshold: float = 36000.0, c: float = 3.0, delta: float = 0.01, null_band: float = 0.0,
                 running: bool = True, window: int = RUN_WINDOW, every: int = RUN_EVERY, n_pre: int = AWAKE_PRE,
                 n_learn: int = AWAKE_LEARN, n_plan: float | None = AWAKE_N_PLAN, up_grid=DEFAULT_UP_GRID,
                 down_fractions=DEFAULT_DOWN_FRACTIONS):
        self.A = float(threshold)
        self.c = c
        self.delta = delta
        self.band = null_band
        self.running = running
        self.window = window
        self.every = every
        self.n_pre = n_pre
        self.n_learn = n_learn
        self.n_plan = n_plan
        self.epochs = set(int(e) for e in awake_epochs(n_learn))
        self.up_grid = np.asarray(up_grid, float)
        self.down_fractions = np.asarray(down_fractions, float)
        self.K = GRID_KMAX - GRID_KMIN + 1
        self.t = -1
        self.events: list[Event] = []
        self.segments: list[dict] = []
        self.recent: deque = deque(maxlen=max(window, n_pre))
        self.calibrated = False

    def _center(self):
        r = np.asarray(self.recent)[-self.n_pre :]
        s = scaled_mad(r)
        if not s > 0:
            nz = np.abs(r[r != 0])
            s = float(np.median(nz) * MAD_NORMAL_SCALE) if nz.size else 0.0
        return self.c * s

    def _start(self):
        center = self._center()
        if not center > 0:
            self.calibrated = False
            return
        self.center = center
        self.levels = grid_levels(np.array([center]))[0]
        self.cs = _LevelsCS(self.levels, self.delta / (2 * self.K), self.n_plan)
        self.snaps_U: list[np.ndarray] = []
        self.snaps_L: list[np.ndarray] = []
        self.S_up = np.zeros((0, self.up_grid.size))
        self.S_dn = np.zeros((0, self.down_fractions.size))
        self.n_seg = 0
        self.seg_returns: list[float] = []
        self.start = self.t + 1
        self.calibrated = True
        self._set_level()
        self.segments.append(dict(start=self.start, center=center))

    def _set_level(self):
        if not self.running:
            self.k = -GRID_KMIN
            return
        s = self.c * scaled_mad(np.asarray(self.recent)[-self.window :])
        k = round(math.log(s / self.center) / math.log(GRID_RATIO)) if s > 0 else GRID_KMIN
        self.k = int(min(max(k, GRID_KMIN), GRID_KMAX)) - GRID_KMIN

    def update(self, r: float) -> Event | None:
        self.t += 1
        if not math.isfinite(r):
            return None
        if not self.calibrated:
            self.recent.append(r)
            if len(self.recent) >= self.n_pre:
                self._start()
            return None
        s = self.n_seg
        if s > 0 and s % self.every == 0:
            self._set_level()
        if s in self.epochs:
            u, l_ = self.cs.bounds()
            self.snaps_U.append(u * (1.0 + self.band))
            self.snaps_L.append(l_ * (1.0 - self.band))
            self.S_up = np.vstack([self.S_up, np.zeros(self.up_grid.size)])
            self.S_dn = np.vstack([self.S_dn, np.zeros(self.down_fractions.size)])
        k = self.k
        tau2 = self.levels[k] ** 2
        x = min(r * r, tau2)
        U = np.array([v[k] for v in self.snaps_U])
        L = np.array([v[k] for v in self.snaps_L])
        fu = 1.0 + self.up_grid[None, :] * (x / U[:, None] - 1.0)
        on = (L > 0) & (tau2 > L)
        Ls = np.where(on, L, 1.0)
        lam = np.where(on[:, None], self.down_fractions[None, :] / (np.where(on, tau2 / Ls, 2.0)[:, None] - 1.0), 0.0)
        fd = 1.0 + lam * (1.0 - (x / Ls)[:, None])
        self.S_up[-1] += 1.0
        self.S_dn[-1] += 1.0
        self.S_up *= fu
        self.S_dn *= fd
        if s < self.n_learn:
            self.cs.update(r)
            self.seg_returns.append(r)
        self.recent.append(r)
        self.n_seg += 1
        mu = float(self.S_up.mean(axis=1).sum())
        md = float(self.S_dn.mean(axis=1).sum())
        ev = None
        if mu >= 2 * self.A or md >= 2 * self.A:
            seg = np.asarray(self.seg_returns)
            scale = float(scaled_mad(seg)) if seg.size >= 30 else self.center / self.c
            ev = Event(self.t, "up" if mu >= 2 * self.A else "down", self.start, scale, max(mu, md))
            self.events.append(ev)
            self._start()
        return ev

    def run(self, returns) -> list[Event]:
        for r in np.asarray(returns, float):
            self.update(r)
        return self.events
