"""Confidence sequences for a bounded mean, and a calibration of the detector baseline that keeps the guarantee.

Why this module exists
----------------------
The e-detectors in ``detectors.py`` are valid when they bet against the TRUE
pre-change mean of ``X = min(r^2, tau^2)``. The monitor estimates that mean
from a calibration window, and a plug-in estimate that comes out too low makes
the increase detector anti-conservative (``docs/math.md``, 5.6 and 6). This
module replaces the plug-in by one-sided confidence bounds that hold
uniformly over time:

* ``U_t``: with probability at least ``1 - delta_upper``, ``U_t`` is above the
  (weighted average) conditional mean of the calibration data at EVERY t;
* ``L_t``: the same from below with ``delta_lower``.

The increase detector then bets against ``E[X] <= U_T`` and the decrease
detector against ``E[X] >= L_T``, where ``T`` is the end of calibration. If
there is no change, the chance of any false alarm within m monitored steps is
at most ``delta_upper + delta_lower + m / A`` (``docs/math.md``, Theorem 6.4).
Because the bounds hold at every t simultaneously, ``T`` may be chosen from
the data, for example "stop calibrating once the bounds are tight enough".

The confidence sequence
-----------------------
Predictable plug-in empirical Bernstein (PrPl-EB) sequence of Waudby-Smith and
Ramdas (2024, JRSS B 86(1), 1 to 27, doi:10.1093/jrsssb/qkad009, Theorem 2),
with one change: the bet cap uses the one-sided range of the centred variable
(``a_t`` below) instead of the full range, which allows much larger bets for
the lower bound of a right-skewed variable such as a truncated squared return.
With ``a_t = 1`` and ``delta_upper = delta_lower = alpha / 2`` the two bounds are
exactly their two-sided interval at level ``alpha``. Validity is proved in full
in ``docs/math.md``, 6.2, using only Ville's inequality.

For ``Z_i = X_i / B`` in ``[0, 1]`` (``B`` = the known bound, ``tau^2`` here):

* ``zhat_t = (1/2 + sum_{i<=t} Z_i) / (t + 1)``,
  ``s2_t = (1/4 + sum_{i<=t} (Z_i - zhat_i)^2) / (t + 1)`` (``zhat_0 = 1/2``, ``s2_0 = 1/4``);
* one-sided ranges ``a^up_t = 1 - zhat_{t-1}`` and ``a^lo_t = zhat_{t-1}``;
* bets ``lam_t = min( sqrt(2 log(1/delta) / (d_t s2_{t-1})), c / a_t )`` with
  ``d_t = n_plan`` (tuned for a planned sample size) or ``d_t = t log(1 + t)``
  (the time-uniform schedule of Waudby-Smith and Ramdas);
* ``psi_a(lam) = (-log(1 - lam a) - lam a) / a^2`` and ``e_t = (Z_t - zhat_{t-1})^2``;
* ``U_t = [sum lam^up Z + log(1/delta_up) + sum psi(lam^up) e] / sum lam^up``,
  ``L_t = [sum lam^lo Z - log(1/delta_lo) - sum psi(lam^lo) e] / sum lam^lo``,
  clipped to ``[0, 1]`` and multiplied by ``B``.

Any predictable choice of ``lam``, ``zhat`` and ``s2`` keeps validity; the
choices only move where the sequence is tightest.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .robust import MAD_NORMAL_SCALE

# ----------------------------------------------------------------------------
# Numerics
# ----------------------------------------------------------------------------


def psi_coef(x):
    """``h(x) = (-log(1 - x) - x) / x^2`` for ``x`` in ``[0, 1)``, with ``h(0) = 1/2``.

    Then ``psi_a(lam) = lam^2 h(lam a)``. A Taylor series is used for small x to
    avoid cancellation (``h(x) = sum_{k>=2} x^(k-2) / k``).
    """
    x = np.asarray(x, float)
    small = x < 1e-3
    xs = np.where(small, 0.5, x)
    direct = (-np.log1p(-xs) - xs) / (xs * xs)
    series = 0.5 + x / 3.0 + x * x / 4.0 + x**3 / 5.0
    return np.where(small, series, direct)


def _psi_coef_scalar(x: float) -> float:
    if x < 1e-3:
        return 0.5 + x / 3.0 + x * x / 4.0 + x**3 / 5.0
    return (-math.log1p(-x) - x) / (x * x)


def _check_params(delta_upper, delta_lower, n_plan, c):
    if not 0 < c < 1:
        raise ValueError("c must lie in (0, 1)")
    for d in (delta_upper, delta_lower):
        if not 0 < d < 1:
            raise ValueError("delta must lie in (0, 1)")
    if n_plan is not None and not (math.isfinite(n_plan) and n_plan > 0):
        raise ValueError("n_plan must be a positive number or None")


# ----------------------------------------------------------------------------
# Vectorised confidence sequence
# ----------------------------------------------------------------------------


@dataclass
class ConfSeq:
    """Bounds after each observation, in the original units of ``x``. Shape (reps, n)."""

    lower: np.ndarray
    upper: np.ndarray
    mean: np.ndarray  # running sample mean (point estimate, not used for validity)
    lam_upper: np.ndarray | None = None  # bets of the upper sequence (normalised units)
    lam_lower: np.ndarray | None = None  # bets of the lower sequence

    def rel_width(self) -> np.ndarray:
        """Relative width ``(upper - lower) / mean``; ``inf`` where the mean is not positive."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.mean > 0, (self.upper - self.lower) / self.mean, np.inf)


def eb_confseq(
    x,
    bound,
    delta_upper: float = 0.005,
    delta_lower: float = 0.005,
    n_plan: float | None = None,
    c: float = 0.75,
) -> ConfSeq:
    """One-sided PrPl-EB confidence sequences for the mean of ``x`` in ``[0, bound]``.

    Parameters
    ----------
    x : array (n,) or (reps, n)
        Observations, each in ``[0, bound]`` (no NaN). Rows are independent paths.
    bound : float or array (reps,)
        Known almost sure upper bound of ``x`` (per row if an array).
    delta_upper, delta_lower : float
        Error levels of the upper and the lower sequence: ``P(exists t: m >= U_t) <= delta_upper``
        and ``P(exists t: m <= L_t) <= delta_lower`` (see the module docstring for which m).
    n_plan : float or None
        Planned sample size for tuning the bets (``d_t = n_plan``). ``None`` uses
        the time-uniform schedule ``d_t = t log(1 + t)``. Validity holds either way.
    c : float in (0, 1)
        Cap on ``lam_t a_t``. Waudby-Smith and Ramdas suggest 1/2 or 3/4.

    Returns ``ConfSeq`` with arrays of shape (reps, n) (or (1, n) for 1-d input).
    """
    _check_params(delta_upper, delta_lower, n_plan, c)
    x = np.atleast_2d(np.asarray(x, float))
    reps, n = x.shape
    B = np.asarray(bound, float)
    B = B.reshape(-1, 1) if B.ndim else B
    if np.any(~np.isfinite(x)):
        raise ValueError("x must be finite (drop missing values first)")
    if np.any(B <= 0):
        raise ValueError("bound must be positive")
    z = x / B
    if np.any(z < 0) or np.any(z > 1 + 1e-12):
        raise ValueError("x must lie in [0, bound]")
    z = np.minimum(z, 1.0)
    t = np.arange(1, n + 1, dtype=float)
    zhat = (0.5 + np.cumsum(z, axis=1)) / (t + 1.0)
    zprev = np.concatenate([np.full((reps, 1), 0.5), zhat[:, :-1]], axis=1)
    s2 = (0.25 + np.cumsum((z - zhat) ** 2, axis=1)) / (t + 1.0)
    s2prev = np.concatenate([np.full((reps, 1), 0.25), s2[:, :-1]], axis=1)
    d = float(n_plan) if n_plan is not None else t * np.log1p(t)
    e = (z - zprev) ** 2

    def side(delta, a, sign):
        logd = math.log(1.0 / delta)
        lam = np.minimum(np.sqrt(2.0 * logd / (d * s2prev)), c / a)
        psi = lam * lam * psi_coef(lam * a)
        sl = np.cumsum(lam, axis=1)
        slz = np.cumsum(lam * z, axis=1)
        spe = np.cumsum(psi * e, axis=1)
        return (slz + sign * (logd + spe)) / sl, lam

    up, lam_u = side(delta_upper, 1.0 - zprev, +1.0)
    lo, lam_l = side(delta_lower, zprev, -1.0)
    up, lo = np.minimum(up, 1.0), np.maximum(lo, 0.0)
    mean = np.cumsum(z, axis=1) / t
    return ConfSeq(lower=lo * B, upper=up * B, mean=mean * B, lam_upper=lam_u, lam_lower=lam_l)


# ----------------------------------------------------------------------------
# Streaming confidence sequence (used by the monitor)
# ----------------------------------------------------------------------------


class StreamingEBConfSeq:
    """The same sequence as ``eb_confseq``, updated one observation at a time."""

    def __init__(self, bound: float, delta_upper: float = 0.005, delta_lower: float = 0.005,
                 n_plan: float | None = None, c: float = 0.75):
        if not (math.isfinite(bound) and bound > 0):
            raise ValueError("bound must be positive")
        _check_params(delta_upper, delta_lower, n_plan, c)
        self.bound = float(bound)
        self.log_u = math.log(1.0 / delta_upper)
        self.log_l = math.log(1.0 / delta_lower)
        self.n_plan = n_plan
        self.c = c
        self.t = 0
        self.sum_z = 0.0
        self.sum_dev2 = 0.0
        self.zhat = 0.5
        self.s2 = 0.25
        self.sl_u = self.slz_u = self.spe_u = 0.0
        self.sl_l = self.slz_l = self.spe_l = 0.0

    def update(self, x: float) -> None:
        """Add one observation ``x`` in ``[0, bound]``."""
        z = float(x) / self.bound
        if not (0.0 <= z <= 1.0 + 1e-12):
            raise ValueError("x must lie in [0, bound]")
        z = min(z, 1.0)
        self.t += 1
        t = self.t
        zprev, s2prev = self.zhat, self.s2
        d = float(self.n_plan) if self.n_plan is not None else t * math.log1p(t)
        e = (z - zprev) ** 2
        a = 1.0 - zprev
        lam = min(math.sqrt(2.0 * self.log_u / (d * s2prev)), self.c / a)
        self.sl_u += lam
        self.slz_u += lam * z
        self.spe_u += lam * lam * _psi_coef_scalar(lam * a) * e
        a = zprev
        lam = min(math.sqrt(2.0 * self.log_l / (d * s2prev)), self.c / a)
        self.sl_l += lam
        self.slz_l += lam * z
        self.spe_l += lam * lam * _psi_coef_scalar(lam * a) * e
        self.sum_z += z
        self.zhat = (0.5 + self.sum_z) / (t + 1.0)
        self.sum_dev2 += (z - self.zhat) ** 2
        self.s2 = (0.25 + self.sum_dev2) / (t + 1.0)

    @property
    def upper(self) -> float:
        """Current upper bound ``U_t`` in the units of ``x`` (``bound`` before any observation)."""
        if self.t == 0:
            return self.bound
        return min((self.slz_u + self.log_u + self.spe_u) / self.sl_u, 1.0) * self.bound

    @property
    def lower(self) -> float:
        """Current lower bound ``L_t`` in the units of ``x`` (0 before any observation)."""
        if self.t == 0:
            return 0.0
        return max((self.slz_l - self.log_l - self.spe_l) / self.sl_l, 0.0) * self.bound

    @property
    def mean(self) -> float:
        """Running sample mean (NaN before any observation)."""
        return self.sum_z / self.t * self.bound if self.t else float("nan")

    @property
    def rel_width(self) -> float:
        """Relative width ``(upper - lower) / mean``; ``inf`` before any observation or when the mean is 0."""
        m = self.mean
        return (self.upper - self.lower) / m if self.t and m > 0 else float("inf")


# ----------------------------------------------------------------------------
# Valid baseline from a calibration window (sample splitting)
# ----------------------------------------------------------------------------


@dataclass
class CSBaseline:
    """Baseline with a sample split: tau from the first ``n_tau`` returns, bounds on E[X] from the rest.

    All fields are arrays of shape (reps,) (one entry per calibration window).
    """

    scale: np.ndarray  # scaled MAD of the first n_tau returns (sets tau)
    tau: np.ndarray
    mean: np.ndarray  # plug-in mean of X over the second part (point estimate)
    upper: np.ndarray  # U_T, upper confidence bound for E[X]
    lower: np.ndarray  # L_T, lower confidence bound for E[X]
    n_tau: int
    n_cs: np.ndarray  # number of observations used by the confidence sequence


def split_tau(r_tau, c: float = 3.0, floor: float = 0.0):
    """``scale`` (scaled MAD about zero) and ``tau = max(c * scale, floor)`` per row."""
    r_tau = np.atleast_2d(np.asarray(r_tau, float))
    scale = np.median(np.abs(r_tau), axis=1) * MAD_NORMAL_SCALE
    tau = np.maximum(c * scale, floor)
    return scale, tau


def cs_baseline(
    r_cal,
    c: float = 3.0,
    delta_upper: float = 0.005,
    delta_lower: float = 0.005,
    tau_frac: float = 0.25,
    n_plan: float | None = None,
    cs_c: float = 0.75,
    floor: float = 0.0,
    stop: np.ndarray | None = None,
) -> CSBaseline:
    """Valid baseline from calibration windows ``r_cal`` (rows = paths, no NaN).

    ``tau`` comes from the first ``n_tau = round(tau_frac * n)`` returns, so it is
    fixed before the confidence sequence starts (``docs/math.md``, 6.4). The
    sequence runs on ``X = min(r^2, tau^2)`` for the remaining returns. If
    ``stop`` is given (per row, number of sequence observations, at least 1),
    the bounds are read at that data-dependent time instead of the end, which is
    allowed because the sequence is time-uniform. ``n_plan`` defaults to the
    number of sequence observations of a full window.
    """
    r_cal = np.atleast_2d(np.asarray(r_cal, float))
    reps, n = r_cal.shape
    n_tau = int(round(tau_frac * n))
    if not 1 <= n_tau < n:
        raise ValueError("tau_frac must leave at least one return on each side of the split")
    scale, tau = split_tau(r_cal[:, :n_tau], c, floor)
    if np.any(tau <= 0):
        raise ValueError("tau is zero on some rows (too many zero returns); pass floor > 0")
    tau2 = tau**2
    x = np.minimum(r_cal[:, n_tau:] ** 2, tau2[:, None])
    plan = n - n_tau if n_plan is None else n_plan
    cs = eb_confseq(x, tau2, delta_upper, delta_lower, n_plan=plan, c=cs_c)
    if stop is None:
        k = np.full(reps, n - n_tau)
    else:
        k = np.asarray(stop, int)
        if np.any(k < 1) or np.any(k > n - n_tau):
            raise ValueError("stop must lie in 1 .. n - n_tau")
    i = np.arange(reps)
    return CSBaseline(scale=scale, tau=tau, mean=cs.mean[i, k - 1], upper=cs.upper[i, k - 1],
                      lower=cs.lower[i, k - 1], n_tau=n_tau, n_cs=k)


def width_stopping_time(cs: ConfSeq, target: float, n_min: int = 1, n_max: int | None = None) -> np.ndarray:
    """First t (1-based count) with ``t >= n_min`` and relative width ``(U - L) / mean <= target``.

    Rows that never reach the target stop at ``n_max`` (default: the last column).
    """
    if n_min < 1:
        raise ValueError("n_min must be at least 1")
    w = cs.rel_width()
    reps, n = w.shape
    n_max = n if n_max is None else min(n_max, n)
    ok = w <= target
    ok[:, : n_min - 1] = False
    ok[:, n_max:] = False
    first = np.where(ok.any(axis=1), ok.argmax(axis=1) + 1, n_max)
    return first
