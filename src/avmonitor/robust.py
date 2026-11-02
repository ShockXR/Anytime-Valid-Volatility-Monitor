"""Realized quarticity estimators, MAD truncation and kernel-weighted local measures.

Conventions used throughout this module
---------------------------------------
* ``r`` holds ``n`` log returns observed on an interval ``[0, T]`` with equal
  spacing ``delta = T / n``. Functions accept a 1-D array (one path) or a 2-D
  array of shape ``(reps, n)`` (one path per row) and work along the last axis.
* Integrated quarticity is ``IQ = int_0^T sigma_t^4 dt``. With ``T = 1`` this is
  the scale of the realized quarticity ``RQ = (n / 3) sum r_i^4`` in eq. (3) of
  Bapat (2026). Setting ``T = n`` gives the "per step" scale, where the target is
  the average of ``sigma_i^4`` for per-step volatilities ``sigma_i``.

Two families of kernel estimators are provided:

* ``*_printed`` functions follow equations (6), (7), (12) and (16) of Bapat
  (2026) literally, including the unnormalised kernel weights. They exist so that
  the published numbers can be compared with the formulas as printed.
* ``local_*`` and ``kernel_iq`` use weights that sum to one and an explicit
  ``1 / delta`` scaling, so they estimate spot quantities and IQ on a fixed,
  documented scale. See ``docs/estimator_issues.md``.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

# Consistency constant that turns the raw MAD into an estimate of the standard
# deviation for normal data: 1 / Phi^{-1}(3/4).
MAD_NORMAL_SCALE = 1.0 / stats.norm.ppf(0.75)  # 1.482602...

# ----------------------------------------------------------------------------
# Kernels (Table 2 of Bapat 2026)
# ----------------------------------------------------------------------------


def _kernel_values(u: np.ndarray, kernel: str) -> np.ndarray:
    a = np.abs(u)
    inside = a <= 1.0
    if kernel == "gaussian":
        return np.exp(-0.5 * u**2) / np.sqrt(2 * np.pi)
    if kernel == "epanechnikov":
        return 0.75 * (1 - u**2) * inside
    if kernel == "uniform":
        return 0.5 * inside
    if kernel == "triangular":
        return (1 - a) * inside
    if kernel == "quartic":
        return 15 / 16 * (1 - u**2) ** 2 * inside
    if kernel == "triweight":
        return 35 / 32 * (1 - u**2) ** 3 * inside
    if kernel == "tricube":
        return 70 / 81 * (1 - a**3) ** 3 * inside
    if kernel == "cosine":
        return np.pi / 4 * np.cos(np.pi / 2 * u) * inside
    if kernel == "logistic":
        return 1.0 / (np.exp(u) + 2 + np.exp(-u))
    if kernel == "sigmoid":
        return 2 / np.pi / (np.exp(u) + np.exp(-u))
    if kernel == "silverman":
        return 0.5 * np.exp(-a / np.sqrt(2)) * np.sin(a / np.sqrt(2) + np.pi / 4)
    raise ValueError(f"unknown kernel {kernel!r}")


KERNELS = (
    "gaussian",
    "epanechnikov",
    "uniform",
    "triangular",
    "quartic",
    "triweight",
    "tricube",
    "cosine",
    "logistic",
    "sigmoid",
    "silverman",
)


def kernel_weights(h: int, kernel: str = "gaussian", normalise: bool = True) -> np.ndarray:
    """One-sided weights for lags ``-h, ..., -1, 0`` (oldest first, newest last).

    The newest return (lag 0) gets ``K(0)``, the largest weight for every kernel
    in Table 2. With ``normalise=False`` the raw values ``K((i - t) / h)`` of
    eq. (6) are returned. Their sum grows roughly linearly in ``h``.
    """
    if h < 1:
        raise ValueError("bandwidth h must be at least 1")
    u = np.arange(-h, 1) / h
    w = _kernel_values(u, kernel).astype(float)
    if normalise:
        w = w / w.sum()
    return w


def kernel_weight_sum(h: int, kernel: str = "gaussian") -> float:
    """Sum of the unnormalised weights ``sum_{k=0}^{h} K(k / h)``."""
    return float(kernel_weights(h, kernel, normalise=False).sum())


def _windowed(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted sums over trailing windows of length ``len(w)``.

    Output position ``k`` uses ``x[..., k : k + len(w)]`` and therefore refers to
    time ``t = k + len(w) - 1``. The output has ``n - len(w) + 1`` entries.
    """
    x = np.asarray(x, float)
    L = len(w)
    m = x.shape[-1] - L + 1
    if m < 1:
        raise ValueError("series shorter than the kernel window")
    out = np.zeros(x.shape[:-1] + (m,))
    for k in range(L):
        out += w[k] * x[..., k : k + m]
    return out


# ----------------------------------------------------------------------------
# MAD and truncation
# ----------------------------------------------------------------------------


def mad(r: np.ndarray, scaled: bool = False, center: str = "median") -> np.ndarray:
    """Median absolute deviation along the last axis.

    ``scaled=True`` multiplies by 1.4826 so that the result estimates the
    standard deviation for normal data (as R's ``mad()`` does by default).
    ``scaled=False`` returns the raw MAD (as ``scipy.stats.median_abs_deviation``
    does by default). Bapat (2026) does not say which one is used.
    ``center="zero"`` uses ``median(|r|)``, which is natural for returns whose
    mean is negligible and is cheaper to update online.
    """
    r = np.asarray(r, dtype=float)
    if center == "median":
        m = np.median(r, axis=-1, keepdims=True)
        out = np.median(np.abs(r - m), axis=-1)
    elif center == "zero":
        out = np.median(np.abs(r), axis=-1)
    else:
        raise ValueError("center must be 'median' or 'zero'")
    return out * MAD_NORMAL_SCALE if scaled else out


def truncate(r: np.ndarray, tau) -> np.ndarray:
    """Eq. (11): clip returns to ``[-tau, tau]``. ``tau`` may be per row or per entry."""
    tau = np.asarray(tau, dtype=float)
    if tau.ndim == 1 and np.ndim(r) == 2:
        tau = tau[:, None]
    return np.clip(r, -tau, tau)


class DegenerateThresholdWarning(UserWarning):
    """The MAD threshold is zero although some returns are not, so truncation erases them."""


def mad_threshold(r: np.ndarray, c: float = 3.0, scaled: bool = False, floor: float = 0.0) -> np.ndarray:
    """``tau = c * MAD(r)`` with an optional floor (useful for tick-size assets).

    When more than half of the returns are exactly zero (coarse price ticks,
    frozen prices), the MAD is zero, every return is clipped to zero and any
    truncated estimator silently returns 0. That case raises a
    ``DegenerateThresholdWarning``; pass ``floor`` (for example
    ``tick_floor(prices)``) to avoid it.
    """
    tau = np.maximum(c * mad(r, scaled=scaled), floor)
    r = np.asarray(r, float)
    nonzero = np.any((r != 0) & np.isfinite(r), axis=-1)
    if np.any((np.asarray(tau) <= 0) & nonzero):
        import warnings

        warnings.warn(
            "MAD threshold is zero (more than half of the returns are exactly zero), so every return is "
            "truncated to zero; pass a positive floor such as tick_floor(prices)",
            DegenerateThresholdWarning,
            stacklevel=2,
        )
    return tau


def tick_floor(prices: np.ndarray, c: float = 3.0, ticks: float = 1.0) -> float:
    """A truncation floor of ``c * ticks`` one-tick log returns.

    The tick is estimated as the smallest positive absolute price change and
    the one-tick log return at the median price. For a coin quoted to 4
    decimals at 0.20 USD this gives about ``c * 5e-4``. Use it as ``floor`` in
    ``mad_threshold`` (and ``scale_floor`` in the monitor, divided by ``c``) for
    assets whose prices move in discrete ticks.
    """
    p = np.asarray(prices, float)
    p = p[np.isfinite(p)]
    d = np.abs(np.diff(p))
    d = d[d > 0]
    if d.size == 0 or p.size == 0:
        return 0.0
    tick = float(np.min(d))
    return float(c * ticks * np.log1p(tick / np.median(p)))


def truncated_moment_factor(a: float, power: int = 4) -> float:
    """``E[clip(Z, -a, a)^p] / E[Z^p]`` for standard normal ``Z`` and p in {2, 4}.

    This is the fraction of the second or fourth moment that survives
    truncation at ``a`` standard deviations. With ``tau = c * MAD``, use
    ``a = c * Phi^{-1}(3/4)`` for the raw MAD and ``a = c`` for the scaled MAD.
    """
    phi = stats.norm.pdf(a)
    tail = 2 * stats.norm.sf(a)  # P(|Z| > a)
    central2 = (1 - tail) - 2 * a * phi  # E[Z^2 1{|Z|<a}]
    if power == 2:
        return central2 + a**2 * tail
    if power == 4:
        central4 = 3 * central2 - 2 * a**3 * phi  # E[Z^4 1{|Z|<a}]
        return (central4 + a**4 * tail) / 3.0
    raise ValueError("power must be 2 or 4")


def effective_cutoff(c: float, scaled: bool) -> float:
    """Truncation point in standard deviations implied by ``tau = c * MAD`` for normal data."""
    return c if scaled else c * stats.norm.ppf(0.75)


# ----------------------------------------------------------------------------
# Classical integrated quarticity estimators (all on the IQ scale)
# ----------------------------------------------------------------------------


def _mu(p: float) -> float:
    """``E|Z|^p`` for standard normal ``Z``."""
    from scipy.special import gamma

    return 2 ** (p / 2) * gamma((p + 1) / 2) / np.sqrt(np.pi)


def _delta(n: int, T: float) -> float:
    return T / n


def realized_quarticity(r, T: float = 1.0):
    """Eq. (3): ``RQ = (1 / (3 delta)) sum r_i^4``; equals ``(n/3) sum r^4`` when T = 1."""
    r = np.asarray(r, float)
    n = r.shape[-1]
    return np.sum(r**4, axis=-1) / (3 * _delta(n, T))


def tripower_quarticity(r, T: float = 1.0):
    """Standard tripower quarticity with p = 4/3 (Barndorff-Nielsen and Shephard 2004)."""
    r = np.abs(np.asarray(r, float))
    n = r.shape[-1]
    p = 4.0 / 3.0
    prod = r[..., 2:] ** p * r[..., 1:-1] ** p * r[..., :-2] ** p
    return (n / (n - 2)) * np.sum(prod, axis=-1) / (_mu(p) ** 3 * _delta(n, T))


def quadpower_quarticity(r, T: float = 1.0):
    """Quad-power quarticity ``mu_1^{-4} sum |r_i||r_{i-1}||r_{i-2}||r_{i-3}|`` scaled to IQ."""
    r = np.abs(np.asarray(r, float))
    n = r.shape[-1]
    prod = r[..., 3:] * r[..., 2:-1] * r[..., 1:-2] * r[..., :-3]
    return (n / (n - 3)) * np.sum(prod, axis=-1) / (_mu(1.0) ** 4 * _delta(n, T))


def bipower_quarticity(r, T: float = 1.0):
    """Eq. (4) form ``sum |r_i| |r_{i-1}|^3`` with the constant and scaling that make it unbiased for IQ.

    For normal increments ``E|Z| E|Z|^3 = mu_1 mu_3``, so the sum is divided by
    ``mu_1 mu_3 delta``. Eq. (4) as printed uses ``mu_1^{-4}`` and no ``1/delta``.
    """
    r = np.abs(np.asarray(r, float))
    n = r.shape[-1]
    prod = r[..., 1:] * r[..., :-1] ** 3
    return (n / (n - 1)) * np.sum(prod, axis=-1) / (_mu(1.0) * _mu(3.0) * _delta(n, T))


def min_rq(r, T: float = 1.0):
    """MinRQ of Andersen, Dobrev and Schaumburg (2012)."""
    r = np.abs(np.asarray(r, float))
    n = r.shape[-1]
    m = np.minimum(r[..., 1:], r[..., :-1])
    return np.pi / (3 * np.pi - 8) * (n / (n - 1)) * np.sum(m**4, axis=-1) / _delta(n, T)


def med_rq(r, T: float = 1.0):
    """MedRQ of Andersen, Dobrev and Schaumburg (2012)."""
    r = np.abs(np.asarray(r, float))
    n = r.shape[-1]
    stack = np.stack([r[..., 2:], r[..., 1:-1], r[..., :-2]], axis=0)
    m = np.median(stack, axis=0)
    const = 3 * np.pi / (9 * np.pi + 72 - 52 * np.sqrt(3))
    return const * (n / (n - 2)) * np.sum(m**4, axis=-1) / _delta(n, T)


def bipower_variation(r):
    """Bipower variation ``mu_1^{-2} sum |r_i||r_{i-1}|``, a jump-robust estimate of ``int sigma^2 dt``."""
    r = np.abs(np.asarray(r, float))
    n = r.shape[-1]
    return (n / (n - 1)) * np.sum(r[..., 1:] * r[..., :-1], axis=-1) / _mu(1.0) ** 2


def threshold_rq(r, T: float = 1.0, c: float = 4.0, varpi: float = 0.49):
    """Threshold realized quarticity (Mancini 2009 style).

    Returns larger than ``u_n = c * sqrt(BV / T) * delta^varpi`` are dropped,
    where BV is bipower variation. Because ``varpi < 1/2`` the threshold shrinks
    more slowly than the diffusive returns, so the truncation bias vanishes as
    n grows. This is the condition the proof of Theorem 1 needs.
    """
    r = np.asarray(r, float)
    n = r.shape[-1]
    d = _delta(n, T)
    bv = bipower_variation(r)
    u = c * np.sqrt(np.maximum(bv, 0) / T) * d**varpi
    if np.ndim(u) > 0:
        u = u[..., None]
    keep = np.abs(r) <= u
    return np.sum((r**4) * keep, axis=-1) / (3 * d)


# ----------------------------------------------------------------------------
# Kernel estimators exactly as printed in Bapat (2026)
# ----------------------------------------------------------------------------


def local_quarticity_printed(r, h: int = 30, kernel: str = "gaussian", reverse_weights: bool = False):
    """Eq. (6)/(12) as printed: ``Q_t = (1/3) sum_{i=t-h}^{t} K((i-t)/h) r_i^4``.

    Weights are not normalised and there is no ``1/delta`` factor. With
    ``reverse_weights=True`` the indexing ``w_{s-t+h}`` of eq. (15) is used,
    which gives the largest weight to the oldest return in the window.
    """
    w = kernel_weights(h, kernel, normalise=False)
    if reverse_weights:
        w = w[::-1]
    return _windowed(np.asarray(r, float) ** 4, w) / 3.0


def kernel_quarticity_printed(r, h: int = 30, kernel: str = "gaussian", reverse_weights: bool = False):
    """Eq. (7)/(16) as printed: the average of eq. (6) over the ``T - h`` admissible ``t``."""
    return local_quarticity_printed(r, h, kernel, reverse_weights).mean(axis=-1)


def robust_kernel_quarticity_printed(
    r, h: int = 30, c: float = 3.0, kernel: str = "gaussian", mad_scaled: bool = False, reverse_weights: bool = False
):
    """Eqs. (11), (12) and (16) as printed: truncate at ``c * MAD`` of the whole sample, then eq. (16)."""
    r = np.asarray(r, float)
    tau = mad_threshold(r, c, scaled=mad_scaled)
    return kernel_quarticity_printed(truncate(r, tau), h, kernel, reverse_weights)


# ----------------------------------------------------------------------------
# Normalised kernel estimators (spot and integrated)
# ----------------------------------------------------------------------------


def local_variance(r, h: int = 30, kernel: str = "gaussian", T: float = 1.0):
    """Spot variance ``sigma_t^2`` from normalised one-sided kernel weights: ``sum w_i r_i^2 / delta``.

    Entry ``k`` of the output is the estimate at time index ``t = k + h``.
    """
    r = np.asarray(r, float)
    d = _delta(r.shape[-1], T)
    return _windowed(r**2, kernel_weights(h, kernel)) / d


def local_quarticity(r, h: int = 30, kernel: str = "gaussian", T: float = 1.0):
    """Spot quarticity ``sigma_t^4`` from normalised weights: ``sum w_i r_i^4 / (3 delta^2)``."""
    r = np.asarray(r, float)
    d = _delta(r.shape[-1], T)
    return _windowed(r**4, kernel_weights(h, kernel)) / (3 * d**2)


def kernel_iq(
    r,
    h: int = 30,
    kernel: str = "gaussian",
    T: float = 1.0,
    c: float | None = None,
    mad_scaled: bool = True,
    bias_correct: bool = False,
):
    """Integrated quarticity from normalised local estimates: ``T * mean_t(sigma_t^4 hat)``.

    ``c=None`` gives the plain kernel estimator. With ``c`` set, returns are
    first truncated at ``c * MAD`` (eq. 11). ``bias_correct=True`` divides by the
    surviving fourth-moment fraction for normal returns (see
    ``truncated_moment_factor``), which removes the first-order truncation bias
    when the diffusive part is locally Gaussian.
    """
    r = np.asarray(r, float)
    if c is not None:
        r = truncate(r, mad_threshold(r, c, scaled=mad_scaled))
    est = T * local_quarticity(r, h, kernel, T).mean(axis=-1)
    if c is not None and bias_correct:
        est = est / truncated_moment_factor(effective_cutoff(c, mad_scaled), 4)
    return est


def truncated_rq_mad(r, c: float = 3.0, T: float = 1.0, mad_scaled: bool = True, bias_correct: bool = False):
    """Realized quarticity of MAD-truncated returns, with no kernel at all.

    ``docs/estimator_issues.md`` shows that the aggregated kernel estimator of
    eq. (16) equals this statistic up to a constant factor and edge effects.
    """
    r = np.asarray(r, float)
    rt = truncate(r, mad_threshold(r, c, scaled=mad_scaled))
    est = realized_quarticity(rt, T)
    if bias_correct:
        est = est / truncated_moment_factor(effective_cutoff(c, mad_scaled), 4)
    return est


# ----------------------------------------------------------------------------
# Online (predictable) scale estimates for streaming use
# ----------------------------------------------------------------------------


def rolling_mad_zero(r: np.ndarray, window: int, min_periods: int | None = None, scaled: bool = True) -> np.ndarray:
    """Predictable rolling ``median(|r|)``: entry ``t`` uses ``r[t-window : t]`` only.

    Entries without enough history are NaN. With ``scaled=True`` the output
    estimates the return standard deviation for normal data.
    """
    import pandas as pd

    s = pd.Series(np.abs(np.asarray(r, float)))
    mp = window if min_periods is None else min_periods
    out = s.rolling(window, min_periods=mp).median().shift(1).to_numpy()
    return out * MAD_NORMAL_SCALE if scaled else out


def ew_quantile_scale(
    r: np.ndarray, halflife: float = 390.0, init: float | None = None, scaled: bool = True, floor: float = 0.0,
    strict: bool = False,
) -> np.ndarray:
    """Predictable exponentially weighted ``median(|r|)`` by stochastic approximation.

    Uses the multiplicative Robbins and Monro update
    ``m_t = m_{t-1} * exp(eta * (1{|r_t| > m_{t-1}} - 1/2))``, which tracks the
    median of ``|r|`` and stays positive. ``eta`` is set from ``halflife``.
    Entry ``t`` uses data before ``t`` only, with one exception: when ``init``
    is None the start value is the median of the first 50 values, so entries
    0 to 49 use those values (a start-up look-ahead, flagged by the independent
    review). Pass ``init`` from earlier data, or ``strict=True`` to return NaN
    for those entries, when strict predictability matters.

    Edge cases (iteration 3 of the log): NaN returns leave the estimate
    unchanged instead of counting as zeros; the estimate never falls below
    ``max(floor, 1e-300)``, so it cannot underflow to exactly zero and get
    stuck there when more than half of the returns are zero; and if the median
    of the first 50 values is zero, the start value is the median of the
    non-zero values instead (otherwise the tracker needs hundreds of thousands
    of steps to climb from 1e-300).
    """
    a = np.abs(np.asarray(r, float))
    n = a.size
    eta = np.log(2) / halflife * 2.0
    out = np.empty(n)
    lo = max(float(floor), 1e-300)
    if init is None:
        head = a[: min(n, 50)]
        head = head[np.isfinite(head)]
        m = float(np.median(head)) if head.size else 0.0
        if m <= 0:
            nz = head[head > 0]
            m = float(np.median(nz)) if nz.size else lo
    else:
        m = float(init)
    m = max(m, lo)
    up, down = np.exp(eta * 0.5), np.exp(-eta * 0.5)
    for t in range(n):
        out[t] = m
        at = a[t]
        if at != at:  # NaN: no information
            continue
        m = m * up if at > m else max(m * down, lo)
    if strict and init is None:
        out[: min(n, 50)] = np.nan
    return out * MAD_NORMAL_SCALE if scaled else out
