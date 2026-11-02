"""Predictable volatility forecasts, used to monitor returns *relative to a forecast*.

Why this module exists (E9, E10): the e-detector's guarantee needs the
conditional mean of the monitored statistic to stay at its baseline. Under
persistent volatility clustering (GARCH-type) that is false for raw returns:
the conditional variance wanders, and the monitor correctly reports those moves
as changes. If the question is instead "has volatility moved away from what my
model predicted?", the statistic should be built from standardised returns
``z_t = r_t / sigma_hat_t``, where ``sigma_hat_t`` uses data up to t - 1 only.
Under the null "the forecast is right up to a constant", ``z_t`` has a constant
conditional variance and the guarantee applies again (exactly with the true
conditional volatility, approximately with an estimate).

All functions treat NaN returns as missing: the forecast is carried forward
unchanged and the standardised value is NaN.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize


def ewma_variance(r: np.ndarray, halflife: float, init: float | None = None, floor: float = 0.0) -> np.ndarray:
    """Predictable EWMA of squared returns: ``v_t = lam v_{t-1} + (1 - lam) r_{t-1}^2``.

    Entry ``t`` uses ``r[:t]`` only. ``lam = 0.5 ** (1 / halflife)``. The start
    value is ``init`` or the mean of the first 100 finite squared returns
    (a small look-ahead at the very start; pass ``init`` to avoid it).
    """
    r = np.asarray(r, float)
    lam = 0.5 ** (1.0 / halflife)
    if init is None:
        head = r[: min(r.size, 100)]
        head = head[np.isfinite(head)]
        init = float(np.mean(head**2)) if head.size else 1.0
    if not np.isfinite(init):
        raise ValueError("init must be finite (a NaN start value would make every forecast NaN)")
    out = np.empty(r.size)
    v = max(float(init), floor)
    for t in range(r.size):
        out[t] = v
        x = r[t]
        if x == x:  # not NaN
            v = max(lam * v + (1 - lam) * x * x, floor)
    return out


def garch11_filter(r: np.ndarray, omega: float, alpha: float, beta: float, init: float | None = None) -> np.ndarray:
    """Predictable GARCH(1,1) variance ``s2_t = omega + alpha r_{t-1}^2 + beta s2_{t-1}``.

    NaN returns leave the variance at its previous value.
    """
    r = np.asarray(r, float)
    s2 = omega / max(1 - alpha - beta, 1e-12) if init is None else float(init)
    out = np.empty(r.size)
    for t in range(r.size):
        out[t] = s2
        x = r[t]
        if x == x:
            s2 = omega + alpha * x * x + beta * s2
    return out


def _garch_nll(theta, r, var_target):
    a, b = theta
    if a < 0 or b < 0 or a + b >= 0.99999:
        return 1e100
    omega = var_target * (1 - a - b)
    s2 = garch11_filter(r, omega, a, b, init=var_target)
    ok = np.isfinite(r)
    return 0.5 * float(np.sum(np.log(s2[ok]) + r[ok] ** 2 / s2[ok]))


def garch11_fit(r: np.ndarray, start=(0.05, 0.90)) -> tuple[float, float, float]:
    """Gaussian quasi maximum likelihood fit of GARCH(1,1) with variance targeting.

    ``omega`` is set so that the unconditional variance equals the sample
    variance (Engle and Mezrich 1996), and ``(alpha, beta)`` are found by
    Nelder-Mead. Returns ``(omega, alpha, beta)``. Meant for a training window,
    after which the parameters are frozen for monitoring.
    """
    r = np.asarray(r, float)
    var_target = float(np.nanvar(r))
    res = optimize.minimize(_garch_nll, np.asarray(start, float), args=(r, var_target), method="Nelder-Mead",
                            options=dict(xatol=1e-5, fatol=1e-3, maxiter=400))
    a, b = (float(v) for v in res.x)
    return var_target * (1 - a - b), a, b


def standardise(r: np.ndarray, variance_forecast: np.ndarray) -> np.ndarray:
    """``z_t = r_t / sqrt(v_t)`` with NaN where the return or the forecast is missing."""
    r = np.asarray(r, float)
    v = np.asarray(variance_forecast, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = r / np.sqrt(v)
    z[~np.isfinite(z)] = np.nan
    return z
