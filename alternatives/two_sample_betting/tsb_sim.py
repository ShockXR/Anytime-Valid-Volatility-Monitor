"""Simulators and the vectorised main-monitor benchmark for the two-sample betting prototype.

All laws have unit variance unless a scale factor says otherwise, so a "change in
shape at fixed variance" really keeps the variance fixed.

* ``draw(spec, shape, rng)``: i.i.d. returns from a spec dict:
  ``{"law": "gauss"}``, ``{"law": "t", "df": 4}`` (variance-matched Student t),
  ``{"law": "jumps", "p": 0.002, "J": 6.0}`` (Gaussian diffusion plus Gaussian
  jumps of sd ``J`` with probability ``p`` per step, diffusion variance
  ``1 - p J^2`` so the total variance is 1), ``{"law": "signskew", "thr": 1.5,
  "p_neg": 0.8}`` (Gaussian magnitudes, but moves larger than ``thr`` in absolute
  value are negative with probability ``p_neg``: the law of ``|u|`` is exactly
  half-normal, only the sign of large moves changes). Every spec accepts
  ``"sd"`` (a multiplier, default 1).
* ``change_paths``: calibration block, pre-change monitored block, post-change block.
* ``garch_paths``: GARCH(1, 1) with unit unconditional variance (volatility clustering).
* ``main_monitor_first_alarm``: the plug-in ``avmonitor.monitor.VolatilityMonitor``
  (calibrate on the first ``cal_len`` returns, two-sided Shiryaev-Roberts at
  ``2 A`` per side) for many paths at once; first alarm only. A test checks it
  against the class itself.
"""

from __future__ import annotations

import numpy as np

from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from avmonitor.robust import MAD_NORMAL_SCALE


def draw(spec: dict, shape, rng: np.random.Generator) -> np.ndarray:
    law = spec.get("law", "gauss")
    sd = float(spec.get("sd", 1.0))
    if law == "gauss":
        x = rng.standard_normal(shape)
    elif law == "t":
        df = float(spec["df"])
        x = rng.standard_t(df, shape) * np.sqrt((df - 2.0) / df)
    elif law == "jumps":
        p, J = float(spec["p"]), float(spec["J"])
        if p * J * J >= 1:
            raise ValueError("jump variance must be below 1")
        x = np.sqrt(1.0 - p * J * J) * rng.standard_normal(shape)
        x += (rng.random(shape) < p) * rng.normal(0.0, J, shape)
    elif law == "signskew":
        z = np.abs(rng.standard_normal(shape))
        big = z >= float(spec["thr"])
        neg = np.where(big, rng.random(shape) < float(spec["p_neg"]), rng.random(shape) < 0.5)
        x = np.where(neg, -z, z)
    else:
        raise ValueError(f"unknown law {law}")
    return sd * x


def change_paths(reps: int, n_cal: int, nu: int, n_post: int, pre: dict, post: dict | None,
                 rng: np.random.Generator) -> np.ndarray:
    """Rows: ``n_cal`` calibration steps and ``nu`` monitored steps from ``pre``, then ``n_post`` from ``post``."""
    a = draw(pre, (reps, n_cal + nu), rng)
    b = draw(pre if post is None else post, (reps, n_post), rng)
    return np.concatenate([a, b], axis=1)


def garch_paths(reps: int, n: int, persistence: float, alpha: float, rng: np.random.Generator,
                burn: int = 3000) -> np.ndarray:
    """GARCH(1, 1) with Gaussian innovations and unconditional variance 1."""
    beta = persistence - alpha
    omega = 1.0 - persistence
    h = np.ones(reps)
    out = np.empty((reps, n))
    z = rng.standard_normal((reps, n + burn))
    r = np.zeros(reps)
    for t in range(n + burn):
        h = omega + alpha * r * r + beta * h
        r = np.sqrt(h) * z[:, t]
        if t >= burn:
            out[:, t - burn] = r
    return out


def main_monitor_first_alarm(r: np.ndarray, threshold: float, cal_len: int = 1440, c: float = 3.0,
                             null_band: float = 0.0, up_grid=DEFAULT_UP_GRID,
                             down_fractions=DEFAULT_DOWN_FRACTIONS):
    """First alarm (monitored step index, -1 if none) and direction of the plug-in VolatilityMonitor.

    Same arithmetic as ``avmonitor.monitor.VolatilityMonitor`` with ``baseline="plugin"``
    on paths without missing values and with a usable calibration window.
    """
    r = np.atleast_2d(np.asarray(r, float))
    reps, n = r.shape
    cal = r[:, :cal_len]
    scale = np.median(np.abs(cal), axis=1) * MAD_NORMAL_SCALE
    tau2 = (c * scale) ** 2
    m0 = np.minimum(cal**2, tau2[:, None]).mean(axis=1)
    m_up = m0 * (1 + null_band)
    m_dn = m0 * (1 - null_band)
    up = np.asarray(up_grid, float)
    bb = tau2 / m_dn
    lam_dn = np.asarray(down_fractions, float)[None, :] / (bb[:, None] - 1.0)
    S_up = np.zeros((reps, up.size))
    S_dn = np.zeros((reps, lam_dn.shape[1]))
    alarm = np.full(reps, -1)
    direction = np.full(reps, "", dtype=object)
    active = np.ones(reps, bool)
    A2 = 2.0 * threshold
    for t in range(cal_len, n):
        x = np.minimum(r[:, t] ** 2, tau2)
        with np.errstate(over="ignore"):
            S_up = (1.0 + up[None, :] * (x / m_up - 1.0)[:, None]) * (S_up + 1.0)
            S_dn = (1.0 + lam_dn * (1.0 - x / m_dn)[:, None]) * (S_dn + 1.0)
        mu, md = S_up.mean(axis=1), S_dn.mean(axis=1)
        hit_up = active & (mu >= A2)
        hit_dn = active & ~hit_up & (md >= A2)
        alarm[hit_up | hit_dn] = t - cal_len
        direction[hit_up] = "up"
        direction[hit_dn] = "down"
        active &= ~(hit_up | hit_dn)
        if not active.any():
            break
    return alarm, direction
