"""Simulators for high-frequency returns with jumps, noise and known volatility changes.

All output is synthetic. Three generators are provided:

* ``paper_returns``: the discrete design of Bapat (2026), eqs. (10) and (13):
  ``r_i = sigma_i * eps_i + B_i * Z_i + alpha_i`` with Bernoulli(p) jumps of size
  ``N(mu_J, sigma_J^2)`` and i.i.d. ``N(0, sigma_alpha^2)`` noise added to returns.
* ``sv_jump_returns``: a continuous-time log-volatility model sampled on a fine
  grid, so that the true integrated quarticity is known, with compound Poisson
  jumps and either return-level noise (as in the paper) or price-level noise
  (the standard microstructure model, which makes returns MA(1)).
* ``change_returns``: piecewise volatility with a known change point, for the
  detection experiments. Optional intraday seasonality, slow stochastic
  volatility, jumps and price-level noise.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def paper_returns(
    reps: int,
    n: int = 390,
    sigma=0.02,
    p_jump: float = 0.01,
    mu_j: float = 0.0,
    sigma_j: float = 0.05,
    noise_var: float = 0.0,
    rng: np.random.Generator | None = None,
):
    """Returns from eq. (13) of Bapat (2026) on the per-step scale.

    ``sigma`` is the per-step diffusive standard deviation (scalar or length-n
    array). The paper does not report the value it used, see
    ``docs/estimator_issues.md`` for what the published tables imply.

    Returns ``(r, target)`` where ``target`` is the per-step quarticity
    ``mean(sigma_i^4)``. Multiply by ``n`` (``T = 1``) to get IQ on the RQ scale.
    """
    rng = np.random.default_rng() if rng is None else rng
    sig = np.broadcast_to(np.asarray(sigma, float), (n,))
    eps = rng.standard_normal((reps, n))
    jumps = (rng.random((reps, n)) < p_jump) * rng.normal(mu_j, sigma_j, (reps, n))
    noise = rng.normal(0.0, np.sqrt(noise_var), (reps, n)) if noise_var > 0 else 0.0
    r = sig * eps + jumps + noise
    return r, float(np.mean(sig**4))


@dataclass
class SVPaths:
    r: np.ndarray  # observed returns, shape (reps, n)
    iq: np.ndarray  # true integrated quarticity per path, shape (reps,)
    iv: np.ndarray  # true integrated variance per path
    spot_var: np.ndarray  # spot variance at each observation time, shape (reps, n)
    jump_count: np.ndarray  # number of jumps per path


def intraday_u_shape(n: int, amplitude: float = 0.6) -> np.ndarray:
    """Multiplicative U-shaped volatility profile with mean square 1 over the day."""
    x = np.linspace(0, 1, n)
    prof = 1 + amplitude * (np.cos(2 * np.pi * x) + 1) / 2
    return prof / np.sqrt(np.mean(prof**2))


def sv_jump_returns(
    reps: int,
    n: int = 390,
    T: float = 1.0,
    sigma_bar: float = 0.2,
    kappa: float = 5.0,
    xi: float = 1.0,
    jump_intensity: float = 0.0,
    jump_size_sd: float = 0.0,
    noise_sd: float = 0.0,
    noise_level: str = "price",
    seasonal_amplitude: float = 0.0,
    substeps: int = 10,
    rng: np.random.Generator | None = None,
) -> SVPaths:
    """Log-OU stochastic volatility with jumps and noise on ``[0, T]``.

    ``log sigma_t^2`` follows ``d x = -kappa (x - m) dt + xi dB``, started from
    its stationary law, and ``m`` is set so that ``E[sigma_t^2] = sigma_bar^2``.
    ``jump_intensity`` is the expected number of jumps on ``[0, T]``. Jump sizes are
    ``N(0, jump_size_sd^2)`` in log-price units. With ``noise_level="price"`` the
    observed log price is ``p + u`` with i.i.d. ``u ~ N(0, noise_sd^2)``, so returns
    pick up ``u_i - u_{i-1}``. With ``noise_level="return"`` i.i.d. noise is added to
    each return, which is eq. (13) of Bapat (2026).
    """
    rng = np.random.default_rng() if rng is None else rng
    N = n * substeps
    dt = T / N
    stat_var = xi**2 / (2 * kappa)
    m = np.log(sigma_bar**2) - stat_var / 2
    x = m + np.sqrt(stat_var) * rng.standard_normal(reps)
    a = np.exp(-kappa * dt)
    sd = np.sqrt(stat_var * (1 - a**2))
    season = intraday_u_shape(N, seasonal_amplitude) if seasonal_amplitude > 0 else np.ones(N)
    var_path = np.empty((reps, N))
    for k in range(N):
        var_path[:, k] = np.exp(x) * season[k] ** 2
        x = m + a * (x - m) + sd * rng.standard_normal(reps)
    dW = rng.standard_normal((reps, N)) * np.sqrt(dt)
    fine = np.sqrt(var_path) * dW
    r = fine.reshape(reps, n, substeps).sum(axis=2)
    iq = np.sum(var_path**2, axis=1) * dt
    iv = np.sum(var_path, axis=1) * dt
    spot = var_path.reshape(reps, n, substeps).mean(axis=2)
    counts = rng.poisson(jump_intensity, reps) if jump_intensity > 0 else np.zeros(reps, int)
    for i in np.nonzero(counts)[0]:
        idx = rng.integers(0, n, counts[i])
        np.add.at(r[i], idx, rng.normal(0, jump_size_sd, counts[i]))
    if noise_sd > 0:
        if noise_level == "price":
            u = rng.normal(0, noise_sd, (reps, n + 1))
            r = r + np.diff(u, axis=1)
        elif noise_level == "return":
            r = r + rng.normal(0, noise_sd, (reps, n))
        else:
            raise ValueError("noise_level must be 'price' or 'return'")
    return SVPaths(r=r, iq=iq, iv=iv, spot_var=spot, jump_count=counts)


@dataclass
class ChangePaths:
    r: np.ndarray  # returns, shape (reps, n)
    sigma: np.ndarray  # per-step diffusive sd, shape (reps, n)
    nu: int | None  # index of the first post-change return (None = no change)


def change_returns(
    reps: int,
    n: int,
    nu: int | None = None,
    ratio: float = 1.5,
    sigma0: float = 1.0,
    p_jump: float = 0.0,
    jump_scale: float = 8.0,
    noise_sd: float = 0.0,
    seasonal_amplitude: float = 0.0,
    period: int = 390,
    sv_xi: float = 0.0,
    sv_halflife: float = 390.0,
    heavy_tail_df: float | None = None,
    rng: np.random.Generator | None = None,
) -> ChangePaths:
    """Per-step returns with a volatility change by factor ``ratio`` at index ``nu``.

    ``sigma_i = sigma0 * season_i * exp(sv_i) * (ratio if i >= nu else 1)``.
    ``ratio`` multiplies the standard deviation, so the variance changes by
    ``ratio**2``. Jumps arrive with probability ``p_jump`` per step with size
    ``N(0, (jump_scale * sigma0)^2)``. ``sv_xi`` is the stationary sd of the slow
    log-volatility AR(1) factor (0 disables it). ``heavy_tail_df`` replaces the
    normal innovations by unit-variance Student t innovations.
    """
    rng = np.random.default_rng() if rng is None else rng
    sig = np.full((reps, n), float(sigma0))
    if seasonal_amplitude > 0:
        prof = intraday_u_shape(period, seasonal_amplitude)
        sig *= np.resize(prof, n)
    if sv_xi > 0:
        phi = 0.5 ** (1.0 / sv_halflife)
        z = np.empty((reps, n))
        z[:, 0] = sv_xi * rng.standard_normal(reps)
        e = rng.standard_normal((reps, n)) * sv_xi * np.sqrt(1 - phi**2)
        for k in range(1, n):
            z[:, k] = phi * z[:, k - 1] + e[:, k]
        sig *= np.exp(z - sv_xi**2)  # E[exp(2 z - 2 xi^2)] = 1, so the mean variance is unchanged
    if nu is not None:
        sig[:, nu:] *= ratio
    if heavy_tail_df is None:
        eps = rng.standard_normal((reps, n))
    else:
        df = heavy_tail_df
        eps = rng.standard_t(df, (reps, n)) * np.sqrt((df - 2) / df)
    r = sig * eps
    if p_jump > 0:
        r += (rng.random((reps, n)) < p_jump) * rng.normal(0, jump_scale * sigma0, (reps, n))
    if noise_sd > 0:
        u = rng.normal(0, noise_sd, (reps, n + 1))
        r += np.diff(u, axis=1)
    return ChangePaths(r=r, sigma=sig, nu=nu)
