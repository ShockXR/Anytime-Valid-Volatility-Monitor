"""Synthetic bivariate returns with a known correlation change, volatility change, jumps or fat tails.

All output is synthetic. ``pair_returns`` draws ``reps`` paths of ``n`` steps of
a pair ``(r1, r2)`` with

    r_k,t = sigma_k,t * eps_k,t + jump_k,t,

where ``(eps_1, eps_2)`` has correlation ``rho`` before step ``nu`` and ``rho1``
from step ``nu`` on, unit variances, and is Gaussian or multivariate Student t
(a common chi-square mixing variable, so the pair has tail dependence, as
crypto pairs do). ``vol1_post`` and ``vol2_post`` multiply the standard
deviation of each asset from step ``nu`` on, which is how a "volatility change
without a correlation change" is produced. Jumps arrive independently in each
asset with probability ``p_jump`` per step and size ``N(0, (jump_scale sd)^2)``;
``p_cojump`` adds simultaneous jumps in both assets with independent signs,
which leaves the sign statistic's mean unchanged (checked in the tests).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class PairPaths:
    r1: np.ndarray  # shape (reps, n)
    r2: np.ndarray
    nu: int | None  # first post-change index (None = no change)
    rho: float
    rho1: float | None


def pair_returns(
    reps: int,
    n: int,
    rho: float = 0.3,
    nu: int | None = None,
    rho1: float | None = None,
    vol1_post: float = 1.0,
    vol2_post: float = 1.0,
    df: float | None = None,
    p_jump: float = 0.0,
    jump_scale: float = 8.0,
    p_cojump: float = 0.0,
    sigma1: float = 1.0,
    sigma2: float = 1.0,
    rng: np.random.Generator | None = None,
) -> PairPaths:
    rng = np.random.default_rng() if rng is None else rng
    e1 = rng.standard_normal((reps, n))
    e2 = rng.standard_normal((reps, n))
    rho_t = np.full(n, float(rho))
    if nu is not None and rho1 is not None:
        rho_t[nu:] = rho1
    z1 = e1
    z2 = rho_t[None, :] * e1 + np.sqrt(1 - rho_t[None, :] ** 2) * e2
    if df is not None:
        if df <= 2:
            raise ValueError("df must exceed 2 for a finite variance")
        w = rng.chisquare(df, (reps, n)) / df
        mix = np.sqrt((df - 2) / df) / np.sqrt(w)  # unit variance Student t scaling, common to both assets
        z1 = z1 * mix
        z2 = z2 * mix
    s1 = np.full(n, float(sigma1))
    s2 = np.full(n, float(sigma2))
    if nu is not None:
        s1[nu:] *= vol1_post
        s2[nu:] *= vol2_post
    r1 = s1[None, :] * z1
    r2 = s2[None, :] * z2
    if p_jump > 0:
        r1 = r1 + (rng.random((reps, n)) < p_jump) * rng.normal(0, jump_scale * sigma1, (reps, n))
        r2 = r2 + (rng.random((reps, n)) < p_jump) * rng.normal(0, jump_scale * sigma2, (reps, n))
    if p_cojump > 0:
        hit = rng.random((reps, n)) < p_cojump
        r1 = r1 + hit * rng.normal(0, jump_scale * sigma1, (reps, n))
        r2 = r2 + hit * rng.normal(0, jump_scale * sigma2, (reps, n))
    return PairPaths(r1, r2, nu, float(rho), None if rho1 is None else float(rho1))
