"""Purely sequential fixed-accuracy estimation of a volatility level (alternative ALT6).

Returns ``r_1, r_2, ...`` arrive one per step. Each procedure decides, after every
step, whether it has seen enough data to report an estimate whose error is at most
``d`` (absolute accuracy) or ``d * theta`` (relative accuracy), with probability at
least ``1 - alpha`` at its stopping time ``N``.

Statistics
----------
* raw: ``X = r^2``, target ``theta = E r^2`` (the variance).
* truncated: ``X = min(r^2, tau^2)`` with ``tau = c * scaled MAD`` of a pilot sample
  of ``m`` returns, target ``theta_tau = E[min(r^2, tau^2) | tau]`` (the truncated
  variance at the realised level). Its second moment ``E[min(r^4, tau^4)]`` is a
  truncated quarticity: the averaged robust kernel estimator of Bapat (2026) up to a
  constant (``docs/estimator_issues.md``, Issue 3).

Procedures (vectorised over independent paths, streamed block by block)
-----------------------------------------------------------------------
* ``CR``: Chow and Robbins (1965), ``N = inf{n >= m : n >= (z/d)^2 (nu_n + 1/n)}`` with
  the sample variance ``nu_n``; relative version (Nadas 1967)
  ``n >= (z/d)^2 (nu_n + 1/n) / Xbar_n^2``. Estimate ``Xbar_N``.
* ``CR-BM``: the same rule with a batch-means long-run variance (batches of fixed
  length, rule checked at batch ends): a fix for serial dependence.
* ``Stein``: Stein (1945) two-stage, ``N = max(m, floor((t_{m-1}/d)^2 s_m^2) + 1)``.
* ``Oracle``: a fixed sample of ``ceil(n*)``, ``n* = z^2 nu / d^2`` with the TRUE nu (a
  reference, not a procedure; nu is the long-run variance for dependent laws).
* ``MinRisk``: Robbins (1959) minimum risk point estimation for the loss
  ``A (Xbar_n - theta)^2 + n``: ``N = inf{n >= m : n >= sqrt(A) (s_n + n^(-gamma))}``.
* ``CS-*``: an anytime-valid confidence sequence for the mean of ``Z = X / tau^2`` in
  ``[0, 1]``, started after the pilot (so ``tau`` is fixed before it starts) and stopped
  the first time its width reaches the target. ``CS-wsr``, ``CS-tuned`` and ``CS-opt``
  are the predictable plug-in empirical Bernstein sequence of Waudby-Smith and Ramdas
  (2024) with the one-sided ranges of ``avmonitor.confseq`` and three bet schedules;
  ``CS-bet`` is a betting (capital process) sequence evaluated on a grid of candidate
  means with an exact concavity bound between grid points.

Every procedure runs on the same simulated paths (common random numbers). Closed-form
predictions for the confidence sequences: ``cs_ratio_limit``, ``cs_n_predict``,
``betting_oracle_n``. Derivations are in this folder's ``README.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize, stats

from avmonitor.confseq import psi_coef
from avmonitor.robust import MAD_NORMAL_SCALE

# ----------------------------------------------------------------------------
# Constants and closed forms
# ----------------------------------------------------------------------------


def zq(alpha: float) -> float:
    """Two-sided normal quantile ``z_{1 - alpha/2}``."""
    return float(stats.norm.ppf(1.0 - alpha / 2.0))


def n_star(nu, d, alpha: float, theta=None):
    """Optimal fixed sample size ``z^2 nu / d^2`` (or ``z^2 nu / (d theta)^2`` if theta is given)."""
    nu = np.asarray(nu, float)
    dd = d if theta is None else d * np.asarray(theta, float)
    return zq(alpha) ** 2 * nu / dd**2


def cs_ratio_limit(alpha: float) -> float:
    """Limit of ``E[N_cs] / n*`` as ``d -> 0`` with bets tuned to the target: ``2 log(2/alpha) / z^2``."""
    return 2.0 * math.log(2.0 / alpha) / zq(alpha) ** 2


def cs_limit_coverage(alpha: float) -> float:
    """Normal-approximation coverage at the stopping time as ``d -> 0``: ``1 - 2 (1 - Phi(sqrt(2 log(2/alpha))))``."""
    return 1.0 - 2.0 * stats.norm.sf(math.sqrt(2.0 * math.log(2.0 / alpha)))


def brownian_limit(alpha: float, kappa: float = 1.0, intersect: bool = True, reps: int = 20_000,
                   du: float = 1e-3, umax: float = 4.0, seed: int = 0, chunk: int = 2000) -> dict:
    """Limit of ``E[N_cs] / n*`` as ``d -> 0`` for a confidence sequence intersected over time.

    With constant bets ``lam = kappa d / var`` the bounds in units of d, at time
    ``u n*`` (``n* = z^2 var / d^2``), converge to
    ``U(u) = W(u) / (z u) + L / (kappa z^2 u) + kappa / 2`` and
    ``L(u) = W(u) / (z u) - L / (kappa z^2 u) - kappa / 2`` with ``L = log(2/alpha)`` and W a
    standard Brownian motion (centred at the true mean). The intersected sequence stops at
    ``u* = inf{u : min_{v<=u} U(v) - max_{v<=u} L(v) <= 2}``. Without intersection
    ``u* = 2 L / (kappa (2 - kappa) z^2)`` (``= 2 log(2/alpha) / z^2`` at kappa = 1). Returns
    the Monte Carlo mean of ``u*``, its standard error and the coverage at ``u*``.
    """
    z = zq(alpha)
    L = math.log(2.0 / alpha)
    if not intersect:
        return dict(ratio=2 * L / (kappa * (2 - kappa) * z * z), se=0.0, coverage=float("nan"))
    rng = np.random.default_rng(seed)
    nst = int(round(umax / du))
    u = du * np.arange(1, nst + 1)
    us, cov = [], []
    for _ in range(max(1, reps // chunk)):
        W = np.cumsum(rng.standard_normal((chunk, nst)) * math.sqrt(du), axis=1)
        Um = np.minimum.accumulate(W / (z * u) + L / (kappa * z * z * u) + kappa / 2, axis=1)
        Lm = np.maximum.accumulate(W / (z * u) - L / (kappa * z * z * u) - kappa / 2, axis=1)
        ok = (Um - Lm) <= 2.0
        hit = ok.any(axis=1)
        idx = ok.argmax(axis=1)
        us.append(np.where(hit, u[idx], np.inf))
        r = np.arange(chunk)
        cov.append((Lm[r, idx] <= 0) & (Um[r, idx] >= 0))
    us = np.concatenate(us)
    return dict(ratio=float(us.mean()), se=float(us.std() / math.sqrt(us.size)),
                coverage=float(np.concatenate(cov).mean()), not_stopped=float(np.isinf(us).mean()))


def kappa2(a):
    """``E[min(Z^2, a^2)]`` for standard normal Z (vectorised)."""
    a = np.asarray(a, float)
    phi = stats.norm.pdf(a)
    tail = 2.0 * stats.norm.sf(a)
    return (1.0 - tail) - 2.0 * a * phi + a * a * tail


def kappa4(a):
    """``E[min(Z^4, a^4)]`` for standard normal Z (vectorised); equals ``3 * truncated_moment_factor(a, 4)``."""
    a = np.asarray(a, float)
    phi = stats.norm.pdf(a)
    tail = 2.0 * stats.norm.sf(a)
    c2 = (1.0 - tail) - 2.0 * a * phi
    return 3.0 * c2 - 2.0 * a**3 * phi + a**4 * tail


def opt_plan(mean_z, var_z, d, alpha: float, cap: float = 0.75, relative: bool = False, d_z=None,
             grid: int = 96, iters: int = 48):
    """Constant bets per side that minimise the predicted stopping time (vectorised over paths).

    For a candidate count n each side's predicted half-width
    ``w_s(lam, n) = log(2/alpha) / (lam n) + lam h(lam a_s) var_z`` is minimised over a
    log grid of bets from 1/30 to 30 times ``d_z / var_z``, capped at ``cap / a_s``; n is then found by bisection as the
    smallest count at which the stopping rule fires. Returns ``(n_hat, lam_up, lam_lo)``
    (``n_hat = inf`` where the rule cannot fire). Absolute rule: ``d_z`` is the target
    half-width in Z units. Relative rule: ``d`` is the relative half-width.
    """
    mz = np.clip(np.atleast_1d(np.asarray(mean_z, float)), 1e-6, 1 - 1e-6)
    vz = np.maximum(np.atleast_1d(np.asarray(var_z, float)), 1e-10)
    logd = math.log(2.0 / alpha)
    a_up, a_lo = 1.0 - mz, mz
    if relative:
        wu, wl, rhs = 1.0 - d, 1.0 + d, 2.0 * d * mz
    else:
        wu, wl = 1.0, 1.0
        rhs = 2.0 * np.broadcast_to(np.asarray(d_z, float), mz.shape)
    # bets on a log grid around the quadratic optimum d_z / var_z, capped at cap / a
    base = (d * mz if relative else rhs / 2.0) / vz
    f = np.geomspace(1.0 / 30.0, 30.0, grid)
    lu = np.minimum(f[None, :] * base[:, None], (cap / a_up)[:, None])
    ll = np.minimum(f[None, :] * base[:, None], (cap / a_lo)[:, None])
    gu = lu * psi_coef(lu * a_up[:, None]) * vz[:, None]
    gl = ll * psi_coef(ll * a_lo[:, None]) * vz[:, None]

    def total(n):
        wU = logd / (lu * n[:, None]) + gu
        wL = logd / (ll * n[:, None]) + gl
        return wu * wU.min(axis=1) + wl * wL.min(axis=1), wU.argmin(axis=1), wL.argmin(axis=1)

    lo = np.ones_like(mz)
    hi = np.full_like(mz, 1e9)
    never = total(hi)[0] > rhs
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        ok = total(mid)[0] <= rhs
        hi = np.where(ok, mid, hi)
        lo = np.where(ok, lo, mid)
    _, iu, il = total(hi)
    rows = np.arange(mz.size)
    n_hat = np.where(never, np.inf, hi)
    return n_hat, lu[rows, iu], ll[rows, il]


def cs_n_predict(mean_z: float, var_z: float, d_z: float, alpha: float, cap: float = 0.75,
                 bets: str = "tuned", relative: bool = False, d_rel: float | None = None) -> float:
    """Predicted number of sequence observations until the empirical Bernstein rule fires.

    Replaces every running estimate by its limit (``zhat -> mean_z``, ``s2 -> var_z``,
    mean of ``e -> var_z``) and keeps the exact bet cap and the exact ``psi`` term, so
    it is the first-order prediction at a practical ``d``. Bets: ``"tuned"``
    (``lam_s = min(d_z / var_z, cap / a_s)``) or ``"opt"`` (``opt_plan``). Absolute rule:
    ``w_up + w_lo <= 2 d_z``. Relative rule (``d_z = d_rel * mean_z``):
    ``(1 - d) w_up + (1 + d) w_lo <= 2 d mean_z``, the linearisation of
    ``(U - L) / (U + L) <= d``. Returns ``inf`` if the rule can never fire.
    """
    logd = math.log(2.0 / alpha)
    a_up, a_lo = 1.0 - mean_z, mean_z
    if bets == "opt":
        n_hat, _, _ = opt_plan(mean_z, var_z, d_rel, alpha, cap, relative, d_z=d_z)
        return float(n_hat[0])
    if bets != "tuned":
        raise ValueError("bets must be 'tuned' or 'opt'")
    wu, wl = (1.0, 1.0) if not relative else (1.0 - d_rel, 1.0 + d_rel)
    rhs = 2.0 * d_z if not relative else 2.0 * d_rel * mean_z
    lu = min(d_z / var_z, cap / a_up)
    ll = min(d_z / var_z, cap / a_lo)
    num = logd * (wu / lu + wl / ll)
    den = rhs - var_z * (wu * lu * float(psi_coef(lu * a_up)) + wl * ll * float(psi_coef(ll * a_lo)))
    return num / den if den > 0 else math.inf


def betting_oracle_n(z_atoms, p, d, alpha: float, relative: bool = False, d_z=None) -> float:
    """Stopping count of the best constant-bet betting sequence for a KNOWN law of Z (a benchmark).

    ``z_atoms``, ``p``: a discretised law of Z on [0, 1]. Side s rejects the candidate
    mean ``mu +- w`` after ``log(2/alpha) / g_s(w)`` steps, where
    ``g_up(w) = max_lam E log(1 + lam (mu + w - Z))`` and
    ``g_lo(w) = max_lam E log(1 + lam (Z - mu + w))`` (the KL-inf rates of a growth-rate
    optimal bettor). Returns the smallest n at which the two half-widths satisfy the
    stopping rule (absolute: ``w_up + w_lo <= 2 d_z``; relative:
    ``(1 - d) w_up + (1 + d) w_lo <= 2 d mu``).
    """
    z = np.asarray(z_atoms, float)
    p = np.asarray(p, float) / np.sum(p)
    mu = float(np.sum(p * z))
    L = math.log(2.0 / alpha)

    def g(w, side):
        if side == "up":
            m = mu + w
            f = lambda lam: -np.sum(p * np.log1p(lam * (m - z)))  # noqa: E731
            ub = 1.0 / (1.0 - m)
        else:
            m = mu - w
            f = lambda lam: -np.sum(p * np.log1p(lam * (z - m)))  # noqa: E731
            ub = 1.0 / m
        r = optimize.minimize_scalar(f, bounds=(0.0, ub * (1 - 1e-9)), method="bounded",
                                     options={"xatol": 1e-10})
        return -r.fun

    def width(n, side):
        hi = (1 - mu) * 0.999 if side == "up" else mu * 0.999
        f = lambda w: n * g(w, side) - L  # noqa: E731
        if f(hi) < 0:
            return math.inf
        return optimize.brentq(f, 1e-12, hi, xtol=1e-12)

    def ok(n):
        wu, wl = width(n, "up"), width(n, "lo")
        if not (math.isfinite(wu) and math.isfinite(wl)):
            return False
        if relative:
            return (1 - d) * wu + (1 + d) * wl <= 2 * d * mu
        return wu + wl <= 2 * d_z

    lo, hi = 1.0, 1e9
    if not ok(hi):
        return math.inf
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if ok(mid):
            hi = mid
        else:
            lo = mid
        if hi / lo < 1 + 1e-6:
            break
    return hi


# ----------------------------------------------------------------------------
# Return laws with exact truncated moments
# ----------------------------------------------------------------------------

_GL_X, _GL_W = np.polynomial.legendre.leggauss(400)


class Law:
    """Base class. ``draw(rng, state, rows, k)`` returns k new returns for the paths in ``rows``."""

    name = "law"
    iid = True

    def start(self, rng, reps):
        return None

    def draw(self, rng, state, rows, k):
        """Return ``(r, cond_var, state)`` for the rows ``rows`` (``cond_var`` may be None)."""
        raise NotImplementedError

    def tail(self, t):
        """``P(|r| > t)``."""
        raise NotImplementedError

    def theta_raw(self) -> float:
        return 1.0

    def m4_raw(self) -> float:
        raise NotImplementedError

    def nu_raw(self) -> float:
        return self.m4_raw() - self.theta_raw() ** 2

    def trunc_moments(self, tau):
        """``(E min(r^2, tau^2), E min(r^4, tau^4))`` from ``E min(Y, b) = int_0^b P(Y > y) dy``."""
        tau = np.atleast_1d(np.asarray(tau, float))
        x = 0.5 * (_GL_X + 1.0)
        w = 0.5 * _GL_W
        t = tau[:, None] * x[None, :]
        p = self.tail(t)
        m2 = (w[None, :] * 2.0 * t * p).sum(axis=1) * tau
        m4 = (w[None, :] * 4.0 * t**3 * p).sum(axis=1) * tau
        return m2, m4

    def median_abs(self) -> float:
        return float(optimize.brentq(lambda t: self.tail(np.array([t]))[0] - 0.5, 1e-9, 100.0))

    def tau_limit(self, c: float = 3.0) -> float:
        """Large-pilot limit of ``tau = c * scaled MAD about zero``."""
        return c * MAD_NORMAL_SCALE * self.median_abs()

    def z_law(self, tau: float, k: int = 4000):
        """Discretised law of ``Z = min(r^2, tau^2) / tau^2`` (atoms, probabilities), for i.i.d. laws."""
        x = np.linspace(0.0, tau, k + 1)
        mid = 0.5 * (x[1:] + x[:-1])
        surv = self.tail(x)
        p = np.append(surv[:-1] - surv[1:], surv[-1])
        z = np.append(mid**2 / tau**2, 1.0)
        return z, p / p.sum()


class Gaussian(Law):
    name = "gauss"

    def __init__(self, sigma: float = 1.0):
        self.sigma = float(sigma)

    def draw(self, rng, state, rows, k):
        return self.sigma * rng.standard_normal((len(rows), k)), None, state

    def tail(self, t):
        return 2.0 * stats.norm.sf(np.asarray(t) / self.sigma)

    def theta_raw(self):
        return self.sigma**2

    def m4_raw(self):
        return 3.0 * self.sigma**4

    def trunc_moments(self, tau):
        a = np.atleast_1d(np.asarray(tau, float)) / self.sigma
        return self.sigma**2 * kappa2(a), self.sigma**4 * kappa4(a)


class StudentT(Law):
    """Unit-variance Student t with ``df`` degrees of freedom (df > 2)."""

    def __init__(self, df: float):
        self.df = float(df)
        self.s = math.sqrt((df - 2.0) / df)
        self.name = f"t{df:g}"

    def draw(self, rng, state, rows, k):
        return rng.standard_t(self.df, (len(rows), k)) * self.s, None, state

    def tail(self, t):
        return 2.0 * stats.t.sf(np.asarray(t) / self.s, self.df)

    def m4_raw(self):
        if self.df <= 4:
            return math.inf
        return 3.0 * (self.df - 2.0) / (self.df - 4.0)


class Jumps(Law):
    """Gaussian diffusion of sd ``sigma`` plus Bernoulli(p) jumps of size N(0, sj^2) (Bapat 2026, eq. 10)."""

    def __init__(self, p: float = 0.01, sj: float = 5.0, sigma: float = 1.0):
        self.p, self.sj, self.sigma = float(p), float(sj), float(sigma)
        self.name = f"jumps(p={p:g},sJ={sj:g})"
        self.sig_j = math.sqrt(sigma**2 + sj**2)

    def draw(self, rng, state, rows, k):
        n = len(rows)
        r = self.sigma * rng.standard_normal((n, k))
        r += (rng.random((n, k)) < self.p) * rng.normal(0.0, self.sj, (n, k))
        return r, None, state

    def tail(self, t):
        t = np.asarray(t)
        return (1 - self.p) * 2 * stats.norm.sf(t / self.sigma) + self.p * 2 * stats.norm.sf(t / self.sig_j)

    def theta_raw(self):
        return self.sigma**2 + self.p * self.sj**2

    def m4_raw(self):
        return 3.0 * ((1 - self.p) * self.sigma**4 + self.p * self.sig_j**4)

    def trunc_moments(self, tau):
        tau = np.atleast_1d(np.asarray(tau, float))
        q, s1, s2 = self.p, self.sigma, self.sig_j
        m2 = (1 - q) * s1**2 * kappa2(tau / s1) + q * s2**2 * kappa2(tau / s2)
        m4 = (1 - q) * s1**4 * kappa4(tau / s1) + q * s2**4 * kappa4(tau / s2)
        return m2, m4


class Garch(Law):
    """GARCH(1, 1) with Gaussian innovations and unconditional variance 1 (Bollerslev 1986)."""

    iid = False

    def __init__(self, a: float = 0.10, b: float = 0.85, burn: int = 3000):
        self.a, self.b, self.burn = float(a), float(b), int(burn)
        self.omega = 1.0 - a - b
        self.name = f"garch(a={a:g},b={b:g})"

    def start(self, rng, reps):
        h = np.ones(reps)
        r = np.zeros(reps)
        for _ in range(self.burn):
            h = self.omega + self.a * r * r + self.b * h
            r = np.sqrt(h) * rng.standard_normal(reps)
        return [h, r]

    def draw(self, rng, state, rows, k):
        h, r = state[0][rows], state[1][rows]
        out = np.empty((len(rows), k))
        hv = np.empty((len(rows), k))
        z = rng.standard_normal((len(rows), k))
        for j in range(k):
            h = self.omega + self.a * r * r + self.b * h
            r = np.sqrt(h) * z[:, j]
            out[:, j] = r
            hv[:, j] = h
        state[0][rows], state[1][rows] = h, r
        return out, hv, state

    def kurtosis(self):
        """``E r^4 / (E r^2)^2 = 3 (1 - (a + b)^2) / (1 - (a + b)^2 - 2 a^2)``."""
        a, b = self.a, self.b
        return 3.0 * (1 - (a + b) ** 2) / (1 - (a + b) ** 2 - 2 * a * a)

    def nu_raw(self):
        """``Var(r^2) = E r^4 - 1`` (finite if ``3 a^2 + 2 a b + b^2 < 1``)."""
        a, b = self.a, self.b
        if 3 * a * a + 2 * a * b + b * b >= 1:
            return math.inf
        return self.kurtosis() - 1.0

    def rho1(self):
        """Lag-one autocorrelation of ``r^2``: ``a (1 - a b - b^2) / (1 - 2 a b - b^2)``; lag k: ``rho1 (a + b)^(k-1)``."""
        a, b = self.a, self.b
        return a * (1 - a * b - b * b) / (1 - 2 * a * b - b * b)

    def lrv_factor_raw(self):
        """Long-run variance of ``r^2`` divided by its variance: ``1 + 2 rho1 / (1 - a - b)``."""
        return 1.0 + 2.0 * self.rho1() / (1.0 - self.a - self.b)


class Bootstrap(Law):
    """Stationary bootstrap (Politis and Romano 1994) of one series, mean block length ``block``."""

    def __init__(self, u: np.ndarray, block: float = 1.0, name: str = "boot"):
        u = np.asarray(u, float)
        self.u = u[np.isfinite(u)]
        self.block = float(block)
        self.name = name

    def start(self, rng, reps):
        return rng.integers(0, self.u.size, reps)

    def draw(self, rng, state, rows, k):
        n = self.u.size
        nr = len(rows)
        if self.block <= 1.0:
            ii = rng.integers(0, n, (nr, k))
            state[rows] = ii[:, -1]
            return self.u[ii], None, state
        idx = state[rows]
        out = np.empty((nr, k))
        jump = rng.random((nr, k)) < 1.0 / self.block
        fresh = rng.integers(0, n, (nr, k))
        for j in range(k):
            idx = np.where(jump[:, j], fresh[:, j], (idx + 1) % n)
            out[:, j] = self.u[idx]
        state[rows] = idx
        return out, None, state

    def theta_raw(self):
        return float(np.mean(self.u**2))

    def m4_raw(self):
        return float(np.mean(self.u**4))

    def trunc_moments(self, tau):
        tau = np.atleast_1d(np.asarray(tau, float))
        u2 = np.sort(self.u**2)
        cs = np.concatenate([[0.0], np.cumsum(u2)])
        cs4 = np.concatenate([[0.0], np.cumsum(u2**2)])
        n = u2.size
        k = np.searchsorted(u2, tau**2, side="right")
        m2 = (cs[k] + (n - k) * tau**2) / n
        m4 = (cs4[k] + (n - k) * tau**4) / n
        return m2, m4


class Windows(Law):
    """Rows of a fixed matrix of returns (real data): path i reads row i from left to right."""

    iid = False

    def __init__(self, R: np.ndarray, name: str = "real"):
        self.R = np.asarray(R, float)
        self.name = name

    def start(self, rng, reps):
        if reps != self.R.shape[0]:
            raise ValueError("reps must equal the number of rows")
        return [0]

    def draw(self, rng, state, rows, k):
        j = state[0]
        if j + k > self.R.shape[1]:
            raise ValueError("ran out of data; lower n_max")
        state[0] = j + k
        return self.R[rows, j: j + k], None, state


# ----------------------------------------------------------------------------
# The streaming engine
# ----------------------------------------------------------------------------


@dataclass
class Setup:
    d: float  # accuracy: absolute (units of theta) or relative
    alpha: float = 0.05
    m: int = 30  # pilot size
    c: float = 3.0  # truncation constant
    relative: bool = False
    n_max: int = 200_000
    block: int = 128
    cs_cap: float = 0.75
    batch: int = 100  # batch length for CR-BM
    minrisk_A: float | None = None  # loss A (Xbar - theta)^2 + n, in the units of the returns
    minrisk_gamma: float = 1.0
    bet_grid: int = 96  # grid points of CS-bet
    bet_span: float = 5.0  # CS-bet grid covers [zhat0 / span, zhat0 * span]
    intersect: bool = True  # CS bounds intersected over time (False reproduces avmonitor.confseq exactly)
    procs: tuple = ("CR-raw", "CR-trunc", "Stein-raw", "Stein-trunc", "Oracle-raw", "Oracle-trunc",
                    "CS-opt")
    extra: dict = field(default_factory=dict)


@dataclass
class Result:
    N: np.ndarray  # stopping times (total observations, pilot included)
    est: np.ndarray  # point estimate at N (original units)
    lo: np.ndarray  # interval at N (original units)
    hi: np.ndarray
    target: np.ndarray  # theta (raw) or theta_tau (trunc) per path, original units
    capped: np.ndarray  # True if the rule had not fired by n_max
    nstar: np.ndarray  # oracle fixed sample size per path (inf if undefined)
    cond_target: np.ndarray | None = None  # average conditional mean over the sample (dependent laws)
    cond_target_w: np.ndarray | None = None  # CS: (lambda-up weighted, lambda-down weighted) averages
    loss: np.ndarray | None = None  # MinRisk loss A (est - target)^2 + N
    var_at: np.ndarray | None = None  # sample variance of X at N (original units squared)


def _first_hit(cond):
    anyh = cond.any(axis=1)
    return anyh, np.where(anyh, cond.argmax(axis=1), -1)


class _Proc:
    """Common bookkeeping: full-size state, updates on the rows that are still running."""

    def __init__(self, reps):
        self.active = np.ones(reps, bool)

    def sel(self, rows):
        """Positions within ``rows`` that are still active, and their global indices."""
        mask = self.active[rows]
        return np.nonzero(mask)[0], rows[mask]


class _MeanVarRule(_Proc):
    """Chow-Robbins, Nadas, batch means, Robbins minimum risk and the fixed oracle on running sums of X."""

    def __init__(self, kind, rule, setup, x_pilot, nstar, cond_pilot=None, A=None):
        reps, m = x_pilot.shape
        super().__init__(reps)
        self.kind, self.rule, self.s = kind, rule, setup
        self.K = x_pilot.mean(axis=1)  # shift for numerical stability
        y = x_pilot - self.K[:, None]
        self.n = m
        self.S1 = y.sum(axis=1)
        self.S2 = (y * y).sum(axis=1)
        self.C = cond_pilot.sum(axis=1) if cond_pilot is not None else None
        self.N = np.full(reps, -1)
        self.mean_at = np.full(reps, np.nan)
        self.var_at = np.full(reps, np.nan)
        self.cond_at = np.full(reps, np.nan)
        self.nstar, self.A = np.asarray(nstar, float), A
        self.z = zq(setup.alpha)
        if rule == "bm":
            B = setup.batch
            nb = m // B
            self.bm_s1 = np.zeros(reps)
            self.bm_s2 = np.zeros(reps)
            for j in range(nb):
                bmean = x_pilot[:, j * B:(j + 1) * B].mean(axis=1)
                self.bm_s1 += bmean
                self.bm_s2 += bmean * bmean
            self.bm_n = nb
            self.bsum = x_pilot[:, nb * B:].sum(axis=1)
            self.bk = m - nb * B
        else:
            n = float(m)
            mean = self.K + self.S1 / n
            var = (self.S2 - self.S1**2 / n) / (n - 1)
            hit = self._cond(np.full(reps, n), mean, var, np.arange(reps))
            self._record(np.nonzero(hit)[0], np.full(hit.sum(), m), mean[hit], var[hit],
                         None if self.C is None else (self.C / n)[hit])

    def _cond(self, n, mean, var, g):
        s, z = self.s, self.z
        if self.rule in ("cr", "bm"):
            with np.errstate(divide="ignore", invalid="ignore"):
                if s.relative:
                    return n >= (z / s.d) ** 2 * (var + 1.0 / n) / mean**2
                return n >= (z / s.d) ** 2 * (var + 1.0 / n)
        if self.rule == "minrisk":
            rA = np.sqrt(self.A[g]).reshape((-1,) + (1,) * (np.ndim(n) - 1))
            return n >= rA * (np.sqrt(np.maximum(var, 0.0)) + n ** (-s.minrisk_gamma))
        if self.rule == "oracle":
            ns = np.maximum(np.ceil(self.nstar[g]), s.m).reshape((-1,) + (1,) * (np.ndim(n) - 1))
            return n >= ns
        raise ValueError(self.rule)

    def _record(self, g, N, mean, var, cond):
        self.N[g] = N
        self.mean_at[g] = mean
        self.var_at[g] = var
        if cond is not None:
            self.cond_at[g] = cond
        self.active[g] = False

    def update(self, x, rows, cond=None):
        pos, g = self.sel(rows)
        k = x.shape[1]
        if g.size == 0:
            self.n += k
            return
        x = x[pos]
        y = x - self.K[g][:, None]
        n = self.n + np.arange(1, k + 1, dtype=float)[None, :]
        s1 = self.S1[g][:, None] + np.cumsum(y, axis=1)
        s2 = self.S2[g][:, None] + np.cumsum(y * y, axis=1)
        cc = None if cond is None else self.C[g][:, None] + np.cumsum(cond[pos], axis=1)
        if self.rule == "bm":
            B = self.s.batch
            hitc = np.zeros((g.size, k), bool)
            var = np.full((g.size, k), np.nan)
            bsum, bk = self.bsum[g].copy(), self.bk
            bm1, bm2, bmn = self.bm_s1[g].copy(), self.bm_s2[g].copy(), self.bm_n
            p0 = 0
            while p0 < k:
                take = min(B - bk, k - p0)
                bsum += x[:, p0:p0 + take].sum(axis=1)
                bk += take
                p0 += take
                if bk == B:
                    bmean = bsum / B
                    bm1 += bmean
                    bm2 += bmean * bmean
                    bmn += 1
                    bsum = np.zeros(g.size)
                    bk = 0
                    if bmn >= 10:
                        v = B * (bm2 - bm1**2 / bmn) / (bmn - 1)
                        j = p0 - 1
                        var[:, j] = v
                        nn = float(self.n + p0)
                        mean_j = self.K[g] + s1[:, j] / nn
                        hitc[:, j] = self._cond(np.full(g.size, nn), mean_j, v, g)
            self.bsum[g], self.bm_s1[g], self.bm_s2[g] = bsum, bm1, bm2
            # bk and bm_n advance identically for every row (time is common)
            self._bk_next, self._bmn_next = bk, bmn
            anyh, idx = _first_hit(hitc)
        else:
            mean = self.K[g][:, None] + s1 / n
            var = (s2 - s1**2 / n) / (n - 1)
            anyh, idx = _first_hit(self._cond(n, mean, var, g))
        h = np.nonzero(anyh)[0]
        j = idx[h]
        nn = self.n + j + 1
        self._record(g[h], nn, self.K[g[h]] + s1[h, j] / nn, var[h, j],
                     None if cc is None else cc[h, j] / nn)
        self.S1[g] = s1[:, -1]
        self.S2[g] = s2[:, -1]
        if cc is not None:
            self.C[g] = cc[:, -1]
        self.n += k
        if self.rule == "bm":
            self.bk, self.bm_n = self._bk_next, self._bmn_next

    def advance(self, k):
        """Time passes for rows that are no longer drawn (all inactive)."""
        self.n += k

    def finish(self):
        a = self.active.copy()
        n = float(self.n)
        g = np.nonzero(a)[0]
        if g.size:
            self.N[g] = self.n
            self.mean_at[g] = self.K[g] + self.S1[g] / n
            self.var_at[g] = (self.S2[g] - self.S1[g] ** 2 / n) / (n - 1)
            if self.C is not None:
                self.cond_at[g] = self.C[g] / n
        return a


class _CS(_Proc):
    """Confidence sequence on Z in [0, 1] (started after the pilot), stopped at the target width.

    ``bets`` in ``{"wsr", "tuned", "opt"}``: predictable plug-in empirical Bernstein
    bounds (``docs/math.md`` 6.2) with one-sided ranges. ``bets="bet"``: capital
    processes ``prod (1 + lam (m - Z))`` (upper) and ``prod (1 + lam (Z - m))`` (lower)
    on a grid of candidate means ``m`` fixed from the pilot (a fine log grid around the
    pilot mean plus coarse extensions to 1e-6 and 1 - 1e-6), with one bet per side shared
    by all grid points (so log-capital is nondecreasing, resp. nonincreasing, and concave
    in m, which makes interpolation between grid points a valid bound; above the fine grid
    the lower bet is capped per point). Bounds are intersected over time.

    ``init="pilot"``: the running centre and variance start from the pilot Z values with
    weight ``m`` (known when the sequence starts, so the bets stay predictable).
    ``init="wsr"``: centre 1/2 and variance 1/4 with weight 1 (Waudby-Smith and Ramdas),
    which reproduces ``avmonitor.confseq.eb_confseq`` exactly.
    """

    def __init__(self, setup, z_pilot, bets="opt", init="pilot", n_plan=None, dz=None):
        reps = z_pilot.shape[0]
        super().__init__(reps)
        s = self.s = setup
        self.bets, self.n_plan = bets, n_plan
        self.dz = dz  # absolute target half-width in Z units (per row) or None (relative)
        self.logd = math.log(2.0 / s.alpha)  # delta_up = delta_lo = alpha / 2
        if init == "pilot":
            self.w0 = float(z_pilot.shape[1])
            self.z0 = np.clip(z_pilot.mean(axis=1), 1e-4, 1 - 1e-4)
            self.v0 = np.maximum(z_pilot.var(axis=1, ddof=1), 1e-6)
        elif init == "wsr":
            self.w0, self.z0, self.v0 = 1.0, np.full(reps, 0.5), np.full(reps, 0.25)
        else:
            raise ValueError(init)
        self.t = 0
        self.sum_z = np.zeros(reps)
        self.sum_dev2 = np.zeros(reps)
        self.zhat = self.z0.copy()
        self.s2 = self.v0.copy()
        self.sl = np.zeros((2, reps))
        self.slz = np.zeros((2, reps))
        self.spe = np.zeros((2, reps))
        self.slm = np.zeros((2, reps))
        self.sm = np.zeros(reps)
        self.U = np.ones(reps)
        self.L = np.zeros(reps)
        self.T = np.full(reps, -1)
        self.L_at = np.full(reps, np.nan)
        self.U_at = np.full(reps, np.nan)
        self.cm_at = np.full((3, reps), np.nan)
        if bets == "bet":
            # core: log grid on [zhat0 / span, span zhat0], where the bounds normally lie; between grid points
            # a second-order bound is used. Extensions: 8 log-spaced points down to 1e-6 and 32 log-spaced
            # points up to 1 - 1e-6, so that a pilot far from the truth cannot leave the true mean outside the
            # grid (the sequence would then never shrink). The bets are chosen for the core; above the core
            # the lower side's bet is capped per point at cs_cap / m (its factors stay positive and its log
            # capital stays nonincreasing in m, but it is no longer concave there, so the bound there is the
            # largest rejected grid point, without interpolation).
            G, NL, NH = s.bet_grid, 8, 32
            f = np.geomspace(1.0 / s.bet_span, s.bet_span, G)
            core = np.clip(self.z0[:, None] * f[None, :], 1e-6, 1 - 1e-6)
            tl = np.arange(NL)[None, :] / NL
            low = np.exp(np.log(1e-6) + tl * (np.log(core[:, :1]) - np.log(1e-6)))
            th = np.arange(1, NH + 1)[None, :] / NH
            high = np.exp(np.log(core[:, -1:]) + th * (np.log(1 - 1e-6) - np.log(core[:, -1:])))
            self.grid = np.concatenate([low, core, high], axis=1)
            self.core_lo, self.core_hi = NL, NL + G - 1  # indices of the first and last core points
            Gt = self.grid.shape[1]
            # per grid point: log capital F, slope sum D = sum lam / x, curvature sum C = sum (lam / x)^2,
            # x = 1 + lam (m - Z) (upper) or 1 + lam (Z - m) (lower)
            self.logK_up = np.zeros((reps, Gt))
            self.logK_lo = np.zeros((reps, Gt))
            self.D_up = np.zeros((reps, Gt))
            self.D_lo = np.zeros((reps, Gt))
            self.C_up = np.zeros((reps, Gt))
            self.C_lo = np.zeros((reps, Gt))
            # upper side: lam (1 - m) <= cs_cap / (1 - core_0) < 1 for every m in [0, 1) since core_0 <= 0.2
            self.cap_up = s.cs_cap / (1.0 - core[:, 0])
            # lower side: lam m <= cs_cap on the core; above it the per-point cap cs_cap / m applies
            self.cap_lo = s.cs_cap / core[:, -1]
            self.lmax_lo = np.zeros(reps)  # largest lower bet so far (where the per-point cap has bitten)
            # binned law of Z for the bets: bins of width zhat0 / 8 up to 10 zhat0, one overflow bin and
            # one bin for Z = 1 (clipped returns); per-bin sums keep the binned mean exact
            self.nbins = 82
            self.hist = np.zeros((reps, self.nbins))
            self.hsum = np.zeros((reps, self.nbins))
            self._add_hist(np.arange(reps), np.clip(z_pilot, 0, 1))

    def _add_hist(self, g, z):
        v = z / self.z0[g][:, None]
        b = np.minimum(np.floor(v * 8.0).astype(int), self.nbins - 2)
        b = np.where(z >= 1.0 - 1e-12, self.nbins - 1, b)
        flat = (np.arange(g.size)[:, None] * self.nbins + b).ravel()
        n = g.size * self.nbins
        self.hist[g] += np.bincount(flat, minlength=n).reshape(g.size, self.nbins)
        self.hsum[g] += np.bincount(flat, weights=z.ravel(), minlength=n).reshape(g.size, self.nbins)

    def _bets(self, g, zprev, s2prev, tt):
        s = self.s
        a_up, a_lo = 1.0 - zprev, zprev
        if self.bets == "tuned":
            dz = s.d * zprev if s.relative else self.dz[g][:, None]
            base = dz / s2prev
            return np.minimum(base, s.cs_cap / a_up), np.minimum(base, s.cs_cap / a_lo), a_up, a_lo
        if self.bets == "wsr":
            dt = tt * np.log1p(tt) if self.n_plan is None else self.n_plan
            base = np.sqrt(2.0 * self.logd / (dt * s2prev))
            return np.minimum(base, s.cs_cap / a_up), np.minimum(base, s.cs_cap / a_lo), a_up, a_lo
        if self.bets == "opt":
            # constant bets for the whole block, re-solved from the plug-in values at its start
            _, lu1, ll1 = opt_plan(zprev[:, 0], s2prev[:, 0], s.d, s.alpha, s.cs_cap, s.relative,
                                   d_z=None if s.relative else self.dz[g])
            lu = np.minimum(np.broadcast_to(lu1[:, None], zprev.shape), s.cs_cap / a_up)
            ll = np.minimum(np.broadcast_to(ll1[:, None], zprev.shape), s.cs_cap / a_lo)
            return lu, ll, a_up, a_lo
        raise ValueError(self.bets)

    def _gro_bets(self, g, zc):
        """Growth-rate optimal bets against ``zhat +- w`` under the binned empirical law of Z seen so far."""
        s = self.s
        w = s.d * zc if s.relative else self.dz[g]
        cnt = self.hist[g]
        p = cnt / cnt.sum(axis=1, keepdims=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            cen = np.where(cnt > 0, self.hsum[g] / cnt, 0.0)
        out = []
        for side, cap in ((+1, self.cap_up[g]), (-1, self.cap_lo[g])):
            lam = np.geomspace(1e-3, 1.0, 48)[None, :] * cap[:, None]
            m = zc + side * w
            xi = (m[:, None] - cen) * side  # up: m - Z; lo: Z - m
            val = np.log1p(np.maximum(lam[:, :, None] * xi[:, None, :], -0.999999))
            gro = (val * p[:, None, :]).sum(axis=2)
            out.append(lam[np.arange(g.size), gro.argmax(axis=1)])
        return out

    def update(self, z, rows, cmean=None):
        pos, g = self.sel(rows)
        k = z.shape[1]
        if g.size == 0:
            self.t += k
            return
        s = self.s
        z = z[pos]
        cmean = None if cmean is None else cmean[pos]
        tt = self.t + np.arange(1, k + 1, dtype=float)[None, :]
        cz = self.sum_z[g][:, None] + np.cumsum(z, axis=1)
        zhat = (self.w0 * self.z0[g][:, None] + cz) / (self.w0 + tt)
        zprev = np.concatenate([self.zhat[g][:, None], zhat[:, :-1]], axis=1)
        cd = self.sum_dev2[g][:, None] + np.cumsum((z - zhat) ** 2, axis=1)
        s2 = (self.w0 * self.v0[g][:, None] + cd) / (self.w0 + tt)
        s2prev = np.maximum(np.concatenate([self.s2[g][:, None], s2[:, :-1]], axis=1), 1e-8)
        zprev = np.clip(zprev, 1e-6, 1 - 1e-6)
        sm = None if cmean is None else self.sm[g][:, None] + np.cumsum(cmean, axis=1)
        sl_prev = self.sl[:, g].copy()
        if self.bets == "bet":
            U, L, lams = self._update_bet(g, z, zprev)
        else:
            e = (z - zprev) ** 2
            lu, ll, a_up, a_lo = self._bets(g, zprev, s2prev, tt)
            lams = (lu, ll)
            res = []
            for side, (lam, a) in enumerate(((lu, a_up), (ll, a_lo))):
                psi = lam * lam * psi_coef(lam * a)
                res.append((self.sl[side, g][:, None] + np.cumsum(lam, axis=1),
                            self.slz[side, g][:, None] + np.cumsum(lam * z, axis=1),
                            self.spe[side, g][:, None] + np.cumsum(psi * e, axis=1)))
            (slu, slzu, speu), (sll, slzl, spel) = res
            U = np.minimum((slzu + self.logd + speu) / slu, 1.0)
            L = np.maximum((slzl - self.logd - spel) / sll, 0.0)
            if s.intersect:  # valid: the bounds hold at all times simultaneously
                U = np.minimum.accumulate(np.concatenate([self.U[g][:, None], U], axis=1), axis=1)[:, 1:]
                L = np.maximum.accumulate(np.concatenate([self.L[g][:, None], L], axis=1), axis=1)[:, 1:]
            for side, (sl_, slz_, spe_) in enumerate(res):
                self.sl[side, g], self.slz[side, g], self.spe[side, g] = sl_[:, -1], slz_[:, -1], spe_[:, -1]
        if s.relative:
            cond = (L > 0) & ((U - L) <= s.d * (U + L))
        else:
            cond = (U - L) <= 2.0 * self.dz[g][:, None]
        anyh, idx = _first_hit(cond)
        h = np.nonzero(anyh)[0]
        j = idx[h]
        gh = g[h]
        self.T[gh] = self.t + j + 1
        self.L_at[gh] = L[h, j]
        self.U_at[gh] = U[h, j]
        if cmean is not None:
            slm = [self.slm[side, g][:, None] + np.cumsum(lams[side] * cmean, axis=1) for side in (0, 1)]
            cl = [sl_prev[side][:, None] + np.cumsum(lams[side], axis=1) for side in (0, 1)]
            self.cm_at[0, gh] = slm[0][h, j] / cl[0][h, j]
            self.cm_at[1, gh] = slm[1][h, j] / cl[1][h, j]
            self.cm_at[2, gh] = sm[h, j] / (self.t + j + 1)
            self.slm[0, g], self.slm[1, g] = slm[0][:, -1], slm[1][:, -1]
            self.sm[g] = sm[:, -1]
        self.active[gh] = False
        self.t += k
        self.sum_z[g] = cz[:, -1]
        self.sum_dev2[g] = cd[:, -1]
        self.zhat[g] = zhat[:, -1]
        self.s2[g] = s2[:, -1]
        self.U[g] = U[:, -1]
        self.L[g] = L[:, -1]

    @staticmethod
    def _quad_root(F0, D0, C0, thr, h):
        """Smallest y in [0, h] with ``F0 + D0 y - C0 y^2 / 2 >= thr`` (``F0 < thr``), or h if none."""
        disc = D0 * D0 - 2.0 * C0 * (thr - F0)
        with np.errstate(invalid="ignore", divide="ignore"):
            y = np.where(C0 > 0, (D0 - np.sqrt(np.maximum(disc, 0.0))) / C0, (thr - F0) / D0)
        ok = (disc >= 0) & (D0 > 0) & (y >= 0) & (y <= h)
        return np.where(ok, y, h)

    def _update_bet(self, g, z, zprev):
        """Betting sequence: per-step loop over the block, vectorised over rows and grid points.

        Bounds between grid points. Upper side: ``F(m) = sum log(1 + lam (m - Z_i))`` is
        nondecreasing and concave, and ``|F''|`` is nonincreasing in m, so on ``[m0, m1]``
        ``F(m) >= max(chord, F(m0) + F'(m0) y - |F''(m0)| y^2 / 2)`` with ``y = m - m0``.
        Lower side: ``F(m) = sum log(1 + lam (Z_i - m))`` is nonincreasing and concave, ``|F''|``
        is nondecreasing in m, so ``F(m) >= max(chord, F(m1) + |F'(m1)| y - |F''(m1)| y^2 / 2)``
        with ``y = m1 - m``. Every candidate mean where a lower bound reaches ``log(2/alpha)``
        is rejected, so the reported bounds are conservative. Above the core of the grid the lower bet is
        capped per point (``cs_cap / m``); there the log capital is still nonincreasing in m but not
        concave, so the bound is the largest rejected grid point, without interpolation.
        """
        k = z.shape[1]
        lu, ll = self._gro_bets(g, zprev[:, 0])  # bets fixed for the block (predictable)
        grid = self.grid[g]
        Fu, Fl = self.logK_up[g], self.logK_lo[g]
        Du, Dl = self.D_up[g], self.D_lo[g]
        Cu, Cl = self.C_up[g], self.C_lo[g]
        Uo, Lo = self.U[g].copy(), self.L[g].copy()
        U = np.empty((g.size, k))
        L = np.empty((g.size, k))
        thr = self.logd
        G = grid.shape[1]
        ar = np.arange(g.size)
        lu2 = lu[:, None]
        ll2 = np.minimum(ll[:, None], self.s.cs_cap / grid)  # per-point cap (binds only above the core)
        self.lmax_lo[g] = np.maximum(self.lmax_lo[g], ll)
        m_conc = self.s.cs_cap / self.lmax_lo[g]  # lower log capital is concave in m on [0, m_conc]
        for j in range(k):
            zj = z[:, j][:, None]
            xu = 1.0 + lu2 * (grid - zj)
            xl = 1.0 + ll2 * (zj - grid)
            Fu = Fu + np.log(xu)
            Fl = Fl + np.log(xl)
            qu = lu2 / xu
            ql = ll2 / xl
            Du = Du + qu
            Dl = Dl + ql
            Cu = Cu + qu * qu
            Cl = Cl + ql * ql
            # upper bound: first grid point that is rejected, then refine inside the cell below it
            ru = Fu >= thr
            anyu = ru.any(axis=1)
            iu = np.where(anyu, ru.argmax(axis=1), G)
            ub = np.ones(g.size)
            first = anyu & (iu == 0)
            ub[first] = grid[first, 0]
            mid = anyu & (iu > 0)
            if mid.any():
                i1 = iu[mid]
                r_ = ar[mid]
                F0, F1 = Fu[r_, i1 - 1], Fu[r_, i1]
                m0, m1 = grid[r_, i1 - 1], grid[r_, i1]
                h = m1 - m0
                y_ch = (thr - F0) / (F1 - F0) * h
                y_q = self._quad_root(F0, Du[r_, i1 - 1], Cu[r_, i1 - 1], thr, h)
                ub[mid] = m0 + np.minimum(y_ch, y_q)
            # lower bound: last grid point that is rejected, then refine inside the cell above it
            rl = Fl >= thr
            anyl = rl.any(axis=1)
            il = np.where(anyl, G - 1 - rl[:, ::-1].argmax(axis=1), -1)
            lb = np.zeros(g.size)
            last = anyl & (il == G - 1)
            lb[last] = grid[last, -1]
            mid = anyl & (il < G - 1)
            if mid.any():
                i0 = il[mid]
                r_ = ar[mid]
                F0, F1 = Fl[r_, i0], Fl[r_, i0 + 1]
                m0, m1 = grid[r_, i0], grid[r_, i0 + 1]
                h = m1 - m0
                y_ch = (thr - F1) / (F0 - F1) * h  # distance from m1 down to the chord crossing
                y_q = self._quad_root(F1, Dl[r_, i0 + 1], Cl[r_, i0 + 1], thr, h)
                # interpolate only where the bet was never capped (concave part); else the rejected grid point
                conc = m1 <= m_conc[r_]
                lb[mid] = np.where(conc, m1 - np.minimum(y_ch, y_q), m0)
            Uo = np.minimum(Uo, ub) if self.s.intersect else ub
            Lo = np.maximum(Lo, lb) if self.s.intersect else lb
            U[:, j] = Uo
            L[:, j] = Lo
        self.logK_up[g], self.logK_lo[g] = Fu, Fl
        self.D_up[g], self.D_lo[g] = Du, Dl
        self.C_up[g], self.C_lo[g] = Cu, Cl
        # update the binned law used by the bets (after the block: predictable for the next block)
        self._add_hist(g, np.clip(z, 0, 1))
        lams = (np.broadcast_to(lu[:, None], z.shape), np.broadcast_to(ll[:, None], z.shape))
        self.sl[0, g] += lu * k
        self.sl[1, g] += ll * k
        return U, L, lams

    def advance(self, k):
        self.t += k

    def finish(self):
        a = self.active.copy()
        g = np.nonzero(a)[0]
        if g.size:
            self.T[g] = self.t
            self.U_at[g] = self.U[g]
            self.L_at[g] = self.L[g]
            if self.slm.any():
                with np.errstate(invalid="ignore", divide="ignore"):
                    self.cm_at[0, g] = self.slm[0, g] / self.sl[0, g]
                    self.cm_at[1, g] = self.slm[1, g] / self.sl[1, g]
                    self.cm_at[2, g] = self.sm[g] / max(self.t, 1)
        return a


def pilot_scale(r_pilot):
    """Scaled MAD about zero of each row (``1.4826 * median |r|``).

    Where more than half of a row is exactly zero (price discreteness on real one-minute
    data) the median is zero; the scaled median of the nonzero ``|r|`` is used instead, as
    in the package's monitor (``VolatilityMonitor._fallback_scale``). Rows with no nonzero
    value get 0.
    """
    a = np.abs(np.asarray(r_pilot, float))
    sc = np.median(a, axis=1) * MAD_NORMAL_SCALE
    bad = sc <= 0
    if np.any(bad):
        ab = a[bad]
        has = (ab > 0).any(axis=1)
        fb = np.zeros(ab.shape[0])
        if np.any(has):
            nz = np.where(ab[has] > 0, ab[has], np.nan)
            fb[has] = np.nanmedian(nz, axis=1) * MAD_NORMAL_SCALE
        sc[bad] = fb
    return sc


def run(law: Law, setup: Setup, reps: int, seed: int, garch_targets=None) -> dict:
    """Run every procedure in ``setup.procs`` on the same ``reps`` paths of ``law``.

    In relative mode returns are divided by the pilot scale first (a predictable,
    scale-free normalisation; Chow and Robbins' ``1/n`` term is not scale-free).
    Returns ``{name: Result}`` with everything in original units, plus ``"_meta"``.
    For dependent laws (GARCH) ``garch_targets`` maps tau to (truncated variance,
    long-run variance of the truncated statistic) in original units.
    """
    s = setup
    rng = np.random.default_rng(seed)
    state = law.start(rng, reps)
    allrows = np.arange(reps)
    r_p, h_p, state = law.draw(rng, state, allrows, s.m)
    scale = pilot_scale(r_p)
    if np.any(scale <= 0):
        raise ValueError("zero pilot scale")
    unit = scale if s.relative else np.ones(reps)  # returns are divided by this
    rp = r_p / unit[:, None]
    tau = s.c * pilot_scale(rp)
    tau2 = tau**2
    if isinstance(law, Garch):
        th_raw = np.full(reps, 1.0) / unit**2
        nu_raw = law.nu_raw() * law.lrv_factor_raw() / unit**4  # long-run variance of r^2
        m2, lrv = garch_targets(tau * unit)
        th_tr = m2 / unit**2
        nu_tr = lrv / unit**4
    elif isinstance(law, Windows):
        th_raw = th_tr = np.full(reps, np.nan)
        nu_raw = nu_tr = np.full(reps, np.nan)
    else:
        th_raw = np.full(reps, law.theta_raw()) / unit**2
        nu_raw = law.nu_raw() / unit**4
        m2, m4 = law.trunc_moments(tau * unit)
        th_tr = m2 / unit**2
        nu_tr = (m4 - m2**2) / unit**4
    z = zq(s.alpha)
    if s.relative:
        ns_raw = z**2 * nu_raw / (s.d * th_raw) ** 2
        ns_tr = z**2 * nu_tr / (s.d * th_tr) ** 2
    else:
        ns_raw = z**2 * nu_raw / s.d**2 * np.ones(reps)
        ns_tr = z**2 * nu_tr / s.d**2
    hp = None if h_p is None else h_p / unit[:, None] ** 2
    xr_p = rp**2
    xt_p = np.minimum(xr_p, tau2[:, None])

    def cond_tr(h, rows):
        return h * kappa2(tau[rows][:, None] / np.sqrt(h))

    A = None if s.minrisk_A is None else np.full(reps, s.minrisk_A) * unit**4
    procs = {}
    for p in s.procs:
        kind = "raw" if p.endswith("-raw") else "trunc"
        xp = xr_p if kind == "raw" else xt_p
        cp = None if hp is None else (hp if kind == "raw" else cond_tr(hp, allrows))
        ns = ns_raw if kind == "raw" else ns_tr
        if p.startswith("CR-BM"):
            procs[p] = _MeanVarRule(kind, "bm", s, xp, ns, cp)
        elif p.startswith("CR"):
            procs[p] = _MeanVarRule(kind, "cr", s, xp, ns, cp)
        elif p.startswith("MinRisk"):
            procs[p] = _MeanVarRule(kind, "minrisk", s, xp, ns, cp, A=A)
        elif p.startswith("Oracle"):
            procs[p] = _MeanVarRule(kind, "oracle", s, xp, ns, cp)
        elif p.startswith("Stein"):
            procs[p] = None
        elif p.startswith("CS"):
            bets = p.split("-")[1]
            procs[p] = _CS(s, xt_p / tau2[:, None], bets=bets, init="wsr" if bets == "wsr" else "pilot",
                           dz=None if s.relative else s.d / tau2)
        else:
            raise ValueError(p)
    stein = [p for p in s.procs if p.startswith("Stein")]
    stein_N, stein_sum = {}, {}
    if stein:
        tq = stats.t.ppf(1 - s.alpha / 2, s.m - 1)
        for p in stein:
            xp = xr_p if p.endswith("-raw") else xt_p
            v = xp.var(axis=1, ddof=1)
            mu = xp.mean(axis=1)
            dd = s.d * mu if s.relative else s.d
            stein_N[p] = np.minimum(np.maximum(s.m, np.floor((tq / dd) ** 2 * v) + 1), s.n_max).astype(int)
            stein_sum[p] = xp.sum(axis=1)
    n = s.m
    while n < s.n_max:
        live = np.zeros(reps, bool)
        for o in procs.values():
            if o is not None:
                live |= o.active
        for p in stein:
            live |= stein_N[p] > n
        if not live.any():
            break
        rows = np.nonzero(live)[0]
        k = min(s.block, s.n_max - n)
        r, h, state = law.draw(rng, state, rows, k)
        r = r / unit[rows][:, None]
        hb = None if h is None else h / unit[rows][:, None] ** 2
        xr = r * r
        xt = np.minimum(xr, tau2[rows][:, None])
        ctr = None if hb is None else cond_tr(hb, rows)
        for p, o in procs.items():
            if o is None:
                continue
            if isinstance(o, _CS):
                o.update(xt / tau2[rows][:, None], rows, None if ctr is None else ctr / tau2[rows][:, None])
            else:
                o.update(xr if o.kind == "raw" else xt, rows, None if hb is None else (hb if o.kind == "raw" else ctr))
        for p in stein:
            x = xr if p.endswith("-raw") else xt
            cols = np.arange(n + 1, n + k + 1)[None, :]
            take = cols <= stein_N[p][rows][:, None]
            stein_sum[p][rows] += np.where(take, x, 0.0).sum(axis=1)
        n += k
    out = {}
    u2 = unit**2
    for p in s.procs:
        kind = "raw" if p.endswith("-raw") else "trunc"
        tgt = th_raw if kind == "raw" else th_tr
        ns = ns_raw if kind == "raw" else ns_tr
        if p.startswith("Stein"):
            N = stein_N[p]
            est = stein_sum[p] / N
            dd = s.d * est if s.relative else s.d
            out[p] = Result(N=N, est=est * u2, lo=(est - dd) * u2, hi=(est + dd) * u2, target=tgt * u2,
                            capped=N >= s.n_max, nstar=ns)
            continue
        o = procs[p]
        capped = o.finish()
        if isinstance(o, _CS):
            Lz, Uz = o.L_at, o.U_at
            if s.relative:
                with np.errstate(divide="ignore", invalid="ignore"):
                    est = np.where(Lz + Uz > 0, 2 * Lz * Uz / (Lz + Uz), 0.0)
            else:
                est = (Lz + Uz) / 2
            cm = o.cm_at * (tau2 * u2)[None, :]
            out[p] = Result(N=o.T + s.m, est=est * tau2 * u2, lo=Lz * tau2 * u2, hi=Uz * tau2 * u2,
                            target=tgt * u2, capped=capped, nstar=ns,
                            cond_target=None if h_p is None else cm[2],
                            cond_target_w=None if h_p is None else cm[:2])
            continue
        est = o.mean_at
        dd = s.d * est if s.relative else s.d
        loss = None
        if o.rule == "minrisk":
            loss = A * (est - tgt) ** 2 + o.N
        out[p] = Result(N=o.N, est=est * u2, lo=(est - dd) * u2, hi=(est + dd) * u2, target=tgt * u2,
                        capped=capped, nstar=ns, cond_target=None if h_p is None else o.cond_at * u2, loss=loss,
                        var_at=o.var_at * u2**2)
    out["_meta"] = dict(tau=tau * unit, unit=unit)
    return out


# ----------------------------------------------------------------------------
# Summaries
# ----------------------------------------------------------------------------


def covered(res: Result, d: float, relative: bool, target=None) -> np.ndarray:
    """Was the target within the prescribed accuracy of the estimate at N? (``|est - theta| <= d`` or ``<= d theta``)."""
    t = res.target if target is None else target
    tol = d * t if relative else d
    return np.abs(res.est - t) <= tol


def summarise(res: Result, d: float, relative: bool, alpha: float | None = None) -> dict:
    cov = covered(res, d, relative)
    N = res.N.astype(float)
    ns = np.asarray(res.nstar, float)
    fin = bool(np.isfinite(ns).all())
    tol = d * res.target if relative else d
    out = dict(
        coverage=float(cov.mean()),
        miss_under=float((res.target - res.est > tol).mean()),
        miss_over=float((res.est - res.target > tol).mean()),
        EN=float(N.mean()), sdN=float(N.std()),
        N_q05=float(np.quantile(N, 0.05)), N_q50=float(np.median(N)), N_q95=float(np.quantile(N, 0.95)),
        nstar=float(np.mean(ns)) if fin else math.inf,
        EN_over_nstar=float(np.mean(N / ns)) if fin else math.nan,
        regret=float(np.mean(N - ns)) if fin else math.nan,
        capped=float(np.mean(res.capped)),
        interval_cover=float(((res.lo <= res.target) & (res.target <= res.hi)).mean()),
        rel_bias=float(np.mean(res.est / res.target - 1.0)),
    )
    if res.var_at is not None and alpha is not None:
        # the studentised interval est +- z sqrt(var / N) (random width), for the decomposition of errors
        out["student_cover"] = float((np.abs(res.est - res.target) <= zq(alpha) * np.sqrt(res.var_at / N)).mean())
    if res.loss is not None:
        out["risk"] = float(np.mean(res.loss))
    return out


# ----------------------------------------------------------------------------
# Mapping the truncated variance back to sigma^2 under a fixed innovation shape
# ----------------------------------------------------------------------------


def sigma2_from_theta_tau(theta_tau, tau, shape: Law | None = None):
    """Solve ``E[min(sigma^2 eps^2, tau^2)] = theta_tau`` for ``sigma^2`` (vectorised, monotone in sigma).

    ``shape`` is a unit-variance law for ``eps`` (default Gaussian, closed form
    ``sigma^2 kappa2(tau / sigma)``). The map is exact if the shape is right.
    """
    theta_tau = np.atleast_1d(np.asarray(theta_tau, float))
    tau = np.broadcast_to(np.asarray(tau, float), theta_tau.shape)
    out = np.empty_like(theta_tau)
    for i, (th, ta) in enumerate(zip(theta_tau, tau)):
        if not (0 < th < ta**2):
            out[i] = math.inf if th >= ta**2 else math.nan
            continue

        def f(ls):
            sg = math.exp(0.5 * ls)
            if shape is None:
                v = sg * sg * float(kappa2(ta / sg))
            else:
                v = float(shape.trunc_moments(np.array([ta / sg]))[0][0]) * sg * sg
            return v - th
        lo, hi = math.log(th) - 1.0, math.log(th) + 20.0
        out[i] = math.exp(optimize.brentq(f, lo, hi, xtol=1e-12))
    return out


def elasticity_gaussian(a: float) -> float:
    """``d log theta_tau / d log sigma^2`` for Gaussian returns at cut-off ``a = tau / sigma``: ``1 - a k2'(a) / (2 k2(a))``."""
    dk = 2.0 * a * 2.0 * stats.norm.sf(a)  # d/da E[min(Z^2, a^2)] = 2 a P(|Z| > a)
    return 1.0 - a * dk / (2.0 * float(kappa2(a)))


def garch_target_table(law: Garch, taus: np.ndarray, n: int = 400_000, reps: int = 8, seed: int = 7,
                       lags: int = 400):
    """Stationary truncated variance ``E min(r^2, tau^2)`` and long-run variance of ``min(r^2, tau^2)`` on a tau grid.

    Monte Carlo from ``reps`` long stationary paths; the long-run variance uses a
    Bartlett window of ``lags`` lags (``lags`` well beyond ``1 / (1 - a - b)``).
    Returns a callable ``tau -> (theta_tau, lrv_tau)`` by linear interpolation.
    """
    rng = np.random.default_rng(seed)
    st = law.start(rng, reps)
    r, _, _ = law.draw(rng, st, np.arange(reps), n)
    m2 = np.empty(taus.size)
    lr = np.empty(taus.size)
    for i, t in enumerate(taus):
        x = np.minimum(r * r, t * t)
        m2[i] = x.mean()
        xc = x - x.mean(axis=1, keepdims=True)
        f = np.fft.rfft(xc, 2 * n, axis=1)
        ac = np.fft.irfft(f * np.conj(f), axis=1)[:, : lags + 1] / n
        ac = ac.mean(axis=0)
        w = 1.0 - np.arange(1, lags + 1) / (lags + 1)
        lr[i] = ac[0] + 2.0 * np.sum(w * ac[1:])

    def f(tau):
        tau = np.asarray(tau, float)
        return np.interp(tau, taus, m2), np.interp(tau, taus, lr)

    return f
