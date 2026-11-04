"""Growth rates, detection delay approximations and bet design for the betting e-detector.

The theory behind every function is in ``docs/theory_delay.md``. Notation:

* ``Y = min(r^2, tau^2) / m0`` is the normalised truncated squared return, with values in ``[0, B]``
  where ``B = tau^2 / m0`` (``B = inf`` without truncation).
* An increase bet is ``f = 1 + lam (Y - 1)`` with ``lam`` in ``[0, 1)``; a decrease bet is
  ``f = 1 + lam (1 - Y)`` with ``lam`` in ``[0, 1 / (B - 1))``.
* The growth rate of a bet under a post-change law ``Q`` is ``g(lam) = E_Q[log f]``.

The module works with *population* quantities: ``baseline`` returns the truncation level and
null mean implied by a return model, and ``y_law`` returns a quadrature representation of the
law of ``Y`` (nodes and weights, including the atom at ``B``), so that every expectation
``E[h(Y)]`` is a weighted sum. Nothing in ``detectors.py`` depends on this module.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, special, stats

from .robust import MAD_NORMAL_SCALE

_GL_X, _GL_W = np.polynomial.legendre.leggauss(64)
EULER_GAMMA = 0.5772156649015329


# ----------------------------------------------------------------------------
# Return models (the ones used in experiment E4)
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class ReturnModel:
    """Symmetric per-step return model with pre-change standard deviation 1.

    ``kind``:
      * ``"gaussian"``: ``r = rho Z``.
      * ``"student_t"``: ``r = rho T sqrt((df - 2) / df)`` with ``T ~ t_df`` (unit variance when rho = 1).
      * ``"jumps"``: ``r = rho Z + I J`` with ``I ~ Bernoulli(p_jump)`` and ``J ~ N(0, jump_sd^2)``.
        The jump part is not scaled by ``rho`` (as in ``simulate.change_returns``).

    ``rho`` is the post-change ratio of standard deviations (1 means no change).
    """

    kind: str = "gaussian"
    df: float = 5.0
    p_jump: float = 0.005
    jump_sd: float = 8.0

    def components(self, rho: float = 1.0):
        """List of ``(weight, df or None, scale)``: |r| is a mixture of scaled |normal| or |t| laws."""
        if self.kind == "gaussian":
            return [(1.0, None, float(rho))]
        if self.kind == "student_t":
            return [(1.0, self.df, float(rho) * np.sqrt((self.df - 2.0) / self.df))]
        if self.kind == "jumps":
            return [(1.0 - self.p_jump, None, float(rho)), (self.p_jump, None, float(np.hypot(rho, self.jump_sd)))]
        raise ValueError(f"unknown model kind {self.kind!r}")

    @staticmethod
    def _dist(df):
        return stats.norm if df is None else stats.t(df)

    def pdf_abs(self, a, rho: float = 1.0):
        """Density of ``|r|`` at ``a >= 0``."""
        a = np.asarray(a, float)
        out = np.zeros_like(a)
        for w, df, s in self.components(rho):
            out = out + w * 2.0 * self._dist(df).pdf(a / s) / s
        return out

    def logpdf_abs(self, a, rho: float = 1.0):
        """Log density of ``|r|`` at ``a >= 0`` (stable in the far tail)."""
        a = np.asarray(a, float)
        terms = [np.log(w) + np.log(2.0) + self._dist(df).logpdf(a / s) - np.log(s)
                 for w, df, s in self.components(rho)]
        return special.logsumexp(np.stack(terms), axis=0)

    def sf_abs(self, a, rho: float = 1.0):
        """``P(|r| > a)``."""
        a = np.asarray(a, float)
        out = np.zeros_like(a)
        for w, df, s in self.components(rho):
            out = out + w * 2.0 * self._dist(df).sf(a / s)
        return out

    def sample(self, size, rho: float = 1.0, rng: np.random.Generator | None = None):
        """Draw returns (same construction as ``simulate.change_returns``)."""
        rng = np.random.default_rng() if rng is None else rng
        if self.kind == "gaussian":
            return rho * rng.standard_normal(size)
        if self.kind == "student_t":
            return rho * rng.standard_t(self.df, size) * np.sqrt((self.df - 2.0) / self.df)
        if self.kind == "jumps":
            r = rho * rng.standard_normal(size)
            return r + (rng.random(size) < self.p_jump) * rng.normal(0.0, self.jump_sd, size)
        raise ValueError(f"unknown model kind {self.kind!r}")


def _abs_nodes(model: ReturnModel, rho: float, upper: float):
    """Gauss-Legendre nodes and weights for integrals over ``|r|`` in ``(0, upper)``."""
    scales = [s for _, _, s in model.components(rho)] + [s for _, _, s in model.components(1.0)]
    s_min, s_max = min(scales), max(scales)
    heavy = model.kind == "student_t"
    top = upper if np.isfinite(upper) else (3000.0 if heavy else 60.0) * s_max
    edges = [0.0]
    e = s_min / 16.0
    while e < top:
        edges.append(e)
        e *= 1.6
    edges.append(top)
    edges = np.unique(np.clip(np.asarray(edges), 0.0, top))
    lo, hi = edges[:-1], edges[1:]
    x = (0.5 * (hi - lo))[:, None] * _GL_X[None, :] + (0.5 * (hi + lo))[:, None]
    w = (0.5 * (hi - lo))[:, None] * _GL_W[None, :]
    return x.ravel(), w.ravel()


def baseline(model: ReturnModel, c: float = 3.0):
    """Population truncation level and null mean: ``tau = c * scaled MAD`` of pre-change returns,
    ``m0 = E[min(r^2, tau^2)]``. ``c = inf`` means no truncation (``m0 = E[r^2]``)."""
    med = optimize.brentq(lambda a: float(model.sf_abs(a, 1.0)) - 0.5, 1e-8, 100.0, xtol=1e-14)
    tau = c * med * MAD_NORMAL_SCALE if np.isfinite(c) else np.inf
    a, w = _abs_nodes(model, 1.0, tau)
    m0 = float(np.sum(w * model.pdf_abs(a, 1.0) * np.minimum(a**2, tau**2)))
    if np.isfinite(tau):
        m0 += float(tau**2 * model.sf_abs(tau, 1.0))
    return float(tau), m0


@dataclass
class YLaw:
    """Discrete quadrature representation of the law of ``Y = min(r^2, tau^2) / m0``."""

    y: np.ndarray
    w: np.ndarray
    B: float  # upper end of the support (tau^2 / m0, or inf)

    def expect(self, h):
        return float(np.sum(self.w * h(self.y)))

    @property
    def mean(self):
        return float(np.sum(self.w * self.y))

    @property
    def var(self):
        mu = self.mean
        return float(np.sum(self.w * (self.y - mu) ** 2))


def y_law(model: ReturnModel, rho: float, tau: float, m0: float) -> YLaw:
    """Law of ``Y`` when returns follow ``model`` with sd ratio ``rho``."""
    a, w = _abs_nodes(model, rho, tau)
    y = np.minimum(a**2, tau**2) / m0
    wt = w * model.pdf_abs(a, rho)
    B = tau**2 / m0 if np.isfinite(tau) else np.inf
    if np.isfinite(tau):
        y = np.append(y, B)
        wt = np.append(wt, float(model.sf_abs(tau, rho)))
    return YLaw(y=y, w=wt, B=float(B))


def sample_y(model: ReturnModel, rho: float, tau: float, m0: float, size, rng: np.random.Generator):
    """Draws of ``Y`` (exact sampling, used for Monte Carlo checks and delay simulations)."""
    r = model.sample(size, rho, rng)
    return np.minimum(r**2, tau**2) / m0


# ----------------------------------------------------------------------------
# Growth rates and information numbers
# ----------------------------------------------------------------------------


def lam_max(direction: str, B: float) -> float:
    """Largest admissible bet (exclusive): 1 for increases, ``1 / (B - 1)`` for decreases."""
    if direction == "up":
        return 1.0
    if direction == "down":
        if not np.isfinite(B):
            return 0.0  # no valid decrease bet without an upper bound (docs/math.md 5.3)
        return 1.0 / (B - 1.0)
    raise ValueError("direction must be 'up' or 'down'")


def _gain(y, direction):
    return (y - 1.0) if direction == "up" else (1.0 - y)


def growth_rate(lams, law: YLaw, direction: str = "up"):
    """``g(lam) = E[log(1 + lam * gain)]`` for each bet in ``lams`` (gain = Y - 1 or 1 - Y)."""
    lams = np.atleast_1d(np.asarray(lams, float))
    gain = _gain(law.y, direction)
    with np.errstate(divide="ignore", invalid="ignore"):
        lf = np.log1p(lams[:, None] * gain[None, :])
    lf = np.where(law.w[None, :] > 0, lf, 0.0)
    return lf @ law.w


def log_factor_moments(lam: float, law: YLaw, direction: str = "up"):
    """Mean and variance of ``log f`` for one bet."""
    gain = _gain(law.y, direction)
    lf = np.log1p(lam * gain)
    mu = float(np.sum(law.w * lf))
    return mu, float(np.sum(law.w * (lf - mu) ** 2))


def best_bet(law: YLaw, direction: str = "up"):
    """Growth-optimal constant bet ``lam*`` and its growth ``KL_inf = max_lam g(lam)``.

    By the duality in docs/theory_delay.md (Lemma 2) the maximum equals the smallest KL divergence
    from the post-change law of Y to the nonparametric null (all laws on [0, B] with mean <= 1 for
    an increase, >= 1 for a decrease).
    """
    top = lam_max(direction, law.B)
    if top <= 0:
        return 0.0, 0.0
    gain = _gain(law.y, direction)
    if np.sum(law.w * gain) <= 0:
        return 0.0, 0.0
    hi = top * (1 - 1e-12)
    res = optimize.minimize_scalar(lambda l: -float(growth_rate([l], law, direction)[0]),
                                   bounds=(0.0, hi), method="bounded", options={"xatol": 1e-10})
    lam = float(res.x)
    return lam, float(growth_rate([lam], law, direction)[0])


def approx_best_bet(law: YLaw, direction: str = "up"):
    """Second-order (aGRAPA type) approximation ``lam ~ d / (v + d^2)`` with ``d = |E Y - 1|``,
    clipped to the admissible range."""
    d = (law.mean - 1.0) if direction == "up" else (1.0 - law.mean)
    lam = max(d, 0.0) / (law.var + d**2)
    return float(min(lam, 0.999 * lam_max(direction, law.B)))


def kl_gaussian_scale(rho):
    """``KL(N(0, rho^2) || N(0, 1)) = (rho^2 - 1) / 2 - log(rho)`` per observation."""
    rho = np.asarray(rho, float)
    return 0.5 * (rho**2 - 1.0) - np.log(rho)


def kl_truncated(model: ReturnModel, rho: float, tau: float):
    """KL divergence between post- and pre-change laws of ``min(|r|, tau)``.

    This is the information per observation in the truncated statistic ``X = min(r^2, tau^2)``
    (``tau = inf`` gives the information in the raw return). It bounds the growth of every
    e-variable that is valid for the true pre-change law of X (Lemma 1 in docs/theory_delay.md).
    """
    a, w = _abs_nodes(model, rho, tau)
    lq = model.logpdf_abs(a, rho)
    lp = model.logpdf_abs(a, 1.0)
    integrand = np.exp(lq) * (lq - lp)
    kl = float(np.sum(w * integrand))
    if np.isfinite(tau):
        qs, ps = float(model.sf_abs(tau, rho)), float(model.sf_abs(tau, 1.0))
        if qs > 0:
            kl += qs * np.log(qs / ps)
    return kl


# ----------------------------------------------------------------------------
# Delay: bounds and approximations
# ----------------------------------------------------------------------------


def delay_first_order(A: float, g: float) -> float:
    """First-order delay ``log(A) / g`` (``inf`` if ``g <= 0``)."""
    return float(np.log(A) / g) if g > 0 else np.inf


def overshoot_bound(lam: float, law: YLaw, direction: str = "up") -> float:
    """Upper bound on the expected overshoot of the log wealth over any level:
    the smaller of the largest one-step log gain and Lorden's (1970) bound ``E[(Z+)^2] / E[Z]``."""
    gain = _gain(law.y, direction)
    z = np.log1p(lam * gain)
    g = float(np.sum(law.w * z))
    sup = float(np.max(z[law.w > 0]))
    if g <= 0:
        return np.inf
    lorden = float(np.sum(law.w * np.maximum(z, 0.0) ** 2) / g)
    return min(sup, lorden)


def wald_delay_bound(A: float, grid, weights, law: YLaw, direction: str = "up"):
    """Rigorous bound on the worst-case conditional expected delay (Theorem 1, docs/theory_delay.md):

    ``sup_nu ess sup E[(T_A - nu + 1)^+ | F_{nu-1}] <= min_k (log(A / w_k) + overshoot_k) / g_k``

    for the SR or CUSUM type mixture with fixed bets ``grid`` and weights ``weights``, when the
    post-change Y are i.i.d. from ``law``. Returns the bound and the index k that attains it.
    """
    grid = np.asarray(grid, float)
    weights = np.asarray(weights, float)
    g = growth_rate(grid, law, direction)
    best, arg = np.inf, -1
    for k in range(len(grid)):
        if g[k] <= 0:
            continue
        val = (np.log(A / weights[k]) + overshoot_bound(grid[k], law, direction)) / g[k]
        if val < best:
            best, arg = float(val), k
    return best, arg


def zeta_brownian(g: float, s2: float) -> float:
    """``E log int_0^inf exp(-(g t + s W_t)) dt = log(2 / s2) - digamma(2 g / s2)`` (Dufresne 1990).

    Brownian approximation of the Shiryaev-Roberts "memory" term ``E log sum_{j>=0} exp(-S_j)``.
    """
    return float(np.log(2.0 / s2) - special.digamma(2.0 * g / s2))


def renewal_constants(lam: float, draw_y, direction: str = "up", n_walks: int = 4000,
                      max_steps: int = 20000, rng: np.random.Generator | None = None):
    """Monte Carlo estimates of the two second-order constants of the SR delay for one bet.

    ``kappa``: limiting mean overshoot of the log-wealth random walk ``S_n = sum log f`` over a high
    level, from the ladder-height formula ``E[H^2] / (2 E[H])``.
    ``zeta``: ``E log V`` with ``V = sum_{j >= 0} exp(-S_j)`` (the SR statistic's memory of all
    start times). ``draw_y(size, rng)`` returns i.i.d. post-change draws of Y.
    Returns a dict with both constants and their Monte Carlo standard errors.
    """
    sign = 1.0 if direction == "up" else -1.0
    return renewal_constants_increments(lambda size, r: np.log1p(lam * sign * (draw_y(size, r) - 1.0)),
                                        n_walks, max_steps, rng)


def renewal_constants_increments(draw_z, n_walks: int = 4000, max_steps: int = 20000,
                                 rng: np.random.Generator | None = None):
    """``renewal_constants`` for a random walk with i.i.d. increments ``draw_z(size, rng)`` (positive mean),
    for example the log likelihood ratio of a parametric detector."""
    rng = np.random.default_rng() if rng is None else rng
    # ladder heights
    H = np.full(n_walks, np.nan)
    S = np.zeros(n_walks)
    alive = np.ones(n_walks, bool)
    for _ in range(max_steps):
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        S[idx] += draw_z(idx.size, rng)
        done = S[idx] > 0
        H[idx[done]] = S[idx[done]]
        alive[idx[done]] = False
    H = H[np.isfinite(H)]
    kappa = float(np.mean(H**2) / (2 * np.mean(H)))
    # crude delta-method standard error
    m1, m2 = np.mean(H), np.mean(H**2)
    cov = np.cov(np.vstack([H, H**2]))
    grad = np.array([-m2 / (2 * m1**2), 1 / (2 * m1)])
    kappa_se = float(np.sqrt(grad @ cov @ grad / H.size))
    # memory term
    S = np.zeros(n_walks)
    V = np.ones(n_walks)
    alive = np.ones(n_walks, bool)
    for _ in range(max_steps):
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        S[idx] += draw_z(idx.size, rng)
        V[idx] += np.exp(-S[idx])
        alive[idx[S[idx] > 40.0]] = False
    logV = np.log(V)
    return dict(kappa=kappa, kappa_se=kappa_se, zeta=float(logV.mean()),
                zeta_se=float(logV.std(ddof=1) / np.sqrt(n_walks)), unfinished=int(alive.sum()))


def delay_second_order(A: float, weight: float, g: float, kappa: float, zeta: float) -> float:
    """Second-order SR delay for a single dominant bet (change at the first monitored step):
    ``(log(A / w) + kappa - zeta) / g``. Heuristic for mixtures, see docs/theory_delay.md."""
    if g <= 0:
        return np.inf
    return float(max(np.log(A / weight) + kappa - zeta, 0.0) / g)


def speed_limit(A: float, factors, weights) -> int:
    """Smallest n at which the SR mixture can reach ``A`` when every bet earns its largest possible
    factor at every step (``factors[k] = sup f``). This mixture, with these bets and weights, cannot
    alarm sooner. For the floor of every valid decrease detector of SR form see
    ``decrease_speed_limit``."""
    factors = np.asarray(factors, float)
    weights = np.asarray(weights, float)
    M = np.zeros_like(factors)
    for n in range(1, 10**7):
        M = factors * (M + 1.0)
        if float(M @ weights) >= A:
            return n
    return -1


def sr_steps_constant_factor(A: float, F: float) -> int:
    """Steps a Shiryaev-Roberts statistic started at 0 needs to reach ``A`` when every step multiplies
    by the same factor ``F > 1``.

    ``M_n = sum_{j=1}^{n} F^(n-j+1) = F (F^n - 1) / (F - 1)``, so ``M_n >= A`` exactly when
    ``F^n >= 1 + A (F - 1) / F`` (docs/theory_delay.md, Corollary 4).
    """
    if not F > 1:
        raise ValueError("F must exceed 1")
    n = int(np.ceil(np.log1p(A * (F - 1.0) / F) / np.log(F)))
    while F * (F**n - 1.0) / (F - 1.0) < A:  # guard against rounding at an exact integer, both ways
        n += 1
    while n > 1 and F * (F ** (n - 1) - 1.0) / (F - 1.0) >= A:
        n -= 1
    return max(n, 1)


def decrease_speed_limit(A: float, B: float) -> int:
    """Fewest steps in which any Shiryaev-Roberts e-detector whose factors are valid for the whole
    nonparametric null can reach ``A`` against a decrease, ``B = tau^2 / m0``.

    Every such factor is at most ``B / (B - 1)`` (Lemma 2a and Corollary 4 of docs/theory_delay.md),
    so ``M_n <= B ((B / (B - 1))^n - 1)`` for every mixture of bets, fixed or chosen predictably.
    The floor is ``ceil(log(1 + A / B) / log(B / (B - 1)))``, reached in the limit by a single bet
    at ``lam_max`` when every return is zero. A grid that puts weight on smaller bets is slower
    (``speed_limit``).
    """
    if not (np.isfinite(B) and B > 1):
        raise ValueError("B must be finite and above 1")
    return sr_steps_constant_factor(A, B / (B - 1.0))


# ----------------------------------------------------------------------------
# Fast simulation of the SR / CUSUM mixture on precomputed Y (for analysis only)
# ----------------------------------------------------------------------------


def run_mixture(y: np.ndarray, grid, weights, direction: str = "up", threshold: float = 1000.0,
                kind: str = "SR", init=None, record_horizons=()):
    """Run the mixture e-detector on rows of ``y`` (already normalised by m0).

    Same recursion as ``detectors.betting_edetector`` (fixed grid, no adaptive bet), but with
    arbitrary weights, an optional initial state ``init`` (shape (reps, K)), and extra outputs:

    * ``alarm``: first index with mixture >= threshold, -1 if none.
    * ``at_alarm``: per-component statistic at the alarm (nan if none), shape (reps, K).
    * ``state``: final per-component statistic, shape (reps, K).
    * ``at_horizon``: dict m -> per-component statistic after m steps (frozen at the alarm).
    """
    y = np.atleast_2d(np.asarray(y, float))
    reps, n = y.shape
    grid = np.asarray(grid, float)
    weights = np.asarray(weights, float)
    sign = 1.0 if direction == "up" else -1.0
    S = np.zeros((reps, len(grid))) if init is None else np.array(init, float)
    alarm = np.full(reps, -1)
    at_alarm = np.full((reps, len(grid)), np.nan)
    active = np.ones(reps, bool)
    horizons = {int(m): None for m in record_horizons}
    for t in range(n):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        f = 1.0 + grid[None, :] * (sign * (y[idx, t] - 1.0))[:, None]
        with np.errstate(over="ignore"):
            S[idx] = f * (S[idx] + 1.0) if kind == "SR" else f * np.maximum(S[idx], 1.0)
        M = S[idx] @ weights
        hit = M >= threshold
        alarm[idx[hit]] = t
        at_alarm[idx[hit]] = S[idx[hit]]
        active[idx[hit]] = False
        if (t + 1) in horizons:
            horizons[t + 1] = S.copy()
    for m in horizons:
        if horizons[m] is None:
            horizons[m] = S.copy()
    return dict(alarm=alarm, at_alarm=at_alarm, state=S, at_horizon=horizons)


# ----------------------------------------------------------------------------
# Exact worst case over the nonparametric null (single bet), by dynamic programming
# ----------------------------------------------------------------------------


def worst_case_false_alarm(lam: float, B: float, A: float, horizon: int, direction: str = "up",
                           n_states: int = 700, n_left: int = 30, n_right: int = 90, compute_upper: bool = True):
    """Bracket ``[lower, upper]`` for ``sup P(T_A <= m)`` over every null law, for m = 1..horizon.

    The supremum is over all adaptive choices of the conditional law of ``Y_t`` on ``[0, B]`` with
    conditional mean <= 1 (increase) or >= 1 (decrease), for the single-bet SR statistic
    ``M_t = f_t (M_{t-1} + 1)``, ``M_0 = 0``, alarm at ``M_t >= A``. With ``u = Y - 1`` (increase)
    or ``u = 1 - Y`` (decrease) the bet is ``1 + lam u``, ``u`` lies in ``[u_lo, u_hi]`` and has
    conditional mean <= 0. The value function is nondecreasing in M, so a mean-zero two-point law
    is optimal at each step (concave envelope). ``upper`` rounds every transition up on a log grid
    of states and evaluates cells at their least favourable corner, so it is a rigorous upper bound
    (up to floating point). ``lower`` is the value of an actual adversary restricted to the grid
    and rounding down, so it is attained (a lower bound on the supremum). Rounding up accumulates
    over many steps, so ``upper`` is only informative for short horizons or fine grids; Theorem B
    (``horizon / A``) is always a valid upper bound. ``compute_upper=False`` skips it (returns nan).
    """
    if direction == "up":
        u_lo, u_hi = -1.0, B - 1.0
    else:
        u_lo, u_hi = 1.0 - B, 1.0
    fmin = 1.0 + lam * u_lo
    if fmin <= 0:
        raise ValueError("bet too large for this B")
    states = np.concatenate([[0.0], np.geomspace(fmin, A, n_states)[:-1]])  # last point < A
    a = np.linspace(u_lo, 0.0, n_left + 1)  # left cells [a_i, a_{i+1}]
    b = np.linspace(0.0, u_hi, n_right + 1)  # right cells [b_j, b_{j+1}]
    S = states[:, None]

    def nxt(u):
        return (1.0 + lam * u)[None, :] * (S + 1.0)

    h_a, h_b = nxt(a), nxt(b)
    # upper: round up to the next state (>= A means absorbed)
    up_a = np.searchsorted(states, h_a, side="left")
    up_b = np.searchsorted(states, h_b, side="left")
    # lower: round down
    lo_a = np.searchsorted(states, h_a, side="right") - 1
    lo_b = np.searchsorted(states, h_b, side="right") - 1
    absorbed_a, absorbed_b = h_a >= A, h_b >= A
    ns = len(states)
    # weights on the right point: p2 = -u1 / (u2 - u1), u1 < 0 <= u2
    a_left = a[:-1]
    p2_up = -a_left[:, None] / (b[None, :-1] - a_left[:, None])  # left ends of both cells
    a_in = a[:-1]
    b_in = b[1:]
    p2_lo = -a_in[:, None] / (b_in[None, :] - a_in[:, None])
    V_up = np.zeros(ns + 1)  # index ns = absorbed (value 1)
    V_lo = np.zeros(ns + 1)
    V_up[ns] = V_lo[ns] = 1.0
    upper, lower = np.empty(horizon), np.empty(horizon)
    ia_up = np.where(absorbed_a, ns, np.minimum(up_a, ns))
    ib_up = np.where(absorbed_b, ns, np.minimum(up_b, ns))
    ia_lo = np.where(absorbed_a, ns, lo_a)
    ib_lo = np.where(absorbed_b, ns, lo_b)
    for k in range(horizon):
        if compute_upper:
            pa = V_up[ia_up]  # (ns, n_left + 1)
            pb = V_up[ib_up]  # (ns, n_right + 1)
            r1 = pa[:, 1:]  # right end of each left cell
            r2 = pb[:, 1:]  # right end of each right cell
            val = r1[:, :, None] + p2_up[None, :, :] * (r2[:, None, :] - r1[:, :, None])
            new_up = val.max(axis=(1, 2))
        else:
            new_up = np.full(ns, np.nan)
        qa = V_lo[ia_lo]
        qb = V_lo[ib_lo]
        l1 = qa[:, :-1]
        l2 = qb[:, 1:]
        val_lo = l1[:, :, None] + p2_lo[None, :, :] * (l2[:, None, :] - l1[:, :, None])
        new_lo = np.maximum(val_lo.max(axis=(1, 2)), qa[:, -1])  # qa[:, -1] is u = 0 (degenerate law)
        V_up[:ns] = np.minimum(new_up, 1.0)
        V_lo[:ns] = new_lo
        upper[k] = V_up[0]
        lower[k] = V_lo[0]
    return lower, upper


# ----------------------------------------------------------------------------
# Bet design
# ----------------------------------------------------------------------------


def sr_delay_approx(A: float, weight: float, g: float, s2: float) -> float:
    """Design approximation of the SR delay of one bet with weight ``w`` inside a mixture:
    ``log(1 + (A / w) exp(-zeta)) / g`` with the Brownian memory term ``zeta`` capped at
    ``log(1 / g) + gamma`` (its value for a growth-optimal bet).

    For large ``A g`` this is the second-order formula ``(log(A / w) - zeta) / g`` without overshoot;
    for a bet with almost no growth it tends to ``A e^{-gamma} / w`` instead of becoming negative.
    The cap keeps poorly matched bets (log-wealth variance much larger than twice the growth),
    where the Brownian memory formula is not reliable, on the pessimistic side.
    """
    if g <= 0:
        return np.inf
    zeta = min(zeta_brownian(g, s2), np.log(1.0 / g) + EULER_GAMMA)
    return float(np.log1p((A / weight) * np.exp(-zeta)) / g)


def oracle_delays(A: float, laws, direction: str = "up", use_memory: bool = True):
    """Approximate delay of the single best bet for each law (``sr_delay_approx`` with weight 1,
    or ``log(A) / g*`` when ``use_memory`` is False)."""
    out = []
    for law in laws:
        lam_o, g_o = best_bet(law, direction)
        _, s2_o = log_factor_moments(lam_o, law, direction)
        out.append(sr_delay_approx(A, 1.0, g_o, s2_o) if use_memory else np.log(A) / g_o)
    return np.asarray(out)


def design_ratio(A: float, grid, weights, laws, direction: str = "up", use_memory: bool = True, oracle=None):
    """Approximate ratio of mixture delay to oracle (single best bet) delay, worst case over ``laws``.

    The mixture delay is approximated by ``min_k sr_delay_approx(A, w_k, g_k, s2_k)`` (each bet
    alone, paying ``log(1 / w_k)``), the oracle delay by ``sr_delay_approx(A, 1, g*, s2*)``. With
    ``use_memory=False`` both use the first-order ``log(A / w) / g``. The approximation ignores the
    overshoot (present in both) and the extra speed a mixture gains from its other bets, so it is
    pessimistic for mixtures. Returns (worst ratio, ratios per law).
    """
    grid = np.asarray(grid, float)
    weights = np.asarray(weights, float)
    if oracle is None:
        oracle = oracle_delays(A, laws, direction, use_memory)
    ratios = []
    for law, d_o in zip(laws, oracle):
        best = np.inf
        g = growth_rate(grid, law, direction)
        for k in range(len(grid)):
            if g[k] <= 0 or weights[k] <= 0:
                continue
            if use_memory:
                _, s2_k = log_factor_moments(grid[k], law, direction)
                d_k = sr_delay_approx(A, weights[k], g[k], s2_k)
            else:
                d_k = np.log(A / weights[k]) / g[k]
            best = min(best, d_k)
        ratios.append(best / d_o)
    ratios = np.asarray(ratios)
    return float(ratios.max()), ratios


def optimise_grid(A: float, laws, K: int, direction: str = "up", lam_cap: float | None = None,
                  free_weights: bool = True, use_memory: bool = True, seed: int = 0):
    """Choose ``K`` bets (and weights) minimising ``design_ratio`` over the given post-change laws.

    Step 1: search geometric grids between ``lam_lo`` and ``lam_hi`` (uniform weights) on a 2-D
    lattice. Step 2: Nelder-Mead polish of free positions (and softmax weights if ``free_weights``)
    from the best lattice point and a few random perturbations. Returns (grid, weights, worst ratio).
    """
    B = laws[0].B
    cap = lam_cap if lam_cap is not None else 0.999 * lam_max(direction, B)
    oracle_d = oracle_delays(A, laws, direction, use_memory)

    def ratio(lam, w):
        return design_ratio(A, lam, w, laws, direction, use_memory, oracle=oracle_d)[0]

    lattice = np.geomspace(cap * 1e-2, cap, 36)
    best_val, best_lam = np.inf, None
    for i, lo in enumerate(lattice):
        for hi in lattice[i:]:
            if K == 1 and hi != lo:
                continue
            lam = np.geomspace(lo, hi, K) if K > 1 else np.array([lo])
            if K > 1 and hi == lo:
                continue
            v = ratio(lam, np.full(K, 1.0 / K))
            if v < best_val:
                best_val, best_lam = v, lam
    best_w = np.full(K, 1.0 / K)

    def unpack(z):
        lam = np.sort(cap / (1.0 + np.exp(-z[:K])))
        if free_weights and K > 1:
            e = np.exp(z[K:] - z[K:].max())
            w = e / e.sum()
        else:
            w = np.full(K, 1.0 / K)
        return lam, w

    def obj(z):
        lam, w = unpack(z)
        return ratio(lam, w)

    rng = np.random.default_rng(seed)
    z_best = np.log(best_lam / (cap - best_lam))
    if free_weights and K > 1:
        z_best = np.concatenate([z_best, np.zeros(K)])
    for trial in range(4):
        z0 = z_best if trial == 0 else z_best + 0.2 * rng.standard_normal(z_best.size)
        res = optimize.minimize(obj, z0, method="Nelder-Mead",
                                options={"maxiter": 3000, "xatol": 1e-5, "fatol": 1e-7})
        if res.fun < best_val - 1e-12:
            best_val = float(res.fun)
            best_lam, best_w = unpack(res.x)
            z_best = res.x
    return np.asarray(best_lam), np.asarray(best_w), float(best_val)
