"""Limits that hold for every valid procedure, the plug-in baseline in the joint regime, and the decrease
bet near its boundary (TH5).

Theory and proofs: ``docs/theory_limits.md``. Experiment: ``experiments/th5_limits.py``. Tests:
``tests/test_limits.py``.

Notation as in ``docs/math.md`` (sections 5 and 6): ``Y_t = min(r_t^2, tau^2) / m0`` in ``[0, B]`` with
``B = tau^2 / m0``. The decrease null class contains every law of ``(Y_t)`` with ``E[Y_t | past] >= 1``; the
increase null class every law with ``E[Y_t | past] <= 1``.

Part 1 (a speed limit for every valid procedure). A stopping rule ``T`` is valid at levels ``alpha_n`` if
``P(T <= n) <= alpha_n`` for every law in the null class and every n. Then on every path ``y`` with
``T(y) <= n``, ``W_n(y) = prod_{s <= n} 1 / pi(y_s) >= 1 / alpha_n``, where ``pi(y)`` is the largest
probability that a law with the null mean constraint can put on the value y (``max_null_mass``):
``(B - 1) / (B - y)`` for ``y < 1`` on the decrease side and ``1 / y`` for ``y > 1`` on the increase side.
For a randomised rule the probability of an alarm by n on ``y`` is at most ``alpha_n W_n(y)``.
``floor_time`` is the first n with ``W_n >= A / n`` (the floor for ``alpha_n = n / A``); ``np_power`` adds
the Neyman-Pearson bound at the least favourable null law ``P*`` (the change of measure behind Lorden's
lower bound, evaluated at finite A).

Part 2 (the plug-in baseline when n, m and A grow together). ``plugin_windows`` simulates calibration
windows and returns the drift ``d(W) = m_F(tau_hat) / m_hat - 1`` of the monitor's plug-in baseline, the
inflated baselines ``m_hat (1 +- z se)``, the confidence-sequence baseline and the quantities used by the
upper bound (``first_order_constant``, ``distribution_free_constant``) and by the lower bounds
(``noise_free_alarm``, ``hold_bold_lower_bound``).

Part 3 (the decrease bet at ``lam_max (1 - eps)``). Its factor is ``eps + (1 - eps) f0(y)`` with
``f0(y) = (B - y) / (B - 1)``; ``boundary_growth`` and ``boundary_terms`` give its growth and the bounds of
the theory note.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate, optimize, special, stats

from . import robustness as R
from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from .robust import MAD_NORMAL_SCALE

# ----------------------------------------------------------------------------
# Part 1: the pathwise floor for every valid procedure
# ----------------------------------------------------------------------------


def max_null_mass(y, B, side: str = "down"):
    """Largest probability a law on ``[0, B]`` with the null mean constraint can put on the value ``y``.

    ``side="down"`` (mean at least 1): ``(B - 1) / (B - y)`` for ``y < 1``, else 1 (attained by the two-point
    law on ``{y, B}`` with mean exactly 1). ``side="up"`` (mean at most 1): ``1 / y`` for ``y > 1``, else 1
    (two-point law on ``{0, y}``). ``side="two-sided"`` (mean exactly 1): both. ``B = inf`` is allowed for
    ``side="up"`` and gives 1 for every ``y`` on the decrease side (no truncation).
    """
    y = np.asarray(y, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        if np.isinf(B):
            down = np.ones_like(y)
        else:
            down = np.where(y < 1, (B - 1.0) / (B - y), 1.0)
        up = np.where(y > 1, 1.0 / y, 1.0)
    if side == "down":
        return down
    if side == "up":
        return up
    if side == "two-sided":
        return np.where(y < 1, down, up)
    raise ValueError("side must be 'down', 'up' or 'two-sided'")


def log_floor_factor(y, B, side: str = "down"):
    """``log(1 / max_null_mass)``: the per-step increment of ``log W_n``. Equals ``log max(1, (B - y) / (B - 1))``
    on the decrease side and ``log max(1, y)`` on the increase side."""
    return -np.log(max_null_mass(y, B, side))


def floor_time(logw, A: float, alpha=None):
    """First n (1-based) with ``sum_{s <= n} logw_s >= log(1 / alpha_n)``; -1 if none within the rows.

    ``alpha=None`` means ``alpha_n = n / A``. ``logw`` has shape ``(paths, horizon)`` (or one row). By the
    floor lemma no deterministic rule valid at these levels alarms on a path before this time.
    """
    S = np.cumsum(np.atleast_2d(np.asarray(logw, float)), axis=-1)
    n = np.arange(1, S.shape[-1] + 1)
    thr = np.log(A / n) if alpha is None else -np.log(np.asarray(alpha, float))
    hit = S >= thr - 1e-12
    first = np.argmax(hit, axis=-1)
    return np.where(hit.any(axis=-1), first + 1, -1)


def collapse_floor(A: float, B: float, side: str = "down", alpha=None) -> int:
    """Floor for the most extreme path: every ``Y = 0`` (decrease side, factor ``B / (B - 1)``) or every
    ``Y = B`` (increase side, factor ``B``). Smallest n with ``n log F >= log(A / n)``."""
    F = B / (B - 1.0) if side == "down" else B
    if not (np.isfinite(F) and F > 1):
        return -1
    for n in range(1, 10**7):
        thr = np.log(A / n) if alpha is None else -np.log(alpha(n))
        if n * np.log(F) >= thr - 1e-12:
            return n
    return -1


def randomised_path_delay(log_w, A: float, horizon: int | None = None):
    """Smallest mean delay on one path over randomised rules valid at levels ``n / A`` (Theorem 2b, sharp form).

    ``log_w[j - 1] = log W_j(y)`` for j = 1, 2, ... (nondecreasing). With ``q_j = P_U(T = j | y)``, the product
    law of Theorem 2c gives the path probability ``1 / W_j`` at every j at once, so validity requires
    ``sum_{j <= k} q_j / W_j <= k / A`` for every k; and a rule that alarms only on this path with these
    ``q_j`` is valid, because every null law gives the path at most ``1 / W_j``. The linear programme
    ``min sum_j j q_j`` subject to these constraints and ``sum_j q_j = 1`` is therefore the exact minimum.
    Returns ``(mean delay, q)``; ``horizon`` defaults to three times the deterministic floor."""
    lw = np.asarray(log_w, float)
    n_det = int(np.argmax(lw + np.log(np.arange(1, lw.size + 1) / A) >= -1e-12)) + 1
    H = int(horizon or min(lw.size, 3 * n_det))
    j = np.arange(1, H + 1)
    inv_w = np.exp(-lw[:H])
    A_ub = np.tril(np.ones((H, H))) * inv_w[None, :]
    res = optimize.linprog(j.astype(float), A_ub=A_ub, b_ub=j / A, A_eq=np.ones((1, H)), b_eq=[1.0],
                           bounds=(0, None), method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    return float(res.fun), res.x


def randomised_floor_cdf(logw, A: float):
    """``E[min(1, (n / A) W_n)]`` over the rows of ``logw``, for n = 1, ..., horizon: an upper bound on
    ``P(T <= n)`` for every randomised rule valid at levels n / A (Corollary of the floor lemma)."""
    S = np.cumsum(np.atleast_2d(np.asarray(logw, float)), axis=-1)
    n = np.arange(1, S.shape[-1] + 1)
    with np.errstate(over="ignore"):
        v = np.minimum(1.0, np.exp(np.minimum(S + np.log(n / A), 50.0)))
    return v.mean(axis=0)


def np_power(s, alpha: float, cap=None) -> float:
    """Neyman-Pearson bound from Monte Carlo draws under Q.

    ``s``: log likelihood ratio ``log dQ/dP*`` of the first n observations, one value per path drawn from Q.
    Returns ``max E_Q[phi]`` over tests ``0 <= phi <= cap`` with ``E_{P*}[phi] <= alpha``, where
    ``P*``-expectations are estimated by the importance weights ``exp(-s) / N``. ``cap=None`` (``cap = 1``):
    every rule, randomised or not. ``cap`` = indicator of ``W_n >= 1 / alpha``: deterministic rules (the
    floor lemma restricts them to that set). ``cap = min(1, alpha W_n)``: randomised rules.
    """
    s = np.asarray(s, float)
    N = s.size
    c = np.ones(N) if cap is None else np.clip(np.asarray(cap, float), 0.0, 1.0)
    order = np.argsort(-s, kind="stable")
    s, c = s[order], c[order]
    w = np.exp(-s) / N
    cost = c * w
    cum = np.cumsum(cost)
    k = int(np.searchsorted(cum, alpha, side="right"))  # first k items fit fully
    power = float(c[:k].sum())
    if k < N and w[k] > 0:
        left = alpha - (cum[k - 1] if k > 0 else 0.0)
        power += float(min(c[k], max(left, 0.0) / w[k]))
    return power / N


def np_cdf(S_star, A: float, logw=None, mode: str = "np", n_max=None):
    """Upper bounds ``beta_n`` on ``P_Q(T <= n)``, n = 1, ..., horizon, for rules valid at levels n / A.

    ``S_star``: cumulative log likelihood ratios ``log dQ/dP*`` per path (shape (paths, horizon));
    ``logw``: per-step floor increments (needed for ``mode`` "det" and "rand"). ``mode="np"``: any rule;
    ``"det"``: deterministic rules (Neyman-Pearson restricted to ``W_n >= A / n``); ``"rand"``: randomised
    rules (test capped at ``min(1, (n / A) W_n)``). Stops early once ``beta_n`` exceeds 0.99999 (later values
    are set to 1).
    """
    S_star = np.atleast_2d(np.asarray(S_star, float))
    H = S_star.shape[1] if n_max is None else min(n_max, S_star.shape[1])
    W = None if logw is None else np.cumsum(np.atleast_2d(np.asarray(logw, float)), axis=1)
    beta = np.ones(H)
    for j in range(H):
        n = j + 1
        alpha = n / A
        cap = None
        if mode == "det":
            cap = (W[:, j] >= np.log(A / n) - 1e-12).astype(float)
        elif mode == "rand":
            cap = np.minimum(1.0, np.exp(np.minimum(W[:, j] + np.log(alpha), 50.0)))
        beta[j] = min(1.0, np_power(S_star[:, j], alpha, cap))
        if beta[j] > 0.99999:
            break
    return beta


def least_favourable_bet(law, side: str = "down"):
    """Growth-optimal bet ``lam*`` and ``KL_inf`` from the first-order condition ``E_Q[U / (1 + lam U)] = 0``.

    ``U = Y - 1`` (increase, ``lam <= 1``) or ``1 - Y`` (decrease, ``lam <= 1 / (B - 1)``). The root is found in
    ``log10(eps)`` with ``lam = lam_max (1 - eps)``, so that a bet within ``1e-8`` of the boundary is resolved:
    the least favourable null law ``P* = Q / (1 + lam* U)`` then has total mass 1 to about ten digits (a
    maximiser of the growth found to the same tolerance in lam can leave P* with mass 0.98, see TH5).
    Returns ``(lam*, KL_inf)``; ``(0, 0)`` if ``E_Q[U] <= 0``; ``lam* = lam_max`` if there is no interior root.
    """
    U = (law.y - 1.0) if side == "up" else (1.0 - law.y)
    top = 1.0 if side == "up" else 1.0 / (law.B - 1.0)
    if float(np.sum(law.w * U)) <= 0:
        return 0.0, 0.0

    def h(le):
        with np.errstate(divide="ignore", invalid="ignore"):
            v = U / (1.0 + top * (1.0 - 10.0**le) * U)
        return float(np.sum(np.where(law.w > 0, law.w * v, 0.0)))

    if h(-40.0) >= 0:
        lam = top
    else:
        le = optimize.brentq(h, -40.0, 0.0, xtol=1e-13, rtol=1e-14)
        lam = top * (1.0 - 10.0**le)
    with np.errstate(divide="ignore"):
        kl = float(np.sum(np.where(law.w > 0, law.w * np.log1p(lam * U), 0.0)))
    return float(lam), kl


def delay_lower_bound(beta) -> float:
    """``E[T] = sum_{n >= 0} P(T > n) >= 1 + sum_{n = 1}^{H - 1} (1 - beta_n)`` for ``beta_n >= P(T <= n)``.
    Truncating at the horizon keeps a valid lower bound."""
    beta = np.asarray(beta, float)
    return float(1.0 + np.sum(1.0 - np.clip(beta[:-1], 0.0, 1.0)))


def quantile_lower_bound(beta, p: float) -> int:
    """Smallest n with ``beta_n >= p``: every valid rule has its p-quantile of the delay at least this n."""
    beta = np.asarray(beta, float)
    idx = np.flatnonzero(beta >= p)
    return int(idx[0] + 1) if idx.size else beta.size + 1


def sr_alarm_times(y, bets, weights, direction: str = "up", thresholds=(7800.0,)):
    """First crossing of each threshold by the Shiryaev-Roberts mixture on rows of ``y`` (change at step 1).

    Returns an int array of shape ``(len(thresholds), paths)`` with 1-based alarm times, -1 if none.
    Same recursion as ``design.run_mixture`` with several thresholds at once.
    """
    y = np.atleast_2d(np.asarray(y, float))
    reps, n = y.shape
    bets = np.asarray(bets, float)
    weights = np.asarray(weights, float)
    thr = np.asarray(thresholds, float)
    sign = 1.0 if direction == "up" else -1.0
    S = np.zeros((reps, bets.size))
    out = np.full((thr.size, reps), -1)
    active = np.arange(reps)
    top = thr.max()
    for t in range(n):
        if active.size == 0:
            break
        f = 1.0 + bets[None, :] * (sign * (y[active, t] - 1.0))[:, None]
        S[active] = np.minimum(f * (S[active] + 1.0), 1e250)
        M = S[active] @ weights
        for i, a in enumerate(thr):
            new = (M >= a) & (out[i, active] < 0)
            out[i, active[new]] = t + 1
        active = active[M < top]
    return out


# ----------------------------------------------------------------------------
# Part 2: the plug-in baseline in the joint regime
# ----------------------------------------------------------------------------


def t_trunc_mean(s, nu: float):
    """``E[min(T^2, s^2)]`` for Student t with ``nu > 2`` (closed form, as in E16)."""
    s = np.asarray(s, float)
    return nu / (nu - 2) * special.betainc(1.5, nu / 2 - 1, s * s / (nu + s * s)) + 2 * s * s * stats.t.sf(s, nu)


def law_draws(law: str, rng, shape):
    """Raw returns (scale does not matter: Y is scale free)."""
    if law == "gaussian":
        return rng.standard_normal(shape)
    if law == "t3":
        return rng.standard_t(3.0, shape)
    raise ValueError(law)


def law_trunc_mean(law: str, s):
    """``m_F(s) = E_F[min(r^2, s^2)]`` for the raw law."""
    if law == "gaussian":
        return R.trunc_mean_normal(s)
    if law == "t3":
        return t_trunc_mean(s, 3.0)
    raise ValueError(law)


def law_abs_cdf(law: str, a):
    """``P(|r| <= a)`` (equivalently the distribution function of ``r^2`` at ``a^2``)."""
    a = np.asarray(a, float)
    if law == "gaussian":
        return 2 * stats.norm.cdf(a) - 1
    if law == "t3":
        return 2 * stats.t.cdf(a, 3.0) - 1
    raise ValueError(law)


def law_population(law: str, c: float = 3.0):
    """Population truncation level, baseline, B and the standard deviation of Y (``sigma_Y``) for the law."""
    dist = stats.norm if law == "gaussian" else stats.t(3.0)
    tau = c * MAD_NORMAL_SCALE * dist.ppf(0.75)
    m0 = float(law_trunc_mean(law, tau))
    pdf = (lambda a: 2 * dist.pdf(a))
    ex2 = integrate.quad(lambda a: min(a * a, tau * tau) ** 2 * pdf(a), 0, tau, limit=200)[0]
    ex2 += tau**4 * 2 * dist.sf(tau)
    var_y = ex2 / m0**2 - 1.0
    return dict(tau=float(tau), m0=m0, B=float(tau**2 / m0), sigma_y=float(np.sqrt(var_y)))


def bet_constants(B, up_grid=None, down_fractions=None):
    """Mean and largest bets of the default detector: ``lam_bar, lam_max, phi_bar, phi_max`` at ``B``."""
    lam = np.asarray(DEFAULT_UP_GRID if up_grid is None else up_grid, float)
    frac = np.asarray(DEFAULT_DOWN_FRACTIONS if down_fractions is None else down_fractions, float)
    return float(lam.mean()), float(lam.max()), float(frac.mean() / (B - 1)), float(frac.max() / (B - 1))


def first_order_constant(sigma_y: float, B: float, up_grid=None, down_fractions=None) -> float:
    """``C_F = (lam_bar + phi_bar) sigma_Y / (4 sqrt(2 pi))``: under an i.i.d. law,
    ``P(T_A <= m) <= (m / A)(1 + C_F (m + 1) / sqrt(n) + o(m / sqrt(n)))`` (Theorem 8 (c))."""
    lb, _, pb, _ = bet_constants(B, up_grid, down_fractions)
    return float((lb + pb) * sigma_y / (4 * np.sqrt(2 * np.pi)))


def b_range(c: float = 3.0):
    """Every calibration window with a positive median of |r| has ``b_lo <= B_hat <= b_hi``:
    ``b_lo = 2 a^2 / (1 + a^2)``, ``b_hi = 2 a^2`` with ``a = c * 1.4826`` (TH4, section 4)."""
    a2 = (c * MAD_NORMAL_SCALE) ** 2
    return 2 * a2 / (1 + a2), 2 * a2


def distribution_free_constants(c: float = 3.0, up_grid=None, down_fractions=None):
    """Constants of Theorem 8 (b): ``|d(W)| <= B_hat D_n`` gives, for every window,
    ``lam_bar d^+ + phi_bar d^- <= k1 D_n`` and ``x* <= k2 D_n`` (``x*`` the largest bet times drift)."""
    lo, hi = b_range(c)
    lam = np.asarray(DEFAULT_UP_GRID if up_grid is None else up_grid, float)
    frac = np.asarray(DEFAULT_DOWN_FRACTIONS if down_fractions is None else down_fractions, float)
    k1 = max(lam.mean() * hi, frac.mean() * lo / (lo - 1))
    k2 = max(lam.max() * hi, frac.max() * lo / (lo - 1))
    return float(k1), float(k2)


def distribution_free_rate_constant(c: float = 3.0) -> float:
    """``C_DF = k1 I0 / 4`` of Theorem 8 (b), with ``I0 = int_0^inf min(1, 2 e^(-2 u^2)) du
    = u0 + sqrt(2 pi) Phi_bar(2 u0)``, ``u0 = sqrt(log 2 / 2)`` (0.888; 1.79 at c = 3)."""
    k1, _ = distribution_free_constants(c)
    u0 = np.sqrt(np.log(2) / 2)
    return float(k1 * (u0 + np.sqrt(2 * np.pi) * stats.norm.sf(2 * u0)) / 4)


def dkw_log_moment(n: int, t: float) -> float:
    """Logarithm of a bound on ``E[D_n exp(t D_n)]`` for the Kolmogorov distance ``D_n`` of n i.i.d. draws
    from any continuous law: ``E[g(D_n)] = int_0^1 g'(x) P(D_n > x) dx`` with ``g(x) = x e^(t x)`` and
    ``P(D_n > x) <= min(1, 2 exp(-2 n x^2))`` (Dvoretzky-Kiefer-Wolfowitz inequality with Massart's
    constant). Computed in log space because ``t`` can be large.

    With ``u = sqrt(n) x`` and ``s = t / sqrt(n)`` the integral is ``n^(-1/2) int_0^sqrt(n) (1 + s u) e^(s u)
    min(1, 2 e^(-2 u^2)) du``. On ``[0, u0]``, ``u0 = sqrt(log 2 / 2)``, it equals ``u0 e^(s u0)`` exactly.
    Beyond ``u0`` the log integrand is concave with its maximum near ``s / 4`` and curvature ``-4``, so the
    integral is taken over ``[u0, min(sqrt(n), max(u0, s / 4) + 30)]`` with breakpoints at the maximum;
    the part beyond is below ``e^(-1800)`` of the peak. (A first version integrated over all of x in
    [0, 1] and under-resolved the tail for n >= 10^6; see docs/iterations/TH5.md.)"""
    rn = np.sqrt(n)
    s = t / rn
    u0 = np.sqrt(np.log(2) / 2)
    head = np.log(u0) + s * u0  # log of int_0^u0 (1 + s u) e^(s u) du; n >= 1 gives sqrt(n) > u0
    peak = min(max(u0, s / 4), rn)
    u_end = min(rn, max(u0, s / 4) + 30.0)

    def h(u):
        return np.log(2) - 2 * u * u + np.log1p(s * u) + s * u

    grid = np.unique(np.concatenate([np.linspace(u0, u_end, 4001), [peak, u_end]]))
    grid = grid[(grid >= u0) & (grid <= u_end)]
    hmax = float(np.max(h(grid)))
    offsets = 2.0 ** -np.arange(0, 24)
    cuts = sorted({u0, u_end, *[p for p in np.concatenate([[peak], peak - offsets, peak + offsets]) if u0 < p < u_end]})
    val = sum(integrate.quad(lambda u: np.exp(h(u) - hmax), a, b, limit=400, epsabs=0, epsrel=1e-9)[0]
              for a, b in zip(cuts[:-1], cuts[1:]))
    tail = hmax + np.log(val) if val > 0 else -np.inf
    return float(np.logaddexp(head, tail) - np.log(rn))


def dkw_moment(n: int, t: float) -> float:
    """``exp(dkw_log_moment(n, t))`` (``inf`` on overflow)."""
    with np.errstate(over="ignore"):
        return float(np.exp(dkw_log_moment(n, t)))


def distribution_free_bound(n: int, m: int, A: float, c: float = 3.0) -> float:
    """Theorem 8 (b): for every law with ``P(r = 0) = 0``,
    ``P(T_A <= m) <= (m / A)(1 + ((m + 1) / 4) k1 E[D_n exp(k2 m D_n)])``. Returned without the cap at 1
    (``inf`` when the bound overflows)."""
    k1, k2 = distribution_free_constants(c)
    with np.errstate(over="ignore"):
        return float(m / A * (1 + (m + 1) / 4 * k1 * dkw_moment(n, k2 * m)))


def kolmogorov_distance(absr_sorted, law: str):
    """``sup_x |F_n(x) - F(x)|`` for the law of |r| (rows sorted ascending)."""
    a = np.atleast_2d(absr_sorted)
    n = a.shape[1]
    u = law_abs_cdf(law, a)
    i = np.arange(1, n + 1)
    return np.maximum((i / n - u).max(axis=1), (u - (i - 1) / n).max(axis=1))


def plugin_windows(law: str, n: int, W: int, rng, c: float = 3.0, zs=(1.0, 2.0, 2.576), cs: bool = False,
                   cs_delta: float = 0.01, kolmogorov: bool = False, chunk_elems: int = 4_000_000):
    """Simulate W calibration windows of n i.i.d. returns and return the monitor's plug-in quantities.

    Keys: ``tau2`` (``tau_hat^2``), ``m_hat``, ``m_true`` (``m_F(tau_hat)``), ``d`` (plug-in drift
    ``m_true / m_hat - 1``), ``B_hat``, ``se`` (``sd(X) / sqrt(n)``), for each z the inflated baselines
    ``up_z = m_hat + z se``, ``dn_z = m_hat - z se``; with ``cs`` the confidence-sequence baseline of
    ``docs/math.md`` section 6 (``cs_tau2``, ``cs_up``, ``cs_dn``, ``cs_mtrue``); with ``kolmogorov`` the
    distance ``D_n``.
    """
    from .confseq import cs_baseline

    out = {k: [] for k in ("tau2", "m_hat", "m_true", "se", "D")}
    out.update({f"cs_{k}": [] for k in ("tau2", "up", "dn", "mtrue")})
    per = max(1, chunk_elems // n)
    for s0 in range(0, W, per):
        p = min(per, W - s0)
        r = law_draws(law, rng, (p, n))
        ar = np.abs(r)
        med = np.median(ar, axis=1)
        tau2 = (c * MAD_NORMAL_SCALE * med) ** 2
        x = np.minimum(r * r, tau2[:, None])
        out["tau2"].append(tau2)
        out["m_hat"].append(x.mean(axis=1))
        out["se"].append(x.std(axis=1, ddof=1) / np.sqrt(n))
        out["m_true"].append(np.asarray(law_trunc_mean(law, np.sqrt(tau2)), float))
        if kolmogorov:
            out["D"].append(kolmogorov_distance(np.sort(ar, axis=1), law))
        if cs:
            b = cs_baseline(r, c=c, delta_upper=cs_delta / 2, delta_lower=cs_delta / 2)
            out["cs_tau2"].append(b.tau**2)
            out["cs_up"].append(b.upper)
            out["cs_dn"].append(b.lower)
            out["cs_mtrue"].append(np.asarray(law_trunc_mean(law, b.tau), float))
    res = {k: np.concatenate(v) for k, v in out.items() if v}
    res["d"] = res["m_true"] / res["m_hat"] - 1.0
    res["B_hat"] = res["tau2"] / res["m_hat"]
    for z in zs:
        res[f"up_{z}"] = res["m_hat"] + z * res["se"]
        res[f"dn_{z}"] = res["m_hat"] - z * res["se"]
    return res


def side_bound(m_true, m_up, m_dn, tau2, m, A, up_grid=None, down_fractions=None):
    """Per-window bound of Theorem 1 (TH4) for a baseline pair ``(m_up, m_dn)``:
    ``min(1, (sum_k w_k S_m(lam_k delta_up) + sum_l v_l S_m(phi_l delta_dn)) / (2 A))``, with
    ``delta_up = (m_true / m_up - 1)^+``, ``delta_dn = (1 - m_true / m_dn)^+`` and decrease bets set by
    ``B_dn = tau2 / m_dn``. A window whose decrease side is off (``B_dn <= 1`` or ``m_dn <= 0``) contributes
    only its increase side."""
    m_true, m_up, m_dn, tau2 = (np.asarray(v, float) for v in (m_true, m_up, m_dn, tau2))
    du = np.maximum(m_true / m_up - 1.0, 0.0)
    on = (m_dn > 0) & (tau2 / np.where(m_dn > 0, m_dn, 1.0) > 1.0 + 1e-12)
    Bd = np.where(on, tau2 / np.where(m_dn > 0, m_dn, 1.0), 2.0)
    dd = np.where(on, np.maximum(1.0 - m_true / np.where(m_dn > 0, m_dn, 1.0), 0.0), 0.0)
    up, dn = R._side_sums(du, dd, m, Bd, up_grid, down_fractions)
    dn = np.where(on, dn, 0.0)
    return np.minimum(1.0, (up + dn) / (2 * A))


def noise_free_alarm(d, B_hat, m, A, up_grid=None, down_fractions=None):
    """Alarm indicator, within m steps, of the noise-free stream ``Y_t = 1 + d`` (per window). This stream
    keeps the conditional truncated mean at ``m_F(tau_hat)``, so it lies in the null class of Corollary 1."""
    d = np.asarray(d, float)
    up, dn = R._side_sums(np.maximum(d, 0), np.maximum(-d, 0), m, B_hat, up_grid, down_fractions)
    return (up >= 2 * A) | (dn >= 2 * A)


def _bold_side(d_abs, B_hat, side, up_grid, down_fractions):
    """Bet, weight, largest factor, loss factor and win probability of the bold phase for one side."""
    if side == "up":
        lam = np.asarray(DEFAULT_UP_GRID if up_grid is None else up_grid, float)
        k = int(np.argmax(lam))
        beta = np.full_like(B_hat, lam[k])
        w = 1.0 / lam.size
        F = 1 + beta * (B_hat - 1)  # factor at Y = B_hat
        Lf = 1 - beta  # factor at Y = 0
        p = (1 + d_abs) / B_hat  # P(Y = B_hat), mean 1 + d
        win_p = p
    else:
        frac = np.asarray(DEFAULT_DOWN_FRACTIONS if down_fractions is None else down_fractions, float)
        k = int(np.argmax(frac))
        beta = frac[k] / (B_hat - 1)
        w = 1.0 / frac.size
        F = 1 + beta  # factor at Y = 0
        Lf = 1 - beta * (B_hat - 1)  # factor at Y = B_hat
        p = (1 - d_abs) / B_hat  # P(Y = B_hat), mean 1 - |d|
        win_p = 1 - p
    return beta, w, F, Lf, win_p


def hold_bold_lower_bound(d, B_hat, m, A, holds=None, up_grid=None, down_fractions=None):
    """Lower bound on ``P(T_A <= m | W)`` for the hold-then-bold stream (Proposition 9), per window.

    The stream holds ``Y_t = 1 + d`` for h steps, then plays the two-point law on ``{0, B_hat}`` with mean
    ``1 + d`` (so the conditional truncated mean is ``m_F(tau_hat)`` at every step). For the side favoured by
    the sign of d, with its largest bet beta, weight w, ``K = 2 A / w``, largest factor F and
    ``rho = E[sqrt(f)] / sqrt(E[f])``:

    ``P(T_A <= m | W) >= max(1{noise-free alarm}, max_h (S_h(beta |d|) - sqrt(K S_h) rho^(m - h)) / (F K))``.

    ``holds``: candidate values of ``m - h`` (bold-phase lengths); default a grid from 0 to m.
    """
    d = np.asarray(d, float)
    B_hat = np.asarray(B_hat, float)
    best = noise_free_alarm(d, B_hat, m, A, up_grid, down_fractions).astype(float)
    L = np.unique(np.clip(np.round(np.geomspace(1, max(m, 1), 60)).astype(int), 1, m)) if holds is None else np.asarray(holds)
    for side, mask in (("up", d > 0), ("down", d < 0)):
        if not mask.any():
            continue
        da = np.abs(d[mask])
        Bh = B_hat[mask]
        beta, w, F, Lf, win_p = _bold_side(da, Bh, side, up_grid, down_fractions)
        K = 2 * A / w
        rho = (win_p * np.sqrt(F) + (1 - win_p) * np.sqrt(np.maximum(Lf, 0.0))) / np.sqrt(1 + beta * da)
        bb = np.zeros_like(da)
        for Lb in L:
            h = m - Lb
            if h < 0:
                continue
            S = R.sr_growth_sum(beta * da, h) if h > 0 else np.zeros_like(da)
            with np.errstate(over="ignore", invalid="ignore"):
                val = (S - np.sqrt(K * S) * rho**Lb) / (F * K)
            val = np.where(S >= K, 1.0, val)
            bb = np.maximum(bb, np.nan_to_num(val, nan=0.0, posinf=1.0))
        best[mask] = np.maximum(best[mask], np.clip(bb, 0.0, 1.0))
    return best


# ----------------------------------------------------------------------------
# Part 3: the decrease bet at lam_max (1 - eps)
# ----------------------------------------------------------------------------


def boundary_growth(law, eps):
    """``g(eps) = E_Q log(eps + (1 - eps) f0(Y))`` with ``f0(y) = (B - y) / (B - 1)``: the growth of the
    decrease bet ``lam = lam_max (1 - eps)``. ``law`` is a ``design.YLaw`` (quadrature nodes, weights, B)."""
    eps = np.atleast_1d(np.asarray(eps, float))
    f0 = (law.B - law.y) / (law.B - 1.0)
    with np.errstate(divide="ignore"):
        v = np.log(eps[:, None] + (1 - eps[:, None]) * f0[None, :])
    return np.where(law.w[None, :] > 0, v, 0.0) @ law.w


def boundary_terms(law):
    """``q_B = Q(Y = B)``, ``G_c = E_Q[log f0(Y); Y < B]`` (the growth at ``lam_max`` without the atom) and
    the growth-optimal ``eps*`` with its growth ``KL_inf``."""
    f0 = (law.B - law.y) / (law.B - 1.0)
    at = np.isclose(law.y, law.B, rtol=0, atol=1e-12 * law.B) | (f0 <= 0)
    qB = float(law.w[at].sum())
    Gc = float(np.sum(law.w[~at] * np.log(f0[~at])))
    res = optimize.minimize_scalar(lambda le: -float(boundary_growth(law, 10.0**le)[0]), bounds=(-16, 0),
                                   method="bounded", options={"xatol": 1e-9})
    eps_star = float(10.0**res.x)
    return dict(q_B=qB, G_c=Gc, eps_star=eps_star, kl_inf=float(boundary_growth(law, eps_star)[0]))


def block_delay_bound(logf, A: float, ns=None):
    """Block bound on the delay of one Shiryaev-Roberts bet (change at the first step, or any change time):
    ``E[T_A] <= n / p_n`` for every n, with ``p_n = P(sum_{i <= n} log f_i >= log A)``, because the statistic
    is at least the product of the factors over any block of n consecutive steps and disjoint blocks are
    independent. ``logf``: log factors on simulated paths (paths, horizon); ``p_n`` is estimated by the
    share of paths. Returns ``(min_n n / p_n, argmin n)`` (a Monte Carlo estimate of a bound)."""
    S = np.cumsum(np.atleast_2d(np.asarray(logf, float)), axis=1)
    H = S.shape[1]
    ns = np.unique(np.clip(np.round(np.geomspace(1, H, 400)).astype(int), 1, H)) if ns is None else np.asarray(ns)
    p = (S[:, ns - 1] >= np.log(A)).mean(axis=0)
    with np.errstate(divide="ignore"):
        v = np.where(p > 0, ns / np.where(p > 0, p, 1.0), np.inf)
    k = int(np.argmin(v))
    return float(v[k]), int(ns[k])


def boundary_growth_bounds(law, eps):
    """Bounds of Proposition 12 (a): ``(1 - eps) G_c + q_B log eps <= g(eps) <= G_c + q_B log eps + E[log(1 +
    eps / f0); Y < B]``. Returns ``(lower, upper)``."""
    eps = np.atleast_1d(np.asarray(eps, float))
    t = boundary_terms(law)
    f0 = (law.B - law.y) / (law.B - 1.0)
    at = np.isclose(law.y, law.B, rtol=0, atol=1e-12 * law.B) | (f0 <= 0)
    extra = np.array([float(np.sum(law.w[~at] * np.log1p(e / f0[~at]))) for e in eps])
    lo = (1 - eps) * t["G_c"] + t["q_B"] * np.log(eps)
    hi = t["G_c"] + t["q_B"] * np.log(eps) + extra
    return lo, hi
