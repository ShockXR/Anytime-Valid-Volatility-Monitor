"""False alarm bounds when the null holds only approximately (TH4, tracker B84).

Theory and proofs: ``docs/theory_robustness.md``. Experiment: ``experiments/th4_robustness.py``.

Set-up (``docs/math.md``, sections 5 and 6). ``X_t = min(r_t^2, tau^2)``, baseline ``m0``,
``Y_t = X_t / m0`` in ``[0, B]`` with ``B = tau^2 / m0``. Increase bets ``1 + lam (Y - 1)`` with
``lam`` in ``[0, 1)``; decrease bets ``1 + phi (1 - Y)`` with ``phi`` in ``[0, 1 / (B - 1))``. Each bet
has a Shiryaev-Roberts statistic ``R_t = f_t (R_{t-1} + 1)``; each side averages its bets and the
two-sided rule alarms when one side average reaches ``2 A``.

The approximate null allows a drift ``delta_t >= 0``, known one step ahead (predictable):

* increase side: ``E[X_t | past] <= m0 (1 + delta_t)``, so ``E[f_t | past] <= 1 + lam delta_t``;
* decrease side: ``E[X_t | past] >= m0 (1 - delta_t)``, so ``E[f_t | past] <= 1 + phi delta_t``.

Main results (all proved in the theory document):

* ``sup_bound``: if ``delta_t <= delta`` always, then for every bet ``E[R_t] <= S_t(lam delta)`` with
  ``S_t(x) = sum_{j=1}^t (1 + x)^j``, and the two-sided default detector has
  ``P(false alarm within m) <= (mean_k S_m(lam_k delta_up) + mean_l S_m(phi_l delta_dn)) / (2 A)``.
  The factor over ``m / A`` is ``sup_inflation``.
* ``exogenous_bound``: if the drift path is fixed in advance or driven by something independent of
  the data's own randomness, the same holds path by path with ``(1 + x)^k`` replaced by the product
  of ``(1 + lam delta_i)`` over the window (``drift_sr_sums``).
* ``threshold_bound`` (any predictable drift, also one that responds to the data, as in GARCH):
  ``P(false alarm within m) <= K m / A + P(D*_m > K)`` for every K, where ``D*`` is the largest
  product of ``(1 + lam d_i)`` or ``(1 - phi d_i)`` over a window, with the signed drift ``d_i``
  (``drift_cusum``), so drift that averages out does not accumulate.
* ``truncated_sup_bound`` and ``mean_drift_bound``: two more bounds for any predictable drift.
* ``plugin_bound``: the plug-in baseline is a drift decided by the calibration window.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize, stats

from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from .robust import MAD_NORMAL_SCALE

DEFAULT_A = 7800.0
DAY = 390
WEEK = 1950
CLIP = 1e250  # cap for statistics that only matter below a few times 2A (avoids overflow)


# ----------------------------------------------------------------------------
# Gaussian reference values and the default bets
# ----------------------------------------------------------------------------


def trunc_mean_normal(k):
    """E[min(Z^2, k^2)] for a standard normal Z."""
    k = np.asarray(k, float)
    return 2 * stats.norm.cdf(k) - 1 - 2 * k * stats.norm.pdf(k) + 2 * k**2 * stats.norm.sf(k)


def gaussian_B(c: float = 3.0) -> float:
    """B = tau^2 / m0 for Gaussian returns with tau = c times the scaled MAD (9.045 for c = 3)."""
    tau = c * MAD_NORMAL_SCALE * stats.norm.ppf(0.75)
    return float(tau**2 / trunc_mean_normal(tau))


def default_bets(B, up_grid=None, down_fractions=None):
    """Increase bets lam_k and decrease bets phi_l = fraction_l / (B - 1) of the default detector.

    ``B`` may be an array (one value per window); then phi has shape ``B.shape + (L,)``.
    """
    lam = np.asarray(DEFAULT_UP_GRID if up_grid is None else up_grid, float)
    frac = np.asarray(DEFAULT_DOWN_FRACTIONS if down_fractions is None else down_fractions, float)
    B = np.asarray(B, float)
    if np.any(B <= 1):
        raise ValueError("B must exceed 1")
    phi = frac / (B[..., None] - 1.0)
    return lam, phi


# ----------------------------------------------------------------------------
# Theorem 1: the sup form
# ----------------------------------------------------------------------------


def sr_growth_sum(x, m):
    """S_m(x) = sum_{j=1}^m (1 + x)^j = (1 + x)((1 + x)^m - 1) / x, with S_m(0) = m (x > -1).

    This is the largest possible E[R_m] of one Shiryaev-Roberts statistic whose factors have
    conditional mean at most 1 + x (Lemma 1 of the theory document). Vectorised; may return inf.
    """
    x = np.asarray(x, float)
    m = np.asarray(m, float)
    if np.any(x <= -1):
        raise ValueError("x must exceed -1")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        safe = np.where(x == 0, 1.0, x)
        val = (1 + x) * np.expm1(m * np.log1p(x)) / safe
    return np.where(x == 0, m + 0 * x, val)


def growth_factor(x, m):
    """G_m(x) = S_m(x) / m, the inflation of one bet's bound over m / A (G_m(0) = 1)."""
    return sr_growth_sum(x, m) / np.asarray(m, float)


def growth_factor_bounds(x, m):
    """Elementary bounds for x >= 0: (1 + x)^((m + 1) / 2) <= G_m(x) <= exp(x (m + 1) / 2 + x^2 m^2 / 24)."""
    x = np.asarray(x, float)
    with np.errstate(over="ignore"):
        lo = np.exp((m + 1) / 2 * np.log1p(x))
        hi = np.exp(x * (m + 1) / 2 + x**2 * m**2 / 24)
    return lo, hi


def _side_sums(delta_up, delta_dn, m, B, up_grid, down_fractions, w_up=None, w_dn=None):
    lam, phi = default_bets(B, up_grid, down_fractions)
    du = np.asarray(delta_up, float)[..., None]
    dd = np.asarray(delta_dn, float)[..., None]
    wu = np.full(lam.size, 1 / lam.size) if w_up is None else np.asarray(w_up, float)
    wd = np.full(phi.shape[-1], 1 / phi.shape[-1]) if w_dn is None else np.asarray(w_dn, float)
    mm = np.asarray(m, float)[..., None]
    up = np.sum(wu * sr_growth_sum(lam * du, mm), axis=-1)
    dn = np.sum(wd * sr_growth_sum(phi * dd, mm), axis=-1)
    return up, dn


def sup_bound(delta_up, delta_dn, m, A=DEFAULT_A, B=None, up_grid=None, down_fractions=None, side="two-sided",
              w_up=None, w_dn=None):
    """Theorem 1: bound on P(false alarm within m steps) when the drift never exceeds the given values.

    ``side="two-sided"``: the default rule (each side at 2A); the bound is
    ``(sum_k w_k S_m(lam_k delta_up) + sum_l v_l S_m(phi_l delta_dn)) / (2 A)``.
    ``side="up"`` or ``"down"``: one side alone at level A; the bound is that side's sum over A.
    Returns the raw bound (it may exceed 1). ``B`` defaults to the Gaussian value for c = 3.
    Arrays broadcast (for example one drift and one B per calibration window).
    """
    B = gaussian_B() if B is None else B
    up, dn = _side_sums(delta_up, delta_dn, m, B, up_grid, down_fractions, w_up, w_dn)
    if side == "two-sided":
        return (up + dn) / (2 * A)
    if side == "up":
        return up / A
    if side == "down":
        return dn / A
    raise ValueError("side must be 'two-sided', 'up' or 'down'")


def sup_inflation(delta_up, delta_dn, m, B=None, up_grid=None, down_fractions=None, side="two-sided"):
    """The factor of Theorem 1 over m / A (1 at zero drift)."""
    return sup_bound(delta_up, delta_dn, m, 1.0, B, up_grid, down_fractions, side) / m


def design_drift(target, m=DAY, A=DEFAULT_A, B=None, direction="both", up_grid=None, down_fractions=None):
    """Drift delta at which the two-sided sup bound of the default detector reaches ``target``.

    ``direction="both"``: |E[Y_t | past] - 1| <= delta (delta_up = delta_dn = delta);
    ``"up"``: only upward drift (delta_dn = 0); ``"down"``: only downward drift (delta_up = 0).
    """
    B = gaussian_B() if B is None else float(B)

    def f(d):
        du = d if direction in ("both", "up") else 0.0
        dd = d if direction in ("both", "down") else 0.0
        return float(sup_bound(du, dd, m, A, B, up_grid, down_fractions)) - target

    if f(0.0) >= 0:
        return 0.0
    hi = 1.0 if direction == "down" else B - 1.0
    if f(hi) < 0:
        return np.nan
    return float(optimize.brentq(f, 0.0, hi, xtol=1e-12))


# ----------------------------------------------------------------------------
# Pathwise forms
# ----------------------------------------------------------------------------


def drift_sr_sums(d, lam, phi, horizons):
    """Pathwise sums sum_{j=1}^m prod_{i=j}^m (1 + lam d_i^+) (increase) and (1 + phi d_i^-) (decrease).

    ``d``: signed drift E[Y_t | past] - 1, shape (P, n). ``phi``: shape (L,) or (P, L). The sum obeys the
    Shiryaev-Roberts recursion s_t = rho_t (s_{t-1} + 1), so it is the detector run on the drift itself.
    Returns {m: (up (P, K), dn (P, L))} for each horizon m.
    """
    d = np.atleast_2d(np.asarray(d, float))
    P, n = d.shape
    lam = np.asarray(lam, float)
    phi = np.broadcast_to(np.asarray(phi, float), (P, np.shape(phi)[-1]))
    su = np.zeros((P, lam.size))
    sd = np.zeros((P, phi.shape[1]))
    out = {}
    hs = set(int(h) for h in horizons)
    for t in range(n):
        dp = np.maximum(d[:, t], 0.0)[:, None]
        dm = np.maximum(-d[:, t], 0.0)[:, None]
        su = np.minimum((1 + lam * dp) * (su + 1), CLIP)
        sd = np.minimum((1 + phi * dm) * (sd + 1), CLIP)
        if t + 1 in hs:
            out[t + 1] = (su.copy(), sd.copy())
    return out


def exogenous_bound(d, m, A=DEFAULT_A, B=None, up_grid=None, down_fractions=None):
    """Theorem 3 (drift fixed in advance or exogenous): per-path conditional bound, unclipped.

    ``(mean_k sum_j prod_{i=j}^m (1 + lam_k d_i^+) + mean_l sum_j prod (1 + phi_l d_i^-)) / (2 A)``.
    The proved bound on P(false alarm within m) is the mean over paths of min(1, this).
    """
    B = gaussian_B() if B is None else B
    lam, phi = default_bets(B, up_grid, down_fractions)
    up, dn = drift_sr_sums(d, lam, phi, [m])[m]
    return (up.mean(axis=1) + dn.mean(axis=1)) / (2 * A)


def drift_cusum(d, lam, phi, horizons):
    """D*_t = largest product over a window ending at t, over every bet and side, of (1 + lam d_i) or (1 - phi d_i).

    Uses the signed drift d (so a calm spell, d < 0, shrinks the increase product). CUSUM recursion
    C_t = rho_t max(C_{t-1}, 1). Returns {m: max_{t <= m} D*_t} (shape (P,)) for each horizon.
    """
    d = np.atleast_2d(np.asarray(d, float))
    P, n = d.shape
    lam = np.asarray(lam, float)
    phi = np.broadcast_to(np.asarray(phi, float), (P, np.shape(phi)[-1]))
    cu = np.zeros((P, lam.size))
    cd = np.zeros((P, phi.shape[1]))
    run = np.zeros(P)
    out = {}
    hs = set(int(h) for h in horizons)
    for t in range(n):
        dt = d[:, t][:, None]
        cu = np.minimum((1 + lam * dt) * np.maximum(cu, 1.0), CLIP)
        cd = np.minimum((1 - phi * dt) * np.maximum(cd, 1.0), CLIP)
        run = np.maximum(run, np.maximum(cu.max(axis=1), cd.max(axis=1)))
        if t + 1 in hs:
            out[t + 1] = run.copy()
    return out


def threshold_bound(dstar, m, A=DEFAULT_A):
    """Theorem 2: min over K of K m / A + P(D* > K), with P estimated by the sample ``dstar`` (one value per path).

    Returns (bound clipped at 1, the minimising K). The minimum over K is attained at a sample value.
    """
    v = np.sort(np.asarray(dstar, float))
    n = v.size
    # P(D > v_i) for each sorted value: number strictly above divided by n
    above = n - np.searchsorted(v, v, side="right")
    cand = v * m / A + above / n
    i = int(np.argmin(cand))
    return float(min(1.0, cand[i])), float(v[i])


def truncated_sup_bound(max_up, max_dn, m, A=DEFAULT_A, B=None, levels=None, up_grid=None, down_fractions=None):
    """Proposition 2b: min over thresholds (a, b) of sup_bound(a, b) + P(max_{t<=m} d_t^+ > a or max d_t^- > b).

    ``max_up``, ``max_dn``: per-path maxima of the positive and negative parts of the drift up to m.
    Thresholds are taken at sample quantiles (``levels``). Returns (bound clipped at 1, a, b).
    """
    B = gaussian_B() if B is None else B
    max_up = np.asarray(max_up, float)
    max_dn = np.asarray(max_dn, float)
    q = np.linspace(0, 1, 81) if levels is None else np.asarray(levels, float)
    a_grid = np.unique(np.concatenate([[0.0], np.quantile(max_up, q)]))
    b_grid = np.unique(np.concatenate([[0.0], np.quantile(max_dn, q)]))
    best = (1.0, np.nan, np.nan)
    sb = sup_bound(a_grid[:, None], b_grid[None, :], m, A, B, up_grid, down_fractions)
    for i, a in enumerate(a_grid):
        exceed = (max_up[:, None] > a) | (max_dn[:, None] > b_grid[None, :])
        tot = sb[i] + exceed.mean(axis=0)
        j = int(np.argmin(tot))
        if tot[j] < best[0]:
            best = (float(tot[j]), float(a), float(b_grid[j]))
    return best


def mean_drift_bound(sum_up, sum_dn, m, A=DEFAULT_A, lam_max=None, phi_max=None, B=None):
    """Proposition 4: m / A + (1 + 1 / (2A)) E[sum_{t <= m} (lam_max d_t^+ + phi_max d_t^-)] for the two-sided rule.

    ``sum_up``, ``sum_dn``: per-path sums of d^+ and d^- over the horizon (or up to the alarm, which is sharper).
    Valid for any predictable drift; the drift term does not shrink as A grows.
    """
    B = gaussian_B() if B is None else B
    lam, phi = default_bets(B)
    lam_max = float(lam.max()) if lam_max is None else lam_max
    phi_max = np.max(phi, axis=-1) if phi_max is None else phi_max
    excess = np.mean(lam_max * np.asarray(sum_up, float) + phi_max * np.asarray(sum_dn, float))
    return float(min(1.0, m / A + (1 + 1 / (2 * A)) * excess))


# ----------------------------------------------------------------------------
# The plug-in baseline (Corollary 1)
# ----------------------------------------------------------------------------


def plugin_drift(m_true, m_hat):
    """Signed drift of the plug-in: the true conditional mean at tau-hat over the estimate, minus 1."""
    return np.asarray(m_true, float) / np.asarray(m_hat, float) - 1.0


def plugin_bound(d, B_hat, m, A=DEFAULT_A, up_grid=None, down_fractions=None, per_window=False):
    """Corollary 1: mean over calibration windows of min(1, sup_bound(d^+, d^-, m, A, B_hat)).

    ``d``: plug-in drift per window (``plugin_drift``); ``B_hat = tau_hat^2 / m_hat`` per window (it sets
    the decrease bets). With ``per_window`` the per-window bounds are returned instead of their mean.
    """
    d = np.asarray(d, float)
    b = np.minimum(1.0, sup_bound(np.maximum(d, 0), np.maximum(-d, 0), m, A, B_hat, up_grid, down_fractions))
    return b if per_window else float(np.mean(b))


# ----------------------------------------------------------------------------
# Counterexamples (section 3 of the theory document)
# ----------------------------------------------------------------------------


def two_step_counterexample(lam: float, c: float):
    """E[f_1 f_2] against E[(1 + lam delta_1)(1 + lam delta_2)] when the drift responds to the data.

    Y_1 is 0 or 2 with probability 1/2 each (mean 1, so delta_1 = 0); delta_2 = c if Y_1 = 2, else 0;
    Y_2 = 1 + delta_2. Returns (E[f_1 f_2], E[product of (1 + lam delta)]); the first exceeds the
    second by lam^2 c / 2.
    """
    e_l = 0.5 * (1 - lam) + 0.5 * (1 + lam) * (1 + lam * c)
    e_d = 1 + lam * c / 2
    return e_l, e_d


def deterministic_alarm_step(y, A=DEFAULT_A, B=None, up_grid=None, down_fractions=None, n_max=100000):
    """First step at which the two-sided default detector alarms when Y_t = y at every step (0 if never by n_max).

    With y = B this is the "jackpot" burst of section 3; with y = 1 + delta it is the deterministic
    constant-drift law at which Theorem 1 is attained.
    """
    B = gaussian_B() if B is None else B
    lam, phi = default_bets(B, up_grid, down_fractions)
    fu = 1 + lam * (y - 1)
    fd = 1 + phi * (1 - y)
    su = np.zeros(lam.size)
    sd = np.zeros(phi.size)
    for t in range(1, n_max + 1):
        su = fu * (su + 1)
        sd = fd * (sd + 1)
        if max(su.mean(), sd.mean()) >= 2 * A:
            return t
    return 0


# ----------------------------------------------------------------------------
# One pass over simulated paths (used by the experiment and the tests)
# ----------------------------------------------------------------------------


def run_paths(y, d, B, A=DEFAULT_A, horizons=(DAY, WEEK), up_grid=None, down_fractions=None, check=False):
    """Run the default two-sided detector on Y = X / m0 and collect everything the bounds need.

    ``y``: (P, n) values of Y in [0, B]; ``d``: (P, n) signed drift E[Y_t | past] - 1 (predictable).
    Returns a dict with ``alarm`` (first alarm, -1 if none; 0-based step), ``alarm_corrected`` (the
    drift-corrected detector, whose factors are divided by 1 + lam d_t or 1 - phi d_t, an e-detector
    by Lemma 2), and per horizon m: ``dstar`` (max_{t<=m} D*_t), ``exo`` (the per-path conditional
    bound of Theorem 3, unclipped), ``max_up``, ``max_dn`` (maxima of d^+, d^-), ``sum_up``,
    ``sum_dn`` (sums of d^+, d^- up to min(alarm, m)). With ``check``, ``ratio_max`` is the largest
    value of N_t / (D*_t N~_t) before the alarm (Lemma 2 says it is at most 1).
    """
    y = np.atleast_2d(np.asarray(y, float))
    d = np.atleast_2d(np.asarray(d, float))
    P, n = y.shape
    lam, phi = default_bets(B, up_grid, down_fractions)
    phi = np.broadcast_to(phi, (P, phi.shape[-1]))
    K, L = lam.size, phi.shape[1]
    Su, Sd = np.zeros((P, K)), np.zeros((P, L))  # plain statistics
    Tu, Td = np.zeros((P, K)), np.zeros((P, L))  # drift-corrected statistics
    Cu, Cd = np.zeros((P, K)), np.zeros((P, L))  # CUSUM of the drift products
    Eu, Ed = np.zeros((P, K)), np.zeros((P, L))  # exogenous sums (positive parts)
    run = np.zeros(P)
    mxu, mxd = np.zeros(P), np.zeros(P)
    smu, smd = np.zeros(P), np.zeros(P)
    alarm = np.full(P, -1)
    alarm_c = np.full(P, -1)
    ratio = np.zeros(P)
    hs = sorted(int(h) for h in horizons)
    out = {m: {} for m in hs}
    thr = 2 * A
    for t in range(n):
        yt = y[:, t][:, None]
        dt = d[:, t][:, None]
        rho_u = 1 + lam * dt
        rho_d = 1 - phi * dt
        fu = 1 + lam * (yt - 1)
        fd = 1 + phi * (1 - yt)
        live = alarm < 0
        Su = np.where(live[:, None], np.minimum(fu * (Su + 1), CLIP), Su)
        Sd = np.where(live[:, None], np.minimum(fd * (Sd + 1), CLIP), Sd)
        Tu = np.minimum(fu / rho_u * (Tu + 1), CLIP)
        Td = np.minimum(fd / rho_d * (Td + 1), CLIP)
        Cu = np.minimum(rho_u * np.maximum(Cu, 1.0), CLIP)
        Cd = np.minimum(rho_d * np.maximum(Cd, 1.0), CLIP)
        dstar_t = np.maximum(Cu.max(axis=1), Cd.max(axis=1))
        run = np.maximum(run, dstar_t)
        dp, dm = np.maximum(dt[:, 0], 0.0), np.maximum(-dt[:, 0], 0.0)
        Eu = np.minimum((1 + lam * dp[:, None]) * (Eu + 1), CLIP)
        Ed = np.minimum((1 + phi * dm[:, None]) * (Ed + 1), CLIP)
        mxu, mxd = np.maximum(mxu, dp), np.maximum(mxd, dm)
        smu = smu + np.where(live, dp, 0.0)
        smd = smd + np.where(live, dm, 0.0)
        mu_, md_ = Su.mean(axis=1), Sd.mean(axis=1)
        if check:
            N = 0.5 * (mu_ + md_)
            Nt = 0.5 * (Tu.mean(axis=1) + Td.mean(axis=1))
            ratio = np.where(live, np.maximum(ratio, N / (dstar_t * Nt)), ratio)
        hit = live & (np.maximum(mu_, md_) >= thr)
        alarm[hit] = t
        hit_c = (alarm_c < 0) & (np.maximum(Tu.mean(axis=1), Td.mean(axis=1)) >= thr)
        alarm_c[hit_c] = t
        if t + 1 in out:
            o = out[t + 1]
            o["dstar"] = run.copy()
            o["exo"] = (Eu.mean(axis=1) + Ed.mean(axis=1)) / (2 * A)
            o["max_up"], o["max_dn"] = mxu.copy(), mxd.copy()
            o["sum_up"], o["sum_dn"] = smu.copy(), smd.copy()
    res = dict(alarm=alarm, alarm_corrected=alarm_c, horizons=out)
    if check:
        res["ratio_max"] = ratio
    return res
