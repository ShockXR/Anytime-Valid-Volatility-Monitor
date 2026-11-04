"""Pre-averaging, then MAD truncation, for noisy high-frequency returns (alternative A1).

Model. The observed log price is ``Y_j = X_j + u_j`` where ``X`` is the efficient
price (``dX = sigma_t dW_t`` plus jumps) and ``u_j`` is i.i.d. microstructure
noise with variance ``omega^2`` (Zhang, Mykland and Ait-Sahalia 2005; Jacod, Li,
Mykland, Podolskij and Vetter 2009). Observed returns ``r_j = Y_j - Y_{j-1}`` then
have variance ``sigma^2 delta + 2 omega^2`` and are MA(1). Bapat (2026) instead
adds noise to each return (eq. 13), which is the same as extra volatility and
cannot be removed by any estimator (``docs/estimator_issues.md``, 4b). We use
the price-noise model, where pre-averaging works.

Pre-averaging. With a window of ``k`` returns and the weight ``g(x) = min(x, 1 - x)``,

    rbar_i = sum_{j=1}^{k-1} g(j / k) r_{i+j},   i = 0, ..., n - k.

The weighted sum of the efficient increments has variance ``sigma^2 delta k psi2_k``
with ``psi2_k = (1/k) sum_j g(j/k)^2`` (about 1/12), and the weighted sum of the
noise differences has variance ``omega^2 psi1_k / k`` with
``psi1_k = k sum_j (g(j/k) - g((j-1)/k))^2`` (exactly 1 for this ``g``). So

    E[rbar_i^2] = k psi2_k delta sigma^2 + psi1_k omega^2 / k.

Dividing by ``k psi2_k delta`` gives a local variance ``V_i`` with mean
``sigma^2 + c`` where ``c = omega^2 psi1_k / (psi2_k k^2 delta)`` is the noise term.
With ``k = theta sqrt(n)`` both pieces are of the same order and the estimator
converges at rate ``n^{1/4}`` (Jacod et al. 2009). For locally Gaussian
``rbar_i`` the local quarticity is

    Q_i = rbar_i^4 / (3 (k psi2_k delta)^2) - 2 c V_i + c^2,   E[Q_i] = sigma_i^4,

because ``E[rbar^4] = 3 (E[rbar^2])^2`` and ``(s + c)^2 - 2c(s + c) + c^2 = s^2``.
The noise variance is estimated by ``omega^2 = -mean(r_j r_{j-1})`` (minus the
first autocovariance, unbiased under the model).

Truncation after pre-averaging. ``rbar_i`` is clipped at ``tau = c_mad x scaled MAD``
of the pre-averaged returns (Bapat's eq. 11 applied to ``rbar``). The clipped
second and fourth moments are then divided by the Gaussian survival factors
``kappa_2(a)`` and ``kappa_4(a)`` (``docs/math.md``, section 3), so the noise
correction above still holds. ``mode="drop"`` zeroes clipped values instead
(a Mancini-style threshold), with the matching factors; a jump that is clipped
rather than dropped keeps ``tau^4`` in every one of the ``k - 1`` windows it enters.

Monitoring. Overlapping pre-averaged returns share ``k - 2`` raw returns with
their neighbours, so ``E[rbar_t^2 | rbar_{t-1}, ...]`` is not constant even under a
constant volatility: the conditional-mean null of the e-detector fails by
construction and false alarms are not controlled. Every ``k``-th window
(``step=k``) uses disjoint price points, so the blocks are independent under the
model and the guarantee holds with ``n / k`` observations per day.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from avmonitor.robust import MAD_NORMAL_SCALE, truncated_moment_factor

# ----------------------------------------------------------------------------
# Weights and constants
# ----------------------------------------------------------------------------


def preavg_weights(k: int) -> np.ndarray:
    """``g(j / k)`` for ``j = 1, ..., k - 1`` with ``g(x) = min(x, 1 - x)``."""
    if k < 2:
        raise ValueError("window k must be at least 2")
    j = np.arange(1, k) / k
    return np.minimum(j, 1.0 - j)


def psi_constants(k: int) -> tuple[float, float]:
    """Discrete ``(psi1_k, psi2_k)``: ``k sum (g(j/k) - g((j-1)/k))^2`` and ``(1/k) sum g(j/k)^2``."""
    g = np.concatenate([[0.0], preavg_weights(k), [0.0]])
    psi1 = k * float(np.sum(np.diff(g) ** 2))
    psi2 = float(np.sum(g**2)) / k
    return psi1, psi2


def choose_k(n: int, theta: float = 0.4) -> int:
    """``k = theta sqrt(n)`` rounded to an even integer of at least 2."""
    k = int(round(theta * np.sqrt(n)))
    return max(2, k + (k % 2))


# ----------------------------------------------------------------------------
# Pre-averaged returns
# ----------------------------------------------------------------------------


def preaverage(r, k: int, step: int = 1) -> np.ndarray:
    """Pre-averaged returns ``rbar_i`` for ``i = 0, step, 2 step, ...`` up to ``n - k``.

    Works along the last axis of ``r``. ``step=1`` gives all overlapping windows,
    ``step=k`` the non-overlapping blocks used by the monitor. Window ``i`` uses
    ``r[..., i + 1 : i + k]`` (prices ``Y_i, ..., Y_{i+k-1}``), so blocks ``i`` and
    ``i + k`` share no price point.
    """
    r = np.asarray(r, float)
    g = preavg_weights(k)
    n = r.shape[-1]
    m = n - k + 1
    if m < 1:
        raise ValueError("series shorter than the pre-averaging window")
    out = np.zeros(r.shape[:-1] + (m,))
    for j, w in enumerate(g, start=1):
        out += w * r[..., j : j + m]
    return out[..., ::step]


def noise_variance(r, c: float | None = 3.0, mad_scaled: bool = True) -> np.ndarray:
    """``omega^2 = -mean(r_j r_{j-1})`` along the last axis, floored at zero.

    With ``c`` set, raw returns are first clipped at ``c x MAD`` so that a jump
    does not swamp the autocovariance.
    """
    r = np.asarray(r, float)
    if c is not None:
        med = np.median(np.abs(r), axis=-1, keepdims=True) * (MAD_NORMAL_SCALE if mad_scaled else 1.0)
        r = np.clip(r, -c * med, c * med)
    cov = np.mean(r[..., 1:] * r[..., :-1], axis=-1)
    return np.maximum(-cov, 0.0)


# ----------------------------------------------------------------------------
# Truncation of pre-averaged returns
# ----------------------------------------------------------------------------


def truncation_factors(a: float, mode: str = "clip") -> tuple[float, float]:
    """Fraction of the second and fourth Gaussian moments kept at cut-off ``a`` sd.

    ``"clip"`` replaces ``|z| > a`` by ``a`` (eq. 11 of Bapat 2026), ``"drop"`` by 0.
    """
    if mode == "clip":
        return truncated_moment_factor(a, 2), truncated_moment_factor(a, 4)
    if mode == "drop":
        tail = 2 * stats.norm.sf(a)
        return truncated_moment_factor(a, 2) - a**2 * tail, truncated_moment_factor(a, 4) - a**4 * tail / 3.0
    raise ValueError("mode must be 'clip' or 'drop'")


def truncate_preavg(rbar, c: float = 3.0, mode: str = "clip", tau=None):
    """Clip (or drop) pre-averaged returns at ``tau = c x scaled MAD`` (per row); returns ``(rt, tau)``."""
    rbar = np.asarray(rbar, float)
    if tau is None:
        tau = c * np.median(np.abs(rbar), axis=-1, keepdims=True) * MAD_NORMAL_SCALE
    tau = np.asarray(tau, float)
    if mode == "clip":
        return np.clip(rbar, -tau, tau), tau
    if mode == "drop":
        return np.where(np.abs(rbar) <= tau, rbar, 0.0), tau
    raise ValueError("mode must be 'clip' or 'drop'")


# ----------------------------------------------------------------------------
# Local and integrated estimators
# ----------------------------------------------------------------------------


def local_moments(
    rbar,
    k: int,
    delta: float,
    omega2,
    c: float | None = 3.0,
    mode: str = "clip",
    noise_correct: bool = True,
):
    """Noise-corrected local variance ``V_i`` and local quarticity ``Q_i`` per window.

    ``rbar`` are (overlapping or block) pre-averaged returns, ``omega2`` the noise
    variance per row (array or scalar). With ``c`` set the returns are first
    truncated at ``c x scaled MAD`` and the moments divided by the matching
    Gaussian factors, so that ``E[V_i] = sigma_i^2`` and ``E[Q_i] = sigma_i^4`` under
    locally constant volatility and no jumps.
    """
    rbar = np.asarray(rbar, float)
    psi1, psi2 = psi_constants(k)
    scale = k * psi2 * delta
    f2 = f4 = 1.0
    if c is not None:
        rbar, _ = truncate_preavg(rbar, c, mode)
        f2, f4 = truncation_factors(c, mode)
    cn = np.asarray(omega2, float) * psi1 / (psi2 * k**2 * delta) if noise_correct else 0.0
    if np.ndim(cn) > 0:
        cn = cn[..., None]
    v_raw = rbar**2 / (f2 * scale)
    q_raw = rbar**4 / (3 * f4 * scale**2)
    V = v_raw - cn
    Q = q_raw - 2 * cn * v_raw + cn**2
    return V, Q


def preavg_estimates(
    r,
    k: int | None = None,
    T: float = 1.0,
    c: float | None = 3.0,
    mode: str = "clip",
    noise_correct: bool = True,
    step: int = 1,
    theta: float = 0.4,
) -> dict:
    """Integrated variance and quarticity from pre-averaged, truncated returns.

    Returns a dict with ``iv``, ``iq`` (per row), ``omega2``, ``k`` and the mean
    noise term ``c_noise`` (in variance units). ``c=None`` skips truncation.
    """
    r = np.asarray(r, float)
    n = r.shape[-1]
    delta = T / n
    k = choose_k(n, theta) if k is None else k
    rbar = preaverage(r, k, step)
    omega2 = noise_variance(r, c)
    V, Q = local_moments(rbar, k, delta, omega2, c, mode, noise_correct)
    psi1, psi2 = psi_constants(k)
    return dict(
        iv=T * V.mean(axis=-1),
        iq=T * Q.mean(axis=-1),
        omega2=omega2,
        k=k,
        c_noise=omega2 * psi1 / (psi2 * k**2 * delta),
    )


def raw_truncated_estimates(r, T: float = 1.0, c: float = 3.0, bias_correct: bool = True) -> dict:
    """The plain route: MAD-clip the raw returns, then realized variance and quarticity.

    This is the normalised, kappa-corrected version of Bapat's eq. (16) (which
    equals truncated RQ up to edge effects, ``docs/math.md`` section 4), applied
    without pre-averaging. Under price noise its mean is inflated by
    ``(1 + 2 gamma^2)`` for IV and ``(1 + 2 gamma^2)^2`` for IQ, where
    ``gamma = omega / (sigma sqrt(delta))`` is the noise-to-signal ratio per return.
    """
    r = np.asarray(r, float)
    n = r.shape[-1]
    delta = T / n
    tau = c * np.median(np.abs(r), axis=-1, keepdims=True) * MAD_NORMAL_SCALE
    rt = np.clip(r, -tau, tau)
    f2, f4 = truncation_factors(c, "clip") if bias_correct else (1.0, 1.0)
    return dict(iv=np.sum(rt**2, axis=-1) / f2, iq=np.sum(rt**4, axis=-1) / (3 * delta * f4))


# ----------------------------------------------------------------------------
# Spot (local) variance for the monitor
# ----------------------------------------------------------------------------


def spot_variance_raw(r, delta: float, h: int, c: float = 3.0) -> np.ndarray:
    """Trailing-window mean of ``min(r^2, tau^2) / (kappa_2 delta)`` over the last ``h`` raw returns.

    Entry ``t`` (for ``t >= h - 1``) uses ``r[..., t - h + 1 : t + 1]``. Earlier
    entries are NaN. ``tau`` is ``c x scaled MAD`` of each row.
    """
    r = np.asarray(r, float)
    tau = c * np.median(np.abs(r), axis=-1, keepdims=True) * MAD_NORMAL_SCALE
    x = np.minimum(r**2, tau**2)
    f2, _ = truncation_factors(c, "clip")
    cs = np.cumsum(x, axis=-1)
    out = np.full(r.shape, np.nan)
    out[..., h - 1 :] = (cs[..., h - 1 :] - np.concatenate([np.zeros(r.shape[:-1] + (1,)), cs[..., :-h]], axis=-1)) / (h * f2 * delta)
    return out


def spot_variance_preavg(r, delta: float, h: int, k: int, c: float = 3.0, step: int = 1) -> np.ndarray:
    """Trailing mean of noise-corrected pre-averaged local variances over the last ``h`` raw steps.

    Window ``i`` of the pre-averaged series ends at raw index ``i + k - 1``. Entry
    ``t`` of the output averages the windows (every ``step``-th one) that end in
    ``(t - h, t]``. Entries with no complete window are NaN.
    """
    r = np.asarray(r, float)
    rbar = preaverage(r, k, 1)
    omega2 = noise_variance(r, c)
    V, _ = local_moments(rbar, k, delta, omega2, c, "clip", True)
    n = r.shape[-1]
    out = np.full(r.shape, np.nan)
    ends = np.arange(V.shape[-1]) + k - 1  # raw index where each window ends
    use = (np.arange(V.shape[-1]) % step) == 0
    for t in range(k - 1 + h - 1, n):
        sel = use & (ends > t - h) & (ends <= t)
        if sel.any():
            out[..., t] = V[..., sel].mean(axis=-1)
    return out


# ----------------------------------------------------------------------------
# Streams for the e-detector
# ----------------------------------------------------------------------------


def monitor_stream(r, k: int, step: int, c: float = 3.0, tau=None):
    """Truncated squared pre-averaged returns ``min(rbar^2, tau^2)`` and the bound ``tau^2``.

    ``tau`` defaults to ``c x scaled MAD`` of the pre-averaged returns of each
    row (in a monitor it must come from the calibration window). ``step=k``
    gives the non-overlapping block version.
    """
    rbar = preaverage(r, k, step)
    if tau is None:
        tau = c * np.median(np.abs(rbar), axis=-1, keepdims=True) * MAD_NORMAL_SCALE
    tau = np.asarray(tau, float)
    return np.minimum(rbar**2, tau**2), tau**2


def block_noise_share(k: int, gamma: float) -> float:
    """``rho_k = psi1_k gamma^2 / (psi2_k k^2)``: noise variance over signal variance in one pre-averaged block.

    The mean of the block monitor's statistic after a sd change by ``v`` is
    about ``(v^2 + rho_k) / (1 + rho_k)`` times its baseline, so ``rho_k`` is the
    dilution of the change. With ``psi1 = 1`` and ``psi2 = 1/12`` it is
    ``12 gamma^2 / k^2``.
    """
    psi1, psi2 = psi_constants(k)
    return psi1 * gamma**2 / (psi2 * k**2)


def rule_of_thumb_k(gamma: float, sd_ratio: float = 1.5, k_min: int = 4) -> int:
    """Block size that minimises the block monitor's delay for a sd change by ``sd_ratio``: ``6 gamma / sd_ratio``.

    Derivation (Gaussian blocks, truncation ignored). After the sd changes by
    ``v`` the block statistic ``Y = Z / m0`` has mean ``(v^2 + rho_k) / (1 + rho_k)``
    and standard deviation ``sqrt(2) (v^2 + rho_k) / (1 + rho_k)``, so the
    betting e-detector grows per block by about ``snr^2 / 2`` with
    ``snr = (v^2 - 1) / (sqrt(2) (v^2 + rho_k))``. The delay in seconds is about
    ``k x 2 log(A) / snr^2``, proportional to ``k (v^2 + rho_k)^2`` with
    ``rho_k = 12 gamma^2 / k^2``; its derivative in ``k`` vanishes where
    ``rho_k = v^2 / 3``, that is ``k = 6 gamma / v``. The curve is flat near the
    optimum (within 10% for ``k`` between ``k* / 1.5`` and ``1.5 k*``), so the
    default target ``v = 1.5`` (``k = 4 gamma``) serves changes of any size.
    ``gamma`` is the noise-to-signal ratio per raw return, estimated from
    ``noise_variance`` divided by the diffusive return variance
    ``mean(r^2) - 2 omega^2``. Q6 of ``run_benchmark.py`` checks the location of
    the minimum against simulated delays.
    """
    return max(k_min, int(np.ceil(6.0 * gamma / sd_ratio)))


# ----------------------------------------------------------------------------
# Benchmarks: modulated bipower variation (Podolskij and Vetter 2009) and
# pre-averaged multipower variation
# ----------------------------------------------------------------------------
#
# Podolskij and Vetter (2009, Bernoulli 15(3), 634 to 658; arXiv 0909.0827)
# average the K-step increments ``Y_{i+K} - Y_i`` inside blocks of ``L + K - 1``
# returns (``K = c1 sqrt(n)``, block length ``c2 sqrt(n)`` observations, so
# ``L = (c2 - c1) sqrt(n)`` increments per block) and form bipower or multipower
# sums of consecutive block averages. The block average ``Ybar_m`` is a
# triangular-trapezoid weighted sum of the returns in the block, so its
# variance is ``n^{-1/2} (nu1 sigma^2 + nu2 omega^2)`` with
#
#     nu1 = c1^2 (3 c2 - 4 c1) / (3 (c2 - c1)^2),   nu2 = 2 c1 / (c2 - c1)^2,
#
# for ``c2 >= 2 c1`` (derived here from the weights: the return ``j`` of the
# block is covered by ``min(j, K, L + K - j)`` increments; the noise part is the
# difference of two averages of ``K`` noise terms). These constants are checked
# by simulation in the tests, not against the printed paper, because the
# journal, arXiv and repository hosts could not be reached from this
# environment; equation numbers are therefore not quoted. Powers below 2
# (bipower 1, 1 for variance, tripower 4/3 for quarticity) make the sums robust
# to finite-activity jumps because a jump enters exactly one block average.


def _mu(p: float) -> float:
    from scipy.special import gamma as G

    return 2 ** (p / 2) * G((p + 1) / 2) / np.sqrt(np.pi)


def mbv_constants(c1: float, c2: float) -> tuple[float, float]:
    """``(nu1, nu2)`` for block length ``c2 sqrt(n)`` and increment length ``c1 sqrt(n)``, ``c2 >= 2 c1``."""
    if c2 < 2 * c1:
        raise ValueError("the closed form needs c2 >= 2 c1")
    nu1 = c1**2 * (3 * c2 - 4 * c1) / (3 * (c2 - c1) ** 2)
    nu2 = 2 * c1 / (c2 - c1) ** 2
    return nu1, nu2


def mbv_block_averages(r, K: int, block: int) -> np.ndarray:
    """Block averages of ``K``-step increments: ``Ybar_m = mean over i in block m of (Y_{i+K} - Y_i)``.

    ``r`` are returns (last axis). Block ``m`` holds observations
    ``m block, ..., (m + 1) block`` and averages the ``L = block - K + 1``
    increments that start inside it. Blocks share only the boundary price.
    """
    r = np.asarray(r, float)
    n = r.shape[-1]
    cs = np.concatenate([np.zeros(r.shape[:-1] + (1,)), np.cumsum(r, axis=-1)], axis=-1)  # prices Y_0..Y_n
    M = n // block
    L = block - K + 1
    if L < 1:
        raise ValueError("block must exceed K")
    out = np.zeros(r.shape[:-1] + (M,))
    for i in range(L):
        starts = np.arange(M) * block + i
        out += cs[..., starts + K] - cs[..., starts]
    return out / L


def mbv_block_weights(K: int, block: int) -> tuple[float, float]:
    """Exact finite-sample factors: ``Var(Ybar_m) = sig_w sigma^2 delta + noise_w omega^2``.

    Return ``j`` of the block (``j = 1, ..., block``) is covered by
    ``min(j, K, L + K - j, L)`` of the ``L = block - K + 1`` increments, so
    ``sig_w = sum_j cover_j^2 / L^2``; the noise part is the difference of two
    averages of ``K`` noise terms, ``noise_w = 2 K / L^2`` (for ``L >= K``; in
    general the exact count of boundary terms). Asymptotically
    ``sig_w ~ nu1 sqrt(n)`` and ``noise_w ~ nu2 / sqrt(n)``.
    """
    L = block - K + 1
    j = np.arange(1, block + 1)
    cover = np.minimum(np.minimum(j, K), np.minimum(L + K - j, L))
    sig_w = float(np.sum(cover.astype(float) ** 2)) / L**2
    # noise weights: the noise term u_j enters increment i with +1 at j = i + K and -1 at j = i
    w = np.zeros(block + 1)
    for i in range(L):
        w[i + K] += 1.0
        w[i] -= 1.0
    noise_w = float(np.sum(w**2)) / L**2
    return sig_w, noise_w


def mbv_estimates(r, c1: float = 0.5, c2: float = 1.0, T: float = 1.0, omega2=None, c: float | None = 3.0) -> dict:
    """Modulated bipower (variance) and tripower (quarticity) estimates with noise correction.

    With ``a = sig_w delta`` and ``b = noise_w omega^2`` from ``mbv_block_weights``
    (so that ``Var(Ybar_m) = a sigma_m^2 + b``, which equals
    ``n^{-1/2} (nu1 sigma^2 + nu2 omega^2)`` asymptotically):
    ``mean_m |Ybar_m| |Ybar_{m+1}| / mu_1^2 -> a IV / T + b`` and
    ``mean_m |Ybar_m Ybar_{m+1} Ybar_{m+2}|^{4/3} / mu_{4/3}^3 -> a^2 IQ / T + 2 a b IV / T + b^2``.
    The dict holds ``iv``, ``iq``, ``iv_se`` (standard error from the independent
    Gaussian block approximation, our derivation, not the paper's theorem),
    the block sizes and the asymptotic constants. ``omega2`` defaults to
    ``noise_variance(r, c)``.
    """
    r = np.asarray(r, float)
    n = r.shape[-1]
    s = np.sqrt(n)
    K = max(1, int(round(c1 * s)))
    block = max(2 * K, int(round(c2 * s)))
    nu1, nu2 = mbv_constants(K / s, block / s)
    sig_w, noise_w = mbv_block_weights(K, block)
    om = noise_variance(r, c) if omega2 is None else np.asarray(omega2, float)
    yb = mbv_block_averages(r, K, block)
    M = yb.shape[-1]
    ab = np.abs(yb)
    a = sig_w * T / n
    b = noise_w * om
    bip = np.mean(ab[..., 1:] * ab[..., :-1], axis=-1) / _mu(1.0) ** 2  # a sigma^2 + b
    iv = T * (bip - b) / a
    p = 4.0 / 3.0
    tri = np.mean(ab[..., 2:] ** p * ab[..., 1:-1] ** p * ab[..., :-2] ** p, axis=-1) / _mu(p) ** 3
    iq = T * (tri - 2 * a * b * (iv / T) - b**2) / a**2
    # CLT standard error for IV from independent Gaussian blocks: Var(sum |Ybar_m||Ybar_{m+1}|)
    # = (M - 1) (mu_2^2 + 2 mu_2 mu_1^2 - 3 mu_1^4) mean(v_m^2), v_m^2 estimated by Ybar_m^2 Ybar_{m+1}^2.
    theta = 1.0 + 2 * _mu(1.0) ** 2 - 3 * _mu(1.0) ** 4
    v2 = np.mean(yb[..., 1:] ** 2 * yb[..., :-1] ** 2, axis=-1)
    iv_se = T * np.sqrt(theta * v2 / (M - 1)) / (_mu(1.0) ** 2 * a)
    return dict(iv=iv, iq=iq, iv_se=iv_se, K=K, block=block, nu1=nu1, nu2=nu2, sig_w=sig_w, noise_w=noise_w, omega2=om)


def preavg_multipower_estimates(r, k: int, T: float = 1.0, omega2=None, c: float | None = 3.0) -> dict:
    """Bipower and tripower variation of pre-averaged returns taken ``k`` apart, noise corrected.

    ``rbar_i`` and ``rbar_{i+k}`` use disjoint returns, so for locally constant
    volatility ``E|rbar_i rbar_{i+k}| = mu_1^2 k psi2 delta (sigma^2 + c)`` and
    ``E|rbar_i rbar_{i+k} rbar_{i+2k}|^{4/3} = mu_{4/3}^3 (k psi2 delta)^2 (sigma^2 + c)^2``,
    with the noise term ``c = omega^2 psi1 / (psi2 k^2 delta)``. Jumps enter
    at most ``k - 1`` consecutive windows, so each product with a jump has at
    least one jump-free factor (Podolskij and Vetter 2009 style robustness).
    """
    r = np.asarray(r, float)
    n = r.shape[-1]
    delta = T / n
    psi1, psi2 = psi_constants(k)
    om = noise_variance(r, c) if omega2 is None else np.asarray(omega2, float)
    cn = om * psi1 / (psi2 * k**2 * delta)
    a = np.abs(preaverage(r, k, 1))
    s = k * psi2 * delta
    bip = np.mean(a[..., k:] * a[..., :-k], axis=-1) / (_mu(1.0) ** 2 * s)  # sigma^2 + c
    iv = T * (bip - cn)
    p = 4.0 / 3.0
    tri = np.mean(a[..., 2 * k :] ** p * a[..., k:-k] ** p * a[..., : -2 * k] ** p, axis=-1) / (_mu(p) ** 3 * s**2)
    iq = T * (tri - 2 * cn * bip + cn**2)
    return dict(iv=iv, iq=iq, omega2=om, c_noise=cn)
