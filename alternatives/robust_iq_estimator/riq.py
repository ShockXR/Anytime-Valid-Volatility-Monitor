"""Robust local plug-in estimators of integrated quarticity (alternative A5).

The estimators here keep the two ingredients of Bapat (2026), truncation at a
multiple of a robust scale (MAD) and local smoothing, but aggregate the local
estimates the way Jacod and Rosenbaum (2013) do: estimate the local variance
``sigma_t^2`` from truncated squared returns, square it, remove the
finite-sample bias of the square, and add up. All functions work on the IQ
scale ``IQ = int_0^T sigma_t^4 dt`` with ``delta = T / n`` and accept one path
(1-D array) or many paths (2-D array, one path per row). NaN returns (missing
minutes) are ignored.

Main pieces
-----------
* ``local_scale``: robust local scale (scaled median of ``|r|`` over a block and
  its neighbours), so the threshold follows stochastic volatility and the
  intraday pattern instead of one daily MAD.
* ``thresholds``: per-return truncation levels ``tau_i = a_n * s_i`` with either
  a fixed cut-off ``a_n = c`` (Bapat's rule, c = 3) or the slowly growing
  "universal" cut-off ``a_n = sqrt(2 log n)``.
* ``block_plugin``: the block plug-in estimator with an exact finite-sample
  correction under local Gaussianity, an estimate of ``int sigma^8`` and a
  standard error for confidence intervals.
* ``kernel_plugin``: Bapat's one-sided Gaussian kernel used for the local
  *variance*, squared and bias corrected (overlapping windows).
* ``truncated_rq_local``: fourth-power (RQ type) estimator with the same local
  threshold, to separate the effect of aggregation from that of the threshold.
* Theory helpers: truncated normal moments, asymptotic variance constants and
  the residual left by a clipped or undetected jump (see README, section 2).
"""

from __future__ import annotations

import math
import warnings

import numpy as np
from scipy import stats
from scipy.special import gamma as _gamma

from avmonitor.robust import MAD_NORMAL_SCALE, kernel_weights, mad

# ----------------------------------------------------------------------------
# Truncated normal moments
# ----------------------------------------------------------------------------


def trunc_even_moments(a: float, qmax: int = 4) -> np.ndarray:
    """``T[q] = E[Z^(2q) 1{|Z| <= a}]`` for q = 0..qmax and standard normal Z.

    Integration by parts gives ``T[q] = (2q - 1) T[q-1] - 2 a^(2q-1) phi(a)``
    with ``T[0] = 2 Phi(a) - 1``. For ``a = inf`` these are the full moments
    ``(2q - 1)!!`` (1, 1, 3, 15, 105).
    """
    T = np.empty(qmax + 1)
    if not np.isfinite(a):
        T[0] = 1.0
        for q in range(1, qmax + 1):
            T[q] = (2 * q - 1) * T[q - 1]
        return T
    phi = stats.norm.pdf(a)
    T[0] = 1.0 - 2.0 * stats.norm.sf(a)
    for q in range(1, qmax + 1):
        T[q] = (2 * q - 1) * T[q - 1] - 2.0 * a ** (2 * q - 1) * phi
    return T


def square_moments(a: float, mode: str = "clip", qmax: int = 4) -> np.ndarray:
    """Raw moments ``E[X^q]``, q = 1..qmax, of the truncated squared standard normal.

    ``mode="clip"``: ``X = min(Z^2, a^2)`` (Bapat's eq. 11 applied to returns).
    ``mode="drop"``: ``X = Z^2`` given ``|Z| <= a`` (Mancini style: large
    returns are removed and the block is re-counted).
    Index 0 of the output is ``E[X^1]``.
    """
    T = trunc_even_moments(a, qmax)
    if not np.isfinite(a):
        return T[1:].copy()
    p = 2.0 * stats.norm.sf(a)
    q = np.arange(1, qmax + 1)
    if mode == "clip":
        return T[1:] + a ** (2.0 * q) * p
    if mode == "drop":
        return T[1:] / (1.0 - p)
    raise ValueError("mode must be 'clip' or 'drop'")


def mean_power_moments(m, mu: np.ndarray):
    """``E[(S/m)^2]`` and ``E[(S/m)^4]`` for ``S`` a sum of ``m`` i.i.d. copies of X.

    ``mu`` holds the raw moments ``E[X], ..., E[X^4]``. The moments of S follow
    from its cumulants, which are ``m`` times those of X. For ``X = Z^2`` this
    gives ``(m + 2) / m`` and ``(m + 2)(m + 4)(m + 6) / m^3``.
    """
    m = np.maximum(np.asarray(m, float), 1.0)
    m1, m2, m3, m4 = mu[:4]
    c1 = m1
    c2 = m2 - m1**2
    c3 = m3 - 3 * m2 * m1 + 2 * m1**3
    c4 = m4 - 4 * m3 * m1 - 3 * m2**2 + 12 * m2 * m1**2 - 6 * m1**4
    K1, K2, K3, K4 = m * c1, m * c2, m * c3, m * c4
    es2 = K2 + K1**2
    es4 = K4 + 4 * K3 * K1 + 3 * K2**2 + 6 * K2 * K1**2 + K1**4
    return es2 / m**2, es4 / m**4


def gaussian_block_relvar(m):
    """Relative variance of the unbiased Gaussian block estimate of sigma^4: ``8 (m + 3) / (m (m + 2))``."""
    m = np.maximum(np.asarray(m, float), 1.0)
    return 8.0 * (m + 3.0) / (m * (m + 2.0))


def universal_cutoff(n: int) -> float:
    """``a_n = sqrt(2 log n)``: the expected number of Gaussian returns beyond it goes to 0."""
    return math.sqrt(2.0 * math.log(n))


# ----------------------------------------------------------------------------
# Blocks and local scale
# ----------------------------------------------------------------------------


def default_block(n: int) -> int:
    """Default block length ``ceil(sqrt(n))`` (20, 38 and 153 for n = 390, 1440, 23400)."""
    return int(math.ceil(math.sqrt(n)))


def block_edges(n: int, k: int) -> np.ndarray:
    """Edges of ``max(1, n // k)`` contiguous blocks whose lengths differ by at most one."""
    B = max(1, n // max(1, k))
    return np.linspace(0, n, B + 1).round().astype(int)


def _gather_index(edges: np.ndarray, span: int = 0):
    """Index matrix (B, L) for block j and ``span`` neighbours on each side, with a validity mask."""
    B = len(edges) - 1
    j = np.arange(B)
    starts = edges[np.clip(j - span, 0, B)]
    ends = edges[np.clip(j + 1 + span, 0, B)]
    L = int((ends - starts).max())
    idx = starts[:, None] + np.arange(L)[None, :]
    valid = idx < ends[:, None]
    return np.where(valid, idx, 0), valid


def _as2d(r):
    r = np.asarray(r, float)
    if r.ndim == 1:
        return r[None, :], True
    return r, False


def _rowwise_nanmedian(x: np.ndarray) -> np.ndarray:
    """Median over the last axis ignoring NaN, via one sort (faster than np.nanmedian with NaN)."""
    s = np.sort(x, axis=-1)  # NaN sorted last
    cnt = np.sum(np.isfinite(x), axis=-1)
    lo = np.clip((cnt - 1) // 2, 0, None)
    hi = np.clip(cnt // 2, 0, None)
    a = np.take_along_axis(s, lo[..., None], -1)[..., 0]
    b = np.take_along_axis(s, hi[..., None], -1)[..., 0]
    out = 0.5 * (a + b)
    return np.where(cnt > 0, out, np.nan)


def local_scale(r, k_s: int | None = None, span: int = 1, floor: float = 0.0, method: str = "trms",
                b: float = 4.0) -> np.ndarray:
    """Robust local standard deviation of each return (same shape as ``r``).

    The day is cut into blocks of about ``k_s`` returns (default ``ceil(sqrt(n))``).
    Each block's scale uses the block and its ``span`` neighbours on each side
    (about ``3 k_s`` returns), which follows volatility through the day.

    ``method="mad"``: ``1.4826 * median(|r|)`` over the window.
    ``method="trms"`` (default): a two-step scale. First the MAD scale s0, then
    the truncated root mean square ``sqrt(mean(r^2 1{|r| <= b s0}) / E[Z^2 1{|Z| <= b}])``
    with b = 4. The MAD has only 37% efficiency, and a noisy scale makes the
    effective cut-off ``a * s_hat / s`` random, which biases any truncation
    correction (iteration 1 in docs/iterations/ALT3.md). The second step is
    about 2.7 times more efficient and still ignores jumps beyond 4 s0.
    If more than half of the window is exactly zero (tick-size assets), the
    plain root mean square of the window is used for s0.
    """
    r2, one = _as2d(r)
    n = r2.shape[-1]
    k_s = default_block(n) if k_s is None else int(k_s)
    edges = block_edges(n, k_s)
    idx, valid = _gather_index(edges, span)
    a = np.abs(r2)[:, idx]
    a = np.where(valid, a, np.nan)
    med = _rowwise_nanmedian(a) * MAD_NORMAL_SCALE
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if np.any(~(med > 0)):
            rms = np.sqrt(np.nanmean(a**2, axis=-1))
            med = np.where(med > 0, med, rms)
        if method == "trms":
            fin = np.isfinite(a)
            keep = fin & (a <= b * med[..., None])
            num = np.sum(np.where(keep, a * a, 0.0), axis=-1)
            cnt = np.sum(fin, axis=-1)
            t1 = trunc_even_moments(b, 1)[1]
            s = np.sqrt(num / np.maximum(cnt, 1) / t1)
            med = np.where(s > 0, s, med)
        elif method != "mad":
            raise ValueError("method must be 'trms' or 'mad'")
    med = np.maximum(med, floor)
    bid = np.repeat(np.arange(len(edges) - 1), np.diff(edges))
    out = med[:, bid]
    return out[0] if one else out


def thresholds(r, rule: str = "universal", c: float = 3.0, scale: str = "local", k_s: int | None = None):
    """Per-return truncation levels ``tau_i = a_n * s_i`` and the cut-off ``a_n``.

    ``rule="fixed"``: ``a_n = c`` (Bapat 2026 uses c = 3 times the MAD).
    ``rule="universal"``: ``a_n = sqrt(2 log n)``, which grows slowly with n so
    that the truncation bias vanishes, while ``tau_n -> 0`` still removes jumps.
    ``scale="local"`` uses ``local_scale``; ``scale="global"`` uses one scaled
    MAD (median centred) for the whole day, as in eq. (11) of the paper.
    """
    r2, one = _as2d(r)
    n = r2.shape[-1]
    if rule == "fixed":
        a = float(c)
    elif rule == "universal":
        a = universal_cutoff(n)
    elif rule == "none":
        a = float("inf")
    else:
        raise ValueError("rule must be 'fixed', 'universal' or 'none'")
    if scale == "local":
        s = local_scale(r2, k_s)
    elif scale == "global":
        x = np.where(np.isfinite(r2), r2, np.nan)
        med = np.nanmedian(x, axis=-1, keepdims=True)
        s = np.nanmedian(np.abs(x - med), axis=-1, keepdims=True) * MAD_NORMAL_SCALE
        s = np.broadcast_to(s, r2.shape)
    else:
        raise ValueError("scale must be 'local' or 'global'")
    tau = a * s if np.isfinite(a) else np.full(r2.shape, np.inf)
    return (tau[0] if one else tau), a


# ----------------------------------------------------------------------------
# Estimators
# ----------------------------------------------------------------------------


def block_plugin(r, tau=None, a: float = float("inf"), k: int | None = None, T: float = 1.0,
                 mode: str = "drop", correct: bool = True) -> dict:
    """Block plug-in estimator of IQ with truncation.

    Returns are cut into blocks of about ``k`` (default ``ceil(sqrt(n))``). In
    block j, with ``m_j`` usable returns, ``S_j = sum of truncated r_i^2``:

    * ``sigma2_j = S_j / (m_j delta)`` (local variance),
    * ``sigma4_j = sigma2_j^2 / E[(S/m)^2]`` and ``sigma8_j = sigma2_j^4 / E[(S/m)^4]``,
      where the expectations are for standard normal increments truncated at
      ``a`` (``mean_power_moments``). Without truncation the first is the
      familiar ``m / (m + 2)`` correction of a squared chi-square.
    * ``IQ = sum_j (block duration) sigma4_j`` and ``IQ8 = int sigma^8`` likewise.
    * ``se`` estimates the standard deviation of IQ from
      ``sum_j duration_j^2 sigma8_j 8 (m_j + 3) / (m_j (m_j + 2))``.

    ``mode="clip"`` clips at ``tau`` (all returns stay in the count);
    ``mode="drop"`` removes returns beyond ``tau`` and re-counts the block, so a
    detected jump leaves no residual. ``tau=None`` means no truncation.
    ``correct=False`` returns the naive ``sum duration_j sigma2_j^2``.
    Output: dict with ``iq``, ``iq8``, ``se`` and ``dropped`` (share of returns
    truncated), each of shape ``(reps,)`` or scalar for 1-D input.
    """
    r2, one = _as2d(r)
    reps, n = r2.shape
    delta = T / n
    k = default_block(n) if k is None else int(k)
    edges = block_edges(n, k)
    idx, valid = _gather_index(edges, 0)
    X = (r2**2)[:, idx]
    ok = valid[None, :, :] & np.isfinite(X)
    if tau is None:
        tau2 = np.full_like(X, np.inf)
        a = float("inf")
    else:
        tau2 = (np.broadcast_to(np.asarray(tau, float), r2.shape) ** 2)[:, idx]
    beyond = ok & (X > tau2)
    if mode == "clip":
        Y = np.where(ok, np.minimum(X, tau2), 0.0)
        use = ok
    elif mode == "drop":
        use = ok & ~beyond
        Y = np.where(use, X, 0.0)
    else:
        raise ValueError("mode must be 'clip' or 'drop'")
    S = Y.sum(axis=-1)
    m = use.sum(axis=-1).astype(float)
    dur = np.diff(edges) * delta
    mu = square_moments(a, mode)
    e2, e4 = mean_power_moments(m, mu)
    with np.errstate(divide="ignore", invalid="ignore"):
        sig2 = S / (np.maximum(m, 1.0) * delta)
        if correct:
            sig4 = sig2**2 / e2
            sig8 = sig2**4 / e4
        else:
            sig4 = sig2**2
            sig8 = sig2**4
    good = m >= 1  # one return already gives an unbiased (if noisy) block estimate
    covered = (dur[None, :] * good).sum(axis=-1)
    scale = T / np.where(covered > 0, covered, np.nan)
    iq = (dur * np.where(good, sig4, 0.0)).sum(axis=-1) * scale
    iq8 = (dur * np.where(good, sig8, 0.0)).sum(axis=-1) * scale
    var = (dur**2 * np.where(good, sig8 * gaussian_block_relvar(m), 0.0)).sum(axis=-1) * scale**2
    out = dict(iq=iq, iq8=iq8, se=np.sqrt(var), dropped=beyond.sum(axis=(-1, -2)) / np.maximum(ok.sum(axis=(-1, -2)), 1))
    if one:
        out = {key: float(v[0]) for key, v in out.items()}
    return out


def kernel_plugin(r, tau=None, a: float = float("inf"), h: int = 30, kernel: str = "gaussian", T: float = 1.0) -> np.ndarray:
    """Bapat's kernel applied to the local *variance*, then squared and corrected.

    ``V_t = sum_k w_k min(r_{t-k}^2, tau^2) / delta`` with the normalised one-sided
    weights of ``avmonitor.robust.kernel_weights`` (lags 0..h). Under locally
    constant volatility and Gaussian increments,
    ``E[V_t^2] = sigma^4 (mu1^2 + (mu2 - mu1^2) sum_k w_k^2)``, where mu are
    the moments of the clipped squared normal. The estimator is
    ``T * mean_t V_t^2 / (mu1^2 + (mu2 - mu1^2) sum w^2)``. Windows that
    contain a missing return are skipped.
    """
    r2, one = _as2d(r)
    n = r2.shape[-1]
    delta = T / n
    X = r2**2
    if tau is not None:
        X = np.minimum(X, np.broadcast_to(np.asarray(tau, float), r2.shape) ** 2)
    else:
        a = float("inf")
    w = kernel_weights(h, kernel)
    L = len(w)
    M = n - L + 1
    V = np.zeros((r2.shape[0], M))
    for j in range(L):
        V += w[j] * X[:, j : j + M]
    V /= delta
    mu = square_moments(a, "clip")
    e2 = mu[0] ** 2 + (mu[1] - mu[0] ** 2) * np.sum(w**2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        iq = T * np.nanmean(V**2, axis=-1) / e2
    return float(iq[0]) if one else iq


def truncated_rq_local(r, tau, a: float, T: float = 1.0, mode: str = "drop") -> np.ndarray:
    """Fourth-power estimator with per-return threshold, corrected for truncation.

    ``IQ = sum g(r_i)^4 / (delta * E[g(Z)^4])`` with ``g`` the clip or drop map
    at ``tau_i`` and ``E[g(Z)^4]`` for a standard normal truncated at ``a``
    (3 without truncation). Uses fourth powers, so it is an RQ-type estimator
    with asymptotic variance at least ``32/3 int sigma^8``.
    """
    r2, one = _as2d(r)
    n = r2.shape[-1]
    delta = T / n
    t4 = np.broadcast_to(np.asarray(tau, float), r2.shape) ** 4
    x4 = r2**4
    ok = np.isfinite(x4)
    if mode == "drop":
        y = np.where(ok & (x4 <= t4), x4, 0.0)
        mom = trunc_even_moments(a, 2)[2]
    else:
        y = np.where(ok, np.minimum(x4, t4), 0.0)
        mom = square_moments(a, "clip")[1]
    nv = ok.sum(axis=-1)
    iq = y.sum(axis=-1) / (delta * mom) * (n / np.maximum(nv, 1))
    return float(iq[0]) if one else iq


# ----------------------------------------------------------------------------
# Asymptotic theory helpers (constant volatility, Gaussian increments)
# ----------------------------------------------------------------------------


def abs_moment(p: float) -> float:
    """``E|Z|^p`` for standard normal Z."""
    return 2 ** (p / 2) * _gamma((p + 1) / 2) / math.sqrt(math.pi)


def avar_block(k) -> float:
    """AVAR of the non-overlapping block plug-in with k returns per block, in units of ``int sigma^8``.

    ``n Var(IQ_hat) -> 8 (k + 3) / (k + 2) int sigma^8`` (T = 1). k = 1 gives
    realized quarticity (32/3) and k -> infinity the efficient bound 8.
    """
    if not np.isfinite(k):
        return 8.0
    return 8.0 * (k + 3.0) / (k + 2.0)


def avar_tpq(p: float = 4.0 / 3.0) -> float:
    """AVAR of tripower quarticity ``sum |r_i r_{i-1} r_{i-2}|^p / mu_p^3`` (3p = 4), units of ``int sigma^8``.

    The summands are 2-dependent, so the long-run variance is
    ``gamma0 + 2 gamma1 + 2 gamma2`` with ``gamma0 = mu_2p^3 - mu_p^6``,
    ``gamma1 = mu_2p^2 mu_p^2 - mu_p^6`` and ``gamma2 = mu_2p mu_p^4 - mu_p^6``.
    """
    a, b = abs_moment(2 * p), abs_moment(p)
    return (a**3 + 2 * a**2 * b**2 + 2 * a * b**4 - 5 * b**6) / b**6


def avar_qpq() -> float:
    """AVAR of quad-power quarticity ``sum |r_i ... r_{i-3}| / mu_1^4`` (3-dependent summands)."""
    b = abs_moment(1.0)
    return (1 + 2 * b**2 + 2 * b**4 + 2 * b**6 - 7 * b**8) / b**8


def avar_truncated_rq(a: float, mode: str = "clip") -> float:
    """AVAR of a fourth-power estimator truncated at a known ``a`` sd and divided by its mean factor.

    ``(E[g(Z)^8] / E[g(Z)^4]^2 - 1)``; for ``a = inf`` this is ``105/9 - 1 = 32/3``.
    """
    mu = square_moments(a, mode)
    return mu[3] / mu[1] ** 2 - 1.0


def avar_block_truncated(a: float, mode: str = "clip") -> float:
    """AVAR of the block plug-in for k -> infinity with an *oracle* (known-scale) cut-off at ``a`` sd.

    ``4 Var(X) / E[X]^2`` with X the truncated squared normal. This is below 8
    for finite ``a`` because the oracle threshold uses the true sigma. With a
    threshold estimated from the data the variance rises back to about 8
    (README, section 2.3), so confidence intervals use the Gaussian factor.
    """
    mu = square_moments(a, mode)
    return 4.0 * (mu[1] - mu[0] ** 2) / mu[0] ** 2


def clipped_jump_residual(a: float, k: float, exact: bool = True) -> float:
    """Expected extra IQ, in units of ``sigma^4 delta``, from one jump clipped to ``tau = a sigma sqrt(delta)``.

    In a block of k returns the clipped jump contributes ``a^2`` (in units of
    ``sigma^2 delta``) in place of a clipped squared normal X with moments
    ``mu1, mu2``. So ``E[S^2]`` rises by ``a^4 + 2 a^2 (k - 1) mu1 - mu2 - 2 (k - 1) mu1^2``,
    and after the block correction (division by ``E[(S/k)^2] k^2``) and the block
    weight ``k delta`` the residual is that rise divided by ``mu2 + (k - 1) mu1^2``.
    ``exact=False`` uses the untruncated moments (1, 3), which gives the simpler
    ``(a^4 + 2 a^2 (k - 1) - 2k - 1) / (k + 2)``, about ``2 (a^2 - 1)`` for large k.
    k = 1 is the fourth-power (truncated RQ) case, ``(a^4 - mu2) / mu2``.
    """
    if exact:
        mu = square_moments(a, "clip")
        mu1, mu2 = mu[0], mu[1]
    else:
        mu1, mu2 = 1.0, 3.0
    rise = a**4 + 2 * a**2 * (k - 1) * mu1 - mu2 - 2 * (k - 1) * mu1**2
    return rise / (mu2 + (k - 1) * mu1**2)


def avar_block_exact(k: int, a: float = float("inf"), mode: str = "clip") -> float:
    """AVAR of the block plug-in with k returns per block and an oracle cut-off at ``a`` sd.

    ``n Var(IQ_hat) / int sigma^8 -> k (E[(S/k)^4] / E[(S/k)^2]^2 - 1)`` for
    constant volatility; for ``a = inf`` this is ``8 (k + 3) / (k + 2)``.
    """
    e2, e4 = mean_power_moments(k, square_moments(a, mode))
    return float(k * (e4 / e2**2 - 1.0))


def undetected_jump_residual(b: float, k: float) -> float:
    """Expected extra IQ, in units of ``sigma^4 delta``, from one jump of ``b`` sd that is *not* truncated.

    ``E[(S + 2 r J + J^2)^2] - E[S^2] = J^2 sigma^2 delta (2k + 4) + J^4``, which after the
    block correction and weight gives ``2 b^2 + b^4 / (k + 2)``. For k = 1 (fourth
    powers) this is ``2 b^2 + b^4 / 3``: the plug-in damps the ``b^4`` term by ``3 / (k + 2)``.
    """
    return 2 * b**2 + b**4 / (k + 2)


def printed_cv_expectation(h: int, a: float = float("inf"), kernel: str = "gaussian") -> tuple[float, float]:
    """Expected CV(h) of eq. (17) under constant volatility, printed weights versus corrected target.

    Units of ``sigma^8 delta^4``. With q1 = E[g(Z)^4] and q2 = E[g(Z)^8] for the
    clipped normal, and leave-one-out weights ``v_k = K(-k/h)``, k = 1..h:

    * printed: ``(q2 - q1^2)(1 + sum v^2 / 9) + q1^2 (1 - sum v / 3)^2``,
    * corrected (weights renormalised to one): ``(q2 - q1^2)(1 + sum w^2)``.

    The printed version has a squared bias that grows with h because the raw
    weights do not sum to one.
    """
    v = kernel_weights(h, kernel, normalise=False)[:-1]  # lags h..1, drop lag 0 (left out)
    mu = square_moments(a, "clip")
    q1, q2 = mu[1], mu[3]
    printed = (q2 - q1**2) * (1 + np.sum(v**2) / 9.0) + q1**2 * (1 - np.sum(v) / 3.0) ** 2
    w = v / v.sum()
    corrected = (q2 - q1**2) * (1 + np.sum(w**2))
    return float(printed), float(corrected)


# ----------------------------------------------------------------------------
# Cross-validation for the bandwidth (paper eqs. 17 and 18)
# ----------------------------------------------------------------------------


def cv_curves(r_trunc: np.ndarray, hs, T: float = 1.0, kernel: str = "gaussian") -> dict:
    """CV scores for each bandwidth h, for one or many paths of already truncated returns.

    * ``printed``: eq. (17) literally, ``mean_{t>h} (r_t^4 - Q_{t,-t})^2`` with
      ``Q_{t,-t} = (1/3) sum_{k=1}^{h} K(-k/h) r_{t-k}^4`` (raw weights), divided by n.
    * ``corrected``: the target ``r_t^4`` against the leave-one-out *normalised*
      weighted mean of ``r^4`` (both have mean ``3 kappa4 sigma^4 delta^2``).
    * ``variance``: the same with squared returns (a less noisy proxy).

    The corrected scores use the common range ``t >= max(hs)`` so that they are
    comparable across h. Returns dict of arrays with shape (reps, len(hs)).
    """
    r2, _ = _as2d(r_trunc)
    reps, n = r2.shape
    x2, x4 = r2**2, r2**4
    hmax = max(hs)
    out = {key: np.empty((reps, len(hs))) for key in ("printed", "corrected", "variance")}
    for j, h in enumerate(hs):
        v = kernel_weights(h, kernel, normalise=False)[:-1]  # oldest (lag h) first, lag 1 last
        M = n - h
        P4 = np.zeros((reps, M))
        P2 = np.zeros((reps, M))
        for i in range(h):
            P4 += v[i] * x4[:, i : i + M]
            P2 += v[i] * x2[:, i : i + M]
        tgt4, tgt2 = x4[:, h:], x2[:, h:]
        out["printed"][:, j] = np.sum((tgt4 - P4 / 3.0) ** 2, axis=-1) / n
        s = v.sum()
        sl = slice(hmax - h, None)
        out["corrected"][:, j] = np.mean((tgt4[:, sl] - P4[:, sl] / s) ** 2, axis=-1)
        out["variance"][:, j] = np.mean((tgt2[:, sl] - P2[:, sl] / s) ** 2, axis=-1)
    return out


def global_threshold(r, c: float = 3.0) -> np.ndarray:
    """Bapat's eq. (11) threshold with the scaled MAD: ``c * 1.4826 * MAD(r)`` per path."""
    return c * mad(np.asarray(r, float), scaled=True)
