"""Robust and econometric benchmarks for the volatility monitor (experiment E21).

The independent referee (``docs/review/referee_report_paper.md``, M5) asked for benchmarks that are not
calibrated on Gaussian data only. This module adds three families; none of them changes any existing code.

1. Classical detectors on the e-detector's own statistic ``X_t = min(r_t^2, tau^2)``, with
   ``tau = c x scaled MAD`` and ``m0`` = mean of ``X`` over the calibration window (the plug-in
   baseline of ``detectors.calibrate``), written ``Y_t = X_t / m0``:

   * Page's two-sided CUSUM, ``G+_t = max(0, G+_{t-1} + Y_t - 1 - k+)`` and
     ``G-_t = max(0, G-_{t-1} + 1 - Y_t - k-)``. The reference values come from the Gaussian
     likelihood ratio of a change of the variance by a factor ``theta``: for chi-square data that ratio is
     an increasing affine function of ``Y - theta log(theta) / (theta - 1)``, so
     ``k = theta log(theta) / (theta - 1) - 1`` for an increase and ``1 - theta log(theta) / (theta - 1)``
     for a decrease (``gaussian_reference_value``).
   * The Shiryaev-Roberts betting statistic of the e-detector itself (same bets, same recursion
     ``M_t = mean_k S_t^k``, ``S_t = f_t (S_{t-1} + 1)``), but with a threshold calibrated by
     simulation or by a bootstrap instead of the level ``2 A`` that Ville's inequality gives.

   Thresholds are calibrated (a) by simulation on Gaussian data, or (b) for each path by a bootstrap of
   its own calibration window (``bootstrap_max_stats``): resample the window's ``X`` values, recompute
   ``m0`` on a resampled calibration window, run the rule on a resampled monitoring stream, and take
   the threshold at the target rate. Neither calibration gives a finite-sample guarantee.

2. Bayesian online change-point detection (Adams and MacKay 2007, Algorithm 1) with a Student t
   observation model of fixed degrees of freedom ``nu`` and a discrete prior on a log-spaced grid of
   variances (``GridBOCPD``). Inside a run the returns are i.i.d. ``t_nu`` with variance ``V`` (scale
   ``sqrt(V (nu - 2) / nu)``); the prior on ``V`` is the inverse gamma prior of
   ``competitors.BOCPDVariance`` (shape ``alpha0``, mean ``prior_var``), discretised on the grid. With a
   finite grid the joint posterior of (run length, ``V``) is updated exactly. Run lengths of at least
   ``w`` (the alarm rule's short-run window) are lumped into one bucket; because the hazard is constant
   and the likelihood does not depend on the run length given ``V``, the lumped joint posterior is
   updated exactly too, so unlike the capped filter of ``competitors`` no approximation enters the
   alarm statistic. ``nu=None`` gives a Gaussian observation model on the same grid (used to check the
   implementation against ``competitors.BOCPDVariance``). This is the Student t version allowed by
   the task, not the beta-divergence BOCPD of Knoblauch, Jewson and Damoulas (2018).

3. A sequential CUSUM monitoring procedure for the variance with a training sample, of the type of
   Horvath, Kokoszka and Zhang (2006, Econometric Theory 22(3)) (``hkz_paths``). The training sample
   is the calibration window ``Z_1, ..., Z_m`` (``Z = r^2``, or ``Z = X`` for the robust variant), and
   the detector is ``Q(m, k) = sum_{i=m+1}^{m+k} Z_i - k Zbar_m``, standardised by ``sigma_m sqrt(m)``
   with ``sigma_m^2`` a Bartlett-kernel long-run variance of the training sample (bandwidth
   ``floor(m^(1/3))``). It alarms when ``|Q(m, k)| / (sigma_m sqrt(m))`` reaches
   ``c g(m, k)`` with the boundary ``g(m, k) = (1 + k / m) (k / (m + k))^gamma``, ``0 <= gamma < 1/2``
   (Horvath, Huskova, Kokoszka and Steinebach 2004). As m grows the probability of an alarm at any
   time under no change tends to ``P(sup_{0<t<1} |W(t)| / t^gamma > c)``. The published procedure of
   Horvath, Kokoszka and Zhang could not be read here (publisher pages blocked); the form above, the
   boundary and its limit were verified from secondary sources only (``docs/report_parts/E21.md``).
   ``boundary="csw"`` uses instead the boundary of the OLS-CUSUM monitoring in the R package
   ``strucchange`` (``R/monitoring.R``, source read): with ``x = (m + k) / m`` the process must stay
   below ``sqrt(x (x - 1) (c^2 + log(x / (x - 1))))``, and ``c`` solves
   ``2 (1 - Phi(c) + c phi(c)) = alpha``.

Every alarm rule here returns, per path, a statistic for each side that is "larger means stronger
evidence"; thresholds are applied with ``first_alarm`` (simulation) or ``restart_events`` (a streaming
series with recalibration after each alarm, as ``monitor.VolatilityMonitor`` does).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats
from scipy.special import gammaln

from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from .robust import MAD_NORMAL_SCALE

# ----------------------------------------------------------------------------
# The truncated statistic with a plug-in baseline
# ----------------------------------------------------------------------------


@dataclass
class TruncatedInputs:
    y: np.ndarray  # X / m0 on the monitored part, (reps, n)
    ycal: np.ndarray  # X / m0 on the calibration window, (reps, n_cal) (NaN where the return is missing)
    B: np.ndarray  # tau^2 / m0 per path, (reps,)
    tau2: np.ndarray  # (reps,)
    m0: np.ndarray  # (reps,)


def truncated_inputs(r: np.ndarray, cal: np.ndarray, c: float = 3.0) -> TruncatedInputs:
    """``Y = min(r^2, tau^2) / m0`` with ``tau = c x`` scaled MAD and ``m0`` from each path's calibration window.

    Same plug-in baseline as ``detectors.calibrate(c=c)`` (band 0, no margin), vectorised over paths.
    Missing returns stay missing (NaN).
    """
    r = np.atleast_2d(np.asarray(r, float))
    cal = np.atleast_2d(np.asarray(cal, float))
    scale = np.nanmedian(np.abs(cal), axis=1) * MAD_NORMAL_SCALE
    tau2 = (c * scale) ** 2
    xc = np.minimum(cal**2, tau2[:, None])
    m0 = np.nanmean(xc, axis=1)
    y = np.minimum(r**2, tau2[:, None]) / m0[:, None]
    return TruncatedInputs(y, xc / m0[:, None], tau2 / m0, tau2, m0)


# ----------------------------------------------------------------------------
# Page's CUSUM on Y
# ----------------------------------------------------------------------------


def gaussian_reference_value(ratio: float) -> float:
    """Reference value k of Page's CUSUM on ``Y = X / m0`` for a change of the sd by ``ratio``.

    For ``r ~ N(0, s^2)`` against ``N(0, theta s^2)``, ``theta = ratio^2``, the log likelihood ratio of one
    return is ``(1 - 1/theta) / 2 (z^2 - theta log(theta) / (theta - 1))`` with ``z = r / s``. Its CUSUM is
    Page's CUSUM on ``z^2`` with reference value ``theta log(theta) / (theta - 1) - 1`` above 1 (increase)
    or ``1 - theta log(theta) / (theta - 1)`` below 1 (decrease). Applied to the truncated ``Y``.
    """
    th = float(ratio) ** 2
    if th == 1.0:
        raise ValueError("ratio must differ from 1")
    v = th * math.log(th) / (th - 1.0) - 1.0
    return v if th > 1 else -v


def page_cusum_paths(y: np.ndarray, k_up: float, k_dn: float):
    """Two-sided Page CUSUM on rows of ``y``: ``(G_up, G_dn)``, both ``(reps, n)``.

    ``G_t = max(0, G_{t-1} + d_t)`` with ``G_0 = 0`` equals ``S_t - min(0, min_{j <= t} S_j)`` for the
    partial sums ``S`` of ``d`` (vectorised). Missing values add nothing (``d = 0``) and are marked NaN
    in the output, so that they never trigger an alarm.
    """
    y = np.atleast_2d(np.asarray(y, float))
    ok = np.isfinite(y)
    full = bool(ok.all())
    out = []
    for sign, k in ((1.0, k_up), (-1.0, k_dn)):
        d = sign * (y - 1.0) - k
        if not full:
            d = np.where(ok, d, 0.0)
        S = np.cumsum(d, axis=1)
        G = S - np.minimum(np.minimum.accumulate(S, axis=1), 0.0)
        out.append(G if full else np.where(ok, G, np.nan))
    return out[0], out[1]


# ----------------------------------------------------------------------------
# The e-detector's Shiryaev-Roberts betting statistic, vectorised (log domain)
# ----------------------------------------------------------------------------


def sr_log_paths(y: np.ndarray, lam: np.ndarray, side: str, chunk: int = 256) -> np.ndarray:
    """``log M_t`` with ``M_t = mean_k S^k_t``, ``S^k_t = f^k_t (S^k_{t-1} + 1)``, ``S_0 = 0``.

    ``f = 1 + lam (Y - 1)`` for ``side="up"`` and ``1 + lam (1 - Y)`` for ``"down"``. ``lam`` is ``(K,)``
    (same bets on every path) or ``(reps, K)`` (per path, as for the decrease side, whose largest bet is
    ``1 / (B - 1)``). Uses ``S_t = sum_{j <= t} exp(L_t - L_{j-1})`` with ``L`` the partial sums of
    ``log f``, accumulated with ``logaddexp`` (no overflow). Missing values (NaN) leave the statistic
    unchanged, as in ``detectors.betting_edetector`` and ``monitor.VolatilityMonitor``; the output is
    NaN there.
    """
    y = np.atleast_2d(np.asarray(y, float))
    reps, n = y.shape
    lam = np.asarray(lam, float)
    lam2 = np.broadcast_to(lam if lam.ndim == 2 else lam[None, :], (reps, lam.shape[-1]))
    K = lam2.shape[1]
    sign = 1.0 if side == "up" else -1.0
    out = np.full((reps, n), np.nan)
    for i0 in range(0, reps, chunk):
        sl = slice(i0, i0 + chunk)
        yy = y[sl]
        ok = np.isfinite(yy)
        g = np.where(ok, sign * (yy - 1.0), 0.0)
        with np.errstate(divide="ignore"):
            lf = np.log1p(lam2[sl][:, None, :] * g[:, :, None])  # (p, n, K)
        L = np.cumsum(lf, axis=1)
        # one "+1" per finite step: the term exp(-L_{j-1}) enters only at steps j with an observation
        a = np.where(ok[:, :, None], -np.concatenate([np.zeros_like(L[:, :1]), L[:, :-1]], axis=1), -np.inf)
        acc = np.logaddexp.accumulate(a, axis=1)
        logS = L + acc
        m = logS.max(axis=2, keepdims=True)
        with np.errstate(invalid="ignore"):
            logM = m[:, :, 0] + np.log(np.exp(logS - m).sum(axis=2)) - math.log(K)
        out[sl] = np.where(ok, logM, np.nan)
    return out


def down_bets(B: np.ndarray, fractions=DEFAULT_DOWN_FRACTIONS) -> np.ndarray:
    """Per-path decrease bets ``fractions / (B - 1)`` (``(reps, K)``); zero bets where ``B <= 1``."""
    B = np.asarray(B, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.asarray(fractions, float)[None, :] / (B[:, None] - 1.0)
    return np.where(B[:, None] > 1.0, lam, 0.0)


# ----------------------------------------------------------------------------
# Alarms from side statistics
# ----------------------------------------------------------------------------


def first_alarm(up: np.ndarray, dn: np.ndarray, h_up, h_dn):
    """First alarm per path (index, -1 if none) and its direction (+1 up, -1 down, 0 none).

    ``h_up`` / ``h_dn`` are scalars or per-path arrays. Up wins a tie at the same step.
    """
    up = np.nan_to_num(np.atleast_2d(up), nan=-np.inf)
    dn = np.nan_to_num(np.atleast_2d(dn), nan=-np.inf)
    hu = np.broadcast_to(np.asarray(h_up, float), (up.shape[0],))[:, None]
    hd = np.broadcast_to(np.asarray(h_dn, float), (dn.shape[0],))[:, None]
    cu, cd = up >= hu, dn >= hd
    n = up.shape[1]
    tu = np.where(cu.any(1), cu.argmax(1), n)
    td = np.where(cd.any(1), cd.argmax(1), n)
    t = np.minimum(tu, td)
    direction = np.where(t == n, 0.0, np.where(tu <= td, 1.0, -1.0))
    return np.where(t == n, -1, t), direction


def side_maxima(up: np.ndarray, dn: np.ndarray, horizon: int):
    """Maximum of each side's statistic over the first ``horizon`` steps (``-inf`` if all missing)."""
    u = np.nan_to_num(np.atleast_2d(up)[:, :horizon], nan=-np.inf).max(1)
    d = np.nan_to_num(np.atleast_2d(dn)[:, :horizon], nan=-np.inf).max(1)
    return u, d


def calibrate_two_sided(max_up: np.ndarray, max_dn: np.ndarray, alpha: float):
    """Per-side thresholds at level ``alpha / 2`` each, from simulated maxima over the target horizon.

    Returns ``(h_up, h_dn, two_sided_rate)``: the ``1 - alpha / 2`` quantiles and the share of the
    simulated paths on which either side reaches its threshold (just below ``alpha`` because the two
    sides rarely alarm on the same path).
    """
    hu = float(np.quantile(max_up, 1 - alpha / 2))
    hd = float(np.quantile(max_dn, 1 - alpha / 2))
    return hu, hd, float(np.mean((max_up >= hu) | (max_dn >= hd)))


# ----------------------------------------------------------------------------
# Bootstrap thresholds from each path's calibration window
# ----------------------------------------------------------------------------


def bootstrap_max_stats(xcal: np.ndarray, tau2, rules: dict, horizon: int, n_boot: int, rng: np.random.Generator,
                        chunk: int = 20):
    """Bootstrap distribution of each side's maximum over ``horizon`` steps, per path and rule.

    ``xcal`` (``(reps, n_cal)``) holds each path's truncated squares ``X`` from its calibration window
    (NaN for missing returns, dropped) and ``tau2`` (``(reps,)``) its truncation level. For every
    bootstrap replicate: draw ``n_cal`` values with replacement and take their mean as ``m0*`` (the
    rule's own calibration), draw ``horizon`` further values as the monitored stream, form
    ``Y* = X* / m0*`` and run each rule on it. ``tau`` is kept at the window's value (``X`` is already
    truncated), so the bootstrap includes the error of ``m0`` but not that of ``tau``. ``rules`` maps a
    name to ``("sr", {})`` (the e-detector's bets, decrease bets from ``B* = tau^2 / m0*``; maxima of
    ``log M``) or ``("cusum", {"k_up": .., "k_dn": ..})`` (Page). All rules see the same replicates.
    Returns ``{name: (max_up, max_dn)}``, each ``(reps, n_boot)``.
    """
    xcal = np.atleast_2d(np.asarray(xcal, float))
    reps = xcal.shape[0]
    out = {nm: (np.empty((reps, n_boot)), np.empty((reps, n_boot))) for nm in rules}
    gu = np.asarray(DEFAULT_UP_GRID, float)
    fd = np.asarray(DEFAULT_DOWN_FRACTIONS, float)
    t2all = np.broadcast_to(np.asarray(tau2, float), (reps,))
    for i0 in range(0, reps, chunk):
        p = min(chunk, reps - i0)
        ys = np.empty((p, n_boot, horizon))
        Bs = np.empty((p, n_boot))
        for j in range(p):
            v = xcal[i0 + j]
            v = v[np.isfinite(v)]
            nc = v.size
            m0s = v[rng.integers(0, nc, (n_boot, nc))].mean(axis=1)
            ys[j] = v[rng.integers(0, nc, (n_boot, horizon))] / m0s[:, None]
            Bs[j] = t2all[i0 + j] / m0s
        for nm, (kind, prm) in rules.items():
            ou, od = out[nm]
            if kind == "sr":
                Su = np.zeros((p, n_boot, gu.size))
                Sd = np.zeros((p, n_boot, fd.size))
                with np.errstate(divide="ignore", invalid="ignore"):
                    lam_d = np.where(Bs[:, :, None] > 1.0, fd[None, None, :] / (Bs[:, :, None] - 1.0), 0.0)
                mu = np.zeros((p, n_boot))
                md = np.zeros((p, n_boot))
                with np.errstate(over="ignore"):
                    for t in range(horizon):
                        yt = ys[:, :, t, None]
                        Su = (1.0 + gu * (yt - 1.0)) * (Su + 1.0)
                        Sd = (1.0 + lam_d * (1.0 - yt)) * (Sd + 1.0)
                        np.maximum(mu, Su.mean(2), out=mu)
                        np.maximum(md, Sd.mean(2), out=md)
                ou[i0 : i0 + p] = np.log(np.maximum(mu, 1e-300))
                od[i0 : i0 + p] = np.log(np.maximum(md, 1e-300))
            elif kind == "cusum":
                Gu, Gd = page_cusum_paths(ys.reshape(p * n_boot, horizon), prm["k_up"], prm["k_dn"])
                ou[i0 : i0 + p] = Gu.max(1).reshape(p, n_boot)
                od[i0 : i0 + p] = Gd.max(1).reshape(p, n_boot)
            else:
                raise ValueError("rule kind must be 'sr' or 'cusum'")
    return out


def bootstrap_thresholds(max_u: np.ndarray, max_d: np.ndarray, alpha: float):
    """Per-path thresholds: the ``1 - alpha / 2`` quantile of each side's bootstrap maxima."""
    return np.quantile(max_u, 1 - alpha / 2, axis=1), np.quantile(max_d, 1 - alpha / 2, axis=1)


# ----------------------------------------------------------------------------
# Sequential CUSUM monitoring with a training sample (Horvath, Kokoszka and Zhang type)
# ----------------------------------------------------------------------------


def bartlett_lrv(z: np.ndarray, q: int | None = None) -> np.ndarray:
    """Bartlett-kernel long-run variance of each row of ``z`` (NaN dropped), bandwidth ``floor(n^(1/3))``.

    ``s^2 = g_0 + 2 sum_{l=1}^{q} (1 - l / (q + 1)) g_l`` with ``g_l`` the lag-l sample autocovariance
    (divisor n, rows centred at their own mean). Floored at ``g_0 / 10`` so that it stays positive.
    """
    z = np.atleast_2d(np.asarray(z, float))
    out = np.empty(z.shape[0])
    for i, row in enumerate(z):
        v = row[np.isfinite(row)]
        n = v.size
        qq = int(math.floor(n ** (1.0 / 3.0))) if q is None else int(q)
        d = v - v.mean()
        g0 = float(d @ d) / n
        s = g0
        for lag in range(1, min(qq, n - 1) + 1):
            s += 2.0 * (1.0 - lag / (qq + 1.0)) * float(d[lag:] @ d[:-lag]) / n
        out[i] = max(s, g0 / 10.0)
    return out


def hkz_paths(zcal: np.ndarray, zmon: np.ndarray, gamma: float = 0.25, boundary: str = "hhks",
              q: int | None = None) -> np.ndarray:
    """Standardised CUSUM monitoring process ``D_k`` (``(reps, n)``); alarm when ``|D_k| >= c``.

    ``Q(m, k) = sum_{i <= k} Z_{m+i} - k Zbar_m`` over the finite monitored values (missing values are
    skipped and do not advance k), ``s = sigma_m sqrt(m)`` with the Bartlett long-run sd ``sigma_m`` of
    the training sample and ``m`` its number of finite values.

    * ``boundary="hhks"``: ``D_k = Q / (s g(m, k))``, ``g = (1 + k/m) (k / (m + k))^gamma``.
    * ``boundary="csw"``: ``D_k = sign(Q) sqrt(max(0, (Q / s)^2 / (x (x - 1)) - log(x / (x - 1))))`` with
      ``x = (m + k) / m``, so that ``|D_k| >= c`` is exactly ``|Q| / s >= sqrt(x (x-1) (c^2 + log(x/(x-1))))``.
    The sign of ``D_k`` gives the direction (positive: the variance has risen).
    """
    zcal = np.atleast_2d(np.asarray(zcal, float))
    zmon = np.atleast_2d(np.asarray(zmon, float))
    if not 0 <= gamma < 0.5:
        raise ValueError("gamma must lie in [0, 1/2)")
    m = np.isfinite(zcal).sum(1).astype(float)
    zbar = np.nanmean(zcal, axis=1)
    sig = np.sqrt(bartlett_lrv(zcal, q))
    ok = np.isfinite(zmon)
    k = np.cumsum(ok, axis=1).astype(float)
    Q = np.cumsum(np.where(ok, zmon, 0.0), axis=1) - k * zbar[:, None]
    s = (sig * np.sqrt(m))[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        if boundary == "hhks":
            g = (1.0 + k / m[:, None]) * (k / (m[:, None] + k)) ** gamma
            D = Q / (s * g)
        elif boundary == "csw":
            x = (m[:, None] + k) / m[:, None]
            v = (Q / s) ** 2 / (x * (x - 1.0)) - np.log(x / (x - 1.0))
            D = np.sign(Q) * np.sqrt(np.maximum(v, 0.0))
        else:
            raise ValueError("boundary must be 'hhks' or 'csw'")
    return np.where(ok & (k > 0), D, np.nan)


def wiener_sup_samples(gamma: float, n_paths: int = 20000, t_max: float = 1.0, rng: np.random.Generator | None = None,
                       n_geom: int = 4000, n_lin: int = 20000, t_min: float = 1e-8, t_switch: float = 1e-2,
                       chunk: int = 500) -> np.ndarray:
    """Samples of ``sup_{0<t<=t_max} |W(t)| / t^gamma`` for a standard Wiener process W.

    W is sampled exactly (independent normal increments) at ``n_geom`` geometrically spaced times from
    ``t_min`` to ``t_switch``, which resolves small t where ``t^gamma`` matters, and at ``n_lin`` equally
    spaced times from ``t_switch`` to ``t_max``. The supremum over the grid slightly underestimates the
    supremum over the interval (for gamma = 0 the simulated 95% point is checked against the exact law of
    ``sup |W|`` in ``tests/test_benchmarks_robust.py``).
    """
    rng = np.random.default_rng(0) if rng is None else rng
    t = np.concatenate([np.geomspace(t_min, t_switch, n_geom, endpoint=False), np.linspace(t_switch, t_max, n_lin)])
    dt = np.diff(np.concatenate([[0.0], t]))
    sups = []
    for i in range(0, n_paths, chunk):
        p = min(chunk, n_paths - i)
        W = np.cumsum(rng.standard_normal((p, t.size)) * np.sqrt(dt), axis=1)
        sups.append((np.abs(W) / t**gamma).max(1))
    return np.concatenate(sups)


def sup_abs_bm_cdf(x: float, terms: int = 50) -> float:
    """Exact ``P(sup_{0<=t<=1} |W(t)| <= x)`` (series of Feller): ``(4/pi) sum_k (-1)^k / (2k+1) exp(-(2k+1)^2 pi^2 / (8 x^2))``."""
    k = np.arange(terms)
    return float(4 / np.pi * np.sum((-1.0) ** k / (2 * k + 1) * np.exp(-((2 * k + 1) ** 2) * np.pi**2 / (8 * x * x))))


def wiener_sup_crit(gamma: float, alpha: float = 0.05, n_paths: int = 20000, rng: np.random.Generator | None = None) -> float:
    """Simulated ``1 - alpha`` quantile of ``sup_{0<t<1} |W(t)| / t^gamma`` (the limit of the HHKS boundary)."""
    return float(np.quantile(wiener_sup_samples(gamma, n_paths, 1.0, rng), 1 - alpha))


def csw_crit(alpha: float = 0.05) -> float:
    """``c`` with ``2 (1 - Phi(c) + c phi(c)) = alpha``: the critical value of strucchange's OLS-CUSUM monitoring."""
    f = lambda c: 2.0 * (1.0 - stats.norm.cdf(c) + c * stats.norm.pdf(c)) - alpha  # noqa: E731
    return float(optimize.brentq(f, 0.5, 10.0))


# ----------------------------------------------------------------------------
# BOCPD with a Student t observation model on a variance grid
# ----------------------------------------------------------------------------


def log_variance_grid(lo: float, hi: float, step: float = 0.1) -> np.ndarray:
    """Equally spaced log-variances from ``log(lo)`` to ``log(hi)`` (inclusive, step in log variance)."""
    a, b = math.log(lo), math.log(hi)
    n = int(round((b - a) / step)) + 1
    return np.linspace(a, b, n)


class GridBOCPD:
    """Adams and MacKay (2007) run-length filter with a variance grid and a Student t (or Gaussian) likelihood.

    State: the joint posterior ``J[r, j] = P(run length r, V = V_j | x_1..x_t)`` for ``r = 0, ..., w - 1``
    (exact run lengths) and ``r = w`` (all runs of at least ``w`` observations), normalised to 1. One
    step with observation x: ``P = J f(x | V_j)``; the new run (no observation yet, prior ``pi``) gets
    ``H sum P``; run r moves to r + 1 with ``(1 - H) P[r]``; runs ``w - 1`` and ``>= w`` both move into the
    lumped bucket. This is the exact recursion for a constant hazard H, a prior supported on the grid
    and the lumped run lengths. Alarm statistic (as ``competitors.BOCPDVariance``): the log-odds that the
    current run is shorter than w; direction: +1 if the posterior mean variance of the runs of 1 to
    ``w - 1`` observations is at least that of the runs of at least w.

    ``nu`` fixed degrees of freedom (``nu > 2``), or ``None`` for a Gaussian likelihood. Prior: inverse
    gamma on V with shape ``alpha0`` and mean ``prior_var`` (per path), discretised on the grid as
    ``pi_j proportional to V_j^(-alpha0) exp(-beta0 / V_j)``, ``beta0 = (alpha0 - 1) prior_var`` (the
    density times V, since the grid is uniform in log V).
    """

    def __init__(self, reps: int, hazard: float, log_var_grid: np.ndarray, short: int, nu: float | None = 3.0,
                 alpha0: float = 2.0, prior_var=1.0):
        if not 0 < hazard < 1:
            raise ValueError("hazard must lie in (0, 1)")
        if nu is not None and nu <= 2:
            raise ValueError("nu must exceed 2 (finite variance)")
        if short < 2:
            raise ValueError("short must be at least 2")
        self.reps, self.w = int(reps), int(short)
        self.H = float(hazard)
        self.logV = np.asarray(log_var_grid, float)
        self.V = np.exp(self.logV)
        self.G = self.V.size
        self.nu = nu
        self.alpha0 = float(alpha0)
        if nu is None:
            self._c0 = -0.5 * np.log(2 * np.pi * self.V)
        else:
            s2 = self.V * (nu - 2.0) / nu
            self._s2nu = nu * s2
            self._c0 = gammaln((nu + 1) / 2) - gammaln(nu / 2) - 0.5 * np.log(np.pi * nu * s2)
        self.set_prior(prior_var)
        self.J = np.zeros((self.reps, self.w + 1, self.G))
        self.J[:, 0, :] = self.prior  # no data yet: a new run

    def set_prior(self, prior_var):
        """Prior of FUTURE new runs (mean of V); may change at every step (a predictable prior)."""
        pv = np.broadcast_to(np.asarray(prior_var, float), (self.reps,))
        beta0 = (self.alpha0 - 1.0) * pv
        lp = -self.alpha0 * self.logV[None, :] - beta0[:, None] / self.V[None, :]
        lp -= lp.max(1, keepdims=True)
        p = np.exp(lp)
        self.prior = p / p.sum(1, keepdims=True)

    def loglik(self, x: np.ndarray) -> np.ndarray:
        """``log f(x | V_j)``, shape ``(reps, G)``."""
        x2 = (np.asarray(x, float) ** 2)[:, None]
        if self.nu is None:
            return self._c0[None, :] - 0.5 * x2 / self.V[None, :]
        return self._c0[None, :] - 0.5 * (self.nu + 1.0) * np.log1p(x2 / self._s2nu[None, :])

    def warm_start(self, xcal: np.ndarray):
        """All mass on the lumped bucket (a run of at least w steps) with the posterior given ``xcal``.

        Encodes "the calibration window has no change", as ``competitors.BOCPDVariance.warm_start``.
        ``xcal`` is ``(reps, n_cal)``; missing values are skipped.
        """
        xcal = np.atleast_2d(np.asarray(xcal, float))
        ll = np.zeros((self.reps, self.G))
        for t in range(xcal.shape[1]):
            x = xcal[:, t]
            ok = np.isfinite(x)
            ll += np.where(ok[:, None], self.loglik(np.where(ok, x, 0.0)), 0.0)
        lp = np.log(np.maximum(self.prior, 1e-300)) + ll
        lp -= lp.max(1, keepdims=True)
        p = np.exp(lp)
        self.J[:] = 0.0
        self.J[:, self.w, :] = p / p.sum(1, keepdims=True)

    def update(self, x):
        x = np.broadcast_to(np.asarray(x, float), (self.reps,))
        ok = np.isfinite(x)
        keep = None if ok.all() else (~ok, self.J[~ok].copy())  # rows without an observation keep their state
        lf = self.loglik(np.where(ok, x, 0.0))
        f = np.exp(lf - lf.max(1, keepdims=True))
        w, H = self.w, self.H
        J = self.J
        J *= f[:, None, :]  # joint of (run, V) and x_t, up to a factor per row
        tot = J.sum(axis=(1, 2))
        last = J[:, w - 1, :] + J[:, w, :]
        J[:, 1:w, :] = J[:, : w - 1, :].copy()  # run r becomes r + 1
        J[:, w, :] = last
        J[:, 0, :] = (H / (1.0 - H)) * tot[:, None] * self.prior  # a new run starts from the prior
        J *= ((1.0 - H) / tot)[:, None, None]
        if keep is not None:
            J[keep[0]] = keep[1]

    def short_log_odds(self) -> np.ndarray:
        s = self.J[:, : self.w, :].sum(axis=(1, 2))
        lg = self.J[:, self.w, :].sum(1)
        with np.errstate(divide="ignore"):
            return np.log(s) - np.log(lg)

    def direction(self) -> np.ndarray:
        ps = self.J[:, 1 : self.w, :].sum(1)  # (reps, G): runs of 1 .. w-1 observations
        pl = self.J[:, self.w, :]
        with np.errstate(invalid="ignore", divide="ignore"):
            vs = (ps * self.V).sum(1) / ps.sum(1)
            vl = (pl * self.V).sum(1) / pl.sum(1)
        return np.where(vs >= vl, 1.0, -1.0)


def grid_bocpd_paths(r, hazard: float, log_var_grid, short: int = 250, nu: float | None = 3.0, cal=None,
                     alpha0: float = 2.0, prior_var=None, prior_var_path=None, chunk: int = 250):
    """Run ``GridBOCPD`` over rows of ``r``; returns ``(L, S)``: short-run log-odds and direction, ``(reps, n)``.

    Same conventions as ``competitors.bocpd_paths``: ``cal`` (optional) is a change-free calibration
    window that sets the prior mean (its mean square unless ``prior_var`` is given) and warm-starts the
    filter; ``prior_var_path`` (``(reps, n)``, predictable) replaces the prior mean of new runs step by
    step where it is finite and positive. Rows are processed in chunks of ``chunk`` paths.
    """
    r = np.atleast_2d(np.asarray(r, float))
    reps, n = r.shape
    L = np.empty((reps, n))
    S = np.empty((reps, n))
    for i0 in range(0, reps, chunk):
        sl = slice(i0, min(reps, i0 + chunk))
        p = sl.stop - sl.start
        pv = None
        cc = None
        if cal is not None:
            cc = np.atleast_2d(np.asarray(cal, float))[sl]
            pv = np.nanmean(cc**2, axis=1)
        if prior_var is not None:
            pv = np.broadcast_to(np.asarray(prior_var, float), (reps,))[sl]
        if pv is None:
            pv = np.ones(p)
        pv = np.maximum(pv, 1e-300)
        f = GridBOCPD(p, hazard, log_var_grid, short, nu, alpha0, pv)
        if cc is not None:
            f.warm_start(cc)
        for t in range(n):
            if prior_var_path is not None:
                q = prior_var_path[sl, t]
                pv = np.where(np.isfinite(q) & (q > 0), q, pv)
                f.set_prior(pv)
            f.update(r[sl, t])
            L[sl, t] = f.short_log_odds()
            S[sl, t] = f.direction()
    return L, S


# ----------------------------------------------------------------------------
# Streaming series with recalibration after each alarm (real data)
# ----------------------------------------------------------------------------


def _calibration(buf: np.ndarray, c: float):
    """Plug-in calibration of ``monitor.VolatilityMonitor`` on a window: (scale, tau2, m0) of its finite values."""
    r = buf[np.isfinite(buf)]
    if r.size == 0:
        return 0.0, 0.0, 0.0, r
    scale = float(np.median(np.abs(r)) * MAD_NORMAL_SCALE)
    tau2 = (c * scale) ** 2
    return scale, tau2, float(np.mean(np.minimum(r**2, tau2))), r


def _usable(*vals):
    return all(math.isfinite(v) and v > 0 for v in vals)


def calibrate_window(u: np.ndarray, s: int, cal_len: int, c: float = 3.0, cal_max: int | None = None):
    """Calibration window starting at ``s`` with the extension and fallback rules of ``VolatilityMonitor``.

    Returns ``(end, scale, tau2, m0, r_finite)``, with monitoring starting at ``end``, or ``None`` if the
    series ends first. A window that is degenerate (scale or m0 not positive) is extended one return at a
    time until it is usable and holds at least ``cal_len // 2`` finite returns, up to ``cal_max``
    (default ``2 cal_len``) returns; at the cap the scale falls back to the scaled median of the nonzero
    absolute returns; if that fails too, calibration restarts after the window.
    """
    cal_max = 2 * cal_len if cal_max is None else cal_max
    n = u.size
    while True:
        e = s + cal_len
        if e > n:
            return None
        extended = False
        while True:
            scale, tau2, m0, r = _calibration(u[s:e], c)
            usable = _usable(scale, tau2, m0)
            if usable and (not extended or r.size >= cal_len // 2):
                return e, scale, tau2, m0, r
            if e - s < cal_max:
                extended = True
                e += 1
                if e > n:
                    return None
                continue
            if not usable:
                nz = np.abs(r[r != 0])
                scale = float(np.median(nz) * MAD_NORMAL_SCALE) if nz.size else 0.0
                tau2 = (c * scale) ** 2
                m0 = float(np.mean(np.minimum(r**2, tau2))) if r.size else 0.0
                if not _usable(scale, tau2, m0):
                    s = e
                    break
            return e, scale, tau2, m0, r


class _SegmentRule:
    """Side statistics of one monitored segment, computed block by block with the state carried over."""

    def __init__(self, kind: str, r_cal: np.ndarray, tau2: float, m0: float, params: dict):
        self.kind = kind
        self.p = params
        self.tau2, self.m0 = tau2, m0
        band = params.get("band", 0.0)
        if kind == "sr":
            self.m_up, self.m_dn = m0 * (1 + band), m0 * (1 - band)
            B = tau2 / self.m_dn if self.m_dn > 0 else math.inf
            self.lam_u = np.asarray(DEFAULT_UP_GRID, float)
            fd = np.asarray(DEFAULT_DOWN_FRACTIONS, float)
            self.down_active = bool(math.isfinite(B) and B > 1)
            self.lam_d = fd / (B - 1.0) if self.down_active else np.zeros(fd.size)
            self.logS_u = np.full(self.lam_u.size, -np.inf)
            self.logS_d = np.full(self.lam_d.size, -np.inf)
        elif kind == "cusum":
            self.Gu = self.Gd = 0.0
        elif kind == "hkz":
            z = r_cal**2 if params.get("on", "x") == "r2" else np.minimum(r_cal**2, tau2)
            self.m = float(z.size)
            self.zbar = float(z.mean())
            self.s = float(np.sqrt(bartlett_lrv(z[None, :])[0]) * np.sqrt(self.m))
            self.k = 0.0
            self.Q = 0.0
        else:
            raise ValueError(kind)

    def block(self, r: np.ndarray):
        """Statistics for the finite returns ``r`` of the next block: ``(up, dn)``, each ``(len(r),)``."""
        x = np.minimum(r**2, self.tau2)
        if self.kind == "sr":
            out = []
            for lam, logS0, y, sign in ((self.lam_u, self.logS_u, x / self.m_up, 1.0),
                                        (self.lam_d, self.logS_d, x / self.m_dn if self.m_dn > 0 else x, -1.0)):
                with np.errstate(divide="ignore"):
                    lf = np.log1p(lam[None, :] * sign * (y[:, None] - 1.0))
                Lc = np.cumsum(lf, axis=0)
                # S_t = exp(Lc_t) (S_0 + sum_{j<=t} exp(-Lc_{j-1})), Lc_0 = 0
                a = np.concatenate([np.zeros((1, lam.size)), -Lc[:-1]], axis=0)
                acc = np.logaddexp(logS0[None, :], np.logaddexp.accumulate(a, axis=0))
                logS = Lc + acc
                m = logS.max(1, keepdims=True)
                out.append((m[:, 0] + np.log(np.exp(logS - m).sum(1)) - math.log(lam.size), logS[-1]))
            (su, self.logS_u), (sd, self.logS_d) = out
            if not self.down_active:
                sd = np.full_like(sd, -np.inf)
            return su, sd
        if self.kind == "cusum":
            y = x / self.m0
            res = []
            for d, G0 in ((y - 1.0 - self.p["k_up"], self.Gu), (1.0 - y - self.p["k_dn"], self.Gd)):
                S = G0 + np.cumsum(d)
                G = S - np.minimum(0.0, np.minimum.accumulate(S))  # max(0, ...) recursion started at G0
                res.append(G)
            self.Gu, self.Gd = float(res[0][-1]), float(res[1][-1])
            return res[0], res[1]
        z = r**2 if self.p.get("on", "x") == "r2" else x
        k = self.k + np.arange(1, z.size + 1)
        Q = self.Q + np.cumsum(z) - (k - self.k) * self.zbar
        self.k, self.Q = float(k[-1]), float(Q[-1])
        m, gam = self.m, self.p.get("gamma", 0.25)
        if self.p.get("boundary", "hhks") == "hhks":
            D = Q / (self.s * (1.0 + k / m) * (k / (m + k)) ** gam)
        else:
            xx = (m + k) / m
            D = np.sign(Q) * np.sqrt(np.maximum((Q / self.s) ** 2 / (xx * (xx - 1.0)) - np.log(xx / (xx - 1.0)), 0.0))
        return D, -D


def restart_events(u: np.ndarray, cal_len: int, kind: str, h_up: float, h_dn: float, params: dict | None = None,
                   c: float = 3.0, block: int = 2880):
    """Alarm times and directions of a calibrate, monitor, alarm, recalibrate loop on one series.

    As ``monitor.VolatilityMonitor``: calibrate on the next ``cal_len`` entries (``calibrate_window``),
    monitor the finite returns after the window, and after an alarm at t start a new calibration at t + 1.
    ``kind``: ``"sr"`` (the e-detector's statistic; ``params["band"]`` widens its null; thresholds in log
    units, so ``h = log(2 A)`` reproduces the e-detector), ``"cusum"`` (``params`` holds ``k_up`` and
    ``k_dn``) or ``"hkz"`` (``params``: ``on`` = ``"x"`` or ``"r2"``, ``gamma``, ``boundary``). Up wins a tie.
    Returns a list of ``(t, "up" | "down")``.
    """
    params = {} if params is None else params
    u = np.asarray(u, float)
    n = u.size
    events = []
    s = 0
    while True:
        cw = calibrate_window(u, s, cal_len, c)
        if cw is None:
            return events
        e, _, tau2, m0, rfin = cw
        rule = _SegmentRule(kind, rfin, tau2, m0, params)
        t0 = e
        alarm = None
        while t0 < n and alarm is None:
            t1 = min(n, t0 + block)
            seg = u[t0:t1]
            idx = np.flatnonzero(np.isfinite(seg))
            if idx.size:
                su, sd = rule.block(seg[idx])
                cu, cd = su >= h_up, sd >= h_dn
                hit = cu | cd
                if hit.any():
                    j = int(hit.argmax())
                    alarm = (t0 + int(idx[j]), "up" if cu[j] else "down")
            t0 = t1
        if alarm is None:
            return events
        events.append(alarm)
        s = alarm[0] + 1


def match_restart_counts(u, cal_len, kind, params, n_up, n_dn, lo_hi, iters: int = 30, rounds: int = 3, c: float = 3.0,
                         refine: int = 0, refine_points: int = 41, refine_span: float = 0.5,
                         grid_points: int = 0, grid_span: float = 2.0, grid_if_error_above: int = 2):
    """Thresholds ``(h_up, h_dn)`` whose numbers of up and down alarms are as close as possible to the targets.

    Coordinate-wise bisection as in ``competitors.match_alarm_counts``, but every evaluation reruns the
    restart loop. The recalibration points depend on the thresholds, so the counts are not monotone in
    either threshold (a small change of one alarm time moves every later calibration window), and the
    bisection can stop several alarms away from the targets. The error of a pair is
    |up count - n_up| + |down count - n_dn|.

    With ``refine > 0``, up to ``refine`` rounds of a local search follow the bisection: each threshold in
    turn is scanned over ``refine_points`` values h exp(s), s in [-refine_span, refine_span], with the other
    held fixed, keeping the value with the smallest error (ties: the smallest |s|); the search stops at an
    exact match or when a round brings no improvement. If the error is still above ``grid_if_error_above``
    and ``grid_points > 0``, a second start is tried: the pair with the smallest error (ties: the closest)
    on the ``grid_points`` x ``grid_points`` grid (h_up exp(s), h_dn exp(t)), s and t in
    [-grid_span, grid_span], around the bisection result, followed by the same local search; the better of
    the two results is kept. Returns ``(h_up, h_dn, events)``.
    """
    lo0, hi0 = lo_hi
    h_up, h_dn = hi0, hi0

    def counts(hu, hd):
        ev = restart_events(u, cal_len, kind, hu, hd, params, c)
        k = sum(1 for e in ev if e[1] == "up")
        return k, len(ev) - k

    def error(hu, hd):
        ku, kd = counts(hu, hd)
        return abs(ku - n_up) + abs(kd - n_dn)

    for _ in range(rounds):
        lo, hi = lo0, hi0
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if counts(mid, h_dn)[0] > n_up else (lo, mid)
        h_up = hi
        lo, hi = lo0, hi0
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if counts(h_up, mid)[1] > n_dn else (lo, mid)
        h_dn = hi

    steps = np.linspace(-refine_span, refine_span, refine_points)
    steps = steps[np.argsort(np.abs(steps), kind="stable")]  # ties go to the smallest change

    def local(hu0, hd0, err):
        for _ in range(refine):
            if err == 0:
                break
            before = err
            for side in (0, 1):
                best = (err, hu0, hd0)
                for st in steps:
                    hu, hd = (hu0 * np.exp(st), hd0) if side == 0 else (hu0, hd0 * np.exp(st))
                    e = error(hu, hd)
                    if e < best[0]:
                        best = (e, hu, hd)
                        if e == 0:
                            break
                err, hu0, hd0 = best
                if err == 0:
                    break
            if err >= before:
                break
        return err, hu0, hd0

    err, hu_a, hd_a = local(h_up, h_dn, error(h_up, h_dn))
    if grid_points > 0 and err > grid_if_error_above:
        g = np.linspace(-grid_span, grid_span, grid_points)
        best = (np.inf, 0.0, h_up, h_dn)
        for s_ in g:
            for t_ in g:
                e = error(h_up * np.exp(s_), h_dn * np.exp(t_))
                if (e, abs(s_) + abs(t_)) < best[:2]:
                    best = (e, abs(s_) + abs(t_), h_up * np.exp(s_), h_dn * np.exp(t_))
        err_b, hu_b, hd_b = local(best[2], best[3], best[0])
        if err_b < err:
            err, hu_a, hd_a = err_b, hu_b, hd_b
    return hu_a, hd_a, restart_events(u, cal_len, kind, hu_a, hd_a, params, c)
