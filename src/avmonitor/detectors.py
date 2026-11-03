"""Sequential detectors for changes in return variance.

The main detector is a nonparametric e-detector (Shin, Ramdas and Rinaldo 2024)
built from bets on the truncated squared return ``X_t = min(r_t^2, tau^2)``.

Validity in one paragraph
-------------------------
Let ``Y_t = X_t / m0``. Under the null for an increase, ``E[Y_t | past] <= 1``
for every t. A bet ``f_t = 1 + lambda_t (Y_t - 1)`` with predictable
``lambda_t in [0, 1)`` is nonnegative because ``Y_t >= 0`` and has conditional
mean at most 1, so the product of bets started at any time ``j`` is a test
supermartingale ``L^(j)``. The Shiryaev-Roberts style statistic
``M_t = sum_{j <= t} L^(j)_t`` obeys ``M_t = f_t (M_{t-1} + 1)``, and ``M_t - t`` is a
supermartingale. Hence the alarm time ``inf{t : M_t >= A}`` has expected value at
least ``A`` under the null (average run length bound), and for any horizon m the
probability of a false alarm within the first m steps is at most ``m / A``.
The CUSUM style statistic ``C_t = f_t max(C_{t-1}, 1)`` never exceeds ``M_t``, so
both guarantees carry over. Averages over several bets are covered too because
the recursions are linear in the starting wealth. For a decrease the bet is
``1 + lambda (1 - Y_t)``, which is nonnegative only when ``Y_t`` is bounded above.
Truncation provides that bound (``Y_t <= tau^2 / m0``).

Benchmarks: a Gaussian likelihood ratio e-detector (valid only if returns are
exactly normal), Page's CUSUM on squared returns, and a classical window
z-test checked at every step (not anytime-valid).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from .robust import MAD_NORMAL_SCALE, kernel_weights


# ----------------------------------------------------------------------------
# Calibration of the pre-change (null) regime
# ----------------------------------------------------------------------------


@dataclass
class Baseline:
    """Pre-change reference values estimated from a calibration window."""

    scale: float  # robust standard deviation of returns (scaled MAD about zero)
    tau: float  # truncation level for returns
    m0: float  # mean of X = min(r^2, tau^2) in the calibration window
    m0_up: float  # conservative null value for the increase detector
    m0_down: float  # conservative null value for the decrease detector
    var_x: float  # variance of X in the calibration window
    n_cal: int


def calibrate(
    r_cal: np.ndarray, c: float = 3.0, margin_z: float = 0.0, floor: float = 0.0, null_band: float = 0.0,
    *, baseline: str = "plugin", delta: float = 0.01, tau_frac: float = 0.25,
) -> Baseline:
    """Estimate the null regime from pre-change returns.

    ``tau = c * scaled MAD`` (about zero), floored at ``floor``. ``margin_z > 0``
    widens ``m0`` by ``margin_z`` standard errors in the conservative direction,
    which protects validity against estimation error in the baseline (the
    standard error ignores serial dependence, so it is optimistic when
    volatility clusters). ``null_band`` widens the null to "the truncated
    variance has moved by less than this fraction", for example 0.2 means the
    increase detector only bets against ``E[X] <= 1.2 m0``.

    ``baseline="cs"`` (keyword only) returns a baseline for which the false
    alarm guarantee holds with estimation error (``docs/math.md``, section 6):
    ``tau`` from the first ``round(tau_frac * n)`` returns, and ``m0_up`` /
    ``m0_down`` are upper / lower confidence bounds for the mean of
    ``min(r^2, tau^2)`` on the remaining returns, each at level ``delta / 2``
    (then times ``1 +/- null_band``). P(any false alarm within m steps) is then
    at most ``delta + m / A``. In that mode ``m0`` is the plug-in mean of the
    second part, ``var_x`` its variance and ``margin_z`` must be 0. Missing
    values are dropped first. A lower bound of 0 is returned as ``m0_down =
    1e-300``; the decrease side should then be switched off (its bets would be
    near 0 and its Shiryaev-Roberts statistic would just count steps).
    """
    if baseline == "cs":
        from .confseq import cs_baseline  # local import: confseq imports this module's dependencies

        if margin_z:
            raise ValueError("margin_z applies to the plug-in baseline only")
        r = np.asarray(r_cal, float)
        r = r[np.isfinite(r)]
        b = cs_baseline(r[None, :], c=c, delta_upper=delta / 2, delta_lower=delta / 2, tau_frac=tau_frac, floor=floor)
        x2 = np.minimum(r[b.n_tau:] ** 2, b.tau[0] ** 2)
        up = float(b.upper[0]) * (1 + null_band)
        down = max(float(b.lower[0]) * (1 - null_band), 1e-300)
        return Baseline(float(b.scale[0]), float(b.tau[0]), float(b.mean[0]), up, down, float(x2.var(ddof=1)), len(r))
    if baseline != "plugin":
        raise ValueError("baseline must be 'plugin' or 'cs'")
    r_cal = np.asarray(r_cal, float)
    r_cal = r_cal[np.isfinite(r_cal)]  # missing values are dropped (review m17), as in the cs branch
    if r_cal.size < 2:
        raise ValueError("calibration window has fewer than 2 finite returns")
    scale = float(np.median(np.abs(r_cal)) * MAD_NORMAL_SCALE)
    tau = max(c * scale, floor)
    x = np.minimum(r_cal**2, tau**2)
    m0 = float(x.mean())
    if not (tau > 0 and m0 > 0):
        raise ValueError("degenerate calibration window (scale or truncated mean 0, for example more than half "
                         "zero returns); pass floor > 0 or use a longer window")
    v = float(x.var(ddof=1))
    se = np.sqrt(v / len(x))
    up = (m0 + margin_z * se) * (1 + null_band)
    down = max((m0 - margin_z * se) * (1 - null_band), 1e-300)
    return Baseline(scale, tau, m0, up, down, v, len(x))


# ----------------------------------------------------------------------------
# Betting e-detector
# ----------------------------------------------------------------------------

DEFAULT_UP_GRID = (0.02, 0.05, 0.1, 0.2, 0.35, 0.5)
DEFAULT_DOWN_FRACTIONS = (0.1, 0.25, 0.5, 0.75, 0.95)


def check_bets(up_grid, down_fractions) -> None:
    """Raise unless every increase bet and every decrease fraction lies in [0, 1), so that every factor is positive."""
    for name, v in (("increase bets", up_grid), ("decrease fractions", down_fractions)):
        v = np.asarray(v, float)
        if v.size == 0 or not np.all(np.isfinite(v)) or np.any(v < 0) or np.any(v >= 1):
            raise ValueError(f"{name} must be finite and lie in [0, 1)")


@dataclass
class DetectorResult:
    """Output of the array detectors: first alarm index per path (-1 if none) and optional paths."""

    alarm: np.ndarray  # index of the first alarm per path, -1 if none
    stat_path: np.ndarray | None = None  # statistic for path 0 (if requested)
    lam_path: np.ndarray | None = None  # adaptive bet for path 0 (if requested)
    info: dict = field(default_factory=dict)


def _predictable_local_moments(y: np.ndarray, h: int, kernel: str = "gaussian"):
    """Kernel-weighted mean and second moment of ``y`` over the ``h + 1`` steps before t.

    These are the normalised one-sided kernel estimators of Bapat (2026) applied
    to ``Y = X / m0`` (local truncated variance and, through ``E[Y^2]``, local
    truncated quarticity), shifted so that entry t uses data up to ``t - 1``.
    Entries without a full window fall back to (1, 1 + v0) style neutral values.
    """
    w = kernel_weights(h, kernel)
    reps, n = y.shape
    m1 = np.full((reps, n), np.nan)
    m2 = np.full((reps, n), np.nan)
    L = h + 1
    if n > L:
        acc1 = np.zeros((reps, n - L + 1))
        acc2 = np.zeros((reps, n - L + 1))
        for k in range(L):
            seg = y[:, k : k + n - L + 1]
            acc1 += w[k] * seg
            acc2 += w[k] * seg**2
        # acc[:, j] uses y[j .. j+h]; it is predictable for time t = j + h + 1
        m1[:, L:] = acc1[:, : n - L]
        m2[:, L:] = acc2[:, : n - L]
    return m1, m2


def betting_edetector(
    x: np.ndarray,
    m0: float,
    direction: str = "up",
    upper: float | None = None,
    threshold: float = 1000.0,
    kind: str = "SR",
    grid=None,
    adaptive_h: int | None = None,
    adaptive_weight: float = 0.5,
    lam_cap: float = 0.5,
    start: int = 0,
    return_path: bool = False,
) -> DetectorResult:
    """Run the betting e-detector on rows of ``x`` (nonnegative statistics).

    Parameters
    ----------
    x : array (reps, n) or (n,)
        Truncated squared returns ``min(r^2, tau^2)`` (any nonnegative statistic works).
    m0 : float
        Null mean of ``x``. For ``direction="up"`` the null is ``E[x_t | past] <= m0``,
        for ``"down"`` it is ``E[x_t | past] >= m0``.
    upper : float
        Almost sure upper bound of ``x`` (``tau^2``). Required for ``"down"``.
    threshold : float
        Alarm level ``A``; average run length under the null is at least ``A``.
    kind : "SR" or "CUSUM"
    grid : sequence of bets. Defaults to ``DEFAULT_UP_GRID`` for "up" and to
        fractions of the largest admissible bet for "down".
    adaptive_h : int or None
        If set, add a plug-in bettor ``lambda_t = (mu - 1) / (v + (mu - 1)^2)``
        (the approximate growth-rate optimal bet, Waudby-Smith and Ramdas 2024),
        where ``mu`` and ``mu^2 + v`` are kernel-weighted local first and second moments
        of ``Y = x / m0`` with bandwidth ``adaptive_h``. The second moment of a
        truncated squared return is the truncated local quarticity, so this is
        where the Bapat (2026) estimator enters the detector.
    adaptive_weight : float
        Weight of the adaptive bettor in the mixture (the grid shares the rest).
    start : int
        Monitoring starts at this index (earlier entries are ignored).

    Missing values (NaN) in ``x`` are skipped: the statistic is left unchanged at
    that step (no bet and no new chain), as in ``VolatilityMonitor`` (valid when
    missingness does not depend on the missing value, see ``docs/math.md`` 5.6).
    Before November 2026 a single NaN silently disabled the detector (independent
    review, m17), and until review round 9 (M1) a NaN step still added 1 to the
    Shiryaev-Roberts statistic.
    """
    if not 0 <= lam_cap <= 1 or not 0 <= adaptive_weight <= 1:
        raise ValueError("lam_cap and adaptive_weight must lie in [0, 1]")
    x = np.atleast_2d(np.asarray(x, float))[:, start:]
    reps, n = x.shape
    y = x / m0
    if direction == "up":
        lam_max = 1.0
        grid = np.asarray(DEFAULT_UP_GRID if grid is None else grid, float)
        sign = 1.0
    elif direction == "down":
        if upper is None:
            raise ValueError("a decrease detector needs an upper bound on x")
        b = upper / m0
        if b <= 1:
            raise ValueError("upper bound must exceed m0")
        if np.nanmax(x, initial=0.0) > upper * (1 + 1e-12):
            raise ValueError("x exceeds the stated upper bound; a decrease bet factor would turn negative")
        lam_max = 1.0 / (b - 1.0)
        grid = lam_max * np.asarray(DEFAULT_DOWN_FRACTIONS if grid is None else grid, float)
        sign = -1.0
    else:
        raise ValueError("direction must be 'up' or 'down'")
    if np.any(grid < 0) or np.any(grid >= lam_max):
        raise ValueError("bets must lie in [0, lam_max)")

    K = len(grid)
    weights = np.full(K, 1.0 / K)
    lam_adapt = None
    if adaptive_h is not None:
        weights = np.append(weights * (1 - adaptive_weight), adaptive_weight)
        m1, m2 = _predictable_local_moments(y, adaptive_h)
        mu = np.nan_to_num(m1, nan=1.0)
        var = np.nan_to_num(m2 - m1**2, nan=1.0)
        edge = sign * (mu - 1.0)
        lam_adapt = np.clip(edge / (np.maximum(var, 1e-12) + edge**2), 0.0, lam_cap * lam_max)
    S = np.zeros((reps, len(weights)))
    alarm = np.full(reps, -1)
    active = np.ones(reps, bool)
    path = np.empty(n) if return_path else None
    lam_path = np.empty(n) if (return_path and lam_adapt is not None) else None
    for t in range(n):
        g = sign * (y[:, t] - 1.0)  # gain per unit bet
        ok = np.isfinite(g)  # missing value: no bet and no new chain
        g = np.where(ok, g, 0.0)
        f = 1.0 + grid[None, :] * g[:, None]
        if lam_adapt is not None:
            f = np.concatenate([f, (1.0 + lam_adapt[:, t] * g)[:, None]], axis=1)
            if lam_path is not None:
                lam_path[t] = lam_adapt[0, t]
        with np.errstate(over="ignore"):
            if kind == "SR":
                S = np.where(ok[:, None], f * (S + 1.0), S)
            elif kind == "CUSUM":
                S = np.where(ok[:, None], f * np.maximum(S, 1.0), S)
            else:
                raise ValueError("kind must be 'SR' or 'CUSUM'")
        M = S @ weights
        if path is not None:
            path[t] = M[0]
        hit = active & (M >= threshold)
        alarm[hit] = t + start
        active &= ~hit
        if not active.any() and path is None:
            break
    return DetectorResult(alarm, path, lam_path, {"grid": grid, "lam_max": lam_max})


def two_sided_alarm(up: DetectorResult, down: DetectorResult) -> np.ndarray:
    """First alarm of either side. Run each side at threshold ``2 A`` to keep ARL >= A overall."""
    a, b = np.asarray(up.alarm, np.int64), np.asarray(down.alarm, np.int64)
    big = np.iinfo(np.int64).max
    out = np.minimum(np.where(a < 0, big, a), np.where(b < 0, big, b))
    return np.where(out == big, -1, out)


# ----------------------------------------------------------------------------
# Benchmarks
# ----------------------------------------------------------------------------


def gaussian_lr_edetector(
    r: np.ndarray, sigma0: float, ratios=(1.25, 1.5, 2.0), threshold: float = 1000.0, start: int = 0
) -> DetectorResult:
    """Shiryaev-Roberts mixture of Gaussian likelihood ratios for ``sd = ratio * sigma0``.

    This is an e-detector only if pre-change returns are exactly ``N(0, sigma0^2)``.
    Jumps and heavy tails break that assumption.
    """
    z2 = (np.atleast_2d(np.asarray(r, float))[:, start:] / sigma0) ** 2
    ratios = np.asarray(ratios, float)
    reps, n = z2.shape
    S = np.zeros((reps, len(ratios)))
    alarm = np.full(reps, -1)
    active = np.ones(reps, bool)
    logc = -np.log(ratios)
    k = 0.5 * (1 - 1 / ratios**2)
    for t in range(n):
        f = np.exp(np.minimum(logc[None, :] + k[None, :] * z2[:, t : t + 1], 700.0))
        with np.errstate(over="ignore"):
            S = f * (S + 1.0)
        M = S.mean(axis=1)
        hit = active & (M >= threshold)
        alarm[hit] = t + start
        active &= ~hit
        if not active.any():
            break
    return DetectorResult(alarm)


def page_cusum(x: np.ndarray, m0: float, drift: float, h: float, start: int = 0) -> DetectorResult:
    """Page's CUSUM on ``x / m0``: ``G_t = max(0, G_{t-1} + x_t/m0 - 1 - drift)``, alarm at ``G_t >= h``.

    The threshold ``h`` has no finite-sample guarantee and must be calibrated by
    simulation under an assumed null model.
    """
    y = np.atleast_2d(np.asarray(x, float))[:, start:] / m0
    reps, n = y.shape
    G = np.zeros(reps)
    alarm = np.full(reps, -1)
    active = np.ones(reps, bool)
    for t in range(n):
        G = np.maximum(0.0, G + y[:, t] - 1.0 - drift)
        hit = active & (G >= h)
        alarm[hit] = t + start
        active &= ~hit
        if not active.any():
            break
    return DetectorResult(alarm)


def repeated_window_ztest(
    r: np.ndarray, sigma0_sq: float, var_r2: float, window: int = 60, level: float = 0.05, start: int = 0
) -> DetectorResult:
    """Classical test checked at every step: alarm when the window variance z-score is large.

    ``z_t = (sum_{window} r^2 - window * sigma0^2) / sqrt(window * var_r2)``, where
    ``var_r2`` is the variance of a squared return (``2 sigma^4`` for normal data,
    which is where quarticity enters). Each single look has size ``level``, but
    looking at every step makes a false alarm certain eventually.
    """
    r2 = np.atleast_2d(np.asarray(r, float))[:, start:] ** 2
    reps, n = r2.shape
    cs = np.cumsum(r2, axis=1)
    win = cs[:, window - 1 :] - np.concatenate([np.zeros((reps, 1)), cs[:, :-window]], axis=1)
    z = (win - window * sigma0_sq) / np.sqrt(window * var_r2)
    crit = stats.norm.ppf(1 - level / 2)
    hit = np.abs(z) >= crit
    first = np.where(hit.any(axis=1), hit.argmax(axis=1) + window - 1 + start, -1)
    return DetectorResult(first)


# ----------------------------------------------------------------------------
# Run-length summaries
# ----------------------------------------------------------------------------


def summarize_alarms(alarm: np.ndarray, nu: int | None, horizon: int) -> dict:
    """False alarm and delay summaries for alarm indices from paths with change at ``nu``."""
    alarm = np.asarray(alarm)
    out = {}
    if nu is None:
        fa = alarm >= 0
        out["false_alarm_rate"] = float(fa.mean())
        out["mean_run_length_censored"] = float(np.where(fa, alarm + 1, horizon).mean())
        return out
    early = (alarm >= 0) & (alarm < nu)
    detected = alarm >= nu
    out["false_alarm_before_change"] = float(early.mean())
    out["detected"] = float(detected.mean())
    delays = alarm[detected] - nu + 1
    out["median_delay"] = float(np.median(delays)) if delays.size else float("nan")
    out["mean_delay"] = float(delays.mean()) if delays.size else float("nan")
    return out
