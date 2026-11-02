"""Anytime-valid monitor for changes in the dependence between two assets (alternative A7).

The main package monitors the volatility of ONE asset through the bounded
statistic ``min(r^2, tau^2)``. This module applies the same recipe (robust
scale, truncation, betting e-detector) to a PAIR of assets and monitors a
bounded statistic whose conditional mean moves with their correlation.

Two statistics are provided. Both are bounded, which is what the e-detector
needs for validity in both directions (``docs/math.md``, section 5.3).

* ``sign``: ``S_t = sign(r1_t) * sign(r2_t)`` in ``{-1, 0, 1}``. Its mean is
  ``P(same sign) - P(opposite signs)``, Blomqvist's beta, a "median
  correlation". For Gaussian (and any elliptical, for example Student t)
  returns with correlation ``rho`` its mean is ``(2 / pi) arcsin(rho)``, which
  increases with ``rho``. It does not change when either asset is multiplied by
  any positive, predictable factor, so a volatility change alone cannot move
  it. A jump changes one ``S_t`` by at most 2.
* ``product``: ``P_t = clip(r1 / s1, -c, c) * clip(r2 / s2, -c, c)`` in
  ``[-c^2, c^2]``, the product of standardised returns truncated at ``c`` robust
  standard deviations (Bapat's rule applied to each asset). Its mean is close
  to the Pearson correlation for Gaussian data. It is robust to jumps, but it
  DOES move when the volatility of either asset changes, so it mistakes a
  volatility shift for a correlation shift. It is kept as the natural
  "truncated covariance" candidate and as a comparison.

Validity (proof sketch; the full statement is in README.md)
-----------------------------------------------------------
Let ``Y_t`` be the chosen statistic, with ``lo <= Y_t <= hi`` almost surely, and
let ``m0`` be the null value of its conditional mean. For an increase the null is
``E[Y_t | F_{t-1}] <= m0`` for all t. With a predictable bet
``lam_t in [0, 1 / (m0 - lo))`` the factor ``f_t = 1 + lam_t (Y_t - m0)`` is
strictly positive and has conditional mean at most 1 under the null, so the
products of factors started at every time j are test supermartingales and the
Shiryaev-Roberts statistic ``M_t = f_t (M_{t-1} + 1)`` is an e-detector:
``E[alarm time] >= A`` and ``P(alarm within m steps) <= m / A`` (Theorems A and
B of ``docs/math.md``, whose proofs only use ``f_t >= 0`` and
``E[f_t | F_{t-1}] <= 1``). For a decrease use ``f_t = 1 + lam_t (m0 - Y_t)`` with
``lam_t in [0, 1 / (hi - m0))``. The null is a statement about one conditional
moment of a bounded statistic: no distribution, independence or moment
assumption is made about the returns themselves.

Benchmarks: a Gaussian likelihood ratio e-detector for the correlation (valid
only under exact bivariate normality with the calibration scales) and a naive
Page CUSUM on the raw, untruncated product of standardised returns, whose
threshold must be tuned by simulation under an assumed model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from avmonitor.detectors import DetectorResult
from avmonitor.robust import MAD_NORMAL_SCALE

STATISTICS = ("sign", "product")
DEFAULT_FRACTIONS = (0.05, 0.1, 0.2, 0.35, 0.5, 0.75)  # bets as fractions of the largest admissible bet
DEFAULT_LR_DELTAS = (0.15, 0.3, 0.45)  # correlation shifts in the Gaussian LR mixture


# ----------------------------------------------------------------------------
# Statistics and their Gaussian means
# ----------------------------------------------------------------------------


def robust_scale(r: np.ndarray) -> float:
    """Scaled ``median(|r|)`` over the finite entries (a robust per-step standard deviation)."""
    r = np.asarray(r, float)
    r = r[np.isfinite(r)]
    return float(np.median(np.abs(r)) * MAD_NORMAL_SCALE) if r.size else 0.0


def sign_concordance(r1, r2) -> np.ndarray:
    """``sign(r1) * sign(r2)``: +1 same sign, -1 opposite signs, 0 if either return is zero, NaN if missing."""
    r1 = np.asarray(r1, float)
    r2 = np.asarray(r2, float)
    out = np.sign(r1) * np.sign(r2)
    out[~(np.isfinite(r1) & np.isfinite(r2))] = np.nan
    return out


def truncated_product(r1, r2, scale1: float, scale2: float, c: float = 3.0) -> np.ndarray:
    """Product of standardised returns, each clipped at ``c`` robust standard deviations."""
    z1 = np.clip(np.asarray(r1, float) / scale1, -c, c)
    z2 = np.clip(np.asarray(r2, float) / scale2, -c, c)
    return z1 * z2


def gaussian_sign_mean(rho) -> np.ndarray:
    """``E[sign(Z1) sign(Z2)] = (2 / pi) arcsin(rho)`` for bivariate normal (and elliptical) pairs."""
    return 2.0 / np.pi * np.arcsin(np.asarray(rho, float))


def rho_from_sign_mean(m) -> np.ndarray:
    """Inverse of ``gaussian_sign_mean``: the correlation implied by a sign-concordance mean."""
    return np.sin(np.pi * np.clip(np.asarray(m, float), -1, 1) / 2.0)


def gaussian_product_moments(rho: float, c: float = 3.0, n: int = 400_000, seed: int = 0) -> tuple[float, float]:
    """Mean and variance of the truncated product for a bivariate normal pair (Monte Carlo, fixed seed)."""
    rng = np.random.default_rng(seed)
    z1 = rng.standard_normal(n)
    z2 = rho * z1 + math.sqrt(1 - rho**2) * rng.standard_normal(n)
    p = np.clip(z1, -c, c) * np.clip(z2, -c, c)
    return float(p.mean()), float(p.var())


# ----------------------------------------------------------------------------
# Calibration
# ----------------------------------------------------------------------------


@dataclass
class PairBaseline:
    stat: str  # "sign" or "product"
    scale1: float  # robust sd of asset 1 in the calibration window
    scale2: float
    c: float  # truncation level in robust sds (product only)
    lo: float  # almost sure lower bound of the statistic
    hi: float  # almost sure upper bound
    m0: float  # calibration mean of the statistic
    m_up: float  # null value used by the increase side (m0 + band + margin)
    m_dn: float  # null value used by the decrease side
    var: float  # calibration variance of the statistic
    n_cal: int
    rho_hat: float  # correlation implied by m0 (sign: sin(pi m0 / 2); product: m0)

    @property
    def lam_max_up(self) -> float:
        return 1.0 / (self.m_up - self.lo) if self.m_up > self.lo else 0.0

    @property
    def lam_max_dn(self) -> float:
        return 1.0 / (self.hi - self.m_dn) if self.hi > self.m_dn else 0.0


def pair_statistic(r1, r2, base: PairBaseline) -> np.ndarray:
    if base.stat == "sign":
        return sign_concordance(r1, r2)
    if base.stat == "product":
        return truncated_product(r1, r2, base.scale1, base.scale2, base.c)
    raise ValueError("stat must be 'sign' or 'product'")


def calibrate_pair(
    r1_cal, r2_cal, stat: str = "sign", c: float = 3.0, band: float = 0.0, margin_z: float = 0.0
) -> PairBaseline:
    """Null values from a calibration window of paired returns.

    ``band`` widens the null additively: the increase side bets against
    ``E[Y] <= m0 + band`` and the decrease side against ``E[Y] >= m0 - band``
    (the statistic lives on a fixed scale, so an additive band is the natural
    one; 0.1 on the sign scale is about 0.15 in correlation near zero).
    ``margin_z`` adds that many standard errors of ``m0`` on top, in the
    conservative direction, to allow for estimation error in the baseline.
    """
    if stat not in STATISTICS:
        raise ValueError("stat must be 'sign' or 'product'")
    r1 = np.asarray(r1_cal, float)
    r2 = np.asarray(r2_cal, float)
    s1, s2 = robust_scale(r1), robust_scale(r2)
    if stat == "sign":
        y = sign_concordance(r1, r2)
        lo, hi = -1.0, 1.0
    else:
        y = truncated_product(r1, r2, s1, s2, c)
        lo, hi = -(c**2), c**2
    y = y[np.isfinite(y)]
    n = int(y.size)
    m0 = float(y.mean()) if n else 0.0
    v = float(y.var(ddof=1)) if n > 1 else float("nan")
    se = math.sqrt(v / n) if n > 1 and np.isfinite(v) else 0.0
    m_up = min(m0 + band + margin_z * se, hi)
    m_dn = max(m0 - band - margin_z * se, lo)
    rho_hat = float(rho_from_sign_mean(m0)) if stat == "sign" else m0
    return PairBaseline(stat, s1, s2, c, lo, hi, m0, m_up, m_dn, v, n, rho_hat)


# ----------------------------------------------------------------------------
# Betting e-detector for a bounded statistic with an additive null
# ----------------------------------------------------------------------------


def bounded_edetector(
    y: np.ndarray,
    m0,
    lo: float,
    hi: float,
    direction: str = "up",
    threshold: float = 1000.0,
    fractions=DEFAULT_FRACTIONS,
    kind: str = "SR",
    start: int = 0,
    return_path: bool = False,
) -> DetectorResult:
    """Shiryaev-Roberts (or CUSUM) betting e-detector on rows of a bounded statistic.

    ``y`` has shape (reps, n) or (n,), with ``lo <= y <= hi``. ``m0`` is the null
    value of the conditional mean (a scalar or one value per row). For
    ``direction="up"`` the null is ``E[y_t | past] <= m0`` and the bets are
    ``fractions * 1 / (m0 - lo)``; for ``"down"`` the null is ``E[y_t | past] >= m0``
    and the bets are ``fractions * 1 / (hi - m0)``. The statistic is the equal
    weight average over the bets (itself an e-detector). NaN entries are skipped
    (factor 1, no step counted). The alarm index is the first t with
    ``M_t >= threshold`` (-1 if none).
    """
    y = np.atleast_2d(np.asarray(y, float))[:, start:]
    reps, n = y.shape
    m0 = np.broadcast_to(np.asarray(m0, float), (reps,)).astype(float)
    fr = np.asarray(fractions, float)
    if np.any(fr < 0) or np.any(fr >= 1):
        raise ValueError("fractions must lie in [0, 1)")
    if direction == "up":
        lam_max = np.where(m0 > lo, 1.0 / np.maximum(m0 - lo, 1e-300), 0.0)
        sign = 1.0
    elif direction == "down":
        lam_max = np.where(hi > m0, 1.0 / np.maximum(hi - m0, 1e-300), 0.0)
        sign = -1.0
    else:
        raise ValueError("direction must be 'up' or 'down'")
    lam = lam_max[:, None] * fr[None, :]  # (reps, K)
    K = fr.size
    S = np.zeros((reps, K))
    alarm = np.full(reps, -1)
    active = np.ones(reps, bool)
    path = np.empty(n) if return_path else None
    for t in range(n):
        g = sign * (y[:, t] - m0)
        ok = np.isfinite(g)
        f = np.where(ok[:, None], 1.0 + lam * np.where(ok, g, 0.0)[:, None], 1.0)
        with np.errstate(over="ignore"):
            if kind == "SR":
                S = np.where(ok[:, None], f * (S + 1.0), S)
            elif kind == "CUSUM":
                S = np.where(ok[:, None], f * np.maximum(S, 1.0), S)
            else:
                raise ValueError("kind must be 'SR' or 'CUSUM'")
        M = S.mean(axis=1)
        if path is not None:
            path[t] = M[0]
        hit = active & (M >= threshold)
        alarm[hit] = t + start
        active &= ~hit
        if not active.any() and path is None:
            break
    return DetectorResult(alarm, path, None, {"lam_max": lam_max, "fractions": fr})


def two_sided_alarm(up: DetectorResult, down: DetectorResult) -> np.ndarray:
    """First alarm of either side (run each side at ``2 A`` to keep the overall guarantee at ``A``)."""
    a, b = up.alarm.copy(), down.alarm.copy()
    big = np.iinfo(a.dtype).max
    a[a < 0] = big
    b[b < 0] = big
    out = np.minimum(a, b)
    out[out == big] = -1
    return out


# ----------------------------------------------------------------------------
# Benchmarks
# ----------------------------------------------------------------------------


def gaussian_corr_lr_factors(z1, z2, rho0, rho1) -> np.ndarray:
    """Likelihood ratio of a bivariate normal pair with unit variances: correlation ``rho1`` against ``rho0``."""
    z1 = np.asarray(z1, float)
    z2 = np.asarray(z2, float)
    rho0 = np.asarray(rho0, float)
    rho1 = np.asarray(rho1, float)

    def q(rho):
        return (z1**2 - 2 * rho * z1 * z2 + z2**2) / (1 - rho**2)

    logf = -0.5 * np.log((1 - rho1**2) / (1 - rho0**2)) - 0.5 * q(rho1) + 0.5 * q(rho0)
    with np.errstate(over="ignore"):
        return np.exp(np.minimum(logf, 700.0))


def gaussian_corr_lr_edetector(
    z1: np.ndarray,
    z2: np.ndarray,
    rho0,
    direction: str = "up",
    deltas=DEFAULT_LR_DELTAS,
    threshold: float = 1000.0,
    rho_cap: float = 0.95,
    start: int = 0,
) -> DetectorResult:
    """Shiryaev-Roberts mixture of Gaussian likelihood ratios for ``rho = rho0 +/- delta``.

    ``z1, z2`` are returns divided by their calibration scales. This is an
    e-detector only if the pre-change pairs are exactly bivariate normal with
    unit variances and correlation ``rho0``. Jumps, fat tails and a change in
    the volatility of either asset all break that assumption.
    """
    z1 = np.atleast_2d(np.asarray(z1, float))[:, start:]
    z2 = np.atleast_2d(np.asarray(z2, float))[:, start:]
    reps, n = z1.shape
    rho0 = np.broadcast_to(np.asarray(rho0, float), (reps,)).astype(float)
    rho0 = np.clip(rho0, -rho_cap, rho_cap)
    sgn = 1.0 if direction == "up" else -1.0
    if direction not in ("up", "down"):
        raise ValueError("direction must be 'up' or 'down'")
    rho1 = np.clip(rho0[:, None] + sgn * np.asarray(deltas, float)[None, :], -rho_cap, rho_cap)  # (reps, K)
    S = np.zeros_like(rho1)
    alarm = np.full(reps, -1)
    active = np.ones(reps, bool)
    for t in range(n):
        a, b = z1[:, t : t + 1], z2[:, t : t + 1]
        ok = np.isfinite(a[:, 0]) & np.isfinite(b[:, 0])
        f = gaussian_corr_lr_factors(np.where(ok[:, None], a, 0.0), np.where(ok[:, None], b, 0.0), rho0[:, None], rho1)
        with np.errstate(over="ignore"):
            S = np.where(ok[:, None], f * (S + 1.0), S)
        M = S.mean(axis=1)
        hit = active & (M >= threshold)
        alarm[hit] = t + start
        active &= ~hit
        if not active.any():
            break
    return DetectorResult(alarm)


def cusum_additive(y: np.ndarray, m0, drift: float, h: float, direction: str = "up", start: int = 0) -> DetectorResult:
    """Page's CUSUM on an unbounded statistic: ``G_t = max(0, G_{t-1} + s (y_t - m0) - drift)``, alarm at ``G_t >= h``.

    ``s = +1`` for an increase and ``-1`` for a decrease. The threshold ``h`` has
    no finite-sample guarantee; ``calibrate_cusum_threshold`` tunes it under an
    assumed null model, which is the weakness of this naive benchmark.
    """
    y = np.atleast_2d(np.asarray(y, float))[:, start:]
    reps, n = y.shape
    m0 = np.broadcast_to(np.asarray(m0, float), (reps,)).astype(float)
    sgn = 1.0 if direction == "up" else -1.0
    G = np.zeros(reps)
    alarm = np.full(reps, -1)
    active = np.ones(reps, bool)
    for t in range(n):
        inc = sgn * (y[:, t] - m0) - drift
        G = np.maximum(0.0, G + np.where(np.isfinite(inc), inc, 0.0))
        hit = active & (G >= h)
        alarm[hit] = t + start
        active &= ~hit
        if not active.any():
            break
    return DetectorResult(alarm)


def calibrate_cusum_threshold(y_null: np.ndarray, m0: float, drift: float, horizon: int, target: float) -> float:
    """Smallest ``h`` on a grid with ``P(one-sided alarm within horizon) <= target`` on the given null paths."""
    best = None
    for h in np.arange(0.5, 60.0, 0.25):
        rate = float((cusum_additive(y_null[:, :horizon], m0, drift, h, "up").alarm >= 0).mean())
        if rate <= target:
            best = float(h)
            break
    return best if best is not None else 60.0


# ----------------------------------------------------------------------------
# Streaming monitor (calibrate, monitor, alarm, recalibrate)
# ----------------------------------------------------------------------------


@dataclass
class PairEvent:
    index: int  # position in the stream
    direction: str  # "up" (dependence increased) or "down"
    monitored_since: int
    m0: float  # calibration mean of the statistic
    rho_hat: float  # implied pre-change correlation
    stat: float  # value of the side's statistic at the alarm
    scale1: float
    scale2: float


class CorrelationMonitor:
    """Streaming two-sided monitor for the dependence between two return series.

    1. Calibrate on ``cal_len`` finite pairs: robust scales, ``m0`` (mean of the
       statistic), null values ``m0 +/- band``.
    2. Monitor with a Shiryaev-Roberts betting e-detector per side, each at
       level ``2 A``, so the first alarm of either side has average run length
       at least ``A`` and P(false alarm within m steps) <= m / A under the null.
    3. On an alarm, record a ``PairEvent`` and recalibrate on fresh data.

    ``stat="sign"`` (default) is invariant to volatility changes; ``"product"``
    is the truncated-covariance version. ``benchmark="gaussian"`` replaces the
    bets by the Gaussian likelihood ratio mixture (not robust, same protocol).
    """

    def __init__(
        self,
        threshold: float = 5760.0,
        stat: str = "sign",
        c: float = 3.0,
        cal_len: int = 1440,
        null_band: float = 0.0,
        fractions=DEFAULT_FRACTIONS,
        margin_z: float = 0.0,
        benchmark: str | None = None,
        lr_deltas=DEFAULT_LR_DELTAS,
    ):
        if stat not in STATISTICS:
            raise ValueError("stat must be 'sign' or 'product'")
        if benchmark not in (None, "gaussian"):
            raise ValueError("benchmark must be None or 'gaussian'")
        self.A = float(threshold)
        self.stat = stat
        self.c = c
        self.cal_len = int(cal_len)
        self.band = float(null_band)
        self.fr = np.asarray(fractions, float)
        self.margin_z = float(margin_z)
        self.benchmark = benchmark
        self.lr_deltas = np.asarray(lr_deltas, float)
        self.t = -1
        self.events: list[PairEvent] = []
        self.segments: list[dict] = []
        self._reset()

    def _reset(self):
        self.buf1: list[float] = []
        self.buf2: list[float] = []
        self.calibrated = False
        K = len(self.lr_deltas) if self.benchmark == "gaussian" else len(self.fr)
        self.S_up = np.zeros(K)
        self.S_dn = np.zeros(K)

    def _calibrate(self):
        base = calibrate_pair(self.buf1, self.buf2, self.stat, self.c, self.band, self.margin_z)
        if not (base.scale1 > 0 and base.scale2 > 0 and np.isfinite(base.m0)):
            self.buf1, self.buf2 = [], []  # degenerate window: start again
            return
        self.base = base
        self.lam_up = base.lam_max_up * self.fr
        self.lam_dn = base.lam_max_dn * self.fr
        if self.benchmark == "gaussian":
            cap = 0.95
            r0 = float(np.clip(base.rho_hat, -cap, cap))
            self.rho_up0 = float(np.clip(r0 + self.band, -cap, cap))
            self.rho_dn0 = float(np.clip(r0 - self.band, -cap, cap))
            self.rho_up1 = np.clip(self.rho_up0 + self.lr_deltas, -cap, cap)
            self.rho_dn1 = np.clip(self.rho_dn0 - self.lr_deltas, -cap, cap)
        self.start = self.t + 1
        self.calibrated = True
        self.segments.append(dict(start=self.start, cal_len=base.n_cal, m0=base.m0, rho_hat=base.rho_hat,
                                  m_up=base.m_up, m_dn=base.m_dn, scale1=base.scale1, scale2=base.scale2))

    def update(self, r1: float, r2: float) -> PairEvent | None:
        """Feed one pair of returns; returns a ``PairEvent`` if an alarm fires at this step."""
        self.t += 1
        if not (np.isfinite(r1) and np.isfinite(r2)):
            return None
        if not self.calibrated:
            self.buf1.append(float(r1))
            self.buf2.append(float(r2))
            if len(self.buf1) >= self.cal_len:
                self._calibrate()
            return None
        b = self.base
        if self.benchmark == "gaussian":
            z1, z2 = r1 / b.scale1, r2 / b.scale2
            fu = gaussian_corr_lr_factors(z1, z2, self.rho_up0, self.rho_up1)
            fd = gaussian_corr_lr_factors(z1, z2, self.rho_dn0, self.rho_dn1)
        else:
            y = float(pair_statistic(np.array([r1]), np.array([r2]), b)[0])
            fu = 1.0 + self.lam_up * (y - b.m_up)
            fd = 1.0 + self.lam_dn * (b.m_dn - y)
        with np.errstate(over="ignore"):
            self.S_up = fu * (self.S_up + 1.0)
            self.S_dn = fd * (self.S_dn + 1.0)
        mu, md = float(self.S_up.mean()), float(self.S_dn.mean())
        ev = None
        if mu >= 2 * self.A:
            ev = PairEvent(self.t, "up", self.start, b.m0, b.rho_hat, mu, b.scale1, b.scale2)
        elif md >= 2 * self.A:
            ev = PairEvent(self.t, "down", self.start, b.m0, b.rho_hat, md, b.scale1, b.scale2)
        if ev is not None:
            self.events.append(ev)
            self._reset()
        return ev

    def run(self, r1, r2) -> list[PairEvent]:
        r1 = np.asarray(r1, float)
        r2 = np.asarray(r2, float)
        if r1.shape != r2.shape:
            raise ValueError("the two return series must have the same length")
        for a, b in zip(r1, r2):
            self.update(float(a), float(b))
        return self.events


# ----------------------------------------------------------------------------
# k assets: basket concordance
# ----------------------------------------------------------------------------


def basket_concordance(r) -> np.ndarray:
    """Average pairwise sign concordance of ``k`` assets, ``2 / (k (k - 1)) sum_{i<j} sign(r_i) sign(r_j)``.

    ``r`` has the assets on its last axis (shape (n, k) or (reps, n, k)). The
    statistic lies in ``[-1 / (k - 1), 1]`` (lower bound from
    ``sum_{i<j} s_i s_j = ((sum s)^2 - sum s^2) / 2 >= -k / 2``), its mean is the
    average Blomqvist beta over the pairs, and it is invariant to any positive
    predictable rescaling of each asset. NaN where any asset is missing.
    """
    r = np.asarray(r, float)
    k = r.shape[-1]
    if k < 2:
        raise ValueError("need at least two assets")
    s = np.sign(r)
    tot = s.sum(axis=-1)
    out = (tot**2 - (s**2).sum(axis=-1)) / (k * (k - 1))
    out[~np.isfinite(r).all(axis=-1)] = np.nan
    return out


def basket_bounds(k: int) -> tuple[float, float]:
    return -1.0 / (k - 1), 1.0


class BasketMonitor:
    """Streaming two-sided monitor for the average dependence inside a basket of ``k`` assets.

    Same protocol as ``CorrelationMonitor`` (calibrate on ``cal_len`` complete
    rows, each side at ``2 A``, recalibrate after an alarm) with the statistic
    ``basket_concordance``. One e-detector covers the whole basket, so there is
    no Bonferroni-type split over the k (k - 1) / 2 pairs.
    """

    def __init__(self, k: int, threshold: float = 5760.0, cal_len: int = 1440, null_band: float = 0.0,
                 fractions=DEFAULT_FRACTIONS):
        self.k = int(k)
        self.lo, self.hi = basket_bounds(self.k)
        self.A = float(threshold)
        self.cal_len = int(cal_len)
        self.band = float(null_band)
        self.fr = np.asarray(fractions, float)
        self.t = -1
        self.events: list[PairEvent] = []
        self.segments: list[dict] = []
        self._reset()

    def _reset(self):
        self.buf: list[float] = []
        self.calibrated = False
        self.S_up = np.zeros(len(self.fr))
        self.S_dn = np.zeros(len(self.fr))

    def _calibrate(self):
        y = np.asarray(self.buf)
        self.m0 = float(y.mean())
        self.m_up = min(self.m0 + self.band, self.hi)
        self.m_dn = max(self.m0 - self.band, self.lo)
        self.lam_up = self.fr / (self.m_up - self.lo) if self.m_up > self.lo else np.zeros(len(self.fr))
        self.lam_dn = self.fr / (self.hi - self.m_dn) if self.hi > self.m_dn else np.zeros(len(self.fr))
        self.start = self.t + 1
        self.calibrated = True
        self.segments.append(dict(start=self.start, cal_len=len(y), m0=self.m0, m_up=self.m_up, m_dn=self.m_dn))

    def update(self, row) -> PairEvent | None:
        self.t += 1
        y = float(basket_concordance(np.asarray(row, float)[None, :])[0])
        if not np.isfinite(y):
            return None
        if not self.calibrated:
            self.buf.append(y)
            if len(self.buf) >= self.cal_len:
                self._calibrate()
            return None
        with np.errstate(over="ignore"):
            self.S_up = (1.0 + self.lam_up * (y - self.m_up)) * (self.S_up + 1.0)
            self.S_dn = (1.0 + self.lam_dn * (self.m_dn - y)) * (self.S_dn + 1.0)
        mu, md = float(self.S_up.mean()), float(self.S_dn.mean())
        ev = None
        if mu >= 2 * self.A:
            ev = PairEvent(self.t, "up", self.start, self.m0, float("nan"), mu, float("nan"), float("nan"))
        elif md >= 2 * self.A:
            ev = PairEvent(self.t, "down", self.start, self.m0, float("nan"), md, float("nan"), float("nan"))
        if ev is not None:
            self.events.append(ev)
            self._reset()
        return ev

    def run(self, r) -> list[PairEvent]:
        for row in np.asarray(r, float):
            self.update(row)
        return self.events


# ----------------------------------------------------------------------------
# Helpers for real data
# ----------------------------------------------------------------------------


def aggregate_returns(logret, bars: int):
    """Sum consecutive 1-minute log returns into non-overlapping ``bars``-minute returns.

    ``logret`` is a pandas Series indexed by UTC minute. A coarse bar is NaN
    unless all of its minutes are present, so gaps never create artificial
    returns. Index = start time of the coarse bar.
    """
    import pandas as pd

    s = pd.Series(np.asarray(logret, float), index=logret.index)
    g = s.groupby(s.index.floor(f"{bars}min"))
    out = g.sum(min_count=bars)
    cnt = g.count()
    out[cnt < bars] = np.nan
    return out


def rolling_sign_rho(r1, r2, window: int):
    """Correlation implied by the rolling mean of the sign statistic (for plots and diagnostics)."""
    import pandas as pd

    s = pd.Series(sign_concordance(r1, r2), index=getattr(r1, "index", None))
    m = s.rolling(window, min_periods=window // 2).mean()
    return pd.Series(rho_from_sign_mean(m.to_numpy()), index=m.index)
