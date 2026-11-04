"""Quarticity-stream e-detector: monitor ``min(r^4, tau^4)`` instead of ``min(r^2, tau^2)``.

The main package monitors the truncated squared return ``X_t = min(r_t^2, tau^2)``,
whose mean moves when the variance moves. This prototype feeds the same betting
e-detector (Shin, Ramdas and Rinaldo 2024) with the truncated fourth power

    Q_t = min(r_t^4, tau^4) = X_t^2,

whose mean is the truncated local quarticity of Bapat (2026). It moves with the
variance too, but also with the shape of the return distribution (tail
thickness) and with the volatility of volatility, while the variance stays put.

Null mean and bet range (the two things that must be derived afresh)
--------------------------------------------------------------------
* Null mean. The increase null is ``E[Q_t | past] <= q0`` and the decrease null is
  ``E[Q_t | past] >= q0``. The model-free choice is ``q0 = mean of min(r^4, tau^4)``
  over the calibration window. For normal returns with sd ``sigma`` and
  ``tau = a sigma`` the closed form is ``q0 = 3 kappa_4(a) sigma^4`` with
  ``kappa_4`` from ``avmonitor.robust.truncated_moment_factor``. Note that
  ``q0`` is NOT ``3 m0^2`` (the Gaussian plug-in from the variance baseline):
  at ``a = 3`` that would be 2.7% too high, at ``a = 2.02`` (raw MAD) 21% too
  high, and an inflated ``q0`` makes the decrease detector anti-conservative.
* Bet range. ``Q_t`` lies in ``[0, tau^4]``, so ``Y_t = Q_t / q0`` lies in
  ``[0, B_q]`` with ``B_q = tau^4 / q0``. The increase bet is ``lambda in [0, 1)``
  as before. The decrease bet ``f = 1 + lambda (1 - Y)`` stays positive only for
  ``lambda < 1 / (B_q - 1)``. For normal returns and ``a = 3``, ``B_q = 28.0``
  against ``B_var = tau^2 / m0 = 9.0`` for the variance stream, so the largest
  admissible decrease bet is 0.037 instead of 0.124. This is why decreases
  are slower to detect on the quarticity stream.

Validity: ``Q_t`` is nonnegative and bounded, so every statement in
``docs/math.md`` section 5 (test supermartingale, ``E[run length] >= A``,
``P(false alarm within m steps) <= m / A``) holds unchanged with ``q0`` in
place of ``m0`` and ``B_q`` in place of ``B``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import integrate, stats

from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DetectorResult, betting_edetector
from avmonitor.monitor import Event, VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE, truncated_moment_factor

# Fourth powers are noisier than squares, so the growth-optimal bet for a given
# relative mean shift is smaller. The grid therefore starts lower than
# ``detectors.DEFAULT_UP_GRID``. The experiments use this grid for both streams.
Q_UP_GRID = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5)


# ----------------------------------------------------------------------------
# Streams and baselines
# ----------------------------------------------------------------------------


def variance_stream(r, tau: float) -> np.ndarray:
    """``X_t = min(r_t^2, tau^2)``, the statistic of the main package."""
    r = np.asarray(r, float)
    return np.minimum(r**2, tau**2)


def quarticity_stream(r, tau: float) -> np.ndarray:
    """``Q_t = min(r_t^4, tau^4) = X_t^2``, bounded by ``tau^4``."""
    r = np.asarray(r, float)
    return np.minimum(r**4, tau**4)


def gaussian_null_moments(sigma: float, a: float) -> tuple[float, float]:
    """Closed-form ``(m0, q0)`` for ``N(0, sigma^2)`` returns clipped at ``tau = a sigma``."""
    return sigma**2 * truncated_moment_factor(a, 2), 3 * sigma**4 * truncated_moment_factor(a, 4)


@dataclass
class StreamBaseline:
    scale: float  # scaled MAD about zero (robust sd) of the calibration returns
    tau: float  # truncation level c * scale
    m0: float  # mean of min(r^2, tau^2)
    q0: float  # mean of min(r^4, tau^4)
    b_var: float  # tau^2 / m0, upper bound of the normalised variance stream
    b_q: float  # tau^4 / q0, upper bound of the normalised quarticity stream
    n_cal: int

    @property
    def lam_max_down_var(self) -> float:
        return 1.0 / (self.b_var - 1.0)

    @property
    def lam_max_down_q(self) -> float:
        return 1.0 / (self.b_q - 1.0)


def calibrate_streams(r_cal, c: float = 3.0, floor: float = 0.0) -> StreamBaseline:
    """Baselines of both streams from a pre-change window (model free)."""
    r = np.asarray(r_cal, float)
    r = r[np.isfinite(r)]
    scale = max(float(np.median(np.abs(r)) * MAD_NORMAL_SCALE), floor)
    tau = c * scale
    x = variance_stream(r, tau)
    m0 = float(x.mean())
    q0 = float((x**2).mean())
    return StreamBaseline(scale, tau, m0, q0, tau**2 / m0, tau**4 / q0, r.size)


def run_stream_detector(
    r,
    base: StreamBaseline,
    stream: str,
    direction: str,
    threshold: float,
    grid=None,
    kind: str = "SR",
    start: int = 0,
    return_path: bool = False,
) -> DetectorResult:
    """Betting e-detector on the chosen stream with the matching null mean and bound.

    ``stream`` is ``"variance"`` (``min(r^2, tau^2)``, null mean ``m0``, bound
    ``tau^2``) or ``"quarticity"`` (``min(r^4, tau^4)``, null mean ``q0``, bound
    ``tau^4``). ``grid`` defaults to ``Q_UP_GRID`` for increases and to the
    fractions ``DEFAULT_DOWN_FRACTIONS`` of the largest admissible bet for
    decreases, for both streams.
    """
    if stream == "variance":
        x, m, upper = variance_stream(r, base.tau), base.m0, base.tau**2
    elif stream == "quarticity":
        x, m, upper = quarticity_stream(r, base.tau), base.q0, base.tau**4
    else:
        raise ValueError("stream must be 'variance' or 'quarticity'")
    if direction == "up" and grid is None:
        grid = Q_UP_GRID
    return betting_edetector(
        x, m, direction, upper=upper, threshold=threshold, grid=grid, kind=kind, start=start, return_path=return_path
    )


# ----------------------------------------------------------------------------
# Theory: how far each stream's mean moves under an alternative distribution
# ----------------------------------------------------------------------------


def _clipped_moment(pdf, p: float, a: float) -> float:
    """``E[min(|Z|^p, a^p)]`` for a symmetric density ``pdf`` on the real line."""
    inner = integrate.quad(lambda z: z**p * pdf(z), 0.0, a, limit=200)[0]
    tail = integrate.quad(pdf, a, np.inf, limit=200)[0]
    return 2.0 * (inner + a**p * tail)


def stream_shift(pdf, a: float, null_pdf=stats.norm.pdf) -> dict:
    """Mean, standard deviation and signal-to-noise ratio of both normalised streams.

    ``pdf`` is the post-change density of the standardised return and
    ``null_pdf`` the pre-change one; both are clipped at ``a`` (in pre-change
    standard deviations). Returns ``E[Y]``, ``sd(Y)`` and
    ``snr = (E[Y] - 1) / sd(Y)`` for ``Y = X / m0`` and ``Y = Q / q0``. The
    signal-to-noise ratio is the per-step edge a bettor can exploit; the
    expected log growth of the best constant bet is about ``snr^2 / 2`` per step.
    """
    m0 = _clipped_moment(null_pdf, 2, a)
    q0 = _clipped_moment(null_pdf, 4, a)
    m1, m2 = _clipped_moment(pdf, 2, a), _clipped_moment(pdf, 4, a)
    q2 = _clipped_moment(pdf, 8, a)
    out = {}
    for name, mean, second, null in (("var", m1, m2, m0), ("q", m2, q2, q0)):
        mu = mean / null
        sd = np.sqrt(max(second - mean**2, 0.0)) / null
        out[f"mean_{name}"] = mu
        out[f"sd_{name}"] = sd
        out[f"snr_{name}"] = (mu - 1.0) / sd
    return out


def student_pdf(df: float):
    """Density of a Student t with ``df`` degrees of freedom scaled to unit variance."""
    s = np.sqrt((df - 2.0) / df)
    return lambda z: stats.t.pdf(z / s, df) / s


def mixture_pdf(d: float):
    """Density of a scale mixture: sd^2 is ``1 - d`` or ``1 + d`` with equal chance (mean variance 1)."""
    lo, hi = np.sqrt(1.0 - d), np.sqrt(1.0 + d)
    return lambda z: 0.5 * stats.norm.pdf(z, scale=lo) + 0.5 * stats.norm.pdf(z, scale=hi)


def scaled_normal_pdf(ratio: float):
    return lambda z: stats.norm.pdf(z, scale=ratio)


# ----------------------------------------------------------------------------
# Simulation of distribution changes with the variance held fixed
# ----------------------------------------------------------------------------


def shape_change_returns(
    reps: int,
    n: int,
    nu: int | None,
    kind: str,
    param: float,
    rng: np.random.Generator,
    p_jump: float = 0.0,
    jump_scale: float = 8.0,
    pre_df: float | None = None,
) -> np.ndarray:
    """Unit-variance returns that switch distribution at index ``nu``.

    Pre-change innovations are standard normal (or unit-variance Student t
    with ``pre_df`` degrees of freedom). From ``nu`` on, ``kind`` selects:

    * ``"sd"``: normal with standard deviation ``param`` (a pure variance change);
    * ``"student"``: Student t with ``param`` degrees of freedom, scaled to unit
      variance (fatter tails, same variance);
    * ``"mixture"``: normal with variance ``1 - param`` or ``1 + param`` chosen
      independently each step with equal chance (vol of vol, same mean variance);
    * ``"sv"``: normal with a slow log-volatility AR(1) factor of stationary sd
      ``param`` and half-life 60 steps, scaled so the mean variance stays 1.

    Jumps (probability ``p_jump`` per step, size ``N(0, jump_scale^2)``) are
    added over the whole path, so they are part of the null as in E4.
    """
    if pre_df is None:
        r = rng.standard_normal((reps, n))
    else:
        r = rng.standard_t(pre_df, (reps, n)) * np.sqrt((pre_df - 2.0) / pre_df)
    if nu is not None:
        m = n - nu
        if kind == "sd":
            r[:, nu:] = rng.standard_normal((reps, m)) * param
        elif kind == "student":
            r[:, nu:] = rng.standard_t(param, (reps, m)) * np.sqrt((param - 2.0) / param)
        elif kind == "mixture":
            s2 = 1.0 + param * rng.choice([-1.0, 1.0], (reps, m))
            r[:, nu:] = rng.standard_normal((reps, m)) * np.sqrt(s2)
        elif kind == "sv":
            phi = 0.5 ** (1.0 / 60.0)
            z = np.empty((reps, m))
            z[:, 0] = param * rng.standard_normal(reps)
            e = rng.standard_normal((reps, m)) * param * np.sqrt(1 - phi**2)
            for k in range(1, m):
                z[:, k] = phi * z[:, k - 1] + e[:, k]
            r[:, nu:] = rng.standard_normal((reps, m)) * np.exp(z - param**2)
        else:
            raise ValueError(f"unknown kind {kind!r}")
    if p_jump > 0:
        r += (rng.random((reps, n)) < p_jump) * rng.normal(0.0, jump_scale, (reps, n))
    return r


# ----------------------------------------------------------------------------
# Streaming monitor on the quarticity stream
# ----------------------------------------------------------------------------


class QuarticityMonitor(VolatilityMonitor):
    """Calibrate, monitor, alarm, recalibrate loop on ``min(r^4, tau^4)``.

    Same protocol and guarantees as ``avmonitor.monitor.VolatilityMonitor``:
    each side runs at level ``2 A`` so that the first alarm of either side has
    average run length at least ``A`` and ``P(false alarm within m steps) <= m / A``
    under the null (with the true ``q0``). ``null_band = b`` makes the sides bet
    against "truncated quarticity within a factor ``1 +/- b`` of the baseline".
    Note that a pure variance change by a factor ``v`` moves the quarticity by
    about ``v^2``, so a band of ``b`` on quarticity is a tighter tolerance on the
    variance than the same band on the variance stream.
    """

    def __init__(
        self,
        threshold: float = 28800.0,
        c: float = 3.0,
        cal_len: int = 1440,
        null_band: float = 0.0,
        up_grid=Q_UP_GRID,
        down_fractions=DEFAULT_DOWN_FRACTIONS,
        scale_floor: float = 0.0,
    ):
        super().__init__(threshold, c, cal_len, null_band, up_grid, down_fractions, scale_floor)

    def _calibrate(self):
        r = np.asarray(self.buffer)
        r = r[np.isfinite(r)]
        self.scale = max(float(np.median(np.abs(r)) * MAD_NORMAL_SCALE), self.scale_floor)
        self.tau2 = (self.c * self.scale) ** 2
        self.tau4 = self.tau2**2
        q0 = float(np.mean(np.minimum(r**4, self.tau4)))
        self.q0 = q0
        self.q_up = q0 * (1 + self.band)
        self.q_dn = q0 * (1 - self.band)
        b = self.tau4 / self.q_dn
        self.lam_dn = self.down_fractions / (b - 1.0)
        self.start = self.t + 1
        self.calibrated = True

    def update(self, r: float) -> Event | None:
        self.t += 1
        if not self.calibrated:
            self.buffer.append(r)
            if len(self.buffer) >= self.cal_len:
                self._calibrate()
            return None
        if not np.isfinite(r):
            return None
        q = min(r**4, self.tau4)
        self.S_up = (1.0 + self.up_grid * (q / self.q_up - 1.0)) * (self.S_up + 1.0)
        self.S_dn = (1.0 + self.lam_dn * (1.0 - q / self.q_dn)) * (self.S_dn + 1.0)
        mu, md = self.S_up.mean(), self.S_dn.mean()
        ev = None
        if mu >= 2 * self.A:
            ev = Event(self.t, "up", self.start, self.scale, float(mu))
        elif md >= 2 * self.A:
            ev = Event(self.t, "down", self.start, self.scale, float(md))
        if ev is not None:
            self.events.append(ev)
            self._reset()
        return ev
