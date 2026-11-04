"""Streaming volatility regime monitor.

``VolatilityMonitor`` runs a calibrate, monitor, alarm, recalibrate loop on a
stream of (deseasonalised) returns:

1. Calibrate: collect ``cal_len`` returns, set ``scale`` = scaled MAD about zero,
   ``tau = c * scale`` and ``m0`` = mean of ``min(r^2, tau^2)``.
2. Monitor: update a Shiryaev-Roberts betting e-detector for an increase and
   one for a decrease of the truncated variance (see ``detectors.py``), each with
   alarm level ``2 A``. Because ``(M_up + M_down) / 2`` is itself an e-detector,
   the first alarm of either side has average run length at least ``A`` under
   the null, and the chance of any false alarm within m steps is at most m / A.
3. On an alarm, record an event and go back to step 1 with fresh data.

``null_band = b`` makes each side bet against "truncated variance within a
factor (1 +/- b) of the baseline" instead of exact equality.

``c_down`` (default: equal to ``c``) gives the decrease side its own truncation
constant: it uses ``X_dn = min(r^2, tau_dn^2)`` with ``tau_dn = c_down * scale``
and its own null mean ``m0_dn`` = mean of ``X_dn``, while the increase side
keeps ``tau = c * scale`` and ``m0``. A smaller clip on the decrease side raises
the cap on the evidence one step can carry against a fall, about ``1 / c_down^2``
nats (``docs/theory_delay.md``, Corollary 4 and section 6.3). Each side is an
e-detector for its own clipped statistic, so the two-sided guarantee needs
only each side's own null: "the mean of ``min(r^2, tau^2)`` has not risen" for
the increase side and "the mean of ``min(r^2, tau_dn^2)`` has not fallen" for
the decrease side. That guarantee is proved when each side's baseline is the
true null mean, or with ``baseline="cs"`` (below); with the plug-in baseline it
is an approximation, as for ``c_down = c``. ``c_down = c`` reproduces the
default exactly. Each entry of ``segments`` records ``tau2_dn`` and ``m0_dn``.

Baselines
---------
``baseline="plugin"`` (default) is the original behaviour described above. Its
guarantee assumes the true baseline. ``margin_z = z`` (plug-in only) moves each
side's null value by z standard errors in the conservative direction.

``baseline="cs"`` makes the guarantee hold with an ESTIMATED baseline
(``docs/math.md``, section 6). The calibration window is split: ``tau`` comes
from its first ``round(tau_frac * cal_len)`` returns, and a time-uniform
confidence sequence for the mean of ``X = min(r^2, tau^2)`` runs on the rest
(``confseq.StreamingEBConfSeq``). The increase side bets against the upper
bound ``U``, the decrease side against the lower bound ``L``. If there is no
change, P(any false alarm within m monitored steps) <= ``delta + m / A``, where
``delta`` is split equally between the two bounds. With ``c_down != c`` there
are two sequences, both started after the same ``tau`` split: the upper bound
``U`` comes from a sequence on ``min(r^2, tau^2)`` and the lower bound ``L``
from a sequence on ``min(r^2, tau_dn^2)`` (bounded by ``tau_dn^2``), each at
level ``delta / 2``. Only one bound of each sequence is used, so the total
error is still ``delta``, and ``docs/math.md`` Theorem 6.4 applies to each side
for its own clipped statistic. Optional adaptive rule: with
``target_rel_width = w``, calibration continues past ``cal_len`` until
``(U - L) / mean <= w`` or the window reaches ``cal_max``; the time-uniform
sequence keeps the guarantee at this data-dependent length. (With two
sequences, both must meet the target, each with its own ``(U - L) / mean``.)

Degenerate windows
------------------
If a calibration window gives ``scale <= 0`` or ``m0 <= 0`` (more than half of
the returns are exactly zero, for example frozen prices or coarse price ticks,
or the window is all NaN), the window is extended, one return at a time, until
it is usable (and, once extended, holds at least ``cal_len / 2`` finite
returns) or reaches ``cal_max`` returns (default ``2 * cal_len``). At the cap
the scale falls back to the scaled median of the NONZERO absolute returns; if
there are none, calibration restarts with fresh data. Counts are kept in
``cal_stats``. Windows with ``scale > 0`` and ``m0 > 0`` are handled exactly as
before, so the default output is unchanged whenever the old code did not fail.
Per-segment calibration details are stored in ``segments``.

Seeing that the monitor is blind
--------------------------------
While it calibrates, the monitor cannot raise an alarm. A stream that is all
NaN, or a price that never moves, keeps it calibrating for ever: every window
reaches ``cal_max`` and is thrown away. ``status()`` returns a small dictionary
that shows this at any moment (steps spent calibrating and monitoring,
calibration restarts in a row, and whether the increase side of the current
segment can gain any evidence at all); ``run`` also warns at its end when the
monitor never finished calibrating. Neither changes any alarm.

Checked arguments
-----------------
``threshold`` must be a finite number above 1: the bound m / A says nothing
for A <= 1, and with A = inf a statistic that overflows to inf still alarms
(inf >= 2 inf). ``c`` and ``c_down`` must be finite and positive. A
calibration window whose returns are so large that ``(c * scale)^2``
overflows gives ``tau^2 = inf``, which the degenerate-window rule above
rejects, instead of an ``OverflowError``.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np

from .confseq import StreamingEBConfSeq
from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID, check_bets
from .robust import MAD_NORMAL_SCALE

def _sq(x) -> float:
    """``x ** 2`` exactly as before for every finite result; ``inf`` instead of ``OverflowError`` (audit R5-B-001)."""
    try:
        return x**2
    except OverflowError:  # a Python float raises here; numpy.float64 already gives inf
        return math.inf


def _as_return(r):
    """One return as a number, for the ``update`` methods of the subclasses (same rule as ``VolatilityMonitor``).

    A Python float or a ``numpy.float64`` (what ``run`` passes) is returned unchanged, so their arithmetic is
    exactly as before; anything ``float()`` cannot convert raises ``TypeError`` at once (audit A4-B-004).
    """
    if type(r) is float or type(r) is np.float64:
        return r
    try:
        return float(r)
    except (TypeError, ValueError):
        raise TypeError(f"a return must be one number, got {type(r).__name__}: {r!r}") from None


# Below this calibration length VolatilityMonitor warns (it still runs). With n returns the scaled MAD has a
# relative standard error of about 1.17 / sqrt(n) for normal data: 12 percent at n = 100, so the baseline
# variance m0 is off by about 23 percent, already as large as the smallest change the default bets target.
MIN_ADVISED_CAL_LEN = 100


@dataclass
class Event:
    """One alarm of a streaming monitor, returned by ``update`` and collected in ``events``."""

    index: int  # position in the stream of the alarm
    direction: str  # "up" or "down"
    monitored_since: int  # first monitored index of this segment
    baseline_scale: float  # robust sd of the calibration window
    stat: float  # value of the side's statistic at the alarm


def _number(v, name: str) -> float:
    """``float(v)``, with a ``ValueError`` that names the argument when ``v`` is not a number."""
    try:
        return float(v)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a finite number, got {v!r}") from None


def _usable(scale: float, tau2: float, m0: float | None = None) -> bool:
    vals = (scale, tau2) if m0 is None else (scale, tau2, m0)
    return all(math.isfinite(v) and v > 0 for v in vals)


class VolatilityMonitor:
    """Calibrate, monitor, alarm, recalibrate loop on a stream of returns (details in the module docstring).

    ``threshold`` is A; each side alarms when its statistic reaches 2 A, so with the true baseline
    P(any false alarm within the first m monitored steps of a segment) <= m / A. Each segment starts
    with a calibration window of ``cal_len`` returns (longer if it is degenerate) that sets ``scale``,
    ``tau = c * scale`` and the null mean ``m0``. With ``baseline="cs"`` the null values are confidence
    bounds instead, and the bound is ``delta + m / A``. Feed returns with ``update`` (one at a time) or
    ``run`` (an array). Alarms are ``Event`` objects in ``events``; calibration details are in
    ``segments``.
    """

    def __init__(
        self,
        threshold: float = 28800.0,
        c: float = 3.0,
        cal_len: int = 1440,
        null_band: float = 0.0,
        up_grid=DEFAULT_UP_GRID,
        down_fractions=DEFAULT_DOWN_FRACTIONS,
        scale_floor: float = 0.0,
        *,
        baseline: str = "plugin",
        delta: float = 0.01,
        margin_z: float = 0.0,
        tau_frac: float = 0.25,
        cs_n_plan: float | None = None,
        cs_c: float = 0.75,
        target_rel_width: float | None = None,
        cal_max: int | None = None,
        c_down: float | None = None,
    ):
        self.A = _number(threshold, "threshold")
        # A <= 1 makes the bound m / A empty, NaN never alarms (audit A1-B-001), and inf alarms once a statistic
        # overflows to inf, because inf >= 2 * inf (audit R5-B-002).
        if not (math.isfinite(self.A) and self.A > 1):
            raise ValueError(f"threshold must be a finite number above 1 (the false alarm bound is m / threshold), "
                             f"got {threshold!r}")
        self.c = c
        c_value = _number(c, "c")
        if not (math.isfinite(c_value) and c_value > 0):
            raise ValueError(f"c must be a finite positive number, got {c!r}")
        # The decrease side may use its own truncation constant (docs/theory_delay.md, 6.3). Each side is an
        # e-detector for its own truncated statistic, so the two-sided guarantee needs only each side's own
        # null; it is proved for a known baseline and for baseline="cs" (docs/math.md, Theorem 6.4).
        self.c_down = float(c) if c_down is None else _number(c_down, "c_down")
        if not (math.isfinite(self.c_down) and self.c_down > 0):
            raise ValueError(f"c_down must be a finite positive number, got {c_down!r}")
        if int(cal_len) < 2:  # a one-return window gives scale = |r| and a meaningless baseline (audit A1-B-002)
            raise ValueError("cal_len must be at least 2")
        if int(cal_len) < MIN_ADVISED_CAL_LEN:
            warnings.warn(
                f"cal_len={int(cal_len)} is short: the scaled MAD of n returns has a relative standard error of "
                f"about 1.17/sqrt(n) ({1.17 / math.sqrt(int(cal_len)):.0%} here), so the baseline variance is "
                f"off by about twice that and the false alarm bound m / threshold no longer describes the "
                f"monitor. Use at least {MIN_ADVISED_CAL_LEN} returns (the experiments use 120 to 1440).",
                UserWarning, stacklevel=2)
        self.cal_len = cal_len
        self.band = null_band
        self.up_grid = np.asarray(up_grid, float)
        self.down_fractions = np.asarray(down_fractions, float)
        check_bets(self.up_grid, self.down_fractions)
        if not null_band >= 0:
            raise ValueError("null_band must be nonnegative")
        if null_band >= 1:
            warnings.warn(f"null_band={null_band} is 1 or more, so the decrease side can never bet and the monitor "
                          "only looks for volatility increases (audit A4-B-003)", UserWarning, stacklevel=2)
        self.scale_floor = scale_floor
        if baseline not in ("plugin", "cs"):
            raise ValueError("baseline must be 'plugin' or 'cs'")
        self.baseline = baseline
        self.delta = float(delta)
        self.margin_z = float(margin_z)
        self.cal_max = int(cal_max) if cal_max is not None else 2 * int(cal_len)
        if self.cal_max < cal_len:
            raise ValueError("cal_max must be at least cal_len")
        if baseline == "cs":
            if not 0 < self.delta < 1:
                raise ValueError("delta must lie in (0, 1)")
            if margin_z:
                raise ValueError("margin_z applies to the plug-in baseline only")
            self.n_tau = int(round(tau_frac * cal_len))
            if not 1 <= self.n_tau < cal_len:
                raise ValueError("tau_frac must leave returns on both sides of the split")
            self.cs_n_plan = cal_len - self.n_tau if cs_n_plan is None else cs_n_plan
        elif target_rel_width is not None:
            raise ValueError("target_rel_width needs baseline='cs'")
        self.tau_frac = tau_frac
        self.cs_c = cs_c
        self.target_rel_width = target_rel_width
        self.t = -1
        self.events: list[Event] = []
        self.segments: list[dict] = []
        self.cal_stats = {"calibrations": 0, "extended": 0, "fallback": 0, "restarted": 0}
        self.cal_steps = 0  # returns fed to calibration so far (no alarm is possible during these steps)
        self.restarts_in_a_row = 0  # calibration restarts since the last finished calibration
        self._reset()

    def _reset(self):
        self.buffer: list[float] = []
        self.calibrated = False
        self.S_up = np.zeros(len(self.up_grid))
        self.S_dn = np.zeros(len(self.down_fractions))
        self._extended = False
        self._fallback = False
        self._cs = None
        self._cs_dn = None  # second sequence (on min(r^2, tau_dn^2)) when c_down != c with baseline="cs"
        self._n_tau_used = 0
        self._n_finite = 0

    # ------------------------------------------------------------------
    # Calibration
    # ------------------------------------------------------------------

    def _plugin_stats(self, r: np.ndarray):
        if r.size == 0:
            return 0.0, 0.0, 0.0
        scale = max(float(np.median(np.abs(r)) * MAD_NORMAL_SCALE), self.scale_floor)
        tau2 = _sq(self.c * scale)
        m0 = float(np.mean(np.minimum(r**2, tau2)))
        return scale, tau2, m0

    def _fallback_scale(self, r: np.ndarray) -> float:
        nz = np.abs(r[r != 0])
        if nz.size == 0:
            return 0.0
        return max(float(np.median(nz) * MAD_NORMAL_SCALE), self.scale_floor)

    def _restart(self):
        self.cal_stats["restarted"] += 1
        self.restarts_in_a_row += 1
        self.buffer = []
        self._cs = None
        self._cs_dn = None
        self._n_tau_used = 0

    def _feed_calibration(self, r: float):
        self.cal_steps += 1
        self.buffer.append(r)
        if self.baseline == "plugin":
            if len(self.buffer) >= self.cal_len:
                self._calibrate()
        else:
            self._feed_cs(r)

    def _calibrate(self):
        """Plug-in calibration on the current buffer (extends the window if it is degenerate)."""
        r = np.asarray(self.buffer)
        r = r[np.isfinite(r)]
        scale, tau2, m0 = self._plugin_stats(r)
        usable = _usable(scale, tau2, m0)
        # A window that is already being extended must also hold cal_len / 2 finite returns.
        if not (usable and (not self._extended or r.size >= self.cal_len // 2)):
            if len(self.buffer) < self.cal_max:
                self._extended = True
                return
            if not usable:
                scale = self._fallback_scale(r)
                tau2 = _sq(self.c * scale)
                m0 = float(np.mean(np.minimum(r**2, tau2))) if r.size else 0.0
                if not _usable(scale, tau2, m0):
                    self._restart()
                    return
                self._fallback = True
        self._n_finite = int(r.size)
        self.scale = scale
        self.tau2 = tau2
        if self.c_down == self.c:
            tau2_dn, m0_dn = tau2, m0
        else:
            tau2_dn = _sq(self.c_down * scale)
            m0_dn = float(np.mean(np.minimum(r**2, tau2_dn))) if r.size else 0.0
        self.tau2_dn = tau2_dn
        if self.margin_z:
            x = np.minimum(r**2, tau2)
            se = float(np.sqrt(x.var(ddof=1) / x.size)) if x.size > 1 else float("inf")
            xd = np.minimum(r**2, tau2_dn)
            se_dn = float(np.sqrt(xd.var(ddof=1) / xd.size)) if xd.size > 1 else float("inf")
            m_up = (m0 + self.margin_z * se) * (1 + self.band)
            m_dn = (m0_dn - self.margin_z * se_dn) * (1 - self.band)
        else:
            m_up = m0 * (1 + self.band)
            m_dn = m0_dn * (1 - self.band)
        self._finish(m0, m_up, m_dn, None, None, m0_dn)

    def _feed_cs(self, r: float):
        if self._cs is None:  # first part: choose tau
            if len(self.buffer) < self.n_tau:
                return
            rr = np.asarray(self.buffer)
            rr = rr[np.isfinite(rr)]
            scale, tau2, _ = self._plugin_stats(rr)
            usable = _usable(scale, tau2)
            if not (usable and rr.size >= max(1, self.n_tau // 2)):
                if len(self.buffer) < self.cal_max - (self.cal_len - self.n_tau):
                    self._extended = True
                    return
                if not usable:
                    scale = self._fallback_scale(rr)
                    tau2 = _sq(self.c * scale)
                    if not _usable(scale, tau2):
                        self._restart()
                        return
                    self._fallback = True
            self.tau_scale = scale
            self.tau2 = tau2
            self.tau2_dn = tau2 if self.c_down == self.c else _sq(self.c_down * scale)
            half = self.delta / 2
            self._cs = StreamingEBConfSeq(tau2, half, half, n_plan=self.cs_n_plan, c=self.cs_c)
            if self.c_down != self.c:
                # The upper bound is read from self._cs and the lower bound from self._cs_dn, each at level
                # delta / 2, so the total error stays at delta (docs/math.md, Theorem 6.4, for each side).
                self._cs_dn = StreamingEBConfSeq(self.tau2_dn, half, half, n_plan=self.cs_n_plan, c=self.cs_c)
            self._n_tau_used = len(self.buffer)
            return
        if np.isfinite(r):  # second part: confidence sequence for E[X] (and for E[X_dn])
            self._cs.update(min(r * r, self.tau2))
            if self._cs_dn is not None:
                self._cs_dn.update(min(r * r, self.tau2_dn))
        n_cs = len(self.buffer) - self._n_tau_used
        if n_cs < self.cal_len - self.n_tau:
            return
        if self.target_rel_width is not None and len(self.buffer) < self.cal_max:
            if not self._cs.rel_width <= self.target_rel_width:
                return
            if self._cs_dn is not None and not self._cs_dn.rel_width <= self.target_rel_width:
                return
        if self._cs.t < max(1, (self.cal_len - self.n_tau) // 2):  # too few finite returns in the second part
            if len(self.buffer) < self.cal_max or self._cs.t == 0:
                # Keep collecting. A sequence that has started is never restarted, so each segment uses exactly
                # one confidence sequence per side and its error probability stays at delta (docs/math.md, 6.3).
                self._extended = True
                return
        rr = np.asarray(self.buffer)
        rr = rr[np.isfinite(rr)]
        self._n_finite = int(rr.size)
        full = float(np.median(np.abs(rr)) * MAD_NORMAL_SCALE)
        self.scale = full if full > 0 else self.tau_scale
        cs_dn = self._cs if self._cs_dn is None else self._cs_dn
        U, L = self._cs.upper, cs_dn.lower
        self._finish(self._cs.mean, U * (1 + self.band), L * (1 - self.band), U, L, cs_dn.mean)

    def _finish(self, m0, m_up, m_dn, upper, lower, m0_dn):
        self.m0 = m0
        self.m0_dn = m0_dn
        self.m_up = m_up
        self.m_dn = m_dn
        b = self.tau2_dn / self.m_dn if self.m_dn > 0 else math.inf
        self.down_active = bool(math.isfinite(b) and b > 1.0)
        if self.down_active:
            self.lam_dn = self.down_fractions / (b - 1.0)
        else:  # no admissible decrease bet: that side's statistic stays at 0 (a valid e-detector)
            self.lam_dn = np.zeros(len(self.down_fractions))
        self.start = self.t + 1
        self.calibrated = True
        self.restarts_in_a_row = 0
        self.cal_stats["calibrations"] += 1
        self.cal_stats["extended"] += int(self._extended)
        self.cal_stats["fallback"] += int(self._fallback)
        self.segments.append(dict(start=self.start, cal_len=len(self.buffer), n_finite=self._n_finite,
                                  baseline=self.baseline,
                                  n_tau=self._n_tau_used if self.baseline == "cs" else None,
                                  scale=self.scale, c=self.c, c_down=self.c_down, tau2=self.tau2,
                                  tau2_dn=self.tau2_dn, m0=m0, m0_dn=m0_dn, m_up=m_up, m_dn=m_dn,
                                  upper=upper, lower=lower, extended=self._extended, fallback=self._fallback,
                                  down_active=self.down_active))

    # ------------------------------------------------------------------
    # Monitoring
    # ------------------------------------------------------------------

    def update(self, r: float) -> Event | None:
        """Feed one return; returns an ``Event`` if an alarm fires at this step."""
        try:  # fail at once on a string, None or an array (audit A4-B-004), not when calibration ends
            r = float(r)
        except (TypeError, ValueError):
            raise TypeError(f"a return must be one number, got {type(r).__name__}: {r!r}") from None
        self.t += 1
        if not self.calibrated:
            self._feed_calibration(r)
            return None
        if not np.isfinite(r):
            return None
        x = min(r * r, self.tau2)
        self.S_up = (1.0 + self.up_grid * (x / self.m_up - 1.0)) * (self.S_up + 1.0)
        if self.down_active:
            xd = x if self.tau2_dn == self.tau2 else min(r * r, self.tau2_dn)
            self.S_dn = (1.0 + self.lam_dn * (1.0 - xd / self.m_dn)) * (self.S_dn + 1.0)
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

    def run(self, returns) -> list[Event]:
        """Feed ``returns`` in order through ``update``; returns ``self.events`` (all alarms so far)."""
        r_all = np.asarray(returns, float)
        for r in r_all:
            self.update(r)
        if r_all.size and not self.segments and not self.calibrated:
            n_fin = int(np.isfinite(r_all).sum())
            warnings.warn(
                f"the monitor never finished calibrating on these {r_all.size} returns ({n_fin} finite), so it "
                f"could not alarm at all: no alarm here does not mean no change. Calibration needs about "
                f"cal_len={self.cal_len} usable returns.", UserWarning, stacklevel=2)
        return self.events

    def status(self) -> dict:
        """Where the monitor stands now, to see whether it is blind (it changes nothing).

        Keys: ``step`` (index of the last return fed, -1 before any), ``state`` ("calibrating" or "monitoring"),
        ``steps_calibrating`` and ``steps_monitoring`` (totals so far; no alarm is possible while calibrating),
        ``window`` (returns in the current calibration window, 0 while monitoring), ``restarts_in_a_row`` and
        ``restarts_total`` (windows thrown away because they were all missing or all zero), ``calibrations``,
        ``alarms``, and for a monitored segment ``tau2_over_m_up`` and ``up_can_gain`` (False when ``tau^2`` is
        not above the increase side's null value, so that side's bets can never win and it cannot alarm) and
        ``down_active``. ``blind`` is True while calibrating after at least one restart in a row, or while
        monitoring with neither side able to gain.
        """
        out = dict(step=self.t, state="monitoring" if self.calibrated else "calibrating",
                   steps_calibrating=self.cal_steps, steps_monitoring=self.t + 1 - self.cal_steps,
                   window=0 if self.calibrated else len(self.buffer), restarts_in_a_row=self.restarts_in_a_row,
                   restarts_total=self.cal_stats["restarted"], calibrations=self.cal_stats["calibrations"],
                   alarms=len(self.events), tau2_over_m_up=None, up_can_gain=None, down_active=None)
        if self.calibrated:
            m_up, tau2 = getattr(self, "m_up", math.nan), getattr(self, "tau2", math.nan)
            ratio = tau2 / m_up if m_up > 0 else math.inf
            out.update(tau2_over_m_up=float(ratio), up_can_gain=bool(ratio > 1.0),
                       down_active=bool(getattr(self, "down_active", False)))
            out["blind"] = not (out["up_can_gain"] or out["down_active"])
        else:
            out["blind"] = self.restarts_in_a_row > 0
        return out


class GaussianLRMonitor(VolatilityMonitor):
    """Same protocol with Gaussian likelihood ratio bets (benchmark, not robust). Plug-in baseline only."""

    def __init__(self, *args, up_ratios=(1.25, 1.5, 2.0), down_ratios=(0.5, 0.7, 0.8), **kw):
        super().__init__(*args, **kw)
        if self.baseline != "plugin":
            raise ValueError("GaussianLRMonitor supports only baseline='plugin'")
        self.up_ratios = np.asarray(up_ratios, float)
        self.down_ratios = np.asarray(down_ratios, float)
        self.S_up = np.zeros(len(self.up_ratios))
        self.S_dn = np.zeros(len(self.down_ratios))

    def _reset(self):
        super()._reset()
        if hasattr(self, "up_ratios"):
            self.S_up = np.zeros(len(self.up_ratios))
            self.S_dn = np.zeros(len(self.down_ratios))

    def update(self, r: float) -> Event | None:
        r = _as_return(r)
        self.t += 1
        if not self.calibrated:
            self._feed_calibration(r)
            return None
        if not np.isfinite(r):
            return None
        with np.errstate(over="ignore"):
            z2 = _sq(r / self.scale)
        with np.errstate(over="ignore"):
            fu = np.exp(np.minimum(-np.log(self.up_ratios) + 0.5 * (1 - 1 / self.up_ratios**2) * z2, 700))
            fd = np.exp(np.minimum(-np.log(self.down_ratios) + 0.5 * (1 - 1 / self.down_ratios**2) * z2, 700))
            self.S_up = fu * (self.S_up + 1.0)
            self.S_dn = fd * (self.S_dn + 1.0)
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
