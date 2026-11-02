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
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .confseq import StreamingEBConfSeq
from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from .robust import MAD_NORMAL_SCALE


@dataclass
class Event:
    index: int  # position in the stream of the alarm
    direction: str  # "up" or "down"
    monitored_since: int  # first monitored index of this segment
    baseline_scale: float  # robust sd of the calibration window
    stat: float  # value of the side's statistic at the alarm


def _usable(scale: float, tau2: float, m0: float | None = None) -> bool:
    vals = (scale, tau2) if m0 is None else (scale, tau2, m0)
    return all(math.isfinite(v) and v > 0 for v in vals)


class VolatilityMonitor:
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
        self.A = float(threshold)
        self.c = c
        # The decrease side may use its own truncation constant (docs/theory_delay.md, 6.3). Each side is an
        # e-detector for its own truncated statistic, so the two-sided guarantee needs only each side's own
        # null; it is proved for a known baseline and for baseline="cs" (docs/math.md, Theorem 6.4).
        self.c_down = float(c) if c_down is None else float(c_down)
        if self.c_down <= 0:
            raise ValueError("c_down must be positive")
        self.cal_len = cal_len
        self.band = null_band
        self.up_grid = np.asarray(up_grid, float)
        self.down_fractions = np.asarray(down_fractions, float)
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
        tau2 = (self.c * scale) ** 2
        m0 = float(np.mean(np.minimum(r**2, tau2)))
        return scale, tau2, m0

    def _fallback_scale(self, r: np.ndarray) -> float:
        nz = np.abs(r[r != 0])
        if nz.size == 0:
            return 0.0
        return max(float(np.median(nz) * MAD_NORMAL_SCALE), self.scale_floor)

    def _restart(self):
        self.cal_stats["restarted"] += 1
        self.buffer = []
        self._cs = None
        self._cs_dn = None
        self._n_tau_used = 0

    def _feed_calibration(self, r: float):
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
                tau2 = (self.c * scale) ** 2
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
            tau2_dn = (self.c_down * scale) ** 2
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
                    tau2 = (self.c * scale) ** 2
                    if not _usable(scale, tau2):
                        self._restart()
                        return
                    self._fallback = True
            self.tau_scale = scale
            self.tau2 = tau2
            self.tau2_dn = tau2 if self.c_down == self.c else (self.c_down * scale) ** 2
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
        for r in np.asarray(returns, float):
            self.update(r)
        return self.events


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
        self.t += 1
        if not self.calibrated:
            self._feed_calibration(r)
            return None
        if not np.isfinite(r):
            return None
        z2 = (r / self.scale) ** 2
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
