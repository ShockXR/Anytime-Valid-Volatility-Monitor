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
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from .robust import MAD_NORMAL_SCALE


@dataclass
class Event:
    index: int  # position in the stream of the alarm
    direction: str  # "up" or "down"
    monitored_since: int  # first monitored index of this segment
    baseline_scale: float  # robust sd of the calibration window
    stat: float  # value of the side's statistic at the alarm


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
    ):
        self.A = float(threshold)
        self.c = c
        self.cal_len = cal_len
        self.band = null_band
        self.up_grid = np.asarray(up_grid, float)
        self.down_fractions = np.asarray(down_fractions, float)
        self.scale_floor = scale_floor
        self.t = -1
        self.events: list[Event] = []
        self._reset()

    def _reset(self):
        self.buffer: list[float] = []
        self.calibrated = False
        self.S_up = np.zeros(len(self.up_grid))
        self.S_dn = np.zeros(len(self.down_fractions))

    def _calibrate(self):
        r = np.asarray(self.buffer)
        r = r[np.isfinite(r)]
        self.scale = max(float(np.median(np.abs(r)) * MAD_NORMAL_SCALE), self.scale_floor)
        self.tau2 = (self.c * self.scale) ** 2
        m0 = float(np.mean(np.minimum(r**2, self.tau2)))
        self.m_up = m0 * (1 + self.band)
        self.m_dn = m0 * (1 - self.band)
        b = self.tau2 / self.m_dn
        self.lam_dn = self.down_fractions / (b - 1.0)
        self.start = self.t + 1
        self.calibrated = True

    def update(self, r: float) -> Event | None:
        """Feed one return; returns an ``Event`` if an alarm fires at this step."""
        self.t += 1
        if not self.calibrated:
            self.buffer.append(r)
            if len(self.buffer) >= self.cal_len:
                self._calibrate()
            return None
        if not np.isfinite(r):
            return None
        x = min(r * r, self.tau2)
        self.S_up = (1.0 + self.up_grid * (x / self.m_up - 1.0)) * (self.S_up + 1.0)
        self.S_dn = (1.0 + self.lam_dn * (1.0 - x / self.m_dn)) * (self.S_dn + 1.0)
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
    """Same protocol with Gaussian likelihood ratio bets (benchmark, not robust)."""

    def __init__(self, *args, up_ratios=(1.25, 1.5, 2.0), down_ratios=(0.5, 0.7, 0.8), **kw):
        super().__init__(*args, **kw)
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
            self.buffer.append(r)
            if len(self.buffer) >= self.cal_len:
                self._calibrate()
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
