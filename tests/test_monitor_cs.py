"""Tests for the monitor: the default output is unchanged, degenerate windows, and the cs baseline.

The legacy classes below are a verbatim copy (renamed) of ``VolatilityMonitor``
and ``GaussianLRMonitor`` as of commit bf67e2b, before the cs baseline and the
degenerate-window fix were added. The regression tests require the new code to
give exactly the same events on streams where the old code worked.
"""

from dataclasses import dataclass

import numpy as np
import pytest

from avmonitor.confseq import cs_baseline
from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE


# ----------------------------------------------------------------------------
# Legacy implementation (verbatim, renamed)
# ----------------------------------------------------------------------------


@dataclass
class _LegacyEvent:
    index: int
    direction: str
    monitored_since: int
    baseline_scale: float
    stat: float


class _LegacyVolatilityMonitor:
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
        self.events: list[_LegacyEvent] = []
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

    def update(self, r: float) -> _LegacyEvent | None:
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
            ev = _LegacyEvent(self.t, "up", self.start, self.scale, float(mu))
        elif md >= 2 * self.A:
            ev = _LegacyEvent(self.t, "down", self.start, self.scale, float(md))
        if ev is not None:
            self.events.append(ev)
            self._reset()
        return ev

    def run(self, returns) -> list[_LegacyEvent]:
        for r in np.asarray(returns, float):
            self.update(r)
        return self.events


class _LegacyGaussianLRMonitor(_LegacyVolatilityMonitor):
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

    def update(self, r: float) -> _LegacyEvent | None:
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
            ev = _LegacyEvent(self.t, "up", self.start, self.scale, float(mu))
        elif md >= 2 * self.A:
            ev = _LegacyEvent(self.t, "down", self.start, self.scale, float(md))
        if ev is not None:
            self.events.append(ev)
            self._reset()
        return ev


# ----------------------------------------------------------------------------
# Regression: default behaviour unchanged
# ----------------------------------------------------------------------------


def _stream(seed):
    rng = np.random.default_rng(seed)
    r = np.concatenate([rng.standard_normal(2000), 2 * rng.standard_normal(1500), 0.5 * rng.standard_normal(1500),
                        rng.standard_t(4, 1000)])
    rng.random(r.size)  # kept so that the stream matches the one used to record the reference values
    jumps = rng.random(r.size) < 0.004
    r[jumps] += rng.normal(0, 8, jumps.sum())
    r[rng.random(r.size) < 0.01] = np.nan
    return r


def _tuples(events):
    return [(e.index, e.direction, e.monitored_since, e.baseline_scale, e.stat) for e in events]


CONFIGS = [
    ("vm", dict(threshold=2000.0, cal_len=500)),
    ("vm", dict(threshold=5000.0, cal_len=390, null_band=0.25, scale_floor=0.1)),
    ("vm", dict(threshold=800.0, cal_len=200, c=2.5, up_grid=(0.1, 0.3), down_fractions=(0.5,))),
    ("glr", dict(threshold=2000.0, cal_len=500)),
    ("glr", dict(threshold=300.0, cal_len=250, null_band=0.5)),
]


@pytest.mark.parametrize("seed", [2024, 7, 8])
@pytest.mark.parametrize("kind,kw", CONFIGS)
def test_default_output_identical_to_legacy(seed, kind, kw):
    new_cls, old_cls = (VolatilityMonitor, _LegacyVolatilityMonitor) if kind == "vm" else (
        GaussianLRMonitor, _LegacyGaussianLRMonitor)
    r = _stream(seed)
    new, old = new_cls(**kw), old_cls(**kw)
    assert _tuples(new.run(r)) == _tuples(old.run(r))  # exact equality, floats included
    assert new.t == old.t and new.calibrated == old.calibrated


def test_default_output_matches_recorded_reference():
    """Values recorded from the legacy code (commit bf67e2b) on _stream(2024)."""
    r = _stream(2024)
    ev = VolatilityMonitor(threshold=2000.0, cal_len=500).run(r)
    ref = [(2010, "up", 500, 0.9487662550788026, 5218.968850763707),
           (2940, "down", 2511, 1.9637944599299129, 4225.249390306718),
           (3551, "down", 3441, 2.049299941495822, 4216.183045318844),
           (4863, "down", 4052, 0.502021941877509, 4168.606480507031),
           (5495, "up", 5364, 0.8641486648404929, 5145.843666188957)]
    got = _tuples(ev)
    assert [g[:3] for g in got] == [x[:3] for x in ref]
    np.testing.assert_allclose([g[3:] for g in got], [x[3:] for x in ref], rtol=1e-9)
    ev = VolatilityMonitor(threshold=5000.0, cal_len=390, null_band=0.25, scale_floor=0.1).run(r)
    assert [(e.index, e.direction, e.monitored_since) for e in ev] == [(2018, "up", 390), (3658, "down", 2409),
                                                                      (5023, "up", 4049)]
    ev = GaussianLRMonitor(threshold=2000.0, cal_len=500).run(r)
    assert [e.index for e in ev] == [914, 2004, 2658, 3201, 3717, 4304, 5007, 5550]


def test_margin_z_zero_is_default_and_positive_margin_is_more_conservative():
    r = _stream(3)
    a = VolatilityMonitor(threshold=2000.0, cal_len=500).run(r)
    b = VolatilityMonitor(threshold=2000.0, cal_len=500, margin_z=0.0).run(r)
    assert _tuples(a) == _tuples(b)
    mon = VolatilityMonitor(threshold=2000.0, cal_len=500, margin_z=1.0)
    mon.run(r[:600])
    seg = mon.segments[0]
    assert seg["m_up"] > seg["m0"] > seg["m_dn"]


# ----------------------------------------------------------------------------
# Degenerate calibration windows (MAD = 0, all NaN)
# ----------------------------------------------------------------------------


def test_many_zero_returns_do_not_crash_and_fall_back():
    rng = np.random.default_rng(40)
    r = rng.standard_normal(20000)
    r[rng.random(r.size) < 0.6] = 0.0  # 60% exact zeros: the MAD is 0 in every window
    for cls in (VolatilityMonitor, GaussianLRMonitor):
        mon = cls(threshold=10000.0, cal_len=500)
        mon.run(r)
        assert mon.cal_stats["calibrations"] >= 1
        assert mon.cal_stats["extended"] == mon.cal_stats["calibrations"]
        assert mon.cal_stats["fallback"] == mon.cal_stats["calibrations"]
        seg = mon.segments[0]
        assert seg["cal_len"] == 1000  # extended to the default cap 2 * cal_len
        assert 0.8 < mon.segments[0]["scale"] < 1.2  # scale of the nonzero returns
        assert seg["m0"] > 0
    mon = VolatilityMonitor(threshold=10000.0, cal_len=500, baseline="cs")
    mon.run(r)
    assert mon.cal_stats["calibrations"] >= 1 and mon.cal_stats["fallback"] >= 1
    assert mon.segments[0]["upper"] > 0


def test_frozen_stretch_extends_the_window():
    rng = np.random.default_rng(41)
    r = np.concatenate([np.zeros(400), rng.standard_normal(5000)])
    mon = VolatilityMonitor(threshold=10000.0, cal_len=390)
    mon.run(r)
    seg = mon.segments[0]
    assert mon.cal_stats["extended"] >= 1
    assert 390 < seg["cal_len"] <= 780 and seg["m0"] > 0 and seg["scale"] > 0
    mon = VolatilityMonitor(threshold=10000.0, cal_len=390, baseline="cs")
    mon.run(r)
    assert mon.segments[0]["scale"] > 0 and mon.cal_stats["extended"] >= 1


def test_all_nan_window_is_extended_or_restarted():
    rng = np.random.default_rng(42)
    r = np.concatenate([np.full(1500, np.nan), rng.standard_normal(3000)])
    for kw in (dict(), dict(baseline="cs")):
        mon = VolatilityMonitor(threshold=10000.0, cal_len=500, **kw)
        mon.run(r)  # no exception, no NaN baseline
        assert mon.cal_stats["calibrations"] >= 1 and mon.cal_stats["restarted"] >= 1
        seg = mon.segments[0]
        assert np.isfinite(seg["m_up"]) and seg["m_up"] > 0 and seg["start"] >= 1500
        assert seg["n_finite"] >= 250  # an extended window needs cal_len / 2 finite returns
    mon = GaussianLRMonitor(threshold=10000.0, cal_len=500)
    mon.run(r)
    assert mon.scale > 0


def test_frozen_window_with_scale_floor_extends():
    mon = VolatilityMonitor(threshold=10000.0, cal_len=100, scale_floor=0.5)
    mon.run(np.concatenate([np.zeros(150), np.random.default_rng(1).standard_normal(500)]))
    assert mon.cal_stats["extended"] == 1 and mon.segments[0]["m0"] > 0


# ----------------------------------------------------------------------------
# Confidence-sequence baseline
# ----------------------------------------------------------------------------


def test_cs_segment_matches_vectorised_baseline():
    rng = np.random.default_rng(43)
    r = rng.standard_normal(1440)
    mon = VolatilityMonitor(threshold=1e9, cal_len=1440, baseline="cs", delta=0.02)
    mon.run(np.append(r, 0.1))
    seg = mon.segments[0]
    b = cs_baseline(r[None, :], c=3.0, delta_upper=0.01, delta_lower=0.01, tau_frac=0.25)
    assert seg["n_tau"] == b.n_tau == 360
    np.testing.assert_allclose(seg["tau2"], b.tau[0] ** 2, rtol=1e-12)
    np.testing.assert_allclose([seg["upper"], seg["lower"], seg["m0"]], [b.upper[0], b.lower[0], b.mean[0]], rtol=1e-9)
    assert seg["lower"] < seg["m0"] < seg["upper"]
    assert mon.scale == pytest.approx(np.median(np.abs(r)) * MAD_NORMAL_SCALE)


def test_cs_monitor_is_conservative_under_null_and_detects():
    rng = np.random.default_rng(44)
    alarms = 0
    for _ in range(40):
        mon = VolatilityMonitor(threshold=10000.0, cal_len=1000, baseline="cs", delta=0.02)
        mon.run(rng.standard_normal(2500))
        alarms += len(mon.events) > 0
    assert alarms / 40 <= 0.02 + 1500 / 10000 + 0.08
    mon = VolatilityMonitor(threshold=3000.0, cal_len=1000, baseline="cs")
    ev = mon.run(np.concatenate([rng.standard_normal(1500), 2 * rng.standard_normal(1500)]))
    ups = [e for e in ev if e.direction == "up" and e.index >= 1500]
    assert ups and ups[0].index < 1800


def test_adaptive_calibration_length():
    rng = np.random.default_rng(45)
    r = rng.standard_normal(20000)
    tight = VolatilityMonitor(threshold=1e9, cal_len=400, baseline="cs", target_rel_width=0.5, cal_max=8000)
    tight.run(r)
    loose = VolatilityMonitor(threshold=1e9, cal_len=400, baseline="cs", target_rel_width=5.0, cal_max=8000)
    loose.run(r)
    never = VolatilityMonitor(threshold=1e9, cal_len=400, baseline="cs", target_rel_width=1e-6, cal_max=3000)
    never.run(r)
    assert loose.segments[0]["cal_len"] == 400  # width target met at the minimum length
    assert 400 < tight.segments[0]["cal_len"] < 8000
    s = tight.segments[0]
    assert (s["upper"] - s["lower"]) / s["m0"] <= 0.5
    assert never.segments[0]["cal_len"] == 3000  # capped


def test_invalid_options():
    with pytest.raises(ValueError):
        VolatilityMonitor(baseline="other")
    with pytest.raises(ValueError):
        VolatilityMonitor(target_rel_width=0.3)
    with pytest.raises(ValueError):
        VolatilityMonitor(baseline="cs", margin_z=1.0)
    with pytest.raises(ValueError):
        GaussianLRMonitor(baseline="cs")
    with pytest.raises(ValueError):
        VolatilityMonitor(cal_len=100, cal_max=50)


def test_recalibration_bound_on_alarm_count():
    """docs/math.md Corollary 6.6: E[alarms in m monitored steps] <= (delta + m / A) / (1 - delta) with restarts."""
    rng = np.random.default_rng(46)
    alarms, monitored, bound = 0, 0, 0.0
    for _ in range(2):
        mon = VolatilityMonitor(threshold=100.0, cal_len=200, baseline="cs", delta=0.02)
        r = rng.standard_normal(50000)
        mon.run(r)
        cal = sum(s["cal_len"] for s in mon.segments) + (0 if mon.calibrated else len(mon.buffer))
        m = r.size - cal
        alarms += len(mon.events)
        monitored += m
        bound += (0.02 + m / 100.0) / (1 - 0.02)
    assert monitored > 50000 and alarms <= bound
