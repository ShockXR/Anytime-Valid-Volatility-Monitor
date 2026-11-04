import numpy as np
import pytest

from avmonitor import detectors as D
from avmonitor import simulate as S
from avmonitor.robust import truncated_moment_factor

TAU = 3.0
M0 = truncated_moment_factor(TAU, 2)


def _x(r):
    return np.minimum(r**2, TAU**2)


def test_null_average_run_length_at_least_threshold():
    rng = np.random.default_rng(11)
    p = S.change_returns(400, 2000, rng=rng)
    A = 50.0
    for direction in ("up", "down"):
        res = D.betting_edetector(_x(p.r), M0, direction, upper=TAU**2, threshold=A)
        rl = np.where(res.alarm >= 0, res.alarm + 1, 2000)
        # E[run length] >= A; allow Monte Carlo error
        assert rl.mean() > 0.85 * A


def test_false_alarm_within_horizon_bound():
    rng = np.random.default_rng(12)
    p = S.change_returns(2000, 20, rng=rng)
    A = 200.0
    res = D.betting_edetector(_x(p.r), M0, "up", threshold=A)
    # P(alarm within m steps) <= m / A
    assert (res.alarm >= 0).mean() <= 20 / A + 0.02


def test_detects_variance_increase_and_decrease():
    rng = np.random.default_rng(13)
    up = S.change_returns(100, 1500, nu=500, ratio=2.0, rng=rng)
    dn = S.change_returns(100, 1500, nu=500, ratio=0.5, rng=rng)
    ru = D.betting_edetector(_x(up.r), M0, "up", threshold=1e4)
    rd = D.betting_edetector(_x(dn.r), M0, "down", upper=TAU**2, threshold=1e4)
    su = D.summarize_alarms(ru.alarm, 500, 1500)
    sd = D.summarize_alarms(rd.alarm, 500, 1500)
    assert su["detected"] > 0.9 and su["median_delay"] < 100
    assert sd["detected"] > 0.9 and sd["median_delay"] < 200


def test_down_detector_needs_bound_and_bets_are_admissible():
    x = np.ones((1, 10))
    with pytest.raises(ValueError):
        D.betting_edetector(x, 1.0, "down")
    res = D.betting_edetector(x * 9.0, 1.0, "down", upper=9.0)
    assert np.all(res.info["grid"] < res.info["lam_max"])
    with pytest.raises(ValueError):
        D.betting_edetector(x, 1.0, "up", grid=[1.5])


def test_adaptive_moments_are_predictable():
    rng = np.random.default_rng(3)
    y = rng.random((1, 300))
    m1, m2 = D._predictable_local_moments(y, 20)
    y2 = y.copy()
    y2[0, 150] = 1e6
    m1b, _ = D._predictable_local_moments(y2, 20)
    assert m1[0, 150] == m1b[0, 150]
    assert m1[0, 151] != m1b[0, 151]


def test_cusum_never_exceeds_sr():
    rng = np.random.default_rng(4)
    p = S.change_returns(1, 3000, nu=1500, ratio=1.5, rng=rng)
    sr = D.betting_edetector(_x(p.r), M0, "up", threshold=np.inf, return_path=True)
    cu = D.betting_edetector(_x(p.r), M0, "up", threshold=np.inf, kind="CUSUM", return_path=True)
    assert np.all(cu.stat_path <= sr.stat_path + 1e-9)


def test_calibrate_band():
    rng = np.random.default_rng(5)
    b = D.calibrate(rng.standard_normal(10000), c=3.0, null_band=0.2)
    assert abs(b.scale - 1) < 0.05
    assert abs(b.m0_up / b.m0 - 1.2) < 1e-9
    assert abs(b.m0_down / b.m0 - 0.8) < 1e-9


def test_streaming_monitor_matches_guarantee_and_detects():
    from avmonitor.monitor import VolatilityMonitor

    rng = np.random.default_rng(21)
    # Null: 40 streams, 1500 monitored steps, A = 10000 so P(any alarm) <= 0.15 with the true m0.
    # The baseline is estimated from 2000 steps, so allow some slack for estimation error.
    alarms = 0
    for _ in range(40):
        mon = VolatilityMonitor(threshold=10000, cal_len=2000)
        mon.run(rng.standard_normal(3500))
        alarms += len(mon.events) > 0
    assert alarms / 40 <= 0.15 + 0.12
    # Change: sd doubles after 1500 steps
    mon = VolatilityMonitor(threshold=3000, cal_len=1000)
    ev = mon.run(np.concatenate([rng.standard_normal(1500), 2 * rng.standard_normal(1500)]))
    ups = [e for e in ev if e.direction == "up" and e.index >= 1500]
    assert ups and ups[0].index < 1700
