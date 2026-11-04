"""Regression tests for independent review round 9 (docs/review/independent_review_round9_code.md).

Proposed as tests/test_review_round9.py together with proposed_fixes.patch. The file name here does not match
pytest's test_*.py pattern on purpose: against the unpatched package every test fails, by design.
"""

import numpy as np
import pandas as pd
import pytest

from avmonitor import detectors as D
from avmonitor.design import sr_steps_constant_factor
from avmonitor.monitor import VolatilityMonitor
from avmonitor.periodic import PeriodicVolatilityMonitor
from avmonitor.robustness import trunc_mean_normal
from avmonitor.running import two_sided_sr

M0 = float(trunc_mean_normal(3.0))


@pytest.mark.parametrize("kind", ["SR", "CUSUM"])
def test_m1_missing_values_leave_the_batch_statistic_unchanged(kind):
    x = np.r_[np.full(50, np.nan), 1.0]  # no data, then one observation with Y = 1
    d = D.betting_edetector(x, 1.0, "up", threshold=10.0, kind=kind, return_path=True)
    assert d.alarm[0] == -1 and d.stat_path[-1] == pytest.approx(1.0)  # one chain with factor 1
    a = D.betting_edetector(x, 1.0, "up", threshold=10.0, adaptive_h=5, return_path=True)
    assert a.alarm[0] == -1


def test_m1_batch_and_streaming_agree_with_missing_returns():
    A = 60.0
    for seed in range(8):
        rng = np.random.default_rng(seed)
        r = np.concatenate([rng.standard_normal(800), 1.6 * rng.standard_normal(400), 0.6 * rng.standard_normal(400)])
        r[rng.random(r.size) < 0.1] = np.nan
        stream = [(e.index, e.direction) for e in
                  PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=3.0, m0=M0).run(r)]
        x = np.minimum(r**2, 9.0)
        batch, s = [], 0
        while s < r.size:
            up = D.betting_edetector(x[s:], M0, "up", threshold=2 * A)
            dn = D.betting_edetector(x[s:], M0, "down", upper=9.0, threshold=2 * A)
            a = int(D.two_sided_alarm(up, dn)[0])
            if a < 0:
                break
            batch.append((s + a, "up" if up.alarm[0] == a else "down"))
            s += a + 1
        assert stream == batch and len(stream) > 3


def test_m1_running_two_sided_sr_skips_a_missing_value():
    rng = np.random.default_rng(0)
    x = np.minimum((3 * rng.standard_normal(500)) ** 2, 9.0)
    ref = two_sided_sr(x[None, 1:], 9.0, M0, M0, [100.0])[100.0]
    x2 = x.copy()
    x2[0] = np.nan
    got = two_sided_sr(x2[None], 9.0, M0, M0, [100.0])[100.0]
    assert ref[0][0] >= 0 and got[0][0] == ref[0][0] + 1 and got[1][0] == (ref[1][0] + 1 if ref[1][0] >= 0 else -1)


def test_m2_calibrate_refuses_a_degenerate_window():
    with pytest.raises(ValueError):
        D.calibrate(np.r_[np.zeros(300), 0.001 * np.ones(200)])
    b = D.calibrate(np.r_[np.zeros(300), 0.001 * np.ones(200)], floor=0.003)
    assert b.tau == 0.003 and b.m0 > 0


def test_m3_bets_that_can_make_a_factor_negative_are_refused():
    with pytest.raises(ValueError):
        VolatilityMonitor(up_grid=(1.5,))
    with pytest.raises(ValueError):
        VolatilityMonitor(down_fractions=(1.2,))
    with pytest.raises(ValueError):
        VolatilityMonitor(null_band=-0.1)
    with pytest.raises(ValueError):
        D.betting_edetector(np.ones(10), 1.0, "up", adaptive_h=3, lam_cap=2.0)
    with pytest.raises(ValueError):
        D.betting_edetector(np.ones(10), 1.0, "up", adaptive_h=3, adaptive_weight=1.5)


def test_m4_known_baseline_refuses_c_down():
    with pytest.raises(ValueError):
        PeriodicVolatilityMonitor(baseline="known", tau=3.0, m0=M0, c_down=1.5)
    mon = PeriodicVolatilityMonitor(baseline="known", tau=3.0, m0=M0, tau_down=1.5, m0_down=float(trunc_mean_normal(1.5)))
    assert mon.tau2_dn == 2.25


def test_m5_bitstamp_filled_minutes_are_missing(tmp_path):
    from avmonitor.data import load_bitstamp_csv

    ts = 1_750_000_000 + 60 * np.arange(6)
    close = [100.0, 101.0, 101.0, 102.0, 101.5, 101.5]
    vol = [1.0, 2.0, 0.0, 3.0, 1.0, 0.5]  # minute 2 is a filled flat candle; minute 5 a real unchanged price
    pd.DataFrame(dict(timestamp=ts, open=close, high=close, low=close, close=close, volume=vol)).to_csv(
        tmp_path / "b.csv", index=False)
    df = load_bitstamp_csv(tmp_path / "b.csv")
    r = df["logret"].to_numpy()
    assert np.isnan(r[2]) and np.isnan(r[3])  # both sides of the filled minute
    assert r[5] == 0.0 and np.isfinite(r[1]) and np.isfinite(r[4])


def test_m6_two_sided_alarm_with_mixed_integer_types():
    a = D.DetectorResult(np.array([-1, 5], np.int64))
    b = D.DetectorResult(np.array([-1, 7], np.int32))
    assert list(D.two_sided_alarm(a, b)) == [-1, 5]


def test_m7_sr_steps_at_an_exact_boundary():
    for F in (1.01, 1.125, 2.0):
        for n in range(1, 40):
            A = F * (F**n - 1.0) / (F - 1.0)  # M_n of the docstring's closed form equals A exactly
            assert sr_steps_constant_factor(A, F) == n
            assert sr_steps_constant_factor(np.nextafter(A, np.inf), F) == n + 1
