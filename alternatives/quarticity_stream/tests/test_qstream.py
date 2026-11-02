"""Tests for the quarticity-stream prototype (run from the repository root with pytest)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import qstream as Q  # noqa: E402
from avmonitor.robust import truncated_moment_factor  # noqa: E402


def test_quarticity_stream_is_square_of_variance_stream_and_bounded():
    rng = np.random.default_rng(0)
    r = rng.standard_normal(1000) * 2.0
    x, q = Q.variance_stream(r, 3.0), Q.quarticity_stream(r, 3.0)
    assert np.allclose(q, x**2)
    assert q.max() <= 81.0 + 1e-12 and q.min() >= 0.0


def test_calibration_matches_gaussian_closed_form_and_q0_is_not_3m0_squared():
    rng = np.random.default_rng(1)
    r = rng.standard_normal(400_000)
    b = Q.calibrate_streams(r, c=3.0)
    m0, q0 = Q.gaussian_null_moments(1.0, 3.0)
    assert abs(b.m0 / m0 - 1) < 0.01
    assert abs(b.q0 / q0 - 1) < 0.02
    # the Gaussian plug-in 3 m0^2 overstates q0 by about 2.7% at a = 3
    assert 3 * m0**2 / q0 > 1.02
    # bounds and admissible decrease bets
    assert abs(b.b_q - 81.0 / q0) / b.b_q < 0.03
    assert b.lam_max_down_q < b.lam_max_down_var


def test_stream_shift_is_neutral_under_null_and_signs_under_alternatives():
    null = Q.stream_shift(Q.scaled_normal_pdf(1.0), 3.0)
    assert abs(null["mean_var"] - 1) < 1e-6 and abs(null["mean_q"] - 1) < 1e-6
    up = Q.stream_shift(Q.scaled_normal_pdf(1.5), 3.0)
    assert up["mean_var"] > 1.5 and up["mean_q"] > up["mean_var"]
    mix = Q.stream_shift(Q.mixture_pdf(0.7), 3.0)
    assert mix["mean_q"] > 1.2 and mix["mean_var"] < 1.0
    # Gaussian kappa_4 check: q0 = 3 kappa_4(a) for sigma = 1
    assert abs(Q._clipped_moment(Q.scaled_normal_pdf(1.0), 4, 3.0) / (3 * truncated_moment_factor(3.0, 4)) - 1) < 1e-6


def test_null_false_alarm_bound_holds_for_quarticity_stream_both_sides():
    rng = np.random.default_rng(2)
    r = Q.shape_change_returns(1500, 60, None, "sd", 1.0, rng)
    base = Q.calibrate_streams(rng.standard_normal(200_000), c=3.0)
    A = 300.0
    for direction in ("up", "down"):
        res = Q.run_stream_detector(r, base, "quarticity", direction, threshold=A)
        # P(alarm within 60 steps) <= 60 / A = 0.2, plus Monte Carlo slack
        assert (res.alarm >= 0).mean() <= 60 / A + 0.03


def test_detects_variance_change_and_vol_of_vol_change():
    rng = np.random.default_rng(3)
    base = Q.calibrate_streams(rng.standard_normal(100_000), c=3.0)
    up = Q.shape_change_returns(60, 1500, 300, "sd", 2.0, rng)
    res = Q.run_stream_detector(up, base, "quarticity", "up", threshold=5000.0)
    det = (res.alarm >= 300).mean()
    assert det > 0.9 and np.median(res.alarm[res.alarm >= 300] - 300) < 150
    mix = Q.shape_change_returns(40, 6000, 300, "mixture", 0.9, rng)
    rq = Q.run_stream_detector(mix, base, "quarticity", "up", threshold=2000.0)
    rv = Q.run_stream_detector(mix, base, "variance", "up", threshold=2000.0)
    assert (rq.alarm >= 300).mean() > 0.5
    assert (rq.alarm >= 300).mean() > (rv.alarm >= 300).mean()


def test_decrease_bet_range_is_enforced():
    rng = np.random.default_rng(4)
    base = Q.calibrate_streams(rng.standard_normal(50_000), c=3.0)
    r = rng.standard_normal((2, 50))
    with pytest.raises(ValueError):
        Q.run_stream_detector(r, base, "quarticity", "down", threshold=10.0, grid=[1.5])
    res = Q.run_stream_detector(r, base, "quarticity", "down", threshold=10.0)
    assert np.all(res.info["grid"] < base.lam_max_down_q + 1e-12)


def test_streaming_monitor_recalibrates_and_alarms_on_vol_of_vol():
    rng = np.random.default_rng(5)
    calm = rng.standard_normal(3000)
    wild = Q.shape_change_returns(1, 6000, 0, "mixture", 0.9, rng)[0]
    mon = Q.QuarticityMonitor(threshold=2000.0, cal_len=2000)
    ev = mon.run(np.concatenate([calm, wild]))
    ups = [e for e in ev if e.direction == "up" and e.index >= 3000]
    assert ups
