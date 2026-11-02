"""Edge cases found while extending the prototype to real data (iteration 3 of ITERATION_LOG.md)."""

import warnings

import numpy as np
import pytest

from avmonitor import robust as R


def _zero_heavy(n=2000, share=0.6, seed=0):
    rng = np.random.default_rng(seed)
    r = rng.standard_normal(n) * 1e-3
    r[rng.random(n) < share] = 0.0
    return r


def test_zero_mad_warns_instead_of_silently_returning_zero():
    r = _zero_heavy()
    with pytest.warns(R.DegenerateThresholdWarning):
        tau = R.mad_threshold(r, 3.0, scaled=True)
    assert tau == 0.0
    with pytest.warns(R.DegenerateThresholdWarning):
        assert R.truncated_rq_mad(r) == 0.0  # the silent failure the warning exists for


def test_floor_fixes_zero_mad():
    r = _zero_heavy()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        tau = R.mad_threshold(r, 3.0, scaled=True, floor=3e-3)
    assert tau == pytest.approx(3e-3)


def test_no_warning_on_ordinary_or_all_zero_data():
    rng = np.random.default_rng(1)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        R.mad_threshold(rng.standard_normal(500), 3.0)
        R.mad_threshold(np.zeros(500), 3.0)  # nothing to erase


def test_tick_floor_recovers_the_tick():
    rng = np.random.default_rng(2)
    ticks = np.cumsum(rng.integers(-1, 2, 5000))
    prices = 0.2 + 1e-4 * (ticks - ticks.min())  # quoted to 4 decimals near 0.2
    f = R.tick_floor(prices, c=3.0)
    assert f == pytest.approx(3 * np.log1p(1e-4 / np.median(prices)), rel=1e-6)


def test_ew_scale_does_not_get_stuck_at_zero():
    """With more than half zeros the true median of |r| is 0, so the tracker drifts down.

    It must never underflow to exactly 0 (it could not recover), and with a
    tick-size floor it must recover within a few thousand steps once ordinary
    returns resume. Without a floor recovery is slow (at most a factor e per
    1 / (eta / 2) = 562 steps for half-life 390), which is a property of the
    multiplicative update, documented here.
    """
    rng = np.random.default_rng(3)
    r = np.concatenate([_zero_heavy(20000, 0.6, 4), rng.standard_normal(4000) * 1e-3])
    s = R.ew_quantile_scale(r, halflife=390, scaled=True)
    assert np.all(s > 0) and np.all(np.isfinite(s))
    s_floor = R.ew_quantile_scale(r, halflife=390, scaled=True, floor=1e-4)
    assert np.all(s_floor >= 1e-4 * R.MAD_NORMAL_SCALE * (1 - 1e-12))
    assert s_floor[-1] == pytest.approx(1e-3, rel=0.25)


def test_ew_scale_ignores_nan():
    rng = np.random.default_rng(5)
    r = rng.standard_normal(3000)
    r_nan = r.copy()
    r_nan[1000:1500] = np.nan
    s, s_nan = R.ew_quantile_scale(r), R.ew_quantile_scale(r_nan)
    # during the gap the estimate is frozen, not dragged toward zero
    assert np.all(s_nan[1001:1501] == s_nan[1000])
    assert s_nan[1500] == pytest.approx(s[1000], rel=1e-12)


def test_ew_scale_unchanged_on_clean_data():
    """Regression: the rewrite gives the same numbers as the original update rule."""
    rng = np.random.default_rng(6)
    a = np.abs(rng.standard_normal(5000))
    eta = np.log(2) / 390 * 2.0
    m = float(np.median(a[:50]))
    ref = np.empty(a.size)
    for t in range(a.size):
        ref[t] = m
        m *= np.exp(eta * ((a[t] > m) - 0.5))
    assert np.allclose(R.ew_quantile_scale(a, scaled=False), ref, rtol=1e-12, atol=0)
