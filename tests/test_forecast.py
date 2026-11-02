import numpy as np
import pytest

from avmonitor import forecast as F


def _garch(n, omega, a, b, seed):
    rng = np.random.default_rng(seed)
    r = np.empty(n)
    s2 = omega / (1 - a - b)
    for t in range(n):
        r[t] = np.sqrt(s2) * rng.standard_normal()
        s2 = omega + a * r[t] ** 2 + b * s2
    return r


def test_forecasts_are_predictable():
    rng = np.random.default_rng(0)
    r = rng.standard_normal(500)
    r2 = r.copy()
    r2[300] = 50.0  # changing r_t must not change the forecast for t, only for t + 1
    for f in (lambda x: F.ewma_variance(x, 30, init=1.0), lambda x: F.garch11_filter(x, 0.05, 0.05, 0.9)):
        v, v2 = f(r), f(r2)
        assert np.array_equal(v[:301], v2[:301])
        assert v2[301] > v[301]


def test_ewma_tracks_a_constant_variance():
    r = np.random.default_rng(1).standard_normal(20000) * 2.0
    v = F.ewma_variance(r, 200)
    assert np.median(v[5000:]) == pytest.approx(4.0, rel=0.05)


def test_nan_returns_carry_the_forecast_forward():
    r = np.random.default_rng(2).standard_normal(300)
    r[100:120] = np.nan
    v = F.ewma_variance(r, 20, init=1.0)
    assert np.all(v[101:121] == v[100])
    z = F.standardise(r, v)
    assert np.all(np.isnan(z[100:120])) and np.all(np.isfinite(z[:100]))


def test_garch_fit_recovers_parameters():
    r = _garch(30000, 0.05, 0.05, 0.90, seed=3)
    omega, a, b = F.garch11_fit(r)
    assert a == pytest.approx(0.05, abs=0.02)
    assert b == pytest.approx(0.90, abs=0.04)
    assert omega / (1 - a - b) == pytest.approx(1.0, rel=0.1)


def test_oracle_standardisation_removes_clustering():
    r = _garch(50000, 0.01, 0.05, 0.94, seed=4)
    v = F.garch11_filter(r, 0.01, 0.05, 0.94)
    z = F.standardise(r, v)
    # squared standardised returns have (almost) no autocorrelation, squared raw returns have a lot
    def ac1(x):
        x = x - x.mean()
        return float(np.dot(x[1:], x[:-1]) / np.dot(x, x))
    assert ac1(r**2) > 0.1
    assert abs(ac1(z**2)) < 0.02
