"""E20 (experiments/e20_inference.py): the helpers of the re-analysis.

* Holm's step-down adjustment on a hand-computed example, with NaN entries left out of the family;
* the circular block bootstrap weights: every replicate has D days, each day is drawn once on average,
  and draws come in runs of L consecutive days (circularly);
* the percentile-inversion p-value at its extremes;
* Driscoll-Kraay standard errors: with one observation per day and bandwidth 0 they are the
  heteroskedasticity-robust sandwich times D / (D - 1); positive autocorrelation of the daily scores
  widens them;
* the QLIKE-optimal scale factor is the mean of y / F.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy import optimize

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def E():
    spec = importlib.util.spec_from_file_location("e20_inference", ROOT / "experiments" / "e20_inference.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_holm_example(E):
    p = np.array([0.01, 0.04, 0.03, 0.005])
    # sorted 0.005, 0.01, 0.03, 0.04 -> 4 x 0.005, 3 x 0.01, 2 x 0.03, max(0.06, 1 x 0.04)
    assert np.allclose(E.holm(p), [0.03, 0.06, 0.06, 0.02])
    q = E.holm([0.01, np.nan, 0.04])
    assert np.isnan(q[1]) and np.allclose(q[[0, 2]], [0.02, 0.04])
    assert np.all(E.holm([0.6, 0.7]) <= 1.0)


def test_cbb_weights(E):
    D, L, B = 23, 5, 4000
    W = E.cbb_weights(D, L, B, seed=1)
    assert W.shape == (B, D)
    assert np.all(W.sum(axis=1) == D)
    assert np.allclose(W.mean(axis=0), 1.0, atol=0.08)
    # with L = D every replicate is a rotation of all days: each day exactly once
    W2 = E.cbb_weights(10, 10, 50, seed=2)
    assert np.all(W2 == 1)
    # with D a multiple of L, a single block per replicate covers L consecutive days (circularly)
    W3 = E.cbb_weights(12, 12 // 2, 200, seed=3)
    for row in W3:
        assert row.sum() == 12


def test_boot_p(E):
    rng = np.random.default_rng(0)
    assert E.boot_p(rng.normal(0, 1, 9999)) > 0.9
    assert E.boot_p(np.full(9999, 0.3)) == pytest.approx(2 / 10000)
    assert E.boot_p(np.full(9999, -0.3)) == pytest.approx(2 / 10000)


def test_ols_dk_matches_sandwich(E):
    rng = np.random.default_rng(4)
    n = 400
    X = np.column_stack([np.ones(n), rng.normal(size=n)])
    y = X @ np.array([1.0, 2.0]) + rng.normal(size=n) * (1 + np.abs(X[:, 1]))
    beta, se, nd = E.ols_dk(y, X, np.arange(n), bandwidth=0)
    b_ls = np.linalg.lstsq(X, y, rcond=None)[0]
    assert np.allclose(beta, b_ls)
    e = y - X @ b_ls
    XtX_inv = np.linalg.inv(X.T @ X)
    V = XtX_inv @ (X.T * e ** 2) @ X @ XtX_inv * n / (n - 1)
    assert np.allclose(se, np.sqrt(np.diag(V)))
    assert nd == n


def test_ols_dk_widens_under_autocorrelation(E):
    rng = np.random.default_rng(5)
    n = 3000
    x = np.zeros(n)
    e = np.zeros(n)
    for t in range(1, n):  # persistent regressor and error
        x[t] = 0.9 * x[t - 1] + rng.normal()
        e[t] = 0.9 * e[t - 1] + rng.normal()
    X = np.column_stack([np.ones(n), x])
    y = X @ np.array([0.0, 1.0]) + e
    _, se0, _ = E.ols_dk(y, X, np.arange(n), bandwidth=0)
    _, se10, _ = E.ols_dk(y, X, np.arange(n), bandwidth=10)
    assert se10[1] > 1.5 * se0[1]
    # several observations on one day are summed before the kernel: duplicating every day's row
    # (perfectly dependent within the day) leaves the estimate and its standard error unchanged
    Xd, yd, dd = np.repeat(X, 2, axis=0), np.repeat(y, 2), np.repeat(np.arange(n), 2)
    b1, s1, _ = E.ols_dk(y, X, np.arange(n), bandwidth=3)
    b2, s2, _ = E.ols_dk(yd, Xd, dd, bandwidth=3)
    assert np.allclose(b1, b2) and np.allclose(s1, s2)


def test_qlike_kappa(E):
    rng = np.random.default_rng(6)
    F = np.exp(rng.normal(size=500))
    y = F * np.exp(rng.normal(0.3, 0.5, size=500))
    k = E.qlike_kappa(y, F)
    res = optimize.minimize_scalar(lambda c: np.mean(E.qlike(y, c * F)), bounds=(0.1, 10), method="bounded")
    assert k == pytest.approx(res.x, rel=1e-4)
    assert np.all(E.qlike(F, F) == 0)
