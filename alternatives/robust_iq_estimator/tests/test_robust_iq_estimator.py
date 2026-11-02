"""Tests for the robust local plug-in IQ estimators (run from the repository root with pytest).

All tests are synthetic, need no data files and run in a few seconds.
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import integrate, stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import riq  # noqa: E402
from avmonitor.robust import realized_quarticity, truncated_moment_factor  # noqa: E402


def test_truncated_moments_match_numerical_integration_and_repo_kappa():
    for a in (2.0, 3.0, 4.5):
        T = riq.trunc_even_moments(a, 4)
        for q in range(5):
            num = integrate.quad(lambda z: z ** (2 * q) * stats.norm.pdf(z), -a, a)[0]
            assert abs(T[q] - num) < 1e-9
        clip = riq.square_moments(a, "clip")
        assert abs(clip[0] - truncated_moment_factor(a, 2)) < 1e-12
        assert abs(clip[1] / 3 - truncated_moment_factor(a, 4)) < 1e-12
    assert np.allclose(riq.square_moments(np.inf), [1, 3, 15, 105])


def test_mean_power_moments_reduce_to_chi_square_for_gaussian():
    m = np.array([1.0, 2.0, 7.0, 40.0])
    e2, e4 = riq.mean_power_moments(m, riq.square_moments(np.inf))
    assert np.allclose(e2, (m + 2) / m)
    assert np.allclose(e4, (m + 2) * (m + 4) * (m + 6) / m**3)
    assert np.allclose(m * (e4 / e2**2 - 1), m * riq.gaussian_block_relvar(m))


def test_avar_constants():
    assert riq.avar_block(1) == pytest.approx(32 / 3)
    assert riq.avar_block_exact(1) == pytest.approx(32 / 3)
    assert riq.avar_block_exact(20) == pytest.approx(riq.avar_block(20))
    assert riq.avar_block(np.inf) == 8.0
    assert riq.avar_truncated_rq(np.inf) == pytest.approx(32 / 3)
    assert riq.avar_tpq() == pytest.approx(13.6495, abs=1e-3)
    assert 14.9 < riq.avar_qpq() < 14.95
    # an oracle cut-off at a known sigma can go below 8; the estimated one cannot (see README 2.3)
    assert riq.avar_block_truncated(3.0) < 8.0


def test_block_plugin_is_unbiased_with_correction_and_biased_without():
    rng = np.random.default_rng(0)
    n, reps = 200, 4000
    r = rng.standard_normal((reps, n)) * 2.0 / np.sqrt(n)  # sigma = 2, IQ = 16
    out = riq.block_plugin(r, k=10)
    assert abs(out["iq"].mean() / 16 - 1) < 0.01
    assert abs(out["iq8"].mean() / 256 - 1) < 0.03
    # n Var / sigma^8 matches 8 (k + 3) / (k + 2) within simulation error (about 2%)
    assert abs(n * out["iq"].var() / 256 / riq.avar_block(10) - 1) < 0.08
    assert abs(out["se"].mean() / out["iq"].std() - 1) < 0.08
    raw = riq.block_plugin(r, k=10, correct=False)
    assert abs(raw["iq"].mean() / 16 - 1.2) < 0.015  # (k + 2) / k
    # the k = 1 plug-in is realized quarticity
    assert np.allclose(riq.block_plugin(r[:3], k=1)["iq"], realized_quarticity(r[:3]))


def test_truncation_correction_exact_with_oracle_threshold():
    rng = np.random.default_rng(1)
    n, reps = 400, 3000
    r = rng.standard_normal((reps, n)) / np.sqrt(n)
    for mode in ("clip", "drop"):
        tau = np.full(r.shape, 3.0 / np.sqrt(n))
        out = riq.block_plugin(r, tau, 3.0, k=20, mode=mode)
        assert abs(out["iq"].mean() - 1) < 0.012, mode


def test_clipped_jump_residual_matches_simulation_and_drop_leaves_none():
    rng = np.random.default_rng(2)
    n, reps, k, a = 400, 6000, 20, 3.0
    z = rng.standard_normal((reps, n)) / np.sqrt(n)
    tau = np.full(z.shape, a / np.sqrt(n))
    zj = z.copy()
    zj[:, 205] += 40.0 / np.sqrt(n)
    d = riq.block_plugin(zj, tau, a, k=k, mode="clip")["iq"] - riq.block_plugin(z, tau, a, k=k, mode="clip")["iq"]
    assert d.mean() * n == pytest.approx(riq.clipped_jump_residual(a, k), rel=0.03)
    d4 = riq.truncated_rq_local(zj, tau, a, mode="clip") - riq.truncated_rq_local(z, tau, a, mode="clip")
    assert d4.mean() * n == pytest.approx(riq.clipped_jump_residual(a, 1), rel=0.03)
    dd = riq.block_plugin(zj, tau, a, k=k, mode="drop")["iq"] - riq.block_plugin(z, tau, a, k=k, mode="drop")["iq"]
    assert abs(dd.mean() * n) < 1.0  # against 17.5 for clipping
    # undetected jump: the plug-in damps the b^4 term
    assert riq.undetected_jump_residual(3.0, 20) < riq.undetected_jump_residual(3.0, 1)


def test_local_scale_follows_a_volatility_step_and_ignores_a_jump():
    rng = np.random.default_rng(3)
    n = 1600
    sd = np.where(np.arange(n) < n // 2, 1.0, 3.0)
    r = rng.standard_normal((200, n)) * sd
    r[:, 300] += 200.0
    s = riq.local_scale(r)
    assert abs(np.median(s[:, 100:600]) - 1.0) < 0.05
    assert abs(np.median(s[:, 1000:1500]) - 3.0) < 0.15
    assert np.median(s[:, 300]) < 1.15


def test_kernel_plugin_and_local_threshold_estimators_are_unbiased_for_constant_vol():
    rng = np.random.default_rng(4)
    n, reps = 600, 1500
    r = rng.standard_normal((reps, n)) / np.sqrt(n)
    assert abs(riq.kernel_plugin(r, h=30).mean() - 1) < 0.015
    tau, a = riq.thresholds(r, rule="fixed", c=3.0)
    assert abs(riq.kernel_plugin(r, tau, a, h=30).mean() - 1) < 0.02
    tau, a = riq.thresholds(r, rule="universal")
    assert a == pytest.approx(np.sqrt(2 * np.log(n)))
    assert abs(riq.block_plugin(r, tau, a, mode="drop")["iq"].mean() - 1) < 0.02


def test_printed_cv_prefers_small_bandwidth_for_scale_reasons():
    printed = [riq.printed_cv_expectation(h, 3.0)[0] for h in range(2, 80)]
    corrected = [riq.printed_cv_expectation(h, 3.0)[1] for h in range(2, 80)]
    h_star = 2 + int(np.argmin(printed))
    assert h_star <= 10  # below the paper's whole grid {10, ..., 50}
    grid = [riq.printed_cv_expectation(h, 3.0)[0] for h in (10, 20, 30, 40, 50)]
    assert np.all(np.diff(grid) > 0)  # so eq. (17) as printed always picks the smallest h of the grid
    assert np.all(np.diff(corrected) <= 1e-12)  # constant volatility: larger h is never worse


def test_cv_curves_printed_matches_direct_formula():
    rng = np.random.default_rng(5)
    x = rng.standard_normal(60)
    out = riq.cv_curves(x[None, :], (5, 10))
    from avmonitor.robust import kernel_weights

    h = 5
    v = kernel_weights(h, normalise=False)[:-1]
    direct = sum((x[t] ** 4 - sum(v[i] * x[t - h + i] ** 4 for i in range(h)) / 3) ** 2 for t in range(h, 60)) / 60
    assert out["printed"][0, 0] == pytest.approx(direct)


def test_nan_and_one_dimensional_input():
    rng = np.random.default_rng(6)
    r = rng.standard_normal(500) / np.sqrt(500)
    r[[10, 11, 250]] = np.nan
    tau, a = riq.thresholds(r, rule="universal")
    out = riq.block_plugin(r, tau, a)
    assert isinstance(out["iq"], float) and np.isfinite(out["iq"]) and out["se"] > 0
    assert np.isfinite(riq.truncated_rq_local(r, tau, a))
    clean = r[np.isfinite(r)]
    assert np.isfinite(realized_quarticity(clean))
