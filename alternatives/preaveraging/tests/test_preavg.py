"""Tests for the pre-averaging prototype (run from the repository root with pytest)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import preavg as P  # noqa: E402
from avmonitor.detectors import betting_edetector  # noqa: E402


def _noisy_paths(reps, n, sigma, gamma, rng, jump=0.0):
    """Constant-volatility Brownian returns with price-level noise of ratio ``gamma``."""
    s = sigma / np.sqrt(n)
    x = rng.standard_normal((reps, n)) * s
    if jump > 0:
        x[:, n // 2] += jump * s
    u = rng.standard_normal((reps, n + 1)) * gamma * s
    return x + np.diff(u, axis=1), s


def test_psi_constants_match_the_continuous_limits():
    psi1, psi2 = P.psi_constants(60)
    assert abs(psi1 - 1.0) < 1e-12
    assert abs(psi2 - 1 / 12) < 1e-3
    assert P.choose_k(23400, 0.4) % 2 == 0


def test_preaveraged_second_moment_formula():
    rng = np.random.default_rng(0)
    n, k, gamma = 20000, 20, 1.5
    r, s = _noisy_paths(20, n, 1.0, gamma, rng)
    rbar = P.preaverage(r, k)
    psi1, psi2 = P.psi_constants(k)
    expected = psi2 * k * s**2 + psi1 * (gamma * s) ** 2 / k
    assert abs(np.mean(rbar**2) / expected - 1) < 0.03
    # blocks are every k-th window and disjoint in price points
    blocks = P.preaverage(r, k, step=k)
    assert blocks.shape[-1] == len(range(0, n - k + 1, k))
    assert np.allclose(blocks, rbar[:, ::k])


def test_noise_variance_estimator_and_noise_corrected_iv_iq():
    rng = np.random.default_rng(1)
    n, gamma = 23400, 2.0
    r, s = _noisy_paths(30, n, 1.0, gamma, rng)
    om = P.noise_variance(r)
    assert abs(om.mean() / (gamma * s) ** 2 - 1) < 0.05
    est = P.preavg_estimates(r, 60, T=1.0, c=3.0)
    assert abs(est["iv"].mean() - 1.0) < 0.05  # true IV = 1
    assert abs(est["iq"].mean() - 1.0) < 0.08  # true IQ = 1
    raw = P.raw_truncated_estimates(r, 1.0, 3.0)
    assert raw["iq"].mean() > 50  # (1 + 2 gamma^2)^2 = 81 without pre-averaging


def test_truncation_factors_and_drop_mode_remove_a_jump():
    f2c, f4c = P.truncation_factors(3.0, "clip")
    f2d, f4d = P.truncation_factors(3.0, "drop")
    assert 0.99 < f2c < 1.0 and 0.95 < f4c < 0.97 and f2d < f2c and f4d < f4c
    rng = np.random.default_rng(2)
    r, s = _noisy_paths(30, 23400, 1.0, 1.0, rng, jump=80.0)
    clip = P.preavg_estimates(r, 60, c=3.0, mode="clip")["iq"].mean()
    drop = P.preavg_estimates(r, 60, c=3.0, mode="drop")["iq"].mean()
    plain = P.preavg_estimates(r, 60, c=None)["iq"].mean()
    assert plain > 2.0  # an 80 sd jump dominates the fourth power
    assert abs(drop - 1.0) < 0.1
    assert drop < clip < plain


def test_spot_variance_raw_is_noise_inflated_and_preavg_is_not():
    rng = np.random.default_rng(3)
    n, gamma = 6000, 1.0
    r, s = _noisy_paths(10, n, 1.0, gamma, rng)
    delta = 1.0 / n
    vr = P.spot_variance_raw(r, delta, 300)
    vp = P.spot_variance_preavg(r, delta, 300, 20)
    assert abs(np.nanmean(vr) / 3.0 - 1) < 0.1  # 1 + 2 gamma^2 = 3
    assert abs(np.nanmean(vp) - 1.0) < 0.15


def test_block_stream_keeps_false_alarm_bound_and_overlap_does_not_matter_for_the_mean():
    rng = np.random.default_rng(4)
    n, k = 6000, 20
    # calibration: ten paths with the same per-step scale, joined end to end
    cal = _noisy_paths(10, n, 1.0, 1.0, rng)[0].reshape(1, -1)
    r, _ = _noisy_paths(400, n, 1.0, 1.0, rng)
    tau2 = P.monitor_stream(cal, k, k)[1]
    xb, _ = P.monitor_stream(r, k, k, tau=np.sqrt(tau2))
    m0 = float(P.monitor_stream(cal, k, k, tau=np.sqrt(tau2))[0].mean())
    A = 500.0
    res = betting_edetector(xb, m0, "up", threshold=A)
    m = xb.shape[1]
    assert (res.alarm >= 0).mean() <= m / A + 0.03


def test_mbv_block_variance_matches_exact_weights_and_estimates_are_unbiased():
    rng = np.random.default_rng(5)
    n, gamma = 23400, 1.0
    r, s = _noisy_paths(60, n, 1.0, gamma, rng)
    e = P.mbv_estimates(r, c1=0.5, c2=1.0)
    yb = P.mbv_block_averages(r, e["K"], e["block"])
    expected = e["sig_w"] * s**2 + e["noise_w"] * (gamma * s) ** 2
    assert abs(np.mean(yb**2) / expected - 1) < 0.04
    # the exact weights agree with the asymptotic constants to a few percent
    assert abs(e["sig_w"] / (e["nu1"] * np.sqrt(n)) - 1) < 0.03
    assert abs(e["noise_w"] * np.sqrt(n) / e["nu2"] - 1) < 0.03
    assert abs(e["iv"].mean() - 1.0) < 0.06
    assert abs(e["iq"].mean() - 1.0) < 0.15
    cover = np.mean(np.abs(e["iv"] - 1.0) <= 1.96 * e["iv_se"])
    assert cover > 0.8


def test_multipower_estimators_keep_a_jump_bias_that_dropping_removes():
    rng = np.random.default_rng(6)
    n, gamma = 23400, 2.0
    r, s = _noisy_paths(40, n, 1.0, gamma, rng, jump=100.0)
    mp = P.preavg_multipower_estimates(r, 60)
    mbv = P.mbv_estimates(r, 0.5, 1.0)
    drop = P.preavg_estimates(r, 60, c=3.0, mode="drop")
    plain = P.preavg_estimates(r, 60, c=None)
    assert abs(mp["iv"].mean() - 1.0) < 0.1  # bipower is nearly jump robust
    assert abs(drop["iq"].mean() - 1.0) < 0.1
    # a 100 sd jump keeps |J|^{4/3} in every product it enters: a visible but bounded bias
    assert 1.0 < mp["iq"].mean() < plain["iq"].mean()
    assert 1.0 < mbv["iq"].mean() < plain["iq"].mean()
    # k* = 6 gamma / v (delay formula), floored at 4: 4 gamma for the default target v = 1.5
    assert P.rule_of_thumb_k(0.0) == 4 and P.rule_of_thumb_k(2.0) == 8 and P.rule_of_thumb_k(4.0, 2.0) == 12
    assert abs(P.block_noise_share(60, 2.0) - 12 * 4 / 3600) < 2e-3
