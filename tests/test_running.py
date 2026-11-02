"""Tests for ``avmonitor.running`` (route A: uniform bands with a running truncation level; route B: ranks).

All tests use synthetic data and run without data/raw or network.
"""

import math

import numpy as np
import pytest
from scipy import stats

from avmonitor import running as R
from avmonitor.detectors import betting_edetector


# ----------------------------------------------------------------------------
# Route A: the bands
# ----------------------------------------------------------------------------


def test_hr_constant_matches_the_closed_form_and_is_monotone():
    """C(alpha) from the authors' optimisation lies just below the closed form 0.8 log(1612 / alpha)."""
    prev = math.inf
    for a in (0.001, 0.01, 0.05, 0.1):
        C = R.hr_dkw_constant(a)
        assert 0.9 * 0.8 * math.log(1612 / a) < C <= 0.8 * math.log(1612 / a)
        assert C < prev
        prev = C
    eps = R.hr_dkw_radius(np.array([5.0, 10.0, 100.0]), 0.05, t_min=10)
    assert np.isinf(eps[0]) and eps[1] > eps[2] > 0
    # at t = t_min the radius is A sqrt(C / t)
    assert eps[1] == pytest.approx(0.85 * math.sqrt(R.hr_dkw_constant(0.05) / 10))


def test_hr_band_is_time_uniform_on_uniform_streams():
    """P(sup_x |F_t - F| > eps_t for some checked t) <= alpha (checked on a log-spaced set of times)."""
    rng = np.random.default_rng(1)
    alpha, reps, n = 0.1, 300, 3000
    times = np.unique(np.geomspace(5, n, 40).astype(int))
    eps = R.hr_dkw_radius(times, alpha)
    u = rng.random((reps, n))
    crossed = np.zeros(reps, bool)
    for t, e in zip(times, eps):
        s = np.sort(u[:, :t], axis=1)
        i = np.arange(1, t + 1)
        d = np.maximum(i / t - s, s - (i - 1) / t).max(axis=1)
        crossed |= d > e
    rate = crossed.mean()
    assert rate <= alpha + 3 * math.sqrt(alpha * (1 - alpha) / reps)


def test_band_integrals_match_brute_force():
    rng = np.random.default_rng(2)
    v = np.sort(rng.standard_t(3, 200) ** 2)
    a = np.array([0.05, 0.5, 1.0, 4.0, 9.0, 50.0])
    for eps in (0.0, 0.03, 0.2):
        up, lo, mh = R.truncated_means_from_band(v, a, eps)
        np.testing.assert_allclose(mh, [np.mean(np.minimum(v, x)) for x in a], rtol=1e-12)
        for k, x in enumerate(a):
            grid = np.linspace(0, x, 200001)
            mid = 0.5 * (grid[1:] + grid[:-1])
            S = 1.0 - np.searchsorted(v, mid, side="right") / v.size
            h = grid[1] - grid[0]
            assert up[k] == pytest.approx(np.sum(np.minimum(1.0, S + eps)) * h, rel=2e-3, abs=1e-6)
            assert lo[k] == pytest.approx(np.sum(np.maximum(0.0, S - eps)) * h, rel=2e-3, abs=1e-6)
        assert np.all(up >= mh - 1e-12) and np.all(lo <= mh + 1e-12)


def _true_m(levels, draw, n=1_000_000, seed=0):
    big = draw(np.random.default_rng(seed), n) ** 2
    return np.array([np.mean(np.minimum(big, s * s)) for s in np.ravel(levels)]).reshape(np.shape(levels))


@pytest.mark.parametrize("method", ["grid", "dkw"])
def test_band_covers_the_truncated_mean_at_every_level(method):
    """With delta = 0.1 the band misses m_F(s_k) at some level on at most about 10% of windows."""
    rng = np.random.default_rng(3)
    reps, n = 200, 600
    r = rng.standard_normal((reps, n))
    b = R.calibrate_band(r, method=method, delta=0.1)
    big = np.sort(np.random.default_rng(99).standard_normal(1_000_000) ** 2)
    csum = np.concatenate([[0.0], np.cumsum(big)])
    s2 = b.levels**2
    k = np.searchsorted(big, s2, side="right")
    truth = (csum[k] + s2 * (big.size - k)) / big.size
    miss = np.any((truth > b.upper) | (truth < b.lower), axis=1)
    assert miss.mean() <= 0.1 + 3 * math.sqrt(0.09 / reps)
    assert np.all(b.upper >= b.lower)
    # the grid band is much tighter than the DKW band at the centre level
    if method == "grid":
        mid = -R.GRID_KMIN
        assert np.median(b.upper[:, mid] / truth[:, mid]) < 1.6


def test_running_level_is_predictable():
    rng = np.random.default_rng(4)
    r = rng.standard_normal((5, 800))
    center = np.full(5, 3.0)
    k1 = R.running_level_index(r, 400, 400, center)
    r2 = r.copy()
    r2[:, 600:] *= 5.0  # change the future after monitored step 200
    k2 = R.running_level_index(r2, 400, 400, center)
    np.testing.assert_array_equal(k1[:, :201], k2[:, :201])
    assert np.any(k1[:, 260:] != k2[:, 260:])
    assert k1.min() >= 0 and k1.max() <= R.GRID_KMAX - R.GRID_KMIN


def test_two_sided_sr_with_fixed_inputs_matches_betting_edetector():
    rng = np.random.default_rng(5)
    m0, tau2, A = 0.99, 9.0, 300.0
    x_up = np.minimum(rng.standard_normal((50, 3000)) ** 2 * 1.4, tau2)
    x_dn = np.minimum(rng.standard_normal((50, 3000)) ** 2 * 0.6, tau2)
    up = betting_edetector(x_up, m0, "up", threshold=A).alarm
    dn = betting_edetector(x_dn, m0, "down", upper=tau2, threshold=A).alarm
    np.testing.assert_array_equal(R.two_sided_sr(x_up, tau2, m0, m0, [A])[A][0], up)
    np.testing.assert_array_equal(R.two_sided_sr(x_dn, tau2, m0, m0, [A])[A][1], dn)
    assert np.all(up >= 0) and np.all(dn >= 0)


def test_route_a_false_alarms_respect_delta_plus_m_over_A():
    """i.i.d. Student t3 returns, running tau, small A: P(alarm within m) <= delta + m / A (two-sided, 2A per side)."""
    rng = np.random.default_rng(7)
    reps, cal, m, A, delta = 600, 400, 200, 2000.0, 0.05
    r = rng.standard_t(3, (reps, cal + m))
    for method in ("grid", "dkw"):
        b = R.calibrate_band(r[:, :cal], method=method, delta=delta)
        tau2, mu, md, _ = R.route_a_inputs(b, r, cal, m)
        x = np.minimum(r[:, cal:] ** 2, tau2)
        au, ad = R.two_sided_sr(x, tau2, mu, md, [2 * A])[2 * A]
        rate = np.mean((au >= 0) | (ad >= 0))
        bound = delta + m / A
        assert rate <= bound + 3 * math.sqrt(bound * (1 - bound) / reps)


# ----------------------------------------------------------------------------
# Route B: ranks
# ----------------------------------------------------------------------------


def _brute_pvalues(a, theta):
    p = np.empty_like(a)
    for t in range(a.size):
        prev = a[: t + 1]
        p[t] = (np.sum(prev > a[t]) + theta[t] * np.sum(prev == a[t])) / (t + 1)
    return p


def test_conformal_pvalues_match_brute_force_with_ties():
    a = np.abs(np.round(np.random.default_rng(8).standard_normal((3, 60)) * 2) / 2)  # many ties
    p = R.conformal_pvalues(a, np.random.default_rng(9))
    theta = np.random.default_rng(9).random(a.shape)
    for i in range(3):
        np.testing.assert_allclose(p[i], _brute_pvalues(a[i], theta[i]), rtol=0, atol=1e-12)


@pytest.mark.parametrize("law", ["t3", "ticks"])
def test_conformal_pvalues_are_uniform_and_independent(law):
    rng = np.random.default_rng(10)
    if law == "t3":
        a = np.abs(rng.standard_t(3, (400, 300)))
    else:
        a = np.abs(np.round(rng.standard_normal((400, 300)) / 0.5))  # about 38% ties at zero
    p = R.conformal_pvalues(a, rng)
    assert stats.kstest(p[:, -1], "uniform").pvalue > 0.001
    assert stats.kstest(p[:, 5], "uniform").pvalue > 0.001
    assert stats.kstest(p.ravel()[::7], "uniform").pvalue > 0.001
    c = np.corrcoef(p[:, :-1].ravel(), p[:, 1:].ravel())[0, 1]
    assert abs(c) < 0.02


def test_betting_functions_have_mean_one():
    u = np.concatenate([np.geomspace(1e-9, 1e-2, 200001), np.linspace(1e-2, 1, 400001)[1:]])
    for f in R.bet_family("lr")[0] + R.bet_family("lr")[1]:
        val = np.trapezoid(f(u), u) + 1e-9 * f(1e-9)
        assert val == pytest.approx(1.0, abs=2e-4)
        assert np.all(f(u) > 0)
    for f in R.bet_family("ind")[0] + R.bet_family("ind")[1]:
        g = np.linspace(0, 1, 1_000_001)
        assert np.mean(f((g[1:] + g[:-1]) / 2)) == pytest.approx(1.0, abs=1e-5)
        assert np.all(f(g) >= 0)


def test_lr_bet_is_bounded_by_the_clip():
    f = R.lr_bet_factory(2.0)
    assert f(1e-12) == pytest.approx(f(f.u0))
    assert f(1e-12) == pytest.approx(0.5 * math.exp(9.0 * 0.75 / 2) / f.norm_const, rel=1e-6)  # g(u0) / C, about 16


@pytest.mark.parametrize("family", ["lr", "ind"])
def test_route_b_exact_false_alarm_bound(family):
    """Exchangeable data with ties and jumps: P(alarm within m) <= m / A with no estimated baseline."""
    rng = np.random.default_rng(11)
    reps, warm, m, A = 2000, 100, 100, 400.0
    r = np.round(rng.standard_t(3, (reps, warm + m)) * 4) / 4
    out, _ = R.conformal_detector(np.abs(r), rng, warm, [2 * A], family=family)
    au, ad = out[2 * A]
    rate = np.mean((au >= 0) | (ad >= 0))
    bound = m / A
    assert rate <= bound + 3 * math.sqrt(bound * (1 - bound) / reps)


def test_route_b_detects_a_doubling_quickly():
    rng = np.random.default_rng(12)
    r = rng.standard_normal((200, 1200))
    r[:, 700:] *= 2.0
    out, _ = R.conformal_detector(np.abs(r), rng, 500, [15600.0])
    au = out[15600.0][0]
    d = au[au >= 700] - 700
    assert np.mean(au >= 700) > 0.95 and np.median(d) < 40


def test_conformal_detector_agrees_with_conformal_pvalues():
    """With warm-up 0 the detector bets on exactly the p-values of ``conformal_pvalues`` (same tie-breakers)."""
    a = np.abs(np.round(np.random.default_rng(13).standard_normal((40, 400)) * 3) / 3)
    p = R.conformal_pvalues(a, np.random.default_rng(14))
    L = 30.0
    out, _ = R.conformal_detector(a, np.random.default_rng(14), 0, [L])
    up, dn = R.bet_family("lr")
    for side, fam in ((0, up), (1, dn)):
        S = np.zeros((40, len(fam)))
        first = np.full(40, -1)
        for t in range(400):
            S = np.column_stack([f(p[:, t]) for f in fam]) * (S + 1.0)
            first[(first < 0) & (S.mean(axis=1) >= L)] = t
        np.testing.assert_array_equal(out[L][side], first)


def test_streaming_conformal_pvalue_matches_brute_force():
    mon = R.ConformalMonitor(cal_len=5)
    vals = [0.5, 0.1, 0.5, 2.0, 0.3]
    for v in vals:
        mon.update(v)
    assert mon.calibrated and mon.ref == sorted(vals)
    mon.rng = np.random.default_rng(0)
    p = mon.pvalue(0.5)
    th = np.random.default_rng(0).random()
    assert p == pytest.approx((1 + th * 3) / 6)


def test_streaming_route_a_matches_the_vectorised_version():
    rng = np.random.default_rng(15)
    cal, n = 600, 2500
    r = rng.standard_normal(cal + n)
    r[cal + 800 :] *= 1.6
    A = 500.0
    mon = R.RunningTauMonitor(threshold=A, cal_len=cal, method="grid", delta=0.01)
    ev = mon.run(r)
    b = R.calibrate_band(r[None, :cal], method="grid", delta=0.01)
    tau2, mu, md, _ = R.route_a_inputs(b, r[None, :], cal, n)
    x = np.minimum(r[None, cal:] ** 2, tau2)
    au, ad = R.two_sided_sr(x, tau2, mu, md, [2 * A])[2 * A]
    first = min(v for v in (au[0], ad[0]) if v >= 0)
    assert ev and ev[0].index == cal + first
    assert ev[0].direction == ("up" if au[0] == first else "down")


def test_streaming_monitors_handle_missing_and_zero_returns():
    rng = np.random.default_rng(16)
    r = rng.standard_normal(6000)
    r[100:900] = 0.0  # frozen stretch longer than half a window
    r[2000:2100] = np.nan
    for mon in (R.RunningTauMonitor(threshold=1e6, cal_len=500, method="grid"),
                R.RunningTauMonitor(threshold=1e6, cal_len=500, method="dkw"),
                R.ConformalMonitor(threshold=1e6, cal_len=500)):
        mon.run(r)
        assert mon.calibrated


# ----------------------------------------------------------------------------
# No blind period
# ----------------------------------------------------------------------------


def test_levels_cs_matches_eb_confseq():
    rng = np.random.default_rng(17)
    r = rng.standard_t(4, 500)
    levels = R.grid_levels(np.array([3.0]))[0]
    cs = R._LevelsCS(levels, 0.001, n_plan=300)
    ups, los = [], []
    for v in r:
        cs.update(v)
        u, l_ = cs.bounds()
        ups.append(u)
        los.append(l_)
    ups, los = np.array(ups).T, np.array(los).T
    x = np.minimum(r[None, :] ** 2, (levels**2)[:, None])
    ref = R.eb_confseq(x, levels**2, 0.001, 0.001, n_plan=300)
    np.testing.assert_allclose(ups, ref.upper, rtol=1e-10)
    np.testing.assert_allclose(los, ref.lower, rtol=1e-10, atol=1e-14)


def test_awake_streaming_matches_vectorised():
    rng = np.random.default_rng(18)
    pre, n = 120, 2500
    r = rng.standard_normal(pre + n)
    r[pre + 300 :] *= 1.8  # change soon after the start, while the band is still learning
    A = 400.0
    mon = R.AwakeRunningTauMonitor(threshold=A, delta=0.01)
    ev = mon.run(r)
    out, _ = R.awake_route_a(r[None, :pre], r[None, pre:], [2 * A], delta=0.01)
    au, ad = out[2 * A]
    first = min(v for v in (au[0], ad[0]) if v >= 0)
    assert ev and ev[0].index == pre + first and ev[0].direction == "up"
    assert mon.segments[0]["start"] == pre


def test_awake_route_a_false_alarms_respect_delta_plus_m_over_A():
    rng = np.random.default_rng(19)
    reps, pre, m, A, delta = 500, 120, 300, 3000.0, 0.05
    r = rng.standard_t(3, (reps, pre + m))
    out, _ = R.awake_route_a(r[:, :pre], r[:, pre:], [2 * A], delta=delta)
    au, ad = out[2 * A]
    rate = np.mean((au >= 0) | (ad >= 0))
    bound = delta + m / A
    assert rate <= bound + 3 * math.sqrt(bound * (1 - bound) / reps)


def test_conformal_monitor_restarts_without_warmup():
    rng = np.random.default_rng(20)
    r = rng.standard_normal(4000)
    r[1500:] *= 3.0
    mon = R.ConformalMonitor(threshold=50000.0, cal_len=500, restart_len=0)
    ev = mon.run(r)
    assert ev and ev[0].direction == "up" and 1500 <= ev[0].index < 1600
    assert len(mon.segments) >= 2 and mon.segments[1]["start"] == ev[0].index + 1
