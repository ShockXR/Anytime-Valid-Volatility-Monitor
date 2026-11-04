"""Checks of the E21 benchmarks (``avmonitor.benchmarks_robust``) against brute-force versions (no data needed)."""

import math

import numpy as np
import pytest
from scipy import optimize, stats

from avmonitor import benchmarks_robust as BR
from avmonitor import competitors as C
from avmonitor import detectors as D
from avmonitor.monitor import VolatilityMonitor


def _paths(seed=1, reps=6, n_cal=500, n=300):
    rng = np.random.default_rng(seed)
    r = rng.standard_normal((reps, n_cal + n)) * np.r_[np.ones(n_cal + n // 2), 1.4 * np.ones(n - n // 2)]
    r[2, n_cal + 50 : n_cal + 55] = np.nan
    r[3, n_cal + 80] = 30.0  # a jump
    return r[:, :n_cal], r[:, n_cal:]


def test_truncated_inputs_match_detectors_calibrate():
    cal, mon = _paths()
    ti = BR.truncated_inputs(mon, cal)
    for i in range(cal.shape[0]):
        b = D.calibrate(cal[i], c=3.0)
        assert ti.tau2[i] == pytest.approx(b.tau**2)
        assert ti.m0[i] == pytest.approx(b.m0)
        assert ti.B[i] == pytest.approx(b.tau**2 / b.m0)
    assert np.nanmax(ti.y - ti.B[:, None]) <= 1e-12


def test_sr_log_paths_equal_the_betting_recursion_with_missing_values():
    cal, mon = _paths()
    ti = BR.truncated_inputs(mon, cal)
    lam_d = BR.down_bets(ti.B)
    lu = BR.sr_log_paths(ti.y, D.DEFAULT_UP_GRID, "up", chunk=4)
    ld = BR.sr_log_paths(ti.y, lam_d, "down", chunk=4)
    gu = np.asarray(D.DEFAULT_UP_GRID)
    for i in range(ti.y.shape[0]):
        Su, Sd = np.zeros(gu.size), np.zeros(lam_d.shape[1])
        for t in range(ti.y.shape[1]):
            y = ti.y[i, t]
            if not np.isfinite(y):
                assert np.isnan(lu[i, t]) and np.isnan(ld[i, t])
                continue
            Su = (1 + gu * (y - 1)) * (Su + 1)
            Sd = (1 + lam_d[i] * (1 - y)) * (Sd + 1)
            assert lu[i, t] == pytest.approx(np.log(Su.mean()), abs=1e-9)
            assert ld[i, t] == pytest.approx(np.log(Sd.mean()), abs=1e-9)
    # the same alarms as detectors.betting_edetector at a common level (one path, its own bets)
    i = 0
    up = D.betting_edetector(ti.y[i], 1.0, "up", threshold=50.0)
    a, _ = BR.first_alarm(lu[i : i + 1], np.full((1, lu.shape[1]), -np.inf), np.log(50.0), np.inf)
    assert a[0] == up.alarm[0]


def test_page_cusum_paths_equal_the_recursion():
    cal, mon = _paths(2)
    ti = BR.truncated_inputs(mon, cal)
    Gu, Gd = BR.page_cusum_paths(ti.y, 0.3, 0.2)
    for i in range(ti.y.shape[0]):
        g1 = g2 = 0.0
        for t in range(ti.y.shape[1]):
            y = ti.y[i, t]
            if not np.isfinite(y):
                assert np.isnan(Gu[i, t])
                continue
            g1, g2 = max(0.0, g1 + y - 1.3), max(0.0, g2 + 1 - y - 0.2)
            assert Gu[i, t] == pytest.approx(g1) and Gd[i, t] == pytest.approx(g2)


def test_gaussian_reference_value_is_the_likelihood_ratio_cusum():
    # increments of the Gaussian log likelihood ratio for a variance factor theta are proportional to z^2 - 1 - k
    for ratio in (1.25, 1.5, 0.8, 0.67):
        th = ratio**2
        k = BR.gaussian_reference_value(ratio)
        z2 = np.linspace(0, 5, 11)
        llr = -0.5 * np.log(th) + 0.5 * (1 - 1 / th) * z2
        lin = (z2 - 1 - k) if th > 1 else (1 - z2 - k)
        assert np.allclose(llr / lin, llr[0] / lin[0])
        assert k > 0
    assert BR.gaussian_reference_value(1.25) == pytest.approx(1.5625 * math.log(1.5625) / 0.5625 - 1)


def test_first_alarm_direction_and_ties():
    up = np.array([[0, 1, 5, 0], [0, 0, 0, 0], [0, 3, 0, 0]], float)
    dn = np.array([[0, 2, 0, 0], [0, 0, 0, 1], [0, 3, 0, 0]], float)
    a, d = BR.first_alarm(up, dn, 3.0, 2.0)
    assert list(a) == [1, -1, 1] and list(d) == [-1.0, 0.0, 1.0]


def test_bartlett_lrv_brute_force():
    rng = np.random.default_rng(3)
    z = rng.standard_normal((2, 400)).cumsum(1) * 0.1 + rng.standard_normal((2, 400))
    q = 5
    for i in range(2):
        d = z[i] - z[i].mean()
        g = [d[l:] @ d[: d.size - l] / d.size for l in range(q + 1)]
        want = g[0] + 2 * sum((1 - l / (q + 1)) * g[l] for l in range(1, q + 1))
        assert BR.bartlett_lrv(z[i : i + 1], q)[0] == pytest.approx(max(want, g[0] / 10))
    # default bandwidth floor(n^(1/3)) = 7 for n = 400
    assert BR.bartlett_lrv(z[:1])[0] == pytest.approx(BR.bartlett_lrv(z[:1], 7)[0])


def test_hkz_paths_brute_force_and_csw_equivalence():
    rng = np.random.default_rng(4)
    zc, zm = rng.standard_normal((1, 300)) ** 2, rng.standard_normal((1, 120)) ** 2
    zm[0, 10] = np.nan
    s = math.sqrt(BR.bartlett_lrv(zc)[0]) * math.sqrt(300)
    Dh = BR.hkz_paths(zc, zm, gamma=0.25)
    Dc = BR.hkz_paths(zc, zm, boundary="csw")
    fin = zm[0][np.isfinite(zm[0])]
    j = 0
    for t in range(zm.shape[1]):
        if not np.isfinite(zm[0, t]):
            assert np.isnan(Dh[0, t])
            continue
        k = j + 1
        Q = fin[:k].sum() - k * zc.mean()
        assert Dh[0, t] == pytest.approx(Q / (s * (1 + k / 300) * (k / (300 + k)) ** 0.25))
        x = (300 + k) / 300
        for c in (0.5, 1.5, 2.8):
            crosses = abs(Q) / s >= math.sqrt(x * (x - 1) * (c * c + math.log(x / (x - 1))))
            assert crosses == (abs(Dc[0, t]) >= c)
        j += 1


def test_critical_values():
    c = BR.csw_crit(0.05)
    assert 2 * (1 - stats.norm.cdf(c) + c * stats.norm.pdf(c)) == pytest.approx(0.05)
    assert c == pytest.approx(2.795, abs=0.001)
    exact = optimize.brentq(lambda x: BR.sup_abs_bm_cdf(x) - 0.95, 1.0, 4.0)
    assert exact == pytest.approx(2.2414, abs=1e-3)
    sim = BR.wiener_sup_crit(0.0, 0.05, n_paths=4000, rng=np.random.default_rng(5))
    assert sim == pytest.approx(exact, abs=0.05)
    # the boundary function makes the limit larger as gamma grows
    s1 = BR.wiener_sup_samples(0.25, 2000, rng=np.random.default_rng(6))
    s2 = BR.wiener_sup_samples(0.45, 2000, rng=np.random.default_rng(6))
    assert np.all(s2 >= s1 - 1e-12)


def test_grid_bocpd_gaussian_matches_the_conjugate_filter_before_the_cap():
    rng = np.random.default_rng(7)
    x = rng.standard_normal((3, 200)) * np.r_[np.ones(100), 2 * np.ones(100)]
    _, S1, L1 = C.bocpd_paths(x, 0.01, alpha0=2.0, prior_var=1.0, r_max=512, short=40, return_log_odds=True)
    L2, S2 = BR.grid_bocpd_paths(x, 0.01, BR.log_variance_grid(1e-3, 1e3, 0.01), 40, nu=None, prior_var=1.0)
    ok = np.isfinite(L1) & np.isfinite(L2)
    assert np.max(np.abs(L1[ok] - L2[ok])) < 1e-4
    assert np.mean(S1 == S2) > 0.99


def test_grid_bocpd_lumping_is_exact():
    """The short-run log-odds with runs of at least w lumped equals the unlumped recursion (w larger than the series)."""
    rng = np.random.default_rng(8)
    x = rng.standard_t(3, (2, 60)) * np.r_[np.ones(30), 1.8 * np.ones(30)]
    x[1, 12] = np.nan
    grid = BR.log_variance_grid(0.05, 20, 0.2)
    w = 6
    fa = BR.GridBOCPD(2, 0.05, grid, w, nu=3.0, prior_var=1.0)
    fb = BR.GridBOCPD(2, 0.05, grid, 80, nu=3.0, prior_var=1.0)
    for t in range(x.shape[1]):
        fa.update(x[:, t])
        fb.update(x[:, t])
        with np.errstate(divide="ignore"):  # before step w the long-run mass is exactly 0
            lo_b = np.log(fb.J[:, :w, :].sum((1, 2))) - np.log(fb.J[:, w:, :].sum((1, 2)))
            lo_a = fa.short_log_odds()
        assert np.array_equal(np.isinf(lo_a), np.isinf(lo_b))
        fin = np.isfinite(lo_a)
        assert np.allclose(lo_a[fin], lo_b[fin], rtol=1e-9, atol=1e-9)
        assert np.allclose(fa.J.sum((1, 2)), 1.0)


def test_grid_bocpd_student_t_ignores_an_isolated_jump():
    rng = np.random.default_rng(9)
    cal = rng.standard_normal((1, 1000))
    x = rng.standard_normal((1, 300))
    x[0, 150] = 25.0
    grid = BR.log_variance_grid(0.04, 25, 0.1)
    Lt, _ = BR.grid_bocpd_paths(x, 1 / 7800, grid, 50, nu=3.0, cal=cal)
    Lg, _ = BR.grid_bocpd_paths(x, 1 / 7800, grid, 50, nu=None, cal=cal)
    assert Lg[0, 150] > 5 and Lt[0, 150] < 0  # Gaussian: a change is near certain; Student t3: not


def test_bootstrap_max_stats_constant_window_and_against_sr_log_paths():
    xcal = np.full((2, 100), 0.7)
    out = BR.bootstrap_max_stats(xcal, np.array([4.0, 4.0]), {"sr": ("sr", {}), "c": ("cusum", dict(k_up=0.2, k_dn=0.2))},
                                 50, 7, np.random.default_rng(0))
    assert np.allclose(out["sr"][0], np.log(50)) and np.allclose(out["sr"][1], np.log(50))  # f = 1: M_t = t
    assert np.allclose(out["c"][0], 0.0) and np.allclose(out["c"][1], 0.0)
    rng = np.random.default_rng(10)
    v = np.minimum(rng.standard_normal(300) ** 2, 9.0)
    out = BR.bootstrap_max_stats(v[None, :], np.array([9.0]), {"sr": ("sr", {})}, 40, 5, np.random.default_rng(11))
    rr = np.random.default_rng(11)  # replay the same draws
    m0s = v[rr.integers(0, 300, (5, 300))].mean(axis=1)
    ys = v[rr.integers(0, 300, (5, 40))] / m0s[:, None]
    lu = BR.sr_log_paths(ys, D.DEFAULT_UP_GRID, "up")
    ld = BR.sr_log_paths(ys, BR.down_bets(9.0 / m0s), "down")
    assert np.allclose(out["sr"][0][0], lu.max(1)) and np.allclose(out["sr"][1][0], ld.max(1))


def _series(seed=12, n=40000):
    rng = np.random.default_rng(seed)
    u = rng.standard_normal(n) * np.repeat(rng.choice([0.6, 1.0, 1.6], n // 2000), 2000)
    u[rng.random(n) < 0.01] = np.nan
    u[5000:5400] = 0.0  # a stretch of zeros: degenerate windows are extended
    return u


@pytest.mark.parametrize("band", [0.0, 0.25])
def test_restart_events_reproduce_volatility_monitor(band):
    u = _series()
    A = 300.0
    ev = [(e.index, e.direction) for e in VolatilityMonitor(threshold=A, cal_len=300, null_band=band).run(u)]
    got = BR.restart_events(u, 300, "sr", math.log(2 * A), math.log(2 * A), dict(band=band), block=777)
    assert got == ev and len(ev) > 5


@pytest.mark.parametrize("kind,params,h", [("cusum", dict(k_up=0.24, k_dn=0.2), 15.0),
                                           ("hkz", dict(on="x", gamma=0.25), 2.0),
                                           ("hkz", dict(on="r2", boundary="csw"), 2.5)])
def test_restart_events_brute_force(kind, params, h):
    u = _series(13, 20000)
    got = BR.restart_events(u, 300, kind, h, h, params, block=500)
    # brute force: recompute each segment from its calibration window with the vectorised path functions
    want, s = [], 0
    while True:
        cw = BR.calibrate_window(u, s, 300)
        if cw is None:
            break
        e, _, tau2, m0, rfin = cw
        mon = u[e:]
        fin = np.flatnonzero(np.isfinite(mon))
        r = mon[fin][None, :]
        if kind == "cusum":
            up, dn = BR.page_cusum_paths(np.minimum(r**2, tau2) / m0, params["k_up"], params["k_dn"])
        else:
            zc = rfin**2 if params["on"] == "r2" else np.minimum(rfin**2, tau2)
            zm = r**2 if params["on"] == "r2" else np.minimum(r**2, tau2)
            Dk = BR.hkz_paths(zc[None, :], zm, gamma=params.get("gamma", 0.25), boundary=params.get("boundary", "hhks"))
            up, dn = Dk, -Dk
        a, d = BR.first_alarm(up, dn, h, h)
        if a[0] < 0:
            break
        t = e + int(fin[a[0]])
        want.append((t, "up" if d[0] > 0 else "down"))
        s = t + 1
    assert got == want and len(want) > 3


def test_match_restart_counts_hits_the_targets():
    u = _series(14, 30000)
    ev = BR.restart_events(u, 300, "cusum", 15.0, 15.0, dict(k_up=0.24, k_dn=0.2))
    nu_ = sum(1 for e in ev if e[1] == "up")
    hu, hd, got = BR.match_restart_counts(u, 300, "cusum", dict(k_up=0.24, k_dn=0.2), nu_, len(ev) - nu_, (0.0, 200.0))
    # with recalibration after each alarm the counts are not exactly monotone in the thresholds: within one alarm
    got_up = sum(1 for e in got if e[1] == "up")
    assert abs(got_up - nu_) <= 1 and abs(len(got) - got_up - (len(ev) - nu_)) <= 1


def test_match_restart_counts_local_and_grid_search_never_worse():
    u = _series(15, 30000)
    prm = dict(on="x", gamma=0.25)
    ev = BR.restart_events(u, 300, "hkz", 1.2, 0.9, prm)
    tu, td = sum(1 for e in ev if e[1] == "up"), sum(1 for e in ev if e[1] == "down")

    def err(got):
        k = sum(1 for e in got if e[1] == "up")
        return abs(k - tu) + abs(len(got) - k - td)

    _, _, plain = BR.match_restart_counts(u, 300, "hkz", prm, tu, td, (0.0, 50.0), iters=20, rounds=2)
    _, _, refined = BR.match_restart_counts(u, 300, "hkz", prm, tu, td, (0.0, 50.0), iters=20, rounds=2, refine=4,
                                            grid_points=11, grid_if_error_above=0)
    assert err(refined) <= err(plain) and err(refined) <= 1
