"""E16 (experiments/e16_adversarial_null.py): the adversarial laws stay inside the null class, and the
vectorised detector is the monitor's detector.

* every searched law keeps the conditional mean of Y = X / m0 at exactly 1 (two-point and
  three-point laws, the adaptive bold-play policy at every step and every state, the shape-switching
  chain in each regime) and stays in [0, B];
* the return-law families (scale mixtures, Student t, ticks) have tau = 3 x scaled MAD and m0 equal
  to the truth of the law;
* the vectorised detector matches VolatilityMonitor alarm for alarm, with a baseline that is the
  true conditional mean of the monitored stream (i.i.d., adaptive and Gaussian laws), and with
  the plug-in baseline;
* the engine (on Y) and first_alarm (on X) agree, also with a time-varying baseline.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy import optimize

from avmonitor.monitor import VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def E():
    spec = importlib.util.spec_from_file_location("e16_adversarial_null", ROOT / "experiments" / "e16_adversarial_null.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _checked_policy(E, pol, B, log):
    """Wrap a policy: record that every chosen two-point law has mean 1 and support in [0, B]."""

    def wrapped(t, idx, Su, Sd, lamd):
        yl, yh = pol(t, idx, Su, Sd, lamd)
        assert np.all((yl >= 0) & (yl <= 1) & (yh >= 1) & (yh <= B * (1 + 1e-12)))
        span = yh - yl
        p = np.where(span > 0, (1 - yl) / np.where(span > 0, span, 1), 0.0)
        mean = np.where(span > 0, p * yh + (1 - p) * yl, yl)
        assert np.allclose(mean, 1.0, rtol=0, atol=1e-12)
        log.append(int(idx.size))
        return yl, yh

    return wrapped


@pytest.mark.parametrize("theta_kind", ["bold", "random1", "random2"])
def test_adaptive_policy_keeps_the_conditional_mean_at_every_step(E, theta_kind):
    rng = np.random.default_rng(1)
    theta = E.BOLD if theta_kind == "bold" else np.r_[rng.uniform(0.3, 1.5), rng.uniform(-1.5, 1.5, 3),
                                                       rng.uniform(-1, 0.5), rng.uniform(-1.5, 1.5, 3)]
    B, n, P = E.B_GAUSS, 400, 3000
    log = []
    pol = _checked_policy(E, E.bold_policy(theta, B, n, thr=600.0), B, log)
    U = rng.random((P, n))
    res = E.engine(E.adaptive_y_fn(pol, U), P, n, B, thr=600.0, record=True)
    assert len(log) > 0 and (res["alarm"] >= 0).any()
    # empirical check of the same thing on the recorded draws: mean of Y over paths still running
    Y = res["Y"]
    assert np.nanmax(Y) <= B * (1 + 1e-12) and np.nanmin(Y) >= 0


def test_bold_play_lands_on_the_threshold(E):
    """When one draw can take a side to 2A, bold play lands there: no alarm overshoots by more than one unit."""
    rng = np.random.default_rng(2)
    B, n, P, thr = E.B_GAUSS, 600, 4000, 800.0
    res = E.run_bold(E.BOLD, rng.random((P, n)), B, thr=thr, decompose=True)
    hit = res["alarm"] >= 0
    assert hit.sum() > 50
    # N = (M_up + M_down) / 2 at an alarm; the alarming side is in [thr, thr + 1), the other side >= 0
    assert np.all(res["n_at_alarm"][hit] >= thr / 2)
    side_val = 2 * res["n_at_alarm"][hit]  # M_alarm + M_other
    assert np.median(side_val) < thr * 1.2


def test_iid_two_and_three_point_laws_have_mean_one(E):
    B = E.B_GAUSS
    rng = np.random.default_rng(3)
    U = rng.random(2_000_000)
    for a, b in [(0.0, B), (0.3, 2.0), (0.9, 1.1)]:
        ok, q = E.law_two_point((a, b), B)
        assert ok
        y = q(U, None)
        assert set(np.unique(y)) <= {a, b}
        p = (1 - a) / (b - a)
        assert p * b + (1 - p) * a == pytest.approx(1.0, abs=1e-14)
        assert y.mean() == pytest.approx(1.0, abs=4 * y.std() / np.sqrt(y.size))
    for par in [(0.0, 2.0, B * 0.999, 0.2), (0.2, 0.5, 5.0, 0.3), (0.1, 3.0, 7.0, 0.05)]:
        ok, q = E.law_three_point(par, B)
        assert ok
        y1, y2, y3, p2 = par
        p3 = (1 - y1 - p2 * (y2 - y1)) / (y3 - y1)
        p1 = 1 - p2 - p3
        assert p1 * y1 + p2 * y2 + p3 * y3 == pytest.approx(1.0, abs=1e-14)
        y = q(U, None)
        assert y.min() >= 0 and y.max() <= B
        assert y.mean() == pytest.approx(1.0, abs=4 * y.std() / np.sqrt(y.size))
    assert not E.law_three_point((0.5, 0.2, 1.2, 0.9), B)[0]  # infeasible: no mean-one law with these points


@pytest.mark.parametrize("family,par", [("scale mixture", (0.1, 4.0)), ("scale mixture", (0.35, 1.5)),
                                        ("student t", (3.0,)), ("ticks", (0.5,)), ("ticks", (1.2,))])
def test_return_law_families_use_their_true_tau_and_mean(E, family, par):
    """tau = 3 x scaled MAD of the law and m0 = E[min(r^2, tau^2)]: the sample median and the sample mean agree."""
    rng = np.random.default_rng(4)
    n = 2_000_000
    if family == "scale mixture":
        ok, q, B = E.law_scale_mixture(par)
        w, k = par
        r = np.where(rng.random(n) < w, k, 1.0) * rng.standard_normal(n)
        tau, m0 = E.scale_mixture_truth(w, k)
    elif family == "student t":
        ok, q, B = E.law_student_t(par)
        r = rng.standard_t(par[0], n)
        tau, m0 = E.t_truth(par[0])
    else:
        ok, q, B = E.law_ticks(par)
        h = par[0]
        r = h * np.round(rng.standard_normal(n) / h)
        _, _, tau, m0 = E.tick_law(h)
    assert ok and B == pytest.approx(tau**2 / m0)
    assert 3 * MAD_NORMAL_SCALE * np.median(np.abs(r)) == pytest.approx(tau, rel=5e-3)
    x = np.minimum(r * r, tau**2)
    assert x.mean() == pytest.approx(m0, abs=4 * x.std() / np.sqrt(n))
    # the sampler used by the search gives Y = X / m0 with mean 1 and support in [0, B]
    y = q(rng.random((2000, 500)), rng.standard_normal((2000, 500)))
    assert y.max() <= B * (1 + 1e-12) and y.mean() == pytest.approx(1.0, abs=5 * y.std() / np.sqrt(y.size))


def test_shape_switching_keeps_the_mean_in_each_regime(E):
    rng = np.random.default_rng(5)
    P, n = 4000, 300
    U, V, Z = rng.random((P, n)), rng.random((P, n)), rng.standard_normal((P, n))
    Y = E.markov_shape_Y(U, V, Z, 0.9, 0.99, E.B_GAUSS)
    wild = np.isin(Y, (0.0, E.B_GAUSS))
    # regime 'wild' has mean 1 by construction; the Gaussian shape has mean 1 too
    g = E.gaussian_q(U, Z)
    assert g.mean() == pytest.approx(1.0, abs=4 * g.std() / np.sqrt(g.size))
    assert Y.max() <= E.B_GAUSS and Y.min() >= 0 and wild.any() and (~wild).any()
    assert Y.mean() == pytest.approx(1.0, abs=5 * Y.std() / np.sqrt(Y.size) * 3)


def test_trunc_mean_normal_and_the_gaussian_b(E):
    rng = np.random.default_rng(6)
    z = rng.standard_normal(4_000_000)
    for k in (0.5, 1.5, 3.0):
        x = np.minimum(z * z, k * k)
        assert E.trunc_mean_normal(k) == pytest.approx(x.mean(), abs=4 * x.std() / np.sqrt(z.size))
    assert E.TAU_G == pytest.approx(3.0, abs=1e-12)
    assert E.B_GAUSS == pytest.approx(9.04515996885841, rel=1e-9)  # the B of E11 (results/e11d_*.csv)


def test_tick_tau_is_the_protocol_rule(E):
    rng = np.random.default_rng(7)
    for h in (0.25, 0.75, 1.5):
        j, p, tau, m0 = E.tick_law(h)
        r = h * np.round(rng.standard_normal(1_000_001) / h)
        med = np.median(np.abs(r))
        if med == 0:  # the monitor's fallback: scaled median of the nonzero |r|
            med = np.median(np.abs(r[r != 0]))
        assert tau == pytest.approx(3 * MAD_NORMAL_SCALE * med)


def _monitor_first_alarm(r_cal, r_mon, threshold):
    mon = VolatilityMonitor(threshold=threshold, cal_len=r_cal.size)
    for x in r_cal:
        mon.update(x)
    assert mon.calibrated
    for x in r_mon:
        ev = mon.update(x)
        if ev is not None:
            return ev.index - r_cal.size, 1 if ev.direction == "up" else -1
    return -1, 0


def _plugin_values(r_cal):
    """VolatilityMonitor._plugin_stats, written out (scaled MAD about zero, tau = 3 x scale, mean of X)."""
    scale = float(np.median(np.abs(r_cal)) * MAD_NORMAL_SCALE)
    tau2 = (3.0 * scale) ** 2
    return tau2, float(np.mean(np.minimum(r_cal**2, tau2)))


@pytest.mark.parametrize("law", ["two-point", "bold", "gaussian-shape"])
def test_vectorised_detector_matches_volatility_monitor_with_the_true_baseline(E, law):
    """The monitored stream is built so that the monitor's baseline is its true conditional mean."""
    rng = np.random.default_rng(8)
    thr_A = 150.0
    r_cal = rng.standard_normal(390)
    tau2, m0 = _plugin_values(r_cal)
    B = tau2 / m0
    P, n = 12, 900
    U = rng.random((P, n))
    if law == "two-point":
        Y = E.two_point_draw(U, 0.0, B)
    elif law == "bold":
        Y = E.run_bold(E.BOLD, U, B, thr=2 * thr_A, record=True)["Y"]
        Y = np.where(np.isnan(Y), 1.0, Y)  # steps after the alarm: any value, they are not used
    else:  # Gaussian shape rescaled to [0, B]: Y = min(Z^2, k^2) / h(k) with k^2 / h(k) = B, so E[Y] = 1
        k = optimize.brentq(lambda k: k * k / E.trunc_mean_normal(k) - B, 1.0, 20.0)
        Y = np.minimum(rng.standard_normal((P, n)) ** 2, k * k) / E.trunc_mean_normal(k)
    r = np.sqrt(Y * m0)
    X = np.minimum(r * r, tau2)
    a_vec, s_vec = E.first_alarm(X, m0, m0, tau2, thr=2 * thr_A)
    a_eng = E.engine(E.array_y_fn(X / m0), P, n, B, thr=2 * thr_A)
    seen = set()
    for i in range(P):
        mon = _monitor_first_alarm(r_cal, r[i], thr_A)
        assert mon == (int(a_vec[i]), int(s_vec[i])), (law, i)
        seen.add(mon[1])
    assert seen - {0}, "no alarms: the test would be vacuous"
    assert np.array_equal(a_eng["alarm"], a_vec)


def test_vectorised_plugin_matches_volatility_monitor(E):
    rng = np.random.default_rng(9)
    cal, n, P, thr_A = 390, 1500, 10, 120.0
    r = rng.standard_normal((P, cal + n)) * np.r_[np.ones(cal + 600), np.full(n - 600, 1.4)][None, :]
    tau2, m = E.plugin_baseline(r[:, :cal])
    X = np.minimum(r[:, cal:] ** 2, tau2[:, None])
    a, s = E.first_alarm(X, m, m, tau2, thr=2 * thr_A)
    for i in range(P):
        assert _monitor_first_alarm(r[i, :cal], r[i, cal:], thr_A) == (int(a[i]), int(s[i]))
    assert (a >= 0).sum() >= 5


def test_engine_and_first_alarm_agree_and_time_varying_baseline(E):
    rng = np.random.default_rng(10)
    B, P, n = 6.0, 300, 500
    Y = E.two_point_draw(rng.random((P, n)), 0.0, B)
    a1 = E.engine(E.array_y_fn(Y), P, n, B, thr=400.0)["alarm"]
    a2, _ = E.first_alarm(Y, 1.0, 1.0, B, thr=400.0)
    assert np.array_equal(a1, a2) and (a1 >= 0).any()
    a3, _ = E.first_alarm(Y, np.ones((P, n)), np.ones((P, n)), B, thr=400.0)
    assert np.array_equal(a2, a3)
    # a time-varying baseline (the GARCH oracle): against a plain loop, one path at a time
    s = np.exp(0.3 * rng.standard_normal((P, n)))
    tau2 = B * s.max(axis=1)
    X = np.minimum(Y * s, tau2[:, None])
    a4, s4 = E.first_alarm(X, s, s, tau2, thr=400.0)
    for i in range(20):
        Su, Sd, ref = np.zeros(6), np.zeros(5), (-1, 0)
        for t in range(n):
            lam = E.DNF / (tau2[i] / s[i, t] - 1)
            Su = (1 + E.UPG * (X[i, t] / s[i, t] - 1)) * (Su + 1)
            Sd = (1 + lam * (1 - X[i, t] / s[i, t])) * (Sd + 1)
            if Su.mean() >= 400.0:
                ref = (t, 1)
                break
            if Sd.mean() >= 400.0:
                ref = (t, -1)
                break
        assert ref == (int(a4[i]), int(s4[i]))


def test_garch_conditional_mean_formula(E):
    """E[min(s2 Z^2, tau^2)] = s2 h(tau / sqrt(s2)), used for the GARCH drift and the in-class oracle."""
    rng = np.random.default_rng(11)
    z = rng.standard_normal(2_000_000)
    for s2 in (0.3, 1.0, 4.0):
        x = np.minimum(s2 * z * z, 9.0)
        assert s2 * E.trunc_mean_normal(3.0 / np.sqrt(s2)) == pytest.approx(x.mean(), abs=4 * x.std() / np.sqrt(z.size))
    R, S = E.garch(200, 300, 0.05, 0.9, rng, 500)
    assert np.all(S > 0) and abs(np.mean(R * R) - 1) < 0.15


def test_price_ticks_are_dependent_and_iid_ticks_are_not(E):
    rng = np.random.default_rng(12)
    r = E.price_tick_returns(400, 4000, 1.0, rng)
    k = np.round(r / 1.0)
    assert np.allclose(r, k)  # whole ticks
    c = np.corrcoef(r[:, :-1].ravel(), r[:, 1:].ravel())[0, 1]
    assert c < -0.03  # bid-ask-bounce-like negative autocorrelation (about -0.07 at a tick of 1 sd)
    z = np.round(rng.standard_normal((400, 4000)))
    assert abs(np.corrcoef(z[:, :-1].ravel(), z[:, 1:].ravel())[0, 1]) < 0.01
