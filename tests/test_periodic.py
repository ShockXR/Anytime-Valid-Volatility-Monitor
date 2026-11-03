"""TH6 (src/avmonitor/periodic.py, docs/theory_per_period.md): per-period false alarm guarantees.

* the closed-form bounds (Theorem 1, Theorem 2, the threshold that makes the window bound 5%);
* the engine against the package's own two-sided detector (Shiryaev-Roberts and CUSUM recursions);
* Proposition 1(b): the constant null law Y = 1 makes the plain detector alarm exactly at step 2A;
* Lemma 1 (time reversal) as an identity on paths and as a tail bound by Monte Carlo;
* Theorem 1 (restart at every period boundary) on the calm-then-noise stream where the plain detector fails;
* Theorem 2 (single bet, i.i.d. null, a late window) by Monte Carlo;
* Proposition 3 (park and strike): no alarm before the strike, every draw has mean 1, the late-window rate is
  above the i.i.d. bound and below the conditional bound of Lemma 2;
* PeriodicVolatilityMonitor: identical to VolatilityMonitor by default, and restarts at period boundaries;
  with baseline="known" it is exactly the engine's detector, and Theorem 1 holds for it on the calm stream.
"""

import math

import numpy as np
import pytest

from avmonitor import periodic as Pm
from avmonitor.detectors import betting_edetector, two_sided_alarm
from avmonitor.monitor import VolatilityMonitor
from avmonitor.robustness import gaussian_B, trunc_mean_normal

B = gaussian_B()
M0 = float(trunc_mean_normal(3.0))


def test_closed_form_bounds():
    assert Pm.entropy(np.full(6, 1 / 6)) == pytest.approx(math.log(6))
    assert Pm.window_bound_two_sided(390, 7800.0) == pytest.approx((390 + math.log(15600) + (math.log(6) + math.log(5)) / 2) / 7800)
    assert Pm.window_bound_two_sided(390, 7800.0) == pytest.approx(0.051456, abs=1e-6)
    assert Pm.window_bound_two_sided(1440, 28800.0) == pytest.approx(0.050440, abs=1e-6)
    assert Pm.window_bound_single(390, 7800.0) == pytest.approx(0.05)
    assert Pm.restarted_sr_bound(390, 7800.0) == pytest.approx(0.05)
    assert Pm.restarted_sr_bound(390, 15600.0, carry=390) == pytest.approx(0.05)
    a = Pm.threshold_for_window(390, 0.05)
    assert Pm.window_bound_two_sided(390, a) == pytest.approx(0.05, abs=1e-9)
    assert 8000 < a < 8050
    # the union-bound form is never better for m >= 1
    w = np.full(4, 0.25)
    assert Pm.window_bound_mixture(10, 1000.0, w) == pytest.approx((10 + math.log(1000) + math.log(4)) / 1000)


@pytest.mark.parametrize("kind", ["SR", "CUSUM"])
def test_engine_matches_package_detector(kind):
    rng = np.random.default_rng(1)
    P, n, A = 200, 1500, 150.0
    x = np.minimum(rng.standard_normal((P, n)) ** 2, 9.0)
    up = betting_edetector(x, M0, "up", threshold=2 * A, kind=kind)
    dn = betting_edetector(x, M0, "down", upper=9.0, threshold=2 * A, kind=kind)
    ref = two_sided_alarm(up, dn)
    y = x / M0
    res = Pm.run_engine(lambda t, e: y[:, t], P, n, B, [Pm.Spec("d", kind, A)], day=100)
    assert np.array_equal(ref, res["d"]["first"])
    assert (ref >= 0).mean() > 0.05  # the comparison covers alarms


def test_constant_null_law_alarms_at_2A():
    # Proposition 1(b): Y_t = 1 is in the null on both sides and makes every Shiryaev-Roberts statistic equal to t,
    # so the two-sided rule (each side at 2A) alarms exactly at step 2A; the restarted and CUSUM versions never do
    P, A, m = 3, 50.0, 10
    specs = [Pm.Spec("sr", "SR", A), Pm.Spec("rsr", "SR", A, period=m), Pm.Spec("cu", "CUSUM", A)]
    res = Pm.run_engine(lambda t, e: np.ones(P), P, 300, B, specs, day=m)
    assert np.all(res["sr"]["first"] == int(2 * A) - 1)  # 0-based index of step 2A
    assert np.all(res["rsr"]["first"] < 0) and np.all(res["cu"]["first"] < 0)


def test_time_reversal_identity():
    rng = np.random.default_rng(2)
    for _ in range(50):
        f = rng.uniform(0.2, 2.5, rng.integers(1, 40))
        c = 0.0
        for v in f:
            c = v * max(c, 1.0)
        assert c == pytest.approx(Pm.reversed_running_max(f))


def test_lemma1_tail_single_bet():
    # i.i.d. two-point law {0, B} with mean 1, bet 0.5: P(C_s >= x) <= 1 / x at a late time s
    rng = np.random.default_rng(3)
    P, s, lam = 40000, 400, 0.5
    C = np.zeros(P)
    for _ in range(s):
        y = np.where(rng.random(P) < 1 / B, B, 0.0)
        C = (1 + lam * (y - 1)) * np.maximum(C, 1.0)
    for x in (2.0, 5.0, 20.0):
        p = np.mean(C >= x)
        assert p <= 1 / x + 3 * math.sqrt(p * (1 - p) / P) + 1e-3


def test_restarted_sr_per_period_bound_on_calm_stream():
    # Y = 1 for 9 periods (inside the null class), then Gaussian: the plain detector carries 9 m into period 10
    rng = np.random.default_rng(4)
    P, m, A = 6000, 10, 100.0
    specs = [Pm.Spec("sr", "SR", A), Pm.Spec("rsr", "SR", A, period=m), Pm.Spec("cu", "CUSUM", A)]

    def y_fn(t, eng):
        if t < 9 * m:
            return np.ones(P)
        return np.minimum(rng.standard_normal(P) ** 2, 9.0) / M0

    res = Pm.run_engine(y_fn, P, 10 * m, B, specs, day=m)
    p_sr = res["sr"]["counts"][:, 9].astype(bool).mean()
    p_rsr = res["rsr"]["counts"][:, 9].astype(bool).mean()
    se = math.sqrt(0.1 * 0.9 / P)
    assert (res["sr"]["first"] < 0).sum() + (res["sr"]["first"] >= 9 * m).sum() == P  # no alarm while calm
    assert p_rsr <= m / A + 3 * se
    assert p_sr > m / A + 5 * se  # Theorem B gives only 10 m / A here
    assert p_sr <= min(1.0, 10 * m / A) + 3 * se


def test_cusum_window_bound_single_bet_iid():
    # Theorem 2(a): one increase bet, one side at level A_single, a late window of m steps
    rng = np.random.default_rng(5)
    P, A_single, m, s = 20000, 40.0, 4, 400
    spec = Pm.Spec("cu", "CUSUM", A_single / 2, on_alarm="continue")  # the engine alarms at 2 x Spec.A
    res = Pm.run_engine(lambda t, e: np.where(rng.random(P) < 1 / B, B, 0.0), P, s + m, B, [spec], day=m,
                        up=(0.5,), down_fractions=())
    p = (res["cu"]["counts"][:, -1] > 0).mean()
    assert p <= m / A_single + 3 * math.sqrt(p * (1 - p) / P)


def test_park_strike_counterexample_small():
    rng = np.random.default_rng(6)
    P, m, A = 4000, 10, 60.0
    spec = Pm.Spec("cu", "CUSUM", A)
    strike = 60 * m
    draws = []

    def y_fn(t, eng):
        y = Pm.park_strike_policy("cu", strike, rng)(t, eng)
        draws.append(y)
        return y

    res = Pm.run_engine(y_fn, P, strike + m, B, [spec], day=m)
    first = res["cu"]["first"]
    assert not np.any((first >= 0) & (first < strike))  # no alarm before the strike
    assert all(np.all((d >= 0) & (d <= B)) for d in draws)
    p = np.mean(first >= strike)
    iid = Pm.window_bound_two_sided(m, A)
    assert p > iid + 5 * math.sqrt(iid * (1 - iid) / P)  # Proposition 3: above the exchangeable bound


def test_two_point_draws_have_mean_one():
    rng = np.random.default_rng(7)
    yl = np.full(200000, 0.3)
    yh = np.full(200000, 4.0)
    y = Pm.two_point(rng.random(200000), yl, yh)
    assert y.mean() == pytest.approx(1.0, abs=0.01)
    assert set(np.unique(y)) == {0.3, 4.0}


def test_conditional_bound_holds_under_martingale_null():
    # Lemma 2 at the start of the strike window, against the park-and-strike stream
    rng = np.random.default_rng(8)
    P, m, A = 3000, 10, 60.0
    spec = Pm.Spec("cu", "CUSUM", A)
    eng = Pm.Engine([spec], P, B)
    pol = Pm.park_strike_policy("cu", 40 * m, rng)
    for t in range(40 * m):
        eng.step(pol(t, eng))
    bound = Pm.conditional_window_bound(eng.Su["cu"], eng.Sd["cu"], m, A)
    res = Pm.run_engine(Pm.bold_policy("cu", rng), P, m, B, [spec], day=m, t0=40 * m,
                        init={"cu": (eng.Su["cu"], eng.Sd["cu"])})
    p = np.mean(res["cu"]["first"] >= 0)
    assert p <= bound.mean() + 3 * math.sqrt(p * (1 - p) / P)


def test_periodic_monitor_default_is_identical():
    rng = np.random.default_rng(9)
    r = np.concatenate([rng.standard_normal(3000), 1.6 * rng.standard_normal(2000), rng.standard_normal(3000)])
    a = VolatilityMonitor(threshold=300, cal_len=200).run(r)
    b = Pm.PeriodicVolatilityMonitor(threshold=300, cal_len=200).run(r)
    assert [(e.index, e.direction, e.stat) for e in a] == [(e.index, e.direction, e.stat) for e in b]
    assert len(a) > 0


def test_periodic_monitor_restarts_and_cusum():
    rng = np.random.default_rng(10)
    r = rng.standard_normal(1500)
    mon = Pm.PeriodicVolatilityMonitor(threshold=1e9, cal_len=100, period=100)
    for t, v in enumerate(r):
        mon.update(v)
        if t >= 100 and t % 100 == 0:
            # this step was the first of a new period: the statistics restarted from 0 before it
            assert np.all(mon.S_up < 6) and np.all(mon.S_dn < 6)
    assert mon.period_resets == 13
    cu = Pm.PeriodicVolatilityMonitor(threshold=1e9, cal_len=100, statistic="cusum")
    sr = Pm.PeriodicVolatilityMonitor(threshold=1e9, cal_len=100)
    checked = 0
    for v in r:
        cu.update(v)
        sr.update(v)
        if sr.events or cu.events:
            break
        if sr.calibrated:  # same calibration window, so the same baseline: C_t <= M_t bet by bet
            assert np.all(cu.S_up <= sr.S_up + 1e-12) and np.all(cu.S_dn <= sr.S_dn + 1e-12)
            checked += 1
    assert checked > 100
    with pytest.raises(ValueError):
        Pm.PeriodicVolatilityMonitor(period=0)


@pytest.mark.parametrize("period", [50, None])
def test_periodic_monitor_known_baseline_matches_engine(period):
    # baseline="known": the streaming monitor is exactly the engine's detector (true baseline, restart at
    # every period boundary and at once after every alarm), so Theorem 1 applies to it as stated
    rng = np.random.default_rng(11)
    A = 40.0
    r = np.concatenate([rng.standard_normal(1500), 1.5 * rng.standard_normal(1500)])
    mon = Pm.PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=3.0, m0=M0, period=period)
    ev = [e.index for e in mon.run(r)]
    eng = Pm.Engine([Pm.Spec("k", "SR", A, period=period)], 1, B)
    y = np.minimum(r**2, 9.0) / M0
    hits = []
    for t in range(r.size):
        eng.boundary(t)
        if eng.step(y[t:t + 1])["k"][0]:
            hits.append(t)
    assert ev == hits and len(ev) > 5
    assert mon.segments[0]["start"] == 0 and mon.baseline == "known"  # no calibration window, no blind period
    assert all(seg["start"] == e + 1 for seg, e in zip(mon.segments[1:], ev))


def test_periodic_monitor_known_baseline_theorem1_on_calm_stream():
    # the referee's stream (Y = 1 for 15 periods, then Gaussian) through the streaming monitor: the plain
    # detector exceeds m / A in period 16, the restarted one does not (Theorem 1)
    rng = np.random.default_rng(12)
    P, m, A = 500, 10, 100.0
    hit_plain, hit_restart = np.zeros(P, bool), np.zeros(P, bool)
    for i in range(P):
        r = np.concatenate([np.full(15 * m, math.sqrt(M0)), rng.standard_normal(m)])
        a = Pm.PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=3.0, m0=M0).run(r)
        b = Pm.PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=3.0, m0=M0, period=m).run(r)
        assert all(e.index >= 15 * m for e in a) and all(e.index >= 15 * m for e in b)  # none while calm
        hit_plain[i], hit_restart[i] = len(a) > 0, len(b) > 0
    se = math.sqrt(0.1 * 0.9 / P)
    assert hit_restart.mean() <= m / A + 3 * se
    assert hit_plain.mean() > m / A + 5 * se


def test_periodic_monitor_known_baseline_errors():
    with pytest.raises(ValueError):
        Pm.PeriodicVolatilityMonitor(baseline="known", tau=3.0)  # m0 missing
    with pytest.raises(ValueError):
        Pm.PeriodicVolatilityMonitor(baseline="known", tau=1.0, m0=2.0)  # m0 above tau^2
    with pytest.raises(ValueError):
        Pm.PeriodicVolatilityMonitor(tau=3.0, m0=M0)  # tau and m0 need baseline="known"
    with pytest.raises(ValueError):
        Pm.PeriodicVolatilityMonitor(baseline="other")


def test_day_table():
    first = np.array([3, -1, 12])
    counts = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 2]])
    tab = Pm.day_table(first, counts)
    assert tab["p_first"].tolist() == pytest.approx([1 / 3, 0, 1 / 3])
    assert tab["at_risk"].tolist() == [3, 2, 2]
    assert tab["hazard"][2] == pytest.approx(0.5)
    assert tab["mean_count"].tolist() == pytest.approx([1 / 3, 0, 2 / 3])
    assert tab["cum_rate"][-1] == pytest.approx(1 / 3)
