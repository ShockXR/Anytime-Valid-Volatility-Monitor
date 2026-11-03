"""Independent audit tests (validation round, November 2026).

Boundary and invalid inputs, invariants of the e-process, causality (no use of data at or after the
monitored step), determinism of replays, and NaN / inf / empty-input behaviour. Every test runs on
small synthetic data with a fixed seed and the whole file takes well under a minute.
"""

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from avmonitor import VolatilityMonitor
from avmonitor import confseq as C
from avmonitor import detectors as D
from avmonitor import robust as R
from avmonitor.data import intraday_scale, load_bitstamp_csv, mark_frozen
from avmonitor.guard import spending_schedule

TAU = 2.5  # fixed truncation level for the null experiments, so that the null mean m0 is known exactly


def _m0_normal(tau):
    """E[min(Z^2, tau^2)] for standard normal Z, in closed form."""
    phi, big_phi = stats.norm.pdf(tau), stats.norm.cdf(tau)
    return (2 * big_phi - 1) - 2 * tau * phi + tau**2 * 2 * (1 - big_phi)


def _null_x(reps, n, seed, tau=TAU):
    z = np.random.default_rng(seed).standard_normal((reps, n))
    return np.minimum(z**2, tau**2)


# --------------------------------------------------------------------------------------------
# e-process invariants
# --------------------------------------------------------------------------------------------


def test_closed_form_null_mean_is_correct():
    x = _null_x(1, 1_500_000, seed=1)  # standard error about 0.001, tolerance 5 times that
    assert x.mean() == pytest.approx(_m0_normal(TAU), rel=5e-3)


@pytest.mark.parametrize("direction", ["up", "down"])
@pytest.mark.parametrize("kind", ["SR", "CUSUM"])
def test_statistic_is_nonnegative_for_any_admissible_input(direction, kind):
    rng = np.random.default_rng(3)
    # heavy tailed, includes the extreme values 0 and tau^2 where the bet factors are smallest
    z = rng.standard_t(2, size=(50, 300))
    x = np.minimum(z**2, TAU**2)
    x[:, ::17] = 0.0
    x[:, 5::17] = TAU**2
    m0 = _m0_normal(TAU)
    res = D.betting_edetector(x, m0, direction, upper=TAU**2, threshold=np.inf, kind=kind, return_path=True)
    assert np.all(np.isfinite(res.stat_path)) and np.all(res.stat_path >= 0)
    res = D.betting_edetector(x, m0, direction, upper=TAU**2, threshold=np.inf, kind=kind, adaptive_h=20,
                              return_path=True)
    assert np.all(res.stat_path >= 0)


def test_sr_path_matches_the_recursion_by_hand():
    x = _null_x(1, 60, seed=4)[0]
    m0 = _m0_normal(TAU)
    grid = np.array(D.DEFAULT_UP_GRID)
    S = np.zeros(len(grid))
    expected = []
    for v in x:
        S = (1 + grid * (v / m0 - 1)) * (S + 1.0)
        expected.append(S.mean())
    got = D.betting_edetector(x, m0, "up", threshold=np.inf, return_path=True).stat_path
    assert np.allclose(got, expected, rtol=1e-12, atol=0)


def test_cusum_never_exceeds_sr():
    x = _null_x(1, 500, seed=5)[0] * np.r_[np.ones(250), 2.0 * np.ones(250)]
    m0 = _m0_normal(TAU)
    sr = D.betting_edetector(x, m0, "up", threshold=np.inf, return_path=True).stat_path
    cu = D.betting_edetector(x, m0, "up", threshold=np.inf, kind="CUSUM", return_path=True).stat_path
    assert np.all(cu <= sr * (1 + 1e-12))


def test_null_mean_of_sr_statistic_is_t_with_the_true_baseline():
    # E[M_t] = t exactly when E[x] = m0, because every bet factor has conditional mean 1 and M_t = f_t (M_{t-1} + 1).
    reps, n = 40_000, 12
    x = _null_x(reps, n, seed=6)
    m0 = _m0_normal(TAU)
    grid = np.array(D.DEFAULT_UP_GRID)
    S = np.zeros((reps, len(grid)))
    for t in range(n):
        S = (1 + grid[None, :] * (x[:, t : t + 1] / m0 - 1)) * (S + 1.0)
        assert S.mean() == pytest.approx(t + 1, rel=0.01)
    # each single bet factor has mean one under the null
    for lam in grid:
        assert (1 + lam * (x / m0 - 1)).mean() == pytest.approx(1.0, abs=5e-3)


def test_false_alarm_rate_respects_the_bound_m_over_a():
    reps, n = 3000, 100
    x = _null_x(reps, n, seed=7)
    m0 = _m0_normal(TAU)
    for A in (400.0, 2000.0):  # bound m / A = 0.25 and 0.05; each side runs at level 2A
        up = D.betting_edetector(x, m0, "up", threshold=2 * A)
        dn = D.betting_edetector(x, m0, "down", upper=TAU**2, threshold=2 * A)
        any_alarm = (up.alarm >= 0) | (dn.alarm >= 0)
        assert any_alarm.mean() <= n / A + 3 * math.sqrt(0.25 / reps)


def test_alarm_is_the_first_index_with_statistic_at_least_the_threshold():
    x = _null_x(1, 400, seed=8)[0] * np.r_[np.ones(200), 3.0 * np.ones(200)]
    m0 = _m0_normal(TAU)
    path = D.betting_edetector(x, m0, "up", threshold=np.inf, return_path=True).stat_path
    for k in (50, 250, 390):
        thr = float(path[k])  # exactly equal to the statistic: >= must fire there (not >)
        res = D.betting_edetector(x, m0, "up", threshold=thr)
        first = int(np.argmax(path >= thr))
        assert res.alarm[0] == first <= k
        assert np.all(path[:first] < thr)
    # one ulp above the maximum: never fires
    assert D.betting_edetector(x, m0, "up", threshold=float(np.nextafter(path.max(), np.inf))).alarm[0] == -1


def test_monitor_alarm_level_is_exactly_two_a():
    rng = np.random.default_rng(9)
    A = 50.0
    r = np.r_[rng.standard_normal(300), 3.0 * rng.standard_normal(300)]
    m = VolatilityMonitor(threshold=A, cal_len=300)
    ev = None
    for v in r:
        ev = m.update(v) or ev
        if ev is not None:
            break
    assert ev is not None and ev.direction == "up"
    assert ev.stat >= 2 * A


def test_down_detector_rejects_input_above_its_upper_bound():
    with pytest.raises(ValueError, match="upper bound"):
        D.betting_edetector(np.array([0.5, 2.0]), 1.0, "down", upper=1.5)
    with pytest.raises(ValueError, match="must exceed m0"):
        D.betting_edetector(np.array([0.5]), 1.0, "down", upper=1.0)


def test_bets_outside_the_admissible_range_are_rejected():
    x = np.ones(10)
    with pytest.raises(ValueError):
        D.betting_edetector(x, 1.0, "up", grid=[0.2, 1.0])  # lam = 1 makes the factor 0 at x = 0
    with pytest.raises(ValueError):
        D.betting_edetector(x, 1.0, "up", grid=[-0.1, 0.2])
    with pytest.raises(ValueError):
        D.check_bets([0.5], [])
    with pytest.raises(ValueError):
        D.check_bets([0.5], [np.nan])


# --------------------------------------------------------------------------------------------
# truncation
# --------------------------------------------------------------------------------------------


def test_truncation_never_exceeds_its_cap():
    rng = np.random.default_rng(10)
    r = rng.standard_t(1.5, size=5000) * 1e-3
    r[::500] = np.array([np.inf, -np.inf, 1e300, -1e300, 5.0, -5.0, 0.0, 1e-300, 7.0, -7.0])
    tau = R.mad_threshold(r[np.isfinite(r)], c=3.0, scaled=True)
    out = R.truncate(r, tau)
    assert np.all(np.abs(out) <= tau)
    assert np.all(np.minimum(r**2, tau**2) <= tau**2)
    assert np.isnan(R.truncate(np.array([np.nan]), 1.0)[0])


def test_monitor_statistic_input_is_capped_at_tau_squared():
    rng = np.random.default_rng(11)
    m = VolatilityMonitor(threshold=1e12, cal_len=200)
    for v in rng.standard_normal(200):
        m.update(v)
    assert m.calibrated
    tau2 = m.tau2
    s_before = m.S_up.copy()
    m.update(1e6)  # absurdly large return
    capped = (1 + m.up_grid * (tau2 / m.m_up - 1.0)) * (s_before + 1.0)
    assert np.allclose(m.S_up, capped, rtol=1e-12)


def test_mad_scale_and_tau_follow_the_definition():
    rng = np.random.default_rng(12)
    r = rng.standard_normal(1001)
    b = D.calibrate(r, c=3.0)
    assert b.scale == pytest.approx(np.median(np.abs(r)) * R.MAD_NORMAL_SCALE)
    assert b.tau == pytest.approx(3.0 * b.scale)
    assert b.m0 == pytest.approx(np.mean(np.minimum(r**2, b.tau**2)))
    wide = D.calibrate(r, c=3.0, margin_z=2.0, null_band=0.1)
    assert wide.m0_down <= b.m0 <= wide.m0_up


# --------------------------------------------------------------------------------------------
# causality: nothing at or after the monitored step may influence it
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("baseline", ["plugin", "cs"])
def test_monitor_events_up_to_t_do_not_depend_on_later_returns(baseline):
    rng = np.random.default_rng(13)
    r = np.r_[rng.standard_normal(250), 2.5 * rng.standard_normal(250), rng.standard_normal(250),
              0.3 * rng.standard_normal(250)]
    cut = 620
    r2 = r.copy()
    r2[cut + 1 :] = rng.standard_normal(r.size - cut - 1) * 10
    kw = dict(threshold=30.0, cal_len=200, baseline=baseline)
    a, b = VolatilityMonitor(**kw), VolatilityMonitor(**kw)
    a.run(r)
    b.run(r2)
    ea = [(e.index, e.direction, e.monitored_since, e.stat) for e in a.events if e.index <= cut]
    eb = [(e.index, e.direction, e.monitored_since, e.stat) for e in b.events if e.index <= cut]
    assert ea == eb
    assert len(ea) >= 1  # the test has teeth: at least one alarm happens before the cut


def test_calibration_window_is_strictly_before_the_monitored_steps():
    rng = np.random.default_rng(14)
    r = np.r_[rng.standard_normal(300), 4 * rng.standard_normal(600)]
    m = VolatilityMonitor(threshold=20.0, cal_len=250)
    m.run(r)
    assert m.segments[0]["start"] == 250  # first monitored index follows the 250 calibration returns
    prev_end = -1
    for seg in m.segments:
        assert seg["start"] > prev_end
        prev_end = seg["start"]
    for ev in m.events:
        assert ev.index >= ev.monitored_since
    # A change placed exactly in the calibration window cannot be alarmed on, only absorbed into the baseline
    m2 = VolatilityMonitor(threshold=20.0, cal_len=250)
    m2.run(r[:250] * 5)
    assert m2.events == []


def test_array_detector_with_adaptive_bet_is_predictable():
    x = _null_x(1, 300, seed=15)[0]
    m0 = _m0_normal(TAU)
    t0 = 150
    x2 = x.copy()
    x2[t0:] *= 4.0
    kw = dict(direction="up", threshold=np.inf, adaptive_h=20, return_path=True)
    a = D.betting_edetector(x, m0, **kw)
    b = D.betting_edetector(x2, m0, **kw)
    assert np.array_equal(a.stat_path[:t0], b.stat_path[:t0])
    # the bet used AT step t0 may only use data before t0
    assert a.lam_path[t0] == b.lam_path[t0]
    assert a.lam_path[t0 + 1] != b.lam_path[t0 + 1]


def test_intraday_scale_uses_only_earlier_days():
    rng = np.random.default_rng(16)
    days = pd.date_range("2024-01-01", periods=40, freq="D", tz="UTC")
    idx = pd.DatetimeIndex([d + pd.Timedelta(minutes=i) for d in days for i in range(120)])
    r = pd.Series(rng.standard_normal(len(idx)) * 1e-3, index=idx)
    base = intraday_scale(r, lookback_days=10, smooth=5, circular=False)
    d0 = days[30]
    r2 = r.copy()
    r2[r2.index >= d0] *= 7.0
    alt = intraday_scale(r2, lookback_days=10, smooth=5, circular=False)
    upto = r.index < d0 + pd.Timedelta(days=1)  # day 30 itself must be unchanged too (current day excluded)
    assert np.allclose(base[upto].to_numpy(), alt[upto].to_numpy(), equal_nan=True)
    assert not np.allclose(base[~upto].to_numpy(), alt[~upto].to_numpy(), equal_nan=True)
    # without enough history the scale is NaN rather than borrowed from the future
    assert base[r.index < days[5]].isna().all()


def test_ew_quantile_scale_is_predictable_when_start_value_is_given():
    rng = np.random.default_rng(17)
    r = rng.standard_normal(400) * 1e-3
    r2 = r.copy()
    r2[200:] *= 5
    a = R.ew_quantile_scale(r, halflife=50, init=1e-3)
    b = R.ew_quantile_scale(r2, halflife=50, init=1e-3)
    assert np.array_equal(a[:201], b[:201])  # entry t uses data before t only
    s1 = R.ew_quantile_scale(r, halflife=50, strict=True)
    assert np.isnan(s1[:50]).all() and np.isfinite(s1[50:]).all()


def test_ew_quantile_scale_default_start_value_looks_ahead_as_documented():
    # Documented start-up look-ahead (init=None): entries 0 to 49 depend on the first 50 returns.
    rng = np.random.default_rng(18)
    r = rng.standard_normal(200) * 1e-3
    r2 = r.copy()
    r2[40:50] *= 100
    assert R.ew_quantile_scale(r, 50)[0] != R.ew_quantile_scale(r2, 50)[0]
    assert R.ew_quantile_scale(r, 50, init=1e-3)[0] == R.ew_quantile_scale(r2, 50, init=1e-3)[0]


def test_rolling_mad_zero_excludes_the_current_step():
    rng = np.random.default_rng(19)
    r = rng.standard_normal(100)
    r2 = r.copy()
    r2[60] = 1e6
    a, b = R.rolling_mad_zero(r, 20), R.rolling_mad_zero(r2, 20)
    assert np.array_equal(a[:61], b[:61], equal_nan=True)
    assert np.isnan(a[:20]).all()


# --------------------------------------------------------------------------------------------
# determinism and state
# --------------------------------------------------------------------------------------------


def _events(seed, **kw):
    rng = np.random.default_rng(seed)
    r = np.r_[rng.standard_normal(400), 2.0 * rng.standard_normal(400), 0.4 * rng.standard_normal(600)]
    m = VolatilityMonitor(threshold=25.0, cal_len=300, **kw)
    m.run(r)
    return [(e.index, e.direction, e.stat) for e in m.events], m.segments


def test_replay_with_same_seed_is_bit_identical_and_other_seed_differs():
    for kw in ({}, {"baseline": "cs"}):
        a, sa = _events(21, **kw)
        b, sb = _events(21, **kw)
        assert a == b and sa == sb
        c, _ = _events(22, **kw)
        assert a != c


def test_streaming_updates_equal_batch_run_and_chunked_run():
    rng = np.random.default_rng(23)
    r = np.r_[rng.standard_normal(300), 3 * rng.standard_normal(300)]
    one = VolatilityMonitor(threshold=15.0, cal_len=200)
    one.run(r)
    two = VolatilityMonitor(threshold=15.0, cal_len=200)
    for v in r:
        two.update(v)
    three = VolatilityMonitor(threshold=15.0, cal_len=200)
    for chunk in np.array_split(r, 7):
        three.run(chunk)
    key = lambda m: [(e.index, e.direction, e.stat) for e in m.events]  # noqa: E731
    assert key(one) == key(two) == key(three) and len(key(one)) >= 1


def test_two_monitors_do_not_share_state():
    rng = np.random.default_rng(24)
    a, b = VolatilityMonitor(threshold=20.0, cal_len=100), VolatilityMonitor(threshold=20.0, cal_len=100)
    a.run(rng.standard_normal(150))
    assert b.t == -1 and b.events == [] and b.segments == [] and not b.calibrated
    assert VolatilityMonitor().up_grid[0] == D.DEFAULT_UP_GRID[0]  # constructing one monitor never edits the defaults
    VolatilityMonitor().up_grid[0] = 0.9
    assert VolatilityMonitor().up_grid[0] == D.DEFAULT_UP_GRID[0]


def test_reset_after_alarm_clears_the_statistics():
    rng = np.random.default_rng(25)
    m = VolatilityMonitor(threshold=10.0, cal_len=100)
    for v in np.r_[rng.standard_normal(150), 6 * rng.standard_normal(100)]:
        if m.update(v) is not None:
            break
    assert len(m.events) == 1 and not m.calibrated and m.buffer == [] and not m.S_up.any() and not m.S_dn.any()


def test_spending_schedule_sums_to_delta():
    for K in (1, 2, 3, 7):
        for kind in ("halving", "equal"):
            assert math.isclose(sum(spending_schedule(0.05, K, kind)), 0.05, rel_tol=1e-12)
    with pytest.raises(ValueError):
        spending_schedule(0.05, 0)
    with pytest.raises(ValueError):
        spending_schedule(1.5, 2)


# --------------------------------------------------------------------------------------------
# confidence sequence
# --------------------------------------------------------------------------------------------


def test_confidence_sequence_covers_and_streaming_equals_vectorised():
    rng = np.random.default_rng(26)
    reps, n, bound = 1500, 200, 4.0
    x = np.minimum(rng.exponential(1.0, (reps, n)), bound)
    true_mean = 1.0 - math.exp(-bound)  # E[min(Exp(1), 4)]
    cs = C.eb_confseq(x, bound, delta_upper=0.025, delta_lower=0.025)
    miss = ((cs.lower > true_mean) | (cs.upper < true_mean)).any(axis=1).mean()
    assert miss <= 0.05 + 3 * math.sqrt(0.05 * 0.95 / reps)
    assert np.all(cs.lower <= cs.upper) and np.all(cs.lower >= 0) and np.all(cs.upper <= bound)
    s = C.StreamingEBConfSeq(bound, 0.025, 0.025)
    for t in range(n):
        s.update(x[0, t])
        assert s.upper == pytest.approx(cs.upper[0, t], rel=1e-9)
        assert s.lower == pytest.approx(cs.lower[0, t], rel=1e-9)


def test_confidence_sequence_input_checks():
    with pytest.raises(ValueError):
        C.eb_confseq(np.array([0.5, np.nan]), 1.0)
    with pytest.raises(ValueError):
        C.eb_confseq(np.array([0.5, 1.5]), 1.0)
    with pytest.raises(ValueError):
        C.eb_confseq(np.array([0.5]), 1.0, delta_upper=0.0)
    with pytest.raises(ValueError):
        C.StreamingEBConfSeq(0.0)
    s = C.StreamingEBConfSeq(2.0)
    assert s.upper == 2.0 and s.lower == 0.0 and math.isnan(s.mean) and s.rel_width == math.inf
    with pytest.raises(ValueError):
        s.update(-0.1)


# --------------------------------------------------------------------------------------------
# NaN, inf, empty and degenerate input
# --------------------------------------------------------------------------------------------


def test_calibrate_rejects_empty_nan_and_degenerate_windows():
    for bad in (np.array([]), np.full(30, np.nan), np.array([0.01])):
        with pytest.raises(ValueError):
            D.calibrate(bad)
    with pytest.raises(ValueError, match="degenerate"):
        D.calibrate(np.zeros(100))
    with pytest.raises(ValueError):
        D.calibrate(np.ones(50), baseline="bogus")
    b = D.calibrate(np.r_[np.random.default_rng(27).standard_normal(100), np.nan, np.inf])
    assert np.isfinite(b.m0) and b.n_cal == 100  # NaN and inf are dropped before calibrating


def test_calibrate_with_infinite_return_drops_it():
    r = np.random.default_rng(28).standard_normal(200)
    with_inf = np.r_[r, np.inf]
    assert D.calibrate(with_inf).m0 == pytest.approx(D.calibrate(r).m0)


def test_nan_return_leaves_the_monitor_statistic_unchanged():
    rng = np.random.default_rng(29)
    m = VolatilityMonitor(threshold=1e9, cal_len=100)
    for v in rng.standard_normal(120):
        m.update(v)
    s_up, s_dn = m.S_up.copy(), m.S_dn.copy()
    assert m.update(np.nan) is None and m.update(np.inf) is None and m.update(-np.inf) is None
    assert np.array_equal(m.S_up, s_up) and np.array_equal(m.S_dn, s_dn)


def test_array_detector_skips_nan_steps_without_adding_to_the_statistic():
    x = _null_x(1, 100, seed=30)[0]
    x_nan = x.copy()
    x_nan[40] = np.nan
    m0 = _m0_normal(TAU)
    a = D.betting_edetector(np.delete(x, 40), m0, "up", threshold=np.inf, return_path=True).stat_path
    b = D.betting_edetector(x_nan, m0, "up", threshold=np.inf, return_path=True).stat_path
    assert np.allclose(np.delete(b, 40), a, rtol=1e-12)
    assert b[40] == b[39]


def test_empty_and_short_inputs_do_not_crash_the_detectors():
    m0 = 1.0
    assert D.betting_edetector(np.empty((3, 0)), m0).alarm.tolist() == [-1, -1, -1]
    assert D.betting_edetector(np.ones(5), m0, start=50).alarm.tolist() == [-1]
    assert D.gaussian_lr_edetector(np.empty((2, 0)), 1.0).alarm.tolist() == [-1, -1]
    assert D.page_cusum(np.empty((2, 0)), 1.0, 0.0, 5.0).alarm.tolist() == [-1, -1]
    assert D.repeated_window_ztest(np.ones((2, 10)), 1.0, 2.0, window=60).alarm.tolist() == [-1, -1]
    out = D.summarize_alarms(np.array([], dtype=int), None, 100)
    assert out["false_alarm_rate"] != out["false_alarm_rate"]  # NaN for an empty set of paths, not an exception


def test_monitor_survives_zero_and_all_nan_streams():
    for stream in (np.zeros(400), np.full(400, np.nan)):
        m = VolatilityMonitor(threshold=100.0, cal_len=50)
        m.run(stream)  # no exception, no calibration, no alarm
        assert len(m.events) == 0 and not m.calibrated


def test_constant_stream_statistic_equals_step_count():
    # With zero variance every bet factor is exactly 1, so M_t = t: the Shiryaev-Roberts statistic grows by one per
    # step and the first alarm is at exactly 2A monitored steps (the average run length bound is attained).
    m = VolatilityMonitor(threshold=100.0, cal_len=50)
    m.run(np.full(400, 1e-3))
    assert [(e.index, e.direction) for e in m.events[:1]] == [(50 + 200 - 1, "up")]


def test_summarize_alarms_delay_arithmetic():
    out = D.summarize_alarms(np.array([-1, 5, 12, 30]), nu=10, horizon=100)
    assert out["false_alarm_before_change"] == 0.25
    assert out["detected"] == 0.5
    assert out["mean_delay"] == pytest.approx(((12 - 10 + 1) + (30 - 10 + 1)) / 2)
    assert D.summarize_alarms(np.array([-1, 3]), None, 50)["mean_run_length_censored"] == (50 + 4) / 2


def test_two_sided_alarm_takes_the_earlier_side():
    up = D.DetectorResult(np.array([5, -1, -1, 9]))
    dn = D.DetectorResult(np.array([7, 4, -1, -1]))
    assert D.two_sided_alarm(up, dn).tolist() == [5, 4, -1, 9]


# --------------------------------------------------------------------------------------------
# data helpers
# --------------------------------------------------------------------------------------------


def test_mark_frozen_run_length_boundary():
    r = pd.Series([0.1, 0, 0, 0, 0.2, 0, 0, 0, 0, 0.3] + [np.nan, 0.0, 0.0, 0.0])
    out = mark_frozen(r, min_run=4)
    assert out.iloc[1:4].notna().all()  # run of 3 is kept
    assert out.iloc[5:9].isna().all()  # run of 4 is blanked
    assert out.iloc[0] == 0.1 and out.iloc[9] == 0.3
    assert out.iloc[-3:].notna().all()  # trailing run of 3 after a NaN is kept


def test_bitstamp_loader_filled_minutes_and_window(tmp_path):
    ts = 1_700_000_000 - 1_700_000_000 % 60
    rows = [(ts + 60 * i, c, 1.0 if i != 3 else 0.0) for i, c in enumerate([100, 101, 102, 102, 104, 105])]
    f = tmp_path / "b.csv"
    pd.DataFrame(rows, columns=["timestamp", "close", "volume"]).to_csv(f, index=False)
    df = load_bitstamp_csv(f)
    assert df["logret"].iloc[0] != df["logret"].iloc[0]  # first return is NaN
    assert df["logret"].iloc[3] != df["logret"].iloc[3] and df["logret"].iloc[4] != df["logret"].iloc[4]
    assert df["logret"].iloc[[1, 2, 5]].notna().all()
    z = load_bitstamp_csv(f, filled="zero")
    assert z["logret"].iloc[3] == 0.0
    sub = load_bitstamp_csv(f, start=pd.Timestamp(ts + 120, unit="s", tz="UTC").isoformat(),
                            end=pd.Timestamp(ts + 240, unit="s", tz="UTC").isoformat())
    assert len(sub) == 2  # start inclusive, end exclusive
    with pytest.raises(ValueError):
        load_bitstamp_csv(f, filled="nope")


# --------------------------------------------------------------------------------------------
# invalid constructor arguments (found by the audit, see reports/CODE_AUDIT.md)
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("bad", [float("nan"), 0.0, -5.0])
def test_monitor_rejects_nonpositive_or_nan_threshold(bad):
    with pytest.raises(ValueError, match="threshold"):
        VolatilityMonitor(threshold=bad, cal_len=50)


@pytest.mark.parametrize("bad", [float("nan"), 0.0, -1.0, float("inf")])
def test_monitor_rejects_invalid_truncation_constant(bad):
    with pytest.raises(ValueError):
        VolatilityMonitor(c=bad, cal_len=50)


def test_monitor_rejects_calibration_window_below_two():
    for bad in (0, 1):
        with pytest.raises(ValueError, match="cal_len"):
            VolatilityMonitor(cal_len=bad)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan")])
def test_array_detector_rejects_invalid_null_mean(bad):
    with pytest.raises(ValueError, match="m0"):
        D.betting_edetector(np.ones(20), bad)


def test_array_detector_rejects_negative_statistic():
    x = np.ones(20)
    x[7] = -0.5
    with pytest.raises(ValueError, match="nonnegative"):
        D.betting_edetector(x, 1.0, "up")


def test_array_detector_rejects_unknown_kind_even_for_empty_input():
    with pytest.raises(ValueError, match="kind"):
        D.betting_edetector(np.empty((1, 0)), 1.0, kind="BOGUS")


def test_benchmark_detectors_go_dead_after_one_nan_documented():
    # Documented limitation (reports/CODE_AUDIT.md, A1-B-008): unlike betting_edetector and VolatilityMonitor, the two
    # benchmarks do not skip missing values, so a single NaN freezes their statistic at NaN and no alarm can follow.
    rng = np.random.default_rng(0)
    r = np.r_[rng.standard_normal(100), 3 * rng.standard_normal(200)]
    r_nan = r.copy()
    r_nan[5] = np.nan
    assert D.gaussian_lr_edetector(r, 1.0, threshold=100).alarm[0] >= 100
    assert D.page_cusum(r**2, 1.0, 0.5, 20.0).alarm[0] >= 100
    assert D.gaussian_lr_edetector(r_nan, 1.0, threshold=100).alarm[0] == -1
    assert D.page_cusum(r_nan**2, 1.0, 0.5, 20.0).alarm[0] == -1
