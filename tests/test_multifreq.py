"""Tests for the multi-frequency monitor (``avmonitor.multifreq``, ``docs/theory_multifreq.md``)."""

import math

import numpy as np
import pytest

from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID, betting_edetector
from avmonitor.monitor import VolatilityMonitor
from avmonitor.multifreq import (
    MultiFrequencyMonitor,
    block_returns,
    calibration_block_returns,
    component_paths,
    first_alarm,
    handicap,
    level_for_budget,
    minute_values,
    plugin_calibration,
    proved_bound,
    resolve_weights,
)
from avmonitor.robust import MAD_NORMAL_SCALE

FREQS = (1, 5, 15)


def true_baseline_gauss(k, gap=False, n=400_000, c=3.0, seed=0):
    """tau2 and m0 of the k-minute block return of i.i.d. N(0, 1) one-minute returns (Monte Carlo)."""
    L = k - 1 if (gap and k > 1) else k
    R = np.random.default_rng(seed).standard_normal(n) * math.sqrt(L)
    tau2 = (c * np.median(np.abs(R)) * MAD_NORMAL_SCALE) ** 2
    return tau2, float(np.mean(np.minimum(R**2, tau2)))


KNOWN = {k: true_baseline_gauss(k, seed=k) for k in FREQS}
KNOWN_GAP = {k: true_baseline_gauss(k, gap=True, seed=k) for k in FREQS}


def vectorised_first_alarm(r, A, freqs, weights, known, combine="average", certified=True, gap=False):
    """First alarm (minute within the monitored segment, direction) from the vectorised building blocks."""
    r = np.atleast_2d(r)
    n = r.shape[1]
    w = resolve_weights(weights, freqs)
    up = np.zeros_like(r) if combine == "average" else np.full_like(r, -np.inf)
    dn = up.copy()
    for wk, k in zip(w, freqs):
        tau2, m0 = known[k]
        x = np.minimum(block_returns(r, k, gap) ** 2, tau2)
        p = component_paths(x, m0, m0, tau2)
        vu = minute_values(p["up_full"], p["up_cert"], k, n, certified)
        vd = minute_values(p["dn_full"], p["dn_cert"], k, n, certified)
        if combine == "average":
            up += wk * vu
            dn += wk * vd
        else:
            up = np.maximum(up, wk * vu)
            dn = np.maximum(dn, wk * vd)
    return first_alarm(up, dn, A)


# ----------------------------------------------------------------------------
# Weights, bounds, levels
# ----------------------------------------------------------------------------


def test_weights_bounds_and_handicap():
    assert np.allclose(resolve_weights("balanced", FREQS), np.array([1, 5, 15]) / 21)
    assert np.allclose(resolve_weights("equal", FREQS), 1 / 3)
    assert np.allclose(resolve_weights("sqrt", FREQS).sum(), 1)
    # balanced weights: every frequency has handicap J = 3
    assert np.allclose(handicap(FREQS, "balanced"), 3.0)
    assert np.allclose(handicap(FREQS, "balanced", m=1440), 3.0)
    # equal weights: k * sum(w / j) / w_k = k * (1 + 1/5 + 1/15)
    assert np.allclose(handicap(FREQS, "equal"), np.array([1, 5, 15]) * (1 + 1 / 5 + 1 / 15))
    # the worst handicap is at least J for any weights (a weighted average of the handicaps equals J)
    rng = np.random.default_rng(1)
    for _ in range(50):
        assert handicap(FREQS, rng.random(3) + 1e-3).max() >= 3 - 1e-9
    # single frequency: the bound is the block-monitor bound ceil(m / k) / A
    assert proved_bound(1440, 1920.0, (15,), [1.0]) == pytest.approx(0.05)
    assert proved_bound(1440, 28800.0, (1,), [1.0]) == pytest.approx(0.05)
    # certified average with balanced weights: 28800 * 3 / 21
    A = level_for_budget(0.05, 1440, FREQS, "balanced")
    assert A == pytest.approx(28800 * 3 / 21)
    assert proved_bound(1440, A, FREQS, "balanced") == pytest.approx(0.05)
    # the weighted maximum has the same bound at m divisible by every k; the stale average a larger one
    assert proved_bound(1440, A, FREQS, "balanced", combine="max") == pytest.approx(0.05)
    assert proved_bound(1440, A, FREQS, "balanced", check="stale") == pytest.approx(0.05 * 3)
    # with the same A as the one-minute monitor, any weight on coarse frequencies lowers the bound
    assert proved_bound(1440, 28800.0, FREQS, "equal") < 0.05


# ----------------------------------------------------------------------------
# k = 1 alone reproduces VolatilityMonitor
# ----------------------------------------------------------------------------


def _stream(seed, n=12000):
    rng = np.random.default_rng(seed)
    sd = np.ones(n)
    sd[3000:5000] = 2.0
    sd[7000:9000] = 0.5
    r = sd * rng.standard_normal(n)
    r[rng.random(n) < 0.01] = np.nan
    r[10000:10250] = 0.0  # degenerate stretch: forces extended windows and possibly the fallback
    r[10500:10520] = 0.0
    return r


@pytest.mark.parametrize("band", [0.0, 0.25])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_k1_reproduces_volatility_monitor(band, seed):
    r = _stream(seed)
    ref = VolatilityMonitor(threshold=300.0, cal_len=200, null_band=band)
    mf = MultiFrequencyMonitor(threshold=300.0, freqs=(1,), weights=[1.0], cal_len=200, null_band=band)
    ev_ref = ref.run(r)
    ev_mf = mf.run(r)
    assert len(ev_ref) > 5
    assert [(e.index, e.direction, e.monitored_since) for e in ev_ref] == [
        (e.index, e.direction, e.monitored_since) for e in ev_mf]
    assert [e.stat for e in ev_ref] == [e.stat for e in ev_mf]
    assert ref.cal_stats == mf.cal_stats
    assert ref.cal_stats["extended"] > 0


def test_single_coarse_frequency_certified_equals_stale():
    r = _stream(3)
    a = MultiFrequencyMonitor(threshold=60.0, freqs=(5,), weights=[1.0], cal_len=300).run(r)
    b = MultiFrequencyMonitor(threshold=60.0, freqs=(5,), weights=[1.0], cal_len=300, check="stale").run(r)
    assert len(a) > 3
    assert [(e.index, e.direction) for e in a] == [(e.index, e.direction) for e in b]
    assert all((e.index - e.monitored_since + 1) % 5 == 0 for e in a)  # a lone 5-minute detector alarms at block ends


# ----------------------------------------------------------------------------
# The weighted-average recursion matches separate detectors
# ----------------------------------------------------------------------------


def test_components_match_separate_detectors():
    rng = np.random.default_rng(5)
    n = 900
    r = rng.standard_normal(n) * np.where(np.arange(n) >= 450, 1.6, 1.0)
    for k in FREQS:
        tau2, m0 = KNOWN[k]
        x = np.minimum(block_returns(r, k) ** 2, tau2)
        p = component_paths(x[None, :], m0, m0, tau2, cap=np.inf)
        up = betting_edetector(x, m0, "up", threshold=np.inf, grid=DEFAULT_UP_GRID, return_path=True)
        dn = betting_edetector(x, m0, "down", upper=tau2, threshold=np.inf, grid=DEFAULT_DOWN_FRACTIONS,
                               return_path=True)
        assert np.allclose(p["up_full"][0], up.stat_path, rtol=1e-12)
        assert np.allclose(p["dn_full"][0], dn.stat_path, rtol=1e-12)


def test_streaming_average_equals_weighted_sum_of_components():
    rng = np.random.default_rng(6)
    cal, n = 100, 900
    r = rng.standard_normal(cal + n) * np.where(np.arange(cal + n) >= cal + 600, 1.3, 1.0)
    r[cal + 37] = np.nan  # one missing return: skipped by the 1-minute detector and by the blocks containing it
    for check in ("certified", "stale"):
        mon = MultiFrequencyMonitor(threshold=1e12, freqs=FREQS, weights="sqrt", cal_len=cal, check=check,
                                    known_baseline=KNOWN, trace=True)
        mon.run(r)
        mu = np.array([row[1] for row in mon.trace_rows])
        md = np.array([row[2] for row in mon.trace_rows])
        w = resolve_weights("sqrt", FREQS)
        exp_u = np.zeros(n)
        exp_d = np.zeros(n)
        for wk, k in zip(w, FREQS):
            tau2, m0 = KNOWN[k]
            x = np.minimum(block_returns(r[cal:], k) ** 2, tau2)
            p = component_paths(x[None, :], m0, m0, tau2, cap=np.inf)
            exp_u += wk * minute_values(p["up_full"], p["up_cert"], k, n, check == "certified")[0]
            exp_d += wk * minute_values(p["dn_full"], p["dn_cert"], k, n, check == "certified")[0]
            if k > 1:  # at block ends the certified value is the detector's own full statistic
                ends = np.arange(k - 1, n, k)
                assert np.array_equal(minute_values(p["up_full"], p["up_cert"], k, n)[0, ends], p["up_full"][0])
        assert np.allclose(mu, exp_u, rtol=1e-10)
        assert np.allclose(md, exp_d, rtol=1e-10)
        # at common block ends (multiples of 15 minutes) the certified and the plain average coincide
        if check == "certified":
            stale = MultiFrequencyMonitor(threshold=1e12, freqs=FREQS, weights="sqrt", cal_len=cal, check="stale",
                                          known_baseline=KNOWN, trace=True)
            stale.run(r)
            su = np.array([row[1] for row in stale.trace_rows])
            common = np.arange(14, n, 15)
            assert np.allclose(mu[common], su[common])
            assert np.all(mu <= su + 1e-9)


@pytest.mark.parametrize("gap", [False, True])
@pytest.mark.parametrize("combine,check", [("average", "certified"), ("average", "stale"), ("max", "stale")])
def test_vectorised_matches_streaming(combine, check, gap):
    rng = np.random.default_rng(7)
    cal, n = 60, 1500
    known = KNOWN_GAP if gap else KNOWN
    for rep in range(6):
        ratio = [1.0, 1.5, 0.6][rep % 3]
        r = rng.standard_normal(cal + n) * np.where(np.arange(cal + n) >= cal + 400, ratio, 1.0)
        mon = MultiFrequencyMonitor(threshold=200.0, freqs=FREQS, weights="balanced", cal_len=cal, combine=combine,
                                    check=check, gap=gap, known_baseline=known)
        mon.run(r)
        idx, d = vectorised_first_alarm(r[cal:], 200.0, FREQS, "balanced", known, combine, check == "certified", gap)
        if mon.events:
            e = mon.events[0]
            assert idx[0] == e.index - cal and d[0] == (1 if e.direction == "up" else -1)
        else:
            assert idx[0] == -1


# ----------------------------------------------------------------------------
# The bound on long i.i.d. null streams
# ----------------------------------------------------------------------------


def test_bound_first_alarm_iid_null():
    rng = np.random.default_rng(8)
    m, A = 300, 300.0  # every bound below is under 1, so the comparison is informative
    r = rng.standard_normal((4000, m))
    for weights in ("balanced", "equal"):
        for combine, check in (("average", "certified"), ("max", "stale"), ("average", "stale")):
            idx, _ = vectorised_first_alarm(r, A, FREQS, weights, KNOWN, combine, check == "certified")
            rate = np.mean(idx >= 0)
            bound = proved_bound(m, A, FREQS, weights, combine, check)
            assert 0 < rate <= bound, (weights, combine, check, rate, bound)
    # one frequency alone at its own level: the block monitor's m / (k A) bound
    for k in FREQS:
        idx, _ = vectorised_first_alarm(r, A / k, (k,), [1.0], KNOWN)
        assert 0 < np.mean(idx >= 0) <= proved_bound(m, A / k, (k,), [1.0])


def test_bound_with_restarts_long_stream():
    """Corollary 4: expected alarms per monitored minute at most about sum_k w_k / k / A (true baseline)."""
    rng = np.random.default_rng(9)
    n, cal, A = 120_000, 60, 25.0
    r = rng.standard_normal(n)
    mon = MultiFrequencyMonitor(threshold=A, freqs=FREQS, weights="balanced", cal_len=cal, known_baseline=KNOWN)
    ev = mon.run(r)
    monitored = n - cal * len(mon.segments)
    w = resolve_weights("balanced", FREQS)
    kappa = float(np.sum(w * (1 - 1 / np.array(FREQS))))
    bound = (monitored * float(np.sum(w / np.array(FREQS))) + kappa) / (A - kappa)
    assert 10 < len(ev) <= bound


# ----------------------------------------------------------------------------
# Why the frozen coarse statistic is not a one-minute supermartingale
# ----------------------------------------------------------------------------


def test_frozen_two_minute_statistic_violates_optional_stopping():
    """Theory note, Proposition 2: stop in the middle of block 2 when |r_3| <= 1, else at its end.

    With i.i.d. N(0, 1) returns, no truncation and one bet lambda on R^2 / 2, E[M_tau] exceeds the expected number
    of completed blocks by 0.484 lambda (exact value lambda E[(Z^2 - 1) 1{|Z| > 1}] = 0.48394 lambda).
    The tolerance is about 4 Monte Carlo standard errors.
    """
    rng = np.random.default_rng(10)
    lam, n = 0.5, 400_000
    z = rng.standard_normal((n, 4))
    f1 = 1 + lam * ((z[:, 0] + z[:, 1]) ** 2 / 2 - 1)
    f2 = 1 + lam * ((z[:, 2] + z[:, 3]) ** 2 / 2 - 1)
    stop_mid = np.abs(z[:, 2]) <= 1
    M_tau = np.where(stop_mid, f1, f2 * (f1 + 1))
    blocks_done = np.where(stop_mid, 1, 2)
    gap = M_tau.mean() - blocks_done.mean()
    assert gap == pytest.approx(0.48394 * lam, abs=0.01)


# ----------------------------------------------------------------------------
# Blocks with a gap under price noise
# ----------------------------------------------------------------------------


def test_gap_blocks_remove_noise_dependence():
    rng = np.random.default_rng(11)
    n, gamma, k = 600_000, 1.0, 5
    u = rng.normal(0, gamma, n + 1)
    r = rng.standard_normal(n) + np.diff(u)
    plain = block_returns(r, k)
    gapped = block_returns(r, k, gap=True)
    rho_plain = np.corrcoef(plain[1:], plain[:-1])[0, 1]
    rho_gap = np.corrcoef(gapped[1:], gapped[:-1])[0, 1]
    assert rho_plain == pytest.approx(-gamma**2 / (k + 2 * gamma**2), abs=0.01)
    assert abs(rho_gap) < 0.01
    # calibration windows have the block length of the monitored statistic
    assert calibration_block_returns(r[:100], k, gap=True).shape[-1] == 100 - (k - 1) + 1
    assert block_returns(np.arange(10.0), 5, gap=True).tolist() == [1 + 2 + 3 + 4, 6 + 7 + 8 + 9]


# ----------------------------------------------------------------------------
# Calibration: missing minutes, vectorised against streaming
# ----------------------------------------------------------------------------


def test_calibration_windows_with_a_missing_minute():
    """Only the windows that contain the missing return are missing (review round 4, M6)."""
    r = np.arange(1.0, 21.0)
    r[3] = np.nan
    R = calibration_block_returns(r, 5)
    assert R.shape == (16,)
    assert np.all(np.isnan(R[:4])) and np.all(np.isfinite(R[4:]))
    assert R[4] == pytest.approx(5 + 6 + 7 + 8 + 9)
    Rg = calibration_block_returns(r, 5, gap=True)  # windows of 4 returns
    assert np.all(np.isnan(Rg[:4])) and np.all(np.isfinite(Rg[4:]))


@pytest.mark.parametrize("pos", [0, 10, 700])
def test_monitor_calibrates_on_time_with_a_missing_minute(pos):
    rng = np.random.default_rng(12)
    cal = 1440
    r = rng.standard_normal(cal + 200)
    ref = MultiFrequencyMonitor(threshold=1e9, freqs=FREQS, cal_len=cal)
    ref.run(r)
    r2 = r.copy()
    r2[pos] = np.nan
    mon = MultiFrequencyMonitor(threshold=1e9, freqs=FREQS, cal_len=cal)
    mon.run(r2)
    seg, seg0 = mon.segments[0], ref.segments[0]
    assert seg["start"] == cal and not seg["extended"] and not seg["fallback"]
    assert mon.cal_stats["restarted"] == 0
    for j in range(len(FREQS)):
        assert seg["scale"][j] == pytest.approx(seg0["scale"][j], rel=0.05)
        assert seg["m0"][j] == pytest.approx(seg0["m0"][j], rel=0.05)


@pytest.mark.parametrize("gap", [False, True])
def test_plugin_calibration_matches_streaming(gap):
    rng = np.random.default_rng(13)
    cal = 300
    r = rng.standard_t(4, (3, cal + 10))
    r[2, 50] = np.nan  # one row with a missing minute
    for i in range(3):
        mon = MultiFrequencyMonitor(threshold=1e9, freqs=FREQS, cal_len=cal, gap=gap, null_band=0.25)
        mon.run(r[i])
        seg = mon.segments[0]
        for j, k in enumerate(FREQS):
            scale, tau2, m0, m_up, m_dn = plugin_calibration(r[i:i + 1, :cal], k, 3.0, gap, null_band=0.25)
            assert scale[0] == pytest.approx(seg["scale"][j], rel=1e-12)
            assert tau2[0] == pytest.approx(seg["tau2"][j], rel=1e-12)
            assert m0[0] == pytest.approx(seg["m0"][j], rel=1e-12)
            assert m_up[0] == pytest.approx(1.25 * seg["m0"][j], rel=1e-12)
            assert m_dn[0] == pytest.approx(0.75 * seg["m0"][j], rel=1e-12)


def test_separate_level_for_the_decrease_side():
    rng = np.random.default_rng(14)
    r = rng.standard_normal(20_000) * np.repeat([1.0, 0.5, 1.6, 0.6, 1.0], 4000)
    base = MultiFrequencyMonitor(threshold=50.0, freqs=FREQS, cal_len=200)
    same = MultiFrequencyMonitor(threshold=50.0, freqs=FREQS, cal_len=200, threshold_down=50.0)
    e1, e2 = base.run(r), same.run(r)
    assert [(e.index, e.direction) for e in e1] == [(e.index, e.direction) for e in e2]
    assert any(e.direction == "down" for e in e1)
    n_down = sum(e.direction == "down" for e in e1)
    high = MultiFrequencyMonitor(threshold=50.0, freqs=FREQS, cal_len=200, threshold_down=1e6).run(r)
    assert sum(e.direction == "down" for e in high) < n_down
    assert any(e.direction == "up" for e in high)
    assert same.proved_bound(1440) == pytest.approx(base.proved_bound(1440))
    # 1 / (1 / (2 A) + 1 / (2 A_down)) replaces A in the bound
    two = MultiFrequencyMonitor(threshold=50.0, freqs=FREQS, threshold_down=150.0)
    assert two.proved_bound(1440) == pytest.approx(proved_bound(1440, 1 / (0.5 / 50 + 0.5 / 150), FREQS, "balanced"))
