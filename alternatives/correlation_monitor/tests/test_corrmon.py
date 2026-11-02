"""Tests for the correlation monitor prototype (synthetic only, no data files, run from the repo root)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import corrmon as C  # noqa: E402
from simulate_pairs import pair_returns  # noqa: E402


def test_sign_mean_matches_arcsin_law_for_gaussian_and_student_t():
    rng = np.random.default_rng(0)
    for df in (None, 4.0):
        p = pair_returns(1, 400_000, rho=0.45, df=df, rng=rng)
        m = np.nanmean(C.sign_concordance(p.r1[0], p.r2[0]))
        assert abs(m - C.gaussian_sign_mean(0.45)) < 0.01
        assert abs(C.rho_from_sign_mean(m) - 0.45) < 0.02


def test_sign_statistic_is_scale_invariant_and_product_is_not():
    rng = np.random.default_rng(1)
    p = pair_returns(1, 20_000, rho=0.3, rng=rng)
    r1, r2 = p.r1[0], p.r2[0]
    base = C.calibrate_pair(r1, r2, "sign")
    s_a = C.pair_statistic(r1, r2, base)
    s_b = C.pair_statistic(3.0 * r1, 0.2 * r2, base)
    assert np.array_equal(s_a, s_b)
    assert -1 <= s_a.min() and s_a.max() <= 1
    pb = C.calibrate_pair(r1, r2, "product", c=3.0)
    prod_a = C.pair_statistic(r1, r2, pb).mean()
    prod_b = C.pair_statistic(2.0 * r1, r2, pb).mean()
    assert prod_b > 1.5 * prod_a  # doubling one asset's volatility roughly doubles the product's mean
    assert abs(prod_a - 0.3) < 0.03 and pb.lo == -9.0 and pb.hi == 9.0


def test_cojumps_with_independent_signs_leave_the_sign_mean_unchanged():
    rng = np.random.default_rng(2)
    clean = pair_returns(1, 300_000, rho=0.3, rng=rng)
    noisy = pair_returns(1, 300_000, rho=0.3, p_cojump=0.02, p_jump=0.01, rng=rng)
    m_clean = np.mean(C.sign_concordance(clean.r1[0], clean.r2[0]))
    m_noisy = np.mean(C.sign_concordance(noisy.r1[0], noisy.r2[0]))
    assert abs(m_noisy - m_clean) < 0.01


def test_calibration_band_and_admissible_bets():
    rng = np.random.default_rng(3)
    p = pair_returns(1, 5000, rho=0.3, rng=rng)
    b = C.calibrate_pair(p.r1[0], p.r2[0], "sign", band=0.1)
    assert b.m_up == pytest.approx(b.m0 + 0.1) and b.m_dn == pytest.approx(b.m0 - 0.1)
    assert 0 < b.lam_max_up < 1 and 0 < b.lam_max_dn
    # the factors are strictly positive for every admissible bet and every value of the statistic
    for y in (-1.0, 0.0, 1.0):
        for frac in C.DEFAULT_FRACTIONS:
            assert 1 + frac * b.lam_max_up * (y - b.m_up) > 0
            assert 1 + frac * b.lam_max_dn * (b.m_dn - y) > 0


def test_false_alarm_bound_holds_for_both_statistics_and_both_sides():
    rng = np.random.default_rng(4)
    p = pair_returns(1500, 60, rho=0.3, p_jump=0.01, rng=rng)
    A = 300.0
    for stat in ("sign", "product"):
        if stat == "sign":
            m0, lo, hi = C.gaussian_sign_mean(0.3), -1.0, 1.0
            y = C.sign_concordance(p.r1, p.r2)
        else:
            m0 = C.gaussian_product_moments(0.3, 3.0)[0]
            lo, hi = -9.0, 9.0
            y = C.truncated_product(p.r1, p.r2, 1.0, 1.0, 3.0)
        for direction in ("up", "down"):
            res = C.bounded_edetector(y, m0, lo, hi, direction, threshold=A)
            assert (res.alarm >= 0).mean() <= 60 / A + 0.03  # bound 0.2 plus Monte Carlo slack


def test_detects_increase_and_decrease_and_ignores_volatility_change():
    rng = np.random.default_rng(5)
    up = pair_returns(200, 2500, rho=0.3, nu=500, rho1=0.8, rng=rng)
    y = C.sign_concordance(up.r1, up.r2)
    res = C.bounded_edetector(y, C.gaussian_sign_mean(0.3), -1, 1, "up", threshold=2 * 7800)
    det = res.alarm >= 500
    assert det.mean() > 0.95 and np.median(res.alarm[det] - 500) < 400
    dn = pair_returns(200, 2500, rho=0.8, nu=500, rho1=0.3, rng=rng)
    y = C.sign_concordance(dn.r1, dn.r2)
    res = C.bounded_edetector(y, C.gaussian_sign_mean(0.8), -1, 1, "down", threshold=2 * 7800)
    det = res.alarm >= 500
    assert det.mean() > 0.95 and np.median(res.alarm[det] - 500) < 400
    vol = pair_returns(200, 2500, rho=0.3, nu=500, vol1_post=3.0, rng=rng)
    y = C.sign_concordance(vol.r1, vol.r2)
    res_up = C.bounded_edetector(y, C.gaussian_sign_mean(0.3), -1, 1, "up", threshold=2 * 7800)
    res_dn = C.bounded_edetector(y, C.gaussian_sign_mean(0.3), -1, 1, "down", threshold=2 * 7800)
    assert (C.two_sided_alarm(res_up, res_dn) >= 0).mean() <= 2500 / 7800 + 0.05
    # the product statistic reacts to the volatility change (this is the known weakness)
    yp = C.truncated_product(vol.r1, vol.r2, 1.0, 1.0, 3.0)
    m0p = C.gaussian_product_moments(0.3, 3.0)[0]
    res_p = C.bounded_edetector(yp, m0p, -9, 9, "up", threshold=2 * 7800)
    assert (res_p.alarm >= 500).mean() > 0.9


def test_gaussian_lr_factor_has_unit_mean_under_the_null_and_is_an_edetector_there():
    rng = np.random.default_rng(6)
    p = pair_returns(1, 400_000, rho=0.3, rng=rng)
    f = C.gaussian_corr_lr_factors(p.r1[0], p.r2[0], 0.3, 0.6)
    assert abs(f.mean() - 1.0) < 0.02
    q = pair_returns(1500, 60, rho=0.3, rng=rng)
    res = C.gaussian_corr_lr_edetector(q.r1, q.r2, 0.3, "up", threshold=300.0)
    assert (res.alarm >= 0).mean() <= 0.2 + 0.03


def test_streaming_monitor_recalibrates_and_handles_nans():
    rng = np.random.default_rng(7)
    p = pair_returns(1, 7000, rho=0.3, nu=4000, rho1=0.85, rng=rng)
    r1, r2 = p.r1[0].copy(), p.r2[0].copy()
    r1[100:110] = np.nan
    # band 0.1: with a plug-in baseline from 1000 steps (standard error about 0.03 on the sign
    # scale) and band 0 the strict null is false by about one standard error, which can give
    # a pre-change alarm (as E4d found for the volatility monitor); the band absorbs it
    mon = C.CorrelationMonitor(threshold=5760.0, stat="sign", cal_len=1000, null_band=0.1)
    events = mon.run(r1, r2)
    assert len(events) == 1
    first = events[0]
    assert first.direction == "up" and 4000 <= first.index < 4600
    assert first.monitored_since == 1000 + 10  # ten NaN steps were skipped during calibration
    assert mon.segments[0]["cal_len"] == 1000
    assert abs(first.rho_hat - 0.3) < 0.15
    # after the alarm a fresh calibration starts and the monitor keeps running
    assert len(mon.segments) == 2 and mon.segments[1]["start"] == first.index + 1001
    # benchmark protocol runs end to end
    g = C.CorrelationMonitor(threshold=5760.0, cal_len=1000, benchmark="gaussian", null_band=0.1)
    ge = g.run(r1, r2)
    assert any(e.direction == "up" and e.index >= 4000 for e in ge)


def test_aggregate_returns_requires_full_bars():
    pd = pytest.importorskip("pandas")
    idx = pd.date_range("2024-01-01", periods=12, freq="1min", tz="UTC")
    s = pd.Series(np.arange(12, dtype=float), index=idx)
    s.iloc[7] = np.nan
    out = C.aggregate_returns(s, 5)
    assert out.iloc[0] == pytest.approx(0 + 1 + 2 + 3 + 4)
    assert np.isnan(out.iloc[1]) and np.isnan(out.iloc[2])  # a missing minute and an incomplete bar
