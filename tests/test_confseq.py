"""Tests for the confidence sequences in ``avmonitor.confseq``."""

import numpy as np
import pytest

from avmonitor.confseq import (
    StreamingEBConfSeq,
    cs_baseline,
    eb_confseq,
    psi_coef,
    width_stopping_time,
)


def test_psi_coef_is_stable_and_matches_the_formula():
    x = np.array([0.0, 1e-8, 1e-4, 9.9e-4, 1.1e-3, 0.1, 0.5, 0.9, 0.99])
    h = psi_coef(x)
    assert h[0] == 0.5 and np.all(np.diff(h) > 0)  # increasing, h(0) = 1/2
    big = x >= 1e-3
    np.testing.assert_allclose(h[big], (-np.log1p(-x[big]) - x[big]) / x[big] ** 2, rtol=1e-12)
    # continuity across the switch between series and direct formula
    np.testing.assert_allclose(psi_coef(1e-3 - 1e-12), psi_coef(1e-3 + 1e-12), rtol=1e-9)


def test_fan_inequality_on_a_grid():
    """log(1 + lam xi) >= lam xi - psi_a(lam) xi^2 for xi >= -a and lam a < 1 (docs/math.md, Lemma 6.1)."""
    for a in (0.05, 0.3, 1.0, 4.0):
        xi = np.linspace(-a, 10 * a, 2001)
        for la in (0.01, 0.3, 0.75, 0.99):
            lam = la / a
            psi = lam**2 * psi_coef(la)
            assert np.all(np.log1p(lam * xi) >= lam * xi - psi * xi**2 - 1e-12)


LAWS = {
    "beta(0.5, 3)": lambda rng, s: rng.beta(0.5, 3.0, s),
    "bernoulli(0.1)": lambda rng, s: (rng.random(s) < 0.1).astype(float),
    "truncated normal square": lambda rng, s: np.minimum(rng.standard_normal(s) ** 2, 9.0) / 9.0,
    "uniform": lambda rng, s: rng.random(s),
}
MEANS = {"beta(0.5, 3)": 0.5 / 3.5, "bernoulli(0.1)": 0.1, "uniform": 0.5}


@pytest.mark.parametrize("law", list(LAWS))
@pytest.mark.parametrize("n_plan", [None, 300])
def test_time_uniform_coverage_on_iid_data(law, n_plan):
    """P(the mean leaves the sequence at ANY t <= 1500) <= delta, estimated from 2500 independent paths."""
    rng = np.random.default_rng(100)
    reps, n, delta = 2500, 1500, 0.05
    x = LAWS[law](rng, (reps, n))
    m = MEANS.get(law)
    if m is None:
        m = float(np.mean(LAWS[law](np.random.default_rng(1), 4_000_000)))
    cs = eb_confseq(x, 1.0, delta, delta, n_plan=n_plan)
    miss_up = np.mean((cs.upper <= m).any(axis=1))
    miss_lo = np.mean((cs.lower >= m).any(axis=1))
    se = np.sqrt(delta * (1 - delta) / reps)
    assert miss_up <= delta + 3 * se
    assert miss_lo <= delta + 3 * se


def test_coverage_of_weighted_average_with_varying_means():
    """With conditional means that change over time the sequence covers the bet-weighted average (Theorem 6.2)."""
    rng = np.random.default_rng(101)
    reps, n, delta = 3000, 1500, 0.05
    t = np.arange(n)
    p = 0.1 + 0.08 * np.sin(2 * np.pi * t / 400)  # deterministic, so conditional means are known
    x = (rng.random((reps, n)) < p).astype(float)
    cs = eb_confseq(x, 1.0, delta, delta, n_plan=500)
    nu_up = np.cumsum(cs.lam_upper * p, axis=1) / np.cumsum(cs.lam_upper, axis=1)
    nu_lo = np.cumsum(cs.lam_lower * p, axis=1) / np.cumsum(cs.lam_lower, axis=1)
    se = np.sqrt(delta * (1 - delta) / reps)
    assert np.mean((cs.upper <= nu_up).any(axis=1)) <= delta + 3 * se
    assert np.mean((cs.lower >= nu_lo).any(axis=1)) <= delta + 3 * se


def test_one_sided_bound_with_means_above_a_constant():
    """If every conditional mean is at least m*, the upper sequence stays above m* (Corollary 6.3)."""
    rng = np.random.default_rng(102)
    reps, n, delta = 2000, 1000, 0.05
    p = np.where(np.arange(n) % 200 < 100, 0.2, 0.35)
    x = (rng.random((reps, n)) < p).astype(float)
    cs = eb_confseq(x, 1.0, delta, delta, n_plan=n)
    assert np.mean((cs.upper <= 0.2).any(axis=1)) <= delta
    assert np.mean((cs.lower >= 0.35).any(axis=1)) <= delta


def test_width_shrinks_with_sample_size():
    rng = np.random.default_rng(103)
    x = np.minimum(rng.standard_normal((200, 20000)) ** 2, 9.0)
    cs = eb_confseq(x, 9.0, 0.01, 0.01)
    w = np.median(cs.upper - cs.lower, axis=0)
    idx = [99, 999, 9999, 19999]
    assert np.all(np.diff(w[idx]) < 0)
    # roughly the 1 / sqrt(n) rate (up to the log factor of the time-uniform schedule)
    assert 2.0 < w[999] / w[9999] < 4.5
    assert np.all(cs.lower <= cs.upper)
    assert np.all((cs.lower >= 0) & (cs.upper <= 9.0))


def test_streaming_equals_vectorised():
    rng = np.random.default_rng(104)
    x = np.minimum(rng.standard_normal(800) ** 2, 6.25)
    for plan in (None, 500):
        cs = eb_confseq(x, 6.25, 0.02, 0.005, n_plan=plan)
        s = StreamingEBConfSeq(6.25, 0.02, 0.005, n_plan=plan)
        for k in range(x.size):
            s.update(x[k])
            if k in (0, 10, 399, 799):
                np.testing.assert_allclose([s.upper, s.lower, s.mean],
                                           [cs.upper[0, k], cs.lower[0, k], cs.mean[0, k]], rtol=1e-10, atol=1e-14)


def test_input_checks():
    with pytest.raises(ValueError):
        eb_confseq([0.5, 1.5], 1.0)
    with pytest.raises(ValueError):
        eb_confseq([0.5, np.nan], 1.0)
    with pytest.raises(ValueError):
        eb_confseq([0.5], 1.0, delta_upper=0.0)
    with pytest.raises(ValueError):
        eb_confseq([0.5], 1.0, c=1.0)
    with pytest.raises(ValueError):
        StreamingEBConfSeq(1.0).update(2.0)


def test_cs_baseline_split_and_stopping():
    rng = np.random.default_rng(105)
    r = rng.standard_normal((500, 1440))
    b = cs_baseline(r, c=3.0, delta_upper=0.01, delta_lower=0.01, tau_frac=0.25)
    assert b.n_tau == 360 and np.all(b.n_cs == 1080)
    np.testing.assert_allclose(b.tau, 3.0 * np.median(np.abs(r[:, :360]), axis=1) * 1.482602218505602, rtol=1e-6)
    assert np.all(b.lower < b.mean) and np.all(b.mean < b.upper)
    # true mean of min(Z^2, tau^2) given tau, by numerical integration per row
    from scipy import stats

    a = b.tau
    m_true = (2 * stats.norm.cdf(a) - 1) - 2 * a * stats.norm.pdf(a) + 2 * a**2 * stats.norm.sf(a)
    assert np.mean(b.upper <= m_true) <= 0.02 and np.mean(b.lower >= m_true) <= 0.02
    x = np.minimum(r[:, 360:] ** 2, (b.tau**2)[:, None])
    cs = eb_confseq(x, b.tau**2, 0.01, 0.01, n_plan=1080)
    k = width_stopping_time(cs, 0.6, n_min=100)
    assert np.all((k >= 100) & (k <= 1080))
    w = cs.rel_width()[np.arange(500), k - 1]
    assert np.all((w <= 0.6) | (k == 1080))
    b2 = cs_baseline(r, delta_upper=0.01, delta_lower=0.01, stop=k)
    np.testing.assert_allclose(b2.upper, cs.upper[np.arange(500), k - 1])


def test_calibrate_default_unchanged_and_cs_mode():
    from avmonitor.detectors import calibrate

    rng = np.random.default_rng(106)
    r = rng.standard_t(5, 1440)
    b = calibrate(r, 3.0, margin_z=1.0, null_band=0.1)
    scale = float(np.median(np.abs(r)) * 1.482602218505602)
    x = np.minimum(r**2, (3 * scale) ** 2)
    se = np.sqrt(x.var(ddof=1) / x.size)
    np.testing.assert_allclose([b.scale, b.m0, b.m0_up, b.m0_down], [scale, x.mean(), (x.mean() + se) * 1.1,
                                                                      (x.mean() - se) * 0.9], rtol=1e-6)
    c = calibrate(r, 3.0, baseline="cs", delta=0.02)
    v = cs_baseline(r[None, :], 3.0, 0.01, 0.01)
    np.testing.assert_allclose([c.m0_up, c.m0_down, c.m0, c.tau], [v.upper[0], v.lower[0], v.mean[0], v.tau[0]])
    assert c.m0_down < c.m0 < c.m0_up and c.n_cal == 1440
    with pytest.raises(ValueError):
        calibrate(r, baseline="cs", margin_z=1.0)
