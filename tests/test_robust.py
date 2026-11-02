import numpy as np
import pytest
from scipy import stats

from avmonitor import robust as R

N, REPS, SIG = 390, 4000, 0.3
IQ = SIG**4


@pytest.fixture(scope="module")
def gauss():
    rng = np.random.default_rng(123)
    return rng.standard_normal((REPS, N)) * SIG / np.sqrt(N)


@pytest.mark.parametrize(
    "est",
    [R.realized_quarticity, R.tripower_quarticity, R.quadpower_quarticity, R.bipower_quarticity, R.min_rq, R.med_rq],
)
def test_classical_estimators_unbiased_constant_vol(gauss, est):
    ratio = est(gauss).mean() / IQ
    assert abs(ratio - 1) < 0.03


def test_normalised_kernel_iq_unbiased(gauss):
    assert abs(R.kernel_iq(gauss, 30).mean() / IQ - 1) < 0.02


def test_printed_kernel_scale_is_weight_sum_over_n_squared(gauss):
    # Eq. (7) as printed estimates W * sigma_step^4 with W the raw kernel weight sum
    for h in (10, 30, 50):
        W = R.kernel_weight_sum(h)
        ratio = R.kernel_quarticity_printed(gauss, h).mean() / IQ
        assert abs(ratio / (W / N**2) - 1) < 0.03


def test_truncated_moment_factor_matches_monte_carlo():
    z = np.random.default_rng(0).standard_normal(2_000_000)
    for a in (1.5, 2.0234692, 3.0):
        mc4 = np.mean(np.clip(z, -a, a) ** 4) / 3
        mc2 = np.mean(np.clip(z, -a, a) ** 2)
        assert abs(R.truncated_moment_factor(a, 4) - mc4) < 0.01
        assert abs(R.truncated_moment_factor(a, 2) - mc2) < 0.005


def test_effective_cutoff_raw_mad():
    assert abs(R.effective_cutoff(3, scaled=False) - 3 * stats.norm.ppf(0.75)) < 1e-12
    assert abs(R.effective_cutoff(3, scaled=True) - 3) < 1e-12


def test_raw_mad_truncation_bias_does_not_vanish_with_n():
    rng = np.random.default_rng(5)
    kappa = R.truncated_moment_factor(R.effective_cutoff(3, False), 4)
    for n in (390, 23400):
        r = rng.standard_normal((200, n)) / np.sqrt(n)
        ratio = R.truncated_rq_mad(r, 3.0, mad_scaled=False).mean()
        assert abs(ratio / kappa - 1) < 0.03
        corrected = R.truncated_rq_mad(r, 3.0, mad_scaled=False, bias_correct=True).mean()
        assert abs(corrected - 1) < 0.03


def test_aggregated_kernel_is_weighted_realized_quarticity(gauss):
    # Eq. (7) equals sum_i c_i r_i^4 / (3 (n - h)) where c_i adds the weights that
    # return i receives over all windows. Interior returns all get c_i = W.
    r = gauss[:3]
    h = 20
    w = R.kernel_weights(h, normalise=False)
    n = r.shape[1]
    c = np.zeros(n)
    for t in range(h, n):
        c[t - h : t + 1] += w
    direct = (r**4 @ c) / (3 * (n - h))
    assert np.allclose(direct, R.kernel_quarticity_printed(r, h))
    assert np.allclose(c[h : n - h], w.sum())


def test_reverse_weights_change_local_but_not_much_aggregate(gauss):
    r = gauss[:50]
    a = R.local_quarticity_printed(r, 30)
    b = R.local_quarticity_printed(r, 30, reverse_weights=True)
    assert not np.allclose(a, b)
    agg_a, agg_b = a.mean(axis=1), b.mean(axis=1)
    assert np.all(np.abs(agg_a / agg_b - 1) < 0.1)


def test_robust_estimator_has_bounded_influence(gauss):
    r = gauss[0].copy()
    base = R.kernel_iq(r, 30, c=3.0)
    tau = R.mad_threshold(r, 3.0, scaled=True)
    r[200] = 1e6
    bumped = R.kernel_iq(r, 30, c=3.0)
    # One return can add at most max weight * tau^4 / (3 delta^2) * T / (#windows) per window it enters
    bound = (tau**4) / (3 * (1 / N) ** 2) * R.kernel_weights(30).max() * 31 / (N - 30) * 1.1
    assert bumped - base < bound
    assert R.realized_quarticity(r) / base > 1e10


def test_threshold_rq_robust_to_single_jump(gauss):
    r = gauss[:200].copy()
    r[:, 100] += 0.2  # about 13 diffusive standard deviations
    assert abs(R.threshold_rq(r).mean() / IQ - 1) < 0.05
    assert R.realized_quarticity(r).mean() / IQ > 10


def test_kernel_weights_shapes_and_newest_largest():
    for k in R.KERNELS:
        w = R.kernel_weights(20, k)
        assert w.shape == (21,)
        assert abs(w.sum() - 1) < 1e-12
        if k != "uniform":
            assert w[-1] >= w[0]


def test_online_scale_estimates_track_sd():
    rng = np.random.default_rng(9)
    r = np.concatenate([rng.standard_normal(5000), 2 * rng.standard_normal(5000)])
    roll = R.rolling_mad_zero(r, 500)
    ew = R.ew_quantile_scale(r, halflife=200)
    assert abs(np.nanmean(roll[4000:5000]) - 1) < 0.1
    assert abs(np.nanmean(roll[9000:]) - 2) < 0.2
    assert abs(ew[4000:5000].mean() - 1) < 0.15
    assert abs(ew[9000:].mean() - 2) < 0.3
    # predictable: entry t must not depend on r[t]
    r2 = r.copy()
    r2[6000] = 1e9
    assert roll[6000] == R.rolling_mad_zero(r2, 500)[6000]
    assert ew[6000] == R.ew_quantile_scale(r2, halflife=200)[6000]
