"""Fast checks of the competitor implementations (no data needed)."""

import numpy as np
import pytest
from scipy import stats
from scipy.special import gammaln, logsumexp

from avmonitor import competitors as C


def test_gumbel_critical_value_matches_simulated_max_of_independent_normals():
    rng = np.random.default_rng(1)
    N = 2000
    m = np.abs(rng.standard_normal((4000, N))).max(axis=1)
    crit = C.gumbel_crit_max_abs_normal(N, 0.05)
    # the extreme value approximation is slightly conservative at this N; it must be close to 5%
    assert 0.02 <= np.mean(m > crit) <= 0.07


def test_mosum_critical_value_formula():
    # symmetric case of mosum::mosum.criticalValue: (b - log(log(1 / sqrt(1 - alpha)))) / a
    n, k, alpha = 1000, 50, 0.1
    x = n / k
    a = np.sqrt(2 * np.log(x))
    b = 2 * np.log(x) + 0.5 * np.log(np.log(x)) + np.log(1.5) - 0.5 * np.log(np.pi)
    assert C.mosum_crit(n, k, alpha) == pytest.approx((b - np.log(np.log(1 / np.sqrt(1 - alpha)))) / a)
    with pytest.raises(ValueError):
        C.mosum_crit(100, 50)


def test_truncated_square_dispersion():
    assert C.truncated_square_dispersion(None) == 2.0
    rng = np.random.default_rng(2)
    z = rng.standard_normal(2_000_000)
    x = np.where(np.abs(z) <= 2.0, z**2, 0.0)
    assert C.truncated_square_dispersion(2.0) == pytest.approx(x.var() / x.mean() ** 2, rel=0.01)


def test_block_statistic_is_standard_normal_and_brute_force():
    rng = np.random.default_rng(3)
    r = rng.standard_normal((3000, 400))
    z, pos = C.bjv_block_statistics(r, 100, trunc_c=None)
    assert z.shape == (3000, 3) and list(pos) == [100, 200, 300]
    assert abs(z.mean()) < 0.05 and abs(z.std() - 1) < 0.05
    # neighbouring comparisons share a block: correlation about -1/2
    assert np.corrcoef(z[:, 0], z[:, 1])[0, 1] == pytest.approx(-0.5, abs=0.05)
    x = r[7]
    L, R = x[100:200] ** 2, x[200:300] ** 2
    expected = (R.mean() - L.mean()) / ((L.mean() + R.mean()) / 2 * np.sqrt(2 * (2 / 100)))
    assert z[7, 1] == pytest.approx(expected)


def test_overlap_statistic_brute_force_with_truncation_and_nan():
    rng = np.random.default_rng(4)
    r = rng.standard_normal((5, 600))
    r[2, 250:255] = np.nan
    r[3, 310] = 40.0  # a jump
    k, c = 50, 4.0
    D = C.bjv_overlap_statistics(r, k, c)
    D_small_chunks = C.bjv_overlap_statistics(r, k, c, chunk=2, max_cells=97)
    assert np.allclose(np.nan_to_num(D, nan=9.0), np.nan_to_num(D_small_chunks, nan=9.0))
    rho = C.truncated_square_dispersion(c)
    for i, b in ((0, 300), (2, 260), (3, 300)):
        Lr, Rr = r[i, b - k : b], r[i, b : b + k]

        def bpv(v):
            p = np.abs(v[1:] * v[:-1])
            ok = np.isfinite(p)
            return np.pi / 2 * p[ok].mean()

        u2 = c**2 * max(bpv(Lr), bpv(Rr))
        kl = np.isfinite(Lr) & (Lr**2 <= u2)
        kr = np.isfinite(Rr) & (Rr**2 <= u2)
        nl, nr = np.isfinite(Lr).sum(), np.isfinite(Rr).sum()
        sl, sr = np.sum(Lr[kl] ** 2), np.sum(Rr[kr] ** 2)
        z = (sr / nr - sl / nl) / ((sl + sr) / (nl + nr) * np.sqrt(rho * (1 / nl + 1 / nr)))
        assert D[i, b] == pytest.approx(z)
    assert np.all(np.isnan(D[:, :k])) and np.all(np.isnan(D[:, 600 - k + 1 :]))
    # the jump is removed by truncation, so the plain statistic reacts much more
    Dp = C.bjv_overlap_statistics(r, k, None)
    assert abs(Dp[3, 300]) > 3 * abs(D[3, 300])


def test_online_statistic_uses_only_the_past():
    rng = np.random.default_rng(5)
    r = rng.standard_normal((2, 500))
    S1 = C.bjv_online_statistic(r, 40, 4.0)
    r2 = r.copy()
    r2[:, 300:] = rng.standard_normal((2, 200)) * 5
    S2 = C.bjv_online_statistic(r2, 40, 4.0)
    assert np.allclose(np.nan_to_num(S1[:, :300]), np.nan_to_num(S2[:, :300]))
    D = C.bjv_overlap_statistics(r, 40, 4.0)
    assert S1[0, 199] == pytest.approx(D[0, 160])


def test_bjv_test_size_power_and_change_point():
    rng = np.random.default_rng(6)
    null = rng.standard_normal((400, 1000))
    for variant in ("blocks", "overlap"):
        res = C.bjv_test(null, 50, 0.05, variant, 4.0)
        assert res.reject.mean() <= 0.08  # the Gumbel critical values are conservative here
    r = rng.standard_normal((200, 1000))
    r[:, 600:] *= 2.0
    # with pooled standardisation |D| <= sqrt(k), so a doubling of the sd needs a window of about 100
    res = C.bjv_test(r, 100, 0.05, "overlap", 4.0)
    assert res.reject.mean() > 0.95
    assert np.median(np.abs(res.cp_index - 600)) <= 20
    assert np.all(res.direction[res.reject] > 0)


def test_bocpd_predictive_is_student_t():
    f = C.BOCPDVariance(1, 0.01, alpha0=3.0, beta0=2.0, r_max=16, short=4)
    # a single bucket at run length 0: the predictive is t with 2 alpha df and scale sqrt(beta / alpha)
    x = 1.3
    a, b = 3.0, 2.0
    lp = (gammaln(a + 0.5) - gammaln(a) - 0.5 * np.log(2 * np.pi * b) - (a + 0.5) * np.log1p(x * x / (2 * b)))
    assert lp == pytest.approx(stats.t.logpdf(x, df=2 * a, scale=np.sqrt(b / a)))
    f.update(x)
    # after one step: P(r=0) = H, P(r=1) = 1 - H, and the run of length 1 has absorbed x
    p = np.exp(f.logR[0])
    assert p[0] == pytest.approx(0.01) and p[1] == pytest.approx(0.99)
    assert f.beta[0, 1] == pytest.approx(b + x * x / 2)


def test_bocpd_matches_exact_recursion_without_cap():
    rng = np.random.default_rng(7)
    x = rng.standard_normal(30) * np.r_[np.ones(15), 3 * np.ones(15)]
    H, a0, b0 = 0.05, 2.0, 1.0
    # exact Adams and MacKay recursion with an unbounded run-length vector
    logR = np.array([0.0])
    alphas, betas = np.array([a0]), np.array([b0])
    f = C.BOCPDVariance(1, H, a0, b0, r_max=64, short=5)
    for v in x:
        lp = stats.t.logpdf(v, df=2 * alphas, scale=np.sqrt(betas / alphas))
        lj = logR + lp
        newR = np.r_[logsumexp(lj) + np.log(H), lj + np.log1p(-H)]
        logR = newR - logsumexp(newR)
        alphas = np.r_[a0, alphas + 0.5]
        betas = np.r_[b0, betas + v * v / 2]
        f.update(v)
    n = logR.size
    assert np.allclose(np.exp(f.logR[0, :n]), np.exp(logR), atol=1e-12)
    assert f.short_mass()[0] == pytest.approx(np.exp(logsumexp(logR[:5])))
    assert f.short_log_odds()[0] == pytest.approx(logsumexp(logR[:5]) - logsumexp(logR[5:]))


def test_bocpd_cap_keeps_a_proper_distribution_and_detects():
    rng = np.random.default_rng(8)
    cal = rng.standard_normal((50, 1000))
    r = rng.standard_normal((50, 600))
    r[:, 300:] *= 2.5
    r[3, 10] = np.nan
    P, S = C.bocpd_paths(r, 1 / 1000, cal=cal, r_max=64, short=20)
    assert np.all((P >= 0) & (P <= 1))
    assert np.quantile(P[:, :300].max(axis=1), 0.5) < 0.5
    a = C.first_crossing(P, 0.5)
    assert np.mean(a >= 300) > 0.9
    assert np.mean(S[np.arange(50), np.maximum(a, 0)][a >= 300] > 0) > 0.9
    with pytest.raises(ValueError):
        C.BOCPDVariance(1, 0.01, alpha0=1.0)


def test_scan_alarms_dead_time_and_matching():
    up = np.zeros(100)
    dn = np.zeros(100)
    up[[10, 12, 40]] = 5.0
    dn[[20, 70]] = 5.0
    ev = C.scan_alarms(up, dn, 1.0, 1.0, start=0, dead=15)
    assert ev == [(10, "up"), (40, "up"), (70, "down")]
    assert C.first_crossing(up[None, :], 1.0, start=11)[0] == 12
    rng = np.random.default_rng(9)
    s = rng.standard_normal(5000)
    hu, hd = C.match_alarm_counts(s, -s, 7, 4, start=0, dead=50, lo_hi=(0.0, 10.0))
    ev = C.scan_alarms(s, -s, hu, hd, 0, 50)
    assert sum(e[1] == "up" for e in ev) == 7 and sum(e[1] == "down" for e in ev) == 4
