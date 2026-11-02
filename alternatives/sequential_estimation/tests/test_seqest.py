"""Fast tests for the sequential estimation prototype (no data, no network, under a minute).

Each closed form used in the README is checked exactly or against Monte Carlo, and each
stopping rule is checked against a brute-force loop on the same numbers.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import seqest as S  # noqa: E402
from avmonitor.confseq import eb_confseq  # noqa: E402
from avmonitor.robust import truncated_moment_factor  # noqa: E402

# ----------------------------------------------------------------------------- closed forms


def test_ratio_limit_and_coverage_limit_values():
    # 2 log(2/alpha) / z^2 at the two levels used in the experiments
    assert S.cs_ratio_limit(0.05) == pytest.approx(2 * math.log(40) / 1.959963984540054**2, rel=1e-12)
    assert S.cs_ratio_limit(0.05) == pytest.approx(1.92056, abs=1e-5)
    assert S.cs_ratio_limit(0.10) == pytest.approx(2.21451, abs=1e-5)
    # coverage of the stopped sequence under a normal approximation: 1 - 2 (1 - Phi(sqrt(2 log 40)))
    assert S.cs_limit_coverage(0.05) == pytest.approx(0.99340, abs=1e-5)


def test_kappa_matches_package_and_monte_carlo():
    for a in (1.0, 2.0234, 3.0, 4.5):
        assert S.kappa2(a) == pytest.approx(truncated_moment_factor(a, 2), rel=1e-12)
        assert S.kappa4(a) == pytest.approx(3 * truncated_moment_factor(a, 4), rel=1e-12)
    z = np.random.default_rng(0).standard_normal(2_000_000)
    assert np.mean(np.minimum(z**2, 9.0)) == pytest.approx(float(S.kappa2(3.0)), rel=3e-3)
    assert np.mean(np.minimum(z**4, 81.0)) == pytest.approx(float(S.kappa4(3.0)), rel=5e-3)


def test_generic_truncated_moments_quadrature():
    g = S.Gaussian(1.3)
    tau = np.array([0.5, 2.0, 3.9])
    m2, m4 = S.Law.trunc_moments(g, tau)  # generic quadrature on the tail function
    c2, c4 = g.trunc_moments(tau)  # closed form
    np.testing.assert_allclose(m2, c2, rtol=1e-9)
    np.testing.assert_allclose(m4, c4, rtol=1e-9)
    rng = np.random.default_rng(1)
    for law in (S.StudentT(3), S.StudentT(5), S.Jumps(0.01, 5.0)):
        r, _, _ = law.draw(rng, None, np.arange(1), 2_000_000)
        r = r.ravel()
        tau = 3.0 * S.pilot_scale(r[None, :])[0]
        m2, m4 = law.trunc_moments(np.array([tau]))
        assert np.mean(np.minimum(r * r, tau * tau)) == pytest.approx(m2[0], rel=5e-3)
        assert np.mean(np.minimum(r**4, tau**4)) == pytest.approx(m4[0], rel=1e-2)
        assert np.var(r) == pytest.approx(law.theta_raw(), rel=0.05)


def test_jump_law_raw_moments_exact():
    j = S.Jumps(0.01, 5.0)
    assert j.theta_raw() == pytest.approx(1.25)
    assert j.m4_raw() == pytest.approx(0.99 * 3 + 0.01 * 3 * 26**2)
    assert S.StudentT(5).m4_raw() == pytest.approx(9.0)
    assert math.isinf(S.StudentT(3).m4_raw())


def test_bootstrap_population_moments_are_empirical():
    u = np.random.default_rng(2).standard_t(4, 1440)
    b = S.Bootstrap(u, block=30)
    for tau in (0.5, 1.7, 10.0):
        m2, m4 = b.trunc_moments(np.array([tau]))
        assert m2[0] == pytest.approx(np.mean(np.minimum(u**2, tau**2)), rel=1e-12)
        assert m4[0] == pytest.approx(np.mean(np.minimum(u**4, tau**4)), rel=1e-12)
    # the stationary bootstrap keeps the empirical law as its marginal
    rng = np.random.default_rng(3)
    st = b.start(rng, 400)
    r, _, _ = b.draw(rng, st, np.arange(400), 3000)
    assert np.mean(r**2) == pytest.approx(np.mean(u**2), rel=0.03)


def test_garch_autocorrelation_and_long_run_variance():
    g = S.Garch(0.10, 0.85)
    rng = np.random.default_rng(4)
    st = g.start(rng, 40)
    r, h, _ = g.draw(rng, st, np.arange(40), 50_000)
    x = r * r
    assert np.mean(x) == pytest.approx(1.0, rel=0.03)
    assert np.mean(x * x) == pytest.approx(g.kurtosis(), rel=0.08)
    xc = x - x.mean(axis=1, keepdims=True)
    rho = [np.mean(xc[:, k:] * xc[:, :-k]) / np.mean(xc * xc) for k in (1, 5)]
    assert rho[0] == pytest.approx(g.rho1(), abs=0.02)
    assert rho[1] == pytest.approx(g.rho1() * 0.95**4, abs=0.02)
    # batch means of length 2000 estimate the long-run variance
    bm = x[:, : 50_000 // 2000 * 2000].reshape(40, -1, 2000).mean(axis=2)
    lrv = 2000 * np.var(bm)
    assert lrv / (g.kurtosis() - 1) == pytest.approx(g.lrv_factor_raw(), rel=0.2)


def test_elasticity_and_sigma_mapping_roundtrip():
    # d log theta_tau / d log sigma^2 by finite differences
    a = 3.0
    eps = 1e-5
    f = lambda s2: s2 * S.kappa2(a / math.sqrt(s2))  # noqa: E731  tau = a at sigma = 1
    num = (math.log(f(1 + eps)) - math.log(f(1 - eps))) / (math.log(1 + eps) - math.log(1 - eps))
    assert S.elasticity_gaussian(a) == pytest.approx(num, rel=1e-6)
    sig2 = np.array([0.5, 1.0, 2.3])
    tau = np.array([3.0, 2.5, 4.0])
    th = sig2 * S.kappa2(tau / np.sqrt(sig2))
    np.testing.assert_allclose(S.sigma2_from_theta_tau(th, tau), sig2, rtol=1e-9)
    t5 = S.StudentT(5)
    th5 = t5.trunc_moments(np.array([2.0]))[0] * 1.7
    tau5 = 2.0 * math.sqrt(1.7)
    assert S.sigma2_from_theta_tau(th5, np.array([tau5]), shape=t5)[0] == pytest.approx(1.7, rel=1e-8)


def test_cs_prediction_tends_to_limit():
    mz = float(S.kappa2(3.0)) / 9
    vz = float(S.kappa4(3.0) - S.kappa2(3.0) ** 2) / 81
    for alpha in (0.05, 0.10):
        ns = S.zq(alpha) ** 2 * vz / (1e-4) ** 2
        r = S.cs_n_predict(mz, vz, 1e-4, alpha, bets="tuned") / ns
        assert r == pytest.approx(S.cs_ratio_limit(alpha), rel=5e-3)
        r2 = S.cs_n_predict(mz, vz, 1e-4, alpha, bets="opt") / ns
        assert r2 == pytest.approx(S.cs_ratio_limit(alpha), rel=1e-2)
    # at a practical d the optimised bets never need more data than the tuned ones (up to the grid)
    for d in (0.4, 0.2, 0.1):
        assert S.cs_n_predict(mz, vz, d / 9, 0.05, bets="opt") <= 1.01 * S.cs_n_predict(mz, vz, d / 9, 0.05)


def test_brownian_limit():
    # without intersection the limit is deterministic: 2 log(2/alpha) / (kappa (2 - kappa) z^2)
    assert S.brownian_limit(0.05, intersect=False)["ratio"] == pytest.approx(S.cs_ratio_limit(0.05), rel=1e-12)
    res = S.brownian_limit(0.05, reps=2000, chunk=1000, du=2e-3)
    # intersecting over time stops earlier and the coverage at the stop stays above 1 - alpha
    assert 1.4 < res["ratio"] < 1.7 and res["coverage"] > 0.95 and res["not_stopped"] == 0


def test_betting_oracle_tends_to_limit():
    z, p = S.Gaussian().z_law(3.0, k=1000)
    mu = float(np.sum(p * z))
    var = float(np.sum(p * (z - mu) ** 2))
    dz = 0.003
    n = S.betting_oracle_n(z, p, None, 0.05, d_z=dz)
    assert n / (S.zq(0.05) ** 2 * var / dz**2) == pytest.approx(S.cs_ratio_limit(0.05), rel=0.01)


# ----------------------------------------------------------------------------- engine against brute force


def _brute_cr(x, d, z, m, relative):
    for n in range(m, x.size + 1):
        xs = x[:n]
        v = xs.var(ddof=1)
        lhs = (z / d) ** 2 * (v + 1.0 / n)
        if relative:
            lhs /= xs.mean() ** 2
        if n >= lhs:
            return n, xs.mean()
    return x.size, x.mean()


@pytest.mark.parametrize("relative", [False, True])
def test_chow_robbins_stopping_matches_brute_force(relative):
    rng = np.random.default_rng(5)
    R = rng.standard_t(5, (25, 3000)) * 0.9
    st = S.Setup(d=0.35, alpha=0.1, m=12, relative=relative, n_max=3000, block=37,
                 procs=("CR-raw", "CR-trunc", "Stein-trunc"))
    out = S.run(S.Windows(R), st, 25, 0)
    z = S.zq(0.1)
    unit = out["_meta"]["unit"]
    tau = out["_meta"]["tau"]
    for i in range(25):
        u = R[i] / unit[i]
        n, mean = _brute_cr(u**2, 0.35, z, 12, relative)
        assert out["CR-raw"].N[i] == n
        assert out["CR-raw"].est[i] == pytest.approx(mean * unit[i] ** 2, rel=1e-9)
        xt = np.minimum(R[i] ** 2, tau[i] ** 2) / unit[i] ** 2
        n, mean = _brute_cr(xt, 0.35, z, 12, relative)
        assert out["CR-trunc"].N[i] == n
        # Stein: N = max(m, floor((t/d)^2 s^2) + 1) from the pilot
        p = xt[:12]
        dd = 0.35 * p.mean() if relative else 0.35
        ns = max(12, math.floor((stats.t.ppf(0.95, 11) / dd) ** 2 * p.var(ddof=1)) + 1)
        assert out["Stein-trunc"].N[i] == min(ns, 3000)
        assert out["Stein-trunc"].est[i] == pytest.approx(xt[: min(ns, 3000)].mean() * unit[i] ** 2, rel=1e-9)


def test_batch_means_rule_matches_brute_force():
    rng = np.random.default_rng(6)
    R = rng.standard_normal((6, 6000))
    st = S.Setup(d=0.15, alpha=0.05, m=30, n_max=6000, block=64, batch=50, procs=("CR-BM-raw",))
    out = S.run(S.Windows(R), st, 6, 0)
    z = S.zq(0.05)
    for i in range(6):
        x = R[i] ** 2
        want = 6000
        for n in range(500, 6001, 50):  # at least 10 complete batches
            b = x[:n].reshape(-1, 50).mean(axis=1)
            v = 50 * b.var(ddof=1)
            if n >= (z / 0.15) ** 2 * (v + 1.0 / n):
                want = n
                break
        assert out["CR-BM-raw"].N[i] == want


def test_eb_sequence_reproduces_package_without_intersection():
    rng = np.random.default_rng(7)
    x = np.minimum(rng.standard_normal((5, 400)) ** 2, 9.0)
    ref = eb_confseq(x, 9.0, 0.025, 0.025, n_plan=None, c=0.75)
    st = S.Setup(d=0.1, alpha=0.05, intersect=False)
    cs = S._CS(st, np.full((5, 3), 0.5), bets="wsr", init="wsr", dz=np.full(5, 1.0))
    cs.active[:] = True
    for b in range(0, 400, 50):  # several blocks: the carried state must be exact
        cs.update(x[:, b:b + 50] / 9.0, np.arange(5))
        cs.active[:] = True
        np.testing.assert_allclose(cs.U * 9.0, ref.upper[:, b + 49], rtol=1e-10)
        np.testing.assert_allclose(cs.L * 9.0, ref.lower[:, b + 49], rtol=1e-10, atol=1e-14)


def test_betting_bounds_are_conservative_and_tight():
    """Grid bounds with interpolation contain the exact capital-process bounds and are close to them."""
    rng = np.random.default_rng(8)
    z = rng.beta(0.6, 4.0, (3, 600))
    st = S.Setup(d=0.02, alpha=0.05, block=600, bet_grid=48, bet_span=5.0, intersect=False)
    cs = S._CS(st, z[:, :40], bets="bet", init="pilot", dz=np.full(3, 0.02))
    lu, ll = cs._gro_bets(np.arange(3), cs.zhat.copy())
    cs.update(z[:, 40:], np.arange(3))
    thr = math.log(40)
    zz = z[:, 40:]
    for i in range(3):
        f_up = lambda m: np.sum(np.log1p(lu[i] * (m - zz[i]))) - thr  # noqa: E731
        f_lo = lambda m: np.sum(np.log1p(ll[i] * (zz[i] - m))) - thr  # noqa: E731
        from scipy.optimize import brentq
        u_exact = brentq(f_up, cs.grid[i, 0], cs.grid[i, -1])
        l_exact = brentq(f_lo, cs.grid[i, 0], cs.grid[i, -1])
        assert cs.U[i] >= u_exact - 1e-12 and cs.U[i] <= u_exact * 1.0005
        assert cs.L[i] <= l_exact + 1e-12 and cs.L[i] >= l_exact * 0.9995


def test_confidence_sequences_cover_at_all_times():
    """Time-uniform coverage of the true mean of a bounded skewed variable (Monte Carlo, loose bound)."""
    rng = np.random.default_rng(9)
    reps, n = 400, 800
    z = rng.beta(0.5, 3.5, (reps, n))
    mu = 0.5 / 4.0
    for bets in ("opt", "bet"):
        st = S.Setup(d=1e-9, alpha=0.10, block=100)  # tiny d: the rule never fires, the bounds run to n
        cs = S._CS(st, z[:, :30], bets=bets, dz=np.full(reps, 1e-9))
        miss = np.zeros(reps, bool)
        for b in range(30, n, 100):
            cs.update(z[:, b:b + 100], np.arange(reps))
            miss |= (cs.U < mu) | (cs.L > mu)  # intersected bounds: checking at block ends covers all times
        assert miss.mean() <= 0.10 + 3 * math.sqrt(0.1 * 0.9 / reps)


def test_engine_smoke_gaussian_coverage_and_ratios():
    st = S.Setup(d=0.25, alpha=0.05, m=30, relative=True,
                 procs=("CR-trunc", "Oracle-trunc", "CS-opt", "CS-bet"))
    out = S.run(S.Gaussian(), st, 600, 11)
    sm = {p: S.summarise(out[p], st.d, True, st.alpha) for p in st.procs}
    assert sm["Oracle-trunc"]["EN_over_nstar"] == pytest.approx(1.0, abs=0.05)
    assert sm["CR-trunc"]["coverage"] > 0.88
    assert sm["CS-bet"]["coverage"] >= 0.95 and sm["CS-opt"]["coverage"] >= 0.95
    assert sm["CS-bet"]["EN"] < sm["CS-opt"]["EN"]
    assert 1.4 < sm["CS-bet"]["EN_over_nstar"] < 3.5


def test_garch_engine_records_conditional_targets():
    g = S.Garch(0.10, 0.85, burn=500)
    taus = np.linspace(0.5, 8, 40)
    gt = S.garch_target_table(g, taus, n=20_000, reps=2, lags=100)
    st = S.Setup(d=0.3, alpha=0.1, m=30, relative=True, procs=("CR-trunc", "CS-bet"))
    out = S.run(g, st, 50, 12, garch_targets=gt)
    for p in st.procs:
        assert np.all(np.isfinite(out[p].cond_target))
    assert np.all(np.isfinite(out["CS-bet"].cond_target_w))


def test_minimum_risk_rule_brute_force():
    rng = np.random.default_rng(13)
    R = rng.standard_normal((8, 4000))
    A = 200.0**2 / 2.0  # n0 = sqrt(A nu) = 200 for Gaussian r^2
    st = S.Setup(d=1.0, m=10, n_max=4000, block=50, minrisk_A=A, procs=("MinRisk-raw",))
    out = S.run(S.Windows(R), st, 8, 0)
    for i in range(8):
        x = R[i] ** 2
        want = 4000
        for n in range(10, 4001):
            if n >= math.sqrt(A) * (x[:n].std(ddof=1) + 1.0 / n):
                want = n
                break
        assert out["MinRisk-raw"].N[i] == want
