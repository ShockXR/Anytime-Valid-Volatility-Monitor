"""Tests for avmonitor.design: growth rates, information numbers, delay bounds and approximations.

Theory: docs/theory_delay.md. All checks are against Monte Carlo or exact identities.
"""

import numpy as np
import pytest

from avmonitor import design as G
from avmonitor import detectors as D
from avmonitor.robust import truncated_moment_factor

MODELS = [G.ReturnModel("gaussian"), G.ReturnModel("jumps"), G.ReturnModel("student_t", df=5.0)]


def test_gaussian_baseline_matches_truncated_moment():
    m = G.ReturnModel("gaussian")
    for c in (2.0, 3.0, 4.0):
        tau, m0 = G.baseline(m, c)
        assert tau == pytest.approx(c, rel=1e-9)  # scaled MAD of N(0, 1) is 1
        assert m0 == pytest.approx(truncated_moment_factor(c, 2), rel=1e-9)
    tau, m0 = G.baseline(m, np.inf)
    assert np.isinf(tau) and m0 == pytest.approx(1.0, rel=1e-9)


@pytest.mark.parametrize("model", MODELS, ids=lambda m: m.kind)
def test_quadrature_law_is_a_probability_law_with_unit_null_mean(model):
    tau, m0 = G.baseline(model, 3.0)
    law = G.y_law(model, 1.0, tau, m0)
    assert law.w.sum() == pytest.approx(1.0, abs=1e-9)
    assert law.mean == pytest.approx(1.0, abs=1e-9)
    assert law.y.max() <= law.B * (1 + 1e-12)


@pytest.mark.parametrize("model", MODELS, ids=lambda m: m.kind)
@pytest.mark.parametrize("direction,rho", [("up", 1.5), ("down", 0.67)])
def test_growth_rate_matches_monte_carlo(model, direction, rho):
    rng = np.random.default_rng(101)
    tau, m0 = G.baseline(model, 3.0)
    law = G.y_law(model, rho, tau, m0)
    lams = G.lam_max(direction, law.B) * np.array([0.1, 0.4, 0.8]) if direction == "down" else np.array([0.05, 0.2, 0.5])
    g = G.growth_rate(lams, law, direction)
    y = G.sample_y(model, rho, tau, m0, 400_000, rng)
    gain = (y - 1.0) if direction == "up" else (1.0 - y)
    z = np.log1p(np.outer(lams, gain))
    mc, se = z.mean(axis=1), z.std(axis=1) / np.sqrt(y.size)
    assert np.all(np.abs(g - mc) < 4 * se + 1e-6)


def test_information_ladder_and_gaussian_kl():
    m = G.ReturnModel("gaussian")
    for rho in (0.5, 1.25, 2.0):
        assert G.kl_truncated(m, rho, np.inf) == pytest.approx(float(G.kl_gaussian_scale(rho)), rel=1e-8)
    for model in MODELS:
        for c in (2.0, 3.0, 5.0):
            tau, m0 = G.baseline(model, c)
            for direction, rho in (("up", 1.25), ("up", 2.0), ("down", 0.5), ("down", 0.8)):
                law = G.y_law(model, rho, tau, m0)
                lam_o, kl_inf = G.best_bet(law, direction)
                grid = np.array(D.DEFAULT_UP_GRID) if direction == "up" else \
                    G.lam_max("down", law.B) * np.array(D.DEFAULT_DOWN_FRACTIONS)
                g_grid = G.growth_rate(grid, law, direction).max()
                kl_x = G.kl_truncated(model, rho, tau)
                kl_r = G.kl_truncated(model, rho, np.inf)
                assert kl_r >= kl_x - 1e-10 >= kl_inf - 2e-10 >= g_grid - 3e-10 > 0
                ceiling = np.log(law.B) if direction == "up" else np.log(law.B / (law.B - 1))
                assert kl_inf <= ceiling


@pytest.mark.parametrize("direction,rho", [("up", 1.5), ("up", 1.1), ("down", 0.9)])
def test_kl_inf_duality_certificate(direction, rho):
    """At an interior optimum, P*(y) = Q(y) / (1 + lam* gain(y)) is a null law with KL(Q || P*) = g*."""
    m = G.ReturnModel("gaussian")
    tau, m0 = G.baseline(m, 3.0)
    law = G.y_law(m, rho, tau, m0)
    lam, g = G.best_bet(law, direction)
    assert 0 < lam < 0.99 * G.lam_max(direction, law.B)
    gain = (law.y - 1.0) if direction == "up" else (1.0 - law.y)
    p_star = law.w / (1.0 + lam * gain)
    assert p_star.sum() == pytest.approx(1.0, abs=1e-6)
    assert np.sum(p_star * law.y) == pytest.approx(1.0, abs=1e-6)  # on the boundary of the null
    kl = np.sum(law.w * np.log(law.w / p_star))
    assert kl == pytest.approx(g, abs=1e-9)


def test_run_mixture_matches_betting_edetector():
    rng = np.random.default_rng(5)
    m = G.ReturnModel("gaussian")
    tau, m0 = G.baseline(m, 3.0)
    y = G.sample_y(m, 1.5, tau, m0, (300, 600), rng)
    a1 = G.run_mixture(y, D.DEFAULT_UP_GRID, np.full(6, 1 / 6), "up", 500.0)["alarm"]
    a2 = D.betting_edetector(y, 1.0, "up", threshold=500.0).alarm
    np.testing.assert_array_equal(a1, a2)
    B = tau**2 / m0
    grid = G.lam_max("down", B) * np.array(D.DEFAULT_DOWN_FRACTIONS)
    y = G.sample_y(m, 0.6, tau, m0, (300, 600), rng)
    a1 = G.run_mixture(y, grid, np.full(5, 0.2), "down", 500.0)["alarm"]
    a2 = D.betting_edetector(y, 1.0, "down", upper=B, threshold=500.0).alarm
    np.testing.assert_array_equal(a1, a2)


def test_delay_bound_and_second_order_approximation_easy_case():
    """Single growth-optimal bet, sd doubles, change at the first step, A = 7800.

    Stated tolerance: the second-order formula (log A + kappa - zeta) / g is within 6% of the
    simulated mean delay; the rigorous Wald bound is never below it; log(A) / g is within 25%."""
    rng = np.random.default_rng(2026)
    m = G.ReturnModel("gaussian")
    tau, m0 = G.baseline(m, 3.0)
    A = 7800.0
    law = G.y_law(m, 2.0, tau, m0)
    lam, g = G.best_bet(law)
    rc = G.renewal_constants(lam, lambda s, r: G.sample_y(m, 2.0, tau, m0, s, r), "up", 3000, rng=rng)
    d2 = G.delay_second_order(A, 1.0, g, rc["kappa"], rc["zeta"])
    ub, _ = G.wald_delay_bound(A, [lam], [1.0], law)
    y = G.sample_y(m, 2.0, tau, m0, (3000, 400), rng)
    alarm = G.run_mixture(y, [lam], [1.0], "up", A)["alarm"]
    assert np.all(alarm >= 0)
    sim = (alarm + 1).mean()
    assert abs(d2 - sim) / sim < 0.06
    assert ub >= sim
    assert abs(G.delay_first_order(A, g) - sim) / sim < 0.25


def test_zeta_brownian_matches_random_walk_memory_for_small_drift():
    g = 0.01
    rng = np.random.default_rng(9)
    rc = G.renewal_constants_increments(lambda s, r: g + np.sqrt(2 * g) * r.standard_normal(s), 3000, 200000, rng)
    assert rc["zeta"] == pytest.approx(G.zeta_brownian(g, 2 * g), abs=0.12)
    # at a growth-optimal bet (variance = 2 g) the memory term is log(1 / g) + Euler's gamma
    assert G.zeta_brownian(g, 2 * g) == pytest.approx(np.log(1 / g) + G.EULER_GAMMA, abs=1e-12)


def test_speed_limit_is_attained_by_extreme_data():
    B = 9.0
    A = 7800.0
    grid = G.lam_max("down", B) * np.array(D.DEFAULT_DOWN_FRACTIONS)
    w = np.full(len(grid), 1 / len(grid))
    n_min = G.speed_limit(A, 1 + grid, w)
    alarm = G.run_mixture(np.zeros((1, 500)), grid, w, "down", A)["alarm"][0]
    assert alarm + 1 == n_min
    up = np.array(D.DEFAULT_UP_GRID)
    n_up = G.speed_limit(A, 1 + up * (B - 1), np.full(6, 1 / 6))
    assert G.run_mixture(np.full((1, 50), B), up, np.full(6, 1 / 6), "up", A)["alarm"][0] + 1 == n_up


def test_worst_case_dp_one_step_and_bound():
    B = 9.0
    lam, A = 0.5, 3.0
    lo, up = G.worst_case_false_alarm(lam, B, A, 1, n_states=400)
    exact = 1.0 / (1.0 + (A - 1.0) / lam)  # Markov's inequality is attained by a two-point law
    assert lo[0] <= exact + 1e-12 <= up[0] + 2e-12
    assert up[0] - lo[0] < 0.01
    lo, up = G.worst_case_false_alarm(0.5, B, 200.0, 40, n_states=300)
    m = np.arange(1, 41)
    assert np.all(lo <= m / 200.0 + 1e-12)
    assert np.all(np.diff(lo) >= -1e-12)
