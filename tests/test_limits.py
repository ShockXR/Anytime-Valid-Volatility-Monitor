"""Tests for src/avmonitor/limits.py (TH5, docs/theory_limits.md)."""

import numpy as np
import pytest
from scipy import optimize, stats

from avmonitor import design as D
from avmonitor import limits as L
from avmonitor import robustness as R
from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID

A = 7800.0


def _gauss(c=3.0):
    model = D.ReturnModel("gaussian")
    tau, m0 = D.baseline(model, c)
    return model, tau, m0, tau**2 / m0


# ----------------------------------------------------------------------------
# Part 1: the floor lemma
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("side", ["down", "up"])
def test_max_null_mass_is_the_linear_programme_optimum(side):
    """Largest P(Y = y) over laws on a grid of [0, B] with the null mean constraint (linear programme)."""
    B = 9.0
    grid = np.linspace(0, B, 91)
    for y in (0.0, 0.3, 0.7, 0.95, 1.0, 1.5, 4.0, 9.0):
        support = np.unique(np.append(grid, y))
        j = int(np.flatnonzero(support == y)[0])
        c = np.zeros(support.size)
        c[j] = -1.0
        sign = -1.0 if side == "down" else 1.0  # down: mean >= 1  <=>  -mean <= -1
        res = optimize.linprog(c, A_ub=[sign * support], b_ub=[sign * 1.0], A_eq=[np.ones(support.size)], b_eq=[1.0],
                               bounds=(0, 1), method="highs")
        assert res.status == 0
        assert -res.fun == pytest.approx(float(L.max_null_mass(y, B, side)), abs=1e-9)


def test_two_point_laws_attain_the_mass_and_lie_in_the_null():
    B = 9.0
    for y in (0.0, 0.4, 0.99):
        p = float(L.max_null_mass(y, B, "down"))
        assert p * y + (1 - p) * B == pytest.approx(1.0)
    for y in (1.01, 3.0, 9.0):
        p = float(L.max_null_mass(y, B, "up"))
        assert p * y == pytest.approx(1.0)


def test_unbounded_decrease_side_gives_no_floor_growth():
    """B = inf: every value has null mass 1 (docs/math.md section 5.3), so W_n = 1 and no rule with
    sup P(T <= n) < 1 can alarm on any path."""
    y = np.array([0.0, 0.1, 0.5, 2.0])
    assert np.all(L.max_null_mass(y, np.inf, "down") == 1.0)
    assert np.all(L.floor_time(L.log_floor_factor(np.zeros((1, 500)), np.inf, "down"), A) == -1)


def test_product_law_probability_of_a_path_is_one_over_W():
    """Under the product of the two-point laws, P(Y_{1:n} = y) = 1 / W_n(y) exactly (Monte Carlo on a short
    path of three values)."""
    B = 9.0
    path = np.array([0.0, 0.5, 0.2])
    p = L.max_null_mass(path, B, "down")
    rng = np.random.default_rng(1)
    N = 400_000
    hit = np.ones(N, bool)
    for s, ys in enumerate(path):
        draw = np.where(rng.random(N) < p[s], ys, B)
        assert abs(np.mean(np.where(draw == ys, ys, B)) - 1.0) < 0.02  # mean 1 at every step
        hit &= draw == ys
    W = np.exp(np.sum(L.log_floor_factor(path, B, "down")))
    assert np.mean(hit) == pytest.approx(1 / W, abs=4 * np.sqrt(1 / W / N))


def test_collapse_floors():
    _, _, _, B3 = _gauss(3.0)
    _, _, _, B2 = _gauss(2.0)
    assert L.collapse_floor(A, B3, "down") == 45
    assert L.collapse_floor(2 * A, B3, "down") == 50
    assert L.collapse_floor(A, B2, "down") == 23
    assert L.collapse_floor(A, B3, "up") == 4
    # the floors of Corollary 4 of docs/theory_delay.md (Shiryaev-Roberts with linear bets, default grid)
    F = B3 / (B3 - 1)
    assert np.log1p(A * (F - 1) / F) / np.log(F) == pytest.approx(57.7, abs=0.05)
    # at Y = 0 the decrease bet phi = fraction / (B - 1) earns the factor 1 + phi
    assert D.speed_limit(A, 1 + np.asarray(DEFAULT_DOWN_FRACTIONS) / (B3 - 1), np.full(5, 0.2)) == 73


def test_randomised_rules_on_the_collapse():
    """Corollary 4 is for deterministic rules. A randomised valid rule alarms by step n on the collapse with
    probability at most min(1, n W_n / A) (0.978 at step 44), and its smallest mean delay there is 43.4
    (linear programme); on the explosion the largest probability at step 3 is 0.285."""
    _, _, _, B = _gauss(3.0)
    F = B / (B - 1)
    assert min(1.0, 44 * F**44 / A) == pytest.approx(0.978, abs=1e-3)
    lw = np.arange(1, 301) * np.log(F)
    mean, q = L.randomised_path_delay(lw, A)
    assert mean == pytest.approx(43.39, abs=0.01) and 43 < mean < 45
    j = np.arange(1, q.size + 1)
    assert q.sum() == pytest.approx(1.0, abs=1e-9) and np.all(q >= -1e-12)
    assert np.all(np.cumsum(q * np.exp(-lw[: q.size])) <= j / A + 1e-12)  # validity on every k
    assert q[41] == pytest.approx(0.738, abs=2e-3)  # most of the mass at step 42
    # the deterministic rule (all mass at the floor) is feasible and slower
    assert mean < L.collapse_floor(A, B, "down")
    assert min(1.0, 3 * B**3 / A) == pytest.approx(0.285, abs=1e-3)
    assert L.randomised_path_delay(np.arange(1, 41) * np.log(B), A)[0] == pytest.approx(3.715, abs=1e-3)


def test_singleton_rule_is_valid_and_attains_the_floor():
    """Sharpness: the rule 'alarm at n iff Y_{1:n} = y' has sup P(T <= k) = 1 / W_n(y) (chain rule), so it is
    valid at levels k / A exactly when W_n(y) >= A / n; it then alarms at the floor time."""
    _, _, _, B = _gauss(3.0)
    path = np.zeros(200)
    n = int(L.floor_time(L.log_floor_factor(path, B, "down"), A)[0])
    W_n = (B / (B - 1)) ** n
    assert W_n >= A / n
    assert 1 / W_n <= n / A  # the singleton rule is valid
    assert (B / (B - 1)) ** (n - 1) < A / (n - 1)  # and nothing valid alarms one step earlier


@pytest.mark.parametrize("model_name,rho", [("gaussian", 0.5), ("gaussian", 2.0), ("student t3", 0.8), ("jumps", 1.25)])
def test_e_detector_never_alarms_before_the_floor(model_name, rho):
    model = {"gaussian": D.ReturnModel("gaussian"), "student t3": D.ReturnModel("student_t", df=3.0),
             "jumps": D.ReturnModel("jumps")}[model_name]
    tau, m0 = D.baseline(model, 3.0)
    B = tau**2 / m0
    side = "up" if rho > 1 else "down"
    rng = np.random.default_rng(7)
    y = D.sample_y(model, rho, tau, m0, (300, 1500), rng)
    lw = L.log_floor_factor(y, B, side)
    law = D.y_law(model, rho, tau, m0)
    lam, _ = D.best_bet(law, side)
    bets = np.asarray(DEFAULT_UP_GRID) if side == "up" else np.asarray(DEFAULT_DOWN_FRACTIONS) / (B - 1)
    for AA in (A, 2 * A):
        nf = L.floor_time(lw, AA)
        for bb, ww in ((bets, np.full(bets.size, 1 / bets.size)), ([lam], [1.0]), ([bets.max()], [1.0])):
            t = L.sr_alarm_times(y, bb, ww, side, (AA,))[0]
            ok = t > 0
            assert np.all(nf[ok] > 0)
            assert np.all(t[ok] >= nf[ok])


def test_randomised_bound_and_neyman_pearson():
    # for a single deterministic path the randomised bound is min(1, (n / A) W_n)
    B = 9.0
    path = np.zeros((1, 80))
    rf = L.randomised_floor_cdf(L.log_floor_factor(path, B, "down"), A)
    n = np.arange(1, 81)
    assert np.allclose(rf, np.minimum(1, n / A * (B / (B - 1)) ** n))
    # Neyman-Pearson against the closed form: Q = N(mu, 1), P* = N(0, 1)
    rng = np.random.default_rng(3)
    mu, alpha = 1.5, 0.05
    x = rng.normal(mu, 1.0, 400_000)
    s = mu * x - mu**2 / 2
    exact = stats.norm.sf(stats.norm.isf(alpha) - mu)
    assert L.np_power(s, alpha) == pytest.approx(exact, abs=0.01)
    # a cap can only lower the power
    cap = (x > 1.0).astype(float)
    assert L.np_power(s, alpha, cap) <= L.np_power(s, alpha) + 1e-12


def test_np_bound_is_below_simulated_alarm_cdf_of_valid_rules():
    """beta_n bounds P(T <= n) for every valid rule: the e-detectors' empirical CDF stays below it."""
    model, tau, m0, B = _gauss(3.0)
    law = D.y_law(model, 0.5, tau, m0)
    lam, _ = L.least_favourable_bet(law, "down")
    rng = np.random.default_rng(11)
    y = D.sample_y(model, 0.5, tau, m0, (2000, 400), rng)
    lw = L.log_floor_factor(y, B, "down")
    S = np.cumsum(np.log1p(lam * (1 - y)), axis=1)
    beta = L.np_cdf(S, A, lw, "det")
    t = L.sr_alarm_times(y, [lam], [1.0], "down", (A,))[0]
    t = np.where(t < 0, 10**9, t)
    emp = np.array([np.mean(t <= n) for n in range(1, beta.size + 1)])
    assert np.all(emp <= beta + 3 * np.sqrt(np.maximum(beta * (1 - beta), 1e-4) / 2000) + 1e-9)
    assert L.delay_lower_bound(beta) <= np.mean(np.minimum(t, 400)) + 1


# ----------------------------------------------------------------------------
# Part 2: the plug-in baseline
# ----------------------------------------------------------------------------


def test_drift_is_bounded_by_B_hat_times_kolmogorov_distance():
    rng = np.random.default_rng(5)
    for law in ("gaussian", "t3"):
        for n in (50, 390):
            win = L.plugin_windows(law, n, 400, rng, 3.0, (), kolmogorov=True)
            assert np.all(np.abs(win["d"]) <= win["B_hat"] * win["D"] + 1e-12)


def test_B_hat_range_for_any_sample():
    lo, hi = L.b_range(3.0)
    rng = np.random.default_rng(9)
    a = 3.0 * 1.4826
    samples = [rng.standard_cauchy((200, 51)), rng.standard_normal((200, 50)), rng.exponential(size=(200, 40)),
               np.where(rng.random((200, 60)) < 0.5, 1.0, 50.0), rng.standard_t(1.5, (200, 391))]
    for r in samples:
        med = np.median(np.abs(r), axis=1)
        tau2 = (a * med) ** 2
        m = np.minimum(r * r, tau2[:, None]).mean(axis=1)
        Bh = tau2 / m
        assert np.all(Bh >= lo - 1e-9) and np.all(Bh <= hi + 1e-9)


def test_growth_sum_bound():
    """G_m(x) = S_m(x) / m <= 1 + x (m + 1) e^(x m) / 2."""
    for m in (1, 5, 98, 390, 1950):
        x = np.geomspace(1e-6, 0.05, 40)
        G = R.sr_growth_sum(x, m) / m
        assert np.all(G <= 1 + x * (m + 1) * np.exp(x * m) / 2 + 1e-12)


def test_distribution_free_bound_is_finite_and_tends_to_nominal():
    k1, k2 = L.distribution_free_constants(3.0)
    assert k1 == pytest.approx(8.045, abs=1e-3) and k2 == pytest.approx(19.783, abs=1e-3)
    r = [L.distribution_free_bound(n, 10, A) / (10 / A) for n in (10**5, 10**6, 10**7, 10**8)]
    assert all(np.isfinite(r)) and r == sorted(r, reverse=True)
    # first-order rate: (ratio - 1) / ((m + 1) / sqrt(n)) -> C_DF = k1 I0 / 4 = 1.79 (m6 of review round 7)
    cdf = L.distribution_free_rate_constant(3.0)
    assert cdf == pytest.approx(1.787, abs=1e-3)
    excess = [(L.distribution_free_bound(n, 10, A) / (10 / A) - 1) / (11 / np.sqrt(n)) for n in (10**8, 10**10, 10**12)]
    assert excess == sorted(excess, reverse=True) and excess[-1] == pytest.approx(cdf, rel=1e-3)
    # Monte Carlo check of the moment bound: E[D exp(t D)] for Gaussian |r| at n = 390
    rng = np.random.default_rng(2)
    a = np.sort(np.abs(rng.standard_normal((4000, 390))), axis=1)
    Dn = L.kolmogorov_distance(a, "gaussian")
    t = 20.0
    assert np.mean(Dn * np.exp(t * Dn)) <= L.dkw_moment(390, t)


@pytest.mark.parametrize("n,m", [(390, 5), (92160, 39), (10**6, 98), (10**7, 98), (10**7, 390), (10**8, 390),
                                 (10**8, 1950)])
def test_dkw_moment_against_closed_form(n, m):
    """With u = sqrt(n) x and s = t / sqrt(n), the integral of (1 + t x) e^(t x) min(1, 2 e^(-2 n x^2)) over
    [0, 1] is n^(-1/2) [u0 e^(s u0) + 2 e^(s^2 / 8) ((1 + s^2 / 4) sqrt(pi / 2) (Phi(2 vb) - Phi(2 va))
    + (s / 4)(e^(-2 va^2) - e^(-2 vb^2)))] with va = u0 - s / 4, vb = sqrt(n) - s / 4. A first version of
    dkw_log_moment under-integrated for n >= 10^6 (review round 7, M2)."""
    _, k2 = L.distribution_free_constants(3.0)
    t = k2 * m
    rn = np.sqrt(n)
    s = t / rn
    u0 = np.sqrt(np.log(2) / 2)
    va, vb = u0 - s / 4, rn - s / 4
    gauss = np.sqrt(np.pi / 2) * (stats.norm.cdf(2 * vb) - stats.norm.cdf(2 * va))
    tail = 2 * np.exp(s * s / 8) * ((1 + s * s / 4) * gauss + (s / 4) * (np.exp(-2 * va**2) - np.exp(-2 * vb**2)))
    exact = np.log(u0 * np.exp(s * u0) + tail) - np.log(rn)
    assert L.dkw_log_moment(n, t) == pytest.approx(exact, rel=1e-8, abs=1e-8)


def test_distribution_free_table_values():
    """Values quoted in docs/theory_limits.md section 4.7 (ratio to m / A)."""
    ratio = lambda n, m: L.distribution_free_bound(n, m, A) / (m / A)
    assert ratio(10**6, 39) == pytest.approx(1.1525, abs=1e-3)
    assert ratio(10**8, 98) == pytest.approx(1.0213, abs=1e-3)
    assert ratio(10**8, 390) == pytest.approx(1.1491, abs=1e-3)
    assert ratio(10**7, 390) == pytest.approx(3.9409, abs=1e-3)


def test_first_order_constant():
    pop = L.law_population("gaussian", 3.0)
    assert pop["sigma_y"] == pytest.approx(1.38585, abs=1e-4)
    assert pop["B"] == pytest.approx(9.04516, abs=1e-4)
    assert L.first_order_constant(pop["sigma_y"], pop["B"]) == pytest.approx(0.0369, abs=2e-4)
    t3 = L.law_population("t3", 3.0)
    assert t3["B"] == pytest.approx(6.36248, abs=1e-4) and t3["sigma_y"] == pytest.approx(1.58267, abs=1e-4)


def test_sqrt_n_drift_is_asymptotically_normal_with_sd_sigma_y():
    rng = np.random.default_rng(4)
    win = L.plugin_windows("gaussian", 5760, 3000, rng, 3.0, ())
    sd = float(np.std(win["d"]) * np.sqrt(5760))
    assert sd == pytest.approx(1.38585, rel=0.06)


def _hold_bold_sim(d, Bh, m, A_, Lb, rng, P):
    """Hold Y = 1 + d for m - Lb steps, then the two-point law on {0, B_hat} with mean 1 + d; default grid,
    each side at 2 A_. Returns the share of paths with an alarm within m."""
    lam = np.asarray(DEFAULT_UP_GRID)
    phi = np.asarray(DEFAULT_DOWN_FRACTIONS) / (Bh - 1)
    Su = np.zeros((P, lam.size))
    Sd = np.zeros((P, phi.size))
    alarm = np.zeros(P, bool)
    for t in range(1, m + 1):
        Y = np.full(P, 1 + d) if t <= m - Lb else np.where(rng.random(P) < (1 + d) / Bh, Bh, 0.0)
        Su = (1 + lam * (Y - 1)[:, None]) * (Su + 1)
        Sd = (1 + phi * (1 - Y)[:, None]) * (Sd + 1)
        alarm |= np.maximum(Su.mean(axis=1), Sd.mean(axis=1)) >= 2 * A_
    return alarm.mean()


@pytest.mark.parametrize("d", [0.04, -0.05])
def test_hold_bold_lower_bound_is_below_simulation_and_upper_bound(d):
    Bh, m, A_ = 9.0, 120, 60.0
    rng = np.random.default_rng(13)
    best, bestL = 0.0, 0
    for Lb in range(1, m + 1):
        v = float(L.hold_bold_lower_bound(np.array([d]), np.array([Bh]), m, A_, holds=[Lb])[0])
        if v > best:
            best, bestL = v, Lb
    assert best > 0.05  # the bound is informative here
    P = 20000
    rate = _hold_bold_sim(d, Bh, m, A_, bestL, rng, P)
    assert best <= rate + 3 * np.sqrt(rate * (1 - rate) / P)
    # and below the upper bound of Corollary 1 (TH4), which covers every stream with this conditional mean
    m_hat = 1.0
    up = L.side_bound(np.array([m_hat * (1 + d)]), np.array([m_hat]), np.array([m_hat]), np.array([Bh * m_hat]), m, A_)
    assert best <= float(up[0]) + 1e-12


def test_inflated_baseline_shares():
    rng = np.random.default_rng(6)
    win = L.plugin_windows("gaussian", 2880, 4000, rng, 3.0, (2.576,))
    bad = (win["m_true"] > win["up_2.576"]) | (win["m_true"] < win["dn_2.576"])
    assert np.mean(bad) < 0.03


# ----------------------------------------------------------------------------
# Part 3: the decrease bet near its boundary
# ----------------------------------------------------------------------------


def test_boundary_terms_gaussian_half():
    model, tau, m0, B = _gauss(3.0)
    law = D.y_law(model, 0.5, tau, m0)
    t = L.boundary_terms(law)
    assert 1e-8 < t["eps_star"] < 5e-8  # the referee's 1 - 2.4e-8
    lam, kl = D.best_bet(law, "down")
    assert t["kl_inf"] == pytest.approx(kl, rel=1e-6)
    eps = np.array([0.5, 0.1, 1e-3, 1e-6, t["eps_star"]])
    g = L.boundary_growth(law, eps)
    lo, hi = L.boundary_growth_bounds(law, eps)
    assert np.all(lo <= g + 1e-12) and np.all(g <= hi + 1e-12)
    assert t["G_c"] - t["q_B"] * (1 + np.log(t["G_c"] / t["q_B"])) <= L.boundary_growth(law, t["q_B"] / t["G_c"])[0] + 1e-12


def test_interior_optimum_iff_mean_of_inverse_factor_exceeds_one():
    """lam* < lam_max exactly when E_Q[(B - 1) / (B - Y)] > 1 (an atom at B makes it infinite)."""
    B = 9.0
    for ys, ws in (([0.0, 0.5], [0.5, 0.5]), ([0.0, 2.0], [0.9, 0.1]), ([0.0, 8.0], [0.5, 0.5]), ([0.2, B], [0.99, 0.01])):
        law = D.YLaw(y=np.array(ys, float), w=np.array(ws, float), B=B)
        with np.errstate(divide="ignore"):
            crit = float(np.sum(law.w * (B - 1) / (B - law.y)))
        grid = np.linspace(0, 1 - 1e-12, 200001)
        g = L.boundary_growth(law, 1 - grid)  # eps = 1 - lam / lam_max
        interior = int(np.argmax(g)) < grid.size - 1
        assert interior == (crit > 1)


def test_block_bound_is_above_the_simulated_delay():
    """E[T_A] <= n / p_n for one Shiryaev-Roberts bet (Proposition 12 (e)), here for the t3 law at sd ratio 0.8
    with a bet so close to the boundary that its growth rate is negative."""
    model = D.ReturnModel("student_t", df=3.0)
    tau, m0 = D.baseline(model, 3.0)
    B = tau**2 / m0
    law = D.y_law(model, 0.8, tau, m0)
    eps = 1e-3
    assert L.boundary_growth(law, eps)[0] < 0
    rng = np.random.default_rng(21)
    y = D.sample_y(model, 0.8, tau, m0, (1500, 4000), rng)
    t = L.sr_alarm_times(y, [(1 - eps) / (B - 1)], [1.0], "down", (A,))[0]
    assert np.all(t > 0)
    bound, _ = L.block_delay_bound(np.log(eps + (1 - eps) * (B - y) / (B - 1)), A)
    assert np.isfinite(bound) and t.mean() <= bound


@pytest.mark.parametrize("model_name,rho", [("gaussian", 0.5), ("gaussian", 0.8), ("student t3", 0.8), ("gaussian", 2.0)])
def test_least_favourable_law_is_a_probability_in_the_null(model_name, rho):
    model = D.ReturnModel("gaussian") if model_name == "gaussian" else D.ReturnModel("student_t", df=3.0)
    tau, m0 = D.baseline(model, 3.0)
    law = D.y_law(model, rho, tau, m0)
    side = "up" if rho > 1 else "down"
    lam, kl = L.least_favourable_bet(law, side)
    U = (law.y - 1) if side == "up" else (1 - law.y)
    p_star = law.w / (1 + lam * U)
    assert p_star.sum() == pytest.approx(1.0, abs=1e-8)
    assert abs(np.sum(p_star * U)) < 1e-8  # mean exactly 1
    assert kl == pytest.approx(D.best_bet(law, side)[1], rel=1e-7)
