"""TH4 (src/avmonitor/robustness.py, docs/theory_robustness.md): bounds when the null holds only approximately.

* the closed form S_m(x) and the elementary bounds on G_m(x) = S_m(x) / m;
* Lemma 1 (expected Shiryaev-Roberts wealth under a drift) by Monte Carlo, also at a stopping time;
* Theorem 1 on simulated streams with a known constant drift (the simulation check), including the
  deterministic stream at which the bound is attained;
* the one-pass engine run_paths agrees with the package's own detector, and Lemma 2 (N_t <= D*_t N~_t) holds
  path by path;
* drift_sr_sums, drift_cusum and threshold_bound agree with brute force;
* Theorems 2 and 3 and Proposition 4 hold on an exogenous switching stream; the counterexamples of section 3;
* Corollary 1 for the plug-in: per window, the capped bound is at least the alarm indicator of the
  deterministic stream that keeps the same conditional mean (the shape adversary of E16 with q = 0).
"""

import numpy as np
import pytest

from avmonitor import robustness as R
from avmonitor.detectors import DEFAULT_UP_GRID, betting_edetector, two_sided_alarm
from avmonitor.robust import MAD_NORMAL_SCALE

B = R.gaussian_B()


def test_gaussian_B():
    assert B == pytest.approx(9.045, abs=5e-4)


def test_growth_sum_closed_form():
    for x in (0.0, 1e-9, 1e-4, 0.01, 0.3, -0.2):
        for m in (1, 2, 7, 390, 1950):
            direct = sum((1 + x) ** j for j in range(1, m + 1))
            assert float(R.sr_growth_sum(x, m)) == pytest.approx(direct, rel=1e-9)
    assert float(R.growth_factor(0.0, 390)) == 1.0
    with pytest.raises(ValueError):
        R.sr_growth_sum(-1.0, 3)


def test_growth_factor_elementary_bounds():
    for x in np.linspace(0, 0.5, 26):
        for m in (1, 2, 10, 98, 390, 1950):
            g = float(R.growth_factor(x, m))
            lo, hi = R.growth_factor_bounds(x, m)
            assert lo <= g * (1 + 1e-12)
            assert g <= hi * (1 + 1e-12)
            assert 1 + x * (m + 1) / 2 <= lo * (1 + 1e-12)


def test_sup_bound_reduces_to_theorem_B():
    for m in (1, 390, 1950):
        assert float(R.sup_bound(0, 0, m, 7800, B)) == pytest.approx(m / 7800)
        assert float(R.sup_bound(0, 0, m, 7800, B, side="up")) == pytest.approx(m / 7800)
        assert float(R.sup_inflation(0, 0, m, B)) == pytest.approx(1.0)
    # increasing in the drift and in the horizon
    v = [float(R.sup_inflation(d, d, 390, B)) for d in (0, 0.001, 0.01, 0.02)]
    assert np.all(np.diff(v) > 0)


def test_first_order_and_upper_bound_of_the_inflation():
    lam, phi = R.default_bets(B)
    for d in (1e-4, 1e-3, 5e-3):
        for m in (98, 390):
            infl = float(R.sup_inflation(d, d, m, B))
            first = 1 + (m + 1) / 4 * (lam.mean() * d + phi.mean() * d)
            xs = max(lam.max(), phi.max()) * d
            assert first <= infl * (1 + 1e-12)
            assert infl <= np.exp(xs * (m + 1) / 2 + xs**2 * m**2 / 24) * (1 + 1e-12)


def test_design_drift_hits_target():
    for direction in ("both", "up", "down"):
        for tgt in (0.10, 0.25):
            d = R.design_drift(tgt, 390, 7800, B, direction)
            du = d if direction != "down" else 0.0
            dd = d if direction != "up" else 0.0
            assert float(R.sup_bound(du, dd, 390, 7800, B)) == pytest.approx(tgt, rel=1e-8)
    # the decrease side tolerates more drift (phi_max = 0.95 / (B - 1) < lam_max = 0.5)
    assert R.design_drift(0.1, direction="down") > 3 * R.design_drift(0.1, direction="up")


def test_lemma1_expected_wealth_and_stopping():
    """E[R_t] = S_t(lam delta) when E[f_t | past] = 1 + lam delta exactly; E[R_{T ^ t}] <= S_t at a stopping time."""
    rng = np.random.default_rng(1)
    lam, delta, n, P = 0.35, 0.03, 40, 200000
    Y = np.where(rng.random((P, n)) < (1 + delta) / B, B, 0.0)
    Rt = np.zeros(P)
    Rs = np.zeros(P)
    stopped = np.zeros(P, bool)
    for t in range(n):
        f = 1 + lam * (Y[:, t] - 1)
        Rt = f * (Rt + 1)
        Rs = np.where(stopped, Rs, f * (Rs + 1))
        stopped |= Rs >= 30  # a stopping time
    S = float(R.sr_growth_sum(lam * delta, n))
    se = Rt.std() / np.sqrt(P)
    assert abs(Rt.mean() - S) < 4 * se
    assert Rs.mean() <= S + 4 * Rs.std() / np.sqrt(P)


def _run_default(Y, mean, A):
    res = R.run_paths(Y, np.full(Y.shape, mean - 1.0), B, A, horizons=(Y.shape[1],))
    return res["alarm"]


@pytest.mark.parametrize("law", ["two-point", "gaussian"])
@pytest.mark.parametrize("direction", ["up", "down"])
def test_simulation_check_constant_drift(law, direction):
    """Theorem 1 on streams with a known constant drift (small A, so that the rates are not tiny)."""
    rng = np.random.default_rng(7 + (law == "gaussian") + 2 * (direction == "down"))
    A, m, P = 150.0, 60, 12000
    d = 0.05 if direction == "up" else 0.15
    mean = 1 + d if direction == "up" else 1 - d
    if law == "two-point":
        Y = np.where(rng.random((P, m)) < mean / B, B, 0.0)
    else:
        m0 = float(R.trunc_mean_normal(3.0))
        from scipy import optimize

        s = optimize.brentq(lambda s: s * s * float(R.trunc_mean_normal(3.0 / s)) - mean * m0, 0.1, 10)
        Y = np.minimum((s * rng.standard_normal((P, m))) ** 2, 9.0) / m0
    alarm = _run_default(Y, mean, A)
    rate = np.mean((alarm >= 0) & (alarm < m))
    du, dd = (d, 0.0) if direction == "up" else (0.0, d)
    bound = float(R.sup_bound(du, dd, m, A, B))
    assert bound < 1  # the check is not vacuous
    assert rate <= bound + 3 * np.sqrt(bound * (1 - bound) / P)
    assert rate > 0


def test_deterministic_stream_attains_the_bound():
    """Y = 1 + delta at every step: one side alone, one bet, level A alarms exactly when S_m(lam delta) >= A."""
    A, lam, d = 500.0, 0.5, 0.02
    m_star = int(np.argmax(R.sr_growth_sum(lam * d, np.arange(1, 5000)) >= A)) + 1
    up = betting_edetector(np.full(5000, 1 + d), 1.0, "up", threshold=A, grid=[lam])
    assert up.alarm[0] + 1 == m_star
    assert float(R.sr_growth_sum(lam * d, m_star - 1)) / A < 1 <= float(R.sr_growth_sum(lam * d, m_star)) / A
    # two-sided default: the bound at the alarm step is at most (1 + lam_max d)(1 + 1/(2A)) + m*/(2A)
    s = R.deterministic_alarm_step(1 + d, 7800, B)
    b = float(R.sup_bound(d, 0, s, 7800, B))
    assert 1 <= b <= (1 + 0.5 * d) * (1 + 1 / 15600) + s / 15600
    assert float(R.sup_bound(d, 0, s - 1, 7800, B)) < 1 + s / 15600


def test_run_paths_matches_package_detector():
    rng = np.random.default_rng(3)
    A = 40.0
    r = rng.standard_normal((50, 300)) * np.where(rng.random((50, 1)) < 0.5, 1.4, 0.8)
    tau2, m0 = 9.0, float(R.trunc_mean_normal(3.0))
    x = np.minimum(r * r, tau2)
    res = R.run_paths(x / m0, np.zeros_like(x), tau2 / m0, A, horizons=(300,))
    up = betting_edetector(x, m0, "up", threshold=2 * A)
    dn = betting_edetector(x, m0, "down", upper=tau2, threshold=2 * A)
    assert np.array_equal(res["alarm"], two_sided_alarm(up, dn))
    assert np.any(res["alarm"] >= 0)


def _random_drift_stream(rng, P, n, eps=0.08, stay=0.95):
    state = np.empty((P, n), bool)
    state[:, 0] = rng.random(P) < 0.5
    flip = rng.random((P, n)) > stay
    for t in range(1, n):
        state[:, t] = state[:, t - 1] ^ flip[:, t]
    mean = np.where(state, 1 + eps, 1 - eps)
    Y = np.where(rng.random((P, n)) < mean / B, B, 0.0)
    return Y, mean - 1


def test_lemma2_pathwise_inequality():
    rng = np.random.default_rng(11)
    Y, d = _random_drift_stream(rng, 400, 400)
    res = R.run_paths(Y, d, B, 50.0, horizons=(400,), check=True)
    assert np.all(res["ratio_max"] <= 1 + 1e-9)
    assert np.any(res["alarm"] >= 0)


def test_cusum_and_sr_sums_brute_force():
    rng = np.random.default_rng(5)
    d = rng.uniform(-0.3, 0.4, (3, 25))
    lam, phi = R.default_bets(B)
    cus = R.drift_cusum(d, lam, phi, [10, 25])
    sr = R.drift_sr_sums(d, lam, phi, [25])[25]
    for p in range(3):
        for m in (10, 25):
            best = -np.inf
            for t in range(m):
                for j in range(t + 1):
                    seg = d[p, j : t + 1]
                    best = max(best, max(np.prod(1 + lam_k * seg) for lam_k in lam),
                               max(np.prod(1 - ph * seg) for ph in phi))
            assert cus[m][p] == pytest.approx(best, rel=1e-12)
        for k, lam_k in enumerate(lam):
            brute = sum(np.prod(1 + lam_k * np.maximum(d[p, j:], 0)) for j in range(25))
            assert sr[0][p, k] == pytest.approx(brute, rel=1e-12)
        for l, ph in enumerate(phi):
            brute = sum(np.prod(1 + ph * np.maximum(-d[p, j:], 0)) for j in range(25))
            assert sr[1][p, l] == pytest.approx(brute, rel=1e-12)


def test_threshold_bound_brute_force():
    rng = np.random.default_rng(2)
    D = np.exp(rng.normal(0.2, 0.5, 300))
    b, K = R.threshold_bound(D, 390, 7800)
    grid = np.concatenate([D, np.linspace(0.01, 50, 20000)])
    brute = min(1.0, min(k * 390 / 7800 + np.mean(D > k) for k in grid))
    assert b == pytest.approx(brute, abs=1e-12)


def test_pathwise_bounds_hold_on_switching_stream():
    """Theorems 2 and 3 and Proposition 4 on an exogenous switching drift (small A)."""
    rng = np.random.default_rng(13)
    A, m, P = 500.0, 150, 6000
    Y, d = _random_drift_stream(rng, P, m, eps=0.08, stay=0.5)
    res = R.run_paths(Y, d, B, A, horizons=(m,))
    rate = np.mean((res["alarm"] >= 0) & (res["alarm"] < m))
    h = res["horizons"][m]
    se = 3 * np.sqrt(rate * (1 - rate) / P) + 3 / P
    thr, _ = R.threshold_bound(h["dstar"], m, A)
    exo = float(np.mean(np.minimum(1, np.minimum(h["exo"], h["dstar"] * m / A))))
    md = R.mean_drift_bound(h["sum_up"], h["sum_dn"], m, A, B=B)
    trs = R.truncated_sup_bound(h["max_up"], h["max_dn"], m, A, B)[0]
    for b in (thr, exo, md, trs):
        assert rate <= b + se
    assert exo < 1 and thr < 1
    # the drift-corrected detector is an e-detector: at most m / A
    rc = np.mean((res["alarm_corrected"] >= 0) & (res["alarm_corrected"] < m))
    assert rc <= m / A + 3 * np.sqrt(m / A / P)


def test_two_step_counterexample():
    lam, c = 0.5, 2.0
    el, ed = R.two_step_counterexample(lam, c)
    assert el == pytest.approx(1.75)
    assert ed == pytest.approx(1.5)
    assert el - ed == pytest.approx(lam**2 * c / 2)
    rng = np.random.default_rng(0)
    y1 = np.where(rng.random(400000) < 0.5, 2.0, 0.0)
    y2 = 1 + c * (y1 == 2)
    mc = np.mean((1 + lam * (y1 - 1)) * (1 + lam * (y2 - 1)))
    assert mc == pytest.approx(el, abs=0.01)


def test_jackpot_counterexample():
    """Drift B - 1 for s* steps on a fraction p of paths: P(alarm) = p, mean drift p (B - 1) s* / m."""
    A = 7800.0
    s = R.deterministic_alarm_step(B, A, B)
    assert s == 7
    p, m = 0.1, 390
    rng = np.random.default_rng(4)
    P = 4000
    jack = rng.random(P) < p
    Y = np.ones((P, m))
    Y[jack, :s] = B
    d = Y - 1.0
    res = R.run_paths(Y, d, B, A, horizons=(m,))
    rate = np.mean((res["alarm"] >= 0) & (res["alarm"] < m))
    assert rate == pytest.approx(jack.mean())
    mean_drift = np.mean(np.abs(d))
    assert mean_drift == pytest.approx(jack.mean() * (B - 1) * s / m)
    h = res["horizons"][m]
    assert R.mean_drift_bound(h["sum_up"], h["sum_dn"], m, A, B=B) >= rate


def test_plugin_bound_dominates_deterministic_stream():
    """Corollary 1 per window: the capped bound is at least the alarm indicator of Y = 1 + d(W) (E16's q = 0)."""
    rng = np.random.default_rng(8)
    n, W, m = 390, 300, 390
    r = rng.standard_normal((W, n))
    scale = np.median(np.abs(r), axis=1) * MAD_NORMAL_SCALE
    tau2 = (3 * scale) ** 2
    m_hat = np.mean(np.minimum(r * r, tau2[:, None]), axis=1)
    m_true = R.trunc_mean_normal(np.sqrt(tau2))
    d = R.plugin_drift(m_true, m_hat)
    Bh = tau2 / m_hat
    b = R.plugin_bound(d, Bh, m, 7800.0, per_window=True)
    hit = 0
    for i in range(W):
        y = np.full((1, m), 1 + d[i])
        a = R.run_paths(y, y - 1, Bh[i], 7800.0, horizons=(m,))["alarm"][0]
        if 0 <= a < m:
            hit += 1
            assert b[i] >= 1 - 1e-12
    assert hit > 0
    assert R.plugin_bound(d, Bh, m, 7800.0) == pytest.approx(b.mean())
    assert np.all(Bh < 2 * (3 * MAD_NORMAL_SCALE) ** 2 + 1e-9)  # B-hat <= 39.6 (section 4)


def test_default_bets_shape():
    lam, phi = R.default_bets(np.array([5.0, 9.0]))
    assert np.allclose(lam, DEFAULT_UP_GRID)
    assert phi.shape == (2, 5)
    assert phi[0, -1] == pytest.approx(0.95 / 4)
    with pytest.raises(ValueError):
        R.default_bets(1.0)
