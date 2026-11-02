"""Tests for the panel monitor (src/avmonitor/panel.py). Synthetic only, no data files, a few seconds."""

import numpy as np
import pytest
from scipy import stats

from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from avmonitor.panel import (PanelEProcess, PanelMonitor, change_weights, ebh, ebh_count, fdp_summary, run_block)

TAU = 3.0


def kappa2(a):
    """E[min(Z^2, a^2)] for a standard normal Z (docs/math.md, section 3)."""
    return (2 * stats.norm.cdf(a) - 1) - 2 * a * stats.norm.pdf(a) + 2 * a * a * stats.norm.sf(a)


M0 = kappa2(TAU)  # exact null mean of min(r^2, 9) for N(0, 1) returns


def test_ebh_matches_definition_and_is_self_consistent():
    rng = np.random.default_rng(0)
    e = np.exp(rng.normal(1.0, 2.0, size=(500, 12)))
    alpha = 0.1
    mask = ebh(e, alpha)
    K = e.shape[1]
    for i in range(e.shape[0]):
        srt = np.sort(e[i])[::-1]
        ok = [r for r in range(1, K + 1) if srt[r - 1] >= K / (alpha * r)]
        R = max(ok) if ok else 0
        assert mask[i].sum() == R
        if R:
            assert set(np.argsort(-e[i])[:R]) == set(np.nonzero(mask[i])[0])
            assert np.all(e[i, mask[i]] >= K / (alpha * R))
    assert np.array_equal(ebh_count(e, alpha), mask.sum(1))
    assert not ebh(np.full((1, 5), 19.0), 0.05).any()  # below 1/alpha nothing is rejected
    assert ebh(np.array([[100.0, 1, 1, 1, 1]]), 0.05)[0, 0]  # K / alpha = 100 is enough on its own


def test_change_weights_sum_to_one():
    for kind, h in (("window", 50), ("pareto", 50), ("pareto", 1)):
        n = np.arange(1, 400)
        w, _ = change_weights(n, kind, h)
        _, tail = change_weights(399, kind, h)
        assert w.sum() + tail == pytest.approx(1.0, abs=1e-12)
        _, tail0 = change_weights(0, kind, h)
        assert tail0 == pytest.approx(1.0)


def test_window_eprocess_is_shiryaev_roberts_over_h_plus_cash():
    """With window weights, E_t = (M_up + M_down) / (2h) + (h - t) / h, M the e-detector of detectors.py."""
    rng = np.random.default_rng(1)
    n, h = 120, 120
    r = rng.normal(0, 1.2, n)
    ep = PanelEProcess(TAU**2, M0, M0, h, "window", reps=1, K=1)
    up = np.asarray(DEFAULT_UP_GRID)
    lam_dn = np.asarray(DEFAULT_DOWN_FRACTIONS) / (TAU**2 / M0 - 1)
    S_up, S_dn = np.zeros(up.size), np.zeros(lam_dn.size)
    for t in range(n):
        ep.step(r[t])
        x = min(r[t] ** 2, TAU**2)
        S_up = (1 + up * (x / M0 - 1)) * (S_up + 1)
        S_dn = (1 + lam_dn * (1 - x / M0)) * (S_dn + 1)
        expect = (S_up.mean() + S_dn.mean()) / (2 * h) + (h - t - 1) / h
        assert ep.live_value[0, 0] == pytest.approx(expect, rel=1e-10)


@pytest.mark.parametrize("weights", ["window", "pareto"])
def test_test_supermartingale_at_fixed_and_stopping_times(weights):
    """Under the exact null E[E_t] <= 1, also at a data-dependent stopping time (optional stopping)."""
    rng = np.random.default_rng(2)
    reps, n = 20000, 150
    r = rng.normal(size=(reps, n))
    ep = PanelEProcess(TAU**2, M0, M0, 50, weights, reps=reps, K=1)
    stopped = np.full(reps, np.nan)
    for t in range(n):
        ep.step(r[:, t : t + 1])
        v = ep.live_value[:, 0]
        hit = np.isnan(stopped) & (v >= 1.5)
        stopped[hit] = v[hit]
    final = ep.live_value[:, 0]
    stopped = np.where(np.isnan(stopped), final, stopped)
    for vals in (final, stopped):
        assert vals.mean() <= 1 + 4 * vals.std() / np.sqrt(reps)


def test_shiryaev_roberts_divided_by_time_is_not_an_e_value_at_a_stopping_time():
    """Theory note, Proposition 1: stop at 1 if f_1 >= 1, else at 2; then E[M_tau / tau] > 1."""
    rng = np.random.default_rng(3)
    lam = np.mean(DEFAULT_UP_GRID)
    y = np.minimum(rng.normal(size=(400000, 2)) ** 2, TAU**2) / M0
    f = 1 + lam * (y - 1)
    m1 = f[:, 0]
    m2 = f[:, 1] * (m1 + 1)
    val = np.where(m1 >= 1, m1, m2 / 2)
    margin = np.mean(np.maximum(1 - m1, 0)) / 2  # the exact excess E[(1 - f_1)^+] / 2
    assert val.mean() > 1 + margin / 2
    assert margin > 0.03


def test_rule_ordering_is_pathwise():
    """Naive flags no later than e-BH, and e-BH no later than Bonferroni, on every path."""
    rng = np.random.default_rng(4)
    reps, K, n = 200, 8, 300
    r = rng.normal(size=(reps, K, n))
    r[:, :4, 100:] *= 1.6
    out = {rule: run_block(r, TAU**2, M0, M0, alpha=0.05, rule=rule).flag_time for rule in ("naive", "ebh", "bonferroni")}
    inf = np.iinfo(np.int64).max
    t = {k: np.where(v < 0, inf, v) for k, v in out.items()}
    assert np.all(t["naive"] <= t["ebh"])
    assert np.all(t["ebh"] <= t["bonferroni"])
    assert (out["ebh"] >= 0).sum() > (out["bonferroni"] >= 0).sum()


def test_sticky_ebh_sets_grow_and_stay_self_consistent():
    rng = np.random.default_rng(5)
    K, n, alpha = 10, 400, 0.1
    r = rng.normal(size=(K, n))
    r[:6, 80:] *= 1.5
    res = run_block(r, TAU**2, M0, M0, alpha=alpha, rule="ebh", return_path=True)
    ft = res.flag_time[0]
    assert (ft >= 0).sum() >= 3
    for t in range(n):
        F = (ft >= 0) & (ft <= t)
        if F.any():
            assert np.all(res.path[t, F] >= K / (alpha * F.sum()) * (1 - 1e-12))
            assert F.sum() == ebh(res.path[t][None], alpha).sum()  # e-BH on stopped values returns this set


def test_global_null_family_wise_error():
    """No change anywhere: P(any discovery) <= alpha for e-BH and Bonferroni; the naive rule is far above."""
    rng = np.random.default_rng(6)
    reps, K, n, alpha = 1500, 10, 200, 0.2
    r = rng.normal(size=(reps, K, n))
    any_flag = {}
    for rule in ("ebh", "bonferroni", "naive"):
        any_flag[rule] = (run_block(r, TAU**2, M0, M0, alpha=alpha, rule=rule).flag_time >= 0).any(1).mean()
    se = np.sqrt(alpha * (1 - alpha) / reps)
    assert any_flag["ebh"] <= alpha + 3 * se
    assert any_flag["bonferroni"] <= alpha + 3 * se
    assert any_flag["naive"] > 2 * any_flag["ebh"]


def test_fdr_with_changes_and_dependence():
    """Common factor and common jumps, half the assets change: FDR at the end and at the first discovery <= alpha."""
    rng = np.random.default_rng(7)
    reps, K, n, alpha, rho = 800, 10, 200, 0.2, 0.7
    f = rng.normal(size=(reps, 1, n))
    r = np.sqrt(rho) * f + np.sqrt(1 - rho) * rng.normal(size=(reps, K, n))
    r[:, :5, 60:] *= 1.5
    nu = np.full((reps, K), n)
    nu[:, :5] = 60
    res = run_block(r, TAU**2, M0, M0, alpha=alpha, rule="ebh")
    s = fdp_summary(res.flag_time, nu, n)
    se = 0.5 / np.sqrt(reps)
    assert s["fdp_end"].mean() <= alpha + 3 * se
    assert np.nanmean(s["fdp_first"]) <= alpha + 3 * se
    assert (res.flag_time[:, :5] >= 60).mean() > 0.5


def test_fdp_summary_by_hand():
    ft = np.array([[3, 1, -1, 1]])
    nu = np.array([[0, 5, 0, 0]])  # asset 1 changes at 5, so its flag at 1 is false
    s = fdp_summary(ft, nu, 10)
    assert s["R"][0] == 3 and s["V"][0] == 1
    assert s["fdp_first"][0] == pytest.approx(0.5)  # step 1: assets 1 (false) and 3 (true)
    assert s["sup_fdp"][0] == pytest.approx(0.5)
    assert s["fdp_end"][0] == pytest.approx(1 / 3)


def test_missing_returns_leave_evalue_and_clock_unchanged():
    ep = PanelEProcess(TAU**2, M0, M0, 30, "pareto", reps=1, K=2)
    ep.step([2.0, 0.3])
    v, n = ep.live_value.copy(), ep.n.copy()
    ep.step([np.nan, 0.5])
    assert ep.live_value[0, 0] == v[0, 0] and ep.n[0, 0] == n[0, 0]
    assert ep.n[0, 1] == n[0, 1] + 1


def test_inactive_and_closed_assets_are_never_discovered():
    rng = np.random.default_rng(8)
    r = rng.normal(size=(1, 4, 300)) * 3.0  # every asset far above its baseline
    r[0, 3] = np.nan  # asset 3 never trades
    res = run_block(r, TAU**2, M0, M0, alpha=0.05, rule="ebh", active=[[True, True, False, True]])
    assert res.flag_time[0, 2] == -1 and res.flag_time[0, 3] == -1
    assert np.all(res.flag_time[0, :2] >= 0)
    assert res.final_value[0, 2] == 1.0 and res.final_value[0, 3] == 1.0


def test_streaming_monitor_blocks_and_consistency_with_run_block():
    rng = np.random.default_rng(9)
    K, n = 4, 1600
    R = rng.normal(size=(n, K))
    R[700:, :2] *= 2.0
    mon = PanelMonitor(["a", "b", "c", "d"], alpha=0.05, rule="ebh", block_len=300, cal_len=400, record=True)
    disc = mon.run(R)
    assert {d.name for d in disc} >= {"a", "b"}
    assert all(d.evalue >= d.threshold for d in disc)
    assert all(d.direction == "up" for d in disc if d.name in ("a", "b"))
    starts = [b.start for b in mon.blocks]
    assert starts[0] == 400 and np.all(np.diff(starts) == 300)
    # each block is the vectorised run_block on that block's returns with that block's baselines
    for b in mon.blocks:
        seg = R[b.start : b.start + 300].T
        res = run_block(seg, b.tau2[None], b.m_up[None], b.m_dn[None], alpha=0.05, rule="ebh", horizon=300,
                        active=b.active[None])
        got = {(d.asset, d.index - b.start) for d in disc if d.block_start == b.start}
        want = {(k, int(t)) for k, t in enumerate(res.flag_time[0]) if t >= 0}
        assert got == want
    assert mon.path.shape == (n, K)
    assert np.all(np.isnan(mon.path[:400]))


def test_streaming_monitor_with_session_asset_and_cs_baseline():
    """An asset that trades only part of the clock (NaN otherwise) gets its own clock and calibration."""
    rng = np.random.default_rng(10)
    n = 3000
    R = rng.normal(size=(n, 3))
    session = (np.arange(n) % 100) < 30  # asset 2 trades 30 of every 100 steps
    R[~session, 2] = np.nan
    mon = PanelMonitor(3, alpha=0.05, rule="ebh", block_len=500, horizon=[500, 500, 150], cal_len=[500, 500, 150],
                       baseline="cs", delta=0.01)
    mon.run(R)
    assert all(b.active.all() for b in mon.blocks)
    assert mon.blocks[0].n_cal[2] == 150
    with pytest.raises(ValueError):
        PanelMonitor(3, rule="holm", block_len=10)
