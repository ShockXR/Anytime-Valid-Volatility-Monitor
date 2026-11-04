"""Tests for the conformal rank betting prototype (no data, no network, a few seconds)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import rankbet as R  # noqa: E402
import tsb_sim as S  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402


def _brute_pvalues(s, n_ref, U):
    reps, n = s.shape
    out = np.empty((reps, n - n_ref))
    for i in range(reps):
        for t in range(n - n_ref):
            x, prev = s[i, n_ref + t], s[i, : n_ref + t]
            out[i, t] = ((prev > x).sum() + U[i, t] * ((prev == x).sum() + 1)) / (n_ref + t + 1)
    return out


def test_pvalues_match_brute_force_and_streaming_with_ties():
    rng = np.random.default_rng(0)
    s = np.round(rng.standard_normal((4, 160)), 1)  # many ties
    p = R.pvalues_against_reference(s, 40, np.random.default_rng(3))
    U = np.random.default_rng(3).random((4, 120))
    assert np.allclose(p, _brute_pvalues(s, 40, U))
    rk = R.SequentialRanker(np.random.default_rng(9))
    for v in s[0, :40]:
        rk.add(v)
    ps = [rk.pvalue(v, u=U[0, t]) for t, v in enumerate(s[0, 40:])]
    assert np.allclose(ps, p[0])


def test_pvalues_are_uniform_and_uncorrelated_for_iid_discrete_data():
    rng = np.random.default_rng(1)
    s = rng.poisson(1.0, (1000, 300)).astype(float)  # heavy ties, like zero returns
    p = R.pvalues_against_reference(s, 50, rng)
    assert abs(p.mean() - 0.5) < 0.005
    assert abs(p.var() - 1 / 12) < 0.002
    assert abs(np.mean(p < 0.05) - 0.05) < 0.004
    assert abs(np.corrcoef(p[:, :-1].ravel(), p[:, 1:].ravel())[0, 1]) < 0.01


def test_gauss_score_has_mean_one_and_known_bounds():
    grid = (np.arange(200_000) + 0.5) / 200_000  # deterministic uniform grid
    for kind in ("abs", "signed"):
        y = R.gauss_score(grid, kind, 3.0)
        assert abs(y.mean() - 1) < 1e-3
        assert y.min() >= 0 and y.max() <= 9 / 0.9946 + 1e-9
    assert abs(R.gauss_score_abs_dev(3.0) - np.abs(R.gauss_score(grid) - 1).mean()) < 2e-3


@pytest.mark.parametrize("band", [0.0, 0.3])
def test_every_bet_is_nonnegative_and_fair_under_the_null(band):
    """For a fixed state of the learners: f >= 0 everywhere, E_U f = 1 when b = 0, and E_h f <= 1 for
    densities h with |h - 1| <= b (checked on h = 1 + b sign(w - median w), the most favourable one)."""
    cfg = R.RankBetConfig(band=band)
    core = R.RankBetCore(cfg, reps=1)
    rng = np.random.default_rng(2)
    for p in rng.beta(0.3, 0.3, 500):  # drive the learners away from uniform
        core.step(np.array([p]))
    grid = (np.arange(100_000) + 0.5) / 100_000
    y = R.gauss_score(grid, cfg.kind, cfg.a) - 1
    fams = {"up": (core.lam_up[:, None] * (y[None, :] - band * core.mY), y),
            "down": (core.lam_dn[:, None] * (-y[None, :] - band * core.mY), -y)}
    tot = core.counts.sum(axis=2, keepdims=True)
    pi = (core.counts + cfg.prior * core.q) / (tot + cfg.prior)
    k = np.clip(np.searchsorted(core.edges, grid, side="right") - 1, 0, core.K - 1)
    wk = (pi / core.q - 1.0)[0][:, k]
    mdev = np.abs(pi - core.q).sum(axis=2)[0][:, None]
    lam = (1 - cfg.eps) / (1 + band * mdev)
    fams["hist"] = (lam * (wk - band * mdev), wk)
    for name, (g, w) in fams.items():
        f = 1 + g
        assert f.min() > 0, name
        if band == 0:
            assert np.allclose(f.mean(axis=1), 1, atol=2e-3), name
        w2 = np.atleast_2d(w)
        h = 1 + band * np.sign(w2 - np.median(w2, axis=1, keepdims=True))
        h = h / h.mean(axis=1, keepdims=True)
        assert np.all((f * h).mean(axis=1) <= 1 + 2e-3), name


def test_shiryaev_roberts_mixture_has_mean_t_under_the_null():
    rng = np.random.default_rng(3)
    reps, m = 4000, 200
    core = R.RankBetCore(R.RankBetConfig(), reps)
    for _ in range(m):
        M = core.step(rng.random(reps))
    assert abs(M.mean() / m - 1) < 0.1  # E[M_t] = t exactly; Monte Carlo error is about 3%


def test_vectorised_main_monitor_matches_volatility_monitor():
    rng = np.random.default_rng(4)
    r = S.change_paths(6, 1440, 300, 1500, {"law": "jumps", "p": 0.003, "J": 6.0}, {"law": "gauss", "sd": 1.4}, rng)
    for band in (0.0, 0.25):
        alarm, direction = S.main_monitor_first_alarm(r, 2000.0, null_band=band)
        for i in range(r.shape[0]):
            mon = VolatilityMonitor(threshold=2000.0, cal_len=1440, null_band=band)
            first = None
            for t, x in enumerate(r[i]):
                ev = mon.update(x)
                if ev is not None:
                    first = ev
                    break
            assert first is not None and alarm[i] == first.index - 1440 and direction[i] == first.direction


def test_rank_monitor_detects_a_tail_change_at_fixed_variance():
    rng = np.random.default_rng(5)
    r = S.change_paths(1, 1440, 200, 3000, {"law": "gauss"}, {"law": "t", "df": 3}, rng)[0]
    mon = R.RankMonitor(threshold=2000.0, cal_len=1440, seed=1)
    ev = mon.run(r)
    assert ev and 1440 + 200 <= ev[0].index < 1440 + 200 + 1500  # first alarm after the change, and soon
    assert abs(ev[0].share_up + ev[0].share_down + ev[0].share_hist - 1) < 1e-9


def test_rank_monitor_skips_nan_and_validates_config():
    mon = R.RankMonitor(threshold=1e9, cal_len=10)
    mon.run([np.nan] * 5 + list(np.random.default_rng(6).standard_normal(30)))
    assert mon.calibrated and mon.n_ref == 10 and not mon.events
    with pytest.raises(ValueError):
        R.RankBetConfig(weights=(0.5, 0.5, 0.5))
    with pytest.raises(ValueError):
        R.RankBetConfig(kind="other")
    with pytest.raises(ValueError):
        R.RankBetConfig(edges=(0.0, 0.5, 0.4, 1.0))
