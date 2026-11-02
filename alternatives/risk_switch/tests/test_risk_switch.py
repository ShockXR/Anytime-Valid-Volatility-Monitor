"""Tests for the risk switch prototype (synthetic data only, run in a few seconds)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import backtest as B  # noqa: E402
import strategies as S  # noqa: E402


def synthetic_frame(n_days=12, sd=5e-4, change_day=None, ratio=1.0, seed=0, spec=S.CRYPTO):
    """Gaussian minute returns with an optional change of standard deviation."""
    rng = np.random.default_rng(seed)
    n = n_days * spec.bars_per_session
    scale = np.full(n, sd)
    if change_day is not None:
        scale[change_day * spec.bars_per_session :] *= ratio
    r = rng.standard_normal(n) * scale
    idx = pd.date_range("2024-01-01", periods=n, freq="1min", tz="UTC")
    close = 100.0 * np.exp(np.cumsum(r))
    df = pd.DataFrame({"close": close}, index=idx)
    df["logret"] = np.log(df["close"]).diff()
    return df


@pytest.fixture(scope="module")
def prepared():
    return S.prepare(synthetic_frame(), S.CRYPTO)


def test_backtest_accounting_matches_hand_calculation():
    P = S.prepare(synthetic_frame(n_days=3), S.CRYPTO)
    w = np.ones(len(P.index))
    bt = B.run_backtest(P, w, cost_bps=10.0)
    # exposure 1 from the first bar: returns equal the asset's, one entry trade of size 1
    assert np.allclose(bt.gross[1:], P.pnl[1:])
    assert bt.dw.sum() == pytest.approx(1.0)
    assert bt.net[0] == pytest.approx(-10.0 / 1e4)
    # a step change of exposure is charged once, at the bar where it happens
    w2 = w.copy()
    w2[1000:] = 0.5
    bt2 = B.run_backtest(P, w2, cost_bps=10.0)
    assert bt2.dw[1000] == pytest.approx(0.5) and bt2.dw[1001] == 0.0
    assert bt2.gross[1001] == pytest.approx(0.5 * P.pnl[1001])


def test_session_returns_and_drawdown():
    sess = np.array([0, 0, 1, 1])
    net = np.array([0.01, 0.01, -0.02, 0.0])
    d = B.session_returns(net, sess)
    assert d == pytest.approx([1.01**2 - 1, -0.02])
    nav = np.array([1.0, 1.2, 0.9, 1.0, 1.3, 1.17])
    assert B.max_drawdown(nav) == pytest.approx(0.9 / 1.2 - 1)


def test_no_lookahead_in_every_rule():
    """Perturbing the data after bar t0 must not change any exposure at or before t0.

    The first 20 days are skipped because ``avmonitor.data.intraday_scale`` gives
    them the overall profile of those days (a warm-up convention). The experiment
    windows start more than 20 days after the first bar for the same reason.
    """
    df = synthetic_frame(n_days=30)
    P = S.prepare(df, S.CRYPTO)
    warm = 21 * 1440
    t0 = 27 * 1440 + 17
    df2 = df.copy()
    df2.loc[df2.index[t0 + 1 :], "logret"] *= 4.0
    df2["close"] = 100.0 * np.exp(df2["logret"].fillna(0.0).cumsum())
    assert not np.allclose(df2["logret"].iloc[t0 + 1 :], df["logret"].iloc[t0 + 1 :])
    P2 = S.prepare(df2, S.CRYPTO)
    k = 1.0
    rules = {
        "ewma_minute": lambda Q: S.ewma_exposure(Q, 0.25, "minute", 0.3, 2.0),
        "ewma_hour": lambda Q: S.ewma_exposure(Q, 0.25, "hour", 0.3, 2.0),
        "ewma_day": lambda Q: S.ewma_exposure(Q, 0.25, "day", 0.3, 2.0),
        "robust_roll": lambda Q: S.robust_rolling_exposure(Q, 240, k, 0.3, 2.0),
        "switch": lambda Q: S.switch_exposure(
            Q, S.monitor_events(Q.u, S.make_monitor("robust", 200.0, 240, 0.0)), k, 0.3, 2.0, fast_window=60
        ),
    }
    for name, rule in rules.items():
        w1, w2 = rule(P), rule(P2)
        assert np.allclose(w1[warm : t0 + 1], w2[warm : t0 + 1]), name
        assert not np.allclose(w1[t0 + 1 :], w2[t0 + 1 :]), name  # the rule does react to the data
    # and the backtest only uses w_t from t + 1 on
    bt1 = B.run_backtest(P, rules["ewma_minute"](P), 5.0)
    bt2 = B.run_backtest(P2, rules["ewma_minute"](P2), 5.0)
    assert np.allclose(bt1.net[warm : t0 + 1], bt2.net[warm : t0 + 1])


def test_ewma_targeting_hits_target_on_constant_volatility():
    df = synthetic_frame(n_days=70, sd=5e-4, seed=1)
    P = S.prepare(df, S.CRYPTO)
    true_ann = 5e-4 * np.sqrt(525600)
    target = 0.5 * true_ann
    w = S.ewma_exposure(P, 2.0, "day", target, 5.0)
    win = np.arange(len(w)) >= 10 * 1440
    bt = B.run_backtest(P, w, 0.0)
    m = B.metrics(P, bt, win, target, 0.0, n_boot=200)
    assert abs(m["mean_exposure"] - 0.5) < 0.1
    # 60 daily returns: the sampling error of a standard deviation is about 9%
    assert abs(m["realised_vol"] / target - 1) < 0.25
    assert m["tracking_rms"] < 0.25
    assert m["sharpe_lo"] <= m["sharpe"] <= m["sharpe_hi"]


def test_update_frequency_limits_trades(prepared):
    P = prepared
    wm = S.ewma_exposure(P, 0.5, "minute", 0.3, 2.0)
    wh = S.ewma_exposure(P, 0.5, "hour", 0.3, 2.0)
    wd = S.ewma_exposure(P, 0.5, "day", 0.3, 2.0)
    changes = lambda w: int((np.abs(np.diff(w)) > 0).sum())  # noqa: E731
    assert changes(wd) < changes(wh) < changes(wm)
    assert changes(wd) <= 12 and changes(wh) <= 12 * 24
    # hourly exposure equals the minute exposure at the top of each hour
    assert np.allclose(wh[P.hour_bar][1:], wm[P.hour_bar][1:])


def test_switch_cuts_exposure_after_volatility_doubles():
    df = synthetic_frame(n_days=14, sd=5e-4, change_day=8, ratio=2.0, seed=2)
    P = S.prepare(df, S.CRYPTO)
    k = S.robust_to_sd_factor(P)
    assert 0.8 < k < 1.2  # Gaussian data: robust scale is already a standard deviation
    mon = S.make_monitor("robust", 1440.0, 720, 0.0)
    ev = S.monitor_events(P.u, mon)
    alarms = [e for e in ev if e["kind"] == "alarm"]
    assert any(e["direction"] == "up" and 8 * 1440 <= e["index"] < 9 * 1440 for e in alarms)
    w = S.switch_exposure(P, ev, k, 0.3, 3.0, fast_window=120)
    before = w[5 * 1440 : 8 * 1440].mean()
    after = w[10 * 1440 :].mean()
    assert 0.35 < after / before < 0.7
    # exposure is piecewise constant: it changes only at events
    n_changes = int((np.abs(np.diff(w)) > 0).sum())
    assert n_changes <= len(ev)


def test_block_bootstrap_indices_are_valid():
    rng = np.random.default_rng(0)
    idx = B.block_bootstrap_indices(37, 10, 50, rng)
    assert idx.shape == (50, 37) and idx.min() >= 0 and idx.max() < 37
    # consecutive entries inside a block step by one (circularly)
    assert np.all((idx[:, 1] - idx[:, 0]) % 37 == 1)
