"""Minute-level backtest, risk metrics and block-bootstrap intervals.

Accounting
----------
Exposure ``w_t`` is decided at the close of bar ``t`` and earns ``w_t * R_{t+1}``,
where ``R`` is the simple close-to-close return of the asset (overnight gaps
included). Changing the exposure from ``w_{t-1}`` to ``w_t`` costs
``cost_bps / 10,000 * |w_t - w_{t-1}|`` of wealth, charged at bar ``t``. The
one-way cost per unit of turnover is the only cost assumption: there is no
financing cost on leverage and no market impact beyond it. Wealth compounds bar
by bar: ``NAV_{t+1} = NAV_t (1 + w_t R_{t+1} - cost_t)``.

Metrics (all on the evaluation window only)
-------------------------------------------
* ``realised_vol``: annualised standard deviation of daily (per session) net returns.
* ``tracking_rms``: root mean square of ``log(RV_week / target)``, where ``RV_week``
  is the annualised realised volatility of net minute returns over a calendar
  week. It measures how far the delivered risk sits from the target. Lower is better.
* ``vol_of_vol``: standard deviation of ``log(RV_week)``, the same dispersion
  without a target (so buy and hold can be compared).
* ``max_drawdown``: largest peak-to-trough fall of wealth, on the minute path.
* ``worst_day``: the worst single session return.
* ``sharpe``: annualised mean over standard deviation of daily net returns, with a
  zero risk-free rate. The 95% interval comes from a circular block bootstrap of
  the daily returns (blocks of 10 sessions keep short-range dependence).
* ``turnover_per_year``: sum of ``|w_t - w_{t-1}|`` per year; ``trades_per_year``
  counts the bars at which exposure changed; ``cost_drag_per_year`` is
  ``cost * turnover``, the annual return lost to trading.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from strategies import MarketSpec, Prepared


@dataclass
class BacktestResult:
    net: np.ndarray  # net return earned at each bar
    gross: np.ndarray  # gross return earned at each bar
    dw: np.ndarray  # |change in exposure| at each bar
    w: np.ndarray  # exposure decided at each bar


def run_backtest(P: Prepared, w: np.ndarray, cost_bps: float) -> BacktestResult:
    w = np.asarray(w, float)
    w_prev = np.r_[0.0, w[:-1]]
    dw = np.abs(w - w_prev)
    gross = np.r_[0.0, w[:-1] * P.pnl[1:]]
    net = gross - cost_bps / 1e4 * dw
    return BacktestResult(net, gross, dw, w)


def _group_sums(x: np.ndarray, codes: np.ndarray):
    """Sums of ``x`` per code (codes need not start at 0); returns (sums, counts) for non-empty codes."""
    c = codes - codes.min()
    sums = np.bincount(c, weights=x)
    counts = np.bincount(c)
    keep = counts > 0
    return sums[keep], counts[keep]


def session_returns(net: np.ndarray, session: np.ndarray) -> np.ndarray:
    """Compound minute returns within each session code."""
    s, _ = _group_sums(np.log1p(np.asarray(net, float)), np.asarray(session))
    return np.expm1(s)


def max_drawdown(nav: np.ndarray) -> float:
    peak = np.maximum.accumulate(nav)
    return float(np.min(nav / peak - 1.0))


def weekly_realised_vol(net: np.ndarray, week: np.ndarray, session: np.ndarray, spec: MarketSpec,
                        min_sessions: int = 3) -> np.ndarray:
    """Annualised realised volatility of net minute returns per calendar week (weeks with few sessions dropped)."""
    sum_r2, _ = _group_sums(np.asarray(net, float) ** 2, week)
    # number of distinct sessions per week
    pair = np.unique(np.stack([week, session], axis=1), axis=0)
    n_sess = np.bincount(pair[:, 0] - week.min())
    n_sess = n_sess[n_sess > 0]
    ok = n_sess >= min_sessions
    return np.sqrt(sum_r2[ok] / n_sess[ok] * spec.sessions_per_year)


def block_bootstrap_indices(n: int, block: int, n_boot: int, rng: np.random.Generator) -> np.ndarray:
    """Circular block bootstrap: ``n_boot`` rows of ``n`` indices into a length-``n`` series."""
    n_blocks = int(np.ceil(n / block))
    starts = rng.integers(0, n, size=(n_boot, n_blocks))
    offs = np.arange(block)
    idx = (starts[:, :, None] + offs[None, None, :]).reshape(n_boot, -1)[:, :n]
    return idx % n


def sharpe(daily: np.ndarray, sessions_per_year: float) -> float:
    sd = daily.std(ddof=1)
    return float(daily.mean() / sd * np.sqrt(sessions_per_year)) if sd > 0 else float("nan")


def sharpe_ci(daily: np.ndarray, sessions_per_year: float, block: int = 10, n_boot: int = 2000, seed: int = 0,
              other: np.ndarray | None = None) -> dict:
    """Bootstrap 95% interval for the Sharpe ratio (and, with ``other``, for the paired difference)."""
    rng = np.random.default_rng(seed)
    idx = block_bootstrap_indices(len(daily), block, n_boot, rng)
    x = daily[idx]
    s = x.mean(axis=1) / x.std(axis=1, ddof=1) * np.sqrt(sessions_per_year)
    lo, hi = np.nanpercentile(s, [2.5, 97.5])
    out = dict(sharpe_lo=float(lo), sharpe_hi=float(hi))
    if other is not None:
        y = other[idx]
        d = s - y.mean(axis=1) / y.std(axis=1, ddof=1) * np.sqrt(sessions_per_year)
        dlo, dhi = np.nanpercentile(d, [2.5, 97.5])
        out.update(diff=float(sharpe(daily, sessions_per_year) - sharpe(other, sessions_per_year)),
                   diff_lo=float(dlo), diff_hi=float(dhi))
    return out


def daily_returns(P: Prepared, bt: BacktestResult, window: np.ndarray) -> np.ndarray:
    return session_returns(bt.net[window], P.session[window])


def metrics(P: Prepared, bt: BacktestResult, window: np.ndarray, target: float | None, cost_bps: float,
            block: int = 10, n_boot: int = 2000, seed: int = 0) -> dict:
    """All metrics on the bars selected by the boolean ``window`` (``n_boot=0`` skips the bootstrap)."""
    net, dw, w = bt.net[window], bt.dw[window], bt.w[window]
    sess, week = P.session[window], P.week[window]
    spec = P.spec
    daily = session_returns(net, sess)
    n_sess = len(daily)
    years = n_sess / spec.sessions_per_year
    nav = np.cumprod(1.0 + net)
    rv = weekly_realised_vol(net, week, sess, spec)
    out = dict(
        sessions=n_sess,
        ann_return=float(daily.mean() * spec.sessions_per_year),
        total_return=float(nav[-1] - 1.0),
        realised_vol=float(daily.std(ddof=1) * np.sqrt(spec.sessions_per_year)),
        tracking_rms=float(np.sqrt(np.mean(np.log(rv / target) ** 2))) if target else float("nan"),
        vol_of_vol=float(np.std(np.log(rv), ddof=1)),
        median_weekly_vol=float(np.median(rv)),
        max_drawdown=max_drawdown(nav),
        worst_day=float(daily.min()),
        worst_week_vol=float(rv.max()),
        sharpe=sharpe(daily, spec.sessions_per_year),
        mean_exposure=float(w.mean()),
        max_exposure=float(w.max()),
        turnover_per_year=float(dw.sum() / years),
        trades_per_year=float((dw > 1e-12).sum() / years),
        cost_drag_per_year=float(cost_bps / 1e4 * dw.sum() / years),
        cost_bps=cost_bps,
    )
    if n_boot > 0:
        out.update(sharpe_ci(daily, spec.sessions_per_year, block, n_boot, seed))
    return out
