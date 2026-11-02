"""Exposure rules for the change-triggered risk switch and its benchmarks.

Every rule returns one exposure per one-minute bar, ``w_t``, computed from
information available at the close of bar ``t`` only. The backtest in
``backtest.py`` applies ``w_t`` to the return of bar ``t + 1``, so there is no
look-ahead by construction.

Rules
-----
* ``buy_and_hold``: ``w_t = 1``.
* ``ewma_exposure``: classic volatility targeting. An exponentially weighted
  moving average (EWMA) of squared minute returns gives a per-minute variance
  forecast, which is annualised and inverted: ``w_t = min(w_max, target / sigma_t)``.
  The exposure is refreshed every minute, every hour or once a day.
* ``robust_rolling_exposure``: the same inversion, but the forecast is the
  scaled median of absolute deseasonalised returns over a trailing window (a
  robust standard deviation), refreshed once a day. It separates "robust
  estimate" from "change-triggered timing".
* ``switch_exposure``: the change-triggered switch. Exposure changes only when
  the monitor (``avmonitor.monitor.VolatilityMonitor`` or the Gaussian
  likelihood-ratio benchmark) raises an alarm, and again when the monitor has
  recalibrated on the new regime. Between those events it is held fixed.

Units
-----
``target`` is an annualised volatility (0.40 means 40% a year). Per-minute
variances are annualised through a ``MarketSpec`` that says how many bars make a
session and how many sessions make a year. For markets that close overnight,
the close-to-open gap is treated as an extra return at the first bar of each
session, and its variance is forecast separately with an EWMA over sessions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from avmonitor.data import intraday_scale
from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE


@dataclass(frozen=True)
class MarketSpec:
    name: str
    bars_per_session: int  # one-minute bars in a trading session (1440 for 24/7 crypto, 375 for Nifty)
    sessions_per_year: float  # 365 for crypto, 252 for Nifty
    overnight: bool  # True if the first bar of a session carries a close-to-open gap


CRYPTO = MarketSpec("crypto", 1440, 365.0, False)
NIFTY = MarketSpec("nifty", 375, 252.0, True)

# Floor on the robust scale of deseasonalised returns (per-minute return units; the
# intraday profile is a unitless multiplier with median 1). It is a numerical guard
# only: 1e-6 per minute is 0.07% a year for a 24/7 market, far below any real regime.
# It bites only when more than half of a calibration window is exact zero returns
# (this happens on Binance BTC on the very quiet weekend of 22 to 23 June 2024),
# which would otherwise make the robust scale exactly 0 and break the monitor.
SCALE_FLOOR = 1e-6


@dataclass
class Prepared:
    """Per-bar arrays derived from a bar frame, all predictable (no look-ahead)."""

    index: pd.DatetimeIndex
    logret: np.ndarray  # log return of the bar, NaN at gaps and at the first bar of a session
    pnl: np.ndarray  # simple close-to-close return earned by holding the asset over the bar (includes overnight gaps)
    s: np.ndarray  # intraday volatility profile (median 1 over the day), from the previous 20 days
    u: np.ndarray  # deseasonalised return logret / s
    s2bar: np.ndarray  # trailing mean of s^2 over one session, turns a deseasonalised variance into a raw one
    overnight_var: np.ndarray  # EWMA forecast of the squared overnight gap (0 for 24/7 markets)
    session: np.ndarray  # integer session code of each bar (increasing)
    week: np.ndarray  # integer calendar-week code of each bar (increasing)
    first_bar: np.ndarray  # True at the first bar of a session
    last_bar: np.ndarray  # True at the last bar of a session
    hour_bar: np.ndarray  # True every 60th bar of a session (first bar included)
    spec: MarketSpec


def prepare(df: pd.DataFrame, spec: MarketSpec, profile_days: int = 20, overnight_halflife_sessions: float = 20.0) -> Prepared:
    """Build the predictable per-bar arrays used by every rule.

    ``df`` is a frame from ``avmonitor.data.load_binance`` or ``load_nifty``
    (after ``regular_sessions``), indexed by bar start time with ``close`` and
    ``logret``; Nifty frames also carry ``open`` and ``session``.
    """
    idx = df.index
    logret = df["logret"].to_numpy(float)
    close = df["close"].ffill()
    pnl = np.array(close.pct_change().to_numpy(float), copy=True)
    pnl[~np.isfinite(pnl)] = 0.0

    s = intraday_scale(df["logret"], lookback_days=profile_days)
    s = s.where(s > 0)
    u = np.array((df["logret"] / s).to_numpy(float), copy=True)
    u[~np.isfinite(u)] = np.nan
    s2bar = (s**2).rolling(spec.bars_per_session, min_periods=max(10, spec.bars_per_session // 4)).mean().to_numpy(float)

    sess = df["session"] if spec.overnight else pd.Series(idx.floor("D"), index=idx)
    codes = pd.factorize(sess)[0]
    naive = idx.tz_localize(None) if idx.tz is not None else idx
    week = pd.factorize(naive.to_period("W"))[0]
    first = np.r_[True, codes[1:] != codes[:-1]]
    last = np.r_[codes[1:] != codes[:-1], True]
    pos = pd.Series(np.ones(len(idx)), index=idx).groupby(codes).cumsum().to_numpy() - 1
    hour = (pos % 60) == 0

    if spec.overnight:
        lo = np.log(df["open"].to_numpy(float))
        lc = np.log(close.to_numpy(float))
        gap = np.full(len(idx), np.nan)
        gap[first] = lo[first] - np.r_[np.nan, lc[:-1]][first]
        g = pd.Series(gap**2)
        ov = g.ewm(halflife=overnight_halflife_sessions, adjust=True, ignore_na=True).mean().ffill().fillna(0.0).to_numpy()
    else:
        ov = np.zeros(len(idx))

    return Prepared(idx, logret, pnl, s.to_numpy(float), u, s2bar, ov, codes, week, first, last, hour, spec)


# ----------------------------------------------------------------------------
# From a variance forecast to an exposure
# ----------------------------------------------------------------------------


def annualised_vol(v_bar, overnight_var, spec: MarketSpec):
    """Annualised volatility implied by a per-bar variance (and a per-session overnight variance)."""
    per_session = np.asarray(v_bar, float) * spec.bars_per_session
    if spec.overnight:
        per_session = per_session + np.asarray(overnight_var, float)
    return np.sqrt(per_session * spec.sessions_per_year)


def exposure_from_vol(sigma_ann, target: float, w_max: float):
    """``min(w_max, target / sigma)``, 0 where the forecast is missing."""
    sig = np.asarray(sigma_ann, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        w = target / sig
    w = np.where(np.isfinite(w), w, 0.0)
    return np.clip(w, 0.0, w_max)


def hold_at(w: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Keep ``w`` only where ``mask`` is True and hold it until the next True (0 before the first)."""
    out = np.where(mask, w, np.nan)
    return pd.Series(out).ffill().fillna(0.0).to_numpy()


# ----------------------------------------------------------------------------
# Benchmarks
# ----------------------------------------------------------------------------


def buy_and_hold(P: Prepared) -> np.ndarray:
    return np.ones(len(P.index))


def ewma_var(logret: np.ndarray, halflife_bars: float) -> np.ndarray:
    """EWMA of squared returns; NaN returns are skipped and the last value is carried forward."""
    r2 = pd.Series(np.asarray(logret, float) ** 2)
    return r2.ewm(halflife=halflife_bars, adjust=True, ignore_na=True).mean().to_numpy(float)


def ewma_exposure(P: Prepared, halflife_days: float, update: str, target: float, w_max: float) -> np.ndarray:
    """Classic volatility targeting with an EWMA variance forecast.

    ``update`` is "minute" (refresh every bar), "hour" (every 60th bar of a
    session) or "day" (at the last bar of each session, so the new exposure
    applies from the next session's first bar).
    """
    v = ewma_var(P.logret, halflife_days * P.spec.bars_per_session)
    w = exposure_from_vol(annualised_vol(v, P.overnight_var, P.spec), target, w_max)
    if update == "minute":
        return w
    if update == "hour":
        return hold_at(w, P.hour_bar)
    if update == "day":
        return hold_at(w, P.last_bar)
    raise ValueError("update must be 'minute', 'hour' or 'day'")


def robust_rolling_exposure(P: Prepared, window_bars: int, robust_to_sd: float, target: float, w_max: float,
                            update: str = "day", scale_floor: float = SCALE_FLOOR) -> np.ndarray:
    """Volatility targeting with a trailing robust scale (scaled median of |u|) instead of an EWMA."""
    a = pd.Series(np.abs(P.u))
    m = a.rolling(window_bars, min_periods=window_bars // 2).median().to_numpy(float) * MAD_NORMAL_SCALE
    m = np.maximum(m, scale_floor)
    v = (m * robust_to_sd) ** 2 * P.s2bar
    w = exposure_from_vol(annualised_vol(v, P.overnight_var, P.spec), target, w_max)
    if update == "minute":
        return w
    return hold_at(w, P.hour_bar if update == "hour" else P.last_bar)


# ----------------------------------------------------------------------------
# The change-triggered switch
# ----------------------------------------------------------------------------


def make_monitor(kind: str, threshold: float, cal_len: int, null_band: float, c: float = 3.0,
                 scale_floor: float = SCALE_FLOOR):
    if kind == "robust":
        return VolatilityMonitor(threshold=threshold, c=c, cal_len=cal_len, null_band=null_band, scale_floor=scale_floor)
    if kind == "gaussian":
        return GaussianLRMonitor(threshold=threshold, c=c, cal_len=cal_len, null_band=null_band, scale_floor=scale_floor)
    raise ValueError("kind must be 'robust' or 'gaussian'")


def monitor_events(u: np.ndarray, monitor) -> list[dict]:
    """Run a monitor over the deseasonalised stream and record alarms and recalibrations.

    Each record has ``kind`` ("alarm" or "calibrated"), ``index`` (bar at whose
    close the event is known), ``direction`` ("up", "down" or "") and ``scale``
    (the robust standard deviation of the new calibration window, in
    deseasonalised units; NaN for alarms).
    """
    events = []
    was_calibrated = monitor.calibrated
    for t, x in enumerate(np.asarray(u, float)):
        ev = monitor.update(x)
        if ev is not None:
            events.append(dict(kind="alarm", index=t, direction=ev.direction, scale=np.nan, stat=ev.stat))
            was_calibrated = False
        elif not was_calibrated and monitor.calibrated:
            events.append(dict(kind="calibrated", index=t, direction="", scale=float(monitor.scale), stat=np.nan))
            was_calibrated = True
    return events


def switch_exposure(P: Prepared, events: list[dict], robust_to_sd: float, target: float, w_max: float,
                    fast_window: int = 120, provisional: bool = True, scale_floor: float = SCALE_FLOOR) -> np.ndarray:
    """Exposure that changes only at monitor events.

    At a recalibration the regime's robust per-minute standard deviation is the
    calibration scale times ``robust_to_sd`` (a fixed conversion from robust scale
    to standard deviation) times the intraday profile factor ``sqrt(s2bar)``. At an
    alarm, if ``provisional`` is True, a provisional exposure is set from the robust
    scale of the last ``fast_window`` deseasonalised returns, which mostly belong to
    the new regime because detection delays are typically shorter than that; it
    bridges the calibration period during which the monitor is collecting data.
    """
    n = len(P.index)
    w = np.full(n, np.nan)
    for e in events:
        t = e["index"]
        if e["kind"] == "alarm":
            if not provisional:
                continue
            seg = P.u[max(0, t - fast_window + 1) : t + 1]
            seg = seg[np.isfinite(seg)]
            if seg.size < fast_window // 2:
                continue
            sc = float(np.median(np.abs(seg))) * MAD_NORMAL_SCALE
        else:
            sc = e["scale"]
        sc = max(sc, scale_floor)
        v = (sc * robust_to_sd) ** 2 * P.s2bar[t]
        w[t] = exposure_from_vol(annualised_vol(v, P.overnight_var[t], P.spec), target, w_max)
    return pd.Series(w).ffill().fillna(0.0).to_numpy()


def robust_to_sd_factor(P: Prepared, mask: np.ndarray | None = None, min_bars: int | None = None) -> float:
    """Median over sessions of (realised sd of returns) / (robust sd implied by the profile).

    Converts the monitor's robust scale, which estimates the standard deviation of
    normal data, into a standard deviation for fat-tailed minute returns. Estimate
    it on the tuning sample only and freeze it.
    """
    n = len(P.index)
    use = np.ones(n, bool) if mask is None else np.asarray(mask, bool)
    min_bars = P.spec.bars_per_session * 9 // 10 if min_bars is None else min_bars
    ratios = []
    for code in np.unique(P.session[use]):
        sel = use & (P.session == code) & np.isfinite(P.u) & np.isfinite(P.logret)
        if sel.sum() < min_bars:
            continue
        sd = np.sqrt(np.mean(P.logret[sel] ** 2))
        rob = np.median(np.abs(P.u[sel])) * MAD_NORMAL_SCALE * np.sqrt(np.mean(P.s[sel] ** 2))
        if rob > 0:
            ratios.append(sd / rob)
    return float(np.median(ratios))
