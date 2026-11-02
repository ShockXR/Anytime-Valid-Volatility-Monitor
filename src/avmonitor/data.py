"""Loaders for 1-minute bars and intraday deseasonalisation.

The first supported source is the public Bitstamp BTC/USD 1-minute file
maintained at https://github.com/ff137/bitstamp-btcusd-minute-data
(columns: timestamp, open, high, low, close, volume; timestamp is the UTC open
time in seconds; missing minutes are filled with flat zero-volume candles).
Raw files are not committed; see ``data/README.md``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .robust import MAD_NORMAL_SCALE

BITSTAMP_URL = (
    "https://raw.githubusercontent.com/ff137/bitstamp-btcusd-minute-data/main/"
    "data/updates/btcusd_bitstamp_1min_latest.csv"
)


def load_bitstamp_csv(path: str | Path, start: str | None = None, end: str | None = None) -> pd.DataFrame:
    """Read the Bitstamp file and return a frame indexed by UTC minute.

    ``start`` is inclusive and ``end`` exclusive (any string pandas parses).
    Adds ``logret`` (close-to-close log return, NaN for the first row) and
    ``zero_volume`` (True for filled or tradeless minutes).
    """
    df = pd.read_csv(path)
    df["time"] = pd.to_datetime(df["timestamp"], unit="s", utc=True)
    df = df.set_index("time").sort_index()
    df = df[~df.index.duplicated()]
    if start is not None:
        df = df[df.index >= pd.Timestamp(start, tz="UTC")]
    if end is not None:
        df = df[df.index < pd.Timestamp(end, tz="UTC")]
    df["logret"] = np.log(df["close"]).diff()
    df["zero_volume"] = df["volume"] <= 0
    return df


def intraday_scale(returns: pd.Series, lookback_days: int = 20, smooth: int = 15) -> pd.Series:
    """Predictable intraday volatility profile, one value per minute.

    For each minute of the day the scale is the scaled median of ``|r|`` at that
    minute over the previous ``lookback_days`` days (the current day is excluded),
    smoothed with a centred circular moving average of ``smooth`` minutes across
    the clock, and normalised so that its median over the day is 1. Medians keep
    jumps from distorting the profile (compare Boudt, Croux and Laurent 2011).
    The first ``lookback_days`` days get the overall profile of those days.
    """
    r = returns.abs()
    day = r.index.floor("D")
    minute = r.index.hour * 60 + r.index.minute
    wide = pd.DataFrame({"a": r.values, "day": day, "minute": minute}).pivot_table(
        index="day", columns="minute", values="a", aggfunc="first"
    )
    prof = wide.rolling(lookback_days, min_periods=max(5, lookback_days // 2)).median().shift(1)
    prof = prof.fillna(wide.iloc[:lookback_days].median())
    if smooth > 1:
        arr = prof.to_numpy()
        pad = smooth // 2
        ext = np.concatenate([arr[:, -pad:], arr, arr[:, :pad]], axis=1)
        ker = np.ones(smooth) / smooth
        arr = np.apply_along_axis(lambda v: np.convolve(v, ker, mode="valid"), 1, ext)
        prof = pd.DataFrame(arr, index=prof.index, columns=prof.columns)
    prof = prof.div(prof.median(axis=1), axis=0)
    long = prof.stack()
    idx = pd.MultiIndex.from_arrays([day, minute])
    out = long.reindex(idx).to_numpy()
    return pd.Series(out, index=returns.index, name="intraday_scale")


def daily_robust_scale(returns: pd.Series) -> pd.Series:
    """Scaled ``median(|r|)`` per UTC day (a robust per-minute standard deviation)."""
    return returns.abs().groupby(returns.index.floor("D")).median() * MAD_NORMAL_SCALE
