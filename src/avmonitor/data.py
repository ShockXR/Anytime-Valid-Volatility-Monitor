"""Loaders for 1-minute bars and intraday deseasonalisation.

Supported sources (all public, fetched over raw.githubusercontent.com because
exchange APIs are blocked in the cloud environment this project was built in):

* Bitstamp BTC/USD, from https://github.com/ff137/bitstamp-btcusd-minute-data
  (columns: timestamp, open, high, low, close, volume; timestamp is the UTC open
  time in seconds; missing minutes are filled with flat zero-volume candles).
* Binance spot 1-minute klines for BTCUSDT, ETHUSDT, XRPUSDT, SOLUSDT and
  DOGEUSDT, from the per-symbol archives at
  https://github.com/parshant231234-eng/<SYMBOL>-Historical-Data (one CSV per
  month, ``1m/<YYYY>/<MM>/<SYMBOL>_1m_<YYYY-MM>.csv``; columns timestamp in
  milliseconds, datetime_utc, open, high, low, close, volume, quote_volume,
  trade_count, taker_buy_volume, taker_buy_quote_volume). The archive is a
  third-party mirror of Binance public kline data, so it is checked on load
  (duplicates, gaps, non-positive prices) and every check is reported.
* Nifty 50 index 1-minute bars (IST, 09:15 to 15:29), from
  https://github.com/technovusin/nifty50-historical-data (one CSV per year for
  2017 to 2025, one per month from 2026; columns Epoch, Timestamp, Open, High,
  Low, Close and, from 2026, Volume).

Raw files are not committed; ``scripts/download_data.py`` fetches them and
writes a manifest with row counts and SHA-256 hashes (``data/MANIFEST.csv``).
"""

from __future__ import annotations

import hashlib
import urllib.request
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


# ----------------------------------------------------------------------------
# Downloading with a local cache
# ----------------------------------------------------------------------------

DATA_ROOT = Path(__file__).resolve().parents[2] / "data" / "raw"

BINANCE_REPO = "https://raw.githubusercontent.com/parshant231234-eng/{symbol}-Historical-Data/main"
NIFTY_REPO = "https://raw.githubusercontent.com/technovusin/nifty50-historical-data/main"
CRYPTO_SYMBOLS = ("BTCUSDT", "ETHUSDT", "XRPUSDT", "SOLUSDT", "DOGEUSDT")


def fetch(url: str, dest: Path, overwrite: bool = False, timeout: float = 120.0) -> Path:
    """Download ``url`` to ``dest`` unless it is already there. Returns ``dest``."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0 and not overwrite:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=timeout) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    tmp.replace(dest)
    return dest


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _months(start: str, end: str) -> list[pd.Period]:
    """Calendar months touched by ``[start, end)``."""
    s = pd.Period(pd.Timestamp(start), "M")
    e = pd.Period(pd.Timestamp(end) - pd.Timedelta(minutes=1), "M")
    return list(pd.period_range(s, e, freq="M"))


def binance_month_path(symbol: str, month: pd.Period, root: Path = DATA_ROOT) -> Path:
    return Path(root) / "crypto" / "binance" / symbol / f"{symbol}_1m_{month.year}-{month.month:02d}.csv"


def binance_month_url(symbol: str, month: pd.Period) -> str:
    base = BINANCE_REPO.format(symbol=symbol)
    return f"{base}/1m/{month.year}/{month.month:02d}/{symbol}_1m_{month.year}-{month.month:02d}.csv"


def download_binance(symbol: str, start: str, end: str, root: Path = DATA_ROOT) -> list[Path]:
    return [fetch(binance_month_url(symbol, m), binance_month_path(symbol, m, root)) for m in _months(start, end)]


def nifty_file_path(key: str, root: Path = DATA_ROOT) -> Path:
    return Path(root) / "nifty" / f"NIFTY50_1min_{key}.csv"


def nifty_file_url(key: str) -> str:
    year = key[:4]
    return f"{NIFTY_REPO}/1min/{year}/NIFTY50_1min_{key}.csv"


def _nifty_keys(start: str, end: str) -> list[str]:
    keys = []
    for m in _months(start, end):
        if m.year <= 2025:
            k = str(m.year)
        else:
            k = f"{m.year}-{m.month:02d}"
        if k not in keys:
            keys.append(k)
    return keys


def download_nifty(start: str, end: str, root: Path = DATA_ROOT) -> list[Path]:
    return [fetch(nifty_file_url(k), nifty_file_path(k, root)) for k in _nifty_keys(start, end)]


# ----------------------------------------------------------------------------
# Loading and checking
# ----------------------------------------------------------------------------


def _ts(x, tz: str) -> pd.Timestamp:
    t = pd.Timestamp(x)
    return t.tz_localize(tz) if t.tzinfo is None else t.tz_convert(tz)


def _finish(df: pd.DataFrame, start, end, tz: str) -> pd.DataFrame:
    df = df.sort_index(kind="mergesort")  # stable, so keep="first" keeps the first row read
    df = df[~df.index.duplicated(keep="first")]
    if start is not None:
        df = df[df.index >= _ts(start, tz)]
    if end is not None:
        df = df[df.index < _ts(end, tz)]
    return df


def load_binance(symbol: str, start: str, end: str, root: Path = DATA_ROOT, download: bool = True) -> pd.DataFrame:
    """Binance spot 1-minute bars for ``symbol`` on ``[start, end)`` (UTC).

    The frame is indexed by the UTC open time of each minute and has ``open``,
    ``high``, ``low``, ``close``, ``volume`` and ``trade_count``. Missing minutes
    are inserted as NaN rows (never forward filled), so ``logret`` is NaN on both
    sides of a gap and gaps cannot create artificial zero returns. The checks
    are stored in ``df.attrs["checks"]``.
    """
    paths = download_binance(symbol, start, end, root) if download else [
        binance_month_path(symbol, m, root) for m in _months(start, end)
    ]
    parts = []
    for p in paths:
        d = pd.read_csv(p, usecols=["timestamp", "open", "high", "low", "close", "volume", "trade_count"])
        parts.append(d)
    raw = pd.concat(parts, ignore_index=True)
    n_raw = len(raw)
    raw["time"] = pd.to_datetime(raw["timestamp"], unit="ms", utc=True)
    df = raw.drop(columns="timestamp").set_index("time")
    n_dup = int(df.index.duplicated().sum())
    df = _finish(df, start, end, "UTC")
    full = pd.date_range(pd.Timestamp(start, tz="UTC"), pd.Timestamp(end, tz="UTC"), freq="1min", inclusive="left")
    missing = int(len(full.difference(df.index)))
    df = df.reindex(full)
    df.index.name = "time"
    bad_price = int((df["close"] <= 0).sum())
    df.loc[df["close"] <= 0, ["open", "high", "low", "close"]] = np.nan
    df["logret"] = np.log(df["close"]).diff()
    df["zero_volume"] = df["volume"].fillna(0) <= 0
    df.attrs["checks"] = dict(source=f"Binance spot {symbol} via github.com/parshant231234-eng mirror",
                              rows_read=n_raw, duplicates=n_dup, missing_minutes=missing,
                              nonpositive_prices=bad_price, first=str(df.index[0]), last=str(df.index[-1]))
    return df


def load_nifty(start: str, end: str, root: Path = DATA_ROOT, download: bool = True) -> pd.DataFrame:
    """Nifty 50 index 1-minute bars on ``[start, end)`` (dates in IST).

    Returns a frame indexed by IST minute with ``open``, ``high``, ``low``,
    ``close``, ``session`` (the trading date) and ``logret``. Returns are
    computed within each session only: the first bar of every day has NaN
    ``logret``, so overnight gaps never enter the intraday statistics. The
    index is the bar start time (09:15 is the first bar of a regular session).
    """
    keys = _nifty_keys(start, end)
    paths = download_nifty(start, end, root) if download else [nifty_file_path(k, root) for k in keys]
    parts = [pd.read_csv(p) for p in paths]
    raw = pd.concat(parts, ignore_index=True)
    raw.columns = [c.strip().lower() for c in raw.columns]
    raw["time"] = pd.to_datetime(raw["epoch"], unit="s", utc=True).dt.tz_convert("Asia/Kolkata")
    df = raw.set_index("time")[["open", "high", "low", "close"]].astype(float)
    n_dup = int(df.index.duplicated().sum())
    df = _finish(df, start, end, "Asia/Kolkata")
    df["session"] = df.index.floor("D")
    lc = np.log(df["close"])
    df["logret"] = lc.groupby(df["session"]).diff()
    # A missing minute inside a session makes the next return span two minutes.
    gap = df.index.to_series().groupby(df["session"]).diff() != pd.Timedelta(minutes=1)
    df.loc[gap.values & df["logret"].notna().values, "logret"] = np.nan
    bars = df.groupby("session").size()
    df.attrs["checks"] = dict(source="Nifty 50 index via github.com/technovusin/nifty50-historical-data",
                              rows=len(df), duplicates=n_dup, sessions=int(bars.size),
                              median_bars_per_session=float(bars.median()),
                              short_sessions=int((bars < 360).sum()),
                              first=str(df.index[0]), last=str(df.index[-1]))
    return df


def regular_sessions(df: pd.DataFrame, min_bars: int = 370) -> pd.DataFrame:
    """Drop Nifty sessions that are short (special or muhurat sessions, outages)."""
    bars = df.groupby("session")["close"].transform("size")
    return df[bars >= min_bars]
