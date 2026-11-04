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


def load_bitstamp_csv(path: str | Path, start: str | None = None, end: str | None = None,
                      filled: str = "missing") -> pd.DataFrame:
    """Read the Bitstamp file and return a frame indexed by UTC minute.

    ``start`` is inclusive and ``end`` exclusive (any string pandas parses).
    Adds ``logret`` (close-to-close log return, NaN for the first row) and
    ``zero_volume`` (True for filled or tradeless minutes).

    ``filled`` sets the returns around a filled minute. ``"missing"`` (default)
    makes the returns on both sides of it NaN. ``"zero"`` is the behaviour
    before review round 9 (m5): the filled minute repeats the previous close,
    so its return is 0 and the next return carries the whole move since the
    last trade. ``"zero"`` breaks the rule of ``docs/math.md`` 5.6 never to
    zero-fill a missing return; it exists only so that the committed results
    of E5, E6b, E7a, E8e and E19, which were computed with it, can be
    reproduced exactly.
    """
    if filled not in ("missing", "zero"):
        raise ValueError("filled must be 'missing' or 'zero'")
    df = pd.read_csv(path)
    n_raw = len(df)
    df["time"] = pd.to_datetime(df["timestamp"], unit="s", utc=True)
    df = df.set_index("time")
    n_dup = int(df.index.duplicated().sum())
    # stable sort, so a duplicated minute always keeps the first row read (the default sort is not stable,
    # so which duplicate survived used to be arbitrary)
    df = _finish(df, start, end, "UTC")
    bad_other = _mask_nonfinite(df)  # first, so that the bars masked for a close <= 0 are not counted twice
    bad_price = _mask_nonpositive(df)
    df["zero_volume"] = df["volume"] <= 0
    # A filled (tradeless) minute repeats the previous close: no price was observed, so the returns on both
    # sides of it are missing, not zero (docs/math.md 5.6, item 5), as for the gaps of load_binance.
    close = df["close"] if filled == "zero" else df["close"].where(~df["zero_volume"])
    df["logret"] = np.log(close).diff()
    df.attrs["checks"] = dict(source="Bitstamp BTC/USD via github.com/ff137 mirror", rows_read=n_raw,
                              duplicates=n_dup, nonpositive_prices=bad_price, nonfinite_prices=bad_other,
                              first=str(df.index[0]) if len(df) else "", last=str(df.index[-1]) if len(df) else "")
    return df


def _mask_nonpositive(df: pd.DataFrame) -> int:
    """Set open, high, low and close to NaN on bars whose close is zero or negative; return how many.

    A price of 0 or below is a feed error, and its log return would be -inf or NaN and could pass as a huge
    move. The returns on both sides of such a bar become NaN (missing), never zero.
    """
    bad = (df["close"] <= 0).to_numpy()
    cols = [c for c in ("open", "high", "low", "close") if c in df.columns]
    if bad.any():
        df.loc[bad, cols] = np.nan
    return int(bad.sum())


def _mask_nonfinite(df: pd.DataFrame) -> int:
    """Set open, high, low and close to NaN on bars whose close is infinite or missing; return how many.

    Added in round 4 next to ``_mask_nonpositive``: an infinite close (a parse or feed error) would give
    returns of +inf and -inf. A missing close already gives NaN returns; it is counted so that the checks show
    it. Call it before any NaN rows for missing minutes are inserted, so that those are not counted.
    """
    close = pd.to_numeric(df["close"], errors="coerce").to_numpy(float)
    bad = ~np.isfinite(close)
    cols = [c for c in ("open", "high", "low", "close") if c in df.columns]
    if bad.any():
        df.loc[bad, cols] = np.nan
    return int(bad.sum())


def intraday_scale(returns: pd.Series, lookback_days: int = 20, smooth: int = 15, start: str = "past",
                   circular: bool = True) -> pd.Series:
    """Predictable intraday volatility profile, one value per minute.

    For each minute of the day the scale is the scaled median of ``|r|`` at that
    minute over the previous ``lookback_days`` days (the current day is excluded),
    smoothed with a centred moving average of ``smooth`` minutes across the clock,
    and normalised so that its median over the day is 1. Medians keep jumps from
    distorting the profile (compare Boudt, Croux and Laurent 2011). Missing
    returns (NaN) are left out of the medians; pass returns with their NaN, do
    not fill them with zeros.

    ``start`` sets what happens before ``max(5, lookback_days // 2)`` earlier days
    exist. ``"past"`` (default) returns NaN there, so every value uses past data
    only. ``"pooled"`` fills those days with the median profile of the first
    ``lookback_days`` days, which uses up to ``lookback_days`` days of future data
    (this was the behaviour before November 2026 and is kept to reproduce old
    results). Load a warm-up period before the monitored period to avoid both.

    ``circular=True`` wraps the smoothing window from the end of the day to its
    start, which is right for 24-hour markets. Use ``circular=False`` for session
    markets (the edges are then averaged over the available minutes only).
    """
    if start not in ("past", "pooled"):
        raise ValueError("start must be 'past' or 'pooled'")
    r = returns.abs()
    day = r.index.floor("D")
    minute = r.index.hour * 60 + r.index.minute
    wide = pd.DataFrame({"a": r.values, "day": day, "minute": minute}).pivot_table(
        index="day", columns="minute", values="a", aggfunc="first", dropna=False
    )
    prof = wide.rolling(lookback_days, min_periods=max(5, lookback_days // 2)).median().shift(1)
    if start == "pooled":
        prof = prof.fillna(wide.iloc[:lookback_days].median())
    if smooth > 1:
        arr = prof.to_numpy()
        pad = smooth // 2
        ker = np.ones(smooth)
        if circular:
            ext = np.concatenate([arr[:, -pad:], arr, arr[:, :pad]], axis=1)
            arr = np.apply_along_axis(lambda v: np.convolve(v, ker / smooth, mode="valid"), 1, ext)
        else:
            ok = np.isfinite(arr)
            num = np.apply_along_axis(lambda v: np.convolve(v, ker, mode="same"), 1, np.where(ok, arr, 0.0))
            den = np.apply_along_axis(lambda v: np.convolve(v, ker, mode="same"), 1, ok.astype(float))
            with np.errstate(invalid="ignore", divide="ignore"):
                arr = np.where(ok, num / den, np.nan)
        prof = pd.DataFrame(arr, index=prof.index, columns=prof.columns)
    prof = prof.div(prof.median(axis=1), axis=0)
    long = prof.stack(future_stack=True)
    idx = pd.MultiIndex.from_arrays([day, minute])
    out = long.reindex(idx).to_numpy()
    return pd.Series(out, index=returns.index, name="intraday_scale")


def mark_frozen(returns: pd.Series, min_run: int = 30) -> pd.Series:
    """Set runs of at least ``min_run`` consecutive exact-zero returns to NaN.

    A price that does not move for half an hour on a liquid market is almost
    always a frozen feed (for example the NSE trading halt of 24 February
    2021), not information about volatility. Left as zeros, such runs look like
    strong evidence for a volatility fall. Shorter runs of zeros (normal for
    coarse ticks) are kept. The rule looks at the whole run, so it is meant for
    cleaning historical data; a live feed should flag staleness at the source.
    """
    r = returns.copy()
    z = (r == 0).to_numpy()
    out = r.to_numpy(copy=True)
    n = z.size
    i = 0
    while i < n:
        if z[i]:
            j = i
            while j < n and z[j]:
                j += 1
            if j - i >= min_run:
                out[i:j] = np.nan
            i = j
        else:
            i += 1
    return pd.Series(out, index=r.index, name=r.name)


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
    bad_other = _mask_nonfinite(df)  # before the NaN rows of missing minutes are inserted
    df = df.reindex(full)
    df.index.name = "time"
    bad_price = _mask_nonpositive(df)
    df["logret"] = np.log(df["close"]).diff()
    df["zero_volume"] = df["volume"].fillna(0) <= 0
    df.attrs["checks"] = dict(source=f"Binance spot {symbol} via github.com/parshant231234-eng mirror",
                              rows_read=n_raw, duplicates=n_dup, missing_minutes=missing,
                              nonpositive_prices=bad_price, nonfinite_prices=bad_other,
                              first=str(df.index[0]), last=str(df.index[-1]))
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
    bad_other = _mask_nonfinite(df)  # first, so that the bars masked for a close <= 0 are not counted twice
    bad_price = _mask_nonpositive(df)
    df["session"] = df.index.floor("D")
    lc = np.log(df["close"])
    df["logret"] = lc.groupby(df["session"]).diff()
    # A missing minute inside a session makes the next return span two minutes.
    gap = df.index.to_series().groupby(df["session"]).diff() != pd.Timedelta(minutes=1)
    df.loc[gap.values & df["logret"].notna().values, "logret"] = np.nan
    bars = df.groupby("session").size()
    df.attrs["checks"] = dict(source="Nifty 50 index via github.com/technovusin/nifty50-historical-data",
                              rows=len(df), duplicates=n_dup, nonpositive_prices=bad_price,
                              nonfinite_prices=bad_other, sessions=int(bars.size),
                              median_bars_per_session=float(bars.median()),
                              short_sessions=int((bars < 360).sum()),
                              first=str(df.index[0]), last=str(df.index[-1]))
    return df


def regular_sessions(df: pd.DataFrame, min_bars: int = 370) -> pd.DataFrame:
    """Drop Nifty sessions that are short (special or muhurat sessions, outages)."""
    bars = df.groupby("session")["close"].transform("size")
    return df[bars >= min_bars]
