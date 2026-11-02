"""Offline tests for the data loaders (no network: small files written to tmp_path)."""

import numpy as np
import pandas as pd

from avmonitor import data as D


def _write_binance(root, symbol, month, minutes, drop=(), dup=()):
    t0 = pd.Timestamp(f"{month}-01", tz="UTC")
    ts = [t0 + pd.Timedelta(minutes=i) for i in range(minutes) if i not in drop]
    ts += [t0 + pd.Timedelta(minutes=i) for i in dup]
    close = 100 * np.exp(np.cumsum(np.full(len(ts), 1e-4)))
    df = pd.DataFrame({"timestamp": [int(t.value // 10**6) for t in ts], "datetime_utc": ts, "open": close,
                       "high": close, "low": close, "close": close, "volume": 1.0, "quote_volume": 1.0,
                       "trade_count": 5, "taker_buy_volume": 0.5, "taker_buy_quote_volume": 0.5})
    p = D.binance_month_path(symbol, pd.Period(month, "M"), root)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False)


def test_binance_gaps_become_nan_not_zero(tmp_path):
    _write_binance(tmp_path, "ETHUSDT", "2025-07", 60, drop=(10, 11), dup=(5,))
    df = D.load_binance("ETHUSDT", "2025-07-01 00:00", "2025-07-01 01:00", root=tmp_path, download=False)
    assert len(df) == 60
    assert df.attrs["checks"]["missing_minutes"] == 2
    assert df.attrs["checks"]["duplicates"] == 1
    # returns next to the gap are NaN, never an artificial zero or a two-minute return
    assert np.isnan(df.logret.iloc[10]) and np.isnan(df.logret.iloc[11]) and np.isnan(df.logret.iloc[12])
    assert np.allclose(df.logret.dropna(), 1e-4)


def test_nifty_returns_stay_within_sessions(tmp_path):
    rows = []
    for day in ["2024-01-01", "2024-01-02"]:
        t0 = pd.Timestamp(f"{day} 09:15", tz="Asia/Kolkata")
        for i in range(375):
            if day == "2024-01-02" and i == 100:
                continue  # one missing minute
            t = t0 + pd.Timedelta(minutes=i)
            px = 20000 + i + (500 if day == "2024-01-02" else 0)
            rows.append(dict(Epoch=int(t.timestamp()), Timestamp=t.strftime("%Y-%m-%dT%H:%M:%S"),
                             Open=px, High=px, Low=px, Close=px))
    p = D.nifty_file_path("2024", tmp_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(p, index=False)
    df = D.load_nifty("2024-01-01", "2024-01-03", root=tmp_path, download=False)
    assert df.attrs["checks"]["sessions"] == 2
    first = df.groupby("session").head(1)
    assert first.logret.isna().all()  # no overnight return
    assert df.logret.max() < 1e-3  # the 500 point overnight gap never appears
    after_gap = df.loc[pd.Timestamp("2024-01-02 09:15", tz="Asia/Kolkata") + pd.Timedelta(minutes=101), "logret"]
    assert np.isnan(after_gap)  # a two-minute return is dropped


def test_month_and_file_keys():
    assert [str(m) for m in D._months("2025-11-15", "2026-02-01")] == ["2025-11", "2025-12", "2026-01"]
    assert D._nifty_keys("2025-06-01", "2026-03-01") == ["2025", "2026-01", "2026-02"]
