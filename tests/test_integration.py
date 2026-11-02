"""End-to-end tests: data file -> loader -> deseasonalisation -> monitor -> alarms, and the CLI.

These run offline on small synthetic files written to tmp_path, so they check that the
pieces fit together, not the statistics (the unit tests and experiments do that).
"""

import subprocess
import sys

import numpy as np
import pandas as pd

from avmonitor import data as D
from avmonitor.cli import main as cli_main
from avmonitor.monitor import VolatilityMonitor
from avmonitor.simulate import change_returns


def _binance_month_with_change(root, symbol="ETHUSDT", month="2025-07", ratio=2.0, seed=0):
    """A month of synthetic minute bars whose volatility doubles on day 20 (with a U-shaped day)."""
    rng = np.random.default_rng(seed)
    t0 = pd.Timestamp(f"{month}-01", tz="UTC")
    n = 31 * 1440
    minute = np.arange(n) % 1440
    season = 1 + 0.5 * np.cos(2 * np.pi * minute / 1440)
    sd = 5e-4 * season * np.where(np.arange(n) >= 19 * 1440, ratio, 1.0)
    r = sd * rng.standard_normal(n)
    close = 2000 * np.exp(np.cumsum(r))
    ts = t0 + pd.to_timedelta(np.arange(n), unit="min")
    ms = (ts - pd.Timestamp(0, tz="UTC")) // pd.Timedelta(milliseconds=1)
    df = pd.DataFrame({"timestamp": ms, "datetime_utc": ts, "open": close, "high": close, "low": close,
                       "close": close, "volume": 1.0, "quote_volume": 1.0, "trade_count": 3,
                       "taker_buy_volume": 0.5, "taker_buy_quote_volume": 0.5})
    p = D.binance_month_path(symbol, pd.Period(month, "M"), root)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False)
    return t0 + pd.Timedelta(days=19)


def test_pipeline_detects_a_volatility_doubling(tmp_path):
    change_time = _binance_month_with_change(tmp_path)
    df = D.load_binance("ETHUSDT", "2025-07-01", "2025-08-01", root=tmp_path, download=False)
    s = D.intraday_scale(df.logret.fillna(0.0), lookback_days=5)
    u = (df.logret / s).to_numpy()
    ev = VolatilityMonitor(threshold=1440 / 0.05, cal_len=1440, null_band=0.25).run(np.nan_to_num(u, nan=0.0))
    ups = [df.index[e.index] for e in ev if e.direction == "up" and df.index[e.index] >= change_time]
    assert ups, "the doubling was not detected"
    assert ups[0] - change_time < pd.Timedelta(hours=6)


def test_monitor_is_quiet_on_a_stationary_stream():
    # 20 days of iid returns, A for a 5% budget per day: expect at most a few alarms in total
    r = change_returns(1, 20 * 1440, None, rng=np.random.default_rng(1)).r[0]
    ev = VolatilityMonitor(threshold=1440 / 0.05, cal_len=1440).run(r)
    assert len(ev) <= 4


def test_cli_simulate_runs(capsys):
    assert cli_main(["simulate", "--seed", "1"]) == 0
    out = capsys.readouterr().out
    assert "[SYNTHETIC]" in out and "robust e-detector" in out


def test_cli_run_on_local_files(tmp_path, monkeypatch, capsys):
    _binance_month_with_change(tmp_path)
    monkeypatch.setattr(D, "DATA_ROOT", tmp_path)
    orig = D.load_binance
    monkeypatch.setattr(D, "load_binance", lambda *a, **k: orig(*a, **{**k, "root": tmp_path}))
    assert cli_main(["run", "--source", "binance", "--symbol", "ETHUSDT", "--start", "2025-07-01",
                     "--end", "2025-08-01", "--lookback-days", "5"]) == 0
    out = capsys.readouterr().out
    assert "[REAL]" in out and "alarms in" in out


def test_module_entry_point():
    res = subprocess.run([sys.executable, "-m", "avmonitor", "--help"], capture_output=True, text=True, timeout=60)
    assert res.returncode == 0 and "simulate" in res.stdout
