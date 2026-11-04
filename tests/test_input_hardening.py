"""Open items of the validation audit (reports/EXECUTION_SUMMARY.md, "Next steps" 1 to 4), fixed on 4 Nov 2026.

1. Loaders: a duplicated minute always keeps the first row read, and nonpositive prices are masked.
2. VolatilityMonitor: warns for a short calibration window and when it never finished calibrating.
3. Benchmark detectors (Gaussian LR, Page CUSUM) skip a missing value instead of going dead.
4. The command line prints a short error message instead of a traceback.

SYNTHETIC data only (small files written to tmp_path, no network).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from avmonitor import data as D
from avmonitor import detectors as Det
from avmonitor.cli import main as cli_main
from avmonitor.monitor import MIN_ADVISED_CAL_LEN, VolatilityMonitor


# ---------------------------------------------------------------------------------------------------------
# 1. loaders
# ---------------------------------------------------------------------------------------------------------


def _bitstamp(tmp_path, close, ts=None, vol=None):
    n = len(close)
    ts = 1_750_000_000 + 60 * np.arange(n) if ts is None else np.asarray(ts)
    vol = np.ones(n) if vol is None else vol
    p = tmp_path / "b.csv"
    pd.DataFrame(dict(timestamp=ts, open=close, high=close, low=close, close=close, volume=vol)).to_csv(p, index=False)
    return p


def test_bitstamp_duplicate_minute_keeps_the_first_row_read(tmp_path):
    base = 1_750_000_000 + 60 * np.arange(200)
    # the same minute appears twice, far apart in the file, with different prices; shuffle the rest so the
    # sort really has to move rows (an unstable sort could then keep either copy)
    rng = np.random.default_rng(3)
    order = rng.permutation(200)
    ts = np.r_[base[order], base[50]]
    close = np.r_[100.0 + order * 0.01, 999.0]
    df = D.load_bitstamp_csv(_bitstamp(tmp_path, close, ts=ts))
    t50 = pd.Timestamp(int(base[50]), unit="s", tz="UTC")
    assert len(df) == 200 and df.index.is_monotonic_increasing
    assert df.loc[t50, "close"] == pytest.approx(100.0 + 50 * 0.01)  # first row read, never the 999 copy
    assert df.attrs["checks"]["duplicates"] == 1


def test_bitstamp_nonpositive_price_is_masked_and_counted(tmp_path):
    close = np.array([100.0, 100.5, 0.0, 101.0, -3.0, 101.2, 101.4])
    df = D.load_bitstamp_csv(_bitstamp(tmp_path, close))
    r = df["logret"].to_numpy()
    assert df.attrs["checks"]["nonpositive_prices"] == 2
    assert np.isnan(df["close"].iloc[2]) and np.isnan(df["close"].iloc[4])
    assert np.isnan(r[2]) and np.isnan(r[3]) and np.isnan(r[4]) and np.isnan(r[5])  # both sides, never -inf
    assert np.isfinite(r[1]) and np.isfinite(r[6])
    assert not np.isinf(r).any()


def test_bitstamp_start_end_still_half_open(tmp_path):
    close = 100 + np.arange(10.0)
    df = D.load_bitstamp_csv(_bitstamp(tmp_path, close), start="2025-06-15 15:08", end="2025-06-15 15:11")
    assert len(df) == 3


def _nifty(tmp_path, bad=None):
    rows = []
    t0 = pd.Timestamp("2024-01-01 09:15", tz="Asia/Kolkata")
    for i in range(375):
        t = t0 + pd.Timedelta(minutes=i)
        px = 20000.0 + i
        if bad is not None and i in bad:
            px = bad[i]
        rows.append(dict(Epoch=int(t.timestamp()), Timestamp=t.strftime("%Y-%m-%dT%H:%M:%S"),
                         Open=px, High=px, Low=px, Close=px))
    p = D.nifty_file_path("2024", tmp_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(p, index=False)


def test_nifty_nonpositive_price_is_masked_and_counted(tmp_path):
    _nifty(tmp_path, bad={100: 0.0, 200: -5.0})
    df = D.load_nifty("2024-01-01", "2024-01-02", root=tmp_path, download=False)
    assert df.attrs["checks"]["nonpositive_prices"] == 2
    r = df["logret"].to_numpy()
    assert not np.isinf(r).any()
    for i in (100, 101, 200, 201):
        assert np.isnan(r[i])
    assert np.nanmax(np.abs(r)) < 1e-3


def test_nifty_clean_file_reports_zero_nonpositive(tmp_path):
    _nifty(tmp_path)
    df = D.load_nifty("2024-01-01", "2024-01-02", root=tmp_path, download=False)
    assert df.attrs["checks"]["nonpositive_prices"] == 0 and df["logret"].notna().sum() == 374


# ---------------------------------------------------------------------------------------------------------
# 2. monitor warnings
# ---------------------------------------------------------------------------------------------------------


def test_short_calibration_window_warns_but_still_runs():
    with pytest.warns(UserWarning, match="cal_len=10 is short"):
        m = VolatilityMonitor(threshold=1000.0, cal_len=10)
    m.run(np.random.default_rng(0).standard_normal(200))
    assert m.segments  # it calibrated


def test_default_and_experiment_calibration_lengths_do_not_warn():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        VolatilityMonitor()
        VolatilityMonitor(cal_len=120)
        VolatilityMonitor(cal_len=MIN_ADVISED_CAL_LEN)


def test_monitor_that_never_calibrates_warns_that_it_was_blind():
    m = VolatilityMonitor(threshold=1000.0, cal_len=200)
    with pytest.warns(UserWarning, match="never finished calibrating"):
        ev = m.run(np.random.default_rng(1).standard_normal(150))
    assert ev == []
    m2 = VolatilityMonitor(threshold=1000.0, cal_len=200)
    r = np.zeros(1000)  # a frozen feed: 1000 finite returns, but a zero scale is never usable
    with pytest.warns(UserWarning, match="1000 finite"):
        m2.run(r)
    m3 = VolatilityMonitor(threshold=1000.0, cal_len=200)
    with pytest.warns(UserWarning, match="0 finite"):
        m3.run(np.full(1000, np.nan))


def test_monitor_that_calibrates_gives_no_blind_warning():
    m = VolatilityMonitor(threshold=1000.0, cal_len=200)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        m.run(np.random.default_rng(2).standard_normal(500))
    assert m.segments


# ---------------------------------------------------------------------------------------------------------
# 3. benchmark detectors skip NaN
# ---------------------------------------------------------------------------------------------------------


def _change_path(seed=0):
    rng = np.random.default_rng(seed)
    return np.r_[rng.standard_normal(100), 3 * rng.standard_normal(200)]


def test_gaussian_lr_skips_a_missing_return():
    r = _change_path()
    clean = Det.gaussian_lr_edetector(r, 1.0, threshold=100).alarm[0]
    r_nan = r.copy()
    r_nan[5] = np.nan
    al = Det.gaussian_lr_edetector(r_nan, 1.0, threshold=100).alarm[0]
    assert clean >= 100 and al >= 100  # still alarms after the change
    # skipping a value is the same as deleting it (alarm index shifts by one when the NaN is before it)
    dropped = Det.gaussian_lr_edetector(np.delete(r, 5), 1.0, threshold=100).alarm[0]
    assert al == dropped + 1


def test_page_cusum_skips_a_missing_value():
    r = _change_path()
    r_nan = r.copy()
    r_nan[5] = np.nan
    al = Det.page_cusum(r_nan**2, 1.0, 0.5, 20.0).alarm[0]
    dropped = Det.page_cusum(np.delete(r, 5) ** 2, 1.0, 0.5, 20.0).alarm[0]
    assert al >= 100 and al == dropped + 1


def test_benchmarks_unchanged_without_nan():
    # the only callers (E4, E11) pass simulated data without NaN; the statistic path must be bit-identical
    rng = np.random.default_rng(7)
    r = rng.standard_normal((50, 400)) * np.r_[np.ones(200), 2 * np.ones(200)]
    a = Det.gaussian_lr_edetector(r, 1.0, threshold=200).alarm
    b = Det.page_cusum(r**2, 1.0, 0.25, 15.0).alarm
    # reference implementation of the old loops (before the NaN skip)
    z2 = r**2
    ratios = np.array([1.25, 1.5, 2.0])
    S = np.zeros((50, 3))
    ref = np.full(50, -1)
    act = np.ones(50, bool)
    for t in range(400):
        f = np.exp(np.minimum(-np.log(ratios)[None, :] + 0.5 * (1 - 1 / ratios**2)[None, :] * z2[:, t:t + 1], 700.0))
        S = f * (S + 1.0)
        hit = act & (S.mean(1) >= 200)
        ref[hit] = t
        act &= ~hit
    assert np.array_equal(a, ref)
    G = np.zeros(50)
    ref2 = np.full(50, -1)
    act = np.ones(50, bool)
    for t in range(400):
        G = np.maximum(0.0, G + z2[:, t] - 1.25)
        hit = act & (G >= 15.0)
        ref2[hit] = t
        act &= ~hit
    assert np.array_equal(b, ref2)


# ---------------------------------------------------------------------------------------------------------
# 4. command line errors
# ---------------------------------------------------------------------------------------------------------


def test_cli_missing_data_gives_a_short_message(capsys, monkeypatch, tmp_path):
    monkeypatch.setattr(D, "DATA_ROOT", tmp_path)
    code = cli_main(["run", "--source", "binance", "--symbol", "ETHUSDT", "--start", "2031-01-01",
                     "--end", "2031-01-02"])
    err = capsys.readouterr().err
    assert code == 2 and "error: data file not found" in err and "--download" in err
    assert "Traceback" not in err


def test_cli_bad_budget_gives_a_short_message(capsys):
    code = cli_main(["simulate", "--budget", "0"])
    err = capsys.readouterr().err
    assert code == 2 and "--budget" in err


def test_cli_end_before_start_gives_a_short_message(capsys):
    code = cli_main(["run", "--source", "binance", "--start", "2025-01-02", "--end", "2025-01-01"])
    err = capsys.readouterr().err
    assert code == 2 and "--end" in err and "Traceback" not in err


@pytest.mark.parametrize("bad", ["abc", None, np.ones(3)])
def test_monitor_update_rejects_non_numbers_at_once(bad):  # audit A4-B-004
    m = VolatilityMonitor(threshold=1000.0, cal_len=200)
    with pytest.raises(TypeError, match="one number"):
        m.update(bad)
    assert m.t == -1 and m.buffer == []  # nothing was stored


def test_null_band_of_one_warns_that_the_decrease_side_is_off():  # audit A4-B-003
    with pytest.warns(UserWarning, match="decrease side"):
        VolatilityMonitor(null_band=1.0)
