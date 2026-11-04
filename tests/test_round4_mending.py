"""Round 4 repairs on top of PR #4 (reports/round4/M/REPORT.md). [SYNTHETIC]

PR #4 (tests/test_input_hardening.py) covers the loaders' duplicate rule and close <= 0, the monitor's
warnings, update() on VolatilityMonitor, the NaN skip of page_cusum and gaussian_lr_edetector and the CLI
errors. The tests here cover only what round 4 adds: the threshold, c and c_down checks (A4-B-003,
R5-B-002), update() on the subclasses (A4-B-004), status() and the blind counters (A4-B-005), the NaN skip of
repeated_window_ztest (A1-B-007), infinite and missing prices in the loaders (A4-B-002), no OverflowError on
huge returns (R5-B-001), and CLI error paths that test_input_hardening.py does not reach. Seeds 9_950_000 up.
"""

from __future__ import annotations

import math
import urllib.error
import warnings

import numpy as np
import pandas as pd
import pytest

from avmonitor import cli
from avmonitor import data as D
from avmonitor import detectors as Det
from avmonitor.guard import GuardedVolatilityMonitor
from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor
from avmonitor.periodic import PeriodicVolatilityMonitor

SEED = 9_950_000

# --------------------------------------------------------------------------------------------
# A4-B-003 and R5-B-002: constructor checks
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("cls", [VolatilityMonitor, GaussianLRMonitor])
@pytest.mark.parametrize("bad", [math.inf, -math.inf, 0.5, 1.0, "abc", None])
def test_threshold_must_be_finite_and_above_one(cls, bad):
    with pytest.raises(ValueError, match="threshold"):
        cls(threshold=bad)


def test_threshold_just_above_one_is_accepted():
    assert VolatilityMonitor(threshold=1.0001).A == 1.0001


def test_infinite_threshold_no_longer_alarms_after_overflow():
    # R5-B-002: threshold=inf was accepted, and inf >= 2 * inf made an alarm fire once a statistic overflowed.
    with pytest.raises(ValueError, match="finite number above 1"):
        VolatilityMonitor(threshold=math.inf, cal_len=100)
    rng = np.random.default_rng(SEED)
    m = VolatilityMonitor(threshold=1e300, cal_len=100)  # the largest finite levels still work
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        m.run(np.r_[rng.standard_normal(100), np.full(500, 1e6)])
    assert all(e.stat >= 2e300 for e in m.events)


@pytest.mark.parametrize("bad", [None, "abc", math.nan, 0.0])
def test_c_errors_name_c(bad):
    with pytest.raises(ValueError, match=r"^c must be a finite"):
        VolatilityMonitor(c=bad)


@pytest.mark.parametrize("bad", [math.inf, math.nan, 0.0, -2.0, "x"])
def test_c_down_must_be_finite_and_positive(bad):
    with pytest.raises(ValueError, match=r"^c_down must be a finite"):
        VolatilityMonitor(c_down=bad)


# --------------------------------------------------------------------------------------------
# A4-B-004: update() of the subclasses
# --------------------------------------------------------------------------------------------


def _subclass_monitors():
    return [GaussianLRMonitor(threshold=100.0, cal_len=200),
            PeriodicVolatilityMonitor(threshold=100.0, cal_len=200, period=100),
            PeriodicVolatilityMonitor(threshold=100.0, cal_len=200, statistic="cusum")]


@pytest.mark.parametrize("bad", ["abc", None, np.ones(3), [0.1, 0.2]])
def test_subclass_update_rejects_non_numbers_at_once(bad):
    for m in _subclass_monitors():
        with pytest.raises(TypeError, match="one number"):
            m.update(bad)
        assert m.t == -1 and m.buffer == [], type(m).__name__


def test_gaussian_lr_monitor_update_with_python_floats_gives_exactly_the_same_alarms():
    rng = np.random.default_rng(SEED + 1)
    r = np.r_[rng.standard_normal(3000), 1.8 * rng.standard_normal(1500), 0.5 * rng.standard_normal(1500)]
    a = GaussianLRMonitor(threshold=500.0, cal_len=500).run(r)
    m = GaussianLRMonitor(threshold=500.0, cal_len=500)
    for x in r.tolist():
        m.update(x)
    assert len(a) >= 2
    assert [(e.index, e.direction, e.stat) for e in a] == [(e.index, e.direction, e.stat) for e in m.events]


# --------------------------------------------------------------------------------------------
# A4-B-005: status() and the blind counters
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("stream", ["nan", "zero"])
def test_status_shows_a_blind_stream_and_its_end(stream):
    x = np.full(2000, np.nan) if stream == "nan" else np.zeros(2000)
    m = VolatilityMonitor(threshold=100.0, cal_len=100)  # cal_max = 200: 10 restarts
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # PR #4's "never finished calibrating" warning
        m.run(x)
    s = m.status()
    assert s["state"] == "calibrating" and s["blind"] is True
    assert s["restarts_in_a_row"] == s["restarts_total"] == 10 == m.cal_stats["restarted"]
    assert s["steps_calibrating"] == m.cal_steps == 2000 and s["steps_monitoring"] == 0 and s["alarms"] == 0
    rng = np.random.default_rng(SEED + 2)
    m.run(rng.standard_normal(300))  # good data ends the blind spell
    s = m.status()
    assert s["restarts_in_a_row"] == 0 and s["restarts_total"] == 10 and s["calibrations"] >= 1
    assert s["steps_calibrating"] + s["steps_monitoring"] == 2300


def test_status_on_ordinary_data_and_before_any_data():
    m = VolatilityMonitor(threshold=1e6, cal_len=1000, null_band=0.25)
    s = m.status()
    assert s["step"] == -1 and s["state"] == "calibrating" and s["blind"] is False and s["window"] == 0
    rng = np.random.default_rng(SEED + 3)
    m.run(rng.standard_normal(300))
    assert m.status()["window"] == 300 and m.status()["steps_calibrating"] == 300
    m.run(rng.standard_normal(2700))
    s = m.status()
    assert s == dict(s, state="monitoring", blind=False, up_can_gain=True, down_active=True, steps_calibrating=1000,
                     steps_monitoring=2000, window=0, restarts_in_a_row=0, alarms=0)
    assert s["tau2_over_m_up"] > 1


def test_status_shows_an_increase_side_that_cannot_gain():
    # The blind confidence-sequence window of tests/test_guard.py: every X equals tau^2, so U >= tau^2.
    rng = np.random.default_rng(SEED + 4)
    cal, nt = 200, 50
    u = np.r_[rng.choice([-1.0, 1.0], nt), rng.choice([-10.0, 10.0], cal - nt)]
    m = VolatilityMonitor(threshold=1000.0, c=3.0, cal_len=cal, null_band=0.25, baseline="cs", delta=0.05)
    m.run(u)
    s = m.status()
    assert s["state"] == "monitoring" and s["up_can_gain"] is False and s["tau2_over_m_up"] <= 1.0


def test_status_works_for_the_subclasses():
    rng = np.random.default_rng(SEED + 5)
    r = rng.standard_normal(600)
    p = PeriodicVolatilityMonitor(threshold=1000.0, baseline="known", tau=3.0, m0=0.97, period=100)
    p.run(r)
    assert p.status()["steps_calibrating"] == 0 and p.status()["state"] == "monitoring"
    g = GuardedVolatilityMonitor(threshold=1000.0, c=3.0, cal_len=200, baseline="cs", delta=0.05)
    g.run(r)
    assert g.status()["steps_calibrating"] == 200 and g.status()["steps_monitoring"] == 400
    lr = GaussianLRMonitor(threshold=1e6, cal_len=100)
    lr.run(r)
    assert lr.status()["steps_calibrating"] == 100


# --------------------------------------------------------------------------------------------
# A1-B-007: repeated_window_ztest skips missing values
# --------------------------------------------------------------------------------------------


def _ztest_reference(r, window, level):
    """The old vectorised rule, written out: first t with |z_t| >= crit over the last ``window`` returns."""
    from scipy import stats

    crit = stats.norm.ppf(1 - level / 2)
    out = []
    for row in np.atleast_2d(r):
        a = -1
        for t in range(window - 1, row.size):
            z = (np.sum(row[t - window + 1:t + 1] ** 2) - window) / np.sqrt(window * 2.0)
            if abs(z) >= crit:
                a = t
                break
        out.append(a)
    return np.array(out)


def _paths(seed, reps=8, n=1200):
    rng = np.random.default_rng(seed)
    r = rng.standard_normal((reps, n))
    r[:, n // 2:] *= 1.6
    return r


def test_ztest_without_missing_values_matches_the_old_rule():
    r = _paths(SEED + 10, reps=6, n=900)
    a = Det.repeated_window_ztest(r, 1.0, 2.0, 60, level=1e-4).alarm
    assert np.array_equal(a, _ztest_reference(r, 60, 1e-4)) and (a >= 0).sum() >= 3


@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_ztest_skips_a_missing_value_as_if_it_were_dropped(bad):
    r = _paths(SEED + 11)
    rng = np.random.default_rng(SEED + 12)
    for i in range(r.shape[0]):
        x = r[i].copy()
        pos = np.sort(rng.choice(np.arange(5, 600), size=4, replace=False))
        x[pos] = bad
        keep = np.setdiff1d(np.arange(x.size), pos)
        a_nan = Det.repeated_window_ztest(x, 1.0, 2.0, 60, level=1e-4).alarm[0]
        a_drop = Det.repeated_window_ztest(r[i, keep], 1.0, 2.0, 60, level=1e-4).alarm[0]
        assert a_drop >= 0 and a_nan == keep[a_drop]
    x = r.copy()  # several rows, each with its own missing positions, in one call
    x[np.arange(8), np.arange(8) * 30 + 10] = bad
    rows = [Det.repeated_window_ztest(x[i], 1.0, 2.0, 60, level=1e-4).alarm[0] for i in range(8)]
    assert Det.repeated_window_ztest(x, 1.0, 2.0, 60, level=1e-4).alarm.tolist() == rows
    assert Det.repeated_window_ztest(x, 1.0, 2.0, 60, level=1e-4, start=100).alarm.min() >= 100


def test_ztest_with_fewer_observed_returns_than_a_window_never_alarms():
    x = np.full(200, np.nan)
    x[:30] = 5.0
    assert Det.repeated_window_ztest(x, 1.0, 2.0, 60).alarm.tolist() == [-1]


# --------------------------------------------------------------------------------------------
# A4-B-002: infinite and missing prices in the loaders
# --------------------------------------------------------------------------------------------


def _bitstamp_file(path, closes):
    ts0 = 1_700_000_000 - 1_700_000_000 % 60
    rows = [(ts0 + 60 * i, c, c, c, c, 1.0) for i, c in enumerate(closes)]
    pd.DataFrame(rows, columns=["timestamp", "open", "high", "low", "close", "volume"]).to_csv(path, index=False)


@pytest.mark.parametrize("filled", ["missing", "zero"])
def test_bitstamp_masks_infinite_and_missing_closes_and_counts_each_bar_once(tmp_path, filled):
    closes = [100, 101, np.inf, 102, 103, -5, 104, 105, np.nan, 106, 107, 0, 108]
    f = tmp_path / "bad.csv"
    _bitstamp_file(f, closes)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # no RuntimeWarning from log(0), log(-5) or inf - inf
        df = D.load_bitstamp_csv(f, filled=filled)
    lr = df["logret"].to_numpy()
    assert not np.isinf(lr).any()
    for i in (2, 5, 8, 11):
        assert np.isnan(lr[i]) and np.isnan(lr[i + 1]) and df[["open", "high", "low", "close"]].iloc[i].isna().all()
    assert np.isfinite(lr[[1, 4, 7, 10]]).all()
    ch = df.attrs["checks"]
    assert ch["nonfinite_prices"] == 2 and ch["nonpositive_prices"] == 2  # inf and NaN; -5 and 0


def test_nifty_masks_an_infinite_close(tmp_path):
    t0 = int(pd.Timestamp("2024-05-02 09:15", tz="Asia/Kolkata").timestamp())
    closes = [22000.0 + i for i in range(20)]
    closes[6] = np.inf
    d = tmp_path / "nifty"
    d.mkdir()
    pd.DataFrame([(t0 + 60 * i, "", c, c, c, c) for i, c in enumerate(closes)],
                 columns=["Epoch", "Timestamp", "Open", "High", "Low", "Close"]).to_csv(d / "NIFTY50_1min_2024.csv",
                                                                                      index=False)
    df = D.load_nifty("2024-05-01", "2024-05-03", root=tmp_path, download=False)
    lr = df["logret"].to_numpy()
    assert not np.isinf(lr).any() and np.isnan(lr[[6, 7]]).all() and np.isfinite(lr[[5, 8]]).all()
    assert df.attrs["checks"]["nonfinite_prices"] == 1 and df.attrs["checks"]["nonpositive_prices"] == 0


def test_binance_masks_an_infinite_close_but_does_not_count_missing_minutes(tmp_path):
    sym = "TESTUSDT"
    t_ms = 1_704_067_200_000 + 60_000 * np.r_[np.arange(15), np.arange(16, 30)]  # 2024-01-01, minute 15 missing
    close = 100.0 + np.arange(29.0)
    close[10] = np.inf
    d = tmp_path / "crypto" / "binance" / sym
    d.mkdir(parents=True)
    pd.DataFrame(dict(timestamp=t_ms, open=close, high=close, low=close, close=close, volume=1.0,
                      trade_count=1)).to_csv(d / f"{sym}_1m_2024-01.csv", index=False)
    df = D.load_binance(sym, "2024-01-01", "2024-01-01 00:30", root=tmp_path, download=False)
    ch = df.attrs["checks"]
    assert ch["nonfinite_prices"] == 1 and ch["missing_minutes"] == 1 and ch["nonpositive_prices"] == 0
    assert not np.isinf(df["logret"].to_numpy()).any() and np.isnan(df["logret"].to_numpy()[[10, 11]]).all()


# --------------------------------------------------------------------------------------------
# R5-B-001: no OverflowError on huge returns
# --------------------------------------------------------------------------------------------


def test_calibration_on_huge_returns_takes_the_degenerate_path_not_an_overflow_error():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        warnings.simplefilter("ignore", UserWarning)
        m = VolatilityMonitor(cal_len=2)
        m.run([1e300, 1e300, 1.0, 1.0])  # tau^2 = inf is not usable: the window is extended, not a crash
        assert m.cal_stats["extended"] >= 1 or m.cal_stats["restarted"] >= 1
        m = VolatilityMonitor(cal_len=4, baseline="cs")
        m.run([1e300] * 4 + [1.0] * 4)
        m = VolatilityMonitor(cal_len=200, c_down=1.5)
        m.run(np.full(400, 1e300))
        assert m.status()["restarts_total"] == 1 and not m.calibrated
        b = Det.calibrate(np.full(10, 1e300))
        assert b.tau > 0 and not math.isfinite(b.m0)
        g = GaussianLRMonitor(cal_len=2)
        g.run([1.0, -1.0])
        g.update(1e200)  # a Python float: (r / scale) ** 2 used to raise


def test_squares_of_ordinary_values_are_exactly_as_before():
    from avmonitor.monitor import _sq

    rng = np.random.default_rng(SEED + 20)
    v = rng.standard_normal(10_000) * 10.0 ** rng.integers(-150, 150, 10_000)
    assert all(_sq(float(x)) == float(x) ** 2 for x in v)
    assert _sq(1e200) == math.inf


# --------------------------------------------------------------------------------------------
# CLI paths that tests/test_input_hardening.py does not reach (the handler is PR #4's)
# --------------------------------------------------------------------------------------------


def _run_args(*extra):
    return ["run", "--source", "binance", "--symbol", "BTCUSDT", "--start", "2025-10-08", "--end", "2025-10-10",
            *extra]


@pytest.mark.parametrize("exc", [
    pd.errors.EmptyDataError("No columns to parse from file"),
    urllib.error.HTTPError("https://example.invalid/x.csv", 404, "Not Found", None, None),
    urllib.error.URLError("offline"),
    KeyError("timestamp"),
])
def test_cli_data_errors_give_one_line_and_exit_code_two(capsys, monkeypatch, exc):
    def boom(*a, **k):
        raise exc

    monkeypatch.setattr(D, "load_binance", boom)
    code = cli.main(_run_args())
    err = capsys.readouterr().err
    assert code == 2 and err.startswith("error: ") and len(err.strip().splitlines()) == 1 and "Traceback" not in err


@pytest.mark.parametrize("argv", [["run", "--start", "not-a-date", "--end", "2025-10-14"],
                                  ["estimators", "--day", "2025-13-45"]])
def test_cli_bad_date_gives_one_line_and_exit_code_two(capsys, argv):
    code = cli.main(argv)
    err = capsys.readouterr().err
    assert code == 2 and err.startswith("error: ") and "Traceback" not in err


def test_cli_monitor_parameter_error_gives_one_line(capsys, monkeypatch):
    rng = np.random.default_rng(SEED + 30)
    idx = pd.date_range("2025-09-28", periods=12 * 1440, freq="1min", tz="UTC")
    lr = 1e-3 * rng.standard_normal(idx.size)
    monkeypatch.setattr(cli, "_load", lambda a: pd.DataFrame(dict(close=100 * np.exp(np.cumsum(lr)), logret=lr),
                                                             index=idx))
    code = cli.main(_run_args("--lookback-days", "5", "--c", "nan"))
    err = capsys.readouterr().err
    assert code == 2 and "c must be a finite" in err and "Traceback" not in err
