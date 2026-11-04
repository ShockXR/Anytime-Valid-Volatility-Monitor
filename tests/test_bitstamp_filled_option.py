"""Review round 9, m5: the Bitstamp loader's option that reproduces the committed results (filled minutes as zeros)."""

import numpy as np
import pandas as pd
import pytest

from avmonitor.data import load_bitstamp_csv


@pytest.fixture
def bitstamp_file(tmp_path):
    ts = 1_750_000_000 + 60 * np.arange(6)
    close = [100.0, 101.0, 101.0, 102.0, 101.5, 101.5]
    vol = [1.0, 2.0, 0.0, 3.0, 1.0, 0.5]  # minute 2 is a filled flat candle; minute 5 a real unchanged price
    frame = {"timestamp": ts, "open": close, "high": close, "low": close, "close": close, "volume": vol}
    pd.DataFrame(frame).to_csv(tmp_path / "b.csv", index=False)
    return tmp_path / "b.csv"


def test_default_is_missing(bitstamp_file):
    a = load_bitstamp_csv(bitstamp_file)["logret"].to_numpy()
    b = load_bitstamp_csv(bitstamp_file, filled="missing")["logret"].to_numpy()
    np.testing.assert_array_equal(a, b)
    assert np.isnan(a[2]) and np.isnan(a[3])


def test_zero_reproduces_the_old_loader(bitstamp_file):
    df = load_bitstamp_csv(bitstamp_file, filled="zero")
    close = np.array([100.0, 101.0, 101.0, 102.0, 101.5, 101.5])
    old = np.r_[np.nan, np.diff(np.log(close))]  # the loader before review round 9: log(close).diff()
    np.testing.assert_array_equal(df["logret"].to_numpy(), old)
    assert df["logret"].iloc[2] == 0.0  # the filled minute
    assert df["logret"].iloc[3] == np.log(102.0) - np.log(101.0)  # the next return carries the whole move
    assert list(df["zero_volume"]) == [False, False, True, False, False, False]


def test_unknown_option_raises(bitstamp_file):
    with pytest.raises(ValueError):
        load_bitstamp_csv(bitstamp_file, filled="ffill")
