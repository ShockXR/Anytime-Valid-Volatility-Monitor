"""B5. Integrity of the REAL Binance BTCUSDT 1-minute files (months given on the command line) and comparison with data/MANIFEST.csv.
Usage: python b5_data_integrity.py <scratch_root> <out.csv> 2024-01 2024-02 ...  (downloads with the package's own fetch into scratch_root)"""
import sys, hashlib
import numpy as np, pandas as pd
from pathlib import Path
from avmonitor import data as DA

root = Path(sys.argv[1]); out = sys.argv[2]; months = sys.argv[3:]
man = pd.read_csv("/home/user/anytime-valid-volatility-monitor/data/MANIFEST.csv")
man["base"] = man["file"].str.split("/").str[-1]
rows = []
for mo in months:
    per = pd.Period(mo, "M")
    p = DA.binance_month_path("BTCUSDT", per, root)
    try:
        DA.fetch(DA.binance_month_url("BTCUSDT", per), p)
    except Exception as e:
        rows.append(dict(month=mo, check="download", value=f"FAILED {type(e).__name__}: {str(e)[:80]}", status="BLOCKED")); continue
    raw = pd.read_csv(p)
    sha = DA.sha256(p); m = man[man.base == p.name]
    exp_sha = m.sha256.iloc[0] if len(m) else None; exp_rows = int(m.rows.iloc[0]) if len(m) else None
    def add(check, value, status): rows.append(dict(month=mo, check=check, value=value, status=status))
    add("file sha256 equals MANIFEST", f"{sha[:12]} vs {str(exp_sha)[:12]}", "VERIFIED" if sha == exp_sha else "FAILED")
    add("row count equals MANIFEST", f"{len(raw)} vs {exp_rows}", "VERIFIED" if len(raw) == exp_rows else "FAILED")
    ts = raw["timestamp"].to_numpy()
    add("timestamps strictly increasing in file order", bool(np.all(np.diff(ts) > 0)), "VERIFIED" if np.all(np.diff(ts) > 0) else "FAILED")
    add("duplicate timestamps", int(raw.timestamp.duplicated().sum()), "VERIFIED" if raw.timestamp.duplicated().sum() == 0 else "FAILED")
    add("timestamps on whole minutes", bool(np.all(ts % 60000 == 0)), "VERIFIED" if np.all(ts % 60000 == 0) else "FAILED")
    t = pd.to_datetime(ts, unit="ms", utc=True)
    full = pd.date_range(per.start_time.tz_localize("UTC"), (per + 1).start_time.tz_localize("UTC"), freq="1min", inclusive="left")
    miss = full.difference(pd.DatetimeIndex(t))
    add("missing minutes in month", int(len(miss)), "VERIFIED" if len(miss) == 0 else f"{len(miss)} missing")
    if len(miss):
        gaps = pd.Series(miss).diff().ne(pd.Timedelta(minutes=1)).cumsum(); add("number of separate gaps", int(gaps.nunique()), "info")
        add("longest gap minutes", int(pd.Series(miss).groupby(gaps.values).size().max()), "info")
    for c in ["open", "high", "low", "close"]:
        add(f"non-positive or NaN {c}", int(((raw[c] <= 0) | raw[c].isna()).sum()), "VERIFIED" if ((raw[c] <= 0) | raw[c].isna()).sum() == 0 else "FAILED")
    bad_ohlc = int(((raw.high < raw.low) | (raw.close > raw.high) | (raw.close < raw.low) | (raw.open > raw.high) | (raw.open < raw.low)).sum())
    add("OHLC inconsistent bars (high<low or open/close outside range)", bad_ohlc, "VERIFIED" if bad_ohlc == 0 else "FAILED")
    zv = int((raw.volume <= 0).sum()); add("zero-volume bars", zv, "info")
    flat = int(((raw.open == raw.high) & (raw.high == raw.low) & (raw.low == raw.close)).sum()); add("flat bars (O=H=L=C)", flat, "info")
    lr = np.log(raw.close).diff().to_numpy()[1:]
    consecutive = (np.diff(ts) == 60000)
    lr = lr[consecutive]
    sd = np.median(np.abs(lr)) * 1.4826
    add("zero returns share", float((lr == 0).mean()), "info")
    add("max |return| / robust sd", float(np.abs(lr).max() / sd), "info")
    add("returns beyond 20 robust sd (count)", int((np.abs(lr) > 20 * sd).sum()), "info")
    add("max abs return", float(np.abs(lr).max()), "info")
    # high/low jumps relative to the close (a fat-finger candle shows up in high/low before the close)
    wick = ((raw.high - raw.low) / raw.close).to_numpy(); add("max (high-low)/close", float(wick.max()), "info")
    # load through the package loader and compare the loader's own checks
    s, e = per.start_time.strftime("%Y-%m-%d"), (per + 1).start_time.strftime("%Y-%m-%d")
    df = DA.load_binance("BTCUSDT", s, e, root=root, download=False)
    add("package loader checks", str(df.attrs["checks"]).replace(",", ";")[:200], "info")
    add("loader NaN returns (gaps+first)", int(df.logret.isna().sum()), "info")
pd.DataFrame(rows).to_csv(out, index=False)
print("EXIT")
