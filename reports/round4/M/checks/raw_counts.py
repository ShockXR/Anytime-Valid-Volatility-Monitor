# Round 4, agent M: real-data counts for reports/round4/M (see the handback report).
"""Count duplicate timestamps, non-monotonic order, non-positive and non-finite prices in the raw files."""
import sys, glob, json
import numpy as np, pandas as pd
root = sys.argv[1]
out = []
def summ(name, ts, close, other=None):
    ts = pd.Series(ts)
    d = dict(file=name, rows=len(ts), duplicated_ts=int(ts.duplicated().sum()),
             monotonic=bool(ts.is_monotonic_increasing), decreasing_steps=int((ts.diff() < 0).sum()),
             nonpositive_close=int((close <= 0).sum()), nonfinite_close=int((~np.isfinite(close)).sum()))
    if other is not None:
        d["nonpositive_ohl"] = int((other <= 0).any(axis=1).sum()); d["nonfinite_ohl"] = int((~np.isfinite(other)).any(axis=1).sum())
    if d["duplicated_ts"]:
        # do duplicated rows differ in content?
        d["dup_rows_differ"] = None
    return d
b = pd.read_csv(f"{root}/crypto/btcusd_bitstamp_1min_latest.csv")
d = summ("bitstamp", b["timestamp"], b["close"].to_numpy(float), b[["open", "high", "low"]].to_numpy(float))
dup = b[b["timestamp"].duplicated(keep=False)]
d["dup_groups_with_different_rows"] = int(dup.groupby("timestamp").apply(lambda g: len(g.drop_duplicates()) > 1).sum()) if len(dup) else 0
d["first"] = str(pd.to_datetime(b["timestamp"].min(), unit="s")); d["last"] = str(pd.to_datetime(b["timestamp"].max(), unit="s"))
d["zero_volume_rows"] = int((b["volume"] <= 0).sum())
out.append(d)
for f in sorted(glob.glob(f"{root}/nifty/*.csv")):
    n = pd.read_csv(f); n.columns = [c.strip().lower() for c in n.columns]
    d = summ(f.split("/")[-1], n["epoch"], n["close"].to_numpy(float), n[["open", "high", "low"]].to_numpy(float))
    out.append(d)
tot = dict(file="binance_all", rows=0, duplicated_ts=0, monotonic=True, decreasing_steps=0, nonpositive_close=0, nonfinite_close=0, nonpositive_ohl=0, nonfinite_ohl=0, files=0)
for f in sorted(glob.glob(f"{root}/crypto/binance/*/*.csv")):
    k = pd.read_csv(f, usecols=["timestamp", "open", "high", "low", "close"])
    d = summ(f.split("/")[-1], k["timestamp"], k["close"].to_numpy(float), k[["open", "high", "low"]].to_numpy(float))
    tot["files"] += 1
    for key in ("rows", "duplicated_ts", "decreasing_steps", "nonpositive_close", "nonfinite_close", "nonpositive_ohl", "nonfinite_ohl"):
        tot[key] += d[key]
    tot["monotonic"] &= d["monotonic"]
out.append(tot)
df = pd.DataFrame(out)
df.to_csv(sys.argv[2], index=False)
print(df.to_string())
