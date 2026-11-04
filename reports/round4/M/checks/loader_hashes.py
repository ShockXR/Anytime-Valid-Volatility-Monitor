"""Hash every column of the frames the loaders return on the real files (run once with the old src, once with the new)."""
import sys, json, hashlib
import numpy as np, pandas as pd
import avmonitor, avmonitor.data as D
root = sys.argv[1]
from pathlib import Path
root = Path(root)

def h(df):
    out = {"index": hashlib.sha256(df.index.asi8.tobytes()).hexdigest(), "shape": list(df.shape), "columns": list(map(str, df.columns))}
    for c in df.columns:
        v = df[c].to_numpy()
        if v.dtype.kind in "fiub":
            out[c] = hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest()
        else:
            out[c] = hashlib.sha256("|".join(map(str, v)).encode()).hexdigest()
    out["checks"] = {k: v for k, v in df.attrs.get("checks", {}).items()}
    return out

res = {"avmonitor_file": avmonitor.__file__}
bs = root / "crypto" / "btcusd_bitstamp_1min_latest.csv"
for filled in ("missing", "zero"):
    res[f"bitstamp_all_{filled}"] = h(D.load_bitstamp_csv(bs, filled=filled))
    res[f"bitstamp_2025H2_{filled}"] = h(D.load_bitstamp_csv(bs, "2025-07-01", "2026-01-01", filled=filled))
res["nifty_all"] = h(D.load_nifty("2017-01-01", "2026-11-01", root=root, download=False))
res["nifty_2024MayJun_regular"] = h(D.regular_sessions(D.load_nifty("2024-05-01", "2024-07-01", root=root, download=False)))
for sym in D.CRYPTO_SYMBOLS:
    res[f"binance_{sym}"] = h(D.load_binance(sym, "2024-01-01", "2026-01-01", root=root, download=False))
json.dump(res, open(sys.argv[2], "w"), indent=1, default=str)
print("done", sys.argv[2])
