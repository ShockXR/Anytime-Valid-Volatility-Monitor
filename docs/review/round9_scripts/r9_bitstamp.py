"""Round 9, m5 [REAL]: filled (zero-volume) minutes in the Bitstamp BTC/USD file and the returns they produce.

Needs the raw file (python scripts/download_data.py); it is not committed. Period: 2025-07-01 to 2026-01-01
(the E8 period).
"""
import numpy as np

from avmonitor.data import DATA_ROOT, load_bitstamp_csv

df = load_bitstamp_csv(DATA_ROOT / "crypto" / "btcusd_bitstamp_1min_latest.csv", "2025-07-01", "2026-01-01")
zv = df["zero_volume"].to_numpy()
r = df["logret"].to_numpy()
flat = ((df.open == df.high) & (df.high == df.low) & (df.low == df.close)).to_numpy()
z = zv.astype(int)
d = np.diff(np.r_[0, z, 0])
runs = np.where(d == -1)[0] - np.where(d == 1)[0]
print(f"[REAL] Bitstamp BTC/USD 2025-07-01 to 2026-01-01: {len(df)} minutes, {zv.sum()} zero-volume "
      f"({zv.mean():.2%}), of which flat candles {np.sum(zv & flat)} and zero returns {np.sum(zv & (r == 0))}; "
      f"zero returns in all {np.sum(r == 0)}; NaN returns {np.sum(np.isnan(r))}; zero-volume runs {runs.size}, "
      f"longest {runs.max()} minutes")
