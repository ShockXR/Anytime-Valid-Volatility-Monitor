"""Download every real dataset used in the experiments and write data/MANIFEST.csv.

Usage:
    python scripts/download_data.py            # crypto Jul to Dec 2025 + 2024, Nifty 2017 to Oct 2026, Bitstamp BTC
    python scripts/download_data.py --quick    # only the paper's window (Jul to Dec 2025) for the five coins

All files are public and are fetched from raw.githubusercontent.com. They are not
committed (data/raw is git-ignored). The manifest records each file's size, row
count and SHA-256 so a later run can confirm it used the same bytes.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from avmonitor.data import (  # noqa: E402
    BITSTAMP_URL,
    CRYPTO_SYMBOLS,
    DATA_ROOT,
    download_binance,
    download_nifty,
    fetch,
    sha256,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    files = []
    crypto_start = "2025-07-01" if a.quick else "2024-01-01"
    for s in CRYPTO_SYMBOLS:
        print("Binance", s, flush=True)
        files += download_binance(s, crypto_start, "2026-01-01")
    if not a.quick:
        print("Nifty 50", flush=True)
        files += download_nifty("2017-01-01", "2026-11-01")
        print("Bitstamp BTC/USD", flush=True)
        files.append(fetch(BITSTAMP_URL, DATA_ROOT / "crypto" / "btcusd_bitstamp_1min_latest.csv"))
    rows = []
    for f in files:
        with open(f, "rb") as fh:
            n = sum(1 for _ in fh) - 1
        rows.append(dict(file=str(Path(f).relative_to(ROOT)), bytes=Path(f).stat().st_size, rows=n, sha256=sha256(f)))
    man = pd.DataFrame(rows)
    out = ROOT / "data" / "MANIFEST.csv"
    if out.exists():
        old = pd.read_csv(out)
        man = pd.concat([old[~old.file.isin(man.file)], man]).sort_values("file")
    man.to_csv(out, index=False)
    print(f"{len(rows)} files, manifest at {out}")


if __name__ == "__main__":
    main()
