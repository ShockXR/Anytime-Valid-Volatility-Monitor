"""E7c. External check of the Nifty monitor with India VIX (NSE implied volatility index).

REAL MARKET DATA. Source: NSE daily index close files
https://archives.nseindia.com/content/indices/ind_close_all_DDMMYYYY.csv
(one file per trading day; the file has an "India VIX" row with open, high,
low and close, and a "Nifty 50" row used here to cross-check the minute data).

Research question
  Q7c. Are up alarms of the streaming monitor followed by India VIX rises, and
       down alarms by falls, more often than random sessions?

Usage
  python -I experiments/e07c_india_vix.py download   # fetch and cache the daily files (polite, 1 s pacing)
  python -I experiments/e07c_india_vix.py            # analysis (needs results/e07b_alarms.csv)
"""

import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/india_vix"
URL = "https://archives.nseindia.com/content/indices/ind_close_all_{d}.csv"


def download(pause: float = 1.0):
    """Fetch one file per Nifty session date (cached; 404 means no file for that day)."""
    from avmonitor.data import load_nifty

    RAW.mkdir(parents=True, exist_ok=True)
    df = load_nifty("2017-01-01", "2026-11-01", download=False)
    days = sorted(pd.DatetimeIndex(df["session"].unique()).strftime("%d%m%Y"))
    log = []
    for i, d in enumerate(days):
        dest = RAW / f"ind_close_all_{d}.csv"
        if dest.exists() and dest.stat().st_size > 0:
            log.append(dict(day=d, status="cached", bytes=dest.stat().st_size))
            continue
        status, size = "error", 0
        for attempt in range(3):
            try:
                with urllib.request.urlopen(URL.format(d=d), timeout=60) as resp:
                    body = resp.read()
                tmp = dest.with_suffix(".part")
                tmp.write_bytes(body)
                tmp.replace(dest)
                status, size = "ok", len(body)
                break
            except urllib.error.HTTPError as e:
                status = f"http {e.code}"
                if e.code == 404:
                    break
            except Exception as e:  # network hiccup: back off and retry
                status = f"error {type(e).__name__}"
            time.sleep(pause * 2 ** (attempt + 1))
        log.append(dict(day=d, status=status, bytes=size))
        if i % 100 == 0:
            print(i, d, status, flush=True)
        time.sleep(pause)
    pd.DataFrame(log).to_csv(RAW / "download_log.csv", index=False)
    print(pd.DataFrame(log).status.value_counts())


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "download":
        download()
