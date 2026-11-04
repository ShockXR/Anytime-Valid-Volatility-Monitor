"""Review round 11 (E25), check 4c: gate 2 cuts only at 00:00 UTC, but the status part's data end can fall inside a
day (E23.data_end is one minute after the last bar in the mirror's file). This reruns S2 on BTC over E23's known
period with data ending at 13:37 UTC on 1 July 2026 (known data, read from disk, no download) and compares its alarms
with the committed results/e23_alarms.csv before the cut, and the release outcome of E23's CPI release of 14 July
2026 is not involved. Usage: nice -n 10 python docs/review/round11_scripts/r11_midday_causality.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT / "src"))
import e25_sequential_events as E25  # noqa: E402

CUTS = ("2026-07-01 13:37", "2026-09-16 18:30")  # the second: half an hour into an FOMC release window
com = pd.read_csv(ROOT / "results/e23_alarms.csv")
com = com[(com.arm == "S2") & (com.coin == "BTC")].reset_index(drop=True)
com["time"] = pd.to_datetime(com.time, utc=True)
for cut in CUTS:
    per = E25.period_e23(cut)
    assert pd.Timestamp(per["end"]) <= pd.Timestamp("2026-10-25")
    r = E25.run_arm(("BTC", "S2", per))
    a = r["alarms"].reset_index(drop=True)
    c = com[com.time < pd.Timestamp(cut, tz="UTC")].reset_index(drop=True)
    same = len(a) == len(c) and np.array_equal(a["index"].to_numpy(), c["index"].to_numpy()) and \
        np.array_equal(a["direction"].to_numpy(), c["direction"].to_numpy())
    print(f"cut {cut} UTC: cut run {len(a)} alarms, committed full run {len(c)} alarms before the cut; "
          f"index and direction identical: {same}; run took {r['seconds']:.0f} s")
