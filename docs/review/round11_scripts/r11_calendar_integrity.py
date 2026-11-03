"""Review round 11 (E25), check 3b [SYNTHETIC calendar, synthetic alarms; no market data]: rule 4 asks that no event-file
row within 30 days of the release still say 'scheduled', but nothing checks that the rows exist. F24a's future part
refuses when a month has no CPI row or a scheduled FOMC row is missing (E24.future_refusals); E25's status part does
not. Demonstration with E25's own functions: the FOMC statement of 28 October 2026 is released, the CPI release of
10 November 2026 is MISSING from the file. At a look with data to 1 December 2026, stratum 1 is complete and is
reported with four controls (11 November, one day from the missing release, is kept). When the row is added later
(status 'released'), the same stratum at a later look has three controls and a different M_1 for the same alarms.
Usage: python docs/review/round11_scripts/r11_calendar_integrity.py"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT / "src"))
import e25_sequential_events as E25  # noqa: E402


def events(rows):
    return pd.DataFrame([dict(date=d, event_type=t, release_time_et="14:00" if t == "fomc" else "08:30", status=s)
                         for d, t, s in rows])


d_end = pd.Timestamp("2026-12-01", tz="UTC")
now = pd.Timestamp("2026-12-02", tz="UTC")
# synthetic up alarms: one in the release window, one in the 11 November control window (19:00 UTC + 1 h)
ups = pd.DatetimeIndex([pd.Timestamp("2026-10-28 19:00", tz="UTC"), pd.Timestamp("2026-11-11 20:00", tz="UTC")])
for label, ev in (("CPI row of 10 November missing", events([("2026-10-28", "fomc", "released")])),
                  ("CPI row of 10 November present", events([("2026-10-28", "fomc", "released"),
                                                             ("2026-11-10", "cpi", "released")]))):
    why = E25.status_refusals(ev, d_end, now=now)
    cands, comp = E25.fresh_plan(ev, d_end)
    strata, path, summ, pend = E25.status_tables(cands, comp, {("S2", "BTC"): ups}, d_end)
    s1 = strata[strata.stratum == 1].iloc[0]
    print(f"{label}: refusals {why or 'none'}; stratum 1 in the test with n = {s1.n} (controls {s1.controls}), "
          f"Z = {s1.Z}, K = {s1.K}, M_1 = {s1.M:.4f}")
