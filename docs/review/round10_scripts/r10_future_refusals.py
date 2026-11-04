"""Review round 10 (E24): does the `future` part of experiments/e24_blind_guard.py refuse to run as stated?

Calls E24.future_refusals with a patched clock and with edited copies of the event file (in a temporary directory;
the committed file is not touched). No market data are read.

Run: python docs/review/round10_scripts/r10_future_refusals.py
"""

import sys
import tempfile
import types
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

import e24_blind_guard as E24  # noqa: E402

ORIG_EVENTS = E24.FUTURE_EVENTS


def with_clock(now: str):
    class FakeTS(pd.Timestamp):
        @classmethod
        def now(cls, tz=None):
            return pd.Timestamp(now, tz=tz)

    shim = types.ModuleType("pandas_shim")
    shim.__dict__.update({k: getattr(pd, k) for k in dir(pd) if not k.startswith("__")})
    shim.Timestamp = FakeTS
    E24.pd = shim


def refusals(ev: pd.DataFrame | None, now: str, data_end=None):
    with_clock(now)
    with tempfile.TemporaryDirectory() as d:
        if ev is None:
            E24.FUTURE_EVENTS = ORIG_EVENTS
        else:
            p = Path(d) / "events.csv"
            ev.to_csv(p, index=False)
            E24.FUTURE_EVENTS = p
        try:
            return E24.future_refusals(E24.FUTURE_END, data_end)
        finally:
            E24.FUTURE_EVENTS = ORIG_EVENTS
            E24.pd = pd


def main():
    ev = pd.read_csv(ORIG_EVENTS, dtype=str, keep_default_na=False)
    end = E24.FUTURE_END
    print("FUTURE_END", end, "| primary", E24.FUTURE_PRIMARY_ARM, E24.FUTURE_PRIMARY_DESIGN)
    r = refusals(None, "2026-11-03 12:00")
    print(f"(a) today, committed file: {len(r)} refusals; first: {r[0]}")
    r = refusals(None, "2027-06-02 12:00")
    print(f"(b) after the end date, committed file: {len(r)} refusals (rows still 'scheduled')")
    upd = ev.copy()
    inside = pd.to_datetime(upd.date) < pd.Timestamp(end)
    upd.loc[inside & (upd.status == "scheduled"), "status"] = "released"
    upd.loc[upd.date == "2027-01-27", "status"] = "no_statement"
    r = refusals(upd, "2027-06-02 12:00")
    print(f"(c) after the end date, every row in the period updated: refusals {r}")
    r = refusals(upd, "2027-05-31 12:00")
    print(f"(c') the same file one day before the end date: refusals {r}")
    drop_fomc = upd[upd.date != "2027-03-17"]
    r = refusals(drop_fomc, "2027-06-02 12:00")
    print(f"(d) the FOMC row of 17 March 2027 deleted instead of updated: refusals {r}  <- not caught")
    drop_cpi = upd[upd.date != "2027-02-11"]
    r = refusals(drop_cpi, "2027-06-02 12:00")
    print(f"(e) the CPI row of 11 February 2027 deleted: refusals {r}")
    moved = upd.copy()
    moved.loc[moved.date == "2027-04-28", ["original_date", "date", "status"]] = ["2027-04-28", "2027-06-02", "rescheduled"]
    r = refusals(moved, "2027-06-02 12:00")
    print(f"(f) the April FOMC statement moved past the end date (rescheduled): refusals {r} (it simply leaves the test)")
    for de in ("2027-05-31 23:59", "2027-06-01 00:00"):
        r = refusals(upd, "2027-06-02 12:00", data_end=de)
        print(f"(g) data end (exclusive) {de}: refusals {r}")
    unsched = pd.concat([upd, pd.DataFrame([dict(date="2027-02-20", event_type="fomc", release_time_et="10:00",
                                                 release_time_utc="2027-02-20 15:00", status="scheduled",
                                                 original_date="", description="hypothetical unscheduled statement",
                                                 source_url="", verification="", access_date="")])])
    unsched = unsched.sort_values("date")
    r = refusals(unsched, "2027-06-02 12:00")
    print(f"(h) a hypothetical extra FOMC row left 'scheduled': refusals {r}")


if __name__ == "__main__":
    main()
