"""Review round 11 (E25), check 5: calendar numbers of docs/report_parts/E25.md recomputed with the independent strata
of r11_strata_lib.py (planned forward calendar to December 2027, and F24a's planning calendar with the 2028
placeholders), under the amended and the original rules: candidates kept, n_i counts, balance, one-sided strata,
strata without controls, overlapping windows of different strata (any kind), close release pairs.
Usage: python docs/review/round11_scripts/r11_report_numbers.py"""
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT / "src"))
import r11_strata_lib as L  # noqa: E402

UTC = timezone.utc
START = datetime(2026, 10, 25, tzinfo=UTC)


def planned():
    ev = pd.read_csv(ROOT / "data/events/us_macro_events_e24_future.csv", dtype=str, keep_default_na=False)
    ev = ev[ev.status.isin(["released", "rescheduled", "scheduled"]) & (ev.date != "2027-01-27")]
    return [(d, t, h) for d, t, h in zip(ev.date, ev.event_type, ev.release_time_et)]


def with_placeholders():
    import e24_blind_guard as E24  # calendar only

    ev = E24.planning_events("2029-01-01")
    return [(d.strftime("%Y-%m-%d"), t, h) for d, t, h in zip(ev.date, ev.event_type, ev.release_time_et)]


def overlaps(st):
    w = []
    for s in st:
        w.append((s["release_utc"], s["stratum"], "release"))
        w += [(a, s["stratum"], "control") for _, a in s["kept"]]
    w.sort()
    pairs = []
    for i, (a, s, k) in enumerate(w):
        for b, s2, k2 in w[i + 1:]:
            if b >= a + L.WIN:
                break
            if s2 != s:
                pairs.append((k, k2))
    return pairs


for label, rels in (("planned to Dec 2027", planned()), ("with 2028 placeholders", with_placeholders())):
    for amended in (True, False):
        st = L.build_strata(rels, [r[0] for r in rels], START, None, amended=amended)
        n = np.array([s["n"] for s in st])
        bal = sum(1 for s in st if s["before"] == s["after"] > 0)
        one = sum(1 for s in st if (s["before"] == 0) != (s["after"] == 0))
        after_only = sum(1 for s in st if s["before"] == 0 and s["after"] > 0)
        ov = overlaps(st)
        print(f"{label} | {'amended' if amended else 'original'}: {len(st)} strata; kept {n.sum()} of {8 * len(st)}; "
              f"n counts 0..8 {np.bincount(n, minlength=9).tolist()}; balanced {bal}; one side only {one} "
              f"(after only {after_only}); no controls {int((n == 0).sum())}; controls before/after "
              f"{sum(s['before'] for s in st)}/{sum(s['after'] for s in st)}; overlapping window pairs {len(ov)} "
              f"(kinds {sorted(set(ov))})")
    if label.startswith("planned"):
        for amended, f in ((True, ""), (False, "_original_rules")):
            st = L.build_strata(rels, [r[0] for r in rels], START, None, amended=amended)
            ref = pd.read_csv(ROOT / f"results/e25_strata_planned{f}.csv")
            ref = ref[ref.kept.astype(bool)]
            mine = sorted((s["stratum"], pd.Timestamp(a)) for s in st for _, a in s["kept"])
            theirs = sorted(zip(ref.stratum, pd.to_datetime(ref.anchor_utc, utc=True)))
            print(f"{label} | {'amended' if amended else 'original'}: kept control anchors (UTC) identical to "
                  f"results/e25_strata_planned{f}.csv: {mine == theirs} ({len(mine)} anchors)")
    days = sorted((L.et_to_utc(L.to_date(d), h), L.to_date(d), t) for d, t, h in rels)
    close = [(str(a[1]), a[2], str(b[1]), b[2]) for a, b in zip(days, days[1:]) if (b[1] - a[1]).days <= 2]
    print(f"{label}: consecutive releases at most 2 Eastern days apart: {len(close)}: {close}")
