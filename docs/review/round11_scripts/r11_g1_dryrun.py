"""Review round 11 (E25), check 5b: the iteration notes (iteration 5) quote a G1 number from a dry run whose outputs are
not committed ("G1 (exploratory) on BTC: 4 of 15, M = 0.55"). This reruns G1 on BTC over E23's known period (data to
24 October 2026 23:59 UTC, on disk, no download) through E25.run_arm and computes the amended-rule e-process with the
independent strata of r11_strata_lib.py. Usage: nice -n 10 python docs/review/round11_scripts/r11_g1_dryrun.py"""
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT / "src"))
import r11_strata_lib as L  # noqa: E402
import e25_sequential_events as E25  # noqa: E402

UTC = timezone.utc
per = E25.period_e23()
assert pd.Timestamp(per["end"]) <= pd.Timestamp("2026-10-25")
r = E25.run_arm(("BTC", "G1", per))
a = r["alarms"]
up = np.sort(pd.DatetimeIndex(a.loc[a.direction == "up", "time"]).as_unit("ns").asi8)
ev = pd.read_csv(ROOT / "data/events/us_macro_events_2026.csv", dtype=str, keep_default_na=False)
relrows = ev[ev.status.isin(["released", "rescheduled"])]
start, end = datetime(2026, 1, 1, tzinfo=UTC), datetime(2026, 10, 25, tzinfo=UTC)
rels = [(d, t, h) for d, t, h in zip(relrows.date, relrows.event_type, relrows.release_time_et)
        if start <= L.et_to_utc(L.to_date(d), h) < end]
st = L.complete_prefix(L.build_strata(rels, list(relrows.date), start, end), ev, end)
rows = [(L.has_up(up, s["release_utc"]), L.has_up(up, s["release_utc"]) + sum(L.has_up(up, x) for _, x in s["kept"]), s["n"])
        for s in st if s["in_test"]]
m = L.e_path(rows)
print(f"G1 BTC, E23 period, amended strata: {len(rows)} strata, releases with an up alarm {sum(z for z, _, _ in rows)}, "
      f"final M {m[-1]:.4f} (run {r['seconds']:.0f} s)")
