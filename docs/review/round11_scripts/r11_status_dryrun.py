"""Review round 11 (E25), check 3: the status part end to end on KNOWN data only.
E25.part_status is pointed at E23's known period (event file us_macro_events_2026.csv, stream from 1 January 2026,
load from 12 December 2025); the downloader is replaced by a guard that only returns files already on disk for months
up to October 2026 (and fails otherwise); outputs go to the scratchpad given on the command line, never to results/.
Only S2 is run (G1 is exploratory). Two looks: cut-off 25 October 2026 (the whole known period) and 10 October 2026
(the stratum of the FOMC statement of 16 September 2026 then lacks its last control window, 14 October 2026). Each look
is compared with the independent strata of r11_strata_lib.py and the committed S2 alarms restricted to before the
look's data end (gate 2 and r11_midday_causality.py show the alarms before a time do not depend on later data).
Usage: nice -n 10 python docs/review/round11_scripts/r11_status_dryrun.py <scratch dir>"""
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
from avmonitor import data as Dm  # noqa: E402

UTC = timezone.utc
scratch = Path(sys.argv[1])
EVF = ROOT / "data/events/us_macro_events_2026.csv"


def guard_fetch(url, dest, overwrite=False, timeout=0):
    dest = Path(dest)
    month = dest.stem.rsplit("_", 1)[-1]
    if month > "2026-10" or not dest.exists():
        raise AssertionError(f"dry run: refused to fetch {url}")
    return dest  # never downloads, never overwrites


Dm.fetch = guard_fetch
E25.FUTURE_EVENTS = EVF
E25.STREAM_START, E25.LOAD_START = "2026-01-01", "2025-12-12"
E25.STATUS_ARMS = ("S2",)
ev = pd.read_csv(EVF, dtype=str, keep_default_na=False)
relrows = ev[ev.status.isin(["released", "rescheduled"])]
al = pd.read_csv(ROOT / "results/e23_alarms.csv")
al = al[(al.arm == "S2") & (al.direction == "up")]
for cut in ("2026-10-25", "2026-10-10"):
    out = scratch / f"look_{cut}"
    out.mkdir(parents=True, exist_ok=True)
    E25.OUT = out
    E25.part_status(2, cut)
    summ = pd.read_csv(out / f"e25_status_{cut}_summary.csv")
    strata = pd.read_csv(out / f"e25_status_{cut}_strata.csv")
    pend = pd.read_csv(out / f"e25_status_{cut}_pending.csv")
    d_end = datetime.fromisoformat(cut).replace(tzinfo=UTC)
    rels = [(d, t, h) for d, t, h in zip(relrows.date, relrows.event_type, relrows.release_time_et)
            if L.et_to_utc(L.to_date(d), h) >= datetime(2026, 1, 1, tzinfo=UTC)]
    st = L.complete_prefix(L.build_strata(rels, list(relrows.date), datetime(2026, 1, 1, tzinfo=UTC), d_end), ev, d_end)
    for coin in ("BTC", "ETH", "SOL"):
        a = al[al.coin == coin]
        t = pd.to_datetime(a.time, utc=True)
        up = np.sort(t[t < pd.Timestamp(d_end)].dt.as_unit("ns").astype("int64").to_numpy())
        rows = [(L.has_up(up, s["release_utc"]), L.has_up(up, s["release_utc"]) + sum(L.has_up(up, x) for _, x in s["kept"]),
                 s["n"]) for s in st if s["in_test"]]
        m = L.e_path(rows)
        mine_strata = [s["stratum"] for s in st if s["in_test"]]
        sub = strata[strata.coin == coin]
        row = summ[(summ.coin == coin) & (summ.cell == "fomc+cpi")].iloc[0]
        ok = (list(sub.stratum) == mine_strata and np.isclose(row.M, m[-1], rtol=1e-12)
              and list(sub.n) == [r[2] for r in rows] and list(sub.Z) == [r[0] for r in rows])
        print(f"look {cut} | {coin}: status strata {list(sub.stratum)[0]}..{list(sub.stratum)[-1]} ({len(sub)}), M "
              f"{row.M:.6f}, role {row.role}, threshold {row.threshold}; independent: {len(mine_strata)} strata, M "
              f"{m[-1]:.6f}; {'MATCH' if ok else 'MISMATCH'}")
    print(f"look {cut}: pending strata {list(pend.stratum)} ({'; '.join(pend.date)}); outcome columns in the pending "
          f"file: {sorted({'Z', 'K', 'S', 'M'} & set(pend.columns))}; files written: "
          f"{sorted(p.name for p in out.iterdir())}")
