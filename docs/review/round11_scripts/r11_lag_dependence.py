"""Review round 11 (E25), check 1c: how strongly does S2's up-alarm outcome in one 6-hour window depend on the outcome
in a window of the previous day? Releases one day apart (FOMC Wednesday 14:00 ET then CPI Thursday 08:30 ET, or CPI
Tuesday 08:30 then FOMC Wednesday 14:00) occur four times on the planned calendar to December 2027, while rule 2 keeps
every control at least 2 days from every release. H0 conditions on the earlier strata, so a day-to-day dependence of
the outcome process is part of what H0 must assume away. [REAL, seen data, descriptive] Committed S2 up alarms (BTC,
ETH, SOL), E22's period (21 January 2024 to 31 December 2025) and E23's (1 January to 24 October 2026); pairs of
windows on consecutive Eastern days, both days at least 2 days from every release day of the period's event file
(so the pairs are free of the releases themselves).
Usage: python docs/review/round11_scripts/r11_lag_dependence.py"""
import sys
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import r11_strata_lib as L  # noqa: E402

UTC = timezone.utc
PER = [(ROOT / "data/events/us_macro_events.csv", date(2024, 1, 22), date(2025, 12, 30), "e22_alarms.csv"),
       (ROOT / "data/events/us_macro_events_2026.csv", date(2026, 1, 2), date(2026, 10, 23), "e23_alarms.csv")]
PAIRS = {"14:00 ET, then 08:30 ET next day (FOMC Wed / CPI Thu)": ("14:00", "08:30"),
         "08:30 ET, then 14:00 ET next day (CPI Tue / FOMC Wed)": ("08:30", "14:00")}

for coin in ("BTC", "ETH", "SOL"):
    for name, (h1, h2) in PAIRS.items():
        tab = np.zeros((2, 2), int)
        for evf, d0, d1, af in PER:
            ev = pd.read_csv(evf, dtype=str, keep_default_na=False)
            rdays = [L.to_date(x) for x in ev.loc[ev.status.isin(["released", "rescheduled"]), "date"]]
            al = pd.read_csv(ROOT / "results" / af)
            al = al[(al.arm == "S2") & (al.coin == coin) & (al.direction == "up")]
            up = np.sort(pd.to_datetime(al.time, utc=True).dt.as_unit("ns").astype("int64").to_numpy())
            d = d0
            while d < d1:
                e = d + timedelta(days=1)
                if d.weekday() < 5 and e.weekday() < 5 and all(abs((x - r).days) >= 2 for x in (d, e) for r in rdays):
                    x = L.has_up(up, L.et_to_utc(d, h1))
                    y = L.has_up(up, L.et_to_utc(e, h2))
                    tab[x, y] += 1
                d = e
        a = tab[0, 1] / tab[0].sum()
        b = tab[1, 1] / tab[1].sum() if tab[1].sum() else np.nan
        # Fisher-type summary: odds ratio with a 0.5 correction
        orr = ((tab[1, 1] + .5) * (tab[0, 0] + .5)) / ((tab[1, 0] + .5) * (tab[0, 1] + .5))
        se = np.sqrt(sum(1 / (v + .5) for v in tab.ravel()))
        print(f"{coin} | {name}: pairs {tab.sum()}; P(second = 1 | first = 0) = {a:.3f} ({tab[0,1]}/{tab[0].sum()}), "
              f"P(second = 1 | first = 1) = {b:.3f} ({tab[1,1]}/{tab[1].sum()}); odds ratio {orr:.2f} "
              f"(log OR {np.log(orr):+.2f} +- {1.96*se:.2f})")
