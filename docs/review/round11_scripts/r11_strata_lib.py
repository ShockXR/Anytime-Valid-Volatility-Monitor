"""Review round 11 (E25): an independent implementation of the strata of docs/preregistration_e25_sequential.md,
section 4 with Amendment 1, and of the e-process of section 5. Written from the pre-registration text only; it imports
nothing from experiments/. Uses zoneinfo (not pandas' tz database calls) for US Eastern time."""

from __future__ import annotations

from datetime import date, datetime, time, timedelta, timezone
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

ET = ZoneInfo("America/New_York")
UTC = timezone.utc
WIN = timedelta(hours=6)
PSI = (2.0, 4.0, 8.0, 16.0)


def et_to_utc(d: date, hhmm: str) -> datetime:
    h, m = (int(x) for x in hhmm.split(":"))
    return datetime.combine(d, time(h, m), tzinfo=ET).astimezone(UTC)


def to_date(s) -> date:
    return pd.Timestamp(s).date()


def build_strata(releases, release_days, stream_start: datetime, data_end: datetime | None, amended: bool = True):
    """releases: list of (date, type, hhmm). Returns a list of strata dicts, in release-time order."""
    rel = sorted(((et_to_utc(to_date(d), h), to_date(d), t, h) for d, t, h in releases), key=lambda x: x[0])
    rdays = sorted({to_date(d) for d in release_days})
    min_anchor = stream_start + timedelta(minutes=1440 + 360)
    used = {}  # anchor -> stratum
    out = []
    for i, (tu, d, typ, hhmm) in enumerate(rel, start=1):
        kept, dropped = [], []
        for k in (-4, -3, -2, -1, 1, 2, 3, 4):
            cd = d + timedelta(days=7 * k)
            a = et_to_utc(cd, hhmm)
            if any(abs((cd - r).days) < 2 for r in rdays):
                dropped.append((k, a, "near release day"))
            elif a < min_anchor:
                dropped.append((k, a, "before start"))
            elif data_end is not None and a + WIN > data_end:
                dropped.append((k, a, "after data end"))
            elif a in used:
                dropped.append((k, a, "used"))
            elif amended and any(abs(a - r[0]) < WIN for r in rel):
                dropped.append((k, a, "overlaps release window"))
            elif amended and any(abs(a - b) < WIN for b in used):
                dropped.append((k, a, "overlaps earlier control"))
            else:
                kept.append((k, a))
                used[a] = i
        need = et_to_utc(d + timedelta(days=28), hhmm) + WIN
        out.append(dict(stratum=i, date=d, type=typ, release_utc=tu, kept=kept, dropped=dropped, n=len(kept),
                        before=sum(1 for k, _ in kept if k < 0), after=sum(1 for k, _ in kept if k > 0),
                        data_needed_until=need))
    return out


def complete_prefix(strata, event_rows: pd.DataFrame, data_end: datetime | None):
    """Rule 4: data reach t + 28 d + 6 h, and no row within 30 days of d is 'scheduled'; strata enter in order."""
    sched = [to_date(x) for x in event_rows.loc[event_rows.status == "scheduled", "date"]]
    ok_prefix = True
    for s in strata:
        dok = data_end is None or data_end >= s["data_needed_until"]
        sok = not any(abs((r - s["date"]).days) <= 30 for r in sched)
        s["complete"] = dok and sok
        ok_prefix = ok_prefix and s["complete"]
        s["in_test"] = ok_prefix
    return strata


def has_up(up_ns: np.ndarray, a: datetime) -> int:
    lo = int(pd.Timestamp(a).value)
    hi = int(pd.Timestamp(a + WIN).value)
    return int(np.searchsorted(up_ns, hi, "left") > np.searchsorted(up_ns, lo, "left"))


def e_val(z, s, n, psi):
    p = s / (n + 1)
    return psi**z / (1 - p + psi * p)


def e_path(rows):
    """rows: list of (Z, S, n). Returns M_j for j = 1..J."""
    prod = np.ones(len(PSI))
    m = []
    for z, s, n in rows:
        prod = prod * np.array([e_val(z, s, n, ps) for ps in PSI])
        m.append(prod.mean())
    return np.array(m)
