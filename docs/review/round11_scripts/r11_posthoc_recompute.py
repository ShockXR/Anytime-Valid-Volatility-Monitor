"""Review round 11 (E25), check 4b: recompute the post hoc e-process (pre-registration section 10) from the committed
S2 alarm files with the independent strata of r11_strata_lib.py, and compare with results/e25_posthoc_strata.csv and
results/e25_posthoc_summary.csv (amended rules) and the *_original_rules.csv files. Reads only committed CSVs.
Usage: python docs/review/round11_scripts/r11_posthoc_recompute.py"""

from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import r11_strata_lib as L  # noqa: E402

UTC = timezone.utc
PERIODS = [
    ("2024 to 2025 (E22 design)", ROOT / "data/events/us_macro_events.csv", datetime(2024, 1, 21, tzinfo=UTC),
     datetime(2026, 1, 1, tzinfo=UTC), "e22_alarms.csv"),
    ("2026 (E23 design)", ROOT / "data/events/us_macro_events_2026.csv", datetime(2026, 1, 1, tzinfo=UTC),
     datetime(2026, 10, 25, tzinfo=UTC), "e23_alarms.csv"),
]


def run(amended: bool):
    sfx = "" if amended else "_original_rules"
    ref = pd.read_csv(ROOT / f"results/e25_posthoc_strata{sfx}.csv")
    refs = pd.read_csv(ROOT / f"results/e25_posthoc_summary{sfx}.csv")
    bad, checked = 0, 0
    lines = []
    for plabel, evf, start, end, af in PERIODS:
        ev = pd.read_csv(evf, dtype=str, keep_default_na=False)
        relrows = ev[ev.status.isin(["released", "rescheduled"])]
        rels = []
        for d, t, h in zip(relrows.date, relrows.event_type, relrows.release_time_et):
            tu = L.et_to_utc(L.to_date(d), h)
            if start <= tu < end:
                rels.append((d, t, h))
        st = L.build_strata(rels, list(relrows.date), start, end, amended=amended)
        st = L.complete_prefix(st, ev, end)
        al = pd.read_csv(ROOT / "results" / af)
        al = al[(al.arm == "S2") & (al.direction == "up")]
        for coin in ("BTC", "ETH", "SOL"):
            up = np.sort(pd.to_datetime(al[al.coin == coin].time, utc=True).dt.as_unit("ns").astype("int64").to_numpy())
            rows = []
            for s in st:
                z = L.has_up(up, s["release_utc"])
                k = sum(L.has_up(up, a) for _, a in s["kept"])
                s_ = dict(stratum=s["stratum"], date=str(s["date"]), type=s["type"], n=s["n"], Z=z, K=k, S=z + k,
                          before=s["before"], after=s["after"], in_test=s["in_test"])
                rows.append(s_)
            mine = pd.DataFrame(rows)
            r = ref[(ref.period == plabel) & (ref.coin == coin)].sort_values("stratum").reset_index(drop=True)
            same_cols = ["stratum", "n", "in_test", "controls_before", "controls_after"]
            m2 = mine.rename(columns=dict(before="controls_before", after="controls_after"))
            eq = (len(r) == len(m2)) and all((r[c].astype(int).to_numpy() == m2[c].astype(int).to_numpy()).all()
                                              for c in same_cols)
            it = m2[m2.in_test]
            ri = r[r.in_test.astype(bool)]
            eq_out = len(it) == len(ri) and all((ri[c].astype(int).to_numpy() == it[c].to_numpy()).all()
                                                for c in ("Z", "K", "S"))
            for cell in ("fomc+cpi", "fomc", "cpi"):
                sub = it if cell == "fomc+cpi" else it[it.type == cell]
                m = L.e_path(list(zip(sub.Z, sub.S, sub.n))) if len(sub) else np.array([1.0])
                first = int(np.argmax(m >= 20)) + 1 if (m >= 20).any() else None
                rs = refs[(refs.period == plabel) & (refs.coin == coin) & (refs.cell == cell)].iloc[0]
                f_ref = None if pd.isna(rs.first_crossing_20_stratum) else int(rs.first_crossing_20_stratum)
                ok = (np.isclose(m[-1], rs.final_M, rtol=1e-12) and np.isclose(m.max(), rs.max_M, rtol=1e-12)
                      and first == f_ref and int(sub.Z.sum()) == int(rs.releases_with_outcome)
                      and int(sub.K.sum()) == int(rs.controls_with_outcome) and int(sub.n.sum()) == int(rs.control_windows)
                      and len(sub) == int(rs.strata_in_test))
                checked += 1
                bad += int(not ok)
                lines.append(f"{'amended' if amended else 'original'} | {plabel} | {coin} | {cell}: strata {len(sub)}, "
                             f"Z {int(sub.Z.sum())}, K {int(sub.K.sum())}/{int(sub.n.sum())}, final M {m[-1]:.6g} "
                             f"(CSV {rs.final_M:.6g}), max M {m.max():.6g}, first >= 20 at {first} (CSV {f_ref}); "
                             f"{'MATCH' if ok else 'MISMATCH'}")
            checked += 2
            bad += int(not eq) + int(not eq_out)
            lines.append(f"  strata table (n, in_test, before, after) equal to CSV: {eq}; outcomes Z, K, S equal: {eq_out}")
            if amended and coin == "BTC":
                lines.append("  BTC strata in the test: " + "; ".join(
                    f"{r_.stratum}:{r_.date} {r_.type} n={r_.n} Z={r_.Z} K={r_.K}" for r_ in it.itertuples()))
    return lines, checked, bad


if __name__ == "__main__":
    tot_c, tot_b = 0, 0
    for amended in (True, False):
        lines, c, b = run(amended)
        print("\n".join(lines))
        tot_c, tot_b = tot_c + c, tot_b + b
    print(f"TOTAL: {tot_c} comparisons, {tot_b} mismatches")
