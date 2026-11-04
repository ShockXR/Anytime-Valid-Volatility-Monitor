"""Review round 10 (E24): every number of docs/report_parts/E24.md and docs/iterations/E24.md that comes from a
committed CSV, recomputed from that CSV. Prints each check with the report's value and the recomputed one.

Run: python docs/review/round10_scripts/r10_report_numbers.py   [REAL and SYNTHETIC CSVs; seconds]
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
R = ROOT / "results"


def show(label, value):
    print(f"{label}: {value}")


def main():
    c = pd.read_csv(R / "e24_real_contrast.csv")
    for per in c.period.unique():
        for arm in ("P1", "G1", "G0", "S2"):
            x = c[(c.period == per) & (c.arm == arm) & (c.coin == "BTC")].set_index("event_type")
            if x.empty:
                continue
            p = x.loc["fomc+cpi"]
            show(f"{per} {arm} BTC pooled", f"{int(p.events_with)}/{int(p.events)}, control {p.rate_control_matched:.3f}, "
                 f"p {p.p_stratified_permutation:.3g}; FOMC {int(x.loc['fomc'].events_with)}/{int(x.loc['fomc'].events)}, "
                 f"CPI {int(x.loc['cpi'].events_with)}/{int(x.loc['cpi'].events)}; watched "
                 f"{int(x.loc['fomc'].watched_at_release + x.loc['cpi'].watched_at_release)}; blind at release "
                 f"{int(x.loc['fomc'].blind_at_release + x.loc['cpi'].blind_at_release)}")
    for coin in ("ETH", "SOL"):
        x = c[(c.period == "2026 (E23 design)") & (c.coin == coin) & (c.event_type == "fomc+cpi")].set_index("arm")
        show(f"2026 {coin} pooled G1 vs P1", f"G1 {int(x.loc['G1'].events_with)}/16 p {x.loc['G1'].p_stratified_permutation:.3f}; "
             f"P1 {int(x.loc['P1'].events_with)}/16")

    # the six BTC releases that differ between G1 and P1
    w = pd.read_csv(R / "e24_real_windows.csv")
    ref = pd.concat([pd.read_csv(R / "e22_windows.csv").assign(period="2024-2025 (E22 design)"),
                     pd.read_csv(R / "e23_windows.csv").assign(period="2026 (E23 design)")])
    g1 = w[(w.arm == "G1") & (w.coin == "BTC")][["period", "date", "event_type", "up_post", "monitoring_at_start"]]
    p1 = ref[(ref.arm == "P1") & (ref.coin == "BTC")][["period", "date", "event_type", "up_post", "monitoring_at_start",
                                                        "tau2_over_m_up_at_release"]]
    m = g1.merge(p1, on=["period", "date", "event_type"], suffixes=("_G1", "_P1"))
    d = m[m.up_post_G1 != m.up_post_P1]
    print("BTC releases where G1 and P1 differ:")
    print(d.to_string(index=False))

    s = pd.read_csv(R / "e24_real_summary.csv")
    g = s[s.arm == "G1"]
    show("G1 period alarm share range", f"{g.period_alarm_share.min():.3f} to {g.period_alarm_share.max():.3f}")
    show("G1 alarms per monitored day range", f"{g.alarms_per_monitored_day.min():.3f} to {g.alarms_per_monitored_day.max():.3f}")
    show("G1 first-day share (first_day_alarms / segments_with_a_day)",
         ", ".join(f"{r.period[:9]} {r.coin} {r.first_day_alarms / r.segments_with_a_day:.3f}" for r in g.itertuples()))
    show("G1 BTC discarded / extra minutes", ", ".join(f"{r.period[:9]}: {r.windows_discarded}, {r.extra_calibration_minutes}, "
                                                   f"days {r.days}" for r in g[g.coin == "BTC"].itertuples()))
    tot_days = g[g.coin == "BTC"].days.sum()
    tot_w = g[g.coin == "BTC"].windows_discarded.sum()
    show("one discarded window per (days)", f"{tot_days / tot_w:.0f} days = {tot_days / tot_w / 30.44:.1f} months")

    e = pd.read_csv(R / "e24_e19_blind.csv")
    e19 = pd.read_csv(R / "e19_blind.csv")
    e19 = e19[(e19.part == "p2") & (e19.baseline == "cs")]
    un = e[e.monitor == "unguarded"].reset_index(drop=True)
    cols = ["alarms", "segments", "segments_up_off", "monitored_minutes", "up_off_minutes", "up_off_share_of_monitored",
            "longest_up_off_days", "down_off_minutes"]
    mm = un.merge(e19, on=["source", "period", "band"], suffixes=("", "_e19"))
    mism = sum(int(not np.isclose(mm[c], mm[c + "_e19"]).all()) for c in cols)
    show("E19 diagnostic: unguarded rows matched / columns with a mismatch", f"{len(mm)} of {len(un)} / {mism}")

    gate = pd.read_csv(R / "e24_future_gate.csv")
    show("future gate rows, columns, mismatches", f"{len(gate)}, {gate.compared.max()}, {int(gate.mismatches.sum())}")

    sel = pd.read_csv(R / "e24_selection.csv")
    show("selection rows (config, K, scheme, p_kept)", "")
    print(sel[["config", "K", "scheme", "p_single_window", "p_kept_window_fails", "formula_1_minus_prod"]].to_string(index=False))


if __name__ == "__main__":
    main()
