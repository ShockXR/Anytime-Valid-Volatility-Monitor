"""Turn r3_null.csv and r3_change.csv into Markdown tables for REPORT.md (prints to stdout)."""
import sys
import numpy as np, pandas as pd

SHORT = {"e-detector, KNOWN baseline (theorem)": "Known baseline (theorem)", "monitor plugin band0.25 (README)": "Plug-in, band 0.25 (README)",
         "monitor plugin band0": "Plug-in, band 0", "monitor cs baseline delta0.01": "Confidence-sequence baseline"}

def cell(r):
    s = f"{r.rate:.3f} [{r.ci_low:.3f},{r.ci_high:.3f}]"
    if r.verdict == "EXCEEDS bound": s += "*"
    return s

def null_tables(path, groups, extra=None):
    df = pd.read_csv(path)
    if extra is not None: df = pd.concat([df, pd.read_csv(extra)], ignore_index=True)
    df = df[~df.horizon_m.astype(str).str.startswith("switch")].copy(); df["horizon_m"] = df.horizon_m.astype(int)
    for title, laws in groups:
        laws = [l for l in laws if l in set(df.law)]
        if not laws: continue
        for m in (390, 1950):
            sub = df[df.horizon_m == m]
            print(f"\n{title}: within the first {m} monitored steps (bound {m/7800:.3f}; cs bound {0.01+m/7800:.3f}). Cell = rate [Wilson 95%]; * = interval entirely above the bound; (u/d) = up/down alarms.\n")
            print("| Method | " + " | ".join(laws) + " |"); print("|---|" + "---|" * len(laws))
            for meth, sh in SHORT.items():
                row = []
                for l in laws:
                    r = sub[(sub.law == l) & (sub.method == meth)]
                    row.append(cell(r.iloc[0]) + f" ({int(r.iloc[0].k_up)}u/{int(r.iloc[0].k_down)}d)" if len(r) else "")
                print(f"| {sh} | " + " | ".join(row) + " |")
            print(f"paths per law: {int(sub[sub.law.isin(laws)].paths.min())} to {int(sub[sub.law.isin(laws)].paths.max())}")
    sw = pd.read_csv(path)
    if extra is not None: sw = pd.concat([sw, pd.read_csv(extra)], ignore_index=True)
    sw = sw[sw.horizon_m.astype(str).str.startswith("switch")]
    if len(sw):
        print("\nAfter the tail switch (monitored step 975), among paths with no alarm before it. Horizon counts from the switch; bound as if a fresh segment.\n")
        print("| Law | Method | Horizon | Alarms / paths | Rate [Wilson 95%] | up / down |"); print("|---|---|---|---|---|---|")
        for r in sw.itertuples():
            print(f"| {r.law} | {SHORT[r.method]} | {r.horizon_m} | {r.k} / {r.paths} | {cell(r)} | {r.k_up} / {r.k_down} |")

def change_tables(path):
    df = pd.read_csv(path)
    def fmt(r, key):
        return f"{r[key]:.3f} {r[key + '_ci']}" if key in r and pd.notna(r[key]) else ""
    print("\n(h) Ramp: variance rises linearly from 1 to 2 over 390 steps, starting at monitored step 50, then stays at 2.\n")
    print("| Method | Alarm before ramp | Detected during ramp (390 steps) | Detected within 780 steps of ramp start | Detected within 5 days | Median delay from ramp start (IQR) | Share up |"); print("|---|---|---|---|---|---|---|")
    for _, r in df[df.cell == "ramp"].iterrows():
        print(f"| {SHORT[r.method]} | {fmt(r,'alarm_before_change')} | {fmt(r,'detected_during_ramp')} | {fmt(r,'detected_within_day_after_ramp_end')} | {fmt(r,'detected_within_5days')} | {r.median_delay_from_ramp_start:.0f} ({r.delay_iqr}) | {r.share_up:.2f} |")
    print("\n(i) Burst: variance x3 for L steps starting at monitored step 200, then back to 1.\n")
    print("| Burst L | Method | Alarm before burst | Inside burst | Within 60 after | Inside or within 60 after | Later than 60 after (5-day horizon) | Median delay from burst start (IQR) | Share up |"); print("|---|---|---|---|---|---|---|---|---|")
    for _, r in df[df.cell.str.startswith("burst")].iterrows():
        print(f"| {r.cell.split('_')[1]} | {SHORT[r.method]} | {fmt(r,'alarm_before_burst')} | {fmt(r,'detected_inside_burst')} | {fmt(r,'detected_within_60_after_burst')} | {fmt(r,'detected_inside_or_within_60_after')} | {fmt(r,'alarm_later_than_60_after_burst')} | {r.median_delay_from_burst_start:.0f} ({r.delay_iqr}) | {r.share_up:.2f} |")
    print("\n(j) Change during calibration: variance x2 from calibration step c of 1,000 onwards (so the last 1000 - c calibration returns already have the new variance).\n")
    print("| c | Method | Alarm within 390 | Alarm within 1950 | Up within 1950 | Down within 1950 | Median delay from monitoring start (IQR) | Median m0hat / true new m0 |"); print("|---|---|---|---|---|---|---|---|")
    for _, r in df[df.cell.str.startswith("incal")].iterrows():
        mm = r.get("median_m0hat_over_true_new_m0", np.nan)
        print(f"| {r.cell.split('_')[1]} | {SHORT[r.method]} | {fmt(r,'alarm_within_390')} | {fmt(r,'alarm_within_1950')} | {fmt(r,'up_alarm_within_1950')} | {fmt(r,'down_alarm_within_1950')} | {r.median_delay_from_monitoring_start:.0f} ({r.delay_iqr}) | {'' if pd.isna(mm) else f'{mm:.2f}'} |")
    print(f"\npaths per cell: {int(df.paths.min())}")

if __name__ == "__main__":
    groups = [("(a) Tail shape switch at constant variance", ["t30_to_t3", "t3_to_t30", "t30_to_t5", "t5_to_t30"]),
              ("(b) MA(1) bid-ask bounce and (d) skewed two-piece normal", ["ma1_-0.1", "ma1_-0.3", "ma1_-0.6", "ma1_0.3", "twopiece_0.6_1.4", "twopiece_0.4_1.6"]),
              ("(c) Prices on a tick grid", ["tick_zero20", "tick_zero40", "tick_zero60"]),
              ("(e) Missing data and (f) estimated intraday profile", ["missing_5pct", "missing_blocks", "seasonal_pipeline", "seasonal_amp1.0"]),
              ("(g) Forecast-relative pipeline with a misspecified GARCH(1,1) fit", ["fc_gjr", "fc_egarch", "fc_lnsv"])]
    import os
    null_tables(sys.argv[1], groups, sys.argv[3] if len(sys.argv) > 3 and os.path.exists(sys.argv[3]) else None)
    if len(sys.argv) > 2: change_tables(sys.argv[2])
    if len(sys.argv) > 4 and os.path.exists(sys.argv[4]):
        df = pd.read_csv(sys.argv[4])
        print("\n(k) Variance x2 or x3 from monitored step 50, carried only by rare jumps (sd 20) that the truncation clips.\n")
        print("| Cell | Method | True variance after | Clipped mean after / before | Alarm within 390 | Alarm within 1950 | Up within 1950 |"); print("|---|---|---|---|---|---|---|")
        for _, r in df.iterrows():
            tv = "" if pd.isna(r.get("true_var_after", np.nan)) else f"{r.true_var_after:.0f}"
            cm = "" if pd.isna(r.get("clipped_mean_ratio_after_over_before", np.nan)) else f"{r.clipped_mean_ratio_after_over_before:.3f}"
            print(f"| {r.cell} | {SHORT[r.method]} | {tv} | {cm} | {r.alarm_within_390:.3f} {r.alarm_within_390_ci} | {r.alarm_within_1950:.3f} {r.alarm_within_1950_ci} | {r.up_alarm_within_1950:.3f} {r.up_alarm_within_1950_ci} |")
        print(f"\npaths per cell: {int(df.paths.min())}")
