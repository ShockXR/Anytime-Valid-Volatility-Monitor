"""R6 step 4: markdown tables from r6_null_rows.csv and r6_delay.csv -> r6_tables.md. Usage: python r6_tables.py"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
null = pd.read_csv(HERE / "r6_null_rows.csv")
delay = pd.read_csv(HERE / "r6_delay.csv")
LABEL = {("plugin", 0.25): "plug-in, band 0.25 (README)", ("plugin", 0.1): "plug-in, band 0.1",
         ("plugin", 0.05): "plug-in, band 0.05", ("plugin", 0.0): "plug-in, band 0",
         ("cs", 0.0): "cs baseline, band 0", ("cs", 0.1): "cs baseline, band 0.1"}
pct = lambda x: f"{100 * x:.2f}%"
out = []


def null_table(ds):
    d = null[null.dataset == ds]
    out.append(f"\n**{ds}** ({int(d.streams.iloc[0]):,} streams per setting)\n")
    out.append("| Setting | A | Alarms in day 1 | Rate day 1 [Wilson 95%] | Bound | Alarms in 5 days | Rate 5 days [Wilson 95%] | Bound | Verdict |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in d.iterrows():
        v1 = "above" if r.one_day_vs_bound == "above" else "holds"
        v5 = "above" if r.five_days_vs_bound == "above" else "holds"
        verdict = "holds" if v1 == v5 == "holds" else f"day 1 {v1}, 5 days {v5}"
        out.append(f"| {LABEL[(r.monitor, r.band)]} | {int(r.A):,} | {int(r.one_day_alarms)} | {pct(r.one_day)} [{pct(r.one_day_lo)}, {pct(r.one_day_hi)}] | {pct(r.one_day_bound)} | "
                   f"{int(r.five_days_alarms)} | {pct(r.five_days)} [{pct(r.five_days_lo)}, {pct(r.five_days_hi)}] | {pct(r.five_days_bound)} | {verdict} |")


out.append("## Null tables [REAL marginal, constructed streams]\n")
out.append("How to read this: each row is one monitor setting run on permuted real calm-day returns, where no alarm is correct. "
           "'Rate day 1' is the share of streams that alarmed within the first monitored day; 'Rate 5 days' within the first five. "
           "A rate is 'above' its bound when the lower end of its Wilson interval is above the bound (E19's rule), and 'holds' otherwise.")
for ds in ["five coins", "BTC", "ETH", "XRP", "SOL", "DOGE", "Nifty 50"]:
    if (null.dataset == ds).any():
        null_table(ds)

out.append("\n## Delay tables [SYNTHETIC]\n")
out.append("How to read this: each row is one setting on simulated paths with 1,000 calibration steps and a variance change at monitored step 400. "
           "'Detected' is the share of paths whose first alarm came after the change and within 1,950 steps (five days of 390). "
           "'Early' is the share whose first alarm came before the change (a false alarm). "
           "The median delay is the middle number of steps from the change to the alarm, among detected paths; the IQR is the middle half.")
for law in ["gaussian", "student_t3"]:
    for f in [2.0, 0.5]:
        d = delay[(delay.law == law) & (delay.var_factor == f)]
        out.append(f"\n**{law}, variance x{f:g}** ({int(d.paths.iloc[0]):,} paths per setting)\n")
        out.append("| Setting | A | Detected [Wilson 95%] | Early | Median delay | IQR | Mean delay |")
        out.append("|---|---|---|---|---|---|---|")
        for _, r in d.iterrows():
            out.append(f"| {LABEL[(r.monitor, r.band)]} | {int(r.A):,} | {pct(r.detect_rate)} [{pct(r.det_ci_low)}, {pct(r.det_ci_high)}] | {pct(r.early_false_alarm)} | "
                       f"{r.median_delay:.0f} | {r.q25_delay:.0f} to {r.q75_delay:.0f} | {r.mean_delay:.0f} |")
(HERE / "r6_tables.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
