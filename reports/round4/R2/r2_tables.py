"""R2 tables (markdown) from r2_rows.csv. Usage: python r2_tables.py <rows.csv> > tables.md"""
import sys
import numpy as np, pandas as pd

d = pd.read_csv(sys.argv[1], keep_default_na=False, na_values=[""])
SHORT = {"README: plug-in, band 0.25": "README (plug-in, band 0.25)", "plug-in, band 0": "Plug-in, band 0",
         "confidence-sequence baseline": "Confidence-sequence baseline", "known baseline (theorem)": "Known baseline (theorem)",
         "Gaussian LR monitor": "Gaussian LR", "Page CUSUM up": "Page CUSUM (up only)", "repeated 60-step z-test": "Repeated 60-step z-test"}
METHODS = list(SHORT)
LAWS = ["gaussian", "student_t3", "jumps_p0.5pct_sd8"]
LAWN = {"gaussian": "Gaussian", "student_t3": "Student t3", "jumps_p0.5pct_sd8": "Gaussian + jumps"}
FACT = [2.0, 1.5, 0.5, 0.6667]
FN = {2.0: "x2", 1.5: "x1.5", 0.5: "x0.5", 0.6667: "x2/3"}
nul = d[d.kind == "nullcheck"]
chg = d[d.kind == "change"]


def fa(law, m, lab, hor=390):
    r = nul[(nul.law == law) & (nul.method == m) & (nul.label == lab) & (nul.horizon == hor)]
    if r.empty:
        return "n/a"
    r = r.iloc[0]
    return f"{100 * r.rate:.1f}% [{100 * r.ci_low:.1f}, {100 * r.ci_high:.1f}]"


def cell(law, m, lab, f, det=True):
    if m == "Page CUSUM up" and f < 1:
        return "n/a"
    r = chg[(chg.law == law) & (chg.method == m) & (chg.label == lab) & (np.isclose(chg.factor, f))]
    if r.empty:
        return "n/a"
    r = r.iloc[0]
    if not np.isfinite(r.median_delay) or r.rate < 0.05:
        return f"{100 * r.rate:.1f}% detected"
    s = f"{r.median_delay:.0f} ({r.q25:.0f}-{r.q75:.0f})"
    return s + (f", {100 * r.rate:.0f}%" if det else "")


def th(law, m, lab):
    r = nul[(nul.law == law) & (nul.method == m) & (nul.label == lab) & (nul.horizon == 390)]
    return r.iloc[0].threshold if not r.empty else np.nan


N_CHECK = int(nul.paths.iloc[0]); N_CELL = int(chg.paths.iloc[0])
print(f"### T1. Oracle thresholds and realised false-alarm rates on fresh null paths ({N_CHECK:,} paths per law)\n")
print("| Law | Method | Oracle threshold (alarm when statistic >= this) | Proved/nominal threshold | False alarms within 390 steps [95% CI] | Within 1,950 steps |")
print("|---|---|---|---|---|---|")
for law in LAWS:
    for m in METHODS:
        h = th(law, m, "oracle")
        nom = "15,600 (= 2A)" if m in METHODS[:5] else "none"
        hs = f"{h:.4g}" if h < 1e6 else f"{h:.2e}"
        print(f"| {LAWN[law]} | {SHORT[m]} | {hs} | {nom} | {fa(law, m, 'oracle')} | {fa(law, m, 'oracle', 1950)} |")
print()
for law in LAWS:
    print(f"### T2-{LAWN[law]}. Oracle tuning (5% per 390 steps on this law): median delay (IQR), detected within 1,950 steps ({N_CELL:,} paths per cell)\n")
    print("| Method | " + " | ".join(f"variance {FN[f]}" for f in FACT) + " |")
    print("|---|" + "---|" * len(FACT))
    for m in METHODS:
        print(f"| {SHORT[m]} | " + " | ".join(cell(law, m, "oracle", f) for f in FACT) + " |")
    print()

print("### T3. Where the README setting's slowness comes from (median delay in steps)\n")
print("| Law | Variance | README at 2A | README at oracle | Band 0 at 2A | Band 0 at oracle | Known baseline at oracle | Best non-e competitor at oracle |")
print("|---|---|---|---|---|---|---|---|")
for law in LAWS:
    for f in FACT:
        def md(m, lab):
            r = chg[(chg.law == law) & (chg.method == m) & (chg.label == lab) & (np.isclose(chg.factor, f))]
            return r.iloc[0].median_delay if not r.empty else np.nan
        comp = [(md(m, "oracle"), SHORT[m]) for m in METHODS[4:] if not (m == "Page CUSUM up" and f < 1)]
        comp = [c for c in comp if np.isfinite(c[0])]
        b = min(comp)
        print(f"| {LAWN[law]} | {FN[f]} | {md(METHODS[0], 'nominal_2A'):.0f} | {md(METHODS[0], 'oracle'):.0f} | {md(METHODS[1], 'nominal_2A'):.0f} | "
              f"{md(METHODS[1], 'oracle'):.0f} | {md(METHODS[3], 'oracle'):.0f} | {b[0]:.0f} ({b[1]}) |")
print()

print("### T4. Practical tuning: competitors tuned on Gaussian data, e-detectors at the proved 2A\n")
print("| Law | Method | Threshold used | False alarms within 390 steps [95% CI] | " + " | ".join(f"variance {FN[f]}" for f in FACT) + " |")
print("|---|---|---|---|" + "---|" * len(FACT))
for law in LAWS[1:]:
    for m in METHODS:
        lab = "nominal_2A" if m in METHODS[:4] else "gauss_tuned"
        print(f"| {LAWN[law]} | {SHORT[m]} | {'2A (proved)' if lab == 'nominal_2A' else 'tuned on Gaussian'} | {fa(law, m, lab)} | "
              + " | ".join(cell(law, m, lab, f) for f in FACT) + " |")
print()
print("### T5. Extra: e-detectors and GLR also tuned on Gaussian data (5% per 390 steps there), used on the other laws\n")
print("| Law | Method | False alarms within 390 steps [95% CI] | " + " | ".join(f"variance {FN[f]}" for f in FACT) + " |")
print("|---|---|---|" + "---|" * len(FACT))
for law in LAWS[1:]:
    for m in METHODS[:4]:
        print(f"| {LAWN[law]} | {SHORT[m]} | {fa(law, m, 'gauss_tuned')} | " + " | ".join(cell(law, m, "gauss_tuned", f) for f in FACT) + " |")
print()
print("### T6. Nominal thresholds for reference (as in B2): e-detectors and GLR at 2A, on each law\n")
print("| Law | Method | False alarms within 390 steps [95% CI] | " + " | ".join(f"variance {FN[f]}" for f in FACT) + " |")
print("|---|---|---|" + "---|" * len(FACT))
for law in LAWS:
    for m in METHODS[:5]:
        print(f"| {LAWN[law]} | {SHORT[m]} | {fa(law, m, 'nominal_2A')} | " + " | ".join(cell(law, m, "nominal_2A", f) for f in FACT) + " |")
print()
print("### T7. Matched over five days instead of one: threshold that gives 25% false alarms within 1,950 steps on this law\n")
print("| Law | Method | False alarms within 390 / 1,950 steps | " + " | ".join(f"variance {FN[f]}" for f in FACT) + " |")
print("|---|---|---|" + "---|" * len(FACT))
for law in LAWS:
    for m in METHODS:
        r1 = nul[(nul.law == law) & (nul.method == m) & (nul.label == "oracle_5d") & (nul.horizon == 390)].iloc[0]
        r5 = nul[(nul.law == law) & (nul.method == m) & (nul.label == "oracle_5d") & (nul.horizon == 1950)].iloc[0]
        print(f"| {LAWN[law]} | {SHORT[m]} | {100 * r1.rate:.1f}% / {100 * r5.rate:.1f}% | " + " | ".join(cell(law, m, "oracle_5d", f, det=False) for f in FACT) + " |")
