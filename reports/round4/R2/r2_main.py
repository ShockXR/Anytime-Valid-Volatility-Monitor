"""R2. Power at EQUAL measured false-alarm rates (SYNTHETIC). Extends reports/stress/b2_power_delay.py.

Usage: python r2_main.py <N_tune> <N_check> <N_cell> <out_prefix>
Seeds (numpy default_rng): truth sample 9_600_000 + law; tuning null paths 9_610_000 + 100 law + chunk;
fresh null paths 9_620_000 + 100 law + chunk; change cells 9_630_000 + 1000 law + 100 factor + chunk.
law index: gaussian 0, student_t3 1, jumps 2. factor index: 2 -> 0, 1.5 -> 1, 0.5 -> 2, 2/3 -> 3.
Calibration 1,000 steps, then monitoring. False alarm = any alarm within the first 390 (and 1,950) monitored
steps of a null path. Change at monitored step 400; detection = first alarm in [400, 400 + 1950); delay = alarm - 400 + 1.
"""
import sys, time, json
import numpy as np, pandas as pd
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/round4/R2")
from r2_lib import (A, PROVED, CAL, METHODS, E_DETECTORS, LAWS, draw, truth, all_stats, first_alarm,
                    oracle_threshold, wilson)

N_TUNE, N_CHECK, N_CELL, OUT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
CHUNK = 500
LAW_NAMES = list(LAWS)
FACTORS = [2.0, 1.5, 0.5, 2 / 3]
NU, H, DAY, FIVE = 400, 1950, 390, 1950
ROC_LEVELS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30]
rng_boot = np.random.default_rng(9_699_999)


def null_maxima(law, li, N, base):
    """Per-path maxima of every statistic over the first 390 and 1,950 monitored steps of null paths."""
    tau2, m0 = TRUTH[law]
    m390 = {m: [] for m in METHODS}
    m1950 = {m: [] for m in METHODS}
    for c in range(int(np.ceil(N / CHUNK))):
        n = min(CHUNK, N - c * CHUNK)
        data = draw(law, n, CAL + FIVE, np.random.default_rng(base + 100 * li + c))
        W = all_stats(data, tau2, m0)
        for m in METHODS:
            m390[m].append(W[m][:, :DAY].max(axis=1))
            m1950[m].append(W[m].max(axis=1))
    return {m: np.concatenate(v) for m, v in m390.items()}, {m: np.concatenate(v) for m, v in m1950.items()}


def boot_median_ci(d, B=400):
    if d.size < 5:
        return np.nan, np.nan
    meds = np.median(rng_boot.choice(d, (B, d.size), replace=True), axis=1)
    return float(np.quantile(meds, 0.025)), float(np.quantile(meds, 0.975))


t_start = time.time()
TRUTH = {law: truth(law, 9_600_000 + li) for li, law in enumerate(LAW_NAMES)}
print("truth (tau^2, m0):", TRUTH, flush=True)

# ---- 1. tuning on null paths of each law: thresholds -------------------------------------------------------------
TH = {}  # (law, method) -> list of (label, h)
tune_max, tune_max5 = {}, {}
for li, law in enumerate(LAW_NAMES):
    t0 = time.time()
    m390, m1950 = null_maxima(law, li, N_TUNE, 9_610_000)
    tune_max[law], tune_max5[law] = m390, m1950
    print(f"tuning {law} done in {time.time() - t0:.0f}s", flush=True)
for law in LAW_NAMES:
    for m in METHODS:
        th = [("oracle", oracle_threshold(tune_max[law][m], 0.05)),
              ("gauss_tuned", oracle_threshold(tune_max["gaussian"][m], 0.05)),
              ("oracle_5d", oracle_threshold(tune_max5[law][m], 0.25))]  # 25% within 1,950 steps (the 5-day bound)
        if m in E_DETECTORS or m == "Gaussian LR monitor":
            th.append(("nominal_2A", PROVED))
        th += [(f"roc_{lv}", oracle_threshold(tune_max[law][m], lv)) for lv in ROC_LEVELS]
        TH[(law, m)] = th
json.dump({f"{k[0]}|{k[1]}": v for k, v in TH.items()}, open(OUT + "_thresholds.json", "w"), indent=1)
np.savez(OUT + "_tune_max.npz", **{f"{law}|{m}|390": tune_max[law][m] for law in LAW_NAMES for m in METHODS},
         **{f"{law}|{m}|1950": tune_max5[law][m] for law in LAW_NAMES for m in METHODS})

rows = []
# ---- 2. fresh null paths: realised false-alarm rates for every threshold -----------------------------------------
for li, law in enumerate(LAW_NAMES):
    t0 = time.time()
    m390, m1950 = null_maxima(law, li, N_CHECK, 9_620_000)
    for m in METHODS:
        for label, h in TH[(law, m)]:
            for horizon, mx in ((DAY, m390[m]), (FIVE, m1950[m])):
                k = int((mx >= h).sum())
                lo, hi = wilson(k, N_CHECK)
                rows.append(dict(kind="nullcheck", law=law, method=m, label=label, threshold=h, horizon=horizon,
                                 factor=1.0, paths=N_CHECK, k=k, rate=k / N_CHECK, ci_low=lo, ci_high=hi))
    print(f"check {law} done in {time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(OUT + "_rows.csv", index=False)

# ---- 3. change cells: detection and delay for every threshold ----------------------------------------------------
for li, law in enumerate(LAW_NAMES):
    tau2, m0 = TRUTH[law]
    for fi, f in enumerate(FACTORS):
        t0 = time.time()
        fa = {(m, lab): [] for m in METHODS for lab, _ in TH[(law, m)]}
        for c in range(int(np.ceil(N_CELL / CHUNK))):
            n = min(CHUNK, N_CELL - c * CHUNK)
            rng = np.random.default_rng(9_630_000 + 1000 * li + 100 * fi + c)
            data = draw(law, n, CAL + NU + H, rng, nu=CAL + NU, f=f)
            W = all_stats(data, tau2, m0)
            for m in METHODS:
                cm = np.maximum.accumulate(W[m], axis=1)
                for lab, h in TH[(law, m)]:
                    fa[(m, lab)].append(first_alarm(cm, h))
        for m in METHODS:
            for lab, h in TH[(law, m)]:
                al = np.concatenate(fa[(m, lab)])
                early = (al >= 0) & (al < NU)
                det = (al >= NU) & (al < NU + H)
                d = (al[det] - NU + 1).astype(float)
                k = int(det.sum())
                lo, hi = wilson(k, N_CELL)
                blo, bhi = boot_median_ci(d)
                rows.append(dict(kind="change", law=law, method=m, label=lab, threshold=h, horizon=H, factor=round(f, 4),
                                 paths=N_CELL, k=k, rate=k / N_CELL, ci_low=lo, ci_high=hi, early=float(early.mean()),
                                 median_delay=float(np.median(d)) if d.size else np.nan,
                                 q25=float(np.quantile(d, .25)) if d.size else np.nan,
                                 q75=float(np.quantile(d, .75)) if d.size else np.nan,
                                 mean_delay=float(d.mean()) if d.size else np.nan, med_ci_low=blo, med_ci_high=bhi))
        print(f"change {law} f={f:.3f} done in {time.time() - t0:.0f}s", flush=True)
        pd.DataFrame(rows).to_csv(OUT + "_rows.csv", index=False)
print(f"EXIT total {time.time() - t_start:.0f}s", flush=True)
