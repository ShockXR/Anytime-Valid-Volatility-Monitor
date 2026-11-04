"""E4. Detecting volatility changes: false alarms and detection delay.

Research questions
  Q4a. Under no change, do the e-detectors keep their guaranteed average run
       length (ARL >= A) when returns have jumps, heavy tails, slow stochastic
       volatility or price noise? Do the benchmarks?
  Q4b. How fast do they detect increases and decreases of different sizes?
  Q4c. Does the quarticity-based adaptive bet help over a fixed grid of bets?
  Q4d. What does truncation buy: compare with the same detector on raw r^2.

Setup (synthetic, per-step units, sigma0 = 1): each path has a 2000-step
calibration window followed by monitoring. The baseline (scale, tau = 3 x
scaled MAD, m0) is estimated from the calibration window of the same path,
which has the same jump, tail and noise features as the monitored data.
E-detectors use A = 7800, so ARL >= 7800 and P(false alarm within a
390-step day) <= 5%. Page's CUSUM is calibrated by simulation to ARL about
7800 on clean Gaussian data. The repeated z-test uses a 60-step window
at level 0.05 per look.
"""

from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import detectors as D
from avmonitor.simulate import change_returns

OUT = Path(__file__).resolve().parents[1] / "results"
A = 7800.0  # P(false alarm within one 390-step day) <= 390 / A = 5%
CAL = 2000

SCENARIOS = {
    "gaussian": dict(),
    "jumps": dict(p_jump=0.005, jump_scale=8.0),
    "student t5": dict(heavy_tail_df=5.0),
    "slow SV": dict(sv_xi=0.15, sv_halflife=390.0),
    "price noise": dict(noise_sd=0.5),
}


def calibrate_page(drift=0.25, reps=300, n=40000, seed=99):
    rng = np.random.default_rng(seed)
    r = change_returns(reps, n, rng=rng).r
    lo, hi = 1.0, 60.0
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        al = D.page_cusum(r**2, 1.0, drift, mid).alarm
        arl = np.where(al >= 0, al + 1, n).mean()
        lo, hi = (mid, hi) if arl < A else (lo, mid)
    return hi


def run_detectors(r, cal, page_h, nu):
    """Return alarm indices (relative to monitoring start) for every detector."""
    reps = r.shape[0]
    out = {}
    bases = [D.calibrate(cal[i], c=3.0) for i in range(reps)]
    tau = np.array([b.tau for b in bases])[:, None]
    m0 = np.array([b.m0 for b in bases])[:, None]
    scale = np.array([b.scale for b in bases])[:, None]
    # Normalise each path by its own baseline so that one m0 = 1 and tau^2 apply.
    y = np.minimum(r**2, tau**2) / m0
    ub = float(np.max(tau**2 / m0))
    # per-path bounds differ slightly; use each path's own bound via max (conservative for bets)
    out["e-det up (grid)"] = D.betting_edetector(y, 1.0, "up", threshold=A).alarm
    out["e-det up (grid + quarticity bet)"] = D.betting_edetector(y, 1.0, "up", threshold=A, adaptive_h=60).alarm
    out["e-det up CUSUM-type"] = D.betting_edetector(y, 1.0, "up", threshold=A, kind="CUSUM").alarm
    out["e-det down (grid)"] = D.betting_edetector(y, 1.0, "down", upper=ub, threshold=A).alarm
    out["e-det down (grid + quarticity bet)"] = D.betting_edetector(y, 1.0, "down", upper=ub, threshold=A, adaptive_h=60).alarm
    raw_m0 = np.array([np.mean(cal[i] ** 2) for i in range(reps)])[:, None]
    out["e-det up on raw r^2 (no truncation)"] = D.betting_edetector(r**2 / raw_m0, 1.0, "up", threshold=A).alarm
    out["Gaussian LR e-det up"] = D.gaussian_lr_edetector(r / scale, 1.0, threshold=A).alarm
    out["Gaussian LR e-det down"] = D.gaussian_lr_edetector(r / scale, 1.0, ratios=(0.5, 0.7, 0.8), threshold=A).alarm
    out["Page CUSUM up (calibrated)"] = D.page_cusum(r**2 / scale**2, 1.0, 0.25, page_h).alarm
    out["repeated 60-step z-test"] = D.repeated_window_ztest(r / scale, 1.0, 2.0, 60).alarm
    return out


def null_job(args):
    name, kw, page_h, seed = args
    rng = np.random.default_rng(seed)
    n = 20000
    p = change_returns(300, CAL + n, rng=rng, **kw)
    al = run_detectors(p.r[:, CAL:], p.r[:, :CAL], page_h, None)
    rows = []
    for det, a in al.items():
        rl = np.where(a >= 0, a + 1, n)
        rows.append(dict(scenario=name, detector=det, false_alarm_within_20000=float((a >= 0).mean()),
                         false_alarm_within_390=float(((a >= 0) & (a < 390)).mean()),
                         mean_run_length_censored=float(rl.mean())))
    return rows


def change_job(args):
    name, kw, ratio, page_h, seed = args
    rng = np.random.default_rng(seed)
    nu, n = 500, 3500
    p = change_returns(300, CAL + n, nu=CAL + nu, ratio=ratio, rng=rng, **kw)
    al = run_detectors(p.r[:, CAL:], p.r[:, :CAL], page_h, nu)
    rows = []
    for det, a in al.items():
        s = D.summarize_alarms(a, nu, n)
        rows.append(dict(scenario=name, ratio=ratio, detector=det, **s))
    return rows


def main():
    page_h = calibrate_page()
    print("Page CUSUM threshold for ARL ~ A:", round(page_h, 3))
    jobs = [(nm, kw, page_h, 100 + i) for i, (nm, kw) in enumerate(SCENARIOS.items())]
    with Pool(4) as pool:
        null_rows = sum(pool.map(null_job, jobs), [])
    nd = pd.DataFrame(null_rows)
    nd.to_csv(OUT / "e04a_null.csv", index=False)
    print(nd.pivot(index="detector", columns="scenario", values="mean_run_length_censored").round(0).to_string())
    print(nd.pivot(index="detector", columns="scenario", values="false_alarm_within_390").round(3).to_string())

    cjobs = []
    for i, nm in enumerate(("gaussian", "jumps", "student t5")):
        for j, ratio in enumerate((1.25, 1.5, 2.0, 0.8, 0.67, 0.5)):
            cjobs.append((nm, SCENARIOS[nm], ratio, page_h, 1000 + 10 * i + j))
    with Pool(4) as pool:
        ch_rows = sum(pool.map(change_job, cjobs), [])
    cd = pd.DataFrame(ch_rows)
    cd.to_csv(OUT / "e04b_change.csv", index=False)
    for nm in ("gaussian", "jumps", "student t5"):
        sub = cd[cd.scenario == nm]
        print(f"\n[{nm}] median delay (detected share) by sd ratio")
        tab = sub.assign(cell=sub.median_delay.round(0).astype(str) + " (" + sub.detected.round(2).astype(str) + ", FA " + sub.false_alarm_before_change.round(2).astype(str) + ")")
        print(tab.pivot(index="detector", columns="ratio", values="cell").to_string())
    print(nd.pivot(index="detector", columns="scenario", values="false_alarm_within_20000").round(3).to_string())


if __name__ == "__main__":
    main()
