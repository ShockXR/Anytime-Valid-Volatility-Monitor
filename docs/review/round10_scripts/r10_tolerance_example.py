"""Review round 10 (E24): Remark 3 of docs/theory_blind_guard.md, first example of the second bullet.

Claim checked: "If volatility rises in the middle of window 1's CS part and stays, the tolerance fails relative to
window 1 but holds relative to window 2 (the guarded monitor is covered, the unguarded is not)". The guarded monitor
is covered by Theorem 2 only if the tolerance holds relative to the KEPT window almost surely, so only if window 1 is
discarded (near-blind) with probability one. This script computes, for P1's window (1440 minutes, tau from the first
360, delta_1 = 0.005, b = 0.25), the share of first windows that the recommended guard (eta = 0.25) would keep when
the standard deviation of i.i.d. Gaussian returns rises by a factor k at minute s of the window. On every kept first
window the tolerance relative to the kept window fails (the rise exceeds the band), so Theorem 2 does not cover the
guarded monitor there either; part 4 of Theorem 6.4 charges P(window 1 kept).

Then a direct run of the guarded P1 (GE) and the unguarded P1 (U) on such paths: share with an alarm on the first
monitored day.

Run: nice -n 10 python docs/review/round10_scripts/r10_tolerance_example.py   [SYNTHETIC; about 2 minutes]
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from avmonitor.confseq import cs_baseline  # noqa: E402
from avmonitor.guard import GuardedPeriodicVolatilityMonitor  # noqa: E402
from avmonitor.periodic import PeriodicVolatilityMonitor  # noqa: E402

CAL, NT, B, ETA = 1440, 360, 0.25, 0.25


def kept_share(k: float, s: int, reps: int, rng) -> tuple[float, float]:
    sd = np.where(np.arange(CAL) < s, 1.0, k)
    r = rng.standard_normal((reps, CAL)) * sd
    cb = cs_baseline(r, c=3.0, delta_upper=0.0025, delta_lower=0.0025, tau_frac=0.25)
    rho = cb.tau**2 / ((1 + B) * cb.upper)
    return float(np.mean(rho > 1 + ETA)), float(np.median(rho))


def first_day_alarm(k: float, s: int, paths: int, rng) -> tuple[float, float]:
    base = dict(threshold=1440 / 0.04, c=3.0, cal_len=CAL, null_band=B, baseline="cs", delta=0.01, period=1440)
    hits_u = hits_g = 0
    for _ in range(paths):
        n = 5 * CAL
        sd = np.where(np.arange(n) < s, 1.0, k)
        u = rng.standard_normal(n) * sd
        a = PeriodicVolatilityMonitor(**base)
        g = GuardedPeriodicVolatilityMonitor(**base, guard_tries=3, guard_eta=ETA, guard_fallback="rank", rank_seed=1)
        a.run(u)
        g.run(u)
        for mon, name in ((a, "u"), (g, "g")):
            s0 = mon.segments[0]["start"]
            hit = any(e.monitored_since == s0 and e.index - s0 < 1440 for e in mon.events)
            if name == "u":
                hits_u += hit
            else:
                hits_g += hit
    return hits_u / paths, hits_g / paths


def main():
    rng = np.random.default_rng(31)
    print("rise by k at minute s of window 1 (tau part = minutes 0 to 359): share of first windows KEPT by the "
          "guard (rho > 1.25), median rho")
    for s in (360, 540, 900, 1260):
        row = []
        for k in (2, 4, 8, 16, 32):
            ks, med = kept_share(k, s, 4000, rng)
            row.append(f"k={k}: kept {ks:.3f} (median rho {med:.2f})")
        print(f"s = {s:4d}: " + "; ".join(row))
    for k in (2.0, 4.0):
        pu, pg = first_day_alarm(k, 900, 100, rng)
        print(f"rise x{k:g} at minute 900 (middle of the CS part), 100 paths: first-day alarm share, unguarded P1 "
              f"{pu:.2f}, guarded P1 (GE) {pg:.2f}")


if __name__ == "__main__":
    main()
