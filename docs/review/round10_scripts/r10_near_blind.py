"""Review round 10 (E24): Lemma 1 and Proposition 2 of docs/theory_blind_guard.md, recomputed independently.

Closed form for the largest increase-side average d steps after a (re)start when every factor is at its maximum
g_k = 1 + lam_k (rho - 1):  a_d(rho) = mean_k g_k (g_k^d - 1) / (g_k - 1)  (a_d = d when rho = 1).
d_min(rho) = min{d : a_d(rho) >= 2A}. Nothing from avmonitor.guard is used for the computation; the package
function near_blind_delay_floor is called only to compare.

Run: nice -n 10 python docs/review/round10_scripts/r10_near_blind.py   [PROVED arithmetic, seconds]
"""

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

LAM = np.array([0.02, 0.05, 0.1, 0.2, 0.35, 0.5])  # DEFAULT_UP_GRID
A, M = 36_000.0, 1440


def a_d(rho: float, d: int) -> float:
    g = 1.0 + LAM * (rho - 1.0)
    if rho == 1.0:
        return float(d)
    with np.errstate(over="ignore"):
        return float(np.mean(g * (np.power(g, d) - 1.0) / (g - 1.0)))


def d_min(rho: float, A: float = A, cap: int = 10**6) -> float:
    if rho < 1:
        return math.inf
    lo, hi = 1, cap
    if a_d(rho, hi) < 2 * A:
        return math.inf
    while lo < hi:  # a_d is increasing in d
        mid = (lo + hi) // 2
        if a_d(rho, mid) >= 2 * A:
            hi = mid
        else:
            lo = mid + 1
    return float(lo)


def main():
    from avmonitor.guard import near_blind_delay_floor

    print("rho, d_min (closed form), d_min (package), max log growth per step log(1 + 0.5 (rho - 1))")
    for rho in (1.0, 1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0):
        print(f"{rho:6.3f} {d_min(rho):8.0f} {near_blind_delay_floor(rho, A, cap=10**6):8.0f} "
              f"{math.log(1 + 0.5 * (rho - 1)):.6f}")
    # threshold rho* with a_M(rho*) = 2A: no increase alarm within one restarted period of M steps iff rho < rho*
    lo, hi = 1.0, 1.1
    for _ in range(100):
        mid = (lo + hi) / 2
        if a_d(mid, M) < 2 * A:
            lo = mid
        else:
            hi = mid
    print(f"rho* (a_1440(rho*) = 72,000): {lo:.10f}")
    for rho in (1.0105, 1.01050, 1.0106):
        print(f"a_1440({rho}) = {a_d(rho, M):.1f}  (2A = 72,000; alarm possible within a day: {a_d(rho, M) >= 2 * A})")
    # Lemma 1 at rho = 1: a_d = d, so the first possible alarm is at step 2A exactly (not 'within 2A steps' impossible)
    print(f"rho = 1: a_(2A - 1) = {a_d(1.0, int(2 * A) - 1):.0f} < 2A, a_(2A) = {a_d(1.0, int(2 * A)):.0f} = 2A -> "
          "alarm possible at step 2A, impossible within 2A - 1 steps")
    # with a carry kappa at each restart the never-alarm condition of Lemma 1 is m + kappa < 2A
    print("P1's carry is 0 (PeriodicVolatilityMonitor default), so m = 1440 < 2A = 72,000 suffices.")
    # the 390-step setting of e24_near_blind.csv
    A2 = 390 / 0.04
    print("390-step setting, A = 9750: d_min at 1.01, 1.25:", d_min(1.01, A2), d_min(1.25, A2))


if __name__ == "__main__":
    main()
