"""Round 9 checks of running.py (public monitors). [SYNTHETIC]

1. Every rank bet of both families integrates to 1 over a uniform p-value and is nonnegative.
2. ConformalMonitor (route B), i.i.d. Student t3 returns: P(false alarm within m) <= m / A (exact).
3. RunningTauMonitor (route A, grid band), i.i.d. Gaussian: P(false alarm within m) <= delta + m / A.
4. AwakeRunningTauMonitor, i.i.d. Gaussian: P(false alarm within m) <= delta + m / A.
Usage: python r9_running.py [paths]
"""
import math
import sys

import numpy as np
from scipy import integrate

from avmonitor.running import AwakeRunningTauMonitor, ConformalMonitor, RunningTauMonitor, bet_family


def integrals():
    worst_int = 0.0
    worst_min = np.inf
    for fam in ("lr", "ind"):
        up, dn = bet_family(fam)
        for f in up + dn:
            pts = [f.u0] if hasattr(f, "u0") else None
            val = integrate.quad(lambda u: float(f(u)), 0.0, 1.0, points=pts, limit=400)[0]
            worst_int = max(worst_int, abs(val - 1.0))
            worst_min = min(worst_min, float(np.min(f(np.linspace(0, 1, 200001)))))
    print(f"1. rank bets: largest |integral - 1| = {worst_int:.2e}; smallest value on [0, 1] = {worst_min:.4f}", flush=True)


def rate(name, make, warm, m, A, P, bound, law="gauss"):
    rng = np.random.default_rng(12)
    k = 0
    for i in range(P):
        mon = make()
        r = rng.standard_normal(warm + m) if law == "gauss" else rng.standard_t(3, warm + m)
        for v in r:
            if mon.update(v) is not None:
                k += 1
                break
    p = k / P
    se = math.sqrt(max(p * (1 - p), 1e-300) / P)
    print(f"[SYNTHETIC] {name}: rate {p:.4%} (se {se:.4%}) over {P} paths, m = {m}, A = {A:g}, bound {bound:.2%}",
          flush=True)


if __name__ == "__main__":
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    integrals()
    rate("2. ConformalMonitor, t3, cal_len 200", lambda: ConformalMonitor(threshold=7800, cal_len=200, seed=1),
         200, 390, 7800, P, 390 / 7800, law="t3")
    rate("3. RunningTauMonitor, Gaussian, cal_len 1440", lambda: RunningTauMonitor(threshold=28800, cal_len=1440),
         1440, 1440, 28800, P // 2, 0.01 + 1440 / 28800)
    rate("4. AwakeRunningTauMonitor, Gaussian, n_pre 120", lambda: AwakeRunningTauMonitor(threshold=28800),
         120, 1440, 28800, P // 2, 0.01 + 1440 / 28800)
