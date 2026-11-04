"""Review round 10 (E24): Proposition 3 of docs/theory_blind_guard.md in numbers, recomputed independently.

Under the construction of Proposition 3 (b = 0.25, eta = 0, so q = 0.8), the confidence-sequence part of every window
is i.i.d. Bernoulli(q) in units of tau^2, and a window is kept exactly when U_T < q, in which case its upper bound
fails. p(delta) = P(U_T < q) for one window's upper sequence at level delta / 2, n_plan = n_cs.

The upper sequence is reimplemented here from docs/math.md 6.2 (remark 2) without importing avmonitor.confseq; it is
checked against avmonitor.confseq.eb_confseq on a small batch first. Then:
  * p(delta) for P1's window (n_cs = 1080) at delta = 0.01, 0.005, 0.0025, and for the illustration (n_cs = 150) at
    delta = 0.3, 0.15, 0.075;
  * the kept window's failure probability with K windows: 1 - prod_k (1 - p(delta_k));
  * K*, the number of windows with the same delta each at which the naive guard first exceeds the nominal delta,
    and the uncapped guard (E19's proposal), for which P(kept window fails) = 1 under the construction.

Run: nice -n 10 python docs/review/round10_scripts/r10_selection.py [--reps=400000]   [SYNTHETIC; about 3 minutes]
"""

import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

Q = 0.8
C_CAP = 0.75


def h(x):
    x = np.asarray(x, float)
    small = x < 1e-3
    xs = np.where(small, 0.5, x)
    return np.where(small, 0.5 + x / 3 + x * x / 4 + x**3 / 5, (-np.log1p(-xs) - xs) / (xs * xs))


def upper_final(z: np.ndarray, delta_u: float, n_plan: float) -> np.ndarray:
    """U_n (clipped at 1) of the predictable plug-in empirical Bernstein upper sequence, rows = replications."""
    reps, n = z.shape
    log_u = math.log(1.0 / delta_u)
    zhat = np.full(reps, 0.5)
    s2 = np.full(reps, 0.25)
    sum_z = np.zeros(reps)
    sum_dev2 = np.zeros(reps)
    sl = np.zeros(reps)
    slz = np.zeros(reps)
    spe = np.zeros(reps)
    for t in range(1, n + 1):
        x = z[:, t - 1]
        a = 1.0 - zhat
        lam = np.minimum(np.sqrt(2.0 * log_u / (n_plan * s2)), C_CAP / a)
        sl += lam
        slz += lam * x
        spe += lam * lam * h(lam * a) * (x - zhat) ** 2
        sum_z += x
        zhat = (0.5 + sum_z) / (t + 1.0)
        sum_dev2 += (x - zhat) ** 2
        s2 = (0.25 + sum_dev2) / (t + 1.0)
    return np.minimum((slz + log_u + spe) / sl, 1.0)


def p_fail(n_cs: int, delta_window: float, reps: int, rng, batch: int = 20_000) -> tuple[float, int]:
    hits = done = 0
    while done < reps:
        b = min(batch, reps - done)
        z = (rng.random((b, n_cs)) < Q).astype(float)
        hits += int((upper_final(z, delta_window / 2, n_cs) < Q).sum())
        done += b
    return hits / reps, hits


def main():
    opts = dict(a.lstrip("-").split("=", 1) for a in sys.argv[1:] if "=" in a)
    reps = int(opts.get("reps", 400_000))
    from avmonitor.confseq import eb_confseq

    rng = np.random.default_rng(20261103)
    z = (rng.random((300, 1080)) < Q).astype(float)
    mine = upper_final(z, 0.005, 1080)
    ref = eb_confseq(z, 1.0, 0.005, 0.005, n_plan=1080).upper[:, -1]
    print(f"check against eb_confseq: max |difference| = {np.max(np.abs(mine - ref)):.2e}")

    out = {}
    for n_cs, deltas, r in ((1080, (0.01, 0.005, 0.0025), reps), (150, (0.3, 0.15, 0.075), 2 * reps)):
        for d in deltas:
            t0 = time.time()
            p, k = p_fail(n_cs, d, r, rng)
            se = math.sqrt(p * (1 - p) / r)
            out[(n_cs, d)] = p
            print(f"n_cs = {n_cs:5d}, delta_window = {d:7.4f}: p = {p:.4g} ({k} of {r}, s.e. {se:.2g}); "
                  f"delta/2 / p = {d / 2 / p:.1f}  [{time.time() - t0:.0f} s]")

    for n_cs, d, half in ((1080, 0.01, (0.005, 0.0025)), (150, 0.3, (0.15, 0.075))):
        p = out[(n_cs, d)]
        same = [1 - (1 - p) ** K for K in (1, 2, 3, 4)]
        halving = 1 - (1 - out[(n_cs, half[0])]) * (1 - out[(n_cs, half[1])]) ** 2
        k_star = math.ceil(math.log(1 - d) / math.log(1 - p))
        print(f"n_cs = {n_cs}, delta = {d}: same delta, K = 1..4: " + ", ".join(f"{x:.3g}" for x in same)
              + f"; halving K = 3: {halving:.3g}; naive guard exceeds delta first at K* = {k_star} windows "
              f"(1 - (1 - p)^K* = {1 - (1 - p) ** k_star:.4f}); uncapped naive guard: P(kept window fails) = 1")


if __name__ == "__main__":
    main()
