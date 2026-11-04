"""Round 9 check: the batch functions and the streaming classes give the same alarms on the same data.

Four comparisons, each on many seeded paths with a volatility change so that alarms occur:
  A. known baseline with restarts: PeriodicVolatilityMonitor(baseline="known") against repeated calls of
     betting_edetector (up and down, each at 2A) + two_sided_alarm, restarted after each alarm;
  B. plug-in baseline: VolatilityMonitor against calibrate() + betting_edetector (first alarm);
  C. confidence-sequence baseline: VolatilityMonitor(baseline="cs") against calibrate(baseline="cs") +
     betting_edetector (first alarm);
  D. as A, with 10% of the returns missing (NaN), for the package's betting_edetector and for the
     proposed fix (patched_betting_edetector below, the M1 diff).
Returns: Gaussian and Student t3 (unit scale).
"""
import math

import numpy as np
from scipy import integrate, stats

from avmonitor import VolatilityMonitor, betting_edetector, calibrate, two_sided_alarm
from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID, DetectorResult
from avmonitor.periodic import PeriodicVolatilityMonitor
from avmonitor.robustness import trunc_mean_normal


def patched_betting_edetector(x, m0, direction="up", upper=None, threshold=1000.0, kind="SR", grid=None):
    """betting_edetector with the M1 fix (fixed grid only): a NaN step leaves every statistic unchanged."""
    x = np.atleast_2d(np.asarray(x, float))
    reps, n = x.shape
    y = x / m0
    if direction == "up":
        grid = np.asarray(DEFAULT_UP_GRID if grid is None else grid, float); sign = 1.0
    else:
        b = upper / m0
        grid = (1.0 / (b - 1.0)) * np.asarray(DEFAULT_DOWN_FRACTIONS if grid is None else grid, float); sign = -1.0
    w = np.full(len(grid), 1.0 / len(grid))
    S = np.zeros((reps, len(grid)))
    alarm = np.full(reps, -1)
    active = np.ones(reps, bool)
    for t in range(n):
        g = sign * (y[:, t] - 1.0)
        ok = np.isfinite(g)
        f = 1.0 + grid[None, :] * np.where(ok, g, 0.0)[:, None]
        new = f * (S + 1.0) if kind == "SR" else f * np.maximum(S, 1.0)
        S = np.where(ok[:, None], new, S)
        hit = active & (S @ w >= threshold)
        alarm[hit] = t
        active &= ~hit
    return DetectorResult(alarm)


def t3_trunc_mean(tau):
    f = lambda v: 2 * stats.t.sf(math.sqrt(v), 3)  # E[min(T^2, tau^2)] = int_0^{tau^2} P(T^2 > v) dv
    return integrate.quad(f, 0, tau * tau, limit=200)[0]


def draw(law, rng, n):
    return rng.standard_normal(n) if law == "gauss" else rng.standard_t(3, n)


def batch_events_known(r, tau, m0, A, det):
    x = np.minimum(r**2, tau**2)
    out, s = [], 0
    while s < r.size:
        up = det(x[s:], m0, "up", threshold=2 * A)
        dn = det(x[s:], m0, "down", upper=tau**2, threshold=2 * A)
        a = int(two_sided_alarm(up, dn)[0])
        if a < 0:
            break
        out.append((s + a, "up" if up.alarm[0] == a else "down"))
        s += a + 1
    return out


def main():
    A = 60.0
    P = 150
    tau_g, tau_t = 3.0, 3.0 * stats.t.ppf(0.75, 3) / stats.norm.ppf(0.75)  # tau = 3 x population scaled MAD
    laws = {"gauss": (tau_g, float(trunc_mean_normal(tau_g))), "t3": (tau_t, t3_trunc_mean(tau_t))}
    for law, (tau, m0) in laws.items():
        nA = nB = nC = nD_pkg = nD_fix = 0
        evA = evB = evC = evD = 0
        for i in range(P):
            rng = np.random.default_rng(1000 + i)
            r = np.concatenate([draw(law, rng, 1500), 1.6 * draw(law, rng, 700), 0.6 * draw(law, rng, 800)])
            # A
            mon = PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=tau, m0=m0)
            s = [(e.index, e.direction) for e in mon.run(r)]
            b = batch_events_known(r, tau, m0, A, betting_edetector)
            nA += s != b; evA += len(s)
            # B (first alarm after a 400-step calibration)
            mon = VolatilityMonitor(threshold=A, cal_len=400)
            mon.run(r)
            s = (mon.events[0].index, mon.events[0].direction) if mon.events else None
            base = calibrate(r[:400])
            x = np.minimum(r[400:] ** 2, base.tau**2)
            up = betting_edetector(x, base.m0_up, "up", threshold=2 * A)
            dn = betting_edetector(x, base.m0_down, "down", upper=base.tau**2, threshold=2 * A)
            a = int(two_sided_alarm(up, dn)[0])
            b = (400 + a, "up" if up.alarm[0] == a else "down") if a >= 0 else None
            nB += s != b; evB += s is not None
            # C
            mon = VolatilityMonitor(threshold=A, cal_len=400, baseline="cs", delta=0.01)
            mon.run(r)
            s = (mon.events[0].index, mon.events[0].direction) if mon.events else None
            base = calibrate(r[:400], baseline="cs", delta=0.01)
            x = np.minimum(r[400:] ** 2, base.tau**2)
            up = betting_edetector(x, base.m0_up, "up", threshold=2 * A)
            dn = betting_edetector(x, base.m0_down, "down", upper=base.tau**2, threshold=2 * A)
            a = int(two_sided_alarm(up, dn)[0])
            b = (400 + a, "up" if up.alarm[0] == a else "down") if a >= 0 else None
            nC += s != b; evC += s is not None
            # D: 10% missing
            rn = r.copy()
            rn[rng.random(r.size) < 0.10] = np.nan
            mon = PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=tau, m0=m0)
            s = [(e.index, e.direction) for e in mon.run(rn)]
            nD_pkg += s != batch_events_known(rn, tau, m0, A, betting_edetector)
            nD_fix += s != batch_events_known(rn, tau, m0, A, patched_betting_edetector)
            evD += len(s)
        print(f"[SYNTHETIC] {law} (tau {tau:.4f}, m0 {m0:.4f}), {P} paths, A = {A:g}:", flush=True)
        print(f"  A known baseline with restarts: {nA} paths differ ({evA} streaming alarms in all)")
        print(f"  B plug-in, first alarm:         {nB} paths differ ({evB} paths with an alarm)")
        print(f"  C cs baseline, first alarm:     {nC} paths differ ({evC} paths with an alarm)")
        print(f"  D 10% NaN, known, restarts:     {nD_pkg} paths differ with betting_edetector, "
              f"{nD_fix} with the M1 fix ({evD} streaming alarms)")


if __name__ == "__main__":
    main()
