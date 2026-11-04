"""Round 9 check: missing values in detectors.betting_edetector versus VolatilityMonitor.

The docstring of betting_edetector says a NaN is skipped ("the bet factor is 1 at that step"), and
benchmarks_robust.sr_log_paths says NaN "leave the statistic unchanged, as in detectors.betting_edetector
and monitor.VolatilityMonitor". This script checks both claims.
"""
import numpy as np

from avmonitor import betting_edetector, two_sided_alarm
from avmonitor.periodic import PeriodicVolatilityMonitor
from avmonitor.benchmarks_robust import sr_log_paths

A = 7800.0
tau, m0 = 3.0, 1.0  # any valid known baseline

# 1. A stream with no observation at all.
x = np.full(2 * int(A) + 10, np.nan)
up = betting_edetector(x, m0, "up", threshold=2 * A)
dn = betting_edetector(x, m0, "down", upper=tau**2, threshold=2 * A)
print("all-NaN stream of", x.size, "steps: batch alarm (up, down, two-sided) =",
      int(up.alarm[0]), int(dn.alarm[0]), int(two_sided_alarm(up, dn)[0]))
mon = PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=tau, m0=m0)
mon.run(np.sqrt(x * m0))
print("all-NaN stream: streaming monitor events =", len(mon.events))

# 2. Leading NaN (asset that starts trading late on a common clock), then ONE calm observation.
x = np.r_[np.full(2 * int(A) - 1, np.nan), 1.0]  # Y = 1 exactly: no evidence of anything
up = betting_edetector(x, m0, "up", threshold=2 * A, return_path=True)
print("2A-1 NaN then one observation Y = 1: batch up alarm at", int(up.alarm[0]),
      "statistic", up.stat_path[-1])
mon = PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=tau, m0=m0)
mon.run(np.sqrt(np.nan_to_num(x, nan=np.nan) * m0))
print("same stream, streaming monitor: events =", len(mon.events), "S_up mean =", mon.S_up.mean())

# 3. sr_log_paths on the same input.
lp = sr_log_paths(np.r_[np.full(5, np.nan), 1.0][None], np.asarray([0.1, 0.2]), "up")
bp = betting_edetector(np.r_[np.full(5, np.nan), 1.0], 1.0, "up", grid=[0.1, 0.2], threshold=1e9, return_path=True)
print("5 NaN then Y = 1: sr_log_paths M_6 =", float(np.exp(lp[0, -1])), " betting_edetector M_6 =", bp.stat_path[-1])

# 4. A session market on a 24-hour minute grid (375 observed minutes, then 1065 NaN, per day), i.i.d. Gaussian
#    returns with the true baseline, A = 375 / 0.05 = 7500 (the CLI's level for Nifty). [SYNTHETIC]
from avmonitor.robustness import trunc_mean_normal

M0 = float(trunc_mean_normal(3.0))
A, D, P = 7500.0, 12, 400
rng = np.random.default_rng(4)
day = np.r_[np.ones(375, bool), np.zeros(1065, bool)]
obs = np.tile(day, D)
r = np.full((P, obs.size), np.nan)
r[:, obs] = rng.standard_normal((P, obs.sum()))
x = np.minimum(r**2, 9.0)
a = two_sided_alarm(betting_edetector(x, M0, "up", threshold=2 * A),
                    betting_edetector(x, M0, "down", upper=9.0, threshold=2 * A))
n_obs_before = np.cumsum(obs)  # observed returns up to and including each grid step
batch_obs = np.where(a >= 0, n_obs_before[np.maximum(a, 0)], -1)
stream_hit = np.zeros(P, bool)
for i in range(P):
    mon = PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=3.0, m0=M0)
    stream_hit[i] = len(mon.run(r[i])) > 0
m_obs = D * 375
print(f"[SYNTHETIC] session market on a 24-hour grid, {D} days ({m_obs} observed returns), A = {A:g}, {P} paths:")
print(f"  batch betting_edetector: P(alarm) = {np.mean(a >= 0):.3f}; median observed returns before the alarm "
      f"{np.median(batch_obs[a >= 0]) if (a >= 0).any() else float('nan'):.0f}")
print(f"  streaming monitor:       P(alarm) = {stream_hit.mean():.3f}")
print(f"  Theorem B bound with m = observed returns: {min(1.0, m_obs / A):.3f}")

# 5. A degenerate calibration window (more than half zero returns) through the documented array API.
import warnings

from avmonitor import calibrate

rng = np.random.default_rng(0)
cal = np.r_[np.zeros(300), 0.001 * rng.standard_normal(200)]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    b = calibrate(cal)
    print("degenerate window: calibrate returns", b)
    x = np.minimum((0.001 * rng.standard_normal(2000)) ** 2, b.tau**2)
    up = betting_edetector(x, b.m0_up, "up", threshold=2 * 50)
    print("  betting_edetector(up, A = 50) alarm at", int(up.alarm[0]), "(= 2A - 1, whatever the data)")
    try:
        betting_edetector(x, b.m0_down, "down", upper=b.tau**2, threshold=2 * 50)
    except ValueError as e:
        print("  down side raises:", e)
