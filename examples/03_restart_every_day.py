"""A monitor that restarts its statistics at every day boundary, for a bound on every day.

[SYNTHETIC] Ten days of 390 one-minute Gaussian returns with sd 1. Volatility doubles at minute 200 of
day 7.

Theorem B bounds the chance of a false alarm within the FIRST m steps after the statistics start by
m / A. It does not bound a later day of the same run: the Shiryaev-Roberts statistic grows like t under
the null, so on day d Ville's inequality gives only d m / A. ``PeriodicVolatilityMonitor(period=390)``
resets both statistics to 0 at every day boundary. What that proves depends on the baseline
(``docs/theory_per_period.md``):

* known baseline (``baseline="known"`` with ``tau`` and ``m0``): P(false alarm on a day | everything
  before it) <= 390 / A for every day (``restarted_sr_bound``). The data here are simulated, so the true
  values are known: tau = 3 and m0 = E[min(Z^2, 9)] for a standard normal Z;
* confidence-sequence baseline (``baseline="cs"``): P(false alarm on a given day of a segment)
  <= delta + 390 / A, not given the past, while the tolerance of the confidence sequence holds;
* plug-in baseline (the default): no per-day bound is proved.

The cost of the restart is that evidence gathered before a boundary is dropped at the boundary. One
synthetic run shows how to use the class; it says nothing about false alarm rates.

Run: python examples/03_restart_every_day.py
"""

import numpy as np

from avmonitor import PeriodicVolatilityMonitor, VolatilityMonitor, restarted_sr_bound
from avmonitor.robustness import trunc_mean_normal

DAY = 390
DAYS = 10
CHANGE = 7 * DAY + 200
DELTA = 0.01
A = DAY / 0.04  # known baseline: 4% per day; confidence-sequence baseline: delta + 4% = 5% per day


def main() -> None:
    rng = np.random.default_rng(3)
    sd = np.where(np.arange(DAYS * DAY) >= CHANGE, 2.0, 1.0)
    r = sd * rng.standard_normal(DAYS * DAY)
    tau, m0 = 3.0, float(trunc_mean_normal(3.0))
    print(f"[SYNTHETIC] {DAYS} days of {DAY} Gaussian returns; sd doubles at day {CHANGE // DAY}, "
          f"minute {CHANGE % DAY}")
    print(f"A = {A:.0f}. Daily restart, known baseline (tau = {tau}, m0 = {m0:.4f}): "
          f"P(false alarm on a day | past) <= {restarted_sr_bound(DAY, A):.0%} for every day.")
    print(f"Daily restart, confidence-sequence baseline (delta = {DELTA}): "
          f"P(false alarm on a given day) <= {DELTA + restarted_sr_bound(DAY, A):.0%}.")
    print(f"No restart, confidence-sequence baseline: only delta + d * {DAY} / A for monitored day d, "
          f"for example {min(1.0, DELTA + 9 * DAY / A):.0%} for day 9.")

    monitors = {
        "daily restart, known baseline": PeriodicVolatilityMonitor(
            threshold=A, period=DAY, offset=0, baseline="known", tau=tau, m0=m0),
        "daily restart, confidence-sequence baseline": PeriodicVolatilityMonitor(
            threshold=A, cal_len=DAY, period=DAY, offset=0, baseline="cs", delta=DELTA),
        "no restart, confidence-sequence baseline": VolatilityMonitor(
            threshold=A, cal_len=DAY, baseline="cs", delta=DELTA),
    }
    for name, mon in monitors.items():
        events = mon.run(r)
        false = [e for e in events if e.index < CHANGE]
        true = [e for e in events if e.index >= CHANGE]
        line = f"{name}: {len(false)} alarm(s) before the change"
        if true:
            e = true[0]
            line += (f"; first alarm after it on day {e.index // DAY}, minute {e.index % DAY} ({e.direction}), "
                     f"{e.index - CHANGE + 1} minutes after the change; {len(true)} alarm(s) after it")
        print(line)
    print("(A known baseline is never re-estimated, so after a lasting change it keeps alarming: each alarm")
    print(" restarts the statistics against the old baseline. The other two recalibrate after an alarm.)")


if __name__ == "__main__":
    main()
