"""Quickstart: calibrate, monitor, alarm.

[SYNTHETIC] Independent Gaussian returns whose standard deviation doubles at step 3000.

``VolatilityMonitor`` calibrates on the first ``cal_len`` returns (robust scale, truncation level
``tau = 3 * scale`` and the null mean ``m0`` of the truncated squared return), then updates two
Shiryaev-Roberts e-detectors, one for an increase and one for a decrease of the truncated variance.
Each side alarms at 2 A. With the true baseline, the chance of any false alarm within the first m
monitored steps is at most m / A; here A = 390 / 0.05, so at most 5% within the first 390 steps. (For a
bound on every later day as well, see 03_restart_every_day.py.) After an alarm the monitor recalibrates
on fresh data.

The second part runs the same detectors on arrays: ``calibrate`` for the baseline, then
``betting_edetector`` for each side and ``two_sided_alarm`` for the first alarm of either.

Run: python examples/01_quickstart.py
"""

import numpy as np

from avmonitor import VolatilityMonitor, betting_edetector, calibrate, two_sided_alarm

CHANGE = 3000  # first step with the new standard deviation
CAL_LEN = 1000
A = 390 / 0.05  # alarm level A = 7800


def main() -> None:
    rng = np.random.default_rng(1)
    r = np.concatenate([rng.standard_normal(CHANGE), 2.0 * rng.standard_normal(1000)])
    print(f"[SYNTHETIC] {r.size} Gaussian returns, sd 1 before step {CHANGE} and 2 from it on")
    print(f"alarm level A = {A:.0f}: with the true baseline, P(false alarm within the first 390 monitored steps) "
          f"<= 390 / A = {390 / A:.0%}; this run estimates the baseline (plug-in), so its level is approximate")

    # 1. Streaming: feed returns one at a time (run() loops over update()).
    mon = VolatilityMonitor(threshold=A, cal_len=CAL_LEN)
    events = mon.run(r)
    seg = mon.segments[0]
    print(f"calibration on steps 0 to {CAL_LEN - 1}: scale {seg['scale']:.3f}, tau {np.sqrt(seg['tau2']):.3f}, "
          f"m0 {seg['m0']:.3f}")
    for e in events:
        when = (f"{e.index - CHANGE + 1} steps after the change" if e.index >= CHANGE
                else "before the change (false alarm)")
        print(f"streaming alarm: step {e.index}, direction {e.direction}, {when}")
    if not events:
        print("streaming monitor: no alarm")

    # 2. Arrays: the same calibration and detectors, one call each.
    base = calibrate(r[:CAL_LEN], c=3.0)
    x = np.minimum(r[CAL_LEN:] ** 2, base.tau**2)  # truncated squared returns of the monitored part
    up = betting_edetector(x, base.m0_up, "up", threshold=2 * A)
    down = betting_edetector(x, base.m0_down, "down", upper=base.tau**2, threshold=2 * A)
    first = int(two_sided_alarm(up, down)[0])
    if first >= 0:
        side = "up" if up.alarm[0] == first else "down"
        print(f"array detectors: first alarm at step {CAL_LEN + first} ({side})")
    else:
        print("array detectors: no alarm")


if __name__ == "__main__":
    main()
