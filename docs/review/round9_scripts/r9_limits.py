"""Round 9 spot checks of limits.py public functions. [SYNTHETIC]

1. sr_alarm_times against betting_edetector (same bets, equal weights, same threshold) on Gaussian-based Y paths.
2. collapse_floor against floor_time on the constant extreme path, and against a direct search.
3. max_null_mass: the two-point law it names has mean exactly 1 and puts that mass on y.
Usage: python r9_limits.py
"""
import math

import numpy as np

from avmonitor import betting_edetector
from avmonitor.limits import collapse_floor, floor_time, max_null_mass, sr_alarm_times
from avmonitor.robustness import trunc_mean_normal

TAU = 3.0
M0 = float(trunc_mean_normal(TAU))
B = TAU * TAU / M0


def check_sr_alarm_times():
    rng = np.random.default_rng(21)
    worst = 0
    n_alarm = 0
    for direction, scale in (("up", 1.6), ("down", 0.55)):
        x = np.minimum((scale * rng.standard_normal((200, 3000))) ** 2, TAU * TAU)
        y = x / M0
        # betting_edetector takes the decrease grid as fractions of 1 / (B - 1); sr_alarm_times takes the bets
        grid = (0.02, 0.05, 0.1, 0.2, 0.35, 0.5) if direction == "up" else (0.1, 0.25, 0.5, 0.75, 0.95)
        bets = np.asarray(grid) if direction == "up" else np.asarray(grid) / (B - 1.0)
        for A in (50.0, 400.0):
            ref = betting_edetector(x, M0, direction, upper=TAU * TAU, threshold=A, grid=grid).alarm
            got = sr_alarm_times(y, bets, np.full(len(grid), 1.0 / len(grid)), direction, thresholds=(A,))[0]
            ref1 = np.where(ref >= 0, ref + 1, -1)
            worst = max(worst, int(np.sum(ref1 != got)))
            n_alarm += int(np.sum(got > 0))
    print(f"1. sr_alarm_times against betting_edetector: {worst} differing paths at worst ({n_alarm} alarms compared)")


def check_collapse_floor():
    bad = 0
    cases = 0
    for A in (50.0, 300.0, 7800.0, 28800.0):
        for Bv in (2.0, 5.0, 9.045, 20.0):
            for side in ("down", "up"):
                F = Bv / (Bv - 1.0) if side == "down" else Bv
                n_cf = collapse_floor(A, Bv, side)
                n_ft = int(floor_time(np.full(10**5, math.log(F)), A)[0])
                n_bf = next(n for n in range(1, 10**6) if n * math.log(F) >= math.log(A / n) - 1e-12)
                cases += 1
                bad += not (n_cf == n_ft == n_bf)
    print(f"2. collapse_floor = floor_time on the constant path = direct search: {cases - bad} of {cases} cases agree")


def check_max_null_mass():
    worst = 0.0
    for Bv in (2.0, 9.045, 50.0):
        for y in np.linspace(0, Bv, 41):
            p_dn = float(max_null_mass(y, Bv, "down"))
            if y < 1:  # two-point law on {y, B}
                worst = max(worst, abs(p_dn * y + (1 - p_dn) * Bv - 1.0))
            p_up = float(max_null_mass(y, Bv, "up"))
            if y > 1:  # two-point law on {0, y}
                worst = max(worst, abs(p_up * y - 1.0))
    print(f"3. max_null_mass: largest |mean of the named two-point law - 1| = {worst:.1e}")


if __name__ == "__main__":
    check_sr_alarm_times()
    check_collapse_floor()
    check_max_null_mass()
