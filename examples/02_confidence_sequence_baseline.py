"""The confidence-sequence baseline: a false alarm guarantee with an estimated baseline.

[SYNTHETIC] Student t returns with 5 degrees of freedom (unit variance); the standard deviation rises
by half at step 3000.

The guarantee P(false alarm within the first m monitored steps) <= m / A assumes the true null mean of
the truncated squared return. The default monitor plugs in the sample mean of the calibration window, and
an estimate that comes out too low makes the increase side fire too early. With ``baseline="cs"`` the
calibration window is split: ``tau`` comes from its first quarter, and a time-uniform confidence sequence
on the rest gives an upper bound U and a lower bound L for the null mean, each at level delta / 2. The
increase side bets against U and the decrease side against L, and then, with no change,
P(any false alarm within the first m monitored steps) <= delta + m / A, with the baseline estimated.
The price is a longer detection delay, because U and L are further from the mean than the plug-in value.

``cs_baseline`` computes the same bounds for calibration windows held in an array.

Run: python examples/02_confidence_sequence_baseline.py
"""

import numpy as np

from avmonitor import VolatilityMonitor, cs_baseline

CHANGE = 3000
CAL_LEN = 1000
DELTA = 0.01
A = 390 / 0.04  # with DELTA: P(false alarm within the first 390 monitored steps) <= 0.01 + 0.04 = 5%


def main() -> None:
    rng = np.random.default_rng(2)
    t5 = rng.standard_t(5, CHANGE + 1500) / np.sqrt(5 / 3)  # unit variance
    r = np.concatenate([t5[:CHANGE], 1.5 * t5[CHANGE:]])
    print(f"[SYNTHETIC] {r.size} Student t5 returns, sd 1 before step {CHANGE} and 1.5 from it on")
    print(f"A = {A:.0f}, delta = {DELTA}: proved bound delta + 390 / A = {DELTA + 390 / A:.0%} "
          "within the first 390 monitored steps")

    # Bounds on the null mean from the calibration window (array version).
    b = cs_baseline(r[:CAL_LEN], c=3.0, delta_upper=DELTA / 2, delta_lower=DELTA / 2, tau_frac=0.25)
    print(f"calibration window: tau {b.tau[0]:.3f} from the first {b.n_tau} returns; "
          f"mean of min(r^2, tau^2) over the other {int(b.n_cs[0])}: {b.mean[0]:.3f}, "
          f"confidence bounds L = {b.lower[0]:.3f}, U = {b.upper[0]:.3f}")

    for baseline in ("plugin", "cs"):
        mon = VolatilityMonitor(threshold=A, cal_len=CAL_LEN, baseline=baseline, delta=DELTA)
        events = mon.run(r)
        seg = mon.segments[0]
        print(f"baseline={baseline!r}: increase side bets against {seg['m_up']:.3f}, "
              f"decrease side against {seg['m_dn']:.3f}")
        if events:
            e = events[0]
            when = (f"{e.index - CHANGE + 1} steps after the change" if e.index >= CHANGE
                    else "before the change (false alarm)")
            print(f"    first alarm at step {e.index} ({e.direction}), {when}")
        else:
            print("    no alarm")


if __name__ == "__main__":
    main()
