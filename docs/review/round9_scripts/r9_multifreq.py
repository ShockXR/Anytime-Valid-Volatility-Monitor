"""Round 9 checks of MultiFrequencyMonitor (certified average, the default).

1. Pathwise: between block ends, the certified value of a coarse frequency is never above that frequency's
   full statistic at the next block end (the property Theorem 1 of docs/theory_multifreq.md rests on),
   with 5% missing returns.
2. freqs=(1,) gives exactly the events of VolatilityMonitor (plug-in), as the docstring says.
3. [SYNTHETIC] Monte Carlo under the null with the known baseline (i.i.d. Gaussian, freqs 1, 5, 15,
   balanced weights, A = level_for_budget(0.05, 390)): P(false alarm within 390 minutes) <= 5%.
"""
import math
import sys

import numpy as np

from avmonitor import VolatilityMonitor
from avmonitor.multifreq import (MultiFrequencyMonitor, block_returns, component_paths, first_alarm,
                                 level_for_budget, minute_values, resolve_weights)
from avmonitor.robustness import trunc_mean_normal

M0 = float(trunc_mean_normal(3.0))
FREQS = (1, 5, 15)
KNOWN = {k: (9.0 * k, k * M0) for k in FREQS}


def pathwise(P=40):
    """Block positions are counted from each segment's start (the monitor restarts its block grid after an alarm)."""
    worst = -np.inf
    checked = 0
    n_alarms = 0
    for i in range(P):
        rng = np.random.default_rng(i)
        r = np.concatenate([rng.standard_normal(1500), 1.4 * rng.standard_normal(3000)])
        r[rng.random(r.size) < 0.05] = np.nan
        mon = MultiFrequencyMonitor(threshold=1e6, freqs=FREQS, cal_len=1500, known_baseline=KNOWN, trace=True)
        mon.run(r)
        n_alarms += len(mon.events)
        t = np.array([row[0] for row in mon.trace_rows])
        vu = np.array([row[3] for row in mon.trace_rows])
        vd = np.array([row[4] for row in mon.trace_rows])
        starts = [seg["start"] for seg in mon.segments]
        ends = [e.index for e in mon.events] + [r.size - 1]
        for s0, e0 in zip(starts, ends):
            sel = (t >= s0) & (t <= e0)
            pos = t[sel] - s0
            assert np.array_equal(pos, np.arange(pos.size))
            for j, k in enumerate(FREQS):
                if k == 1:
                    continue
                for v in (vu[sel, j], vd[sel, j]):
                    for b0 in range(0, v.size - k + 1, k):
                        end = b0 + k - 1
                        if end == v.size - 1 and e0 != r.size - 1:
                            break  # the segment ended with an alarm at this minute: the block end is the alarm
                        inside = v[b0:end]
                        if inside.size:
                            worst = max(worst, float(np.max(inside - v[end])))
                            checked += inside.size
    print(f"1. certified value minus full value at the next block end: largest {worst:.3g} over {checked} minutes "
          f"in {P} paths with {n_alarms} alarms (must be <= 0)")


def single_freq():
    diffs = 0
    tot = 0
    for i in range(30):
        rng = np.random.default_rng(100 + i)
        r = np.concatenate([rng.standard_normal(2000), 1.7 * rng.standard_normal(800), 0.6 * rng.standard_normal(800)])
        r[rng.random(r.size) < 0.03] = np.nan
        a = VolatilityMonitor(threshold=100, cal_len=500).run(r)
        b = MultiFrequencyMonitor(threshold=100, freqs=(1,), cal_len=500).run(r)
        diffs += [(e.index, e.direction) for e in a] != [(e.index, e.direction) for e in b]
        tot += len(a)
    print(f"2. freqs=(1,) against VolatilityMonitor: {diffs} of 30 paths differ ({tot} alarms in all)")


def vectorised_first_alarm(r, A, freqs=FREQS, weights="balanced", known=KNOWN):
    """First alarm of the certified average from the vectorised building blocks (as in tests/test_multifreq.py)."""
    r = np.atleast_2d(r)
    n = r.shape[1]
    w = resolve_weights(weights, freqs)
    up = np.zeros_like(r)
    dn = np.zeros_like(r)
    for wk, k in zip(w, freqs):
        tau2, m0 = known[k]
        x = np.minimum(block_returns(r, k) ** 2, tau2)
        p = component_paths(x, m0, m0, tau2)
        up += wk * minute_values(p["up_full"], p["up_cert"], k, n, True)
        dn += wk * minute_values(p["dn_full"], p["dn_cert"], k, n, True)
    return first_alarm(up, dn, A)


def vec_vs_stream_nan():
    rng = np.random.default_rng(21)
    cal, n, A = 60, 1500, 200.0
    diffs = 0
    hits = 0
    for rep in range(30):
        ratio = [1.0, 1.5, 0.6][rep % 3]
        r = rng.standard_normal(cal + n) * np.where(np.arange(cal + n) >= cal + 400, ratio, 1.0)
        r[cal:][rng.random(n) < 0.05] = np.nan
        mon = MultiFrequencyMonitor(threshold=A, freqs=FREQS, cal_len=cal, known_baseline=KNOWN)
        mon.run(r)
        idx, d = vectorised_first_alarm(r[cal:], A)
        s_ = (mon.events[0].index - cal, 1 if mon.events[0].direction == "up" else -1) if mon.events else (-1, 0)
        diffs += s_ != (int(idx[0]), int(d[0]))
        hits += idx[0] >= 0
    print(f"2b. vectorised building blocks against MultiFrequencyMonitor with 5% NaN: {diffs} of 30 paths differ "
          f"({hits} with an alarm)")


def mc(P):
    m = 390
    A = level_for_budget(0.05, m, FREQS)
    rng = np.random.default_rng(5)
    k = 0
    for i in range(P):
        mon = MultiFrequencyMonitor(threshold=A, freqs=FREQS, cal_len=30, known_baseline=KNOWN)
        r = rng.standard_normal(30 + m)
        for v in r:
            if mon.update(v) is not None:
                k += 1
                break
    p = k / P
    se = math.sqrt(max(p * (1 - p), 1e-300) / P)
    print(f"3. [SYNTHETIC] multi-frequency (streaming), known baseline, Gaussian, A = {A:.0f}, m = {m}: "
          f"rate {p:.4%} (se {se:.4%}) over {P} paths, proved bound {mon.proved_bound(m):.2%}")
    Pv = 10 * P
    idx, _ = vectorised_first_alarm(np.random.default_rng(6).standard_normal((Pv, m)), A)
    p = float(np.mean(idx >= 0))
    se = math.sqrt(max(p * (1 - p), 1e-300) / Pv)
    print(f"   [SYNTHETIC] the same with the vectorised building blocks: rate {p:.4%} (se {se:.4%}) over {Pv} paths")


if __name__ == "__main__":
    pathwise()
    single_freq()
    vec_vs_stream_nan()
    mc(int(sys.argv[1]) if len(sys.argv) > 1 else 4000)
