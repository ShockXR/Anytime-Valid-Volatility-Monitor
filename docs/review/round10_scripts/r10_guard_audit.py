"""Review round 10 (E24): does src/avmonitor/guard.py implement what Theorems 1 to 3 assume?

Checks (all [SYNTHETIC], seconds to a few minutes):
  1. budgets: the k-th window's sequences are built at level delta_k / 2 per side (also with c_down != c);
  2. the discard decision is taken at the end of each window, windows are consecutive and disjoint;
  3. the rank fallback restarts its statistics at every period boundary (and not at its first step);
  4. guard_tries=1 is the unguarded monitor only with the "keep" fallback (docstring says "the unguarded one");
  5. after a run, mon.delta holds the current window's budget, not the total delta;
  6. Monte Carlo, Theorem 3 part 2 on the rank branch: per-period alarm rate after a fallback, exchangeable data.

Run: nice -n 10 python docs/review/round10_scripts/r10_guard_audit.py
"""

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from avmonitor.guard import GuardedPeriodicVolatilityMonitor, GuardedVolatilityMonitor  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402
from avmonitor.robust import MAD_NORMAL_SCALE  # noqa: E402

CAL, NT = 200, 50
SMALL = dict(c=3.0, cal_len=CAL, null_band=0.25, baseline="cs")
TAU = 3.0 * MAD_NORMAL_SCALE


def blind_window(rng, big=10.0):
    return np.r_[rng.choice([-1.0, 1.0], NT), rng.choice([-big, big], CAL - NT)]


def record_budgets(mon):
    """Wrap _finish to record, at the end of each window, (try index, delta attribute, delta_u and delta_l used)."""
    rec = []
    orig = mon._finish

    def wrapped(*a):
        cs_dn = mon._cs if mon._cs_dn is None else mon._cs_dn
        rec.append(dict(t=mon.t, try_index=mon._try, delta_attr=mon.delta, delta_u=math.exp(-mon._cs.log_u),
                        delta_l=math.exp(-cs_dn.log_l), buffer=len(mon.buffer)))
        return orig(*a)

    mon._finish = wrapped
    return rec


def check_budgets():
    rng = np.random.default_rng(1)
    u = np.r_[blind_window(rng), blind_window(rng, 20.0), rng.standard_normal(CAL) * 10, rng.standard_normal(300) * 10]
    for kw in ({}, dict(c_down=4.0)):
        g = GuardedVolatilityMonitor(threshold=1e6, delta=0.04, guard_tries=3, **SMALL, **kw)
        rec = record_budgets(g)
        g.run(u)
        print("1. budgets", kw or "(c_down = c)", "schedule", g.guard_deltas)
        for r in rec:
            print(f"   window {r['try_index']}: ends at t = {r['t']}, length {r['buffer']}, delta attribute "
                  f"{r['delta_attr']:.4f}, delta_u = {r['delta_u']:.5f}, delta_l = {r['delta_l']:.5f}")
        ok = all(abs(r["delta_u"] - g.guard_deltas[r["try_index"] - 1] / 2) < 1e-12 for r in rec)
        ok &= all(abs(r["delta_l"] - g.guard_deltas[r["try_index"] - 1] / 2) < 1e-12 for r in rec)
        print("   every window at delta_k / 2 per side:", ok)
        starts = [d["window_start"] for d in g.discarded] + [g.segments[0]["start"] - g.segments[0]["cal_len"]]
        ends = [d["window_end"] for d in g.discarded] + [g.segments[0]["start"] - 1]
        print("2. windows (start, end):", list(zip(starts, ends)), "consecutive and disjoint:",
              all(starts[i + 1] == ends[i] + 1 for i in range(len(ends) - 1)),
              "| decisions at", [r["t"] for r in rec], "= window ends:", [r["t"] for r in rec] == ends)


def check_rank_restart():
    """At every step of the fallback segment, record the rank statistics as rm.update receives them: they must be 0
    exactly at the period boundaries after the fallback's first step, and equal to the previous step's output
    elsewhere."""
    rng = np.random.default_rng(2)
    m = 97  # a period that does not divide the window length, so boundaries fall at arbitrary segment steps
    g = GuardedPeriodicVolatilityMonitor(threshold=1e9, delta=0.05, guard_tries=2, guard_fallback="rank", rank_seed=3,
                                         period=m, offset=5, **SMALL)
    u = np.r_[blind_window(rng), blind_window(rng), rng.standard_normal(1000)]
    entries, exits = [], []
    wrapped = False
    for x in u:
        if g._rank is not None and not wrapped:
            rm = g._rank
            orig = rm.update

            def upd(r, _orig=orig, _rm=rm):
                entries.append((g.t, _rm.S_up.copy(), _rm.S_dn.copy()))
                out = _orig(r)
                exits.append((_rm.S_up.copy(), _rm.S_dn.copy()))
                return out

            rm.update = upd
            wrapped = True
            start = g.start
            resets0 = g.period_resets
        g.update(x)
    ok_b = ok_c = True
    n_b = 0
    for i, (t, su, sd) in enumerate(entries):
        boundary = (t - g.offset) % m == 0 and t > start
        if boundary:
            n_b += 1
            ok_b &= bool(np.all(su == 0) and np.all(sd == 0))
        elif i > 0:
            ok_c &= bool(np.array_equal(su, exits[i - 1][0]) and np.array_equal(sd, exits[i - 1][1]))
    print(f"3. rank fallback from t = {start} (first step {entries[0][0]}): {n_b} boundaries, statistics 0 at each: "
          f"{ok_b}; carried unchanged elsewhere: {ok_c}; period_resets added {g.period_resets - resets0}; "
          f"first step at a boundary: {(start - g.offset) % m == 0}")


def check_identity_rank_k1():
    rng = np.random.default_rng(4)
    u = np.r_[blind_window(rng), rng.standard_normal(400) * 10]
    a = VolatilityMonitor(threshold=50.0, delta=0.05, **SMALL)
    b = GuardedVolatilityMonitor(threshold=50.0, delta=0.05, guard_tries=1, guard_fallback="rank", rank_seed=1, **SMALL)
    a.run(u)
    b.run(u)
    print("4. guard_tries=1, rank fallback, blind first window: unguarded segment", a.segments[0]["tau2"] /
          a.segments[0]["m_up"] <= 1, "(blind kept), up alarms", sum(e.direction == "up" for e in a.events),
          "| guarded:", b.segments[0]["guard"], ", up alarms", sum(e.direction == "up" for e in b.events),
          "-> not the unguarded monitor")


def check_delta_attribute():
    rng = np.random.default_rng(5)
    u = np.r_[blind_window(rng), rng.standard_normal(CAL), rng.standard_normal(100)]
    g = GuardedVolatilityMonitor(threshold=1e6, delta=0.01, guard_tries=3, **SMALL)
    g.run(u)
    print(f"5. after a run: mon.delta = {g.delta} (window {g._try}'s budget), mon.delta_total = {g.delta_total}")


def mc_rank_period(paths=400, seed=6):
    """Theorem 3 part 2 on the rank branch: K = 1, rank fallback, period m; i.i.d. data after the fallback.
    Per-period alarm rate of the fallback segment for its first periods, against m / A."""
    rng = np.random.default_rng(seed)
    m, A = 100, 500.0
    n_per = 4
    hits = np.zeros(n_per)
    alive = np.zeros(n_per)
    falls = 0
    for _ in range(paths):
        u = rng.standard_normal(CAL + (n_per + 1) * m)
        u[:NT] /= 8.0  # quiet tau part: the first window is blind, a constant law from then on
        g = GuardedPeriodicVolatilityMonitor(threshold=A, delta=0.05, guard_tries=1, guard_fallback="rank",
                                             rank_seed=rng, period=m, **SMALL)
        g.run(u)
        if not g.segments or g.segments[0]["guard"] != "rank fallback":
            continue
        falls += 1
        s = g.segments[0]["start"]
        a = next((e.index for e in g.events if e.monitored_since == s), None)
        b0 = s + (m - s % m)  # first boundary strictly after s (P1 restarts at t > start)
        bounds = [s] + [b0 + i * m for i in range(n_per)]
        for i in range(n_per):
            lo, hi = bounds[i], bounds[i + 1]
            if hi > len(u) or lo == hi:
                continue
            if a is None or a >= lo:
                alive[i] += 1
            if a is not None and lo <= a < hi:
                hits[i] += 1
    print(f"6. rank branch, {falls} of {paths} paths fell back; per-period alarm rate (marginal, all fallback paths): "
          + ", ".join(f"{h / falls:.3f}" for h in hits) + f"; bound m / A = {m / A:.3f}")


if __name__ == "__main__":
    check_budgets()
    check_rank_restart()
    check_identity_rank_k1()
    check_delta_attribute()
    mc_rank_period()
