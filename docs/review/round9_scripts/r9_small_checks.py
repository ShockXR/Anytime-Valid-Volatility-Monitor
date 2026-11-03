"""Round 9 small deterministic checks (each prints PASS/FAIL style lines with numbers).

1. e-BH (panel.ebh) against a brute-force implementation of Wang and Ramdas (2022).
2. Predictability: the bet at step t does not change when Y_t changes (adaptive bet of betting_edetector,
   both bets of eb_confseq).
3. eb_confseq (vectorised) against StreamingEBConfSeq (the monitor's), and calibrate(baseline="cs") against
   the monitor's cs segment.
4. Boundaries of Y: r = 0 (Y = 0) and r = tau (Y = B) in the streaming monitor; every factor is positive
   and the decrease factor at Y = B is 1 - fraction.
5. Off by one: the calm law Y = 1 makes the known-baseline monitor alarm at the 2A-th monitored step
   (0-based index 2A - 1), and never with period < 2A.
6. Input validation: bets that make a factor negative are accepted by VolatilityMonitor and by
   betting_edetector's lam_cap / adaptive_weight.
7. two_sided_alarm with alarm arrays of different integer dtypes.
8. PeriodicVolatilityMonitor(baseline="known", c_down=...) silently ignores c_down.
9. running.two_sided_sr with one missing value.
10. design.sr_steps_constant_factor at exact boundaries, and decrease_speed_limit at the reported settings.
"""
import itertools
import math

import numpy as np

from avmonitor import VolatilityMonitor, betting_edetector, calibrate, two_sided_alarm
from avmonitor.confseq import StreamingEBConfSeq, eb_confseq
from avmonitor.detectors import DetectorResult
from avmonitor.panel import ebh
from avmonitor.periodic import PeriodicVolatilityMonitor
from avmonitor.robustness import trunc_mean_normal

M0 = float(trunc_mean_normal(3.0))


def ebh_brute(e, alpha):
    K = len(e)
    best = 0
    for r in range(1, K + 1):
        if np.sort(e)[::-1][r - 1] >= K / (alpha * r):
            best = r
    if best == 0:
        return np.zeros(K, bool)
    order = np.argsort(-e, kind="stable")
    mask = np.zeros(K, bool)
    mask[order[:best]] = True
    return mask


def check_ebh():
    """Mismatches are counted per kind of input: 0 continuous values, 1 values computed as K / (alpha r),
    which sit exactly on the e-BH boundaries, so floating-point rounding decides; 2 ties, zeros and inf."""
    rng = np.random.default_rng(1)
    bad = {0: 0, 1: 0, 2: 0}
    cnt = {0: 0, 1: 0, 2: 0}
    n = 0
    for trial in range(20000):
        K = int(rng.integers(1, 12))
        alpha = float(rng.choice([0.05, 0.1, 0.2]))
        kind = trial % 3
        if kind == 0:
            e = rng.exponential(1.0, K) * rng.choice([1, 10, 100], K)
        elif kind == 1:  # values exactly on the e-BH boundaries K / (alpha r)
            e = np.array([K / (alpha * int(rng.integers(1, K + 1))) for _ in range(K)])
        else:  # ties and zeros
            e = rng.choice([0.0, 1.0, K / alpha, 2 * K / alpha, np.inf], K)
        with np.errstate(over="ignore"):
            a = ebh(e, alpha)
        b = ebh_brute(e, alpha)
        n += 1
        cnt[kind] += 1
        if a.sum() != b.sum() or not np.array_equal(np.sort(e[a]), np.sort(e[b])):
            bad[kind] += 1
    print(f"1. e-BH against brute force, {n} random cases: mismatches {bad[0]} of {cnt[0]} continuous, "
          f"{bad[2]} of {cnt[2]} with ties, zeros and inf, {bad[1]} of {cnt[1]} with every value computed as "
          f"K / (alpha r) (exactly on a boundary, where rounding of e r against K / alpha decides)")


def check_predictable():
    rng = np.random.default_rng(2)
    x = np.minimum(rng.standard_normal(400) ** 2, 9.0)
    base = betting_edetector(x, M0, "up", threshold=1e30, adaptive_h=20, return_path=True)
    worst = 0.0
    for t in (25, 100, 399):
        x2 = x.copy()
        x2[t] = 9.0 if x[t] < 1 else 0.0
        alt = betting_edetector(x2, M0, "up", threshold=1e30, adaptive_h=20, return_path=True)
        worst = max(worst, float(np.max(np.abs(alt.lam_path[: t + 1] - base.lam_path[: t + 1]))))
    cs0 = eb_confseq(x[None], 9.0, n_plan=300)
    worst_cs = 0.0
    for t in (0, 10, 399):
        x2 = x.copy()
        x2[t] = 9.0 - x[t]
        cs1 = eb_confseq(x2[None], 9.0, n_plan=300)
        worst_cs = max(worst_cs, float(np.max(np.abs(cs1.lam_upper[0, : t + 1] - cs0.lam_upper[0, : t + 1]))),
                       float(np.max(np.abs(cs1.lam_lower[0, : t + 1] - cs0.lam_lower[0, : t + 1]))))
    print(f"2. predictability: largest change of a bet at or before t when Y_t changes: adaptive {worst:.3g}, "
          f"confidence sequence {worst_cs:.3g}")


def check_cs_equal():
    rng = np.random.default_rng(3)
    x = np.minimum(rng.standard_t(3, 2000) ** 2, 12.0)
    cs = eb_confseq(x[None], 12.0, 0.005, 0.005, n_plan=1000)
    s = StreamingEBConfSeq(12.0, 0.005, 0.005, n_plan=1000)
    du = dl = 0.0
    for t, v in enumerate(x):
        s.update(v)
        du = max(du, abs(s.upper - cs.upper[0, t]))
        dl = max(dl, abs(s.lower - cs.lower[0, t]))
    # monitor's cs segment against calibrate(baseline="cs")
    r = rng.standard_normal(1440)
    mon = VolatilityMonitor(cal_len=1440, baseline="cs", delta=0.01)
    for v in r:
        mon.update(v)
    b = calibrate(r, baseline="cs", delta=0.01)
    seg = mon.segments[0]
    print(f"3. eb_confseq vs StreamingEBConfSeq over 2000 steps: max |dU| {du:.3g}, max |dL| {dl:.3g}; "
          f"monitor vs calibrate(cs): tau2 {seg['tau2'] - b.tau**2:.3g}, m_up {seg['m_up'] - b.m0_up:.3g}, "
          f"m_dn {seg['m_dn'] - b.m0_down:.3g}")


def check_boundaries():
    mon = PeriodicVolatilityMonitor(threshold=1e12, baseline="known", tau=3.0, m0=M0)
    B = 9.0 / M0
    out = []
    for r in (0.0, 3.0, -3.0, 0.0, 1e6, 3.0 * (1 + 1e-15)):
        su, sd = mon.S_up.copy(), mon.S_dn.copy()
        mon.update(r)
        fu = mon.S_up / (su + 1.0)
        fd = mon.S_dn / (sd + 1.0)
        out.append((r, float(fu.min()), float(fd.min())))
    ok = all(a > 0 and b > 0 for _, a, b in out)
    at_B = [b for r, _, b in out if abs(r) >= 3.0]
    print(f"4. Y at 0 and at B = {B:.4f}: smallest factors (up, down) per step {[(r, round(a, 6), round(b, 6)) for r, a, b in out]}; "
          f"all positive: {ok}; decrease factor at Y = B: {min(at_B):.15f} (1 - 0.95 = 0.05)")


def check_calm():
    A = 50.0
    r = np.full(400, math.sqrt(M0))
    a = PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=3.0, m0=M0).run(r)
    b = PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=3.0, m0=M0, period=99).run(r)
    print(f"5. calm law Y = 1, A = {A:g}: alarms at {[e.index for e in a]} (expected 99, 199, 299, 399), "
          f"with period 99 < 2A: {len(b)} alarms")


def check_validation():
    msgs = []
    try:
        mon = PeriodicVolatilityMonitor(threshold=1e9, baseline="known", tau=3.0, m0=M0, up_grid=(1.5,),
                                        down_fractions=(1.5,))
        mon.update(0.0)
        mon.update(3.0)
        msgs.append(f"VolatilityMonitor accepts up_grid=(1.5,), down_fractions=(1.5,): S_up {mon.S_up}, S_dn {mon.S_dn}")
    except ValueError as e:
        msgs.append(f"VolatilityMonitor rejects bad bets: {e}")
    x = np.r_[np.full(30, 1.5), 0.0, 0.0]  # local mean 1.5, local variance 0: aGRAPA bet 0.5 / 0.25 = 2
    try:
        d = betting_edetector(x, 1.0, "up", threshold=1e9, adaptive_h=5, lam_cap=3.0, return_path=True)
        msgs.append(f"betting_edetector accepts lam_cap=3: min statistic {np.min(d.stat_path):.3g}, "
                    f"largest adaptive bet {np.max(d.lam_path):.3g}")
    except ValueError as e:
        msgs.append(f"lam_cap rejected: {e}")
    try:
        d = betting_edetector(np.full(10, 0.0), 1.0, "up", threshold=1e9, adaptive_h=5, adaptive_weight=1.5,
                              return_path=True)
        msgs.append(f"betting_edetector accepts adaptive_weight=1.5: final statistic {d.stat_path[-1]:.3g} (negative weight)")
    except ValueError as e:
        msgs.append(f"adaptive_weight rejected: {e}")
    print("6. " + "\n   ".join(msgs))


def check_dtype():
    a = DetectorResult(np.array([-1, 5], np.int64))
    b = DetectorResult(np.array([-1, 7], np.int32))
    print(f"7. two_sided_alarm(int64 [-1, 5], int32 [-1, 7]) = {two_sided_alarm(a, b)} (expected [-1, 5])")


def check_known_cdown():
    mon = PeriodicVolatilityMonitor(threshold=1e9, baseline="known", tau=3.0, m0=M0, c_down=1.5)
    seg = mon.segments[0]
    print(f"8. known baseline with c_down=1.5: segment says c_down={seg['c_down']}, but tau2_dn={seg['tau2_dn']} "
          f"(= tau^2 = 9; c_down is ignored without an error)")


def check_two_sided_sr():
    from avmonitor.running import two_sided_sr

    rng = np.random.default_rng(0)
    x = np.minimum((3 * rng.standard_normal(500)) ** 2, 9.0)
    a = two_sided_sr(x[None], 9.0, M0, M0, [100.0])[100.0]
    x2 = x.copy()
    x2[0] = np.nan
    b = two_sided_sr(x2[None], 9.0, M0, M0, [100.0])[100.0]
    print(f"9. two_sided_sr first (up, down) alarm: complete data ({int(a[0][0])}, {int(a[1][0])}); "
          f"with x[0] = NaN ({int(b[0][0])}, {int(b[1][0])})")


def check_speed_limit():
    from avmonitor.design import decrease_speed_limit, sr_steps_constant_factor
    from avmonitor.robustness import gaussian_B

    bad = tot = 0
    for F in (1.01, 1.125, 2.0, 1.5, 10 / 9):
        for n in range(1, 40):
            A = F * (F**n - 1.0) / (F - 1.0)  # the docstring's closed form M_n equals A exactly
            tot += 1
            bad += sr_steps_constant_factor(A, F) != n
    B = gaussian_B()
    rows = []
    for A in (300.0, 7800.0, 15600.0, 28800.0, 57600.0):
        F, M, k = B / (B - 1.0), 0.0, 0
        while M < A:
            M, k = F * (M + 1.0), k + 1
        rows.append((A, decrease_speed_limit(A, B), k))
    print(f"10. sr_steps_constant_factor one step too large at {bad} of {tot} exact boundaries; "
          f"decrease_speed_limit vs brute force at B = {B:.3f}: {rows}")


if __name__ == "__main__":
    check_ebh()
    check_predictable()
    check_cs_equal()
    check_boundaries()
    check_calm()
    check_validation()
    check_dtype()
    check_known_cdown()
    check_two_sided_sr()
    check_speed_limit()
