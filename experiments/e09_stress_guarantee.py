"""E9. Stress test of the false alarm guarantee at scale, and runtime profiling.

All data here are SYNTHETIC. Seeds are fixed.

Research questions
  Q9a. Tightness. With the true baseline, Theorem B (docs/math.md, 5.4) says the
       two-sided monitor (each side at level 2A) raises a false alarm within m
       steps with probability at most m / A. How close to the bound is the actual
       rate, for horizons from a quarter of a day to a week?
  Q9b. Distribution-free validity. The theorem needs only E[X_t | past] = m0 and
       a bounded X. Does the empirical rate stay below m / A for heavy tails with
       infinite fourth moment (Student t3), rare large jumps, and coarse price
       ticks (a tick of half a standard deviation, about 20% zero returns)?
  Q9c. Serial dependence, where the theorem's null is FALSE. Under GARCH(1,1)
       volatility clustering the unconditional mean of X equals m0 but the
       conditional mean moves around it, so the guarantee does not apply. How
       large must the null band b be to bring the one-day false alarm rate back to
       the 5% budget, as a function of persistence (alpha + beta)?
  Q9d. Profiling. Time per step of the streaming monitor and of the vectorised
       detector, and of the kernel estimators, so users know what scales.
"""

import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import robust as R
from avmonitor.detectors import betting_edetector, two_sided_alarm
from avmonitor.monitor import VolatilityMonitor

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
SEED = 20261103


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return c - h, c + h


# ----------------------------------------------------------------------------
# Return generators with a known truncation level and a known null mean
# ----------------------------------------------------------------------------


def draw(kind, size, rng):
    if kind == "gaussian":
        return rng.standard_normal(size)
    if kind == "student_t3":
        return rng.standard_t(3, size) / np.sqrt(3.0)  # unit variance, infinite fourth moment
    if kind == "jumps":
        r = rng.standard_normal(size)
        return r + (rng.random(size) < 0.005) * rng.normal(0, 10, size)
    if kind == "ticks":
        return np.round(rng.standard_normal(size) / 0.5) * 0.5  # tick of 0.5 sd: about 20% zero returns
    raise ValueError(kind)


def truth(kind, rng, n=20_000_000):
    """tau = 3 x scaled MAD of the distribution (floored at 3 ticks for 'ticks') and m0 = E[min(r^2, tau^2)]."""
    big = draw(kind, n, rng)
    scale = np.median(np.abs(big)) * R.MAD_NORMAL_SCALE
    tau = 3 * scale
    if kind == "ticks":
        tau = max(tau, 3 * 0.5)
    x = np.minimum(big**2, tau**2)
    return tau, float(x.mean()), float(x.std() / np.sqrt(n)), float(np.mean(big == 0))


def two_sided(x, m0, tau2, A, band=0.0):
    up = betting_edetector(x, m0 * (1 + band), "up", threshold=2 * A)
    dn = betting_edetector(x, m0 * (1 - band), "down", upper=tau2, threshold=2 * A)
    return two_sided_alarm(up, dn)


def q9ab(reps=4000):
    rng = np.random.default_rng(SEED)
    A = 7800.0  # 390 / 0.05: at most 5% chance of a false alarm in a 390-step day
    horizons = [98, 390, 1950, 3900, 7800]
    n = max(horizons)
    rows = []
    for kind in ["gaussian", "student_t3", "jumps", "ticks"]:
        tau, m0, se, zero = truth(kind, rng)
        alarms = []
        for _ in range(reps // 1000):
            r = draw(kind, (1000, n), rng)
            x = np.minimum(r**2, tau**2)
            alarms.append(two_sided(x, m0, tau**2, A))
        a = np.concatenate(alarms)
        for m in horizons:
            k = int(np.sum((a >= 0) & (a < m)))
            lo, hi = wilson(k, a.size)
            rows.append(dict(distribution=kind, tau=tau, m0=m0, m0_mc_se=se, zero_share=zero, A=A, horizon=m,
                             bound=m / A, false_alarm_rate=k / a.size, ci_low=lo, ci_high=hi, paths=a.size,
                             rate_over_bound=(k / a.size) / (m / A)))
        print(kind, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e09ab_bound_tightness.csv", index=False)
    print(df.pivot(index="horizon", columns="distribution", values="false_alarm_rate").round(4).to_string())
    print("bound m/A:", [round(m / A, 4) for m in horizons])
    return df


def garch(reps, n, alpha, beta, rng, burn=5000):
    omega = 1 - alpha - beta  # unconditional variance 1
    s2 = np.ones(reps)
    r_prev = rng.standard_normal(reps)
    out = np.empty((reps, n))
    for t in range(burn + n):
        s2 = omega + alpha * r_prev**2 + beta * s2
        r_prev = np.sqrt(s2) * rng.standard_normal(reps)
        if t >= burn:
            out[:, t - burn] = r_prev
    return out


def q9c(reps=2000):
    rng = np.random.default_rng(SEED + 1)
    A = 28800.0  # 1440 / 0.05, minute data
    n = 1440
    rows = []
    for alpha, beta in [(0.05, 0.90), (0.05, 0.94), (0.03, 0.969), (0.02, 0.9799)]:
        pers = round(alpha + beta, 6)
        long = garch(200, 50_000, alpha, beta, rng)
        scale = np.median(np.abs(long)) * R.MAD_NORMAL_SCALE
        tau = 3 * scale
        m0 = float(np.minimum(long**2, tau**2).mean())
        r = garch(reps, n, alpha, beta, rng)
        x = np.minimum(r**2, tau**2)
        for band in [0.0, 0.1, 0.25, 0.5, 1.0]:
            up = betting_edetector(x, m0 * (1 + band), "up", threshold=2 * A)
            if band < 1.0:
                dn = betting_edetector(x, m0 * (1 - band), "down", upper=tau**2, threshold=2 * A)
                a = two_sided_alarm(up, dn)
                dn_rate = float(np.mean(dn.alarm >= 0))
            else:
                a, dn_rate = up.alarm, 0.0
            k = int(np.sum(a >= 0))
            lo, hi = wilson(k, reps)
            rows.append(dict(alpha=alpha, beta=beta, persistence=pers, half_life_steps=np.log(0.5) / np.log(pers),
                             band=band, one_day_false_alarm=k / reps, ci_low=lo, ci_high=hi,
                             up_side=float(np.mean(up.alarm >= 0)), down_side=dn_rate, budget=n / A))
        print("garch", pers, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e09c_garch_band.csv", index=False)
    print(df.pivot(index="band", columns="persistence", values="one_day_false_alarm").round(3).to_string())
    return df


def q9d():
    rng = np.random.default_rng(SEED + 2)
    rows = []
    r = rng.standard_normal(200_000)
    t = time.perf_counter()
    VolatilityMonitor(threshold=28800, cal_len=1440, null_band=0.25).run(r)
    rows.append(dict(task="VolatilityMonitor.run (streaming, one series)", size=r.size,
                     seconds=time.perf_counter() - t))
    for reps in [1, 100, 1000, 10000]:
        x = np.minimum(rng.standard_normal((reps, 2000)) ** 2, 9.0)
        t = time.perf_counter()
        betting_edetector(x, 0.9946, "up", threshold=1e12)
        rows.append(dict(task=f"betting_edetector, {reps} paths x 2000 steps", size=reps * 2000,
                         seconds=time.perf_counter() - t))
    for n in [390, 1440, 23400]:
        x = rng.standard_normal((200, n))
        t = time.perf_counter()
        R.kernel_iq(x, 30, c=3)
        rows.append(dict(task=f"kernel_iq robust, 200 days of {n} returns", size=200 * n, seconds=time.perf_counter() - t))
        t = time.perf_counter()
        R.med_rq(x)
        rows.append(dict(task=f"med_rq, 200 days of {n} returns", size=200 * n, seconds=time.perf_counter() - t))
    df = pd.DataFrame(rows)
    df["microseconds_per_unit"] = df.seconds / df["size"] * 1e6
    df["machine"] = f"{platform.processor() or platform.machine()}, Python {platform.python_version()}, shared with other jobs"
    df.to_csv(OUT / "e09d_profile.csv", index=False)
    print(df[["task", "size", "seconds", "microseconds_per_unit"]].round(3).to_string(index=False))


def figure():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ab = pd.read_csv(OUT / "e09ab_bound_tightness.csv")
    c = pd.read_csv(OUT / "e09c_garch_band.csv")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for kind, g in ab.groupby("distribution"):
        ax[0].errorbar(g.horizon, g.false_alarm_rate, yerr=[g.false_alarm_rate - g.ci_low, g.ci_high - g.false_alarm_rate],
                       marker="o", ms=3, lw=1, capsize=2, label=kind)
    hs = np.sort(ab.horizon.unique())
    ax[0].plot(hs, hs / ab.A.iloc[0], "k--", lw=1, label="bound m / A")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("horizon m (steps)")
    ax[0].set_ylabel("P(false alarm within m)")
    ax[0].set_title("Q9a, Q9b [synthetic]: true baseline, i.i.d. returns")
    ax[0].legend(fontsize=8)
    for p, g in c.groupby("persistence"):
        ax[1].plot(g.band, g.one_day_false_alarm, marker="o", label=f"alpha + beta = {p:g}")
    ax[1].axhline(0.05, color="k", ls="--", lw=1, label="5% budget")
    ax[1].set_xlabel("null band b")
    ax[1].set_ylabel("one-day false alarm rate")
    ax[1].set_title("Q9c [synthetic]: GARCH clustering (null false)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "figures/e09_stress.png", dpi=130)


def main(parts=("ab", "c", "d")):
    if "ab" in parts:
        q9ab()
    if "c" in parts:
        q9c()
    if "d" in parts:
        q9d()
    figure()


if __name__ == "__main__":
    import sys

    main(tuple(sys.argv[1:]) or ("ab", "c", "d"))
