"""Round 9 Monte Carlo: false alarm rates of the streaming classes under the null, as a user calls them.

Every path is a fresh monitor fed one return at a time with update(). [SYNTHETIC] throughout.
Jobs (python r9_mc_null.py <job> [paths]):
  iid      : PeriodicVolatilityMonitor(baseline="known") (period=None, so VolatilityMonitor.update), i.i.d.
             Gaussian and Student t3 returns, (m, A) = (390, 7800) and (1440, 28800). Bound m / A = 5%.
  bold     : bold play (the adaptive two-point null law of E16 with exact landing) against the same streaming
             monitor, default grid and a single bet per side (0.5 up, 0.95 of 1/(B-1) down). Bound 5%.
  periodic : PeriodicVolatilityMonitor(baseline="known", period=390): per-day alarm rates over 5 days for
             i.i.d. Gaussian returns and for bold play (single bets) on every day. Bound 5% per day given the past.
  calm     : the stream that is calm (Y = 1) for 9 days and Gaussian on day 10, with and without period=390.
  cs       : VolatilityMonitor(baseline="cs", cal_len=1440, threshold=28800, delta=0.01), i.i.d. Gaussian and t3:
             P(false alarm within the first 1440 monitored steps). Bound delta + m / A = 6%.
"""
import csv
import math
import os
import sys
import tempfile
import time

import numpy as np
from scipy import integrate, stats

from avmonitor import VolatilityMonitor
from avmonitor.periodic import PeriodicVolatilityMonitor
from avmonitor.robustness import trunc_mean_normal

OUT = os.environ.get("R9_OUT", tempfile.gettempdir())  # where mc_<job>.csv is written
LAND = 1e-9


def t3_trunc_mean(tau):
    return integrate.quad(lambda v: stats.t.sf(math.sqrt(v), 3) * 2 if v > 0 else 1.0, 0, tau * tau, limit=200)[0]


TAU_G = 3.0
M0_G = float(trunc_mean_normal(TAU_G))
TAU_T = 3.0 * stats.t.ppf(0.75, 3) / stats.norm.ppf(0.75)  # 3 x the population scaled MAD of t3
M0_T = t3_trunc_mean(TAU_T)


def draw(law, rng, n):
    return rng.standard_normal(n) if law == "gauss" else rng.standard_t(3, n)


def row(job, case, m, A, P, k, bound):
    p = k / P
    se = math.sqrt(max(p * (1 - p), 1e-300) / P)
    d = dict(job=job, case=case, m=m, A=A, paths=P, alarms=k, rate=p, se=se, bound=bound,
             z_vs_bound=(p - bound) / se if se > 0 else float("nan"))
    print(f"[SYNTHETIC] {job:9s} {case:40s} m={m:5d} A={A:7.0f} paths={P:6d} rate={p:.4%} (se {se:.4%}) bound {bound:.2%}",
          flush=True)
    return d


def known(law, A, **kw):
    tau, m0 = (TAU_G, M0_G) if law == "gauss" else (TAU_T, M0_T)
    return PeriodicVolatilityMonitor(threshold=A, baseline="known", tau=tau, m0=m0, **kw), tau, m0


def bold_y(mon, A, B, u, reset=False):
    """One draw of E16's bold play against the monitor's current state (Shiryaev-Roberts, two sides).

    ``reset=True``: the next update() starts with a period boundary, so play against statistics at 0."""
    tgt = 2 * A * (1 + LAND)
    S_up = np.zeros_like(mon.S_up) if reset else mon.S_up
    S_dn = np.zeros_like(mon.S_dn) if reset else mon.S_dn
    gu = S_up + 1.0
    au, bu = gu.mean(), (gu * mon.up_grid).mean()
    yh = B
    if bu > 0:
        ys = 1.0 + (tgt - au) / bu
        if ys <= B:
            yh = max(ys, 1.0)
    yl = 0.0
    if mon.down_active:
        gd = S_dn + 1.0
        ad, bd = gd.mean(), (gd * mon.lam_dn).mean()
        if bd > 0:
            ys = 1.0 - (tgt - ad) / bd
            if ys >= 0:
                yl = min(ys, 1.0)
    span = yh - yl
    p = (1.0 - yl) / span if span > 0 else 0.0
    return yh if u < p else yl


def job_iid(P):
    rows = []
    for m, A in ((390, 7800.0), (1440, 28800.0)):
        for law in ("gauss", "t3"):
            rng = np.random.default_rng(91 + m)
            k = 0
            for i in range(P if m == 390 else P // 2):
                mon, _, _ = known(law, A)
                for v in draw(law, rng, m):
                    if mon.update(v) is not None:
                        k += 1
                        break
            rows.append(row("iid", f"{law}, known baseline", m, A, P if m == 390 else P // 2, k, m / A))
    return rows


def job_bold(P):
    rows = []
    m, A = 390, 7800.0
    for name, kw in (("default grid", {}), ("single bets 0.5 / 0.95", dict(up_grid=(0.5,), down_fractions=(0.95,)))):
        rng = np.random.default_rng(7)
        k = 0
        for i in range(P):
            mon, tau, m0 = known("gauss", A, **kw)
            B = tau * tau / m0
            u = rng.random(m)
            for t in range(m):
                y = bold_y(mon, A, B, u[t])
                if mon.update(math.sqrt(y * m0)) is not None:
                    k += 1
                    break
        rows.append(row("bold", f"bold play, {name}", m, A, P, k, m / A))
    return rows


def job_periodic(P):
    rows = []
    m, A, D = 390, 7800.0, 5
    for name in ("gauss iid", "bold play, single bets"):
        rng = np.random.default_rng(8)
        hits = np.zeros((P, D), bool)
        for i in range(P):
            if name == "gauss iid":
                mon, tau, m0 = known("gauss", A, period=m)
                r = rng.standard_normal(m * D)
                for t in range(m * D):
                    if mon.update(r[t]) is not None:
                        hits[i, t // m] = True
            else:
                mon, tau, m0 = known("gauss", A, period=m, up_grid=(0.5,), down_fractions=(0.95,))
                B = tau * tau / m0
                u = rng.random(m * D)
                for t in range(m * D):
                    # update() first restarts the statistics at a period boundary; the adversary plays against that
                    reset = t % m == 0 and t > mon.start
                    y = bold_y(mon, A, B, u[t], reset)
                    if mon.update(math.sqrt(y * m0)) is not None:
                        hits[i, t // m] = True
        for d in range(D):
            rows.append(row("periodic", f"{name}, P(alarm on day {d + 1})", m, A, P, int(hits[:, d].sum()), m / A))
    return rows


def job_calm(P):
    rows = []
    m, A, d = 390, 7800.0, 10
    rng = np.random.default_rng(9)
    k_plain = k_rest = 0
    calm = np.full((d - 1) * m, math.sqrt(M0_G))
    for i in range(P):
        r = np.concatenate([calm, rng.standard_normal(m)])
        a, _, _ = known("gauss", A)
        b, _, _ = known("gauss", A, period=m)
        ea = [e.index for e in a.run(r)]
        eb = [e.index for e in b.run(r)]
        assert all(e >= (d - 1) * m for e in ea + eb)
        k_plain += any(e >= (d - 1) * m for e in ea)
        k_rest += any(e >= (d - 1) * m for e in eb)
    rows.append(row("calm", "calm 9 days, Gaussian day 10, no period", m, A, P, k_plain, m / A))
    rows.append(row("calm", "calm 9 days, Gaussian day 10, period=390", m, A, P, k_rest, m / A))
    return rows


def job_cs(P):
    rows = []
    m, A, cal = 1440, 28800.0, 1440
    for law in ("gauss", "t3"):
        rng = np.random.default_rng(10)
        k = 0
        for i in range(P):
            mon = VolatilityMonitor(threshold=A, cal_len=cal, baseline="cs", delta=0.01)
            r = draw(law, rng, cal + m)
            for v in r:
                ev = mon.update(v)
                if ev is not None:
                    assert ev.index < cal + m
                    k += 1
                    break
            assert mon.segments and mon.segments[0]["start"] == cal
        rows.append(row("cs", f"{law}, cs baseline (delta 0.01)", m, A, P, k, 0.01 + m / A))
    return rows


if __name__ == "__main__":
    job = sys.argv[1]
    P = int(sys.argv[2])
    t0 = time.time()
    rows = {"iid": job_iid, "bold": job_bold, "periodic": job_periodic, "calm": job_calm, "cs": job_cs}[job](P)
    with open(f"{OUT}/mc_{job}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{job}: {time.time() - t0:.0f} s", flush=True)
