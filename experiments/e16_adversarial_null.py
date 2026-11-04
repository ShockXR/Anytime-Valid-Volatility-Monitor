"""E16. Adversarial null search: how close can a legitimate adversary push the default two-sided
detector to its false alarm bound, and where exactly does it break outside the theorem's conditions?

All results are SYNTHETIC. Seeds are fixed (listed below). Theory: docs/math.md, sections 5 and 6;
docs/theory_delay.md, section 8 (why the bound looks loose); docs/theory_running_tau.md,
Proposition A.6 (the plug-in baseline has no bound of the form C m / A).

The detector is the monitor's default: X_t = min(r_t^2, tau^2), increase bets 0.02 to 0.5, decrease
bets 0.1 to 0.95 times 1 / (B - 1) with B = tau^2 / m0, Shiryaev-Roberts statistics, each side
alarming at 2A, A = 7800 (so P(false alarm within m steps) <= m / A: 5% for one 390-step day and 25%
for five days, Theorem B, whenever E[X_t | F_{t-1}] = m0 with the true m0).

Part a: inside the null class (true baseline, conditional mean exactly m0)
  a1  i.i.d. laws. Searched on Y = X / m0 in [0, B] with mean 1: two-point and three-point laws
      (B fixed at the Gaussian value 9.045, and B free in [2, 40]); and return laws with tau = 3 x
      scaled MAD and m0 computed from the law itself: scale mixtures of two normals, Student t,
      returns rounded to a tick. Grid search, then Nelder-Mead, with common random numbers (CRN);
      every reported rate comes from fresh seeds.
  a2  Dependent laws. (i) Exogenous shape switching: a hidden Markov chain switches the shape of
      the law between a two-point law {0, B} and the Gaussian shape, the conditional mean staying
      at m0 (grid search). (ii) Adaptive laws: at each step the adversary sees the detector's
      state and picks a two-point law with mean 1 (two points suffice at each step, docs/
      theory_delay.md 8.3). "Bold play": spread {0, B}, except that when one draw can take a side
      exactly to 2A it lands there (no overshoot). Policy search over an 8-parameter family that
      makes the spread depend on both statistics and on the time left (random search, then
      Nelder-Mead, CRN), and a "hold, then bold" family. Fresh-seed evaluation at one and five days,
      the whole horizon curve, and a scan over B.
  a3  Where the rest of the slack sits: the optional-stopping decomposition of
      docs/theory_delay.md 8.1 under the worst adversary (stopping, wealth stranded on paths that
      never alarm, by bet, and overshoot); and a check of the bold policy on single bets against
      E11's dynamic programme (results/e11d_worst_case_dp.csv; its values are lower bounds on the
      worst case).
  a4  B below the searched range (B = 1.2 to 2): bold play and the i.i.d. law {0, B}. Under the theorem
      tau is any level fixed in advance, so B can be any number above 1.
  a5  B above the searched range (B = 100, 150 and 1,000): bold play. When more than half of the
      returns are zero, the monitor's fallback for tau (median of the nonzero |r|) allows B above
      39.6; the price ticks of b3 reach B = 149.
  a6  What can be proved, next to what was found: (i) the floor bound. On a path without an alarm
      by step m, N_m = (M_up + M_down) / 2 is at least c_m, its value when every draw goes against
      every bet (Y = 0 for the increase bets, Y = B for the decrease bets), so optional stopping
      gives P(false alarm within m) <= (m - c_m) / (A - c_m) for every null law. (ii) E11's dynamic
      programme for the single bet 0.05 on finer wealth grids (its values are lower bounds that
      rise as the grid is refined).

Part b: outside the class (do not read these as failures of the theorem)
  b1  GARCH(1,1) returns, unit unconditional variance, tau and m0 at their unconditional values:
      grid over persistence and alpha. Failure mapped against the size of the conditional mean
      drift E|mu_t / m0 - 1| (mu_t = E[X_t | F_{t-1}], exact for Gaussian innovations), and
      checked by the same detector betting against the true mu_t (inside the class again).
  b2  Plug-in baselines, calibration 390 and 1440, i.i.d. returns (the null holds exactly):
      Gaussian, scale mixtures, Student t, return ticks; three detectors on the same paths: plug-in
      (tau-hat, m-hat), the in-class "tau-hat oracle" (tau-hat with the true mean m(tau-hat)) and the
      full oracle (true tau and m0). Plus an adversary inside the null class that only changes the
      shape of the law after calibration (constant |r| after calibration, same conditional mean).
  b3  Price-level ticks: observed price = tick x round(efficient price / tick), so returns are
      dependent (bid-ask-bounce-like). Tick sizes 0.05 to 6 sd, true unconditional tau and m0.
      In-class control: i.i.d. draws from the marginal law of the price-rounded returns (pooled
      from the long paths that give tau and m0), with the same tau and m0, so only the dependence
      is removed. The older "twin" rounds each return independently to the same tick: it is
      inside the class too, but it is a different law (far more zero returns, a much larger B),
      so it checks the tick grid, not the values. Conditional mean drift measured as
      E[X_t | previous return] / m0.

Usage (each part writes its own CSV files; run with at most 2 processes on a shared machine):
  python experiments/e16_adversarial_null.py a1
  python experiments/e16_adversarial_null.py a2
  python experiments/e16_adversarial_null.py a3
  python experiments/e16_adversarial_null.py a4
  python experiments/e16_adversarial_null.py a5
  python experiments/e16_adversarial_null.py a6
  python experiments/e16_adversarial_null.py b1
  python experiments/e16_adversarial_null.py b2
  python experiments/e16_adversarial_null.py b3
  python experiments/e16_adversarial_null.py figures
  python experiments/e16_adversarial_null.py a        # a1 a2 a3 a4 a5 a6
  python experiments/e16_adversarial_null.py b        # b1 b2 b3

Seeds (SEED = 16000): search (common random numbers) SEED + 101..135; fresh evaluation SEED + 900..988 (part a)
and SEED + 991..1016 (part b2); part a5 SEED + 600..602; part b1 SEED + 201..264; b2 search SEED + 300..436;
b3 SEED + 400..410 (price-rounded paths and the twin) and SEED + 500..510 (the control with the same marginal
law); a6 draws nothing. No reported rate uses the seeds its search used. Overlaps: b3 reuses the seeds of the
b2 search with a 1440-step calibration (SEED + 400..410, laws 0 to 10); the two parts are unrelated, so this
only makes them dependent. SEED + 911 serves a1 (ticks, five days) and a2 (shape switching, one day), and
SEED + 921 and 925 serve all three adaptive families of a2, which makes those comparisons paired.
"""

from __future__ import annotations

import ast
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from avmonitor.robust import MAD_NORMAL_SCALE

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("E16_OUT", ROOT / "results"))  # E16_OUT and E16_SCALE: smoke runs only
FIG = OUT / "figures"
SCALE = float(os.environ.get("E16_SCALE", "1"))
SEED = 16000
A = 7800.0
THR = 2 * A  # each side alarms at 2A
DAY, WEEK = 390, 1950
HORIZONS = (DAY, WEEK)
C = 3.0
UPG = np.asarray(DEFAULT_UP_GRID, float)
DNF = np.asarray(DEFAULT_DOWN_FRACTIONS, float)
LAND = 1e-9  # bold play aims at 2A (1 + LAND), so floating point never falls just short


def _n(P):
    """Number of paths (or evaluations), scaled down only in smoke runs."""
    return int(max(200, round(P * SCALE))) if SCALE != 1 else int(P)


def _f(k):
    return int(max(3, round(k * SCALE))) if SCALE != 1 else int(k)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(c - h, 0.0), min(c + h, 1.0)


def rate_row(alarm, m):
    k = int(np.sum((alarm >= 0) & (alarm < m)))
    n = alarm.size
    lo, hi = wilson(k, n)
    return dict(horizon=m, bound=m / A, rate=k / n, ci_low=lo, ci_high=hi, paths=n, rate_over_bound=(k / n) / (m / A),
                above_bound=bool(lo > m / A))


# ----------------------------------------------------------------------------
# Gaussian helpers and the true baseline
# ----------------------------------------------------------------------------


def trunc_mean_normal(k):
    """E[min(Z^2, k^2)] for a standard normal Z."""
    k = np.asarray(k, float)
    return 2 * stats.norm.cdf(k) - 1 - 2 * k * stats.norm.pdf(k) + 2 * k**2 * stats.norm.sf(k)


TAU_G = C * MAD_NORMAL_SCALE * stats.norm.ppf(0.75)  # = 3 for the standard normal
B_GAUSS = float(TAU_G**2 / trunc_mean_normal(TAU_G))  # 9.045


# ----------------------------------------------------------------------------
# The detector, vectorised
# ----------------------------------------------------------------------------


def first_alarm(X, m_up, m_dn, tau2, thr=THR, up=UPG, dnf=DNF):
    """Two-sided default detector on rows of X (clipped squared returns, X <= tau2).

    m_up, m_dn: the baselines the two sides bet against (scalar, one per path, or one per path and
    step for a time-varying baseline). tau2: truncation level squared (scalar or one per path); it
    sets the decrease bets dnf / (tau2 / m_dn - 1). Same recursions and the same order of checks as
    VolatilityMonitor.update (increase side first). Returns (alarm index or -1, side +1 / -1 / 0).
    """
    X = np.atleast_2d(np.asarray(X, float))
    P, n = X.shape
    Xt = np.ascontiguousarray(X.T)

    def arr(v):
        v = np.asarray(v, float)
        return np.full(P, float(v)) if v.ndim == 0 else v

    m_up, m_dn, tau2 = arr(m_up), arr(m_dn), arr(tau2)
    tv_up, tv_dn = m_up.ndim == 2, m_dn.ndim == 2
    mu_t = np.ascontiguousarray(m_up.T) if tv_up else None
    md_t = (mu_t if m_dn is m_up else np.ascontiguousarray(m_dn.T)) if tv_dn else None
    idx = np.arange(P)
    Su = np.zeros((P, up.size))
    Sd = np.zeros((P, dnf.size))
    alarm = np.full(P, -1)
    side = np.zeros(P, np.int8)
    mu_p, md_p, t2_p = (None if tv_up else m_up.copy()), (None if tv_dn else m_dn.copy()), tau2.copy()
    lamd = None if tv_dn else dnf[None, :] / (t2_p / md_p - 1.0)[:, None]
    for t in range(n):
        x = Xt[t, idx]
        mu = mu_t[t, idx] if tv_up else mu_p
        md = md_t[t, idx] if tv_dn else md_p
        Su = (1.0 + up[None, :] * (x / mu - 1.0)[:, None]) * (Su + 1.0)
        ld = dnf[None, :] / (t2_p / md - 1.0)[:, None] if tv_dn else lamd
        Sd = (1.0 + ld * (1.0 - x / md)[:, None]) * (Sd + 1.0)
        hu = Su.mean(axis=1) >= thr
        hd = (~hu) & (Sd.mean(axis=1) >= thr)
        hit = hu | hd
        if hit.any():
            alarm[idx[hit]] = t
            side[idx[hu]] = 1
            side[idx[hd]] = -1
            keep = ~hit
            idx, Su, Sd, t2_p = idx[keep], Su[keep], Sd[keep], t2_p[keep]
            if not tv_up:
                mu_p = mu_p[keep]
            if not tv_dn:
                md_p = md_p[keep]
                lamd = lamd[keep]
            if idx.size == 0:
                break
    return alarm, side


def engine(y_fn, P, n, B, thr=THR, up=UPG, dnf=DNF, record=False, decompose=False):
    """The same detector on Y = X / m0 in [0, B] (true baseline m0, so m_up = m_dn = 1, tau2 = B).

    ``y_fn(t, idx, Su, Sd, lamd)`` returns the draws Y_t for the paths still running (``idx``); it
    may use the detector's state (an adaptive adversary). Returns a dict with alarm, side and,
    with ``record``, the Y paths; with ``decompose``, what the optional-stopping identity needs
    (N = (M_up + M_down) / 2 at the alarm, the survivors' end states) and ``m_alarm``, the
    statistic of the side that alarmed.
    """
    lamd = dnf / (B - 1.0)
    idx = np.arange(P)
    Su = np.zeros((P, up.size))
    Sd = np.zeros((P, dnf.size))
    alarm = np.full(P, -1)
    side = np.zeros(P, np.int8)
    Yrec = np.full((P, n), np.nan) if record else None
    n_at_alarm = np.full(P, np.nan) if decompose else None
    m_alarm = np.full(P, np.nan) if decompose else None
    for t in range(n):
        y = y_fn(t, idx, Su, Sd, lamd)
        if record:
            Yrec[idx, t] = y
        Su = (1.0 + up[None, :] * (y - 1.0)[:, None]) * (Su + 1.0)
        Sd = (1.0 + lamd[None, :] * (1.0 - y)[:, None]) * (Sd + 1.0)
        Mu = Su.mean(axis=1) if up.size else np.zeros(idx.size)
        Md = Sd.mean(axis=1) if dnf.size else np.zeros(idx.size)
        hu = Mu >= thr
        hd = (~hu) & (Md >= thr)
        hit = hu | hd
        if hit.any():
            alarm[idx[hit]] = t
            side[idx[hu]] = 1
            side[idx[hd]] = -1
            if decompose:
                n_at_alarm[idx[hit]] = 0.5 * (Mu[hit] + Md[hit])
                m_alarm[idx[hit]] = np.where(hu, Mu, Md)[hit]
            keep = ~hit
            idx, Su, Sd = idx[keep], Su[keep], Sd[keep]
            if idx.size == 0:
                break
    out = dict(alarm=alarm, side=side, Y=Yrec)
    if decompose:
        out.update(n_at_alarm=n_at_alarm, m_alarm=m_alarm, survivors=idx, Su_end=Su, Sd_end=Sd)
    return out


# ----------------------------------------------------------------------------
# Adversaries inside the null class (all draws have conditional mean exactly 1)
# ----------------------------------------------------------------------------


def two_point_draw(u, yl, yh):
    """Draw from the two-point law on {yl, yh} with mean 1 (yl <= 1 <= yh); degenerate at 1 if yl = yh."""
    span = yh - yl
    p = np.where(span > 0, (1.0 - yl) / np.where(span > 0, span, 1.0), 0.0)
    return np.where(u < p, yh, yl)


BOLD = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])  # far law {0, B}, landing on


def bold_policy(theta, B, n, up=UPG, thr=THR, hold=0):
    """State-dependent two-point law. Far from the thresholds: high point 1 + (B - 1) fh and low
    point fl, both clipped linear functions of zu, zd (log statistics over log 2A) and the share
    of time left. Near a threshold: land exactly (high point that takes the increase side to
    2A, low point that takes the decrease side to 2A). theta = BOLD gives the spread {0, B}.
    ``hold`` steps of Y = 1 first (the "hold, then bold" family)."""
    th = np.asarray(theta, float)
    tgt = thr * (1 + LAND)
    L = np.log(thr)

    def pol(t, idx, Su, Sd, lamd):
        au = Su.mean(axis=1) + 1.0 if up.size else np.ones(idx.size)
        ad = Sd.mean(axis=1) + 1.0 if lamd.size else np.ones(idx.size)
        if t < hold:
            return np.ones(idx.size), np.ones(idx.size)
        zu, zd, tl = np.log(au) / L, np.log(ad) / L, (n - t) / n
        fh = np.clip(th[0] + th[1] * zu + th[2] * zd + th[3] * tl, 0.0, 1.0)
        fl = np.clip(th[4] + th[5] * zu + th[6] * zd + th[7] * tl, 0.0, 1.0)
        yh, yl = 1.0 + (B - 1.0) * fh, fl
        if up.size:
            bu = (Su + 1.0) @ up / up.size
            ysu = 1.0 + (tgt - au) / bu
            yh = np.where(ysu <= B, np.maximum(ysu, 1.0), yh)
        if lamd.size:
            bd = (Sd + 1.0) @ lamd / lamd.size
            ysd = 1.0 - (tgt - ad) / bd
            yl = np.where(ysd >= 0.0, np.minimum(ysd, 1.0), yl)
        return yl, yh

    return pol


def adaptive_y_fn(pol, U):
    Ut = np.ascontiguousarray(U.T)

    def y_fn(t, idx, Su, Sd, lamd):
        yl, yh = pol(t, idx, Su, Sd, lamd)
        return two_point_draw(Ut[t, idx], yl, yh)

    return y_fn


def array_y_fn(Y):
    Yt = np.ascontiguousarray(Y.T)
    return lambda t, idx, Su, Sd, lamd: Yt[t, idx]


def run_bold(theta, U, B, hold=0, **kw):
    n = U.shape[1]
    return engine(adaptive_y_fn(bold_policy(theta, B, n, hold=hold, up=kw.get("up", UPG), thr=kw.get("thr", THR)), U),
                  U.shape[0], n, B, **kw)


# --- i.i.d. families: quantile functions of Y given uniforms (CRN), plus B --------------------


def law_two_point(par, B):
    a, b = par
    ok = 0 <= a < 1 < b <= B
    return ok, (lambda U, Z: two_point_draw(U, np.full(U.shape, a), np.full(U.shape, b)))


def law_three_point(par, B):
    y1, y2, y3, p2 = par
    if not (0 <= y1 < 1 < y3 <= B and 0 <= y2 <= B and 0 <= p2 < 1):
        return False, None
    p3 = (1 - y1 - p2 * (y2 - y1)) / (y3 - y1)
    p1 = 1 - p2 - p3
    if p1 < 0 or p3 < 0:
        return False, None
    pts = np.array([y1, y2, y3])
    pr = np.array([p1, p2, p3])
    o = np.argsort(pts)
    pts, cum = pts[o], np.cumsum(pr[o])

    def q(U, Z):
        return pts[np.minimum(np.searchsorted(cum, U, side="right"), 2)]

    return True, q


def scale_mixture_truth(w, k):
    """Mixture of N(0, 1) (weight 1 - w) and N(0, k^2): tau = 3 x scaled MAD of the law, m0, B."""
    f = lambda q: (1 - w) * (2 * stats.norm.cdf(q) - 1) + w * (2 * stats.norm.cdf(q / k) - 1) - 0.5
    med = optimize.brentq(f, 1e-9, 50 * max(k, 1))
    tau = C * MAD_NORMAL_SCALE * med
    m0 = (1 - w) * trunc_mean_normal(tau) + w * k**2 * trunc_mean_normal(tau / k)
    return tau, float(m0)


def law_scale_mixture(par):
    w, k = par
    if not (0 < w < 1 and k >= 1):
        return False, None, np.nan
    tau, m0 = scale_mixture_truth(w, k)

    def q(U, Z):
        s = np.where(U < w, k, 1.0)
        return np.minimum((s * Z) ** 2, tau**2) / m0

    return True, q, tau**2 / m0


def t_trunc_mean(s, nu):
    """E[min(T^2, s^2)] for Student t with nu > 2 degrees of freedom, in closed form: T^2 / (nu + T^2) is
    Beta(1/2, nu/2), which gives E[T^2 ; |T| < s] = nu / (nu - 2) I(s^2 / (nu + s^2); 3/2, nu/2 - 1)."""
    s = np.asarray(s, float)
    return nu / (nu - 2) * special.betainc(1.5, nu / 2 - 1, s * s / (nu + s * s)) + 2 * s * s * stats.t.sf(s, nu)


def t_truth(nu):
    tau = C * MAD_NORMAL_SCALE * stats.t.ppf(0.75, nu)
    return tau, float(t_trunc_mean(tau, nu))


def law_student_t(par):
    (nu,) = par
    if not 2.0 < nu <= 200:
        return False, None, np.nan
    tau, m0 = t_truth(nu)

    def q(U, Z):
        T = stats.t.ppf(np.clip(U, 1e-12, 1 - 1e-12), nu)
        return np.minimum(T * T, tau**2) / m0

    return True, q, tau**2 / m0


def tick_law(h, kmax=None):
    """Returns rounded to a tick h (in sd units): values j h, probabilities, tau (protocol rule), m0."""
    kmax = int(np.ceil(12 / h)) + 1 if kmax is None else kmax
    j = np.arange(-kmax, kmax + 1)
    p = stats.norm.cdf((j + 0.5) * h) - stats.norm.cdf((j - 0.5) * h)
    tau = tick_tau(np.abs(j) * h, p)
    m0 = float(np.sum(p * np.minimum((j * h) ** 2, tau**2)))
    return j, p, tau, m0


def tick_tau(absval, p):
    """tau = 3 x scaled median of |r| for a discrete law; if the median is 0 (more than half of the
    returns are zero), the monitor's fallback: the scaled median of the nonzero |r|."""
    v = np.unique(absval)
    cdf = np.array([p[absval <= x].sum() for x in v])
    med = v[np.searchsorted(cdf, 0.5 - 1e-15)]
    if med == 0:
        nz = absval > 0
        pz = p[nz] / p[nz].sum()
        v2 = np.unique(absval[nz])
        cdf2 = np.array([pz[absval[nz] <= x].sum() for x in v2])
        med = v2[np.searchsorted(cdf2, 0.5 - 1e-15)]
    return float(C * MAD_NORMAL_SCALE * med)


def law_ticks(par):
    (h,) = par
    if not 0.01 <= h <= 1.3:
        return False, None, np.nan
    _, _, tau, m0 = tick_law(h)

    def q(U, Z):
        r = h * np.round(Z / h)
        return np.minimum(r * r, tau**2) / m0

    return True, q, tau**2 / m0


def gaussian_q(U, Z):
    return np.minimum(Z * Z, TAU_G**2) / trunc_mean_normal(TAU_G)


def iid_alarms(q, B, U, Z, chunk=None):
    chunk = _chunk(U.shape[1]) if chunk is None else chunk
    out = []
    for s in range(0, U.shape[0], chunk):
        Y = q(U[s : s + chunk], Z[s : s + chunk] if Z is not None else None)
        out.append(engine(array_y_fn(Y), Y.shape[0], Y.shape[1], B)["alarm"])
    return np.concatenate(out)


def crn(P, n, seed, normals=False):
    rng = np.random.default_rng(seed)
    U = rng.random((P, n), dtype=np.float32).astype(float)
    Z = rng.standard_normal((P, n), dtype=np.float32).astype(float) if normals else None
    return U, Z


def _chunk(n, budget=8e6):
    """Paths per chunk so that one (paths x steps) array holds at most ``budget`` numbers."""
    return int(max(500, min(20000, budget // n)))


def fresh_iid(q, B, P, n, seed, normals, chunk=None):
    """Fresh-seed evaluation in chunks (memory), never reusing training seeds."""
    chunk = _chunk(n) if chunk is None else chunk
    rng = np.random.default_rng(seed)
    out = []
    for s in range(0, P, chunk):
        p = min(chunk, P - s)
        U = rng.random((p, n))
        Z = rng.standard_normal((p, n)) if normals else None
        out.append(engine(array_y_fn(q(U, Z)), p, n, B)["alarm"])
    return np.concatenate(out)


def fresh_bold(theta, B, P, n, seed, hold=0, chunk=None):
    chunk = _chunk(n) if chunk is None else chunk
    rng = np.random.default_rng(seed)
    out = []
    for s in range(0, P, chunk):
        p = min(chunk, P - s)
        out.append(run_bold(theta, rng.random((p, n)), B, hold=hold)["alarm"])
    return np.concatenate(out)


# ----------------------------------------------------------------------------
# Part a1: i.i.d. laws
# ----------------------------------------------------------------------------


def _search(obj, starts, maxfev, tag, trace):
    best = None
    for x0 in starts:
        res = optimize.minimize(obj, np.asarray(x0, float), method="Nelder-Mead",
                                options=dict(maxfev=maxfev, xatol=1e-3, fatol=1e-6, adaptive=True))
        if best is None or res.fun < best.fun:
            best = res
        log(tag, "NM from", np.round(x0, 3), "->", round(-res.fun, 5), np.round(res.x, 4))
    return best


def part_a1():
    rows, trace = [], []
    specs = []  # (family, B setting, horizon, best params, train rate, B, quantile maker, normals)
    for hz in HORIZONS:
        P_train = _n(20000 if hz == DAY else 6000)
        U, Z = crn(P_train, hz, SEED + 100 + hz // DAY, normals=True)
        fe = [0]

        def make_obj(maker, Bfun, family):
            def obj(x):
                ok, q, B = maker(x)
                if not ok:
                    return 1.0
                a = iid_alarms(q, B, U, Z)
                v = float(np.mean(a >= 0))
                fe[0] += 1
                trace.append(dict(part="a1", family=family, horizon=hz, eval=fe[0], params=str(np.round(x, 5).tolist()),
                                  train_rate=v))
                return -v
            return obj

        # two-point, B fixed
        mk2 = lambda x: (*law_two_point(x, B_GAUSS), B_GAUSS)
        obj = make_obj(mk2, None, "two-point, B = 9.045")
        grid = [(a, b) for a in (0.0, 0.1, 0.25, 0.4, 0.6, 0.8) for b in np.geomspace(1.15, B_GAUSS, 9)]
        vals = [obj(g) for g in grid]
        order = np.argsort(vals)
        best = _search(obj, [grid[i] for i in order[:2]], _f(50 if hz == DAY else 25), f"two-point B fixed {hz}", trace)
        specs.append(("two-point, B = 9.045", hz, best.x, -best.fun, mk2, True))
        # two-point, B free: x = (a, fraction of B, log B)
        mk2f = lambda x: ((lambda B: (*law_two_point((x[0], 1 + x[1] * (B - 1)), B), B))(float(np.exp(x[2])))
                          if (np.log(2) <= x[2] <= np.log(40) and 0 < x[1] <= 1) else (False, None, np.nan))
        obj = make_obj(mk2f, None, "two-point, B free")
        grid = [(a, f, np.log(B)) for a in (0.0, 0.25) for f in (0.5, 1.0) for B in (2.0, 3.0, 5.0, 9.045, 15.0, 25.0, 39.0)]
        vals = [obj(g) for g in grid]
        order = np.argsort(vals)
        best = _search(obj, [grid[i] for i in order[:2]], _f(50 if hz == DAY else 25), f"two-point B free {hz}", trace)
        specs.append(("two-point, B free in [2, 40]", hz, best.x, -best.fun, mk2f, True))
        # three-point, B fixed: start from the best two-point law plus a middle point
        bp = [s for s in specs if s[0] == "two-point, B = 9.045" and s[1] == hz][0][2]
        mk3 = lambda x: (*law_three_point(x, B_GAUSS), B_GAUSS)
        obj = make_obj(mk3, None, "three-point, B = 9.045")
        starts = [(bp[0], 1.5, bp[1], 0.1), (0.0, 2.0, B_GAUSS * 0.999, 0.2), (0.0, 0.5, B_GAUSS * 0.999, 0.3)]
        best = _search(obj, starts, _f(60 if hz == DAY else 30), f"three-point {hz}", trace)
        specs.append(("three-point, B = 9.045", hz, best.x, -best.fun, mk3, True))
        # scale mixtures of two normals
        obj = make_obj(law_scale_mixture, None, "scale mixture of two normals")
        grid = [(w, k) for w in (0.01, 0.03, 0.1, 0.2, 0.35, 0.5) for k in (1.5, 2.5, 4.0, 7.0, 12.0)]
        vals = [obj(g) for g in grid]
        order = np.argsort(vals)
        best = _search(obj, [grid[i] for i in order[:2]], _f(40 if hz == DAY else 20), f"scale mixture {hz}", trace)
        specs.append(("scale mixture of two normals", hz, best.x, -best.fun, law_scale_mixture, True))
        # Student t
        obj = make_obj(law_student_t, None, "Student t")
        grid = [(nu,) for nu in (2.1, 2.5, 3.0, 4.0, 6.0, 10.0, 30.0)]
        vals = [obj(g) for g in grid]
        i = int(np.argmin(vals))
        specs.append(("Student t", hz, np.asarray(grid[i]), -vals[i], law_student_t, False))
        # return ticks
        obj = make_obj(law_ticks, None, "returns rounded to a tick")
        grid = [(h,) for h in (0.05, 0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.2, 1.3)]
        vals = [obj(g) for g in grid]
        i = int(np.argmin(vals))
        specs.append(("returns rounded to a tick", hz, np.asarray(grid[i]), -vals[i], law_ticks, True))
        log("a1 search done for horizon", hz)

    # fresh-seed evaluation of every family's best law, plus the Gaussian reference
    P_fresh = {DAY: _n(200000), WEEK: _n(60000)}
    for k, (family, hz, x, tr, maker, normals) in enumerate(specs):
        ok, q, B = maker(x)
        a = fresh_iid(q, B, P_fresh[hz], hz, SEED + 900 + k, normals=True)
        r = rate_row(a, hz)
        rows.append(dict(part="a1", family=family, B=B, params=str(np.round(x, 5).tolist()), train_rate=tr, **r))
        log("fresh", family, hz, r["rate"], r["ci_low"], r["ci_high"])
    for hz in HORIZONS:
        a = fresh_iid(gaussian_q, B_GAUSS, P_fresh[hz], hz, SEED + 950 + hz // DAY, normals=True)
        rows.append(dict(part="a1", family="Gaussian (reference, as in E9)", B=B_GAUSS, params="", train_rate=np.nan,
                         **rate_row(a, hz)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e16a_iid.csv", index=False)
    pd.DataFrame(trace).to_csv(OUT / "e16a_search_trace_iid.csv", index=False)
    print(df[["family", "horizon", "B", "params", "train_rate", "rate", "ci_low", "ci_high", "rate_over_bound"]]
          .round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Part a2: dependent and adaptive laws
# ----------------------------------------------------------------------------


def markov_shape_Y(U, V, Z, s_wild, s_calm, B):
    """Shape switching with constant conditional mean: regime 'wild' draws Y from {0, B}, regime
    'calm' from the Gaussian shape (both mean 1, same B = 9.045). Regime chain with stay
    probabilities s_wild, s_calm, started from its stationary law; independent of the returns."""
    P, n = U.shape
    pw = (1 - s_calm) / ((1 - s_wild) + (1 - s_calm))
    reg = np.empty((P, n), bool)
    reg[:, 0] = V[:, 0] < pw
    for t in range(1, n):
        reg[:, t] = np.where(reg[:, t - 1], V[:, t] < s_wild, V[:, t] >= s_calm)
    wild = two_point_draw(U, np.zeros_like(U), np.full_like(U, B))
    return np.where(reg, wild, gaussian_q(U, Z))


def part_a2():
    rows, trace = [], []
    B = B_GAUSS
    # (i) exogenous shape switching: grid search on CRN, then fresh seeds
    for hz in HORIZONS:
        P = _n(20000 if hz == DAY else 6000)
        rng = np.random.default_rng(SEED + 110 + hz // DAY)
        U, V, Z = rng.random((P, hz)), rng.random((P, hz)), rng.standard_normal((P, hz))
        best = None
        for s_w in (0.5, 0.9, 0.99, 0.999):
            for s_c in (0.5, 0.9, 0.99, 0.999):
                a = engine(array_y_fn(markov_shape_Y(U, V, Z, s_w, s_c, B)), P, hz, B)["alarm"]
                v = float(np.mean(a >= 0))
                trace.append(dict(part="a2", family="shape switching", horizon=hz, params=str([s_w, s_c]), train_rate=v))
                if best is None or v > best[0]:
                    best = (v, s_w, s_c)
        Pf = _n(100000 if hz == DAY else 40000)
        rng = np.random.default_rng(SEED + 910 + hz // DAY)
        al = []
        ch = min(_chunk(hz, 4e6), Pf)
        # full chunks only: 10 x 10,256 = 102,560 paths at one day and 20 x 2,051 = 41,020 at five days
        for s in range(0, Pf, ch):
            U, V, Z = rng.random((ch, hz)), rng.random((ch, hz)), rng.standard_normal((ch, hz))
            al.append(engine(array_y_fn(markov_shape_Y(U, V, Z, best[1], best[2], B)), ch, hz, B)["alarm"])
        r = rate_row(np.concatenate(al), hz)
        rows.append(dict(part="a2", family="exogenous shape switching ({0, B} and Gaussian)", B=B,
                         params=str([best[1], best[2]]), train_rate=best[0], **r))
        log("shape switching", hz, best, r["rate"])

    # (ii) adaptive: bold play, policy search over the 8-parameter family, hold-then-bold
    for hz in HORIZONS:
        P = _n(20000 if hz == DAY else 8000)
        U, _ = crn(P, hz, SEED + 120 + hz // DAY)
        fe = [0]

        def obj(th, U=U, hz=hz):
            a = run_bold(th, U, B)["alarm"]
            v = float(np.mean(a >= 0))
            fe[0] += 1
            trace.append(dict(part="a2", family="adaptive (8-parameter family)", horizon=hz, eval=fe[0],
                              params=str(np.round(th, 4).tolist()), train_rate=v))
            return -v

        v_bold = -obj(BOLD)
        # random search on a quarter of the paths, then Nelder-Mead from the two best and from bold play
        rs = np.random.default_rng(SEED + 130 + hz // DAY)
        cand = []
        Uq = U[: P // 4]
        for i in range(_f(30 if hz == DAY else 16)):
            th = rs.uniform(-1.5, 1.5, 8)
            th[0] = rs.uniform(0.3, 1.5)
            th[4] = rs.uniform(-1.0, 0.5)
            a = run_bold(th, Uq, B)["alarm"]
            cand.append((float(np.mean(a >= 0)), th))
        cand.sort(key=lambda z: -z[0])
        starts = [c[1] for c in cand[: (2 if hz == DAY else 1)]]
        best = _search(obj, starts, _f(120 if hz == DAY else 50), f"adaptive {hz}", trace)
        theta_opt = best.x if -best.fun > v_bold else BOLD
        # hold, then bold
        holds = [0, 50, 100, 150, 200, 250, 300] if hz == DAY else [0, 250, 500, 1000, 1300, 1600]
        hv = []
        for hd in holds:
            a = run_bold(BOLD, U, B, hold=hd)["alarm"]
            hv.append(float(np.mean(a >= 0)))
            trace.append(dict(part="a2", family="hold, then bold", horizon=hz, params=str([hd]), train_rate=hv[-1]))
        hold_best = holds[int(np.argmax(hv))]
        log("adaptive", hz, "bold", v_bold, "opt", -best.fun, "hold", dict(zip(holds, hv)))
        Pf = _n(200000 if hz == DAY else 60000)
        for fam, th, hd, tr in (("adaptive: bold play (spread {0, B}, exact landing)", BOLD, 0, v_bold),
                                ("adaptive: best of the 8-parameter family", theta_opt, 0, max(-best.fun, v_bold)),
                                ("adaptive: hold, then bold", BOLD, hold_best, max(hv))):
            a = fresh_bold(th, B, Pf, hz, SEED + 920 + hz // DAY, hold=hd)
            r = rate_row(a, hz)
            rows.append(dict(part="a2", family=fam, B=B, params=str(np.round(th, 4).tolist()) + (f", hold {hd}" if hd else ""),
                             train_rate=tr, **r))
            log("fresh", fam, hz, r["rate"], r["ci_low"], r["ci_high"])

    # (iii) bold play over B, and the whole horizon curve at the Gaussian B (fresh seeds)
    curve = []
    for i, Bv in enumerate((2.0, 3.0, 5.0, B_GAUSS, 15.0, 25.0, 39.0)):
        a = fresh_bold(BOLD, Bv, _n(60000), WEEK, SEED + 930 + i)
        for hz in HORIZONS:
            curve.append(dict(B=Bv, curve="B scan", **rate_row(a, hz)))
        log("B scan", Bv, [round(float(np.mean((a >= 0) & (a < h))), 4) for h in HORIZONS])
    a = fresh_bold(BOLD, B_GAUSS, _n(30000), 7800, SEED + 940)
    g = fresh_iid(gaussian_q, B_GAUSS, _n(30000), 7800, SEED + 941, normals=True)
    tp = fresh_iid(law_two_point((0.0, B_GAUSS), B_GAUSS)[1], B_GAUSS, _n(30000), 7800, SEED + 942, normals=False)
    for m in (98, 195, 390, 780, 1170, 1950, 2925, 3900, 5850, 7800):
        curve.append(dict(B=B_GAUSS, curve="adaptive: bold play", **rate_row(a, m)))
        curve.append(dict(B=B_GAUSS, curve="Gaussian", **rate_row(g, m)))
        curve.append(dict(B=B_GAUSS, curve="i.i.d. two-point {0, B}", **rate_row(tp, m)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e16a_dependent.csv", index=False)
    pd.DataFrame(curve).to_csv(OUT / "e16a_curves.csv", index=False)
    pd.DataFrame(trace).to_csv(OUT / "e16a_search_trace_dependent.csv", index=False)
    print(df[["family", "horizon", "params", "train_rate", "rate", "ci_low", "ci_high", "rate_over_bound"]]
          .round(4).to_string())
    print(pd.DataFrame(curve).pivot_table(index=["curve", "B"], columns="horizon", values="rate").round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Part a3: where the remaining slack sits, and a check against the exact DP of E11
# ----------------------------------------------------------------------------


def decomposition(chunks, m, up=UPG, dnf=DNF):
    """Optional-stopping identity (docs/theory_delay.md 8.1), exact under an exact null:
    P(sigma <= m) E[N_sigma | sigma <= m] + E[N_m ; sigma > m] = E[sigma ^ m], N = (M_up + M_down) / 2.
    ``chunks``: engine results run with ``decompose=True`` on independent batches of paths."""
    w_up = 0.5 / up.size
    w_dn = 0.5 / dnf.size
    P, n_hit, sum_stop, sum_n_alarm = 0, 0, 0.0, 0.0
    stranded = np.zeros(up.size + dnf.size)
    for res in chunks:
        al = res["alarm"]
        hit = (al >= 0) & (al < m)
        P += al.size
        n_hit += int(hit.sum())
        sum_stop += float(np.where(hit, al + 1, m).sum())
        sum_n_alarm += float(res["n_at_alarm"][hit].sum())
        stranded += np.concatenate([w_up * res["Su_end"], w_dn * res["Sd_end"]], axis=1).sum(axis=0)
    stranded /= P
    E_unc = stranded.sum()
    E_stop = sum_stop / P
    E_N_alarm = sum_n_alarm / n_hit
    rate = n_hit / P
    lo, hi = wilson(n_hit, P)
    out = dict(horizon=m, rate=rate, ci_low=lo, ci_high=hi, bound=m / A, rate_over_bound=rate / (m / A),
               F1_stopping=E_stop / m, F2_crossed=1 - E_unc / E_stop, F3_overshoot=A / E_N_alarm,
               mean_N_at_alarm_over_A=E_N_alarm / A, identity_lhs=rate * E_N_alarm + E_unc, identity_rhs=E_stop, paths=P)
    out["product"] = out["F1_stopping"] * out["F2_crossed"] * out["F3_overshoot"]
    labels = [f"up {g:g}" for g in up] + [f"down {f:g} x max" for f in dnf]
    for lab, v in zip(labels, stranded):
        out[f"stranded share: {lab}"] = v / E_unc
    return out


def part_a3():
    rows = []
    m = DAY
    P, ch = _n(100000), min(_n(100000), 20000)
    for law, seed in (("adaptive: bold play", SEED + 960), ("Gaussian", SEED + 961), ("i.i.d. two-point {0, B}", SEED + 962)):
        rng = np.random.default_rng(seed)
        chunks = []
        for _ in range(P // ch):
            U = rng.random((ch, m))
            if law == "adaptive: bold play":
                chunks.append(run_bold(BOLD, U, B_GAUSS, decompose=True))
                continue
            Y = gaussian_q(U, rng.standard_normal((ch, m))) if law == "Gaussian" else two_point_draw(U, 0.0, B_GAUSS)
            chunks.append(engine(array_y_fn(Y), ch, m, B_GAUSS, decompose=True))
        rows.append(dict(law=law, **decomposition(chunks, m)))
        log("decomposition", law, round(rows[-1]["rate"], 4))
    dec = pd.DataFrame(rows)
    dec.to_csv(OUT / "e16a_slack.csv", index=False)
    print(dec.T.to_string())

    # bold play on single bets and on each side alone against the E11 dynamic programme
    dp = pd.read_csv(ROOT / "results" / "e11d_worst_case_dp.csv")
    chk = []
    lam_max = 1.0 / (B_GAUSS - 1.0)

    def run_side(up, dn, thr, n, P, seed):
        rng = np.random.default_rng(seed)
        ch = _chunk(n)
        out = []
        for _ in range(int(np.ceil(P / ch))):
            U = rng.random((ch, n))
            pol = bold_policy(BOLD, B_GAUSS, n, up=up, thr=thr)
            out.append(engine(adaptive_y_fn(pol, U), ch, n, B_GAUSS, thr=thr, up=up, dnf=dn)["alarm"])
        return np.concatenate(out)

    for i, (direction, lam) in enumerate((("up", 0.5), ("up", 0.2), ("up", 0.05), ("down", 0.95 * lam_max))):
        up = np.array([lam]) if direction == "up" else np.zeros(0)
        dn = np.array([lam / lam_max]) if direction == "down" else np.zeros(0)
        al = run_side(up, dn, A, 3900, _n(40000), SEED + 970 + i)
        for h in (98, 390, 1950, 3900):
            r = rate_row(al, h)
            d = dp[(dp.direction == direction) & np.isclose(dp.lam, lam) & (dp.horizon == h)]
            chk.append(dict(detector=f"single bet {direction} {lam:.4g} at A", **r,
                            dp_worst_case_lower=float(d.worst_case_lower.iloc[0]) if len(d) else np.nan,
                            dp_over_bound=float(d.lower_over_bound.iloc[0]) if len(d) else np.nan))
        log("single bet", direction, lam, [round(c["rate_over_bound"], 3) for c in chk[-4:]])
    # each side of the default grid alone, at 2A (at most half the two-sided budget each)
    for i, (name, up, dn) in enumerate((("increase side alone (default grid) at 2A", UPG, np.zeros(0)),
                                        ("decrease side alone (default grid) at 2A", np.zeros(0), DNF))):
        al = run_side(up, dn, THR, WEEK, _n(40000), SEED + 980 + i)
        for h in HORIZONS:
            r = rate_row(al, h)
            r["bound"] = h / (2 * A)
            r["rate_over_bound"] = r["rate"] / r["bound"]
            r["above_bound"] = bool(r["ci_low"] > r["bound"])
            chk.append(dict(detector=name, **r, dp_worst_case_lower=np.nan, dp_over_bound=np.nan))
        log(name, [round(c["rate"], 4) for c in chk[-2:]])
    chk = pd.DataFrame(chk)
    chk.to_csv(OUT / "e16a_single_bet_check.csv", index=False)
    print(chk.round(4).to_string())


def part_a4():
    """B below the searched range. Under the theorem tau is any level fixed in advance, so B = tau^2 / m0 can be
    any number above 1; when tau is 3 x the scaled MAD of the law itself, B lies in [1.90, 39.6] (report, set-up;
    when the median absolute return is positive; part a5 covers larger B).
    Bold play and the i.i.d. law {0, B} at B = 1.2, 1.5, 1.9 and 2, fresh seeds."""
    rows = []
    for i, Bv in enumerate((1.2, 1.5, 1.9, 2.0)):
        a = fresh_bold(BOLD, Bv, _n(60000), WEEK, SEED + 985 + i)
        t = fresh_iid(law_two_point((0.0, Bv), Bv)[1], Bv, _n(60000), WEEK, SEED + 975 + i, normals=False)
        for hz in HORIZONS:
            rows.append(dict(B=Bv, law="adaptive: bold play", **rate_row(a, hz)))
            rows.append(dict(B=Bv, law="i.i.d. two-point {0, B}", **rate_row(t, hz)))
        log("low B", Bv, [round(r["rate"], 4) for r in rows[-4:]])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e16a_low_B.csv", index=False)
    print(df.round(4).to_string())


def part_a5():
    """B above the searched range. With more than half of the returns at zero the monitor's fallback (median of
    the nonzero |r|) gives B above 39.6: 47.6 to 149 for the price ticks of 2 to 6 sd in part b3. Bold play at
    B = 100, 150 and 1,000, 60,000 fresh paths of five days each."""
    rows = []
    for i, Bv in enumerate((100.0, 150.0, 1000.0)):
        a = fresh_bold(BOLD, Bv, _n(60000), WEEK, SEED + 600 + i)
        for hz in HORIZONS:
            rows.append(dict(B=Bv, law="adaptive: bold play", **rate_row(a, hz)))
        log("high B", Bv, [round(r["rate"], 4) for r in rows[-2:]])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e16a_high_B.csv", index=False)
    print(df.round(4).to_string())


def floor_constant(m, up=UPG, dnf=DNF):
    """Smallest value of N_m = (M_up + M_down) / 2 after m steps over every path of draws in [0, B].

    Each bet's statistic S_t = f_t (S_{t-1} + 1) is nondecreasing in its factor f_t and in S_{t-1}, so it
    is smallest when every draw goes against that bet: Y = 0 for an increase bet (f = 1 - lam) and Y = B
    for a decrease bet (f = 1 - phi (B - 1) = 1 - frac, whatever B). Then S_m = (1 - g)(1 - (1 - g)^m) / g
    with g = lam or frac. The bound holds bet by bet on every path, so also for the averages."""
    su = (1 - up) * (1 - (1 - up) ** m) / up
    sd = (1 - dnf) * (1 - (1 - dnf) ** m) / dnf
    return 0.5 * (su.mean() + sd.mean())


def part_a6():
    """What can be proved, next to what was found (no random draws).

    (i) Floor bound: under an exact null N_t - t is a martingale, so E[N at sigma ^ m] = E[sigma ^ m] <= m. At
    an alarm N >= A; on a path without an alarm by step m, N_m >= c_m (floor_constant). Hence
    A P + c_m (1 - P) <= m, that is P(false alarm within m) <= (m - c_m) / (A - c_m), for every null law.
    (ii) E11's dynamic programme (design.worst_case_false_alarm, the call of experiments/e11_delay_theory.py)
    for the single increase bet 0.05 at level A on finer log grids of wealth states. Its values are the
    false alarm probabilities of an actual adversary restricted to the grid, so they are lower bounds on
    the worst case; a coarse grid rounds the step M -> M + 1 (a draw at Y = 1) back down to M once M is
    large, which is how a small bet builds its wealth."""
    from avmonitor import design

    rows = [dict(kind="floor bound (proved)", horizon=m, bound=m / A, c_m=floor_constant(m),
                 value=(m - floor_constant(m)) / (A - floor_constant(m)))
            for m in (98, 195, 390, 780, 1950, 3900)]
    for ns, H in ((700, 3900), (1400, 1950), (2800, 3900), (5600, 390)):
        t0 = time.time()
        lo, _ = design.worst_case_false_alarm(0.05, B_GAUSS, A, H, "up", n_states=ns, n_left=30, n_right=90,
                                              compute_upper=False)
        for m in (98, 390, 1950, 3900):
            if m <= H:
                rows.append(dict(kind="dynamic programme, single bet up 0.05 at A (lower bound)", n_states=ns,
                                 horizon=m, bound=m / A, value=float(lo[m - 1])))
        log("dp", ns, H, [round(r["value"] / r["bound"], 4) for r in rows if r.get("n_states") == ns],
            f"{time.time() - t0:.0f} s")
    df = pd.DataFrame(rows)
    df["value_over_bound"] = df.value / df.bound
    df.to_csv(OUT / "e16a_floor_and_dp.csv", index=False)
    print(df.round(4).to_string())


# ----------------------------------------------------------------------------
# Part b1: GARCH(1,1)
# ----------------------------------------------------------------------------


def garch(P, n, alpha, beta, rng, burn):
    """GARCH(1,1) with unit unconditional variance; returns r and the predictable variance s2."""
    omega = 1.0 - alpha - beta
    s2 = np.ones(P)
    r = rng.standard_normal(P)
    R = np.empty((P, n))
    S = np.empty((P, n))
    for t in range(burn + n):
        s2 = omega + alpha * r * r + beta * s2
        r = np.sqrt(s2) * rng.standard_normal(P)
        if t >= burn:
            R[:, t - burn] = r
            S[:, t - burn] = s2
    return R, S


def part_b1():
    rows = []
    pers_grid = (0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999)
    alpha_grid = (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3)
    k = 0
    for pers in pers_grid:
        for alpha in alpha_grid:
            k += 1
            beta = pers - alpha
            burn = int(max(1000, 7 * np.log(2) / -np.log(pers)))
            rng = np.random.default_rng(SEED + 200 + k)
            R, S = garch(_n(400), 10000, alpha, beta, rng, burn)
            tau = C * MAD_NORMAL_SCALE * float(np.median(np.abs(R)))
            mu = S * trunc_mean_normal(tau / np.sqrt(S))  # E[X_t | F_{t-1}] exactly (Gaussian innovations)
            m0 = float(mu.mean())
            dev = mu / m0 - 1
            drift = dict(drift_mean_abs=float(np.mean(np.abs(dev))), drift_sd=float(dev.std()),
                         drift_q01=float(np.quantile(dev, 0.01)), drift_q99=float(np.quantile(dev, 0.99)))
            R, S = garch(_n(6000), WEEK, alpha, beta, rng, burn)
            X = np.minimum(R * R, tau**2)
            del R
            a_unc, s_unc = first_alarm(X, m0, m0, tau**2)
            mu = S * trunc_mean_normal(tau / np.sqrt(S))
            del S
            a_orc, _ = first_alarm(X, mu, mu, tau**2)  # the same detector betting against the true E[X_t | F_{t-1}]
            del X, mu
            for hz in HORIZONS:
                ru, ro = rate_row(a_unc, hz), rate_row(a_orc, hz)
                rows.append(dict(persistence=pers, alpha=alpha, beta=beta, half_life=np.log(0.5) / np.log(pers),
                                 tau=tau, m0=m0, **drift, horizon=hz, bound=hz / A,
                                 rate=ru["rate"], ci_low=ru["ci_low"], ci_high=ru["ci_high"],
                                 rate_over_bound=ru["rate_over_bound"], above_bound=ru["above_bound"],
                                 up_share=float(np.mean((a_unc >= 0) & (a_unc < hz) & (s_unc == 1)) / max(ru["rate"], 1e-12)),
                                 oracle_rate=ro["rate"], oracle_ci_low=ro["ci_low"], oracle_ci_high=ro["ci_high"],
                                 oracle_above_bound=ro["above_bound"], paths=ru["paths"]))
            log("garch", pers, alpha, "drift", round(drift["drift_mean_abs"], 4),
                [(r["rate"], r["oracle_rate"]) for r in rows[-2:]])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e16b_garch.csv", index=False)
    print(df.pivot_table(index="persistence", columns="alpha", values="rate_over_bound", aggfunc="first")
          .loc[:, :].round(2).to_string())
    # boundary: largest alpha on the grid with the rate at or below the bound, per persistence and horizon
    return df


# ----------------------------------------------------------------------------
# Part b2: plug-in baselines
# ----------------------------------------------------------------------------


def law_returns(kind, par, rng, shape):
    """Return draws (unit scale) and the function s -> m(s) = E[min(r^2, s^2)] of the law."""
    if kind == "gaussian":
        return rng.standard_normal(shape), lambda s: trunc_mean_normal(s)
    if kind == "scale mixture":
        w, k = par
        sd = np.where(rng.random(shape) < w, k, 1.0)
        return sd * rng.standard_normal(shape), lambda s: (1 - w) * trunc_mean_normal(s) + w * k**2 * trunc_mean_normal(s / k)
    if kind == "student t":
        (nu,) = par
        return rng.standard_t(nu, shape), lambda s: t_trunc_mean(s, nu)
    if kind == "ticks":
        (h,) = par
        j, p, _, _ = tick_law(h)
        r = h * np.round(rng.standard_normal(shape) / h)
        return r, lambda s: np.sum(p[None, :] * np.minimum((j[None, :] * h) ** 2, np.asarray(s)[:, None] ** 2), axis=1)
    raise ValueError(kind)


def law_true_tau(kind, par):
    if kind == "gaussian":
        return TAU_G
    if kind == "scale mixture":
        return scale_mixture_truth(*par)[0]
    if kind == "student t":
        return C * MAD_NORMAL_SCALE * stats.t.ppf(0.75, par[0])
    if kind == "ticks":
        return tick_law(par[0])[2]
    raise ValueError(kind)


def plugin_baseline(r_cal):
    """Exactly VolatilityMonitor's plug-in calibration: scaled MAD about zero, tau = c x scale, mean of X."""
    scale = np.median(np.abs(r_cal), axis=1) * MAD_NORMAL_SCALE
    tau2 = (C * scale) ** 2
    m = np.mean(np.minimum(r_cal**2, tau2[:, None]), axis=1)
    return tau2, m


def plugin_run(kind, par, cal, P, seed, n=WEEK, q_shape=1.0):
    """Three detectors on the same paths: plug-in (tau-hat, m-hat), tau-hat oracle (tau-hat, m(tau-hat),
    inside the class) and the full oracle (true tau and m0). q_shape < 1: after calibration, each
    monitored X is the law's own draw with probability q_shape and the constant m(tau-hat)
    otherwise (same conditional mean, so still inside the null class for the tau-hat oracle)."""
    rng = np.random.default_rng(seed)
    tau = law_true_tau(kind, par)
    out = {k: [] for k in ("plugin", "plugin_side", "tau_hat_oracle", "oracle", "err")}
    ch = _chunk(cal + n, 1.2e7)
    for s0 in range(0, P, ch):
        p = min(ch, P - s0)
        r, mfun = law_returns(kind, par, rng, (p, cal + n))
        m0 = float(np.atleast_1d(mfun(np.atleast_1d(tau)))[0])
        tau2_hat, m_hat = plugin_baseline(r[:, :cal])
        assert np.all(tau2_hat > 0) and np.all(m_hat > 0), "degenerate calibration window"
        m_true_hat = np.asarray(mfun(np.sqrt(tau2_hat)), float)
        X = np.minimum(r[:, cal:] ** 2, tau2_hat[:, None])
        Xo = np.minimum(r[:, cal:] ** 2, tau**2)
        del r
        if q_shape < 1:
            keep = rng.random((p, n)) < q_shape
            X = np.where(keep, X, m_true_hat[:, None])
            Xo = np.where(keep, Xo, m0)
        a, sd = first_alarm(X, m_hat, m_hat, tau2_hat)
        out["plugin"].append(a)
        out["plugin_side"].append(sd)
        out["tau_hat_oracle"].append(first_alarm(X, m_true_hat, m_true_hat, tau2_hat)[0])
        out["oracle"].append(first_alarm(Xo, m0, m0, tau**2)[0])
        out["err"].append(m_hat / m_true_hat - 1)
    return {k: np.concatenate(v) for k, v in out.items()}


def part_b2():
    laws = [("gaussian", ())]
    laws += [("scale mixture", (w, k)) for w in (0.02, 0.05, 0.1, 0.2, 0.35) for k in (1.5, 2.0, 3.0, 5.0, 8.0)]
    laws += [("student t", (nu,)) for nu in (2.2, 2.5, 3.0, 4.0, 6.0, 10.0)]
    laws += [("ticks", (h,)) for h in (0.1, 0.25, 0.5, 0.75, 0.9)]  # 0.9: 35% zeros, so the window median is never 0
    search = []
    for cal in (390, 1440):
        for i, (kind, par) in enumerate(laws):
            res = plugin_run(kind, par, cal, _n(6000), SEED + 300 + i + 100 * (cal == 1440))
            for hz in HORIZONS:
                search.append(dict(cal=cal, law=kind, params=str(list(par)), horizon=hz, bound=hz / A,
                                   plugin=float(np.mean((res["plugin"] >= 0) & (res["plugin"] < hz))),
                                   tau_hat_oracle=float(np.mean((res["tau_hat_oracle"] >= 0) & (res["tau_hat_oracle"] < hz))),
                                   oracle=float(np.mean((res["oracle"] >= 0) & (res["oracle"] < hz))),
                                   err_sd=float(res["err"].std())))
        log("b2 search cal", cal, "done")
    search = pd.DataFrame(search)
    search.to_csv(OUT / "e16b_plugin_search.csv", index=False)
    # fresh seeds for the worst law of each family and the Gaussian, per calibration and horizon
    rows = []
    k = 0
    for cal in (390, 1440):
        for hz in HORIZONS:
            sub = search[(search.cal == cal) & (search.horizon == hz)]
            picks = sub.loc[sub.groupby("law").plugin.idxmax()]
            for _, pk in picks.iterrows():
                par = tuple(ast.literal_eval(pk.params))
                k += 1
                res = plugin_run(pk.law, par, cal, _n(30000), SEED + 990 + k)
                d = dict(cal=cal, law=pk.law, params=pk.params, selected_for_horizon=hz, search_rate=pk.plugin)
                for name in ("plugin", "tau_hat_oracle", "oracle"):
                    r = rate_row(res[name], hz)
                    d.update({f"{name}_rate": r["rate"], f"{name}_ci_low": r["ci_low"], f"{name}_ci_high": r["ci_high"],
                              f"{name}_above_bound": r["above_bound"]})
                d.update(horizon=hz, bound=hz / A, paths=r["paths"], err_sd=float(res["err"].std()))
                rows.append(d)
                log("b2 fresh", cal, hz, pk.law, pk.params, round(d["plugin_rate"], 4), round(d["tau_hat_oracle_rate"], 4))
    # the shape adversary (constant |r| after calibration with probability 1 - q), Gaussian calibration
    for cal in (390, 1440):
        for q in (1.0, 0.5, 0.2, 0.05, 0.0):
            k += 1
            res = plugin_run("gaussian", (), cal, _n(30000), SEED + 990 + k, q_shape=q)
            for hz in HORIZONS:
                d = dict(cal=cal, law="gaussian calibration, shape adversary", params=f"q = {q}",
                         selected_for_horizon=np.nan, search_rate=np.nan)
                for name in ("plugin", "tau_hat_oracle", "oracle"):
                    r = rate_row(res[name], hz)
                    d.update({f"{name}_rate": r["rate"], f"{name}_ci_low": r["ci_low"], f"{name}_ci_high": r["ci_high"],
                              f"{name}_above_bound": r["above_bound"]})
                d.update(horizon=hz, bound=hz / A, paths=r["paths"], err_sd=float(res["err"].std()))
                rows.append(d)
            # rate by size of the baseline error (the assumption that fails)
            if q in (1.0, 0.0):
                e = res["err"]
                bins = np.quantile(e, np.linspace(0, 1, 11))
                b = np.clip(np.searchsorted(bins, e, side="right") - 1, 0, 9)
                for j in range(10):
                    sel = b == j
                    for hz in HORIZONS:
                        rows.append(dict(cal=cal, law="gaussian calibration, rate by baseline error decile",
                                         params=f"q = {q}, decile {j + 1}", selected_for_horizon=np.nan,
                                         search_rate=np.nan, horizon=hz, bound=hz / A, paths=int(sel.sum()),
                                         err_low=float(bins[j]), err_high=float(bins[j + 1]),
                                         plugin_rate=float(np.mean((res["plugin"][sel] >= 0) & (res["plugin"][sel] < hz))),
                                         tau_hat_oracle_rate=float(np.mean((res["tau_hat_oracle"][sel] >= 0)
                                                                           & (res["tau_hat_oracle"][sel] < hz)))))
            log("b2 shape adversary", cal, q, [round(rw.get("plugin_rate", np.nan), 4) for rw in rows[-2:]])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e16b_plugin.csv", index=False)
    print(df[["cal", "law", "params", "horizon", "plugin_rate", "plugin_ci_low", "plugin_ci_high", "tau_hat_oracle_rate",
              "oracle_rate", "err_sd"]].round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Part b3: price-level ticks
# ----------------------------------------------------------------------------


def price_tick_returns(P, n, h, rng):
    """Efficient log price: Gaussian random walk (sd 1 per step), started at a uniform position
    within a tick; observed price h x round(price / h); returns are its first differences."""
    eff = np.cumsum(rng.standard_normal((P, n + 1)), axis=1)
    eff += h * (rng.random(P) - 0.5)[:, None]
    eff /= h
    np.round(eff, out=eff)
    return h * np.diff(eff, axis=1)


def marginal_iid_draws(k, h, shape, rng):
    """I.i.d. returns from the empirical marginal law of the price-rounded returns h k (k: whole ticks, pooled
    over all paths and steps). The dependence between steps is removed, the values and their frequencies kept."""
    vals, cnt = np.unique(np.asarray(k).ravel(), return_counts=True)
    return h * rng.choice(vals, size=shape, p=cnt / cnt.sum()).astype(float)


def part_b3():
    rows = []
    for i, h in enumerate((0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0)):
        rng = np.random.default_rng(SEED + 400 + i)
        long = price_tick_returns(_n(1000), 10000, h, rng)
        k = np.round(long / h).astype(int)
        absval = np.abs(k).ravel()
        vals, cnt = np.unique(absval, return_counts=True)
        tau = tick_tau(vals * h, cnt / cnt.sum())
        X = np.minimum(long**2, tau**2)
        m0 = float(X.mean())
        prev = np.abs(k[:, :-1])
        nxt = X[:, 1:]
        cond = {f"cond_mean_after_{lab}": float(nxt[sel].mean() / m0)
                for lab, sel in (("zero", prev == 0), ("one_tick", prev == 1), ("two_or_more", prev >= 2))
                if sel.any()}
        ac = float(np.corrcoef(long[:, :-1].ravel(), long[:, 1:].ravel())[0, 1])
        acx = float(np.corrcoef(X[:, :-1].ravel(), X[:, 1:].ravel())[0, 1])
        r = price_tick_returns(_n(10000), WEEK, h, rng)
        a_price, _ = first_alarm(np.minimum(r * r, tau**2), m0, m0, tau**2)
        del r
        # twin, inside the class but a different law: each return rounded independently (its own exact tau and m0)
        j, p, tau_i, m0_i = tick_law(h, kmax=int(np.ceil(12 / h)) + 1)
        z = h * np.round(rng.standard_normal((_n(10000), WEEK)) / h)
        a_iid, _ = first_alarm(np.minimum(z * z, tau_i**2), m0_i, m0_i, tau_i**2)
        del z
        # control inside the class with the same marginal law: i.i.d. draws from the pooled price-rounded returns of
        # the long paths, with the same tau and m0 (m0 is exactly the mean of min(r^2, tau^2) under this law)
        ctrl = marginal_iid_draws(k, h, (_n(10000), WEEK), np.random.default_rng(SEED + 500 + i))
        xc = np.minimum(ctrl * ctrl, tau**2)
        ctrl_zero, ctrl_mean = float(np.mean(ctrl == 0)), float(xc.mean() / m0)
        del ctrl
        a_ctrl, _ = first_alarm(xc, m0, m0, tau**2)
        del xc
        for hz in HORIZONS:
            rp, ri, rc = rate_row(a_price, hz), rate_row(a_iid, hz), rate_row(a_ctrl, hz)
            rows.append(dict(tick=h, tau=tau, m0=m0, B=tau**2 / m0, zero_share=float(np.mean(k == 0)),
                             lag1_corr_returns=ac, lag1_corr_X=acx, **cond, horizon=hz, bound=hz / A,
                             price_rounding_rate=rp["rate"], price_ci_low=rp["ci_low"], price_ci_high=rp["ci_high"],
                             price_rate_over_bound=rp["rate_over_bound"], price_above_bound=rp["above_bound"],
                             marginal_iid_rate=rc["rate"], marginal_iid_ci_low=rc["ci_low"],
                             marginal_iid_ci_high=rc["ci_high"], marginal_iid_rate_over_bound=rc["rate_over_bound"],
                             marginal_iid_above_bound=rc["above_bound"], marginal_iid_zero_share=ctrl_zero,
                             marginal_iid_mean_X_over_m0=ctrl_mean,
                             iid_rounding_rate=ri["rate"], iid_ci_low=ri["ci_low"], iid_ci_high=ri["ci_high"],
                             iid_above_bound=ri["above_bound"], twin_tau=tau_i, twin_m0=m0_i, twin_B=tau_i**2 / m0_i,
                             twin_zero_share=float(p[j == 0].sum()), paths=rp["paths"]))
        log("price ticks", h, cond, [(rw["price_rounding_rate"], rw["marginal_iid_rate"], rw["iid_rounding_rate"])
                                     for rw in rows[-2:]])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e16b_price_ticks.csv", index=False)
    print(df.round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

INK, MUTED, GRIDC = "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical slots 1 to 4 (as in E11)


def _style(ax):
    ax.grid(True, color=GRIDC, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def figures():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    from matplotlib.ticker import FuncFormatter, NullFormatter

    plain = FuncFormatter(lambda v, _: f"{v:g}")

    FIG.mkdir(exist_ok=True)
    # Figure 1: inside the null class
    iid = pd.read_csv(OUT / "e16a_iid.csv")
    dep = pd.read_csv(OUT / "e16a_dependent.csv")
    cur = pd.read_csv(OUT / "e16a_curves.csv")
    allr = pd.concat([iid, dep], ignore_index=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0), gridspec_kw=dict(width_ratios=[1.3, 1]))
    ax = axes[0]
    fams = list(dict.fromkeys(allr.family))
    yy = np.arange(len(fams))
    for k, (hz, col, mk) in enumerate(((DAY, SERIES[0], "o"), (WEEK, SERIES[1], "s"))):
        sub = allr[allr.horizon == hz].set_index("family").reindex(fams)
        b = hz / A
        ax.errorbar(sub.rate / b, yy + (k - 0.5) * 0.28, xerr=[(sub.rate - sub.ci_low) / b, (sub.ci_high - sub.rate) / b],
                    fmt=mk, ms=5, color=col, lw=1, capsize=2, label=f"within {hz} steps (bound {b:.0%})")
    ax.axvline(1.0, color=INK, lw=1)
    ax.set_yticks(yy)
    ax.set_yticklabels(fams, fontsize=7.5, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("false alarm rate / bound m / A (fresh seeds, 95% Wilson intervals)", color=INK, fontsize=9)
    ax.set_title("Worst law found in each family (true baseline; the line is the bound)", fontsize=9.5, color=INK, loc="left")
    _style(ax)
    ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.12), ncol=2)
    ax = axes[1]
    for (name, col) in (("adaptive: bold play", SERIES[0]), ("i.i.d. two-point {0, B}", SERIES[2]), ("Gaussian", SERIES[1])):
        s_ = cur[(cur.curve == name) & (cur.rate > 0)].sort_values("horizon")
        ax.plot(s_.horizon, s_.rate, color=col, lw=2, marker="o", ms=4, label=name)
    hs = np.array(sorted(cur[cur.curve == "Gaussian"].horizon.unique()))
    ax.plot(hs, hs / A, color=INK, lw=1, ls="--", label="bound m / A")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("horizon m (steps; 390 = one day, 7800 = A)", color=INK, fontsize=9)
    ax.set_ylabel("P(false alarm within m)", color=INK, fontsize=9)
    ax.set_title("Horizon curves, B = 9.045 (zeros not shown)", fontsize=9.5, color=INK, loc="left")
    _style(ax)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.suptitle("E16a: adversarial laws inside the null class, against the bound (synthetic, A = 7800, both sides at 2A)",
                 fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "e16_inside_class.png", dpi=140)
    plt.close(fig)

    # Figure 2: outside the class
    g = pd.read_csv(OUT / "e16b_garch.csv")
    pl = pd.read_csv(OUT / "e16b_plugin.csv")
    tk = pd.read_csv(OUT / "e16b_price_ticks.csv")
    fig = plt.figure(figsize=(13, 8.6))
    gs = fig.add_gridspec(2, 6)
    cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])  # diverging: blue, gray, red
    for col0, hz in ((0, DAY), (3, WEEK)):
        ax = fig.add_subplot(gs[0, col0 : col0 + 3])
        sub = g[g.horizon == hz].pivot_table(index="persistence", columns="alpha", values="rate_over_bound")
        vals = np.clip(sub.values, 0.05, 20)
        im = ax.imshow(vals, cmap=cmap, norm=LogNorm(vmin=0.05, vmax=20), aspect="auto", origin="lower")
        for i in range(vals.shape[0]):
            for j in range(vals.shape[1]):
                ax.text(j, i, f"{sub.values[i, j]:.2g}", ha="center", va="center", fontsize=6.5, color=INK)
        ax.set_xticks(range(sub.shape[1]))
        ax.set_xticklabels([f"{a:g}" for a in sub.columns], fontsize=7.5, color=MUTED)
        ax.set_yticks(range(sub.shape[0]))
        ax.set_yticklabels([f"{p:g}" for p in sub.index], fontsize=7.5, color=MUTED)
        ax.set_xlabel("alpha (weight of the last squared return)", color=INK, fontsize=9)
        ax.set_ylabel("persistence alpha + beta", color=INK, fontsize=9)
        extra = "; 4 = every path" if hz == WEEK else ""
        ax.set_title(f"GARCH(1,1): rate / bound within {hz} steps (red: above 1{extra})", fontsize=9.5, color=INK, loc="left")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02).ax.tick_params(labelsize=7, colors=MUTED)
    floor = 0.02
    ax = fig.add_subplot(gs[1, 0:2])
    sub = g[g.horizon == DAY]
    ax.scatter(sub.drift_mean_abs, np.maximum(sub.rate_over_bound, floor), s=22, color=SERIES[1], edgecolor="white",
               lw=0.6, label="betting against m0 (outside the class)", zorder=3)
    ax.scatter(sub.drift_mean_abs, np.maximum(sub.oracle_rate / sub.bound, floor), s=22, color=SERIES[0],
               edgecolor="white", lw=0.6, label="betting against the true μ_t (inside)", zorder=3)
    ax.axhline(1.0, color=INK, lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("conditional mean drift E|μ_t / m0 - 1|", color=INK, fontsize=9)
    ax.set_ylabel(f"one-day rate / bound (0 shown at {floor})", color=INK, fontsize=9)
    ax.set_title("GARCH: failure follows the drift", fontsize=9.5, color=INK, loc="left")
    _style(ax)
    ax.yaxis.set_major_formatter(plain)
    ax.xaxis.set_major_formatter(plain)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax = fig.add_subplot(gs[1, 2:4])
    sh = pl[pl.law == "gaussian calibration, shape adversary"].copy()
    sh["q"] = sh.params.str.replace("q = ", "").astype(float)
    for (cal, hz), col, ls in (((390, DAY), SERIES[0], "-"), ((390, WEEK), SERIES[0], "--"),
                               ((1440, DAY), SERIES[1], "-"), ((1440, WEEK), SERIES[1], "--")):
        s_ = sh[(sh.cal == cal) & (sh.horizon == hz)].sort_values("q")
        ax.plot(s_.q, s_.plugin_rate / (hz / A), color=col, ls=ls, lw=2, marker="o", ms=4,
                label=f"plug-in, cal. {cal}, {hz} steps")
    s_ = sh[(sh.cal == 390) & (sh.horizon == WEEK)].sort_values("q")
    ax.plot(s_.q, np.maximum(s_.tau_hat_oracle_rate / (WEEK / A), floor), color=SERIES[2], lw=2, marker="s", ms=4,
            label="true mean at τ-hat, cal. 390, 1950 steps (0 shown at 0.02)")
    ax.axhline(1.0, color=INK, lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("q: share of monitored draws from the calibration law", color=INK, fontsize=9)
    ax.set_ylabel("rate / bound", color=INK, fontsize=9)
    ax.set_title("Plug-in baseline, i.i.d. (q = 1) to constant |r| (q = 0)", fontsize=9.5, color=INK, loc="left")
    _style(ax)
    ax.yaxis.set_major_formatter(plain)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax = fig.add_subplot(gs[1, 4:6])
    for hz, ls in ((DAY, "-"), (WEEK, "--")):
        s_ = tk[tk.horizon == hz].sort_values("tick")
        ax.plot(s_.tick, np.maximum(s_.price_rounding_rate / (hz / A), floor), color=SERIES[2], ls=ls, lw=2, marker="s",
                ms=4, label=f"price rounded to the tick, {hz} steps")
        ax.plot(s_.tick, np.maximum(s_.marginal_iid_rate / (hz / A), floor), color=SERIES[3], ls=ls, lw=2, marker="o",
                ms=4, label=f"same values drawn i.i.d. (inside), {hz} steps")
    ax.axhline(1.0, color=INK, lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("tick size (standard deviations of one return)", color=INK, fontsize=9)
    ax.set_ylabel("rate / bound", color=INK, fontsize=9)
    ax.set_title("Price-level ticks (true unconditional baseline)", fontsize=9.5, color=INK, loc="left")
    _style(ax)
    ax.set_yticks([0.2, 0.5, 1, 2, 3])
    ax.yaxis.set_major_formatter(plain)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(plain)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle("E16b: outside the null class, each failure matches a broken assumption (synthetic, A = 7800)",
                 fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "e16_outside_class.png", dpi=140)
    plt.close(fig)
    log("figures written")


if __name__ == "__main__":
    parts = sys.argv[1:] or ["a", "b"]
    expand = {"a": ["a1", "a2", "a3", "a4", "a5", "a6"], "b": ["b1", "b2", "b3"]}
    todo = [p for x in parts for p in expand.get(x, [x])]
    fns = dict(a1=part_a1, a2=part_a2, a3=part_a3, a4=part_a4, a5=part_a5, a6=part_a6, b1=part_b1, b2=part_b2,
               b3=part_b3, figures=figures)
    for p in todo:
        t0 = time.time()
        log("start", p)
        fns[p]()
        log("done", p, f"({time.time() - t0:.0f} s)")
