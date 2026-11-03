"""TH6 (referee M2 and M4). What Theorem B bounds after the first day, and two constructions with a proved
guarantee for every period.

All results are SYNTHETIC. Seeds are fixed (listed below). Theory: docs/theory_per_period.md (Theorems 1 and 2,
Proposition 3); docs/math.md, Theorem B and Corollary C (section 5.4) and Theorem 6.4 (Theorem D in the referee report).
Code: src/avmonitor/periodic.py. The detector is the monitor's default with the TRUE baseline: X_t = min(r_t^2,
tau^2), Y_t = X_t / m0 in [0, B], increase bets 0.02 to 0.5, decrease bets 0.1 to 0.95 times 1 / (B - 1), each
side the average of its bets, each side alarming at 2A. A = 7800 for 390-step days and A = 28,800 for 1440-step
days, so m / A = 5% for one day.

Research questions
  Q6a. Under i.i.d. Gaussian, Student t3, jump and tick returns (the four laws of E9, true tau and m0 computed
       exactly), for days d = 1 to 20: P(first false alarm on day d), the conditional hazard P(alarm on day d |
       none before), and, with a restart after every alarm (recalibration with the true baseline), the expected
       number of false alarms on day d and per day over days 1 to d (Corollary C bounds the latter by m / A).
  Q6b. The same for null streams that are not i.i.d.: the referee's stream (|r_t| = sqrt(m0), so Y_t = 1, for
       the first days, then Gaussian), as one stream (calm for 9 days) and as a family (calm for d - 1 days, then
       one Gaussian or two-point day d); and E16's bold play run over 20 days.
  Q6c. Two constructions with a per-period guarantee, on the same streams: Shiryaev-Roberts restarted at every
       day boundary (Theorem 1: at most m / A in every day) and the CUSUM-type e-detector (Theorem 2: at most
       5.146% (390-step days) or 5.044% (1440-step days) in any window of one day under an i.i.d. null). The
       counterexample of Proposition 3: a martingale-difference null that parks the CUSUM statistic and strikes
       on a chosen day, at full scale and scaled down (A = 400, 20-step periods) to show the limit.
  Q6d. Delay cost. A change of the return standard deviation by a factor 1.25, 1.5, 2, 0.8, 0.67 or 0.5 at a
       uniformly random step of day 3 (Gaussian returns). Shiryaev-Roberts at A (first-day guarantee only),
       restarted Shiryaev-Roberts at A (5% per day), restarted with a carry of m at 2A (5% per day), and the
       e-CUSUM at the A that makes its window bound 5% (8028 and 29,053). Delay of the first alarm at or after
       the change, by the position of the change within the day.
  Q6e. (M4) What the tolerance hypothesis of Theorem D excludes: under GARCH(1,1) (the four cells of E9) and a
       slow log-normal stochastic volatility, the smallest b for which sup_{monitoring} mu_t <= (1 + b)
       min_{calibration} mu_i and inf_{monitoring} mu_t >= (1 - b) max_{calibration} mu_i hold over a
       calibration window and one monitored day, with mu_t = E[min(r_t^2, tau^2) | past] exact.

Parts (each writes its own CSV files; run at most 2 processes on a shared machine, with nice -n 10):
  laws        exact tau, m0 and B of the four i.i.d. laws (th6_laws.csv)
  hazard390   Q6a and the calm-9 stream for 390-step days (th6_hazard_390.csv)
  hazard1440  the same for 1440-step days (th6_hazard_1440.csv)
  family      the referee's family, Q6b (th6_calm_family.csv)
  adversary   bold play over 20 days and the park-and-strike counterexample (th6_adversary.csv,
              th6_park_strike.csv, th6_park_scaled.csv)
  delay       Q6d (th6_delay.csv, th6_delay_by_position.csv)
  renewal     Q6d, Proposition 2: the renewal bound on the restarted detector's delay (th6_delay_renewal.csv;
              needs th6_delay.csv; added after the delay results were seen)
  tolerance   Q6e (th6_tolerance.csv)
  figures     results/figures/th6_hazard.png, th6_delay.png
  python experiments/th6_per_period.py laws hazard390 hazard1440 family adversary delay renewal tolerance figures
  TH6_SCALE=0.05 TH6_OUT=/tmp/th6 python experiments/th6_per_period.py ...   # smoke run

Seeds (SEED = 26060): hazard SEED + 1..10 (stream index, plus 5 for 1440-step days); family SEED + 21..100;
adversary SEED + 100..152 (SEED + 100 is also used by one family cell, for an unrelated stream); delay
SEED + 200..215; renewal SEED + 400..415; tolerance SEED + 300..309.
"""

from __future__ import annotations

import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate, optimize, stats

from avmonitor import periodic as Pm
from avmonitor.robust import MAD_NORMAL_SCALE

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("TH6_OUT", ROOT / "results"))
FIG = OUT / "figures"
SCALE = float(os.environ.get("TH6_SCALE", "1"))
SEED = 26060
DAYS = 20
ALPHA = 0.05
SETTINGS = {390: 7800.0, 1440: 28800.0}
LAWS = ("gaussian", "student_t3", "jumps", "ticks")


def _n(P):
    return int(max(200, round(P * SCALE))) if SCALE != 1 else int(P)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ----------------------------------------------------------------------------
# The four i.i.d. laws of E9, with exact tau and m0
# ----------------------------------------------------------------------------


def trunc_mean_normal(k):
    k = np.asarray(k, float)
    return 2 * stats.norm.cdf(k) - 1 - 2 * k * stats.norm.pdf(k) + 2 * k**2 * stats.norm.sf(k)


def law_truth(kind):
    """tau = 3 x scaled MAD of the law (floored at 3 ticks for 'ticks', as in E9) and m0 = E[min(r^2, tau^2)]."""
    if kind == "gaussian":
        tau = 3 * MAD_NORMAL_SCALE * stats.norm.ppf(0.75)
        return tau, float(trunc_mean_normal(tau)), 0.0
    if kind == "student_t3":
        s = math.sqrt(3.0)
        med = stats.t.ppf(0.75, 3) / s
        tau = 3 * MAD_NORMAL_SCALE * med
        dens = lambda x: s * stats.t.pdf(s * x, 3)  # noqa: E731
        inner = 2 * integrate.quad(lambda x: x * x * dens(x), 0, tau, epsabs=1e-13, epsrel=1e-12)[0]
        return tau, float(inner + tau**2 * 2 * stats.t.sf(s * tau, 3)), 0.0
    if kind == "jumps":
        sj = math.sqrt(1 + 100.0)
        F = lambda x: 0.995 * (2 * stats.norm.cdf(x) - 1) + 0.005 * (2 * stats.norm.cdf(x / sj) - 1)  # noqa: E731
        med = optimize.brentq(lambda x: F(x) - 0.5, 0.1, 2.0, xtol=1e-14)
        tau = 3 * MAD_NORMAL_SCALE * med
        m0 = 0.995 * trunc_mean_normal(tau) + 0.005 * sj**2 * trunc_mean_normal(tau / sj)
        return tau, float(m0), 0.0
    if kind == "ticks":
        h = 0.5
        j = np.arange(-40, 41)
        p = stats.norm.cdf(h * j + h / 2) - stats.norm.cdf(h * j - h / 2)
        a = np.abs(j) * h
        v = np.unique(a)
        cdf = np.array([p[a <= x].sum() for x in v])
        med = v[np.searchsorted(cdf, 0.5)]
        tau = max(3 * MAD_NORMAL_SCALE * med, 3 * h)
        return tau, float(np.sum(p * np.minimum(a * a, tau * tau))), float(p[j == 0][0])
    raise ValueError(kind)


TRUTH = {k: law_truth(k) for k in LAWS}


def draw_returns(kind, size, rng):
    if kind == "gaussian":
        return rng.standard_normal(size)
    if kind == "student_t3":
        return rng.standard_t(3, size) / math.sqrt(3.0)
    if kind == "jumps":
        r = rng.standard_normal(size)
        return r + (rng.random(size) < 0.005) * rng.normal(0, 10, size)
    if kind == "ticks":
        return np.round(rng.standard_normal(size) / 0.5) * 0.5
    raise ValueError(kind)


def law_Y(kind, size, rng):
    tau, m0, _ = TRUTH[kind]
    r = draw_returns(kind, size, rng)
    return np.minimum(r * r, tau * tau) / m0


def law_B(kind):
    tau, m0, _ = TRUTH[kind]
    return tau * tau / m0


class DaySampler:
    """y_fn for run_engine: i.i.d. draws of one law, generated one day at a time; Y = 1 before ``calm`` steps."""

    def __init__(self, kind, P, day, rng, calm=0):
        self.kind, self.P, self.day, self.rng, self.calm = kind, P, day, rng, calm
        self.block, self.start = None, None

    def __call__(self, t, eng):
        if t < self.calm:
            return np.ones(self.P)
        if self.block is None or t >= self.start + self.day:
            self.start = t
            self.block = np.ascontiguousarray(law_Y(self.kind, (self.day, self.P), self.rng))
        return self.block[t - self.start]


def part_laws():
    rows = []
    for k in LAWS:
        tau, m0, zero = TRUTH[k]
        rows.append(dict(law=k, tau=tau, m0=m0, B=tau * tau / m0, zero_share=zero))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th6_laws.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Q6a, Q6b, Q6c: per-day tables
# ----------------------------------------------------------------------------


def specs_for(day, A):
    return [Pm.Spec("SR", "SR", A),
            Pm.Spec("SR_restarted_daily", "SR", A, period=day),
            Pm.Spec("CUSUM", "CUSUM", A),
            Pm.Spec("CUSUM_no_restart", "CUSUM", A, on_alarm="continue")]


def table_rows(res, stream, day, A, P, extra=None):
    """Rows of per-day summaries for every detector in ``res``."""
    rows = []
    for name, r in res.items():
        tab = Pm.day_table(r["first"], r["counts"])
        cont = name.endswith("no_restart")  # counts are steps at or above 2A, not alarms: only p_any is used
        cum = np.cumsum(r["counts"], axis=1) / np.arange(1, r["counts"].shape[1] + 1)
        for d in range(r["counts"].shape[1]):
            pf, h, n_r = tab["p_first"][d], tab["hazard"][d], tab["at_risk"][d]
            an, mc = tab["any"][d], tab["mean_count"][d]
            rows.append(dict(
                stream=stream, day_len=day, A=A, detector=name, day=d + 1, paths=P,
                p_first=pf, p_first_se=math.sqrt(pf * (1 - pf) / P),
                at_risk=int(n_r), hazard=h, hazard_se=math.sqrt(h * (1 - h) / n_r) if n_r else np.nan,
                p_any=an, p_any_se=math.sqrt(an * (1 - an) / P),
                mean_count=np.nan if cont else mc, mean_count_se=np.nan if cont else tab["count_sd"][d] / math.sqrt(P),
                cum_rate=np.nan if cont else tab["cum_rate"][d],
                cum_rate_se=np.nan if cont else cum[:, d].std(ddof=1) / math.sqrt(P),
                bound_first_day=day / A, bound_ville_day_d=min(1.0, (d + 1) * day / A),
                bound_corollary_c=day / A, bound_restarted=Pm.restarted_sr_bound(day, A),
                bound_cusum_window=Pm.window_bound_two_sided(day, A), **(extra or {})))
    return rows


def part_hazard(day):
    A = SETTINGS[day]
    P = _n(40000)
    rows = []
    for i, kind in enumerate(LAWS + ("calm9_gaussian",)):
        law = "gaussian" if kind == "calm9_gaussian" else kind
        calm = 9 * day if kind == "calm9_gaussian" else 0
        rng = np.random.default_rng(SEED + 1 + i + (5 if day == 1440 else 0))
        t0 = time.time()
        res = Pm.run_engine(DaySampler(law, P, day, rng, calm=calm), P, DAYS * day, law_B(law), specs_for(day, A), day)
        rows += table_rows(res, kind, day, A, P)
        log(f"hazard {day} {kind}: {time.time() - t0:.0f} s")
        pd.DataFrame(rows).to_csv(OUT / f"th6_hazard_{day}.csv", index=False)
    df = pd.DataFrame(rows)
    show = df[df.day.isin([1, 2, 5, 10, 15, 20])]
    print(show.pivot_table(index=["stream", "day"], columns="detector", values="hazard").round(4).to_string())


def part_family():
    """The referee's family: calm (Y = 1) for d - 1 days, then one day of Gaussian or two-point {0, B} draws.

    After the calm days every Shiryaev-Roberts bet sits at (d - 1) m and every CUSUM bet at 1, deterministically,
    so only day d is simulated from that state. The restarted detector is back at 0 on day d.
    """
    rows = []
    P = _n(40000)
    for day, A in SETTINGS.items():
        B = law_B("gaussian")
        for li, law in enumerate(("gaussian", "two_point")):
            for d in range(1, DAYS + 1):
                rng = np.random.default_rng(SEED + 20 + 40 * li + d + (20 if day == 1440 else 0))
                if law == "gaussian":
                    y_fn = DaySampler("gaussian", P, day, rng)
                else:
                    y_fn = lambda t, eng, rng=rng: np.where(rng.random(P) < 1 / B, B, 0.0)  # noqa: E731
                specs = [Pm.Spec("SR", "SR", A), Pm.Spec("SR_restarted_daily", "SR", A, period=day),
                         Pm.Spec("CUSUM", "CUSUM", A)]
                s0 = float((d - 1) * day)
                init = {"SR": (s0, s0), "SR_restarted_daily": (0.0, 0.0), "CUSUM": (1.0 if d > 1 else 0.0,) * 2}
                res = Pm.run_engine(y_fn, P, day, B, specs, day, t0=(d - 1) * day, init=init)
                for name, r in res.items():
                    p = float(np.mean(r["first"] >= 0))
                    rows.append(dict(day_len=day, A=A, law_on_day_d=law, day=d, detector=name, paths=P,
                                     p_alarm_day_d=p, se=math.sqrt(p * (1 - p) / P), bound_ville=min(1, d * day / A),
                                     bound_restarted=day / A, bound_cusum_window=Pm.window_bound_two_sided(day, A)))
            log(f"family {day} {law} done")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th6_calm_family.csv", index=False)
    print(df[df.day.isin([1, 2, 5, 10, 15, 19, 20])].pivot_table(index=["day_len", "law_on_day_d", "day"],
                                                               columns="detector", values="p_alarm_day_d").round(4))


def park_strike(day, A, P, rng, n_days=DAYS, B=None):
    """Proposition 3 against the default two-sided e-CUSUM: hunt with {0, B}, park with Y = 1 once one draw could
    reach 2A, strike with bold play on day d. Returns rows for target days d = 1..n_days.

    One hunting run serves every target day: the state at the start of day d is the same for every target
    d' >= d, and the strike day is then simulated from that state with fresh draws.
    """
    B = law_B("gaussian") if B is None else B
    spec = Pm.Spec("CUSUM", "CUSUM", A)
    eng = Pm.Engine([spec], P, B)
    snaps, parked_at = [], np.full(P, -1)
    for t in range(n_days * day):
        if t % day == 0:
            snaps.append((eng.Su["CUSUM"].copy(), eng.Sd["CUSUM"].copy()))
        yl, yh, uo, do = eng.landing("CUSUM")
        park = uo | do
        parked_at[park & (parked_at < 0)] = t
        y = np.where(park, 1.0, Pm.two_point(rng.random(P), np.zeros(P), np.full(P, B)))
        if eng.step(y)["CUSUM"].any():
            raise RuntimeError("an alarm during the hunt: the parking rule is wrong")
    rows = []
    for d in range(1, n_days + 1):
        su, sd = snaps[d - 1]
        res = Pm.run_engine(Pm.bold_policy("CUSUM", rng), P, day, B, [spec], day, t0=(d - 1) * day,
                            init={"CUSUM": (su, sd)})
        p = float(np.mean(res["CUSUM"]["first"] >= 0))
        cond = Pm.conditional_window_bound(su, sd, day, A).mean()
        rows.append(dict(day_len=day, A=A, B=B, day=d, paths=P, parked_before_day=float(np.mean((parked_at >= 0) & (parked_at < (d - 1) * day))),
                         p_alarm_day_d=p, se=math.sqrt(p * (1 - p) / P), mean_conditional_bound=cond,
                         bound_cusum_window_iid=Pm.window_bound_two_sided(day, A), one_over_B=1 / B))
    return rows


def part_adversary():
    rows = []
    P = _n(30000)
    for day, A in SETTINGS.items():
        B = law_B("gaussian")
        off = 0 if day == 390 else 50
        for j, (name, spec) in enumerate((("SR", Pm.Spec("SR", "SR", A)),
                                          ("SR_restarted_daily", Pm.Spec("SR_restarted_daily", "SR", A, period=day)),
                                          ("CUSUM", Pm.Spec("CUSUM", "CUSUM", A)))):
            t0 = time.time()
            rng = np.random.default_rng(SEED + 100 + off + j)
            res = Pm.run_engine(Pm.bold_policy(name, rng), P, DAYS * day, B, [spec], day)
            rows += table_rows(res, "bold_play", day, A, P)
            log(f"bold {day} {name}: {time.time() - t0:.0f} s")
            pd.DataFrame(rows).to_csv(OUT / "th6_adversary.csv", index=False)
    df = pd.DataFrame(rows)
    print(df[df.day.isin([1, 2, 5, 10, 20])].pivot_table(index=["day_len", "day"], columns="detector",
                                                         values=["hazard", "mean_count", "cum_rate"]).round(4))
    prow = []
    for day, A in SETTINGS.items():
        t0 = time.time()
        rng = np.random.default_rng(SEED + 120 + (1 if day == 1440 else 0))
        prow += park_strike(day, A, _n(30000), rng)
        log(f"park and strike {day}: {time.time() - t0:.0f} s")
        pd.DataFrame(prow).to_csv(OUT / "th6_park_strike.csv", index=False)
    print(pd.DataFrame(prow)[["day_len", "day", "parked_before_day", "p_alarm_day_d", "se", "bound_cusum_window_iid"]]
          .round(4).to_string())
    # Scaled down, to show the limit as d grows: periods of 20 steps, A = 400 (m / A = 5%), 120 periods
    rng = np.random.default_rng(SEED + 130)
    sc = park_strike(20, 400.0, _n(20000), rng, n_days=120)
    pd.DataFrame(sc).to_csv(OUT / "th6_park_scaled.csv", index=False)
    print(pd.DataFrame(sc)[["day", "parked_before_day", "p_alarm_day_d", "bound_cusum_window_iid"]].iloc[::10].round(4))


# ----------------------------------------------------------------------------
# Q6d: delays
# ----------------------------------------------------------------------------

RATIOS = (1.25, 1.5, 2.0, 0.8, 0.67, 0.5)


def part_delay():
    rows, prow = [], []
    for day, A in SETTINGS.items():
        A_cusum = Pm.threshold_for_window(day, ALPHA)
        specs = [Pm.Spec("SR", "SR", A), Pm.Spec("SR_restarted_daily", "SR", A, period=day),
                 Pm.Spec("SR_restarted_carry_m", "SR", 2 * A, period=day, carry=float(day)),
                 Pm.Spec("CUSUM_certified", "CUSUM", A_cusum)]
        proved = {"SR": "first day after (re)calibration only: m / A = 5%",
                  "SR_restarted_daily": "every day: m / A = 5% (Theorem 1)",
                  "SR_restarted_carry_m": "every day: (m + m) / (2A) = 5% (Theorem 1)",
                  "CUSUM_certified": "any window of m steps, i.i.d. null: 5% (Theorem 2)"}
        n_days = 12
        n = n_days * day
        tau, m0, _ = TRUTH["gaussian"]
        B = law_B("gaussian")
        for i, rho in enumerate(RATIOS):
            P = _n(20000)
            rng = np.random.default_rng(SEED + 200 + i + (10 if day == 1440 else 0))
            nu = 2 * day + rng.integers(0, day, P)  # change at a uniform step of day 3 (0-based index)
            sd = lambda t: np.where(t >= nu, rho, 1.0)  # noqa: E731

            def y_fn(t, eng, rng=rng, sd=sd):
                z = rng.standard_normal(P) * sd(t)
                return np.minimum(z * z, tau * tau) / m0

            t0 = time.time()
            res = Pm.run_engine(y_fn, P, n, B, specs, day, after=nu)
            pos = (nu - 2 * day) / day
            for s in specs:
                r = res[s.name]
                fa = (r["first"] >= 0) & (r["first"] < nu)
                det = r["first_after"] >= 0
                delay = np.where(det, r["first_after"] - nu + 1, n - nu)  # censored at the horizon
                rows.append(dict(day_len=day, ratio=rho, detector=s.name, A=s.A, proved=proved[s.name], paths=P,
                                 false_alarm_before_change=fa.mean(), detected_within_horizon=det.mean(),
                                 mean_delay=delay.mean(), mean_delay_se=delay.std(ddof=1) / math.sqrt(P),
                                 median_delay=float(np.median(delay)), q90_delay=float(np.quantile(delay, 0.9)),
                                 horizon_days_after_change_min=(n - 3 * day) / day))
                for lo, hi, lab in ((0.0, 0.5, "first half"), (0.5, 0.9, "50% to 90%"), (0.9, 1.0, "last 10%")):
                    sel = (pos >= lo) & (pos < hi)
                    prow.append(dict(day_len=day, ratio=rho, detector=s.name, position=lab, paths=int(sel.sum()),
                                     mean_delay=delay[sel].mean(), mean_delay_se=delay[sel].std(ddof=1) / math.sqrt(sel.sum()),
                                     median_delay=float(np.median(delay[sel])), detected_within_horizon=det[sel].mean()))
            log(f"delay {day} ratio {rho}: {time.time() - t0:.0f} s")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th6_delay.csv", index=False)
    pd.DataFrame(prow).to_csv(OUT / "th6_delay_by_position.csv", index=False)
    print(df.pivot_table(index=["day_len", "ratio"], columns="detector", values="mean_delay").round(1))
    print(pd.DataFrame(prow).pivot_table(index=["day_len", "ratio", "position"], columns="detector",
                                         values="mean_delay").round(1))


def part_renewal():
    """Proposition 2 of docs/theory_per_period.md: the delay of the detector restarted at every day boundary is at
    most E[min(D0, m - k)] + P(D0 > m - k) E[min(D0, m)] / P(D0 <= m) for a change at step k of a day (k = 0..m-1),
    with D0 the delay of the detector started at 0 at the change. D0 needs to be simulated over one day only.
    The bound is averaged over k uniform, as in part delay, and set beside the measured delays of th6_delay.csv.
    Added after the delay results were seen (post hoc), to explain them.
    """
    rows = []
    meas = pd.read_csv(OUT / "th6_delay.csv")
    tau, m0, _ = TRUTH["gaussian"]
    B = law_B("gaussian")
    for day, A in SETTINGS.items():
        for i, rho in enumerate(RATIOS):
            P = _n(20000)
            rng = np.random.default_rng(SEED + 400 + i + (10 if day == 1440 else 0))

            def y_fn(t, eng, rng=rng):
                z = rng.standard_normal(P) * rho
                return np.minimum(z * z, tau * tau) / m0

            res = Pm.run_engine(y_fn, P, day, B, [Pm.Spec("SR", "SR", A)], day)
            first = res["SR"]["first"]
            d0 = np.where(first >= 0, first + 1, day + 1).astype(float)  # day + 1 stands for "after the day"
            p = float(np.mean(d0 <= day))
            e_dr0 = float(np.mean(np.minimum(d0, day))) / p if p > 0 else math.inf
            ks = np.arange(day)
            rem = day - ks
            e_min = np.array([np.mean(np.minimum(d0, r)) for r in rem])
            tail = np.array([np.mean(d0 > r) for r in rem])
            bound_k = e_min + tail * e_dr0
            mm = meas[(meas.day_len == day) & (meas.ratio == rho)].set_index("detector")
            rows.append(dict(day_len=day, ratio=rho, paths=P, p_D0_within_day=p, mean_D0_given_within=float(np.mean(d0[d0 <= day])),
                             E_restarted_change_at_day_start=e_dr0, bound_mean_delay_uniform=float(bound_k.mean()),
                             bound_last_10pct=float(bound_k[ks >= 0.9 * day].mean()),
                             bound_first_half=float(bound_k[ks < 0.5 * day].mean()),
                             measured_restarted_mean=float(mm.loc["SR_restarted_daily", "mean_delay"]),
                             measured_SR_mean=float(mm.loc["SR", "mean_delay"])))
            log(f"renewal {day} ratio {rho}: p = {p:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th6_delay_renewal.csv", index=False)
    print(df.round(3).to_string())


# ----------------------------------------------------------------------------
# Q6e: what the tolerance hypothesis of Theorem D excludes
# ----------------------------------------------------------------------------


def part_tolerance():
    rows = []
    P = _n(2000)
    cells = [("garch", 0.05, 0.90), ("garch", 0.05, 0.94), ("garch", 0.03, 0.969), ("garch", 0.02, 0.9799),
             ("slow_sv", 0.15, 390.0)]
    for cal in (390, 1440):
        n_tau = int(round(0.25 * cal))
        n = cal + cal  # calibration, then one monitored day of the same length
        for i, (kind, a, b) in enumerate(cells):
            rng = np.random.default_rng(SEED + 300 + i + (5 if cal == 1440 else 0))
            sig = np.empty((P, n))
            r = np.empty((P, n))
            if kind == "garch":
                # the returns that drive the recursion are the returns monitored, so sigma_t is a function of
                # r_1, ..., r_{t-1} and mu_t is the conditional mean given the returns
                omega = 1 - a - b
                s2 = np.ones(P)
                rprev = rng.standard_normal(P)
                for t in range(5000 + n):  # burn-in, then keep sigma_t^2 (known at t - 1)
                    s2 = omega + a * rprev**2 + b * s2
                    rprev = np.sqrt(s2) * rng.standard_normal(P)
                    if t >= 5000:
                        sig[:, t - 5000] = np.sqrt(s2)
                        r[:, t - 5000] = rprev
                label = f"GARCH(1,1), alpha {a}, beta {b} (persistence {a + b:.4g})"
            else:
                # mu_t is the conditional mean given the returns and the volatility path up to t (a filtration
                # to which the returns are adapted, so Theorem 6.4 applies with it)
                phi = 0.5 ** (1 / b)
                h = rng.normal(0, a, P)
                for t in range(n):
                    h = phi * h + math.sqrt(1 - phi * phi) * a * rng.standard_normal(P)
                    sig[:, t] = np.exp(h)
                r[:] = sig * rng.standard_normal((P, n))
                label = f"log-normal SV, log-sd {a}, half-life {b:g} steps"

            tau = 3 * MAD_NORMAL_SCALE * np.median(np.abs(r[:, :n_tau]), axis=1)
            mu = sig**2 * trunc_mean_normal(tau[:, None] / sig)
            calp, mon = mu[:, n_tau:cal], mu[:, cal:]
            up = mon.max(axis=1) / calp.min(axis=1) - 1
            dn = 1 - mon.min(axis=1) / calp.max(axis=1)
            need = np.maximum(up, dn)
            row = dict(model=label, cal_len=cal, monitored=cal, paths=P,
                       b_up_median=np.median(up), b_up_q90=np.quantile(up, 0.9),
                       b_down_median=np.median(dn), b_down_q90=np.quantile(dn, 0.9),
                       b_two_sided_median=np.median(need), b_two_sided_q90=np.quantile(need, 0.9),
                       b_two_sided_q95=np.quantile(need, 0.95))
            for bb in (0.1, 0.25, 0.5, 1.0):
                row[f"p_holds_b{bb:g}"] = float(np.mean(need <= bb))
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th6_tolerance.csv", index=False)
    print(df.round(3).to_string())


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

INK, MUTED, GRIDC = "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a5cc2"]


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

    FIG.mkdir(exist_ok=True)
    hz = pd.read_csv(OUT / "th6_hazard_390.csv")
    adv = pd.read_csv(OUT / "th6_adversary.csv")
    fam = pd.read_csv(OUT / "th6_calm_family.csv")
    ps = pd.read_csv(OUT / "th6_park_strike.csv")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    A, m = 7800.0, 390
    # Panel 1: hazard of the plain Shiryaev-Roberts detector
    ax = axes[0]
    for (stream, lab), col in zip((("gaussian", "Gaussian"), ("student_t3", "Student t3"), ("jumps", "jumps"),
                                   ("ticks", "ticks")), SERIES):
        s = hz[(hz.stream == stream) & (hz.detector == "SR")]
        ax.plot(s.day, 100 * s.hazard, color=col, lw=1.3, marker="o", ms=3, label=lab)
    s = adv[(adv.day_len == m) & (adv.detector == "SR")]
    ax.plot(s.day, 100 * s.hazard, color=SERIES[4], lw=1.3, marker="^", ms=3.5, label="bold play (E16)")
    ax.axhline(100 * m / A, color=INK, lw=1, ls="--", label="m / A (day 1 only)")
    ax.set_title("Shiryaev-Roberts: P(alarm on day d | none before)", fontsize=9.5, loc="left", color=INK)
    ax.set_xlabel("day d", fontsize=9, color=INK)
    ax.set_ylabel("percent", fontsize=9, color=INK)
    ax.set_ylim(0, 7)
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    _style(ax)
    # Panel 2: the referee's family against the Ville line and the two per-period constructions
    ax = axes[1]
    f = fam[(fam.day_len == m) & (fam.law_on_day_d == "gaussian")]
    for (det, lab), col in zip((("SR", "Shiryaev-Roberts"), ("SR_restarted_daily", "restarted daily (Theorem 1)"),
                                ("CUSUM", "e-CUSUM (Theorem 2)")), SERIES):
        s = f[f.detector == det]
        ax.plot(s.day, 100 * s.p_alarm_day_d, color=col, lw=1.3, marker="o", ms=3, label=lab)
    ax.plot(f.day.unique(), 100 * np.minimum(1, f.day.unique() * m / A), color=MUTED, lw=1, ls=":", label="d m / A (Ville)")
    ax.axhline(100 * m / A, color=INK, lw=1, ls="--", label="m / A")
    ax.set_title("Calm for d - 1 days, Gaussian on day d", fontsize=9.5, loc="left", color=INK)
    ax.set_xlabel("day d", fontsize=9, color=INK)
    ax.set_ylabel("P(false alarm on day d), percent", fontsize=9, color=INK)
    ax.set_ylim(0, 70)
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    _style(ax)
    # Panel 3: e-CUSUM, i.i.d. against the park-and-strike stream
    ax = axes[2]
    for (stream, lab), col in zip((("gaussian", "Gaussian, i.i.d."), ("ticks", "ticks, i.i.d.")), SERIES):
        s = hz[(hz.stream == stream) & (hz.detector == "CUSUM_no_restart")]
        ax.plot(s.day, 100 * s.p_any, color=col, lw=1.3, marker="o", ms=3, label=lab)
    s = adv[(adv.day_len == m) & (adv.detector == "CUSUM")]
    ax.plot(s.day, 100 * s.p_any, color=SERIES[2], lw=1.3, marker="^", ms=3.5, label="bold play, restarts")
    s = ps[ps.day_len == m]
    ax.plot(s.day, 100 * s.p_alarm_day_d, color=SERIES[4], lw=1.3, marker="s", ms=3, label="park, then strike on day d")
    ax.axhline(100 * Pm.window_bound_two_sided(m, A), color=INK, lw=1, ls="--", label="Theorem 2 bound (i.i.d.)")
    ax.set_title("e-CUSUM: P(an alarm on day d)", fontsize=9.5, loc="left", color=INK)
    ax.set_xlabel("day d", fontsize=9, color=INK)
    ax.set_ylabel("percent", fontsize=9, color=INK)
    ax.set_ylim(0, 7)
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.legend(fontsize=7.5, frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.08))
    _style(ax)
    fig.suptitle("TH6: false alarms after the first day, default two-sided detector, true baseline, A = 7800, "
                 "390-step days (synthetic)", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "th6_hazard.png", dpi=140)
    plt.close(fig)

    dl = pd.read_csv(OUT / "th6_delay.csv")
    dp = pd.read_csv(OUT / "th6_delay_by_position.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    labs = (("SR_restarted_daily", "restarted daily, A"), ("SR_restarted_carry_m", "restarted, carry m, 2A"),
            ("CUSUM_certified", "e-CUSUM, certified A"))
    ax = axes[0]
    x = np.arange(len(RATIOS))
    for k, ((det, lab), col) in enumerate(zip(labs, SERIES)):
        for day, mk, ls in ((390, "o", "-"), (1440, "s", "--")):
            base = dl[(dl.day_len == day) & (dl.detector == "SR")].set_index("ratio").mean_delay
            s = dl[(dl.day_len == day) & (dl.detector == det)].set_index("ratio").mean_delay
            ax.plot(x, (s / base).loc[list(RATIOS)].values, color=col, marker=mk, ls=ls, lw=1.2, ms=4,
                    label=f"{lab}, {day}-step days")
    ax.axhline(1, color=INK, lw=0.8)
    ax.set_xticks(x, [f"x{r:g}" for r in RATIOS])
    ax.set_xlabel("change of the return standard deviation", fontsize=9, color=INK)
    ax.set_ylabel("mean delay / mean delay of Shiryaev-Roberts", fontsize=9, color=INK)
    ax.set_title("Delay cost of a per-period guarantee (change on day 3)", fontsize=9.5, loc="left", color=INK)
    ax.set_ylim(0.95, 2.3)
    ax.legend(fontsize=7, frameon=False, loc="upper center", ncol=2)
    _style(ax)
    ax = axes[1]
    s = dp[(dp.day_len == 390)]
    positions = ["first half", "50% to 90%", "last 10%"]
    for (det, lab), col in zip((("SR", "Shiryaev-Roberts, A"),) + labs, [INK] + SERIES):
        v = [s[(s.detector == det) & (s.position == p) & (s.ratio == 1.5)].mean_delay.iloc[0] for p in positions]
        ax.plot(positions, v, color=col, marker="o", lw=1.2, ms=4, label=lab)
    ax.set_ylabel("mean delay (steps)", fontsize=9, color=INK)
    ax.set_xlabel("position of the change within the day", fontsize=9, color=INK)
    ax.set_title("sd x1.5, 390-step days, by position of the change", fontsize=9.5, loc="left", color=INK)
    ax.legend(fontsize=7.5, frameon=False)
    _style(ax)
    fig.suptitle("TH6: delay cost of a proved bound for every day (synthetic, Gaussian returns, true baseline)",
                 fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG / "th6_delay.png", dpi=140)
    plt.close(fig)


PARTS = {"laws": part_laws, "hazard390": lambda: part_hazard(390), "hazard1440": lambda: part_hazard(1440),
         "family": part_family, "adversary": part_adversary, "delay": part_delay, "renewal": part_renewal,
         "tolerance": part_tolerance, "figures": figures}

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for p in sys.argv[1:] or ["laws"]:
        log("part", p)
        PARTS[p]()
