"""Per-period false alarm guarantees (TH6; theory and proofs in ``docs/theory_per_period.md``).

Theorem B (``docs/math.md``, 5.4) bounds the probability of a false alarm within the FIRST m steps after
the statistics start (at the start of monitoring or after a recalibration) by m / A. It does not bound
the probability of a false alarm on a later day of the same run: under the null the Shiryaev-Roberts
statistic M_t has E[M_t] <= t, so on day d Ville's inequality gives only d m / A. Corollary C bounds the
expected number of false alarms per monitored step by 1 / A when the monitor restarts after every alarm.

This module gives two constructions with a guarantee for every period, and the tools to check them:

* **Restarted Shiryaev-Roberts** (``period`` in ``Spec`` or ``PeriodicVolatilityMonitor``): both
  statistics restart at 0 at every period boundary. With a known baseline (the engine below, or
  ``PeriodicVolatilityMonitor(baseline="known", tau=..., m0=...)``), P(false alarm in a period | past)
  <= m / A for every period of m steps, under every null law of Theorem B (``restarted_sr_bound``).
  Keeping at most ``carry`` of each side's average across a boundary gives (m + carry) / A instead. With
  the confidence-sequence baseline only the marginal bound delta + m / A holds for each period of a
  segment, and only while the tolerance of Theorem 6.4 holds; with the plug-in baseline nothing per
  period is proved.
* **CUSUM-type e-detector** (e-CUSUM; ``kind="CUSUM"``): C_t = f_t max(C_{t-1}, 1), the maximum over start
  times of the product of bets. Under an i.i.d. (or exchangeable) null with fixed bets, time reversal
  gives P(C_s >= x) <= 1 / x at every s, so a single bet has P(alarm in any window of m steps) <= m / A
  (``window_bound_single``; an unconditional bound, not one given the past), and the default average
  over bets has
  (m + log A + H(w)) / A, H the entropy of the bet weights (``window_bound_mixture``,
  ``window_bound_two_sided``). Under a martingale-difference null the reversal fails and no such bound
  holds (``park_strike_policy`` builds the counterexample). ``conditional_window_bound`` is the part
  that holds under every null of Theorem B.

The simulation engine ``run_engine`` runs several detectors on the same paths of Y = X / m0 (true
baseline, so the null mean of Y is 1 and Y lies in [0, B]) and records the first alarm and the number of
alarms per day. Policies (``bold_policy``, ``park_strike_policy``) are adaptive null streams that watch
one detector's state; every draw has conditional mean exactly 1.

Nothing here changes ``VolatilityMonitor`` or any existing function.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import optimize

from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from .monitor import Event, VolatilityMonitor

LAND = 1e-9  # landing draws aim at 2A (1 + LAND), so floating point never falls just short (as in E16)


# ----------------------------------------------------------------------------
# Closed-form bounds
# ----------------------------------------------------------------------------


def entropy(weights) -> float:
    """Shannon entropy -sum w log w (natural log) of bet weights that sum to 1; log K for K equal weights."""
    w = np.asarray(weights, float)
    w = w[w > 0]
    if not np.isclose(w.sum(), 1.0):
        raise ValueError("weights must sum to 1")
    return float(-(w * np.log(w)).sum())


def restarted_sr_bound(m: int, A: float, carry: float = 0.0) -> float:
    """Theorem 1: per-period bound of the restarted Shiryaev-Roberts detector, (m + carry) / A, capped at 1.

    ``carry`` is the largest two-sided average N = (M_up + M_down) / 2 kept across a boundary (0 = restart).
    Holds conditionally on the past at the boundary, for every null law of Theorem B.
    """
    return min(1.0, (m + carry) / A)


def window_bound_single(m: int, A: float) -> float:
    """Theorem 2(a): one fixed bet, CUSUM statistic, alarm at A, exchangeable null: any window of m steps."""
    return min(1.0, m / A)


def window_bound_mixture(m: int, A: float, weights) -> float:
    """Theorem 2(b): average of per-bet CUSUM statistics with weights w, alarm at A, exchangeable null.

    P(an alarm at some step of any window of m steps) <= (m + log A + H(w)) / A. A union bound over the
    m steps of per-step tails gives m (1 + log A + H(w)) / A, never better for m >= 1; the smaller is
    returned (capped at 1).
    """
    h = entropy(weights)
    a = (m + math.log(A) + h) / A
    b = m * (1 + math.log(A) + h) / A
    return min(1.0, a, b)


def window_bound_two_sided(m: int, A: float, n_up: int = len(DEFAULT_UP_GRID), n_down: int = len(DEFAULT_DOWN_FRACTIONS),
                           w_up=None, w_down=None) -> float:
    """Theorem 2(c): default two-sided e-CUSUM, each side the average of its bets, each side alarming at 2A.

    Union over the two sides of Theorem 2(b) at level 2A:
    (m + log(2A) + (H_up + H_down) / 2) / A. For the default grid (6 + 5 equal weights) and A = 7800,
    m = 390 this is 5.146%; for A = 28,800, m = 1440 it is 5.044%.
    """
    w_up = np.full(n_up, 1.0 / n_up) if w_up is None else np.asarray(w_up, float)
    w_down = np.full(n_down, 1.0 / n_down) if w_down is None else np.asarray(w_down, float)
    return min(1.0, window_bound_mixture(m, 2 * A, w_up) + window_bound_mixture(m, 2 * A, w_down))


def threshold_for_window(m: int, alpha: float, bound=window_bound_two_sided, **kw) -> float:
    """Smallest A with bound(m, A, **kw) <= alpha (the bound is decreasing in A)."""
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0, 1)")
    f = lambda A: bound(m, A, **kw) - alpha  # noqa: E731
    lo, hi = m / alpha * 0.5, m / alpha * 4.0 + 100.0
    while f(hi) > 0:
        hi *= 2
    return float(optimize.brentq(f, lo, hi, xtol=1e-9))


def conditional_window_bound(c_up, c_down, m: int, A: float, w_up=None, w_down=None) -> np.ndarray:
    """Lemma 2: P(an alarm in the next m steps | F_s) for the two-sided e-CUSUM (each side at 2A).

    ``c_up`` (paths x K_up) and ``c_down`` are the per-bet CUSUM statistics at time s. Holds under every
    null of Theorem B (martingale differences), with no exchangeability:
    sum over sides of min(1, (sum_k w_k max(C_k, 1) + m - 1) / (2 A)), capped at 1.
    """
    cu = np.atleast_2d(np.asarray(c_up, float))
    cd = np.atleast_2d(np.asarray(c_down, float))
    wu = np.full(cu.shape[1], 1.0 / cu.shape[1]) if w_up is None else np.asarray(w_up, float)
    wd = np.full(cd.shape[1], 1.0 / cd.shape[1]) if w_down is None else np.asarray(w_down, float)
    bu = np.minimum(1.0, (np.maximum(cu, 1.0) @ wu + m - 1) / (2 * A))
    bd = np.minimum(1.0, (np.maximum(cd, 1.0) @ wd + m - 1) / (2 * A))
    return np.minimum(1.0, bu + bd)


def reversed_running_max(f: np.ndarray) -> np.ndarray:
    """Lemma 1 as an identity on one path: C_s = max_{j <= s} prod_{i=j}^{s} f_i equals the running maximum
    max_{k <= s} prod_{i=1}^{k} f_{s+1-i} of the forward partial products of the reversed sequence.
    Returns that maximum for s = len(f) (used in tests)."""
    g = np.asarray(f, float)[::-1]
    return float(np.max(np.cumprod(g)))


# ----------------------------------------------------------------------------
# Vectorised engine on Y = X / m0 (true baseline)
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Spec:
    """One detector in ``run_engine``.

    kind: "SR" (Shiryaev-Roberts, R_t = f_t (R_{t-1} + 1)) or "CUSUM" (C_t = f_t max(C_{t-1}, 1)).
    A: each side alarms when the average over its bets reaches 2A (the monitor's two-sided rule).
    period: restart both sides at every multiple of ``period`` on the stream clock (None: never).
    carry: at a boundary, scale each side's bets down so that the side average is at most ``carry``
        (0 restarts at 0).
    on_alarm: "restart" (statistics back to 0, as a recalibration with the true baseline) or "continue"
        (the statistic keeps running; counts then record steps at or above the threshold).
    """

    name: str
    kind: str = "SR"
    A: float = 7800.0
    period: int | None = None
    carry: float = 0.0
    on_alarm: str = "restart"


def default_bets(B: float, up=DEFAULT_UP_GRID, down_fractions=DEFAULT_DOWN_FRACTIONS):
    """The monitor's default increase bets and decrease bets (fractions of 1 / (B - 1))."""
    if B <= 1:
        raise ValueError("B must exceed 1")
    return np.asarray(up, float), np.asarray(down_fractions, float) / (B - 1.0)


class Engine:
    """State of several detectors on P paths; see ``run_engine``."""

    def __init__(self, specs, P: int, B: float, up=DEFAULT_UP_GRID, down_fractions=DEFAULT_DOWN_FRACTIONS):
        self.specs = {s.name: s for s in specs}
        if len(self.specs) != len(specs):
            raise ValueError("detector names must be unique")
        for s in specs:
            if s.kind not in ("SR", "CUSUM"):
                raise ValueError("kind must be 'SR' or 'CUSUM'")
            if s.on_alarm not in ("restart", "continue"):
                raise ValueError("on_alarm must be 'restart' or 'continue'")
        self.P, self.B = int(P), float(B)
        self.up, self.dn = default_bets(B, up, down_fractions)
        self.wu = np.full(self.up.size, 1.0 / max(self.up.size, 1))  # an empty side has statistic 0
        self.wd = np.full(self.dn.size, 1.0 / max(self.dn.size, 1))
        self.Su = {s.name: np.zeros((self.P, self.up.size)) for s in specs}
        self.Sd = {s.name: np.zeros((self.P, self.dn.size)) for s in specs}

    def boundary(self, t: int):
        for s in self.specs.values():
            if s.period and t > 0 and t % s.period == 0:
                for S, w in ((self.Su[s.name], self.wu), (self.Sd[s.name], self.wd)):
                    if s.carry <= 0:
                        S[:] = 0.0
                    else:
                        avg = S @ w
                        big = avg > s.carry
                        if big.any():
                            S[big] *= (s.carry / avg[big])[:, None]

    def landing(self, name: str):
        """Bold play against detector ``name``: the high draw that takes the increase side exactly to 2A
        (or B if it cannot reach) and the low draw that takes the decrease side to 2A (or 0)."""
        s = self.specs[name]
        tgt = 2 * s.A * (1 + LAND)
        gu = self.Su[name] + 1.0 if s.kind == "SR" else np.maximum(self.Su[name], 1.0)
        gd = self.Sd[name] + 1.0 if s.kind == "SR" else np.maximum(self.Sd[name], 1.0)
        au, bu = gu @ self.wu, gu @ (self.wu * self.up)
        ad, bd = gd @ self.wd, gd @ (self.wd * self.dn)
        with np.errstate(divide="ignore", invalid="ignore"):
            ysu = 1.0 + (tgt - au) / bu
            ysd = 1.0 - (tgt - ad) / bd
        up_ok = (ysu <= self.B) & (bu > 0)
        dn_ok = (ysd >= 0.0) & (bd > 0)
        yh = np.where(up_ok, np.maximum(ysu, 1.0), self.B)
        yl = np.where(dn_ok, np.minimum(ysd, 1.0), 0.0)
        return yl, yh, up_ok, dn_ok

    def step(self, y: np.ndarray):
        """Update every detector with the draws y (P,); returns {name: alarm mask}."""
        g = y - 1.0
        fu = 1.0 + self.up[None, :] * g[:, None]
        fd = 1.0 - self.dn[None, :] * g[:, None]
        out = {}
        for s in self.specs.values():
            Su, Sd = self.Su[s.name], self.Sd[s.name]
            if s.kind == "SR":
                Su += 1.0
                Sd += 1.0
            else:
                np.maximum(Su, 1.0, out=Su)
                np.maximum(Sd, 1.0, out=Sd)
            Su *= fu
            Sd *= fd
            thr = 2 * s.A
            hit = (Su @ self.wu >= thr) | (Sd @ self.wd >= thr)
            if s.on_alarm == "restart" and hit.any():
                Su[hit] = 0.0
                Sd[hit] = 0.0
            out[s.name] = hit
        return out


def run_engine(y_fn, P: int, n: int, B: float, specs, day: int, t0: int = 0, init=None, after=None,
               up=DEFAULT_UP_GRID, down_fractions=DEFAULT_DOWN_FRACTIONS):
    """Run detectors ``specs`` on P paths for steps t0, ..., t0 + n - 1 of the stream clock.

    ``y_fn(t, eng)`` returns the draws Y_t (P,), mean-1 law under the null; it may read ``eng`` (an adaptive
    stream). ``day`` is the period length used for the per-day counts. ``init`` maps a detector name to
    (S_up, S_down) initial per-bet statistics (for example the deterministic state after calm days).
    ``after`` (P,) optional: also record, per path, the first alarm at a step t >= after[path] ("first_after",
    for delays after a change). Returns {name: {"first": first alarm step or -1, "counts": alarms per day
    (P, days)}}; days are counted from t0 // day.
    """
    eng = Engine(specs, P, B, up, down_fractions)
    if init:
        for name, (su, sd) in init.items():
            eng.Su[name][:] = su
            eng.Sd[name][:] = sd
    d0 = t0 // day
    n_days = (t0 + n - 1) // day - d0 + 1
    res = {s.name: {"first": np.full(P, -1, np.int64), "counts": np.zeros((P, n_days), np.int32)} for s in specs}
    if after is not None:
        after = np.asarray(after)
        for s in specs:
            res[s.name]["first_after"] = np.full(P, -1, np.int64)
    for t in range(t0, t0 + n):
        eng.boundary(t)
        y = np.asarray(y_fn(t, eng), float)
        hits = eng.step(y)
        di = t // day - d0
        for name, hit in hits.items():
            if hit.any():
                r = res[name]
                r["counts"][hit, di] += 1
                newly = hit & (r["first"] < 0)
                r["first"][newly] = t
                if after is not None:
                    newly = hit & (r["first_after"] < 0) & (after <= t)
                    r["first_after"][newly] = t
    return res


def two_point(u: np.ndarray, yl: np.ndarray, yh: np.ndarray) -> np.ndarray:
    """Draw from the law on {yl, yh} with mean 1 (yl <= 1 <= yh); degenerate at 1 when yl = yh = 1."""
    span = yh - yl
    p = np.where(span > 0, (1.0 - yl) / np.where(span > 0, span, 1.0), 0.0)
    return np.where(u < p, yh, yl)


def bold_policy(name: str, rng: np.random.Generator):
    """E16's bold play against detector ``name``: spread {0, B}, landing exactly on 2A when one draw can."""

    def y_fn(t, eng):
        yl, yh, _, _ = eng.landing(name)
        return two_point(rng.random(eng.P), yl, yh)

    return y_fn


def park_strike_policy(name: str, strike_start: int, rng: np.random.Generator):
    """Counterexample of Proposition 3 against a CUSUM detector ``name`` (a martingale-difference null).

    Before ``strike_start``: draw from {0, B} (mean 1) until one draw could take a side to 2A, then play
    Y = 1, which leaves every CUSUM statistic at max(C, 1) (it parks). From ``strike_start`` on: bold
    play with landing. Every draw has conditional mean exactly 1; no alarm can occur before the strike.
    """

    def y_fn(t, eng):
        yl, yh, up_ok, dn_ok = eng.landing(name)
        u = rng.random(eng.P)
        if t >= strike_start:
            return two_point(u, yl, yh)
        y = two_point(u, np.zeros(eng.P), np.full(eng.P, eng.B))
        return np.where(up_ok | dn_ok, 1.0, y)

    return y_fn


def day_table(first: np.ndarray, counts: np.ndarray) -> dict:
    """Per-day summaries from ``run_engine`` output for one detector.

    p_first[d]: P(first alarm on day d); hazard[d]: P(alarm on day d | none before), with the number at risk;
    any[d]: P(at least one alarm on day d) (with restarts after alarms); mean_count[d]: expected number of
    alarms on day d; cum_rate[d]: expected alarms per day over days 1..d (Corollary C bounds it by m / A).
    """
    P, D = counts.shape
    day_of_first = np.full(P, -1)
    has = counts.sum(axis=1) > 0
    day_of_first[has] = np.argmax(counts[has] > 0, axis=1)
    p_first = np.array([np.mean(day_of_first == d) for d in range(D)])
    at_risk = np.array([np.sum((day_of_first < 0) | (day_of_first >= d)) for d in range(D)])
    hazard = np.array([np.sum(day_of_first == d) / at_risk[d] if at_risk[d] else np.nan for d in range(D)])
    anyd = (counts > 0).mean(axis=0)
    mean_count = counts.mean(axis=0)
    cum_rate = np.cumsum(mean_count) / np.arange(1, D + 1)
    return dict(p_first=p_first, at_risk=at_risk, hazard=hazard, any=anyd, mean_count=mean_count,
                count_sd=counts.std(axis=0, ddof=1) if P > 1 else np.zeros(D), cum_rate=cum_rate)


# ----------------------------------------------------------------------------
# A streaming monitor with a per-period option
# ----------------------------------------------------------------------------


class PeriodicVolatilityMonitor(VolatilityMonitor):
    """``VolatilityMonitor`` with options for a guarantee in every period (``docs/theory_per_period.md``).

    period : int or None
        Restart both statistics at 0 at every step t of the stream with (t - offset) % period == 0 (for
        example the first minute of each trading day). ``carry > 0`` keeps at most ``carry`` of each side's
        average across a boundary instead (bound (period + carry) / A below). What this buys depends on
        the baseline:

        * ``baseline="known"`` (with ``tau`` and ``m0``): P(false alarm in a period | everything before it)
          <= period / A for every period, under the null of Theorem B for that tau and m0 (Theorem 1). The
          monitor then never calibrates: it starts at the first return and restarts at once after an alarm.
        * ``baseline="cs"``: P(false alarm in a given period of a segment) <= delta + period / A, marginally
          (not given the past), if the tolerance of Theorem 6.4 holds up to the end of that period
          (Corollary 1 of the theory note).
        * ``baseline="plugin"`` (the default): no per-period bound is proved.
    statistic : "sr" (default) or "cusum"
        "cusum" replaces the Shiryaev-Roberts recursion by C_t = f_t max(C_{t-1}, 1). With a known baseline,
        fixed bets and an i.i.d. or exchangeable null, the probability of a false alarm in any window of m
        steps is at most ``window_bound_two_sided(m, A)`` (Theorem 2; unconditional, not given the past).
        Under a martingale-difference null it has no such bound (Proposition 3).
    baseline : "plugin" (default), "cs" or "known"
        "known" needs ``tau`` (the truncation level) and ``m0`` (the null mean of min(r^2, tau^2)), and
        optionally ``tau_down`` and ``m0_down`` for the decrease side (default: the same values).

    With ``period=None``, ``statistic="sr"`` and a plug-in or confidence-sequence baseline the output is
    identical to ``VolatilityMonitor``. Calibration and the two-sided rule are inherited unchanged.
    """

    def __init__(self, *args, period: int | None = None, offset: int = 0, carry: float = 0.0,
                 statistic: str = "sr", baseline: str = "plugin", tau: float | None = None,
                 m0: float | None = None, tau_down: float | None = None, m0_down: float | None = None, **kw):
        self._known = None
        if baseline == "known":
            if tau is None or m0 is None:
                raise ValueError("baseline='known' needs tau and m0")
            if kw.get("margin_z") or kw.get("target_rel_width") is not None:
                raise ValueError("margin_z and target_rel_width do not apply to a known baseline")
            td = float(tau) if tau_down is None else float(tau_down)
            md = float(m0) if m0_down is None else float(m0_down)
            for name, (tt, mm) in (("increase", (float(tau), float(m0))), ("decrease", (td, md))):
                if not (math.isfinite(tt) and tt > 0 and math.isfinite(mm) and 0 < mm <= tt * tt):
                    raise ValueError(f"{name} side: need tau > 0 and 0 < m0 <= tau^2")
            known = (float(tau), float(m0), td, md)
            super().__init__(*args, baseline="plugin", **kw)
            self.baseline = "known"
            self._known = known
            self._reset()
        else:
            if any(v is not None for v in (tau, m0, tau_down, m0_down)):
                raise ValueError("tau, m0, tau_down and m0_down need baseline='known'")
            super().__init__(*args, baseline=baseline, **kw)
        if period is not None and int(period) < 1:
            raise ValueError("period must be a positive integer")
        if statistic not in ("sr", "cusum"):
            raise ValueError("statistic must be 'sr' or 'cusum'")
        if carry < 0:
            raise ValueError("carry must be nonnegative")
        self.period = None if period is None else int(period)
        self.offset = int(offset)
        self.carry = float(carry)
        self.statistic = statistic
        self.period_resets = 0

    def _reset(self):
        super()._reset()
        if getattr(self, "_known", None) is not None:  # known baseline: monitor from the next return on
            tau, m0, td, md = self._known
            self.scale = tau / self.c
            self.tau2 = tau * tau
            self.tau2_dn = td * td
            self._n_finite = 0
            self._finish(m0, m0 * (1 + self.band), md * (1 - self.band), None, None, md)

    def _boundary(self):
        if self.period is None or (self.t - self.offset) % self.period != 0 or self.t <= self.start:
            return
        self.period_resets += 1
        for attr in ("S_up", "S_dn"):
            S = getattr(self, attr)
            if self.carry <= 0:
                setattr(self, attr, np.zeros_like(S))
            else:
                avg = S.mean()
                if avg > self.carry:
                    setattr(self, attr, S * (self.carry / avg))

    def update(self, r: float) -> Event | None:
        if self.period is None and self.statistic == "sr":
            return super().update(r)
        self.t += 1
        if not self.calibrated:
            self._feed_calibration(r)
            return None
        self._boundary()
        if not np.isfinite(r):
            return None
        x = min(r * r, self.tau2)
        fu = 1.0 + self.up_grid * (x / self.m_up - 1.0)
        if self.statistic == "sr":
            self.S_up = fu * (self.S_up + 1.0)
        else:
            self.S_up = fu * np.maximum(self.S_up, 1.0)
        if self.down_active:
            xd = x if self.tau2_dn == self.tau2 else min(r * r, self.tau2_dn)
            fd = 1.0 + self.lam_dn * (1.0 - xd / self.m_dn)
            if self.statistic == "sr":
                self.S_dn = fd * (self.S_dn + 1.0)
            else:
                self.S_dn = fd * np.maximum(self.S_dn, 1.0)
        mu, md = self.S_up.mean(), self.S_dn.mean()
        ev = None
        if mu >= 2 * self.A:
            ev = Event(self.t, "up", self.start, self.scale, float(mu))
        elif md >= 2 * self.A:
            ev = Event(self.t, "down", self.start, self.scale, float(md))
        if ev is not None:
            self.events.append(ev)
            self._reset()
        return ev
