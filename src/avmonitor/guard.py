"""A retry guard against blind segments of the confidence-sequence baseline (E24).

Theory and proofs: ``docs/theory_blind_guard.md``. Nothing in ``monitor.py`` or ``periodic.py`` is changed; the
guarded monitors are subclasses that override only the end of calibration, the reset after an alarm and (for the
rank fallback) the monitoring step.

The blind condition
-------------------
With ``baseline="cs"`` the truncation level tau is fixed from the first quarter of the calibration window and the
upper bound U from the rest; the increase side bets against ``m_up = (1 + b) U``. If ``tau^2 <= m_up``, every
increase factor ``1 + lam (X / m_up - 1)`` with ``X = min(r^2, tau^2)`` is at most 1, so the increase side cannot
alarm for the life of the segment (with the daily restart of P1, never; without it, not within 2A - 1 monitored steps).
Both tau and U are known when calibration ends, so the condition is decided at a stopping time. The guard uses the
"near-blind" form ``tau^2 <= (1 + eta) m_up`` with ``eta >= 0`` fixed in advance (``eta = 0`` is E19's criterion).

The guard
---------
If the window just calibrated is (near-)blind and fewer than ``guard_tries = K`` windows have been tried in this
segment, the window is discarded and a fresh calibration window starts at the next return. Window k uses its own
error budget ``delta_k`` (split equally between U and L, as in the unguarded monitor), with
``delta_1 + ... + delta_K <= delta``. After K blind windows the fallback is either ``"keep"`` (monitor with the
K-th window, as the unguarded monitor would) or ``"rank"`` (E14's awake rank monitor, awake-B, for the rest of the
segment, at the same A and band; after its alarm the next segment calibrates again with window 1).

Guarantee (``docs/theory_blind_guard.md``, Theorems 1 to 3), with ``delta`` the total budget:

* P(the bounds of the window kept in a segment fail) <= delta_1 + ... + delta_K <= delta, given everything before
  the segment. Not per window: given that a window is kept its bounds can fail with probability up to 1
  (Proposition 3), so each window tried must have its own share of delta.
* Hence ``docs/math.md`` Theorem 6.4 and Corollary 6.6, and Corollary 1 of ``docs/theory_per_period.md`` (with the
  daily restart), hold unchanged, with T the end of the kept window and the tolerance (H-up), (H-down) relative to
  the kept window's confidence-sequence part (sufficient, and independent of the selection: the tolerance relative
  to each of the K windows, docs/theory_blind_guard.md (H*)): P(an alarm within the first m monitored steps of a
  segment) <= delta + m / A; with ``period=m``, P(an alarm in a given period of a segment) <=
  delta + (m + carry) / A (marginal); expected alarms in the first M monitored steps <= (delta + M / A) / (1 - delta).
  Calibration steps, including discarded windows, are not monitored steps.
* Rank fallback: the same bounds hold if, in a segment that falls back, the returns after the fallback time are
  exchangeable given everything up to it (Theorem 3; this replaces the tolerance for that segment only).

Which window is kept depends on delta_k, the bets and eta, so the null "tolerance relative to the kept window"
depends on the guard's tuning. A sufficient condition that does not depend on which window is kept is the
tolerance relative to each of the K windows ((H*) of the theory note), which implies both the guarded and the
unguarded hypotheses; it holds under a constant conditional law of the returns from the end of the first window's
tau part on (``docs/math.md``, Corollary 6.5), for example i.i.d. returns with any law, possibly after a step in
volatility inside the first window's tau part.

With ``guard_tries=1`` and ``guard_fallback="keep"`` the guarded monitor is the unguarded one, step for step.

The constructors' defaults (``guard_eta=0``, ``guard_fallback="keep"``) are E19's criterion (E24's arm G0), not the
recommended guard, which is ``guard_eta=0.25, guard_fallback="rank"`` (``docs/theory_blind_guard.md`` 6.1).
"""

from __future__ import annotations

import math

import numpy as np

from .monitor import Event, VolatilityMonitor
from .periodic import PeriodicVolatilityMonitor


def spending_schedule(delta: float, tries: int, kind: str = "halving") -> tuple[float, ...]:
    """Error budgets delta_1, ..., delta_K for the K windows a segment may try; they sum to ``delta``.

    ``"halving"``: delta / 2, delta / 4, ..., delta / 2^(K-1), and the last window gets the same as the one before
    it (so the sum is exactly delta; K = 1 gives delta, the unguarded monitor). ``"equal"``: delta / K each.
    """
    K = int(tries)
    if K < 1:
        raise ValueError("tries must be at least 1")
    if not 0 < delta < 1:
        raise ValueError("delta must lie in (0, 1)")
    if kind == "halving":
        if K == 1:
            return (float(delta),)
        d = [delta * 2.0 ** -(k + 1) for k in range(K - 1)]
        return tuple(d + [d[-1]])
    if kind == "equal":
        return tuple([delta / K] * K)
    raise ValueError("kind must be 'halving' or 'equal'")


def near_blind_delay_floor(rho: float, A: float, up_grid=None, cap: int = 10**7) -> float:
    """Proposition 2(b): the fewest monitored steps after a (re)start before the increase side can alarm, over all data.

    ``rho = tau^2 / m_up``. Every increase factor is at most ``1 + lam (rho - 1)``, attained when ``X = tau^2``, and the
    Shiryaev-Roberts recursion is increasing in each factor, so the largest side average after d steps is
    ``mean_k sum_{i=1}^{d} (1 + lam_k (rho - 1))^i``. Returns the first d at which it reaches 2A (``inf`` when
    ``rho < 1``; exactly 2A steps when ``rho = 1``; capped at ``cap``)."""
    from .detectors import DEFAULT_UP_GRID

    lam = np.asarray(DEFAULT_UP_GRID if up_grid is None else up_grid, float)
    if rho < 1:
        return math.inf
    g = 1.0 + lam * (rho - 1.0)
    R = np.zeros_like(g)
    for d in range(1, cap + 1):
        R = g * (R + 1.0)
        if R.mean() >= 2 * A:
            return float(d)
    return math.inf


class _GuardMixin:
    """The guard, shared by the plain and the periodic monitor (see the module docstring)."""

    def _guard_init(self, guard_tries, guard_deltas, guard_eta, guard_fallback, guard_spending, rank_seed):
        if self.baseline != "cs":
            raise ValueError("the guard needs baseline='cs'")
        K = int(guard_tries)
        if guard_deltas is None:
            deltas = spending_schedule(self.delta, K, guard_spending)
        else:
            deltas = tuple(float(d) for d in guard_deltas)
            if len(deltas) != K:
                raise ValueError("guard_deltas must have one entry per try")
        if any(not 0 < d < 1 for d in deltas):
            raise ValueError("every delta_k must lie in (0, 1)")
        if sum(deltas) > self.delta * (1 + 1e-12):
            raise ValueError("the budgets delta_k must sum to at most delta")
        if not (math.isfinite(guard_eta) and guard_eta >= 0):
            raise ValueError("guard_eta must be a nonnegative number")
        if guard_fallback not in ("keep", "rank"):
            raise ValueError("guard_fallback must be 'keep' or 'rank'")
        self.delta_total = float(self.delta)
        self.guard_tries = K
        self.guard_deltas = deltas
        self.guard_eta = float(guard_eta)
        self.guard_fallback = guard_fallback
        self._rng = np.random.default_rng(rank_seed)
        self.discarded: list[dict] = []
        self.guard_stats = {"windows_tried": 0, "windows_discarded": 0, "fallback_keep": 0, "fallback_rank": 0}
        self._rank = None
        self._try = 1
        self._seg_cal_start = self.t + 1
        self.delta = deltas[0]

    # A segment starts with window 1 and the budget delta_1 (after an alarm, and at the start of the stream).
    def _reset(self):
        super()._reset()
        if getattr(self, "guard_deltas", None) is not None:
            self._try = 1
            self.delta = self.guard_deltas[0]
            self._rank = None
            self._seg_cal_start = self.t + 1

    def _blind(self, m_up: float) -> bool:
        return bool(m_up > 0 and self.tau2 / m_up <= 1.0 + self.guard_eta)

    def _finish(self, m0, m_up, m_dn, upper, lower, m0_dn):
        self.guard_stats["windows_tried"] += 1
        ratio = self.tau2 / m_up if m_up > 0 else math.inf
        blind = self._blind(m_up)
        info = dict(try_index=self._try, delta_window=self.delta, tau2=self.tau2, m_up=m_up, tau2_over_m_up=ratio,
                    window_start=self.t + 1 - len(self.buffer), window_end=self.t, window_len=len(self.buffer))
        if blind and self._try < self.guard_tries:
            # discard: a fresh window starts at the next return, with the next budget (Theorem 1)
            self.guard_stats["windows_discarded"] += 1
            self.discarded.append(info)
            k = self._try + 1
            VolatilityMonitor._reset(self)  # clears the buffer and both sequences; keeps t, events and segments
            self._try = k
            self.delta = self.guard_deltas[k - 1]
            return
        if blind and self.guard_fallback == "rank":
            self.guard_stats["fallback_rank"] += 1
            self._start_rank(info)
            return
        if blind:
            self.guard_stats["fallback_keep"] += 1
        super()._finish(m0, m_up, m_dn, upper, lower, m0_dn)
        self.segments[-1].update(guard="blind kept" if blind else "kept", try_index=self._try,
                                 delta_window=self.delta, guard_cal_steps=self.t + 1 - self._seg_cal_start)

    # ------------------------------------------------------------------
    # Rank fallback (E14's awake rank monitor, awake-B, started empty at the fallback time)
    # ------------------------------------------------------------------

    def _start_rank(self, info):
        from .running import ConformalMonitor

        rm = ConformalMonitor(threshold=self.A, cal_len=0, family="lr", band=self.band, seed=self._rng, restart_len=0)
        rm.t = self.t
        rm.start = self.t + 1
        rm.segments[-1]["start"] = rm.start
        self._rank = rm
        self.calibrated = True
        self.start = self.t + 1
        self.scale = info["tau2"] ** 0.5 / self.c
        self.cal_stats["calibrations"] += 1
        nan = float("nan")
        self.segments.append(dict(start=self.start, cal_len=len(self.buffer), n_finite=self._n_finite, baseline="rank",
                                  n_tau=None, scale=self.scale, c=self.c, c_down=self.c_down, tau2=nan, tau2_dn=nan,
                                  m0=nan, m0_dn=nan, m_up=nan, m_dn=nan, upper=None, lower=None,
                                  extended=self._extended, fallback=self._fallback, down_active=True,
                                  guard="rank fallback", try_index=self._try, delta_window=0.0,
                                  guard_cal_steps=self.t + 1 - self._seg_cal_start, blind_tau2_over_m_up=info["tau2_over_m_up"]))

    def _rank_boundary(self) -> bool:
        period = getattr(self, "period", None)
        if period is None or (self.t - self.offset) % period != 0 or self.t <= self.start:
            return False
        return True

    def update(self, r: float) -> Event | None:
        if self._rank is None:
            return super().update(r)
        self.t += 1
        rm = self._rank
        if self._rank_boundary():  # the daily restart of P1 applies to the fallback's statistics too
            self.period_resets += 1
            rm.S_up = np.zeros_like(rm.S_up)
            rm.S_dn = np.zeros_like(rm.S_dn)
        rm.t = self.t - 1
        ev = rm.update(r)
        if ev is None:
            return None
        ev = Event(self.t, ev.direction, self.start, ev.baseline_scale, ev.stat)
        self.events.append(ev)
        self._reset()
        return ev

    @property
    def extra_calibration_steps(self) -> int:
        """Calibration steps spent on discarded windows so far (the time cost of the guard)."""
        return int(sum(d["window_len"] for d in self.discarded))


class GuardedVolatilityMonitor(_GuardMixin, VolatilityMonitor):
    """``VolatilityMonitor(baseline="cs")`` with the retry guard against blind segments. The defaults
    (guard_eta=0, guard_fallback="keep") are E19's criterion (E24's arm G0); the recommended guard is
    guard_eta=0.25, guard_fallback="rank" (docs/theory_blind_guard.md 6.1).

    Extra keyword arguments (all fixed in advance):

    guard_tries : int, default 3
        K, the most calibration windows a segment may try. 1 with ``guard_fallback="keep"`` reproduces the unguarded
        monitor exactly.
    guard_deltas : tuple of K floats or None
        Error budget of each window (U and L each get half). Default ``spending_schedule(delta, K, guard_spending)``.
        Must sum to at most ``delta``.
    guard_spending : "halving" (default) or "equal"
    guard_eta : float, default 0.0
        A window is (near-)blind when ``tau^2 <= (1 + eta) m_up``; 0 is the blind condition itself.
    guard_fallback : "keep" (default) or "rank"
        After K blind windows: monitor with the K-th window, or switch to the awake rank monitor for the segment.
    rank_seed : seed of the rank fallback's tie-breaking draws.

    Guarantee (``docs/theory_blind_guard.md``): if the tolerance (H-up), (H-down) of ``docs/math.md`` Theorem 6.4
    holds for the first m monitored steps relative to the kept window (sufficient, and independent of the selection:
    the tolerance relative to each of the K windows, docs/theory_blind_guard.md (H*)) (and, with the rank fallback,
    the returns after a fallback are exchangeable given the past), P(an alarm within the first m monitored steps of
    a segment) <= delta + m / A; over the first M monitored steps the expected number of alarms is at most
    (delta + M / A) / (1 - delta) when the tolerance holds in every segment up to its alarm (h = infinity).

    ``delta`` holds the current window's budget delta_k during and after a run; the total budget is ``delta_total``.

    Records: ``segments`` (kept windows and rank segments, with ``guard``, ``try_index``, ``delta_window`` and
    ``guard_cal_steps``), ``discarded`` (one dict per discarded window), ``guard_stats``.
    """

    def __init__(self, *args, guard_tries: int = 3, guard_deltas=None, guard_spending: str = "halving",
                 guard_eta: float = 0.0, guard_fallback: str = "keep", rank_seed=None, **kw):
        self.guard_deltas = None
        super().__init__(*args, **kw)
        self._guard_init(guard_tries, guard_deltas, guard_eta, guard_fallback, guard_spending, rank_seed)
        self.period = None
        self.period_resets = 0


class GuardedPeriodicVolatilityMonitor(_GuardMixin, PeriodicVolatilityMonitor):
    """P1 with the retry guard: ``PeriodicVolatilityMonitor(baseline="cs", period=m)`` plus the guard. The defaults
    (guard_eta=0, guard_fallback="keep") are E19's criterion (E24's arm G0); the recommended guard is
    guard_eta=0.25, guard_fallback="rank" (docs/theory_blind_guard.md 6.1).

    The guard arguments are those of ``GuardedVolatilityMonitor``; ``period``, ``offset`` and ``carry`` are those of
    ``PeriodicVolatilityMonitor`` (``statistic`` must be "sr"). In rank-fallback segments the rank statistics restart
    at the same period boundaries.

    Guarantee (Corollary 1 of ``docs/theory_per_period.md`` with Theorem 1 of ``docs/theory_blind_guard.md``): for
    every period (b, b'] of a segment with b' - b <= m, if the tolerance (H-up), (H-down) holds relative to the kept
    window from its end to b' (sufficient, and independent of the selection: the tolerance relative to each of the K
    windows, docs/theory_blind_guard.md (H*)), P(the segment raises an alarm in (b, b']) <= delta + (m + carry) / A,
    marginally (not given the past). In a rank-fallback segment the bound is m / A (its statistics restart at 0
    whatever the carry), given the p-values before b, if the returns after the fallback are exchangeable given the
    past. The bound is per segment and period; for a fixed calendar day see Remark 6 of the theory note; the long-run
    share of alarmed periods is at most (m + carry) / (A (1 - delta)) under the tolerance in every segment (its
    Corollary 2).

    ``delta`` holds the current window's budget delta_k during and after a run; the total budget is ``delta_total``.
    """

    def __init__(self, *args, guard_tries: int = 3, guard_deltas=None, guard_spending: str = "halving",
                 guard_eta: float = 0.0, guard_fallback: str = "keep", rank_seed=None, **kw):
        if kw.get("baseline", "cs") != "cs":
            raise ValueError("the guard needs baseline='cs'")
        if kw.get("statistic", "sr") != "sr":
            raise ValueError("the guarded periodic monitor supports statistic='sr' only")
        kw["baseline"] = "cs"
        self.guard_deltas = None
        super().__init__(*args, **kw)
        self._guard_init(guard_tries, guard_deltas, guard_eta, guard_fallback, guard_spending, rank_seed)
