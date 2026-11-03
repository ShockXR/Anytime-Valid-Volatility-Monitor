"""Tests of the retry guard against blind segments (avmonitor.guard; docs/theory_blind_guard.md). [SYNTHETIC]"""

import math

import numpy as np
import pytest

from avmonitor.confseq import cs_baseline, eb_confseq
from avmonitor.guard import (
    GuardedPeriodicVolatilityMonitor,
    GuardedVolatilityMonitor,
    near_blind_delay_floor,
    spending_schedule,
)
from avmonitor.monitor import VolatilityMonitor
from avmonitor.periodic import PeriodicVolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE

CAL, NT = 200, 50  # small calibration window: tau from the first 50 returns, the sequence on the next 150
SMALL = dict(c=3.0, cal_len=CAL, null_band=0.25, baseline="cs")


def blind_window(rng, big=10.0):
    """tau part: +-1 (so tau = 3 * 1.4826 exactly); sequence part: +-big, every X equal to tau^2: blind."""
    return np.r_[rng.choice([-1.0, 1.0], NT), rng.choice([-big, big], CAL - NT)]


def calm_window(rng):
    return np.r_[rng.choice([-1.0, 1.0], NT), rng.standard_normal(CAL - NT)]


# ----------------------------------------------------------------------------
# The blind condition
# ----------------------------------------------------------------------------


def test_blind_condition_on_constructed_data():
    rng = np.random.default_rng(0)
    u = np.r_[blind_window(rng), rng.choice([-50.0, 50.0], 600)]  # then huge returns for six periods
    mon = PeriodicVolatilityMonitor(threshold=1000.0, delta=0.05, period=100, **SMALL)
    mon.run(u)
    seg = mon.segments[0]
    tau = 3.0 * MAD_NORMAL_SCALE
    assert seg["tau2"] == pytest.approx(tau * tau)
    assert seg["tau2"] / seg["m_up"] <= 1.0  # blind (Lemma 1)
    assert not any(e.direction == "up" for e in mon.events)  # huge returns, yet no increase alarm can come
    # Lemma 1: the increase statistic grows by at most 1 per step after a restart
    m = PeriodicVolatilityMonitor(threshold=1e9, delta=0.05, period=10**6, **SMALL)
    for r in u[:CAL]:
        m.update(r)
    for k in range(1, 300):
        m.update(50.0)
        assert m.S_up.mean() <= k + 1e-9


def test_guard_discards_a_blind_window_and_keeps_the_next():
    rng = np.random.default_rng(1)
    u = np.r_[blind_window(rng), rng.standard_normal(CAL) * 10.0, rng.standard_normal(300) * 10.0]
    g = GuardedVolatilityMonitor(threshold=1000.0, delta=0.05, guard_tries=2, **SMALL)
    g.run(u)
    assert len(g.discarded) == 1 and g.discarded[0]["tau2_over_m_up"] <= 1.0
    assert g.discarded[0]["window_start"] == 0 and g.discarded[0]["window_len"] == CAL
    seg = g.segments[0]
    assert seg["guard"] == "kept" and seg["try_index"] == 2 and seg["start"] == 2 * CAL
    assert seg["tau2"] / seg["m_up"] > 1.0
    assert seg["guard_cal_steps"] == 2 * CAL and g.extra_calibration_steps == CAL


def test_near_blind_eta():
    """A window with ratio rho is discarded exactly when rho <= 1 + eta (eta fixed in advance)."""
    rng = np.random.default_rng(2)
    tau = 3.0 * MAD_NORMAL_SCALE
    z = np.where(rng.random(CAL - NT) < 0.55, tau * 1.5, 0.0) * rng.choice([-1.0, 1.0], CAL - NT)
    u = np.r_[rng.choice([-1.0, 1.0], NT), z, rng.standard_normal(2 * CAL)]
    ref = VolatilityMonitor(threshold=1000.0, delta=0.025, **SMALL)  # = window 1 of the guard below (delta_1)
    ref.run(u[:CAL])
    rho = ref.segments[0]["tau2"] / ref.segments[0]["m_up"]
    assert 1.0 < rho < 2.0
    for eta, discarded in ((rho - 1 + 0.01, 1), (rho - 1 - 0.01, 0)):
        g = GuardedVolatilityMonitor(threshold=1000.0, delta=0.05, guard_tries=2, guard_eta=eta, **SMALL)
        g.run(u)
        assert len(g.discarded) == discarded


def test_near_blind_delay_floor():
    """Proposition 2(b): exactly 2A steps at rho = 1, attained by X = tau^2, and decreasing in rho."""
    assert near_blind_delay_floor(1.0, 50.0) == 100.0
    assert math.isinf(near_blind_delay_floor(0.99, 50.0))
    A = 36_000.0
    d = [near_blind_delay_floor(r, A) for r in (1.01, 1.05, 1.25, 2.0)]
    assert d == sorted(d, reverse=True) and d[0] == 1503.0 and d[-1] == 30.0
    # attained: a monitor whose window has rho = tau^2 / m_up fed X = tau^2 alarms upward exactly then
    rng = np.random.default_rng(3)
    tau = 3.0 * MAD_NORMAL_SCALE
    u = np.r_[rng.choice([-1.0, 1.0], NT), np.where(rng.random(CAL - NT) < 0.4, tau, 0.0)]
    mon = VolatilityMonitor(threshold=200.0, delta=0.05, **SMALL)
    mon.run(u)
    seg = mon.segments[0]
    rho = seg["tau2"] / seg["m_up"]
    floor = near_blind_delay_floor(rho, 200.0)
    mon.run(np.full(int(floor) + 5, 2 * tau))
    up = [e for e in mon.events if e.direction == "up"]
    assert up and up[0].index - seg["start"] + 1 == floor


# ----------------------------------------------------------------------------
# Spending and bookkeeping
# ----------------------------------------------------------------------------


def test_spending_schedule():
    assert spending_schedule(0.01, 1) == (0.01,)
    assert spending_schedule(0.01, 3) == pytest.approx((0.005, 0.0025, 0.0025))
    assert sum(spending_schedule(0.01, 5)) == pytest.approx(0.01)
    assert spending_schedule(0.01, 4, "equal") == pytest.approx((0.0025,) * 4)
    with pytest.raises(ValueError):
        spending_schedule(0.01, 0)
    with pytest.raises(ValueError):
        GuardedVolatilityMonitor(threshold=100.0, delta=0.01, guard_tries=2, guard_deltas=(0.008, 0.004), **SMALL)
    with pytest.raises(ValueError):
        GuardedVolatilityMonitor(threshold=100.0, delta=0.01, guard_tries=3, guard_deltas=(0.005, 0.005), **SMALL)
    with pytest.raises(ValueError):
        GuardedVolatilityMonitor(threshold=100.0, delta=0.01, c=3.0, cal_len=CAL)  # plug-in baseline


def test_retry_bookkeeping_and_budgets():
    """Each window's bounds are exactly a fresh sequence at level delta_k / 2 on that window alone."""
    rng = np.random.default_rng(4)
    w1, w2 = blind_window(rng), blind_window(rng, big=20.0)
    w3 = np.r_[rng.choice([-1.0, 1.0], NT), rng.standard_normal(CAL - NT) * 1.2]
    u = np.r_[w1, w2, w3, rng.standard_normal(100) * 1.2]
    deltas = (0.02, 0.01, 0.02)
    g = GuardedVolatilityMonitor(threshold=1e6, delta=0.05, guard_tries=3, guard_deltas=deltas, **SMALL)
    g.run(u)
    assert [d["try_index"] for d in g.discarded] == [1, 2]
    assert [d["delta_window"] for d in g.discarded] == [0.02, 0.01]
    seg = g.segments[0]
    assert seg["try_index"] == 3 and seg["delta_window"] == 0.02 and seg["start"] == 3 * CAL
    ref = cs_baseline(w3[None, :], c=3.0, delta_upper=0.01, delta_lower=0.01, tau_frac=NT / CAL)
    assert seg["upper"] == pytest.approx(float(ref.upper[0]), rel=1e-12)
    assert seg["lower"] == pytest.approx(float(ref.lower[0]), rel=1e-12)
    assert g.guard_stats == {"windows_tried": 3, "windows_discarded": 2, "fallback_keep": 0, "fallback_rank": 0}


def test_budget_resets_after_an_alarm_and_fallbacks():
    rng = np.random.default_rng(5)
    # segment 1: window 1 blind, window 2 kept; then a big rise -> up alarm; segment 2 starts again at try 1
    u = np.r_[blind_window(rng), rng.standard_normal(CAL), rng.standard_normal(200) * 6.0,
              rng.standard_normal(CAL) * 6.0, rng.standard_normal(200) * 6.0]
    g = GuardedVolatilityMonitor(threshold=50.0, delta=0.05, guard_tries=2, **SMALL)
    g.run(u)
    assert g.events and g.events[0].direction == "up"
    assert g.segments[0]["try_index"] == 2 and g.segments[1]["try_index"] == 1
    assert g.segments[1]["delta_window"] == g.guard_deltas[0]
    # fallback "keep": all K windows blind, the K-th is monitored as the unguarded monitor would
    v = np.r_[blind_window(rng), blind_window(rng), rng.standard_normal(100) * 10.0]
    k = GuardedVolatilityMonitor(threshold=1e6, delta=0.05, guard_tries=2, **SMALL)
    k.run(v)
    assert k.segments[0]["guard"] == "blind kept" and k.guard_stats["fallback_keep"] == 1
    # fallback "rank": the awake rank monitor takes over at once and is never blind
    r = GuardedPeriodicVolatilityMonitor(threshold=50.0, delta=0.05, guard_tries=2, guard_fallback="rank",
                                         rank_seed=1, period=10**6, **SMALL)
    w = np.r_[blind_window(rng), blind_window(rng), rng.standard_normal(300), rng.standard_normal(300) * 8.0]
    flags = []
    for x in w:
        flags.append(r.calibrated)
        ev = r.update(x)
        if ev is not None:
            break
    assert r.segments[0]["guard"] == "rank fallback" and r.segments[0]["start"] == 2 * CAL
    assert flags[2 * CAL] and not flags[2 * CAL - 1]
    assert ev is not None and ev.direction == "up" and ev.index >= 2 * CAL + 300 and ev.monitored_since == 2 * CAL
    assert r.guard_stats["fallback_rank"] == 1
    # after its alarm the next segment calibrates again with window 1 and the budget delta_1
    assert not r.calibrated and r._rank is None and r._try == 1 and r.delta == r.guard_deltas[0]


# ----------------------------------------------------------------------------
# Identity with the unguarded monitor
# ----------------------------------------------------------------------------


def changing_stream(rng, n, law="gauss"):
    z = rng.standard_normal(n) if law == "gauss" else rng.standard_t(3, n)
    s = np.ones(n)
    for a in range(300, n, 450):  # volatility moves now and then (alarms, several segments), never 6x in a window
        s[a:] = s[a - 1] * rng.choice([0.6, 1.0, 1.7])
    return z * s


@pytest.mark.parametrize("law", ["gauss", "t3"])
def test_identity_with_one_try(law):
    rng = np.random.default_rng(6)
    for _ in range(4):
        u = changing_stream(rng, 3000, law)
        for cls, gcls, kw in ((VolatilityMonitor, GuardedVolatilityMonitor, {}),
                              (PeriodicVolatilityMonitor, GuardedPeriodicVolatilityMonitor, dict(period=100))):
            a = cls(threshold=300.0, delta=0.05, **SMALL, **kw)
            b = gcls(threshold=300.0, delta=0.05, guard_tries=1, **SMALL, **kw)
            a.run(u)
            b.run(u)
            assert [(e.index, e.direction, e.monitored_since) for e in a.events] == \
                   [(e.index, e.direction, e.monitored_since) for e in b.events]
            assert [s["upper"] for s in a.segments] == [s["upper"] for s in b.segments]


def test_identity_when_no_window_is_blind():
    """K = 3 with no blind window equals the unguarded monitor run at delta_1."""
    rng = np.random.default_rng(7)
    n_alarms = 0
    for _ in range(4):
        u = changing_stream(rng, 3000)
        a = PeriodicVolatilityMonitor(threshold=300.0, delta=0.025, period=100, **SMALL)
        b = GuardedPeriodicVolatilityMonitor(threshold=300.0, delta=0.05, guard_tries=3, period=100, **SMALL)
        a.run(u)
        b.run(u)
        assert b.guard_stats["windows_discarded"] == 0
        assert [(e.index, e.direction) for e in a.events] == [(e.index, e.direction) for e in b.events]
        assert [s["lower"] for s in a.segments] == [s["lower"] for s in b.segments]
        n_alarms += len(a.events)
    assert n_alarms >= 4  # several segments were compared


# ----------------------------------------------------------------------------
# Proposition 3: a kept window that is not blind has a failing upper bound
# ----------------------------------------------------------------------------


def test_selection_construction():
    """Proposition 3: under the construction a window that is kept and not blind always has a failing upper bound,
    and with the same budget for every window the kept window fails about K times as often as one window."""
    rng = np.random.default_rng(8)
    tau = 3.0 * MAD_NORMAL_SCALE
    q = 1.0 / 1.25  # eta = 0, b = 0.25
    P, K, d = 600, 3, 0.3
    kept = 0
    for _ in range(P):
        parts = []
        for _k in range(K):
            parts += [rng.choice([-1.0, 1.0], NT), np.where(rng.random(CAL - NT) < q, tau, 0.0)]
        g = GuardedVolatilityMonitor(threshold=1e6, delta=0.95, guard_tries=K, guard_deltas=(d,) * K, **SMALL)
        g.run(np.concatenate(parts))
        seg = g.segments[0]
        if seg["guard"] == "kept":
            kept += 1
            assert seg["upper"] < q * tau * tau  # Proposition 3(b): kept and not blind => the bound fails
        else:
            assert seg["upper"] >= q * tau * tau * (1 - 1e-12)
    z = (rng.random((40000, CAL - NT)) < q).astype(float)
    p1 = float(np.mean(eb_confseq(z, 1.0, d / 2, d / 2, n_plan=CAL - NT).upper[:, -1] < q))  # one window
    expected = 1 - (1 - p1) ** K  # Proposition 3(c)
    rate = kept / P
    assert abs(rate - expected) < 3.5 * math.sqrt(expected * (1 - expected) / P)
    assert rate > 1.8 * p1


# ----------------------------------------------------------------------------
# Monte Carlo null checks of the proved bounds
# ----------------------------------------------------------------------------


def null_paths(rng, law, n_paths, n):
    z = rng.standard_normal((n_paths, n)) if law == "gauss" else rng.standard_t(3, (n_paths, n))
    s = np.where(np.arange(n) < NT, 1.0, 8.0)  # a step inside the first tau part: blind first windows,
    return z * s  # and a constant law from the end of that tau part on (Corollary 6.5)


@pytest.mark.parametrize("law", ["gauss", "t3"])
def test_null_bounds_first_m_and_per_period(law):
    rng = np.random.default_rng(9 if law == "gauss" else 10)
    delta, m, A = 0.05, 100, 100 / 0.15  # bound delta + m / A = 0.2
    P, n_per = 260, 3
    n = 3 * CAL + (n_per + 1) * m
    paths = null_paths(rng, law, P, n)
    first_m = per = 0
    per_counts = np.zeros(n_per)
    retried = 0
    for u in paths:
        g = GuardedPeriodicVolatilityMonitor(threshold=A, delta=delta, guard_tries=3, period=m, **SMALL)
        g.run(u)
        retried += g.guard_stats["windows_discarded"] > 0
        if not g.segments:
            continue
        s = g.segments[0]["start"]
        e = [ev for ev in g.events if ev.monitored_since == s]
        a = e[0].index if e else None
        first_m += a is not None and a - s < m
        b0 = s + (-s) % m  # first boundary at or after s; periods (b0 - m, b0], (b0, b0 + m], ...
        bounds = [s] + [b0 + i * m for i in range(n_per)]
        for i in range(n_per):
            lo, hi = bounds[i], bounds[i + 1]
            if hi <= n and a is not None and lo <= a < hi:
                per_counts[i] += 1
        per += 1
    bound = delta + m / A
    se = math.sqrt(bound * (1 - bound) / P)
    assert retried > P / 2  # the guard was exercised
    assert first_m / P <= bound + 3 * se
    assert np.all(per_counts / max(per, 1) <= bound + 3 * se)


def test_null_bound_rank_fallback():
    """Theorem 3 on the rank branch: K = 1 and the rank fallback, i.i.d. returns after the fallback."""
    rng = np.random.default_rng(11)
    m, A, P = 100, 500.0, 250
    hits = falls = 0
    for u in null_paths(rng, "gauss", P, CAL + 2 * m):
        g = GuardedVolatilityMonitor(threshold=A, delta=0.05, guard_tries=1, guard_fallback="rank", rank_seed=rng, **SMALL)
        g.run(u)
        if g.segments and g.segments[0]["guard"] == "rank fallback":
            falls += 1
            s = g.segments[0]["start"]
            hits += any(ev.monitored_since == s and ev.index - s < m for ev in g.events)
    bound = 0.05 + m / A
    assert falls > P / 2
    assert hits / P <= bound + 3 * math.sqrt(bound * (1 - bound) / P)


# ----------------------------------------------------------------------------
# The forward pre-registration's guards (experiments/e24_blind_guard.py; no market data)
# ----------------------------------------------------------------------------


def _e24():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
    import e24_blind_guard as E24

    return E24


def test_future_part_refuses_early_and_with_unverified_releases():
    E24 = _e24()
    assert E24.future_refusals(None)  # no end date fixed
    import pandas as pd

    why = E24.future_refusals("2027-07-01")
    if pd.Timestamp.now(tz="UTC") < pd.Timestamp("2027-07-01", tz="UTC"):
        assert any("before the end date" in w for w in why)
    assert any("status 'scheduled'" in w for w in why)  # every 2027 row of the event file is still scheduled
    assert any("before the end date 2027-07-01" in w for w in E24.future_refusals("2027-07-01", data_end="2027-06-30 23:59"))


def test_event_file_and_planning_calendar():
    import pandas as pd

    E24 = _e24()
    ev = pd.read_csv(E24.FUTURE_EVENTS, dtype=str, keep_default_na=False)
    assert list(ev.columns) == ["date", "event_type", "release_time_et", "release_time_utc", "status", "original_date",
                                "description", "source_url", "verification", "access_date"]
    assert ev.date.is_monotonic_increasing
    assert set(ev.status) <= {"released", "rescheduled", "cancelled", "no_statement", "scheduled"}
    et = pd.to_datetime(ev.date + " " + ev.release_time_et).dt.tz_localize("America/New_York").dt.tz_convert("UTC")
    assert (et.dt.strftime("%Y-%m-%d %H:%M") == ev.release_time_utc).all()  # UTC column agrees with the time zone rules
    assert ((ev.event_type == "fomc") <= (ev.release_time_et == "14:00")).all()
    assert ((ev.event_type == "cpi") <= (ev.release_time_et == "08:30")).all()
    assert ev.verification.str.len().min() > 50 and (ev.access_date == "2026-11-03").all()
    # every month from November 2026 to December 2027 has exactly one CPI row
    m = pd.to_datetime(ev.loc[ev.event_type == "cpi", "date"]).dt.to_period("M")
    assert list(m) == list(pd.period_range("2026-11", "2027-12", freq="M"))
    # the plan: the no-statement meeting dropped, no placeholder before 2028, every 6-hour window inside the period
    pe = E24.planning_events("2027-12-01")
    assert (pe.kind == "event file").all() and "2027-01-27" not in set(pe.date.dt.strftime("%Y-%m-%d"))
    assert pe.release_utc.min() >= pd.Timestamp(E24.FUTURE_MONITOR_START, tz="UTC")
    assert (pe.release_utc + pd.Timedelta(hours=6) <= pd.Timestamp("2027-12-01", tz="UTC")).all()
    pe28 = E24.planning_events("2029-01-01")
    ph = pe28[pe28.kind != "event file"]
    assert (ph.date.dt.year == 2028).all() and (ph.event_type == "fomc").sum() == 6 and (ph.event_type == "cpi").sum() == 12
