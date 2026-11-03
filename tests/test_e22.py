"""Tests for experiments/e22_guaranteed_events.py on synthetic data only [SYNTHETIC]."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import e22_guaranteed_events as E22  # noqa: E402


def minute_index(start: str, days: int) -> pd.DatetimeIndex:
    return pd.date_range(pd.Timestamp(start, tz="UTC"), periods=days * 1440, freq="1min")


def test_utc_day_offset():
    assert E22.utc_day_offset(minute_index("2024-01-21 00:00", 2)) == 0
    assert E22.utc_day_offset(minute_index("2024-01-21 05:00", 2)) == 1140
    idx = minute_index("2024-01-21 00:00", 2)
    with pytest.raises(ValueError):
        E22.utc_day_offset(idx.delete(100))  # a missing minute shifts every later boundary


def test_arms_are_the_registered_configurations():
    p1, s1, p0 = E22.make_arm("P1", 7), E22.make_arm("S1"), E22.make_arm("P0")
    assert (p1.A, p1.baseline, p1.delta, p1.period, p1.offset, p1.band, p1.cal_len) == (36000.0, "cs", 0.01, 1440, 7, 0.25, 1440)
    assert (s1.A, s1.baseline, s1.delta, s1.band) == (36000.0, "cs", 0.01, 0.25)
    assert not hasattr(s1, "period")
    assert (p0.A, p0.baseline, p0.band) == (28800.0, "plugin", 0.25)
    s2 = E22.make_arm("S2")
    assert (s2.A, s2.band, s2.restart_len, s2.cal_len) == (28800.0, 0.25, 0, 120)
    # every per-day budget is 5%: delta + m / A for the cs arms, m / A for the rank and plug-in arms
    assert np.isclose(E22.DELTA + 1440 / p1.A, 0.05) and np.isclose(1440 / s2.A, 0.05)


def test_p1_restarts_at_every_utc_midnight_after_calibration():
    idx = minute_index("2024-01-21 05:00", 6)
    off = E22.utc_day_offset(idx)
    u = np.random.default_rng(1).standard_normal(idx.size)
    mon = E22.make_arm("P1", off)
    monitoring = E22.run_monitor(mon, u)
    assert not mon.events  # i.i.d. Gaussian, six days: no alarm with this seed
    start = mon.segments[0]["start"]
    midnights = np.flatnonzero((idx.hour == 0) & (idx.minute == 0))
    assert mon.period_resets == int(np.sum(midnights > start))
    assert monitoring[start:].all() and not monitoring[:start].any()
    # without the restart (S1) the statistics are the same up to the first boundary after calibration
    s1 = E22.make_arm("S1")
    E22.run_monitor(s1, u[: midnights[midnights > start][0]])
    p1 = E22.make_arm("P1", off)
    E22.run_monitor(p1, u[: midnights[midnights > start][0]])
    assert np.array_equal(s1.S_up, p1.S_up) and np.array_equal(s1.S_dn, p1.S_dn)


def test_blind_segment_is_not_watched_and_cannot_raise_an_up_alarm():
    rng = np.random.default_rng(2)
    n = 4 * 1440
    u = rng.standard_normal(n)
    u[:360] *= 0.05  # the quarter of the calibration window that sets tau is very quiet: tau^2 <= m_up
    u[2000:2600] *= 4  # a burst that a watching monitor would flag
    for arm in ("P1", "S1"):
        mon = E22.make_arm(arm, 0)
        monitoring = E22.run_monitor(mon, u)
        segs, seg_id = E22.segment_table(mon, n, x_based=True)
        assert np.array_equal(monitoring, seg_id >= 0)
        assert segs.loc[0, "tau2_over_m_up"] <= 1.0 and not segs.loc[0, "up_active"]
        watched = E22.watched_flags(monitoring, seg_id, segs)
        first = slice(segs.loc[0, "start"], segs.loc[0, "end"] + 1)
        assert monitoring[first].all() and not watched[first].any()
        assert all(e.direction != "up" for e in mon.events if e.monitored_since == segs.loc[0, "start"])


def test_segment_table_rank_monitor_is_never_blind():
    rng = np.random.default_rng(3)
    u = rng.standard_normal(3 * 1440)
    u[1500:1800] *= 5
    mon = E22.make_arm("S2")
    monitoring = E22.run_monitor(mon, u)
    segs, seg_id = E22.segment_table(mon, u.size, x_based=False)
    assert len(mon.events) >= 1
    assert np.array_equal(monitoring, seg_id >= 0)
    assert segs["up_active"].all()
    assert np.array_equal(E22.watched_flags(monitoring, seg_id, segs), monitoring)
    assert not monitoring[:120].any() and monitoring[120:].all()  # awake after the first 120 returns, never blind
    # consecutive segments: a new one starts at the step after each alarm
    for k in range(1, len(segs)):
        assert segs.loc[k, "start"] == segs.loc[k - 1, "end"] + 1


def test_period_counts_and_first_day():
    n = 5 * 1440
    monitoring = np.zeros(n, bool)
    seg_id = np.full(n, -1, np.int32)
    monitoring[1440:3000] = True  # segment 0: day 1 and part of day 2, alarm at 2999
    seg_id[1440:3000] = 0
    monitoring[4440:n] = True  # segment 1: from day 3 to the end, no alarm
    seg_id[4440:n] = 1
    pc = E22.period_counts(monitoring, seg_id, 0, [2999])
    assert pc == dict(periods=4, periods_with_alarm=1)  # segment 0 on days 1 and 2, segment 1 on days 3 and 4
    segs = pd.DataFrame(dict(start=[1440, 4440]))

    class Ev:
        def __init__(self, i, s):
            self.index, self.monitored_since = i, s

    fd = E22.first_day(segs, [Ev(2999, 1440)], n)
    assert fd == dict(segments_with_a_day=2, first_day_alarms=0)  # 2999 - 1440 = 1559 >= 1440
    fd = E22.first_day(segs, [Ev(2000, 1440)], n)
    assert fd == dict(segments_with_a_day=2, first_day_alarms=1)
    fd = E22.first_day(pd.DataFrame(dict(start=[n - 100])), [], n)
    assert fd == dict(segments_with_a_day=0, first_day_alarms=0)


def test_long_run_bounds():
    M = 400 * 1440
    b = E22.long_run_bound_per_monitored_day("P1", M)
    assert np.isclose(b, (0.01 + M / 36000) / 0.99 / 400)
    assert 0.0404 < b < 0.0405
    assert E22.long_run_bound_per_monitored_day("S1", M) == b
    assert np.isclose(E22.long_run_bound_per_monitored_day("S2", M), 0.05)
    assert np.isclose(E22.long_run_bound_per_monitored_day("P0", M), 0.05)


def test_families_and_holm():
    rows = []
    rng = np.random.default_rng(4)
    for arm in E22.ARMS:
        for coin in ("BTC", "ETH", "SOL"):
            for etype in ("fomc", "cpi"):
                for outcome, sample in (("up_post", "all"), ("up_post", "monitoring at start"), ("up_post_1h", "all")):
                    rows.append(dict(arm=arm, coin=coin, event_type=etype, outcome=outcome, sample=sample,
                                     p_stratified_permutation=float(rng.uniform(0.001, 0.2))))
    t = E22.add_families(pd.DataFrame(rows))
    sizes = t[t.family != ""].groupby("family").size().to_dict()
    assert sizes == {"F22a": 2, "F22b": 4, "F22c": 6, "F22d": 6, "F22e": 2}
    a = t[t.family == "F22a"].set_index("event_type")
    p = a["p_stratified_permutation"]
    lo, hi = p.idxmin(), p.idxmax()
    assert np.isclose(a.loc[lo, "p_holm_family"], min(1, 2 * p[lo]))
    assert np.isclose(a.loc[hi, "p_holm_family"], max(min(1, 2 * p[lo]), p[hi]))
    assert t["primary_test"].sum() == 1
    prim = t[t["primary_test"]].iloc[0]
    assert (prim["arm"], prim["coin"], prim["event_type"]) == ("P1", "BTC", "fomc")
    assert (t[t.family == ""]["status"] != "pre-registered").all()


def test_compare_detects_any_difference():
    a = pd.DataFrame(dict(coin=["BTC", "ETH"], x=[0.1, np.nan], y=[True, False]))
    ok = pd.DataFrame(E22._compare(a, a.copy(), ["coin"], ["x", "y"], "same"))
    assert (ok.mismatches == 0).all()
    b = a.copy()
    b.loc[0, "x"] = 0.1 + 1e-15
    bad = pd.DataFrame(E22._compare(a, b, ["coin"], ["x", "y"], "diff"))
    assert bad.set_index("column").loc["x", "mismatches"] == 1
    c = a.iloc[:1]
    keys = pd.DataFrame(E22._compare(a, c, ["coin"], ["x"], "rows"))
    assert keys.set_index("column").loc["(row keys)", "mismatches"] == 1
