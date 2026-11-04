"""Tests for experiments/e23_fresh_releases.py on synthetic data only [SYNTHETIC] (and on the committed event file)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import e12_crypto_events as E12  # noqa: E402
import e23_fresh_releases as E23  # noqa: E402

COLS_2024 = list(pd.read_csv(ROOT / "data/events/us_macro_events.csv", nrows=0).columns)


def write_events(path: Path, rows) -> Path:
    """A synthetic event file in the format of data/events/us_macro_events_2026.csv."""
    recs = []
    for date, etype in rows:
        tet = "14:00" if etype == "fomc" else "08:30"
        utc = pd.Timestamp(f"{date} {tet}", tz="America/New_York").tz_convert("UTC").strftime("%Y-%m-%d %H:%M")
        recs.append(dict(date=date, event_type=etype, release_time_et=tet, release_time_utc=utc, status="released",
                         original_date="", description="synthetic", source_url="https://example.org", verification="",
                         access_date="2026-11-03"))
    pd.DataFrame(recs).to_csv(path, index=False)
    return path


# ----------------------------------------------------------------------------
# The committed event file
# ----------------------------------------------------------------------------


def test_event_file_schema():
    ev = pd.read_csv(E23.EVENTS_2026, dtype=str, keep_default_na=False)
    assert list(ev.columns) == COLS_2024 + ["access_date"]
    assert set(ev.event_type) == {"fomc", "cpi"}
    assert set(ev.status) <= {"released", "rescheduled"}
    assert (ev.loc[ev.event_type == "fomc", "release_time_et"] == "14:00").all()
    assert (ev.loc[ev.event_type == "cpi", "release_time_et"] == "08:30").all()
    assert (ev.access_date == "2026-11-03").all()
    assert ev.source_url.str.startswith("https://www.federalreserve.gov/").sum() == (ev.event_type == "fomc").sum()
    assert ev.source_url.str.startswith("https://www.bls.gov/").sum() == (ev.event_type == "cpi").sum()
    assert ev.verification.str.startswith("web search 2026-11-03").all()
    assert (ev.loc[ev.status == "rescheduled", "original_date"] != "").all()
    d = pd.to_datetime(ev.date)
    assert d.is_monotonic_increasing and not d.duplicated().any()
    assert (d.dt.year == 2026).all()
    assert (d[ev.event_type == "fomc"].dt.weekday == 2).all()  # every 2026 statement is on a Wednesday
    assert (ev.event_type == "fomc").sum() == 7 and (ev.event_type == "cpi").sum() == 10
    # E12's reader recomputes every UTC time from the Eastern time and stops on a disagreement
    with E23.e12_period(E23.EVENTS_2026, E23.LOAD_START, E23.MONITOR_START, "2026-11-01"):
        used, raw, excluded = E12.read_events()
    assert len(used) == 17 and len(excluded) == 0
    assert set(used.eastern_zone) == {"EST", "EDT"}


# ----------------------------------------------------------------------------
# The period switch and the inclusion rule
# ----------------------------------------------------------------------------


def test_period_switch_restores_e12():
    before = (E12.EVENTS, E12.START, E12.MONITOR_START, E12.END, E12.read_events)
    with pytest.raises(RuntimeError):
        with E23.e12_period("x.csv", "2025-12-12", "2026-01-01", "2026-02-01"):
            assert E12.END == "2026-02-01" and E12.read_events is E23.read_events
            raise RuntimeError
    assert (E12.EVENTS, E12.START, E12.MONITOR_START, E12.END, E12.read_events) == before
    assert E23.E12_PERIOD["end"] == "2026-01-01" and E23._E12_READ_EVENTS is before[4]


def test_inclusion_rule(tmp_path):
    # data end 2026-02-01 00:00 UTC; a CPI release at 2026-01-31 08:30 EST = 13:30 UTC has its window end at 19:30:
    # included. A window ending exactly at the end is included; one ending a minute later is not.
    path = write_events(tmp_path / "ev.csv", [("2026-01-13", "cpi"), ("2026-01-28", "fomc"), ("2026-01-30", "cpi")])
    with E23.e12_period(path, "2025-12-12", "2026-01-01", "2026-01-30 19:30"):
        used, _, excluded = E12.read_events()
    assert list(used.date.dt.strftime("%Y-%m-%d")) == ["2026-01-13", "2026-01-28", "2026-01-30"]
    assert len(excluded) == 0
    with E23.e12_period(path, "2025-12-12", "2026-01-01", "2026-01-30 19:29"):
        used, _, excluded = E12.read_events()
    assert list(used.date.dt.strftime("%Y-%m-%d")) == ["2026-01-13", "2026-01-28"]
    assert list(excluded.date.dt.strftime("%Y-%m-%d")) == ["2026-01-30"]
    assert "6-hour window" in excluded.reason.iloc[0]
    # a release after the end is dropped by E12's own rule (release < END) before the 6-hour rule
    with E23.e12_period(path, "2025-12-12", "2026-01-01", "2026-01-20"):
        used, _, excluded = E12.read_events()
    assert list(used.date.dt.strftime("%Y-%m-%d")) == ["2026-01-13"]
    assert list(excluded.date.dt.strftime("%Y-%m-%d")) == ["2026-01-28", "2026-01-30"]
    assert excluded.reason.str.startswith("at or after the data end").all()


def test_controls_exclude_event_days_and_follow_the_rules(tmp_path):
    path = write_events(tmp_path / "ev.csv", [("2026-01-13", "cpi"), ("2026-01-28", "fomc"), ("2026-02-13", "cpi"),
                                              ("2026-03-11", "cpi"), ("2026-03-18", "fomc")])
    end = "2026-04-01"
    index = pd.date_range(pd.Timestamp("2026-01-01", tz="UTC"), pd.Timestamp(end, tz="UTC"), freq="1min", inclusive="left")
    with E23.e12_period(path, "2025-12-12", "2026-01-01", end):
        events, _, _ = E12.read_events()
        days = events["date"].to_list()
        for _, r in events.iterrows():
            anchors = E12.control_anchors(index, days, r["date"].weekday(), r["release_time_et"], index.size)
            assert anchors
            t_et = index[anchors].tz_convert("America/New_York")
            assert (t_et.weekday == r["date"].weekday()).all()
            assert (t_et.strftime("%H:%M") == r["release_time_et"]).all()
            d = t_et.tz_localize(None).normalize()
            assert all(min(abs((x - e).days) for e in days) >= 2 for x in d)  # at least 2 days from every event day
            assert min(anchors) >= E12.CAL + E12.PRE and max(anchors) + E12.DOWN_TO <= index.size
    # the Wednesday 14:00 ET controls skip the event Wednesdays (28 January, 11 and 18 March)
    with E23.e12_period(path, "2025-12-12", "2026-01-01", end):
        anchors = E12.control_anchors(index, days, 2, "14:00", index.size)
    wed = set(index[anchors].tz_convert("America/New_York").strftime("%Y-%m-%d"))
    assert {"2026-01-28", "2026-03-11", "2026-03-18"}.isdisjoint(wed) and "2026-01-21" in wed


# ----------------------------------------------------------------------------
# The pooled test
# ----------------------------------------------------------------------------


def synthetic_strata(seed=0):
    """Releases of two types on a few weekdays, with E12-style estats, cstats and a controls function."""
    rng = np.random.default_rng(seed)
    rows = [("2026-01-07", "fomc"), ("2026-01-14", "fomc"), ("2026-02-04", "fomc"), ("2026-01-13", "cpi"),
            ("2026-01-20", "cpi"), ("2026-01-16", "cpi"), ("2026-01-21", "cpi")]
    ev = pd.DataFrame(rows, columns=["date", "event_type"])
    ev["date"] = pd.to_datetime(ev["date"])
    ev["release_time_et"] = np.where(ev.event_type == "fomc", "14:00", "08:30")
    ev["i0"] = np.arange(len(ev)) * 10_000
    ev["weekday"] = ev["date"].dt.day_name()
    groups = {}
    nxt = 500_000
    for _, r in ev.iterrows():
        key = (r["date"].weekday(), r["release_time_et"])
        if key not in groups:
            n = int(rng.integers(5, 30))
            groups[key] = list(range(nxt, nxt + n))
            nxt += n
    keys = ["up_post", "up_post_1h", "up_pre", "down_follow", "any_post", "monitoring_at_start"]
    estats = {int(i): {k: bool(rng.random() < 0.4) for k in keys} for i in ev["i0"]}
    cstats = {i: {k: bool(rng.random() < 0.15) for k in keys} for g in groups.values() for i in g}

    def controls(row):
        return groups[(row["date"].weekday(), row["release_time_et"])]

    return ev, estats, cstats, controls


def windows_from(ev, estats, cstats, controls):
    rows = []
    for _, r in ev.iterrows():
        cs = controls(r)
        rows.append(dict(arm="P1", coin="BTC", detector="x", event_type=r["event_type"], weekday=r["weekday"],
                         n_controls=len(cs), up_post=estats[int(r["i0"])]["up_post"],
                         up_post_1h=estats[int(r["i0"])]["up_post_1h"],
                         control_up_post_rate=float(np.mean([cstats[i]["up_post"] for i in cs])),
                         control_up_post_1h_rate=float(np.mean([cstats[i]["up_post_1h"] for i in cs]))))
    return pd.DataFrame(rows)


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_pooled_function_reproduces_e12_per_type(seed):
    ev, estats, cstats, controls = synthetic_strata(seed)
    ref = pd.DataFrame(E12.stratified_tests("BTC", "x", "streaming", ev, estats, cstats, controls))
    w = windows_from(ev, estats, cstats, controls)
    for etype in ("fomc", "cpi"):
        for outcome in ("up_post", "up_post_1h"):
            mine = E23.stratified_test(E23.strata_from_windows(w[w.event_type == etype], outcome))
            r = ref[(ref.event_type == etype) & (ref.outcome == outcome) & (ref["sample"] == "all")].iloc[0]
            for c in ("events", "events_with", "rate_control_matched", "expected_under_null", "p_stratified_permutation",
                      "p_stratified_permutation_fewer", "p_stratified_fisher"):
                assert mine[c] == r[c], (etype, outcome, c)
    # pooled: E12's test with every release relabelled as one type gives the same numbers
    pooled = E23.stratified_test(E23.strata_from_windows(w))
    ref_all = pd.DataFrame(E12.stratified_tests("BTC", "x", "streaming", ev.assign(event_type="fomc"), estats, cstats,
                                                controls))
    r = ref_all[(ref_all.event_type == "fomc") & (ref_all.outcome == "up_post") & (ref_all["sample"] == "all")].iloc[0]
    assert pooled["events"] == len(ev)
    for c in ("events_with", "expected_under_null", "p_stratified_permutation", "p_stratified_fisher"):
        assert pooled[c] == r[c], c
    chk = E23.check_pooled_against_e12(pd.DataFrame(E12.stratified_tests("BTC", "x", "streaming", ev, estats, cstats,
                                                                         controls)).assign(arm="P1"), w)
    assert (chk.mismatches == 0).all() and len(chk) == 4


def test_pooled_p_value_by_hand():
    # two releases, each alone in its group: y = (1, 0); (n, k) = (9, 0) and (4, 1): p_i = 1/10 and 1/5; T = 1
    w = pd.DataFrame(dict(event_type=["fomc", "cpi"], weekday=["Wednesday", "Tuesday"], n_controls=[9, 4],
                          up_post=[True, False], control_up_post_rate=[0.0, 0.25]))
    t = E23.stratified_test(E23.strata_from_windows(w))
    assert t["events_with"] == 1 and np.isclose(t["expected_under_null"], 0.3)
    assert np.isclose(t["p_stratified_permutation"], 1 - 0.9 * 0.8)
    assert np.isclose(t["rate_control_matched"], 0.125)
    bad = w.assign(control_up_post_rate=[0.0, 0.3])  # 0.3 * 4 is not an integer
    with pytest.raises(ValueError):
        E23.strata_from_windows(bad)


# ----------------------------------------------------------------------------
# Power on toy cases
# ----------------------------------------------------------------------------


def toy_calendar(groups):
    rows = []
    for g, (etype, n, m) in enumerate(groups):
        rows += [dict(event_type=etype, group=(g, "x"), n=n)] * m
    return pd.DataFrame(rows)


def test_power_toy_cases():
    cal = toy_calendar([("fomc", 10, 3)])
    sure = E23.simulate_power(cal, {"fomc": (1.0, 0.0)}, reps=200, seed=1, designs=("fomc",))
    assert sure["fomc"]["power"] == 1.0 and sure["fomc"]["mean_T"] == 3
    never = E23.simulate_power(cal, {"fomc": (0.0, 0.0)}, reps=200, seed=1, designs=("fomc",))
    assert never["fomc"]["power"] == 0.0
    # one release, 19 controls that never have the outcome: y = 1 gives p = 1/20 = alpha, so power = P(y = 1) = 0.3
    one = toy_calendar([("fomc", 19, 1)])
    res = E23.simulate_power(one, {"fomc": (0.3, 0.0)}, reps=20_000, seed=2, designs=("fomc",))["fomc"]
    assert abs(res["power"] - 0.3) < 4 * res["mc_se"]
    # pooled over two types equals the type-specific design when only one type is present
    both = E23.simulate_power(cal, {"fomc": (0.5, 0.1)}, reps=500, seed=3, designs=("fomc", E23.POOLED))
    assert both["fomc"]["power"] == both[E23.POOLED]["power"]
    # common random numbers and a fixed seed: the result is reproducible
    a = E23.simulate_power(toy_calendar([("fomc", 29, 2), ("cpi", 40, 3)]), {"fomc": (0.5, 0.1), "cpi": (0.2, 0.06)},
                           reps=300, seed=4)
    b = E23.simulate_power(toy_calendar([("fomc", 29, 2), ("cpi", 40, 3)]), {"fomc": (0.5, 0.1), "cpi": (0.2, 0.06)},
                           reps=300, seed=4)
    assert a == b


def test_power_increases_with_effect():
    cal = toy_calendar([("fomc", 29, 7), ("cpi", 29, 6), ("cpi", 40, 4)])
    lo = E23.simulate_power(cal, {"fomc": (0.2, 0.1), "cpi": (0.08, 0.06)}, reps=2000, seed=5)
    hi = E23.simulate_power(cal, {"fomc": (0.6, 0.1), "cpi": (0.3, 0.06)}, reps=2000, seed=5)
    for d in E23.DESIGNS:
        assert hi[d]["power"] > lo[d]["power"]


# ----------------------------------------------------------------------------
# Families, data end
# ----------------------------------------------------------------------------


def test_families():
    rows = []
    for arm in E23.ARMS:
        for coin in E23.COINS:
            for et in ("fomc", "cpi", E23.POOLED):
                for outcome, sample in (("up_post", "all"), ("up_post_1h", "all"), ("up_post", "watched at start")):
                    rows.append(dict(arm=arm, coin=coin, event_type=et, outcome=outcome, sample=sample,
                                     p_stratified_permutation=0.01))
    t = E23.add_families(pd.DataFrame(rows))
    sizes = t[t.family != ""].groupby("family").size().to_dict()
    assert sizes == {"F23a": 1, "F23b": 2, "F23c": 6, "F23d": 6, "F23e": 6, "F23f": 6}
    prim = t[t.primary_test]
    assert len(prim) == 1 and prim.iloc[0][["arm", "coin", "event_type"]].tolist() == ["P1", "BTC", E23.POOLED]
    assert prim.iloc[0]["p_holm_family"] == 0.01 and bool(prim.iloc[0]["primary_passes"])
    assert np.allclose(t.loc[t.family == "F23c", "p_holm_family"], 0.06)


def test_data_end(tmp_path):
    paths = {}
    for coin, last in (("BTC", "2026-10-31 23:59"), ("ETH", "2026-10-31 23:58"), ("SOL", "2026-10-31 23:59")):
        ts = pd.date_range("2026-10-31 23:50", last, freq="1min", tz="UTC")
        p = tmp_path / f"{coin}.csv"
        ms = (ts - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta(milliseconds=1)
        pd.DataFrame(dict(timestamp=ms)).to_csv(p, index=False)
        paths[coin] = [p]
    assert E23.data_end(paths) == "2026-10-31 23:59"
