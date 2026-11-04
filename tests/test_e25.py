"""Tests for experiments/e25_sequential_events.py [SYNTHETIC]: the e-value and the e-process, the strata of the
pre-registration (section 4), completion, the status part's refusals, its event-file checks and ledger (review round
11, M1) and the missing minutes per window (m10). No market data are read; two tests read the committed forward event
file only."""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import e12_crypto_events as E12  # noqa: E402
import e25_sequential_events as E25  # noqa: E402


def rel(rows):
    """Synthetic releases: (date, type) with the type's Eastern clock time."""
    return pd.DataFrame([dict(date=d, event_type=t, release_time_et="14:00" if t == "fomc" else "08:30") for d, t in rows])


def events(rows):
    """A synthetic event file: (date, type, status), each row verified on its own date (access date = date)."""
    return pd.DataFrame([dict(date=d, event_type=t, release_time_et="14:00" if t == "fomc" else "08:30", status=s,
                              original_date="", access_date=d) for d, t, s in rows],
                        columns=["date", "event_type", "release_time_et", "status", "original_date", "access_date"])


# ----------------------------------------------------------------------------
# The e-value (Lemma 1) and the e-process
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("psi", E25.PSI + (1.5, 100.0))
def test_lemma1_identity(psi):
    for n in range(0, 9):
        for s in range(0, n + 2):
            p = s / (n + 1)
            mean = p * E25.e_value(1, s, n, psi) + (1 - p) * E25.e_value(0, s, n, psi) if 0 < s < n + 1 else (
                E25.e_value(1 if s else 0, s, n, psi))
            assert abs(float(mean) - 1.0) < 1e-12
            q, pf = Fraction(psi).limit_denominator(10**6), Fraction(s, n + 1)
            exact = pf * q / (1 - pf + q * pf) + (1 - pf) / (1 - pf + q * pf)
            assert exact == 1


def test_identity_table_gate1():
    t = E25.identity_table()
    assert len(t) == sum(n + 2 for n in range(9)) * len(E25.PSI)
    assert t.mean_exact_is_1.all() and t.abs_deviation_float.max() <= 1e-12


@pytest.mark.parametrize("psi", E25.PSI)
def test_e_is_one_when_uninformative(psi):
    for n in range(0, 9):
        assert E25.e_value(0, 0, n, psi) == 1.0  # no window has the outcome
        assert E25.e_value(1, n + 1, n, psi) == 1.0  # every window has it
    assert E25.e_value(0, 0, 0, psi) == 1.0 and E25.e_value(1, 1, 0, psi) == 1.0  # no control


def test_e_process_hand_computed():
    # strata (Z, S, n): (1, 1, 3), (0, 1, 2), (1, 2, 4)
    z, s, n = np.array([1, 0, 1]), np.array([1, 1, 2]), np.array([3, 2, 4])
    _, m = E25.e_process(E25.stratum_e(z, s, n))
    exp = []
    prod = {ps: Fraction(1) for ps in (2, 4, 8, 16)}
    for zi, si, ni in zip(z, s, n):
        for ps in prod:
            p = Fraction(int(si), int(ni) + 1)
            prod[ps] *= Fraction(ps) ** int(zi) / (1 - p + ps * p)
        exp.append(sum(prod.values()) / 4)
    assert np.allclose(m, [float(x) for x in exp], rtol=1e-14, atol=0)
    # by hand: M_1 = (8/5 + 16/7 + 32/11 + 64/19) / 4
    assert m[0] == pytest.approx((8 / 5 + 16 / 7 + 32 / 11 + 64 / 19) / 4, rel=1e-14)
    assert m[1] == pytest.approx((6 / 5 + 8 / 7 + 48 / 55 + 32 / 57) / 4, rel=1e-14)


def test_e_process_batch_shape_and_threshold():
    z = np.ones((5, 3), int)
    s = np.ones((5, 3), int)
    n = np.array([8, 8, 8])
    prod, m = E25.e_process(E25.stratum_e(z, s, n))
    assert prod.shape == (5, 3, 4) and m.shape == (5, 3)
    # three releases with the outcome and no control with it: M = mean over psi of (9 psi / (8 + psi))^3
    want = np.mean([(9 * ps / (8 + ps)) ** 3 for ps in E25.PSI])
    assert np.allclose(m[:, -1], want)
    assert E25.THRESHOLD == 20.0 and E25.THRESHOLD_E25B == 40.0


def test_oracle_e_value_mean_one():
    for p0 in (0.05, 0.165, 0.5):
        for ps in E25.PSI:
            assert p0 * E25.e_value_oracle(1, p0, ps) + (1 - p0) * E25.e_value_oracle(0, p0, ps) == pytest.approx(1, abs=1e-14)


def test_simulated_validity_small():
    rng = np.random.default_rng(1)
    ns = np.array([3, 5, 8, 0, 6])
    r = E25.simulate_h0(np.resize(ns, 40), "all or nothing", 0.3, 500, rng)
    assert r["share_M_ever_ge_20"] == 0.0 and r["mean_M_final"] == 1.0
    r = E25.simulate_h0(np.resize(ns, 40), "exchangeable Bernoulli", 0.2, 4000, rng)
    assert r["share_M_ever_ge_20"] <= 0.05 + 3 * np.sqrt(0.05 * 0.95 / 4000)
    assert abs(r["mean_M_1"] - 1) < 4 * r["mean_M_1_se"] + 1e-9


# ----------------------------------------------------------------------------
# Strata (section 4, rules 1 to 3)
# ----------------------------------------------------------------------------

START = "2026-10-25"


def test_candidates_same_weekday_and_eastern_clock_time_across_dst():
    # CPI on Wednesday 10 March 2027, 08:30 EST; US daylight saving time starts on 14 March 2027
    c = E25.assign_strata(rel([("2027-03-10", "cpi")]), ["2027-03-10"], START, None)
    assert list(c.k) == list(E25.WEEKS)
    et = pd.DatetimeIndex(c.anchor_utc).tz_convert(E12.TZ)
    assert (et.hour == 8).all() and (et.minute == 30).all() and (et.dayofweek == 2).all()
    utc_hours = dict(zip(c.candidate_day, pd.DatetimeIndex(c.anchor_utc).hour))
    assert utc_hours["2027-03-03"] == 13 and utc_hours["2027-03-17"] == 12  # EST before, EDT after
    assert set(c.anchor_eastern_zone) == {"EST", "EDT"}
    # DST end: a FOMC statement on 28 October 2026 (EDT), controls in November at 14:00 EST = 19:00 UTC
    c = E25.assign_strata(rel([("2026-10-28", "fomc")]), ["2026-10-28"], START, None)
    kept = c[c.kept]
    assert list(pd.DatetimeIndex(kept.anchor_utc).hour) == [19] * len(kept)


def test_two_day_rule():
    # candidates of a Wednesday release; event days one and two days from two of them
    r = rel([("2027-06-16", "fomc")])
    c = E25.assign_strata(r, ["2027-06-16", "2027-06-24", "2027-07-09"], START, None)
    by_day = c.set_index("candidate_day")
    assert not by_day.loc["2027-06-23", "kept"] and "rule 2: within 2 days" in by_day.loc["2027-06-23", "reason"]
    assert by_day.loc["2027-07-07", "kept"]  # 2 days from 9 July: not less than 2 days, kept
    assert by_day.loc["2027-06-09", "kept"]
    c2 = E25.assign_strata(r, ["2027-06-16", "2027-07-08"], START, None)
    assert not c2.set_index("candidate_day").loc["2027-07-07", "kept"]  # one day away


def test_start_rule_cal_plus_pre():
    r = rel([("2026-11-18", "fomc")])  # candidate k = -1: 11 November 2026, 14:00 EST = 19:00 UTC
    anchor = pd.Timestamp("2026-11-11 19:00", tz="UTC")
    ok = E25.assign_strata(r, ["2026-11-18"], anchor - pd.Timedelta(minutes=1800), None).set_index("k")
    assert ok.loc[-1, "kept"]
    late = E25.assign_strata(r, ["2026-11-18"], anchor - pd.Timedelta(minutes=1799), None).set_index("k")
    assert not late.loc[-1, "kept"] and "1800 minutes" in late.loc[-1, "reason"]


def test_data_end_rule():
    r = rel([("2026-11-18", "fomc")])
    last = pd.Timestamp("2026-12-16 19:00", tz="UTC")  # k = +4
    c = E25.assign_strata(r, ["2026-11-18"], START, last + pd.Timedelta(minutes=360)).set_index("k")
    assert c.loc[4, "kept"]
    c = E25.assign_strata(r, ["2026-11-18"], START, last + pd.Timedelta(minutes=359)).set_index("k")
    assert not c.loc[4, "kept"] and c.loc[4, "reason"] == "rule 2: window not inside the data"


def test_rule3_disjoint_and_order():
    # two Wednesday FOMC-time releases six weeks apart share two candidates; the earlier release keeps them
    r = rel([("2027-04-28", "fomc"), ("2027-03-17", "fomc")])  # given out of order: strata follow release time
    days = ["2027-03-17", "2027-04-28"]
    c = E25.assign_strata(r, days, START, None)
    assert list(c.groupby("stratum").date.first()) == ["2027-03-17", "2027-04-28"]
    kept = c[c.kept]
    assert not kept.anchor_utc.duplicated().any()
    late = c[c.stratum == 2].set_index("candidate_day")
    for d in ("2027-03-31", "2027-04-07"):
        assert not late.loc[d, "kept"] and late.loc[d, "reason"] == "rule 3: used by stratum 1"
        assert c[(c.stratum == 1) & (c.candidate_day == d)].kept.iloc[0]
    # no control is a release window, and the assignment is a pure function of its inputs
    assert not set(kept.anchor_utc) & set(c.release_utc)
    again = E25.assign_strata(r.iloc[::-1], days[::-1], START, None)
    pd.testing.assert_frame_equal(c, again)


def test_candidate_dropped_by_rule2_is_not_used():
    # a candidate beyond the data end for the first release is not 'used'; the same anchor is also beyond the data
    # for the second release, so neither keeps it
    r = rel([("2026-12-02", "fomc"), ("2026-12-16", "fomc")])
    end = pd.Timestamp("2026-12-24", tz="UTC")
    c = E25.assign_strata(r, ["2026-12-02", "2026-12-16"], START, end)
    assert not c[c.candidate_day == "2026-12-30"].kept.any()
    assert (c[c.candidate_day == "2026-12-30"].reason == "rule 2: window not inside the data").all()


def test_window_overlaps_cpi_and_fomc_controls_original_rules():
    # under the rules as first committed, a Wednesday CPI release and a Wednesday FOMC statement share control days:
    # windows 08:30 to 14:30 and 14:00 to 20:00 ET overlap by 30 minutes
    r = rel([("2027-01-13", "cpi"), ("2026-12-09", "fomc")])
    c = E25.assign_strata(r, ["2026-12-09", "2027-01-13"], START, None, rules="original")
    ov = E25.window_overlaps(c)
    assert len(ov) == 4 and (ov.overlap_minutes == 30).all()
    assert set(ov.kind_1) | set(ov.kind_2) == {"control"}


def test_amended_rule3_same_day_cpi_and_fomc_controls():
    """Amendment 1: the CPI release's 08:30 ET candidates on the Wednesdays already kept (14:00 ET) by the earlier FOMC
    statement are dropped; no two windows of different strata overlap."""
    r = rel([("2027-01-13", "cpi"), ("2026-12-09", "fomc")])
    c = E25.assign_strata(r, ["2026-12-09", "2027-01-13"], START, None)  # amended is the default
    assert E25.window_overlaps(c).empty
    cpi = c[c.stratum == 2].set_index("candidate_day")
    for d in ("2026-12-16", "2026-12-23", "2026-12-30", "2027-01-06"):
        assert not cpi.loc[d, "kept"] and cpi.loc[d, "reason"] == "rule 3 (Amendment 1): overlaps a control of stratum 1"
        assert c[(c.stratum == 1) & (c.candidate_day == d)].kept.iloc[0]  # the earlier stratum keeps its control
    t = E25.strata_table(c).set_index("stratum")
    assert t.loc[2, "n"] == 4 and t.loc[2, "controls_before"] == 0 and t.loc[2, "controls_after"] == 4
    assert not t.loc[2, "balanced"] and t.loc[1, "n"] == 8 and t.loc[1, "balanced"]


def test_amended_rule3_control_overlapping_a_release_window():
    """A candidate whose window overlaps the window of a later (or earlier) release is dropped by the amended rule 3,
    even when that release's day is not among the event days of rule 2 (with it, rule 2 drops the candidate first)."""
    r = rel([("2027-03-17", "fomc"), ("2027-03-24", "cpi")])
    c = E25.assign_strata(r, [], START, None)  # no event days: only rule 3 can see the release windows
    f = c[c.stratum == 1].set_index("candidate_day")
    assert not f.loc["2027-03-24", "kept"]
    assert f.loc["2027-03-24", "reason"] == "rule 3 (Amendment 1): overlaps the release window of stratum 2"
    cp = c[c.stratum == 2].set_index("candidate_day")
    assert cp.loc["2027-03-17", "reason"] == "rule 3 (Amendment 1): overlaps the release window of stratum 1"
    o = E25.assign_strata(r, [], START, None, rules="original").set_index(["stratum", "candidate_day"])
    assert o.loc[(1, "2027-03-24"), "kept"] and o.loc[(2, "2027-03-17"), "kept"]
    c2 = E25.assign_strata(r, ["2027-03-17", "2027-03-24"], START, None).set_index(["stratum", "candidate_day"])
    assert c2.loc[(1, "2027-03-24"), "reason"].startswith("rule 2: within 2 days")


def test_release_windows_never_dropped_overlap_reported():
    # a CPI release (08:30 ET) and an FOMC statement (14:00 ET) on the same day, as on 12 June 2024
    r = rel([("2024-06-12", "cpi"), ("2024-06-12", "fomc")])
    c = E25.assign_strata(r, ["2024-06-12"], "2024-01-21", "2026-01-01")
    t = E25.strata_table(c)
    assert list(t.event_type) == ["cpi", "fomc"]
    assert t.release_window_overlaps.iloc[0] == "stratum 2 (2024-06-12 18:00)"
    assert t.release_window_overlaps.iloc[1] == "stratum 1 (2024-06-12 12:30)"
    ov = E25.window_overlaps(c)
    assert len(ov) == 1 and ov.kind_1.iloc[0] == ov.kind_2.iloc[0] == "release" and ov.overlap_minutes.iloc[0] == 30
    assert not c[c.kept].anchor_utc.duplicated().any()


def test_balance_columns():
    r = rel([("2027-09-14", "cpi")])
    t = E25.strata_table(E25.assign_strata(r, ["2027-09-14"], START, None))
    assert t.controls_before.iloc[0] == 4 and t.controls_after.iloc[0] == 4 and t.balanced.iloc[0]
    t = E25.strata_table(E25.assign_strata(rel([("2026-10-28", "fomc")]), ["2026-10-28"], START, None))
    assert (t.controls_before.iloc[0], t.controls_after.iloc[0], t.balanced.iloc[0]) == (0, 4, False)
    assert E25._one_sided(t).iloc[0]
    assert (t.controls_before + t.controls_after == t.n).all()
    # no controls at all (earlier candidates before the stream, later ones after the data end): neither
    # balanced nor one-sided, since its e-value is 1 whatever happens
    t = E25.strata_table(E25.assign_strata(rel([("2026-10-28", "fomc")]), ["2026-10-28"], START, "2026-11-02"))
    assert t.n.iloc[0] == 0 and not t.balanced.iloc[0] and not E25._one_sided(t).iloc[0]


def test_planned_calendar_first_stratum():
    """The committed forward event file: the FOMC statement of 28 October 2026 keeps three controls (4, 18 and 25
    November; 11 November is one day from the CPI release of 10 November) and completes on 26 November 01:00 UTC."""
    pr = E25.planned_releases()
    c = E25.assign_strata(pr, pr["date"], E25.STREAM_START, None)
    t = E25.strata_table(c)
    first = t.iloc[0]
    assert first.date == "2026-10-28" and first.n == 3
    assert first.controls == "2026-11-04 19:00;2026-11-18 19:00;2026-11-25 19:00"
    assert first.data_needed_until == pd.Timestamp("2026-11-26 01:00", tz="UTC")
    assert "2027-01-27" not in set(t.date)
    assert not c[c.kept].anchor_utc.duplicated().any()


# ----------------------------------------------------------------------------
# Completion (rule 4) and the status part's refusals
# ----------------------------------------------------------------------------


def test_window_outcomes_half_open():
    a = pd.DatetimeIndex([pd.Timestamp("2026-11-04 19:00", tz="UTC")])
    assert E25.window_outcomes(a, a)[0]
    assert not E25.window_outcomes(a + pd.Timedelta(minutes=360), a)[0]
    assert E25.window_outcomes(a + pd.Timedelta(minutes=359), a)[0]
    assert not E25.window_outcomes(a - pd.Timedelta(minutes=1), a)[0]
    assert not E25.window_outcomes(pd.DatetimeIndex([], tz="UTC"), a)[0]


def test_window_outcomes_equal_e12_window_stats():
    """On a gap-free minute index, an up alarm in [anchor, anchor + 6 h) is E12.window_stats' up_post."""
    rng = np.random.default_rng(7)
    idx = pd.date_range("2026-01-01", periods=20_000, freq="1min", tz="UTC")
    pos = np.sort(rng.choice(np.arange(2000, 19_000), 60, replace=False))
    ev = pd.DataFrame(dict(index=pos, direction=rng.choice(["up", "down"], pos.size), max_abs_recent_over_scale=1.0,
                           post_over_baseline=1.0))
    u = rng.standard_normal(idx.size)
    mon = np.ones(idx.size, bool)
    anchors = rng.choice(np.arange(500, 19_000), 300, replace=False)
    want = [E12.window_stats(ev, mon, u, int(i))["up_post"] for i in anchors]
    up = idx[ev.loc[ev.direction == "up", "index"].to_numpy()]
    got = E25.window_outcomes(up, idx[anchors])
    assert list(got) == want


def test_completion_rule_and_prefix():
    ev = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "released"),
                 ("2026-12-09", "fomc", "scheduled")])
    r = ev[ev.status == "released"]
    c = E25.assign_strata(r, r["date"], START, None)
    t = E25.strata_table(c)
    need1 = pd.Timestamp("2026-11-26 01:00", tz="UTC")  # 25 November 14:00 EST + 6 h (DST ended on 1 November)
    assert t.data_needed_until.iloc[0] == need1
    comp = E25.completion(t, ev, need1)
    assert comp.complete.iloc[0] and comp.in_test.iloc[0]
    comp = E25.completion(t, ev, need1 - pd.Timedelta(minutes=1))
    assert not comp.complete.iloc[0] and not comp.in_test.any()
    # stratum 2 (10 November CPI) has the scheduled 9 December row within 30 days: not complete even with the data
    comp = E25.completion(t, ev, "2027-03-01")
    assert list(comp.complete) == [True, False] and "2026-12-09" in comp.not_complete_because.iloc[1]
    # a row 31 days away does not block
    ev2 = events([("2026-10-28", "fomc", "released"), ("2026-11-28", "cpi", "scheduled")])
    comp = E25.completion(E25.strata_table(E25.assign_strata(ev2.iloc[:1], ["2026-10-28"], START, None)), ev2, "2027-01-01")
    assert comp.complete.iloc[0]
    # prefix: a complete stratum after an incomplete one does not enter the test
    t2 = t.copy()
    ev3 = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "released"),
                  ("2026-09-30", "cpi", "scheduled")])
    comp = E25.completion(t2, ev3, "2027-03-01")
    assert list(comp.complete) == [False, True] and not comp.in_test.any()


def test_status_refuses_when_nothing_complete():
    ev = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "scheduled"),
                 ("2026-12-09", "fomc", "scheduled"), ("2026-12-10", "cpi", "scheduled")])
    now = pd.Timestamp("2026-11-20", tz="UTC")
    why = E25.status_refusals(ev, "2026-11-20", now=now)
    assert why and "no stratum is complete" in why[0]
    assert E25.status_refusals(ev, "2026-12-01", now=now)  # cut-off after the time of the run
    late = pd.Timestamp("2026-12-01", tz="UTC")
    why = E25.status_refusals(ev, "2026-11-26 01:00", now=late)  # data complete, but 10 November still 'scheduled'
    assert why and "scheduled rows within 30 days: 2026-11-10" in why[0]
    ev.loc[1, "status"] = "released"
    assert E25.status_refusals(ev, "2026-11-26 01:00", now=late) == []
    assert E25.status_refusals(ev, "2026-11-26 00:59", now=late)
    assert E25.status_refusals(events([]).reindex(columns=["date", "event_type", "release_time_et", "status"]),
                               "2026-12-01", now=pd.Timestamp("2026-12-02", tz="UTC"))


def test_status_part_refuses_today_before_any_download(monkeypatch, tmp_path):
    """On the committed event file the first stratum completes on 26 November 2026; the status part must refuse
    before downloading or reading anything, and write nothing."""
    from avmonitor import data as Dm

    def no_fetch(*a, **k):
        raise AssertionError("a download was attempted")

    monkeypatch.setattr(Dm, "fetch", no_fetch)
    monkeypatch.setattr(E25, "OUT", tmp_path)
    with pytest.raises(SystemExit) as err:
        E25.part_status(1, "2026-11-03")
    assert "no stratum is complete" in str(err.value)
    assert not list(tmp_path.iterdir())
    why = E25.status_refusals(E25.future_event_rows(), "2026-11-26 00:59", now=pd.Timestamp("2026-12-01", tz="UTC"))
    assert why and "data end 2026-11-26 00:59" in why[0]


def test_status_tables_report_only_complete_strata():
    ev = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "released"),
                 ("2026-12-09", "fomc", "released"), ("2026-12-10", "cpi", "scheduled")])
    d_end = pd.Timestamp("2026-12-01", tz="UTC")  # stratum 1 complete, stratum 2 needs 8 December
    cands, comp = E25.fresh_plan(ev, d_end)
    assert list(comp.in_test) == [True, False, False]
    s1 = cands[(cands.stratum == 1) & cands.kept]
    base = pd.DatetimeIndex([pd.Timestamp("2026-10-28 19:30", tz="UTC"), s1.anchor_utc.iloc[1] + pd.Timedelta(hours=1)])
    # alarms in the windows of stratum 2 (release and a control) and after the data end must not change anything
    s2 = cands[(cands.stratum == 2) & cands.kept]
    extra = pd.DatetimeIndex([pd.Timestamp("2026-11-10 14:00", tz="UTC"), s2.anchor_utc.iloc[0] + pd.Timedelta(minutes=5),
                              pd.Timestamp("2026-12-05", tz="UTC")])
    a = E25.status_tables(cands, comp, {("S2", "BTC"): base}, d_end)
    b = E25.status_tables(cands, comp, {("S2", "BTC"): base.append(extra).sort_values()}, d_end)
    for x, y in zip(a[:3], b[:3]):
        pd.testing.assert_frame_equal(x, y)
    strata, path, summ, pend = a
    assert list(strata.stratum) == [1] and list(pend.stratum) == [2, 3]
    assert not {"Z", "K", "S", "M"} & set(pend.columns)
    assert strata.Z.iloc[0] == 1 and strata.K.iloc[0] == 1 and strata.n.iloc[0] == 3
    m = np.mean([ps / (1 - 0.5 + ps * 0.5) for ps in E25.PSI])  # p = 2 / 4
    assert strata.M.iloc[0] == pytest.approx(m, rel=1e-14)
    row = summ[(summ.cell == "fomc+cpi")].iloc[0]
    assert row.role == "E25a primary" and row.threshold == 20.0 and not row.crossed
    # a window of a stratum in the test may never end after the data end
    with pytest.raises(RuntimeError):
        E25.status_tables(cands, comp.assign(in_test=True), {("S2", "BTC"): base}, d_end)


# ----------------------------------------------------------------------------
# Review round 11: event-file checks (M1 i), the ledger (M1 ii), missing minutes (m10), provenance and tags (m11)
# ----------------------------------------------------------------------------


def test_calendar_integrity_missing_cpi_row_blocks_the_stratum():
    """The reviewer's case: with the CPI row of 10 November 2026 missing, the FOMC stratum of 28 October would be
    complete with 11 November kept (n = 4); the event-file check now keeps it out of the test."""
    ev = events([("2026-10-28", "fomc", "released"), ("2026-12-09", "fomc", "released")])
    cands, comp = E25.fresh_plan(ev, "2026-12-01")
    assert comp.n.iloc[0] == 4 and comp.data_complete.iloc[0] and comp.status_complete.iloc[0]
    assert not comp.calendar_complete.iloc[0] and not comp.in_test.any()
    assert "no CPI row for 2026-11" in comp.not_complete_because.iloc[0]
    why = E25.status_refusals(ev, "2026-12-01", now=pd.Timestamp("2026-12-02", tz="UTC"))
    assert why and "no CPI row for 2026-11" in why[0]
    # the range is within 30 days of the stratum's windows, [d - 58, d + 58]: December 2026 needs its row as well
    assert "no CPI row for 2026-12" in comp.not_complete_because.iloc[0]
    ev2 = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "released"), ("2026-12-09", "fomc", "released"),
                  ("2026-12-10", "cpi", "scheduled")])
    cands, comp = E25.fresh_plan(ev2, "2026-12-01")
    assert comp.in_test.iloc[0] and comp.n.iloc[0] == 3 and comp.calendar_complete.iloc[0]
    # a cancelled CPI row is a row: the month is accounted for
    ev3 = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "cancelled"), ("2026-12-09", "fomc", "released"),
                  ("2026-12-10", "cpi", "scheduled")])
    assert E25.fresh_plan(ev3, "2026-12-01")[1].calendar_complete.iloc[0]


def test_calendar_integrity_fomc_rows_access_dates_and_schedule_coverage():
    base = [("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "released"), ("2026-12-10", "cpi", "scheduled")]
    tab = E25.strata_table(E25.assign_strata(rel([("2026-10-28", "fomc")]), ["2026-10-28", "2026-11-10"], START, None))
    # a scheduled meeting inside [d - 30, d + 30] without a row
    why = E25.calendar_integrity(tab, events(base), START, ("2026-10-28", "2026-11-18"))
    assert "no FOMC row for the scheduled meeting of 2026-11-18" in why[0]
    # present as the original date of a rescheduled row: accounted for
    ev = events(base + [("2026-11-19", "fomc", "rescheduled")])
    ev.loc[3, "original_date"] = "2026-11-18"
    assert E25.calendar_integrity(tab, ev, START, ("2026-10-28", "2026-11-18")) == [""]
    # a meeting outside the range [d - 58, d + 58] (to 25 December) is not required
    assert E25.calendar_integrity(tab, events(base), START, ("2026-10-28", "2026-12-26")) == [""]
    assert E25.calendar_integrity(tab, events(base), START, ("2026-10-28", "2026-12-25")) != [""]
    # 'released' with an access date before the release
    ev = events(base)
    ev.loc[1, "access_date"] = "2026-11-09"
    assert "2026-11-10 cpi: 'released' with an access date before the release" in E25.calendar_integrity(
        tab, ev, START, ("2026-10-28",))[0]
    ev.loc[1, "access_date"] = ""
    assert E25.calendar_integrity(tab, ev, START, ("2026-10-28",))[0]
    # the schedule must be known to d + 30: a December release needs the next year's schedule
    tab2 = E25.strata_table(E25.assign_strata(rel([("2026-12-10", "cpi")]), ["2026-12-10"], START, None))
    ev = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "released"), ("2026-12-10", "cpi", "released"),
                 ("2027-01-13", "cpi", "released"), ("2027-02-10", "cpi", "scheduled")])
    why = E25.calendar_integrity(tab2, ev, START, ("2026-10-28",))
    assert "known only to 2026-12-31; extend it to 2027-02-06" in why[0]
    assert E25.calendar_integrity(tab2, ev, START, ("2026-10-28", "2027-06-16")) == [""]
    # the committed schedule reaches December 2027, so the last planned strata wait for the 2028 schedule
    pr = E25.planned_releases()
    t = E25.strata_table(E25.assign_strata(pr, pr["date"], E25.STREAM_START, None))
    why = E25.calendar_integrity(t, E25.future_event_rows(), E25.STREAM_START)
    assert all(w == "" for w in why[:20]) and all("known only to 2027-12-31" in w for w in why[20:])


def test_calendar_integrity_blocks_later_strata_too():
    ev = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "released"),
                 ("2026-12-09", "fomc", "released"), ("2026-12-10", "cpi", "released"),
                 ("2027-01-13", "cpi", "released"), ("2027-03-10", "cpi", "released")])  # no row for February 2027
    cands, comp = E25.fresh_plan(ev, "2027-06-01", fomc_scheduled=("2026-10-28", "2026-12-09", "2027-03-17"))
    # the range of stratum 3 (9 December 2026) reaches 5 February 2027; stratum 6 (10 March) also misses April, May
    # and the March meeting
    assert list(comp.calendar_complete) == [True, True, False, False, False, False]
    assert "no CPI row for 2027-02" in comp.not_complete_because.iloc[2]
    assert "no FOMC row for the scheduled meeting of 2027-03-17" in comp.not_complete_because.iloc[5]
    assert list(comp.in_test) == [True, True, False, False, False, False]


def test_missing_minutes_in_windows():
    a = pd.DatetimeIndex([pd.Timestamp("2026-11-04 19:00", tz="UTC"), pd.Timestamp("2026-11-18 19:00", tz="UTC")])
    gaps = pd.DataFrame(dict(coin=["BTC", "BTC", "ETH"], first_missing=["2026-11-04 18:55:00+00:00",
                                                                         "2026-11-05 00:58:00+00:00",
                                                                         "2026-11-18 20:00:00+00:00"],
                             minutes=[10, 5, 3]))
    assert list(E25.missing_in_windows(gaps, "BTC", a)) == [5 + 2, 0]  # 19:00 to 19:05; 00:58 and 00:59
    assert list(E25.missing_in_windows(gaps, "ETH", a)) == [0, 3]
    assert list(E25.missing_in_windows(None, "BTC", a)) == [0, 0]
    idx = pd.date_range("2026-11-04 18:50", periods=20, freq="1min", tz="UTC")
    miss = np.zeros(20, bool)
    miss[5:15] = True
    miss[19] = True
    g = E25.gap_runs(idx, miss, "BTC")
    assert list(g.minutes) == [10, 1] and g.first_missing.iloc[0] == str(idx[5])


def test_ledger_rows_and_differences():
    strata = pd.DataFrame(dict(arm="S2", coin="BTC", stratum=[1, 2], date=["2026-10-28", "2026-11-10"],
                               event_type=["fomc", "cpi"],
                               release_utc=pd.to_datetime(["2026-10-28 18:00", "2026-11-10 13:30"], utc=True),
                               controls=["2026-11-04 19:00;2026-11-18 19:00;2026-11-25 19:00", "2026-11-03 13:30"],
                               n=[3, 1], Z=[1, 0], K=[1, 0], S=[2, 0]))
    strata = E25.add_e_process(strata)
    led = E25.ledger_rows(strata, "2026-12-09 00:00", "2026-12-09 00:00", "abc123")
    assert list(led.stratum) == ["1", "2"] and led.first_reported_head.iloc[0] == "abc123"
    assert E25.ledger_differences(None, led) == []
    assert E25.ledger_differences(led, led) == []
    # round trip through the CSV
    import io
    buf = io.StringIO()
    led.to_csv(buf, index=False)
    back = pd.read_csv(io.StringIO(buf.getvalue()), dtype=str, keep_default_na=False)
    later = E25.ledger_rows(strata, "2027-01-09 00:00", "2027-01-09 00:00", "def456")
    assert E25.ledger_differences(back, later) == []
    # a changed outcome, a changed control set, a stratum that left the test
    s2 = strata.copy()
    s2.loc[0, "K"] = 0
    s2.loc[0, "S"] = 1
    d = E25.ledger_differences(back, E25.ledger_rows(E25.add_e_process(s2), "2027-01-09", "2027-01-09", "x"))
    assert len(d) == 2 and "K 1 then, 0 now" in d[0] and "stratum 2" in d[1] and "M " in d[1]
    s3 = strata.copy()
    s3.loc[1, "controls"] = "2026-11-24 13:30"
    d = E25.ledger_differences(back, E25.ledger_rows(s3, "2027-01-09", "2027-01-09", "x"))
    assert len(d) == 1 and "controls" in d[0]
    d = E25.ledger_differences(back, E25.ledger_rows(strata.iloc[:1], "2027-01-09", "2027-01-09", "x"))
    assert len(d) == 1 and "not in the test now" in d[0]


def _fake_status(monkeypatch, tmp_path, ev, alarms, gaps=None):
    """The status part on a synthetic calendar and synthetic alarms: no download, no market data, outputs to tmp."""
    import e23_fresh_releases as E23
    import e24_blind_guard as E24
    from avmonitor import data as Dm

    monkeypatch.setattr(E25, "OUT", tmp_path)
    monkeypatch.setattr(E25, "STREAM_START", "2026-01-01")
    monkeypatch.setattr(E25, "LOAD_START", "2025-12-12")
    monkeypatch.setattr(E25, "FOMC_SCHEDULED", ("2026-01-28", "2026-03-18"))
    monkeypatch.setattr(E25, "future_event_rows", lambda: ev.copy())
    monkeypatch.setattr(Dm, "fetch", lambda url, dest, overwrite=False: dest)
    monkeypatch.setattr(E24, "_month_complete", lambda dest, m: True)
    monkeypatch.setattr(E23, "data_end", lambda paths: "2030-01-01 00:00")
    g = pd.DataFrame(columns=["coin", "first_missing", "minutes", "in_monitored_period"]) if gaps is None else gaps
    monkeypatch.setattr(E23, "data_checks", lambda paths, end: (pd.DataFrame(dict(file=["x"])),
                                                                 pd.DataFrame(dict(coin=["BTC"])), g))

    def runs(jobs, procs):
        return [dict(arm=a, coin=c, alarms=pd.DataFrame(dict(direction="up", time=pd.DatetimeIndex(alarms))))
                for c, a, _ in jobs]

    monkeypatch.setattr(E25, "_run_many", runs)


def test_status_part_ledger_tags_and_provenance(monkeypatch, tmp_path):
    ev = events([("2026-01-13", "cpi", "released"), ("2026-01-28", "fomc", "released"),
                 ("2026-02-11", "cpi", "released"), ("2026-03-11", "cpi", "released"),
                 ("2026-03-18", "fomc", "released"), ("2026-04-10", "cpi", "released")])
    alarms = [pd.Timestamp("2026-01-13 14:00", tz="UTC"), pd.Timestamp("2026-02-11 15:00", tz="UTC")]
    gaps = pd.DataFrame(dict(coin=["BTC"], first_missing=["2026-01-13 19:00:00+00:00"], minutes=[45],
                             in_monitored_period=[True]))
    _fake_status(monkeypatch, tmp_path, ev, alarms, gaps)
    # refused before any download: nothing written, no provenance
    with pytest.raises(SystemExit):
        E25.part_status(1, "2026-01-20")
    assert not list(tmp_path.iterdir())
    E25.part_status(1, "2026-03-01")
    first = sorted(p.name for p in tmp_path.iterdir())
    assert "e25_status_2026-03-01_0000_provenance.txt" in first and E25.LEDGER in first
    prov = (tmp_path / "e25_status_2026-03-01_0000_provenance.txt").read_text()
    assert "git status --porcelain" in prov and "cut-off 2026-03-01 00:00 UTC" in prov
    led = pd.read_csv(tmp_path / E25.LEDGER, dtype=str)
    assert len(led) == 2 * len(E25.STATUS_ARMS) * len(E12.COINS) and set(led.stratum) == {"1", "2"}
    st = pd.read_csv(tmp_path / "e25_status_2026-03-01_0000_strata.csv")
    btc = st[(st.coin == "BTC") & (st.arm == "S2")]
    # the gap 19:00 to 19:45 UTC on 13 January covers the last 30 minutes of the CPI window (13:30 to 19:30 UTC)
    assert btc.missing_minutes_release_window.tolist() == [30, 0] and btc.missing_minutes_total.tolist() == [30, 0]
    sm = pd.read_csv(tmp_path / "e25_status_2026-03-01_0000_summary.csv")
    assert sm[(sm.cell == "fomc+cpi")].set_index(["arm", "coin"]).missing_minutes_in_windows.to_dict() == {
        (a, c): 30 if c == "BTC" else 0 for a in E25.STATUS_ARMS for c in E12.COINS}
    # a second look at a later cut-off, with the same data: the ledger is verified, new strata added, nothing
    # overwritten
    E25.part_status(1, "2026-04-01 06:30")
    assert "e25_status_2026-04-01_0630_summary.csv" in {p.name for p in tmp_path.iterdir()}
    led2 = pd.read_csv(tmp_path / E25.LEDGER, dtype=str)
    assert set(led2.stratum) == {"1", "2", "3"} and len(led2) == 3 * len(E25.STATUS_ARMS) * len(E12.COINS)
    assert (led2[led2.stratum == "1"].first_reported_cutoff == "2026-03-01 00:00").all()
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    # a revision of the data (an alarm of stratum 1's release window disappears): refused, nothing written
    _fake_status(monkeypatch, tmp_path, ev, alarms[1:], gaps)
    with pytest.raises(SystemExit) as err:
        E25.part_status(1, "2026-04-02")
    assert "recomputed differently" in str(err.value) and "stratum 1" in str(err.value)
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before
    # a release added later inside a reported stratum's range changes its controls: refused as well
    ev2 = pd.concat([ev, events([("2026-02-02", "cpi", "rescheduled")])], ignore_index=True)
    _fake_status(monkeypatch, tmp_path, ev2, alarms, gaps)
    with pytest.raises(SystemExit) as err:
        E25.part_status(1, "2026-04-02")
    assert "recomputed differently" in str(err.value)
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


def test_status_part_reports_missing_minutes(monkeypatch, tmp_path):
    ev = events([("2026-01-13", "cpi", "released"), ("2026-01-28", "fomc", "released"),
                 ("2026-02-11", "cpi", "released"), ("2026-03-11", "cpi", "released")])
    cands, comp = E25.fresh_plan(ev, "2026-03-01", stream_start="2026-01-01", fomc_scheduled=("2026-01-28",))
    a = cands[(cands.stratum == 1) & cands.kept].anchor_utc.iloc[0]
    gaps = pd.DataFrame(dict(coin=["BTC"], first_missing=[str(a + pd.Timedelta(minutes=350))], minutes=[30],
                             in_monitored_period=[True]))
    _fake_status(monkeypatch, tmp_path, ev, [], gaps)
    E25.part_status(1, "2026-03-01")
    st = pd.read_csv(tmp_path / "e25_status_2026-03-01_0000_strata.csv", keep_default_na=False)
    btc = st[(st.coin == "BTC") & (st.arm == "S2") & (st.stratum == 1)].iloc[0]
    assert btc.missing_minutes_total == 10 and btc.windows_with_missing_minutes == 1
    assert str(btc.missing_minutes_controls).split(";")[0] == "10"
    eth = st[(st.coin == "ETH") & (st.arm == "S2")]
    assert (eth.missing_minutes_total == 0).all()
