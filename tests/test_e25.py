"""Tests for experiments/e25_sequential_events.py [SYNTHETIC]: the e-value and the e-process, the strata of the
pre-registration (section 4), completion, and the status part's refusals. No market data are read; two tests read the
committed forward event file only."""

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
    """A synthetic event file: (date, type, status)."""
    return pd.DataFrame([dict(date=d, event_type=t, release_time_et="14:00" if t == "fomc" else "08:30", status=s)
                         for d, t, s in rows])


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


def test_window_overlaps_cpi_and_fomc_controls():
    # a Wednesday CPI release and a Wednesday FOMC statement share control days: windows 08:30 to 14:30 and
    # 14:00 to 20:00 ET overlap by 30 minutes
    r = rel([("2027-01-13", "cpi"), ("2026-12-09", "fomc")])
    c = E25.assign_strata(r, ["2026-12-09", "2027-01-13"], START, None)
    ov = E25.window_overlaps(c)
    assert len(ov) > 0 and (ov.overlap_minutes == 30).all()
    assert set(ov.kind_1) | set(ov.kind_2) == {"control"}


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
    ev = events([("2026-10-28", "fomc", "released"), ("2026-11-10", "cpi", "scheduled")])
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
                 ("2026-12-09", "fomc", "released")])
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
