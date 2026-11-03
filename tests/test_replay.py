"""scripts/build_replay.py: the pure functions on a short synthetic stream (no raw data needed)."""

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from avmonitor.monitor import VolatilityMonitor

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("build_replay", ROOT / "scripts" / "build_replay.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["build_replay"] = mod
    spec.loader.exec_module(mod)
    return mod


br = _load()

A, CAL, BAND = 150.0, 200, 0.25


@pytest.fixture(scope="module")
def stream():
    rng = np.random.default_rng(97)
    u = np.concatenate([rng.standard_normal(1200), 2.0 * rng.standard_normal(500), 0.6 * rng.standard_normal(900)])
    u[[50, 700, 1300]] = np.nan  # missing returns during calibration and monitoring
    u[400] = 25.0  # one jump: clipped at tau
    return u


@pytest.fixture(scope="module")
def trace(stream):
    return br.edetector_trace(stream, A, CAL, c=3.0, band=BAND)


def test_trace_matches_the_monitor(stream, trace):
    ref = VolatilityMonitor(threshold=A, c=3.0, cal_len=CAL, null_band=BAND)
    ref.run(stream)
    assert trace["events"] == [(e.index, e.direction) for e in ref.events]
    assert {d for _, d in trace["events"]} == {"up", "down"}
    assert len(trace["segments"]) == len(ref.segments)
    # state when each return arrives: calibrating for the first CAL returns and after every alarm
    mon = trace["monitoring"]
    assert not mon[:CAL].any() and mon[CAL]
    for t, _ in trace["events"]:
        assert mon[t] and not mon[t + 1 : t + 1 + CAL].any()
        assert max(trace["wu"][t], trace["wd"][t]) >= 2 * A  # the crossing value is kept
    assert np.isnan(trace["wu"][~mon]).all() and np.isfinite(trace["wu"][mon]).all()
    assert np.isnan(trace["wd"][~mon]).all() and np.isfinite(trace["wd"][mon]).all()
    below = mon.copy()
    below[[t for t, _ in trace["events"]]] = False
    assert (np.fmax(trace["wu"], trace["wd"])[below] < 2 * A).all()


def test_wealth_recomputed_from_the_clipped_statistic(stream, trace):
    """The page's Y_t and B reproduce both wealth paths with the bets of docs/math.md (null band included)."""
    Y, B, seg_of = br.clipped_statistic(stream, trace["segments"])
    lam_up = np.asarray(DEFAULT_UP_GRID)
    fr = np.asarray(DEFAULT_DOWN_FRACTIONS)
    mon = trace["monitoring"]
    su = sd = None
    for t in range(stream.size):
        if not mon[t]:
            su = sd = None
            continue
        if su is None:
            su, sd = np.zeros(lam_up.size), np.zeros(fr.size)
        if np.isfinite(Y[t]):
            lam_dn = fr / (B[t] / (1 - BAND) - 1)
            su = (1 + lam_up * (Y[t] / (1 + BAND) - 1)) * (su + 1)
            sd = (1 + lam_dn * (1 - Y[t] / (1 - BAND))) * (sd + 1)
        assert su.mean() == pytest.approx(trace["wu"][t], rel=1e-9)
        assert sd.mean() == pytest.approx(trace["wd"][t], rel=1e-9)
    ok = np.isfinite(Y)
    assert (Y[ok] <= B[ok] * (1 + 1e-12)).all() and (Y[ok] >= 0).all()
    assert np.isnan(Y[:CAL]).all() and (seg_of[:CAL] == -1).all()
    assert Y[400] == pytest.approx(B[400])  # the jump is clipped


def test_intervals_and_pauses():
    assert br.intervals([0, 1, 1, 0, 1]) == [[1, 3], [4, 5]]
    assert br.intervals([1, 1]) == [[0, 2]] and br.intervals([0, 0]) == []
    m = br.pause_mask([(2, "up"), (10, "down")], dead=3, n=12)
    assert br.intervals(m) == [[3, 6], [11, 12]]


def test_times_round_trip():
    m = np.r_[np.arange(100, 110), np.arange(1500, 1503), 3000]
    t0, rle = br.encode_times(m)
    assert t0 == 100 and rle == [[1, 9], [1391, 1], [1, 2], [1498, 1]]
    assert (br.decode_times(t0, rle) == m).all()


def _when(times, t):
    return pd.Timestamp(times[t]).tz_convert("UTC")


def test_alarm_list_check(trace):
    times = pd.date_range("2025-10-01", periods=3000, freq="min", tz="UTC")
    got = {"e-detector": br.window_alarms(trace["events"], times, 0, 3000), "other": []}
    lines = br.check_alarm_lists(got, {"e-detector": list(got["e-detector"])}, "toy")
    assert len(lines) == 2 and all(x.endswith("identical") for x in lines)
    t, d = trace["events"][0]
    with pytest.raises(AssertionError, match="only in CSV"):
        br.check_alarm_lists(got, {"e-detector": got["e-detector"], "other": [(_when(times, t), d)]}, "toy")
    flipped = [(x, "down" if y == "up" else "up") for x, y in got["e-detector"]]
    with pytest.raises(AssertionError, match="e-detector"):
        br.check_alarm_lists(got, {"e-detector": flipped}, "toy")
    with pytest.raises(AssertionError, match="not rerun"):
        br.check_alarm_lists(got, {"e-detector": got["e-detector"], "BOCPD, w=9": []}, "toy")


def test_payload_and_page(stream, trace):
    """The page holds derived statistics only, valid JSON, no NaN and no en or em dash."""
    times = pd.date_range("2025-10-01", periods=stream.size, freq="min", tz="UTC")
    Y, B, seg_of = br.clipped_statistic(stream, trace["segments"])
    comp = [(1500, "up"), (2400, "down")]
    alarms = {"e-detector": trace["events"], "BOCPD, w=250": comp}
    a, b, t_ref = 1000, 2600, 1600
    pay = br.window_payload(trace, Y, B, seg_of, alarms, CAL, times, a, b, t_ref, None, 2 * A)
    assert pay["n"] == b - a == len(pay["Y"]) == len(pay["wu"]) == len(pay["wd"]) and pay["T"] == t_ref - a
    assert (br.decode_times(pay["t0"], pay["dt"]) == np.arange(b - a) + pay["t0"]).all()
    ed = pay["methods"][0]
    assert ed["alarms"] == [[t - a, 1 if d == "up" else -1] for t, d in trace["events"] if a <= t < b]
    assert ed["paused"] == pay["recal"] and pay["methods"][1]["paused"] == [[501, 701], [1401, 1600]]
    assert set(pay) == {"t0", "dt", "n", "T", "latest", "Y", "wu", "wd", "seg", "recal", "thr", "methods"}
    pay.update(key="X1", event="X1", name="toy", title="Toy stream", subtitle="synthetic", asset="toy", market="crypto",
               tz="UTC", tzOffset=0, A=A, day=CAL, cal=CAL, band=BAND, unit="minute", facts=["A toy day."],
               init=dict(cursor=pay["T"], width=500, speed=5), marks=[["start", 0], ["run-up", 100], ["T", pay["T"]]],
               checks=["toy: identical"], n_methods=2)
    for m in pay["methods"]:
        m.update(h=None, state="monitoring", resumes_txt="", first_up_txt="", delay=None, delay_latest=None,
                 n7=[0, 0], n24=[0, 0])
    quotes = dict(quote="Quoted sentence.", header=["a", "b"], rows=[["1", "2"]])
    page = br.render_html([pay], quotes)
    assert "__DATA__" not in page and "__QUOTE__" not in page
    assert not re.search(r"[\u2013\u2014]", page)
    blob = re.search(r'<script id="replay-data" type="application/json">(.*?)</script>', page, re.S).group(1)
    assert "NaN" not in blob and "Infinity" not in blob
    data = json.loads(blob)
    assert data["days"][0]["Y"][:3] == pay["Y"][:3]
    assert "http://" not in page and "https://" not in page  # no external resource
