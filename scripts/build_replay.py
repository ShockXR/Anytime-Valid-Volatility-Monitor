"""Minute-by-minute replay of E18's landmark case studies as one self-contained HTML page (tracker B97). [REAL]

What it does. For each chosen day it reruns E18's primary design (experiments/e18_landmark_events.py) on the same
stream, with the same frozen settings: the e-detector (VolatilityMonitor, A = 7,500 per Nifty session or 28,800 per
crypto day, c = 3, calibration 750 bars or 1,440 minutes, null band 0.25) and every competitor configuration E18 used
on that asset (BOCPD and BOCPD on truncated returns, w = 250 and 375 on Nifty, 250, 720 and 1440 on crypto; the
two-window variance-ratio scan, k = 120 and 375, or 120, 480 and 1440), with the up and down thresholds read from
results/e18_thresholds.csv (matching "outside event windows (primary)") and E18's pause after each alarm. It records
the e-detector's clipped statistic and its two wealth paths, then checks that every method's alarms from 7 days before
to 24 hours after the reference time T are exactly the rows of results/e18_alarms.csv (design "primary"), and that the
pre-registered outcomes (state at T, first up alarm, counts) equal results/e18_outcomes.csv. Any mismatch stops the
build with an error. Nothing is tuned here.

Days. N1: Nifty 50, 4 June 2024 (Indian election results, T = 09:15 IST), the main day. C3: BTC, 10 October 2025
(crypto liquidation cascade, T = 20:00 UTC), where the e-detector was recalibrating at T after an earlier alarm.

Data rule. The page embeds derived statistics only: Y_t = min(u_t^2, tau^2) / m0 (u_t the deseasonalised return),
the two wealth paths, thresholds, alarm times, pause and recalibration intervals and states. It embeds no price level
and no raw return. Raw data stay in data/raw (git-ignored); fetch them with python scripts/download_data.py.

Usage
  nice -n 10 python scripts/build_replay.py                 # N1 and C3: about 21 minutes with 2 processes (BOCPD)
  nice -n 10 python scripts/build_replay.py --days N1       # Nifty only, about 6 minutes
  python scripts/build_replay.py --cache DIR                # keep the BOCPD arrays between builds (never commit DIR);
                                                            # a build from a filled cache takes about 12 seconds
Output: docs/demo/replay.html
"""

from __future__ import annotations

import argparse
import html
import json
import math
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from avmonitor.monitor import VolatilityMonitor  # noqa: E402

RESULTS = ROOT / "results"
OUT_HTML = ROOT / "docs" / "demo" / "replay.html"
PRIMARY_MATCHING = "outside event windows (primary)"

# The days on the page. "event" and "asset" pick E18's rows; everything else comes from E18's own definitions.
DAYS = {
    "N1": dict(event="N1", asset="Nifty 50", tz="IST", tz_offset=330, title="Nifty 50, 4 June 2024",
               subtitle="Indian general election results; T = 09:15 IST, the open",
               run_up=("3 Jun 14:45 IST, before monitoring resumed", pd.Timestamp("2024-06-03 14:45", tz="Asia/Kolkata")),
               cursor_after_T=20, width=375, speed=1),
    "C3": dict(event="C3", asset="BTC", tz="UTC", tz_offset=0, title="BTC, 10 October 2025",
               subtitle="crypto liquidation cascade; T = 20:00 UTC, the US close (the tariff post came by 20:57 UTC)",
               run_up=("10 Oct 14:30 UTC, before the e-detector's alarm", pd.Timestamp("2025-10-10 14:30", tz="UTC")),
               cursor_after_T=90, width=1440, speed=5, latest_label="20:57, latest time of the post"),
}

SHORT = {"e-detector": "e-detector"}  # canvas row labels (full names are in the table under the panels)


def short_label(method: str) -> str:
    if method in SHORT:
        return SHORT[method]
    w = method.split("=")[-1]
    if method.startswith("BOCPD on truncated"):
        return f"BOCPD clip {w}"
    if method.startswith("BOCPD"):
        return f"BOCPD {w}"
    return f"scan {w}"


def long_label(method: str) -> str:
    if method == "e-detector":
        return "e-detector (VolatilityMonitor, plug-in baseline)"
    w = method.split("=")[-1]
    if method.startswith("BOCPD on truncated"):
        return f"BOCPD on returns clipped at 3 robust sd, w = {w}"
    if method.startswith("BOCPD"):
        return f"BOCPD, w = {w}"
    return f"two-window variance-ratio scan, k = {w} (CSV label: {method})"


def family(method: str) -> str:
    if method == "e-detector":
        return "edet"
    if method.startswith("BOCPD on truncated"):
        return "bocpdc"
    if method.startswith("BOCPD"):
        return "bocpd"
    return "scan"


# ----------------------------------------------------------------------------
# Pure functions (tested in tests/test_replay.py on a synthetic stream)
# ----------------------------------------------------------------------------


class TracingMonitor(VolatilityMonitor):
    """``VolatilityMonitor`` that also keeps the two wealths at the step an alarm resets them.

    Only ``_reset`` is extended, and only to read ``S_up`` and ``S_dn`` before the parent clears them; every
    decision (calibration, bets, alarms) is the parent's code, unchanged.
    """

    def _reset(self):
        if hasattr(self, "S_up"):
            self.wealth_at_reset = (float(self.S_up.mean()), float(self.S_dn.mean()))
        super()._reset()


def edetector_trace(u, threshold: float, cal_len: int, c: float = 3.0, band: float = 0.25) -> dict:
    """Stream the e-detector over ``u`` and record its state and both wealths at every step.

    ``monitoring[t]`` is the state when the return at t arrives (as in E18's ``run_edetector``); ``wu[t]`` and
    ``wd[t]`` are the increase-side and decrease-side wealths after that return (the Shiryaev-Roberts statistics
    averaged over the bet grid, the quantities compared with 2A), NaN while calibrating. At an alarm step they hold
    the values that crossed, before the reset.
    """
    u = np.asarray(u, float)
    mon = TracingMonitor(threshold=threshold, c=c, cal_len=cal_len, null_band=band)
    n = u.size
    monitoring = np.zeros(n, bool)
    wu = np.full(n, np.nan)
    wd = np.full(n, np.nan)
    for t in range(n):
        monitoring[t] = mon.calibrated
        ev = mon.update(u[t])
        if ev is not None:
            wu[t], wd[t] = mon.wealth_at_reset
        elif monitoring[t]:
            wu[t], wd[t] = float(mon.S_up.mean()), float(mon.S_dn.mean())
    return dict(events=[(e.index, e.direction) for e in mon.events], monitoring=monitoring, wu=wu, wd=wd,
                segments=list(mon.segments), starts=np.array([s["start"] for s in mon.segments], int))


def clipped_statistic(u, segments) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``Y_t = min(u_t^2, tau^2) / m0`` and ``B_t = tau^2 / m0`` with the baseline of the last segment started by t.

    While monitoring this is the statistic the bets use (divided by m0 instead of m0 (1 +/- band)). While
    recalibrating it is the same quantity against the previous baseline, which the monitor no longer bets on (the page
    draws it in grey). NaN before the first calibration ends and where the return is missing.
    """
    u = np.asarray(u, float)
    n = u.size
    starts = np.array([s["start"] for s in segments], int)
    tau2 = np.array([s["tau2"] for s in segments], float)
    m0 = np.array([s["m0"] for s in segments], float)
    k = np.searchsorted(starts, np.arange(n), side="right") - 1
    Y = np.full(n, np.nan)
    B = np.full(n, np.nan)
    ok = k >= 0
    Y[ok] = np.minimum(u[ok] ** 2, tau2[k[ok]]) / m0[k[ok]]
    B[ok] = tau2[k[ok]] / m0[k[ok]]
    return Y, B, k


def intervals(mask) -> list[list[int]]:
    """Runs of True in a boolean array as half-open ``[start, end)`` pairs."""
    m = np.asarray(mask, bool).astype(np.int8)
    d = np.diff(np.r_[0, m, 0])
    return [[int(a), int(b)] for a, b in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1))]


def pause_mask(alarms, dead: int, n: int) -> np.ndarray:
    """Steps at which a competitor is paused: ``t + 1 .. t + dead`` after each alarm at t (E13's ``scan_alarms``)."""
    m = np.zeros(n, bool)
    for t, _ in alarms:
        m[t + 1 : min(n, t + dead + 1)] = True
    return m


def window_alarms(alarms, times, a: int, b: int) -> list[tuple[pd.Timestamp, str]]:
    """Alarms with position in ``[a, b)``, as (UTC time, direction)."""
    return [(pd.Timestamp(times[t]).tz_convert("UTC"), d) for t, d in alarms if a <= t < b]


def check_alarm_lists(computed: dict, expected: dict, label: str) -> list[str]:
    """Compare per-method alarm lists; raise ``AssertionError`` naming every difference.

    ``computed`` and ``expected`` map a method name to a list of (UTC timestamp, direction) in time order. A method
    missing from ``expected`` must have no alarm (E18's CSV has rows only for alarms). Returns one line per method.
    """
    unknown = sorted(set(expected) - set(computed))
    problems = [f"{label}: method {m!r} is in the CSV but was not rerun" for m in unknown]
    lines = []
    for m, got in computed.items():
        want = expected.get(m, [])
        got_s = [(str(pd.Timestamp(t).tz_convert("UTC")), d) for t, d in got]
        want_s = [(str(pd.Timestamp(t).tz_convert("UTC")), d) for t, d in want]
        if got_s != want_s:
            extra = [x for x in got_s if x not in want_s]
            missing = [x for x in want_s if x not in got_s]
            problems.append(f"{label}, {m}: rerun has {len(got_s)} alarms, CSV {len(want_s)}; "
                            f"only in rerun {extra}, only in CSV {missing}")
        else:
            lines.append(f"{label}, {m}: {len(got_s)} alarms, identical")
    if problems:
        raise AssertionError("alarm lists differ from results/e18_alarms.csv:\n  " + "\n  ".join(problems))
    return lines


def sig(x, digits: int) -> list:
    """Round to ``digits`` significant digits for the page; NaN and infinities become None (JSON null)."""
    return [float(f"{v:.{digits}g}") if math.isfinite(v) else None for v in np.asarray(x, float)]


def encode_times(minutes) -> tuple[int, list[list[int]]]:
    """Epoch minutes as (first value, run-length coded steps [[step, count], ...])."""
    m = np.asarray(minutes, np.int64)
    d = np.diff(m)
    rle = []
    for v in d.tolist():
        if rle and rle[-1][0] == v:
            rle[-1][1] += 1
        else:
            rle.append([int(v), 1])
    return int(m[0]), rle


def decode_times(t0: int, rle) -> np.ndarray:
    steps = np.concatenate([np.full(c, s, np.int64) for s, c in rle]) if rle else np.zeros(0, np.int64)
    return np.r_[t0, t0 + np.cumsum(steps)]


def window_payload(trace: dict, Y, B, seg_of, alarms: dict, dead: int, times, a: int, b: int, t_ref: int,
                   latest: int | None, threshold: float) -> dict:
    """Derived statistics on ``[a, b)`` for the page (positions are relative to ``a``). No return enters it."""
    n_all = len(Y)
    sl = slice(a, b)
    epoch = pd.Timestamp("1970-01-01", tz="UTC")
    minutes = np.asarray((pd.DatetimeIndex(times[sl]).tz_convert("UTC") - epoch) // pd.Timedelta(minutes=1), np.int64)
    t0, rle = encode_times(minutes)
    seg = []
    for i in range(a, b):
        if seg_of[i] >= 0 and (i == a or seg_of[i] != seg_of[i - 1]):
            seg.append([i - a, float(f"{B[i]:.4g}")])
    mon = trace["monitoring"]
    recal = [[max(s - a, 0), min(e - a, b - a)] for s, e in intervals(~mon) if e > a and s < b]
    methods = []
    for m, al in alarms.items():
        if m == "e-detector":
            paused = recal
        else:
            pm = pause_mask(al, dead, n_all)
            paused = [[max(s - a, 0), min(e - a, b - a)] for s, e in intervals(pm) if e > a and s < b]
        methods.append(dict(name=m, label=short_label(m), long=long_label(m), family=family(m), paused=paused,
                            alarms=[[t - a, 1 if d == "up" else -1] for t, d in al if a <= t < b]))
    return dict(t0=t0, dt=rle, n=b - a, T=t_ref - a, latest=None if latest is None else latest - a,
                Y=sig(Y[sl], 3), wu=sig(trace["wu"][sl], 4), wd=sig(trace["wd"][sl], 4), seg=seg, recal=recal,
                thr=threshold, methods=methods)


# ----------------------------------------------------------------------------
# Rerun of E18's primary design on one stream
# ----------------------------------------------------------------------------


def _e18():
    import e18_landmark_events as E18

    return E18


def required_files(event: dict, asset: str) -> list[Path]:
    from avmonitor import data as D

    if event["market"] == "nifty":
        return [D.nifty_file_path(k) for k in D._nifty_keys("2023-10-01", "2026-01-01")]
    E18 = _e18()
    if event["stream"] != "H2 2025":
        raise ValueError(f"stream {event['stream']} is not wired into the replay")
    return [D.binance_month_path(E18.COINS[asset], m) for m in D._months("2024-01-01", "2026-01-01")]


def load_stream(event: dict, asset: str):
    """The stream E18 ran for this event and asset: (index, deseasonalised returns, first monitored position)."""
    E18 = _e18()
    if event["stream"] == "Nifty 2024-2025":
        return E18.nifty_stream()[event["stream"]][asset]
    if event["stream"] == "H2 2025":
        return E18.prepare(asset, "Jul-Dec 2025")  # E18.crypto_streams()["H2 2025"][asset]
    raise ValueError(f"stream {event['stream']} is not wired into the replay")


def bocpd_arrays(market: str, stream: str, asset: str, u, i0: int, cache: Path | None, workers: int) -> dict:
    """E18's BOCPD arrays (S, L) for one series: same inputs, hazard, calibration, caps and windows as ``bocpd_all``."""
    E18 = _e18()
    E13 = E18.E13
    cfg = E18.MARKETS[market]
    um = u[i0:][None, :]
    prior = E13.trailing_robust_var(u, cfg["day"])[i0:][None, :]
    out, todo = {}, []
    for w in sorted(cfg["ws"], reverse=True):
        for tr in (False, True):
            f = cache / f"bocpd_{stream}_{asset}_w{w}_t{int(tr)}.npz".replace(" ", "_") if cache else None
            if f is not None and f.exists():
                z = np.load(f)
                out[(w, tr)] = (z["S"], z["L"])
            else:
                todo.append(((um, prior, 1 / cfg["A"], cfg["cal"], tr, w), f))
    if todo:
        with Pool(max(1, min(workers, len(todo)))) as pool:
            for (w, tr, S, L), (_, f) in zip(pool.imap(E18.bocpd_job, [j for j, _ in todo]), todo):
                out[(w, tr)] = (S[0], L[0])
                if f is not None:
                    np.savez(f, S=S[0], L=L[0])
    return out


def read_thresholds(market: str, stream: str, asset: str) -> dict:
    t = pd.read_csv(RESULTS / "e18_thresholds.csv", float_precision="round_trip")
    t = t[(t.market == market) & (t.stream == stream) & (t.asset == asset) & (t.matching == PRIMARY_MATCHING)]
    if t.method.duplicated().any():
        raise AssertionError(f"duplicate threshold rows for {stream} {asset}")
    return {r.method: (float(r.h_up), float(r.h_down)) for r in t.itertuples()}


def expected_alarms(event: str, asset: str) -> dict:
    a = pd.read_csv(RESULTS / "e18_alarms.csv")
    a = a[(a.design == "primary") & (a.event == event) & (a.asset == asset)]
    out: dict = {}
    for r in a.itertuples():
        out.setdefault(r.method, []).append((pd.Timestamp(r.time_utc), r.direction))
    return out


OUTCOME_FIELDS = ("state_at_T", "resumes_at", "blinded_by_alarm_at", "n_up_7d_before", "n_down_7d_before",
                  "last_up_before_T", "n_up_24h_after", "n_down_24h_after", "first_alarm_after_T", "first_up_after_T",
                  "delay_min", "delay_from_latest_min")


def _same(x, y) -> bool:
    def norm(v):
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return ""
        if isinstance(v, (int, float, np.integer, np.floating)):
            return float(v)
        return str(v)

    a, b = norm(x), norm(y)
    if isinstance(a, float) and isinstance(b, float):
        return abs(a - b) < 1e-9
    return a == b


def check_outcomes(recs: dict, event: str, asset: str, label: str) -> None:
    o = pd.read_csv(RESULTS / "e18_outcomes.csv")
    o = o[(o.design == "primary") & (o.event == event) & (o.asset == asset)].set_index("method")
    problems = []
    for m, rec in recs.items():
        if m not in o.index:
            problems.append(f"{label}, {m}: no row in results/e18_outcomes.csv")
            continue
        for f in OUTCOME_FIELDS:
            if not _same(rec.get(f), o.loc[m, f]):
                problems.append(f"{label}, {m}, {f}: rerun {rec.get(f)!r}, CSV {o.loc[m, f]!r}")
    if problems:
        raise AssertionError("outcomes differ from results/e18_outcomes.csv:\n  " + "\n  ".join(problems))


def rerun_day(key: str, cache: Path | None, workers: int) -> dict:
    """Rerun E18's primary design for one day; check it against E18's CSV files; return the page payload."""
    import avmonitor.competitors as C

    E18 = _e18()
    spec = DAYS[key]
    ev = next(e for e in E18.EVENTS if e["event"] == spec["event"])
    market, stream, asset = ev["market"], ev["stream"], spec["asset"]
    cfg = E18.MARKETS[market]
    missing = [str(p) for p in required_files(ev, asset) if not p.exists()]
    if missing:
        raise SystemExit("raw data missing (run python scripts/download_data.py):\n  " + "\n  ".join(missing))
    idx, u, i0 = load_stream(ev, asset)
    um, tt = u[i0:], idx[i0:]
    print(f"{key}: stream {stream}, {asset}: {um.size} monitored steps", flush=True)

    trace = edetector_trace(um, cfg["A"], cfg["cal"], c=3.0, band=E18.BAND)
    Y, B, seg_of = clipped_statistic(um, trace["segments"])
    alarms = {"e-detector": trace["events"]}
    bocpd = bocpd_arrays(market, stream, asset, u, i0, cache, workers)
    stats_ = E18.competitor_stats(market, u, i0, bocpd)
    thr = read_thresholds(market, stream, asset)
    names = E18.method_names(market)
    if set(thr) != set(names[1:]):
        raise AssertionError(f"{key}: thresholds CSV has methods {sorted(thr)}, E18 runs {names[1:]}")
    for m in names[1:]:
        su, sd, _ = stats_[m]
        alarms[m] = C.scan_alarms(su, sd, thr[m][0], thr[m][1], cfg["cal"], cfg["cal"])

    T, latest = ev["T"], ev["latest"]
    a = int(tt.searchsorted(T - E18.D7))
    b = int(tt.searchsorted(T + E18.H24))
    t_ref = int(tt.searchsorted(T))
    computed = {m: window_alarms(al, tt, a, b) for m, al in alarms.items()}
    check_lines = check_alarm_lists(computed, expected_alarms(ev["event"], asset), key)
    recs = {}
    for m, al in alarms.items():
        kind = "e-detector" if m == "e-detector" else "competitor"
        recs[m] = E18.window_outcomes(al, tt, ev, kind, cfg["cal"], trace["monitoring"], trace["starts"])
    check_outcomes(recs, ev["event"], asset, key)
    print("\n".join(check_lines), flush=True)

    latest_i = int(tt.searchsorted(latest)) if latest is not None else None
    pay = window_payload(trace, Y, B, seg_of, alarms, cfg["cal"], tt, a, b, t_ref, latest_i, 2 * cfg["A"])
    for m in pay["methods"]:
        nm = m["name"]
        m["h"] = None if nm == "e-detector" else [float(f"{x:.6g}") for x in thr[nm]]
        r = recs[nm]
        m["state"] = r["state_at_T"]
        m["resumes_txt"] = fmt_time(pd.Timestamp(r["resumes_at"]), spec) if r["resumes_at"] else ""
        m["first_up_txt"] = "" if r["first_up_after_T"] == "missed" else fmt_time(pd.Timestamp(r["first_up_after_T"]), spec)
        m["delay"] = None if not math.isfinite(r["delay_min"]) else r["delay_min"]
        m["delay_latest"] = None if not math.isfinite(r["delay_from_latest_min"]) else r["delay_from_latest_min"]
        m["n7"] = [r["n_up_7d_before"], r["n_down_7d_before"]]
        m["n24"] = [r["n_up_24h_after"], r["n_down_24h_after"]]
    run_up = int(tt.searchsorted(spec["run_up"][1])) - a
    s1 = pd.read_csv(RESULTS / "e18_outcomes.csv")
    s1 = s1[(s1.design == "S1") & (s1.event == ev["event"]) & (s1.asset == asset) & (s1.method == "e-detector")]
    facts = day_facts(spec, ev, tt, trace, alarms, recs, a, b, t_ref, cfg)
    if len(s1) == 1 and math.isfinite(float(s1.delay_min.iloc[0])):
        facts.append(f"For comparison, E18's secondary design S1 (not a rerun here; results/e18_outcomes.csv) started "
                     f"the e-detector fresh, calibrated on the {cfg['cal']:,} {'bars' if market == 'nifty' else 'minutes'} "
                     f"just before T: it flagged this case {float(s1.delay_min.iloc[0]):.0f} minutes after T.")
    pay.update(
        latestLabel=spec.get("latest_label"), key=key, event=ev["event"], name=ev["name"], title=spec["title"],
        subtitle=spec["subtitle"], asset=asset, market=market, tz=spec["tz"], tzOffset=spec["tz_offset"], A=cfg["A"], day=cfg["day"], cal=cfg["cal"],
        band=E18.BAND, unit="bar" if market == "nifty" else "minute",
        init=dict(cursor=min(pay["T"] + spec["cursor_after_T"], pay["n"] - 1), width=spec["width"], speed=spec["speed"]),
        marks=[["7 days before T", 0], [f"Run-up: {spec['run_up'][0]}", run_up], ["T, the reference time", pay["T"]]]
        + [[f"e-detector {d} alarm, {fmt_time(tt[t], spec)}", t - a] for t, d in trace["events"] if a <= t < b],
        facts=facts,
        checks=check_lines, n_methods=len(names),
    )
    return pay


def fmt_time(x, spec, date=True) -> str:
    x = pd.Timestamp(x).tz_convert("Asia/Kolkata" if spec["tz"] == "IST" else "UTC")
    s = f"{x.day} {x:%b} {x:%H:%M}" if date else f"{x:%H:%M}"
    return f"{s} {spec['tz']}"


def _minutes(m: float) -> str:
    m = float(m)
    if abs(m) < 120:
        return f"{m:.0f} minutes"
    return f"{m / 60:.1f} hours"


def day_facts(spec, ev, tt, trace, alarms, recs, a, b, t_ref, cfg) -> list[str]:
    """Plain sentences about the day, built from the rerun (whose numbers were just checked against E18's CSVs)."""
    T = ev["T"]
    ed = recs["e-detector"]
    before = [(t, d) for t, d in trace["events"] if t < t_ref]
    out = []
    if ed["state_at_T"] == "monitoring":
        seg = trace["starts"][trace["starts"] <= t_ref][-1]
        unit = "bars" if ev["market"] == "nifty" else "minutes"
        if before:
            out.append(f"At T the e-detector was monitoring: after its {before[-1][1]} alarm of "
                       f"{fmt_time(tt[before[-1][0]], spec)} it had recalibrated on the next {cfg['cal']:,} {unit} and "
                       f"resumed at {fmt_time(tt[seg], spec)}.")
        else:
            out.append(f"At T the e-detector was monitoring, since {fmt_time(tt[seg], spec)}.")
    else:
        t_al, d_al = before[-1]
        lead = (T - tt[t_al]) / pd.Timedelta(minutes=1)
        out.append(f"At T the e-detector was recalibrating: it had raised {'an' if d_al == 'up' else 'a'} {d_al} alarm at "
                   f"{fmt_time(tt[t_al], spec)}, {_minutes(lead)} before T, and could not alarm again until "
                   f"{fmt_time(pd.Timestamp(ed['resumes_at']), spec)}. It was blind at T.")
    if ed["first_up_after_T"] != "missed":
        out.append(f"Its first up alarm after T came at {fmt_time(pd.Timestamp(ed['first_up_after_T']), spec)}, "
                   f"{ed['delay_min']:.0f} minutes after T.")
    else:
        first = ed["first_alarm_after_T"]
        nxt = f" Its next alarm was a {first.split('(')[1].rstrip(')')} alarm at {fmt_time(pd.Timestamp(first.split(' (')[0]), spec)}." if first else ""
        out.append(f"It raised no up alarm in the 24 hours after T.{nxt}")
    comp = [m for m in recs if m != "e-detector"]
    watching = [m for m in comp if recs[m]["state_at_T"] == "watching"]
    caught = [m for m in comp if recs[m]["first_up_after_T"] != "missed"]
    w_txt = ", ".join(short_label(m) for m in watching) or "none"
    out.append(f"Competitor configurations watching at T: {len(watching)} of {len(comp)} ({w_txt}); the others were "
               f"paused after their own alarms.")
    if caught:
        parts = []
        for m in caught:
            r = recs[m]
            extra = ""
            if ev["latest"] is not None and math.isfinite(r["delay_from_latest_min"]):
                extra = f", {r['delay_from_latest_min']:.0f} after 20:57"
            first_minute = r["state_at_T"] == "paused" and pd.Timestamp(r["first_up_after_T"]) == pd.Timestamp(r["resumes_at"])
            parts.append(f"{short_label(m)} at {fmt_time(pd.Timestamp(r['first_up_after_T']), spec)} "
                         f"({r['delay_min']:.0f} minutes{extra}{'; first minute after its pause' if first_minute else ''})")
        out.append(f"With an up alarm within 24 hours of T: {len(caught)} of {len(comp)}: " + "; ".join(parts) + ".")
    else:
        out.append("No competitor configuration raised an up alarm within 24 hours of T.")
    return out


# ----------------------------------------------------------------------------
# Page
# ----------------------------------------------------------------------------

E18_QUOTE = (
    "In four pre-registered case studies, the streaming e-detector flagged all three Nifty shocks, 11 to 172 minutes "
    "after the open, but none of the fifteen crypto coin and event pairs within 24 hours, because on each crypto event "
    "it had raised an up alarm 4.5 to 6.9 hours before the reference time, as the turbulence began, and was still "
    "recalibrating; the competitors, at the same alarm rates, were mostly paused for the same reason. With a fresh "
    "calibration just before the reference time it flagged all 18 cases within 15 to 207 minutes."
)
E18_ROWS = {
    "N1": "| N1, Indian election results, Nifty (09:15 IST, 4 Jun 2024) | monitoring | 11 minutes | 2 of 6 | 2 of 6 (0) |",
    "C3": "| C3, crypto liquidation cascade, five coins (20:00 UTC, 10 Oct 2025; the tariff post came by 20:57 UTC) | "
          "recalibrating on all five, after up alarms at 15:23 to 15:30 UTC that day (4.5 to 4.6 hours before T) | "
          "none on any coin | 5 of 45 (variance-ratio scan, k = 1440, on all five coins) | 10 of 45 (1); the five "
          "watching configurations within 17 to 73 minutes |",
}
E18_HEADER = ["Event (reference time T)", "e-detector at T", "e-detector: first up alarm within 24 hours of T",
              "Competitor configurations watching at T",
              "Competitor configurations with an up alarm within 24 hours (of which in the first minute after a pause ended)"]


def e18_quotes(keys) -> dict:
    """The E18 headline sentence and table rows, checked to appear verbatim in docs/report_parts/E18.md."""
    text = (ROOT / "docs" / "report_parts" / "E18.md").read_text(encoding="utf-8")
    missing = [k for k in keys if E18_ROWS[k] not in text]
    if E18_QUOTE not in text or "| " + " | ".join(E18_HEADER) + " |" not in text or missing:
        raise AssertionError("the E18 headline quoted on the page no longer appears verbatim in "
                             f"docs/report_parts/E18.md (rows {missing}); update E18_QUOTE / E18_ROWS")
    rows = [[c.strip() for c in E18_ROWS[k].strip("|").split("|")] for k in keys]
    return dict(quote=E18_QUOTE, header=E18_HEADER, rows=rows)


def render_html(days: list[dict], quotes: dict) -> str:
    """The self-contained page: CSS, the data as JSON and plain JavaScript drawing on canvas. No external resource."""
    data = json.dumps(dict(days=days), separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    head = "".join(f"<th>{html.escape(c)}</th>" for c in quotes["header"])
    rows = "".join("<tr>" + "".join(f'<td class="lead">{html.escape(c)}</td>' if j == 0 else
                                    f'<td data-label="{html.escape(h)}">{html.escape(c)}</td>'
                                    for j, (h, c) in enumerate(zip(quotes["header"], r))) + "</tr>" for r in quotes["rows"])
    checks = "".join(f"<li>{html.escape(line)}</li>" for d in days for line in d["checks"])
    out = (TEMPLATE.replace("__QUOTE__", html.escape(quotes["quote"])).replace("__HEAD__", head)
           .replace("__ROWS__", rows).replace("__CHECKS__", checks).replace("__DATA__", data))
    bad = [ch for ch in ("\u2013", "\u2014") if ch in out]
    if bad:
        raise AssertionError("the page contains an en or em dash")
    return out


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Volatility Monitor Replay</title>
<meta name="description" content="Minute-by-minute replay of pre-registered case studies of the anytime-valid volatility monitor (experiment E18, primary design, real data).">
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --edet: #2a78d6; --edetdn: #184f95; --bocpd: #eb6834; --bocpdc: #1baf7a; --scan: #6250d6;
  --bar: #7c8797; --barclip: #2f3a4a; --barcal: #dcdad3; --recal: rgba(137,135,129,0.32);
  --blind: rgba(250,178,25,0.22); --blindline: #b07800; --band: rgba(42,120,214,0.09);
  --latest: rgba(137,135,129,0.14); --pill-mon: #e3eefb; --pill-recal: #ecebe6; --pill-blind: #fde9bd;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --edet: #3987e5; --edetdn: #86b6ef; --bocpd: #d95926; --bocpdc: #199e70; --scan: #9085e9;
    --bar: #8e99aa; --barclip: #dfe5ee; --barcal: #3a3935; --recal: rgba(137,135,129,0.38);
    --blind: rgba(250,178,25,0.20); --blindline: #fab219; --band: rgba(57,135,229,0.14);
    --latest: rgba(137,135,129,0.18); --pill-mon: #16304f; --pill-recal: #2c2c2a; --pill-blind: #4a3608;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --edet: #3987e5; --edetdn: #86b6ef; --bocpd: #d95926; --bocpdc: #199e70; --scan: #9085e9;
  --bar: #8e99aa; --barclip: #dfe5ee; --barcal: #3a3935; --recal: rgba(137,135,129,0.38);
  --blind: rgba(250,178,25,0.20); --blindline: #fab219; --band: rgba(57,135,229,0.14);
  --latest: rgba(137,135,129,0.18); --pill-mon: #16304f; --pill-recal: #2c2c2a; --pill-blind: #4a3608;
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body { margin: 0; background: var(--page); color: var(--ink); font: 16px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1120px; margin: 0 auto; padding: 20px 16px 48px; }
h1 { font-size: 1.6rem; line-height: 1.25; margin: 0 0 4px; }
h2 { font-size: 1.2rem; margin: 32px 0 8px; }
p, li { color: var(--ink); }
.sub { color: var(--ink2); margin: 0; }
.real { display: inline-block; font-size: 0.75rem; font-weight: 600; letter-spacing: 0.02em; padding: 1px 7px; border-radius: 999px; border: 1px solid var(--border); background: var(--surface); color: var(--ink2); vertical-align: 2px; }
.top { display: flex; gap: 12px; justify-content: space-between; align-items: flex-start; }
.btn, select { font: inherit; font-size: 0.9rem; color: var(--ink); background: var(--surface); border: 1px solid var(--border); border-radius: 8px; padding: 6px 10px; min-height: 36px; }
.btn { cursor: pointer; }
.btn:hover, select:hover { border-color: var(--axis); }
.btn:focus-visible, select:focus-visible, canvas:focus-visible { outline: 2px solid var(--edet); outline-offset: 2px; }
.btn.primary { background: var(--ink); color: var(--surface); border-color: var(--ink); min-width: 92px; }
.lede { max-width: 72ch; }
.card { background: var(--surface); border: 1px solid var(--border); border-radius: 12px; padding: 14px; margin-top: 16px; }
.days { display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 10px; }
.days .btn[aria-pressed="true"] { border-color: var(--edet); box-shadow: inset 0 0 0 1px var(--edet); }
.daytitle { margin: 2px 0 10px; }
.daytitle b { font-size: 1.05rem; }
.daytitle span { color: var(--ink2); }
.controls { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
.controls label { font-size: 0.85rem; color: var(--ink2); display: inline-flex; gap: 6px; align-items: center; max-width: 100%; }
.controls select { max-width: 100%; min-width: 0; }
#jump { width: 15em; }
.caption { margin: 12px 0 6px; padding: 10px 12px; border-radius: 10px; background: var(--page); border: 1px solid var(--border); font-size: 0.92rem; min-height: 7.2em; }
.caption .when { font-weight: 600; }
.caption .rel { color: var(--ink2); font-variant-numeric: tabular-nums; margin-left: 6px; }
.caption div + div { margin-top: 3px; }
.pill { display: inline-block; padding: 0 8px; border-radius: 999px; font-size: 0.82rem; font-weight: 600; margin-right: 4px; }
.pill.mon { background: var(--pill-mon); }
.pill.recal { background: var(--pill-recal); }
.pill.blind { background: var(--pill-blind); }
.pill.alarm { background: var(--ink); color: var(--surface); }
.meter { display: inline-block; width: 110px; height: 8px; border-radius: 4px; background: var(--grid); vertical-align: 1px; margin: 0 6px; overflow: hidden; }
.meter span { display: block; height: 100%; background: var(--edet); }
.meter.dn span { background: var(--edetdn); }
.num { font-variant-numeric: tabular-nums; }
.panel { margin-top: 10px; }
.ptitle { display: flex; flex-wrap: wrap; gap: 4px 14px; align-items: baseline; font-size: 0.85rem; color: var(--ink2); margin-bottom: 2px; }
.ptitle b { color: var(--ink); font-weight: 600; }
.key { display: inline-flex; align-items: center; gap: 5px; white-space: nowrap; }
.sw { display: inline-block; width: 18px; height: 0; border-top: 2px solid var(--edet); }
.sw.dash { border-top: 2px dashed var(--edetdn); }
.sw.thr { border-top: 1.5px dotted var(--ink2); }
.sw.box { width: 12px; height: 10px; border: 0; background: var(--recal); opacity: 0.6; }
.sw.blind { width: 12px; height: 10px; border: 1px solid var(--blindline); background: var(--blind); }
.sw.band { width: 12px; height: 10px; border: 0; background: var(--band); }
.sw.bar { width: 6px; height: 10px; border: 0; background: var(--bar); }
.sw.barcal { width: 6px; height: 10px; border: 0; background: var(--barcal); }
canvas { display: block; width: 100%; touch-action: pan-y; cursor: crosshair; }
#tip { position: fixed; z-index: 5; pointer-events: none; background: var(--surface); color: var(--ink); border: 1px solid var(--border); border-radius: 8px; padding: 6px 9px; font-size: 0.8rem; line-height: 1.4; box-shadow: 0 2px 10px rgba(0,0,0,0.12); display: none; max-width: 260px; }
.tablewrap { overflow-x: auto; margin-top: 12px; }
table { border-collapse: collapse; font-size: 0.85rem; width: 100%; }
th, td { text-align: left; padding: 6px 8px; border-bottom: 1px solid var(--grid); vertical-align: top; }
th { color: var(--ink2); font-weight: 600; }
td.num, th.num { text-align: right; }
td.lead .note { font-weight: 400; font-size: 0.8rem; }
.dot { display: inline-block; width: 9px; height: 9px; border-radius: 2px; margin-right: 6px; vertical-align: 0; }
ul.facts li + li { margin-top: 4px; }
blockquote { margin: 10px 0; padding: 8px 14px; border-left: 3px solid var(--axis); color: var(--ink); background: var(--surface); }
.note { color: var(--ink2); font-size: 0.9rem; }
code { font-size: 0.88em; background: var(--surface); border: 1px solid var(--border); border-radius: 4px; padding: 0 4px; overflow-wrap: anywhere; }
.math { font-family: "Cambria Math", "STIX Two Math", Georgia, serif; white-space: nowrap; }
dl { margin: 0; }
dt { font-weight: 600; margin-top: 10px; }
dd { margin: 2px 0 0; }
footer { margin-top: 36px; font-size: 0.85rem; color: var(--ink2); }
footer ul { padding-left: 18px; }
@media (max-width: 600px) {
  body { font-size: 15px; }
  .controls label.jump { flex: 1 1 100%; }
  #jump { width: auto; flex: 1 1 auto; }
  table.stack thead { display: none; }
  table.stack tr { display: block; padding: 8px 0; border-bottom: 1px solid var(--grid); }
  table.stack td { display: flex; gap: 10px; justify-content: space-between; border: 0; padding: 2px 0; text-align: right; }
  table.stack td::before { content: attr(data-label); color: var(--ink2); text-align: left; flex: 0 0 42%; }
  table.stack td.lead { display: block; text-align: left; font-weight: 600; }
  table.stack td.lead::before { content: none; }
  h1 { font-size: 1.3rem; }
  .card { padding: 10px; border-radius: 10px; }
  .btn, select { padding: 5px 8px; font-size: 0.85rem; }
  .btn.primary { min-width: 76px; }
  .meter { width: 70px; }
  .caption { font-size: 0.86rem; min-height: 9.5em; }
}
</style>
</head>
<body>
<main>
<div class="top">
  <div>
    <h1>Replay of a landmark market day</h1>
    <p class="sub">Anytime-valid volatility monitor, experiment E18 (pre-registered case studies), primary design <span class="real">[REAL]</span></p>
  </div>
  <button class="btn" id="theme" type="button" aria-label="Switch between light and dark theme">Dark</button>
</div>

<p class="lede">This page replays, minute by minute, what the streaming volatility monitor (the <i>e-detector</i>) and the competing change-point methods saw on a day of a large, well-remembered market shock. The numbers come from a rerun of E18's primary design on market data <span class="real">[REAL]</span>; the alarm lists were checked against E18's saved results when the page was built. This is <b>one pre-registered case study, not evidence that the method is better</b>: E18 has four events, no pooled test, and does not rank the methods.</p>

<section class="card" aria-label="Replay">
  <div class="days" id="days" role="group" aria-label="Choose a day"></div>
  <p class="daytitle"><b id="dtitle"></b> <span id="dsub"></span></p>
  <div class="controls">
    <button class="btn" id="back" type="button" aria-label="Step back one minute">&#8249; Step</button>
    <button class="btn primary" id="play" type="button">&#9654; Play</button>
    <button class="btn" id="fwd" type="button" aria-label="Step forward one minute">Step &#8250;</button>
    <label>Speed <select id="speed" aria-label="Minutes per frame"></select></label>
    <label>Window <select id="width" aria-label="Width of the time window"></select></label>
    <label class="jump">Jump to <select id="jump" aria-label="Jump to a moment"></select></label>
  </div>
  <div class="caption" id="caption" aria-live="polite"></div>

  <div class="panel">
    <div class="ptitle"><b>1. Clipped, scaled squared return Y<sub>t</sub></b>
      <span class="key"><span class="sw bar"></span>monitoring</span>
      <span class="key"><span class="sw barcal"></span>recalibrating (not used)</span>
      <span class="key"><span class="sw band"></span>null band 0.75 to 1.25</span>
      <span class="key"><span class="sw thr"></span>clip level B</span></div>
    <canvas id="c1" data-h="130" tabindex="0" aria-label="Panel 1: clipped scaled statistic"></canvas>
  </div>
  <div class="panel">
    <div class="ptitle"><b>2. e-process wealth, log scale</b>
      <span class="key"><span class="sw"></span>increase side</span>
      <span class="key"><span class="sw dash"></span>decrease side</span>
      <span class="key"><span class="sw thr"></span>alarm threshold 2A</span>
      <span class="key"><span class="sw box"></span>recalibrating</span>
      <span class="key"><span class="sw blind"></span>blind at T</span></div>
    <canvas id="c2" data-h="170" tabindex="0" aria-label="Panel 2: wealth of the two sides"></canvas>
  </div>
  <div class="panel">
    <div class="ptitle"><b>3. Alarms, one row per method</b>
      <span class="key">&#9650; up alarm</span><span class="key">&#9661; down alarm</span>
      <span class="key"><span class="sw box"></span>paused or recalibrating</span></div>
    <canvas id="c3" data-h="180" tabindex="0" aria-label="Panel 3: alarms of every method"></canvas>
  </div>
  <div class="panel">
    <div class="ptitle"><b>Context: 7 days before T to 24 hours after</b><span>increase-side wealth; the box is the window above; click or tap to move the cursor</span></div>
    <canvas id="c4" data-h="64" aria-label="Context strip of the whole replay"></canvas>
  </div>
  <div class="tablewrap">
    <p class="note">Methods on this asset, as configured in E18 (thresholds from <code>results/e18_thresholds.csv</code>; outcomes as in <code>results/e18_outcomes.csv</code>) <span class="real">[REAL]</span></p>
    <table id="mtable" class="stack"><thead></thead><tbody></tbody></table>
  </div>
</section>

<h2>What happened on this day <span class="real">[REAL]</span></h2>
<ul class="facts" id="facts"></ul>

<h2>How to read the panels</h2>
<dl>
  <dt>1. Clipped statistic</dt>
  <dd>Each bar is one minute (one bar of the Nifty session): <span class="math">Y<sub>t</sub> = min(u<sub>t</sub><sup>2</sup>, &tau;<sup>2</sup>) / m<sub>0</sub></span>, where u<sub>t</sub> is the one-minute log return divided by a predictable intraday volatility profile, &tau; is 3 robust standard deviations of the calibration window, and m<sub>0</sub> is the mean of min(u<sup>2</sup>, &tau;<sup>2</sup>) on that window. If volatility has not changed, Y<sub>t</sub> averages about 1. No bar can exceed <span class="math">B = &tau;<sup>2</sup> / m<sub>0</sub></span> (dotted line), so one jump adds only a bounded amount of evidence. Grey bars fall in a recalibration: the monitor is not betting, and they show the same quantity against the previous baseline.</dd>
  <dt>2. e-process</dt>
  <dd>The wealth of the increase side and of the decrease side, on a log scale. Each side is a Shiryaev-Roberts e-process: a sum, over every possible start time, of the wealth of a bettor who started then. The increase side bets that Y<sub>t</sub> averages above 1.25, the decrease side that it averages below 0.75 (the null band, shaded in panel 1). An alarm is raised when either side reaches 2A. Hatched: the monitor is recalibrating after an alarm; amber: the recalibration that covered the reference time T, when the monitor was blind.</dd>
  <dt>3. Alarms</dt>
  <dd>One row per method, with exactly the alarms that E18 recorded for this asset from 7 days before to 24 hours after T (<code>results/e18_alarms.csv</code>; the build stops if any differs). Filled triangles are up alarms (volatility rose), open triangles down alarms. A tinted span means the method was paused after its own alarm (750 bars or 1,440 minutes, as in E13 and E18), a hatched span on the e-detector's row that it was recalibrating.</dd>
  <dt>Caption</dt>
  <dd>The clock time at the cursor (IST for Nifty, UTC for crypto), the monitor's state, and each side's wealth against the threshold 2A (the bars show it on a log scale from 1 to 2A). Everything after the cursor is drawn faintly. Keys: space plays or pauses, the arrow keys step.</dd>
</dl>

<h2>E18 headline, quoted from <code>docs/report_parts/E18.md</code> <span class="real">[REAL]</span></h2>
<blockquote>__QUOTE__</blockquote>
<div class="tablewrap"><table class="stack"><thead><tr>__HEAD__</tr></thead><tbody>__ROWS__</tbody></table></div>
<p class="note">Case studies, not evidence: with four events, and five coins that move together, E18 has no pooled statistic, no p-value, and does not rank the methods. The thresholds of the competitors were matched to the e-detector's alarm counts using the whole half-year (or two years) outside the event windows, which a user could not have done in real time; that favours the competitors.</p>

<h2>The monitor in brief</h2>
<dl>
  <dt>Calibration</dt>
  <dd>On n returns (750 bars, two Nifty sessions; 1,440 minutes for crypto) set the robust scale s = 1.4826 &times; median |u|, the clip &tau; = 3s, and m<sub>0</sub> = mean of min(u<sup>2</sup>, &tau;<sup>2</sup>).</dd>
  <dt>Increase side</dt>
  <dd><span class="math">M<sub>t</sub><sup>up</sup> = (1/6) &Sigma;<sub>&lambda;</sub> &Sigma;<sub>j&le;t</sub> &Pi;<sub>i=j..t</sub> (1 + &lambda; (Y<sub>i</sub>/1.25 &minus; 1))</span>, with bets &lambda; in {0.02, 0.05, 0.1, 0.2, 0.35, 0.5}. If the conditional mean of Y<sub>i</sub> is at most 1.25, each product (the wealth of a bettor who starts at j with 1) is a nonnegative supermartingale, so a large M<sup>up</sup> is evidence that the mean of Y has risen above 1.25.</dd>
  <dt>Decrease side</dt>
  <dd>The same with factors 1 + &lambda; (1 &minus; Y<sub>i</sub>/0.75), &lambda; = f / (B/0.75 &minus; 1), f in {0.1, 0.25, 0.5, 0.75, 0.95}. Because Y is bounded by B, a fall in volatility can be detected at all: without a bound, no nonparametric e-test can detect a fall in the variance (<code>docs/math.md</code>, section 5.3).</dd>
  <dt>Alarm and guarantee</dt>
  <dd>Alarm when either side reaches 2A, with A = (steps per day) / 0.05: 7,500 per 375-bar Nifty session, 28,800 per 1,440-minute crypto day. If the clipped squared returns keep the baseline's conditional mean and the baseline is the true one, the chance of any false alarm within the first day after a calibration is at most m/A = 5% (Theorem B in <code>docs/math.md</code>). Here the baseline is estimated on the calibration window (plug-in), for which the bound is not proved, and real returns do not meet the conditions exactly: on real data the false alarm rate is measured, not guaranteed.</dd>
  <dt>After an alarm</dt>
  <dd>The monitor recalibrates on the next n returns and cannot alarm meanwhile. This blind period decided the crypto cases of E18.</dd>
  <dt>Competitors</dt>
  <dd>BOCPD (Bayesian online change point detection, Adams and MacKay 2007) for the variance of Gaussian returns, whose alarm statistic is the posterior log odds that the current run is shorter than w steps, signed by whether the short runs have the higher variance; run on the returns and on returns clipped at 3 robust standard deviations of the previous day; and a two-window variance-ratio scan that compares the realised variance of the last k returns with that of the k before (inspired by Bibinger, Jirak and Vetter 2017; not their test). Each pauses for n steps after its own alarm, and its up and down thresholds were set so that it raises as many up and as many down alarms as the e-detector on the same stream, counting only alarms outside the windows from 8 days before to 2 days after each event.</dd>
</dl>

<footer>
  <p>Built by <code>scripts/build_replay.py</code>, which reruns E18's primary design (<code>experiments/e18_landmark_events.py</code>) on the stream of each day with the frozen settings and the thresholds of <code>results/e18_thresholds.csv</code>, and stops if any alarm differs from <code>results/e18_alarms.csv</code> or any pre-registered outcome from <code>results/e18_outcomes.csv</code>. Pre-registration: <code>docs/preregistration_e18_events.md</code>. The page holds derived statistics only (Y<sub>t</sub>, the two wealth paths, thresholds, alarm times and states): no price level and no raw return. Checks at build time:</p>
  <ul>__CHECKS__</ul>
</footer>
</main>
<div id="tip" role="status"></div>
<script id="replay-data" type="application/json">__DATA__</script>
<script>
(function () {
'use strict';
var DATA = JSON.parse(document.getElementById('replay-data').textContent);
var MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
var WDAY = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'];
var FAM = { edet: '--edet', bocpd: '--bocpd', bocpdc: '--bocpdc', scan: '--scan' };
var TICK_MS = 110;
var D = null, cur = 0, width = 0, speed = 1, timer = null, v0 = 0, drag = null, hatch = null;
var COL = {};
function $(id) { return document.getElementById(id); }
function readColors() {
  var cs = getComputedStyle(document.documentElement);
  ['page', 'surface', 'ink', 'ink2', 'muted', 'grid', 'axis', 'edet', 'edetdn', 'bocpd', 'bocpdc', 'scan', 'bar', 'barclip',
   'barcal', 'recal', 'blind', 'blindline', 'band', 'latest'].forEach(function (k) { COL[k] = cs.getPropertyValue('--' + k).trim(); });
  var p = document.createElement('canvas'); p.width = 6; p.height = 6;
  var g = p.getContext('2d'); g.strokeStyle = COL.recal; g.lineWidth = 1;
  g.beginPath(); g.moveTo(-1, 7); g.lineTo(7, -1); g.moveTo(5, 7); g.lineTo(7, 5); g.moveTo(-1, 1); g.lineTo(1, -1); g.stroke();
  hatch = p;
}
function pad2(x) { return (x < 10 ? '0' : '') + x; }
function local(i) { return new Date((D.t[i] + D.tzOffset) * 60000); }
function clock(i) { var d = local(i); return pad2(d.getUTCHours()) + ':' + pad2(d.getUTCMinutes()); }
function dlabel(i) { var d = local(i); return d.getUTCDate() + ' ' + MON[d.getUTCMonth()]; }
function short(i) { return dlabel(i) + ' ' + clock(i) + ' ' + D.tz; }
function full(i) { var d = local(i); return WDAY[d.getUTCDay()] + ' ' + d.getUTCDate() + ' ' + MON[d.getUTCMonth()] + ' ' + d.getUTCFullYear() + ', ' + clock(i) + ' ' + D.tz; }
function fmtNum(v) {
  if (v === null || v === undefined) return 'none';
  if (v >= 1000) return Math.round(v).toLocaleString('en-US');
  if (v >= 100) return v.toFixed(0);
  if (v >= 10) return v.toFixed(1);
  return v.toFixed(2);
}
function relT(i) {
  var m = D.t[i] - D.t[D.T];
  if (m === 0) return 'at T';
  var s = m > 0 ? 'T + ' : 'T - ', a = Math.abs(m);
  if (a < 120) return s + a + ' min';
  if (a < 2880) return s + Math.floor(a / 60) + ' h ' + pad2(a % 60) + ' min';
  return s + Math.floor(a / 1440) + ' d ' + Math.floor((a % 1440) / 60) + ' h';
}
function prep(d) {
  if (d.t) return;
  var t = new Float64Array(d.n), k = 1, j;
  t[0] = d.t0;
  d.dt.forEach(function (sc) { for (j = 0; j < sc[1]; j++) { t[k] = t[k - 1] + sc[0]; k++; } });
  d.t = t;
  d.B = new Float64Array(d.n).fill(NaN);
  d.seg.forEach(function (s, q) { var e = q + 1 < d.seg.length ? d.seg[q + 1][0] : d.n; for (j = s[0]; j < e; j++) d.B[j] = s[1]; });
  d.mon = new Uint8Array(d.n).fill(1);
  d.recal.forEach(function (r) { for (j = r[0]; j < r[1]; j++) d.mon[j] = 0; });
  d.blind = null;
  d.recal.forEach(function (r) { if (r[0] <= d.T && d.T < r[1]) d.blind = r; });
  d.sess = [];
  for (j = 0; j < d.n; j++) if (j === 0 ? d.market === 'nifty' : t[j] - t[j - 1] > 1) d.sess.push(j);
  var bmax = 0, wmax = d.thr;
  for (j = 0; j < d.n; j++) {
    if (d.B[j] > bmax) bmax = d.B[j];
    if (d.wu[j] !== null && d.wu[j] > wmax) wmax = d.wu[j];
    if (d.wd[j] !== null && d.wd[j] > wmax) wmax = d.wd[j];
  }
  d.ymax = bmax * 1.06;
  d.whi = Math.ceil(Math.log10(wmax) + 0.35);
  d.edet = d.methods[0];
  d.alarmAt = {};
  d.methods.forEach(function (m, q) { m.alarms.forEach(function (a) { (d.alarmAt[a[0]] = d.alarmAt[a[0]] || []).push([q, a[1]]); }); });
}
function inIv(ivs, i) { for (var q = 0; q < ivs.length; q++) if (ivs[q][0] <= i && i < ivs[q][1]) return ivs[q]; return null; }

/* geometry */
function gut(W) {
  var c = $('c3').getContext('2d'), wmax = 0, f = (W < 560 ? '10px' : '11.5px') + ' system-ui, sans-serif';
  c.font = '600 ' + f;
  D.methods.forEach(function (m) { wmax = Math.max(wmax, c.measureText(m.label).width); });
  c.font = '600 11px system-ui, sans-serif';
  wmax = Math.max(wmax, c.measureText(W < 560 ? '2A' : '2A = 15,000').width - 14);
  return Math.ceil(wmax) + 24;
}
var G = { L: 106, R: 10, pw: 1, w: 1 };
function xEdge(i) { return G.L + (i - v0) * G.pw / G.w; }
function xOf(i) { return xEdge(i + 0.5); }
function setup(cv, h) {
  var dpr = window.devicePixelRatio || 1, W = cv.clientWidth, H = h || +cv.dataset.h;
  cv.style.height = H + 'px';
  if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
  var ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = COL.surface; ctx.fillRect(0, 0, W, H);
  return { ctx: ctx, W: W, H: H };
}
function twoPass(ctx, W, H, fn) {
  var xc = xOf(cur) + 0.5;
  ctx.save(); ctx.beginPath(); ctx.rect(G.L, 0, Math.max(0, xc - G.L), H); ctx.clip(); fn(); ctx.restore();
  ctx.save(); ctx.beginPath(); ctx.rect(xc, 0, W, H); ctx.clip(); ctx.globalAlpha = 0.2; fn(); ctx.restore();
}
function clipPlot(ctx, W, H, fn) { ctx.save(); ctx.beginPath(); ctx.rect(G.L, 0, W - G.L - G.R, H); ctx.clip(); fn(); ctx.restore(); }
function shade(ctx, ivs, y0, y1, style) {
  ivs.forEach(function (r) {
    var a = Math.max(r[0], v0), b = Math.min(r[1], v0 + G.w);
    if (b <= a) return;
    ctx.fillStyle = style; ctx.fillRect(xEdge(a), y0, xEdge(b) - xEdge(a), y1 - y0);
  });
}
function backdrop(ctx, y0, y1, withTicks) {
  var pat = ctx.createPattern(hatch, 'repeat');
  shade(ctx, D.recal, y0, y1, pat);
  if (D.blind) {
    shade(ctx, [D.blind], y0, y1, COL.blind);
  }
  if (D.latest !== null) shade(ctx, [[D.T, D.latest]], y0, y1, COL.latest);
  ctx.strokeStyle = COL.axis; ctx.lineWidth = 1;
  D.sess.forEach(function (i) { if (i > v0 && i < v0 + G.w) { var x = Math.round(xEdge(i)) + 0.5; ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y1); ctx.stroke(); } });
  if (withTicks) {
    ctx.strokeStyle = COL.grid;
    ticks().forEach(function (tk) { var x = Math.round(xEdge(tk.i)) + 0.5; ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y1); ctx.stroke(); });
  }
}
function markers(ctx, y0, y1) {
  if (D.T >= v0 && D.T < v0 + G.w) {
    var x = Math.round(xEdge(D.T)) + 0.5;
    ctx.save(); ctx.strokeStyle = COL.ink; ctx.setLineDash([4, 3]); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y1); ctx.stroke(); ctx.restore();
  }
  var xc = Math.round(xOf(cur)) + 0.5;
  ctx.save(); ctx.strokeStyle = COL.ink2; ctx.lineWidth = 1.5; ctx.globalAlpha = 0.85;
  ctx.beginPath(); ctx.moveTo(xc, y0); ctx.lineTo(xc, y1); ctx.stroke(); ctx.restore();
}
function ticks() {
  var px = G.pw / G.w, steps = [1, 2, 5, 10, 15, 30, 60, 120, 180, 360, 720, 1440, 2880], s = 2880, q, i;
  for (q = 0; q < steps.length; q++) if (steps[q] * px >= 62) { s = steps[q]; break; }
  var out = [], last = -1e9, isS = {};
  D.sess.forEach(function (j) { isS[j] = 1; });
  for (i = v0; i < v0 + G.w && i < D.n; i++) {
    var lt = D.t[i] + D.tzOffset, mid = lt % 1440 === 0, ses = !!isS[i];
    if (ses || mid || (s < 1440 && lt % s === 0)) {
      var x = xEdge(i);
      if (x - last < 44 && !ses) continue;
      if (ses && out.length && x - last < 44) out.pop();
      out.push({ i: i, label: (ses || mid) ? dlabel(i) : clock(i), major: ses || mid });
      last = x;
    }
  }
  return out;
}
function tickLabels(ctx, y) {
  ctx.font = '11px system-ui, sans-serif'; ctx.textBaseline = 'top'; ctx.textAlign = 'left';
  ticks().forEach(function (tk) {
    ctx.fillStyle = tk.major ? COL.ink2 : COL.muted;
    ctx.font = (tk.major ? '600 ' : '') + '11px system-ui, sans-serif';
    ctx.fillText(tk.label, xEdge(tk.i) + 3, y);
  });
}
function gutterText(ctx, txt, y, color, bold) {
  ctx.fillStyle = color || COL.muted; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  ctx.font = (bold ? '600 ' : '') + '11px system-ui, sans-serif';
  ctx.fillText(txt, G.L - 6, y);
}

/* panel 1 */
function drawY() {
  var o = setup($('c1')), ctx = o.ctx, W = o.W, H = o.H, top = 6, bot = 4, ph = H - top - bot;
  var yOf = function (v) { return top + ph * (1 - v / D.ymax); };
  clipPlot(ctx, W, H, function () {
    backdrop(ctx, top, H - bot, true);
    ctx.fillStyle = COL.band; ctx.fillRect(G.L, yOf(1 + D.band), W - G.L - G.R, yOf(1 - D.band) - yOf(1 + D.band));
  });
  ctx.strokeStyle = COL.grid; ctx.lineWidth = 1;
  [0, 1].forEach(function (v) { var y = Math.round(yOf(v)) + 0.5; ctx.beginPath(); ctx.moveTo(G.L, y); ctx.lineTo(W - G.R, y); ctx.stroke(); });
  gutterText(ctx, '0', yOf(0) - 5); gutterText(ctx, '1', yOf(1));
  var bw = Math.max(1, G.pw / G.w - (G.pw / G.w > 3 ? 1 : 0));
  twoPass(ctx, W, H, function () {
    for (var i = v0; i < v0 + G.w && i < D.n; i++) {
      var v = D.Y[i]; if (v === null) continue;
      ctx.fillStyle = !D.mon[i] ? COL.barcal : (v >= D.B[i] * 0.999 ? COL.barclip : COL.bar);
      var y = yOf(v); ctx.fillRect(xEdge(i), y, bw, yOf(0) - y);
    }
    ctx.save(); ctx.strokeStyle = COL.ink2; ctx.setLineDash([2, 3]); ctx.lineWidth = 1.2;
    D.seg.forEach(function (s, q) {
      var a = Math.max(s[0], v0), b = Math.min(q + 1 < D.seg.length ? D.seg[q + 1][0] : D.n, v0 + G.w);
      if (b <= a) return; var y = Math.round(yOf(s[1])) + 0.5;
      ctx.beginPath(); ctx.moveTo(xEdge(a), y); ctx.lineTo(xEdge(b), y); ctx.stroke();
    });
    ctx.restore();
  });
  var bNow = D.B[Math.min(Math.max(cur, 0), D.n - 1)];
  if (isFinite(bNow)) gutterText(ctx, 'B = ' + bNow.toFixed(1), yOf(bNow), COL.ink2, true);
  markers(ctx, top, H - bot);
  ctx.font = '600 11px system-ui, sans-serif'; ctx.textBaseline = 'top'; ctx.textAlign = 'left'; ctx.fillStyle = COL.ink;
  if (D.T >= v0 && D.T < v0 + G.w) ctx.fillText('T', xEdge(D.T) + 3, top + 12);
  if (D.latest !== null && D.latest > v0 && D.latest < v0 + G.w && D.latestLabel) {
    ctx.font = '11px system-ui, sans-serif'; ctx.fillStyle = COL.ink2;
    var xl = Math.round(xEdge(D.latest)) + 0.5, lw = ctx.measureText(D.latestLabel).width;
    if (xl + 3 + lw > W - G.R) { ctx.textAlign = 'right'; ctx.fillText(D.latestLabel, xl - 3, top + 26); ctx.textAlign = 'left'; }
    else ctx.fillText(D.latestLabel, xl + 3, top + 26);
    ctx.strokeStyle = COL.ink2; ctx.setLineDash([1, 2]);
    ctx.beginPath(); ctx.moveTo(xl, top); ctx.lineTo(xl, H - bot); ctx.stroke(); ctx.setLineDash([]);
  }
}

/* panel 2 */
function drawW() {
  var o = setup($('c2')), ctx = o.ctx, W = o.W, H = o.H, top = 8, bot = 6, ph = H - top - bot, lo = -1, hi = D.whi;
  var yOf = function (v) { return top + ph * (1 - (Math.log10(Math.max(v, 0.1)) - lo) / (hi - lo)); };
  clipPlot(ctx, W, H, function () { backdrop(ctx, top, H - bot, true); });
  ctx.lineWidth = 1;
  var yt = Math.round(yOf(D.thr)) + 0.5;
  for (var e = lo; e <= hi; e++) {
    var y = Math.round(yOf(Math.pow(10, e))) + 0.5;
    ctx.strokeStyle = COL.grid; ctx.beginPath(); ctx.moveTo(G.L, y); ctx.lineTo(W - G.R, y); ctx.stroke();
    var lab = e < 0 ? '0.1' : e < 3 ? String(Math.pow(10, e)) : e < 6 ? Math.pow(10, e - 3) + 'k' : Math.pow(10, e - 6) + 'M';
    if (Math.abs(y - yt) > 11) gutterText(ctx, lab, y);
  }
  ctx.save(); ctx.strokeStyle = COL.ink2; ctx.setLineDash([1.5, 3]); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(G.L, yt); ctx.lineTo(W - G.R, yt); ctx.stroke(); ctx.restore();
  gutterText(ctx, W < 560 ? '2A' : '2A = ' + Math.round(D.thr).toLocaleString('en-US'), yt, COL.ink2, true);
  function line(arr, color, dash, lw) {
    ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = lw; ctx.setLineDash(dash); ctx.lineJoin = 'round';
    ctx.beginPath(); var pen = false;
    for (var i = Math.max(v0 - 1, 0); i < v0 + G.w + 1 && i < D.n; i++) {
      var v = arr[i]; if (v === null) { pen = false; continue; }
      var x = xOf(i), y = yOf(v);
      if (pen) ctx.lineTo(x, y); else { ctx.moveTo(x, y); pen = true; }
    }
    ctx.stroke(); ctx.restore();
  }
  twoPass(ctx, W, H, function () {
    line(D.wd, COL.edetdn, [5, 3], 1.5);
    line(D.wu, COL.edet, [], 2);
    D.edet.alarms.forEach(function (a) {
      var i = a[0]; if (i < v0 || i >= v0 + G.w) return;
      var v = a[1] > 0 ? D.wu[i] : D.wd[i];
      ctx.beginPath(); ctx.arc(xOf(i), yOf(v), 4.5, 0, 2 * Math.PI);
      ctx.fillStyle = a[1] > 0 ? COL.edet : COL.surface; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = a[1] > 0 ? COL.surface : COL.edetdn; ctx.stroke();
    });
  });
  if (D.blind && D.blind[1] > v0 && D.blind[0] < v0 + G.w) {
    var xa = Math.max(xEdge(D.blind[0]), G.L) + 4;
    ctx.fillStyle = COL.blindline; ctx.font = '600 11px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillText('blind at T: recalibrating', Math.min(xa, W - G.R - 150), top + 2);
  }
  markers(ctx, top, H - bot);
}

/* panel 3 */
function rowH(W) { return W < 560 ? 19 : 21; }
function tri(ctx, x, y, up, color, r) {
  ctx.beginPath();
  if (up) { ctx.moveTo(x, y - r); ctx.lineTo(x + r, y + r * 0.8); ctx.lineTo(x - r, y + r * 0.8); }
  else { ctx.moveTo(x, y + r); ctx.lineTo(x + r, y - r * 0.8); ctx.lineTo(x - r, y - r * 0.8); }
  ctx.closePath();
  if (up) { ctx.fillStyle = color; ctx.fill(); ctx.lineWidth = 1.5; ctx.strokeStyle = COL.surface; ctx.stroke(); }
  else { ctx.fillStyle = COL.surface; ctx.fill(); ctx.lineWidth = 1.6; ctx.strokeStyle = color; ctx.stroke(); }
}
function drawA() {
  var cv = $('c3'), Wc = cv.clientWidth, rh = rowH(Wc), top = 6, bot = 22, nm = D.methods.length;
  var o = setup(cv, top + nm * rh + bot), ctx = o.ctx, W = o.W, H = o.H;
  var pat = ctx.createPattern(hatch, 'repeat');
  clipPlot(ctx, W, H, function () {
    backdrop(ctx, top, top + nm * rh, true);
  });
  D.methods.forEach(function (m, q) {
    var y0 = top + q * rh, yc = y0 + rh / 2, col = COL[m.family === 'edet' ? 'edet' : m.family];
    clipPlot(ctx, W, H, function () {
      if (m.family === 'edet') shade(ctx, m.paused, y0 + 2, y0 + rh - 2, pat);
      else { ctx.save(); ctx.globalAlpha = 0.16; shade(ctx, m.paused, y0 + 3, y0 + rh - 3, col); ctx.restore(); }
      ctx.strokeStyle = COL.grid; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(G.L, Math.round(yc) + 0.5); ctx.lineTo(W - G.R, Math.round(yc) + 0.5); ctx.stroke();
    });
    ctx.fillStyle = col; ctx.fillRect(4, yc - 4, 8, 8);
    ctx.fillStyle = m.family === 'edet' ? COL.ink : COL.ink2; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    ctx.font = (m.family === 'edet' ? '600 ' : '') + (W < 560 ? '10px' : '11.5px') + ' system-ui, sans-serif';
    ctx.fillText(m.label, 16, yc);
  });
  twoPass(ctx, W, H, function () {
    D.methods.forEach(function (m, q) {
      var yc = top + q * rh + rh / 2, col = COL[m.family === 'edet' ? 'edet' : m.family];
      m.alarms.forEach(function (a) { if (a[0] >= v0 && a[0] < v0 + G.w) tri(ctx, xOf(a[0]), yc, a[1] > 0, col, rh < 20 ? 5 : 5.5); });
    });
  });
  markers(ctx, top, top + nm * rh);
  tickLabels(ctx, top + nm * rh + 5);
  if (D.T >= v0 && D.T < v0 + G.w) {
    ctx.fillStyle = COL.ink; ctx.font = '600 11px system-ui, sans-serif'; ctx.textAlign = 'right'; ctx.textBaseline = 'top';
    ctx.fillText('T', xEdge(D.T) - 3, top + nm * rh + 5);
  }
}

/* context strip */
function drawC() {
  var o = setup($('c4')), ctx = o.ctx, W = o.W, H = o.H, top = 4, bot = 18, ph = H - top - bot, L = G.L, pw = W - L - G.R, n = D.n;
  var xe = function (i) { return L + i * pw / n; };
  var lo = -1, hi = D.whi, yOf = function (v) { return top + ph * (1 - (Math.log10(Math.max(v, 0.1)) - lo) / (hi - lo)); };
  var pat = ctx.createPattern(hatch, 'repeat');
  D.recal.forEach(function (r) { ctx.fillStyle = pat; ctx.fillRect(xe(r[0]), top, xe(r[1]) - xe(r[0]), ph); });
  if (D.blind) { ctx.fillStyle = COL.blind; ctx.fillRect(xe(D.blind[0]), top, xe(D.blind[1]) - xe(D.blind[0]), ph); }
  ctx.strokeStyle = COL.ink2; ctx.setLineDash([1.5, 3]); ctx.lineWidth = 1;
  var yt = Math.round(yOf(D.thr)) + 0.5; ctx.beginPath(); ctx.moveTo(L, yt); ctx.lineTo(W - G.R, yt); ctx.stroke(); ctx.setLineDash([]);
  ctx.strokeStyle = COL.edet; ctx.lineWidth = 1.2; ctx.beginPath();
  var pen = false, px = -1, best = null;
  for (var i = 0; i <= n; i++) {
    var x = i < n ? Math.floor(xe(i + 0.5)) : -2, v = i < n ? D.wu[i] : null;
    if (x !== px || i === n) {
      if (best !== null) { var yy = yOf(best); if (pen) ctx.lineTo(px + 0.5, yy); else { ctx.moveTo(px + 0.5, yy); pen = true; } } else pen = false;
      px = x; best = null;
    }
    if (v !== null && (best === null || v > best)) best = v;
  }
  ctx.stroke();
  D.edet.alarms.forEach(function (a) { tri(ctx, xe(a[0] + 0.5), top + 5, a[1] > 0, COL.edet, 4); });
  ctx.strokeStyle = COL.ink; ctx.setLineDash([3, 2]);
  var xT = Math.round(xe(D.T)) + 0.5; ctx.beginPath(); ctx.moveTo(xT, top); ctx.lineTo(xT, top + ph); ctx.stroke(); ctx.setLineDash([]);
  ctx.strokeStyle = COL.ink2; ctx.lineWidth = 1.2; ctx.fillStyle = COL.band;
  ctx.fillRect(xe(v0), top, xe(v0 + G.w) - xe(v0), ph); ctx.strokeRect(Math.round(xe(v0)) + 0.5, top + 0.5, Math.max(2, xe(v0 + G.w) - xe(v0)), ph - 1);
  var xc = Math.round(xe(cur + 0.5)) + 0.5; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(xc, top); ctx.lineTo(xc, top + ph); ctx.stroke();
  ctx.font = '11px system-ui, sans-serif'; ctx.textBaseline = 'top'; ctx.textAlign = 'left'; ctx.fillStyle = COL.muted;
  var last = -1e9;
  for (i = 0; i < n; i++) {
    var lt = D.t[i] + D.tzOffset, start = D.market === 'nifty' ? (i === 0 || D.t[i] - D.t[i - 1] > 1) : lt % 1440 === 0;
    if (start) { var xx = xe(i); if (xx - last >= 48) { ctx.fillText(dlabel(i), xx + 2, top + ph + 3); last = xx; } }
  }
  gutterText(ctx, 'whole replay', top + ph / 2, COL.ink2);
}

/* caption, table, controls */
function caption() {
  var i = cur, h = '<div><span class="when">' + full(i) + '</span><span class="rel">' + relT(i) + '</span></div>';
  var rc = inIv(D.recal, i);
  if (rc) {
    var al = null; D.edet.alarms.forEach(function (a) { if (a[0] === rc[0] - 1) al = a[1] > 0 ? 'up' : 'down'; });
    var why = al ? 'after its ' + al + ' alarm at ' + short(rc[0] - 1) : 'after an alarm before this window';
    var until = rc[1] < D.n ? 'can alarm again from ' + short(rc[1]) : 'resumes after the end of this replay';
    var blind = D.blind === rc;
    h += '<div><span class="pill ' + (blind ? 'blind' : 'recal') + '">' + (blind ? '&#9888; Recalibrating, blind at T' : 'Recalibrating') + '</span>' +
      why + '; collecting ' + D.cal.toLocaleString('en-US') + ' ' + D.unit + 's for a new baseline, cannot alarm; ' + until + '.</div>';
    h += '<div>No wealth: the monitor is not betting. Y<sub>t</sub> against the previous baseline: <span class="num">' + fmtNum(D.Y[i]) + '</span>.</div>';
  } else {
    var prev = null; D.recal.forEach(function (r) { if (r[1] <= i) prev = r; });
    h += '<div><span class="pill mon">Monitoring</span>' + (prev ? 'since ' + short(prev[1]) : 'since before this window') + '. Y<sub>t</sub> = <span class="num">' +
      fmtNum(D.Y[i]) + '</span> (clip B = ' + D.B[i].toFixed(1) + ').</div>';
    var up = D.wu[i], dn = D.wd[i], lt = Math.log10(D.thr);
    var f = function (v) { return Math.max(0, Math.min(1, Math.log10(Math.max(v || 1, 1)) / lt)) * 100; };
    h += '<div>Increase side <b class="num">' + fmtNum(up) + '</b><span class="meter"><span style="width:' + f(up).toFixed(1) + '%"></span></span>' +
      ' decrease side <b class="num">' + fmtNum(dn) + '</b><span class="meter dn"><span style="width:' + f(dn).toFixed(1) + '%"></span></span>' +
      ' threshold 2A = <span class="num">' + Math.round(D.thr).toLocaleString('en-US') + '</span></div>';
  }
  var watching = 0, nc = D.methods.length - 1;
  D.methods.forEach(function (m, q) { if (q > 0 && !inIv(m.paused, i)) watching++; });
  var here = (D.alarmAt[i] || []).map(function (qa) {
    return '<span class="pill alarm">' + (qa[1] > 0 ? '&#9650; UP' : '&#9661; DOWN') + '</span>' + D.methods[qa[0]].label;
  });
  h += '<div>Competitor configurations watching: <b class="num">' + watching + ' of ' + nc + '</b>' + (here.length ? '. Alarm at this ' + D.unit + ': ' + here.join(', ') : '') + '.</div>';
  $('caption').innerHTML = h;
}
function table() {
  var th = '<tr><th>Method</th><th class="num">Thresholds (up / down)</th><th class="num">Alarms in 7 days before T (up / down)</th><th>State at T</th><th>First up alarm within 24 h of T</th></tr>';
  $('mtable').querySelector('thead').innerHTML = th;
  var rows = D.methods.map(function (m) {
    var col = COL[m.family === 'edet' ? 'edet' : m.family];
    var thr = m.h ? m.h[0] + ' / ' + m.h[1] : 'wealth 2A = ' + Math.round(D.thr).toLocaleString('en-US');
    var st = m.state + (m.resumes_txt ? ' (back ' + m.resumes_txt + ')' : '');
    var fu = m.delay === null ? 'missed' : m.first_up_txt + ', ' + m.delay.toFixed(0) + ' min' + (m.delay_latest !== null ? ' (' + m.delay_latest.toFixed(0) + ' from 20:57)' : '');
    return '<tr><td class="lead"><span class="dot" style="background:' + col + '"></span>' + m.label + '<br><span class="note">' + m.long + '</span></td>' +
      '<td class="num" data-label="Thresholds (up / down)">' + thr + '</td><td class="num" data-label="Alarms in 7 days before T (up / down)">' + m.n7[0] + ' / ' + m.n7[1] +
      '</td><td data-label="State at T">' + st + '</td><td data-label="First up alarm within 24 h of T">' + fu + '</td></tr>';
  });
  $('mtable').querySelector('tbody').innerHTML = rows.join('');
}
function facts() {
  var ul = $('facts'); ul.innerHTML = '';
  D.facts.forEach(function (s) { var li = document.createElement('li'); li.textContent = s; ul.appendChild(li); });
}
function clampV(v) { return Math.max(0, Math.min(v, D.n - G.w)); }
function follow(frac) {
  if (frac !== undefined) { v0 = clampV(cur - Math.round(G.w * frac)); return; }
  if (cur < v0 + G.w * 0.05 && v0 > 0) v0 = clampV(cur - Math.round(G.w * 0.7));
  else if (cur >= v0 + G.w * 0.92) v0 = clampV(cur - Math.round(G.w * 0.3));
}
function draw() {
  var W = $('c1').clientWidth;
  G.L = gut(W); G.R = 10; G.pw = W - G.L - G.R; G.w = Math.min(width, D.n);
  if (!drag) follow();
  drawY(); drawW(); drawA(); drawC(); caption();
}
var raf = 0;
function redraw() { if (!raf) raf = requestAnimationFrame(function () { raf = 0; draw(); }); }
function setCur(i, frac) { cur = Math.max(0, Math.min(D.n - 1, i)); if (frac !== undefined) { G.w = Math.min(width, D.n); follow(frac); } redraw(); }
function stop() { if (timer) { clearInterval(timer); timer = null; } $('play').innerHTML = '&#9654; Play'; }
function play() {
  if (timer) { stop(); return; }
  if (cur >= D.n - 1) setCur(D.marks[1][1], 0.3);
  $('play').innerHTML = '&#10074;&#10074; Pause';
  timer = setInterval(function () { if (cur >= D.n - 1) { stop(); return; } setCur(cur + speed); }, TICK_MS);
}
function options(sel, list, value) {
  sel.innerHTML = '';
  list.forEach(function (o) { var e = document.createElement('option'); e.value = o[1]; e.textContent = o[0]; if (+o[1] === +value) e.selected = true; sel.appendChild(e); });
}
function selectDay(k) {
  stop();
  D = DATA.days[k]; prep(D);
  Array.prototype.forEach.call($('days').children, function (b, q) { b.setAttribute('aria-pressed', q === k ? 'true' : 'false'); });
  $('dtitle').textContent = D.title + ' (' + D.event + ').';
  $('dsub').textContent = D.subtitle.charAt(0).toUpperCase() + D.subtitle.slice(1) + '.';
  cur = D.init.cursor; width = D.init.width; speed = D.init.speed;
  var u = D.unit + 's';
  options($('speed'), [1, 2, 5, 15, 60].map(function (s) { return [s + ' ' + (s === 1 ? D.unit : u) + ' per frame', s]; }), speed);
  options($('width'), D.market === 'nifty'
    ? [['1 session', 375], ['2 sessions', 750], ['whole replay', D.n]]
    : [['6 hours', 360], ['1 day', 1440], ['2 days', 2880], ['whole replay', D.n]], width);
  options($('jump'), [['choose a moment', -1]].concat(D.marks), -1);
  G.w = Math.min(width, D.n); follow(0.7);
  table(); facts(); redraw();
}
function indexAt(cv, ev) {
  var r = cv.getBoundingClientRect(), x = ev.clientX - r.left;
  if (cv.id === 'c4') return Math.floor((x - G.L) / (r.width - G.L - G.R) * D.n);
  return v0 + Math.floor((x - G.L) / G.pw * G.w);
}
function tip(ev, i) {
  var t = $('tip');
  if (i === null || i < 0 || i >= D.n) { t.style.display = 'none'; return; }
  var rc = inIv(D.recal, i);
  var s = '<b>' + short(i) + '</b> (' + relT(i) + ')<br>' + (rc ? 'recalibrating' : 'monitoring') + '<br>Y<sub>t</sub> ' + fmtNum(D.Y[i]);
  if (!rc) s += '<br>wealth up ' + fmtNum(D.wu[i]) + ', down ' + fmtNum(D.wd[i]);
  (D.alarmAt[i] || []).forEach(function (qa) { s += '<br>' + (qa[1] > 0 ? '&#9650; up' : '&#9661; down') + ' alarm: ' + D.methods[qa[0]].label; });
  t.innerHTML = s; t.style.display = 'block';
  var x = ev.clientX + 14, y = ev.clientY + 14, tw = t.offsetWidth, thh = t.offsetHeight;
  if (x + tw > window.innerWidth - 8) x = ev.clientX - tw - 14;
  if (y + thh > window.innerHeight - 8) y = ev.clientY - thh - 14;
  t.style.left = Math.max(4, x) + 'px'; t.style.top = Math.max(4, y) + 'px';
}
function hookCanvas(cv) {
  cv.addEventListener('pointerdown', function (ev) {
    var i = indexAt(cv, ev); if (i < 0 || i >= D.n) return;
    stop(); drag = cv.id;
    if (cv.id === 'c4') { drag = null; setCur(i, 0.6); return; }
    try { cv.setPointerCapture(ev.pointerId); } catch (e) { /* ignore */ }
    setCur(i);
  });
  cv.addEventListener('pointermove', function (ev) {
    var i = indexAt(cv, ev);
    if (drag === cv.id) { setCur(Math.max(v0, Math.min(v0 + G.w - 1, i))); }
    if (ev.pointerType === 'mouse') tip(ev, cv.id === 'c4' ? null : i);
  });
  cv.addEventListener('pointerup', function () { drag = null; redraw(); });
  cv.addEventListener('pointercancel', function () { drag = null; });
  cv.addEventListener('pointerleave', function () { $('tip').style.display = 'none'; });
}
function init() {
  readColors();
  var dark = matchMedia('(prefers-color-scheme: dark)').matches;
  $('theme').textContent = dark ? 'Light' : 'Dark';
  $('theme').addEventListener('click', function () {
    var isDark = document.documentElement.getAttribute('data-theme') === 'dark' || (!document.documentElement.getAttribute('data-theme') && matchMedia('(prefers-color-scheme: dark)').matches);
    document.documentElement.setAttribute('data-theme', isDark ? 'light' : 'dark');
    $('theme').textContent = isDark ? 'Dark' : 'Light';
    readColors(); table(); redraw();
  });
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', function () { readColors(); table(); redraw(); });
  DATA.days.forEach(function (d, k) {
    var b = document.createElement('button'); b.className = 'btn'; b.type = 'button';
    b.textContent = d.event + ': ' + d.title; b.addEventListener('click', function () { selectDay(k); }); $('days').appendChild(b);
  });
  $('play').addEventListener('click', play);
  $('back').addEventListener('click', function () { stop(); setCur(cur - 1); });
  $('fwd').addEventListener('click', function () { stop(); setCur(cur + 1); });
  $('speed').addEventListener('change', function () { speed = +this.value; });
  $('width').addEventListener('change', function () { width = +this.value; G.w = Math.min(width, D.n); follow(0.7); redraw(); });
  $('jump').addEventListener('change', function () { var v = +this.value; if (v >= 0) { stop(); setCur(v, 0.3); } this.value = -1; });
  ['c1', 'c2', 'c3', 'c4'].forEach(function (id) { hookCanvas($(id)); });
  document.addEventListener('keydown', function (ev) {
    var tg = ev.target.tagName; if (tg === 'SELECT' || tg === 'INPUT' || tg === 'BUTTON' && ev.key === ' ') return;
    if (ev.key === ' ') { ev.preventDefault(); play(); }
    else if (ev.key === 'ArrowRight') { stop(); setCur(cur + (ev.shiftKey ? 10 : 1)); }
    else if (ev.key === 'ArrowLeft') { stop(); setCur(cur - (ev.shiftKey ? 10 : 1)); }
  });
  window.addEventListener('resize', redraw);
  selectDay(0);
}
init();
})();
</script>
</body>
</html>
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--days", nargs="+", default=list(DAYS), choices=list(DAYS))
    ap.add_argument("--cache", type=Path, default=None, help="folder for BOCPD arrays (derived, never committed)")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", type=Path, default=OUT_HTML)
    args = ap.parse_args(argv)
    if args.cache:
        args.cache.mkdir(parents=True, exist_ok=True)
    quotes = e18_quotes(args.days)
    payloads = [rerun_day(k, args.cache, args.workers) for k in args.days]
    html = render_html(payloads, quotes)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html, encoding="utf-8")
    print(f"wrote {args.out} ({args.out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
