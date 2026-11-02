"""A pre-registered live run of the monitor on one Nifty 50 session.

This is a demonstration, not evidence: one session cannot test anything. Its point is to show
the monitor working on a session it has never seen, with every setting and every reported
quantity fixed in a commit made before the session opened.

Pre-registration (see live/PREREGISTRATION_2026-11-03.md, committed before 09:15 IST on
3 November 2026):

* Data: Nifty 50 one-minute bars from the public mirror used in E7 (`avmonitor.data.load_nifty`),
  cleaned exactly as in E7 (`e07_nifty.load_clean`: regular sessions, frozen-price runs of 5 or
  more zero returns set to missing).
* Deseasonalising: E7's predictable session profile (median |r| per minute over the previous 20
  sessions, 5-minute centred smoothing, current session excluded).
* Monitor: E7's pre-registered primary configuration. `VolatilityMonitor` with A = 7500
  (375 / 0.05, each side alarming at 2A), c = 3, calibration 750 returns, null band 0.25,
  plug-in baseline, recalibration after each alarm.
* Stream: every regular session from 1 October 2026 to the target session, in order, so the
  monitor's state at the open of the target session is set by the preceding sessions only.
  The intraday profile uses sessions from 1 August 2026 onwards.
* Reported for the target session, and nothing else:
  1. every alarm in the session: time (IST), direction, and the baseline scale it was measured
     against;
  2. whether the monitor was monitoring or calibrating at 09:15, and the share of the session's
     minutes it spent monitoring;
  3. the session's robust volatility (scaled median |u| of the deseasonalised returns) divided by
     the baseline scale in force at the open (or at the first moment the monitor was monitoring);
  4. for each alarm, the robust volatility of the rest of the session after the alarm divided by
     the alarm's baseline (above 1 after an up alarm, below 1 after a down alarm would agree with
     the alarm);
  5. the SHA-256 of the data file used, so the run can be repeated byte for byte.

Usage:
    python experiments/live_nifty_preregistered.py config                   # print the frozen settings and their hash
    python experiments/live_nifty_preregistered.py run --day 2026-11-03     # run on a session once its bars are out
    python experiments/live_nifty_preregistered.py run --day 2026-11-02 --dry-run   # pipeline check on a past session
Outputs: results/live/<day>/ (alarms.csv, summary.json), or results/live/dryrun_<day>/ with --dry-run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from avmonitor import data as D  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402
from avmonitor.robust import MAD_NORMAL_SCALE  # noqa: E402

import e07_nifty as E7  # noqa: E402

CONFIG = dict(
    data="Nifty 50 one-minute bars, public GitHub mirror (avmonitor.data.load_nifty), cleaned as in E7 load_clean",
    profile=dict(kind="E7 session_profile", lookback_sessions=E7.LOOKBACK, smooth_minutes=E7.SMOOTH,
                 profile_data_from="2026-08-01"),
    monitor=dict(cls="VolatilityMonitor", threshold_A=E7.A, c=E7.C, cal_len=E7.PRIMARY[0], null_band=E7.PRIMARY[1],
                 baseline="plugin", c_down=None, recalibrate_after_alarm=True),
    stream_from="2026-10-01",
    steps_per_session=E7.STEPS_PER_SESSION,
    reported=["alarms in the target session (time, direction, baseline scale)",
              "state at 09:15 and share of the session monitored",
              "session robust volatility / baseline scale at the open (or first monitored minute)",
              "rest-of-session robust volatility / alarm baseline, per alarm",
              "SHA-256 of the data files used"],
)


def config_hash() -> str:
    return hashlib.sha256(json.dumps(CONFIG, sort_keys=True).encode()).hexdigest()


def robust_scale(u: np.ndarray) -> float:
    u = u[np.isfinite(u)]
    return float(np.median(np.abs(u)) * MAD_NORMAL_SCALE) if u.size else float("nan")


def run(day: str, dry_run: bool) -> dict:
    target = pd.Timestamp(day).date()
    end = (pd.Timestamp(day) + pd.Timedelta(days=1)).strftime("%Y-%m-%d")
    # make sure the month files are present (the current month's file is updated daily upstream)
    D.download_nifty(CONFIG["profile"]["profile_data_from"], end)
    df, info, _ = E7.load_clean(start=CONFIG["profile"]["profile_data_from"], end=end)
    sessions = sorted({s.date() for s in pd.to_datetime(df["session"].unique())})
    if target not in sessions:
        raise SystemExit(f"session {target} not in the data yet (last session: {sessions[-1]})")
    u_all = E7.deseasonalise(df)
    keep = df.index >= pd.Timestamp(CONFIG["stream_from"], tz=E7.TZ)
    df_s, u = df[keep], u_all[keep].to_numpy()
    sess = pd.to_datetime(df_s["session"]).dt.date.to_numpy()
    in_target = sess == target

    m = CONFIG["monitor"]
    mon = VolatilityMonitor(threshold=m["threshold_A"], c=m["c"], cal_len=m["cal_len"], null_band=m["null_band"])
    monitoring = np.zeros(u.size, bool)
    base = np.full(u.size, np.nan)
    for t, x in enumerate(u):
        mon.update(x)
        if mon.calibrated:
            monitoring[t] = True
            base[t] = mon.scale
    idx_t = np.flatnonzero(in_target)
    times = df_s.index
    alarms = []
    for e in mon.events:
        if in_target[e.index]:
            rest = u[e.index + 1: idx_t[-1] + 1]
            alarms.append(dict(time_ist=times[e.index].strftime("%Y-%m-%d %H:%M"), direction=e.direction,
                               baseline_scale=float(e.baseline_scale),
                               rest_of_session_minutes=int(np.isfinite(rest).sum()),
                               rest_of_session_ratio=robust_scale(rest) / float(e.baseline_scale)))
    first_mon = idx_t[monitoring[idx_t]]
    open_base = float(base[first_mon[0]]) if first_mon.size else float("nan")
    files = sorted({str(D.nifty_file_path(k)) for k in D._nifty_keys(CONFIG["profile"]["profile_data_from"], end)})
    summary = dict(
        session=str(target), dry_run=dry_run, config_sha256=config_hash(), config=CONFIG,
        returns_in_session=int(np.isfinite(u[idx_t]).sum()),
        state_at_open="monitoring" if monitoring[idx_t[0]] else "calibrating",
        share_of_session_monitored=float(monitoring[idx_t].mean()),
        first_monitored_minute=times[first_mon[0]].strftime("%H:%M") if first_mon.size else None,
        session_robust_vol_over_open_baseline=robust_scale(u[idx_t]) / open_base if first_mon.size else None,
        alarms_in_session=len(alarms),
        alarms=alarms,
        alarms_in_stream_before_session=int(sum(e.index < idx_t[0] for e in mon.events)),
        stream_sessions=int(len(set(sess))),
        data_files={Path(f).name: hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in files},
        cleaning=dict(sessions=info["sessions"], frozen_sessions_dropped=info["frozen_sessions_dropped"],
                      flagged_returns=info["flagged_returns"]),
    )
    out = ROOT / "results" / "live" / (f"dryrun_{target}" if dry_run else str(target))
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(alarms, columns=["time_ist", "direction", "baseline_scale", "rest_of_session_minutes",
                                  "rest_of_session_ratio"]).to_csv(out / "alarms.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    return summary


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("config")
    p = sub.add_parser("run")
    p.add_argument("--day", required=True)
    p.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "config":
        print(json.dumps(CONFIG, indent=2, sort_keys=True))
        print("sha256:", config_hash())
    else:
        s = run(a.day, a.dry_run)
        print(json.dumps({k: v for k, v in s.items() if k != "config"}, indent=2, default=str))


if __name__ == "__main__":
    main()
