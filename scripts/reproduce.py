"""One-command reproduction of the repository's results, with hash checks.

Usage:
    python scripts/reproduce.py                 # data check, tests, quick tier, comparison (about 15 minutes)
    python scripts/reproduce.py --tier full     # every experiment (several hours on one machine)
    python scripts/reproduce.py --steps data,tests
    python scripts/reproduce.py --steps manifest   # rewrite results/MANIFEST.csv from the files on disk
    python scripts/reproduce.py --tier full --steps run,compare --from E15   # resume a full run at E15
    python scripts/reproduce.py --tier full --steps run,compare --from E13 --to E14   # only the E13 and E14 entries

Steps (run in this order; pick a subset with --steps):

data      Check every file listed in data/MANIFEST.csv: present, same size, same SHA-256.
          Missing files are fetched with scripts/download_data.py when --download is given.
tests     Run the test suite (no data needed).
run       Run the experiments of the chosen tier. Each experiment script writes its CSV files
          to results/ (or alternatives/<name>/results/), overwriting the committed copies.
compare   Compare every result file against results/MANIFEST.csv: identical bytes, numerically
          equal within a tolerance (CSV only), equal except timing columns (CSV only; the columns
          listed in TIMING_COLUMNS hold run times, which no rerun reproduces), different, or
          missing. With git available the committed copy can also be restored with
          `git checkout -- results`.
report    Rebuild docs/experiment_report.md from docs/report_parts/.
manifest  Write results/MANIFEST.csv (path, bytes, SHA-256) for every result file on disk.

Exit status is 0 when every checked file matched, 1 otherwise.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable

# (label, command relative to ROOT, needs real data, tier). Tier "quick" runs in minutes on one core;
# "full" adds the long simulations and every real-data study. Approximate run times are in
# results/MANIFEST_RUNTIMES.csv after a full run.
EXPERIMENTS = [
    ("E1 reproduce tables", ["experiments/e01_reproduce_tables.py"], False, "full"),
    ("E2 truncation bias", ["experiments/e02_truncation_and_fair_comparison.py"], False, "full"),
    ("E3 online truncation", ["experiments/e03_online_truncation.py"], False, "quick"),
    ("E4 detection", ["experiments/e04_detection.py"], False, "full"),
    ("E4c null band", ["experiments/e04c_null_band.py"], False, "full"),
    ("E4d calibration length", ["experiments/e04d_calibration_length.py"], False, "full"),
    ("E4e review checks", ["experiments/e04e_review_checks.py"], False, "full"),
    ("E5 Bitstamp BTC", ["experiments/e05_real_btc.py"], True, "full"),
    ("E6a five coins tables", ["experiments/e06a_five_coins_tables.py"], True, "full"),
    ("E6b five coins monitor", ["experiments/e06b_five_coins_monitor.py"], True, "full"),
    ("E7 Nifty", ["experiments/e07_nifty.py"], True, "full"),
    ("E7c India VIX", ["experiments/e07c_india_vix.py"], True, "full"),
    ("E7d Nifty events", ["experiments/e07d_events.py"], True, "full"),
    ("E8 confidence-sequence baseline", ["experiments/e08_cs_baseline.py"], True, "full"),
    ("E9 stress tests", ["experiments/e09_stress_guarantee.py"], False, "full"),
    ("E10 forecast-relative", ["experiments/e10_forecast_relative.py"], True, "full"),
    ("E10d signal cap", ["experiments/e10_forecast_relative.py", "d"], False, "quick"),
    ("E11 delay theory", ["experiments/e11_delay_theory.py"], False, "full"),
    ("E11f decrease-side clip", ["experiments/e11f_split_c.py"], True, "full"),
    ("E12 crypto events", ["experiments/e12_crypto_events.py"], True, "full"),
    ("E13 offline size", ["experiments/e13_benchmarks.py", "offline-null"], False, "full"),
    ("E13 offline power", ["experiments/e13_benchmarks.py", "offline-change"], False, "full"),
    ("E13 online false alarms", ["experiments/e13_benchmarks.py", "online-null"], False, "full"),
    ("E13 online detection, part 0", ["experiments/e13_benchmarks.py", "online-change", "0"], False, "full"),
    ("E13 online detection, part 1", ["experiments/e13_benchmarks.py", "online-change", "1"], False, "full"),
    ("E13 real-data thresholds", ["experiments/e13_benchmarks.py", "real-calib"], False, "full"),
    ("E13 five coins, Jul to Dec 2025", ["experiments/e13_benchmarks.py", "real-crypto"], True, "full"),
    ("E13 five coins, Jan to Jun 2025 (held out)", ["experiments/e13_benchmarks.py", "real-crypto-h1-2025"], True, "full"),
    ("E13 Nifty", ["experiments/e13_benchmarks.py", "real-nifty"], True, "full"),
    ("E13 cost per step", ["experiments/e13_benchmarks.py", "cost"], False, "full"),
    # Review round 3 of E13, in the order of the script's docstring: these parts add rows to, or rewrite rows of, the
    # files written above (e13a_*, e13b_*, e13c_*_part0/1, e13d_real_*), so they must run after them; the figures last.
    ("E13 review round 3: VR scan with k = 250 and 500 offline", ["experiments/e13_benchmarks.py", "offline-long"], False, "full"),
    ("E13 review round 3: rolling VR scan with k = 250 and 500 online", ["experiments/e13_benchmarks.py", "online-long"], False, "full"),
    ("E13 review round 3: nominal rolling VR scan rows with the right series length", ["experiments/e13_benchmarks.py", "fix-nominal"], False, "full"),
    ("E13 review round 3: capped BOCPD against the exact recursion", ["experiments/e13_benchmarks.py", "cap-check"], False, "full"),
    ("E13 review round 3: longer competitor windows, Nifty", ["experiments/e13_benchmarks.py", "real-long", "nifty"], True, "full"),
    ("E13 review round 3: longer competitor windows, five coins, both half-years", ["experiments/e13_benchmarks.py", "real-long", "crypto"], True, "full"),
    ("E13 review round 3: pooled table with clustered tests", ["experiments/e13_benchmarks.py", "long-pooled"], False, "full"),
    ("E13 figures", ["experiments/e13_benchmarks.py", "figures"], True, "full"),
    # Parts a and b of E14 were committed one calibration length or one law at a time (e14a_validity_cal<n>.csv,
    # e14b_delay_<law>.csv); the seeds do not depend on the split, so these commands rewrite exactly those files.
    ("E14 running truncation, validity, calibration 390", ["experiments/e14_running_tau.py", "a", "--cal=390"], False, "full"),
    ("E14 running truncation, validity, calibration 1440", ["experiments/e14_running_tau.py", "a", "--cal=1440"], False, "full"),
    ("E14 running truncation, delay, Gaussian", ["experiments/e14_running_tau.py", "b", "--law=gaussian"], False, "full"),
    ("E14 running truncation, delay, Student t3", ["experiments/e14_running_tau.py", "b", "--law=student_t3"], False, "full"),
    ("E14 running truncation, delay, jumps", ["experiments/e14_running_tau.py", "b", "--law=jumps"], False, "full"),
    ("E14 running truncation, synthetic", ["experiments/e14_running_tau.py", "c", "e", "f", "h", "w"], False, "full"),
    ("E14 running truncation, real data and figures", ["experiments/e14_running_tau.py", "d", "g", "figures"], True, "full"),
    ("E15 panel, synthetic", ["experiments/e15_panel.py", "synthetic"], False, "full"),
    ("E15 panel, horizon", ["experiments/e15_panel.py", "horizon"], False, "full"),
    ("E15 panel, confidence-sequence baseline", ["experiments/e15_panel.py", "cs"], False, "full"),
    ("E15 panel, real data", ["experiments/e15_panel.py", "real"], True, "full"),
    ("E15 panel, figures", ["experiments/e15_panel.py", "figures"], True, "full"),
    ("E16 adversarial null search, inside the class", ["experiments/e16_adversarial_null.py", "a"], False, "full"),
    ("E16 adversarial null search, outside the class", ["experiments/e16_adversarial_null.py", "b"], False, "full"),
    ("E16 figures", ["experiments/e16_adversarial_null.py", "figures"], False, "full"),
    ("E17 multi-frequency, synthetic", ["experiments/e17_multifrequency.py", "baselines", "levels", "null", "change"], False, "full"),
    ("E17 multi-frequency, real data and figures", ["experiments/e17_multifrequency.py", "real", "ratematch", "figures"], True, "full"),
    ("E18 landmark events, pre-registered case studies", ["experiments/e18_landmark_events.py", "run", "figure", "context", "awake"], True, "full"),
    ("E21 benchmarks not calibrated on Gaussian data only, synthetic", ["experiments/e21_benchmarks.py", "calib"], False, "full"),
    ("E21 null", ["experiments/e21_benchmarks.py", "null"], False, "full"),
    ("E21 changes, part 0", ["experiments/e21_benchmarks.py", "change", "0"], False, "full"),
    ("E21 changes, part 1", ["experiments/e21_benchmarks.py", "change", "1"], False, "full"),
    ("E21 changes, merge", ["experiments/e21_benchmarks.py", "change-merge"], False, "full"),
    ("E21 on E13 paths, calibration", ["experiments/e21_benchmarks.py", "e13-calib"], False, "full"),
    ("E21 on E13 paths, null", ["experiments/e21_benchmarks.py", "e13-null"], False, "full"),
    ("E21 on E13 paths, changes part 0", ["experiments/e21_benchmarks.py", "e13-change", "0"], False, "full"),
    ("E21 on E13 paths, changes part 1", ["experiments/e21_benchmarks.py", "e13-change", "1"], False, "full"),
    ("E21 on E13 paths, merge", ["experiments/e21_benchmarks.py", "e13-change-merge"], False, "full"),
    ("E21 real data, calibration", ["experiments/e21_benchmarks.py", "real-calib"], True, "full"),
    ("E21 real data, BOCPD-t w = 250", ["experiments/e21_benchmarks.py", "real-bocpd", "250"], True, "full"),
    ("E21 real data, BOCPD-t w = 720", ["experiments/e21_benchmarks.py", "real-bocpd", "720"], True, "full"),
    ("E21 real data, BOCPD-t w = 1440", ["experiments/e21_benchmarks.py", "real-bocpd", "1440"], True, "full"),
    ("E21 real data, h2", ["experiments/e21_benchmarks.py", "real", "h2"], True, "full"),
    ("E21 real data, h1", ["experiments/e21_benchmarks.py", "real", "h1"], True, "full"),
    ("E21 real data, pooled", ["experiments/e21_benchmarks.py", "real-pooled"], True, "full"),
    ("E21 figures", ["experiments/e21_benchmarks.py", "figures"], True, "full"),
    ("E19 the guarantee on real data (pre-registered parts, then post hoc parts)", ["experiments/e19_real_null.py", "pools", "p1", "p2", "p3", "p4", "blind", "p1std", "p2boot", "figs"], True, "full"),
    ("E20 re-analysis of the real-data inference", ["experiments/e20_inference.py", "checks", "bootstrap", "holm", "events", "forecast", "classify"], True, "full"),
    ("E22 FOMC and CPI test with monitors that carry a proved per-day bound (pre-registered; checks E12 and E14 first)", ["experiments/e22_guaranteed_events.py", "run"], True, "full"),
    ("E23 power calculation of the pre-registration", ["experiments/e23_fresh_releases.py", "power"], False, "full"),
    ("E23 the same test on fresh releases (FOMC and CPI of 2026; pre-registered; gates against E12 and E22 first)", ["experiments/e23_fresh_releases.py", "run", "--procs=2"], True, "full"),
    ("E24 guard against blind segments: closed forms", ["experiments/e24_blind_guard.py", "theory"], False, "full"),
    ("E24 guard: the selection effect of naive retries", ["experiments/e24_blind_guard.py", "selection"], False, "full"),
    ("E24 guard, synthetic designs (800 paths, about an hour)", ["experiments/e24_blind_guard.py", "synthetic", "--procs=2", "--paths=800"], False, "full"),
    ("E24 guard on real data (post hoc)", ["experiments/e24_blind_guard.py", "real", "--procs=2"], True, "full"),
    ("E24 E19 blind diagnostic with the guard", ["experiments/e24_blind_guard.py", "e19blind", "--procs=2"], True, "full"),
    ("E24 pipeline gate on E23's period", ["experiments/e24_blind_guard.py", "future_gate", "--procs=2"], True, "full"),
    ("E24 forward pre-registration, power", ["experiments/e24_blind_guard.py", "power", "--procs=2"], False, "full"),
    ("E24 figures", ["experiments/e24_blind_guard.py", "figs"], False, "full"),
    ("E25 anytime-valid release test: exact identity and simulated validity", ["experiments/e25_sequential_events.py", "theory"], False, "full"),
    ("E25 planned strata", ["experiments/e25_sequential_events.py", "strata"], False, "full"),
    ("E25 causality and pipeline gates on known data", ["experiments/e25_sequential_events.py", "gates", "--procs=2"], True, "full"),
    ("E25 power", ["experiments/e25_sequential_events.py", "power"], False, "full"),
    ("E25 post hoc illustration on seen releases (reads the gates and the committed alarm files)", ["experiments/e25_sequential_events.py", "posthoc"], False, "full"),
    ("E25 rules as first committed (comparison): validity", ["experiments/e25_sequential_events.py", "theory", "--rules=original"], False, "full"),
    ("E25 rules as first committed (comparison): strata", ["experiments/e25_sequential_events.py", "strata", "--rules=original"], False, "full"),
    ("E25 rules as first committed (comparison): power", ["experiments/e25_sequential_events.py", "power", "--rules=original"], False, "full"),
    ("E25 rules as first committed (comparison): post hoc", ["experiments/e25_sequential_events.py", "posthoc", "--rules=original"], False, "full"),
    ("E25 figures", ["experiments/e25_sequential_events.py", "figs"], False, "full"),
    ("Replay page of two landmark days (rebuilds docs/demo/replay.html and checks it against E18)", ["scripts/build_replay.py"], True, "full"),
    ("TH4 robustness to drift in the conditional mean", ["experiments/th4_robustness.py", "sup", "constant", "switching", "garch", "plugin", "ticks", "counterexamples", "figures"], False, "full"),
    ("TH6 false alarms after the first day, per-day guarantees", ["experiments/th6_per_period.py", "laws", "hazard390", "hazard1440", "family", "adversary", "delay", "renewal", "tolerance", "figures"], False, "full"),
    ("TH5 limits for every valid procedure and the plug-in regime", ["experiments/th5_limits.py", "collapse", "floors", "asymptotic", "boundary", "plugin_bound", "plugin_sim", "holdbold", "figures"], False, "full"),
]

# Columns that hold measured run times (wall-clock or CPU time on the machine of the run), by result file. Compare
# leaves them out and reports such a file as "equal except timing columns" when every other column matches. The E9,
# E13 and E14 columns were checked against a full rerun (only these columns differed); the E22 to E24 "seconds"
# column is the wall time of one monitor run (experiments/e22_guaranteed_events.py, also used by E23 and E24).
TIMING_COLUMNS = {
    "results/e09d_profile.csv": ("seconds", "microseconds_per_unit"),
    "results/e13e_cost.csv": ("seconds", "microseconds_per_step"),
    "results/e14d_real_summary.csv": ("seconds", "us_per_step"),
    "results/e14e_cost.csv": ("seconds", "us_per_step"),
    "results/e14h_streaming.csv": ("us_per_step",),
    "results/e22_monitor_summary.csv": ("seconds",),
    "results/e23_monitor_summary.csv": ("seconds",),
    "results/e24_real_summary.csv": ("seconds",),
}

RESULT_DIRS = ["results"] + [str(p.relative_to(ROOT)) for p in sorted((ROOT / "alternatives").glob("*/results"))]
RESULT_SUFFIXES = {".csv", ".txt", ".md", ".json"}
MANIFEST = ROOT / "results" / "MANIFEST.csv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def result_files() -> list[Path]:
    out = []
    for d in RESULT_DIRS:
        for p in sorted((ROOT / d).rglob("*")):
            if p.is_file() and p.suffix in RESULT_SUFFIXES and p.name not in {"MANIFEST.csv", "MANIFEST_RUNTIMES.csv"}:
                out.append(p)
    return out


def step_data(download: bool) -> bool:
    import csv

    man = ROOT / "data" / "MANIFEST.csv"
    if not man.exists():
        print("data: no data/MANIFEST.csv; run scripts/download_data.py first")
        return False
    rows = list(csv.DictReader(open(man)))
    missing = [r for r in rows if not (ROOT / r["file"]).exists()]
    if missing and download:
        print(f"data: {len(missing)} files missing, downloading")
        subprocess.run([PY, str(ROOT / "scripts" / "download_data.py")], cwd=ROOT, check=False)
    ok = bad = miss = changing = 0
    for r in rows:
        p = ROOT / r["file"]
        if not p.exists():
            miss += 1
            continue
        if p.stat().st_size == int(r["bytes"]) and sha256(p) == r["sha256"]:
            ok += 1
        elif "latest" in p.name:
            changing += 1  # the Bitstamp file is updated daily upstream; experiments cut it to fixed dates
        else:
            bad += 1
            print(f"  hash mismatch: {r['file']}")
    print(f"data: {ok} match, {bad} differ, {miss} missing, {changing} daily-updated file(s) changed upstream "
          f"(of {len(rows)} in data/MANIFEST.csv)")
    return bad == 0 and miss == 0


def step_tests() -> bool:
    res = subprocess.run([PY, "-m", "pytest", "-q"], cwd=ROOT)
    return res.returncode == 0


def step_run(tier: str, have_data: bool, start: str | None = None, only: str | None = None,
             stop: str | None = None) -> bool:
    import csv

    rows, ok = [], True
    started = start is None
    last = max((i for i, e in enumerate(EXPERIMENTS) if e[0].startswith(stop)), default=-1) if stop else len(EXPERIMENTS)
    if stop and last < 0:
        print(f"run: no experiment label begins with {stop!r}")
        return False
    for i, (label, cmd, needs_data, t) in enumerate(EXPERIMENTS):
        if i > last:
            break
        if tier == "quick" and t != "quick":
            continue
        started = started or label.startswith(start)
        if not started or (only is not None and only not in label):
            continue
        if needs_data and not have_data:
            print(f"run: skip {label} (needs data)")
            continue
        script = ROOT / cmd[0]
        if not script.exists():
            print(f"run: skip {label} ({cmd[0]} not found)")
            continue
        print(f"run: {label}: {' '.join(cmd)}", flush=True)
        t0 = time.time()
        res = subprocess.run([PY, *cmd], cwd=ROOT)
        dt = time.time() - t0
        rows.append(dict(experiment=label, command=" ".join(cmd), seconds=round(dt), status=res.returncode))
        ok &= res.returncode == 0
    if tier == "full" and rows and start is None and only is None and stop is None:
        with open(ROOT / "results" / "MANIFEST_RUNTIMES.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return ok


def _numeric_equal(a: Path, b_text: str, rtol: float, ignore: tuple = ()) -> bool:
    """Same columns and shape, numeric columns equal within rtol, other columns equal as text; columns in ``ignore``
    are left out of the comparison (they must still exist in both files)."""
    import io

    import numpy as np
    import pandas as pd

    try:
        x = pd.read_csv(a)
        y = pd.read_csv(io.StringIO(b_text))
    except Exception:
        return False
    if list(x.columns) != list(y.columns) or x.shape != y.shape:
        return False
    for c in x.columns:
        if c in ignore:
            continue
        xa, ya = x[c], y[c]
        if pd.api.types.is_numeric_dtype(xa) and pd.api.types.is_numeric_dtype(ya):
            xv, yv = xa.to_numpy(float), ya.to_numpy(float)
            same_nan = np.isnan(xv) == np.isnan(yv)
            fin = ~np.isnan(xv) & ~np.isnan(yv)
            if not same_nan.all() or not np.allclose(xv[fin], yv[fin], rtol=rtol, atol=1e-12):
                return False
        elif not ((xa.astype(str) == ya.astype(str)) | (xa.isna() & ya.isna())).all():
            return False  # (pandas reads the text "null" as missing; missing equals missing)
    return True


def step_compare(rtol: float) -> bool:
    import csv

    if not MANIFEST.exists():
        print("compare: no results/MANIFEST.csv")
        return False
    ref = {r["file"]: r for r in csv.DictReader(open(MANIFEST))}
    counts = dict(identical=0, close=0, timing=0, different=0, missing=0)
    for f, r in ref.items():
        p = ROOT / f
        if not p.exists():
            counts["missing"] += 1
            print(f"  missing: {f}")
            continue
        if sha256(p) == r["sha256"]:
            counts["identical"] += 1
            continue
        committed = None
        try:
            committed = subprocess.run(["git", "show", f"HEAD:{f}"], cwd=ROOT, capture_output=True, text=True,
                                       check=True).stdout
        except Exception:
            pass
        if committed is not None and p.suffix == ".csv" and _numeric_equal(p, committed, rtol):
            counts["close"] += 1
        elif (committed is not None and p.suffix == ".csv" and f in TIMING_COLUMNS
              and _numeric_equal(p, committed, rtol, ignore=TIMING_COLUMNS[f])):
            counts["timing"] += 1
            print(f"  equal except timing columns ({', '.join(TIMING_COLUMNS[f])}): {f}")
        else:
            counts["different"] += 1
            print(f"  different: {f}")
    new = [p for p in result_files() if str(p.relative_to(ROOT)) not in ref]
    print(f"compare: {counts['identical']} identical, {counts['close']} equal within rtol {rtol:g}, "
          f"{counts['timing']} equal except timing columns, {counts['different']} different, "
          f"{counts['missing']} missing; {len(new)} files not in the manifest")
    return counts["different"] == 0 and counts["missing"] == 0


def step_report() -> bool:
    return subprocess.run([PY, str(ROOT / "scripts" / "build_report.py")], cwd=ROOT).returncode == 0


def step_manifest() -> bool:
    import csv

    files = result_files()
    with open(MANIFEST, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "bytes", "sha256"])
        for p in files:
            w.writerow([str(p.relative_to(ROOT)), p.stat().st_size, sha256(p)])
    print(f"manifest: {len(files)} result files written to {MANIFEST.relative_to(ROOT)}")
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", choices=["quick", "full"], default="quick")
    ap.add_argument("--steps", default="data,tests,run,compare")
    ap.add_argument("--download", action="store_true", help="fetch missing data files")
    ap.add_argument("--from", dest="start", default=None,
                    help="run step: start at the first experiment whose label begins with this text (to resume)")
    ap.add_argument("--only", default=None, help="run step: run only experiments whose label contains this text")
    ap.add_argument("--to", dest="stop", default=None,
                    help="run step: stop after the last experiment whose label begins with this text")
    ap.add_argument("--rtol", type=float, default=1e-9, help="relative tolerance for numeric comparison")
    a = ap.parse_args(argv)
    steps = [s.strip() for s in a.steps.split(",") if s.strip()]
    status = {}
    have_data = True
    for s in steps:
        print(f"== {s}", flush=True)
        if s == "data":
            status[s] = have_data = step_data(a.download)
        elif s == "tests":
            status[s] = step_tests()
        elif s == "run":
            status[s] = step_run(a.tier, have_data, a.start, a.only, a.stop)
        elif s == "compare":
            status[s] = step_compare(a.rtol)
        elif s == "report":
            status[s] = step_report()
        elif s == "manifest":
            status[s] = step_manifest()
        else:
            raise SystemExit(f"unknown step {s}")
    print("summary:", ", ".join(f"{k} {'ok' if v else 'FAILED'}" for k, v in status.items()))
    return 0 if all(status.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
