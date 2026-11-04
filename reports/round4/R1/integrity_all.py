"""R1 Part A. Integrity of every real data file: the 140 files of the committed data/MANIFEST.csv, the daily-updated
Bitstamp file (cut to the committed last bar) and the 2026 Binance months that E23 uses (checked against
results/e23_data_files.csv). Extends reports/stress/b5_data_integrity.py (same checks, every file, plus a look at the
10 largest absolute returns of each file).

Usage: python integrity_all.py <out_dir>
Writes <out_dir>/integrity_files.csv (one row per file), <out_dir>/integrity_top10.csv (10 largest |returns| per
file with a verdict) and <out_dir>/integrity_nifty_sessions.csv (per session day: bars present, missing minutes).
"""
import hashlib
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/user/anytime-valid-volatility-monitor")
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
BITSTAMP_COMMITTED_ROWS = 956_327
BITSTAMP_COMMITTED_LAST = "2026-11-02"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def committed_manifest() -> pd.DataFrame:
    txt = subprocess.run(["git", "show", "HEAD:data/MANIFEST.csv"], cwd=ROOT, capture_output=True, text=True,
                         check=True).stdout
    from io import StringIO
    m = pd.read_csv(StringIO(txt))
    m["source"] = "committed data/MANIFEST.csv"
    return m


def e23_manifest() -> pd.DataFrame:
    m = pd.read_csv(ROOT / "results" / "e23_data_files.csv")
    m["source"] = "results/e23_data_files.csv (E23 fresh months)"
    return m


def read_any(path: Path) -> tuple[pd.DataFrame, str]:
    """Return a frame with columns ts (UTC, ns), open, high, low, close, volume (NaN if absent) and the kind."""
    raw = pd.read_csv(path)
    name = path.name
    if "binance" in str(path):
        kind = "binance"
        ts = pd.to_datetime(raw["timestamp"].to_numpy(), unit="ms", utc=True)
        df = pd.DataFrame(dict(ts=ts, open=raw.open, high=raw.high, low=raw.low, close=raw.close, volume=raw.volume))
    elif "bitstamp" in name:
        kind = "bitstamp"
        ts = pd.to_datetime(raw["timestamp"].to_numpy(), unit="s", utc=True)
        df = pd.DataFrame(dict(ts=ts, open=raw.open, high=raw.high, low=raw.low, close=raw.close, volume=raw.volume))
    else:
        kind = "nifty"
        ts = pd.to_datetime(raw["Epoch"].to_numpy(), unit="s", utc=True)
        vol = raw["Volume"] if "Volume" in raw.columns else pd.Series(np.nan, index=raw.index)
        df = pd.DataFrame(dict(ts=ts, open=raw.Open, high=raw.High, low=raw.Low, close=raw.Close, volume=vol))
        df.attrs["has_volume"] = "Volume" in raw.columns
        # the text timestamp should be the epoch in IST
        txt = pd.to_datetime(raw["Timestamp"]).dt.tz_localize("Asia/Kolkata")
        df.attrs["text_ts_mismatch"] = int((txt.dt.tz_convert("UTC").to_numpy() != ts.to_numpy()).sum())
    df.attrs["nrows_raw"] = len(raw)
    return df, kind


def check_file(path: Path, exp_rows, exp_sha, tag: str) -> tuple[dict, list[dict], list[dict]]:
    t0 = time.time()
    row = dict(file=str(path.relative_to(ROOT)), group=tag)
    sha = sha256(path)
    df, kind = read_any(path)
    n = len(df)
    row["kind"] = kind
    row["rows"] = n
    row["rows_manifest"] = exp_rows
    row["sha_match"] = (sha == exp_sha) if exp_sha is not None else None
    row["rows_match"] = (n == exp_rows) if exp_rows is not None else None
    ts = (df.ts.astype("datetime64[ns, UTC]").astype("int64") // 1_000_000_000).to_numpy()
    row["duplicate_ts"] = int(pd.Series(ts).duplicated().sum())
    d = np.diff(ts)
    row["strictly_increasing"] = bool(np.all(d > 0))
    row["whole_minutes"] = bool(np.all(ts % 60 == 0))
    row["first_bar"] = str(df.ts.iloc[0])[:16]
    row["last_bar"] = str(df.ts.iloc[-1])[:16]
    sessions = []
    if kind == "binance":
        per = pd.Period(df.ts.iloc[0].strftime("%Y-%m"), "M")
        full = pd.date_range(per.start_time.tz_localize("UTC"), (per + 1).start_time.tz_localize("UTC"), freq="1min",
                             inclusive="left")
        miss = full.difference(pd.DatetimeIndex(df.ts.drop_duplicates()))
        row["expected_minutes"] = len(full)
    elif kind == "bitstamp":
        full = pd.date_range(df.ts.iloc[0], df.ts.iloc[-1], freq="1min")
        miss = full.difference(pd.DatetimeIndex(df.ts.drop_duplicates()))
        row["expected_minutes"] = len(full)
    else:
        # regular NSE session 09:15 to 15:29 IST, 375 one-minute bars per day; count inside the days present
        ist = df.ts.dt.tz_convert("Asia/Kolkata")
        day = ist.dt.date
        hm = ist.dt.hour * 60 + ist.dt.minute
        inside = (hm >= 9 * 60 + 15) & (hm <= 15 * 60 + 29)
        row["bars_outside_0915_1529"] = int((~inside).sum())
        miss_list = []
        for dday, g in pd.Series(hm[inside].to_numpy(), index=day[inside].to_numpy()).groupby(level=0):
            present = set(g.to_numpy())
            m = 375 - len(present)
            sessions.append(dict(file=path.name, day=str(dday), bars=len(present), missing=m,
                                 first=f"{min(present)//60:02d}:{min(present)%60:02d}",
                                 last=f"{max(present)//60:02d}:{max(present)%60:02d}"))
            miss_list.append(m)
        miss = np.array(miss_list)
        row["sessions"] = len(sessions)
        row["short_sessions"] = int((miss > 0).sum())
        row["sessions_missing_over_60min"] = int((miss > 60).sum())
        row["expected_minutes"] = 375 * len(sessions)
        row["text_timestamp_mismatch"] = df.attrs.get("text_ts_mismatch")
        row["has_volume_column"] = df.attrs.get("has_volume")
    if kind in ("binance", "bitstamp"):
        row["missing_minutes"] = int(len(miss))
        if len(miss):
            gaps = pd.Series(miss).diff().ne(pd.Timedelta(minutes=1)).cumsum()
            sizes = pd.Series(miss).groupby(gaps.values).size()
            row["gaps"] = int(sizes.size)
            row["longest_gap_min"] = int(sizes.max())
            row["gaps_over_60min"] = int((sizes > 60).sum())
        else:
            row["gaps"] = 0
            row["longest_gap_min"] = 0
            row["gaps_over_60min"] = 0
    else:
        row["missing_minutes"] = int(miss.sum())
    for c in ["open", "high", "low", "close"]:
        row[f"bad_{c}"] = int(((df[c] <= 0) | df[c].isna()).sum())
    row["ohlc_inconsistent"] = int(((df.high < df.low) | (df.close > df.high) | (df.close < df.low)
                                    | (df.open > df.high) | (df.open < df.low)).sum())
    row["zero_volume_bars"] = int((df.volume <= 0).sum()) if df.volume.notna().any() else None
    row["flat_bars"] = int(((df.open == df.high) & (df.high == df.low) & (df.low == df.close)).sum())
    close = df.close.to_numpy(float)
    lr_all = np.diff(np.log(np.where(close > 0, close, np.nan)))
    consecutive = d == 60
    lr = lr_all[consecutive]
    lr = lr[np.isfinite(lr)]
    sd = np.median(np.abs(lr)) * 1.4826
    if lr.size == 0:
        lr = np.array([np.nan])
    row["n_returns"] = int(np.isfinite(lr).sum())
    row["robust_sd"] = float(sd)
    row["zero_return_share"] = float((lr == 0).mean())
    row["beyond_20sd"] = int((np.abs(lr) > 20 * sd).sum())
    row["max_abs_return"] = float(np.abs(lr).max())
    row["max_return_over_sd"] = float(np.abs(lr).max() / sd)
    wick = ((df.high - df.low) / df.close).to_numpy(float)
    row["max_range_over_close"] = float(np.nanmax(wick))
    # 10 largest |returns| over consecutive minutes: real move or bad print?
    idx = np.where(consecutive & np.isfinite(lr_all))[0]  # return i is from bar i to bar i+1
    order = idx[np.argsort(-np.abs(lr_all[idx]))[:10]]
    top = []
    n_suspect = 0
    for i in order:
        r = lr_all[i]
        j = i + 1
        hi = max(df.high.iloc[i], df.high.iloc[j])
        lo = min(df.low.iloc[i], df.low.iloc[j])
        move = abs(close[j] - close[i])
        covered = (hi - lo) >= move * (1 - 1e-9)  # the two bars' high-low span covers the move
        nxt = lr_all[j] if j < len(lr_all) and d[j] == 60 else np.nan
        reversal = bool(np.isfinite(nxt) and np.sign(nxt) == -np.sign(r) and abs(nxt) >= 0.5 * abs(r))
        # was the move already visible before (close of bar i within 10% of its own open)? not used for verdict
        bar_range_j = (df.high.iloc[j] - df.low.iloc[j]) / close[j]
        if not covered:
            verdict = "SUSPECT: move larger than the high-low span of the two bars"
        elif reversal:
            verdict = "SUSPECT: at least half reversed within one minute"
        else:
            verdict = "real move (range covers it, no instant reversal)"
        n_suspect += verdict.startswith("SUSPECT")
        top.append(dict(file=path.name, group=tag, time=str(df.ts.iloc[j])[:19], logret=float(r),
                        ret_over_sd=float(abs(r) / sd), close_before=float(close[i]), close_after=float(close[j]),
                        bar_range_over_close=float(bar_range_j), next_logret=float(nxt) if np.isfinite(nxt) else np.nan, reversal_fraction=float(-nxt / r) if np.isfinite(nxt) else np.nan,
                        volume=float(df.volume.iloc[j]) if df.volume.notna().any() else np.nan, verdict=verdict))
    row["top10_suspect"] = n_suspect
    row["seconds"] = round(time.time() - t0, 1)
    return row, top, sessions


def main():
    man = committed_manifest()
    e23 = e23_manifest()
    done = set(man.file)
    e23 = e23[~e23.file.isin(done)]
    rows, tops, sess = [], [], []
    for frame, tag in ((man, "manifest"), (e23, "e23_2026")):
        for r in frame.itertuples():
            p = ROOT / r.file
            if not p.exists():
                rows.append(dict(file=r.file, group=tag, kind="missing", sha_match=False, rows_match=False))
                print("MISSING", r.file, flush=True)
                continue
            exp_rows, exp_sha = int(r.rows), r.sha256
            if "latest" in p.name:
                # daily-updated file: compare the rows up to the committed last bar instead of the hash
                raw = pd.read_csv(p)
                last = pd.to_datetime(raw.timestamp, unit="s", utc=True)
                keep = last < pd.Timestamp(BITSTAMP_COMMITTED_LAST, tz="UTC") + pd.Timedelta(days=1)
                cut = OUT / "bitstamp_cut_to_committed_last_bar.csv"
                raw[keep].to_csv(cut, index=False)
                row, top, s = check_file(cut, exp_rows, exp_sha, tag)
                row["file"] = r.file + " (rows with timestamp <= 2026-11-02 23:59 UTC)"
                row["rows_on_disk_today"] = len(raw)
                row["sha_on_disk_today_matches"] = sha256(p) == exp_sha
                row["note"] = ("daily-updated upstream; hash mismatch expected; row count compared on the rows up to "
                               "the committed last bar")
                cut.unlink()
            else:
                row, top, s = check_file(p, exp_rows, exp_sha, tag)
            rows.append(row)
            tops += top
            sess += s
            print(f"{row['file']}: sha {row.get('sha_match')} rows {row.get('rows_match')} dup {row.get('duplicate_ts')} "
                  f"miss {row.get('missing_minutes')} >20sd {row.get('beyond_20sd')} suspect {row.get('top10_suspect')} "
                  f"({row.get('seconds')}s)", flush=True)
            pd.DataFrame(rows).to_csv(OUT / "integrity_files.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "integrity_files.csv", index=False)
    pd.DataFrame(tops).to_csv(OUT / "integrity_top10.csv", index=False)
    pd.DataFrame(sess).to_csv(OUT / "integrity_nifty_sessions.csv", index=False)
    print("INTEGRITY_DONE")


if __name__ == "__main__":
    main()
