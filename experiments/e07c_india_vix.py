"""E7c. External check of the Nifty monitor with India VIX (NSE implied volatility index).

REAL MARKET DATA. Source: NSE daily index close files
https://archives.nseindia.com/content/indices/ind_close_all_DDMMYYYY.csv
(one file per trading day; the file has an "India VIX" row with open, high,
low and close, and a "Nifty 50" row used here to cross-check the minute data).
India VIX is computed by NSE from Nifty option prices, so it measures the volatility
the options market expects over the next 30 days. It is a forward-looking quantity
that the monitor never sees.

Research question
  Q7c. Are up alarms of the streaming monitor followed by India VIX rises, and
       down alarms by falls, more often than random sessions?

Method (session level)
  For each session d, VIX change on the alarm day = log VIX close(d) - log VIX close(d-1)
  and on the next day = log VIX close(d+1) - log VIX close(d), with d-1 and d+1 the
  neighbouring sessions that have a VIX file. A session counts as an "up alarm session"
  if it has at least one up alarm (same for down). The control is all sessions of the
  same period. Alarms come from results/e07b_alarms.csv (run e07_nifty.py q7b first).

Usage
  python -I experiments/e07c_india_vix.py download   # fetch and cache the daily files (polite, 1 s pacing)
  python -I experiments/e07c_india_vix.py            # analysis
"""

import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/india_vix"
OUT = ROOT / "results"
URL = "https://archives.nseindia.com/content/indices/ind_close_all_{d}.csv"
DEV_END = pd.Timestamp("2022-01-01")


def download(pause: float = 1.0):
    """Fetch one file per Nifty session date (cached; 404 means no file for that day)."""
    from avmonitor.data import load_nifty

    RAW.mkdir(parents=True, exist_ok=True)
    df = load_nifty("2017-01-01", "2026-11-01", download=False)
    days = sorted(pd.DatetimeIndex(df["session"].unique()).strftime("%d%m%Y"))
    log = []
    for i, d in enumerate(days):
        dest = RAW / f"ind_close_all_{d}.csv"
        if dest.exists() and dest.stat().st_size > 0:
            log.append(dict(day=d, status="cached", bytes=dest.stat().st_size))
            continue
        status, size = "error", 0
        for attempt in range(3):
            try:
                with urllib.request.urlopen(URL.format(d=d), timeout=60) as resp:
                    body = resp.read()
                tmp = dest.with_suffix(".part")
                tmp.write_bytes(body)
                tmp.replace(dest)
                status, size = "ok", len(body)
                break
            except urllib.error.HTTPError as e:
                status = f"http {e.code}"
                if e.code == 404:
                    break
            except Exception as e:  # network hiccup: back off and retry
                status = f"error {type(e).__name__}"
            time.sleep(pause * 2 ** (attempt + 1))
        log.append(dict(day=d, status=status, bytes=size))
        if i % 100 == 0:
            print(i, d, status, flush=True)
        time.sleep(pause)
    pd.DataFrame(log).to_csv(RAW / "download_log.csv", index=False)
    print(pd.DataFrame(log).status.value_counts())


def parse() -> pd.DataFrame:
    """One row per file: India VIX OHLC and Nifty 50 OHLC (as published by NSE)."""
    rows, bad = [], []
    for p in sorted(RAW.glob("ind_close_all_*.csv")):
        try:
            f = pd.read_csv(p, dtype=str)
            f.columns = [c.strip() for c in f.columns]
            name = f["Index Name"].str.strip().str.lower()
            vix, nifty = f[name == "india vix"], f[name == "nifty 50"]
            if len(vix) != 1 or len(nifty) != 1:
                bad.append(p.name)
                continue
            num = lambda r, c: float(r[c].iloc[0])  # noqa: E731
            rows.append(dict(date=pd.to_datetime(p.stem[-8:], format="%d%m%Y"),
                             file_date=pd.to_datetime(vix["Index Date"].iloc[0], dayfirst=True),
                             vix_open=num(vix, "Open Index Value"), vix_high=num(vix, "High Index Value"),
                             vix_low=num(vix, "Low Index Value"), vix_close=num(vix, "Closing Index Value"),
                             nifty_open=num(nifty, "Open Index Value"), nifty_high=num(nifty, "High Index Value"),
                             nifty_low=num(nifty, "Low Index Value"), nifty_close=num(nifty, "Closing Index Value")))
        except Exception:
            bad.append(p.name)
    d = pd.DataFrame(rows).sort_values("date").reset_index(drop=True)
    d.attrs["unparsed"] = bad
    return d


def wilson(k, n, z=1.959964):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return (c - h, c + h)


def crosscheck(v: pd.DataFrame) -> pd.DataFrame:
    """Compare NSE's official Nifty 50 daily values with the minute data."""
    from avmonitor.data import load_nifty

    m = load_nifty("2017-01-01", "2026-11-01", download=False)
    g = m.groupby("session").agg(close=("close", "last"), high=("high", "max"), low=("low", "min"), bars=("close", "size"))
    g.index = g.index.tz_localize(None)
    j = v.set_index("date").join(g, how="inner")
    out = []
    for per, mask in (("all", np.ones(len(j), bool)), ("before Aug 2026", j.index < "2026-08-01"),
                      ("Aug to Oct 2026", j.index >= "2026-08-01")):
        x = j[mask]
        for col in ("close", "high", "low"):
            rel = (x[col] / x[f"nifty_{col}"] - 1).abs()
            out.append(dict(period=per, field=col, days=len(x), median_abs_rel_diff=float(rel.median()),
                            share_within_0_01pct=float((rel <= 1e-4).mean()),
                            share_within_0_1pct=float((rel <= 1e-3).mean()), max_abs_rel_diff=float(rel.max())))
    return pd.DataFrame(out)


def analyse(v: pd.DataFrame):
    al = pd.read_csv(OUT / "e07b_alarms.csv", parse_dates=["session"])
    sess = pd.read_csv(OUT / "e07b_session_vol.csv", parse_dates=["session"])["session"]
    v = v.set_index("date")
    lv = np.log(v["vix_close"])
    d_same = lv.diff()  # previous VIX day to alarm day
    d_next = lv.shift(-1) - lv  # alarm day to next VIX day
    s = pd.DataFrame({"session": sess})
    s["same"] = s["session"].map(d_same)
    s["next"] = s["session"].map(d_next)
    s["period"] = np.where(s["session"] < DEV_END, "dev", "held-out")
    rows = []
    for (det, cal, band), a in al.groupby(["detector", "cal_len", "band"]):
        for per in ("dev", "held-out", "all"):
            base = s if per == "all" else s[s["period"] == per]
            for d, sign in (("up", 1), ("down", -1)):
                ss = set(a.loc[a["direction"] == d, "session"])
                hit = base[base["session"].isin(ss)]
                rest = base[~base["session"].isin(ss)]
                for horizon in ("same", "next"):
                    x = hit[horizon].dropna()
                    c = base[horizon].dropna()
                    r = rest[horizon].dropna()
                    k = int((sign * x > 0).sum())
                    kc = int((sign * c > 0).sum())
                    kr = int((sign * r > 0).sum())
                    lo, hi = wilson(k, len(x))
                    clo, chi = wilson(kc, len(c))
                    _, p = stats.fisher_exact([[k, len(x) - k], [kr, len(r) - kr]]) if len(x) else (np.nan, np.nan)
                    rows.append(dict(detector=det, cal_len=cal, band=band, period=per, direction=d, horizon=horizon,
                                     alarm_sessions=len(x), share_vix_moved_same_way=k / len(x) if len(x) else np.nan,
                                     lo=lo, hi=hi, control_sessions=len(c), control_share=kc / len(c),
                                     control_lo=clo, control_hi=chi, fisher_p_vs_other_sessions=p,
                                     median_log_vix_change=float(x.median()) if len(x) else np.nan,
                                     control_median_log_vix_change=float(c.median())))
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "e07c_vix_alarm_check.csv", index=False)
    return res


def main():
    v = parse()
    print(f"REAL DATA: {len(v)} NSE daily files parsed, {len(v.attrs['unparsed'])} without a VIX and Nifty row:",
          v.attrs["unparsed"][:10])
    print("file date equals requested date:", bool((v["date"] == v["file_date"]).all()))
    v.to_csv(RAW / "india_vix_daily.csv", index=False)
    log = RAW / "download_log.csv"
    if log.exists():
        print(pd.read_csv(log)["status"].value_counts().to_dict())
    cc = crosscheck(v)
    cc.to_csv(OUT / "e07c_nse_crosscheck.csv", index=False)
    print(cc.to_string(index=False))
    desc = v["vix_close"].describe()
    print("India VIX close:", desc.round(2).to_dict())
    res = analyse(v)
    show = res[(res.period != "all")]
    show = show[((show.detector == "robust e-detector") & (show.cal_len == 750) & (show.band == 0.25))
                | ((show.detector == "Gaussian LR") & (show.cal_len == 750))]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(show.round(3).to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "download":
        download()
    else:
        main()
