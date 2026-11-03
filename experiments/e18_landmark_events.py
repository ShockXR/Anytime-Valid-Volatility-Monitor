"""E18. Landmark market shocks as pre-registered case studies (tracker B85).

REAL MARKET DATA. Binance spot one-minute bars (BTC, ETH, XRP, SOL, DOGE) and Nifty 50 one-minute bars, as in E6, E7
and E13.

Pre-registration: docs/preregistration_e18_events.md, committed on its own (commit b8ee0f7, 2026-11-03 00:37 UTC)
before this script existed and before any return in the event windows was loaded for this study. Everything below
follows it; nothing was tuned on the event windows.

Question. When a large, well-remembered market shock happens, how quickly do the e-detector and its published
competitors (E13) raise an up alarm, were they watching at that moment, and how many alarms did they raise in the week
before? Four events, up to six assets each: case studies, not evidence. No pooled test.

Events and reference times T (verified by web search before any returns were looked at; sources in the
pre-registration and in results/e18_events.csv):
  N1  4 Jun 2024, Indian general election results: Nifty open 09:15 IST (03:45 UTC).
  C1  5 Aug 2024, global sell-off: crypto 00:00 UTC (Tokyo open, 09:00 JST); N2 Nifty open 09:15 IST.
  C2  7 Apr 2025, tariff sell-off: crypto 00:00 UTC; N3 Nifty open 09:15 IST.
  C3  10 Oct 2025, crypto liquidation cascade: 20:00 UTC (US close), marked: the 100% tariff post came after the close
      and no later than 20:57 UTC, so delays are also given from 20:57 UTC.

Streams (as a user would have run them). Crypto: three half-years (H2 2024 for C1, H1 2025 for C2, H2 2025 for C3),
returns divided by the predictable intraday profile from data starting 1 Jan 2024; H1 2025 and H2 2025 are E13's two
samples. Nifty: E13's stream (2024 to 2025, E7's session profile).

Methods. e-detector with the frozen E6 / E7 settings (A = 28,800 or 7,500, c = 3, calibration 1440 or 750, band 0.25),
streaming from the stream start. E13's competitors with the first-run windows and the round-3 windows: BOCPD and BOCPD
on truncated returns (w = 250, 720, 1440 on crypto; 250, 375 on Nifty) and the truncated Gaussian-normalised
two-window variance-ratio scan (k = 120, 480, 1440; 120, 375), inspired by Bibinger, Jirak and Vetter (2017); the paper
could not be obtained, so this is not their test (earlier called the "BJV-type" test; the code and the CSV files keep
E13's internal label "BJV rolling"). Each competitor's up and down thresholds are matched to the e-detector's
up and down alarm counts on the same stream, counting only alarms outside the exclusion windows [T - 8 d, T + 2 d).

Outcomes per event, asset and method (primary): state at T (monitoring / calibrating, watching / paused), delay from T
to the first up alarm within 24 hours, alarms in the 7 days before T (with the number expected at the method's average
rate), alarms in the 24 hours after T. Secondary: S1 fresh start at T, S2 the e-detector streamed from another start
date (crypto from 21 Jan 2024 as in E12; Nifty from 2017 as in E7), S3 thresholds matched on the whole sample (E13's
rule).

Usage
  python experiments/e18_landmark_events.py run      # all streams and methods (about 40 minutes with 2 processes)
  python experiments/e18_landmark_events.py figure   # results/figures/e18_events.png from the saved CSV files
  python experiments/e18_landmark_events.py context  # results/e18_context.csv: descriptive return sizes around the
                                                     # e-detector's alarms on the crypto events (NOT pre-registered)
  python experiments/e18_landmark_events.py awake    # results/e18_awake.csv: E14's monitors with no blind period on the
                                                     # same events (exploratory, NOT pre-registered, post hoc)
  Optional: E18_CACHE=<folder> keeps the BOCPD arrays and finished series between attempts (not result files).

Outputs: results/e18_events.csv, e18_outcomes.csv, e18_alarms.csv, e18_thresholds.csv, e18_rates.csv,
results/figures/e18_events.png; results/e18_context.csv (descriptive, not pre-registered); results/e18_awake.csv
(exploratory, not pre-registered).
"""

from __future__ import annotations

import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
sys.path.insert(0, str(Path(__file__).resolve().parent))

import e13_benchmarks as E13  # noqa: E402  (real-data competitor code of E13, reused unchanged)
from e06b_five_coins_monitor import COINS, prepare  # noqa: E402
from e07_nifty import deseasonalise, load_clean  # noqa: E402

from avmonitor import competitors as C  # noqa: E402
from avmonitor.data import intraday_scale, load_binance  # noqa: E402
from avmonitor.monitor import VolatilityMonitor  # noqa: E402

# ----------------------------------------------------------------------------
# Fixed in the pre-registration
# ----------------------------------------------------------------------------
BAND = 0.25
H24 = pd.Timedelta(hours=24)
D7 = pd.Timedelta(days=7)
EXCL_BEFORE, EXCL_AFTER = pd.Timedelta(days=8), pd.Timedelta(days=2)
MARKETS = {
    "crypto": dict(day=1440, A=28800.0, cal=1440, ws=(250, 720, 1440), ks=(120, 480, 1440)),
    "nifty": dict(day=375, A=7500.0, cal=750, ws=(250, 375), ks=(120, 375)),
}

SRC_N1 = ("Election Commission via Deccan Herald (PTI): counting of votes from 8 am on June 4 (postal ballots first, "
          "machines from 08:30 IST), https://www.deccanherald.com/elections/india/lok-sabha-elections-2024-counting-of-"
          "votes-from-8-am-on-june-4-says-ec-3047984 | Business Standard live blog 4 Jun 2024: Sensex crashes 4,390 pts, "
          "Nifty gives up 22k on sub-par BJP tally, https://www.business-standard.com/markets/news/stock-market-live-"
          "updates-june-4-sensex-nifty-lok-sabha-election-results-2024-psu-stocks-nse-124060400081_1.html")
SRC_AUG = ("Reuters via Kontan: Nikkei fell 12.4% on Monday 5 Aug 2024, biggest rout since 1987, after weak US jobs data "
           "on Friday 2 Aug, https://english.kontan.co.id/news/japans-nikkei-sees-biggest-rout-since-1987-black-monday | "
           "JPX: Tokyo trading hours 9:00 to 11:30 and 12:30 to 15:00 (15:30 from 5 Nov 2024), https://www.jpx.co.jp/"
           "english/corporate/news/news-releases/1030/uorii50000002f2a-att/pressrelease_extension_of_trading_hours_en.pdf"
           " | Business Standard: indices witness a gap-down opening mirroring global markets, https://www.business-"
           "standard.com/amp/markets/capital-market-news/indices-witness-a-gap-down-opening-mirroring-global-mkts-"
           "124080500173_1.html | Benzinga: bitcoin's low came during Asian trading hours, https://www.benzinga.com/"
           "markets/cryptocurrency/24/08/40153940/bitcoin-crashes-below-50k-briefly-crypto-market-records-over-1b-in-"
           "liquidations-in-24-hour")
SRC_APR = ("Bloomberg: Bitcoin plunges into Asia open as crypto feels trade war pain (7 Apr 2025), https://www.bloomberg."
           "com/news/articles/2025-04-07/bitcoin-btc-plunges-into-asia-open-as-crypto-feels-trade-war-pain | CNBC: bitcoin "
           "drops Sunday evening as cryptocurrencies join global market rout (6 Apr 2025), https://www.cnbc.com/2025/04/"
           "06/bitcoin-drops-sunday-evening-as-cryptocurrencies-join-global-market-rout.html | Business Standard: Nikkei "
           "hits lower limit (7 Apr 2025), https://www.business-standard.com/amp/markets/news/stock-market-trading-guide-"
           "apr-7-nikkei-hits-lower-limit-us-tariff-woes-125040700067_1.html | JPX trading hours as for 5 Aug 2024")
SRC_OCT = ("ABC News (AP): 100% tariff announced after the close of trading on Friday, https://abc13.com/post/stock-market-"
           "today-wall-street-drops-worst-day-april-trumps-threats-tariffs-shatter-calm/17980388/ | CNBC article on the "
           "post published 4:57 PM EDT (20:57 UTC), https://www.cnbc.com/2025/10/10/trump-trade-tariffs-china-software."
           "html | CoinDesk (21:51 UTC): post came late Friday afternoon, https://www.coindesk.com/markets/2025/10/10/"
           "bitcoin-crashes-below-usd110k-cryptos-in-freefall-on-further-trump-tariff-on-china | text archived at "
           "https://www.presidency.ucsb.edu/documents/truth-social-posts-october-10-2025 | earlier threat shortly before "
           "11 a.m. EDT, Bloomberg first story 15:07 UTC, https://www.bloomberg.com/news/articles/2025-10-10/trump-says-"
           "he-sees-no-reason-to-meet-xi-threatens-new-tariffs")

EVENTS = [
    dict(event="N1", date="2024-06-04", name="Indian general election results", market="nifty", stream="Nifty 2024-2025",
         T="2024-06-04 03:45", latest=None, local="09:15 IST (Nifty open)",
         rule="Counting began 08:00 IST (02:30 UTC) before the open, so T is the open", status="verified", sources=SRC_N1),
    dict(event="C1", date="2024-08-05", name="Global sell-off", market="crypto", stream="H2 2024",
         T="2024-08-05 00:00", latest=None, local="09:00 JST (Tokyo open)",
         rule="Monday sell-off that started in Asia: T is the Tokyo open", status="verified", sources=SRC_AUG),
    dict(event="N2", date="2024-08-05", name="Global sell-off", market="nifty", stream="Nifty 2024-2025",
         T="2024-08-05 03:45", latest=None, local="09:15 IST (Nifty open)",
         rule="News before the open: T is the open", status="verified", sources=SRC_AUG),
    dict(event="C2", date="2025-04-07", name="Tariff sell-off", market="crypto", stream="H1 2025",
         T="2025-04-07 00:00", latest=None, local="09:00 JST (Tokyo open)",
         rule="Monday sell-off that started in Asia: T is the Tokyo open", status="verified", sources=SRC_APR),
    dict(event="N3", date="2025-04-07", name="Tariff sell-off", market="nifty", stream="Nifty 2024-2025",
         T="2025-04-07 03:45", latest=None, local="09:15 IST (Nifty open)",
         rule="News before the open: T is the open", status="verified", sources=SRC_APR),
    dict(event="C3", date="2025-10-10", name="Crypto liquidation cascade", market="crypto", stream="H2 2025",
         T="2025-10-10 20:00", latest="2025-10-10 20:57", local="16:00 EDT (US stock close)",
         rule="100% tariff post after the US close; latest possible time 20:57 UTC (CNBC article)",
         status="marked: time of the post known only to within an hour (20:00 to 20:57 UTC)", sources=SRC_OCT),
]
for _e in EVENTS:
    _e["T"] = pd.Timestamp(_e["T"], tz="UTC")
    _e["latest"] = pd.Timestamp(_e["latest"], tz="UTC") if _e["latest"] else None

ASSETS = {"crypto": list(COINS), "nifty": ["Nifty 50"]}


def method_names(market):
    cfg = MARKETS[market]
    names = ["e-detector"]
    names += [f"BOCPD, w={w}" for w in cfg["ws"]]
    names += [f"BOCPD on truncated returns, w={w}" for w in cfg["ws"]]
    names += [f"BJV rolling, k={k}" for k in cfg["ks"]]
    return names


# ----------------------------------------------------------------------------
# Streams
# ----------------------------------------------------------------------------


def prepare_h2_2024():
    """H2 2024 crypto stream, built like E13's H1 2025 one: profile from 1 Jan 2024, monitoring from 1 Jul 2024."""
    prepared = {}
    for c, sym in COINS.items():
        df = load_binance(sym, "2024-01-01", "2025-01-01", download=False)
        s = intraday_scale(df.logret, lookback_days=20)
        u = (df.logret / s).to_numpy(copy=True)
        u[~np.isfinite(u)] = np.nan
        i0 = int(np.searchsorted(df.index, pd.Timestamp("2024-07-01", tz="UTC")))
        prepared[c] = (df.index, u, i0)
    return prepared


def crypto_streams():
    return {"H2 2024": prepare_h2_2024(), "H1 2025": E13.prepare_h1_2025(),
            "H2 2025": {c: prepare(c, "Jul-Dec 2025") for c in COINS}}


def nifty_stream():
    df, _, _ = load_clean("2023-10-01", "2026-01-01")
    u = deseasonalise(df).to_numpy()
    i0 = int(np.searchsorted(df.index, pd.Timestamp("2024-01-01", tz="Asia/Kolkata")))
    return {"Nifty 2024-2025": {"Nifty 50": (df.index, u, i0)}}


# ----------------------------------------------------------------------------
# Methods
# ----------------------------------------------------------------------------


def bocpd_job(args):
    um, prior, hazard, cal_len, truncated, w = args
    t0 = time.time()
    r_max = E13.R_MAX if w == E13.W_SHORT else 2 * w  # first run: cap 512; round 3: cap 2w (as in E13)
    _, S, L = E13.real_bocpd(um, prior, hazard, cal_len, truncated, r_max=r_max, short=w)
    print(f"BOCPD w={w} truncated={truncated}: {time.time() - t0:.0f} s for {um.shape}", flush=True)
    return w, truncated, S, L


def bocpd_all(market, series):
    """BOCPD arrays (S, L) for every series of a market, in one padded batch per window and version (E13's code)."""
    cfg = MARKETS[market]
    n_max = max(u.size - i0 for _, (_, u, i0) in series)
    um = np.full((len(series), n_max), np.nan)
    prior = np.full((len(series), n_max), np.nan)
    for j, (_, (_, u, i0)) in enumerate(series):
        um[j, : u.size - i0] = u[i0:]
        prior[j, : u.size - i0] = E13.trailing_robust_var(u, cfg["day"])[i0:]
    jobs = [(um, prior, 1 / cfg["A"], cfg["cal"], tr, w) for w in sorted(cfg["ws"], reverse=True) for tr in (False, True)]
    cache = os.environ.get("E18_CACHE")
    cfile = Path(cache) / f"bocpd_{market}.npz" if cache else None
    if cfile is not None and cfile.exists():  # all arrays of a market in one file
        z = np.load(cfile)
        return {(int(w), bool(tr)): (z[f"S_{w}_{tr}"], z[f"L_{w}_{tr}"]) for w, tr in z["keys"]}
    out = {}

    def jfile(w, tr):  # one file per job, written as soon as the job ends (survives an interrupted run)
        return Path(cache) / f"bocpd_{market}_w{w}_t{int(tr)}.npz" if cache else None

    todo = []
    for job in jobs:
        f = jfile(job[5], job[4])
        if f is not None and f.exists():
            z = np.load(f)
            out[(job[5], job[4])] = (z["S"], z["L"])
        else:
            todo.append(job)
    if todo:
        with Pool(min(2, len(todo))) as pool:
            for w, tr, S, L in pool.imap_unordered(bocpd_job, todo):
                out[(w, tr)] = (S, L)
                if cache:
                    np.savez(jfile(w, tr), S=S, L=L)
    return out


def run_edetector(um, cfg):
    """Stream the e-detector; ``monitoring[t]`` is its state when the return at t arrives (as in E12)."""
    mon = VolatilityMonitor(threshold=cfg["A"], c=3.0, cal_len=cfg["cal"], null_band=BAND)
    monitoring = np.zeros(um.size, bool)
    for t in range(um.size):
        monitoring[t] = mon.calibrated
        mon.update(um[t])
    starts = np.array([s["start"] for s in mon.segments], int)
    return [(e.index, e.direction) for e in mon.events], monitoring, starts


def match_masked(su, sd, n_up, n_dn, start, dead, keep, lo_hi, iters=50):
    """E13's bisection (``competitors.match_alarm_counts``) with alarms counted only where ``keep`` is True."""
    lo0, hi0 = lo_hi
    h_up, h_dn = hi0, hi0

    def counts(hu, hd):
        ev = C.scan_alarms(su, sd, hu, hd, start, dead)
        u = sum(1 for t, d in ev if d == "up" and keep[t])
        return u, sum(1 for t, d in ev if d == "down" and keep[t])

    for _ in range(4):
        lo, hi = lo0, hi0
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if counts(mid, h_dn)[0] > n_up else (lo, mid)
        h_up = hi
        lo, hi = lo0, hi0
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if counts(h_up, mid)[1] > n_dn else (lo, mid)
        h_dn = hi
    return h_up, h_dn


def competitor_stats(market, u, i0, bocpd_row):
    cfg = MARKETS[market]
    um = u[i0:]
    out = {}
    for w in cfg["ws"]:
        for tr in (False, True):
            S, L = bocpd_row[(w, tr)]
            nm = ("BOCPD on truncated returns" if tr else "BOCPD") + f", w={w}"
            out[nm] = (np.where(S > 0, L, -np.inf), np.where(S < 0, L, -np.inf), (-50.0, 1000.0))
    for k in cfg["ks"]:
        if k == E13.K_ONLINE:  # first run: statistic from the monitored returns only (E13 run_real_series)
            Z = C.bjv_online_statistic(um[None, :], k, E13.TRUNC)[0]
        else:  # round 3: history before the monitored period (E13 run_real_series_long)
            assert i0 >= 2 * k
            Z = C.bjv_online_statistic(u[None, i0 - 2 * k :], k, E13.TRUNC, norm="gaussian", method="sums")[0, 2 * k :]
        out[f"BJV rolling, k={k}"] = (Z, -Z, (0.0, 1000.0))
    return out


# ----------------------------------------------------------------------------
# Outcomes
# ----------------------------------------------------------------------------


def window_outcomes(alarms, tt, ev, kind, dead, monitoring=None, starts=None, rate=None):
    """Pre-registered outcomes 1 to 4 for one alarm list (positions into tt) around event ev."""
    T, latest = ev["T"], ev["latest"]
    t_ref = int(tt.searchsorted(T))
    pos = np.array([a[0] for a in alarms], int)
    dirs = np.array([a[1] for a in alarms])
    times = tt[pos] if pos.size else pd.DatetimeIndex([], tz=tt.tz)
    rec = dict(reference_bar=str(tt[t_ref]))
    before = pos[pos < t_ref]
    if kind == "e-detector":
        awake = bool(monitoring[t_ref])
        rec["state_at_T"] = "monitoring" if awake else "calibrating"
        later = starts[starts > t_ref]
        rec["resumes_at"] = "" if awake or not later.size else str(tt[later[0]])
    else:
        awake = not (before.size and t_ref <= before[-1] + dead)
        rec["state_at_T"] = "watching" if awake else "paused"
        rec["resumes_at"] = "" if awake else str(tt[before[-1] + dead + 1])
    rec["blinded_by_alarm_at"] = "" if awake or not before.size else f"{str(tt[before[-1]])} ({dirs[pos < t_ref][-1]})"
    pre = (times >= T - D7) & (times < T) if pos.size else np.zeros(0, bool)
    post = (times >= T) & (times < T + H24) if pos.size else np.zeros(0, bool)
    up = dirs == "up"
    rec["n_up_7d_before"] = int(np.sum(pre & up))
    rec["n_down_7d_before"] = int(np.sum(pre & ~up))
    if np.any(pre & up):
        last = times[pre & up][-1]
        rec["last_up_before_T"] = str(last)
        rec["last_up_lead_min"] = float((T - last) / pd.Timedelta(minutes=1))
    else:
        rec["last_up_before_T"], rec["last_up_lead_min"] = "", np.nan
    if rate is not None:
        steps = int(np.sum((tt >= T - D7) & (tt < T)))
        rec["expected_alarms_7d"] = rate * steps
    rec["n_up_24h_after"] = int(np.sum(post & up))
    rec["n_down_24h_after"] = int(np.sum(post & ~up))
    rec["first_alarm_after_T"] = f"{str(times[post][0])} ({dirs[post][0]})" if np.any(post) else ""
    if np.any(post & up):
        first = times[post & up][0]
        rec["first_up_after_T"] = str(first)
        rec["delay_min"] = float((first - T) / pd.Timedelta(minutes=1))
        rec["delay_from_latest_min"] = float((first - latest) / pd.Timedelta(minutes=1)) if latest is not None else np.nan
    else:
        rec["first_up_after_T"], rec["delay_min"], rec["delay_from_latest_min"] = "missed", np.nan, np.nan
    return rec


def alarm_rows(alarms, tt, ev, design, asset, method):
    T = ev["T"]
    rows = []
    for t, d in alarms:
        x = tt[t]
        if T - D7 <= x < T + H24:
            rows.append(dict(design=design, event=ev["event"], date=ev["date"], asset=asset, method=method,
                             time_utc=str(x.tz_convert("UTC")), time_local=str(x), direction=d,
                             minutes_from_T=float((x - T) / pd.Timedelta(minutes=1)),
                             window="7 days before" if x < T else "24 hours after"))
    return rows


def fresh_outcome(kind, ev, tt, um, cfg, su=None, h_up=None):
    """S1: an e-detector calibrated on the cal returns just before T and monitoring from T; for a competitor, its first
    up crossing at or after T ignoring any pause. Delay within 24 hours."""
    T = ev["T"]
    t_ref = int(tt.searchsorted(T))
    w = int(np.sum((tt >= T) & (tt < T + H24)))
    if kind == "e-detector":
        cal = cfg["cal"]
        mon = VolatilityMonitor(threshold=cfg["A"], c=3.0, cal_len=cal, null_band=BAND)
        mon.run(um[t_ref - cal : t_ref + w])
        ups = [e.index - cal for e in mon.events if e.direction == "up" and e.index >= cal]
        extra = dict(fresh_monitoring_from=str(tt[t_ref - cal + mon.segments[0]["start"]]) if mon.segments else "never")
        first = ups[0] if ups else None
    else:
        hit = np.flatnonzero(np.nan_to_num(su[t_ref : t_ref + w], nan=-np.inf) >= h_up)
        first = int(hit[0]) if hit.size else None
        extra = {}
    if first is None:
        return dict(first_up_after_T="missed", delay_min=np.nan, delay_from_latest_min=np.nan, **extra)
    x = tt[t_ref + first]
    lt = ev["latest"]
    return dict(first_up_after_T=str(x), delay_min=float((x - T) / pd.Timedelta(minutes=1)),
                delay_from_latest_min=float((x - lt) / pd.Timedelta(minutes=1)) if lt is not None else np.nan, **extra)


# ----------------------------------------------------------------------------
# One series: all methods, both matchings, all outcomes
# ----------------------------------------------------------------------------


def run_series(market, stream, asset, idx, u, i0, bocpd_row):
    cfg = MARKETS[market]
    um = u[i0:]
    tt = idx[i0:]
    n = um.size
    start = dead = cfg["cal"]
    evs = [e for e in EVENTS if e["stream"] == stream]
    keep = np.ones(n, bool)
    for e in evs:
        keep &= ~((tt >= e["T"] - EXCL_BEFORE) & (tt < e["T"] + EXCL_AFTER))
    eligible = keep & (np.arange(n) >= start)
    ed, monitoring, starts = run_edetector(um, cfg)
    nu_out = sum(1 for t, d in ed if d == "up" and keep[t])
    nd_out = sum(1 for t, d in ed if d == "down" and keep[t])
    nu_all = sum(1 for t, d in ed if d == "up")
    nd_all = len(ed) - nu_all
    stats_ = competitor_stats(market, u, i0, bocpd_row)
    alarms = {"primary": {"e-detector": ed}, "S3": {}}
    thr, rates, outc, alr = [], [], [], []
    for nm, (su, sd, lo_hi) in stats_.items():
        hu, hd = match_masked(su, sd, nu_out, nd_out, start, dead, keep, lo_hi)
        a = C.scan_alarms(su, sd, hu, hd, start, dead)
        alarms["primary"][nm] = a
        au = sum(1 for t, d in a if d == "up" and keep[t])
        thr.append(dict(market=market, stream=stream, asset=asset, method=nm, matching="outside event windows (primary)",
                        h_up=hu, h_down=hd, target_up=nu_out, target_down=nd_out, achieved_up=au,
                        achieved_down=sum(1 for t, d in a if d == "down" and keep[t])))
        hu3, hd3 = C.match_alarm_counts(su, sd, nu_all, nd_all, start, dead, lo_hi=lo_hi, iters=50)
        a3 = C.scan_alarms(su, sd, hu3, hd3, start, dead)
        alarms["S3"][nm] = a3
        thr.append(dict(market=market, stream=stream, asset=asset, method=nm, matching="whole sample (E13 rule, S3)",
                        h_up=hu3, h_down=hd3, target_up=nu_all, target_down=nd_all,
                        achieved_up=sum(1 for t, d in a3 if d == "up"), achieved_down=sum(1 for t, d in a3 if d == "down")))
    for design in ("primary", "S3"):
        for nm, a in alarms[design].items():
            n_out = sum(1 for t, _ in a if keep[t])
            rate = n_out / max(int(eligible.sum()), 1)
            if design == "primary":
                rates.append(dict(market=market, stream=stream, asset=asset, method=nm, monitored_steps=n,
                                  steps_outside_windows=int(eligible.sum()), alarms=len(a), alarms_outside_windows=n_out,
                                  alarms_per_day_outside_windows=rate * cfg["day"]))
            kind = "e-detector" if nm == "e-detector" else "competitor"
            for e in evs:
                rec = window_outcomes(a, tt, e, kind, dead, monitoring, starts, rate)
                outc.append(dict(design=design, event=e["event"], date=e["date"], market=market, stream=stream,
                                 asset=asset, method=nm, **rec))
                alr += alarm_rows(a, tt, e, design, asset, nm)
    for e in evs:  # S1: fresh start at T
        rec = fresh_outcome("e-detector", e, tt, um, cfg)
        outc.append(dict(design="S1", event=e["event"], date=e["date"], market=market, stream=stream, asset=asset,
                         method="e-detector", **rec))
        for nm, (su, _, _) in stats_.items():
            h_up = next(r["h_up"] for r in thr if r["method"] == nm and r["matching"].startswith("outside"))
            rec = fresh_outcome("competitor", e, tt, um, cfg, su=su, h_up=h_up)
            outc.append(dict(design="S1", event=e["event"], date=e["date"], market=market, stream=stream, asset=asset,
                             method=nm, **rec))
    check = dict(stream=stream, asset=asset, edetector_alarms=[(tt[t], d) for t, d in ed])
    return outc, alr, thr, rates, check


# ----------------------------------------------------------------------------
# S2: the e-detector streamed from another start date
# ----------------------------------------------------------------------------


def s2_long_streams():
    rows, alr = [], []
    jobs = []
    for c, sym in COINS.items():
        df = load_binance(sym, "2024-01-01", "2026-01-01", download=False)  # E12's stream design
        s = intraday_scale(df.logret, lookback_days=20)
        u = (df.logret / s).to_numpy(copy=True)
        u[~np.isfinite(u)] = np.nan
        i0 = int(np.searchsorted(df.index, pd.Timestamp("2024-01-21", tz="UTC")))
        jobs.append(("crypto", c, df.index, u, i0, "crypto, 21 Jan 2024 to 31 Dec 2025 (E12 stream)"))
    df, _, _ = load_clean()  # E7's stream: every session from April 2017
    jobs.append(("nifty", "Nifty 50", df.index, deseasonalise(df).to_numpy(), 0, "Nifty, April 2017 onwards (E7 stream)"))
    for market, asset, idx, u, i0, label in jobs:
        cfg = MARKETS[market]
        um, tt = u[i0:], idx[i0:]
        ed, monitoring, starts = run_edetector(um, cfg)
        evs = [e for e in EVENTS if e["market"] == market]
        keep = np.ones(um.size, bool)
        for e in evs:
            keep &= ~((tt >= e["T"] - EXCL_BEFORE) & (tt < e["T"] + EXCL_AFTER))
        eligible = keep & (np.arange(um.size) >= cfg["cal"])
        rate = sum(1 for t, _ in ed if keep[t]) / max(int(eligible.sum()), 1)
        for e in evs:
            rec = window_outcomes(ed, tt, e, "e-detector", cfg["cal"], monitoring, starts, rate)
            rows.append(dict(design="S2", event=e["event"], date=e["date"], market=market, stream=label, asset=asset,
                             method="e-detector", **rec))
            alr += alarm_rows(ed, tt, e, "S2", asset, "e-detector")
        print(f"S2 done {asset}", flush=True)
    return rows, alr


# ----------------------------------------------------------------------------
# Reproduction checks against E13's saved files (inclusive matching = E13's rule)
# ----------------------------------------------------------------------------

E13_FIRST = {"BOCPD, w=250": "BOCPD (log-odds)", "BOCPD on truncated returns, w=250": "BOCPD on truncated returns (log-odds)",
             "BJV rolling, k=120": "BJV rolling"}


def e13_reference(thr):
    """Add E13's band-0.25 thresholds (where E13 ran the same stream) next to our whole-sample (S3) thresholds."""
    first = {"H1 2025": "e13d_real_crypto_h1_2025_matched_thresholds.csv", "H2 2025": "e13d_real_crypto_matched_thresholds.csv",
             "Nifty 2024-2025": "e13d_real_nifty_matched_thresholds.csv"}
    longm = {"H1 2025": "five Binance coins, Jan to Jun 2025", "H2 2025": "five Binance coins, Jul to Dec 2025",
             "Nifty 2024-2025": "Nifty 50"}
    lw = pd.concat([E13._read(OUT / "e13d_real_long_windows_crypto_matched_thresholds.csv"),
                    E13._read(OUT / "e13d_real_long_windows_nifty_matched_thresholds.csv")], ignore_index=True)
    lw = lw[lw.band == 0.25]
    e13_up, e13_dn = [], []
    for _, r in thr.iterrows():
        hu = hd = np.nan
        if r.stream in first and r.matching.startswith("whole"):
            if r.method in E13_FIRST:
                f = E13._read(OUT / first[r.stream])
                m = f[(f.series == r.asset) & (f.band == 0.25) & (f.method == E13_FIRST[r.method])]
            else:
                m = lw[(lw.market == longm[r.stream]) & (lw.series == r.asset) & (lw.method == r.method)]
            if len(m) == 1:
                hu, hd = float(m.h_up.iloc[0]), float(m.h_down.iloc[0])
        e13_up.append(hu)
        e13_dn.append(hd)
    thr["e13_h_up"], thr["e13_h_down"] = e13_up, e13_dn
    return thr


def e13_edetector_check(checks):
    """Are the e-detector's alarms on E13's streams exactly E13's saved band-0.25 alarms?"""
    files = {"H1 2025": "e13d_real_crypto_h1_2025_events.csv", "H2 2025": "e13d_real_crypto_events.csv",
             "Nifty 2024-2025": "e13d_real_nifty_events.csv"}
    out = []
    for ch in checks:
        if ch["stream"] not in files:
            continue
        ev = E13._read(OUT / files[ch["stream"]])
        ev = ev[(ev.series == ch["asset"]) & (ev.method == "e-detector, band 0.25")]
        mine = ch["edetector_alarms"]
        same = len(ev) == len(mine) and all(pd.Timestamp(a) == b and d0 == d1
                                            for a, d0, (b, d1) in zip(ev.time, ev.direction, mine))
        out.append(dict(stream=ch["stream"], asset=ch["asset"], e13_alarms=len(ev), e18_alarms=len(mine),
                        identical_times_and_directions=same))
    return out


# ----------------------------------------------------------------------------
# Main run
# ----------------------------------------------------------------------------


def run():
    t0 = time.time()
    ev_tab = pd.DataFrame([dict(event=e["event"], date=e["date"], name=e["name"], market=e["market"],
                                assets=", ".join(ASSETS[e["market"]]), stream=e["stream"],
                                reference_utc=str(e["T"]), reference_local=e["local"],
                                latest_reference_utc=str(e["latest"]) if e["latest"] is not None else "",
                                rule=e["rule"], status=e["status"], sources=e["sources"]) for e in EVENTS])
    ev_tab.to_csv(OUT / "e18_events.csv", index=False)
    import pickle

    cache = os.environ.get("E18_CACHE")
    outc, alr, thr, rates, checks = [], [], [], [], []
    for market, streams in (("nifty", nifty_stream()), ("crypto", crypto_streams())):
        series = [((st, a), streams[st][a]) for st in streams for a in streams[st]]
        bocpd = bocpd_all(market, series)
        print(f"{market}: BOCPD done after {time.time() - t0:.0f} s", flush=True)
        for j, ((st, a), (idx, u, i0)) in enumerate(series):
            pfile = Path(cache) / f"series_{st}_{a}.pkl".replace(" ", "_") if cache else None
            if pfile is not None and pfile.exists():  # finished in an earlier, interrupted attempt
                with open(pfile, "rb") as fh:
                    o, al, th, ra, ch = pickle.load(fh)
            else:
                n = u.size - i0
                row = {key: (S[j, :n], L[j, :n]) for key, (S, L) in bocpd.items()}
                o, al, th, ra, ch = run_series(market, st, a, idx, u, i0, row)
                if pfile is not None:
                    with open(pfile, "wb") as fh:
                        pickle.dump((o, al, th, ra, ch), fh)
            outc += o
            alr += al
            thr += th
            rates += ra
            checks.append(ch)
            print(f"done {st} {a} after {time.time() - t0:.0f} s", flush=True)
        del bocpd
    o, al = s2_long_streams()
    outc += o
    alr += al
    pd.DataFrame(outc).to_csv(OUT / "e18_outcomes.csv", index=False)
    pd.DataFrame(alr).to_csv(OUT / "e18_alarms.csv", index=False)
    thr = e13_reference(pd.DataFrame(thr))
    thr.to_csv(OUT / "e18_thresholds.csv", index=False)
    pd.DataFrame(rates).to_csv(OUT / "e18_rates.csv", index=False)
    chk = pd.DataFrame(e13_edetector_check(checks))
    print(chk.to_string(index=False))
    s3 = thr[thr.matching.str.startswith("whole") & thr.e13_h_up.notna()]
    rel = np.abs(np.r_[s3.h_up - s3.e13_h_up, s3.h_down - s3.e13_h_down]) / np.maximum(
        np.abs(np.r_[s3.e13_h_up, s3.e13_h_down]), 1e-12)
    print(f"S3 thresholds against E13's: {len(s3)} rows, largest relative difference {rel.max():.2e}")
    print(f"total {time.time() - t0:.0f} s")


# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------

FAMILY = {"e-detector": ("e-detector", "#2a78d6"), "BOCPD,": ("BOCPD", "#eb6834"),
          "BOCPD on": ("BOCPD on truncated returns", "#1baf7a"), "BJV": ("two-window variance-ratio scan", "#6250d6")}


def _family(m):
    for key, (lab, col) in FAMILY.items():
        if m.startswith(key):
            return lab, col
    raise KeyError(m)


def _short(m):
    if m == "e-detector":
        return "e-detector"
    w = m.split("=")[-1]
    if m.startswith("BOCPD on"):
        return f"BOCPD trunc. {w}"
    if m.startswith("BOCPD"):
        return f"BOCPD {w}"
    return f"scan {w}"


def raw_returns(market, asset, T):
    """Raw one-minute log returns over [T - 7 d, T + 24 h) from data/raw (not saved anywhere)."""
    if market == "crypto":
        a = (T - D7).strftime("%Y-%m-%d")
        b = (T + H24 + pd.Timedelta(days=1)).strftime("%Y-%m-%d")
        df = load_binance(COINS[asset], a, b, download=False)
        r = df.logret
    else:
        df, _, _ = load_clean((T - D7 - pd.Timedelta(days=1)).strftime("%Y-%m-%d"),
                              (T + H24 + pd.Timedelta(days=1)).strftime("%Y-%m-%d"))
        r = df.logret
    return r[(r.index >= T - D7) & (r.index < T + H24)]


def figure():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    al = pd.read_csv(OUT / "e18_alarms.csv")
    al = al[al.design == "primary"]
    dates = {"crypto": ["2024-08-05", "2025-04-07", "2025-10-10"], "nifty": ["2024-06-04", "2024-08-05", "2025-04-07"]}
    rows = list(COINS) + ["Nifty 50"]
    ink, muted, grid = "#0b0b0b", "#52514e", "#e6e5e1"
    fig = plt.figure(figsize=(16, 23), facecolor="#fcfcfb")
    outer = fig.add_gridspec(len(rows), 3, hspace=0.42, wspace=0.13, left=0.085, right=0.99, top=0.94, bottom=0.04)
    for i, asset in enumerate(rows):
        market = "nifty" if asset == "Nifty 50" else "crypto"
        names = method_names(market)
        for j, d in enumerate(dates[market]):
            ev = next(e for e in EVENTS if e["date"] == d and e["market"] == market)
            inner = outer[i, j].subgridspec(2, 1, height_ratios=[1.3, 1.0], hspace=0.05)
            ax, lx = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            for a_ in (ax, lx):
                a_.set_facecolor("#fcfcfb")
                for s in ("top", "right"):
                    a_.spines[s].set_visible(False)
                a_.tick_params(labelsize=7, colors=muted)
            T = ev["T"]
            r = raw_returns(market, asset, T)
            if market == "crypto":
                x = (r.index - T) / pd.Timedelta(hours=1)
                xa = lambda s: (pd.to_datetime(s, utc=True) - T) / pd.Timedelta(hours=1)  # noqa: E731
                xlab = "hours from T"
            else:
                pos = np.arange(len(r))
                ref = int(r.index.searchsorted(T))
                x = (pos - ref) / 375.0
                where = pd.Series(x, index=r.index.tz_convert("UTC"))
                xa = lambda s: where.reindex(pd.to_datetime(s, utc=True)).to_numpy()  # noqa: E731
                xlab = "sessions from T (trading minutes / 375)"
                for b in np.flatnonzero(np.r_[True, np.diff(r.index.normalize()) > pd.Timedelta(0)]):
                    for a_ in (ax, lx):
                        a_.axvline(x[b], color=grid, lw=0.6, zorder=0)
            ax.plot(np.asarray(x), 100 * r.to_numpy(), color="#8a8984", lw=0.35)
            ax.set_ylabel("return, %", fontsize=7, color=muted)
            ax.set_xticklabels([])
            ax.set_title(f"{asset}, {pd.Timestamp(d):%d %b %Y} ({ev['event']})", fontsize=9, color=ink, loc="left")
            for a_ in (ax, lx):
                a_.axvline(0, color=ink, ls="--", lw=0.8)
                if ev["latest"] is not None:
                    a_.axvspan(0, (ev["latest"] - T) / pd.Timedelta(hours=1), color="#d9d8d3", lw=0, zorder=0)
            sub = al[(al.event == ev["event"]) & (al.asset == asset)]
            for k, m in enumerate(names):
                _, col = _family(m)
                y = len(names) - 1 - k
                s = sub[sub.method == m]
                for direc, mk in (("up", "^"), ("down", "v")):
                    xs = xa(s[s.direction == direc].time_utc)
                    lx.scatter(xs, np.full(len(xs), y), marker=mk, s=22 if direc == "up" else 16,
                               color=col if direc == "up" else "none", edgecolors=col, linewidths=0.9, zorder=3)
                lx.axhline(y, color=grid, lw=0.5, zorder=0)
            lx.set_ylim(-0.7, len(names) - 0.3)
            lx.set_yticks(range(len(names)))
            first_col = j == 0
            lx.set_yticklabels([_short(m) for m in names[::-1]] if first_col else [], fontsize=6.5, color=muted)
            lx.tick_params(axis="y", length=0)
            lx.set_xlabel(xlab, fontsize=7, color=muted)
            lo = float(np.nanmin(np.asarray(x)))
            hi = float(np.nanmax(np.asarray(x)))
            for a_ in (ax, lx):
                a_.set_xlim(lo, hi)
    handles = [Line2D([], [], color=c, marker="s", ls="", ms=7, label=lab) for lab, c in FAMILY.values()]
    handles += [Line2D([], [], color=muted, marker="^", ls="", ms=7, label="up alarm"),
                Line2D([], [], color=muted, marker="v", ls="", ms=7, mfc="none", label="down alarm"),
                Line2D([], [], color=ink, ls="--", lw=0.8, label="reference time T")]
    fig.legend(handles=handles, loc="upper center", ncol=7, fontsize=9, frameon=False, bbox_to_anchor=(0.53, 0.975))
    fig.suptitle("E18 [REAL], case studies: one-minute log returns and alarms, 7 days before to 24 hours after the "
                 "reference time T (grey band on 10 Oct 2025: the post came between 20:00 and 20:57 UTC)",
                 fontsize=11, color=ink, y=0.99)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "e18_events.png", dpi=110, facecolor=fig.get_facecolor())
    print("wrote", FIG / "e18_events.png")


def context():
    """Descriptive context, NOT pre-registered (added after the run): how large the raw minute returns were around the
    e-detector's last up alarm before T on the crypto events. Size = mean absolute one-minute log return in a span,
    divided by its mean over the first five of the seven days before T. Price change = exp(sum of log returns) - 1."""
    o = pd.read_csv(OUT / "e18_outcomes.csv")
    h = pd.Timedelta(hours=1)
    rows = []
    for ev in EVENTS:
        if ev["market"] != "crypto":
            continue
        T = pd.Timestamp(ev["T"])
        for a in COINS:
            r = raw_returns("crypto", a, T)
            x = r.abs()
            sel = o[(o.event == ev["event"]) & (o.asset == a) & (o.method == "e-detector") & (o.design == "primary")]
            tA = pd.Timestamp(sel.iloc[0].last_up_before_T)
            base = x[(x.index >= T - D7) & (x.index < T - D7 + pd.Timedelta(days=5))].mean()

            def size(lo, hi):
                return float(x[(x.index >= lo) & (x.index < hi)].mean() / base)

            def move(lo, hi):
                return float(100 * (np.exp(r[(r.index >= lo) & (r.index < hi)].sum()) - 1))

            rows.append(dict(event=ev["event"], asset=a, edetector_last_up_before_T=tA,
                             size_3h_before_alarm=size(tA - 3 * h, tA), size_alarm_to_T=size(tA, T),
                             size_T_to_T6h=size(T, T + 6 * h), size_T6h_to_T24h=size(T + 6 * h, T + 24 * h),
                             price_change_pct_3h_before_alarm=move(tA - 3 * h, tA),
                             price_change_pct_alarm_to_T=move(tA, T), price_change_pct_T_to_T24h=move(T, T + 24 * h)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "e18_context.csv", index=False, float_format="%.4g")
    print(out.round({c: 2 for c in out.columns if c.startswith(("size", "price"))}).to_string(index=False))
    print("wrote", OUT / "e18_context.csv")


# ----------------------------------------------------------------------------
# Exploratory, NOT pre-registered: E14's monitors with no blind period on the same events
# ----------------------------------------------------------------------------
# Decided after seeing that the e-detector was recalibrating at every crypto event. Settings exactly as committed in
# E14 (experiments/e14_running_tau.py, make_monitor and real_job): awake-A = route A with a learning band and no
# blind period (AwakeRunningTauMonitor, A = day / (0.05 - 0.01), delta = 0.01, c = 3); awake-B = the rank (conformal)
# detector restarted at once after each alarm (ConformalMonitor, A = 28,800 or 7,500, 120 returns of warm-up at the
# stream start only, family "lr", seed 14000). Both with E14's two null bands, 0.25 (E12's frozen setting, the one the
# e-detector uses here) and 0. Same streams, events and reference times as the primary design.

E14_SEED = 14_000
DELTA_AWAKE = 0.01
AWAKE_MONITORS = [(m, b) for m in ("awake-A", "awake-B") for b in (0.25, 0.0)]


def make_awake(method, band, cfg):
    from avmonitor import running as RT

    if method == "awake-A":
        return RT.AwakeRunningTauMonitor(threshold=cfg["day"] / (0.05 - DELTA_AWAKE), c=3.0, delta=DELTA_AWAKE,
                                         null_band=band)
    return RT.ConformalMonitor(threshold=cfg["A"], cal_len=RT.AWAKE_PRE, family="lr", band=band, seed=E14_SEED,
                               restart_len=0)


def awake_job(job):
    market, stream, asset, tt, um = job
    cfg = MARKETS[market]
    evs = [e for e in EVENTS if e["stream"] == stream and tt[0] <= e["T"] and e["T"] + H24 <= tt[-1]]
    keep = np.ones(um.size, bool)
    for e in evs:
        keep &= ~((tt >= e["T"] - EXCL_BEFORE) & (tt < e["T"] + EXCL_AFTER))
    rows, lists = [], {}
    for method, band in AWAKE_MONITORS:
        t0 = time.time()
        mon = make_awake(method, band, cfg)
        watching = np.zeros(um.size, bool)
        for t in range(um.size):
            watching[t] = mon.calibrated
            mon.update(um[t])
        alarms = [(e.index, e.direction) for e in mon.events]
        lists[(method, band)] = alarms
        starts = np.array([s["start"] for s in mon.segments], int)
        n_out = sum(1 for t, _ in alarms if keep[t])
        rate = n_out / max(int(keep.sum()), 1)
        for e in evs:
            rec = window_outcomes(alarms, tt, e, "e-detector", None, watching, starts, rate)
            T = e["T"]
            pos = np.array([a[0] for a in alarms], int)
            before = pos[pos < int(tt.searchsorted(T))]
            steps7 = int(np.sum((tt >= T - D7) & (tt < T)))
            n7 = rec["n_up_7d_before"] + rec["n_down_7d_before"]
            rows.append(dict(
                event=e["event"], date=e["date"], market=market, stream=stream, asset=asset, monitor=method,
                null_band=band, A=cfg["day"] / (0.05 - DELTA_AWAKE) if method == "awake-A" else cfg["A"],
                watching_at_T=rec["state_at_T"] == "monitoring",
                last_alarm_before_T=str(tt[before[-1]]) if before.size else "",
                last_alarm_direction=dict(alarms)[int(before[-1])] if before.size else "",
                last_alarm_lead_min=float((T - tt[before[-1]]) / pd.Timedelta(minutes=1)) if before.size else np.nan,
                n_up_7d_before=rec["n_up_7d_before"], n_down_7d_before=rec["n_down_7d_before"],
                alarms_per_day_7d_before=n7 / (steps7 / cfg["day"]),
                first_up_after_T=rec["first_up_after_T"], delay_min=rec["delay_min"],
                delay_from_latest_min=rec["delay_from_latest_min"],
                n_up_24h_after=rec["n_up_24h_after"], n_down_24h_after=rec["n_down_24h_after"],
                first_alarm_after_T=rec["first_alarm_after_T"],
                stream_alarms_per_day=len(alarms) / (um.size / cfg["day"]),
                stream_alarms_per_day_outside_windows=rate * cfg["day"], watching_share=float(watching.mean())))
        print(f"awake {stream} {asset} {method} band {band}: {len(alarms)} alarms, {time.time() - t0:.0f} s", flush=True)
    return rows, (stream, asset, lists)


def e14_awake_check(lists):
    """Compare with E14's saved real-data alarm lists (H2 2025 coins and Nifty 2024 to 2025, same steps)."""
    f = OUT / "e14d_real_alarms.csv"
    if not f.exists():
        return {}
    ref = pd.read_csv(f)
    out = {}
    for stream, asset, by_mon in lists:
        src = COINS.get(asset, asset)
        if stream not in ("H2 2025", "Nifty 2024-2025"):
            continue
        for (method, band), alarms in by_mon.items():
            r = ref[(ref.source == src) & (ref.method == method) & (ref.band == band)]
            same = [(int(i), d) for i, d in zip(r["index"], r["direction"])] == [(int(i), d) for i, d in alarms]
            out[(stream, asset, method, band)] = same
            print(f"E14 check {stream} {asset} {method} band {band}: {'identical' if same else 'DIFFERENT'} "
                  f"({len(alarms)} alarms here, {len(r)} in E14)", flush=True)
    return out


def awake():
    """Exploratory, NOT pre-registered: results/e18_awake.csv (see the block comment above)."""
    t0 = time.time()
    jobs = []
    for market, streams in (("nifty", nifty_stream()), ("crypto", crypto_streams())):
        for st in streams:
            for a, (idx, u, i0) in streams[st].items():
                jobs.append((market, st, a, idx[i0:], u[i0:]))
    with Pool(2) as pool:
        res = pool.map(awake_job, jobs, chunksize=1)
    rows = [r for rr, _ in res for r in rr]
    check = e14_awake_check([lst for _, lst in res])
    out = pd.DataFrame(rows)
    out["e14_alarm_list_check"] = [
        ("identical" if check[k] else "different") if k in check else "no E14 run"
        for k in zip(out.stream, out.asset, out.monitor, out.null_band)]
    o = pd.read_csv(OUT / "e18_outcomes.csv")
    ed = o[(o.design == "primary") & (o.method == "e-detector")].set_index(["event", "asset"])
    key = list(zip(out.event, out.asset))
    out["edetector_n_7d_before"] = [int(ed.loc[k, "n_up_7d_before"] + ed.loc[k, "n_down_7d_before"]) for k in key]
    out["edetector_first_up_after_T"] = [ed.loc[k, "first_up_after_T"] for k in key]
    out.to_csv(OUT / "e18_awake.csv", index=False)
    cols = ["event", "asset", "monitor", "null_band", "watching_at_T", "last_alarm_before_T", "last_alarm_direction",
            "delay_min", "n_up_7d_before", "n_down_7d_before", "alarms_per_day_7d_before", "stream_alarms_per_day"]
    print(out[cols].round(2).to_string(index=False))
    print(f"wrote {OUT / 'e18_awake.csv'} after {time.time() - t0:.0f} s")


if __name__ == "__main__":
    parts = sys.argv[1:] or ["run", "figure"]
    for p in parts:
        {"run": run, "figure": figure, "context": context, "awake": awake}[p]()
