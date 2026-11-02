"""E6b. The streaming monitor on five Binance coins, in and out of sample, across venues and against forecasts.

REAL MARKET DATA: Binance spot 1-minute klines (BTCUSDT, ETHUSDT, XRPUSDT,
SOLUSDT, DOGEUSDT) via ``avmonitor.data.load_binance``, and Bitstamp BTC/USD
1-minute bars (data/raw/crypto/btcusd_bitstamp_1min_latest.csv) for the
cross-venue check.

Settings are FROZEN from E5 (chosen on Bitstamp BTC, Jul to Dec 2025) and not
tuned here: A = 28,800, each side alarming at 2A (at most 5% chance of any false alarm per day
when the null holds), c = 3, cal_len = 1440 minutes, null_band 0 / 0.25 / 0.5,
returns divided by the predictable intraday profile ``intraday_scale`` (20-day
lookback). The Gaussian likelihood ratio monitor (band 0) is the benchmark.

Research questions
  Q6d. Alarms per day and the E5c validation ("followed by higher / lower
       volatility over the next 1440 minutes than the calibration baseline")
       with Wilson 95% intervals, against random times, for all five coins in
       Jul to Dec 2025 and OUT OF SAMPLE in 21 Jan 2024 to 30 Jun 2025 (same
       frozen settings). Cross-venue consistency: how many Binance BTCUSDT
       alarms have a Bitstamp BTC/USD alarm within 60 minutes, against
       randomly shifted alarm times?
  Q6e. Does an alarm carry information beyond a volatility forecast? At each
       alarm, forecast the robust volatility of the next 1440 minutes with
       (i) an EWMA of |deseasonalised return| (half-life 1440 minutes, scaled
       by a factor estimated on the previous 30 days only) and (ii) a HAR
       regression on log robust volatility (1-day, 5-day and 22-day averages),
       refitted each day on past data only. Is realised volatility above the
       forecast after "up" alarms and below it after "down" alarms more often
       than at random times? Rate-matched competitor: a non-sequential rule that
       alarms when short-window (60 or 240 minutes) robust volatility over the
       preceding 1440-minute robust volatility crosses thresholds chosen to give
       the same numbers of up and down alarms as the e-detector, with the same
       1440-minute dead time after each alarm. This validates usefulness only;
       it says nothing about the false alarm guarantee.

Gaps: missing minutes are NaN rows in the loaders, so the return before and
after a gap is NaN (no multi-minute returns). The monitor skips NaN returns
(statistic unchanged); NaN values in a calibration window are dropped. Robust
volatilities are medians over the finite returns of a window and need at least
1000 of 1440 (the Binance mirror has no missing minutes in either period, so
this only matters for robustness of the code).

Outputs: results/e06d_monitor_summary.csv, e06d_events.csv,
e06d_cross_venue.csv, e06d_cross_coin.csv, e06e_forecast_hits.csv,
e06e_competitor.csv, e06e_pooled.csv, e06d_pooled.csv, results/figures/e06_alarm_timelines.png,
results/figures/e06_hit_rates.png. Seeds: 6 (random times), 66 (shifts).
"""

from __future__ import annotations

import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from avmonitor.data import intraday_scale, load_binance, load_bitstamp_csv
from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
BITSTAMP = ROOT / "data/raw/crypto/btcusd_bitstamp_1min_latest.csv"
COINS = {"BTC": "BTCUSDT", "ETH": "ETHUSDT", "XRP": "XRPUSDT", "SOL": "SOLUSDT", "DOGE": "DOGEUSDT"}
# label -> (monitor start, end, history start used for the intraday profile and the forecasts)
PERIODS = {
    # History from Jan 2024 (all of it before the monitored period) so that the HAR regression has a long
    # training sample; with only 3 months of history its median log(realised / forecast) at random times was
    # +0.18 for BTC (iteration 5 in docs/iterations/E6.md).
    "Jul-Dec 2025": ("2025-07-01", "2026-01-01", "2024-01-01"),
    "Jan 2024-Jun 2025 (out of sample)": ("2024-01-21", "2025-07-01", "2024-01-01"),
}
A, C, CAL = 28800.0, 3.0, 1440
CONFIGS = [("robust e-detector", 0.0), ("robust e-detector", 0.25), ("robust e-detector", 0.5), ("Gaussian LR", 0.0)]
DAY = 1440
MIN_FINITE = 1000
N_RANDOM = 3000
EWMA_HALFLIFE = 1440
SHORT_WINDOWS = (60, 240)
MATCH_MINUTES = 60


# ----------------------------------------------------------------------------
# Data and monitor runs
# ----------------------------------------------------------------------------


def prepare(source: str, period: str):
    start, end, hist = PERIODS[period]
    if source == "Bitstamp BTC":
        df = load_bitstamp_csv(BITSTAMP, hist, end)
    else:
        df = load_binance(COINS[source], hist, end, download=False)
    s = intraday_scale(df.logret, lookback_days=20)  # NaN returns are skipped by the rolling medians
    u = (df.logret / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    i0 = int(np.searchsorted(df.index, pd.Timestamp(start, tz="UTC")))
    return df.index, u, i0


def run_job(job):
    source, period = job
    idx, u, i0 = prepare(source, period)
    out = {}
    for name, band in CONFIGS:
        cls = VolatilityMonitor if name.startswith("robust") else GaussianLRMonitor
        mon = cls(threshold=A, c=C, cal_len=CAL, null_band=band)
        ev = mon.run(u[i0:])  # NaN passed through: skipped by the monitor
        out[(name, band)] = [(e.index + i0, e.direction, e.monitored_since + i0, e.baseline_scale) for e in ev]
    return source, period, out


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def wilson(k: int, n: int, z: float = 1.959964):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, centre - half, centre + half


def rolling_robust_sd(u: np.ndarray, w: int, min_periods: int) -> np.ndarray:
    """Entry t: scaled median of |u| over u[t-w+1 .. t] (finite values only)."""
    return pd.Series(np.abs(u)).rolling(w, min_periods=min_periods).median().to_numpy() * MAD_NORMAL_SCALE


class Series:
    """Precomputed robust volatilities and forecasts for one deseasonalised return series."""

    def __init__(self, idx, u, i0):
        self.idx, self.u, self.i0, self.n = idx, u, i0, u.size
        self.v1 = rolling_robust_sd(u, DAY, MIN_FINITE)  # (t-1439 .. t]
        self.y = np.full(self.n, np.nan)  # realised robust vol over (t, t+1440]
        self.y[: self.n - DAY] = self.v1[DAY:]
        self.prev = np.full(self.n, np.nan)  # (t-1440 .. t-1]
        self.prev[1:] = self.v1[:-1]
        self.ewma = pd.Series(np.abs(u)).ewm(halflife=EWMA_HALFLIFE, ignore_na=True).mean().to_numpy()
        self.rho = {}
        for S in SHORT_WINDOWS:
            vs = rolling_robust_sd(u, S, int(0.8 * S))
            rho = np.full(self.n, np.nan)
            rho[S:] = vs[S:] / self.v1[:-S]  # short window ending at t over the 1440 minutes before it
            self.rho[S] = rho
        self._fit_har()

    # HAR on log robust volatility: features at t use data up to t; refit once per day on samples whose
    # target window ended before that day started.
    def har_x(self, t):
        t = np.atleast_1d(t)
        lags = [self.v1[np.clip(t - DAY * j, 0, None)] for j in range(22)]
        lags = np.array(lags)
        bad = (t - DAY * 21) < DAY
        x1 = lags[0]
        x5 = lags[:5].mean(axis=0)
        x22 = lags.mean(axis=0)
        X = np.column_stack([np.ones(t.size), np.log(x1), np.log(x5), np.log(x22)])
        X[bad] = np.nan
        return X

    def _fit_har(self):
        grid = np.arange(22 * DAY, self.n - DAY, 60)
        X = self.har_x(grid)
        yy = np.log(self.y[grid])
        ok = np.all(np.isfinite(X), axis=1) & np.isfinite(yy)
        grid, X, yy = grid[ok], X[ok], yy[ok]
        days = np.arange(0, self.n, DAY)  # the series starts at midnight UTC, so these are day boundaries
        self.har_days = days
        self.har_coef = np.full((days.size, 4), np.nan)
        for k, d0 in enumerate(days):
            m = grid + DAY <= d0
            if m.sum() >= 240:  # at least 10 days of hourly samples
                self.har_coef[k] = np.linalg.lstsq(X[m], yy[m], rcond=None)[0]

    def forecast_har(self, t):
        t = np.atleast_1d(t)
        k = np.searchsorted(self.har_days, t, side="right") - 1
        coef = np.where(k[:, None] >= 0, self.har_coef[np.clip(k, 0, None)], np.nan)
        return np.exp(np.sum(self.har_x(t) * coef, axis=1))

    def forecast_ewma(self, t):
        t = np.atleast_1d(t)
        js = np.arange(1, 31)
        s = t[:, None] - DAY * js[None, :]
        valid = s >= 0
        s = np.clip(s, 0, None)
        ratio = np.where(valid, self.y[s] / self.ewma[s], np.nan)
        with np.errstate(all="ignore"):
            k = np.nanmedian(np.where(np.isfinite(ratio), ratio, np.nan), axis=1)
        enough = np.sum(np.isfinite(ratio), axis=1) >= 10
        return np.where(enough, k * self.ewma[t], np.nan)


def window_sd(u, a, b):
    seg = u[max(a, 0) : b]
    seg = seg[np.isfinite(seg)]
    return np.median(np.abs(seg)) * MAD_NORMAL_SCALE if seg.size >= MIN_FINITE else np.nan


# ----------------------------------------------------------------------------
# Q6d: summary with Wilson intervals, random-time controls
# ----------------------------------------------------------------------------


def summarise(source, period, S: Series, events, rng):
    rows, ev_rows = [], []
    days = (S.n - S.i0) / DAY
    lo, hi = S.i0 + DAY, S.n - DAY - 1
    rand_t = rng.integers(lo, hi, N_RANDOM)
    ctrl_ratio = S.y[rand_t] / S.prev[rand_t]
    ctrl_ratio = ctrl_ratio[np.isfinite(ctrl_ratio)]
    for (name, band), ev in events.items():
        t = np.array([e[0] for e in ev], int)
        dirs = np.array([e[1] for e in ev])
        since = np.array([e[2] for e in ev], int)
        base = np.array([e[3] for e in ev], float)
        post = S.y[t] / base if t.size else np.array([])
        # matched-delay control: random times, baseline from a 1440-minute window that ended as long before
        # as the alarm's calibration window did
        mc_up, mc_dn = [], []
        for d, dr in zip(t - since, dirs):
            if max(lo, d + DAY) >= hi:
                continue
            for tt in rng.integers(max(lo, d + DAY), hi, 5):
                b0 = window_sd(S.u, tt - d - DAY, tt - d)
                val = S.y[tt] / b0
                if np.isfinite(val):
                    (mc_up if dr == "up" else mc_dn).append(val)
        res = dict(source=source, period=period, detector=name, band=band, days=days, alarms=t.size,
                   alarms_per_day=t.size / days, up=int(np.sum(dirs == "up")), down=int(np.sum(dirs == "down")))
        for dname, cmp in (("up", lambda x: x > 1), ("down", lambda x: x < 1)):
            p = post[(dirs == dname) & np.isfinite(post)] if t.size else np.array([])
            k = int(np.sum(cmp(p)))
            share, l95, h95 = wilson(k, p.size)
            res[f"{dname}_validated_n"] = p.size
            res[f"{dname}_median_post_over_baseline"] = float(np.median(p)) if p.size else np.nan
            res[f"{dname}_share_in_flagged_direction"] = share
            res[f"{dname}_wilson_lo"] = l95
            res[f"{dname}_wilson_hi"] = h95
            mc = np.array(mc_up if dname == "up" else mc_dn)
            res[f"{dname}_matched_delay_control_share"] = float(np.mean(cmp(mc))) if mc.size else np.nan
        rows.append(res)
        for i in range(t.size):
            ev_rows.append(dict(source=source, period=period, detector=name, band=band, time=S.idx[t[i]], index=int(t[i]),
                                direction=dirs[i], delay_minutes=int(t[i] - since[i]), baseline_scale=base[i],
                                post_over_baseline=post[i]))
    rows.append(dict(source=source, period=period, detector="control: random times (next day / previous day)", band=np.nan,
                     days=days, alarms=ctrl_ratio.size, up_share_in_flagged_direction=float(np.mean(ctrl_ratio > 1)),
                     down_share_in_flagged_direction=float(np.mean(ctrl_ratio < 1)),
                     up_median_post_over_baseline=float(np.median(ctrl_ratio)), down_median_post_over_baseline=float(np.median(ctrl_ratio))))
    return rows, ev_rows


# ----------------------------------------------------------------------------
# Q6e: forecasts and the rate-matched competitor
# ----------------------------------------------------------------------------


def ratio_rule(S: Series, short: int, k_up: float, k_dn: float):
    """Alarm when short-window robust vol / previous 1440-minute robust vol leaves [k_dn, k_up]; dead time 1440."""
    rho = S.rho[short]
    with np.errstate(invalid="ignore"):
        cand = np.flatnonzero((rho >= k_up) | (rho <= k_dn))
    out, nxt = [], S.i0 + DAY  # like the monitor: the first day is calibration
    while True:
        j = np.searchsorted(cand, nxt)
        if j >= cand.size:
            break
        t = int(cand[j])
        out.append((t, "up" if rho[t] >= k_up else "down", S.v1[t - short]))
        nxt = t + DAY + 1
    return out


def match_rule(S: Series, short: int, n_up: int, n_dn: int):
    """Thresholds that reproduce the e-detector's numbers of up and down alarms (coordinate-wise bisection)."""
    k_up, k_dn = 2.0, 0.5
    for _ in range(4):
        lo, hi = 1.0, 20.0
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            nu = sum(1 for e in ratio_rule(S, short, mid, k_dn) if e[1] == "up")
            lo, hi = (mid, hi) if nu > n_up else (lo, mid)
        k_up = hi
        lo, hi = 0.01, 1.0
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            nd = sum(1 for e in ratio_rule(S, short, k_up, mid) if e[1] == "down")
            lo, hi = (lo, mid) if nd > n_dn else (mid, hi)
        k_dn = lo
    return k_up, k_dn


def forecast_hits(source, period, S: Series, events, rng):
    rows = []
    lo, hi = S.i0 + DAY, S.n - DAY - 1
    rand_t = rng.integers(lo, hi, N_RANDOM)
    fc = {"EWMA": S.forecast_ewma, "HAR": S.forecast_har}
    rand = {}
    for fname, f in fc.items():
        F = f(rand_t)
        y = S.y[rand_t]
        ok = np.isfinite(F) & np.isfinite(y)
        rand[fname] = (y[ok], F[ok])
        rows.append(dict(source=source, period=period, rule="random times", band=np.nan, short_window=np.nan, forecast=fname,
                         n_up=int(ok.sum()), up_hits=int(np.sum(y[ok] > F[ok])), n_down=int(ok.sum()),
                         down_hits=int(np.sum(y[ok] < F[ok])),
                         up_median_log_realised_over_forecast=float(np.median(np.log(y[ok] / F[ok]))),
                         down_median_log_realised_over_forecast=float(np.median(np.log(y[ok] / F[ok])))))

    def score(rule, band, short, ev, extra):
        t = np.array([e[0] for e in ev], int)
        d = np.array([e[1] for e in ev])
        for fname, f in fc.items():
            F = f(t) if t.size else np.array([])
            y = S.y[t] if t.size else np.array([])
            ok = np.isfinite(F) & np.isfinite(y)
            rec = dict(source=source, period=period, rule=rule, band=band, short_window=short, forecast=fname, **extra)
            for dn, cmp in (("up", np.greater), ("down", np.less)):
                m = ok & (d == dn)
                rec[f"n_{dn}"] = int(m.sum())
                rec[f"{dn}_hits"] = int(np.sum(cmp(y[m], F[m])))
                rec[f"{dn}_median_log_realised_over_forecast"] = float(np.median(np.log(y[m] / F[m]))) if m.any() else np.nan
            rows.append(rec)

    comp_rows = []
    for (name, band), ev in events.items():
        score(name, band, np.nan, ev, {})
        if not name.startswith("robust"):
            continue
        n_up = sum(1 for e in ev if e[1] == "up")
        n_dn = len(ev) - n_up
        for short in SHORT_WINDOWS:
            k_up, k_dn = match_rule(S, short, n_up, n_dn)
            cev = ratio_rule(S, short, k_up, k_dn)
            score(f"ratio rule matched to band {band}", band, short, cev, dict(k_up=k_up, k_down=k_dn))
            t = np.array([e[0] for e in cev], int)
            d = np.array([e[1] for e in cev])
            base = np.array([e[2] for e in cev])
            post = S.y[t] / base if t.size else np.array([])
            rec = dict(source=source, period=period, band=band, short_window=short, k_up=k_up, k_down=k_dn,
                       edetector_up=n_up, edetector_down=n_dn, rule_up=int(np.sum(d == "up")), rule_down=int(np.sum(d == "down")))
            for dn, cmp in (("up", np.greater), ("down", np.less)):
                p = post[(d == dn) & np.isfinite(post)] if t.size else np.array([])
                rec[f"{dn}_share_post_over_baseline_in_direction"] = float(np.mean(cmp(p, 1))) if p.size else np.nan
            # Do the two rules flag the same events? Share of e-detector alarms with a same-direction rule alarm
            # within 60 or 360 minutes, and the median of (e-detector time - rule time) for those within 360.
            te = np.array([e[0] for e in ev], int)
            de = np.array([e[1] for e in ev])
            for win in (60, 360):
                hit, lead = [], []
                for ti, di in zip(te, de):
                    m = (d == di) & (np.abs(t - ti) <= win) if t.size else np.array([], bool)
                    hit.append(bool(np.any(m)))
                    if win == 360 and np.any(m):
                        j = np.argmin(np.where(m, np.abs(t - ti), np.inf))
                        lead.append(ti - t[j])
                rec[f"share_edetector_alarms_with_rule_alarm_within_{win}"] = float(np.mean(hit)) if hit else np.nan
            rec["median_edetector_minus_rule_minutes"] = float(np.median(lead)) if lead else np.nan
            comp_rows.append(rec)
    out = pd.DataFrame(rows)
    for dn in ("up", "down"):
        w = [wilson(k, n) for k, n in zip(out[f"{dn}_hits"], out[f"n_{dn}"])]
        out[f"{dn}_hit_rate"] = [x[0] for x in w]
        out[f"{dn}_wilson_lo"] = [x[1] for x in w]
        out[f"{dn}_wilson_hi"] = [x[2] for x in w]
        # two-proportion z-test against random times with the same forecast (independence assumed)
        pz = []
        for _, r in out.iterrows():
            ref = out[(out.rule == "random times") & (out.forecast == r.forecast)].iloc[0]
            k1, n1, k2, n2 = r[f"{dn}_hits"], r[f"n_{dn}"], ref[f"{dn}_hits"], ref[f"n_{dn}"]
            if r.rule == "random times" or n1 == 0:
                pz.append(np.nan)
                continue
            p = (k1 + k2) / (n1 + n2)
            se = np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
            pz.append(float(2 * stats.norm.sf(abs(k1 / n1 - k2 / n2) / se)) if se > 0 else np.nan)
        out[f"{dn}_p_vs_random"] = pz
    return out, pd.DataFrame(comp_rows)


# ----------------------------------------------------------------------------
# Cross-venue and cross-coin coincidence
# ----------------------------------------------------------------------------


def coincidence(ta, da, tb, db, span_minutes, rng, n_shift=1000):
    """Share of alarms in A with an alarm in B within MATCH_MINUTES (any / same direction), and a shift null."""
    ta, tb = np.asarray(ta, float), np.asarray(tb, float)

    def share(tb_, same):
        if ta.size == 0 or tb_.size == 0:
            return np.nan
        hit = []
        for t, d in zip(ta, da):
            m = np.abs(tb_ - t) <= MATCH_MINUTES
            hit.append(bool(np.any(m & (np.asarray(db) == d))) if same else bool(np.any(m)))
        return float(np.mean(hit))

    obs_any, obs_same = share(tb, False), share(tb, True)
    null_any, null_same = [], []
    for _ in range(n_shift):
        off = rng.uniform(0, span_minutes)
        tb_ = (tb + off) % span_minutes
        null_any.append(share(tb_, False))
        null_same.append(share(tb_, True))
    null_any, null_same = np.array(null_any), np.array(null_same)
    return dict(share_any=obs_any, share_same_direction=obs_same, null_mean_any=float(null_any.mean()),
                null_mean_same=float(null_same.mean()), p_any=float((1 + np.sum(null_any >= obs_any)) / (1 + n_shift)),
                p_same=float((1 + np.sum(null_same >= obs_same)) / (1 + n_shift)))


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

UP_COLOR, DOWN_COLOR = "#eb6834", "#2a78d6"
BAR_COLORS = {"robust e-detector (band 0.25)": "#2a78d6", "ratio rule, 60 min": "#eb6834",
              "ratio rule, 240 min": "#1baf7a", "random times": "#9a9893"}


def fig_timelines(series, events):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    srcs = ["BTC", "Bitstamp BTC", "ETH", "XRP", "SOL", "DOGE"]
    fig, axes = plt.subplots(len(srcs), 1, figsize=(11, 12), sharex=True)
    for ax, src in zip(axes, srcs):
        S = series[(src, "Jul-Dec 2025")]
        idx = S.idx[S.i0 :]
        hv = pd.Series(S.v1[S.i0 :], index=idx).iloc[::60]
        ax.plot(hv.index, hv.values, lw=0.8, color="#52514e")
        ax.set_yscale("log")
        ax.set_ylabel(f"{src}\n1-day robust vol", fontsize=8)
        ax.grid(axis="y", color="#e6e5e1", lw=0.5)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ev = events[(src, "Jul-Dec 2025")][("robust e-detector", 0.25)]
        nu = nd = 0
        for t, d, _, _ in ev:
            col = UP_COLOR if d == "up" else DOWN_COLOR
            ax.axvline(S.idx[t], color=col, lw=0.6, alpha=0.7)
            nu += d == "up"
            nd += d == "down"
        ax.text(0.005, 0.92, f"{nu} up, {nd} down", transform=ax.transAxes, fontsize=8, va="top", color="#0b0b0b")
    axes[0].set_title("Robust e-detector alarms, band 0.25, frozen settings [REAL]. Orange: increase, blue: decrease.\n"
                      "Line: robust volatility of the previous 1440 deseasonalised returns (log scale). Binance unless noted.",
                      fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "figures/e06_alarm_timelines.png", dpi=130)
    plt.close(fig)


def fig_hits(hits):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h = hits[hits.forecast == "EWMA"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.4), sharey=True)
    coins = list(COINS)
    groups = [("robust e-detector (band 0.25)", lambda d: (d.rule == "robust e-detector") & (d.band == 0.25)),
              ("ratio rule, 60 min", lambda d: (d.rule == "ratio rule matched to band 0.25") & (d.short_window == 60)),
              ("ratio rule, 240 min", lambda d: (d.rule == "ratio rule matched to band 0.25") & (d.short_window == 240)),
              ("random times", lambda d: d.rule == "random times")]
    for i, period in enumerate(PERIODS):
        for j, dn in enumerate(("up", "down")):
            ax = axes[i, j]
            hp = h[h.period == period]
            x = np.arange(len(coins))
            wdt = 0.2
            for g, (lab, sel) in enumerate(groups):
                vals, lo, hi = [], [], []
                for c in coins:
                    r = hp[(hp.source == c) & sel(hp)]
                    if len(r):
                        vals.append(r[f"{dn}_hit_rate"].iloc[0])
                        lo.append(r[f"{dn}_wilson_lo"].iloc[0])
                        hi.append(r[f"{dn}_wilson_hi"].iloc[0])
                    else:
                        vals.append(np.nan)
                        lo.append(np.nan)
                        hi.append(np.nan)
                vals, lo, hi = map(np.array, (vals, lo, hi))
                xs = x + (g - 1.5) * (wdt + 0.02)
                ax.bar(xs, vals, width=wdt, color=BAR_COLORS[lab], label=lab, edgecolor="white", linewidth=1)
                ax.errorbar(xs, vals, yerr=[vals - lo, hi - vals], fmt="none", ecolor="#0b0b0b", elinewidth=0.8, capsize=2)
            ax.axhline(0.5, color="#9a9893", lw=0.8, ls="--")
            ax.set_xticks(x)
            ax.set_xticklabels(coins)
            ax.set_ylim(0, 1)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            word = "above" if dn == "up" else "below"
            per = "Jul to Dec 2025" if period.startswith("Jul") else "Jan 2024 to Jun 2025 (out of sample)"
            ax.set_title(f"{per}\nafter '{dn}' alarms: next-day vol {word} EWMA forecast", fontsize=9, loc="left")
            if j == 0:
                ax.set_ylabel("share of alarms (Wilson 95%)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("Do alarms beat a volatility forecast? [REAL] Five Binance coins; ratio rules are matched on alarm counts",
                 fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "figures/e06_hit_rates.png", dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------


def pooled_monitor_summary():
    """Pool the E5c-style validation over the five Binance coins, per period, detector and band."""
    m = pd.read_csv(OUT / "e06d_monitor_summary.csv")
    m = m[m.source.isin(list(COINS)) & ~m.detector.str.startswith("control")]
    rows = []
    for (period, det, band), g in m.groupby(["period", "detector", "band"], sort=False):
        rec = dict(period=period, detector=det, band=band, coins=len(g), alarms_per_day_min=g.alarms_per_day.min(),
                   alarms_per_day_max=g.alarms_per_day.max(), up=int(g.up.sum()), down=int(g.down.sum()))
        for dn in ("up", "down"):
            n = int(g[f"{dn}_validated_n"].sum())
            k = int(np.round(g[f"{dn}_validated_n"] * g[f"{dn}_share_in_flagged_direction"]).sum())
            p_, lo, hi = wilson(k, n)
            rec.update({f"{dn}_n": n, f"{dn}_share_in_flagged_direction": p_, f"{dn}_wilson_lo": lo, f"{dn}_wilson_hi": hi,
                        f"{dn}_matched_delay_control_mean": g[f"{dn}_matched_delay_control_share"].mean(),
                        f"{dn}_median_post_over_baseline_range": f"{g[f'{dn}_median_post_over_baseline'].min():.2f} to "
                                                                f"{g[f'{dn}_median_post_over_baseline'].max():.2f}"})
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "e06d_pooled.csv", index=False)
    print(out.round(3).to_string(index=False))
    return out


def pooled_summary():
    """Pool the forecast hit counts over the five Binance coins, per period, rule and forecast."""
    h = pd.read_csv(OUT / "e06e_forecast_hits.csv")
    h = h[h.source.isin(list(COINS))]
    h["short_window"] = h.short_window.fillna(0).astype(int)
    h["band"] = h.band.fillna(-1)
    g = h.groupby(["period", "forecast", "rule", "band", "short_window"], sort=False)[["up_hits", "n_up", "down_hits", "n_down"]].sum()
    g = g.reset_index()
    for dn in ("up", "down"):
        w = [wilson(k, n) for k, n in zip(g[f"{dn}_hits"], g[f"n_{dn}"])]
        g[f"{dn}_hit_rate"] = [x[0] for x in w]
        g[f"{dn}_wilson_lo"] = [x[1] for x in w]
        g[f"{dn}_wilson_hi"] = [x[2] for x in w]
        # difference to the e-detector with the same band, period and forecast (two-proportion z-test)
        diff, pv = [], []
        for _, r in g.iterrows():
            ref = g[(g.period == r.period) & (g.forecast == r.forecast) & (g.rule == "robust e-detector") & (g.band == r.band)]
            if r.rule == "robust e-detector" or ref.empty or r[f"n_{dn}"] == 0:
                diff.append(np.nan)
                pv.append(np.nan)
                continue
            ref = ref.iloc[0]
            k1, n1, k2, n2 = ref[f"{dn}_hits"], ref[f"n_{dn}"], r[f"{dn}_hits"], r[f"n_{dn}"]
            pp = (k1 + k2) / (n1 + n2)
            se = np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
            diff.append(k1 / n1 - k2 / n2)
            pv.append(float(2 * stats.norm.sf(abs(k1 / n1 - k2 / n2) / se)) if se > 0 else np.nan)
        g[f"{dn}_edetector_minus_this"] = diff
        g[f"{dn}_p_edetector_vs_this"] = pv
    g["band"] = g.band.replace(-1, np.nan)
    g.to_csv(OUT / "e06e_pooled.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["period", "forecast", "rule", "band", "short_window", "n_up", "up_hit_rate", "up_wilson_lo", "up_wilson_hi",
            "up_edetector_minus_this", "up_p_edetector_vs_this", "n_down", "down_hit_rate", "down_wilson_lo", "down_wilson_hi",
            "down_edetector_minus_this", "down_p_edetector_vs_this"]
    print(g[cols].round(3).to_string(index=False))
    return g


def main():
    if "--summaries-only" in sys.argv:  # rebuild the pooled table and the hit-rate figure from saved results
        pooled_monitor_summary()
        pooled_summary()
        fig_hits(pd.read_csv(OUT / "e06e_forecast_hits.csv"))
        return
    jobs = [(c, p) for p in PERIODS for c in COINS] + [("Bitstamp BTC", "Jul-Dec 2025")]
    with Pool(2) as pool:
        res = pool.map(run_job, jobs)
    events = {(s, p): o for s, p, o in res}
    rng = np.random.default_rng(6)
    series, idx_only, summ, evs, hits, comps = {}, {}, [], [], [], []
    for (src, period) in jobs:
        idx, u, i0 = prepare(src, period)
        S = Series(idx, u, i0)
        r, e = summarise(src, period, S, events[(src, period)], rng)
        summ += r
        evs += e
        h, cp = forecast_hits(src, period, S, events[(src, period)], rng)
        hits.append(h)
        comps.append(cp)
        print(f"done {src} {period}: {len(u) - i0} minutes, NaN returns {int(np.sum(~np.isfinite(u[i0:])))}", flush=True)
        if period == "Jul-Dec 2025":
            series[(src, period)] = S
        else:
            idx_only[(src, period)] = idx
        del S
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "e06d_monitor_summary.csv", index=False)
    pd.DataFrame(evs).to_csv(OUT / "e06d_events.csv", index=False)
    hits = pd.concat(hits, ignore_index=True)
    hits.to_csv(OUT / "e06e_forecast_hits.csv", index=False)
    comps = pd.concat(comps, ignore_index=True)
    comps.to_csv(OUT / "e06e_competitor.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["source", "period", "detector", "band", "alarms_per_day", "up", "down", "up_share_in_flagged_direction",
            "up_wilson_lo", "up_wilson_hi", "up_matched_delay_control_share", "down_share_in_flagged_direction",
            "down_wilson_lo", "down_wilson_hi", "down_matched_delay_control_share"]
    print(summ[cols].round(3).to_string(index=False))
    hc = ["source", "period", "rule", "band", "short_window", "forecast", "n_up", "up_hit_rate", "up_p_vs_random", "n_down",
          "down_hit_rate", "down_p_vs_random"]
    print(hits[hc].round(3).to_string(index=False))
    print(comps.round(3).to_string(index=False))

    # Cross-venue: Binance BTCUSDT against Bitstamp BTC/USD, Jul to Dec 2025
    cv = []
    rng2 = np.random.default_rng(66)
    Sb, Ss = series[("BTC", "Jul-Dec 2025")], series[("Bitstamp BTC", "Jul-Dec 2025")]
    t0 = pd.Timestamp(PERIODS["Jul-Dec 2025"][0], tz="UTC")
    span = (pd.Timestamp(PERIODS["Jul-Dec 2025"][1], tz="UTC") - t0) / pd.Timedelta(minutes=1)
    to_min = lambda S, ev: [(S.idx[t] - t0) / pd.Timedelta(minutes=1) for t, *_ in ev]
    for cfg in CONFIGS:
        eb, es = events[("BTC", "Jul-Dec 2025")][cfg], events[("Bitstamp BTC", "Jul-Dec 2025")][cfg]
        for a_name, (Sa, ea), b_name, (Sb_, eb_) in (("Binance", (Sb, eb), "Bitstamp", (Ss, es)), ("Bitstamp", (Ss, es), "Binance", (Sb, eb))):
            r = coincidence(to_min(Sa, ea), [e[1] for e in ea], to_min(Sb_, eb_), [e[1] for e in eb_], span, rng2)
            cv.append(dict(detector=cfg[0], band=cfg[1], alarms_of=a_name, matched_against=b_name, n_alarms=len(ea),
                           n_other=len(eb_), window_minutes=MATCH_MINUTES, **r))
    cv = pd.DataFrame(cv)
    cv.to_csv(OUT / "e06d_cross_venue.csv", index=False)
    print(cv.round(3).to_string(index=False))
    # Cross-coin: share of each coin's alarms within 60 minutes of a BTC alarm (band 0.25), both periods
    cc = []
    for period in PERIODS:
        p0 = pd.Timestamp(PERIODS[period][0], tz="UTC")
        sp = (pd.Timestamp(PERIODS[period][1], tz="UTC") - p0) / pd.Timedelta(minutes=1)
        ix = {c: (series[(c, period)].idx if (c, period) in series else idx_only[(c, period)]) for c in COINS}
        tm = lambda c, ev: [(ix[c][t] - p0) / pd.Timedelta(minutes=1) for t, *_ in ev]
        eb = events[("BTC", period)][("robust e-detector", 0.25)]
        for c in ("ETH", "XRP", "SOL", "DOGE"):
            ec = events[(c, period)][("robust e-detector", 0.25)]
            r = coincidence(tm(c, ec), [e[1] for e in ec], tm("BTC", eb),
                            [e[1] for e in eb], sp, rng2, n_shift=500)
            cc.append(dict(period=period, coin=c, n_alarms=len(ec), n_btc=len(eb), **r))
    cc = pd.DataFrame(cc)
    cc.to_csv(OUT / "e06d_cross_coin.csv", index=False)
    print(cc.round(3).to_string(index=False))

    pooled_monitor_summary()
    pooled_summary()
    fig_timelines(series, events)
    fig_hits(hits)


if __name__ == "__main__":
    sys.exit(main())
