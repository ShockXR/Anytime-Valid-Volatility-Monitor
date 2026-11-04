"""E20. A re-analysis of the real-data inference: Holm within families, a block bootstrap over calendar days,
and forecast evaluation with realised variance.

REAL MARKET DATA [REAL]: Binance spot one-minute bars of BTC, ETH, XRP, SOL and DOGE against USDT
(``avmonitor.data.load_binance``, download=False) and the Nifty 50 one-minute bars of E7 (``e07_nifty.load_clean``).

The data and the alarms analysed here had been seen before (E6, E7, E12, E13, E14, E15, E17), so E20 is a
re-analysis whose methods were fixed in advance of their own results, not a confirmatory study. The methods are
in docs/preregistration_e20.md, committed (bbc8f87) before any adjusted result of this script was computed.
Nothing here changes a monitor setting, and no committed experiment script is modified: E6's and E17's
functions are imported and called as they are.

Parts (``python experiments/e20_inference.py <part> ...``)
  checks    Reload or rerun the alarm lists and check them against the committed results (pre-registration 3.2).
            E6: ``summarise`` and ``forecast_hits`` are called in E6's job order with E6's seed, so the random times,
            the matched-delay control draws and the ratio-rule alarms are regenerated exactly; E17: the matched
            one-minute monitor is rerun at its stored levels. Writes results/e20_checks.csv and the per-item files
            results/e20_items_e6.csv.gz and results/e20_items_e17.csv.gz (one row per scored alarm, random time or
            control draw).
  bootstrap Families F1 to F4, F6 to F9, F14 to F16: circular block bootstrap over calendar days (7 days) or Nifty
            sessions (5 sessions), days resampled jointly across coins and methods, B = 9999, seed 20; percentile
            intervals and p-values, Holm within families. Sensitivity with 14 days and 10 sessions.
            -> results/e20_bootstrap.csv
  holm      Families F5, F10a, F10c, F11, F12a to F12d, F13: the published p-values with Holm within families.
            -> results/e20_holm_existing.csv
  events    F10b: circular-shift test of the Nifty event calendar against the alarm series.
            -> results/e20_nifty_event_shift.csv
  forecast  Part 3: next-day and next 6-hour log realised variance (and bipower variation) regressed on the log
            forecast (EWMA, HAR, HAR with the last hour) plus up- and down-alarm indicators, Driscoll-Kraay standard
            errors by calendar day; QLIKE of each forecast against the forecast scaled at alarms.
            -> results/e20_forecast_choice.csv, e20_forecast_regressions.csv, e20_forecast_qlike.csv,
               e20_forecast_alarm_counts.csv
  classify  The classification of every real-data test of the paper (paper at commit 8f71a2e), its adjusted status,
            and the primary endpoint of each data set. -> results/e20_classification.csv,
            results/e20_primary_endpoints.csv

Runtime: checks about 2 minutes on one core (E6 series, E17 reruns), bootstrap about 1 minute, forecast about
1 minute. Run with ``nice -n 10``.
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
sys.path.insert(0, str(ROOT / "experiments"))

DAY = 1440
B_BOOT = 9999
SEED = 20
BLOCK_CRYPTO, BLOCK_NIFTY = 7, 5
DK_BANDWIDTH = 5
EQUIV_MARGIN = 0.05
COINS = ["BTC", "ETH", "XRP", "SOL", "DOGE"]
P1, P2 = "Jul-Dec 2025", "Jan 2024-Jun 2025 (out of sample)"
PERIOD_DAYS = {P1: ("2025-07-01", "2026-01-01"), P2: ("2024-01-21", "2025-07-01"), "Jan-Jun 2025": ("2025-01-01", "2025-07-01")}


# ----------------------------------------------------------------------------
# Helpers (tested in tests/test_e20_inference.py)
# ----------------------------------------------------------------------------


def holm(p) -> np.ndarray:
    """Holm step-down adjusted p-values; NaN entries are left out of the family and stay NaN."""
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    m = ok.size
    running = 0.0
    for rank, i in enumerate(ok[np.argsort(p[ok], kind="stable")]):
        running = max(running, min(1.0, (m - rank) * p[i]))
        out[i] = running
    return out


def cbb_weights(D: int, L: int, B: int, seed: int) -> np.ndarray:
    """Circular block bootstrap: B x D matrix of how often each of D consecutive days is drawn."""
    rng = np.random.default_rng(seed)
    nb = -(-D // L)
    starts = rng.integers(0, D, size=(B, nb))
    idx = ((starts[:, :, None] + np.arange(L)) % D).reshape(B, nb * L)[:, :D]
    flat = (np.arange(B)[:, None] * D + idx).ravel()
    return np.bincount(flat, minlength=B * D).reshape(B, D).astype(float)


def boot_p(delta_star: np.ndarray) -> float:
    """Two-sided p-value for H0: delta = 0 by inverting the percentile interval."""
    d = np.asarray(delta_star, float)
    d = d[np.isfinite(d)]
    b = d.size
    return float(min(1.0, 2.0 * min(1 + np.sum(d <= 0), 1 + np.sum(d >= 0)) / (b + 1)))


def ols_dk(y, X, day, bandwidth: int):
    """OLS with Driscoll-Kraay standard errors: scores summed within each calendar day (all units), Bartlett
    weights up to ``bandwidth`` days, factor D / (D - 1). ``day`` is an integer day number; days without
    observations count as zero scores. Returns (beta, se, number of days)."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    e = y - X @ beta
    g = X * e[:, None]
    d0 = int(np.min(day))
    dd = np.asarray(day, int) - d0
    D = int(dd.max()) + 1
    h = np.zeros((D, X.shape[1]))
    np.add.at(h, dd, g)
    S = h.T @ h
    for lag in range(1, bandwidth + 1):
        w = 1.0 - lag / (bandwidth + 1.0)
        G = h[lag:].T @ h[:-lag]
        S += w * (G + G.T)
    n_days = int(np.unique(dd).size)
    V = XtX_inv @ S @ XtX_inv * n_days / max(n_days - 1, 1)
    return beta, np.sqrt(np.diag(V)), n_days


def qlike(y, F):
    """QLIKE loss y / F - log(y / F) - 1 (zero at y = F)."""
    q = np.asarray(y, float) / np.asarray(F, float)
    return q - np.log(q) - 1.0


def qlike_kappa(y, F) -> float:
    """The constant k minimising the mean QLIKE of k F against y: the mean of y / F."""
    return float(np.mean(np.asarray(y, float) / np.asarray(F, float)))


def norm_p(z) -> float:
    return float(math.erfc(abs(z) / math.sqrt(2.0)))


# ----------------------------------------------------------------------------
# Part: checks (alarm lists, regenerated controls) for E6 and E17
# ----------------------------------------------------------------------------


def _close(a, b, rel=1e-12):
    a, b = np.asarray(a, float), np.asarray(b, float)
    both_nan = np.isnan(a) & np.isnan(b)
    return bool(np.all(both_nan | (np.abs(a - b) <= rel * np.maximum(1.0, np.maximum(np.abs(a), np.abs(b))))))


def _compare_frames(new: pd.DataFrame, old: pd.DataFrame, keys: list[str], label: str, checks: list, rel=1e-12):
    """Compare numeric columns of two frames row by row after sorting on keys; record the result."""
    new = new.copy()
    old = old.copy()
    for k in keys:  # map(str), not astype(str): pandas 3 keeps NaN as missing under astype(str)
        new[k] = new[k].map(str)
        old[k] = old[k].map(str)
    new = new.sort_values(keys).reset_index(drop=True)
    old = old.sort_values(keys).reset_index(drop=True)
    ok = len(new) == len(old) and bool((new[keys].to_numpy() == old[keys].to_numpy()).all())
    bad = []
    if ok:
        for c in old.columns:
            if c in keys or c not in new.columns:
                continue
            if pd.api.types.is_numeric_dtype(old[c]):
                if not _close(pd.to_numeric(new[c], errors="coerce"), old[c], rel):
                    bad.append(c)
    checks.append(dict(check=label, rows=len(old), passed=ok and not bad, detail="; ".join(bad) if bad else ""))
    return ok and not bad


def _mirror_controls(S, events, rng):
    """Mirror of the random draws of e06b ``summarise``: the random times and, per alarm, the 5 matched-delay
    control draws, with their times (``summarise`` does not return them). Same generator calls in the same order."""
    from e06b_five_coins_monitor import window_sd

    lo, hi = S.i0 + DAY, S.n - DAY - 1
    rng.integers(lo, hi, 3000)  # the random times of summarise (E5c control, not used by E20)
    out = []
    for (name, band), ev in events.items():
        t = np.array([e[0] for e in ev], int)
        dirs = np.array([e[1] for e in ev])
        since = np.array([e[2] for e in ev], int)
        for d, dr in zip(t - since, dirs):
            if max(lo, d + DAY) >= hi:
                continue
            for tt in rng.integers(max(lo, d + DAY), hi, 5):
                b0 = window_sd(S.u, tt - d - DAY, tt - d)
                val = S.y[tt] / b0
                if np.isfinite(val):
                    out.append((name, band, dr, int(tt), float(val)))
    return out


def part_checks():
    import e06b_five_coins_monitor as E6

    checks = []
    ev_all = pd.read_csv(OUT / "e06d_events.csv")
    summ_old = pd.read_csv(OUT / "e06d_monitor_summary.csv")
    hits_old = pd.read_csv(OUT / "e06e_forecast_hits.csv")
    comp_old = pd.read_csv(OUT / "e06e_competitor.csv")
    rm_levels = pd.read_csv(OUT / "e17_ratematch_levels.csv")
    e17_ev = pd.read_csv(OUT / "e17_real_events.csv")
    rng = np.random.default_rng(6)  # E6's generator, used in E6's job order
    items, items17 = [], []
    for period in E6.PERIODS:
        for coin in E6.COINS:
            t0 = time.time()
            idx, u, i0 = E6.prepare(coin, period)
            S = E6.Series(idx, u, i0)
            pstart = pd.Timestamp(PERIOD_DAYS[period][0], tz="UTC")
            dayof = lambda tt: ((S.idx[np.asarray(tt, int)] - pstart).days).to_numpy() if np.size(tt) else np.array([], int)
            # events in the order of E6's CONFIGS, alarms in time order, as run_job returned them
            events = {}
            for name, band in E6.CONFIGS:
                e = ev_all[(ev_all.source == coin) & (ev_all.period == period) & (ev_all.detector == name) & (ev_all.band == band)]
                events[(name, band)] = [(int(r.index), r.direction, int(r.index) - int(r.delay_minutes), float(r.baseline_scale))
                                        for r in e.itertuples()]
            # summarise, with a mirror of its random draws from the same generator state
            state = rng.bit_generator.state
            rows, ev_rows = E6.summarise(coin, period, S, events, rng)
            mirror = np.random.default_rng()
            mirror.bit_generator.state = state
            ctrl = _mirror_controls(S, events, mirror)
            if mirror.bit_generator.state != rng.bit_generator.state:
                checks.append(dict(check=f"E6 {coin} {period}: mirror of summarise draws", rows=0, passed=False, detail="generator state differs"))
            new = pd.DataFrame(rows)
            old = summ_old[(summ_old.source == coin) & (summ_old.period == period)]
            _compare_frames(new, old, ["detector", "band"], f"E6 {coin} {period}: e06d_monitor_summary rows", checks)
            # mirror against the summary's matched-delay control shares
            ok = True
            for (name, band) in events:
                for dn, cmp in (("up", lambda x: x > 1), ("down", lambda x: x < 1)):
                    vals = np.array([c[4] for c in ctrl if c[0] == name and c[1] == band and c[2] == dn])
                    share = float(np.mean(cmp(vals))) if vals.size else np.nan
                    ref = new[(new.detector == name) & (new.band == band)][f"{dn}_matched_delay_control_share"].iloc[0]
                    ok &= _close(share, ref)
            checks.append(dict(check=f"E6 {coin} {period}: mirrored matched-delay control draws", rows=len(ctrl), passed=ok, detail=""))
            # forecast_hits (random times, e-detector and ratio-rule hits), random times regenerated
            state = rng.bit_generator.state
            hits, comps = E6.forecast_hits(coin, period, S, events, rng)
            mirror = np.random.default_rng()
            mirror.bit_generator.state = state
            lo, hi = S.i0 + DAY, S.n - DAY - 1
            rand_t = mirror.integers(lo, hi, E6.N_RANDOM)
            hold = hits_old[(hits_old.source == coin) & (hits_old.period == period)]
            hits_k = hits.copy()
            hits_k["short_window"] = hits_k.short_window.fillna(0)
            hold = hold.copy()
            hold["short_window"] = hold.short_window.fillna(0)
            hits_k["band"] = hits_k.band.fillna(-1)
            hold["band"] = hold.band.fillna(-1)
            _compare_frames(hits_k, hold, ["rule", "band", "short_window", "forecast"], f"E6 {coin} {period}: e06e_forecast_hits rows", checks)
            cold = comp_old[(comp_old.source == coin) & (comp_old.period == period)]
            _compare_frames(comps, cold, ["band", "short_window"], f"E6 {coin} {period}: e06e_competitor rows (thresholds and counts)", checks)
            # per-item outcomes
            fe = {"ewma": S.forecast_ewma, "har": S.forecast_har}

            def add(kind, rule, band, short, tt, dirs, post):
                tt = np.asarray(tt, int)
                if tt.size == 0:
                    return
                y = S.y[tt]
                rec = dict(period=period, coin=coin, kind=kind, rule=rule, band=band, short=short, index=tt,
                           day=dayof(tt), direction=dirs, post_over_baseline=post)
                for k, f in fe.items():
                    rec[f"y_over_{k}"] = y / f(tt)
                items.append(pd.DataFrame(rec))

            add("random", "random times", np.nan, 0, rand_t, "none", np.nan)
            for (name, band), ev in events.items():
                tt = np.array([e[0] for e in ev], int)
                base = np.array([e[3] for e in ev], float)
                add("alarm", name, band, 0, tt, np.array([e[1] for e in ev]), S.y[tt] / base if tt.size else [])
            for r in comps.itertuples():
                cev = E6.ratio_rule(S, int(r.short_window), r.k_up, r.k_down)
                tt = np.array([e[0] for e in cev], int)
                add("alarm", "ratio rule", r.band, int(r.short_window), tt, np.array([e[1] for e in cev]), np.nan)
            c25 = [c for c in ctrl if c[0] == "robust e-detector" and c[1] == 0.25]
            if c25:
                tt = np.array([c[3] for c in c25], int)
                items.append(pd.DataFrame(dict(period=period, coin=coin, kind="control", rule="matched-delay control",
                                               band=0.25, short=0, index=tt, day=dayof(tt), direction=[c[2] for c in c25],
                                               post_over_baseline=[c[4] for c in c25], y_over_ewma=np.nan, y_over_har=np.nan)))
            # E17 Q17g: coarse detectors (stored) and the matched one-minute monitor (rerun at its stored levels)
            if period == P1:
                items17.extend(_e17_items(coin, S, i0, u, rm_levels, e17_ev, dayof, checks))
            print(f"checks {coin} {period}: {time.time() - t0:.0f} s", flush=True)
            del S
    it = pd.concat(items, ignore_index=True)
    it.to_csv(OUT / "e20_items_e6.csv.gz", index=False, float_format="%.10g")
    pd.DataFrame(items17).to_csv(OUT / "e20_items_e17.csv.gz", index=False, float_format="%.10g")
    checks += _other_checks(it, pd.DataFrame(items17))
    ch = pd.DataFrame(checks)
    ch.to_csv(OUT / "e20_checks.csv", index=False)
    print(ch.to_string())
    print("all checks passed" if ch.passed.all() else "SOME CHECKS FAILED")


def _e17_items(coin, S, i0, u, rm_levels, e17_ev, dayof, checks):
    from e17_multifrequency import _one_minute_events

    out = []
    for target in ("5 min", "15 min"):
        lv = rm_levels[(rm_levels.source == coin) & (rm_levels.target == target)].iloc[0]
        ev1 = _one_minute_events(u[i0:], float(lv.A_up), float(lv.A_down))
        nu = sum(e.direction == "up" for e in ev1)
        nd = len(ev1) - nu
        checks.append(dict(check=f"E17 {coin}: one-minute monitor matched to {target}, counts at stored levels",
                           rows=len(ev1), passed=(nu == lv.matched_up and nd == lv.matched_down),
                           detail=f"{nu} / {nd} against {lv.matched_up} / {lv.matched_down}"))
        cev = e17_ev[(e17_ev.source == coin) & (e17_ev.monitor == target)]
        nu_c = int((cev.direction == "up").sum())
        checks.append(dict(check=f"E17 {coin}: stored {target} detector counts against the rate-matching targets",
                           rows=len(cev), passed=(nu_c == lv.target_up and len(cev) - nu_c == lv.target_down),
                           detail=f"{nu_c} / {len(cev) - nu_c} against {lv.target_up} / {lv.target_down}"))
        for mon, tt, dirs, base in ((f"1 min matched to {target}", np.array([e.index + i0 for e in ev1], int),
                                     [e.direction for e in ev1], np.array([e.baseline_scale for e in ev1], float)),
                                    (target, cev["index"].to_numpy(int), cev.direction.tolist(), cev.baseline_scale.to_numpy(float))):
            y = S.y[tt]
            fe, fh = S.forecast_ewma(tt), S.forecast_har(tt)
            for j in range(tt.size):
                out.append(dict(coin=coin, comparison=target, monitor=mon, index=int(tt[j]), day=int(dayof(tt[j:j + 1])[0]),
                                direction=dirs[j], post_over_baseline=y[j] / base[j], y_over_ewma=y[j] / fe[j], y_over_har=y[j] / fh[j]))
    return out


def _hits(df, col, thresh=1.0):
    """Hit indicator for the flagged direction (NaN when the outcome cannot be scored)."""
    v = df[col].to_numpy(float)
    up = df.direction.to_numpy() == "up"
    h = np.where(up, v > thresh, v < thresh).astype(float)
    h[~np.isfinite(v)] = np.nan
    return h


def _other_checks(it, it17):
    """Checks of the stored alarm lists of E7, E13, E14, E15 and of the E17 items against the committed summaries."""
    checks = []
    # E6 pooled shares (Q6d) from the items, against e06d_pooled
    pooled = pd.read_csv(OUT / "e06d_pooled.csv")
    for r in pooled.itertuples():
        a = it[(it.kind == "alarm") & (it.rule == r.detector) & (it.band == r.band) & (it.period == r.period)]
        ok = True
        for dn in ("up", "down"):
            h = _hits(a[a.direction == dn], "post_over_baseline")
            h = h[np.isfinite(h)]
            ok &= h.size == getattr(r, f"{dn}_n") and _close(h.mean(), getattr(r, f"{dn}_share_in_flagged_direction"))
        checks.append(dict(check=f"E6 pooled Q6d {r.period} {r.detector} {r.band}", rows=len(a), passed=ok, detail=""))
    # E17 pooled
    p17 = pd.read_csv(OUT / "e17_ratematch_pooled.csv")
    for r in p17[~p17.monitor.str.startswith("p-value")].itertuples():
        a = it17[(it17.comparison == r.comparison) & (it17.monitor == r.monitor)]
        ok = True
        for dn in ("up", "down"):
            ad = a[a.direction == dn]
            for col, key in (("post_over_baseline", "confirmed_share"), ("y_over_ewma", "ewma_hit_rate"), ("y_over_har", "har_hit_rate")):
                h = _hits(ad, col)
                ok &= _close(np.nanmean(h), getattr(r, f"{dn}_{key}"), 1e-9)
        checks.append(dict(check=f"E17 pooled {r.comparison}: {r.monitor}", rows=len(a), passed=ok, detail=""))
    # E13
    lp = pd.read_csv(OUT / "e13d_real_long_windows_pooled.csv")
    al = pd.concat([pd.read_csv(OUT / "e13d_real_long_windows_crypto_alarms.csv"), pd.read_csv(OUT / "e13d_real_long_windows_nifty_alarms.csv")])
    ok_all = True
    for r in lp.itertuples():
        a = al[(al.market == r.market) & (al.method == r.method)]
        for dn in ("up", "down"):
            h = a[a.direction == dn].hit_next_day.dropna()
            ok_all &= h.size == getattr(r, f"n_{dn}_day") and _close(h.mean(), getattr(r, f"{dn}_rate_day"))
    checks.append(dict(check="E13 long-window alarm lists against e13d_real_long_windows_pooled", rows=len(lp), passed=ok_all, detail=""))
    # E7 Q7b
    a7 = pd.read_csv(OUT / "e07b_alarms.csv")
    s7 = pd.read_csv(OUT / "e07b_monitor_summary.csv")
    ok_all = True
    for r in s7.itertuples():
        a = a7[(a7.detector == r.detector) & (a7.cal_len == r.cal_len) & (a7.band == r.band) & (a7.period == r.period)]
        ok_all &= len(a) == r.alarms
        for dn in ("up", "down"):
            ad = a[a.direction == dn]
            h = np.where(dn == "up", ad.next_vol > ad.baseline, ad.next_vol < ad.baseline).astype(float)
            h[~np.isfinite(ad.next_vol.to_numpy())] = np.nan
            ok_all &= _close(np.nanmean(h), getattr(r, f"{dn}_next_vs_baseline"))
    checks.append(dict(check="E7 alarm list against e07b_monitor_summary", rows=len(s7), passed=ok_all, detail=""))
    # E14 Q14d
    a14 = pd.read_csv(OUT / "e14d_real_alarms.csv")
    s14 = pd.read_csv(OUT / "e14d_real_summary.csv")
    ok_all = True
    for r in s14.itertuples():
        if "alarms" not in s14.columns:
            break
        a = a14[(a14.source == r.source) & (a14.method == r.method) & (a14.band == r.band)]
        ok_all &= len(a) == r.alarms
    checks.append(dict(check="E14 alarm list counts against e14d_real_summary", rows=len(s14), passed=bool(ok_all), detail=""))
    # E15
    d15 = pd.read_csv(OUT / "e15d_discoveries.csv.gz")
    t15 = pd.read_csv(OUT / "e15d_rate_matched_tests.csv")
    ok_all = True
    for r in t15.itertuples():
        g = d15[(d15.panel == r.panel) & (d15.period == r.period) & (d15.band == r.band) & (d15.direction == r.direction)]
        for rule, n, hit in (("e-BH", r.ebh_n, r.ebh_hit), ("rate-matched threshold (hindsight)", r.matched_n, r.matched_hit)):
            h = _hits(g[g.rule == rule], "post_ratio")
            h = h[np.isfinite(h)]
            ok_all &= h.size == n and _close(h.mean(), hit)
    checks.append(dict(check="E15 discoveries against e15d_rate_matched_tests", rows=len(t15), passed=bool(ok_all), detail=""))
    return checks


# ----------------------------------------------------------------------------
# Part: bootstrap families
# ----------------------------------------------------------------------------


class Boot:
    """Circular block bootstrap over the D days of one period; one weight matrix shared by every statistic."""

    def __init__(self, D: int, L: int):
        self.D, self.L = D, L
        self.W = cbb_weights(D, L, B_BOOT, SEED)

    def share(self, day, hit):
        day = np.asarray(day, int)
        hit = np.asarray(hit, float)
        ok = np.isfinite(hit)
        day, hit = day[ok], hit[ok]
        if np.any((day < 0) | (day >= self.D)):
            raise ValueError("day outside the period")
        n = np.bincount(day, minlength=self.D).astype(float)
        h = np.bincount(day, weights=hit, minlength=self.D)
        with np.errstate(invalid="ignore", divide="ignore"):
            star = (self.W @ h) / (self.W @ n)
        return int(n.sum()), float(h.sum() / n.sum()) if n.sum() else np.nan, star


def _row(family, test, period, direction, measure, name_a, name_b, bt: Boot, a, b=None, fixed_b=None, orig=None):
    """One bootstrap row: share of a, and if given the difference a - b (b another item set or a fixed number)."""
    na, sa, sa_star = bt.share(*a)
    rec = dict(family=family, test=test, period=period, direction=direction, measure=measure, a=name_a, b=name_b,
               block=bt.L, days=bt.D, n_a=na, share_a=sa, share_a_lo=np.nanpercentile(sa_star, 2.5),
               share_a_hi=np.nanpercentile(sa_star, 97.5))
    if b is not None:
        nb, sb, sb_star = bt.share(*b)
        dstar = sa_star - sb_star
        rec.update(n_b=nb, share_b=sb, share_b_lo=np.nanpercentile(sb_star, 2.5), share_b_hi=np.nanpercentile(sb_star, 97.5))
    elif fixed_b is not None:
        sb = fixed_b
        dstar = sa_star - sb
        rec.update(n_b=np.nan, share_b=sb)
    else:
        dstar = None
    if dstar is not None:
        rec.update(diff=sa - sb, diff_se=float(np.nanstd(dstar, ddof=1)), diff_lo=np.nanpercentile(dstar, 2.5),
                   diff_hi=np.nanpercentile(dstar, 97.5), p_boot=boot_p(dstar))
    rec.update(orig or {})
    return rec


def _items_e6():
    it = pd.read_csv(OUT / "e20_items_e6.csv.gz")
    it["short"] = it["short"].fillna(0).astype(int)
    return it


def _days(period):
    a, b = PERIOD_DAYS[period]
    return (pd.Timestamp(b) - pd.Timestamp(a)).days


def families_crypto_e6(L):
    it = _items_e6()
    pooled = pd.read_csv(OUT / "e06e_pooled.csv")
    dpool = pd.read_csv(OUT / "e06d_pooled.csv")
    pooled["short_window"] = pooled.short_window.fillna(0).astype(int)
    rows = []
    for period in (P1, P2):
        bt = Boot(_days(period), L)
        ip = it[it.period == period]

        def sel(rule, band=None, short=0, kind="alarm", direction=None):
            m = (ip.kind == kind) & (ip.rule == rule) & (ip["short"] == short)
            if band is not None:
                m &= np.isclose(ip.band, band)
            if direction is not None:
                m &= ip.direction == direction
            return ip[m]

        # F1: Q6d, e-detector band 0.25 against the matched-delay control
        for dn in ("up", "down"):
            a = sel("robust e-detector", 0.25, direction=dn)
            c = sel("matched-delay control", 0.25, kind="control", direction=dn)
            dp = dpool[(dpool.period == period) & (dpool.detector == "robust e-detector") & (dpool.band == 0.25)].iloc[0]
            rows.append(_row("F1", "e-detector (band 0.25) minus matched-delay control", period, dn, "next-day robust vol vs baseline",
                             "robust e-detector, band 0.25", "matched-delay control", bt,
                             (a.day, _hits(a, "post_over_baseline")), (c.day, _hits(c, "post_over_baseline")),
                             orig=dict(orig_share_a=dp[f"{dn}_share_in_flagged_direction"], orig_lo=dp[f"{dn}_wilson_lo"],
                                       orig_hi=dp[f"{dn}_wilson_hi"], orig_share_b=dp[f"{dn}_matched_delay_control_mean"])))
        # F2: forecast hits of the e-detector (band 0.25) against random times
        for fc in ("EWMA", "HAR"):
            col = f"y_over_{fc.lower()}"
            r = sel("random times", kind="random")
            for dn in ("up", "down"):
                a = sel("robust e-detector", 0.25, direction=dn)
                rr = r.assign(direction=dn)
                po = pooled[(pooled.period == period) & (pooled.forecast == fc) & (pooled.rule == "robust e-detector") & (pooled.band == 0.25)].iloc[0]
                pr = pooled[(pooled.period == period) & (pooled.forecast == fc) & (pooled.rule == "random times")].iloc[0]
                rows.append(_row("F2", f"e-detector (band 0.25) minus random times, {fc}", period, dn, f"next-day robust vol vs {fc}",
                                 "robust e-detector, band 0.25", "random times", bt, (a.day, _hits(a, col)), (rr.day, _hits(rr, col)),
                                 orig=dict(orig_share_a=po[f"{dn}_hit_rate"], orig_lo=po[f"{dn}_wilson_lo"], orig_hi=po[f"{dn}_wilson_hi"],
                                           orig_share_b=pr[f"{dn}_hit_rate"])))
        # F3: e-detector against the ratio rules, all bands
        for band in (0.0, 0.25, 0.5):
            for fc in ("EWMA", "HAR"):
                col = f"y_over_{fc.lower()}"
                for short in (60, 240):
                    for dn in ("up", "down"):
                        a = sel("robust e-detector", band, direction=dn)
                        b = sel("ratio rule", band, short=short, direction=dn)
                        po = pooled[(pooled.period == period) & (pooled.forecast == fc) & (pooled.rule == f"ratio rule matched to band {band}")
                                    & (pooled.short_window == short)].iloc[0]
                        rows.append(_row("F3", f"e-detector minus ratio rule {short} min, band {band}, {fc}", period, dn,
                                         f"next-day robust vol vs {fc}", f"robust e-detector, band {band}", f"ratio rule {short} min", bt,
                                         (a.day, _hits(a, col)), (b.day, _hits(b, col)),
                                         orig=dict(orig_share_b=po[f"{dn}_hit_rate"], orig_diff=po[f"{dn}_edetector_minus_this"],
                                                   orig_p=po[f"{dn}_p_edetector_vs_this"], orig_p_kind="two-proportion z-test")))
        # F4: robust band 0 against the Gaussian LR monitor (EWMA)
        for dn in ("up", "down"):
            a = sel("robust e-detector", 0.0, direction=dn)
            b = sel("Gaussian LR", 0.0, direction=dn)
            po = pooled[(pooled.period == period) & (pooled.forecast == "EWMA") & (pooled.rule == "Gaussian LR")].iloc[0]
            rows.append(_row("F4", "robust e-detector (band 0) minus Gaussian LR, EWMA", period, dn, "next-day robust vol vs EWMA",
                             "robust e-detector, band 0", "Gaussian LR", bt, (a.day, _hits(a, "y_over_ewma")), (b.day, _hits(b, "y_over_ewma")),
                             orig=dict(orig_diff=po[f"{dn}_edetector_minus_this"], orig_p=po[f"{dn}_p_edetector_vs_this"],
                                       orig_p_kind="two-proportion z-test")))
    return rows


def _first_run(method):
    return ("k=" not in method) and ("w=" not in method)


def families_e13(L_crypto, L_nifty, nifty_sessions):
    al_c = pd.read_csv(OUT / "e13d_real_long_windows_crypto_alarms.csv")
    al_n = pd.read_csv(OUT / "e13d_real_long_windows_nifty_alarms.csv")
    lp = pd.read_csv(OUT / "e13d_real_long_windows_pooled.csv")
    rows = []
    markets = [("five Binance coins, Jul to Dec 2025", P1, al_c, "crypto"), ("five Binance coins, Jan to Jun 2025", "Jan-Jun 2025", al_c, "crypto"),
               ("Nifty 50", "2024 to 2025", al_n, "nifty")]
    for market, period, al, kind in markets:
        a = al[al.market == market].copy()
        if kind == "crypto":
            pstart = pd.Timestamp(PERIOD_DAYS[period][0], tz="UTC")
            a["day"] = (pd.to_datetime(a.time, utc=True).dt.floor("D") - pstart).dt.days
            bt = Boot(_days(period), L_crypto)
        else:
            sess = pd.to_datetime(a.time).map(lambda x: x.date())
            pos = {d: i for i, d in enumerate(nifty_sessions)}
            a["day"] = sess.map(pos)
            bt = Boot(len(nifty_sessions), L_nifty)
        e = a[a.method == "e-detector, band 0.25"]
        comps = [m for m in a.method.unique() if m.endswith("matched to band 0.25")]
        for m in comps:
            fam = ("F6" if _first_run(m) else "F7") if kind == "crypto" else ("F8a" if _first_run(m) else "F8b")
            c = a[a.method == m]
            o = lp[(lp.market == market) & (lp.method == m)].iloc[0]
            for dn in ("up", "down"):
                ea, ca = e[e.direction == dn], c[c.direction == dn]
                rows.append(_row(fam, f"e-detector minus {m.replace(', matched to band 0.25', '')}", period, dn, "next-day robust vol vs EWMA",
                                 "e-detector, band 0.25", m, bt, (ea.day, ea.hit_next_day), (ca.day, ca.hit_next_day),
                                 orig=dict(orig_share_b=getattr(o, f"{dn}_rate_day"), orig_p=getattr(o, f"{dn}_day_p_ztest"),
                                           orig_p_kind="two-proportion z-test",
                                           orig_p_week_clustered=getattr(o, f"{dn}_day_p_week_clustered"),
                                           orig_p_week_clustered_unpaired=getattr(o, f"{dn}_day_p_week_clustered_unpaired"))))
    return rows


def families_e7(L, sessions_by_period):
    a7 = pd.read_csv(OUT / "e07b_alarms.csv")
    s7 = pd.read_csv(OUT / "e07b_monitor_summary.csv")
    rows = []
    for period, sess in sessions_by_period.items():
        pos = {d: i for i, d in enumerate(sess)}
        bt = Boot(len(sess), L)
        a = a7[a7.period == period].copy()
        a["day"] = pd.to_datetime(a.session).map(lambda x: pos[x.date()])
        rob = a[(a.detector == "robust e-detector") & (a.cal_len == 750) & (a.band == 0.25)]
        gau = a[(a.detector == "Gaussian LR") & (a.cal_len == 750)]
        sr = s7[(s7.detector == "robust e-detector") & (s7.cal_len == 750) & (s7.band == 0.25) & (s7.period == period)].iloc[0]
        for dn in ("up", "down"):
            def hit(df):
                df = df[df.direction == dn]
                h = np.where(dn == "up", df.next_vol > df.baseline, df.next_vol < df.baseline).astype(float)
                h[~np.isfinite(df.next_vol.to_numpy())] = np.nan
                return df.day.to_numpy(), h
            ctrl = float(getattr(sr, f"control_{dn}_next_vs_baseline"))
            rows.append(_row("F9a", "Nifty primary configuration minus random-time control (fixed)", period, dn, "next-session robust vol vs baseline",
                             "robust e-detector, cal 750, band 0.25", "random-time control (fixed)", bt, hit(rob), fixed_b=ctrl,
                             orig=dict(orig_share_a=getattr(sr, f"{dn}_next_vs_baseline"), orig_lo=getattr(sr, f"{dn}_next_vs_baseline_lo"),
                                       orig_hi=getattr(sr, f"{dn}_next_vs_baseline_hi"), orig_share_b=ctrl)))
            rows.append(_row("F9b", "Nifty robust (primary) minus Gaussian LR (cal 750)", period, dn, "next-session robust vol vs baseline",
                             "robust e-detector, cal 750, band 0.25", "Gaussian LR, cal 750", bt, hit(rob), hit(gau)))
    return rows


def families_e14(L):
    a = pd.read_csv(OUT / "e14d_real_alarms.csv")
    a = a[a.source.str.endswith("USDT") & (a.band == 0.25)].copy()
    a["day"] = a["index"] // DAY
    bt = Boot(_days(P1), L)
    rows = []
    pl = a[a.method == "plugin"]
    for m in ("awake-A", "awake-B"):
        b = a[a.method == m]
        for dn in ("up", "down"):
            x, y = b[b.direction == dn], pl[pl.direction == dn]
            rows.append(_row("F14", f"{m} minus plug-in, five coins, band 0.25", P1, dn, "next-day robust vol vs baseline", m, "plugin", bt,
                             (x.day, _hits(x, "post_over_baseline")), (y.day, _hits(y, "post_over_baseline"))))
    return rows


def families_e15(L):
    d = pd.read_csv(OUT / "e15d_discoveries.csv.gz")
    t15 = pd.read_csv(OUT / "e15d_rate_matched_tests.csv")
    rows = []
    for (panel, period), g in d.groupby(["panel", "period"], sort=False):
        days = pd.to_datetime(g.day)
        d0, d1 = days.min(), days.max()
        bt = Boot((d1 - d0).days + 1, L)
        g = g.assign(dayi=(days - d0).dt.days)
        for r in t15[(t15.panel == panel) & (t15.period == period)].itertuples():
            gg = g[(g.band == r.band) & (g.direction == r.direction)]
            e, m = gg[gg.rule == "e-BH"], gg[gg.rule == "rate-matched threshold (hindsight)"]
            rows.append(_row("F15", f"e-BH minus rate-matched threshold, {panel}, band {r.band}", period, r.direction,
                             "next-day robust vol vs baseline", "e-BH", "rate-matched threshold (hindsight)", bt,
                             (e.dayi, _hits(e, "post_ratio")), (m.dayi, _hits(m, "post_ratio")),
                             orig=dict(orig_p=r.p_value, orig_p_kind="two-proportion z-test")))
    return rows


def families_e17(L):
    it = pd.read_csv(OUT / "e20_items_e17.csv.gz")
    p17 = pd.read_csv(OUT / "e17_ratematch_pooled.csv")
    bt = Boot(_days(P1), L)
    rows = []
    for target in ("5 min", "15 min"):
        pr = p17[(p17.comparison == target) & p17.monitor.str.startswith("p-value")].iloc[0]
        for col, key, meas in (("post_over_baseline", "confirmed_share", "next-day robust vol vs baseline"),
                               ("y_over_ewma", "ewma_hit_rate", "next-day robust vol vs EWMA"), ("y_over_har", "har_hit_rate", "next-day robust vol vs HAR")):
            for dn in ("up", "down"):
                a = it[(it.comparison == target) & (it.monitor == f"1 min matched to {target}") & (it.direction == dn)]
                b = it[(it.comparison == target) & (it.monitor == target) & (it.direction == dn)]
                rows.append(_row("F16", f"one-minute monitor matched to {target} minus the {target} detector", P1, dn, meas,
                                 f"1 min matched to {target}", target, bt, (a.day, _hits(a, col)), (b.day, _hits(b, col)),
                                 orig=dict(orig_p=getattr(pr, f"{dn}_{key}"), orig_p_kind="two-proportion z-test")))
    return rows


def _nifty_sessions():
    from e07_nifty import DEV_END, load_clean

    df, _, _ = load_clean()
    sess = pd.Index(df["session"].unique())
    dates = [s.date() for s in sess]
    dev = [d for s, d in zip(sess, dates) if s < DEV_END]
    held = [d for s, d in zip(sess, dates) if s >= DEV_END]
    y2425 = [d for d in dates if 2024 <= d.year <= 2025]
    return {"dev": dev, "held-out": held}, y2425


def part_bootstrap():
    t0 = time.time()
    by_period, y2425 = _nifty_sessions()
    all_rows = []
    for lc, ln, tag in ((BLOCK_CRYPTO, BLOCK_NIFTY, "primary"), (2 * BLOCK_CRYPTO, 2 * BLOCK_NIFTY, "sensitivity")):
        rows = families_crypto_e6(lc) + families_e13(lc, ln, y2425) + families_e7(ln, by_period) + families_e14(lc) + families_e15(lc) + families_e17(lc)
        df = pd.DataFrame(rows)
        df["analysis"] = tag
        df["p_holm"] = np.nan
        for fam, g in df.groupby("family"):
            df.loc[g.index, "p_holm"] = holm(g.p_boot.to_numpy())
        if "orig_p" in df:
            df["orig_p_holm"] = np.nan
            for fam, g in df.groupby("family"):
                df.loc[g.index, "orig_p_holm"] = holm(g.orig_p.to_numpy())
        df["equivalence_within_5_points"] = (df.diff_lo >= -EQUIV_MARGIN) & (df.diff_hi <= EQUIV_MARGIN)
        all_rows.append(df)
    out = pd.concat(all_rows, ignore_index=True)
    out.to_csv(OUT / "e20_bootstrap.csv", index=False, float_format="%.6g")
    pd.set_option("display.width", 250)
    prim = out[out.analysis == "primary"]
    print(prim[["family", "test", "period", "direction", "n_a", "share_a", "diff", "diff_lo", "diff_hi", "p_boot", "p_holm", "orig_p"]]
          .round(4).to_string())
    print(f"bootstrap: {time.time() - t0:.0f} s")


# ----------------------------------------------------------------------------
# Part: Holm on published p-values (F5, F10a, F10c, F11, F12, F13)
# ----------------------------------------------------------------------------


def part_holm():
    rows = []

    def fam(name, df, pcol, label_cols, note=""):
        df = df.copy()
        df["family"] = name
        df["p_original"] = df[pcol].astype(float)
        df["p_holm"] = holm(df["p_original"].to_numpy())
        df["test"] = df[label_cols].astype(str).agg(", ".join, axis=1)
        df["note"] = note
        rows.append(df[["family", "test", "p_original", "p_holm", "note"]])

    cv = pd.read_csv(OUT / "e06d_cross_venue.csv")
    fam("F5", cv, "p_same", ["detector", "band", "alarms_of", "matched_against"], "circular-shift p-value, same-direction alarm within 60 minutes")
    er = pd.read_csv(OUT / "e07d_event_alarm_rates.csv")
    fam("F10a", er[(er.outcome == "any alarm") & (er.sessions_used == "all")], "fisher_p", ["event_type", "outcome", "sessions_used"], "Fisher exact test")
    fam("F10c", er[(er.outcome == "any alarm") & (er.sessions_used == "monitoring at open")], "fisher_p", ["event_type", "outcome", "sessions_used"],
        "Fisher exact test, sessions monitoring at the open")
    vx = pd.read_csv(OUT / "e07c_vix_alarm_check.csv")
    vx = vx[(vx.detector.str.startswith("robust")) & (vx.cal_len == 750) & (vx.band == 0.25) & vx.period.isin(["dev", "held-out"])]
    fam("F11", vx, "fisher_p_vs_other_sessions", ["period", "direction", "horizon"], "Fisher exact test, primary configuration")
    t = pd.read_csv(OUT / "e12_tests.csv")
    rob = t[(t.detector == "robust e-detector") & (t.outcome == "up_post")]
    s_all = rob[(rob.design == "streaming") & (rob["sample"] == "all")]
    fam("F12a", s_all[s_all.coin == "BTC"], "p_stratified_permutation", ["coin", "design", "event_type"], "pre-registered primary family (Holm as published)")
    fam("F12b", s_all[s_all.coin != "BTC"], "p_stratified_permutation", ["coin", "design", "event_type"], "replications")
    fam("F12c", rob[(rob.design == "fresh calibration") & (rob["sample"] == "all")], "p_stratified_permutation", ["coin", "design", "event_type"],
        "secondary design added after the first results")
    fam("F12d", rob[(rob.design == "streaming") & (rob["sample"] == "monitoring at start")], "p_stratified_permutation",
        ["coin", "design", "event_type", "sample"], "post hoc restriction to releases with the monitor awake")
    g = pd.read_csv(OUT / "e14g_e12_tests.csv")
    fam("F13", g[g.method != "plugin"], "stratified_p", ["coin", "method", "event_type"], "E14 Q14g, restarted and rank monitors")
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT / "e20_holm_existing.csv", index=False, float_format="%.6g")
    pd.set_option("display.width", 250)
    print(out.round(5).to_string())


# ----------------------------------------------------------------------------
# Part: circular-shift test of the Nifty event calendar (F10b)
# ----------------------------------------------------------------------------


def part_events():
    from e07_nifty import load_clean
    from e07d_events import TYPES, tag_events

    df, _, _ = load_clean()
    sess = pd.Index(df["session"].unique())
    s = pd.DataFrame({"session": sess})
    s["date"] = s["session"].dt.date
    a7 = pd.read_csv(OUT / "e07b_alarms.csv")
    a = a7[(a7.detector == "robust e-detector") & (a7.cal_len == 750) & (a7.band == 0.25)]
    alarm_dates = set(pd.to_datetime(a.session).map(lambda x: x.date()))
    s["n_alarms"] = s["date"].isin(alarm_dates).astype(int)
    s, _, _ = tag_events(s)
    x = s["n_alarms"].to_numpy() > 0
    N = x.size
    ref = pd.read_csv(OUT / "e07d_event_alarm_rates.csv")
    ref = ref[(ref.outcome == "any alarm") & (ref.sessions_used == "all")]
    rows = []
    for t in TYPES:
        ev = np.flatnonzero(s[t].to_numpy())
        T = np.array([x[(ev + k) % N].sum() for k in range(N)])
        dev = np.abs(T - T.mean())
        p = float(np.mean(dev >= dev[0] - 1e-9))
        r = ref[ref.event_type == t].iloc[0]
        o = s["ordinary"].to_numpy()
        rows.append(dict(event_type=t, event_sessions=ev.size, with_alarm=int(T[0]), shift_mean=float(T.mean()),
                         shift_q025=float(np.percentile(T, 2.5)), shift_q975=float(np.percentile(T, 97.5)),
                         p_shift_two_sided=p, ordinary_with=int(x[o].sum()), ordinary_sessions=int(o.sum()),
                         check_against_e07d=bool(T[0] == r.event_sessions_with and ev.size == r.event_sessions
                                                 and int(x[o].sum()) == r.ordinary_sessions_with),
                         fisher_p_original=r.fisher_p))
    out = pd.DataFrame(rows)
    out["p_shift_holm"] = holm(out.p_shift_two_sided)
    out["fisher_p_holm"] = holm(out.fisher_p_original)
    out.to_csv(OUT / "e20_nifty_event_shift.csv", index=False, float_format="%.6g")
    print(out.round(5).to_string())


# ----------------------------------------------------------------------------
# Part 3: forecast regressions and QLIKE
# ----------------------------------------------------------------------------

HORIZONS = (1440, 360)
MEASURES = ("RV", "BV")
MODELS = ("EWMA", "HAR", "HAR-I")
MIN_TRAIN = 720
EWMA_HALFLIFE = 1440


def _coin_grid(coin):
    """Hourly grid of one coin with outcomes, regressors and past-only forecasts (pre-registration 3.3 to 3.5)."""
    from avmonitor.data import load_binance

    from e06b_five_coins_monitor import COINS as SYM

    df = load_binance(SYM[coin], "2024-01-01", "2026-01-01", download=False)
    r = df["logret"].to_numpy(float)
    n = r.size
    miss = ~np.isfinite(r)
    r0 = np.where(miss, 0.0, r)
    C = np.concatenate([[0.0], np.cumsum(r0 ** 2)])  # C[k] = sum of r^2 over indices < k
    M = np.concatenate([[0], np.cumsum(miss)])
    ab = np.zeros(n)
    ab[1:] = np.abs(r0[1:]) * np.abs(r0[:-1])
    abm = np.zeros(n, bool)
    abm[1:] = miss[1:] | miss[:-1]
    A = np.concatenate([[0.0], np.cumsum(ab)])
    AM = np.concatenate([[0], np.cumsum(abm)])
    g = np.arange(0, n, 60)
    times = df.index[g]
    out = pd.DataFrame({"coin": coin, "t": g, "time": times})
    out["dayn"] = g // DAY

    def win(Cs, Ms, a, b):  # sum over indices a..b inclusive, NaN if any missing or out of range
        a = np.asarray(a)
        b = np.asarray(b)
        ok = (a >= 0) & (b < n)
        aa, bb = np.clip(a, 0, n - 1), np.clip(b, 0, n - 1)
        s = Cs[bb + 1] - Cs[aa]
        bad = (Ms[bb + 1] - Ms[aa]) > 0
        return np.where(ok & ~bad, s, np.nan)

    for H in HORIZONS:
        out[f"RV{H}"] = win(C, M, g + 1, g + H)
        out[f"BV{H}"] = (np.pi / 2) * (H / (H - 1)) * win(A, AM, g + 2, g + H)
    rvd = lambda tt: win(C, M, tt - DAY + 1, tt)
    out["RVd"] = rvd(g)
    out["RVw"] = np.mean([rvd(g - DAY * j) for j in range(5)], axis=0)
    out["RVm"] = np.mean([rvd(g - DAY * j) for j in range(22)], axis=0)
    out["RVh"] = np.maximum(win(C, M, g - 59, g), 0.01 * out["RVd"] / 24)
    ew = pd.Series(np.where(miss, np.nan, r ** 2)).ewm(halflife=EWMA_HALFLIFE, ignore_na=True).mean().to_numpy()
    out["E"] = ew[g]
    out["hour"] = times.hour
    out["dow"] = times.dayofweek
    return out


def _design(df, model):
    cols = [np.ones(len(df))]
    if model == "EWMA":
        cols.append(np.log(df["E"].to_numpy()))
    else:
        cols += [np.log(df[c].to_numpy()) for c in ("RVd", "RVw", "RVm")]
        if model == "HAR-I":
            cols.append(np.log(df["RVh"].to_numpy()))
    for h in range(1, 24):
        cols.append((df["hour"].to_numpy() == h).astype(float))
    for d in range(1, 7):
        cols.append((df["dow"].to_numpy() == d).astype(float))
    return np.column_stack(cols)


def _fit_forecasts(df):
    """Daily refitted expanding-window OLS of log outcome on each model's regressors (past data only)."""
    ndays = int(df["dayn"].max()) + 2
    for H in HORIZONS:
        for meas in MEASURES:
            y = np.log(df[f"{meas}{H}"].to_numpy())
            known = (df["t"].to_numpy() + H) // DAY + 1  # first day on which the outcome is fully observed
            for model in MODELS:
                X = _design(df, model)
                ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
                k = X.shape[1]
                XtX = np.zeros((ndays + 1, k, k))
                Xty = np.zeros((ndays + 1, k))
                yy = np.zeros(ndays + 1)
                cnt = np.zeros(ndays + 1)
                kd = np.minimum(known[ok], ndays)
                np.add.at(XtX, kd, X[ok][:, :, None] * X[ok][:, None, :])
                np.add.at(Xty, kd, X[ok] * y[ok][:, None])
                np.add.at(yy, kd, y[ok] ** 2)
                np.add.at(cnt, kd, 1)
                XtX, Xty, yy, cnt = (np.cumsum(a, axis=0) for a in (XtX, Xty, yy, cnt))
                beta = np.full((ndays + 1, k), np.nan)
                s2 = np.full(ndays + 1, np.nan)
                for d in range(ndays + 1):
                    if cnt[d] >= MIN_TRAIN:
                        try:
                            b = np.linalg.solve(XtX[d], Xty[d])
                        except np.linalg.LinAlgError:
                            continue
                        beta[d] = b
                        s2[d] = (yy[d] - 2 * b @ Xty[d] + b @ XtX[d] @ b) / (cnt[d] - k)
                dn = df["dayn"].to_numpy()
                logF = np.sum(X * beta[dn], axis=1)
                df[f"logF_{model}_{meas}{H}"] = logF
                df[f"F_{model}_{meas}{H}"] = np.exp(logF + s2[dn] / 2)
    return df


def _alarm_indicators(df, coin, items, tally):
    """U and D indicators of each rule: an alarm in (t - 60, t] for grid time t (indices from 2024-01-01 00:00).

    ``tally`` collects, per indicator and period, the number of alarms of the five coins, the number that landed on a
    grid point, and the number whose grid point lies past the end of the data (an alarm in the last hour of 2025),
    which therefore enter no indicator; two alarms of one rule, coin and direction in the same hour would count once.
    """
    rules = {"e-detector": ("robust e-detector", 0), "ratio 60": ("ratio rule", 60), "ratio 240": ("ratio rule", 240)}
    tpos = {t: i for i, t in enumerate(df["t"].to_numpy())}
    for band in (0.25, 0.5):
        for key, (rule, short) in rules.items():
            a = items[(items.coin == coin) & (items.kind == "alarm") & (items.rule == rule) & (items["short"] == short)
                      & np.isclose(items.band, band)]
            for dn, lab in (("up", "U"), ("down", "D")):
                col = f"{lab}_{key}_{band}"
                v = np.zeros(len(df))
                ad = a[a.direction == dn]
                for per, ai in zip(ad["period"].to_numpy(), ad["index"].to_numpy(int)):
                    tg = -(-ai // 60) * 60
                    rec = tally.setdefault((col, per), dict(alarms=0, on_grid=0, past_end=0))
                    rec["alarms"] += 1
                    if tg in tpos:
                        v[tpos[tg]] += 1.0
                        rec["on_grid"] += 1
                    else:
                        rec["past_end"] += 1
                df[col] = np.minimum(v, 1.0)
                tally.setdefault(("same hour", col), dict(alarms=0, on_grid=0, past_end=0))["alarms"] += int(np.sum(v > 1))
    return df


def part_forecast():
    t0 = time.time()
    items = _items_e6()
    grids, tally = [], {}
    for coin in COINS:
        g = _coin_grid(coin)
        g = _fit_forecasts(g)
        g = _alarm_indicators(g, coin, items, tally)
        grids.append(g)
        print(f"forecast grid {coin}: {time.time() - t0:.0f} s", flush=True)
    G = pd.concat(grids, ignore_index=True)
    # periods: grid time inside the period and the outcome window ending inside it
    G["period"] = ""
    for p in (P1, P2):
        a, b = (pd.Timestamp(x, tz="UTC") for x in PERIOD_DAYS[p])
        G.loc[(G.time >= a) & (G.time < b), "period"] = p
    end = {p: pd.Timestamp(PERIOD_DAYS[p][1], tz="UTC") for p in (P1, P2)}
    G["gday"] = (G.time.dt.floor("D") - pd.Timestamp("2024-01-01", tz="UTC")).dt.days
    # pre-registration 3.4: every alarm of the e-detector and ratio rules must land on a grid time of its own period
    chk = []
    for c in [c for c in G.columns if c.startswith(("U_", "D_"))]:
        t1, t2 = (tally.get((c, p), dict(alarms=0, on_grid=0, past_end=0)) for p in (P1, P2))
        chk.append(dict(indicator=c, alarms_P1=t1["alarms"], alarms_P2=t2["alarms"], P1=int(G.loc[G.period == P1, c].sum()),
                        P2=int(G.loc[G.period == P2, c].sum()), outside=int(G.loc[G.period == "", c].sum()),
                        past_end_of_data=t1["past_end"] + t2["past_end"], same_hour_pairs=tally[("same hour", c)]["alarms"]))
    choice, regs, ql = [], [], []

    def sample(p, H, meas):
        m = (G.period == p) & (G.time + pd.Timedelta(minutes=H) < end[p]) & np.isfinite(np.log(G[f"{meas}{H}"]))
        for model in MODELS:
            m &= np.isfinite(G[f"logF_{model}_{meas}{H}"]) & np.isfinite(G[f"F_{model}_{meas}{H}"])
        return G[m]

    # strongest forecast: smallest mean QLIKE of next-day RV on the development period
    for p in (P1, P2):
        for H in HORIZONS:
            for meas in MEASURES:
                s = sample(p, H, meas)
                for model in MODELS:
                    q = qlike(s[f"{meas}{H}"], s[f"F_{model}_{meas}{H}"])
                    e = np.log(s[f"{meas}{H}"]) - s[f"logF_{model}_{meas}{H}"]
                    choice.append(dict(period=p, horizon=H, measure=meas, model=model, n=len(s), days=s.gday.nunique(),
                                       first=str(s.time.min()), last=str(s.time.max()),
                                       mean_qlike=float(q.mean()), mse_log=float(np.mean(e ** 2)), mean_log_error=float(e.mean())))
    ch = pd.DataFrame(choice)
    dev = ch[(ch.period == P1) & (ch.horizon == 1440) & (ch.measure == "RV")]
    best = dev.sort_values("mean_qlike").model.iloc[0]
    ch["strongest"] = ch.model == best
    rules = ("e-detector", "ratio 60", "ratio 240")
    for p in (P1, P2):
        for H in HORIZONS:
            for meas in MEASURES:
                s = sample(p, H, meas)
                y = np.log(s[f"{meas}{H}"].to_numpy())
                coin_d = np.column_stack([(s.coin == c).to_numpy(float) for c in COINS])
                for model in MODELS:
                    lf = s[f"logF_{model}_{meas}{H}"].to_numpy()
                    for band in (0.25, 0.5):
                        specs = [(r,) for r in rules] + [("e-detector", "ratio 60"), ("e-detector", "ratio 240")]
                        for spec in specs:
                            cols, names = [coin_d, lf[:, None]], ["b_logF"]
                            for r in spec:
                                cols += [s[f"U_{r}_{band}"].to_numpy()[:, None], s[f"D_{r}_{band}"].to_numpy()[:, None]]
                                names += [f"c_up[{r}]", f"c_dn[{r}]"]
                            X = np.hstack(cols)
                            for bw in (DK_BANDWIDTH, 0):
                                beta, se, nd = ols_dk(y, X, s.gday.to_numpy(), bw)
                                for j, nm in enumerate(names):
                                    b, sj = beta[len(COINS) + j], se[len(COINS) + j]
                                    rule = nm[nm.find("[") + 1:-1] if "[" in nm else ""
                                    regs.append(dict(period=p, horizon=H, measure=meas, model=model, band=band, spec=" + ".join(spec),
                                                     coef=nm.split("[")[0], rule=rule, estimate=b, se=sj, z=b / sj, p=norm_p(b / sj),
                                                     se_kind="Driscoll-Kraay, 5 days" if bw else "clustered by day",
                                                     n=len(s), days=nd, alarms=int(s[f"{'U' if 'up' in nm else 'D'}_{rule}_{band}"].sum()) if rule else np.nan,
                                                     strongest=model == best))
    R = pd.DataFrame(regs)
    R["family"] = ""
    prim = (R.period == P2) & (R.horizon == 1440) & (R.measure == "RV") & (R.strongest) & (R.band == 0.25) & (R.spec.isin(rules)) \
        & (R.coef != "b_logF") & (R.se_kind.str.startswith("Driscoll"))
    R.loc[prim, "family"] = "F17"
    R["p_holm"] = np.nan
    R.loc[prim, "p_holm"] = holm(R.loc[prim, "p"].to_numpy())
    blk = R[(R.coef != "b_logF") & R.spec.isin(rules)]
    for _, g in blk.groupby(["period", "horizon", "measure", "model", "band", "se_kind"]):
        if g.family.eq("F17").all():
            continue
        R.loc[g.index, "p_holm"] = holm(g.p.to_numpy())
    # QLIKE: forecast against the forecast scaled at alarms (kappas from one period, evaluated on the other)
    for fit_p, ev_p in ((P1, P2), (P2, P1)):
        for H in HORIZONS:
            meas = "RV"
            sf, se_ = sample(fit_p, H, meas), sample(ev_p, H, meas)
            for model in MODELS:
                for band in (0.25, 0.5):
                    for r in rules:
                        kap = {}
                        for lab in ("U", "D"):
                            m = sf[f"{lab}_{r}_{band}"] == 1
                            kap[lab] = qlike_kappa(sf.loc[m, f"{meas}{H}"], sf.loc[m, f"F_{model}_{meas}{H}"])
                        y = se_[f"{meas}{H}"].to_numpy()
                        F = se_[f"F_{model}_{meas}{H}"].to_numpy()
                        U, Dn = se_[f"U_{r}_{band}"].to_numpy(), se_[f"D_{r}_{band}"].to_numpy()
                        F2 = F * kap["U"] ** U * kap["D"] ** Dn
                        dl = qlike(y, F2) - qlike(y, F)
                        beta, se, nd = ols_dk(dl, np.ones((dl.size, 1)), se_.gday.to_numpy(), DK_BANDWIDTH)
                        al = (U + Dn) > 0
                        ql.append(dict(fit_period=fit_p, eval_period=ev_p, horizon=H, measure=meas, model=model, band=band, rule=r,
                                       kappa_up=kap["U"], kappa_down=kap["D"], n=dl.size, alarm_points=int(al.sum()),
                                       mean_qlike_forecast=float(qlike(y, F).mean()), mean_diff=float(beta[0]), se=float(se[0]),
                                       z=float(beta[0] / se[0]), p=norm_p(beta[0] / se[0]),
                                       mean_diff_at_alarms=float(dl[al].mean()) if al.any() else np.nan,
                                       relative_change=float(beta[0] / qlike(y, F).mean()),
                                       up_mean_y_over_F_eval=float(np.mean(y[U == 1] / F[U == 1])) if (U == 1).any() else np.nan,
                                       down_mean_y_over_F_eval=float(np.mean(y[Dn == 1] / F[Dn == 1])) if (Dn == 1).any() else np.nan,
                                       strongest=model == best))
    Q = pd.DataFrame(ql)
    Q["family"] = ""
    primq = (Q.fit_period == P1) & (Q.horizon == 1440) & (Q.strongest) & (Q.band == 0.25)
    Q.loc[primq, "family"] = "F18"
    Q["p_holm"] = np.nan
    Q.loc[primq, "p_holm"] = holm(Q.loc[primq, "p"].to_numpy())
    for _, g in Q[~primq].groupby(["fit_period", "horizon", "model", "band"]):
        Q.loc[g.index, "p_holm"] = holm(g.p.to_numpy())
    ch.to_csv(OUT / "e20_forecast_choice.csv", index=False, float_format="%.6g")
    R.to_csv(OUT / "e20_forecast_regressions.csv", index=False, float_format="%.6g")
    Q.to_csv(OUT / "e20_forecast_qlike.csv", index=False, float_format="%.6g")
    pd.DataFrame(chk).to_csv(OUT / "e20_forecast_alarm_counts.csv", index=False)
    pd.set_option("display.width", 250)
    print(ch.round(4).to_string())
    print(pd.DataFrame(chk).to_string())
    print(R[R.family == "F17"].round(4).to_string())
    print(Q[Q.family == "F18"].round(5).to_string())
    print(f"forecast: {time.time() - t0:.0f} s")


# ----------------------------------------------------------------------------
# Part: classification table
# ----------------------------------------------------------------------------

# Every real-data claim of paper/short_paper.md at commit 8f71a2e (line numbers of that version). Classes (pre-registration
# Part 1): P pre-specified primary, S pre-specified secondary, X exploratory or post hoc; "verifiable" when the describing
# commit is strictly earlier than the first commit with the results, "stated" when the description is committed together
# with the results.
CLASSIFICATION = [
    # (claim_id, claim, line, source, data set, class, verifiable, evidence)
    ("C1", "On five cryptocurrencies the alarms carry information beyond EWMA and HAR volatility forecasts", "9, 404 to 419",
     "E6 Q6e", "five coins (E6)", "S", "stated",
     "Q6e is in the docstring of e06b_five_coins_monitor.py, first committed together with its results (d0caeb1); the monitor "
     "settings were frozen in E5 (bf67e2b), which precedes E6 (verifiable). No test is named primary. The HAR forecast's training "
     "history was changed after the first full run (docs/iterations/E6.md, iteration 5), so the HAR rows are post hoc in their "
     "forecast definition."),
    ("C2", "On the Nifty 50 index and around FOMC statements the alarms follow scheduled events far more often than control times",
     "9", "E7 Q7d (Budget); E12 (FOMC, BTC)", "Nifty events; FOMC and CPI", "S (Budget); P (FOMC on BTC)", "stated",
     "Budget: one of five event types in the docstring of e07d_events.py, committed with the results (c20ea73); no primary event "
     "type. FOMC on BTC: member of the primary family stated in the docstring of e12_crypto_events.py, committed with the results "
     "(47a5982)."),
    ("C3", "A ratio rule tuned with hindsight is about as precise", "9, 419", "E6 Q6e", "five coins (E6)", "S", "stated",
     "The rate-matched ratio rule is in the E6b docstring (d0caeb1, with results); windows 60 and 240 minutes 'fixed before "
     "running' (docs/iterations/E6.md, iteration 6). The claim is an equivalence claim; no margin was set in advance."),
    ("C4", "Published change-point methods given day-long windows come close", "9, 441", "E13 Q13d, longer windows",
     "five coins (E13)", "X", "verifiable",
     "The longer windows were added in review round 3 (fca74a5) after the first real results (53f9eba); docstring: 'Added after "
     "the first runs (review round 3)'."),
    ("C5", "Restarted monitors were confirmed by the next day's volatility about as often (five coins, band 0.25)", "379, 457, 788",
     "E14 Q14d", "five coins (E14)", "S", "verifiable",
     "Q14d is in the E14 docstring at d9871f4 (2 Nov 21:20), before the results (e955821, 21:52). An equivalence claim with no "
     "margin set in advance."),
    ("C6", "Restarted monitors: the BTC CPI cell gave p = 0.003 and 0.0003 instead of 0.11", "379, 449, 792", "E14 Q14g",
     "FOMC and CPI", "X", "verifiable",
     "Q14g was described before its own run (d9871f4, 21:20) but designed after E12's results, in which the BTC CPI test failed, "
     "were committed (47a5982, 21:01); it is outside E12's pre-registered primary family and cannot rescue it."),
    ("C7", "At matched alarm counts the 5- and 15-minute detectors alone were no more precise than the one-minute monitor",
     "381, 834", "E17 Q17g", "five coins (E17)", "X", "verifiable",
     "part_ratematch is marked 'POST HOC (review round 4, M5)' (ef65698), after the real results (5da9ac2)."),
    ("C8", "75% to 84% of increase alarms and 83% to 93% of decrease alarms were followed by a day on the flagged side of the "
     "baseline, against about 50% at random times with matched delays", "402", "E6 Q6d", "five coins (E6)", "S", "stated",
     "Q6d is in the E6b docstring, committed with its results (d0caeb1); no primary endpoint."),
    ("C9", "Out of sample these shares fell by only 1 to 4 points", "402", "E6 Q6d", "five coins (E6)", "S", "stated",
     "Descriptive comparison of the two periods named in the docstring (d0caeb1)."),
    ("C10", "Bitstamp had a same-direction alarm within 60 minutes of 72% of Binance BTC alarms (73% the other way), against "
     "about 3% for shifted times", "402", "E6 Q6d cross-venue", "five coins (E6)", "S", "stated",
     "Named in the E6b docstring (d0caeb1), committed with the results."),
    ("C11", "After an increase alarm next-day volatility exceeded both forecasts 76% to 81% of the time, against 48% to 51% at "
     "random times; after a decrease alarm it fell short 71% to 81% of the time", "404 to 419", "E6 Q6e", "five coins (E6)",
     "S (EWMA); X (HAR)", "stated", "As C1."),
    ("C12", "The rate-matched ratio rule was within 0 to 5 points at band 0.25, and none of these differences was significant",
     "419", "E6 Q6e", "five coins (E6)", "S", "stated", "As C3; 48 cells (3 bands x 2 periods x 2 forecasts x 2 windows x 2 "
     "directions) were compared, with no multiplicity plan."),
    ("C13", "Only at band 0.5, out of sample, EWMA: the e-detector better by about 4 points for increases (p = 0.03 to 0.06) and "
     "10 points for decreases against the 60-minute rule (p < 0.001)", "419", "E6 Q6e", "five coins (E6)", "X", "verifiable",
     "One subgroup of the 48 cells, highlighted after the results were seen (docs/iterations/E6.md, iteration 6; referee item 17)."),
    ("C14", "Increase alarms of the Gaussian likelihood-ratio monitor were 7 to 14 points less informative", "419", "E6 Q6e",
     "five coins (E6)", "S", "stated", "The Gaussian LR benchmark is in E6's CONFIGS and docstring (d0caeb1)."),
    ("C15", "First-run windows: e-detector alarms followed by the announced volatility more often than every competitor's "
     "(p <= 0.02; week-clustered p <= 0.010), Jul to Dec 2025", "441", "E13 Q13d", "five coins (E13)", "S", "verifiable",
     "Q13d is in the E13 docstring at 6cfe39c (21:24), before the real results (53f9eba, 22:14); no primary endpoint. The "
     "week-clustered tests were added in review round 3 (fca74a5): X."),
    ("C16", "Held-out half-year Jan to Jun 2025: the same comparisons, except truncated BOCPD down alarms (67.0% against 69.7%, "
     "p = 0.47)", "441", "E13 Q13d held-out", "five coins (E13)", "X", "verifiable",
     "Added after the first real results (iteration 7, 9c782d5) as a replication with every setting frozen."),
    ("C17", "With day-long windows the competitors' down alarms were as informative as the e-detector's (p >= 0.28)", "441",
     "E13 longer windows", "five coins (E13)", "X", "verifiable", "As C4."),
    ("C18", "No competitor reached the e-detector's up alarms; the best longer-window competitor, picked after seeing the "
     "results, was about 6 points lower (p = 0.065 / 0.047; 0.054 / 0.028; unpaired 0.21 / 0.19)", "441", "E13 longer windows",
     "five coins (E13)", "X", "verifiable", "Comparator chosen after the results; three tests without a primary one (referee item 18)."),
    ("C19", "Scored on the next four hours, the short-window competitors did as well as the e-detector or better", "441",
     "E13 four-hour outcome", "five coins (E13)", "X", "verifiable", "Outcome added in review round 3 (fca74a5)."),
    ("C20", "Nifty 50 (2024 and 2025): all five methods scored 69% to 83% with the first-run windows and no difference was "
     "significant; a session-long BJV window did worse", "441", "E13 Q13d Nifty", "Nifty 50 (E13)",
     "S (first-run); X (session-long windows)", "verifiable", "First run as C15; session-long windows as C4."),
    ("C21", "Nifty held-out: 78% of increase alarms followed by a more volatile session than the baseline (47% at random "
     "times), 76% of decrease alarms by a calmer one (53%)", "445, 842", "E7 Q7b", "Nifty 50 monitor (E7)", "S", "stated",
     "The 'pre-registered primary configuration' (calibration 750, band 0.25) is in the docstring of e07_nifty.py, committed "
     "with the results (38856a3); a configuration, not an endpoint, was named."),
    ("C22", "Against 61% for the increase alarms of the Gaussian likelihood-ratio monitor", "445, 842", "E7 Q7b",
     "Nifty 50 monitor (E7)", "S", "stated", "Benchmark named in the E7 docstring (38856a3)."),
    ("C23", "India VIX agreed with the alarms on the same day but did not move in their direction on the next day more often",
     "445, 842", "E7 Q7c", "Nifty 50 and India VIX (E7)", "S", "verifiable",
     "Q7c is in the docstring of e07c_india_vix.py at 0627d14 (20:29), before the VIX results (83ce131, 21:12), but after the "
     "alarms were known (38856a3, 20:16). No primary cell among the eight."),
    ("C24", "9 of 11 Union Budget sessions had an alarm, against 25% of ordinary sessions (Fisher p = 0.0001)", "445, 842",
     "E7 Q7d", "Nifty events (E7)", "S", "stated", "As C2 (Budget)."),
    ("C25", "All 7 budget sessions at which the monitor was watching at the open had an alarm", "842", "E7 Q7d",
     "Nifty events (E7)", "X", "stated", "The 'monitoring at open' filter replaced a flawed first filter after the results were "
     "seen (docs/iterations/E7.md, iteration 5)."),
    ("C26", "Reserve Bank of India policy days did not differ (32% against 25%, p = 0.36)", "445, 842", "E7 Q7d",
     "Nifty events (E7)", "S", "stated", "As C24; a 'no difference' statement."),
    ("C27", "Monthly expiry days had more alarms (38%, p = 0.002), which may be an artefact of the profile", "842", "E7 Q7d",
     "Nifty events (E7)", "S (test); X (interpretation)", "stated", "As C24."),
    ("C28", "In March 2020 only 6 of 21 sessions had an alarm", "445, 842", "E7 Q7d", "Nifty events (E7)", "S", "stated",
     "Descriptive; the COVID window is one of the five event types (c20ea73)."),
    ("C29", "Primary test (BTC, streaming, Holm): FOMC 9 of 16 against 14% of controls (p = 0.0002); CPI not significant "
     "(p = 0.11)", "449, 846 to 848", "E12 Q12a", "FOMC and CPI (E12)", "P", "stated",
     "'the PRIMARY family, fixed before the first run' in the docstring of e12_crypto_events.py at its first commit, which also "
     "holds the results (47a5982). The event calendar was committed earlier (08fb1bf) but holds no analysis plan."),
    ("C30", "The CPI test replicated on ETH (p = 0.003) and SOL (p = 0.00004)", "449, 848", "E12 Q12a", "FOMC and CPI (E12)",
     "S", "stated", "'ETH and SOL are replications' in the same docstring (47a5982)."),
    ("C31", "At the FOMC statements at which the streaming monitor was monitoring, it alarmed every time", "449, 848", "E12 Q12a",
     "FOMC and CPI (E12)", "X", "stated", "Restriction to releases with the monitor awake, introduced after the first run "
     "(docs/iterations/EV1.md, iteration 1, remaining issues)."),
    ("C32", "Fresh calibration: 14 to 15 of 16 FOMC statements against 29% to 32% of controls; CPI 59% to 73% against 23% to 28%",
     "449, 848", "E12 fresh calibration", "FOMC and CPI (E12)", "X", "stated",
     "Added in iteration 2 after the streaming results (docstring: 'added in iteration 2 after the streaming run'); marked secondary."),
    ("C33", "Gaussian LR (fresh): alarms at almost every release but also in 62% to 76% of control windows; robust monitor at band 0: "
     "39% to 51%", "449, 848", "E12 Q12d", "FOMC and CPI (E12)", "X", "stated",
     "Fresh design (iteration 2) and band-0 robust monitor (review round 2a, 990c289) added after the first results."),
    ("C34", "Panel: discoveries pointed the right way (up 78% to 83% against 46% to 50% at matched random times; down 85% to 97% "
     "against 49% to 55%)", "453, 852", "E15 Q15d", "five coins and Nifty panel (E15)", "S", "verifiable",
     "Q15d is in the docstring of e15_panel.py at 2eaa2e3 (21:26), before the real results (d8cf4ae, 21:38)."),
    ("C35", "A per-asset threshold chosen with hindsight had its discoveries confirmed at the same rate (no significant "
     "difference in any of 12 cells)", "852", "E15 Q15d", "five coins and Nifty panel (E15)", "S", "verifiable",
     "The rate-matched threshold is in the docstring at 2eaa2e3, before the results. An equivalence claim with no margin set in advance."),
    ("C36", "5- and 15-minute detectors alone were confirmed more often (up 90% and 95%, down 95% and 96%) than the averages "
     "and the one-minute monitor", "834", "E17 Q17f", "five coins (E17)", "S", "verifiable",
     "Real settings frozen in the script (34c6155) before the real run (5da9ac2); a comparison at unequal alarm rates."),
]

PRIMARY_ENDPOINTS = [
    ("five coins (E6)", "none was pre-specified", "Q6d and Q6e are research questions; no test is named primary."),
    ("five coins (E13)", "none was pre-specified", "Q13d (verifiable, 6cfe39c) names the comparison, not a primary test."),
    ("five coins (E14, E15, E17)", "none was pre-specified", "Q14d, Q15d and Q17f name analyses, not a primary test."),
    ("Nifty 50 monitor (E7)", "none was pre-specified", "A primary configuration (calibration 750, band 0.25) was stated, "
     "in the same commit as the results (38856a3), but no endpoint or test."),
    ("Nifty events (E7)", "none was pre-specified", "Five event types, no primary type or multiplicity plan (c20ea73)."),
    ("Nifty 50 and India VIX (E7)", "none was pre-specified", "Eight cells, no primary cell (0627d14)."),
    ("Nifty 50 (E13)", "none was pre-specified", "As E13 crypto."),
    ("FOMC and CPI (E12)", "BTC, streaming robust monitor (band 0.25), up alarm within 6 hours, FOMC and CPI, exact stratified "
     "permutation test, Holm at 5%", "Stated as fixed before the first run in the docstring committed with the results (47a5982); "
     "not verifiable from the commit order."),
]

# Adjusted status of each claim (pre-registration Part 4), with the families whose rows support it. Every number
# below is a row of results/e20_bootstrap.csv, e20_holm_existing.csv, e20_nifty_event_shift.csv,
# e20_forecast_regressions.csv or e20_forecast_qlike.csv.
STATUS = {
    "C1": ("F2, F17, F18", "partly holds: hit rates against EWMA and HAR forecasts hold (F2); against the strongest forecast "
           "(HAR with the last hour) only decrease alarms add information (F17: e-detector down Holm p = 0.046, up 0.32); "
           "scaling that forecast at alarms did not lower held-out QLIKE (F18)"),
    "C2": ("F10a, F10b, F12a", "holds for Union Budget sessions (F10a, F10b) and FOMC on BTC (F12a); does not hold for CPI on "
           "BTC, the other member of the primary family (F12a, p = 0.11)"),
    "C3": ("F3", "no difference significant after Holm at band 0.25; equivalence within 5 points shown in 7 of 8 held-out "
           "cells and in none of the 8 Jul to Dec 2025 cells"),
    "C4": ("F7", "exploratory; down alarms: no significant difference, equivalence not shown; up alarms: e-detector higher "
           "after Holm in 5 of 8 day-long cells"),
    "C5": ("F14", "does not hold for the rank version's increase alarms (6.7 points fewer, Holm p = 0.006); band version: "
           "no significant difference, equivalence not shown"),
    "C6": ("F13", "exploratory; Holm p = 0.011 and 0.0027 within F13; cannot replace E12's primary CPI test"),
    "C7": ("F16", "exploratory; no cell favours the coarse detector, 2 of 12 favour the one-minute monitor after Holm; "
           "equivalence shown in 1 of 12 cells"),
    "C8": ("F1", "holds (all four Holm p = 0.0008)"),
    "C9": ("F1", "descriptive, not tested; pooled shares fell by 1.6 (up) and 3.8 (down) points"),
    "C10": ("F5", "holds (all eight Holm p = 0.008)"),
    "C11": ("F2", "holds (all eight Holm p = 0.0016); HAR rows exploratory"),
    "C12": ("F3", "holds as 'not significant' (all 16 band 0.25 cells Holm p = 1); equivalence shown out of sample only"),
    "C13": ("F3", "exploratory; decreases against the 60-minute rule hold after Holm (p = 0.0096); increases hold only "
            "unadjusted (Holm p = 0.11 and 0.24)"),
    "C14": ("F4", "holds for increase alarms (Holm p <= 0.001); for decrease alarms the Gaussian monitor was better out of "
            "sample (5.5 points, Holm p = 0.003)"),
    "C15": ("F6", "holds (15 of 16 cells Holm p <= 0.021; the exception is the one the paper names)"),
    "C16": ("F6", "exploratory (held-out replication, settings frozen); holds as C15"),
    "C17": ("F7", "exploratory; no significant difference; equivalence not shown in any of the 8 cells"),
    "C18": ("F7", "exploratory; holds only unadjusted (Holm p = 0.53 and 0.44)"),
    "C19": ("", "not tested by E20"),
    "C20": ("F8a, F8b", "first-run windows: no significant difference, equivalence not shown; session-long BJV worse for "
            "decrease alarms only (Holm p = 0.021 and 0.035)"),
    "C21": ("F9a", "holds (Holm p <= 0.0008; control share treated as fixed)"),
    "C22": ("F9b", "holds (held-out increases, Holm p = 0.0008)"),
    "C23": ("F11", "same day: holds in development for both directions and for held-out decreases, not for held-out "
            "increases (p = 0.062); next day: no significant difference"),
    "C24": ("F10a, F10b", "holds (Fisher Holm p = 0.0006; circular-shift Holm p = 0.0021)"),
    "C25": ("F10c", "exploratory; Holm p = 0.011"),
    "C26": ("F10a, F10b", "no difference detected (Fisher Holm p = 0.72, circular-shift p = 0.50); not an equivalence"),
    "C27": ("F10a, F10b", "holds (Fisher Holm p = 0.0095; circular-shift Holm p = 0.015)"),
    "C28": ("F10a", "descriptive; the March 2020 window did not differ from ordinary sessions (Holm p = 0.72)"),
    "C29": ("F12a", "FOMC holds (Holm p = 0.00025); CPI does not hold (p = 0.11)"),
    "C30": ("F12b", "holds (Holm p = 0.0032 and 0.00015)"),
    "C31": ("F12d", "exploratory; FOMC Holm p <= 0.0011 on all three coins (BTC CPI at awake releases p = 0.27)"),
    "C32": ("F12c", "exploratory; holds (Holm p <= 0.0004)"),
    "C33": ("", "descriptive, not tested by E20"),
    "C34": ("", "not tested by E20 (no family registered)"),
    "C35": ("F15", "holds; equivalence within 5 points shown in all 12 cells"),
    "C36": ("", "descriptive at unequal alarm rates, not tested by E20"),
}


def part_classify():
    df = pd.DataFrame(CLASSIFICATION, columns=["claim_id", "claim", "line", "source", "data_set", "class", "verifiable", "evidence"])
    assert set(df.claim_id) == set(STATUS)
    df["families"] = df.claim_id.map(lambda k: STATUS[k][0])
    df["adjusted_status"] = df.claim_id.map(lambda k: STATUS[k][1])
    df.to_csv(OUT / "e20_classification.csv", index=False)
    pe = pd.DataFrame(PRIMARY_ENDPOINTS, columns=["data_set", "primary_endpoint", "evidence"])
    pe.to_csv(OUT / "e20_primary_endpoints.csv", index=False)
    print(df[["claim_id", "line", "source", "class", "verifiable", "families"]].to_string())
    print(pe[["data_set", "primary_endpoint"]].to_string())


PARTS = {"checks": part_checks, "bootstrap": part_bootstrap, "holm": part_holm, "events": part_events,
         "forecast": part_forecast, "classify": part_classify}

if __name__ == "__main__":
    for name in sys.argv[1:] or ["classify"]:
        PARTS[name]()
