"""E6a. Tables 5 and 6 of Bapat (2026) for all five coins on Binance spot, and daily IQ.

REAL MARKET DATA: Binance spot 1-minute klines (BTCUSDT, ETHUSDT, XRPUSDT,
SOLUSDT, DOGEUSDT), 1 July to 31 December 2025, UTC, via the public mirror read
by ``avmonitor.data.load_binance``. Published numbers come from
``paper_tables.py`` ([PUBLISHED]). The paper does not name its venue, so every
statement below is "consistent with" or "inconsistent with", never "is".

Research questions
  Q6a. Table 5 side by side with the published rows. Which statistics match?
       Falsifiable hypotheses for the mismatches:
       H1 DOGE prices rounded to d decimals (or a coarser tick) give Q1 = Q3 = 0;
       H2 a stale (sparsely updated) DOGE price gives Q1 = Q3 = 0;
       H3 the ETH mean has the wrong sign because of a sign or ordering slip
          (the mean is fixed by the first and last prices);
       H4 a different calendar window (time zone) explains the small mean gaps;
       H4b the window has about 2,000 fewer returns than a complete half year
           (the length implied by Kernel / RQ in Q6b);
       H5 the extremes come from one crash (10 Oct 2025) and from venue
          microstructure, so gaps or multi-minute returns could move them.
  Q6b. Table 6 for all five coins: which printed formula or alternative gives
       each published magnitude, per coin? Is Kernel / RQ = W(30) / (n - 30)?
       Raw against scaled MAD for the robust kernel, also after adjusting for
       the MAD implied by the paper's own quartiles. Which definition of the
       "RMSE" column reproduces the published values?
  Q6c. Daily IQ on a common, documented scale (T = 1 day) per coin: ranking of
       estimators, and the share of the fourth moment removed by MAD truncation
       (raw returns, deseasonalised returns, Gaussian theory).

Outputs: results/e06a_table5.csv, e06a_hypotheses.csv, e06a_window_offsets.csv, e06a_window_length.csv,
e06a_extremes.csv, e06b_table6_formulas.csv, e06b_power_scan.csv,
e06b_kernel_rq_ratio.csv, e06b_robust_mad.csv, e06b_rmse_definitions.csv,
e06c_daily_iq.csv, e06c_summary.csv, e06c_truncation_share.csv.
Seed: 0 (only H2 uses random numbers).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from paper_tables import TABLE5_ALL, TABLE6_ALL, TABLE6_TRUE_QUARTICITY  # noqa: E402

from avmonitor import robust as R  # noqa: E402
from avmonitor.data import intraday_scale, load_binance  # noqa: E402

warnings.simplefilter("ignore", category=getattr(R, "DegenerateThresholdWarning", UserWarning))

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
COINS = {"BTC": "BTCUSDT", "ETH": "ETHUSDT", "XRP": "XRPUSDT", "SOL": "SOLUSDT", "DOGE": "DOGEUSDT"}
START, END = "2025-07-01", "2026-01-01"
THETA = TABLE6_TRUE_QUARTICITY
H = 30
mu = R._mu
CRASH = (pd.Timestamp("2025-10-10 21:00", tz="UTC"), pd.Timestamp("2025-10-10 23:00", tz="UTC"))


def load_all(start=START, end=END):
    return {c: load_binance(s, start, end, download=False) for c, s in COINS.items()}


def table5_row(r: np.ndarray) -> dict:
    q = np.quantile(r, [0, 0.25, 0.5, 0.75, 1])
    return dict(min=q[0], q1=q[1], median=q[2], mean=float(np.mean(r)), q3=q[3], max=q[4])


def printed_half_unit(x: float) -> float:
    """Half a unit in the last printed digit (3 significant figures, or 4 decimals for 0)."""
    if x == 0:
        return 5e-5
    e = np.floor(np.log10(abs(x)))
    return 0.5 * 10 ** (e - 2)


# ----------------------------------------------------------------------------
# Q6a
# ----------------------------------------------------------------------------


def q6a(data):
    rows = []
    for c, df in data.items():
        r = df.logret.dropna().to_numpy()
        ours = table5_row(r)
        for k, pub in TABLE5_ALL[c].items():
            v = ours[k]
            tol = printed_half_unit(pub)
            rows.append(dict(coin=c, stat=k, binance=v, published=pub, ratio=(v / pub) if pub != 0 else np.nan,
                             within_print_precision=bool(abs(v - pub) <= tol),
                             within_5pct=bool(pub != 0 and abs(v / pub - 1) <= 0.05) or (pub == 0 and abs(v) < 5e-5)))
        rows.append(dict(coin=c, stat="zero_share", binance=float(np.mean(r == 0)), published=np.nan, ratio=np.nan,
                         within_print_precision=np.nan, within_5pct=np.nan))
        rows.append(dict(coin=c, stat="n_returns", binance=r.size, published=np.nan, ratio=np.nan,
                         within_print_precision=np.nan, within_5pct=np.nan))
    t5 = pd.DataFrame(rows)
    t5.to_csv(OUT / "e06a_table5.csv", index=False)
    print("\n[Q6a] Table 5, Binance against published")
    print(t5.pivot(index="coin", columns="stat", values="ratio")[["min", "q1", "mean", "q3", "max"]].round(3))

    hyp = []
    # H1 and H2: DOGE price discreteness or staleness
    doge = data["DOGE"].close.to_numpy()
    def stats_for(p, label, param):
        rr = np.diff(np.log(p))
        rr = rr[np.isfinite(rr)]
        q = np.quantile(rr, [0, 0.25, 0.75, 1])
        mad_raw = R.mad(rr, scaled=False)
        rk = float(R.robust_kernel_quarticity_printed(rr, H, 3.0, mad_scaled=False))
        return dict(hypothesis=label, param=param, zero_share=float(np.mean(rr == 0)), q1=q[1], q3=q[2], min=q[0], max=q[3],
                    mean=float(rr.mean()), mad_raw=float(mad_raw), robust_kernel_raw=rk,
                    RQ=float(np.sum(rr**4) / 3), MinRQ=float(R.min_rq(rr) / rr.size), MedRQ=float(R.med_rq(rr) / rr.size),
                    q1_q3_zero=bool(q[1] == 0 and q[2] == 0))
    hyp.append(stats_for(doge, "DOGE as on Binance (tick 1e-5)", "none"))
    for d in (4, 3):
        hyp.append(stats_for(np.round(doge, d), "H1 DOGE closes rounded to d decimals", f"d={d}"))
    for tick in (2e-4, 2.5e-4, 3e-4, 4e-4, 5e-4):
        hyp.append(stats_for(np.round(doge / tick) * tick, "H1 DOGE closes on a coarser tick", f"tick={tick:g}"))
    rng = np.random.default_rng(0)
    for q_upd in (0.5, 0.4, 0.3):
        keep = rng.random(doge.size) < q_upd
        keep[0] = True
        idx = np.maximum.accumulate(np.where(keep, np.arange(doge.size), 0))
        hyp.append(stats_for(doge[idx], "H2 DOGE price refreshed only in a share q of minutes", f"q={q_upd}"))
    pub = TABLE5_ALL["DOGE"]
    p6 = TABLE6_ALL["DOGE"]
    hyp.append(dict(hypothesis="published DOGE (Tables 5, 6)", param="", zero_share=np.nan, q1=pub["q1"], q3=pub["q3"],
                    min=pub["min"], max=pub["max"], mean=pub["mean"], mad_raw=np.nan, robust_kernel_raw=p6["Robust kernel"][0],
                    RQ=p6["RQ"][0], MinRQ=p6["MinRQ"][0], MedRQ=p6["MedRQ"][0], q1_q3_zero=True))
    # Rounding check on the other four coins (a common 4-decimal rule must leave them unchanged)
    for c in ("BTC", "ETH", "XRP", "SOL"):
        p = data[c].close.to_numpy()
        hyp.append(stats_for(np.round(p, 4), "H1 control: other coin rounded to 4 decimals", c))
    hdf = pd.DataFrame(hyp)
    # H3: ETH sign; mean fixed by first and last price
    eth = data["ETH"].close
    n = int(data["ETH"].logret.notna().sum())
    implied = np.log(eth.iloc[-1] / eth.iloc[0])
    h3 = []
    for c, df in data.items():
        p = df.close
        nn = int(df.logret.notna().sum())
        h3.append(dict(coin=c, first_close=float(p.iloc[0]), last_close=float(p.iloc[-1]), n_returns=nn,
                       log_change_binance=float(np.log(p.iloc[-1] / p.iloc[0])),
                       log_change_implied_by_published_mean=TABLE5_ALL[c]["mean"] * nn,
                       ratio=TABLE5_ALL[c]["mean"] * nn / float(np.log(p.iloc[-1] / p.iloc[0]))))
    h3 = pd.DataFrame(h3)
    print("\n[Q6a H3] log change first to last close, Binance against implied by the published mean")
    print(h3.round(4).to_string(index=False))
    rr = data["ETH"].logret.dropna().to_numpy()
    print(f"ETH under a sign flip: predicted min {-rr.max():.4f} (published {TABLE5_ALL['ETH']['min']}), "
          f"predicted max {-rr.min():.4f} (published {TABLE5_ALL['ETH']['max']}); implied ETH log change {implied:.4f}, n {n}")
    hdf.to_csv(OUT / "e06a_hypotheses.csv", index=False)
    h3.to_csv(OUT / "e06a_mean_log_change.csv", index=False)
    print("\n[Q6a H1, H2] DOGE")
    print(hdf.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    # H4: calendar window shifted by a time-zone offset
    closes = {c: load_binance(s, "2025-06-29", END, download=False).close for c, s in COINS.items()}
    rows = []
    for off in np.arange(0, 30.01, 0.5):
        st = pd.Timestamp(START, tz="UTC") - pd.Timedelta(hours=off)
        en = pd.Timestamp(END, tz="UTC") - pd.Timedelta(hours=off) - pd.Timedelta(minutes=1)
        nmin = int((en - st) / pd.Timedelta(minutes=1))
        rec = dict(offset_hours_before_utc=off)
        for c in COINS:
            m = np.log(closes[c][en] / closes[c][st]) / nmin
            pubm = TABLE5_ALL[c]["mean"] * (-1 if c == "ETH" else 1)
            rec[f"{c}_ratio"] = m / pubm
        rec["max_abs_err_excl_ETH"] = max(abs(rec[f"{c}_ratio"] - 1) for c in ("BTC", "XRP", "SOL", "DOGE"))
        rows.append(rec)
    w = pd.DataFrame(rows)
    w.to_csv(OUT / "e06a_window_offsets.csv", index=False)
    print("\n[Q6a H4] best calendar offsets (ETH sign flipped); 0 = UTC, 5.5 = IST, 8 = UTC+8")
    print(w.sort_values("max_abs_err_excl_ETH").head(5).round(4).to_string(index=False))
    print(w[w.offset_hours_before_utc.isin([0, 5.5, 8])].round(4).to_string(index=False))

    # H4b: the window length implied by Kernel / RQ (Q6b) is about 262,950 returns, about 2,000 fewer than a
    # complete UTC half year. Three ways to lose them: the window ends early, starts late, or has missing minutes
    # (same first and last price, fewer returns, returns spanning the gaps).
    n_t = 262950
    t_start, t_last = pd.Timestamp(START, tz="UTC"), pd.Timestamp(END, tz="UTC") - pd.Timedelta(minutes=1)
    cases = {"same window, about 2,000 missing minutes": (t_start, t_last),
             "window ends early (30 Dec 14:30 UTC)": (t_start, t_start + pd.Timedelta(minutes=n_t)),
             "window starts late (2 Jul 09:29 UTC)": (t_last - pd.Timedelta(minutes=n_t), t_last)}
    rows = []
    for lab, (st, en) in cases.items():
        rec = dict(case=lab, n_returns=n_t, start=str(st), end=str(en))
        for c in COINS:
            m = np.log(closes[c][en] / closes[c][st]) / n_t
            rec[f"{c}_ratio"] = m / (TABLE5_ALL[c]["mean"] * (-1 if c == "ETH" else 1))
        rec["max_abs_err_all_five_ETH_flipped"] = max(abs(rec[f"{c}_ratio"] - 1) for c in COINS)
        rows.append(rec)
    wl = pd.DataFrame(rows)
    wl.to_csv(OUT / "e06a_window_length.csv", index=False)
    print("\n[Q6a H4b] means for a sample of 262,950 returns (ratio Binance / published, ETH sign flipped)")
    print(wl.round(4).to_string(index=False))

    # H5: extremes, when they happen, and multi-minute returns
    rows = []
    for c, df in data.items():
        lp = np.log(df.close)
        r1 = df.logret
        crash = (r1.index >= CRASH[0]) & (r1.index < CRASH[1])
        for k in (1, 2, 3, 5, 10):
            rk = (lp - lp.shift(k)).dropna()
            rows.append(dict(coin=c, k_minutes=k, min=rk.min(), min_time=str(rk.idxmin()), max=rk.max(), max_time=str(rk.idxmax()),
                             published_min=TABLE5_ALL[c]["min"], published_max=TABLE5_ALL[c]["max"]))
        ex = r1[~crash].dropna()
        rows.append(dict(coin=c, k_minutes="1, without 10 Oct 21:00 to 23:00", min=ex.min(), min_time=str(ex.idxmin()), max=ex.max(),
                         max_time=str(ex.idxmax()), published_min=TABLE5_ALL[c]["min"], published_max=TABLE5_ALL[c]["max"]))
    e = pd.DataFrame(rows)
    e.to_csv(OUT / "e06a_extremes.csv", index=False)
    print("\n[Q6a H5] extremes")
    print(e.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    return t5


# ----------------------------------------------------------------------------
# Q6b
# ----------------------------------------------------------------------------


def candidate_formulas(r: np.ndarray) -> dict:
    a = np.abs(r)
    n = r.size
    return {
        ("RQ", "(1/3) sum r^4 (no n)"): np.sum(r**4) / 3,
        ("RQ", "(n/3) sum r^4 (eq. 3, T = 1)"): n / 3 * np.sum(r**4),
        ("BPQ", "eq. (4) as printed: mu1^-4 sum |r_i||r_i-1|^3"): mu(1) ** -4 * np.sum(a[1:] * a[:-1] ** 3),
        ("BPQ", "sum |r_i||r_i-1| (no constant)"): np.sum(a[1:] * a[:-1]),
        ("BPQ", "mu1^-2 sum |r_i||r_i-1| (bipower variation)"): mu(1) ** -2 * np.sum(a[1:] * a[:-1]),
        ("BPQ", "(1/3) mu1^-2 sum |r_i||r_i-1|"): mu(1) ** -2 * np.sum(a[1:] * a[:-1]) / 3,
        ("TPQ", "eq. (5) as printed, p = 1/2"): mu(0.5) ** -3 * np.sum(a[2:] ** 0.5 * a[1:-1] ** 0.5 * a[:-2] ** 0.5),
        ("TPQ", "sum prod |r|^(1/2) (no constant)"): np.sum(a[2:] ** 0.5 * a[1:-1] ** 0.5 * a[:-2] ** 0.5),
        ("TPQ", "(1/3) eq. (5), p = 1/2"): mu(0.5) ** -3 * np.sum(a[2:] ** 0.5 * a[1:-1] ** 0.5 * a[:-2] ** 0.5) / 3,
        ("TPQ", "eq. (5) as printed, p = 2/3"): mu(2 / 3) ** -3 * np.sum(a[2:] ** (2 / 3) * a[1:-1] ** (2 / 3) * a[:-2] ** (2 / 3)),
        ("TPQ", "eq. (5), p = 4/3 (standard TPQ, no n)"): mu(4 / 3) ** -3 * np.sum(a[2:] ** (4 / 3) * a[1:-1] ** (4 / 3) * a[:-2] ** (4 / 3)),
        ("MinRQ", "ADS constant, no n"): R.min_rq(r) / n,
        ("MinRQ", "ADS constant, with n"): R.min_rq(r),
        ("MedRQ", "ADS constant, no n"): R.med_rq(r) / n,
        ("MedRQ", "ADS constant, with n"): R.med_rq(r),
        ("Kernel", "eq. (7) as printed, h = 30"): R.kernel_quarticity_printed(r, H),
        ("Robust kernel", "eqs. (11), (16) as printed, raw MAD"): R.robust_kernel_quarticity_printed(r, H, 3.0, mad_scaled=False),
        ("Robust kernel", "eqs. (11), (16) as printed, scaled MAD"): R.robust_kernel_quarticity_printed(r, H, 3.0, mad_scaled=True),
    }


def q6b(data):
    W = R.kernel_weight_sum(H)
    rows, rets = [], {}
    for c, df in data.items():
        r = df.logret.dropna().to_numpy()
        rets[c] = r
        for (est, f), v in candidate_formulas(r).items():
            pub = TABLE6_ALL[c][est][0]
            rows.append(dict(coin=c, estimator=est, formula=f, binance=float(v), published=pub, ratio=float(v) / pub))
        # Same formulas with the 10 Oct crash window removed (a venue-independent "bulk" check)
        ix = df.logret.dropna().index
        keep = ~((ix >= CRASH[0]) & (ix < CRASH[1]))
        for (est, f), v in candidate_formulas(r[keep]).items():
            pub = TABLE6_ALL[c][est][0]
            rows.append(dict(coin=c, estimator=est, formula=f + " [10 Oct 21:00 to 23:00 UTC removed]", binance=float(v),
                             published=pub, ratio=float(v) / pub))
    f = pd.DataFrame(rows)
    f.to_csv(OUT / "e06b_table6_formulas.csv", index=False)
    piv = f.pivot_table(index=["estimator", "formula"], columns="coin", values="ratio")[list(COINS)]
    piv["sd_log_ratio_4coins"] = np.log(piv[["BTC", "ETH", "XRP", "SOL"]]).std(axis=1)
    print("\n[Q6b] Binance / published, per coin (sd of log ratio over BTC, ETH, XRP, SOL measures cross-coin fit)")
    print(piv.to_string(float_format=lambda v: f"{v:.3g}"))

    # Power scan: which total power in |r| reproduces the cross-coin pattern of BPQ and TPQ?
    scan = []
    for est, lags in (("BPQ", 2), ("TPQ", 3)):
        for p in np.round(np.arange(1.0, 4.01, 0.25), 2):
            logs = {}
            for c, r in rets.items():
                a = np.abs(r)
                prod = np.ones(a.size - lags + 1)
                for j in range(lags):
                    prod = prod * a[j : a.size - lags + 1 + j] ** (p / lags)
                logs[c] = np.log(TABLE6_ALL[c][est][0] / prod.sum())
            v4 = np.array([logs[c] for c in ("BTC", "ETH", "XRP", "SOL")])
            const = float(np.exp(v4.mean()))
            scan.append(dict(estimator=est, total_power=p, sd_log_ratio_4coins=float(v4.std(ddof=1)),
                             implied_constant=const, mu_constant=mu(p / lags) ** -lags,
                             **{f"log_pub_over_sum_{c}": logs[c] for c in COINS}))
    scan = pd.DataFrame(scan)
    scan.to_csv(OUT / "e06b_power_scan.csv", index=False)
    print("\n[Q6b] power scan for BPQ (2 lags) and TPQ (3 lags)")
    print(scan[["estimator", "total_power", "sd_log_ratio_4coins", "implied_constant", "mu_constant"]].round(4).to_string(index=False))

    # Kernel / RQ ratio
    kr = []
    for c, r in rets.items():
        n = r.size
        k = R.kernel_quarticity_printed(r, H)
        rq = np.sum(r**4) / 3
        kr.append(dict(coin=c, n=n, W30_over_n_minus_30=W / (n - H), binance_kernel_over_rq=float(k / rq),
                       published_kernel_over_rq=TABLE6_ALL[c]["Kernel"][0] / TABLE6_ALL[c]["RQ"][0],
                       published_kernel_over_rq_lo=(TABLE6_ALL[c]["Kernel"][0] - printed_half_unit(TABLE6_ALL[c]["Kernel"][0]))
                       / (TABLE6_ALL[c]["RQ"][0] + printed_half_unit(TABLE6_ALL[c]["RQ"][0])),
                       published_kernel_over_rq_hi=(TABLE6_ALL[c]["Kernel"][0] + printed_half_unit(TABLE6_ALL[c]["Kernel"][0]))
                       / (TABLE6_ALL[c]["RQ"][0] - printed_half_unit(TABLE6_ALL[c]["RQ"][0]))))
    kr = pd.DataFrame(kr)
    kr["W30_inside_published_rounding_interval"] = (kr.W30_over_n_minus_30 >= kr.published_kernel_over_rq_lo) & (
        kr.W30_over_n_minus_30 <= kr.published_kernel_over_rq_hi)
    # If Kernel / RQ = W(30) / (n - 30) exactly, the published ratio pins the sample length n.
    kr["n_implied_lo"] = W / kr.published_kernel_over_rq_hi + H
    kr["n_implied_hi"] = W / kr.published_kernel_over_rq_lo + H
    kr["n_implied_mid"] = W / kr.published_kernel_over_rq + H
    lo, hi = kr.n_implied_lo.max(), kr.n_implied_hi.min()
    print(f"W(30) = {W:.4f}. Sample length n consistent with all five published ratios: {lo:.0f} to {hi:.0f} "
          f"(Binance has {kr.n.iloc[0]} returns)" if lo <= hi else "No common n fits all five published ratios")
    kr.to_csv(OUT / "e06b_kernel_rq_ratio.csv", index=False)
    print("\n[Q6b] Kernel / RQ against W(30) / (n - 30)")
    print(kr.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # Robust kernel: raw against scaled MAD, and adjusted to the paper's own quartiles
    rm = []
    for c, r in rets.items():
        q1, q3 = np.quantile(r, [0.25, 0.75])
        half_iqr_ours = (q3 - q1) / 2
        half_iqr_pub = (TABLE5_ALL[c]["q3"] - TABLE5_ALL[c]["q1"]) / 2
        kern = R.kernel_quarticity_printed(r, H)
        for scaled in (False, True):
            v = float(R.robust_kernel_quarticity_printed(r, H, 3.0, mad_scaled=scaled))
            adj = v * (half_iqr_pub / half_iqr_ours) ** 4
            pub = TABLE6_ALL[c]["Robust kernel"][0]
            rm.append(dict(coin=c, mad="scaled" if scaled else "raw", binance=v, published=pub, ratio=v / pub,
                           half_iqr_binance=half_iqr_ours, half_iqr_published=half_iqr_pub,
                           ratio_after_quartile_adjustment=adj / pub,
                           robust_over_kernel_binance=v / kern,
                           robust_over_kernel_published=pub / TABLE6_ALL[c]["Kernel"][0]))
    rm = pd.DataFrame(rm)
    rm.to_csv(OUT / "e06b_robust_mad.csv", index=False)
    print("\n[Q6b] Robust kernel, raw against scaled MAD")
    print(rm.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    # RMSE definitions
    k4s = R.truncated_moment_factor(3.0, 4)
    k4r = R.truncated_moment_factor(R.effective_cutoff(3.0, False), 4)
    rows = []
    for c, df in data.items():
        r = rets[c]
        days = [g.dropna().to_numpy() for _, g in df.logret.groupby(df.index.floor("D"))]
        days = [d for d in days if d.size >= 1400]
        cand_day = [candidate_formulas(d) for d in days]
        local_k = R.local_quarticity_printed(r, H)
        local_rk = R.local_quarticity_printed(R.truncate(r, R.mad_threshold(r, 3.0, scaled=False)), H)
        main = {"RQ": ("RQ", "(1/3) sum r^4 (no n)"), "BPQ": ("BPQ", "mu1^-2 sum |r_i||r_i-1| (bipower variation)"),
                "TPQ": ("TPQ", "eq. (5) as printed, p = 1/2"), "MinRQ": ("MinRQ", "ADS constant, no n"),
                "MedRQ": ("MedRQ", "ADS constant, no n"), "Kernel": ("Kernel", "eq. (7) as printed, h = 30"),
                "Robust kernel": ("Robust kernel", "eqs. (11), (16) as printed, raw MAD")}
        full = candidate_formulas(r)
        for est, key in main.items():
            q_pub, rmse_pub = TABLE6_ALL[c][est]
            dvals = np.array([cd[key] for cd in cand_day])
            rec = dict(coin=c, estimator=est, published_estimate=q_pub, published_rmse=rmse_pub,
                       D1_abs_error_published_estimate=abs(q_pub - THETA),
                       D1_abs_error_binance_estimate=abs(full[key] - THETA),
                       D2_rmse_over_daily_estimates_binance=float(np.sqrt(np.mean((dvals - THETA) ** 2))),
                       D3_rmse_over_local_estimates_binance=np.nan,
                       D4_abs_error_vs_kappa4_scaled_theta=abs(q_pub - k4s * THETA),
                       D4_abs_error_vs_kappa4_raw_theta=abs(q_pub - k4r * THETA))
            if est == "Kernel":
                rec["D3_rmse_over_local_estimates_binance"] = float(np.sqrt(np.mean((local_k - THETA) ** 2)))
            if est == "Robust kernel":
                rec["D3_rmse_over_local_estimates_binance"] = float(np.sqrt(np.mean((local_rk - THETA) ** 2)))
            rows.append(rec)
    d = pd.DataFrame(rows)

    def consistent(target_factor):
        # The published estimate is rounded to 3 significant figures, so |Q - f theta| lies in an interval.
        # The published RMSE R is taken as consistent if it equals that interval after rounding OR truncating
        # to 3 significant figures, i.e. if [R - half unit, R + 2 half units] overlaps it.
        hq = d.published_estimate.map(printed_half_unit)
        a = np.abs(d.published_estimate - hq - target_factor * THETA)
        b = np.abs(d.published_estimate + hq - target_factor * THETA)
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        hr = d.published_rmse.map(printed_half_unit)
        return (d.published_rmse - hr <= hi) & (d.published_rmse + 2 * hr >= lo)

    d["D1_consistent_with_published"] = consistent(1.0)
    d["D4_scaled_consistent_with_published"] = consistent(k4s)
    d["D4_raw_consistent_with_published"] = consistent(k4r)
    # Any RMSE that averages squared errors of estimates Q_j is at least |mean(Q_j) - theta| (Jensen).
    # To print R < theta, the estimates would need a mean of at least theta - R.
    d["mean_estimate_needed_for_averaged_rmse"] = np.where(d.published_rmse < THETA, THETA - d.published_rmse, np.nan)
    d.to_csv(OUT / "e06b_rmse_definitions.csv", index=False)
    print(d.groupby("estimator", sort=False)[["D1_consistent_with_published", "D4_scaled_consistent_with_published",
                                                "D4_raw_consistent_with_published"]].sum())
    print(f"\n[Q6b] RMSE definitions (theta = {THETA}, kappa4 scaled = {k4s:.4f}, raw = {k4r:.4f})")
    print(d.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return f


# ----------------------------------------------------------------------------
# Q6c
# ----------------------------------------------------------------------------


def q6c():
    k4s = R.truncated_moment_factor(3.0, 4)
    k4r = R.truncated_moment_factor(R.effective_cutoff(3.0, False), 4)
    rows, trunc = [], []
    for c, s in COINS.items():
        df = load_binance(s, "2025-06-01", END, download=False)  # June only warms up the intraday profile
        sc = intraday_scale(df.logret, lookback_days=20)
        df = df[df.index >= pd.Timestamp(START, tz="UTC")]
        sc = sc[sc.index >= pd.Timestamp(START, tz="UTC")]
        u = df.logret / sc
        for day, g in df.groupby(df.index.floor("D")):
            r = g.logret.dropna().to_numpy()
            if r.size < 1400:
                continue
            ud = u.loc[g.index].dropna().to_numpy()
            rec = dict(coin=c, day=day, n=r.size, zero_share=float(np.mean(r == 0)), RV=float(np.sum(r**2)),
                       RQ=R.realized_quarticity(r), TPQ=R.tripower_quarticity(r), QPQ=R.quadpower_quarticity(r),
                       MinRQ=R.min_rq(r), MedRQ=R.med_rq(r), TRQ=R.threshold_rq(r), Kernel=R.kernel_iq(r, H),
                       RobustKernel_scaled_bc=R.kernel_iq(r, H, c=3, mad_scaled=True, bias_correct=True),
                       RobustKernel_scaled=R.kernel_iq(r, H, c=3, mad_scaled=True),
                       RobustKernel_raw=R.kernel_iq(r, H, c=3, mad_scaled=False))
            rows.append(rec)
            for lab, x in (("raw returns", r), ("deseasonalised returns", ud)):
                s4 = np.sum(x**4)
                for scaled in (False, True):
                    xt = R.truncate(x, R.mad_threshold(x, 3.0, scaled=scaled))
                    trunc.append(dict(coin=c, day=day, returns=lab, mad="scaled" if scaled else "raw", mad_window="one day",
                                      share_4th_kept=float(np.sum(xt**4) / s4), share_clipped=float(np.mean(np.abs(x) > np.max(np.abs(xt))))))
        # Full-sample MAD over the six months (as in the paper)
        r = df.logret.dropna().to_numpy()
        ud = u.dropna().to_numpy()
        for lab, x in (("raw returns", r), ("deseasonalised returns", ud)):
            for scaled in (False, True):
                xt = R.truncate(x, R.mad_threshold(x, 3.0, scaled=scaled))
                trunc.append(dict(coin=c, day=pd.NaT, returns=lab, mad="scaled" if scaled else "raw", mad_window="six months",
                                  share_4th_kept=float(np.sum(xt**4) / np.sum(x**4)),
                                  share_clipped=float(np.mean(np.abs(x) > np.max(np.abs(xt))))))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "e06c_daily_iq.csv", index=False)
    cols = ["RQ", "Kernel", "TPQ", "QPQ", "MinRQ", "MedRQ", "TRQ", "RobustKernel_scaled_bc", "RobustKernel_scaled", "RobustKernel_raw"]
    summ = []
    for c, g in d.groupby("coin", sort=False):
        for k in cols:
            summ.append(dict(coin=c, estimator=k, days=len(g), median_daily_iq=float(g[k].median()),
                             median_ratio_to_MedRQ=float((g[k] / g.MedRQ).median()),
                             q25_ratio_to_MedRQ=float((g[k] / g.MedRQ).quantile(0.25)),
                             q75_ratio_to_MedRQ=float((g[k] / g.MedRQ).quantile(0.75)),
                             median_ratio_to_RQ=float((g[k] / g.RQ).median())))
    s = pd.DataFrame(summ)
    s["rank_within_coin"] = s.groupby("coin").median_ratio_to_MedRQ.rank(ascending=False).astype(int)
    s.to_csv(OUT / "e06c_summary.csv", index=False)
    print("\n[Q6c] median ratio of daily IQ to MedRQ")
    print(s.pivot(index="estimator", columns="coin", values="median_ratio_to_MedRQ")[list(COINS)].loc[cols].round(3))
    t = pd.DataFrame(trunc)
    t.to_csv(OUT / "e06c_truncation_share.csv", index=False)
    ts = t.groupby(["coin", "returns", "mad", "mad_window"], sort=False).share_4th_kept.median().unstack(["returns", "mad_window"])
    print(f"\n[Q6c] median share of the fourth moment kept after clipping at 3 MAD (Gaussian theory: raw {k4r:.3f}, scaled {k4s:.3f})")
    print(ts.round(3))
    return d, s, t


def main():
    data = load_all()
    for c, df in data.items():
        ch = df.attrs["checks"]
        print(f"REAL DATA {c}: {ch['rows_read']} rows, missing minutes {ch['missing_minutes']}, duplicates {ch['duplicates']}")
    q6a(data)
    q6b(data)
    q6c()


if __name__ == "__main__":
    main()
