"""Risk switch experiment: change-triggered exposure against volatility targeting.

REAL MARKET DATA [REAL]. Binance spot one-minute bars for BTCUSDT, ETHUSDT,
XRPUSDT, SOLUSDT and DOGEUSDT (January 2024 to December 2025) and Nifty 50
index one-minute bars (April 2017 to October 2026), loaded with
``avmonitor.data`` (``download=False``; see ``scripts/download_data.py``).

Research questions
  RQ1 (risk control). Does a long position whose size is reset only when the
       anytime-valid monitor raises an alarm (exposure = target volatility over
       the new regime's robust volatility) deliver risk closer to the target,
       smaller drawdowns and milder worst days than buy and hold and than classic
       EWMA volatility targeting refreshed every minute, hour or day?
  RQ2 (costs). After transaction costs, how do turnover and the cost drag of the
       switch compare with the EWMA rules, and how sensitive is the ranking to
       the cost assumption (0 to 20 basis points per unit of turnover for crypto,
       0 to 10 for Nifty)?
  RQ3 (driver). Does the robust e-detector make a better trigger than the
       Gaussian likelihood-ratio monitor, which alarms on jumps?
  RQ4 (out of sample). With every setting frozen on 2024 crypto, do the answers
       hold on 2025 crypto and on nine years of Nifty 50?

Protocol, fixed before the out-of-sample runs
  * Exposure decided at the close of minute t applies to minute t + 1.
  * Targets: 40% a year for crypto, 15% for Nifty. Exposure capped at 2.
  * Base cost: 5 bp per unit of turnover for crypto, 2.5 bp for Nifty.
  * Tuning on 2024 crypto only (window 1 Feb to 31 Dec 2024, after a warm-up),
    one rule for every family: among the variants whose cost drag at the base
    cost is at most 2% a year, take the one with the lowest mean tracking RMS
    across the five coins; if no variant of a family is within budget, take the
    lowest tracking RMS anyway and flag it as over budget. The conversion from
    robust scale to standard deviation is the median daily ratio on the same
    window. Then everything is frozen.
  * Monitor: A = bars per session / (false alarm budget per session), with
    budgets 5% and 1% in the grid; calibration 1440 bars; c = 3.
  * Sharpe intervals: circular block bootstrap, blocks of 10 sessions, 2000
    resamples, seed 0.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import backtest as B  # noqa: E402
import strategies as S  # noqa: E402
from avmonitor.data import load_binance, load_nifty, regular_sessions  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
FIG = OUT / "figures"

CONFIG = dict(
    symbols=["BTCUSDT", "ETHUSDT", "XRPUSDT", "SOLUSDT", "DOGEUSDT"],
    crypto_load=("2024-01-01", "2026-01-01"),
    is_window=("2024-02-01", "2025-01-01"),
    oos_window=("2025-01-01", "2026-01-01"),
    nifty_load=("2017-04-01", "2026-11-01"),
    nifty_window=("2017-06-01", "2026-11-01"),
    target=dict(crypto=0.40, nifty=0.15),
    w_max=2.0,
    cost_base=dict(crypto=5.0, nifty=2.5),
    cost_grid=dict(crypto=[0.0, 2.0, 5.0, 10.0, 20.0], nifty=[0.0, 1.0, 2.5, 5.0, 10.0]),
    ewma_halflives_days=[0.25, 1.0, 3.0, 10.0, 30.0],
    ewma_updates=["minute", "hour", "day"],
    robust_windows=[1440, 4320],
    bands=[0.25, 0.5, 0.75],
    false_alarm_budgets=[0.05, 0.01],
    cal_len=1440,
    fast_window=120,
    max_cost_drag=0.02,
    block=10,
    n_boot=2000,
    n_boot_grid=500,
    seed=0,
    profile_days=20,
)

PALETTE = dict(switch="#2a78d6", ewma_day="#eb6834", ewma_hour="#1baf7a", ewma_minute="#eda100",
               gaussian="#e87ba4", robust_roll="#6250d6", buy_and_hold="#777777", target="#333333")
LABELS = dict(buy_and_hold="buy and hold", ewma_minute="EWMA, minute", ewma_hour="EWMA, hour", ewma_day="EWMA, day",
              robust_roll="robust rolling, day", switch_robust="switch, robust monitor", switch_gaussian="switch, Gaussian LR")
COLORS = dict(buy_and_hold=PALETTE["buy_and_hold"], ewma_minute=PALETTE["ewma_minute"], ewma_hour=PALETTE["ewma_hour"],
              ewma_day=PALETTE["ewma_day"], robust_roll=PALETTE["robust_roll"], switch_robust=PALETTE["switch"],
              switch_gaussian=PALETTE["gaussian"])


def window_mask(P: S.Prepared, start: str, end: str) -> np.ndarray:
    tz = P.index.tz
    lo, hi = pd.Timestamp(start, tz=tz), pd.Timestamp(end, tz=tz)
    return np.asarray((P.index >= lo) & (P.index < hi))


def build_exposures(P: S.Prepared, k: float, target: float, cfg: dict) -> tuple[dict, dict]:
    """Every exposure rule in the grid. Returns ``{name: (w, meta)}`` and the monitor events."""
    ex = {}
    w_max = cfg["w_max"]
    ex["buy_and_hold"] = (S.buy_and_hold(P), dict(family="buy_and_hold"))
    for upd in cfg["ewma_updates"]:
        for hl in cfg["ewma_halflives_days"]:
            ex[f"ewma_{upd}_hl{hl:g}d"] = (S.ewma_exposure(P, hl, upd, target, w_max),
                                           dict(family="ewma", update=upd, halflife_days=hl))
    for win in cfg["robust_windows"]:
        ex[f"robust_roll_day_{win}"] = (S.robust_rolling_exposure(P, win, k, target, w_max, "day"),
                                        dict(family="robust_roll", update="day", window_bars=win))
    events = {}
    for fa in cfg["false_alarm_budgets"]:
        A = P.spec.bars_per_session / fa
        for band in cfg["bands"]:
            ev = S.monitor_events(P.u, S.make_monitor("robust", A, cfg["cal_len"], band))
            events[("robust", fa, band)] = ev
            for prov in (True, False):
                name = f"switch_robust_fa{fa:g}_b{band:g}_{'prov' if prov else 'hold'}"
                ex[name] = (S.switch_exposure(P, ev, k, target, w_max, cfg["fast_window"], prov),
                            dict(family="switch", kind="robust", false_alarm=fa, band=band, provisional=prov))
        ev = S.monitor_events(P.u, S.make_monitor("gaussian", A, cfg["cal_len"], 0.0))
        events[("gaussian", fa, 0.0)] = ev
        for prov in (True, False):
            name = f"switch_gaussian_fa{fa:g}_{'prov' if prov else 'hold'}"
            ex[name] = (S.switch_exposure(P, ev, k, target, w_max, cfg["fast_window"], prov),
                        dict(family="switch", kind="gaussian", false_alarm=fa, band=np.nan, provisional=prov))
    return ex, events


def evaluate(P, ex, names, window, target, cost, n_boot, cfg, asset, period) -> pd.DataFrame:
    rows = []
    for name in names:
        w, meta = ex[name]
        bt = B.run_backtest(P, w, cost)
        tgt = None if meta["family"] == "buy_and_hold" else target
        m = B.metrics(P, bt, window, tgt, cost, cfg["block"], n_boot, cfg["seed"])
        rows.append(dict(asset=asset, period=period, strategy=name, **meta, **m))
    return pd.DataFrame(rows)


def select(is_df: pd.DataFrame, cfg: dict) -> tuple[dict, pd.DataFrame]:
    """Apply the pre-stated selection rule to the 2024 crypto grid (mean over coins)."""
    agg = is_df.groupby("strategy").agg(
        family=("family", "first"), update=("update", "first"), halflife_days=("halflife_days", "first"),
        window_bars=("window_bars", "first"), kind=("kind", "first"), false_alarm=("false_alarm", "first"),
        band=("band", "first"), provisional=("provisional", "first"),
        tracking_rms=("tracking_rms", "mean"), cost_drag_per_year=("cost_drag_per_year", "mean"),
        turnover_per_year=("turnover_per_year", "mean"), trades_per_year=("trades_per_year", "mean"),
        sharpe=("sharpe", "mean"), max_drawdown=("max_drawdown", "mean"), worst_day=("worst_day", "mean"),
        realised_vol=("realised_vol", "mean"), mean_exposure=("mean_exposure", "mean"),
    )
    agg["eligible"] = agg["cost_drag_per_year"] <= cfg["max_cost_drag"]
    chosen, over_budget = {}, {}

    def best(key, mask):
        c = agg[mask & agg["eligible"]]
        over_budget[key] = c.empty
        if c.empty:
            c = agg[mask]
        return c["tracking_rms"].idxmin()

    for upd in cfg["ewma_updates"]:
        chosen[f"ewma_{upd}"] = best(f"ewma_{upd}", (agg["family"] == "ewma") & (agg["update"] == upd))
    chosen["robust_roll"] = best("robust_roll", agg["family"] == "robust_roll")
    chosen["switch_robust"] = best("switch_robust", (agg["family"] == "switch") & (agg["kind"] == "robust"))
    row = agg.loc[chosen["switch_robust"]]
    prov = bool(row["provisional"])
    chosen["switch_gaussian"] = f"switch_gaussian_fa{row['false_alarm']:g}_{'prov' if prov else 'hold'}"
    over_budget["switch_gaussian"] = not bool(agg.loc[chosen["switch_gaussian"], "eligible"])
    chosen["buy_and_hold"] = "buy_and_hold"
    over_budget["buy_and_hold"] = False
    agg["chosen"] = agg.index.isin(list(chosen.values()))
    return chosen, over_budget, agg.reset_index()


def alarm_stats(events: dict, P: S.Prepared, window: np.ndarray, asset: str, period: str) -> list[dict]:
    idx = np.flatnonzero(window)
    lo, hi = idx[0], idx[-1]
    n_sess = len(np.unique(P.session[window]))
    years = n_sess / P.spec.sessions_per_year
    rows = []
    for (kind, fa, band), ev in events.items():
        al = [e for e in ev if e["kind"] == "alarm" and lo <= e["index"] <= hi]
        rows.append(dict(asset=asset, period=period, monitor=kind, false_alarm=fa, band=band, alarms=len(al),
                         up=sum(e["direction"] == "up" for e in al), down=sum(e["direction"] == "down" for e in al),
                         alarms_per_year=len(al) / years))
    return rows


def paired_differences(P, ex, chosen, ref_key, window, target, cost, cfg, n_boot, asset, period) -> list[dict]:
    ref = B.daily_returns(P, B.run_backtest(P, ex[chosen[ref_key]][0], cost), window)
    rows = []
    for key, name in chosen.items():
        d = B.daily_returns(P, B.run_backtest(P, ex[name][0], cost), window)
        ci = B.sharpe_ci(d, P.spec.sessions_per_year, cfg["block"], n_boot, cfg["seed"], other=ref)
        rows.append(dict(asset=asset, period=period, strategy=name, reference=chosen[ref_key],
                         sharpe=B.sharpe(d, P.spec.sessions_per_year), **ci))
    return rows


def ew_portfolio(Ps: list, exs: list, name: str, cost: float) -> B.BacktestResult:
    """Equal-weight portfolio of one strategy across coins (same minute index for every coin)."""
    bts = [B.run_backtest(P, ex[name][0], cost) for P, ex in zip(Ps, exs)]
    return B.BacktestResult(np.mean([b.net for b in bts], axis=0), np.mean([b.gross for b in bts], axis=0),
                            np.mean([b.dw for b in bts], axis=0), np.mean([b.w for b in bts], axis=0))


def run_crypto(cfg: dict, quick: bool):
    symbols = cfg["symbols"][:1] if quick else cfg["symbols"]
    target = cfg["target"]["crypto"]
    cost = cfg["cost_base"]["crypto"]
    Ps, checks = {}, {}
    t0 = time.time()
    for sym in symbols:
        df = load_binance(sym, *cfg["crypto_load"], download=False)
        Ps[sym] = S.prepare(df, S.CRYPTO, cfg["profile_days"])
        checks[sym] = df.attrs["checks"]
        print(f"{sym}: {df.attrs['checks']}", flush=True)
    is_win = {s: window_mask(P, *cfg["is_window"]) for s, P in Ps.items()}
    oos_win = {s: window_mask(P, *cfg["oos_window"]) for s, P in Ps.items()}
    ks = {s: S.robust_to_sd_factor(P, is_win[s]) for s, P in Ps.items()}
    k = float(np.median(list(ks.values())))
    print(f"robust_to_sd per coin {ks}, pooled {k:.4f}", flush=True)

    exs, events = {}, {}
    for sym, P in Ps.items():
        exs[sym], events[sym] = build_exposures(P, k, target, cfg)
        print(f"{sym}: exposures built ({time.time() - t0:.0f} s)", flush=True)

    # 1. In-sample grid on 2024
    n_boot_grid = 100 if quick else cfg["n_boot_grid"]
    is_rows = [evaluate(Ps[s], exs[s], list(exs[s]), is_win[s], target, cost, n_boot_grid, cfg, s, "2024_IS")
               for s in symbols]
    is_df = pd.concat(is_rows, ignore_index=True)
    is_df.to_csv(OUT / "is_grid_crypto_2024.csv", index=False)
    chosen, over_budget, agg = select(is_df, cfg)
    agg.to_csv(OUT / "is_selection_crypto_2024.csv", index=False)
    print("chosen:", chosen, "over budget:", over_budget, flush=True)

    # 2. Out of sample 2025 with the frozen set, per coin and equal-weight portfolio
    n_boot = 300 if quick else cfg["n_boot"]
    names = list(chosen.values())
    oos_rows = [evaluate(Ps[s], exs[s], names, oos_win[s], target, cost, n_boot, cfg, s, "2025_OOS") for s in symbols]
    is_rows_frozen = [evaluate(Ps[s], exs[s], names, is_win[s], target, cost, n_boot, cfg, s, "2024_IS") for s in symbols]
    P0 = Ps[symbols[0]]
    for period, win in (("2025_OOS", oos_win[symbols[0]]), ("2024_IS", is_win[symbols[0]])):
        rows = []
        for key, name in chosen.items():
            bt = ew_portfolio(list(Ps.values()), list(exs.values()), name, cost)
            meta = exs[symbols[0]][name][1]
            m = B.metrics(P0, bt, win, None if meta["family"] == "buy_and_hold" else target, cost, cfg["block"], n_boot, cfg["seed"])
            rows.append(dict(asset="EW5", period=period, strategy=name, **meta, **m))
        (oos_rows if period == "2025_OOS" else is_rows_frozen).append(pd.DataFrame(rows))
    frozen = pd.concat(is_rows_frozen + oos_rows, ignore_index=True)
    frozen["role"] = frozen["strategy"].map({v: kk for kk, v in chosen.items()})
    frozen["over_budget"] = frozen["role"].map(over_budget)
    frozen.to_csv(OUT / "frozen_crypto.csv", index=False)

    # 3. Paired Sharpe differences against daily EWMA targeting (out of sample)
    diffs = []
    for s in symbols:
        diffs += paired_differences(Ps[s], exs[s], chosen, "ewma_day", oos_win[s], target, cost, cfg, n_boot, s, "2025_OOS")
    ref = B.daily_returns(P0, ew_portfolio(list(Ps.values()), list(exs.values()), chosen["ewma_day"], cost), oos_win[symbols[0]])
    for key, name in chosen.items():
        d = B.daily_returns(P0, ew_portfolio(list(Ps.values()), list(exs.values()), name, cost), oos_win[symbols[0]])
        ci = B.sharpe_ci(d, S.CRYPTO.sessions_per_year, cfg["block"], n_boot, cfg["seed"], other=ref)
        diffs.append(dict(asset="EW5", period="2025_OOS", strategy=name, reference=chosen["ewma_day"],
                          sharpe=B.sharpe(d, S.CRYPTO.sessions_per_year), **ci))
    pd.DataFrame(diffs).to_csv(OUT / "sharpe_diff_crypto_2025.csv", index=False)

    # 4. Cost sensitivity (equal-weight portfolio, out of sample)
    rows = []
    for c in cfg["cost_grid"]["crypto"]:
        for key, name in chosen.items():
            bt = ew_portfolio(list(Ps.values()), list(exs.values()), name, c)
            meta = exs[symbols[0]][name][1]
            m = B.metrics(P0, bt, oos_win[symbols[0]], None if meta["family"] == "buy_and_hold" else target, c, cfg["block"], 0, cfg["seed"])
            rows.append(dict(asset="EW5", period="2025_OOS", strategy=name, role=key, family=meta["family"], **m))
    pd.DataFrame(rows).to_csv(OUT / "cost_sensitivity_crypto_2025.csv", index=False)

    # 5. Alarm counts
    al = []
    for s in symbols:
        al += alarm_stats(events[s], Ps[s], is_win[s], s, "2024_IS") + alarm_stats(events[s], Ps[s], oos_win[s], s, "2025_OOS")
    pd.DataFrame(al).to_csv(OUT / "alarms_crypto.csv", index=False)

    # 6. Per-minute exposure paths of the frozen set for BTC (figures) and the frozen parameters
    btc = symbols[0]
    meta = exs[btc][chosen["switch_robust"]][1]
    ev_btc = [e for e in events[btc][("robust", meta["false_alarm"], meta["band"])] if e["kind"] == "alarm"]
    paths = pd.DataFrame({key: exs[btc][name][0] for key, name in chosen.items()}, index=Ps[btc].index)
    paths["close"] = load_binance(btc, *cfg["crypto_load"], download=False)["close"].to_numpy()
    paths = paths[oos_win[btc]]
    sel = dict(chosen=chosen, over_budget=over_budget, robust_to_sd=k, robust_to_sd_per_coin=ks,
               data_checks=checks, config=cfg)
    with open(OUT / "frozen_parameters.json", "w") as fh:
        json.dump(sel, fh, indent=2, default=str)
    return chosen, over_budget, k, paths, [Ps[btc].index[e["index"]] for e in ev_btc]


def run_nifty(cfg: dict, chosen: dict, over_budget: dict, k: float, quick: bool):
    target = cfg["target"]["nifty"]
    cost = cfg["cost_base"]["nifty"]
    load = ("2024-01-01", "2026-11-01") if quick else cfg["nifty_load"]
    win = ("2024-03-01", "2026-11-01") if quick else cfg["nifty_window"]
    df = regular_sessions(load_nifty(*load, download=False))
    print(f"NIFTY: {df.attrs['checks']}", flush=True)
    P = S.prepare(df, S.NIFTY, cfg["profile_days"])
    window = window_mask(P, *win)
    ex, events = build_exposures(P, k, target, cfg)
    names = list(chosen.values())
    n_boot = 300 if quick else cfg["n_boot"]
    res = evaluate(P, ex, names, window, target, cost, n_boot, cfg, "NIFTY50", "2017_2026_OOS")
    res["role"] = res["strategy"].map({v: kk for kk, v in chosen.items()})
    res["over_budget"] = res["role"].map(over_budget)
    res.to_csv(OUT / "frozen_nifty.csv", index=False)
    full = evaluate(P, ex, list(ex), window, target, cost, 0, cfg, "NIFTY50", "2017_2026_OOS")
    full.to_csv(OUT / "grid_nifty_not_used_for_selection.csv", index=False)
    diffs = paired_differences(P, ex, chosen, "ewma_day", window, target, cost, cfg, n_boot, "NIFTY50", "2017_2026_OOS")
    pd.DataFrame(diffs).to_csv(OUT / "sharpe_diff_nifty.csv", index=False)
    rows = []
    for c in cfg["cost_grid"]["nifty"]:
        for key, name in chosen.items():
            bt = B.run_backtest(P, ex[name][0], c)
            meta = ex[name][1]
            m = B.metrics(P, bt, window, None if meta["family"] == "buy_and_hold" else target, c, cfg["block"], 0, cfg["seed"])
            rows.append(dict(asset="NIFTY50", period="2017_2026_OOS", strategy=name, role=key, family=meta["family"], **m))
    pd.DataFrame(rows).to_csv(OUT / "cost_sensitivity_nifty.csv", index=False)
    pd.DataFrame(alarm_stats(events, P, window, "NIFTY50", "2017_2026_OOS")).to_csv(OUT / "alarms_nifty.csv", index=False)
    # yearly breakdown of the frozen set (risk metrics only, no bootstrap)
    years = []
    for y in range(int(win[0][:4]), int(win[1][:4]) + 1):
        wy = window & window_mask(P, f"{y}-01-01", f"{y + 1}-01-01")
        if wy.sum() < 60 * P.spec.bars_per_session:
            continue
        yr = evaluate(P, ex, names, wy, target, cost, 0, cfg, "NIFTY50", str(y))
        yr["role"] = yr["strategy"].map({v: kk for kk, v in chosen.items()})
        years.append(yr)
    pd.concat(years, ignore_index=True).to_csv(OUT / "frozen_nifty_by_year.csv", index=False)
    # weekly realised vol paths for the figure
    weekly = {}
    for key in ("buy_and_hold", "switch_robust", "ewma_day"):
        bt = B.run_backtest(P, ex[chosen[key]][0], cost)
        weekly[key] = B.weekly_realised_vol(bt.net[window], P.week[window], P.session[window], P.spec)
    return weekly, P, window


def figures(cfg, chosen, paths, btc_alarms, nifty_weekly, nifty_P, nifty_window):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5})
    FIG.mkdir(parents=True, exist_ok=True)

    # Figure 1: BTC 2025 exposure paths
    fig, ax = plt.subplots(3, 1, figsize=(11, 7.5), sharex=True)
    hourly = paths.resample("1h").last()
    ax[0].plot(hourly.index, hourly["close"], color=PALETTE["buy_and_hold"], lw=1.0)
    ax[0].set_ylabel("BTC/USDT close")
    for t in btc_alarms:
        if paths.index[0] <= t <= paths.index[-1]:
            ax[0].axvline(t, color=PALETTE["switch"], lw=0.4, alpha=0.35)
    ax[0].set_title("BTC 2025 [REAL]: price with robust-monitor alarms (thin blue lines), then exposures")
    ax[1].step(hourly.index, hourly["switch_robust"], color=PALETTE["switch"], lw=1.2, where="post", label="change-triggered switch (robust monitor)")
    ax[1].step(hourly.index, hourly["switch_gaussian"], color=PALETTE["gaussian"], lw=0.8, where="post", alpha=0.8, label="switch driven by Gaussian LR monitor")
    ax[1].set_ylabel("exposure")
    ax[1].legend(loc="upper right", frameon=False)
    ax[2].step(hourly.index, hourly["ewma_day"], color=PALETTE["ewma_day"], lw=1.2, where="post", label="EWMA targeting, daily refresh")
    ax[2].step(hourly.index, hourly["ewma_hour"], color=PALETTE["ewma_hour"], lw=0.8, where="post", alpha=0.8, label="EWMA targeting, hourly refresh")
    ax[2].set_ylabel("exposure")
    ax[2].legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_btc_2025_exposures.png", dpi=130)
    plt.close(fig)

    # Figure 2: weekly realised volatility against the target, Nifty
    fig, ax = plt.subplots(figsize=(11, 4))
    wk_sessions = pd.Series(nifty_P.session[nifty_window]).groupby(nifty_P.week[nifty_window]).nunique()
    wk_idx = pd.Series(nifty_P.index[nifty_window]).groupby(nifty_P.week[nifty_window]).first()
    keep = (wk_sessions >= 3).to_numpy()
    x = wk_idx.to_numpy()[keep]
    ax.plot(x, nifty_weekly["buy_and_hold"], color=PALETTE["buy_and_hold"], lw=0.8, label="buy and hold")
    ax.plot(x, nifty_weekly["ewma_day"], color=PALETTE["ewma_day"], lw=0.9, label="EWMA targeting, daily refresh")
    ax.plot(x, nifty_weekly["switch_robust"], color=PALETTE["switch"], lw=0.9, label="change-triggered switch")
    ax.axhline(cfg["target"]["nifty"], color=PALETTE["target"], lw=1.0, ls="--", label="target 15%")
    ax.set_yscale("log")
    ticks = [0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{100 * v:.0f}%"))
    ax.set_ylabel("weekly realised volatility (annualised)")
    ax.set_title("Nifty 50 [REAL]: delivered risk week by week (parameters frozen on 2024 crypto)")
    ax.legend(loc="upper right", frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_nifty_weekly_vol.png", dpi=130)
    plt.close(fig)

    # Figure 3: cost sensitivity
    cs_c = pd.read_csv(OUT / "cost_sensitivity_crypto_2025.csv")
    cs_n = pd.read_csv(OUT / "cost_sensitivity_nifty.csv")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for a, cs, ttl in ((ax[0], cs_c, "Crypto, equal-weight 5 coins, 2025 (out of sample)"),
                       (ax[1], cs_n, "Nifty 50, 2017 to 2026 (out of sample)")):
        for role in LABELS:
            g = cs[cs.role == role].sort_values("cost_bps")
            a.plot(g.cost_bps, g.sharpe, marker="o", ms=4, lw=1.2, color=COLORS[role], label=LABELS[role])
        a.set_xlabel("one-way cost per unit of turnover (basis points)")
        a.set_ylabel("net Sharpe ratio")
        a.set_title(ttl, fontsize=9)
    ax[1].legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_cost_sensitivity.png", dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="one coin, short Nifty sample, fewer bootstrap draws")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    chosen, over_budget, k, paths, btc_alarms = run_crypto(CONFIG, args.quick)
    print(f"crypto done in {time.time() - t0:.0f} s", flush=True)
    nifty_weekly, P, window = run_nifty(CONFIG, chosen, over_budget, k, args.quick)
    print(f"nifty done in {time.time() - t0:.0f} s", flush=True)
    figures(CONFIG, chosen, paths, btc_alarms, nifty_weekly, P, window)
    print(f"all done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
