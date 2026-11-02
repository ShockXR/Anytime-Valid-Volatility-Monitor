"""E15. Monitoring a panel of assets with false discovery control (tracker B74).

Theory and proofs: docs/theory_panel.md. Code: src/avmonitor/panel.py.

Research questions
  Q15a [SYNTHETIC] Validity. K = 5, 10, 20, 50 assets, one block of 390 steps (a
       "session"), returns either independent or strongly dependent (a common factor
       with correlation 0.7 between every pair, and jumps that hit every asset at the
       same minute), no change or a variance increase in 20% or 60% of the assets.
       Is the false discovery rate (FDR) of the e-BH panel monitor at most alpha =
       5%: at the end of the block, at the first discovery, and in the worst false
       discovery proportion over the block (proved bound alpha (1 + log(1/alpha)),
       about 20%)? What do per-asset Bonferroni monitors (Shiryaev-Roberts at
       A x K) and naive per-asset monitors (at A) give?
  Q15b [SYNTHETIC] Power and delay. How many more changes does e-BH find than
       Bonferroni, how much earlier, and how does this depend on K, on the share of
       assets that change and on whether they change at the same minute?
  Q15c [SYNTHETIC] Horizon. Over ten blocks (3,900 steps): one "no change since the
       start" e-process per asset (Pareto weights over the change time, one e-BH
       guarantee for the whole horizon) against daily blocks restarted every 390 steps
       (one guarantee per block). Family-wise error over the horizon, FDR, and the
       delay as a function of how late the change comes (the log(1/w) cost).
  Q15d [REAL] Five Binance coins (BTC, ETH, XRP, SOL, DOGE), one-minute returns
       divided by the predictable intraday profile, one block per UTC day, each coin
       calibrated on its previous 1,440 returns, c = 3, band 0.25 (the E5/E6 primary
       setting) and 0.5 (the other band of the E6 grid), alpha = 5% per day. Jul to Dec 2025 (the E6 period) and Jan 2024 to
       Jun 2025 (the E6 out-of-sample period). Discoveries per rule, discoveries that
       only e-BH's adaptivity allowed, and validation as in E6: is volatility over
       the next 1,440 returns above (up) or below (down) the calibration baseline,
       against matched random times; rule sets nested (Bonferroni, then e-BH, then
       naive) and a rate-matched per-asset threshold chosen with hindsight.
  Q15e [REAL] A mixed panel on the UTC minute clock: the five coins plus the Nifty 50
       index (regular sessions only, 03:45 to 10:00 UTC; its returns are missing
       outside the session, and it has its own clock: horizon 374 returns per
       session, calibration on its previous 750 returns as in E7), Feb 2024 to Dec
       2025. Same questions as Q15d.

Settings are fixed before any run: alpha = 0.05, c = 3, the bet grids of
detectors.py, window weights over one block. Real data: bands 0.25 and 0.5 (the
E6 grid without band 0), plug-in baselines (as E6). Synthetic: the exact null mean (oracle
baseline), computed in closed form, so that the multiple testing layer is tested
on its own; one robustness row uses the confidence-sequence baseline.

Seeds: 15 (Q15a/b), 151 (Q15c), 152 (real-data controls). Usage:
  python experiments/e15_panel.py synthetic   # Q15a, Q15b  (about 4 minutes on 2 cores)
  python experiments/e15_panel.py horizon     # Q15c        (about 1 minute)
  python experiments/e15_panel.py real        # Q15d, Q15e  (needs data/raw; about 4 minutes)
  python experiments/e15_panel.py figures
"""

from __future__ import annotations

import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from avmonitor.detectors import betting_edetector, two_sided_alarm
from avmonitor.panel import calibrate_panel, ebh, fdp_summary, run_block
from avmonitor.robust import MAD_NORMAL_SCALE

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"

ALPHA = 0.05
TAU = 3.0  # synthetic: unit scale, c = 3
M = 390  # synthetic block length
A = M / ALPHA  # per-asset Shiryaev-Roberts level for alpha per block (Theorem B of docs/math.md)
P_JUMP, S_JUMP = 0.002, 4.0  # jump probability per step and jump sd (in units of the diffusive sd)
RHO = 0.7  # factor correlation in the dependent design
NU = 100  # change time in the "common" design
CHUNK = 200


def kappa2(a):
    """E[min(Z^2, a^2)] for a standard normal Z (docs/math.md, section 3)."""
    return (2 * stats.norm.cdf(a) - 1) - 2 * a * stats.norm.pdf(a) + 2 * a * a * stats.norm.sf(a)


def null_mean(tau=TAU, p=P_JUMP, sj=S_JUMP):
    """Exact E[min(r^2, tau^2)] for r = Z + B J, B ~ Bernoulli(p), J ~ N(0, sj^2), Z ~ N(0, 1)."""
    s2 = 1 + sj**2
    return (1 - p) * kappa2(tau) + p * s2 * kappa2(tau / np.sqrt(s2))


M0 = null_mean()


# ----------------------------------------------------------------------------
# Synthetic panels
# ----------------------------------------------------------------------------


def simulate(rng, reps, K, n, dep, k1, ratio, timing):
    """Returns (reps, K, n) and change times (reps, K); n means "no change"."""
    if dep == "independent":
        z = rng.standard_normal((reps, K, n))
        jump = (rng.random((reps, K, n)) < P_JUMP) * rng.standard_normal((reps, K, n)) * S_JUMP
    elif dep == "factor+cojumps":
        f = rng.standard_normal((reps, 1, n))
        z = np.sqrt(RHO) * f + np.sqrt(1 - RHO) * rng.standard_normal((reps, K, n))
        jump = (rng.random((reps, 1, n)) < P_JUMP) * rng.standard_normal((reps, 1, n)) * S_JUMP
    else:
        raise ValueError(dep)
    nu = np.full((reps, K), n)
    if k1:
        if timing == "common":
            nu[:, :k1] = NU
        else:
            nu[:, :k1] = rng.integers(50, 301, size=(reps, k1))
        z = z * np.where(np.arange(n)[None, None, :] >= nu[:, :, None], ratio, 1.0)
    return z + jump, nu


def sr_per_asset(r, level):
    """First two-sided alarm of the per-asset Shiryaev-Roberts e-detector (each side at 2 x level)."""
    reps, K, n = r.shape
    x = np.minimum(r.reshape(reps * K, n) ** 2, TAU**2)
    up = betting_edetector(x, M0, "up", threshold=2 * level)
    dn = betting_edetector(x, M0, "down", upper=TAU**2, threshold=2 * level)
    return two_sided_alarm(up, dn).reshape(reps, K)


def metrics(ft, nu, n, k1):
    s = fdp_summary(ft, nu, n)
    reps = ft.shape[0]
    true = (ft >= 0) & (ft >= nu)
    out = dict(
        fwer=float((s["V"] >= 1).mean()),
        fdr_end=float(s["fdp_end"].mean()),
        fdr_end_se=float(s["fdp_end"].std(ddof=1) / np.sqrt(reps)),
        fdr_first=float(np.nan_to_num(s["fdp_first"], nan=0.0).mean()),
        sup_fdp=float(s["sup_fdp"].mean()),
        mean_V=float(s["V"].mean()),
        mean_R=float(s["R"].mean()),
    )
    if k1:
        d = (ft - nu + 1)[:, :k1][true[:, :k1]]
        out.update(power=float(true[:, :k1].mean()), median_delay=float(np.median(d)) if d.size else np.nan,
                   mean_delay=float(d.mean()) if d.size else np.nan)
    else:
        out.update(power=np.nan, median_delay=np.nan, mean_delay=np.nan)
    return out


def synthetic_job(job):
    K, dep, pi1, timing, ratio, reps, seed = job
    rng = np.random.default_rng(seed)
    k1 = int(round(pi1 * K))
    fts = {m: [] for m in ("naive SR (A)", "Bonferroni SR (A x K)", "panel Bonferroni", "panel e-BH")}
    nus = []
    for start in range(0, reps, CHUNK):
        b = min(CHUNK, reps - start)
        r, nu = simulate(rng, b, K, M, dep, k1, ratio, timing)
        nus.append(nu)
        fts["naive SR (A)"].append(sr_per_asset(r, A))
        fts["Bonferroni SR (A x K)"].append(sr_per_asset(r, A * K))
        fts["panel Bonferroni"].append(run_block(r, TAU**2, M0, M0, ALPHA, "bonferroni").flag_time)
        fts["panel e-BH"].append(run_block(r, TAU**2, M0, M0, ALPHA, "ebh").flag_time)
    nu = np.concatenate(nus)
    rows = []
    for name, parts in fts.items():
        ft = np.concatenate(parts)
        rows.append(dict(K=K, dependence=dep, share_changed=pi1, changed=k1, timing=timing if k1 else "none",
                         sd_ratio=ratio if k1 else np.nan, reps=reps, method=name, **metrics(ft, nu, M, k1)))
    return rows


def run_synthetic():
    jobs, seed = [], 15
    designs = [(0.0, "none", 1.0), (0.2, "common", 1.25), (0.2, "staggered", 1.25), (0.6, "common", 1.25),
               (0.6, "staggered", 1.25), (0.2, "common", 1.5), (0.6, "common", 1.5)]
    for K in (5, 10, 20, 50):
        for dep in ("independent", "factor+cojumps"):
            for pi1, timing, ratio in designs:
                seed += 1
                jobs.append((K, dep, pi1, timing, ratio, 2000 if pi1 == 0 else 600, seed))
    with Pool(2) as pool:
        rows = [r for res in pool.map(synthetic_job, jobs) for r in res]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e15a_synthetic.csv", index=False)
    print(df.to_string(max_rows=400))
    robustness_cs()


def robustness_cs():
    """One design with ESTIMATED baselines (confidence-sequence baseline, delta = 0.01 split over K)."""
    from avmonitor.panel import calibrate_panel as cp

    rng = np.random.default_rng(1500)
    rows = []
    for K, dep, pi1 in ((20, "factor+cojumps", 0.0), (20, "factor+cojumps", 0.2), (20, "independent", 0.2)):
        k1 = int(round(pi1 * K))
        reps, cal = 400, 1440
        ft_e, ft_b, nus = [], [], []
        for _ in range(reps // CHUNK):
            rc, _ = simulate(rng, CHUNK, K, cal, dep, 0, 1.0, "none")
            r, nu = simulate(rng, CHUNK, K, M, dep, k1, 1.25, "common")
            tau2 = np.ones((CHUNK, K))
            m_up = np.ones((CHUNK, K))
            m_dn = np.ones((CHUNK, K))
            act = np.zeros((CHUNK, K), bool)
            for i in range(CHUNK):
                b = cp(list(rc[i]), cal, c=3.0, null_band=0.0, baseline="cs", delta=0.01)
                tau2[i], m_up[i], m_dn[i], act[i] = b.tau2, b.m_up, b.m_dn, b.active
            ft_e.append(run_block(r, tau2, m_up, m_dn, ALPHA, "ebh", active=act).flag_time)
            ft_b.append(run_block(r, tau2, m_up, m_dn, ALPHA, "bonferroni", active=act).flag_time)
            nus.append(nu)
        nu = np.concatenate(nus)
        for name, ft in (("panel e-BH, cs baseline", ft_e), ("panel Bonferroni, cs baseline", ft_b)):
            rows.append(dict(K=K, dependence=dep, share_changed=pi1, changed=k1, timing="common" if k1 else "none",
                             sd_ratio=1.25 if k1 else np.nan, reps=reps, method=name,
                             **metrics(np.concatenate(ft), nu, M, k1)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e15a_cs_baseline.csv", index=False)
    print(df.to_string())


# ----------------------------------------------------------------------------
# Q15c: one guarantee for the whole horizon against one per block
# ----------------------------------------------------------------------------


def run_horizon():
    rng = np.random.default_rng(151)
    K, n_blocks, reps = 20, 10, 400
    N = n_blocks * M
    rows, delay_rows = [], []
    for k1 in (0, 6):
        fts = {m: [] for m in ("daily blocks, e-BH", "daily blocks, Bonferroni", "since start, Pareto h = 390, e-BH",
                               "since start, window h = 3900, e-BH")}
        nus = []
        for _ in range(reps // CHUNK):
            r, _ = simulate(rng, CHUNK, K, N, "independent", 0, 1.0, "none")
            nu = np.full((CHUNK, K), N)
            if k1:
                nu[:, :k1] = rng.integers(1, N, size=(CHUNK, k1))
                r = np.where(np.arange(N)[None, None, :] >= nu[:, :, None], 1.5 * r, r)
            nus.append(nu)
            for rule, name in (("ebh", "daily blocks, e-BH"), ("bonferroni", "daily blocks, Bonferroni")):
                first = np.full((CHUNK, K), -1)
                for b in range(n_blocks):
                    ft = run_block(r[:, :, b * M : (b + 1) * M], TAU**2, M0, M0, ALPHA, rule).flag_time
                    new = (first < 0) & (ft >= 0)
                    first[new] = ft[new] + b * M
                fts[name].append(first)
            fts["since start, Pareto h = 390, e-BH"].append(
                run_block(r, TAU**2, M0, M0, ALPHA, "ebh", horizon=M, weights="pareto").flag_time)
            fts["since start, window h = 3900, e-BH"].append(
                run_block(r, TAU**2, M0, M0, ALPHA, "ebh", horizon=N, weights="window").flag_time)
        nu = np.concatenate(nus)
        for name, parts in fts.items():
            ft = np.concatenate(parts)
            rows.append(dict(changed=k1, method=name, reps=reps, **metrics(ft, nu, N, k1)))
            if k1:
                ch = nu[:, :k1].ravel()
                f = ft[:, :k1].ravel()
                for lo, hi in ((0, 390), (390, 1300), (1300, 2600), (2600, N)):
                    sel = (ch >= lo) & (ch < hi)
                    ok = sel & (f >= ch)
                    d = (f - ch + 1)[ok]
                    delay_rows.append(dict(method=name, change_from=lo, change_to=hi, n_changes=int(sel.sum()),
                                           detected=float(ok.sum() / max(sel.sum(), 1)),
                                           median_delay=float(np.median(d)) if d.size else np.nan))
    pd.DataFrame(rows).to_csv(OUT / "e15c_horizon.csv", index=False)
    pd.DataFrame(delay_rows).to_csv(OUT / "e15c_delay_by_change_time.csv", index=False)
    print(pd.DataFrame(rows).to_string())
    print(pd.DataFrame(delay_rows).to_string())


# ----------------------------------------------------------------------------
# Real data
# ----------------------------------------------------------------------------

COINS = {"BTC": "BTCUSDT", "ETH": "ETHUSDT", "XRP": "XRPUSDT", "SOL": "SOLUSDT", "DOGE": "DOGEUSDT"}
CRYPTO_H, CRYPTO_CAL = 1440, 1440
NIFTY_H, NIFTY_CAL = 374, 750  # 374 returns per regular session; calibration as E7's primary setting
NIFTY_LOOKBACK, NIFTY_SMOOTH = 20, 5  # E7's session profile settings, but strictly past-only here
N_CONTROL = 5


def crypto_u(symbol: str, start: str, end: str) -> pd.Series:
    from avmonitor.data import intraday_scale, load_binance

    df = load_binance(symbol, start, end, download=False)
    s = intraday_scale(df.logret, lookback_days=20)  # past-only by default
    u = df.logret / s
    return u.where(np.isfinite(u))


def nifty_u(start: str, end: str) -> pd.Series:
    """Nifty returns divided by a past-only session profile (E7's profile without its start-up look-ahead)."""
    from avmonitor.data import load_nifty, regular_sessions

    df = regular_sessions(load_nifty(start, end, download=False))
    minute = df.index.hour * 60 + df.index.minute
    a = df["logret"].abs()
    wide = pd.DataFrame({"a": a.values, "s": df["session"].values, "m": minute}).pivot_table(
        index="s", columns="m", values="a", aggfunc="first", dropna=False)
    prof = wide.rolling(NIFTY_LOOKBACK, min_periods=NIFTY_LOOKBACK // 2).median().shift(1)
    prof = prof.T.rolling(NIFTY_SMOOTH, center=True, min_periods=1).mean().T
    prof = prof.div(prof.median(axis=1), axis=0)
    long = prof.stack(future_stack=True)
    s = long.reindex(pd.MultiIndex.from_arrays([df["session"].values, minute])).to_numpy()
    u = pd.Series(df["logret"].to_numpy() / s, index=df.index.tz_convert("UTC"))
    return u.where(np.isfinite(u))


def panel_array(series: dict, start: str, end: str):
    idx = pd.date_range(pd.Timestamp(start, tz="UTC"), pd.Timestamp(end, tz="UTC"), freq="1min", inclusive="left")
    U = np.column_stack([s.reindex(idx).to_numpy() for s in series.values()])
    return idx, U


def offline_sticky_ebh(path: np.ndarray, dir_path: np.ndarray, alpha: float):  # check of run_block
    """Sticky e-BH from the unstopped e-value path of one block: (flag step, direction) per asset."""
    n, K = path.shape
    frozen = np.zeros(K, bool)
    fval = np.ones(K)
    ft = np.full(K, -1)
    fd = np.zeros(K, int)
    for t in range(n):
        v = np.where(frozen, fval, path[t])
        if v.max() < 1 / alpha:  # nothing can clear 1 / alpha: skip the sort
            continue
        new = ebh(v, alpha) & ~frozen
        if new.any():
            ft[new], fd[new], fval[new] = t, dir_path[t, new], v[new]
            frozen |= new
    return ft, fd


def run_panel(names, U, idx, block_from: str, cal_len, horizon, band: float, n_check: int = 40):
    """Daily UTC blocks from ``block_from``; calibration on each asset's previous ``cal_len`` finite returns.

    All blocks are run at once (blocks are independent given their baselines). Returns a frame of
    block-level records per asset, the unstopped e-value paths (blocks, 1440, K) and directions. The
    sticky e-BH flags of the vectorised run are checked against a step-by-step offline e-BH on the
    unstopped paths for the first ``n_check`` blocks.
    """
    K = len(names)
    cal_len = np.broadcast_to(np.asarray(cal_len, int), (K,))
    horizon = np.broadcast_to(np.asarray(horizon, float), (K,))
    days = pd.date_range(pd.Timestamp(block_from, tz="UTC"), idx[-1].floor("D"), freq="D")
    starts = idx.get_indexer(days)
    starts = starts[(starts >= 0) & (starts + 1440 <= len(idx))]
    fin = [np.nonzero(np.isfinite(U[:, k]))[0] for k in range(K)]
    B = len(starts)
    tau2, m_up, m_dn, scale = np.ones((B, K)), np.ones((B, K)), np.ones((B, K)), np.zeros((B, K))
    act = np.zeros((B, K), bool)
    for b, s in enumerate(starts):
        windows = []
        for k in range(K):
            j = np.searchsorted(fin[k], s)
            windows.append(U[fin[k][max(0, j - cal_len[k]) : j], k])
        base = calibrate_panel(windows, cal_len, c=3.0, null_band=band)
        tau2[b], m_up[b], m_dn[b], scale[b], act[b] = base.tau2, base.m_up, base.m_dn, base.scale, base.active
    seg = np.stack([U[s : s + 1440].T for s in starts])  # (B, K, 1440)
    res = run_block(seg, tau2, m_up, m_dn, ALPHA, "none", horizon=horizon[None], active=act, return_path="all")
    eb = run_block(seg, tau2, m_up, m_dn, ALPHA, "ebh", horizon=horizon[None], active=act)
    for b in range(min(n_check, B)):
        ft, fd = offline_sticky_ebh(res.path[b], res.dir_path[b], ALPHA)
        if np.any(ft != eb.flag_time[b]) or np.any(fd[ft >= 0] != eb.direction[b][ft >= 0]):
            raise RuntimeError(f"offline e-BH differs from the vectorised run in block {b}")
    recs = []
    for b, s in enumerate(starts):
        for k in range(K):
            recs.append(dict(block=b, start=s, end=s + 1440, day=idx[s].date(), asset=names[k], k=k,
                             active=bool(act[b, k]), scale=scale[b, k], n_obs=int(np.isfinite(seg[b, k]).sum()),
                             block_max=float(np.nanmax(res.path[b, :, k])), ebh_t=int(eb.flag_time[b, k]),
                             ebh_dir=int(eb.direction[b, k])))
    return pd.DataFrame(recs), res.path, res.dir_path, starts


def robust_sd(x):
    x = x[np.isfinite(x)]
    return float(np.median(np.abs(x)) * MAD_NORMAL_SCALE) if x.size else np.nan


def post_ratio(U, fin_k, k, t, h, scale):
    """Robust sd of the next h finite returns of asset k after step t, over the block's baseline scale."""
    j = np.searchsorted(fin_k, t, side="right")
    if j + h > fin_k.size or not scale > 0:
        return np.nan
    return robust_sd(U[fin_k[j : j + h], k]) / scale


def wilson(k, n, z=1.959964):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def validate(panel_name, period, band, names, U, recs, path, dirp, horizon, rng):
    """Discoveries per rule, nested categories, timing, validation against matched random times."""
    K = len(names)
    fin = [np.nonzero(np.isfinite(U[:, k]))[0] for k in range(K)]
    active = recs[recs.active].reset_index(drop=True)
    bi, ki = active.block.to_numpy(), active.k.to_numpy()
    P = path[bi, :, ki]  # (n_active, 1440) unstopped e-values of each active asset-day

    def first(level):
        hit = P >= level
        return np.where(hit.any(1), hit.argmax(1), -1)

    rules = {"e-BH": active.ebh_t.to_numpy(), "naive (1/alpha)": first(1 / ALPHA),
             "Bonferroni (K/alpha)": first(K / ALPHA)}
    n_ebh = int((rules["e-BH"] >= 0).sum())
    bm = np.sort(active.block_max.to_numpy())[::-1]
    c_star = float(bm[n_ebh - 1]) if n_ebh else np.inf
    rm_name = "rate-matched threshold (hindsight)"
    rules[rm_name] = first(c_star)
    disc = []
    for rule, tt in rules.items():
        for i in np.nonzero(tt >= 0)[0]:
            row = active.iloc[i]
            gt = int(row.start + tt[i])
            d = int(row.ebh_dir) if rule == "e-BH" else int(dirp[row.block, tt[i], row.k])
            ratio = post_ratio(U, fin[row.k], row.k, gt, int(horizon[row.k]), row.scale)
            disc.append(dict(panel=panel_name, period=period, band=band, rule=rule, asset=row.asset, day=row.day,
                             step=gt, offset=int(tt[i]), direction="up" if d > 0 else "down", post_ratio=ratio,
                             evalue=float(path[row.block, tt[i], row.k])))
    disc = pd.DataFrame(disc)
    # matched random times: same asset, same offset in a random active block of that asset
    ctrl = []
    by_asset = {a: g for a, g in active.groupby("asset")}
    for row in disc[disc.rule == "e-BH"].itertuples():
        g = by_asset[row.asset]
        k = names.index(row.asset)
        for _ in range(N_CONTROL):
            rb = g.iloc[rng.integers(len(g))]
            ctrl.append(dict(asset=row.asset, direction=row.direction,
                             post_ratio=post_ratio(U, fin[k], k, int(rb.start + row.offset), int(horizon[k]),
                                                   rb.scale)))
    ctrl = pd.DataFrame(ctrl)
    te, tb, tn = rules["e-BH"], rules["Bonferroni (K/alpha)"], rules["naive (1/alpha)"]
    in_b, in_e = tb >= 0, te >= 0
    disc["category"] = ""
    key = {(a, d): (eb_, bb_) for a, d, eb_, bb_ in zip(active.asset, active.day, in_e, in_b)}
    flags = np.array([key.get((a, d), (False, False)) for a, d in zip(disc.asset, disc.day)], bool)
    is_e = (disc.rule == "e-BH").to_numpy()
    disc.loc[is_e & flags[:, 1], "category"] = "e-BH and Bonferroni"
    disc.loc[is_e & ~flags[:, 1], "category"] = "e-BH only (not Bonferroni)"
    disc.loc[(disc.rule == "naive (1/alpha)").to_numpy() & ~flags[:, 0], "category"] = "naive only (not e-BH)"
    summ = []
    groups = [(r, disc[disc.rule == r]) for r in rules] + [
        (c, disc[disc.category == c]) for c in ("e-BH and Bonferroni", "e-BH only (not Bonferroni)",
                                                 "naive only (not e-BH)")]
    n_days = int(recs.block.nunique())
    for label, g in groups:
        row = dict(panel=panel_name, period=period, band=band, set=label, days=n_days, discoveries=len(g),
                   per_day=len(g) / n_days)
        for d in ("up", "down"):
            gd = g[(g.direction == d) & g.post_ratio.notna()]
            hits = int((gd.post_ratio > 1).sum() if d == "up" else (gd.post_ratio < 1).sum())
            lo, hi = wilson(hits, len(gd))
            row.update({f"{d}_n": len(gd), f"{d}_hits": hits, f"{d}_hit": hits / len(gd) if len(gd) else np.nan,
                        f"{d}_lo": lo, f"{d}_hi": hi})
            if label == "e-BH" and len(ctrl):
                cd = ctrl[(ctrl.direction == d) & ctrl.post_ratio.notna()]
                row[f"{d}_control"] = float(((cd.post_ratio > 1) if d == "up" else (cd.post_ratio < 1)).mean())
        summ.append(row)
    summ = pd.DataFrame(summ)
    eb = disc[disc.rule == "e-BH"]
    per_day = eb.groupby("day").size()
    both_eb = in_e & in_b
    both_ne = in_e & (tn >= 0)
    struct = dict(panel=panel_name, period=period, band=band, days=n_days, active_asset_days=len(active),
                  inactive_asset_days=int((~recs.active).sum()), ebh=int(in_e.sum()), bonferroni=int(in_b.sum()),
                  naive=int((tn >= 0).sum()), ebh_only=int((in_e & ~in_b).sum()),
                  naive_only=int(((tn >= 0) & ~in_e).sum()),
                  ebh_below_K_over_alpha=int((eb.evalue < K / ALPHA).sum()),
                  median_minutes_ebh_before_bonferroni=float(np.median((tb - te)[both_eb])) if both_eb.any() else np.nan,
                  mean_minutes_ebh_before_bonferroni=float(np.mean((tb - te)[both_eb])) if both_eb.any() else np.nan,
                  share_ebh_strictly_earlier=float(np.mean((tb - te)[both_eb] > 0)) if both_eb.any() else np.nan,
                  median_minutes_naive_before_ebh=float(np.median((te - tn)[both_ne])) if both_ne.any() else np.nan,
                  days_with_ebh_discovery=int(per_day.size), days_with_2plus=int((per_day >= 2).sum()),
                  max_per_day=int(per_day.max()) if per_day.size else 0, rate_matched_threshold=c_star)
    return disc, summ, struct


def two_prop(k1, n1, k2, n2):
    if min(n1, n2) == 0:
        return np.nan
    p = (k1 + k2) / (n1 + n2)
    se = np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    return float(2 * stats.norm.sf(abs(k1 / n1 - k2 / n2) / se)) if se > 0 else np.nan


def in_session(t: pd.DatetimeIndex) -> np.ndarray:
    """True for UTC minutes that carry a Nifty return (bars 03:46 to 09:59 UTC, that is 09:16 to 15:29 IST)."""
    m = t.hour * 60 + t.minute
    return np.asarray((m >= 226) & (m <= 599))


def run_real():
    rng = np.random.default_rng(152)
    crypto = {c: crypto_u(sym, "2024-01-01", "2026-01-01") for c, sym in COINS.items()}
    nifty = nifty_u("2023-11-01", "2026-01-01")
    mixed = dict(crypto)
    mixed["NIFTY"] = nifty
    idx_c, U_c = panel_array(crypto, "2024-01-01", "2026-01-01")
    idx_m, U_m = panel_array(mixed, "2024-01-01", "2026-01-01")
    runs = [("five coins", "Jul-Dec 2025", idx_c, U_c, "2025-07-01", "2026-01-01", 1440, 1440),
            ("five coins", "Jan 2024-Jun 2025 (out of sample)", idx_c, U_c, "2024-01-22", "2025-07-01", 1440, 1440),
            ("five coins + Nifty", "Feb 2024-Dec 2025", idx_m, U_m, "2024-02-01", "2026-01-01",
             np.array([CRYPTO_CAL] * 5 + [NIFTY_CAL]), np.array([CRYPTO_H] * 5 + [NIFTY_H]))]
    all_disc, all_summ, all_struct = [], [], []
    for band in (0.25, 0.5):  # the E6 grid; 0.25 is the E5/E6 primary setting
        for panel, period, idx, U, b0, b1, cal, hz in runs:
            sel = idx < pd.Timestamp(b1, tz="UTC")
            names = list(mixed) if panel.endswith("Nifty") else list(crypto)
            hz = np.broadcast_to(np.asarray(hz), (len(names),))
            recs, path, dirp, starts = run_panel(names, U[sel], idx[sel], b0, cal, hz, band)
            disc, summ, struct = validate(panel, period, band, names, U[sel], recs, path, dirp, hz, rng)
            all_disc.append(disc), all_summ.append(summ), all_struct.append(struct)
            print(summ.round(3).to_string())
            print(struct, flush=True)
    disc = pd.concat(all_disc, ignore_index=True)
    summ = pd.concat(all_summ, ignore_index=True)
    struct = pd.DataFrame(all_struct)
    nifty_rows = []
    for (band, rule), g in disc[disc.asset == "NIFTY"].groupby(["band", "rule"]):
        for d in ("up", "down"):
            gd = g[(g.direction == d) & g.post_ratio.notna()]
            hits = int((gd.post_ratio > 1).sum() if d == "up" else (gd.post_ratio < 1).sum())
            lo, hi = wilson(hits, len(gd))
            nifty_rows.append(dict(band=band, rule=rule, direction=d, n=len(gd), hits=hits,
                                   hit=hits / len(gd) if len(gd) else np.nan, lo=lo, hi=hi,
                                   discovered_outside_session=int(((g.direction == d).to_numpy() & ~in_session(idx_m[g.step])).sum())))
    pd.DataFrame(nifty_rows).to_csv(OUT / "e15e_nifty_discoveries.csv", index=False)
    print(pd.DataFrame(nifty_rows).to_string())
    eb = disc[(disc.rule == "e-BH") & (disc.panel == "five coins + Nifty")]
    co = (eb.groupby(["band", "day"]).asset.apply(lambda s: "+".join(sorted(s)))
          .groupby(level=0).value_counts().rename("days"))
    co.to_csv(OUT / "e15e_codiscoveries.csv")
    tests = []
    for (panel, period, band), g in summ.groupby(["panel", "period", "band"]):
        e = g[g.set == "e-BH"].iloc[0]
        rm = g[g.set.str.startswith("rate-matched")].iloc[0]
        for d in ("up", "down"):
            tests.append(dict(panel=panel, period=period, band=band, direction=d, ebh_n=e[f"{d}_n"],
                              ebh_hit=e[f"{d}_hit"], matched_n=rm[f"{d}_n"], matched_hit=rm[f"{d}_hit"],
                              p_value=two_prop(e[f"{d}_hits"], e[f"{d}_n"], rm[f"{d}_hits"], rm[f"{d}_n"])))
    disc.to_csv(OUT / "e15d_discoveries.csv", index=False)
    summ.to_csv(OUT / "e15d_validation.csv", index=False)
    struct.to_csv(OUT / "e15d_structure.csv", index=False)
    pd.DataFrame(tests).to_csv(OUT / "e15d_rate_matched_tests.csv", index=False)
    print(pd.DataFrame(tests).to_string())


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------


def figures():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(OUT / "e15a_synthetic.csv")
    colors = {"naive SR (A)": "#9e9e9e", "Bonferroni SR (A x K)": "#1f77b4", "panel Bonferroni": "#6baed6",
              "panel e-BH": "#d62728"}
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
    panels = [("factor+cojumps", 0.0, "none", "fwer", "P(any false discovery), no change"),
              ("factor+cojumps", 0.6, "common", "fdr_end", "FDR at block end, 60% change together"),
              ("factor+cojumps", 0.6, "common", "sup_fdp", "E[max over block of FDP], 60% together"),
              ("factor+cojumps", 0.2, "common", "power", "Power, 20% change together"),
              ("factor+cojumps", 0.6, "common", "power", "Power, 60% change together"),
              ("factor+cojumps", 0.6, "common", "median_delay", "Median delay (steps), 60% together")]
    for ax, (dep, pi1, timing, col, title) in zip(axes.ravel(), panels):
        sel = (df.dependence == dep) & np.isclose(df.share_changed, pi1) & (df.timing == timing)
        if pi1 > 0:
            sel &= np.isclose(df.sd_ratio, 1.25)
        for m, c in colors.items():
            g = df[sel & (df.method == m)].sort_values("K")
            ax.plot(g.K, g[col], marker="o", color=c, label=m, lw=1.8)
        if col in ("fwer", "fdr_end"):
            ax.axhline(ALPHA, color="k", ls="--", lw=1)
        if col == "sup_fdp":
            ax.axhline(ALPHA * (1 + np.log(1 / ALPHA)), color="k", ls=":", lw=1)
            ax.axhline(ALPHA, color="k", ls="--", lw=1)
        ax.set_title(title, fontsize=10)
        ax.set_xscale("log")
        ax.set_xticks([5, 10, 20, 50], ["5", "10", "20", "50"])
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("number of assets K")
    axes[0, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle("E15 [SYNTHETIC]: panel of K assets, common factor (correlation 0.7) and common jumps, "
                 "sd x 1.25 at the change, alpha = 5% per 390-step block", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "e15_synthetic.png", dpi=130)
    plt.close(fig)

    disc = pd.read_csv(OUT / "e15d_discoveries.csv", parse_dates=["day"])
    d = disc[(disc.panel == "five coins + Nifty") & (disc.rule == "e-BH")]
    names = ["BTC", "ETH", "XRP", "SOL", "DOGE", "NIFTY"]
    fig, ax = plt.subplots(figsize=(13, 3.6))
    for i, a in enumerate(names):
        g = d[d.asset == a]
        below = g.evalue < len(names) / ALPHA
        for sub, mk, lab in ((g[~below], "|", "e-value >= K / alpha (Bonferroni would also flag)"),
                             (g[below], "x", "e-value < K / alpha (only e-BH flags)")):
            up = sub.direction == "up"
            ax.scatter(sub.day[up], np.full(up.sum(), i + 0.12), marker=mk, color="#d62728", s=30,
                       label=f"up, {lab}" if i == 0 else None)
            ax.scatter(sub.day[~up], np.full((~up).sum(), i - 0.12), marker=mk, color="#1f77b4", s=30,
                       label=f"down, {lab}" if i == 0 else None)
    ax.set_yticks(range(len(names)), names)
    ax.set_ylim(-0.6, len(names) - 0.4)
    ax.grid(axis="x", alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    ax.set_title("E15 [REAL]: e-BH discoveries in the mixed panel (five Binance coins and Nifty 50), one block per "
                 "UTC day, alpha = 5% per day", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "e15_real_discoveries.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "synthetic"
    {"synthetic": run_synthetic, "horizon": run_horizon, "real": run_real, "figures": figures,
     "cs": robustness_cs}[what]()
