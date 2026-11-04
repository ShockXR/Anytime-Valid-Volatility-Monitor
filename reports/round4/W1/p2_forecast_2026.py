"""W1 part 2: does a plug-in alarm improve a HAR-I forecast of next-hour or next-day realised variance on
1 January to 24 October 2026? See PLAN.md. E20's helpers (ols_dk, qlike, cbb_weights) are imported unchanged;
E20's grid, design and expanding-window fit are copied here and adapted (horizons 60 and one day, Nifty sessions,
alarm dummies as extra regressors). Nothing after 24 October 2026 is read; no event times are used."""
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments"))
from e20_inference import ols_dk, qlike, cbb_weights  # noqa: E402

SCR = Path(__file__).resolve().parents[1] / "R4" / "runs"  # R4's run arrays, rebuilt by R4/run_monitors.py
OUT = Path(__file__).resolve().parent
END = "2026-10-25"
SYM = {"BTC": "BTCUSDT", "ETH": "ETHUSDT", "SOL": "SOLUSDT", "XRP": "XRPUSDT", "DOGE": "DOGEUSDT"}
MIN_TRAIN = 720
B = 1999
EVAL0 = pd.Timestamp("2026-01-01", tz="UTC")
EVAL1 = pd.Timestamp(END, tz="UTC")


def winsum(r, a, b, min_frac):
    """Sum of r^2 over indices a..b inclusive; missing returns allowed up to 1 - min_frac, sum rescaled."""
    n = r.size
    miss = ~np.isfinite(r)
    C = np.concatenate([[0.0], np.cumsum(np.where(miss, 0.0, r) ** 2)])
    M = np.concatenate([[0], np.cumsum(miss)])
    a = np.asarray(a); b = np.asarray(b)
    ok = (a >= 0) & (b < n)
    aa, bb = np.clip(a, 0, n - 1), np.clip(b, 0, n - 1)
    L = bb - aa + 1
    nm = M[bb + 1] - M[aa]
    s = (C[bb + 1] - C[aa]) * L / np.maximum(L - nm, 1)
    good = ok & ((L - nm) >= min_frac * L) & (s > 0)
    return np.where(good, s, np.nan)


def grid_crypto(coin):
    from avmonitor.data import load_binance
    df = load_binance(SYM[coin], "2024-01-01", END, download=False)
    r = df["logret"].to_numpy(float)
    n = r.size
    g = np.arange(0, n, 60)
    times = df.index[g]
    out = pd.DataFrame({"asset": coin, "t": g, "time": times})
    out["dayn"] = g // 1440
    DAYL, mf = 1440, 1.0  # E20 rule: a window with any missing return is dropped
    out["RV60"] = winsum(r, g + 1, g + 60, mf)
    out["RVday"] = winsum(r, g + 1, g + DAYL, mf)
    rvd = lambda tt: winsum(r, tt - DAYL + 1, tt, mf)
    out["RVd"] = rvd(g)
    out["RVw"] = np.mean([rvd(g - DAYL * j) for j in range(5)], axis=0)
    out["RVm"] = np.mean([rvd(g - DAYL * j) for j in range(22)], axis=0)
    out["RVh"] = np.maximum(winsum(r, g - 59, g, mf), 0.01 * out["RVd"] / 24)
    out["slot"] = times.hour
    out["dow"] = times.dayofweek
    known = {60: (g + 60) // DAYL + 1, "day": (g + DAYL) // DAYL + 1}
    # plug-in alarms of R4, by timestamp
    z = np.load(SCR / f"{coin}_plugin.npz")
    ai = np.searchsorted(df.index.asi8, z["ts"][z["ev_idx"]])
    assert np.all(df.index.asi8[ai] == z["ts"][z["ev_idx"]])
    tg = -(-ai // 60) * 60
    pos = tg // 60
    for lab, m in (("U", z["ev_up"]), ("D", ~z["ev_up"])):
        v = np.zeros(len(out))
        p = pos[m]; p = p[p < len(out)]
        np.add.at(v, p, 1.0)
        out[lab] = np.minimum(v, 1.0)
    out["gday"] = (out.time.dt.floor("D") - pd.Timestamp("2024-01-01", tz="UTC")).dt.days
    return out, known, 23, 6


def grid_nifty():
    import e07_nifty as E7
    df, _, _ = E7.load_clean("2017-01-01", "2026-11-01")
    df = df[df.index < pd.Timestamp(END, tz="UTC")]
    r = df["logret"].to_numpy(float)
    pos_in = df["minute"].to_numpy() - 555
    sess = pd.factorize(df["session"])[0]
    slots = [59, 119, 179, 239, 299, 359, 374]
    g = np.flatnonzero(np.isin(pos_in, slots))
    out = pd.DataFrame({"asset": "NIFTY", "t": g, "time": df.index[g]})
    out["dayn"] = sess[g]
    DAYL, mf = 375, 0.9  # Nifty: up to 10% missing (frozen stretches) allowed, sum rescaled (declared adaptation)
    out["RV60"] = winsum(r, g + 1, g + 60, mf)
    out["RVday"] = winsum(r, g + 1, g + DAYL, mf)
    rvd = lambda tt: winsum(r, tt - DAYL + 1, tt, mf)
    out["RVd"] = rvd(g)
    out["RVw"] = np.mean([rvd(g - DAYL * j) for j in range(5)], axis=0)
    out["RVm"] = np.mean([rvd(g - DAYL * j) for j in range(22)], axis=0)
    out["RVh"] = np.maximum(winsum(r, g - 59, g, mf), 0.01 * out["RVd"] / 6.25)
    out["slot"] = np.searchsorted(slots, pos_in[g])
    out["dow"] = pd.DatetimeIndex(df["session"].to_numpy()[g]).dayofweek
    n = r.size
    known = {60: sess[np.minimum(g + 60, n - 1)] + 1, "day": sess[np.minimum(g + DAYL, n - 1)] + 1}
    z = np.load(SCR / "NIFTY_plugin.npz")
    assert z["u"].size == n and np.all(z["ts"] == df.index.asi8)
    ev, up = z["ev_idx"], z["ev_up"]
    gi = np.searchsorted(g, ev, side="left")  # first grid row at or after the alarm: alarm in (previous grid, grid]
    for lab, m in (("U", up), ("D", ~up)):
        v = np.zeros(len(out))
        p = gi[m]; p = p[p < len(out)]
        np.add.at(v, p, 1.0)
        out[lab] = np.minimum(v, 1.0)
    out["gday"] = sess[g]  # consecutive session number (blocks of 7 sessions)
    return out, known, 6, 4


def design(df, model, nslot, ndow):
    cols = [np.ones(len(df))] + [np.log(df[c].to_numpy()) for c in ("RVd", "RVw", "RVm", "RVh")]
    if model == "HAR-I+A":
        cols += [df["U"].to_numpy(float), df["D"].to_numpy(float)]
    for h in range(1, nslot + 1):
        cols.append((df["slot"].to_numpy() == h).astype(float))
    for d in range(1, ndow + 1):
        cols.append((df["dow"].to_numpy() == d).astype(float))
    return np.column_stack(cols)


def fit(df, known, nslot, ndow):
    """E20's daily refitted expanding-window OLS (past data only), for HAR-I and HAR-I plus alarm dummies."""
    ndays = int(df["dayn"].max()) + 2
    dn = df["dayn"].to_numpy()
    for H, col in ((60, "RV60"), ("day", "RVday")):
        y = np.log(df[col].to_numpy())
        for model in ("HAR-I", "HAR-I+A"):
            X = design(df, model, nslot, ndow)
            ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
            k = X.shape[1]
            XtX = np.zeros((ndays + 1, k, k)); Xty = np.zeros((ndays + 1, k)); yy = np.zeros(ndays + 1); cnt = np.zeros(ndays + 1)
            kd = np.minimum(known[H][ok], ndays)
            np.add.at(XtX, kd, X[ok][:, :, None] * X[ok][:, None, :]); np.add.at(Xty, kd, X[ok] * y[ok][:, None])
            np.add.at(yy, kd, y[ok] ** 2); np.add.at(cnt, kd, 1)
            XtX, Xty, yy, cnt = (np.cumsum(a, axis=0) for a in (XtX, Xty, yy, cnt))
            beta = np.full((ndays + 1, k), np.nan); s2 = np.full(ndays + 1, np.nan)
            for d in range(ndays + 1):
                if cnt[d] >= MIN_TRAIN:
                    try:
                        b = np.linalg.solve(XtX[d], Xty[d])
                    except np.linalg.LinAlgError:
                        continue
                    beta[d] = b; s2[d] = (yy[d] - 2 * b @ Xty[d] + b @ XtX[d] @ b) / (cnt[d] - k)
            logF = np.sum(X * beta[dn], axis=1)
            df[f"logF_{model}_{H}"] = logF
            df[f"F_{model}_{H}"] = np.exp(logF + s2[dn] / 2)
            if model == "HAR-I+A":
                last = beta[dn][-1]
                df.attrs[f"alarm_coef_last_{H}"] = (float(last[5]), float(last[6]))
    return df


def boot_reg(y, X, day, W, names, jcols):
    """Day-block bootstrap of OLS coefficients: per-day cross products weighted by the bootstrap day counts."""
    D = W.shape[1]
    XtXd = np.zeros((D, X.shape[1], X.shape[1])); Xtyd = np.zeros((D, X.shape[1]))
    np.add.at(XtXd, day, X[:, :, None] * X[:, None, :]); np.add.at(Xtyd, day, X * y[:, None])
    A = np.einsum("bd,dij->bij", W, XtXd); c = np.einsum("bd,di->bi", W, Xtyd)
    bs = np.linalg.solve(A, c[:, :, None])[:, :, 0]
    return {nm: np.percentile(bs[:, j], [2.5, 97.5]) for nm, j in zip(names, jcols)}


def evaluate(G, label, seed):
    rows = []
    for H, col in ((60, "RV60"), ("day", "RVday")):
        m = (G.time >= EVAL0) & (G.time < EVAL1) & np.isfinite(np.log(G[col]))
        for mod in ("HAR-I", "HAR-I+A"):
            m &= np.isfinite(G[f"logF_{mod}_{H}"]) & np.isfinite(G[f"F_{mod}_{H}"])
        # outcome window must end by 24 October 2026 23:59 (crypto: time + H minutes; Nifty: rows exist by construction)
        hmin = 60 if H == 60 else (1440 if label != "NIFTY" else 0)
        if label != "NIFTY":
            m &= (G.time + pd.Timedelta(minutes=hmin)) < EVAL1
        s = G[m]
        assets = sorted(s.asset.unique())
        day = pd.factorize(s.gday)[0]
        order = np.argsort(s.gday.to_numpy()); udays = np.unique(s.gday.to_numpy())
        dayc = np.searchsorted(udays, s.gday.to_numpy())
        y = np.log(s[col].to_numpy())
        A = np.column_stack([(s.asset == a).to_numpy(float) for a in assets])
        lf = s[f"logF_HAR-I_{H}"].to_numpy()
        X = np.hstack([A, lf[:, None], s[["U", "D"]].to_numpy(float)])
        beta, se, nd = ols_dk(y, X, s.gday.to_numpy(), 5)
        Wt = cbb_weights(len(udays), 7, B, seed + (0 if H == 60 else 1))
        k0 = len(assets)
        ci = boot_reg(y, X, dayc, Wt, ["b", "cU", "cD"], [k0, k0 + 1, k0 + 2])
        q0 = qlike(s[col].to_numpy(), s[f"F_HAR-I_{H}"].to_numpy()); q1 = qlike(s[col].to_numpy(), s[f"F_HAR-I+A_{H}"].to_numpy())
        dq = q1 - q0
        Sd = np.bincount(dayc, weights=dq, minlength=len(udays)); Nd = np.bincount(dayc, minlength=len(udays)).astype(float)
        rep = (Wt @ Sd) / (Wt @ Nd)
        # QLIKE only at grid points with an alarm (where the two forecasts can differ most)
        al = (s.U.to_numpy() + s.D.to_numpy()) > 0
        rows.append(dict(asset=label, horizon=H, n=len(s), days=len(udays), up_alarm_points=int(s.U.sum()),
                         down_alarm_points=int(s.D.sum()),
                         b=beta[k0], cU=beta[k0 + 1], cU_lo=ci["cU"][0], cU_hi=ci["cU"][1], cU_se_dk=se[k0 + 1],
                         cD=beta[k0 + 2], cD_lo=ci["cD"][0], cD_hi=ci["cD"][1], cD_se_dk=se[k0 + 2],
                         qlike_HARI=q0.mean(), qlike_HARI_A=q1.mean(), dqlike=dq.mean(),
                         dqlike_lo=np.percentile(rep, 2.5), dqlike_hi=np.percentile(rep, 97.5),
                         dqlike_rel_pct=100 * dq.mean() / q0.mean(),
                         dqlike_at_alarm_points=dq[al].mean() if al.any() else np.nan,
                         seed=seed + (0 if H == 60 else 1)))
    return rows


if __name__ == "__main__":
    t0 = time.time()
    which = sys.argv[1:] or ["BTC", "ETH", "NIFTY", "SOL", "XRP", "DOGE"]
    seeds = {"BTC": 9_989_000, "ETH": 9_989_010, "NIFTY": 9_989_020, "SOL": 9_989_030, "XRP": 9_989_040,
             "DOGE": 9_989_050, "CRYPTO_POOLED": 9_989_060, "BTC+ETH": 9_989_070}
    rows, grids = [], {}
    for a in which:
        if a == "NIFTY":
            G, known, ns, nd = grid_nifty()
        else:
            G, known, ns, nd = grid_crypto(a)
        G = fit(G, known, ns, nd)
        grids[a] = G
        rows += evaluate(G, a, seeds[a])
        print(a, f"{time.time() - t0:.0f} s", G.attrs, flush=True)
        pd.DataFrame(rows).to_csv(OUT / "p2_forecast_2026.csv", index=False)
    cr = [a for a in which if a != "NIFTY"]
    if {"BTC", "ETH"} <= set(cr):
        rows += evaluate(pd.concat([grids["BTC"], grids["ETH"]], ignore_index=True), "BTC+ETH", seeds["BTC+ETH"])
    if len(cr) == 5:
        rows += evaluate(pd.concat([grids[a] for a in cr], ignore_index=True), "CRYPTO_POOLED", seeds["CRYPTO_POOLED"])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "p2_forecast_2026.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.round(4).to_string())
