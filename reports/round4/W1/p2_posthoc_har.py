"""W1 part 2, POST HOC (not in PLAN.md): the same 2026 regression as p2_forecast_2026.py, but against the plain HAR
forecast (no last-hour term), to see whether the alarm dummies matter once the last hour is left out.
Seeds: 9,989,100 + 10 x case + horizon (0 next hour, 1 next day)."""
import numpy as np
import pandas as pd
import p2_forecast_2026 as P

_orig = P.design


def design_har(df, model, nslot, ndow):
    X = _orig(df, "HAR-I", nslot, ndow)
    return np.delete(X, 4, axis=1)  # drop log RVh (column 4: const, RVd, RVw, RVm, RVh, slot..., dow...)


rows = []
cases = ["BTC", "ETH", "NIFTY"]
for ci, a in enumerate(cases):
    G, known, ns, nd = P.grid_nifty() if a == "NIFTY" else P.grid_crypto(a)
    P.design = lambda df, model, nslot, ndow: design_har(df, model, nslot, ndow)
    G = P.fit(G, known, ns, nd)  # both columns now hold the plain HAR fit; only "HAR-I" ones are used below
    P.design = _orig
    for H, col in ((60, "RV60"), ("day", "RVday")):
        m = (G.time >= P.EVAL0) & (G.time < P.EVAL1) & np.isfinite(np.log(G[col])) & np.isfinite(G[f"logF_HAR-I_{H}"])
        if a != "NIFTY":
            m &= (G.time + pd.Timedelta(minutes=60 if H == 60 else 1440)) < P.EVAL1
        s = G[m]
        udays = np.unique(s.gday.to_numpy()); dayc = np.searchsorted(udays, s.gday.to_numpy())
        y = np.log(s[col].to_numpy())
        X = np.column_stack([np.ones(len(s)), s[f"logF_HAR-I_{H}"].to_numpy(), s.U.to_numpy(float), s.D.to_numpy(float)])
        beta, se, _ = P.ols_dk(y, X, s.gday.to_numpy(), 5)
        seed = 9_989_100 + 10 * ci + (0 if H == 60 else 1)
        W = P.cbb_weights(len(udays), 7, P.B, seed)
        cis = P.boot_reg(y, X, dayc, W, ["cU", "cD"], [2, 3])
        rows.append(dict(asset=a, horizon=H, forecast="HAR (no last hour)", n=len(s), b=beta[1], cU=beta[2],
                         cU_lo=cis["cU"][0], cU_hi=cis["cU"][1], cD=beta[3], cD_lo=cis["cD"][0], cD_hi=cis["cD"][1], seed=seed))
    print(a, flush=True)
df = pd.DataFrame(rows)
df.to_csv(P.OUT / "p2_posthoc_har.csv", index=False)
pd.set_option("display.width", 250); print(df.round(4).to_string())
