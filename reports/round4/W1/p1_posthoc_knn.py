"""W1 part 1, POST HOC (not in PLAN.md; written after the declared results were seen): nearest-neighbour matching on
recent volatility over several look-backs. For each alarm, the control is drawn at random from the K = 200 pool minutes
closest in (log x60), in (log x60, log x360), or in (log x60, log x360, log x1440), where xW = RV of the last W minutes
over the calibration RV. Nifty: also the same session hour. Same pool, hit rule and test as p1_matched_placebo.py.
Variant 3 (added after variants 0 to 2 were seen): (log x10, log x60, log x360, log x1440), x10 = last 10 minutes,
to test whether the alarm's extra information is about the very last minutes.
Seeds: 9,987,000 + 500 x variant + 50 x asset + 2 x (year - 2017) + direction (9,987,000 to 9,988,819)."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from p1_matched_placebo import SCR, ASSETS, windows

OUT = Path(__file__).resolve().parent
K, NSET = 200, 2000


def pre_w(u, Wn, minfrac=0.8):
    fin = np.isfinite(u); c2 = np.concatenate([[0.0], np.cumsum(np.where(fin, u * u, 0.0))])
    cn = np.concatenate([[0], np.cumsum(fin)]); n = u.size; t = np.arange(n)
    out = np.full(n, np.nan); ok = t >= Wn - 1; tt = t[ok]
    m = cn[tt + 1] - cn[tt + 1 - Wn]; v = (c2[tt + 1] - c2[tt + 1 - Wn]) / np.maximum(m, 1)
    g = (m >= minfrac * Wn) & (v > 0); out[tt[g]] = v[g]
    return out


rows = []
for ai, asset in enumerate(ASSETS):
    z = np.load(SCR / f"{asset}_plugin.npz")
    fin, pre, post, base = windows(z)
    u, calib, year, tod = z["u"], z["calib"], z["year"], z["tod"]
    long_w = (360, 1440) if asset != "NIFTY" else (180, 375)
    p6, p24 = pre_w(u, long_w[0]), pre_w(u, long_w[1])
    ev, evup = z["ev_idx"], z["ev_up"]
    isal = np.zeros(year.size, bool); isal[ev] = True
    valid = ~calib & fin & np.isfinite(pre) & np.isfinite(post) & np.isfinite(base) & (base > 0) & np.isfinite(p6) & np.isfinite(p24)
    bb = np.where(base > 0, base, 1.0)
    p10 = pre_w(u, 10, 0.8)
    valid &= np.isfinite(p10)
    F = np.column_stack([np.log(pre / bb), np.log(p6 / bb), np.log(p24 / bb), np.log(p10 / bb)])
    hb = np.clip((tod - 555) // 60, 0, 6) if asset == "NIFTY" else np.zeros(year.size, int)
    for y in np.unique(year):
        for d, dname in ((0, "up"), (1, "down")):
            sel = ev[(evup == (d == 0)) & (year[ev] == y)]; sel = sel[valid[sel]]
            if sel.size < 10:
                continue
            hit = (lambda a: post[a] > base[a]) if d == 0 else (lambda a: post[a] < base[a])
            for var, ncol in enumerate((1, 2, 3, 4)):
                nbrs = np.empty((sel.size, K), int)
                for h in np.unique(hb[sel]):
                    pool = np.flatnonzero(valid & (year == y) & ~isal & (hb == h))
                    tree = cKDTree(F[pool, :ncol])
                    m = hb[sel] == h
                    _, ix = tree.query(F[sel[m], :ncol], k=K)
                    nbrs[m] = pool[ix.reshape(m.sum(), K)]
                seed = 9_987_000 + 500 * var + 50 * ai + 2 * (int(y) - 2017) + d
                rng = np.random.default_rng(seed)
                Hs = np.array([hit(nbrs[np.arange(sel.size), rng.integers(0, K, sel.size)]).mean() for _ in range(NSET)])
                ha = hit(sel).mean()
                rows.append(dict(asset=asset, year=int(y), direction=dname,
                                 match=["last hour", "last hour + 6h", "last hour + 6h + day", "+ last 10 min"][var], n_alarms=sel.size,
                                 hit_alarm=ha, hit_control=Hs.mean(), hit_diff=ha - Hs.mean(),
                                 p_hit=(1 + np.sum(Hs >= ha)) / (NSET + 1), seed=seed))
    print(asset, flush=True)
df = pd.DataFrame(rows)
df.to_csv(OUT / "p1_posthoc_knn.csv", index=False)
s = []
for (m, dr), g in df.groupby(["match", "direction"]):
    w = g.n_alarms
    s.append(dict(match=m, direction=dr, cells=len(g), cells_p_lt_005=int((g.p_hit < 0.05).sum()),
                  hit_alarm=np.average(g.hit_alarm, weights=w), hit_control=np.average(g.hit_control, weights=w),
                  diff_points=100 * (np.average(g.hit_alarm, weights=w) - np.average(g.hit_control, weights=w)),
                  min_diff=100 * g.hit_diff.min(), max_diff=100 * g.hit_diff.max()))
s = pd.DataFrame(s); s.to_csv(OUT / "p1_posthoc_knn_summary.csv", index=False)
pd.set_option("display.width", 250); print(s.round(3).to_string())
