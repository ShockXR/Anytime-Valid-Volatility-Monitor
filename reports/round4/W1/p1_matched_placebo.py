"""W1 part 1: volatility-matched placebo for the plug-in (README) monitor, every asset and year. See PLAN.md.
Uses the R4 run arrays (same definitions as R4 and V1: 60-minute windows, >= 50 finite returns, calibration base)."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

SCR = Path(__file__).resolve().parents[1] / "R4" / "runs"  # R4's run arrays, rebuilt by R4/run_monitors.py
OUT = Path(__file__).resolve().parent
ASSETS = ["NIFTY", "BTC", "ETH", "XRP", "SOL", "DOGE", "BITSTAMP"]
W, NSET = 60, 2000


def windows(z):
    u, calib = z["u"], z["calib"]
    n = u.size
    fin = np.isfinite(u); u2 = np.where(fin, u * u, 0.0)
    c2 = np.concatenate([[0.0], np.cumsum(u2)]); cn = np.concatenate([[0], np.cumsum(fin)])
    t = np.arange(n)
    pre = np.full(n, np.nan); post = np.full(n, np.nan)
    ok = (t >= W - 1) & (t + W < n); tt = t[ok]
    npre = cn[tt + 1] - cn[tt + 1 - W]; npost = cn[tt + 1 + W] - cn[tt + 1]
    vp = (c2[tt + 1] - c2[tt + 1 - W]) / np.maximum(npre, 1); vq = (c2[tt + 1 + W] - c2[tt + 1]) / np.maximum(npost, 1)
    g = (npre >= 50) & (npost >= 50) & (vp > 0) & (vq > 0)
    pre[tt[g]] = vp[g]; post[tt[g]] = vq[g]
    ss, cl = z["seg_start"], z["seg_cal_len"]
    base_seg = (c2[ss] - c2[ss - cl]) / np.maximum(cn[ss] - cn[ss - cl], 1)
    sid = np.searchsorted(ss, t, side="right") - 1
    base = np.where(~calib & (sid >= 0), base_seg[np.maximum(sid, 0)], np.nan)
    return fin, pre, post, base


def run(nbins_list=(10, 100)):
    rows = []
    for ai, asset in enumerate(ASSETS):
        z = np.load(SCR / f"{asset}_plugin.npz")
        fin, pre, post, base = windows(z)
        calib, year, tod = z["calib"], z["year"], z["tod"]
        ev, evup = z["ev_idx"], z["ev_up"]
        isal = np.zeros(year.size, bool); isal[ev] = True
        valid = ~calib & fin & np.isfinite(pre) & np.isfinite(post) & np.isfinite(base) & (base > 0)
        x = np.where(valid, pre / np.where(base > 0, base, 1.0), np.nan)
        lr = np.where(valid, np.log(post / np.where(base > 0, base, 1.0)), np.nan)
        hb = np.clip((tod - 555) // 60, 0, 6) if asset == "NIFTY" else np.zeros(year.size, int)
        for y in np.unique(year):
            pool = np.flatnonzero(valid & (year == y) & ~isal)
            if pool.size == 0:
                continue
            for d, dname in ((0, "up"), (1, "down")):
                sel = ev[(evup == (d == 0)) & (year[ev] == y)]
                sel = sel[valid[sel]]
                if sel.size == 0:
                    continue
                hit_fn = (lambda a: post[a] > base[a]) if d == 0 else (lambda a: post[a] < base[a])
                for nb in nbins_list:
                    edges = np.quantile(x[pool], np.linspace(0, 1, nb + 1)[1:-1])
                    bp = np.searchsorted(edges, x[pool], side="right")
                    ba = np.searchsorted(edges, x[sel], side="right")
                    key_p = bp * 7 + hb[pool]; key_a = ba * 7 + hb[sel]
                    order = np.argsort(key_p, kind="stable"); kp = key_p[order]; pool_o = pool[order]
                    lo = np.searchsorted(kp, key_a, side="left"); hi = np.searchsorted(kp, key_a, side="right")
                    have = hi > lo
                    # Nifty: fall back to the decile alone if a stratum with the hour is empty (counted)
                    nfall = int((~have).sum())
                    if nfall:
                        key_p2 = bp; o2 = np.argsort(key_p2, kind="stable"); k2 = key_p2[o2]; pool2 = pool[o2]
                    seed = 9_980_000 + 1000 * ai + 10 * (int(y) - 2017) + d + (2 if nb == 100 else 0)
                    rng = np.random.default_rng(seed)
                    hit_a = hit_fn(sel).mean(); lr_a = lr[sel].mean()
                    H = np.empty(NSET); L = np.empty(NSET)
                    for s in range(NSET):
                        pick = np.empty(sel.size, int)
                        r = rng.random(sel.size)
                        pick[have] = pool_o[lo[have] + (r[have] * (hi[have] - lo[have])).astype(int)]
                        if nfall:
                            l2 = np.searchsorted(k2, ba[~have], "left"); h2 = np.searchsorted(k2, ba[~have], "right")
                            pick[~have] = pool2[l2 + (r[~have] * (h2 - l2)).astype(int)]
                        H[s] = hit_fn(pick).mean(); L[s] = lr[pick].mean()
                    if d == 0:
                        p_hit = (1 + np.sum(H >= hit_a)) / (NSET + 1); p_lr = (1 + np.sum(L >= lr_a)) / (NSET + 1)
                    else:
                        p_hit = (1 + np.sum(H >= hit_a)) / (NSET + 1); p_lr = (1 + np.sum(L <= lr_a)) / (NSET + 1)
                    rnd_hit = hit_fn(pool).mean()
                    rows.append(dict(asset=asset, year=int(y), direction=dname, bins=nb, n_alarms=int(sel.size),
                                     pool=int(pool.size), median_x_alarm=float(np.median(x[sel])),
                                     median_x_control=float(np.median(x[pick])), hit_alarm=hit_a,
                                     hit_control=float(H.mean()), hit_random=float(rnd_hit), hit_diff=hit_a - H.mean(),
                                     p_hit=p_hit, rv_ratio_alarm_over_control=float(np.exp(lr_a - L.mean())),
                                     p_rv=p_lr, nifty_hour_fallbacks=nfall, seed=seed))
            print(asset, y, flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = run()
    df.to_csv(OUT / "p1_matched_placebo.csv", index=False)
    print(df.to_string())
