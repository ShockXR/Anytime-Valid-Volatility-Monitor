"""V1 probe of the R4 hit-share placebo. Same definitions as R4 (60-minute windows, >= 50 finite returns,
calibration-window baseline), recomputed here. Compares the hit share of up alarms with (a) all watched minutes
(R4's placebo pool) and (b) watched minutes whose previous 60 minutes were as volatile, relative to their
baseline, as before the alarms (pre/base at or above the 25th percentile of the alarms' pre/base; Nifty also
by time-of-day hour). Uses the run files in the scratchpad, which v1_r4_check.py showed equal to a fresh run."""
import numpy as np
from pathlib import Path
SCR = Path("/tmp/claude-0/-home-user/413dde72-c15d-5c36-9682-3b5d711ff95d/scratchpad/R4_runs")
W = 60
for asset, y in (("BTC", 2025), ("NIFTY", 2024)):
    z = np.load(SCR / f"{asset}_plugin.npz")
    u, calib, year, tod = z["u"], z["calib"], z["year"], z["tod"]
    n = u.size
    fin = np.isfinite(u); u2 = np.where(fin, u * u, 0.0)
    c2 = np.concatenate([[0.0], np.cumsum(u2)]); cn = np.concatenate([[0], np.cumsum(fin)])
    t = np.arange(n)
    ok = (t >= W - 1) & (t + W < n)
    pre = np.full(n, np.nan); post = np.full(n, np.nan)
    tt = t[ok]
    npre = cn[tt + 1] - cn[tt + 1 - W]; npost = cn[tt + 1 + W] - cn[tt + 1]
    vp = (c2[tt + 1] - c2[tt + 1 - W]) / np.maximum(npre, 1); vq = (c2[tt + 1 + W] - c2[tt + 1]) / np.maximum(npost, 1)
    g = (npre >= 50) & (npost >= 50) & (vp > 0) & (vq > 0)
    pre[tt[g]] = vp[g]; post[tt[g]] = vq[g]
    ss, cl = z["seg_start"], z["seg_cal_len"]
    base_seg = (c2[ss] - c2[ss - cl]) / np.maximum(cn[ss] - cn[ss - cl], 1)
    sid = np.searchsorted(ss, t, side="right") - 1
    base = np.where(~calib & (sid >= 0), base_seg[np.maximum(sid, 0)], np.nan)
    ev, up = z["ev_idx"], z["ev_up"]
    sel = ev[up & (year[ev] == y)]
    sel = sel[np.isfinite(pre[sel]) & np.isfinite(post[sel])]
    hit_al = (post[sel] > base[sel]).mean()
    cand = np.flatnonzero(~calib & fin & (year == y) & np.isfinite(pre) & np.isfinite(post) & np.isfinite(base))
    hit_c = (post[cand] > base[cand]).mean()
    q = np.quantile(pre[sel] / base[sel], 0.25)
    hi = cand[pre[cand] / base[cand] >= q]
    print(f"{asset} {y}: up alarms {sel.size}, hit {hit_al:.3f} | all watched minutes {hit_c:.3f} (R4 placebo pool) | "
          f"watched minutes with pre/base >= {q:.2f} (alarms' 25th pct): {hi.size} minutes, hit {(post[hi] > base[hi]).mean():.3f}; "
          f"median pre/base at alarms {np.median(pre[sel] / base[sel]):.2f}")
    if asset == "NIFTY":
        # time-of-day check: alarms' hours vs pool hours
        print("   alarm tod quartiles", np.percentile(tod[sel], [25, 50, 75]), "pool tod quartiles", np.percentile(tod[cand], [25, 50, 75]))
