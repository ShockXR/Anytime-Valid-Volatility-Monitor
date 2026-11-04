"""R4: per asset and calendar year metrics from the run files of run_monitors.py (definitions in PLAN.md)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
ASSETS = ["NIFTY", "BTC", "ETH", "XRP", "SOL", "DOGE", "BITSTAMP"]
MONS = ["plugin", "cs"]
K = 2000  # placebo sets
W = 60  # window length (steps)
MIN_FIN = 50
SEED0 = 9_800_000


def seed_for(ai, mi, year, d):
    return SEED0 + 1000 * ai + 100 * mi + 10 * (year - 2017) + d  # stays in 9,800,000 to 9,806,191


def runs_of(mask):
    """Start and length of each maximal run of True."""
    m = np.concatenate([[False], mask, [False]]).astype(np.int8)
    d = np.diff(m)
    s = np.flatnonzero(d == 1)
    e = np.flatnonzero(d == -1)
    return s, e - s


def to_ts(v):
    v = int(v)
    unit = "ns" if v > 1e17 else "us" if v > 1e14 else "ms" if v > 1e11 else "s"
    return pd.Timestamp(v, unit=unit, tz="UTC")


def load(asset, mon):
    z = np.load(RUNS / f"{asset}_{mon}.npz")
    return {k: z[k] for k in z.files}


def derived(z, day):
    u = z["u"]
    n = u.size
    fin = np.isfinite(u)
    u2 = np.where(fin, u * u, 0.0)
    c2 = np.concatenate([[0.0], np.cumsum(u2)])
    cn = np.concatenate([[0], np.cumsum(fin)])
    t = np.arange(n)
    pre = np.full(n, np.nan)
    post = np.full(n, np.nan)
    ok = (t - W + 1 >= 0) & (t + W < n)
    tt = t[ok]
    npre = cn[tt + 1] - cn[tt + 1 - W]
    npost = cn[tt + 1 + W] - cn[tt + 1]
    vpre = (c2[tt + 1] - c2[tt + 1 - W]) / np.maximum(npre, 1)
    vpost = (c2[tt + 1 + W] - c2[tt + 1]) / np.maximum(npost, 1)
    good = (npre >= MIN_FIN) & (npost >= MIN_FIN) & (vpre > 0) & (vpost > 0)
    pre[tt[good]] = vpre[good]
    post[tt[good]] = vpost[good]
    # segment of each monitored step, baseline RV of each segment's calibration window
    ss, cl = z["seg_start"], z["seg_cal_len"]
    a = np.maximum(ss - cl, 0)
    base = (c2[ss] - c2[a]) / np.maximum(cn[ss] - cn[a], 1)
    seg_id = np.searchsorted(ss, t, side="right") - 1
    calib = z["calib"]
    seg_id[calib] = -1
    up_off = np.zeros(n, bool)
    off = z["seg_tau2"] <= z["seg_m_up"]
    mon_steps = ~calib & (seg_id >= 0)
    up_off[mon_steps] = off[seg_id[mon_steps]]
    base_t = np.full(n, np.nan)
    base_t[mon_steps] = base[seg_id[mon_steps]]
    return dict(fin=fin, pre=pre, post=post, base_t=base_t, up_off=up_off, watched=~calib & fin, seg_base=base)


def placebo_test(rng, lr_al, hit_al, cand_lr, cand_hit, groups_al=None, groups_c=None, up=True):
    """Matched Monte Carlo placebo: one random watched minute per alarm, K sets."""
    n = lr_al.size
    if n == 0 or cand_lr.size == 0:
        return dict(pl_ratio=np.nan, p_ratio=np.nan, pl_hit=np.nan, p_hit=np.nan)
    if groups_al is None:
        draw = rng.integers(0, cand_lr.size, size=(K, n))
    else:
        draw = np.empty((K, n), np.int64)
        order = np.argsort(groups_c, kind="stable")
        gs = groups_c[order]
        for j, g in enumerate(groups_al):
            lo, hi = np.searchsorted(gs, g, "left"), np.searchsorted(gs, g, "right")
            if hi == lo:  # no candidate at this minute: fall back to all minutes
                draw[:, j] = rng.integers(0, cand_lr.size, size=K)
            else:
                draw[:, j] = order[rng.integers(lo, hi, size=K)]
    med = np.median(cand_lr[draw], axis=1)
    hs = cand_hit[draw].mean(axis=1)
    m_al, h_al = np.median(lr_al), hit_al.mean()
    if up:
        p_r = (1 + (med >= m_al).sum()) / (1 + K)
    else:
        p_r = (1 + (med <= m_al).sum()) / (1 + K)
    p_h = (1 + (hs >= h_al).sum()) / (1 + K)
    return dict(pl_ratio=float(np.exp(np.median(med))), p_ratio=float(p_r), pl_hit=float(hs.mean()), p_hit=float(p_h))


def main():
    rows, alarm_rows, stretch_rows = [], [], []
    for ai, asset in enumerate(ASSETS):
        day = 375 if asset == "NIFTY" else 1440
        for mi, mon in enumerate(MONS):
            f = RUNS / f"{asset}_{mon}.npz"
            if not f.exists():
                print("missing", f)
                continue
            z = load(asset, mon)
            D = derived(z, day)
            year, tod, calib = z["year"], z["tod"], z["calib"]
            ev, up, since = z["ev_idx"], z["ev_up"], z["ev_since"]
            # blind stretches
            s, L = runs_of(calib)
            for a_, l_ in zip(s, L):
                stretch_rows.append(dict(asset=asset, monitor=mon, start=int(a_), year=int(year[a_]), steps=int(l_)))
            sy, sl = year[s], L
            # alarm-level table
            lr_all = np.log(D["post"][ev] / D["pre"][ev])
            seg_k = np.searchsorted(z["seg_start"], since)
            base_al = D["seg_base"][seg_k]
            post_al = D["post"][ev]
            hit_all = np.where(up, post_al > base_al, post_al < base_al)
            for j in range(ev.size):
                alarm_rows.append(dict(asset=asset, monitor=mon, idx=int(ev[j]), ts=str(to_ts(z["ts"][ev[j]])),
                                       year=int(year[ev[j]]), direction="up" if up[j] else "down",
                                       steps_since_cal=int(ev[j] - since[j] + 1), log_post_pre=float(lr_all[j]),
                                       post_rv=float(post_al[j]), pre_rv=float(D["pre"][ev[j]]), base_rv=float(base_al[j]),
                                       hit=bool(hit_all[j]) if np.isfinite(post_al[j]) else None, tod=int(tod[ev[j]])))
            for y in np.unique(year):
                my = year == y
                steps = int(my.sum())
                watched = int((D["watched"] & my).sum())
                ey = year[ev] == y
                row = dict(asset=asset, monitor=mon, year=int(y), steps=steps, days=round(steps / day, 1), watched=watched,
                           blind_share=float(calib[my].mean()),
                           upoff_share_of_steps=float(D["up_off"][my].mean()),
                           n_up=int((ey & up).sum()), n_down=int((ey & ~up).sum()))
                row["blind_for_up_share"] = row["blind_share"] + row["upoff_share_of_steps"]
                row["up_per_1000"] = 1000 * row["n_up"] / max(watched, 1)
                row["down_per_1000"] = 1000 * row["n_down"] / max(watched, 1)
                row["alarms_per_1000"] = row["up_per_1000"] + row["down_per_1000"]
                ly = sl[sy == y]
                row.update(n_stretches=int(ly.size), median_stretch=float(np.median(ly)) if ly.size else np.nan,
                           longest_stretch=int(ly.max()) if ly.size else 0)
                row["median_stretch_days"] = row["median_stretch"] / day
                row["longest_stretch_days"] = row["longest_stretch"] / day
                # placebo candidates: watched minutes of this year with valid windows
                cand = np.flatnonzero(D["watched"] & my & np.isfinite(D["pre"]) & np.isfinite(D["post"]) & np.isfinite(D["base_t"]))
                c_lr = np.log(D["post"][cand] / D["pre"][cand])
                for d, is_up in ((1, True), (0, False)):
                    sel = ey & (up == is_up)
                    lr = lr_all[sel]
                    hit = hit_all[sel]
                    v = np.isfinite(lr)
                    lr, hit = lr[v], hit[v]
                    c_hit = D["post"][cand] > D["base_t"][cand] if is_up else D["post"][cand] < D["base_t"][cand]
                    sd = seed_for(ai, mi, int(y), d)
                    rng = np.random.default_rng(sd)
                    if asset == "NIFTY":
                        res = placebo_test(rng, lr, hit, c_lr, c_hit, tod[ev[sel]][v], tod[cand], up=is_up)
                    else:
                        res = placebo_test(rng, lr, hit, c_lr, c_hit, up=is_up)
                    k = "up" if is_up else "down"
                    row[f"{k}_n_valid"] = int(lr.size)
                    row[f"{k}_ratio"] = float(np.exp(np.median(lr))) if lr.size else np.nan
                    row[f"{k}_ratio_placebo"] = res["pl_ratio"]
                    row[f"{k}_p_ratio"] = res["p_ratio"]
                    row[f"{k}_hit"] = float(hit.mean()) if hit.size else np.nan
                    row[f"{k}_hit_placebo"] = res["pl_hit"]
                    row[f"{k}_p_hit"] = res["p_hit"]
                    row[f"{k}_seed"] = sd
                rows.append(row)
            print(asset, mon, "done", flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(HERE / "asset_year.csv", index=False)
    pd.DataFrame(alarm_rows).to_csv(HERE / "alarms.csv", index=False)
    pd.DataFrame(stretch_rows).to_csv(HERE / "blind_stretches.csv", index=False)
    # whole-run summary per asset and monitor
    S = []
    for (a, m), g in T.groupby(["asset", "monitor"], sort=False):
        st = pd.read_csv(HERE / "blind_stretches.csv") if False else None
        S.append(dict(asset=a, monitor=m, years=f"{g.year.min()}-{g.year.max()}", steps=int(g.steps.sum()),
                      blind_share=float((g.blind_share * g.steps).sum() / g.steps.sum()),
                      blind_for_up_share=float((g.blind_for_up_share * g.steps).sum() / g.steps.sum()),
                      up_per_1000=1000 * g.n_up.sum() / g.watched.sum(), down_per_1000=1000 * g.n_down.sum() / g.watched.sum(),
                      n_up=int(g.n_up.sum()), n_down=int(g.n_down.sum())))
    pd.DataFrame(S).to_csv(HERE / "asset_summary.csv", index=False)
    print(pd.DataFrame(S).round(3).to_string())


if __name__ == "__main__":
    main()
