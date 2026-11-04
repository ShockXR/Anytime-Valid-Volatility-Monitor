"""R4: apply the declared pass / fail lines of PLAN.md to asset_year.csv and alarms.csv; write compact tables."""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
T = pd.read_csv(HERE / "asset_year.csv")
AL = pd.read_csv(HERE / "alarms.csv")
CRYPTO = ["BTC", "ETH", "XRP", "SOL", "DOGE", "BITSTAMP"]
BINANCE = ["BTC", "ETH", "XRP", "SOL", "DOGE"]
out = []


def md(df, floatfmt=3):
    return df.round(floatfmt).to_string(index=False)


# ---------------- Q1: stability across years -----------------
q1 = []
for (a, m), g in T.groupby(["asset", "monitor"], sort=False):
    med = g.alarms_per_1000.median()
    stable = ((g.alarms_per_1000 - med).abs() <= 0.3 * med)
    use_up = g.up_p_hit < 0.05
    use_dn = g.down_p_hit < 0.05
    q1.append(dict(asset=a, monitor=m, years=len(g), rate_min=g.alarms_per_1000.min(), rate_max=g.alarms_per_1000.max(),
                   rate_median=med, years_rate_stable=int(stable.sum()),
                   up_hit_min=g.up_hit.min(), up_hit_max=g.up_hit.max(), up_placebo_max=g.up_hit_placebo.max(),
                   years_up_useful=int(use_up.sum()),
                   down_hit_min=g.down_hit.min(), down_hit_max=g.down_hit.max(), down_placebo_max=g.down_hit_placebo.max(),
                   years_down_useful=int(use_dn.sum()),
                   up_ratio_min=g.up_ratio.min(), up_ratio_max=g.up_ratio.max(),
                   years_up_ratio_beats_placebo=int((g.up_p_ratio < 0.05).sum()),
                   down_ratio_min=g.down_ratio.min(), down_ratio_max=g.down_ratio.max(),
                   years_down_ratio_beats_placebo=int((g.down_p_ratio < 0.05).sum())))
Q1 = pd.DataFrame(q1)
Q1.to_csv(HERE / "q1_stability.csv", index=False)
for m in ["plugin", "cs"]:
    g = T[T.monitor == m]
    n = len(g)
    fu, fd = (g.up_p_hit < 0.05).mean(), (g.down_p_hit < 0.05).mean()
    q = Q1[Q1.monitor == m]
    allstable = bool((q.years_rate_stable == q.years).all())
    verdict = "PASS" if fu >= 0.8 and fd >= 0.8 and allstable else ("PARTIAL" if min(fu, fd) >= 0.5 else "FAILED")
    out.append(f"Q1 {m}: asset-years {n}; up useful {fu:.2%}; down useful {fd:.2%}; every asset's rate stable: {allstable} -> {verdict}")
    out.append(f"   ratio metric: up post/pre beats placebo in {(g.up_p_ratio < 0.05).sum()} of {n}; down in {(g.down_p_ratio < 0.05).sum()} of {n}")
    out.append(f"   unstable asset-years: " + ", ".join(f"{r.asset} {r.year} ({r.alarms_per_1000:.2f})" for _, r in g.iterrows()
               if abs(r.alarms_per_1000 - q.set_index('asset').loc[r.asset, 'rate_median']) > 0.3 * q.set_index('asset').loc[r.asset, 'rate_median']))

# ---------------- Q2: 2026 against earlier years -----------------
q2 = []
for (a, m), g in T.groupby(["asset", "monitor"], sort=False):
    g = g.set_index("year")
    if 2026 not in g.index:
        continue
    ref_years = [2025] if a == "BITSTAMP" else ([2024, 2025] if a in BINANCE else list(range(2018, 2026)))
    ref = g.loc[ref_years]
    r26 = g.loc[2026]
    w = ref.steps
    ref_rate = ref.alarms_per_1000.mean()
    d = dict(asset=a, monitor=m, ref_years=f"{min(ref_years)}-{max(ref_years)}",
             rate_ref=ref_rate, rate_2026=r26.alarms_per_1000, rate_change=r26.alarms_per_1000 / ref_rate - 1,
             blind_ref=ref.blind_share.mean(), blind_2026=r26.blind_share,
             up_hit_ref=ref.up_hit.mean(), up_hit_2026=r26.up_hit, down_hit_ref=ref.down_hit.mean(), down_hit_2026=r26.down_hit,
             up_ratio_ref=ref.up_ratio.mean(), up_ratio_2026=r26.up_ratio,
             down_per_1000_ref=ref.down_per_1000.mean(), down_per_1000_2026=r26.down_per_1000)
    d["ok_rate"] = abs(d["rate_change"]) <= 0.3
    d["ok_blind"] = abs(d["blind_2026"] - d["blind_ref"]) <= 0.10
    d["ok_up"] = abs(d["up_hit_2026"] - d["up_hit_ref"]) <= 0.10
    d["ok_down"] = abs(d["down_hit_2026"] - d["down_hit_ref"]) <= 0.10
    d["all_ok"] = d["ok_rate"] and d["ok_blind"] and d["ok_up"] and d["ok_down"]
    q2.append(d)
Q2 = pd.DataFrame(q2)
Q2.to_csv(HERE / "q2_2026.csv", index=False)
for m in ["plugin", "cs"]:
    q = Q2[Q2.monitor == m]
    k = int(q.all_ok.sum())
    out.append(f"Q2 {m}: assets with all four checks within limits: {k} of {len(q)} -> {'PASS' if k >= 5 else 'FAILED'}")

# ---------------- Q3: blind share -----------------
p = T[T.monitor == "plugin"]
c = T[T.monitor == "cs"].set_index(["asset", "year"])
q3 = p[["asset", "year", "blind_share"]].copy()
q3["in_60_80"] = q3.blind_share.between(0.6, 0.8)
q3["cs_blind_share"] = [c.loc[(a, y), "blind_share"] for a, y in zip(q3.asset, q3.year)]
q3["cs_blind_for_up"] = [c.loc[(a, y), "blind_for_up_share"] for a, y in zip(q3.asset, q3.year)]
q3.to_csv(HERE / "q3_blind.csv", index=False)
cr = q3[q3.asset.isin(CRYPTO)]
ni = q3[q3.asset == "NIFTY"]
out.append(f"Q3 crypto asset-years in [60%, 80%]: {int(cr.in_60_80.sum())} of {len(cr)} (range {cr.blind_share.min():.3f} to {cr.blind_share.max():.3f}) -> "
           f"{'VERIFIED' if cr.in_60_80.mean() >= 0.8 else 'FAILED'}")
out.append(f"Q3 Nifty years in [60%, 80%]: {int(ni.in_60_80.sum())} of {len(ni)} (range {ni.blind_share.min():.3f} to {ni.blind_share.max():.3f})")
out.append(f"Q3 cs crypto calibrating share range {cr.cs_blind_share.min():.3f} to {cr.cs_blind_share.max():.3f}; blind for up {cr.cs_blind_for_up.min():.3f} to {cr.cs_blind_for_up.max():.3f}")
out.append(f"Q3 cs Nifty calibrating share range {ni.cs_blind_share.min():.3f} to {ni.cs_blind_share.max():.3f}")

# ---------------- Q4: early against late alarms -----------------
rng = np.random.default_rng(9_800_201)
AL["fam"] = np.where(AL.asset == "NIFTY", "Nifty", np.where(AL.asset == "BITSTAMP", "Bitstamp BTC", "Binance (5 coins)"))
day = np.where(AL.asset == "NIFTY", 375, 1440)
AL["bin"] = np.where(AL.steps_since_cal <= 60, "early (<= 60 min)", np.where(AL.steps_since_cal <= day, "middle (61 min to 1 day)", "late (> 1 day)"))
AL["post_base"] = AL.post_rv / AL.base_rv
q4 = []
tests = []
for (fam, mon, dirn), g in AL.groupby(["fam", "monitor", "direction"], sort=False):
    g = g[np.isfinite(g.log_post_pre)]
    for b, h in g.groupby("bin"):
        q4.append(dict(family=fam, monitor=mon, direction=dirn, bin=b, n=len(h), share_of_alarms=len(h) / len(g),
                       hit=h.hit.astype(float).mean(), ratio_post_pre=float(np.exp(np.median(h.log_post_pre))),
                       ratio_post_base=float(np.median(h.post_base))))
    e = g[g.bin.str.startswith("early")].hit.astype(float).to_numpy()
    l = g[g.bin.str.startswith("late")].hit.astype(float).to_numpy()
    if e.size >= 5 and l.size >= 5:
        obs = e.mean() - l.mean()
        allv = np.concatenate([e, l])
        perm = np.empty(2000)
        for i in range(2000):
            rng.shuffle(allv)
            perm[i] = allv[: e.size].mean() - allv[e.size:].mean()
        pval = (1 + (perm <= obs).sum()) / 2001
        el = np.exp(np.median(g[g.bin.str.startswith("early")].log_post_pre))
        ll = np.exp(np.median(g[g.bin.str.startswith("late")].log_post_pre))
        tests.append(dict(family=fam, monitor=mon, direction=dirn, n_early=e.size, n_late=l.size, hit_early=e.mean(),
                          hit_late=l.mean(), diff_points=100 * obs, p_one_sided=pval, ratio_early=el, ratio_late=ll,
                          less_useful=bool(obs <= -0.10 and pval < 0.05)))
Q4 = pd.DataFrame(q4)
Q4.to_csv(HERE / "q4_restart_bins.csv", index=False)
Q4t = pd.DataFrame(tests)
Q4t.to_csv(HERE / "q4_restart_tests.csv", index=False)

(HERE / "verdicts.txt").write_text("\n".join(out) + "\n")
print("\n".join(out))
print(md(Q1))
print(md(Q2))
print(md(Q4))
print(md(Q4t))
