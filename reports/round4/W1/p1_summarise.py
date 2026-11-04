"""W1 part 1 summary: verdict counts per direction and match (PLAN.md lines), Holm over cells, pooled differences."""
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
d = pd.read_csv(OUT / "p1_matched_placebo.csv")


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = p.size
    adj = np.maximum.accumulate((m - np.arange(m)) * p[o]); out = np.empty(m); out[o] = np.minimum(adj, 1)
    return out


rows = []
for (nb, dr), s in d.groupby(["bins", "direction"]):
    s = s[s.n_alarms >= 10]
    w = s.n_alarms
    pooled_a = np.average(s.hit_alarm, weights=w); pooled_c = np.average(s.hit_control, weights=w)
    pooled_r = np.average(s.hit_random, weights=w)
    share = (s.p_hit < 0.05).mean(); diff = 100 * (pooled_a - pooled_c)
    verdict = "VERIFIED" if (share >= 0.8 and diff >= 5) else ("FAILED" if (share <= 0.2 or diff < 2) else "PARTIAL")
    rows.append(dict(bins=nb, direction=dr, cells=len(s), cells_p_lt_005=int((s.p_hit < 0.05).sum()),
                     cells_holm_lt_005=int((holm(s.p_hit) < 0.05).sum()), cells_rv_p_lt_005=int((s.p_rv < 0.05).sum()),
                     alarms=int(w.sum()), pooled_hit_alarm=pooled_a, pooled_hit_matched=pooled_c, pooled_hit_random=pooled_r,
                     pooled_diff_points=diff, min_diff=100 * s.hit_diff.min(), max_diff=100 * s.hit_diff.max(),
                     rv_ratio_median=s.rv_ratio_alarm_over_control.median(), rv_ratio_min=s.rv_ratio_alarm_over_control.min(),
                     rv_ratio_max=s.rv_ratio_alarm_over_control.max(), verdict=verdict))
v = pd.DataFrame(rows)
v.to_csv(OUT / "p1_verdicts.csv", index=False)
pd.set_option("display.width", 250)
print(v.round(3).to_string())
# per asset (pooled over years), percentile match
pa = []
for (nb, dr, a), s in d.groupby(["bins", "direction", "asset"]):
    w = s.n_alarms
    pa.append(dict(bins=nb, direction=dr, asset=a, years=len(s), alarms=int(w.sum()),
                   hit_alarm=np.average(s.hit_alarm, weights=w), hit_matched=np.average(s.hit_control, weights=w),
                   hit_random=np.average(s.hit_random, weights=w), years_p_lt_005=int((s.p_hit < 0.05).sum()),
                   rv_ratio_range=f"{s.rv_ratio_alarm_over_control.min():.2f} to {s.rv_ratio_alarm_over_control.max():.2f}"))
pa = pd.DataFrame(pa)
pa.to_csv(OUT / "p1_by_asset.csv", index=False)
print(pa.round(3).to_string())
