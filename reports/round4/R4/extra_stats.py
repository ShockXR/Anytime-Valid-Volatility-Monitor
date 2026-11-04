"""R4: extra descriptive numbers used in REPORT.md (blind-for-increase stretches of cs, pooled post/pre shares)."""
from pathlib import Path
import numpy as np
import pandas as pd
import analyse as AN

HERE = Path(__file__).resolve().parent
rows = []
for asset in AN.ASSETS:
    day = 375 if asset == "NIFTY" else 1440
    for mon in AN.MONS:
        z = AN.load(asset, mon)
        D = AN.derived(z, day)
        blind_up = z["calib"] | D["up_off"]
        s, L = AN.runs_of(blind_up)
        yr = z["year"][s]
        for y in np.unique(z["year"]):
            l = L[yr == y]
            rows.append(dict(asset=asset, monitor=mon, year=int(y), n_blind_for_up_stretches=int(l.size),
                             median_days=float(np.median(l)) / day, longest_days=float(l.max()) / day,
                             share_of_steps_in_stretches_over_2_days=float(l[l > 2 * day].sum() / (z["year"] == y).sum())))
B = pd.DataFrame(rows)
B.to_csv(HERE / "blind_for_up_stretches.csv", index=False)
print(B[B.monitor == "cs"].round(3).to_string(index=False))
print(B[(B.monitor == "plugin")].groupby("asset")[["median_days", "longest_days"]].max())
AL = pd.read_csv(HERE / "alarms.csv")
AL = AL[np.isfinite(AL.log_post_pre)]
AL["fam"] = np.where(AL.asset == "NIFTY", "Nifty", "crypto")
g = AL.groupby(["fam", "monitor", "direction"])
print(pd.DataFrame(dict(n=g.size(), share_post_lt_pre=g.log_post_pre.apply(lambda x: (x < 0).mean()),
                        median_ratio=g.log_post_pre.apply(lambda x: np.exp(np.median(x))),
                        hit=g.hit.apply(lambda x: x.astype(float).mean()),
                        median_post_over_base=g.apply(lambda h: np.median(h.post_rv / h.base_rv), include_groups=False))).round(3))
