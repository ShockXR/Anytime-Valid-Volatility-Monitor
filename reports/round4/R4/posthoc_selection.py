"""R4 post hoc [SYNTHETIC]: is post / pre < 1 for up alarms expected even when the change is persistent?
Plug-in monitor (A = 28,800, cal 1,440, band 0.25), i.i.d. returns, sd multiplied by 1.5 (or 0.6) for ever
at step 1,440 + 300. Post / pre as in the real-data metric. Seeds 9,800,301 to 9,800,304."""
import numpy as np
import pandas as pd
from avmonitor.monitor import VolatilityMonitor

rows = []
for k, (law, ratio) in enumerate([("gauss", 1.5), ("t3", 1.5), ("gauss", 0.6), ("t3", 0.6)]):
    rng = np.random.default_rng(9_800_301 + k)
    lr, hit = [], []
    for p in range(300):
        n = 1440 + 3000
        x = rng.standard_normal(n) if law == "gauss" else rng.standard_t(3, n) / np.sqrt(3.0)
        cp = 1440 + 300
        x[cp:] *= ratio
        mon = VolatilityMonitor(threshold=28800.0, c=3.0, cal_len=1440, null_band=0.25)
        for t in range(n):
            ev = mon.update(x[t])
            if ev is not None:
                if t >= cp + 60 and t + 60 < n and ((ev.direction == "up") == (ratio > 1)):
                    pre = np.mean(x[t - 59:t + 1] ** 2)
                    post = np.mean(x[t + 1:t + 61] ** 2)
                    lr.append(np.log(post / pre))
                break
    lr = np.array(lr)
    rows.append(dict(law=law, sd_ratio=ratio, alarms_after_change=lr.size, median_post_pre=float(np.exp(np.median(lr))),
                     share_post_below_pre=float((lr < 0).mean()), seed=9_800_301 + k))
R = pd.DataFrame(rows)
R.to_csv("posthoc_selection.csv", index=False)
print(R.round(3).to_string(index=False))
