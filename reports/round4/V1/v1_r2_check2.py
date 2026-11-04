"""V1: two more R2 cells, streaming package classes on Student t3 (own generator, seeds 9_960_010..13).
D: GaussianLRMonitor at R2's t3 oracle 2.834e122 (R2: null FA 4.8%; x2 249 (167-343), 95%).
E: README VolatilityMonitor at R2's t3 oracle (R2: null FA 4.9%; x2 80 (48-130), 95%)."""
import json, time
import numpy as np
from avmonitor.monitor import VolatilityMonitor, GaussianLRMonitor
TH = json.load(open("/home/user/anytime-valid-volatility-monitor/reports/round4/R2/r2_thresholds.json"))
h_glr = dict(TH["student_t3|Gaussian LR monitor"])["oracle"]
h_rd = dict(TH["student_t3|README: plug-in, band 0.25"])["oracle"]
CAL, NU, H = 1000, 400, 1950
t3 = lambda rng, n, L: rng.standard_t(3, (n, L)) / np.sqrt(3.0)
def first(mk, x):
    m = mk()
    for r in x:
        ev = m.update(float(r))
        if ev is not None:
            return ev.index - ev.monitored_since
    return -1
t0 = time.time()
for tag, mk, h, s in (("D GLR t3 oracle", lambda: GaussianLRMonitor(threshold=h_glr / 2, cal_len=CAL), h_glr, 9_960_010),
                      ("E README t3 oracle", lambda: VolatilityMonitor(threshold=h_rd / 2, cal_len=CAL, null_band=0.25), h_rd, 9_960_012)):
    al = np.array([first(mk, x) for x in t3(np.random.default_rng(s), 2000, CAL + 390)])
    p = (al >= 0).mean()
    print(f"{tag} h={h:.4g}: null FA within 390 = {p:.4f} +- {1.96*np.sqrt(p*(1-p)/al.size):.4f}", flush=True)
    x = t3(np.random.default_rng(s + 1), 1000, CAL + NU + H); x[:, CAL + NU:] *= np.sqrt(2.0)
    al = np.array([first(mk, r) for r in x])
    d = al[al >= NU] - NU + 1
    print(f"{tag} x2: early {(al >= 0).mean() - (al >= NU).mean():.3f} detected {(al >= NU).mean():.3f} "
          f"median {np.median(d):.0f} IQR ({np.percentile(d, 25):.0f}-{np.percentile(d, 75):.0f}) time {time.time()-t0:.0f}s", flush=True)
