"""R2 check: the vectorised statistics in r2_lib.py give the same first alarms as the package's own code.
Seeds 9_690_000 + cell. Usage: python r2_validate.py <paths per cell>"""
import sys, time
import numpy as np
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/round4/R2")
from r2_lib import A, CAL, METHODS, draw, truth, all_stats, first_alarm
from avmonitor import detectors as D
from avmonitor.robust import MAD_NORMAL_SCALE
from avmonitor.monitor import VolatilityMonitor, GaussianLRMonitor
from scipy import stats


def stream_first(make, data):
    out = np.full(data.shape[0], -1)
    for i in range(data.shape[0]):
        m = make()
        for r in data[i]:
            ev = m.update(r)
            if ev is not None:
                out[i] = ev.index - ev.monitored_since
                break
    return out


N = int(sys.argv[1])
cells = [("gaussian", 2.0), ("student_t3", 0.5), ("jumps_p0.5pct_sd8", 1.0), ("gaussian", 1.0)]
tot_bad = 0
for ci, (law, f) in enumerate(cells):
    t0 = time.time()
    tau2_t, m0_t = truth(law, 9_600_000 + ci, big=1_000_000)
    rng = np.random.default_rng(9_690_000 + ci)
    nu = 400 if f != 1.0 else None
    data = draw(law, N, CAL + 1200, rng, nu=None if nu is None else CAL + nu, f=f)
    W = all_stats(data, tau2_t, m0_t)
    cal, mon = data[:, :CAL], data[:, CAL:]
    scale = (np.median(np.abs(cal), axis=1) * MAD_NORMAL_SCALE)[:, None]
    for meth in METHODS:
        cm = np.maximum.accumulate(W[meth], axis=1)
        # two thresholds: the proved/nominal one and a low one hit by about half the paths
        h_low = float(np.median(cm[:, 389]))
        if meth == METHODS[6]:
            h_low = min(h_low, 6.0)  # beyond about 8 the level 2(1 - Phi(h)) underflows in the package's ppf call
        for h in ([2 * A, h_low] if meth in METHODS[:5] else [h_low]):
            mine = first_alarm(cm, h)
            if meth == METHODS[0]:
                ref = stream_first(lambda: VolatilityMonitor(threshold=h / 2, cal_len=CAL, null_band=0.25), data)
            elif meth == METHODS[1]:
                ref = stream_first(lambda: VolatilityMonitor(threshold=h / 2, cal_len=CAL), data)
            elif meth == METHODS[2]:
                ref = stream_first(lambda: VolatilityMonitor(threshold=h / 2, cal_len=CAL, baseline="cs", delta=0.01), data)
            elif meth == METHODS[3]:
                x = np.minimum(mon**2, tau2_t)
                up = D.betting_edetector(x, m0_t, "up", threshold=h)
                dn = D.betting_edetector(x, m0_t, "down", upper=tau2_t, threshold=h)
                ref = D.two_sided_alarm(up, dn)
            elif meth == METHODS[4]:
                ref = stream_first(lambda: GaussianLRMonitor(threshold=h / 2, cal_len=CAL), data)
            elif meth == METHODS[5]:
                ref = D.page_cusum(mon**2 / scale**2, 1.0, 0.25, h).alarm
            else:
                level = 2 * (1 - stats.norm.cdf(h))
                ref = D.repeated_window_ztest(mon / scale, 1.0, 2.0, 60, level=level).alarm
            bad = int((ref != mine).sum())
            tot_bad += bad
            print(f"{law:18s} f={f:<4} {meth:30s} h={h:10.4g} alarms={int((mine >= 0).sum()):3d}/{N} mismatches={bad}", flush=True)
    print(law, f, "cell done in", round(time.time() - t0), "s", flush=True)
print("TOTAL MISMATCHES", tot_bad)
