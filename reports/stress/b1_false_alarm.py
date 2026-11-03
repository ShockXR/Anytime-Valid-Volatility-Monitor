"""B1. Fresh-seed Monte Carlo of the false-alarm guarantee (SYNTHETIC). Not a rerun of experiments/.
Usage: python b1_false_alarm.py <n_paths> <outcsv> <law> [<law>...]
Seeds: law i uses default_rng(9_100_000 + 1000*i) for the monitored paths and 9_150_000 + 1000*i for the 'truth' sample."""
import sys, time
import numpy as np, pandas as pd
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/stress")
from common import wilson, draw, LAWS
from avmonitor import detectors as D, robust as R
from avmonitor.monitor import VolatilityMonitor, GaussianLRMonitor

A = 7800.0           # 390/0.05
CAL = 1000
HORIZ = [98, 390, 1950]
NMON = max(HORIZ)

def truth(law, rng, big=400_000):
    """tau = 3 x scaled MAD of the marginal and m0 = E min(r^2, tau^2): marginal from a long independent sample."""
    if law.startswith("garch"):
        r = draw(law, rng, big // 4, 1)  # one draw after the burn-in per chain: stationary marginal
        r = np.concatenate([r.ravel(), draw(law, rng, big // 4, 1).ravel(), draw(law, rng, big // 4, 1).ravel(), draw(law, rng, big // 4, 1).ravel()])
    else:
        r = draw(law, rng, 1, 5_000_000).ravel()
    scale = np.median(np.abs(r)) * R.MAD_NORMAL_SCALE
    tau = 3 * scale
    return tau, float(np.minimum(r**2, tau**2).mean())

def calibrate_page(seed, drift=0.25, reps=300, n=40000):
    rng = np.random.default_rng(seed)
    r = rng.standard_normal((reps, n))
    lo, hi = 1.0, 60.0
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        al = D.page_cusum(r**2, 1.0, drift, mid).alarm
        arl = np.where(al >= 0, al + 1, n).mean()
        lo, hi = (mid, hi) if arl < A else (lo, mid)
    return hi

def first_alarm_stream(make, data):
    out = np.full(data.shape[0], -1)
    nev = 0
    for i in range(data.shape[0]):
        m = make()
        for r in data[i]:
            ev = m.update(r)
            if ev is not None:
                out[i] = ev.index - ev.monitored_since
                break
    return out

def rows_for(law, method, bound_fn, alarm, N, note):
    rows = []
    for m in HORIZ:
        k = int(((alarm >= 0) & (alarm < m)).sum())
        lo, hi = wilson(k, N)
        b = bound_fn(m)
        verdict = "n/a (no bound)" if b is None else ("within bound" if hi <= b else ("EXCEEDS bound" if lo > b else "inconclusive"))
        rows.append(dict(law=law, method=method, horizon_m=m, k=k, paths=N, rate=k / N, ci_low=lo, ci_high=hi,
                         bound=b, rate_over_bound=None if b is None else (k / N) / b, verdict=verdict, note=note))
    return rows

def run_law(args):
    li, law, N, page_h = args
    t0 = time.time()
    rng = np.random.default_rng(9_100_000 + 1000 * li)
    data = draw(law, rng, N, CAL + NMON)
    cal, mon = data[:, :CAL], data[:, CAL:]
    rows = []
    # (1) known baseline, the setting of the theorem (E[min(r^2,tau^2)] known from the marginal)
    tau, m0 = truth(law, np.random.default_rng(9_150_000 + 1000 * li))
    x = np.minimum(mon**2, tau**2)
    up = D.betting_edetector(x, m0, "up", threshold=2 * A)
    dn = D.betting_edetector(x, m0, "down", upper=tau**2, threshold=2 * A)
    rows += rows_for(law, "e-detector, KNOWN baseline (theorem)", lambda m: m / A, D.two_sided_alarm(up, dn), N,
                     f"tau={tau:.4f} m0={m0:.5f}; proved bound only if conditional mean = m0")
    # (2) streaming monitors from the package, calibration then monitoring
    mk = {
        "monitor plugin band0.25 (README)": (lambda: VolatilityMonitor(threshold=A, cal_len=CAL, null_band=0.25), lambda m: m / A, "approximate: plug-in baseline"),
        "monitor plugin band0": (lambda: VolatilityMonitor(threshold=A, cal_len=CAL), lambda m: m / A, "approximate: plug-in baseline"),
        "monitor cs baseline delta0.01": (lambda: VolatilityMonitor(threshold=A, cal_len=CAL, baseline="cs", delta=0.01), lambda m: 0.01 + m / A, "proved bound delta+m/A if iid-type null"),
        "Gaussian LR monitor (baseline)": (lambda: GaussianLRMonitor(threshold=A, cal_len=CAL), lambda m: m / A, "valid only for exactly normal returns"),
    }
    for name, (make, bf, note) in mk.items():
        rows += rows_for(law, name, bf, first_alarm_stream(make, data), N, note)
    # (3) vectorised classical baselines, scale from the calibration window
    scale = (np.median(np.abs(cal), axis=1) * R.MAD_NORMAL_SCALE)[:, None]
    a_page = D.page_cusum(mon**2 / scale**2, 1.0, 0.25, page_h).alarm
    rows += rows_for(law, "Page CUSUM up (Gaussian-calibrated)", lambda m: None, a_page, N, "no finite-sample guarantee; threshold set for ARL=A on Gaussian")
    a_z = D.repeated_window_ztest(mon / scale, 1.0, 2.0, 60).alarm
    rows += rows_for(law, "repeated 60-step z-test", lambda m: None, a_z, N, "not anytime-valid by design")
    print(law, "done in", round(time.time() - t0), "s", flush=True)
    return rows

if __name__ == "__main__":
    N = int(sys.argv[1]); out = sys.argv[2]; laws = sys.argv[3:] or LAWS
    page_h = calibrate_page(9_100_099)
    print("page_h", page_h, flush=True)
    from multiprocessing import Pool
    jobs = [(LAWS.index(l), l, N, page_h) for l in laws]
    with Pool(2) as p:
        res = p.map(run_law, jobs, chunksize=1)
    df = pd.DataFrame([r for rr in res for r in rr])
    df.to_csv(out, index=False)
    print("EXIT", flush=True)
