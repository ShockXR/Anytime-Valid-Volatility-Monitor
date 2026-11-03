"""B2. Power and detection delay under a variance change (SYNTHETIC). Seeds: 9_200_000 + 100*cell.
Usage: python b2_power_delay.py <paths> <out.csv>
Change factor f multiplies the VARIANCE (sd ratio sqrt(f)) at monitored step nu (0-based, counted after calibration).
Detection = first alarm in [nu, nu+H), H=1950 steps (5 days). Early false alarm = alarm before nu. Delay = alarm - nu + 1."""
import sys, time, itertools
import numpy as np, pandas as pd
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/stress")
from common import wilson
from avmonitor import detectors as D, robust as R
from avmonitor.simulate import change_returns
from avmonitor.monitor import VolatilityMonitor, GaussianLRMonitor
from scipy import stats

A, CAL, H = 7800.0, 1000, 1950
LAWS = {"gaussian": dict(), "jumps_p0.5pct_sd8": dict(p_jump=0.005, jump_scale=8.0)}
FACTORS = [1.5, 2.0, 3.0, 1 / 1.5, 0.5, 1 / 3]
NUS = [50, 400, 1500]

def page_h_cal(seed, drift=0.25, reps=300, n=40000):
    rng = np.random.default_rng(seed); r = rng.standard_normal((reps, n)); lo, hi = 1.0, 60.0
    for _ in range(18):
        mid = 0.5 * (lo + hi); al = D.page_cusum(r**2, 1.0, drift, mid).alarm
        lo, hi = (mid, hi) if np.where(al >= 0, al + 1, n).mean() < A else (lo, mid)
    return hi

def z_level_cal(seed, target=0.05, window=60, reps=3000, n=390):
    """repeated window z-test level such that P(false alarm within 390 steps) = 5 percent on Gaussian data."""
    rng = np.random.default_rng(seed); r = rng.standard_normal((reps, n)); lo, hi = 1e-9, 0.05
    for _ in range(30):
        mid = np.sqrt(lo * hi); fa = (D.repeated_window_ztest(r, 1.0, 2.0, window, level=mid).alarm >= 0).mean()
        lo, hi = (lo, mid) if fa > target else (mid, hi)
    return np.sqrt(lo * hi)

def first_alarm_stream(make, data):
    out = np.full(data.shape[0], -1)
    for i in range(data.shape[0]):
        m = make()
        for r in data[i]:
            ev = m.update(r)
            if ev is not None:
                out[i] = ev.index - ev.monitored_since; break
    return out

def summ(law, direction, f, nu, method, alarm, valid, N):
    early = (alarm >= 0) & (alarm < nu); det = (alarm >= nu) & (alarm < nu + H)
    delays = alarm[det] - nu + 1
    lo, hi = wilson(int(det.sum()), N)
    return dict(law=law, direction=direction, var_factor=round(f, 4), nu=nu, method=method, paths=N, detect_rate=det.mean(), det_ci_low=lo, det_ci_high=hi,
                early_false_alarm=early.mean(), median_delay=float(np.median(delays)) if delays.size else np.nan,
                q25_delay=float(np.quantile(delays, .25)) if delays.size else np.nan, q75_delay=float(np.quantile(delays, .75)) if delays.size else np.nan,
                valid_level_controlled=valid)

def job(args):
    ci, law, f, nu, N, page_h, z_level = args
    rng = np.random.default_rng(9_200_000 + 100 * ci)
    n = CAL + nu + H
    p = change_returns(N, n, nu=CAL + nu, ratio=float(np.sqrt(f)), rng=rng, **LAWS[law])
    data = p.r; cal, mon = data[:, :CAL], data[:, CAL:]
    direction = "increase" if f > 1 else "decrease"
    rows = []
    for name, make, valid in [
        ("e-monitor plugin band0.25 (README)", lambda: VolatilityMonitor(threshold=A, cal_len=CAL, null_band=0.25), "approximate (plug-in)"),
        ("e-monitor plugin band0", lambda: VolatilityMonitor(threshold=A, cal_len=CAL), "approximate (plug-in)"),
        ("e-monitor cs baseline", lambda: VolatilityMonitor(threshold=A, cal_len=CAL, baseline="cs", delta=0.01), "proved (delta + m/A)"),
        ("Gaussian LR monitor", lambda: GaussianLRMonitor(threshold=A, cal_len=CAL), "only for exactly normal returns")]:
        rows.append(summ(law, direction, f, nu, name, first_alarm_stream(make, data), valid, N))
    scale = (np.median(np.abs(cal), axis=1) * R.MAD_NORMAL_SCALE)[:, None]
    if f > 1:
        rows.append(summ(law, direction, f, nu, "Page CUSUM up (ARL=A on Gaussian)", D.page_cusum(mon**2 / scale**2, 1.0, 0.25, page_h).alarm, "no finite-sample guarantee", N))
    rows.append(summ(law, direction, f, nu, "repeated 60-step z, level tuned to 5%/390 on Gaussian", D.repeated_window_ztest(mon / scale, 1.0, 2.0, 60, level=z_level).alarm, "not valid (tuned for Gaussian only)", N))
    rows.append(summ(law, direction, f, nu, "repeated 60-step z, level 0.05 (naive)", D.repeated_window_ztest(mon / scale, 1.0, 2.0, 60, level=0.05).alarm, "not valid (certain false alarm)", N))
    return rows

if __name__ == "__main__":
    N = int(sys.argv[1]); out = sys.argv[2]
    page_h = page_h_cal(9_200_099); z_level = z_level_cal(9_200_098)
    print("page_h", page_h, "z_level", z_level, flush=True)
    jobs = [(i, law, f, nu, N, page_h, z_level) for i, (law, f, nu) in enumerate(itertools.product(LAWS, FACTORS, NUS))]
    from multiprocessing import Pool
    with Pool(2) as pool:
        res = pool.map(job, jobs, chunksize=1)
    pd.DataFrame([r for rr in res for r in rr]).to_csv(out, index=False)
    print("EXIT", flush=True)
