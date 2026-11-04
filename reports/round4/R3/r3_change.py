"""R3 realistic change shapes (SYNTHETIC), Gaussian returns. Usage: python r3_change.py <n_paths> <outcsv>
Cell j uses default_rng(9_760_000 + 100 j). Methods: README plug-in band 0.25, plug-in band 0, cs baseline delta 0.01."""
import sys, time
import numpy as np, pandas as pd
from r3_common import *

def sd_profile(kind, n_cal=CAL, n_mon=NMON):
    """Standard deviation at each step of calibration + monitoring for one change shape."""
    v = np.ones(n_cal + n_mon)
    if kind.startswith("ramp"):            # variance 1 -> 2 linearly over 390 steps, starting at monitored step 50
        s = n_cal + 50
        v[s:s + 390] = np.linspace(1, 2, 390); v[s + 390:] = 2.0
    elif kind.startswith("burst"):          # x3 variance for L steps starting at monitored step 200
        L = int(kind.split("_")[1]); s = n_cal + 200
        v[s:s + L] = 3.0
    elif kind.startswith("incal"):          # x2 variance from calibration step c onwards, forever
        c = int(kind.split("_")[1]); v[c:] = 2.0
    return np.sqrt(v)

CELLS = ["ramp", "burst_30", "burst_60", "incal_100", "incal_500", "incal_900"]

def summarise(kind, method, a, d, N, m0hat=None):
    row = dict(cell=kind, method=method, paths=N)
    def rate(mask, key):
        k = int(mask.sum()); lo, hi = wilson(k, N)
        row[key] = k / N; row[key + "_ci"] = f"[{lo:.3f},{hi:.3f}]"; row[key + "_k"] = k
    det = a >= 0
    if kind.startswith("ramp"):
        s = 50
        rate(det & (a < s), "alarm_before_change")
        rate(det & (a >= s) & (a < s + 390), "detected_during_ramp")
        rate(det & (a >= s) & (a < s + 390 + 390), "detected_within_day_after_ramp_end")
        rate(det & (a >= s), "detected_within_5days")
        dl = a[det & (a >= s)] - s + 1
        row["median_delay_from_ramp_start"] = float(np.median(dl)) if dl.size else np.nan
        row["delay_iqr"] = f"{np.percentile(dl,25):.0f} to {np.percentile(dl,75):.0f}" if dl.size else ""
        row["share_up"] = float(np.mean(d[det & (a >= s)] > 0)) if dl.size else np.nan
    elif kind.startswith("burst"):
        L = int(kind.split("_")[1]); s = 200
        rate(det & (a < s), "alarm_before_burst")
        rate(det & (a >= s) & (a < s + L), "detected_inside_burst")
        rate(det & (a >= s + L) & (a < s + L + 60), "detected_within_60_after_burst")
        rate(det & (a >= s) & (a < s + L + 60), "detected_inside_or_within_60_after")
        rate(det & (a >= s + L + 60), "alarm_later_than_60_after_burst")
        dl = a[det & (a >= s)] - s + 1
        row["median_delay_from_burst_start"] = float(np.median(dl)) if dl.size else np.nan
        row["delay_iqr"] = f"{np.percentile(dl,25):.0f} to {np.percentile(dl,75):.0f}" if dl.size else ""
        row["share_up"] = float(np.mean(d[det & (a >= s)] > 0)) if dl.size else np.nan
    else:
        rate(det & (a < 390), "alarm_within_390")
        rate(det & (a < 1950), "alarm_within_1950")
        rate(det & (a < 1950) & (d > 0), "up_alarm_within_1950")
        rate(det & (a < 1950) & (d < 0), "down_alarm_within_1950")
        dl = a[det] + 1
        row["median_delay_from_monitoring_start"] = float(np.median(dl)) if dl.size else np.nan
        row["delay_iqr"] = f"{np.percentile(dl,25):.0f} to {np.percentile(dl,75):.0f}" if dl.size else ""
        if m0hat is not None:
            row["median_m0hat_over_true_new_m0"] = float(np.median(m0hat))
    return row

if __name__ == "__main__":
    N = int(sys.argv[1]); out = sys.argv[2]
    rows = []
    # true clipped mean of the new (x2) regime, for the in-calibration cells: Gaussian sd sqrt 2 with tau = 3 x sd
    g = np.random.default_rng(9_769_999).standard_normal(2_000_000) * np.sqrt(2)
    _, m0_new = truth(g)
    for j, kind in enumerate(CELLS):
        t0 = time.time()
        rng = np.random.default_rng(9_760_000 + 100 * j)
        data = rng.standard_normal((N, CAL + NMON)) * sd_profile(kind)[None, :]
        cal, mon = data[:, :CAL], data[:, CAL:]
        scale, tau2, m0p = plugin_stats(cal)
        for band, name in ((0.25, "monitor plugin band0.25 (README)"), (0.0, "monitor plugin band0")):
            a, d = vec_plugin(mon, tau2, m0p, band)
            rows.append(summarise(kind, name, a, d, N, m0p / m0_new if kind.startswith("incal") else None))
        a, d, seg = stream(lambda: VolatilityMonitor(threshold=A, cal_len=CAL, baseline="cs", delta=0.01), data)
        rows.append(summarise(kind, "monitor cs baseline delta0.01", a, d, N))
        pd.DataFrame(rows).to_csv(out, index=False)
        print(kind, "done in", round(time.time() - t0), "s", flush=True)
    print("EXIT", flush=True)
