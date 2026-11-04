"""R3 null laws (SYNTHETIC). Usage: python r3_null.py <n_paths> <outcsv> [law ...]
Law i (index in LAWS): paths default_rng(9_700_000 + 1000 i); truth sample default_rng(9_750_000 + 1000 i);
GARCH(1,1) fit for the fc_ laws: default_rng(9_750_500 + 1000 i)."""
import sys, time
import numpy as np, pandas as pd
from r3_common import *

def truth_sample(law, li):
    rng = np.random.default_rng(9_750_000 + 1000 * li)
    if "_to_t" in law and law.startswith("t"): return t_unit(rng, int(law[1:].split("_to_t")[0]), 3_000_000)
    if law == "seasonal_pipeline": return seasonal_pipeline(rng, 150).ravel()
    if law.startswith("seasonal_amp"): return seasonal_pipeline(rng, 150, amplitude=float(law[12:])).ravel()
    if law.startswith("fc_"): return None  # set in main (needs the fit)
    r = draw(law, rng, 300, 10000).ravel()
    return r[np.isfinite(r)]

if __name__ == "__main__":
    N = int(sys.argv[1]); out = sys.argv[2]; laws = sys.argv[3:] or LAWS
    calibrate_ticks(); print("ticks", TICK, flush=True)
    allrows = []
    for law in laws:
        li = LAWS.index(law); t0 = time.time()
        fit = None; note = ""
        if law.startswith("fc_"):
            fit = fit_for(law[3:], 9_750_500 + 1000 * li); note = f"garch11 fit (omega,a,b)={tuple(round(v,4) for v in fit)}"
            ts = forecast_pipeline(law[3:], np.random.default_rng(9_750_000 + 1000 * li), 200, 10000, fit).ravel()
        else:
            ts = truth_sample(law, li)
        tau, m0 = truth(ts[np.isfinite(ts)])
        rng = np.random.default_rng(9_700_000 + 1000 * li)
        data = draw(law, rng, N, CAL + NMON + (1000 if law.startswith("tick") else 0), fit)
        if law.startswith("tick"): note = f"tick={TICK[int(law[9:])/100]:.3f} sd; zero share={np.mean(data==0):.3f}"
        if law.startswith("missing"): note = f"NaN share={np.mean(~np.isfinite(data)):.3f}"
        rows = run_methods(law, data, (tau, m0), N, extra_note=note, use_stream=law.startswith("tick"))
        if "_to_t" in law and law.startswith("t"):
            # extra: alarms after the switch at monitored step 975, among paths silent before it
            half = NMON // 2
            cal, mon = data[:, :CAL], data[:, CAL:]
            scale, tau2, m0p = plugin_stats(cal)
            for name, (a, d) in {"e-detector, KNOWN baseline (theorem)": known_baseline(mon, tau, m0),
                                 "monitor plugin band0.25 (README)": vec_plugin(mon, tau2, m0p, 0.25),
                                 "monitor plugin band0": vec_plugin(mon, tau2, m0p, 0.0)}.items():
                silent = (a < 0) | (a >= half)
                for m in (390, 975):
                    hit = silent & (a >= half) & (a < half + m)
                    k = int(hit.sum()); n = int(silent.sum()); lo, hi = wilson(k, n)
                    allrows.append(dict(law=law, method=name, horizon_m=f"switch+{m}", k=k, paths=n, rate=k / n, ci_low=lo, ci_high=hi, bound=m / A,
                                        k_up=int((hit & (d > 0)).sum()), k_down=int((hit & (d < 0)).sum()),
                                        verdict="EXCEEDS bound" if lo > m / A else ("within bound" if hi <= m / A else "inconclusive"),
                                        note="paths silent before the switch; bound as if a fresh segment"))
        allrows += rows
        pd.DataFrame(allrows).to_csv(out, index=False)
        print(law, "done in", round(time.time() - t0), "s", flush=True)
    print("EXIT", flush=True)
