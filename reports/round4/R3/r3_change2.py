"""R3 extra change cell (SYNTHETIC): the variance doubles or triples from monitored step 50 onwards, but only through rare
large jumps (probability p per step, jump sd 20, so p x 400 = 1 or 2 added variance), which the truncation clips.
Usage: python r3_change2.py <n_paths> <outcsv>. Cell j uses default_rng(9_761_000 + 100 j)."""
import sys, time
import numpy as np, pandas as pd
from r3_common import *
from r3_change import summarise

if __name__ == "__main__":
    N = int(sys.argv[1]); out = sys.argv[2]; rows = []
    for j, (kind, p) in enumerate((("jumpvar_x2", 0.0025), ("jumpvar_x3", 0.005))):
        rng = np.random.default_rng(9_761_000 + 100 * j)
        data = rng.standard_normal((N, CAL + NMON))
        s = CAL + 50
        data[:, s:] += (rng.random((N, NMON - 50)) < p) * rng.normal(0, 20.0, (N, NMON - 50))
        cal, mon = data[:, :CAL], data[:, CAL:]
        scale, tau2, m0p = plugin_stats(cal)
        for band, name in ((0.25, "monitor plugin band0.25 (README)"), (0.0, "monitor plugin band0")):
            a, d = vec_plugin(mon, tau2, m0p, band)
            r = summarise("incal_0", name, a, d, N); r["cell"] = kind; r["true_var_after"] = 1 + 400 * p
            r["clipped_mean_ratio_after_over_before"] = float(np.mean(np.minimum(mon[:, 50:]**2, tau2[:, None])) / np.mean(np.minimum(mon[:, :50]**2, tau2[:, None])))
            rows.append(r)
        a, d, seg = stream(lambda: VolatilityMonitor(threshold=A, cal_len=CAL, baseline="cs", delta=0.01), data)
        r = summarise("incal_0", "monitor cs baseline delta0.01", a, d, N); r["cell"] = kind; rows.append(r)
        pd.DataFrame(rows).to_csv(out, index=False); print(kind, "done", flush=True)
    print("EXIT", flush=True)
