"""Round 9 check of the panel monitor under the global null (known baselines). [SYNTHETIC]

panel.run_block with K = 5 assets, true baselines, one block of 390 steps, e-BH at alpha = 0.05, window
weights (h = 390): P(any discovery in the block) <= alpha (docs/theory_panel.md). Gaussian and t3 returns,
plus 10% missing returns. Also PanelMonitor (streaming) against run_block on the same data.
"""
import math

import numpy as np
from scipy import integrate, stats

from avmonitor.panel import PanelMonitor, run_block
from avmonitor.robustness import trunc_mean_normal

M0 = float(trunc_mean_normal(3.0))
TAU_T = 3.0 * stats.t.ppf(0.75, 3) / stats.norm.ppf(0.75)
M0_T = integrate.quad(lambda v: 2 * stats.t.sf(math.sqrt(v), 3), 0, TAU_T**2, limit=200)[0]


def mc(law, nan_frac, reps=20000, K=5, n=390, alpha=0.05):
    rng = np.random.default_rng(3)
    if law == "gauss":
        r, tau2, m0 = rng.standard_normal((reps, K, n)), 9.0, M0
    else:
        r, tau2, m0 = rng.standard_t(3, (reps, K, n)), TAU_T**2, M0_T
    if nan_frac:
        r[rng.random(r.shape) < nan_frac] = np.nan
    res = run_block(r, tau2, m0, m0, alpha=alpha, rule="ebh")
    p = float((res.flag_time >= 0).any(axis=1).mean())
    se = math.sqrt(max(p * (1 - p), 1e-300) / reps)
    print(f"[SYNTHETIC] panel e-BH, K={K}, n={n}, {law}, {nan_frac:.0%} missing: P(any discovery) = {p:.4%} "
          f"(se {se:.4%}) over {reps} panels, bound {alpha:.0%}", flush=True)


def stream_vs_block():
    rng = np.random.default_rng(4)
    K, cal, n = 4, 300, 600
    diffs = 0
    tot = 0
    for i in range(20):
        r = rng.standard_normal((cal + n, K))
        r[cal + 200:, 0] *= 2.0
        r[cal + 100:, 2] *= 0.4
        r[cal:][rng.random((n, K)) < 0.05] = np.nan  # missing returns in the monitored part only
        mon = PanelMonitor(K, alpha=0.1, block_len=n, cal_len=cal, warmup=cal)
        disc = mon.run(r)
        b = mon.blocks[0]
        res = run_block(r[cal:].T, b.tau2, b.m_up, b.m_dn, alpha=0.1, rule="ebh", horizon=n, active=b.active)
        a = sorted((d.asset, d.index - cal) for d in disc)
        c = sorted((k, int(res.flag_time[0, k])) for k in range(K) if res.flag_time[0, k] >= 0)
        diffs += a != c
        tot += len(a)
        assert b.active.all()
    print(f"PanelMonitor against run_block on the same data: {diffs} of 20 panels differ ({tot} discoveries)")


if __name__ == "__main__":
    stream_vs_block()
    for law in ("gauss", "t3"):
        mc(law, 0.0)
    mc("gauss", 0.10)
