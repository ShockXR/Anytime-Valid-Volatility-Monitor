"""Two-minute demo on synthetic data (no download needed).

Simulates one-minute returns with jumps whose volatility rises 50% at step
3000 and falls back at step 6000, runs the streaming monitor, and saves a
figure to results/figures/demo_monitor.png. All data here are synthetic.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor  # noqa: E402

rng = np.random.default_rng(42)
n = 9000
sd = np.ones(n)
sd[3000:6000] = 1.5
r = sd * rng.standard_normal(n)
r += (rng.random(n) < 0.003) * rng.normal(0, 10, n)  # rare jumps, same rate throughout

robust = VolatilityMonitor(threshold=7800, cal_len=1000).run(r)
gauss = GaussianLRMonitor(threshold=7800, cal_len=1000).run(r)

print("True changes: up at 3000, down at 6000")
print("Robust e-detector alarms:", [(e.index, e.direction) for e in robust])
print("Gaussian LR alarms:      ", [(e.index, e.direction) for e in gauss])

fig, ax = plt.subplots(figsize=(11, 3.5))
ax.plot(np.abs(r), lw=0.3, color="#999")
ax.plot(sd, color="k", lw=1.2, label="true sd")
for e in robust:
    ax.axvline(e.index, color="#c0392b" if e.direction == "up" else "#27ae60", lw=1.5)
for e in gauss:
    ax.axvline(e.index, color="#4a6fa5", lw=0.8, ls=":")
ax.set_ylim(0, 6)
ax.set_title("Demo [synthetic]: robust alarms (red up, green down) vs Gaussian LR alarms (blue dotted)")
ax.set_xlabel("minute")
ax.legend(loc="upper left")
fig.tight_layout()
out = Path(__file__).resolve().parents[1] / "results/figures/demo_monitor.png"
fig.savefig(out, dpi=130)
print("saved", out)
