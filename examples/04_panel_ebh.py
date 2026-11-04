"""A panel of assets monitored together, with discoveries controlled by e-BH.

[SYNTHETIC] Eight assets with independent Gaussian returns, two blocks ("days") of 390 steps after a
390-step warm-up. In the second block the volatility of asset_2 doubles at step 100 of the block and the
volatility of asset_5 falls by half at step 200; the other six assets never change.

``PanelMonitor`` gives every asset a two-sided mixture e-process, recalibrated at the start of each block
on the asset's own last ``cal_len`` returns. At every step it applies e-BH to the e-values, each stopped
when its asset is discovered. Within a block the false discovery rate is then at most alpha at every
stopping time, under any dependence between the assets (with the true baselines; ``docs/theory_panel.md``).
``rule="bonferroni"`` instead flags an asset when its e-process reaches K / alpha (family-wise error at
most alpha per block); path by path e-BH flags no later than Bonferroni. The false discovery rate is an
expectation, so a single run can still contain a false discovery. The baselines here are plug-in
estimates from 390 returns, which makes false discoveries more likely than with the true baselines;
``baseline="cs"`` keeps the guarantee with estimated baselines, at the cost of slower discoveries.

Run: python examples/04_panel_ebh.py
"""

import numpy as np

from avmonitor import PanelMonitor, ebh

K = 8
BLOCK = 390
ALPHA = 0.1


def simulate(rng: np.random.Generator) -> tuple[np.ndarray, dict]:
    n = BLOCK + 2 * BLOCK  # warm-up, then two blocks
    sd = np.ones((n, K))
    up_at, down_at = 2 * BLOCK + 100, 2 * BLOCK + 200  # steps 100 and 200 of the second block
    sd[up_at:, 2] = 2.0
    sd[down_at:, 5] = 0.5
    return sd * rng.standard_normal((n, K)), {2: up_at, 5: down_at}


def main() -> None:
    returns, changes = simulate(np.random.default_rng(4))
    names = [f"asset_{k}" for k in range(K)]
    print(f"[SYNTHETIC] {K} assets, {returns.shape[0]} steps; blocks of {BLOCK} steps after a {BLOCK}-step warm-up")
    print("changes: " + ", ".join(f"{names[k]} at step {t}" for k, t in changes.items()) + "; no other asset changes")
    e = [50.0, 25.0, 1.0, 1.0]
    print(f"e-BH on the e-values {e} at alpha = {ALPHA} rejects {ebh(e, ALPHA).tolist()}: the second one is "
          f"below the Bonferroni level {len(e) / ALPHA:.0f} but above the e-BH level {len(e) / (ALPHA * 2):.0f}")

    for rule in ("ebh", "bonferroni"):
        mon = PanelMonitor(names, alpha=ALPHA, rule=rule, block_len=BLOCK, cal_len=BLOCK)
        found = mon.run(returns)
        print(f"rule={rule!r}: {len(found)} discoveries")
        for d in found:
            if d.asset in changes and d.index >= changes[d.asset]:
                note = f"true discovery, {d.index - changes[d.asset] + 1} steps after the change"
            else:
                note = "false discovery (no change)"
            print(f"    {d.name}: step {d.index} (block {d.block}), direction {d.direction}, "
                  f"e-value {d.evalue:.0f} >= level {d.threshold:.0f}; {note}")


if __name__ == "__main__":
    main()
