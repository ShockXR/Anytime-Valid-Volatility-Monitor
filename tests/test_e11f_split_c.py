"""Review 2b: the E11f script (experiments/e11f_split_c.py).

* m14: the vectorised first alarms match VolatilityMonitor, alarm for alarm, on four laws.
* C2: the tick law rounds after the change of scale, so the tick is fixed in units of the pre-change sd.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_e11f():
    spec = importlib.util.spec_from_file_location("e11f_split_c", ROOT / "experiments" / "e11f_split_c.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_e11f_vectorised_first_alarms_match_the_class():
    """Four laws, sd ratios 1, 0.8 and 0.5, c_down 3, 2 and 1.5, 6 paths each: 216 runs, first alarm and side."""
    E = _load_e11f()
    cal, horizon, thr = 300, 1200, 1500.0
    rng = np.random.default_rng(2026)
    kinds = set()
    for law in ("gaussian", "t4", "jumps", "ticks"):
        for ratio in (1.0, 0.8, 0.5):
            if law == "ticks":
                r = E.tick_paths(rng.standard_normal((6, cal + horizon)), ratio, 0.5, "fixed", cal=cal)
            else:
                r = E.draw(law, (6, cal + horizon), rng) * np.r_[np.ones(cal), np.full(horizon, ratio)][None, :]
            for cd in (3.0, 2.0, 1.5):
                t_up, t_dn = E.first_alarms(r, cd, cal=cal, threshold=thr)
                for i in range(r.shape[0]):
                    mine = E.vectorised_first_alarm(t_up, t_dn, i)
                    assert E.class_first_alarm(r[i], cd, cal=cal, threshold=thr) == mine, (law, ratio, cd, i)
                    kinds.add(None if mine is None else mine[1])
    assert kinds == {None, "up", "down"}  # the check covers paths with an up alarm, a down alarm and no alarm


def test_e11f_tick_law_is_fixed_in_units_of_the_pre_change_sd():
    """Review 2b, C2: after a fall the share of zero returns rises (fixed tick); the old law kept it constant."""
    E = _load_e11f()
    z = np.random.default_rng(3).standard_normal((200, 2 * E.CAL))
    fixed = E.tick_paths(z, 0.5, 0.5, "fixed")
    old = E.tick_paths(z, 0.5, 0.5, "scaled")
    # zero shares: 2 Phi(0.25) - 1 = 0.197 before; 2 Phi(0.5) - 1 = 0.383 after with a fixed tick
    assert np.mean(fixed[:, :E.CAL] == 0) == pytest.approx(0.197, abs=0.01)
    assert np.mean(fixed[:, E.CAL:] == 0) == pytest.approx(0.383, abs=0.01)
    assert np.mean(old[:, E.CAL:] == 0) == pytest.approx(0.197, abs=0.01)
    after = np.unique(np.round(fixed[:, E.CAL:] / 0.5, 9))
    assert np.allclose(after, np.round(after))  # still on the grid of multiples of 0.5
    np.testing.assert_array_equal(fixed[:, :E.CAL], old[:, :E.CAL])
