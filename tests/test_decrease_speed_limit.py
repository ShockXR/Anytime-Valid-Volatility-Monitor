"""Review 2b, C1: the decrease speed limit for every valid bet and for the default grid.

The floor of 73 steps for a total collapse of volatility at c = 3 and A = 7800 belongs to the default
five-bet decrease grid. Any Shiryaev-Roberts e-detector with factors valid for the whole null needs at least
ceil(log(1 + A / B) / log(B / (B - 1))) steps, which is 58 at c = 3.
"""

import numpy as np
import pytest

from avmonitor import design as G
from avmonitor import detectors as D


def _gaussian_B(c):
    tau, m0 = G.baseline(G.ReturnModel("gaussian"), c)
    return tau**2 / m0


# c: (default grid at A, at 2A, any valid bet at A, at 2A), A = 7800 (docs/theory_delay.md, Corollary 4)
FLOORS = {1.5: (23, 25, 19, 21), 2.0: (36, 39, 29, 32), 3.0: (73, 79, 58, 64), 4.0: (122, 134, 96, 107),
          5.0: (180, 199, 141, 158)}


@pytest.mark.parametrize("c", sorted(FLOORS))
def test_decrease_floors_default_grid_and_any_valid_bet(c):
    B = _gaussian_B(c)
    A = 7800.0
    grid = 1 + np.asarray(D.DEFAULT_DOWN_FRACTIONS) / (B - 1)
    w = np.full(grid.size, 1 / grid.size)
    got = (G.speed_limit(A, grid, w), G.speed_limit(2 * A, grid, w),
           G.decrease_speed_limit(A, B), G.decrease_speed_limit(2 * A, B))
    assert got == FLOORS[c]
    # closed form against the recursion, and the floor is below the first-order log A / log(B / (B - 1))
    F = B / (B - 1)
    assert G.decrease_speed_limit(A, B) == G.speed_limit(A, [F], [1.0])
    assert G.decrease_speed_limit(A, B) < np.log(A) / np.log(F) + 1


def test_c3_numbers_quoted_in_the_text():
    B = _gaussian_B(3.0)
    assert B == pytest.approx(9.04516, abs=1e-5)
    F = B / (B - 1)
    assert np.log(F) == pytest.approx(0.11716, abs=1e-5)
    assert np.log1p(7800 / B) / np.log(F) == pytest.approx(57.7, abs=0.05)
    assert np.log1p(15600 / B) / np.log(F) == pytest.approx(63.6, abs=0.05)
    assert np.log(7800) / np.log(F) == pytest.approx(76.5, abs=0.05)


def test_single_bet_near_lam_max_attains_the_floor_on_zero_returns():
    B, A = _gaussian_B(3.0), 7800.0
    lam = 0.999 * G.lam_max("down", B)
    alarm = G.run_mixture(np.zeros((1, 200)), [lam], [1.0], "down", A)["alarm"][0]
    assert alarm + 1 == 58 == G.decrease_speed_limit(A, B)
    grid = G.lam_max("down", B) * np.asarray(D.DEFAULT_DOWN_FRACTIONS)
    w = np.full(grid.size, 1 / grid.size)
    assert G.run_mixture(np.zeros((1, 200)), grid, w, "down", A)["alarm"][0] + 1 == 73


def test_sr_steps_constant_factor_matches_the_recursion():
    for F in (1.01, 1.1, 1.5, 3.0):
        for A in (2.0, 50.0, 7800.0):
            assert G.sr_steps_constant_factor(A, F) == G.speed_limit(A, [F], [1.0])
    with pytest.raises(ValueError):
        G.sr_steps_constant_factor(10.0, 1.0)
    with pytest.raises(ValueError):
        G.decrease_speed_limit(10.0, np.inf)
