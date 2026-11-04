"""E19 (experiments/e19_real_null.py): the vectorised monitors are the package's monitors, and the streams are
what the pre-registration says.

* ``x_engine`` matches ``VolatilityMonitor`` alarm for alarm (plug-in and cs, bands 0 and 0.25), on streams with
  volatility shifts and long stretches of zero returns, so that the extension and fallback rules for degenerate
  calibration windows are used;
* ``rank_engine`` matches ``ConformalMonitor(restart_len=0)`` (awake-B) alarm for alarm, with the same tie-breakers;
* the known-baseline mode matches a direct Shiryaev-Roberts recursion with restarts at once;
* permutation streams are permutations (no index twice), the same seed gives the same streams, the stationary
  bootstrap has the stated mean block length;
* the calm-day rule uses only the day itself and the 20 days before it;
* the vectorised confidence-sequence split clip of Part 4 matches the class;
* (post hoc) when a cs calibration gives tau^2 <= m_up the increase side cannot alarm, and ``side_activity``
  reports the segment as switched off.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from avmonitor.monitor import VolatilityMonitor
from avmonitor.running import ConformalMonitor

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def E():
    spec = importlib.util.spec_from_file_location("e19_real_null", ROOT / "experiments" / "e19_real_null.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _streams(seed=1, n=8, L=4000):
    rng = np.random.default_rng(seed)
    r = rng.standard_normal((n, L)) * np.repeat(rng.choice([0.6, 1.0, 1.6], (n, L // 400)), 400, axis=1)
    r[:, 700:1000] = np.where(rng.random((n, 300)) < 0.7, 0.0, r[:, 700:1000])  # degenerate windows
    r[2] = np.where(rng.random(L) < 0.55, 0.0, r[2])  # MAD = 0 almost everywhere: fallback scale
    return r


def _summary(update, calibrated, r, m):
    """(first alarm as monitored count, its side, up alarms, down alarms within m, monitored steps)."""
    mcount = first = fdir = nu = nd = 0
    for x in r:
        was = calibrated()
        ev = update(x)
        if was:
            mcount += 1
            if ev is not None:
                if first == 0:
                    first, fdir = mcount, (1 if ev.direction == "up" else -1)
                if mcount <= m:
                    nu += ev.direction == "up"
                    nd += ev.direction == "down"
            if mcount >= m:
                break
    return first, fdir, nu, nd, mcount


@pytest.mark.parametrize("mode", ["plugin", "cs"])
@pytest.mark.parametrize("band", [0.0, 0.25])
def test_x_engine_matches_the_class(E, mode, band):
    r = _streams()
    cal, A, m = 160, 300.0, 1800
    R = E.x_engine(None, E.ArrayGen(r), mode, band, A, cal, m)
    stats = {"extended": 0, "fallback": 0}
    alarms = 0
    for i in range(r.shape[0]):
        mon = VolatilityMonitor(threshold=A, cal_len=cal, null_band=band, baseline=mode, delta=0.01)
        ref = _summary(mon.update, lambda: mon.calibrated, r[i], m)
        got = tuple(int(R[k][i]) for k in ("first_t", "first_dir", "n_up", "n_dn", "mon_steps"))
        assert got == ref, (mode, band, i)
        stats["extended"] += mon.cal_stats["extended"]
        stats["fallback"] += mon.cal_stats["fallback"]
        alarms += ref[2] + ref[3]
    assert alarms > 10  # the comparison covers many alarms and recalibrations
    assert stats["extended"] > 0 and stats["fallback"] > 0  # and the degenerate-window rules


@pytest.mark.parametrize("band", [0.0, 0.25])
def test_rank_engine_matches_the_awake_rank_monitor(E, band):
    r = _streams(seed=3)
    cal, A, m = 160, 300.0, 1800
    n, L = r.shape
    _, inv = np.unique(np.abs(r).ravel(), return_inverse=True)
    ranks = (inv + 1).reshape(r.shape)
    theta = np.zeros((n, L))
    for i in range(n):
        theta[i, cal:] = np.random.default_rng(500 + i).random(L - cal)
    R = E.rank_engine(ranks, E.ArrayGen(r), band, A, cal, m, theta=theta)
    alarms = 0
    for i in range(n):
        mon = ConformalMonitor(threshold=A, cal_len=cal, band=band, seed=500 + i, restart_len=0)
        ref = _summary(mon.update, lambda: mon.calibrated, r[i], m)
        got = tuple(int(R[k][i]) for k in ("first_t", "first_dir", "n_up", "n_dn", "mon_steps"))
        assert got == ref, (band, i)
        alarms += ref[2] + ref[3]
    assert alarms > 10


def test_known_baseline_matches_a_direct_recursion(E):
    rng = np.random.default_rng(7)
    n, L, cal, A, m = 6, 3000, 100, 200.0, 2500
    r = rng.standard_t(4, (n, L)) * np.repeat(rng.choice([0.7, 1.0, 1.4], (n, L // 300)), 300, axis=1)
    tau2, m0 = 9.0 * 0.9, 0.85
    R = E.x_engine(None, E.ArrayGen(r), "known", 0.25, A, cal, m, known=(tau2, m0))
    up = np.asarray(DEFAULT_UP_GRID)
    lam = np.asarray(DEFAULT_DOWN_FRACTIONS) / (tau2 / (0.75 * m0) - 1.0)
    total = 0
    for i in range(n):
        Su, Sd = np.zeros(up.size), np.zeros(lam.size)
        k = first = fdir = nu = nd = 0
        for x in r[i, cal:]:
            k += 1
            X = min(x * x, tau2)
            Su = (1 + up * (X / (1.25 * m0) - 1)) * (Su + 1)
            Sd = (1 + lam * (1 - X / (0.75 * m0))) * (Sd + 1)
            side = 1 if Su.mean() >= 2 * A else (-1 if Sd.mean() >= 2 * A else 0)
            if side:
                first, fdir = (k, side) if first == 0 else (first, fdir)
                nu += side == 1
                nd += side == -1
                Su[:], Sd[:] = 0.0, 0.0  # restart at once
            if k >= m:
                break
        got = tuple(int(R[c][i]) for c in ("first_t", "first_dir", "n_up", "n_dn", "mon_steps"))
        assert got == (first, fdir, nu, nd, k)
        total += nu + nd
    assert total > 5


def test_streams(E):
    N, n = 5000, 4
    g1 = E.StreamGen("perm", N, n, 3000, np.random.default_rng(1))
    g2 = E.StreamGen("perm", N, n, 3000, np.random.default_rng(1))
    a = np.array([g1.next() for _ in range(3000)]).T
    b = np.array([g2.next() for _ in range(3000)]).T
    assert np.array_equal(a, b)
    assert all(np.unique(row).size == row.size for row in a)
    assert a.min() >= 0 and a.max() < N
    g = E.StreamGen("perm", 50, 2, 1000, np.random.default_rng(2))  # a stream that runs out of pool
    c = np.array([g.next() for _ in range(60)])
    assert (c[:50] >= 0).all() and (c[50:] == -1).all()
    g = E.StreamGen("sb390", N, 50, 10**9, np.random.default_rng(3))
    d = np.array([g.next() for _ in range(20000)])
    jumps = np.sum(d[1:] != (d[:-1] + 1) % N)
    assert 0.8 < (d.size - d.shape[1]) / jumps / 390 < 1.25  # mean block length 390 (new blocks may continue by chance)


def test_calm_days_use_only_past_days(E):
    rng = np.random.default_rng(4)
    days = 60
    lev = np.ones(days)
    lev[40] = 2.0  # one volatile day
    keys = np.repeat(np.arange(days), 1440)
    u = rng.standard_normal(days * 1440) * np.repeat(lev, 1440)
    in_period = keys >= 30
    tab, mask = E.calm_flags(keys, u, in_period, 1000)
    assert list(tab.day) == list(range(30, 60))
    assert not tab.set_index("day").loc[40, "calm"]
    assert tab.calm.sum() >= 25  # i.i.d. days with the same level are almost all calm
    u2 = u.copy()
    u2[50 * 1440:] *= 5.0  # later days do not change earlier flags
    tab2, _ = E.calm_flags(keys, u2, in_period, 1000)
    assert tab2.set_index("day").loc[30:49, "calm"].equals(tab.set_index("day").loc[30:49, "calm"])


def test_split_clip_cs_matches_the_class(E):
    rng = np.random.default_rng(11)
    cal = 400
    r = rng.standard_normal((6, cal + 1500)) * np.r_[np.ones(cal), np.full(1500, 0.6)][None, :]
    for cd in (3.0, 2.0, 1.5):
        tu, td = E.first_alarms_cs(r, cd, cal=cal, threshold=2000.0)
        for i in range(r.shape[0]):
            ev = VolatilityMonitor(threshold=2000.0, cal_len=cal, baseline="cs", delta=0.01, c_down=cd).run(r[i])
            ref = (ev[0].index - cal, ev[0].direction) if ev else None
            got = (int(tu[i]), "up") if tu[i] >= 0 else ((int(td[i]), "down") if td[i] >= 0 else None)
            assert got == ref, (cd, i)


def test_stacked_configurations_equal_single_runs(E):
    """x_engine_multi runs several configurations in one pass; each must equal its own single run."""
    rng = np.random.default_rng(21)
    pool = rng.standard_t(4, 20000) * np.repeat(rng.lognormal(0, 0.3, 50), 400)
    pool[rng.random(pool.size) < 0.04] = 0.0
    tau2, m0 = 9.0 * np.median(np.abs(pool)) ** 2 * 1.4826**2, None
    m0 = float(np.mean(np.minimum(pool**2, tau2)))
    configs = [("known", 0.0, 300.0), ("plugin", 0.0, 300.0), ("cs", 0.25, 400.0), ("plugin", 0.25, 300.0)]
    cal, m = 150, 1500
    multi = E.x_engine_multi(pool, E.StreamGen("sb30", pool.size, 7, 10**6, np.random.default_rng(3)), configs, cal, m,
                             known=(tau2, m0))
    for cfg, R in zip(configs, multi):
        single = E.x_engine(pool, E.StreamGen("sb30", pool.size, 7, 10**6, np.random.default_rng(3)), *cfg, cal, m,
                            known=(tau2, m0))
        for key in ("first_t", "first_dir", "n_up", "n_dn", "mon_steps", "used", "calibrations"):
            assert np.array_equal(R[key], single[key]), (cfg, key)
    assert sum(int((R["n_up"] + R["n_dn"]).sum()) for R in multi) > 5


def test_increase_side_switched_off_when_tau2_below_m_up(E):
    """Post hoc diagnostic: a cs calibration whose first quarter is quiet gives tau^2 <= U (1 + b); the increase
    statistic then never exceeds the number of monitored steps (each factor is at most 1), so no up alarm comes, and
    side_activity counts the whole monitored stretch as switched off."""
    rng = np.random.default_rng(5)
    cal = 400
    r = np.r_[0.2 * rng.standard_normal(cal // 4), rng.standard_normal(3 * cal // 4), 3.0 * rng.standard_normal(3000)]
    mon = VolatilityMonitor(threshold=36000.0, cal_len=cal, null_band=0.5, baseline="cs", delta=0.01)
    ev = mon.run(r)
    seg = mon.segments[0]
    assert seg["tau2"] <= seg["m_up"]
    assert not any(e.direction == "up" for e in ev)
    assert np.all(mon.S_up <= len(r) - seg["start"] + 1)
    act = E.side_activity(mon, len(r))
    assert act["segments_up_off"] == act["segments"] == 1
    assert act["up_off_minutes"] == act["monitored_minutes"] == len(r) - seg["start"]
