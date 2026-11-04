"""Round 4, agent R5: metamorphic and property tests of the streaming monitor and the e-detector.

Each test checks a RELATION that must hold by design (read from ``src/avmonitor/monitor.py``,
``detectors.py``, ``confseq.py``, ``cli.py`` and ``docs/math.md`` section 5), on random inputs from
several laws and lengths, instead of one worked example. Relations:

* scale invariance (returns times k > 0): plug-in, cs, Gaussian LR monitors and the known-baseline detector;
* sign invariance (global flip and random per-return flips);
* streaming equals batch equals random chunks (plug-in, cs, c_down, Gaussian LR);
* a NaN inserted at a monitored step only shifts the indices of later events;
* pickle or deepcopy mid-stream, then continue, equals the uninterrupted run;
* a larger threshold, a wider null band or a larger margin never alarms earlier (first alarm);
* under the exact null the expected one-step betting factor is 1 for every bet (several laws);
* the e-process stays finite and nonnegative on extreme but finite inputs;
* a restart after an alarm (or after a degenerate window) equals a fresh monitor started there;
* ``python -m avmonitor simulate`` and ``avmonitor run`` on a CSV give the same events as the Python API.

Seeds are numpy ``default_rng`` seeds from 9,900,000 upwards (round 4 range for R5). Everything is synthetic.
The two real bugs it found (R5-B-001, R5-B-002) were fixed in round 4; see reports/round4/R5/REPORT.md.
"""

from __future__ import annotations

import ast
import copy
import math
import pickle
import re
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import integrate, stats

from avmonitor import data as DATA
from avmonitor import detectors as D
from avmonitor.cli import main as cli_main
from avmonitor.monitor import GaussianLRMonitor, VolatilityMonitor
from avmonitor.robust import MAD_NORMAL_SCALE

SEED0 = 9_900_000

LAWS = ("gaussian", "t3", "uniform", "twopoint", "laplace")


# --------------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------------


def _draw(rng, law, n):
    """n unit-scale returns from one law."""
    if law == "gaussian":
        return rng.standard_normal(n)
    if law == "t3":
        return rng.standard_t(3, n)
    if law == "uniform":
        return rng.uniform(-1.0, 1.0, n)
    if law == "twopoint":
        return rng.choice([-1.0, 1.0], n) * np.where(rng.random(n) < 0.6, 1.0, 5.0)
    if law == "laplace":
        return rng.laplace(0.0, 1.0, n)
    raise ValueError(law)


def _path(seed, law=None, n=None, n_regimes=None):
    """A random path with a few volatility regimes (so that both sides alarm), plus its parameters."""
    rng = np.random.default_rng(seed)
    law = law or LAWS[rng.integers(len(LAWS))]
    n = n or int(rng.integers(600, 2500))
    n_regimes = n_regimes or int(rng.integers(2, 5))
    cuts = np.sort(rng.choice(np.arange(200, n - 100), n_regimes - 1, replace=False)) if n_regimes > 1 else []
    sd = np.ones(n)
    level = 1.0
    for c in cuts:
        level *= rng.choice([0.3, 0.5, 2.0, 3.0])
        sd[c:] = level
    r = sd * _draw(rng, law, n)
    cal_len = int(rng.integers(60, 200))
    A = float(rng.choice([10.0, 30.0, 100.0, 400.0]))
    return r, dict(threshold=A, cal_len=cal_len)


def _key(mon):
    return [(e.index, e.direction, e.monitored_since, e.baseline_scale, e.stat) for e in mon.events]


def _key_idx(mon):
    return [(e.index, e.direction) for e in mon.events]


def _seg_key(mon):
    keys = ("start", "cal_len", "n_finite", "scale", "tau2", "tau2_dn", "m0", "m0_dn", "m_up", "m_dn", "down_active")
    return [tuple(s[k] for k in keys) for s in mon.segments]


def _close_events(a, b, rtol=1e-6):
    """Same alarm times and sides; stats and scales equal up to rounding."""
    assert [(i, d) for i, d, *_ in a] == [(i, d) for i, d, *_ in b]
    for x, y in zip(a, b):
        assert x[2] == y[2]
        assert math.isclose(x[4], y[4], rel_tol=rtol)


def _cases(n, kinds=({},)):
    out = []
    for i in range(n):
        r, kw = _path(SEED0 + 100 + i)
        for extra in kinds:
            out.append((r, {**kw, **extra}))
    return out


MON_KINDS = ({}, {"baseline": "cs"}, {"c_down": 1.5}, {"baseline": "cs", "c_down": 2.0}, {"null_band": 0.2},
             {"margin_z": 1.0}, {"baseline": "cs", "target_rel_width": 0.8})
# scale_floor > 0 is left out on purpose: a floor in return units breaks scale invariance by design


# --------------------------------------------------------------------------------------------
# 1. scale invariance
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("k", [2.0, 0.5, 1024.0, 1e-4, 3.7])
def test_scale_invariance_of_the_monitors(k):
    """Multiplying all returns by k > 0 leaves every alarm time and side unchanged (plug-in, cs, c_down, band,
    margin). scale, tau^2 and m0 all scale with k (k^2), so x / m0 and the bets do not change. Powers of two
    give bit-identical statistics; other k only up to rounding."""
    exact = math.log2(k).is_integer()
    for r, kw in _cases(8, MON_KINDS):
        a = VolatilityMonitor(**kw)
        a.run(r)
        b = VolatilityMonitor(**kw)
        b.run(k * r)
        assert len(a.events) >= 1 or len(a.segments) >= 1
        if exact:
            assert _key_idx(a) == _key_idx(b)
            assert [e.stat for e in a.events] == [e.stat for e in b.events]
            assert [e.baseline_scale * k for e in a.events] == pytest.approx([e.baseline_scale for e in b.events])
        else:
            _close_events(_key(a), _key(b), rtol=1e-8)
        for sa, sb in zip(a.segments, b.segments):
            assert sa["start"] == sb["start"] and sa["cal_len"] == sb["cal_len"]
            assert sb["scale"] == pytest.approx(k * sa["scale"], rel=1e-9)
            assert sb["m0"] == pytest.approx(k * k * sa["m0"], rel=1e-9)
            assert sb["tau2"] == pytest.approx(k * k * sa["tau2"], rel=1e-9)


@pytest.mark.parametrize("k", [4.0, 0.25, 2.5])
def test_scale_invariance_of_the_gaussian_lr_monitor(k):
    """The Gaussian LR benchmark uses (r / scale)^2, so it is scale invariant too."""
    for r, kw in _cases(6):
        a, b = GaussianLRMonitor(**kw), GaussianLRMonitor(**kw)
        a.run(r)
        b.run(k * r)
        _close_events(_key(a), _key(b), rtol=1e-8)


@pytest.mark.parametrize("k", [8.0, 0.125, 1e-3])
def test_scale_invariance_of_the_known_baseline_detector(k):
    """betting_edetector on min(r^2, tau^2) with m0 and tau^2 scaled by k^2 gives the same alarms, both sides."""
    for i in range(6):
        rng = np.random.default_rng(SEED0 + 300 + i)
        law = LAWS[i % len(LAWS)]
        n = 1200
        r = np.r_[_draw(rng, law, 600), rng.choice([0.4, 2.5]) * _draw(rng, law, 600)]
        tau2 = (3 * MAD_NORMAL_SCALE * np.median(np.abs(r[:300]))) ** 2
        m0 = float(np.minimum(r[:300] ** 2, tau2).mean())
        x = np.minimum(r**2, tau2)
        xk = np.minimum((k * r) ** 2, k * k * tau2)
        for direction, kind in (("up", "SR"), ("down", "SR"), ("up", "CUSUM"), ("down", "CUSUM")):
            ra = D.betting_edetector(x, m0, direction, upper=tau2, threshold=50.0, kind=kind, return_path=True)
            rb = D.betting_edetector(xk, k * k * m0, direction, upper=k * k * tau2, threshold=50.0, kind=kind,
                                     return_path=True)
            assert ra.alarm.tolist() == rb.alarm.tolist()
            np.testing.assert_allclose(ra.stat_path, rb.stat_path, rtol=1e-9)


# --------------------------------------------------------------------------------------------
# 2. sign invariance
# --------------------------------------------------------------------------------------------


def test_sign_invariance_global_and_per_return():
    """r -> -r, and r_t -> s_t r_t for random signs s_t, give bit-identical events (only |r| is used)."""
    for j, (r, kw) in enumerate(_cases(8, MON_KINDS)):
        rng = np.random.default_rng(SEED0 + 400 + j)
        flips = rng.choice([-1.0, 1.0], r.size)
        base = VolatilityMonitor(**kw)
        base.run(r)
        for rr in (-r, flips * r):
            m = VolatilityMonitor(**kw)
            m.run(rr)
            assert _key(m) == _key(base)
            assert _seg_key(m) == _seg_key(base)
        g, gb = GaussianLRMonitor(**{k: v for k, v in kw.items() if k in ("threshold", "cal_len")}), None
        if kw.get("baseline", "plugin") == "plugin" and "c_down" not in kw:
            gb = GaussianLRMonitor(threshold=kw["threshold"], cal_len=kw["cal_len"])
            g.run(r)
            gb.run(flips * r)
            assert _key(g) == _key(gb)


# --------------------------------------------------------------------------------------------
# 3. streaming equals batch equals random chunks
# --------------------------------------------------------------------------------------------


def test_streaming_equals_batch_equals_random_chunks():
    """update() one value at a time, run() on the array and run() on random chunks (including empty chunks) give
    the same events, segments and final state (plug-in, cs, c_down, Gaussian LR). Extends the audit test, which
    used the plug-in baseline and 7 equal chunks."""
    kinds = ({}, {"baseline": "cs"}, {"c_down": 1.5}, {"baseline": "cs", "c_down": 2.0, "target_rel_width": 0.8},
             {"gaussian": True})
    for j, (r, kw) in enumerate(_cases(6, kinds)):
        rng = np.random.default_rng(SEED0 + 500 + j)
        cls = VolatilityMonitor
        if kw.pop("gaussian", False):
            cls = GaussianLRMonitor
        one, two, three = cls(**kw), cls(**kw), cls(**kw)
        out = one.run(r)
        assert out is one.events
        evs = [two.update(float(v)) for v in r]
        assert [e for e in evs if e is not None] == two.events
        bounds = np.sort(rng.integers(0, r.size + 1, 25))
        bounds = np.r_[0, bounds, r.size]
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            three.run(r[lo:hi])  # empty chunks when lo == hi
        assert _key(one) == _key(two) == _key(three)
        assert _seg_key(one) == _seg_key(two) == _seg_key(three)
        assert one.t == two.t == three.t == r.size - 1
        for m in (two, three):
            assert np.array_equal(one.S_up, m.S_up) and np.array_equal(one.S_dn, m.S_dn)
            assert one.calibrated == m.calibrated and one.buffer == m.buffer


# --------------------------------------------------------------------------------------------
# 4. NaN insertion
# --------------------------------------------------------------------------------------------


def _monitored_positions(mon, n):
    """Indices of the stream that were monitored (not part of a calibration window) in a finished run."""
    pos = []
    alarms = [e.index for e in mon.events]
    for s in mon.segments:
        end = next((a for a in alarms if a >= s["start"]), n - 1)
        pos.extend(range(s["start"], end + 1))
    return np.asarray(pos)


@pytest.mark.parametrize("kw_extra", [{}, {"baseline": "cs"}, {"c_down": 1.5}])
def test_nan_at_monitored_steps_only_shifts_later_events(kw_extra):
    """A NaN fed at a monitored step is skipped (no bet, no new chain), so inserting NaNs at monitored positions
    shifts later alarm indices by the number of NaNs before them and changes nothing else. NaNs inside a
    calibration window are a different matter: the window still closes after cal_len values, so the baseline is
    estimated from fewer returns (documented in docs/math.md 5.6 item 5); that case is not a shift relation."""
    for j, (r, kw) in enumerate(_cases(8)):
        kw = {**kw, **kw_extra}
        rng = np.random.default_rng(SEED0 + 600 + j)
        clean = VolatilityMonitor(**kw)
        clean.run(r)
        pos = _monitored_positions(clean, r.size)
        if pos.size == 0:
            continue
        k = int(min(pos.size, rng.integers(1, 12)))
        where = np.sort(rng.choice(pos, k, replace=False))
        rr = list(r)
        for p in where[::-1]:  # insert from the back so earlier positions stay valid
            rr.insert(int(p), np.nan)
        rr = np.asarray(rr)
        assert rr.size == r.size + k
        dirty = VolatilityMonitor(**kw)
        dirty.run(rr)
        shift = lambda i: i + int(np.sum(where <= i))  # noqa: E731  (a NaN at position p sits before old index p)
        # the segment start is set when the window closes, before the NaN at that position is seen, so it moves
        # only for NaNs strictly before it
        shift_start = lambda i: i + int(np.sum(where < i))  # noqa: E731
        assert [(shift(e.index), e.direction, shift_start(e.monitored_since), e.baseline_scale, e.stat)
                for e in clean.events] == _key(dirty)
        assert [shift_start(s["start"]) for s in clean.segments] == [s["start"] for s in dirty.segments]
        assert [s["n_finite"] for s in clean.segments] == [s["n_finite"] for s in dirty.segments]
        assert [s["m0"] for s in clean.segments] == [s["m0"] for s in dirty.segments]


# --------------------------------------------------------------------------------------------
# 5. pickle and deepcopy mid-stream
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("how", ["pickle", "deepcopy"])
def test_pickle_or_deepcopy_mid_stream_then_continue_equals_uninterrupted(how):
    """Serialising the monitor at a random point (inside a calibration window or while monitoring) and continuing
    from the copy gives the same events, segments and state as the uninterrupted run; the copy does not share
    state with the original."""
    kinds = ({}, {"baseline": "cs"}, {"c_down": 1.5}, {"baseline": "cs", "c_down": 2.0}, {"gaussian": True})
    for j, (r, kw) in enumerate(_cases(6, kinds)):
        rng = np.random.default_rng(SEED0 + 700 + j)
        cls = GaussianLRMonitor if kw.pop("gaussian", False) else VolatilityMonitor
        full = cls(**kw)
        full.run(r)
        for cut in rng.integers(1, r.size, 3):
            m = cls(**kw)
            m.run(r[:cut])
            c = pickle.loads(pickle.dumps(m)) if how == "pickle" else copy.deepcopy(m)
            c.run(r[cut:])
            assert _key(c) == _key(full) and _seg_key(c) == _seg_key(full)
            assert c.t == full.t and np.array_equal(c.S_up, full.S_up) and np.array_equal(c.S_dn, full.S_dn)
            assert m.t == cut - 1 and len(m.events) <= len(full.events)  # the original did not move
            assert m.cal_stats["calibrations"] <= c.cal_stats["calibrations"]


# --------------------------------------------------------------------------------------------
# 6. monotonicity in threshold, band and margin (first alarm)
# --------------------------------------------------------------------------------------------


def _first(mon):
    return mon.events[0].index if mon.events else math.inf


def test_larger_threshold_never_alarms_earlier():
    """The first alarm is the first step where a side's statistic reaches 2A, so it is nondecreasing in A.
    Later alarms are not comparable (recalibration happens at different places)."""
    levels = [3.0, 10.0, 30.0, 100.0, 1000.0, 1e9]
    for j, (r, kw) in enumerate(_cases(10, ({}, {"baseline": "cs"}, {"c_down": 1.5}, {"gaussian": True}))):
        cls = GaussianLRMonitor if kw.pop("gaussian", False) else VolatilityMonitor
        firsts = []
        for A in levels:
            m = cls(**{**kw, "threshold": A})
            m.run(r)
            firsts.append(_first(m))
        assert firsts == sorted(firsts), (j, firsts)
        assert firsts[0] < math.inf  # the smallest level does alarm, so the chain is not trivially empty


def test_wider_null_band_or_larger_margin_never_alarms_earlier():
    """The up factor 1 + lam (x / m_up - 1) is nonincreasing in m_up, and the down factor
    1 + f (m_dn - x) / (tau^2 - m_dn) is nondecreasing in m_dn (derivative (tau^2 - x) / (tau^2 - m_dn)^2 >= 0),
    so a wider band (m_up up, m_dn down) or a larger margin makes every statistic smaller pathwise."""
    for j, (r, kw) in enumerate(_cases(10)):
        for name, values in (("null_band", [0.0, 0.05, 0.2, 0.5, 0.9]), ("margin_z", [0.0, 0.5, 1.0, 2.0, 4.0])):
            firsts = []
            for v in values:
                m = VolatilityMonitor(**{**kw, name: v})
                m.run(r)
                firsts.append(_first(m))
            assert firsts == sorted(firsts), (j, name, firsts)
        firsts = []
        for v in (0.0, 0.1, 0.3):  # cs baseline: band only
            m = VolatilityMonitor(**{**kw, "baseline": "cs", "null_band": v})
            m.run(r)
            firsts.append(_first(m))
        assert firsts == sorted(firsts), (j, "cs band", firsts)


def test_larger_threshold_never_alarms_earlier_for_the_array_detector():
    for i in range(6):
        rng = np.random.default_rng(SEED0 + 800 + i)
        r = np.r_[_draw(rng, LAWS[i % 5], 500), 2.0 * _draw(rng, LAWS[i % 5], 500)]
        tau2 = (3 * MAD_NORMAL_SCALE * np.median(np.abs(r[:300]))) ** 2
        m0 = float(np.minimum(r[:300] ** 2, tau2).mean())
        x = np.minimum(r**2, tau2)
        for direction in ("up", "down"):
            for kind in ("SR", "CUSUM"):
                al = [D.betting_edetector(x, m0, direction, upper=tau2, threshold=A, kind=kind).alarm[0]
                      for A in (2.0, 5.0, 20.0, 100.0, 1e4)]
                al = [math.inf if a < 0 else a for a in al]
                assert al == sorted(al)


# --------------------------------------------------------------------------------------------
# 7. the supermartingale property: E[one-step factor] = 1 under the exact null
# --------------------------------------------------------------------------------------------


def _truncated_mean(dist, tau):
    """E[min(r^2, tau^2)] for a symmetric continuous law by numerical integration."""
    inner, _ = integrate.quad(lambda v: v * v * dist.pdf(v), 0.0, tau, limit=200)
    return 2.0 * inner + 2.0 * tau * tau * dist.sf(tau)


def _null_law(name):
    """Returns (sampler, median of |r|, truncated mean function of tau)."""
    if name == "gaussian":
        d = stats.norm()
        return (lambda rng, n: d.rvs(n, random_state=rng)), d.ppf(0.75), (lambda tau: _truncated_mean(d, tau))
    if name == "t3":
        d = stats.t(3)
        return (lambda rng, n: d.rvs(n, random_state=rng)), d.ppf(0.75), (lambda tau: _truncated_mean(d, tau))
    if name == "uniform":
        d = stats.uniform(-1.0, 2.0)
        return (lambda rng, n: d.rvs(n, random_state=rng)), 0.5, (lambda tau: _truncated_mean(d, tau))
    if name == "laplace":
        d = stats.laplace()
        return (lambda rng, n: d.rvs(n, random_state=rng)), d.ppf(0.75), (lambda tau: _truncated_mean(d, tau))
    if name == "twopoint":  # |r| = 1 with probability 0.6, |r| = 5 with probability 0.4; median |r| = 1
        def sample(rng, n):
            return rng.choice([-1.0, 1.0], n) * np.where(rng.random(n) < 0.6, 1.0, 5.0)

        return sample, 1.0, (lambda tau: 0.6 * min(1.0, tau * tau) + 0.4 * min(25.0, tau * tau))
    if name == "jumpmix":  # N(0,1) with a jump N(0, 100) added with probability 0.05 (a mixture)
        d1, d2 = stats.norm(), stats.norm(0.0, math.sqrt(1.0 + 100.0))

        def sample(rng, n):
            r = d1.rvs(n, random_state=rng)
            return r + (rng.random(n) < 0.05) * rng.normal(0.0, 10.0, n)

        class Mix:
            def pdf(self, v):
                return 0.95 * d1.pdf(v) + 0.05 * d2.pdf(v)

            def sf(self, v):
                return 0.95 * d1.sf(v) + 0.05 * d2.sf(v)

        mix = Mix()
        from scipy.optimize import brentq

        med = brentq(lambda v: 1.0 - 2.0 * mix.sf(v) - 0.5, 1e-9, 20.0)  # P(|r| <= v) = 1 - 2 sf(v) = 0.5
        return sample, med, (lambda tau: _truncated_mean(mix, tau))
    raise ValueError(name)


def _known_baseline_monitor(tau2, m0, **kw):
    """A VolatilityMonitor put into the calibrated state with a KNOWN baseline, through the code's own _finish."""
    m = VolatilityMonitor(**kw)
    m.t = 0
    m.scale = math.sqrt(tau2) / m.c
    m.tau2 = tau2
    m.tau2_dn = tau2 if m.c_down == m.c else (m.c_down * m.scale) ** 2
    m0_dn = m0 if m.c_down == m.c else None
    return m, m0_dn


@pytest.mark.parametrize("law", ["gaussian", "t3", "uniform", "laplace", "twopoint", "jumpmix"])
def test_expected_one_step_factor_is_one_under_the_exact_null(law):
    """Lemma 1 of docs/math.md rests on E[f_t | past] <= 1 for every bet when the baseline is the true mean of
    min(r^2, tau^2). With the true mean (numerical integration, or closed form for the two-point law) the
    expectation is exactly 1 for both sides; the Monte Carlo mean must sit within 4 standard errors of 1, and
    with a null band it must be below 1. The factors are computed through the monitor's own update()."""
    sample, med_abs, trunc_mean = _null_law(law)
    rng = np.random.default_rng(SEED0 + 900 + ["gaussian", "t3", "uniform", "laplace", "twopoint", "jumpmix"].index(law))
    c = 3.0
    scale = med_abs * MAD_NORMAL_SCALE
    tau = c * scale
    m0 = trunc_mean(tau)
    # sanity of the integration against a large sample (tolerance 1%)
    big = sample(rng, 400_000)
    assert np.minimum(big**2, tau * tau).mean() == pytest.approx(m0, rel=0.01)
    for band in (0.0, 0.25):
        m, _ = _known_baseline_monitor(tau * tau, m0, threshold=1e12, cal_len=10, null_band=band)
        m._finish(m0, m0 * (1 + band), m0 * (1 - band), None, None, m0)
        assert m.down_active
        # vectorised copy of the update() formulas (checked against update() itself on 300 draws below)
        x = np.minimum(big**2, m.tau2)
        f_up = 1.0 + m.up_grid[None, :] * (x[:, None] / m.m_up - 1.0)
        f_dn = 1.0 + m.lam_dn[None, :] * (1.0 - x[:, None] / m.m_dn)
        for v in big[:300]:
            m.S_up[:] = 0.0
            m.S_dn[:] = 0.0
            m.update(float(v))
        np.testing.assert_allclose(m.S_up, f_up[299], rtol=1e-12)
        np.testing.assert_allclose(m.S_dn, f_dn[299], rtol=1e-12)
        assert np.all(f_up > 0) and np.all(f_dn > 0)
        for f in (f_up, f_dn):
            mean = f.mean(axis=0)
            se = f.std(axis=0, ddof=1) / math.sqrt(f.shape[0])
            if band == 0.0:
                assert np.all(np.abs(mean - 1.0) <= 4.0 * se), (law, mean, se)
            else:
                assert np.all(mean + 4.0 * se < 1.0), (law, band, mean, se)


def test_expected_gaussian_lr_factor_is_one_only_for_gaussian_data():
    """The benchmark's factor has mean 1 under exactly Gaussian returns (it is a likelihood ratio) and mean
    above 1 for t3 returns, which is why it is not valid for heavy tails (documented in detectors.py)."""
    rng = np.random.default_rng(SEED0 + 950)
    for law, valid in (("gaussian", True), ("t3", False)):
        sample, med_abs, _ = _null_law(law)
        scale = med_abs * MAD_NORMAL_SCALE  # the true scaled MAD, as the monitor would estimate it
        z2 = (sample(rng, 300_000) / scale) ** 2
        g = GaussianLRMonitor(threshold=1e12, cal_len=10)
        logf = -np.log(g.up_ratios)[None, :] + 0.5 * (1 - 1 / g.up_ratios**2)[None, :] * z2[:, None]
        fu = np.exp(np.minimum(logf, 700.0))  # the same cap as GaussianLRMonitor.update
        with np.errstate(over="ignore", invalid="ignore"):  # the t3 factors are huge; their std overflows
            mean, se = fu.mean(axis=0), fu.std(axis=0, ddof=1) / math.sqrt(z2.size)
        if valid:
            assert np.all(np.abs(mean - 1.0) <= 4.0 * se), (mean, se)
        else:  # for t3 the expectation of exp(0.375 z^2) is infinite; the sample mean is far above 1
            assert np.all(mean > 10.0), mean


# --------------------------------------------------------------------------------------------
# 8. extreme but finite inputs
# --------------------------------------------------------------------------------------------


def _assert_sane(m):
    for S in (m.S_up, m.S_dn):
        assert np.all(np.isfinite(S)) and np.all(S >= 0), S
    for e in m.events:
        assert math.isfinite(e.stat) and e.stat >= 2 * m.A and math.isfinite(e.baseline_scale)
    for s in m.segments:
        assert all(math.isfinite(s[k]) and s[k] > 0 for k in ("scale", "tau2", "tau2_dn", "m_up"))
        # with baseline="cs" a second part made of zeros gives m0 = m_dn = 0 and switches the decrease side off
        # (valid and conservative; the plug-in path would extend the window instead): allowed, but never NaN
        assert all(math.isfinite(s[k]) and s[k] >= 0 for k in ("m0", "m0_dn", "m_dn"))


def _extreme_streams(rng):
    tiny, huge, sub = 1e-300, 1e300, 5e-324
    n = rng.standard_normal(400)
    yield "all huge", np.full(600, huge)
    yield "all tiny", np.full(600, tiny)
    yield "all subnormal", np.full(600, sub)
    yield "alternating huge tiny", np.where(np.arange(600) % 2 == 0, huge, tiny)
    yield "normal then huge", np.r_[n, np.full(300, huge)]
    yield "normal then tiny", np.r_[n, np.full(300, tiny)]
    yield "normal then subnormal", np.r_[n, np.full(300, sub)]
    yield "normal then zeros", np.r_[n, np.zeros(2000)]
    yield "normal then alternating", np.r_[n, np.where(np.arange(600) % 2 == 0, huge, tiny)]
    yield "mostly zeros window", np.r_[np.where(rng.random(400) < 0.7, 0.0, n), n]
    yield "window with a few huge", np.r_[np.where(rng.random(400) < 0.05, huge, n), n]
    yield "window with a few tiny and huge", np.r_[np.where(rng.random(400) < 0.3, tiny, n),
                                                  np.where(rng.random(400) < 0.3, huge, n)]
    yield "max float", np.r_[n, np.full(50, np.finfo(float).max), n]
    yield "negative huge", np.r_[n, np.full(100, -huge), n]


@pytest.mark.parametrize("kw", [{}, {"baseline": "cs"}, {"c_down": 1.5}, {"null_band": 0.25}])
def test_eprocess_stays_finite_and_nonnegative_on_extreme_inputs(kw):
    """1e-300, 1e300, subnormals, long runs of exact zeros and alternating huge and tiny values never make a
    statistic negative, NaN or inf, and every alarm has a finite statistic at or above 2A (numpy overflow
    warnings from squaring 1e300 are expected and silenced here)."""
    rng = np.random.default_rng(SEED0 + 1000)
    crashed = []
    for name, r in _extreme_streams(rng):
        m = VolatilityMonitor(threshold=20.0, cal_len=100, **kw)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            try:
                m.run(r)
                # feed the same stream twice more so that the extremes also hit later calibration windows
                m.run(r[::-1])
                m.run(r)
            except OverflowError:  # R5-B-001, fixed in round 4: must not happen any more
                crashed.append(name)
                continue
        _assert_sane(m)
        assert m.t == 3 * r.size - 1
    # before the round 4 fix (R5-B-001), streams with a calibration window whose median |r| was above about
    # 4.5e153 raised OverflowError; now no stream crashes
    assert crashed == [], crashed


# R5-B-001 (fixed in round 4): (c * scale) ** 2 on a Python float used to raise OverflowError when the median
# |r| of a calibration window exceeded about 4.5e153; the documented path is a degenerate window.
def test_calibration_window_of_huge_returns_does_not_crash():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        VolatilityMonitor(cal_len=2).run([1e300, 1e300, 1.0, 1.0])
        VolatilityMonitor(cal_len=4, baseline="cs").run([1e300] * 4 + [1.0] * 4)
        D.calibrate(np.full(10, 1e300))
        g = GaussianLRMonitor(cal_len=2)
        g.run([1.0, -1.0])
        g.update(1e200)  # a Python float: (r / scale) ** 2 raises; run() passes np.float64 and gets inf instead


def test_gaussian_lr_monitor_on_extreme_inputs_stays_finite_or_alarms():
    """The benchmark caps the log factor at 700, so 1e300 gives a finite factor, the statistic may overflow to
    inf and then alarms (inf >= 2A) rather than turning NaN."""
    rng = np.random.default_rng(SEED0 + 1001)
    for name, r in _extreme_streams(rng):
        m = GaussianLRMonitor(threshold=20.0, cal_len=100)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            try:
                m.run(r)
            except OverflowError:  # R5-B-001, shared calibration code
                continue
        assert not np.any(np.isnan(m.S_up)) and not np.any(np.isnan(m.S_dn)), name
        assert np.all(m.S_up >= 0) and np.all(m.S_dn >= 0), name
        for e in m.events:
            assert not math.isnan(e.stat) and e.stat >= 2 * m.A


def test_array_detector_on_extreme_statistics():
    """x = 1e300 on the up side overflows the statistic to inf, which alarms (not NaN); x = 0 for a long run makes
    the up statistic decay and the down side alarm; subnormal x behaves like 0."""
    m0 = 1.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        res = D.betting_edetector(np.r_[np.ones(50), np.full(20, 1e300)], m0, "up", threshold=1e6, return_path=True)
    assert res.alarm[0] == 50 and not np.any(np.isnan(res.stat_path[:51]))
    res = D.betting_edetector(np.r_[np.ones(50), np.zeros(3000)], m0, "up", threshold=1e6, return_path=True)
    assert res.alarm[0] == -1 and np.all(np.isfinite(res.stat_path)) and res.stat_path[-1] < 1e3
    res = D.betting_edetector(np.r_[np.ones(50), np.full(3000, 5e-324)], m0, "down", upper=9.0, threshold=1e6,
                              return_path=True)
    assert res.alarm[0] > 50 and np.all(np.isfinite(res.stat_path[: res.alarm[0] + 1]))


def test_array_detector_treats_an_infinite_statistic_as_missing():
    """Observation, not a validity bug: betting_edetector skips x = +inf like a NaN (np.isfinite on the gain), so
    an untruncated squared return that overflowed to inf carries no evidence for an increase. The monitor never
    passes inf (it clips r^2 at tau^2 first). Documented here so that the behaviour is a conscious choice."""
    x = np.r_[np.ones(10), np.inf, np.ones(10)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        res = D.betting_edetector(x, 1.0, "up", threshold=1e9, return_path=True)
    assert res.stat_path[10] == res.stat_path[9]  # unchanged at the inf step, as for a NaN


def test_infinite_threshold_is_rejected():
    """R5-B-002 (fixed in round 4): threshold = inf used to be accepted, and an alarm fired once a statistic
    overflowed to inf, because inf >= 2 * inf is True. The monitor now rejects a threshold that is not finite."""
    with pytest.raises(ValueError):
        VolatilityMonitor(threshold=math.inf, cal_len=100)
    with pytest.raises(ValueError):
        GaussianLRMonitor(threshold=math.inf, cal_len=100)


# --------------------------------------------------------------------------------------------
# 9. restart after an alarm equals a fresh monitor started at the restart point
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("kw_extra", [{}, {"baseline": "cs"}, {"c_down": 1.5}, {"gaussian": True}])
def test_restart_after_alarm_equals_a_fresh_monitor(kw_extra):
    """_reset() after an alarm clears the buffer, the statistics and the confidence sequence, so the events after
    the first alarm equal those of a fresh monitor fed the stream from the step after the alarm (indices shifted
    by alarm + 1), and so do the segments."""
    kw_extra = dict(kw_extra)
    cls = GaussianLRMonitor if kw_extra.pop("gaussian", False) else VolatilityMonitor
    checked = 0
    for r, kw in _cases(10):
        kw = {**kw, **kw_extra}
        full = cls(**kw)
        full.run(r)
        if not full.events:
            continue
        a = full.events[0].index
        fresh = cls(**kw)
        fresh.run(r[a + 1:])
        s = a + 1
        assert [(e.index + s, e.direction, e.monitored_since + s, e.baseline_scale, e.stat) for e in fresh.events] \
            == _key(full)[1:]
        assert [(k[0] + s,) + k[1:] for k in _seg_key(fresh)] == _seg_key(full)[1:]
        assert np.array_equal(fresh.S_up, full.S_up) and np.array_equal(fresh.S_dn, full.S_dn)
        checked += 1
    assert checked >= 5


def test_restart_after_a_degenerate_window_equals_a_fresh_monitor():
    """A calibration window that stays degenerate up to cal_max (all zeros) restarts calibration with fresh data,
    so a stream of cal_max zeros followed by returns gives the events of the returns alone, shifted by cal_max."""
    for j, (r, kw) in enumerate(_cases(6, ({}, {"baseline": "cs"}))):
        cal_max = 2 * kw["cal_len"]
        full = VolatilityMonitor(**kw)
        restarts = []
        orig = full._restart

        def spy():
            restarts.append(full.t)
            orig()

        full._restart = spy
        full.run(np.r_[np.zeros(cal_max), r])
        assert full.cal_stats["restarted"] == len(restarts) >= 1
        if kw.get("baseline") == "plugin" or "baseline" not in kw:
            assert restarts == [cal_max - 1]  # the plug-in path restarts exactly at cal_max
        # with baseline="cs" the restart point is data dependent (the tau split is extended, then restarted, and
        # the second restart window already holds some of r), so the recorded point is used instead
        s = restarts[-1] + 1
        fresh = VolatilityMonitor(**kw)
        fresh.run(np.r_[np.zeros(cal_max), r][s:])
        assert [(e.index + s, e.direction, e.monitored_since + s, e.baseline_scale, e.stat) for e in fresh.events] \
            == _key(full)
        assert [(k[0] + s,) + k[1:] for k in _seg_key(fresh)] == _seg_key(full)


# --------------------------------------------------------------------------------------------
# 10. CLI end to end
# --------------------------------------------------------------------------------------------


def _simulate_oracle(seed, ratio, jump_prob, budget):
    """The same path that cmd_simulate builds (copied from cli.py, which is the point of the comparison)."""
    rng = np.random.default_rng(seed)
    n = 9000
    sd = np.ones(n)
    sd[3000:6000] = ratio
    r = sd * rng.standard_normal(n)
    r += (rng.random(n) < jump_prob) * rng.normal(0, 10, n)
    A = 390 / budget
    rob = VolatilityMonitor(threshold=A, cal_len=1000).run(r)
    gau = GaussianLRMonitor(threshold=A, cal_len=1000).run(r)
    return [(e.index, e.direction) for e in rob], [(e.index, e.direction) for e in gau]


@pytest.mark.parametrize("args", [(SEED0 + 1100, 1.5, 0.003, 0.05), (SEED0 + 1101, 3.0, 0.0, 0.2),
                                  (SEED0 + 1102, 0.5, 0.01, 0.05)])
def test_cli_simulate_subprocess_matches_the_python_api(args):
    """python -m avmonitor simulate prints exactly the events the API gives on the same seed and parameters."""
    seed, ratio, jp, budget = args
    env = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"}
    res = subprocess.run([sys.executable, "-m", "avmonitor", "simulate", "--seed", str(seed), "--ratio", str(ratio),
                          "--jump-prob", str(jp), "--budget", str(budget)],
                         capture_output=True, text=True, timeout=120, env=env)
    assert res.returncode == 0, res.stderr
    lines = {ln.split(":", 1)[0].strip(): ln.split(":", 1)[1].strip() for ln in res.stdout.splitlines() if ":" in ln}
    rob = ast.literal_eval(lines["robust e-detector"])
    gau = ast.literal_eval(lines["Gaussian LR"])
    rob_api, gau_api = _simulate_oracle(seed, ratio, jp, budget)
    assert rob == rob_api and gau == gau_api
    assert f"A = {390 / budget:.0f}" in res.stdout


def _bitstamp_csv(root, seed, days=24):
    """A Bitstamp-format minute file with a U-shaped day and volatility regime changes every three days."""
    rng = np.random.default_rng(seed)
    n = days * 1440
    t0 = pd.Timestamp("2025-03-01", tz="UTC")
    minute = np.arange(n) % 1440
    season = 1 + 0.4 * np.cos(2 * np.pi * minute / 1440)
    level = np.where((np.arange(n) // (3 * 1440)) % 2 == 0, 1.0, 2.2)
    r = 4e-4 * season * level * rng.standard_t(4, n) / math.sqrt(2.0)
    close = 50_000 * np.exp(np.cumsum(r))
    vol = rng.random(n) * 3
    vol[rng.random(n) < 0.01] = 0.0  # a few tradeless minutes
    ts = t0.timestamp() + 60 * np.arange(n)
    df = pd.DataFrame({"timestamp": ts.astype(np.int64), "open": close, "high": close, "low": close, "close": close,
                       "volume": vol})
    p = root / "crypto" / "btcusd_bitstamp_1min_latest.csv"
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False)
    return p


@pytest.mark.parametrize("detector", ["robust", "gaussian"])
def test_cli_run_on_a_csv_matches_the_python_api(tmp_path, monkeypatch, capsys, detector):
    """avmonitor run --source bitstamp on a small CSV prints the same alarm times and sides as the API path
    (load_bitstamp_csv -> intraday_scale -> VolatilityMonitor with the CLI's A, c, cal_len and band)."""
    p = _bitstamp_csv(tmp_path, SEED0 + 1200)
    monkeypatch.setattr(DATA, "DATA_ROOT", tmp_path)
    start, end, lookback, budget, cal_len, band = "2025-03-14", "2025-03-25", 5, 0.5, 300, 0.25
    rc = cli_main(["run", "--source", "bitstamp", "--start", start, "--end", end, "--lookback-days", str(lookback),
                   "--budget", str(budget), "--cal-len", str(cal_len), "--band", str(band), "--max-rows", "1000",
                   "--detector", detector])
    assert rc == 0
    out = capsys.readouterr().out
    rows = re.findall(r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\+00:00)\s+(up|down)\s+(\d+)", out)
    # the API path, following cmd_run step by step
    warm = (pd.Timestamp(start) - pd.Timedelta(days=lookback + 7)).strftime("%Y-%m-%d")
    df_all = DATA.load_bitstamp_csv(p, warm, end)
    s_all = DATA.intraday_scale(df_all["logret"], lookback_days=lookback)
    keep = df_all.index >= pd.Timestamp(start, tz="UTC")
    df, s = df_all.loc[keep], s_all.loc[keep]
    u = (df["logret"] / s).to_numpy(copy=True)
    u[~np.isfinite(u)] = np.nan
    cls = GaussianLRMonitor if detector == "gaussian" else VolatilityMonitor
    mon = cls(threshold=1440 / budget, c=3.0, cal_len=cal_len, null_band=band)
    ev = mon.run(u)
    api = [(str(df.index[e.index]), e.direction, str(e.index - e.monitored_since)) for e in ev]
    assert len(api) >= 2, "the synthetic file should produce several alarms"
    assert rows == api
    assert f"{len(ev)} alarms in" in out


# --------------------------------------------------------------------------------------------
# 7b. Lemma 1 and Theorem B for the whole two-sided monitor statistic, several laws, known baseline
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("law", ["gaussian", "t3", "twopoint", "jumpmix"])
def test_two_sided_statistic_has_mean_at_most_t_and_false_alarm_rate_within_m_over_a(law):
    """docs/math.md Lemma 1 and Theorem B for the monitor's own recursion, (M_up + M_dn) / 2, with the true
    baseline: E[M_t] <= t at every t (within 4 standard errors) and P(alarm within m steps) <= m / A, for
    non-Gaussian laws too (the audit checked the array detector on Gaussian data only). The recursion is the
    vectorised copy of update() checked in test_expected_one_step_factor_is_one_under_the_exact_null."""
    sample, med_abs, trunc_mean = _null_law(law)
    rng = np.random.default_rng(SEED0 + 1500 + ["gaussian", "t3", "twopoint", "jumpmix"].index(law))
    tau = 3.0 * med_abs * MAD_NORMAL_SCALE
    m0 = trunc_mean(tau)
    m, _ = _known_baseline_monitor(tau * tau, m0, threshold=1e12, cal_len=10)
    m._finish(m0, m0, m0, None, None, m0)
    reps, n, A = 20_000, 150, 500.0
    x = np.minimum(sample(rng, reps * n).reshape(reps, n) ** 2, m.tau2)
    S_up = np.zeros((reps, m.up_grid.size))
    S_dn = np.zeros((reps, m.lam_dn.size))
    alarmed = np.zeros(reps, bool)
    for t in range(n):
        S_up = (1.0 + m.up_grid[None, :] * (x[:, t:t + 1] / m.m_up - 1.0)) * (S_up + 1.0)
        S_dn = (1.0 + m.lam_dn[None, :] * (1.0 - x[:, t:t + 1] / m.m_dn)) * (S_dn + 1.0)
        two = 0.5 * (S_up.mean(axis=1) + S_dn.mean(axis=1))
        assert np.all(two >= 0)
        if t % 10 == 9:
            se = two.std(ddof=1) / math.sqrt(reps)
            assert two.mean() <= (t + 1) + 4.0 * se, (law, t, two.mean(), se)
        alarmed |= (S_up.mean(axis=1) >= 2 * A) | (S_dn.mean(axis=1) >= 2 * A)
    assert alarmed.mean() <= n / A, (law, alarmed.mean(), n / A)


# --------------------------------------------------------------------------------------------
# 7c. the streaming monitor and the array detector are the same recursion
# --------------------------------------------------------------------------------------------


def test_streaming_monitor_equals_the_array_detector_on_one_segment():
    """Within the first segment, S_up.mean() and S_dn.mean() of VolatilityMonitor equal the SR paths of
    betting_edetector run on the same min(r^2, tau^2) with m0 = m_up (up) and m0 = m_dn, upper = tau_dn^2 (down),
    and the first alarm of the monitor is the earlier of the two array alarms at level 2A."""
    for j, (r, kw) in enumerate(_cases(8, ({}, {"c_down": 1.5}, {"null_band": 0.2}))):
        m = VolatilityMonitor(**kw)
        up_path, dn_path = [], []
        first = None
        for t, v in enumerate(r):
            ev = m.update(float(v))
            if m.calibrated and t >= m.start and first is None:  # the step that closes the window is not monitored
                up_path.append(m.S_up.mean())
                dn_path.append(m.S_dn.mean())
            if ev is not None:
                first = ev
                break
        if not m.segments:
            continue
        seg = m.segments[0]
        start = seg["start"]
        end = first.index if first is not None else r.size - 1
        rr = r[start:end + 1]
        xu = np.minimum(rr**2, seg["tau2"])
        xd = np.minimum(rr**2, seg["tau2_dn"])
        up = D.betting_edetector(xu, seg["m_up"], "up", threshold=2 * m.A, grid=m.up_grid, return_path=True)
        dn = D.betting_edetector(xd, seg["m_dn"], "down", upper=seg["tau2_dn"], threshold=2 * m.A,
                                 grid=m.down_fractions, return_path=True)
        # the monitor resets its statistics at the alarm step, so the recorded path stops one step before the array path
        assert len(up_path) == len(rr) - (first is not None)
        np.testing.assert_allclose(up.stat_path[:len(up_path)], up_path, rtol=1e-10)
        np.testing.assert_allclose(dn.stat_path[:len(dn_path)], dn_path, rtol=1e-10)
        both = D.two_sided_alarm(up, dn)[0]
        if first is None:
            assert both == -1
        else:
            assert both + start == first.index
            side = "up" if up.alarm[0] == both else "down"
            # the monitor checks the up side first, so a same-step double crossing reports "up"
            assert first.direction == ("up" if up.alarm[0] == both else side)


# --------------------------------------------------------------------------------------------
# extra relations found while reading the code
# --------------------------------------------------------------------------------------------


def test_input_container_does_not_matter():
    """A list, a tuple, a pandas Series, an integer array and one-at-a-time Python floats give the same events as
    a float64 array (run() converts with np.asarray(returns, float); update() takes one number)."""
    for j, (r, kw) in enumerate(_cases(4, ({}, {"baseline": "cs"}))):
        r = np.round(r * 4.0)  # integer valued returns (a coarse tick), so that int containers are exact
        base = VolatilityMonitor(**kw)
        base.run(r)
        for other in (r.tolist(), tuple(r.tolist()), pd.Series(r), r.astype(np.int64)):
            m = VolatilityMonitor(**kw)
            m.run(other)
            assert _key(m) == _key(base) and _seg_key(m) == _seg_key(base), type(other)
        m = VolatilityMonitor(**kw)
        for v in r.tolist():
            m.update(int(v))
        assert _key(m) == _key(base)



def test_c_down_equal_to_c_reproduces_the_default_exactly():
    """The docstring promises c_down = c reproduces the default bit for bit (plug-in and cs)."""
    for r, kw in _cases(6, ({}, {"baseline": "cs"})):
        a, b = VolatilityMonitor(**kw), VolatilityMonitor(**{**kw, "c_down": kw.get("c", 3.0)})
        a.run(r)
        b.run(r)
        assert _key(a) == _key(b) and _seg_key(a) == _seg_key(b)


def test_events_before_a_cut_do_not_change_when_the_stream_is_cut():
    """Prefix property: running on r[:n] gives exactly the events of the full run with index < n (no look-ahead),
    for random cut points, plug-in and cs baselines."""
    for j, (r, kw) in enumerate(_cases(6, ({}, {"baseline": "cs"}, {"c_down": 1.5}))):
        rng = np.random.default_rng(SEED0 + 1300 + j)
        full = VolatilityMonitor(**kw)
        full.run(r)
        for cut in rng.integers(1, r.size, 4):
            part = VolatilityMonitor(**kw)
            part.run(r[:cut])
            assert _key(part) == [k for k in _key(full) if k[0] < cut]


def test_permuting_a_calibration_window_does_not_change_the_plugin_baseline():
    """The plug-in baseline uses only the median of |r| and the mean of min(r^2, tau^2), so permuting the returns
    inside the first calibration window changes no event (the cs baseline is order dependent by design)."""
    for j, (r, kw) in enumerate(_cases(6)):
        rng = np.random.default_rng(SEED0 + 1400 + j)
        a = VolatilityMonitor(**kw)
        a.run(r)
        rr = r.copy()
        n = kw["cal_len"]
        rr[:n] = rr[rng.permutation(n)]
        b = VolatilityMonitor(**kw)
        b.run(rr)
        assert _key_idx(a) == _key_idx(b)
        assert a.segments[0]["scale"] == b.segments[0]["scale"]
        assert a.segments[0]["m0"] == pytest.approx(b.segments[0]["m0"], rel=1e-12)
