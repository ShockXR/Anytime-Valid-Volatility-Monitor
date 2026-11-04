"""R3 shared code (SYNTHETIC only). New null laws and change shapes not covered by B1/B2.
Seeds: block 9_700_000 to 9_799_999 only. Law i: paths default_rng(9_700_000 + 1000 i), truth 9_750_000 + 1000 i.
Change cell j: default_rng(9_760_000 + 100 j)."""
import sys, time
import numpy as np, pandas as pd
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/stress")
from common import wilson
from avmonitor import detectors as D, robust as R, forecast as F
from avmonitor.monitor import VolatilityMonitor
from avmonitor.simulate import intraday_u_shape
from avmonitor.data import intraday_scale

A = 7800.0
CAL = 1000
HORIZ = [98, 390, 1950]
NMON = 1950
DAY = 390

# ----------------------------------------------------------------------------- laws
def t_unit(rng, df, shape):
    return rng.standard_t(df, shape) * np.sqrt((df - 2) / df)

def ma1(rng, theta, N, n):
    e = rng.standard_normal((N, n + 1))
    return (e[:, 1:] + theta * e[:, :-1]) / np.sqrt(1 + theta**2)

def tick_prices(rng, tick, N, n):
    """Random walk price with unit-variance steps, rounded to a grid of size ``tick``; returns = price differences."""
    p = 1000.0 + np.cumsum(rng.standard_normal((N, n + 1)), axis=1)
    q = np.round(p / tick) * tick
    return np.diff(q, axis=1)

TICK = {0.2: None, 0.4: None, 0.6: None}  # filled by calibrate_ticks

def calibrate_ticks():
    """Find the tick size (in units of the step sd) that gives the wanted share of exact-zero returns."""
    rng = np.random.default_rng(9_799_999)
    for target in TICK:
        lo, hi = 0.05, 5.0
        for _ in range(25):
            mid = 0.5 * (lo + hi)
            z = float(np.mean(tick_prices(rng, mid, 200, 2000) == 0))
            lo, hi = (mid, hi) if z < target else (lo, mid)
        TICK[target] = 0.5 * (lo + hi)
    return TICK

def two_piece(rng, s_left, s_right, shape):
    """Two-piece normal (Fechner): N(0, s_left) on the left, N(0, s_right) on the right, standardised to mean 0 var 1."""
    z = np.abs(rng.standard_normal(shape))
    side = rng.random(shape) < s_left / (s_left + s_right)
    x = np.where(side, -s_left * z, s_right * z)
    mean = np.sqrt(2 / np.pi) * (s_right - s_left)
    # exact: E[x^2] = (s_l^3 + s_r^3) / (s_l + s_r)
    ex2 = (s_left**3 + s_right**3) / (s_left + s_right)
    return (x - mean) / np.sqrt(ex2 - mean**2)

def with_missing_random(r, rng, p):
    r = r.copy()
    r[rng.random(r.shape) < p] = np.nan
    return r

def with_missing_blocks(r, rng, per_steps=390):
    """About one gap per 390 steps, each 30 to 120 steps long, at a random place (calibration and monitoring)."""
    r = r.copy()
    N, n = r.shape
    for i in range(N):
        k = rng.poisson(n / per_steps)
        for _ in range(k):
            L = int(rng.integers(30, 121))
            s = int(rng.integers(0, n - L))
            r[i, s:s + L] = np.nan
    return r

def seasonal_pipeline(rng, N, n_days=32, amplitude=0.6, lookback=20):
    """iid Gaussian x U-shape over n_days sessions of 390 minutes, deseasonalised with intraday_scale(start='past').
    Returns the deseasonalised series after the warm-up days (first 10 days are NaN by construction) as (N, m)."""
    prof = intraday_u_shape(DAY, amplitude)
    days = pd.date_range("2026-01-01", periods=n_days, freq="D")
    idx = pd.DatetimeIndex(np.concatenate([pd.date_range(d + pd.Timedelta(hours=9, minutes=15), periods=DAY, freq="min").to_numpy() for d in days]))
    out = []
    for i in range(N):
        r = np.tile(prof, n_days) * rng.standard_normal(n_days * DAY)
        s = intraday_scale(pd.Series(r, index=idx), lookback_days=lookback, smooth=15, start="past", circular=False)
        z = (r / s.to_numpy())
        out.append(z[np.isfinite(z)])
    m = min(len(z) for z in out)
    return np.stack([z[:m] for z in out])

# --- misspecified conditional-variance models for the forecast-relative pipeline (g)
def gjr(rng, N, n, omega, a, b, gam, burn=2000):
    h = np.full(N, omega / (1 - a - b - gam / 2)); out = np.empty((N, n))
    for t in range(n + burn):
        z = rng.standard_normal(N); r = np.sqrt(h) * z
        if t >= burn: out[:, t - burn] = r
        h = omega + (a + gam * (r < 0)) * r**2 + b * h
    return out

def egarch(rng, N, n, omega, a, gam, b, burn=2000):
    lh = np.full(N, omega / (1 - b)); out = np.empty((N, n))
    for t in range(n + burn):
        z = rng.standard_normal(N); r = np.exp(0.5 * lh) * z
        if t >= burn: out[:, t - burn] = r
        lh = omega + a * (np.abs(z) - np.sqrt(2 / np.pi)) + gam * z + b * lh
    return out

def lnsv(rng, N, n, phi, sig_eta, burn=2000):
    lh = np.zeros(N); out = np.empty((N, n))
    for t in range(n + burn):
        r = np.exp(0.5 * lh) * rng.standard_normal(N)
        if t >= burn: out[:, t - burn] = r
        lh = phi * lh + sig_eta * rng.standard_normal(N)
    return out

MISSPEC = {"gjr": lambda rng, N, n: gjr(rng, N, n, 0.02, 0.02, 0.92, 0.08),
           "egarch": lambda rng, N, n: egarch(rng, N, n, 0.0, 0.10, -0.05, 0.97),
           "lnsv": lambda rng, N, n: lnsv(rng, N, n, 0.98, 0.15)}

def forecast_pipeline(law, rng, N, n, fit):
    r = MISSPEC[law](rng, N, n + 2000)
    om, a, b = fit
    z = np.stack([F.standardise(r[i], F.garch11_filter(r[i], om, a, b)) for i in range(N)])
    return z[:, 2000:]

def fit_for(law, seed):
    rng = np.random.default_rng(seed)
    train = MISSPEC[law](rng, 1, 30000)[0]
    return F.garch11_fit(train)

LAWS = ["t30_to_t3", "t3_to_t30", "ma1_-0.3", "ma1_-0.6", "tick_zero20", "tick_zero40", "tick_zero60",
        "twopiece_0.6_1.4", "twopiece_0.4_1.6", "missing_5pct", "missing_blocks", "seasonal_pipeline",
        "fc_gjr", "fc_egarch", "fc_lnsv",
        "t30_to_t5", "t5_to_t30", "ma1_-0.1", "ma1_0.3", "seasonal_amp1.0"]  # extra laws, indices 15 to 19

def draw(law, rng, N, n, fit=None):
    """Returns (N, n) with n = CAL + NMON. The first CAL columns are calibration."""
    if "_to_t" in law and law.startswith("t"):
        a, b = (int(v) for v in law[1:].split("_to_t"))
        half = CAL + NMON // 2
        out = np.empty((N, n))
        out[:, :half] = t_unit(rng, a, (N, half))
        out[:, half:] = t_unit(rng, b, (N, n - half))
        return out
    if law.startswith("ma1_"):
        return ma1(rng, float(law[4:]), N, n)
    if law.startswith("tick_zero"):
        tgt = int(law[9:]) / 100
        if TICK[tgt] is None: calibrate_ticks()
        return tick_prices(rng, TICK[tgt], N, n)
    if law.startswith("twopiece_"):
        _, sl, sr = law.split("_")
        return two_piece(rng, float(sl), float(sr), (N, n))
    if law == "missing_5pct":
        return with_missing_random(rng.standard_normal((N, n)), rng, 0.05)
    if law == "missing_blocks":
        return with_missing_blocks(rng.standard_normal((N, n)), rng)
    if law == "seasonal_pipeline":
        z = seasonal_pipeline(rng, N)
        return z[:, :n]
    if law.startswith("seasonal_amp"):
        z = seasonal_pipeline(rng, N, amplitude=float(law[12:]))
        return z[:, :n]
    if law.startswith("fc_"):
        return forecast_pipeline(law[3:], rng, N, n, fit)
    raise ValueError(law)

# ----------------------------------------------------------------------------- detectors
def plugin_stats(cal):
    """Per-path scale, tau2, m0 exactly as VolatilityMonitor._plugin_stats (NaN dropped)."""
    N = cal.shape[0]
    scale = np.empty(N); m0 = np.empty(N)
    for i in range(N):
        r = cal[i][np.isfinite(cal[i])]
        scale[i] = np.median(np.abs(r)) * R.MAD_NORMAL_SCALE
        m0[i] = np.mean(np.minimum(r**2, (3 * scale[i])**2))
    return scale, (3 * scale)**2, m0

def vec_plugin(mon, tau2, m0, band, A=A):
    """Vectorised copy of VolatilityMonitor.update (plug-in baseline, c_down = c). Returns first alarm step (-1 none)
    and direction code (+1 up, -1 down, 0 none). NaN steps are skipped as in the monitor."""
    N, n = mon.shape
    up_grid = np.asarray(D.DEFAULT_UP_GRID, float); dfr = np.asarray(D.DEFAULT_DOWN_FRACTIONS, float)
    m_up = m0 * (1 + band); m_dn = m0 * (1 - band)
    b = tau2 / m_dn
    active = np.isfinite(b) & (b > 1)
    lam_dn = np.where(active[:, None], dfr[None, :] / np.where(active, b - 1, 1.0)[:, None], 0.0)
    S_up = np.zeros((N, len(up_grid))); S_dn = np.zeros((N, len(dfr)))
    alarm = np.full(N, -1); direc = np.zeros(N, int)
    live = np.ones(N, bool)
    for t in range(n):
        r = mon[:, t]
        ok = np.isfinite(r) & live
        x = np.minimum(r * r, tau2)
        nu = (1 + up_grid[None, :] * (x / m_up - 1)[:, None]) * (S_up + 1)
        nd = (1 + lam_dn * (1 - x / m_dn)[:, None]) * (S_dn + 1)
        S_up = np.where(ok[:, None], nu, S_up); S_dn = np.where(ok[:, None], nd, S_dn)
        mu = S_up.mean(1); md = S_dn.mean(1)
        hit_up = ok & (mu >= 2 * A); hit_dn = ok & ~hit_up & (md >= 2 * A)
        hit = hit_up | hit_dn
        alarm[hit] = t; direc[hit_up] = 1; direc[hit_dn] = -1
        live &= ~hit
    return alarm, direc

def truth(sample):
    scale = np.median(np.abs(sample)) * R.MAD_NORMAL_SCALE
    if not scale > 0:  # more than half the returns are exactly zero: same fallback as the monitor (median of the nonzero)
        scale = np.median(np.abs(sample[sample != 0])) * R.MAD_NORMAL_SCALE
    tau = 3 * scale
    return tau, float(np.minimum(sample**2, tau**2).mean())

def known_baseline(mon, tau, m0):
    x = np.minimum(mon**2, tau**2)
    up = D.betting_edetector(x, m0, "up", threshold=2 * A)
    dn = D.betting_edetector(x, m0, "down", upper=tau**2, threshold=2 * A)
    a = D.two_sided_alarm(up, dn)
    au = np.asarray(up.alarm); ad = np.asarray(dn.alarm)
    direc = np.where(a < 0, 0, np.where((au >= 0) & ((ad < 0) | (au <= ad)), 1, -1))
    return a, direc

def stream(make, data):
    out = np.full(data.shape[0], -1); direc = np.zeros(data.shape[0], int); seg = []
    for i in range(data.shape[0]):
        m = make()
        for r in data[i]:
            ev = m.update(r)
            if ev is not None:
                out[i] = ev.index - ev.monitored_since; direc[i] = 1 if ev.direction == "up" else -1
                break
        seg.append(m.segments[0] if m.segments else {})
    return out, direc, seg

def rows_for(law, method, bound_fn, alarm, direc, N, note, horiz=HORIZ):
    rows = []
    for m in horiz:
        hit = (alarm >= 0) & (alarm < m)
        k = int(hit.sum()); lo, hi = wilson(k, N); b = bound_fn(m)
        verdict = "n/a" if b is None else ("within bound" if hi <= b else ("EXCEEDS bound" if lo > b else "inconclusive"))
        rows.append(dict(law=law, method=method, horizon_m=m, k=k, paths=N, rate=k / N, ci_low=lo, ci_high=hi, bound=b,
                         k_up=int((hit & (direc > 0)).sum()), k_down=int((hit & (direc < 0)).sum()), verdict=verdict, note=note))
    return rows

METHODS = ["known", "plugin025", "plugin0", "cs"]

def run_methods(law, data, tau_m0, N, do_cs=True, horiz=HORIZ, extra_note="", use_stream=False):
    cal, mon = data[:, :CAL], data[:, CAL:]
    rows = []
    if tau_m0 is not None:
        tau, m0 = tau_m0
        a, d = known_baseline(mon, tau, m0)
        rows += rows_for(law, "e-detector, KNOWN baseline (theorem)", lambda m: m / A, a, d, N, f"tau={tau:.4f} m0={m0:.5f} {extra_note}", horiz)
    scale, tau2, m0p = plugin_stats(cal)
    for band, name in ((0.25, "monitor plugin band0.25 (README)"), (0.0, "monitor plugin band0")):
        if use_stream:  # the monitor may extend a degenerate window (tick grids): run the real stream
            a, d, seg = stream(lambda: VolatilityMonitor(threshold=A, cal_len=CAL, null_band=band), data)
            ext = np.mean([s.get("extended", False) for s in seg]); fb = np.mean([s.get("fallback", False) for s in seg])
            rows += rows_for(law, name, lambda m: m / A, a, d, N, f"plug-in stream; window extended on {ext:.2f}, fallback scale on {fb:.2f} of paths {extra_note}", horiz)
        else:
            a, d = vec_plugin(mon, tau2, m0p, band)
            rows += rows_for(law, name, lambda m: m / A, a, d, N, f"plug-in; median m0hat={np.median(m0p):.4f} {extra_note}", horiz)
    if do_cs:
        a, d, seg = stream(lambda: VolatilityMonitor(threshold=A, cal_len=CAL, baseline="cs", delta=0.01), data)
        rows += rows_for(law, "monitor cs baseline delta0.01", lambda m: 0.01 + m / A, a, d, N, f"cs {extra_note}", horiz)
    return rows
