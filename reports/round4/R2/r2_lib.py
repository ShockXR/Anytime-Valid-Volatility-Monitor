"""R2 helpers: vectorised running statistics for every method (SYNTHETIC data only).

Each function returns W, shape (N, T): the method's alarm statistic at each monitored step, computed WITHOUT
stopping. The method alarms at the first step with W >= h, so the first alarm for ANY threshold h is the first
passage of the running maximum of W over h. This lets one simulation serve every threshold (tuning, checking,
ROC). Large values are capped at 1e300, which changes nothing for thresholds below 1e300.

The recursions copy src/avmonitor/monitor.py (VolatilityMonitor, GaussianLRMonitor) and detectors.py
(betting_edetector, page_cusum, repeated_window_ztest) line by line; r2_validate.py checks that the first alarms
agree with the package's own code.
"""
import sys
import numpy as np
from avmonitor.detectors import DEFAULT_UP_GRID, DEFAULT_DOWN_FRACTIONS
from avmonitor.robust import MAD_NORMAL_SCALE
from avmonitor.confseq import cs_baseline
from avmonitor.simulate import change_returns

sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/stress")
from common import wilson  # noqa: E402

A = 7800.0          # 390 / 0.05
PROVED = 2 * A      # each side alarms at 2A, so W >= 2A is the package's alarm rule
CAL = 1000
CAP = 1e300
LAWS = {"gaussian": dict(), "student_t3": dict(heavy_tail_df=3.0), "jumps_p0.5pct_sd8": dict(p_jump=0.005, jump_scale=8.0)}
METHODS = ["README: plug-in, band 0.25", "plug-in, band 0", "confidence-sequence baseline", "known baseline (theorem)",
           "Gaussian LR monitor", "Page CUSUM up", "repeated 60-step z-test"]
E_DETECTORS = METHODS[:4]


def draw(law, N, n, rng, nu=None, f=1.0):
    """Returns of one law; variance multiplied by f from index nu on (the factor applies to the diffusive part)."""
    return change_returns(N, n, nu=nu, ratio=float(np.sqrt(f)), rng=rng, **LAWS[law]).r


def truth(law, seed, big=5_000_000):
    """tau = 3 x scaled MAD of the marginal and m0 = E min(r^2, tau^2), from one long sample (as in B1)."""
    r = draw(law, 1, big, np.random.default_rng(seed)).ravel()
    tau2 = (3.0 * np.median(np.abs(r)) * MAD_NORMAL_SCALE) ** 2
    return tau2, float(np.minimum(r**2, tau2).mean())


def sr_stat(mon, tau2, m_up, m_dn):
    """Two-sided betting e-detector of VolatilityMonitor: W_t = max(mean S_up, mean S_dn). Per-path tau2, m_up, m_dn."""
    N, T = mon.shape
    up = np.asarray(DEFAULT_UP_GRID, float)[None, :]
    fr = np.asarray(DEFAULT_DOWN_FRACTIONS, float)[None, :]
    tau2, m_up, m_dn = (np.broadcast_to(np.asarray(v, float), (N,)) for v in (tau2, m_up, m_dn))
    b = np.where(m_dn > 0, tau2 / np.where(m_dn > 0, m_dn, 1.0), np.inf)
    act = np.isfinite(b) & (b > 1.0)
    lam_dn = np.where(act[:, None], fr / (np.where(act, b, 2.0) - 1.0)[:, None], 0.0)
    x = np.minimum(mon**2, tau2[:, None])
    Su = np.zeros((N, up.shape[1]))
    Sd = np.zeros((N, fr.shape[1]))
    W = np.empty((N, T))
    with np.errstate(over="ignore", invalid="ignore"):
        for t in range(T):
            xt = x[:, t:t + 1]
            Su = np.minimum((1.0 + up * (xt / m_up[:, None] - 1.0)) * (Su + 1.0), CAP)
            Sd = np.where(act[:, None], np.minimum((1.0 + lam_dn * (1.0 - xt / m_dn[:, None])) * (Sd + 1.0), CAP), 0.0)
            W[:, t] = np.maximum(Su.mean(axis=1), Sd.mean(axis=1))
    return W


def glr_stat(mon, scale, up_ratios=(1.25, 1.5, 2.0), down_ratios=(0.5, 0.7, 0.8)):
    """GaussianLRMonitor: Shiryaev-Roberts mixture of Gaussian likelihood ratios, plug-in scale."""
    N, T = mon.shape
    ur, dr = np.asarray(up_ratios, float)[None, :], np.asarray(down_ratios, float)[None, :]
    z2 = (mon / scale[:, None]) ** 2
    Su = np.zeros((N, ur.shape[1]))
    Sd = np.zeros((N, dr.shape[1]))
    W = np.empty((N, T))
    with np.errstate(over="ignore"):
        for t in range(T):
            zt = z2[:, t:t + 1]
            fu = np.exp(np.minimum(-np.log(ur) + 0.5 * (1 - 1 / ur**2) * zt, 700))
            fd = np.exp(np.minimum(-np.log(dr) + 0.5 * (1 - 1 / dr**2) * zt, 700))
            Su = np.minimum(fu * (Su + 1.0), CAP)
            Sd = np.minimum(fd * (Sd + 1.0), CAP)
            W[:, t] = np.maximum(Su.mean(axis=1), Sd.mean(axis=1))
    return W


def page_stat(mon, scale, drift=0.25):
    """detectors.page_cusum on mon^2 / scale^2 with m0 = 1: G_t = max(0, G_{t-1} + y_t - 1 - drift)."""
    y = (mon**2 / scale[:, None] ** 2) / 1.0
    N, T = y.shape
    G = np.zeros(N)
    W = np.empty((N, T))
    for t in range(T):
        G = np.maximum(0.0, G + y[:, t] - 1.0 - drift)
        W[:, t] = G
    return W


def z_stat(mon, scale, window=60):
    """detectors.repeated_window_ztest on mon / scale (sigma0^2 = 1, var_r2 = 2): W_t = |z_t|, 0 before a full window.
    The test with level a alarms when W >= norm.ppf(1 - a / 2)."""
    r2 = (mon / scale[:, None]) ** 2
    N, T = r2.shape
    cs = np.cumsum(r2, axis=1)
    win = cs[:, window - 1:] - np.concatenate([np.zeros((N, 1)), cs[:, :-window]], axis=1)
    W = np.zeros((N, T))
    W[:, window - 1:] = np.abs((win - window * 1.0) / np.sqrt(window * 2.0))
    return W


def all_stats(data, tau2_true, m0_true):
    """Every method's statistic on paths data (N, CAL + T): calibration on the first CAL returns of each path."""
    cal, mon = data[:, :CAL], data[:, CAL:]
    scale = np.median(np.abs(cal), axis=1) * MAD_NORMAL_SCALE
    tau2 = (3.0 * scale) ** 2
    m0 = np.minimum(cal**2, tau2[:, None]).mean(axis=1)
    cs = cs_baseline(cal, c=3.0, delta_upper=0.005, delta_lower=0.005, tau_frac=0.25, n_plan=CAL - 250, cs_c=0.75)
    out = {}
    out[METHODS[0]] = sr_stat(mon, tau2, m0 * 1.25, m0 * 0.75)
    out[METHODS[1]] = sr_stat(mon, tau2, m0, m0)
    out[METHODS[2]] = sr_stat(mon, cs.tau**2, cs.upper, cs.lower)
    out[METHODS[3]] = sr_stat(mon, tau2_true, m0_true, m0_true)
    out[METHODS[4]] = glr_stat(mon, scale)
    out[METHODS[5]] = page_stat(mon, scale)
    out[METHODS[6]] = z_stat(mon, scale)
    return out


def first_alarm(cm, h):
    """First monitored index with running max >= h (-1 if none). cm = np.maximum.accumulate(W, axis=1)."""
    hit = cm >= h
    return np.where(hit[:, -1], hit.argmax(axis=1), -1)


def oracle_threshold(maxstat, target=0.05):
    """Threshold h with exactly round(target N) of the N tuning paths reaching it (midpoint of two order statistics)."""
    s = np.sort(maxstat)
    N = s.size
    k = int(round(target * N))
    return 0.5 * (s[N - k - 1] + s[N - k])
