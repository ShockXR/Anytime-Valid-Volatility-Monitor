"""V1: independent recomputation of three R2 cells with the package's own streaming classes / detectors.
Own data generation (numpy directly, not change_returns), seeds 9_960_000 + cell. Plain loops, first alarm only.
Usage: python v1_r2_check.py <N_change> <N_null>"""
import sys, time
import numpy as np
from avmonitor.monitor import VolatilityMonitor, GaussianLRMonitor
from avmonitor import detectors as D
from avmonitor.robust import MAD_NORMAL_SCALE

NC, NN = int(sys.argv[1]), int(sys.argv[2])
CAL, NU, H = 1000, 400, 1950
H_README = 349.55965624105613      # R2 oracle, Gaussian, README
H_PAGE_T3 = 936.6199108655383      # R2 oracle, t3, Page CUSUM
H_GLR_G = 21232.618046337833       # R2 oracle on Gaussian, GLR (used on t3 in T4)


def gen(law, n, length, rng, change_at=None, f=1.0):
    if law == "gaussian":
        e = rng.standard_normal((n, length))
    else:
        e = rng.standard_t(3, (n, length)) / np.sqrt(3.0)
    if change_at is not None:
        e[:, change_at:] *= np.sqrt(f)
    return e


def first_event(make, path):
    m = make()
    for r in path:
        ev = m.update(float(r))
        if ev is not None:
            return ev.index - ev.monitored_since, ev.direction
    return -1, None


def summarise(name, al, nu=NU):
    al = np.asarray(al)
    early = (al >= 0) & (al < nu)
    det = al >= nu
    d = al[det] - nu + 1
    q = np.percentile(d, [25, 50, 75]) if d.size else [np.nan] * 3
    print(f"{name}: N={al.size} early={early.mean():.3f} detected={det.mean():.3f} "
          f"median={q[1]:.0f} IQR=({q[0]:.0f}-{q[2]:.0f})  bootSE(median)~{boot_se(d):.1f}", flush=True)


def boot_se(d, B=500):
    if d.size < 5:
        return np.nan
    rb = np.random.default_rng(1)
    return float(np.std(np.median(rb.choice(d, (B, d.size)), axis=1)))


def fa(name, al, horizon=390):
    al = np.asarray(al)
    k = int(((al >= 0) & (al < horizon)).sum())
    p = k / al.size
    print(f"{name}: false alarm within {horizon}: {k}/{al.size} = {p:.4f} (+-{1.96*np.sqrt(p*(1-p)/al.size):.4f})", flush=True)


t0 = time.time()
# ---- Cell A: README monitor, Gaussian, var x2, threshold 349.56 (R2: 37 (24-58), 95%; null FA 4.7%)
mk = lambda: VolatilityMonitor(threshold=H_README / 2, cal_len=CAL, null_band=0.25)
X = gen("gaussian", NN, CAL + 390, np.random.default_rng(9_960_000))
fa("A README gaussian null @349.6", [first_event(mk, x)[0] for x in X])
X = gen("gaussian", NC, CAL + NU + H, np.random.default_rng(9_960_001), CAL + NU, 2.0)
summarise("A README gaussian x2 @349.6", [first_event(mk, x)[0] for x in X])
print("time", round(time.time() - t0), flush=True)

# ---- Cell B: Page CUSUM up, t3, var x2, threshold 936.62 (R2: 203 (134-269), 94%; null FA 4.8%)
def page_alarm(X, h):
    cal, mon = X[:, :CAL], X[:, CAL:]
    out = np.empty(X.shape[0], int)
    for i in range(X.shape[0]):   # one path at a time, own scale
        s = float(np.median(np.abs(cal[i])) * MAD_NORMAL_SCALE)
        out[i] = D.page_cusum(mon[i] ** 2 / s**2, 1.0, 0.25, h).alarm[0]
    return out
X = gen("t3", 4 * NN, CAL + 390, np.random.default_rng(9_960_002))
fa("B Page t3 null @936.6", page_alarm(X, H_PAGE_T3))
X = gen("t3", 2 * NC, CAL + NU + H, np.random.default_rng(9_960_003), CAL + NU, 2.0)
summarise("B Page t3 x2 @936.6", page_alarm(X, H_PAGE_T3))
# fairness probe (not an R2 cell): CUSUM normalised by the calibration MEAN of r^2 instead of MAD scale^2,
# threshold tuned here to 5% within 390 steps on t3 null paths (seeds 9_960_006/7)
def page_mean_stat(X, steps):
    cal, mon = X[:, :CAL], X[:, CAL:CAL + steps]
    m0 = (cal**2).mean(axis=1)[:, None]
    y = mon**2 / m0
    G = np.zeros(X.shape[0]); W = np.empty_like(y)
    for t in range(y.shape[1]):
        G = np.maximum(0.0, G + y[:, t] - 1.25); W[:, t] = G
    return np.maximum.accumulate(W, axis=1)
Xn = gen("t3", 10000, CAL + 390, np.random.default_rng(9_960_006))
mx = page_mean_stat(Xn, 390)[:, -1]
h_mean = float(np.quantile(mx, 0.95))
Xc = gen("t3", 2 * NC, CAL + NU + H, np.random.default_rng(9_960_007), CAL + NU, 2.0)
cm = page_mean_stat(Xc, NU + H)
al = np.where((cm >= h_mean)[:, -1], (cm >= h_mean).argmax(axis=1), -1)
print(f"B' probe: Page CUSUM on r^2 / mean(cal r^2), t3, h={h_mean:.2f}", end=" ")
summarise("x2", al)
print("time", round(time.time() - t0), flush=True)

# ---- Cell C: Gaussian LR monitor tuned on Gaussian (21233), run on t3 (R2 T4: FA 99.9%, x2 0.1% detected)
mk = lambda: GaussianLRMonitor(threshold=H_GLR_G / 2, cal_len=CAL)
n_c = max(300, NN // 3)
X = gen("t3", n_c, CAL + 390, np.random.default_rng(9_960_004))
fa("C GLR(gauss-tuned) t3 null @21233", [first_event(mk, x)[0] for x in X])
X = gen("t3", n_c, CAL + NU + H, np.random.default_rng(9_960_005), CAL + NU, 2.0)
summarise("C GLR(gauss-tuned) t3 x2 @21233", [first_event(mk, x)[0] for x in X])
# same monitor on Gaussian null, to check the tuning point itself (R2: 4.9%)
X = gen("gaussian", NN, CAL + 390, np.random.default_rng(9_960_008))
fa("C0 GLR gaussian null @21233", [first_event(mk, x)[0] for x in X])
print("EXIT time", round(time.time() - t0), flush=True)
