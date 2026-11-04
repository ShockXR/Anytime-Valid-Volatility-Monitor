"""Shared helpers for the Agent 4 stress scripts (SYNTHETIC data only). Seeds: block 9_100_000 to 9_999_999,
which no file in experiments/ or src/ uses (largest seed constant found there: 26060 and the date-like 20261103/20261111)."""
import numpy as np

def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return c - h, c + h

def garch(rng, N, n, alpha, beta, burn=500, df=None):
    """Unit-variance GARCH(1,1) returns, shape (N, n). Innovations Gaussian (or unit-variance t if df)."""
    omega = 1.0 - alpha - beta
    h = np.ones(N)
    out = np.empty((N, n))
    for t in range(n + burn):
        z = rng.standard_normal(N) if df is None else rng.standard_t(df, N) * np.sqrt((df - 2) / df)
        r = np.sqrt(h) * z
        if t >= burn:
            out[:, t - burn] = r
        h = omega + alpha * r**2 + beta * h
    return out

def seasonal_profile(period=390, amplitude=0.6):
    from avmonitor.simulate import intraday_u_shape
    return intraday_u_shape(period, amplitude)

LAWS = ["gaussian", "student_t3", "student_t5", "jumps_p1pct_sd8", "pareto_alpha1.5", "garch_0.05_0.90",
        "garch_0.08_0.90", "unremoved_seasonality"]

def draw(law, rng, N, n):
    if law == "gaussian":
        return rng.standard_normal((N, n))
    if law == "student_t3":
        return rng.standard_t(3, (N, n)) / np.sqrt(3.0)
    if law == "student_t5":
        return rng.standard_t(5, (N, n)) * np.sqrt(3.0 / 5.0)
    if law == "jumps_p1pct_sd8":
        return rng.standard_normal((N, n)) + (rng.random((N, n)) < 0.01) * rng.normal(0, 8.0, (N, n))
    if law == "pareto_alpha1.5":  # symmetric Pareto, tail index 1.5: infinite variance, finite MAD
        u = rng.random((N, n))
        return np.where(rng.random((N, n)) < 0.5, 1.0, -1.0) * (1.0 - u) ** (-1.0 / 1.5)
    if law == "garch_0.05_0.90":
        return garch(rng, N, n, 0.05, 0.90)
    if law == "garch_0.08_0.90":
        return garch(rng, N, n, 0.08, 0.90)
    if law == "unremoved_seasonality":  # iid Gaussian times a U-shaped intraday profile, random phase, NOT removed
        prof = seasonal_profile()
        ph = rng.integers(0, 390, N)
        idx = (ph[:, None] + np.arange(n)[None, :]) % 390
        return prof[idx] * rng.standard_normal((N, n))
    raise ValueError(law)
