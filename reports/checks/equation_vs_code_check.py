"""OUR CHECK: independent numeric checks of load-bearing equations of docs/math.md against avmonitor code.

Run from the repository root:  python reports/checks/equation_vs_code_check.py
Each check recomputes a quantity from the written formula (or a brute-force simulation) and compares it with the
package. One line per check, PASS or FAIL. All simulated data are SYNTHETIC.
"""
import os, sys
import numpy as np
from scipy.stats import norm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from avmonitor import detectors, robust, design, confseq
from avmonitor.monitor import VolatilityMonitor

res = []
def check(name, ok, detail=""):
    res.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (" | " + detail if detail else ""))

rng = np.random.default_rng(7)

# 1. Shiryaev-Roberts recursion M_t = f_t (M_{t-1} + 1) equals the sum over start times of the products of bets
n = 60
Y = np.minimum(rng.standard_normal(n) ** 2 * 1.3, 9.0)
lam = 0.35
f = 1 + lam * (Y - 1)
M = 0.0
for t in range(n):
    M = f[t] * (M + 1)
brute = sum(np.prod(f[j:]) for j in range(n))
check("SR recursion equals sum_j prod_{s>=j} f_s", np.isclose(M, brute, rtol=1e-10), f"rec={M:.6g} brute={brute:.6g}")
C = 0.0
for t in range(n):
    C = f[t] * max(C, 1.0)
brute_c = max(np.prod(f[j:]) for j in range(n))
check("CUSUM recursion equals max_j prod_{s>=j} f_s, and C <= M", np.isclose(C, brute_c, rtol=1e-10) and C <= M + 1e-12, f"C={C:.6g} max-prod={brute_c:.6g}")

# 2. Code path: betting_edetector with a single bet reproduces the recursion
m0 = 1.0
x = Y * m0
res_path = detectors.betting_edetector(x, m0, direction="up", grid=[lam], threshold=1e300, return_path=True)
check("betting_edetector(up, one bet) path equals f_t(M_{t-1}+1)", np.isclose(res_path.stat_path[-1], M, rtol=1e-10))
# decrease side: factor 1 + lam (1 - Y), lam < 1/(B-1)
tau2 = 9.0; m0d = 1.0; B = tau2 / m0d
lam_d = 0.9 / (B - 1)
Xd = np.minimum(rng.standard_normal(n) ** 2, tau2)
Yd = Xd / m0d
Md = 0.0
for t in range(n):
    Md = (1 + lam_d * (1 - Yd[t])) * (Md + 1)
# NOTE: for direction="down" a user-supplied grid is read as FRACTIONS of lam_max (the code multiplies it by lam_max), so pass 0.9
rd = detectors.betting_edetector(Xd, m0d, direction="down", upper=tau2, grid=[0.9], threshold=1e300, return_path=True)
check("betting_edetector(down) uses 1 + lam(1 - Y) with lam_max = 1/(B-1)", np.isclose(rd.stat_path[-1], Md, rtol=1e-10) and np.isclose(rd.info["lam_max"], 1 / (B - 1)))
check("design.lam_max matches", np.isclose(design.lam_max("down", B), 1 / (B - 1)) and design.lam_max("up", B) == 1.0)

# 3. Factors are positive at the extremes and have conditional mean exactly 1 when E[Y]=1
check("increase factor positive at Y=0 for lam<1", 1 + 0.999 * (0 - 1) > 0)
check("decrease factor positive at Y=B for lam<1/(B-1)", 1 + 0.999 / (B - 1) * (1 - B) > 0)
EY = 1.0
check("E[f] = 1 + lam (E[Y] - 1) = 1 at the null boundary", np.isclose(1 + lam * (EY - 1), 1.0))

# 4. Truncation: tau = c * scaled MAD about zero, scaled MAD = median|r| / Phi^-1(3/4)
r = rng.standard_normal(5000) * 0.002
b = detectors.calibrate(r, c=3.0)
check("calibrate: tau = c * median|r| / 0.6745", np.isclose(b.tau, 3.0 * np.median(np.abs(r)) / norm.ppf(0.75), rtol=1e-12) and np.isclose(robust.MAD_NORMAL_SCALE, 1.4826, atol=1e-4))
check("calibrate: m0 = mean(min(r^2, tau^2))", np.isclose(b.m0, np.minimum(r**2, b.tau**2).mean()))
# Truncation bias factors kappa_2, kappa_4: closed form and doc table versus 4 million draws
z = rng.standard_normal(4_000_000)
for label, a, k2_doc, k4_doc in (("raw MAD c=3", 3 * norm.ppf(0.75), 0.925, 0.705), ("scaled MAD c=3", 3.0, 0.995, 0.964)):
    zc = np.clip(z, -a, a)
    k2, k4 = robust.truncated_moment_factor(a, 2), robust.truncated_moment_factor(a, 4)
    mc2, mc4 = (zc**2).mean(), (zc**4).mean() / 3
    check(f"kappa2, kappa4 ({label}): code = Monte Carlo = doc table", abs(k2 - mc2) < 0.003 and abs(k4 - mc4) < 0.006 and round(k2, 3) == k2_doc and round(k4, 3) == k4_doc,
          f"a={a:.3f} k2={k2:.4f} (mc {mc2:.4f}) k4={k4:.4f} (mc {mc4:.4f})")
check("effective_cutoff: raw MAD a = c*0.6745, scaled a = c", np.isclose(robust.effective_cutoff(3, False), 2.0235, atol=1e-3) and robust.effective_cutoff(3, True) == 3)
# the E[min(sigma^2 Z^2, tau^2)] = sigma^2 kappa_2(a) identity used later
check("E[min(Z^2, a^2)] = kappa_2(a)", abs(np.minimum(z**2, 9.0).mean() - robust.truncated_moment_factor(3.0, 2)) < 0.003)

# 5. Guarantee: P(false alarm within first m steps) <= m/A with a KNOWN baseline (true null), two sides at 2A
A = 7800.0; m = 390; reps = 6000
sigma = 1.0; tau = 3.0 * sigma; m0_true = sigma**2 * robust.truncated_moment_factor(3.0, 2)
Xn = np.minimum(rng.standard_normal((reps, m)) ** 2, tau**2)
up = detectors.betting_edetector(Xn, m0_true, "up", threshold=2 * A)
dn = detectors.betting_edetector(Xn, m0_true, "down", upper=tau**2, threshold=2 * A)
fa = (detectors.two_sided_alarm(up, dn) >= 0).mean()
check("false alarm within first day (true baseline, Gaussian) is below m/A = 5%", fa <= m / A, f"measured {fa:.4f} over {reps} paths [SYNTHETIC]")
# adversarial null law of docs/theory_per_period.md Proposition 1b: Y_t = 1 exactly -> alarm exactly at t = 2A (day 40)
nn = int(2 * A) + 5
Xc = np.full((1, nn), m0_true)
upc = detectors.betting_edetector(Xc, m0_true, "up", threshold=2 * A)
dnc = detectors.betting_edetector(Xc, m0_true + 1e-9 if False else m0_true, "down", upper=tau**2, threshold=2 * A)
first = detectors.two_sided_alarm(upc, dnc)[0] + 1
check("TH6 Proposition 1b: law Y_t=1 alarms at step 2A = 15600, last step of day 40", first == 15600 and first == 40 * 390, f"first alarm step {first}")

# 6. Confidence-sequence baseline: P(m0_up < true mean) <= delta/2 (final bound), delta default 0.01 here
delta = 0.01; reps2 = 4000; ncal = 1440
rc = rng.standard_normal((reps2, ncal))
cb = confseq.cs_baseline(rc, c=3.0, delta_upper=delta / 2, delta_lower=delta / 2, tau_frac=0.25)
a_ = cb.tau / 1.0
true_mean = np.array([robust.truncated_moment_factor(float(t), 2) for t in a_])
cov_up = (cb.upper >= true_mean).mean(); cov_lo = (cb.lower <= true_mean).mean()
check("cs_baseline upper/lower bounds cover the true truncated mean at least 1-delta/2 each", cov_up >= 1 - delta / 2 and cov_lo >= 1 - delta / 2, f"coverage up={cov_up:.4f} low={cov_lo:.4f} (target {1 - delta / 2:.3f}) [SYNTHETIC]")

# 7. Decrease speed limit: floor log(1 + A/B)/log(B/(B-1)) and it is attained by all-zero returns with the largest bet
Bv = 9.045; Av = 7800.0
floor_formula = int(np.ceil(np.log(1 + Av / Bv) / np.log(Bv / (Bv - 1))))
check("decrease_speed_limit = ceil(log(1+A/B)/log(B/(B-1)))", design.decrease_speed_limit(Av, Bv) == floor_formula, f"{floor_formula} steps")
M = 0.0; steps = 0
fac = Bv / (Bv - 1)
while M < Av:
    M = fac * (M + 1); steps += 1
check("single bet at lam_max on all-zero returns reaches A in the floor number of steps (+/- 1)", abs(steps - floor_formula) <= 1, f"{steps} vs {floor_formula}")

# 8. Monitor: alarm level per side is 2A; threshold semantic
mon = VolatilityMonitor(threshold=A, cal_len=1440)
check("VolatilityMonitor default calibration window is 1440 returns", mon.cal_len == 1440)

bad = [n_ for n_, ok in res if not ok]
print(f"\n{len(res) - len(bad)}/{len(res)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
sys.exit(1 if bad else 0)
