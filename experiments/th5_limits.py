"""TH5. Limits for every valid procedure, the plug-in baseline when n, m and A grow together, and the
decrease bet near its boundary.

All results are SYNTHETIC. Seeds are fixed (listed below). Theory: docs/theory_limits.md. Functions:
src/avmonitor/limits.py. True baseline unless stated: tau = 3 x the population scaled MAD, m0 = E[min(r^2,
tau^2)], c = 3, A = 7800 (one 390-step day has m / A = 5%).

Parts
  collapse     Closed-form floors for the most extreme paths (complete collapse, Y = 0 at every step, and
               complete explosion, Y = B at every step) for c = 1.5 to 5, at A and 2A: every deterministic valid
               procedure (floor lemma), the smallest mean delay of a randomised valid procedure on the path (a
               linear programme) and its largest alarm probability one step before the floor, every Shiryaev-Roberts e-detector with linear bets (Corollary 4 of docs/theory_delay.md), the default
               grid (th5_collapse.csv). No simulation.
  floors       18 post-change laws (Gaussian, Student t3, jumps; sd ratios 0.5, 0.67, 0.8, 1.25, 1.5, 2), change
               at the first monitored step, 10,000 paths of 4,000 steps each: the pathwise floor of every
               deterministic valid rule at A and 2A, the Neyman-Pearson bound at the least favourable null law
               (any rule, deterministic rules, randomised rules), the delays of the default grid and of the
               growth-optimal single bet, and a pathwise check that no simulated e-detector alarm comes before
               the floor (th5_floors.csv, th5_floor_cdf.csv).
  asymptotic   Six of these laws at A = 10^3, 10^4, 10^6, 10^8 (4,000 paths): floor and Neyman-Pearson bound
               divided by log A / KL_inf and by log A / gamma (th5_asymptotic.csv).
  boundary     The decrease bet lam_max (1 - eps): atom q_B, growth G_c without the atom, eps*, KL_inf, the bounds
               of Proposition 12, the Wald and block delay bounds and simulated delays for eps from 0.5 to 1e-6
               (nine laws, 4,000 paths of 6,000 steps) (th5_boundary.csv).
  plugin_bound Calibration windows only (no monitoring): for n = 390 to 92,160 and m = 98 to 7800, Corollary 1 of
               TH4 averaged over windows, the first-order formula of Theorem 8 (c), the explicit bound of
               Theorem 8 (a), the distribution-free bound of Theorem 8 (b), the noise-free stream and the hold-then-bold
               lower bound (Proposition 9), the inflated plug-in bound, the confidence-sequence baseline's
               effective z, and the share of windows whose plug-in has a positive growth rate under the i.i.d.
               law (Proposition A.6) (th5_plugin_bound.csv, th5_plugin_drift.csv); also runs part dfree.
  dfree        The distribution-free bound of Theorem 8 (b) for n = 390 to 10^8 and m = 5 to 1950, and the
               rate constant C_DF (th5_distribution_free.csv). No simulation.
  plugin_sim   Monitoring simulation: Gaussian and Student t3 returns, n = 390 to 23,040, 10,000 paths per cell,
               each with its own calibration window, 7800 monitored steps of the same law; the default
               two-sided detector with the plug-in baseline, the inflated plug-in (z = 1, 2, 2.576), the
               confidence-sequence baseline (delta = 0.01) and the true baseline on the same monitored returns
               (th5_plugin_sim.csv).
  holdbold     Simulation check of Proposition 9: the hold-then-bold stream run on fresh windows (Gaussian,
               n = 390 and 2880, m = 390 and 1950, 4,000 windows) against the per-window proved lower bound
               (th5_holdbold.csv).
  figures      results/figures/th5_floors.png and th5_plugin_regime.png.

Usage (each part writes its own CSV files; run with at most 2 processes on a shared machine):
  python experiments/th5_limits.py collapse floors asymptotic boundary plugin_bound plugin_sim holdbold figures
  TH5_SCALE=0.05 TH5_OUT=/tmp/th5 python experiments/th5_limits.py ...   # smoke run

Seeds (SEED = 19500): floors SEED + 1..18; asymptotic SEED + 100..105; boundary SEED + 200..208; plugin_bound
SEED + 300..315; plugin_sim SEED + 400..411 (windows and monitored returns); holdbold SEED + 500..503.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from avmonitor import design as D
from avmonitor import limits as L
from avmonitor import robustness as R
from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("TH5_OUT", ROOT / "results"))
FIG = OUT / "figures"
SCALE = float(os.environ.get("TH5_SCALE", "1"))
SEED = 19500
A = 7800.0
C = 3.0
LAM_UP = np.asarray(DEFAULT_UP_GRID, float)
FRAC_DN = np.asarray(DEFAULT_DOWN_FRACTIONS, float)
MODELS = {"gaussian": D.ReturnModel("gaussian"), "student t3": D.ReturnModel("student_t", df=3.0),
          "jumps": D.ReturnModel("jumps")}
RATIOS = (0.5, 0.67, 0.8, 1.25, 1.5, 2.0)


def _n(P):
    return int(max(200, round(P * SCALE))) if SCALE != 1 else int(P)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(c - h, 0.0), min(c + h, 1.0)


def _setup(name, rho):
    model = MODELS[name]
    tau, m0 = D.baseline(model, C)
    B = tau**2 / m0
    side = "up" if rho > 1 else "down"
    law = D.y_law(model, rho, tau, m0)
    lam, kl = L.least_favourable_bet(law, side)
    return model, tau, m0, B, side, law, lam, kl


def _grid(side, B):
    return (LAM_UP, np.full(6, 1 / 6)) if side == "up" else (FRAC_DN / (B - 1), np.full(5, 1 / 5))


def _gamma(law, B, side):
    """E_Q[log(1 / max_null_mass(Y))]: the growth rate of log W_n."""
    return float(np.sum(law.w * L.log_floor_factor(law.y, B, side)))


# ----------------------------------------------------------------------------
# Part collapse
# ----------------------------------------------------------------------------


def part_collapse():
    rows = []
    for c in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0):
        tau, m0 = D.baseline(MODELS["gaussian"], c)
        B = tau**2 / m0
        for level, AA in (("A", A), ("2A", 2 * A)):
            for side in ("down", "up"):
                F = B / (B - 1) if side == "down" else B
                bets, w = _grid(side, B)
                fac = 1 + bets * (1.0 if side == "down" else (B - 1))
                n_any = L.collapse_floor(AA, B, side)
                from scipy import optimize

                n_frac = optimize.brentq(lambda n: n * np.log(F) - np.log(AA / n), 1e-6, 1e6)
                rows.append(dict(c=c, B=B, side=side, level=level, A=AA, factor_per_step=F, log_factor=np.log(F),
                                 floor_any_procedure=n_any, floor_any_procedure_fractional=n_frac,
                                 floor_sr_linear_bets=D.sr_steps_constant_factor(AA, F),
                                 floor_sr_linear_bets_fractional=np.log1p(AA * (F - 1) / F) / np.log(F),
                                 floor_default_grid=D.speed_limit(AA, fac, w),
                                 first_order_log_A_over_log_F=np.log(AA) / np.log(F),
                                 rand_max_prob_by_floor_minus_1=min(1.0, (n_any - 1) * F ** (n_any - 1) / AA),
                                 rand_min_mean_delay=L.randomised_path_delay(np.arange(1, 4 * n_any + 1) * np.log(F), AA)[0]))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th5_collapse.csv", index=False)
    print(df[(df.side == "down")][["c", "B", "level", "floor_any_procedure", "floor_any_procedure_fractional",
                                   "rand_min_mean_delay", "rand_max_prob_by_floor_minus_1",
                                   "floor_sr_linear_bets_fractional", "floor_default_grid"]].to_string(index=False))
    print(df[(df.side == "up")][["c", "B", "level", "floor_any_procedure", "rand_min_mean_delay",
                                 "rand_max_prob_by_floor_minus_1", "floor_sr_linear_bets",
                                 "floor_default_grid"]].to_string(index=False))


# ----------------------------------------------------------------------------
# Part floors
# ----------------------------------------------------------------------------


def _summ(x, pref):
    x = np.asarray(x, float)
    return {f"{pref}_mean": float(np.mean(x)), f"{pref}_median": float(np.median(x)),
            f"{pref}_q10": float(np.quantile(x, 0.1)), f"{pref}_q90": float(np.quantile(x, 0.9))}


def floors_cell(name, rho, N, H, seed, levels=(A, 2 * A), keep_cdf=False):
    model, tau, m0, B, side, law, lam, kl = _setup(name, rho)
    rng = np.random.default_rng(seed)
    y = D.sample_y(model, rho, tau, m0, (N, H), rng)
    lw = L.log_floor_factor(y, B, side)
    u = (y - 1.0) if side == "up" else (1.0 - y)
    S_star = np.cumsum(np.log1p(lam * u), axis=1)
    bets, w = _grid(side, B)
    grid = L.sr_alarm_times(y, bets, w, side, levels)
    orac = L.sr_alarm_times(y, [lam], [1.0], side, levels)
    pstar_mass = float(np.sum(law.w / (1 + lam * ((law.y - 1) if side == "up" else (1 - law.y)))))
    row = dict(model=name, ratio=rho, side=side, B=B, lam_star=lam, lam_star_over_lam_max=lam * (1 if side == "up" else B - 1),
               kl_inf=kl, gamma=_gamma(law, B, side), q_B=float(law.w[-1]), p_star_mass=pstar_mass, paths=N, horizon=H)
    cdfs = []
    for i, AA in enumerate(levels):
        tag = "A" if AA == A else "2A"
        nf = L.floor_time(lw, AA)
        miss = int(np.sum(nf < 0))
        nf = np.where(nf < 0, H, nf)
        g = grid[i]
        o = orac[i]
        row.update({f"floor_missing_{tag}": miss, f"grid_missing_{tag}": int(np.sum(g < 0)),
                    f"oracle_missing_{tag}": int(np.sum(o < 0))})
        row.update(_summ(nf, f"floor_{tag}"))
        row.update(_summ(np.where(g < 0, H, g), f"grid_{tag}"))
        row.update(_summ(np.where(o < 0, H, o), f"oracle_{tag}"))
        both = (g > 0)
        row[f"pathwise_violations_grid_{tag}"] = int(np.sum(g[both] < nf[both]))
        row[f"pathwise_violations_oracle_{tag}"] = int(np.sum((o > 0) & (o < nf)))
        b_np = L.np_cdf(S_star, AA)
        b_det = L.np_cdf(S_star, AA, lw, "det")
        b_rand = L.np_cdf(S_star, AA, lw, "rand")
        rf = L.randomised_floor_cdf(lw, AA)
        row.update({f"np_any_{tag}": L.delay_lower_bound(b_np), f"np_det_{tag}": L.delay_lower_bound(b_det),
                    f"np_rand_{tag}": L.delay_lower_bound(b_rand), f"floor_rand_{tag}": L.delay_lower_bound(rf),
                    f"np_det_median_{tag}": L.quantile_lower_bound(b_det, 0.5),
                    f"np_any_median_{tag}": L.quantile_lower_bound(b_np, 0.5)})
        lb_det = max(row[f"floor_{tag}_mean"], row[f"np_det_{tag}"])
        row[f"lower_bound_det_{tag}"] = lb_det
        row[f"efficiency_grid_{tag}"] = row[f"grid_{tag}_mean"] / lb_det
        row[f"efficiency_oracle_{tag}"] = row[f"oracle_{tag}_mean"] / lb_det
        lb_rand = max(row[f"floor_rand_{tag}"], row[f"np_rand_{tag}"])
        row[f"lower_bound_rand_{tag}"] = lb_rand
        row[f"efficiency_grid_rand_{tag}"] = row[f"grid_{tag}_mean"] / lb_rand
        row[f"efficiency_oracle_rand_{tag}"] = row[f"oracle_{tag}_mean"] / lb_rand
        if AA == A:
            # batch-means Monte Carlo standard errors (10 batches of the same paths; no new random numbers)
            gm, om = np.where(g < 0, H, g), np.where(o < 0, H, o)
            est = []
            for ii in np.array_split(np.arange(N), 10):
                lb_b = max(float(np.mean(nf[ii])), L.delay_lower_bound(L.np_cdf(S_star[ii], AA, lw[ii], "det")))
                est.append((lb_b, np.mean(gm[ii]) / lb_b, np.mean(om[ii]) / lb_b))
            se = np.std(np.asarray(est), axis=0, ddof=1) / np.sqrt(10)
            row.update({"se_lower_bound_det_A": float(se[0]), "se_efficiency_grid_A": float(se[1]),
                        "se_efficiency_oracle_A": float(se[2])})
        if keep_cdf:
            nmax = int(min(H, max(np.quantile(np.where(g < 0, H, g), 0.999), 50)))
            ns = np.arange(1, nmax + 1)
            ecdf = lambda t: np.searchsorted(np.sort(np.where(t < 0, H + 1, t)), ns, side="right") / t.size
            nfr = L.floor_time(lw, AA)
            cdfs.append(pd.DataFrame(dict(model=name, ratio=rho, level=tag, n=ns, floor_det=ecdf(nfr),
                                          np_det=np.pad(b_det, (0, max(0, nmax - b_det.size)), constant_values=1)[:nmax],
                                          np_any=np.pad(b_np, (0, max(0, nmax - b_np.size)), constant_values=1)[:nmax],
                                          floor_rand=rf[:nmax], grid=ecdf(g), oracle=ecdf(o))))
    if len(levels) == 2:
        # each side of the two-sided monitor runs at 2A, but the monitor is valid at n / A over H=: compare
        # the side at 2A with the bound for rules valid at n / A
        row["efficiency_grid_2A_vs_A"] = row["grid_2A_mean"] / row["lower_bound_det_A"]
        row["efficiency_oracle_2A_vs_A"] = row["oracle_2A_mean"] / row["lower_bound_det_A"]
    return row, cdfs


def part_floors():
    rows, cdfs = [], []
    k = 0
    for name in MODELS:
        for rho in RATIOS:
            k += 1
            t0 = time.time()
            row, cd = floors_cell(name, rho, _n(10000), 4000, SEED + k, keep_cdf=True)
            rows.append(row)
            cdfs += cd
            log("floors", name, rho, f"floor {row['floor_A_mean']:.1f}", f"np_det {row['np_det_A']:.1f}",
                f"grid {row['grid_A_mean']:.1f}", f"oracle {row['oracle_A_mean']:.1f}",
                f"eff {row['efficiency_grid_A']:.2f}", f"viol {row['pathwise_violations_grid_A']}",
                f"({time.time() - t0:.0f} s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th5_floors.csv", index=False)
    pd.concat(cdfs).to_csv(OUT / "th5_floor_cdf.csv", index=False)
    print(df[["model", "ratio", "kl_inf", "gamma", "floor_A_mean", "np_any_A", "np_det_A", "np_rand_A",
              "grid_A_mean", "oracle_A_mean", "efficiency_grid_A", "efficiency_oracle_A"]].round(3).to_string(index=False))


# ----------------------------------------------------------------------------
# Part asymptotic
# ----------------------------------------------------------------------------


def _root(rate, AA):
    """Solution of n * rate = log(A / n): the first-order delay with the log n correction that comes from
    alpha_n = n / A."""
    from scipy import optimize

    return float(optimize.brentq(lambda n: n * rate - np.log(AA / n), 1e-9, AA))


def part_asymptotic():
    rows = []
    cells = [("gaussian", 0.5), ("gaussian", 0.8), ("gaussian", 1.5), ("student t3", 0.5), ("student t3", 0.8),
             ("student t3", 1.5)]
    for k, (name, rho) in enumerate(cells):
        model, tau, m0, B, side, law, lam, kl = _setup(name, rho)
        gam = _gamma(law, B, side)
        rng = np.random.default_rng(SEED + 100 + k)
        Ab = (1e3, 1e4, 1e6, 1e8)
        H = int(min(20000, 1.6 * np.log(1e8) / kl + 400))
        Np = _n(4000)
        y = D.sample_y(model, rho, tau, m0, (Np, H), rng)
        lw = L.log_floor_factor(y, B, side)
        u = (y - 1.0) if side == "up" else (1.0 - y)
        S_star = np.cumsum(np.log1p(lam * u), axis=1)
        orac = L.sr_alarm_times(y, [lam], [1.0], side, Ab)
        for i, AA in enumerate(Ab):
            nf = L.floor_time(lw, AA)
            b_det = L.np_cdf(S_star, AA, lw, "det")
            o = orac[i]
            rows.append(dict(model=name, ratio=rho, side=side, A=AA, kl_inf=kl, gamma=gam, paths=Np, horizon=H,
                             floor_mean=float(np.mean(np.where(nf < 0, H, nf))), floor_missing=int(np.sum(nf < 0)),
                             np_det=L.delay_lower_bound(b_det),
                             oracle_mean=float(np.mean(np.where(o < 0, H, o))), oracle_missing=int(np.sum(o < 0)),
                             log_A_over_kl=np.log(AA) / kl, log_A_over_gamma=np.log(AA) / gam,
                             refined_gamma=_root(gam, AA), refined_kl=_root(kl, AA)))
        log("asymptotic", name, rho, H)
    df = pd.DataFrame(rows)
    df["floor_over_logA_gamma"] = df.floor_mean / df.log_A_over_gamma
    df["np_det_over_logA_kl"] = df.np_det / df.log_A_over_kl
    df["oracle_over_logA_kl"] = df.oracle_mean / df.log_A_over_kl
    df["floor_over_logA_kl"] = df.floor_mean / df.log_A_over_kl
    df["floor_over_refined_gamma"] = df.floor_mean / df.refined_gamma
    df["np_det_over_refined_kl"] = df.np_det / df.refined_kl
    df["oracle_over_refined_kl"] = df.oracle_mean / df.refined_kl
    df.to_csv(OUT / "th5_asymptotic.csv", index=False)
    print(df[["model", "ratio", "A", "floor_over_logA_gamma", "floor_over_refined_gamma", "floor_over_logA_kl",
              "np_det_over_logA_kl", "np_det_over_refined_kl", "oracle_over_logA_kl"]].round(3).to_string(index=False))


# ----------------------------------------------------------------------------
# Part boundary
# ----------------------------------------------------------------------------


def part_boundary():
    rows = []
    EPS = (0.5, 0.15, 0.05, 0.01, 1e-3, 1e-6)
    k = 0
    for name in MODELS:
        for rho in (0.5, 0.67, 0.8):
            model, tau, m0, B, side, law, lam, kl = _setup(name, rho)
            t = L.boundary_terms(law)
            lam_max = 1 / (B - 1)
            rng = np.random.default_rng(SEED + 200 + k)
            k += 1
            Np = _n(4000)
            H = 6000
            y = D.sample_y(model, rho, tau, m0, (Np, H), rng)
            eps_list = list(EPS) + [t["eps_star"], max(t["q_B"] / t["G_c"], 1e-300) if t["G_c"] > 0 else np.nan]
            labels = [f"{e:g}" for e in EPS] + ["eps_star", "q_B / G_c"]
            for lab, e in zip(labels, eps_list):
                if not np.isfinite(e) or e <= 0 or e >= 1:
                    continue
                g = float(L.boundary_growth(law, e)[0])
                lo, hi = L.boundary_growth_bounds(law, e)
                al = L.sr_alarm_times(y, [lam_max * (1 - e)], [1.0], "down", (A,))[0]
                wald = (np.log(A) + np.log(B / (B - 1))) / g if g > 0 else np.inf
                blk, blk_n = L.block_delay_bound(np.log(e + (1 - e) * (B - y) / (B - 1)), A)
                rows.append(dict(model=name, ratio=rho, B=B, q_B=t["q_B"], G_c=t["G_c"], eps_star=t["eps_star"],
                                 kl_inf=t["kl_inf"], lam_star_over_lam_max=1 - t["eps_star"], eps_label=lab, eps=e,
                                 growth=g, growth_lower=float(lo[0]), growth_upper=float(hi[0]),
                                 atom_cost=t["q_B"] * np.log(1 / e), growth_over_kl_inf=g / t["kl_inf"],
                                 sim_mean_delay=float(np.mean(np.where(al < 0, H, al))), sim_missing=int(np.sum(al < 0)),
                                 wald_bound=wald, block_bound=blk, block_n=blk_n,
                                 first_order=np.log(A) / g if g > 0 else np.inf, paths=Np))
            log("boundary", name, rho, f"q_B {t['q_B']:.2e} eps* {t['eps_star']:.2e} G_c {t['G_c']:.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th5_boundary.csv", index=False)
    print(df[["model", "ratio", "q_B", "G_c", "eps_star", "kl_inf", "eps_label", "growth", "growth_lower",
              "growth_upper", "sim_mean_delay", "wald_bound", "block_bound"]].to_string(index=False))


# ----------------------------------------------------------------------------
# Part plugin_bound (windows only)
# ----------------------------------------------------------------------------

ZS = (1.0, 2.0, 2.576)
HORIZONS = (98, 390, 1950, 7800)


def _growth_rates(law, tau2, m_hat, nodes=96):
    """Largest growth rate of each side's bets given the window, under the i.i.d. law (quadrature in u)."""
    dist = stats.norm if law == "gaussian" else stats.t(3.0)
    x, w = np.polynomial.legendre.leggauss(nodes)
    x, w = (x + 1) / 2, w / 2
    tau = np.sqrt(tau2)
    Ft = L.law_abs_cdf(law, tau)  # P(|r| <= tau)
    u = Ft[:, None] * x[None, :]
    a = dist.ppf(0.5 + u / 2)
    Y = (a * a) / m_hat[:, None]
    wq = Ft[:, None] * w[None, :]
    Bh = tau2 / m_hat
    pt = 1 - Ft
    gu = np.full(tau2.size, -np.inf)
    for lam in LAM_UP:
        g = np.sum(wq * np.log1p(lam * (Y - 1)), axis=1) + pt * np.log1p(lam * (Bh - 1))
        gu = np.maximum(gu, g)
    gd = np.full(tau2.size, -np.inf)
    for fr in FRAC_DN:
        phi = fr / (Bh - 1)
        g = np.sum(wq * np.log1p(phi[:, None] * (1 - Y)), axis=1) + pt * np.log1p(phi * (1 - Bh))
        gd = np.maximum(gd, g)
    return gu, gd


def part_plugin_bound():
    rows, drows = [], []
    ns = (390, 1440, 2880, 5760, 11520, 23040, 46080, 92160)
    k = 0
    for law in ("gaussian", "t3"):
        pop = L.law_population(law, C)
        CF = L.first_order_constant(pop["sigma_y"], pop["B"])
        for n in ns:
            Wn = _n(20000 if n <= 23040 else 5000)
            rng = np.random.default_rng(SEED + 300 + k)
            k += 1
            t0 = time.time()
            win = L.plugin_windows(law, n, Wn, rng, C, ZS, cs=(n <= 23040 and Wn <= 20000), kolmogorov=True)
            d, Bh = win["d"], win["B_hat"]
            lb_, lmax, pb, pmax = L.bet_constants(pop["B"])
            # check of |d| <= B_hat D_n (deterministic inequality)
            ratio_check = float(np.max(np.abs(d) / (Bh * win["D"])))
            gu, gd = _growth_rates(law, win["tau2"], win["m_hat"])
            drows.append(dict(law=law, n=n, windows=Wn, sd_drift=float(d.std()), sqrt_n_sd_drift=float(d.std() * np.sqrt(n)),
                              sigma_y=pop["sigma_y"], mean_abs_drift=float(np.mean(np.abs(d))),
                              sqrt_n_mean_pos_drift=float(np.mean(np.maximum(d, 0)) * np.sqrt(n)),
                              normal_limit_mean_pos=pop["sigma_y"] / np.sqrt(2 * np.pi),
                              share_growth_up=float(np.mean(gu > 0)), share_growth_down=float(np.mean(gd > 0)),
                              share_growth_either=float(np.mean((gu > 0) | (gd > 0))),
                              max_abs_d_over_Bhat_D=ratio_check, B_hat_min=float(Bh.min()), B_hat_max=float(Bh.max()),
                              cs_z_eff_up=float(np.median((win["cs_up"] / win["cs_mtrue"] - 1) * np.sqrt(n) / pop["sigma_y"]))
                              if "cs_up" in win else np.nan,
                              cs_z_eff_down=float(np.median((1 - win["cs_dn"] / win["cs_mtrue"]) * np.sqrt(n) / pop["sigma_y"]))
                              if "cs_up" in win else np.nan))
            for m in HORIZONS:
                cor1 = L.side_bound(win["m_true"], win["m_hat"], win["m_hat"], win["tau2"], m, A)
                # explicit form of Theorem 8 (a): G_m(x) <= 1 + x (m + 1) e^(x m) / 2
                lam_bar_d = np.where(d > 0, lb_ * d, (0.51 / (Bh - 1)) * (-d))
                xs = np.where(d > 0, lmax * d, (0.95 / (Bh - 1)) * (-d))
                with np.errstate(over="ignore"):
                    expl = np.minimum(1.0, m / A * (1 + (m + 1) / 4 * lam_bar_d * np.exp(np.minimum(xs * m, 700))))
                    expl = expl / (m / A)
                nf = L.noise_free_alarm(d, Bh, m, A)
                hb = L.hold_bold_lower_bound(d, Bh, m, A)
                row = dict(law=law, n=n, m=m, A=A, windows=Wn, m_over_sqrt_n=m / np.sqrt(n), nominal=m / A,
                           corollary1=float(np.mean(cor1)), corollary1_ratio=float(np.mean(cor1)) / (m / A),
                           first_order_ratio=1 + CF * (m + 1) / np.sqrt(n), C_F=CF,
                           explicit_ratio=float(np.mean(expl)),
                           distribution_free_ratio=min(L.distribution_free_bound(n, m, A, C), 1.0) / (m / A),
                           noise_free=float(np.mean(nf)), hold_bold_lower=float(np.mean(hb)),
                           hold_bold_lower_ratio=float(np.mean(hb)) / (m / A))
                for z in ZS:
                    bz = L.side_bound(win["m_true"], win[f"up_{z}"], win[f"dn_{z}"], win["tau2"], m, A)
                    bad = (win["m_true"] > win[f"up_{z}"]) | (win["m_true"] < win[f"dn_{z}"])
                    row[f"inflated_{z}_bound"] = float(np.mean(bz))
                    row[f"inflated_{z}_share_bad"] = float(np.mean(bad))
                    row[f"inflated_{z}_normal_limit"] = 2 * stats.norm.sf(z)
                rows.append(row)
            log("plugin_bound", law, n, f"sd(d) sqrt(n) {d.std() * np.sqrt(n):.3f}", f"growth share {drows[-1]['share_growth_either']:.3f}",
                f"({time.time() - t0:.0f} s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th5_plugin_bound.csv", index=False)
    part_dfree()
    pd.DataFrame(drows).to_csv(OUT / "th5_plugin_drift.csv", index=False)
    print(pd.DataFrame(drows).round(4).to_string(index=False))
    print(df[["law", "n", "m", "m_over_sqrt_n", "corollary1_ratio", "first_order_ratio", "explicit_ratio",
              "distribution_free_ratio", "noise_free", "hold_bold_lower_ratio", "inflated_2.576_bound"]].round(3).to_string(index=False))


def part_dfree():
    """The distribution-free bound of Theorem 8 (b) on a grid of n and m (no randomness)."""
    k1, k2 = L.distribution_free_constants(C)
    lo_b, hi_b = L.b_range(C)
    cdf = L.distribution_free_rate_constant(C)
    dfree = [dict(n=n, m=m, A=A, k1=k1, k2=k2, b_lo=lo_b, b_hi=hi_b, C_DF=cdf, nominal=m / A,
                  bound=min(1.0, L.distribution_free_bound(n, m, A, C)),
                  ratio=L.distribution_free_bound(n, m, A, C) / (m / A),
                  first_order_ratio=1 + cdf * (m + 1) / np.sqrt(n))
             for n in (390, 1440, 5760, 23040, 92160, 10**6, 10**7, 10**8) for m in (5, 10, 39, 98, 390, 1950)]
    pd.DataFrame(dfree).to_csv(OUT / "th5_distribution_free.csv", index=False)
    log("dfree", f"C_DF {cdf:.4f}")


# ----------------------------------------------------------------------------
# Part plugin_sim (monitoring)
# ----------------------------------------------------------------------------


def monitor(rng, law, tau2, m_up, m_dn, T, chunk_t=100):
    """Two-sided default detector (each side at 2A) for P paths and V baselines on the same monitored returns.

    ``tau2``, ``m_up``, ``m_dn``: arrays of shape (P, V). Returns 1-based alarm times (P, V), -1 if none."""
    P, V = tau2.shape
    on = (m_dn > 0) & (tau2 > m_dn)
    Bd = np.where(on, tau2 / np.where(m_dn > 0, m_dn, 1.0), 2.0)
    phi = np.where(on[..., None], FRAC_DN[None, None, :] / (Bd[..., None] - 1), 0.0)
    mdn = np.where(on, m_dn, 1.0)
    Su = np.zeros((P, V, LAM_UP.size))
    Sd = np.zeros((P, V, FRAC_DN.size))
    alarm = np.full((P, V), -1)
    thr = 2 * A
    for t0 in range(0, T, chunk_t):
        r = L.law_draws(law, rng, (P, min(chunk_t, T - t0)))
        r2 = r * r
        for j in range(r.shape[1]):
            X = np.minimum(r2[:, j][:, None], tau2)
            Su = np.minimum((1 + LAM_UP * (X / m_up - 1)[..., None]) * (Su + 1), 1e200)
            Sd = np.minimum((1 + phi * (1 - X / mdn)[..., None]) * (Sd + 1), 1e200)
            hit = (np.maximum(Su.mean(axis=2), Sd.mean(axis=2)) >= thr) & (alarm < 0)
            if hit.any():
                alarm[hit] = t0 + j + 1
    return alarm


def part_plugin_sim():
    rows = []
    ns = (390, 1440, 2880, 5760, 11520, 23040)
    k = 0
    for law in ("gaussian", "t3"):
        pop = L.law_population(law, C)
        for n in ns:
            P = _n(10000)
            rng = np.random.default_rng(SEED + 400 + k)
            k += 1
            t0 = time.time()
            win = L.plugin_windows(law, n, P, rng, C, ZS, cs=True)
            names = ["plug-in"] + [f"inflated z = {z:g}" for z in ZS] + ["confidence sequence", "true baseline"]
            tau2 = np.column_stack([win["tau2"]] * 4 + [win["cs_tau2"], np.full(P, pop["tau"] ** 2)])
            mup = np.column_stack([win["m_hat"]] + [win[f"up_{z}"] for z in ZS] + [win["cs_up"], np.full(P, pop["m0"])])
            mdn = np.column_stack([win["m_hat"]] + [win[f"dn_{z}"] for z in ZS] + [win["cs_dn"], np.full(P, pop["m0"])])
            mtrue = np.column_stack([win["m_true"]] * 4 + [win["cs_mtrue"], np.full(P, pop["m0"])])
            alarm = monitor(rng, law, tau2, mup, mdn, 7800)
            for v, name in enumerate(names):
                for m in HORIZONS:
                    kk = int(np.sum((alarm[:, v] > 0) & (alarm[:, v] <= m)))
                    lo, hi = wilson(kk, P)
                    bnd = L.side_bound(mtrue[:, v], mup[:, v], mdn[:, v], tau2[:, v], m, A)
                    rows.append(dict(law=law, n=n, baseline=name, m=m, A=A, paths=P, m_over_sqrt_n=m / np.sqrt(n),
                                     nominal=m / A, rate=kk / P, ci_low=lo, ci_high=hi, rate_over_nominal=kk / P / (m / A),
                                     corollary1_bound=float(np.mean(bnd)),
                                     share_drift_positive=float(np.mean((mtrue[:, v] > mup[:, v]) | (mtrue[:, v] < mdn[:, v]))),
                                     decrease_side_off=float(np.mean(~((mdn[:, v] > 0) & (tau2[:, v] > mdn[:, v]))))))
            log("plugin_sim", law, n, {nm: round(float(np.mean((alarm[:, v] > 0) & (alarm[:, v] <= 390))), 4) for v, nm in enumerate(names)},
                f"({time.time() - t0:.0f} s)")
            pd.DataFrame(rows).to_csv(OUT / "th5_plugin_sim.csv", index=False)
    df = pd.DataFrame(rows)
    df["rate_above_bound"] = df.ci_low > df.corollary1_bound
    df.to_csv(OUT / "th5_plugin_sim.csv", index=False)
    print(df.pivot_table(index=["law", "n", "m"], columns="baseline", values="rate").round(4).to_string())


# ----------------------------------------------------------------------------
# Part holdbold (simulation check of Proposition 9)
# ----------------------------------------------------------------------------


def holdbold_run(d, Bh, m, rng, L_bold):
    """Run the hold-then-bold stream (hold Y = 1 + d for m - L steps, then the two-point law on {0, B_hat}
    with mean 1 + d) through the two-sided default detector. Returns the alarm indicator within m steps."""
    P = d.size
    phi = FRAC_DN[None, :] / (Bh[:, None] - 1)
    Su = np.zeros((P, LAM_UP.size))
    Sd = np.zeros((P, FRAC_DN.size))
    alarm = np.zeros(P, bool)
    p_top = (1 + d) / Bh
    for t in range(1, m + 1):
        hold = t <= m - L_bold
        Y = np.where(hold, 1 + d, np.where(rng.random(P) < p_top, Bh, 0.0))
        Su = np.minimum((1 + LAM_UP * (Y - 1)[:, None]) * (Su + 1), 1e200)
        Sd = np.minimum((1 + phi * (1 - Y)[:, None]) * (Sd + 1), 1e200)
        alarm |= np.maximum(Su.mean(axis=1), Sd.mean(axis=1)) >= 2 * A
    return alarm


def part_holdbold():
    rows = []
    k = 0
    for n in (390, 2880):
        for m in (390, 1950):
            rng = np.random.default_rng(SEED + 500 + k)
            k += 1
            Wn = _n(4000)
            win = L.plugin_windows("gaussian", n, Wn, rng, C, ())
            d, Bh = win["d"], win["B_hat"]
            # one bold length per window: the one that maximises the proved bound
            Ls = np.unique(np.clip(np.round(np.geomspace(1, m, 60)).astype(int), 1, m))
            best_val = L.noise_free_alarm(d, Bh, m, A).astype(float)
            best_L = np.zeros(Wn, int)
            for Lb in Ls:
                v = L.hold_bold_lower_bound(d, Bh, m, A, holds=[Lb])
                better = v > best_val + 1e-15
                best_val = np.where(better, v, best_val)
                best_L = np.where(better, Lb, best_L)
            hit = np.zeros(Wn, bool)
            for Lb in np.unique(best_L):
                sel = best_L == Lb
                hit[sel] = holdbold_run(d[sel], Bh[sel], m, rng, int(Lb))
            lo, hi = wilson(int(hit.sum()), Wn)
            nf = L.noise_free_alarm(d, Bh, m, A)
            cor1 = L.side_bound(win["m_true"], win["m_hat"], win["m_hat"], win["tau2"], m, A)
            rows.append(dict(law="gaussian", n=n, m=m, windows=Wn, nominal=m / A, proved_lower_bound=float(best_val.mean()),
                             simulated_rate=float(hit.mean()), ci_low=lo, ci_high=hi, noise_free=float(nf.mean()),
                             corollary1_upper=float(cor1.mean()), lower_below_rate=bool(best_val.mean() <= hi)))
            log("holdbold", n, m, rows[-1])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "th5_holdbold.csv", index=False)
    print(df.round(4).to_string(index=False))


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

INK, MUTED, GRIDC = "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a5cc2"]


def _style(ax):
    ax.grid(True, color=GRIDC, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def figures():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(exist_ok=True)
    cdf = pd.read_csv(OUT / "th5_floor_cdf.csv")
    cells = [("gaussian", 0.5), ("gaussian", 0.8), ("student t3", 0.8), ("gaussian", 1.25), ("gaussian", 2.0),
             ("jumps", 0.5)]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4))
    for ax, (mod, rho) in zip(axes.ravel(), cells):
        s = cdf[(cdf.model == mod) & (cdf.ratio == rho) & (cdf.level == "A")]
        ax.plot(s.n, s.floor_det, color=SERIES[1], lw=1.4, label="pathwise floor (every deterministic valid rule)")
        ax.plot(s.n, s.np_det, color=INK, lw=1.6, label="floor and Neyman-Pearson at P* (upper bound on P(T <= n))")
        ax.plot(s.n, s.oracle, color=SERIES[2], lw=1.2, ls="--", label="growth-optimal single bet (simulated)")
        ax.plot(s.n, s.grid, color=SERIES[0], lw=1.2, ls="-.", label="default grid (simulated)")
        ax.set_title(f"{mod}, sd ratio {rho:g}", fontsize=9, color=INK)
        ax.set_xlabel("steps after the change", fontsize=8, color=INK)
        ax.set_ylabel("P(alarm by step n)", fontsize=8, color=INK)
        _style(ax)
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.suptitle("Every deterministic valid rule has its alarm-time CDF below the black curve; e-detectors simulated, A = 7800 [synthetic]",
                 fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    fig.savefig(FIG / "th5_floors.png", dpi=150)
    plt.close(fig)

    pb = pd.read_csv(OUT / "th5_plugin_bound.csv")
    ps = pd.read_csv(OUT / "th5_plugin_sim.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    for ax, law in zip(axes, ("gaussian", "t3")):
        b = pb[pb.law == law].sort_values("m_over_sqrt_n")
        ax.scatter(b.m_over_sqrt_n, b.corollary1_ratio, s=16, color=INK, label="upper bound: Corollary 1 (TH4), all streams in the class")
        xs = np.geomspace(b.m_over_sqrt_n.min(), b.m_over_sqrt_n.max(), 100)
        CF = float(b.C_F.iloc[0])
        ax.plot(xs, 1 + CF * xs, color=MUTED, lw=1.0, ls="--", label="first-order formula 1 + C_F m / sqrt(n)")
        ax.scatter(b.m_over_sqrt_n, np.maximum(b.hold_bold_lower_ratio, 1e-2), s=16, marker="v", color=SERIES[1],
                   label="lower bound: hold-then-bold stream (proved)")
        s = ps[(ps.law == law) & (ps.baseline == "plug-in")]
        ax.scatter(s.m_over_sqrt_n, np.maximum(s.rate_over_nominal, 1e-2), s=22, marker="s", facecolor="none",
                   color=SERIES[0], label="i.i.d. returns, plug-in (simulated)")
        s = ps[(ps.law == law) & (ps.baseline == "true baseline")]
        ax.scatter(s.m_over_sqrt_n, np.maximum(s.rate_over_nominal, 1e-2), s=22, marker="o", facecolor="none",
                   color=SERIES[2], label="i.i.d. returns, true baseline (simulated)")
        ax.axhline(1, color=MUTED, lw=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("m / sqrt(n)", fontsize=9, color=INK)
        ax.set_title(f"{law}: P(false alarm within m) divided by m / A, A = 7800", fontsize=9, color=INK)
        _style(ax)
    axes[0].set_ylabel("ratio to the nominal m / A (values below 0.01 drawn at 0.01)", fontsize=8, color=INK)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(FIG / "th5_plugin_regime.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    parts = sys.argv[1:] or ["collapse", "floors", "asymptotic", "boundary", "plugin_bound", "plugin_sim", "holdbold",
                             "figures"]
    fns = dict(collapse=part_collapse, floors=part_floors, asymptotic=part_asymptotic, boundary=part_boundary,
               plugin_bound=part_plugin_bound, dfree=part_dfree, plugin_sim=part_plugin_sim, holdbold=part_holdbold,
               figures=figures)
    for p in parts:
        t0 = time.time()
        log("start", p)
        fns[p]()
        log("done", p, f"({time.time() - t0:.0f} s)")
