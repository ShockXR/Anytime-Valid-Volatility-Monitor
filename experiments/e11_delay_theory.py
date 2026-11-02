"""E11. Detection delay theory for the betting e-detector: growth rates, price of robustness,
bet design, and the slack in the false alarm bound.

All results are SYNTHETIC (simulations and numerical integration). Seeds are fixed.
Theory: docs/theory_delay.md. Helper functions: src/avmonitor/design.py.

Research questions
  Q11a (part A). For a fixed bet lam, the log wealth grows at g(lam) = E_Q log(1 + lam (Y - 1)) after
        the change. How well do (i) the first-order delay log(A) / g, (ii) the rigorous Wald bound
        and (iii) a second-order Shiryaev-Roberts approximation predict the simulated delays of the
        E4 scenarios (Gaussian, jumps, Student t5; sd ratios 1.25, 1.5, 2 up and 0.8, 0.67, 0.5 down)?
        Where does the first-order approximation fail?
  Q11b (part B). Price of robustness: how does the growth of the robust detector compare with the
        Kullback-Leibler information that drives the Gaussian likelihood ratio detector? Split the
        ratio into truncation loss, nonparametric loss and grid loss, and study the truncation
        constant c in {2, 3, 4, 5, none}. Does the ratio explain the E4 delay penalties?
  Q11c (part C). Bet design: which grid of bets and weights minimises the worst-case ratio of the
        mixture's delay to the oracle (single best bet) delay over a prior range of change sizes
        (sd ratio 1.2 to 2 up, 0.5 to 0.85 down)? How much does it gain over the defaults, by
        simulation (design on formulas, validation on fresh simulations)?
  Q11d (part D). Why is the false alarm bound conservative by about 4 at one day (E9)? Exact
        decomposition of P(false alarm within m) into stopping, uncrossed-mass and overshoot factors,
        and the exact worst case over the whole nonparametric null for a single bet (dynamic
        programming), to separate what is provable from what is only empirical.

Usage: python experiments/e11_delay_theory.py [A|B|C|D ...]   (default: all parts)
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from avmonitor import design as G
from avmonitor.detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
A = 7800.0  # one-sided alarm level used in E4 (P(false alarm in 390 steps) <= 5%)
SEED = 20261111

MODELS = {
    "gaussian": G.ReturnModel("gaussian"),
    "jumps": G.ReturnModel("jumps", p_jump=0.005, jump_sd=8.0),
    "student t5": G.ReturnModel("student_t", df=5.0),
}


def default_grid(direction, B):
    if direction == "up":
        grid = np.asarray(DEFAULT_UP_GRID, float)
    else:
        grid = G.lam_max("down", B) * np.asarray(DEFAULT_DOWN_FRACTIONS, float)
    return grid, np.full(len(grid), 1.0 / len(grid))


def simulate_delays(model, rho, tau, m0, grid, weights, direction, threshold, nu, reps, horizon, rng,
                    batch=1000):
    """Delays (alarm - nu + 1) of the mixture SR with change at index nu (0 = first monitored step).

    Returns (delays of paths that alarmed at or after nu, share of false alarms before nu,
    share of paths with no alarm by the horizon)."""
    delays, early, missing = [], 0, 0
    done = 0
    while done < reps:
        b = min(batch, reps - done)
        pre = G.sample_y(model, 1.0, tau, m0, (b, nu), rng) if nu > 0 else np.zeros((b, 0))
        post = G.sample_y(model, rho, tau, m0, (b, horizon), rng)
        out = G.run_mixture(np.concatenate([pre, post], axis=1), grid, weights, direction, threshold)
        a = out["alarm"]
        early += int(np.sum((a >= 0) & (a < nu)))
        missing += int(np.sum(a < 0))
        ok = a >= nu
        delays.append(a[ok] - nu + 1)
        done += b
    return np.concatenate(delays).astype(float), early / reps, missing / reps


def draw_fn(model, rho, tau, m0):
    return lambda size, rng: G.sample_y(model, rho, tau, m0, size, rng)


# ----------------------------------------------------------------------------
# Part A: growth rates and delay predictions against simulation
# ----------------------------------------------------------------------------


def part_a():
    rng = np.random.default_rng(SEED)
    e4 = pd.read_csv(OUT / "e04b_change.csv")
    cells = [(nm, "up", r) for nm in MODELS for r in (1.25, 1.5, 2.0)]
    cells += [(nm, "down", r) for nm in MODELS for r in (0.8, 0.67, 0.5)]
    cells += [("gaussian", "up", 1.1), ("gaussian", "down", 0.9)]
    grows, rows = [], []
    for nm, direction, rho in cells:
        t0 = time.time()
        model = MODELS[nm]
        tau, m0 = G.baseline(model, 3.0)
        law = G.y_law(model, rho, tau, m0)
        grid, w = default_grid(direction, law.B)
        g = G.growth_rate(grid, law, direction)
        for lam, gk in zip(grid, g):
            grows.append(dict(scenario=nm, direction=direction, ratio=rho, lam=lam, frac_of_max=lam / G.lam_max(direction, law.B),
                              growth=gk))
        k = int(np.argmax(g))
        lam_o, g_o = G.best_bet(law, direction)
        lam_ag = G.approx_best_bet(law, direction)
        g_ag = float(G.growth_rate([lam_ag], law, direction)[0])
        rc_k = G.renewal_constants(grid[k], draw_fn(model, rho, tau, m0), direction, 4000, rng=rng)
        rc_o = G.renewal_constants(lam_o, draw_fn(model, rho, tau, m0), direction, 4000, rng=rng)
        _, s2_k = G.log_factor_moments(grid[k], law, direction)
        _, s2_o = G.log_factor_moments(lam_o, law, direction)
        ub, k_ub = G.wald_delay_bound(A, grid, w, law, direction)
        d1_grid = G.delay_first_order(A, g[k])
        d1_or = G.delay_first_order(A, g_o)
        d2_grid = G.delay_second_order(A, w[k], g[k], rc_k["kappa"], rc_k["zeta"])
        d2_or = G.delay_second_order(A, 1.0, g_o, rc_o["kappa"], rc_o["zeta"])
        d2_or_bm = G.delay_second_order(A, 1.0, g_o, 0.0, G.zeta_brownian(g_o, s2_o))
        horizon = int(min(12000, max(400, 4 * d1_grid)))
        d0, fa0, miss0 = simulate_delays(model, rho, tau, m0, grid, w, direction, A, 0, 3000, horizon, rng)
        do, _, misso = simulate_delays(model, rho, tau, m0, [lam_o], [1.0], direction, A, 0, 3000, horizon, rng)
        d5, fa5, miss5 = simulate_delays(model, rho, tau, m0, grid, w, direction, A, 500, 2000, horizon, rng)
        det = "e-det up (grid)" if direction == "up" else "e-det down (grid)"
        e4row = e4[(e4.scenario == nm) & (np.isclose(e4.ratio, rho)) & (e4.detector == det)]
        rows.append(dict(
            scenario=nm, direction=direction, ratio=rho, tau=tau, m0=m0, B=law.B, mean_Y=law.mean, var_Y=law.var,
            lam_oracle=lam_o, kl_inf=g_o, lam_agrapa=lam_ag, growth_agrapa=g_ag,
            best_grid_lam=grid[k], growth_grid=g[k], kl_trunc=G.kl_truncated(model, rho, tau),
            kl_raw=G.kl_truncated(model, rho, np.inf),
            kappa_grid=rc_k["kappa"], zeta_grid=rc_k["zeta"], zeta_grid_bm=G.zeta_brownian(g[k], s2_k),
            kappa_oracle=rc_o["kappa"], zeta_oracle=rc_o["zeta"], zeta_oracle_bm=G.zeta_brownian(g_o, s2_o),
            pred_first_order_grid=d1_grid, pred_first_order_oracle=d1_or, pred_wald_bound_grid=ub,
            pred_second_order_grid=d2_grid, pred_second_order_oracle=d2_or, pred_second_order_oracle_bm=d2_or_bm,
            sim_nu0_grid_mean=d0.mean(), sim_nu0_grid_se=d0.std(ddof=1) / np.sqrt(d0.size),
            sim_nu0_grid_median=np.median(d0), sim_nu0_grid_missing=miss0,
            sim_nu0_oracle_mean=do.mean(), sim_nu0_oracle_se=do.std(ddof=1) / np.sqrt(do.size), sim_nu0_oracle_missing=misso,
            sim_nu500_grid_mean=d5.mean(), sim_nu500_grid_median=np.median(d5), sim_nu500_false_alarm=fa5,
            sim_nu500_missing=miss5,
            e4_median=float(e4row.median_delay.iloc[0]) if len(e4row) else np.nan,
            e4_mean=float(e4row.mean_delay.iloc[0]) if len(e4row) else np.nan,
            horizon=horizon,
        ))
        print(f"{nm:10s} {direction:4s} {rho:4.2f} g*={g_o:.4f} ggrid={g[k]:.4f} D1={d1_grid:7.1f} D2or={d2_or:6.1f} "
              f"sim_or={do.mean():6.1f} UB={ub:7.1f} sim0={d0.mean():6.1f} sim500={d5.mean():6.1f} "
              f"E4={rows[-1]['e4_mean']:6.1f} ({time.time() - t0:.1f}s)", flush=True)
    pd.DataFrame(grows).to_csv(OUT / "e11a_growth.csv", index=False)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e11a_delay.csv", index=False)
    return df


# ----------------------------------------------------------------------------
# Part B: price of robustness
# ----------------------------------------------------------------------------

C_VALUES = (2.0, 3.0, 4.0, 5.0, np.inf)


def info_row(model, nm, c, direction, rho, rng=None):
    tau, m0 = G.baseline(model, c)
    law = G.y_law(model, rho, tau, m0)
    if direction == "down" and not np.isfinite(law.B):
        return None
    grid, w = default_grid(direction, law.B)
    g = G.growth_rate(grid, law, direction)
    k = int(np.argmax(g))
    lam_o, g_o = G.best_bet(law, direction)
    kl_raw = G.kl_truncated(model, rho, np.inf)
    kl_x = G.kl_truncated(model, rho, tau)
    _, s2_o = G.log_factor_moments(lam_o, law, direction) if g_o > 0 else (0, np.nan)
    _, s2_k = G.log_factor_moments(grid[k], law, direction)
    if direction == "up":
        cap = 1 + grid * (law.B - 1) if np.isfinite(law.B) else np.full(len(grid), np.inf)
        ceiling = np.log(law.B) if np.isfinite(law.B) else np.inf
    else:
        cap = 1 + grid
        ceiling = np.log(law.B / (law.B - 1))
    limit = G.speed_limit(A, cap, w) if np.all(np.isfinite(cap)) else 1
    # second-order delays (Monte Carlo overshoot and memory constants) for the oracle bet and for the
    # dominant grid bet with its weight 1/K (the latter ignores help from the other grid bets)
    d2_rob = d2_grid = np.inf
    if rng is not None and g_o > 0:
        rc_o = G.renewal_constants(lam_o, draw_fn(model, rho, tau, m0), direction, 2000, rng=rng)
        rc_k = G.renewal_constants(grid[k], draw_fn(model, rho, tau, m0), direction, 2000, rng=rng)
        d2_rob = G.delay_second_order(A, 1.0, g_o, rc_o["kappa"], rc_o["zeta"])
        d2_grid = G.delay_second_order(A, w[k], g[k], rc_k["kappa"], rc_k["zeta"])
    return dict(scenario=nm, c=c, direction=direction, ratio=rho, tau=tau, m0=m0, B=law.B,
                clip_share_pre=float(model.sf_abs(tau, 1.0)) if np.isfinite(tau) else 0.0,
                clip_share_post=float(model.sf_abs(tau, rho)) if np.isfinite(tau) else 0.0,
                mean_Y=law.mean, kl_raw=kl_raw, kl_trunc=kl_x, kl_inf=g_o, lam_oracle=lam_o,
                growth_grid=g[k], best_grid_lam=grid[k],
                loss_truncation=kl_raw / kl_x if kl_x > 0 else np.inf,
                loss_nonparametric=kl_x / g_o if g_o > 0 else np.inf,
                loss_grid=g_o / g[k] if g[k] > 0 else np.inf,
                kl_ratio_total=kl_raw / g[k] if g[k] > 0 else np.inf,
                growth_ceiling=ceiling, speed_limit_grid=limit,
                max_log_factor_one_step=float(np.log(cap.max())) if np.all(np.isfinite(cap)) else np.inf,
                pred2_oracle_robust=d2_rob, pred2_grid_robust=d2_grid)


def run_lr_sr(r, ratios, threshold):
    """SR mixture of Gaussian likelihood ratios for sd = ratio (sigma0 = 1), alarm index per path."""
    from avmonitor.detectors import gaussian_lr_edetector

    return gaussian_lr_edetector(r, 1.0, ratios=ratios, threshold=threshold).alarm


def part_b():
    rng = np.random.default_rng(SEED + 1)
    rows = []
    for nm, model in MODELS.items():
        for c in C_VALUES:
            for rho in (1.1, 1.25, 1.5, 2.0, 3.0):
                rows.append(info_row(model, nm, c, "up", rho, rng))
            for rho in (0.5, 0.67, 0.8, 0.9):
                row = info_row(model, nm, c, "down", rho, rng)
                if row is not None:
                    rows.append(row)
    info = pd.DataFrame(rows)
    # parametric reference: second-order delay of the Gaussian LR oracle (Gaussian data)
    lr_rows = []
    for rho in (1.1, 1.25, 1.5, 2.0, 3.0, 0.5, 0.67, 0.8, 0.9):
        kl = float(G.kl_gaussian_scale(rho))
        a = 0.5 * (1 - 1 / rho**2)
        rc = G.renewal_constants_increments(lambda size, r_: -np.log(rho) + a * r_.standard_normal(size) ** 2 * rho**2,
                                            4000, rng=rng)
        s2 = (rho**2 - 1) ** 2 / 2
        lr_rows.append(dict(ratio=rho, kl=kl, kappa_lr=rc["kappa"], zeta_lr=rc["zeta"],
                            zeta_lr_bm=G.zeta_brownian(kl, s2),
                            pred2_lr=G.delay_second_order(A, 1.0, kl, rc["kappa"], rc["zeta"]),
                            pred1_lr=G.delay_first_order(A, kl)))
    lr = pd.DataFrame(lr_rows)
    info = info.merge(lr[["ratio", "kl", "pred2_lr"]].rename(columns={"kl": "kl_gauss"}), on="ratio", how="left")
    info.to_csv(OUT / "e11b_information.csv", index=False)
    lr.to_csv(OUT / "e11b_lr_reference.csv", index=False)
    show = info[(info.scenario == "gaussian")][["c", "direction", "ratio", "kl_raw", "kl_trunc", "kl_inf", "growth_grid",
                                                 "loss_truncation", "loss_nonparametric", "loss_grid", "kl_ratio_total",
                                                 "speed_limit_grid"]]
    print(show.round(3).to_string())

    # simulated price on clean Gaussian data, common random numbers across c and detectors
    sims = []
    model = MODELS["gaussian"]
    for direction, ratios, horizons in (("up", (1.25, 1.5, 2.0), (2500, 800, 400)),
                                        ("down", (0.8, 0.67, 0.5), (2500, 1200, 1200))):
        for rho, horizon in zip(ratios, horizons):
            reps = 3000
            r = rho * rng.standard_normal((reps, horizon))
            a_or = run_lr_sr(r, (rho,), A)
            mix = (1.25, 1.5, 2.0) if direction == "up" else (0.5, 0.7, 0.8)
            a_mix = run_lr_sr(r, mix, A)
            base = dict(direction=direction, ratio=rho, reps=reps, horizon=horizon)
            for det, a in (("Gaussian LR oracle", a_or), ("Gaussian LR mixture (E4)", a_mix)):
                d = a[a >= 0] + 1.0
                sims.append(dict(base, detector=det, c=np.nan, mean_delay=d.mean(), se=d.std(ddof=1) / np.sqrt(d.size),
                                 median_delay=np.median(d), missing=float(np.mean(a < 0))))
            for c in C_VALUES:
                if direction == "down" and not np.isfinite(c):
                    continue
                tau, m0 = G.baseline(model, c)
                y = np.minimum(r**2, tau**2) / m0
                law = G.y_law(model, rho, tau, m0)
                grid, w = default_grid(direction, law.B)
                lam_o, _ = G.best_bet(law, direction)
                for det, gr, ww in (("robust grid", grid, w), ("robust oracle bet", [lam_o], [1.0])):
                    a = G.run_mixture(y, gr, ww, direction, A)["alarm"]
                    d = a[a >= 0] + 1.0
                    sims.append(dict(base, detector=det, c=c, mean_delay=d.mean(), se=d.std(ddof=1) / np.sqrt(d.size),
                                     median_delay=np.median(d), missing=float(np.mean(a < 0))))
            print(direction, rho, "done", flush=True)
    sim = pd.DataFrame(sims)
    sim.to_csv(OUT / "e11b_price_simulated.csv", index=False)
    print(sim.round(2).to_string())

    # jump bursts: sd unchanged, jump probability rises from 0.005 to 0.02 (calibrated on 0.005)
    burst = []
    calm = MODELS["jumps"]
    hot = G.ReturnModel("jumps", p_jump=0.02, jump_sd=8.0)
    for c in C_VALUES:
        tau, m0 = G.baseline(calm, c)
        law = G.y_law(hot, 1.0, tau, m0)
        grid, w = default_grid("up", law.B)
        g = G.growth_rate(grid, law, "up")
        lam_o, g_o = G.best_bet(law, "up")
        law_ref = G.y_law(calm, 1.25, tau, m0)
        g_ref = G.growth_rate(grid, law_ref, "up").max()
        burst.append(dict(c=c, mean_Y_burst=law.mean, growth_grid_burst=g.max(), kl_inf_burst=g_o,
                          growth_grid_sd125=g_ref, burst_over_sd125=g.max() / g_ref))
    burst = pd.DataFrame(burst)
    burst.to_csv(OUT / "e11b_jump_burst.csv", index=False)
    print(burst.round(4).to_string())

    # fine curve in c: efficiency of the robust detector relative to the full-information KL
    curve = []
    for nm, model in MODELS.items():
        for c in list(np.round(np.arange(1.0, 8.01, 0.25), 2)) + [np.inf]:
            for direction, ratios in (("up", (1.25, 1.5, 2.0)), ("down", (0.5, 0.67, 0.8, 0.9))):
                for rho in ratios:
                    row = info_row(model, nm, c, direction, rho)
                    if row is not None:
                        curve.append({k_: row[k_] for k_ in ("scenario", "c", "direction", "ratio", "B", "kl_raw",
                                                             "kl_trunc", "kl_inf", "growth_grid", "speed_limit_grid")})
    curve = pd.DataFrame(curve)
    curve["eff_grid"] = curve.growth_grid / curve.kl_raw
    curve["eff_oracle"] = curve.kl_inf / curve.kl_raw
    curve.to_csv(OUT / "e11b_c_curve.csv", index=False)
    best = curve.loc[curve.groupby(["scenario", "direction", "ratio"]).growth_grid.idxmax()]
    print(best[["scenario", "direction", "ratio", "c", "growth_grid", "eff_grid"]].round(3).to_string())
    return info, sim, burst


PARTS = {"A": part_a, "B": part_b}

if __name__ == "__main__":
    which = sys.argv[1:] or list(PARTS)
    for p in which:
        print(f"=== part {p} ===", flush=True)
        PARTS[p]()
