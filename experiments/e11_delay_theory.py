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
  Q11e (parts E, E2; tracker B73). Does the delay attain the nonparametric lower bound
        log(A) / KL_inf as A grows? delay / (log A / KL_inf) for A = 10^2 to 10^8 under Gaussian,
        Student t3 and jump post-change laws, increases and decreases, for the oracle bet, the
        default grid, a fine (continuous) mixture and the designed grids; second-order reading of
        the gap; KL_inf / KL as the share of the likelihood ratio detector's speed that is kept.

Parts: A (Q11a), B (Q11b), C, C2, C3, C4, C4V (Q11c: Gaussian design, variants, validation,
design pooled over tail shapes, its validation), D, D2 (Q11d), E, E2 (Q11e), FIG (figures).
Usage: python experiments/e11_delay_theory.py [PART ...]   (default: all parts)
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


# ----------------------------------------------------------------------------
# Part C: bet design
# ----------------------------------------------------------------------------

UP_RANGE = (1.2, 2.0)
DOWN_RANGE = (0.5, 0.85)


def mean_delay(model, rho, tau, m0, grid, w, direction, reps, rng, nu=0):
    law = G.y_law(model, rho, tau, m0)
    g = G.growth_rate(grid, law, direction).max()
    horizon = int(min(15000, max(400, 6 * np.log(A) / max(g, 1e-6))))
    d, fa, miss = simulate_delays(model, rho, tau, m0, grid, w, direction, A, nu, reps, horizon, rng)
    return d.mean(), d.std(ddof=1) / np.sqrt(d.size), miss, fa


def part_c():
    model = MODELS["gaussian"]
    tau, m0 = G.baseline(model, 3.0)
    B = tau**2 / m0
    lm = G.lam_max("down", B)
    designs = []
    cands = {"up": {}, "down": {}}
    # 1. design on formulas only
    for direction, rng_, n in (("up", UP_RANGE, 17), ("down", DOWN_RANGE, 15)):
        laws = [G.y_law(model, r, tau, m0) for r in np.linspace(*rng_, n)]
        grid0, w0 = default_grid(direction, B)
        cands[direction]["default"] = (grid0, w0)
        if direction == "up":
            cands[direction]["default without 0.02, 0.05"] = (np.array([0.1, 0.2, 0.35, 0.5]), np.full(4, 0.25))
            options = [(K, fw) for K in (1, 2, 3, 4) for fw in (False, True)]
            cap = 0.95
        else:
            cands[direction]["single 0.95 x max"] = (np.array([0.95 * lm]), np.array([1.0]))
            cands[direction]["0.9 and 0.99 x max"] = (lm * np.array([0.9, 0.99]), np.full(2, 0.5))
            options = [(1, False), (2, False), (3, False)]
            cap = 0.99 * lm
        for K, fw in options:
            lam, w, v = G.optimise_grid(A, laws, K, direction, lam_cap=cap, free_weights=fw)
            cands[direction][f"optimised K={K}{' weighted' if fw else ''}"] = (lam, w)
        for name, (lam, w) in cands[direction].items():
            v1, r1 = G.design_ratio(A, lam, w, laws, direction, use_memory=False)
            v2, r2 = G.design_ratio(A, lam, w, laws, direction, use_memory=True)
            designs.append(dict(direction=direction, design=name, K=len(lam), bets=" ".join(f"{x:.4f}" for x in lam),
                                bets_frac_of_max=" ".join(f"{x / G.lam_max(direction, B):.4f}" for x in lam),
                                weights=" ".join(f"{x:.3f}" for x in w), pred_worst_first_order=v1,
                                pred_worst_second_order=v2))
            print(designs[-1], flush=True)
    pd.DataFrame(designs).to_csv(OUT / "e11c_design.csv", index=False)

    # 2. selection by simulation (seed S1), oracle single bet as the reference
    sel = []
    rng = np.random.default_rng(SEED + 2)
    sel_rhos = {"up": (1.2, 1.35, 1.5, 1.75, 2.0), "down": (0.5, 0.6, 0.7, 0.85)}
    for direction in ("up", "down"):
        for rho in sel_rhos[direction]:
            law = G.y_law(model, rho, tau, m0)
            lam_o, _ = G.best_bet(law, direction)
            seed = int(rng.integers(1 << 31))
            d_or = mean_delay(model, rho, tau, m0, [lam_o], [1.0], direction, 2000, np.random.default_rng(seed))[0]
            for name, (lam, w) in cands[direction].items():
                d = mean_delay(model, rho, tau, m0, lam, w, direction, 2000, np.random.default_rng(seed))
                sel.append(dict(stage="selection", direction=direction, ratio=rho, design=name, mean_delay=d[0], se=d[1],
                                missing=d[2], oracle_delay=d_or, ratio_to_oracle=d[0] / d_or))
    sel = pd.DataFrame(sel)
    worst = sel.groupby(["direction", "design"]).ratio_to_oracle.max().reset_index()
    print(worst.sort_values(["direction", "ratio_to_oracle"]).round(3).to_string(), flush=True)
    chosen = {}
    for direction in ("up", "down"):
        wd = worst[worst.direction == direction].copy()
        wd["K"] = [len(cands[direction][n][0]) for n in wd.design]
        best = wd.ratio_to_oracle.min()
        # fewest bets among designs within 1% of the best simulated worst case
        pick = wd[wd.ratio_to_oracle <= best + 0.01].sort_values(["K", "ratio_to_oracle"]).iloc[0]
        chosen[direction] = pick.design
    print("chosen:", chosen, flush=True)

    # 3. validation on fresh simulations (seed S2)
    val = []
    rng = np.random.default_rng(SEED + 3)
    val_rhos = {"up": (1.1, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6, 1.75, 2.0, 3.0),
                "down": (0.4, 0.5, 0.6, 0.67, 0.75, 0.8, 0.85, 0.9)}
    for nm, mdl in MODELS.items():
        t_, m_ = G.baseline(mdl, 3.0)
        B_ = t_**2 / m_
        lm_ = G.lam_max("down", B_)
        for direction in ("up", "down"):
            rhos = val_rhos[direction] if nm == "gaussian" else ((1.25, 1.5, 2.0) if direction == "up" else (0.5, 0.67, 0.8))
            for rho in rhos:
                law = G.y_law(mdl, rho, t_, m_)
                lam_o, _ = G.best_bet(law, direction)
                for nu in ((0, 500) if (nm == "gaussian" and rho in (1.25, 1.5, 2.0, 0.5, 0.67, 0.8)) else (0,)):
                    seed = int(rng.integers(1 << 31))
                    d_or = mean_delay(mdl, rho, t_, m_, [lam_o], [1.0], direction, 4000, np.random.default_rng(seed), nu)[0]
                    for name in ("default", chosen[direction]):
                        lam, w = cands[direction][name]
                        if direction == "down":  # same fractions of the model's own largest bet
                            lam = lam / lm * lm_
                        d = mean_delay(mdl, rho, t_, m_, lam, w, direction, 4000, np.random.default_rng(seed), nu)
                        val.append(dict(stage="validation", scenario=nm, direction=direction, ratio=rho, nu=nu,
                                        in_design_range=bool((UP_RANGE[0] <= rho <= UP_RANGE[1]) if direction == "up"
                                                             else (DOWN_RANGE[0] <= rho <= DOWN_RANGE[1])),
                                        design=name, mean_delay=d[0], se=d[1], missing=d[2], false_alarm_before=d[3],
                                        oracle_delay=d_or, ratio_to_oracle=d[0] / d_or))
            print(nm, direction, "validated", flush=True)
    val = pd.DataFrame(val)
    # null check: one-day false alarm rate (validity holds for any fixed grid; reported for completeness)
    nul = []
    rng = np.random.default_rng(SEED + 4)
    y = G.sample_y(model, 1.0, tau, m0, (20000, 390), rng)
    for direction in ("up", "down"):
        for name in ("default", chosen[direction]):
            lam, w = cands[direction][name]
            a = G.run_mixture(y, lam, w, direction, A)["alarm"]
            nul.append(dict(direction=direction, design=name, one_day_false_alarm=float(np.mean(a >= 0)), bound=390 / A))
    nul = pd.DataFrame(nul)
    print(nul.to_string())
    pd.concat([sel, val], ignore_index=True).to_csv(OUT / "e11c_design_sim.csv", index=False)
    nul.to_csv(OUT / "e11c_null_check.csv", index=False)
    summ = val[val.nu == 0].pivot_table(index=["scenario", "direction", "ratio"], columns="design",
                                         values="ratio_to_oracle").round(3)
    print(summ.to_string())
    return designs, sel, val


def part_c2():
    """Robustness check of a few simple variants of the chosen grids (fresh seed), reported as such.

    The single decrease bet at 0.99 of the largest admissible bet won the Gaussian design, but a bet
    so close to the edge loses log(1 - 0.99) at every capped return, so it is fragile when the
    post-change law has mass at the cap (fat tails). Growth at the cap is computed exactly, then
    simulated with a new seed for all three models."""
    rng = np.random.default_rng(SEED + 7)
    up_c = {"default": (np.asarray(DEFAULT_UP_GRID), np.full(6, 1 / 6)),
            "optimised K=2 weighted (0.2565, 0.4893; 0.569, 0.431)": (np.array([0.2565, 0.4893]), np.array([0.569, 0.431])),
            "rounded (0.25, 0.5) equal weights": (np.array([0.25, 0.5]), np.array([0.5, 0.5])),
            "wider (0.1, 0.25, 0.5) equal weights": (np.array([0.1, 0.25, 0.5]), np.full(3, 1 / 3))}
    dn_f = {"default": (np.asarray(DEFAULT_DOWN_FRACTIONS), np.full(5, 0.2)),
            "single 0.99": (np.array([0.99]), np.array([1.0])),
            "single 0.95": (np.array([0.95]), np.array([1.0])),
            "0.9 and 0.99 equal weights": (np.array([0.9, 0.99]), np.array([0.5, 0.5]))}
    rows, grow = [], []
    for nm, mdl in MODELS.items():
        tau, m0 = G.baseline(mdl, 3.0)
        lm = G.lam_max("down", tau**2 / m0)
        for direction, rhos, cands in (("up", (1.1, 1.2, 1.35, 1.5, 1.75, 2.0, 3.0), up_c),
                                       ("down", (0.5, 0.67, 0.8, 0.85, 0.9), dn_f)):
            for rho in rhos:
                law = G.y_law(mdl, rho, tau, m0)
                lam_o, g_o = G.best_bet(law, direction)
                seed = int(rng.integers(1 << 31))
                d_or = mean_delay(mdl, rho, tau, m0, [lam_o], [1.0], direction, 4000, np.random.default_rng(seed))[0]
                for name, (lam, w) in cands.items():
                    lam_ = lam * lm if direction == "down" else lam
                    g = G.growth_rate(lam_, law, direction)
                    d = mean_delay(mdl, rho, tau, m0, lam_, w, direction, 4000, np.random.default_rng(seed))
                    rows.append(dict(scenario=nm, direction=direction, ratio=rho, design=name, mean_delay=d[0], se=d[1],
                                     missing=d[2], oracle_delay=d_or, ratio_to_oracle=d[0] / d_or,
                                     growth_best=float(g.max()), kl_inf=g_o, cap_mass=float(mdl.sf_abs(tau, rho))))
            print(nm, direction, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e11c2_variants.csv", index=False)
    print(df.pivot_table(index=["scenario", "direction", "ratio"], columns="design", values="ratio_to_oracle").round(3).to_string())
    return df


FINAL_UP = (np.array([0.25, 0.5]), np.array([0.5, 0.5]))
FINAL_DOWN_FRACTIONS = (np.array([0.9, 0.99]), np.array([0.5, 0.5]))


def part_c3():
    """Fresh-seed validation of the final recommendation (chosen after parts C and C2):
    increase bets (0.25, 0.5) with equal weights; decrease bets 0.9 and 0.99 of the largest
    admissible bet with equal weights. Delays with change at step 1 and at step 501 (E4 setting),
    3000 paths, all three models, sd ratios inside the design ranges."""
    rng = np.random.default_rng(SEED + 9)
    rows = []
    for nm, mdl in MODELS.items():
        tau, m0 = G.baseline(mdl, 3.0)
        lm = G.lam_max("down", tau**2 / m0)
        for direction, rhos in (("up", (1.2, 1.25, 1.5, 1.75, 2.0)), ("down", (0.5, 0.6, 0.67, 0.8, 0.85))):
            grid0, w0 = default_grid(direction, tau**2 / m0)
            fin = FINAL_UP if direction == "up" else (FINAL_DOWN_FRACTIONS[0] * lm, FINAL_DOWN_FRACTIONS[1])
            for rho in rhos:
                law = G.y_law(mdl, rho, tau, m0)
                lam_o, _ = G.best_bet(law, direction)
                for nu in (0, 500):
                    seed = int(rng.integers(1 << 31))
                    d_or = mean_delay(mdl, rho, tau, m0, [lam_o], [1.0], direction, 3000, np.random.default_rng(seed), nu)
                    d_def = mean_delay(mdl, rho, tau, m0, grid0, w0, direction, 3000, np.random.default_rng(seed), nu)
                    d_fin = mean_delay(mdl, rho, tau, m0, fin[0], fin[1], direction, 3000, np.random.default_rng(seed), nu)
                    rows.append(dict(scenario=nm, direction=direction, ratio=rho, nu=nu, oracle=d_or[0], default=d_def[0],
                                     default_se=d_def[1], final=d_fin[0], final_se=d_fin[1],
                                     default_over_oracle=d_def[0] / d_or[0], final_over_oracle=d_fin[0] / d_or[0],
                                     final_over_default=d_fin[0] / d_def[0], false_alarm_default=d_def[3],
                                     false_alarm_final=d_fin[3]))
            print(nm, direction, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e11c3_final_validation.csv", index=False)
    print(df.round(3).to_string())
    return df


# ----------------------------------------------------------------------------
# Part C4: bet design pooled over tail shapes (added after C3 and E showed that the Gaussian-tuned
# grids lose speed under Student t3 and t5 data)
# ----------------------------------------------------------------------------

POOL_MODELS = {**MODELS, "student t3": G.ReturnModel("student_t", df=3.0)}
GAUSS_TUNED_UP = FINAL_UP
GAUSS_TUNED_DOWN_FRACTIONS = FINAL_DOWN_FRACTIONS


def _pooled(direction, rhos):
    out = []
    for nm, mdl in POOL_MODELS.items():
        tau, m0 = G.baseline(mdl, 3.0)
        laws = [G.y_law(mdl, r, tau, m0) for r in rhos]
        scale = G.lam_max("down", tau**2 / m0) if direction == "down" else 1.0
        out.append((nm, scale, laws, G.oracle_delays(A, laws, direction)))
    return out


def pooled_ratio(spec, w, pooled, direction):
    """Worst predicted delay / oracle delay over all models and change sizes. ``spec`` is the
    grid of bets (increase) or of fractions of each model's largest admissible bet (decrease)."""
    per = [G.design_ratio(A, np.asarray(spec) * s, w, laws, direction, oracle=orc)[0] for _, s, laws, orc in pooled]
    return max(per), per


def optimise_pooled(pooled, K, direction, cap):
    from scipy import optimize

    lattice = np.geomspace(0.02 * cap, cap, 30)
    best_v, best = np.inf, None
    w = np.full(K, 1.0 / K)
    for i, lo in enumerate(lattice):
        for hi in lattice[i:]:
            if (K == 1) != (hi == lo):
                continue
            spec = np.geomspace(lo, hi, K) if K > 1 else np.array([lo])
            v = pooled_ratio(spec, w, pooled, direction)[0]
            if v < best_v:
                best_v, best = v, spec

    def obj(z):
        return pooled_ratio(np.sort(cap / (1.0 + np.exp(-z))), w, pooled, direction)[0]

    start = np.minimum(best, cap * (1 - 1e-6))
    res = optimize.minimize(obj, np.log(start / (cap - start)), method="Nelder-Mead",
                            options={"maxiter": 1500, "xatol": 1e-4, "fatol": 1e-6})
    if res.fun < best_v:
        best_v, best = float(res.fun), np.sort(cap / (1.0 + np.exp(-res.x)))
    return best, w, best_v


def part_c4():
    """Design step (formulas only, no simulation). Pre-registered rule: among the candidates, pick
    the smallest predicted worst ratio over the pooled models (Gaussian, jumps, t5, t3) and change
    sizes (sd ratio 1.2 to 2 up, 0.5 to 0.85 down); designs within 1% of the best count as ties and
    the one with fewer bets wins; optimised bets are then rounded to two digits (increase) or to
    the nearest 0.05 (decrease fractions) and the rounded grid is the one validated in part C4V."""
    rows, chosen = [], {}
    for direction, rhos, cap in (("up", np.linspace(*UP_RANGE, 9), 0.95), ("down", np.linspace(*DOWN_RANGE, 8), 0.99)):
        pooled = _pooled(direction, rhos)
        if direction == "up":
            cands = {"default": (np.asarray(DEFAULT_UP_GRID), np.full(6, 1 / 6)),
                     "Gaussian-tuned (0.25, 0.5)": GAUSS_TUNED_UP,
                     "(0.1, 0.25, 0.5)": (np.array([0.1, 0.25, 0.5]), np.full(3, 1 / 3))}
        else:
            cands = {"default": (np.asarray(DEFAULT_DOWN_FRACTIONS), np.full(5, 0.2)),
                     "Gaussian-tuned (0.9, 0.99)": GAUSS_TUNED_DOWN_FRACTIONS,
                     "(0.5, 0.75, 0.95)": (np.array([0.5, 0.75, 0.95]), np.full(3, 1 / 3))}
        for K in (1, 2, 3, 4):
            spec, w, v = optimise_pooled(pooled, K, direction, cap)
            cands[f"pooled optimum K={K}"] = (spec, w)
            step = 100.0 if direction == "up" else 20.0
            rounded = np.unique(np.clip(np.round(spec * step) / step, 1 / step, cap))
            cands[f"pooled optimum K={K}, rounded"] = (rounded, np.full(len(rounded), 1 / len(rounded)))
            print(direction, K, np.round(spec, 4), round(v, 4), flush=True)
        for name, (spec, w) in cands.items():
            v, per = pooled_ratio(spec, w, pooled, direction)
            rows.append(dict(direction=direction, design=name, K=len(spec), bets=" ".join(f"{x:.4f}" for x in spec),
                             weights=" ".join(f"{x:.3f}" for x in w), pred_worst_pooled=v,
                             **{f"pred_worst_{nm}": p for (nm, _, _, _), p in zip(pooled, per)}))
        df = pd.DataFrame([r for r in rows if r["direction"] == direction])
        ok = df[df.design.str.endswith("rounded") | ~df.design.str.startswith("pooled")]
        best = ok.pred_worst_pooled.min()
        pick = ok[ok.pred_worst_pooled <= best * 1.01].sort_values(["K", "pred_worst_pooled"]).iloc[0]
        chosen[direction] = (pick.design, pick.bets)
    df = pd.DataFrame(rows)
    df["chosen"] = [chosen[d][0] == n for d, n in zip(df.direction, df.design)]
    df.to_csv(OUT / "e11c4_pooled_design.csv", index=False)
    print(df.round(3).to_string())
    print("chosen:", chosen)
    return df


# chosen by the design step of part C4 (see results/e11c4_pooled_design.csv); validated in part C4V
ROBUST_UP = (np.array([0.16, 0.39]), np.array([0.5, 0.5]))
ROBUST_DOWN_FRACTIONS = (np.array([0.5, 0.85]), np.array([0.5, 0.5]))


def part_c4v():
    """Fresh-seed validation of the pooled (robust) grids against the default and the Gaussian-tuned
    grids, all four models, change at step 1 and at step 501, 2000 paths per cell."""
    rng = np.random.default_rng(SEED + 12)
    rows = []
    for nm, mdl in POOL_MODELS.items():
        tau, m0 = G.baseline(mdl, 3.0)
        B = tau**2 / m0
        lm = G.lam_max("down", B)
        for direction, rhos in (("up", (1.1, 1.2, 1.5, 2.0, 3.0)), ("down", (0.4, 0.5, 0.67, 0.85, 0.9))):
            dets = {"default": default_grid(direction, B)}
            if direction == "up":
                dets["Gaussian-tuned"] = GAUSS_TUNED_UP
                dets["robust"] = ROBUST_UP
            else:
                dets["Gaussian-tuned"] = (GAUSS_TUNED_DOWN_FRACTIONS[0] * lm, GAUSS_TUNED_DOWN_FRACTIONS[1])
                dets["robust"] = (ROBUST_DOWN_FRACTIONS[0] * lm, ROBUST_DOWN_FRACTIONS[1])
            for rho in rhos:
                law = G.y_law(mdl, rho, tau, m0)
                lam_o, g_o = G.best_bet(law, direction)
                inside = bool(UP_RANGE[0] <= rho <= UP_RANGE[1]) if direction == "up" else bool(DOWN_RANGE[0] <= rho <= DOWN_RANGE[1])
                for nu in (0, 500):
                    seed = int(rng.integers(1 << 31))
                    d_or = mean_delay(mdl, rho, tau, m0, [lam_o], [1.0], direction, 2000, np.random.default_rng(seed), nu)
                    for det, (gr, w) in dets.items():
                        d = mean_delay(mdl, rho, tau, m0, gr, w, direction, 2000, np.random.default_rng(seed), nu)
                        rows.append(dict(scenario=nm, direction=direction, ratio=rho, in_design_range=inside, nu=nu,
                                         detector=det, mean_delay=d[0], se=d[1], missing=d[2], false_alarm_before=d[3],
                                         oracle_delay=d_or[0], ratio_to_oracle=d[0] / d_or[0],
                                         growth_best=float(G.growth_rate(np.asarray(gr), law, direction).max()),
                                         kl_inf=g_o))
            print(nm, direction, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e11c4_robust_validation.csv", index=False)
    print(df.pivot_table(index=["scenario", "direction", "ratio", "nu"], columns="detector",
                         values="ratio_to_oracle").round(3).to_string())
    return df


# ----------------------------------------------------------------------------
# Part D: why the false alarm bound is conservative
# ----------------------------------------------------------------------------

HORIZONS = (98, 390, 1950, 3900, 7800)


def decomposition(y, sides, A_total, horizons=HORIZONS):
    """Exact optional-stopping decomposition of P(alarm within m) for an e-detector under an exact null.

    ``sides`` is a list of (direction, grid, weights, side_threshold) and the combined statistic is
    ``N_t = mean over sides of M_side``; the alarm is the first time any side reaches its threshold.
    Under the exact null N_t - t is a martingale, so for every horizon m
        P(sigma <= m) E[N_sigma | sigma <= m] + E[N_m; sigma > m] = E[min(sigma, m)].
    Returns one row per horizon with every term, the three factors and a per-bet split of the
    uncrossed mass."""
    reps, n = y.shape
    S = [np.zeros((reps, len(gr))) for _, gr, _, _ in sides]
    alarm = np.full(reps, -1)
    N_alarm = np.full(reps, np.nan)
    own_alarm = np.full(reps, np.nan)
    active = np.ones(reps, bool)
    snap = {}
    hset = set(horizons)
    for t in range(n):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        Ms = []
        hit = np.zeros(idx.size, bool)
        for j, (direction, gr, w, thr) in enumerate(sides):
            sign = 1.0 if direction == "up" else -1.0
            f = 1.0 + np.asarray(gr)[None, :] * (sign * (y[idx, t] - 1.0))[:, None]
            S[j][idx] = f * (S[j][idx] + 1.0)
            M = S[j][idx] @ np.asarray(w)
            Ms.append(M)
            hit |= M >= thr
        Ms = np.array(Ms)
        alarm[idx[hit]] = t
        N_alarm[idx[hit]] = Ms[:, hit].mean(axis=0)
        ratio_to_thr = Ms[:, hit] / np.array([thr for *_, thr in sides])[:, None]
        own_alarm[idx[hit]] = Ms[:, hit][np.argmax(ratio_to_thr, axis=0), np.arange(hit.sum())] / len(sides)
        active[idx[hit]] = False
        if (t + 1) in hset:
            snap[t + 1] = [Sj.copy() for Sj in S]
    rows = []
    for m in horizons:
        crossed = (alarm >= 0) & (alarm < m)
        p = crossed.mean()
        stop = np.where(crossed, alarm + 1, m).mean()  # E[min(sigma, m)]
        Sm = snap[m]
        N_m = np.mean([Sj @ np.asarray(w) for Sj, (_, _, w, _) in zip(Sm, sides)], axis=0)
        uncrossed = np.mean(np.where(crossed, 0.0, N_m))
        e_over = float(np.nanmean(N_alarm[crossed])) if p > 0 else np.nan
        row = dict(horizon=m, bound=m / A_total, false_alarm=p, E_stop=stop, E_uncrossed=uncrossed,
                   E_N_at_alarm=e_over, E_own_side_at_alarm=float(np.nanmean(own_alarm[crossed])) if p > 0 else np.nan,
                   identity_lhs=p * (e_over if p > 0 else 0.0) + uncrossed,
                   F1_stopping=stop / m, F2_uncrossed=1 - uncrossed / stop,
                   F3_overshoot=A_total / e_over if p > 0 else np.nan)
        row["product"] = row["F1_stopping"] * row["F2_uncrossed"] * row["F3_overshoot"]
        row["rate_over_bound"] = p / (m / A_total)
        # share of the uncrossed mass held by each bet (averaged over sides)
        for (direction, gr, w, _), Sj in zip(sides, Sm):
            for lam, wk, col in zip(gr, w, Sj.T):
                row[f"uncrossed_{direction}_{lam:.4f}"] = float(np.mean(np.where(crossed, 0.0, wk * col))) / len(sides) / max(uncrossed, 1e-12)
        rows.append(row)
    return rows


def part_d():
    model = MODELS["gaussian"]
    tau, m0 = G.baseline(model, 3.0)
    B = tau**2 / m0
    up, wu = default_grid("up", B)
    dn, wd = default_grid("down", B)
    rng = np.random.default_rng(SEED + 5)
    reps, n = 10000, max(HORIZONS)
    configs = {
        "two-sided default grids (E9)": [("up", up, wu, 2 * A), ("down", dn, wd, 2 * A)],
        "increase side only, default grid": [("up", up, wu, A)],
        "decrease side only, default grid": [("down", dn, wd, A)],
        "increase side only, single bet 0.5": [("up", [0.5], [1.0], A)],
        "increase side only, single bet 0.2": [("up", [0.2], [1.0], A)],
        "increase side only, single bet 0.05": [("up", [0.05], [1.0], A)],
    }
    rows = []
    for name, sides in configs.items():
        t0 = time.time()
        acc = []
        for b in range(reps // 2500):
            y = G.sample_y(model, 1.0, tau, m0, (2500, n), rng)
            acc.append(pd.DataFrame(decomposition(y, sides, A)))
        # average the batch-level expectations (equal batch sizes), recompute ratios
        df = pd.concat(acc).groupby("horizon").mean().reset_index()
        df["F1_stopping"] = df.E_stop / df.horizon
        df["F2_uncrossed"] = 1 - df.E_uncrossed / df.E_stop
        df["F3_overshoot"] = A / df.E_N_at_alarm
        df["product"] = df.F1_stopping * df.F2_uncrossed * df.F3_overshoot
        df["rate_over_bound"] = df.false_alarm / df.bound
        df.insert(0, "detector", name)
        df["paths"] = reps
        rows.append(df)
        print(name, f"({time.time() - t0:.0f}s)", flush=True)
        print(df[["horizon", "false_alarm", "bound", "rate_over_bound", "F1_stopping", "F2_uncrossed",
                  "F3_overshoot", "product", "E_stop", "identity_lhs"]].round(4).to_string(), flush=True)
    dec = pd.concat(rows, ignore_index=True)
    dec.to_csv(OUT / "e11d_decomposition.csv", index=False)

    # exact worst case over the nonparametric null for single bets (lower bound by DP)
    dp = []
    for direction, lam in (("up", 0.5), ("up", 0.2), ("up", 0.05), ("down", 0.95 * G.lam_max("down", B))):
        t0 = time.time()
        hmax = 3900  # at m = A the bound is attained exactly by Y = 1 (M_t = t), see docs/theory_delay.md
        lo, upb = G.worst_case_false_alarm(lam, B, A, hmax, direction, n_states=700, n_left=30, n_right=90,
                                           compute_upper=False)
        for m in [h for h in HORIZONS if h <= hmax]:
            dp.append(dict(direction=direction, lam=lam, B=B, A=A, horizon=m, bound=m / A,
                           worst_case_lower=lo[m - 1], worst_case_upper=m / A,
                           lower_over_bound=lo[m - 1] / (m / A)))
        print(direction, lam, f"DP done ({time.time() - t0:.0f}s)", [round(lo[m - 1] / (m / A), 3) for m in HORIZONS if m <= hmax], flush=True)
    dp = pd.DataFrame(dp)
    dp.to_csv(OUT / "e11d_worst_case_dp.csv", index=False)

    # i.i.d. two-point adversaries (lower bounds on the worst case for the actual grids)
    adv = []
    rng = np.random.default_rng(SEED + 6)
    m = 390
    for a_ in (0.0, 0.25, 0.5, 0.75, 0.9):
        for b_ in (1.25, 1.5, 2.0, 3.0, 5.0, B):
            p_b = (1 - a_) / (b_ - a_)
            u = rng.random((20000, m))
            y = np.where(u < p_b, b_, a_)
            for name, (direction, gr, w) in {"up default grid": ("up", up, wu), "up single 0.5": ("up", [0.5], [1.0]),
                                             "down default grid": ("down", dn, wd)}.items():
                al = G.run_mixture(y, gr, w, direction, A)["alarm"]
                adv.append(dict(detector=name, low=a_, high=b_, p_high=p_b, horizon=m, false_alarm=float(np.mean(al >= 0)),
                                bound=m / A))
    adv = pd.DataFrame(adv)
    adv.to_csv(OUT / "e11d_two_point_adversary.csv", index=False)
    print(adv.loc[adv.groupby("detector").false_alarm.idxmax()].to_string())
    return dec, dp, adv


def part_d2():
    """Two-sided split and grid choice at one day (Gaussian null, true baseline, 40,000 paths)."""
    model = MODELS["gaussian"]
    tau, m0 = G.baseline(model, 3.0)
    B = tau**2 / m0
    lm = G.lam_max("down", B)
    up, wu = default_grid("up", B)
    dn, wd = default_grid("down", B)
    designs = {
        "default grids": ((up, wu), (dn, wd)),
        "E11 recommended grids (up 0.25/0.5, down 0.9/0.99 x max)": (FINAL_UP, (FINAL_DOWN_FRACTIONS[0] * lm,
                                                                                 FINAL_DOWN_FRACTIONS[1])),
        "single bets (up 0.5, down 0.95 x max)": ((np.array([0.5]), np.array([1.0])),
                                                  (np.array([0.95 * lm]), np.array([1.0]))),
    }
    rng = np.random.default_rng(SEED + 8)
    rows = []
    m = 390
    y = G.sample_y(model, 1.0, tau, m0, (40000, m), rng)
    for name, ((gu, wu_), (gd, wd_)) in designs.items():
        res = {}
        for thr_name, thr in (("A", A), ("2A", 2 * A)):
            res[("up", thr_name)] = G.run_mixture(y, gu, wu_, "up", thr)["alarm"]
            res[("down", thr_name)] = G.run_mixture(y, gd, wd_, "down", thr)["alarm"]
        a_up, a_dn = res[("up", "2A")], res[("down", "2A")]
        two = np.where((a_up >= 0) | (a_dn >= 0), 1, 0)
        rows.append(dict(design=name, horizon=m, bound_total=m / A,
                         up_at_A=float(np.mean(res[("up", "A")] >= 0)), down_at_A=float(np.mean(res[("down", "A")] >= 0)),
                         up_at_2A=float(np.mean(a_up >= 0)), down_at_2A=float(np.mean(a_dn >= 0)),
                         both_at_2A=float(np.mean((a_up >= 0) & (a_dn >= 0))), two_sided=float(two.mean()),
                         paths=y.shape[0]))
    df = pd.DataFrame(rows)
    df["two_sided_over_bound"] = df.two_sided / df.bound_total
    df.to_csv(OUT / "e11d_two_sided_split.csv", index=False)
    print(df.round(4).to_string())
    return df


# ----------------------------------------------------------------------------
# Part E: first-order optimality as A grows (tracker B73)
# ----------------------------------------------------------------------------

A_LIST = 10.0 ** np.arange(2, 9)


def first_passages(y, grid, weights, direction, thresholds):
    """First index at which the SR mixture reaches each threshold (one run, all thresholds), -1 if never."""
    reps, n = y.shape
    grid = np.asarray(grid, float)
    weights = np.asarray(weights, float)
    sign = 1.0 if direction == "up" else -1.0
    thresholds = np.asarray(thresholds, float)
    S = np.zeros((reps, len(grid)))
    out = np.full((reps, len(thresholds)), -1)
    active = np.ones(reps, bool)
    for t in range(n):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        f = 1.0 + grid[None, :] * (sign * (y[idx, t] - 1.0))[:, None]
        S[idx] = f * (S[idx] + 1.0)
        M = S[idx] @ weights
        newly = (M[:, None] >= thresholds[None, :]) & (out[idx] < 0)
        rows, cols = np.nonzero(newly)
        out[idx[rows], cols] = t
        active[idx[out[idx, -1] >= 0]] = False
    return out


def part_e():
    rng = np.random.default_rng(SEED + 10)
    models = {"gaussian": MODELS["gaussian"], "student t3": G.ReturnModel("student_t", df=3.0), "jumps": MODELS["jumps"]}
    cells = [("up", 1.25), ("up", 1.5), ("up", 2.0), ("down", 0.5), ("down", 0.8)]
    rows = []
    for nm, mdl in models.items():
        tau, m0 = G.baseline(mdl, 3.0)
        B = tau**2 / m0
        lm = G.lam_max("down", B)
        for direction, rho in cells:
            law = G.y_law(mdl, rho, tau, m0)
            lam_o, kl_inf = G.best_bet(law, direction)
            kl = G.kl_truncated(mdl, rho, np.inf)
            grid0, w0 = default_grid(direction, B)
            cap = 0.95 if direction == "up" else 0.99 * lm
            cont = (np.arange(100) + 0.5) / 100 * cap  # fine uniform grid: a discretised continuous mixture
            dets = {"oracle bet": ([lam_o], [1.0]), "default grid": (grid0, w0),
                    "continuous mixture (100 bets, uniform on (0, cap))": (cont, np.full(100, 0.01)),
                    "Gaussian-tuned grid": FINAL_UP if direction == "up" else (FINAL_DOWN_FRACTIONS[0] * lm,
                                                                              FINAL_DOWN_FRACTIONS[1]),
                    "robust grid": ROBUST_UP if direction == "up" else (ROBUST_DOWN_FRACTIONS[0] * lm,
                                                                        ROBUST_DOWN_FRACTIONS[1])}
            g_grid = float(G.growth_rate(grid0, law, direction).max())
            horizon = int(min(6000, 2.5 * np.log(A_LIST[-1]) / (0.8 * g_grid) + 50))
            y = G.sample_y(mdl, rho, tau, m0, (1000, horizon), rng)
            for det, (gr, w) in dets.items():
                fp = first_passages(y, gr, w, direction, A_LIST)
                g_det = float(G.growth_rate(np.asarray(gr), law, direction).max())
                for i, a_ in enumerate(A_LIST):
                    d = fp[:, i][fp[:, i] >= 0] + 1.0
                    rows.append(dict(scenario=nm, direction=direction, ratio=rho, detector=det, A=a_,
                                     mean_delay=d.mean(), se=d.std(ddof=1) / np.sqrt(d.size), missing=float(np.mean(fp[:, i] < 0)),
                                     kl_inf=kl_inf, kl=kl, lam_oracle=lam_o, growth_detector=g_det,
                                     lower_bound_first_order=np.log(a_) / kl_inf,
                                     delay_over_bound=d.mean() / (np.log(a_) / kl_inf),
                                     limit_fixed_grid=kl_inf / g_det, kl_inf_over_kl=kl_inf / kl))
            print(nm, direction, rho, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e11e_asymptotic_optimality.csv", index=False)
    piv = df.pivot_table(index=["scenario", "direction", "ratio", "detector"], columns="A", values="delay_over_bound").round(3)
    print(piv.to_string())
    return df


def part_e2():
    """Second-order reading of part E, and a cross-check with E10 (Q10e).

    1. For every part-E cell, Monte Carlo renewal constants (kappa, zeta) of the oracle bet give the
       predicted ratio delay / (log A / KL_inf) = 1 + (kappa - zeta) / log A, compared with the
       simulated ratio; and the excess of the continuous mixture and of the default grid over the
       oracle bet in nats, (delay - oracle delay) * KL_inf, which theory says grows like
       (1/2) log log A + O(1) (mixture, interior optimum) or linearly in log A (fixed grid that
       misses the optimum).
    2. E10 observed delays of about 0.6 times the step at which the cumulative KL reaches ln A. In an
       i.i.d. Gaussian version of E10's steady state (variance ratio s1/s0 at t = 3000 from
       results/e10e_transient.csv, true baseline, increase side, c = 3, A = 28,800), compare
       ln A / KL with the second-order prediction and with simulated mean and median delays."""
    rng = np.random.default_rng(SEED + 14)
    e = pd.read_csv(OUT / "e11e_asymptotic_optimality.csv")
    models = {"gaussian": MODELS["gaussian"], "student t3": G.ReturnModel("student_t", df=3.0), "jumps": MODELS["jumps"]}
    mix = "continuous mixture (100 bets, uniform on (0, cap))"
    rows = []
    for (nm, direction, rho), grp in e.groupby(["scenario", "direction", "ratio"], sort=False):
        mdl = models[nm]
        tau, m0 = G.baseline(mdl, 3.0)
        lam_o = float(grp.lam_oracle.iloc[0])
        kl_inf = float(grp.kl_inf.iloc[0])
        rc = G.renewal_constants(lam_o, draw_fn(mdl, rho, tau, m0), direction, 4000, rng=rng)
        piv = grp.pivot_table(index="A", columns="detector", values="mean_delay")
        for a_, r in piv.iterrows():
            la = np.log(a_)
            rows.append(dict(scenario=nm, direction=direction, ratio=rho, A=a_, kappa=rc["kappa"], zeta=rc["zeta"],
                             zeta_se=rc["zeta_se"], oracle_ratio_simulated=r["oracle bet"] / (la / kl_inf),
                             oracle_ratio_second_order=max(la + rc["kappa"] - rc["zeta"], 0.0) / la,
                             mixture_excess_nats=(r[mix] - r["oracle bet"]) * kl_inf,
                             default_excess_nats=(r["default grid"] - r["oracle bet"]) * kl_inf,
                             half_loglog_A=0.5 * np.log(la)))
        print(nm, direction, rho, "renewal constants done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e11e2_second_order_check.csv", index=False)
    print(df.round(3).to_string())

    tr = pd.read_csv(OUT / "e10e_transient.csv")
    A10 = 28800.0
    model = MODELS["gaussian"]
    tau, m0 = G.baseline(model, 3.0)
    xs = []
    for k in (2.0, 4.0):
        q = float(tr[(tr.break_kind == "scale") & (tr.k == k)].mean_ratio_t3000.iloc[0])
        rho = float(np.sqrt(q))
        law = G.y_law(model, rho, tau, m0)
        lam_o, kl_inf = G.best_bet(law, "up")
        kl = float(G.kl_gaussian_scale(rho))
        _, s2 = G.log_factor_moments(lam_o, law, "up")
        rc = G.renewal_constants(lam_o, draw_fn(model, rho, tau, m0), "up", 3000, rng=rng)
        grid0, w0 = default_grid("up", tau**2 / m0)
        g_grid = float(G.growth_rate(grid0, law, "up").max())
        horizon = int(6 * np.log(A10) / g_grid)
        for det, (gr, w) in (("oracle bet", ([lam_o], [1.0])), ("default grid", (grid0, w0))):
            d, fa, miss = simulate_delays(model, rho, tau, m0, gr, w, "up", A10, 0, 1000, horizon, rng, batch=250)
            xs.append(dict(k=k, variance_ratio=q, sd_ratio=rho, detector=det, kl=kl, kl_inf=kl_inf,
                           lnA_over_kl=np.log(A10) / kl, lnA_over_kl_inf=np.log(A10) / kl_inf,
                           second_order_oracle=G.delay_second_order(A10, 1.0, kl_inf, rc["kappa"], rc["zeta"]),
                           zeta=rc["zeta"], zeta_brownian=G.zeta_brownian(kl_inf, s2), kappa=rc["kappa"],
                           mean_delay=d.mean(), median_delay=float(np.median(d)), se=d.std(ddof=1) / np.sqrt(d.size),
                           missing=miss, mean_over_lnA_kl=d.mean() / (np.log(A10) / kl),
                           median_over_lnA_kl=float(np.median(d)) / (np.log(A10) / kl)))
            print(xs[-1], flush=True)
    pd.DataFrame(xs).to_csv(OUT / "e11e2_e10_crosscheck.csv", index=False)
    return df


# ----------------------------------------------------------------------------
# Figures (from the CSV files written by parts A, B and D)
# ----------------------------------------------------------------------------

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical slots 1 to 4
INK, MUTED, GRIDC = "#0b0b0b", "#52514e", "#e4e3df"


def _style(ax):
    ax.grid(True, color=GRIDC, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def part_fig():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(exist_ok=True)
    # 1. prediction quality
    a = pd.read_csv(OUT / "e11a_delay.csv")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    specs = [("first order log(A)/g, grid", a.pred_first_order_grid / a.sim_nu0_grid_mean, a.sim_nu0_grid_mean, "o"),
             ("Wald bound (proved), grid", a.pred_wald_bound_grid / a.sim_nu0_grid_mean, a.sim_nu0_grid_mean, "s"),
             ("second order, single best bet", a.pred_second_order_oracle / a.sim_nu0_oracle_mean, a.sim_nu0_oracle_mean, "D")]
    for (lab, ratio, x, mk), col in zip(specs, SERIES):
        ax.scatter(x, ratio, s=28, marker=mk, color=col, edgecolor="white", linewidth=0.8, label=lab, zorder=3)
    ax.axhline(1.0, color=INK, lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("simulated mean delay (steps, change at the first monitored step)", color=INK, fontsize=9)
    ax.set_ylabel("predicted / simulated delay", color=INK, fontsize=9)
    ax.set_title("E11a: the first-order delay misses the Shiryaev-Roberts memory term (synthetic, A = 7800)",
                 fontsize=9.5, color=INK, loc="left")
    _style(ax)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "e11_delay_predictions.png", dpi=140)
    plt.close(fig)

    # 2. price of robustness against c
    c = pd.read_csv(OUT / "e11b_c_curve.csv")
    c = c[np.isfinite(c.c) & (c.scenario == "gaussian")]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), sharey=True)
    for ax, direction, ratios, title in ((axes[0], "up", (1.25, 1.5, 2.0), "Increases (sd ratio)"),
                                          (axes[1], "down", (0.5, 0.67, 0.8, 0.9), "Decreases (sd ratio)")):
        for rho, col in zip(ratios, SERIES):
            s_ = c[(c.direction == direction) & np.isclose(c.ratio, rho)].sort_values("c")
            ax.plot(s_.c, s_.eff_grid, color=col, lw=2, label=f"{rho:g}")
        ax.set_xlabel("truncation constant c (tau = c x scaled MAD)", color=INK, fontsize=9)
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.set_ylim(0, 1.02)
        _style(ax)
        ax.legend(fontsize=8, frameon=False, title="sd ratio", title_fontsize=8)
    axes[0].set_ylabel("efficiency g_grid / KL (Gaussian LR = 1)", color=INK, fontsize=9)
    fig.suptitle("E11b: growth of the robust detector relative to the Gaussian likelihood ratio, clean Gaussian data",
                 fontsize=9.5, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "e11_price_of_robustness.png", dpi=140)
    plt.close(fig)

    # 3. slack in the false alarm bound
    d = pd.read_csv(OUT / "e11d_decomposition.csv")
    dp = pd.read_csv(OUT / "e11d_worst_case_dp.csv")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for (lab, sub), col in zip((("two-sided default grids, Gaussian null (E9)", d[d.detector == "two-sided default grids (E9)"]),
                                ("increase side, default grid, Gaussian null", d[d.detector == "increase side only, default grid"]),
                                ("increase side, single bet 0.5, Gaussian null", d[d.detector == "increase side only, single bet 0.5"])),
                               SERIES):
        ax.plot(sub.horizon, sub.rate_over_bound, color=col, lw=2, marker="o", ms=5, label=lab)
    w = dp[(dp.direction == "up") & np.isclose(dp.lam, 0.5)]
    ax.plot(list(w.horizon) + [A], list(w.lower_over_bound) + [1.0], color=SERIES[3], lw=2, ls="--", marker="s", ms=5,
            label="single bet 0.5, worst null (DP lower bound; exact at m = A)")
    ax.axhline(1.0, color=INK, lw=1)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("horizon m (steps)", color=INK, fontsize=9)
    ax.set_ylabel("P(false alarm within m) / (m / A)", color=INK, fontsize=9)
    ax.set_title("E11d: the bound m/A is nearly attained by an adversarial null, not by Gaussian data",
                 fontsize=9.5, color=INK, loc="left")
    _style(ax)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "e11_false_alarm_slack.png", dpi=140)
    plt.close(fig)

    # 4. first-order optimality as A grows (tracker B73)
    e = pd.read_csv(OUT / "e11e_asymptotic_optimality.csv")
    e2 = pd.read_csv(OUT / "e11e2_second_order_check.csv")
    cells = [("up", 1.25), ("up", 1.5), ("up", 2.0), ("down", 0.8), ("down", 0.5)]
    fig, axes = plt.subplots(1, 5, figsize=(12.5, 3.4), sharey=True)
    mix = "continuous mixture (100 bets, uniform on (0, cap))"
    for ax, (direction, rho) in zip(axes, cells):
        sub = e[(e.scenario == "gaussian") & (e.direction == direction) & np.isclose(e.ratio, rho)]
        pr = e2[(e2.scenario == "gaussian") & (e2.direction == direction) & np.isclose(e2.ratio, rho)].sort_values("A")
        for (det, lab), col in zip((("oracle bet", "oracle bet"), (mix, "fine mixture"), ("default grid", "default grid")), SERIES):
            s_ = sub[sub.detector == det].sort_values("A")
            ax.plot(s_.A, s_.delay_over_bound, color=col, lw=2, marker="o", ms=4, label=lab)
        ax.plot(pr.A, pr.oracle_ratio_second_order, color=SERIES[0], lw=1.5, ls="--", label="oracle, second-order formula")
        ax.axhline(1.0, color=INK, lw=1)
        ax.set_xscale("log")
        ax.set_ylim(0, 1.25)
        ax.set_title(f"{'increase' if direction == 'up' else 'decrease'}, sd ratio {rho:g}", fontsize=9, color=INK, loc="left")
        ax.set_xlabel("alarm level A", color=INK, fontsize=8.5)
        _style(ax)
    axes[0].set_ylabel("mean delay / (log A / KL_inf)", color=INK, fontsize=9)
    axes[0].legend(fontsize=7.5, frameon=False, loc="lower right")
    fig.suptitle("E11e: delay relative to the nonparametric lower bound log(A)/KL_inf, Gaussian data, change at step 1 (synthetic)",
                 fontsize=9.5, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "e11_asymptotic_optimality.png", dpi=140)
    plt.close(fig)


PARTS = {"A": part_a, "B": part_b, "C": part_c, "C2": part_c2, "C3": part_c3, "C4": part_c4, "C4V": part_c4v, "D": part_d, "D2": part_d2, "E": part_e, "E2": part_e2, "FIG": part_fig}

if __name__ == "__main__":
    which = sys.argv[1:] or list(PARTS)
    for p in which:
        print(f"=== part {p} ===", flush=True)
        PARTS[p]()
