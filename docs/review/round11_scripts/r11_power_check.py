"""Review round 11 (E25), check 4a: the power simulation (pre-registration section 9) and gate 5.
(1) Independent strata on F24a's planning calendar (E24.planning_events supplies the calendar only) against
    results/e25_power_calendar.csv: n_i, balance, entry times.
(2) Exact reproduction of results/e25_power.csv with E25's own simulate_power and seed 25,000 (nothing is written).
(3) An independent simulation (own code, seed 1111) of the same model: differences in units of the Monte Carlo s.e.
(4) Exact expected log e-values (own code) against results/e25_power_growth.csv.
(5) The dates at which E25's crossing probability first reaches F24a's power (results/e25_power_curve.csv).
(6) Gate 5: one law rerun with E25's simulate_h0 and the stated seed, against results/e25_validity.csv.
Usage: nice -n 10 python docs/review/round11_scripts/r11_power_check.py"""

from __future__ import annotations

import sys
from datetime import datetime, timezone
from math import comb, log
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT / "src"))
import r11_strata_lib as L  # noqa: E402

UTC = timezone.utc


def own_planning_table():
    import e24_blind_guard as E24  # calendar only

    ev = E24.planning_events("2029-01-01")
    rels = [(d.strftime("%Y-%m-%d"), t, h) for d, t, h in zip(ev.date, ev.event_type, ev.release_time_et)]
    st = L.build_strata(rels, [r[0] for r in rels], datetime(2026, 10, 25, tzinfo=UTC), None, amended=True)
    rows = pd.read_csv(ROOT / "data/events/us_macro_events_e24_future.csv", dtype=str, keep_default_na=False)
    settle = [(L.to_date(d), L.et_to_utc(L.to_date(d), h)) for d, h in zip(rows.date, rows.release_time_et)]
    settle += [(L.to_date(d), L.et_to_utc(L.to_date(d), h)) for d, h, k in
               zip(ev.date.dt.strftime("%Y-%m-%d"), ev.release_time_et, ev.kind) if k != "event file"]
    run = None
    for s in st:
        e = max([s["data_needed_until"]] + [t for d, t in settle if abs((d - s["date"]).days) <= 30])
        run = e if run is None else max(run, e)
        s["enters"] = run
    return st


def own_sim(st, rates, reps, seed):
    rng = np.random.default_rng(seed)
    n = np.array([s["n"] for s in st])
    typ = [s["type"] for s in st]
    p1 = np.array([rates[t][0] for t in typ])
    p0 = np.array([rates[t][1] for t in typ])
    J = len(st)
    first = np.zeros(reps, int)
    for b0 in range(0, reps, 2000):
        b = min(2000, reps - b0)
        z = (rng.random((b, J)) < p1).astype(float)
        k = rng.binomial(n[None, :].repeat(b, 0), p0[None, :].repeat(b, 0)).astype(float)
        p = (z + k) / (n + 1)
        m = np.zeros((b, J))
        for ps in L.PSI:
            m += np.cumprod(ps**z / (1 - p + ps * p), axis=1) / len(L.PSI)
        hit = m >= 20
        first[b0:b0 + b] = np.where(hit.any(1), hit.argmax(1) + 1, 0)
    return first


def main():
    import e24_blind_guard as E24
    import e25_sequential_events as E25

    out = []
    st = own_planning_table()
    ref = pd.read_csv(ROOT / "results/e25_power_calendar.csv")
    same_n = [s["n"] for s in st] == ref.n.tolist()
    same_bal = [s["before"] for s in st] == ref.controls_before.tolist() and [s["after"] for s in st] == ref.controls_after.tolist()
    ent_ref = pd.to_datetime(ref.enters_test, utc=True)
    same_ent = all(pd.Timestamp(s["enters"]) == e for s, e in zip(st, ent_ref))
    out.append(f"(1) planning calendar: {len(st)} strata (CSV {len(ref)}); n_i equal: {same_n}; before/after equal: "
               f"{same_bal}; entry times equal: {same_ent}")
    out.append(f"    n_i counts: {np.bincount([s['n'] for s in st], minlength=9).tolist()}; balanced "
               f"{sum(1 for s in st if s['before'] == s['after'] > 0)}; one side only "
               f"{sum(1 for s in st if (s['before'] == 0) != (s['after'] == 0))}")
    for d in ("2027-06-01", "2028-01-01", "2029-01-01"):
        cut = datetime.fromisoformat(d).replace(tzinfo=UTC)
        out.append(f"    strata entered by {d}: {sum(1 for s in st if s['enters'] <= cut)}")
    # (2) exact reproduction with E25's own simulator
    tab = E25.planning_table()
    inp = pd.read_csv(ROOT / "results/e24_future_planning_inputs.csv")
    scen = E24.planning_scenarios(inp, "S2")
    pw = pd.read_csv(ROOT / "results/e25_power.csv")
    entry_ns = E25._ns(tab["enters_test"])
    for sname in E25.SCENARIOS:
        res = E25.simulate_power(tab, scen[sname])
        r = pw[(pw.scenario == sname) & pw.procedure.str.startswith("E25")].iloc[0]
        f = res["first"]
        vals = []
        for d in ("2027-06-01", "2028-01-01", "2029-01-01"):
            cd = np.where(f > 0, entry_ns[np.maximum(f, 1) - 1], np.iinfo(np.int64).max)
            vals.append(float((cd <= pd.Timestamp(d, tz="UTC").value).mean()))
        refv = [r[f"P_cross_by_{d}"] for d in ("2027-06-01", "2028-01-01", "2029-01-01")]
        out.append(f"(2) {sname}: rerun {vals} CSV {refv} identical: {np.allclose(vals, refv, rtol=0, atol=0)}")
        # (3) independent simulation
        fo = own_sim(st, scen[sname], 20000, 1111)
        ent = np.array([pd.Timestamp(s["enters"]).value for s in st])
        cd = np.where(fo > 0, ent[np.maximum(fo, 1) - 1], np.iinfo(np.int64).max)
        mine = [float((cd <= pd.Timestamp(d, tz="UTC").value).mean()) for d in ("2027-06-01", "2028-01-01", "2029-01-01")]
        z = [(a - b) / np.sqrt(2 * max(b * (1 - b), 1e-9) / 20000) for a, b in zip(mine, refv)]
        out.append(f"    independent (seed 1111): {np.round(mine, 4).tolist()}; difference in s.e. of the difference: "
                   f"{np.round(z, 2).tolist()}")
    # (4) exact expected log growth, own code, psi = 8 pooled and psi = 4 2026 rates, over strata entered by 2029
    gr = pd.read_csv(ROOT / "results/e25_power_growth.csv")
    cut = datetime(2029, 1, 1, tzinfo=UTC)
    sub = [s for s in st if s["enters"] <= cut]
    for sname, ps in (("pooled 2024 to 2026 (planning)", 8.0), ("2026 only (E23 design)", 4.0),
                      ("half of the pooled excess", 4.0)):
        rates = scen[sname]
        g = []
        for s in sub:
            p1, p0 = rates[s["type"]]
            n = s["n"]
            tot = 0.0
            for z in (0, 1):
                pz = p1 if z else 1 - p1
                for k in range(n + 1):
                    tot += pz * comb(n, k) * p0**k * (1 - p0) ** (n - k) * log(L.e_val(z, z + k, n, ps))
            g.append(tot)
        r = gr[(gr.scenario == sname) & gr.procedure.str.startswith("E25") & (gr.psi == ps)].iloc[0]
        out.append(f"(4) {sname}, psi {ps}: own mean E[log e] {np.mean(g):.6f} over {len(g)} strata, CSV "
                   f"{r.mean_expected_log_e_per_release:.6f}")
    # (5) the date E25 first reaches F24a's power
    cur = pd.read_csv(ROOT / "results/e25_power_curve.csv")
    for sname in E25.SCENARIOS[:3]:
        f24 = pw[(pw.scenario == sname) & pw.procedure.str.startswith("E25")]["F24a_power_at_2027-06-01"].iloc[0]
        c = cur[(cur.scenario == sname) & cur.procedure.str.startswith("E25")]
        hit = c[c.P_cross_by_date >= f24].iloc[0]
        out.append(f"(5) {sname}: first date with P >= F24a's {f24:.4f}: {hit.date} ({hit.P_cross_by_date:.4f})")
    # (6) gate 5, one law, with the stated seed
    cal = E25.planning_table()
    ns = np.resize(cal["n"].to_numpy(int), 200)
    val = pd.read_csv(ROOT / "results/e25_validity.csv")
    for i, (law, rate) in ((7, ("S chosen from the earlier strata", np.nan)), (0, ("exchangeable Bernoulli", 0.05))):
        res = E25.simulate_h0(ns, law, rate, 20000, np.random.default_rng([25000, i]))
        r = val.iloc[i]
        out.append(f"(6) gate 5 law {i} ({law}): rerun P(M ever >= 20) {res['share_M_ever_ge_20']:.5f}, CSV "
                   f"{r.share_M_ever_ge_20:.5f}; identical: {res['share_M_ever_ge_20'] == r.share_M_ever_ge_20}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
