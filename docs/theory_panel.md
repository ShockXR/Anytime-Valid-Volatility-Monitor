# Monitoring a panel of assets with false discovery control

This note extends the single-asset monitor (`docs/math.md`, sections 5 and 6) to K assets watched at the same time. It states the problem precisely, explains why the multiple testing tool for e-values (e-BH) cannot be applied to the Shiryaev-Roberts statistic directly, builds e-values that it can be applied to, proves what the resulting panel monitor controls, and says plainly what it does not control. Code: `src/avmonitor/panel.py`. Experiments: E15 (`docs/report_parts/E15.md`).

Every theorem below is proved here, except where a result is explicitly cited. Section 9 lists which results are standard, which are cited, and how each reference was checked.

**Short version for a reader in a hurry.**

| Rule (one block, for example one day) | What is guaranteed, under any dependence between assets | What is not guaranteed |
|---|---|---|
| Naive: each asset alarms at its own level (e-value 1/α, or Shiryaev-Roberts at A = m/α) | each asset: P(false alarm in the block) ≤ α | the panel: P(at least one false alarm) can reach Kα (25% for five assets at 5%) |
| Bonferroni: e-value K/α (Shiryaev-Roberts at A·K) | P(any false discovery in the block) ≤ α, at all times of the block at once | power: every asset pays log K |
| **Sticky e-BH on mixture e-processes (Theorem 1)** | false discovery rate ≤ α at every stopping time; E[max over the block of the false discovery proportion] ≤ α(1 + log(1/α)); P(any discovery) ≤ α if nothing changes | the probability of at least one false discovery when some assets do change; the expected number of false discoveries (up to αK); the false discovery rate of all discoveries pooled over many blocks |

## 1. The problem

**Assets, clock, information.** Assets k = 1, ..., K are observed on one common clock t = 1, 2, ... (one-minute steps in UTC). At step t asset k may have a return `r_{k,t}` or none (its market is closed, or the minute is missing); `O_{k,t} = 1` if it is observed. `F_t` is everything observed on the whole panel up to step t (the **panel filtration**). This matters: below, "conditional on the past" always means conditional on the past of **all** assets.

**Statistic and null for one asset (as in `docs/math.md`, 5.1).** For the current block (section 6), asset k has a truncation level `τ_k` and null values `m_k^↑ ≥ m_k^↓ > 0`, all fixed before the block starts (they come from data strictly before the block). When asset k is observed, `X_{k,t} = min(r_{k,t}², τ_k²) ∈ [0, τ_k²]`, and

`μ_{k,t} = E[X_{k,t} | F_{t-1}]` (on `O_{k,t} = 1`), the conditional truncated variance given the whole panel's past.

**Assumption (M), observation.** Whether `r_{k,t}` is observed is known at step t - 1 (`O_{k,t}` is `F_{t-1}`-measurable). Exchange calendars and trading hours satisfy this. (Random data gaps can be handled as in `docs/math.md` 6.7 by assuming that, given the past, being observed says nothing about the size of the return; every proof below goes through with the conditional means taken given `O_{k,t} = 1`.)

**Change time.** The null for asset k at step s says `m_k^↓ ≤ μ_{k,s} ≤ m_k^↑`. Let

`ν_k = inf{ s ≥ s_0 : O_{k,s} = 1 and μ_{k,s} ∉ [m_k^↓, m_k^↑] }` (∞ if this never happens),

the first step of the block (which starts at `s_0`) at which asset k's null fails. Because `μ_{k,s}` and `O_{k,s}` are known at step s - 1, the event `{ν_k = s}` is known at s - 1: `ν_k` is a **predictable** time, and `ν_k - 1` is a stopping time. Under the model of `docs/math.md` 5.1 (returns `σ_t ε_t` with a fixed innovation law), `ν_k` is the first step at which the volatility of asset k leaves its calibration level. Nothing is assumed about how assets depend on each other.

**Discoveries and errors.** A panel monitor declares, at random times, that some assets have changed. Let `τ_k^D` be the step at which asset k is discovered (∞ if never) and `D_t = {k : τ_k^D ≤ t}`. A discovery of asset k is **false** if it is made before the asset changed, `τ_k^D < ν_k`. With `V_t` the number of false discoveries in `D_t` and `R_t = |D_t|`:

* false discovery proportion `FDP_t = V_t / max(R_t, 1)`;
* false discovery rate at a stopping time σ: `FDR(σ) = E[FDP_σ]` (Benjamini and Hochberg 1995 at a fixed time);
* family-wise error: P(V ≥ 1);
* the time-uniform quantity `E[sup_t FDP_t]`.

Indexing the nulls by time is what makes this a change detection problem and not only a test of "no change ever": an asset that changes late in the block is a true null until it changes, and a discovery before that is counted as false.

## 2. e-values and e-BH

An **e-value** for a null is a nonnegative random variable e with `E[e] ≤ 1` under the null. Large values are evidence against it (Markov: `P(e ≥ 1/α) ≤ α`). The **e-BH procedure** of Wang and Ramdas (2022) takes e-values `e_1, ..., e_K` and rejects the R largest, where

`R = max{ r : e_(r) ≥ K / (α r) }`, `e_(r)` the r-th largest (R = 0 if no r qualifies).

A set S of hypotheses is **self-consistent** if every k in S has `e_k ≥ K / (α |S|)`. The e-BH set is self-consistent, and it contains every self-consistent set (if S is self-consistent, the |S| largest e-values are all at least `K/(α|S|)`, so R ≥ |S| and every member of S clears `K / (α R)`).

**Lemma 1 (Wang and Ramdas 2022).** If S is self-consistent, then `FDP(S) ≤ (α / K) Σ_{k null} e_k`. Hence, whatever the dependence between the e-values, `FDR ≤ α K_0 / K` (K_0 the number of true nulls).

*Proof.* For a null k in S, `1 / |S| ≤ α e_k / K`. Sum over the nulls in S: `V / |S| ≤ (α / K) Σ_{k null in S} e_k ≤ (α / K) Σ_{k null} e_k`. Take expectations. ∎

Two features matter here. First, only `E[e_k] ≤ 1` for each null k is used, nothing about the joint law: e-BH needs no independence or positive dependence. (The p-value BH procedure needs positive regression dependence or a log K correction under arbitrary dependence; Benjamini and Yekutieli 2001.) Second, the e-values must be e-values **at the moment e-BH is applied**. That is the difficulty in monitoring.

## 3. Why a Shiryaev-Roberts statistic is not an e-value

The single-asset detector (`docs/math.md`, 5.4) uses `M_t = Σ_{j ≤ t} L_t^{(j)}`, where `L^{(j)}` is the wealth of a bettor who starts betting at step j. Under the null `M_t - t` is a supermartingale, so `E[M_t] ≤ t` and, by optional stopping, `E[M_σ] ≤ E[σ]` for stopping times. This is what makes M an **e-detector** (Shin, Ramdas and Rinaldo 2024): it controls the average run length and the false alarm rate per unit time. It is not an e-value: its expectation grows with time.

Dividing by elapsed time does not fix this.

**Proposition 1.** Suppose `E[f_t | F_{t-1}] = 1` for every t (the null holds with equality), and `P(f_1 < 1) > 0`. Let σ = 1 if `f_1 ≥ 1` and σ = 2 otherwise. Then `E[M_σ / σ] = 1 + E[(1 - f_1)⁺] / 2 > 1`. So `M_t / t`, which has expectation at most 1 at every fixed t, is not an e-value at a stopping time.

*Proof.* `M_1 = f_1` and `M_2 = f_2 (f_1 + 1)`. Since σ = 2 only on `{f_1 < 1}`, which is known at step 1, `E[M_σ/σ] = E[f_1 1{f_1 ≥ 1}] + E[1{f_1 < 1} (f_1 + 1) E[f_2 | F_1]] / 2 = E[f_1 1{f_1 ≥ 1}] + E[1{f_1 < 1} (f_1 + 1)] / 2`. Write `f_1 = f_1 1{f_1 ≥ 1} + f_1 1{f_1 < 1}` and use `E[f_1] = 1`: the expression equals `1 + E[1{f_1 < 1} (1 - f_1)] / 2`. ∎

(`tests/test_panel.py::test_shiryaev_roberts_divided_by_time_is_not_an_e_value_at_a_stopping_time` checks this by simulation for the bets of `detectors.py`.)

At a fixed time t, `M_t / t` is an e-value, so e-BH could be applied once a day at a fixed clock time. But a monitor should act when the evidence arrives. The deeper obstacle is a theorem of Dandapanthula and Ramdas (2025): any multi-stream procedure with a finite average run length under the global null has worst-case false discovery rate (and family-wise error) equal to one at some stopping time. Their answer is to divide the error at a stopping time by the expected waiting time ("error over patience") and to control that with e-BH applied to e-detectors. This note takes the other route: give up a finite average run length over an unbounded horizon, and recover genuine false discovery control at every stopping time, either with an e-process that pays more for late changes (section 4 and Proposition 2) or with blocks that restart on a fixed schedule (section 6). Both routes are legitimate; they control different things.

## 4. e-processes for "no change since the block started"

Fix one asset and drop k. Index its observed returns in the block by its own clock n = 1, 2, .... For a bet λ from the increase grid, the factor is `f_n = 1 + λ (X_n / m^↑ - 1)`; for the decrease grid, `f_n = 1 + λ (1 - X_n / m^↓)` with `λ < 1 / (B - 1)`, `B = τ² / m^↓` (if `B ≤ 1` the decrease side is switched off). Both factors are positive, and under the null their conditional mean is at most 1 (`docs/math.md`, 5.2). `L_n^{(j)} = Π_{i=j}^{n} f_i` is the wealth of a bettor who starts at own observation j.

Choose **prior weights over the change time**, `w_j ≥ 0` with `Σ_j w_j = 1`, and define for each bet

`S_n = Σ_{j ≤ n} w_j L_n^{(j)}`, which satisfies `S_n = f_n (S_{n-1} + w_n)`, `S_0 = 0`,

and the **mixture e-process** (average over the bets of each side, then over the two sides)

`E_n = ½ (mean over the increase grid of S_n + mean over the decrease grid of S_n) + T_n`, with `T_n = Σ_{j > n} w_j`.

The term `T_n` is the weight of bettors who have not started yet; they still hold their stake. On the panel clock, `E` stays constant at steps where the asset is not observed.

**Lemma 2.** Let `Z_t = E_{t ∧ (ν - 1)}` (the e-process on the panel clock, stopped just before the change). Under assumption (M), Z is a nonnegative supermartingale for the panel filtration with `Z_0 = 1` (a *test supermartingale*). In particular `E[Z_σ] ≤ 1` for every stopping time σ.

*Proof.* `E_0 = T_0 = 1` and `E ≥ 0`. Fix t. Because ν is predictable, the three events `{ν ≤ t - 1}`, `{ν = t}` and `{ν > t}` are all known at step t - 1. On the first two, `Z_t = Z_{t-1}` (the process is stopped at `ν - 1 ≤ t - 1`). On `{ν > t}` the null holds at step t. If the asset is not observed at t, nothing moves. If it is, say as its n-th observation (n is known at t - 1 by assumption (M)), then for each bet `E[S_n | F_{t-1}] = E[f_n | F_{t-1}] (S_{n-1} + w_n) ≤ S_{n-1} + w_n`, while `T_n = T_{n-1} - w_n`. Averaging over bets and sides gives `E[Z_t | F_{t-1}] ≤ Z_{t-1}`. ∎

So `E_σ` is an e-value at every stopping time σ of the panel, as long as the asset has not changed. This is the object e-BH needs.

**Two choices of weights** (`panel.change_weights`).

* *Window*: `w_j = 1/h` for `j ≤ h`, 0 after. For `n ≤ h`, `E_n = M_n / h + (h - n) / h`, where `M = (M^↑ + M^↓)/2` is exactly the two-sided Shiryaev-Roberts statistic of `docs/math.md` 5.4 (checked to 10 digits in `tests/test_panel.py`). This is the supermartingale `S_t` used to prove Theorem B there, now used as an e-process. Consequence: an e-value threshold c corresponds to a Shiryaev-Roberts threshold of about `c h`; the naive rule `E ≥ 1/α` with h = one block of m steps is the single-asset monitor with `A = m / α`, and `E ≥ K/α` is that monitor at `A·K` (Bonferroni). Evidence from a change after the first h observations is not added (later bettors get weight 0), which is why the window is used with blocks.
* *Pareto*: `w_j = h / ((j + h - 1)(j + h))`, with tail `T_n = h / (n + h)`. The weights are about 1/h for j up to h and decay like `h / j²` after. This covers an unbounded horizon.

**The price of a late change.** If the asset changes at own observation ν, then `E_n ≥ w_ν L_n^{(ν)} / (2G)` for every bet of the right side (G the grid size), so the alarm level c is reached once `log L_n^{(ν)} ≥ log c + log(1/w_ν) + log(2G)`. The post-change growth rate of `log L` decides how fast that happens. For the window this is `log(c h)` (the Shiryaev-Roberts cost at `A = c h`). For Pareto weights `log(1/w_ν) ≈ log h + 2 log(1 + ν/h)`: a change ten blocks in (ν = 10h) costs about `2 log 11 ≈ 4.8` nats more than a change right at the start. The next result shows that some such price is unavoidable.

**Proposition 2 (no e-process detects late changes with bounded delay).** Take one asset with i.i.d. observations from a law P (a member of the null) and, for a change at ν, observations from a law Q from step ν on. Let `D = KL(Q ‖ P) < ∞`. Let E be any e-process for P (`E_P[E_σ] ≤ 1` for all stopping times σ of the data), let c > 1 and `τ_c = inf{t : E_t ≥ c}`. Suppose there are an integer d ≥ 1 and β < 1 with `(1 - β) c > 1` such that, for every change time `ν ∈ {1, ..., N}`, `Q_ν(ν ≤ τ_c < ν + d) ≥ 1 - β`. Then

`d · D ≥ (1 - β) log(c ⌊N / d⌋) - log 2`.

*Proof.* First, `P(τ_c < ∞) ≤ 1/c`: on `{τ_c ≤ n}`, `E_{τ_c ∧ n} ≥ c`, so `P(τ_c ≤ n) ≤ E_P[E_{τ_c ∧ n}] / c ≤ 1/c`; let n grow. Take the J = ⌊N/d⌋ disjoint windows starting at `ν_i = 1 + (i - 1) d` and the events `A_i = {ν_i ≤ τ_c < ν_i + d}`, which depend only on the first `ν_i + d - 1` observations. Their P-probabilities `p_i` add up to at most 1/c. Under `Q_{ν_i}` and under P the first `ν_i + d - 1` observations differ only in the last d, so their Kullback-Leibler divergence is `d D`, and by the data processing inequality the binary divergence satisfies `kl(Q_{ν_i}(A_i), p_i) ≤ d D`. If `p_i ≥ 1 - β` for some i, then `1 - β ≤ 1/c`, which is excluded. So `p_i < 1 - β ≤ Q_{ν_i}(A_i)`; `kl(x, p)` increases in x for x ≥ p, and `kl(x, p) ≥ x log(1/p) - log 2`. Hence `(1 - β) log(1/p_i) - log 2 ≤ d D`, that is `p_i ≥ exp(-(d D + log 2)/(1 - β))`. Summing over the J windows, `J exp(-(d D + log 2)/(1 - β)) ≤ 1/c`, which rearranges to the claim. ∎

So the worst-case delay over change times in a horizon of N observations is at least about `(1 - β) log(c N) / D` (up to the `log d` inside): it grows without bound as N grows. The window e-process with `h = N` and level c is the Shiryaev-Roberts detector at `A ≈ c N` (the window identity of this section). With likelihood-ratio bets its delay is about `log(c N) / D` (Pollak 1985; Lai 1998 proves the matching lower bound `log A / D` under an average run length constraint), so the order `log N` is attained. With the nonparametric bets of this project, D is replaced by the growth rate of the bets, which is smaller (E11). In words: **to have an e-value at every stopping time, the monitor must treat its whole horizon as if it were the average run length.** This is the quantitative counterpart, for one stream, of the impossibility in Dandapanthula and Ramdas (2025). The argument is the standard change-of-measure argument of the change detection lower bounds (Lorden 1971, Lai 1998), written out for e-processes; we make no claim that it is new.

## 5. The panel monitor and its guarantee

**Sticky e-BH.** Within a block, give each asset its mixture e-process `E_k`. At every step t:

1. stop the e-process of every asset already discovered at its value at discovery: `Ẽ_k(t) = E_k(t ∧ τ_k^D)`;
2. compute the e-BH set of `Ẽ(t)` at level α; its new members are discovered at t.

(`panel.run_block` and `panel.PanelMonitor` with `rule="ebh"`. The code freezes a discovered asset, which is the same thing.)

**Theorem 1.** Under assumption (M), with block parameters fixed before the block, and whatever the dependence between assets:

1. (*Monotone and self-consistent.*) `D_t` only grows, `D_t = {k : τ_k^D ≤ t}`, and every `D_t` is self-consistent for `Ẽ(t)`.
2. (*False discovery rate at stopping times.*) For every stopping time σ of the panel filtration (bounded or not), `FDR(σ) ≤ α`.
3. (*Time-uniform.*) `P(sup_t FDP_t ≥ x) ≤ α / x` for every x in (0, 1], and `E[sup_t FDP_t] ≤ α (1 + log(1/α))` (0.20 for α = 0.05).
4. (*No change anywhere.*) If no asset changes in the block, the probability of any discovery in the block is at most α.

*Proof.* (1) By induction on t. Suppose `D_s = {k : τ_k^D ≤ s}` and these sets grow for s < t. A k in `D_{t-1}` was discovered at `τ_k^D ≤ t - 1` with `E_k(τ_k^D) ≥ K / (α |D_{τ_k^D}|) ≥ K / (α |D_{t-1}|)`, and its stopped value stays at `E_k(τ_k^D)`. So `D_{t-1}` is self-consistent for `Ẽ(t)` and is contained in the e-BH set (section 2). Hence nothing discovered is ever lost, the new members are exactly the assets with `τ_k^D = t`, and e-BH sets are self-consistent by construction.

(2) Let `Z_k(t) = E_k(t ∧ τ_k^D ∧ (ν_k - 1))`. Since `τ_k^D` and `ν_k - 1` are stopping times of the panel filtration, Lemma 2 and the stopping theorem make each `Z_k` a test supermartingale, and so is the average `Z̄(t) = (1/K) Σ_k Z_k(t)` (all are supermartingales for the same filtration). If k is a false discovery in `D_t`, then `τ_k^D ≤ t` and `τ_k^D ≤ ν_k - 1`, so `Ẽ_k(t) = E_k(τ_k^D) = Z_k(t)`. Self-consistency gives `1 / R_t ≤ α Z_k(t) / K`, and summing over false discoveries,

`FDP_t ≤ (α / K) Σ_k Z_k(t) = α Z̄(t)` for every t, on every path.

For a stopping time σ, `E[FDP_σ] ≤ α E[Z̄_σ] ≤ α`. (For unbounded σ use `σ ∧ n`, Fatou's lemma, and the fact that `FDP_t` converges because `D_t` is monotone.)

(3) `sup_t FDP_t ≤ min(1, α sup_t Z̄_t)`, and Ville's inequality gives `P(sup_t Z̄_t ≥ x / α) ≤ α / x`. Then `E[sup_t FDP_t] = ∫_0^1 P(sup_t FDP_t ≥ x) dx ≤ ∫_0^1 min(1, α/x) dx = α + α log(1/α)`.

(4) If nothing changes, every discovery is false, so "some discovery" means `sup_t FDP_t = 1`, which has probability at most α by (3). ∎

**Remarks.**

1. *What filtration.* The proof needs each `Z_k` to be a supermartingale for the **panel** filtration, because discovery times depend on all assets. This is why the null conditions on the whole panel's past. If assets are independent, an asset's own-past null implies the panel null. If they are not (for example volatility that spills over from one asset to another with a lag), the panel null is stronger. At a fixed, data-independent time only `E[E_k(t)] ≤ 1` is needed, and the own-past null suffices. Wang, Dandapanthula and Ramdas (2025) make this local versus global filtration point for e-BH applied to stopped e-processes, give a causal condition under which local e-processes are global, and show that e-BH with global (compound) e-processes is necessary and sufficient for FDR control at stopping times. Xu, Wang and Ramdas (2021) used e-BH on e-processes at stopping times in a bandit setting.
2. *Why stop discovered assets.* Plain e-BH on the live e-values at every step has the same guarantee (2) at each stopping time, but its set can shrink when evidence fades, so an alarm could be withdrawn. Stopping at discovery makes the set monotone at no cost to the guarantee. Running maxima would not do: e-BH applied to running suprema of e-processes does not control the false discovery rate (Tavyrikov, Goeman and de Heide 2026, who restore control for suprema with "adjusters").
3. *Bonferroni and naive rules on the same e-processes.* Bonferroni (`E_k ≥ K/α`): the probability of any false discovery in the block is at most `Σ_k P(sup_t Z_k ≥ K/α) ≤ α` (Ville per asset, union bound), at all times at once. Naive (`E_k ≥ 1/α`): each asset alone is at α. Path by path, a naive rule flags an asset no later than e-BH (e-BH needs `Ẽ ≥ K/(αR) ≥ 1/α`), and e-BH no later than Bonferroni (an e-value at least K/α is self-consistent on its own). This ordering is checked in the tests.
4. *Estimated baselines.* With the confidence-sequence baseline of `docs/math.md` section 6, each asset's bounds fail with probability at most `δ_k`, on an event `G_k^c` known before the block starts. Replacing `Z_k` by `1{G_k} Z_k` (still a nonnegative supermartingale started at most at 1) and bounding the false discoveries of assets with failed bounds by 1 gives `FDR(σ) ≤ α + Σ_k δ_k` and `P(sup_t FDP_t ≥ x) ≤ α / x + Σ_k δ_k`. `PanelMonitor(baseline="cs", delta=δ)` uses `δ_k = δ / K`, so the bound is `α + δ`.
5. *Fewer than K nulls.* If some assets change at the first step of the block, they never contribute to the sum in (2), and the bound improves to `α K_0 / K` with `K_0` the number of assets that do not change at the first step.
6. *What is not controlled.* The probability of at least one false discovery is not controlled once some assets change (that is the point of a false discovery rate). The number of false discoveries is not either: e-BH needs `Ẽ_k ≥ 1/α` at least, so `E[V] ≤ Σ_k P(sup_t Z_k ≥ 1/α) ≤ α K`, the same worst case as the naive rule. For a target on `E[sup_t FDP_t]` instead of the FDR, run e-BH at the α' that solves `α'(1 + log(1/α')) = α` (α' ≈ 0.0087 for α = 0.05).

## 6. Genuine change detection over a long horizon

Theorem 1 already handles one change per asset: nulls are indexed by time, and a discovery is false only if it precedes the change. Two issues remain for a monitor that runs for months: an asset can change several times, and Proposition 2 says that one e-process per asset for the whole horizon gets slower and slower for late changes.

**(a) One block for the whole horizon.** Use Pareto weights (or a window as long as the horizon) and never restart. Theorem 1 then holds for the whole horizon: FDR at most α at every stopping time, and at most α probability of any discovery if nothing changes, ever. The price is the delay of Proposition 2, growing like `2 log(ν/h)` for a change at ν. After an asset is discovered it is frozen; recalibrating it and starting a new hypothesis for it would add hypotheses over time, and Theorem 1 would no longer cover the whole set (this is the territory of online multiple testing, where hypotheses arrive over time; for e-values see Xu and Ramdas 2024).

**(b) Blocks on a fixed schedule (the construction used in E15).** Fix block start times `s_1 < s_2 < ...` in advance (for example every UTC midnight). At each `s_b`, recalibrate every asset on its own last `cal_len` finite returns (data strictly before `s_b`, so the block parameters are known at `s_b - 1`), start fresh e-processes with window weights over one block, and run sticky e-BH until the block ends. Theorem 1 applies to every block, conditionally on `F_{s_b - 1}`. Hence:

* for every block b and every stopping time σ in it, `E[FDP_σ^{(b)} | F_{s_b - 1}] ≤ α`, and `E[sup_{t in b} FDP_t | F_{s_b - 1}] ≤ α (1 + log(1/α))`;
* the average over blocks of the per-block false discovery rate is at most α;
* if no asset changes in any block, the expected number of blocks with a discovery, among the first B, is at most `α B` (α per block, so `α / m` per step for blocks of m steps, which is the per-unit-time guarantee of an e-detector);
* with the Bonferroni rule, the expected number of false discoveries over B blocks is at most `α B`.

Each block's null is relative to that block's baseline, which is the recent past: a block-b discovery means "this asset's conditional truncated variance in block b left the band around its level in the calibration window before b". A change that happened before the block started has already been absorbed into the baseline. This is a different null from section 6.3 of `docs/math.md` only in the choice of calibration window; the per-block guarantee does not depend on how the baseline was chosen, as long as it is fixed before the block.

What blocks do **not** give: a false discovery rate for the pool of all discoveries over many blocks, or a family-wise statement for the whole horizon. These are the same limits as for the single-asset monitor (`docs/limitations.md`, section 1, "recalibration after an alarm").

**(c) Relation to existing notions for multi-stream change detection.** Chen, Zhang and Poor (2020) define the false discovery rate of parallel sequential change detection as the expected share of falsely declared streams among declared streams, give a sequential procedure that controls it, quantify its delay in asymptotic regimes, and find that FDR-oriented procedures react much faster than family-wise ones. Chen and Li (2023) keep pre-change streams active as long as possible while controlling, at every time, the expected share of post-change streams among the active ones (a false non-discovery type criterion), and show their procedure is optimal among those that do so. Dandapanthula and Ramdas (2025) work with e-detectors and nonparametric nulls, prove the impossibility recalled in section 3, and control error over patience with e-d-BH and related procedures. Wang, Wang and Zhang (2026) propose sequential resetting procedures built from e-values and test supermartingales that locate intervals in which some null hypotheses fail, with explicit FDR bounds. Our construction uses the same two ingredients as Dandapanthula and Ramdas (e-BH and evidence built from Shiryaev-Roberts bet chains) and differs in the error criterion: an ordinary false discovery rate (not divided by expected time) at every stopping time of a block, plus the time-uniform bound of Theorem 1(3), bought by turning the chains into e-processes (section 4) and by restarting on a fixed calendar so that each block's statement needs no assumption on earlier blocks. As we read the abstract of Wang, Wang and Zhang, their hypotheses are indexed by time within a sequence and the output is a set of intervals; here the hypotheses are assets. We make no claim of priority over either.

## 7. Assets with different trading hours

The panel runs on the UTC minute clock. A crypto asset trades in every minute; the Nifty 50 index only in its regular session, 09:15 to 15:30 IST, which is 03:45 to 10:00 UTC (374 one-minute returns per session in the data used here). Outside its session the Nifty return is missing.

* **Own clock.** Each asset's weights are indexed by its own observations, so a window of h = 374 covers one Nifty session while h = 1,440 covers one crypto day. Lemma 2 needs only that the observation pattern is known one step ahead (assumption (M)), which a published trading calendar provides.
* **Frozen evidence when closed.** Outside its session, an asset's e-process does not move. It can still be discovered then: e-BH thresholds depend on how many assets are discovered, so if several crypto assets are discovered at 14:00 UTC, the e-BH level drops and the Nifty evidence collected during its session may clear it. Theorem 1 covers this (the discovery set stays self-consistent), and it is the cross-asset borrowing of strength that e-BH is designed for. A user should read such a discovery as "Nifty's session today already looked different, and the market-wide picture makes it count".
* **Blocks.** One block per UTC day contains exactly one Nifty session, so the per-block guarantee is per crypto day and per Nifty session at the same time. Calibration windows are counted on each asset's own clock: the last 1,440 crypto returns, the last 750 Nifty returns (about two sessions, the primary setting of E7).
* **No common time needed for returns.** Returns are never synchronised or interpolated; each asset only uses its own returns. This avoids the bias that asynchronous trading creates in correlation estimates (the Epps effect discussed in `alternatives/correlation_monitor/README.md`), which matters for statistics that multiply returns of different assets, not for this panel.

The panel null for Nifty conditions on crypto prices before and during the Nifty session. If crypto volatility predicted Nifty volatility, the panel null would be stronger than Nifty's own-past null (Remark 1). E15 does not test this.

## 8. Practical recipe

1. Choose a block (one day, or one session) and α per block.
2. For each asset, set the window h to its expected number of returns per block, and a calibration length on its own clock.
3. At each block start, calibrate (plug-in with a band, or the confidence-sequence baseline with δ for a proved bound α + δ).
4. At every step, update the e-processes and apply sticky e-BH. Report each discovery with its direction (the side with more evidence), its e-value and the level it cleared (`K / (α R)`).
5. Read the guarantee as: in every block, the expected share of false alarms among the alarms raised is at most α, whenever one looks.

## 9. Status and references

* **Proved here:** Lemma 1 (as in Wang and Ramdas 2022), Proposition 1, Lemma 2, Proposition 2, Theorem 1 and Remarks 3 to 6, the block statements of 6(b).
* **Cited, not re-proved:** Ville's inequality and optional stopping for nonnegative supermartingales (see `docs/math.md`, 6.7); the impossibility theorem of Dandapanthula and Ramdas (2025); the results of Wang, Dandapanthula and Ramdas (2025) and of Tavyrikov, Goeman and de Heide (2026) quoted in Remarks 1 and 2.
* **Checked numerically:** `tests/test_panel.py` (e-BH against its definition; the identity with the Shiryaev-Roberts statistic; E[E] ≤ 1 at fixed and stopping times; Proposition 1; the pathwise ordering of the three rules; monotone self-consistent sets; family-wise error and false discovery rate under a common factor; streaming against vectorised code). E15 measures FDR, the time-uniform FDP and power in 56 synthetic designs.
* **Assumed, not proved:** the nulls on real data (as in `docs/math.md` 6.5: under volatility clustering the strict null is often false, which is why the real-data runs use a band); assumption (M) on real data; that the intraday profile used to deseasonalise is known.

**How the references were checked.** The network of this environment blocks arxiv.org, the Semantic Scholar API, Project Euclid and most publisher sites for direct fetching, so every reference below was checked with web searches on 2 November 2026: title, authors, venue and year against the publisher, journal or arXiv listing shown in the search results, and the specific statements we use against the abstract or summaries returned by the search. For Dandapanthula and Ramdas (2025), Wang, Dandapanthula and Ramdas (2025), Tavyrikov, Goeman and de Heide (2026), Chen, Zhang and Poor (2020), Chen and Li (2023) and Wang, Wang and Zhang (2026) we could not read the full text, so we describe only what their abstracts state. References already listed in `docs/math.md` or `docs/literature_review.md` (Ville 1939; Shin, Ramdas and Rinaldo 2024; Pollak 1985; Lorden 1971) are not repeated in the checks.

* Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS B* 57(1), 289 to 300. doi:10.1111/j.2517-6161.1995.tb02031.x. Checked: https://www.dcscience.net/Benjamini-Hochberg-1995-FDR.pdf
* Benjamini, Y. and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. *Annals of Statistics* 29(4), 1165 to 1188. doi:10.1214/aos/1013699998. Checked: https://projecteuclid.org/euclid.aos/1013699998
* Chen, J., Zhang, W. and Poor, H. V. (2020). A false discovery rate oriented approach to parallel sequential change detection problems. *IEEE Transactions on Signal Processing* 68, 1823 onwards (DOI not confirmed). Checked: https://collaborate.princeton.edu/en/publications/a-false-discovery-rate-oriented-approach-to-parallel-sequential-c/ and https://par.nsf.gov/biblio/10194694
* Chen, Y. and Li, X. (2023). Compound sequential change-point detection in parallel data streams. *Statistica Sinica* 33, 453 to 474. doi:10.5705/ss.202020.0508. Checked: https://www3.stat.sinica.edu.tw/statistica/J33N1/J33N120/J33N120.html
* Dandapanthula, S. and Ramdas, A. (2025). Multiple testing in multi-stream sequential change detection. arXiv:2501.04130. Checked: https://arxiv.org/abs/2501.04130 (listing seen through search results; statements used: the impossibility of finite average run length with nontrivial worst-case FDR and family-wise error, the error-over-patience criterion, and the e-d-BH and e-d-Bonferroni procedures).
* Katsevich, E. and Ramdas, A. (2020). Simultaneous high-probability bounds on the false discovery proportion in structured, regression and online settings. *Annals of Statistics* 48(6), 3465 to 3487. doi:10.1214/19-AOS1938. Checked: https://projecteuclid.org/journals/annals-of-statistics/volume-48/issue-6/Simultaneous-high-probability-bounds-on-the-false-discovery-proportion-in/10.1214/19-AOS1938.full (related work on time-uniform FDP bounds by martingale arguments; not used in a proof).
* Lai, T. L. (1998). Information bounds and quick detection of parameter changes in stochastic systems. *IEEE Transactions on Information Theory* 44(7), 2917 to 2929. Checked: https://ieeexplore.ieee.org/abstract/document/737522
* Ramdas, A., Grünwald, P., Vovk, V. and Shafer, G. (2023). Game-theoretic statistics and safe anytime-valid inference. *Statistical Science* 38(4), 576 to 601. doi:10.1214/23-STS894. Checked: https://projecteuclid.org/journals/statistical-science/volume-38/issue-4/Game-Theoretic-Statistics-and-Safe-Anytime-Valid-Inference/10.1214/23-STS894.full (background on test supermartingales and e-processes).
* Tavyrikov, Y., Goeman, J. J. and de Heide, R. (2026). Carefree multiple testing with e-processes. *Electronic Journal of Statistics* 20, 3178 to 3189. doi:10.1214/26-EJS2546. Checked: https://ir.cwi.nl/pub/36622 and https://arxiv.org/abs/2501.19360
* Wang, H., Dandapanthula, S. and Ramdas, A. (2025). Anytime-valid FDR control with the stopped e-BH procedure. *Statistics & Probability Letters* (special issue on e-values and multiple testing; volume and pages not confirmed). arXiv:2502.08539. Checked: https://www.sciencedirect.com/science/article/abs/pii/S0167715225001579 and https://arxiv.org/abs/2502.08539
* Wang, Q., Wang, R. and Zhang, Z. (2026). Sequential resetting procedures and false discovery rate. arXiv:2610.09339 (submitted 7 October 2026). Checked: https://arxiv.org/abs/2610.09339
* Wang, R. and Ramdas, A. (2022). False discovery rate control with e-values. *JRSS B* 84(3), 822 to 852. doi:10.1111/rssb.12489. Checked: https://sas.uwaterloo.ca/~wang/papers/2021Wang-Ramdas-JRSSB.pdf and https://ideas.repec.org/a/bla/jorssb/v84y2022i3p822-852.html
* Xu, Z. and Ramdas, A. (2024). Online multiple testing with e-values. *Proceedings of AISTATS 2024*, PMLR 238. Checked: https://proceedings.mlr.press/v238/xu24a/xu24a.pdf
* Xu, Z., Wang, R. and Ramdas, A. (2021). A unified framework for bandit multiple testing. *NeurIPS 2021*. arXiv:2107.07322. Checked: https://proceedings.neurips.cc/paper_files/paper/2021/file/8c460674cd61bf189e62b4da4bd9d7c1-Paper.pdf
