# Independent review, round 10: E24, the guard against blind segments of the cs baseline

> *Note added 3 November 2026: the "independent reviewer" or "referee" below was a separate AI agent (Claude, by Anthropic) instructed to act as a sceptical referee. It was not a human expert, and the review has not been confirmed by a human.*

Scope: commits `4c294e4`, `09866ec`, `22f1988`, `d83e5f0` and `5df43c8`, reviewed at `95a4ed5` (no file in scope changed after `5df43c8`). Files: `docs/theory_blind_guard.md`, `src/avmonitor/guard.py`, `tests/test_guard.py`, `docs/report_parts/E24.md`, `docs/iterations/E24.md`, `docs/preregistration_e24_future.md`, `data/events/us_macro_events_e24_future.csv`, and the parts of `experiments/e24_blind_guard.py` they rest on. Every number below was recomputed with the scripts in `docs/review/round10_scripts/` (outputs in `docs/review/round10_scripts/outputs/`). Labels: [PROVED] re-derived by hand, [NUMERICAL] exact arithmetic, [SYNTHETIC] simulation, [REAL] recomputed from committed real-data CSVs, [PUBLISHED] official schedules found by web search. No file outside `docs/review/` was edited; fixes are given as replacement text.

## Recommendation

**Accept after revision. No experiment needs to be rerun.**

Lemma 1, Lemma 2, Propositions 1, 2 and 4, Theorems 1, 2 and 3 and the upper bound of Remark 6 are correct as proved. The code implements what they assume: one budget δ_k/2 per side per window, the decision at each window's end, consecutive windows, and a rank fallback that restarts at every period boundary. Every number in the E24 report part that comes from a committed CSV reproduces. The primary forward test is fully specified, and its decision rule was committed before the committed power run.

Two statements must be corrected before the note is cited:

* **M1.** Proposition 3 and the report read more into the construction than it shows. Part (a) of the proposition is false as stated. With three windows the naive guard stays 13 times below δ on this law.
* **M2.** The null of Theorem 2 depends on the guard's tuning through the selection. One example in Remark 3 is false. A selection-free hypothesis fixes both.

Eleven minor points (m1 to m11) need text or small code fixes. Two of them protect the forward test: the future part does not catch a deleted FOMC row, and it does not rerun its own gate (m8). Two tests will fail at the first lawful update of the event file (m9).

## Major findings

### M1. Proposition 3: part (a) is false, and the numbers do not show that the union bound "cannot be improved" or that the naive guard breaks δ

**Where.**
* `docs/theory_blind_guard.md:11` (summary item 3).
* `docs/theory_blind_guard.md:75` (part (a)) and `:79` (its proof).
* `docs/theory_blind_guard.md:77` (last sentence of part (c)).
* `docs/theory_blind_guard.md:81`.
* `docs/report_parts/E24.md:8` and `:59`.
* The same reading appears in `ITERATION_LOG.md:611` ("the naive guard does not"). That file is not mine to edit.

**Problem 1: part (a) is false as stated.** It says the tolerance holds "relative to every window, for every horizon". After a discarded window k, the next n_τ steps are window k + 1's τ part. There the returns are ±s, so μ_t = s² = τ² / (1.4826 c)² = 0.051 τ² for c = 3. (H↓) relative to window k needs μ_t ≥ (1 - b) q τ² = 0.6 τ². So the tolerance fails relative to every discarded window. Part (a) holds relative to the kept window, up to the segment's alarm.

This does not affect parts (b) and (c), which are coverage statements and need no tolerance. It does matter for how the law is read: the law is inside the guarded monitor's null (the kept window) and inside the unguarded monitor's null (window 1 is always kept there). It is not inside the selection-free null (H*) proposed in M2.

**Problem 2: the numbers show less than the text says.**

[SYNTHETIC, `r10_selection.py`, 1,000,000 replications per level for n_cs = 1080, 2,000,000 for n_cs = 150]

| Window | p(δ) for one window | 3 windows, δ each | 3 windows, halving | Committed |
|---|---|---|---|---|
| P1's window (n_cs = 1080, δ = 0.01) | 2.59 × 10⁻⁴ (s.e. 1.6 × 10⁻⁵) | 7.77 × 10⁻⁴ | 1.88 × 10⁻⁴ | 2.525 × 10⁻⁴, 7.57 × 10⁻⁴, 1.87 × 10⁻⁴ |
| Illustration (n_cs = 150, δ = 0.3) | 0.01398 (s.e. 8.3 × 10⁻⁵) | 0.0413 | 0.00644 | 0.0139, 0.0411, 0.0063 |

The committed numbers reproduce. Their reading does not:

* With three windows and δ per window, the naive guard's failure probability on this law is 7.6 × 10⁻⁴. That is 13 times below δ = 0.01, so the naive capped guard keeps P(Gᶜ) ≤ δ here. In the illustration it is 0.041 against δ = 0.3.
* The naive guard first exceeds δ at K* = 40 windows with the committed p, or 39 with mine (1 - (1 - p)^K ≥ δ). In the illustration that happens at K* = 26.
* The budget δ/2 is larger than p(δ) by a factor of 19 (δ = 0.01), 23 (0.005) and 32 (0.0025), and by 11, 19 and 30 for n_cs = 150. The bound Σ_k δ_k of Theorem 1 is therefore attained by this law only up to that Ville gap. Line 77 says it "cannot be improved beyond the gap between p(δ) and δ/2", which is true but says nothing. Line 11 says it "cannot be improved in general", which is not shown.

What the construction does show sharply is threefold:
1. Given that a window before the K-th is kept, its bound fails with probability 1.
2. E19's guard with no cap (discard and recalibrate until a window is not blind) keeps a failing window with probability 1. The windows are i.i.d. and each is kept exactly when its bound fails, so the first kept window always fails.
3. With a cap K and δ per window, the failure probability grows about K-fold.

E24.md:59 ("The selection multiplies the failure probability by the number of windows tried ...; spending keeps it below δ") implies that the naive guard does not keep it below δ. With K = 3 on this law, it does.

**Reproduction.** `nice -n 10 python docs/review/round10_scripts/r10_selection.py`. The output is in `outputs/r10_selection.txt`, with the last two lines for K* and the uncapped guard. The script checks its vectorised sequence against `eb_confseq` (largest difference 1.1 × 10⁻¹⁶).

**Fix (replacement text).**

`docs/theory_blind_guard.md:11`:

> 3. **Choosing which window to keep does break per-window coverage** [PROVED, Proposition 3]. There is a law (adaptive; inside the tolerance relative to the kept window, not relative to the discarded ones) under which a window before the K-th is kept exactly when its upper bound fails: given that such a window is kept, its bound fails with probability 1. With no cap on the number of windows (E19's proposal) the kept window's bound then fails with probability 1. With a cap K and the same budget δ for every window it fails with probability `1 - Π_k (1 - p(δ))`, about K times the single-window failure probability p(δ). For P1's window p(δ) is about 1/20 of δ/2, so on this law the naive guard with K = 3 stays below δ (7.6 × 10⁻⁴ against 0.01) and exceeds δ only from K = 40. Spending δ over the windows makes the bound hold for every law and every K.

`docs/theory_blind_guard.md:75`:

> (a) the tolerance (H↑), (H↓) holds relative to the kept window of every segment, for every horizon up to the segment's alarm and every b ≥ 0 (μ_t = q τ² in every CS part and during monitoring). It fails relative to a discarded window k: the next n_τ steps are window k + 1's τ part, where `μ_t = s² = τ² / (1.4826 c)²`, below `(1 - b) q τ²` (for c = 3 and b = 0.25: 0.051 τ² against 0.6 τ²). So the law is in the null of Theorem 2, but not in the selection-free null (H*) of Remark 3;

`docs/theory_blind_guard.md:79`, proof of (a):

> (a) `μ_t = E[min(r_t², τ²) | F_{t-1}] = q τ²` outside the τ parts; after the kept window every step up to the segment's alarm is a monitoring step, so no τ part follows it within the segment.

`docs/theory_blind_guard.md:77`, last sentence:

> With the same budget δ for every window (no spending) this is about `K p(δ)`, K times the single-window failure probability; it exceeds δ once K is about `δ / p(δ)`. The bound `Σ_k δ_k` of Theorem 1 is attained by this law only up to the factor `p(δ) / (δ / 2)` (1/19 to 1/32 for P1's window).

`docs/theory_blind_guard.md:81`, first two sentences:

> So E19's worry is right in this sense: conditioning on being kept destroys the per-window coverage; without a cap the kept window's bound fails with probability 1, and with a cap only the unconditional union over the windows tried survives. [NUMERICAL, `results/e24_selection.csv`] For P1's window (n_cs = 1080, δ = 0.01, b = 0.25, η = 0, so q = 0.8), one window's upper sequence fails with probability p(δ) = 2.5 × 10⁻⁴ under the construction (400,000 replications per level; the budget δ/2 = 0.005 is conservative by a factor of 20 here); with three windows and δ for each, the kept window's bound fails with probability 7.6 × 10⁻⁴, three times as often but still 13 times below δ (the naive guard exceeds δ on this law from K = 40); with the halving budgets of Theorem 1, 1.9 × 10⁻⁴.

`docs/report_parts/E24.md:8`, the last sentence:

> [NUMERICAL] For P1's window the construction triples the failure probability with three windows and no spending (7.6 × 10⁻⁴ against 2.5 × 10⁻⁴; `e24_selection.csv`), which is still below δ = 0.01; the naive guard exceeds δ on this law only from about 40 windows, and with no cap on the number of windows the kept window's bound fails with probability 1.

`docs/report_parts/E24.md:59`:

> The selection multiplies the failure probability by the number of windows tried. On this law the naive guard with three windows still stays below δ (7.6 × 10⁻⁴ against 0.01; 0.041 against 0.3), because one window's sequence fails far less often than its budget δ/2; it exceeds δ from about 40 windows (26 in the illustration), and with no cap the kept window's bound fails with probability 1. Spending keeps the failure probability below δ for every law and every K.

### M2. The null of Theorem 2 depends on the guard's tuning, and Remark 3's first example is false

**Where.**
* `docs/theory_blind_guard.md:12` (summary item 4).
* `docs/theory_blind_guard.md:83` to `:89` (Theorem 2) and `:98` (Remark 3, second bullet).
* `src/avmonitor/guard.py:31` to `:32`, `:241` and `:267` to `:268`.
* `docs/report_parts/E24.md:8` ("with the tolerance taken relative to the kept window").

**Problem.** Theorem 2 is correct as proved. Its hypothesis is almost sure and is checked path by path, and items (i) to (iv) of the proof hold; see Q2. But which window is kept depends on U_T, so it depends on δ_k, on the predictable bets and on η. The null of Theorem 2 is therefore a set of laws that changes with the guard's tuning. Two guards with different δ_k or η can have different nulls on the same law.

`docs/math.md` 6.3 states that the hypotheses do not involve the bets, and the remark after Theorem 6.4 rejected a bet-weighted tolerance for that reason. Theorem 2 brings that dependence back through the selection. Summary item 4 says "neither null implies the other". It does not mention that a common, selection-free sufficient condition exists and covers both monitors.

**The example at line 98 is false.** It says: "If volatility rises in the middle of window 1's CS part and stays, the tolerance fails relative to window 1 but holds relative to window 2 (the guarded monitor is covered ...)". Theorem 2 covers the guarded monitor only on paths where window 1 is discarded, which needs ρ_1 ≤ 1 + η. A rise inside the CS part raises U_T, but a rise in the middle of the CS part leaves ρ_1 far above 1.25, so window 1 is kept.

[SYNTHETIC, `r10_tolerance_example.py`, P1's window, i.i.d. Gaussian returns, standard deviation multiplied by k from minute s, 4000 paths per cell]. The table gives the share of first windows kept by the recommended guard (η = 0.25):

| Rise at minute s | k = 2 | k = 4 | k = 8 | k = 16 | k = 32 |
|---|---|---|---|---|---|
| 360 (start of the CS part) | 1.000 | 0.133 | 0 | 0 | 0 |
| 540 | 1.000 | 1.000 | 0.756 | 0.079 | 0 |
| 900 (middle of the CS part) | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 1260 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

The median ρ_1 is 1.98 to 3.33 at s = 900. In the remark's own scenario, then, the guarded monitor is outside the null of Theorem 2 on every path, just like the unguarded one. A direct run with the rise at minute 900 (100 paths) gives these first-day alarm shares:

| Rise | Unguarded P1 | Guarded P1 (GE) |
|---|---|---|
| k = 2 | 0.93 | 0.89 |
| k = 4 | 1.00 | 1.00 |

Part 4 of Theorem 6.4 still applies, with P(V_m) ≥ P(window 1 kept).

**Reproduction.** `nice -n 10 python docs/review/round10_scripts/r10_tolerance_example.py`. The output is in `outputs/r10_tolerance_example.txt`.

**Fix (replacement and added text).**

After `docs/theory_blind_guard.md:25` (end of the tolerance definition), add:

> The **tolerance relative to every window** (H*) (horizon h): for every k ≤ K, (H↑) and (H↓) hold relative to window k (T = T_k, the minimum and maximum over window k's CS part, μ computed with τ_k) for `T_k < t ≤ T_k + h`, almost surely, where window k is the k-th window the guard would build from S (defined on every path, whether or not it is tried). (H*) depends on the law, c, b, K and the window lengths, not on δ_k, η or the bets. It implies the tolerance relative to the kept window (κ ≤ K), and its k = 1 part is the hypothesis of the unguarded monitor.

Add to Theorem 2 (after `:87`):

> 4. Parts 1 to 3 hold in particular under (H*), which does not depend on which window is kept.

Replace the first two sentences of `docs/theory_blind_guard.md:98`:

> * The guarded and unguarded monitors test different nulls, and neither implies the other, but (H*) implies both. If volatility rises in the CS part of window 1 and stays, the tolerance fails relative to window 1 and holds relative to window 2. The guarded monitor is covered by Theorem 2 only on paths where window 1 is discarded (ρ_1 ≤ 1 + η): for P1 with η = 0.25, a rise by a factor of 8 or more at the start of the CS part does that on every path, but a rise in the middle of the CS part leaves ρ_1 near 2 to 3, and window 1 is kept on every path [SYNTHETIC, round 10 review]. There both monitors are outside the null, and part 4 of Theorem 6.4 charges the guarded monitor `P(V_m) ≥ P(window 1 kept)`.

Replace `docs/theory_blind_guard.md:12`, second and third sentences:

> The tolerance is an almost-sure hypothesis about conditional means, so the selection does not enter the proofs; but which window is kept depends on the sequence's output (δ_k, the bets, η), so this null depends on the guard's tuning. A sufficient condition that does not is the tolerance relative to each of the K windows the guard may try ((H*), section 1), which implies both the guarded and the unguarded hypotheses. It holds under a constant conditional law from the end of the first window's τ part on (Corollary 6.5 of `docs/math.md`), which covers a step in volatility inside that τ part.

In `src/avmonitor/guard.py:31` to `:32`, `:241` and `:267` to `:268`, after "relative to the kept window", add: "(sufficient, and independent of the selection: the tolerance relative to each of the K windows, docs/theory_blind_guard.md (H*))".

## Minor findings

### m1. Off by one step, and the carry, in the blind condition

**Where.** `docs/theory_blind_guard.md:9` and `docs/report_parts/E24.md:20` say "cannot alarm within 2A monitored steps". `docs/theory_blind_guard.md:29` (Lemma 1) and `:47` (Proposition 2(b)) say "with a restart every m < 2A steps it can never alarm".

**Problem.** At ρ = 1 the data X ≡ τ² give R_t = t - t₀ exactly. The increase side can therefore alarm at step 2A, and the correct statement is "within 2A - 1", as Lemma 1 itself says. With a carry κ the statistic restarts at up to κ, not 0, so "never" needs m + κ < 2A, and Proposition 2(b)'s floor assumes carry 0.

**Reproduction.** `outputs/r10_near_blind.txt`, line "rho = 1: a_(2A - 1) = 71999 < 2A, a_(2A) = 72000 = 2A".

**Fix.**
* In `:9` and `E24.md:20`, replace "within 2A monitored steps" with "within 2A - 1 monitored steps".
* In `:29`, replace "with a restart every m < 2A steps (P1: m = 1440, 2A = 72,000) it can never alarm" with "with a restart every m steps and a carry of at most κ it can never alarm when m + κ < 2A (P1: m = 1440, κ = 0, 2A = 72,000)".
* In `:47`, replace "With a restart every m steps" with "With a restart every m steps and carry 0".

The same wording sits in `docs/report_parts/E19.md:135`, which is outside this scope.

### m2. Two estimates of the same p in one paragraph

**Where.** `docs/theory_blind_guard.md:81` ("the formula gives 0.0139 for one window and 0.041 for three ... against 4.26% from the formula") and `docs/report_parts/E24.md:57` ("formula 0.0426").

**Problem.** The 4.26% uses a second 200,000-replication estimate, p = 0.0144 (`experiments/e24_blind_guard.py:186`). With the paragraph's own p = 0.0139 the formula gives 4.11%. My 2,000,000-replication estimate is p = 0.01398 (s.e. 8.3 × 10⁻⁵), which gives 4.13%. Every version agrees with the direct 4.3% (s.e. 0.32 points), so only the text is inconsistent.

**Fix.**
* In `:81`, replace "against 4.26% from the formula" with "against 4.1% from the formula with p = 0.0139 (4.26% with that run's own estimate, p = 0.0144 from another 200,000 replications)".
* In `E24.md:57`, replace "formula 0.0426" with "formula 0.041; 0.0426 with that run's own p = 0.0144".

### m3. The filtration, what Theorem 1 uses, and a missing Proposition 5

**Where.** `docs/theory_blind_guard.md:19`, `:10`, `:71` and `:13`.

**Problem.**
* Line 19 puts "any independent randomisation the monitor uses" in F_0. Theorem 3 (line 129) needs the tie-breaking draws to be independent of F_{T_K}, so that the p-values are i.i.d. uniform given F_{T_K}. If the draws are in F_0 ⊆ F_{T_K}, the p-values are not uniform given F_{T_K} whenever returns tie, and BTC's τ parts hold 9% to 35% exact zeros.
* Lines 10 and 71 say the proof of Theorem 1 uses that the windows are disjoint. The proof (line 69) uses W_k ∈ F_{T_{k-1}} ⊆ F_{σ_k}, that is T_{k-1} ≤ σ_k, together with one stopping time and one budget per sequence. Disjointness is not used.
* Line 13 cites "Propositions 4 and 5". There is no Proposition 5; the fallback is Theorem 3.

**Fix.**
* In `:19`, replace "(with any independent randomisation the monitor uses added to `F_0`)" with "(with the monitor's tie-breaking draw U_t added to F_t at the step t where it is used; the draws are i.i.d. uniform and independent of the returns, so conditional means of the returns are unchanged)".
* In `:10`, replace the last sentence with "The proof uses only that whether window k is tried is decided by the time its sequence starts (T_{k-1} ≤ σ_k), and that each sequence starts at a stopping time with its own budget δ_k."
* In `:71`, replace "the windows do not overlap, so window k's sequence is a fresh sequence started at its own stopping time; whether window k is tried is known before its sequence starts" with "whether window k is tried is known before its sequence starts (W_k ∈ F_{σ_k}), and window k's sequence is a fresh sequence started at its own stopping time; that the windows do not overlap is not needed".
* In `:13`, replace "[PROVED, Propositions 4 and 5, Theorem 3]" with "[PROVED, Proposition 4, Theorem 3]".

### m4. Remark 6: the sketch's mechanism is not available inside the null, and the bound can be made explicit

**Where.** `docs/theory_blind_guard.md:111` (sketch) and `:14` (summary item 6).

**Problem.** The sketch ends the other segments "quickly with alarms". Inside the null that cannot happen: a good segment alarms in any period with probability at most (m + κ)/A (Corollary 1). Under the tolerance with h = ∞ in every segment, Corollary 6.6 (Theorem 2, part 2) also bounds the number of earlier segments. Combining the two gives an explicit bound. [PROVED]

E[#{segments whose calibration starts before b'}] ≤ 1 + E[N_{b'}] ≤ 1 + (δ + b'/A)/(1 - δ), because at most b' steps before b' are monitored. So for the d-th period of the stream clock (b' = d m):

`P(an alarm in it) ≤ δ (1 + (δ + d m / A) / (1 - δ)) + (m + κ) / A`.

For P1 this is 5.05% on day 1, 6.02% on day 25 and 9.05% on day 100 [NUMERICAL]. The excess over δ + m/A grows by at most δ m / (A (1 - δ)) = 0.04 points per day. Whether some adaptive null attains more than δ + (m + κ)/A on a given day is open.

**Fix.** Replace the sketch at `:111` with:

> The first term grows with the number of earlier segments, but inside the null it grows slowly. Under the tolerance with h = ∞ in every segment, Theorem 2 part 2 gives `E[#{segments whose calibration starts before b'}] ≤ 1 + (δ + b' / A) / (1 - δ)`, so for the d-th period of the stream clock `P(an alarm in it) ≤ δ (1 + (δ + d m / A) / (1 - δ)) + (m + κ) / A` [PROVED]; for P1, 5.05% on day 1, 6.02% on day 25 and 9.05% on day 100. A null cannot end good segments quickly, since each period of a good segment alarms with probability at most `(m + κ) / A`. Whether an adaptive null attains more than `δ + (m + κ) / A` on some day is open; we have not constructed one. So "δ + m / A for every day" is a statement about the periods of each segment, not about calendar days.

In `:14`, replace "and an adaptive null can push the probability above `δ + m / A`" with "which is at most `δ (1 + (δ + d m / A) / (1 - δ)) + (m + κ) / A` for the d-th day under the tolerance in every segment (5.05% on day 1 and 6.02% on day 25 for P1); whether a null attains more than `δ + m / A` is open".

### m5. The per-day share of E22 and E23 is not "the quantity Corollary 1 bounds"

**Where.** `docs/report_parts/E22.md:48`, `docs/report_parts/E23.md:68` and `docs/preregistration_e24_future.md:103`.

**Problem.** Corollary 1 bounds the probability of an alarm in a given period of a given segment. The share of monitored periods with an alarm pools periods that exist only because their segment has not alarmed yet, so it is not that probability. The right bound comes from counting.

**Counting corollary [PROVED].** Let 𝒫_M be the set of periods with at least one of the first M monitored steps, and suppose the tolerance holds with h = ∞ in every segment, relative to its kept window. Then

`E[N_M] ≤ (δ + ((m + κ) / A) E|𝒫_M|) / (1 - δ)`.

*Proof.*
1. Each alarm ends its segment, so N_M is also the number of alarmed periods.
2. For period p of segment j, starting at b_p (the boundary, or T_j for the first period), the event {p ∈ 𝒫_M} is in F_{b_p}. Corollary 1 gives `P(G_j, alarm in p | F_{b_p}) ≤ 1{G_j, p ∈ 𝒫_M} (m + κ) / A`.
3. Summing over all (j, p) and taking expectations bounds the alarms of good segments by `((m + κ) / A) E|𝒫_M|`.
4. Bad segments add at most `δ E[#{j : C_j}] ≤ δ (1 + E[N_M])`, by Theorem 1 at each start, as in Theorem 3 part 3. ∎

So the long-run share of alarmed periods, as a ratio of expectations, is at most `(m + κ) / (A (1 - δ)) = 4.04%` for P1. Comparing E22's and E23's shares (about one half) with 5% is therefore lenient to the monitor, and their conclusions stand.

**Fix.** In `E22.md:48` and `E23.md:68`, replace "the quantity Corollary 1 bounds" with "whose long-run expectation is at most (m + κ)/(A (1 - δ)) = 4.04% under the tolerance in every segment (each period has probability at most δ + m/A = 5% by Corollary 1; the share pools periods selected by earlier alarms, so its bound comes from counting, round 10 review m5)". The pre-registration must not change after its commit, so state the same reading in the future report part when the test runs.

### m6. Docstrings in `src/avmonitor/guard.py`

**Problems and fixes.**
* **`:43`.** "With guard_tries=1 the guarded monitor is the unguarded one, step for step" holds only with `guard_fallback="keep"`. With "rank", a blind single window falls back to the rank monitor (`outputs/r10_guard_audit.txt`, item 4: "rank fallback" against "blind kept"). Replace with "With guard_tries=1 and guard_fallback="keep" the guarded monitor is the unguarded one, step for step."
* **`:33`.** "with ``period=m``, P(an alarm in a given period of a segment) <= delta + m / A" omits the carry. Replace "delta + m / A (marginal)" with "delta + (m + carry) / A (marginal)". The same omission is in `docs/theory_blind_guard.md:132` (Theorem 3 part 2): there the cs branch has `δ + (m + κ)/A`, and the rank branch has m/A because its statistics restart at 0 (audit item 3).
* **`:243` to `:244`.** "when the hypothesis holds in every segment" points back to the first-m tolerance (h = m), but Corollary 6.6 needs h = ∞. Replace with "when the tolerance holds in every segment up to its alarm (h = infinity)".
* **`:250` to `:251` and `:274` to `:275`.** The defaults `guard_eta=0.0, guard_fallback="keep"` are E19's criterion (arm G0), not the recommended guard, and `docs/api.md` shows them as the API. Add to both class docstrings: "The defaults (guard_eta=0, guard_fallback="keep") are E19's criterion (E24's arm G0); the recommended guard is guard_eta=0.25, guard_fallback="rank" (docs/theory_blind_guard.md 6.1)." Changing the defaults would change the public API, so that is the owner's call.
* **`:130` and `:137`.** `self.delta` is overwritten with the current window's budget. After a run with δ = 0.01, `mon.delta` is 0.0025 and `mon.delta_total` is 0.01 (audit item 5). Nothing in the repository reads `mon.delta` after a run, but a user computing δ + m/A from it gets 4.25% instead of 5%. Add to both class docstrings: "``delta`` holds the current window's budget δ_k during and after a run; the total budget is ``delta_total``."

### m7. Theorem 3 part 3: "exchangeability of the whole remainder" is not defined

**Where.** `docs/theory_blind_guard.md:133`.

**Problem.** The remainder of a segment has random length, and exchangeability of a random-length block is not a standard notion. The proof (optional stopping at min(alarm, T_j + R_j) with R_j ≤ M) needs that, given F_{T_K}, the returns r_{T_K+1}, ..., r_{T_K+M} are exchangeable for the M in question. Exchangeability of the infinite sequence after T_K is enough for every M.

**Fix.** Replace "(h = ∞; exchangeability of the whole remainder of a fallback segment)" with "(h = ∞; in a fallback segment, the returns r_{T_K + 1}, r_{T_K + 2}, ... exchangeable given F_{T_K}, or at least the first M of them)".

### m8. The future part: a deleted FOMC row passes, and the gate is not rerun at run time

**Where.** `experiments/e24_blind_guard.py:930` to `:953` (`future_refusals`) and `:956` to `:990` (`part_future`). The claim is in `docs/preregistration_e24_future.md:3` and `docs/report_parts/E24.md:13`.

**Problem.** [`r10_future_refusals.py`, clock patched, edited copies of the event file in a temporary directory]

| Case | Result | Expected? |
|---|---|---|
| (a) Committed file, today | 12 refusals | yes |
| (b) Committed file, after the end date | 11 refusals | yes |
| (c) Every row updated, after the end date | none | yes |
| (c') The same file, one day before the end date | the date refusal | yes |
| (e) The CPI row of 11 February 2027 deleted | caught | yes |
| (g) Data ending at 31 May 23:59 | refused | yes |
| (g) Data ending at 1 June 00:00 | passes | yes |
| (h) An extra row left "scheduled" | caught | yes |
| **(d) The FOMC row of 17 March 2027 deleted instead of updated** | **no refusal** | **no** |
| (f) A row moved past the end date as "rescheduled" | leaves the test | yes, by design |

The update rule forbids deletion, but the code does not enforce it.

Separately, `part_future` never reruns `part_future_gate`. If E12, E22, E23, `running.py` or `periodic.py` change before June 2027, the pipeline can drift from the registered one without any refusal. These modules are shared, and an audit was changing some of them during E24 (`docs/iterations/E24.md:3`). Last, "unverified" is checked through `status` alone.

**Fix.** These checks are additions that only refuse more often and change nothing in the test. Commit them as a disclosed amendment and note it in the future report part.

```python
# experiments/e24_blind_guard.py, next to FUTURE_EVENTS
FUTURE_FOMC_SCHEDULED = ("2026-10-28", "2026-12-09", "2027-01-27", "2027-03-17", "2027-04-28", "2027-06-09",
                         "2027-07-28", "2027-09-15", "2027-10-27", "2027-12-08")  # second meeting days, Fed schedule

# in future_refusals, after the status loop
    fomc = ev[ev.event_type == "fomc"]
    seen = set(fomc.date) | set(fomc.original_date)
    for d0 in FUTURE_FOMC_SCHEDULED:
        if pd.Timestamp(d0) < pd.Timestamp(end) and d0 not in seen:
            why.append(f"no FOMC row for the scheduled meeting of {d0} (update its status; rows are never deleted)")
    rel = ev[inside & (ev.status == "released")]
    for _, r in rel[rel.access_date < rel.date].iterrows():
        why.append(f"{r.date} {r.event_type}: 'released' with an access date before the release; record the verification")
```

In `part_future`, before any download:

```python
    gate = part_future_gate(procs)          # make part_future_gate return its DataFrame g
    if int(gate.mismatches.sum()):
        raise SystemExit("future part refused: the pipeline no longer reproduces E23's committed rows")
```

Also write `git rev-parse HEAD` and `git diff --stat 22f1988 -- experiments/e12_crypto_events.py experiments/e22_guaranteed_events.py experiments/e23_fresh_releases.py src/avmonitor` into an output file of the future part.

### m9. `tests/test_guard.py`: two tests will fail at the first lawful update, and one checks an empty period

**Where and fix.**
* **`:359`.** The test asserts `(ev.access_date == "2026-11-03").all()`. The update rule requires a new access date with each update, so this fails at the first lawful update. Replace with `(ev.access_date >= "2026-11-03").all()`.
* **`:342`.** The test asserts that the live event file still has rows with status "scheduled" before 2027-07-01. It will fail once the file is fully updated, which must happen before the test runs. Run it on a copy instead:
  ```python
  def test_future_part_refuses_early_and_with_unverified_releases(tmp_path, monkeypatch):
      E24 = _e24()
      import pandas as pd
      ev = pd.read_csv(E24.FUTURE_EVENTS, dtype=str, keep_default_na=False)
      ev.loc[ev.date >= "2026-11-01", "status"] = "scheduled"
      p = tmp_path / "events.csv"
      ev.to_csv(p, index=False)
      monkeypatch.setattr(E24, "FUTURE_EVENTS", p)
      ...  # the existing assertions unchanged
  ```
* **`:288` to `:293`.** `b0 = s + (-s) % m` equals s when s is on a period boundary. Here s is always 200, 400 or 600 with m = 100, so the first "period" [s, s) is empty, `per_counts[0]` is always 0, and only one period beyond the first-m check is tested. Replace with `b0 = (s // m + 1) * m  # first boundary strictly after s (restarts happen at t > start)`. With that change my audit (item 6) runs the same check on the rank branch and finds per-period rates of 0.000, 0.003, 0.000 and 0.000 against the bound 0.2.

### m10. The timing of the power rule and of η, and the post hoc labels

**Where.** `docs/report_parts/E24.md:13`, `:82` and `:101`; `docs/theory_blind_guard.md:15` and `:152`; `docs/preregistration_e24_future.md:62`.

**What the record shows.**
* The decision rule is in the docstring of `part_power` in `09866ec` (08:30:35). It is unchanged in `22f1988` (08:46:32) apart from the added sensitivity block. `FUTURE_END` and the primary arm are `None` in `09866ec` and set in `22f1988`. The run at 08:39 rests on the notes but is consistent with the commit order.
* `docs/iterations/E24.md:113` records a scratch power check of S2 alone (2,000 replications; "S2 pooled about 0.88 by July 2027"). It was made in iteration 4, which was "written while the synthetic part ran" (`:99`, 07:21 to 08:24, `:55`). That is before the rule was committed, so the rule's candidate set, which includes S2, was written after a look at S2's power. No data from the test period were involved, so the test is not biased, but "fixed before the power calculation ran" needs this qualification.
* η = 0.25, K = 3 and the monitor list are in `4c294e4` (07:14:57), before the synthetic run (07:21), as `E24.md:101` says. They came after iteration 1's scratch look at P1's 2026 BTC windows with ρ ≤ 1.2 (`docs/iterations/E24.md:25` and `:29`), and the margin 1.25 covers those four windows. The synthetic experiments chose GE among G, GR and GE; they did not choose η.
* The design `frozen` (disclosed at `E24.md:103`) and the null-segment filter `seg1_start` were added after `4c294e4`. They were committed only with the results, in `09866ec`.

**Fix.**
* In `E24.md:13`, replace "By a rule fixed before the power calculation ran" with "By a rule committed before the power part ran (a scratch power check of S2 alone, 2,000 replications, preceded it; `docs/iterations/E24.md`, iteration 4)".
* In `E24.md:82`, replace "Rule (committed before the power part first ran)" with "Rule (committed before the power part first ran, after a scratch power check of S2 alone)".
* In `E24.md:101`, use: "The near-blind margin η = 0.25 and K = 3 were fixed before the synthetic run (commit `4c294e4`), after a 3-path smoke run and after iteration 1's scratch look at P1's 2026 BTC windows with ρ ≤ 1.2; other values were not explored. The design `frozen` and the null-segment filter were committed with the results."
* In `theory:15`, replace "fixed from the synthetic experiments before any real-data rerun" with "the monitor chosen from the synthetic experiments, before any real-data rerun; η = 0.25 was set before those experiments, after a scratch look at P1's real BTC windows".
* In `theory:152`, add: "η = 0.25 itself was set before the synthetic run, after iteration 1's scratch look at real BTC windows with ρ ≤ 1.2."
* The pre-registration stays as committed; state the scratch check in the future report part.

### m11. Three numbers in the report text

**Where and fix.**
* **"About one day per three months"** (`E24.md:76`, `docs/iterations/E24.md:93`). The guard discarded 8 windows in 711 + 297 = 1008 BTC days, one per 126 days, which is 4.1 months [REAL]. Replace with "about one day per four months (8 windows in 1008 days)".
* **"390 to 1010 times larger"** (`E24.md:76`, `docs/iterations/E24.md:89`). These are ratios of deseasonalised scales (`cs_part_scale / tau_part_scale`: 390, 629 and 1012), while "1.48 price ticks" in the same sentence is a raw-return scale. The raw ratios (`cs_part_scale_logret / tau_part_scale_logret`) are 542, 617 and 1121 [REAL]. Replace with "390 to 1010 times larger (deseasonalised returns, which the monitor sees; 540 to 1120 in raw log returns)".
* **"About twice as fast (25 to 30 minutes against 50 to 64)"** (`E24.md:47`). This holds for step6, step8, frozen and drift8 (GE 48 to 64 minutes). In step2 the comparison is 25 against 144, and in step4 27 against 49. Under Student t3 returns the rank monitor is slower (127 against 115). Replace with "The rank monitor detects faster in eight of the nine designs (median 25 to 30 minutes against 48 to 64 for GE in step6, step8, frozen and drift8; 25 against 144 in step2), and slower under Student t3 returns (127 against 115)".

## Answers to Q1 to Q5

**Q1, Theorem 1's conditioning: correct.**
* The statement is conditional on F_S, as Corollary 6.6 needs (P(G_jᶜ | F_{S_j}) ≤ δ).
* The proof applies Lemma 2 at σ_k with A = W_k ∩ A'. Here A' ∈ F_S ⊆ F_{σ_k} because S ≤ σ_k, and W_k ∈ F_{T_{k-1}} ⊆ F_{σ_k}.
* Lemma 2 itself is correct. 1_A K_{σ+j} is a nonnegative supermartingale started at 1_A, and Ville's inequality gives δ_u P(A ∩ {σ < ∞}).
* G_kᶜ is contained in "the sequence fails at some t > σ_k", because T_k > σ_k.
* The union over k adds the budgets of the windows that are actually tried, so the factor P(window k is tried | F_S) is right.
* The theorem does not give coverage given that a window is kept, and the note says so (line 71). The F_0 point in m3 does not affect Theorem 1.

**Q2, Theorem 2's hypothesis: the proof is correct; the hypothesis needs M2.**
* (i) {T = t} = ∪_k {κ = k, T_k = t} ∈ F_t, with {κ = k} = W_k ∩ B_kᶜ ∈ F_{T_k}.
* (ii) On {κ = k}, the quantities are window k's.
* (iii) G ∩ {T = t} = ∪_k {κ = k} ∩ G_k ∩ {T_k = t} ∈ F_t.
* (iv) Step 2 of Theorem 6.4 applies window by window.
* The issue is not validity but what the null means. It depends on δ_k, η and the bets through the selection, and the motivating example of Remark 3 is not covered (M2). (H*) restores a selection-free null that implies both monitors' hypotheses.

**Q3, Proposition 3: (b) and (c) are correct; (a) is false as stated; the numbers reproduce, but their reading is overstated (M1).**
* (b): in normalised units ν_i = q for every CS observation, so near-blind ⟺ min(U_T, 1) ≥ q and failure ⟺ U_T ≤ q. A non-near-blind window therefore fails.
* (c): the CS data of the windows are independent Bernoulli(q) sequences, and the first window with U < q is kept and fails.
* The numbers: 2.59 × 10⁻⁴ against 2.525 × 10⁻⁴ (s.e. 1.6 × 10⁻⁵), and 7.77 × 10⁻⁴ against 7.57 × 10⁻⁴.
* 7.6 × 10⁻⁴ against 2.5 × 10⁻⁴ is a tripling of the coverage failure, but it is 13 times below δ. The naive guard exceeds δ on this law only with about 40 windows, or with no cap.

**Q4, Theorem 3's optional stopping: correct.**
* R_j is known when segment j starts monitoring, so min(alarm, T_j + R_j) is a bounded stopping time: in F on the cs branch, and in H given F_{T_j} on the rank branch.
* 1_E (N_t - (t - T_j)) is a supermartingale with N_{T_j} = 0 and E in F_{T_j}. Restarts only lower N; the rank statistics restart at 0 (audit item 3).
* Optional stopping gives A P(E, alarm among the first M monitored steps) ≤ E[1_E ℓ_j], and Σ_j ℓ_j ≤ M.
* The bad branch costs δ P(C_j), and Σ_j P(C_j) ≤ 1 + E[N_M]. Summing per-segment expectations avoids the mixed-filtration supermartingale correctly.
* Two wording fixes are needed: m7 (what is exchangeable) and m3 (the draws enter at their time of use).
* [SYNTHETIC] Rank branch, 395 of 400 paths fell back: per-period alarm rates 0.000, 0.003, 0.000 and 0.000 against m/A = 0.2.

**Q5, Remark 6 and the per-day share: the upper bound is correct, the sketch is not, and the E22 and E23 label is imprecise.**
* The upper bound is correct. With calibration at least m long, at most one segment monitors in a period, and G_Jᶜ ⊆ ∪_j G_jᶜ over the segments started before b'.
* The sketch's mechanism (ending good segments quickly with alarms) is not available inside the null. The explicit bound in m4 limits the excess over δ + m/A to at most 0.04 points per day for P1.
* "Per-day share" in E22 and E23 is not the quantity Corollary 1 bounds. The counting corollary in m5 bounds its long-run expectation by 4.04% for P1, so the conclusions of E22 and E23 (shares near one half) stand.

## Checked and correct

**[PROVED], re-derived.**
* Lemma 1, with the correct "2A - 1" in its own statement.
* Proposition 1: B ∈ F_T, and T_B is a stopping time.
* Proposition 2 (a) to (c). The recursion is monotone in f_t and R_{t-1}, and Jensen's inequality gives (c).
* Lemma 2 and Theorem 1 (Q1).
* Theorem 2, items (i) to (iv) (Q2).
* Proposition 3 (b) and (c) (Q3).
* Proposition 4. #{k : C_k} ≤ 1 + N_M + ⌊M/(D m)⌋, and the per-period limit is (δ/D + m/A)/(1 - δ). With D = 7, δ/D = 0.0014.
* Theorem 3, parts 1 to 3 (Q4).
* The upper bound of Remark 6 (Q5).
* Remark 1: √(log 400 / log 200) = 1.063.
* Remark 2.

**[NUMERICAL]** (`r10_near_blind.py`, a closed form for a_d(ρ), checked against the package).
* All twelve d_min values for P1 match `results/e24_near_blind.csv` and the note: 72000, 9854, 2703, 1503, 827, 371, 202, 91, 51, 30, 18 and 12.
* ρ* = 1.0105072094 (CSV `eta_star` 0.0105072). a_1440(1.0105) = 71,701 < 72,000 and a_1440(1.0106) = 75,968, so "ρ ≤ 1.0105" is correct.
* The 390-step setting: d_min = 1217 at ρ = 1.01 and 80 at 1.25, as in the CSV.

**[SYNTHETIC], code audit** (`r10_guard_audit.py`, `outputs/r10_guard_audit.txt`).
* Every window's sequences use δ_k/2 per side, also with c_down = 4.
* The windows are consecutive and disjoint, and the decisions fall at the window ends (199, 399 and 599).
* The rank fallback's statistics are 0 at all 10 boundaries and carried unchanged elsewhere.
* The halving schedule sums to δ (audit item 1), and a new segment starts again with window 1 and δ_1 (code reading, `guard.py:133` to `:139`).

**[REAL], report numbers** (`r10_report_numbers.py`, `outputs/r10_report_numbers.txt`).
* The contrast table for P1, G1, G0 and S2 in both periods: counts, control rates and p-values, watched releases (26, 28, 28, 38 and 11, 12, 12, 16) and blind at release (P1 1 and 3, the others 0).
* ETH and SOL G1 5 of 16, p = 0.023 and 0.032.
* The six BTC releases where G1 and P1 differ, including ρ = 0.948 (29 January 2025) and 0.980 (16 September 2026).
* The G1 period alarm share, 0.456 to 0.525; alarms per monitored day, 0.839 to 1.107; first-day shares, 0.61 to 0.69.
* 3 and 5 discarded BTC windows (4320 and 7200 minutes), all starting on a weekend. ETH and SOL smallest ρ 1.59 and 1.76. Three tick windows; the other five with 9% to 27% zeros and ratios 4.9 to 26.
* E19's blind diagnostic: 33 of 33 unguarded rows match `e19_blind.csv` (part p2, cs) in 8 columns.
* The future gate: 192 rows, 14 columns, 0 mismatches.
* The synthetic table matches `e24_synthetic.csv`. The selection table matches `e24_selection.csv`.
* The power table and the candidates match. S2 pooled reaches 0.80 first at June 2027 (0.798 at May, 0.820 at June in the grid). In the sensitivity grid, the 2026 rates reach 0.80 at 1 October 2027 (0.844), and half the excess reaches 0.80 in August 2028 (0.801).

**The pre-registration.**
* The primary test is fully specified: the S2 constructor with every argument; BTC; FOMC and CPI pooled; inclusion and control rules; the outcome `up_post` within [0, 6 h); a one-sided exact stratified permutation test with its p_i; pass at p ≤ 0.05; one test, so no multiplicity; and the end date 1 June 2027 00:00 UTC, exclusive.
* The rule precedes its result in the commit record (m10 for the scratch check).
* The refusals work except for case (d) (m8).
* The event file's UTC times follow the US daylight saving rules (test at `tests/test_guard.py:356`).

**[PUBLISHED], event dates checked by web search on 3 November 2026.**
* The Federal Reserve's tentative 2027 FOMC schedule (press release of 5 September 2025): 26 to 27 January, 27 to 28 April, 8 to 9 June, 27 to 28 July, 14 to 15 September, 26 to 27 October and 7 to 8 December. These agree with the file. The 16 to 17 March meeting was not shown in the snippets I saw.
* BLS CPI dates 10 December 2026, 13 January 2027, 11 February 2027 and 10 March 2027 agree with the file.
* **Not verified by me:** the CPI dates of 10 November 2026, 13 April 2027 and 12 May 2027, the FOMC statement of 9 December 2026, and the report (28 October 2026) that one 2027 meeting will have no statement. The file marks every row after 28 October 2026 as `scheduled` with its search basis, and the update rule requires verification before the test runs, so this is acceptable for planning.
* The official hosts were not fetched (one fetch of a third-party calendar was refused by the environment, and I did not retry it or route around the block).

## Scripts, commands and sources

All run from the repository root with `nice -n 10`, at most 2 workers. Outputs are in `docs/review/round10_scripts/outputs/`.

| Script | What it checks | Time |
|---|---|---|
| `r10_near_blind.py` | d_min, ρ*, the off by one step | seconds |
| `r10_selection.py` | Proposition 3: p(δ), K-window failure, K*, the Ville gap | 6 minutes |
| `r10_guard_audit.py` | budgets, windows, rank restarts, guard_tries=1, `mon.delta`, the rank-branch Monte Carlo | minutes |
| `r10_tolerance_example.py` | M2: the share of first windows kept, first-day alarms | 2 minutes |
| `r10_future_refusals.py` | m8: cases (a) to (h), no market data | seconds |
| `r10_report_numbers.py` | m11 and the report's numbers from the CSVs | seconds |

`python -m pytest -q`: all tests pass at the reviewed commit.

Sources:
* https://www.federalreserve.gov/newsevents/pressreleases/monetary20250905a.htm (seen through search results)
* https://www.mnimarkets.com/articles/mni-federal-reserve-sets-2027-fomc-meeting-schedule-1757093400287
* https://www.bls.gov/schedule/news_release/cpi.htm
* https://www.bls.gov/schedule/2027/home.htm
