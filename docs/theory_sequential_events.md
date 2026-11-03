# An anytime-valid release test (E25)

This note proves the guarantee behind E25 (`docs/preregistration_e25_sequential.md`): a test of "release windows look like their matched control windows" that may be checked after every release and run without an end date. It uses only the conditional structure inside each stratum, no independence between strata, and no model for returns.

## 1. Setting

Releases (FOMC statements, CPI releases) are taken in time order, i = 1, 2, .... Each release i has a stratum: its own window and n_i control windows at the same weekday and Eastern clock time, one to four weeks before and after it. The control windows are chosen from the calendar only (the release dates and the data's coverage), never from outcomes, no window belongs to two strata, and no control window overlaps in time a window of an earlier stratum (pre-registration, section 4 and Amendment 1). A monitor (E14's awake rank monitor, S2) runs over the whole stream, and each window gets a binary outcome: 1 if the monitor raised an up alarm within six hours of the window's start.

For stratum i write Z_i for the release window's outcome, C_{i1}, ..., C_{in_i} for the controls' outcomes, S_i = Z_i + Σ_l C_{il} and p_i = S_i / (n_i + 1). Let G_j be the information in the outcomes of strata 1 to j (G_0 trivial). The numbers n_i are fixed in advance given the calendar and the data's coverage, which we treat as given.

**Null hypothesis H0.** For every i, conditionally on G_{i-1}, the vector (Z_i, C_{i1}, ..., C_{in_i}) is exchangeable.

H0 says that, once the earlier strata are known, nothing distinguishes the release window from its controls as far as the outcome goes. It does not require the strata to be independent of each other, the outcomes within a stratum to be independent, or the windows' outcome probability to stay constant over time; it does require that whatever moves the outcome over the four weeks either side of a release treats the release window like the others.

## 2. The e-value of a stratum

For ψ > 0 define

    e_i(ψ) = ψ^{Z_i} / (1 - p_i + ψ p_i).

**Lemma 1.** Under H0, E[e_i(ψ) | G_{i-1}, S_i] = 1 for every ψ > 0.

*Proof.* Conditionally on G_{i-1}, the n_i + 1 outcomes are exchangeable; conditionally also on their sum S_i, each of the C(n_i + 1, S_i) arrangements with S_i ones is equally likely, so the release window is one of the windows with outcome 1 with probability S_i / (n_i + 1) = p_i. Given (G_{i-1}, S_i), p_i is fixed and Z_i is Bernoulli(p_i), so

    E[e_i(ψ) | G_{i-1}, S_i] = p_i ψ / (1 - p_i + ψ p_i) + (1 - p_i) / (1 - p_i + ψ p_i) = 1.  ∎

e_i(ψ) is nonnegative, equals 1 when p_i is 0 or 1 (a stratum with no window, or every window, with the outcome) or when n_i = 0, and is a likelihood ratio: for the conditional law of Z_i given S_i it compares "the release window has odds ψ times those of a control" (the noncentral hypergeometric law with one release slot, whose probability of Z_i = 1 is ψ p_i / (1 - p_i + ψ p_i)) with "the release window is a random one of the n_i + 1".

## 3. The e-process and its guarantee

For a finite grid Ψ of odds ratios with weights w_ψ ≥ 0 summing to 1, let

    M_j = Σ_{ψ in Ψ} w_ψ Π_{i ≤ j} e_i(ψ),   M_0 = 1.

**Theorem 1.** Under H0, (M_j) is a nonnegative martingale with respect to (G_j) with M_0 = 1, and for every α in (0, 1),

    P(M_j ≥ 1/α for some j ≥ 1) ≤ α.

Consequently, for any stopping rule τ, possibly unbounded and possibly using information beyond G (for example, a look at the data whenever someone wants), P(M_τ ≥ 1/α, τ < ∞) ≤ α.

*Proof.* e_i(ψ) is a function of (Z_i, S_i), hence G_i-measurable. By Lemma 1 and the tower property, E[e_i(ψ) | G_{i-1}] = E[E[e_i(ψ) | G_{i-1}, S_i] | G_{i-1}] = 1, so for each ψ the product Π_{i ≤ j} e_i(ψ) is a nonnegative martingale with initial value 1, and so is the weighted average M_j. Ville's inequality for nonnegative supermartingales gives the bound on the supremum. The event {M_τ ≥ 1/α, τ < ∞} is contained in {sup_j M_j ≥ 1/α} whatever τ is, which gives the second statement. ∎

With Ψ = {2, 4, 8, 16}, equal weights and α = 0.05, this is the rule of the pre-registration: reject H0 the first time M reaches 20.

**Remarks.**

1. *Why the strata must be disjoint.* If a window belonged to strata j < i, its outcome would be in G_{i-1}, and H0 would ask for the release window of stratum i to be exchangeable with a window whose outcome is already fixed given G_{i-1}, which fails as soon as the release window's outcome is random. The disjoint assignment makes H0 a coherent statement.
2. *No independence between strata.* The fixed-sample stratified permutation test of E12 combines strata through a Poisson-binomial law, which treats them as independent; strata that share controls can push its level above nominal (5.8% for seven FOMC statements sharing controls, E23). The martingale argument needs only the conditional statement within each stratum.
3. *Optional continuation.* The bound holds over an unbounded number of releases, so the test can continue past any date, and an interim look carries no penalty. This is what a fixed-sample test such as F24a cannot offer: extending its sample after a failure would spend more than its level.
4. *What the evidence grows like.* If the release window has outcome probability q_i and its controls p0_i, independently, M grows roughly like exp(Σ_i E[log e_i(ψ*)]) for the best ψ* in the grid, less a constant log(1/w_ψ*) = log 4; a stratum adds at most log(ψ*(n_i + 1)/(n_i + ψ*)), less than log ψ*, when the release has the outcome and no control does. With few controls p_i is a noisy estimate of the control rate, which costs growth compared with an oracle that knows p0_i; that is the price of a null that needs no model for the control rate.
5. *What H0 does not cover.* H0 is a statement about the monitor's outcomes in windows, not about returns, and no property of S2 implies it. A volatility trend over the eight weeks around a release, or an alarm rate that differs between weeks for reasons unrelated to releases, can make it false; controls balanced around the release cancel a linear trend to first order, but after the exclusions most strata are not balanced (pre-registration, Amendment 1). A rejection says that release windows differ from their controls under this monitor, not why.

## 4. Relation to the monitor's own guarantee

The monitor's e-detector applies Ville's inequality to minutes; E25 applies it to releases. The two guarantees are separate: E25's H0 can hold while the monitor's own null (exchangeable returns for S2) fails on most real days, as E14 and E19 found, and vice versa. A crossing of E25 is evidence about releases, not about the monitor's false alarm rate.
