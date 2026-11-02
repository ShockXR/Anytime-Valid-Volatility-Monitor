---
title: "Anytime-valid detection of volatility changes in high-frequency returns with robust truncation"
author: "Shreyan Thanki"
date: "Draft, November 2026"
---

## Abstract

We study the online detection of changes in the volatility of one-minute returns when the returns contain jumps and fat tails and the monitor is inspected after every observation. Following Bapat (2026), each return is clipped at $c$ times a robust scale, the median absolute deviation (MAD), which we estimate from past data only. We monitor the clipped squared return with Shiryaev-Roberts betting e-detectors (Shin, Ramdas and Rinaldo 2024), one for increases and one for decreases. With a known baseline, the average run length to a false alarm is at least $A$ and the probability of any false alarm within $m$ steps is at most $m/A$, for every return distribution and every form of dependence that keeps the conditional mean of the clipped squared return on the null side. A one-sided form of the impossibility result of Bahadur and Savage (1956) shows that no nontrivial e-process exists for the null "the variance has not fallen"; a bounded statistic such as the clipped squared return is needed to detect decreases, and the clipping level then sets a speed limit for them. With a baseline estimated by a confidence sequence the guarantee becomes $\delta + m/A$. We give a non-asymptotic bound on the detection delay, a decomposition of the delay cost of robustness into a truncation loss, a nonparametric loss and a grid loss, and a closed-form cap on the signal that remains after dividing returns by a GARCH forecast. In simulations the bound holds under jumps, Student $t$ returns and price ticks (1.1% to 1.4% false alarms per day against a bound of 5%), while Gaussian likelihood-ratio and CUSUM detectors fail. On five cryptocurrencies, the Nifty 50 index and around FOMC and CPI releases, alarms carry information beyond volatility forecasts, but a ratio rule tuned with hindsight is about as precise: on real data the contribution is a false alarm guarantee with settings fixed in advance, not extra precision.

**Keywords.** Sequential change detection; e-processes; testing by betting; high-frequency volatility; truncation; jumps; anytime validity.

**Status of this draft.** Every number below comes from an experiment of the accompanying repository, cited by its label (for example "(E9)"). Each result is marked **[SYNTHETIC]** (our simulation), **[REAL]** (our run on market data) or **[PUBLISHED]** (printed in Bapat 2026). Lines marked "[PENDING ...]" are placeholders for experiments that were still running when this draft was written.

## 1 Introduction

Volatility in high-frequency markets moves between regimes: calm, stressed, calm again. Risk limits, position sizes and model refits should react when it moves. This paper asks for an online monitor that reads a stream of one-minute returns and raises an alarm when volatility shifts up or down, with a proven limit on false alarms that holds however often the monitor is checked. Two obstacles make this hard.

**Jumps and fat tails.** A few extreme returns dominate sums of squared returns, and even more so sums of fourth powers. A detector built on a Gaussian likelihood treats a single jump as strong evidence of higher volatility. Jump-robust high-frequency estimators address this for a fixed sample: bipower variation (Barndorff-Nielsen and Shephard 2004), nearest-neighbour truncation (Andersen, Dobrev and Schaumburg 2012), threshold estimators (Mancini 2009), and the kernel quarticity estimator of Bapat (2026), which clips every return at $c$ times its MAD before smoothing, so that no single return has unbounded influence.

**Repeated looks.** A fixed-sample test applied at every minute will eventually reject by chance. In our simulations a 60-step $z$-test repeated at every step raised a false alarm within one 390-step day on 71% of Gaussian paths (E4) [SYNTHETIC]. Classical sequential change detection (Page 1954; Shiryaev 1963; Roberts 1966; Lorden 1971; Pollak 1985) controls the average run length to a false alarm, but it needs the pre-change law, usually in parametric form. E-detectors (Shin, Ramdas and Rinaldo 2024) give the same control for composite, nonparametric nulls, using e-processes built by betting (Shafer 2021; Ramdas, Grünwald, Vovk and Shafer 2023; Waudby-Smith and Ramdas 2024).

The robust estimators solve the first problem for a fixed sample and do not address the second. Anytime-valid methods solve the second, but they need a statistic whose conditional mean has a meaning under the null and, for the direction "volatility went down", a bounded one. The clipping step of Bapat (2026) supplies exactly such a statistic. This paper puts the two together.

**What is known and what is ours.** The building blocks are known: Ville's (1939) inequality, test supermartingales and testing by betting (Shafer 2021; Waudby-Smith and Ramdas 2024), e-detectors and their run-length guarantee (Shin, Ramdas and Rinaldo 2024), empirical Bernstein confidence sequences (Howard, Ramdas, McAuliffe and Sekhon 2021; Waudby-Smith and Ramdas 2024), the impossibility result of Bahadur and Savage (1956), and Lorden's (1971) and Lai's (1998) delay lower bounds. Our contributions are the following.

* **A robust statistic and null for online volatility monitoring.** We monitor $X_t = \min(r_t^2, \tau^2)$, the clipped squared return, with a truncation level $\tau = c\,\hat s$ computed from past data only ($\hat s$ is the scaled MAD). The nulls restrict only the conditional mean of $X_t$. Under a fixed innovation shape this mean is a nondecreasing function of the conditional volatility (Lemma 1), so the nulls have a volatility meaning without a model for the return distribution.
* **Validity, with exact conditions.** For Shiryaev-Roberts betting e-detectors on $X_t$ we state and prove the run-length bound $\mathbb{E}[T_A] \ge A$ (Theorem A, an instance of Shin, Ramdas and Rinaldo 2024), the finite-horizon bound $\mathbb{P}(T_A \le m) \le m/A$ (Theorem B), a two-sided rule, and a bound on the expected number of false alarms under recalibration after every alarm (Corollary C).
* **Why a bounded statistic.** A one-sided form of the Bahadur and Savage (1956) impossibility (Proposition 2, Corollary 3): every e-process for the sequential null "the conditional variance has not fallen" is bounded by 1, so no e-detector can detect a fall in the variance itself. A bounded functional is needed; the truncated variance is one with a direct volatility meaning, and truncation then imposes an explicit speed limit on decrease detection (Corollary 7). The argument is elementary and we do not claim it is new; its role in this application is the point.
* **An estimated baseline.** We measure how much a plug-in baseline inflates the false alarm rate, and we prove $\mathbb{P}(\text{false alarm within } m) \le \delta + m/A$ for a baseline estimated by a confidence sequence (Theorem D), with a version for recalibration after every alarm. The confidence sequence is the predictable plug-in empirical Bernstein construction of Waudby-Smith and Ramdas (2024) with one-sided ranges; its use as the baseline of a detector, with this bound, is what we add.
* **Detection delay.** A non-asymptotic bound on the worst-case delay (Theorem 8), first-order optimality of the oracle bet among all rules that control the run length over the nonparametric null (Theorem 9), and a decomposition of the delay cost of robustness into a truncation loss, a nonparametric loss and a grid loss (Corollary 6). Under GARCH clustering, monitoring returns divided by a forecast restores validity but caps the visible signal of a variance break at $(1-\beta)/\alpha$ (Proposition 10, closed form).
* **Evidence.** Simulations of validity and delay against Gaussian likelihood-ratio, CUSUM and repeated-test benchmarks, and of the failure modes (clustering, estimated baselines). Real-data studies on the five cryptocurrencies of Bapat (2026), the Nifty 50 index with India VIX, and FOMC and CPI releases, with random-time, forecast and rate-matched controls.

We also report, briefly and as questions, what an independent reproduction of the published tables of Bapat (2026) found about the estimator (Section 7). We do not claim to have improved the estimator.

**Related work.** Retrospective tests for volatility breaks include the cumulative sums of squares of Inclán and Tiao (1994) and the high-frequency, nonparametric and minimax-optimal test of Bibinger, Jirak and Vetter (2017). Online monitoring for a switch to explosive GARCH volatility on daily data is studied by Horváth, Trapani and Wang (2024). Andersen, Todorov and Zhou (2025) monitor high-frequency prices in real time for local departures from the semimartingale model, a different target. Our literature search (November 2026, not exhaustive) did not find a method that combines jump-robust high-frequency volatility measures with anytime-valid sequential change detection; we state this with that caveat.

**Organisation.** Section 2 sets up the statistic. Section 3 gives the detectors and their guarantees, the impossibility result and the estimated baseline. Section 4 studies the detection delay. Sections 5 and 6 report simulation and real-data evidence. Section 7 summarises the reproduction of the estimator. Section 8 discusses limitations and open problems. Longer proofs and settings are in the appendix.

## 2 Setting and the truncated statistic

**Returns and information.** Let $r_t$ be the one-minute log return ending at minute $t$. Monitoring starts after a calibration window $r_{-n+1}, \dots, r_0$; we write $\mathcal{F}_t = \sigma(r_{-n+1}, \dots, r_t)$ for the information up to time $t$. A quantity is *predictable* if its value used at time $t$ is $\mathcal{F}_{t-1}$-measurable. Everything the monitor uses (truncation level, baseline, bets, intraday profile) must be predictable.

**The truncation of Bapat (2026).** Bapat (2026, eq. 11) clips returns at $\tau = c\,\mathrm{MAD}(r)$ with $c = 3$, $\tilde r_i = \max(-\tau, \min(r_i, \tau))$, where the MAD is computed over the whole sample, and then forms the kernel-weighted local quarticity $\hat Q_t = \tfrac13 \sum_{i=t-h}^{t} K((i-t)/h)\, \tilde r_i^4$ (eq. 12) and its average over $t$ (eq. 16). We keep the clipping and change three things.

1. *Scaled MAD.* We use $\hat s = 1.4826\, \mathrm{median}_i \lvert r_i \rvert$, the MAD about zero scaled to estimate a Gaussian standard deviation (the median of one-minute returns is 0 in all our data sets, E6, E7 [REAL]). With $\tau = c\,\hat s$ and Gaussian returns the cut-off sits at $c$ standard deviations. With the raw MAD it sits at $0.6745\,c$ standard deviations.
2. *Past data only.* $\hat s$ and $\tau$ are computed from the calibration window, so $\tau$ is $\mathcal{F}_0$-measurable and predictable for every monitored $t$. A rolling or exponentially weighted MAD keeps the jump protection of the full-sample MAD in simulation (E3) [SYNTHETIC]: the 95th percentile of the absolute log error of a local variance estimate falls from about 1.4 without truncation to about 0.52 with an online MAD.
3. *The squared, not the fourth-power, statistic.* We monitor
$$
X_t = \min(r_t^2, \tau^2) = \tilde r_t^{\,2} \in [0, \tau^2],
$$
the summand of truncated realised variance, which is the second-moment analogue of eq. (12) with bandwidth $h = 0$. The same detectors applied to $\min(r_t^4, \tau^4)$ are valid as well; for pure variance changes they were 1% to 28% faster for increases and 1.8 to 2.6 times slower for decreases (ALT1) [SYNTHETIC]. The local truncated quarticity enters the variance monitor only through the variance of $X_t$, which governs bet sizing (Section 4).

**Truncated moments.** For $Z \sim N(0,1)$ and a cut-off $a > 0$, write $\varphi, \Phi$ for its density and distribution function and
$$
\kappa_2(a) = \mathbb{E}[\min(Z^2, a^2)] = 2\Phi(a) - 1 - 2a\varphi(a) + 2a^2(1 - \Phi(a)),
$$
$$
\kappa_4(a) = \tfrac13 \mathbb{E}[\min(Z^4, a^4)] = \tfrac13\left[3\big(2\Phi(a) - 1 - 2a\varphi(a)\big) - 2a^3\varphi(a) + 2a^4(1 - \Phi(a))\right].
$$
For $c = 3$: with the scaled MAD ($a = 3$), 0.27% of Gaussian returns are clipped, $\kappa_2 = 0.995$ and $\kappa_4 = 0.964$; with the raw MAD ($a = 2.023$), 4.3% are clipped, $\kappa_2 = 0.925$ and $\kappa_4 = 0.705$ (E2) [SYNTHETIC agreement with the closed form]. For Gaussian returns with standard deviation $\sigma_0$ and the scaled MAD, the null mean of $X_t$ is $m_0 = \kappa_2(c)\sigma_0^2$ and the ratio $B = \tau^2/m_0 = c^2/\kappa_2(c)$ equals 9.045 for $c = 3$. The constant $B$, the largest possible value of $X_t/m_0$, governs the decrease detector.

**Intraday seasonality and missing returns.** Intraday volatility follows a strong pattern (Andersen and Bollerslev 1997). Before monitoring, each return is divided by a predictable profile: for crypto, the median of $\lvert r \rvert$ at the same minute of the day over the previous 20 days, smoothed across neighbouring minutes; for the Nifty 50, a session-aware version that does not wrap the close into the open (E7). Medians limit the effect of jumps on the profile (Boudt, Croux and Laurent 2011). The profile is treated as known; errors in it act as a small violation of the null that no experiment isolates (Section 8). A missing return is skipped. This is a zero bet chosen after seeing that $r_t$ is missing, and it keeps validity when, given the past, whether a return is missing says nothing about its size. Missing returns must not be replaced by zeros, which look like evidence for lower volatility.

## 3 Betting e-detectors on the truncated statistic

### 3.1 Nulls and their volatility meaning

Let $\mu_t = \mathbb{E}[X_t \mid \mathcal{F}_{t-1}]$ be the conditional truncated variance, and let $m_0 > 0$ be a baseline. For a band $b \in [0, 1)$ fixed in advance, the nulls are
$$
H_0^{\uparrow}:\ \mu_t \le (1+b)\,m_0 \ \text{ for all } t \ge 1, \qquad H_0^{\downarrow}:\ \mu_t \ge (1-b)\,m_0 \ \text{ for all } t \ge 1.
$$
They restrict one conditional moment and allow any return distribution, any jump process and any dependence that keeps that moment on the null side. "Any dependence" is a statement about the shape of the joint law, not a licence for any data: the conditional mean must stay on the null side at every step. Persistent volatility clustering and price-level microstructure noise both move it, even when the unconditional variance is constant (Sections 5.4 and 8).

**Lemma 1 (volatility meaning).** Suppose $r_t = \sigma_t \varepsilon_t$, where $\sigma_t > 0$ is $\mathcal{F}_{t-1}$-measurable and $\varepsilon_t$ is independent of $\mathcal{F}_{t-1}$ with a law $G$ that does not depend on $t$. Let $\tau$ be $\mathcal{F}_0$-measurable and $g_\tau(s) = \mathbb{E}_G[\min(s^2\varepsilon^2, \tau^2)]$. Then $\mu_t = g_\tau(\sigma_t)$ and $g_\tau$ is nondecreasing on $(0, \infty)$. Hence, with $m_0 = g_\tau(\sigma_0)$ and $b = 0$, $H_0^{\uparrow}$ holds whenever $\sigma_t \le \sigma_0$ for all $t$, and $H_0^{\downarrow}$ whenever $\sigma_t \ge \sigma_0$ for all $t$.

*Proof.* For each fixed $\varepsilon$, $s \mapsto \min(s^2\varepsilon^2, \tau^2)$ is nondecreasing. Since $\sigma_t$ and $\tau$ are $\mathcal{F}_{t-1}$-measurable and $\varepsilon_t$ is independent of $\mathcal{F}_{t-1}$, $\mu_t = g_\tau(\sigma_t)$, and monotonicity passes through the expectation. $\square$

Lemma 1 holds for every innovation law $G$: only the size of the change in $\mu_t$ for a given change in $\sigma_t$ depends on $G$. Two caveats follow from it. A change in the tail shape alone, with $\sigma_t$ fixed, also moves $\mu_t$ and can trigger an alarm (with $c = 3$, a switch from Gaussian to Student $t$ innovations at the same variance lowers $\mu_t$ and is seen as a volatility decrease, ALT1 [SYNTHETIC]). And without the independence of $\varepsilon_t$ and $\sigma_t$ (for example under a contemporaneous leverage effect) $\mu_t$ need not be monotone in $\sigma_t$.

### 3.2 Bets, Shiryaev-Roberts e-detectors and the false alarm guarantees

Write $Y_t = X_t/m^{\uparrow}$ with $m^{\uparrow} = (1+b)m_0$ for the increase side, and $Y_t = X_t/m^{\downarrow}$ with $m^{\downarrow} = (1-b)m_0$ and $B = \tau^2/m^{\downarrow} > 1$ for the decrease side, so that $Y_t \in [0, B]$. With a predictable bet $\lambda_t$ define the betting factors
$$
f_t^{\uparrow} = 1 + \lambda_t (Y_t - 1), \ \ \lambda_t \in [0, 1), \qquad f_t^{\downarrow} = 1 + \lambda_t (1 - Y_t), \ \ \lambda_t \in \big[0, 1/(B-1)\big).
$$
Since $Y_t \ge 0$, $f_t^{\uparrow} \ge 1 - \lambda_t > 0$, and since $Y_t \le B$, $f_t^{\downarrow} \ge 1 - \lambda_t(B-1) > 0$. Under the corresponding null, $\mathbb{E}[f_t \mid \mathcal{F}_{t-1}] = 1 + \lambda_t(\mathbb{E}[Y_t \mid \mathcal{F}_{t-1}] - 1) \le 1$ for the increase side, and symmetrically for the decrease side. Hence for every start time $j$ the wealth $L_t^{(j)} = \prod_{s=j}^{t} f_s$ (with $L_{j-1}^{(j)} = 1$) is a nonnegative supermartingale with initial value 1, a *test supermartingale*.

**Definition (Shiryaev-Roberts e-detector).** Fix bets $\lambda_1, \dots, \lambda_K$ in the admissible range and weights $w_k > 0$ with $\sum_k w_k = 1$. For each $k$ let $M_t^{(k)} = f_t^{(k)}(M_{t-1}^{(k)} + 1)$ with $M_0^{(k)} = 0$, so that $M_t^{(k)} = \sum_{j=1}^{t} \prod_{s=j}^{t} f_s^{(k)}$ is the total wealth of bettors who start at every time $j \le t$. Let $M_t = \sum_k w_k M_t^{(k)}$ and raise an alarm at $T_A = \inf\lbrace t \ge 1 : M_t \ge A \rbrace$.

The default grids are $\lambda \in \lbrace 0.02, 0.05, 0.1, 0.2, 0.35, 0.5 \rbrace$ for increases and $\lambda/\lambda_{\max} \in \lbrace 0.1, 0.25, 0.5, 0.75, 0.95 \rbrace$ with $\lambda_{\max} = 1/(B-1)$ for decreases, with equal weights. The updates cost $O(K)$ per step (about 10 microseconds per step for the whole streaming monitor in plain Python and NumPy, E9) [SYNTHETIC].

**Lemma 2.** Under the null, $\mathbb{E}[M_t \mid \mathcal{F}_{t-1}] \le M_{t-1} + 1$ for every $t \ge 1$, so $M_t - t$ is a supermartingale. The same holds for any predictable bets, and for any average of such statistics with fixed weights.

*Proof.* $\mathbb{E}[M_t^{(k)} \mid \mathcal{F}_{t-1}] = \mathbb{E}[f_t^{(k)} \mid \mathcal{F}_{t-1}](M_{t-1}^{(k)} + 1) \le M_{t-1}^{(k)} + 1$; average over $k$. $\square$

**Theorem A (average run length).** Under the null, $\mathbb{E}[T_A] \ge A$ for every $A > 0$.

*Proof.* If $\mathbb{E}[T_A] = \infty$ there is nothing to show; otherwise $T_A < \infty$ almost surely. For fixed $t$, optional stopping of the supermartingale $M_t - t$ (Lemma 2) at the bounded time $T_A \wedge t$ gives $\mathbb{E}[M_{T_A \wedge t}] \le \mathbb{E}[T_A \wedge t] \le \mathbb{E}[T_A]$. As $t \to \infty$, $M_{T_A \wedge t} \to M_{T_A} \ge A$ almost surely, and since $M \ge 0$, Fatou's lemma gives $A \le \mathbb{E}[M_{T_A}] \le \liminf_t \mathbb{E}[M_{T_A \wedge t}] \le \mathbb{E}[T_A]$. $\square$

Theorem A is the run-length guarantee of e-detectors (Shin, Ramdas and Rinaldo 2024), who define an e-detector by $\mathbb{E}[M_\rho] \le \mathbb{E}[\rho]$ for all stopping times $\rho$; we include the proof to make the conditions explicit. A practitioner usually needs a statement for a finite horizon, such as one trading day.

**Theorem B (false alarm within a horizon).** Under the null, $\mathbb{P}(T_A \le m) \le m/A$ for every integer $m \ge 1$.

*Proof.* For $0 \le t \le m$ let $S_t = M_t + (m - t)$. By Lemma 2, $S$ is a supermartingale on $\lbrace 0, \dots, m \rbrace$; it is nonnegative and $S_0 = m$. Ville's inequality gives $\mathbb{P}(\max_{t \le m} S_t \ge A) \le m/A$, and $M_t \ge A$ implies $S_t \ge A$. $\square$

The argument is elementary and we do not claim it is new. Choosing $A = 7800$ for a 390-step trading day, or $A = 28{,}800$ for a 1440-minute crypto day, bounds the probability of any false alarm per day by 5%. The CUSUM-type statistic $C_t = f_t \max(C_{t-1}, 1)$ satisfies $C_t \le M_t$, so it alarms no earlier and inherits both guarantees.

**Proposition 1 (two-sided monitoring).** Under the two-sided null $(1-b)m_0 \le \mu_t \le (1+b)m_0$ for all $t$, the average $\bar M_t = (M_t^{\uparrow} + M_t^{\downarrow})/2$ satisfies Lemma 2, and $\max(M_t^{\uparrow}, M_t^{\downarrow}) \ge 2A$ implies $\bar M_t \ge A$. Hence the rule that alarms at the first $t$ with $\max(M_t^{\uparrow}, M_t^{\downarrow}) \ge 2A$ satisfies Theorems A and B, and the side that crosses gives the direction.

*Proof.* Both sides satisfy Lemma 2 under the two-sided null, so their average does; the inclusion is immediate since both statistics are nonnegative. $\square$

**Corollary C (restarts after alarms).** Suppose the monitor restarts after every alarm, with a new truncation level and baseline that are measurable at the start of the new segment, and that the two-sided null holds in every segment relative to that segment's baseline. Let $N_m$ be the number of alarms within the first $m$ monitored steps (calibration steps are not counted). Then $\mathbb{E}[N_m] \le m/A$.

*Proof sketch.* Let $W_t$ be the sum of $\bar M$ at all past alarms, plus the current $\bar M$, minus the number of monitored steps so far. By Lemma 2 applied in each segment, $W$ is a supermartingale with $W_0 = 0$ (at an alarm the statistic moves into the sum unchanged, and calibration steps leave $W$ unchanged). Each alarm adds at least $A$ to the sum. Optional stopping at the bounded time $\rho \wedge n$, with $\rho$ the time of the $m$-th monitored step, and monotone convergence as $n \to \infty$ give $A\,\mathbb{E}[N_m] - m \le 0$. $\square$

Numerically (E4e) [SYNTHETIC]: with $A = 300$ the true baseline gave 26.5 alarms per 10,000 monitored steps and the streaming monitor with a plug-in baseline recalibrated on 390 steps gave 30.1, against the bound of 33.3. Corollary C is the bookkeeping behind "false alarms per day" over a long series; it is not a single anytime-valid statement for the whole series, and it assumes the true baseline in every segment (Corollary F below covers an estimated one).

### 3.3 Why a decrease detector needs a bounded statistic

**Proposition 2.** Let $\mathcal{H}$ be the set of all probability laws on $[0, \infty)$ with mean at least 1. If $E: [0, \infty) \to [0, \infty]$ is measurable and $\mathbb{E}_P[E(Y)] \le 1$ for every $P \in \mathcal{H}$, then $E(y) \le 1$ for every $y \ge 0$.

*Proof.* For $y \ge 1$ the point mass $\delta_y$ lies in $\mathcal{H}$, so $E(y) \le 1$. For $y < 1$ and $\varepsilon \in (0, 1)$, let $z = (1 - (1-\varepsilon)y)/\varepsilon$, so that $P_\varepsilon = (1-\varepsilon)\delta_y + \varepsilon \delta_z$ has mean exactly 1 and lies in $\mathcal{H}$. Then $(1-\varepsilon)E(y) + \varepsilon E(z) \le 1$, and since $E \ge 0$, $E(y) \le 1/(1-\varepsilon)$. Let $\varepsilon \to 0$. $\square$

The sequential version is what matters for monitoring.

**Corollary 3 (no e-process for "the variance has not fallen").** Let $\mathcal{P}^{\downarrow}$ be the set of laws of real sequences $(r_t)_{t \ge 1}$ with $\mathbb{E}[r_t^2 \mid r_1, \dots, r_{t-1}] \ge m_0$ for every $t$. Let $W_t = W_t(r_1, \dots, r_t) \ge 0$ satisfy $\mathbb{E}_P[W_t] \le 1$ for every $P \in \mathcal{P}^{\downarrow}$ and every fixed $t$ (a requirement weaker than being an e-process). Then $W_t \le 1$ for every $t$ and every path. Consequently a Shiryaev-Roberts statistic built from such processes started at each time satisfies $M_t \le t$ on every path, so it cannot alarm before time $A$, whatever the data.

*Proof.* Fix $t$, a path $(x_1, \dots, x_t)$ and $\varepsilon \in (0, 1)$. Build $P_\varepsilon \in \mathcal{P}^{\downarrow}$ step by step. If the past agrees with $(x_1, \dots, x_{s-1})$ and $x_s^2 \ge m_0$, draw $r_s = x_s$; if $x_s^2 < m_0$, draw $r_s = x_s$ with probability $1 - \varepsilon$ and $r_s = z_s$ otherwise, with $z_s^2 = (m_0 - (1-\varepsilon)x_s^2)/\varepsilon$, so that the conditional mean of $r_s^2$ is exactly $m_0$. Off the path, and after time $t$, draw $r_s = \sqrt{m_0}$. Then $P_\varepsilon(r_1 = x_1, \dots, r_t = x_t) \ge (1-\varepsilon)^t$ and, since $W_t \ge 0$, $1 \ge \mathbb{E}_{P_\varepsilon}[W_t] \ge (1-\varepsilon)^t\, W_t(x_1, \dots, x_t)$. Let $\varepsilon \to 0$. For the last claim, each started process is at most 1, so $M_t \le t$. $\square$

This is a one-sided relative of the classical result of Bahadur and Savage (1956), who showed that without restrictions on the tails the mean of a distribution admits no nontrivial test. The asymmetry is what matters here. The *increase* null "mean at most 1" on $[0, \infty)$ does admit nontrivial e-variables, for example $E(y) = y$, which is why an increase detector works on raw squared returns (Section 5). A fall in the variance itself, however, cannot be detected by any nonparametric e-test, whatever statistic it computes. A decrease detector must therefore change its target to a bounded functional of the returns. For a statistic with values in $[0, B]$ the null "conditional mean at least 1" does admit nontrivial e-variables, namely the bets $1 + \lambda(1 - y)$ with $\lambda \in [0, 1/(B-1)]$, and up to domination these are the only ones (Lemma 5 in Section 4). Truncation is one way to obtain a bounded functional, and the truncated variance $\mathbb{E}[\min(r^2, \tau^2)]$ falls exactly when $\sigma$ falls under the fixed-shape assumption of Lemma 1. It is not the only way: exceedance indicators $1\lbrace \lvert r_t \rvert > q \rbrace$, ranks, or any bounded increasing transform of $r_t^2$ would also do (a rank-based detector, exactly valid for exchangeable returns, is studied in ALT5). Truncation is the one with a direct volatility meaning and a known bias (Section 2). The price is a speed limit: one step can multiply the wealth of a decrease bet by at most $B/(B-1)$ (Corollary 7). For increases, clipping is not needed for validity, but it bounds the effect of a single jump on the wealth, $f_t^{\uparrow} \le 1 + \lambda(B-1)$, which is the betting counterpart of the bounded influence of the clipped estimator.

### 3.4 An estimated baseline

The guarantees above assume that the baseline $m_0$ and the level $\tau$ are fixed in advance, or at least that the null holds with whatever values the monitor uses. In practice both are estimated from a calibration window. The *plug-in* baseline (the mean of $X$ over the window, with $\tau$ from the same window) is not covered, and its effect is measurable [SYNTHETIC]:

* With 390 calibration steps, the one-day false alarm rate of the increase side was 4.58%, 4.68% and 5.10% (Gaussian, Student $t_3$, jumps), against 1.75%, 1.73% and 1.62% with the true baseline; with 1440 steps it was 2.75%, 2.55% and 2.00% (4000 paths, Monte Carlo standard error 0.34 points, E4e). The plug-in roughly triples the rate at short calibration, and only the slack of Ville's inequality keeps it near the bound.
* The excess grows with the horizon, because a baseline that is slightly too low gives the log wealth a small positive drift that the bets eventually find. With exactly i.i.d. standardised returns and a one-day plug-in baseline, the false alarm rate over five days was 44% (95% interval 39% to 50%), where Theorem B with the true baseline guarantees 25% (E10, Q10b2). Over seven days at band 0 it was 51% to 57% for Gaussian, Student $t_4$ and jump returns, against a guaranteed 35% (ALT5).

A confidence sequence restores a proof. The construction uses a sample split, which is needed because a $\tau$ computed from the whole window makes each $X_i$ in the window depend on later returns.

**Construction.** The first $n_\tau$ returns of the calibration window fix $\tau = c\,\hat s$, so $\tau$ is $\mathcal{F}_{n_\tau}$-measurable (the default uses 25% of the window). For the later returns, indexed $i = 1, 2, \dots$, let $Z_i = \min(r_i^2, \tau^2)/\tau^2 \in [0, 1]$ and $\eta_i = \mathbb{E}[Z_i \mid \mathcal{F}_{i-1}]$. Choose predictably a centre $\hat z_{i-1} \in (0, 1)$ and bets $\lambda_i^{\uparrow}, \lambda_i^{\downarrow} > 0$ with $\lambda_i^{\uparrow}(1 - \hat z_{i-1}) \le c_0$ and $\lambda_i^{\downarrow}\hat z_{i-1} \le c_0$ for a fixed $c_0 \in (0, 1)$. With $\psi_a(\lambda) = (-\log(1 - \lambda a) - \lambda a)/a^2$, $e_i = (Z_i - \hat z_{i-1})^2$, $\psi_i^{\uparrow} = \psi_{1-\hat z_{i-1}}(\lambda_i^{\uparrow})$, $\psi_i^{\downarrow} = \psi_{\hat z_{i-1}}(\lambda_i^{\downarrow})$ and error levels $\delta_u, \delta_l \in (0, 1)$, put
$$
U_t = \frac{\sum_{i \le t} \lambda_i^{\uparrow} Z_i + \log(1/\delta_u) + \sum_{i \le t} \psi_i^{\uparrow} e_i}{\sum_{i \le t} \lambda_i^{\uparrow}}, \qquad
L_t = \frac{\sum_{i \le t} \lambda_i^{\downarrow} Z_i - \log(1/\delta_l) - \sum_{i \le t} \psi_i^{\downarrow} e_i}{\sum_{i \le t} \lambda_i^{\downarrow}}.
$$
They bound the bet-weighted averages $\bar\eta_t^{\uparrow} = \sum_{i \le t}\lambda_i^{\uparrow}\eta_i / \sum_{i \le t}\lambda_i^{\uparrow}$ and $\bar\eta_t^{\downarrow}$ (defined in the same way with $\lambda^{\downarrow}$) at all times at once: $\mathbb{P}(\exists t: \bar\eta_t^{\uparrow} \ge U_t) \le \delta_u$ and $\mathbb{P}(\exists t: \bar\eta_t^{\downarrow} \le L_t) \le \delta_l$, with no assumption on dependence or on how $\eta_i$ varies (Theorem D.1 in the appendix). With $1 - \hat z = \hat z = 1$ these are the predictable plug-in empirical Bernstein bounds of Waudby-Smith and Ramdas (2024, Theorem 2); the one-sided ranges let the lower bound use bets up to $c_0/\hat z$, which matters because $\hat z \approx 0.11$ for a clipped squared Gaussian return.

**Monitoring rule.** Calibration ends at a stopping time $T_0 > n_\tau$ (a fixed length, or a rule that stops when the relative width of the bounds falls below a target). Set $m^{\uparrow} = (1+b)\,\tau^2 \min(U_{T_0}, 1)$ and $m^{\downarrow} = (1-b)\,\tau^2 \max(L_{T_0}, 0)$, run the two Shiryaev-Roberts statistics of Section 3.2 from $T_0$ with these baselines and bets fixed at $T_0$ (the decrease side is switched off when $m^{\downarrow} = 0$), and alarm at the first $t > T_0$ with $\max(M_t^{\uparrow}, M_t^{\downarrow}) \ge 2A$.

**Theorem D (guarantee with an estimated baseline).** Let $\delta = \delta_u + \delta_l$ and assume
$$
(H^{\uparrow})\ \ \mu_t \le (1+b)\,\tau^2\,\bar\eta_{T_0}^{\uparrow} \quad\text{and}\quad (H^{\downarrow})\ \ \mu_t \ge (1-b)\,\tau^2\,\bar\eta_{T_0}^{\downarrow} \qquad \text{for all } t > T_0,
$$
with nothing assumed about the calibration window beyond $0 \le X \le \tau^2$. Then

1. for every integer $m \ge 1$, $\mathbb{P}(\text{an alarm at one of } T_0 + 1, \dots, T_0 + m) \le \delta + m/A$, and more sharply $\le \delta + (1-\delta)\,m/A$ when $m < A$;
2. if $T_0 < \infty$ almost surely, the expected number of monitored steps before the first alarm is at least $(1-\delta)A$;
3. one side alone at level $A$ satisfies the bound with $\delta_u$ (increase, under $H^{\uparrow}$) or $\delta_l$ (decrease, under $H^{\downarrow}$) in place of $\delta$.

*Proof sketch* (full proof in Appendix A). Let $G$ be the event that both bounds hold at $T_0$; it is known at $T_0$, and because the coverage of Theorem D.1 is uniform in time, $\mathbb{P}(G^c) \le \delta$ for every stopping time $T_0$. On $G$, $(H^{\uparrow})$ and $(H^{\downarrow})$ imply $\mu_t \le m^{\uparrow}$ and $\mu_t \ge m^{\downarrow}$, so on $G \cap \lbrace t > T_0 \rbrace$ both statistics satisfy Lemma 2. Applying Ville's inequality to $1_{G \cap \lbrace T_0 = s\rbrace}(\bar M_t + m - (t - s))$ for each $s$ and summing over $s$ gives $\mathbb{P}(G \cap \text{alarm within } m) \le \mathbb{P}(G)\,m/A$. Add $\mathbb{P}(G^c) \le \delta$. $\square$

**Corollary E (a constant conditional law).** If from step $n_\tau + 1$ until the change the conditional law of $r_t$ given $\mathcal{F}_{t-1}$ is one fixed law $P_0$ (for example i.i.d. returns with jumps and fat tails), then $\mu_t = \mathbb{E}_{P_0}[\min(r^2, \tau^2)]$ is constant, $(H^{\uparrow})$ and $(H^{\downarrow})$ hold with $b = 0$, and with no change $\mathbb{P}(\text{false alarm within } m \text{ monitored steps}) \le \delta + m/A$ with the baseline estimated from the data.

**Corollary F (recalibration after every alarm).** If every segment has its own split, bounds and error budget $\delta$, and $(H^{\uparrow})$, $(H^{\downarrow})$ hold in every segment relative to its own calibration averages, then the number $N_m$ of alarms within the first $m$ monitored steps satisfies $\mathbb{E}[N_m] \le (\delta + m/A)/(1-\delta)$. With $\delta = 0$ this is Corollary C. (Proof in Appendix A.)

**Choosing $A$ and $\delta$.** For a target $\alpha$ per day of $m$ steps, choose $\delta < \alpha$ and $A = m/(\alpha - \delta)$. With $\delta = 0.01$ and $\alpha = 5\%$, $A = 9750$ for a 390-step day and $A = 36{,}000$ for a 1440-minute day. Because the delay grows like $\log A$, the larger $A$ costs little ($\log 9750/\log 7800 \approx 1.025$); the real cost is the width of the bounds (Section 5.3). Two remarks. Under volatility clustering the bounds still cover the bet-weighted average of the conditional means within the calibration window, but the monitoring null "never above that average" is then simply false (Section 5.4): the confidence sequence protects against estimation error, not against this model mismatch. And because the coverage is uniform in time, the calibration may end at any stopping time, which a fixed-$n$ interval would not allow.

