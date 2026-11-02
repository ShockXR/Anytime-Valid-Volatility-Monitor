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

## 4 Detection delay

Section 3 controls false alarms. This section asks how fast the detector reacts to a real change, what robustness costs in speed, and what limits the speed of decrease detection. Proofs not given here are in Appendix B. Numerical statements are from E10 and E11 and are [SYNTHETIC].

### 4.1 Growth rate and information

**Post-change model.** From a change time $\nu$ on, $Y_\nu, Y_{\nu+1}, \dots$ are i.i.d. with law $Q$ and independent of $\mathcal{F}_{\nu-1}$. Write $V = Y - 1$ for the increase side and $V = 1 - Y$ for the decrease side; $V$ takes values in $[v_{\mathrm{lo}}, v_{\mathrm{hi}}]$, equal to $[-1, B-1]$ or $[1-B, 1]$, and the admissible bets are $\lambda \in [0, \lambda_{\max})$ with $\lambda_{\max} = 1/\lvert v_{\mathrm{lo}}\rvert$. A change means $\mathbb{E}_Q[V] > 0$. The worst-case delay of a stopping rule $T$ is Lorden's (1971)
$$
\mathrm{WADD}(T) = \sup_{\nu \ge 1}\ \operatorname*{ess\,sup}\ \mathbb{E}_\nu\big[(T - \nu + 1)^+ \mid \mathcal{F}_{\nu-1}\big].
$$
The *growth rate* of a bet is $g(\lambda) = \mathbb{E}_Q[\log(1 + \lambda V)]$. It is strictly concave with $g(0) = 0$ and $g'(0) = \mathbb{E}_Q[V] > 0$, so it has a unique maximiser $\lambda^\ast(Q)$ in $[0, \lambda_{\max}]$. Let $\mathcal{H} = \lbrace P \text{ on } [0, B] : \mathbb{E}_P[V] \le 0 \rbrace$ be the nonparametric null class of one step, and
$$
\mathrm{KL}_{\inf}(Q) = \inf_{P \in \mathcal{H}} \mathrm{KL}(Q \,\Vert\, P).
$$

**Lemma 4 (growth is bounded by information).** If $E \ge 0$ and $\mathbb{E}_P[E(Y)] \le 1$ for a law $P$, then $\mathbb{E}_Q[\log E(Y)] \le \mathrm{KL}(Q \,\Vert\, P)$.

**Lemma 5 (the best nonparametric growth).** Let $B < \infty$. (a) Every e-variable for the whole class $\mathcal{H}$ is dominated by a bet: there is $\lambda \in [0, \lambda_{\max}]$ with $E(y) \le 1 + \lambda v(y)$ for all $y$. (b) If $0 < \lambda^\ast < \lambda_{\max}$, then $\max_\lambda g(\lambda) = g(\lambda^\ast) = \mathrm{KL}_{\inf}(Q)$, and the infimum is attained by $P^\ast(dy) = Q(dy)/(1 + \lambda^\ast v(y))$.

Part (b) is the dual form, for a bounded statistic, of the information number $K_{\inf}$ of Honda and Takemura (2010), and $P^\ast$ is the reverse information projection of $Q$ onto the null (Grünwald, de Heide and Koolen 2024; Larsson, Ramdas and Ruf 2024). Since every bet is an e-variable for the true pre-change law $P_0$ (when $P_0 \in \mathcal{H}$), and since $X$ and $Y$ are functions of $r$, Lemmas 4 and 5 and the data processing inequality give:

**Corollary 6 (information ladder).** With $g_{\mathrm{grid}} = \max_k g(\lambda_k)$,
$$
\mathrm{KL}(Q_r \,\Vert\, P_{0,r}) \ \ge\ \mathrm{KL}(Q_X \,\Vert\, P_{0,X}) \ \ge\ \mathrm{KL}_{\inf}(Q) \ \ge\ g_{\mathrm{grid}}.
$$
The three ratios between neighbouring terms are the *truncation loss*, the *nonparametric loss* (the price of a null that restricts only the mean of $Y$, which is what makes the detector valid for every return law) and the *grid loss*. The ratio $\mathrm{KL}(Q_r \Vert P_{0,r})/g_{\mathrm{grid}}$ is the first-order delay ratio of the robust detector against the likelihood-ratio detector that knows both laws.

**Corollary 7 (speed limits).** Every bet satisfies $f_\lambda \le 1 + \lambda_{\max} v_{\mathrm{hi}}$, so by Lemma 5(a) every e-variable for $\mathcal{H}$ satisfies $\log E \le \log B$ on the increase side and $\log E \le \log(B/(B-1))$ on the decrease side, and $\mathrm{KL}_{\inf}(Q)$ obeys the same bounds however large the change. For the detector itself the statement is deterministic: $M_n \le \sum_k w_k \sum_{j=1}^{n} F_k^{\,n-j+1}$ with $F_k = 1 + \lambda_k v_{\mathrm{hi}}$, so $T_A$ is never smaller than the first $n$ at which this sum reaches $A$.

With $c = 3$ and Gaussian returns, $\log(B/(B-1)) = 0.117$ nats per step, and for the default decrease grid at $A = 7800$ the floor is 36, 73, 122 and 180 steps for $c = 2, 3, 4$ and 5 (E11). A complete collapse of volatility (all returns zero) is detected after exactly 73 steps with $c = 3$. This is the quantitative form of the impossibility of Section 3.3: the bound $B$ that makes a decrease detectable also caps the evidence that one step can provide.

### 4.2 Delay bounds and first-order optimality

**Theorem 8 (non-asymptotic delay bound).** Under the post-change model, for the Shiryaev-Roberts mixture or its CUSUM-type version:

(a) for every $A$, every change time $\nu$ and every pre-change history,
$$
\mathbb{E}\big[(T_A - \nu + 1)^+ \mid \mathcal{F}_{\nu-1}\big] \ \le\ \min_{k:\, g_k > 0}\ \frac{\log(A/w_k) + o_k}{g_k},
$$
where $g_k = g(\lambda_k)$, $Z_k = \log(1 + \lambda_k V)$ with $V$ drawn from $Q$, and $o_k = \min\big(\sup Z_k,\ \mathbb{E}[(Z_k^+)^2]/\mathbb{E}[Z_k]\big)$ bounds the mean overshoot;

(b) $\mathrm{WADD}(T_A) = \mathbb{E}_Q[T_A^0]$, where $T_A^0$ is the alarm time of the detector started from $M_0 = 0$ on i.i.d. data from $Q$: the worst case is a change at the first monitored step.

*Proof.* (a) For $t \ge \nu$ all terms of the Shiryaev-Roberts sum are nonnegative, so $M_t \ge w_k L_t^{(\nu, k)}$. Let $S_n = \sum_{i=0}^{n-1} Z_{k, \nu+i}$ and $N = \inf\lbrace n \ge 1 : S_n \ge \log(A/w_k)\rbrace$; then $(T_A - \nu + 1)^+ \le N$. $N$ depends only on post-change data, which are independent of $\mathcal{F}_{\nu-1}$. $S$ is a random walk with drift $g_k > 0$ and bounded increments, so Wald's identity gives $g_k\,\mathbb{E}[N] = \log(A/w_k) + \mathbb{E}[S_N - \log(A/w_k)]$, and the overshoot is bounded by $\sup Z_k$ and, in mean, by Lorden's (1970) inequality. (b) Each recursion $x \mapsto f(x + 1)$ (or $f\max(x, 1)$) is nondecreasing in $x$, so from any state $M_{\nu-1} \ge 0$ the statistic after $\nu$ dominates the one started from 0 on the same post-change data. $\square$

Part (a) is the elementary version, for linear bets on a bounded statistic, of the delay bounds that Shin, Ramdas and Rinaldo (2024) give for mixtures of exponential baseline e-detectors. It pays the full $\log(1/w_k) = \log K$ for not knowing the right bet and ignores the help of the other bets and start times; in simulation it is 1.3 to 2.7 times the true delay (E11). Part (b) means that simulations with the change at the first monitored step measure the worst-case delay exactly, and any pre-change history can only shorten it.

**Theorem 9 (first-order asymptotics).** Fix $Q$ with $\mathbb{E}_Q[V] > 0$ and an interior maximiser $0 < \lambda^\ast < \lambda_{\max}$, so that $g(\lambda^\ast) = \mathrm{KL}_{\inf}(Q)$. As $A \to \infty$:

(a) *Lower bound.* Let $T$ be any stopping rule based on the returns whose average run length is at least $A$ under every i.i.d. return law $P$ with $\mathbb{E}_P[\min(r^2, \tau^2)] \le m^{\uparrow}$ (increase side; symmetrically for the decrease side). Then $\mathrm{WADD}_Q(T) \ge (\log A / \mathrm{KL}_{\inf}(Q))(1 + o(1))$.

(b) *Fixed grid.* For a fixed grid with $g_{\mathrm{grid}} > 0$, $\mathrm{WADD}(T_A) \big/ (\log A/\mathrm{KL}_{\inf}) \to \mathrm{KL}_{\inf}/g_{\mathrm{grid}} \ge 1$.

(c) *Refined grid.* Fix $0 < a < b' < \lambda_{\max}$ with $\lambda^\ast \in [a, b']$, and use $K_A = 1 + \lceil (b' - a)\sqrt{\log A}\,\rceil$ equally spaced bets on $[a, b']$ with equal weights. Then $\mathrm{WADD}(T_A) \le (\log A + \tfrac12 \log\log A + O(1))/\mathrm{KL}_{\inf}$.

(d) *Continuous mixture.* For a mixing law on $[0, b']$ with a density at least $p_0 > 0$ on $[a, b']$ and $\lambda^\ast \in (a, b')$, $\mathrm{WADD}(T_A) \le (\log A + \log\log A + O(1))/\mathrm{KL}_{\inf}$.

Hence the Shiryaev-Roberts e-detector with the oracle bet, a refined grid or a continuous mixture is first-order minimax optimal among all rules that control the run length over the nonparametric null, and a fixed grid has asymptotic efficiency $g_{\mathrm{grid}}/\mathrm{KL}_{\inf}$.

*Proof sketch* (details in Appendix B). (a) The law $P_r^\ast(dr) = Q_r(dr)/(1 + \lambda^\ast v(r))$ is a probability law (its mass is 1 by the first-order condition $\mathbb{E}_Q[V/(1 + \lambda^\ast V)] = 0$), it satisfies the null with equality, and $\mathrm{KL}(Q_r \Vert P_r^\ast) = g(\lambda^\ast)$. The rule has run length at least $A$ under $P_r^\ast$, and Lorden's (1971) lower bound for a change from $P_r^\ast$ to $Q_r$ gives the claim. (b) The upper half is Theorem 8(a) with the best bet; the lower half uses $\log M_n \le \log n + \max_{k} \max_{i < n}(S_n^k - S_i^k)$ and the strong law of large numbers. (c) Some bet lies within half a grid step of $\lambda^\ast$, where $g$ is within $O(h^2)$ of its maximum (Taylor), and Theorem 8(a) with $w_k = 1/K_A$ gives the bound. (d) Concavity of $\lambda \mapsto \log(1 + \lambda v)$ bounds the mixture from below by its value on a window of width $2/\log A$ around $\lambda^\ast$, and the Wald argument applies to the resulting random walk. $\square$

Ram and Ramdas (2026, preprint) prove a universal lower bound of the form (a) for composite pre- and post-change classes and a matching upper bound for bounded means, which covers the truncated statistic. We give (a) to (d) because for linear bets on a bounded statistic the arguments are short, and because the comparison between fixed grids, refined grids and mixtures is what matters for the choice of bets.

### 4.3 Second-order behaviour at practical alarm levels [HEURISTIC, checked numerically]

The first-order formula $\log A / g$ is a poor guide at $A = 7800$. For one bet and the change at the first step, $M_n = e^{S_n} V_n$ with $S_n = \sum_{i \le n} \log f_i$ and $V_n = \sum_{j=0}^{n-1} e^{-S_j} \uparrow V_\infty < \infty$, so $\log M_n$ is a random walk plus a slowly changing term. Nonlinear renewal theory (Woodroofe 1982; Pollak 1987 for the Shiryaev-Roberts statistic) suggests
$$
\mathbb{E}[T_A] \approx \frac{\log A + \kappa - \zeta}{g}, \qquad \kappa = \frac{\mathbb{E}[H^2]}{2\,\mathbb{E}[H]}, \quad \zeta = \mathbb{E}[\log V_\infty],
$$
with $H$ the first ascending ladder height of $S$ (so $\kappa$ is the limiting mean overshoot). We have not checked the regularity conditions of those theorems for our law of $Y$ (a continuous part and an atom at the cap), so we call the formula heuristic. Approximating $S$ by a Brownian motion with drift $g$ and variance $s^2$ per step, Dufresne's (1990) identity for exponential functionals gives $\zeta \approx \log(2/s^2) - \psi(2g/s^2)$, with $\psi$ the digamma function. For a growth-optimal bet and a small change $s^2 \approx 2g$, so $\zeta \approx \log(1/g) + \gamma_E$ with Euler's constant $\gamma_E = 0.5772$, and
$$
\mathbb{E}[\text{delay}] \approx \frac{\log(A g) - \gamma_E + \kappa}{g}.
$$
The memory term $\zeta$ is the reason the first-order formula fails: for a 25% increase in standard deviation $g \approx 0.05$, so $\zeta \approx 3.6$ while $\log A = 8.96$.

Table 1 compares the predictions with simulation (true baseline, $c = 3$, $A = 7800$, change at the first monitored step, 3000 paths per case, mean delays, E11) [SYNTHETIC].

| Case | $g_{\mathrm{grid}}$ | first order $\log A/g_{\mathrm{grid}}$ | proved bound (Thm. 8a) | simulated, default grid | second order, best single bet | simulated, best single bet |
|---|---|---|---|---|---|---|
| Gaussian, sd ratio 1.1 | 0.0093 | 964 | 1221 | 455 | 405 | 417 |
| Gaussian, sd ratio 1.25 | 0.0502 | 179 | 233 | 122 | 109 | 107 |
| Gaussian, sd ratio 1.5 | 0.161 | 56 | 75 | 48 | 42 | 42 |
| Gaussian, sd ratio 2 | 0.413 | 22 | 30 | 23 | 18.4 | 18.7 |
| Jumps, sd ratio 1.25 | 0.0473 | 190 | 248 | 128 | 115 | 115 |
| Student $t_5$, sd ratio 1.25 | 0.0320 | 280 | 362 | 178 | 158 | 159 |
| Gaussian, sd ratio 0.8 | 0.0340 | 264 | 314 | 195 | 157 | 156 |
| Gaussian, sd ratio 0.67 | 0.0598 | 150 | 179 | 124 | 99 | 99 |
| Gaussian, sd ratio 0.5 | 0.0839 | 107 | 127 | 93 | 74 | 74 |

*Table 1. Predicted and simulated mean delays in steps (E11, Q11a). Standard errors are below 1.5% of the mean.*

In all 20 cases of E11 the second-order prediction for a single bet is within 3% of its simulated mean delay, the first-order formula overstates the delay of small changes by up to 112%, and the proved bound is 1.3 to 2.7 times the simulated delay. The default grid costs 9% to 26% against the best single bet. With the change after 500 pre-change steps, as in E4, the delays are 17% to 34% shorter, because the Shiryaev-Roberts statistic enters the change with a head start, and they match E4 to within about 10% (E11). The asymptotic regime arrives very slowly: at $A = 10^8$ even the best single bet is 5% to 24% faster than $\log A/\mathrm{KL}_{\inf}$ in 15 cases (Gaussian, Student $t_3$ and jump returns), and the second-order formula explains the gap to within 0.024 in ratio for every $A \ge 10^4$ (E11, Q11e; Figure `e11_asymptotic_optimality.png`).

![Figure 1. Ratio of predicted to simulated delay for the first-order formula, the proved bound of Theorem 8(a) and the second-order formula (E11) [SYNTHETIC].](../results/figures/e11_delay_predictions.png)

### 4.4 The price of robustness and the truncation constant

Table 2 splits the first-order delay ratio against the Gaussian likelihood-ratio detector into the three factors of Corollary 6, for clean Gaussian data and $c = 3$, and compares it with simulated mean delays (E11, Q11b) [SYNTHETIC].

| sd ratio | 1.25 | 1.5 | 2 | 0.8 | 0.67 | 0.5 |
|---|---|---|---|---|---|---|
| truncation loss | 1.03 | 1.08 | 1.27 | 1.00 | 1.00 | 1.00 |
| nonparametric loss | 1.10 | 1.23 | 1.48 | 1.25 | 2.00 | 3.61 |
| grid loss | 1.02 | 1.02 | 1.04 | 1.02 | 1.05 | 1.05 |
| predicted penalty $\mathrm{KL}/g_{\mathrm{grid}}$ | 1.16 | 1.37 | 1.95 | 1.27 | 2.09 | 3.79 |
| simulated penalty against the Gaussian LR mixture | 1.14 | 1.35 | 1.73 | 1.48 | 2.19 | 3.45 |

*Table 2. The price of robustness on clean Gaussian data, $c = 3$ (E11). Penalties are ratios of mean delays, change at the first step, 3000 paths with common random numbers.*

The information ratio explains the 10% to 70% delay penalty measured in E4 (Section 5.3). The largest factor is the nonparametric loss, not truncation. For decreases the nonparametric loss is the speed limit of Corollary 7: $\mathrm{KL}_{\inf}$ cannot exceed 0.117 nats per step, while the Gaussian information keeps growing with the size of the change. Asymptotically, the robust detector keeps $\mathrm{KL}_{\inf}/\mathrm{KL}$ = 53% to 88% of the speed of the Gaussian likelihood-ratio detector for increases (sd ratio 2 to 1.25) and 28% to 80% for decreases (0.5 to 0.8); under Student $t_3$ returns it keeps 94% to 98% of the speed of the correctly specified detector for increases, because heavy tails carry little information about scale (E11, Q11e). In exchange the robust detector keeps its false alarm guarantee where the Gaussian one loses it (Section 5.1).

**The truncation constant.** For increases on clean data the delay barely changes for $c \ge 4$, and $c = 3$ costs 3%, 6% and 12% at sd ratios 1.25, 1.5 and 2 against no truncation. With jumps (probability 0.005 per step, standard deviation 8 times the diffusive one) the growth at sd ratio 1.25 peaks at $c$ between 3.25 and 4.25 depending on the change size, and without truncation it falls to 29% of the $c = 3$ value. A burst of jumps at constant volatility (jump probability raised from 0.005 to 0.02) moves the $c = 3$ detector only 3.5% as much as a real 25% volatility increase, against 39% without truncation (E11). For decreases the speed limit dominates and a smaller $c$ is faster: with $c = 2$ the simulated decrease delays fall by 26% to 49% relative to $c = 3$. The two sides need not share $c$, since each is a valid e-detector for its own truncated statistic and the two-sided rule of Proposition 1 needs only each side's guarantee: $c$ of about 3 to 4 for increases and 1.5 to 2 for decreases is a legitimate design. With a plug-in baseline and four return laws, $c = 2$ on the decrease side shortened the median decrease delay in all 12 simulated cells, by 4% to 51%, while $c = 1.5$ was slower than $c = 3$ for a 20% fall when returns were rounded to price ticks, because zero and one-tick returns then carry most of the statistic (E11f) [SYNTHETIC]. On the five cryptocurrencies of Section 6 the choice $c = 2$ for decreases was neutral to slightly positive for the confirmation rate of decrease alarms, within sampling noise (E11f) [REAL].

### 4.5 Clustering, forecast-relative monitoring and a signal cap

Under persistent volatility clustering the conditional mean of $X_t$ wanders even when the unconditional variance is constant, so the nulls of Section 3 are false and the raw monitor alarms on most days (Section 5.4). A natural fix is to monitor $z_t = r_t/\hat\sigma_t$, with $\hat\sigma_t^2$ a predictable volatility forecast; the null then reads "the forecast is right up to a constant". The following exact calculation shows what this costs when the forecast model is GARCH(1,1) (Bollerslev 1986).

**Proposition 10 (signal cap).** Let $\alpha, \beta > 0$ with $\alpha + \beta < 1$ and $\omega > 0$. Before a break the returns follow GARCH$(\omega, \alpha, \beta)$, and the forecast is the pre-break filter $s_t = \omega + \alpha r_{t-1}^2 + \beta s_{t-1}$, run without change after the break. From the break on, the returns follow GARCH$(k\omega, \alpha, \beta)$ with $k > 0$: $r_t = \sigma_t \varepsilon_t$, $\sigma_t^2 = k\omega + \alpha r_{t-1}^2 + \beta \sigma_{t-1}^2$, with i.i.d. $\varepsilon_t$ of mean 0 and variance 1. Then $d_t = \sigma_t^2 - s_t$ satisfies $d_t = (k-1)\omega + \beta d_{t-1}$ deterministically, so $d_t \to (k-1)\omega/(1-\beta)$ geometrically at rate $\beta$, and in the post-break stationary regime
$$
\rho(k) := \frac{\mathbb{E}[r_t^2]}{\mathbb{E}[s_t]} = \frac{k(1-\beta)}{1 - \alpha - \beta + \alpha k},
$$
which does not depend on $\omega$, increases in $k$, and satisfies $\rho(k) < (1-\beta)/\alpha$ for every $k$, however large the break.

*Proof.* Subtracting the two recursions, the terms $\alpha r_{t-1}^2$ cancel, which gives the recursion for $d_t$. In the post-break stationary regime $\mathbb{E}[r_t^2] = \mathbb{E}[\sigma_t^2] = k\omega/(1 - \alpha - \beta)$ and $\mathbb{E}[s_t] = \mathbb{E}[\sigma_t^2] - (k-1)\omega/(1-\beta)$. Dividing and simplifying gives $\rho(k)$. Its derivative in $k$ has the sign of $(1-\beta)(1-\alpha-\beta) > 0$, and $\rho(k) \to (1-\beta)/\alpha$ as $k \to \infty$. $\square$

Two exact facts sit behind this. A scale break ($r_t \mapsto \sqrt{k}\, r_t$) and an intercept break ($\omega \mapsto k\omega$) lead to the same post-break model and differ only in the transient right after the break. And the filter absorbs most of the break through its $\alpha r^2$ term: at persistence $\alpha + \beta = 0.99$ with $\alpha = 0.05$, a doubling of the long-run variance appears after standardisation as an increase of the ratio of means to $\rho(2) = 1.091$, and no break can push it above $(1 - \beta)/\alpha = 1.2$. The monitor sees the mean of the ratio, $\mathbb{E}[r_t^2/s_t]$, which is a little larger by Jensen's inequality (1.116 for $k = 2$ in simulation) and which stayed below 1.28 for $k$ up to 1000 at this persistence; we have a proof of the cap only for the ratio of means. Figure 2 compares the closed form with simulation (E10, Q10d) [SYNTHETIC].

![Figure 2. The ratio of means $\rho(k)$ (lines, closed form; filled circles, simulation) and the mean of the ratio that the monitor sees (open triangles), with the caps $(1-\beta)/\alpha$ (dashed), for three persistences. 200 paths of 18,000 post-break steps per point, standard errors below 0.002 (E10) [SYNTHETIC].](figures/signal_cap.png)

**An information rate.** With Gaussian innovations the per-step Kullback-Leibler divergence between the post-break and pre-break conditional laws of $r_t$ is $\tfrac12(q_t - 1 - \ln q_t)$ with $q_t = \sigma_t^2/s_t$. Let $I$ be its average over the post-break path. By Lai (1998), any detector whose average run length to a false alarm is at least $A$ needs, as $A \to \infty$ and under his condition that the post-change log-likelihood ratio grows at rate $I$ uniformly in the change time, a worst-case delay of at least $(\ln A/I)(1 + o(1))$. At $A = 28{,}800$ this is an order-of-magnitude guide, not an exact bound. Table 3 gives the numbers.

| $\alpha + \beta$ | $k$ | $\rho(k)$ | simulated mean of ratio | $I$ (nats per step) | $\ln A / I$ (days of 1440 steps) |
|---|---|---|---|---|---|
| 0.95 | 2 | 1.333 | 1.361 | 0.0276 | 0.26 |
| 0.99 | 1.5 | 1.059 | 1.074 | 0.0015 | 4.7 |
| 0.99 | 2 | 1.091 | 1.116 | 0.0037 | 1.9 |
| 0.99 | 4 | 1.143 | 1.187 | 0.0093 | 0.77 |
| 0.99 | 1000 | 1.200 | 1.271 | 0.0187 | 0.38 |
| 0.999 | 2 | 1.016 | 1.031 | 0.0003 | 21 |

*Table 3. Signal left after a frozen GARCH(1,1) forecast and the information rate of a break (E10, Q10d) [SYNTHETIC]. $A = 28{,}800$.*

At persistence 0.99 no detector with this false alarm budget can be expected to see a doubling of the long-run variance much faster than within about two days, and at persistence 0.999 it is practically invisible for weeks. This is a property of the question, not of the monitor: a highly persistent volatility model can explain a level shift as its own slow wandering, which is the mirror image of the known finding that unmodelled level shifts make fitted GARCH models look nearly integrated (Diebold 1986; Lamoureux and Lastrapes 1990; Mikosch and Starica 2004). The observed median delays of the forecast-relative monitors (1523 to 1791 steps for a doubling, Section 5.4) are about 0.6 times the step at which the cumulative information reaches $\ln A$. The memory term of Section 4.3 can explain a factor of this size on its own: in an i.i.d. version of the same steady state with the true baseline, $\zeta$ is 54% to 62% of $\ln A$ and the detector alarms at 0.41 to 0.58 times $\ln A/\mathrm{KL}$ (E11, cross-check with E10). We did not rerun the GARCH setting itself to separate this from the effect of the plug-in baseline.

## 5 Simulation evidence

All results in this section are **[SYNTHETIC]**: returns are simulated in per-step units with fixed seeds, and the settings of each experiment are collected in Appendix C. Unless stated otherwise the monitor uses $c = 3$ with the scaled MAD, the default grids of Section 3.2 and $A = 7800$, so that with a known baseline Theorem B bounds the probability of a false alarm within one 390-step day by 5%. "Jumps" adds to Gaussian returns, with probability 0.005 per step, an independent normal jump whose standard deviation is 8 times (E4, E8, E11) or 10 times (E9) the diffusive one. Student $t_\nu$ returns are scaled to unit variance. "Price ticks" are Gaussian returns rounded to multiples of half a standard deviation (about 20% zero returns).

### 5.1 Validity where classical detectors fail

E4 compares the two sides of the robust monitor with an e-detector on raw squared returns, a Shiryaev-Roberts mixture of Gaussian likelihood ratios for standard deviation ratios 1.25, 1.5 and 2 (increase side) at the same $A$, Page's (1954) CUSUM on $r_t^2$ with its threshold calibrated by simulation to an average run length of about 7800 on Gaussian data, and a two-sided $z$-test of the variance over the last 60 steps at level 0.05, repeated at every step. The baseline of every method is estimated from a 2000-step calibration window drawn from the same law as the monitored data, so jumps and fat tails are part of "normal".

| Detector | Gaussian | jumps | Student $t_5$ | price noise | slow SV |
|---|---|---|---|---|---|
| Robust e-detector, increase | 1.7% | 2.3% | 0.7% | 3.7% | 26.3% |
| Robust e-detector, decrease | 0.7% | 0.7% | 0.3% | 0.3% | 27.0% |
| E-detector on raw $r_t^2$, increase | 1.7% | 2.3% | 1.3% | 3.0% | 26.3% |
| Gaussian likelihood-ratio e-detector, increase | 2.7% | **70.7%** | **96.3%** | 6.3% | 29.0% |
| Page CUSUM, calibrated on Gaussian data | 5.0% | **70.3%** | **97.3%** | 9.7% | 32.3% |
| Repeated 60-step $z$-test | **70.7%** | **90.0%** | **100%** | **78.0%** | **95.0%** |

*Table 4. Share of paths with a false alarm within one 390-step day, no change, 300 paths per cell (Monte Carlo standard error about 0.8 points at 2% and 1.3 points at 5%). "Price noise" adds i.i.d. noise with standard deviation 0.5 per-step return standard deviations to log prices; "slow SV" is a log-normal stochastic volatility with log-sd 0.15 and half-life 390 steps (E4, Q4a) [SYNTHETIC].*

In the three columns where the null holds, all e-detectors stay below the bound, including the increase detector on raw squared returns, as Section 3.3 predicts. The Gaussian likelihood-ratio and CUSUM detectors, whose thresholds are correct only for Gaussian data, alarm on 70% to 97% of days under jumps or fat tails, and the repeated $z$-test fails everywhere. The last two columns are not cases where the null holds. Price noise makes the conditional mean of $r_t^2$ depend on the previous return (Section 5.4), so the 3.7% shows robustness to a mild violation over one day. Under slow stochastic volatility the conditional variance genuinely drifts, and every method alarms on a quarter or more of days.

### 5.2 How tight is the bound

E9 runs the two-sided monitor with the true baseline ($m_0$ from 20 million draws, $\tau$ floored at three ticks for the tick law) on 4000 paths per law (Table 5, Figure 3).

| Horizon $m$ (steps) | bound $m/A$ | Gaussian | Student $t_3$ | jumps | price ticks |
|---|---|---|---|---|---|
| 98 | 1.3% | 0.00% | 0.18% | 0.12% | 0.05% |
| 390 | 5.0% | 1.1% | 1.4% | 1.3% | 1.4% |
| 1950 | 25% | 12.8% | 14.3% | 13.2% | 13.7% |
| 3900 | 50% | 29.6% | 30.6% | 29.7% | 31.2% |
| 7800 | 100% | 55.9% | 55.9% | 56.4% | 56.9% |

*Table 5. False alarm rate within $m$ steps with the true baseline, two-sided rule, $A = 7800$ (E9, Q9a and Q9b) [SYNTHETIC]. Student $t_3$ has an infinite fourth moment. 95% intervals at one day are about $\pm 0.3$ points.*

![Figure 3. Left: false alarm rate within $m$ steps against the bound $m/A$ for four return laws, true baseline. Right: one-day false alarm rate under GARCH(1,1) clustering, where the null is false, as a function of the band $b$ (E9) [SYNTHETIC].](../results/figures/e09_stress.png)

The bound holds for every law, and the four laws give almost the same rates, which is what a nonparametric guarantee should look like once the statistic is bounded. It is conservative by a factor of about 4 at one day and about 2 at $m = A$. E11 locates the slack with an exact optional stopping identity that writes the ratio of the false alarm rate to $m/A$ as a product of a stopping factor, the share $F_2$ of expected wealth held on paths that crossed $A$, and an overshoot factor $F_3$ (E11, Q11d). For the default two-sided monitor (40,000 paths, rate 1.03%) the stopping factor is close to 1, $F_2 = 0.24$ and $F_3 = 0.80$: the slack is overshoot (a factor of 1.25) times wealth stranded on paths that never alarm (a factor of about 4). Of the stranded wealth, 60% is held by the decrease side and 31% by the three smallest increase bets, which cannot build wealth within a day. The split of the budget between the two sides costs almost nothing (both sides alarm on the same day in 0.00% of days). The bound itself is nearly sharp: for a single increase bet of 0.5, dynamic programming over the whole nonparametric null finds a law with a one-day false alarm rate of 4.73%, 95% of the bound, against 3.52% for Gaussian data. For the default mixture the worst law found reaches 42% of the bound; a certified sharper bound for mixtures is open. Provable gains are possible by splitting the budget unevenly between the sides and by dropping bets that cannot use it.

### 5.3 Detection delay, bet design and the cost of an estimated baseline

**Against the benchmarks.** Table 6 reports median delays after a change 500 steps into monitoring (E4, Q4b).

| Detector | Gaussian | jumps | Student $t_5$ |
|---|---|---|---|
| *Increases, sd ratio 1.25, 1.5, 2* | | | |
| Robust e-detector (default grid) | 87, 35, 17 ($\le$ 4.3%) | 102, 36, 18 ($\le$ 4%) | 122, 49, 22 (3%) |
| Robust e-detector, CUSUM type | 154, 51, 23 ($\le$ 1%) | 171, 52, 23 ($\le$ 0.3%) | 215, 72, 32 ($\le$ 0.3%) |
| E-detector on raw $r_t^2$ | 84, 33, 15 ($\le$ 4.3%) | 246, 60, 23 ($\le$ 5.3%) | 134, 51, 21 ($\le$ 4%) |
| Gaussian likelihood-ratio e-detector | 78, 26, 10 ($\le$ 5.3%) | 61, 20, 9 (**74 to 78%**) | 46, 10, 7 (**96 to 98%**) |
| Page CUSUM | 64, 24, 10 (6 to 10%) | 54, 20, 11 (**74 to 79%**) | 50, 9, 6 (**97 to 99%**) |
| *Decreases, sd ratio 0.5, 0.67, 0.8* | | | |
| Robust e-detector (default grid) | 71, 92, 147 ($\le$ 2%) | 72, 98, 160 (1%) | 68, 96, 168 ($\le$ 3%) |
| Gaussian likelihood-ratio e-detector | 23, 46, 104 (6 to 9%) | 23, 49, 121 (3 to 7%) | 29, 81, 282 ($\le$ 2%) |

*Table 6. Median delay in steps after a change at step 500 of monitoring; in brackets, the share of paths with a false alarm before the change. 300 paths per cell (E4, Q4b) [SYNTHETIC]. The decrease detector on raw $r_t^2$ does not exist (Section 3.3).*

On clean Gaussian data the parametric detectors are faster, and the robust detector pays 10% to 70% extra delay for increases, the penalty that the information ladder of Table 2 predicts. With jumps or fat tails the parametric detectors alarm before the change on 74% to 99% of paths, so their delays have no meaning. Truncation matters for increases under jumps: at sd ratio 1.25 the median delay is 102 steps with truncation and 246 without. The CUSUM-type version of the robust detector is slower, as $C_t \le M_t$ implies. For decreases the Gaussian likelihood-ratio detector was faster and close to its nominal false alarm rate in these designs, because jumps only add evidence against a decrease; this rests on normality, and real one-minute returns contain many exact zeros (5.7% of Bitstamp BTC/USD minutes, E5 [REAL]), which a Gaussian likelihood reads as evidence of low volatility.

**Bet design.** The quarticity estimator can enter the detector through the bet. A plug-in bettor with $\lambda_t = (\hat\mu_t - 1)/(\hat v_t + (\hat\mu_t - 1)^2)$, where $\hat\mu_t$ and $\hat v_t$ are kernel-weighted local mean and variance of $Y$ over the last 60 steps (so that the local second moment is a truncated local quarticity), added with weight one half to the grid, was slower than the grid alone in 14 of 18 cells, by up to 9% (E4, Q4c). E11 explains this negative result: the quadratic approximation behind the bet makes it 1.4 to 4.4 times too small for increases, and it keeps only 49% to 94% of the attainable growth (E11, Q11a). A grid designed by minimising the worst-case ratio to the best single bet over a range of changes (sd ratio 1.2 to 2 up, 0.5 to 0.85 down) and four tail shapes, $\lambda \in \lbrace 0.16, 0.39 \rbrace$ up and $\lambda/\lambda_{\max} \in \lbrace 0.5, 0.85 \rbrace$ down, was at least as fast as the default in all 48 validation cases in its design range, by 0.4% to 10.5% (median 5%), but up to 17% slower for small increases (sd ratio 1.1) under fat tails. A grid tuned on Gaussian data alone was up to 18% faster on Gaussian or jump data and up to 55% slower under Student $t_3$ data (E11, Q11c). The gains from bet design are modest, and tuning to one tail shape is risky.

**The cost of an estimated baseline.** E8 compares the plug-in baseline with the confidence-sequence baseline of Theorem D ($\delta = 0.01$ per side and $A = 9750$, so that $\delta + m/A = 5\%$ per day). Without a change the confidence-sequence monitor had at most 0.27% false alarms per day in nine settings (three laws, calibration of 390, 1440 and 5000 steps, 10,000 paths each), and the confidence bound failed on at most 0.08% of paths against an allowed 1% (E8, Q8a). The plug-in was at 4.5% to 4.9% on the increase side at 390 calibration steps, at the bound but not beyond it within Monte Carlo error at one day; Section 3.4 shows that it does exceed the bound over longer horizons. The price of the proof is delay (Table 7).

| Returns, sd ratio | 390 steps: plug-in / CS | 1440 steps: plug-in / CS | 5000 steps: plug-in / CS |
|---|---|---|---|
| Gaussian, 1.25 | 90 / 645 (76%) | 88 / 222 | 88 / 134 |
| Gaussian, 1.5 | 36 / 103 | 35 / 59 | 35 / 46 |
| Gaussian, 2 | 17 / 31 | 17 / 24 | 17 / 20 |
| Gaussian, 0.8 | 143 / 633 (85%) | 137 / 269 | 139 / 199 |
| Gaussian, 0.67 | 89 / 207 | 87 / 135 | 87 / 114 |
| Gaussian, 0.5 | 64 / 122 | 65 / 92 | 64 / 80 |
| Student $t_5$, 1.25 | 128 / 904 (35%) | 131 / 391 | 130 / 219 |

*Table 7. Median delay in steps after a change at step 500, plug-in against confidence-sequence (CS) baseline, by calibration length; in brackets, the share detected within 2500 steps when below 90%. 2000 paths per cell (E8, Q8b) [SYNTHETIC]. With 390 calibration steps the plug-in detects only 92% to 94% of increases in this sense, because it alarms before the change on 6% to 7.5% of paths.*

The proof is expensive with short calibration and cheap with long calibration: at 390 steps the confidence-sequence monitor needs 1.8 to 7.2 times the plug-in delay, at 1440 steps about 1.4 times for a doubling or halving and 1.5 to 3.0 times for smaller changes, and at 5000 steps 1.2 to 1.7 times. Most of the width is intrinsic rather than a weakness of the bound: from 292 observations, the median upper bound sits 37% above the true mean for Gaussian returns, while the best bound that any constant bet could give for the true law (an oracle) sits 29% above it and a normal approximation without a finite-sample guarantee 19% above it; at 1440 steps and more our bound is within 3.5 points of the oracle (E8, Q8d). The truncated squared return has a standard deviation 1.4 to 1.6 times its mean, so a mean is hard to pin down from a few hundred minutes. The choice of $\delta$ matters little between 0.01 and 0.03 (delays within 5%, E8, Q8g), and a calibration that stops when the relative width of the bounds reaches a target behaves almost like a fixed length chosen for that width (E8, Q8c).

### 5.4 Clustering, forecasts and price noise: when the null is false

**Volatility clustering.** Under GARCH(1,1) returns with constant unconditional variance, the conditional mean of $X_t$ moves above and below $m_0$, so the nulls of Section 3 are false. With the true unconditional baseline and $A = 28{,}800$ for a 1440-step day, the share of days with an alarm is 78% at persistence $\alpha + \beta = 0.95$ (half-life 14 steps) and 100% at 0.99, 0.999 and 0.9999 (E9, Q9c; Figure 3, right). A band $b = 0.5$ absorbs short-lived clustering (3.6% at persistence 0.95), but no reasonable static band absorbs persistent clustering: even $b = 1$ on the increase side leaves 12% to 20% of days with an alarm at persistence 0.99 or more. For slow log-normal stochastic volatility, a band of 0.5 brought the daily rates to 1.7% (increase) and 0.0% (decrease) while a doubling or halving of volatility was still detected on 97% and 100% of paths (E4, E4c); the band is a design choice, not a guarantee. These alarms are not errors of the proof: persistent GARCH clusters are changes in conditional volatility, and the monitor reports them.

**Forecast-relative monitoring.** Monitoring $r_t/\hat\sigma_t$ with a GARCH(1,1) forecast fitted on a separate path and then frozen restores the one-day false alarm rate to 3.0% to 6.2% at all three persistences (600 paths, plug-in baseline from one day; the oracle with the true conditional volatility gives 4.7% to 5.3%), while an EWMA forecast is only as valid as its match to the clustering (half-life 360 steps: 72% to 100% at band 0) (E10, Q10a). The price is the signal cap of Proposition 10. At persistence 0.99, a doubling of the long-run variance was detected within one day on 34% to 43% of paths and within five days on 93% to 94%, with median delays of 1523 to 1791 steps; without a break the same monitors alarmed within one day on 4% to 6% and within five days on 38% to 44% of paths, the latter far above the 25% that Theorem B gives with the true baseline (E10, Q10b and Q10b2). A fourfold increase of the long-run variance was detected within a day on 89% to 90% of paths (median 612 to 712 steps). For the oracle, whose standardised returns are exactly i.i.d., the excess over five days is the plug-in effect of Section 3.4, not a property of the forecast.

**Price-level microstructure noise.** With i.i.d. noise on log prices, $\mathbb{E}[r_t^2 \mid \mathcal{F}_{t-1}]$ depends on the last noise value, so squared returns are positively autocorrelated and the null fails without any change in volatility. At one-second sampling (23,400 steps per day, $A$ set to give 5% per day) the one-day false alarm rate of the increase side on raw returns rose from 5.0% without noise to 8.5%, 18.0%, 31.0% and 41.5% at noise-to-signal ratios 0.5, 1, 2 and 4 (200 days per cell); overlapping pre-averaged returns alarmed on every day, because consecutive windows share most of their prices. Pre-averaging over non-overlapping blocks of prices, skipping the one return between blocks, gives each block fresh noise values, so the null holds for the block statistic whenever the spot variance and the noise law are unchanged, and Theorems A and B apply to it for every block length $k \ge 2$; the block monitors were below the bound or within Monte Carlo error of it in all 58 cells tested, and with $k = 10$ were about 8 times slower than raw returns without noise and faster than them from a noise ratio of about 2 (ALT1) [SYNTHETIC]. One-minute crypto returns were not tested for noise of this kind; the effect is largest at the highest frequencies.

### 5.5 Experiments in progress

[PENDING E13: benchmark against published change-point methods, the offline test of Bibinger, Jirak and Vetter (2017) and an online Bayesian change-point detector, under the designs of Table 4 and Table 6.]

[PENDING E14: false alarm and delay results for a running, data-driven truncation level with a guarantee, and for a restart after an alarm without a blind recalibration period.]

## 6 Real data

All results in this section are **[REAL]**. Real data cannot test the false alarm guarantee, which concerns periods in which volatility does not change, because such periods cannot be identified. What real data can show is how often the monitor alarms, whether its alarms are followed by the volatility they announce, and whether they carry information beyond a forecast. In all runs the settings were fixed before the data were examined: they were frozen after one earlier run on Bitstamp BTC/USD for July to December 2025 (E5), or pre-registered (E7). Returns are divided by the predictable intraday profile of Section 2, the truncation is $c = 3$ scaled MADs of the calibration window, the baseline is the plug-in mean over a calibration window of one day (two sessions for the Nifty 50), missing returns are skipped, and the monitor recalibrates after every alarm. The total alarm level is $A = 28{,}800$ for a 1440-minute crypto day (each side at $2A$, Proposition 1) and $A = 7500$ for a 375-bar Nifty session, so that with a known baseline the probability of any false alarm is at most 5% per day or session. Because the baseline is a plug-in, this bound is approximate (Section 3.4). The validation measure is a robust volatility, $1.4826$ times the median of $\lvert u \rvert$ over the next day or session, where $u$ is the deseasonalised return; it is not a quantity a trader can trade.

### 6.1 Five cryptocurrencies

**Data.** Binance spot one-minute bars for BTC, ETH, XRP, SOL and DOGE against USDT, the coins of Bapat (2026), from a public mirror of the exchange's files with no missing minutes, for 1 July to 31 December 2025 (264,959 returns per coin) and, out of sample, 21 January 2024 to 30 June 2025 (527 days); Bitstamp BTC/USD for a cross-venue check (E6).

**Alarms are followed by the announced volatility.** With band $b = 0.25$ the monitor raised 0.66 to 0.71 alarms per day for every coin in both periods. In both periods, 75% to 84% of increase alarms and 83% to 93% of decrease alarms were followed by a day whose robust volatility was on the flagged side of the calibration baseline, against about 50% at random times whose baseline window ended as long before as the alarm's (E6, Q6d). Out of sample these shares fell by only 1 to 4 points. Within 60 minutes of a Binance BTC alarm, Bitstamp had a same-direction alarm 72% of the time, and 73% of Bitstamp alarms had one on Binance, against about 3% for randomly shifted alarm times; across coins about a quarter of alarms coincided with a BTC alarm.

**Alarms carry information beyond a forecast, but a tuned ratio rule does as well.** At each alarm we forecast the robust volatility of the next 1440 minutes with an EWMA (half-life one day) and with a HAR regression (Corsi 2009) refitted daily, using data up to the alarm only. A rate-matched competitor alarms when the robust volatility of the last $S$ minutes ($S = 60$ or 240) divided by that of the preceding day crosses thresholds that were chosen on the same data to give the same numbers of up and down alarms as the e-detector, with the same one-day dead time after each alarm. This hindsight favours the competitor. Table 8 pools the five coins.

| Rule | Period | Up alarms: next day above forecast | Down alarms: next day below forecast |
|---|---|---|---|
| Random times (EWMA / HAR) | Jul to Dec 2025 | 49% / 51% | 51% / 49% |
| Robust e-detector (EWMA) | Jul to Dec 2025 | 79% (74 to 82), $n = 367$ | 81% (75 to 85), $n = 256$ |
| Robust e-detector (HAR) | Jul to Dec 2025 | 81% (76 to 84) | 78% (73 to 83) |
| Ratio rule, 60 / 240 min (EWMA) | Jul to Dec 2025 | 74% / 74% | 78% / 82% |
| Random times (EWMA / HAR) | out of sample | 50% / 48% | 50% / 52% |
| Robust e-detector (EWMA) | out of sample | 77% (74 to 80), $n = 951$ | 71% (68 to 74), $n = 870$ |
| Robust e-detector (HAR) | out of sample | 76% (73 to 79) | 74% (71 to 77) |
| Ratio rule, 60 / 240 min (EWMA) | out of sample | 77% / 76% | 67% / 71% |

*Table 8. Share of alarms after which the next day's robust volatility was on the flagged side of a forecast made at the alarm, pooled over five coins, band 0.25, Wilson 95% intervals in brackets (E6, Q6e) [REAL]. Alarms across coins are correlated, so pooled intervals are optimistic.*

After an increase alarm the next day's volatility exceeded both forecasts 76% to 81% of the time, against 48% to 51% at random times, and after a decrease alarm it fell short 71% to 81% of the time. The rate-matched ratio rule was within 0 to 5 points at band 0.25, and none of these differences was significant. Only at the strictest band (0.5), out of sample and against the EWMA forecast, was the e-detector better, by about 4 points for increases ($p = 0.03$ to 0.06) and 10 points for decreases against the 60-minute rule ($p < 0.001$); against the HAR forecast the differences were smaller and not significant. The increase alarms of the Gaussian likelihood-ratio monitor were 7 to 14 points less informative than those of the robust monitor. On real data, then, the e-detector does not buy precision over a simple rule tuned with hindsight. What it adds is a false alarm guarantee with settings fixed in advance, which the ratio rule does not have.

### 6.2 The Nifty 50 index and India VIX

**Data.** Nifty 50 one-minute bars from April 2017 to October 2026, sessions of 375 bars, overnight returns excluded. Stretches of five or more minutes with a frozen index value (exchange outages, the two market-wide halts of March 2020, and the closing auction minutes from 3 August 2026) were set to missing, leaving 2,361 sessions and 882,001 returns. The intraday profile is session-aware, so that the close is not smoothed into the open. Pre-registered configuration: calibration 750 bars, band 0.25; development period 2017 to 2021, held-out period 2022 to October 2026, with nothing changed after the held-out period was seen (E7).

**Results.** The monitor raised 0.26 alarms per session in both periods. In the held-out period, 78% of increase alarms were followed by a session more volatile than the baseline (47% at random times) and 76% of decrease alarms by a calmer one (53%); the Gaussian likelihood-ratio monitor's increase alarms were followed by higher volatility only 61% of the time. India VIX agreed with the alarms on the same day (in the development years VIX closed higher on 68% of days with an increase alarm, against 47% of all days) but did not move in the alarm's direction on the next day more often than on other days (45% to 48% against 46% to 47%): the agreement is concurrent, not a lead. Among event days, 9 of 11 Union Budget sessions had an alarm, against 25% of ordinary sessions (Fisher $p = 0.0001$), and all 7 budget sessions at which the monitor was already watching at the open; Reserve Bank of India policy days were not different (32% against 25%, $p = 0.36$), consistent with their volatility being close to the recent norm (median 1.06 times the previous 20 sessions). Monthly derivative expiry days had more alarms (38%, $p = 0.002$), which may be an artefact of a profile that pools all weekdays. In March 2020 volatility was 2.7 times that of the previous 20 sessions, yet only 6 of 21 sessions had an alarm: after each alarm the monitor recalibrates on the new, already high level. It answers "has volatility changed from the recent baseline?", not "is volatility high?".

### 6.3 FOMC statements and CPI releases

**Design.** Binance BTC, ETH and SOL from 21 January 2024 to 31 December 2025; the 16 FOMC rate decisions (statement at 14:00 US Eastern) and 22 US CPI releases (08:30 US Eastern) of 2024 and 2025 inside the monitored period, with every date checked against a published schedule. Control windows are all days at least two days away from any event day, with the same weekday, at the same Eastern clock time (134 for FOMC, 299 for CPI). The outcome is an increase alarm within six hours of the anchor, and the test is an exact stratified permutation test. The primary family, fixed before the first run, is FOMC and CPI on BTC in the streaming design (the monitor as a user would have run it), with a Holm correction; ETH and SOL are replications, and a design in which a fresh monitor is calibrated on the day before every anchor was added after the streaming results were seen and is secondary (E12).

**Results.** In the streaming design, 9 of 16 FOMC statements on BTC were followed by an increase alarm within six hours, against 14% of control windows (Holm-adjusted $p = 0.0002$); the CPI test was not significant on BTC (5 of 22 against 11%, $p = 0.11$) but replicated on ETH ($p = 0.003$) and SOL ($p = 0.00004$). The streaming monitor was recalibrating after its own alarms about 70% of the time and was awake at only 9 of the 16 FOMC statements on BTC; at those it alarmed every time (and at 8 of 8 on ETH and 9 of 9 on SOL). With a fresh calibration, 14 to 15 of the 16 FOMC statements raised an alarm within six hours, against 29% to 31% of controls, and 59% to 73% of CPI releases, against 25% to 27%. The alarms came after the release (median 23 to 34 minutes, never before 14 minutes on BTC and SOL and 10 minutes on ETH), and the hours before a release had fewer alarms than ordinary hours. The Gaussian likelihood-ratio monitor with a fresh calibration also alarmed at almost every release, but with a median delay of 0 to 2 minutes, on a single return of 8 to 9 robust scale units (median), and in 64% to 77% of control windows: one large minute is enough to trigger it. The truncated detector cannot be moved by one minute (one clipped observation multiplies its statistic by a bounded factor), so its alarms mark a short volatility regime rather than a jump.

### 6.4 What the real data say about the method

Three lessons qualify the simulation results. First, alarm rates on real data are far above the 5% budget: 0.66 to 0.71 per day on crypto at band 0.25, and 0.14 to 0.57 per session on the Nifty 50 across configurations. This does not contradict the guarantee, which concerns periods in which the null holds; one-minute volatility relative to a one-day baseline simply moves often. Dividing returns by an EWMA forecast with a half-life of 6 or 24 hours changed the crypto alarm rate by between +2% and $-17\%$ (E10, Q10c). Second, a monitor that recalibrates for a full day after every alarm is blind most of the time when alarms come every 1.4 days (about 70% in E12, and 85% with band 0 on Bitstamp BTC/USD, E8). A restart without a blind period is studied in E14 (Section 5.5). Third, on Bitstamp BTC/USD the confidence-sequence baseline of Theorem D reduced the number of alarms by 16% to 20% at every band, almost all of it in decrease alarms, and behaved like a data-driven band of about 0.16 to 0.18 (E8, Q8e). Its value on real data is the proof, not better alarms. Finally, our crypto and Nifty data come from public mirrors of exchange data rather than from the exchanges, and Bapat (2026) does not name the source of its crypto data, so the comparison with its tables in Section 7 uses a data set that is close to, but not the same as, the published one.

## 7 What a reproduction found about the estimator

Before building the monitor we reimplemented the estimators of Bapat (2026) as printed and compared them with the published simulation and real-data tables. We report four findings, as questions for the author rather than corrections. Where we say what the code behind the tables probably does, this is an inference from published numbers, not a fact. None of the findings touches the central idea, truncation at a robust scale followed by kernel smoothing, and the monitor of this paper uses the scaled MAD and the squared statistic, so none of them affects Sections 3 to 6. We do not claim to have improved the estimator.

**Question 1: are the kernel weights normalised?** The local estimator $\hat Q_t = \tfrac13 \sum_{i=t-h}^{t} K((i-t)/h)\, r_i^4$ (eqs. 6 and 12) uses raw kernel weights and no factor for the step length $\Delta$. Under locally constant volatility $\mathbb{E}[\hat Q_t] = W(h)\,\sigma_t^4\Delta^2$ with $W(h) = \sum_{k=0}^{h} K(k/h)$, about $0.341h + 0.32$ for the Gaussian kernel ($W(30) = 10.6$), so the estimate grows with the bandwidth. The real-data tables give a clean check: if the kernel estimate is the average of $\hat Q_t$ over the sample and RQ is $\tfrac13\sum r_i^4$, their ratio is close to $W(30)/(n - 30) \approx 3.99 \times 10^{-5}$ whatever the data. The published ratios are 4.00 to $4.03 \times 10^{-5}$ for all five coins [PUBLISHED], and our Binance data give $3.986 \times 10^{-5}$ for all five [REAL] (E1, E6). Once both kernels are put on the same scale, the comparison is fair and the robust idea holds up: with the scaled MAD the normalised robust kernel has a smaller RMSE than the normalised plain kernel in all 45 cells of the simulation design, and with the raw MAD it wins whenever jumps are moderate or large (E1) [SYNTHETIC].

**Question 2: is the MAD scaled by 1.4826?** Because the MAD shrinks at the same rate as the returns, $\tau/(\sigma\sqrt\Delta)$ does not change with $n$, the share of clipped returns stays fixed, and the truncated estimator converges to $\kappa_4(a)$ times integrated quarticity rather than to integrated quarticity: $\kappa_4 = 0.705$ with the raw MAD ($c = 3$ clips at 2.02 standard deviations) and 0.964 with the scaled MAD. In simulation the ratio to the truth was 0.705 at $n = 390$, 1440 and 23,400 (E2) [SYNTHETIC]. Thresholds of the form $c\Delta^{\varpi}$ with $\varpi < 1/2$ (Mancini 2009) avoid this, as does division by $\kappa_4$. Two independent checks point to the raw MAD in the published tables: the simulated robust RMSE matches the published one with a median ratio of 0.90 for the raw MAD against 1.29 for the scaled MAD (E1), and on four coins, after adjusting for differences in quartiles between the data sets, the robust kernel matches the published values with ratios of 0.81 to 0.96 for the raw MAD against 1.86 to 2.48 for the scaled MAD (E6) [PUBLISHED against SYNTHETIC and REAL]. For a monitor this constant is a matter of definition, since the null concerns the truncated variance itself; for an estimator of integrated quarticity it has to be divided out.

**Question 3: is the averaged estimator meant to equal truncated realised quarticity?** Exchanging the two sums in eq. (16), every interior return collects each kernel weight exactly once, so apart from the first and last $h$ returns the averaged estimator equals $W(h)/(3(n-h))\sum_i \tilde r_i^{\,4}$, a constant times the realised quarticity of the truncated returns. In simulation the normalised averaged robust kernel and plain truncated realised quarticity had correlation 0.98 to 0.99 and mean ratio within 0.4% of 1 in seven scenarios (E2) [SYNTHETIC]. For integrated quarticity, therefore, the kernel and the bandwidth mostly change a constant and the edges. The kernel's real contribution is the local estimate $\hat Q_t$, which is what an online procedure needs.

**Question 4: is one scale per day or per sample enough when volatility moves within the day?** Under stochastic volatility with an intraday pattern, a single daily MAD left a bias of about $-14\%$ at $n = 23{,}400$ even after dividing by $\kappa_4$, and its RMSE fell like $n^{-0.1}$, while estimators truncating against a local scale (a few dozen neighbouring returns) converged at the usual rate (fitted slopes $-0.50$ to $-0.68$) (ALT3) [SYNTHETIC]. On one-minute crypto data, 10 to 15 times more returns are clipped within a day than for Gaussian returns, and with one MAD over six months, as in the published crypto table, the robust kernel keeps under 2% of the fourth moment (E6) [REAL]. Large one-minute returns come in bursts (44% to 65% have a large neighbour), and on BTC and DOGE the scale-truncated estimators come out at 0.33 to 0.44 times the median realised quarticity of Andersen, Dobrev and Schaumburg (2012), local plug-ins at 0.66 to 0.73 (ALT3) [REAL]. Was the MAD in the crypto application computed per day, over the whole sample, or locally?

The accompanying documentation lists smaller points in the same spirit (the weight order in eq. 15, the noise model, the constants of the bipower and tripower quarticities, how the RMSE of the real-data table was computed, and the DOGE row, whose zero quartiles imply a zero MAD and hence a zero robust kernel, while a positive value is printed).

## 8 Discussion, limitations and open problems

**What is new and what is not.** The probabilistic tools (Ville's inequality, test supermartingales, Shiryaev-Roberts e-detectors, empirical Bernstein confidence sequences, the Bahadur-Savage phenomenon, Lorden's and Lai's lower bounds) are known, and Theorems A and B are instances of known results. What this paper adds is their combination for jump-robust online volatility monitoring: the choice of the clipped squared return with a predictable robust truncation level as the monitored statistic, the conditional-mean nulls and their volatility meaning, the explanation of why a decrease detector needs a bounded statistic and what speed limit truncation then imposes, the confidence-sequence baseline with its proof, the decomposition of the delay cost of robustness, the signal cap under forecast-relative monitoring, and the evidence of Sections 5 and 6. Sections 3.3 and 4 are where we would most value a check by a specialist in sequential analysis.

**Limitations.**

* *The null is a statement about conditional means.* Persistent volatility clustering, price-level microstructure noise, estimated baselines and errors in the intraday profile all move the conditional mean of $X_t$ without a structural change. Each has a remedy (forecast-relative monitoring, non-overlapping pre-averaged blocks, the confidence-sequence baseline), and each remedy costs speed. The effect of an estimated intraday profile was not isolated in any experiment.
* *Tail shape.* By Lemma 1, a change in the shape of the return distribution at constant volatility, or a leverage effect, can move the truncated variance and raise an alarm that is not a volatility change.
* *Estimated baselines.* The plug-in baseline that the real-data runs use has no proof and exceeds the bound over horizons of several days (Section 3.4). The confidence-sequence baseline has a proof but needs long calibration to be fast (Table 7).
* *Conservatism.* With the default grid the bound is about four times larger than the observed rate at one day (Section 5.2).
* *Restarts and blindness.* Corollaries C and F bound the expected number of false alarms when the monitor recalibrates after every alarm; they are not a single anytime-valid statement for a whole series, and on real data the monitor is blind most of the time (Section 6.4).
* *Several assets.* Five coins monitored at 5% per day each give, by the union bound, up to 25% per day for at least one false alarm in the group. [PENDING E15: panel monitoring of many assets with false discovery control.]
* *Delay theory.* The delay results assume i.i.d. post-change data and a known baseline; the second-order formula is a heuristic whose regularity conditions we did not check; Theorem 9 is asymptotic and the asymptotic regime arrives slowly (Section 4.3).
* *Data.* One-minute bars only, from public mirrors of exchange data; no tick data; real-data validation has no ground truth and uses indirect checks (later volatility, forecasts, India VIX, scheduled events).

**Open problems.**

1. *A certified sharper threshold for mixtures.* For single bets the bound $m/A$ is nearly attained by an adversarial null law, but for the default mixture the worst law found reaches only 42% of it (Section 5.2). A sharper bound that is still provable would translate directly into speed.
2. *A running truncation level and a restart without a blind period, with a proof that costs little.* Real data want $\tau_t$ to follow the current scale. A band that holds for all truncation levels at once, rank-based (conformal) bets that need no baseline, and restarts that learn the baseline while monitoring are candidates (E14, Section 5.5).
3. *A narrower baseline.* A betting confidence sequence could recover at most the gap to the oracle bound of Section 5.3, 7.6 to 10.6 points at 390 calibration steps.
4. *Theory under clustering.* A proof of the cap for the mean of the ratio $\mathbb{E}[r_t^2/s_t]$ that the monitor actually sees, delay bounds for forecast-relative monitoring with estimated forecast parameters, and a principled choice of band for clustered data.
5. *What should the monitor watch?* The same detectors applied to the clipped fourth power were faster for increases and slower for decreases (ALT1). Whether to monitor the truncated variance, the truncated quarticity or both, and how to combine them with one false alarm budget, is a question on which the author of the quarticity estimator is better placed than we are.
6. *Estimation.* A central limit theorem for truncation against a local scale, which the evidence of Section 7 (Question 4) suggests is the right object for one-minute data, and a comparison of classical purely sequential fixed-accuracy estimation of volatility, with a robust quarticity estimate as the variance plug-in, against anytime-valid confidence sequences.

## References

Andersen, T. G. and Bollerslev, T. (1997). Intraday periodicity and volatility persistence in financial markets. *Journal of Empirical Finance* 4(2-3), 115-158. doi:10.1016/S0927-5398(97)00004-2

Andersen, T. G., Dobrev, D. and Schaumburg, E. (2012). Jump-robust volatility estimation using nearest neighbor truncation. *Journal of Econometrics* 169(1), 75-93. doi:10.1016/j.jeconom.2012.01.011

Andersen, T. G., Todorov, V. and Zhou, B. (2025). Real-time detection of local no-arbitrage violations. *Quantitative Economics* 16(2), 459-495. doi:10.3982/QE2585

Bahadur, R. R. and Savage, L. J. (1956). The nonexistence of certain statistical procedures in nonparametric problems. *Annals of Mathematical Statistics* 27(4), 1115-1122.

Bapat, S. R. (2026). A jump and noise robust kernel based estimator of quarticity for high frequency data. *International Review of Economics and Finance* 109, 105482. doi:10.1016/j.iref.2026.105482

Barndorff-Nielsen, O. E. and Shephard, N. (2004). Power and bipower variation with stochastic volatility and jumps. *Journal of Financial Econometrics* 2(1), 1-37. doi:10.1093/jjfinec/nbh001

Bibinger, M., Jirak, M. and Vetter, M. (2017). Nonparametric change-point analysis of volatility. *Annals of Statistics* 45, 1542-1578. doi:10.1214/16-AOS1499

Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. *Journal of Econometrics* 31(3), 307-327. doi:10.1016/0304-4076(86)90063-1

Boudt, K., Croux, C. and Laurent, S. (2011). Robust estimation of intraweek periodicity in volatility and jump detection. *Journal of Empirical Finance* 18(2), 353-367. doi:10.1016/j.jempfin.2010.11.005

Corsi, F. (2009). A simple approximate long-memory model of realized volatility. *Journal of Financial Econometrics* 7(2), 174-196. doi:10.1093/jjfinec/nbp001

Diebold, F. X. (1986). Modeling the persistence of conditional variances: a comment. *Econometric Reviews* 5(1), 51-56.

Dufresne, D. (1990). The distribution of a perpetuity, with applications to risk theory and pension funding. *Scandinavian Actuarial Journal* 1990(1), 39-79.

Fan, X., Grama, I. and Liu, Q. (2015). Exponential inequalities for martingales with applications. *Electronic Journal of Probability* 20, paper 1, 1-22. doi:10.1214/EJP.v20-3496

Grünwald, P., de Heide, R. and Koolen, W. (2024). Safe testing. *Journal of the Royal Statistical Society Series B* 86(5), 1091-1128. doi:10.1093/jrsssb/qkae011

Honda, J. and Takemura, A. (2010). An asymptotically optimal bandit algorithm for bounded support models. *Proceedings of the 23rd Conference on Learning Theory (COLT 2010)*.

Horváth, L., Trapani, L. and Wang, S. (2024). Sequential monitoring for explosive volatility regimes. arXiv:2404.17885.

Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2021). Time-uniform, nonparametric, nonasymptotic confidence sequences. *Annals of Statistics* 49(2), 1055-1080. doi:10.1214/20-AOS1991

Inclán, C. and Tiao, G. C. (1994). Use of cumulative sums of squares for retrospective detection of changes of variance. *Journal of the American Statistical Association* 89(427), 913-923.

Lai, T. L. (1998). Information bounds and quick detection of parameter changes in stochastic systems. *IEEE Transactions on Information Theory* 44(7), 2917-2929.

Lamoureux, C. G. and Lastrapes, W. D. (1990). Persistence in variance, structural change, and the GARCH model. *Journal of Business and Economic Statistics* 8(2), 225-234.

Larsson, M., Ramdas, A. and Ruf, J. (2024). The numeraire e-variable and reverse information projection. arXiv:2402.18810.

Lorden, G. (1970). On excess over the boundary. *Annals of Mathematical Statistics* 41(2), 520-527.

Lorden, G. (1971). Procedures for reacting to a change in distribution. *Annals of Mathematical Statistics* 42(6), 1897-1908.

Mancini, C. (2009). Non-parametric threshold estimation for models with stochastic diffusion coefficient and jumps. *Scandinavian Journal of Statistics* 36, 270-296.

Mikosch, T. and Starica, C. (2004). Nonstationarities in financial time series, the long-range dependence, and the IGARCH effects. *Review of Economics and Statistics* 86(1), 378-390. doi:10.1162/003465304323023886

Page, E. S. (1954). Continuous inspection schemes. *Biometrika* 41(1-2), 100-115. doi:10.1093/biomet/41.1-2.100

Pollak, M. (1985). Optimal detection of a change in distribution. *Annals of Statistics* 13(1), 206-227. doi:10.1214/aos/1176346587

Pollak, M. (1987). Average run lengths of an optimal method of detecting a change in distribution. *Annals of Statistics* 15(2), 749-779. doi:10.1214/aos/1176350373

Ram, A. and Ramdas, A. (2026). Asymptotically optimal sequential change detection for bounded means. arXiv:2602.05272.

Ramdas, A., Grünwald, P., Vovk, V. and Shafer, G. (2023). Game-theoretic statistics and safe anytime-valid inference. *Statistical Science* 38(4), 576-601. doi:10.1214/23-STS894

Roberts, S. W. (1966). A comparison of some control chart procedures. *Technometrics* 8(3), 411-430.

Shafer, G. (2021). Testing by betting: a strategy for statistical and scientific communication. *Journal of the Royal Statistical Society Series A* 184(2), 407-431. doi:10.1111/rssa.12647

Shin, J., Ramdas, A. and Rinaldo, A. (2024). E-detectors: a nonparametric framework for sequential change detection. *New England Journal of Statistics in Data Science* 2(2), 229-260. doi:10.51387/23-NEJSDS51

Shiryaev, A. N. (1963). On optimum methods in quickest detection problems. *Theory of Probability and Its Applications* 8(1), 22-46. doi:10.1137/1108002

Ville, J. (1939). *Étude critique de la notion de collectif*. Thèse, Paris.

Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded random variables by betting. *Journal of the Royal Statistical Society Series B* 86(1), 1-27. doi:10.1093/jrsssb/qkad009

Woodroofe, M. (1982). *Nonlinear Renewal Theory in Sequential Analysis*. CBMS-NSF Regional Conference Series in Applied Mathematics 39. SIAM, Philadelphia.
