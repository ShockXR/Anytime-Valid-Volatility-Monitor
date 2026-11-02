"""Strong published competitors for the volatility monitor (experiment E13, tracker B67).

Two families of methods are implemented here so that the robust e-detector can be
benchmarked against them on equal terms.

1. Adjacent-window tests for a jump in volatility, in the spirit of Bibinger, Jirak
   and Vetter (2017), "Nonparametric change-point analysis of volatility", Annals of
   Statistics 45(4), 1542 to 1578, doi:10.1214/16-AOS1499, arXiv:1502.00043.
   Local realised variances over two adjacent windows of ``k`` returns are compared,
   the comparison is standardised so that it is asymptotically standard normal when
   volatility does not jump, and the maximum over the sample is compared with an
   extreme value (Gumbel-type) critical value. A jump-robust variant drops returns
   above a threshold (truncated realised variance, Mancini 2009) before comparing,
   and the location of the maximum estimates the change point.

   IMPORTANT, read before citing equation numbers. The text of the paper could not be
   retrieved in the environment where this module was written (arxiv.org, its
   ar5iv mirror and the Humboldt University repository edoc.hu-berlin.de were blocked
   by the network policy). What was verified through web search: the paper builds its test
   from rescaled differences of local spot volatility estimates on adjacent
   windows of ``k_n`` observations, its main statistic uses non-overlapping blocks and
   a second statistic uses overlapping windows (reported to be more powerful), the
   null limit is a Gumbel-type extreme value law, and its simulations use
   ``k_n = 125``. The exact normalising constants and the published equation numbers
   were NOT checked. The constants below are derived here from standard extreme
   value theory for maxima of (weakly dependent) Gaussian variables, and for the
   overlapping version from the moving-sum (MOSUM) limit used by the CRAN package
   ``mosum`` (function ``mosum.criticalValue``, Eichinger and Kirch 2018, Bernoulli
   24(1)). Because convergence to these limits is slow, every experiment also
   reports critical values calibrated by simulation. Treat these functions as a
   faithful implementation of the published IDEA, not as a line-by-line
   reproduction, until the equations are checked against the paper.

2. Bayesian online change-point detection (Adams and MacKay 2007, arXiv:0710.3742,
   Algorithm 1) for a zero-mean Gaussian with unknown variance and an inverse gamma
   prior on the variance (conjugate), a constant hazard, and a capped run-length
   posterior. The alarm rule is "posterior probability that the current run is
   shorter than ``short`` steps is at least a threshold". There is no frequentist
   false alarm guarantee; the threshold has to be calibrated by simulation.

Units: returns are per-step returns (any scale). Arrays can be one path ``(n,)`` or
many paths ``(reps, n)``. Missing returns (NaN) are skipped where stated.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy import stats
from scipy.special import gammaln, logsumexp

# ----------------------------------------------------------------------------
# Extreme value critical values
# ----------------------------------------------------------------------------


def gumbel_crit_max_abs_normal(n_terms: int, alpha: float = 0.05) -> float:
    """Asymptotic critical value for ``max_{j <= N} |Z_j|`` with ``Z_j`` standard normal.

    With ``a = sqrt(2 log 2N)`` and ``b = a - (log log 2N + log 4 pi) / (2 a)``, the
    classical result for maxima of independent normals (applied to ``2N`` one-sided
    tails, since ``P(|Z| > u) = 2 (1 - Phi(u))``) is ``P(a (max |Z_j| - b) <= x) -> exp(-e^{-x})``.
    The same limit holds for stationary Gaussian sequences whose correlation dies out
    (Berman's condition), which covers the adjacent-block statistics below (lag-one
    correlation -1/2, independent beyond). Returns ``b - log(-log(1 - alpha)) / a``.
    """
    if n_terms < 2:
        raise ValueError("need at least two terms")
    N2 = 2.0 * n_terms
    a = np.sqrt(2.0 * np.log(N2))
    b = a - (np.log(np.log(N2)) + np.log(4.0 * np.pi)) / (2.0 * a)
    return float(b - np.log(-np.log(1.0 - alpha)) / a)


def mosum_crit(n: int, k: int, alpha: float = 0.05) -> float:
    """Asymptotic critical value for the maximum of a standardised two-window moving sum.

    For windows of ``k`` on each side and ``x = n / k``: ``a = sqrt(2 log x)``,
    ``b = 2 log x + log log x / 2 + log(3/2) - log(pi) / 2`` and
    ``P(a max_i |D_i| - b <= y) -> exp(-2 e^{-y})``. This is the symmetric case of
    ``mosum.asymptoticA`` / ``mosum.asymptoticB`` / ``mosum.criticalValue`` in the CRAN
    package ``mosum`` (source verified at
    https://raw.githubusercontent.com/cran/mosum/master/R/mosum_test.R).
    """
    x = n / k
    if x <= np.e:
        raise ValueError("n / k must exceed e")
    a = np.sqrt(2.0 * np.log(x))
    b = 2.0 * np.log(x) + 0.5 * np.log(np.log(x)) + np.log(1.5) - 0.5 * np.log(np.pi)
    return float((b - np.log(np.log(1.0 / np.sqrt(1.0 - alpha)))) / a)


def truncated_square_dispersion(c: float | None) -> float:
    """``Var(Z^2 1{|Z| <= c}) / E[Z^2 1{|Z| <= c}]^2`` for standard normal Z (2 without truncation).

    Uses ``E[Z^2 1] = (2 Phi(c) - 1) - 2 c phi(c)`` and ``E[Z^4 1] = 3 E[Z^2 1] - 2 c^3 phi(c)``
    (``docs/math.md``, section 3).
    """
    if c is None or not np.isfinite(c):
        return 2.0
    phi, Phi = stats.norm.pdf(c), stats.norm.cdf(c)
    m2 = (2 * Phi - 1) - 2 * c * phi
    m4 = 3 * m2 - 2 * c**3 * phi
    return float(m4 / m2**2 - 1.0)


# ----------------------------------------------------------------------------
# Adjacent-window volatility jump statistics (Bibinger, Jirak and Vetter type)
# ----------------------------------------------------------------------------

BPV_CONST = np.pi / 2.0  # E|Z| ^ -2 for standard normal Z: bipower variation (pi/2) sum |r_i||r_{i-1}| estimates IV


def _pair_z(S1L, S1R, S2L, S2R, NL, NR, rho, norm):
    """Standardised difference of the mean (truncated) squared return on the right and left windows."""
    with np.errstate(invalid="ignore", divide="ignore"):
        mL, mR = S1L / NL, S1R / NR
        if norm == "gaussian":
            m = (S1L + S1R) / (NL + NR)
            z = (mR - mL) / (m * np.sqrt(rho * (1.0 / NL + 1.0 / NR)))
        elif norm == "empirical":
            vL = np.maximum(S2L / NL - mL**2, 0.0)
            vR = np.maximum(S2R / NR - mR**2, 0.0)
            z = (mR - mL) / np.sqrt(vL / NL + vR / NR)
        else:
            raise ValueError("norm must be 'gaussian' or 'empirical'")
    return z


@dataclass
class BJVResult:
    stat: np.ndarray  # max_j |Z_j| per path
    crit: float  # critical value used
    reject: np.ndarray  # stat > crit
    cp_index: np.ndarray  # estimated change point: index of the first return after the change
    direction: np.ndarray  # +1 increase, -1 decrease (sign of Z at the maximum)
    z: np.ndarray  # all standardised comparisons, shape (reps, m)
    positions: np.ndarray  # change-point index each comparison refers to


def bjv_block_statistics(r, k: int, trunc_c: float | None = 4.0, norm: str = "gaussian"):
    """Standardised comparisons of adjacent NON-OVERLAPPING blocks of ``k`` returns.

    Block j holds returns ``j k, ..., (j + 1) k - 1`` (a tail shorter than k is dropped).
    For each pair of neighbouring blocks (j, j + 1) the statistic is

    ``Z_j = (RV_{j+1} - RV_j) / ( RV_pool * sqrt(rho (1/N_j + 1/N_{j+1})) )``

    where ``RV`` is the mean squared return in a block, ``RV_pool`` the mean over both
    blocks, ``N`` the number of finite returns and ``rho = Var(r^2) / E[r^2]^2`` (2 for
    normal returns). With locally normal returns and no volatility jump, ``Z_j`` is
    asymptotically standard normal, consecutive ``Z_j`` have correlation -1/2 and
    non-neighbours are independent. ``norm="empirical"`` uses the sample variance of
    the squared returns in each block instead of ``rho RV^2`` (robust to fat-tailed
    innovations at the price of a noisier denominator).

    Jump-robust variant (``trunc_c`` not None): in each pair, returns with
    ``|r| > trunc_c * sqrt(BPV)`` are dropped, where BPV is the larger of the two
    blocks' bipower variations ``(pi/2) mean |r_i| |r_{i-1}|`` (a jump-robust local
    variance). ``rho`` is then the truncated normal value. Returns ``(Z, positions)``
    where ``positions[j] = (j + 1) k`` is the first index of the right block.
    """
    r = np.atleast_2d(np.asarray(r, float))
    reps, n = r.shape
    m = n // k
    if m < 3:
        raise ValueError("need at least three blocks")
    B = r[:, : m * k].reshape(reps, m, k)
    fin = np.isfinite(B)
    Bz = np.where(fin, B, 0.0)
    rho = truncated_square_dispersion(trunc_c)
    L, R = Bz[:, :-1, :], Bz[:, 1:, :]
    fL, fR = fin[:, :-1, :], fin[:, 1:, :]
    if trunc_c is not None:
        prod = np.abs(Bz[:, :, 1:] * Bz[:, :, :-1])
        pok = fin[:, :, 1:] & fin[:, :, :-1]
        bpv = BPV_CONST * prod.sum(axis=2) / np.maximum(pok.sum(axis=2), 1)
        u2 = trunc_c**2 * np.maximum(bpv[:, :-1], bpv[:, 1:])[:, :, None]
        kL, kR = (L**2 <= u2) & fL, (R**2 <= u2) & fR
    else:
        kL, kR = fL, fR
    xL, xR = np.where(kL, L**2, 0.0), np.where(kR, R**2, 0.0)
    z = _pair_z(xL.sum(2), xR.sum(2), (xL**2).sum(2), (xR**2).sum(2), fL.sum(2), fR.sum(2), rho, norm)
    positions = (np.arange(1, m)) * k
    return z, positions


def bjv_overlap_statistics(r, k: int, trunc_c: float | None = 4.0, norm: str = "gaussian", chunk: int = 16,
                           max_cells: int = 40000):
    """Standardised comparisons of OVERLAPPING adjacent windows, one per possible change point.

    ``D[:, b]`` compares the left window ``r[b-k : b]`` with the right window
    ``r[b : b+k]`` exactly as in ``bjv_block_statistics`` (same normalisation, same
    pairwise truncation rule). It is defined for ``k <= b <= n - k`` and NaN
    elsewhere. Online use: ``D[:, b]`` is known at time ``b + k - 1``.
    """
    r = np.atleast_2d(np.asarray(r, float))
    reps, n = r.shape
    if n < 2 * k + 1:
        raise ValueError("series too short for two windows")
    rho = truncated_square_dispersion(trunc_c)
    D = np.full((reps, n), np.nan)
    fin_all = np.isfinite(r)
    rz_all = np.where(fin_all, r, 0.0)
    if trunc_c is None:  # cumulative sums: O(n)
        x = rz_all**2
        C1 = np.concatenate([np.zeros((reps, 1)), np.cumsum(x, axis=1)], axis=1)
        C2 = np.concatenate([np.zeros((reps, 1)), np.cumsum(x**2, axis=1)], axis=1)
        Cn = np.concatenate([np.zeros((reps, 1)), np.cumsum(fin_all, axis=1)], axis=1)
        b = np.arange(k, n - k + 1)
        S1L, S1R = C1[:, b] - C1[:, b - k], C1[:, b + k] - C1[:, b]
        S2L, S2R = C2[:, b] - C2[:, b - k], C2[:, b + k] - C2[:, b]
        NL, NR = Cn[:, b] - Cn[:, b - k], Cn[:, b + k] - Cn[:, b]
        D[:, b] = _pair_z(S1L, S1R, S2L, S2R, NL, NR, rho, norm)
        return D
    n_win = n - 2 * k + 1  # window w covers r[w : w + 2k] and has its boundary at b = w + k
    for i0 in range(0, reps, chunk):
        rows = slice(i0, i0 + chunk)
        n_rows = rz_all[rows].shape[0]
        step = max(1, max_cells // n_rows)
        for w0 in range(0, n_win, step):
            w1 = min(n_win, w0 + step)
            W = sliding_window_view(rz_all[rows, w0 : w1 + 2 * k - 1], 2 * k, axis=1)  # (rows, w1 - w0, 2k)
            F = sliding_window_view(fin_all[rows, w0 : w1 + 2 * k - 1], 2 * k, axis=1)
            prod = np.abs(W[:, :, 1:] * W[:, :, :-1])
            pok = F[:, :, 1:] & F[:, :, :-1]
            # bipower within each half (pairs fully inside the half)
            bl = BPV_CONST * prod[:, :, : k - 1].sum(2) / np.maximum(pok[:, :, : k - 1].sum(2), 1)
            br = BPV_CONST * prod[:, :, k:].sum(2) / np.maximum(pok[:, :, k:].sum(2), 1)
            u2 = trunc_c**2 * np.maximum(bl, br)[:, :, None]
            X = W**2
            Xk = np.where((X <= u2) & F, X, 0.0)
            S1L, S1R = Xk[:, :, :k].sum(2), Xk[:, :, k:].sum(2)
            S2L, S2R = (Xk[:, :, :k] ** 2).sum(2), (Xk[:, :, k:] ** 2).sum(2)
            NL, NR = F[:, :, :k].sum(2), F[:, :, k:].sum(2)
            D[rows, w0 + k : w1 + k] = _pair_z(S1L, S1R, S2L, S2R, NL, NR, rho, norm)
    return D


def bjv_test(
    r, k: int, alpha: float = 0.05, variant: str = "blocks", trunc_c: float | None = 4.0, norm: str = "gaussian",
    crit: float | None = None,
) -> BJVResult:
    """Offline (a posteriori) test for a volatility jump in a fixed window, with a change-point estimate.

    ``variant="blocks"``: maximum of ``|Z_j|`` over the ``m - 1`` neighbouring pairs of
    non-overlapping blocks, critical value ``gumbel_crit_max_abs_normal(m - 1, alpha)``.
    ``variant="overlap"``: maximum of ``|D_b|`` over all boundaries, critical value
    ``mosum_crit(n, k, alpha)``. ``crit`` overrides the asymptotic value (for example
    a value calibrated by simulation). The change-point estimate is the boundary
    with the largest ``|Z|`` (the first index after the estimated change) and the
    direction is the sign of Z there.
    """
    r = np.atleast_2d(np.asarray(r, float))
    n = r.shape[1]
    if variant == "blocks":
        z, pos = bjv_block_statistics(r, k, trunc_c, norm)
        c = gumbel_crit_max_abs_normal(z.shape[1], alpha) if crit is None else crit
    elif variant == "overlap":
        D = bjv_overlap_statistics(r, k, trunc_c, norm)
        pos = np.arange(k, n - k + 1)
        z = D[:, pos]
        c = mosum_crit(n, k, alpha) if crit is None else crit
    else:
        raise ValueError("variant must be 'blocks' or 'overlap'")
    az = np.where(np.isfinite(z), np.abs(z), -np.inf)
    j = np.argmax(az, axis=1)
    stat = az[np.arange(len(j)), j]
    direction = np.sign(z[np.arange(len(j)), j])
    return BJVResult(stat, float(c), stat > c, pos[j], direction, z, pos)


def bjv_online_statistic(r, k: int, trunc_c: float | None = 4.0, norm: str = "gaussian") -> np.ndarray:
    """Online version: entry t compares the last k returns ``r[t-k+1 .. t]`` with the k before them.

    ``out[:, t] = D[:, t - k + 1]`` from ``bjv_overlap_statistics`` (NaN for t < 2k - 1). It uses
    data up to time t only, so ``|out[:, t]| >= h`` is a valid online alarm rule. It is a
    repeated significance test with no anytime guarantee: its false alarm rate per day
    must be measured, not assumed.
    """
    D = bjv_overlap_statistics(r, k, trunc_c, norm)
    out = np.full_like(D, np.nan)
    out[:, k - 1 :] = D[:, : D.shape[1] - k + 1]
    return out


# ----------------------------------------------------------------------------
# Bayesian online change-point detection (Adams and MacKay 2007)
# ----------------------------------------------------------------------------


class BOCPDVariance:
    """Run-length filter of Adams and MacKay (2007, Algorithm 1) for changes in variance.

    Model inside a run: ``r_t | s2 ~ N(0, s2)`` with ``s2 ~ InvGamma(alpha0, beta0)``. The
    posterior after observations ``x_1..x_m`` of the run is ``InvGamma(alpha0 + m/2,
    beta0 + sum x^2 / 2)`` and the predictive density of the next return is Student t with
    ``2 alpha`` degrees of freedom and scale ``sqrt(beta / alpha)``:

    ``log p(x) = lgamma(alpha + 1/2) - lgamma(alpha) - log(2 pi beta) / 2 - (alpha + 1/2) log(1 + x^2 / (2 beta))``.

    Recursion (Algorithm 1 of the paper, constant hazard ``H``): with ``pi_r`` the
    predictive of the new return under run length r,
    growth ``P(r_t = r + 1, x_{1:t}) = P(r_{t-1} = r, x_{1:t-1}) pi_r (1 - H)`` and
    change ``P(r_t = 0, x_{1:t}) = sum_r P(r_{t-1} = r, x_{1:t-1}) pi_r H``, then normalise.
    The new run (r_t = 0) starts from the prior; the other runs absorb ``x_t``.

    Pruning: run lengths are capped at ``r_max - 1``. The last bucket means "a run of at
    least r_max - 1 steps"; when the run of length ``r_max - 2`` grows into it, the two are
    merged (probabilities added) and the merged bucket keeps the sufficient statistics of
    whichever of the two had more posterior mass. This only touches very long runs and
    never the short runs the alarm rule looks at.

    Works on ``reps`` independent paths at once. ``update(x)`` takes one return per path
    (NaN means "no observation", the state is unchanged) and returns
    ``(p_short, direction)``: the posterior probability that the current run is shorter
    than ``short`` steps, and +1 / -1 if the variance of the short runs is above / below
    that of the long runs (posterior means).
    """

    def __init__(self, reps: int, hazard: float, alpha0: float = 2.0, beta0=1.0, r_max: int = 256, short: int = 60):
        if not 0 < hazard < 1:
            raise ValueError("hazard must lie in (0, 1)")
        if alpha0 <= 1:
            raise ValueError("alpha0 must exceed 1 (finite prior mean of the variance)")
        if not 2 <= short < r_max - 1:
            raise ValueError("need 2 <= short < r_max - 1")
        self.reps, self.R, self.w = int(reps), int(r_max), int(short)
        self.logH, self.log1mH = np.log(hazard), np.log1p(-hazard)
        self.alpha0 = float(alpha0)
        self.beta0 = np.broadcast_to(np.asarray(beta0, float), (self.reps,)).copy()
        a = self.alpha0 + 0.5 * np.arange(self.R - 1)  # regular buckets: alpha depends only on the run length
        self.alpha_reg = a
        self.const_reg = gammaln(a + 0.5) - gammaln(a)
        self.logR = np.full((self.reps, self.R), -np.inf)
        self.logR[:, 0] = 0.0
        self.beta = np.tile(self.beta0[:, None], (1, self.R))
        self.alpha_last = np.full(self.reps, self.alpha0)
        self.t = 0

    def warm_start(self, n0, sum_x2):
        """Start with all mass on one long run that has already seen ``n0`` returns with sum of squares ``sum_x2``.

        This encodes "the calibration window has no change", the same assumption the
        e-detector makes when it estimates its baseline.
        """
        n0 = np.broadcast_to(np.asarray(n0, float), (self.reps,))
        s = np.broadcast_to(np.asarray(sum_x2, float), (self.reps,))
        self.logR[:] = -np.inf
        self.logR[:, -1] = 0.0
        self.alpha_last = self.alpha0 + 0.5 * n0
        self.beta[:, -1] = self.beta0 + 0.5 * s

    def set_prior_scale(self, beta0):
        """Change the prior of FUTURE new runs (an empirical-Bayes, predictable prior)."""
        self.beta0 = np.broadcast_to(np.asarray(beta0, float), (self.reps,)).copy()

    def update(self, x):
        x = np.broadcast_to(np.asarray(x, float), (self.reps,))
        ok = np.isfinite(x)
        x2 = np.where(ok, x * x, 0.0)
        R = self.R
        beta = self.beta
        lp = np.empty((self.reps, R))
        lp[:, :-1] = (self.const_reg[None, :] - 0.5 * np.log(2 * np.pi * beta[:, :-1])
                      - (self.alpha_reg[None, :] + 0.5) * np.log1p(x2[:, None] / (2 * beta[:, :-1])))
        aL = self.alpha_last
        lp[:, -1] = (gammaln(aL + 0.5) - gammaln(aL) - 0.5 * np.log(2 * np.pi * beta[:, -1])
                     - (aL + 0.5) * np.log1p(x2 / (2 * beta[:, -1])))
        lj = self.logR + lp
        cp = logsumexp(lj, axis=1) + self.logH
        g = lj + self.log1mH
        newR = np.empty_like(lj)
        newR[:, 0] = cp
        newR[:, 1:-1] = g[:, : R - 2]
        newR[:, -1] = np.logaddexp(g[:, R - 2], g[:, R - 1])
        newR -= logsumexp(newR, axis=1)[:, None]
        take_new = g[:, R - 2] > g[:, R - 1]  # merged bucket keeps the statistics of the heavier run
        newB = np.empty_like(beta)
        newB[:, 0] = self.beta0
        newB[:, 1:-1] = beta[:, : R - 2] + 0.5 * x2[:, None]
        newB[:, -1] = np.where(take_new, beta[:, R - 2], beta[:, R - 1]) + 0.5 * x2
        new_aL = np.where(take_new, self.alpha_reg[R - 2], aL) + 0.5
        # rows without an observation keep their state
        self.logR = np.where(ok[:, None], newR, self.logR)
        self.beta = np.where(ok[:, None], newB, beta)
        self.alpha_last = np.where(ok, new_aL, aL)
        self.t += 1
        return self.short_mass(), self.direction()

    def short_mass(self) -> np.ndarray:
        """Posterior probability that the current run is shorter than ``short`` steps."""
        return np.exp(logsumexp(self.logR[:, : self.w], axis=1))

    def short_log_odds(self) -> np.ndarray:
        """``log(p_short / (1 - p_short))`` computed in the log domain (no rounding to 0 or 1)."""
        return logsumexp(self.logR[:, : self.w], axis=1) - logsumexp(self.logR[:, self.w :], axis=1)

    def direction(self) -> np.ndarray:
        """+1 if short runs (1 to short - 1 steps) imply a higher variance than long runs, else -1."""
        a = np.concatenate([np.broadcast_to(self.alpha_reg, (self.reps, self.R - 1)), self.alpha_last[:, None]], axis=1)
        v = self.beta / (a - 1.0)  # posterior mean of the variance per run length
        p = np.exp(self.logR)
        ps, pl = p[:, 1 : self.w], p[:, self.w :]
        with np.errstate(invalid="ignore", divide="ignore"):
            vs = (ps * v[:, 1 : self.w]).sum(1) / ps.sum(1)
            vl = (pl * v[:, self.w :]).sum(1) / pl.sum(1)
        return np.where(vs >= vl, 1.0, -1.0)


def bocpd_paths(
    r, hazard: float, cal=None, alpha0: float = 2.0, prior_var=None, r_max: int = 256, short: int = 60,
    prior_var_path=None, return_log_odds: bool = False,
):
    """Run ``BOCPDVariance`` over rows of ``r`` and return ``(p_short, direction)``, both ``(reps, n)``.

    With ``return_log_odds=True`` a third array holds ``log(p_short / (1 - p_short))``,
    which keeps its resolution when ``p_short`` rounds to 1 (useful for matching alarm counts).

    ``cal`` (optional, ``(reps, n_cal)``) is a change-free calibration window: it sets the
    prior mean of the variance (``prior_var`` = its mean square unless given) and warm-starts
    the filter with one long run that has seen it. ``prior_var_path`` (``(reps, n)``)
    replaces the prior mean of new runs step by step; it must be predictable (use only
    data before t).
    """
    r = np.atleast_2d(np.asarray(r, float))
    reps, n = r.shape
    if cal is not None:
        cal = np.atleast_2d(np.asarray(cal, float))
        fin = np.isfinite(cal)
        s2 = np.where(fin, cal**2, 0.0).sum(1)
        nc = fin.sum(1)
        if prior_var is None:
            prior_var = s2 / np.maximum(nc, 1)
    if prior_var is None:
        prior_var = 1.0
    prior_var = np.maximum(np.asarray(prior_var, float), 1e-300)
    beta0 = (alpha0 - 1.0) * np.broadcast_to(np.asarray(prior_var, float), (reps,))
    f = BOCPDVariance(reps, hazard, alpha0, beta0, r_max, short)
    if cal is not None:
        f.warm_start(nc, s2)
    P = np.empty((reps, n))
    S = np.empty((reps, n))
    L = np.empty((reps, n)) if return_log_odds else None
    for t in range(n):
        if prior_var_path is not None:
            pv = prior_var_path[:, t]
            f.set_prior_scale(np.where(np.isfinite(pv) & (pv > 0), (alpha0 - 1.0) * pv, f.beta0))
        P[:, t], S[:, t] = f.update(r[:, t])
        if L is not None:
            L[:, t] = f.short_log_odds()
    return (P, S, L) if return_log_odds else (P, S)


# ----------------------------------------------------------------------------
# Alarm bookkeeping shared by the online competitors
# ----------------------------------------------------------------------------


def first_crossing(stat, h: float, start: int = 0) -> np.ndarray:
    """First index ``t >= start`` with ``stat[:, t] >= h`` per row, -1 if none (NaN never crosses)."""
    s = np.atleast_2d(np.asarray(stat, float))[:, start:]
    hit = np.nan_to_num(s, nan=-np.inf) >= h
    return np.where(hit.any(axis=1), hit.argmax(axis=1) + start, -1)


def scan_alarms(up_stat, dn_stat, h_up: float, h_dn: float, start: int = 0, dead: int = 0):
    """Alarm times and directions on one series with a dead time after each alarm.

    An "up" alarm fires at t when ``up_stat[t] >= h_up``, a "down" alarm when
    ``dn_stat[t] >= h_dn`` (up wins a tie). After an alarm at t the next one can fire
    at ``t + dead + 1`` at the earliest (the same pause as the e-detector's
    recalibration). Returns a list of ``(t, "up" | "down")``.
    """
    up = np.nan_to_num(np.asarray(up_stat, float), nan=-np.inf)
    dn = np.nan_to_num(np.asarray(dn_stat, float), nan=-np.inf)
    cu, cd = up >= h_up, dn >= h_dn
    cand = np.flatnonzero(cu | cd)
    out, nxt = [], start
    while True:
        j = np.searchsorted(cand, nxt)
        if j >= cand.size:
            break
        t = int(cand[j])
        out.append((t, "up" if cu[t] else "down"))
        nxt = t + dead + 1
    return out


def match_alarm_counts(up_stat, dn_stat, n_up: int, n_dn: int, start: int, dead: int, lo_hi=(0.0, 50.0), iters: int = 40):
    """Thresholds ``(h_up, h_dn)`` whose alarm counts are as close as possible to ``(n_up, n_dn)``.

    Coordinate-wise bisection, as for the rate-matched ratio rule of E6b. Both
    statistics must be "larger means stronger evidence".
    """
    lo0, hi0 = lo_hi
    h_up, h_dn = hi0, hi0

    def counts(hu, hd):
        ev = scan_alarms(up_stat, dn_stat, hu, hd, start, dead)
        u = sum(1 for e in ev if e[1] == "up")
        return u, len(ev) - u

    for _ in range(4):
        lo, hi = lo0, hi0
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if counts(mid, h_dn)[0] > n_up else (lo, mid)
        h_up = hi
        lo, hi = lo0, hi0
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if counts(h_up, mid)[1] > n_dn else (lo, mid)
        h_dn = hi
    return h_up, h_dn
