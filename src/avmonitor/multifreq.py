"""Multi-frequency volatility monitor: e-detectors at several sampling frequencies combined into one.

A k-minute detector bets on the clipped squared k-minute block return
``X_k = min(R_k^2, tau_k^2)``, where ``R_k`` is the sum of k one-minute returns
over non-overlapping blocks (optionally without the first return of each block,
``gap=True``, so that consecutive blocks share no price; this is the
non-overlapping block design of the pre-averaging alternative, which keeps the
null under i.i.d. price noise). Each frequency has its own truncation level
``tau_k = c * scaled MAD`` of k-minute returns and its own baseline ``m_k``, both
estimated on the SAME calibration window, and its own increase and decrease
Shiryaev-Roberts statistics (``docs/math.md`` 5.4). The monitor combines the
frequencies with weights ``w_k`` (summing to 1) and alarms when the combination
of either side reaches ``2 A``.

Why the combination is not simply "average the statistics every minute"
-----------------------------------------------------------------------
A k-minute statistic changes only at the end of a block, and its bet at that
moment depends on the k returns of the block. On the one-minute filtration it is
therefore NOT a supermartingale: one minute before the block ends, most of the
block is already known, so the conditional mean of the bet is not at most 1.
Averaging the frozen ("stale") statistics every minute is a legitimate monitor,
but Lemma 1 of ``docs/math.md`` does not apply to it (``docs/theory_multifreq.md``,
Proposition 2). Three combinations are implemented, each with a proved bound on
P(false alarm within m monitored minutes) under each frequency's null with the
true baseline:

``combine="average", check="certified"`` (default)
    Every minute, each coarse frequency enters with a *certified* value: its
    full statistic at the minute its block ends, and between block ends the
    smallest value its statistic can take at the next block end,
    ``mean_j min(S_j, l_j (S_j + 1))`` with ``l_j = 1 - lambda_j`` (increase side)
    or ``1 - lambda_j (B_k - 1)`` (decrease side). These values are below a
    supermartingale on the one-minute filtration (Theorem 1), so
    ``P <= sum_k w_k ceil(m / k) / A``. At minutes where every block ends at once
    (multiples of the coarsest frequency when the frequencies divide it) the
    certified average equals the plain average.
``combine="average", check="stale"``
    The plain weighted average of the frozen statistics, every minute. Proved
    only through a union bound: ``P <= (sum_k sqrt(w_k floor(m / k)))^2 / A``
    (Theorem 3), which is larger.
``combine="max"``
    Weighted Bonferroni: alarm when ``max_k w_k M_k >= 2 A``. Each frequency is
    a valid e-detector on its own block filtration, so
    ``P <= sum_k w_k floor(m / k) / A`` (Theorem 3).

``level_for_budget`` returns the level A that makes the proved bound equal a
target, so that monitors with different frequencies and weights can be compared
at the same proved false alarm budget.

Weights: ``"equal"`` (1 / J each), ``"balanced"`` (proportional to k, which
gives every frequency the same handicap J against its own single-frequency
monitor at the same budget; ``handicap``), ``"sqrt"`` (proportional to sqrt(k)),
or any nonnegative sequence.

Missing returns: a block in which one of the returns it sums is missing is
skipped (its statistic is unchanged), as the one-minute monitor skips a missing
return; the certified value covers this case. With ``gap=True`` a missing return
in the left-out first minute of a block does not skip the block. In a
calibration window, only the windows that contain a missing return are dropped.

Separate levels for the two sides: ``threshold_down`` (default: the same as
``threshold``) lets the decrease side alarm at ``2 A_down`` while the increase
side alarms at ``2 A``. Every bound below then holds with A replaced by
``1 / (1 / (2 A) + 1 / (2 A_down))`` (``docs/theory_multifreq.md``, remarks after
Theorem 1); this is used only for rate-matched comparisons in E17.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID
from .robust import MAD_NORMAL_SCALE

DEFAULT_FREQS = (1, 5, 15)


# ----------------------------------------------------------------------------
# Weights, proved bounds and levels
# ----------------------------------------------------------------------------


def resolve_weights(weights, freqs) -> np.ndarray:
    """Weights as an array summing to 1. Strings: "equal", "balanced" (proportional to k), "sqrt"."""
    k = np.asarray(freqs, float)
    if isinstance(weights, str):
        if weights == "equal":
            w = np.ones_like(k)
        elif weights == "balanced":
            w = k.copy()
        elif weights == "sqrt":
            w = np.sqrt(k)
        else:
            raise ValueError("weights must be 'equal', 'balanced', 'sqrt' or a sequence")
    else:
        w = np.asarray(weights, float)
        if w.shape != k.shape:
            raise ValueError("one weight per frequency is needed")
    if np.any(w < 0) or not np.all(np.isfinite(w)) or w.sum() <= 0:
        raise ValueError("weights must be finite, nonnegative and not all zero")
    return w / w.sum()


def _check_freqs(freqs) -> tuple:
    f = tuple(int(k) for k in freqs)
    if not f or any(k < 1 for k in f) or len(set(f)) != len(f):
        raise ValueError("frequencies must be distinct positive integers")
    return f


def proved_bound(m: int, threshold: float, freqs=DEFAULT_FREQS, weights="balanced", combine: str = "average",
                 check: str = "certified") -> float:
    """Proved bound on P(false alarm within m monitored minutes) for the two-sided monitor at level ``threshold``.

    Each side alarms at ``2 * threshold``. Under each frequency's two-sided null with the true baseline
    (``docs/theory_multifreq.md``): certified average ``sum_k w_k ceil(m / k) / A`` (Theorem 1); weighted
    maximum ``sum_k w_k floor(m / k) / A``; stale average ``(sum_k sqrt(w_k floor(m / k)))^2 / A``
    (Theorem 3). A frequency with weight 0 is left out.
    """
    f = _check_freqs(freqs)
    w = resolve_weights(weights, f)
    k = np.asarray(f, float)
    if combine == "max":
        return float(np.sum(w * np.floor(m / k)) / threshold)
    if combine != "average":
        raise ValueError("combine must be 'average' or 'max'")
    if check == "certified":
        return float(np.sum(w * np.ceil(m / k)) / threshold)
    if check == "stale":
        return float(np.sum(np.sqrt(w * np.floor(m / k))) ** 2 / threshold)
    raise ValueError("check must be 'certified' or 'stale'")


def level_for_budget(alpha: float, m: int, freqs=DEFAULT_FREQS, weights="balanced", combine: str = "average",
                     check: str = "certified") -> float:
    """The level A for which ``proved_bound(m, A, ...)`` equals ``alpha`` (each side then alarms at 2A)."""
    if not 0 < alpha:
        raise ValueError("alpha must be positive")
    return proved_bound(m, 1.0, freqs, weights, combine, check) / alpha


def handicap(freqs=DEFAULT_FREQS, weights="balanced", m: int | None = None) -> np.ndarray:
    """Factor by which frequency k alone must exceed its own single-frequency level in the certified average.

    With A set by ``level_for_budget`` for the average and ``A_k = ceil(m / k) / alpha`` for frequency k alone,
    the average alarms on evidence from k only when ``M_k >= h_k * 2 A_k`` with
    ``h_k = sum_j w_j ceil(m / j) / (w_k ceil(m / k))``; for m a multiple of every k this is
    ``k * sum_j (w_j / j) / w_k``. Balanced weights give ``h_k = J`` for every k, the smallest possible
    worst case (``docs/theory_multifreq.md``, section 5). ``m=None`` uses the formula for m divisible by all k.
    """
    f = _check_freqs(freqs)
    w = resolve_weights(weights, f)
    k = np.asarray(f, float)
    n = (1.0 / k) if m is None else np.ceil(m / k)
    with np.errstate(divide="ignore"):
        return np.sum(w * n) / (w * n)


# ----------------------------------------------------------------------------
# Block statistics and vectorised paths (used by the experiments and the tests)
# ----------------------------------------------------------------------------


def block_length(k: int, gap: bool) -> int:
    """Number of one-minute returns summed in a k-minute block (k - 1 with a gap, never for k = 1)."""
    return k - 1 if (gap and k > 1) else k


def block_returns(r, k: int, gap: bool = False) -> np.ndarray:
    """Non-overlapping k-minute block returns along the last axis.

    Block b sums ``r[..., b k : (b + 1) k]``; with ``gap=True`` (k >= 2) the first return of each block is left
    out, so consecutive blocks share no price. An incomplete last block is dropped; NaN propagates.
    """
    r = np.asarray(r, float)
    nb = r.shape[-1] // k
    blk = r[..., : nb * k].reshape(r.shape[:-1] + (nb, k))
    if gap and k > 1:
        blk = blk[..., 1:]
    return blk.sum(axis=-1)


def calibration_block_returns(r, k: int, gap: bool = False, overlap: bool = True) -> np.ndarray:
    """k-minute returns of a calibration window along the last axis.

    ``overlap=True`` uses every window of ``block_length(k, gap)`` consecutive returns (all phases of the block
    grid; more data for the plug-in estimates), otherwise the non-overlapping blocks of ``block_returns``.
    """
    r = np.asarray(r, float)
    L = block_length(k, gap)
    if not overlap or L == 1:
        return block_returns(r, k, gap) if not overlap else r
    # A missing return makes only the windows that contain it missing (a plain cumulative sum would make every
    # later window missing).
    fin = np.isfinite(r)
    z = np.zeros(r.shape[:-1] + (1,))
    c = np.cumsum(np.concatenate([z, np.where(fin, r, 0.0)], axis=-1), axis=-1)
    nbad = np.cumsum(np.concatenate([z, (~fin).astype(float)], axis=-1), axis=-1)
    out = c[..., L:] - c[..., :-L]
    return np.where(nbad[..., L:] - nbad[..., :-L] > 0, np.nan, out)


def plugin_calibration(r_cal, k: int, c: float = 3.0, gap: bool = False, overlap: bool = True,
                       null_band: float = 0.0):
    """Vectorised plug-in baseline of frequency k for each row of ``r_cal``.

    Returns ``(scale, tau2, m0, m_up, m_dn)`` with ``scale`` = scaled median |R_k|, ``tau2 = (c scale)^2``,
    ``m0`` = mean of ``min(R_k^2, tau2)``, ``m_up = m0 (1 + b)``, ``m_dn = m0 (1 - b)``. Windows with a missing
    return are left out, as in the streaming monitor (no extension or fallback here).
    """
    R = calibration_block_returns(r_cal, k, gap, overlap)
    med, mean = (np.nanmedian, np.nanmean) if np.isnan(R).any() else (np.median, np.mean)
    scale = med(np.abs(R), axis=-1) * MAD_NORMAL_SCALE
    tau2 = (c * scale) ** 2
    m0 = mean(np.minimum(R**2, tau2[..., None]), axis=-1)
    return scale, tau2, m0, m0 * (1 + null_band), m0 * (1 - null_band)


def component_paths(x, m_up, m_dn, tau2, up_grid=DEFAULT_UP_GRID, down_fractions=DEFAULT_DOWN_FRACTIONS,
                    cap: float = 1e12) -> dict:
    """Increase and decrease Shiryaev-Roberts mixtures of one frequency, at each block end.

    ``x`` (reps, nb): clipped squared block returns (NaN = skipped block). ``m_up``, ``m_dn``, ``tau2``: scalars
    or (reps,). Returns arrays of shape (reps, nb): ``up_full[:, b]`` is the grid mean of the increase statistics
    after block b, ``up_cert[:, b]`` the certified value used strictly inside block b + 1,
    ``mean_j min(S_j, (1 - lambda_j)(S_j + 1))``, and the same for the decrease side with
    ``1 - lambda_j (B - 1)``. Statistics are capped at ``cap`` (only to avoid overflow long after an alarm).
    A path whose decrease side has no admissible bet (``tau2 <= m_dn``) keeps a decrease statistic of 0.
    """
    x = np.atleast_2d(np.asarray(x, float))
    reps, nb = x.shape
    m_up = np.broadcast_to(np.asarray(m_up, float), (reps,))
    m_dn = np.broadcast_to(np.asarray(m_dn, float), (reps,))
    tau2 = np.broadcast_to(np.asarray(tau2, float), (reps,))
    gu = np.asarray(up_grid, float)
    fd = np.asarray(down_fractions, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        B = tau2 / m_dn
    active = np.isfinite(B) & (B > 1)
    lam_dn = np.where(active[:, None], fd[None, :] / np.where(active, B - 1, 1.0)[:, None], 0.0)
    l_up = 1.0 - gu
    l_dn = 1.0 - fd
    Su = np.zeros((reps, gu.size))
    Sd = np.zeros((reps, fd.size))
    out = {key: np.empty((reps, nb)) for key in ("up_full", "up_cert", "dn_full", "dn_cert")}
    for b in range(nb):
        xb = x[:, b]
        ok = np.isfinite(xb)
        xb = np.where(ok, xb, 0.0)
        fu = 1.0 + gu[None, :] * (xb / m_up - 1.0)[:, None]
        fdn = 1.0 + lam_dn * (1.0 - xb / m_dn)[:, None]
        Su = np.where(ok[:, None], np.minimum(fu * (Su + 1.0), cap), Su)
        Sd = np.where((ok & active)[:, None], np.minimum(fdn * (Sd + 1.0), cap), Sd)
        out["up_full"][:, b] = Su.mean(axis=1)
        out["dn_full"][:, b] = Sd.mean(axis=1)
        out["up_cert"][:, b] = np.minimum(Su, l_up * (Su + 1.0)).mean(axis=1)
        out["dn_cert"][:, b] = np.where(active, np.minimum(Sd, l_dn * (Sd + 1.0)).mean(axis=1), 0.0)
    return out


def minute_values(full, cert, k: int, n: int, certified: bool = True) -> np.ndarray:
    """Per-minute value of one frequency's statistic for minutes 0 .. n-1 of a monitored segment.

    Minute i ends block ``(i + 1) // k - 1`` when ``(i + 1) % k == 0``; there the value is the full statistic.
    Elsewhere it is the certified value after the last completed block (``certified=True``) or the frozen full
    value (``certified=False``). Before the first block ends the value is 0.
    """
    full = np.atleast_2d(full)
    cert = np.atleast_2d(cert)
    reps = full.shape[0]
    i = np.arange(n)
    nb = (i + 1) // k
    pad_full = np.concatenate([np.zeros((reps, 1)), full], axis=1)
    pad_cert = np.concatenate([np.zeros((reps, 1)), cert], axis=1)
    stale = pad_full[:, nb]
    if k == 1 or not certified:
        return stale
    at_end = (i + 1) % k == 0
    return np.where(at_end[None, :], stale, pad_cert[:, nb])


def first_alarm(up, dn, threshold: float):
    """First minute at which ``up`` or ``dn`` (reps, n) reaches ``2 * threshold``; -1 if none.

    Returns ``(index, direction)`` with direction 1 for up, -1 for down, 0 for none (up wins ties, as in
    ``VolatilityMonitor``).
    """
    hu = up >= 2 * threshold
    hd = dn >= 2 * threshold
    hit = hu | hd
    any_ = hit.any(axis=1)
    idx = np.where(any_, hit.argmax(axis=1), -1)
    rows = np.arange(up.shape[0])
    safe = np.maximum(idx, 0)
    direction = np.where(any_, np.where(hu[rows, safe], 1, -1), 0)
    return idx, direction


# ----------------------------------------------------------------------------
# Streaming monitor
# ----------------------------------------------------------------------------


@dataclass
class MFEvent:
    index: int  # position in the stream of the alarm
    direction: str  # "up" or "down"
    monitored_since: int  # first monitored index of this segment
    baseline_scale: float  # robust sd of the ONE-minute returns of the calibration window
    stat: float  # combined statistic of the alarming side
    shares: dict = field(default_factory=dict)  # share of each frequency in the combined statistic at the alarm


def _usable(*vals) -> bool:
    return all(math.isfinite(v) and v > 0 for v in vals)


class MultiFrequencyMonitor:
    """Calibrate, monitor, alarm, recalibrate loop with e-detectors at several frequencies.

    Parameters follow ``VolatilityMonitor`` (plug-in baseline only). ``threshold`` is A; each side alarms when
    its combined statistic reaches 2A (the decrease side at ``2 * threshold_down`` if that is given). ``freqs`` are block lengths in minutes; ``weights`` see
    ``resolve_weights``; ``combine`` and ``check`` see the module docstring. ``gap=True`` leaves out the first
    return of every block of length k >= 2. ``cal_overlap=True`` estimates each coarse frequency's MAD and
    baseline from all windows of the block length in the calibration window. ``known_baseline`` (for
    simulations): a dict ``{k: (tau2, m0)}`` used instead of the plug-in estimates (the calibration window
    is still consumed). With ``freqs=(1,)`` the monitor gives exactly the events of
    ``VolatilityMonitor(threshold, c, cal_len, null_band, ...)`` with the plug-in baseline.
    """

    def __init__(self, threshold: float = 28800.0, freqs=DEFAULT_FREQS, weights="balanced", c: float = 3.0,
                 cal_len: int = 1440, null_band: float = 0.0, up_grid=DEFAULT_UP_GRID,
                 down_fractions=DEFAULT_DOWN_FRACTIONS, scale_floor: float = 0.0, *, combine: str = "average",
                 check: str = "certified", gap: bool = False, cal_overlap: bool = True, cal_max: int | None = None,
                 known_baseline: dict | None = None, trace: bool = False, threshold_down: float | None = None):
        self.A = float(threshold)
        self.A_down = self.A if threshold_down is None else float(threshold_down)
        self.freqs = _check_freqs(freqs)
        self.w = resolve_weights(weights, self.freqs)
        if combine not in ("average", "max"):
            raise ValueError("combine must be 'average' or 'max'")
        if check not in ("certified", "stale"):
            raise ValueError("check must be 'certified' or 'stale'")
        self.combine = combine
        self.check = check
        self.c = float(c)
        self.cal_len = int(cal_len)
        self.cal_max = int(cal_max) if cal_max is not None else 2 * int(cal_len)
        if self.cal_max < self.cal_len:
            raise ValueError("cal_max must be at least cal_len")
        self.band = float(null_band)
        self.up_grid = np.asarray(up_grid, float)
        self.down_fractions = np.asarray(down_fractions, float)
        self.l_up = 1.0 - self.up_grid
        self.l_dn = 1.0 - self.down_fractions
        self.scale_floor = float(scale_floor)
        self.gap = bool(gap)
        self.cal_overlap = bool(cal_overlap)
        self.known = None if known_baseline is None else {int(k): tuple(map(float, v)) for k, v in known_baseline.items()}
        if self.known is not None and set(self.known) != set(self.freqs):
            raise ValueError("known_baseline needs one (tau2, m0) per frequency")
        self.L = [block_length(k, self.gap) for k in self.freqs]
        self.t = -1
        self.events: list[MFEvent] = []
        self.segments: list[dict] = []
        self.cal_stats = {"calibrations": 0, "extended": 0, "fallback": 0, "restarted": 0}
        self.trace = trace
        self.trace_rows: list[tuple] = []
        self._reset()

    def proved_bound(self, m: int) -> float:
        """Proved bound on P(false alarm within m monitored minutes) for this configuration (true baseline)."""
        a_eff = 1.0 / (0.5 / self.A + 0.5 / self.A_down)  # equals A when both sides use the same level
        return proved_bound(m, a_eff, self.freqs, self.w, self.combine, self.check)

    def _reset(self):
        self.buffer: list[float] = []
        self.calibrated = False
        J = len(self.freqs)
        self.S_up = [np.zeros(self.up_grid.size) for _ in range(J)]
        self.S_dn = [np.zeros(self.down_fractions.size) for _ in range(J)]
        self.full_up = np.zeros(J)
        self.full_dn = np.zeros(J)
        self.cert_up = np.zeros(J)
        self.cert_dn = np.zeros(J)
        self.acc = np.zeros(J)
        self.bad = np.zeros(J, bool)
        self.pos = 0
        self._extended = False
        self._fallback = False

    # ------------------------------------------------------------------
    # Calibration (plug-in, the same window for every frequency)
    # ------------------------------------------------------------------

    def _block_cal(self, r: np.ndarray, j: int) -> np.ndarray:
        k = self.freqs[j]
        R = r if k == 1 else calibration_block_returns(r, k, self.gap, self.cal_overlap)
        return R[np.isfinite(R)]

    def _stats(self, R: np.ndarray, j: int):
        if R.size == 0:
            return 0.0, 0.0, 0.0
        floor = self.scale_floor * math.sqrt(self.L[j])
        scale = max(float(np.median(np.abs(R)) * MAD_NORMAL_SCALE), floor)
        tau2 = (self.c * scale) ** 2
        m0 = float(np.mean(np.minimum(R**2, tau2)))
        return scale, tau2, m0

    def _fallback_scale(self, R: np.ndarray, j: int) -> float:
        nz = np.abs(R[R != 0])
        if nz.size == 0:
            return 0.0
        return max(float(np.median(nz) * MAD_NORMAL_SCALE), self.scale_floor * math.sqrt(self.L[j]))

    def _calibrate(self):
        rb = np.asarray(self.buffer, float)
        rf = rb[np.isfinite(rb)]
        J = len(self.freqs)
        Rs = [self._block_cal(rb, j) for j in range(J)]
        st = [self._stats(Rs[j], j) for j in range(J)]
        usable = all(_usable(*s) for s in st)
        if not (usable and (not self._extended or rf.size >= self.cal_len // 2)):
            if len(self.buffer) < self.cal_max:
                self._extended = True
                return
            if not usable:
                for j in range(J):
                    if not _usable(*st[j]):
                        sc = self._fallback_scale(Rs[j], j)
                        tau2 = (self.c * sc) ** 2
                        m0 = float(np.mean(np.minimum(Rs[j] ** 2, tau2))) if Rs[j].size else 0.0
                        if not _usable(sc, tau2, m0):
                            self.cal_stats["restarted"] += 1
                            self.buffer = []
                            return
                        st[j] = (sc, tau2, m0)
                self._fallback = True
        if self.known is not None:
            st = [(math.sqrt(self.known[k][0]) / self.c, self.known[k][0], self.known[k][1]) for k in self.freqs]
        self.scale1 = float(np.median(np.abs(rf)) * MAD_NORMAL_SCALE) if rf.size else float("nan")
        self.tau2 = np.array([s[1] for s in st])
        self.m0 = np.array([s[2] for s in st])
        self.m_up = self.m0 * (1 + self.band)
        self.m_dn = self.m0 * (1 - self.band)
        self.lam_dn, self.down_active = [], []
        for j in range(J):
            b = self.tau2[j] / self.m_dn[j] if self.m_dn[j] > 0 else math.inf
            act = bool(math.isfinite(b) and b > 1.0)
            self.down_active.append(act)
            self.lam_dn.append(self.down_fractions / (b - 1.0) if act else np.zeros(self.down_fractions.size))
        self.start = self.t + 1
        self.calibrated = True
        self.cal_stats["calibrations"] += 1
        self.cal_stats["extended"] += int(self._extended)
        self.cal_stats["fallback"] += int(self._fallback)
        self.segments.append(dict(start=self.start, cal_len=len(self.buffer), n_finite=int(rf.size),
                                  scale1=self.scale1, freqs=self.freqs, scale=[s[0] for s in st],
                                  tau2=self.tau2.tolist(), m0=self.m0.tolist(), down_active=list(self.down_active),
                                  extended=self._extended, fallback=self._fallback))

    # ------------------------------------------------------------------
    # Monitoring
    # ------------------------------------------------------------------

    def _bet(self, j: int, x: float):
        self.S_up[j] = (1.0 + self.up_grid * (x / self.m_up[j] - 1.0)) * (self.S_up[j] + 1.0)
        if self.down_active[j]:
            self.S_dn[j] = (1.0 + self.lam_dn[j] * (1.0 - x / self.m_dn[j])) * (self.S_dn[j] + 1.0)

    def _end_block(self, j: int):
        su, sd = self.S_up[j], self.S_dn[j]
        self.full_up[j] = su.mean()
        self.full_dn[j] = sd.mean()
        self.cert_up[j] = np.minimum(su, self.l_up * (su + 1.0)).mean()
        self.cert_dn[j] = np.minimum(sd, self.l_dn * (sd + 1.0)).mean() if self.down_active[j] else 0.0

    def update(self, r: float) -> MFEvent | None:
        """Feed one one-minute return; returns an ``MFEvent`` if an alarm fires at this minute."""
        self.t += 1
        if not self.calibrated:
            self.buffer.append(r)
            if len(self.buffer) >= self.cal_len:
                self._calibrate()
            return None
        i = self.pos
        self.pos += 1
        fin = math.isfinite(r)
        J = len(self.freqs)
        vu = np.empty(J)
        vd = np.empty(J)
        for j, k in enumerate(self.freqs):
            if k == 1:
                if fin:
                    self._bet(j, min(r * r, self.tau2[j]))
                    self.full_up[j] = self.S_up[j].mean()
                    self.full_dn[j] = self.S_dn[j].mean()
                vu[j], vd[j] = self.full_up[j], self.full_dn[j]
                continue
            off = i % k
            if off == 0:
                self.acc[j] = 0.0
                self.bad[j] = False
            if not (self.gap and off == 0):
                if fin:
                    self.acc[j] += r
                else:
                    self.bad[j] = True
            if off == k - 1:
                if not self.bad[j]:
                    a = self.acc[j]
                    self._bet(j, min(a * a, self.tau2[j]))
                self._end_block(j)
                vu[j], vd[j] = self.full_up[j], self.full_dn[j]
            elif self.combine == "average" and self.check == "certified":
                vu[j], vd[j] = self.cert_up[j], self.cert_dn[j]
            else:
                vu[j], vd[j] = self.full_up[j], self.full_dn[j]
        if self.combine == "average":
            mu, md = float(self.w @ vu), float(self.w @ vd)
        else:
            mu, md = float(np.max(self.w * vu)), float(np.max(self.w * vd))
        if self.trace:
            self.trace_rows.append((self.t, mu, md, vu.copy(), vd.copy()))
        ev = None
        if mu >= 2 * self.A:
            ev = self._event("up", mu, self.w * vu)
        elif md >= 2 * self.A_down:
            ev = self._event("down", md, self.w * vd)
        if ev is not None:
            self.events.append(ev)
            self._reset()
        return ev

    def _event(self, direction, stat, parts):
        tot = parts.sum()
        shares = {k: (float(p / tot) if tot > 0 else float("nan")) for k, p in zip(self.freqs, parts)}
        return MFEvent(self.t, direction, self.start, self.scale1, float(stat), shares)

    def run(self, returns) -> list[MFEvent]:
        for r in np.asarray(returns, float):
            self.update(float(r))
        return self.events
