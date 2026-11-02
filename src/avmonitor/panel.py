"""Monitoring a panel of K assets with false discovery control.

Theory, proofs and references: ``docs/theory_panel.md``. Summary:

* Each asset k gets a two-sided **mixture e-process** for the null "the
  conditional truncated variance of asset k has not moved away from its
  baseline since the block started" (``docs/math.md`` 5.1, with the
  conditioning on the information of the WHOLE panel). It is built from the
  same bet chains as the Shiryaev-Roberts e-detector of ``detectors.py``, but
  the chain that starts at step j gets a prior weight ``w_j`` with
  ``sum_j w_j = 1``, and the weight of chains that have not started yet is
  kept as cash:

      E_t = 1/2 (mean over the up grid + mean over the down grid) of S_t + T_t,
      S_t = f_t (S_{t-1} + w_t),   T_t = sum_{j > t} w_j.

  E is a nonnegative supermartingale that starts at 1 (a test supermartingale)
  as long as the null holds, so ``E_tau`` is an e-value at EVERY stopping time
  tau (optional stopping). The Shiryaev-Roberts statistic itself is not: its
  expectation grows like t.
* Two weight families: ``"window"`` (``w_j = 1/h`` for the first h own
  observations of the block, so ``E_t = M_t / h + (h - t) / h`` with the
  e-detector M: the process in the proof of Theorem B of ``docs/math.md``),
  and ``"pareto"`` (``w_j = h / ((j + h - 1)(j + h))``, an infinite horizon;
  a change after j steps costs ``log(1 / w_j)``, about ``log h + 2 log(j / h)``).
* At every step the panel applies **e-BH** (Wang and Ramdas 2022) to the
  e-processes, each STOPPED at the time its asset was discovered (``rule="ebh"``).
  Stopped test supermartingales are still test supermartingales, the
  discovery set only grows within a block, and every set it reports is
  self-consistent. Consequences (proved in the theory note): at every
  stopping time the false discovery rate is at most alpha, under any
  dependence between assets; ``P(sup_t FDP_t >= c) <= alpha / c``;
  ``E[sup_t FDP_t] <= alpha (1 + log(1 / alpha))``; with no change anywhere,
  P(any discovery in the block) <= alpha.
* ``rule="bonferroni"`` flags an asset when its e-process reaches ``K / alpha``
  (family-wise error at most alpha per block) and ``rule="naive"`` at
  ``1 / alpha`` (alpha per asset, no correction). Path by path, the naive rule
  flags no later than e-BH, and e-BH no later than Bonferroni.
* Blocks: the guarantee is per block. At the start of a block every asset is
  recalibrated on its own last ``cal_len`` finite returns (data strictly
  before the block, so the baseline is predictable) and gets a fresh
  e-process. Missing returns (NaN) leave an asset's e-process unchanged and
  do not advance its own clock, which is how assets with different trading
  hours share one panel clock (valid when the trading calendar, or more
  generally whether a return is observed, is known one step ahead).

Vectorised core: ``PanelEProcess`` (many independent panels at once, used by
the experiments) and ``run_block``. Streaming API: ``PanelMonitor``.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field

import numpy as np

from .detectors import DEFAULT_DOWN_FRACTIONS, DEFAULT_UP_GRID, calibrate

RULES = ("ebh", "bonferroni", "naive", "none")
WEIGHTS = ("window", "pareto")
_CAP = 1e250  # numerical cap for the chain sums (far above every alarm level used)


# ----------------------------------------------------------------------------
# e-BH
# ----------------------------------------------------------------------------


def ebh(e, alpha: float) -> np.ndarray:
    """e-BH rejection mask along the last axis (Wang and Ramdas 2022).

    With K = ``e.shape[-1]`` hypotheses, ``R = max{r : e_(r) >= K / (alpha r)}``
    (``e_(r)`` the r-th largest e-value, R = 0 if no r qualifies), and the R
    hypotheses with the largest e-values are rejected. The returned mask is
    ``e >= K / (alpha R)`` (exactly R entries, also with ties). NaN counts as 0.
    """
    e = np.nan_to_num(np.asarray(e, float), nan=0.0)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0, 1)")
    K = e.shape[-1]
    srt = -np.sort(-e, axis=-1)
    r = np.arange(1, K + 1)
    ok = srt * r >= K / alpha
    R = np.where(ok.any(axis=-1), K - np.argmax(ok[..., ::-1], axis=-1), 0)
    with np.errstate(divide="ignore"):
        thr = np.where(R > 0, K / (alpha * np.maximum(R, 1)), np.inf)
    return e >= thr[..., None]


def ebh_count(e, alpha: float) -> np.ndarray:
    """Number of e-BH rejections along the last axis."""
    return ebh(e, alpha).sum(axis=-1)


# ----------------------------------------------------------------------------
# Prior weights over the change time
# ----------------------------------------------------------------------------


def change_weights(n, kind: str, horizon):
    """Weight ``w_n`` of the chain that starts at own observation n (n >= 1) and the tail ``sum_{j > n} w_j``.

    ``kind="window"``: ``w_n = 1/h`` for ``n <= h`` and 0 after; tail ``max(h - n, 0) / h``.
    ``kind="pareto"``: ``w_n = h / ((n + h - 1)(n + h))``; tail ``h / (n + h)``. Both sum to 1.
    ``n`` and ``horizon`` broadcast. Returns ``(w, tail)`` as float arrays.
    """
    n = np.asarray(n, float)
    h = np.asarray(horizon, float)
    if np.any(h <= 0):
        raise ValueError("horizon must be positive")
    if kind == "window":
        w = np.where((n >= 1) & (n <= h), 1.0 / h, 0.0)
        tail = np.maximum(h - n, 0.0) / h
    elif kind == "pareto":
        nn = np.maximum(n, 1.0)
        w = np.where(n >= 1, h / ((nn + h - 1.0) * (nn + h)), 0.0)
        tail = h / (np.maximum(n, 0.0) + h)
    else:
        raise ValueError(f"weights must be one of {WEIGHTS}")
    return w, tail


# ----------------------------------------------------------------------------
# Vectorised mixture e-processes
# ----------------------------------------------------------------------------


class PanelEProcess:
    """Two-sided mixture e-processes for K assets on ``reps`` independent panels.

    Parameters broadcast to shape (reps, K): ``tau2`` (truncation level squared),
    ``m_up`` (null value of the increase side: ``E[X | past] <= m_up``), ``m_dn``
    (decrease side: ``E[X | past] >= m_dn``), ``horizon`` (h of the weights, in the
    asset's own observations) and ``active`` (False: the asset's e-process is the
    constant 1 for this block, for example when it could not be calibrated).
    The decrease side is switched off (contributes the constant 1/2) when
    ``tau2 / m_dn <= 1``, where no decrease bet is admissible.
    """

    def __init__(self, tau2, m_up, m_dn, horizon, weights: str = "window", reps: int | None = None,
                 K: int | None = None, active=None, up_grid=DEFAULT_UP_GRID, down_fractions=DEFAULT_DOWN_FRACTIONS):
        arrs = [np.asarray(a, float) for a in (tau2, m_up, m_dn, horizon)]
        shape = np.broadcast_shapes(*(a.shape for a in arrs), (reps or 1, K or 1))
        if len(shape) != 2:
            raise ValueError("parameters must broadcast to (reps, K)")
        tau2, m_up, m_dn, horizon = (np.broadcast_to(a, shape).copy() for a in arrs)
        act = np.ones(shape, bool) if active is None else np.broadcast_to(np.asarray(active, bool), shape).copy()
        act &= np.isfinite(tau2) & (tau2 > 0) & np.isfinite(m_up) & (m_up > 0)
        if weights not in WEIGHTS:
            raise ValueError(f"weights must be one of {WEIGHTS}")
        self.shape = shape
        self.weights = weights
        self.tau2 = np.where(act, tau2, 1.0)
        self.m_up = np.where(act, m_up, 1.0)
        m_dn = np.where(act & np.isfinite(m_dn) & (m_dn > 0), m_dn, np.inf)
        b = self.tau2 / m_dn
        self.down_active = act & (b > 1.0)
        self.m_dn = np.where(self.down_active, m_dn, 1.0)
        self.horizon = np.where(horizon > 0, horizon, 1.0)
        self.active = act
        self.lam_up = np.asarray(up_grid, float)
        if np.any(self.lam_up < 0) or np.any(self.lam_up >= 1):
            raise ValueError("increase bets must lie in [0, 1)")
        fr = np.asarray(down_fractions, float)
        if np.any(fr < 0) or np.any(fr >= 1):
            raise ValueError("decrease fractions must lie in [0, 1)")
        with np.errstate(divide="ignore", invalid="ignore"):
            lam_max_dn = np.where(self.down_active, 1.0 / (b - 1.0), 0.0)
        self.lam_dn = lam_max_dn[..., None] * fr  # (reps, K, Gd)
        self.S_up = np.zeros(shape + (self.lam_up.size,))
        self.S_dn = np.zeros(shape + (fr.size,))
        self.n = np.zeros(shape, np.int64)  # own clock: observed returns since the block started
        self.tail = np.ones(shape)
        self.frozen = np.zeros(shape, bool)
        self.frozen_value = np.ones(shape)

    def step(self, r) -> None:
        """Feed one return per asset and panel (array broadcastable to (reps, K)); NaN = not observed."""
        r = np.broadcast_to(np.asarray(r, float), self.shape)
        obs = np.isfinite(r) & self.active & ~self.frozen
        if not obs.any():
            return
        x = np.minimum(np.where(obs, r, 0.0) ** 2, self.tau2)
        self.n = self.n + obs
        w, tail = change_weights(self.n, self.weights, self.horizon)
        w = np.where(obs, w, 0.0)
        g_up = np.where(obs, x / self.m_up - 1.0, 0.0)
        g_dn = np.where(obs & self.down_active, 1.0 - x / self.m_dn, 0.0)
        f_up = 1.0 + self.lam_up * g_up[..., None]
        f_dn = 1.0 + self.lam_dn * g_dn[..., None]
        self.S_up = np.minimum(f_up * (self.S_up + w[..., None]), _CAP)
        self.S_dn = np.minimum(f_dn * (self.S_dn + w[..., None]), _CAP)
        self.tail = np.where(obs, tail, self.tail)

    @property
    def up_part(self) -> np.ndarray:
        return self.S_up.mean(axis=-1)

    @property
    def down_part(self) -> np.ndarray:
        return self.S_dn.mean(axis=-1)

    @property
    def live_value(self) -> np.ndarray:
        """Current e-value of every asset, ignoring freezing. Inactive assets: 1."""
        e = 0.5 * (self.up_part + self.down_part) + self.tail
        return np.where(self.active, e, 1.0)

    @property
    def value(self) -> np.ndarray:
        """e-values with discovered (frozen) assets held at their value at discovery."""
        return np.where(self.frozen, self.frozen_value, self.live_value)

    @property
    def direction(self) -> np.ndarray:
        """+1 where the increase side carries more evidence, -1 otherwise."""
        return np.where(self.up_part >= self.down_part, 1, -1)

    def freeze(self, mask) -> None:
        mask = np.asarray(mask, bool) & ~self.frozen
        self.frozen_value = np.where(mask, self.live_value, self.frozen_value)
        self.frozen |= mask


def flags(values: np.ndarray, frozen: np.ndarray, alpha: float, rule: str) -> np.ndarray:
    """New discoveries (mask) at this step for one of ``RULES``, given the (stopped) e-values."""
    K = values.shape[-1]
    if rule == "ebh":
        new = ebh(values, alpha)
    elif rule == "bonferroni":
        new = values >= K / alpha
    elif rule == "naive":
        new = values >= 1.0 / alpha
    elif rule == "none":
        return np.zeros(values.shape, bool)
    else:
        raise ValueError(f"rule must be one of {RULES}")
    return new & ~frozen


@dataclass
class BlockResult:
    """Outcome of ``run_block`` for (reps, K) asset-panels; index -1 means never discovered."""

    flag_time: np.ndarray  # step of discovery
    direction: np.ndarray  # +1 up, -1 down, 0 none
    value_at_flag: np.ndarray
    set_size_at_flag: np.ndarray  # number of discoveries in the panel right after this one was made
    final_value: np.ndarray  # stopped e-values at the end of the block
    path: np.ndarray | None = None  # stopped e-values at every step: (n, K) for panel 0, or (reps, n, K)
    dir_path: np.ndarray | None = None  # direction (+1 / -1) at every step, same shape as path


def run_block(r, tau2, m_up, m_dn, alpha: float = 0.05, rule: str = "ebh", horizon=None, weights: str = "window",
              active=None, up_grid=DEFAULT_UP_GRID, down_fractions=DEFAULT_DOWN_FRACTIONS,
              return_path: bool | str = False) -> BlockResult:
    """Run one block of the panel monitor on returns ``r`` of shape (reps, K, n) or (K, n).

    Parameters as in ``PanelEProcess``; ``horizon`` defaults to n (window weights
    spread over the block). Discovered assets are frozen (their e-process is
    stopped), so with ``rule="ebh"`` the discovery set only grows.
    ``return_path=True`` keeps the stopped e-values and directions of panel 0
    (shape (n, K)); ``return_path="all"`` keeps them for every panel (shape (reps, n, K)).
    """
    r = np.asarray(r, float)
    if r.ndim == 2:
        r = r[None]
    reps, K, n = r.shape
    horizon = n if horizon is None else horizon
    ep = PanelEProcess(tau2, m_up, m_dn, horizon, weights, reps=reps, K=K, active=active,
                       up_grid=up_grid, down_fractions=down_fractions)
    flag_time = np.full((reps, K), -1)
    direction = np.zeros((reps, K), int)
    vflag = np.full((reps, K), np.nan)
    size = np.zeros((reps, K), int)
    keep_all = isinstance(return_path, str)
    if keep_all and return_path != "all":
        raise ValueError("return_path must be True, False or 'all'")
    pshape = (reps, n, K) if keep_all else (n, K)
    path = np.empty(pshape) if return_path else None
    dir_path = np.empty(pshape, np.int8) if return_path else None
    for t in range(n):
        ep.step(r[:, :, t])
        v = ep.value
        new = flags(v, ep.frozen, alpha, rule)
        if new.any():
            flag_time[new] = t
            direction[new] = ep.direction[new]
            vflag[new] = v[new]
            ep.freeze(new)
            size = np.where(new, ep.frozen.sum(axis=1, keepdims=True), size)
        if path is not None and keep_all:
            path[:, t] = ep.value
            dir_path[:, t] = ep.direction
        elif path is not None:
            path[t] = ep.value[0]
            dir_path[t] = ep.direction[0]
    return BlockResult(flag_time, direction, vflag, size, ep.value, path, dir_path)


def fdp_summary(flag_time: np.ndarray, null_until: np.ndarray, n: int) -> dict:
    """False discovery summaries for sticky discovery sets.

    ``null_until[i, k]`` is the change time of asset k on panel i (the first step at which
    its null fails; use n or larger for an asset that never changes). A discovery at step
    t is false when ``t < null_until``. Returns per-panel arrays: ``V``, ``R`` at the end,
    ``fdp_end``, ``sup_fdp`` (maximum over time of the FDP of the discovery set) and
    ``fdp_first`` (FDP at the first step with a discovery, NaN if none).
    """
    ft = np.asarray(flag_time)
    nu = np.asarray(null_until)
    reps, K = ft.shape
    flagged = ft >= 0
    false = flagged & (ft < nu)
    V = false.sum(1)
    R = flagged.sum(1)
    fdp_end = np.where(R > 0, V / np.maximum(R, 1), 0.0)
    sup_fdp = np.zeros(reps)
    fdp_first = np.full(reps, np.nan)
    order_t = np.where(flagged, ft, np.iinfo(np.int64).max)
    idx = np.argsort(order_t, axis=1, kind="stable")
    ts = np.take_along_axis(order_t, idx, 1)
    fs = np.take_along_axis(false, idx, 1)
    for i in np.nonzero(R > 0)[0]:
        t_i, f_i = ts[i, : R[i]], fs[i, : R[i]]
        cum_r = np.arange(1, R[i] + 1)
        cum_v = np.cumsum(f_i)
        last = np.r_[t_i[1:] != t_i[:-1], True]  # FDP is read after all discoveries of the same step
        fdp = cum_v[last] / cum_r[last]
        sup_fdp[i] = fdp.max()
        fdp_first[i] = fdp[0]
    return dict(V=V, R=R, fdp_end=fdp_end, sup_fdp=sup_fdp, fdp_first=fdp_first)


# ----------------------------------------------------------------------------
# Calibration of a panel
# ----------------------------------------------------------------------------


@dataclass
class PanelBaseline:
    """Per-asset baselines for one block (arrays of length K)."""

    active: np.ndarray
    tau2: np.ndarray
    m_up: np.ndarray
    m_dn: np.ndarray
    scale: np.ndarray  # scaled MAD of the calibration window (robust sd)
    n_cal: np.ndarray


def calibrate_panel(windows, cal_len, c: float = 3.0, null_band: float = 0.0, baseline: str = "plugin",
                    delta: float = 0.01) -> PanelBaseline:
    """Calibrate each asset on its own window of recent finite returns (``detectors.calibrate``).

    ``windows`` is a list of K arrays (NaN are dropped). An asset is inactive when it has
    fewer than ``cal_len`` finite returns (int or one per asset) or its window is
    degenerate (zero scale or zero truncated variance). With ``baseline="cs"`` every asset
    uses error level ``delta / K``.
    """
    K = len(windows)
    cal_len = np.broadcast_to(np.asarray(cal_len, int), (K,))
    out = PanelBaseline(np.zeros(K, bool), np.ones(K), np.ones(K), np.full(K, np.inf), np.zeros(K), np.zeros(K, int))
    for k, w in enumerate(windows):
        r = np.asarray(w, float)
        r = r[np.isfinite(r)]
        out.n_cal[k] = r.size
        if r.size < cal_len[k] or r.size < 4:
            continue
        try:
            if baseline == "cs":
                b = calibrate(r, c=c, null_band=null_band, baseline="cs", delta=delta / K)
            else:
                b = calibrate(r, c=c, null_band=null_band)
        except ValueError:
            continue
        if not (math.isfinite(b.tau) and b.tau > 0 and math.isfinite(b.m0_up) and b.m0_up > 0):
            continue
        out.active[k] = True
        out.tau2[k], out.m_up[k], out.m_dn[k], out.scale[k] = b.tau**2, b.m0_up, b.m0_down, b.scale
    return out


# ----------------------------------------------------------------------------
# Streaming panel monitor
# ----------------------------------------------------------------------------


@dataclass
class PanelDiscovery:
    index: int  # step of the panel clock at which the asset was discovered
    asset: int
    name: str
    direction: str  # "up" or "down"
    block: int  # block number (0, 1, ...)
    block_start: int  # first step of the block
    evalue: float  # the asset's e-value at discovery
    threshold: float  # level it had to clear: K / (alpha R) for e-BH, K / alpha or 1 / alpha otherwise
    set_size: int  # number of discoveries in this block so far, including this one
    rule: str


@dataclass
class _BlockInfo:
    start: int
    active: np.ndarray
    tau2: np.ndarray
    m_up: np.ndarray
    m_dn: np.ndarray
    scale: np.ndarray
    n_cal: np.ndarray
    extra: dict = field(default_factory=dict)


class PanelMonitor:
    """Streaming panel monitor: feed one vector of K returns per step, get discoveries.

    Parameters
    ----------
    names : list of K asset names (or an int K).
    alpha : target level (false discovery rate per block for ``rule="ebh"``,
        family-wise error per block for ``"bonferroni"``, per asset for ``"naive"``).
    rule : one of ``RULES``.
    block_len : None or int. If an int, a new block starts automatically every
        ``block_len`` steps, the first one once ``warmup`` steps have been seen
        (default: the largest ``cal_len``). If None, blocks start only when
        ``update(..., block_start=True)`` is called (for example at each UTC
        midnight), and a single block that is never restarted gives the
        "no change since the start" monitor.
    horizon : h of the weights, in each asset's own observations (int or one per
        asset). Defaults to ``block_len``. With ``weights="window"`` evidence from
        changes after the first h observations of a block is no longer added.
    weights : ``"window"`` or ``"pareto"`` (see ``change_weights``).
    cal_len : number of the asset's own most recent finite returns used to
        calibrate it at a block start (int or one per asset). An asset with fewer
        is inactive for that block (its e-value is 1 and it cannot be discovered).
    c, null_band, baseline, delta : as in ``detectors.calibrate``. With
        ``baseline="cs"`` each asset gets error level ``delta / K`` for its
        confidence bounds, so all K baselines hold together with probability at
        least ``1 - delta`` and the false discovery rate bound becomes
        ``alpha + delta`` (``docs/theory_panel.md``, Remark 4).
    record : keep the stopped e-values and directions of every step
        (``self.path``, ``self.path_dir``) and the live, unstopped e-values
        (``self.path_live``).
    """

    def __init__(self, names, alpha: float = 0.05, rule: str = "ebh", block_len: int | None = None,
                 horizon=None, weights: str = "window", cal_len=1440, c: float = 3.0, null_band: float = 0.0,
                 baseline: str = "plugin", delta: float = 0.01, up_grid=DEFAULT_UP_GRID,
                 down_fractions=DEFAULT_DOWN_FRACTIONS, warmup: int | None = None, record: bool = False):
        if isinstance(names, (int, np.integer)):
            names = [f"asset{k}" for k in range(int(names))]
        self.names = list(names)
        self.K = K = len(self.names)
        if K < 1:
            raise ValueError("need at least one asset")
        if not 0 < alpha < 1:
            raise ValueError("alpha must lie in (0, 1)")
        if rule not in RULES:
            raise ValueError(f"rule must be one of {RULES}")
        if weights not in WEIGHTS:
            raise ValueError(f"weights must be one of {WEIGHTS}")
        if baseline not in ("plugin", "cs"):
            raise ValueError("baseline must be 'plugin' or 'cs'")
        self.alpha, self.rule, self.weights = float(alpha), rule, weights
        self.block_len = block_len
        hz = block_len if horizon is None else horizon
        if hz is None:
            raise ValueError("give horizon (or block_len)")
        self.horizon = np.broadcast_to(np.asarray(hz, float), (K,)).copy()
        self.cal_len = np.broadcast_to(np.asarray(cal_len, int), (K,)).copy()
        if np.any(self.cal_len < 4):
            raise ValueError("cal_len must be at least 4")
        self.c, self.band, self.baseline, self.delta = c, null_band, baseline, delta
        self.up_grid, self.down_fractions = up_grid, down_fractions
        self.warmup = int(self.cal_len.max()) if warmup is None else int(warmup)
        self.buffers = [deque(maxlen=int(n)) for n in self.cal_len]
        self.t = -1
        self.block = -1
        self.ep: PanelEProcess | None = None
        self.blocks: list[_BlockInfo] = []
        self.discoveries: list[PanelDiscovery] = []
        self.record = record
        self._path: list[np.ndarray] = []
        self._path_live: list[np.ndarray] = []
        self._path_dir: list[np.ndarray] = []

    # ------------------------------------------------------------------
    def start_block(self) -> None:
        """Recalibrate every asset on its own recent returns and start fresh e-processes."""
        b = calibrate_panel([np.asarray(q, float) for q in self.buffers], self.cal_len, self.c, self.band,
                            self.baseline, self.delta)
        self.block += 1
        self.ep = PanelEProcess(b.tau2[None], b.m_up[None], b.m_dn[None], self.horizon[None], self.weights,
                                active=b.active[None], up_grid=self.up_grid, down_fractions=self.down_fractions)
        self.blocks.append(_BlockInfo(self.t + 1, b.active, b.tau2, b.m_up, b.m_dn, b.scale, b.n_cal))

    def _threshold(self, values: np.ndarray) -> float:
        if self.rule == "ebh":
            R = int(ebh(values, self.alpha).sum())
            return self.K / (self.alpha * R) if R else math.inf
        if self.rule == "bonferroni":
            return self.K / self.alpha
        if self.rule == "naive":
            return 1.0 / self.alpha
        return math.inf

    def update(self, r, block_start: bool = False) -> list[PanelDiscovery]:
        """Feed the K returns of the next step (NaN = no observation). Returns new discoveries.

        ``block_start=True`` starts a new block BEFORE this step's returns are used.
        """
        r = np.asarray(r, float).reshape(-1)
        if r.size != self.K:
            raise ValueError(f"expected {self.K} returns")
        auto = (self.block_len is not None and self.t + 1 >= self.warmup
                and (self.ep is None or (self.t + 1 - self.blocks[-1].start) >= self.block_len))
        if block_start or auto:
            self.start_block()
        self.t += 1
        out: list[PanelDiscovery] = []
        if self.ep is not None:
            ep = self.ep
            ep.step(r[None])
            v = ep.value
            new = flags(v, ep.frozen, self.alpha, self.rule)[0]
            if new.any():
                thr = self._threshold(v[0])
                dirs = ep.direction[0]
                ep.freeze(new[None])
                size = int(ep.frozen.sum())
                for k in np.nonzero(new)[0]:
                    out.append(PanelDiscovery(self.t, int(k), self.names[k], "up" if dirs[k] > 0 else "down",
                                              self.block, self.blocks[-1].start, float(v[0, k]), float(thr),
                                              size, self.rule))
                self.discoveries.extend(out)
            if self.record:
                self._path.append(ep.value[0].copy())
                self._path_live.append(ep.live_value[0].copy())
                self._path_dir.append(ep.direction[0].copy())
        elif self.record:
            self._path.append(np.full(self.K, np.nan))
            self._path_live.append(np.full(self.K, np.nan))
            self._path_dir.append(np.zeros(self.K, int))
        for k in range(self.K):  # buffers get this step's returns only after it was monitored
            if np.isfinite(r[k]):
                self.buffers[k].append(float(r[k]))
        return out

    def run(self, returns, block_starts=None) -> list[PanelDiscovery]:
        """Feed a (n, K) array. ``block_starts``: boolean array of length n (or indices) marking block starts."""
        R = np.asarray(returns, float)
        if R.ndim != 2 or R.shape[1] != self.K:
            raise ValueError("returns must have shape (n, K)")
        bs = np.zeros(R.shape[0], bool)
        if block_starts is not None:
            b = np.asarray(block_starts)
            if b.dtype == bool:
                bs[:] = b
            else:
                bs[b] = True
        for i in range(R.shape[0]):
            self.update(R[i], block_start=bool(bs[i]))
        return self.discoveries

    @property
    def path(self) -> np.ndarray:
        return np.asarray(self._path)

    @property
    def path_live(self) -> np.ndarray:
        return np.asarray(self._path_live)

    @property
    def path_dir(self) -> np.ndarray:
        return np.asarray(self._path_dir)
