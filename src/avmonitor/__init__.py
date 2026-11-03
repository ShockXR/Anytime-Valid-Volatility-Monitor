"""Anytime-valid online detection of volatility regime changes.

The main entry points are importable from the top level, for example
``from avmonitor import VolatilityMonitor``. Every name below is also available
from its own module (``avmonitor.monitor``, ``avmonitor.detectors`` and so on),
so imports written for earlier versions keep working. The reference for every
module is ``docs/api.md``.
"""

__version__ = "0.3.0"

from .confseq import ConfSeq, CSBaseline, StreamingEBConfSeq, cs_baseline, eb_confseq
from .detectors import (
    DEFAULT_DOWN_FRACTIONS,
    DEFAULT_UP_GRID,
    Baseline,
    DetectorResult,
    betting_edetector,
    calibrate,
    gaussian_lr_edetector,
    page_cusum,
    repeated_window_ztest,
    summarize_alarms,
    two_sided_alarm,
)
from .monitor import Event, GaussianLRMonitor, VolatilityMonitor
from .multifreq import MFEvent, MultiFrequencyMonitor, level_for_budget
from .panel import PanelDiscovery, PanelMonitor, ebh
from .periodic import PeriodicVolatilityMonitor, restarted_sr_bound
from .robust import (
    MAD_NORMAL_SCALE,
    bipower_quarticity,
    bipower_variation,
    ew_quantile_scale,
    kernel_iq,
    local_quarticity,
    local_variance,
    mad,
    mad_threshold,
    med_rq,
    min_rq,
    quadpower_quarticity,
    realized_quarticity,
    threshold_rq,
    tripower_quarticity,
    truncate,
    truncated_rq_mad,
)

__all__ = [
    # streaming monitors
    "VolatilityMonitor",
    "GaussianLRMonitor",
    "Event",
    "PeriodicVolatilityMonitor",
    "restarted_sr_bound",
    "MultiFrequencyMonitor",
    "MFEvent",
    "level_for_budget",
    "PanelMonitor",
    "PanelDiscovery",
    "ebh",
    # detectors on arrays
    "calibrate",
    "Baseline",
    "betting_edetector",
    "two_sided_alarm",
    "DetectorResult",
    "gaussian_lr_edetector",
    "page_cusum",
    "repeated_window_ztest",
    "summarize_alarms",
    "DEFAULT_UP_GRID",
    "DEFAULT_DOWN_FRACTIONS",
    # confidence-sequence baseline
    "cs_baseline",
    "CSBaseline",
    "eb_confseq",
    "ConfSeq",
    "StreamingEBConfSeq",
    # robust estimators
    "MAD_NORMAL_SCALE",
    "mad",
    "mad_threshold",
    "truncate",
    "ew_quantile_scale",
    "realized_quarticity",
    "truncated_rq_mad",
    "kernel_iq",
    "local_variance",
    "local_quarticity",
    "min_rq",
    "med_rq",
    "tripower_quarticity",
    "quadpower_quarticity",
    "bipower_quarticity",
    "bipower_variation",
    "threshold_rq",
]
