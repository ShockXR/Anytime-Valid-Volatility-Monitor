"""The top-level API: every exported name resolves, old import paths still work, the version is single-sourced."""

import re
import subprocess
import sys
from pathlib import Path

import pytest

import avmonitor

ROOT = Path(__file__).resolve().parents[1]


def test_all_names_resolve():
    assert len(avmonitor.__all__) == len(set(avmonitor.__all__))
    for name in avmonitor.__all__:
        assert hasattr(avmonitor, name), name


def test_top_level_names_are_the_module_objects():
    from avmonitor import confseq, detectors, monitor, multifreq, panel, periodic, robust

    assert avmonitor.VolatilityMonitor is monitor.VolatilityMonitor
    assert avmonitor.PeriodicVolatilityMonitor is periodic.PeriodicVolatilityMonitor
    assert avmonitor.MultiFrequencyMonitor is multifreq.MultiFrequencyMonitor
    assert avmonitor.PanelMonitor is panel.PanelMonitor
    assert avmonitor.betting_edetector is detectors.betting_edetector
    assert avmonitor.cs_baseline is confseq.cs_baseline
    assert avmonitor.kernel_iq is robust.kernel_iq


def test_submodules_still_importable_as_modules():
    import types

    from avmonitor import detectors, robust

    assert isinstance(robust, types.ModuleType) and isinstance(detectors, types.ModuleType)


def test_version_is_single_sourced():
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert re.search(r'^dynamic\s*=\s*\["version"\]', text, flags=re.M)
    assert not re.search(r'^version\s*=\s*"', text, flags=re.M)  # no static version in [project]
    assert 'attr = "avmonitor.__version__"' in text
    assert re.fullmatch(r"\d+\.\d+\.\d+", avmonitor.__version__)


def test_api_docs_builder_runs(tmp_path):
    script = ROOT / "scripts" / "build_api_docs.py"
    if not script.exists():  # a source distribution ships without scripts/
        pytest.skip("scripts/build_api_docs.py is not in this tree")
    out = tmp_path / "api.md"
    proc = subprocess.run([sys.executable, str(script), "-o", str(out)],
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr
    text = out.read_text(encoding="utf-8")
    for name in avmonitor.__all__:
        assert f"`{name}`" in text, name
