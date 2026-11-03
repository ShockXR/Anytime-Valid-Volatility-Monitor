"""Every script in examples/ runs to completion in a fresh interpreter (synthetic data only)."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = sorted((ROOT / "examples").glob("*.py"))


def test_examples_exist():
    assert len(EXAMPLES) >= 3


@pytest.mark.parametrize("script", EXAMPLES, ids=lambda p: p.name)
def test_example_runs(script, tmp_path):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(ROOT / "src"), env.get("PYTHONPATH")]))
    env["MPLBACKEND"] = "Agg"
    proc = subprocess.run([sys.executable, str(script)], cwd=tmp_path, env=env, capture_output=True, text=True,
                          timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert "[SYNTHETIC]" in proc.stdout
    assert "Traceback" not in proc.stderr
