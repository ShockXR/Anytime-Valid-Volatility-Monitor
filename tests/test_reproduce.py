"""scripts/reproduce.py: order of the E13 and E14 entries, the subset options of the run step, and the comparison."""

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("reproduce", ROOT / "scripts" / "reproduce.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["reproduce"] = mod
    spec.loader.exec_module(mod)
    return mod


R = _load()


def _commands(prefix):
    return [e[1][1:] for e in R.EXPERIMENTS if e[0].startswith(prefix)]


def test_e13_round3_parts_run_after_the_first_parts_and_figures_last():
    cmds = [" ".join(c) for c in _commands("E13")]
    first = ["offline-null", "offline-change", "online-null", "online-change 0", "online-change 1", "real-calib",
             "real-crypto", "real-crypto-h1-2025", "real-nifty", "cost"]
    round3 = ["offline-long", "online-long", "fix-nominal", "cap-check", "real-long nifty", "real-long crypto",
              "long-pooled"]
    assert cmds == first + round3 + ["figures"]


def test_e14_parts_a_and_b_write_the_committed_split_files():
    cmds = _commands("E14")
    assert ["a", "--cal=390"] in cmds and ["a", "--cal=1440"] in cmds
    assert {c[1] for c in cmds if c[0] == "b"} == {"--law=gaussian", "--law=student_t3", "--law=jumps"}
    assert all("a" not in c and "b" not in c for c in cmds if c[0] not in ("a", "b"))  # no unsuffixed a or b run
    assert cmds[-1][-1] == "figures"
    for name in ("e14a_validity_cal390", "e14a_validity_cal1440", "e14b_delay_gaussian", "e14b_delay_student_t3",
                 "e14b_delay_jumps"):
        assert (ROOT / "results" / f"{name}.csv").exists()


def test_run_subset_options(monkeypatch):
    ran = []

    class Done:
        returncode = 0

    monkeypatch.setattr(R.subprocess, "run", lambda cmd, cwd=None: ran.append(cmd[1:]) or Done())
    assert R.step_run("full", True, start="E13", stop="E14")
    labels = [e[0] for e in R.EXPERIMENTS]
    expected = [e[1] for e in R.EXPERIMENTS if e[0].startswith(("E13", "E14"))]
    assert ran == expected
    assert labels.index("E15 panel, synthetic") > max(i for i, s in enumerate(labels) if s.startswith("E14"))
    ran.clear()
    assert R.step_run("full", True, only="E9 stress")
    assert ran == [["experiments/e09_stress_guarantee.py"]]
    ran.clear()
    assert not R.step_run("full", True, stop="no such label")
    assert ran == []


def test_numeric_equal_ignores_only_the_listed_columns(tmp_path):
    committed = 'task,steps,seconds,note\n"a, b",10,1.5,null\nc,20,2.5,x\n'
    rerun = tmp_path / "f.csv"
    rerun.write_text('task,steps,seconds,note\n"a, b",10,1.7,null\nc,20,2.4,x\n')
    assert not R._numeric_equal(rerun, committed, 1e-9)
    assert R._numeric_equal(rerun, committed, 1e-9, ignore=("seconds",))  # "null" is read as missing on both sides
    rerun.write_text('task,steps,seconds,note\n"a, b",11,1.7,null\nc,20,2.4,x\n')
    assert not R._numeric_equal(rerun, committed, 1e-9, ignore=("seconds",))
    rerun.write_text('task,steps,seconds,note\n"a, b",10,1.7,y\nc,20,2.4,x\n')
    assert not R._numeric_equal(rerun, committed, 1e-9, ignore=("seconds",))


def test_timing_columns_exist_in_the_committed_files():
    import pandas as pd

    for f, cols in R.TIMING_COLUMNS.items():
        header = list(pd.read_csv(ROOT / f, nrows=0).columns)
        assert set(cols) <= set(header), f
