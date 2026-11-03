"""scripts/trace_numbers.py: extraction rules, the rounding rule and a small end-to-end trace."""

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("trace_numbers", ROOT / "scripts" / "trace_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["trace_numbers"] = mod
    spec.loader.exec_module(mod)
    return mod


tn = _load()


def values(line, **kw):
    return [(t.mantissa, t.decimals, t.exponent, t.percent) for t in tn.extract_line(line, 1, **kw)]


def texts(line, **kw):
    return [t.text for t in tn.extract_line(line, 1, **kw)]


@pytest.mark.parametrize("line, expected", [
    ("with 1,940 returns and $A = 28{,}800$", [(1940, 0, 0, False), (28800, 0, 0, False)]),
    ("rates of 93.6% and 5.0\\% at p = 0.010", [(93.6, 1, 0, True), (5.0, 1, 0, True), (0.010, 3, 0, False)]),
    ("slopes $-0.48$ to $-0.68$ and a bias of $-17\\%$", [(-0.48, 2, 0, False), (-0.68, 2, 0, False),
                                                         (-17, 0, 0, True)]),
    ("74 to 78% of alarms", [(74, 0, 0, True), (78, 0, 0, True)]),
    ("ratios 4.00 to $4.03 \\times 10^{-5}$", [(4.00, 2, -5, False), (4.03, 2, -5, False)]),
    ("at $A = 10^8$ and from 20 million draws", [(1, 0, 8, False), (20, 0, 6, False)]),
    ("the 95th percentile and a 1440-minute day", [(95, 0, 0, False), (1440, 0, 0, False)]),
])
def test_extraction_spellings(line, expected):
    got = values(line)
    assert len(got) == len(expected)
    for (m, d, e, p), (m2, d2, e2, p2) in zip(got, expected):
        assert m == pytest.approx(m2) and (d, e, p) == (d2, e2, p2)


@pytest.mark.parametrize("line", [
    "Following Bapat (2026) and Page 1954; see Section 5.6 and Sections 5 and 6.",
    "Proposition A.6, Lemma A.7, Theorem 6.4, eq. (12), Table 7, Figure 2 and Supplement S4.1.",
    "Labels E14, Q14h and ALT1; GARCH(1,1) and N(0,1) returns; Student $t_5$.",
    "From 1 July to 31 December 2025 and 21 January 2024; statements at 14:00 on 2026-11-03.",
    "Run `python scripts/reproduce.py --quick` or read `results/e14d_real_alarms.csv`.",
    "The Nifty 50 index; $r_t^2$, $\\kappa_2$ and $\\tfrac13$; the 16 FOMC statements of 2024 and 2025.",
    "See [the report](docs/report_parts/E12.md) and ![fig](../results/figures/e09_stress.png).",
])
def test_non_results_are_skipped(line):
    assert [t for t in texts(line) if t != "16"] == []


def test_years_and_counts():
    assert texts("within 2000 steps and a 2000-step window, 2000 paths") == ["2000", "2000", "2000"]
    assert texts("| 1950 | 25% |") == ["1950", "25%"]
    assert texts("development 2017 to 2021, held-out 2022") == []


def test_small_integers_flag():
    assert texts("on 9 of 11 days and 3 times") == ["11"]
    assert texts("on 9 of 11 days and 3 times", keep_small_ints=True) == ["9", "11", "3"]
    assert texts("fell by 1 to 4 points and from 2 to 1.25") == ["1", "4", "2", "1.25"]


def test_document_skips_front_matter_references_and_code():
    doc = "\n".join([
        "---", 'title: "x"', 'date: "Draft 2026, 12 items"', "---", "",
        "We saw 26.5 alarms.", "", "```bash", "python run.py --n 3000   # 3000 steps", "```", "",
        "## References", "", "Bapat, S. (2026). Journal 12, 105482.", "", "## Appendix", "", "And 33.3."])
    assert [t.text for t in tn.extract_numbers(doc)] == ["26.5", "33.3"]


@pytest.mark.parametrize("mant, dec, exp, percent, data, source, expected", [
    (0.94, 2, 0, False, False, 0.9363, "rounded"),
    (93.6, 1, 0, True, False, 0.9363, "rounded, source x 100"),
    (94, 0, 0, True, False, 0.9363, "rounded, source x 100"),
    (0.93, 2, 0, False, False, 0.9363, None),
    (93.6, 1, 0, False, False, 0.9363, None),          # x 100 only for percentages
    (0.94, 2, 0, False, True, 93.63, "rounded, source / 100"),
    (0.94, 2, 0, False, False, 93.63, None),           # / 100 only against data files
    (3.99, 2, -5, False, False, 3.986e-5, "rounded"),
    (1791, 0, 0, False, False, 1790.6, "rounded"),
])
def test_rounding_rule(mant, dec, exp, percent, data, source, expected):
    arr = np.array([source])
    assert tn.value_match(arr, mant, dec, exp, percent=percent, data=data) == expected


def test_end_to_end(tmp_path):
    (tmp_path / "docs" / "report_parts").mkdir(parents=True)
    (tmp_path / "results").mkdir()
    (tmp_path / "docs" / "report_parts" / "E4.md").write_text(
        "The true baseline gave 26.5 alarms against a bound of 33.3.\n\n"
        "| law | a | b | c | d |\n|---|---|---|---|---|\n| Gaussian | 1.11 | 2.22 | 3.33 | 4.44 |\n")
    (tmp_path / "docs" / "report_parts" / "E9.md").write_text("Unrelated numbers 77.7 and 12.5.\n")
    (tmp_path / "results" / "e04a_null.csv").write_text("rate,paths\n0.9363,300\n")
    target = tmp_path / "paper.md"
    target.write_text(
        "With the true baseline it gave 26.5 alarms against 33.3 (E4). The share was 93.6% (E4). "
        "A stray 77.7% (E4). A missing 61.3% (E4).\n\n"
        "| law | a | b | c | d |\n|---|---|---|---|---|\n| Gaussian | 1.11 | 2.22 | 3.33 | 4.45 |\n\n"
        "*Table 1. Delays (E4).*\n")
    out = tmp_path / "trace.csv"
    assert tn.main(["--root", str(tmp_path), "--target", str(target), "--out", str(out), "--quiet"]) == 0
    rows = {r["number"]: r for r in csv.DictReader(out.open())}
    assert rows["26.5"]["status"] == "traced_md"
    assert rows["26.5"]["source"] == "docs/report_parts/E4.md"
    assert rows["26.5"]["cooccur"] == "docs/report_parts/E4.md"
    assert rows["93.6%"]["status"] == "traced_csv"
    assert rows["93.6%"]["source"] == "results/e04a_null.csv"
    assert (rows["77.7%"]["status"], rows["77.7%"]["in_scope"]) == ("traced_md", "no")  # E9 only, not weak
    assert rows["61.3%"]["status"] == "untraced"
    assert rows["4.45"]["row_check"].startswith("not in docs/report_parts/E4.md:5")
    assert rows["1.11"]["row_check"] == ""
    assert tn.main(["--root", str(tmp_path), "--target", str(target), "--out", str(out), "--quiet",
                    "--fail-on-untraced"]) == 1
