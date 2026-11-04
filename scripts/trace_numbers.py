"""Trace every number in a Markdown document to a committed source.

Run from the repository root:

    python scripts/trace_numbers.py                        # paper/short_paper.md
    python scripts/trace_numbers.py --target README.md
    python scripts/trace_numbers.py --target PITCH_BAPAT.md --keep-small-ints

The script writes one row per number to ``results/trace/<name>_trace.csv``
(``paper_trace.csv`` for the paper, ``readme_trace.csv`` for README.md) and
prints a summary and the list of untraced numbers. It reads files only; it
never edits the target.

Extraction rules
----------------
A number is a run of ASCII digits, with these spellings understood:

* decimals ``0.94``; the number of decimals shown is kept (``0.010`` has 3);
* percentages ``93.6%`` and LaTeX ``93.6\\%``; ``points`` after a number also
  marks it as a percentage (only used for the small-integer rule below);
* negative numbers ``-0.48``, ``$-17\\%$``, ``\u22120.48`` (Unicode minus), when the
  sign follows a space, an opening bracket, ``$``, ``=`` or the start of the
  line (``B-1``, ``(1-b)`` are subtractions, not negative numbers);
* thousands separators ``1,940``, ``264,959`` and LaTeX ``28{,}800``: a comma
  followed by exactly three digits (a comma and a space separate list items);
* scientific notation ``3.99 \\times 10^{-5}``, ``2 \u00d7 10\u207b\u2076``, ``1e-6``, bare
  powers ``10^8`` and ``20 million``; the precision is that of the mantissa;
* ordinals ``95th`` (the number 95);
* ranges ``a to b``: each end is a number. When only the right end carries a
  percent sign or a power of ten (``74 to 78%``, ``4.00 to 4.03 \\times 10^{-5}``)
  the left end inherits it;
* numbers inside inline and display LaTeX math are read with the same rules.

Skipped, because they are not results:

1. YAML front matter and the References or Bibliography section.
2. Link and image targets ``](...)``, bare URLs, ``doi:`` strings, HTML comments.
3. Fenced code blocks (examples, not results), and inline code spans that look
   like a path or a command: they contain ``/`` or a file extension (``.py``,
   ``.csv``, ...), start with a command word (``python``, ``pytest``, ``git``,
   ...) or contain an option such as ``--quick``. Other code spans are read
   like text.
4. Labels and identifiers: a word that starts with a letter and contains a
   digit (``E14``, ``Q14h``, ``ALT1``, ``S4.1``, ``t5``, ``H1``, ``x86_64``,
   ``e09_stress``); theorem-style labels such as ``A.6`` or ``B.2``; model
   orders and Landau symbols written as a word followed by digits in brackets
   (``GARCH(1,1)``, ``N(0,1)``, ``O(1)``); index names (``Nifty 50``, ``S&P 500``).
5. Cross-references: numbers right after Section, Table, Figure, eq.,
   equation, Theorem, Lemma, Proposition, Corollary, Definition, Remark,
   Assumption, Appendix, Supplement, Step, item or round (plural forms take lists,
   as in "Sections 5 and 6"), and the leading section number of a heading.
6. LaTeX indices and fixed fractions: sub- and superscripts (``r_t^2``,
   ``\\kappa_2``, ``M_t^{(k)}``) and ``\\frac12``, ``\\tfrac13``.
7. Dates and times: ISO dates (``2026-11-03``), clock times (``14:00``), day
   numbers next to a month name (``1 July``, ``July 1``), years after a month,
   a half or a quarter (``December 2025``, ``H1 2025``), and every other
   integer from 1900 to 2099 used as a year: in a citation (``Bapat (2026)``,
   ``Page 1954``), or continuing a list or range of years (``2024 and 2025``,
   ``2017 to 2021``). Such an integer counts as a number, not a year, when a
   hyphen or a lower-case word other than and/or/to follows it (``2000-step``,
   ``2000 paths``, ``1950 steps``), when it fills a table cell on its own, or
   when it follows ``=``, ``<``, ``>``, ``\\le`` or ``\\ge``.
8. Small integers 0 to 10 that stand alone: no decimals, no percent sign, no
   power of ten, and not the end of a range whose other end is kept. Use
   ``--keep-small-ints`` to keep them.

Sources
-------
* Markdown: every ``docs/**/*.md`` except the target itself and
  ``docs/experiment_report.md`` (it is assembled from ``docs/report_parts``, so it
  would trace everything trivially). Nothing under ``paper/`` is a source.
  Markdown sources come in three tiers: *primary* (``docs/report_parts``,
  ``docs/theory*.md``, ``docs/math.md``, ``docs/preregistration*.md``), *other*
  (the rest of ``docs/``, e.g. ``docs/iterations``) and *review*
  (``docs/review/`` and ``docs/paper_review.md``, which may quote the paper).
* Data: every CSV (also ``.csv.gz``), JSON and plain-text file under ``results/``,
  including subfolders, except ``results/trace/`` (this script's own output).
  Numeric cells, JSON numbers and numeric strings, and numbers in text files
  are read as values.

Scope
-----
Experiment labels are read from the sentence that holds the number (the
narrow scope); when it has none, from its paragraph or table (for a table also
its caption and the paragraph before it); when that has none, from the whole
section. ``E4``, ``E4e`` and ``Q4b`` all mean E4; ``ALT1`` means ALT1. The wide
scope adds every label of the paragraph and the section. In-scope data files
are ``results/e04*`` for E4 and so on; in-scope Markdown files are
``docs/report_parts/E4.md`` (or ``E01_E05.md``), ``docs/iterations/E4.md``,
``docs/report_parts/ALT1.md`` and the theory note of the experiment
(``theory_running_tau`` for E14, ``theory_panel`` for E15, ``theory_multifreq``
for E17, ``theory_delay`` for E11). A number whose section has no label is
matched against all sources; ``--no-scope`` does that for every number.

Tracing
-------
Each number is tried in this order; the first hit decides ``status`` and
``source``.

1. ``traced_md``, literal, narrow scope: the number as written appears as a
   number in an in-scope Markdown source (sign, digits and decimals must agree;
   thousands separators, ``{,}``, Unicode minus, ``minus`` spelled out and the
   percent sign are not compared; scientific notation is compared as
   ``<mantissa>e<exponent>``). Primary tier first, review tier last.
2. ``traced_md``, literal, wide scope (``in_scope`` = paragraph).
3. ``traced_md``, rounded: a number in an in-scope Markdown source (primary or
   other tier) satisfies the rounding rule.
4. ``traced_csv``: a value in an in-scope data file satisfies the rounding rule.
5. ``traced_md``, literal, any Markdown source out of scope (``in_scope`` = no),
   only when the match is not weak (see below).

Rounding rule: a source value v traces a number p shown with d decimals if v
lies within half a unit of the last shown digit of p (for scientific notation,
of the mantissa); for a percentage also 100 v (a share stored as a fraction),
and for a non-percentage between 0 and 1 matched against a data file also
v / 100 (a share stored in percent). So 0.9363 traces 0.94, 93.6% and 94%, and
does not trace 0.93. Rounding matches are only accepted in the narrow scope,
because 25 MB of results contain a value near almost any short number.

How much a trace says
---------------------
Because the sources are dense in numbers, a trace on its own is weak evidence:
when every number of the paper is moved by three units of its last digit, about
90% of the wrong values still trace. Three extra columns help:

* ``weak``: four wrong values (the last digit moved by -7, -3, +3, +7) are tried
  against the file that traced the number; ``yes`` when two or more match too.
* ``cooccur``: the first in-scope Markdown file (narrow, then wide scope) in
  which one sentence or table row holds the number together with another
  number of its own sentence or table row; ``single`` when the sentence has no
  other number; ``common`` when the number appears in at least five sentences
  of the in-scope reports (a setting such as 390 steps); ``no`` otherwise. In
  the same perturbation test 34% to 40% of wrong values co-occur, against 94%
  of the paper's own numbers. Data rows are not used for this: a CSV
  with thousands of rows holds two short numbers in one row by chance too often.
* ``row_check``: each table row of the target with at least three distinct
  numbers is compared with the in-scope Markdown table row that holds most of
  them; when that row holds at least 75% of them, every number of the target
  row missing from it is reported with the file and line of the source row.
  It catches about half of single wrong cells and flags rows that merge
  several source rows only rarely.

Output columns: ``line``, ``number`` (as written), ``context``, ``status``,
``source`` (first source file that traced it), then ``how``, ``tier``,
``in_scope``, ``weak``, ``cooccur``, ``row_check``, ``scope`` (narrow labels)
and ``hint`` (for untraced numbers: a file under ``alternatives/*/results``,
which is not a declared source, or a data or Markdown file outside the scope
that would match, as a lead for checking by hand).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

MONTHS = (r"(?:January|February|March|April|May|June|July|August|September|October|November|December|"
          r"Jan|Feb|Mar|Apr|Jun|Jul|Aug|Sept|Sep|Oct|Nov|Dec)")
SUPERSCRIPTS = str.maketrans("\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b\u207a",
                             "0123456789-+")
SUP_CLASS = "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079"
MINUS = "\u2212"

CODE_PATH_OR_COMMAND = re.compile(
    r"/|\\\\|\.(?:py|md|csv|json|png|pdf|txt|gz|tex|bib|toml|ya?ml|ipynb|sh|cfg|ini|html|parquet|pkl)\b"
    r"|^\s*(?:python3?|pip|pytest|git|make|pandoc|bash|sh|cd|ls|cat|grep|curl|wget|conda|uv|avmonitor|jupyter)\b"
    r"|(?:^|\s)--?[A-Za-z]")

SINGULAR_REF = r"Section|section|Table|table|Figure|figure|Fig\.|fig\.|[Ee]q\.|[Ee]quation|Theorem|Lemma|" \
               r"Proposition|Corollary|Definition|Remark|Assumption|Appendix|Supplement|Step|[Ii]tem|[Rr]ound|Part|Chapter|§"
PLURAL_REF = r"Sections|sections|Tables|tables|Figures|figures|Figs\.|[Ee]qs\.|[Ee]quations|Theorems|Lemmas|" \
             r"Propositions|Corollaries|Definitions|Remarks|Assumptions|Appendices|Supplements|Steps|[Ii]tems|[Rr]ounds|Parts"
REF_ITEM = r"(?:[A-Z]?[0-9]+(?:\.[0-9]+)*[a-z]?|[A-Z](?:\.[0-9]+)*)(?![0-9A-Za-z])"
REF_SINGLE_RE = re.compile(rf"(?<![A-Za-z])(?:{SINGULAR_REF})\s*\(?\s*({REF_ITEM})")
REF_LIST_RE = re.compile(
    rf"(?<![A-Za-z])(?:{PLURAL_REF})\s*\(?\s*({REF_ITEM}(?:\s*(?:,|and|or|to|&|-)\s*{REF_ITEM})*)")

# Scientific notation (found first, then blanked out of the working line).
SCI_TIMES_RE = re.compile(
    rf"(?<![A-Za-z0-9_.])(?P<sign>[-+{MINUS}])?(?P<mant>[0-9]+(?:\.[0-9]+)?)\s*(?:\\times|\u00d7|\\cdot|x|\*)\s*"
    rf"10\s*(?:\^\s*\{{\s*(?P<e1>[-+{MINUS}]?[0-9]+)\s*\}}|\^\s*(?P<e2>[-+{MINUS}]?[0-9]+)|(?P<e3>[\u207b\u207a]?[{SUP_CLASS}]+))")
SCI_POWER_RE = re.compile(
    rf"(?<![A-Za-z0-9_.])10\s*(?:\^\s*\{{\s*(?P<e1>[-+{MINUS}]?[0-9]+)\s*\}}|\^\s*(?P<e2>[-+{MINUS}]?[0-9]+)"
    rf"|(?P<e3>[\u207b\u207a]?[{SUP_CLASS}]+))")
SCI_E_RE = re.compile(rf"(?<![A-Za-z0-9_.])(?P<sign>[-+{MINUS}])?(?P<mant>[0-9]+(?:\.[0-9]+)?)[eE](?P<e1>[-+]?[0-9]+)(?![A-Za-z0-9])")

SCI_WORD_RE = re.compile(rf"(?<![A-Za-z0-9_.])(?P<sign>[-+{MINUS}])?(?P<mant>[0-9]+(?:\.[0-9]+)?)\s+(?P<word>million|billion)\b")
WORD_EXP = {"million": 6, "billion": 9}
NUM_RE = re.compile(
    rf"(?<![A-Za-z0-9_.])(?P<sign>[-+{MINUS}])?"
    r"(?P<int>[0-9]{1,3}(?:(?:,|\{,\})[0-9]{3})+(?![0-9])|[0-9]+)"
    r"(?P<frac>\.[0-9]+)?(?![0-9])(?!\.[0-9])"
    r"(?P<ord>st|nd|rd|th)?(?![A-Za-z]*[0-9])"
    r"(?P<pct>\s?\\?%)?")

INDEX_NAME_RE = re.compile(r"\b(?:Nifty|NIFTY|S&P|FTSE|DAX|Stoxx|STOXX|Russell|Nikkei|Sensex)\s+([0-9]+)\b")
LABEL_RE = re.compile(r"(?<![A-Za-z0-9_])(?:(?P<eq>[EQ])(?P<n>[0-9]{1,2})[a-z]?|ALT(?P<alt>[0-9])|TH(?P<th>[0-9]))(?![0-9])")
SENTENCE_BREAK_RE = re.compile(r"(?<=[.!?])\s+(?=[A-Z(\[*$\\])")
YEAR_CONNECTORS = {"and", "or", "to", "through", "until"}


@dataclass
class NumberToken:
    """One number found in the target."""

    line: int                 # 1-based line number
    start: int                # column in the line
    end: int
    text: str                 # as written
    mantissa: float           # value without the power of ten
    decimals: int             # decimals shown in the mantissa
    exponent: int = 0         # power of ten
    percent: bool = False
    has_sep: bool = False     # thousands separator used
    ordinal: bool = False
    sci: bool = False
    context: str = ""
    unit: int = 0             # sentence (or table row) of the line that holds the number

    @property
    def value(self) -> float:
        return self.mantissa * 10.0 ** self.exponent

    def literal_key(self) -> str:
        """Normalised spelling used for literal matching."""
        if self.sci:
            return sci_key(self.mantissa, self.decimals, self.exponent)
        return format_fixed(self.mantissa, self.decimals)


def format_fixed(value: float, decimals: int) -> str:
    s = f"{value:.{decimals}f}"
    return s.lstrip("-") if float(s) == 0 else s


def sci_key(mantissa: float, decimals: int, exponent: int) -> str:
    return f"{mantissa:.{decimals}f}e{exponent}"


def _blank(chars: list[str], start: int, end: int) -> None:
    for k in range(max(start, 0), min(end, len(chars))):
        chars[k] = " "


def _apply(chars: list[str], pattern: re.Pattern, group: int | str = 0, keep=None) -> None:
    """Blank every match (or one group of it) of pattern in the working line."""
    text = "".join(chars)
    for m in pattern.finditer(text):
        if keep is not None and keep(m):
            continue
        if m.group(group) is None:
            continue
        _blank(chars, m.start(group), m.end(group))


def _exp_value(m: re.Match) -> int:
    gd = m.groupdict()
    raw = gd.get("e1") or gd.get("e2")
    if raw is None:
        raw = (gd.get("e3") or "").translate(SUPERSCRIPTS)
    return int(raw.replace(MINUS, "-"))


def _decimals(mant: str) -> int:
    return len(mant.split(".")[1]) if "." in mant else 0


def _sign(sign: str | None, line: str, pos: int) -> int:
    """+1 or -1; a minus right after a letter, digit or closing bracket is a subtraction."""
    if not sign or sign == "+":
        return 1
    before = line[pos - 1] if pos > 0 else " "
    if before.isalnum() or before in ")]}":
        return 1
    return -1


def is_code_path_or_command(code: str) -> bool:
    return bool(CODE_PATH_OR_COMMAND.search(code))


def _context(line: str, start: int, end: int, width: int = 45) -> str:
    a, b = max(0, start - width), min(len(line), end + width)
    snippet = re.sub(r"\s+", " ", line[a:b]).strip()
    return ("..." if a > 0 else "") + snippet + ("..." if b < len(line) else "")


def _skip_regions(lines: list[str]) -> tuple[set[int], set[int]]:
    """Indices of lines to skip entirely, and indices of lines inside fenced code."""
    skip, fenced = set(), set()
    if lines and lines[0].strip() == "---":
        for k in range(1, len(lines)):
            if lines[k].strip() in ("---", "..."):
                skip.update(range(0, k + 1))
                break
    in_fence, fence = False, ""
    ref_level = None
    for k, line in enumerate(lines):
        stripped = line.strip()
        if in_fence:
            if stripped.startswith(fence):
                in_fence = False
                skip.add(k)
            else:
                fenced.add(k)
            continue
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_fence, fence = True, stripped[:3]
            skip.add(k)
            continue
        h = re.match(r"^(#{1,6})\s+(.*)$", line)
        if h:
            level = len(h.group(1))
            if ref_level is not None and level <= ref_level:
                ref_level = None
            if re.match(r"(?i)(references|bibliography|literature cited)\s*$", h.group(2).strip()):
                ref_level = level
                skip.add(k)
                continue
        if ref_level is not None:
            skip.add(k)
    return skip, fenced


def _is_table_cell_alone(line: str, start: int, end: int) -> bool:
    if not line.lstrip().startswith("|"):
        return False
    left = line[:start].rstrip()
    right = line[end:].lstrip()
    return left.endswith("|") and right.startswith("|")


def _year_like(tok: NumberToken) -> bool:
    return (tok.decimals == 0 and not tok.percent and not tok.has_sep and not tok.sci and not tok.ordinal
            and tok.exponent == 0 and 1900 <= tok.mantissa <= 2099)


def _between(line: str, a: NumberToken, b: NumberToken) -> str:
    return line[a.end:b.start]


def _is_range_gap(gap: str) -> bool:
    return bool(re.fullmatch(r"\s*\$?\s*(?:to|and|or|-)\s*\$?\s*", gap))


def extract_line(line: str, lineno: int, *, keep_small_ints: bool = False) -> list[NumberToken]:
    """Extract the result-like numbers of one line of prose (see the module docstring for the rules)."""
    chars = list(line)

    # 2. link targets, URLs, doi strings, HTML comments
    _apply(chars, re.compile(r"\]\(([^)]*)\)"), 1)
    _apply(chars, re.compile(r"https?://[^\s)>\]]+|doi:\S+|<!--.*?-->"))
    # 3. code spans
    text = "".join(chars)
    for m in re.finditer(r"`([^`]*)`", text):
        if is_code_path_or_command(m.group(1)):
            _blank(chars, m.start(), m.end())
        else:
            chars[m.start()] = " "
            chars[m.end() - 1] = " "

    tokens: list[NumberToken] = []

    # scientific notation first, so that its exponent is not taken for a superscript
    text = "".join(chars)
    for pattern in (SCI_TIMES_RE, SCI_E_RE, SCI_WORD_RE):
        for m in pattern.finditer(text):
            mant = m.group("mant")
            sgn = _sign(m.group("sign"), line, m.start())
            exp = WORD_EXP[m.group("word")] if pattern is SCI_WORD_RE else _exp_value(m)
            tokens.append(NumberToken(lineno, m.start(), m.end(), line[m.start():m.end()], sgn * float(mant),
                                      _decimals(mant), exp, sci=True))
            _blank(chars, m.start(), m.end())
        text = "".join(chars)
    for m in SCI_POWER_RE.finditer(text):
        tokens.append(NumberToken(lineno, m.start(), m.end(), line[m.start():m.end()], 1.0, 0, _exp_value(m),
                                  sci=True))
        _blank(chars, m.start(), m.end())

    # 7. dates and times
    _apply(chars, re.compile(r"(?<![0-9])[0-9]{4}-[0-9]{2}(?:-[0-9]{2})?(?![0-9])"))
    _apply(chars, re.compile(r"(?<![0-9.,])[0-9]{1,2}:[0-9]{2}(?::[0-9]{2})?(?![0-9])"))
    _apply(chars, re.compile(rf"(?<![0-9.,])([0-9]{{1,2}})(?:st|nd|rd|th)?\s+(?:of\s+)?{MONTHS}\b"), 1)
    _apply(chars, re.compile(rf"\b{MONTHS}\.?\s+([0-9]{{1,2}})(?:st|nd|rd|th)?\b(?!\s*%)(?![0-9,])"), 1)
    _apply(chars, re.compile(rf"\b{MONTHS}\.?\s+(?:[0-9]{{1,2}},?\s+)?([0-9]{{4}})\b"), 1)
    _apply(chars, re.compile(r"\b[HQ][1-4]\s+([0-9]{4})\b"), 1)

    # 5. cross-references and heading numbers
    _apply(chars, REF_LIST_RE, 1)
    _apply(chars, REF_SINGLE_RE, 1)
    _apply(chars, re.compile(r"^(#{1,6}\s+)([A-Z]?[0-9]+(?:\.[0-9]+)*\.?)"), 2)

    # 6. LaTeX fractions, sub- and superscripts
    _apply(chars, re.compile(r"\\[td]?frac\s*(?:[0-9]\s*[0-9]|\{[0-9]+\}\s*\{[0-9]+\})"))
    _apply(chars, re.compile(r"[_^]\s*\{(?:[^{}]|\{[^{}]*\})*\}"),
           keep=lambda m: re.fullmatch(rf"\^\s*\{{\s*[-{MINUS}]?[0-9]*\.[0-9]+\s*\}}", m.group(0)) is not None)
    _apply(chars, re.compile(r"[_^]\s*[0-9A-Za-z]"))

    # 4. identifiers, theorem-style labels, model orders and Landau symbols
    _apply(chars, re.compile(r"(?<![\w])[^\W\d_]\w*?[0-9][\w]*(?:\.[0-9]+[A-Za-z]?)*"))
    _apply(chars, re.compile(r"(?<![\w.])[A-Z]\.[0-9]+(?:\.[0-9]+)*(?![0-9])"))
    _apply(chars, re.compile(r"(?<=[A-Za-z])\(\s*[0-9]+(?:\s*,\s*[0-9]+)*\s*\)"))
    _apply(chars, INDEX_NAME_RE, 1)

    work = "".join(chars)
    for m in NUM_RE.finditer(work):
        sgn = _sign(m.group("sign"), line, m.start())
        start = m.start() if sgn < 0 or m.group("sign") == "+" else m.start("int")
        int_part = m.group("int")
        frac = m.group("frac") or ""
        mant_text = int_part.replace("{,}", "").replace(",", "") + frac
        pct = bool(m.group("pct"))
        end = m.end()
        if not pct and re.match(r"\s*(?:percentage\s+)?points?\b", work[end:]):
            pct = True
        tokens.append(NumberToken(lineno, start, end, line[start:end].strip(), sgn * float(mant_text),
                                  len(frac) - 1 if frac else 0, percent=pct,
                                  has_sep=("," in int_part), ordinal=bool(m.group("ord"))))

    tokens.sort(key=lambda t: t.start)

    # ranges: the left end inherits a percent sign or a power of ten from the right end
    for a, b in zip(tokens, tokens[1:]):
        gap = _between(line, a, b)
        if re.fullmatch(r"\s*\$?\s*to\s*\$?\s*", gap) or (re.fullmatch(r"\s*\$?\s*(?:and|-)\s*\$?\s*", gap)
                                                         and (b.percent or b.sci) and a.decimals == b.decimals):
            if b.percent and not a.percent:
                a.percent = True
            if b.sci and not a.sci and b.exponent != 0:
                a.sci, a.exponent = True, b.exponent

    # years
    kept: list[NumberToken] = []
    year_flags: list[bool] = []
    for k, tok in enumerate(tokens):
        is_year = False
        if _year_like(tok):
            after = work[tok.end:]
            before = line[:tok.start]
            prev_year = (k > 0 and year_flags[k - 1] and _is_range_gap(_between(line, tokens[k - 1], tok)))
            if re.match(r"\s*-\s*[A-Za-z]", after):
                is_year = False
            elif _is_table_cell_alone(line, tok.start, tok.end):
                is_year = False
            elif re.search(r"(?:[=<>\u2264\u2265]|\\le|\\ge|\\leq|\\geq)\s*\$?\s*$", before):
                is_year = False
            elif prev_year:
                is_year = True
            else:
                w = re.match(r"\s+([a-z][a-z-]*)", after)
                is_year = not (w and w.group(1) not in YEAR_CONNECTORS)
        year_flags.append(is_year)

    def small(t: NumberToken) -> bool:
        return (t.decimals == 0 and not t.percent and not t.sci and not t.has_sep and t.exponent == 0
                and abs(t.mantissa) <= 10)

    breaks = [] if line.lstrip().startswith("|") else [m.end() for m in SENTENCE_BREAK_RE.finditer(line)]
    for k, tok in enumerate(tokens):
        if year_flags[k]:
            continue
        if not keep_small_ints and small(tok):
            partners = []
            if k > 0 and _is_range_gap(_between(line, tokens[k - 1], tok)):
                partners.append(tokens[k - 1])
            if k + 1 < len(tokens) and _is_range_gap(_between(line, tok, tokens[k + 1])):
                partners.append(tokens[k + 1])
            if not any(not small(p) and not year_flags[tokens.index(p)] for p in partners):
                continue
        tok.context = _context(line, tok.start, tok.end)
        tok.unit = sum(1 for b in breaks if b <= tok.start)
        kept.append(tok)
    return kept


def extract_numbers(text: str, keep_small_ints: bool = False) -> list[NumberToken]:
    """All result-like numbers of a Markdown document, in reading order."""
    lines = text.splitlines()
    skip, fenced = _skip_regions(lines)
    out: list[NumberToken] = []
    in_comment = False
    for k, line in enumerate(lines):
        if k in skip:
            continue
        if in_comment:
            if "-->" in line:
                in_comment = False
                line = " " * (line.index("-->") + 3) + line[line.index("-->") + 3:]
            else:
                continue
        if "<!--" in line and "-->" not in line[line.index("<!--"):]:
            in_comment = True
            line = line[:line.index("<!--")]
        if k in fenced:
            continue
        out.extend(extract_line(line, k + 1, keep_small_ints=keep_small_ints))
    return out


# ---------------------------------------------------------------------------
# Scope: experiment labels of each line
# ---------------------------------------------------------------------------

def labels_in(text: str) -> set[str]:
    out = set()
    for m in LABEL_RE.finditer(text):
        if m.group("alt"):
            out.add(f"ALT{int(m.group('alt'))}")
        elif m.group("th"):
            out.add(f"TH{int(m.group('th'))}")
        else:
            out.add(f"E{int(m.group('n'))}")
    return out


def sentence_units(line: str) -> list[str]:
    """The sentences of a line of prose (a table row is one unit), as used for NumberToken.unit."""
    if line.lstrip().startswith("|"):
        return [line]
    cuts = [0] + [m.end() for m in SENTENCE_BREAK_RE.finditer(line)] + [len(line)]
    return [line[a:b] for a, b in zip(cuts, cuts[1:])]


def line_scopes(text: str) -> dict[tuple[int, int], tuple[set[str], set[str]]]:
    """Scopes of each sentence, keyed by (line, unit) with 1-based lines: (narrow, wide).

    narrow: the labels of the sentence itself; if it has none, those of its paragraph (for a table: the
    table, its caption and the paragraph before it); if that has none, those of its section.
    wide: the labels of the paragraph and of the section. Both are empty when the section has none."""
    lines = text.splitlines()
    n = len(lines)
    # blocks: runs of non-blank lines; sections: runs between headings
    block_of = [-1] * n
    blocks: list[list[int]] = []
    for k, line in enumerate(lines):
        if line.strip() == "" or re.match(r"^#{1,6}\s", line):
            continue
        if k > 0 and block_of[k - 1] >= 0 and lines[k - 1].strip() != "":
            block_of[k] = block_of[k - 1]
            blocks[-1].append(k)
        else:
            block_of[k] = len(blocks)
            blocks.append([k])
    section_of = [0] * n
    sec = 0
    for k, line in enumerate(lines):
        if re.match(r"^#{1,6}\s", line):
            sec += 1
        section_of[k] = sec
    sec_labels: dict[int, set[str]] = {}
    for k, line in enumerate(lines):
        sec_labels.setdefault(section_of[k], set()).update(labels_in(line))
    block_labels = [set().union(*(labels_in(lines[i]) for i in b)) for b in blocks]

    out: dict[tuple[int, int], tuple[set[str], set[str]]] = {}
    for k in range(n):
        b = block_of[k]
        sec = set(sec_labels.get(section_of[k], set()))
        if b < 0:
            out[(k + 1, 0)] = (sec, sec)
            continue
        labs = set(block_labels[b])
        is_table = lines[blocks[b][0]].lstrip().startswith("|")
        if is_table:
            if b + 1 < len(blocks) and re.match(r"^\s*[*_]?(Table|Figure)\b", lines[blocks[b + 1][0]]):
                labs |= block_labels[b + 1]
            if b > 0:
                labs |= block_labels[b - 1]
        wide = labs | sec
        if not labs:
            labs = sec
        for u, sentence in enumerate(sentence_units(lines[k])):
            own = set() if is_table else labels_in(sentence)
            out[(k + 1, u)] = (own or labs, wide)
    return out


# ---------------------------------------------------------------------------
# Corpus
# ---------------------------------------------------------------------------


THEORY_SCOPE = {"theory_running_tau": {"E14"}, "theory_panel": {"E15"}, "theory_multifreq": {"E17"},
                "theory_delay": {"E11"}, "theory_robustness": {"TH4"}, "theory_limits": {"TH5"},
                "theory_per_period": {"TH6"}}
SRC_NUM_RE = re.compile(r"(?<![A-Za-z0-9_.])(-?)([0-9]{1,3}(?:,[0-9]{3})+(?![0-9])|[0-9]+)(\.[0-9]+)?(?![0-9])(?!\.[0-9])")
COMMON_HITS = 5000   # a value matching more cells than this in one file says nothing about co-occurrence
COMMON_SEGMENTS = 5  # a number in this many sentences of the cited experiment's reports is one of its settings
ROW_MIN = 3          # table rows with fewer distinct numbers are not compared row by row
ROW_SHARE = 0.75     # share of a row's numbers the best source row must hold before a mismatch is reported


@dataclass
class Source:
    """One source file: its numbers, and the segment (Markdown line, CSV row) each one sits in."""

    path: str                       # repository-relative, POSIX
    kind: str                       # "md" or "data"
    tier: str                       # primary, other, review, results or alternatives
    labels: frozenset = frozenset()
    literals: dict = field(default_factory=dict)                     # literal key -> segment ids
    seg_lines: list = field(default_factory=list)                    # segment id -> 1-based line (md)
    row_segs: frozenset = frozenset()                                # segment ids that are table rows (md)
    values: np.ndarray = field(default_factory=lambda: np.empty(0))  # sorted values
    segs: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))  # segment of each value


def _md_labels(path: Path) -> frozenset:
    stem = path.stem
    if stem in THEORY_SCOPE:
        return frozenset(THEORY_SCOPE[stem])
    m = re.fullmatch(r"E(\d+)_E(\d+)", stem)
    if m:
        return frozenset(f"E{i}" for i in range(int(m.group(1)), int(m.group(2)) + 1))
    m = re.fullmatch(r"E(\d+)", stem)
    if m:
        return frozenset({f"E{int(m.group(1))}"})
    m = re.fullmatch(r"ALT(\d+)", stem)
    if m:
        return frozenset({f"ALT{int(m.group(1))}"})
    m = re.fullmatch(r"TH(\d+)", stem)
    if m:
        return frozenset({f"TH{int(m.group(1))}"})
    return frozenset()


def _md_tier(rel: str) -> str:
    if rel.startswith("docs/review/") or rel == "docs/paper_review.md":
        return "review"
    if (rel.startswith("docs/report_parts/") or re.fullmatch(r"docs/theory[^/]*\.md", rel)
            or rel == "docs/math.md" or re.fullmatch(r"docs/preregistration[^/]*\.md", rel)):
        return "primary"
    return "other"


def source_segments(text: str) -> list[tuple[int, str]]:
    """(1-based line, sentence or table row) of a source text, after normalising number spellings."""
    text = text.replace("{,}", ",").replace(MINUS, "-").replace("\\%", "%")
    text = re.sub(r"\bminus\s+(?=[0-9])", "-", text)
    return [(i + 1, u) for i, line in enumerate(text.splitlines()) for u in sentence_units(line)]


def text_numbers(text: str):
    """Yield (segment, literal key, value) for every number of a source text.

    A segment is a sentence, or a table row (see source_segments)."""
    for seg, (_, line) in enumerate(source_segments(text)):
        for pattern in (SCI_TIMES_RE, SCI_E_RE, SCI_WORD_RE):
            for m in pattern.finditer(line):
                mant = m.group("mant")
                e = WORD_EXP[m.group("word")] if pattern is SCI_WORD_RE else _exp_value(m)
                sgn = -1 if m.group("sign") == "-" else 1
                yield seg, sci_key(sgn * float(mant), _decimals(mant), e), sgn * float(mant) * 10.0 ** e
        for m in SCI_POWER_RE.finditer(line):
            e = _exp_value(m)
            yield seg, sci_key(1.0, 0, e), 10.0 ** e
        for m in SRC_NUM_RE.finditer(line):
            core = m.group(2).replace(",", "") + (m.group(3) or "")
            value = float(core)
            yield seg, core, value
            if m.group(1):
                yield seg, "-" + core, -value


def _cell_value(cell: str) -> float | None:
    s = cell.strip()
    if not s or len(s) > 40:
        return None
    s = s.rstrip("%").strip()
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def _json_values(obj, out: list[float]) -> None:
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        if math.isfinite(obj):
            out.append(float(obj))
    elif isinstance(obj, str):
        v = _cell_value(obj)
        if v is not None:
            out.append(v)
    elif isinstance(obj, dict):
        for v in obj.values():
            _json_values(v, out)
    elif isinstance(obj, list):
        for v in obj:
            _json_values(v, out)


def _sorted_arrays(pairs: list[tuple[int, float]]) -> tuple[np.ndarray, np.ndarray]:
    if not pairs:
        return np.empty(0), np.empty(0, dtype=np.int64)
    segs = np.fromiter((p[0] for p in pairs), dtype=np.int64, count=len(pairs))
    vals = np.fromiter((p[1] for p in pairs), dtype=float, count=len(pairs))
    order = np.argsort(vals, kind="stable")
    return vals[order], segs[order]


def md_source(path: Path, rel: str, tier: str, labels: frozenset) -> Source:
    literals: dict[str, set[int]] = {}
    pairs: list[tuple[int, float]] = []
    text = path.read_text(encoding="utf-8", errors="replace")
    for seg, key, value in text_numbers(text):
        literals.setdefault(key, set()).add(seg)
        pairs.append((seg, value))
    vals, segs = _sorted_arrays(pairs)
    segments = source_segments(text)
    return Source(rel, "md", tier, labels, literals=literals, values=vals, segs=segs,
                  seg_lines=[ln for ln, _ in segments],
                  row_segs=frozenset(i for i, (_, u) in enumerate(segments) if u.lstrip().startswith("|")))


def data_source(path: Path, rel: str, tier: str, labels: frozenset) -> Source:
    name = path.name
    pairs: list[tuple[int, float]] = []
    if name.endswith((".csv", ".csv.gz")):
        opener = gzip.open if name.endswith(".gz") else open
        with opener(path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
            for seg, row in enumerate(csv.reader(fh)):
                for cell in row:
                    v = _cell_value(cell)
                    if v is not None:
                        pairs.append((seg, v))
    elif name.endswith(".json"):
        vals: list[float] = []
        with open(path, encoding="utf-8") as fh:
            _json_values(json.load(fh), vals)
        pairs = [(0, v) for v in vals]
    elif name.endswith(".txt"):
        pairs = [(seg, v) for seg, _, v in text_numbers(path.read_text(encoding="utf-8", errors="replace"))]
    vals_arr, segs = _sorted_arrays(pairs)
    return Source(rel, "data", tier, labels, values=vals_arr, segs=segs)


def _data_labels(path: Path) -> frozenset:
    m = re.match(r"e(\d{2})", path.name)
    if m:
        return frozenset({f"E{int(m.group(1))}"})
    m = re.match(r"th(\d)_", path.name)
    return frozenset({f"TH{int(m.group(1))}"}) if m else frozenset()


DATA_SUFFIXES = (".csv", ".csv.gz", ".json", ".txt")


@dataclass
class Corpus:
    md: list[Source]
    data: list[Source]
    extra: list[Source] = field(default_factory=list)   # alternatives/*/results: leads only
    cache: dict = field(default_factory=dict)


def build_corpus(root: Path, target: Path, with_alternatives: bool = False) -> Corpus:
    """Markdown sources under docs/ and data files under results/ (see the module docstring)."""
    target = target.resolve()
    md: list[Source] = []
    for p in sorted((root / "docs").rglob("*.md")):
        rel = p.relative_to(root).as_posix()
        if p.resolve() == target or rel == "docs/experiment_report.md" or rel.startswith("paper/"):
            continue
        md.append(md_source(p, rel, _md_tier(rel), _md_labels(p)))
    data: list[Source] = []
    res = root / "results"
    if res.is_dir():
        for p in sorted(res.rglob("*")):
            rel = p.relative_to(root).as_posix()
            if p.is_file() and not rel.startswith("results/trace/") and p.name.endswith(DATA_SUFFIXES):
                data.append(data_source(p, rel, "results", _data_labels(p)))
    # alternatives/<folder>/results, labelled with the ALT part that names the folder
    alt_of: dict[str, set[str]] = {}
    for part in sorted((root / "docs" / "report_parts").glob("ALT*.md")):
        for folder in set(re.findall(r"alternatives/([A-Za-z0-9_]+)", part.read_text(encoding="utf-8"))):
            alt_of.setdefault(folder, set()).add(part.stem)
    extra: list[Source] = []
    alt_root = root / "alternatives"
    if alt_root.is_dir():
        for p in sorted(alt_root.glob("*/results/**/*")):
            rel = p.relative_to(root).as_posix()
            if p.is_file() and p.name.endswith(DATA_SUFFIXES):
                labels = frozenset(alt_of.get(rel.split("/")[1], set()))
                extra.append(data_source(p, rel, "alternatives", labels))
    if with_alternatives:
        data, extra = data + extra, []
    return Corpus(md, data, extra)


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def _intervals(tok: NumberToken, data: bool) -> list[tuple[float, float, str]]:
    """Value intervals that trace tok (see the rounding rule in the module docstring)."""
    half = 0.5 * 10.0 ** (tok.exponent - tok.decimals)
    centre = tok.mantissa * 10.0 ** tok.exponent
    lo, hi = centre - half, centre + half
    out = [(lo, hi, "rounded")]
    if tok.percent:
        out.append((lo / 100, hi / 100, "rounded, source x 100"))
    elif data and 0 < centre <= 1:
        out.append((lo * 100, hi * 100, "rounded, source / 100"))
    return out


def value_match(values: np.ndarray, mantissa: float, decimals: int, exponent: int = 0, *,
                percent: bool = False, data: bool = False) -> str | None:
    """How a sorted array of source values traces a number shown with `decimals` decimals, or None.

    The source value itself is always tried. A percentage is also compared with 100 times the source
    (a share stored as a fraction); a number between 0 and 1 that is not a percentage is also compared
    with the source divided by 100 when the source is a data file (a share stored in percent).
    """
    if values.size == 0:
        return None
    tok = NumberToken(0, 0, 0, "", mantissa, decimals, exponent, percent)
    for a, b, how in _intervals(tok, data):
        tol = 1e-9 * max(abs(a), abs(b)) + 1e-15
        if np.searchsorted(values, b + tol, "right") > np.searchsorted(values, a - tol, "left"):
            return how
    return None


def token_match(src: Source, tok: NumberToken) -> str | None:
    return value_match(src.values, tok.mantissa, tok.decimals, tok.exponent, percent=tok.percent,
                       data=src.kind == "data")


def literal_match(literals, tok: NumberToken) -> bool:
    return tok.literal_key() in literals


def _value_segs(src: Source, tok: NumberToken, cache: dict) -> set[int] | None:
    """Segments of src holding a value that traces tok; None when there are none or too many."""
    key = (src.path, tok.mantissa, tok.decimals, tok.exponent, tok.percent)
    if key in cache:
        return cache[key]
    out: set[int] = set()
    n = 0
    for a, b, _ in _intervals(tok, src.kind == "data"):
        tol = 1e-9 * max(abs(a), abs(b)) + 1e-15
        i = int(np.searchsorted(src.values, a - tol, "left"))
        j = int(np.searchsorted(src.values, b + tol, "right"))
        n += j - i
        if n > COMMON_HITS:
            cache[key] = None
            return None
        out.update(src.segs[i:j].tolist())
    cache[key] = out or None
    return cache[key]


def cooccurs(src: Source, tok: NumberToken, others: list[NumberToken], cache: dict, need: int = 1) -> bool:
    """True if one segment of src (a Markdown sentence or table row, a CSV row) holds tok and at least
    `need` other numbers of the same paper sentence or table row."""
    distinct = list({o.literal_key(): o for o in others if o.literal_key() != tok.literal_key()}.values())
    if len(distinct) < need:
        return False
    if src.kind == "md":
        segs = src.literals.get(tok.literal_key())
        sets = [src.literals.get(o.literal_key(), set()) for o in distinct]
    else:
        segs = _value_segs(src, tok, cache)
        sets = [_value_segs(src, o, cache) or set() for o in distinct]
    if not segs:
        return False
    if need == 1:
        return any(segs & s for s in sets)
    counts: dict[int, int] = {}
    for s in sets:
        for seg in segs & s:
            counts[seg] = counts.get(seg, 0) + 1
    return bool(counts) and max(counts.values()) >= need


def _perturbed(tok: NumberToken) -> list[NumberToken]:
    step = 10.0 ** (-tok.decimals)
    out = []
    for k in (-7, -3, 3, 7):
        m = round(tok.mantissa + k * step, tok.decimals)
        out.append(NumberToken(tok.line, tok.start, tok.end, "", m, tok.decimals, tok.exponent, tok.percent,
                               sci=tok.sci))
    return out


def _weak(src: Source, tok: NumberToken, literal: bool) -> bool:
    hits = 0
    for p in _perturbed(tok):
        hits += literal_match(src.literals, p) if literal else token_match(src, p) is not None
    return hits >= 2


@dataclass
class Trace:
    token: NumberToken
    status: str
    source: str = ""
    how: str = ""
    tier: str = ""
    in_scope: str = ""
    weak: str = ""
    cooccur: str = ""
    scope: str = ""
    hint: str = ""
    row_check: str = ""


def _in_scope(src: Source, scope: set[str] | None) -> bool:
    return scope is None or bool(src.labels & scope)


def _md_order(corpus: Corpus, scope: set[str] | None, inside: bool) -> list[Source]:
    """Markdown sources in or out of scope, primary tier first and review tier last."""
    rank = {"primary": 0, "other": 1, "review": 2}
    pick = [s for s in corpus.md if _in_scope(s, scope) == inside]
    return sorted(pick, key=lambda s: rank[s.tier])


def _scope_text(scope: set[str] | None) -> str:
    if scope is None:
        return "all"
    return " ".join(sorted(scope, key=lambda s: (s[0], int(re.sub(r"\D", "", s)))))


def find_cooccurrence(tok: NumberToken, others: list[NumberToken], scope: set[str] | None,
                      corpus: Corpus) -> str:
    """The first in-scope Markdown source (not the review tier) in which one sentence or table row holds
    tok and another number of tok's own sentence or table row. Returns the path, "single" when the
    sentence has no other number, or "" when there is no such source.

    When there is none, "common" means that the number appears in at least COMMON_SEGMENTS sentences
    of the in-scope reports (a setting of the experiment, such as 390 steps). Data rows are not used:
    a CSV with thousands of rows holds two given short numbers in one row by chance too often for the
    test to mean anything."""
    others = [o for o in others if o.literal_key() != tok.literal_key()]
    if not others:
        return "single"
    in_scope = [s for s in _md_order(corpus, scope, True) if s.tier != "review"]
    for src in in_scope:
        if cooccurs(src, tok, others, corpus.cache):
            return src.path
    if scope is not None and sum(len(s.literals.get(tok.literal_key(), ())) for s in in_scope) >= COMMON_SEGMENTS:
        return "common"
    return ""


def trace_token(tok: NumberToken, scope: set[str] | None, corpus: Corpus,
                others: list[NumberToken] | None = None, wide: set[str] | None = None) -> Trace:
    """Trace one number. scope: the narrow scope (None: every source is in scope); wide: the labels of
    the paragraph and section, used for literal Markdown matches only."""
    tr = _trace_status(tok, scope, corpus, wide)
    tr.scope = _scope_text(scope)
    if others is not None:
        tr.cooccur = find_cooccurrence(tok, others, scope, corpus) or "no"
        if tr.cooccur == "no" and scope is not None and wide and wide - scope:
            tr.cooccur = find_cooccurrence(tok, others, wide, corpus) or "no"
    return tr


def _trace_status(tok: NumberToken, scope: set[str] | None, corpus: Corpus,
                  wide: set[str] | None = None) -> Trace:
    yes = lambda flag: "yes" if flag else ""  # noqa: E731
    # 1. literal, Markdown in scope
    for src in _md_order(corpus, scope, True):
        if literal_match(src.literals, tok):
            return Trace(tok, "traced_md", src.path, "literal", src.tier, "yes", yes(_weak(src, tok, True)))
    # 2. literal, Markdown of an experiment named elsewhere in the paragraph or section
    if scope is not None and wide and wide - scope:
        for src in _md_order(corpus, wide - scope, True):
            if src.tier != "review" and literal_match(src.literals, tok):
                return Trace(tok, "traced_md", src.path, "literal", src.tier, "paragraph",
                             yes(_weak(src, tok, True)))
    # 3. rounded, Markdown in scope (not the review tier)
    for src in _md_order(corpus, scope, True):
        if src.tier == "review":
            continue
        how = token_match(src, tok)
        if how:
            return Trace(tok, "traced_md", src.path, how, src.tier, "yes", yes(_weak(src, tok, False)))
    # 4. rounded, data in scope
    for src in corpus.data:
        if _in_scope(src, scope):
            how = token_match(src, tok)
            if how:
                return Trace(tok, "traced_csv", src.path, how, src.tier, "yes", yes(_weak(src, tok, False)))
    # 5. literal, Markdown out of scope: accepted only when the match is not weak
    hint = ""
    for src in _md_order(corpus, scope, False):
        if literal_match(src.literals, tok):
            if not _weak(src, tok, True):
                return Trace(tok, "traced_md", src.path, "literal", src.tier, "no", "")
            hint = hint or f"{src.path} (literal, weak)"
    # untraced: leave a lead, first among the alternatives' results, then outside the scope
    for src in corpus.extra:
        how = token_match(src, tok)
        if how and (scope is None or _in_scope(src, scope)):
            hint = f"{src.path} ({how})" + (f"; {hint}" if hint else "")
            break
    if not hint and scope is not None:
        for src in corpus.data + [s for s in corpus.md if s.tier != "review"]:
            if not _in_scope(src, scope):
                how = token_match(src, tok)
                if how:
                    hint = f"{src.path} ({how})"
                    break
    return Trace(tok, "untraced", hint=hint)


def row_check(tokens: list[NumberToken], scope: set[str] | None, corpus: Corpus) -> dict[int, str]:
    """Compare a table row of the target with the source table row that shares most of its numbers.

    Only rows with at least ROW_MIN distinct numbers are checked, against Markdown table rows in scope
    (outside the review tier). When the best source row holds at least ROW_MIN of them and at least
    ROW_SHARE of them, every number of the row missing from it is reported (key: index in tokens): such
    a number was rounded differently, comes from another run, or is wrong. Rows that merge several
    source rows rarely reach ROW_SHARE, so they are not reported."""
    keys = {t.literal_key() for t in tokens}
    if len(keys) < ROW_MIN:
        return {}
    best = (0, None, -1)
    for src in corpus.md:
        if src.tier == "review" or not _in_scope(src, scope):
            continue
        counts: dict[int, int] = {}
        for key in keys:
            for seg in src.literals.get(key, ()):
                if seg in src.row_segs:
                    counts[seg] = counts.get(seg, 0) + 1
        for seg, c in counts.items():
            if c > best[0]:
                best = (c, src, seg)
    c, src, seg = best
    if src is None or c < ROW_MIN or c < ROW_SHARE * len(keys):
        return {}
    where = f"{src.path}:{src.seg_lines[seg]}"
    held = {k for k in keys if seg in src.literals.get(k, ())}
    return {i: f"not in {where} ({c} of {len(keys)} numbers of the row are)"
            for i, t in enumerate(tokens) if t.literal_key() not in held}


def trace_document(target: Path, corpus: Corpus, *, keep_small_ints: bool = False,
                   use_scope: bool = True) -> list[Trace]:
    text = target.read_text(encoding="utf-8")
    tokens = extract_numbers(text, keep_small_ints=keep_small_ints)
    scopes = line_scopes(text)
    units: dict[tuple[int, int], list[NumberToken]] = {}
    for tok in tokens:
        units.setdefault((tok.line, tok.unit), []).append(tok)
    out = []
    for tok in tokens:
        narrow, wide = scopes.get((tok.line, tok.unit), (set(), set()))
        scope = (narrow or None) if use_scope else None
        others = [o for o in units[(tok.line, tok.unit)] if o is not tok]
        out.append(trace_token(tok, scope, corpus, others, wide if use_scope else None))
    lines = text.splitlines()
    by_token = {id(t.token): t for t in out}
    for (line, unit), toks in units.items():
        if not lines[line - 1].lstrip().startswith("|"):
            continue
        narrow, wide = scopes.get((line, unit), (set(), set()))
        scope = ((narrow | wide) or None) if use_scope else None
        for i, msg in row_check(toks, scope, corpus).items():
            by_token[id(toks[i])].row_check = msg
    return out


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

COLUMNS = ["line", "number", "context", "status", "source", "how", "tier", "in_scope", "weak", "cooccur",
           "row_check", "scope", "hint"]


def default_output(target: Path, root: Path) -> Path:
    name = "paper" if target.resolve() == (root / "paper" / "short_paper.md").resolve() else target.stem.lower()
    return root / "results" / "trace" / f"{name}_trace.csv"


def write_csv(traces: list[Trace], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(COLUMNS)
        for t in traces:
            w.writerow([t.token.line, t.token.text, t.token.context, t.status, t.source, t.how, t.tier,
                        t.in_scope, t.weak, t.cooccur, t.row_check, t.scope, t.hint])


def summarise(traces: list[Trace], target: str) -> str:
    n = len(traces)
    traced = [t for t in traces if t.status != "untraced"]
    md_lit = sum(t.status == "traced_md" and t.how == "literal" for t in traces)
    md_rnd = sum(t.status == "traced_md" and t.how != "literal" for t in traces)
    out_scope = sum(t.in_scope == "no" for t in traced)
    para = sum(t.in_scope == "paragraph" for t in traced)
    csv_n = sum(t.status == "traced_csv" for t in traces)
    un = n - len(traced)
    weak = sum(t.weak == "yes" for t in traced)
    co = sum(t.cooccur not in ("", "no", "single", "common") for t in traced)
    single = sum(t.cooccur == "single" for t in traced)
    common = sum(t.cooccur == "common" for t in traced)
    alone = sum(t.cooccur == "no" for t in traced)
    return "\n".join([
        f"{target}: {n} numbers",
        f"  traced      {len(traced)}",
        f"    traced_md   {md_lit + md_rnd} (literal {md_lit}, rounded {md_rnd}; in a report the paragraph cites {para};"
        f" outside the cited experiments {out_scope})",
        f"    traced_csv  {csv_n}",
        f"  untraced    {un}",
        f"  in context (an in-scope Markdown sentence or row also holds another number of the same sentence):",
        f"    {co}; a setting of the experiment (common) {common}; only number of its sentence {single};"
        f" not in context {alone}",
        f"  table cells missing from the best-matching source row: {sum(bool(t.row_check) for t in traces)}",
        f"  weak matches (nearby wrong values also match the same file): {weak}"])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--target", default="paper/short_paper.md", help="Markdown file to check")
    ap.add_argument("--out", default=None, help="output CSV (default results/trace/<name>_trace.csv)")
    ap.add_argument("--keep-small-ints", action="store_true", help="keep stand-alone integers 0 to 10")
    ap.add_argument("--no-scope", action="store_true", help="match every number against every source")
    ap.add_argument("--with-alternatives", action="store_true",
                    help="also trace against alternatives/*/results (otherwise they only give hints)")
    ap.add_argument("--root", default=str(ROOT), help="repository root")
    ap.add_argument("--quiet", action="store_true", help="do not list the untraced numbers")
    ap.add_argument("--list-alone", action="store_true",
                    help="also list traced numbers that share no source line or row with their sentence")
    ap.add_argument("--fail-on-untraced", action="store_true", help="exit with status 1 if a number is untraced")
    args = ap.parse_args(argv)

    root = Path(args.root).resolve()
    target = Path(args.target)
    if not target.is_absolute():
        target = root / target
    corpus = build_corpus(root, target, with_alternatives=args.with_alternatives)
    traces = trace_document(target, corpus, keep_small_ints=args.keep_small_ints, use_scope=not args.no_scope)
    out = Path(args.out) if args.out else default_output(target, root)
    if not out.is_absolute():
        out = root / out
    write_csv(traces, out)
    rel = target.relative_to(root).as_posix() if target.is_relative_to(root) else str(target)
    print(summarise(traces, rel))
    print(f"  sources: {len(corpus.md)} Markdown files, {len(corpus.data)} data files"
          + (f" ({len(corpus.extra)} files under alternatives/ used for hints only)" if corpus.extra else ""))
    print(f"wrote {out.relative_to(root) if out.is_relative_to(root) else out}")
    if not args.quiet:
        groups = [("Untraced numbers", lambda t: t.status == "untraced")]
        if args.list_alone:
            groups.append(("Traced, but sharing no source line or row with their sentence",
                           lambda t: t.status != "untraced" and t.cooccur == "no"))
        groups.append(("Table cells missing from the best-matching source row", lambda t: bool(t.row_check)))
        for title, keep in groups:
            rows = [t for t in traces if keep(t)]
            if rows:
                print(f"\n{title} (line, number, context, hint or source):")
                for t in rows:
                    note = t.row_check if title.startswith("Table") else (t.hint or t.source)
                    print(f"  {t.token.line:4d}  {t.token.text:<14} {t.token.context}" + (f"  [{note}]" if note else ""))
    if args.fail_on_untraced and any(t.status == "untraced" for t in traces):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
