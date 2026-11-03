"""Build paper/short_paper.pdf from paper/short_paper.md without LaTeX.

Steps: pandoc turns the Markdown into a standalone HTML page with MathML formulas
(paper/build/short_paper.html, with the figures embedded), then a headless Chromium, driven by
Node's Playwright, prints that page to A4 PDF. If Node or Playwright is missing, the HTML page is
still written and can be printed to PDF from any browser.

Usage:
    python scripts/build_paper_pdf.py            # writes paper/build/short_paper.html and paper/short_paper.pdf
    python scripts/build_paper_pdf.py --html-only

Needs pandoc (2.19 or later). The LaTeX version, paper/short_paper.tex, is generated separately
with pandoc and is not used here.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
SRC = PAPER / "short_paper.md"
BUILD = PAPER / "build"
HTML = BUILD / "short_paper.html"
PDF = PAPER / "short_paper.pdf"

CSS = """
body { font-family: "DejaVu Serif", Georgia, serif; font-size: 10.5pt; line-height: 1.42;
       max-width: 46em; margin: 0 auto; color: #111; }
h1.title { font-size: 17pt; text-align: center; margin-bottom: 0.2em; }
p.author, p.date { text-align: center; margin: 0.1em; }
h2 { font-size: 13pt; margin-top: 1.4em; border-bottom: 1px solid #bbb; }
h3 { font-size: 11.5pt; margin-top: 1.1em; }
h4 { font-size: 10.5pt; }
table { border-collapse: collapse; margin: 0.8em auto; font-size: 8.8pt; }
th, td { border-top: 1px solid #999; border-bottom: 1px solid #999; padding: 2px 6px;
         vertical-align: top; }
thead th { border-bottom: 1.5px solid #333; }
figure { margin: 1em auto; text-align: center; page-break-inside: avoid; }
figure img { max-width: 100%; }
figcaption { font-size: 9pt; text-align: left; }
code { font-size: 9pt; }
math[display="block"] { margin: 0.5em 0; overflow-x: auto; }
tr, img { page-break-inside: avoid; }
"""

PRINT_JS = r"""
const { chromium } = require('playwright');
(async () => {
  const [html, pdf] = process.argv.slice(2);
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto('file://' + html, { waitUntil: 'load' });
  await page.pdf({ path: pdf, format: 'A4', printBackground: true,
                   margin: { top: '18mm', bottom: '18mm', left: '16mm', right: '16mm' },
                   displayHeaderFooter: true, headerTemplate: '<span></span>',
                   footerTemplate: '<div style="font-size:8px;width:100%;text-align:center;">' +
                                   '<span class="pageNumber"></span> / <span class="totalPages"></span></div>' });
  await browser.close();
})().catch(e => { console.error(e); process.exit(1); });
"""


def build_html() -> Path:
    if shutil.which("pandoc") is None:
        raise SystemExit("pandoc not found; install it from https://pandoc.org")
    BUILD.mkdir(parents=True, exist_ok=True)
    css = BUILD / "paper.css"
    css.write_text(CSS)
    cmd = ["pandoc", str(SRC), "--standalone", "--mathml", "--embed-resources",
           f"--resource-path={PAPER}", f"--css={css}", "--metadata", "lang=en",
           "-o", str(HTML)]
    subprocess.run(cmd, check=True, cwd=PAPER)
    return HTML


def print_pdf(html: Path) -> bool:
    if shutil.which("node") is None:
        print("node not found: open", html, "in a browser and print it to PDF")
        return False
    with tempfile.TemporaryDirectory() as tmp:
        js = Path(tmp) / "print.js"
        js.write_text(PRINT_JS)
        r = subprocess.run(["node", str(js), str(html), str(PDF)], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr.strip()[-2000:], file=sys.stderr)
        print("printing failed: open", html, "in a browser and print it to PDF")
        return False
    return True


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--html-only", action="store_true")
    a = ap.parse_args(argv)
    html = build_html()
    print("wrote", html.relative_to(ROOT))
    if not a.html_only and print_pdf(html):
        print("wrote", PDF.relative_to(ROOT))


if __name__ == "__main__":
    main()
