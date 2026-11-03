# Short paper draft

Draft of a short paper on the monitor in this repository: *Anytime-valid detection of volatility changes in high-frequency returns with robust truncation* (Shreyan Thanki, draft of November 2026).

## Files and how they relate

| File | What it is | Edit it? |
|---|---|---|
| `short_paper.md` | The paper. Markdown with `$...$` and `$$...$$` LaTeX math, readable on GitHub. This is the source of truth. | Yes |
| `short_paper.tex` | Standalone LaTeX generated from `short_paper.md` by pandoc. | No: regenerate it |
| `references.bib` | The same references as the hand-written list at the end of `short_paper.md`, in BibTeX, for a journal template. Every entry was checked against a publisher, journal, arXiv or index listing; authors are given with initials. | Keep in step with the Markdown list |
| `make_paper_figures.py` | Makes the one figure that has no counterpart in `results/figures/` (Figure 1, the signal cap, from `results/e10d_signal_cap.csv`; the simulated ratio of means is the pooled column `ratio_of_means_pooled`, not the per-path average). | Yes |
| `figures/signal_cap.png` | Output of `make_paper_figures.py`. | No: regenerate it |

The other two figures are read from `results/figures/` (`e09_stress.png` for Figure 2 and `e11_delay_predictions.png` for Figure S1); the experiment scripts in `experiments/` make them. Image paths in the paper are relative to this folder (`figures/...` and `../results/figures/...`), so the Markdown renders on GitHub and the LaTeX compiles from this folder.

## Structure

Main text: abstract; 1 introduction; 2 setting and the truncated statistic; 3 betting e-detectors (Theorems A and B, two-sided rule, Corollary C, the one-sided Bahadur-Savage argument, the confidence-sequence baseline with Theorem D and Corollaries E and F); 4 detection delay (information ladder, speed limits, Theorems 8 and 9, the price of robustness, the signal cap of Proposition 10); 5 simulation evidence; 6 real data; 7 questions from the reproduction of the estimator; 8 discussion, limitations and open problems. Then references, Appendices A (proofs for the estimated baseline), B (proofs for the delay theory), C (experiment settings), D (where each result lives in the repository), and a supplement S1 to S6 with the detailed versions of results that the main text summarises.

Conventions: every number cites its experiment label (E1 to E17, ALT1 to ALT5, as in `docs/experiment_report.md`) and is marked **[SYNTHETIC]**, **[REAL]** or **[PUBLISHED]**. All experiments cited in the draft are finished; there are no placeholders left. **[HEURISTIC]** marks the one approximation that is not proved.

## Regenerating

From the repository root:

```bash
python paper/make_paper_figures.py          # writes paper/figures/signal_cap.png
cd paper
pandoc -s --shift-heading-level-by=-1 short_paper.md -o short_paper.tex
```

`--shift-heading-level-by=-1` turns the Markdown `##` headings into LaTeX `\section`. The section numbers are written in the headings, and the pandoc template switches LaTeX numbering off, so they are not doubled. The `.tex` was produced with pandoc 3.1.3.

## Compiling the LaTeX elsewhere

No LaTeX engine was available where this draft was written, so `short_paper.tex` has not been compiled yet; a first run may show small layout problems, most likely in the wide tables. With a TeX distribution that has the usual packages of the pandoc template (`amsmath`, `amssymb`, `graphicx`, `longtable`, `booktabs`, `hyperref`, `xcolor`), run from this folder:

```bash
cd paper
pdflatex short_paper.tex && pdflatex short_paper.tex
# or: xelatex short_paper.tex (or lualatex), twice
```

or go straight from the Markdown:

```bash
cd paper
pandoc -s --shift-heading-level-by=-1 short_paper.md -o short_paper.pdf --pdf-engine=xelatex
```

For a journal template, copy the body of `short_paper.tex` (between `\begin{document}` and `\end{document}`) into the template, replace the hand-written reference list by `\bibliography{references}` with the template's author-year style (for example natbib), and turn the in-text citations into `\citet`/`\citep` commands. The keys in `references.bib` are of the form `surnameYEARword` (for example `shin2024edetectors`).

## Length

The main text (abstract to Section 8) has about 17,000 words including inline mathematics, plus tables and displayed equations, longer than a 10 to 14 page short paper. Results that the main text only summarises are already in the supplement. The next candidates for moving out of the main text are the construction of the confidence sequence in Section 3.4 (to Appendix A), the information rate and Table 2 in Section 4.5, and Sections 6.2 to 6.4.
