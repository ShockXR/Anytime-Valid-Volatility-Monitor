# Paper draft

Draft of a paper on the monitor in this repository: *Anytime-valid detection of volatility changes in high-frequency returns with robust truncation* (Shreyan Thanki, draft of November 2026), with a supplement and a companion note on extensions. This version follows the restructuring asked for in the referee report `docs/review/referee_report_paper.md` (Section 5 and items M1 to M8, 1 to 21 and m1 to m35).

## Files and how they relate

| File | What it is | Edit it? |
|---|---|---|
| `short_paper.md` | The paper: main text, references, then the supplement (S1 to S9). Markdown with `$...$` and `$$...$$` LaTeX math, readable on GitHub. This is the source of truth. | Yes |
| `extensions.md` | Companion note: running truncation level and restarts, many assets (e-BH), several sampling frequencies, adversarial null streams. The main text keeps one short paragraph on each (Sections 6.1 and 6.5). | Yes |
| `short_paper.tex` | Standalone LaTeX generated from `short_paper.md` by pandoc. | No: regenerate it |
| `references.bib` | The references of the paper and of the companion note in BibTeX, for a journal template. | Keep in step with the Markdown lists |
| `make_paper_figures.py` | Draws every figure of the paper from committed result files; it runs no simulation. | Yes |
| `figures/validity.png` | Figure 1 (from `results/e09ab_bound_tightness.csv` and `results/e09c_garch_band.csv`). | No: regenerate it |
| `figures/delay_predictions.png` | Figure S1 (from `results/e11a_delay.csv`). | No: regenerate it |
| `figures/signal_cap.png` | Figure S2 (from `results/e10d_signal_cap.csv`; the simulated ratio of means is the pooled column `ratio_of_means_pooled`). | No: regenerate it |

All figures now live in `paper/figures/`, with no experiment labels in their titles; the paper no longer reads anything from `results/figures/`.

## Structure

Main text: Abstract (no citations); 1 Introduction (three contributions and a related-work paragraph); 2 Model, statistic and nulls; 3 Validity (3.1 e-detectors and the guarantees, including their reading over many days; 3.2 a small drift of the conditional mean; 3.3 why a decrease detector needs a bounded statistic); 4 Estimated baseline (4.1 plug-in; 4.2 confidence sequence); 5 Detection delay; 6 Simulations; 7 Application (7.1 data and what real data can test; 7.2 five cryptocurrencies; 7.3 FOMC statements as the pre-specified validation; 7.4 other data and lessons); 8 Discussion; Reproducibility; References.

Supplement: S1 proofs and details for Sections 2 and 3 (Corollary 6; Theorem 7 with its design numbers, Table S1 and the sharpness result; drift that varies; Corollary 9; the short proofs of Lemma 1, Theorem 3 and Propositions 5 and 8); S2 proofs for Section 4; S3 proofs for Section 5; S4 second-order delay at practical levels (Table S2, Figure S1); S5 truncation constant, slack of the bound and bet design (Table S3); S6 clustering, forecasts, the GARCH bound and price noise (Proposition S8, Figure S2, Table S4); S7 real-data details (Nifty 50, FOMC and CPI, competitors with Table S5); S8 settings of every experiment (Tables S6 and S7); S9 where each table and figure comes from (script and result file).

Numbering: one counter for all results of the main text (Lemma 1 to Theorem 19); results stated only in the supplement are Lemma S1 to Proposition S8. Map from the previous draft (commit 8f71a2e): Theorem A to Theorem 3, Theorem B to Theorem 4, Proposition 1 to Proposition 5, Corollary C to Corollary 6, Proposition 2 to Proposition 8, Corollary 3 to Corollary 9, Theorem D to Theorem 11, Corollaries E and F to 12 and 13, Lemmas 4 and 5 to Lemmas 14 and 15, Corollaries 6 and 7 to 16 and 17, Theorems 8 and 9 to 18 and 19, Proposition 10 to Proposition S8, Lemma A.1 to Lemma S6, Theorem D.1 to Lemma S7. New: Theorem 7 (drift of the conditional mean, from the TH4 note `docs/theory_robustness.md`) and Corollary 10 (plug-in baseline through Theorem 7). Proposition A.6 and Lemma A.7 of the running-truncation note are not used in the paper; they are cited only in the companion note. Section 7 of the previous draft (the reproduction of the quarticity estimator) and Appendix D were removed.

## Length

Words are counted as whitespace-separated tokens after removing HTML comments and displayed equations, with each inline formula counted as one word; tables are included unless stated.

| Part | Words |
|---|---|
| Abstract | 196 |
| Main text, Section 1 to Reproducibility | 10,250 (9,409 without tables) |
| Introduction | 949 |
| References | about 1,170 |
| Supplement | about 10,700 |
| Companion note `extensions.md` | about 5,700 |

By section: Introduction 952, Section 2 765, Section 3 1,385, Section 4 1,033, Section 5 1,692, Section 6 2,474, Section 7 1,379, Section 8 455, Reproducibility 115 (these include the section headings). The open slots below will add to the main text; a figure of about 10,000 leaves little room, so filling them should move other material to the supplement.

## Provenance of numbers

Experiment labels (E1 to E21, ALT1 to ALT5, TH4, and their sub-labels such as Q4a) no longer appear in the text of the paper. Each paragraph, table caption or figure caption that reports results ends with a hidden HTML comment such as `<!-- provenance: E4 Q4a -->` naming the sections of `docs/experiment_report.md` behind it. Pandoc drops these comments, so they do not appear in the LaTeX or a PDF, but `scripts/trace_numbers.py` reads them to set the scope of each number. Supplement S9 maps every table and figure to its script and result file. The companion note keeps the labels in the text, since it is a technical note on the repository's experiments.

Both files were checked with

```bash
python scripts/trace_numbers.py --target paper/short_paper.md --out <scratch>/trace_paper.csv
python scripts/trace_numbers.py --target paper/extensions.md --out <scratch>/trace_ext.csv
```

with no untraced number in either. Captions say whether a table or figure shows simulated or market data.

## Open slots

HTML comments of the form `<!-- SLOT: ... -->` mark where pending results go. Each says what is expected.

| Slot | Where | What |
|---|---|---|
| TH5, universal speed limit | end of Section 3.3 | impossibility and speed limit for every stopping rule with $P(T \le n) \le n/A$, not only e-detectors |
| TH5, item 10 | end of Section 5.1 | replace the e-detector floors by the bound for every rule (45 steps for a complete collapse at $c = 3$, $A = 7800$) |
| TH5, joint regime | end of Section 4.1 | the plug-in, or a plug-in inflated by $z$ standard errors, as a function of calibration length, horizon and level |
| TH6, per-period guarantees | Section 3.1, after Theorem 4 | a construction with a proved bound per period and its cost |
| TH6, Theorem 11 | Section 4.2, after Theorem 11 | the algorithm-free tolerance hypothesis; the statement already uses it, TH6 may sharpen it |
| E19 | Section 7.1 | real-data null streams; the confidence-sequence baseline as the primary configuration |
| E20 | Section 7.2, after Table 6 | Holm within families, day-block bootstrap intervals, forecast regressions or QLIKE, a realised-variance outcome |
| E21 | Section 6.2 | robust CUSUM and Shiryaev-Roberts rules on the clipped statistic, a robust BOCPD, an econometric variance monitor |
| E18 | Section 7.3 | landmark market events |

TH4 is filled: Theorem 7 in Section 3.2, Corollary 10 in Section 4.1, proofs and the pathwise results in Supplement S1.2 and S1.3.

## References

Every entry of the reference lists was checked against a publisher, journal, arXiv or index listing, either by the referee (`docs/review/referee_report_paper.md`, Sections 2.2 and 2.3) or by our own web search on 3 November 2026. New in this version: Berkes et al. (2004), Chu, Stinchcombe and White (1996), Choe and Ramdas (2024), Clerico (2025), Dandapanthula and Ramdas (2025), Dette and Gösmann (2020), Fearnhead and Rigaill (2019), Horváth, Kokoszka and Zhang (2006), Huber (1965), Knoblauch, Jewson and Damoulas (2018), Mei (2008), Mineiro (2022), Shekhar and Ramdas (2024), Tartakovsky, Nikiforov and Basseville (2014), Unnikrishnan, Veeravalli and Meyn (2009), Vovk et al. (2021), Wang and Ramdas (2023) and Wang, Wang and Ziegel (2026). Horváth, Trapani and Wang is now cited as the *Journal of Time Series Analysis* article of 2025.

Not fully verified, to check before submission:

* Horváth, Trapani and Wang (2025): the article exists (IDEAS/RePEc listing, vol. 46(5), 981-996, doi 10.1111/jtsa.12824). Whether it is the published version of arXiv:2404.17885 ("Sequential monitoring for explosive volatility regimes", cited in the previous draft) is not confirmed; the search describes the arXiv paper as covering similar ground.
* Unnikrishnan, Veeravalli and Meyn: only the arXiv version (0911.2551) is confirmed; the journal version (often cited as *IEEE Transactions on Information Theory*, 2011) is not.
* Mei (2008): volume and issue confirmed (*Sequential Analysis* 27(4)); page range and DOI not found, so they are left out.
* Wang, Wang and Ziegel (2026): DOI confirmed; volume 72(6) and pages 4952-4973 come from a single secondary listing.
* Choe and Ramdas (2024): arXiv:2402.09698; its latest arXiv version says it was accepted by the *Journal of the Royal Statistical Society Series B*, but no published citation was found, so the arXiv reference is kept.
* Larsson, Ramdas and Ruf (2024): one listing says *Annals of Statistics* (2025), without volume or pages; the arXiv reference is kept.
* Lorden (1970): DOI 10.1214/aoms/1177697092 from a search summary that cites Project Euclid and Wikipedia; not opened directly.
* Inclán and Tiao (1994): DOI 10.1080/01621459.1994.10476824 from a software package citation, not from the publisher.
* Bahadur and Savage (1956): DOI 10.1214/aoms/1177728077 from a repository listing; volume 27 confirmed by several sources (one later paper cites volume 25).
* Mancini (2009): volume and pages confirmed; issue and DOI not found, so they are left out.
* Chu, Stinchcombe and White (1996): no DOI was confirmed (the JSTOR stable identifier suggests one, which we did not add).
* Waudby-Smith and Ramdas (2024): the paper is confirmed; the theorem number of the predictable plug-in empirical Bernstein bound could not be checked, so the paper cites it without a theorem number.
* Ram and Ramdas (2026): exists (arXiv:2602.05272); its text could not be opened, so the statement that Theorem 19 is an instance of it rests on the abstract.
* Bapat (2026): details as checked against the PDF in an earlier round; volume and DOI not confirmed online.

## Regenerating

From the repository root:

```bash
python paper/make_paper_figures.py          # writes the three PNG files in paper/figures/
cd paper
pandoc -s --shift-heading-level-by=-1 short_paper.md -o short_paper.tex
```

`--shift-heading-level-by=-1` turns the Markdown `##` headings into LaTeX `\section`. The section numbers are written in the headings, and the pandoc template switches LaTeX numbering off, so they are not doubled. The `.tex` was produced with pandoc 3.1.3. The PDF was not rebuilt in this round.

## Compiling the LaTeX elsewhere

With a TeX distribution that has the usual packages of the pandoc template (`amsmath`, `amssymb`, `graphicx`, `longtable`, `booktabs`, `hyperref`, `xcolor`), run from this folder:

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

For a journal template, copy the body of `short_paper.tex` (between `\begin{document}` and `\end{document}`) into the template, replace the hand-written reference list by `\bibliography{references}` with the template's author-year style (for example natbib), and turn the in-text citations into `\citet`/`\citep` commands. The keys in `references.bib` are of the form `surnameYEARword` (for example `shin2024edetectors`). The supplement is the part of `short_paper.md` after the heading "Supplementary material"; a journal will usually want it as a separate file.
