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

Main text: Abstract (no citations); 1 Introduction (three contributions and a related-work paragraph); 2 Model, statistic and nulls; 3 Validity (3.1 e-detectors and the guarantees, their reading over many days (Proposition 7) and the restart at every day boundary that bounds every day (Theorem 8); 3.2 a small drift of the conditional mean; 3.3 why a decrease detector needs a bounded statistic, and the pathwise floor for every valid rule); 4 Estimated baseline (4.1 plug-in, and when it keeps its level as $n$, $m$ and $A$ grow; 4.2 confidence sequence, with the null stated as a tolerance with a horizon, part 4 with no assumption, and what the tolerance excludes); 5 Detection delay (5.1 including the collapse floors and the delay floor for every valid rule); 6 Simulations; 7 Application (7.1 data, which results are confirmatory, and what real data can test; 7.2 five cryptocurrencies, with intervals from a block bootstrap over days, Holm adjustment and forecast regressions; 7.3 scheduled announcements and landmark shocks: the FOMC and CPI test, the one pre-specified primary test, and four case studies; 7.4 other data and lessons); 8 Discussion; Reproducibility; References.

Supplement: S1 proofs and details for Sections 2 and 3 (Corollary 6; Theorem 9 with its design numbers, Table S1 and the sharpness result; drift that varies; Corollary 11; the short proofs of Lemma 1, Theorem 3 and Propositions 5 and 10; S1.6 the pathwise floor, Theorem 12 and Corollary 13; S1.7 later days: proofs of Proposition 7 and Theorem 8, Proposition S6 on the delay of the restart, Theorem S7 on the CUSUM-type statistic, Proposition S8 on its failure under martingale differences, and the numbers); S2 proofs for Section 4 (including Theorem 15 and its numbers, and part 4 of Theorem 16); S3 proofs for Section 5 (including Theorem 23, Proposition S11 on the growth of the floor, the numbers for Theorem 23, and Proposition S12 on the decrease bet at the boundary); S4 second-order delay at practical levels (Table S2, Figure S1); S5 truncation constant, slack of the bound, bet design and further comparisons (Table S3; offline scan and cost; width of the confidence-sequence baseline); S6 clustering, forecasts, the GARCH bound and price noise (Proposition S13, Figure S2, Table S4); S7 real-data details (S7.1 Nifty 50 and India VIX; S7.2 FOMC and CPI; S7.3 competitors with Table S5; S7.4 the re-analysis of the real-data inference with Table S6; S7.5 the landmark shocks); S8 settings of every experiment (Tables S7 and S8); S9 where each table and figure comes from (script and result file).

Numbering: one counter for all results of the main text (Lemma 1 to Theorem 25); results stated only in the supplement are Lemma S1 to Proposition S13. Tables are 1 to 6 in the main text and S1 to S8 in the supplement. New with TH6 (`docs/theory_per_period.md`, referee M2 and M4): Proposition 7 (what Theorem 4 gives after the first day) and Theorem 8 (restart at every day boundary) in Section 3.1, Propositions S6 and S8 and Theorem S7 in Supplement S1.7; Theorem 16 is restated with the tolerance over a horizon $h$ and a part 4 with no assumption, as in `docs/math.md`, Theorem 6.4. Map from the previous version of this draft (commit 126588d): Lemma 1 to Corollary 6 and Lemma S1 to Theorem S5 are unchanged; every main-text result from Theorem 7 on moves up by 2 (Theorem 7 to 9, Proposition 8 to 10, Corollary 9 to 11, Theorem 10 to 12, Corollary 11 to 13, Corollary 12 to 14, Theorem 13 to 15, Theorem 14 to 16, Corollaries 15 and 16 to 17 and 18, Lemmas 17 and 18 to 19 and 20, Corollaries 19 and 20 to 21 and 22, Theorems 21, 22 and 23 to 23, 24 and 25); Lemmas S6 and S7 become S9 and S10, and Propositions S8, S9 and S10 become S11, S12 and S13. New with TH5 (`docs/theory_limits.md`): Theorem 12 (pathwise floor), Corollary 13 (no truncation; complete collapse), Theorem 15 (the plug-in as $n$, $m$ and $A$ grow), Theorem 23 (delay floor), Propositions S11 and S12. Map from the draft before the referee round (commit 8f71a2e) to the current numbers: Theorem A to Theorem 3, Theorem B to Theorem 4, Proposition 1 to Proposition 5, Corollary C to Corollary 6, Proposition 2 to Proposition 10, Corollary 3 to Corollary 11, Theorem D to Theorem 16, Corollaries E and F to 17 and 18, Lemmas 4 and 5 to Lemmas 19 and 20, Corollaries 6 and 7 to 21 and 22, Theorems 8 and 9 to 24 and 25, Proposition 10 to Proposition S13, Lemma A.1 to Lemma S9, Theorem D.1 to Lemma S10. Theorem 9 (drift of the conditional mean) and Corollary 14 (plug-in baseline through Theorem 9) come from the TH4 note `docs/theory_robustness.md`. Proposition A.6 and Lemma A.7 of the running-truncation note are not used in the paper; they are cited only in the companion note. Section 7 of the previous draft (the reproduction of the quarticity estimator) and Appendix D were removed. The theory notes and report parts name the paper's results by the numbers they had when they were written (for example `docs/report_parts/TH6.md` and `docs/math.md` call Theorem 16 "Theorem 14" or "Theorem D", and Theorem 4 "Theorem B").

## Length

Words are counted as whitespace-separated tokens after removing HTML comments and displayed equations, with each inline formula counted as one word; tables are included unless stated.

| Part | Words |
|---|---|
| Abstract | 197 |
| Main text, Section 1 to Reproducibility | 12,725 |
| Introduction | 1,038 |
| References | about 1,290 |
| Supplement | about 17,000 |
| Companion note `extensions.md` | about 6,300 |

By section: Introduction 1,038, Section 2 765, Section 3 2,338, Section 4 1,647, Section 5 2,102, Section 6 2,329, Section 7 1,812, Section 8 579, Reproducibility 115 (these include the section headings). TH5 added about 1,140 words to the main text, and TH6 with the corrections of the TH5 review about 1,030 more (Proposition 7, Theorem 8 and their numbers in Section 3.1; the tolerance with a horizon, part 4 and what the tolerance excludes in Section 4.2; short changes to the abstract, the contributions and the Discussion); the proofs and most numbers are in Supplement S1.6, S1.7, S2 and S3. The main text is now about 2,700 words above 10,000, and the open slots below (E19, E21) will add to it. We suggest one pass that moves material to the supplement once they are filled, rather than cutting piecemeal; candidates are the measured later-day numbers after Proposition 7 and the delay and CUSUM sentences after Theorem 8 (all already repeated in Supplement S1.7), the plug-in numbers of Section 4.1 after Corollary 14, the discussion after Theorem 15, "What the tolerance excludes" in Section 4.2, the paragraph of Section 5.1 that compares the collapse floors, and parts of Section 6.

## Provenance of numbers

Experiment labels (E1 to E21, ALT1 to ALT5, TH4 to TH6, and their sub-labels such as Q4a) no longer appear in the text of the paper. Each paragraph, table caption or figure caption that reports results ends with a hidden HTML comment such as `<!-- provenance: E4 Q4a -->` naming the sections of `docs/experiment_report.md` behind it. Pandoc drops these comments, so they do not appear in the LaTeX or a PDF, but `scripts/trace_numbers.py` reads them to set the scope of each number. Supplement S9 maps every table and figure to its script and result file. The companion note keeps the labels in the text, since it is a technical note on the repository's experiments.

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
| E19 | Section 7.1 | real-data null streams; the confidence-sequence baseline as the primary configuration |
| E21 | Section 6.2 | robust CUSUM and Shiryaev-Roberts rules on the clipped statistic, a robust BOCPD, an econometric variance monitor |

Filled: TH4 (Theorem 9 in Section 3.2, Corollary 14 in Section 4.1, proofs and the pathwise results in Supplement S1.2 and S1.3); E20 (Sections 7.1 to 7.4, Table 6, Supplement S7.1 to S7.4 with Table S6, and the restart, panel and frequency paragraphs of the companion note); E18 (Section 7.3, Supplement S7.5, and its exploratory run of the restarted monitors in the companion note, Section 1); TH5 (Theorem 12 and Corollary 13 at the end of Section 3.3, Theorem 15 at the end of Section 4.1, the floors and Theorem 23 at the end of Section 5.1, with proofs and numbers in Supplement S1.6, S2 and S3, Propositions S11 and S12, and a row of Table S7), corrected after the independent review of TH5 (`docs/review/independent_review_round7_th5.md`: the 45-step floor is for deterministic rules, with 0.98 and 43.4 for randomised ones; the delay ratios are upper bounds on the inefficiency against the best valid deterministic procedure, with the randomised ratios added; the break of the plug-in needs $\log A$ small compared with $m$; the worst case of the plug-in is about 33% to 42% from 20,000 windows; 0.72 to 0.86 at $A = 10^8$); TH6 (Proposition 7 and Theorem 8 in Section 3.1, the restated Theorem 16 in Section 4.2, Supplement S1.7, a row of Table S7 and the result files in S9), with the per-day readings of Theorem 4 corrected in the abstract, Sections 3.1, 6 and 7.4 and the companion note (referee M2).

## References

Every entry of the reference lists was checked against a publisher, journal, arXiv or index listing, either by the referee (`docs/review/referee_report_paper.md`, Sections 2.2 and 2.3) or by our own web search on 3 November 2026. New in this version: Berkes et al. (2004), Chu, Stinchcombe and White (1996), Choe and Ramdas (2024), Clerico (2025), Dandapanthula and Ramdas (2025), Dette and Gösmann (2020), Fearnhead and Rigaill (2019), Horváth, Kokoszka and Zhang (2006), Huber (1965), Knoblauch, Jewson and Damoulas (2018), Mei (2008), Mineiro (2022), Shekhar and Ramdas (2024), Tartakovsky, Nikiforov and Basseville (2014), Unnikrishnan, Veeravalli and Meyn (2009), Vovk et al. (2021), Wang and Ramdas (2023) and Wang, Wang and Ziegel (2026). Horváth, Trapani and Wang is now cited as the *Journal of Time Series Analysis* article of 2025. Added with the re-analysis of the real-data inference: Holm (1979), Driscoll and Kraay (1998), Patton (2011) and Diebold and Mariano (1995), checked by our own web search on 3 November 2026 against index and library listings. Added with TH5: Dvoretzky, Kiefer and Wolfowitz (1956) and Massart (1990), used in the proof of Theorem 15, checked by our own web search on 3 November 2026 against journal and index listings.

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
* Ram and Ramdas (2026): exists (arXiv:2602.05272); its text could not be opened, so the statement that Theorem 25 is an instance of it rests on the abstract.
* Bapat (2026): details as checked against the PDF in an earlier round; volume and DOI not confirmed online.
* Bibinger, Jirak and Vetter (2017): the reference is confirmed, but we could not obtain the paper. The two-window variance-ratio scan of Sections 6.2 and 7.2 is therefore called "inspired by" their work and is not presented as their procedure (it was called a BJV-type test in earlier drafts and in the experiment reports).
* Holm (1979): volume, issue and pages confirmed by several listings; the JSTOR stable identifier (10.2307/4615733, from a software package's documentation) was not added as a DOI.
* Driscoll and Kraay (1998): DOI 10.1162/003465398557825 from two secondary listings, not from the publisher.
* Diebold and Mariano (1995): volume, issue and pages confirmed; no DOI was confirmed, so none is given.
* Clerico (2025): cited for a related optimality of betting for a bounded mean, next to Theorem 12. The TH5 note checked it from search listings only; the paper itself was not opened.
* Waudby-Smith, Sandoval and Jordan (2025, arXiv:2504.02818) and Kesten (1973), cited in the TH5 note, are not cited in the paper: the first was checked from search listings only, and the second supports only a heuristic that the paper states without a reference.
* Politis and Romano (1992), the circular block bootstrap, in *Exploring the Limits of Bootstrap* (Wiley): not cited, because neither the chapter title nor the pages could be confirmed; the paper describes the resampling scheme without a citation.

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
