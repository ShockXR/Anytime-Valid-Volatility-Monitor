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

Main text: Abstract (no citations); 1 Introduction (three contributions and a related-work paragraph); 2 Model, statistic and nulls; 3 Validity (3.1 e-detectors and the guarantees, their reading over many days (Proposition 7) and the restart at every day boundary that bounds every day with a known baseline (Theorem 8); 3.2 a small drift of the conditional mean; 3.3 why a decrease detector needs a bounded statistic, and the pathwise floor for every valid rule); 4 Estimated baseline (4.1 plug-in, and when it keeps its level as $n$, $m$ and $A$ grow; 4.2 confidence sequence, with the null stated as a tolerance with a horizon, part 4 with no assumption, and what the tolerance excludes); 5 Detection delay (5.1 including the collapse floors and the delay floor for every valid rule; 5.2 delay bounds; 5.3 practical levels and the price of robustness); 6 Simulations (6.1 validity; 6.2 delay against parametric, robust and calibrated classical competitors; 6.3 the cost of an estimated baseline; 6.4 when the null is false, which now holds the former Section 5.4 on clustering and forecast-relative monitoring; 6.5 extensions); 7 Application (7.1 data, which results are confirmatory, and null streams built from real returns; 7.2 five cryptocurrencies, with the confidence-sequence baseline as the pre-registered primary configuration next to the plug-in, forecasts and a rule tuned with hindsight, forecast regressions and competitors at matched alarm counts; 7.3 scheduled announcements and landmark shocks; 7.4 other data and lessons); 8 Discussion; Reproducibility; References. The main text has two tables (Table 1, the price of robustness; Table 2, false alarms of the simulated detectors); the other tables of the previous draft are now in the supplement.

Supplement: S1 proofs and details for Sections 2 and 3 (S1.1 Corollary 6 with its numbers; S1.2 Theorem 9 with its design numbers, Table S1 and the sharpness result; S1.3 drift that varies; S1.4 Corollary 11; S1.5 the short proofs of Lemma 1, Theorem 3 and Propositions 5 and 10; S1.6 the pathwise floor, Theorem 12 and Corollary 13; S1.7 later days: proofs of Proposition 7 and Theorem 8, the window across a boundary, estimated baselines, Proposition S6 on the delay of the restart, Theorem S7 on the CUSUM-type statistic, Proposition S8 on its failure under martingale differences, and the numbers; S1.8 fourth powers and zero returns); S2 proofs for Section 4 (the confidence-sequence bounds, Lemmas S9 and S10, Theorem 16 and its corollaries, the bet-weighted condition, Theorem 15 with the corrected distribution-free constant, the plug-in simulations and the numbers); S3 proofs for Section 5 (including Theorem 23, Proposition S11, the numbers for Theorem 23 and Proposition S12); S4 second-order delay at practical levels (Table S2, Figure S1); S5 details for Section 6 (truncation constant and split clip, with the confidence-sequence baseline for the Gaussian law, Table S3; slack of the bound and the false alarm rates with the true baseline, Table S4; delays in the designs of Table 2, Table S5; BOCPD and the scan; the competitors calibrated for the data at hand, with their definitions, false alarms and delays; the Student $t$ BOCPD; width of the confidence-sequence baseline and its delays, Table S6); S6 clustering, forecasts, the GARCH bound and price noise (Proposition S13, Figure S2, Table S7); S7 real-data details (S7.1 Nifty 50 and India VIX; S7.2 FOMC and CPI; S7.3 forecasts and competitors, with Tables S8 and S9 and the classical rules on real data; S7.4 the re-analysis of the real-data inference with Table S10; S7.5 the landmark shocks; S7.6 null streams from real returns, the confidence-sequence baseline on the five coins, blind segments and the forecast-relative null); S8 settings of every experiment (Tables S11 and S12); S9 where each table and figure comes from (script and result file).

Numbering: one counter for all results of the main text (Lemma 1 to Theorem 25); results stated only in the supplement are Lemma S1 to Proposition S13. Tables are 1 and 2 in the main text and S1 to S12 in the supplement. Map of the tables from the previous version: Tables 1 and 2 are unchanged; Table 3 is now Table S4, Table 4 is S5, Table 5 is S6 and Table 6 is S8; Tables S1 to S3 are unchanged; Table S4 is now S7, S5 is S9, S6 is S10, S7 is S11 and S8 is S12. The result numbers (Lemma 1 to Theorem 25, Lemma S1 to Proposition S13) did not change in this round. New with TH6 (`docs/theory_per_period.md`, referee M2 and M4): Proposition 7 (what Theorem 4 gives after the first day) and Theorem 8 (restart at every day boundary) in Section 3.1, Propositions S6 and S8 and Theorem S7 in Supplement S1.7; Theorem 16 is restated with the tolerance over a horizon $h$ and a part 4 with no assumption, as in `docs/math.md`, Theorem 6.4. Map from the previous version of this draft (commit 126588d): Lemma 1 to Corollary 6 and Lemma S1 to Theorem S5 are unchanged; every main-text result from Theorem 7 on moves up by 2 (Theorem 7 to 9, Proposition 8 to 10, Corollary 9 to 11, Theorem 10 to 12, Corollary 11 to 13, Corollary 12 to 14, Theorem 13 to 15, Theorem 14 to 16, Corollaries 15 and 16 to 17 and 18, Lemmas 17 and 18 to 19 and 20, Corollaries 19 and 20 to 21 and 22, Theorems 21, 22 and 23 to 23, 24 and 25); Lemmas S6 and S7 become S9 and S10, and Propositions S8, S9 and S10 become S11, S12 and S13. New with TH5 (`docs/theory_limits.md`): Theorem 12 (pathwise floor), Corollary 13 (no truncation; complete collapse), Theorem 15 (the plug-in as $n$, $m$ and $A$ grow), Theorem 23 (delay floor), Propositions S11 and S12. Map from the draft before the referee round (commit 8f71a2e) to the current numbers: Theorem A to Theorem 3, Theorem B to Theorem 4, Proposition 1 to Proposition 5, Corollary C to Corollary 6, Proposition 2 to Proposition 10, Corollary 3 to Corollary 11, Theorem D to Theorem 16, Corollaries E and F to 17 and 18, Lemmas 4 and 5 to Lemmas 19 and 20, Corollaries 6 and 7 to 21 and 22, Theorems 8 and 9 to 24 and 25, Proposition 10 to Proposition S13, Lemma A.1 to Lemma S9, Theorem D.1 to Lemma S10. Theorem 9 (drift of the conditional mean) and Corollary 14 (plug-in baseline through Theorem 9) come from the TH4 note `docs/theory_robustness.md`. Proposition A.6 and Lemma A.7 of the running-truncation note are not used in the paper; they are cited only in the companion note. Section 7 of the previous draft (the reproduction of the quarticity estimator) and Appendix D were removed. The theory notes and report parts name the paper's results by the numbers they had when they were written (for example `docs/report_parts/TH6.md` and `docs/math.md` call Theorem 16 "Theorem 14" or "Theorem D", and Theorem 4 "Theorem B").

## Length

Words are counted as whitespace-separated tokens after removing HTML comments and displayed equations, with each inline formula counted as one word; tables are included unless stated.

| Part | Words |
|---|---|
| Abstract | 216 |
| Main text, Section 1 to Reproducibility | 10,225 |
| Introduction | 837 |
| References | about 1,340 |
| Supplement | about 21,900 |
| Companion note `extensions.md` | about 6,300 |

By section: Introduction 837, Section 2 621, Section 3 1,903, Section 4 1,347, Section 5 1,517, Section 6 1,632, Section 7 1,776, Section 8 510, Reproducibility 82 (these include the section headings). The length pass of this round brought the main text from 12,725 words to 10,225 while adding the results of E19 and E21: Tables 3 to 6 moved to the supplement (as Tables S4, S5, S6 and S8), the former Section 5.4 was merged into Section 6.4, the formulas of the confidence-sequence bounds moved to Supplement S2, and the discussions after most results were shortened; every number cut from the main text is now in the supplement (S1.1, S1.2, S1.6, S1.7, S1.8, S2, S3, S5, S7.3, S7.6) or, for the exploratory restarts, in the companion note. A further cut of about 225 words would need to move a result or a table out of the main text.

## Provenance of numbers

Experiment labels (E1 to E21, ALT1 to ALT5, TH4 to TH6, and their sub-labels such as Q4a) no longer appear in the text of the paper. Each paragraph, table caption or figure caption that reports results ends with a hidden HTML comment such as `<!-- provenance: E4 Q4a -->` naming the sections of `docs/experiment_report.md` behind it. Pandoc drops these comments, so they do not appear in the LaTeX or a PDF, but `scripts/trace_numbers.py` reads them to set the scope of each number. Supplement S9 maps every table and figure to its script and result file. The companion note keeps the labels in the text, since it is a technical note on the repository's experiments.

Both files were checked with

```bash
python scripts/trace_numbers.py --target paper/short_paper.md --out <scratch>/trace_paper.csv
python scripts/trace_numbers.py --target paper/extensions.md --out <scratch>/trace_ext.csv
```

with no untraced number in either. Captions say whether a table or figure shows simulated or market data.

## Open slots

There are no open slots. HTML comments of the form `<!-- SLOT: ... -->` would mark where pending results go.

Filled: TH4 (Theorem 9 in Section 3.2, Corollary 14 in Section 4.1, proofs and the pathwise results in Supplement S1.2 and S1.3); E20 (Sections 7.1 to 7.4, Supplement S7.1 to S7.4 with Tables S8 and S10, and the restart, panel and frequency paragraphs of the companion note); E18 (Section 7.3, Supplement S7.5, and its exploratory run of the restarted monitors in the companion note, Section 1); TH5 (Theorem 12 and Corollary 13 at the end of Section 3.3, Theorem 15 at the end of Section 4.1, the floors and Theorem 23 at the end of Section 5.1, with proofs and numbers in Supplement S1.6, S2 and S3, Propositions S11 and S12, and a row of Table S11), corrected after the independent review of TH5 (`docs/review/independent_review_round7_th5.md`); TH6 (Proposition 7 and Theorem 8 in Section 3.1, the restated Theorem 16 in Section 4.2, Supplement S1.7, a row of Table S11 and the result files in S9), with the per-day readings of Theorem 4 corrected (referee M2).

New in this round:

* E21 (referee M5): the paragraph "Competitors calibrated for the data at hand" in Section 6.2 [simulated], the paragraph "Competitors at matched alarm counts" in Section 7.2 [market data], the definitions, false alarms, delays and the Student $t$ BOCPD in Supplement S5, the classical rules on real data in Supplement S7.3, rows of Tables S11 and S12, and the result files in S9. The scan is called "two-window variance-ratio scan inspired by Bibinger, Jirak and Vetter (2017)" everywhere. Not done, because the papers could not be read: the exact procedures of Horváth, Kokoszka and Zhang (2006) and Chu, Stinchcombe and White (1996), and a $\beta$-divergence BOCPD.
* E19 (referee M1): the paragraph "Null streams built from real returns" in Section 7.1 [constructed streams with the real marginal law], the paragraph "The confidence-sequence baseline" in Section 7.2 [market data], the limitations in Section 8, Supplement S7.6, the confidence-sequence run of the split clip in Supplement S5 (m32; the Gaussian row of Table S3), rows of Tables S11 and S12, and the result files in S9. The confidence-sequence baseline is reported as the pre-registered primary configuration next to the plug-in; the forecast comparisons, regressions and competitor analyses of Section 7.2 are still those of the plug-in. Post hoc items are named as such (day standardisation of the bootstrap streams, the blind segments).
* The TH5 correction of the distribution-free constant: $C = 1.79$ (not 2.52) in Theorem 15(a) and Supplement S2, with $C = k_1 I_0/4$, 48 times the Gaussian constant (not 68), and its bound informative only for $m$ far below $\sqrt n$.
* The TH6 corrections of the review of round 8 (`docs/review/independent_review_round8_th6.md`): Theorem 8 is stated for a known baseline; with the confidence-sequence baseline the restart gives $\delta + m/A$ per day only marginally and while the tolerance holds to the end of the day, and with the plug-in nothing per day is proved; Theorem S7 is an unconditional bound; a window across one boundary has $(2m + \kappa)/A$; the sharpness of $dm/A$ is no longer claimed (35.7% against 50% on day 10, 73.2% against 100% on day 20); the parking stream passes the bound on day 16 (5.30%) with 1440-step days; the GARCH bands are per path, while parts 1 to 3 of Theorem 16 need the tolerance almost surely; bold play is 5.23% pooled (standard error about 0.06) against 5.25% as the average of the daily hazards; the restarted detector's largest probability is 2.66% (390) and its largest hazard 2.75%. The paper had no 4.6% figure for Mei's criterion; it says 4.7%.

## References

Every entry of the reference lists was checked against a publisher, journal, arXiv or index listing, either by the referee (`docs/review/referee_report_paper.md`, Sections 2.2 and 2.3) or by our own web search on 3 November 2026. New in this version: Berkes et al. (2004), Chu, Stinchcombe and White (1996), Choe and Ramdas (2024), Clerico (2025), Dandapanthula and Ramdas (2025), Dette and Gösmann (2020), Fearnhead and Rigaill (2019), Horváth, Kokoszka and Zhang (2006), Huber (1965), Knoblauch, Jewson and Damoulas (2018), Mei (2008), Mineiro (2022), Shekhar and Ramdas (2024), Tartakovsky, Nikiforov and Basseville (2014), Unnikrishnan, Veeravalli and Meyn (2009), Vovk et al. (2021), Wang and Ramdas (2023) and Wang, Wang and Ziegel (2026). Horváth, Trapani and Wang is now cited as the *Journal of Time Series Analysis* article of 2025. Added with the re-analysis of the real-data inference: Holm (1979), Driscoll and Kraay (1998), Patton (2011) and Diebold and Mariano (1995), checked by our own web search on 3 November 2026 against index and library listings. Added with TH5: Dvoretzky, Kiefer and Wolfowitz (1956) and Massart (1990), used in the proof of Theorem 15, checked by our own web search on 3 November 2026 against journal and index listings. Added with E21: Horváth, Hušková, Kokoszka and Steinebach (2004), for the boundary of the training-sample rules, and Zeileis, Leisch, Hornik and Kleiber (2002), the R package strucchange, whose OLS-CUSUM monitoring boundary was used as a check.

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
* Horváth, Hušková, Kokoszka and Steinebach (2004): journal, volume, year and pages confirmed by listings and by the E21 report; the issue (given as 1 in the E21 report) and the DOI were not confirmed, so they are left out. The paper itself could not be read; the boundary was taken from search excerpts and later papers.
* Zeileis, Leisch, Hornik and Kleiber (2002): the package paper is cited for strucchange. The monitoring functions of the package may have been described in a later paper, which we did not check, and the package version used by E21 (1.6-0, read from a GitHub mirror of CRAN) differs from the version that CRAN listed (1.5-4); no version is cited in the paper.
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
