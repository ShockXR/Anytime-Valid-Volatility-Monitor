# Anytime-Valid Volatility Monitor

[![tests](https://github.com/ShockXR/Anytime-Valid-Volatility-Monitor/actions/workflows/tests.yml/badge.svg)](https://github.com/ShockXR/Anytime-Valid-Volatility-Monitor/actions/workflows/tests.yml)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776ab)
[![License: MIT](https://img.shields.io/badge/license-MIT-5b6475)](LICENSE)
![Status: research prototype](https://img.shields.io/badge/status-research%20prototype-1f4e8c)

**Can we detect a change in volatility as one-minute returns arrive, while limiting the influence of isolated price jumps?**

This prototype squares each return and caps it at a cutoff learned during calibration. Two sequential detectors then accumulate evidence that the average of this capped quantity has increased or decreased. After an alarm, the default monitor pauses to learn a new baseline.

The clipping idea is drawn from [Bapat's quarticity estimator (2026)](https://doi.org/10.1016/j.iref.2026.105482); the detection framework comes from [Shin, Ramdas and Rinaldo's e-detectors (2024)](https://doi.org/10.51387/23-NEJSDS51). This project explores their combination through implementation, calibration comparisons, simulations and market replays. The main monitor tracks **clipped squared returns**, not integrated quarticity itself.

**Research status:** preliminary. The default estimated baseline does not have the general false-alarm guarantee available under the known-baseline assumptions. The mathematical arguments and code have not received human expert review. See [scope](#scope-and-guarantees) and [development provenance](#development-provenance).

## Start here

No installation is needed to read the documents or the executed demonstration.

| To investigate…                                   | Start with…                                                                                                                                              |
| ------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| What is being measured and what “no change” means | [Statistic and null hypothesis](docs/math.md#51-the-statistic-and-the-null)                                                                              |
| The false-alarm argument and its conditions       | [Guarantees](docs/math.md#54-e-detectors-and-the-false-alarm-guarantees) and [what they do not cover](docs/math.md#56-what-is-not-covered-by-the-proofs) |
| A baseline learned from data                      | [Estimated-baseline formulation](docs/math.md#6-keeping-the-guarantee-with-an-estimated-baseline)                                                        |
| Simulation evidence                               | [Detection, false alarms and delay](docs/experiment_report.md#e4-detecting-changes-false-alarms-and-delay)                                               |
| Later checks and negative findings                | [Latest validation summary](reports/round4/ROUND4_SUMMARY.md) and [matched-control / forecasting checks](reports/round4/W1/REPORT.md)                    |
| Practical and statistical limitations             | [Limitations](docs/limitations.md)                                                                                                                       |
| A demonstration without running code              | [Executed notebook](notebooks/demo.ipynb)                                                                                                                |
| The implementation                                | [Streaming monitor](src/avmonitor/monitor.py), [detectors](src/avmonitor/detectors.py), [confidence sequences](src/avmonitor/confseq.py)                 |
| The longer research formulation                   | [Working paper and supplement](paper/short_paper.md), not an expert-reviewed or submission-ready manuscript                                              |
| Sources and their relevance                       | [Annotated literature review](docs/literature_review.md)                                                                                                 |

## Scope and guarantees

An alarm indicates evidence of a change in the **conditional mean of the capped squared return**. It does not by itself identify a change in the underlying variance: changes in the tails or jump distribution can affect this statistic too. Nor is an alarm a forecast that volatility will continue rising.

| Configuration or claim       | What the current mathematical notes establish, subject to their assumptions                                                                                                                                                    |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Known baseline               | If the conditional mean stays at the specified baseline, the probability of any false alarm within the first `m` monitored steps after a fresh start is at most `m/A`. The two-sided monitor uses threshold `2A` on each side. |
| Default plug-in baseline     | The baseline is estimated from a calibration window. The known-baseline guarantee does **not** automatically apply; estimation error can increase false alarms.                                                                |
| Confidence-sequence baseline | The notes give a bound of `delta + m/A`, including calibration uncertainty, under the stated conditions linking calibration and monitoring. This does not remove the need for assumptions about the data.                      |
| Repeated days and restarts   | A fresh-start bound is not a bound for every later day. Daily restarts have separate results and conditions.                                                                                                                   |

For illustration, `m = 390` monitored steps and `A = 7,800` give a fresh-start bound of 5% under the known-baseline null. This is a numerical illustration, not a certified 5% rate for the default market replay; Nifty sessions and data-cleaning rules use different step counts.

Details: [mathematical notes, sections 5–6](docs/math.md#5-anytime-valid-detection-of-variance-changes), [period-by-period guarantees](docs/theory_per_period.md), and [limitations](docs/limitations.md#1-what-the-false-alarm-guarantee-covers).

## Selected evidence

### Simulations with jumps

In experiment E4, with rare jumps and no regime change, the monitor's increase side raised false alarms on **2.3%** of simulated paths, compared with **70.7%** for a Gaussian likelihood-ratio e-detector and **70.3%** for Page's CUSUM calibrated on normal returns. Settings: 2,000 calibration steps, `A = 7,800`, a 390-step evaluation horizon and 300 paths. These are measured results with an estimated baseline, not a general guarantee.

When the standard deviation rose by 50% at step 500, the capped monitor's median delay was 36 steps. This comparison illustrates failure of Gaussian-calibrated methods under that jump model; it does not establish superiority over every classical detector. Later comparisons include classical rules on the same clipped statistic and matched false-alarm rates.

Sources: [E4 tables and settings](docs/experiment_report.md#e4-detecting-changes-false-alarms-and-delay), [later benchmark checks](reports/round4/ROUND4_SUMMARY.md).

### Market data: evidence of a current burst, not forecasting superiority

Across 27 asset-years of crypto and Nifty data, 91.4% of up alarms were followed by an hour above the calibrated level, compared with 74.6% of control minutes matched on previous-hour volatility. Matching on the last 10 minutes as well, in a later exploratory check, raised the control rate to 83.5%, narrowing the difference. On fresh 2026 data, adding alarms did not improve a HAR forecast already using the last hour for BTC, ETH or Nifty 50.

The research focus is false-alarm control under explicit assumptions; these market replays do not establish forecasting gains. See [the matched-control and forecasting report](reports/round4/W1/REPORT.md), which also documents the weaker evidence for decrease alarms.

### One historical example

<p align="center"><a href="docs/figures/nifty_2024-06-04_alarm.png"><img src="docs/figures/nifty_2024-06-04_alarm.png" width="560" alt="Nifty 50 increase-side evidence crosses the alarm threshold at 09:26 IST on 4 June 2024, followed by a period of recalibration without monitoring."></a></p>

On **4 June 2024**, the monitor raised an up alarm at **09:26 IST**, 11 minutes after the Nifty open. This is one historical case (E18), using the default plug-in baseline, which the known-baseline guarantee does not cover. The threshold here is `2A = 15,000`. Time since the open is not a measured detection delay from a known volatility changepoint.

After this alarm, the monitor spends 750 minutes—two sessions—recalibrating without monitoring. See the [two-panel figure](docs/figures/nifty_2024-06-04_alarm_detail.png), [executed notebook](notebooks/demo.ipynb), or download [the interactive replay](docs/demo/replay.html) and open it in a browser.

## Questions for further investigation

1. **Learning the baseline:** which assumptions are defensible for market data, and how much calibration uncertainty must the monitor allow?
2. **Interpreting the target:** how should changes in the clipped statistic be distinguished from changes in tails, jumps or microstructure noise?
3. **Keeping the monitor useful:** can calibration uncertainty and time spent recalibrating be reduced while retaining clearly stated guarantees?

The [limitations](docs/limitations.md) and [latest validation summary](reports/round4/ROUND4_SUMMARY.md) describe the existing evidence behind these questions. Further work should first select a focused question and explicit assumptions, then evaluate the proposed method against the existing benchmarks.

## Run and reproduce

<details>
<summary><b>Installation and a simulated demonstration</b></summary>

Requires Python 3.10 or newer. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
# Windows: .venv\Scripts\activate
pip install -e ".[dev]"
python -m avmonitor simulate --ratio 2
```

The example stream contains rare jumps; its standard deviation doubles at step 3000 and returns to its original level at step 6000. The default monitor recalibrates for 1,000 steps after each alarm.

The CLI prints a nominal false-alarm bound. Interpret that statement using the assumptions in [Scope and guarantees](#scope-and-guarantees): it is not a proved bound for the default plug-in calibration.

</details>

<details>
<summary><b>Tests, result tracing and experiment reproduction</b></summary>

```bash
# Software tests; no market data needed
python -m pytest -q

# Match numerical statements to committed result records
python scripts/trace_numbers.py --target README.md

# Data hashes, tests and experiment reruns
python scripts/reproduce.py
```

Software tests run through GitHub Actions. Passing tests does not establish that the mathematical assumptions hold on market data. Number tracing checks for matching records; it does not independently validate the numbers or their interpretation.

Raw market data are not committed. See [data sources and download instructions](data/README.md); `python scripts/download_data.py` retrieves data from public sources and records hashes. Reproducing market experiments requires those data and additional computation.

</details>

## Repository guide

| Location                           | Contents                                                                 |
| ---------------------------------- | ------------------------------------------------------------------------ |
| [src/avmonitor/](src/avmonitor/)   | Package implementation                                                   |
| [tests/](tests/)                   | Software and statistical simulation checks                               |
| [experiments/](experiments/)       | Experiment scripts                                                       |
| [results/](results/)               | Committed result records                                                 |
| [docs/](docs/)                     | Mathematical notes, experiment report, limitations and literature        |
| [reports/round4/](reports/round4/) | Later validation reports; begin with the summary                         |
| [paper/](paper/)                   | Working manuscript, supplement and extensions                            |
| [alternatives/](alternatives/)     | Exploratory alternative prototypes; not all are part of the main monitor |

[Estimator replication questions](docs/estimator_issues.md) are documented separately from the monitoring results. They include interpretations that require confirmation from the original author's implementation.

For the current overview, use this README together with the limitations and latest validation summary. [PITCH\_BAPAT.md](PITCH_BAPAT.md) is an earlier discussion brief and [README\_FULL.md](README_FULL.md) preserves the previous, fuller presentation; some framing may precede later qualifications. The working paper and mathematical notes remain material for scrutiny, not independent confirmation of the claims.

## Development provenance

Shreyan Thanki directed the project. The code, experiments, mathematical arguments and write-ups were produced with extensive AI assistance (Claude, by Anthropic). Reviews in [docs/review/](docs/review/) were generated by separate AI agents instructed to act as sceptical referees; they are **not human peer review**. No statistician has yet checked the proofs. The code is released under the [MIT license](LICENSE).

## References and supporting reading

**Direct sources and benchmark**

* Bapat, S. R. (2026). *A jump and noise robust kernel based estimator of quarticity for high frequency data.* International Review of Economics and Finance, 109, 105482. [DOI](https://doi.org/10.1016/j.iref.2026.105482). Source of the clipping idea.

* Shin, J., Ramdas, A. and Rinaldo, A. (2024). *E-detectors: a nonparametric framework for sequential change detection.* New England Journal of Statistics in Data Science, 2(2), 229–260. [DOI](https://doi.org/10.51387/23-NEJSDS51). Detection framework used here.

* Page, E. S. (1954). *Continuous inspection schemes.* Biometrika, 41(1/2), 100–115. [DOI](https://doi.org/10.1093/biomet/41.1-2.100). Source of the CUSUM benchmark.

**Background and possible directions**

The [annotated literature review](docs/literature_review.md#7-change-detection-and-testing-with-e-values-anytime-valid-inference) includes SAVI (Ramdas, Grünwald, Vovk and Shafer, 2023), confidence-sequence methods (Howard et al., 2021; Waudby-Smith and Ramdas, 2024), and change detection through sequential estimation (Shekhar and Ramdas, 2024). These support understanding of the theory and investigation of calibration alternatives; citing them does not certify this implementation or establish that an extension is valid. The review marks bibliographic details and full texts that remain unverified.
