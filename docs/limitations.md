# Limitations

What the prototype does not establish, in one place. Each item says how large the problem is, what evidence we have, and what would fix it. Items are ordered by how much they matter for the claims made in `README.md` and `PITCH_BAPAT.md`.

## 1. What the false alarm guarantee covers

* **The null is a statement about conditional means.** Theorems A and B (`docs/math.md`, 5.4) hold when the conditional mean of the truncated squared return stays at (or below, for the increase side) its baseline. Under persistent volatility clustering this is false even when the unconditional variance never changes, and the monitor then alarms on most days (E9, Q9c: 12% to 20% of days even with a band of 1 when the GARCH half-life is an hour or more). These alarms are real changes in conditional volatility, but a user who only cares about structural breaks will count them as false. Remedy: monitor returns divided by a predictable volatility forecast (E10), which restores validity at the cost of power against breaks the forecast absorbs.
* **The baseline is estimated.** With a plug-in baseline the bound holds only approximately (E4d; E10a found about 5% per day with a one-day calibration window against a bound that is about 1% with the true baseline). The confidence-sequence baseline (E8, `docs/math.md` section 6) makes the bound `delta + m / A` hold with an estimated baseline, under the stated assumptions.
* **Recalibration after an alarm.** Each monitored segment has its own guarantee. Statements about "alarms per day" over a long series are averages over segments, not a single anytime-valid statement for the whole series.
* **Several assets.** Running the monitor on five coins with a 5% budget each gives up to about 25% chance of at least one false alarm per day across the group (union bound). A family-wise statement needs the budget split across assets.
* **The bound is conservative.** At one day it is about 4 times larger than the observed rate with the true baseline (E9). This costs detection speed; a sharper threshold is studied in E11 but the default keeps the proven level.

## 2. Real-data validation has no ground truth

* Whether a volatility regime "really" changed at a given minute is not observable. The real-data checks (E5, E6, E7, E10) are indirect: realised volatility after the alarm, comparison with forecasts, India VIX, and scheduled events. Volatility clusters, so any detector that reacts to a burst will look somewhat predictive; the random-time and forecast controls reduce but do not remove this.
* Data come from public mirrors (Binance klines via a third-party GitHub archive, a cleaned Nifty 50 archive, a Bitstamp archive), not from the exchanges directly. The loaders report gaps, duplicates and bad prices, and `data/MANIFEST.csv` pins the exact bytes used.
* The paper does not name its crypto venue, so our reproductions of Tables 5 and 6 compare against a likely different data set (E6 discusses which features match).

## 3. Estimation

* Truncation does not remove microstructure noise (`docs/estimator_issues.md`, 4b). At one-minute sampling the noise is small for the assets used here; at higher frequencies a noise-robust first stage would be needed (alternative `alternatives/preaveraging/`).
* Clipping keeps τ⁴ for each large jump, so very frequent large jumps still bias truncated quarticity upward (E2c). Dropping instead of clipping (threshold RQ) removes more of this.
* The κ₄ bias correction is exact only for locally Gaussian diffusive returns.
* Our reading of how the paper's tables were computed (unnormalised weights, raw MAD, no factor n in RQ) is an inference from the published numbers. It needs the author's code to confirm, and is presented as questions, not corrections.

## 4. Scope

* One-minute bars only. Tick data, other venues and other asset classes are untested.
* The risk switch and correlation monitor are minimal prototypes, kept in `alternatives/`, and their results are not claims about trading performance.
