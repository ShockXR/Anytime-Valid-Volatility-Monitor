## Null tables [REAL marginal, constructed streams]

How to read this: each row is one monitor setting run on permuted real calm-day returns, where no alarm is correct. 'Rate day 1' is the share of streams that alarmed within the first monitored day; 'Rate 5 days' within the first five. A rate is 'above' its bound when the lower end of its Wilson interval is above the bound (E19's rule), and 'holds' otherwise.

**five coins** (40,000 streams per setting)

| Setting | A | Alarms in day 1 | Rate day 1 [Wilson 95%] | Bound | Alarms in 5 days | Rate 5 days [Wilson 95%] | Bound | Verdict |
|---|---|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 28,800 | 1 | 0.00% [0.00%, 0.01%] | 5.00% | 9 | 0.02% [0.01%, 0.04%] | 25.00% | holds |
| plug-in, band 0.1 | 28,800 | 106 | 0.27% [0.22%, 0.32%] | 5.00% | 563 | 1.41% [1.30%, 1.53%] | 25.00% | holds |
| plug-in, band 0.05 | 28,800 | 457 | 1.14% [1.04%, 1.25%] | 5.00% | 3494 | 8.73% [8.46%, 9.02%] | 25.00% | holds |
| plug-in, band 0 | 28,800 | 2245 | 5.61% [5.39%, 5.84%] | 5.00% | 16777 | 41.94% [41.46%, 42.43%] | 25.00% | day 1 above, 5 days above |
| cs baseline, band 0.1 | 36,000 | 0 | 0.00% [0.00%, 0.01%] | 5.00% | 3 | 0.01% [0.00%, 0.02%] | 21.00% | holds |
| cs baseline, band 0 | 36,000 | 10 | 0.03% [0.01%, 0.05%] | 5.00% | 51 | 0.13% [0.10%, 0.17%] | 21.00% | holds |

**BTC** (8,000 streams per setting)

| Setting | A | Alarms in day 1 | Rate day 1 [Wilson 95%] | Bound | Alarms in 5 days | Rate 5 days [Wilson 95%] | Bound | Verdict |
|---|---|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 28,800 | 1 | 0.01% [0.00%, 0.07%] | 5.00% | 2 | 0.03% [0.01%, 0.09%] | 25.00% | holds |
| plug-in, band 0.1 | 28,800 | 29 | 0.36% [0.25%, 0.52%] | 5.00% | 110 | 1.38% [1.14%, 1.65%] | 25.00% | holds |
| plug-in, band 0.05 | 28,800 | 99 | 1.24% [1.02%, 1.50%] | 5.00% | 684 | 8.55% [7.96%, 9.18%] | 25.00% | holds |
| plug-in, band 0 | 28,800 | 437 | 5.46% [4.99%, 5.98%] | 5.00% | 3355 | 41.94% [40.86%, 43.02%] | 25.00% | day 1 holds, 5 days above |
| cs baseline, band 0.1 | 36,000 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 1 | 0.01% [0.00%, 0.07%] | 21.00% | holds |
| cs baseline, band 0 | 36,000 | 4 | 0.05% [0.02%, 0.13%] | 5.00% | 10 | 0.12% [0.07%, 0.23%] | 21.00% | holds |

**ETH** (8,000 streams per setting)

| Setting | A | Alarms in day 1 | Rate day 1 [Wilson 95%] | Bound | Alarms in 5 days | Rate 5 days [Wilson 95%] | Bound | Verdict |
|---|---|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 28,800 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 2 | 0.03% [0.01%, 0.09%] | 25.00% | holds |
| plug-in, band 0.1 | 28,800 | 20 | 0.25% [0.16%, 0.39%] | 5.00% | 104 | 1.30% [1.07%, 1.57%] | 25.00% | holds |
| plug-in, band 0.05 | 28,800 | 87 | 1.09% [0.88%, 1.34%] | 5.00% | 713 | 8.91% [8.31%, 9.56%] | 25.00% | holds |
| plug-in, band 0 | 28,800 | 438 | 5.47% [5.00%, 6.00%] | 5.00% | 3316 | 41.45% [40.37%, 42.53%] | 25.00% | day 1 holds, 5 days above |
| cs baseline, band 0.1 | 36,000 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 1 | 0.01% [0.00%, 0.07%] | 21.00% | holds |
| cs baseline, band 0 | 36,000 | 1 | 0.01% [0.00%, 0.07%] | 5.00% | 8 | 0.10% [0.05%, 0.20%] | 21.00% | holds |

**XRP** (8,000 streams per setting)

| Setting | A | Alarms in day 1 | Rate day 1 [Wilson 95%] | Bound | Alarms in 5 days | Rate 5 days [Wilson 95%] | Bound | Verdict |
|---|---|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 28,800 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 2 | 0.03% [0.01%, 0.09%] | 25.00% | holds |
| plug-in, band 0.1 | 28,800 | 18 | 0.22% [0.14%, 0.36%] | 5.00% | 119 | 1.49% [1.24%, 1.78%] | 25.00% | holds |
| plug-in, band 0.05 | 28,800 | 95 | 1.19% [0.97%, 1.45%] | 5.00% | 741 | 9.26% [8.65%, 9.92%] | 25.00% | holds |
| plug-in, band 0 | 28,800 | 479 | 5.99% [5.49%, 6.53%] | 5.00% | 3403 | 42.54% [41.46%, 43.62%] | 25.00% | day 1 above, 5 days above |
| cs baseline, band 0.1 | 36,000 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 0 | 0.00% [0.00%, 0.05%] | 21.00% | holds |
| cs baseline, band 0 | 36,000 | 1 | 0.01% [0.00%, 0.07%] | 5.00% | 10 | 0.12% [0.07%, 0.23%] | 21.00% | holds |

**SOL** (8,000 streams per setting)

| Setting | A | Alarms in day 1 | Rate day 1 [Wilson 95%] | Bound | Alarms in 5 days | Rate 5 days [Wilson 95%] | Bound | Verdict |
|---|---|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 28,800 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 1 | 0.01% [0.00%, 0.07%] | 25.00% | holds |
| plug-in, band 0.1 | 28,800 | 18 | 0.22% [0.14%, 0.36%] | 5.00% | 117 | 1.46% [1.22%, 1.75%] | 25.00% | holds |
| plug-in, band 0.05 | 28,800 | 89 | 1.11% [0.90%, 1.37%] | 5.00% | 712 | 8.90% [8.30%, 9.54%] | 25.00% | holds |
| plug-in, band 0 | 28,800 | 458 | 5.73% [5.24%, 6.26%] | 5.00% | 3375 | 42.19% [41.11%, 43.27%] | 25.00% | day 1 above, 5 days above |
| cs baseline, band 0.1 | 36,000 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 1 | 0.01% [0.00%, 0.07%] | 21.00% | holds |
| cs baseline, band 0 | 36,000 | 2 | 0.03% [0.01%, 0.09%] | 5.00% | 15 | 0.19% [0.11%, 0.31%] | 21.00% | holds |

**DOGE** (8,000 streams per setting)

| Setting | A | Alarms in day 1 | Rate day 1 [Wilson 95%] | Bound | Alarms in 5 days | Rate 5 days [Wilson 95%] | Bound | Verdict |
|---|---|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 28,800 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 2 | 0.03% [0.01%, 0.09%] | 25.00% | holds |
| plug-in, band 0.1 | 28,800 | 21 | 0.26% [0.17%, 0.40%] | 5.00% | 113 | 1.41% [1.18%, 1.70%] | 25.00% | holds |
| plug-in, band 0.05 | 28,800 | 87 | 1.09% [0.88%, 1.34%] | 5.00% | 644 | 8.05% [7.47%, 8.67%] | 25.00% | holds |
| plug-in, band 0 | 28,800 | 433 | 5.41% [4.94%, 5.93%] | 5.00% | 3328 | 41.60% [40.52%, 42.68%] | 25.00% | day 1 holds, 5 days above |
| cs baseline, band 0.1 | 36,000 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 0 | 0.00% [0.00%, 0.05%] | 21.00% | holds |
| cs baseline, band 0 | 36,000 | 2 | 0.03% [0.01%, 0.09%] | 5.00% | 8 | 0.10% [0.05%, 0.20%] | 21.00% | holds |

**Nifty 50** (8,000 streams per setting)

| Setting | A | Alarms in day 1 | Rate day 1 [Wilson 95%] | Bound | Alarms in 5 days | Rate 5 days [Wilson 95%] | Bound | Verdict |
|---|---|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 7,500 | 1 | 0.01% [0.00%, 0.07%] | 5.00% | 5 | 0.06% [0.03%, 0.15%] | 25.00% | holds |
| plug-in, band 0.1 | 7,500 | 27 | 0.34% [0.23%, 0.49%] | 5.00% | 225 | 2.81% [2.47%, 3.20%] | 25.00% | holds |
| plug-in, band 0.05 | 7,500 | 86 | 1.07% [0.87%, 1.33%] | 5.00% | 853 | 10.66% [10.00%, 11.36%] | 25.00% | holds |
| plug-in, band 0 | 7,500 | 230 | 2.88% [2.53%, 3.26%] | 5.00% | 2508 | 31.35% [30.34%, 32.38%] | 25.00% | day 1 holds, 5 days above |
| cs baseline, band 0.1 | 9,375 | 0 | 0.00% [0.00%, 0.05%] | 5.00% | 0 | 0.00% [0.00%, 0.05%] | 21.00% | holds |
| cs baseline, band 0 | 9,375 | 1 | 0.01% [0.00%, 0.07%] | 5.00% | 4 | 0.05% [0.02%, 0.13%] | 21.00% | holds |

## Delay tables [SYNTHETIC]

How to read this: each row is one setting on simulated paths with 1,000 calibration steps and a variance change at monitored step 400. 'Detected' is the share of paths whose first alarm came after the change and within 1,950 steps (five days of 390). 'Early' is the share whose first alarm came before the change (a false alarm). The median delay is the middle number of steps from the change to the alarm, among detected paths; the IQR is the middle half.

**gaussian, variance x2** (2,000 paths per setting)

| Setting | A | Detected [Wilson 95%] | Early | Median delay | IQR | Mean delay |
|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 7,800 | 100.00% [99.81%, 100.00%] | 0.00% | 114 | 77 to 163 | 128 |
| plug-in, band 0.1 | 7,800 | 99.70% [99.35%, 99.86%] | 0.30% | 70 | 48 to 98 | 78 |
| plug-in, band 0.05 | 7,800 | 99.10% [98.58%, 99.43%] | 0.90% | 59 | 41 to 85 | 66 |
| plug-in, band 0 | 7,800 | 97.85% [97.12%, 98.40%] | 2.15% | 50 | 35 to 71 | 56 |
| cs baseline, band 0 | 7,800 | 100.00% [99.81%, 100.00%] | 0.00% | 114 | 75 to 163 | 128 |
| cs baseline, band 0.1 | 7,800 | 100.00% [99.81%, 100.00%] | 0.00% | 177 | 116 to 261 | 202 |
| cs baseline, band 0 | 9,750 | 100.00% [99.81%, 100.00%] | 0.00% | 118 | 80 to 168 | 133 |
| cs baseline, band 0.1 | 9,750 | 100.00% [99.81%, 100.00%] | 0.00% | 189 | 120 to 267 | 210 |

**gaussian, variance x0.5** (2,000 paths per setting)

| Setting | A | Detected [Wilson 95%] | Early | Median delay | IQR | Mean delay |
|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 7,800 | 100.00% [99.81%, 100.00%] | 0.00% | 262 | 225 to 311 | 275 |
| plug-in, band 0.1 | 7,800 | 99.90% [99.64%, 99.97%] | 0.10% | 154 | 133 to 178 | 157 |
| plug-in, band 0.05 | 7,800 | 99.30% [98.83%, 99.58%] | 0.70% | 132 | 112 to 153 | 132 |
| plug-in, band 0 | 7,800 | 97.85% [97.12%, 98.40%] | 2.15% | 113 | 93 to 132 | 112 |
| cs baseline, band 0 | 7,800 | 100.00% [99.81%, 100.00%] | 0.00% | 195 | 163 to 235 | 202 |
| cs baseline, band 0.1 | 7,800 | 100.00% [99.81%, 100.00%] | 0.00% | 270 | 223 to 330 | 284 |
| cs baseline, band 0 | 9,750 | 100.00% [99.81%, 100.00%] | 0.00% | 203 | 169 to 243 | 210 |
| cs baseline, band 0.1 | 9,750 | 100.00% [99.81%, 100.00%] | 0.00% | 280 | 231 to 340 | 295 |

**student_t3, variance x2** (2,000 paths per setting)

| Setting | A | Detected [Wilson 95%] | Early | Median delay | IQR | Mean delay |
|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 7,800 | 99.95% [99.72%, 99.99%] | 0.00% | 281 | 173 to 433 | 334 |
| plug-in, band 0.1 | 7,800 | 99.60% [99.21%, 99.80%] | 0.40% | 140 | 87 to 209 | 160 |
| plug-in, band 0.05 | 7,800 | 98.75% [98.16%, 99.15%] | 1.25% | 112 | 72 to 167 | 128 |
| plug-in, band 0 | 7,800 | 96.65% [95.77%, 97.35%] | 3.35% | 89 | 57 to 134 | 103 |
| cs baseline, band 0 | 7,800 | 100.00% [99.81%, 100.00%] | 0.00% | 262 | 160 to 411 | 314 |
| cs baseline, band 0.1 | 7,800 | 93.85% [92.71%, 94.82%] | 0.00% | 516 | 304 to 804 | 606 |
| cs baseline, band 0 | 9,750 | 100.00% [99.81%, 100.00%] | 0.00% | 275 | 166 to 429 | 329 |
| cs baseline, band 0.1 | 9,750 | 93.20% [92.01%, 94.22%] | 0.00% | 540 | 322 to 842 | 631 |

**student_t3, variance x0.5** (2,000 paths per setting)

| Setting | A | Detected [Wilson 95%] | Early | Median delay | IQR | Mean delay |
|---|---|---|---|---|---|---|
| plug-in, band 0.25 (README) | 7,800 | 98.25% [97.58%, 98.74%] | 0.00% | 513 | 323 to 765 | 593 |
| plug-in, band 0.1 | 7,800 | 99.75% [99.42%, 99.89%] | 0.25% | 214 | 148 to 298 | 239 |
| plug-in, band 0.05 | 7,800 | 98.60% [97.98%, 99.03%] | 1.40% | 171 | 121 to 235 | 190 |
| plug-in, band 0 | 7,800 | 96.90% [96.05%, 97.57%] | 3.10% | 139 | 98 to 194 | 154 |
| cs baseline, band 0 | 7,800 | 99.70% [99.35%, 99.86%] | 0.00% | 341 | 232 to 530 | 416 |
| cs baseline, band 0.1 | 7,800 | 93.45% [92.28%, 94.45%] | 0.00% | 601 | 372 to 917 | 695 |
| cs baseline, band 0 | 9,750 | 99.70% [99.35%, 99.86%] | 0.00% | 364 | 243 to 556 | 438 |
| cs baseline, band 0.1 | 9,750 | 92.20% [90.94%, 93.30%] | 0.00% | 637 | 399 to 962 | 725 |
