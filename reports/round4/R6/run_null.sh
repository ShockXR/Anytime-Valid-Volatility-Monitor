#!/bin/bash
# R6: 8 batches of 1000 permuted streams per data set (seeds 9_970_000 + 100 d + b)
for ds in BTC ETH XRP SOL DOGE "Nifty 50"; do
  for b in 0 1 2 3 4 5 6 7; do
    python r6_null.py "$ds" $b 1000
  done
done
echo EXIT
