# Data

Raw and processed data are not committed to this repository.

## Bitstamp BTC/USD 1-minute bars (used in E5)

Public file maintained at <https://github.com/ff137/bitstamp-btcusd-minute-data>. The daily-updated file covers 7 January 2025 onwards. Columns: `timestamp` (UTC open time, seconds), `open`, `high`, `low`, `close`, `volume`. Missing minutes are filled with flat zero-volume candles.

```bash
mkdir -p data/raw/crypto
curl -L -o data/raw/crypto/btcusd_bitstamp_1min_latest.csv \
  https://raw.githubusercontent.com/ff137/bitstamp-btcusd-minute-data/main/data/updates/btcusd_bitstamp_1min_latest.csv
```

The file used for the committed results was downloaded on 2 November 2026 (956,327 rows, last bar 2026-11-02). Bapat (2026) does not name the exchange behind its crypto data, so Bitstamp is a stand-in, not the same dataset.

## Expected layout

```
data/
├── raw/
│   ├── crypto/        1-minute bars (Bitstamp BTC/USD now; ETH, XRP, SOL, DOGE later)
│   ├── nifty/         1-minute Nifty 50 bars (planned)
│   └── india_vix/     daily India VIX (planned)
└── processed/         cleaned returns and detector outputs
```
