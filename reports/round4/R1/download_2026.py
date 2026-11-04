"""Download the 2026 Binance months that E23 uses (BTC, ETH, SOL, 2025-12 to 2026-10), with the package's own function."""
import sys
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/experiments")
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor")
import e23_fresh_releases as E23
import e12_crypto_events as E12
for coin, sym in E12.COINS.items():
    ps = E23.download_months(sym)
    print(coin, len(ps), ps[-1], flush=True)
print("DONE2026")
