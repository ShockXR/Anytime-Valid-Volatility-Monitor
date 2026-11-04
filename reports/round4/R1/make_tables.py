"""Markdown tables for REPORT.md from integrity_files.csv and integrity_top10.csv."""
import pandas as pd, numpy as np
f = pd.read_csv("integrity_files.csv"); t = pd.read_csv("integrity_top10.csv")
out = []
def pct(x): return f"{100*x:.1f}%"
# Binance per coin (manifest months 2024-01 to 2025-12, and the E23 2026 months separately)
b = f[f.kind == "binance"].copy(); b["coin"] = b.file.str.split("/").str[4]
for grp, label in (("manifest", "Binance, committed manifest months (2024-01 to 2025-12)"), ("e23_2026", "Binance, E23 fresh months (2026-01 to 2026-10), checked against results/e23_data_files.csv")):
    g = b[b.group == grp]
    out.append(f"\n**{label}**\n")
    out.append("| Coin | Files | SHA-256 match | Rows match | Rows | Dup. ts | Missing minutes | Bad prices | OHLC inconsistent | Zero-vol bars | Flat bars | Zero-return share (min to max) | Returns > 20 robust sd | Largest return / sd | Largest abs. return | Top-10 flagged |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for coin, c in g.groupby("coin"):
        out.append(f"| {coin} | {len(c)} | {int(c.sha_match.sum())}/{len(c)} | {int(c.rows_match.sum())}/{len(c)} | {int(c.rows.sum()):,} | {int(c.duplicate_ts.sum())} | {int(c.missing_minutes.sum()):,} | {int(c[['bad_open','bad_high','bad_low','bad_close']].sum().sum())} | {int(c.ohlc_inconsistent.sum())} | {int(c.zero_volume_bars.sum())} | {int(c.flat_bars.sum())} | {pct(c.zero_return_share.min())} to {pct(c.zero_return_share.max())} | {int(c.beyond_20sd.sum())} | {c.max_return_over_sd.max():.0f} | {pct(c.max_abs_return.max())} | {int(c.top10_suspect.sum())}/{10*len(c)} |")
# Nifty per file
n = f[f.kind == "nifty"]
out.append("\n**Nifty 50 (data/raw/nifty, one file per year, then per month in 2026)**\n")
out.append("| File | SHA-256 | Rows (= manifest) | Sessions | Short sessions | Missing minutes in sessions | Bars outside 09:15 to 15:29 | Flat bars | Zero-return share | Returns > 20 sd | Largest return / sd | Largest abs. return | Top-10 flagged |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in n.itertuples():
    out.append(f"| {r.file.split('/')[-1]} | {'yes' if r.sha_match else 'NO'} | {r.rows:,} ({'yes' if r.rows_match else 'NO'}) | {int(r.sessions)} | {int(r.short_sessions)} | {int(r.missing_minutes)} | {int(r.bars_outside_0915_1529)} | {int(r.flat_bars)} | {pct(r.zero_return_share)} | {int(r.beyond_20sd)} | {r.max_return_over_sd:.0f} | {pct(r.max_abs_return)} | {int(r.top10_suspect)}/10 |")
# Bitstamp
bs = f[f.kind == "bitstamp"].iloc[0]
out.append("\n**Bitstamp BTC/USD (daily-updated file, rows up to the committed last bar 2026-11-02 02:47 UTC)**\n")
out.append("| Check | Value |\n|---|---|")
for k, v in [("Rows on disk today", f"{int(bs.rows_on_disk_today):,}"), ("Rows up to 2026-11-02 23:59 UTC", f"{int(bs.rows):,}"), ("Committed manifest rows", f"{int(bs.rows_manifest):,}"),
             ("Last bar of the first 956,327 rows", "2026-11-02 02:47 UTC"), ("SHA-256 of header + first 956,327 rows", "810a84b5... = committed manifest (byte-exact prefix, 53,449,859 bytes)"),
             ("First bar", bs.first_bar), ("Duplicate timestamps", int(bs.duplicate_ts)), ("Missing minutes (gaps)", f"{int(bs.missing_minutes)} ({int(bs.gaps)} gaps)"), ("Strictly increasing, whole minutes", f"{bs.strictly_increasing}, {bs.whole_minutes}"),
             ("Bad prices, OHLC inconsistent", f"{int(bs[['bad_open','bad_high','bad_low','bad_close']].sum())}, {int(bs.ohlc_inconsistent)}"), ("Zero-volume bars (filled minutes)", f"{int(bs.zero_volume_bars):,}"), ("Flat bars", f"{int(bs.flat_bars):,}"),
             ("Zero-return share", pct(bs.zero_return_share)), ("Returns > 20 robust sd", int(bs.beyond_20sd)), ("Largest return / sd, largest abs. return", f"{bs.max_return_over_sd:.0f}, {pct(bs.max_abs_return)}"), ("Top-10 flagged", f"{int(bs.top10_suspect)}/10")]:
    out.append(f"| {k} | {v} |")
# top10 verdict summary
t["src"] = np.where(t.file.str.contains("NIFTY"), "Nifty 50", np.where(t.file.str.contains("bitstamp"), "Bitstamp", t.file.str.split("_").str[0]))
v = t.groupby("src").agg(n=("verdict", "size"), real=("verdict", lambda s: s.str.startswith("real").sum()), rev=("verdict", lambda s: s.str.contains("reversed").sum()), big=("verdict", lambda s: s.str.contains("larger").sum()), max_sd=("ret_over_sd", "max"))
out.append("\n**The 10 largest absolute one-minute returns of every file (1,700 returns in all)**\n")
out.append("| Source | Returns looked at | Range covers the move and no instant reversal | At least half reversed in the next minute | Move larger than the two bars' high-low span | Largest return / robust sd |\n|---|---|---|---|---|---|")
for s, r in v.iterrows():
    out.append(f"| {s} | {int(r.n)} | {int(r.real)} | {int(r.rev)} | {int(r.big)} | {r.max_sd:.0f} |")
open("tables.md", "w").write("\n".join(out) + "\n"); print("\n".join(out))
