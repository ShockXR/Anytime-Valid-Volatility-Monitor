"""R6 step 3: summarise the null batches (r6_null.py) into r6_null_rows.csv, with Wilson 95% intervals and
E19's bounds (E19.bounds). One row per data set, monitor and band, plus 'five coins' pooled rows.
Usage: python r6_summarise.py"""
import os, sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
os.environ.setdefault("E19_OUT", str(HERE / "pools"))
os.environ.setdefault("E19_SCRATCH", str(HERE / "pools"))
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/experiments")
import e19_real_null as E19  # noqa: E402

parts = []
for f in sorted((HERE / "null_batches").glob("*.pkl")):
    parts.extend(pd.read_pickle(f))
cat = {}
for r in parts:
    cat.setdefault((r["dataset"], r["monitor"], r["band"]), []).append(r)


def row(label, ps, name_for_bounds):
    day = E19.DATASETS[name_for_bounds][0]
    first = np.concatenate([p["first_t"] for p in ps])
    fdir = np.concatenate([p["first_dir"] for p in ps])
    ended = np.concatenate([p["ended"] for p in ps])
    n = first.size
    k1 = int(np.sum((first > 0) & (first <= day)))
    k5 = int(np.sum((first > 0) & (first <= 5 * day)))
    b1, b5, _ = E19.bounds(name_for_bounds, label[1])
    p1, l1, h1 = E19.wilson(k1, n)
    p5, l5, h5 = E19.wilson(k5, n)
    return dict(dataset=label[0], monitor=label[1], band=label[2], A=ps[0]["A"], streams=n, batches=len(ps),
                seeds=",".join(str(p["seed"]) for p in ps),
                one_day=p1, one_day_lo=l1, one_day_hi=h1, one_day_bound=b1, one_day_alarms=k1,
                five_days=p5, five_days_lo=l5, five_days_hi=h5, five_days_bound=b5, five_days_alarms=k5,
                first_alarm_up_share=float(np.mean(fdir[first > 0] == 1)) if np.any(first > 0) else np.nan,
                streams_ended_early=int(ended.sum()))


rows = []
order = lambda k: (list(E19.DATASETS).index(k[0]), ["plugin", "cs"].index(k[1]), -k[2])
for key in sorted(cat, key=order):
    rows.append(row(key, cat[key], key[0]))
for mon in ("plugin", "cs"):
    for band in sorted({k[2] for k in cat if k[1] == mon}, reverse=True):
        ps = [p for c in E19.CRYPTO for p in cat.get((c, mon, band), [])]
        if ps:
            rows.append(row(("five coins", mon, band), ps, "BTC"))
df = pd.DataFrame(rows)
for col in ("one_day", "five_days"):
    df[f"{col}_vs_bound"] = np.where(df[f"{col}_lo"] > df[f"{col}_bound"], "above", "at or below")
df.to_csv(HERE / "r6_null_rows.csv", index=False)
show = df[["dataset", "monitor", "band", "A", "streams", "one_day", "one_day_lo", "one_day_hi", "one_day_bound",
           "five_days", "five_days_lo", "five_days_hi", "five_days_bound", "one_day_vs_bound", "five_days_vs_bound"]]
pd.set_option("display.width", 250)
print(show.to_string(float_format=lambda x: f"{x:.4f}"))
