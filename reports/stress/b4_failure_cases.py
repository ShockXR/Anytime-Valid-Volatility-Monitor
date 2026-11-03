"""B4. Failure simulation for the monitor and the loaders (SYNTHETIC inputs). Seeds: 9_500_000 + case number.
Writes one row per case with the observed behaviour; the class column is assigned in b4_classify.py after review."""
import copy, io, os, pickle, sys, tempfile, warnings
import numpy as np, pandas as pd
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/reports/stress")
from avmonitor.monitor import VolatilityMonitor, GaussianLRMonitor
from avmonitor import detectors as D, data as DA

ROWS = []
def case(name, fn):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            obs = fn()
            exc = ""
        except BaseException as e:  # noqa
            obs, exc = "", f"{type(e).__name__}: {str(e)[:110]}"
        wn = sorted({f"{x.category.__name__}: {str(x.message)[:60]}" for x in w})
    ROWS.append(dict(case=name, observed=str(obs)[:230], exception=exc, warnings=" | ".join(wn)[:200]))

def base(n=3000, seed=9_500_000):
    return np.random.default_rng(seed).standard_normal(n)

def M(**kw):
    kw.setdefault("threshold", 7800.0); kw.setdefault("cal_len", 500)
    return VolatilityMonitor(**kw)

def summ(m):
    return f"events={[(e.index, e.direction) for e in m.events][:4]} n_events={len(m.events)} calibrated={m.calibrated} cal_stats={m.cal_stats} t={m.t}"

# ---- data-value failures fed to the streaming monitor
def run_with(r, **kw):
    m = M(**kw); m.run(r); return summ(m)
def nan_mon():
    r = base(); r[1000] = np.nan; return run_with(r)
case("NaN single point in monitoring", nan_mon)
def nan_block():
    r = base(); r[800:1100] = np.nan; return run_with(r)
case("NaN block of 300 in monitoring", nan_block)
def nan_cal():
    r = base(); r[10:60] = np.nan; return run_with(r)
case("NaN 50 points inside calibration window", nan_cal)
def all_nan():
    return run_with(np.full(3000, np.nan))
case("all NaN stream", all_nan)
def inf_mon():
    r = base(); r[1000] = np.inf; r[1001] = -np.inf; return run_with(r)
case("+inf and -inf in monitoring", inf_mon)
def inf_cal():
    r = base(); r[100] = np.inf; return run_with(r)
case("inf inside calibration window", inf_cal)
def const_zero():
    return run_with(np.zeros(5000))
case("all-zero returns (constant price) 5000 steps", const_zero)
def const_nonzero():
    return run_with(np.full(3000, 0.001))
case("constant nonzero return 0.001", const_nonzero)
def mostly_zero():
    r = base(); r[np.random.default_rng(1).random(3000) < 0.6] = 0.0; return run_with(r)
case("60 percent zero returns (frozen/tick data)", mostly_zero)
def empty():
    m = M(); m.run(np.array([])); return summ(m)
case("empty input to run()", empty)
def one_pt():
    m = M(); m.run([0.01]); return summ(m)
case("one-point input to run()", one_pt)
case("calibrate() on empty window", lambda: D.calibrate(np.array([])))
case("calibrate() on one point", lambda: D.calibrate(np.array([0.01])))
case("calibrate() on all-zero window", lambda: D.calibrate(np.zeros(100)))
def spike_mon():
    r = base(); r[1500] = 1e6; return run_with(r)
case("fat finger: one return of 1e6 sd in monitoring", spike_mon)
def spike_pair():
    r = base(); r[1500] = 0.5; r[1501] = -0.5 ; r = r * 1.0; return run_with(r)  # +50 sd then -50 sd (bounce back)
case("fat finger: +50 sd then -50 sd bounce", spike_pair)
def spike_cal():
    r = base(); r[100:110] = 200.0; return run_with(r)
case("10 consecutive 200 sd returns inside calibration", spike_cal)
def huge():
    r = base(); r[1500] = 1e300; return run_with(r)
case("return 1e300 (square overflows)", huge)
case("string input to update()", lambda: M().update("abc"))
case("None input to update()", lambda: M().update(None))
case("2-D array to run()", lambda: M().run(np.ones((3, 3))))
# ---- parameters
for name, kw in [("cal_len=1", dict(cal_len=1)), ("cal_len=2", dict(cal_len=2)), ("cal_len=5", dict(cal_len=5)),
                 ("cal_len=0", dict(cal_len=0)), ("cal_len=-5", dict(cal_len=-5)), ("cal_len=10", dict(cal_len=10)),
                 ("threshold=0", dict(threshold=0.0)), ("threshold=-5", dict(threshold=-5.0)),
                 ("threshold=1 (A<2)", dict(threshold=1.0)), ("threshold=NaN", dict(threshold=float("nan"))),
                 ("threshold=inf", dict(threshold=float("inf"))),
                 ("c=0", dict(c=0.0)), ("c=-3", dict(c=-3.0)), ("c=NaN", dict(c=float("nan"))), ("c=1e6", dict(c=1e6)),
                 ("c=0.1", dict(c=0.1)),
                 ("null_band=-0.1", dict(null_band=-0.1)), ("null_band=2 (decrease side dead)", dict(null_band=2.0)),
                 ("cs delta=0", dict(baseline="cs", delta=0.0)), ("cs delta=1", dict(baseline="cs", delta=1.0)),
                 ("cs delta=-0.1", dict(baseline="cs", delta=-0.1)), ("cs delta=NaN", dict(baseline="cs", delta=float("nan"))),
                 ("plugin delta=0 (ignored)", dict(delta=0.0)),
                 ("cs tau_frac=0", dict(baseline="cs", tau_frac=0.0)), ("cs tau_frac=1", dict(baseline="cs", tau_frac=1.0)),
                 ("c_down=0", dict(c_down=0.0)), ("c_down=-1", dict(c_down=-1.0)),
                 ("up_grid containing 1.0", dict(up_grid=(0.5, 1.0))), ("up_grid empty", dict(up_grid=())),
                 ("up_grid negative", dict(up_grid=(-0.1,))), ("down_fractions containing 1.0", dict(down_fractions=(1.0,))),
                 ("margin_z=2 plugin", dict(margin_z=2.0)), ("cal_max<cal_len", dict(cal_max=10)),
                 ("scale_floor=1e6", dict(scale_floor=1e6)), ("baseline='foo'", dict(baseline="foo"))]:
    def f(kw=kw):
        kw = dict(kw); kw.setdefault("cal_len", 500)
        m = VolatilityMonitor(**{"threshold": 7800.0, **kw}); m.run(base(3000)); return summ(m)
    case("param " + name, f)
def cs_len_small():
    m = VolatilityMonitor(threshold=7800.0, cal_len=3, baseline="cs"); m.run(base(300)); return summ(m)
case("cs baseline with cal_len=3", cs_len_small)
# ---- interrupt and restart
def interrupt():
    r = np.concatenate([base(3000, 9_500_100), 2.0 * base(3000, 9_500_101)])
    full = M(); full.run(r)
    a = M(); a.run(r[:3700]); blob = pickle.dumps(a); b = pickle.loads(blob); b.run(r[3700:])
    c = M(); c.run(r[:3700]); d = copy.deepcopy(c); d.run(r[3700:])
    fresh = M(); fresh.run(r[3700:])
    ev = lambda m: [(e.index, e.direction) for e in m.events]
    return (f"uninterrupted={ev(full)} pickle_resume_identical={ev(full)==ev(b)} deepcopy_resume_identical={ev(full)==ev(d)} "
            f"cold_restart_on_remaining_data(indices relative to resume)={[(e.index+3700,e.direction) for e in fresh.events]} "
            f"lost_calibration_steps={fresh.cal_stats}")
case("interrupt at 3700 then resume (pickle, deepcopy, cold restart)", interrupt)
def restart_gap():
    r = base(4000, 9_500_102); r2 = np.concatenate([r[:2000], np.full(200, np.nan), r[2000:]]); a = M(); a.run(r); b = M(); b.run(r2)
    return f"identical events with 200 NaN gap inserted: {[ (e.index,e.direction) for e in a.events]==[(e.index-200 if e.index>2000 else e.index,e.direction) for e in b.events]}; {summ(b)}"
case("200 missing returns (NaN) inserted mid-stream", restart_gap)
# ---- loaders
def tmpcsv(text):
    f = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False); f.write(text); f.close(); return f.name
def bitstamp_csv(rows, header="timestamp,open,high,low,close,volume\n"):
    return tmpcsv(header + "\n".join(rows) + "\n")
t0 = 1_700_000_000
good = [f"{t0+60*i},100,101,99,{100+0.01*i},5" for i in range(10)]
def load_ok(rows, header="timestamp,open,high,low,close,volume\n", **kw):
    df = DA.load_bitstamp_csv(bitstamp_csv(rows, header), **kw)
    return f"rows={len(df)} index_monotonic={df.index.is_monotonic_increasing} n_logret_nan={int(df.logret.isna().sum())} min_logret={df.logret.min():.4g} max_logret={df.logret.max():.4g} n_inf={int(np.isinf(df.logret).sum())}"
case("bitstamp loader: good file (control)", lambda: load_ok(good))
case("bitstamp loader: zero close price", lambda: load_ok(good[:5] + [f"{t0+60*5},100,101,99,0,5"] + good[6:]))
case("bitstamp loader: negative close price", lambda: load_ok(good[:5] + [f"{t0+60*5},100,101,99,-5,5"] + good[6:]))
case("bitstamp loader: NaN close", lambda: load_ok(good[:5] + [f"{t0+60*5},100,101,99,,5"] + good[6:]))
case("bitstamp loader: unsorted rows", lambda: load_ok(good[::-1]))
case("bitstamp loader: duplicate timestamps with different closes", lambda: load_ok(good + [f"{t0+60*3},100,101,99,150,5"]))
case("bitstamp loader: missing minutes (3 rows removed, no zero-volume fill)", lambda: load_ok(good[:4] + good[7:]))
case("bitstamp loader: zero-volume filled minute", lambda: load_ok(good[:4] + [f"{t0+60*4},100,100,100,{100+0.01*3},0"] + good[5:]))
case("bitstamp loader: empty file (no header)", lambda: DA.load_bitstamp_csv(tmpcsv("")))
case("bitstamp loader: header only", lambda: load_ok([]))
case("bitstamp loader: missing close column", lambda: load_ok([f"{t0+60*i},100,101,99,5" for i in range(5)], header="timestamp,open,high,low,volume\n"))
case("bitstamp loader: text in price column", lambda: load_ok(good[:5] + [f"{t0+60*5},100,101,99,abc,5"] + good[6:]))
case("bitstamp loader: timestamps in ms instead of s", lambda: load_ok([f"{(t0+60*i)*1000},100,101,99,{100+0.01*i},5" for i in range(10)]))
case("bitstamp loader: truncated last row", lambda: load_ok(good + [f"{t0+60*10},100,101"]))
case("bitstamp loader: extra whitespace/quotes", lambda: load_ok([f'{t0+60*i}, 100, 101, 99, "{100+0.01*i}", 5' for i in range(10)]))
case("bitstamp loader: UTF-8 BOM header", lambda: load_ok(good, header="﻿timestamp,open,high,low,close,volume\n"))
case("bitstamp loader: all identical closes (flat)", lambda: load_ok([f"{t0+60*i},100,100,100,100,5" for i in range(10)]))
def nifty_root(rows, header="Epoch,Timestamp,Open,High,Low,Close\n", year="2024"):
    root = tempfile.mkdtemp(); p = DA.nifty_file_path(DA._nifty_keys(f"{year}-01-01", f"{int(year)+1}-01-01")[0], root)
    os.makedirs(os.path.dirname(p), exist_ok=True) if os.path.dirname(p) else None
    open(p, "w").write(header + "\n".join(rows) + "\n"); return root
ist0 = 1704421500  # 2024-01-05 09:15 IST
nrows = [f"{ist0+60*i},x,100,101,99,{100+0.01*i}" for i in range(20)]
def nif(rows, **kw):
    root = nifty_root(rows, **kw)
    df = DA.load_nifty("2024-01-01", "2025-01-01", root=root, download=False)
    return f"rows={len(df)} logret_nan={int(df.logret.isna().sum())} n_inf={int(np.isinf(df.logret).sum())} checks_dups={df.attrs['checks']['duplicates']}"
case("nifty loader: good file (control)", lambda: nif(nrows))
case("nifty loader: zero close", lambda: nif(nrows[:8] + [f"{ist0+60*8},x,100,101,99,0"] + nrows[9:]))
case("nifty loader: negative close", lambda: nif(nrows[:8] + [f"{ist0+60*8},x,100,101,99,-3"] + nrows[9:]))
case("nifty loader: duplicate minute", lambda: nif(nrows + [nrows[5]]))
case("nifty loader: missing minute inside session", lambda: nif(nrows[:8] + nrows[9:]))
case("nifty loader: unsorted", lambda: nif(nrows[::-1]))
case("nifty loader: missing Close column", lambda: nif([f"{ist0+60*i},x,100,101,99" for i in range(5)], header="Epoch,Timestamp,Open,High,Low\n"))
case("nifty loader: empty body", lambda: nif([]))
# ---- deseasonalisation helpers
def scale_short():
    s = pd.Series(base(500), index=pd.date_range("2024-01-01 09:15", periods=500, freq="1min", tz="Asia/Kolkata"))
    out = DA.intraday_scale(s); return f"len={len(out)} nan={int(np.isnan(out).sum())}"
case("intraday_scale on a series shorter than one lookback", scale_short)
def mark_frozen_nan():
    s = pd.Series(np.concatenate([np.zeros(100), base(100)]), index=pd.date_range("2024-01-01", periods=200, freq="1min", tz="UTC"))
    out = DA.mark_frozen(s); return f"flagged={int(out.isna().sum()) if hasattr(out,'isna') else out.sum()}"
case("mark_frozen on 100 zero returns", mark_frozen_nan)

# ---- extra cases added after the first review of the raw output
def str_delayed():
    m = M(cal_len=50); m.update("abc")
    for x in base(60): m.update(x)
    return summ(m)
case("string fed to update() during calibration, then 60 normal returns", str_delayed)
def twod_delayed():
    m = M(cal_len=50); m.run(np.ones((3, 3)))
    for x in base(60): m.update(x)
    return summ(m)
case("2-D array fed to run(), then 60 normal returns", twod_delayed)
def tiny_spike():
    rng = np.random.default_rng(9_500_300); fails = 0; n = 300
    for i in range(n):
        r = 0.0005 * rng.standard_normal(4000); r[2500] = 0.05; r[2501] = -0.05  # a 5 percent fat-finger bar that reverts
        m = M(cal_len=500, null_band=0.25); m.run(r); fails += len(m.events) > 0
    return f"alarm in {fails}/{n} paths (spike +-100 sd at step 2500, iid N(0,0.0005^2) otherwise)"
case("fat finger +5 percent then -5 percent on 0.05 percent minute returns, README setting, 300 paths", tiny_spike)
def dup_stability():
    n = 2000; rows = [f"{t0+60*i},100,101,99,{100+0.01*i},5" for i in range(n)]
    rows += [f"{t0+60*i},100,101,99,{900+i},5" for i in range(0, n, 7)]  # later duplicate rows with absurd closes
    rng = np.random.default_rng(9_500_301); order = rng.permutation(len(rows)); shuffled = [rows[i] for i in order]
    df = DA.load_bitstamp_csv(bitstamp_csv(shuffled)); kept_bad = int((df.close > 500).sum())
    return f"rows={len(df)} rows_where_the_absurd_duplicate_was_kept={kept_bad} of {len(range(0, n, 7))} duplicated minutes (duplicates appear after the original in file order only before shuffling)"
case("bitstamp loader: 286 duplicated minutes with absurd closes, file shuffled", dup_stability)
def dup_stability_sorted():
    n = 2000; rows = [f"{t0+60*i},100,101,99,{100+0.01*i},5" for i in range(n)]
    rows += [f"{t0+60*i},100,101,99,{900+i},5" for i in range(0, n, 7)]
    df = DA.load_bitstamp_csv(bitstamp_csv(rows)); kept_bad = int((df.close > 500).sum())
    return f"rows={len(df)} absurd_duplicates_kept={kept_bad} (file order: original first, duplicate appended at the end)"
case("bitstamp loader: 286 duplicated minutes appended at end of sorted file", dup_stability_sorted)
def scale_unsorted():
    idx = pd.date_range("2024-01-01", periods=60 * 24 * 30, freq="1min", tz="UTC")
    s = pd.Series(np.random.default_rng(9_500_302).standard_normal(len(idx)), index=idx)
    a = DA.intraday_scale(s); b = DA.intraday_scale(s.iloc[::-1]).iloc[::-1]
    return f"intraday_scale(sorted) vs intraday_scale(reversed order) max abs diff = {np.nanmax(np.abs(a.values-b.values)) if np.isfinite(a.values).any() else 'all NaN'}; n_finite_sorted={int(np.isfinite(a.values).sum())} n_finite_reversed={int(np.isfinite(b.values).sum())}"
case("intraday_scale on reverse-ordered index (30 days)", scale_unsorted)
def cal_contam():
    n = 300; res = {}
    for k in (0, 10, 25):
        al = 0
        for i in range(n):
            rng = np.random.default_rng(9_500_400 + i); r = rng.standard_normal(3000)
            if k: r[rng.choice(500, k, replace=False)] = 15.0 * rng.choice([-1, 1], k)
            m = M(cal_len=500, null_band=0.25); m.run(r); al += len(m.events) > 0
        res[k] = al / n
    return f"alarm share within 2500 monitored steps by number of 15-sd outliers inside the 500-step calibration window: {res} (n=300 each)"
case("outliers inside calibration window, README setting with cal_len=500", cal_contam)

pd.DataFrame(ROWS).to_csv(sys.argv[1] if len(sys.argv) > 1 else "b4_raw_cases.csv", index=False)
print(len(ROWS), "cases written")
