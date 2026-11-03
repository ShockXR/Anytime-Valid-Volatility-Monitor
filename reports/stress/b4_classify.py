"""Assign the failure class to each raw case of b4_raw_cases.csv (reviewed by hand against the observed column).
DETECTED = the package raises its own error naming the problem. SAFELY HANDLED = continues with a sensible, documented result.
SILENTLY IGNORED = accepted with no error, wrong or degenerate result and no warning. CRASHED = loud failure from library internals
(no package message) or a delayed failure far from the bad input. Run: python b4_classify.py"""
import pandas as pd
d = pd.read_csv("b4_raw_cases.csv", keep_default_na=False)
D_, S_, I_, C_ = "DETECTED", "SAFELY HANDLED", "SILENTLY IGNORED", "CRASHED"
K = {  # row index -> (class, root cause, note)
 0:(S_,"by design","NaN skipped as missing, state untouched"), 1:(S_,"by design","300 NaN skipped"), 2:(S_,"by design","NaN dropped from calibration window"),
 3:(I_,"implementation gap","no alarm and no warning, calibration restarts forever; only cal_stats['restarted'] shows it"),
 4:(S_,"by design","non-finite returns treated as missing"), 5:(S_,"by design","non-finite dropped from calibration"),
 6:(I_,"implementation gap","constant price: monitor blind forever, no warning, only cal_stats['restarted']=5"),
 7:(S_,"by design",""), 8:(S_,"by design","documented fallback scale path, counted in cal_stats"), 9:(S_,"by design","no-op"), 10:(S_,"by design","no-op"),
 11:(D_,"by design",""), 12:(D_,"by design",""), 13:(D_,"by design",""),
 14:(S_,"by design","truncation bounds the influence of one huge return"), 15:(S_,"by design",""),
 16:(I_,"invalid assumption","10 consecutive 200-sd returns inside calibration (null_band=0) inflate the baseline and cause a later false down alarm; see sensitivity table"),
 17:(S_,"by design","RuntimeWarning only"),
 18:(C_,"implementation gap","accepted silently during calibration, TypeError from numpy at the end of calibration (delayed)"),
 19:(C_,"implementation gap","same as the string: accepted silently during calibration, TypeError when calibration ends, or at once while monitoring"),
 20:(C_,"implementation gap","2-D input iterates rows, accepted during calibration, fails later"),
 21:(I_,"implementation gap","cal_len=1 accepted: many false alarms (12 in 3000 null steps, expected 0.4); no minimum cal_len check"),
 22:(I_,"implementation gap","cal_len=2 accepted, same"), 23:(I_,"implementation gap","cal_len=5 accepted, plug-in baseline very noisy"),
 24:(I_,"implementation gap","cal_len=0 behaves as 1"), 25:(D_,"accidental","message is about cal_max, not cal_len"),
 26:(I_,"implementation gap","cal_len=10 accepted: 8 alarms in 3000 null steps"),
 27:(I_,"implementation gap","threshold<=0 accepted: alarm at every first monitored step"), 28:(I_,"implementation gap","threshold<0 accepted"),
 29:(I_,"implementation gap","threshold=1 accepted although the bound m/A exceeds 1"),
 30:(I_,"implementation gap","threshold=NaN accepted: never alarms"), 31:(I_,"by design","threshold=inf never alarms"),
 32:(D_,"accidental","raised by c_down check, message names c_down not c"), 33:(D_,"accidental","same"),
 34:(I_,"implementation gap","c=NaN accepted: calibration never completes, monitor blind forever"),
 35:(S_,"by design","c=1e6 is a valid but non-robust setting (no truncation)"), 36:(S_,"by design",""),
 37:(D_,"by design",""), 38:(I_,"implementation gap","null_band>=1 silently switches the decrease side off (down_active False)"),
 39:(D_,"by design",""),40:(D_,"by design",""),41:(D_,"by design",""),42:(D_,"by design",""), 43:(S_,"by design","delta ignored for plug-in"),
 44:(D_,"by design",""),45:(D_,"by design",""),46:(D_,"by design",""),47:(D_,"by design",""),48:(D_,"by design",""),49:(D_,"by design",""),
 50:(D_,"by design",""),51:(D_,"by design",""),52:(S_,"by design",""),53:(D_,"by design",""),54:(S_,"by design","floor above data scale means no alarms below the floor"),
 55:(D_,"by design",""), 56:(S_,"by design",""),
 57:(S_,"by design","pickle and deepcopy resume gives identical alarms; there is no save/load API; a cold restart re-calibrates and cannot see a shift that predates its calibration"),
 58:(S_,"by design","NaN gap leaves later alarms unchanged apart from the index shift"),
 59:(S_,"control",""),
 60:(I_,"implementation gap","zero close gives +-inf returns, RuntimeWarning only (Binance loader masks non-positive prices, this one does not)"),
 61:(I_,"implementation gap","negative close becomes NaN returns, RuntimeWarning only"),
 62:(S_,"by design",""), 63:(S_,"by design","sorted silently"),
 64:(I_,"data problem","duplicates dropped silently, number not reported; see rows 89 and 90 for which duplicate is kept"),
 65:(I_,"implementation gap","gap of 3 minutes gives a 4-minute return treated as a 1-minute return; the source file normally fills gaps"),
 66:(S_,"by design",""),67:(C_,"environment/library","pandas EmptyDataError"),68:(C_,"implementation gap","obscure TypeError on a header-only file"),
 69:(C_,"library","KeyError 'close'"),70:(C_,"library","TypeError from numpy"),
 71:(I_,"implementation gap","millisecond timestamps read as seconds give dates around year 55000, no check"),
 72:(S_,"by design","truncated row becomes NaN"),73:(C_,"library","space before a quote makes the price a string: TypeError"),
 74:(S_,"by design",""),75:(S_,"by design","flat bars give zero returns, handled by the calibration extension"),
 76:(S_,"control",""),77:(I_,"implementation gap","zero close gives +-inf returns, RuntimeWarning only"),78:(I_,"implementation gap","negative close gives NaN, RuntimeWarning only"),
 79:(D_,"by design","reported in attrs['checks']['duplicates']"),80:(S_,"by design","next return is NaN"),81:(S_,"by design",""),
 82:(C_,"library","KeyError from pandas"),83:(C_,"implementation gap","IndexError on an empty file body"),
 84:(S_,"by design","documented NaN start; a monitor fed all-NaN scale is blind (see row 3)"),85:(S_,"by design",""),
 86:(C_,"implementation gap","delayed TypeError"),87:(C_,"implementation gap","delayed ValueError"),
 88:(S_,"by design","0 of 300 paths alarm after a reverting 5 percent bar"),
 89:(I_,"implementation gap","sort_index() is not stable, so which duplicate is kept is arbitrary (Nifty loader uses a stable sort)"),
 90:(I_,"implementation gap","original-first file: the later absurd duplicate was kept for about half of the minutes"),
 91:(S_,"by design","intraday_scale does not depend on row order"),
 92:(I_,"invalid assumption","25 outliers of 15 sd in a 500-step calibration window give 69 percent false alarms within 2500 steps; 10 outliers 0.7 percent"),
}
assert set(K) == set(range(len(d))), set(range(len(d))) ^ set(K)
d.insert(0, "id", [f"B4-{i+1:02d}" for i in range(len(d))])
d["class"] = [K[i][0] for i in range(len(d))]; d["root_cause"] = [K[i][1] for i in range(len(d))]; d["note"] = [K[i][2] for i in range(len(d))]
d.to_csv("b4_failure_cases.csv", index=False)
print(d["class"].value_counts().to_string())
