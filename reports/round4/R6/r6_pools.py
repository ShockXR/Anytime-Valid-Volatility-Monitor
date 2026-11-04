"""R6 step 0: build E19's calm-day pools with E19's own code (prepare_pools), written to reports/round4/R6/pools.
E19_OUT and E19_SCRATCH must point to that folder so that results/ is not touched."""
import sys, time
sys.path.insert(0, "/home/user/anytime-valid-volatility-monitor/experiments")
import e19_real_null as E19
t0 = time.time()
df = E19.prepare_pools(list(E19.DATASETS))
print(df.to_string())
print("pools built in %.0f s" % (time.time() - t0), flush=True)
