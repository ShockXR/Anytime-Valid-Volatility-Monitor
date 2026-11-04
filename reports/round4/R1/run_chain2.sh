#!/bin/bash
# waits for E18 until 02:31 UTC (kills it then), runs E25 gates and E19 pools, one core each
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=/home/user/wt-repro/src
source /home/user/.venv-bapat/bin/activate
cd /home/user/wt-repro
LOGS=/home/user/anytime-valid-volatility-monitor/reports/round4/R1/logs
while kill -0 2987 2>/dev/null; do
  if [ $(date -u +%H%M) -ge 0231 ]; then echo "E18 killed at $(date -u) (over the 25 minute budget)" >> $LOGS/chain.log; pkill -f e18_landmark_events.py; sleep 3; break; fi
  sleep 10
done
echo "E18 ended at $(date -u)" >> $LOGS/chain.log
run() { name=$1; shift; echo "START $name $(date -u)" >> $LOGS/chain.log; t0=$(date +%s); python "$@" > $LOGS/$name.log 2>&1; rc=$?; echo "END $name rc=$rc wall=$(( $(date +%s) - t0 ))s $(date -u)" >> $LOGS/chain.log; }
run e25_gates experiments/e25_sequential_events.py gates --procs=1
run e19_pools experiments/e19_real_null.py pools
echo CHAIN2_DONE >> $LOGS/chain.log
