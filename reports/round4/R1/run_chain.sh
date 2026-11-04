#!/bin/bash
# runs E23, E18, E25 gates, E19 pools p1 one after another in the worktree, one core each, timing each
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=/home/user/wt-repro/src
source /home/user/.venv-bapat/bin/activate
cd /home/user/wt-repro
LOGS=/home/user/anytime-valid-volatility-monitor/reports/round4/R1/logs
echo "chain start $(date -u)" > $LOGS/chain.log
run() { name=$1; shift; echo "START $name $(date -u)" >> $LOGS/chain.log; t0=$(date +%s); python "$@" > $LOGS/$name.log 2>&1; rc=$?; echo "END $name rc=$rc wall=$(( $(date +%s) - t0 ))s $(date -u)" >> $LOGS/chain.log; }
run e23_run experiments/e23_fresh_releases.py run --procs=1
run e18_run experiments/e18_landmark_events.py run figure context awake
run e25_gates experiments/e25_sequential_events.py gates --procs=1
run e19_pools_p1 experiments/e19_real_null.py pools p1
echo CHAIN_DONE >> $LOGS/chain.log
