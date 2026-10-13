#!/bin/bash
# Runs every examples/*.DATA through OPM Flow into build/opm/<ID>/ (reference results for compare_opm).
cd "$(dirname "$0")/.."
for f in ${@:-examples/*.DATA}; do
  id=$(basename $f .DATA)
  d=build/opm/$id; rm -rf $d; mkdir -p $d; cp $f $d/
  start=$(date +%s%N)
  (cd $d && OMP_NUM_THREADS=1 timeout 900 flow $id.DATA --output-dir=. --threads-per-process=1 > flow.log 2>&1)
  rc=$?
  t=$(( ($(date +%s%N) - start) / 1000000 ))
  echo "$id rc=$rc ${t}ms $(grep -E '^Error|aborted' $d/flow.log | head -2 | tr '\n' ' ')"
done
