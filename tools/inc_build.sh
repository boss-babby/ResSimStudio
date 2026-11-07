#!/bin/bash
# incremental native build: objects older than their source or the newest header (1 job at a time to fit 2 GB)
cd /home/user/ressim
LOG=/tmp/sd/build_err.log
H=$(ls -t engine/include/ressim/*.h | head -1)
J=${JOBS:-2}
for f in engine/src/*.cpp; do b=$(basename $f .cpp); o=build/native/$b.o; if [ ! -f $o ] || [ $o -ot $H ] || [ $o -ot $f ]; then echo $f; fi; done > /tmp/sd/stale.txt
cat /tmp/sd/stale.txt | xargs -P $J -I{} bash -c 'f={}; b=$(basename $f .cpp); g++ -std=c++17 -O2 -fopenmp -fPIC -Wall -Wno-unused-function -Iengine/include -c $f -o build/native/$b.o.tmp 2>>'"$LOG"' && mv build/native/$b.o.tmp build/native/$b.o || echo "FAIL $b" >> '"$LOG"
echo OBJDONE >> $LOG
g++ -shared -fopenmp -o build/native/libressim.so build/native/*.o -lpthread -ldl 2>>$LOG || echo "FAIL libressim.so" >> $LOG
echo SODONE >> $LOG
