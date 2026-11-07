#!/bin/bash
# Native (Linux) build of all engine objects into build/native (2 parallel jobs); writes OBJDONE / FAIL lines to $LOG
cd /home/user/ressim
LOG=${LOG:-/tmp/sd/build_err.log}
mkdir -p build/native
ls engine/src/*.cpp | xargs -P 2 -I{} bash -c 'f={}; b=$(basename $f .cpp); g++ -std=c++17 -O2 -fopenmp -fPIC -Wall -Wno-unused-function -Iengine/include -c $f -o build/native/$b.o 2>>'"$LOG"' || echo "FAIL $b" >> '"$LOG"
echo OBJDONE >> $LOG
# C33: shared library for the Python package (ctypes) - all objects are position independent
g++ -shared -fopenmp -o build/native/libressim.so build/native/*.o -lpthread -ldl 2>>$LOG || echo "FAIL libressim.so" >> $LOG
echo SODONE >> $LOG
