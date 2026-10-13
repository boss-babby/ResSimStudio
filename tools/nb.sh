#!/bin/bash
# Native (Linux) build of the engine objects -> build/native/*.o (incremental). Prints "NB OK".
cd "$(dirname "$0")/.."
mkdir -p build/native
fail=0
for f in engine/src/*.cpp; do
  o=build/native/$(basename ${f%.cpp}).o
  if [ -f "$o" ] && [ "$o" -nt "$f" ] && [ -z "$(find engine/include engine/src -name '*.h' -newer "$o" -o -name '*.inc' -newer "$o" | head -1)" ]; then continue; fi
  echo "  CXX $f"
  g++ -std=c++17 -O2 -fopenmp -Wall -Wno-unused-function -Iengine/include -c $f -o $o || fail=1
done
[ $fail = 0 ] && echo "NB OK" || echo "NB FAILED"
