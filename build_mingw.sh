#!/bin/bash
# Cross-build ResSimStudio.exe with MinGW-w64 GCC (posix threads, OpenMP) on Linux.
#   sudo apt install g++-mingw-w64-x86-64-posix
#   ./build_mingw.sh            -> build/mingw/ResSimStudio.exe (fully static, no DLLs needed)
# Env: JOBS (default 2), CXX (default x86_64-w64-mingw32-g++-posix)
set -u
cd "$(dirname "$0")"
CXX=${CXX:-x86_64-w64-mingw32-g++-posix}
JOBS=${JOBS:-2}
OUT=build/mingw
mkdir -p $OUT/obj
INC="-Iengine/include -Ithird_party/imgui -Ithird_party/imgui/backends -Ithird_party/implot -Ithird_party/textedit -Iapp/src"
DEFS="-DNOMINMAX -DWIN32_LEAN_AND_MEAN -D_USE_MATH_DEFINES"
BASE="-std=c++17 -pipe -march=x86-64-v2 -mtune=generic $DEFS $INC"

flags_for() {
  case "$1" in
    engine/*)                 echo "-O3 -fopenmp -ffp-contract=off" ;;
    third_party/implot/implot_items.cpp) echo "-O1 --param ggc-min-expand=20 --param ggc-min-heapsize=65536" ;;
    third_party/*)            echo "-O2" ;;
    app/*)                    echo "-O2 -fopenmp -Wall -Wno-unknown-pragmas" ;;
  esac
}

SRCS="$(ls engine/src/*.cpp) $(ls app/src/*.cpp) $(ls third_party/imgui/*.cpp) $(ls third_party/imgui/backends/*.cpp) third_party/implot/implot.cpp third_party/textedit/TextEditor.cpp"
fail=0
compile() {
  local f=$1 o=$OUT/obj/$(echo $1 | tr '/' '_' | sed 's/\.cpp$/.o/')
  if [ -f "$o" ] && [ "$o" -nt "$f" ] && [ -z "$(find app/src engine/include -newer "$o" -name '*.h' | head -1)" ]; then return 0; fi
  echo "  CXX $f"
  $CXX $BASE $(flags_for $f) -c "$f" -o "$o" || { echo "FAILED: $f"; return 1; }
}
running=0
for f in $SRCS; do
  compile $f &
  running=$((running+1))
  if [ $running -ge $JOBS ]; then wait -n || fail=1; running=$((running-1)); fi
done
while [ $running -gt 0 ]; do wait -n || fail=1; running=$((running-1)); done
# the heaviest translation unit is compiled alone to keep memory low
compile third_party/implot/implot_items.cpp || fail=1
[ $fail = 0 ] || { echo "BUILD FAILED"; exit 1; }

# embedded resources: app icon, version info, keyword manual, icon font
WINDRES=${WINDRES:-x86_64-w64-mingw32-windres}
if [ ! -f $OUT/obj/zz_resources.o ] || [ -n "$(find app/res -newer $OUT/obj/zz_resources.o | head -1)" ]; then
  echo "  RC  app/res/ressim.rc"
  $WINDRES --include-dir app/res -O coff app/res/ressim.rc -o $OUT/obj/zz_resources.o || { echo "RC FAILED"; exit 1; }
fi
echo "  LINK $OUT/ResSimStudio.exe"
$CXX -o $OUT/ResSimStudio.exe $OUT/obj/*.o -mwindows -static -static-libgcc -static-libstdc++ -fopenmp -s \
  -ld3d11 -ldxgi -ld3dcompiler_47 -ld2d1 -ldwrite -ldwmapi -lcomdlg32 -lshell32 -lole32 -lwindowscodecs -luuid -luser32 -lgdi32 -limm32 \
  || { echo "LINK FAILED"; exit 1; }

# console test runner (engine regression tests)
if [ "${TESTS:-1}" = 1 ]; then
  echo "  LINK $OUT/ressim_tests.exe"
  $CXX $BASE -O2 -fopenmp tests/run_tests.cpp $OUT/obj/engine_src_*.o -o $OUT/ressim_tests.exe -static -s || echo "tests link failed"
fi
ls -la $OUT/*.exe
echo "BUILD OK"
