#!/bin/bash
# C36 headless check of the Geomodel 3-D scene builders / seismic interpretation UI (Linux, no D3D).
# Needs the native engine objects (tools/rebuild_native.sh -> build/native/*.o). Run from the repository root.
set -e
B=${TMPDIR:-/tmp}/geo3d_headless
mkdir -p "$B"
cp tests/geo3d_headless/app.h tests/geo3d_headless/geoscene3d.h app/src/geo_state.h app/src/ui_geo3d.cpp "$B/"
if [ ! -f "$B/imgui_all.a" ]; then
  for f in third_party/imgui/imgui.cpp third_party/imgui/imgui_draw.cpp third_party/imgui/imgui_tables.cpp third_party/imgui/imgui_widgets.cpp \
           third_party/implot/implot.cpp third_party/implot/implot_items.cpp app/src/colormaps.cpp; do
    g++ -std=c++17 -O1 -Ithird_party/imgui -Ithird_party/implot -Iapp/src -c "$f" -o "$B/$(basename "$f" .cpp).o"
  done
  ar rcs "$B/imgui_all.a" "$B"/imgui*.o "$B"/implot*.o "$B"/colormaps.o
fi
g++ -std=c++17 -O1 -fopenmp -Wall -Wno-unused-function -I"$B" -Iapp/src -Iengine/include -Ithird_party/imgui -Ithird_party/implot \
    tests/geo3d_headless/driver.cpp $(ls build/native/*.o) "$B/imgui_all.a" -o "$B/geo3d_headless" -lpthread -ldl
OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} "$B/geo3d_headless"
