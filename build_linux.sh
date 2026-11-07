#!/bin/bash
# C39 native build of ResSim Studio on Linux / macOS (CMake + GLFW + OpenGL 3.3 core).
#   Debian / Ubuntu : sudo apt install build-essential cmake libglfw3-dev libgl-dev zenity
#   Fedora          : sudo dnf install gcc-c++ cmake glfw-devel mesa-libGL-devel zenity
#   Arch            : sudo pacman -S base-devel cmake glfw mesa zenity
#   macOS (Homebrew): brew install cmake glfw libomp   (OpenGL.framework is part of the OS)
#   ./build_linux.sh            -> build/linux/ResSimStudio, ResGeo (symlink), ResSimCLI, ressim_tests
# Env: JOBS (default: number of cores), TYPE (Release), AVX2 (ON)
set -eu
cd "$(dirname "$0")"
JOBS=${JOBS:-$( (nproc || sysctl -n hw.ncpu) 2>/dev/null || echo 2)}
TYPE=${TYPE:-Release}
AVX2=${AVX2:-ON}
EXTRA=()
if [ "$(uname)" = "Darwin" ]; then
  AVX2=OFF   # universal / Apple silicon builds
  if command -v brew >/dev/null && [ -d "$(brew --prefix libomp 2>/dev/null)" ]; then
    OMP=$(brew --prefix libomp)
    EXTRA+=(-DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$OMP/include" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=$OMP/lib/libomp.dylib)
  fi
fi
cmake -S . -B build/linux -DCMAKE_BUILD_TYPE=$TYPE -DRESSIM_AVX2=$AVX2 "${EXTRA[@]}"
cmake --build build/linux -j "$JOBS"
echo "BUILD OK: build/linux/ResSimStudio  (file dialogs use zenity / kdialog / osascript when available)"
