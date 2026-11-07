#!/bin/bash
# GUI smoke test under Xvfb + wine: smoke_app.sh [deck]
export WINEPREFIX=$HOME/.cache/wp WINEDEBUG=${WDBG:--all}
Xvfb :7 -screen 0 1600x1000x24 >/dev/null 2>&1 &
sleep 2
export DISPLAY=:7
cd ~/ressim
if [ -n "$1" ]; then exec wine explorer /desktop=rs,1580x960 build/mingw/ResSimStudio.exe "$1"; else exec wine explorer /desktop=rs,1580x960 build/mingw/ResSimStudio.exe; fi
