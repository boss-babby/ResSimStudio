#!/bin/bash
# xclick.sh x y [out.png]  - slow click in the Xvfb :7 smoke session, then screenshot
export DISPLAY=:7
xdotool mousemove $1 $2; sleep 0.6; xdotool mousedown 1; sleep 0.3; xdotool mouseup 1; sleep ${WAIT:-2}
[ -n "$3" ] && import -window root "$3" && sync && sleep 1
