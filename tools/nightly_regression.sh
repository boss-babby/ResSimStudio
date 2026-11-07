#!/bin/bash
# C33 nightly deck regression (Linux / cluster). crontab:  0 2 * * *  /path/to/ressim/tools/nightly_regression.sh
set -u
HOME_DIR=${RESSIM_HOME:-$(cd "$(dirname "$0")/.." && pwd)}
TESTS=${RESSIM_TESTS:-$HOME_DIR/opm-tests}
CLI=${RESSIM_CLI:-$HOME_DIR/ResSimCLI}
STAMP=$(date +%Y%m%d)
mkdir -p "$HOME_DIR/regression/$STAMP"
python3 "$HOME_DIR/tools/fetch_opm_tests.py" -o "$TESTS" --update
"$CLI" --regress "$TESTS/manifest.txt" -o "$HOME_DIR/regression/$STAMP" --history "$HOME_DIR/regression/history.tsv" | tee "$HOME_DIR/regression/$STAMP/run.log"
