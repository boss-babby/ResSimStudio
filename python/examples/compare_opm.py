"""Run a deck with ResSim and compare with an OPM Flow / ECLIPSE summary of the same deck."""
import sys
import ressim

if len(sys.argv) < 3:
    sys.exit("usage: compare_opm.py CASE.DATA OPM_RUN/CASE.SMSPEC")
s = ressim.Session(sys.argv[1])
s.run()
ref = ressim.open_summary(sys.argv[2])
for d in ressim.compare(ref, s, ["FOPR", "FOPT", "FWCT", "FGOR", "FPR"]):
    print(f"{d['key']:10s} ref {d['ref_final']:12.5g}  ressim {d['cmp_final']:12.5g}  rel RMS {100 * d['rel_rms']:6.2f} %")
s.report("comparison.html")
