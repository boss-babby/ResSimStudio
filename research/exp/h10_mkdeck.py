#!/usr/bin/env python3
"""H10 (Research Cycle 9) deck generator: CO2 storage (ResSim CO2STORE core) or a water-injection control,
2-D vertical section with a fault monitored by FAULTMECH (Dieterich seismicity, one-way, analytic poroelastic stress).

  mkdeck(case, u, dt_report=10.0) -> deck text
    case 'N'  CO2, fault at x = 150 m (plume reaches the fault)
    case 'F'  CO2, fault at x = 300 m (pressure only)
    case 'W'  water injection (black oil, initially water-filled), fault at x = 150 m: the near-linear control
    case 'WF' water control with the fault at x = 300 m
    case 'NN' / 'WNN' (H10c) fault at x = 60 m; 'NH' / 'WNH' (H10c) fault at x = 150 m, 3x rate
    u         K = 4 rate multipliers (mean 1) for 90 d intervals, followed by 360 d of shut-in
"""
import sys

NX, NZ, DX, DY, DZ = 40, 10, 10.0, 100.0, 5.0
Q_CO2 = 40000.0      # sm3/d mean CO2 rate
Q_W = 150.0          # sm3/d mean water rate (~ the CO2 reservoir volume rate at ~110 bar / 50 C)
TK, TSHUT = 90.0, 360.0


def _multpv():
    row = ["%d*1" % (NX - 1), "1e4"]
    return " ".join(row * NZ)


def _schedule(u, dt, phase, q):
    out = []
    for uk in u:
        n = int(round(TK / dt))
        out.append("WCONINJE\n 'INJ' '%s' 'OPEN' 'RATE' %.6f 1* 400 /\n/\nTSTEP\n %d*%g /" % (phase, uk * q, n, dt))
    n = int(round(TSHUT / dt))
    out.append("WCONINJE\n 'INJ' '%s' 'SHUT' 'RATE' 0 1* 400 /\n/\nTSTEP\n %d*%g /" % (phase, n, dt))
    return "\n".join(out)


CASES = {  # case -> (fault column, water?, rate multiplier)
    "N": (15, False, 1.0), "F": (30, False, 1.0), "W": (15, True, 1.0), "WF": (30, True, 1.0),
    "NN": (6, False, 1.0), "WNN": (6, True, 1.0), "NH": (15, False, 3.0), "WNH": (15, True, 3.0),   # H10c exploratory
    # H10e (Cycle 10): sealing fault (MULTFLT 0.02) reopened by two-way coupling (FAULTMECH permFactor 50 -> T x 1.0)
    "NS": (15, False, 1.0, True), "NHS": (15, False, 3.0, True), "WNS": (15, True, 1.0, True), "WNHS": (15, True, 3.0, True)}
# H10e amendment A1 (pre-check, u = 1 only): the default fault is 30-37 bar from failure and never reactivates. A uniform shear
# pre-load c (bar, entered as cohesion -c; slip is off so cohesion only shifts CFF) puts the most-stressed face at CFF = +2 bar at u = 1.
PRELOAD = {"NS": 32.0, "NHS": 23.0, "WNS": 27.0, "WNHS": 14.0, "NSC": 32.0, "WNSC": 27.0}
# H10e exploratory controls: same sealing fault and pre-load, but no permeability switching (one-way, permFactor 1)
CASES.update({"NSC": (15, False, 1.0, True), "WNSC": (15, True, 1.0, True)})
NOSWITCH = {"NSC", "WNSC"}
# Cycle 11 (CYCLE11_PLAN.md): H10f-a smoothed ramp (FAULTMECH item 10 = 0.5 instead of 0.05), H10f-b reactivation constraint
# (H10e switching, pre-load 3 bar below A1 so that u = 1 stays below failure)
CASES.update({"NSR": (15, False, 1.0, True), "WNSR": (15, True, 1.0, True), "NSB": (15, False, 1.0, True), "WNSB": (15, True, 1.0, True)})
PRELOAD.update({"NSR": 32.0, "WNSR": 27.0, "NSB": 29.0, "WNSB": 24.0})
RAMP = {"NSR": 0.5, "WNSR": 0.5}
# Cycle 11 EXPLORATORY (after H10f-a): dose-response in the ramp width (log-slope of the T multiplier vs CFF)
CASES.update({"NSR15": (15, False, 1.0, True), "WNSR15": (15, True, 1.0, True)})
PRELOAD.update({"NSR15": 32.0, "WNSR15": 27.0})
RAMP.update({"NSR15": 0.15, "WNSR15": 0.15})
CASES.update({"NSR08": (15, False, 1.0, True), "WNSR08": (15, True, 1.0, True)})
PRELOAD.update({"NSR08": 32.0, "WNSR08": 27.0})
RAMP.update({"NSR08": 0.08, "WNSR08": 0.08})


def _register(case):
    """Cycle 12 (H10g) generic case 'X:<fault col>:<water 0/1>:<pre-load bar>:<ramp w>:<permFactor>' (sealing fault).
    permFactor 1 = no switching (one-way control). Existing named cases are unchanged."""
    if case in CASES or not case.startswith("X:"):
        return
    _, col, wat, pre, w, pf = case.split(":")
    CASES[case] = (int(col), wat == "1", 1.0, True)
    PRELOAD[case] = float(pre)
    if float(pf) == 1.0:
        NOSWITCH.add(case)
    else:
        RAMP[case] = float(w)


def mkdeck(case, u, dt_report=10.0):
    _register(case)
    fx, water, qm = CASES[case][:3]
    seal = len(CASES[case]) > 3 and CASES[case][3]
    n = NX * NZ
    head = ["-- H10 (Research Cycle 9) case %s, u = %s" % (case, " ".join("%.4f" % x for x in u)),
            "RUNSPEC", "DIMENS\n %d 1 %d /" % (NX, NZ)]
    if water:
        head += ["OIL", "WATER"]
    else:
        head += ["OIL", "GAS", "CO2STORE"]
    head += ["METRIC", "START\n 1 'JAN' 2025 /", "WELLDIMS\n 1 10 1 1 /", "FAULTDIM\n 1 /", "TABDIMS\n 1 1 40 20 /",
             "GRID",
             "DX\n %d*%g /" % (n, DX), "DY\n %d*%g /" % (n, DY), "DZ\n %d*%g /" % (n, DZ), "TOPS\n %d*1000 /" % NX,
             "PORO\n %d*0.25 /" % n, "PERMX\n %d*50 /" % n, "PERMY\n %d*50 /" % n, "PERMZ\n %d*5 /" % n,
             "MULTPV\n %s /" % _multpv(),
             "FAULTS\n 'F1' %d %d 1 1 1 %d 'X' /\n/" % (fx, fx, NZ),
             "-- RESSIM: STRESSFIELD 0.226 0.7 0.95 30 0.5",
             ("MULTFLT\n 'F1' 0.02 /\n/" if seal else "-- (fault not sealing)"),
             "-- RESSIM: FAULTMECH F1 0.6 %g 1.0 100 5 1.0 -1.0" % (-PRELOAD.get(case, 0.0) if seal else 0.0) + (" %g" % (float(case.split(":")[5]) if case.startswith("X:") else 50) if seal and case not in NOSWITCH else "") + (" %g" % RAMP[case] if case in RAMP else ""),
             "PROPS", "ROCK\n 110 5e-5 /"]
    if water:
        head += ["SWOF\n 0 0 1 0\n 1 1 0 0\n/",
                 "PVTW\n 110 1.0 4.5e-5 0.5 0 /",
                 "PVDO\n 50 1.0 1.0\n 400 0.99 1.0\n/",
                 "DENSITY\n 800 1030 1 /",
                 "SOLUTION", "PRESSURE\n %d*110 /" % n, "SWAT\n %d*1 /" % n]
    else:
        head += ["SGOF\n 0 0 1 0\n 0.05 0 0.85 0\n 0.2 0.03 0.55 0\n 0.4 0.15 0.25 0\n 0.6 0.35 0.07 0\n 0.8 0.6 0 0\n 1.0 1.0 0 0\n/",
                 "SALINITY\n 0.7 /",
                 "SOLUTION", "EQUIL\n 1000 100 1100 0 900 0 1 1 0 /",
                 "RTEMPVD\n 1000 50\n 1100 50\n/", "RSVD\n 1000 0\n 1100 0\n/"]
    head += ["SUMMARY", "FPR", "SCHEDULE",
             "WELSPECS\n 'INJ' 'G1' 1 1 1* '%s' /\n/" % ("WATER" if water else "GAS"),
             "COMPDAT\n 'INJ' 1 1 8 10 'OPEN' 1* 1* 0.2 /\n/",
             _schedule(u, dt_report, "WATER" if water else "GAS", qm * (Q_W if water else Q_CO2)), "END", ""]
    return "\n".join(head)


if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "N"
    u = [float(x) for x in sys.argv[2:6]] if len(sys.argv) > 5 else [1, 1, 1, 1]
    sys.stdout.write(mkdeck(case, u))
