#!/usr/bin/env python3
"""H12 (Research Cycle 13, CYCLE13_PLAN.md) deck generator: seasonal cyclic storage with signed rates.

  mkdeck(case, u, dt_report=5.0, asig=1.0) -> deck text
    case 'SW'  water storage: initially water-filled 2-D section (as H10 'W'), injector INJ and producer PRD in column 1
    case 'SG'  depleted dry-gas field storage: gas-filled (Sg = 1 - Swc) black-oil section, gas INJ / PRD in column 1
    u          8 signed multipliers for 45-d intervals (u > 0 inject, u < 0 withdraw, u = 0 idle), then 180 d rest
    unit       'pulse:k' gives u = e_k (used for the H12b response matrix)
Fault F1 at column 15 (x = 150 m), not sealing, FAULTMECH one-way (no permeability switching), Dieterich Asigma = asig bar.
"""
import sys

NX, NZ, DX, DY, DZ = 40, 10, 10.0, 100.0, 5.0
Q_W = 150.0          # sm3/d water (as H10)
Q_G = 12000.0        # sm3/d gas ~ 150 rm3/d at ~110 bar (Bg ~ 0.0121 rm3/sm3)
TK, TREST, K = 45.0, 180.0, 8


def _multpv():
    return " ".join(["%d*1" % (NX - 1), "1e4"] * NZ)


def _schedule(u, dt, gas):
    q, ph = (Q_G, "GAS") if gas else (Q_W, "WATER")
    ctl = "GRAT" if gas else "WRAT"
    out = []
    for uk in list(u) + [None]:
        n = int(round((TK if uk is not None else TREST) / dt))
        if uk is not None and uk > 0:
            out.append("WCONINJE\n 'INJ' '%s' 'OPEN' 'RATE' %.6f 1* 400 /\n/\nWCONPROD\n 'PRD' 'SHUT' '%s' 3* 0 /\n/" % (ph, uk * q, ctl))
        elif uk is not None and uk < 0:
            rates = ("1* 1* %.6f" % (-uk * q)) if gas else ("1* %.6f 1*" % (-uk * q))
            out.append("WCONINJE\n 'INJ' '%s' 'SHUT' 'RATE' 0 1* 400 /\n/\nWCONPROD\n 'PRD' 'OPEN' '%s' %s 1* 1* 20 /\n/" % (ph, ctl, rates))
        else:
            out.append("WCONINJE\n 'INJ' '%s' 'SHUT' 'RATE' 0 1* 400 /\n/\nWCONPROD\n 'PRD' 'SHUT' '%s' 3* 0 /\n/" % (ph, ctl))
        out.append("TSTEP\n %d*%g /" % (n, dt))
    return "\n".join(out)


def mkdeck(case, u, dt_report=5.0, asig=1.0):
    gas = case == "SG"
    if isinstance(u, str) and u.startswith("pulse:"):
        k = int(u.split(":")[1]); u = [1.0 if i == k else 0.0 for i in range(K)]
    n = NX * NZ
    d = ["-- H12 (Research Cycle 13) case %s, u = %s" % (case, " ".join("%.4f" % x for x in u)),
         "RUNSPEC", "DIMENS\n %d 1 %d /" % (NX, NZ), "OIL", "WATER"] + (["GAS"] if gas else []) + [
         "METRIC", "START\n 1 'JAN' 2025 /", "WELLDIMS\n 2 10 1 2 /", "FAULTDIM\n 1 /", "TABDIMS\n 1 1 40 20 /",
         "GRID", "DX\n %d*%g /" % (n, DX), "DY\n %d*%g /" % (n, DY), "DZ\n %d*%g /" % (n, DZ), "TOPS\n %d*1000 /" % NX,
         "PORO\n %d*0.25 /" % n, "PERMX\n %d*50 /" % n, "PERMY\n %d*50 /" % n, "PERMZ\n %d*5 /" % n,
         "MULTPV\n %s /" % _multpv(),
         "FAULTS\n 'F1' 15 15 1 1 1 %d 'X' /\n/" % NZ,
         "-- RESSIM: STRESSFIELD 0.226 0.7 0.95 30 0.5",
         "-- RESSIM: FAULTMECH F1 0.6 0 %g 100 5 1.0 -1.0" % asig,
         "PROPS", "ROCK\n 110 5e-5 /",
         "SWOF\n 0.2 0 1 0\n 1 1 0 0\n/" if gas else "SWOF\n 0 0 1 0\n 1 1 0 0\n/",
         "PVTW\n 110 1.0 4.5e-5 0.5 0 /", "PVDO\n 50 1.0 1.0\n 400 0.99 1.0\n/"]
    if gas:
        d += ["PVDG\n 50 0.0250 0.013\n 100 0.0125 0.015\n 200 0.0062 0.018\n 300 0.0045 0.022\n 500 0.0030 0.028 /",
              "SGOF\n 0 0 1 0\n 0.1 0.02 0.7 0\n 0.4 0.3 0.2 0\n 0.8 0.9 0 0 /",
              "DENSITY\n 800 1030 0.8 /",
              "SOLUTION", "PRESSURE\n %d*110 /" % n, "SWAT\n %d*0.2 /" % n, "SGAS\n %d*0.8 /" % n]
    else:
        d += ["DENSITY\n 800 1030 1 /", "SOLUTION", "PRESSURE\n %d*110 /" % n, "SWAT\n %d*1 /" % n]
    ph = "GAS" if gas else "WATER"
    d += ["SUMMARY", "FPR", "SCHEDULE",
          "WELSPECS\n 'INJ' 'G1' 1 1 1* '%s' /\n 'PRD' 'G1' 1 1 1* '%s' /\n/" % (ph, ph),
          "COMPDAT\n 'INJ' 1 1 8 10 'OPEN' 1* 1* 0.2 /\n 'PRD' 1 1 8 10 'OPEN' 1* 1* 0.2 /\n/",
          _schedule(u, dt_report, gas), "END", ""]
    return "\n".join(d)


if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "SW"
    u = [float(x) for x in sys.argv[2:10]] if len(sys.argv) > 9 else [1, 1, 1, 1, -1, -1, -1, -1]
    sys.stdout.write(mkdeck(case, u))
