#!/usr/bin/env python3
"""C25 examples: FE geomechanics + fracture growth + fault mechanics in black oil (corner point), CO2 storage
(compositional) and EGS (thermal). Writes examples/GEOMECH_CP.DATA, CO2_GEOMECH.DATA, EGS_FAULT.DATA."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, '..', 'examples')


def wrap(vals, per=8):
    out, line = [], []
    for v in vals:
        line.append(v)
        if len(line) == per:
            out.append(' ' + ' '.join(line)); line = []
    if line: out.append(' ' + ' '.join(line))
    return '\n'.join(out)


def geomech_cp():
    nx, ny, nz, dx, dy, dz = 16, 10, 4, 50.0, 50.0, 10.0
    ifault, throw, dip = 8, 8.0, 0.04          # columns i >= 8 (0-based) are down-thrown by 8 m; 2.3 deg dip to +x

    def top(x, i):
        return 2000.0 + dip * x + (throw if i >= ifault else 0.0)
    coord = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x, y = i * dx, j * dy
            coord += ['%g %g %g %g %g %g' % (x, y, 1900, x, y, 2300)]
    z = []
    for k in range(nz):
        for side in (0, 1):                     # top / bottom of layer k
            for j in range(ny):
                for jy in (0, 1):
                    for i in range(nx):
                        for ix in (0, 1):
                            x = (i + ix) * dx
                            z.append('%.4f' % (top(x, i) + (k + side) * dz))
    n = nx * ny * nz
    deck = f"""-- C25 example: FE geomechanics on a CORNER-POINT grid in the black-oil engine, with in-run fracture growth and fault
-- mechanics driven by the finite-element stress.
--   grid      16 x 10 x 4 corner-point (COORD / ZCORN), 2.3 deg dip, fault F1 with an 8 m throw between columns 8 and 9
--   GEOMECH   two-way sequential FE geomechanics (isoparametric hexahedra on the corner-point nodes, 2 overburden, 1 under-
--             burden and 1 side-burden padding layers); stress-dependent permeability exp(-permStress dsigma'_m)
--   HFGROW    the hydraulic fracture HF1 through the injector grows by LEFM against the FE minimum principal stress
--             (FESTRESS) instead of a prescribed gradient
--   FAULTMECH Coulomb stress / seismicity on F1 with the FE stress change (STRESSFIELD ... FE), two-way fault permeability
-- Results: 3D view (SIGMIN, SIGEFFM, UPLIFT, STRAINV, FAILURE, KDAMAGE ...), Results > Faults & seismicity, run log.
-- OPM Flow / ECLIPSE ignore the '-- RESSIM:' lines (the deck runs there as a plain corner-point waterflood).
RUNSPEC
TITLE
Corner-point FE geomechanics + fracture growth + fault mechanics (C25, black oil)
DIMENS
 {nx} {ny} {nz} /
OIL
WATER
METRIC
TABDIMS
 1 1 20 20 /
WELLDIMS
 2 5 1 2 /
FAULTDIM
 1 /
START
 1 JAN 2026 /
GRID
SPECGRID
 {nx} {ny} {nz} 1 F /
COORD
{chr(10).join(' ' + c for c in coord)}
/
ZCORN
{wrap(z)}
/
PORO
 {n}*0.18 /
PERMX
 {n}*80 /
PERMY
 {n}*80 /
PERMZ
 {n}*8 /
FAULTS
 'F1' 8 8 1 10 1 4 'X' /
/
MULTFLT
 'F1' 0.1 /
/
-- RESSIM: GEOMECH ON 1 1e+10 0.25 0.7 1e-05 0 2 1 1 600 300 300 8e+09 0.25 1e-05 1 23000 0.38 2e+06 30 1e+06 4 0.0001 2e-08
-- RESSIM: HFRAC HF1 225 30 225 470 -1 -1 0.003 0 0.5 0 0
-- RESSIM: HFGROW HF1 1.0 380 0.1 60 1 FESTRESS
-- RESSIM: STRESSFIELD 0.23 0.65 1.05 60 0.5 0 FE
-- RESSIM: FAULTMECH F1 0.6 0 0 100 5 1.0 -1.0 10
PROPS
PVCDO
 200 1.02 5e-5 2.0 0 /
PVTW
 200 1.0 4.5e-5 0.4 0 /
DENSITY
 850 1020 1 /
ROCK
 200 5e-5 /
SWOF
 0.20 0.00 1.00 0
 0.30 0.02 0.70 0
 0.40 0.06 0.45 0
 0.50 0.13 0.26 0
 0.60 0.22 0.12 0
 0.70 0.34 0.04 0
 0.80 0.50 0.00 0
 1.00 1.00 0.00 0 /
SOLUTION
EQUIL
 2020 205 2400 0 /
SCHEDULE
WELSPECS
 INJ G 5 5 2010 WATER /
 PROD G 14 5 2040 OIL /
/
COMPDAT
 INJ 5 5 1 4 OPEN 1* 1* 0.2 /
 PROD 14 5 1 4 OPEN 1* 1* 0.2 /
/
WCONINJE
 INJ WATER OPEN RATE 1500 1* 420 /
/
WCONPROD
 PROD OPEN BHP 5* 150 /
/
TSTEP
 6*5 6*20 6*50 /
END
"""
    return deck


def co2_geomech():
    src = open(os.path.join(EX, 'CO2STORE.DATA')).read()
    head = """-- C25 example: CO2 storage (compositional CO2STORE core) with FE geomechanics, fault mechanics and fracture growth.
--   GEOMECH   two-way sequential FE geomechanics: pore volume and permeability of the CO2STORE cells follow the FE
--             porosity / stress after every report step (mole-conservative); stress, strain, uplift and failure maps
--   FAULTMECH a vertical fault F1 at x = 50 m (between columns 10 and 11) with the FE stress change (STRESSFIELD ... FE):
--             Coulomb stress, Dieterich seismicity, two-way fault permeability
--   HFGROW    the injection fracture HF1 (bottom layers) grows towards the fault by LEFM against the FE minimum stress\n--             (tight 2 mD formation, 40000 sm3/d, BHP limit 230 bar; poro-elastic stressing arrests the growth)
-- OPM Flow ignores the '-- RESSIM:' lines.
"""
    src = src.replace('-- CO2 storage in a brine aquifer (OPM Flow CO2STORE format, after opm-tests co2store/CO2STORE.DATA)',
                      head + '-- base: CO2 storage in a brine aquifer (OPM Flow CO2STORE format, after opm-tests co2store/CO2STORE.DATA)', 1)
    src = src.replace('CO2 storage - vertical section', 'CO2 storage + FE geomechanics + fault + fracture growth (C25)', 1)
    src = src.replace('WELLDIMS\n 1 20 1 1 /', 'WELLDIMS\n 1 20 1 1 /\nFAULTDIM\n 1 /', 1)
    extra = """FAULTS
 'F1' 10 10 1 1 1 20 'X' /
/
-- RESSIM: GEOMECH ON 1 5e+09 0.25 1 1e-05 0 2 1 1 500 300 300 5e+09 0.25 1e-05 1 22600 0.5 1e+06 30 5e+05 4 0.0001 1e-08
-- RESSIM: HFRAC HF1 0.5 50 60 50 1085 1100 0.002 0 0.5 0 0
-- RESSIM: HFGROW HF1 0.5 160 0.1 6 1 FESTRESS
-- RESSIM: STRESSFIELD 0.226 0.7 0.95 30 0.5 0 FE
-- RESSIM: FAULTMECH F1 0.6 0 0 100 5 1.0 -1.0 5

PROPS"""
    src = src.replace('\nPROPS', '\n' + extra, 1)
    # k0 0.5 -> S_h ~ 177 bar at 1090 m; the injector (BHP limit 230 bar) exceeds it, poro-elastic stressing raises S_h with p
    n = src.count("'INJ' 'GAS' 'OPEN' 'RATE' 20000 1* 200 /")
    assert n == 1, n
    src = src.replace("'INJ' 'GAS' 'OPEN' 'RATE' 20000 1* 200 /", "'INJ' 'GAS' 'OPEN' 'RATE' 40000 1* 230 /")
    # tight saline formation (2 mD): the injection pressure reaches the FE closure stress near the well
    for a, b in (('PERMX\n 400*2000 /', 'PERMX\n 400*2 /'), ('PERMY\n 400*2000 /', 'PERMY\n 400*2 /'), ('PERMZ\n 400*200 /', 'PERMZ\n 400*0.2 /')):
        assert src.count(a) == 1, a
        src = src.replace(a, b)
    src = src.replace("'INJ' 'GAS' 'SHUT' 'RATE' 0 1* 200 /", "'INJ' 'GAS' 'SHUT' 'RATE' 0 1* 230 /")
    return src


def egs_fault():
    src = open(os.path.join(EX, 'EGS_EDFM.DATA')).read()
    head = """-- C25 example: EGS doublet (thermal core) with THM geomechanics and fault mechanics on a fault between the wells.
--   GEOMECH   thermo-poro-elastic FE mechanics (fixed-stress, iterative) of the thermal engine
--   FAULTMECH fault F1 (x-faces between columns 11 and 12) with the FE stress change of the THM mechanics
--             (STRESSFIELD ... FE: cooling contracts the rock -> the horizontal stress drops -> Coulomb stress rises);
--             TH E alpha nu: thermo-elastic stress used when the analytic source is selected
--   Coulomb stress, Dieterich seismicity, catalogue, two-way fault permeability, post-run rate-and-state trace (FAULTRSF)
"""
    src = head + src
    src = src.replace('Thermal core: EGS doublet in granite with embedded discrete fractures (EDFM + Barton-Bandis + hydro-shear)',
                      'EGS doublet + THM geomechanics + fault mechanics (C25, thermal)', 1)
    src = src.replace('\nMETRIC\n', '\nMETRIC\nFAULTDIM\n 1 /\n', 1)
    extra = """INIT

FAULTS
 'F1' 11 11 1 11 1 3 'X' /
/
-- RESSIM: GEOMECH ON 1 3e+10 0.25 0.7 8e-06 0 1 1 0 500 500 300 3e+10 0.25 8e-06 1 26000 0.7 5e+06 32 2e+06 4 0.0001
-- RESSIM: STRESSFIELD 0.26 0.62 0.98 35 0.5 0 FE TH 30 8e-06 0.25
-- RESSIM: FAULTMECH F1 0.6 0 0 100 2 1.0 -1.0 5
-- RESSIM: FAULTRSF F1 2 0.010 0.015 1.0
"""
    assert '\nINIT\n' in src
    src = src.replace('\nINIT\n', '\n' + extra, 1)
    return src


if __name__ == '__main__':
    for name, fn in (('GEOMECH_CP.DATA', geomech_cp), ('CO2_GEOMECH.DATA', co2_geomech), ('EGS_FAULT.DATA', egs_fault)):
        with open(os.path.join(EX, name), 'w') as f:
            f.write(fn())
        print('wrote', name)
