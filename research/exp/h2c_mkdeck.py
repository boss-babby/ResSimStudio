#!/usr/bin/env python3
"""H2c helper: HFGROW_INJECTION at along-fracture cell size dx (10 / 5 / 2.5 m) in two new settings
   GAS : dry-gas injection (GAS phase, PVDG / SGOF; same reservoir-volume rate x gasMult)
   CP  : corner-point grid (COORD / ZCORN) with a 2 degree dip along y, water injection
-> build/h3/HC_<variant>_<dx>[_SG|_SGL].DATA   usage: h2c_mkdeck.py dx [SG|SGL|SGLC|-] GAS|CP [gasRate]"""
import sys, math, os
dx = float(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else '-'; var = sys.argv[3] if len(sys.argv) > 3 else 'GAS'
gasRate = float(sys.argv[4]) if len(sys.argv) > 4 else 6000
sub = mode in ('SG', 'SGL', 'SGLC'); leak = mode in ('SGL', 'SGLC')   # C48 H2e: SGLC = gas-aware sub-cell storage
src = open('examples/HFGROW_INJECTION.DATA').read()
nx = int(round(810 / dx)); nx += (nx + 1) % 2
ic = nx // 2 + 1; xc = (ic - 0.5) * dx; n = nx * 31
dyTxt = src.split('DY\n')[1].split('/')[0]
dyRow = []
for tok in dyTxt.split():
    c, v = tok.split('*'); dyRow += [float(v)] * (int(c) // 81)
assert len(dyRow) == 31
out = src.replace(' 81 31 1 /', f' {nx} 31 1 /').replace('2511*10 /\nDY', f'{n}*{dx:g} /\nDY', 1)
out = out.split('DY\n')[0] + 'DY\n' + dyTxt.replace('81*', f'{nx}*') + '/' + out.split('DY\n')[1].split('/', 1)[1]
out = out.replace('2511*', f'{n}*')
out = out.replace('HFRAC F1 305 688.793902 505 688.793902', f'HFRAC F1 {xc - 100:g} 688.793902 {xc + 100:g} 688.793902')
out = out.replace('HFGROW F1 1.0 340 0 20', 'HFGROW F1 1.0 340 0 20 1' + (' SUBGRID' if sub else '') + (' LEAK' if leak else '') + ('C' if mode == 'SGLC' else ''))
out = out.replace(' P G 41 16 ', f' P G {ic} 16 ').replace(' P 41 16 1 1 ', f' P {ic} 16 1 1 ')
if var == 'GAS':
    out = out.replace('OIL\nWATER\n', 'OIL\nWATER\nGAS\n', 1)
    out = out.replace('DENSITY\n 800 1000 1 /', 'DENSITY\n 800 1000 0.9 /\nPVDG\n 100 0.0125 0.015\n 300 0.0045 0.022\n 500 0.0030 0.028 /\n'
                      'SGOF\n 0 0 1 0\n 0.1 0.02 0.7 0\n 0.4 0.3 0.2 0\n 0.8 0.9 0 0 /')
    out = out.replace('WCONINJE\n P WATER OPEN RATE 2 1* 700 /', f'WCONINJE\n P GAS OPEN RATE {gasRate:g} 1* 700 /')
    out = out.replace('EQUIL\n 2005 300 3000 0 /', 'EQUIL\n 2005 300 3000 0 1000 0 /')   # GOC above the reservoir: oil zone
elif var == 'CP':
    # pillars at the cell corners; tops 2000 m + y tan(2 deg); one 10 m layer
    xs = [i * dx for i in range(nx + 1)]
    ys = [0.0]
    for d in dyRow: ys.append(ys[-1] + d)
    t = math.tan(math.radians(2.0))
    coord = []
    for y in ys:
        for x in xs:
            z0 = 2000 + y * t
            coord += [x, y, z0, x, y, z0 + 10]
    zc = []
    for kk in range(2):          # top then bottom of the single layer
        for j in range(31):
            for jj in range(2):
                for i in range(nx):
                    for ii in range(2):
                        zc.append(2000 + ys[j + jj] * t + 10 * kk)
    def block(vals):
        s = ''; line = []
        for v in vals:
            line.append(f'{v:.6f}'.rstrip('0').rstrip('.'))
            if len(line) == 6: s += ' ' + ' '.join(line) + '\n'; line = []
        if line: s += ' ' + ' '.join(line) + '\n'
        return s
    g0 = out.index('DX\n'); g1 = out.index('PORO\n')
    out = out[:g0] + 'COORD\n' + block(coord) + '/\nZCORN\n' + block(zc) + '/\n' + out[g1:]
    out = out.replace('WELSPECS\n P G %d 16 2005 OIL /' % ic, 'WELSPECS\n P G %d 16 1* OIL /' % ic)
    out = out.replace('EQUIL\n 2005 300 3000 0 /', 'EQUIL\n 2010 300 3000 0 /')
os.makedirs('build/h3', exist_ok=True)
name = f'build/h3/HC_{var}_{dx:g}{("_SGL" if leak else "_SG") if sub else ""}.DATA'
open(name, 'w').write(out); print(name, nx, ic, xc)
