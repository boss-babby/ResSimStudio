#!/usr/bin/env python3
"""H2 helper: HFGROW_INJECTION example at along-fracture cell size dx (10 / 5 / 2.5 m) -> build/h3/HF_<dx>[_SG].DATA.
nx is odd so the injection cell is centred; the fracture trace is centre +- 100 m; everything else as the example."""
import sys
dx = float(sys.argv[1]); sub = len(sys.argv) > 2 and sys.argv[2] in ('SG', 'SGL', 'SGLC'); leak = len(sys.argv) > 2 and sys.argv[2] in ('SGL', 'SGLC')
src = open('examples/HFGROW_INJECTION.DATA').read()
nx = int(round(810 / dx)); nx += (nx + 1) % 2
ic = nx // 2 + 1; xc = (ic - 0.5) * dx; n = nx * 31
dy = src.split('DY\n')[1].split('/')[0].replace('81*', f'{nx}*')
out = src.replace(' 81 31 1 /', f' {nx} 31 1 /').replace('2511*10 /\nDY', f'{n}*{dx:g} /\nDY', 1)
out = out.split('DY\n')[0] + 'DY\n' + dy + '/' + out.split('DY\n')[1].split('/', 1)[1]
out = out.replace('2511*', f'{n}*')
out = out.replace('HFRAC F1 305 688.793902 505 688.793902', f'HFRAC F1 {xc - 100:g} 688.793902 {xc + 100:g} 688.793902')
out = out.replace('HFGROW F1 1.0 340 0 20', 'HFGROW F1 1.0 340 0 20 1' + (' SUBGRID' if sub else '') + (' LEAK' if leak else '') + ('C' if len(sys.argv) > 2 and sys.argv[2] == 'SGLC' else ''))
out = out.replace(' P G 41 16 ', f' P G {ic} 16 ').replace(' P 41 16 1 1 ', f' P {ic} 16 1 1 ')
name = f'build/h3/HF_{dx:g}{("_SGL" if leak else "_SG") if sub else ""}.DATA'
open(name, 'w').write(out); print(name, nx, ic, xc)
