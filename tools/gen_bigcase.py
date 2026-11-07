#!/usr/bin/env python3
"""C37 performance round: synthetic black-oil benchmark deck of any size (SPE1 PVT / rel-perm, heterogeneous layered
permeability, 4 corner producers + 1 central water injector, 3-phase with gas liberation below the bubble point).

    python3 tools/gen_bigcase.py NX NY NZ out.DATA [--months M]
"""
import math, sys, os

def main():
    if len(sys.argv) < 5:
        print(__doc__); sys.exit(1)
    nx, ny, nz = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    out = sys.argv[4]
    months = int(sys.argv[sys.argv.index('--months') + 1]) if '--months' in sys.argv else 6
    here = os.path.dirname(os.path.abspath(__file__))
    spe1 = open(os.path.join(here, '..', 'examples', 'SPE1.DATA')).read()
    props = spe1[spe1.index('\nPROPS') + 1: spe1.index('\nSOLUTION') + 1]
    n = nx * ny * nz
    def arr(f):
        vals = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    vals.append(f(i, j, k))
        lines, row = [], []
        for v in vals:
            row.append('%.4g' % v)
            if len(row) == 12:
                lines.append(' '.join(row)); row = []
        if row: lines.append(' '.join(row))
        return '\n'.join(lines)
    def lperm(i, j, k):
        g = 0.9 * math.sin(i / 7.3 + 1.7 * k) * math.cos(j / 9.1 - 0.6 * k) + 0.5 * math.sin((i + 2 * j) / 23.0) + 0.35 * math.cos(1.3 * k)
        return 10 ** (2.0 + g)
    d = []
    d.append('-- synthetic benchmark %dx%dx%d = %d cells (tools/gen_bigcase.py)\nRUNSPEC\nTITLE\nSynthetic %d-cell black-oil benchmark\n' % (nx, ny, nz, n, n))
    d.append('DIMENS\n %d %d %d /\nOIL\nWATER\nGAS\nDISGAS\nFIELD\nTABDIMS\n 1 1 20 20 /\nEQLDIMS\n 1 1* 20 /\nWELLDIMS\n 5 %d 1 5 /\nSTART\n 1 \'JAN\' 2015 /\nUNIFOUT\n' % (nx, ny, nz, nz))
    d.append('GRID\nDX\n %d*100 /\nDY\n %d*100 /\nDZ\n %d*10 /\nTOPS\n %d*8325 /\n' % (n, n, n, nx * ny))
    d.append('PERMX\n' + arr(lperm) + '\n/\nCOPY\n PERMX PERMY /\n/\nPERMZ\n' + arr(lambda i, j, k: 0.1 * lperm(i, j, k)) + '\n/\n')
    d.append('PORO\n' + arr(lambda i, j, k: 0.12 + 0.05 * math.log10(lperm(i, j, k)) / 2.0) + '\n/\n')
    d.append(props)
    d.append('SOLUTION\nEQUIL\n 8400 4800 8600 0 8300 0 1 1* 0 /\nRSVD\n 7300 1.27\n 9450 1.27\n/\n')
    d.append('SUMMARY\nFOPR\nFWPR\nFGPR\nFWIR\nFPR\nFOPT\n')
    wells = [('P1', 1, 1), ('P2', nx, 1), ('P3', 1, ny), ('P4', nx, ny)]
    s = 'SCHEDULE\nWELSPECS\n'
    for w, i, j in wells: s += " '%s' 'G' %d %d 1* OIL /\n" % (w, i, j)
    s += " 'INJ' 'G' %d %d 1* WATER /\n/\nCOMPDAT\n" % ((nx + 1) // 2, (ny + 1) // 2)
    for w, i, j in wells: s += " '%s' %d %d 1 %d OPEN 1* 1* 0.5 /\n" % (w, i, j, nz)
    s += " 'INJ' %d %d 1 %d OPEN 1* 1* 0.5 /\n/\nWCONPROD\n" % ((nx + 1) // 2, (ny + 1) // 2, nz)
    rate = max(500.0, n * 0.1)                 # total liquid rate (stb/d): ~6 % of the pore volume per half year
    for w, _, _ in wells: s += " '%s' OPEN LRAT 3* %.0f 1* 2500 /\n" % (w, rate / 4)
    s += "/\nWCONINJE\n 'INJ' WATER OPEN RATE %.0f 1* 7000 /\n/\nTSTEP\n %d*30.4375 /\nEND\n" % (1.1 * rate, months)
    d.append(s)
    open(out, 'w').write(''.join(d))
    print('wrote %s: %d cells' % (out, n))

main()
