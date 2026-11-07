#!/usr/bin/env python3
"""Rasterise the official SPE11A geometry (11thSPE-CSP geometries/spe11a.geo, MIT licence, D. Glaeser) onto a
Cartesian nx x nz grid: facies of the cell centre (point in polygon over the Plane Surfaces). Output: one line per row
from the TOP row down, nx integers (1-7, 7 = impermeable / inactive). Usage: spe11_raster.py spe11a.geo nx nz out.txt"""
import re
import sys


def parse(path):
    txt = open(path).read()
    pts, lines, loops, surf_fac = {}, {}, {}, {}
    for m in re.finditer(r'^Point\((\d+)\)\s*=\s*\{([^}]*)\}', txt, re.M):
        v = [s.strip() for s in m.group(2).split(',')]
        pts[int(m.group(1))] = (float(v[0]), float(v[1]))
    for m in re.finditer(r'^Line\((\d+)\)\s*=\s*\{\s*(\d+)\s*,\s*(\d+)\s*\}', txt, re.M):
        lines[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
    for m in re.finditer(r'^Curve Loop\((\d+)\)\s*=\s*\{([^}]*)\}', txt, re.M):
        loops[int(m.group(1))] = [int(s) for s in m.group(2).split(',')]
    for m in re.finditer(r'Physical Surface\("Facies (\d)",\s*\d+\)\s*=\s*\{([^}]*)\}', txt):
        for s in m.group(2).split(','):
            surf_fac[int(s)] = int(m.group(1))
    polys = []
    for sid, fac in surf_fac.items():
        poly = []
        for e in loops[sid]:
            a, b = lines[abs(e)]
            if e < 0:
                a, b = b, a
            if not poly:
                poly.append(pts[a])
            poly.append(pts[b])
        polys.append((fac, poly))
    return polys


def inside(x, y, poly):
    c = False
    n = len(poly)
    for i in range(n - 1):
        (x1, y1), (x2, y2) = poly[i], poly[i + 1]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            c = not c
    return c


def main():
    geo, nx, nz, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    polys = parse(geo)
    xmax = max(max(p[0] for p in poly) for _, poly in polys)
    ymax = max(max(p[1] for p in poly) for _, poly in polys)
    dx, dz = xmax / nx, ymax / nz
    rows = []
    miss = 0
    for k in range(nz):
        y = ymax - (k + 0.5) * dz
        row = []
        for i in range(nx):
            x = (i + 0.5) * dx
            f = 0
            for fac, poly in polys:
                if inside(x, y, poly):
                    f = fac
                    break
            if f == 0:   # on a polygon edge: nudge
                for fac, poly in polys:
                    if inside(x + 1e-6, y + 1e-6, poly):
                        f = fac
                        break
            if f == 0:
                miss += 1
                f = row[-1] if row else 3
            row.append(f)
        rows.append(row)
    with open(out, 'w') as fh:
        fh.write('# SPE11A facies, %d x %d cells (%.4f x %.4f m), top row first; from spe11a.geo (11thSPE-CSP, MIT)\n' % (nx, nz, dx, dz))
        for r in rows:
            fh.write(' '.join(str(v) for v in r) + '\n')
    cnt = {}
    for r in rows:
        for v in r:
            cnt[v] = cnt.get(v, 0) + 1
    print('domain %.3f x %.3f m, facies counts %s, unresolved %d' % (xmax, ymax, dict(sorted(cnt.items())), miss))


if __name__ == '__main__':
    main()
