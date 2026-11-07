#!/usr/bin/env python3
"""Print summary vectors of an ECLIPSE/OPM binary summary (SMSPEC + UNSMRY).
usage: smry.py CASE.SMSPEC KW[:WGNAME][:NUM] ...   (e.g. SOFR:PROD01:3 WBHP:PROD01 FGOR)"""
import struct, sys
def blocks(path):
    d = open(path, 'rb').read(); i = 0
    while i < len(d):
        n = struct.unpack('>i', d[i:i+4])[0]; hdr = d[i+4:i+4+n]; i += 8 + n
        kw = hdr[:8].decode().strip(); cnt = struct.unpack('>i', hdr[8:12])[0]; typ = hdr[12:16].decode()
        sz = {'INTE': 4, 'REAL': 4, 'DOUB': 8, 'LOGI': 4, 'CHAR': 8, 'MESS': 0}.get(typ, 8)
        data = b''; rem = cnt * sz
        while rem > 0:
            m = struct.unpack('>i', d[i:i+4])[0]; data += d[i+4:i+4+m]; i += 8 + m; rem -= m
        yield kw, typ, cnt, data
sp = sys.argv[1]
b = {}
for k, t, c, d in blocks(sp):
    b.setdefault(k, (t, c, d))
n = b['KEYWORDS'][1]
kws = [b['KEYWORDS'][2][i*8:i*8+8].decode().strip() for i in range(n)]
wg = [b['WGNAMES' if 'WGNAMES' in b else 'NAMES'][2][i*8:i*8+8].decode().strip() for i in range(n)]
nums = struct.unpack('>%di' % n, b['NUMS'][2]) if 'NUMS' in b else [0] * n
def idx(spec):
    p = spec.split(':')
    for i in range(n):
        if kws[i] != p[0]: continue
        if len(p) > 1 and p[1] and wg[i] != p[1]: continue
        if len(p) > 2 and nums[i] != int(p[2]): continue
        return i
    raise SystemExit('not found: ' + spec)
cols = [idx(s) for s in sys.argv[2:]]
ti = kws.index('TIME')
print('%8s ' % 'TIME' + ' '.join('%14s' % s for s in sys.argv[2:]))
for k, t, c, d in blocks(sp.replace('.SMSPEC', '.UNSMRY')):
    if k == 'PARAMS':
        v = struct.unpack('>%df' % c, d)
        print('%8.1f ' % v[ti] + ' '.join('%14.5g' % v[j] for j in cols))
