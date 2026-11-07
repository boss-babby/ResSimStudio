#!/usr/bin/env python3
"""C47 developer tool: time-aligned difference of two ECLIPSE/OPM summaries (e.g. two OPM versions of a reference).
usage: smrydiff.py A.SMSPEC B.SMSPEC [N=25] [KW,KW..]   -> vectors ranked by max |A - B(t_A)| / max|A| (B linearly interpolated)"""
import struct, sys, bisect
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
def load(sp):
    kws, wg, nums = [], [], []
    for kw, t, c, d in blocks(sp):
        if kw == 'KEYWORDS': kws = [d[j*8:j*8+8].decode().strip() for j in range(c)]
        if kw in ('WGNAMES', 'NAMES'): wg = [d[j*8:j*8+8].decode().strip() for j in range(c)]
        if kw == 'NUMS': nums = list(struct.unpack('>%di' % c, d))
    rows = [struct.unpack('>%df' % c, d) for kw, t, c, d in blocks(sp.replace('SMSPEC', 'UNSMRY')) if kw == 'PARAMS']
    keys = ['%s:%s' % (k, w) if w and w != ':+:+:+:+' else ('%s:%d' % (k, n) if k.startswith(('B', 'C', 'R')) else k) for k, w, n in zip(kws, wg, nums or [0] * len(kws))]
    return keys, rows
if __name__ == '__main__':
    ka, ra = load(sys.argv[1]); kb, rb = load(sys.argv[2])
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    only = sys.argv[4].split(',') if len(sys.argv) > 4 else None
    ta = [r[ka.index('TIME')] for r in ra]; tb = [r[kb.index('TIME')] for r in rb]
    def at(col, t):
        i = bisect.bisect_left(tb, t)
        if i <= 0: return rb[0][col]
        if i >= len(tb): return rb[-1][col]
        w = (t - tb[i-1]) / (tb[i] - tb[i-1]); return (1 - w) * rb[i-1][col] + w * rb[i][col]
    out = []
    for j, k in enumerate(ka):
        if k not in kb or k == 'TIME' or (only and k.split(':')[0] not in only): continue
        jb = kb.index(k)
        a = [r[j] for r in ra]; mx = max(abs(x) for x in a)
        if mx == 0: continue
        dev = [(abs(x - at(jb, t)), t) for x, t in zip(a, ta) if t <= tb[-1]]
        d, t = max(dev)
        out.append((d / mx, k, t))
    out.sort(reverse=True)
    print('steps A %d (t_end %.0f), B %d (t_end %.0f)' % (len(ra), ta[-1], len(rb), tb[-1]))
    for d, k, t in out[:N]: print('%7.2f %%  %-22s at t = %.0f' % (100 * d, k, t))
