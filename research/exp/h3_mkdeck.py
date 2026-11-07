#!/usr/bin/env python3
"""H3 helper: scale examples/CHANNELS.DATA (60x60x3) vertically to 60x60x(3*f) by repeating the layer stack -> build/h3/."""
import re, sys
f = int(sys.argv[1]) if len(sys.argv) > 1 else 5
src = open('examples/CHANNELS.DATA').read().split('\n')
out, i, n0 = [], 0, 10800
while i < len(src):
    L = src[i]
    if L.strip() == '60 60 3 /':
        out.append(f' 60 60 {3 * f} /'); i += 1; continue
    if L.strip() in ('PORO', 'PERMX', 'PERMY', 'PERMZ'):
        j = i + 1; body = []
        while src[j].strip() != '/':
            body.append(src[j]); j += 1
        out.append(L); out += body * f; out.append('/'); i = j + 1; continue
    out.append(L.replace('10800*', f'{n0 * f}*')); i += 1
open(f'build/h3/CH{3 * f}.DATA', 'w').write('\n'.join(out))
print(f'build/h3/CH{3 * f}.DATA  cells {n0 * f}')
