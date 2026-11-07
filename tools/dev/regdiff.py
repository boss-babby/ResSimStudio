#!/usr/bin/env python3
"""Compare two ResSim deck-regression results (HTML report or console log) case by case.
usage: regdiff.py BASE NEW   (each: report .html or --regress console log)"""
import re, sys, html

def load(path):
    s = open(path, encoding='utf-8', errors='replace').read()
    out = {}
    if '<html' in s[:2000].lower() or '<table' in s:
        t = html.unescape(re.sub(r'<[^>]+>', ' ', s))
        for m in re.finditer(r'(\S+)\s+(PASS|FAIL|ERROR)\s+(\S+)\s+([\d.]+)\s+\d+\s+\d+\s+[\d.]+\s+(?:[\d.]+\s+)?(?:[^m]*?material worst (\S+) ([\d.]+) %)?', t):
            if m.group(1) in out: continue
            out[m.group(1)] = (m.group(2), m.group(3), float(m.group(4)), m.group(5), float(m.group(6)) if m.group(6) else None)
    else:
        rx = re.compile(r'(PASS|FAIL|ERROR)\s+(\S+)\s+worst (\S+) ([\d.]+) %.*?(?:\[material worst (\S+) ([\d.]+) %\])?$')
        for line in s.splitlines():
            m = rx.search(line.strip())
            if m: out[m.group(2)] = (m.group(1), m.group(3), float(m.group(4)), m.group(5), float(m.group(6)) if m.group(6) else None)
    return out

a, b = load(sys.argv[1]), load(sys.argv[2])
print(f"{'case':48s} {'base':>14s} {'new':>14s}  change")
for k in sorted(set(a) | set(b)):
    x, y = a.get(k), b.get(k)
    f = lambda r: '-' if r is None else f"{r[0][0]} {r[2]:.3g}" + (f"[{r[4]:.3g}]" if r[4] is not None else '')
    d = '' if (x is None or y is None) else f"{y[2]-x[2]:+.2f}"
    print(f"{k:48s} {f(x):>14s} {f(y):>14s}  {d}")
