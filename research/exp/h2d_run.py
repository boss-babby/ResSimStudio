#!/usr/bin/env python3
"""H2d (Research Cycle 9, pre-registered in CYCLE9_PLAN.md): H2c GAS at 300 sm3/d (growth over weeks, not hours).
Pre-check: SG t60 at dx 10 m < 10 d -> rate 100 sm3/d once. Criterion (unchanged from H2c): SGL t60 spread < 1 d over
dx 10 / 5 / 2.5 m and smaller than the SG spread -> SUPPORTED.   usage (from ressim/): h2d_run.py DRIVER > research/exp/h2d.log"""
import subprocess, sys, os
drv = sys.argv[1]
os.makedirs('build/h3', exist_ok=True)

def one(dx, mode, rate):
    name = subprocess.run([sys.executable, 'research/exp/h2c_mkdeck.py', str(dx), mode, 'GAS', str(rate)], capture_output=True, text=True, check=True).stdout.split()[0]
    out = subprocess.run([drv, name], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS='1'), timeout=7200).stdout
    ev = [tuple(map(float, l.split()[1:])) for l in out.splitlines() if l.startswith('EV ')]
    def tL(L):
        ts = [e[0] for e in ev if e[2] >= L]
        return min(ts) if ts else float('nan')
    done = [l for l in out.splitlines() if l.startswith('DONE')]
    t60, t90 = tL(60), tL(90)
    print('h2d_GAS%g_%g_%s t60 %.2f t90 %.2f | %s' % (rate, dx, mode, t60, t90, done[0] if done else 'NO DONE'), flush=True)
    return t60

rate = float(sys.argv[2]) if len(sys.argv) > 2 else 300     # H2d' (exploratory addendum): 600
pre = one(10, 'SG', rate)
if len(sys.argv) <= 2 and pre < 10:          # (pre-registered fallback only for growth that is too fast; nan = no growth)
    print('PRECHECK SG t60 at 10 m = %.2f d < 10 d -> rate 100 sm3/d (pre-registered fallback)' % pre, flush=True)
    rate = 100
    pre = one(10, 'SG', rate)
res = {('SG', 10): pre}
for dx in (10, 5, 2.5):
    for mode in ('SG', 'SGL'):
        if (mode, dx) in res:
            continue
        res[(mode, dx)] = one(dx, mode, rate)
sp = {m: max(res[(m, d)] for d in (10, 5, 2.5)) - min(res[(m, d)] for d in (10, 5, 2.5)) for m in ('SG', 'SGL')}
import math
v = ('INCONCLUSIVE (no growth to 60 m)' if any(math.isnan(x) for x in res.values()) else
     'SUPPORTED' if sp['SGL'] < 1.0 and sp['SGL'] < sp['SG'] else ('NON-DISCRIMINATING' if res[('SG', 10)] < 10 else 'REFUTED'))
print('SUMMARY rate %g: t60 spread SGL %.2f d, SG %.2f d; SGL %s; SG %s' % (rate, sp['SGL'], sp['SG'],
      ' / '.join('%.2f' % res[('SGL', d)] for d in (10, 5, 2.5)), ' / '.join('%.2f' % res[('SG', d)] for d in (10, 5, 2.5))), flush=True)
print('VERDICT H2d%s %s' % ("'" if len(sys.argv) > 2 else '', v), flush=True)
print('H2DDONE', flush=True)
