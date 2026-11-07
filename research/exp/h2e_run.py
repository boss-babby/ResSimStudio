#!/usr/bin/env python3
"""H2e (Research Cycle 10, pre-registered in CYCLE10_PLAN.md): gas-aware sub-cell storage (HFGROW ... SUBGRID LEAKC).
H2d' GAS deck at 600 sm3/d, dx 10 / 5 / 2.5 m, SGLC vs the F30 SGL / SG values; water guard (F21 deck, SGLC vs SGL).
usage (from ressim/): h2e_run.py DRIVER > research/exp/h2e.log"""
import subprocess, sys, os, json
drv = sys.argv[1]
os.makedirs('build/h3', exist_ok=True)
F30 = {'SG': {10: (18.80, 40.78), 5: (17.80, 39.32), 2.5: (17.55, 38.56)},
       'SGL': {10: (16.60, 36.21), 5: (16.60, 37.30), 2.5: (17.20, 38.40)}}

def times(out):
    ev = [tuple(map(float, l.split()[1:])) for l in out.splitlines() if l.startswith('EV ')]
    tL = lambda L: min([e[0] for e in ev if e[2] >= L] or [float('nan')])
    done = [l for l in out.splitlines() if l.startswith('DONE')]
    return tL(60), tL(90), done[0] if done else 'NO DONE'

def run(gen, dx, mode, extra=()):
    name = subprocess.run([sys.executable, gen, str(dx), mode] + list(extra), capture_output=True, text=True, check=True).stdout.split()[0]
    out = subprocess.run([drv, name], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS='1'), timeout=7200).stdout
    return times(out)

r = run('research/exp/h2c_mkdeck.py', 10, 'SGL', ['GAS', '600'])   # reproduction check of F30 with the current engine
print('h2e_GAS600_10_SGL (F30 repro: 16.60 / 36.21) t60 %.2f t90 %.2f | %s' % r, flush=True)
res = {}
for dx in (10, 5, 2.5):
    t60, t90, d = run('research/exp/h2c_mkdeck.py', dx, 'SGLC', ['GAS', '600'])
    res[dx] = (t60, t90)
    print('h2e_GAS600_%g_SGLC t60 %.2f t90 %.2f | %s' % (dx, t60, t90, d), flush=True)
sp = lambda i, tab: max(tab[d][i] for d in (10, 5, 2.5)) - min(tab[d][i] for d in (10, 5, 2.5))
s60, s90 = sp(0, res), sp(1, res)
g = {}
for mode in ('SGL', 'SGLC'):
    g[mode] = run('research/exp/h2_mkdeck.py', 10, mode)
    print('h2e_WATER_10_%s t60 %.2f t90 %.2f | %s' % ((mode,) + g[mode]), flush=True)
guard = abs(g['SGLC'][0] - g['SGL'][0]) <= 0.3 and abs(g['SGLC'][1] - g['SGL'][1]) <= 0.3
if s60 < 1 and s90 < 1 and s90 < sp(1, F30['SG']):
    v = 'SUPPORTED'
elif s90 >= sp(1, F30['SGL']):
    v = 'REFUTED'
else:
    v = 'INCONCLUSIVE'
print('SUMMARY SGLC t60 %s (spread %.2f) t90 %s (spread %.2f) | SGL t90 spread %.2f, SG t90 spread %.2f | water guard %s' % (
    ' / '.join('%.2f' % res[d][0] for d in (10, 5, 2.5)), s60, ' / '.join('%.2f' % res[d][1] for d in (10, 5, 2.5)), s90,
    sp(1, F30['SGL']), sp(1, F30['SG']), 'OK' if guard else 'FAIL'), flush=True)
print('VERDICT H2e %s%s' % (v, '' if guard else ' (guard failed: result void)'), flush=True)
json.dump({'SGLC': {str(k): v_ for k, v_ in res.items()}, 'water': g, 'spread60': s60, 'spread90': s90, 'guard': guard, 'verdict': v},
          open('research/exp/h2e_results.json', 'w'), indent=1)
