import numpy as np
from scipy.special import erfc
from h1_seis_opt import *
p = dict(P); p['Th_fac'] = 20.0; g = build(p); S = schedules(g, p)
r = p['rpatch'][0]; Th = g['Th']
# analytic dose within horizon for impulse at time u: f*mu*V/(4 pi k r) * erfc(r/(2 sqrt(c (Th-u)))) ~ averaged
full = p['f']*p['mu']*p['V']/(4*np.pi*p['k']*r)
qo, _ = optimise(g, p); S['OPTIMAL'] = qo
print(f"schedule-independent dose limit f*mu*V/(4 pi k r) = {full/1e6/DAY:.4f} MPa*day")
for k, q in S.items():
    dS = g['A'][0] @ q; dose = np.sum(g['w']*dS)
    print(f"  {k:16s} dose(0..20Ti)={dose/1e6/DAY:.4f} MPa*day  ({100*dose/full:.1f}% of limit)  max dS={dS.max()/1e6:.3f}")
print("Pi = dSmax/a scaling of schedule sensitivity (Th=5Ti):")
for a in [0.4e6, 0.2e6, 0.1e6, 0.05e6, 0.03e6, 0.02e6, 0.01e6]:
    pp = dict(P); pp['a'] = a; gg = build(pp); SS = schedules(gg, pp)
    E = {k: excess(q, gg) for k, q in SS.items()}; qo, _ = optimise(gg, pp); Eo = excess(qo, gg)
    Pi = (gg['A'][0] @ SS['constant']).max()/a
    print(f"  a={a/1e6:5.2f} MPa  Pi={Pi:5.2f}  opt gain vs const={100*(1-Eo/E['constant']):5.1f}%  "
          f"worst std schedule={100*(max(E.values())/E['constant']-1):6.1f}%")
