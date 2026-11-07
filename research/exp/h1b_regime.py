import numpy as np, json
import h1_seis_opt as H
from h1b_poro import build_poro, optimise_general, poro_step
DAY = 86400.0; P = dict(H.P); n = P['n']; out = {}
print("T7: strongly coupled regime (alpha=1, nu=0.15 -> eta=0.41), does the optimum taper?")
for name, geom in [('stabilising nf||n', dict(c1=1.0, c2=0.0, eta=0.41)),
                   ('stab. 60deg', dict(c1=np.cos(np.pi/6), c2=np.sin(np.pi/6), eta=0.41)),
                   ('pressure-dominated', dict(c1=np.sqrt(.5), c2=-np.sqrt(.5), eta=0.41))]:
    g = build_poro(P, [(300.0, geom)]); S = H.schedules(g, P); dt = g['dt']; V = P['V']; qmax = 3*V/g['Ti']
    st = poro_step(np.array([1, 3, 10, 30])*DAY, 300.0, P, geom)
    S['OPT events'] = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V)
    S['OPT peak'] = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V, obj='peak')
    b = H.excess(S['constant'], g); pk = H.rate_stats(S['constant'], g)[0]
    print(f" -- {name}: step kernel dCFS at 1,3,10,30 d = {np.round(st/1e3, 2)} kPa per m3/s")
    res = {}
    for k in ['constant', 'taper-last30%', 'step-down(4)', 'OPT events', 'OPT peak']:
        s = H.rate_stats(S[k], g)[0]; E = H.excess(S[k], g)
        last = S[k][-n//5:].mean()/S[k].mean()
        res[k] = dict(E=E/b-1, Rpk=s['Rpk']/pk['Rpk']-1, tpk=s['tpk_d'], last20=last)
        print(f"    {k:14s} events {100*(E/b-1):+6.1f}%  peak {100*(s['Rpk']/pk['Rpk']-1):+6.1f}% @ {s['tpk_d']:5.2f} d   "
              f"rate in last 20% = {last:.2f}x mean")
    out[name] = res
json.dump(out, open('h1b_regime.json', 'w'), default=float)
