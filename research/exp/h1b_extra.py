import numpy as np, json
import h1_seis_opt as H
from h1b_poro import build_poro, optimise_general, table
from scipy.optimize import minimize
DAY = 86400.0; P = dict(H.P); out = {}
geom = dict(c1=np.sqrt(.5), c2=-np.sqrt(.5), eta=0.233)
print("== T1b truly destabilising poroelastic geometry (45 deg, slip sense reversed)")
g = build_poro(P, [(300.0, geom)]); S = H.schedules(g, P); n = g['n']; dt = g['dt']; V = P['V']; qmax = 3*V/g['Ti']
S['OPTIMAL(events)'] = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V)
S['OPTIMAL(peak LP)'] = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V, obj='peak')
out['T1b'] = table(g, S, 'destab')
print("== T6 Pareto front events vs peak dS (pore-only, r=300): min E s.t. max dS <= z")
g = H.build(P); M = g['A'][0]
qE = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V)
qP = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V, obj='peak')
zmin, zmax = (M@qP).max(), (M@qE).max(); front = []
for z in np.linspace(zmin*1.0005, zmax, 6):
    sc = V/g['Ti']; E0 = H.excess(qE, g)
    cons = [{'type': 'eq', 'fun': lambda y: y.sum()*dt*sc/V-1, 'jac': lambda y: np.full(n, dt*sc/V)},
            {'type': 'ineq', 'fun': lambda y, z=z: (z - M@(y*sc))/1e4, 'jac': lambda y: -M*sc/1e4}]
    r = minimize(lambda y: H.excess(y*sc, g)/E0, qP/sc, jac=lambda y: H.excess(y*sc, g, True)[1]*sc/E0,
                 method='SLSQP', bounds=[(0, 3)]*n, constraints=cons, options=dict(maxiter=500, ftol=1e-12))
    q = r.x*sc; st = H.rate_stats(q, g)[0]
    front.append((float(z/1e6), H.excess(q, g)/DAY, st['Rpk'])); print(f"   max dS<={z/1e6:.4f} MPa -> E={front[-1][1]:.4g} d, Rpk={st['Rpk']:.3g}")
out['T6'] = front
print("== T3b EGS stimulation: injection fixed (constant 10 d), optimal flowback after shut-in with budget B")
P3 = dict(P); P3['Ti_d'] = 20.0; P3['Th_fac'] = 2.5; P3['V'] = V
g = H.build(P3); n3 = g['n']; dt3 = g['dt']; t = (np.arange(n3)+0.5)*dt3; inj = t < 10*DAY
q_inj = np.where(inj, V/(10*DAY), 0.0); E_ref = H.excess(q_inj, g); res = []
for B in [0.05, 0.1, 0.2, 0.4]:
    sc = V/(10*DAY); fb_max = 0.5*sc
    def full(y): q = q_inj.copy(); q[~inj] = -y; return q
    m = (~inj).sum(); y0 = np.full(m, B*V/(m*dt3))
    def f(y):
        e, gr = H.excess(full(y), g, True); return e/E_ref, -gr[~inj]/E_ref
    r = minimize(lambda y: f(y)[0], y0, jac=lambda y: f(y)[1], method='SLSQP', bounds=[(0, fb_max)]*m,
                 constraints=[{'type': 'eq', 'fun': lambda y: y.sum()*dt3/(B*V)-1}], options=dict(maxiter=400))
    E = H.excess(full(r.x), g); yc = np.zeros(m); k = int(np.ceil(B*V/(fb_max*dt3)))
    yc[:k] = fb_max; yc[k-1] -= (yc.sum()*dt3-B*V)/dt3; Eimm = H.excess(full(yc), g)
    first = r.x[:max(1, m//5)].sum()/r.x.sum()
    print(f"   flowback {B*100:4.0f}% of V: optimal E={100*(E/E_ref-1):+6.1f}%  immediate-max-rate E={100*(Eimm/E_ref-1):+6.1f}%  "
          f"(optimal puts {first*100:.0f}% of flowback in first 20% of post window)")
    res.append(dict(B=B, opt=E/E_ref-1, immediate=Eimm/E_ref-1, front_frac=float(first)))
out['T3b'] = res
json.dump(out, open('h1b_extra.json', 'w'), default=float)
