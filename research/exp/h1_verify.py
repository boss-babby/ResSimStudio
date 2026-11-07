import numpy as np, json
from scipy.integrate import solve_ivp
from h1_seis_opt import *
p = dict(P); p['ta_yr'] = 0.05     # short t_a so depletion matters (stringent test)
g = build(p); S = schedules(g, p)
ta, a = g['ta'], g['a']; taur = a/ta
print("V1: identity check vs original Dieterich gamma-ODE (r=1)")
for name in ['constant', 'cyclic(6x,50%)', 'step-up(4)']:
    q = S[name]; M = g['A'][0]; s = g['s']
    dS = M @ q
    dSdt = lambda t: np.interp(t, s, np.gradient(dS, s))
    # gamma-ODE: dg/dt = (1 - g*Sdot)/a, R = 1/(g*taur); use log-gamma for stiffness
    def rhs(t, y):
        lg, N = y; gam = np.exp(lg)
        return [(1/gam - (taur + dSdt(t)))/a, 1/(gam*taur)]
    sol = solve_ivp(rhs, [0, g['Th']], [np.log(1/taur), 0.0], method='LSODA', rtol=1e-10, atol=1e-12,
                    max_step=g['Ti']/2000)
    Eode = sol.y[1, -1] - g['Th']
    Eint = excess(q, g)
    print(f"  {name:16s} E_ode={Eode/DAY:.6f} d  E_integral={Eint/DAY:.6f} d  rel={abs(Eode-Eint)/abs(Eint):.2e}")
print("V2: midpoint convexity on 2000 random pairs (3 patches)")
p2 = dict(P); p2['rpatch'] = [150., 300., 600.]; g2 = build(p2); rng = np.random.default_rng(7); viol = 0; worst = 0
for i in range(2000):
    q1 = rng.random(g2['n'])**4*0.05; q2 = rng.random(g2['n'])**4*0.05; lam = rng.random()
    lhs = excess(lam*q1+(1-lam)*q2, g2); rhs_ = lam*excess(q1, g2)+(1-lam)*excess(q2, g2)
    gap = (lhs-rhs_)/abs(rhs_); worst = max(worst, gap); viol += gap > 1e-10
print(f"  violations={viol}  worst (lhs-rhs)/rhs={worst:.2e}  (<=0 expected)")
print("V3: gradient vs finite difference")
q = S['cyclic(6x,50%)']; E, gr = excess(q, g2, True); j = 37; h = 1e-6*q.mean()
e = np.zeros_like(q); e[j] = h
print(f"  dE/dq[{j}] adj={gr[j]:.6e} fd={(excess(q+e, g2)-excess(q-e, g2))/(2*h):.6e}")
print("V4: KKT + shape of optimum (baseline, Th=5Ti and Th=1Ti)")
for th in [5.0, 1.0]:
    p3 = dict(P); p3['Th_fac'] = th; g3 = build(p3); qo, res = optimise(g3, p3)
    E, gr = excess(qo, g3, True); qmax = p3['qmax_fac']*p3['V']/g3['Ti']
    inner = (qo > 1e-4*qmax) & (qo < qmax*(1-1e-4))
    print(f"  Th={th}Ti: grad on interior support rel-spread={np.ptp(gr[inner])/np.abs(gr[inner]).mean():.2e} "
          f"n_int={inner.sum()}, n_zero={(qo<=1e-4*qmax).sum()}, n_max={(qo>=qmax*(1-1e-4)).sum()}")
    # unique fault-stress history check & coarse shape (10 bins, relative to constant)
    qb = qo.reshape(10, -1).mean(1)/(p3['V']/g3['Ti'])
    print("   10-bin q/q_const:", np.round(qb, 2))
    dS = g3['A'][0] @ qo; i = np.searchsorted(g3['s'], [2, 4, 6, 8, 10, 12, 20, 50]*np.array(DAY))
    print("   dS(MPa) at 2,4,6,8,10,12,20,50 d:", np.round(dS[np.minimum(i, len(dS)-1)]/1e6, 4))
