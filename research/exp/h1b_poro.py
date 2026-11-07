"""Cycle 2 (H1b): poroelastic Coulomb kernel, taper/shut-in, peak-rate objective, flowback,
window-flip robustness (Monte Carlo), two wells.  Builds on h1_seis_opt.py (same objective)."""
import numpy as np, json, sys, time
from scipy.special import erfc, logsumexp
from scipy.optimize import minimize, linprog
from scipy.integrate import cumulative_trapezoid
import h1_seis_opt as H
DAY = 86400.0

def poro_step(t, r, p, geom):
    """dCFS (Pa) at distance r for a unit step rate: pore pressure + Rudnicki-type irrotational
    poroelastic stress (exact for an infinite homogeneous medium):
      sigma_rr = -4 eta M/r^3,  sigma_tt = 2 eta (M/r^3 - p),  M = int_0^r s^2 p(s) ds   (tension +)."""
    c1, c2, eta, f = geom['c1'], geom['c2'], geom['eta'], p['f']
    P0 = p['mu']/(4*np.pi*p['k'])
    t = np.maximum(np.atleast_1d(t), 0.0); out = np.zeros_like(t)
    sg = np.linspace(0, r, 401)[1:]
    for i, ti in enumerate(t):
        if ti <= 0: continue
        pr = P0/r*erfc(r/(2*np.sqrt(p['c']*ti)))
        ps = P0/sg*erfc(sg/(2*np.sqrt(p['c']*ti)))
        M = np.trapezoid(np.concatenate([[0], sg**2*ps]), np.concatenate([[0], sg]))
        srr = -4*eta*M/r**3; stt = 2*eta*(M/r**3 - pr)
        sn = stt + (srr-stt)*c1**2; tau = (srr-stt)*c1*c2
        out[i] = tau + f*(sn + pr)
    return out

def build_poro(p, geoms):
    g = H.build(dict(p, rpatch=[p['rpatch'][0]]))
    Ti = g['Ti']; n = g['n']; tj = np.linspace(0, Ti, n+1); s = g['s']
    A = []
    for (r, geom) in geoms:
        # step response on a fine time table, then interpolate (fast)
        tt = np.concatenate([[0], np.geomspace(60.0, g['Th']+Ti, 900)])
        st = poro_step(tt, r, p, geom)
        stepf = lambda x: np.interp(np.maximum(x, 0), tt, st)
        M = np.empty((len(s), n))
        for j in range(n): M[:, j] = stepf(s-tj[j]) - stepf(s-tj[j+1])
        A.append(M)
    g['A'] = A
    return g

def optimise_general(g, p, lb, ub, Aeq_row, beq, obj='events', q0=None):
    """events: convex excess-count;  peak: min max_t dS (LP, the t<<ta limit of min peak rate)."""
    n = g['n']; sc = abs(beq)/(g['Ti'])
    if obj == 'peak':
        M = g['A'][0]; m = M.shape[0]
        # variables [q (n), z]; min z s.t. M q - z <= 0
        c = np.zeros(n+1); c[-1] = 1
        Aub = np.hstack([M/1e6, -np.ones((m, 1))]); bub = np.zeros(m)
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=np.concatenate([Aeq_row, [0]])[None], b_eq=[beq],
                      bounds=[(l, u) for l, u in zip(lb, ub)]+[(None, None)], method='highs')
        return res.x[:n]
    if q0 is None: q0 = np.clip(np.full(n, beq/Aeq_row.sum()), lb, ub)
    E0 = H.excess(q0, g)
    def fg(y):
        e, gr = H.excess(y*sc, g, True); return e/abs(E0), gr*sc/abs(E0)
    res = minimize(lambda y: fg(y)[0], q0/sc, jac=lambda y: fg(y)[1], method='SLSQP',
                   bounds=[(l/sc, u/sc) for l, u in zip(lb, ub)],
                   constraints=[{'type': 'eq', 'fun': lambda y: Aeq_row@y*sc/beq-1, 'jac': lambda y: Aeq_row*sc/beq}],
                   options=dict(maxiter=500, ftol=1e-12))
    return res.x*sc

def table(g, S, label):
    base = H.excess(S['constant'], g)
    rows = {}
    for k, q in S.items():
        st = H.rate_stats(q, g)[0]; E = H.excess(q, g)
        rows[k] = dict(E_d=E/DAY, rel=100*(E/base-1), Rpk=st['Rpk'], tpk_d=st['tpk_d'], post_d=(st['Etot_d']-st['Eco_d']),
                       dSmax=st['dSmax_MPa'])
        print(f"   {k:18s} E={E/DAY:9.4g} d ({rows[k]['rel']:+6.1f}%)  Rpk={st['Rpk']:8.3g} @{st['tpk_d']:6.2f}d  "
              f"post-shut-in={rows[k]['post_d']:8.3g} d  dSmax={st['dSmax_MPa']:.3f}")
    return rows

if __name__ == '__main__':
    out = {}
    P = dict(H.P); P['Th_fac'] = 5.0
    # fault geometries: c1 = n.nf, c2 = n.s (n = unit vector source->patch)
    geoms = {'pore-only': dict(c1=0.0, c2=0.0, eta=0.0),
             'poro-stabilising (nf || n)': dict(c1=1.0, c2=0.0, eta=0.233),
             'poro-destabilising (45deg)': dict(c1=np.sqrt(.5), c2=np.sqrt(.5), eta=0.233)}
    print("== T1  taper / shut-in with poroelastic kernel (r=300 m, a=0.05 MPa)")
    out['T1'] = {}
    for gname, geom in geoms.items():
        g = build_poro(P, [(300.0, geom)]); S = H.schedules(g, P)
        n = g['n']; dt = g['dt']; V = P['V']; qmax = P['qmax_fac']*V/g['Ti']
        S['OPTIMAL(events)'] = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V)
        S['OPTIMAL(peak LP)'] = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V, obj='peak')
        print(f"  -- {gname}")
        out['T1'][gname] = table(g, S, gname)
    print("== T2  peak-rate: exact min max R(t) multistart (convexity not guaranteed) vs LP surrogate")
    g = build_poro(P, [(300.0, geoms['poro-destabilising (45deg)'])]); n = g['n']; dt = g['dt']; V = P['V']; qmax = P['qmax_fac']*V/g['Ti']
    def peakR(q):
        return H.rate_stats(q, g)[0]['Rpk']
    rng = np.random.default_rng(3); vals = []
    from scipy.optimize import minimize as mz
    for i in range(4):
        q0 = rng.random(n)+0.2; q0 *= V/(q0.sum()*dt); q0 = np.minimum(q0, qmax*0.99)
        # smooth max via log-sum-exp of log R (temperature 50) on the eval grid
        def sm(q):
            M = g['A'][0]; dS = M@q; s = g['s']; ta = g['ta']; a = g['a']
            K = np.exp(s/ta+dS/a); I = np.concatenate([[0], np.cumsum((K[1:]+K[:-1])/2*np.diff(s))])
            lR = np.log(K) - np.log1p(I/ta); return logsumexp(50*lR)/50
        res = mz(lambda y: sm(y*V/g['Ti']), q0*g['Ti']/V, method='SLSQP', bounds=[(0, qmax*g['Ti']/V)]*n,
                 constraints=[{'type': 'eq', 'fun': lambda y: y.sum()*dt/g['Ti']-1}], options=dict(maxiter=300))
        vals.append(peakR(res.x*V/g['Ti']))
    qlp = optimise_general(g, P, np.zeros(n), np.full(n, qmax), np.full(n, dt), V, obj='peak')
    print(f"   exact-smooth multistart peak R/r: {np.round(vals, 4)}   LP-surrogate: {peakR(qlp):.4f}   constant: {peakR(np.full(n, V/g['Ti'])):.4f}")
    out['T2'] = dict(multistart=vals, lp=peakR(qlp), const=peakR(np.full(n, V/g['Ti'])))
    print("== T3  flowback allowed after injection (net volume fixed = V), window [0, 1.5 Ti] controls")
    P3 = dict(P); P3['Ti_d'] = 15.0; P3['Th_fac'] = 10/3*1.0; P3['V'] = P['V']
    g = build_poro(P3, [(300.0, geoms['pore-only'])]); n = g['n']; dt = g['dt']; qmax = 3*P['V']/(10*DAY)
    t = (np.arange(n)+0.5)*dt; inj = t < 10*DAY
    lb0 = np.zeros(n); ub0 = np.where(inj, qmax, 0.0)
    q_noflow = optimise_general(g, P3, lb0, ub0, np.full(n, dt), P['V'])
    for fb in [0.1, 0.25]:
        lb = np.where(inj, 0.0, -fb*qmax); Vinj = P['V']*(1+fb)   # inject more, flow some back -> same net
        qf = optimise_general(g, P3, lb, ub0, np.full(n, dt), P['V'])
        print(f"   flowback cap {fb*100:.0f}% qmax: E={H.excess(qf, g)/DAY:.4g} d vs no-flowback optimum {H.excess(q_noflow, g)/DAY:.4g} d "
              f"({100*(H.excess(qf, g)/H.excess(q_noflow, g)-1):+.1f}%), produced back {-(qf[qf<0].sum()*dt):.0f} m3")
    qc = np.where(inj, P['V']/(10*DAY), 0.0)
    print(f"   constant-then-shut-in E={H.excess(qc, g)/DAY:.4g} d")
    out['T3'] = dict(E_const=H.excess(qc, g)/DAY, E_opt=H.excess(q_noflow, g)/DAY)
    print("== T4  window-flip robustness: Monte Carlo over r, c, a, Ti (pore-pressure kernel)")
    rng = np.random.default_rng(11); flips = 0; cnt = 0; rec = []
    for i in range(150):
        pp = dict(H.P); pp['rpatch'] = [float(10**rng.uniform(np.log10(100), np.log10(1000)))]
        pp['c'] = float(10**rng.uniform(-2, 0)); pp['a'] = float(10**rng.uniform(np.log10(0.01e6), np.log10(0.3e6)))
        pp['Ti_d'] = float(10**rng.uniform(0, 2)); pp['V'] = 8640.0*pp['Ti_d']/10; pp['n'] = 40; pp['neval'] = 800
        for th in [1.0, 5.0]:
            pp['Th_fac'] = th
        g1 = H.build(dict(pp, Th_fac=1.0)); g5 = H.build(dict(pp, Th_fac=5.0)); S = H.schedules(g1, pp)
        d1 = H.excess(S['step-up(4)'], g1)/H.excess(S['constant'], g1)-1
        d5 = H.excess(S['step-up(4)'], g5)/H.excess(S['constant'], g5)-1
        if abs(d1) < 1e-3 and abs(d5) < 1e-3: continue
        cnt += 1; flips += (d1 < 0) and (d5 > 0); rec.append((d1, d5))
    rec = np.array(rec)
    print(f"   informative samples {cnt}: step-up better during injection in {np.mean(rec[:,0]<0)*100:.0f}%, "
          f"better incl. tail in {np.mean(rec[:,1]<0)*100:.0f}%, flip (better->worse) in {100*flips/cnt:.0f}%")
    out['T4'] = dict(n=cnt, co_better=float(np.mean(rec[:, 0] < 0)), tail_better=float(np.mean(rec[:, 1] < 0)), flip=flips/cnt)
    print("== T5  two wells (W1 at 250 m, W2 at 700 m from the fault patch; second patch 400 m from W2)")
    P5 = dict(H.P); n = P5['n']
    gA = H.build(dict(P5, rpatch=[250.0, 700.0])); gB = H.build(dict(P5, rpatch=[900.0, 400.0]))
    # patch1 sees W1 at 250, W2 at 900; patch2 sees W1 at 700, W2 at 400 -> block matrices
    A1 = np.hstack([gA['A'][0], gB['A'][0]]); A2 = np.hstack([gA['A'][1], gB['A'][1]])
    g2 = dict(gA); g2['A'] = [A1, A2]; g2['n'] = 2*n
    dt = gA['dt']; V = 2*P5['V']; qmax = P5['qmax_fac']*P5['V']/gA['Ti']
    qe = np.full(2*n, P5['V']/gA['Ti'])
    qo = optimise_general(g2, P5, np.zeros(2*n), np.full(2*n, qmax), np.full(2*n, dt), V, q0=qe)
    print(f"   equal constant split E={H.excess(qe, g2)/DAY:.4g} d; global optimum E={H.excess(qo, g2)/DAY:.4g} d "
          f"({100*(H.excess(qo, g2)/H.excess(qe, g2)-1):+.1f}%); W1 share {qo[:n].sum()/qo.sum()*100:.0f}%")
    out['T5'] = dict(E_equal=H.excess(qe, g2)/DAY, E_opt=H.excess(qo, g2)/DAY, w1_share=float(qo[:n].sum()/qo.sum()))
    json.dump(out, open('h1b_results.json', 'w'), default=float)
