"""H1: optimal injection schedule under Dieterich(1994) seismicity + pore-pressure diffusion.
Key identity (Heimisson & Segall 2018 / Heimisson 2019 ODE form):
  N(T)/r = t_a ln(1 + (1/t_a) int_0^T exp(S(s)/a) ds),   S = tau_r s + dS(s),  a = A sigma0, t_a = a/tau_r
=> excess E(T) = N(T)/r - T = t_a ln( e^{-T/t_a} + (1/t_a) int e^{(s-T)/t_a} e^{dS/a} ds )
dS linear in rates q  =>  ln(...) is log-sum-exp of affine maps  =>  E convex in q (any #wells, #patches).
"""
import numpy as np, json, sys, time
from scipy.special import erfc, logsumexp
from scipy.optimize import minimize
DAY = 86400.0
P = dict(k=1e-14, mu=1e-3, c=0.1, f=0.6, a=0.05e6, ta_yr=100.0, Ti_d=10.0, V=8640.0,
         rpatch=[300.0], qmax_fac=3.0, n=100, Th_fac=5.0, neval=3000)

def pstep(t, r, p):            # pressure (Pa) at distance r for unit step rate (1 m3/s), 3D point source
    t = np.maximum(t, 0.0)
    out = np.zeros_like(t); m = t > 0
    out[m] = p['mu']/(4*np.pi*p['k']*r)*erfc(r/(2*np.sqrt(p['c']*t[m])))
    return out

def build(p):
    Ti = p['Ti_d']*DAY; Th = p['Th_fac']*Ti; n = p['n']
    tj = np.linspace(0, Ti, n+1); dt = Ti/n
    # eval grid: dense where things change (geometric after Ti)
    s1 = np.linspace(0, Ti, p['neval']//2+1)
    s2 = Ti + (Th-Ti)*np.linspace(0, 1, p['neval']//2+1)[1:]**2 if Th > Ti else np.array([])
    s = np.concatenate([s1, s2])
    w = np.zeros_like(s); ds = np.diff(s); w[:-1] += ds/2; w[1:] += ds/2   # trapezoid
    A = []
    for r in p['rpatch']:
        M = np.empty((len(s), n))
        for j in range(n):
            M[:, j] = p['f']*(pstep(s-tj[j], r, p) - pstep(s-tj[j+1], r, p))
        A.append(M)
    return dict(Ti=Ti, Th=Th, dt=dt, s=s, w=w, A=A, ta=p['ta_yr']*365.25*DAY, a=p['a'], n=n)

def excess(q, g, grad=False):
    """sum over patches of E (seconds of background-equivalent events); optional gradient."""
    ta, a, Th, s, w = g['ta'], g['a'], g['Th'], g['s'], g['w']
    tot, gr = 0.0, np.zeros_like(q)
    logc = -Th/ta
    with np.errstate(divide='ignore'):
        logw = np.log(w/ta) + (s-Th)/ta
    for M in g['A']:
        dS = M @ q
        terms = np.concatenate([[logc], logw + dS/a])
        L = logsumexp(terms)
        tot += ta*L
        if grad:
            sm = np.exp(terms[1:]-L)          # softmax weights (patch-wise)
            gr += ta*(M.T @ sm)/a
    return (tot, gr) if grad else tot

def rate_stats(q, g):
    """peak normalised rate R/r and split co/post-injection excess, max fault dS (MPa)."""
    ta, a, s, w = g['ta'], g['a'], g['s'], g['w']
    out = []
    for M in g['A']:
        dS = M @ q
        K = np.exp(s/ta + dS/a)
        I = np.concatenate([[0], np.cumsum((K[1:]+K[:-1])/2*np.diff(s))])
        R = K/(1+I/ta)
        N = ta*np.log1p(I/ta) - s
        iTi = np.searchsorted(s, g['Ti'])
        out.append(dict(Rpk=float(R.max()), tpk_d=float(s[R.argmax()]/DAY), Eco_d=float(N[iTi]/DAY),
                        Etot_d=float(N[-1]/DAY), dSmax_MPa=float(dS.max()/1e6)))
    return out

def schedules(g, p):
    n, dt, V, Ti = g['n'], g['dt'], p['V'], g['Ti']
    q0 = V/Ti; x = (np.arange(n)+0.5)/n
    S = {}
    S['constant'] = np.full(n, q0)
    S['step-up(4)'] = np.floor(x*4)+1.0
    S['step-down(4)'] = 4.0-np.floor(x*4)
    S['cyclic(6x,50%)'] = (np.floor(x*12) % 2 == 0)*1.0
    S['taper-last30%'] = np.where(x < 0.7, 1.0, np.maximum(0, (1-x)/0.3))
    S['ramp-up-linear'] = x
    for k in S: S[k] = S[k]*V/(S[k].sum()*dt)
    return S

def optimise(g, p, q_init=None, maxiter=400):
    n, dt, V, Ti = g['n'], g['dt'], p['V'], g['Ti']
    qmax = p['qmax_fac']*V/Ti
    sc = V/Ti
    if q_init is None: q_init = np.full(n, sc)
    E0 = excess(q_init, g)
    fun = lambda y: tuple(np.array(v)/abs(E0) for v in ((lambda e: (e[0], e[1]*sc))(excess(y*sc, g, True))))
    cons = [{'type': 'eq', 'fun': lambda y: y.sum()*dt*sc/V - 1.0, 'jac': lambda y: np.full(n, dt*sc/V)}]
    res = minimize(lambda y: fun(y)[0], q_init/sc, jac=lambda y: fun(y)[1], method='SLSQP',
                   bounds=[(0, qmax/sc)]*n, constraints=cons, options=dict(maxiter=maxiter, ftol=1e-12))
    return res.x*sc, res

def run_case(p, tag, nstarts=0, verbose=True):
    g = build(p); S = schedules(g, p)
    rows = {}
    for k, q in S.items():
        rows[k] = dict(E_d=excess(q, g)/DAY, stats=rate_stats(q, g))
    t0 = time.time(); qo, res = optimise(g, p); topt = time.time()-t0
    rows['OPTIMAL'] = dict(E_d=excess(qo, g)/DAY, stats=rate_stats(qo, g), success=bool(res.success),
                           nit=int(res.nit), sec=topt)
    starts = []
    rng = np.random.default_rng(1)
    for i in range(nstarts):
        qi = rng.random(g['n'])**3; qi *= p['V']/(qi.sum()*g['dt'])
        qi = np.minimum(qi, p['qmax_fac']*p['V']/g['Ti']*0.999); qi *= p['V']/(qi.sum()*g['dt'])
        qr, rr = optimise(g, p, qi)
        starts.append(dict(E_d=excess(qr, g)/DAY, dq_rel=float(np.abs(qr-qo).max()/qo.max())))
    base = rows['constant']['E_d']
    if verbose:
        print(f"== {tag}")
        for k, v in rows.items():
            st = v['stats'][0]
            print(f"  {k:16s} E={v['E_d']:10.4g} d  ({100*(v['E_d']/base-1):+6.1f}% vs const)  "
                  f"Rpk={st['Rpk']:9.3g} @{st['tpk_d']:6.2f}d  co={st['Eco_d']:9.3g} dSmax={st['dSmax_MPa']:.3f}")
        if starts: print("  multistart:", [f"{s['E_d']:.6g}/{s['dq_rel']:.1e}" for s in starts])
    return dict(tag=tag, params={k: v for k, v in p.items()}, rows=rows, q_opt=qo.tolist(),
                starts=starts, t_d=((np.arange(g['n'])+0.5)*g['dt']/DAY).tolist())

if __name__ == '__main__':
    out = []
    out.append(run_case(dict(P), 'baseline Th=5Ti', nstarts=4))
    for th in [1.0, 2.0, 20.0]:
        q = dict(P); q['Th_fac'] = th; out.append(run_case(q, f'Th={th:g}Ti'))
    for a in [0.02e6, 0.2e6]:
        q = dict(P); q['a'] = a; out.append(run_case(q, f'a={a/1e6:g}MPa'))
    q = dict(P); q['rpatch'] = [150.0, 300.0, 600.0]; out.append(run_case(q, 'patches 150/300/600 m', nstarts=2))
    q = dict(P); q['ta_yr'] = 0.05; out.append(run_case(q, 'ta=0.05yr (strong depletion)'))
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else 'h1_results.json', 'w'))
