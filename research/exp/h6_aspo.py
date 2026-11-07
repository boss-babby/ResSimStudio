#!/usr/bin/env python3
"""Research Cycle 5 - H6a / H6b / H1h on the Aspo HRL 2015 hydraulic-fracturing data (pre-registered in
CYCLE5_PLAN.md). numpy only. Writes h6_results.json and h6_fig.png (if matplotlib is present).
Data: research/data/aspo (1 s means of GFZ doi 10.5880/GFZ.2.6.2023.004 + Kwiatek 2018 catalogue)."""
import json, math, os, glob, datetime as dt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data', 'aspo')
TA = 1e8           # Dieterich aftershock duration (s) - irrelevant for t << ta
QMIN, GAP, DMIN, POST, PMIN = 0.3, 5, 10, 600, 1.0

def load_test(name):
    parts = []
    for f in sorted(glob.glob(os.path.join(DATA, name + '_1s_*.csv'))):
        parts.append(np.genfromtxt(f, delimiter=';', skip_header=1))
    return parts

def stages(d):
    t, q = d[:, 0], d[:, 1]
    on = q > QMIN
    st = []
    i = 0
    n = len(t)
    while i < n:
        if not on[i]: i += 1; continue
        j = i
        while j + 1 < n and (on[j + 1] or (t[min(n - 1, j + GAP)] - t[j] <= GAP and on[j + 1:j + 1 + GAP].any())):
            j += 1
        # deviation from the pre-registration (documented in CYCLE5_REPORT): flow with interval pressure < 1 MPa
        # is packer inflation / line filling, not injection into the test interval -> not a stage
        if t[j] - t[i] >= DMIN and np.median(d[i:j + 1, 2]) >= PMIN: st.append((i, j))
        i = j + 1
    return st

def events():
    out = []
    with open(os.path.join(DATA, 'aspo_catalog_kwiatek2018niemz2020_match.csv')) as f:
        next(f)
        for line in f:
            s = line.split(';')
            tt = dt.datetime.strptime(s[0], '%Y-%m-%d_%H:%M:%S.%f').replace(tzinfo=dt.timezone.utc).timestamp()
            out.append((tt, float(s[4])))
    return np.array(out)

def dieterich_rate(t, p, Asig):
    """R/r for S(t) = p(t) (MPa), exact per-step update of the Dieterich state gamma."""
    tr = Asig / TA
    g = 1.0 / tr
    R = np.empty(len(t))
    R[0] = 1.0
    for k in range(1, len(t)):
        h = t[k] - t[k - 1]
        sd = (p[k] - p[k - 1]) / h
        if abs(sd) < 1e-12: g = g + h / Asig
        else:
            x = -sd * h / Asig
            g = (g - 1.0 / sd) * math.exp(min(x, 700.0)) + 1.0 / sd
            if g <= 0: g = 1e-300
        R[k] = 1.0 / (g * tr)
    return R

TESTS = ['HF1', 'HF2', 'HF3', 'HF4', 'HF5', 'HF6']
SCHEME = {'HF1': 'conv AG', 'HF2': 'conv AG', 'HF3': 'progressive AG', 'HF4': 'conv fgDG', 'HF5': 'pulse fgDG', 'HF6': 'conv fgG'}

def build():
    ev = events()
    S = []   # stage records
    raw = {}
    for name in TESTS:
        raw[name] = load_test(name)
        for part, d in enumerate(raw[name]):
            st = stages(d)
            for k, (i, j) in enumerate(st):
                t, q, p = d[:, 0], d[:, 1], d[:, 2]
                t0, t1 = t[i], t[j]
                tnext = t[st[k + 1][0]] if k + 1 < len(st) else t[-1] + POST
                tw = min(tnext, t1 + POST)
                sl = slice(i, j + 1)
                V = float(np.sum(q[sl]) / 60.0)                   # L
                E = float(np.sum(q[sl] * p[sl]) / 60.0)          # L*MPa = kJ
                co = int(np.sum((ev[:, 0] >= t0) & (ev[:, 0] < t1)))
                post = int(np.sum((ev[:, 0] >= t1) & (ev[:, 0] < tw)))
                S.append(dict(test=name, part=part, k=k, t0=t0, t1=t1, tw=tw, i=i, j=j, dur=float(t1 - t0), V=V, E=E,
                              pmed=float(np.median(p[sl])), pmax=float(p[sl].max()), qmean=float(np.mean(q[sl])), co=co, post=post))
    return ev, S, raw

def kernel_shapes(S, raw, kind, Asig=None):
    """unnormalised expected counts per stage (co+post), and the co-only part."""
    out, co = [], []
    if kind == 'K1':
        cache = {}
        for s in S:
            key = (s['test'], s['part'])
            if key not in cache:
                d = raw[s['test']][s['part']]
                cache[key] = (d[:, 0], dieterich_rate(d[:, 0], np.maximum(d[:, 2], 0.0), Asig))
            t, R = cache[key]
            w = (t >= s['t0']) & (t < s['tw'])
            c = (t >= s['t0']) & (t < s['t1'])
            out.append(float(R[w].sum())); co.append(float(R[c].sum()))
        return np.array(out), np.array(co)
    if kind == 'K2': x = np.array([s['E'] for s in S])
    elif kind == 'K3': x = np.array([s['V'] for s in S])
    elif kind == 'K4':
        x = []
        mx = {}
        for s in S:
            m = mx.get(s['test'], 0.0)
            x.append(max(0.0, s['E'] - m)); mx[s['test']] = max(m, s['E'])
        x = np.array(x)
        x = np.maximum(x, 0.01 * x.mean())        # fixed 1 % floor (zero expectations would be infinitely penalised)
    return x, None

def poisson_fit(S, shape, tests):
    """per-test MLE rate factor; returns loglik, factors."""
    ll, fac = 0.0, {}
    for tn in tests:
        idx = [i for i, s in enumerate(S) if s['test'] == tn]
        n = np.array([S[i]['co'] + S[i]['post'] for i in idx], float)
        m = shape[idx]
        a = n.sum() / max(m.sum(), 1e-300)
        fac[tn] = a
        mu = np.maximum(a * m, 1e-300)
        ll += float(np.sum(n * np.log(mu) - mu - [math.lgamma(v + 1) for v in n]))
    return ll, fac

def deviance(S, shape, tn):
    idx = [i for i, s in enumerate(S) if s['test'] == tn]
    n = np.array([S[i]['co'] + S[i]['post'] for i in idx], float)
    m = shape[idx]; mu = np.maximum(n.sum() * m / max(m.sum(), 1e-300), 1e-300)
    return float(2 * np.sum(np.where(n > 0, n * np.log(np.maximum(n, 1e-300) / mu), 0) - (n - mu))), int(len(idx))

def main():
    ev, S, raw = build()
    res = {'stages': S}
    print('stage table:')
    for s in S:
        print(f"  {s['test']}.{s['part']}.{s['k']}  dur {s['dur']:5.0f}s  V {s['V']:6.2f} L  E {s['E']:7.1f} kJ  p_med {s['pmed']:5.2f} p_max {s['pmax']:5.2f}  co {s['co']:3d} post {s['post']:3d}")
    tot = {tn: sum(s['co'] + s['post'] for s in S if s['test'] == tn) for tn in TESTS}
    print('events assigned per test', tot, 'of', len(ev))
    res['assigned'] = tot
    train = ['HF1', 'HF2']
    # --- H6a
    out = {}
    grid = np.exp(np.linspace(math.log(0.02), math.log(20), 61))
    best = None
    prof = []
    for A in grid:
        sh, _ = kernel_shapes(S, raw, 'K1', A)
        ll, fac = poisson_fit(S, sh, train)
        prof.append((float(A), ll))
        if best is None or ll > best[1]: best = (A, ll, fac)
    # refine (golden section on log A)
    lo, hi = math.log(best[0]) - 0.12, math.log(best[0]) + 0.12
    f = lambda la: -poisson_fit(S, kernel_shapes(S, raw, 'K1', math.exp(la))[0], train)[0]
    for _ in range(25):
        m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
        if f(m1) < f(m2): hi = m2
        else: lo = m1
    Abest = math.exp(0.5 * (lo + hi))
    sh1, co1 = kernel_shapes(S, raw, 'K1', Abest)
    ll1, fac1 = poisson_fit(S, sh1, train)
    out['K1'] = dict(Asig=Abest, ll=ll1, k=len(train) + 1, aic=2 * (len(train) + 1) - 2 * ll1, profile=prof)
    shapes = {'K1': sh1}
    for K in ['K2', 'K3', 'K4']:
        sh, _ = kernel_shapes(S, raw, K)
        ll, fac = poisson_fit(S, sh, train)
        out[K] = dict(ll=ll, k=len(train), aic=2 * len(train) - 2 * ll)
        shapes[K] = sh
    # null: constant per stage
    sh0 = np.ones(len(S)); ll0, _ = poisson_fit(S, sh0, train)
    out['K0_const'] = dict(ll=ll0, k=len(train), aic=2 * len(train) - 2 * ll0)
    shd = np.array([s['dur'] for s in S]); lld, _ = poisson_fit(S, shd, train)
    out['K5_duration'] = dict(ll=lld, k=len(train), aic=2 * len(train) - 2 * lld)
    shapes['K0_const'] = sh0; shapes['K5_duration'] = shd
    # EXPLORATORY (not pre-registered): volume x mean-rate^beta - does the rate matter at fixed volume?
    Vv = np.array([s['V'] for s in S]); qv = np.array([s['qmean'] for s in S])
    bprof = []
    for b in np.linspace(-1.5, 2.5, 81):
        bprof.append((float(b), poisson_fit(S, Vv * qv ** b, train)[0]))
    bb = max(bprof, key=lambda x: x[1])
    llb = bb[1]
    ci = [x[0] for x in bprof if x[1] >= llb - 1.92]
    out['K6_V_qbeta'] = dict(beta=bb[0], ll=llb, k=len(train) + 1, aic=2 * (len(train) + 1) - 2 * llb, beta_ci95=[min(ci), max(ci)])
    shapes['K6_V_qbeta'] = Vv * qv ** bb[0]
    # confound check: later stages have both higher rate and a larger pre-existing fracture (cumulative volume)
    vc, acc = [], {}
    for s_ in S:
        vc.append(acc.get(s_['test'], 0.0)); acc[s_['test']] = acc.get(s_['test'], 0.0) + s_['V']
    vc = np.array(vc)
    gprof = [(float(g), poisson_fit(S, Vv * (1 + vc) ** g, train)[0]) for g in np.linspace(-1, 2.5, 71)]
    gg = max(gprof, key=lambda x: x[1])
    out['K7_V_cumV'] = dict(beta=gg[0], ll=gg[1], k=len(train) + 1, aic=2 * (len(train) + 1) - 2 * gg[1],
                            beta_ci95=[min(x[0] for x in gprof if x[1] >= gg[1] - 1.92), max(x[0] for x in gprof if x[1] >= gg[1] - 1.92)])
    shapes['K7_V_cumV'] = Vv * (1 + vc) ** gg[0]
    # both together: rate exponent at fixed cumulative-volume exponent (profile over a 2-D grid)
    best2 = None
    for b in np.linspace(-1, 2, 31):
        for g in np.linspace(-1, 2.5, 36):
            l_ = poisson_fit(S, Vv * qv ** b * (1 + vc) ** g, train)[0]
            if best2 is None or l_ > best2[2]: best2 = (float(b), float(g), l_)
    profb = {}
    for b in np.linspace(-2, 3, 101):
        profb[round(float(b), 3)] = max(poisson_fit(S, Vv * qv ** b * (1 + vc) ** g, train)[0] for g in np.linspace(-2, 3, 101))
    LJ = max(profb.values())
    bci = [min(b for b, l in profb.items() if l >= LJ - 1.92), max(b for b, l in profb.items() if l >= LJ - 1.92)]
    out['K8_V_q_cumV'] = dict(beta=best2[0], gamma=best2[1], ll=best2[2], k=len(train) + 2, aic=2 * (len(train) + 2) - 2 * best2[2], beta_ci95=bci)
    print('  rate/cumulative-volume correlation (train stages): %.2f' % np.corrcoef(qv[[i for i, s_ in enumerate(S) if s_['test'] in train]], vc[[i for i, s_ in enumerate(S) if s_['test'] in train]])[0, 1])
    print('\nH6a train HF1+HF2:  (Dieterich A*sigma = %.3f MPa)' % Abest)
    for K, v in out.items(): print(f"  {K:12s} logL {v['ll']:8.2f}  AIC {v['aic']:7.2f}" + (f"  beta {v['beta']:.2f} CI {v['beta_ci95']}" if 'beta' in v else '') + (f" gamma {v['gamma']:.2f}" if 'gamma' in v else ''))
    aics = {K: v['aic'] for K, v in out.items() if K[:2] not in ('K6', 'K7', 'K8')}   # ranking over the pre-registered set
    order = sorted(aics, key=aics.get)
    pref = order[0] if aics[order[1]] - aics[order[0]] >= 2 else None
    print('  ranking', order, '-> preferred', pref)
    oos = {}
    for K in shapes:
        dv = {tn: deviance(S, shapes[K], tn) for tn in ['HF6', 'HF3', 'HF4']}
        oos[K] = dv
        print(f"  OOS deviance {K:12s} HF6 {dv['HF6'][0]:6.2f} (n={dv['HF6'][1]})  HF3 {dv['HF3'][0]:6.2f}")
    res['H6a'] = dict(models={K: {k: v for k, v in d.items() if k != 'profile'} for K, d in out.items()}, K1_profile=out['K1']['profile'],
                      ranking=order, preferred=pref, oos=oos)
    # --- H6b: per-test rate factors under the preferred (or best) kernel
    K = order[0]
    shK = shapes[K]
    allfac = poisson_fit(S, shK, TESTS)[1]
    conv = ['HF1', 'HF2', 'HF4', 'HF6']
    lf = np.log([max(allfac[t], 1e-12) for t in conv])
    ag = np.log([allfac['HF1'], allfac['HF2']])
    sig_all = float(np.std(lf, ddof=1))
    mu_ag = float(np.mean(ag))
    # predictive (t, df = n-1) for a new conventional test in AG, spread from all conventional tests
    from math import sqrt
    tq = {3: 2.353}[len(conv) - 1]       # one-sided 95 % t quantile, df 3
    lo5 = mu_ag - tq * sig_all * sqrt(1 + 1 / 2)
    z = (math.log(max(allfac['HF3'], 1e-12)) - mu_ag) / (sig_all * sqrt(1 + 1 / 2))
    print(f"\nH6b kernel {K}: factors " + ', '.join(f"{t} {allfac[t]:.3g}" for t in TESTS))
    print(f"  conventional log-spread sigma {sig_all:.2f};  HF3 / AG-mean = {math.exp(math.log(allfac['HF3']) - mu_ag):.3f};  t-score {z:.2f}; 5% bound {math.exp(lo5 - mu_ag):.3f} x AG mean")
    res['H6b'] = dict(kernel=K, factors=allfac, sigma=sig_all, ratio_HF3=math.exp(math.log(allfac['HF3']) - mu_ag),
                      t=z, lower5_ratio=math.exp(lo5 - mu_ag), significant=bool(math.log(allfac['HF3']) < lo5))
    # raw-count version (no kernel normalisation)
    # --- H1h: post-shut-in fraction, observed vs Dieterich
    po = sum(s['post'] for s in S if s['test'] in train); to = sum(s['co'] + s['post'] for s in S if s['test'] in train)
    idx = [i for i, s in enumerate(S) if s['test'] in train]
    fpred = float(np.sum((sh1 - co1)[idx] * np.array([fac1[S[i]['test']] for i in idx])) / np.sum(sh1[idx] * np.array([fac1[S[i]['test']] for i in idx])))
    print(f"\nH1h post-shut-in fraction train: observed {po}/{to} = {po/to:.3f}; Dieterich (A*sigma {Abest:.3f}) predicts {fpred:.3f}")
    # binomial z
    pz = (po / to - fpred) / math.sqrt(max(fpred * (1 - fpred) / to, 1e-12))
    res['H1h'] = dict(observed=po / to, n=to, post=po, dieterich=fpred, z=pz)
    # all tests
    po2 = sum(s['post'] for s in S); to2 = sum(s['co'] + s['post'] for s in S)
    res['H1h']['all_tests'] = dict(post=po2, n=to2, observed=po2 / max(to2, 1))
    json.dump(res, open(os.path.join(HERE, 'h6_results.json'), 'w'), indent=1, default=float)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
        P = np.array(out['K1']['profile']); ax[0].semilogx(P[:, 0], P[:, 1] - P[:, 1].max(), 'k-')
        ax[0].set_xlabel(r'$A\sigma$ (MPa)'); ax[0].set_ylabel(r'$\Delta$ log L (HF1+HF2)'); ax[0].set_ylim(-30, 1); ax[0].grid(alpha=.3)
        n = np.array([s['co'] + s['post'] for s in S]); m = [s['test'] in train for s in S]
        for KK, c in [('K1', 'C3'), ('K2', 'C0'), ('K3', 'C2')]:
            f_ = poisson_fit(S, shapes[KK], TESTS)[1]
            mu = np.array([f_[s['test']] * shapes[KK][i] for i, s in enumerate(S)])
            ax[1].plot(mu[m], n[m], 'o', color=c, ms=4, label=KK)
        ax[1].plot([0, 60], [0, 60], 'k:'); ax[1].set_xlabel('expected stage count'); ax[1].set_ylabel('observed'); ax[1].legend(); ax[1].grid(alpha=.3)
        ax[2].bar(TESTS, [allfac[t] for t in TESTS], color=['C0', 'C0', 'C3', 'C0', 'C1', 'C0'])
        ax[2].set_yscale('log'); ax[2].set_ylabel(f'rate factor ({K})'); ax[2].grid(alpha=.3, axis='y')
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'h6_fig.png'), dpi=110)
    except Exception as e:
        print('no figure:', e)

if __name__ == '__main__':
    main()
