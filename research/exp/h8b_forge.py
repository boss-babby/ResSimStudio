#!/usr/bin/env python3
"""Research Cycle 8 - H8b: Utah FORGE 16A(78)-32 April 2022 stimulation (3 stages). Pre-registered in CYCLE8_PLAN.md.
Kernels for 5-min binned event counts (Poisson, one rate factor per stage): volume, energy, volume with exponential
memory (tau profiled), Dieterich (A sigma profiled, S = surface treating pressure), constant. numpy only.
usage: h8b_forge.py [shift_min=-60]   -> prints the table and writes h8b_results.json (primary shift only)"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data', 'forge')
SHIFT = float(sys.argv[1]) if len(sys.argv) > 1 else -60.0
MC, BIN, POST, QON = -1.5, 300.0, 3 * 3600.0, 1.0
PSI = 0.00689476          # MPa
TA = 1e9                  # Dieterich aftershock time (s): t << ta

def dieterich(t, S, Asig):
    tr = Asig / TA; g = 1.0 / tr; R = np.empty(len(t)); R[0] = 1.0
    for k in range(1, len(t)):
        h = t[k] - t[k - 1]
        sd = (S[k] - S[k - 1]) / h
        if abs(sd) < 1e-12: g = g + h / Asig
        else:
            x = -sd * h / Asig
            g = (g - 1.0 / sd) * math.exp(min(x, 700.0)) + 1.0 / sd
            if not math.isfinite(g) or g > 1e300: g = 1e300
            if g <= 0: g = 1e-300
        R[k] = 1.0 / (g * tr)
    return R

def load():
    cat = np.genfromtxt(os.path.join(DATA, 'forge16a_catalogue_ges_v2.csv'), delimiter=';', skip_header=1)
    st = []
    for i in (1, 2, 3):
        pm = np.genfromtxt(os.path.join(DATA, f'forge16a_stage{i}_pump.csv'), delimiter=';', skip_header=1)
        t = pm[:, 0] + SHIFT * 60.0
        q = np.maximum(pm[:, 2], 0.0) * 0.158987 / 60.0          # bpm -> m3/s
        p = pm[:, 1] * PSI
        on = np.where(pm[:, 2] > QON)[0]
        t0, t1 = t[on[0]], t[on[-1]]
        ev = cat[(cat[:, 3] == i) & (cat[:, 2] == 2) & (cat[:, 1] >= MC), 0]
        tend = min(t1 + POST, ev.max() if len(ev) else t1 + POST)
        edges = np.arange(t0, tend + 1e-6, BIN)
        n, _ = np.histogram(ev, bins=edges)
        # 1 s series extended with zero rate after the record (pressure held at the last value)
        dt1 = 1.0
        tt = np.arange(t0, edges[-1], dt1)
        qq = np.interp(tt, t, q, left=0.0, right=0.0)
        pp = np.interp(tt, t, p, left=p[0], right=p[-1])
        st.append(dict(stage=i, t0=t0, t1=t1, edges=edges, n=n.astype(float), tt=tt, qq=qq, pp=pp,
                       pre=int(np.sum((ev >= t0 - 3600) & (ev < t0))), co=int(np.sum((ev >= t0) & (ev <= t1))),
                       post=int(np.sum((ev > t1) & (ev < tend))), V=float(np.sum(qq) * dt1)))
    return st

def binsum(s, x):
    idx = np.minimum(((s['tt'] - s['edges'][0]) // BIN).astype(int), len(s['n']) - 1)
    return np.bincount(idx, weights=x, minlength=len(s['n']))

def ll_stage(n, m):
    a = n.sum() / max(m.sum(), 1e-300)
    mu = np.maximum(a * m, 1e-300)
    return float(np.sum(n * np.log(mu) - mu) - sum(math.lgamma(v + 1) for v in n))

def floor(m): return np.maximum(m, 0.01 * m.mean())

def shapes(s, kind, par=None):
    if kind == 'KC': return np.ones(len(s['n']))
    if kind == 'KV': return floor(binsum(s, s['qq']))
    if kind == 'KE': return floor(binsum(s, s['qq'] * s['pp']))
    if kind == 'KVtau':
        a = math.exp(-1.0 / par); m = np.empty(len(s['qq'])); acc = 0.0
        for k, v in enumerate(s['qq']): acc = acc * a + v; m[k] = acc
        return np.maximum(binsum(s, m), 1e-12)
    if kind == 'KD':
        step = 10
        tt, pp = s['tt'][::step], s['pp'][::step]
        R = dieterich(tt, pp, par)
        Rf = np.repeat(R, step)[:len(s['tt'])]
        return binsum(s, Rf)

def main():
    st = load()
    out = dict(shift_min=SHIFT, Mc=MC, bin_s=BIN, stages=[dict(stage=s['stage'], pre_1h=s['pre'], co=s['co'], post=s['post'],
               bins=len(s['n']), V_m3=round(s['V'], 1)) for s in st], kernels={})
    for s in st: print(f"stage {s['stage']}: {len(s['n'])} bins, events pre-1h {s['pre']}, co {s['co']}, post {s['post']}, V {s['V']:.0f} m3")
    for K in ('KC', 'KV', 'KE'):
        ll = sum(ll_stage(s['n'], shapes(s, K)) for s in st)
        out['kernels'][K] = dict(ll=ll, k=3, aic=2 * 3 - 2 * ll)
    for K, grid in (('KVtau', [60, 120, 300, 600, 1200, 1800, 3600, 5400, 7200, 10800, 14400, 21600, 43200]),
                    ('KD', [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 300.0])):
        prof = []
        for g in grid:
            ll = sum(ll_stage(s['n'], shapes(s, K, g)) for s in st)
            prof.append((g, ll))
        g, ll = max(prof, key=lambda x: x[1])
        out['kernels'][K] = dict(par=g, ll=ll, k=4, aic=2 * 4 - 2 * ll, profile=prof)
    ks = out['kernels']
    for K, v in ks.items(): print(f"  {K:6s} logL {v['ll']:11.2f}  AIC {v['aic']:10.2f}" + (f"  par {v['par']}" if 'par' in v else ''))
    best_vol = min(('KV', 'KE', 'KVtau'), key=lambda K: ks[K]['aic'])
    d = ks['KD']['aic'] - ks[best_vol]['aic']
    verdict = 'SUPPORTED (volume-type)' if d >= 2 else 'REFUTED (Dieterich)' if d <= -2 else 'INCONCLUSIVE'
    out['best_volume_kernel'] = best_vol; out['dAIC_D_minus_bestvol'] = d; out['verdict'] = verdict
    # post-shut-in fraction: observed vs predicted by each fitted kernel (as in F16)
    for K in ('KV', 'KE', 'KVtau', 'KD'):
        fr = []
        for s in st:
            m = shapes(s, K, ks[K].get('par'))
            post = s['edges'][:-1] >= s['t1']
            fr.append(float(m[post].sum() / m.sum()))
        out['kernels'][K]['post_fraction_pred'] = fr
    out['post_fraction_obs'] = [float(s['n'][s['edges'][:-1] >= s['t1']].sum() / max(s['n'].sum(), 1)) for s in st]
    print(f"  best volume-type {best_vol}; AIC(KD) - AIC({best_vol}) = {d:.1f} -> {verdict}")
    print('  post-shut-in fraction obs', [round(x, 3) for x in out['post_fraction_obs']],
          '| KD', [round(x, 3) for x in ks['KD']['post_fraction_pred']], '| KVtau', [round(x, 3) for x in ks['KVtau']['post_fraction_pred']])
    if SHIFT == -60.0:
        json.dump(out, open(os.path.join(HERE, 'h8b_results.json'), 'w'), indent=1)

if __name__ == '__main__':
    main()
