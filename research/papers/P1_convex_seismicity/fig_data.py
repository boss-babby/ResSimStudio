"""Curve data for the P1 figures (synthetic experiments). The figures themselves are rendered by ResSim Studio's
publication figure engine (app/src/figure.cpp) through tools/paper_figs.cpp - this script only computes numbers.
Runs in ~1 min; writes fig_synth.json next to this file."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "exp"))
import h1_seis_opt as H

DAY = H.DAY
t0 = time.time()
p = dict(H.P)
g = H.build(p)
S = H.schedules(g, p)
qo, res = H.optimise(g, p)
S["optimal"] = qo


def history(q):
    """time (d), cumulative excess events E(t) = N/r - t (days of background), normalised rate R/r, fault dS (MPa)"""
    ta, a, s = g["ta"], g["a"], g["s"]
    dS = g["A"][0] @ q
    K = np.exp(s / ta + dS / a)
    I = np.concatenate([[0], np.cumsum((K[1:] + K[:-1]) / 2 * np.diff(s))])
    R = K / (1 + I / ta)
    N = ta * np.log1p(I / ta) - s
    return s / DAY, N / DAY, R, dS / 1e6


out = dict(Ti_d=g["Ti"] / DAY, Th_d=g["Th"] / DAY, t_ctrl_d=((np.arange(g["n"]) + 0.5) * g["dt"] / DAY).tolist(), sched={})
iTi = int(np.searchsorted(g["s"], g["Ti"]))
rank = {}
for k in ["constant", "step-up(4)", "taper-last30%", "optimal"]:
    t, N, R, dS = history(S[k])
    st = max(1, len(t) // 1500)
    out["sched"][k] = dict(q=(S[k] * 86400).tolist(), t=t[::st].tolist(), E=N[::st].tolist(), R=R[::st].tolist(), dS=dS[::st].tolist(),
                           E_shutin=float(N[iTi]), E_total=float(N[-1]))
    rank[k] = (float(N[iTi]), float(N[-1]))
out["rank"] = rank

# convexity along segments between random feasible schedules (Theorem 1)
rng = np.random.default_rng(7)
segs = []
qmax = p["qmax_fac"] * p["V"] / g["Ti"]
for i in range(3):
    qs = []
    for _ in range(2):
        q = rng.random(g["n"]) ** 3
        q *= p["V"] / (q.sum() * g["dt"]); q = np.minimum(q, qmax * 0.999); q *= p["V"] / (q.sum() * g["dt"])
        qs.append(q)
    lam = np.linspace(0, 1, 41)
    E = [H.excess((1 - l) * qs[0] + l * qs[1], g) / DAY for l in lam]
    segs.append(dict(lam=lam.tolist(), E=E, chord=[(1 - l) * E[0] + l * E[-1] for l in lam]))
out["segments"] = segs
# 2000 random midpoint tests (as in exp/h1_verify.py)
viol = 0; worst = 0.0
for i in range(2000):
    a_ = rng.random(g["n"]) ** 2; b_ = rng.random(g["n"]) ** 2
    a_ *= p["V"] / (a_.sum() * g["dt"]); b_ *= p["V"] / (b_.sum() * g["dt"])
    Ea, Eb, Em = H.excess(a_, g), H.excess(b_, g), H.excess(0.5 * (a_ + b_), g)
    gap = (0.5 * (Ea + Eb) - Em) / max(1e-300, abs(0.5 * (Ea + Eb)))
    worst = min(worst, gap); viol += gap < -1e-12
out["midpoint"] = dict(tests=2000, violations=int(viol), worst_rel_gap=float(worst))
out["runtime_s"] = time.time() - t0
json.dump(out, open(os.path.join(HERE, "fig_synth.json"), "w"))
print(json.dumps(dict(rank=rank, midpoint=out["midpoint"], runtime=out["runtime_s"]), indent=1))
