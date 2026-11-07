#!/usr/bin/env python3
"""H1d (pre-registered in CYCLE3_PLAN.md): magnitude-aware objectives.

D1  Constant b (per patch): Lambda_m(q) = sum_k N_k(q) 10^{-b_k (m - Mc_k)} is a non-negative weighted sum of convex functions
    -> convex; P(M >= m) = 1 - exp(-Lambda_m) is monotone in Lambda_m -> quasi-convex with the same minimiser; a chance
    constraint P <= alpha  <=>  Lambda_m <= -ln(1 - alpha) is a convex constraint.  Numerical check: midpoint test.
D2  Stress-dependent b: b(t) = b0 - beta * dS(t)  ->  exceedance density  lambda_m = R(t) 10^{-b0 dm} exp(gamma dS),
    gamma = beta dm ln10.  Is Lambda_m still convex in q?  2000 random midpoint pairs per configuration, incl. pulse pairs.
Kernel: Cycle-1 3-D diffusion point source (exp/h1_seis_opt.py), one patch at 300 m, 10-day injection, horizon 5x.
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h1_seis_opt import P, build, DAY

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(7)
t0 = time.time()


def lam_exceed(q, g, gamma, ta):
    """expected exceedances / (r 10^{-b0 dm}) over the horizon: int R/r * exp(gamma dS) dt, Dieterich with depletion."""
    s, w, a = g["s"], g["w"], g["a"]
    dS = g["A"][0] @ q
    x = s / ta + dS / a                              # ln of the Dieterich numerator (tau_r t / a = t / t_a)
    m = x.max()
    e = np.exp(x - m)
    I = np.concatenate([[0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * np.diff(s))]) * np.exp(m)    # int_0^t e^x
    Rr = np.exp(x) / (1.0 + I / ta)                  # R / r
    return float(np.sum(w * Rr * np.exp(gamma * dS)))


def feasible(n, V, dt, qmax, kind):
    if kind == "random":
        q = rng.gamma(0.5, 1.0, n)
    else:                                            # concentrated pulse somewhere
        q = np.zeros(n); j = rng.integers(0, n); w = rng.integers(1, max(2, n // 10)); q[j:j + w] = 1.0
    q *= V / (q.sum() * dt)
    return np.minimum(q, qmax)


res = {"D1": {}, "D2": []}
p = dict(P); p["n"] = 60; p["neval"] = 1200
for ta_yr in (100.0, 0.05):                          # t << t_a  and strong population depletion
    p["ta_yr"] = ta_yr
    g = build(p)
    ta = g["ta"]
    V, dt = p["V"], g["dt"]
    qmax = p["qmax_fac"] * V / g["Ti"]
    S_const = (g["A"][0] @ np.full(g["n"], V / g["Ti"])).max()
    for beta in (0.0, 0.1, 0.25, 0.5):               # per MPa
        for dm in (1.0, 2.0, 3.0):
            gamma = beta * dm * np.log(10) / 1e6     # per Pa
            viol, worst, npairs = 0, 0.0, 0
            for it in range(2000):
                kind = "pulse" if it % 2 else "random"
                q1, q2 = feasible(g["n"], V, dt, qmax, kind), feasible(g["n"], V, dt, qmax, "random" if it % 3 else "pulse")
                f1, f2, fm = lam_exceed(q1, g, gamma, ta), lam_exceed(q2, g, gamma, ta), lam_exceed(0.5 * (q1 + q2), g, gamma, ta)
                gap = (fm - 0.5 * (f1 + f2)) / max(abs(0.5 * (f1 + f2)), 1e-300)
                npairs += 1
                if gap > 1e-9:
                    viol += 1; worst = max(worst, gap)
            row = dict(ta_yr=ta_yr, beta_per_MPa=beta, dm=dm, gamma_dSmax=gamma * S_const, pairs=npairs, violations=viol, worst_rel_gap=worst)
            res["D2"].append(row)
            print(f"t_a {ta_yr:6.2f} yr  beta {beta:4.2f}/MPa  dm {dm:.0f}  gamma*dSmax {gamma * S_const:6.2f}: violations {viol}/{npairs}  worst {worst:.2e}")
# D2b synthetic stress test beyond the physical range: gamma * dSmax = 1 ... 30
res["D2b_synthetic"] = []
p = dict(P); p["n"] = 60; p["neval"] = 1200
for ta_yr in (100.0, 0.05, 0.002):
    p["ta_yr"] = ta_yr; g = build(p); ta = g["ta"]; V, dt = p["V"], g["dt"]; qmax = p["qmax_fac"] * V / g["Ti"]
    S_const = (g["A"][0] @ np.full(g["n"], V / g["Ti"])).max()
    for gds in (1.0, 3.0, 10.0, 30.0):
        gamma = gds / S_const
        viol, worst = 0, 0.0
        for it in range(2000):
            q1 = feasible(g["n"], V, dt, qmax, "pulse" if it % 2 else "random"); q2 = feasible(g["n"], V, dt, qmax, "random" if it % 3 else "pulse")
            f1, f2, fm = lam_exceed(q1, g, gamma, ta), lam_exceed(q2, g, gamma, ta), lam_exceed(0.5 * (q1 + q2), g, gamma, ta)
            gap = (fm - 0.5 * (f1 + f2)) / max(abs(0.5 * (f1 + f2)), 1e-300)
            if gap > 1e-9: viol += 1; worst = max(worst, gap)
        res["D2b_synthetic"].append(dict(ta_yr=ta_yr, gamma_dSmax=gds, violations=viol, worst_rel_gap=worst))
        print(f"[synthetic] t_a {ta_yr:6.3f} yr  gamma*dSmax {gds:5.1f}: violations {viol}/2000  worst {worst:.2e}")
# D1 numerical: two patches with different b and Mc -> weighted sum, plus quasi-convexity of P
p["ta_yr"] = 100.0; p["rpatch"] = [150.0, 400.0]
g = build(p); V, dt = p["V"], g["dt"]; qmax = p["qmax_fac"] * V / g["Ti"]
from h1_seis_opt import excess
wts = [10 ** (-1.0 * (2.0 - 0.0)), 10 ** (-1.4 * (2.0 - 0.5))]


def Lam(q):
    tot = 0.0
    for k, M in enumerate(g["A"]):
        gk = dict(g); gk["A"] = [M]
        tot += wts[k] * excess(q, gk)
    return tot


viol = 0
for it in range(2000):
    q1, q2 = feasible(g["n"], V, dt, qmax, "random"), feasible(g["n"], V, dt, qmax, "pulse")
    L1, L2, Lm = Lam(q1), Lam(q2), Lam(0.5 * (q1 + q2))
    if Lm > 0.5 * (L1 + L2) * (1 + 1e-12):
        viol += 1
    # quasi-convexity of P = 1 - exp(-L): P(mid) <= max(P1, P2)
    if 1 - np.exp(-Lm) > max(1 - np.exp(-L1), 1 - np.exp(-L2)) + 1e-15:
        viol += 1
res["D1"] = dict(pairs=2000, violations=viol)
print(f"D1 weighted two-patch exceedance + quasi-convexity of P: violations {viol}/4000")
res["runtime_s"] = time.time() - t0
json.dump(res, open(os.path.join(HERE, "h1d_results.json"), "w"), indent=1)
print(f"done in {res['runtime_s']:.0f} s")
