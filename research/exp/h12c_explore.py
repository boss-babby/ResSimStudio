#!/usr/bin/env python3
"""EXPLORATORY (after H12b; not pre-registered), CYCLE13_PLAN.md.
E13a: the Kaiser LP optimum is non-unique (withdrawal timing hardly changes new maxima). Lexicographic check: among Kaiser-optimal
      schedules (N_K <= N_K*(1 + 1e-4)), minimise N_D, giving the fairest X_D.
E13b: robust design across memory models (NEW_FIELDS N2): minimise max(N_D/N_D*, N_K/N_K*), convex as a max of convex functions.
usage: h12c_explore.py CLI  -> exp/h12c_results.json
"""
import json, os, random, sys
sys.argv = [sys.argv[0], sys.argv[1], "none"]
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "h12_storage.py")).read().split("calibrate()\nif")[0]
exec(src)
import numpy as np
from scipy.optimize import minimize, linprog

G = [run("SW", "pulse:%d" % k, hist=True)[1] for k in range(8)]
_, H0 = run("SW", [0.0] * 8, hist=True)
faces = sorted(G[0]); times = [t for (t, _, _) in G[0][faces[0]]]; area = G[0][faces[0]][0][2]
G0 = np.array([[H0[f][j][1] for j in range(len(times))] for f in faces])
Gm = np.array([[[G[k][f][j][1] for j in range(len(times))] for f in faces] for k in range(8)]) - G0[None]
nf, nt = len(faces), len(times)


def D(u):
    return G0 + np.tensordot(np.asarray(u, float), Gm, axes=1)


def ND(u, a):
    Dm = D(u)
    return dieterich({f: [(times[j], Dm[i, j], area) for j in range(nt)] for i, f in enumerate(faces)}, a)


def NK(u):
    return float(np.sum(np.maximum(0.0, D(u).max(axis=1))))


# Kaiser LP optimum
c = np.r_[np.zeros(8), np.ones(nf)]
Aub = []; bub = []
for i in range(nf):
    for j in range(nt):
        row = np.zeros(8 + nf); row[:8] = Gm[:, i, j]; row[8 + i] = -1; Aub.append(row); bub.append(-G0[i, j])
Aub = np.array(Aub); bub = np.array(bub)
Aeq = np.zeros((2, 8 + nf)); Aeq[0, :4] = 1; Aeq[1, 4:8] = 1
bounds = [(0, 2)] * 4 + [(-2, 0)] * 4 + [(0, None)] * nf
lp = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[4, -4], bounds=bounds, method="highs")
NKs = lp.fun
BASE8 = [1, 1, 1, 1, -1, -1, -1, -1]
out = {"NK_star": NKs, "levels": {}}
for a in (1.0, 10.0):
    # Dieterich optimum (as H12b)
    cons_eq = [{"type": "eq", "fun": lambda z: sum(z[:4]) - 4}, {"type": "eq", "fun": lambda z: sum(z[4:8]) + 4}]
    sc = ND(BASE8, a)
    best = None
    for s in [BASE8, list(lp.x[:8])] + [sched(random.Random(7 + i)) for i in range(4)]:
        o = minimize(lambda u: ND(u, a) / sc, s, method="SLSQP", bounds=[(0, 2)] * 4 + [(-2, 0)] * 4, constraints=cons_eq,
                     options={"maxiter": 200, "ftol": 1e-12})
        best = o if best is None or o.fun < best.fun else best
    uD = list(best.x); NDs = ND(uD, a)
    # E13a lexicographic: variables z = (u, s); min ND(u) s.t. Aub z <= bub, sum s <= NK*(1+1e-4)
    zc = [{"type": "ineq", "fun": lambda z: bub - Aub @ z},
          {"type": "ineq", "fun": lambda z: NKs * (1 + 1e-4) - np.sum(z[8:])}] + \
         [{"type": "eq", "fun": lambda z: np.sum(z[:4]) - 4}, {"type": "eq", "fun": lambda z: np.sum(z[4:8]) + 4}]
    z0 = np.array(lp.x)
    o = minimize(lambda z: ND(z[:8], a) / sc, z0, method="SLSQP", bounds=bounds, constraints=zc,
                 options={"maxiter": 300, "ftol": 1e-12})
    uKlex = list(o.x[:8])
    XD_lex = ND(uKlex, a) / NDs - 1
    # E13b robust: min t s.t. ND(u)/NDs <= t, NK(u)/NKs <= t  (NK via epigraph s)
    rc = [{"type": "ineq", "fun": lambda z: z[-1] - ND(z[:8], a) / NDs},
          {"type": "ineq", "fun": lambda z: z[-1] - np.sum(z[8:8 + nf]) / NKs},
          {"type": "ineq", "fun": lambda z: bub - Aub @ z[:-1]},
          {"type": "eq", "fun": lambda z: np.sum(z[:4]) - 4}, {"type": "eq", "fun": lambda z: np.sum(z[4:8]) + 4}]
    z0 = np.r_[np.array(uKlex), np.maximum(0, D(uKlex).max(axis=1)), 1.1]
    o2 = minimize(lambda z: z[-1], z0, method="SLSQP", bounds=bounds + [(1, None)], constraints=rc,
                  options={"maxiter": 300, "ftol": 1e-12})
    uR = list(o2.x[:8])
    lvl = {"uD": uD, "uK_lp": list(lp.x[:8]), "uK_lex": uKlex, "XD_lp": ND(list(lp.x[:8]), a) / NDs - 1, "XD_lex": XD_lex,
           "NK_uKlex_rel": NK(uKlex) / NKs - 1, "XK": NK(uD) / NKs - 1, "uR": uR, "robust_ND_rel": ND(uR, a) / NDs - 1,
           "robust_NK_rel": NK(uR) / NKs - 1, "base_ND_rel": ND(BASE8, a) / NDs - 1, "base_NK_rel": NK(BASE8) / NKs - 1,
           "lex_ok": bool(o.success), "robust_ok": bool(o2.success)}
    out["levels"][str(a)] = lvl
    print("E13 Asig %g | XD lp %+.4f lex %+.4f (NK at lex %+.1e) | XK %+.4f | robust ND %+.4f NK %+.4f | base ND %+.4f NK %+.4f | "
          "uKlex %s | uR %s | ok %s %s" % (a, lvl["XD_lp"], XD_lex, lvl["NK_uKlex_rel"], lvl["XK"], lvl["robust_ND_rel"],
                                         lvl["robust_NK_rel"], lvl["base_ND_rel"], lvl["base_NK_rel"],
                                         " ".join("%.3f" % x for x in uKlex), " ".join("%.3f" % x for x in uR), o.success, o2.success), flush=True)
json.dump(out, open(os.path.join(HERE, "h12c_results.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("E13DONE", flush=True)
