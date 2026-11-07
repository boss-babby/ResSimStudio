#!/usr/bin/env python3
"""H10b (Research Cycle 9, pre-registered in CYCLE9_PLAN.md): sequential convex programming (SCP) driven by the simulator
for minimising the Dieterich expected event count of a fault during CO2 injection (ResSim CO2STORE, FAULTMECH).

  At iterate u_j: one run with the per-face fault history (--fault-csv) + K = 4 pulse runs (u_j + h e_k) give the
  linearised Coulomb history dCFS(u) ~ dCFS(u_j) + J (u - u_j). The Dieterich count of the linearised history is convex
  in u (F1: convex in the stress history; linear map), and is minimised on {sum u = 4, 0.25 <= u <= 1.75,
  |u - u_j|_inf <= Delta} (SLSQP, several starts). The step is accepted on the TRUE simulated count (trust-region ratio).
  Reference: exhaustive grid search, u_k in {0.25, 0.5, ..., 1.75}, sum 4.

usage: h10b_scp.py CLI [cases=N,F] [iters=6]
"""
import csv, itertools, json, math, os, subprocess, sys, tempfile
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h10_mkdeck import mkdeck

CLI = sys.argv[1]
CASES = (sys.argv[2] if len(sys.argv) > 2 else "N,F").split(",")
ITERS = int(sys.argv[3]) if len(sys.argv) > 3 else 6
WORK = tempfile.mkdtemp(prefix="h10b_", dir="/tmp/sd" if os.path.isdir("/tmp/sd") else None)
TA, R0 = 100.0, 5.0          # FAULTMECH ta (years), r0 (events / km2 / yr) of h10_mkdeck
LO, HI, H = 0.25, 1.75, 0.1
NRUNS = 0


def run(case, u, history=False):
    global NRUNS
    NRUNS += 1
    path = os.path.join(WORK, "d.DATA")
    with open(path, "w") as f:
        f.write(mkdeck(case, list(u)))
    cmd = [CLI, path]
    fc = os.path.join(WORK, "f.csv")
    if history:
        cmd += ["--fault-csv", fc]
    out = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900).stdout
    i = out.find("expected events ")
    if i < 0 or "DONE" not in out:
        raise RuntimeError(out[-600:])
    n = float(out[i + 16:].split(",")[0])
    if not history:
        return n, None
    t, faces = [], {}
    with open(fc) as f:
        for r in csv.DictReader(f):
            tt = float(r["t"])
            if not t or t[-1] != tt:
                t.append(tt)
            faces.setdefault(int(r["face"]), {"area": float(r["area"]), "aSig": float(r["aSig"]), "d": []})["d"].append(float(r["dcfs"]))
    ids = sorted(faces)
    D = np.array([faces[i]["d"] for i in ids])            # faces x steps
    meta = (np.array([faces[i]["area"] for i in ids]), np.array([faces[i]["aSig"] for i in ids]), np.array(t))
    return n, (D, meta)


def dieterich(D, meta):
    """Engine algorithm (faultmech.cpp): per face, constant stressing rate per report step, 16 sub-steps, Simpson."""
    area, asig, t = meta
    Sr = asig / TA
    r = R0 * area * 1e-6
    g = TA / asig
    tot = 0.0
    for j in range(1, len(t)):
        dt = (t[j] - t[j - 1]) / 365.25
        if dt <= 0:
            continue
        Sdot = Sr + (D[:, j] - D[:, j - 1]) / dt
        h = dt / 16

        def adv(gg, hh):
            x = Sdot * hh / asig
            e = np.exp(np.clip(-x, -700, 700))
            small = np.abs(x) < 1e-9
            return np.where(small, gg + hh / asig * (1 - gg * Sdot), (gg - 1 / Sdot) * e + 1 / Sdot)
        for _ in range(16):
            gm = np.maximum(1e-300, adv(g, 0.5 * h)); g1 = np.maximum(1e-300, adv(g, h))
            tot += float(np.sum(h / 6 * (r / (np.maximum(1e-300, g) * Sr) + 4 * r / (gm * Sr) + r / (g1 * Sr))))
            g = g1
    return tot


def grid_reference(case):
    vals = [LO + 0.25 * i for i in range(7)]
    best, bu, n = math.inf, None, 0
    for u in itertools.product(vals, repeat=4):
        if abs(sum(u) - 4) > 1e-9:
            continue
        v, _ = run(case, u); n += 1
        if v < best:
            best, bu = v, u
    return best, bu, n


def scp(case, u0, log):
    u = np.array(u0, float)
    Delta = 0.5
    Nj, (D0, meta) = run(case, u, True)
    log("  start u %s N %.8g (replicated Dieterich %.8g)" % (np.round(u, 3).tolist(), Nj, dieterich(D0, meta)))
    hist = [(0, u.tolist(), Nj, 5 * 0 + 1)]
    runs = 1
    for it in range(1, ITERS + 1):
        J = []
        for k in range(4):
            up = u.copy(); up[k] += H
            _, (Dk, _) = run(case, up, True); runs += 1
            J.append((Dk - D0) / H)
        J = np.array(J)                                           # K x faces x steps
        model = lambda v: dieterich(D0 + np.tensordot(v - u, J, axes=1), meta)
        m0 = model(u)
        lo = np.maximum(LO, u - Delta); hi = np.minimum(HI, u + Delta)
        cons = [{"type": "eq", "fun": lambda v: np.sum(v) - 4.0}]
        best = None
        starts = [u, np.clip(np.full(4, 1.0), lo, hi)] + [np.clip(u + Delta * s, lo, hi) for s in (np.array([1, -1, 1, -1]) * 0.5, np.array([-1, 1, -1, 1]) * 0.5, np.array([-1, -1, 1, 1]) * 0.5, np.array([1, 1, -1, -1]) * 0.5)]
        for s0 in starts:
            r = minimize(model, s0, method="SLSQP", bounds=list(zip(lo, hi)), constraints=cons, options={"maxiter": 200, "ftol": 1e-12})
            if r.success or r.status in (8, 9):
                v = np.clip(r.x, lo, hi)
                if abs(v.sum() - 4) < 1e-6 and (best is None or r.fun < best[1]):
                    best = (v, r.fun)
        if best is None:
            log("  iter %d: subproblem failed" % it); break
        v, mv = best
        Nv, hv = run(case, v, True); runs += 1
        pred, act = m0 - mv, Nj - Nv
        rho = act / pred if pred > 1e-15 else 0.0
        acc = rho > 0.1 and act > 0
        log("  iter %d Delta %.3f model %.6g -> %.6g | true %.8g -> %.8g | rho %.3f %s | u %s" % (
            it, Delta, m0, mv, Nj, Nv, rho, "accept" if acc else "reject", np.round(v, 3).tolist()))
        if acc:
            u, Nj, (D0, meta) = v, Nv, hv
        Delta = min(0.75, 2 * Delta) if rho > 0.75 else (0.5 * Delta if rho < 0.25 else Delta)
        hist.append((it, u.tolist(), Nj, runs))
        if Delta < 1e-3:
            break
    return u, Nj, hist, runs


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    out = {}
    for case in CASES:
        log = lambda s: print(s, flush=True)
        log("CASE %s" % case)
        g, gu, ng = grid_reference(case)
        log("GRID case %s best %.8g at %s (%d feasible grid points)" % (case, g, gu, ng))
        res = {"grid_best": g, "grid_u": gu, "grid_n": ng}
        for name, u0 in (("primary", [1.75, 1.75, 0.25, 0.25]), ("constant", [1, 1, 1, 1])):
            log(" SCP start %s" % name)
            u, N, hist, runs = scp(case, u0, log)
            gap = (N - g) / g
            within6 = None
            for (it, uu, nn, rr) in hist:
                if it <= 6 and (nn - g) / g <= 0.01:
                    within6 = (it, rr); break
            log(" RESULT case %s start %s N %.8g gap to grid %+.3f %% | first <= 1 %%: %s | runs %d" % (case, name, N, 100 * gap, within6, runs))
            res[name] = {"u": list(map(float, u)), "N": N, "gap": gap, "first_within_1pct": within6, "runs": runs, "hist": hist}
        p = res["primary"]
        v = "SUPPORTED" if p["first_within_1pct"] else ("REFUTED" if p["gap"] > 0.05 else "INCONCLUSIVE")
        res["verdict"] = v
        log("VERDICT case %s %s" % (case, v))
        out[case] = res
        with open(os.path.join(here, "h10b_results.json"), "w") as f:
            json.dump(out, f, indent=1)
    print("H10BDONE total runs %d" % NRUNS, flush=True)
