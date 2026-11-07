#!/usr/bin/env python3
"""Research Cycle 13 (CYCLE13_PLAN.md, pre-registered): H12a convexity of signed storage schedules, H12b Dieterich vs Kaiser.
usage: h12_storage.py CLI [a|b|ab]  -> exp/h12_results.json (log on stdout)
"""
import csv, json, math, os, random, re, statistics, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h12_mkdeck as M

CLI = sys.argv[1]
WHAT = sys.argv[2] if len(sys.argv) > 2 else "ab"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "h12_results.json")
WORK = tempfile.mkdtemp(prefix="h12_", dir="/tmp/sd" if os.path.isdir("/tmp/sd") else None)
RX = re.compile(r"fault F1: max dCFS ([-0-9.eE+]+), expected events ([-0-9.eE+]+)")
BASE = [1, 1, 1, 1, -1, -1, -1, -1]
res = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    json.dump(res, open(OUT, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


def run(case, u, dt=5.0, asig=1.0, hist=False):
    d = os.path.join(WORK, "d.DATA"); c = os.path.join(WORK, "d.csv")
    open(d, "w").write(M.mkdeck(case, u, dt, asig))
    cmd = [CLI, d] + (["--fault-csv", c] if hist else [])
    o = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900).stdout
    m = RX.search(o)
    if not m or "DONE" not in o:
        raise RuntimeError("run failed %s %s\n%s" % (case, u, o[-800:]))
    if not hist:
        return float(m.group(2)), float(m.group(1))
    H = {}
    for r in csv.DictReader(open(c)):
        H.setdefault(int(r["face"]), []).append((float(r["t"]), float(r["dcfs"]), float(r["area"])))
    return float(m.group(2)), H


def draw(rng):
    while True:
        u = [rng.random() for _ in range(4)]; s = sum(u); u = [4 * x / s for x in u]
        if min(u) >= 0.25 and max(u) <= 1.75:
            return u


def sched(rng):
    return draw(rng) + [-x for x in draw(rng)]


# ---------------------------------------------------------------- calibration (rule G)
def calibrate():
    if "Q_G" in res:
        M.Q_G = res["Q_G"]; print("CAL Q_G (cached) %.6g" % M.Q_G, flush=True); return
    target = run("SW", BASE)[1]
    lo, hi = 360000.0, 600000.0
    for it in range(30):
        mid = 0.5 * (lo + hi); M.Q_G = mid
        d = run("SG", BASE)[1]
        print("CAL it %d Q_G %.6g maxdCFS %.5f target %.5f" % (it, mid, d, target), flush=True)
        if abs(d - target) <= 0.05 * target * 0.2:      # well inside the +-5 % band
            break
        lo, hi = (mid, hi) if d < target else (lo, mid)
    res["Q_G"] = M.Q_G; res["target_dcfs"] = target; save()


# ---------------------------------------------------------------- H12a
def h12a():
    rng = random.Random(31)
    pairs = [(sched(rng), sched(rng)) for _ in range(40)]
    res.setdefault("H12a", {})
    for case in ("SW", "SG"):
        n5, d5 = run(case, BASE, 5.0); n25, _ = run(case, BASE, 2.5)
        eps = 3 * abs(n5 - n25) / n5
        print("EPS %s N5 %.10g N2.5 %.10g eps %.3e maxdCFS %.4f" % (case, n5, n25, eps, d5), flush=True)
        viol, deltas, rows = 0, [], []
        for p, (A, B) in enumerate(pairs):
            C = [(a + b) / 2 for a, b in zip(A, B)]
            (na, da), (nb, db), (nc, dc) = run(case, A), run(case, B), run(case, C)
            avg = (na + nb) / 2; delta = (nc - avg) / nc; v = nc > avg + eps * nc
            viol += v; deltas.append(delta)
            rows.append({"A": A, "B": B, "N": [na, nb, nc], "dCFS": [da, db, dc], "delta": delta, "viol": bool(v)})
            print("PAIR %s %2d N %.6g %.6g mid %.6g delta %+.3e%s" % (case, p, na, nb, nc, delta, " VIOL" if v else ""), flush=True)
        vd = [r["delta"] for r in rows if r["viol"]]
        med_v = statistics.median(vd) if vd else 0.0
        verdict = ("SUPPORTED" if viol <= 2 else "REFUTED" if (viol > 6 and med_v > 0.01) else "INCONCLUSIVE")
        print("SUMMARY %s violations %d/40 median delta %+.3e median violator %+.3e max %+.3e" % (
            case, viol, statistics.median(deltas), med_v, max(deltas)), flush=True)
        print("VERDICT H12a %s %s" % (case, verdict), flush=True)
        res["H12a"][case] = {"eps": eps, "N_base": n5, "maxdcfs_base": d5, "violations": viol, "median_delta": statistics.median(deltas),
                             "median_violator_delta": med_v, "max_delta": max(deltas), "verdict": verdict, "pairs": rows}
        save()


# ---------------------------------------------------------------- H12b
def dieterich(H, asig, ta=100.0, r0=5.0):
    """Python copy of faultmech.cpp: constant stressing rate per report step, 16 Simpson sub-steps on R = r/(gamma Sr)."""
    tot = 0.0
    Sr = asig / ta
    for f, ser in H.items():
        g = ta / asig
        r = r0 * ser[0][2] * 1e-6
        tp, dp = 0.0, 0.0
        for (t, d, _) in ser:
            dtYr = (t - tp) / 365.25
            if dtYr <= 0:
                tp, dp = t, d; continue
            Sd = Sr + (d - dp) / dtYr
            h = dtYr / 16

            def adv(gg, hh):
                if abs(Sd) * hh / asig < 1e-9:
                    return gg + hh / asig * (1 - gg * Sd)
                e = math.exp(max(-700, min(700, -Sd * hh / asig)))
                return (gg - 1 / Sd) * e + 1 / Sd
            for _ in range(16):
                gm, g1 = max(1e-300, adv(g, 0.5 * h)), max(1e-300, adv(g, h))
                tot += h / 6 * (r / (max(1e-300, g) * Sr) + 4 * r / (gm * Sr) + r / (g1 * Sr))
                g = g1
            tp, dp = t, d
    return tot


def h12b():
    import numpy as np
    from scipy.optimize import minimize, linprog
    G = []
    for k in range(8):
        _, H = run("SW", "pulse:%d" % k, hist=True)
        G.append(H)
    _, H0 = run("SW", [0.0] * 8, hist=True)                     # amendment A1: affine response (zero-rate drift)
    faces = sorted(G[0]); times = [t for (t, _, _) in G[0][faces[0]]]; area = G[0][faces[0]][0][2]
    G0 = np.array([[H0[f][j][1] for j in range(len(times))] for f in faces])                        # f, t
    Gm = np.array([[[G[k][f][j][1] for j in range(len(times))] for f in faces] for k in range(8)]) - G0[None]   # k, f, t

    def H_of(u):
        D = G0 + np.tensordot(np.asarray(u, float), Gm, axes=1)                                     # f, t
        return {f: [(times[j], D[i, j], area) for j in range(len(times))] for i, f in enumerate(faces)}, D

    def NK(u):
        D = H_of(u)[1]
        return float(np.sum(np.maximum(0.0, D.max(axis=1))))

    # linearity / implementation check against the engine (5 random H12a schedules)
    rng = random.Random(31); errs = []
    for _ in range(5):
        u = sched(rng)
        ne = run("SW", u)[0]; npy = dieterich(H_of(u)[0], 1.0)
        errs.append(abs(npy - ne) / ne)
        print("H12b CHECK engine %.8g matrix %.8g rel %.3e" % (ne, npy, errs[-1]), flush=True)
    med_err = statistics.median(errs)
    valid = med_err <= 0.02
    print("H12b LINEARITY median rel err %.3e -> %s" % (med_err, "VALID" if valid else "NON-DISCRIMINATING"), flush=True)
    out = {"linearity_median_err": med_err, "valid": valid, "levels": {}}
    # Kaiser optimum: LP in (u, s_f); min sum s_f, s_f >= D_f(t) = sum_k u_k G_kft, s_f >= 0
    nf, nt = len(faces), len(times)
    c = np.r_[np.zeros(8), np.ones(nf)]
    Aub, bub = [], []
    for i in range(nf):
        for j in range(nt):
            row = np.zeros(8 + nf); row[:8] = Gm[:, i, j]; row[8 + i] = -1; Aub.append(row); bub.append(-G0[i, j])
    Aeq = np.zeros((2, 8 + nf)); Aeq[0, :4] = 1; Aeq[1, 4:8] = 1; beq = [4.0, -4.0]
    bounds = [(0, 2)] * 4 + [(-2, 0)] * 4 + [(0, None)] * nf
    lp = linprog(c, A_ub=np.array(Aub), b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method="highs")
    uK = list(lp.x[:8])
    print("H12b uK %s NK %.6g (LP status %d)" % (" ".join("%.3f" % x for x in uK), NK(uK), lp.status), flush=True)
    for asig in (1.0, 10.0):
        f = lambda u: dieterich(H_of(u)[0], asig)
        cons = [{"type": "eq", "fun": lambda u: sum(u[:4]) - 4}, {"type": "eq", "fun": lambda u: sum(u[4:]) + 4}]
        best = None; r2 = random.Random(7)
        starts = [BASE, uK] + [sched(r2) for _ in range(6)]
        for s in starts:
            sc = f(BASE)
            o = minimize(lambda u: f(u) / sc, s, method="SLSQP", bounds=[(0, 2)] * 4 + [(-2, 0)] * 4, constraints=cons,
                         options={"maxiter": 200, "ftol": 1e-12})
            if best is None or o.fun < best.fun:
                best = o
        uD = list(best.x)
        L1 = sum(abs(a - b) for a, b in zip(uD, uK)) / sum(abs(b) for b in uK)
        XD = f(uK) / f(uD) - 1; XK = NK(uD) / NK(uK) - 1
        Pi = float(np.max(H_of(BASE)[1])) / asig
        lvl = {"Pi": Pi, "uD": uD, "uK": uK, "ND_uD": f(uD), "ND_uK": f(uK), "ND_base": f(BASE), "NK_uD": NK(uD), "NK_uK": NK(uK),
               "NK_base": NK(BASE), "L1": L1, "XD": XD, "XK": XK}
        if asig == 1.0:
            ok = XD < 0.05 and XK < 0.05; lvl["verdict"] = "H12b-1 " + ("SUPPORTED" if ok else "REFUTED")
        else:
            ok = L1 > 0.10 and max(XD, XK) > 0.05; lvl["verdict"] = "H12b-2 " + ("SUPPORTED" if ok else "REFUTED")
        if not valid:
            lvl["verdict"] += " (NOT INTERPRETED: linearity check failed)"
        print("H12b Asig %g Pi %.3f uD %s | ND base %.6g uD %.6g uK %.6g | NK base %.6g uD %.6g uK %.6g | L1 %.3f XD %+.3f XK %+.3f | %s" % (
            asig, Pi, " ".join("%.3f" % x for x in uD), lvl["ND_base"], lvl["ND_uD"], lvl["ND_uK"], lvl["NK_base"], lvl["NK_uD"],
            lvl["NK_uK"], L1, XD, XK, lvl["verdict"]), flush=True)
        out["levels"][str(asig)] = lvl
    # engine verification of the Aσ = 1 designs
    l1 = out["levels"]["1.0"]
    eD, eK, eB = run("SW", l1["uD"])[0], run("SW", l1["uK"])[0], run("SW", BASE)[0]
    out["engine_check"] = {"N_uD": eD, "N_uK": eK, "N_base": eB, "ordering_ok": eD <= eK}
    print("H12b ENGINE N(base) %.6g N(uD) %.6g N(uK) %.6g ordering %s" % (eB, eD, eK, "OK" if eD <= eK else "FAILED"), flush=True)
    res["H12b"] = out; save()


calibrate()
if "a" in WHAT:
    h12a()
if "b" in WHAT:
    h12b()
print("H12DONE", flush=True)
