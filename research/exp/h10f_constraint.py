#!/usr/bin/env python3
"""H10f-b (Research Cycle 11, pre-registered in CYCLE11_PLAN.md): reactivation constraint.

usage: h10f_constraint.py CLI [M=40] [cases=NSB,WNSB] [seed=13]
A design is feasible when no fault face reaches CFF >= 0 at any report time (--fault-csv column cff), so it never switches
the fault permeability. For the H10e pairs (seed 13) with A and B feasible: (1) is C = (A + B)/2 feasible (set convexity)?
(2) for A, B, C feasible: midpoint-convexity violation of the expected event count as in H10e (eps from 10 d / 5 d).
Output: h10f_b.log lines PRE / EPS / PAIR / SUMMARY / VERDICT and h10f_b_results.json.
Helpers are copied from h10_co2_convexity.py (that script parses argv at import).
"""
import csv, json, os, random, re, subprocess, sys, tempfile, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h10_mkdeck import mkdeck

CLI = sys.argv[1]
M = int(sys.argv[2]) if len(sys.argv) > 2 else 40
CASES = (sys.argv[3] if len(sys.argv) > 3 else "NSB,WNSB").split(",")
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 13
WORK = tempfile.mkdtemp(prefix="h10f_", dir="/tmp/sd" if os.path.isdir("/tmp/sd") else None)
RX = re.compile(r"fault F1: max dCFS ([-0-9.eE+]+), expected events ([-0-9.eE+]+)")


def run(case, u, dt=10.0):
    tag = "%s_%s_%g" % (case, "_".join("%.5f" % x for x in u), dt)
    path = os.path.join(WORK, "d%d.DATA" % (abs(hash(tag)) % 10**9))
    fc = path + ".csv"
    with open(path, "w") as f:
        f.write(mkdeck(case, u, dt))
    out = subprocess.run([CLI, path, "--fault-csv", fc], capture_output=True, text=True,
                         env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900).stdout
    os.remove(path)
    m = RX.search(out)
    if not m or "DONE" not in out:
        raise RuntimeError("run failed: %s %s\n%s" % (case, u, out[-800:]))
    cmax = -1e300
    with open(fc) as f:
        for r in csv.DictReader(f):
            c = float(r["cff"])
            if c == c:
                cmax = max(cmax, c)
    os.remove(fc)
    return float(m.group(2)), cmax


def random_schedule(rng):
    while True:
        u = [rng.random() for _ in range(4)]
        s = sum(u)
        u = [4 * x / s for x in u]
        if min(u) >= 0.25 and max(u) <= 1.75:
            return u


def verdict_obj(frac, med_v, n):
    if n == 0:
        return "NON-DISCRIMINATING (no pair with A, B, C feasible)"
    if frac <= 0.05:
        return "SUPPORTED (<= 5 %)"
    if frac > 0.15 and med_v > 0.01:
        return "REFUTED (> 15 %, median violator delta > 1 %)"
    return "INCONCLUSIVE"


if __name__ == "__main__":
    rng = random.Random(SEED)
    pairs = [(random_schedule(rng), random_schedule(rng)) for _ in range(M)]
    res = {"M": M, "seed": SEED, "cases": {}}
    for case in CASES:
        n10, c1 = run(case, [1, 1, 1, 1], 10.0)
        n5, _ = run(case, [1, 1, 1, 1], 5.0)
        eps = 3 * abs(n10 - n5) / n10
        print("PRE case %s u=1 max CFF %.4g bar (%s)" % (case, c1, "below failure" if c1 < 0 else "REACTIVATES - non-discriminating"), flush=True)
        print("EPS case %s N(dt 10) %.10g N(dt 5) %.10g eps %.3e (relative)" % (case, n10, n5, eps), flush=True)
        rows, designs = [], []
        for p, (A, B) in enumerate(pairs):
            na, ca = run(case, A)
            nb, cb = run(case, B)
            designs += [(na, ca < 0), (nb, cb < 0)]
            row = {"A": A, "B": B, "N": [na, nb, None], "cff": [ca, cb, None], "feasAB": ca < 0 and cb < 0}
            if row["feasAB"]:
                C = [(a + b) / 2 for a, b in zip(A, B)]
                nc, cc = run(case, C)
                row["N"][2], row["cff"][2] = nc, cc
                row["feasC"] = cc < 0
                avg = (na + nb) / 2
                row["delta"] = (nc - avg) / nc
                row["viol"] = bool(cc < 0 and nc > avg + eps * nc)
            rows.append(row)
            print("PAIR %s %2d A [%s] B [%s] maxCFF %.3g %.3g%s N %.6g %.6g%s" % (
                case, p, " ".join("%.2f" % x for x in A), " ".join("%.2f" % x for x in B), ca, cb,
                (" mid %.3g %s" % (row["cff"][2], "feasible" if row["feasC"] else "INFEASIBLE")) if row["feasAB"] else " (not both feasible)",
                na, nb, (" mid %.6g delta %+.3e%s" % (row["N"][2], row["delta"], " VIOL" if row["viol"] else "")) if row["feasAB"] else ""), flush=True)
        ab = [r for r in rows if r["feasAB"]]
        abc = [r for r in ab if r["feasC"]]
        nab = len(ab)
        set_frac = (len(abc) / nab) if nab else float("nan")
        if nab < 10:
            set_v = "NON-DISCRIMINATING (%d < 10 pairs with A and B feasible)" % nab
        elif set_frac >= 0.95:
            set_v = "SUPPORTED (>= 95 % of midpoints feasible)"
        elif set_frac < 0.85:
            set_v = "REFUTED (< 85 % of midpoints feasible)"
        else:
            set_v = "INCONCLUSIVE"
        viol = sum(r["viol"] for r in abc)
        vfrac = viol / len(abc) if abc else float("nan")
        vd = [r["delta"] for r in abc if r["viol"]]
        med_v = statistics.median(vd) if vd else 0.0
        obj_v = verdict_obj(vfrac if abc else 0.0, med_v, len(abc))
        feas = [n for n, f in designs if f]
        allN = [n for n, _ in designs]
        print("SUMMARY case %s pairs A,B feasible %d/%d | midpoints feasible %d/%d (%.0f %%) | objective violations %d/%d | "
              "median delta %+.3e | feasible designs %d/%d | best feasible N %.6g vs best N %.6g" % (
                  case, nab, M, len(abc), nab, 100 * set_frac if nab else float("nan"), viol, len(abc),
                  statistics.median([r["delta"] for r in abc]) if abc else float("nan"), len(feas), len(allN),
                  min(feas) if feas else float("nan"), min(allN)), flush=True)
        print("VERDICT case %s set %s | objective %s" % (case, set_v, obj_v), flush=True)
        res["cases"][case] = {"eps": eps, "N_const": n10, "cff_const": c1, "pairs_AB_feasible": nab, "midpoints_feasible": len(abc),
                              "set_verdict": set_v, "violations": viol, "objective_verdict": obj_v,
                              "feasible_designs": len(feas), "designs": len(allN),
                              "best_feasible_N": min(feas) if feas else None, "best_N": min(allN), "pairs": rows}
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "h10f_b_results.json"), "w") as f:
            json.dump(res, f, indent=1)
    print("H10FDONE", flush=True)
