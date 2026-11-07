#!/usr/bin/env python3
"""H10 (Research Cycle 9, pre-registered in CYCLE9_PLAN.md): midpoint convexity of the Dieterich expected event count
in the injection-rate history when the fault pressure comes from the full nonlinear ResSim CO2STORE engine.

usage: h10_co2_convexity.py CLI [M=40] [cases=N,F,W,WF] [seed=11]
  1. epsilon: N(u = 1) with report steps of 10 d vs 5 d per case; eps = 3 x |difference| (relative)
  2. M random pairs (A, B) shared by all cases (common random numbers), C = (A + B) / 2,
     u_k in [0.25, 1.75], mean 1; violation if N(C) > (N(A) + N(B)) / 2 + eps * N(C)
Output: lines 'EPS', 'PAIR', 'SUMMARY', 'VERDICT' (h10.log) and h10_results.json
"""
import json, os, random, re, subprocess, sys, tempfile, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h10_mkdeck import mkdeck

CLI = sys.argv[1]
M = int(sys.argv[2]) if len(sys.argv) > 2 else 40
CASES = (sys.argv[3] if len(sys.argv) > 3 else "N,F,W,WF").split(",")
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 11
WORK = tempfile.mkdtemp(prefix="h10_", dir="/tmp/sd" if os.path.isdir("/tmp/sd") else None)
RX = re.compile(r"fault F1: max dCFS ([-0-9.eE+]+), expected events ([-0-9.eE+]+)")


HIST = os.environ.get("H10_HIST", "0") == "1"      # H10c: history-norm nonlinearity of dCFS (needs --fault-csv)
OUTJSON = os.environ.get("H10_OUT", "h10_results.json")


def history(csvpath):
    import csv as _csv
    rows = {}
    with open(csvpath) as f:
        for r in _csv.DictReader(f):
            rows[(float(r["t"]), int(r["face"]))] = float(r["dcfs"])
    os.remove(csvpath)
    return [rows[k] for k in sorted(rows)]


def run(case, u, dt=10.0, csv=None):
    tag = "%s_%s_%g" % (case, "_".join("%.5f" % x for x in u), dt)
    path = os.path.join(WORK, "d%d.DATA" % (abs(hash(tag)) % 10**9))
    with open(path, "w") as f:
        f.write(mkdeck(case, u, dt))
    cmd = [CLI, path] + (["--fault-csv", csv] if csv else [])
    out = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900).stdout
    os.remove(path)
    m = RX.search(out)
    if not m or "DONE" not in out:
        raise RuntimeError("run failed: %s %s\n%s" % (case, u, out[-800:]))
    if csv:
        return float(m.group(2)), history(csv)
    return float(m.group(2)), float(m.group(1))


def random_schedule(rng):
    while True:
        u = [rng.random() for _ in range(4)]
        s = sum(u)
        u = [4 * x / s for x in u]
        if min(u) >= 0.25 and max(u) <= 1.75:
            return u


if __name__ == "__main__":
    rng = random.Random(SEED)
    pairs = [(random_schedule(rng), random_schedule(rng)) for _ in range(M)]
    res = {"M": M, "seed": SEED, "cases": {}}
    for case in CASES:
        n10, d10 = run(case, [1, 1, 1, 1], 10.0)
        n5, _ = run(case, [1, 1, 1, 1], 5.0)
        eps = 3 * abs(n10 - n5) / n10
        print("EPS case %s N(dt 10) %.10g N(dt 5) %.10g eps %.3e (relative)" % (case, n10, n5, eps), flush=True)
        viol, deltas, nl, rows = 0, [], [], []
        for p, (A, B) in enumerate(pairs):
            C = [(a + b) / 2 for a, b in zip(A, B)]
            if HIST:
                hp = os.path.join(WORK, "h.csv")
                (na, ha), (nb, hb), (nc, hc) = run(case, A, csv=hp), run(case, B, csv=hp), run(case, C, csv=hp)
                num = sum((c - (a + b) / 2) ** 2 for a, b, c in zip(ha, hb, hc)) ** 0.5
                den = sum(c * c for c in hc) ** 0.5
                hnl = num / den if den else 0.0
                da, db, dc = max(ha), max(hb), max(hc)
            else:
                (na, da), (nb, db), (nc, dc) = run(case, A), run(case, B), run(case, C)
                hnl = None
            avg = (na + nb) / 2
            delta = (nc - avg) / nc
            v = nc > avg + eps * nc
            viol += v
            deltas.append(delta)
            nl.append(hnl if hnl is not None else (abs(dc - (da + db) / 2) / dc if dc else 0.0))
            rows.append({"A": A, "B": B, "N": [na, nb, nc], "dCFS": [da, db, dc], "delta": delta, "viol": v})
            print("PAIR %s %2d A [%s] B [%s] N %.6g %.6g mid %.6g delta %+.3e%s | maxdCFS %.4g %.4g mid %.4g" % (
                case, p, " ".join("%.2f" % x for x in A), " ".join("%.2f" % x for x in B), na, nb, nc, delta,
                " VIOL" if v else "", da, db, dc), flush=True)
        vd = [d for d, r in zip(deltas, rows) if r["viol"]]
        frac = viol / M
        med_v = statistics.median(vd) if vd else 0.0
        if frac <= 0.05:
            verdict = "SUPPORTED (<= 5 %)"
        elif frac > 0.15 and med_v > 0.01:
            verdict = "REFUTED (> 15 %, median violator delta > 1 %)"
        else:
            verdict = "INCONCLUSIVE"
        print("SUMMARY case %s violations %d/%d (%.0f %%) | median delta all %+.3e | median delta violators %+.3e | max delta %+.3e | "
              "median nonlinearity (%s) %.3e" % (case, viol, M, 100 * frac, statistics.median(deltas), med_v, max(deltas),
                                                     "history norm" if HIST else "max dCFS", statistics.median(nl)), flush=True)
        print("VERDICT case %s %s" % (case, verdict), flush=True)
        res["cases"][case] = {"eps": eps, "N_const": n10, "violations": viol, "frac": frac, "median_delta": statistics.median(deltas),
                              "median_violator_delta": med_v, "max_delta": max(deltas), "median_nonlin_dcfs": statistics.median(nl),
                              "verdict": verdict, "pairs": rows}
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), OUTJSON), "w") as f:
            json.dump(res, f, indent=1)
    print("H10DONE", flush=True)
