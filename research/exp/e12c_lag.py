#!/usr/bin/env python3
"""E12c (exploratory, after the H10g verdict): is the two-way convexity loss an artefact of the one-report-step coupling lag?
Re-run every violating H10g HIGH pair with report steps 5 d and 2.5 d (lag halved / quartered) and recompute delta.
usage: e12c_lag.py CLI   -> e12c.log, e12c_results.json"""
import json, os, re, subprocess, sys, tempfile, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h10_mkdeck import mkdeck
CLI = sys.argv[1]; WORK = tempfile.mkdtemp(prefix="e12c_", dir="/tmp/sd")
RX = re.compile(r"fault F1: max dCFS ([-0-9.eE+]+), expected events ([-0-9.eE+]+)")
def N(case, u, dt):
    p = os.path.join(WORK, "d.DATA"); open(p, "w").write(mkdeck(case, u, dt))
    out = subprocess.run([CLI, p], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900).stdout
    m = RX.search(out)
    if not m or "DONE" not in out: raise RuntimeError(out[-500:])
    return float(m.group(2))
R = json.load(open("h10g_results.json")); res = {}
for case in ("X:22:0:32.68:0.03761:50", "X:22:1:27.06:0.04438:50"):
    c = R["cases"][case]; rows = []
    for p in c["pairs"]:
        if not p["viol"]: continue
        A, B = p["A"], p["B"]; C = [(a + b) / 2 for a, b in zip(A, B)]; r = {"delta10": p["delta"]}
        for dt in (5.0, 2.5):
            na, nb, nc = N(case, A, dt), N(case, B, dt), N(case, C, dt); r["delta%g" % dt] = (nc - (na + nb) / 2) / nc
        rows.append(r); print(case, " ".join("%s %+.4f" % kv for kv in r.items()), flush=True)
    keep = lambda k: sum(1 for r in rows if r[k] > c["eps"])
    res[case] = {"eps": c["eps"], "n_viol_dt10": len(rows), "still_viol_dt5": keep("delta5"), "still_viol_dt2.5": keep("delta2.5"),
                 "median_delta": {k: statistics.median(r[k] for r in rows) for k in ("delta10", "delta5", "delta2.5")}, "rows": rows}
    print("SUMMARY", case, json.dumps({k: v for k, v in res[case].items() if k != "rows"}), flush=True)
json.dump(res, open("e12c_results.json", "w"), indent=1)
