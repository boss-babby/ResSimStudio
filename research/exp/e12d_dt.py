#!/usr/bin/env python3
"""E12d (exploratory, after H10g): unbiased lag check. All 40 H10g pairs of a case at report step DT, with its own epsilon
eps = 3 |N(u=1, DT) - N(u=1, DT/2)| / N(DT).  usage: e12d_dt.py CLI CASE DT [subset=all|viol] OUT.json"""
import json, os, re, subprocess, sys, tempfile, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h10_mkdeck import mkdeck
CLI, CASE, DT = sys.argv[1], sys.argv[2], float(sys.argv[3]); SUB = sys.argv[4]; OUT = sys.argv[5]
WORK = tempfile.mkdtemp(prefix="e12d_", dir="/tmp/sd"); RX = re.compile(r"expected events ([-0-9.eE+]+)")
def N(u, dt):
    p = os.path.join(WORK, "d.DATA"); open(p, "w").write(mkdeck(CASE, u, dt))
    out = subprocess.run([CLI, p], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=1800).stdout
    m = RX.search(out)
    if not m or "DONE" not in out: raise RuntimeError(out[-500:])
    return float(m.group(1))
n1, n2 = N([1, 1, 1, 1], DT), N([1, 1, 1, 1], DT / 2); eps = 3 * abs(n1 - n2) / n1
print("EPS", CASE, DT, n1, n2, eps, flush=True)
pairs = json.load(open("h10g_results.json"))["cases"][CASE]["pairs"]; rows = []
for i, p in enumerate(pairs):
    if SUB == "viol" and not p["viol"]: continue
    A, B = p["A"], p["B"]; C = [(a + b) / 2 for a, b in zip(A, B)]
    na, nb, nc = N(A, DT), N(B, DT), N(C, DT); d = (nc - (na + nb) / 2) / nc
    rows.append({"i": i, "delta10": p["delta"], "delta": d, "viol": d > eps}); print("PAIR", i, "%+.4f %+.4f" % (p["delta"], d), d > eps, flush=True)
v = [r for r in rows if r["viol"]]
s = {"case": CASE, "dt": DT, "subset": SUB, "eps": eps, "n": len(rows), "violations": len(v),
     "median_delta_violators": statistics.median([r["delta"] for r in v]) if v else 0.0, "median_delta_all": statistics.median([r["delta"] for r in rows])}
print("SUMMARY", json.dumps(s), flush=True); s["rows"] = rows; json.dump(s, open(OUT, "w"), indent=1)
