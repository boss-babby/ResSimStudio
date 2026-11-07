#!/usr/bin/env python3
"""H10g calibration (CYCLE12_PLAN.md): u = 1 runs only, fixed rules, no pair information.
usage: h10g_calib.py CLI [col=22]   -> exp/h10g_calib.json
  rule 1: pre-load c = 2 - max CFF0 (sealing fault, no pre-load, no switching); require >= 3 faces CFF > 0 with switching
  rule 2: per fluid, bisection on ramp w for m1 = 1.40 +- 0.10 (LOW) and m1 = 3.50 +- 0.30 (HIGH)
"""
import csv, json, math, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h10_mkdeck import mkdeck

CLI = sys.argv[1]
COL = int(sys.argv[2]) if len(sys.argv) > 2 else 22
WORK = tempfile.mkdtemp(prefix="h10g_", dir="/tmp/sd" if os.path.isdir("/tmp/sd") else None)


def run1(case):
    d, c = os.path.join(WORK, "c.DATA"), os.path.join(WORK, "c.csv")
    open(d, "w").write(mkdeck(case, [1, 1, 1, 1]))
    out = subprocess.run([CLI, d, "--fault-csv", c], capture_output=True, text=True,
                         env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900).stdout
    if "DONE" not in out:
        raise RuntimeError(out[-600:])
    cff, mult = {}, 1.0
    for r in csv.DictReader(open(c)):
        f = int(r["face"]); cff[f] = max(cff.get(f, -1e9), float(r["cff"])); mult = max(mult, float(r["mult"]))
    return max(cff.values()), sum(1 for v in cff.values() if v > 0), mult


def case(col, water, pre, w, pf=50):
    return "X:%d:%d:%g:%g:%g" % (col, water, pre, w, pf)


def bisect(col, water, pre, target, tol):
    lo, hi = 0.01, 2.0                       # m1 decreases with w
    m_lo = run1(case(col, water, pre, lo))[2]
    if m_lo < target - tol:
        return lo, m_lo, "cap"
    for _ in range(30):
        mid = math.sqrt(lo * hi)
        m = run1(case(col, water, pre, round(mid, 5)))[2]
        if abs(m - target) <= tol:
            return round(mid, 5), m, "ok"
        lo, hi = (mid, hi) if m > target else (lo, mid)
    return round(mid, 5), m, "maxiter"


if __name__ == "__main__":
    # sanity: reproduce G1 (F35) multipliers
    chk = {"NSR08 (F35: 2.13)": run1("NSR08")[2], "WNSR08 (F35: 2.92)": run1("WNSR08")[2]}
    print("G1 check", chk, flush=True)
    res = {"col": COL, "g1_check": chk}
    for water, name in ((0, "CO2"), (1, "water")):
        cff0 = run1(case(COL, water, 0, 0.05, 1))[0]
        pre = round(2.0 - cff0, 2)
        while True:
            mx, nre, m = run1(case(COL, water, pre, 0.05))
            if nre >= 3 or pre > 60:
                break
            pre += 1.0
        lw, lm, ls = bisect(COL, water, pre, 1.40, 0.10)
        hw, hm, hs = bisect(COL, water, pre, 3.50, 0.30)
        res[name] = {"cff0": cff0, "preload": pre, "faces_reactivated_w005": nre, "m1_w005": m,
                     "LOW": {"w": lw, "m1": lm, "status": ls, "case": case(COL, water, pre, lw)},
                     "HIGH": {"w": hw, "m1": hm, "status": hs, "case": case(COL, water, pre, hw)}}
        print(name, json.dumps(res[name]), flush=True)
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "h10g_calib.json"), "w"), indent=1)
