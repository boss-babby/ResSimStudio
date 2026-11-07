#!/usr/bin/env python3
"""H10d (Research Cycle 9, exploratory, pre-stated in CYCLE9_PLAN.md): exact decomposition of the convexity gap
   delta = R(D_C) - (R(D_A)+R(D_B))/2 = NL + CV,  NL = R(D_C) - R(Dbar),  CV = R(Dbar) - (R(D_A)+R(D_B))/2 <= 0
for midpoint pairs of 4-period injection schedules (R = Dieterich functional of the per-face Coulomb history, replicated
from --fault-csv in h10b_scp.py).

usage: h10d_decomp.py CLI CASES M SEED > h10d.log
"""
import json, os, random, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
cli, cases, M, seed = sys.argv[1], sys.argv[2].split(","), int(sys.argv[3]), int(sys.argv[4])
sys.argv = [sys.argv[0], cli, "N", "1"]          # h10b_scp reads its own argv at import
sys.path.insert(0, here)
import h10b_scp as H


def random_schedule(rng):          # identical to h10_co2_convexity.random_schedule (not imported: it parses argv at import)
    while True:
        u = [rng.random() for _ in range(4)]
        t = sum(u)
        u = [4 * x / t for x in u]
        if min(u) >= 0.25 and max(u) <= 1.75:
            return u


out = {}
for case in cases:
    rng = random.Random(seed)
    rows = []
    print("CASE %s" % case, flush=True)
    for m in range(M):
        A, B = random_schedule(rng), random_schedule(rng)
        C = [(a + b) / 2 for a, b in zip(A, B)]
        nA, (DA, meta) = H.run(case, A, True)
        nB, (DB, _) = H.run(case, B, True)
        nC, (DC, _) = H.run(case, C, True)
        Db = 0.5 * (DA + DB)
        RA, RB, RC, Rb = (H.dieterich(D, meta) for D in (DA, DB, DC, Db))
        delta = RC - 0.5 * (RA + RB)
        NL, CV = RC - Rb, Rb - 0.5 * (RA + RB)
        s = float(np.sum(DC - Db) / max(1e-30, np.sum(np.abs(DC))))
        rep = max(abs(RA - nA) / nA, abs(RB - nB) / nB, abs(RC - nC) / nC)
        rows.append(dict(A=A, B=B, nA=nA, nB=nB, nC=nC, delta=delta, NL=NL, CV=CV, s=s, replica_err=rep))
        print("PAIR %s %2d delta %+.4e = NL %+.4e + CV %+.4e | NL/|CV| %+.3f | s %+.3e | replica %.1e" % (
            case, m, delta / nC, NL / nC, CV / nC, NL / abs(CV) if CV else float("nan"), s, rep), flush=True)
    nl = np.array([r["NL"] / abs(r["CV"]) for r in rows]); ss = np.array([r["s"] for r in rows])
    med = float(np.median(nl))
    print("SUMMARY case %s median NL/|CV| %+.3f (min %+.3f max %+.3f) | NL > 0 in %d/%d | median s %+.3e | violations %d/%d%s" % (
        case, med, nl.min(), nl.max(), int(np.sum(nl > 0)), M, float(np.median(ss)), sum(r["delta"] > 0 for r in rows), M,
        " | NEAR THE EDGE" if med > 0.5 else ""), flush=True)
    out[case] = rows
    with open(os.path.join(here, os.environ.get("H10D_OUT", "h10d_results.json")), "w") as f:
        json.dump(out, f, indent=1)
print("H10DDONE runs %d" % H.NRUNS, flush=True)
