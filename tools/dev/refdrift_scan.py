#!/usr/bin/env python3
"""C47 developer tool: OPM's own version-to-version drift of opm-tests reference summaries (the noise floor of a
regression case). For each case, the current master reference is compared with the reference at older opm-tests commits
using the regression metric (refdrift = ecl::compareSummaries, relative RMS / max|ref|, + material worst).

usage: refdrift_scan.py REFDRIFT_EXE OUTDIR CASEPATH.SMSPEC ...   (paths relative to the opm-tests root)
"""
import os, re, subprocess, sys, urllib.request
RAW = "https://raw.githubusercontent.com/OPM/opm-tests"
# repository-wide commits (automatic reference updates), newest first
SHAS = {"2026-07": "d07f6a5c5128591cd913066934ce751df67154a6", "2025-05": "da782d48cd", "2024-09": "67db106e68",
        "2023-04": "d5a01f9c5f", "2021-01": "66a7fbe991"}


def full(sha):
    if len(sha) == 40:
        return sha
    with urllib.request.urlopen("https://api.github.com/repos/OPM/opm-tests/commits/" + sha) as r:
        return re.search(r'"sha":\s*"([0-9a-f]{40})"', r.read().decode()).group(1)


def fetch(sha, path, out):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    for ext in ("SMSPEC", "UNSMRY"):
        try:
            with urllib.request.urlopen("%s/%s/%s" % (RAW, sha, path[:-6] + ext)) as r, open(out[:-6] + ext, "wb") as f:
                f.write(r.read())
        except Exception:
            return False
    return True


if __name__ == "__main__":
    exe, outdir, cases = sys.argv[1], sys.argv[2], sys.argv[3:]
    shas = {k: full(v) for k, v in SHAS.items()}
    print("| case | " + " | ".join(shas) + " | max |")
    for p in cases:
        name = os.path.basename(p)[:-7]
        cur = os.path.join(outdir, "master", name + ".SMSPEC")
        if not fetch("master", p, cur):
            print("| %s | missing |" % name); continue
        cells, mx = [], 0.0
        for k, s in shas.items():
            old = os.path.join(outdir, k, name + ".SMSPEC")
            if not fetch(s, p, old):
                cells.append("–"); continue
            o = subprocess.run([exe, cur, old], capture_output=True, text=True).stdout
            m = re.search(r"worst (\S+) ([0-9.]+) %", o)
            mm = re.search(r"material worst (\S+) ([0-9.]+) %", o)
            if not m:
                cells.append("?"); continue
            v = float(m.group(2)); mx = max(mx, v)
            cells.append("%.2f%s" % (v, " [%.2f]" % float(mm.group(2)) if mm else ""))
        print("| %s | %s | %.2f |" % (name, " | ".join(cells), mx), flush=True)
