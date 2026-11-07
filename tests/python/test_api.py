"""C33 Python API checks (run: python tests/python/test_api.py; needs RESSIM_LIB or a built library)."""
import math
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import ressim  # noqa: E402

fails = 0


def check(cond, what):
    global fails
    print(("  [PASS] " if cond else "  [FAIL] ") + what)
    fails += 0 if cond else 1


tmp = tempfile.mkdtemp(prefix="ressim_py_")
print("library:", ressim.version())
check("qfs" in ressim.examples(), "example list")

base = ressim.Session.example("qfs")
info = base.info()
check(info["dims"] == [25, 25, 1] and len(info["wells"]) == 2, "info: dims and wells")
r = base.run()
check(r["ok"] and r["timesteps"] > 0, f"whole run ok ({r['timesteps']} steps, {r['seconds']:.2f} s)")
t = base.summary("TIME")
fopt = base.summary("FOPT")
check(len(t) == len(fopt) == base.summary_length if hasattr(base, "summary_length") else len(t) == len(fopt), "summary vector lengths")
check(abs(fopt[-1] - sum(a * b for a, b in zip(base.summary("FOPR"), [t[0]] + [t[i] - t[i - 1] for i in range(1, len(t))]))) < 1e-3 * fopt[-1], "FOPT = sum FOPR dt")
prod = [w["name"] for w in info["wells"] if w["type"] == "producer"][0]
wopt = base.summary("WOPT:" + prod)
check(abs(wopt[-1] - fopt[-1]) < 1e-4 * fopt[-1], "WOPT of the only producer = FOPT")
p0, p1 = base.cell("PRESSURE", 0), base.cell("PRESSURE", -1)
check(len(p0) == 625 and len(p1) == 625, "cell arrays per active cell")
full = base.cell("SWAT", -1, full=True)
check(len(full) == 625 and not any(math.isnan(x) for x in full), "full-grid array")
check(base.num_snapshots == len(base.snapshot_times()) and base.num_snapshots >= 2, "snapshots")

# step-wise run without steering = whole run
s = ressim.Session.example("qfs")
n = 0
for step, tt in s.steps():
    n += 1
check(abs(s.summary("FOPT")[-1] - fopt[-1]) <= 1e-9 * fopt[-1], f"step-wise run = whole run ({n} report steps)")

# steering: shut the producer half way -> less oil, no production afterwards
s = ressim.Session.example("qfs")
s.begin()
half = s.num_report_steps // 2
while s.report_step < half and s.step():
    pass
msgs = s.apply_schedule(f"WELOPEN\n '{prod}' SHUT /\n/\n")
t_shut = s.time
while s.step():
    pass
f2 = s.summary("FOPT")
check(f2[-1] < fopt[-1] * 0.999, f"steering: shut at {t_shut:.0f} d reduces FOPT ({f2[-1]:.4g} < {fopt[-1]:.4g})")
fopr = s.summary("FOPR")
tt = s.summary("TIME")
check(all(q == 0 for q, x in zip(fopr, tt) if x > t_shut + 1e-6), "no oil after the shut-in")

# deck edit
s = ressim.Session.example("qfs")
kx = s.static("PERMX")
s.set("PERMX", " ".join(f"{v * 2:.6g}" for v in kx))
check(abs(s.static("PERMX")[0] - 2 * kx[0]) < 1e-9 * kx[0], "set PERMX re-imports the deck")

# reports and ECLIPSE output
for ext in ("docx", "pptx", "html", "md"):
    path = os.path.join(tmp, "rep." + ext)
    base.report(path, author="test", maxwells=2)
    check(os.path.getsize(path) > 2000, f"report .{ext} ({os.path.getsize(path)} bytes)")
basep = os.path.join(tmp, "QFS")
base.write_eclipse(basep)
ref = ressim.open_summary(basep + ".SMSPEC")
d = {x["key"]: x for x in ressim.compare(ref, base, ["FOPT", "FPR", "FWCT"])}
check(set(d) == {"FOPT", "FPR", "FWCT"} and max(x["rel_rms"] for x in d.values()) < 1e-6, "write_eclipse -> open_summary -> compare (identical)")
try:
    base.summary("NOPE")
    check(False, "error on an unknown vector")
except ressim.ResSimError as e:
    check("NOPE" in str(e), "error on an unknown vector")
print("ALL PASS" if fails == 0 else f"{fails} FAILED")
sys.exit(1 if fails else 0)
