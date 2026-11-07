"""Run a built-in example, print key results and write a Word report."""
import ressim

print(ressim.version(), "-", ressim.library_path if ressim.library_path else "")
s = ressim.Session.example("spe1")
print(s.info()["title"], s.dims)
info = s.run()
print("steps", info["timesteps"], "newton", info["newton"], "seconds", round(info["seconds"], 2))
t, fopt = s.summary("TIME"), s.summary("FOPT")
print(f"FOPT after {t[-1]:.0f} days: {fopt[-1]:.4g}")
p = s.cell("PRESSURE", step=-1)
print("final pressure range", min(p), max(p))
s.report("quickstart_report.docx", author="ResSim Python API")
print("wrote quickstart_report.docx")
