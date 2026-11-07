"""Permeability multiplier sweep: edit the deck, run, collect FOPT (any deck with PERMX)."""
import sys
import ressim

deck = sys.argv[1] if len(sys.argv) > 1 else "example:spe1"
base = ressim.Session(deck)
nx, ny, nz = base.dims
kx = base.static("PERMX")
for mult in (0.5, 1.0, 2.0):
    s = ressim.Session(deck)
    vals = " ".join(f"{v * mult:.6g}" for v in kx)
    s.set("PERMX", vals)
    s.run()
    print(f"PERMX x {mult}: FOPT = {s.summary('FOPT')[-1]:.5g}")
