"""PYACTION-style steering: step through the schedule and act on the state from Python.
Every producer whose water cut exceeds a limit is shut (default deck: the qfs quarter five-spot waterflood, whose
producer passes 60 % water cut after ~4 years; spe1 has no water production, so nothing would be shut there)."""
import sys
import ressim

LIMIT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6

deck = sys.argv[1] if len(sys.argv) > 1 else None
s = ressim.Session(deck) if deck else ressim.Session.example("qfs")
producers = [w["name"] for w in s.info()["wells"] if w["type"] == "producer"]
shut = set()
for step, t in s.steps():
    for w in producers:
        if w in shut:
            continue
        try:
            wct = s.summary("WWCT:" + w)[-1]
        except ressim.ResSimError:
            continue
        if wct > LIMIT:
            msgs = s.apply_schedule(f"WELOPEN\n '{w}' SHUT /\n/\n")
            shut.add(w)
            print(f"t = {t:.0f} d: {w} water cut {wct:.2f} -> shut", msgs or "")
print("FOPT", s.summary("FOPT")[-1], "shut:", sorted(shut))
