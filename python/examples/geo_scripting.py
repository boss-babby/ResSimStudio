"""C43 geomodel scripting example (Petrel-parity P10): demo project -> surface calculator -> polygon clip ->
seismic calculator / crop -> files. Run:  python geo_scripting.py [output_dir]"""
import sys
from ressim.geo import Geomodel

out = sys.argv[1] if len(sys.argv) > 1 else "geo_demo"
with Geomodel() as g:
    g.demo(out)
    print("surfaces:", g.surfaces())
    for w in g.wells():
        print(f"  {w['name']:4s} tops: " + ", ".join(f"{t['name']} {t['z']:.1f}" for t in w["tops"]))
    # isochore between two horizons, limited to a lease polygon
    g.calc("BaseRes - TopRes", onto="TopRes", name="Isochore")
    lease = g.add_polygon("lease", [(600, 500), (2600, 400), (3000, 1800), (1400, 2500), (500, 1600)])
    g.clip("Isochore", lease)
    iso = g.surface("Isochore")
    vals = [z for z in iso["z"] if z == z]
    print(f"isochore inside the lease: {len(vals)} nodes, mean {sum(vals) / len(vals):.2f}")
    g.write_surface("Isochore", f"{out}/isochore_lease.irap")
    g.write_polygons(f"{out}/lease.txt")
    # seismic: envelope below 1700 ms on a decimated sub-volume
    g.seismic_calc("ENVELOPE * (T > 1700)")
    g.crop(i=(5, 35), j=(5, 25), step=(1, 1, 1))
    g.write_segy(f"{out}/envelope_crop.sgy")
    print("seismic:", g.seismic())
