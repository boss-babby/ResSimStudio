# ResSim Python API (C33)

```
pip install -e python            # or just put python/ on PYTHONPATH
set RESSIM_LIB=C:\ResSim\ressim.dll   (only needed when the library is not next to the package / in release/)
```

The package is a ctypes wrapper of the engine C API (`engine/include/ressim/capi.h`), so it works with any Python 3.8+
without a compiler. numpy / pandas are used when installed.

| Task | Code |
|---|---|
| open | `s = ressim.Session("CASE.DATA")`, `ressim.Session.example("spe1")`, `ressim.Session(text=deck)` |
| edit | `s.set("PERMX", "300*200")`, `s.set("RESSIM:FAULTSLIP", "20 0.25 0.1")` |
| run | `s.run()` (any engine) or `s.begin(); while s.step(): ...` (black oil, step-wise) |
| steer | `s.apply_schedule("WCONPROD\n 'P1' OPEN ORAT 500 4* 1000 /\n/\n")` - any SCHEDULE keywords, like an ACTIONX body |
| summary | `s.summary_keys()`, `s.summary("FOPR")`, `s.summary("WBHP:P1")`, `s.to_dataframe()` |
| cells | `s.cell("PRESSURE", step=-1)`, `s.cell("SWAT", 0, full=True)`, `s.static("PERMX")`, `s.active_index()` |
| output | `s.write_eclipse("OUT/CASE")`, `s.export_deck("x.DATA")`, `s.report("case.docx")` (.pptx / .html / .md) |
| compare | `ressim.compare(ressim.open_summary("OPM/CASE.SMSPEC"), s, ["FOPR", "FPR"])` |

Examples: `examples/quickstart.py`, `examples/steering.py` (PYACTION-style water-cut control),
`examples/sweep.py` (permeability multiplier sweep), `examples/compare_opm.py`.

## Geomodel scripting (C43)

`from ressim.geo import Geomodel` (or `ressim.Geomodel`) wraps the `rs_geo_*` C API:

| Task | Code |
|---|---|
| project | `g = Geomodel()`, `g.demo(dir)` |
| surfaces | `g.load_surface(path)`, `g.surface("TopRes")` (dict with lattice + z array), `g.set_surface(...)`, `g.sample(s, x, y)`, `g.calc("TopRes + 25", onto="TopRes", name="T25")`, `g.grid_points(xyz, name)`, `g.write_surface(s, path)` |
| wells | `g.load_well_heads(p)`, `g.load_deviation(well, p)`, `g.load_tops(p)`, `g.load_las(p)`, `g.wells()`, `g.well(i)`, `g.log(well, "PHIE")`, `g.position(well, md)` |
| polygons | `g.load_polygons(p)`, `g.add_polygon(name, pts)`, `g.polygons()`, `g.write_polygons(p)`, `g.clip(s, poly, keep_inside=True)` |
| velocity | `g.layer_cake(time_surfaces, tops, modes, velocities, datum, residual)` |
| seismic | `g.load_segy(p)`, `g.seismic()`, `g.trace(i, j)`, `g.seismic_calc("ENVELOPE * (T > 1500)")`, `g.crop(i=(0, 20), step=(2, 2, 1))`, `g.write_segy(p)` |

Example: `examples/geo_scripting.py`.
