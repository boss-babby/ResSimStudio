"""ResSim geomodel scripting (C43, Petrel-parity P10) - ctypes wrapper of the rs_geo_* C API.

    from ressim.geo import Geomodel
    g = Geomodel()
    g.demo("/tmp/geo_demo")                         # or load_surface / load_well_heads / load_tops / load_las / load_segy
    top = g.surface("TopRes")                       # dict: name, nx, ny, x0, y0, dx, dy, rot, z (numpy array or list)
    g.calc("TopRes + 25", onto="TopRes", name="TopRes_shift")
    poly = g.add_polygon("lease", [(x0, y0), (x1, y0), (x1, y1)])
    g.clip("TopRes_shift", poly)
    g.write_surface("TopRes_shift", "top_shift.irap")
    g.seismic_calc("ENVELOPE * (T > 1500)")
    g.write_segy("envelope_below_1500.sgy")

Surfaces are nx*ny arrays, i fastest, NaN undefined; positions use the project coordinates (depth positive down).
"""
from __future__ import annotations

import ctypes as _C
import json as _json

from . import _load, _b, _s, _check, _err, ResSimError

try:
    import numpy as _np
except Exception:  # pragma: no cover
    _np = None

__all__ = ["Geomodel"]

_G = _C.c_void_p
_D = _C.POINTER(_C.c_double)
_I = _C.POINTER(_C.c_int)
_F = _C.POINTER(_C.c_float)
_SIG = {
    "rs_geo_new": ([], _G), "rs_geo_free": ([_G], None), "rs_geo_demo": ([_G, _C.c_char_p, _C.c_int], _C.c_int),
    "rs_geo_load_surface": ([_G, _C.c_char_p], _C.c_int), "rs_geo_num_surfaces": ([_G], _C.c_int),
    "rs_geo_surface_json": ([_G, _C.c_int], _C.c_char_p), "rs_geo_surface_values": ([_G, _C.c_int, _D, _C.c_int], _C.c_int),
    "rs_geo_set_surface": ([_G, _C.c_char_p, _C.c_int, _C.c_int] + [_C.c_double] * 5 + [_D], _C.c_int),
    "rs_geo_surface_sample": ([_G, _C.c_int, _D, _D, _C.c_int, _D], _C.c_int),
    "rs_geo_surface_calc": ([_G, _C.c_char_p, _C.c_int, _C.c_char_p], _C.c_int),
    "rs_geo_grid_points": ([_G, _D, _C.c_int, _C.c_int, _C.c_double, _C.c_char_p], _C.c_int),
    "rs_geo_write_surface": ([_G, _C.c_int, _C.c_char_p], _C.c_int),
    "rs_geo_load_well_heads": ([_G, _C.c_char_p], _C.c_int), "rs_geo_load_deviation": ([_G, _C.c_char_p, _C.c_char_p], _C.c_int),
    "rs_geo_load_tops": ([_G, _C.c_char_p], _C.c_int), "rs_geo_load_las": ([_G, _C.c_char_p, _C.c_char_p], _C.c_int),
    "rs_geo_num_wells": ([_G], _C.c_int), "rs_geo_well_json": ([_G, _C.c_int], _C.c_char_p),
    "rs_geo_well_log": ([_G, _C.c_int, _C.c_char_p, _D, _D, _C.c_int], _C.c_int),
    "rs_geo_well_position": ([_G, _C.c_int, _D, _C.c_int, _D], _C.c_int),
    "rs_geo_load_polygons": ([_G, _C.c_char_p], _C.c_int), "rs_geo_add_polygon": ([_G, _C.c_char_p, _C.c_int, _D, _C.c_int], _C.c_int),
    "rs_geo_num_polygons": ([_G], _C.c_int), "rs_geo_polygon_json": ([_G, _C.c_int], _C.c_char_p),
    "rs_geo_write_polygons": ([_G, _C.c_char_p], _C.c_int), "rs_geo_clip_surface": ([_G, _C.c_int, _C.c_int, _C.c_int], _C.c_int),
    "rs_geo_layer_cake": ([_G, _I, _C.c_char_p, _I, _D, _C.c_int, _C.c_double, _C.c_int], _C.c_char_p),
    "rs_geo_load_segy": ([_G, _C.c_char_p], _C.c_int), "rs_geo_seismic_json": ([_G], _C.c_char_p),
    "rs_geo_seismic_trace": ([_G, _C.c_int, _C.c_int, _F, _C.c_int], _C.c_int),
    "rs_geo_seismic_calc": ([_G, _C.c_char_p, _C.c_int], _C.c_int),
    "rs_geo_seismic_crop": ([_G] + [_C.c_int] * 9, _C.c_int), "rs_geo_write_segy": ([_G, _C.c_char_p], _C.c_int),
}
_ready = False


def _lib():
    global _ready
    L = _load()
    if not _ready:
        for name, (args, res) in _SIG.items():
            try:
                f = getattr(L, name)
            except AttributeError as e:  # older library
                raise ResSimError(f"the ResSim library has no {name} (needs C43 or later)") from e
            f.argtypes = args
            f.restype = res
        _ready = True
    return L


def _dbuf(values):
    vals = [float(v) for v in values]
    return (_C.c_double * max(1, len(vals)))(*vals), len(vals)


def _out(buf, n):
    if _np is not None:
        return _np.ctypeslib.as_array(buf)[:n].copy()
    return list(buf)[:n]


def _js(p, what):
    if not p:
        _err(what)
    return _json.loads(_s(p))


class Geomodel:
    """A geomodel project: surfaces, wells (heads, trajectories, tops, logs), polygons and a seismic volume."""

    def __init__(self):
        self._h = _lib().rs_geo_new()

    def close(self):
        if self._h:
            _lib().rs_geo_free(self._h)
            self._h = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def demo(self, directory: str, field_units: bool = False) -> int:
        """Write the synthetic demo project (3 horizons, 6 wells, SEG-Y) to `directory` and load it."""
        return _check(_lib().rs_geo_demo(self._h, _b(directory), int(field_units)), "demo")

    # ---- surfaces
    def _si(self, s) -> int:
        if isinstance(s, int):
            return s
        for i in range(self.num_surfaces()):
            if _js(_lib().rs_geo_surface_json(self._h, i), "surface")["name"] == s:
                return i
        raise ResSimError(f"no surface {s!r}")

    def num_surfaces(self) -> int:
        return _lib().rs_geo_num_surfaces(self._h)

    def surfaces(self) -> list:
        return [_js(_lib().rs_geo_surface_json(self._h, i), "surface")["name"] for i in range(self.num_surfaces())]

    def load_surface(self, path: str) -> int:
        return _check(_lib().rs_geo_load_surface(self._h, _b(path)), "load_surface")

    def surface(self, s) -> dict:
        i = self._si(s)
        info = _js(_lib().rs_geo_surface_json(self._h, i), "surface")
        n = info["nx"] * info["ny"]
        buf = (_C.c_double * n)()
        _check(_lib().rs_geo_surface_values(self._h, i, buf, n), "surface")
        info["z"] = _out(buf, n)
        return info

    def set_surface(self, name: str, nx: int, ny: int, x0: float, y0: float, dx: float, dy: float, z, rot: float = 0.0) -> int:
        buf, n = _dbuf(list(z))
        if n != nx * ny:
            raise ResSimError(f"set_surface: {n} values for {nx} x {ny}")
        return _check(_lib().rs_geo_set_surface(self._h, _b(name), nx, ny, x0, y0, dx, dy, rot, buf), "set_surface")

    def sample(self, s, x, y):
        xb, n = _dbuf(x)
        yb, _ = _dbuf(y)
        out = (_C.c_double * max(1, n))()
        _check(_lib().rs_geo_surface_sample(self._h, self._si(s), xb, yb, n, out), "sample")
        return _out(out, n)

    def calc(self, expr: str, onto, name: str) -> int:
        """Surface calculator: names of the surfaces (non-alphanumerics as _), X, Y, + - * / ^, if(), min/max/abs/sqrt/log ..."""
        return _check(_lib().rs_geo_surface_calc(self._h, _b(expr), self._si(onto), _b(name)), "calc")

    def grid_points(self, xyz, name: str, lattice_of=None, variogram_range: float = 0.0) -> int:
        flat = [float(v) for p in xyz for v in p]
        buf, n = _dbuf(flat)
        lat = -1 if lattice_of is None else self._si(lattice_of)
        return _check(_lib().rs_geo_grid_points(self._h, buf, n // 3, lat, variogram_range, _b(name)), "grid_points")

    def write_surface(self, s, path: str) -> None:
        _check(_lib().rs_geo_write_surface(self._h, self._si(s), _b(path)), "write_surface")

    # ---- wells
    def load_well_heads(self, path: str) -> int:
        return _check(_lib().rs_geo_load_well_heads(self._h, _b(path)), "load_well_heads")

    def load_deviation(self, well: str, path: str) -> int:
        return _check(_lib().rs_geo_load_deviation(self._h, _b(well), _b(path)), "load_deviation")

    def load_tops(self, path: str) -> int:
        return _check(_lib().rs_geo_load_tops(self._h, _b(path)), "load_tops")

    def load_las(self, path: str, well: str | None = None) -> int:
        return _check(_lib().rs_geo_load_las(self._h, _b(well or ""), _b(path)), "load_las")

    def num_wells(self) -> int:
        return _lib().rs_geo_num_wells(self._h)

    def well(self, i) -> dict:
        if not isinstance(i, int):
            i = [w["name"] for w in self.wells()].index(i)
        return _js(_lib().rs_geo_well_json(self._h, i), "well")

    def wells(self) -> list:
        return [_js(_lib().rs_geo_well_json(self._h, i), "well") for i in range(self.num_wells())]

    def _wi(self, w) -> int:
        return w if isinstance(w, int) else [x["name"] for x in self.wells()].index(w)

    def log(self, well, mnem: str):
        """(md, values) of a log curve; nulls as NaN."""
        i = self._wi(well)
        n = _check(_lib().rs_geo_well_log(self._h, i, _b(mnem), None, None, 0), "log")
        md, v = (_C.c_double * max(1, n))(), (_C.c_double * max(1, n))()
        _lib().rs_geo_well_log(self._h, i, _b(mnem), md, v, n)
        return _out(md, n), _out(v, n)

    def position(self, well, md):
        """x, y, tvd at the measured depths (list of triples / n x 3 array)."""
        mb, n = _dbuf(md)
        out = (_C.c_double * max(1, 3 * n))()
        _check(_lib().rs_geo_well_position(self._h, self._wi(well), mb, n, out), "position")
        a = _out(out, 3 * n)
        return a.reshape(n, 3) if _np is not None else [a[3 * k:3 * k + 3] for k in range(n)]

    # ---- polygons
    def load_polygons(self, path: str) -> int:
        return _check(_lib().rs_geo_load_polygons(self._h, _b(path)), "load_polygons")

    def add_polygon(self, name: str, points, closed: bool = True) -> int:
        buf, n = _dbuf([float(v) for p in points for v in p])
        return _check(_lib().rs_geo_add_polygon(self._h, _b(name), int(closed), buf, n // 2), "add_polygon")

    def polygons(self) -> list:
        return [_js(_lib().rs_geo_polygon_json(self._h, i), "polygon") for i in range(_lib().rs_geo_num_polygons(self._h))]

    def write_polygons(self, path: str) -> int:
        return _check(_lib().rs_geo_write_polygons(self._h, _b(path)), "write_polygons")

    def clip(self, s, polygon: int, keep_inside: bool = True) -> int:
        """Undefine the nodes outside (or inside) a closed polygon; returns the defined nodes left."""
        return _check(_lib().rs_geo_clip_surface(self._h, self._si(s), polygon, int(keep_inside)), "clip")

    # ---- layer cake
    def layer_cake(self, time_surfaces, tops, modes=None, velocities=None, datum: float = 0.0, residual: bool = True) -> dict:
        """Depth surfaces from time surfaces (TWT ms) and well tops; modes per layer 0 mean / 1 map / 2 user velocity."""
        idx = [self._si(s) for s in time_surfaces]
        n = len(idx)
        ib = (_C.c_int * n)(*idx)
        mb = (_C.c_int * n)(*(modes or [0] * n))
        vb = (_C.c_double * n)(*(velocities or [2500.0] * n))
        r = _js(_lib().rs_geo_layer_cake(self._h, ib, _b(",".join(tops)), mb, vb, n, datum, int(residual)), "layer_cake")
        r["surface_names"] = [_js(_lib().rs_geo_surface_json(self._h, i), "surface")["name"] for i in r["surfaces"]]
        return r

    # ---- seismic
    def load_segy(self, path: str) -> int:
        return _check(_lib().rs_geo_load_segy(self._h, _b(path)), "load_segy")

    def seismic(self) -> dict:
        return _js(_lib().rs_geo_seismic_json(self._h), "seismic")

    def trace(self, i: int, j: int):
        n = _check(_lib().rs_geo_seismic_trace(self._h, i, j, None, 0), "trace")
        buf = (_C.c_float * max(1, n))()
        _lib().rs_geo_seismic_trace(self._h, i, j, buf, n)
        return _out(buf, n)

    def seismic_calc(self, expr: str, window: int = 6) -> None:
        """Replace the volume: AMP, ENVELOPE, PHASE, FREQUENCY, COSPHASE, SWEETNESS, RMS, RELIMPEDANCE, COHERENCE, VARIANCE, T, I, J."""
        _check(_lib().rs_geo_seismic_calc(self._h, _b(expr), window), "seismic_calc")

    def crop(self, i=(0, -1), j=(0, -1), s=(0, -1), step=(1, 1, 1)) -> int:
        info = self.seismic()
        fix = lambda r, n: (r[0], n - 1 if r[1] < 0 else r[1])
        i0, i1 = fix(i, info["ni"]); j0, j1 = fix(j, info["nj"]); s0, s1 = fix(s, info["ns"])
        return _check(_lib().rs_geo_seismic_crop(self._h, i0, i1, j0, j1, s0, s1, *step), "crop")

    def write_segy(self, path: str) -> None:
        _check(_lib().rs_geo_write_segy(self._h, _b(path)), "write_segy")
