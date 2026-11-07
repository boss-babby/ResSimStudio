"""ResSim Python API (C33).

Thin ctypes wrapper of the ResSim engine C API (ressim.dll on Windows, libressim.so on Linux) - no compiler, no
pybind11 build and no particular Python version needed. numpy / pandas are used when installed (arrays / DataFrames),
plain lists otherwise.

    import ressim
    s = ressim.Session("SPE1CASE1.DATA")        # ECLIPSE deck, .rsim project or ressim.Session.example("spe1")
    s.set("PERMX", "300*200")                   # edit the deck (ECLIPSE keyword data up to '/')
    s.run()                                     # whole schedule, any engine
    t, fopr = s.summary("TIME"), s.summary("FOPR")
    p = s.cell("PRESSURE", step=-1)             # per active cell (full=True -> nx*ny*nz with NaN for inactive)
    s.report("spe1.docx")                       # Word / PowerPoint / HTML report

Step-wise runs with steering (black-oil engine):

    s.begin()
    while s.step():
        if s.summary("FWCT")[-1] > 0.8:
            s.apply_schedule("WELOPEN\\n 'PROD' SHUT /\\n/\\n")

Geomodel scripting (C43): ``from ressim.geo import Geomodel`` (also ``ressim.Geomodel``) - surfaces, wells, logs, polygons,
layer-cake depth conversion, seismic calculator / crop, IRAP and SEG-Y output; see ``python/examples/geo_scripting.py``.

The library is searched in: $RESSIM_LIB, the package folder, the folders above it (release layout), PATH.
"""
from __future__ import annotations

import ctypes as _C
import json as _json
import math as _math
import os as _os
import sys as _sys

try:  # optional
    import numpy as _np
except Exception:  # pragma: no cover
    _np = None

__all__ = ["Session", "Summary", "open_summary", "compare", "examples", "version", "set_threads", "ResSimError", "library_path", "Geomodel"]
__version__ = "1.0.0"


class ResSimError(RuntimeError):
    pass


# ------------------------------------------------------------------------------------------------ library
def _candidates():
    names = ["ressim.dll"] if _sys.platform.startswith("win") else (["libressim.dylib"] if _sys.platform == "darwin" else ["libressim.so"])
    env = _os.environ.get("RESSIM_LIB")
    if env:
        yield env
    here = _os.path.dirname(_os.path.abspath(__file__))
    d = here
    for _ in range(5):
        for n in names:
            yield _os.path.join(d, n)
            yield _os.path.join(d, "build", "native", n)
            yield _os.path.join(d, "build", "mingw", n)
            yield _os.path.join(d, "release", n)
        d = _os.path.dirname(d)
    for n in names:
        yield n


_lib = None
library_path = None


def _load():
    global _lib, library_path
    if _lib is not None:
        return _lib
    errors = []
    for p in _candidates():
        if _os.path.sep in p and not _os.path.exists(p):
            continue
        try:
            if _sys.platform.startswith("win") and _os.path.sep in p and hasattr(_os, "add_dll_directory"):
                _os.add_dll_directory(_os.path.dirname(p))
            _lib = _C.CDLL(p)
            library_path = p
            break
        except OSError as e:
            errors.append(f"{p}: {e}")
    if _lib is None:
        raise ResSimError("ResSim library not found - set RESSIM_LIB to ressim.dll / libressim.so\n" + "\n".join(errors[-4:]))
    L = _lib
    S = _C.c_void_p
    sig = {
        "rs_version": ([], _C.c_char_p), "rs_last_error": ([], _C.c_char_p), "rs_set_threads": ([_C.c_int], None),
        "rs_open": ([_C.c_char_p], S), "rs_open_text": ([_C.c_char_p, _C.c_char_p], S), "rs_open_summary": ([_C.c_char_p], S),
        "rs_close": ([S], None), "rs_set_keyword": ([S, _C.c_char_p], _C.c_int), "rs_num_warnings": ([S], _C.c_int),
        "rs_warning": ([S, _C.c_int], _C.c_char_p), "rs_info_json": ([S], _C.c_char_p), "rs_list_examples": ([_C.c_char_p, _C.c_int], _C.c_int),
        "rs_run": ([S], _C.c_int), "rs_begin": ([S], _C.c_int), "rs_step": ([S], _C.c_int), "rs_apply_schedule": ([S, _C.c_char_p], _C.c_int),
        "rs_time": ([S], _C.c_double), "rs_report_step": ([S], _C.c_int), "rs_num_report_steps": ([S], _C.c_int), "rs_run_json": ([S], _C.c_char_p),
        "rs_summary_count": ([S], _C.c_int), "rs_summary_key": ([S, _C.c_int], _C.c_char_p), "rs_summary_unit": ([S, _C.c_int], _C.c_char_p),
        "rs_summary_length": ([S], _C.c_int), "rs_summary_get": ([S, _C.c_char_p, _C.POINTER(_C.c_double), _C.c_int], _C.c_int),
        "rs_grid_dims": ([S] + [_C.POINTER(_C.c_int)] * 4, _C.c_int),
        "rs_static_array": ([S, _C.c_char_p, _C.POINTER(_C.c_double), _C.c_int], _C.c_int),
        "rs_active_index": ([S, _C.POINTER(_C.c_int), _C.c_int], _C.c_int), "rs_num_snapshots": ([S], _C.c_int),
        "rs_snapshot_time": ([S, _C.c_int], _C.c_double), "rs_cell_array_names": ([S, _C.c_int], _C.c_char_p),
        "rs_cell_array": ([S, _C.c_int, _C.c_char_p, _C.POINTER(_C.c_float), _C.c_int], _C.c_int),
        "rs_write_eclipse": ([S, _C.c_char_p], _C.c_int), "rs_export_deck": ([S, _C.c_char_p], _C.c_int),
        "rs_report": ([S, _C.c_char_p, _C.c_char_p], _C.c_int), "rs_compare_json": ([S, S, _C.c_char_p], _C.c_char_p),
    }
    for name, (args, res) in sig.items():
        f = getattr(L, name)
        f.argtypes = args
        f.restype = res
    return L


def _b(s):
    return None if s is None else (s if isinstance(s, bytes) else str(s).encode("utf-8"))


def _s(b):
    return "" if b is None else b.decode("utf-8", "replace")


def _err(what=""):
    msg = _s(_load().rs_last_error())
    raise ResSimError(f"{what}: {msg}" if what else msg)


def _check(rc, what):
    if rc < 0:
        _err(what)
    return rc


def _arr(values, kind="d"):
    if _np is not None:
        return _np.ctypeslib.as_array(values).copy()
    return list(values)


def version() -> str:
    return _s(_load().rs_version())


def set_threads(n: int) -> None:
    _load().rs_set_threads(int(n))


def examples() -> list:
    L = _load()
    n = L.rs_list_examples(None, 0)
    buf = _C.create_string_buffer(n)
    L.rs_list_examples(buf, n)
    return [x for x in _s(buf.value).split("\n") if x]


# ------------------------------------------------------------------------------------------------ summaries
class _SummaryMixin:
    _h = None

    def summary_keys(self) -> list:
        L = _load()
        return [_s(L.rs_summary_key(self._h, i)) for i in range(max(0, L.rs_summary_count(self._h)))]

    def summary_unit(self, key: str) -> str:
        L = _load()
        for i in range(max(0, L.rs_summary_count(self._h))):
            if _s(L.rs_summary_key(self._h, i)) == key:
                return _s(L.rs_summary_unit(self._h, i))
        return ""

    def summary(self, key: str):
        """Vector by ECLIPSE name: 'TIME' (days), 'FOPR', 'WBHP:PROD', ... (one value per time step)."""
        L = _load()
        n = _check(L.rs_summary_get(self._h, _b(key), None, 0), "summary " + key)
        buf = (_C.c_double * n)()
        L.rs_summary_get(self._h, _b(key), buf, n)
        return _arr(buf)

    def to_dataframe(self, keys=None):
        """pandas DataFrame indexed by TIME (days) with the requested (default: all field) vectors."""
        import pandas as pd  # optional dependency
        keys = keys or [k for k in self.summary_keys() if ":" not in k and k != "TIME"]
        return pd.DataFrame({k: self.summary(k) for k in keys}, index=pd.Index(self.summary("TIME"), name="TIME"))


class Summary(_SummaryMixin):
    """Read-only ECLIPSE / OPM Flow summary (SMSPEC + UNSMRY / Snnnn)."""

    def __init__(self, smspec: str):
        self._h = _load().rs_open_summary(_b(smspec))
        if not self._h:
            _err("open summary " + smspec)
        self.path = smspec

    def info(self) -> dict:
        return _json.loads(_s(_load().rs_info_json(self._h)))

    def close(self):
        if self._h:
            _load().rs_close(self._h)
            self._h = None

    __del__ = close

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


def open_summary(smspec: str) -> Summary:
    return Summary(smspec)


def compare(reference, run, keys=None) -> list:
    """Per-vector deviations of `run` against `reference` (Session or Summary) at the reference report times:
    list of dicts key / ref_final / cmp_final / max_abs / rms / rel_rms (rms / max |ref|)."""
    k = ",".join(keys) if keys else ""
    return _json.loads(_s(_load().rs_compare_json(reference._h, run._h, _b(k))))


# ------------------------------------------------------------------------------------------------ session
class Session(_SummaryMixin):
    """A ResSim case: open, edit, run (whole or step-wise with steering) and read results."""

    def __init__(self, path: str | None = None, *, text: str | None = None, base_dir: str | None = None):
        L = _load()
        if text is not None:
            self._h = L.rs_open_text(_b(text), _b(base_dir or ""))
        elif path is not None:
            self._h = L.rs_open(_b(path))
        else:
            raise ValueError("give a deck path or text=")
        if not self._h:
            _err("open")

    @classmethod
    def example(cls, example_id: str) -> "Session":
        return cls("example:" + example_id)

    # ---- lifetime
    def close(self):
        if getattr(self, "_h", None):
            _load().rs_close(self._h)
            self._h = None

    __del__ = close

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    # ---- model
    def info(self) -> dict:
        return _json.loads(_s(_load().rs_info_json(self._h)))

    @property
    def warnings(self) -> list:
        L = _load()
        return [_s(L.rs_warning(self._h, i)) for i in range(L.rs_num_warnings(self._h))]

    def set(self, key: str, value=None) -> None:
        """Edit the deck and re-import it: set('PERMX', '300*200'), set('RESSIM:FAULTSLIP', '20 0.25 0.1')
        or set('PORO=300*0.25')."""
        kv = key if value is None else f"{key}={value}"
        _check(_load().rs_set_keyword(self._h, _b(kv)), "set " + kv)

    @property
    def dims(self):
        a = [_C.c_int() for _ in range(4)]
        _check(_load().rs_grid_dims(self._h, *[_C.byref(x) for x in a]), "dims")
        return tuple(x.value for x in a[:3])

    def static(self, name: str):
        """Input array per global cell (ECLIPSE order): PORO PERMX PERMY PERMZ NTG DX DY DZ TOPS ACTNUM SATNUM ..."""
        L = _load()
        n = _check(L.rs_static_array(self._h, _b(name), None, 0), "static " + name)
        buf = (_C.c_double * n)()
        L.rs_static_array(self._h, _b(name), buf, n)
        return _arr(buf)

    # ---- run
    def run(self) -> dict:
        _check(_load().rs_run(self._h), "run")
        return self.run_info()

    def begin(self) -> None:
        _check(_load().rs_begin(self._h), "begin")

    def step(self) -> bool:
        """Run one report step; True while more steps remain."""
        return _check(_load().rs_step(self._h), "step") > 0

    def steps(self):
        """Generator: begin() if needed, then yields (report_step, time_days) after every report step."""
        if self.num_snapshots == 0:
            self.begin()
        while True:
            more = self.step()
            yield self.report_step, self.time
            if not more:
                break

    def apply_schedule(self, text: str) -> list:
        """SCHEDULE keywords (deck text) acting from the current report step on, like an ACTIONX body (PYACTION-style
        steering). Returns the importer messages (keywords that cannot act at run time)."""
        n = _check(_load().rs_apply_schedule(self._h, _b(text)), "apply_schedule")
        return [m for m in _s(_load().rs_last_error()).split("\n") if m] if n > 0 else []

    @property
    def time(self) -> float:
        return _load().rs_time(self._h)

    @property
    def report_step(self) -> int:
        return _load().rs_report_step(self._h)

    @property
    def num_report_steps(self) -> int:
        return _load().rs_num_report_steps(self._h)

    def run_info(self) -> dict:
        return _json.loads(_s(_load().rs_run_json(self._h)))

    # ---- cell arrays
    @property
    def num_snapshots(self) -> int:
        return max(0, _load().rs_num_snapshots(self._h))

    def snapshot_times(self) -> list:
        L = _load()
        return [L.rs_snapshot_time(self._h, i) for i in range(self.num_snapshots)]

    def cell_names(self, step: int = -1) -> list:
        return [x for x in _s(_load().rs_cell_array_names(self._h, step)).split("\n") if x]

    def active_index(self):
        L = _load()
        n = _check(L.rs_active_index(self._h, None, 0), "active_index")
        buf = (_C.c_int * n)()
        L.rs_active_index(self._h, buf, n)
        return _arr(buf)

    def cell(self, name: str, step: int = -1, full: bool = False):
        """Dynamic array (PRESSURE SWAT SOIL SGAS RS RV ...) of snapshot `step` (0 = initial, -1 = last) per active
        cell, or per global cell (nx*ny*nz, NaN for inactive cells) with full=True."""
        L = _load()
        n = _check(L.rs_cell_array(self._h, int(step), _b(name), None, 0), "cell " + name)
        buf = (_C.c_float * n)()
        L.rs_cell_array(self._h, int(step), _b(name), buf, n)
        v = _arr(buf)
        if not full:
            return v
        nx, ny, nz = self.dims
        idx = self.active_index()
        if _np is not None:
            out = _np.full(nx * ny * nz, _np.nan, dtype=_np.float32)
            idx = _np.asarray(idx)
            ok = idx < nx * ny * nz
            out[idx[ok]] = _np.asarray(v)[ok]
            return out
        out = [_math.nan] * (nx * ny * nz)
        for a, g in enumerate(idx):
            if g < len(out):
                out[g] = v[a]
        return out

    # ---- output
    def write_eclipse(self, base: str) -> None:
        """BASE.EGRID / INIT / SMSPEC / UNSMRY / UNRST (ResInsight, resdata, OPM tools)."""
        _check(_load().rs_write_eclipse(self._h, _b(base)), "write_eclipse")

    def export_deck(self, path: str) -> None:
        _check(_load().rs_export_deck(self._h, _b(path)), "export_deck")

    def report(self, path: str, **options) -> str:
        """Word (.docx), PowerPoint (.pptx), HTML or Markdown report of the run. Options: title, author, maps (bool),
        wells (bool), layer (1-based), maxwells, inputs, log, notes."""
        opt = ";".join(f"{k}={int(v) if isinstance(v, bool) else v}" for k, v in options.items())
        _check(_load().rs_report(self._h, _b(path), _b(opt)), "report")
        return path


def __getattr__(name):  # lazy: ressim.Geomodel without importing the geo module up front (PEP 562)
    if name == "Geomodel":
        from .geo import Geomodel
        return Geomodel
    raise AttributeError(f"module 'ressim' has no attribute {name!r}")
