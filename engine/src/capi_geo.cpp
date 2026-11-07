// C43 C API of the geomodel (capi.h, rs_geo_*): Petrel-parity P10 scripting (Python package `ressim.geo`, MATLAB, Julia ...)
#define RESSIM_CAPI_BUILD
#include "ressim/capi.h"
#include "ressim/geomodel.h"
#include "ressim/geo_tools.h"
#include "ressim/geo_analysis.h"
#include "ressim/geo_import.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

using namespace ressim;
using namespace ressim::geo;

void rs_capi_set_error(const std::string& e);   // capi.cpp

struct rs_geo {
    std::vector<Surface> hz;
    std::vector<GeoWell> wells;
    std::vector<ressim::geo::Polygon> polys;
    SeismicCube cube;
    bool hasCube = false;
    std::string json;
};

namespace {
int gfail(const std::string& e) { rs_capi_set_error(e); return -1; }
const char* gfailS(rs_geo* g, const std::string& e) { rs_capi_set_error(e); if (g) { g->json.clear(); return nullptr; } return nullptr; }
std::string esc(const std::string& s) {
    std::string o;
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += char(c); }
        else if (c == '\n') o += "\\n";
        else if (c < 32) o += ' ';
        else o += char(c);
    }
    return o;
}
std::string num(double v) {
    if (!std::isfinite(v)) return "null";
    char b[40]; std::snprintf(b, sizeof b, "%.12g", v); return b;
}
std::string baseNoExt(const std::string& p) {
    const size_t s = p.find_last_of("/\\");
    std::string b = p.substr(s == std::string::npos ? 0 : s + 1);
    const size_t d = b.rfind('.');
    return d == std::string::npos ? b : b.substr(0, d);
}
int putSurface(rs_geo* g, const Surface& s) {
    for (size_t k = 0; k < g->hz.size(); ++k) if (g->hz[k].name == s.name) { g->hz[k] = s; return int(k); }
    g->hz.push_back(s);
    return int(g->hz.size()) - 1;
}
int wellIndex(rs_geo* g, const std::string& name) {
    for (size_t k = 0; k < g->wells.size(); ++k) if (g->wells[k].name == name) return int(k);
    return -1;
}
#define NEED_G(g) if (!(g)) return gfail("null geomodel handle")
#define NEED_SURF(g, i) NEED_G(g); if ((i) < 0 || (i) >= int((g)->hz.size())) return gfail("surface index out of range")
}  // namespace

extern "C" {

rs_geo* rs_geo_new(void) { return new rs_geo(); }
void rs_geo_free(rs_geo* g) { delete g; }

int rs_geo_demo(rs_geo* g, const char* dir, int field) {
    NEED_G(g);
    if (!dir) return gfail("null directory");
    std::string e;
    const auto files = writeDemoProject(dir, field != 0, e);
    if (files.empty()) return gfail("demo project: " + e);
    const std::string d(dir);
    for (const char* h : {"TopRes", "MidRes", "BaseRes"}) { Surface s = readIrap(d + "/" + h + ".irap", e); s.name = h; putSurface(g, s); }
    auto ws = readWellHeads(d + "/well_heads.txt", e);
    for (auto& w : ws) {
        std::FILE* f = std::fopen((d + "/" + w.name + "_survey.txt").c_str(), "rb");
        if (f) { std::fclose(f); readDeviation(d + "/" + w.name + "_survey.txt", w, e); }
        w.logs.push_back(readLas(d + "/" + w.name + ".las", e));
    }
    for (const auto& t : readTops(d + "/well_tops.txt", e)) for (auto& w : ws) if (w.name == t.first) w.tops.push_back(t.second);
    g->wells = ws;
    g->cube = readSegy(d + "/seismic.segy", e);
    g->hasCube = g->cube.ni > 0;
    return int(files.size());
}

// ---------------------------------------------------------------- surfaces
int rs_geo_load_surface(rs_geo* g, const char* path) {
    NEED_G(g);
    if (!path) return gfail("null path");
    std::string e; Surface s = readIrap(path, e);
    if (!s.valid()) { e.clear(); s = readZmap(path, e); }
    if (!s.valid()) return gfail(e.empty() ? std::string("not an IRAP / ZMAP surface: ") + path : e);
    s.name = baseNoExt(path);
    return putSurface(g, s);
}
int rs_geo_num_surfaces(rs_geo* g) { NEED_G(g); return int(g->hz.size()); }
const char* rs_geo_surface_json(rs_geo* g, int i) {
    if (!g || i < 0 || i >= int(g->hz.size())) return gfailS(g, "surface index out of range");
    const Surface& s = g->hz[size_t(i)];
    double lo, hi; s.range(lo, hi);
    int def = 0; for (double z : s.z) def += std::isfinite(z);
    std::ostringstream o;
    o << "{\"name\":\"" << esc(s.name) << "\",\"nx\":" << s.nx << ",\"ny\":" << s.ny << ",\"x0\":" << num(s.x0) << ",\"y0\":" << num(s.y0)
      << ",\"dx\":" << num(s.dx) << ",\"dy\":" << num(s.dy) << ",\"rot\":" << num(s.rotDeg) << ",\"min\":" << num(def ? lo : NAN)
      << ",\"max\":" << num(def ? hi : NAN) << ",\"defined\":" << def << "}";
    g->json = o.str();
    return g->json.c_str();
}
int rs_geo_surface_values(rs_geo* g, int i, double* buf, int len) {
    NEED_SURF(g, i);
    const Surface& s = g->hz[size_t(i)];
    const int n = int(s.z.size());
    if (buf) std::copy(s.z.begin(), s.z.begin() + std::min(n, std::max(0, len)), buf);
    return n;
}
int rs_geo_set_surface(rs_geo* g, const char* name, int nx, int ny, double x0, double y0, double dx, double dy, double rot, const double* z) {
    NEED_G(g);
    if (!name || !z || nx < 2 || ny < 2 || !(dx > 0) || !(dy > 0)) return gfail("set_surface: need a name, nx, ny >= 2, dx, dy > 0 and values");
    Surface s; s.name = name; s.nx = nx; s.ny = ny; s.x0 = x0; s.y0 = y0; s.dx = dx; s.dy = dy; s.rotDeg = rot;
    s.z.assign(z, z + size_t(nx) * size_t(ny));
    return putSurface(g, s);
}
int rs_geo_surface_sample(rs_geo* g, int i, const double* x, const double* y, int n, double* out) {
    NEED_SURF(g, i);
    if (!x || !y || !out || n < 0) return gfail("null buffer");
    for (int k = 0; k < n; ++k) out[k] = g->hz[size_t(i)].at(x[k], y[k]);
    return n;
}
int rs_geo_surface_calc(rs_geo* g, const char* expr, int onto, const char* name) {
    NEED_SURF(g, onto);
    if (!expr) return gfail("null expression");
    Surface out; std::string e;
    if (!surfaceCalculator(g->hz, onto, expr, out, e)) return gfail(e);
    out.name = name && *name ? name : std::string("calc");
    return putSurface(g, out);
}
int rs_geo_grid_points(rs_geo* g, const double* xyz, int n, int latticeOf, double range, const char* name) {
    NEED_G(g);
    if (!xyz || n < 3) return gfail("grid_points: need at least 3 points");
    std::vector<P3> pts(static_cast<size_t>(n));
    double x0 = 1e300, x1 = -1e300;
    for (int k = 0; k < n; ++k) { pts[size_t(k)] = {xyz[3 * k], xyz[3 * k + 1], xyz[3 * k + 2]}; x0 = std::min(x0, xyz[3 * k]); x1 = std::max(x1, xyz[3 * k]); }
    Surface lat = latticeOf >= 0 && latticeOf < int(g->hz.size()) ? g->hz[size_t(latticeOf)] : latticeFor(pts, 81);
    geostat::Variogram vg; vg.model = geostat::VarioModel::Gaussian;
    vg.rangeMax = vg.rangeMin = range > 0 ? range : std::max(1.0, 0.35 * (x1 - x0));
    Surface s = gridPoints(pts, lat, vg);
    s.name = name && *name ? name : std::string("gridded");
    return putSurface(g, s);
}
int rs_geo_write_surface(rs_geo* g, int i, const char* path) {
    NEED_SURF(g, i);
    if (!path || !writeIrap(path, g->hz[size_t(i)])) return gfail(std::string("cannot write ") + (path ? path : "(null)"));
    return 0;
}

// ---------------------------------------------------------------- wells
int rs_geo_load_well_heads(rs_geo* g, const char* path) {
    NEED_G(g);
    if (!path) return gfail("null path");
    std::string e; auto ws = readWellHeads(path, e);
    if (ws.empty()) return gfail(e.empty() ? "no wells in the file" : e);
    for (auto& w : ws) { const int k = wellIndex(g, w.name); if (k >= 0) { g->wells[size_t(k)].x = w.x; g->wells[size_t(k)].y = w.y; g->wells[size_t(k)].kb = w.kb; } else g->wells.push_back(w); }
    return int(ws.size());
}
int rs_geo_load_deviation(rs_geo* g, const char* well, const char* path) {
    NEED_G(g);
    if (!well || !path) return gfail("null argument");
    const int k = wellIndex(g, well);
    if (k < 0) return gfail(std::string("no well ") + well);
    std::string e;
    if (!readDeviation(path, g->wells[size_t(k)], e)) return gfail(e);
    return int(g->wells[size_t(k)].traj.size());
}
int rs_geo_load_tops(rs_geo* g, const char* path) {
    NEED_G(g);
    if (!path) return gfail("null path");
    std::string e; int n = 0;
    for (const auto& t : readTops(path, e)) { const int k = wellIndex(g, t.first); if (k >= 0) { g->wells[size_t(k)].tops.push_back(t.second); ++n; } }
    if (!n && !e.empty()) return gfail(e);
    return n;
}
int rs_geo_load_las(rs_geo* g, const char* well, const char* path) {
    NEED_G(g);
    if (!path) return gfail("null path");
    std::string e; WellLog L = readLas(path, e);
    if (L.md.empty()) return gfail(e.empty() ? "empty LAS" : e);
    const std::string wn = well && *well ? std::string(well) : (L.well.empty() ? baseNoExt(path) : L.well);
    int k = wellIndex(g, wn);
    if (k < 0) { GeoWell w; w.name = wn; g->wells.push_back(w); k = int(g->wells.size()) - 1; }
    g->wells[size_t(k)].logs.push_back(L);
    return int(L.curves.size());
}
int rs_geo_num_wells(rs_geo* g) { NEED_G(g); return int(g->wells.size()); }
const char* rs_geo_well_json(rs_geo* g, int i) {
    if (!g || i < 0 || i >= int(g->wells.size())) return gfailS(g, "well index out of range");
    const GeoWell& w = g->wells[size_t(i)];
    std::ostringstream o;
    o << "{\"name\":\"" << esc(w.name) << "\",\"x\":" << num(w.x) << ",\"y\":" << num(w.y) << ",\"kb\":" << num(w.kb)
      << ",\"td\":" << num(std::max(w.td, w.totalMd())) << ",\"tops\":[";
    for (size_t t = 0; t < w.tops.size(); ++t) {
        const P3 p = w.at(w.tops[t].md);
        o << (t ? "," : "") << "{\"name\":\"" << esc(w.tops[t].surface) << "\",\"md\":" << num(w.tops[t].md) << ",\"x\":" << num(p[0])
          << ",\"y\":" << num(p[1]) << ",\"z\":" << num(p[2]) << "}";
    }
    o << "],\"logs\":[";
    bool first = true;
    for (const auto& L : w.logs) for (const auto& c : L.curves) { o << (first ? "" : ",") << "\"" << esc(c.mnem) << "\""; first = false; }
    o << "]}";
    g->json = o.str();
    return g->json.c_str();
}
int rs_geo_well_log(rs_geo* g, int i, const char* mnem, double* md, double* v, int len) {
    NEED_G(g);
    if (i < 0 || i >= int(g->wells.size()) || !mnem) return gfail("well index out of range / null mnemonic");
    for (const auto& L : g->wells[size_t(i)].logs) {
        const LogCurve* c = L.curve(mnem);
        if (!c) continue;
        const int n = int(std::min(L.md.size(), c->v.size()));
        for (int k = 0; k < std::min(n, std::max(0, len)); ++k) {
            if (md) md[k] = L.md[size_t(k)];
            if (v) v[k] = (c->v[size_t(k)] == L.nullValue) ? NAN : c->v[size_t(k)];
        }
        return n;
    }
    return gfail(std::string("no log ") + mnem);
}
int rs_geo_well_position(rs_geo* g, int i, const double* md, int n, double* xyz) {
    NEED_G(g);
    if (i < 0 || i >= int(g->wells.size()) || !md || !xyz || n < 0) return gfail("well index out of range / null buffer");
    for (int k = 0; k < n; ++k) { const P3 p = g->wells[size_t(i)].at(md[k]); xyz[3 * k] = p[0]; xyz[3 * k + 1] = p[1]; xyz[3 * k + 2] = p[2]; }
    return n;
}

// ---------------------------------------------------------------- polygons
int rs_geo_load_polygons(rs_geo* g, const char* path) {
    NEED_G(g);
    if (!path) return gfail("null path");
    std::string e; auto v = readPolygons(path, e);
    if (v.empty()) return gfail(e.empty() ? "no polygons in the file" : e);
    if (v.size() == 1 && (v[0].name.empty() || v[0].name.rfind("Polygon", 0) == 0)) v[0].name = baseNoExt(path);
    for (auto& p : v) g->polys.push_back(p);
    return int(v.size());
}
int rs_geo_add_polygon(rs_geo* g, const char* name, int closed, const double* xy, int n) {
    NEED_G(g);
    if (!xy || n < 2) return gfail("polygon: need at least 2 points");
    ressim::geo::Polygon p; p.name = name && *name ? name : "Polygon " + std::to_string(g->polys.size() + 1); p.closed = closed != 0 && n >= 3;
    for (int k = 0; k < n; ++k) p.pts.push_back({xy[2 * k], xy[2 * k + 1]});
    g->polys.push_back(p);
    return int(g->polys.size()) - 1;
}
int rs_geo_num_polygons(rs_geo* g) { NEED_G(g); return int(g->polys.size()); }
const char* rs_geo_polygon_json(rs_geo* g, int i) {
    if (!g || i < 0 || i >= int(g->polys.size())) return gfailS(g, "polygon index out of range");
    const auto& p = g->polys[size_t(i)];
    std::ostringstream o;
    o << "{\"name\":\"" << esc(p.name) << "\",\"closed\":" << (p.closed ? "true" : "false") << ",\"area\":" << num(polygonArea(p))
      << ",\"length\":" << num(polylineLength(p)) << ",\"points\":[";
    for (size_t k = 0; k < p.pts.size(); ++k) o << (k ? "," : "") << "[" << num(p.pts[k][0]) << "," << num(p.pts[k][1]) << "]";
    o << "]}";
    g->json = o.str();
    return g->json.c_str();
}
int rs_geo_write_polygons(rs_geo* g, const char* path) {
    NEED_G(g);
    if (!path || !writePolygons(path, g->polys)) return gfail(std::string("cannot write ") + (path ? path : "(null)"));
    return int(g->polys.size());
}
int rs_geo_clip_surface(rs_geo* g, int surface, int polygon, int keepInside) {
    NEED_SURF(g, surface);
    if (polygon < 0 || polygon >= int(g->polys.size()) || !g->polys[size_t(polygon)].closed) return gfail("need a closed polygon");
    const std::string nm = g->hz[size_t(surface)].name;
    g->hz[size_t(surface)] = clipSurface(g->hz[size_t(surface)], g->polys[size_t(polygon)], keepInside != 0);
    g->hz[size_t(surface)].name = nm;
    int def = 0; for (double z : g->hz[size_t(surface)].z) def += std::isfinite(z);
    return def;
}

// ---------------------------------------------------------------- layer cake
const char* rs_geo_layer_cake(rs_geo* g, const int* ts, const char* topsCsv, const int* mode, const double* vUser, int n, double datum, int residual) {
    if (!g || !ts || n < 1) return gfailS(g, "layer cake: need at least one interface");
    LayerCakeInput in; in.datumZ = datum; in.residual = residual != 0;
    std::vector<std::string> tops;
    if (topsCsv) { std::string t; std::istringstream is(topsCsv); while (std::getline(is, t, ',')) { while (!t.empty() && t.front() == ' ') t.erase(0, 1); while (!t.empty() && t.back() == ' ') t.pop_back(); tops.push_back(t); } }
    for (int k = 0; k < n; ++k) {
        if (ts[k] < 0 || ts[k] >= int(g->hz.size())) return gfailS(g, "layer cake: time surface index out of range");
        in.twt.push_back(g->hz[size_t(ts[k])]);
        in.topName.push_back(size_t(k) < tops.size() ? tops[size_t(k)] : std::string());
        in.mode.push_back(mode ? mode[k] : 0);
        in.vUser.push_back(vUser ? vUser[k] : 2500.0);
    }
    const LayerCakeResult R = layerCakeDepth(in, g->wells);
    std::ostringstream o;
    o << "{\"surfaces\":[";
    for (size_t k = 0; k < R.depth.size(); ++k) { const int idx = putSurface(g, R.depth[k]); o << (k ? "," : "") << idx; }
    o << "],\"vmean\":[";
    for (size_t k = 0; k < R.vMean.size(); ++k) o << (k ? "," : "") << num(R.vMean[k]);
    o << "],\"nwells\":[";
    for (size_t k = 0; k < R.nWells.size(); ++k) o << (k ? "," : "") << R.nWells[k];
    o << "],\"rms_before\":" << num(R.rmsBefore) << ",\"rms_after\":" << num(R.rmsAfter) << ",\"warnings\":[";
    for (size_t k = 0; k < R.warnings.size(); ++k) o << (k ? "," : "") << "\"" << esc(R.warnings[k]) << "\"";
    o << "]}";
    g->json = o.str();
    return g->json.c_str();
}

// ---------------------------------------------------------------- seismic
int rs_geo_load_segy(rs_geo* g, const char* path) {
    NEED_G(g);
    if (!path) return gfail("null path");
    std::string e; SeismicCube c = readSegy(path, e);
    if (!e.empty() || c.ni <= 0) return gfail(e.empty() ? "empty SEG-Y" : e);
    g->cube = c; g->hasCube = true;
    return c.ni * c.nj;
}
const char* rs_geo_seismic_json(rs_geo* g) {
    if (!g || !g->hasCube) return gfailS(g, "no seismic volume");
    const auto& c = g->cube;
    std::ostringstream o;
    o << "{\"name\":\"" << esc(c.name) << "\",\"ni\":" << c.ni << ",\"nj\":" << c.nj << ",\"ns\":" << c.ns << ",\"t0\":" << num(c.t0) << ",\"dt\":" << num(c.dt) << "}";
    g->json = o.str();
    return g->json.c_str();
}
int rs_geo_seismic_trace(rs_geo* g, int i, int j, float* buf, int len) {
    NEED_G(g);
    if (!g->hasCube) return gfail("no seismic volume");
    const auto& c = g->cube;
    if (i < 0 || j < 0 || i >= c.ni || j >= c.nj) return gfail("trace index out of range");
    if (buf) for (int s = 0; s < std::min(c.ns, std::max(0, len)); ++s) buf[s] = c.at(i, j, s);
    return c.ns;
}
int rs_geo_seismic_calc(rs_geo* g, const char* expr, int window) {
    NEED_G(g);
    if (!g->hasCube) return gfail("no seismic volume");
    if (!expr) return gfail("null expression");
    SeismicCube o; std::string e;
    if (!seismicCalculator(g->cube, expr, window > 0 ? window : 6, o, e)) return gfail(e);
    g->cube = o;
    return 0;
}
int rs_geo_seismic_crop(rs_geo* g, int i0, int i1, int j0, int j1, int s0, int s1, int di, int dj, int ds) {
    NEED_G(g);
    if (!g->hasCube) return gfail("no seismic volume");
    g->cube = cropCube(g->cube, i0, i1, j0, j1, s0, s1, std::max(1, di), std::max(1, dj), std::max(1, ds));
    return g->cube.ni * g->cube.nj;
}
int rs_geo_write_segy(rs_geo* g, const char* path) {
    NEED_G(g);
    if (!g->hasCube) return gfail("no seismic volume");
    if (!path || !writeSegy(path, g->cube)) return gfail(std::string("cannot write ") + (path ? path : "(null)"));
    return 0;
}

}  // extern "C"
