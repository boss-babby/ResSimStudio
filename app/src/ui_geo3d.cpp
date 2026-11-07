// C36 Geomodel 3-D scene + seismic interpretation views (geoscene3d.h renderer):
//   * seismic inline / crossline / time slices at full resolution (out-of-core SEG-Y or in-memory cube), any attribute,
//     GPU ray-marched volume of the proxy / attribute cube, depth (velocity model) or time domain
//   * horizons (contoured, coloured by depth), the tracked time horizon, fault surfaces from sticks, well tubes coloured
//     by a log with tops, the corner-point geomodel grid (exposed faces of an I/J/K box filter, cell edges)
//   * interpretation in 3-D and on the 2-D section: tracking seeds and fault-stick digitising; flattened sections
// Positions are sent to the renderer relative to a local origin (float precision of UTM coordinates).
#include "i18n.h"
#include "geo_state.h"
#include "ui_geoimport.h"
#include "geoscene3d.h"
#include "widgets.h"
#include "implot.h"
#include "icons.h"
#include "theme.h"
#include "colormaps.h"
#include "plotkit.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>

// ============================================================================ seismic source helpers
bool SeisGeom::ij(double x, double y, double& i, double& j) const {
    const double det = ix * jy - iy * jx;
    if (std::fabs(det) < 1e-30) return false;
    const double dx = x - x0, dy = y - y0;
    i = (dx * jy - dy * jx) / det; j = (ix * dy - iy * dx) / det;
    return i >= -0.5 && i <= ni - 0.5 && j >= -0.5 && j <= nj - 0.5;
}
SeisGeom cubeGeom(const SeismicCube& c) {
    SeisGeom g;
    g.ni = c.ni; g.nj = c.nj; g.ns = c.ns; g.il0 = c.il0; g.xl0 = c.xl0; g.dil = c.dil; g.dxl = c.dxl;
    g.t0 = c.t0; g.dt = c.dt; g.x0 = c.x0; g.y0 = c.y0; g.ix = c.ix; g.iy = c.iy; g.jx = c.jx; g.jy = c.jy;
    return g;
}
bool seisGeom(const App::GeoState& S, SeisGeom& g) {
    if (S.segy && S.segy->isOpen()) {
        const SegyVolume& v = *S.segy;
        g.ni = v.ni; g.nj = v.nj; g.ns = v.ns; g.il0 = v.il0; g.xl0 = v.xl0; g.dil = v.dil; g.dxl = v.dxl;
        g.t0 = v.t0; g.dt = v.dt; g.x0 = v.x0; g.y0 = v.y0; g.ix = v.ix; g.iy = v.iy; g.jx = v.jx; g.jy = v.jy;
        return g.ni > 0 && g.nj > 0 && g.ns > 1;
    }
    if (!S.hasCube || S.cube.ni < 1) return false;
    g = cubeGeom(S.cube);
    return g.nj > 0 && g.ns > 1;
}
std::vector<float> seisSlice(const App::GeoState& S, int kind, int idx) {
    if (S.segy && S.segy->isOpen()) {
        const SegyVolume& v = *S.segy;
        if (kind == 0) return v.inlineSlice(std::clamp(idx, 0, v.ni - 1));
        if (kind == 1) return v.xlineSlice(std::clamp(idx, 0, v.nj - 1));
        return v.timeSlice(std::clamp(idx, 0, v.ns - 1));
    }
    const SeismicCube& c = S.cube;
    std::vector<float> out;
    if (!S.hasCube || c.ni < 1 || c.a.size() < size_t(c.ni) * size_t(c.nj) * size_t(c.ns)) return out;
    if (kind == 0) {
        const int i = std::clamp(idx, 0, c.ni - 1);
        out.assign(c.a.begin() + long(size_t(i) * size_t(c.nj) * size_t(c.ns)), c.a.begin() + long(size_t(i + 1) * size_t(c.nj) * size_t(c.ns)));
    } else if (kind == 1) {
        const int j = std::clamp(idx, 0, c.nj - 1);
        out.resize(size_t(c.ni) * size_t(c.ns));
        for (int i = 0; i < c.ni; ++i) std::copy_n(&c.a[(size_t(i) * size_t(c.nj) + size_t(j)) * size_t(c.ns)], size_t(c.ns), &out[size_t(i) * size_t(c.ns)]);
    } else {
        const int s = std::clamp(idx, 0, c.ns - 1);
        out.resize(size_t(c.ni) * size_t(c.nj));
        for (size_t t = 0; t < out.size(); ++t) out[t] = c.a[t * size_t(c.ns) + size_t(s)];
    }
    return out;
}
std::vector<float> seisAttrSection(const App::GeoState& S, int kind, int idx, SeisAttr a) {
    SeisGeom g;
    if (!seisGeom(S, g)) return {};
    std::vector<float> sec = seisSlice(S, kind, idx);
    if (a == SeisAttr::Amplitude || sec.empty()) return sec;
    const int ntr = kind == 0 ? g.nj : g.ni, nmax = kind == 0 ? g.ni : g.nj;
    if (seisAttrMultiTrace(a)) {
        std::vector<float> p, n;
        if (idx > 0) p = seisSlice(S, kind, idx - 1);
        if (idx + 1 < nmax) n = seisSlice(S, kind, idx + 1);
        return sectionAttribute(a, sec, ntr, g.ns, g.dt, S.attrWin, p.empty() ? nullptr : &p, n.empty() ? nullptr : &n);
    }
    return sectionAttribute(a, sec, ntr, g.ns, g.dt, S.attrWin);
}
void seisAttrRange(SeisAttr a, const std::vector<float>& v, double clipPct, float& lo, float& hi) {
    const float c = amplitudeClip(v, clipPct);
    switch (a) {
        case SeisAttr::Phase: lo = -180; hi = 180; break;
        case SeisAttr::CosPhase: lo = -1; hi = 1; break;
        case SeisAttr::Coherence: lo = 0; hi = 1; break;
        case SeisAttr::Amplitude: case SeisAttr::RelImpedance: lo = -c; hi = c; break;
        default: lo = 0; hi = c; break;
    }
    if (!(hi > lo)) { lo = 0; hi = 1; }
}
double trackedTimeAt(const App::GeoState& S, double x, double y) {
    if (!S.hasTh || !S.hasCube || S.th.ni != S.cube.ni || S.th.nj != S.cube.nj) return std::nan("");
    double fi, fj;
    if (!S.cube.ij(x, y, fi, fj)) return std::nan("");
    const int i = std::clamp(int(std::lround(fi)), 0, S.th.ni - 1), j = std::clamp(int(std::lround(fj)), 0, S.th.nj - 1);
    return S.th.t[size_t(i) * size_t(S.th.nj) + size_t(j)];
}

// ============================================================================ 3-D view state
struct Geo3DView {
    GeoScene3D sc;
    bool tried = false, ok = false;
    double ox = 0, oy = 0;                 // local origin
    bool originSet = false;
    // display
    bool domainTime = false;
    bool showIL = true, showXL = true, showTS = false, showVol = false;
    bool showLines = true;                 // C42 2-D line curtains + picks
    float seisOpacity = 1.0f; int seisCmap = int(Colormap::Seismic); bool seisRev = false;
    float volOpacity = 0.55f, volThreshold = 0.35f; int volRes = 160;
    int vi0 = 0, vi1 = 1 << 20, vj0 = 0, vj1 = 1 << 20;   // volume box (proxy indices)
    std::vector<char> hzOn; float hzOpacity = 0.85f; bool hzContours = true; int hzCmap = int(Colormap::Turbo);
    bool showTh = true, showSeeds = true;
    bool showFaults = true, faultSticks = true; float faultOpacity = 0.75f;
    bool showWells = true, wellNames = true, wellTops = true; int wellLog = -1;
    bool showGrid = false, gridEdges = true; int gridProp = 2; float gridOpacity = 1.0f;
    bool vFilter = false; float vLo = 0, vHi = 1;   // C42 value filter: hide cells outside [vLo, vHi] of the displayed property
    int gi[2] = {0, 1 << 20}, gj[2] = {0, 1 << 20}, gk[2] = {0, 1 << 20};
    bool panel = true;
    // caches
    std::string kSlices, kMeshes, kSolids, kVol, kBounds;
    std::shared_ptr<SeismicCube> attrCube; std::string attrKey, attrJobKey;
    std::future<std::shared_ptr<SeismicCube>> attrJob;
    float rLo = -1, rHi = 1; std::string rKey;
    // last slices (hover read-out)
    std::vector<float> secIL, secXL, secTS; int curIL = -1, curXL = -1, curTS = -1; bool tsFromAttr = false;
    std::string legend; float legLo = 0, legHi = 1; Colormap legCm = Colormap::Seismic; bool legRev = false;
    std::vector<GeoLabel> labels;
    ~Geo3DView() { sc.shutdown(); }
};

namespace {
struct Key {
    std::ostringstream o;
    template <class T> Key& operator<<(const T& v) { o << v << '|'; return *this; }
    std::string str() const { return o.str(); }
};
float tScaleOf(const App::GeoState& S) { return S.vm.v0 > 0 ? float(S.vm.v0 / 2000.0) : 1.0f; }
// display z of a TWT (ms) / of a depth
double zOfT(const App::GeoState& S, const Geo3DView& V, double t) { return V.domainTime ? t * tScaleOf(S) : timeToDepth(S.vm, t); }
double zOfDepth(const App::GeoState& S, const Geo3DView& V, double d) { return V.domainTime ? depthToTime(S.vm, d) * tScaleOf(S) : d; }
double tOfZ(const App::GeoState& S, const Geo3DView& V, double z) { return V.domainTime ? z / tScaleOf(S) : depthToTime(S.vm, z); }
[[maybe_unused]] double depthOfZ(const App::GeoState& S, const Geo3DView& V, double z) { return V.domainTime ? timeToDepth(S.vm, z / tScaleOf(S)) : z; }
void rgbOf(Colormap m, float t, float rgb[3]) { colormapRGB(m, std::clamp(t, 0.0f, 1.0f), rgb); }
void colOf(ImU32 c, float a, float out[4]) { const ImVec4 v = ImGui::ColorConvertU32ToFloat4(c); out[0] = v.x; out[1] = v.y; out[2] = v.z; out[3] = a; }
void addSeg(std::vector<GeoLineVtx>& L, float x0, float y0, float z0, float x1, float y1, float z1, const float c[4]) {
    L.push_back({x0, y0, z0, c[0], c[1], c[2], c[3]}); L.push_back({x1, y1, z1, c[0], c[1], c[2], c[3]});
}
// tube along a polyline (n sides), colour per point
void addTube(GeoSolidDesc& D, const std::vector<std::array<float, 3>>& p, const std::vector<std::array<float, 4>>& c, float r, float zs, int sides = 8) {
    if (p.size() < 2) return;
    const unsigned base = unsigned(D.pos.size() / 3);
    for (size_t k = 0; k < p.size(); ++k) {
        // tangent in display space (z exaggerated) so the tube looks round on screen
        const size_t a = k ? k - 1 : 0, b = std::min(p.size() - 1, k + 1);
        float t[3] = {p[b][0] - p[a][0], p[b][1] - p[a][1], (p[b][2] - p[a][2]) * zs};
        float tl = std::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]); if (tl < 1e-12f) { t[0] = 0; t[1] = 0; t[2] = 1; tl = 1; }
        for (float& q : t) q /= tl;
        float u[3] = {0, 0, 1};
        if (std::fabs(t[2]) > 0.9f) { u[0] = 1; u[2] = 0; }
        float n1[3] = {t[1] * u[2] - t[2] * u[1], t[2] * u[0] - t[0] * u[2], t[0] * u[1] - t[1] * u[0]};
        float l1 = std::sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]); for (float& q : n1) q /= l1;
        float n2[3] = {t[1] * n1[2] - t[2] * n1[1], t[2] * n1[0] - t[0] * n1[2], t[0] * n1[1] - t[1] * n1[0]};
        for (int s = 0; s < sides; ++s) {
            const float ang = 6.2831853f * float(s) / float(sides), cs = std::cos(ang), sn = std::sin(ang);
            const float d[3] = {cs * n1[0] + sn * n2[0], cs * n1[1] + sn * n2[1], cs * n1[2] + sn * n2[2]};
            // back to model space: display z offset / zs
            D.pos.insert(D.pos.end(), {p[k][0] + r * d[0], p[k][1] + r * d[1], p[k][2] + r * d[2] / zs});
            D.nrm.insert(D.nrm.end(), {d[0], d[1], d[2] * zs});   // model-space normal (the shader divides z by zScale)
            D.col.insert(D.col.end(), {c[k][0], c[k][1], c[k][2], c[k][3]});
        }
    }
    for (size_t k = 0; k + 1 < p.size(); ++k)
        for (int s = 0; s < sides; ++s) {
            const unsigned a0 = base + unsigned(k * size_t(sides) + size_t(s)), a1 = base + unsigned(k * size_t(sides) + size_t((s + 1) % sides));
            const unsigned b0 = a0 + unsigned(sides), b1 = a1 + unsigned(sides);
            D.idx.insert(D.idx.end(), {a0, b0, a1, a1, b0, b1});
        }
}
// resample a stick by arc length to m points
std::vector<P3> resampleStick(const std::vector<P3>& s, int m) {
    std::vector<P3> out;
    if (s.empty()) return out;
    if (s.size() == 1) { out.assign(size_t(m), s[0]); return out; }
    std::vector<double> L(s.size(), 0.0);
    for (size_t k = 1; k < s.size(); ++k) L[k] = L[k - 1] + std::sqrt(std::pow(s[k][0] - s[k - 1][0], 2) + std::pow(s[k][1] - s[k - 1][1], 2) + std::pow(s[k][2] - s[k - 1][2], 2));
    for (int q = 0; q < m; ++q) {
        const double l = L.back() * q / double(m - 1);
        size_t k = 1; while (k + 1 < s.size() && L[k] < l) ++k;
        const double f = L[k] > L[k - 1] ? std::clamp((l - L[k - 1]) / (L[k] - L[k - 1]), 0.0, 1.0) : 0.0;
        out.push_back({s[k - 1][0] + f * (s[k][0] - s[k - 1][0]), s[k - 1][1] + f * (s[k][1] - s[k - 1][1]), s[k - 1][2] + f * (s[k][2] - s[k - 1][2])});
    }
    return out;
}
std::vector<std::string> logCurves(const App::GeoState& S) {
    std::vector<std::string> c;
    for (const auto& w : S.wells) for (const auto& l : w.logs) for (const auto& x : l.curves) if (std::find(c.begin(), c.end(), x.mnem) == c.end()) c.push_back(x.mnem);
    return c;
}
double niceStep(double range) {
    if (!(range > 0)) return 1;
    const double r = range / 10, p = std::pow(10.0, std::floor(std::log10(r))), m = r / p;
    return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}
const char* kAttrNames[] = {"Amplitude", "Envelope", "Inst. phase", "Inst. frequency", "Cosine of phase", "Sweetness", "RMS amplitude", "Relative impedance", "Coherence (semblance)", "Variance"};
}  // namespace

// ============================================================================ scene building
static void buildSlices(App::GeoState& S, Geo3DView& V, const SeisGeom& g, bool hasSeis) {
    std::vector<GeoSliceDesc> out;
    V.secIL.clear(); V.secXL.clear(); V.secTS.clear(); V.curIL = V.curXL = V.curTS = -1;
    const SeisAttr a = SeisAttr(S.secAttr);
    const Colormap cm = Colormap(V.seisCmap);
    const int nvz = V.domainTime ? 2 : 48;
    auto setStyle = [&](GeoSliceDesc& D) { D.lo = V.rLo; D.hi = V.rHi; D.cmap = cm; D.reverse = V.seisRev; D.opacity = V.seisOpacity; D.alphaMode = 0; };
    auto vertical = [&](int kind, int idx, const std::vector<float>& sec) {
        const int ntr = kind == 0 ? g.nj : g.ni;
        if (sec.size() < size_t(ntr) * size_t(g.ns) || ntr < 2) return;
        GeoSliceDesc D;
        D.nu = 2; D.nv = nvz;
        double xa, ya, xb, yb;
        if (kind == 0) { g.xy(idx, 0, xa, ya); g.xy(idx, g.nj - 1, xb, yb); }
        else { g.xy(0, idx, xa, ya); g.xy(g.ni - 1, idx, xb, yb); }
        for (int b = 0; b < D.nv; ++b) {
            const double t = g.t0 + (g.ns - 1) * g.dt * b / double(D.nv - 1), z = zOfT(S, V, t);
            D.pos.insert(D.pos.end(), {float(xa - V.ox), float(ya - V.oy), float(z), float(xb - V.ox), float(yb - V.oy), float(z)});
        }
        D.tw = ntr; D.th = g.ns; D.tex.resize(size_t(ntr) * size_t(g.ns));
        for (int tr = 0; tr < ntr; ++tr) for (int s = 0; s < g.ns; ++s) D.tex[size_t(s) * size_t(ntr) + size_t(tr)] = sec[size_t(tr) * size_t(g.ns) + size_t(s)];
        setStyle(D);
        out.push_back(std::move(D));
    };
    if (hasSeis && V.showIL) { V.curIL = std::clamp(S.il, 0, g.ni - 1); V.secIL = seisAttrSection(S, 0, V.curIL, a); vertical(0, V.curIL, V.secIL); }
    if (hasSeis && V.showXL) { V.curXL = std::clamp(S.xl, 0, g.nj - 1); V.secXL = seisAttrSection(S, 1, V.curXL, a); vertical(1, V.curXL, V.secXL); }
    if (hasSeis && V.showTS) {
        // amplitude: full resolution from the source; attributes: from the attribute cube (proxy resolution)
        SeisGeom tg = g; std::vector<float> ts;
        const double t = S.tSlice < 0 ? 0.5 * (g.t0 + g.tEnd()) : S.tSlice;
        V.tsFromAttr = false;
        if (a == SeisAttr::Amplitude) { V.curTS = std::clamp(int(std::lround((t - g.t0) / g.dt)), 0, g.ns - 1); ts = seisSlice(S, 2, V.curTS); }
        else if (V.attrCube && V.attrKey == V.attrJobKey) {
            const SeismicCube& c = *V.attrCube; tg = cubeGeom(c);
            V.curTS = std::clamp(int(std::lround((t - c.t0) / c.dt)), 0, c.ns - 1);
            ts.resize(size_t(c.ni) * size_t(c.nj));
            for (size_t q = 0; q < ts.size(); ++q) ts[q] = c.a[q * size_t(c.ns) + size_t(V.curTS)];
            V.tsFromAttr = true;
        }
        if (ts.size() == size_t(tg.ni) * size_t(tg.nj) && tg.ni > 1 && tg.nj > 1) {
            GeoSliceDesc D; D.nu = 2; D.nv = 2;
            const double z = zOfT(S, V, tg.t0 + V.curTS * tg.dt);
            for (int b = 0; b < 2; ++b) for (int u = 0; u < 2; ++u) {
                double x, y; tg.xy(u * (tg.ni - 1), b * (tg.nj - 1), x, y);
                D.pos.insert(D.pos.end(), {float(x - V.ox), float(y - V.oy), float(z)});
            }
            D.tw = tg.ni; D.th = tg.nj; D.tex.resize(ts.size());
            for (int i = 0; i < tg.ni; ++i) for (int j = 0; j < tg.nj; ++j) D.tex[size_t(j) * size_t(tg.ni) + size_t(i)] = ts[size_t(i) * size_t(tg.nj) + size_t(j)];
            setStyle(D);
            V.secTS = ts;
            out.push_back(std::move(D));
        }
    }
    // C42 2-D lines as curtains (traces decimated to <= 2048 x 1024 texels), colour range per line from the clip percentile
    if (V.showLines)
        for (size_t l = 0; l < S.lines.size(); ++l) {
            if (l < S.lineOn.size() && !S.lineOn[l]) continue;
            const SeismicLine& L = S.lines[l];
            if (L.ntr() < 2 || L.ns < 2) continue;
            const int du = std::max(1, (L.ntr() + 2047) / 2048), ds = std::max(1, (L.ns + 1023) / 1024);
            const int tw = (L.ntr() + du - 1) / du, th = (L.ns + ds - 1) / ds;
            GeoSliceDesc D;
            D.nu = std::max(2, std::min(tw, 256)); D.nv = nvz;
            for (int b = 0; b < D.nv; ++b) {
                const double t = L.t0 + (L.ns - 1) * L.dt * b / double(D.nv - 1), z = zOfT(S, V, t);
                for (int u = 0; u < D.nu; ++u) {
                    const double f = double(u) / (D.nu - 1) * (L.ntr() - 1);
                    const int i0 = std::min(int(f), L.ntr() - 2); const double w = f - i0;
                    const double x = L.x[size_t(i0)] * (1 - w) + L.x[size_t(i0 + 1)] * w, y = L.y[size_t(i0)] * (1 - w) + L.y[size_t(i0 + 1)] * w;
                    D.pos.insert(D.pos.end(), {float(x - V.ox), float(y - V.oy), float(z)});
                }
            }
            D.tw = tw; D.th = th; D.tex.resize(size_t(tw) * size_t(th));
            std::vector<float> mag; mag.reserve(D.tex.size() / 5 + 1);
            for (int r = 0; r < th; ++r) for (int u = 0; u < tw; ++u) {
                const float a = L.at(std::min(L.ntr() - 1, u * du), std::min(L.ns - 1, r * ds));
                D.tex[size_t(r) * size_t(tw) + size_t(u)] = a;
                if (((r * tw + u) % 5) == 0) mag.push_back(std::fabs(a));
            }
            float hi = 1;
            if (!mag.empty()) { const size_t q = std::min(mag.size() - 1, size_t(S.lineClip * double(mag.size() - 1))); std::nth_element(mag.begin(), mag.begin() + long(q), mag.end()); hi = std::max(1e-12f, mag[q]); }
            D.lo = -hi; D.hi = hi; D.cmap = cm; D.reverse = V.seisRev; D.opacity = V.seisOpacity; D.alphaMode = 0;
            out.push_back(std::move(D));
        }
    V.sc.setSlices(out);
}

static void buildMeshes(App::GeoState& S, Geo3DView& V) {
    std::vector<GeoMeshDesc> out;
    // ---- horizons (coloured by depth, contoured)
    V.hzOn.resize(S.hz.size(), 1);
    for (size_t h = 0; h < S.hz.size(); ++h) {
        if (!V.hzOn[h] || !S.hz[h].valid()) continue;
        const Surface& s = S.hz[h];
        const int st = std::max(1, std::max(s.nx, s.ny) / 300);
        const int nu = (s.nx - 1) / st + 1, nv = (s.ny - 1) / st + 1;
        GeoMeshDesc D;
        double lo, hi; s.range(lo, hi);
        D.lo = float(lo); D.hi = float(hi); D.cmap = Colormap(V.hzCmap); D.reverse = true; D.opacity = V.hzOpacity;
        D.contour = V.hzContours ? float(niceStep(hi - lo)) : 0.0f;
        std::vector<char> okv(size_t(nu) * size_t(nv), 0);
        for (int b = 0; b < nv; ++b)
            for (int a = 0; a < nu; ++a) {
                const int i = std::min(s.nx - 1, a * st), j = std::min(s.ny - 1, b * st);
                double x, y; s.nodeXY(i, j, x, y);
                const double z = s.node(i, j);
                okv[size_t(b) * size_t(nu) + size_t(a)] = std::isfinite(z);
                D.pos.insert(D.pos.end(), {float(x - V.ox), float(y - V.oy), float(std::isfinite(z) ? zOfDepth(S, V, z) : 0.0)});
                D.val.push_back(std::isfinite(z) ? float(z) : NAN);
            }
        for (int b = 0; b + 1 < nv; ++b)
            for (int a = 0; a + 1 < nu; ++a) {
                const unsigned p = unsigned(b * nu + a), q = p + 1, r = p + unsigned(nu), t = r + 1;
                const bool P = okv[p], Q = okv[q], R = okv[r], T = okv[t];
                if (P && Q && R) D.idx.insert(D.idx.end(), {p, q, r});
                if (Q && T && R) D.idx.insert(D.idx.end(), {q, t, r});
                if (!R && P && Q && T) D.idx.insert(D.idx.end(), {p, q, t});
                if (!Q && P && T && R) D.idx.insert(D.idx.end(), {p, t, r});
            }
        out.push_back(std::move(D));
    }
    // ---- tracked time horizon on the cube lattice
    if (V.showTh && S.hasTh && S.hasCube && S.th.ni == S.cube.ni && S.th.nj == S.cube.nj && S.th.ni > 1 && S.th.nj > 1) {
        const SeismicCube& c = S.cube;
        GeoMeshDesc D;
        double lo = 1e300, hi = -1e300;
        for (double t : S.th.t) if (std::isfinite(t)) { lo = std::min(lo, t); hi = std::max(hi, t); }
        if (lo <= hi) {
            D.lo = float(lo); D.hi = float(hi + 1e-6); D.cmap = Colormap::Viridis; D.reverse = true; D.opacity = 1.0f; D.contour = V.hzContours ? float(niceStep(hi - lo)) : 0.0f;
            for (int i = 0; i < c.ni; ++i)
                for (int j = 0; j < c.nj; ++j) {
                    const double t = S.th.t[size_t(i) * size_t(c.nj) + size_t(j)];
                    double x, y; c.xy(i, j, x, y);
                    D.pos.insert(D.pos.end(), {float(x - V.ox), float(y - V.oy), float(std::isfinite(t) ? zOfT(S, V, t) : 0.0)});
                    D.val.push_back(std::isfinite(t) ? float(t) : NAN);
                }
            auto ok = [&](int i, int j) { return std::isfinite(S.th.t[size_t(i) * size_t(c.nj) + size_t(j)]); };
            for (int i = 0; i + 1 < c.ni; ++i)
                for (int j = 0; j + 1 < c.nj; ++j) {
                    const unsigned p = unsigned(i * c.nj + j), q = p + 1, r = p + unsigned(c.nj), t = r + 1;
                    if (ok(i, j) && ok(i, j + 1) && ok(i + 1, j)) D.idx.insert(D.idx.end(), {p, q, r});
                    if (ok(i, j + 1) && ok(i + 1, j + 1) && ok(i + 1, j)) D.idx.insert(D.idx.end(), {q, t, r});
                }
            out.push_back(std::move(D));
        }
    }
    // ---- geomodel grid: exposed faces of the I/J/K box filter
    if (V.showGrid && S.hasGrid) {
        const GeoGrid& G = S.g;
        const int nx = G.grid.nx, ny = G.grid.ny, nz = G.grid.nz;
        auto prop = geoGridProperty(S, V.gridProp);
        const auto* tpl = geoTemplate(S, geoGridPropName(S, V.gridProp));   // C43 P2 template (colour map, range, log)
        const int i0 = std::clamp(V.gi[0], 0, nx - 1), i1 = std::clamp(V.gi[1], i0, nx - 1);
        const int j0 = std::clamp(V.gj[0], 0, ny - 1), j1 = std::clamp(V.gj[1], j0, ny - 1);
        const int k0 = std::clamp(V.gk[0], 0, nz - 1), k1 = std::clamp(V.gk[1], k0, nz - 1);
        const bool hasCoord = G.grid.coord.size() == size_t(6) * size_t(nx + 1) * size_t(ny + 1);
        auto active = [&](int i, int j, int k) {
            if (i < i0 || i > i1 || j < j0 || j > j1 || k < k0 || k > k1) return false;
            const size_t c = size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k));
            if (V.vFilter && !(prop[c] >= double(V.vLo) && prop[c] <= double(V.vHi))) return false;
            if (!geoCellVisible(S, c)) return false;   // C43 P7 shared zone / block / facies filter
            return G.grid.actnum.empty() || G.grid.actnum[c] != 0;
        };
        auto corner = [&](int i, int j, int k, int a, int b, int f, float out3[3]) {
            const double z = G.grid.zcorn[size_t(((2 * k + f) * 2 * ny + 2 * j + b) * 2 * nx + 2 * i + a)];
            double x, y;
            if (hasCoord) {
                const double* p = &G.grid.coord[size_t(6) * size_t((i + a) + (nx + 1) * (j + b))];
                const double w = std::fabs(p[5] - p[2]) > 1e-12 ? (z - p[2]) / (p[5] - p[2]) : 0.0;
                x = p[0] + w * (p[3] - p[0]); y = p[1] + w * (p[4] - p[1]);
            } else G.cornerXY(i + a, j + b, x, y);
            out3[0] = float(x - V.ox); out3[1] = float(y - V.oy); out3[2] = float(zOfDepth(S, V, z));
        };
        double lo = 1e300, hi = -1e300;
        for (double v : prop) if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
        if (lo > hi) { lo = 0; hi = 1; }
        geoTemplateRange(tpl, lo, hi);
        GeoMeshDesc D;
        D.lo = float(lo); D.hi = float(hi > lo ? hi : lo + 1); D.opacity = V.gridOpacity; D.edges = V.gridEdges;
        D.cmap = tpl ? Colormap(std::clamp(tpl->cmap, 0, int(Colormap::Count) - 1)) : Colormap::Turbo; D.reverse = tpl && tpl->reverse;
        // faces: corner (a, b, f) lists, counter-clockwise seen from outside (orientation does not matter: two-sided lighting)
        static const int F[6][4][3] = {{{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}}, {{1, 0, 0}, {1, 0, 1}, {1, 1, 1}, {1, 1, 0}},
                                       {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 0, 0}}, {{0, 1, 0}, {1, 1, 0}, {1, 1, 1}, {0, 1, 1}},
                                       {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {{0, 0, 1}, {0, 1, 1}, {1, 1, 1}, {1, 0, 1}}};
        static const int N[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
        static const float UV[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        size_t faces = 0;
        for (int k = k0; k <= k1 && faces < 3000000; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    if (!active(i, j, k)) continue;
                    const double v = geoTemplateValue(tpl, prop[size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k))]);
                    for (int f = 0; f < 6; ++f) {
                        if (active(i + N[f][0], j + N[f][1], k + N[f][2])) {
                            // faulted neighbour: the shared face is only partially covered - draw it when the depths differ
                            if (f >= 4) continue;
                            bool same = true;
                            const int ni2 = i + N[f][0], nj2 = j + N[f][1];
                            for (int q = 0; q < 4 && same; ++q) {
                                float p1[3], p2[3];
                                corner(i, j, k, F[f][q][0], F[f][q][1], F[f][q][2], p1);
                                corner(ni2, nj2, k, F[f][q][0] - N[f][0], F[f][q][1] - N[f][1], F[f][q][2], p2);
                                same = std::fabs(p1[2] - p2[2]) < 1e-3f * (1.0f + std::fabs(p1[2]));
                            }
                            if (same) continue;
                        }
                        const unsigned base = unsigned(D.pos.size() / 3);
                        for (int q = 0; q < 4; ++q) {
                            float p[3]; corner(i, j, k, F[f][q][0], F[f][q][1], F[f][q][2], p);
                            D.pos.insert(D.pos.end(), {p[0], p[1], p[2]});
                            D.val.push_back(std::isfinite(v) ? float(v) : NAN);
                            D.uv.insert(D.uv.end(), {UV[q][0], UV[q][1]});
                        }
                        D.idx.insert(D.idx.end(), {base, base + 1, base + 2, base, base + 2, base + 3});
                        ++faces;
                    }
                }
        if (!D.idx.empty()) out.push_back(std::move(D));
    }
    V.sc.setMeshes(out);
}

static void buildSolidsLines(App::GeoState& S, Geo3DView& V, const float lo[3], const float hi[3]) {
    std::vector<GeoSolidDesc> solids;
    std::vector<GeoLineVtx> lines;
    V.labels.clear();
    const float ext = std::max({hi[0] - lo[0], hi[1] - lo[1], 1.0f});
    const float zs = V.sc.zScale;
    // ---- faults: surfaces between consecutive sticks + stick lines
    if (V.showFaults) {
        auto addFault = [&](const Fault& f, ImU32 col, bool draft) {
            GeoSolidDesc D; D.transparent = V.faultOpacity < 0.99f;
            float c[4]; colOf(col, draft ? 0.6f : V.faultOpacity, c);
            const int M = 16;
            std::vector<std::vector<P3>> rs;
            for (const auto& s : f.sticks) if (s.size() >= 2) rs.push_back(resampleStick(s, M));
            for (size_t a = 0; a + 1 < rs.size(); ++a) {
                const unsigned base = unsigned(D.pos.size() / 3);
                for (size_t q = 0; q < 2; ++q)
                    for (int m = 0; m < M; ++m) {
                        const P3& p = rs[a + q][size_t(m)];
                        D.pos.insert(D.pos.end(), {float(p[0] - V.ox), float(p[1] - V.oy), float(zOfDepth(S, V, p[2]))});
                        D.col.insert(D.col.end(), {c[0], c[1], c[2], c[3]});
                    }
                for (int m = 0; m + 1 < M; ++m) {
                    const unsigned p0 = base + unsigned(m), p1 = p0 + 1, q0 = base + unsigned(M + m), q1 = q0 + 1;
                    D.idx.insert(D.idx.end(), {p0, q0, p1, p1, q0, q1});
                }
            }
            if (!D.idx.empty()) solids.push_back(std::move(D));
            if (V.faultSticks || draft) {
                float lc[4] = {1, 1, 1, 0.9f}; if (draft) { lc[0] = 1; lc[1] = 0.85f; lc[2] = 0.1f; }
                for (const auto& s : f.sticks)
                    for (size_t k = 0; k + 1 < s.size(); ++k)
                        addSeg(lines, float(s[k][0] - V.ox), float(s[k][1] - V.oy), float(zOfDepth(S, V, s[k][2])), float(s[k + 1][0] - V.ox), float(s[k + 1][1] - V.oy), float(zOfDepth(S, V, s[k + 1][2])), lc);
            }
            if (!f.sticks.empty() && !f.sticks[0].empty() && !f.name.empty()) {
                const P3& p = f.sticks[0][0];
                GeoLabel L; L.p[0] = float(p[0] - V.ox); L.p[1] = float(p[1] - V.oy); L.p[2] = float(zOfDepth(S, V, p[2])); L.text = f.name;
                L.rgba[0] = 0.7f; L.rgba[1] = 0.15f; L.rgba[2] = 0.12f; L.small = true;
                V.labels.push_back(L);
            }
        };
        for (size_t f = 0; f < S.faults.size(); ++f) addFault(S.faults[f], seriesColor(int(f) + 3), false);
        if (!S.draft.sticks.empty()) addFault(S.draft, IM_COL32(255, 210, 30, 255), true);
    }
    // ---- fault stick being digitised
    if (!S.stick.empty()) {
        const float c[4] = {1, 0.9f, 0.1f, 1};
        for (size_t k = 0; k < S.stick.size(); ++k) {
            const float x = float(S.stick[k][0] - V.ox), y = float(S.stick[k][1] - V.oy), z = float(zOfDepth(S, V, S.stick[k][2]));
            const float d = ext * 0.006f;
            addSeg(lines, x - d, y, z, x + d, y, z, c); addSeg(lines, x, y - d, z, x, y + d, z, c);
            if (k + 1 < S.stick.size()) addSeg(lines, x, y, z, float(S.stick[k + 1][0] - V.ox), float(S.stick[k + 1][1] - V.oy), float(zOfDepth(S, V, S.stick[k + 1][2])), c);
        }
    }
    // ---- tracking seeds
    if (V.showSeeds && S.hasCube) {
        const float c[4] = {0.15f, 1, 0.35f, 1};
        for (const auto& sd : S.seeds) {
            double x, y; S.cube.xy(sd[0], sd[1], x, y);
            const float X = float(x - V.ox), Y = float(y - V.oy), Z = float(zOfT(S, V, sd[2])), d = ext * 0.008f;
            addSeg(lines, X - d, Y, Z, X + d, Y, Z, c); addSeg(lines, X, Y - d, Z, X, Y + d, Z, c);
            addSeg(lines, X, Y, Z - d / zs, X, Y, Z + d / zs, c);
        }
    }
    // ---- wells: tubes (coloured by a log), tops, names
    if (V.showWells && !S.wells.empty()) {
        const auto curves = logCurves(S);
        const std::string curve = V.wellLog >= 0 && V.wellLog < int(curves.size()) ? curves[size_t(V.wellLog)] : std::string();
        double clo = 1e300, chi = -1e300;
        if (!curve.empty())
            for (const auto& w : S.wells) for (const auto& l : w.logs) if (const LogCurve* c = l.curve(curve)) for (double v : c->v) if (std::isfinite(v)) { clo = std::min(clo, v); chi = std::max(chi, v); }
        const bool logCol = clo < chi;
        // display window: from a little above the data to the bottom of the well
        const double zTopDisp = lo[2] - 0.05 * (hi[2] - lo[2]);
        GeoSolidDesc D;
        const float r = ext * 0.0035f;
        for (size_t wi = 0; wi < S.wells.size(); ++wi) {
            const GeoWell& w = S.wells[wi];
            double mdEnd = w.totalMd();
            for (const auto& l : w.logs) if (!l.md.empty()) mdEnd = std::max(mdEnd, l.md.back());
            for (const auto& t : w.tops) mdEnd = std::max(mdEnd, t.md + 20);
            if (!(mdEnd > 0)) continue;
            std::vector<std::array<float, 3>> p; std::vector<std::array<float, 4>> c;
            const int n = 160;
            ImVec4 base = ImGui::ColorConvertU32ToFloat4(seriesColor(int(wi)));
            for (int q = 0; q <= n; ++q) {
                const double md = mdEnd * q / n;
                const P3 P = w.at(md);
                const double z = zOfDepth(S, V, P[2]);
                if (z < zTopDisp && q < n) continue;
                p.push_back({float(P[0] - V.ox), float(P[1] - V.oy), float(z)});
                std::array<float, 4> cc = {base.x * 0.5f + 0.5f, base.y * 0.5f + 0.5f, base.z * 0.5f + 0.5f, 1};
                if (logCol) {
                    double v = std::nan("");
                    for (const auto& l : w.logs) if (const LogCurve* lc = l.curve(curve)) { v = l.valueAt(*lc, md); if (std::isfinite(v)) break; }
                    if (std::isfinite(v)) { float rgb[3]; rgbOf(Colormap::Viridis, float((v - clo) / (chi - clo)), rgb); cc = {rgb[0], rgb[1], rgb[2], 1}; }
                }
                c.push_back(cc);
            }
            if (p.size() >= 2) {
                addTube(D, p, c, logCol ? r * 1.8f : r, zs);
                if (V.wellNames) {
                    GeoLabel L; L.p[0] = p[0][0]; L.p[1] = p[0][1]; L.p[2] = p[0][2]; L.text = w.name;
                    colOf(seriesColor(int(wi)), 0.92f, L.rgba);
                    V.labels.push_back(L);
                }
            }
            if (V.wellTops)
                for (size_t t = 0; t < w.tops.size(); ++t) {
                    const P3 P = w.at(w.tops[t].md);
                    const float X = float(P[0] - V.ox), Y = float(P[1] - V.oy), Z = float(zOfDepth(S, V, P[2])), d = ext * 0.012f;
                    float cc[4]; colOf(seriesColor(int(t) + 1), 1, cc);
                    addSeg(lines, X - d, Y, Z, X + d, Y, Z, cc); addSeg(lines, X, Y - d, Z, X, Y + d, Z, cc);
                    GeoLabel L; L.p[0] = X + d; L.p[1] = Y; L.p[2] = Z; L.text = w.tops[t].surface; L.small = true;
                    L.rgba[0] = cc[0] * 0.6f; L.rgba[1] = cc[1] * 0.6f; L.rgba[2] = cc[2] * 0.6f; L.rgba[3] = 0.85f;
                    V.labels.push_back(L);
                }
        }
        if (!D.idx.empty()) solids.push_back(std::move(D));
    }
    // ---- slice labels
    SeisGeom g;
    if (seisGeom(S, g)) {
        auto lab = [&](double i, double j, double t, const std::string& s) {
            double x, y; g.xy(i, j, x, y);
            GeoLabel L; L.p[0] = float(x - V.ox); L.p[1] = float(y - V.oy); L.p[2] = float(zOfT(S, V, t)); L.text = s; L.small = true;
            L.rgba[0] = 0.12f; L.rgba[1] = 0.12f; L.rgba[2] = 0.14f; L.rgba[3] = 0.8f;
            V.labels.push_back(L);
        };
        if (V.curIL >= 0) lab(V.curIL, 0, g.t0, "IL " + std::to_string(g.il0 + V.curIL * g.dil));
        if (V.curXL >= 0) lab(0, V.curXL, g.t0, "XL " + std::to_string(g.xl0 + V.curXL * g.dxl));
        if (V.curTS >= 0 && V.showTS) {
            char b[48]; const double t = V.tsFromAttr && V.attrCube ? V.attrCube->t0 + V.curTS * V.attrCube->dt : g.t0 + V.curTS * g.dt;
            std::snprintf(b, sizeof b, "%.0f ms", t); lab(g.ni - 1, g.nj - 1, t, b);
        }
    }
    if (V.showLines)   // C42 2-D line picks + line names
        for (size_t l = 0; l < S.lines.size(); ++l) {
            if (l < S.lineOn.size() && !S.lineOn[l]) continue;
            const SeismicLine& L = S.lines[l];
            if (L.ntr() < 2) continue;
            for (size_t h = 0; h < S.lineHz.size(); ++h) {
                if (l >= S.lineHz[h].t.size()) continue;
                const auto& t = S.lineHz[h].t[l];
                static const float hc[6][4] = {{1, 0.9f, 0.1f, 1}, {0.2f, 1, 0.4f, 1}, {0.3f, 0.8f, 1, 1}, {1, 0.45f, 0.9f, 1}, {1, 0.6f, 0.2f, 1}, {0.7f, 0.6f, 1, 1}};
                for (int k = 0; k + 1 < L.ntr() && k + 1 < int(t.size()); ++k)
                    if (std::isfinite(t[size_t(k)]) && std::isfinite(t[size_t(k + 1)]))
                        addSeg(lines, float(L.x[size_t(k)] - V.ox), float(L.y[size_t(k)] - V.oy), float(zOfT(S, V, t[size_t(k)])),
                               float(L.x[size_t(k + 1)] - V.ox), float(L.y[size_t(k + 1)] - V.oy), float(zOfT(S, V, t[size_t(k + 1)])), hc[h % 6]);
            }
            GeoLabel G; G.p[0] = float(L.x[0] - V.ox); G.p[1] = float(L.y[0] - V.oy); G.p[2] = float(zOfT(S, V, L.t0)); G.text = L.name; G.small = true;
            V.labels.push_back(G);
        }
    V.sc.setSolids(solids);
    V.sc.setLines(lines);
    V.sc.setLabels(V.labels);
}

static void buildVolume(App::GeoState& S, Geo3DView& V) {
    const SeisAttr a = SeisAttr(S.secAttr);
    const SeismicCube* src = nullptr;
    if (a == SeisAttr::Amplitude) src = S.hasCube ? &S.cube : nullptr;
    else if (V.attrCube && V.attrKey == V.attrJobKey) src = V.attrCube.get();
    if (!V.showVol || !src || src->ni < 2 || src->nj < 2 || src->ns < 2) { V.sc.setVolume(nullptr); return; }
    const SeismicCube& c = *src;
    const int si = std::max(1, (c.ni + V.volRes - 1) / V.volRes), sj = std::max(1, (c.nj + V.volRes - 1) / V.volRes), ss = std::max(1, (c.ns + 383) / 384);
    GeoVolumeDesc D;
    D.ni = (c.ni - 1) / si + 1; D.nj = (c.nj - 1) / sj + 1; D.ns = (c.ns - 1) / ss + 1;
    D.data.resize(size_t(D.ni) * size_t(D.nj) * size_t(D.ns));
    for (int i = 0; i < D.ni; ++i)
        for (int j = 0; j < D.nj; ++j)
            for (int s = 0; s < D.ns; ++s) D.data[(size_t(i) * size_t(D.nj) + size_t(j)) * size_t(D.ns) + size_t(s)] = c.at(i * si, j * sj, s * ss);
    D.x0 = c.x0 - V.ox; D.y0 = c.y0 - V.oy; D.ix = c.ix * si; D.iy = c.iy * si; D.jx = c.jx * sj; D.jy = c.jy * sj;
    D.t0 = c.t0; D.dt = c.dt * ss;
    D.zMode = V.domainTime ? 0 : 1; D.zScale = tScaleOf(S); D.v0 = S.vm.v0; D.k = S.vm.k; D.datum = S.vm.datum;
    D.lo = V.rLo; D.hi = V.rHi; D.cmap = Colormap(V.seisCmap); D.reverse = V.seisRev; D.opacity = V.volOpacity; D.threshold = V.volThreshold;
    D.iMin = V.vi0 / si; D.iMax = V.vi1 / si; D.jMin = V.vj0 / sj; D.jMax = V.vj1 / sj;
    V.sc.setVolume(&D);
}

static void sceneBounds(App::GeoState& S, Geo3DView& V, float lo[3], float hi[3]) {
    double L[3] = {1e300, 1e300, 1e300}, H[3] = {-1e300, -1e300, -1e300};
    auto add = [&](double x, double y, double z) {
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return;
        const double p[3] = {x - V.ox, y - V.oy, z};
        for (int m = 0; m < 3; ++m) { L[m] = std::min(L[m], p[m]); H[m] = std::max(H[m], p[m]); }
    };
    SeisGeom g;
    if (seisGeom(S, g))
        for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c) {
            double x, y; g.xy(a * (g.ni - 1), b * (g.nj - 1), x, y); add(x, y, zOfT(S, V, c ? g.tEnd() : g.t0));
        }
    for (const auto& h : S.hz) {
        if (!h.valid()) continue;
        for (int j = 0; j < h.ny; j += std::max(1, h.ny / 40)) for (int i = 0; i < h.nx; i += std::max(1, h.nx / 40)) {
            const double z = h.node(i, j); if (!std::isfinite(z)) continue;
            double x, y; h.nodeXY(i, j, x, y); add(x, y, zOfDepth(S, V, z));
        }
    }
    if (S.hasGrid) for (size_t c = 0; c < S.g.center.size(); c += std::max<size_t>(1, S.g.center.size() / 4000)) add(S.g.center[c][0], S.g.center[c][1], zOfDepth(S, V, S.g.center[c][2]));
    for (const auto& f : S.faults) for (const auto& s : f.sticks) for (const auto& p : s) add(p[0], p[1], zOfDepth(S, V, p[2]));
    for (const auto& ln : S.lines) for (int k = 0; k < ln.ntr(); k += std::max(1, ln.ntr() / 50)) { add(ln.x[size_t(k)], ln.y[size_t(k)], zOfT(S, V, ln.t0)); add(ln.x[size_t(k)], ln.y[size_t(k)], zOfT(S, V, ln.t0 + (ln.ns - 1) * ln.dt)); }
    const bool any = L[0] <= H[0];
    for (const auto& w : S.wells) {   // plan extent only when other data exist (well tubes start at the data top)
        if (any) { add(w.x, w.y, L[2]); for (const auto& st : w.traj) add(st.x, st.y, L[2]); }
        else { add(w.x, w.y, 0); for (const auto& t : w.tops) { const P3 P = w.at(t.md); add(P[0], P[1], zOfDepth(S, V, P[2])); } }
    }
    if (L[0] > H[0]) { for (int m = 0; m < 3; ++m) { lo[m] = 0; hi[m] = 1; } return; }
    for (int m = 0; m < 3; ++m) { if (H[m] - L[m] < 1e-6) { L[m] -= 1; H[m] += 1; } lo[m] = float(L[m]); hi[m] = float(H[m]); }
}

// ray -> (slice kind, i, j, t) on the visible seismic slices; false when nothing is hit
static bool rayHitSlices(const App::GeoState& S, const Geo3DView& V, const SeisGeom& g, const double o[3], const double d[3],
                         int& kind, double& fi, double& fj, double& t, double& zHit) {
    const double det = g.ix * g.jy - g.iy * g.jx;
    if (std::fabs(det) < 1e-30) return false;
    const double X0 = g.x0 - V.ox, Y0 = g.y0 - V.oy;
    auto Iof = [&](double x, double y) { return ((x - X0) * g.jy - (y - Y0) * g.jx) / det; };
    auto Jof = [&](double x, double y) { return (g.ix * (y - Y0) - g.iy * (x - X0)) / det; };
    const double zA = zOfT(S, V, g.t0), zB = zOfT(S, V, g.tEnd());
    double best = 1e300; bool hit = false;
    auto test = [&](int k, double s) {
        if (!(s > 0) || s >= best) return;
        const double x = o[0] + s * d[0], y = o[1] + s * d[1], z = o[2] + s * d[2];
        const double I = Iof(x, y), J = Jof(x, y);
        if (I < -0.5 || I > g.ni - 0.5 || J < -0.5 || J > g.nj - 0.5 || z < std::min(zA, zB) - 1e-6 || z > std::max(zA, zB) + 1e-6) return;
        best = s; hit = true; kind = k; fi = I; fj = J; zHit = z; t = tOfZ(S, V, z);
    };
    const double dI = (d[0] * g.jy - d[1] * g.jx) / det, dJ = (g.ix * d[1] - g.iy * d[0]) / det;
    if (V.curIL >= 0 && std::fabs(dI) > 1e-30) test(0, (V.curIL - Iof(o[0], o[1])) / dI);
    if (V.curXL >= 0 && std::fabs(dJ) > 1e-30) test(1, (V.curXL - Jof(o[0], o[1])) / dJ);
    if (V.curTS >= 0 && V.showTS && std::fabs(d[2]) > 1e-30) {
        const double tt = V.tsFromAttr && V.attrCube ? V.attrCube->t0 + V.curTS * V.attrCube->dt : g.t0 + V.curTS * g.dt;
        test(2, (zOfT(S, V, tt) - o[2]) / d[2]);
    }
    return hit;
}

// ============================================================================ interpretation helpers
static void addSeedAt(App::GeoState& S, const SeisGeom& g, double fi, double fj, double t) {
    double x, y; g.xy(fi, fj, x, y);
    double ci, cj;
    if (!S.hasCube || !S.cube.ij(x, y, ci, cj)) { S.msg("seed outside the in-memory cube"); return; }
    S.seeds.push_back({std::clamp(ci, 0.0, double(S.cube.ni - 1)), std::clamp(cj, 0.0, double(S.cube.nj - 1)), t});
    char b[160]; std::snprintf(b, sizeof b, "seed: IL %.0f XL %.0f, %.1f ms", g.il0 + fi * g.dil, g.xl0 + fj * g.dxl, t); S.msg(b);
}
static void addStickPoint(App::GeoState& S, const SeisGeom& g, double fi, double fj, double t) {
    double x, y; g.xy(fi, fj, x, y);
    S.stick.push_back({x, y, timeToDepth(S.vm, t)});
}
static void finishStick(App::GeoState& S) {
    if (S.stick.size() >= 2) {
        std::sort(S.stick.begin(), S.stick.end(), [](const P3& a, const P3& b) { return a[2] < b[2]; });
        S.draft.sticks.push_back(S.stick);
        S.msg("fault stick " + std::to_string(S.draft.sticks.size()) + ": " + std::to_string(S.stick.size()) + " points");
    }
    S.stick.clear();
}
// ============================================================================ C36-4 horizon brush
static const double kZero = 0.0, kOne = 1.0;
void geoLatticeXY(const Surface& h, double fi, double fj, double& x, double& y) {
    double x0, y0, x1, y1, x2, y2;
    h.nodeXY(0, 0, x0, y0); h.nodeXY(1, 0, x1, y1); h.nodeXY(0, 1, x2, y2);
    x = x0 + fi * (x1 - x0) + fj * (x2 - x0); y = y0 + fi * (y1 - y0) + fj * (y2 - y0);
}
int geoApplyBrush(App::GeoState& S, double x, double y) {
    if (S.editHz < 0 || S.editHz >= int(S.hz.size()) || !S.hz[size_t(S.editHz)].valid()) return 0;
    Surface& h = S.hz[size_t(S.editHz)];
    double amount = S.brushAmt;
    if (S.brushMode == 1) amount = std::clamp(S.brushStrength, 0.0, 1.0);
    if (S.brushMode == 2) {
        if (S.flatDepth == 0) { const double z = h.at(x, y); if (!std::isfinite(z)) return 0; S.flatDepth = z; }
        amount = S.flatDepth;
    }
    S.pushUndo();
    const int n = brushSurface(h, x, y, S.brushR, amount, S.brushMode);
    if (n == 0) { S.hzUndo.pop_back(); return 0; }
    ++S.hzRev;
    static const char* mn[] = {"raise / lower", "smooth", "flatten"};
    char b[200]; std::snprintf(b, sizeof b, "brush (%s) on %s at %.0f, %.0f: %d nodes%s", mn[std::clamp(S.brushMode, 0, 2)], h.name.c_str(), x, y, n, S.hasGrid ? " - rebuild the grid to apply" : "");
    S.msg(b);
    return n;
}
// first intersection of a model-space ray with horizon h (bounding-box clip, march + bisection); world x, y and depth
static bool rayHitHorizon(const App::GeoState& S, const Geo3DView& V, int h, const double o[3], const double d[3], double& X, double& Y, double& D) {
    if (h < 0 || h >= int(S.hz.size()) || !S.hz[size_t(h)].valid()) return false;
    const Surface& sf = S.hz[size_t(h)];
    double lo[3] = {1e300, 1e300, 0}, hi[3] = {-1e300, -1e300, 0};
    const int ci[4] = {0, sf.nx - 1, 0, sf.nx - 1}, cj[4] = {0, 0, sf.ny - 1, sf.ny - 1};
    for (int q = 0; q < 4; ++q) { double x, y; sf.nodeXY(ci[q], cj[q], x, y); x -= V.ox; y -= V.oy; lo[0] = std::min(lo[0], x); hi[0] = std::max(hi[0], x); lo[1] = std::min(lo[1], y); hi[1] = std::max(hi[1], y); }
    double zlo, zhi; sf.range(zlo, zhi);
    const double za = zOfDepth(S, V, zlo - 1), zb = zOfDepth(S, V, zhi + 1);
    lo[2] = std::min(za, zb); hi[2] = std::max(za, zb);
    double s0 = 0, s1 = 1e300;
    for (int a = 0; a < 3; ++a) {
        if (std::fabs(d[a]) < 1e-15) { if (o[a] < lo[a] || o[a] > hi[a]) return false; continue; }
        double t0 = (lo[a] - o[a]) / d[a], t1 = (hi[a] - o[a]) / d[a];
        if (t0 > t1) std::swap(t0, t1);
        s0 = std::max(s0, t0); s1 = std::min(s1, t1);
    }
    if (!(s1 > s0)) return false;
    auto f = [&](double t) {
        const double x = o[0] + t * d[0] + V.ox, y = o[1] + t * d[1] + V.oy, z = o[2] + t * d[2];
        const double hz = sf.at(x, y);
        return std::isfinite(hz) ? z - zOfDepth(S, V, hz) : std::nan("");
    };
    const int N = 600;
    double tp = s0, fp = f(s0);
    for (int q = 1; q <= N; ++q) {
        const double t = s0 + (s1 - s0) * q / N, ft = f(t);
        if (std::isfinite(fp) && std::isfinite(ft) && (fp <= 0) != (ft <= 0)) {
            double a = tp, b = t, fa = fp;
            for (int it = 0; it < 50; ++it) { const double m = 0.5 * (a + b), fm = f(m); if (!std::isfinite(fm)) break; if ((fm <= 0) == (fa <= 0)) { a = m; fa = fm; } else b = m; }
            const double tm = 0.5 * (a + b);
            X = o[0] + tm * d[0] + V.ox; Y = o[1] + tm * d[1] + V.oy; D = sf.at(X, Y);
            return std::isfinite(D);
        }
        tp = t; fp = ft;
    }
    return false;
}
static void brushControls(App::GeoState& S) {
    const float fs = ImGui::GetFontSize();
    if (S.hz.empty()) { ImGui::SameLine(); ImGui::TextDisabled("| no horizons"); return; }
    S.editHz = std::clamp(S.editHz, 0, int(S.hz.size()) - 1);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 8);
    if (ImGui::BeginCombo("##ehz", S.hz[size_t(S.editHz)].name.c_str())) {
        for (size_t h = 0; h < S.hz.size(); ++h) if (ImGui::Selectable((S.hz[h].name + "##e" + std::to_string(h)).c_str(), int(h) == S.editHz)) { S.editHz = int(h); S.flatDepth = 0; }
        ImGui::EndCombo();
    }
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7);
    if (ImGui::Combo("##bm", &S.brushMode, "raise / lower\0smooth\0flatten\0")) S.flatDepth = 0;
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::InputDouble("R##br", &S.brushR, 0, 0, "%.0f"); S.brushR = std::max(S.brushR, 1e-3);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5);
    if (S.brushMode == 0) ImGui::InputDouble("dz##ba", &S.brushAmt, 0, 0, "%+.1f");
    else if (S.brushMode == 1) ImGui::SliderScalar("##bs", ImGuiDataType_Double, &S.brushStrength, &kZero, &kOne, "strength %.2f");
    else { ImGui::InputDouble("depth##bf", &S.flatDepth, 0, 0, "%.1f"); }
    ui::hint("Click on the edited horizon (map view or 3-D scene). Raise / lower: dz > 0 deepens; smooth: towards the local mean; flatten: to the "
             "depth (0 = taken from the first click). Gaussian taper exp(-3 r^2/R^2). Ctrl+Z undoes the last edit; rebuild the grid afterwards.");
    ImGui::SameLine();
    ImGui::BeginDisabled(S.hzUndo.empty());
    if (ImGui::SmallButton(ICON_UNDO_2 "##hzu") || (!S.hzUndo.empty() && ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_Z) && !ImGui::GetIO().WantTextInput)) {
        S.hz = S.hzUndo.back(); S.hzUndo.pop_back(); ++S.hzRev; S.msg("horizon edit undone");
    }
    ImGui::EndDisabled();
}
static void interpretationBar(App::GeoState& S) {
    ImGui::TextDisabled("Pick:"); ImGui::SameLine();
    ImGui::RadioButton(TR("navigate"), &S.pickMode, 0); ImGui::SameLine();
    ImGui::RadioButton(TR("tracking seeds"), &S.pickMode, 1); ImGui::SameLine();
    ImGui::RadioButton(TR("fault sticks"), &S.pickMode, 2); ImGui::SameLine();
    ImGui::RadioButton(TR("edit horizon"), &S.pickMode, 3);
    ui::hint("Seeds: click on a seismic section / slice (then Auto-track in the Seismic step). Fault sticks: click points along the fault "
             "trace on successive sections; Enter finishes a stick, Backspace removes the last point, 'Save fault' adds it to the project");
    if (S.pickMode == 2) {
        ImGui::SameLine();
        ImGui::Text("| stick %zu pts, fault %zu sticks", S.stick.size(), S.draft.sticks.size());
        ImGui::SameLine(); if (ImGui::SmallButton(TR("Finish stick"))) finishStick(S);
        ImGui::SameLine(); if (ImGui::SmallButton(TR("Undo"))) { if (!S.stick.empty()) S.stick.pop_back(); else if (!S.draft.sticks.empty()) S.draft.sticks.pop_back(); }
        ImGui::SameLine();
        if (ImGui::SmallButton(TR("Save fault"))) {
            finishStick(S);
            if (S.draft.sticks.size() >= 2) {
                Fault f = S.draft;
                if (f.name.empty()) f.name = "F" + std::to_string(S.faults.size() + 1);
                S.faults.push_back(f); S.msg("fault " + f.name + " saved: " + std::to_string(f.sticks.size()) + " sticks (rebuild the grid to use it)");
                S.draft = Fault();
            } else S.msg("a fault needs at least two sticks");
        }
        ImGui::SameLine(); if (ImGui::SmallButton(TR("Discard"))) { S.stick.clear(); S.draft = Fault(); }
        if (!ImGui::GetIO().WantTextInput) {
            if (ImGui::IsKeyPressed(ImGuiKey_Enter) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter)) finishStick(S);
            if (ImGui::IsKeyPressed(ImGuiKey_Backspace) && !S.stick.empty()) S.stick.pop_back();
        }
    } else if (S.pickMode == 3) {
        brushControls(S);
    } else if (S.pickMode == 1) {
        ImGui::SameLine(); ImGui::Text("| %zu seed(s)", S.seeds.size());
        ImGui::SameLine(); if (ImGui::SmallButton(TR("Clear seeds"))) S.seeds.clear();
    }
}

// ============================================================================ 2-D seismic section (view 1)
void App::geoSeismicSection(GeoState& S, ImVec2 size) {
    if (S.curLine >= 0 && S.curLine < int(S.lines.size())) { geoLineSection(S, size); return; }   // C42 2-D line
    SeisGeom g;
    if (!seisGeom(S, g)) {
        if (!S.lines.empty()) { S.curLine = 0; geoLineSection(S, size); return; }
        ImGui::TextDisabled("No seismic loaded - import SEG-Y files (2-D lines or a 3-D volume) in the Data step (or load the demo project)."); return;
    }
    const float fs = ImGui::GetFontSize();
    if (!S.lines.empty()) { if (ImGui::SmallButton(TR("2-D lines"))) { S.curLine = 0; return; } ImGui::SameLine(); }
    ImGui::SetNextItemWidth(fs * 7); ImGui::Combo("##or", &S.secOrient, "Inline\0Crossline\0");
    ImGui::SameLine();
    const int nmax = S.secOrient == 0 ? g.ni : g.nj;
    int& idx = S.secOrient == 0 ? S.il : S.xl;
    idx = std::clamp(idx, 0, nmax - 1);
    ImGui::SetNextItemWidth(fs * 12); ImGui::SliderInt("##idx", &idx, 0, nmax - 1, S.secOrient == 0 ? "IL index %d" : "XL index %d");
    ImGui::SameLine(); ImGui::Text("%s %d", S.secOrient == 0 ? "inline" : "crossline", S.secOrient == 0 ? g.il0 + idx * g.dil : g.xl0 + idx * g.dxl);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 10); ImGui::Combo("##attr", &S.secAttr, kAttrNames, int(SeisAttr::Count));
    ImGui::SameLine(); ImGui::BeginDisabled(!S.hasTh); ImGui::Checkbox(TR("flatten"), &S.flatten); ImGui::EndDisabled();
    if (!S.hasTh && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("Track a horizon first (Seismic step)");
    if (S.segy) { ImGui::SameLine(); ImGui::TextDisabled("| out-of-core %s", S.segy->name.c_str()); }
    interpretationBar(S);
    // section (cached per orientation / index / attribute / window)
    static std::string cKey; static std::vector<float> cSec; static float cLo = -1, cHi = 1;
    Key k; k << &S << S.segy.get() << S.cube.a.data() << S.cubeRev << S.cube.ni << S.secOrient << idx << S.secAttr << S.attrWin << S.clipPct;
    if (k.str() != cKey) {
        cKey = k.str(); cSec = seisAttrSection(S, S.secOrient, idx, SeisAttr(S.secAttr));
        seisAttrRange(SeisAttr(S.secAttr), cSec, S.clipPct, cLo, cHi);
    }
    const int ntr = S.secOrient == 0 ? g.nj : g.ni;
    if (cSec.size() < size_t(ntr) * size_t(g.ns)) return;
    // flattening on the tracked horizon (per trace sample of the horizon, reference = its mean)
    std::vector<double> hzS(size_t(ntr), std::nan(""));
    double ref = 0; int nref = 0;
    if (S.hasTh) {
        for (int tr = 0; tr < ntr; ++tr) {
            double x, y; if (S.secOrient == 0) g.xy(idx, tr, x, y); else g.xy(tr, idx, x, y);
            const double t = trackedTimeAt(S, x, y);
            if (std::isfinite(t)) { hzS[size_t(tr)] = (t - g.t0) / g.dt; ref += hzS[size_t(tr)]; ++nref; }
        }
        if (nref) ref /= nref;
    }
    const bool flat = S.flatten && nref > 0;
    const std::vector<float> disp = flat ? flattenSection(cSec, ntr, g.ns, hzS, ref) : cSec;
    std::vector<float> img(size_t(g.ns) * size_t(ntr));
    for (int r = 0; r < g.ns; ++r) for (int tr = 0; tr < ntr; ++tr) img[size_t(r) * size_t(ntr) + size_t(tr)] = disp[size_t(tr) * size_t(g.ns) + size_t(g.ns - 1 - r)];
    const SeisAttr a = SeisAttr(S.secAttr);
    const float plotH = size.y - ImGui::GetFrameHeightWithSpacing() * 2;
    if (ImPlot::BeginPlot("##sec", ImVec2(size.x, std::max(100.0f, plotH)))) {
        ImPlot::SetupAxes(S.secOrient == 0 ? "crossline index" : "inline index", flat ? "flattened TWT (ms)" : "TWT (ms)", 0, ImPlotAxisFlags_Invert);
        ImPlot::PushColormap(a == SeisAttr::Amplitude || a == SeisAttr::RelImpedance ? ImPlotColormap_RdBu : a == SeisAttr::Phase ? ImPlotColormap_Twilight : a == SeisAttr::Coherence ? ImPlotColormap_Greys : ImPlotColormap_Viridis);
        ImPlot::PlotHeatmap("##amp", img.data(), g.ns, ntr, cLo, cHi, nullptr, ImPlotPoint(0, g.t0), ImPlotPoint(ntr - 1, g.tEnd()));
        ImPlot::PopColormap();
        auto dispT = [&](int tr, double t) { return flat && std::isfinite(hzS[size_t(tr)]) ? t - (hzS[size_t(tr)] - ref) * g.dt : t; };
        auto realT = [&](int tr, double t) { return flat && std::isfinite(hzS[size_t(tr)]) ? t + (hzS[size_t(tr)] - ref) * g.dt : t; };
        if (S.hasTh) {
            std::vector<double> xs, ys;
            for (int tr = 0; tr < ntr; ++tr) if (std::isfinite(hzS[size_t(tr)])) { xs.push_back(tr); ys.push_back(dispT(tr, g.t0 + hzS[size_t(tr)] * g.dt)); }
            ImPlotSpec sp; sp.LineColor = ImVec4(1, 0.9f, 0.1f, 1); sp.LineWeight = 2;
            ImPlot::PlotLine(S.thName.c_str(), xs.data(), ys.data(), int(xs.size()), sp);
        }
        // seeds near this section
        for (const auto& sd : S.seeds) {
            double x, y; S.cube.xy(sd[0], sd[1], x, y);
            double fi, fj; g.ij(x, y, fi, fj);
            const double along = S.secOrient == 0 ? fi : fj, tr = S.secOrient == 0 ? fj : fi;
            if (std::fabs(along - idx) > 0.75) continue;
            const int tri = std::clamp(int(std::lround(tr)), 0, ntr - 1);
            const double yy = dispT(tri, sd[2]);
            ImPlotSpec sp; sp.Marker = ImPlotMarker_Diamond; sp.MarkerSize = 7; sp.MarkerFillColor = sp.MarkerLineColor = ImVec4(0.1f, 1, 0.3f, 1);
            ImPlot::PlotScatter("seed", &tr, &yy, 1, sp);
        }
        // fault sticks crossing this section (points within half a trace spacing) + the stick being digitised
        auto drawStick = [&](const std::vector<P3>& st, ImVec4 col, const char* id) {
            std::vector<double> xs, ys;
            for (const auto& p : st) {
                double fi, fj; g.ij(p[0], p[1], fi, fj);
                const double along = S.secOrient == 0 ? fi : fj, tr = S.secOrient == 0 ? fj : fi;
                if (std::fabs(along - idx) > 1.5) continue;
                const int tri = std::clamp(int(std::lround(tr)), 0, ntr - 1);
                xs.push_back(tr); ys.push_back(dispT(tri, depthToTime(S.vm, p[2])));
            }
            ImPlotSpec sp; sp.LineColor = col; sp.LineWeight = 2.5f; sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 4; sp.MarkerFillColor = sp.MarkerLineColor = col;
            if (!xs.empty()) ImPlot::PlotLine(id, xs.data(), ys.data(), int(xs.size()), sp);
        };
        for (const auto& f : S.faults) for (const auto& st : f.sticks) drawStick(st, ImVec4(1, 0.3f, 0.25f, 1), ("##f" + f.name).c_str());
        for (const auto& st : S.draft.sticks) drawStick(st, ImVec4(1, 0.75f, 0.1f, 1), "##draft");
        if (!S.stick.empty()) drawStick(S.stick, ImVec4(1, 0.95f, 0.2f, 1), "##stick");
        if (ImPlot::IsPlotHovered()) {
            const ImPlotPoint p = ImPlot::GetPlotMousePos();
            const int tr = std::clamp(int(std::lround(p.x)), 0, ntr - 1);
            const double t = realT(tr, p.y);
            const int s = std::clamp(int(std::lround((t - g.t0) / g.dt)), 0, g.ns - 1);
            ImGui::SetTooltip("%s %d  %s %d\nTWT %.1f ms  depth %.1f\n%s %.4g", S.secOrient == 0 ? "IL" : "XL", S.secOrient == 0 ? g.il0 + idx * g.dil : g.xl0 + idx * g.dxl,
                              S.secOrient == 0 ? "XL" : "IL", S.secOrient == 0 ? g.xl0 + tr * g.dxl : g.il0 + tr * g.dil, t, timeToDepth(S.vm, t), kAttrNames[S.secAttr], cSec[size_t(tr) * size_t(g.ns) + size_t(s)]);
            if (ImGui::IsMouseClicked(0) && !ImGui::GetIO().KeyCtrl && S.pickMode > 0) {
                const double fi = S.secOrient == 0 ? idx : p.x, fj = S.secOrient == 0 ? p.x : idx;
                if (S.pickMode == 1) addSeedAt(S, g, fi, fj, t);
                else addStickPoint(S, g, fi, fj, t);
            }
        }
        ImPlot::EndPlot();
    }
    ImGui::TextDisabled("%s range %.4g .. %.4g %s | clip %.1f %% | window %d samples", kAttrNames[S.secAttr], cLo, cHi, seisAttrUnit(a), S.clipPct * 100, S.attrWin);
}

// ============================================================================ 3-D scene (view 4)
void App::geomodel3D(GeoState& S, ImVec2 size) {
    if (!S.v3) S.v3 = std::make_shared<Geo3DView>();
    Geo3DView& V = *S.v3;
    if (!V.tried) { V.tried = true; V.ok = V.sc.init(); if (!V.ok) S.msg("Error: the 3-D scene could not be initialised (D3D11 shaders)"); }
    if (!V.ok) { ImGui::TextDisabled("3-D scene unavailable (graphics initialisation failed - see the log)."); return; }
    const float fs = ImGui::GetFontSize();
    SeisGeom g; const bool hasSeis = seisGeom(S, g);
    const SeisAttr attr = SeisAttr(S.secAttr);

    // ---- local origin (first data set)
    {
        double ox = 0, oy = 0; bool set = false;
        if (hasSeis) { ox = g.x0; oy = g.y0; set = true; }
        else if (!S.hz.empty()) { ox = S.hz[0].x0; oy = S.hz[0].y0; set = true; }
        else if (S.hasGrid) { ox = S.g.spec.x0; oy = S.g.spec.y0; set = true; }
        else if (!S.wells.empty()) { ox = S.wells[0].x; oy = S.wells[0].y; set = true; }
        else if (!S.lines.empty() && S.lines[0].ntr() > 0) { ox = S.lines[0].x[0]; oy = S.lines[0].y[0]; set = true; }
        if (set && (!V.originSet || std::fabs(ox - V.ox) > 1e5 || std::fabs(oy - V.oy) > 1e5)) { V.ox = std::round(ox); V.oy = std::round(oy); V.originSet = true; V.kBounds.clear(); }
    }
    // ---- attribute cube (background) for time slices / volume of an attribute
    if (V.attrJob.valid() && V.attrJob.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        try { V.attrCube = V.attrJob.get(); V.attrKey = V.attrJobKey; S.msg(std::string("attribute cube ready: ") + seisAttrName(attr)); }
        catch (const std::exception& e) { S.msg(std::string("Error: attribute cube: ") + e.what()); }
        V.kSlices.clear(); V.kVol.clear(); V.rKey.clear();
    }
    if (attr != SeisAttr::Amplitude && S.hasCube && (V.showTS || V.showVol)) {
        Key k; k << S.cube.a.data() << S.cubeRev << S.cube.ni << S.cube.nj << S.cube.ns << S.secAttr << S.attrWin;
        if (k.str() != V.attrJobKey && !V.attrJob.valid()) {
            V.attrJobKey = k.str();
            auto cube = std::make_shared<SeismicCube>(S.cube); const int win = S.attrWin;
            V.attrJob = std::async(std::launch::async, [cube, attr, win]() { return std::make_shared<SeismicCube>(cubeAttribute(*cube, attr, win)); });
            S.msg(std::string("computing the ") + seisAttrName(attr) + " cube ...");
        }
    }
    // ---- colour range of the seismic display (whole proxy / attribute cube: stable between slices)
    {
        Key k; k << S.cube.a.data() << S.cubeRev << S.secAttr << S.clipPct << (V.attrCube ? V.attrCube->a.data() : nullptr) << V.attrKey << S.attrWin;
        if (k.str() != V.rKey && hasSeis) {
            V.rKey = k.str();
            if (attr == SeisAttr::Amplitude && S.hasCube) seisAttrRange(attr, S.cube.a, S.clipPct, V.rLo, V.rHi);
            else if (V.attrCube && V.attrKey == V.attrJobKey) seisAttrRange(attr, V.attrCube->a, S.clipPct, V.rLo, V.rHi);
            else { const auto sec = seisAttrSection(S, 0, std::clamp(S.il, 0, g.ni - 1), attr); seisAttrRange(attr, sec, S.clipPct, V.rLo, V.rHi); V.rKey.clear(); }
            V.kSlices.clear(); V.kVol.clear();
        }
    }

    // ================================================================ layout: scene panel | canvas
    const float panelW = V.panel ? fs * 16.5f : 0.0f;
    if (V.panel) {
        ImGui::BeginChild("##scenepanel", ImVec2(panelW, size.y), ImGuiChildFlags_Borders);
        ImGui::TextDisabled("SCENE");
        if (ImGui::CollapsingHeader(TR(ICON_ACTIVITY " Seismic"), ImGuiTreeNodeFlags_DefaultOpen)) {
            if (!S.lines.empty()) { ImGui::Checkbox(TR("2-D lines (curtains, picks)"), &V.showLines); if (!hasSeis) { ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("opacity##s2", &V.seisOpacity, 0.1f, 1.0f, "%.2f"); } }
            if (!hasSeis) { if (S.lines.empty()) ImGui::TextDisabled("no seismic loaded"); }
            else {
                ImGui::Checkbox(TR("Inline"), &V.showIL); ImGui::SameLine(fs * 5.5f); ImGui::SetNextItemWidth(-1); ImGui::SliderInt("##il3", &S.il, 0, g.ni - 1, "%d");
                ImGui::Checkbox(TR("Crossline"), &V.showXL); ImGui::SameLine(fs * 5.5f); ImGui::SetNextItemWidth(-1); ImGui::SliderInt("##xl3", &S.xl, 0, g.nj - 1, "%d");
                float ts = float(S.tSlice < 0 ? 0.5 * (g.t0 + g.tEnd()) : S.tSlice);
                ImGui::Checkbox(TR("Time"), &V.showTS); ImGui::SameLine(fs * 5.5f); ImGui::SetNextItemWidth(-1);
                if (ImGui::SliderFloat("##ts3", &ts, float(g.t0), float(g.tEnd()), "%.0f ms")) S.tSlice = ts;
                ImGui::SetNextItemWidth(-1); ImGui::Combo("##attr3", &S.secAttr, kAttrNames, int(SeisAttr::Count));
                if (attr != SeisAttr::Amplitude) { ImGui::SetNextItemWidth(fs * 6); ImGui::SliderInt("window", &S.attrWin, 1, 32); }
                std::vector<std::string> cms; for (int c = 0; c < int(Colormap::Count); ++c) cms.push_back(colormapName(Colormap(c)));
                ImGui::SetNextItemWidth(fs * 7); ImGui::Combo("##cm", &V.seisCmap, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(size_t(i)).c_str(); }, &cms, int(cms.size()));
                ImGui::SameLine(); ImGui::Checkbox(TR("rev."), &V.seisRev);
                float clip = float(S.clipPct * 100); ImGui::SetNextItemWidth(fs * 7);
                if (ImGui::SliderFloat("clip %", &clip, 80, 100, "%.1f")) S.clipPct = clip / 100;
                ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("opacity##s", &V.seisOpacity, 0.1f, 1.0f, "%.2f");
                ImGui::Checkbox(TR("Volume rendering"), &V.showVol);
                if (V.showVol) {
                    ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("opacity##v", &V.volOpacity, 0.02f, 1.0f, "%.2f");
                    ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("threshold", &V.volThreshold, 0.0f, 0.95f, "%.2f");
                    ui::hint("Transfer function: samples with |value| below the threshold (fraction of the clip range) are transparent - bright spots and strong reflectors remain");
                    if (S.hasCube) {
                        ImGui::SetNextItemWidth(fs * 9); ImGui::DragIntRange2("IL box", &V.vi0, &V.vi1, 0.2f, 0, S.cube.ni - 1);
                        ImGui::SetNextItemWidth(fs * 9); ImGui::DragIntRange2("XL box", &V.vj0, &V.vj1, 0.2f, 0, S.cube.nj - 1);
                    }
                    ImGui::SetNextItemWidth(fs * 7); ImGui::SliderInt("resolution", &V.volRes, 48, 320);
                }
                if (S.segy) ImGui::TextDisabled("slices: out-of-core, full res.\nvolume: proxy %dx%dx%d", S.cube.ni, S.cube.nj, S.cube.ns);
                if (V.attrJob.valid()) ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(Pal::Warn), "computing attribute cube ...");
            }
        }
        if (ImGui::CollapsingHeader(TR(ICON_MOUNTAIN " Horizons"), ImGuiTreeNodeFlags_DefaultOpen)) {
            V.hzOn.resize(S.hz.size(), 1);
            for (size_t h = 0; h < S.hz.size(); ++h) { bool on = V.hzOn[h] != 0; if (ImGui::Checkbox((S.hz[h].name + "##hz" + std::to_string(h)).c_str(), &on)) V.hzOn[h] = on; }
            if (S.hasTh) ImGui::Checkbox(("tracked: " + S.thName).c_str(), &V.showTh);
            ImGui::Checkbox(TR("contours"), &V.hzContours);
            ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("opacity##h", &V.hzOpacity, 0.1f, 1.0f, "%.2f");
        }
        if (ImGui::CollapsingHeader(TR(ICON_SPLIT " Faults"), ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox(TR("fault surfaces"), &V.showFaults); ImGui::SameLine(); ImGui::Checkbox(TR("sticks"), &V.faultSticks);
            ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("opacity##f", &V.faultOpacity, 0.1f, 1.0f, "%.2f");
            ImGui::Checkbox(TR("seeds"), &V.showSeeds);
        }
        if (ImGui::CollapsingHeader(TR(ICON_PIN " Wells"), ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox(TR("wells"), &V.showWells); ImGui::SameLine(); ImGui::Checkbox(TR("names"), &V.wellNames); ImGui::SameLine(); ImGui::Checkbox(TR("tops"), &V.wellTops);
            auto curves = logCurves(S); curves.insert(curves.begin(), "(no log)");
            int sel = V.wellLog + 1; ImGui::SetNextItemWidth(fs * 9);
            if (ImGui::Combo("log##wl", &sel, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(size_t(i)).c_str(); }, &curves, int(curves.size()))) V.wellLog = sel - 1;
        }
        if (ImGui::CollapsingHeader(TR(ICON_LAYERS " Grid"))) {
            if (!S.hasGrid) ImGui::TextDisabled("build the grid first");
            else {
                ImGui::Checkbox(TR("show grid"), &V.showGrid); ImGui::SameLine(); ImGui::Checkbox(TR("cell edges"), &V.gridEdges);
                ImGui::SetNextItemWidth(fs * 10); geoGridPropCombo("##gp3", &V.gridProp, S);
                V.gi[1] = std::min(V.gi[1], S.g.grid.nx - 1); V.gj[1] = std::min(V.gj[1], S.g.grid.ny - 1); V.gk[1] = std::min(V.gk[1], S.g.grid.nz - 1);
                ImGui::SetNextItemWidth(fs * 9); ImGui::DragIntRange2("I", &V.gi[0], &V.gi[1], 0.2f, 0, S.g.grid.nx - 1);
                ImGui::SetNextItemWidth(fs * 9); ImGui::DragIntRange2("J", &V.gj[0], &V.gj[1], 0.2f, 0, S.g.grid.ny - 1);
                ImGui::SetNextItemWidth(fs * 9); ImGui::DragIntRange2("K", &V.gk[0], &V.gk[1], 0.2f, 0, S.g.grid.nz - 1);
                ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("opacity##g", &V.gridOpacity, 0.1f, 1.0f, "%.2f");
                // C42 value filter (Petrel "value filter"): keeps the cells whose displayed value lies in the range
                {
                    static int rangeFor = -1; static std::string rangeKey; static float plo = 0, phi = 1;
                    const std::string rk = std::to_string(V.gridProp) + "|" + std::to_string(S.propRev) + "|" + std::to_string(size_t(S.phi.v.data())) +
                                           "|" + std::to_string(size_t(S.fr.facies.data())) + "|" + std::to_string(S.userProps.size());
                    if (rangeFor != V.gridProp || rangeKey != rk) {
                        const auto pv = geoGridProperty(S, V.gridProp);
                        double a = 1e300, b = -1e300;
                        for (double x : pv) if (std::isfinite(x)) { a = std::min(a, x); b = std::max(b, x); }
                        if (a > b) { a = 0; b = 1; }
                        plo = float(a); phi = float(b > a ? b : a + 1);
                        if (rangeFor != V.gridProp) { V.vLo = plo; V.vHi = phi; }
                        rangeFor = V.gridProp; rangeKey = rk;
                    }
                    ImGui::Checkbox(TR("value filter"), &V.vFilter);
                    ui::hint("Shows only the cells whose value of the displayed property lies in the range (inclusive); the faces of the cells "
                             "next to hidden ones are drawn, so the filtered body looks solid. Use it with a calculated property "
                             "(Analysis > Property calculator) for conditions such as PORO > 0.2 && ZONE == 2.");
                    if (V.vFilter) {
                        ImGui::SetNextItemWidth(fs * 12);
                        ImGui::DragFloatRange2("range##vf", &V.vLo, &V.vHi, (phi - plo) / 300.0f, plo, phi, "%.4g", "%.4g", ImGuiSliderFlags_AlwaysClamp);
                    }
                }
            }
        }
        if (ImGui::CollapsingHeader(TR(ICON_AXIS_3D " Display"), ImGuiTreeNodeFlags_DefaultOpen)) {
            int dom = V.domainTime ? 1 : 0;
            if (ImGui::Combo("domain", &dom, "Depth (velocity model)\0Time (TWT)\0")) V.domainTime = dom == 1;
            ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("vertical exag.", &V.sc.zScale, 0.5f, 20.0f, "x%.1f", ImGuiSliderFlags_Logarithmic);
            ImGui::Checkbox(TR("legend"), &V.sc.showLegend); ImGui::SameLine(); ImGui::Checkbox(TR("box"), &V.sc.showBox); ImGui::SameLine(); ImGui::Checkbox(TR("labels"), &V.sc.showLabels);
            ImGui::ColorEdit3("background", V.sc.background, ImGuiColorEditFlags_NoInputs);
        }
        ImGui::EndChild();
        ImGui::SameLine();
    }

    // ================================================================ rebuild what changed
    float lo[3], hi[3];
    {
        Key kb; kb << V.ox << V.oy << V.domainTime << S.vm.v0 << S.vm.k << S.vm.datum << S.hz.size() << S.faults.size() << S.wells.size() << S.hasGrid << S.g.center.size()
                   << (S.segy ? S.segy.get() : nullptr) << S.cube.a.data() << S.cube.ni << S.lines.size();
        for (const auto& h : S.hz) kb << h.z.data();
        static float blo[3] = {0, 0, 0}, bhi[3] = {1, 1, 1};
        if (kb.str() != V.kBounds) { sceneBounds(S, V, blo, bhi); V.sc.setBounds(blo, bhi); V.kBounds = kb.str(); V.kSlices.clear(); V.kMeshes.clear(); V.kSolids.clear(); V.kVol.clear(); }
        std::copy(blo, blo + 3, lo); std::copy(bhi, bhi + 3, hi);
    }
    if (hasSeis || (V.showLines && !S.lines.empty())) {
        Key k; k << V.kBounds << V.showIL << V.showXL << V.showTS << S.il << S.xl << S.tSlice << S.secAttr << S.attrWin << V.rLo << V.rHi << V.seisCmap << V.seisRev << V.seisOpacity
                 << V.attrKey << V.attrJobKey << (V.attrCube ? V.attrCube->a.data() : nullptr) << V.showLines << S.lines.size() << S.lineClip;
        for (char c : S.lineOn) k << int(c);
        if (k.str() != V.kSlices) { buildSlices(S, V, g, hasSeis); V.kSlices = k.str(); V.kSolids.clear(); }
    } else if (!V.kSlices.empty() || V.curIL >= 0 || V.curXL >= 0) { V.sc.setSlices({}); V.kSlices.clear(); V.curIL = V.curXL = V.curTS = -1; V.secIL.clear(); V.secXL.clear(); V.secTS.clear(); }
    {
        Key k; k << V.kBounds << V.hzOpacity << V.hzContours << V.hzCmap << V.showTh << S.hasTh << S.th.t.data() << V.showGrid << V.gridProp << V.gridEdges << V.gridOpacity
                 << V.gi[0] << V.gi[1] << V.gj[0] << V.gj[1] << V.gk[0] << V.gk[1] << S.hasPetro << S.phi.v.data() << S.hasFacies << S.fr.facies.data() << S.g.grid.zcorn.data() << S.hzRev << S.hz.size()
                 << S.propRev << S.userProps.size() << S.kx.data() << V.vFilter << V.vLo << V.vHi;
        for (char c : V.hzOn) k << int(c);
        k << S.cf.on;   // C43 shared filter + grid-property template
        for (char c : S.cf.zoneOff) k << int(c);
        for (char c : S.cf.blockOff) k << int(c);
        for (char c : S.cf.faciesOff) k << int(c);
        if (const auto* tp = geoTemplate(S, geoGridPropName(S, V.gridProp))) k << tp->cmap << tp->reverse << tp->fixed << tp->logScale << tp->lo << tp->hi;
        if (k.str() != V.kMeshes) { buildMeshes(S, V); V.kMeshes = k.str(); }
    }
    {
        Key k; k << V.kBounds << V.showFaults << V.faultSticks << V.faultOpacity << V.showWells << V.wellNames << V.wellTops << V.wellLog << V.showSeeds << S.seeds.size()
                 << S.stick.size() << S.draft.sticks.size() << V.sc.zScale << V.curIL << V.curXL << V.curTS << V.showTS << V.showLines << S.lineRev << S.lines.size();
        for (const auto& f : S.faults) k << f.sticks.size() << f.name;
        if (!S.stick.empty()) k << S.stick.back()[0] << S.stick.back()[2];
        if (k.str() != V.kSolids) { buildSolidsLines(S, V, lo, hi); V.kSolids = k.str(); }
    }
    {
        Key k; k << V.kBounds << V.showVol << S.secAttr << V.attrKey << V.attrJobKey << V.volRes << (V.attrCube ? V.attrCube->a.data() : nullptr);
        if (k.str() != V.kVol) { buildVolume(S, V); V.kVol = k.str(); }
        if (V.showVol) {   // style without re-upload
            GeoVolumeDesc st; st.zMode = V.domainTime ? 0 : 1; st.zScale = tScaleOf(S); st.v0 = S.vm.v0; st.k = S.vm.k; st.datum = S.vm.datum;
            st.lo = V.rLo; st.hi = V.rHi; st.cmap = Colormap(V.seisCmap); st.reverse = V.seisRev; st.opacity = V.volOpacity; st.threshold = V.volThreshold;
            const SeismicCube* src = attr == SeisAttr::Amplitude ? (S.hasCube ? &S.cube : nullptr) : V.attrCube.get();
            if (src) {
                const int si = std::max(1, (src->ni + V.volRes - 1) / V.volRes), sj = std::max(1, (src->nj + V.volRes - 1) / V.volRes);
                st.iMin = V.vi0 / si; st.iMax = V.vi1 / si; st.jMin = V.vj0 / sj; st.jMax = V.vj1 / sj;
            }
            V.sc.setVolumeStyle(st);
        }
    }
    // legend: seismic when shown, else the first horizon, else the grid property
    {
        std::string title = "Geomodel";
        if (hasSeis && (V.showIL || V.showXL || V.showTS || V.showVol)) {
            title = std::string(seisAttrName(attr)) + (S.segy ? "  (" + S.segy->name + ")" : std::string());
            V.sc.setTitle(title, std::string(kAttrNames[S.secAttr]) + (seisAttrUnit(attr)[0] ? std::string(" (") + seisAttrUnit(attr) + ")" : ""), V.rLo, V.rHi, Colormap(V.seisCmap), V.seisRev);
        } else if (V.showGrid && S.hasGrid) {
            const auto v = geoGridProperty(S, V.gridProp); double a = 1e300, b = -1e300; for (double x : v) if (std::isfinite(x)) { a = std::min(a, x); b = std::max(b, x); }
            V.sc.setTitle("Geomodel grid", geoGridPropName(S, V.gridProp).c_str(), float(a), float(b), Colormap::Turbo, false);
        } else if (!S.hz.empty()) {
            double a, b; S.hz[0].range(a, b);
            V.sc.setTitle("Structure", S.hz[0].name + " depth", float(a), float(b), Colormap(V.hzCmap), true);
        } else V.sc.setTitle(title, "", 0, 1, Colormap::Turbo, false);
        V.sc.setCaption(V.domainTime ? "time domain (TWT x v0/2000)" : "depth domain (linear velocity model)");
    }

    // ================================================================ toolbar + canvas
    ImGui::BeginGroup();
    if (ImGui::Button(V.panel ? ICON_PANEL_LEFT "##hide" : ICON_PANEL_LEFT "##show")) V.panel = !V.panel;
    ImGui::SameLine(); if (ImGui::Button(TR("Top"))) V.sc.viewPreset(0);
    ImGui::SameLine(); if (ImGui::Button(TR("South"))) V.sc.viewPreset(1);
    ImGui::SameLine(); if (ImGui::Button(TR("East"))) V.sc.viewPreset(2);
    ImGui::SameLine(); if (ImGui::Button(TR("3D"))) V.sc.viewPreset(3);
    ImGui::SameLine(); if (ImGui::Button(ICON_MAXIMIZE "##fit3")) V.sc.resetCamera();
    ImGui::SameLine();
    if (ImGui::Button(ICON_CAMERA "##png3")) ImGui::OpenPopup("snap3");
    if (ImGui::BeginPopup("snap3")) {
        static int scale = 2;
        ImGui::TextUnformatted("Save the 3-D view as PNG");
        ImGui::RadioButton(TR("1x"), &scale, 1); ImGui::SameLine(); ImGui::RadioButton(TR("2x"), &scale, 2); ImGui::SameLine(); ImGui::RadioButton(TR("3x (print)"), &scale, 3);
        ImGui::TextDisabled("%d x %d pixels", V.sc.width() * scale, V.sc.height() * scale);
        std::string path = "geomodel_3d.png";
        if (ImGui::Button(TR(ICON_DOWNLOAD " Save...")) && fileDialog(true, "PNG image (*.png)\0*.png\0", "png", path)) {
            std::string err;
            if (V.sc.saveImage(path, scale, &err)) notify("3-D view saved: " + path); else notify("Could not save the image: " + err);
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
    ImGui::SameLine(); ImGui::TextDisabled("|"); ImGui::SameLine();
    interpretationBar(S);
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const int w = std::max(32, int(avail.x)), h = std::max(32, int(std::min(avail.y, size.y - ImGui::GetFrameHeightWithSpacing())));
    const ImTextureID srv = V.sc.render(w, h);
    ImGui::InvisibleButton("##geo3d", ImVec2(float(w), float(h)), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle);
    if (srv) ImGui::GetWindowDrawList()->AddImage(srv, p0, ImVec2(p0.x + w, p0.y + h));
    ImGuiIO& io = ImGui::GetIO();
    const float mx = io.MousePos.x - p0.x, my = io.MousePos.y - p0.y;
    if (ImGui::IsItemActive()) {
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left, 3.0f) && !io.KeyShift) V.sc.orbit(io.MouseDelta.x, io.MouseDelta.y);
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Right, 2.0f) || ImGui::IsMouseDragging(ImGuiMouseButton_Middle, 0.0f) || (ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f) && io.KeyShift))
            V.sc.pan(io.MouseDelta.x, io.MouseDelta.y);
    }
    if (ImGui::IsItemHovered()) {
        if (io.MouseWheel != 0) V.sc.zoom(io.MouseWheel);
        if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) V.sc.resetCamera();
        double o[3], d[3];
        if (S.pickMode == 3 && V.sc.pixelRay(mx, my, o, d)) {   // C36-4 horizon brush
            double X, Y, Dp;
            if (rayHitHorizon(S, V, S.editHz, o, d, X, Y, Dp)) {
                // brush footprint draped on the horizon
                ImDrawList* dl = ImGui::GetWindowDrawList();
                ImVec2 prev; bool hp = false;
                const Surface& sf = S.hz[size_t(S.editHz)];
                for (int q = 0; q <= 48; ++q) {
                    const double a = 2 * 3.14159265358979 * q / 48, x = X + S.brushR * std::cos(a), y = Y + S.brushR * std::sin(a), z = sf.at(x, y);
                    float sx, sy;
                    const bool ok = std::isfinite(z) && V.sc.project(float(x - V.ox), float(y - V.oy), float(zOfDepth(S, V, z)), sx, sy);
                    const ImVec2 pt(p0.x + sx, p0.y + sy);
                    if (ok && hp) dl->AddLine(prev, pt, IM_COL32(255, 255, 255, 220), 2.0f);
                    prev = pt; hp = ok;
                }
                ImGui::SetTooltip("%s\nx %.1f  y %.1f\ndepth %.2f\nclick: %s (R %.0f)", sf.name.c_str(), X, Y, Dp, S.brushMode == 0 ? "raise / lower" : S.brushMode == 1 ? "smooth" : "flatten", S.brushR);
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && io.MouseDragMaxDistanceSqr[0] < 16.0f && !io.KeyShift) geoApplyBrush(S, X, Y);
            }
        } else if (hasSeis && V.sc.pixelRay(mx, my, o, d)) {
            int kind; double fi, fj, t, z;
            if (rayHitSlices(S, V, g, o, d, kind, fi, fj, t, z)) {
                const int I = std::clamp(int(std::lround(fi)), 0, g.ni - 1), J = std::clamp(int(std::lround(fj)), 0, g.nj - 1), Sm = std::clamp(int(std::lround((t - g.t0) / g.dt)), 0, g.ns - 1);
                float val = NAN;
                if (kind == 0 && V.secIL.size() >= size_t(g.nj) * size_t(g.ns)) val = V.secIL[size_t(J) * size_t(g.ns) + size_t(Sm)];
                else if (kind == 1 && V.secXL.size() >= size_t(g.ni) * size_t(g.ns)) val = V.secXL[size_t(I) * size_t(g.ns) + size_t(Sm)];
                else if (kind == 2 && !V.tsFromAttr && V.secTS.size() >= size_t(g.ni) * size_t(g.nj)) val = V.secTS[size_t(I) * size_t(g.nj) + size_t(J)];
                double x, y; g.xy(fi, fj, x, y);
                ImGui::SetTooltip("IL %d  XL %d\nx %.1f  y %.1f\nTWT %.1f ms  depth %.1f\n%s %.4g", g.il0 + I * g.dil, g.xl0 + J * g.dxl, x, y, t, timeToDepth(S.vm, t), kAttrNames[S.secAttr], val);
                const bool click = ImGui::IsMouseReleased(ImGuiMouseButton_Left) && io.MouseDragMaxDistanceSqr[0] < 16.0f && !io.KeyShift;
                if (click && kind < 2 && S.pickMode == 1) addSeedAt(S, g, fi, fj, t);
                if (click && kind < 2 && S.pickMode == 2) addStickPoint(S, g, fi, fj, t);
                if (click && kind == 2 && S.pickMode == 1) addSeedAt(S, g, fi, fj, t);
            }
        }
        // keyboard: slice stepping (page up / down: inline, arrows with Ctrl: crossline)
        if (hasSeis && !io.WantTextInput) {
            if (ImGui::IsKeyPressed(ImGuiKey_PageUp)) S.il = std::min(g.ni - 1, S.il + 1);
            if (ImGui::IsKeyPressed(ImGuiKey_PageDown)) S.il = std::max(0, S.il - 1);
            if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_RightArrow)) S.xl = std::min(g.nj - 1, S.xl + 1);
            if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_LeftArrow)) S.xl = std::max(0, S.xl - 1);
        }
    }
    ImGui::EndGroup();
}
