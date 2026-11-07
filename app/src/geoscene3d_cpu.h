// C39: CPU side of GeoScene3D shared by the Direct3D 11 (geoscene3d.cpp) and OpenGL (geoscene3d_gl.cpp) back ends:
// vertex layouts, mesh / slice / solid vertex building, the bounding box and the volume ray-march parameters.
#pragma once
#include "geoscene3d.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace geo3d {

struct VtxM { float p[3]; float n[3]; float v; float uv[2]; };
struct VtxT { float p[3]; float uv[2]; };
struct VtxC { float p[3]; float n[3]; float c[4]; };
// per-object constants (D3D cbuffer Obj layout; the GL back end uploads the same values as uniforms)
struct ObjCB { float range[4]; float flags[4]; float volN[4]; float volP[4]; float volA[4]; float volZ[4]; float volV[4]; float volB0[4]; float volB1[4]; };

inline void norm3(float v[3]) { const float l = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); if (l > 0) { v[0] /= l; v[1] /= l; v[2] /= l; } }
inline void cross(const float a[3], const float b[3], float o[3]) { o[0] = a[1] * b[2] - a[2] * b[1]; o[1] = a[2] * b[0] - a[0] * b[2]; o[2] = a[0] * b[1] - a[1] * b[0]; }
bool invert4(const float m[16], float out[16]);

// slice lattice -> vertices + triangle indices; false if the description is invalid
inline bool sliceGeometry(const GeoSliceDesc& D, std::vector<VtxT>& v, std::vector<unsigned>& idx) {
    if (D.nu < 2 || D.nv < 2 || D.pos.size() < size_t(D.nu) * size_t(D.nv) * 3 || D.tw < 1 || D.th < 1 || D.tex.size() < size_t(D.tw) * size_t(D.th)) return false;
    v.assign(size_t(D.nu) * size_t(D.nv), VtxT{});
    for (int b = 0; b < D.nv; ++b)
        for (int a = 0; a < D.nu; ++a) {
            const size_t k = size_t(b) * size_t(D.nu) + size_t(a);
            v[k] = VtxT{{D.pos[3 * k], D.pos[3 * k + 1], D.pos[3 * k + 2]}, {float(a) / float(D.nu - 1), float(b) / float(D.nv - 1)}};
        }
    idx.clear();
    for (int b = 0; b + 1 < D.nv; ++b)
        for (int a = 0; a + 1 < D.nu; ++a) {
            const unsigned p = unsigned(b * D.nu + a);
            idx.insert(idx.end(), {p, p + 1, p + unsigned(D.nu), p + 1, p + unsigned(D.nu) + 1, p + unsigned(D.nu)});
        }
    return true;
}

// value mesh -> vertices (area-weighted normals when none are given) and centroid; false if invalid
inline bool meshGeometry(const GeoMeshDesc& D, std::vector<VtxM>& v, float cz[3]) {
    const size_t n = D.pos.size() / 3;
    if (n < 3 || D.idx.size() < 3) return false;
    v.assign(n, VtxM{});
    for (size_t k = 0; k < n; ++k) {
        v[k] = VtxM{{D.pos[3 * k], D.pos[3 * k + 1], D.pos[3 * k + 2]}, {0, 0, 0}, k < D.val.size() ? D.val[k] : NAN, {0, 0}};
        if (D.nrm.size() >= 3 * n) { v[k].n[0] = D.nrm[3 * k]; v[k].n[1] = D.nrm[3 * k + 1]; v[k].n[2] = D.nrm[3 * k + 2]; }
        if (D.uv.size() >= 2 * n) { v[k].uv[0] = D.uv[2 * k]; v[k].uv[1] = D.uv[2 * k + 1]; }
    }
    if (D.nrm.size() < 3 * n) {   // area-weighted vertex normals
        for (size_t t = 0; t + 2 < D.idx.size(); t += 3) {
            const unsigned a = D.idx[t], b = D.idx[t + 1], c = D.idx[t + 2];
            if (a >= n || b >= n || c >= n) continue;
            const float e1[3] = {v[b].p[0] - v[a].p[0], v[b].p[1] - v[a].p[1], v[b].p[2] - v[a].p[2]};
            const float e2[3] = {v[c].p[0] - v[a].p[0], v[c].p[1] - v[a].p[1], v[c].p[2] - v[a].p[2]};
            float nn[3]; cross(e1, e2, nn);
            for (unsigned q : {a, b, c}) for (int m = 0; m < 3; ++m) v[q].n[m] += nn[m];
        }
        for (auto& x : v) norm3(x.n);
    }
    cz[0] = cz[1] = cz[2] = 0;
    for (size_t k = 0; k < n; ++k) for (int m = 0; m < 3; ++m) cz[m] += v[k].p[m] / float(n);
    return true;
}

// lit solid -> vertices; false if invalid
inline bool solidGeometry(const GeoSolidDesc& D, std::vector<VtxC>& v) {
    const size_t n = D.pos.size() / 3;
    if (n < 3 || D.idx.size() < 3) return false;
    v.assign(n, VtxC{});
    for (size_t k = 0; k < n; ++k) {
        v[k] = VtxC{{D.pos[3 * k], D.pos[3 * k + 1], D.pos[3 * k + 2]}, {0, 0, 0}, {0.7f, 0.7f, 0.7f, 1}};
        if (D.col.size() >= 4 * n) for (int m = 0; m < 4; ++m) v[k].c[m] = D.col[4 * k + size_t(m)];
        if (D.nrm.size() >= 3 * n) for (int m = 0; m < 3; ++m) v[k].n[m] = D.nrm[3 * k + size_t(m)];
    }
    if (D.nrm.size() < 3 * n) {
        for (size_t t = 0; t + 2 < D.idx.size(); t += 3) {
            const unsigned a = D.idx[t], b = D.idx[t + 1], c = D.idx[t + 2];
            if (a >= n || b >= n || c >= n) continue;
            const float e1[3] = {v[b].p[0] - v[a].p[0], v[b].p[1] - v[a].p[1], v[b].p[2] - v[a].p[2]};
            const float e2[3] = {v[c].p[0] - v[a].p[0], v[c].p[1] - v[a].p[1], v[c].p[2] - v[a].p[2]};
            float nn[3]; cross(e1, e2, nn);
            for (unsigned q : {a, b, c}) for (int m = 0; m < 3; ++m) v[q].n[m] += nn[m];
        }
        for (auto& x : v) { norm3(x.n); if (x.n[0] == 0 && x.n[1] == 0 && x.n[2] == 0) x.n[2] = 1e-3f; }
    }
    return true;
}

inline std::vector<VtxC> lineVertices(const std::vector<GeoLineVtx>& l) {
    std::vector<VtxC> v(l.size());
    for (size_t k = 0; k < l.size(); ++k) v[k] = VtxC{{l[k].x, l[k].y, l[k].z}, {0, 0, 0}, {l[k].r, l[k].g, l[k].b, l[k].a}};
    return v;
}

// bounding box lines of the model box
inline std::vector<VtxC> boxLines(const float l[3], const float u[3]) {
    std::vector<VtxC> b;
    auto L = [&](float x0, float y0, float z0, float x1, float y1, float z1, float a) {
        b.push_back(VtxC{{x0, y0, z0}, {0, 0, 0}, {0.75f, 0.78f, 0.85f, a}}); b.push_back(VtxC{{x1, y1, z1}, {0, 0, 0}, {0.75f, 0.78f, 0.85f, a}});
    };
    for (int k = 0; k < 2; ++k) {
        const float z = k ? u[2] : l[2];
        const float a = k ? 0.55f : 0.35f;
        L(l[0], l[1], z, u[0], l[1], z, a); L(u[0], l[1], z, u[0], u[1], z, a); L(u[0], u[1], z, l[0], u[1], z, a); L(l[0], u[1], z, l[0], l[1], z, a);
    }
    L(l[0], l[1], l[2], l[0], l[1], u[2], 0.45f); L(u[0], l[1], l[2], u[0], l[1], u[2], 0.45f); L(u[0], u[1], l[2], u[0], u[1], u[2], 0.45f); L(l[0], u[1], l[2], l[0], u[1], u[2], 0.45f);
    return b;
}

// ray-march constants of the seismic volume
inline ObjCB volumeConstants(const GeoVolumeDesc& vol) {
    ObjCB o = {};
    o.range[0] = vol.lo; o.range[1] = vol.hi > vol.lo ? 1.0f / (vol.hi - vol.lo) : 0.0f; o.range[2] = vol.opacity;
    o.volN[0] = float(vol.ni); o.volN[1] = float(vol.nj); o.volN[2] = float(vol.ns); o.volN[3] = float(std::clamp(vol.steps, 32, 2048));
    o.volP[0] = float(vol.x0); o.volP[1] = float(vol.y0); o.volP[2] = std::clamp(vol.threshold, 0.0f, 0.99f); o.volP[3] = float(vol.zMode);
    const double det = vol.ix * vol.jy - vol.iy * vol.jx;
    if (std::fabs(det) > 1e-12) { o.volA[0] = float(vol.jy / det); o.volA[1] = float(-vol.jx / det); o.volA[2] = float(-vol.iy / det); o.volA[3] = float(vol.ix / det); }
    o.volZ[0] = float(vol.t0); o.volZ[1] = float(vol.dt); o.volZ[2] = float(vol.zScale);
    o.volV[0] = float(vol.v0); o.volV[1] = float(vol.k); o.volV[2] = float(vol.datum);
    auto zOfS = [&](double s) {
        const double t = vol.t0 + s * vol.dt;
        if (vol.zMode == 0) return t * vol.zScale;
        const double tt = t / 2000.0;
        return vol.datum + (std::fabs(vol.k) > 1e-12 ? vol.v0 * (std::exp(vol.k * tt) - 1) / vol.k : vol.v0 * tt);
    };
    const int i0 = std::clamp(vol.iMin, 0, vol.ni - 1), i1 = std::clamp(vol.iMax, i0, vol.ni - 1);
    const int j0 = std::clamp(vol.jMin, 0, vol.nj - 1), j1 = std::clamp(vol.jMax, j0, vol.nj - 1);
    const int s0 = std::clamp(vol.sMin, 0, vol.ns - 1), s1 = std::clamp(vol.sMax, s0, vol.ns - 1);
    o.volB0[0] = float(i0); o.volB0[1] = float(j0); o.volB0[2] = float(zOfS(s0));
    o.volB1[0] = float(i1); o.volB1[1] = float(j1); o.volB1[2] = float(zOfS(s1));
    return o;
}

}  // namespace geo3d
