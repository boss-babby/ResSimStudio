// C39 OpenGL 3.3 back end of GeoScene3D (Linux / macOS). Camera, content preparation (geoscene3d_cpu.h) and the overlay
// layout are shared with the Direct3D 11 back end (geoscene3d.cpp).
#include "geoscene3d.h"
#ifdef RESSIM_GL
#include "geoscene3d_cpu.h"
#include "gl_api.h"
#include "gl_view.h"
#include <algorithm>
#include <cmath>
#include <cstring>

using namespace geo3d;

namespace {
const char* kPrelude = R"(
uniform mat4 VP;
uniform mat4 invVP;
uniform vec4 center;      // xyz model centre, w = vertical exaggeration
uniform vec4 eye;         // xyz camera (world)
uniform vec4 light;       // xyz light direction (world)
uniform vec4 vpx;         // width, height, 1/width, 1/height
uniform vec4 range;       // lo, 1/(hi-lo), opacity, contour interval
uniform vec4 flags;       // alpha mode, edges
vec3 toWorld(vec3 p) { return vec3(p.x - center.x, -(p.z - center.z) * center.w, p.y - center.y); }
vec3 nrmWorld(vec3 n) { return normalize(vec3(n.x, -n.z / center.w, n.y)); }
vec4 clipOf(vec3 w) { vec4 c = VP * vec4(w, 1.0); return vec4(c.x, -c.y, 2.0 * c.z - c.w, c.w); }
bool isNaNf(float v) { return (floatBitsToUint(v) & 0x7fffffffu) > 0x7f800000u; }
vec3 shade(vec3 base, vec3 n, vec3 w) {
    vec3 v = normalize(eye.xyz - w);
    if (dot(n, v) < 0.0) n = -n;
    float d1 = clamp(dot(n, normalize(light.xyz)), 0.0, 1.0);
    float d2 = clamp(dot(n, v), 0.0, 1.0);
    vec3 h = normalize(normalize(light.xyz) + v);
    float sp = pow(clamp(dot(n, h), 0.0, 1.0), 40.0) * 0.15;
    return base * (0.32 + 0.45 * d1 + 0.28 * d2) + sp;
}
)";
// ---- value meshes
const char* kMeshVS = R"(
layout(location = 0) in vec3 aPos; layout(location = 1) in vec3 aNrm; layout(location = 2) in float aVal; layout(location = 3) in vec2 aUV;
out vec3 vNrm; out float vVal; out vec2 vUV; out vec3 vW;
void main() { vec3 w = toWorld(aPos); gl_Position = clipOf(w); vNrm = nrmWorld(aNrm); vVal = aVal; vUV = aUV; vW = w; }
)";
const char* kMeshFS = R"(
in vec3 vNrm; in float vVal; in vec2 vUV; in vec3 vW;
uniform sampler1D cmap;
out vec4 outColor;
void main() {
    vec3 base = vec3(0.5, 0.52, 0.55);
    bool nan = isNaNf(vVal);
    if (!nan) base = textureLod(cmap, clamp((vVal - range.x) * range.y, 0.0, 1.0), 0.0).rgb;
    vec3 c = shade(base, normalize(vNrm), vW);
    float f = nan ? 0.0 : vVal / max(range.w, 1e-30);
    float fwf = max(fwidth(f), 1e-6);
    vec2 fw = max(fwidth(vUV), vec2(1e-6));
    if (range.w > 0.0 && !nan) {
        float d = abs(fract(f - 0.5) - 0.5) / fwf;
        c = mix(c, c * 0.30, (1.0 - clamp(d - 0.4, 0.0, 1.0)) * 0.85);
    }
    if (flags.y > 0.5) {
        vec2 a = smoothstep(vec2(0.0), fw * 1.2, vUV) * smoothstep(vec2(0.0), fw * 1.2, 1.0 - vUV);
        c = mix(c, c * 0.3, (1.0 - min(a.x, a.y)) * 0.8);
    }
    outColor = vec4(c, range.z);
}
)";
// ---- textured slices
const char* kTexVS = R"(
layout(location = 0) in vec3 aPos; layout(location = 1) in vec2 aUV;
out vec2 vUV;
void main() { gl_Position = clipOf(toWorld(aPos)); vUV = aUV; }
)";
const char* kTexFS = R"(
in vec2 vUV;
uniform sampler2D stex;
uniform sampler1D cmap;
out vec4 outColor;
void main() {
    float v = texture(stex, vUV).r;
    if (isNaNf(v)) discard;
    float t = clamp((v - range.x) * range.y, 0.0, 1.0);
    vec3 c = textureLod(cmap, t, 0.0).rgb;
    float a = range.z;
    if (flags.x > 1.5) a *= clamp((1.0 - t) * 1.6, 0.0, 1.0);
    else if (flags.x > 0.5) a *= clamp(abs(t - 0.5) * 3.0, 0.0, 1.0);
    if (a < 0.004) discard;
    outColor = vec4(c, a);
}
)";
// ---- lit coloured solids and lines
const char* kColVS = R"(
layout(location = 0) in vec3 aPos; layout(location = 1) in vec3 aNrm; layout(location = 2) in vec4 aCol;
out vec3 vNrm; out vec4 vCol; out vec3 vW;
void main() { vec3 w = toWorld(aPos); gl_Position = clipOf(w); vNrm = (dot(aNrm, aNrm) > 0.01) ? nrmWorld(aNrm) : vec3(0.0); vCol = aCol; vW = w; }
)";
const char* kColFS = R"(
in vec3 vNrm; in vec4 vCol; in vec3 vW;
out vec4 outColor;
void main() {
    if (dot(vNrm, vNrm) < 0.01) { outColor = vCol; return; }
    outColor = vec4(shade(vCol.rgb, normalize(vNrm), vW), vCol.a);
}
)";
// ---- ray-marched volume (full-screen triangle, premultiplied output)
const char* kQuadVS = R"(
void main() {
    vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
    gl_Position = vec4(p * vec2(2.0, -2.0) + vec2(-1.0, 1.0), 0.0, 1.0);
}
)";
const char* kVolFS = R"(
uniform vec4 volN;        // ni, nj, ns, steps
uniform vec4 volP;        // x0, y0, threshold, zMode
uniform vec4 volA;        // inverse plan matrix
uniform vec4 volZ;        // t0, dt, zScale
uniform vec4 volV;        // v0, k, datum
uniform vec4 volB0;       // iMin, jMin, zMin
uniform vec4 volB1;       // iMax, jMax, zMax
uniform sampler1D cmap;
uniform sampler3D vol;
#if MSDEPTH
uniform sampler2DMS sceneDepth;
#else
uniform sampler2D sceneDepth;
#endif
out vec4 outColor;
float zToS(float z) {
    float t;
    if (volP.w < 0.5) t = z / max(volZ.z, 1e-9);
    else {
        float dz = z - volV.z;
        t = abs(volV.y) > 1e-9 ? 2000.0 * log(max(1e-6, 1.0 + volV.y * dz / volV.x)) / volV.y : 2000.0 * dz / volV.x;
    }
    return (t - volZ.x) / volZ.y;
}
void main() {
    vec2 pos = gl_FragCoord.xy;    // rows are stored top-first (flipped clip space), as SV_Position in Direct3D
    vec2 ndc = vec2(pos.x * vpx.z * 2.0 - 1.0, 1.0 - pos.y * vpx.w * 2.0);
    vec4 pn = invVP * vec4(ndc, 0.0, 1.0); pn /= pn.w;
    vec4 pf = invVP * vec4(ndc, 1.0, 1.0); pf /= pf.w;
#if MSDEPTH
    float dz = texelFetch(sceneDepth, ivec2(pos), 0).r;
#else
    float dz = texelFetch(sceneDepth, ivec2(pos), 0).r;
#endif
    vec3 ro = pn.xyz, rd = pf.xyz - pn.xyz;
    float tEnd = 1.0;
    if (dz < 1.0) { vec4 ps = invVP * vec4(ndc, dz, 1.0); ps /= ps.w; tEnd = dot(ps.xyz - ro, rd) / dot(rd, rd); }
    vec3 mo = vec3(ro.x + center.x, ro.z + center.y, center.z - ro.y / center.w);
    vec3 md = vec3(rd.x, rd.z, -rd.y / center.w);
    vec2 q = mo.xy - volP.xy;
    float i0 = volA.x * q.x + volA.y * q.y, di = volA.x * md.x + volA.y * md.y;
    float j0 = volA.z * q.x + volA.w * q.y, dj = volA.z * md.x + volA.w * md.y;
    float t0 = 0.0, t1 = tEnd;
    vec3 o3 = vec3(i0, j0, mo.z), d3 = vec3(di, dj, md.z);
    for (int a = 0; a < 3; ++a) {
        float lo = volB0[a], hi = volB1[a];
        if (abs(d3[a]) < 1e-12) { if (o3[a] < lo || o3[a] > hi) { outColor = vec4(0.0); return; } }
        else {
            float ta = (lo - o3[a]) / d3[a], tb = (hi - o3[a]) / d3[a];
            t0 = max(t0, min(ta, tb)); t1 = min(t1, max(ta, tb));
        }
    }
    if (t1 <= t0) { outColor = vec4(0.0); return; }
    int N = int(volN.w);
    float dt = (t1 - t0) / float(N);
    vec3 C = vec3(0.0); float A = 0.0;
    float corr = 128.0 / max(volN.w, 1.0);
    float jit = fract(sin(dot(pos, vec2(12.9898, 78.233))) * 43758.5453);
    for (int k = 0; k < N; ++k) {
        float t = t0 + (float(k) + jit) * dt;
        vec3 p = o3 + t * d3;
        float s = zToS(p.z);
        vec3 uvw = vec3((s + 0.5) / volN.z, (p.y + 0.5) / volN.y, (p.x + 0.5) / volN.x);
        float v = textureLod(vol, uvw, 0.0).r;
        float tc = clamp((v - range.x) * range.y, 0.0, 1.0);
        float strength = abs(tc - 0.5) * 2.0;
        float al = clamp((strength - volP.z) / max(1.0 - volP.z, 1e-3), 0.0, 1.0);
        al = 1.0 - pow(clamp(1.0 - al * al * range.z, 0.0, 1.0), corr);
        if (al > 0.001) {
            vec3 c = textureLod(cmap, tc, 0.0).rgb;
            C += (1.0 - A) * al * c; A += (1.0 - A) * al;
            if (A > 0.985) break;
        }
    }
    outColor = vec4(C, A);
}
)";

unsigned makeColormapTex(Colormap cm, bool rev) {
    float px[256 * 4];
    for (int i = 0; i < 256; ++i) { float rgb[3]; colormapRGB(cm, i / 255.0f, rgb, rev); px[4 * i] = rgb[0]; px[4 * i + 1] = rgb[1]; px[4 * i + 2] = rgb[2]; px[4 * i + 3] = 1; }
    GLuint t = 0;
    glGenTextures(1, &t);
    glBindTexture(GL_TEXTURE_1D, t);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA32F, 256, 0, GL_RGBA, GL_FLOAT, px);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_1D, 0);
    return t;
}
// VAO with vertex + optional index buffer; layout 0 = VtxM, 1 = VtxT, 2 = VtxC
template <class V> void makeVAO(unsigned& vao, unsigned& vb, unsigned* ib, const std::vector<V>& v, const std::vector<unsigned>* idx, int layout) {
    glGenVertexArrays(1, &vao); glGenBuffers(1, &vb);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vb);
    glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(v.size() * sizeof(V)), v.data(), GL_STATIC_DRAW);
    const GLsizei st = sizeof(V);
    if (layout == 0) {
        glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, st, (void*)0);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, st, (void*)12);
        glEnableVertexAttribArray(2); glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, st, (void*)24);
        glEnableVertexAttribArray(3); glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, st, (void*)28);
    } else if (layout == 1) {
        glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, st, (void*)0);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, st, (void*)12);
    } else {
        glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, st, (void*)0);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, st, (void*)12);
        glEnableVertexAttribArray(2); glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, st, (void*)24);
    }
    if (ib && idx) {
        glGenBuffers(1, ib);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *ib);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, GLsizeiptr(idx->size() * sizeof(unsigned)), idx->data(), GL_STATIC_DRAW);
    }
    glBindVertexArray(0);
}
void del(unsigned& vao, unsigned& vb, unsigned* ib = nullptr) {
    if (ib && *ib) { glDeleteBuffers(1, ib); *ib = 0; }
    if (vb) glDeleteBuffers(1, &vb);
    if (vao) glDeleteVertexArrays(1, &vao);
    vao = vb = 0;
}
void delTex(unsigned& t) { if (t) glDeleteTextures(1, &t); t = 0; }
void setSamplers(unsigned p) {
    glUseProgram(p);
    const int l0 = glv::loc(p, "stex"), l1 = glv::loc(p, "cmap"), l2 = glv::loc(p, "vol"), l3 = glv::loc(p, "sceneDepth");
    if (l0 >= 0) glUniform1i(l0, 0);
    if (l1 >= 0) glUniform1i(l1, 1);
    if (l2 >= 0) glUniform1i(l2, 2);
    if (l3 >= 0) glUniform1i(l3, 3);
}
}  // namespace

// ============================================================================ init / shutdown
bool GeoScene3D::init() {
    msaa_ = glv::maxSamples(4);
    const char* def = msaa_ > 1 ? "#define MSDEPTH 1\n" : "#define MSDEPTH 0\n";
    progMesh_ = glv::program(kMeshVS, kMeshFS, kPrelude, def);
    progTex_ = glv::program(kTexVS, kTexFS, kPrelude, def);
    progC_ = glv::program(kColVS, kColFS, kPrelude, def);
    progVol_ = glv::program(kQuadVS, kVolFS, kPrelude, def);
    for (unsigned p : {progMesh_, progTex_, progC_, progVol_}) if (p) setSamplers(p);
    glUseProgram(0);
    glGenVertexArrays(1, &emptyVao_);
    ready_ = progMesh_ && progTex_ && progC_ && progVol_;
    return ready_;
}

void GeoScene3D::releaseObjects() {
    for (Slice& s : slices_) { del(s.vao, s.vb, &s.ib); delTex(s.tex); delTex(s.cm); }
    for (Mesh& m : meshes_) { del(m.vao, m.vb, &m.ib); delTex(m.cm); }
    for (Solid& s : solids_) del(s.vao, s.vb, &s.ib);
    slices_.clear(); meshes_.clear(); solids_.clear();
    del(lvao_, lvb_); lineCount_ = 0;
}
void GeoScene3D::releaseVolume() { delTex(volTex_); delTex(volCm_); hasVol_ = false; volCmapBuilt_ = Colormap(-1); }

void GeoScene3D::shutdown() {
    releaseObjects(); releaseVolume(); releaseTargets();
    del(boxVao_, boxVb_);
    for (unsigned* p : {&progMesh_, &progTex_, &progC_, &progVol_}) { if (*p) glDeleteProgram(*p); *p = 0; }
    if (emptyVao_) glDeleteVertexArrays(1, &emptyVao_);
    emptyVao_ = 0;
    ready_ = false;
}

void GeoScene3D::clear() { releaseObjects(); releaseVolume(); labels_.clear(); }

void GeoScene3D::releaseTargets() { glv::releaseTargets(tg_); w_ = h_ = 0; }
void GeoScene3D::createTargets(int w, int h) {
    glv::createTargets(tg_, w, h, msaa_, true, true);
    w_ = w; h_ = h;
}

// ============================================================================ content
void GeoScene3D::setSlices(const std::vector<GeoSliceDesc>& in) {
    for (Slice& s : slices_) { del(s.vao, s.vb, &s.ib); delTex(s.tex); delTex(s.cm); }
    slices_.clear();
    for (const GeoSliceDesc& D : in) {
        std::vector<VtxT> v; std::vector<unsigned> idx;
        if (!sliceGeometry(D, v, idx)) continue;
        Slice S;
        makeVAO(S.vao, S.vb, &S.ib, v, &idx, 1);
        S.n = unsigned(idx.size());
        glGenTextures(1, &S.tex);
        glBindTexture(GL_TEXTURE_2D, S.tex);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, D.tw, D.th, 0, GL_RED, GL_FLOAT, D.tex.data());
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
        S.cm = makeColormapTex(D.cmap, D.reverse);
        S.lo = D.lo; S.hi = D.hi; S.opacity = D.opacity; S.alphaMode = D.alphaMode;
        slices_.push_back(S);
    }
}

void GeoScene3D::setMeshes(const std::vector<GeoMeshDesc>& in) {
    for (Mesh& m : meshes_) { del(m.vao, m.vb, &m.ib); delTex(m.cm); }
    meshes_.clear();
    for (const GeoMeshDesc& D : in) {
        std::vector<VtxM> v; float cz[3];
        if (!meshGeometry(D, v, cz)) continue;
        Mesh M;
        makeVAO(M.vao, M.vb, &M.ib, v, &D.idx, 0);
        M.n = unsigned(D.idx.size());
        M.cm = makeColormapTex(D.cmap, D.reverse);
        M.lo = D.lo; M.hi = D.hi; M.opacity = D.opacity; M.contour = D.contour; M.edges = D.edges && D.uv.size() >= 2 * v.size();
        for (int m = 0; m < 3; ++m) M.cz[m] = cz[m];
        meshes_.push_back(M);
    }
}

void GeoScene3D::setSolids(const std::vector<GeoSolidDesc>& in) {
    for (Solid& s : solids_) del(s.vao, s.vb, &s.ib);
    solids_.clear();
    for (const GeoSolidDesc& D : in) {
        std::vector<VtxC> v;
        if (!solidGeometry(D, v)) continue;
        Solid S;
        makeVAO(S.vao, S.vb, &S.ib, v, &D.idx, 2);
        S.n = unsigned(D.idx.size()); S.transparent = D.transparent;
        solids_.push_back(S);
    }
}

void GeoScene3D::setLines(const std::vector<GeoLineVtx>& l) {
    del(lvao_, lvb_); lineCount_ = 0;
    if (l.size() < 2) return;
    const std::vector<VtxC> v = lineVertices(l);
    makeVAO(lvao_, lvb_, nullptr, v, nullptr, 2);
    lineCount_ = unsigned(v.size() & ~size_t(1));
}

void GeoScene3D::setVolume(const GeoVolumeDesc* v) {
    releaseVolume();
    if (!v || v->ni < 2 || v->nj < 2 || v->ns < 2 || v->data.size() < size_t(v->ni) * size_t(v->nj) * size_t(v->ns)) return;
    while (glGetError() != GL_NO_ERROR) {}
    glGenTextures(1, &volTex_);
    glBindTexture(GL_TEXTURE_3D, volTex_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, v->ns, v->nj, v->ni, 0, GL_RED, GL_FLOAT, v->data.data());
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    for (GLenum w : {GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T, GL_TEXTURE_WRAP_R}) glTexParameteri(GL_TEXTURE_3D, w, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_3D, 0);
    if (glGetError() != GL_NO_ERROR) { delTex(volTex_); return; }   // e.g. larger than GL_MAX_3D_TEXTURE_SIZE
    vol_ = *v;
    vol_.data.clear(); vol_.data.shrink_to_fit();   // the GPU keeps the samples
    hasVol_ = volTex_ != 0;
}

// ============================================================================ render
ImTextureID GeoScene3D::render(int w, int h) {
    w = std::max(16, w); h = std::max(16, h);
    if (w != w_ || h != h_ || !tg_.fbo) createTargets(w, h);
    if (!tg_.msFbo || !ready_) return ImTextureID(0);
    extent_ = std::max({hi_[0] - lo_[0], hi_[1] - lo_[1], (hi_[2] - lo_[2]) * zScale, 1.0f});
    if (builtBoxZ_ != zScale || !boxVao_) {
        del(boxVao_, boxVb_); boxCount_ = 0;
        const std::vector<VtxC> b = boxLines(lo_, hi_);
        makeVAO(boxVao_, boxVb_, nullptr, b, nullptr, 2);
        boxCount_ = unsigned(b.size());
        builtBoxZ_ = zScale;
    }
    glBindFramebuffer(GL_FRAMEBUFFER, tg_.msFbo);
    glViewport(0, 0, w, h);
    glDisable(GL_SCISSOR_TEST); glDisable(GL_CULL_FACE); glDisable(GL_STENCIL_TEST);
    glClearColor(background[0], background[1], background[2], 1.0f);
    glClearDepth(1.0);
    glDepthMask(GL_TRUE);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST); glDepthFunc(GL_LEQUAL);
    glEnable(GL_BLEND); glBlendEquation(GL_FUNC_ADD);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

    const auto M = viewProj();
    float inv[16];
    if (!invert4(M.data(), inv)) std::memset(inv, 0, sizeof(inv));
    const float cen[4] = {center_[0], center_[1], center_[2], zScale};
    float eye[4] = {0, 0, 0, 0}; eyePos(eye);
    const float light[4] = {0.4f, 0.85f, 0.35f, 0};
    const float vpx[4] = {float(w), float(h), 1.0f / float(w), 1.0f / float(h)};
    for (unsigned p : {progMesh_, progTex_, progC_, progVol_}) {
        glUseProgram(p);
        glUniformMatrix4fv(glv::loc(p, "VP"), 1, GL_TRUE, M.data());
        glUniformMatrix4fv(glv::loc(p, "invVP"), 1, GL_TRUE, inv);
        glUniform4fv(glv::loc(p, "center"), 1, cen);
        glUniform4fv(glv::loc(p, "eye"), 1, eye);
        glUniform4fv(glv::loc(p, "light"), 1, light);
        glUniform4fv(glv::loc(p, "vpx"), 1, vpx);
    }
    auto setObj = [&](unsigned p, const float range[4], const float flags[4]) {
        glUniform4fv(glv::loc(p, "range"), 1, range);
        glUniform4fv(glv::loc(p, "flags"), 1, flags);
    };
    auto drawMesh = [&](const Mesh& m) {
        glUseProgram(progMesh_);
        const float r[4] = {m.lo, m.hi > m.lo ? 1.0f / (m.hi - m.lo) : 0.0f, m.opacity, m.contour}, f[4] = {0, m.edges ? 1.0f : 0.0f, 0, 0};
        setObj(progMesh_, r, f);
        glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_1D, m.cm);
        glBindVertexArray(m.vao);
        glDrawElements(GL_TRIANGLES, GLsizei(m.n), GL_UNSIGNED_INT, nullptr);
    };
    auto drawSlice = [&](const Slice& s) {
        glUseProgram(progTex_);
        const float r[4] = {s.lo, s.hi > s.lo ? 1.0f / (s.hi - s.lo) : 0.0f, s.opacity, 0}, f[4] = {float(s.alphaMode), 0, 0, 0};
        setObj(progTex_, r, f);
        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, s.tex);
        glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_1D, s.cm);
        glBindVertexArray(s.vao);
        glDrawElements(GL_TRIANGLES, GLsizei(s.n), GL_UNSIGNED_INT, nullptr);
    };
    auto drawSolid = [&](const Solid& s) {
        glUseProgram(progC_);
        glBindVertexArray(s.vao);
        glDrawElements(GL_TRIANGLES, GLsizei(s.n), GL_UNSIGNED_INT, nullptr);
    };
    auto drawLines = [&](unsigned vao, unsigned n) {
        if (!vao || n < 2) return;
        glUseProgram(progC_);
        glBindVertexArray(vao);
        glDrawArrays(GL_LINES, 0, GLsizei(n));
    };
    // 1) opaque geometry
    for (const Mesh& m : meshes_) if (m.opacity >= 0.999f && m.vao) drawMesh(m);
    for (const Solid& s : solids_) if (!s.transparent && s.vao) drawSolid(s);
    for (const Slice& s : slices_) if (s.opacity >= 0.999f && s.alphaMode == 0 && s.vao && s.tex) drawSlice(s);
    drawLines(lvao_, lineCount_);
    if (showBox) drawLines(boxVao_, boxCount_);
    // 2) volume: ray-marched against the opaque depth (colour-only framebuffer: the depth texture is sampled)
    if (hasVol_ && volTex_ && tg_.msDepth) {
        if (volCmapBuilt_ != vol_.cmap || volRevBuilt_ != vol_.reverse || !volCm_) {
            delTex(volCm_);
            volCm_ = makeColormapTex(vol_.cmap, vol_.reverse);
            volCmapBuilt_ = vol_.cmap; volRevBuilt_ = vol_.reverse;
        }
        const ObjCB o = volumeConstants(vol_);
        glBindFramebuffer(GL_FRAMEBUFFER, tg_.msFboColorOnly);
        glDisable(GL_DEPTH_TEST);
        glBlendFuncSeparate(GL_ONE, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
        glUseProgram(progVol_);
        glUniform4fv(glv::loc(progVol_, "range"), 1, o.range);
        glUniform4fv(glv::loc(progVol_, "volN"), 1, o.volN);
        glUniform4fv(glv::loc(progVol_, "volP"), 1, o.volP);
        glUniform4fv(glv::loc(progVol_, "volA"), 1, o.volA);
        glUniform4fv(glv::loc(progVol_, "volZ"), 1, o.volZ);
        glUniform4fv(glv::loc(progVol_, "volV"), 1, o.volV);
        glUniform4fv(glv::loc(progVol_, "volB0"), 1, o.volB0);
        glUniform4fv(glv::loc(progVol_, "volB1"), 1, o.volB1);
        glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_1D, volCm_);
        glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_3D, volTex_);
        glActiveTexture(GL_TEXTURE3); glBindTexture(tg_.msaa > 1 ? GL_TEXTURE_2D_MULTISAMPLE : GL_TEXTURE_2D, tg_.msDepth);
        glBindVertexArray(emptyVao_);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glActiveTexture(GL_TEXTURE3); glBindTexture(tg_.msaa > 1 ? GL_TEXTURE_2D_MULTISAMPLE : GL_TEXTURE_2D, 0);
        glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_3D, 0);
        glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_DEPTH_TEST);
        glBindFramebuffer(GL_FRAMEBUFFER, tg_.msFbo);
    }
    // 3) transparent geometry (depth tested, not written), far objects first
    glDepthMask(GL_FALSE);
    std::vector<std::pair<float, int>> order;
    for (int k = 0; k < int(meshes_.size()); ++k) if (meshes_[size_t(k)].opacity < 0.999f && meshes_[size_t(k)].vao) {
        float wc[3]; worldOf(meshes_[size_t(k)].cz[0], meshes_[size_t(k)].cz[1], meshes_[size_t(k)].cz[2], wc);
        order.push_back({-((wc[0] - eye[0]) * (wc[0] - eye[0]) + (wc[1] - eye[1]) * (wc[1] - eye[1]) + (wc[2] - eye[2]) * (wc[2] - eye[2])), k});
    }
    std::sort(order.begin(), order.end());
    for (const auto& pr : order) drawMesh(meshes_[size_t(pr.second)]);
    for (const Slice& s : slices_) if ((s.opacity < 0.999f || s.alphaMode != 0) && s.vao && s.tex) drawSlice(s);
    for (const Solid& s : solids_) if (s.transparent && s.vao) drawSolid(s);
    glDepthMask(GL_TRUE);
    glBindVertexArray(0);
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_1D, 0);
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, 0);
    glv::resolve(tg_);
    glv::endView();
    drawOverlay();
    glv::endView();
    return (ImTextureID)(intptr_t)tg_.color;
}

void GeoScene3D::drawOverlay() {
    const float s = g_gfx.dpiScale;
    glv::renderOverlay(tg_, overlayScale_, 11.5f * s, 12.5f * s, 16.0f * s, [this](OverlayPainter& p) { paintOverlay(p); });
}

bool GeoScene3D::saveImage(const std::string& path, int scale, std::string* err) {
    const int w0 = w_, h0 = h_;
    if (w0 <= 0 || h0 <= 0) { if (err) *err = "view not rendered yet"; return false; }
    scale = std::max(1, std::min(4, scale));
    overlayScale_ = float(scale);
    render(w0 * scale, h0 * scale);
    overlayScale_ = 1.0f;
    const bool ok = glv::savePng(tg_, path, err);
    render(w0, h0);
    return ok;
}
#endif  // RESSIM_GL
