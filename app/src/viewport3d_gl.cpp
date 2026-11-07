// C39 OpenGL 3.3 back end of Viewport3D (Linux / macOS). Geometry, camera, picking, per-frame preparation and the overlay
// layout are shared with the Direct3D 11 back end (viewport3d.cpp); this file holds the GPU objects, shaders and drawing.
#include "viewport3d.h"
#include "theme.h"
#ifdef RESSIM_GL
#include "gl_api.h"
#include "gl_view.h"
#include <algorithm>
#include <cmath>
#include <cstring>

namespace {
// GLSL port of the HLSL viewport shader (viewport3d.cpp). Clip space: y flipped (top row first, as Direct3D) and the
// Direct3D depth range [0, 1] mapped to OpenGL's [-1, 1].
const char* kPrelude = R"(
uniform mat4 VP;
uniform vec4 center;      // xyz model centre, w = vertical exaggeration
uniform vec4 range;       // x = vmin, y = 1/(vmax-vmin), z = edge strength, w = log scale flag
uniform vec4 lightDir;    // xyz light direction, w = discrete colour bands (0 = continuous)
uniform vec4 eye;         // xyz camera position (world), w = cell opacity
vec3 toWorld(vec3 p) { return vec3(p.x - center.x, -(p.z - center.z) * center.w, p.y - center.y); }
vec3 nrmWorld(vec3 n) { return normalize(vec3(n.x, -n.z / center.w, n.y)); }
vec4 clipOf(vec3 w) { vec4 c = VP * vec4(w, 1.0); return vec4(c.x, -c.y, 2.0 * c.z - c.w, c.w); }
vec3 shade(vec3 base, vec3 n, vec3 w) {
    vec3 v = normalize(eye.xyz - w);
    if (dot(n, v) < 0.0) n = -n;
    float d1 = clamp(dot(n, normalize(lightDir.xyz)), 0.0, 1.0);
    float d2 = clamp(dot(n, v), 0.0, 1.0);
    vec3 h = normalize(normalize(lightDir.xyz) + v);
    float sp = pow(clamp(dot(n, h), 0.0, 1.0), 48.0) * 0.12;
    return base * (0.30 + 0.45 * d1 + 0.30 * d2) + sp;
}
)";
const char* kCellVS = R"(
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNrm;
layout(location = 2) in vec2 aUV;
layout(location = 3) in uint aCell;
uniform samplerBuffer vals;
out vec3 vNrm; out vec2 vUV; out float vVal; out vec3 vW;
void main() {
    vec3 w = toWorld(aPos);
    gl_Position = clipOf(w);
    vNrm = nrmWorld(aNrm); vUV = aUV; vVal = texelFetch(vals, int(aCell)).r; vW = w;
}
)";
const char* kCellFS = R"(
in vec3 vNrm; in vec2 vUV; in float vVal; in vec3 vW;
uniform sampler1D cmap;
out vec4 outColor;
void main() {
    vec3 base;
    if (vVal < -1e29) base = vec3(0.45, 0.47, 0.50);
    else {
        float v = vVal;
        if (range.w > 0.5) v = log(max(v, 1e-12)) / log(10.0);
        float t = clamp((v - range.x) * range.y, 0.0, 1.0);
        if (lightDir.w > 1.5) t = (min(floor(t * lightDir.w), lightDir.w - 1.0) + 0.5) / lightDir.w;
        base = textureLod(cmap, t, 0.0).rgb;
    }
    vec3 c = shade(base, normalize(vNrm), vW);
    vec2 fw = max(fwidth(vUV), vec2(1e-6));
    vec2 a = smoothstep(vec2(0.0), fw * 1.3, vUV) * smoothstep(vec2(0.0), fw * 1.3, 1.0 - vUV);
    float edge = (1.0 - min(a.x, a.y)) * range.z;
    c = mix(c, c * 0.25, edge);
    outColor = vec4(c, clamp(eye.w + edge * 0.35 * (1.0 - eye.w), 0.0, 1.0));
}
)";
const char* kColVS = R"(
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNrm;
layout(location = 2) in vec4 aCol;
out vec3 vNrm; out vec4 vCol; out vec3 vW;
void main() {
    vec3 w = toWorld(aPos);
    gl_Position = clipOf(w);
    vNrm = (dot(aNrm, aNrm) > 0.01) ? nrmWorld(aNrm) : vec3(0.0);
    vCol = aCol; vW = w;
}
)";
const char* kColFS = R"(
in vec3 vNrm; in vec4 vCol; in vec3 vW;
out vec4 outColor;
void main() {
    if (dot(vNrm, vNrm) < 0.01) { outColor = vCol; return; }
    outColor = vec4(shade(vCol.rgb, normalize(vNrm), vW), vCol.a);
}
)";

template <class V> void uploadVB(unsigned& vao, unsigned& vb, const std::vector<V>& v, bool colourLayout) {
    if (!vao) glGenVertexArrays(1, &vao);
    if (!vb) glGenBuffers(1, &vb);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vb);
    glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(v.size() * sizeof(V)), v.data(), GL_STATIC_DRAW);
    if (colourLayout) {   // VtxC: p3 n3 c4
        glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(V), (void*)0);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(V), (void*)12);
        glEnableVertexAttribArray(2); glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(V), (void*)24);
    }
}
void delVAO(unsigned& vao, unsigned& vb) {
    if (vb) glDeleteBuffers(1, &vb);
    if (vao) glDeleteVertexArrays(1, &vao);
    vao = vb = 0;
}
void delBuf(unsigned& b) { if (b) glDeleteBuffers(1, &b); b = 0; }
}  // namespace

bool Viewport3D::init() {
    prog_ = glv::program(kCellVS, kCellFS, kPrelude);
    progC_ = glv::program(kColVS, kColFS, kPrelude);
    msaa_ = glv::maxSamples(4);
    if (prog_) {
        glUseProgram(prog_);
        glUniform1i(glv::loc(prog_, "vals"), 0);
        glUniform1i(glv::loc(prog_, "cmap"), 1);
        glUseProgram(0);
    }
    return prog_ && progC_;
}

void Viewport3D::shutdown() {
    releaseTargets();
    gpuReleaseCells(); gpuReleaseWells(); gpuReleasePanels();
    delVAO(slvao_, slvb_); slCount_ = 0;
    if (valTex_) glDeleteTextures(1, &valTex_);
    delBuf(valBuf_);
    if (cmapTex_) glDeleteTextures(1, &cmapTex_);
    valTex_ = cmapTex_ = 0; valCap_ = 0;
    if (prog_) glDeleteProgram(prog_);
    if (progC_) glDeleteProgram(progC_);
    prog_ = progC_ = 0;
}

void Viewport3D::releaseTargets() { glv::releaseTargets(tg_); w_ = h_ = 0; }
void Viewport3D::createTargets(int w, int h) {
    glv::createTargets(tg_, w, h, msaa_, false, false);
    w_ = w; h_ = h;
}

// ----------------------------------------------------------------------------- uploads
void Viewport3D::gpuReleaseCells() {
    if (ib_) glDeleteBuffers(1, &ib_);
    ib_ = 0;
    delVAO(vao_, vb_); delVAO(lvao_, lvb_);
    indexCount_ = 0; lineCount_ = 0;
}
void Viewport3D::gpuCells(const std::vector<Vtx>& vtx, const std::vector<unsigned>& idx, const std::vector<VtxC>& lv) {
    if (!vtx.empty()) {
        glGenVertexArrays(1, &vao_); glGenBuffers(1, &vb_); glGenBuffers(1, &ib_);
        glBindVertexArray(vao_);
        glBindBuffer(GL_ARRAY_BUFFER, vb_);
        glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(vtx.size() * sizeof(Vtx)), vtx.data(), GL_STATIC_DRAW);
        glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vtx), (void*)0);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vtx), (void*)12);
        glEnableVertexAttribArray(2); glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vtx), (void*)24);
        glEnableVertexAttribArray(3); glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(Vtx), (void*)32);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ib_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, GLsizeiptr(idx.size() * sizeof(unsigned)), idx.data(), GL_DYNAMIC_DRAW);
        glBindVertexArray(0);
        indexCount_ = unsigned(idx.size());
    }
    if (!lv.empty()) {
        uploadVB(lvao_, lvb_, lv, true);
        glBindVertexArray(0);
        lineCount_ = unsigned(lv.size());
    }
}
void Viewport3D::gpuReleaseWells() {
    if (wib_) glDeleteBuffers(1, &wib_);
    wib_ = 0;
    delVAO(wvao_, wvb_);
    wIndexCount_ = 0;
}
void Viewport3D::gpuWells(const std::vector<VtxC>& v, const std::vector<unsigned>& id) {
    uploadVB(wvao_, wvb_, v, true);
    glGenBuffers(1, &wib_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, wib_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, GLsizeiptr(id.size() * sizeof(unsigned)), id.data(), GL_STATIC_DRAW);
    glBindVertexArray(0);
    wIndexCount_ = unsigned(id.size());
}
void Viewport3D::gpuReleasePanels() { delVAO(pvao_, pvb_); delVAO(polvao_, polvb_); pCount_ = polCount_ = 0; }
void Viewport3D::gpuPanels(const std::vector<VtxC>& tv, const std::vector<VtxC>& lv) {
    if (!tv.empty()) { uploadVB(pvao_, pvb_, tv, true); pCount_ = unsigned(tv.size()); }
    if (!lv.empty()) { uploadVB(polvao_, polvb_, lv, true); polCount_ = unsigned(lv.size()); }
    glBindVertexArray(0);
}
void Viewport3D::gpuStreamlines(const std::vector<VtxC>& lv) {
    delVAO(slvao_, slvb_); slCount_ = 0;
    if (lv.empty()) return;
    uploadVB(slvao_, slvb_, lv, true);
    glBindVertexArray(0);
    slCount_ = unsigned(lv.size());
}
void Viewport3D::gpuValues(const std::vector<float>& v) {
    const int na = int(v.size());
    if (!valBuf_) { glGenBuffers(1, &valBuf_); glGenTextures(1, &valTex_); }
    glBindBuffer(GL_TEXTURE_BUFFER, valBuf_);
    if (na > valCap_) {
        glBufferData(GL_TEXTURE_BUFFER, GLsizeiptr(std::max(1, na) * sizeof(float)), nullptr, GL_DYNAMIC_DRAW);
        valCap_ = std::max(1, na);
        glBindTexture(GL_TEXTURE_BUFFER, valTex_);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, valBuf_);
        glBindTexture(GL_TEXTURE_BUFFER, 0);
    }
    if (na > 0) glBufferSubData(GL_TEXTURE_BUFFER, 0, GLsizeiptr(na * sizeof(float)), v.data());
    glBindBuffer(GL_TEXTURE_BUFFER, 0);
}
bool Viewport3D::gpuHasColormap() const { return cmapTex_ != 0; }
void Viewport3D::gpuColormap(const unsigned char* px) {
    if (!cmapTex_) glGenTextures(1, &cmapTex_);
    glBindTexture(GL_TEXTURE_1D, cmapTex_);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA8, 256, 0, GL_RGBA, GL_UNSIGNED_BYTE, px);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_1D, 0);
}
bool Viewport3D::gpuHasIndexBuffer() const { return ib_ != 0; }
void Viewport3D::gpuIndexUpdate() {
    glBindVertexArray(vao_);   // the element buffer binding is VAO state
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, GLsizeiptr(sortedIdx_.size() * sizeof(unsigned)), sortedIdx_.data());
    glBindVertexArray(0);
}

// ----------------------------------------------------------------------------- rendering
ImTextureID Viewport3D::render(int w, int h) {
    w = std::max(16, w); h = std::max(16, h);
    if (w != w_ || h != h_ || !tg_.fbo) createTargets(w, h);
    if (!tg_.fbo || !prog_) return ImTextureID(0);
    prepareFrame();

    glBindFramebuffer(GL_FRAMEBUFFER, tg_.msFbo);
    glViewport(0, 0, w, h);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_STENCIL_TEST);
    { const float* b3 = background3d(); glClearColor(b3[0], b3[1], b3[2], 1.0f); }   // C40 setting
    glClearDepth(1.0);
    glDepthMask(GL_TRUE);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

    // constants (same values as the D3D constant buffer)
    const auto M = viewProj();
    float lo = vmin_, hi = vmax_;
    if (set_.logScale) { lo = std::log10(std::max(vmin_, 1e-12f)); hi = std::log10(std::max(vmax_, 1e-12f)); }
    const float range[4] = {lo, (hi > lo) ? 1.0f / (hi - lo) : 0.0f, set_.edges ? 0.85f : 0.0f, set_.logScale ? 1.0f : 0.0f};
    const float light[4] = {0.4f, 0.85f, 0.35f, float(std::max(0, set_.bands))};
    const float d = dist_ * extent_;
    const float opacity = std::min(1.0f, std::max(0.05f, set_.opacity));
    const float eye[4] = {panX_ + d * std::cos(pitch_) * std::sin(yaw_), panY_ + d * std::sin(pitch_), panZ_ + d * std::cos(pitch_) * std::cos(yaw_), opacity};
    const float cen[4] = {center_[0], center_[1], center_[2], set_.zScale};
    for (unsigned p : {prog_, progC_}) {
        glUseProgram(p);
        glUniformMatrix4fv(glv::loc(p, "VP"), 1, GL_TRUE, M.data());
        glUniform4fv(glv::loc(p, "center"), 1, cen);
        glUniform4fv(glv::loc(p, "range"), 1, range);
        glUniform4fv(glv::loc(p, "lightDir"), 1, light);
        glUniform4fv(glv::loc(p, "eye"), 1, eye);
    }
    const bool transparent = opacity < 0.999f;
    auto drawCells = [&]() {
        if (!(set_.showCells && vao_ && indexCount_)) return;
        if (transparent) { sortFaces(false); glDepthMask(GL_FALSE); }
        glUseProgram(prog_);
        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_BUFFER, valTex_);
        glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_1D, cmapTex_);
        glBindVertexArray(vao_);
        glDrawElements(GL_TRIANGLES, GLsizei(indexCount_), GL_UNSIGNED_INT, nullptr);
        glDepthMask(GL_TRUE);
    };
    auto noDepth = [&](bool off) { if (off) glDisable(GL_DEPTH_TEST); else glEnable(GL_DEPTH_TEST); };
    auto drawWells = [&]() {
        if (!(set_.wells && wvao_ && wIndexCount_)) return;
        glUseProgram(progC_);
        noDepth(set_.wellsXray && !transparent);
        glBindVertexArray(wvao_);
        glDrawElements(GL_TRIANGLES, GLsizei(wIndexCount_), GL_UNSIGNED_INT, nullptr);
        noDepth(false);
    };
    auto drawLines = [&]() {
        glUseProgram(progC_);
        if (set_.bbox && lvao_ && lineCount_) { glBindVertexArray(lvao_); glDrawArrays(GL_LINES, 0, GLsizei(lineCount_)); }
        if (slvao_ && slCount_) {
            noDepth(set_.linesXray);
            glBindVertexArray(slvao_); glDrawArrays(GL_LINES, 0, GLsizei(slCount_));
            noDepth(false);
        }
    };
    auto drawPanels = [&]() {
        if (!(pvao_ && pCount_) || !set_.fractures) return;
        glUseProgram(progC_);
        noDepth(set_.fracXray);
        glBindVertexArray(pvao_); glDrawArrays(GL_TRIANGLES, 0, GLsizei(pCount_));
        if (polvao_ && polCount_) { glBindVertexArray(polvao_); glDrawArrays(GL_LINES, 0, GLsizei(polCount_)); }
        noDepth(false);
    };
    if (transparent) {
        drawWells(); drawPanels(); drawLines(); drawCells();
    } else {
        drawCells(); drawPanels();
        if (!set_.wellsXray) drawWells();
        drawLines();
        if (set_.wellsXray) drawWells();
    }
    glBindVertexArray(0);
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_1D, 0);
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_BUFFER, 0);
    glv::resolve(tg_);
    glv::endView();
    drawOverlay();
    glv::endView();
    return (ImTextureID)(intptr_t)tg_.color;
}

void Viewport3D::drawOverlay() {
    const float s = g_gfx.dpiScale;
    glv::renderOverlay(tg_, overlayScale_, 12.0f * s, 13.0f * s, 17.0f * s, [this](OverlayPainter& p) { paintOverlay(p); });
}

bool Viewport3D::saveImage(const std::string& path, int scale, std::string* err) {
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
