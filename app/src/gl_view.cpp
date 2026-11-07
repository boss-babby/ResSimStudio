// C39 OpenGL 3.3 helpers for the portable 3-D views - see gl_view.h
#include "gl_view.h"
#ifdef RESSIM_GL
#include "gl_api.h"
#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "theme.h"
#include "plutovg.h"
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace glv {

static unsigned compileStage(GLenum type, const std::string& src) {
    const GLuint s = glCreateShader(type);
    const char* p = src.c_str();
    glShaderSource(s, 1, &p, nullptr);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[4096]; GLsizei n = 0;
        glGetShaderInfoLog(s, sizeof(log), &n, log);
        std::fprintf(stderr, "[ResSim GL] %s shader error:\n%.*s\n", type == GL_VERTEX_SHADER ? "vertex" : "fragment", int(n), log);
        glDeleteShader(s);
        return 0;
    }
    return s;
}

unsigned program(const char* vs, const char* fs, const char* prelude, const char* defines) {
    const std::string head = std::string("#version 330 core\n") + defines + prelude + "\n";
    const GLuint v = compileStage(GL_VERTEX_SHADER, head + vs), f = compileStage(GL_FRAGMENT_SHADER, head + fs);
    if (!v || !f) { if (v) glDeleteShader(v); if (f) glDeleteShader(f); return 0; }
    const GLuint p = glCreateProgram();
    glAttachShader(p, v); glAttachShader(p, f);
    glLinkProgram(p);
    glDeleteShader(v); glDeleteShader(f);
    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[4096]; GLsizei n = 0;
        glGetProgramInfoLog(p, sizeof(log), &n, log);
        std::fprintf(stderr, "[ResSim GL] link error:\n%.*s\n", int(n), log);
        glDeleteProgram(p);
        return 0;
    }
    return p;
}

int loc(unsigned prog, const char* name) { return glGetUniformLocation(prog, name); }

int maxSamples(int wanted) {
    GLint m = 1;
    glGetIntegerv(GL_MAX_SAMPLES, &m);
    return std::max(1, std::min(wanted, int(m)));
}

void releaseTargets(Targets& t) {
    if (t.msFbo) glDeleteFramebuffers(1, &t.msFbo);
    if (t.msFboColorOnly) glDeleteFramebuffers(1, &t.msFboColorOnly);
    if (t.fbo) glDeleteFramebuffers(1, &t.fbo);
    if (t.msColor) glDeleteRenderbuffers(1, &t.msColor);
    if (t.msDepth) { if (t.depthTex) glDeleteTextures(1, &t.msDepth); else glDeleteRenderbuffers(1, &t.msDepth); }
    if (t.color) glDeleteTextures(1, &t.color);
    t = Targets{};
}

bool createTargets(Targets& t, int w, int h, int msaa, bool depthTex, bool colorOnlyFbo) {
    releaseTargets(t);
    t.w = w; t.h = h; t.msaa = std::max(1, msaa); t.depthTex = depthTex;
    const bool ms = t.msaa > 1;
    glGenRenderbuffers(1, &t.msColor);
    glBindRenderbuffer(GL_RENDERBUFFER, t.msColor);
    if (ms) glRenderbufferStorageMultisample(GL_RENDERBUFFER, t.msaa, GL_RGBA8, w, h);
    else glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, w, h);
    if (depthTex) {
        glGenTextures(1, &t.msDepth);
        if (ms) {
            glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, t.msDepth);
            glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, t.msaa, GL_DEPTH_COMPONENT32F, w, h, GL_TRUE);
            glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
        } else {
            glBindTexture(GL_TEXTURE_2D, t.msDepth);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32F, w, h, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
            glBindTexture(GL_TEXTURE_2D, 0);
        }
    } else {
        glGenRenderbuffers(1, &t.msDepth);
        glBindRenderbuffer(GL_RENDERBUFFER, t.msDepth);
        if (ms) glRenderbufferStorageMultisample(GL_RENDERBUFFER, t.msaa, GL_DEPTH_COMPONENT32F, w, h);
        else glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, w, h);
    }
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glGenFramebuffers(1, &t.msFbo);
    glBindFramebuffer(GL_FRAMEBUFFER, t.msFbo);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, t.msColor);
    if (depthTex) glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, ms ? GL_TEXTURE_2D_MULTISAMPLE : GL_TEXTURE_2D, t.msDepth, 0);
    else glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, t.msDepth);
    bool ok = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
    if (colorOnlyFbo) {
        glGenFramebuffers(1, &t.msFboColorOnly);
        glBindFramebuffer(GL_FRAMEBUFFER, t.msFboColorOnly);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, t.msColor);
        ok = ok && glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
    }
    glGenTextures(1, &t.color);
    glBindTexture(GL_TEXTURE_2D, t.color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);
    glGenFramebuffers(1, &t.fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, t.color, 0);
    ok = ok && glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    if (!ok) std::fprintf(stderr, "[ResSim GL] incomplete offscreen framebuffer (%dx%d, %d samples)\n", w, h, t.msaa);
    return ok;
}

void resolve(const Targets& t) {
    glBindFramebuffer(GL_READ_FRAMEBUFFER, t.msFbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, t.fbo);
    glBlitFramebuffer(0, 0, t.w, t.h, 0, 0, t.w, t.h, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void endView() {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glActiveTexture(GL_TEXTURE0);
}

// ---------------------------------------------------------------------------- ImGui draw-list painter
namespace {
inline ImU32 col32(OvColor c) {
    auto b = [](float v) { return int(std::max(0.0f, std::min(1.0f, v)) * 255.0f + 0.5f); };
    return IM_COL32(b(c.r), b(c.g), b(c.b), b(c.a));
}
class ImPainter final : public OverlayPainter {
public:
    ImPainter(ImDrawList* dl, float scale, const float sizes[3]) : dl_(dl), S_(scale) {
        for (int k = 0; k < 3; ++k) size_[k] = sizes[k];
        font_[0] = g_fonts.ui ? g_fonts.ui : ImGui::GetFont();
        font_[1] = g_fonts.bold ? g_fonts.bold : font_[0];
        font_[2] = font_[1];
    }
    void fillRoundRect(float x0, float y0, float x1, float y1, float r, OvColor c) override { dl_->AddRectFilled(P(x0, y0), P(x1, y1), col32(c), r * S_); }
    void fillRect(float x0, float y0, float x1, float y1, OvColor c) override { dl_->AddRectFilled(P(x0, y0), P(x1, y1), col32(c)); }
    void strokeRect(float x0, float y0, float x1, float y1, OvColor c, float w) override { dl_->AddRect(P(x0, y0), P(x1, y1), col32(c), 0.0f, 0, w * S_); }
    void line(float x0, float y0, float x1, float y1, OvColor c, float w) override { dl_->AddLine(P(x0, y0), P(x1, y1), col32(c), w * S_); }
    void fillCircle(float cx, float cy, float r, OvColor c) override { dl_->AddCircleFilled(P(cx, cy), r * S_, col32(c), 0); }
    float textWidth(const std::string& t, OvFont f) override {
        if (t.empty()) return 0;
        return font_[int(f)]->CalcTextSizeA(size_[int(f)], FLT_MAX, 0.0f, t.c_str(), t.c_str() + t.size()).x;
    }
    void text(const std::string& t, OvFont f, float x, float y, OvColor c) override {
        if (t.empty()) return;
        // DirectWrite puts the line box top at y with a slightly taller line than ImGui: nudge down to keep the layout
        const float sz = size_[int(f)];
        dl_->AddText(font_[int(f)], sz * S_, P(x, y + 0.08f * sz), col32(c), t.c_str(), t.c_str() + t.size());
    }

private:
    ImVec2 P(float x, float y) const { return ImVec2(x * S_, y * S_); }
    ImDrawList* dl_;
    float S_;
    float size_[3];
    ImFont* font_[3];
};
}  // namespace

void renderOverlay(const Targets& t, float scale, float fontSmall, float fontLabel, float fontTitle,
                   const std::function<void(OverlayPainter&)>& paint) {
    if (!t.fbo || !ImGui::GetCurrentContext()) return;
    ImDrawList dl(ImGui::GetDrawListSharedData());
    dl._ResetForNewFrame();
    dl.PushTexture(ImGui::GetIO().Fonts->TexRef);
    dl.PushClipRect(ImVec2(0, 0), ImVec2(float(t.w), float(t.h)), false);
    const float sizes[3] = {fontSmall, fontLabel, fontTitle};
    {
        ImPainter p(&dl, scale, sizes);
        paint(p);
    }
    if (dl.VtxBuffer.Size == 0) return;
    // flip to the top-row-first layout of the view texture
    const float H = float(t.h);
    for (ImDrawVert& v : dl.VtxBuffer) v.pos.y = H - v.pos.y;
    for (ImDrawCmd& c : dl.CmdBuffer) { const float y0 = c.ClipRect.y, y1 = c.ClipRect.w; c.ClipRect.y = H - y1; c.ClipRect.w = H - y0; }
    // glyphs baked for new sizes during this frame must reach the GPU atlas before the overlay is drawn
    // (a texture created by an atlas resize in this frame is only in the atlas list until ImGui::Render registers it)
    auto upload = [](ImTextureData* tex) {
        if (tex && (tex->Status == ImTextureStatus_WantCreate || tex->Status == ImTextureStatus_WantUpdates)) ImGui_ImplOpenGL3_UpdateTexture(tex);
    };
    for (ImTextureData* tex : ImGui::GetIO().Fonts->TexList) upload(tex);
    for (ImTextureData* tex : ImGui::GetPlatformIO().Textures) upload(tex);
    ImDrawData dd;
    dd.Valid = true;
    dd.DisplayPos = ImVec2(0, 0);
    dd.DisplaySize = ImVec2(float(t.w), float(t.h));
    dd.FramebufferScale = ImVec2(1, 1);
    dd.Textures = nullptr;
    dd.AddDrawList(&dl);
    glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    ImGui_ImplOpenGL3_RenderDrawData(&dd);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

// ---------------------------------------------------------------------------- PNG
bool writePngRGBA(const std::string& path, int w, int h, const unsigned char* rgba) {
    if (w <= 0 || h <= 0 || !rgba) return false;
    // plutovg surfaces are premultiplied ARGB32 (native endian: B, G, R, A bytes on little-endian); alpha is 255 here
    std::vector<unsigned char> argb(size_t(w) * size_t(h) * 4);
    for (size_t i = 0, n = size_t(w) * size_t(h); i < n; ++i) {
        const uint32_t px = (uint32_t(255) << 24) | (uint32_t(rgba[4 * i]) << 16) | (uint32_t(rgba[4 * i + 1]) << 8) | uint32_t(rgba[4 * i + 2]);
        std::memcpy(&argb[4 * i], &px, 4);
    }
    plutovg_surface_t* s = plutovg_surface_create_for_data(argb.data(), w, h, w * 4);
    if (!s) return false;
    const bool ok = plutovg_surface_write_to_png(s, path.c_str());
    plutovg_surface_destroy(s);
    return ok;
}

bool savePng(const Targets& t, const std::string& path, std::string* err) {
    if (!t.fbo) { if (err) *err = "render target"; return false; }
    std::vector<unsigned char> px(size_t(t.w) * size_t(t.h) * 4);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, t.fbo);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, t.w, t.h, GL_RGBA, GL_UNSIGNED_BYTE, px.data());
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    if (!writePngRGBA(path, t.w, t.h, px.data())) { if (err) *err = "PNG encoding failed"; return false; }
    return true;
}

}  // namespace glv
#endif
