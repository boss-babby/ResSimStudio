// C39 OpenGL 3.3 helpers shared by the portable 3-D views (viewport3d_gl.cpp, geoscene3d_gl.cpp):
// shader programs, MSAA offscreen targets, the ImGui draw-list overlay painter and PNG output.
//
// Orientation: the views render with the clip-space y flipped (and the overlay vertices flipped), so the resolved texture
// stores the TOP image row first - the same layout as the Direct3D path. ImGui::Image shows it with the default UVs and
// glReadPixels returns rows top-down for the PNG snapshots.
#pragma once
#include "gfx.h"
#ifdef RESSIM_GL
#include "overlay_paint.h"
#include <functional>
#include <string>
#include <vector>

namespace glv {

// compiles + links a program from GLSL 330 sources (a common prelude is prepended to both stages); 0 on failure
// (the log goes to stderr). 'defines' is inserted after the #version line, e.g. "#define MSDEPTH 1\n".
unsigned program(const char* vs, const char* fs, const char* prelude = "", const char* defines = "");
int loc(unsigned prog, const char* name);

// offscreen targets: multisampled colour (+ depth: renderbuffer, or a multisampled depth TEXTURE when depthTex) and a
// single-sample RGBA8 texture that receives the resolve. msFboColorOnly (optional) = same colour attachment without depth.
struct Targets {
    unsigned msFbo = 0, msColor = 0, msDepth = 0, fbo = 0, color = 0, msFboColorOnly = 0;
    int w = 0, h = 0, msaa = 1;
    bool depthTex = false;
};
bool createTargets(Targets& t, int w, int h, int msaa, bool depthTex, bool colorOnlyFbo);
void releaseTargets(Targets& t);
void resolve(const Targets& t);                     // MSAA colour -> single-sample texture
// supported sample count (<= wanted)
int maxSamples(int wanted);

// paints an overlay (logical units * scale) into the resolved texture of 't'
void renderOverlay(const Targets& t, float scale, float fontSmall, float fontLabel, float fontTitle,
                   const std::function<void(OverlayPainter&)>& paint);

// reads the resolved texture (top row first) and writes a PNG (alpha forced to 255)
bool savePng(const Targets& t, const std::string& path, std::string* err);
// RGBA8 rows (top first) -> PNG
bool writePngRGBA(const std::string& path, int w, int h, const unsigned char* rgba);

// restores the GL state ImGui's renderer expects after a view was drawn mid-frame
void endView();

}  // namespace glv
#endif
