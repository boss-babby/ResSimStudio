// 3D reservoir viewer: 4x MSAA rendering of the cell exterior faces coloured by a property, wells as cylinders, then a
// 2-D overlay (legend, labels, axis triad) on the same texture.
//   Windows: Direct3D 11 + Direct2D / DirectWrite (viewport3d.cpp)   Linux / macOS: OpenGL 3.3 + ImGui draw list (viewport3d_gl.cpp)
// Geometry building, camera, picking and the overlay layout are shared (viewport3d.cpp).
#pragma once
#include "gfx.h"
#include "imgui.h"
class OverlayPainter;
#ifndef RESSIM_D3D11
#include "gl_view.h"
#endif
#include "ressim/grid.h"
#include "colormaps.h"
#include <vector>
#include <string>
#include <array>
#include <cmath>

struct WellGlyph {
    std::string name;
    // trajectory in model coordinates (x, y, depth): wellhead first, then the connection cells in completion order
    // (deviated / horizontal wells follow their real path)
    std::vector<std::array<float, 3>> path;
    // perforated segments, each through one connection cell along its penetration direction (X, Y or Z)
    struct Perf { std::array<float, 3> a{}, b{}; bool open = true; };
    std::vector<Perf> perfs;
    bool injector = false;
    bool open = true;
};


struct ViewSettings {
    float zScale = 5.0f;
    bool edges = true, wells = true, labels = true, legend = true, axes = true, bbox = true;
    Colormap cmap = Colormap::Turbo;
    bool cmapReverse = false;     // C30: reversed colour map
    int iMin = 0, iMax = 1 << 30, jMin = 0, jMax = 1 << 30, kMin = 0, kMax = 1 << 30;
    bool filterValue = false;
    float filterMin = 0, filterMax = 1;
    bool logScale = false;
    bool showCells = true;        // draw the cell faces (off: only lines / wells, e.g. for streamlines)
    bool linesXray = false;       // draw overlay lines without depth test
    float opacity = 1.0f;         // cell opacity (< 1: see-through grid, wells and streamlines visible inside)
    float wellRadius = 1.0f;      // well tube radius multiplier
    int focusWell = -1;           // index into the wells: highlighted, others dimmed (-1 = none)
    int bands = 0;                // discrete colour bands (0 = continuous)
    bool wellsXray = false;       // wells drawn on top of the cells
    bool fractures = true;        // C13: draw the overlay panels (EDFM fracture planes)
    bool fracXray = true;         // panels drawn on top of the cells (they lie inside the grid)
};

struct LineVertex { float x, y, z; float r, g, b, a; };   // model coordinates (x, y, depth)
// C13: flat quads drawn inside the grid, e.g. EDFM fracture planes (model coordinates, corners in order); value = NaN
// -> fixed colour 'rgba', else coloured with the current colour map and range of the cell property ('byValue')
struct PanelQuad { std::array<float, 3> p[4]; float rgba[4] = {1.0f, 0.55f, 0.15f, 0.85f}; float value = NAN; std::string label; };

class Viewport3D {
public:
    bool init();
    void shutdown();

    // Geometry source. Call again whenever the grid changes (pointer must stay valid until next call).
    void setGrid(const ressim::Grid* grid);
    // Per active cell values; NaN = undefined (drawn grey).
    void setValues(const std::vector<float>& v, float vmin, float vmax, const std::string& title, const std::string& units);
    void setWells(const std::vector<WellGlyph>& w) { wells_ = w; wellsDirty_ = true; }
    void setCaption(const std::string& s) { caption_ = s; }
    // Overlay polylines as a line list (pairs of vertices), e.g. streamlines. Empty = none.
    void setLines(std::vector<LineVertex> v) { lines_ = std::move(v); linesDirty_ = true; }
    // C13: overlay panels (EDFM fractures): solid quads + outline, drawn with the cells; labels at the first panel of a name
    void setPanels(std::vector<PanelQuad> q) { panels_ = std::move(q); panelsDirty_ = true; }
    bool hasPanels() const { return !panels_.empty(); }
    ViewSettings& settings() { return set_; }
    const ViewSettings& settings() const { return set_; }
    void markGeometryDirty() { geomDirty_ = true; }

    // Renders into the offscreen texture (size in pixels). Returns the texture to display with ImGui::Image (0 = none).
    ImTextureID render(int w, int h);

    // Interaction (pixels)
    void orbit(float dx, float dy);
    void pan(float dx, float dy);
    void zoom(float wheel);
    void resetCamera();
    void viewPreset(int which);   // 0 top, 1 front (south), 2 side (east), 3 iso
    // Returns active cell index under the pixel or -1.
    int pick(float px, float py) const;
    // Index of the well whose trajectory passes within 'tol' pixels of (px, py), or -1.
    int pickWell(float px, float py, float tol = 8.0f) const;
    // Centres the camera on well 'idx' (keeps the view direction) and frames its trajectory.
    void focusWellCamera(int idx);
    const std::vector<WellGlyph>& wells() const { return wells_; }
    // Saves the current view (cells, wells, overlay) as PNG; scale renders at a multiple of the view size.
    bool saveImage(const std::string& path, int scale, std::string* err = nullptr);
    int width() const { return w_; }
    int height() const { return h_; }
    float valueOf(int active) const { return (active >= 0 && active < int(values_.size())) ? values_[active] : 0.0f; }

private:
    struct Vtx { float p[3]; float n[3]; float uv[2]; unsigned cell; };
    struct VtxC { float p[3]; float n[3]; float c[4]; };
    void createTargets(int w, int h);
    void releaseTargets();
    void buildGeometry();
    void buildWellGeometry();
    void sortFaces(bool force);
    void uploadValues();
    void updateColormap();
    void drawOverlay();
    void paintOverlay(OverlayPainter& p);
    // back-end uploads (D3D11 in viewport3d.cpp, OpenGL in viewport3d_gl.cpp)
    void gpuReleaseCells();
    void gpuCells(const std::vector<Vtx>& vtx, const std::vector<unsigned>& idx, const std::vector<VtxC>& box);
    void gpuReleaseWells();
    void gpuWells(const std::vector<VtxC>& v, const std::vector<unsigned>& idx);
    void gpuReleasePanels();
    void gpuPanels(const std::vector<VtxC>& tri, const std::vector<VtxC>& lines);
    void gpuValues(const std::vector<float>& v);
    bool gpuHasColormap() const;
    void gpuColormap(const unsigned char* rgba256);
    bool gpuHasIndexBuffer() const;
    void gpuIndexUpdate();
    void gpuStreamlines(const std::vector<VtxC>& lv);
    void prepareFrame();          // rebuilds dirty geometry / values / colour map / lines / panels before drawing
    std::array<float, 16> viewProj() const;
    bool project(float x, float y, float z, float& sx, float& sy) const;   // model -> pixels
    void worldOf(float x, float y, float z, float out[3]) const;

    // GPU object counts (all back ends)
    unsigned indexCount_ = 0, wIndexCount_ = 0, lineCount_ = 0, slCount_ = 0, pCount_ = 0, polCount_ = 0;
    int valCap_ = 0;
    std::vector<float> faceCenters_;          // model-space centre per face (x, y, z) for back-to-front sorting
    std::vector<unsigned> sortedIdx_;
    float sortedCam_[7] = {1e30f, 0, 0, 0, 0, 0, 0};
    int builtFocus_ = -2; float builtRadius_ = -1;
    float builtWellZ_ = -1, overlayScale_ = 1.0f;
    std::vector<LineVertex> lines_;
    bool linesDirty_ = false;
    std::vector<PanelQuad> panels_;
    bool panelsDirty_ = false;
    void buildPanels();
    float builtPanelsZ_ = -1, builtPanelsLo_ = 0, builtPanelsHi_ = 0; Colormap builtPanelsCmap_ = Colormap(-1); bool builtPanelsLog_ = false;
    int w_ = 0, h_ = 0;
#ifdef RESSIM_D3D11
    // D3D resources
    ID3D11VertexShader* vs_ = nullptr; ID3D11PixelShader* ps_ = nullptr;
    ID3D11VertexShader* vsC_ = nullptr; ID3D11PixelShader* psC_ = nullptr;
    ID3D11InputLayout* il_ = nullptr; ID3D11InputLayout* ilC_ = nullptr;
    ID3D11Buffer* cb_ = nullptr;
    ID3D11Buffer* vb_ = nullptr; ID3D11Buffer* ib_ = nullptr;
    ID3D11Buffer* wvb_ = nullptr; ID3D11Buffer* wib_ = nullptr;
    ID3D11Buffer* lvb_ = nullptr;
    ID3D11Buffer* slvb_ = nullptr;
    ID3D11DepthStencilState* dssOff_ = nullptr;
    ID3D11DepthStencilState* dssNoWrite_ = nullptr;
    ID3D11Buffer* pvb_ = nullptr;      // panel triangles
    ID3D11Buffer* polvb_ = nullptr;    // panel outlines
    ID3D11Buffer* valBuf_ = nullptr; ID3D11ShaderResourceView* valSRV_ = nullptr;
    ID3D11Texture1D* cmapTex_ = nullptr; ID3D11ShaderResourceView* cmapSRV_ = nullptr;
    ID3D11SamplerState* samp_ = nullptr;
    ID3D11RasterizerState* rs_ = nullptr;
    ID3D11DepthStencilState* dss_ = nullptr;
    ID3D11BlendState* bs_ = nullptr;
    // targets
    ID3D11Texture2D* msColor_ = nullptr; ID3D11RenderTargetView* msRTV_ = nullptr;
    ID3D11Texture2D* msDepth_ = nullptr; ID3D11DepthStencilView* msDSV_ = nullptr;
    ID3D11Texture2D* color_ = nullptr; ID3D11ShaderResourceView* colorSRV_ = nullptr;
    ID2D1Bitmap1* d2dTarget_ = nullptr;
    IDWriteTextFormat* tfSmall_ = nullptr; IDWriteTextFormat* tfLabel_ = nullptr; IDWriteTextFormat* tfTitle_ = nullptr;
    UINT msaa_ = 4;
#else
    // OpenGL 3.3 resources (names are GLuint)
    unsigned prog_ = 0, progC_ = 0;                       // cell shader, coloured shader
    unsigned vao_ = 0, vb_ = 0, ib_ = 0;                  // cells
    unsigned wvao_ = 0, wvb_ = 0, wib_ = 0;               // wells
    unsigned lvao_ = 0, lvb_ = 0;                         // bounding box
    unsigned slvao_ = 0, slvb_ = 0;                       // streamlines
    unsigned pvao_ = 0, pvb_ = 0, polvao_ = 0, polvb_ = 0;  // panels
    unsigned valBuf_ = 0, valTex_ = 0;                    // per-cell values (texture buffer)
    unsigned cmapTex_ = 0;
    glv::Targets tg_;
    int msaa_ = 4;
#endif
    // data
    const ressim::Grid* grid_ = nullptr;
    std::vector<float> values_;
    float vmin_ = 0, vmax_ = 1;
    std::string title_, units_, caption_;
    std::vector<WellGlyph> wells_;
    ViewSettings set_;
    ViewSettings built_;          // settings used for the current geometry
    float builtZ_ = -1;
    Colormap builtCmap_ = Colormap(-1);
    bool builtRev_ = false;
    bool geomDirty_ = true, wellsDirty_ = true, valuesDirty_ = true;
    float center_[3] = {0, 0, 0}, extent_ = 1;
    // camera
    float yaw_ = 0.7f, pitch_ = 0.55f, dist_ = 2.5f, panX_ = 0, panY_ = 0, panZ_ = 0;
};

// Colormap evaluation shared with 2D plots/legends (t in [0,1]) -> RGB 0..1
