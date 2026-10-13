// 3D reservoir viewer: D3D11 (4x MSAA) rendering of the cell exterior faces coloured by a property,
// wells as cylinders, then a Direct2D/DirectWrite overlay (legend, labels, axis triad) on the same texture.
#pragma once
#include "gfx.h"
#include "ressim/grid.h"
#include <vector>
#include <string>
#include <array>

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

enum class Colormap { Turbo, Viridis, Jet, CoolWarm, Grey };

struct ViewSettings {
    float zScale = 5.0f;
    bool edges = true, wells = true, labels = true, legend = true, axes = true, bbox = true;
    Colormap cmap = Colormap::Turbo;
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
};

struct LineVertex { float x, y, z; float r, g, b, a; };   // model coordinates (x, y, depth)

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
    ViewSettings& settings() { return set_; }
    const ViewSettings& settings() const { return set_; }
    void markGeometryDirty() { geomDirty_ = true; }

    // Renders into the offscreen texture (size in pixels). Returns the SRV to display with ImGui::Image.
    ID3D11ShaderResourceView* render(int w, int h);

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
    std::array<float, 16> viewProj() const;
    bool project(float x, float y, float z, float& sx, float& sy) const;   // model -> pixels
    void worldOf(float x, float y, float z, float out[3]) const;

    // D3D resources
    ID3D11VertexShader* vs_ = nullptr; ID3D11PixelShader* ps_ = nullptr;
    ID3D11VertexShader* vsC_ = nullptr; ID3D11PixelShader* psC_ = nullptr;
    ID3D11InputLayout* il_ = nullptr; ID3D11InputLayout* ilC_ = nullptr;
    ID3D11Buffer* cb_ = nullptr;
    ID3D11Buffer* vb_ = nullptr; ID3D11Buffer* ib_ = nullptr; unsigned indexCount_ = 0;
    ID3D11Buffer* wvb_ = nullptr; ID3D11Buffer* wib_ = nullptr; unsigned wIndexCount_ = 0;
    ID3D11Buffer* lvb_ = nullptr; unsigned lineCount_ = 0;
    ID3D11Buffer* slvb_ = nullptr; unsigned slCount_ = 0;
    ID3D11DepthStencilState* dssOff_ = nullptr;
    ID3D11DepthStencilState* dssNoWrite_ = nullptr;
    std::vector<float> faceCenters_;          // model-space centre per face (x, y, z) for back-to-front sorting
    std::vector<unsigned> sortedIdx_;
    float sortedCam_[7] = {1e30f, 0, 0, 0, 0, 0, 0};
    int builtFocus_ = -2; float builtRadius_ = -1;
    float builtWellZ_ = -1, overlayScale_ = 1.0f;
    std::vector<LineVertex> lines_;
    bool linesDirty_ = false;
    ID3D11Buffer* valBuf_ = nullptr; ID3D11ShaderResourceView* valSRV_ = nullptr; int valCap_ = 0;
    ID3D11Texture1D* cmapTex_ = nullptr; ID3D11ShaderResourceView* cmapSRV_ = nullptr;
    ID3D11SamplerState* samp_ = nullptr;
    ID3D11RasterizerState* rs_ = nullptr;
    ID3D11DepthStencilState* dss_ = nullptr;
    ID3D11BlendState* bs_ = nullptr;
    // targets
    int w_ = 0, h_ = 0;
    ID3D11Texture2D* msColor_ = nullptr; ID3D11RenderTargetView* msRTV_ = nullptr;
    ID3D11Texture2D* msDepth_ = nullptr; ID3D11DepthStencilView* msDSV_ = nullptr;
    ID3D11Texture2D* color_ = nullptr; ID3D11ShaderResourceView* colorSRV_ = nullptr;
    ID2D1Bitmap1* d2dTarget_ = nullptr;
    IDWriteTextFormat* tfSmall_ = nullptr; IDWriteTextFormat* tfLabel_ = nullptr; IDWriteTextFormat* tfTitle_ = nullptr;
    UINT msaa_ = 4;

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
    bool geomDirty_ = true, wellsDirty_ = true, valuesDirty_ = true;
    float center_[3] = {0, 0, 0}, extent_ = 1;
    // camera
    float yaw_ = 0.7f, pitch_ = 0.55f, dist_ = 2.5f, panX_ = 0, panY_ = 0, panZ_ = 0;
};

// Colormap evaluation shared with 2D plots/legends (t in [0,1]) -> RGB 0..1
void colormapRGB(Colormap m, float t, float rgb[3]);
