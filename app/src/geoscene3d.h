// C36 3-D geomodel / seismic scene (Windows: Direct3D 11 + Direct2D / DirectWrite overlay; C39 Linux / macOS: OpenGL 3.3 +
// ImGui draw-list overlay in geoscene3d_gl.cpp), used by the Geomodel workspace and
// ResGeo. Independent of the simulation grid viewer (viewport3d): the scene holds a list of objects in model
// coordinates (x east, y north, z = display depth, positive down - true depth or two-way time x a scale):
//   * slices   - textured quads / curved sheets (seismic inline / crossline / time slices, attribute sections, maps);
//                a float texture coloured on the GPU with a colour map, clip range, opacity, optional value-alpha
//                (e.g. coherence over amplitude), NaN = transparent
//   * meshes   - triangle meshes with a value per vertex (horizons coloured by depth / property with GPU contour lines,
//                grid shells and sections with flat cell colours)
//   * solids   - lit triangle meshes with a colour per vertex (fault surfaces, well tubes, top markers)
//   * lines    - line lists (fault sticks, picks, survey outline, cell edges)
//   * volume   - GPU ray-marched seismic volume (3-D float texture of a proxy cube, amplitude-dependent opacity
//                transfer function, depth-composited with the opaque scene, time or depth domain through the
//                v0 + k z velocity model)
// Labels (wells, tops, slice annotations), title, legend and an axis triad are drawn with Direct2D on top and go into
// the PNG snapshots.
#pragma once
#include "gfx.h"
#include "imgui.h"
class OverlayPainter;
#ifndef RESSIM_D3D11
#include "gl_view.h"
#endif
#include "colormaps.h"
#include <array>
#include <string>
#include <vector>

struct GeoLabel { float p[3]; std::string text; float rgba[4] = {0.15f, 0.45f, 0.75f, 0.92f}; bool small = false; };

struct GeoSliceDesc {
    int nu = 2, nv = 2;                  // vertex lattice (nu along the texture width, nv along the height)
    std::vector<float> pos;              // nu*nv*3 model coordinates, u fastest
    int tw = 0, th = 0;                  // texture size
    std::vector<float> tex;              // tw*th values, row-major (u fastest); NaN = transparent
    float lo = -1, hi = 1;               // colour range
    Colormap cmap = Colormap::Seismic; bool reverse = false;
    float opacity = 1.0f;
    int alphaMode = 0;                   // 0 constant, 1 alpha grows with |t - 0.5| (amplitudes), 2 alpha grows as t falls (coherence: edges opaque)
};
struct GeoMeshDesc {                     // value mesh (vertex values coloured on the GPU)
    std::vector<float> pos;              // n*3
    std::vector<float> nrm;              // n*3 (empty: computed per triangle, flat)
    std::vector<float> val;              // n (NaN = grey)
    std::vector<unsigned> idx;           // triangles
    float lo = 0, hi = 1;
    Colormap cmap = Colormap::Turbo; bool reverse = false;
    float opacity = 1.0f;
    float contour = 0;                   // contour interval in value units (0 = none)
    bool edges = false;                  // quad-cell edge lines from 'uv' (per vertex 0/1 corner coordinates)
    std::vector<float> uv;               // n*2 (optional, for edges)
};
struct GeoSolidDesc {                    // lit mesh with a colour per vertex
    std::vector<float> pos, nrm, col;    // n*3, n*3 (empty: flat), n*4
    std::vector<unsigned> idx;
    bool transparent = false;
};
struct GeoLineVtx { float x, y, z, r, g, b, a; };
struct GeoVolumeDesc {
    int ni = 0, nj = 0, ns = 0;
    std::vector<float> data;             // [i][j][s] as SeismicCube
    double x0 = 0, y0 = 0, ix = 1, iy = 0, jx = 0, jy = 1;   // plan geometry of the lattice
    double t0 = 0, dt = 4;               // sample axis (ms or depth units)
    int zMode = 0;                       // display z of sample s: 0 z = zScale * (t0 + s dt)   1 depth from TWT by the velocity model
    double zScale = 1.0;                 // mode 0
    double v0 = 2000, k = 0, datum = 0;  // mode 1 (geo::VelocityModel)
    float lo = -1, hi = 1;               // amplitude range (colour + transfer function)
    Colormap cmap = Colormap::Seismic; bool reverse = false;
    float opacity = 0.6f;                // overall opacity scale
    float threshold = 0.25f;             // |t - 0.5| * 2 below this is fully transparent
    int steps = 384;                     // ray-march samples across the volume
    int iMin = 0, iMax = 1 << 30, jMin = 0, jMax = 1 << 30, sMin = 0, sMax = 1 << 30;   // sub-volume (indices)
};

class GeoScene3D {
public:
    bool init();
    void shutdown();
    bool ok() const { return ready_; }

    // scene content (replace whole categories; ids are positions in the lists)
    void clear();
    void setSlices(const std::vector<GeoSliceDesc>& s);
    void setMeshes(const std::vector<GeoMeshDesc>& m);
    void setSolids(const std::vector<GeoSolidDesc>& s);
    void setLines(const std::vector<GeoLineVtx>& l);
    void setVolume(const GeoVolumeDesc* v);          // nullptr = none
    void setVolumeStyle(const GeoVolumeDesc& v);     // colour / opacity / threshold / sub-volume without re-uploading
    void setLabels(const std::vector<GeoLabel>& l) { labels_ = l; }
    void setBounds(const float lo[3], const float hi[3]);   // model box (frames the camera)
    void setTitle(const std::string& t, const std::string& legendTitle, float lo, float hi, Colormap cm, bool rev) {
        title_ = t; legTitle_ = legendTitle; legLo_ = lo; legHi_ = hi; legCmap_ = cm; legRev_ = rev;
    }
    void setCaption(const std::string& s) { caption_ = s; }

    float zScale = 3.0f;                 // vertical exaggeration
    bool showLegend = true, showAxes = true, showLabels = true, showBox = true;
    float background[3] = {0.075f, 0.08f, 0.095f};

    ImTextureID render(int w, int h);   // texture for ImGui::Image (0 = none)
    void orbit(float dx, float dy);
    void pan(float dx, float dy);
    void zoom(float wheel);
    void resetCamera();
    void viewPreset(int which);          // 0 top, 1 south, 2 east, 3 iso
    // ray through a pixel in model coordinates (origin, unit direction); false before the first render
    bool pixelRay(float px, float py, double o[3], double d[3]) const;
    bool project(float x, float y, float z, float& sx, float& sy) const;
    bool saveImage(const std::string& path, int scale, std::string* err = nullptr);
    int width() const { return w_; }
    int height() const { return h_; }

private:
    void releaseObjects();
    void releaseVolume();
    void createTargets(int w, int h);
    void releaseTargets();
    std::array<float, 16> viewProj() const;
    void eyePos(float e[3]) const;
    void worldOf(float x, float y, float z, float o[3]) const;
    void drawOverlay();
    void paintOverlay(OverlayPainter& p);
    bool ready_ = false;
    int w_ = 0, h_ = 0;
    float overlayScale_ = 1.0f;
    std::vector<GeoLabel> labels_;
    unsigned lineCount_ = 0, boxCount_ = 0;
    GeoVolumeDesc vol_; bool hasVol_ = false; Colormap volCmapBuilt_ = Colormap(-1); bool volRevBuilt_ = false;
#ifdef RESSIM_D3D11
    struct Slice { ID3D11Buffer* vb = nullptr; ID3D11Buffer* ib = nullptr; unsigned n = 0; ID3D11Texture2D* tex = nullptr; ID3D11ShaderResourceView* srv = nullptr;
                   ID3D11Texture1D* cm = nullptr; ID3D11ShaderResourceView* cmSrv = nullptr; float lo = 0, hi = 1, opacity = 1; int alphaMode = 0; };
    struct Mesh { ID3D11Buffer* vb = nullptr; ID3D11Buffer* ib = nullptr; unsigned n = 0; ID3D11Texture1D* cm = nullptr; ID3D11ShaderResourceView* cmSrv = nullptr;
                  float lo = 0, hi = 1, opacity = 1, contour = 0; bool edges = false; float cz[3] = {0, 0, 0}; };
    struct Solid { ID3D11Buffer* vb = nullptr; ID3D11Buffer* ib = nullptr; unsigned n = 0; bool transparent = false; };
    ID3D11ShaderResourceView* makeColormap(Colormap cm, bool rev, ID3D11Texture1D** tex);

    // pipeline
    ID3D11VertexShader* vs_ = nullptr; ID3D11PixelShader* psMesh_ = nullptr; ID3D11InputLayout* il_ = nullptr;      // value meshes
    ID3D11VertexShader* vsT_ = nullptr; ID3D11PixelShader* psT_ = nullptr; ID3D11InputLayout* ilT_ = nullptr;       // slices
    ID3D11VertexShader* vsC_ = nullptr; ID3D11PixelShader* psC_ = nullptr; ID3D11InputLayout* ilC_ = nullptr;       // solids / lines
    ID3D11VertexShader* vsQ_ = nullptr; ID3D11PixelShader* psV_ = nullptr;                                          // volume (full-screen)
    ID3D11Buffer* cbFrame_ = nullptr; ID3D11Buffer* cbObj_ = nullptr;
    ID3D11SamplerState* samp_ = nullptr; ID3D11SamplerState* sampPt_ = nullptr;
    ID3D11RasterizerState* rs_ = nullptr;
    ID3D11DepthStencilState* dss_ = nullptr; ID3D11DepthStencilState* dssNoWrite_ = nullptr;
    ID3D11BlendState* bs_ = nullptr; ID3D11BlendState* bsPre_ = nullptr;   // straight alpha / premultiplied (volume)
    IDWriteTextFormat* tfSmall_ = nullptr; IDWriteTextFormat* tfLabel_ = nullptr; IDWriteTextFormat* tfTitle_ = nullptr;
    UINT msaa_ = 4;
    // targets
    ID3D11Texture2D* msColor_ = nullptr; ID3D11RenderTargetView* msRTV_ = nullptr;
    ID3D11Texture2D* msDepth_ = nullptr; ID3D11DepthStencilView* msDSV_ = nullptr; ID3D11ShaderResourceView* depthSRV_ = nullptr;
    ID3D11Texture2D* color_ = nullptr; ID3D11ShaderResourceView* colorSRV_ = nullptr;
    ID2D1Bitmap1* d2dTarget_ = nullptr;
    // objects
    ID3D11Buffer* lvb_ = nullptr;
    ID3D11Buffer* boxVb_ = nullptr;
    // volume
    ID3D11Texture3D* volTex_ = nullptr; ID3D11ShaderResourceView* volSrv_ = nullptr;
    ID3D11Texture1D* volCm_ = nullptr; ID3D11ShaderResourceView* volCmSrv_ = nullptr;
#else
    // OpenGL 3.3 (names are GLuint)
    struct Slice { unsigned vao = 0, vb = 0, ib = 0, n = 0, tex = 0, cm = 0; float lo = 0, hi = 1, opacity = 1; int alphaMode = 0; };
    struct Mesh { unsigned vao = 0, vb = 0, ib = 0, n = 0, cm = 0; float lo = 0, hi = 1, opacity = 1, contour = 0; bool edges = false; float cz[3] = {0, 0, 0}; };
    struct Solid { unsigned vao = 0, vb = 0, ib = 0, n = 0; bool transparent = false; };
    unsigned progMesh_ = 0, progTex_ = 0, progC_ = 0, progVol_ = 0, emptyVao_ = 0;
    unsigned lvao_ = 0, lvb_ = 0, boxVao_ = 0, boxVb_ = 0;
    unsigned volTex_ = 0, volCm_ = 0;
    glv::Targets tg_;
    int msaa_ = 4;
#endif
    std::vector<Slice> slices_;
    std::vector<Mesh> meshes_;
    std::vector<Solid> solids_;
    // view
    float lo_[3] = {0, 0, 0}, hi_[3] = {1, 1, 1}, center_[3] = {0.5f, 0.5f, 0.5f}, extent_ = 1;
    float builtBoxZ_ = -1;
    float yaw_ = 0.7f, pitch_ = 0.5f, dist_ = 2.0f, panX_ = 0, panY_ = 0, panZ_ = 0;
    std::string title_, legTitle_, caption_;
    float legLo_ = 0, legHi_ = 1; Colormap legCmap_ = Colormap::Seismic; bool legRev_ = false;
};
