// C36 3-D geomodel / seismic scene (Direct3D 11 + Direct2D / DirectWrite overlay), used by the Geomodel workspace and
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

struct ID3D11ShaderResourceView;
// recording mock of the D3D11 renderer
class GeoScene3D {
public:
    bool init() { return true; }
    void shutdown() {}
    bool ok() const { return true; }
    void clear() {}
    std::vector<GeoSliceDesc> slices; std::vector<GeoMeshDesc> meshes; std::vector<GeoSolidDesc> solids; std::vector<GeoLineVtx> lines;
    std::vector<GeoLabel> labels; bool hasVol = false; GeoVolumeDesc vol; float blo[3] = {0, 0, 0}, bhi[3] = {1, 1, 1};
    int nSlices = 0, nMeshes = 0, nSolids = 0, nVol = 0;
    void setSlices(const std::vector<GeoSliceDesc>& s) { slices = s; ++nSlices; }
    void setMeshes(const std::vector<GeoMeshDesc>& m) { meshes = m; ++nMeshes; }
    void setSolids(const std::vector<GeoSolidDesc>& s) { solids = s; ++nSolids; }
    void setLines(const std::vector<GeoLineVtx>& l) { lines = l; }
    void setVolume(const GeoVolumeDesc* v) { hasVol = v != nullptr; if (v) vol = *v; ++nVol; }
    void setVolumeStyle(const GeoVolumeDesc&) {}
    void setLabels(const std::vector<GeoLabel>& l) { labels = l; }
    void setBounds(const float lo[3], const float hi[3]) { for (int m = 0; m < 3; ++m) { blo[m] = lo[m]; bhi[m] = hi[m]; } }
    void setTitle(const std::string&, const std::string&, float, float, Colormap, bool) {}
    void setCaption(const std::string&) {}
    float zScale = 3.0f;
    bool showLegend = true, showAxes = true, showLabels = true, showBox = true;
    float background[3] = {0, 0, 0};
    ID3D11ShaderResourceView* render(int, int) { return nullptr; }
    void orbit(float, float) {} void pan(float, float) {} void zoom(float) {} void resetCamera() {} void viewPreset(int) {}
    bool pixelRay(float, float, double*, double*) const { return false; }
    bool project(float, float, float, float&, float&) const { return false; }
    bool saveImage(const std::string&, int, std::string* = nullptr) { return false; }
    int width() const { return 0; }
    int height() const { return 0; }
};
