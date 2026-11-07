// C36: Geomodel workspace state shared by ui_geomodel.cpp (workflow, 2-D views) and ui_geo3d.cpp (3-D scene).
#pragma once
#include <deque>
#include <functional>
#include "app.h"
#include "ressim/geomodel.h"
#include "ressim/seismic3d.h"
#include "ressim/geo_import.h"
#include "ressim/geo_tools.h"
#include "ressim/geo_project.h"
#include <array>
#include <chrono>
#include <functional>
#include <future>
#include <memory>
#include <string>
#include <vector>

struct Geo3DView;   // ui_geo3d.cpp

// C42 import wizard (column assignment) state
struct GeoImportWiz {
    bool open = false, requestOpen = false;
    ressim::geo::ImportKind kind = ressim::geo::ImportKind::WellHeads;
    std::vector<std::string> files;      // empty + clip = clipboard / pasted text
    std::string text;                    // pasted text source
    bool fromText = false;
    int previewFile = 0;
    ressim::geo::TableFormat fmt;
    ressim::geo::TextTable preview;
    ressim::geo::ColumnMap map;
    ressim::geo::ImportOptions opt;
    std::string err, nullsText = "-999.25, -999, -9999, 9999900, NaN", targetWell;   // targetWell: import into this well only
    int wellFrom = 0;                    // deviation / tops / LAS without a well column: 0 file name, 1 Petrel header, 2 fixed
    char fixedWell[64] = "";
    bool createMissing = true, gridPoints = true, dirty = true, mapTouched = false;
    std::vector<int> unitSel;            // per column: 0 as is, 1 ft->m, 2 m->ft, 3 custom
};
// C42 well manager (manual heads, per-well trajectory and tops editors)
struct GeoWellEdit {
    bool open = false;
    int sel = -1;
    int mode = 0;                        // trajectory table: 0 MD / INC / AZI, 1 MD / X / Y / TVDSS
    std::string err;
    bool live = true;
};

using namespace ressim;
using namespace ressim::geo;

struct App::GeoState {
    bool field = false;
    std::string dir;
    std::vector<Surface> hz;
    std::vector<Fault> faults;
    // C41 fault network: detected junctions / truncation rules, tolerance, faults before truncation (revert), throw QC
    std::vector<FaultJunction> fnJ;
    double fnTol = 0;                             // 0 = 2 grid cells
    std::vector<Fault> fnBackup;
    std::string fnMsg;
    int fnThrowF = 0, fnThrowH = 0;
    ThrowProfile fnThrow; bool fnThrowOk = false;
    std::vector<GeoWell> wells;
    // C42 import / manual wells / 2-D seismic
    GeoImportWiz wiz;
    GeoWellEdit we;
    SegyByteMap segyMap;
    std::vector<SeismicLine> lines;               // 2-D lines
    std::vector<LineHorizon> lineHz;              // horizons interpreted on the 2-D lines
    std::vector<char> lineOn;                     // per line: shown in the 3-D scene / map
    int curLine = -1, curLineHz = 0, linePick = 1;   // section source (-1 = 3-D volume), pick mode 0 off 1 track 2 manual 3 erase 4 fault sticks
    std::vector<std::pair<int, double>> manualPts;   // manual picks of the segment being drawn (trace, TWT)
    int lineRev = 0;                              // bumped on every line / pick edit (3-D cache key)
    double lineClip = 0.99;
    SeismicCube cube;
    bool hasCube = false;
    std::vector<std::string> log;
    int step = 0, view = 0;
    // seismic
    VelocityModel vm;
    int il = 0, trackMode = 0;
    double window = 12, maxJump = 12;
    std::vector<std::array<double, 3>> seeds;
    TimeHorizon th;
    bool hasTh = false;
    std::string thName = "TopRes";
    // grid
    GridSpec gs;
    GeoGrid g;
    bool hasGrid = false;
    std::vector<Blocked> blocked;
    FaciesSpec fs;
    FaciesResult fr;
    bool hasFacies = false;
    PetroSpec ps;
    PetroResult phi;
    std::vector<PoroPerm> fits;
    std::vector<double> kx, kz, ntg;
    double kvkh = 0.1, scatter = 0.2;
    bool useSeis = false, ntgFromFacies = true;
    double corr = 0.6;
    bool hasPetro = false;
    VolSpec vs;
    VolResult vol;
    Percentiles pct;
    int nReal = 5;
    bool hasVol = false;
    UpscaleSpec us;
    Upscaled up;
    bool hasUp = false;
    // display
    int mapSurface = 0, gridProp = 2, gridLayer = 0, sectionRow = 0, corrCurve = 0;
    // worker
    std::future<std::function<void(GeoState&)>> job;
    std::string jobName;
    std::deque<std::pair<std::string, std::function<std::function<void(GeoState&)>()>>> jobQueue;   // C43: jobs started while busy
    std::chrono::steady_clock::time_point jobStart;
    // C36: out-of-core seismic, attributes, interpretation, 3-D scene
    std::shared_ptr<SegyVolume> segy;             // large SEG-Y opened out-of-core (full resolution); 'cube' is its proxy
    int xl = 0;                                   // crossline index (full resolution of the active source)
    double tSlice = -1;                           // time slice TWT (ms; <0: middle)
    int secOrient = 0;                            // seismic section view: 0 inline, 1 crossline
    int secAttr = 0, attrWin = 6;                 // SeisAttr shown on sections / slices, attribute window (samples)
    double clipPct = 0.99;                        // amplitude clip percentile
    bool flatten = false;                         // flatten the section on the tracked horizon
    int pickMode = 0;                             // 0 navigate, 1 tracking seeds, 2 fault stick digitising, 3 horizon brush
    std::vector<P3> stick;                        // fault stick being digitised (x, y, depth)
    Fault draft;                                  // fault being interpreted (finished sticks)
    std::shared_ptr<Geo3DView> v3;                // 3-D scene (ui_geo3d.cpp)
    // C36-4: horizon editing (pickMode 3 = brush on the map / 3-D scene), MPS facies, uncertainty loop
    int hzRev = 0;                                // bumped on every horizon edit (3-D cache key)
    // C42 property calculator results (one value per cell, NaN = undefined) and a revision counter for edits in place
    std::vector<std::pair<std::string, std::vector<double>>> userProps;
    int propRev = 0;
    // C43 polygons (pickMode 4 = digitising a polygon on the map, 5 = an intersection path), volumetrics limits,
    // layer-cake velocity model, mis-ties
    std::vector<ressim::geo::Polygon> polys;   // qualified: windows.h declares a GDI Polygon()
    ressim::geo::Polygon polyDraft;
    int curPoly = -1, volPoly = -1, isecPath = -1;   // isecPath: polygon index of the intersection path (Intersection view)
    bool volPerBlock = false;
    std::vector<Contacts> blockContacts;
    struct LcRow { int src = 0; std::string top; int mode = 0; double v = 2500; };   // src: 0 tracked 3-D horizon, 1+n line horizon n, 1000+n surface n (time grid)
    std::vector<LcRow> lc;
    LayerCakeResult lcRes; bool lcHas = false;
    MistieFit mistie; bool mistieHas = false; int mistieRef = 0;
    // C43 P2 property templates (name match -> colour map, range, log scale) used by the map, sections, intersection and
    // 3-D scene; P7 cell filter (zones / fault blocks / facies switched off) shared by every window and the volumetrics
    struct PropTemplate { std::string name, match, unit; int cmap = 1; bool reverse = false, fixed = false, logScale = false; double lo = 0, hi = 1; };
    std::vector<PropTemplate> templates; bool templatesInit = false, showTemplates = false;
    struct CellFilter { bool on = false; std::vector<char> zoneOff, blockOff, faciesOff; } cf;
    bool volUseFilter = true;
    // C43 P5 well section: tracks (curve, log scale), wells hidden, flattening top, top picking
    struct WsTrack { std::string curve; bool logX = false; float width = 1.0f; bool useTpl = true; };   // C44 width, template
    std::vector<WsTrack> wsTracks; std::vector<std::string> wsHidden;
    std::vector<std::string> wsOrder;   // C44: well order in the section (names; wells not listed follow in project order)
    // C46 object-tree folders and the project file (.rgeo); cubePath: file of the 3-D volume (re-opened with the project)
    ressim::geo::ObjFolders folders;
    std::string projectPath, cubePath;
    bool treeExpand = false;            // C46: expand every tree node once (smoke tests / "Expand all")
    std::string wsFlat, wsPickTop; bool wsPick = false;
    int cubeRev = 0;     // bumped whenever the in-memory cube is replaced in place (section / attribute cache keys)
    int openTool = -1;   // smoke tests / links: open this C43 panel once (0 polygons, 1 mis-tie, 2 seismic calculator, 3 layer cake)
    int editHz = 0, brushMode = 0;                // brush: 0 raise / lower, 1 smooth, 2 flatten
    double brushR = 300, brushAmt = 5, brushStrength = 0.5, flatDepth = 0;   // flatDepth 0: depth under the first click
    int smoothPasses = 2;
    double minThick = 0.5, isoThick = 20;
    std::vector<std::vector<Surface>> hzUndo;     // whole horizon set before each edit (max 30)
    void pushUndo() { hzUndo.push_back(hz); if (hzUndo.size() > 30) hzUndo.erase(hzUndo.begin()); }
    std::shared_ptr<TrainingImage> ti;            // MPS training image
    std::string tiName;
    int tiNx = 200, tiNy = 200; unsigned tiSeed = 7;
    double tiWidth = 6, tiWave = 80, tiAmp = 12, tiFrac = 0.3;
    MpsSpec mps;
    UncSpec unc;
    UncInputs uncIn;                              // inputs of the last loop (realisations are rebuilt from them)
    UncResult ur;
    bool hasUnc = false;
    PoroPerm userFit; bool useUserFit = false;    // C40 poro-perm transform from the Data-analysis cross-plot (replaces the per-facies fit)
    void msg(const std::string& m) { log.push_back(m); if (log.size() > 400) log.erase(log.begin()); }
    double zMid() const {
        double s = 0; int n = 0;
        for (const auto& h : hz) for (double v : h.z) if (std::isfinite(v)) { s += v; ++n; }
        return n ? s / n : 0;
    }
};

// ---- seismic source: the out-of-core SEG-Y at full resolution when open, else the in-memory cube
struct SeisGeom {
    int ni = 0, nj = 0, ns = 0, il0 = 1, xl0 = 1, dil = 1, dxl = 1;
    double t0 = 0, dt = 4, x0 = 0, y0 = 0, ix = 1, iy = 0, jx = 0, jy = 1;
    void xy(double i, double j, double& x, double& y) const { x = x0 + i * ix + j * jx; y = y0 + i * iy + j * jy; }
    bool ij(double x, double y, double& i, double& j) const;
    double tEnd() const { return t0 + (ns - 1) * dt; }
};
bool seisGeom(const App::GeoState& S, SeisGeom& g);                       // false: no seismic
SeisGeom cubeGeom(const SeismicCube& c);
std::vector<float> seisSlice(const App::GeoState& S, int kind, int idx);  // kind 0 inline [j][s], 1 crossline [i][s], 2 time sample [i][j]
// attribute section (kind 0 / 1) with the neighbouring slices for the multi-trace attributes
std::vector<float> seisAttrSection(const App::GeoState& S, int kind, int idx, SeisAttr a);
void seisAttrRange(SeisAttr a, const std::vector<float>& v, double clipPct, float& lo, float& hi);   // colour range
double trackedTimeAt(const App::GeoState& S, double x, double y);       // tracked horizon TWT at a plan point (NaN)
std::vector<double> geoGridProperty(const App::GeoState& S, int which);
// C43: background job (queued behind a running one); the result is applied on the UI thread
void geoLaunch(App::GeoState& S, const std::string& name, std::function<std::function<void(App::GeoState&)>()> f);
// C36-4: apply the horizon brush at a map point (undo entry, cache invalidation); returns the nodes changed
int geoApplyBrush(App::GeoState& S, double x, double y);
// world coordinates of a fractional lattice position of a surface (rotation aware)
void geoLatticeXY(const Surface& h, double fi, double fj, double& x, double& y);
extern const char* kGeoGridProps[7];
// C42: built-in properties (kGeoGridProps) followed by the calculator's user properties
int geoGridPropCount(const App::GeoState& S);
std::string geoGridPropName(const App::GeoState& S, int which);
bool geoGridPropCombo(const char* id, int* which, const App::GeoState& S);
