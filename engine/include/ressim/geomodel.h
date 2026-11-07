// C32 geomodelling workflow (Petrel / RMS-style, the ResGeo workspace): data -> seismic -> structure -> 3-D grid ->
// well blocking -> facies -> petrophysics -> volumetrics -> upscaling -> simulation case.
//
// Data formats (industry exchange formats, read and written):
//   wells      : LAS 2.0 logs (~V ~W ~C ~A sections, NULL value, wrapped / unwrapped), well heads (CSV / whitespace:
//                name x y kb), deviation surveys (Petrel well-path ASCII 'MD X Y Z TVD DX DY AZIM INCL ...' with a
//                '# WELL NAME:' header, or MD INC AZI columns - minimum curvature), well tops (well surface MD)
//   surfaces   : IRAP classic ASCII grid (ROXAR / Petrel), ZMAP+ grid, scattered XYZ points (gridded by kriging)
//   faults     : fault sticks as 'x y z stick [fault]' columns or x y z blocks separated by blank / -999 lines
//   seismic    : SEG-Y rev 1 (IBM or IEEE float, inline / crossline / CDP X / CDP Y from the trace headers, coordinate
//                scalar) - 3-D post-stack volumes
// Depths are positive down (TVD SS), x east, y north, in the project length unit (m or ft).
//
// Structure: horizons (top -> base) + faults -> pillar grid with vertical pillars on a regular (optionally rotated)
// plan lattice; faults are stair-stepped along cell faces: every column belongs to a fault block (side of each fault
// trace), the corner depths of a column are sampled inside its own block (no smearing across the fault), so the ZCORN
// split nodes carry the throw; the faces between blocks become FAULTS entries (MULTFLT). Zones between horizons are
// layered proportionally, parallel to the top or parallel to the base. Horizons can be tied to the well tops (kriged
// residual correction).
// Properties: well-log blocking (arithmetic / geometric / harmonic / most-of), facies by sequential indicator simulation,
// truncated Gaussian simulation or object-based channels; petrophysics by SGS per facies with optional collocated
// co-kriging to a seismic attribute (Markov-1), permeability from a fitted poro-perm transform or co-simulated.
// Volumetrics: bulk / net / pore / hydrocarbon pore volume and STOIIP / GIIP per zone, facies, fault block with contacts
// and a saturation-height function, over realisations (P90 / P50 / P10). Upscaling to a coarser simulation grid:
// pore-volume weighted porosity, bulk-weighted NTG, permeability arithmetic / harmonic / geometric or flow based
// (local single-phase TPFA solve per coarse block, no-flow sides).
#pragma once
#include "ressim/case.h"
#include "ressim/geostat.h"
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <functional>
#include <memory>

namespace ressim {
namespace geo {

using P3 = std::array<double, 3>;

// ---------------------------------------------------------------- surfaces
struct Surface {
    std::string name;
    int nx = 0, ny = 0;
    double x0 = 0, y0 = 0, dx = 1, dy = 1, rotDeg = 0;   // origin = node (0,0); rotation counter-clockwise about it
    std::vector<double> z;                                // nx*ny, i fastest; NaN undefined
    double node(int i, int j) const { return z[size_t(i + nx * j)]; }
    void nodeXY(int i, int j, double& x, double& y) const;
    double at(double x, double y) const;                  // bilinear; NaN outside / undefined
    double atBlock(double x, double y, double tx, double ty) const;   // sample, falling back toward (tx, ty) if undefined
    bool valid() const { return nx > 1 && ny > 1 && z.size() == size_t(nx) * size_t(ny); }
    void range(double& lo, double& hi) const;
};
Surface readIrap(const std::string& path, std::string& err);
bool writeIrap(const std::string& path, const Surface& s);
Surface readZmap(const std::string& path, std::string& err);
std::vector<P3> readXyz(const std::string& path, std::string& err);
// grid scattered points on the lattice of 'tmpl' (ordinary kriging with the variogram, or inverse distance when vg.rangeMax <= 0)
Surface gridPoints(const std::vector<P3>& pts, const Surface& tmpl, const geostat::Variogram& vg, int maxNbr = 16);

// ---------------------------------------------------------------- wells
struct LogCurve { std::string mnem, unit, descr; std::vector<double> v; bool discrete = false; };
struct WellLog {
    std::string well, file;
    double nullValue = -999.25;
    std::vector<double> md;                       // index curve (DEPT / MD)
    std::vector<LogCurve> curves;                 // excluding the index curve
    const LogCurve* curve(const std::string& mnem) const;
    double valueAt(const LogCurve& c, double md) const;   // linear between samples, NaN for nulls / outside
};
WellLog readLas(const std::string& path, std::string& err);
bool writeLas(const std::string& path, const WellLog& log);
struct Station { double md = 0, inc = 0, azi = 0, x = 0, y = 0, tvd = 0; };   // inc / azi degrees (azimuth from north, clockwise)
struct Top { std::string surface; double md = 0; };
struct GeoWell {
    std::string name;
    double x = 0, y = 0, kb = 0;                  // head position, kelly-bushing elevation (TVD = TVDKB - kb)
    std::vector<Station> traj;                    // computed x, y, tvd (minimum curvature); empty = vertical
    std::vector<WellLog> logs;
    std::vector<Top> tops;
    int trajSrc = 0;                              // C42: 0 MD/INC/AZI (recomputed from the head), 1 absolute X/Y/TVDSS survey
    double td = 0;                                // C42: total depth from the head file / manual entry (0 = unknown)
    P3 at(double md) const;                       // position (x, y, tvd) at a measured depth
    double mdAtTvd(double tvd) const;             // first crossing (NaN if never reached)
    double totalMd() const;
};
void minimumCurvature(GeoWell& w);               // fills x, y, tvd of traj from md / inc / azi
std::vector<GeoWell> readWellHeads(const std::string& path, std::string& err);
// Petrel well-path file or MD INC AZI columns; well name from the file header or 'defaultName'
bool readDeviation(const std::string& path, GeoWell& w, std::string& err);
std::vector<std::pair<std::string, Top>> readTops(const std::string& path, std::string& err);

// ---------------------------------------------------------------- faults
struct Fault {
    std::string name;
    std::vector<std::vector<P3>> sticks;          // each stick: points from top to bottom
    double throwHint = 0;                         // added to horizons on the hanging-wall side when > 0 (smooth horizons)
    double mult = 1.0;                            // MULTFLT of the generated fault faces
    std::vector<std::array<double, 2>> traceAt(double depth) const;   // plan polyline at a depth (stick interpolation)
};
std::vector<Fault> readFaultSticks(const std::string& path, std::string& err);

// ---- C41 fault network (Petrel-style fault modelling): fault-fault junctions and truncation rules.
// detectJunctions compares the plan traces at depth 'zRef': X = the traces cross with both faults continuing beyond
// the crossing; Y = a fault tip lies on (or within 'tol' of) the other trace - the shorter / ending fault is the minor
// one; a tip short of the major by <= tol is a gap that truncation closes (lambda junctions are Y junctions whose
// minor fault is drawn from the major's side). truncateFaults clips every minor fault at its major on one side
// ('keepSide' +1 / -1 relative to the major's trace direction, 0 = the side holding most of the minor's trace):
// sticks on the other side are removed and a new terminal stick is built from the trace intersections at a series
// of depths, so the minor ends exactly on the major surface at every depth (closed fault blocks in buildGrid).
struct FaultJunction {
    int major = -1, minor = -1;
    char type = 'Y';                  // 'Y' (minor ends at major, incl. lambda), 'X' (crossing)
    double x = 0, y = 0;              // junction position at zRef
    double gap = 0;                   // distance from the minor tip to the major trace (0 = touching / crossing)
    int keepSide = 0;                 // truncation side of the minor relative to the major (0 = automatic)
    bool truncate = true;             // apply this rule in truncateFaults
};
std::vector<FaultJunction> detectJunctions(const std::vector<Fault>& faults, double zRef, double tol);
struct FaultNetReport { int truncated = 0, sticksRemoved = 0, sticksAdded = 0; std::vector<std::string> notes; };
std::vector<Fault> truncateFaults(const std::vector<Fault>& faults, const std::vector<FaultJunction>& rules, double tol, FaultNetReport* rep = nullptr);
// horizon-fault contacts: throw (hanging wall minus footwall depth, sampled 'offset' each side, positive = normal
// throw where the hanging wall is deeper) at 'n' points along the trace at the horizon's depth
struct ThrowProfile { std::vector<double> s, throwV, x, y; double maxThrow = 0; bool reverses = false, tipsTaper = true; };
ThrowProfile faultThrowProfile(const Fault& f, const Surface& horizon, double offset, int n = 25);

// ---------------------------------------------------------------- seismic
struct SeismicCube {
    std::string name, textHeader;
    int ni = 0, nj = 0, ns = 0;                   // inlines, crosslines, samples
    int il0 = 1, xl0 = 1, dil = 1, dxl = 1;
    double t0 = 0, dt = 4;                        // ms (two-way time) - or depth units for depth-migrated data
    double x0 = 0, y0 = 0;                        // position of (il0, xl0)
    double ix = 1, iy = 0, jx = 0, jy = 1;        // plan step per inline index / per crossline index
    int format = 5;                               // 1 IBM float, 5 IEEE float (as read)
    std::vector<float> a;                         // [i][j][s]
    float at(int i, int j, int s) const { return a[(size_t(i) * size_t(nj) + size_t(j)) * size_t(ns) + size_t(s)]; }
    void xy(double i, double j, double& x, double& y) const { x = x0 + i * ix + j * jx; y = y0 + i * iy + j * jy; }
    bool ij(double x, double y, double& i, double& j) const;   // fractional indices of a plan point
    float sampleT(double i, double j, double t) const;         // trilinear (NaN outside)
};
SeismicCube readSegy(const std::string& path, std::string& err);
bool writeSegy(const std::string& path, const SeismicCube& c);   // IEEE float, rev 1 headers
struct TimeHorizon { int ni = 0, nj = 0; std::vector<double> t; };   // NaN = not picked
// seeded auto-tracking of an extremum (mode 0 peak, 1 trough, 2 zero crossing + to -) within +-window ms, flood fill
TimeHorizon autoTrack(const SeismicCube& c, const std::vector<std::array<double, 3>>& seeds /* i j t */, int mode, double window, double maxJump);
struct VelocityModel { double v0 = 2000, k = 0.5, datum = 0; };   // interval velocity v(z) = v0 + k (z - datum); TWT ms <-> depth
double timeToDepth(const VelocityModel& v, double twtMs);
double depthToTime(const VelocityModel& v, double z);
Surface horizonToDepth(const SeismicCube& c, const TimeHorizon& h, const VelocityModel& v, const std::string& name);
// attribute per target point: 0 amplitude at the point, 1 RMS amplitude in +-win ms, 2 mean absolute amplitude
std::vector<double> extractAttribute(const SeismicCube& c, const VelocityModel& v, const std::vector<P3>& pts, int kind, double winMs);

// ---------------------------------------------------------------- structural grid
struct Zone { std::string name; int top = 0, base = 1; int nLayers = 5; int mode = 0; };   // horizon indices; mode 0 proportional, 1 parallel to top, 2 parallel to base
struct GridSpec {
    double x0 = 0, y0 = 0, rotDeg = 0, dx = 100, dy = 100;
    int nx = 20, ny = 20;
    std::vector<Zone> zones;
    double minThickness = 0.0;                    // cells thinner than this are inactive (pinch-out)
    bool tieToWells = true;                       // horizons corrected to the well tops
    double tieRange = 0;                          // variogram range of the tie correction (0 = 10 cells)
    // C36 pillars: 0 vertical, 1 fault-parallel (sloped near the faults: the pillar dip follows the nearest fault plane
    // with a Gaussian taper exp(-(d / range)^2) of the plan distance, vertical far away)
    int pillars = 0;
    double pillarRange = 0;                       // taper distance (0 = 4 cells)
};
struct GeoGrid {
    GridInput grid;                               // corner-point (coord / zcorn / actnum) + faults
    GridSpec spec;
    std::vector<int> layerZone;                   // zone per layer
    std::vector<int> block;                       // fault block per column (nx*ny)
    int nBlocks = 1;
    std::vector<P3> center;                       // cell centres (nx*ny*nz)
    std::vector<double> bulk;                     // cell bulk volumes
    std::vector<double> top, bot;                 // average top / bottom depth per cell
    std::vector<std::string> tieReport;           // well-tie residuals
    int cellOf(double x, double y, double z) const;   // natural index or -1
    void columnOf(double x, double y, double& fi, double& fj) const;   // fractional column indices
    void cornerXY(int i, int j, double& x, double& y) const;            // plan position of pillar (i, j) on the lattice
    void pillarXY(int i, int j, double z, double& x, double& y) const;  // C36: position of pillar (i, j) at depth z (COORD line; sloped pillars)
    double maxPillarSlope() const;                                      // C36: max |dxy/dz| of the pillars (0 = vertical)
};
GeoGrid buildGrid(const GridSpec& spec, std::vector<Surface> horizons, const std::vector<Fault>& faults,
                  const std::vector<GeoWell>& wells, std::string& err);
void computeGeometry(GeoGrid& g);                 // centres / volumes / tops (after zcorn changes)

// ---------------------------------------------------------------- blocking + property modelling
struct Blocked { int cell = -1; std::map<std::string, double> v; std::string well; double md = 0; };
// average every log curve over the trajectory length inside each cell; 'methods': mnem -> 0 arithmetic, 1 geometric,
// 2 harmonic, 3 most-of (discrete); unknown mnemonics: arithmetic (PERM*: geometric, FACIES / LITH*: most-of)
std::vector<Blocked> blockWells(const GeoGrid& g, const std::vector<GeoWell>& wells, const std::map<std::string, int>& methods = {});

struct FaciesSpec {
    int method = 0;                               // 0 SIS, 1 truncated Gaussian, 2 object channels
    std::vector<int> codes;                       // facies codes (from the blocked curve if empty)
    std::vector<double> proportions;              // target proportions (from the blocked data if empty)
    geostat::Variogram vg;                        // indicator / Gaussian variogram (sill normalised)
    std::string curve = "FACIES";
    unsigned seed = 1;
    // channels (method 2): background = codes[0], channel = codes[1]
    double chWidth = 200, chThick = 6, chAzimuth = 0, chAmplitude = 300, chWavelength = 2000;
    int maxNbr = 12;
};
struct FaciesResult { std::vector<int> facies; std::vector<double> proportions; std::vector<int> codes; int conditioned = 0, mismatches = 0; };
FaciesResult modelFacies(const GeoGrid& g, const std::vector<Blocked>& data, const FaciesSpec& spec);

struct PetroSpec {
    std::string curve = "PHIE";                   // blocked curve modelled
    int method = 1;                               // 0 kriging, 1 SGS
    geostat::Variogram vg;
    bool perFacies = true;
    std::vector<double> secondary;                // optional secondary per cell (seismic attribute) - collocated co-kriging
    double corrSecondary = 0;                     // correlation of the normal scores with the secondary
    double minV = 0, maxV = 1e300;
    bool logTransform = false;                    // model log10 (permeability)
    unsigned seed = 7;
    int maxNbr = 16;
};
struct PetroResult { std::vector<double> v; double mean = 0, sd = 0; int data = 0; std::string note; };
PetroResult modelProperty(const GeoGrid& g, const std::vector<Blocked>& data, const std::vector<int>* facies, const PetroSpec& spec);
// log10 k = a + b phi (least squares on the blocked pairs, per facies if given); returns per-cell permeability
struct PoroPerm { double a = 0, b = 0, r2 = 0; int n = 0; };
std::vector<PoroPerm> fitPoroPerm(const std::vector<Blocked>& data, const std::string& phi, const std::string& perm, const std::vector<int>* facies,
                                  const std::vector<int>& codes);
std::vector<double> permFromPoro(const std::vector<double>& phi, const std::vector<PoroPerm>& fits, const std::vector<int>* facies,
                                 const std::vector<int>& codes, double scatterSd, unsigned seed);

// ---------------------------------------------------------------- volumetrics
struct Contacts { double owc = 1e30, goc = -1e30; };   // depths; defaults: no water / no gas
struct SatHeight { double swirr = 0.2, pd = 2.0, lambda = 2.0, drhoOW = 0.2, drhoGO = 0.6; };   // Brooks-Corey: Sw = swirr + (1-swirr)(Pd / Pc)^lambda, Pc = drho h (pressure / length units)
struct VolSpec {
    Contacts contacts; SatHeight sh; bool satHeight = true; double swConst = 0.25;
    double bo = 1.2, bg = 0.005;                  // rb/stb | rm3/sm3 ; rcf/scf | rm3/sm3
    bool field = false;                           // field units: STOIIP in stb (bulk ft3 -> rb factor), GIIP in scf
    std::vector<char> mask;                       // C43: per cell, 0 = excluded (polygon boundary); empty = every cell
    std::vector<Contacts> blockContacts;          // C43: contacts per fault block (index = block); missing -> 'contacts'
};
struct VolRow { std::string key; double bulk = 0, net = 0, pore = 0, hcpvO = 0, hcpvG = 0, stoiip = 0, giip = 0; };
struct VolResult { std::vector<VolRow> byZone, byFacies, byBlock; VolRow total; };
VolResult volumetrics(const GeoGrid& g, const std::vector<double>& poro, const std::vector<double>& ntg, const std::vector<int>* facies,
                      const VolSpec& v);
struct Percentiles { double p90 = 0, p50 = 0, p10 = 0, mean = 0; };   // P90 = low (exceeded with 90 % probability)
Percentiles percentiles(std::vector<double> v);

// ---------------------------------------------------------------- upscaling + export
struct UpscaleSpec { int fx = 2, fy = 2, fz = 2; int permMethod = 3; };   // 0 arithmetic, 1 harmonic, 2 geometric, 3 flow based
struct Upscaled { GeoGrid g; std::vector<double> poro, ntg, kx, ky, kz; std::vector<int> facies; };
Upscaled upscale(const GeoGrid& fine, const std::vector<double>& poro, const std::vector<double>& ntg, const std::vector<double>& kx,
                 const std::vector<double>& ky, const std::vector<double>& kz, const UpscaleSpec& s);
// flow-based effective permeability of a box of fine cells (TPFA, unit pressure drop along 'dir', no-flow sides)
double flowBasedPerm(int nx, int ny, int nz, const std::vector<double>& k, double dx, double dy, double dz, int dir);
// geomodel -> simulation case grid (+ properties, faults); wells with trajectories become WELSPECS / COMPDAT on the
// cells they cross (producers by default)
void toCase(Case& c, const GeoGrid& g, const std::vector<double>& poro, const std::vector<double>& ntg, const std::vector<double>& kx,
            const std::vector<double>& ky, const std::vector<double>& kz, const std::vector<GeoWell>& wells);
bool writeGrdecl(const std::string& path, const GeoGrid& g, const std::map<std::string, const std::vector<double>*>& props);

// ---------------------------------------------------------------- demo project
// writes a synthetic but complete project (faulted anticline, 6 wells with LAS logs / surveys / tops, IRAP horizons,
// fault sticks, a SEG-Y volume) into 'dir'; returns the list of files
std::vector<std::string> writeDemoProject(const std::string& dir, bool field, std::string& err);
// headless workflow on the demo project (ResSimCLI --geo-demo): read -> grid (tied) -> blocking -> SIS facies -> SGS
// porosity per facies -> poro-perm -> volumetrics (+ 5 realisations) -> flow-based upscaling; writes
// GEOMODEL_FINE.GRDECL and GEOMODEL_COARSE.GRDECL into 'dir' and returns a text report (empty + err on failure)
std::string runDemoWorkflow(const std::string& dir, bool field, std::string& err);

// ================================================================ C36: structural editing, MPS, uncertainty loop
// ---- surface editing (interactive brush + whole-surface operations)
// Gaussian-tapered brush, weight w = exp(-3 r^2 / R^2) inside R: mode 0 raise / lower by 'amount' (depth units, +
// deeper), 1 smooth toward the 3x3 local mean (amount = strength 0..1), 2 set to the depth 'amount' (flatten).
// Returns the number of nodes changed.
int brushSurface(Surface& s, double x, double y, double radius, double amount, int mode);
void smoothSurface(Surface& s, int passes);                       // 3x3 binomial filter over the defined nodes
int fillSurfaceHoles(Surface& s);                                 // undefined nodes in-painted (Laplace); returns the count
int enforceThickness(std::vector<Surface>& hz, double minThick);  // horizon h >= h-1 + minThick (top-down); nodes moved
Surface isochoreSurface(const Surface& ref, double thickness, const std::string& name);   // conformable copy
bool surfaceIJ(const Surface& s, double x, double y, double& fi, double& fj);            // fractional node indices (rotation aware)

// ---- multiple-point statistics facies: Direct Sampling (Mariethoz, Renard & Straubhaar 2010, WRR 46 W11536)
struct TrainingImage { int nx = 0, ny = 0, nz = 1; std::vector<int> v; int ncat() const; int at(int i, int j, int k) const { return v[size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k))]; } };
// sinuous channels (code 1) in a background (0) along TI +x: width / wavelength / amplitude in TI cells, target fraction
TrainingImage channelTrainingImage(int nx, int ny, int nz, double width, double wavelength, double amplitude, double fraction, unsigned seed);
TrainingImage readTrainingImage(const std::string& path, std::string& err);   // GSLIB / SGeMS ASCII (title line 'nx ny nz', nvar, names, values)
bool writeTrainingImage(const std::string& path, const TrainingImage& ti);
struct MpsSpec {
    int maxNbr = 30;                 // conditioning neighbours (data + simulated)
    double threshold = 0.05;         // accepted mismatch fraction
    double maxScan = 0.3;            // max fraction of the TI scanned per node
    int maxScanAbs = 25000;          // absolute cap
    double servo = 0.5;              // servosystem strength (0 = plain DS): favour under-represented categories vs the target
    std::vector<double> target;      // target proportions per category (empty = TI marginal)
    double azimuthDeg = 0;           // map direction of the TI +x axis (degrees clockwise from north)
    double tiCell = 0;               // plan size of a TI pixel (0 = grid dx)
    double tiLayer = 1;              // grid layers per TI layer (3-D TIs)
    double radius = 15;              // search radius in TI cells (plan); vertical: radiusK grid layers
    int radiusK = 3;
    unsigned seed = 31;
};
// TI categories c = 0..K-1 map to the facies codes 'spec.codes[c]' (or c when no codes); conditioned on the blocked
// 'spec.curve' data
FaciesResult modelFaciesMps(const GeoGrid& g, const std::vector<Blocked>& data, const FaciesSpec& spec, const TrainingImage& ti, const MpsSpec& m);
// indicator connectivity: P(I(x)=1 and I(x+h)=1) / P(I=1) along a plan lag (di, dj) in index space, for category 'cat'
double indicatorContinuity(const std::vector<int>& v, int nx, int ny, int nz, int di, int dj, int cat);

// ---- uncertainty loop (volumetric Monte Carlo with structural, contact, petrophysical and stochastic uncertainty)
struct Tri {                         // triangular distribution (fixed when hi <= lo)
    double lo = 0, mode = 0, hi = 0;
    double sample(double u) const;   // inverse CDF
    bool fixed() const { return !(hi > lo); }
};
struct UncSpec {
    int nReal = 30;
    unsigned seed = 2026;
    double zSd = 0, zRange = 2000;   // structural: Gaussian random field (sd, practical range) added to all horizons, zero at the well tops
    Tri owc, goc;                    // contact depths (fixed: VolSpec value)
    Tri poroMult{1, 1, 1}, ntgMult{1, 1, 1}, sw, bo;   // multipliers; Sw (constant Sw or Swirr of the saturation-height function); Bo
    bool stochastic = true;          // new facies / porosity seeds per realisation
    bool latinHypercube = true;      // stratified sampling of the Tri inputs (lower Monte-Carlo noise for small nReal)
    int tornadoStochastic = 6;       // realisations of the stochastic tornado bar (0 = none)
};
struct UncInputs {
    GridSpec gs; std::vector<Surface> hz; std::vector<Fault> faults; std::vector<GeoWell> wells;
    bool facies = true; FaciesSpec fs; PetroSpec ps; VolSpec vs; std::vector<double> ntgConst;   // ntg per cell (empty = 1)
    std::shared_ptr<const TrainingImage> ti; MpsSpec mps;   // fs.method == 3: MPS facies from this training image
    bool ntgFromFacies = false;      // cells of facies code 0 (shale) get NTG 0 (the workspace's "NTG from facies")
};
struct UncSample { double owc = 0, goc = 0, poroMult = 1, ntgMult = 1, sw = 0, bo = 0, zShift = 0, zSd = 0, zRange = 2000; unsigned fieldSeed = 0, faciesSeed = 0, petroSeed = 0; bool field = false; };
struct UncReal { UncSample s; double stoiip = 0, giip = 0, pore = 0, bulk = 0; };
struct TornadoBar { std::string name; double low = 0, high = 0; };
struct UncResult {
    std::vector<UncReal> real; Percentiles stoiip, giip; double base = 0;
    std::vector<TornadoBar> tornado;   // sorted by swing (largest first)
    int iP90 = -1, iP50 = -1, iP10 = -1;   // realisations closest to the percentiles
    std::string report, err;
};
UncResult uncertaintyLoop(const UncInputs& in, const UncSpec& u, const std::function<bool(double)>& progress = {});
struct UncModel { GeoGrid g; std::vector<double> poro, ntg; FaciesResult fr; VolResult vol; std::vector<Surface> hz; };
UncModel realise(const UncInputs& in, const UncSample& s, std::string& err);   // deterministic rebuild of one realisation
// structural perturbation of a set of horizons: spectral Gaussian random field (sd, practical range), conditioned to 0 at the points
std::vector<Surface> perturbHorizons(const std::vector<Surface>& hz, double sd, double range, unsigned seed, const std::vector<std::array<double, 2>>& zeroAt, double shift = 0);

}  // namespace geo
}  // namespace ressim
