// C43 geomodel tools (Petrel parity items P4 / P6 / P8 / P9 of docs/PETREL_PARITY.md):
//   * 2-D mis-tie correction - least-squares bulk time shift per line from the line intersections
//   * polygons / boundaries  - objects, point-in-polygon, import / export, surface clipping, cell masks
//   * layer-cake velocity model - time surfaces + well tops -> interval velocities per layer (constant or map),
//                                 depth surfaces with residual correction at the wells; checkshot fitting
//   * seismic calculator / cropping - expressions over the volume and its attributes, sub-volumes
#pragma once
#include "ressim/geomodel.h"
#include "ressim/geo_import.h"
#include <array>
#include <string>
#include <vector>

namespace ressim {
namespace geo {

// ---------------------------------------------------------------- P6 mis-tie correction
struct MistieFit {
    std::vector<double> shift;           // per line, ms (add to the line's times)
    double rmsBefore = 0, rmsAfter = 0;  // over all tie / horizon pairs with picks on both lines
    int pairs = 0;                       // tie x horizon pairs used
    std::vector<int> component;          // connected component of each line (lines without ties: -1)
};
// minimises sum (t_a + s_a - t_b - s_b)^2 over the intersections of every horizon picked on both lines; the reference
// line of each connected component keeps s = 0 (refLine if it is in the component, else its first line)
MistieFit fitMisties(const std::vector<SeismicLine>& lines, const std::vector<LineTie>& ties, const std::vector<LineHorizon>& hz, int refLine);
// shifts the data (t0) and every pick of the lines
void applyLineShifts(std::vector<SeismicLine>& lines, std::vector<LineHorizon>& hz, const std::vector<double>& shift);

// ---------------------------------------------------------------- P4 polygons
struct Polygon {
    std::string name;
    std::vector<std::array<double, 2>> pts;
    bool closed = true;                  // boundary (closed) or open polyline (e.g. an intersection path)
};
bool insidePolygon(const Polygon& p, double x, double y);       // even-odd rule; open polylines: false
double polygonArea(const Polygon& p);                           // |shoelace|
double polylineLength(const Polygon& p);
// text: "x y [z]" rows; a new polygon starts after a blank line, a "999" / ">" separator or a "# name" line
// (the comment gives the name); also Petrel-style 4-column "x y z id" files (id change = new polygon)
std::vector<Polygon> readPolygons(const std::string& path, std::string& err);
std::vector<Polygon> parsePolygons(const std::string& text);
bool writePolygons(const std::string& path, const std::vector<Polygon>& polys);
// nodes outside (keepInside) or inside (!keepInside) the polygon become undefined
Surface clipSurface(const Surface& s, const Polygon& p, bool keepInside = true);
// per-cell mask of a grid (cell centre inside the polygon)
std::vector<char> polygonMask(const GeoGrid& g, const Polygon& p);

// ---------------------------------------------------------------- P8 layer-cake velocity model
struct LayerCakeInput {
    std::vector<Surface> twt;            // interface time surfaces (TWT ms), top down
    std::vector<std::string> topName;    // well-top name of each interface (calibration)
    std::vector<int> mode;               // per layer: 0 constant interval velocity (mean of the wells), 1 map (inverse distance), 2 user value
    std::vector<double> vUser;           // per layer user interval velocity (mode 2, or fallback without wells)
    double datumZ = 0;                   // depth of TWT 0 (seismic reference datum, depth positive down)
    bool residual = true;                // correct the depth surfaces to the well tops (inverse-distance residual maps)
};
struct LayerCakeWell { std::string well; int iface = 0; double x = 0, y = 0, zTop = 0, zBefore = 0, zAfter = 0; };
struct LayerCakeResult {
    std::vector<Surface> depth;          // per interface
    std::vector<Surface> vint;           // per layer interval velocity map (m/s or ft/s)
    std::vector<double> vMean;           // per layer mean interval velocity at the wells
    std::vector<int> nWells;             // per layer wells used
    std::vector<LayerCakeWell> ties;     // per well and interface: top depth, model depth before / after the residual
    double rmsBefore = 0, rmsAfter = 0;
    std::vector<std::string> warnings;
};
LayerCakeResult layerCakeDepth(const LayerCakeInput& in, const std::vector<GeoWell>& wells);
// time surface on the lattice of a cube (TimeHorizon -> Surface, z = TWT)
Surface timeSurface(const SeismicCube& c, const TimeHorizon& h, const std::string& name);
// picks of a line horizon as (x, y, TWT) points
std::vector<P3> lineHorizonTimePoints(const std::vector<SeismicLine>& lines, const LineHorizon& h, int step = 1);
// checkshot / time-depth pairs (depth, TWT ms) -> linear velocity function v = v0 + k (z - datum) (least squares on the
// time residuals, golden-section search on k); rms in ms
VelocityModel fitCheckshots(const std::vector<double>& z, const std::vector<double>& twt, double datum, double* rmsMs = nullptr);

// ---------------------------------------------------------------- P9 seismic calculator / cropping
// variables: AMP (the volume), ENVELOPE PHASE FREQUENCY COSPHASE SWEETNESS RMS RELIMPEDANCE COHERENCE VARIANCE (computed
// on demand with 'winSamples'), T (sample time), I J (0-based inline / crossline index); same operators as the log calculator
bool seismicCalculator(const SeismicCube& c, const std::string& expr, int winSamples, SeismicCube& out, std::string& err);
std::vector<std::string> seismicCalculatorNames();
// sub-volume [i0, i1] x [j0, j1] x [s0, s1] (inclusive) with decimation steps
SeismicCube cropCube(const SeismicCube& c, int i0, int i1, int j0, int j1, int s0, int s1, int di = 1, int dj = 1, int ds = 1);

// C44: where a fault cuts a plan path (intersection / 2-D line): at nz+1 depth levels in [z0, z1] (clipped to the stick
// depths) the fault's interpolated plan trace is intersected with the path; returns branches of (distance along path, depth)
std::vector<std::vector<std::array<double, 2>>> faultCutOnPath(const Fault& f, const Polygon& path, double z0, double z1, int nz);

}  // namespace geo
}  // namespace ressim
