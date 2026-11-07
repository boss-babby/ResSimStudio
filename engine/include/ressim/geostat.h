// C28 geostatistics (CMG Builder / Petrel-style property modelling, GSLIB conventions - Deutsch & Journel 1998):
//   variogram models (spherical / exponential / Gaussian, nugget, geometric anisotropy: major / minor horizontal range at
//   an azimuth + vertical range), experimental variograms, ordinary / simple kriging with a moving neighbourhood, and
//   sequential Gaussian simulation (normal-score transform, random path, simple kriging in Gaussian space on the data +
//   previously simulated nodes, back transform). applyToCase fills PORO / PERMX(Y,Z) / NTG of a case from hard data at
//   the well completions (the current values there) or an explicit point list.
#pragma once
#include "ressim/case.h"
#include <array>
#include <string>
#include <vector>

namespace ressim {
namespace geostat {

enum class VarioModel { Spherical = 0, Exponential = 1, Gaussian = 2 };
struct Variogram {
    VarioModel model = VarioModel::Spherical;
    double nugget = 0, sill = 1;                 // total sill = nugget + sill
    double rangeMax = 1000, rangeMin = 1000;     // horizontal ranges (practical range for exponential / Gaussian), case length units
    double rangeVert = 20;
    double azimuthDeg = 0;                       // direction of rangeMax, degrees clockwise from north (+y)
};
// semivariance for the separation (dx east, dy north, dz down)
double gamma(const Variogram& v, double dx, double dy, double dz);
inline double covariance(const Variogram& v, double dx, double dy, double dz) { return v.nugget + v.sill - gamma(v, dx, dy, dz); }

struct Point { double x = 0, y = 0, z = 0, v = 0; };
using Xyz = std::array<double, 3>;

struct Options {
    Variogram vg;
    int maxNbr = 16;                             // conditioning points per estimate (data + simulated nodes)
    double searchRanges = 2.0;                   // search radius in units of the (anisotropic) range
    bool simpleKriging = false; double skMean = 0;   // kriging only (SGS always uses simple kriging in Gaussian space)
    unsigned seed = 1;
    bool normalScore = true;                     // SGS: normal-score transform of the data (else the values are taken as Gaussian)
};
struct KrigeOut { std::vector<double> est, var; int noData = 0; };
KrigeOut krige(const std::vector<Xyz>& targets, const std::vector<Point>& data, const Options& opt);
// one SGS realisation at the targets (data located on a target node are honoured exactly)
std::vector<double> sgs(const std::vector<Xyz>& targets, const std::vector<Point>& data, const Options& opt);

struct ExpVario { std::vector<double> lag, gamma; std::vector<long> pairs; };
// experimental semivariogram of point values; direction 0 = horizontal (|dz| < lag/2, any azimuth), 1 = vertical,
// 2 = horizontal along azimuthDeg (+-22.5 deg)
ExpVario experimental(const std::vector<Point>& pts, double lagWidth, int nLags, int direction = 0, double azimuthDeg = 0, long maxPairs = 4000000);

struct ApplyInput {
    int property = 0;                            // 0 PORO, 1 PERMX (PERMY = PERMX, PERMZ = kvkh PERMX), 2 NTG
    int method = 1;                              // 0 kriging, 1 SGS
    Options opt;
    bool fromWells = true;                       // hard data = current values at the completed cells of every well
    std::vector<Point> data;                     // explicit hard data (added to the well data)
    bool logPerm = true;                         // PERMX modelled as log10 k
    double kvkh = 0.1;
    double minV = -1e300, maxV = 1e300;          // clamp of the result (property units)
};
struct ApplyReport { bool ok = false; std::string message; int cells = 0, dataPoints = 0; double mean = 0, stdDev = 0, minV = 0, maxV = 0; };
ApplyReport applyToCase(Case& c, const ApplyInput& in);
// C40: run a stored property model of the case (Case::propModels) - with the poro-perm transform for PORO models
ApplyReport applyPropModel(Case& c, const PropModelInput& m);
// C41: FACIES property model (property 3) by truncated Gaussian simulation - see PropModelInput
ApplyReport applyFaciesModel(Case& c, const PropModelInput& m);

}  // namespace geostat
}  // namespace ressim
