// C21 geomechanics G1-c: EDFM flow <-> embedded-fracture contact mechanics (staggered, lagged by one report step).
//
// Each coupling step takes the fracture-piece pressures of one grid layer (2-D plane-strain slice), solves the EFEM contact
// problem of efem.h with the along-fracture pressure profile (slip history carried over), and turns the contact state of
// every piece into a hydraulic aperture:
//   closed (effective normal stress s'_n = -t_n > 0): Barton-Bandis hyperbolic closure  a = a_res + V_m - V_m s'_n / (K_n0 V_m + s'_n)
//   open   (normal jump w_n > 0):                      a = a_res + V_m + w_n
//   shear dilation:                                    + min(|w_t| tan(psi), a_dil,max)
// The map-view opening is corrected for a finite fracture height H: w_n *= min(1, H / L_f) (PKN / KGD width ratio).
// and the fracture conductance multiplier (cubic law) m = (a / a_ref)^3, a_ref = the HFRAC width, applied to the
// fracture-fracture transmissibilities through Simulator::setHfConductance (harmonic mean of neighbouring pieces). The same
// multiplier is used for all layers of a piece column (plane strain).
// Limitations (documented, ENGINES_AND_MATH.md section 6.13): total-stress elasticity with fluid pressure on the fracture
// faces only (no poroelastic matrix coupling), one report-step lag (explicit staggering), 2-D slice.
#pragma once
#include "ressim/efem.h"
#include <atomic>
#include <functional>
#include <string>
#include <vector>

namespace ressim {
class Simulator;
namespace geocouple {

struct BartonBandis {
    double aRes = 2e-5;        // residual hydraulic aperture (m)
    double Vm = 3e-4;          // maximum closure (m)
    double Kn0 = 2e10;         // initial normal stiffness (Pa/m)
    double dilationDeg = 0;    // shear dilation angle psi (deg)
    double maxDilation = 5e-4; // cap of the dilation aperture (m)
};
// hydraulic aperture (m) from the effective normal stress (Pa, compression +), normal opening and slip (m)
double aperture(const BartonBandis& bb, double sigmaNeff, double opening, double slip);

struct Input {
    int layer = 0;                       // 0-based k of the slice whose pressures load the contact problem
    double E = 20e9, nu = 0.25;          // Pa, -
    double sxx = -40e6, syy = -30e6;     // far-field total in-plane stresses (Pa, tension +): grid x and y directions
    double friction = 0.6, cohesion = 0; // fracture Coulomb parameters
    int nx = 100, ny = 100;              // EFEM elements over the grid's x-y extent
    BartonBandis bb;
    double aRef = 0;                     // reference aperture (m); 0 = the HFRAC width of each fracture
    double penalty = 1e3;
    // PKN height correction of the plane-strain (map-view) opening: w *= min(1, H / L_f) (w_PKN / w_KGD = H / L_f).
    // height < 0: H from the fracture pieces (sum of the piece heights of one column); 0: off (infinite-height map view)
    double height = -1;
};

struct StepOut {
    double t = 0;                        // days
    // per fracture piece (all layers, same order as the hydraulic-fracture cells)
    std::vector<double> p, wn, wt, tn, a, mult;
    double maxOpening = 0, maxSlip = 0, minMult = 1, maxMult = 1;
    int openPieces = 0, slipPieces = 0, newton = 0;
};
struct Output {
    bool ok = false; std::string error;
    std::vector<StepOut> steps;
    std::vector<double> s;               // per piece: signed distance from its fracture midpoint (m)
    std::vector<int> frac;               // per piece: fracture index
    double efemSeconds = 0;
};

// One contact solve at the simulator's current state (pressures); updates prev (slip history) and the conductances.
bool coupleOnce(Simulator& sim, const Input& in, efem::Result& prev, StepOut& st, Output* geometry, std::string& err);
// Whole run: couple at t = 0, then after every report step until the schedule ends (sim must be set up + initialised).
bool runCoupled(Simulator& sim, const Input& in, Output& out, std::atomic<bool>* cancel = nullptr,
                const std::function<void(const std::string&)>& log = {});

}  // namespace geocouple
}  // namespace ressim
