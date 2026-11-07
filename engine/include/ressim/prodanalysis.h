// C28 production analysis (CMG CoFlow / IPSM-style engineering tools on the case and its results):
//   decline-curve analysis   Arps (1945) exponential / hyperbolic / harmonic fit of a rate history (log-rate least
//                            squares, b in [0, bMax]), forecast, time to the economic limit and EUR
//   nodal analysis           composite Vogel IPR (straight line above the bubble point), vertical lift performance of
//                            the tubing by a homogeneous no-slip multiphase model on the case PVT (hydrostatic + Haaland
//                            friction, marched in segments), operating point = IPR / VLP intersection, and generation of a
//                            VFPPROD table from the same model (flow x THP x WCT x GOR) for the simulator
// Limits (stated in the UI): the no-slip model has no liquid hold-up / flow-pattern map (Hagedorn-Brown / Beggs-Brill
// give higher BHP at low rates in gassy wells); no surface network solver (the group controls of the simulator are used).
#pragma once
#include "ressim/case.h"
#include "ressim/simulator.h"
#include <string>
#include <vector>

namespace ressim {
namespace prod {

struct ArpsFit {
    bool ok = false; std::string error;
    double qi = 0, Di = 0, b = 0;           // rate at t0 (rate units), nominal decline (1/day), exponent
    double t0 = 0;                          // start of the fitted period (days)
    double r2 = 0, rmseLog = 0;             // fit quality (rates / log rates)
    int points = 0;
};
// fit q(t) = qi / (1 + b Di (t - t0))^(1/b) to the points from the rate peak on (t in days, q > 0); fixB >= 0 fixes b
ArpsFit fitArps(const std::vector<double>& t, const std::vector<double>& q, double bMax = 2.0, double fixB = -1);
double arpsRate(const ArpsFit& f, double t);
double arpsCum(const ArpsFit& f, double t);                  // cumulative from t0
// time (days after t0) at which the rate falls to qEcon (capped at tMax) and the cumulative until then
double arpsEUR(const ArpsFit& f, double qEcon, double tMax, double* tEnd = nullptr);

// ---- nodal analysis (case units: pressure, length, surface rates, GOR in the case's Rs unit)
struct Ipr { double pr = 0, pb = 0, J = 0; };                // reservoir pressure, bubble point, PI above pb (rate / pressure)
double iprRate(const Ipr& ipr, double pwf);
double iprAof(const Ipr& ipr);                               // absolute open flow (pwf = 0)
struct Tubing {
    double depthTVD = 2000, lengthMD = 2000;                 // length units of the case
    double diameter = 0.0762;                                // inner diameter (case length units)
    double roughness = 1.5e-5;                               // absolute roughness (case length units)
    double tempGradientIgnored = 0;                          // isothermal (the PVT tables are isothermal)
    int segments = 60;
};
// bottom-hole pressure for an oil rate (surface) at the given THP, water cut and producing GOR (case PVT, region 0)
double vlpBhp(const Case& c, const Tubing& tb, double qOil, double thp, double wct, double gor);
struct OperatingPoint { bool ok = false; double q = 0, pwf = 0; std::string note; };
OperatingPoint operatingPoint(const Case& c, const Ipr& ipr, const Tubing& tb, double thp, double wct, double gor);
// VFPPROD table from the tubing model (flo = OIL, wfr = WCT, gfr = GOR, no ALQ); datum = depthTVD
VfpTable makeVfpProd(const Case& c, const Tubing& tb, int id, const std::vector<double>& rates, const std::vector<double>& thps,
                     const std::vector<double>& wcts, const std::vector<double>& gors);
// IPR of a well from a simulation state: PI from the last summary point (rate / (FPR - BHP)), bubble point from the PVT
Ipr iprFromResults(const Case& c, const std::vector<SummaryPoint>& summary, int well, std::string* note = nullptr);

}  // namespace prod
}  // namespace ressim
