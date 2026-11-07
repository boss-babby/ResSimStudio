// C14: hydraulic fracture propagation (stimulation sizing). Planar fluid-driven fracture from a wellbore, one wing
// symmetric, Newtonian fluid, lubrication flow in the fracture, linear elastic rock:
//   PKN (Perkins-Kern-Nordgren, fixed height H, L >> H): local elasticity  p_net = E' w / (2 H)  (w = maximum width of
//        the elliptical cross-section), area (pi/4) H w, flow q = -(pi H w^3 / 64 mu) dp/dx
//   KGD (plane strain, H >> L): non-local elasticity by constant displacement-discontinuity elements (both wings),
//        p_i = sum_j C_ij w_j,  C_ij = E'/(4 pi) [1/(x_i - x_j + a) - 1/(x_i - x_j - a)] (+ mirror wing), flow per height
//        q = -(w^3 / 12 mu) dp/dx
// Mass balance  d(A)/dt + dq/dx + q_L = 0 with Carter leak-off q_L = 2 h C_L / sqrt(t - tau(x)) (tau: exposure time),
// solved implicitly (Newton, backward Euler) on a fixed grid; the front is the degenerate-diffusion front (w -> 0) of
// the zero-toughness (viscosity-dominated) solution. Toughness: the next cell opens only when the stress intensity at
// the tip reaches K_IC (PKN: K_I = p_net sqrt(pi H / 2); KGD: displacement correlation w = 8 K_I sqrt(r) / (E' sqrt(2 pi))).
// Verified against Nordgren (1972) PKN and the Adachi-Detournay (2002) KGD viscosity-dominated similarity solutions.
#pragma once
#include "ressim/case.h"
#include <array>
#include <string>
#include <vector>

namespace ressim {
namespace hf {

struct Input {
    int model = 0;              // 0 PKN, 1 KGD, 2 P3D (C20: PKN-type lateral flow with equilibrium height growth)
    double E = 30e9, nu = 0.25; // Pa, -
    double height = 30;         // m: PKN fracture height; KGD: height over which the rate is spread
    double mu = 1e-3;           // fluid viscosity, Pa s
    double rate = 0.05;         // total injection rate into both wings, m3/s
    double cl = 0;              // Carter leak-off coefficient, m / sqrt(s)
    double kic = 0;             // fracture toughness, Pa sqrt(m)
    double sigmaMin = 0;        // minimum principal stress (Pa) - reported bottom-hole pressure = sigmaMin + p_net(0)
    double tEnd = 600;          // pumping time, s
    double maxLength = 500;     // computational domain per wing, m
    int cells = 200;            // grid cells per wing
    int outputs = 20;           // snapshots (log-spaced in time)
    // C20 P3D (model 2): 'height' is the pay-zone thickness h; symmetric bounding layers with closure stress
    // sigmaMin + stressContrast above and below. Each cross-section takes the equilibrium height H >= h at which the
    // stress intensity at both vertical tips equals K_IC (Simonson et al. 1978, Fung et al. 1987):
    //   p_net = K_IC / sqrt(pi H / 2) + (2 dsigma / pi) acos(h / H)                       (H > h; H = h while contained)
    //   A     = (pi H^2 / (2 E')) [p_net - (2 dsigma / pi) (acos r - r sqrt(1 - r^2))],  r = h / H
    // (Betti reciprocity with the opening of a uniformly loaded crack); lateral flow q = -(A^3 / (pi^2 H^2 mu)) dp/dx
    // (elliptic-equivalent section), Carter leak-off over the full height. Large dsigma recovers PKN exactly.
    double stressContrast = 0;  // Pa (model 2)
    double maxHeight = 0;       // m, cap of the height growth (0 = 20 h); reaching it flags runaway height growth
};

struct Snapshot {
    double t = 0, L = 0, w0 = 0, pnet0 = 0, eff = 1;   // half-length, wellbore width, net pressure, fluid efficiency
    std::vector<double> x, w, p;                        // profiles along the open part (m, m, Pa)
    std::vector<double> H;                              // P3D: local fracture height along the open part (m)
};

struct Result {
    std::vector<Snapshot> snaps;
    double volInjected = 0, volFracture = 0, volLeak = 0;
    int steps = 0, newton = 0;
    double height = 0;          // fracture height used (m); P3D: length-averaged height of the last snapshot
    double heightWell = 0;      // P3D: height at the wellbore (m) at the end
    bool heightRunaway = false; // P3D: the height reached maxHeight (weak barriers: the containment assumption failed)
    std::string error;
    const Snapshot& last() const { return snaps.back(); }
    double meanWidth() const;   // fracture volume / (2 L H) of the last snapshot (PKN: (pi/4) of the mean max width)
};

bool simulate(const Input& in, Result& out);
// same with the domain length chosen automatically (in.maxLength ignored; in.cells <= 0 picks a default)
bool simulateAdaptive(const Input& in, Result& out);

// C20 P3D section relations (SI): equilibrium height and net pressure of a cross-section with area A
struct P3DSection { double H = 0, pnet = 0, dpdA = 0, dHdA = 0; bool capped = false; };   // with the analytic slopes along the curve
P3DSection p3dSection(double A, double h, double dsigma, double kic, double Ep, double Hmax);
double p3dArea(double H, double h, double dsigma, double kic, double Ep);   // A on the equilibrium curve (H >= h)

// similarity solutions (zero toughness, no leak-off)
double nordgrenLength(const Input& in, double t);   // PKN half-length
double nordgrenWidth(const Input& in, double t);    // PKN maximum width at the wellbore
double kgdLength(const Input& in, double t);        // KGD half-length (Adachi & Detournay 2002, gamma = 0.6152)
double kgdWidth(const Input& in, double t);         // KGD wellbore opening (w0^2 = 2.1 sqrt(12 mu Q / E') L)

// ---- C15 multistage hydraulic fractures for the black-oil / compositional grid (GridInput::hfracs, EDFM)
struct StageLayout {
    std::string well;                   // well whose completions define the lateral (heel -> toe in completion order)
    int stages = 10;
    double spacing = 0;                 // case length units; <= 0: stages spread evenly over the completed length
    double toeOffset = 0;               // distance of the first stage from the toe (spacing > 0)
    double halfLength = 100, height = 30;   // case length units
    double width = 0.003;               // propped width, case length units
    double perm = 0, poro = 0.35;       // proppant pack permeability (mD; <= 0: cubic law) and porosity
    double gamma = 0;                   // stress-dependent conductivity (1 / pressure units)
    double strikeDeg = -1;              // fracture strike (degrees clockwise from +y); < 0: perpendicular to the well
    std::string prefix = "HF";
};
// trajectory points (cell centres of the completions) of a well; empty when the well has no completion in the grid
std::vector<std::array<double, 3>> wellPath(const Case& c, const std::string& well);
// HFRAC records of a multistage design (throws std::runtime_error with a readable message)
std::vector<GridInput::HFrac> layoutStages(const Case& c, const StageLayout& L);
// propped width from the pumped proppant mass (kg) over both wings (2 L_p H): w_p = M / (rho_p (1 - phi_p) 2 L_p H)  (m)
inline double proppedWidth(double massKg, double rhoP, double phiP, double proppedHalfLengthM, double heightM) {
    return massKg / (rhoP * (1 - phiP) * 2 * proppedHalfLengthM * heightM);
}

// ---- C17 simultaneous multi-cluster growth with stress shadowing (one stage, PKN clusters, shared pump rate)
// Each cluster k is a PKN fracture (zero toughness, Carter leak-off) fed by q_k from one wellbore pressure p_w:
//   p_w = S_hmin + ds_k + dsShadow_k + p_net,k(0) + K_perf q_k^2,   sum_k q_k = Q   (fully implicit with the fractures)
// Stress shadow (Sneddon 2D crack under uniform pressure, along the normal through its centre):
//   dsShadow_k = sum_{j != k} pbar_j [1 - xi^3 / (1 + xi^2)^(3/2)],  xi = |s_k - s_j| / c_j,  c_j = min(H, 2 L_j) / 2,
// with pbar_j the mean net pressure of cluster j (explicit, start of the step). Interior clusters see the shadow of both
// neighbours, take less fluid and stay shorter; limited-entry perforation friction evens the distribution out.
struct ClusterInput {
    Input base;                          // E, nu, height, mu, rate (total stage rate), cl, tEnd, maxLength, cells (model must be PKN)
    int clusters = 5;
    double spacing = 10;                 // m between neighbouring clusters
    std::vector<double> stressOffset;    // Pa per cluster added to S_hmin (heterogeneity), empty = 0
    double perfFriction = 0;             // K_perf, Pa / (m3/s)^2 per cluster (perfFrictionCoef)
    bool shadow = true;
};
struct ClusterResult {
    std::vector<double> t, pw;                          // output times (s), wellbore pressure above S_hmin (Pa)
    std::vector<std::vector<double>> q, L;              // [output][cluster] rate (m3/s) and half-length (m)
    std::vector<double> Lf, w0f, qf, shadowf, volf;     // final per cluster: half-length, wellbore width, rate, shadow stress (Pa), stored volume (m3)
    double volInjected = 0, volFracture = 0, volLeak = 0;
    int steps = 0, newton = 0;
    std::string error;
};
bool simulateClusters(const ClusterInput& in, ClusterResult& out);
// perforation friction coefficient: dp = K q^2, K = 8 rho / (pi^2 n^2 d^4 Cd^2)  (rho kg/m3, n perforations, d m)
inline double perfFrictionCoef(double rho, int nPerf, double dPerf, double cd = 0.6) {
    const double pi = 3.14159265358979323846, n = double(nPerf > 0 ? nPerf : 1);
    return 8 * rho / (pi * pi * n * n * dPerf * dPerf * dPerf * dPerf * cd * cd);
}

} // namespace hf
} // namespace ressim
