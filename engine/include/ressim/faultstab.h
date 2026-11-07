// C15: fault / fracture stability and fracture movement (geomechanical analysis of planes in a 3-D stress field).
//
// Stress state (compression positive, x East, y North, z down): principal stresses S_v (vertical), S_H (maximum
// horizontal, azimuth theta clockwise from North) and S_h; pore pressure p, Biot coefficient alpha:
//   sigma = S_H a a^T + S_h b b^T + S_v e_z e_z^T,  a = (sin theta, cos theta, 0), b = (cos theta, -sin theta, 0)
// A plane (strike phi clockwise from North, dip delta to the right of the strike, right-hand rule) has the normal n,
// traction t = sigma n, normal stress sigma_n = n.t, shear stress tau = |t - sigma_n n|, effective sigma_n' = sigma_n - alpha p.
//   slip tendency       T_s = tau / sigma_n'                         (Morris et al. 1996)
//   dilation tendency   T_d = (sigma_1 - sigma_n) / (sigma_1 - sigma_3)  (Ferrill et al. 1999)
//   Coulomb failure     CFF = tau - mu sigma_n' - c ;  critical pore pressure  p_c = (sigma_n - (tau - c) / mu) / alpha
// Poroelastic stress path of a laterally extensive reservoir (Engelder & Fischer 1994; Segall 1989):
//   d S_h = d S_H = A dp,  A = alpha (1 - 2 nu) / (1 - nu),  d S_v = 0;  thermal  d S_h = E alpha_T dT / (1 - nu)
// Fracture movement: Barton-Bandis normal closure  dv_n = sigma_n' V_m / (K_n0 V_m + sigma_n'), shear dilation
// w_dil = min(delta tan(psi), w_dil,max), aperture w = w_res + V_m - dv_n + w_dil, permeability multiplier (w / w_ref)^3; and
// quasi-dynamic spring-slider slip with rate-and-state friction (Dieterich 1979, Ruina 1983 aging law, radiation damping
// eta = G / (2 v_s), Rice 1993):
//   tau_0 + d tau(t) - k delta - eta V = sigma_n'(t) [ f_0 + a ln(V / V_0) + b ln(V_0 theta / D_c) ],  d theta/dt = 1 - V theta / D_c
// with the patch stiffness k = 7 pi G / (16 r) (circular patch of radius r, Eshelby). Instability when k < k_c =
// sigma_n' (b - a) / D_c; seismic events (V > V_seis) report slip, seismic moment M_0 = G pi r^2 d and
// M_w = (2/3)(log10 M_0 - 9.1).
#pragma once
#include <string>
#include <vector>

namespace ressim {
namespace stab {

struct StressState {
    double Sv = 0, SH = 0, Sh = 0;      // Pa (compression positive)
    double azimuthDeg = 0;              // S_H azimuth, degrees clockwise from North (+y)
    double p = 0;                       // pore pressure, Pa
    double biot = 1;
};
struct Plane { double strikeDeg = 0, dipDeg = 60; };
struct Strength { double mu = 0.6, cohesion = 0; };

struct Resolved {
    double sn = 0, tau = 0, snEff = 0;  // Pa
    double slipTend = 0, dilTend = 0;   // -
    double cff = 0;                     // Pa (> 0: slip)
    double pCrit = 0, dpCrit = 0;       // critical pore pressure and its margin from p (Pa; < 0: already critical)
    double rakeDeg = 0;                 // slip direction in the plane (degrees from the strike, + up-dip... Aki-Richards)
    double n[3] = {0, 0, 0};            // unit normal (x E, y N, z down)
};
Resolved resolve(const StressState& s, const Plane& pl, const Strength& st);
// stress state at depth from gradients (Pa/m) and the S_H azimuth
StressState fromGradients(double depth, double svGrad, double sHGrad, double shGrad, double azimuthDeg, double p, double biot);
// poroelastic (and thermal) change of the horizontal stresses for a pressure / temperature change
void stressPath(StressState& s, double dp, double nu, double dT = 0, double E = 0, double alphaT = 0);
inline double stressPathCoeff(double biot, double nu) { return biot * (1 - 2 * nu) / (1 - nu); }
// slip tendency of all orientations (lower-hemisphere poles): strike 0..360 by dStrike, dip 0..90 by dDip;
// row-major [dip][strike]
std::vector<double> slipTendencyGrid(const StressState& s, double dStrike, double dDip, int& nStrike, int& nDip);

// ---- Barton-Bandis + dilation aperture of a fracture under effective normal stress and slip
struct BBParams { double w0 = 0.05e-3, Kn0 = 20e9, Vm = 0.4e-3, psiDeg = 5, dilMax = 2e-3; };   // w0 = residual aperture (m), Pa/m, m, deg, m
double bbAperture(const BBParams& b, double snEff, double slip);

// ---- rate-and-state spring slider
struct RsfInput {
    double a = 0.010, b = 0.015, Dc = 1e-4, f0 = 0.6, V0 = 1e-6;   // friction (-, -, m, -, m/s)
    double G = 20e9, vs = 3000;         // shear modulus (Pa), shear-wave speed (m/s) -> radiation damping G / (2 vs)
    double radius = 100;                // slipping patch radius (m)
    double Vinit = 1e-15;               // initial slip rate (m/s): locked fault
    double Vseis = 1e-3;                // seismic slip-rate threshold (m/s)
    double tau0 = 0, sn0 = 0;           // initial shear / effective normal stress on the plane (Pa)
    // loading history (piecewise linear in time, s): d tau(t) and the effective normal stress change d sigma_n'(t)
    std::vector<double> t, dTau, dSn;
    double tEnd = 0;                    // <= 0: last history time
    int maxSteps = 400000;
    BBParams bb;                        // aperture evolution of the plane
};
struct RsfSample { double t, slip, V, tau, snEff, theta, mu, aperture; };
struct RsfEvent { double tStart = 0, tEnd = 0, slip = 0, Vmax = 0, M0 = 0, Mw = 0; };
struct RsfResult {
    std::vector<RsfSample> samples;
    std::vector<RsfEvent> events;
    double kPatch = 0, kCrit = 0;       // patch stiffness and critical stiffness at the initial state (Pa/m)
    int steps = 0;
    std::string error;
};
bool simulateRsf(const RsfInput& in, RsfResult& out);

} // namespace stab
} // namespace ressim
