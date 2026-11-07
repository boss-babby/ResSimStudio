// C23 geomechanics G3: rate-and-state friction on an embedded fault - quasi-dynamic earthquake cycles and injection-induced
// event catalogues (2-D plane strain, in-plane shear on a straight fault in an elastic full space).
//
// Fault discretised into n equal elements (length dx). Slip delta_i, slip rate V_i, state psi_i (Rice et al. 2001 form):
//   friction     tau_f = sigma' a asinh( V / (2 V0) exp(psi / a) )                       (regularised Dieterich-Ruina law)
//   aging law    dpsi/dt = (b V0 / Dc) (exp((f0 - psi) / b) - V / V0) - alpha dsigma'/dt / sigma'   (Linker-Dieterich term)
//   steady state psi_ss(V) = f0 - b ln(V / V0),  f_ss = f0 + (a - b) ln(V / V0)
//   elasticity   tau_i = tau0_i + sum_j K_ij (delta_j - Vpl t) + dtau_i(t)                (back-slip loading at Vpl)
//   quasi-dynamic balance  tau_i - eta V_i = tau_f(V_i, psi_i, sigma'_i),  eta = G / (2 c_s)  (radiation damping, Rice 1993)
//   effective normal stress sigma'_i(t) = sigma_n,i - p_i(t)  (compression +, floored at 0.1 MPa)
// Elastic kernel of a constant-slip edge-dislocation segment [xi-, xi+] on the fault line (plane strain, mode II):
//   K_ij = G / (2 pi (1 - nu)) [ 1 / (x_i - xi+_j) - 1 / (x_i - xi-_j) ]   (self term -G / (pi (1 - nu) dx/2))
// The balance is solved per element for the frictional stress s = tau - eta V with a bracketed Newton iteration (monotone),
// the ODE system (delta, psi) by an adaptive embedded Runge-Kutta 4(5) (Dormand-Prince) with error control on psi and slip.
// Events: max V above vEvent (default 1 mm/s) opens an event, below closes it; moment M0 = G * sum(dslip dx) * W with the
// out-of-plane width W = rupture length (square-rupture assumption of the 2-D model), Mw = 2/3 (log10 M0 - 9.1).
// Pore pressure / stress forcing: built-in line-source injection (Theis, superposed shut-in), a table (e.g. from the G1-c
// coupled EDFM run via fromCoupled) or a user callback. Initial normal / shear stress either uniform or resolved from an
// EFEM contact solve of efem.h along one embedded fracture (fromEfem) - the "embedded fault" of G1.
// Literature: Dieterich 1979, Ruina 1983, Rice 1993, Lapusta et al. 2000, Rubin & Ampuero 2005 (nucleation length),
// Linker & Dieterich 1992, Luo & Ampuero QDYN; injection-induced: McClure & Horne 2011, Dublanchet 2019.
#pragma once
#include <atomic>
#include <functional>
#include <string>
#include <vector>

namespace ressim {
namespace efem { struct Input; struct Result; }
namespace geocouple { struct Output; }
namespace rsf {

struct Input {
    // geometry (m); the fault is the segment (x0,y0)-(x1,y1), n elements
    double x0 = 0, y0 = 0, x1 = 2000, y1 = 0;
    int n = 128;
    // elasticity
    double G = 30e9, nu = 0.25, rho = 2700;        // Pa, -, kg/m3 (shear-wave speed c_s = sqrt(G / rho))
    // friction (uniform; a velocity-weakening patch between s0 and s1 (fractions of the length) uses aVW, the rest aVS)
    double aVW = 0.010, aVS = 0.025, b = 0.015;   // -, a < b weakening, a > b strengthening
    double vwFrom = 0.2, vwTo = 0.8;               // fraction of the fault length
    double Dc = 1e-3, f0 = 0.6, V0 = 1e-6;         // m, -, m/s
    double alpha = 0.0;                            // Linker-Dieterich normal-stress coupling (0 .. 0.5)
    std::vector<double> aProf;                     // optional per-element a (overrides aVW / aVS)
    // stresses (Pa, compression +)
    double sigmaN = 50e6;                          // total normal stress
    double p0 = 0;                                 // initial pore pressure
    double tauRatio = -1;                          // initial tau0 / sigma'; < 0: steady state at vInit (tau0 = f_ss(vInit) sigma')
    std::vector<double> sigmaNProf, tau0Prof;      // optional per-element values (e.g. from an EFEM contact solve)
    double Vpl = 1e-9;                             // back-slip loading rate (m/s); 0 = no tectonic loading
    double vInit = 1e-9;                           // initial slip rate for the steady state (m/s)
    double tauRate = 0;                            // extra uniform shear-stress rate (Pa/s)
    // pore-pressure forcing (all superposed)
    struct Injection {                             // line source (2-D Theis) at distance 'offset' from the fault, at fault fraction 'at'
        bool on = false;
        double at = 0.5, offset = 50;              // -, m
        double pStar = 1e6;                        // q mu / (4 pi k h) (Pa)
        double diffusivity = 0.05;                 // c = k / (mu phi c_t) (m2/s)
        double tStart = 0, tStop = 1e300;          // s (shut-in superposed with a negative source)
    } inj;
    std::vector<double> forceT;                    // table: times (s), per-time per-element dp (Pa) and dtau (Pa), linear in time
    std::vector<std::vector<double>> forceP, forceTau;
    std::function<void(double t, std::vector<double>& dp, std::vector<double>& dtau)> forcing;   // optional callback (adds)
    // run control
    double tEnd = 3.15e9;                          // s
    double vEvent = 1e-3;                          // seismic threshold (m/s)
    double rtol = 1e-7;
    double dtMax = 3.15e7;                         // s
    long maxSteps = 4000000;
    int snapshots = 400;                           // stored log10 V profiles (uniform in step count during events + in time)
};

struct Event {
    double tStart = 0, tEnd = 0;                   // s
    int hypo = 0;                                  // element where V first exceeded vEvent
    double sHypo = 0;                              // m along the fault
    double slipMax = 0, slipMean = 0;              // m (over the ruptured elements)
    double length = 0;                             // m (elements with > 1 % of the max slip)
    double moment = 0, Mw = 0;                     // N m, -
    double vMax = 0;                               // m/s
    double stressDrop = 0;                         // Pa (mean over the ruptured elements)
};

struct Output {
    bool ok = false; std::string error;
    std::vector<double> s;                         // element centres (m from (x0,y0))
    std::vector<double> a;                         // per-element a
    std::vector<double> t, vmax, pMax;             // per accepted step (decimated to <= 20000): time, max V, max dp
    std::vector<double> snapT;                     // snapshot times
    std::vector<std::vector<double>> snapLogV;     // per snapshot log10 V per element
    std::vector<std::vector<double>> snapSlip;     // per snapshot slip per element (m)
    std::vector<Event> events;
    std::vector<double> finalSlip, finalV;
    double hStar = 0, Lb = 0, dx = 0;              // nucleation length (Rubin-Ampuero) and process-zone scale Lb = G' Dc / (b sigma'), element size
    double kc = 0, kSelf = 0;                      // critical stiffness (b - a) sigma' / Dc and the self stiffness of one element
    long steps = 0, rejected = 0;
    double seconds = 0;
};

bool run(const Input& in, Output& out, std::atomic<bool>* cancel = nullptr);

// Initial stresses + geometry from an EFEM contact solve: the fault is embedded fracture 'frac' of 'efin'; sigma_n and |tau|
// per element are interpolated from the contact tractions of the cut elements along the fracture (tn, tt). tau0 is capped
// at the steady-state strength at vInit times capFrac (a fault that the static solve shows slipping starts close to failure).
bool fromEfem(const efem::Input& efin, const efem::Result& efres, int frac, Input& in, std::string& err, double capFrac = 0.98);
// Pore-pressure table from a G1-c coupled flow run (piece pressures per report step, relative to the first step), mapped
// onto the rsf elements by the signed distance along fracture 'frac' (pieces of all layers averaged).
bool fromCoupled(const geocouple::Output& co, int frac, Input& in, std::string& err);

// helpers
double nucleationLength(double G, double nu, double a, double b, double Dc, double sigmaEff);   // h* = 2 G' Dc b / (pi sigma' (b-a)^2)

}  // namespace rsf
}  // namespace ressim
