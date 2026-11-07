// C23 geomechanics G4: phase-field model of a fluid-driven fracture (2-D plane strain; lab / near-wellbore scale).
//
// Variational fracture (Francfort-Marigo, Bourdin et al. 2000) with the AT2 regularisation and fluid pressure work
// (Bourdin, Chukwudozie & Yoshioka 2012, SPE 159154; Wheeler, Wick & Wollner 2014):
//   Pi(u, v) = int (v^2 + k) psi(eps(u)) dx  -  p int u . grad v dx  +  Gc int [ (1 - v)^2 / (4 l) + l |grad v|^2 ] dx
// v = 1 intact, v = 0 broken; the crack volume is V = int u . grad v dx. Irreversibility by the history field
// H = max_t psi+(eps) (Miehe et al. 2010) with the volumetric-deviatoric split psi+ = K/2 <tr eps>+^2 + mu |dev eps|^2
// (Amor et al. 2009) in the hybrid formulation (Ambati et al. 2015): isotropic degradation in the momentum equation, split
// history in the phase-field equation - both sub-problems are linear and are solved alternately (staggered) per load step:
//   momentum   int (v^2 + k) eps(du) : C : eps(u) = int p grad v . du
//   phase      int (2 H + Gc/(2 l)) v dv + 2 l Gc grad v . grad dv = int Gc/(2 l) dv + p u . grad dv,  then v <- min(v, v_prev)
// Rectilinear mesh, uniform (size h) in a core box around the crack line and geometrically graded outside; outer boundary
// fixed (an infinite plate for a domain >> crack). Momentum: Q4 + Wilson incompatible modes (condensed per element; plain Q4
// shear-locks on the elongated elements where the graded rows meet the core and halves the opening). The initial crack is
// three node rows (|y| <= h) with v = 0, so the elements next to the crack line are fully broken. Sparse systems by BiCGSTAB
// + smoothed-aggregation AMG (ILU0 fallback). Diagnostics: env PF_DEBUG=1 prints linear-solver statistics, PF_FROZEN=1
// freezes v at the analytic AT2 profile 1 - exp(-d / 2l) (elasticity-only check).
// Verification (tests/c23_tests.inc, testPhaseFieldSneddonGriffith): like every diffuse model the crack-opening volume and the
// critical pressure carry an O(l / a) error; with l = 2h: Sneddon ratio V / (2 pi p a^2 / E') 1.80 -> 1.38 and
// p_c / p_Griffith 0.875 -> 0.925 (load-step midpoints) for l/a = 0.2 -> 0.1; linear extrapolation l -> 0 gives 0.96 and 0.975 (Griffith
// p_c = sqrt(E' Gc / (pi a)), Sneddon 1950). Use l/a <= 0.1 for quantitative work. Experimental (2-D, horizontal crack, pressure
// control: growth after p_c is unstable and is not followed - the run stops when the crack leaves the core box).
// C26 (C24-b): any crack orientation (angleDeg; the uniform core box is rotated to contain it), anisotropic in-situ stress
// (sxx, syy, sxy, compression positive: the crack faces carry the traction (p I - S) n, i.e. the scalar pressure of the
// Bourdin terms becomes the tensor p I - S in the momentum load int grad v . T . du and in the phase-field load u . T . grad dv,
// which gives the net pressure p - sigma_n and the shear traction tau on the faces and leaves the far field unloaded - the
// displacement is the perturbation of the in-situ state), and VOLUME control (control = 1): the momentum problem is linear
// in p for a fixed phase field, u = p u_1 + u_S, so every staggered iteration sets p = (V_target - V(u_S)) / V(u_1) exactly
// - growth after the peak pressure is stable and followed (toughness regime: p - sigma_n -> K_IC / sqrt(pi a) as a grows).
// Griffith reference with stress: p_c = sigma_n + sqrt(E' Gc / (pi a) - tau^2) (mixed-mode energy release rate).
#pragma once
#include <atomic>
#include <functional>
#include <string>
#include <vector>

namespace ressim {
namespace phasefield {

struct Step;
// C30 field output of a completed load step for the viewers (local frame: crack through the origin along +x for angle 0)
struct Field {
    std::vector<double> x, y;            // rectilinear node coordinates (m)
    std::vector<float> v;                // phase field per node i + nx1 j (1 intact, 0 broken)
    std::vector<float> ux, uy;           // displacement per node (m; perturbation of the in-situ state)
    std::vector<float> sxx, syy, sxy;    // stress change per element i + nx j at its centre (Pa, tension +), degraded (v^2 + k)
    std::vector<float> codX, cod;        // crack-opening displacement along the crack line: element-column centres x (m) and
                                         // opening w(x) = int u . grad v dy (m); empty for an inclined crack
    double p = 0, volume = 0, half = 0, angleDeg = 0;   // pressure of the step (Pa), crack volume (m2), half-length (m)
};
struct Input {
    double crackA = 0.1;            // initial half-length a (m), crack along x through the domain centre
    double domain = 20;             // domain half-width / a (square domain of side 2 domain a)
    double h = 0.005;               // core element size (m)
    double ell = 0;                 // regularisation length l (m); 0 = 2 h
    double coreX = 3.0, coreY = 0.5;// core box half-sizes / a (uniform h inside)
    double grade = 1.12;            // geometric growth of the element size outside the core
    double E = 10e9, nu = 0.25;     // Pa, -
    double Gc = 100;                // fracture toughness (J/m2); K_IC = sqrt(Gc E')
    double kres = 1e-6;             // residual stiffness
    // loading: pressure ramp 0 .. pMax in 'steps' increments (pMax <= 0: 1.3 x the Griffith critical pressure)
    double pMax = 0; int steps = 26;
    int maxStagger = 300; double tol = 1e-4;   // staggered iteration: max |dv| per step
    double stopLength = 0.9;        // stop when the crack half-length reaches this fraction of the core box
    // C26 (C24-b)
    double angleDeg = 0;            // initial crack direction from the x axis (degrees)
    double sxx = 0, syy = 0, sxy = 0;   // in-situ total stress (Pa, compression positive)
    int control = 0;                // 0 pressure ramp (pMax), 1 injected-volume ramp (vMax): stable growth followed
    double vMax = 0;                // injected volume per m of thickness (m2); <= 0: 6 x the Sneddon volume at the Griffith pressure
    // C28 in-run coupling (Simulator, HFGROW ... PHASEFIELD): with volume control, the target of load step st comes from
    // nextVolume (return false = stop) instead of the ramp and the number of steps is open-ended; onStep reports every
    // completed load step (state, direction of the + tip). The run keeps its state between the calls (one worker thread).
    std::function<bool(int, double&)> nextVolume;
    std::function<void(const Step&, double)> onStep;
    std::function<void(const Field&)> onField;   // C30: field of every completed load step (optional, O(nodes) per step)
};

struct Step { double p = 0; double volume = 0, lengthHalf = 0; int stagger = 0; bool converged = true; };

struct Output {
    bool ok = false; std::string error;
    std::vector<double> x, y;       // node coordinates (rectilinear)
    std::vector<double> v;          // final phase field per node (i + nx1 j)
    std::vector<double> vAtPc;      // phase field at the first step with growth
    std::vector<Step> steps;
    double ell = 0, h = 0, GcEff = 0, Eprime = 0;
    double pcTheory = 0;            // Griffith sqrt(E' Gc / (pi a)), physical Gc (GcEff = Gc (1 + 3h/8l) for reference only)
    double pcNumeric = 0;           // midpoint of the load step in which the crack jumps (unstable growth), Pa
    double pcBracketLo = 0, pcBracketHi = 0;   // last pressure without / first pressure with the jump
    double sneddonRatio = 0;        // V / (2 pi p a^2 / E') at the first (linear-regime) step; -> 1 + O(l/a) (see ENGINES 6.16)
    double sigmaN = 0, tauS = 0;    // C26: normal (compression +) and shear stress of the in-situ state on the initial crack, Pa
    double tipAngleDeg = 0;         // C26: direction of the furthest damaged point (v < 0.5) of the + tip from the centre, degrees
    int controlUsed = 0;
    int nodes = 0; double seconds = 0;
};

bool run(const Input& in, Output& out, std::atomic<bool>* cancel = nullptr,
         const std::function<void(const std::string&)>& log = {});

}  // namespace phasefield
}  // namespace ressim
