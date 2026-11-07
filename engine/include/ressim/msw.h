// C32: multisegment wells (ECLIPSE / OPM WELSEGS, COMPSEGS, WSEGVALV, WSEGSICD, WSEGAICD) in the ResSim black-oil engine.
//
// The deck text stays verbatim in Case::extraKeywords (exact round trip, OPM Flow runs the same deck); the simulator parses
// the definitions in force at every report step.
//
// Model (per well, per time step, from the state at the start of the step - the same explicit treatment ResSim uses for the
// standard-well hydrostatic heads, with the device drops balanced locally):
//   * segment tree: segment 1 = top node at the BHP reference; every other segment has an outlet segment
//     (branches / laterals by branch number); node depth and length along the tubing from WELSEGS (INC or ABS)
//   * connection -> segment from COMPSEGS (the segment of its branch whose node length is nearest to the centre of the
//     perforated interval - ECLIPSE rule); connections not listed: the nearest node by depth (warning)
//   * segment flow = reservoir inflow of all connections in the segment's subtree; mixture density / viscosity weighted
//     by those inflows (no flow: the connection cells' mobility-weighted values)
//   * pressure: p(seg) = p(outlet) + rho_mix g (z_seg - z_out) + sgn(q) [friction + device]
//       friction: Darcy-Weisbach, f = 64 / Re (Re < 2300), Haaland (1983) otherwise, over |L_seg - L_out| (+ WSEGVALV
//       additional length); WELSEGS item 6 'H--' switches friction off
//       devices on a segment:  WSEGVALV  dp = rho v^2 / (2 Cv^2),  v = q / A_c   (SI)
//                              WSEGSICD  dp = K (rho / rho_cal)^(3/4) (mu / mu_cal)^(1/4) q^2   (OPM Flow form)
//                              WSEGAICD  dp = a (rho / rho_cal) (mu_cal / mu)^y q^x   (mixture exponents 1)
//       a device drop is solved with the subtree inflow (q = A (D0 - dp), A = productivity of the subtree, D0 = drawdown
//       without the device): dp - c (A (D0 - dp))^x = 0 on [0, D0] - stable for strong devices (no explicit feedback)
//   * connection pressure = p(seg) + rho_seg g (z_conn - z_seg); friction + device parts under-relaxed 1/2 with the
//     previous step (restart state)
//   * a SHUT valve (WSEGVALV / WSEGSICD / WSEGAICD status) closes every connection of its subtree
//   * C36: WSEGLABY labyrinth device  dp = rho v^2 / 2 [g_in + g_out + (n_ch - 1) g_lab + f(Re, R / D) L_ch / D],
//       v = q / (N_config A_ch), L_ch = L1 + (n_ch - 1) L2 (loss-coefficient reconstruction from the keyword items -
//       the ECLIPSE formula is not public, OPM Flow parses but does not simulate WSEGLABY)
//   * C36: acceleration (WELSEGS item 6 'HFA'): p(seg) - p(out) += (M_seg - sum M_child) / A_seg, momentum flux
//       M = rho q^2 / A (homogeneous) or (rho_g j_g^2 / alpha + rho_l j_l^2 / (1 - alpha)) A (drift flux); inflow
//       through the perforations carries no axial momentum
//   * C36: drift-flux slip (WELSEGS item 7 'DF', WSEGDFPA parameters): Shi et al. (2005) gas-liquid model -
//       v_g = C0 j + v_d, C0 = A / (1 + (A - 1) gamma^2), gamma = (beta - B) / (1 - B), beta = max(alpha, Fv alpha |j| / v_sgf),
//       v_d = (1 - alpha C0) C0 K(alpha) V_c m(theta) / (alpha C0 sqrt(rho_g / rho_l) + 1 - alpha C0),
//       V_c = (sigma g drho / rho_l^2)^(1/4), K = 1.53 / C0 (alpha < a1) .. Ku (alpha > a2), Ku from the Wallis-Makkenchery
//       bond-number fit, m(theta) = sqrt(cos theta) (1 + sin theta)^1.2 (WSEGDFMD default 'H-K'); the in-situ holdup sets
//       the hydrostatic density of co-current upward / downward flow toward the outlet (counter-current and injection:
//       homogeneous); oil and water travel as one liquid (no oil-water slip)
//   * C37 tubing flash: with the connections' surface component rates (oil, water, gas incl. dissolved / vaporised) the
//       segment fluid is re-flashed with the black-oil PVT at the segment pressure (Rs / Rv capped by the saturated
//       values, Bo / Bw / Bg, phase densities and viscosities) - the free gas released up the tubing lightens the head and
//       raises the local volumetric rate (friction); 2 fixed-point passes on the segment pressures
//   * C37 oil-water slip (WSEGDFPA items 7-11, item 7 > 0): Shi et al. (2005) oil-water drift flux inside the liquid -
//       v_o = C0' j_l + v_d', C0' = A / (1 + (A - 1) gamma^2), gamma = (alpha_o - B1) / (B2 - B1) in [0, 1],
//       v_d' = 1.53 V_c' (1 - alpha_o)^N m(theta) * factor, V_c' = (sigma_ow g (rho_w - rho_o) / rho_w^2)^(1/4)
//   * C37 segment storage (transient wellbore): every segment holds its fluid volume V (WELSEGS item 10, default A L);
//       the head density / holdup move toward the steady value with the volume swept in the time step,
//       x_new = x_old + w (x_ss - x_old), w = |q| dt / (V + |q| dt) - a shut-in or slow segment keeps its column
//   * C37 WSEGFLIM flow-limiting valve: per limited phase (OIL / WAT / GAS / LIQ) the subtree's potential phase rate Q0
//       is cut smoothly to the limit, conductance factor f = (1 + (Q0 / Qlim)^a)^(-1/a) (a = item A, sharp for large a)
//   * C37 WSEGPULL downhole pump: pressure boost across the segment for producing flow,
//       dp = P * max(0, DP - G * max(0, q_o - q_lim)) (q_o subtree oil surface rate) - included in the relaxed friction part
//   (WSEGFLIM / WSEGPULL / oil-water formulas are reconstructions from the keyword items: the ECLIPSE manual pages were
//   not accessible; OPM Flow parses but does not simulate them)
// Not modelled: the emulsion / critical-water-fraction terms of SICD / AICD, the rate scaling by device length,
// gravity segregation inside a shut-in segment column.
#pragma once
#include <functional>
#include <string>
#include <vector>

namespace ressim {
struct ExtraKeyword;
namespace msw {

struct Segment {
    int id = 1, branch = 1, outlet = 0;
    double length = 0, depth = 0;         // node length along the tubing and node depth (case units, absolute)
    double diam = 0.1, rough = 1e-5;      // case length units
    int device = 0;                       // 0 none, 1 WSEGVALV, 2 WSEGSICD, 3 WSEGAICD, 4 WSEGLABY
    bool shut = false;
    double cv = 0, area = 0;              // VALV (area in case length^2)
    double addLen = 0, vDiam = 0, vRough = 0;   // VALV additional pipe
    double strength = 0, rhoCal = 0, muCal = 0, x = 2, y = 1;   // SICD / AICD (case units; rhoCal case density units, mu cP)
    // C36 WSEGLABY items 4-13 (lengths / areas in case units)
    int labN = 1, labCh = 2;
    double labA = 0, labL1 = 0, labL2 = 0, labD = 0, labR = 0, gIn = 0.35, gOut = 0.5, gLab = 0.5;
    double volume = 0;                    // C37 WELSEGS item 10 (case length^3; 0: area x length)
    // C37 WSEGFLIM: limited phase per limit (0 none, 1 OIL, 2 WAT, 3 GAS, 4 LIQ), limit (case surface rate), exponent a
    int flimPh[2] = {0, 0};
    double flimQ[2] = {0, 0}, flimA[2] = {8, 8};
    // C37 WSEGPULL pump: scale P, boost DP (case pressure), oil-rate limit (case surface rate), gradient (pressure / rate)
    bool pump = false;
    double pumpP = 1, pumpDp = 0, pumpQo = 0, pumpG = 0;
};
struct DriftParams {                      // WSEGDFPA (Shi et al. 2005 gas-liquid model) + surface tension used for V_c
    double vdFactor = 1.0, a1 = 0.2, a2 = 0.4, A = 1.2, B = 0.3, Fv = 1.0, sigma = 0.03;   // sigma N/m
    // C37 oil-water slip (items 7-11): factor (<= 0: no oil-water slip), profile A, B1, B2, exponent N; sigma_ow N/m
    double owFactor = -1.0, owA = 1.2, owB1 = 0.4, owB2 = 0.7, owN = 2.0, owSigma = 0.025;
};
struct Conn { int i = 0, j = 0, k = 0, branch = 1; double d1 = 0, d2 = 0; int seg = -1; };   // zero-based i j k; seg = Segment id
struct Well {
    std::string name;
    double topDepth = 0, topLength = 0;
    bool frictionOn = true;               // WELSEGS item 6: HFA / HF- (true), H-- (false)
    bool accelOn = false;                 // C36: WELSEGS item 6 'HFA'
    bool driftFlux = false;               // C36: WELSEGS item 7 'DF'
    DriftParams df;
    std::vector<Segment> segs;            // in the order of segment ids (index = position)
    std::vector<Conn> conns;
    std::string warn;
    int indexOf(int id) const;            // -1 if missing
};

// definitions in force at report step 'step' (latest WELSEGS / COMPSEGS per well at or before it, then the device
// keywords in deck order); messages appended to 'err'
std::vector<Well> parse(const std::vector<ExtraKeyword>& kws, int step, bool field, std::string& err);

// ---------------------------------------------------------------- per-time-step pressure solve
struct Inflow {                           // per segment (index = position in Well::segs): connections allocated to it
    std::vector<double> q;                // reservoir rate, + = into the well (case units: rb/d | rm3/d)
    std::vector<double> rho, mu;          // inflow mixture density (case units) / viscosity (cP)
    std::vector<double> a;                // productivity: d q / d (p_res - p_well), case rate / case pressure
    std::vector<double> pres;             // productivity-weighted reservoir pressure of the connections
    // C36 drift flux (optional - empty: homogeneous): free-gas part of q (same sign convention) and the phase densities
    std::vector<double> qg, rhoG, rhoL;
    // C37 tubing flash (optional): surface component rates of the inflow, 3 per segment (oil, water, gas incl. dissolved;
    // + = toward the top), and the PVT flash at a pressure (case units) of |surface rates| -> local phase rates / properties
    std::vector<double> qs;
    struct Flash { double qo = 0, qw = 0, qg = 0, rhoO = 0, rhoW = 0, rhoG = 0, muO = 1, muW = 0.5, muG = 0.02; };
    std::function<Flash(double p, const double* qs)> flash;
    bool flashAtSegment = true;           // false: properties at the subtree reservoir pressure (no tubing flash)
    void resize(size_t n) { q.assign(n, 0); rho.assign(n, 0); mu.assign(n, 1); a.assign(n, 0); pres.assign(n, 0); qg.assign(n, 0); rhoG.assign(n, 0); rhoL.assign(n, 0); }
};
struct SolveOpts {                        // C37
    double dt = 0;                        // time step (days) for the segment storage
    const std::vector<double>* sto = nullptr;   // committed storage state: head density (n) then gas holdup (n); null: none
    const std::vector<double>* flimPrev = nullptr;   // WSEGFLIM potential (subtree reservoir rate) used for the factors in force (n)
    const std::vector<double>* xfSweep = nullptr;    // C45 cross-flow reservoir rate carried down through each segment (n; >= 0)
};
struct Units {
    double g = 0;                         // case pressure per (case density * case length)
    double toM = 1, toRho = 1, toQ = 1.0 / 86400.0, toPa = 1e5;   // case length / density / rate / pressure -> SI
};
struct Solution {
    std::vector<double> p;                // node pressure (case units)
    std::vector<double> qSeg, rhoSeg, muSeg;   // segment flow (+ toward the outlet), mixture density / viscosity
    std::vector<double> dHydro, dFric, dDev;   // drops across each segment (case pressure; signed: p(seg) - p(out) parts)
    std::vector<double> nonHydro;         // accumulated friction from the top node to each node (devices excluded)
    // C32 device coupling as a connection-conductance factor: the drop dp of a device at the linearised subtree inflow
    // q = A (D0 - dp) is represented by scaling the subtree's connection factors with (D0 - dp) / D0 instead of a lagged
    // pressure shift (a lagged shift of the order of the drawdown blocks / floods connections when the BHP moves);
    // product over the devices on the path to the top; 1 without devices
    std::vector<double> devFac;
    // C36: acceleration drop per segment (signed, included in p and nonHydro); in-situ gas holdup and the no-slip gas
    // fraction (drift flux; holdup = no-slip fraction for homogeneous flow); density used for the hydrostatic head
    std::vector<double> dAcc, holdup, noSlip, rhoHead;
    // C37: oil fraction of the in-situ liquid (oil-water slip; no-slip value otherwise), pump boost (signed part of p),
    // flow-limiter factor, local phase rates after the tubing flash (case reservoir rate) and the free-gas surface rate
    std::vector<double> oilHold, dPump, flimFac, qoLoc, qwLoc, qgLoc;
    std::vector<double> flimQ0;           // C37 WSEGFLIM potential used (subtree reservoir rate; 0 none) - state for the next solve
};
// pTop: pressure at the top node; nonHydroPrev (optional): accumulated non-hydrostatic part of the previous step for the
// 1/2 under-relaxation (empty: none)
Solution solve(const Well& W, double pTop, const Inflow& in, const Units& u, const std::vector<double>& nonHydroPrev = {},
               const SolveOpts& opts = SolveOpts());
double frictionFactor(double Re, double relRough);   // Darcy f
// C36: labyrinth device drop (Pa) at a total flow qSI (m3/s, >= 0) with mixture density / viscosity (SI)
double labyrinthDp(const Segment& g, double qSI, double rhoSI, double muSI, double toM);
// C36: Shi et al. (2005) in-situ gas fraction for superficial velocities jg, jl >= 0 (m/s, co-current along the flow),
// densities kg/m3, diameter m, theta = angle of the flow direction from the upward vertical (0 = up, pi/2 horizontal,
// pi = down); optional outputs C0 and v_d (m/s)
double driftHoldup(double jg, double jl, double rhoG, double rhoL, double D, double theta, const DriftParams& P,
                   double* C0 = nullptr, double* vd = nullptr);
// C37: Shi et al. (2005) oil-water model: in-situ oil fraction of the liquid for superficial jo, jw >= 0 (m/s, co-current)
double owHoldup(double jo, double jw, double rhoO, double rhoW, double theta, const DriftParams& P, double* C0 = nullptr, double* vd = nullptr);
// C37 WSEGFLIM smooth limiter factor for a potential rate q0 >= 0 and limit qlim > 0
double flowLimitFactor(double q0, double qlim, double a);
std::vector<std::vector<int>> children(const Well& W);   // per segment position: child positions
std::vector<int> topoOrder(const Well& W);               // parents before children (positions)

}  // namespace msw
}  // namespace ressim
