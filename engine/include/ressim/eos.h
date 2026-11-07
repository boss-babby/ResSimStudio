// Cubic equation-of-state thermodynamics for compositional simulation (E300 class).
//
//  * Peng-Robinson (1976, and the 1978 m(omega) correction for omega > 0.49) and Soave-Redlich-Kwong, van der Waals
//    one-fluid mixing with binary interaction coefficients, Peneloux volume shift (dimensionless SSHIFT s_i = c_i / b_i).
//  * Fugacity coefficients with exact analytic derivatives with respect to mole numbers, pressure and temperature, from
//    the reduced residual Helmholtz function F(T, V, n) of Michelsen & Mollerup ("Thermodynamic Models: Fundamentals and
//    Computational Aspects", 2007, ch. 3); residual enthalpy for thermal use.
//  * Phase stability (Michelsen tangent-plane distance, vapour- and liquid-like trial phases, successive substitution
//    followed by the second-order method in alpha = 2 sqrt(W) variables).
//  * Two-phase PT flash: Rachford-Rice with Leibovici-Neoschil/Nielsen-Lia safeguarded Newton, successive substitution
//    with GDEM acceleration, then Newton on the vapour mole numbers with the full Jacobian (quadratic convergence).
//  * Saturation (bubble / dew) pressure, single-phase labelling (Li critical-temperature rule), Lohrenz-Bray-Clark
//    viscosity, Wilson K-values.
// All internal units are SI: Pa, K, mol, m3, kg, Pa s. Conversion helpers for FIELD / METRIC are at the end.
#pragma once
#include <string>
#include <vector>

namespace ressim {
namespace eos {

constexpr double Rgas = 8.31446261815324;   // J / (mol K)

enum class EosType { PR, SRK };

struct Component {
    std::string name;
    double tc = 0, pc = 0;        // K, Pa
    double vc = 0, zc = 0;        // m3/mol (0 = from zc), critical Z
    double omega = 0;             // acentric factor
    double mw = 0;                // kg/mol
    double sshift = 0;            // Peneloux shift, dimensionless (c = s b)
    double parachor = 0;          // for interfacial tension (optional)
    double vcVisc = 0;            // VCRITVIS for LBC (0 = vc)
    double omegaA = -1, omegaB = -1;   // per-component OMEGAA / OMEGAB (< 0 = EOS default)
};

struct Model {
    EosType type = EosType::PR;
    bool pr78 = true;                    // PR: use the 1978 m(omega) for omega > 0.49 (E300 PRCORR)
    std::vector<Component> comp;
    std::vector<double> kij;             // nc x nc symmetric binary interaction coefficients (empty = 0)
    double lbc[5] = {0.1023, 0.023364, 0.058533, -0.040758, 0.0093324};   // LBCCOEF
    int nc() const { return int(comp.size()); }
    double k(int i, int j) const { return kij.empty() ? 0.0 : kij[size_t(i) * comp.size() + size_t(j)]; }
    void setK(int i, int j, double v) {
        if (kij.empty()) kij.assign(comp.size() * comp.size(), 0.0);
        kij[size_t(i) * comp.size() + size_t(j)] = kij[size_t(j) * comp.size() + size_t(i)] = v;
    }
    double delta1() const;               // PR: 1 + sqrt 2, SRK: 1
    double delta2() const;               // PR: 1 - sqrt 2, SRK: 0
};

// temperature-dependent pure and mixture parameters (computed once per temperature)
struct TParams {
    double T = 0;
    std::vector<double> a, dadT, d2adT2, b;     // per component
    std::vector<double> aij, daijdT, d2aijdT2;  // nc x nc
};
TParams atTemperature(const Model& m, double T);

enum class Root { Liquid, Vapour, MinGibbs };

// Fugacity coefficients of a phase with mole numbers n (any scale) at T, P.
struct PhaseEval {
    double Z = 0;                        // compressibility factor (unshifted)
    double vm = 0;                       // molar volume incl. the Peneloux shift, m3/mol
    double rho = 0;                      // mass density, kg/m3
    double mw = 0;                       // kg/mol
    std::vector<double> lnphi;           // incl. the shift term -c_i P / RT
    std::vector<double> dlnphi_dn;       // nc x nc, n_T * d lnphi_i / d n_j at T, P (scaled: independent of the total)
    std::vector<double> dlnphi_dP;       // 1/Pa
    std::vector<double> dlnphi_dT;       // 1/K
    double hres = 0;                     // residual molar enthalpy, J/mol (H - H_ideal gas)
    double dvm_dP = 0, dvm_dT = 0;       // molar volume derivatives
    std::vector<double> dvm_dx;          // d vm / d n_j * n_T (partial molar volume minus vm)
};
struct EvalOptions { bool dn = true, dP = true, dT = false, hres = false; Root root = Root::MinGibbs; };
PhaseEval evalPhase(const Model& m, const TParams& tp, double P, const std::vector<double>& n, const EvalOptions& o = {});
// all real roots of the cubic in Z (ascending), for diagnostics
std::vector<double> cubicRoots(const Model& m, const TParams& tp, double P, const std::vector<double>& x);

std::vector<double> wilsonK(const Model& m, double T, double P);

// ---- stability
struct StabilityResult {
    bool stable = true;
    double tpd = 0;                      // minimum modified tangent-plane distance found
    std::vector<double> K;               // K estimate from the unstable trial phase (y/z or z/x), for the flash
    int iterations = 0;
};
StabilityResult stability(const Model& m, const TParams& tp, double P, const std::vector<double>& z);

// ---- flash
struct FlashResult {
    int phases = 1;                      // 1 or 2
    bool converged = true;
    double beta = 0;                     // vapour mole fraction (1 phase: 0 = liquid-like, 1 = vapour-like)
    std::vector<double> x, y, K;         // liquid / vapour mole fractions, K = y / x
    PhaseEval L, V;                      // phase properties (only the existing phase is filled when single phase)
    int iterations = 0;
    bool vapourLabel = false;            // single phase: labelled vapour
};
struct FlashOptions {
    double tol = 1e-11;                  // on max |ln f_V - ln f_L|
    int maxSS = 60, maxNewton = 30;
    const std::vector<double>* Kinit = nullptr;   // previous K (skips the stability test when it gives two phases)
    bool skipStability = false;
};
FlashResult flash(const Model& m, const TParams& tp, double P, const std::vector<double>& z, const FlashOptions& o = {});
// Rachford-Rice: vapour fraction beta for given K (may lie outside [0,1] when it is a negative flash)
double rachfordRice(const std::vector<double>& z, const std::vector<double>& K, bool clampToPhysical, int* iterations = nullptr);

// Derivatives of a converged two-phase flash (implicit-function theorem on ln f_V = ln f_L, as used by the
// compositional Jacobian). z is treated as mole numbers of a unit feed: d/dz_j at fixed other z_k.
struct FlashSens {
    double dbeta_dP = 0, dbeta_dT = 0;
    std::vector<double> dx_dP, dy_dP, dx_dT, dy_dT;   // nc
    std::vector<double> dbeta_dz;                     // nc
    std::vector<double> dx_dz, dy_dz;                 // nc x nc, [i * nc + j] = d x_i / d z_j
};
// withT = false: no temperature derivatives (isothermal simulation) - reuses the phase derivatives of the flash result
bool flashSensitivities(const Model& m, const TParams& tp, double P, const std::vector<double>& z, const FlashResult& f, FlashSens& s, bool withT = true);

// Single-phase label (Li's mixture critical temperature): true = vapour-like
bool isVapourLike(const Model& m, const std::vector<double>& z, double T);

// Saturation pressure at T: bubble (z is the liquid) or dew (z is the vapour). Returns < 0 if none found in [pLo, pHi].
double saturationPressure(const Model& m, double T, const std::vector<double>& z, bool bubble, double pGuess = 0,
                          double pLo = 1e4, double pHi = 1e8, std::vector<double>* incipient = nullptr);

// Lohrenz-Bray-Clark viscosity of a phase (Pa s) from its composition, temperature and molar density (mol/m3)
double lbcViscosity(const Model& m, const std::vector<double>& x, double T, double molarDensity);
// same, with the analytic derivatives d mu / d x_i (x treated as independent) and d mu / d molarDensity
double lbcViscosity(const Model& m, const std::vector<double>& x, double T, double molarDensity, double* dmu_dx, double* dmu_dxi);

// ---- E300-style compositional fluid description (PROPS / SOLUTION keywords), SI internally
struct ZmfTable { std::vector<double> depth; std::vector<std::vector<double>> z; };   // depth (m), z per row
struct CompFluid {
    bool present = false;            // COMPS / NCOMPS found
    std::vector<Model> eosRegions;   // one model per EOS region (NTEOS); [0] = reservoir EOS
    std::vector<std::string> names;
    double resT = 0;                 // RTEMP / RTEMPA (K); 0 = not given
    std::vector<std::pair<double, double>> tempvd;   // TEMPVD (depth m, T K) of table 1
    std::vector<ZmfTable> zmfvd;     // per equilibration region
    std::vector<double> zi;          // ZI of region 1 (empty = none)
    double stdT = 288.71, stdP = 101325.0;            // STCOND (K, Pa)
    bool co2store = false, water = false;
    bool co2sol = false;             // CO2SOL (E300): the CO2 component also dissolves in the aqueous phase
    double salinity = 0;             // CO2STORE / CO2SOL: SALINITY (mol NaCl / kg H2O)
    // CO2SOL: SOLUBILI table (E300; table 1), SI: p (Pa), Rsw (sm3 CO2 / sm3 water), Bw,sat (rm3/sm3), cw,sat (1/Pa),
    // mu_w,sat (Pa s). Empty: default solubility (ResSim: Spycher-Pruess / Duan-Sun at SALINITY).
    struct SolRow { double p = 0, rsw = 0, bw = 1, cw = 0, muw = 5e-4; };
    std::vector<SolRow> solubili;
    bool pr78 = false;               // PRCORR
};
}  // namespace eos
namespace deck { struct RawDeck; }
namespace eos {
// Builds the fluid from the deck keywords (NCOMPS, EOS, PRCORR, CNAMES, TCRIT, PCRIT, VCRIT, ZCRIT, ACF, MW, BIC,
// SSHIFT, PARACHOR, VCRITVIS, ZCRITVIS, LBCCOEF, OMEGAA, OMEGAB, RTEMP, RTEMPA, TEMPVD, ZMFVD, ZI, STCOND, CO2STORE).
// fieldUnits: deck unit system (FIELD) - otherwise METRIC. Problems are appended to 'messages'. Returns present.
bool fluidFromDeck(const deck::RawDeck& d, bool fieldUnits, CompFluid& out, std::vector<std::string>* messages = nullptr);
// Library critical properties of common pure components (Poling-Prausnitz-O'Connell / NIST values)
Component standardComponent(const std::string& name);   // C1, C2, C3, IC4, NC4, IC5, NC5, C6, N2, CO2, H2S, H2O, H2

// ---- units
inline double psiaToPa(double p) { return p * 6894.757293168; }
inline double paToPsia(double p) { return p / 6894.757293168; }
inline double barToPa(double p) { return p * 1e5; }
inline double rankineToK(double t) { return t * 5.0 / 9.0; }
inline double fahrenheitToK(double t) { return (t + 459.67) * 5.0 / 9.0; }
inline double celsiusToK(double t) { return t + 273.15; }

}  // namespace eos
}  // namespace ressim
