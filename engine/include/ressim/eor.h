// CO2 / gas EOR laboratory tools on the EOS fluid model:
//  - swelling test and P-X (pseudo-binary) diagram: saturation pressure, swelling factor and first-contact
//    miscibility pressure (FCMP = cricondenbar of the oil / injection-gas pseudo-binary)
//  - multi-contact minimum miscibility pressure by the mixing-cell method (Ahmadi & Johns, SPE J. 2011):
//    combined vaporising / condensing drive, key tie lines tracked to zero length
//  - MMP correlations for cross-checks (Yellig-Metcalfe 1980, Cronquist 1978)
//  - WAG schedule generator (alternating gas-stream / water injection cycles on an injector)
#pragma once
#include "ressim/eos.h"
#include <string>
#include <vector>

namespace ressim {
struct Case;
namespace eor {

// ---- swelling test / P-X diagram
struct SwellRow {
    double alpha = 0;        // mole fraction of injection gas in the mixture
    double gasOilRatio = 0;  // mol gas added per mol of original oil = alpha / (1 - alpha)
    double gasAddedStd = 0;  // standard volume of gas added per saturated volume of original oil (sm3/m3 = scf/scf)
    double psat = -1;        // saturation pressure (Pa); < 0: none found (supercritical / single phase)
    bool dew = false;        // true: dew point (vapour-like mixture), false: bubble point
    double swelling = 0;     // V_sat(mixture) / V_sat(original oil) per mole of original oil
    double rhoSat = 0;       // density of the saturated mixture (kg/m3)
};
struct SwellingResult {
    bool ok = false;
    std::string message;
    double psatOil = -1;             // bubble point of the original oil (Pa)
    double fcmp = -1, fcmpAlpha = 0; // first-contact miscibility pressure (maximum saturation pressure of the P-X curve)
    std::vector<SwellRow> rows;
};
// alphas: mole fractions of gas in the mixture (0 <= alpha < 1). stdT / stdP: standard conditions for the gas volume.
SwellingResult swellingTest(const eos::Model& m, double T, const std::vector<double>& oil, const std::vector<double>& gas,
                            const std::vector<double>& alphas, double stdT = 288.71, double stdP = 101325.0);

// ---- mixing-cell MMP (Ahmadi & Johns 2011)
struct MmpOptions {
    double pLo = 5e6, pHi = 5e7;     // pressure search interval (Pa)
    int nP = 16;                     // pressures of the scan
    int contacts = 60;               // mixing contacts per pressure
    double ratio = 0.5;              // mixing ratio (mole fraction of the upstream vapour in each new cell)
};
struct MmpResult {
    bool ok = false;
    std::string message;
    std::vector<double> p, tlMin;    // scan: pressure (Pa) and minimum key tie-line length after the last contact
    double mmp = -1;                 // MMP: zero intercept of the minimum tie-line length (Pa)
    double mmpBracket = -1;          // first scanned pressure at which the minimum tie line vanished (Pa)
    std::vector<double> profileTL;   // tie-line length of every cell at the highest immiscible pressure (diagnostic)
    double profileP = -1;
    int flashes = 0;
};
// minimum tie-line length of the mixing-cell profile at pressure P (0 = a cell crossed the critical region);
// tls: optional tie-line length of each cell after the last contact (single-phase cells: 0)
double mixingCellMinTieLine(const eos::Model& m, const eos::TParams& tp, double P, const std::vector<double>& oil,
                            const std::vector<double>& gas, int contacts, double ratio, std::vector<double>* tls = nullptr,
                            int* flashes = nullptr);
MmpResult mixingCellMmp(const eos::Model& m, double T, const std::vector<double>& oil, const std::vector<double>& gas,
                        const MmpOptions& o = {});

// ---- slim-tube MMP: 1-D compositional displacement (ResSim EOS core) at each pressure; recovery of the stock-tank oil
// in place after 'pvInjected' hydrocarbon pore volumes of gas. MMP = break-over of the recovery-pressure curve
// (intersection of the immiscible and the miscible straight lines) and the pressure where recovery reaches 'rfMiscible'.
struct SlimTubeOptions {
    int cells = 100;
    double length = 12.0;           // m (40 ft)
    double pvInjected = 1.2;        // hydrocarbon pore volumes
    double days = 1.0;              // injection time
    double rfMiscible = 0.90;
    double corey = 2.0;             // Corey exponent of oil and gas kr (no residual saturations)
};
struct SlimTubeResult {
    bool ok = false;
    std::string message;
    std::vector<double> p, rf;      // pressure (Pa), oil recovery fraction at pvInjected
    double mmpBreak = -1, mmpRf = -1;   // Pa
    std::vector<std::vector<double>> profiles;   // per pressure: recovery vs injected PV (10 points)
};
SlimTubeResult slimTubeMmp(const eos::Model& m, double T, const std::vector<double>& oil, const std::vector<double>& gas,
                           const std::vector<double>& pressures, const SlimTubeOptions& o = {});

// ---- correlations (pure CO2 injection), T in K, result in Pa
double mmpYelligMetcalfe(double T);
// Cronquist: mwC5plus = molecular weight of the C5+ fraction (g/mol), yLight = mole fraction of C1 + N2 in the oil
double mmpCronquist(double T, double mwC5plus, double yLight);

// ---- WAG schedule generator
struct WagSpec {
    std::string well;                // injector (any injector type; the phase is set per half-cycle)
    std::string stream = "CO2";      // WELLSTRE stream injected in the gas half-cycles (empty: injector's default gas)
    double start = 0;                // days from the start of the run
    double gasDays = 90, waterDays = 90;
    int cycles = 6;
    bool gasFirst = true;
    double gasRate = 0, waterRate = 0;   // surface rates (case units)
    double bhpMax = 0;               // injection BHP limit (case units)
    double taper = 1.0;              // gas half-cycle length multiplier per cycle (tapered WAG; 1 = constant)
    double chaseWaterDays = 0;       // final water chase (0 = none)
    bool shutAfter = false;          // shut the injector after the last half-cycle / chase
};
// Splits report steps at every half-cycle boundary (appending steps if the schedule is shorter), removes the
// well's existing controls / WINJGAS from 'start' on and writes the alternating injection controls.
// Returns the number of half-cycles written; message receives a summary or the error.
int applyWagSchedule(Case& c, const WagSpec& w, std::string* message = nullptr);

}  // namespace eor
}  // namespace ressim
