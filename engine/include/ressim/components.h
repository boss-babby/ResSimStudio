// Multi-component product analysis (compositional and thermal engines): component metadata, standard heating values,
// produced-stream composition and product identification (H2, syngas, LHV / HHV, CO2 / O2 slip).
#pragma once
#include <string>
#include <vector>

namespace ressim {

struct SummaryPoint;

// kind: 0 hydrocarbon / K-value component, 1 non-condensable gas, 2 solid (never produced), 3 water,
//       4 passive tracer of the black-oil engine (C27: 'mol' = tracer amount unit, mw = 1, species = carrier phase 0 water / 1 oil / 2 gas)
struct ComponentInfo {
    std::string name;
    double mw = 0;                      // kg/mol
    int kind = 0;
    int species = -1;                   // recognised species (Species enum) or -1
    double hhv = 0, lhv = 0;            // standard heating values at 25 degC, J/mol (0 for inerts / unknown)
};

// species recognised by name (CNAMES conventions: H2, CO, CO2, CH4 / C1, C2 / C2H6, C3, IC4 / NC4 / C4, IC5 / NC5 / C5,
// C6, N2, O2, H2S, H2O / WATER / STEAM, NH3, AR)
enum Species { SP_H2, SP_CO, SP_CO2, SP_CH4, SP_C2H6, SP_C3H8, SP_C4H10, SP_C5H12, SP_C6H14, SP_N2, SP_O2, SP_H2S, SP_H2O,
               SP_NH3, SP_AR, SP_COUNT };
const char* speciesFormula(int sp);
int recogniseSpecies(const std::string& componentName);
// fills species / hhv / lhv from the names (heavy pseudo components with mw >= 0.072 kg/mol get the Dulong-style
// paraffin correlation LHV = 44.0 MJ/kg x mw, HHV = 47.3 MJ/kg x mw; solids and water 0)
void assignHeatingValues(std::vector<ComponentInfo>& comps);

// analysis of one produced stream (a well or the field) from component molar rates (mol/s, production > 0)
struct StreamAnalysis {
    double total = 0;                   // mol/s of all produced components (incl. water)
    double dry = 0;                     // mol/s without water
    double gas = 0;                     // mol/s of the non-condensable gas components + light species (H2..C4, N2, O2, ...)
    double mass = 0;                    // kg/s (all)
    std::vector<double> yDry;           // mole fraction per component, dry basis (water = 0)
    std::vector<double> yGas;           // mole fraction per component in the produced gas (gas-kind / light species only)
    double h2 = 0, co = 0, co2 = 0, ch4 = 0, o2 = 0, n2 = 0, h2s = 0;   // mol/s of recognised species
    double syngas = 0;                  // CO + H2 (mol/s)
    double h2co = 0;                    // H2 / CO molar ratio (0 when no CO)
    double lhvRate = 0, hhvRate = 0;    // W (heating power of the produced stream, dry)
    double lhvGas = 0, hhvGas = 0;      // MJ/sm3 of produced gas (ideal gas at 15 degC, 1 atm: 42.29 mol/sm3)
    double h2Mass = 0;                  // kg/s of hydrogen
    double co2Slip = 0, o2Slip = 0;     // fraction of the injected CO2 / O2 that is produced (0 when nothing injected)
};
// rates: production molar rates; inj (optional): injection molar rates (positive) for the slip fractions
StreamAnalysis analyseStream(const std::vector<ComponentInfo>& comps, const std::vector<double>& rates,
                             const std::vector<double>* inj = nullptr);

// accumulates component cumulative totals (fcompPT / fcompIT, well compCumP / compCumI) along a summary in place
// (rates must be filled; cumulatives integrate rate x dt of each point, piecewise constant like ECLIPSE totals)
void accumulateComponentTotals(std::vector<SummaryPoint>& summary);
// incremental form for streaming runners: p's cumulatives = prev's + rate x p.dt (prev == nullptr: start from zero)
void accumulateComponentPoint(const SummaryPoint* prev, SummaryPoint& p);

constexpr double kMolPerSm3 = 101325.0 / (8.314462618 * 288.15);   // ideal gas at 15 degC, 1 atm: 42.29 mol/sm3

} // namespace ressim
