// Gas / solute property library for the unified multi-component framework (C10): critical constants for the
// Peng-Robinson gas phase, aqueous solubility (IAPWS G7-04 Henry's constants, partial molar volumes at infinite
// dilution, Sechenov salting-out coefficients), ideal-gas heat capacities, element formulas and element-balance checks
// for user-defined reactions.
//
// Sources: critical constants - Poling, Prausnitz & O'Connell (5th ed.) / NIST WebBook; Henry's constants - IAPWS G7-04
// (Fernandez-Prini, Alvarez & Harvey, J. Phys. Chem. Ref. Data 32 (2003) 903), valid from the triple point to near the
// critical point of water; partial molar volumes - Brelvi & O'Connell type values at 25 degC (CO2: Garcia 2001
// temperature fit); Sechenov coefficients (NaCl, natural-log basis) - Weisenberger & Schumpe (1996) at 25 degC;
// heat capacities - linear fits of NIST JANAF data between 300 and 1000 K.
#pragma once
#include <map>
#include <string>
#include <vector>

namespace ressim {
namespace gasdata {

struct GasInfo {
    const char* name;                   // canonical name (CNAMES convention)
    const char* formula;                // element formula
    double mw;                          // kg/mol
    double tc, pc, omega;               // K, Pa, - (Peng-Robinson)
    double hA, hB, hC;                  // IAPWS G7-04 ln(kH / p*) = A / Tr + B tau^0.355 / Tr + C Tr^-0.41 exp(tau); 0 = none
    double vInf;                        // partial molar volume at infinite dilution in water, m3/mol (< 0: Garcia CO2 fit)
    double ksalt;                       // Sechenov coefficient, ln basis, kg/mol: ln(x_pure / x_brine) = ksalt m_NaCl
    double cp300, cpSlope;              // ideal-gas heat capacity cp = cp300 + cpSlope (T - 300 K), J/(mol K), J/(mol K2)
};

// library entry by component name (case-insensitive; aliases C1 = CH4, C2 = C2H6, STEAM / WATER = H2O, ...); nullptr if
// unknown
const GasInfo* find(const std::string& name);
const std::vector<GasInfo>& all();

// IAPWS G7-04 Henry's constant k_H = lim f / x (Pa) of a solute on the water saturation curve at T (K); T is clamped to
// [273.16, 0.999 Tc] (the formulation is singular at the critical point)
double henryIapws(double A, double B, double C, double T);
// d ln k_H / dT (1/K), analytic
double henryIapwsDlnDT(double A, double B, double C, double T);
// partial molar volume of CO2 at infinite dilution (Garcia 2001), m3/mol, T in K
double vInfCO2(double T);
// water vapour pressure from the Wagner-Pruss auxiliary equation (Pa) used by G7-04 (consistent with IF97 to < 0.02 %)
double psatWagnerPruss(double T);

// ---- element formulas
// parses "CaCO3", "C2H6", "Ca(OH)2", "CH1.6" (fractional counts allowed); false on syntax errors
bool parseFormula(const std::string& formula, std::map<std::string, double>& elements);
// mol CO2 stored per mol of a carbonate mineral (number of CO3 groups of the formula; 0 if none)
double carbonateCO2(const std::string& formula);
// element balance of a reaction: sum_i (prod_i - sto_i) elements_i for every element (components without a formula are
// skipped and listed in 'unknown'); returns the largest |imbalance| (mol element per mol of reaction extent)
double elementImbalance(const std::vector<std::string>& formulas, const std::vector<double>& sto, const std::vector<double>& prod,
                        std::map<std::string, double>* perElement = nullptr, std::vector<int>* unknown = nullptr);

} // namespace gasdata
} // namespace ressim
