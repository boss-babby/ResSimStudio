// Water / steam properties for the thermal simulator: IAPWS-IF97 regions 1 (liquid), 2 + 5 (vapour, to 2273 K) and 4
// (saturation line) and the IAPWS 2008 viscosity correlation (without the critical enhancement). SI units throughout: Pa, K,
// m3/kg, J/kg, Pa s. Outside their nominal ranges the region equations are smooth extrapolations (metastable liquid /
// vapour), which the thermal flash uses for derivatives near the saturation line.
#pragma once

namespace ressim {
namespace steam {

constexpr double kMw = 0.018015268;   // kg/mol
constexpr double kR = 461.526;        // J/(kg K)
constexpr double kTc = 647.096, kPc = 22.064e6, kRhoc = 322.0;

double psat(double T);                // saturation pressure (Pa); clamped to [273.15 K, Tc]
double tsat(double p);                // saturation temperature (K); clamped to [611.213 Pa, Pc]
double dpsatdT(double T);             // Pa/K

struct Props {
    double v = 0;                     // specific volume, m3/kg
    double h = 0;                     // specific enthalpy, J/kg
    double u = 0;                     // specific internal energy, J/kg
    double cp = 0;                    // J/(kg K)
};
Props liquid(double T, double p);     // region 1
Props vapour(double T, double p);     // region 2, region 5 above 1073 K (p -> 0: ideal gas)
double vapourZ(double T, double p);   // compressibility factor p v / (R T) of region 2 (finite at p = 0)
double vapourMixCap(double T);    // enthalpy-evaluation cap (Pa) for water in gas mixtures: p on the 0.75 rho_c isochore
double viscosity(double T, double rho);   // IAPWS 2008, Pa s

} // namespace steam
} // namespace ressim
