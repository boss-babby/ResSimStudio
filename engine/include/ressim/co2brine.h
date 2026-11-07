// CO2 - brine thermodynamics for CO2 storage (CO2STORE, milestone C4).
//
// Components: CO2 and H2O; NaCl is a fixed salinity (molality m, mol / kg H2O) carried by the aqueous phase.
// Phases: aqueous (brine with dissolved CO2) and a CO2-rich phase (gas / supercritical, with vaporised water).
//
//   Mutual solubility   Spycher, Pruess & Ennis-King (2003), GCA 67, 3015 (12-100 C, up to 600 bar): Redlich-Kwong
//                       CO2 with a(T) = 7.54e7 - 4.13e4 T, b = 27.80 (CO2), 18.18 (H2O) cm3/mol, a(H2O-CO2) = 7.89e7,
//                       true equilibrium constants K(H2O), K(CO2,g) with Poynting corrections (V = 18.1 / 32.6 cm3/mol).
//                       Salinity: Spycher & Pruess (2005), GCA 69, 3309: y_H2O = A (1 - x_CO2 - x_salt), x_CO2 = B' (1 - y_H2O),
//                       B' = B / gamma'_CO2 with the Duan & Sun (2003) activity coefficient ln gamma = 2 lambda m + zeta m^2
//                       (converted to the mole-fraction scale), solved by a short fixed-point loop on the salt mole fraction.
//   CO2-rich phase      molar volume from the same Redlich-Kwong equation (fitted to pure CO2 PVT data by Spycher et al.),
//                       viscosity Fenghour, Wakeham & Vesovic (1998).
//   Aqueous phase       brine density / viscosity Batzle & Wang (1992); dissolved CO2 by the apparent molar volume of
//                       Garcia (2001):  v_aq = x_w M_bw / rho_brine + x_CO2 V_phi  per mole of (CO2 + H2O).
// All quantities SI (Pa, K, m3, kg, Pa s) unless stated otherwise.
#pragma once
#include <string>
#include <utility>
#include <vector>

namespace ressim {
namespace co2 {

constexpr double MW_CO2 = 44.0095e-3, MW_H2O = 18.01528e-3, MW_NACL = 58.4428e-3;   // kg/mol

struct Fluid {
    double T = 323.15;          // K
    double molality = 0;        // NaCl, mol / kg H2O
    // derived (init)
    double tC = 50, aCO2 = 0, lnK0w = 0, lnK0c = 0, salt = 0, vphi = 0, muBrine = 0;
    double mwBrinePerWater = MW_H2O;   // kg of brine per mol of H2O (water + its salt)
};
Fluid makeFluid(double T, double molality);
double molalityFromPpm(double ppmNaCl);                     // mass ppm of NaCl in brine -> molality
double ppmFromMolality(double m);

// Equilibrium compositions at (T, P): CO2 mole fraction of the aqueous phase on a salt-free basis (CO2 + H2O = 1) and
// H2O mole fraction of the CO2-rich phase. Optional pressure derivatives (central differences of the closed form).
void solubility(const Fluid& f, double P, double& xCO2, double& yH2O);
void solubility(const Fluid& f, double P, double& xCO2, double& yH2O, double& dx_dP, double& dy_dP);
// CO2 solubility in mol / kg H2O (reporting, tests)
double co2Molality(const Fluid& f, double P);

// CO2-rich phase molar volume (m3/mol) and dV/dP (Redlich-Kwong, pure CO2 volumetrics)
double gasMolarVolume(const Fluid& f, double P, double* dVdP = nullptr);
double co2Density(const Fluid& f, double P);               // kg/m3 (pure CO2)
// Fenghour et al. (1998) CO2 viscosity (Pa s) at density rho (kg/m3) with d mu / d rho
double co2Viscosity(double T, double rho, double* dmu_drho = nullptr);

// brine (CO2-free) density, kg/m3, and its pressure derivative (Batzle-Wang)
double brineDensity(const Fluid& f, double P, double* drho_dP = nullptr);
double brineViscosity(const Fluid& f);                      // Pa s (Batzle-Wang, pressure independent)
// aqueous phase molar volume per mole of (CO2 + H2O) with derivatives in P and the salt-free CO2 fraction x
double aqueousMolarVolume(const Fluid& f, double P, double x, double* dP = nullptr, double* dx = nullptr);
double aqueousDensity(const Fluid& f, double P, double x);   // kg/m3 incl. salt and dissolved CO2

// PT flash of the binary (salt-free overall CO2 mole fraction z)
struct FlashResult {
    int phases = 1;             // 1 or 2
    bool gas = false;           // single phase: CO2-rich (true) or aqueous (false)
    double beta = 0;            // CO2-rich phase mole fraction
    double x = 0, y = 1;        // CO2 mole fraction in the aqueous / CO2-rich phase
    double dx_dP = 0, dy_dP = 0;
};
FlashResult flash(const Fluid& f, double P, double z);

// standard-condition volumes per mole: CO2 gas (sm3/mol) and brine (sm3 per mol H2O)
double gasSurfaceVolume(double Tstd, double Pstd);
double brineSurfaceVolume(const Fluid& f, double Tstd, double Pstd);

std::string describe(const Fluid& f);   // one-line summary for logs / UI

// ---------------------------------------------------------------------------------------------------------------------
// CO2 trapping mechanisms and their interactions (milestone C5) - options shared by the case (UI / deck) and the core.
//   structural / stratigraphic  sealing caprock = SATNUM region with a capillary entry pressure and / or low permeability
//                               (inventory: mobile CO2 at rest beneath a seal, see structTol)
//   residual (capillary)        gas relative-permeability hysteresis, Killough or Carlson / Land (SATOPTS HYSTER, EHYSTR,
//                               IMBNUM imbibition tables, or a Land maximum trapped gas saturation sgtMax)
//   solubility                  Spycher-Pruess / Duan-Sun dissolution (co2brine), density-driven convective mixing (the
//                               Garcia volume makes CO2-rich brine heavier) and molecular diffusion (DIFFUSE, DIFFCWAT,
//                               DIFFCGAS): J = - phi S xi D grad x  (OPM / E300 form; optional Millington-Quirk tortuosity)
//   mineral                     kinetic carbonation of reactive silicates (Lasaga transition-state rate law with the
//                               Palandri & Kharaka (2004) neutral-mechanism constants), lumped as
//                               silicate + nu CO2(aq) + nu_w H2O -> carbonate (+ clay / quartz)
//                               r = A k25 exp(-Ea/R (1/T - 1/298.15)) Sw (1 - Q/K),  Q/K = mEq / m_CO2(aq)
//   interactions                salt precipitation on dry-out (PRECSALT, VAPWAT, SALTSOL), porosity loss from salt and
//                               minerals -> permeability (PERMFACT table or Verma-Pruess tubes-in-series) -> injectivity;
//                               trapped gas keeps dissolving; dissolution feeds mineralisation.
struct MineralReaction {
    std::string name = "Anorthite -> calcite + kaolinite";
    double logk25 = -9.12;          // log10 rate constant at 25 C, mol / (m2 s) (Palandri & Kharaka 2004, neutral)
    double Ea = 17.8e3;             // activation energy, J/mol
    double area = 1000.0;           // reactive surface area, m2 per m3 of pore volume
    double nuCO2 = 1.0, nuH2O = 2.0;   // moles consumed per mole of reaction
    double dVsolid = 35.66e-6;      // net solid volume change per mole of reaction (m3/mol; > 0 = porosity loss)
    double reactant0 = 0.05;        // initial reactive silicate, volume fraction of the pore volume (inventory limit)
    double vmReactant = 100.79e-6;  // molar volume of the reactive silicate, m3/mol
    double mEq = 0.01;              // dissolved-CO2 molality at which the overall reaction is at equilibrium (mol/kg)
};
struct TrappingOptions {
    // residual trapping: maximum trapped gas saturation of Land's model (< 0: from the IMBNUM imbibition table)
    double sgtMax = -1;
    // molecular diffusion (m2/s); DIFFUSE switches it on
    bool diffusion = false;
    double dCo2Water = 2.0e-9, dH2oWater = 2.0e-9, dCo2Gas = 1.0e-7, dH2oGas = 1.0e-7;
    int tortuosity = 0;             // 0: phi S (OPM / E300), 1: Millington-Quirk phi^(4/3) S^(10/3)
    // salt precipitation (dry-out)
    bool saltPrecip = false;
    double saltSolubility = 6.15;   // NaCl solubility limit, mol / kg H2O (SALTSOL)
    double rhoSalt = 2165.0;        // kg/m3 (halite)
    std::vector<std::pair<double, double>> permFact;   // PERMFACT rows (phi / phi0, k / k0); empty = Verma-Pruess
    double vpPhiR = 0.8, vpGamma = 0.8;                // Verma-Pruess: fraction of phi0 at which k = 0, tube parameter
    // mineral trapping
    std::vector<MineralReaction> minerals;
    // inventory: mobile CO2 counts as structurally trapped when its gas outflow is below structTol x (cell gas volume) / year
    double structTol = 0.01;
    bool any() const { return diffusion || saltPrecip || !minerals.empty() || sgtMax > 0; }
};
// presets (Palandri & Kharaka 2004 neutral-mechanism rate constants; approximate molar volumes)
std::vector<MineralReaction> mineralPresets();
// permeability factor k / k0 for a porosity factor phi / phi0 (PERMFACT table or Verma-Pruess)
double permFactor(const TrappingOptions& t, double phiRatio);
}  // namespace co2
}  // namespace ressim
