// Thermal reservoir simulator core (milestone C7) - STARS-style steam / hot-water processes.
//
// Formulation (fully implicit, mass + energy conservative):
//   primary unknowns per cell   p, N_1 .. N_nc (hydrocarbon / solvent component moles), N_w (water moles), U (internal
//                               energy of fluids + rock in the cell, J)
//   equations per cell          component balances  N_k - N_k^n + dt (sum_faces F_k + q_k) = 0       (k = 1..nc, w)
//                               energy balance      U - U^n + dt (sum_faces [F_E + F_cond] + q_E + Q_loss) = 0
//                               volume balance      PV(p, T) - V_W - V_O - V_G = 0
//   F_k = T_ab sum_alpha [kr/mu xi x_k]_up dPhi_alpha, F_E = T_ab sum_alpha [kr/mu xi h]_up dPhi_alpha,
//   F_cond = (A/d) lambda_eff (T_a - T_b).
// Using the energy U (not T) as primary variable makes the state unique also for univariant systems (pure water /
// steam with a non-volatile oil): the temperature and the steam quality follow from a (p, N, U) flash - an energy
// search T -> E(T) = U over PT flashes; on the saturation line (E jumps by the latent heat) T = Tsat(p) and the steam
// fraction is set by the energy. Local cell properties are differentiated numerically (forward differences of the
// local (p, N, U) state), the face, well and heat-loss terms by the chain rule; the block system (nb = nc + 3) is solved
// with GMRES + BlockCPR.
//
// Phase behaviour (unified with the compositional EOS core in later milestones): phases aqueous (W: water), oil (O:
// hydrocarbon liquid) and gas (G: hydrocarbon vapour + steam). Water / steam equilibrium y_w p = Psat(T) (IAPWS-IF97);
// hydrocarbon K-values K = (kv1 + kv2 / p + kv3 p) exp(kv4 / (T - kv5)) (STARS / E300 KVCR form). Free-water three-phase
// Rachford-Rice: with K'_i = K_i / (1 - y_w) the oil / gas split is a standard Rachford-Rice problem; without an
// aqueous phase water is a gas-only component. Liquid densities rho = rho_ref exp(c_p dp - ct1 dT - ct2 / 2 d(T^2)),
// liquid viscosities mu = a exp(b / T) or a table vs T (log mixing), gas ideal + IF97 steam compressibility (Amagat),
// enthalpies from heat capacities and a vaporisation enthalpy (Watson form), water / steam from IF97.
//
// Rock: PV = PV_0 (1 + c_r (p - p_ref) + c_T (T - T_ref)), volumetric heat capacity, thermal conductivity
// lambda_eff = phi sum S_alpha lambda_alpha + (1 - phi) lambda_rock. Heat loss to over / underburden: Vinsome-Westerveld
// semi-analytical profile (theta + p z + q z^2) exp(-z / d), d = sqrt(kappa t) / 2.
// Wells: producers on BHP or a surface-rate target (oil, water, liquid; BHP eliminated analytically, all connection cells
// coupled), injectors on a molar rate (or cold-water-equivalent steam rate) with a BHP limit, injected fluid of any
// composition at a temperature or as wet steam of a given quality at the wellbore saturation temperature.
#pragma once
#include "ressim/bsr.h"
#include "ressim/geomech.h"
#include <functional>
#include <string>
#include <vector>

namespace ressim {
namespace thermal {

// Hydrocarbon / solvent component (K-value mode)
struct Component {
    std::string name = "OIL";
    double mw = 0.4;                    // kg/mol
    double kv1 = 0, kv2 = 0, kv3 = 0, kv4 = 0, kv5 = 0;   // K = (kv1 + kv2/p + kv3 p) exp(kv4 / (T - kv5)); all 0 = non-volatile
    bool dead = false;                  // non-volatile (E300 CVTYPE DEAD / TCRIT >= 1e4) - overrides the K-value coefficients
    // component type: 0 K-value (oil / gas), 1 non-condensable gas (gas phase only: O2, N2, CO2, CO, H2, CH4, ...),
    // 2 solid (coke, char, coal, ash: immobile, occupies pore space; density rhoRef / cpL / ct1, heat capacity cpl1 / cpl2)
    int type = 0;
    bool solid() const { return type == 2; }
    // oil-phase density (E300 DREF / CREF / THERMEX1): rho = rhoRef / ((1 - cpL (p - pRef)) (1 + ct1 dT + ct2 dT^2)), dT = T - tRef
    double rhoRef = 950;                // kg/m3
    double cpL = 0;                     // 1/Pa
    double ct1 = 0, ct2 = 0;            // 1/K, 1/K^2
    double avisc = 1e-3, bvisc = 0;     // oil-phase viscosity a exp(b / T) (Andrade), Pa s / K
    std::vector<double> viscT, viscV;   // optional oil-phase viscosity table (E300 OILVISCT): T (K), mu (Pa s), log-linear in mu
    double avg = 1e-5, bvg = 0;         // gas-phase viscosity avg T^bvg (E300 GASVISCF), Pa s
    double z0 = 1, z1 = 0;              // gas Z-factor Z = z0 - z1 p / (R T) (E300 ZFACTOR / ZFACT1; z1 in m3/mol)
    // enthalpies (E300 rules, molar): oil  H = cpl1 dT + cpl2 dT^2 / 2;  gas  H = hvaps + cpg1 dT + cpg2 dT^2 / 2 - z1 (p - pst);
    // a phase without heat capacities follows from the other one and the vaporisation enthalpy
    // dH(T) = hvapA (1 - T / tcrit)^hvapB (0 above tcrit; constant hvapA when tcrit <= 0); above tcrit H_oil = H_gas
    double cpl1 = 0, cpl2 = 0;          // J/(mol K), J/(mol K^2)  (E300 SPECHA / SPECHB x MW)
    double cpg1 = 0, cpg2 = 0;          // J/(mol K), J/(mol K^2)  (E300 SPECHG / SPECHH x MW)
    double hvaps = 0;                   // gas enthalpy at tRef, J/mol (E300 HEATVAPS x MW)
    double hvapA = 0, hvapB = 0.38, tcrit = 0;   // J/mol (HEATVAP x MW), - (HEATVAPE), K (TCRIT)
    bool volatileK() const { return type == 1 || (type == 0 && !dead && (kv1 != 0 || kv2 != 0 || kv3 != 0)); }
    // ---- C10 unified multi-component framework (same physics as the compositional / CO2 cores)
    // aqueous solubility: Henry's law with the Krichevsky-Kasarnovsky (Poynting) correction and Sechenov salting-out
    //   y_i phi_i p = x_aq,i kH_i(T) exp(vInf_i (p - psat(T)) / (R T)) exp(ksalt_i m_NaCl)
    // aqModel 0 insoluble, 1 IAPWS G7-04 form (hA, hB, hC), 2 van't Hoff kH = kh298 exp(-hdT (1/T - 1/298.15)) (Pa, K).
    // Dissolution enthalpy h_aq - h_ig = -R T^2 d ln kH / dT; aqueous volume x_w v_w + sum x_i vInf_i.
    int aqModel = 0;
    double hA = 0, hB = 0, hC = 0, kh298 = 0, hdT = 0;
    double vInf = 0;                    // m3/mol (< 0: Garcia 2001 CO2 temperature fit)
    double ksalt = 0;                   // kg/mol, natural-log basis
    // Peng-Robinson gas phase (Case::gasPR): critical constants; tcG <= 0: not in the PR mixture (Z-factor volume)
    double tcG = 0, pcG = 0, omegaG = 0;
    std::string formula;                // element formula ("CO2", "CaCO3", "CH1.6"): element balance of reactions
    double trackEq = 0;                 // solids: mol of the tracked component fixed per mol (carbonates: CO2 -> mineral trapping)
    bool soluble() const { return aqModel > 0 && type != 2; }
};

// Chemical reaction (STARS-style kinetics). Rate per cell (mol of reaction extent / s):
//   r = V_b A exp(-Ea / (R T)) prod_i C_i^{n_i}               (T >= tMin, else 0)
//   reversible: r = V_b A exp(-Ea / (R T)) [prod_reac C_i^{n_i} - prod_prod C_j^{m_j} / Keq(T)], ln Keq = keqA + keqB / T
// Concentrations per bulk volume (mol/m3): fluid component i in phase j: (moles of i in phase j) / V_b; solid: N_i / V_b;
// the component 'o2pp' enters as its gas partial pressure (kPa, STARS O2PP). Component sources (p_i - s_i) r, heat
// release dH r (dH > 0 exothermic, J per mol of reaction at tRef - the component enthalpies carry the T dependence).
struct Reaction {
    std::string name;
    std::vector<double> sto, prod;      // reactant / product stoichiometric coefficients (nc + 1, water last)
    std::vector<double> order;          // concentration exponents of the reactants (nc + 1; empty = sto)
    std::vector<int> phase;             // phase of each reactant concentration: -1 auto (solid / water: aqueous /
                                        // non-condensable: gas / other: oil), 0 aqueous, 1 oil, 2 gas (nc + 1; empty = auto)
    double A = 0;                       // frequency factor (units of the concentrations)
    double Ea = 0;                      // activation energy, J/mol
    double dH = 0;                      // heat of reaction, J/mol (> 0 exothermic)
    int o2pp = -1;                      // component whose concentration is its gas partial pressure in kPa (-1: none)
    double tMin = 0;                    // K
    bool reversible = false;
    double keqA = 0, keqB = 0;          // ln Keq = keqA + keqB / T (product concentrations with exponents = prod)
};

struct Well {
    std::string name;
    bool injector = false;
    bool open = true;
    std::vector<int> cells;
    std::vector<double> wi;             // well index per connection, m3 (q_vol = WI kr / mu dp)
    double refDepth = -1;               // BHP reference depth (< 0: first connection)
    // producer: control 0 BHP, 1 surface oil (sm3/s), 2 surface water (sm3/s), 3 surface liquid (sm3/s)
    int control = 0;
    double target = 0;
    double bhp = 1e5;                   // producer BHP target / lower limit (Pa)
    // producer steam-trap (subcool) constraint, STARS OPERATE STEAMTRAP semantics: the produced fluid stays subcool K
    // below saturation, i.e. a dynamic lower BHP limit  bhp >= Psat(T_hot + subcool)  (T_hot: hottest connection cell),
    // treated implicitly (d bhp / d X of that cell in the Jacobian); <= 0: off
    double subcool = 0;
    // C17: producer steam-rate limit (STARS OPERATE MAX 'STEAM'): produced steam, cold-water equivalent m3/s (<= 0: off).
    // A further rate constraint, implicit in Newton: when it binds, the BHP is eliminated from q_steam = steamMax with the
    // vapour-phase water mobilities (the most restrictive of target / steam limit / BHP floor controls the well)
    double steamMax = 0;
    // injector: molar rate (mol/s) with an upper BHP limit, or BHP control
    double rate = 0;
    double bhpMax = 0;                  // <= 0: none
    bool bhpControl = false;
    std::vector<double> zinj;           // injected composition (nc + 1, water last)
    // injected enthalpy (E300 WINJTEMP semantics):
    //   quality >= 0 and tinj > 0: saturated at tinj; quality >= 0 and pinj > 0: saturated at Tsat(pinj);
    //   quality >= 0, neither: saturated at Tsat(BHP) (lagged); quality < 0: single phase at (tinj, pinj or BHP);
    //   hinj > 0 (J/kg of injected fluid): given enthalpy
    double tinj = 0;                    // K
    double pinj = 0;                    // Pa
    double hinj = 0;                    // J/kg
    double quality = -1;                // steam quality 0..1
};

struct Case {
    // ---- geometry: Cartesian box ...
    int nx = 1, ny = 1, nz = 1;
    std::vector<double> dx, dy, dz;     // m (size 1 = uniform / per cell)
    std::vector<double> poro, permx, permy, permz;   // -, m2 (size 1 = uniform / per cell)
    double top = 0;
    // ... or general connections (cellPV / cellDepth / cellBulk non-empty)
    struct Conn { int a, b; double T; double G; };   // T: transmissibility (m3), G: conduction geometry A / d (m)
    std::vector<double> cellPV, cellDepth, cellBulk;
    std::vector<Conn> conns;
    // heat-loss faces (cell, area m2); Cartesian boxes: top / bottom layer faces are added when heatLoss is on
    struct LossFace { int cell; double area; };
    std::vector<LossFace> lossFaces;
    bool heatLoss = false;
    double obCv = 2.35e6, obK = 1.73;   // over / underburden volumetric heat capacity (J/m3/K) and conductivity (W/m/K)
    double obT = 0;                     // initial over / underburden temperature (K); <= 0: initial temperature of the cell
    // ---- rock
    double cr = 0, ct = 0, pref = 1e5;  // PV compressibility (1/Pa) and thermal expansion (1/K)
    double rockCv = 2.35e6;             // rock (solid) volumetric heat capacity, J/(m3 K) (E300 HEATCR, per rock volume)
    double rockCvT = 0;                 // its temperature coefficient, J/(m3 K^2) (E300 HEATCRT): H_r = Cv dT + CvT dT^2 / 2
    std::vector<double> rockCvCell, rockCvTCell;   // optional per-cell values (override rockCv / rockCvT)
    // thermal conductivity (E300 THCONR option): cell conductivity THCONR (W/m/K) x (1 - thconsf S_g); when thconr is empty
    // THCONR = phi (lamO + lamW) / 2 + (1 - phi) lamRock and thconsf (< 0) = phi ((lamO + lamW) / 2 - lamG) / THCONR
    double lamRock = 1.7, lamW = 0.6, lamO = 0.13, lamG = 0.04;   // W/(m K)
    std::vector<double> thconr;         // size 1 = uniform / per cell
    double thconsf = -1;
    std::vector<double> thconsfCell;    // optional per-cell THCONSF (overrides thconsf)
    // water / steam viscosity: optional water table (E300 WATVISCT; else IAPWS-2008), steam Ejiogu-Fiori
    // mu = A + B Tc + C Pp^D (cP, deg C, MPa; E300 THSVC) or IAPWS-2008
    std::vector<double> watViscT, watViscV;   // K, Pa s
    double steamVisc[4] = {4.9402e-3, 5.0956e-5, 2.9223e-6, 2.5077};
    bool steamViscIapws = false;
    // ---- saturation functions (Corey; ECLIPSE default three-phase oil kr)
    double swc = 0.2, sorw = 0.2, sorg = 0.1, sgc = 0.05;
    double krw0 = 0.3, kro0 = 1.0, krg0 = 0.5, nw = 2, now = 2, nog = 2, ng = 2;
    std::vector<double> swT, krwT, krowT, slT, krgT, krogT;   // optional SWT / SLT style tables (liquid sat for SLT)
    // capillary pressure (p = oil-phase pressure; p_w = p - Pcow(Sw), p_g = p + Pcgo(Sg)): tables aligned with swT / slT
    // (Pa), else linear Corey-mode curves Pcow = pcw0 (1 - Swn), Pcgo = pcg0 Sgn; wells use the oil-phase pressure
    std::vector<double> pcowT, pcgoT;
    double pcw0 = 0, pcg0 = 0;
    // ---- chemistry (reactions, solids) and permeability alteration
    std::vector<Reaction> reactions;
    std::vector<std::vector<double>> solid0;   // initial solid concentrations, mol per m3 of pore volume (nc; size 1 = uniform)
    // permeability vs fluid porosity phi_f = (PV - V_solid) / V_b (solids, geomechanics), lagged per time step:
    // 0 off, 1 Kozeny-Carman (phi/phi0)^3 ((1 - phi0) / (1 - phi))^2, 2 exponential exp(permC (phi - phi0))
    int permLaw = 0;
    double permC = 20, permMultMax = 1e5;
    // ---- geomechanics (THM / THMC; Cartesian boxes): linear thermo-poroelastic FEM, fixed-stress coupling
    geomech::Input geo;
    // general-connection grids (cellPV non-empty): the Cartesian box of the mechanics mesh and the cell -> box index map
    int geoNx = 0, geoNy = 0, geoNz = 0;
    std::vector<double> geoDx, geoDy, geoDz;   // per box cell or size 1 (m)
    double geoTop = 0;
    std::vector<int> geoCell;
    std::vector<double> geoNodes;   // C25 corner-point / irregular grids: box node coordinates (m) -> isoparametric hexahedra
    // ---- fluid
    std::vector<Component> comps;       // hydrocarbon / solvent components (water is implicit, last)
    double tRef = 288.15, pRef = 101325.0;   // enthalpy / density reference
    double stdT = 288.71, stdP = 101325.0;
    // ---- initial state (size 1 = uniform)
    std::vector<double> p0, T0, sw0, sg0;
    std::vector<std::vector<double>> x0, y0;   // initial oil / gas compositions (nc); y0 empty = x0
    std::vector<Well> wells;
    bool gravity = true;
    // ---- C10: unified multi-component framework
    double salinity = 0;                // NaCl molality of the aqueous phase (salting-out of dissolved gases), mol/kg
    bool gasPR = false;                 // Peng-Robinson gas phase (volume, fugacity, residual enthalpy, LBC viscosity)
    // gas relative-permeability hysteresis (residual trapping), Land trapping with the Carlson free-gas krg (as the
    // compositional core): S_gt = S_gc + (S_gmax - S_gc) / (1 + C (S_gmax - S_gc)), C = 1/(sgtMax - S_gc) - 1/(1 - S_wc - S_gc)
    struct GasHyst { bool on = false; double sgtMax = 0.3; } hyst;
    int track = -1;                     // tracked component of the storage / trapping inventory (-1: none)
    double structTol = 0.01;            // inventory: mobile gas is structurally held when its outflow < structTol x volume / year
    // ---- C12: embedded discrete fractures (filled by addFractures): fracture cells with Barton-Bandis / Coulomb mechanics
    struct FracCell {
        int cell = -1, frac = -1, host = -1;   // fracture cell, fracture index, host matrix cell
        double w0 = 1e-3;               // initial hydraulic aperture (m)
        double segLen = 0, segH = 0, kf = 0;   // trace length and height of the segment in its host (m), permeability (m2)
        double sn0 = 0, tau0 = 0;       // initial total normal / shear stress on the plane (Pa, compression +)
        double p0 = 0, T0 = 0;          // initial pressure / temperature (Pa, K)
        double wOpen = 1e-3;            // Barton-Bandis aperture at zero effective stress (m)
        double kn0 = 0, vm = 0;         // Pa/m, m (kn0 <= 0: rigid)
        double mu = 0.6, ks = 0, tanPsi = 0;   // friction, shear stiffness (Pa/m; 0: no hydro-shear), tan(dilation)
        double thermStress = 0;         // d sigma_n / dT (Pa/K): E alpha / (1 - nu)
        double wShear = 0;              // accumulated (irreversible) shear dilation (m)
        double xa = 0, ya = 0, xb = 0, yb = 0, zt = 0, zb = 0;   // C13: segment geometry (m from the box corner, depth) for display
    };
    std::vector<FracCell> fracCells;
    std::vector<int> fracFF;            // indices of fracture-fracture connections in conns (scaled by the cubic law)
    int size() const { return cellPV.empty() ? nx * ny * nz : int(cellPV.size()); }
    bool aqueousSolutes() const { for (const Component& q : comps) if (q.soluble()) return true; return false; }
};

// storage / trapping inventory of the tracked component (kg), same classes as the CO2 cores (Co2Inventory)
struct Inventory {
    double dissolvedAq = 0;             // solubility trapping (aqueous phase)
    double dissolvedOil = 0;            // dissolved in the oil phase
    double residual = 0;                // gas phase up to the trapped saturation (hysteresis S_gt, else S_gc)
    double immobileKr = 0;              // gas phase in cells with krg = 0
    double structural = 0;              // mobile gas at rest beneath a seal
    double mobile = 0;                  // migrating gas
    double mineral = 0;                 // fixed in solids (net change since t = 0, trackEq of the solid components)
    double total() const { return dissolvedAq + dissolvedOil + residual + structural + mobile + mineral; }
};

struct Stats { int steps = 0, newton = 0, linear = 0, cuts = 0, flashes = 0, mechSolves = 0, mechIters = 0, mechCg = 0; };

struct WellReport {
    double bhp = 0;
    bool onBhp = false;
    bool trapActive = false;            // BHP held by the steam-trap (subcool) constraint
    bool steamLimited = false;          // C17: rate held by the steam-rate limit
    double qo = 0, qw = 0, qg = 0;      // surface production rates (sm3/s)
    double qinj = 0;                    // injection (mol/s); for water: cold-water equivalent qcwe (m3/s)
    double qcwe = 0;
    double hRate = 0;                   // enthalpy rate (W): injection < 0, production > 0
    double tFluid = 0;                  // C13: well-stream temperature (K): producers molar-rate weighted over the
                                        // connections, injectors the injection temperature (0 = no flow)
    std::vector<double> molar;          // component rates (mol/s, production > 0), nc + 1
    double cumO = 0, cumW = 0, cumG = 0, cumCwe = 0, cumH = 0;
    double qsteam = 0, cumSteam = 0;    // C16: produced steam (water leaving in the vapour phase), cold-water equivalent m3/s, m3
    std::vector<double> cumMolP, cumMolI;   // cumulative produced / injected moles per component (nc + 1, water last)
};

// C18 discrete adjoint (Simulator::adjointGradient). Observation qty: 0 surface oil, 1 surface gas, 2 surface water
// (sm3/s, production), 3 BHP (Pa), 4 enthalpy rate (W, production > 0), 6 well-stream temperature (K), 10 + k molar
// rate of component k (water = nc; mol/s, production > 0)
struct AdjObs { int well = 0; int qty = 0; double t = 0, value = 0, sigma = 1; };
// parameter: apply(case, m) scales a property by the multiplier m (base 1); gradients are dF/dm at m = 1
struct AdjParam { std::string name; std::function<void(Case&, double)> apply; double h = 1e-4; };
struct AdjResult { bool ok = false; std::string error; double F = 0; std::vector<double> grad; int steps = 0, linear = 0; double seconds = 0; };
AdjParam paramPermMult(const std::string& name, const std::vector<int>& cells);     // transmissibilities (box: k; grid: T)
AdjParam paramPoroMult(const std::string& name, const std::vector<int>& cells);     // pore volume (box: phi; grid: PV)
AdjParam paramWellIndexMult(const std::string& name, int well);
AdjParam paramConductivityMult(const std::string& name);                           // rock / fluid thermal conductivities

// Rachford-Rice free-water flash at (p, T) (exposed for tests / tools)
struct FlashState {
    double nW = 0, nO = 0, nG = 0;      // phase moles (per mole of feed when called with z)
    std::vector<double> x, y;           // oil (nc) and gas (nc + 1, steam last) mole fractions
    std::vector<double> a;              // aqueous mole fractions (nc + 1, water last; dissolved gases with Component::aqModel)
};
// general three-phase K-value flash with aqueous solutes (used by flashPT when Case::aqueousSolutes()): convex
// minimisation of Q(beta) = sum_j beta_j - sum_i z_i ln(sum_j beta_j a_ij) over beta >= 0 (Michelsen 1994), affinities
// a_ij = 1 / H_ij from the gas-referenced fugacity constants; outer successive substitution on the gas fugacity
// coefficients (Peng-Robinson). Exposed to cross-check it against the free-water flash.
FlashState flashGeneral(const Case& c, double p, double T, const std::vector<double>& z);
FlashState flashPT(const Case& c, double p, double T, const std::vector<double>& z);

class Simulator {
public:
    bool setup(const Case& c, std::string* error = nullptr);
    bool advanceTo(double t, double dtMax);   // seconds
    double time() const { return t_; }
    int size() const { return n_; }
    int nc() const { return nc_; }
    std::vector<Well>& wells() { return c_.wells; }
    const Case& model() const { return c_; }
    // cell state
    const std::vector<double>& pressure() const { return p_; }
    const std::vector<double>& temperature() const { return T_; }
    const std::vector<double>& sw() const { return sw_; }
    const std::vector<double>& so() const { return so_; }
    const std::vector<double>& sg() const { return sg_; }
    const std::vector<double>& moles() const { return X_; }   // state vector (n * nb)
    const std::vector<WellReport>& wellReports() const { return rep_; }
    const Stats& stats() const { return st_; }
    // field balances
    double fluidMoles(int k) const;     // in place (k = nc: water)
    double energyInPlace() const;       // J (fluids + rock, relative to the reference states)
    double cumHeatLoss() const { return cumLoss_; }   // J
    double cumEnthalpyIn() const { return cumHin_; }  // J (injected - produced)
    double cumReactionHeat() const { return cumReact_; }   // J released by reactions
    double reactionHeatRate() const { return reactRate_; } // W (last step)
    // porosities: total pore volume (void + solids) and fluid (void) fraction of the bulk volume
    std::vector<double> porosity() const;
    std::vector<double> fluidPorosity() const;
    double bulkVolume(int q) const { return bulk_[size_t(q)]; }   // m3
    const std::vector<double>& permMultiplier() const { return permMult_; }
    // C25: transmissibility multipliers (absolute w.r.t. the base) of the connections between active-cell pairs (fault
    // mechanics two-way coupling); composed with the porosity / damage permeability law
    void setConnectionMultipliersByCells(const std::vector<std::pair<int, int>>& cells, const std::vector<double>& mult);
    std::vector<double> fractureApertures() const;   // current hydraulic aperture per Case::fracCells entry (m)
    double solidMoles(int k) const { return fluidMoles(k); }
    const geomech::Solver* mechanics() const { return geo_.ready() ? &geo_ : nullptr; }
    // C10: phase compositions per cell (x oil, y gas, a aqueous; each n * (nc + 1), water last), gas hysteresis state,
    // storage / trapping inventory of the tracked component (Case::track)
    void phaseCompositions(std::vector<double>& x, std::vector<double>& y, std::vector<double>& a) const;
    const std::vector<double>& sgMax() const { return sgHist_; }
    std::vector<double> trappedGasSaturation() const;
    Inventory trackedInventory() const;
    const std::vector<int>& mechCellMap() const { return gmap_; }   // cell -> mechanics box cell
    double stockTankOilInPlace() const; // sm3
    double dtLast() const { return dtLast_; }
    void limitNextDt(double dt) { dt_ = std::max(1.0, std::min(dt_, dt)); }   // seconds (well-control changes)
    // diagnostics: max relative difference between the assembled Jacobian and central differences of the residual at
    // the current state for a step dt (columns of all cells; rows scaled by the row norm)
    double jacobianCheck(double dt, int* worstRow = nullptr, int* worstCol = nullptr, int onlyCell = -1);   // onlyCell >= 0: columns of that cell
    // C14 restart: complete dynamic state (moles / energy, heat-loss memory, fracture apertures and shear dilation,
    // mechanics, permeability / pore-volume laws, well reports, lagged yields, hysteresis, cumulatives, step size);
    // loadState() follows setup() of the same model (same cell / well / fracture counts). Throws on a mismatch.
    std::string saveState() const;
    void loadState(const std::string& blob);
    // ---- C18 adjoint: record converged states (call after setup, before the first advanceTo)
    void setAdjointRecording(bool on);
    int adjointSteps() const { return adjRec_.empty() ? 0 : int(adjRec_.size()) - 1; }
    std::vector<double> stepSizes() const;
    void setForcedSteps(const std::vector<double>& dts) { forced_ = dts; forcedIdx_ = 0; }
    // Newton tolerances: component residual / cell moles, energy residual / heat capacity (K), volume residual / PV
    void setNewtonTolerance(double tolM, double tolE, double tolV) { tolM_ = tolM; tolE_ = tolE; tolV_ = tolV; }
    double adjointMisfit(const std::vector<AdjObs>& obs);
    AdjResult adjointGradient(const std::vector<AdjObs>& obs, const std::vector<AdjParam>& params);
private:
    struct LossState;
    struct AdjRec;
    void adjRecord(double t0, double dt);
    void adjEval(const AdjRec& prev, const AdjRec& cur, const Simulator& base, std::vector<double>& R, BSR* J);
    double adjObsValue(const AdjObs& o, const WellReport& r) const;
    std::vector<AdjRec> adjRec_;
    bool adjOn_ = false;
    std::vector<double> forced_; size_t forcedIdx_ = 0;
    double tolM_ = 1e-6, tolE_ = 1e-3, tolV_ = 1e-6;
    template <class Ar> void stateIO(Ar& ar);
    struct LossState { double I = 0, thetaN = 0, Ti = 0; };
    struct AdjRec {
        double t0 = 0, t = 0, dt = 0;
        std::vector<double> X, T, sgHist, permMult, pvA, pvCp, pvP, pvCt, pvT, connT, fracW, wbRho, wellRamp;
        std::vector<LossState> loss; std::vector<std::vector<double>> surfYield; std::vector<std::vector<char>> wellAct;
        std::vector<int> wellMode, steamMode; std::vector<WellReport> rep; std::vector<Well> wells;
    };
    Case c_;
    int n_ = 0, nc_ = 0, nb_ = 0, m_ = 0;
    double t_ = 0, dt_ = 8640, dtLast_ = 0, cumLoss_ = 0, cumHin_ = 0;
    std::vector<double> pv0_, depth_, bulk_, crock_, crock1_;   // reference PV, depth, bulk volume, rock heat capacity (J/K, J/K2)
    std::vector<double> thc_, thcSf_;   // cell conductivity THCONR (W/m/K) and its gas-saturation factor
    // pore volume PV(p, T) = pvA + pvCp (p - pvP) + pvCt (T - pvT): rock compressibility / thermal expansion, or the
    // fixed-stress form centred at the last mechanics state (geomechanics)
    std::vector<double> pvA_, pvCp_, pvP_, pvCt_, pvT_;
    double poreVolume(int q, double p, double T) const { const size_t i = size_t(q); return pvA_[i] + pvCp_[i] * (p - pvP_[i]) + pvCt_[i] * (T - pvT_[i]); }
    std::vector<double> fracW_;         // C12: current fracture apertures (Case::fracCells order)
    void updateFractures(bool stepEnd);
    std::vector<double> connX_;    // C25 external connection multipliers (empty = 1)
    std::vector<double> phiF0_, permMult_, connT0_;   // initial fluid porosity, permeability multiplier, base transmissibility
    std::vector<double> p0_, T0init_;   // initial state (mechanics reference)
    geomech::Solver geo_;
    std::vector<int> gmap_;             // cell -> mechanics box cell
    int gn_ = 0;                        // mechanics box cells
    double cumReact_ = 0, reactRate_ = 0;
    std::vector<double> sgHist_;        // historical maximum gas saturation (hysteresis, lagged one step)
    std::vector<double> solidInit_;     // initial solid moles per cell (n * nc) - mineral inventory
    double trappedGas(int q, double sgHist) const;   // Land trapped gas saturation of the scanning curve from sgHist
    int rbase_ = 0, stride_ = 0;        // Y layout: reaction rows, phase block stride
    bool hasSolid_ = false, hasPc_ = false;
    double solidVolume(int q, const double* X, double p, double T) const;
    bool solveStep(double dt, int& its);
    bool mechanicsUpdate(bool iterate, double& change);
    void updatePermeability();
    double rockEnergy(int q, double T) const { const double d = T - c_.tRef; return crock_[size_t(q)] * d + 0.5 * crock1_[size_t(q)] * d * d; }
    std::vector<Case::Conn> conns_;
    std::vector<double> X_, Xn_;        // state (n * nb)
    std::vector<double> p_, T_, sw_, so_, sg_;
    std::vector<double> Y_, dY_;        // local outputs (n * m) and derivatives (n * m * nb)
    std::vector<LossState> loss_;
    std::vector<WellReport> rep_;
    std::vector<std::vector<double>> surfYield_;   // per well: surface volumes per mole of component (oil, water) lagged
    std::vector<double> wbRho_;         // wellbore density per well (lagged)
    Stats st_;
    BSR J_;                             // Jacobian (pattern built on the first Newton iteration)
    bool assembleOnly_ = false;
    int newtonIt_ = 0;                  // Newton iteration within the current attempt (relaxed tolerances, oscillation control)
    double relax_ = 1.0;                // update relaxation (OPM-style: reduced when the residual oscillates / stagnates)
    std::vector<double> resHist_;       // residual norms of the current attempt
    std::vector<std::vector<char>> wellAct_;   // active (flowing) connections per well, previous iteration (hysteresis)
    int sinceCut_ = 100;                       // accepted steps since the last time-step cut
    std::vector<int> steamMode_, steamSw_;     // C17: per well within a time step: steam limit binding (0/1), switches
    std::vector<int> wellMode_, wellSw_;       // per well within a time step: control mode (-1 / 0 rate / 1 BHP), switches
    std::vector<double> wellRamp_;             // per injector: width of the C1 drawdown ramp (Pa), lagged per attempt
    std::vector<double> lastR_, lastJ_; // unscaled residual / Jacobian of the last assembly
    struct Local;
    bool evalCell(int c, const double* x, double* y, Local* L, double Tguess) const;
    bool evalAll();
    bool newtonStep(double dt, bool& converged, double& resNorm);
    void finishStep(double dt);
    void surfaceVolumes(const std::vector<double>& molar, double& vo, double& vw, double& vg) const;
};

// ---- C12: embedded discrete fracture model (EDFM; Li & Lee 2008, Moinfar et al. 2014) for Cartesian boxes
// Vertical planar fractures (trace (x0, y0)-(x1, y1) in the box plane, depth range) are cut by the matrix columns; each
// (column, layer) piece becomes a fracture cell (PV = w A) connected to its host matrix cell by T_fm = 2 A k_n / <d>
// (<d>: mean distance of the host cell from the plane), to its neighbours in the same fracture (half-transmissibility
// harmonic, k_f w h / (L/2)) and vertically. Fracture-fracture intersections share a connection through the common cell
// pieces when two fractures cut the same host cell (T = harmonic of the two half transmissibilities over the
// intersection length). Mechanics: Barton-Bandis closure w = w_open - Vm s' / (Kn0 Vm + s') with s' = s_n - p and a
// thermoelastic change of s_n (E alpha / (1 - nu) (T - T0): cooling unloads the plane), Coulomb hydro-shear: slip when
// tau > mu s', irreversible dilation (tau - mu s') / ks tan(psi). The aperture scales the fracture transmissibilities
// (cubic law) after each step and the fracture pore volume implicitly (linearised in p).
struct Fracture {
    std::string name = "F";
    double x0 = 0, y0 = 0, x1 = 0, y1 = 0;   // m from the box corner (x along i, y along j)
    double zTop = -1, zBot = -1;        // depth range (m, absolute as Case::top); < 0: the whole box thickness
    double aperture = 1e-3;             // hydraulic aperture at the initial state (m)
    double perm = 0;                    // m2 (0: cubic law w^2 / 12)
    bool mech = false;                  // Barton-Bandis / Coulomb aperture update
    double kn0 = 20e9;                  // initial normal stiffness (Pa/m)
    double vm = 0.4e-3;                 // maximum closure (m)
    double mu = 0.6;                    // friction coefficient
    double ks = 0;                      // shear stiffness (Pa/m); 0 = no hydro-shear
    double dilationDeg = 5;             // shear dilation angle
};
struct FractureStress {
    double sHmaxGrad = 0, sHminGrad = 0;   // horizontal principal total stresses (Pa per m of depth)
    double azimuthDeg = 0;              // direction of SHmax measured from the x axis (degrees)
    double E = 30e9, nu = 0.25, alphaT = 0;   // thermoelastic unloading of the plane (alphaT = 0: off)
};
// Converts a Cartesian case to the general-connection form and appends the fracture cells / connections
// (per-cell arrays are extended with the host values; geomechanics keeps the matrix cells, fracture cells unmapped)
bool addFractures(Case& c, const std::vector<Fracture>& fr, const FractureStress& st, std::string* error = nullptr);
// Well - fracture connections (EDFM, Moinfar et al. 2013): every completion in a host cell of fracture segments also
// connects to those fracture cells with WI_f = 2 pi k_f w / ln(r_e / r_w), r_e = 0.14 sqrt(L_s^2 + H_s^2). Called by
// addFractures for c.wells; call it for well lists defined later (deck schedules). Returns the number of connections added.
int connectWellsToFractures(const Case& c, std::vector<Well>& wells, double rw = 0.1);

// ---- C12: Langmuir sorption (gas adsorption on coal / char / shale, solute retention: nanoparticles, polymer)
// Mass-action form on a shared site species (immobile solid 'site'):  C + Site <=> C_ads, rate per bulk volume
//   r = A [c_C [Site] - [C_ads] / b(T)],   b(T) = bRef exp(dH / R (1/T - 1/tRef))   (van 't Hoff, dH > 0 exothermic)
// At equilibrium with [Site] + sum_j [C_ads,j] = capacity this is the (extended, competitive) Langmuir isotherm
//   [C_ads,i] = capacity b_i c_i / (1 + sum_j b_j c_j)
// c_C: gas partial pressure in kPa (phase 2 / default for non-condensable gases) or the aqueous moles per m3 of bulk
// volume (phase 0). A = b / (2 tau) gives the relaxation time tau (s) at c = 1 / b. The sorbed species keeps the
// component's molar mass; adsRho > 0 lets it occupy pore space (kg/m3: deposition, permeability via PERMLAW),
// 0 = no volume. Site and sorbed species are appended as solids before water; a site name already present is shared.
struct Sorption {
    std::string comp;                   // sorbing component (name)
    double capacity = 0;                // site capacity, mol per m3 of bulk rock (Langmuir V_L rho_rock / v_STP)
    double bRef = 0;                    // Langmuir constant at tRef (1/kPa for gas, m3 bulk / mol for aqueous)
    double tRef = 298.15;               // K
    double dH = 0;                      // heat of adsorption, J/mol (> 0 exothermic; also released by the reaction)
    double tau = 3600;                  // relaxation time, s (kinetic Langmuir; small = equilibrium)
    int phase = -1;                     // -1 auto (gas for non-condensable / K-value, aqueous for soluble-only), 0 aq, 2 gas
    double adsRho = 0;                  // sorbed-phase density, kg/m3 (0: no pore-volume occupancy)
    std::string site = "SITES";
    bool initEquilibrium = true;        // initial sorbed loading in equilibrium with the initial gas (p0 y0)
};
// Appends the site / sorbed solids and the reaction to the case (resizes every per-component array, wells included)
bool addSorption(Case& c, const Sorption& s, std::string* error = nullptr);

} // namespace thermal
} // namespace ressim
