// Compositional (EOS) reservoir simulator core - milestones C2 / C3.
//
// Formulation (ECLIPSE 300 FIM / Acs-Coats-Watts volume balance):
//   primary unknowns per cell   p (oil-phase pressure), N_1 .. N_nc (component moles) [, N_w water moles]
//   equations per cell          component balances        N_i - N_i^n + dt * (sum_faces F_i + q_i) = 0
//                               volume balance            PV(p) - V_L - V_V [- V_W] = 0
//   F_i = T * sum_alpha [ kr_alpha / mu_alpha * xi_alpha * x_i,alpha ]_upstream * dPhi_alpha,
//   dPhi_alpha = (p_a + Pc_alpha,a) - (p_b + Pc_alpha,b) + rho_alpha g (D_b - D_a)   (phase-wise upwinding),
//   Pc_water = -Pcow(Sw), Pc_gas = +Pcgo(Sg), Pc_oil = 0.
// Phase behaviour comes from a PT flash in every cell and Newton iteration (eos::flash, K-value warm start, stability
// test for single-phase cells). The Jacobian is analytic: phase fractions / compositions are differentiated with the
// implicit-function theorem (eos::flashSensitivities), molar volumes with the EOS partial molar volumes, viscosities
// with the analytic LBC derivatives, saturation functions by their table slopes. The block system (nb = nc + 1 or
// nc + 2) is solved with GMRES + BlockCPR (true-IMPES pressure weights).
//
// Water: either an immobile connate saturation (water = false), or a mobile immiscible aqueous phase (water = true):
// PVTW-style water (Bw = Bw,ref / (1 + X + X^2/2), X = cw (p - pref,w)).
// Saturation functions: Corey curves or ECLIPSE SWFN / SGFN / SOF3 tables (incl. Pcow / Pcgo); three-phase oil kr by the
// ECLIPSE default model kro = (Sg krog + (Sw - Swc) krow) / (Sg + Sw - Swc) or Stone II.
//
// Wells (standard well model, no crossflow): connection flow = WI * molar mobility * (Phi_conn - p_wellbore,conn) with
// p_wellbore,conn = BHP + rho_well g (D_conn - D_ref) (wellbore density lagged one time step, ECLIPSE practice). Under a
// rate target the BHP is eliminated analytically (exact Schur complement of the well equation: all connection cells of
// the well are coupled in the Jacobian) and the well switches to its BHP limit inside Newton. Producer targets: BHP,
// hydrocarbon molar rate, reservoir volume, stock-tank oil, surface gas, surface water, liquid; surface rates come from
// the separator train (SEPCOND / FIELDSEP; per-mole conversion factors lagged one time step, E300 practice). Injectors:
// molar rate of any composition (WELLSTRE / WINJGAS) or water, total-mobility injectivity, BHP limit.
#pragma once
#include "ressim/eos.h"
#include "ressim/co2brine.h"
#include <functional>
#include <string>
#include <vector>

namespace ressim {
struct BSR;
namespace comp {

struct CompWell {
    std::string name;
    bool injector = false;
    bool open = true;
    std::vector<int> cells;             // perforated cells (cell index of the simulator)
    std::vector<double> wi;             // well index per connection, m3 (q_vol = WI * kr / mu * dp)
    double refDepth = -1;               // BHP reference depth (m); < 0 = depth of the first connection
    // producer: control 0 = BHP, 1 = hydrocarbon molar rate (mol/s), 2 = reservoir volume of all phases (m3/s),
    // 3 = stock-tank oil (sm3/s), 4 = surface gas (sm3/s), 5 = surface water (sm3/s), 6 = surface liquid (sm3/s)
    int control = 0;
    double target = 0;
    double bhp = 0;                     // producer: BHP target (control 0) / lower BHP limit (rate controls), Pa
    // injector: molar rate target (mol/s) with an upper BHP limit, or BHP control
    double rate = 0;
    double bhpMax = 0;                  // injector BHP limit (<= 0: none) / BHP target when injBhpControl
    bool injBhpControl = false;
    bool water = false;                 // injector: injects water (requires CompCase::water) instead of zinj
    std::vector<double> zinj;           // injector: injected composition
    // producer: further rate limits (control code as above, value in the same units as target); the most restrictive
    // constraint (highest BHP) is active, chosen inside Newton (E300 / ECLIPSE semantics)
    std::vector<std::pair<int, double>> limits;
    int sep = 0;                        // separator train (WSEPCOND): 0 = field separator, k > 0 = CompCase::sepTrains[k-1]
};

// ECLIPSE saturation-function tables (family II): any of them may be empty (Corey parameters are used instead)
struct SatTables {
    std::vector<double> sw, krw, pcow;  // SWFN (Pcow in Pa)
    std::vector<double> sg, krg, pcog;  // SGFN (Pcgo in Pa)
    std::vector<double> so, krow, krog; // SOF3
    bool empty() const { return sw.empty() && sg.empty() && so.empty(); }
};

struct CompCase {
    // ---- geometry: either a Cartesian box ...
    int nx = 1, ny = 1, nz = 1;
    std::vector<double> dx, dy, dz;     // m, per cell (size 1 = uniform)
    std::vector<double> poro, permx, permy, permz;   // -, m2 (size 1 = uniform)
    double top = 0;                     // depth of the top of layer 1 (m)
    // ... or a general unstructured cell list (corner point, faults, NNC, LGR, from the black-oil grid)
    struct Conn { int a, b; double T; double Td = 0; double Tdisp = 0; };   // T: transmissibility (m3); Td: diffusive phi A / d (m), 0 = from T;
                                                                           // Tdisp: dispersivity / centre distance (-), DISPERC
    double dispersivity = 0;            // box grid: uniform mechanical dispersivity (m), alpha_L = alpha_T
    std::vector<double> cellPV, cellDepth; // reference pore volume (m3) and centre depth (m) per cell
    std::vector<Conn> conns;
    // ---- rock and saturation functions
    double cr = 0, pref = 1e5;          // rock compressibility (1/Pa) and reference pressure
    double swc = 0;                     // connate water saturation (SWFN: its first saturation)
    double sorg = 0, sgc = 0, kro0 = 1, krg0 = 1, no = 2, ng = 2;   // Corey (oil in gas, gas)
    bool water = false;                 // mobile water phase
    // water vaporisation (OPM VAPWAT semantics in EOS mode): water vapour in the gas phase, free-water Raoult law
    // y_w = Psat(T) / p (IAPWS-IF97); the vapour takes y_w / (1 - y_w) moles of water per hydrocarbon vapour mole while
    // the aqueous phase lasts (dry-out: all water vaporised). Needs the mobile water phase.
    bool waterVap = false;
    double rhoWsc = 1000, bwRef = 1.0, cw = 4.5e-10, pwRef = 1e5, muw = 5e-4;   // kg/m3, rm3/sm3, 1/Pa, Pa, Pa s
    double sorw = 0, krw0 = 1, nw = 2, now = 2;   // Corey (water, oil in water)
    SatTables sat;                      // SATNUM region 1
    std::vector<SatTables> satRegions;  // SATNUM regions 2..N
    std::vector<int> satnum, imbnum;    // per cell: 0-based table index into {sat, satRegions} (empty = 0 / = satnum)
    // gas relative-permeability hysteresis (residual trapping): model 0 = Carlson / Land free gas, 1 = Killough
    struct Hysteresis { bool on = false; int model = 1; } hyst;
    std::vector<double> cellPoro;       // reference porosity per cell (Millington-Quirk tortuosity; empty = 0.2)
    std::vector<double> diffOil, diffGas;   // EOS mode (DIFFC): component diffusion coefficients in oil / gas, m2/s
    int kroModel = 0;                   // 0 ECLIPSE default three-phase oil kr, 1 Stone II
    // ---- fluid
    eos::Model fluid;
    double T = 350;                     // K (isothermal)
    std::vector<std::pair<double, double>> separator;   // stages (T K, p Pa); E300 semantics: the last stage is the stock tank (empty = one flash at standard conditions)
    std::vector<std::vector<std::pair<double, double>>> sepTrains;   // SEPCOND trains used by WSEPCOND (CompWell::sep)
    double stdT = 288.71, stdP = 101325.0;
    // ---- initial state
    std::vector<double> p0;             // initial pressure (size 1 = uniform), Pa
    std::vector<std::vector<double>> z0;   // initial composition (size 1 = uniform)
    std::vector<double> sw0;            // initial water saturation (size 1 = uniform; empty = swc)
    std::vector<CompWell> wells;
    bool gravity = true;
    // ---- CO2STORE (milestone C4): components {CO2, H2O}, phases aqueous (liquid slot) + CO2-rich (vapour slot), NaCl
    // as a fixed molality; thermodynamics from co2brine.h (fluid is only used for names / molar masses). Saturation
    // functions: brine krw(Sw) (sat.sw / krw), CO2 krg(Sg) (sat.sg / krg), Pc_gas = Pcow(Sw) + Pcgo(Sg); water = false.
    bool co2store = false;
    double salinity = 0;                // NaCl molality (mol / kg H2O)
    co2::TrappingOptions trap;          // diffusion, salt precipitation, minerals, Land sgtMax, structural tolerance
    // ---- CO2SOL (E300 style, C6): CO2 is an ordinary EOS component that also dissolves in the aqueous phase (water =
    // true). Aqueous CO2 from fugacity equality f_CO2(HC phases, EOS) = x_aq H(p), H(p) = f_CO2^pure(p) / x_sat(p) with
    // x_sat from Spycher-Pruess / Duan-Sun at the salinity (co2brine); aqueous volume N_w v_w(p) + s V_phi (Garcia),
    // viscosity mu_w (1 + 4.65 x^1.0134) (Islam-Carlson). solSalinity: NaCl molality.
    bool co2sol = false;
    int co2Comp = -1;                   // index of CO2 in the EOS component list (resolved by name when < 0)
    double solSalinity = 0;
    std::vector<double> solTabP, solTabX; // SOLUBILI: pressure (Pa) and saturated aqueous CO2 mole fraction (empty = default)
    std::vector<double> solLnP, solLnH; // derived by the simulator setup: ln p (Pa) grid and ln H (Pa)
    int size() const { return cellPV.empty() ? nx * ny * nz : int(cellPV.size()); }
};

// CO2 inventory by trapping mechanism (kg CO2) and the precipitated salt (kg)
struct Co2Inventory {
    double dissolved = 0;               // solubility trapping: CO2 in the brine
    double residual = 0;                // residual trapping: gas-phase CO2 up to the trapped saturation (hysteresis Sgt / Sgcr)
    double immobileKr = 0;              // gas-phase CO2 in cells with krg = 0 (strictly immobile now)
    double structural = 0;              // mobile gas-phase CO2 at rest beneath a seal (outflow < structTol per year)
    double mobile = 0;                  // remaining (migrating) gas-phase CO2
    double mineral = 0;                 // mineral trapping: CO2 fixed in carbonates
    double saltSolid = 0;               // precipitated NaCl, kg
    double total() const { return dissolved + residual + structural + mobile + mineral; }
};

struct CompStats { int steps = 0, newton = 0, linear = 0, cuts = 0, flashes = 0, stabilityTests = 0, relaxed = 0; };

// well results of the last converged time step (production > 0)
struct CompWellReport {
    double bhp = 0;
    bool onBhp = false;
    int control = 0;                    // active control code (rate mode) / 0 on BHP
    double qo = 0, qg = 0, qw = 0;      // surface production rates: stock-tank oil, gas, water (sm3/s)
    double qgi = 0, qwi = 0;            // surface injection rates (sm3/s)
    double resv = 0;                    // reservoir volume rate (m3/s, production > 0)
    std::vector<double> molar;          // component rates (mol/s, production > 0, injection < 0), ncw
    double cumO = 0, cumG = 0, cumW = 0, cumGI = 0, cumWI = 0;   // sm3
    std::vector<double> cumMolP, cumMolI;                         // cumulative produced / injected moles per component (ncw)
};

// group rate target (GCONPROD / GCONINJE RATE): shared among the member wells under group control by guide rates
// (productivity x drawdown to the well's BHP limit x the well's last surface yield of the target phase); members that
// cannot deliver their share are fixed at their BHP-limited rate and the rest is re-distributed (step re-solved)
struct CompGroup {
    std::string name;
    std::vector<int> wells;             // member wells (indices into the well list) under group control
    bool injector = false;
    int control = 3;                    // producer: rate control code (2 RESV, 3 oil, 4 gas, 5 water, 6 liquid)
    double target = 0;                  // producer: same units as CompWell::target (m3/s); injector: mol/s
    // GUIDERAT (ECLIPSE): GR = POT^A / (B + C R1^D + E R2^F) of the nominated phase (0 OIL, 1 LIQ, 2 GAS, 3 RES;
    // -1 = guide rates from the well potentials). OIL / RES: R1 = WOR, R2 = GOR; LIQ: R1 = WCT, R2 = GLR; GAS: R1 = WGR,
    // R2 = OGR. potScale / r1Scale / r2Scale convert SI (sm3/s, sm3/sm3) to the deck units of the formula.
    int grPhase = -1;
    double grA = 0, grB = 0, grC = 0, grD = 0, grE = 0, grF = 0;
    double potScale = 86400.0, r1Scale = 1.0, r2Scale = 1.0;
};

// gas relative permeability with hysteresis (plots / checks): drainage table `table` (0 = sat, k = satRegions[k-1]),
// imbibition table imbTable (-1 = none: Land sgtMax / drainage), scanning curve from the historical maximum sgMax
// (sg >= sgMax: drainage); trapped = Land / Killough trapped gas saturation of that scanning curve
double gasRelPerm(const CompCase& c, int table, int imbTable, double sg, double sgMax, double* trapped = nullptr);

// C18 discrete adjoint of the compositional FIM step sequence (CompSimulator::adjointGradient)
// observation: well `well`, quantity qty = 0 stock-tank oil, 1 surface gas, 2 surface water (sm3/s, production),
// 3 BHP (Pa), 4 gas injection, 5 water injection (sm3/s), 10 + i molar rate of component i (mol/s, production > 0)
struct CompAdjObs { int well = 0; int qty = 0; double t = 0, value = 0, sigma = 1; };
// parameter: apply(case, m) scales the property by the multiplier m (base m = 1); gradients are dF/dm at m = 1
struct CompAdjParam { std::string name; std::function<void(CompCase&, double)> apply; double h = 1e-4; };
struct CompAdjResult { bool ok = false; std::string error; double F = 0; std::vector<double> grad; int steps = 0, linear = 0; double seconds = 0; };
// ready-made parameters (cells: 0-based simulator cells; empty = all)
CompAdjParam compParamPermMult(const std::string& name, const std::vector<int>& cells);   // transmissibilities touching the cells
CompAdjParam compParamPoroMult(const std::string& name, const std::vector<int>& cells);   // pore volume (and initial moles)
CompAdjParam compParamWellIndexMult(const std::string& name, int well);

class CompSimulator {
public:
    // ---- C18 adjoint: record converged states from now on (call after setup, before the first advanceTo)
    void setAdjointRecording(bool on);
    int adjointSteps() const { return adjRec_.empty() ? 0 : int(adjRec_.size()) - 1; }
    std::vector<double> stepSizes() const;          // recorded dt sequence
    // replay a fixed step sequence (finite-difference checks); a step that does not converge fails advanceTo
    void setForcedSteps(const std::vector<double>& dts) { forced_ = dts; forcedIdx_ = 0; }
    // Newton tolerances (component / volume residual relative to cell moles / pore volume), relaxed acceptance, GMRES tol
    void setNewtonTolerance(double tolN, double tolV, bool allowRelaxed, double linTol) { tolN_ = tolN; tolV_ = tolV; relaxedOk_ = allowRelaxed; linTol_ = linTol; }
    double adjointMisfit(const std::vector<CompAdjObs>& obs);                 // F from the recorded states
    CompAdjResult adjointGradient(const std::vector<CompAdjObs>& obs, const std::vector<CompAdjParam>& params);

    void setup(const CompCase& c);
    // advance to time t (s) with adaptive steps; returns false if the step size collapsed
    bool advanceTo(double t, double dtMax = 86400.0 * 30);
    double time() const { return t_; }
    const CompStats& stats() const { return st_; }
    int nc() const { return nc_; }
    int ncw() const { return ncw_; }                                // nc (+1 with mobile water: index nc = water)
    int cells() const { return n_; }
    const CompCase& model() const { return c_; }
    // schedule: wells may be changed between advanceTo calls (controls, targets, open / shut, connections)
    std::vector<CompWell>& wells() { return c_.wells; }
    std::vector<CompGroup>& groups() { return groups_; }
    // surface rates of a well report (m3/s): stock-tank oil, gas, water; injected moles (mol/s, > 0)
    void surfaceRates(size_t w, const CompWellReport& X, double& qo, double& qg, double& qw, double& qinj) const;
    // state
    const std::vector<double>& pressure() const { return p_; }
    const std::vector<double>& moles() const { return N_; }        // n * ncw
    std::vector<double> gasSaturation() const;
    std::vector<double> waterSaturation() const;
    double componentInPlace(int i) const;                           // mol
    double vaporisedWater() const;                                  // mol of water in the gas phase (waterVap)
    double injected(int i) const { return cumInj_[size_t(i)]; }     // mol
    double produced(int i) const { return cumProd_[size_t(i)]; }    // mol
    double hcpvPressure() const;                                    // hydrocarbon-pore-volume weighted pressure (FPR)
    void inPlaceSurface(double& oil, double& gas, double& wat) const;   // FOIP / FGIP / FWIP (sm3) via the separator
    // per mole of hydrocarbon stream z: stock-tank oil (m3) and surface gas (sm3) from the separator train
    void surfaceVolumes(const std::vector<double>& z, double& vo, double& vg, int train = 0) const;
    // per-component yields of the separator train for feed z: stock-tank oil (m3) and surface gas (sm3) per mole of component i
    // (sum_i z_i fo_i = vo); used as lagged, composition-linear surface-rate weights (E300 practice)
    void componentYields(const std::vector<double>& z, std::vector<double>& fo, std::vector<double>& fg, int train = 0) const;
    double waterSurfaceVolumePerMole() const;                       // sm3/mol (CO2STORE: brine per mol H2O)
    // CO2STORE inventory (kg CO2): free CO2 in mobile gas, free CO2 where Sg <= critical gas saturation (krg = 0,
    // residually trapped), CO2 dissolved in brine
    void co2Inventory(double& mobile, double& immobile, double& dissolved) const;
    Co2Inventory co2Trapping() const;                               // all trapping mechanisms (C5)
    // per-cell state of the trapping / interaction models
    const std::vector<double>& sgMax() const { return sgHist_; }    // historical maximum gas saturation (hysteresis)
    std::vector<double> trappedGasSaturation() const;               // Sgt of the current scanning curve
    std::vector<double> saltSaturation() const;                     // solid salt volume / reference pore volume
    std::vector<double> mineralCO2() const;                         // CO2 fixed in minerals per cell, kg / m3 pore volume
    void saltInPlace(double& dissolved, double& solid) const;       // NaCl in place (mol): dissolved / precipitated (salt model)
    const std::vector<double>& porosityFactor() const { return pvMult_; }   // phi / phi0
    const std::vector<double>& permFactorCells() const { return kMult_; }   // k / k0
    // C25 sequential geomechanics / fault mechanics coupling (applied between report steps; moles are conserved, so a pore
    // volume change acts on the pressure of the next step): pore-volume and permeability factors per cell (empty = 1) and
    // transmissibility multipliers of the connections between cell pairs (absolute w.r.t. the base transmissibility)
    void setGeoFactors(const std::vector<double>& pvFactor, const std::vector<double>& permFactor);
    void setConnectionMultipliersByCells(const std::vector<std::pair<int, int>>& cells, const std::vector<double>& mult);
    const std::vector<double>& geoPvFactor() const { return geoPv_; }
    // liquid / vapour / overall mole fractions per cell (n * nc; absent phase = 0)
    void compositions(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z) const;
    // CO2SOL: moles of CO2 dissolved in the brine of cell q; field split of CO2 (kg) between hydrocarbon phases and brine
    double aqueousCo2Moles(int q) const;
    void co2SolSplit(double& hcKg, double& aqKg) const;
    std::vector<double> dissolvedCO2() const;                       // CO2 mole fraction of the brine per cell (CO2STORE)
    const co2::Fluid& co2Fluid() const { return cf_; }
    double gasSurfaceVolumePerMole() const;                         // ideal gas at standard conditions, sm3/mol
    const std::vector<CompWellReport>& wellReports() const { return rep_; }
    // diagnostics: residual and Jacobian at the current iterate (dense, for tests on small grids)
    void residualAndJacobian(double dt, std::vector<double>& R, std::vector<double>* Jdense);
    void setState(const std::vector<double>& p, const std::vector<double>& N) { p_ = p; N_ = N; }
    void storeOld() { pOld_ = p_; NOld_ = N_; beginStep(); }
    std::vector<double> lastK(int cell) const { return K_[size_t(cell)]; }
    int phasesIn(int cell) const { return nph_[size_t(cell)]; }
    double wellBhp(int w) const { return size_t(w) < rep_.size() ? rep_[size_t(w)].bhp : 0.0; }
    bool wellOnBhp(int w) const { return size_t(w) < rep_.size() && rep_[size_t(w)].onBhp; }
    // C14 restart: complete dynamic state (pressure / moles, K warm starts, cumulatives, well reports and lagged yields,
    // hysteresis / salt / mineral / rock-factor history, step size); loadState() follows setup() of the same model
    std::string saveState() const;
    void loadState(const std::string& blob);
    template <class Ar> void stateIO(Ar& ar);

    // implementation details (public for the assembly helpers in compsim.cpp)
    struct Conn { int a, b; double T, dD, Td, Tdisp = 0; };
    struct WellAux { std::vector<double> H; std::vector<double> fo, fg, fw; };   // head per connection (Pa); per-component separator yields
    struct AdjRec {
        double t = 0, dt = 0;
        std::vector<double> p, N, pvEff, kMult, sgHist, connT;
        std::vector<std::vector<double>> K; std::vector<int> nph;
        std::vector<WellAux> aux; std::vector<CompWell> wells;
    };
private:
    void adjRecord(double dt);
    // residual / Jacobian / well reports at a recorded step (state cur, old state prev); the model's own wells keep their
    // connection data (perturbed models), controls come from the record
    void adjEval(const AdjRec& prev, const AdjRec& cur, const CompSimulator& base, std::vector<double>& R, BSR* J,
                 std::vector<CompWellReport>* wr);
    double adjObsValue(const CompAdjObs& o, const CompWellReport& X) const;
    std::vector<AdjRec> adjRec_;
    bool adjOn_ = false;
    std::vector<double> forced_; size_t forcedIdx_ = 0;
    double tolN_ = 1e-7, tolV_ = 1e-7; bool relaxedOk_ = true;
    void beginStep();                   // lagged well quantities (wellbore density, separator factors)
    void endStep(double dt);            // hysteresis state, salt transport / precipitation, mineral inventory
    void updateRockFactors();           // porosity / permeability factors from salt + minerals (lagged one step)
    void allocateGroups(const std::vector<char>& fixed, const std::vector<double>& fixRate);
    double groupQuantity(size_t w, const CompWellReport& X, const CompGroup& G) const;
    std::vector<CompGroup> groups_;
    bool dbg_ = false;
    double linTol_ = 1e-5;
    double dpTarget_ = 5e6, dzTarget_ = 0.2;   // target changes per time step (Pa, overall mole fraction)
    bool cutLast_ = false;
    std::vector<double> sAq_;   // CO2SOL: aqueous fraction of each cell's CO2 (warm start of the split solve)
public:
    double tAsm_ = 0, tSetup_ = 0, tSolve_ = 0;   // wall time split (s): assembly (incl. flash), CPR setup, GMRES
private:
    CompCase c_;
    int n_ = 0, nc_ = 0, ncw_ = 0, nv_ = 0;
    double t_ = 0;
    eos::TParams tp_;
    co2::Fluid cf_;
    double co2Vg_ = 0, co2Vw_ = 0;      // CO2STORE surface volumes: sm3 per mol CO2 / per mol H2O (brine)
    std::vector<std::vector<eos::TParams>> sepTp_;   // per separator train (0 = field), stages (last = stock tank)
    std::vector<std::vector<std::pair<double, double>>> sepSt_;   // stage (T, p) per train
    std::vector<double> T0_;            // connection transmissibilities before rock-factor changes
    std::vector<double> geoPv_, geoK_, connX_;   // C25 external (geomechanics / fault) factors; empty = 1
    void composeFactors();
    std::vector<double> pvEff_;         // pv0 * porosity factor
    std::vector<double> sgHist_, saltDiss_, saltSolid_, minProd_, minReact_, minRate_, pvMult_, kMult_;
    std::vector<double> pv0_, depth_;
    std::vector<Conn> conns_;
    std::vector<double> p_, N_, pOld_, NOld_;
    std::vector<std::vector<double>> K_;     // last converged K per cell (empty = single phase)
    std::vector<int> nph_;
    std::vector<double> cumInj_, cumProd_;
    std::vector<WellAux> aux_;
    std::vector<CompWellReport> rep_;
    CompStats st_;
    double dtLast_ = 86400.0;
};

}  // namespace comp
}  // namespace ressim
