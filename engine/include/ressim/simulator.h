// Fully-implicit black-oil simulator (water / oil / gas with dissolved gas).
//   Primary variables per cell: Po, Sw, and Sg (free gas present) or Rs (undersaturated)
//   -> variable switching. Wells are fully coupled with BHP as an extra unknown.
#pragma once
#include "ressim/case.h"
#include "ressim/fluid.h"
#include "ressim/grid.h"
#include "ressim/linsolve.h"
#include "ressim/bsr.h"
#include <cstdlib>
#include "ressim/msw.h"
#include <vector>
#include <unordered_map>
#include <array>
#include <map>
#include <set>
#include <string>
#include <functional>
#include <atomic>
#include <cstdint>
#include <memory>

namespace ressim {
namespace hffield { class Field; struct Stats; }
namespace phasefield { struct Field; }
namespace udq { class Set; struct Action; struct Source; }
namespace pyact { class Host; struct Action; }

enum Phase { PW = 0, PO = 1, PG = 2 };

struct WellRates {
    double oil = 0, water = 0, gas = 0;        // surface production rates (>0)
    double waterInj = 0, gasInj = 0;           // surface injection rates (>0)
    double resv = 0;                           // reservoir voidage rate
    double eff = 1.0;                          // C33 WEFAC efficiency factor: rates above are 'while flowing'; group / field
                                               // rates, cumulatives and the reservoir source use eff x rate (ECLIPSE / OPM)
    double polyInj = 0, polyProd = 0;          // polymer injection / production rate (mass / day, POLYMER option)
    double solvInj = 0, solvProd = 0;          // C35 SOLVENT surface injection / production rate (gas units / day)
    double lift = 0;                           // C35 lift gas injection rate (gas units / day; WGLIR)
    double bhp = 0;
    bool stopped = false;                      // C49 well with STOP status (no surface flow): bhp is its zero-inflow BHP,
                                               // reported as WBHP like OPM / ECLIPSE (SHUT wells report 0)
    double thpNet = 0;                         // C29 THP limit imposed by the surface network (0: none)
    double dpWellbore = 0;                     // C29 friction + ICD pressure drop at the deepest connection (case units, + = more drawdown)
    double oilH = -1, waterH = -1, gasH = -1;  // observed (history) rates, -1 if none
    Control mode = Control::BHP;
    bool open = false;
    // multi-component engines (compositional / thermal): component molar rates (mol/s; production > 0, injection < 0)
    // and cumulative produced / injected moles, one entry per component of the run (RunComponents order)
    std::vector<double> comp, compCumP, compCumI;
    // C43 well performance events of the time step (OPM WPWE0-7): 0 drilled (always 0: no drilling queue), 1 connections
    // opened, 2 connections closed (not those swept by a '+CON' reach), 3 closed to bottom, 4 stopped / 7 shut by the
    // simulator, 5 injector -> producer, 6 producer -> injector
    int pwe[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // C46 productivity / injectivity index per phase (surface volume per pressure unit; OPM WPIO / WPIW / WPIG; WPIL = oil +
    // water, WPI = preferred phase); producers: oil incl. Rv x gas, gas incl. Rs x oil; injectors: the injected phase
    double pi[3] = {0, 0, 0};
    double piPref = 0;
};

struct SummaryPoint {
    double time = 0, dt = 0;
    int newton = 0, linear = 0;
    double fopr = 0, fwpr = 0, fgpr = 0, fwir = 0, fgir = 0;
    double fopt = 0, fwpt = 0, fgpt = 0, fwit = 0, fgit = 0;
    double fpr = 0;                            // average pressure, hydrocarbon-pore-volume weighted (ECLIPSE FPR)
    double fprp = 0;                           // pore-volume weighted average pressure (ECLIPSE FPRP)
    double fwct = 0, fgor = 0;
    double oip = 0, wip = 0, gip = 0;          // fluids in place (surface units)
    double faqr = 0, faqt = 0;                 // aquifer influx rate / cumulative (reservoir volumes)
    double faqrS = 0, faqtS = 0;               // same in surface water volumes
    double fcir = 0, fcpr = 0, fcit = 0, fcpt = 0;   // polymer injection / production rate and totals (mass)
    double fcip = 0, fcad = 0;                 // polymer in solution / adsorbed on the rock (mass)
    double fnir = 0, fnpr = 0, fnit = 0, fnpt = 0, fnip = 0;   // C35 SOLVENT injection / production rates, totals, in place (gas units)
    double fglir = 0, fglit = 0;               // C35 lift gas injection rate and total (FGLIR, FGLIT)
    double fgcr = 0, fgimr = 0, fgsr = 0;      // C35 gas consumption (GCONSUMP), import, sales = produced - injected - consumed + imported
    double fgct = 0, fgimt = 0, fgst = 0;      //     and their totals (FGCT, FGIMT, FGST); satellite (GSATPROD) rates are in FOPR / FWPR / FGPR
    double co2Mob = 0, co2Trap = 0, co2Diss = 0;   // CO2STORE: CO2 mass (kg) mobile / residually trapped / dissolved in brine
    double co2Struct = 0, co2Min = 0, co2KrImmob = 0, saltSolid = 0;   // structurally held part of co2Mob / mineralised / krg = 0 / NaCl solid
    double co2Oil = 0;                         // tracked component dissolved in the oil phase (thermal / compositional), kg
    // thermal (STARS-like core, deck units: degC / degF, kJ / Btu): PV-average / max temperature, steam PV fraction,
    // cumulative heat loss to over / underburden, net injected enthalpy, cumulative steam-oil ratio (CWE)
    double favgT = 0, fmaxT = 0, fsteamPV = 0, fhlt = 0, fenergyIn = 0, csor = 0;
    double freactRate = 0, freactCum = 0, fuplift = 0, fupliftRes = 0, fmaxFail = 0;   // THMC: reaction heat (W, J), heave (m), failure index
    // multi-component engines: field component production / injection rates (mol/s, both >= 0) and cumulatives (mol)
    std::vector<double> fcompP, fcompI, fcompPT, fcompIT;
    std::vector<WellRates> wells;
    // C31 user-defined quantities after this step: FU* values and WU*:WELL / GU*:GROUP entries (NaN = undefined)
    std::vector<std::pair<std::string, double>> udq;
};

struct Snapshot {                              // per active cell, at a report time
    double time = 0;
    std::vector<float> pressure, sw, sg, so, rs;
    std::vector<float> rv;                      // vaporised oil (only with VAPOIL)
    std::vector<float> poly, polyAds;           // polymer concentration / adsorbed polymer (only with POLYMER)
    std::vector<float> xco2;                    // CO2 mole fraction in the brine (CO2STORE) / of the cell (compositional: z of component 1)
    // further named arrays (compositional / CO2STORE): ZMFk / XMFk / YMFk (k = 1..nc), SGTRAP, SGMAX, SSALT, MINCO2,
    // PORFAC, PERMFAC
    std::vector<std::pair<std::string, std::vector<float>>> extra;
    const std::vector<float>& field(const std::string& n) const {
        if (n == "PRESSURE") return pressure;
        if (n == "SWAT") return sw;
        if (n == "SOIL") return so;
        if (n == "SGAS") return sg;
        if (n == "RS") return rs;
        if (n == "RV") return rv;
        if (n == "POLYMER") return poly;
        if (n == "PADS") return polyAds;
        if (n == "XCO2") return xco2;
        for (const auto& e : extra) if (e.first == n) return e.second;
        static const std::vector<float> none;
        return none;
    }
};

struct InPlace { double oil = 0, water = 0, gasFree = 0, gasDissolved = 0, oilVaporized = 0, gasAdsorbed = 0; };   // gasAdsorbed: Langmuir (C15)   // oil = liquid oil

struct AquiferState { int id = 0; bool fetkovich = true; int cells = 0; double pressure = 0, rate = 0, cumulative = 0; };

struct RunStats { int timesteps = 0, cuts = 0, newton = 0, linear = 0; double seconds = 0; bool cpr = false; int amgLevels = 0; int polyOuter = 0; int relaxedIts = 0;
                  double secAssemble = 0, secPrecond = 0, secSolve = 0, secUpdate = 0;   // C12 wall-clock profile of the Newton loop
                  int gpuSolves = 0, gpuFallbacks = 0;                                     // C16 device linear solves
                  double gpuWallSec = 0, gpuDeviceSec = -1;                                // C18 device solve wall / GPU-clock busy time (-1 n/a)
                  int ddSolves = 0, ddFallbacks = 0; double ddSec = 0;                     // C18 MPI domain-decomposed linear solves
                  double lastDP = 0, lastDS = 0;                                           // C27 largest p / S change of the last step
                  int econShut = 0, econTests = 0; };                                      // C27 WECON shut-ins / WTEST re-openings

// C27: '/'-terminated records of a verbatim deck keyword body -> tokens (quotes stripped, n* expanded, comments removed)
std::vector<std::vector<std::string>> verbatimRecords(const std::string& body);
double verbatimNum(const std::vector<std::string>& r, size_t i, double def);              // item i of a verbatim record (number)
std::string verbatimStr(const std::vector<std::string>& r, size_t i, const std::string& def);   // item i (upper case string)

// C28 (C21-3) distributed assembly for MPI: the state is replicated on every rank (all ranks run the same Newton loop in
// lockstep), but each rank evaluates the cell properties of / assembles only its own contiguous range of cell rows (plus
// halo and well-completion cells), the convergence norms are reduced over the ranks, the linear system is solved by the
// collective domain-decomposed solver straight from the owned rows (no gather / scatter of the matrix) and the Newton
// update is all-gathered. Every decision (convergence, chops, well switching) is taken on identical global numbers.
struct DistAssembly {
    int a0 = 0, a1 = 0;                                   // owned cells [a0, a1)
    std::vector<int> starts;                              // cell partition of all ranks (parts + 1 offsets)
    std::function<void(double*, int)> sumN, maxN;         // in-place global reductions
    // collective solve of the Schur-reduced block system (only the owned rows of S / b are valid); x gets ALL entries
    std::function<bool(const BSR&, const std::vector<double>&, std::vector<double>&, double, int, int&)> solve;
};

class Simulator {
public:
    // Throws std::runtime_error with a readable message on invalid input.
    void setup(const Case& c);
    // C28: distributed assembly (after setup, before initialize); cellPartition = nnz-balanced contiguous cell ranges
    void setDistributedAssembly(const DistAssembly& d) { da_ = d; daOn_ = d.solve != nullptr; }
    bool distributedAssembly() const { return daOn_; }
    std::vector<int> cellPartition(int parts) const;
    void initialize();                          // hydrostatic equilibrium (EQUIL)
    bool runReportStep();                       // false when the schedule is finished or cancelled
    bool finished() const { return step_ >= int(case_.schedule.size()) || econStop_; }
    int  reportStep() const { return step_; }
    double time() const { return time_; }

    const Case& caseData() const { return case_; }
    const Grid& grid() const { return grid_; }
    const FluidModel& fluid(int pvtRegion = 0) const { return fluids_[size_t(pvtRegion)]; }
    int numPvtRegions() const { return int(fluids_.size()); }
    const std::vector<SummaryPoint>& summary() const { return summary_; }
    // well connections after setup (open completions): global cell (0-based i,j,k), well index (CF), depth
    struct ConnInfo { std::string well; int i = 0, j = 0, k = 0; double wi = 0, depth = 0; };
    std::vector<ConnInfo> connections() const;
    // per open connection at the current state: surface rates (+ = production, injection < 0), block / connection
    // pressure (ECLIPSE PRT well-connection report)
    struct ConnRate { std::string well; int i = 0, j = 0, k = 0; double oil = 0, water = 0, gas = 0, resv = 0, pBlock = 0, pConn = 0, wi = 0, sw = 0, sg = 0; };
    std::vector<ConnRate> connectionRates() const;
    // constant initial Rs used per EQLNUM region by the equilibration (-1: RSVD table / not live / region empty)
    const std::vector<double>& initialRsByRegion() const { return initRs_; }
    const std::vector<double>& initialRvByRegion() const { return initRv_; }
    const std::vector<Snapshot>& snapshots() const { return snaps_; }
    InPlace inPlace() const;
    const RunStats& stats() const { return stats_; }
    double totalSimTime() const;
    std::vector<AquiferState> aquiferStates() const;
    // per active cell: total phase flux through each connection at the last converged state
    // (reservoir volume / day, positive a->b), used for streamlines
    void connectionFluxes(std::vector<double>& qw, std::vector<double>& qo, std::vector<double>& qg) const;
    // per well: active cell, signed reservoir rate (+ out of the reservoir) at the last state
    struct WellCellRate { int cell; double rate; int well; };
    std::vector<WellCellRate> wellCellRates() const;
    const std::vector<std::string>& wellNames() const { return wellNames_; }
    // C27 passive tracers (ECLIPSE TRACER / WTRACER + RESSIM TRACERPROP): names, phase (PW / PO / PG), per-cell concentration
    // (tracer amount per surface volume of the phase) - the snapshots also carry them as TR_<name>
    struct TracerInfo { std::string name; int phase = 0; double halfLife = 0, kd = 0, kPart = 0, alphaL = 0, dm = 0; int partner = -1, parent = -1; };
    const std::vector<TracerInfo>& tracers() const { return trInfo_; }
    const std::vector<double>& tracerConc() const { return trC_; }   // [tracer * numActive + cell]
    double tracerInPlace(size_t t) const;      // free-phase + adsorbed tracer amount at the current state
    double tracerDecayed(size_t t) const { return t < trDecayCum_.size() ? trDecayCum_[t] : 0.0; }
    // restore the state of a saved report step (pressure, saturations, well BHPs) for post-processing
    void restoreReport(int snapIndex);
    // set the state from an external snapshot (e.g. an ECLIPSE / OPM restart step mapped to this grid's active cells)
    // with well BHPs in wellNames() order (<= 0: shut) - wells are put on BHP control; used for streamlines of
    // results computed by another simulator. Throws if the snapshot size does not match.
    void restoreExternal(const Snapshot& sn, int reportStep, const std::vector<double>& bhp);

    // ---- discrete adjoint / forward sensitivities of an observation misfit (history-matching gradients, C9)
    // Enable recording BEFORE initialize(); after the forward run adjointGradient() replays the recorded converged
    // states: backward pass  J_n^T lam_n = -dF/dx_n - (dR_{n+1}/dx_n)^T lam_{n+1}  (one transposed linear solve per
    // time step), gradient dF/dtheta = dF/dtheta|explicit + sum_n lam_n^T dR_n/dtheta, where dR_n/dtheta comes from the
    // residual of two perturbed models evaluated at the SAME recorded states (central difference, no extra simulation).
    // forwardMode = true instead solves J_n s_n = -dR_n/dtheta - dR_n/dx_{n-1} s_{n-1} per parameter and returns the
    // full sensitivity matrix dy/dtheta (Gauss-Newton / LM) - plus the gradient for cross-checking.
    // The explicit well-head update d_n = H(x_{n-1}, bhp, d_{n-1}) is differentiated as an augmented lagged state (own
    // multipliers); other explicit quantities - group targets, THP limits, RESV conversion, hysteresis history, DRSDT
    // limits, aquifer state - are frozen at their recorded values. Control / phase switching events make the true
    // misfit piecewise smooth: the gradient is exact on the current branch.
    struct AdjObs { ObsQty qty = ObsQty::FOPR; int well = -1; double time = 0, value = 0, sigma = 1; };
    struct AdjParam { const Case* plus = nullptr; const Case* minus = nullptr; double h = 1e-4; };   // theta +- h
    // field parameters with analytic dR/dtheta (C11): one reverse pass gives dF/d ln(theta) for EVERY cell / face /
    // connection - cost independent of the parameter count (no perturbed replays)
    struct AdjFieldReq { bool perm = false, poro = false, trans = false, wi = false; };
    struct AdjResult {
        bool ok = false; std::string error;
        double F = 0;                                  // sum ((y - d) / sigma)^2
        std::vector<double> y;                         // simulated value per observation
        std::vector<double> grad;                      // dF/dtheta per parameter
        std::vector<std::vector<double>> dydp;         // forward mode: dy_i / dtheta_j  [obs][param]
        int steps = 0, linearIts = 0;
        double maxLinRes = 0, replayDev = 0, seconds = 0;   // worst linear residual; replay vs forward rate deviation
        // field gradients (AdjFieldReq): dF/d ln kx, ln ky, ln kz, ln phi per active cell; dF/d ln T per grid connection;
        // dF/d ln WI per well and connection
        std::vector<double> gLogKx, gLogKy, gLogKz, gLogPoro, gLogTrans;
        std::vector<std::vector<double>> gLogWI;
        double secField = 0;
    };
    void setAdjointRecording(bool on) { adjRec_ = on; }
    bool adjointRecording() const { return adjRec_; }
    AdjResult adjointGradient(const std::vector<AdjObs>& obs, const std::vector<AdjParam>& params, bool forwardMode = false,
                              const AdjFieldReq* field = nullptr);

    // ---- C14 restart: complete dynamic state at a report-step boundary (primary variables, well / aquifer / hysteresis /
    // VAPPARS / polymer / DRSDT history, cumulatives, time-step size). loadState() follows setup() + initialize() of the
    // SAME grid and well list (schedule items after the restart step may differ, as in ECLIPSE RESTART); the summary /
    // snapshot lists start empty and the run continues with runReportStep(). Throws on a mismatch.
    std::string saveState() const;
    void loadState(const std::string& blob);
    // ---- C18 in-run hydraulic-fracture growth (HFRAC with HFGROW / kic > 0)
    struct HfEvent { double t = 0; int frac = 0, side = 0; double length = 0, KI = 0, pnet = 0; };   // side -1 / +1 wing, length = wing length (case units), KI Pa m^0.5, pnet case pressure
    const std::vector<HfEvent>& hfGrowthEvents() const { return hfEvents_; }
    // C28 in-run phase-field growth (HFGROW ... PHASEFIELD): one record per report step and phase-field fracture
    struct PfRecord { double t = 0; int frac = 0; double volume = 0, pressure = 0, halfLength = 0;   // m2 per m height, Pa, m
                      double pFlow = 0, sigma = 0; int stagger = 0; bool active = true; std::string note;
                      std::shared_ptr<const phasefield::Field> field; };   // case units; C30: field of the step (local frame)
    const std::vector<PfRecord>& phaseFieldRecords() const { return pfRec_; }
    // C30 fracture field output (hf_field.h): statistics (nullptr when inactive) and a probe of the stress (Pa, tension +)
    // / displacement (m) change at a point (m, z up) for the openings of the last report
    const hffield::Stats* hfFieldStats() const;
    bool hfFieldAt(const double x[3], double S[6], double u[3]) const;
    // C30 phase-field frame of fracture f: origin (injection point) and unit strike in plan (case length units)
    bool phaseFieldFrame(int frac, double& ox, double& oy, double& tx, double& ty) const;
    // C31 UDQ / UDA / ACTIONX (udq.h): the quantity set (nullptr without UDQ), fired actions
    const udq::Set* udqSet() const { return udq_.get(); }
    struct ActionEvent { double t = 0; int step = 0; std::string name, wells, notes; };
    const std::vector<ActionEvent>& actionEvents() const { return actEvents_; }
    int numActions() const;
    double wellRefDepth(int w) const { return w >= 0 && w < int(wells_.size()) ? wells_[size_t(w)].refDepth : 0.0; }   // C48: BHP reference depth
    double groupRateTarget(const std::string& g, const char* q) const;   // C48: current GCONPROD rate (ORAT WRAT GRAT LRAT RESV; -1 none)
    // C37 PYACTION (pyaction.h): number of Python actions, interpreter version ("" until started), last error
    int numPyActions() const;
    std::string pythonVersion() const;
    // C33 scripting (Python API / api::Session): applies SCHEDULE keywords (deck text) at the current report-step
    // boundary like an ACTIONX body (UDQ records extend the quantity set, COMPDAT may add connections); recorded for the
    // restart replay. Returns the importer messages (keywords that cannot act at run time).
    std::vector<std::string> applySchedule(const std::string& scheduleText);
    // C32 multisegment wells (msw.h): definitions in force, and the segment state per accepted time step
    const std::vector<msw::Well>& mswWells() const { return msw_; }
    struct MswRecord { double t = 0; int msw = 0; std::vector<double> p, q, rho, dFric, dDev, dAcc, holdup, noSlip;   // per segment (Well::segs order), case units; C36 acceleration / gas holdup
                       std::vector<double> oilHold, dPump, flimFac, rhoHead; };   // C37 oil fraction of the liquid, pump boost, limiter factor, head density
    const std::vector<MswRecord>& mswRecords() const { return mswRec_; }
    int mswOfWell(int w) const { return w >= 0 && size_t(w) < mswOf_.size() ? mswOf_[size_t(w)] : -1; }
    // C45 producer cross-flow in multisegment wells: MSWOPTS third token or RESSIM_MSW_XFLOW=1 (experimental, default off)
    bool mswXflow() const { static const bool env = std::getenv("RESSIM_MSW_XFLOW") != nullptr; return case_.mswOpt.xflow || env; }
    // C49 converged (Picard-coupled) multisegment connection heads; RESSIM_MSW_HEADS=0/1 overrides the deck option
    bool mswConvergedHeads() const {
        static const char* env = std::getenv("RESSIM_MSW_HEADS");
        if (env && *env) return *env == '1';
        return case_.mswOpt.heads != 0;   // C49 A/B on 24 OPM MSW cases: never worse (TEST_WSEGSICD 93.4 -> 34.0 %), so the default
    }
    bool headsMswOnly_ = false;               // C49 computeHeads(): update multisegment wells only (Picard head iterations)
    int mswPicardIts_ = 0, mswPicardSteps_ = 0, mswPicardUnconv_ = 0;   // C49 statistics
    // C29 surface network (NODEPROP / BRANPROP extended network, balanced explicitly at the start of every time step)
    struct NetNode {
        std::string name, parent; double p = 0, qo = 0, qw = 0, qg = 0; int vfp = 0; bool fixed = false;
        // C35 facilities: branch ALQ (BRANPROP / GRUPNET item 4; alqEq 1 / 2 = surface oil / gas density), no-loss branch
        // (VFP 9999), lift gas added to the node's flow, automatic choke (NODEPROP AS_CHOKE, CHOKE_GROUP), NEFAC
        double alq = 0; int alqEq = 0; bool noLoss = false, addLift = false, choke = false; std::string chokeGroup; double nefac = 1.0;
        double dpChoke = 0;   // current pressure drop of an automatic choke (node pressure - uptree pressure)
        // C36 automatic compressor on the branch to the parent (NETCOMPA): VFP table / ALQ when on, gas consumption taken
        // at the extraction group's node, switched on (one sequence number per time step) while the group misses its
        // production target in 'compPhase'; multi-level: ALQ / consumption linear in the level; COMPOFF switches all off
        int compVfp = -1, compLevels = 1, compLevel = 0, compSeq = 0;   // compVfp -1: no compressor, 0: keep the branch table
        double compAlq = 0, compAlq1 = 0, compGas = 0;
        std::string compGroup, compPhase = "GAS", compExtract;
    };
    const std::vector<NetNode>& networkNodes() const { return netNodes_; }
    struct NetRecord { double t = 0; std::vector<double> p; };   // node pressures (networkNodes order) per accepted step
    const std::vector<NetRecord>& networkHistory() const { return netHist_; }
    // C29 salt precipitation (SALT): solid fraction of the pore volume and permeability factor per active cell
    const std::vector<double>& saltSolid() const { return saltS_; }
    const std::vector<double>& saltPermFactor() const { return saltK_; }
    double saltSolidMass() const;             // kg of precipitated NaCl
    double saltEvaporated() const { return saltEvapCum_; }   // cumulative evaporated water (surface volume)
    // open the pieces of both wings of fracture 'frac' out to the wing length L (case length units; growth only)
    bool setHfWingLength(int frac, double L);
    double hfWingLength(int frac, int side) const;     // current open wing length (case length units)
    bool hfCellOpen(int hfCell) const { return hfOpen_.empty() || hfOpen_[size_t(hfCell)] != 0; }
    // C21 G1-c: per hydraulic-fracture piece conductance multiplier from contact mechanics (cubic law (a / a_ref)^3);
    // fracture-fracture connections take the harmonic mean of their two pieces; empty vector = none
    void setHfConductance(const std::vector<double>& mult);
    // C24 two-way fault coupling: transmissibility of the given connections = base value x mult (base = value at the first
    // call, i.e. after MULTFLT / MULTX ...); repeated calls replace the previous multipliers of those connections
    void setConnectionMultipliers(const std::vector<int>& conns, const std::vector<double>& mult);
    // C25 sequential geomechanics: pore-volume multipliers per active cell w.r.t. the initial reference pore volume (empty =
    // reset; mass is conserved: the old-time accumulation keeps its mass, the pressure responds in the next step) and the
    // smallest principal total stress per active cell (case pressure units) used by HFGROW ... FESTRESS
    void setPoreVolumeMultipliers(const std::vector<double>& mult);
    void setMinStress(const std::vector<double>& sMin) { sMin_ = sMin; }
    const std::vector<double>& minStress() const { return sMin_; }
    const std::vector<double>& hfConductance() const { return hfCond_; }
    double hfCellPressure(int hfCell) const;            // current pressure of a fracture piece (case units)
    double hfCellOpenFraction(int hfCell) const {        // C21: 0..1 (sub-grid tips), else 0 / 1
        if (!hfFrac_.empty()) return hfFrac_[size_t(hfCell)];
        return hfCellOpen(hfCell) ? 1.0 : 0.0;
    }
    std::function<void(const std::string&)> log;              // optional message sink
    std::function<void(const SummaryPoint&)> onTimestep;      // called after each accepted step
    std::atomic<bool>* cancel = nullptr;

private:
    bool dumped_ = false, dumpedFail_ = false;
    // C18 HF growth state
    std::vector<uint8_t> hfOpen_;                 // per HF cell (empty: no growth)
    std::vector<double> hfS_, hfHalf_;            // signed distance of the piece centre from the injection point, half piece length
    std::vector<std::vector<int>> hfConnsOf_;     // per HF cell: its fracture-fracture connections (indices into grid_.conns)
    std::vector<double> hfT0_;                    // original transmissibility per connection (only connections of HF cells used)
    std::vector<HfEvent> hfEvents_;
    struct PfLink;                                // C28: persistent phase-field model per PHASEFIELD fracture (worker thread)
    std::vector<std::shared_ptr<PfLink>> pfLinks_;
    std::vector<PfRecord> pfRec_;
    std::shared_ptr<hffield::Field> hff_;         // C30 fracture field output
    void pfFieldArrays(Snapshot& s) const;      // C30 phase-field damage / stress / opening on the grid
    void pfReportStep();
    // C21 sub-grid tip (HFGROW ... SUBGRID): per wing (2 f + (side > 0)) continuous length and its maximum, per piece the
    // open fraction (non-sub-grid pieces mirror hfOpen_)
    std::vector<double> hfTip_, hfTipMax_, hfFrac_;
    bool hfSubLeak_ = false;
    std::vector<double> hfPrest_, hfLeakM_, hfSg0_;   // + C48 H2e: host Sg at first exposure (LEAKC)
    double cellCt(int c, double sg, double p) const;   // C48 H2e: total compressibility of cell c at saturation sg, pressure p
              // C28 H2b sub-grid leak-off: host pressure at first exposure, leak-off multiplier
    std::unordered_map<int, double> connT0_;      // C24: base transmissibility of multiplied connections
    std::vector<double> hfCond_;                  // C21 G1-c conductance multipliers per piece
    std::vector<double> pvRef0_, sMin_;           // C25 initial reference pore volumes (set on first use), FE minimum stress
    double hfPieceSigma(int f, size_t k, double z) const;
    struct HfFF { int c, a, b; };
    std::vector<HfFF> hfFF_;                      // fracture-fracture connections (index, piece a, piece b)
    void buildHfConns();
    std::vector<std::vector<int>> hfHostConnsOf_;   // C21: per HF cell its host-fracture connections (sub-grid leak-off area)
    bool hfSubgrid_ = false;
    int stateVer_ = 3;
    void setHfTips();                             // hfTip_ -> hfFrac_ / hfOpen_ of the sub-grid fractures + transmissibilities
    bool hfTipKI(int frac, int side, double& KI, double& pnet) const;   // K_I (Pa m^0.5) at the current tip
    void setupHfGrowth();
    void applyHfTrans();
    bool growFractures();
    // ---- adjoint recording / replay (adjoint.cpp)
    struct AdjWell {
        double bhp = 0; Control mode = Control::BHP; WellType type = WellType::Producer; bool open = false;
        WellControl ctrl, own; std::vector<double> dpConn; double headRho = 0, resvC[3] = {0, 0, 0}; bool resvOk = false;
        double oilH = -1, waterH = -1, gasH = -1;
    };
    struct AdjRec {
        int report = 0; double t0 = 0, dt = 0;
        std::vector<double> p, sw, x; std::vector<uint8_t> hg;
        std::vector<AdjWell> w;
        std::vector<double> aqPa, aqWe;
        std::vector<double> rsLim, sgHy, soHy, dG, dO, soMax, fRv, fRs;
        std::vector<double> hbhp; std::vector<std::vector<double>> dOld;   // inputs of the explicit well-head update
        std::vector<uint8_t> hfo;   // C20: open hydraulic-fracture pieces during this step (HFGROW frozen-topology adjoint)
        std::vector<double> hff;    // C21: open fractions (sub-grid tips) during this step
    };
    std::vector<std::pair<double, std::vector<double>>> adjHeadIn_;
    int headsOnly_ = -1;                      // computeHeads(): update only this well (adjoint finite differences)
    void adjHeads(const AdjRec& prev, const AdjRec& cur);   // recompute the well heads of step 'cur' at its inputs
    struct AdjConnD { int well = 0, k = 0, cell = 0; double dc[3] = {0, 0, 0}; double dw = 0;    // dR/d(dpConn_k)
                      double vc[3] = {0, 0, 0}; double vw = 0; };   // the connection's own terms (cell rows unscaled, well row) = dR/d ln WI_k
    std::vector<double> adjFace_;             // per grid connection and equation: flux a -> b (unscaled) = dR_a/d ln T
    std::vector<AdjConnD> adjCd_;
    bool adjRec_ = false;
    std::vector<AdjRec> adj_;                 // [0] = initial state, [n] = converged state of accepted time step n
    void adjCapture(double t0, double dt);
    void adjSetState(const AdjRec& r);
    void adjAssembleStep(const AdjRec& prev, const AdjRec& cur);
    SummaryPoint adjPoint() const;            // rates / FPR at the current state (after assemble)
    bool adjCollect_ = false;
    struct AdjRateD { int well = 0, cell = 0, k = 0; double d[5][4]; double v[5] = {0, 0, 0, 0, 0}; };   // d(oil, water, gas, winj, ginj)/d(cell vars, bhp)
    std::vector<AdjRateD> adjRd_;
    bool quasiFallback_ = false;
    mutable double cnvPh_[3] = {0, 0, 0};   // per-phase CNV of the last converged() call (oscillation detection)     // true-IMPES CPR failed once: quasi-IMPES for the rest of the run
    void restoreCells(const Snapshot& sn);
    struct CellFlux {                       // compact per-cell properties used by flux terms
        AD3 P[3], mob[3], rho[3], mobRs, mobRv, acc[3];
    };
    struct CellFull : CellFlux {
        AD3 mobS, rhoS;                     // C35 SOLVENT: solvent mobility (bs kr / mu) and density (0 without SOLVENT)
        AD3 b[3], S[3], rs, rv, ss, bs;
    };
    struct WellRT {                         // runtime well data
        int input = -1;
        WellType type = WellType::Producer;
        bool lockBhp = false;               // (reset when the next report step applies its controls)
        int stallIts = 0;                   // consecutive Newton iterations with a stalled rate-target residual               // moved to BHP by the in-Newton rule: no switch back within this time step
        std::vector<int> cell;              // active cell per open completion
        std::vector<double> wi, cdepth;
        std::vector<double> wiBase;         // connection factor when the connection is open
        std::vector<std::array<double, 3>> dlwi;   // d ln WI / d ln (kx, ky, kz) per connection (0 for a given CF / KH)
        std::vector<int> cFrom, cUntil;     // report-step window of each connection (late COMPDAT / COMPDAT SHUT)
        std::vector<int> cDef;              // C46: step the connection exists from (defined SHUT: before cFrom) - WELPI scope
        std::vector<double> cMult, cWelpi;  // C46: Completion::mult baked into wiBase; cumulative run-time WELPI factor
        std::vector<int> cComp;             // C31: completion (WellInput::completions index) of each connection
        int nActive = 0;                    // connections open in the current report step
        bool xflow = true;                  // C41: WELSPECS item 10 (crossflow allowed; ECLIPSE / OPM default YES)
        double xmix[3] = {0, 0, 0};         // C41: wellbore mixture (surface oil, water, gas fractions) of the last assembly
        double xrec[3] = {0, 0, 0};         // C44 injector crossflow: back-produced surface oil / water / gas of the last assembly,
        double xinjTot = 0;                 //   re-injected through the injecting connections (share of xinjTot), not produced
        std::vector<double> xmixConn;       // C45 MSW producer: wellbore mixture at each connection's segment (3 per connection:
                                            //   surface oil, water, gas fractions of the gross inflow of the segment's subtree)
        std::vector<double> xfRes;          // C45 MSW producer: cross-flow reservoir rate per connection at the head update
                                            //   (transient, recomputed with the heads; sweeps the tubing storage downward)
        std::vector<int> cpos;              // position of completion cell in well row cell list
    std::vector<double> qwConn;         // surface water rate per connection at the last assembly (+ out of the reservoir)
    std::vector<double> rw, re, plen;   // wellbore radius, area-equivalent block radius, perforation length per connection
    std::vector<double> wmulConn, zwConn;   // polymer: connection water mobility multiplier and shear factor in it
        std::vector<int> cellWellIdx;       // index of this well in the cell's well list
        std::vector<int> rowCells;          // unique sorted cells for the well row
        double refDepth = 0, headRho = 0;
        std::vector<double> dpConn;         // wellbore pressure difference bhp -> connection (OPM/ECLIPSE segregated densities)
        std::vector<double> dpFric;         // C29 WELLBORE: lagged friction + ICD part of dpConn (case units)
        std::vector<double> dpF;            // C29 WELLBORE: relaxed friction part (accumulated from the datum, magnitude)
        double thpNet = -1;                 // C29 network THP limit (case units, < 0 none)
        int bore = -1;                      // C29 index into Case::wellbores (-1: none)
        double pConn(size_t k, double g) const { return bhp + (k < dpConn.size() ? dpConn[k] : headRho * g * (cdepth[k] - refDepth)); }
        double bhp = 0;
        // reservoir-voidage conversion of the surface rates (ECLIPSE / OPM: field hydrocarbon-pore-volume averaged
        // pressure, Rs, Rv at the start of the time step): resv = cW qw + cO qo + cG qg (qo / qg incl. dissolved
        // / vaporised). resvOk = false -> in-situ voidage at connection conditions.
        double resvC[3] = {0, 0, 0};        // [W, O, G]
        bool resvOk = false;
        bool resvHist = false;              // WCONHIST RESV: target derived from the observed rates each time step
        WellControl own;                    // the well's own control record (before group targets / THP limits)
        double bhpDeck = -1;                // BHP limit of the record (the THP limit may raise / lower ctrl.bhp)
        bool grpCtl = false;                // a group target is applied in this time step
        WellControl ctrl;                   // control record in force
        Control mode = Control::BHP;        // active mode (after switching)
        bool open = false;
        WellRates rates;
        double eff = 1.0;                   // C33 WEFAC efficiency factor in force
        double effBase = 1.0;               // C37 WEFAC factor before the WCYCLE start-up ramp
        bool cycSched = false;              // C37 WCYCLE: the on-period started with a schedule opening (no start-up ramp)
        double cycOpen = -1e300, cycClose = -1e300;   // C37 WCYCLE: time the well last opened / was shut by its cycle (-1e300 = none)
        double guide = -1, guideT = -1e30;  // C34 GUIDERAT: last guide rate (nominated phase, damped) and when it was computed
        // C27 economic limits (WECON) / well testing (WTEST)
        int econ = 0;                       // 0 = producing normally, 1 = shut by an economic limit, 2 = re-opened by WTEST
        double econNextTest = -1;           // time of the next WTEST re-opening (-1 = none)
        int econTests = 0;                  // WTEST re-openings so far
        // C43 WPWE tracker: baseline of the last accepted time step (open connections by cell, well type) and the
        // simulator decisions of the current step ('+CON' swept cells, reach past the offender, 4 stop / 7 shut)
        bool pweBase = false, pweInj = false, pweReach = false;
        std::vector<int> pweOpenCells;
        std::set<int> pweSwept;
        int pweSim = 0;
        bool autoStop = false;              // C44 WELSPECS item 9 = STOP: simulator closures stop (WPWE4) instead of shut (WPWE7)
        std::map<int, int> econConns;       // C43: cells closed by economic connection workovers -> their cUntil before
        int pweLast[8] = {0, 0, 0, 0, 0, 0, 0, 0};   // events of the last accepted step (UDQ / ACTIONX)
        double alq = 0;                     // C35 artificial lift quantity in force (VFPPROD ALQ; lift gas rate when alqLift)
        bool alqLift = false;               // C35 alq is a lift gas rate from the optimisation / WLIFTOPT (converted for IGLR / TGLR tables)
    };
    struct EconLimit {                      // WECON record in force for a well
        bool on = false;
        double minOil = 0, minGas = 0, maxWct = 0, maxGor = 0, maxWgr = 0, minLiq = 0;
        std::string action = "NONE";        // workover for ratio limits: NONE / WELL / CON / +CON / PLUG (C34: connection workovers, COMPLUMP)
        std::string action2 = "NONE";       // C34 workover for the secondary water-cut limit (WECON items 11 / 12)
        double maxWct2 = 0;
        bool potn = false;                  // C34 WECON item 10 = POTN: rate limits on the well potentials
        bool endRun = false;
    };
    struct WellTest { double interval = 0; std::string reasons; int maxTests = 0; };
    std::vector<EconLimit> econLim_;        // per runtime well
    struct ConnEcon {                       // C41 CECON limits of one connection (0 = no limit)
        bool on = false, checkStopped = false;
        double maxWct = 0, maxGor = 0, maxWgr = 0;
        double minOil = 0, minGas = 0;      // C42 items 11 / 12: minimum completion surface rates (not in OPM)
        std::string workover = "CON";       // CON / +CON / WELL / PLUG (C42: depth-ordered plug-back, see ceconCheck)
    };
    std::vector<std::vector<ConnEcon>> cecon_;   // per runtime well, per connection
    void ceconCheck();                      // C41 connection economic limits (OPM WellTest::updateWellTestStateCECON)
    void pweBaseline(WellRT& w);            // C43 WPWE: snapshot of the current connection / type state
    void pweTrack();                        // C43 WPWE: events of the accepted step -> last summary point, new baseline
    bool pweExists(const WellRT& w) const;  // C44 WPWE: the well has entered the schedule
    std::vector<WellTest> wtest_;           // per runtime well (interval <= 0 = none)
    bool econAny_ = false, econStop_ = false;
    void econControls(int step);            // WECON / WTEST records of this report step (kept as extra keywords)
    void econCheck();                       // end of every converged time step: limits -> shut-in / end of run
    void econTest();                        // start of every time step: WTEST re-openings

    struct AquiferRT {
        int id = 0;
        AquiferType type = AquiferType::Fetkovich;
        std::vector<int> cell;
        std::vector<double> alpha, cdepth;
        double datum = 0, p0 = 0, pa = 0;   // datum, initial and current (Fetkovich) pressure
        double ctV = 0, J = 0;              // Fetkovich
        double Tc = 1, beta = 0;            // Carter-Tracy time constant (days) and influx constant
        double We = 0, rate = 0;            // cumulative influx and last rate (reservoir volumes)
        double WeS = 0, rateS = 0;          // surface volumes
        bool autoP0 = true;
        double muw = 0.5;                   // water viscosity (outflow mobility)
        std::vector<double> tdTab, pdTab;   // C34 AQUTAB influence table (empty: built-in)
    };
    void log_(const std::string& s) const { if (log) log(s); }
    template <class Ar> void stateIO(Ar& ar);   // C14 restart member list (simulator.cpp)
    void updateResvConversion();            // per time step: field averages -> WellRT::resvC, history RESV targets
    void buildAquifers();
    void initAquifers();
    std::vector<AquiferRT> aquifers_;
    std::vector<std::string> wellNames_;
    void buildWells();
    void buildPattern();
    void applyControls(int step);
    void applyWcycle(int step);             // C37 WCYCLE records in force (schedule extra keywords, incl. ACTIONX-added ones)
    void wcycleStep(double& dt, bool& clipped);   // C37 start of every time step: cycle status, start-up ramp, step size
    struct WCycleRec { std::string well; double on = 0, off = 0, start = 0, maxDt = 0; bool ctrl = false; };
    std::vector<WCycleRec> wcyc_;
    void applyWefac(int step);              // C33 WEFAC from the schedule extra keywords (incl. ACTIONX-added ones)
    void updateConnections(int step);       // connection factors for the report step (0 when not yet / no longer open)
    CellFull evalCell(int c, int satOverride = -1) const;   // satOverride: saturation table (C34 directional KRNUM)
    void evalAll(bool restricted = false);   // restricted: only the cells this rank needs (distributed assembly)
    DistAssembly da_; bool daOn_ = false; std::vector<char> daNeed_; void daNeedCells();
    void assemble(double dt);
    bool converged(double dt, double& cnvMax, double& mbMax) const;
    void updateState(const std::vector<double>& dx);
    bool newton(double dt, int& its, int& lin);
    bool checkWellSwitching();
    void computeHeads();
    void storeOld();
    void record(double dt, int its, int lin);
    void snapshot();
    double wellRateOf(const WellRT& w, Control m) const;

    Case case_;
    Units U_{};
    Grid grid_;
    std::vector<FluidModel> fluids_;          // one per PVTNUM region
    std::vector<SatFunctions> sats_;          // one per SATNUM region
    std::vector<int> pvtReg_, satReg_, eqlReg_;   // per active cell, 0-based
    const FluidModel& fl(int a) const { return fluids_[size_t(pvtReg_[size_t(a)])]; }
    const SatFunctions& st(int a) const { return sats_[size_t(satReg_[size_t(a)])]; }
    std::vector<SatFunctions::EndPoints> eps_;  // per active cell end points (ENDSCALE), empty = no scaling
    const SatFunctions::EndPoints* ep(int a) const { return eps_.empty() ? nullptr : &eps_[size_t(a)]; }
    void setupEndPointScaling(const Case& c);

    // ---- relative permeability hysteresis (SATOPTS HYSTER + EHYSTR, Carlson non-wetting phase: gas in gas-oil, oil in
    // oil-water; wetting phases follow drainage). History (maximum non-wetting saturations) and the Carlson shifts are
    // updated after each accepted time step (explicit, as in ECLIPSE / OPM Flow).
    bool hyst_ = false;
    bool hystWetImb_ = false;   // C44 EHYSTR model 1 / 3: wetting phases (krw, krog) on the imbibition curve (as OPM)
    std::vector<int> imbReg_;                       // IMBNUM region per active cell (0-based)
    std::vector<SatFunctions::EndPoints> ieps_;     // imbibition end points per cell (empty = table)
    std::vector<double> sgHy_, soHy_;               // historical maximum Sg (gas-oil) and oil saturation (oil-water)
    std::vector<double> dG_, dO_;                   // Carlson shifts (NaN = drainage curve)
    const SatFunctions& sti(int a) const { return sats_[size_t(imbReg_[size_t(a)])]; }
    const SatFunctions::EndPoints* iep(int a) const { return ieps_.empty() ? nullptr : &ieps_[size_t(a)]; }
    void setupHysteresis(const Case& c);
    // ---- groups (GRUPTREE / GCONPROD / GCONINJE / GCONPRI) and THP limits (VFPPROD / VFPINJ)
    std::map<std::string, std::string> gParent_;
    std::map<std::string, GroupProdControl> gprod_;
    std::map<std::string, GroupInjControl> ginj_[2];      // water, gas
    std::vector<std::vector<std::string>> wellChain_;     // per runtime well: its group, the parent, ..., FIELD
    bool groupsOn_ = false;
    std::set<std::string> vfpMissing_;
    void groupControls(int step);
    void applyGroupTargets();
    GuideRateSpec guide_;                    // C34 GUIDERAT in force (phase NONE: potentials)
    double guideRateOf(WellRT& w, Control m);   // C34 guide rate of a well for the controlled quantity m
    double wellPotential(const WellRT& w, Control m) const;
    void updateThp();
    void updateHysteresis();
    double soOW(int a, double sw, double sg) const;  // oil saturation argument of krow (model dependent)
    // ---- VAPPARS: saturated Rv / Rs reduced where So fell below its historical maximum (explicit factors)
    double vap1_ = 0, vap2_ = 0;
    std::vector<double> soMaxV_, fRv_, fRs_;
    void updateVappars();
    double rvSatC(int a, double p) const { return fl(a).rvSat(p) * (fRv_.empty() ? 1.0 : fRv_[size_t(a)]); }
    double rsSatC(int a, double p) const { return fl(a).rsSat(p) * (fRs_.empty() ? 1.0 : fRs_[size_t(a)]); }

    // ---- polymer flooding (POLYMER): sequential-implicit component transport coupled to the flow solve by outer
    // iterations inside each time step until the water mobility multiplier no longer changes (fully coupled limit)
    bool poly_ = false;
    std::vector<double> pc_;                // polymer concentration per active cell
    std::vector<double> pcAdsMax_;          // historical maximum adsorption (irreversible adsorption)
    std::vector<double> pcMassOld_;         // polymer mass at the old time level
    std::vector<double> wmul_;              // water mobility multiplier (viscosity correction / resistance factor)
    std::vector<double> pcInj_;             // injected concentration per well (WPOLYMER)
    std::vector<double> pcFw_;              // water surface flux per connection (scratch)
    std::vector<int> mixReg_;               // PLMIXNUM region per active cell (0-based)
    bool pcShear_ = false;                  // PLYSHLOG shear effects
    std::vector<double> pcArea_, pcZw_, pcZp_;   // per connection: face area, water / polymer shear factors (lagged)
    std::vector<std::vector<double>> shLogV_, shM_;   // per PVT region: log reference velocities, reference multipliers
    double shearFactor(int pvtReg, double c, double v0) const;
    CSR pcJ_; ILU0 pcIlu_;
    std::vector<int> pcSelf_, pcAdj_;       // positions in pcJ_ for the diagonal / adjacency entries
    double pcCum_[2] = {0, 0};
    int pcOuter_ = 0;                       // outer coupling iterations (statistics)
    struct PolyEval { double wvc, dwvc, pvc, dpvc, ads, dads, rk; };
    PolyEval polyEval(int a, double c) const;
    double polyMass(int a, double c, double* dM) const;
    void polymerSetup();
    // ---- C35 chemical EOR (chem.cpp): surfactant (SURFACT), alkaline (ALKALINE) and foam (FOAM) as sequential-implicit
    // transported components on the converged flow fluxes, coupled to the flow solve by the polymer outer iterations
    struct ChemSp {
        int kind = 0;                           // 0 surfactant, 1 alkaline, 2 foam
        std::string name;                       // SURFACT / ALKALINE / FOAM
        int phase = 0;                          // carrier phase (PW, or PG for foam with FOAMOPTS GAS)
        std::vector<std::vector<double>> ads;   // per SATNUM: (c, adsorbed mass / rock mass)
        std::vector<double> rho;                // per SATNUM: rock density
        std::vector<int> adsIdx;                // per SATNUM: 1 desorption allowed, 2 irreversible
        std::vector<std::vector<double>> dcyW, dcyO;   // foam: half-life vs Sw / So per SATNUM
        std::vector<double> c, adsMax, massOld, adsMul, vol, lam;   // per cell
        std::vector<double> inj;                // per well
        double scale = 0, cumInj = 0, cumProd = 0, cumDecay = 0, decayPending = 0;
    };
    std::vector<ChemSp> chem_;
    int chSurf_ = -1, chAlk_ = -1, chFoam_ = -1;
    std::vector<double> chWmul_, chGmul_, chMisc_, chPcMul_, chLogNc_, chVg_, polyAdsMul_;   // per cell (lagged)
    std::vector<int> chMisSat_;                 // SURFNUM table (0-based) per cell
    std::vector<double> chArea_;                // per connection face area
    std::vector<std::vector<double>> surfVisc_, surfSt_, surfCapd_, alSurfSt_, alSurfAd_, alPolAds_, foamMob_, foamMobP_, foamMobS_;
    struct ChemBak { std::vector<std::vector<double>> c; std::vector<double> w, g, m, pc, pa; } chBak_;
    int chOuterMax_ = 0;
    CSR chJ_; ILU0 chIlu_;
    std::vector<int> chSelf_, chAdj_;
    double chMaxRelPc_ = 0;
    // ---- C35 SOLVENT (solvent.cpp): the fourth (solvent) component. Flow is fully implicit in the total gas-like
    // component G+S (gas incl. dissolved + solvent) with the solvent mass fraction z of G+S lagged; the solvent
    // equation is then solved for z on the converged fluxes and iterated with the flow solve (outer coupling)
    bool solv_ = false, solvMisc_ = false;
    std::vector<double> zs_, zsBak_, ntOld_;     // per cell: z = Ns / (Ns + Ng), backup, old-time solvent amount z N_T
    std::vector<double> svInjF_;                // per well: WSOLVENT solvent fraction of the injected gas
    std::vector<Table1D> svB_, svBMu_;          // per PVTNUM: PVDS 1/Bs, 1/(Bs mu_s)
    std::vector<double> svRho_;                 // per PVTNUM: SDENSITY surface density
    std::vector<Table1D> ssfnG_, ssfnS_;        // per SATNUM: SSFN gas / solvent kr multipliers vs phase fraction
    std::vector<Table1D> misc_, pmisc_, msfnG_, msfnO_, sorwmis_, sgcwmis_, tlpmix_, sof2_;   // miscible tables (per MISCNUM / SATNUM)
    std::vector<double> tlmix_;                 // per MISCNUM: TLMIXPAR viscosity mixing parameter
    std::vector<int> miscReg_;                  // per cell MISCNUM (0-based)
    std::vector<double> svF_;                   // per connection: converged solvent surface flux (foam carrier)
    std::vector<std::vector<std::array<double, 2>>> svWq_;   // per well: (cell, solvent surface rate + produced / - injected)
    double svCum_[2] = {0, 0};                  // FNIT, FNPT
    int svOuter_ = 0;
    CSR svJ_; ILU0 svIlu_;
    std::vector<int> svSelf_, svAdj_;           // positions in svJ_ (diagonal / adjacency)
    void solventSetup();
    void solventControls(int step);
    // ---- C35 gas-lift optimisation (LIFTOPT / WLIFTOPT / GLIFTOPT; gaslift.cpp)
    struct GasLiftWell { bool on = false; double fixed = -1, maxRate = -1, weight = 1, minRate = 0, gasW = 0; };
    struct GasLiftGroup { double maxLift = -1, maxTotal = -1; };
    double gloInc_ = 0, gloEco_ = 0, gloInterval_ = 0, gloLast_ = -1e30, gloCum_ = 0;
    // C35 facilities: satellite groups (GSATPROD: oil, water, gas, RESV, lift gas), gas consumption / import (GCONSUMP)
    struct SatProd { double o = 0, w = 0, g = 0, resv = 0, lift = 0; };
    std::map<std::string, SatProd> gsat_;
    struct GasConsump { double cons = 0, imp = 0; std::string node; };
    std::map<std::string, GasConsump> gcons_;
    struct GasSale { double target = -1, maxRate = -1, minRate = -1; std::string proc; };
    std::map<std::string, GasSale> gsale_;      // C35 GCONSALE: gas sales target / limits of a group (met by its gas injection)
    double groupGasAvail(const std::string& g) const;   // produced (incl. satellite) - consumed + imported gas below a group
    double gasBalCum_[3] = {0, 0, 0};          // cumulative consumption, import, sales
    bool netStd_ = false;                       // standard network (GRUPNET) instead of the extended one
    std::vector<GasLiftWell> glw_;              // per runtime well
    std::map<std::string, GasLiftGroup> glg_;
    int gloRuns_ = 0;
    std::set<std::string> gloNotes_;
    void gasLiftControls(int step);             // LIFTOPT / WLIFTOPT / GLIFTOPT records of the report step; ALQ of the other wells
    bool gasLiftOptimize();                     // start of a time step: lift gas allocation (true = an ALQ changed)
    const VfpTable* wellVfp(const WellRT& w) const;
    bool thpOperatingPoint(const WellRT& w, double alq, bool liftRate, double& bhp, double* q) const;
    double liftGasRate(const WellRT& w) const;  // lift gas injection rate of a producer (0 without gas lift)
    double solventTransport(double dt);         // returns the largest z change (< 0: failure)
    void solventAccept();
    void solventFixStates();                    // gas phase appears where the transport brought solvent into an undersaturated cell
    double solventTotalNT(int a) const;         // N_T = Ng + Ns of a cell at the current state (surface gas units)
public:
    bool hasSolvent() const { return solv_; }
    const std::vector<double>& solventFraction() const { return zs_; }
    double solventInPlace() const;
    double solventCum(int k) const { return svCum_[k & 1]; }   // 0 injected, 1 produced
    // C35 gas lift: ALQ in force per runtime well (lift gas rate for optimised wells), optimisations run so far
    double wellAlq(size_t i) const { return i < wells_.size() ? wells_[i].alq : 0.0; }
    int gasLiftRuns() const { return gloRuns_; }
private:
    void chemSetup();
    void chemControls(int step);
    double chemTransport(double dt);            // returns the largest change of the coupling multipliers (< 0: failure)
    bool chemSolve(ChemSp& S, double dt, const std::vector<double>& F, const std::vector<std::vector<std::array<double, 2>>>& wq);
    double chemEffects();                       // concentrations -> mobility / IFT / miscibility multipliers (returns the change)
    void chemAccept();
    void chemBackup();
    void chemRestore();
    double chemMass(const ChemSp& S, int a, double c, double* dM) const;
public:
    // C35 chemical components: names / per-cell concentrations / in-place mass (free + adsorbed)
    size_t numChem() const { return chem_.size(); }
    const std::string& chemName(size_t i) const { return chem_[i].name; }
    const std::vector<double>& chemConc(size_t i) const { return chem_[i].c; }
    double chemInPlace(size_t i, double* adsorbed = nullptr) const;
    double chemDecayed(size_t i) const { return i < chem_.size() ? chem_[i].cumDecay : 0.0; }   // cumulative decayed amount (foam)
private:
    // C27 passive tracers (tracer.cpp): one linear implicit upwind transport solve per tracer after each converged step
    std::vector<TracerInfo> trInfo_;
    std::vector<double> trC_, trVolOld_;      // concentrations [tr * na + a]; free-phase surface volume [tr * na + a] at the old time
    std::vector<std::vector<double>> trInj_;  // [tracer][well] injected concentration (WTRACER)
    std::vector<double> trDecayCum_;          // cumulative decayed amount per tracer
    std::vector<double> trB_, trS_;           // C28 b and S per phase and cell at the new state [3 * a + p] (partitioning / dispersion)
    std::vector<int> trOrder_;                // C28 solve order (decay-chain parents first)
    CSR trJ_; ILU0 trIlu_; std::vector<int> trSelf_, trAdj_;
    void tracerSetup();                       // after the initial equilibrium (initialize)
    void tracerControls(int step);            // WTRACER records of the report step
    void tracerVolumes();                     // trVolOld_ from the current state
    // C29 salt / dry-out, surface network, wellbore friction
    std::vector<double> saltS_, saltK_, saltEvap_, saltPvF_, saltTF_;   // solid fraction, k factor, evaporation sink rate, applied PV / T factors
    double saltEvapCum_ = 0, saltOmega_ = 0, saltSolidPrev_ = 0;
    int saltTr_ = -1;                         // index of the internal NACL tracer
    void saltSetup();
    void saltStep(double dt);
    void saltEquilibrium();
    void saltApply();
    std::vector<double> saltPv0_, saltT0_;   // start-of-run pore-volume references / transmissibilities
                 // precipitation / dissolution, evaporation rates, porosity / permeability update
    // C31 UDQ / UDA / ACTIONX (sim_udq.cpp)
    std::shared_ptr<udq::Set> udq_;
    std::shared_ptr<std::vector<udq::Action>> actions_;   // (shared_ptr: udq.h stays out of this header)
    std::vector<double> wCum_;                // per well: oil / water / gas produced, water / gas injected (surface)
    std::vector<ActionEvent> actEvents_;
    std::vector<msw::Well> msw_;                // C32
    std::vector<int> mswOf_;                    // per well: index into msw_ or -1
    std::vector<std::vector<int>> mswSegOf_;    // per well, per connection: segment position (-1 none)
    std::vector<msw::Solution> mswLast_;        // last solve per msw_ entry
    std::vector<MswRecord> mswRec_;
    std::string mswSig_;                        // definitions signature (re-setup only when they change)
    void mswSetup(int step);
    void mswMixture(WellRT& w, const std::vector<double>& gIn);   // C45 segment mixtures for producer cross-flow
    void mswHeads(WellRT& w, const std::vector<double>& qres, const std::vector<double>& rho, const std::vector<double>& mu,
                  const std::vector<double>& pCell, const std::vector<double>& aRes,
                  const std::vector<double>* gasFrac = nullptr, const std::vector<double>* rhoGas = nullptr, const std::vector<double>* rhoLiq = nullptr,
                  const std::vector<double>* qSurf = nullptr);   // C37 connection surface component rates (3 per connection)
    double headDt_ = 0;                         // C37 time step of the current attempt (segment storage)
    void mswRecord();
    std::vector<int> actStep_; std::vector<std::string> actText_;   // applied action bodies (restart replay)
    std::vector<ReportStep> baseSchedule_; std::vector<WellInput> baseWells_;   // schedule before any action
    std::vector<std::string> udqGroups_;
    void udqSetup();
    void udqEvaluate(const SummaryPoint* s, SummaryPoint* out, bool assignsOnly = false);
    void udqSource(const SummaryPoint& s, udq::Source& src) const;
    void applyUda();
    void actionsEndOfStep();
    void pyActionSetup();                       // C37 PYACTION (sim_udq.cpp)
    void pyActionsEndOfStep();
    std::shared_ptr<pyact::Host> pyHost_;
    std::shared_ptr<std::vector<pyact::Action>> pyActions_;
    bool pyFailed_ = false;
    void actionApply(int step, const std::string& body, bool quiet);
    void refreshConnectionWindows(bool quiet);
    bool rebuildActionConnections(bool quiet);   // C32: connections added by ACTIONX COMPDAT during the run
    std::vector<NetNode> netNodes_;
    std::vector<NetRecord> netHist_;
    int netStep_ = -1;
    void networkSetup(int step);              // NODEPROP / BRANPROP of the report step (latest definitions)
    void facilityControls(int step);          // C35 GSATPROD / GCONSUMP records of the report step
    bool groupBelow(const std::string& g, const std::string& top) const;   // g == top or a subgroup of top (GRUPTREE)
    double satRate(const std::string& top, int what) const;   // satellite rate below a group: 0 oil 1 water 2 gas 3 RESV 4 lift
    void networkBalance();                    // node pressures from the latest group rates -> well THP limits
    bool networkPressures(std::vector<double>* p, std::vector<int>* wnode);   // C30: node pressures from the current rates
    double chokePressure(size_t node, double pUp, const std::vector<int>& wnode) const;   // C35 automatic choke: THP meeting the group target
    bool chokeGroup(const std::string& g) const;   // C35: the group's rate target is met by an automatic network choke
    void networkApply(const std::vector<double>& p, const std::vector<int>& wnode);   // C30: node pressures -> THP limits
    // C36 NETBALAN: interval (days; 0 every time step, < 0 NUPCOL = every time step here), node pressure tolerance (case
    // pressure), max iterations, max balance error (case pressure; exceeded -> the time step is cut); given = keyword present
    struct NetBalan { bool given = false; double interval = 0, ptol = 0, maxErr = -1, targetErr = -1; int maxIter = 10; } netBal_;
    double netLastBal_ = -1e300;
    void networkCompressors();
    double compressorGas() const;             // current gas consumption of the compressors that are on (surface gas rate)
    int netIterMax_ = 12; long netIters_ = 0, netResolves_ = 0; std::set<int> netWarn_;   // C30 implicit network balance statistics
    void wellboreFriction(WellRT& w, const std::vector<double>& qres, const std::vector<double>& rhoSeg, const std::vector<double>& pCell, const std::vector<double>& aRes);
    void tracerTransport(double dt);          // converged step -> new concentrations, well tracer rates (WellRates::comp)
    void polymerControls(int step);
    // DRSDT: limit on the growth of Rs per time step (OPM: lastRs + DRSDT dt)
    double drsdt_ = -1; bool drsdtFree_ = false;
    std::vector<double> lastRs_, rsLim_;
    void drsdtControls(int step);
    void drsdtAccept();
    void polymerMultipliers();
    double polymerTransport(double dt);     // returns the largest relative change of the water mobility multiplier (< 0: failed)
    void polymerAccept();
    bool W_ = false, O_ = false, G_ = false, live_ = false, vap_ = false;
    int nv_ = 1, iSw_ = -1, iX_ = -1;
    int eq_[3] = {-1, -1, -1};             // phase -> equation offset

    // state
    std::vector<double> p_, sw_, x_;        // x = Sg (hasGas) or Rs
    std::vector<uint8_t> hasGas_;           // cell state: 0 = no free gas (x = Rs), 1 = oil + gas (x = Sg), 2 = no oil (x = Rv)
    std::vector<double> accOld_;            // na * 3
    // C15 unconventional: Langmuir adsorbed gas coefficient per active cell (surface gas volume at p >> P_L) and P_L;
    // stress-dependent permeability exponent per active cell and the initial pressure it is measured from
    std::vector<double> langK_, sGam_, pI_;
    double langPL_ = 0;
public:
    double adsorbedGas() const;                       // total adsorbed gas (surface volume)
    int numHydraulicFractureCells() const { return grid_.numHfCells(); }
private:
    std::vector<CellFlux> cf_;
    std::vector<int> krReg_[3];               // C34 KRNUMX/Y/Z (0-based tables per active cell; empty = off)
    std::vector<CellFlux> cfDir_[3];          // C34 directional mobilities (flux terms of X / Y / Z connections)
    // C45 memory: the solvent mobility / density live outside CellFlux (64 B per cell less without SOLVENT)
    std::vector<std::array<AD3, 2>> cfS_;     // per cell {mobS, rhoS} (SOLVENT only)
    std::vector<AD3> cfDirS_[3];              // directional solvent mobility (KRNUMX/Y/Z with SOLVENT)
    void storeCell(int c, const CellFull& f) {
        cf_[size_t(c)] = static_cast<const CellFlux&>(f);
        if (!cfS_.empty()) { cfS_[size_t(c)][0] = f.mobS; cfS_[size_t(c)][1] = f.rhoS; }
    }
    const AD3& mobSDir(int kd, int c) const { return kd >= 0 && !cfDirS_[kd].empty() ? cfDirS_[kd][size_t(c)] : cfS_[size_t(c)][0]; }
    bool dirKr_ = false;
    void evalDirKr(bool restricted);
    void applyWelpi(int step);                // C34 WELPI
    void applySchedMultflt(int step);         // C46 MULTFLT in SCHEDULE (cumulative, from the event's step on)
    std::vector<char> multfltDone_;
    std::array<double, 3> wellPIs(size_t wi) const;   // C46 surface-volume PI per phase {oil, water, gas} (OPM convention)
    std::vector<int> completionLumps(size_t wi) const;                  // C34 COMPLUMP completion of each connection
    std::vector<std::array<double, 3>> connRatesOf(size_t wi) const;    // C34 per-connection surface rates (o, w, g)
    // C34 rock compaction (ROCKCOMP): region, overburden and history per active cell
    bool rockOn_ = false, rockTm_ = false;
    std::vector<int> rockReg_;
    std::vector<double> rockPob_, rockPmin_, rockSwMax_, rockSw0_, rockP0_;
    void setupRockComp();
    void updateRockHistory();
    AD3 rockMult(int c, const AD3& P, const AD3& Sw, bool trans) const;
    std::vector<WellRT> wells_;
    std::vector<std::vector<int>> cellWells_;

    // linear system
    // C45 block assembly (blkAsm_): with the block solver the cell-cell Jacobian blocks are assembled straight into the
    // Schur matrix (schur_), J_ then keeps only the well columns of the cell rows and the well rows - no scalar copy of
    // the cell Jacobian (-0.8 kB per cell). Off for the scalar solver, adjoint recording and the matrix dump tools.
    bool blkAsm_ = false, forceScalarJ_ = false;
    double* jrow(int a, int e, int slot) {
        return blkAsm_ ? schur_.cellRow(size_t(blkStart_[size_t(a)] + slot), e) : &J_.val[size_t(J_.rowPtr[size_t(a) * nv_ + e] + slot * nv_)];
    }
    int jwOff(int a) const { return blkAsm_ ? 0 : (blkStart_[size_t(a) + 1] - blkStart_[size_t(a)]) * nv_; }
    CSR fullJacobian() const;
    bool leanLinear() const {           // C45 low-memory BILU0 (TUNING leanMemory; auto from 500 000 active cells)
        static const char* env = std::getenv("RESSIM_LEAN");
        if (env) return std::atoi(env) != 0;
        return case_.tuning.leanMemory > 0 || (case_.tuning.leanMemory < 0 && grid_.numActive() >= 500000);
    }          // scalar CSR of the whole system in block-assembly mode (before the Schur step)
    void ensureScalarJacobian() { if (blkAsm_) { forceScalarJ_ = true; buildPattern(); } }
    CSR J_;
    ILU0 ilu_;
    CPR cpr_;
    WellSchur schur_;                // block path: wells Schur-eliminated, BSR cell system
    BILU0 bilu_;
    BILU0f biluf_;                    // C16 mixed precision
    BlockCPR bcpr_;
    int amgAge_ = 0; double amgBaseIts_ = 0; bool amgStale_ = true;   // C26 AMG reuse across time steps (quality monitor)
    long amgCrossReuse_ = 0;
    std::vector<double> bc_;
    bool useCPR_ = false;
    std::vector<double> accD_;
    std::vector<double> pcScale_;
    std::vector<double> thp_;               // threshold pressure per connection (THPRES), empty = none
    // GRAVDR (dual porosity gravity drainage): matrix block height per connection (0 = none) and, per cell, the normalised
    // water and gas column heights (segregated fluids in the block) - see setupGravityDrainage()
    std::vector<double> gdZ_;
    std::vector<AD3> gdH_;
    bool gravdr_ = false;
    void setupGravityDrainage(const Case& c);
    void evalGravityHeights(int c, const CellFull& f);
    void setupThresholds();           // per-cell water-oil capillary pressure multiplier (SWATINIT), empty = 1              // true-IMPES accumulation blocks (na * nv * nv)
    std::vector<double> R_;
    std::vector<int> blkStart_, blkCells_, adjPos_, selfPos_;
    std::vector<double> bRef_;              // average b per phase (for convergence norms)

    double time_ = 0, dtNext_ = 1;
    int step_ = 0;
    double cum_[5] = {0, 0, 0, 0, 0};
    std::vector<double> initRs_, initRv_;
    std::vector<SummaryPoint> summary_;
    std::vector<Snapshot> snaps_;
    RunStats stats_;
};

} // namespace ressim
