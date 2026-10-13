// Fully-implicit black-oil simulator (water / oil / gas with dissolved gas).
//   Primary variables per cell: Po, Sw, and Sg (free gas present) or Rs (undersaturated)
//   -> variable switching. Wells are fully coupled with BHP as an extra unknown.
#pragma once
#include "ressim/case.h"
#include "ressim/fluid.h"
#include "ressim/grid.h"
#include "ressim/linsolve.h"
#include <vector>
#include <string>
#include <functional>
#include <atomic>
#include <cstdint>

namespace ressim {

enum Phase { PW = 0, PO = 1, PG = 2 };

struct WellRates {
    double oil = 0, water = 0, gas = 0;        // surface production rates (>0)
    double waterInj = 0, gasInj = 0;           // surface injection rates (>0)
    double resv = 0;                           // reservoir voidage rate
    double polyInj = 0, polyProd = 0;          // polymer injection / production rate (mass / day, POLYMER option)
    double bhp = 0;
    double oilH = -1, waterH = -1, gasH = -1;  // observed (history) rates, -1 if none
    Control mode = Control::BHP;
    bool open = false;
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
    std::vector<WellRates> wells;
};

struct Snapshot {                              // per active cell, at a report time
    double time = 0;
    std::vector<float> pressure, sw, sg, so, rs;
    std::vector<float> rv;                      // vaporised oil (only with VAPOIL)
    std::vector<float> poly, polyAds;           // polymer concentration / adsorbed polymer (only with POLYMER)
};

struct InPlace { double oil = 0, water = 0, gasFree = 0, gasDissolved = 0, oilVaporized = 0; };   // oil = liquid oil

struct AquiferState { int id = 0; bool fetkovich = true; int cells = 0; double pressure = 0, rate = 0, cumulative = 0; };

struct RunStats { int timesteps = 0, cuts = 0, newton = 0, linear = 0; double seconds = 0; bool cpr = false; int amgLevels = 0; int polyOuter = 0; };

class Simulator {
public:
    // Throws std::runtime_error with a readable message on invalid input.
    void setup(const Case& c);
    void initialize();                          // hydrostatic equilibrium (EQUIL)
    bool runReportStep();                       // false when the schedule is finished or cancelled
    bool finished() const { return step_ >= int(case_.schedule.size()); }
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
    // restore the state of a saved report step (pressure, saturations, well BHPs) for post-processing
    void restoreReport(int snapIndex);
    // set the state from an external snapshot (e.g. an ECLIPSE / OPM restart step mapped to this grid's active cells)
    // with well BHPs in wellNames() order (<= 0: shut) - wells are put on BHP control; used for streamlines of
    // results computed by another simulator. Throws if the snapshot size does not match.
    void restoreExternal(const Snapshot& sn, int reportStep, const std::vector<double>& bhp);

    std::function<void(const std::string&)> log;              // optional message sink
    std::function<void(const SummaryPoint&)> onTimestep;      // called after each accepted step
    std::atomic<bool>* cancel = nullptr;

private:
    void restoreCells(const Snapshot& sn);
    struct CellFlux {                       // compact per-cell properties used by flux terms
        AD3 P[3], mob[3], rho[3], mobRs, mobRv, acc[3];
    };
    struct CellFull : CellFlux { AD3 b[3], S[3], rs, rv; };
    struct WellRT {                         // runtime well data
        int input = -1;
        WellType type = WellType::Producer;
        std::vector<int> cell;              // active cell per open completion
        std::vector<double> wi, cdepth;
        std::vector<double> wiBase;         // connection factor when the connection is open
        std::vector<int> cFrom, cUntil;     // report-step window of each connection (late COMPDAT / COMPDAT SHUT)
        int nActive = 0;                    // connections open in the current report step
        std::vector<int> cpos;              // position of completion cell in well row cell list
    std::vector<double> qwConn;         // surface water rate per connection at the last assembly (+ out of the reservoir)
    std::vector<double> rw, re, plen;   // wellbore radius, area-equivalent block radius, perforation length per connection
    std::vector<double> wmulConn, zwConn;   // polymer: connection water mobility multiplier and shear factor in it
        std::vector<int> cellWellIdx;       // index of this well in the cell's well list
        std::vector<int> rowCells;          // unique sorted cells for the well row
        double refDepth = 0, headRho = 0;
        std::vector<double> dpConn;         // wellbore pressure difference bhp -> connection (OPM/ECLIPSE segregated densities)
        double pConn(size_t k, double g) const { return bhp + (k < dpConn.size() ? dpConn[k] : headRho * g * (cdepth[k] - refDepth)); }
        double bhp = 0;
        WellControl ctrl;                   // control record in force
        Control mode = Control::BHP;        // active mode (after switching)
        bool open = false;
        WellRates rates;
    };

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
    };
    void log_(const std::string& s) const { if (log) log(s); }
    void buildAquifers();
    void initAquifers();
    std::vector<AquiferRT> aquifers_;
    std::vector<std::string> wellNames_;
    void buildWells();
    void buildPattern();
    void applyControls(int step);
    void updateConnections(int step);       // connection factors for the report step (0 when not yet / no longer open)
    CellFull evalCell(int c) const;
    void evalAll();
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
    std::vector<CellFlux> cf_;
    std::vector<WellRT> wells_;
    std::vector<std::vector<int>> cellWells_;

    // linear system
    CSR J_;
    ILU0 ilu_;
    CPR cpr_;
    bool useCPR_ = false;
    std::vector<double> accD_;
    std::vector<double> pcScale_;
    std::vector<double> thp_;               // threshold pressure per connection (THPRES), empty = none
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
