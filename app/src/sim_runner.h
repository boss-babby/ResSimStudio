// Runs the simulator on a worker thread; results are shared with the UI under a mutex.
#pragma once
#include "ressim/simulator.h"
#include "ressim/components.h"
#include "ressim/results_file.h"
#include "ressim/faultmech.h"
#include "ressim/geocouple.h"
#include "ressim/hf_ddm.h"
#include "opm_engine.h"
#include <thread>
#include <mutex>
#include <atomic>
#include <memory>
#include <string>
#include <vector>

struct RunResults {
    std::vector<ressim::SummaryPoint> summary;
    std::vector<ressim::Snapshot> snaps;
    std::vector<std::string> wellNames;
    std::vector<ressim::WellType> wellTypes;
    std::shared_ptr<const ressim::Grid> grid;   // grid the results belong to
    ressim::Case input;                         // copy of the case that was run
    ressim::InPlace initial;
    ressim::RunStats stats;
    double totalTime = 0;
    bool valid = false;
    std::shared_ptr<ressim::Simulator> sim;     // kept after the run finishes (streamlines / restore of report steps)
    std::string engine = "ResSim";              // "ResSim", "ResSim compositional", "OPM Flow" or "ECLIPSE"
    std::string opmDeck;                        // OPM run: exported deck (results next to it)
    // multi-component runs (compositional / CO2STORE / thermal): component list of WellRates::comp, fcompP / fcompI
    std::vector<ressim::ComponentInfo> components;
    std::string trackName;                      // storage / trapping inventory component ("CO2"; thermal: Case::track)
    // ---- C13: geothermal / EDFM analysis of thermal runs (per report step, parallel to snaps)
    struct FracSeg { int frac = 0; int host = -1; float xa = 0, ya = 0, xb = 0, yb = 0, zt = 0, zb = 0; };   // model coordinates (case units)
    std::vector<std::string> fracNames;         // fracture names (index FracSeg::frac)
    std::vector<FracSeg> fracSegs;              // one per EDFM fracture cell
    struct FracState { std::vector<float> ap, p, T, cond; };   // per segment: aperture (mm), pressure / temperature (case units), conductivity (mD m)
    std::vector<FracState> fracSnap;            // per report step
    std::vector<std::vector<float>> wellT;      // [report][well] well-stream temperature (case units; 0 = no flow)
    std::vector<std::vector<float>> wellMW;     // [report][well] enthalpy rate (MW; production > 0, injection < 0)
    std::vector<float> cooledFrac;              // [report] fraction of the matrix bulk volume cooled by more than 10 K
    double heatInPlace = 0;                     // J above the injection temperature (rock + water, matrix cells, t = 0)
    double tInjRef = 0;                         // reference (injection) temperature, case units
    std::vector<std::pair<std::string, std::string>> wellFracLinks;   // well - fracture connections (names)
    // ---- C14: engine restart states saved during the run (exact continuation from these report steps)
    std::vector<ressim::ResultsFile::Restart> restarts;
    int restartedFrom = -1;                     // > 0: this run continued from that report step of a base run
    std::string sourceFile;                     // .rsres the results were loaded from / saved to
    // ---- C24: coupled fault mechanics / induced seismicity and fracture contact (G1-c) of the main run
    std::shared_ptr<const ressim::faultmech::Output> faultMech;   // per fault face and report step (null: feature off)
    std::shared_ptr<const ressim::geocouple::Output> geoCouple;   // per hydraulic-fracture piece and report step
    std::shared_ptr<const ressim::hfddm::Output> hfDdm;            // C26 HFGROW ... DDMWIDTH: DDM apertures per piece and report step
    std::shared_ptr<const std::vector<ressim::Simulator::PfRecord>> hfPf;   // C28 HFGROW ... PHASEFIELD: in-run phase-field crack per report step
    std::vector<ressim::Simulator::ActionEvent> actions;   // C31 ACTIONX blocks that fired (time, report step, wells)
    std::vector<ressim::msw::Well> msw;                     // C32 multisegment wells in force (last definitions)
    std::shared_ptr<const std::vector<ressim::Simulator::MswRecord>> mswRec;   // C32 segment state per accepted time step
    std::vector<ressim::Simulator::NetNode> netNodes;                            // C35-4b surface network at the end of the run
    std::shared_ptr<const std::vector<ressim::Simulator::NetRecord>> netRec;    // node pressures per accepted step (node order = netNodes prefix)
};

// C14: continue a run from a saved report step (restart). The base results up to that step are copied into the new
// results; the case may have a different schedule after the restart step (same grid / wells / components).
struct RestartRequest {
    int step = 0;                               // report step (number of completed report steps of the base run)
    std::string state;                          // engine blob (ResultsFile::Restart::state)
    std::string kind;                           // "BO", "COMP", "THERMAL"
    std::shared_ptr<const RunResults> base;     // base results (summary / snapshots / analysis up to the step)
};

// C18: a ResSim run in an external process (mpiexec + ResSimCLI --mpi-dd) whose results come back as a .rsres file
struct CliRunSpec { std::string exe; std::vector<std::string> args; std::string workDir, resultsFile, engineName; };

class SimRunner {
public:
    void startCli(const ressim::Case& c, const CliRunSpec& s);
    ~SimRunner() { cancel(); join(); }
    void start(const ressim::Case& c);
    void startOpm(const ressim::Case& c, const OpmOptions& o);    // external OPM Flow run
    void startComp(const ressim::Case& c, const std::string& deckPath);   // E300 deck on the compositional core
    void startThermal(const ressim::Case& c);                              // E300 THERMAL case on the thermal core
    // C14 restart control: states saved every N report steps (0 = none) and / or at the end of the run; a pending
    // restart request is consumed by the next start / startComp / startThermal
    void setRestartSaving(int everyN, bool last) { rstEvery_ = everyN; rstLast_ = last; }
    void requestRestart(std::shared_ptr<RestartRequest> rq) { rq_ = std::move(rq); }
    void cancel() { cancel_ = true; std::lock_guard<std::mutex> l(procM_); if (proc_) proc_->kill(); }
    bool running() const { return running_; }
    void join() { if (th_.joinable()) th_.join(); }
    void setResults(RunResults r) { if (running_) return; join(); { std::lock_guard<std::mutex> l(m_); r_ = std::move(r); error_.clear(); } ++version_; }
    void clear() { if (running_) return; join(); { std::lock_guard<std::mutex> l(m_); r_ = RunResults{}; error_.clear(); } ++version_; }

    // Access to shared state (lock while reading)
    std::mutex& mutex() { return m_; }
    RunResults& results() { return r_; }
    std::vector<std::string>& log() { return log_; }
    std::string error() { std::lock_guard<std::mutex> l(m_); return error_; }
    double progress() const { return progress_; }
    double wallSeconds() const { return wall_; }
    unsigned version() const { return version_; }     // incremented when new results arrive
    void touch() { ++version_; }                      // C30: results changed in place (post-processing)
    void addLog(const std::string& s) { std::lock_guard<std::mutex> l(m_); log_.push_back(s); }

private:
    void work(ressim::Case c);
    void workOpm(ressim::Case c, OpmOptions o);
    void workComp(ressim::Case c, std::string deckPath);
    void workThermal(ressim::Case c);
    void workCli(ressim::Case c, CliRunSpec s);
    int rstEvery_ = 0; bool rstLast_ = true;
    std::shared_ptr<RestartRequest> rq_, rqActive_;
    bool wantsRestart(int step, bool last) const { return (rstLast_ && last) || (rstEvery_ > 0 && step % rstEvery_ == 0); }
    void prefillFromBase(const RestartRequest& rq);   // copy the base results up to the restart step (under m_)
    std::mutex procM_;
    ChildProcess* proc_ = nullptr;
    std::thread th_;
    std::mutex m_;
    RunResults r_;
    std::vector<std::string> log_;
    std::string error_;
    std::atomic<bool> cancel_{false}, running_{false};
    std::atomic<double> progress_{0}, wall_{0};
    std::atomic<unsigned> version_{0};
};
