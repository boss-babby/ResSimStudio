// Runs the simulator on a worker thread; results are shared with the UI under a mutex.
#pragma once
#include "ressim/simulator.h"
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
    std::string engine = "ResSim";              // "ResSim" or "OPM Flow"
    std::string opmDeck;                        // OPM run: exported deck (results next to it)
};

class SimRunner {
public:
    ~SimRunner() { cancel(); join(); }
    void start(const ressim::Case& c);
    void startOpm(const ressim::Case& c, const OpmOptions& o);    // external OPM Flow run
    void cancel() { cancel_ = true; std::lock_guard<std::mutex> l(procM_); if (proc_) proc_->kill(); }
    bool running() const { return running_; }
    void join() { if (th_.joinable()) th_.join(); }
    void clear() { if (running_) return; join(); { std::lock_guard<std::mutex> l(m_); r_ = RunResults{}; error_.clear(); } ++version_; }

    // Access to shared state (lock while reading)
    std::mutex& mutex() { return m_; }
    RunResults& results() { return r_; }
    std::vector<std::string>& log() { return log_; }
    std::string error() { std::lock_guard<std::mutex> l(m_); return error_; }
    double progress() const { return progress_; }
    double wallSeconds() const { return wall_; }
    unsigned version() const { return version_; }     // incremented when new results arrive
    void addLog(const std::string& s) { std::lock_guard<std::mutex> l(m_); log_.push_back(s); }

private:
    void work(ressim::Case c);
    void workOpm(ressim::Case c, OpmOptions o);
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
