// C33 scripting / automation layer shared by the C API (capi.h -> ressim.dll / libressim.so -> Python package
// `ressim`), the CLI (--regress, --report) and the regression / report generators.
//
//   Session s; s.open("SPE1.DATA");            // ECLIPSE deck (.DATA / any text), .rsim project or "example:ID"
//   s.setKeyword("PERMX=300*200");              // edit the deck before the run (same rules as ResSimCLI --set)
//   s.run();                                    // any engine (black oil incl. coupled mechanics, compositional, thermal)
//   or: s.begin(); while (s.step() > 0) { ... s.applySchedule("WCONPROD\n 'P1' OPEN ORAT 500 /\n/\n"); }
//   s.summary() -> ecl::Summary (FOPR, WBHP:P1, ... as in ECLIPSE), s.output().snaps -> per-cell arrays
#pragma once
#include "ressim/case.h"
#include "ressim/ecl.h"
#include "ressim/simulator.h"
#include <atomic>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ressim {
namespace api {

struct RunOutput {
    bool ok = false;
    std::string error, engine;                 // engine: "ResSim black oil", "ResSim compositional", ...
    std::vector<SummaryPoint> summary;
    std::vector<Snapshot> snaps;               // black oil: initial state + one per report step (per active cell)
    std::vector<int> globalOf;                 // active -> global cell (black oil, after begin / run)
    std::vector<double> depth, porv;           // per active cell (black oil)
    RunStats stats;
    double seconds = 0;
    std::vector<std::string> wellNames;
    std::vector<std::pair<std::string, double>> kpi;   // coupled mechanics scalars
    InPlace initial;
    std::vector<std::string> log;              // engine messages (ACTIONX, warnings, ...)
};

// "KEY=VALUE" deck edit (ECLIPSE keyword: data up to its '/' replaced; RESSIM:KEY: the '-- RESSIM: KEY' line)
bool applySet(std::string& deckText, const std::string& kv, std::string& err);
// engine that runs this case
std::string engineName(const Case& c);
// runs any engine headless (one thread per run unless c.tuning.threads); progress(t_days) per report step
RunOutput runCase(const Case& c, const std::function<void(double)>& progress = {}, std::atomic<bool>* cancel = nullptr);
// run -> ECLIPSE-style summary vectors (FOPR, FWCT, FPR, WOPR:NAME, WBHP:NAME, ...)
ecl::Summary toSummary(const RunOutput& r, const Case& c);

class Session {
public:
    Session();
    ~Session();
    bool open(const std::string& path);                       // .DATA (or any deck text file), .rsim, "example:ID"
    bool openText(const std::string& deckText, const std::string& baseDir = "");
    bool setKeyword(const std::string& kv);                   // re-imports the edited deck
    bool run(const std::function<void(double)>& progress = {}, std::atomic<bool>* cancel = nullptr);
    // step-wise black-oil run with steering from a script
    bool begin();                                             // setup + equilibration
    int step();                                               // one report step: 1 more to come, 0 finished, -1 error
    bool applySchedule(const std::string& scheduleText, std::vector<std::string>* messages = nullptr);
    double time() const;
    int reportStep() const;
    int numReportSteps() const { return int(c_.schedule.size()); }

    const Case& caseData() const { return c_; }
    Case& mutableCase() { return c_; }                        // edits before run() / begin()
    const RunOutput& output() const { return out_; }
    const ecl::Summary& summary();                            // built lazily from the run
    const std::vector<std::string>& warnings() const { return warnings_; }
    const std::string& error() const { return error_; }
    const std::string& deckText() const { return text_; }
    bool hasResults() const { return !out_.summary.empty(); }
    // static per-cell input arrays (global cells, ECLIPSE order): PORO PERMX PERMY PERMZ NTG DX DY DZ TOPS ACTNUM SATNUM
    // PVTNUM EQLNUM FIPNUM; empty when unknown
    std::vector<double> staticArray(const std::string& name) const;
    // dynamic per-active-cell array of snapshot s (PRESSURE SWAT SOIL SGAS RS RV ... and Snapshot::extra names)
    std::vector<float> cellArray(int snap, const std::string& name) const;
    std::vector<std::string> cellArrayNames(int snap) const;
    bool writeEclipse(const std::string& base, std::string* err = nullptr) const;   // EGRID/INIT/SMSPEC/UNSMRY/UNRST
    bool exportDeck(const std::string& path, std::string* err = nullptr) const;

private:
    void finishStepping();
    Case c_;
    std::string text_, baseDir_, error_;
    bool fromProject_ = false;
    std::vector<std::string> warnings_;
    std::unique_ptr<Simulator> sim_;
    RunOutput out_;
    ecl::Summary sum_;
    bool sumValid_ = false;
    double t0_ = 0;
};

}  // namespace api
}  // namespace ressim
