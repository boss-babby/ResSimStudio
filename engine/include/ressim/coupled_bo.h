// C27: headless coupled black-oil run (moved from the CLI so that studies, the CLI and scripts share one runner).
// Same physics and order as the GUI runner: FE geomechanics (GEOMECH, explicit or iterative coupling) -> fault
// mechanics (FAULTMECH / FAULTSLIP) -> DDM fracture apertures (HFGROW ... DDMWIDTH), one report-step lag for the
// explicit couplings. Results are the flow summary plus scalar key performance indicators (seismicity, uplift, ...).
#pragma once
#include "ressim/case.h"
#include <atomic>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace ressim {
class Simulator;
namespace faultmech { class Monitor; }
namespace coupledbo {

// true when the case asks for coupled mechanics the plain Simulator does not run on its own
bool needed(const Case& c);

struct Output {
    std::vector<std::pair<std::string, double>> kpi;   // e.g. fault_F1_expected_events, surface_uplift_m, ...
    std::vector<std::string> report;                   // human-readable end-of-run lines
    int iterSteps = 0, iterPasses = 0;                 // iterative geomechanics statistics
    bool finished = false;
};

// Runs the whole schedule on a simulator that is set up and initialised for 'c'. log: progress messages (may be empty);
// progress: per report step callback (time in days). Returns false when cancelled / the run stopped early.
bool run(Simulator& sim, const Case& c, Output& out, const std::function<void(const std::string&)>& log = {},
         std::atomic<bool>* cancel = nullptr, const std::function<void(double)>& progress = {});

// C28: end-of-run fault lines + KPIs of a fault monitor (shared by the black-oil and the compositional / thermal runners)
void faultReport(const faultmech::Monitor& fm, Output& out);

// aggregate KPIs used by the study responses (NaN when the quantity does not exist in the run)
double kpiSum(const std::vector<std::pair<std::string, double>>& kpi, const std::string& suffix);   // sum over faults
double kpiMax(const std::vector<std::pair<std::string, double>>& kpi, const std::string& suffix);   // max over faults

}  // namespace coupledbo
}  // namespace ressim
