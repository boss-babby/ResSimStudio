// Multi-run studies: parallel batch runner, sensitivity analysis (tornado + Latin-hypercube Monte Carlo)
// and history matching (manual parameter sets, Levenberg-Marquardt, ES-MDA).
#pragma once
#include "ressim/simulator.h"
#include <vector>
#include <string>
#include <functional>
#include <atomic>

namespace ressim {

namespace ecl { struct Summary; }   // ecl.h

// ---------------------------------------------------------------- uncertain parameters
const char* toString(ParamKind k);
bool paramIsMultiplier(ParamKind k);
bool paramUsesBox(ParamKind k);

// StudyParam / ParamKind live in case.h (persisted with the case)

std::string paramLabel(const StudyParam& p, const Case& c);

// value <-> normalised coordinate u in [0,1]
double paramToUnit(const StudyParam& p, double v);
double paramFromUnit(const StudyParam& p, double u);
// apply parameter values (one per entry of params, inactive ones use base) to a copy of the case
Case applyParams(const Case& base, const std::vector<StudyParam>& params, const std::vector<double>& values);

// ---------------------------------------------------------------- runs
struct RunResult {
    bool ok = false;
    std::string error;
    std::vector<SummaryPoint> summary;
    RunStats stats;
    double oip0 = 0;                   // initial oil in place (surface)
};
// Runs all cases on 'workers' threads (0 = hardware concurrency; each run is single-threaded).
std::vector<RunResult> runBatch(const std::vector<Case>& cases, int workers, std::atomic<bool>* cancel = nullptr,
                                const std::function<void(int done, int total)>& progress = {});
RunResult runSingle(const Case& c, std::atomic<bool>* cancel = nullptr);

// ---------------------------------------------------------------- observations and misfit
// value of an observed quantity in a run (linear interpolation in time); NaN if the well is unknown
double simulatedValue(const RunResult& r, const Case& c, ObsQty q, const std::string& well, double t);
std::vector<double> autoSigma(const std::vector<Observation>& obs);       // per observation
// normalised residuals (sim - obs) / sigma; failed run -> large values
std::vector<double> residuals(const RunResult& r, const Case& c, const std::vector<Observation>& obs,
                              const std::vector<double>& sigma);
double misfit(const RunResult& r, const Case& c, const std::vector<Observation>& obs, const std::vector<double>& sigma);
// observations from WCONHIST rates (per step) and synthetic "truth" data with noise
std::vector<Observation> observationsFromHistory(const Case& c);
// observations from an ECLIPSE / OPM summary file (SMSPEC/UNSMRY): the history vectors (WOPRH, WWPRH, WGPRH, WWIRH,
// WBHPH, WWCTH, WGORH, FOPRH, ...) when history = true, else the simulated vectors (another run used as the truth).
// Sampled at the report steps, at least minSpacingDays apart; times are shifted to the start date of the case.
struct SummaryObsOptions { bool history = true; double minSpacingDays = 0; bool skipShut = true; unsigned qtyMask = 0xFFFFFFFFu; };
std::vector<Observation> observationsFromSummary(const ecl::Summary& s, const Case& c, const SummaryObsOptions& opt,
                                                 std::string* note = nullptr);
std::vector<Observation> syntheticObservations(const RunResult& truth, const Case& c, const std::vector<ObsQty>& qtys,
                                               double relNoise, unsigned seed);

// ---------------------------------------------------------------- multi-objective misfit
// Observations grouped into objective terms; misfit of a term = mean squared normalised residual of its observations
// (observation weights are folded into sigma by autoSigma: sigma / sqrt(weight)).
enum class ObjectiveGrouping { Quantity, Well, QuantityWell, PhaseClass };
struct ObjectiveTerm { std::string name; std::vector<int> idx; };
std::vector<ObjectiveTerm> objectiveTerms(const std::vector<Observation>& obs, ObjectiveGrouping g);
std::vector<double> misfitTerms(const RunResult& r, const Case& c, const std::vector<Observation>& obs,
                                const std::vector<double>& sigma, const std::vector<ObjectiveTerm>& terms);
// indices of the non-dominated points (minimisation of every column of f; rows with NaN are skipped)
std::vector<int> paretoFront(const std::vector<std::vector<double>>& f);

// ---------------------------------------------------------------- responses for sensitivity analysis
enum class Response { FOPT, RF, FWPT, FGPT, FWIT, FPR, FWCT, NPV, Misfit };
const char* toString(Response r);
struct Economics { double oilPrice = 70, gasPrice = 3, waterCost = 5, injCost = 2, discount = 0.10; };
double responseValue(const RunResult& r, const Case& c, Response resp, const Economics& eco,
                     const std::vector<Observation>& obs = {}, const std::vector<double>& sigma = {});

// ---------------------------------------------------------------- study callbacks
struct StudyCallbacks {
    std::atomic<bool>* cancel = nullptr;
    std::function<void(const std::string&)> log;
    std::function<void(int done, int total)> progress;          // runs finished in the current batch
    int workers = 0;
    // called (on the study thread) for every simulation a study performs: tag ("base", "tornado", "MC", "LM",
    // "ES-MDA", "design", "proxy"), parameter values (one per StudyParam) and the run - used for the run archive
    std::function<void(const char* tag, const std::vector<double>& values, const RunResult& run)> onRun;
};

// ---- sensitivity: one-at-a-time tornado
struct TornadoBar { int param = 0; double low = 0, high = 0; bool okLow = false, okHigh = false; };
struct TornadoResult {
    double base = 0; bool ok = false;
    std::vector<TornadoBar> bars;                                // sorted by swing (largest first)
    RunResult baseRun;
};
TornadoResult runTornado(const Case& base, const std::vector<StudyParam>& params, Response resp, const Economics& eco,
                         const StudyCallbacks& cb);

// ---- sensitivity: Latin-hypercube Monte Carlo
struct MonteCarloResult {
    std::vector<std::vector<double>> samples;                    // [run][param] values
    std::vector<double> response;                                // per run (NaN if failed)
    std::vector<double> spearman;                                // per param
    double p10 = 0, p50 = 0, p90 = 0, mean = 0;
    std::vector<double> times;                                   // fan chart of FOPT
    std::vector<double> fanP10, fanP50, fanP90;
    int failed = 0;
};
MonteCarloResult runMonteCarlo(const Case& base, const std::vector<StudyParam>& params, int samples, unsigned seed,
                               Response resp, const Economics& eco, const StudyCallbacks& cb);

// ---- history matching
struct HMIteration { int iter = 0; double misfit = 0; double lambda = 0; std::vector<double> values; double spread = 0; };
struct HMResult {
    bool ok = false;
    std::vector<double> best;                                    // best parameter values
    double bestMisfit = 0, initialMisfit = 0;
    std::vector<HMIteration> history;
    RunResult bestRun, initialRun;
    int runs = 0;
    // ES-MDA: final ensemble
    std::vector<std::vector<double>> ensemble;                   // [member][param]
    std::vector<double> ensembleMisfit;
    std::vector<RunResult> ensembleRuns;
};
struct LMOptions { int maxIter = 12; double fdStep = 0.04; double lambda0 = 0.01; double tol = 1e-3; };
HMResult runLevenbergMarquardt(const Case& base, const std::vector<StudyParam>& params, const LMOptions& opt,
                               const StudyCallbacks& cb);
struct ESMDAOptions { int ensemble = 32; int assimilations = 4; double priorStd = 1.0; unsigned seed = 1; };
HMResult runESMDA(const Case& base, const std::vector<StudyParam>& params, const ESMDAOptions& opt, const StudyCallbacks& cb);

// ---- experimental design (values per run, one per StudyParam; inactive parameters stay at their base value)
enum class DesignKind { LatinHypercube, FullFactorial2, CentralComposite, PlackettBurman, Random };
const char* toString(DesignKind k);
std::vector<std::vector<double>> designPoints(const std::vector<StudyParam>& params, DesignKind kind, int samples, unsigned seed);
// runs the given parameter sets in parallel (reports every run through cb.onRun with the given tag)
std::vector<RunResult> runDesign(const Case& base, const std::vector<StudyParam>& params,
                                 const std::vector<std::vector<double>>& values, const char* tag, const StudyCallbacks& cb);

// ---- response-surface proxy (polynomial regression in normalised parameter space)
struct ProxyModel {
    bool ok = false, quadratic = true;
    std::vector<int> act;                  // active parameter indices
    std::vector<double> coef;              // 1, u_k, u_k*u_l (k<=l)
    double r2 = 0, q2 = 0, rmse = 0;       // fit R^2, leave-one-out Q^2, RMS error
    int samples = 0;
    std::vector<double> importance;        // per active parameter: |linear| + |quadratic| effect over [0,1]
};
ProxyModel fitProxy(const std::vector<StudyParam>& params, const std::vector<std::vector<double>>& values,
                    const std::vector<double>& y, bool quadratic);
double proxyPredict(const ProxyModel& m, const std::vector<StudyParam>& params, const std::vector<double>& values);
// minimise (or maximise) the proxy by dense random search + coordinate refinement; returns parameter values
std::vector<double> proxyOptimum(const ProxyModel& m, const std::vector<StudyParam>& params, bool maximise, unsigned seed,
                                 double* predicted = nullptr);

} // namespace ressim
