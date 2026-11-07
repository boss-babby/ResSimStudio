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
    // C27: scalar results of the coupled mechanics (fault_<name>_expected_events / _max_mw / _max_dcfs / _max_slip_m,
    // surface_uplift_m, max_failure_index, hf_max_aperture_m ...) - empty for plain flow runs
    std::vector<std::pair<std::string, double>> kpi;
    // C49: well names of the per-well rates in summary[].wells for the compositional / CO2STORE / thermal engines (the
    // black-oil path takes them from the Simulator); empty = no per-well data
    std::vector<std::string> wellNames;
};
// Runs all cases on 'workers' threads (0 = hardware concurrency; each run is single-threaded).
std::vector<RunResult> runBatch(const std::vector<Case>& cases, int workers, std::atomic<bool>* cancel = nullptr,
                                const std::function<void(int done, int total)>& progress = {});
RunResult runSingle(const Case& c, std::atomic<bool>* cancel = nullptr);
// C16: optional external batch runner (MPI cluster / job scheduler). When set, runBatch delegates to it; it must
// return one result per case (failed runs ok = false).
using BatchRunner = std::function<std::vector<RunResult>(const std::vector<Case>&, std::atomic<bool>*, const std::function<void(int, int)>&)>;
void setBatchRunner(BatchRunner r);
bool hasBatchRunner();

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
// C27: SeisEvents .. MaxFailure need the coupled mechanics in the case (FAULTMECH / GEOMECH), else NaN
enum class Response { FOPT, RF, FWPT, FGPT, FWIT, FPR, FWCT, NPV, Misfit, SeisEvents, MaxMw, MaxDCFS, Uplift, MaxFailure };
constexpr int kNumResponses = 14;
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
// one accepted iteration of a history-matching run. C13: LM / adjoint record the gradient norm, cumulative
// simulations and wall time, the gradient of the misfit with respect to the normalised parameters (per parameter,
// NaN = inactive) and the run itself (summary vectors only) so every iterate's responses can be inspected.
struct HMIteration { int iter = 0; double misfit = 0; double lambda = 0; std::vector<double> values; double spread = 0;
                     double gradNorm = 0; int runs = 0; double seconds = 0; std::vector<double> grad; RunResult run; };
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
// Jacobian / gradient source: FiniteDifference = one extra simulation per parameter; Analytic = one recorded forward
// run + discrete forward sensitivities (exact Jacobian of the frozen-coefficient discrete model, no extra simulations)
enum class GradientMethod { FiniteDifference, Analytic };
const char* toString(GradientMethod g);
struct LMOptions { int maxIter = 12; double fdStep = 0.04; double lambda0 = 0.01; double tol = 1e-3;
                   GradientMethod gradient = GradientMethod::FiniteDifference; double resStep = 1e-4; };
HMResult runLevenbergMarquardt(const Case& base, const std::vector<StudyParam>& params, const LMOptions& opt,
                               const StudyCallbacks& cb);
// ---- adjoint-gradient history matching: projected L-BFGS (box [0,1] in normalised parameters) driven by the discrete
// adjoint gradient of the misfit - one forward + one backward pass per iteration, independent of the parameter count
struct AdjointHMOptions { int maxIter = 25; double tol = 1e-5; int memory = 6; double resStep = 1e-4; double maxStep = 0.25; };
HMResult runAdjointLBFGS(const Case& base, const std::vector<StudyParam>& params, const AdjointHMOptions& opt,
                         const StudyCallbacks& cb);
// gradient of the misfit with respect to the normalised active parameters (u in [0,1]) at their base values, by any
// ---- C11: per-cell field history matching (JutulDarcy-style): projected L-BFGS on ln-multipliers of the permeability
// (all directions) and / or porosity of EVERY active cell, gradients from one adjoint pass with analytic dR/dtheta.
// Objective = misfit / nd + regWeight * mean over grid faces of (x_a - x_b)^2 (smoothness of the log-multipliers).
struct FieldHMOptions { bool perm = true, poro = false; int maxIter = 20; double regWeight = 1e-3; double lo = 0.1, hi = 10;
                        double maxLogStep = 0.7; double tol = 1e-6; int memory = 8; };
struct FieldHMResult {
    bool ok = false; std::string error;
    double initialMisfit = 0, bestMisfit = 0, bestObjective = 0;
    int runs = 0, iterations = 0;
    std::vector<double> multK, multPhi;                          // per global cell (1 = unchanged / inactive)
    std::vector<double> history;                                 // misfit per accepted iteration
    Case best;                                                   // the base case with the multipliers applied
    RunResult bestRun, initialRun;
    double secAdjoint = 0, secForward = 0;
    // C13: full record of every accepted iterate (iteration 0 = prior). ln-multiplier snapshots and the data-misfit
    // gradient are stored per GLOBAL cell as float (0 / NaN for inactive cells); 'run' keeps the summary vectors.
    struct Iter {
        int iter = 0, runs = 0;
        double misfit = 0, objective = 0, reg = 0, gradNorm = 0, step = 0, seconds = 0;
        double kMin = 1, kMax = 1, kMean = 1, phiMin = 1, phiMax = 1, phiMean = 1;   // multiplier statistics (active cells)
        std::vector<float> lnK, lnPhi;                               // per global cell (empty when not matched)
        RunResult run;
    };
    std::vector<Iter> iters;
    std::vector<float> gradK0, gradPhi0, gradK, gradPhi;         // d(misfit)/d ln m per global cell: prior and final
    int nx = 0, ny = 0, nz = 0;                                  // grid dimensions of the per-cell arrays
};
FieldHMResult runAdjointFieldHM(const Case& base, const FieldHMOptions& opt, const StudyCallbacks& cb);
Case applyFieldMultipliers(const Case& base, const std::vector<double>& multK, const std::vector<double>& multPhi);

// combination of methods (mask: 1 = adjoint, 2 = forward sensitivities, 4 = central finite differences of full runs)
struct GradientCheck {
    bool ok = false; std::string error;
    std::vector<int> act;                                          // active parameter indices
    std::vector<double> adjoint, forward, fd;                      // per active parameter (empty when not requested)
    double misfit = 0; double secAdjoint = 0, secForward = 0, secFd = 0; int fdRuns = 0;
    double replayDev = 0, maxLinRes = 0;
};
GradientCheck misfitGradient(const Case& base, const std::vector<StudyParam>& params, int methods, double fdStep,
                             const StudyCallbacks& cb, double resStep = 1e-4);
struct ESMDAOptions { int ensemble = 32; int assimilations = 4; double priorStd = 1.0; unsigned seed = 1; };
HMResult runESMDA(const Case& base, const std::vector<StudyParam>& params, const ESMDAOptions& opt, const StudyCallbacks& cb);

// ---- C27: iTOUGH2-style linearised inversion statistics at a parameter set (normally the history-match optimum).
// Jacobian of the weighted residuals r_i = (sim_i - obs_i) / sigma_i by forward differences in the search space
// (log10 for log parameters), then (Finsterle, iTOUGH2; Hill & Tiedeman 2007):
//   s0^2 = r^T r / (m - n),  C_pp = s0^2 (J^T J)^-1,  correlation, eigen-analysis, composite scaled sensitivities,
//   leverage h_ii = [J (J^T J)^-1 J^T]_ii, studentized residuals w_i = r_i / (s0 sqrt(1 - h_ii)) (|w| > 1.96 -> outlier),
//   prediction standard deviation of every observed quantity sigma_i s0 sqrt(h_ii), objective contribution per data set.
struct HMStatsOptions {
    double fdStep = 0.02;                          // forward step as a fraction of the search range
    // C28 FOSM prediction uncertainty + data worth (iTOUGH2): predictions are study responses of the same runs (no extra
    // simulations); candidate data = well / field vectors not yet observed, at the observation times, automatic sigma
    std::vector<Response> predict{Response::FOPT, Response::NPV};
    Economics eco;
    bool dataWorth = true;
};
struct HMStatistics {
    bool ok = false; std::string error;
    int nObs = 0, nPar = 0, runs = 0;
    std::vector<int> param;                        // index into the StudyParam list of every analysed (active) parameter
    std::vector<std::string> name;
    std::vector<double> value;                     // physical value at the analysis point
    std::vector<bool> logSpace;                    // search coordinate log10(value) (else value)
    double objective = 0, s02 = 1;                 // weighted sum of squared residuals, estimated error variance
    std::vector<double> stdDev;                    // a-posteriori standard deviation in the search coordinate
    std::vector<double> ci95lo, ci95hi;            // 95 % confidence interval in physical units (linearised)
    std::vector<std::vector<double>> cov, corr;    // parameter covariance / correlation (search coordinates)
    std::vector<double> eigVal;                    // covariance eigenvalues (descending: worst-determined direction first)
    std::vector<std::vector<double>> eigVec;       // [k][param]
    double condition = 0;                          // condition number of J^T J
    std::vector<double> css;                       // composite scaled sensitivity per parameter
    std::vector<double> residual, studentized, leverage, predStd;   // per observation
    double resMean = 0, resStd = 0, resSkew = 0, resKurt = 0;
    std::vector<int> outliers;                     // observation indices with |studentized| > 1.96
    std::vector<std::pair<std::string, double>> contribution;       // objective share [%] per data set (quantity + well)
    // C28 FOSM: linearised prediction uncertainty sigma_z^2 = g^T C_pp g (g = dz/dx from the Jacobian runs)
    struct Prediction { std::string name; double value = 0, stdDev = 0, lo = 0, hi = 0; };
    std::vector<Prediction> pred;
    // data worth: change of every prediction's standard deviation [%] when an observed data set is removed (existing,
    // > 0 = the set was worth that much) or a candidate data set is added (< 0 = expected reduction)
    struct DataWorth { std::string set; bool existing = true; int points = 0; std::vector<double> dStdPct; };
    std::vector<DataWorth> worth;
};
HMStatistics hmStatistics(const Case& base, const std::vector<StudyParam>& params, const std::vector<double>& values,
                          const HMStatsOptions& opt, const StudyCallbacks& cb);

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
    // C28 Gaussian-process (kriging) proxy - CMOST-style ML proxy: ARD squared-exponential kernel in the unit cube,
    // hyper-parameters by maximum marginal likelihood, exact leave-one-out Q^2, predictive standard deviation
    int kind = 0;                          // 0 polynomial regression, 1 Gaussian process
    std::vector<std::vector<double>> gpX;  // training points (unit cube, active parameters)
    std::vector<double> gpAlpha, gpL;      // K^-1 (y - mean) / scale, Cholesky factor of K (row-major lower)
    std::vector<double> gpEll;             // length scale per active parameter (unit cube)
    double gpSf2 = 1, gpSn2 = 1e-6, gpMean = 0, gpScale = 1, gpLogLik = 0;
};
// C28: Gaussian-process proxy (same inputs as fitProxy)
ProxyModel fitProxyGP(const std::vector<StudyParam>& params, const std::vector<std::vector<double>>& values, const std::vector<double>& y);
// predictive standard deviation of a Gaussian-process proxy (NaN for polynomial proxies)
double proxyStd(const ProxyModel& m, const std::vector<StudyParam>& params, const std::vector<double>& values);
ProxyModel fitProxy(const std::vector<StudyParam>& params, const std::vector<std::vector<double>>& values,
                    const std::vector<double>& y, bool quadratic);
double proxyPredict(const ProxyModel& m, const std::vector<StudyParam>& params, const std::vector<double>& values);
// minimise (or maximise) the proxy by dense random search + coordinate refinement; returns parameter values
std::vector<double> proxyOptimum(const ProxyModel& m, const std::vector<StudyParam>& params, bool maximise, unsigned seed,
                                 double* predicted = nullptr);

} // namespace ressim
