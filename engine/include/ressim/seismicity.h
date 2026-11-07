// C20: rate-and-state (Dieterich 1994) induced-seismicity forecasting and the convex minimum-seismicity injection designer.
//
// Theory (research/CYCLE1_REPORT.md, research/papers/P1_convex_seismicity):
//   R/r = exp((S - tau_r N/r)/a),  a = A*sigma (MPa),  t_a = a/tau_r,  S = tau_r t + dS(t)
//   N(T)/r = t_a ln(1 + I/t_a),   I = int_0^T exp(s/t_a + dS(s)/a) ds          (Heimisson & Segall 2018, exact)
//   excess events E = N - rT = r t_a [ln(1 + I/t_a) - T/t_a]
// dS is integrated exactly for piecewise-linear dS (exp of a linear function per segment).
// When dS = G q is linear in the injection rates, E(q) is a log-sum-exp of affine maps -> CONVEX: the designer finds the
// global optimum (FISTA with exact projection onto bounds + volume equalities).
#pragma once
#include "ressim/case.h"
#include <vector>
#include <string>
#include <atomic>
#include <functional>

namespace ressim {
namespace seis {

struct Patch {                       // a fault patch represented by one grid cell
    std::string name = "F1";
    int i = 0, j = 0, k = 0;         // zero-based global cell
    double friction = 0.6;           // f in dCFS = f * dp * cfsFactor  (pore-pressure Coulomb change)
    double cfsFactor = 1.0;          // stress-path factor (1 = pore pressure only; <1 poroelastic stabilisation)
    double aSigma = 0.05;            // A*sigma_eff (MPa)
    double taYears = 100.0;          // aftershock duration t_a = A sigma / tau_r (years)
    double rate = 1.0;               // background rate r (events / year)
    double weight = 1.0;             // weight in the design objective
};

struct Forecast {                    // per patch
    std::vector<double> t;           // days
    std::vector<double> dS;          // Coulomb stress change (MPa)
    std::vector<double> R;           // seismicity rate (events / year)
    std::vector<double> N;           // cumulative events (incl. background)
    double excess = 0;               // expected induced (excess) events at the last time
    double excessAtShutIn = 0;       // excess events at tShutIn (if inside the window, else = excess)
    double peakR = 0, tPeak = 0;     // events / year, days
    double maxDS = 0;                // MPa
};
// Exact Dieterich / Heimisson-Segall forecast for piecewise-linear dS(t) (t ascending, days; dS MPa)
Forecast forecast(const std::vector<double>& tDays, const std::vector<double>& dSMPa, const Patch& p, double tShutIn = -1);

// ---- convex design on a linear response model  dS_p(t_i) = dS0_p(t_i) + sum_v G_p[i][v] x_v
enum class Objective { Events, Peak, Blend };
const char* toString(Objective o);
struct LinearModel {
    std::vector<double> t;                                 // eval times (days)
    std::vector<Patch> patches;
    std::vector<std::vector<double>> dS0;                  // [patch][i]  stress with all design controls at zero
    std::vector<std::vector<double>> G;                    // [patch][i * nv + v]
    int nv = 0;
};
struct Constraints {
    std::vector<double> lb, ub;                            // per variable
    std::vector<int> group;                                // per variable: equality group index
    std::vector<double> weight;                            // per variable: coefficient in its group's sum (e.g. dt)
    std::vector<double> target;                            // per group: sum_v weight_v x_v = target_g
};
struct OptOptions { Objective obj = Objective::Events; double blend = 0.5; double peakSharpness = 40; int maxIter = 4000;
                    double tol = 1e-9; };
struct OptResult { std::vector<double> x; double F = 0, F0 = 0; double events = 0, events0 = 0; double peak = 0, peak0 = 0;
                   int iters = 0; bool converged = false; double kkt = 0; };
// total weighted excess events  sum_p w_p E_p(x)  and its gradient
double modelEvents(const LinearModel& m, const std::vector<double>& x, std::vector<double>* grad = nullptr);
// smooth maximum of dS_p/a_p over time and patches (log-sum-exp with sharpness beta) and its gradient
double modelPeak(const LinearModel& m, const std::vector<double>& x, double beta, std::vector<double>* grad = nullptr);
std::vector<double> modelStress(const LinearModel& m, const std::vector<double>& x, int patch);
// projection onto {lb <= x <= ub, sum_{v in g} w_v x_v = target_g}
void project(const Constraints& c, std::vector<double>& x);
OptResult optimize(const LinearModel& m, const Constraints& c, const std::vector<double>& x0, const OptOptions& o);

// ---- engine-driven designer: Green's functions from the ResSim simulator, convex optimisation, nonlinear verification
struct DesignSetup {
    std::vector<std::string> wells;        // injectors whose schedules are designed
    std::vector<Patch> patches;
    int nCtrl = 10;                        // control intervals in the injection window
    double injDays = 30, postDays = 60;    // injection window and post-shut-in horizon
    int subSteps = 2, postSteps = 16;      // report steps per control interval / in the post window (geometric)
    std::vector<double> wellVolume;        // per well surface volume (case units); empty -> totalVolume over all wells
    double totalVolume = 0;
    double maxRateFactor = 3.0;            // per-well rate bound = factor * mean rate
    double lastRateFactor = 0;             // C22 (H1f): rate bound of the LAST control interval / mean rate (0 = same as
                                           // maxRateFactor). The event-count optimum ends with a max-rate burst that raises
                                           // the stress at shut-in; a bound <= 1 removes it (post-shut-in safety)
    int correctionIters = 3;               // sequential convex programming: add the (simulated - predicted) stress of the
                                           // last verified schedule to the model, re-optimise, re-verify (1 run each)
    double correctionTol = 0.01;           // stop when the verified excess events agree with the model to this fraction
    OptOptions opt;
};
struct DesignOutput {
    bool ok = false; std::string error;
    std::vector<std::string> logLines;
    std::vector<double> ctrlStart;                         // days, nCtrl + 1 entries
    std::vector<std::vector<double>> qConst, qOpt;         // [well][ctrl] surface rate (case units)
    std::vector<Forecast> predConst, predOpt;              // per patch, linear model
    std::vector<Forecast> simConst, simOpt;                // per patch, nonlinear forward runs
    double eventsConst = 0, eventsOpt = 0, simEventsConst = 0, simEventsOpt = 0;
    double peakConst = 0, peakOpt = 0;                     // max R/r over patches (linear model)
    double linearityError = 0;                             // first-pass max relative |sim - pred| of excess events
    double finalModelError = 0;                            // |sim - pred| / |pred| of the returned schedule after correction
    std::vector<double> corrPred, corrSim;                 // per correction pass: predicted / simulated excess events
    int runs = 0, iters = 0, corrections = 0; double seconds = 0;
    // C22 (research H1f, F12): design potential of the case from the closed form - Pi = max dCFS/(A sigma) of the constant
    // schedule, theta = stress-response time / injection window (63 % rise of the constant-rate stress) on the governing patch
    double Pi = 0, theta = 0, gainTheory = 0, gainTheoryNoBurst = 0;
    bool endsWithBurst = false;            // the designed schedule ends at its rate bound (shut-in burst)
    bool gainTheoryIsBound = false;        // outside the closed-form range: gainTheory is a lower bound of the first-order optimum
    bool farField = false;                 // governing patch: stress keeps rising after shut-in (> 5 %) -> closed form not applicable
};
// C22 (H1f): closed-form gain of the optimal schedule over constant rate for one exponential stress kernel
// (Dieterich, background-dominated): bang-singular-bang = max-rate ramp, constant-stress plateau, max-rate burst at
// shut-in. Pi = max stress / (A sigma) of the constant schedule, theta = kernel time / injection time, uMax = max/mean rate.
struct GainEstimate {
    double gain = 0, gainNoBurst = 0, burstFraction = 0 /* burst volume / total */, sPlateau = 1, sEnd = 1;
    bool exact = false;   // inside the validated range (theta <= 0.1, or uMax <= 3 and theta <= 0.4); otherwise a lower bound (a coast arc appears)
};
GainEstimate scheduleGain(double Pi, double theta, double uMax);

// C23 (research Cycle 5, finding F14): which injection "kernel" do observed stage-wise event counts follow? Poisson fits
// with one free rate factor per test (site / interval effect) and the stage measure as the kernel:
//   constant, duration, volume V, hydraulic energy E = int p q dt, V q^beta, V (1 + V_cum)^gamma, V q^beta (1 + V_cum)^gamma
// (V_cum = volume injected earlier in the same test = growing fracture). AIC ranks them; the profile-likelihood 95 %
// interval of beta in the joint model says whether a lower rate reduces events at fixed volume - the premise of every
// rate-scheduling design (research/CYCLE5_REPORT.md: at Aspo HRL it does not, beta = -0.1 [-1.25, 0.75]).
struct StageObs {
    int test = 0;                // test / interval id (own rate factor)
    double volume = 0;           // injected volume of the stage (any unit)
    double energy = 0;           // hydraulic energy (any unit; <= 0: not used)
    double rate = 0;             // mean injection rate (any unit; <= 0: rate models skipped)
    double duration = 0;         // s (<= 0: not used)
    double cumVolume = -1;       // volume injected before this stage in the same test; < 0: accumulated in the given order
    int count = 0;               // observed events
};
struct KernelModel { std::string name; double logL = 0, aic = 0; int params = 0; double exponent = 0, exponent2 = 0; };
struct KernelFit {
    std::vector<KernelModel> models;   // sorted by AIC
    std::string best;
    double rateExponent = 0, rateLo = 0, rateHi = 0;   // beta of the joint model + profile 95 % interval
    double rateExponentAlone = 0;                      // beta of V q^beta alone (confounded by stage order in step-rate protocols)
    double cumExponent = 0;                            // gamma of the joint model
    double rateCumCorrelation = 0;                     // correlation of log q and log(1 + V_cum) over the stages
    std::string verdict;
};
KernelFit fitStageKernels(const std::vector<StageObs>& stages);
// a runnable case with the designed schedule (q[well][ctrl]); other wells keep their first control record
Case designCase(const Case& base, const DesignSetup& s, const std::vector<std::vector<double>>& q);
// SCHEDULE section text (WCONINJE + TSTEP per control interval) for export to ECLIPSE / OPM decks
std::string designScheduleDeck(const Case& base, const DesignSetup& s, const std::vector<std::vector<double>>& q);
DesignOutput designInjection(const Case& base, const DesignSetup& s, std::atomic<bool>* cancel = nullptr,
                             const std::function<void(const std::string&)>& log = {});

} // namespace seis
} // namespace ressim
