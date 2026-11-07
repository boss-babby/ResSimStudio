// C40 geomodel data analysis (Petrel "Data analysis" / "Well tie" level):
//   * log calculator   - vectorised expressions over the curves of a well log
//   * variogram fit    - weighted least-squares fit of a variogram model to an experimental variogram
//   * regression       - linear / log-linear regression for cross-plots (poro-perm transforms)
//   * well tie         - synthetic seismogram from sonic + density, time-depth relation, correlation with the seismic
//                        trace at the well (bulk shift search) and a velocity model calibrated on the tie
#pragma once
#include "ressim/geomodel.h"
#include "ressim/geostat.h"
#include <map>
#include <string>
#include <vector>

namespace ressim {
namespace geo {

// ---------------------------------------------------------------- log calculator
// Expression over named sample vectors (all the same length). Operators + - * / ^, comparisons < <= > >= == !=,
// && || ! ; functions abs sqrt ln log10 exp pow(a,b) min(a,b) max(a,b) clamp(x,lo,hi) if(c,a,b) isnan(x) and constants
// pi, nan. Names are case-insensitive. NaN (null samples) propagates, also through comparisons and if(); && / || only
// short-circuit on a definite operand (0 && null = 0, 1 || null = 1); isnan(x) tests for nulls explicitly.
// Returns false with a message (position of the error) on syntax errors or unknown names.
bool evalExpression(const std::string& expr, const std::map<std::string, const std::vector<double>*>& vars, size_t n,
                    std::vector<double>& out, std::string& err);
// convenience: all curves of a log plus MD (index) and TVD (trajectory of 'w', if given) as variables;
// nulls (log.nullValue) become NaN. The result is a new curve 'mnem' added to (or replacing it in) the log.
bool addCalculatedCurve(WellLog& log, const GeoWell* w, const std::string& mnem, const std::string& unit, const std::string& expr,
                        std::string& err);
// names usable in expressions for this log (curves + MD + TVD)
std::vector<std::string> expressionNames(const WellLog& log);

// ---------------------------------------------------------------- C42 surface / property calculators, statistics
// Name of a surface / property as an expression variable: upper case, every character outside [A-Z0-9_] -> '_',
// a leading digit gets a '_' prefix ("Top Ness 2" -> TOP_NESS_2).
std::string exprName(const std::string& name);
// Surface calculator: every valid surface in 'hz' is resampled (bilinear, NaN outside) onto the lattice of hz[onto] and
// available under exprName(name); X and Y are the node coordinates. 'out' gets the lattice of hz[onto] and the result.
bool surfaceCalculator(const std::vector<Surface>& hz, int onto, const std::string& expr, Surface& out, std::string& err);
// Descriptive statistics over the finite values of v
struct Stats { size_t n = 0, nNull = 0; double min = 0, max = 0, mean = 0, sd = 0, p10 = 0, p50 = 0, p90 = 0, sum = 0; };
Stats describe(const std::vector<double>& v);

// petrophysical presets (expressions written with the given curve names)
struct CalcPreset { std::string name, mnem, unit, expr, help; };
std::vector<CalcPreset> calcPresets();

// ---------------------------------------------------------------- regression
struct Regression { double a = 0, b = 0, r2 = 0; int n = 0; bool logY = false; };   // y = a + b x  (logY: log10 y = a + b x)
Regression regress(const std::vector<double>& x, const std::vector<double>& y, bool logY);

// ---------------------------------------------------------------- variogram fitting
struct VarioFit { double nugget = 0, sill = 1, range = 1, rmse = 0; bool ok = false; };
// unit-sill structure function of the model at h / range (practical range for exponential and Gaussian)
double structure(geostat::VarioModel m, double hr);
// weighted least squares (Cressie weights N(h) / gamma_model(h)^2, iterated) over the range by a log grid + golden
// refinement, nugget and sill by non-negative linear least squares for each trial range. fixNugget >= 0 keeps it fixed.
VarioFit fitVariogram(const geostat::ExpVario& ev, geostat::VarioModel model, double fixNugget = -1, double maxRange = 0);

// ---------------------------------------------------------------- well tie (synthetic seismogram)
struct TieInput {
    std::string dt = "DT", rhob = "RHOB";       // sonic (slowness) and bulk density curves
    int dtUnit = 0;                              // 0 automatic, 1 us/ft, 2 us/m
    double rhoConst = 2.3;                       // used when there is no density curve (g/cc)
    double freq = 25;                            // Ricker peak frequency (Hz)
    double dtMs = 1.0;                           // sampling of the synthetic (ms)
    double maxShift = 80;                        // bulk shift search window (+- ms)
    bool field = false;                          // project lengths in feet
};
struct TieResult {
    bool ok = false;
    std::string msg;
    std::vector<double> twt, ai, rc, synth;      // on the synthetic time axis (TWT ms, after the best shift)
    std::vector<double> traceT, trace;           // seismic trace at the well (cube time axis)
    std::vector<double> tdZ, tdT;                // time-depth pairs (TVD, TWT ms) from the sonic integration (shifted)
    double shift = 0, corr0 = 0, corr = 0;       // best bulk shift (ms), correlation before / after
    VelocityModel vmFit;                         // linear-gradient velocity model fitted to the time-depth pairs
    double vmRmsMs = 0;                          // misfit of the fitted model (ms)
};
TieResult wellTie(const GeoWell& w, const WellLog& log, const SeismicCube& cube, const VelocityModel& vm, const TieInput& in);
// time-depth fit for v(z) = v0 + k (z - datum)
VelocityModel fitVelocityModel(const std::vector<double>& z, const std::vector<double>& twtMs, double datum, double* rmsMs = nullptr);

}  // namespace geo
}  // namespace ressim
