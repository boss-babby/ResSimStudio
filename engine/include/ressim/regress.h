// C33 deck regression harness: runs a list of decks with the ResSim engine and compares the summary vectors with
// reference results (OPM Flow / ECLIPSE SMSPEC + UNSMRY, e.g. the opm-simulation-reference files of OPM/opm-tests),
// tracks wall-clock time against the previous run (performance regression) and writes HTML / DOCX / JUnit reports.
//
// Manifest (text, one directive per line, '#' comments, paths relative to the manifest):
//   tol 0.05                                   default tolerance: relative RMS deviation (RMS / max|ref|) per vector
//   keys FOPR,FOPT,FWCT,FGOR,FPR,WBHP:*        default vectors ('*' = every well / group present in both)
//   case NAME deck=spe1/SPE1CASE1.DATA ref=spe1/opm-simulation-reference/flow/SPE1CASE1.SMSPEC [tol=0.1] [keys=...] [group=spe1]
//   skip NAME reason text                      listed in the report, not run
#pragma once
#include "ressim/ecl.h"
#include "ressim/report.h"
#include <atomic>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace ressim {
namespace regress {

struct CaseSpec {
    std::string name, deck, reference, group;
    std::vector<std::string> keys;    // empty = manifest default
    double tol = -1;                  // < 0 = manifest default
    std::string skipReason;           // non-empty: not run
};
struct Manifest {
    std::string dir;                  // directory of the manifest (paths are relative to it)
    double tol = 0.05;
    std::vector<std::string> keys = {"FOPR", "FOPT", "FWPR", "FWPT", "FGPR", "FGPT", "FWIR", "FWIT", "FGIR", "FGIT", "FPR", "FWCT", "FGOR", "WBHP:*", "WOPR:*"};
    std::vector<CaseSpec> cases;
};
bool readManifest(const std::string& path, Manifest& m, std::string& err);

struct VectorResult { std::string key; double relRms = 0, refFinal = 0, cmpFinal = 0; bool pass = true; };
struct CaseResult {
    std::string name, group, deck;
    std::string status;               // PASS, FAIL, ERROR (did not run / crashed), SKIP (not supported / listed as skip)
    std::string message, engine;
    double seconds = 0, prevSeconds = -1, tol = 0.05;
    bool slow = false;                // > slowFactor x the previous wall time
    int steps = 0, newton = 0;
    std::vector<VectorResult> vectors;
    double worst = 0;
    std::string worstKey;
    ecl::Summary run, ref;            // kept for the charts
};
struct Options {
    std::string outDir = "regression_out";
    std::string history;              // TSV with previous runs (default outDir/history.tsv)
    double slowFactor = 1.10;
    std::string filter;               // run only cases whose name contains this text
    double tolOverride = -1;
    double maxSeconds = 0;            // > 0: a run is cancelled after this wall time (status ERROR "timeout")
    int chartCases = 12;              // cases with overlay charts in the report (failures first)
    std::function<void(const std::string&)> log;
    const std::atomic<bool>* cancel = nullptr;   // C36: set to true to stop (the running case is aborted, the rest is not run)
};
std::vector<CaseResult> runAll(const Manifest& m, const Options& o);
// compares one run with its reference (exposed for tests)
void compare(CaseResult& r, const std::vector<std::string>& keys, double tol);
report::Doc reportDoc(const std::vector<CaseResult>& r, const Options& o);
bool writeJUnit(const std::string& path, const std::vector<CaseResult>& r, std::string& err);
// history: date \t name \t status \t seconds \t worst
std::map<std::string, double> readHistory(const std::string& path);
bool appendHistory(const std::string& path, const std::vector<CaseResult>& r);
// writes report.html, report.docx, junit.xml, appends history; returns the number of FAIL + ERROR cases (-1: I/O error)
int writeOutputs(const std::vector<CaseResult>& r, const Options& o, std::string& err);

}  // namespace regress
}  // namespace ressim
