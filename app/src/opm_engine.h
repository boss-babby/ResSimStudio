// OPM Flow as an external simulation engine: locate flow(.exe), export the case as an ECLIPSE deck, run it as a child
// process with live progress from its standard output, then read the ECLIPSE result files back.
// No ImGui / Direct3D dependency (also builds natively on Linux for testing).
#pragma once
#include "ressim/case.h"
#include <atomic>
#include <functional>
#include <string>
#include <vector>

// Options of an external engine run. kind 0 = OPM Flow (flowPath = flow.exe), kind 1 = ECLIPSE (flowPath = eclrun.exe
// or eclipse.exe; eclVersion / eclProgram select the installed version and E100 / E300).
struct OpmOptions {
    int kind = 0;                // 0 OPM Flow, 1 ECLIPSE
    std::string eclVersion;      // ECLIPSE: "-v" of eclrun (empty = latest installed)
    std::string eclProgram = "eclipse";   // ECLIPSE: "eclipse" (E100 black oil) or "e300"
    std::string flowPath;        // flow.exe / flow_blackoil.exe (empty = auto-detect)
    std::string sourceDeck;      // non-empty: run this original deck unchanged (imported, not edited) instead of exporting
    std::string runDir;          // directory for the deck and the OPM output (created)
    int threads = 0;             // --threads-per-process (0 = all hardware threads, at most 8)
    double maxDtDays = 0;        // --solver-max-time-step-in-days (0 = OPM default)
    bool enableTuning = false;   // --enable-tuning=true (OPM ignores TUNING otherwise)
    bool keepWellsOpen = true;   // --shut-unsolvable-wells=false + no minimum step before shutting problematic wells
    std::string extraArgs;       // appended verbatim
};

// Searches (in order): 'hint' (file or folder), FLOW_EXE / OPM_FLOW environment variables, PATH, the installer and
// portable-zip locations (Program Files\OPM Flow\bin, %LOCALAPPDATA%\Programs\OPM*, ...). Empty when not found.
std::string detectFlowExecutable(const std::string& hint = std::string(), std::vector<std::string>* searched = nullptr);
// Problems that would stop flow from starting (e.g. MS-MPI runtime missing on Windows). Empty = none found.
std::string flowPrerequisiteProblem(const std::string& flowPath);
// Modelling differences between ResSim and OPM Flow for this case (shown before an OPM run).
std::vector<std::string> opmCompatibilityNotes(const ressim::Case& c);
// "flow 2024.10" style version string (runs flow --version), empty if it cannot be run.
std::string flowVersion(const std::string& flowPath);

// Child process with merged stdout/stderr delivered line by line.
class ChildProcess {
public:
    ~ChildProcess();
    // returns false (and err) if the process could not be started
    // stdinText (optional) is written to the child's standard input, which is then closed (ECLIPSE asks for the case name)
    bool start(const std::string& exe, const std::vector<std::string>& args, const std::string& cwd, std::string* err,
               const std::string& stdinText = std::string());
    // blocks, calling onLine for every output line; returns the exit code (or -1)
    int pump(const std::function<void(const std::string&)>& onLine);
    void kill();
private:
#ifdef _WIN32
    void* proc_ = nullptr; void* thread_ = nullptr; void* out_ = nullptr; void* job_ = nullptr;   // job: kills child trees (eclrun -> eclipse, mpiexec)
#else
    int pid_ = -1; int fd_ = -1;
#endif
    std::atomic<bool> killed_{false};
};

// Progress parsed from flow's terminal output.
struct OpmProgress {
    int reportStep = 0, reportSteps = 0;
    double day = 0, totalDays = 0;
    int timesteps = 0, newton = 0, linear = 0;     // end-of-run statistics
    bool finished = false;                          // "End of simulation" seen
    std::string lastError;                          // last "Error:" line
    bool parse(const std::string& line);            // true if progress changed
};

// run folder next to the project (<project>_OPM / <project>_ECL) or in %TEMP%\ResSimStudio\<tag>\<title>
std::string runDirFor(const std::string& projectPath, const std::string& title, const std::string& tag = "OPM");
std::string caseBaseName(const ressim::Case& c);
