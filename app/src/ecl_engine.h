// Schlumberger ECLIPSE as an external simulation engine (same workflow as OPM Flow: export the case as a deck, run it
// with live progress, read the ECLIPSE result files back). ResSim does not ship ECLIPSE: a licensed local installation
// is located automatically (eclrun.exe launcher preferred, eclipse.exe directly as a fallback) or chosen by the user.
// No ImGui / Direct3D dependency (also builds natively on Linux for testing).
#pragma once
#include <string>
#include <vector>

// Searches (in order): 'hint' (file or folder), %ECLIPSE_EXE% / %ECLRUN% / %ECLPATH%, PATH, <drive>:\ecl for drives
// C-F, Program Files\Schlumberger and \SLB (one launcher folder per suite). Empty when not found.
std::string detectEclipseExecutable(const std::string& hint = std::string(), std::vector<std::string>* searched = nullptr);
// true for the eclrun launcher (eclrun.exe / eclrun), false for a simulator executable (eclipse.exe, e300.exe)
bool isEclrunLauncher(const std::string& exe);
// Installed ECLIPSE versions, newest first: "eclrun --report-versions eclipse" for the launcher, the version folder
// (<root>\2023.1\bin\...) for a direct simulator executable. Empty if unknown.
std::vector<std::string> eclipseVersions(const std::string& exe, const std::string& program = "eclipse");

// Progress parsed from the ECLIPSE terminal output (E100 "STEP n TIME= t DAYS (... n ITS)" lines and @-messages).
struct EclProgress {
    double day = 0;
    int timesteps = 0, newton = 0;
    int errors = 0, problems = 0, bugs = 0;    // counted @--ERROR / @--PROBLEM / @--BUG messages
    int reportedErrors = -1;                   // "Errors  n" from the end-of-run summary (-1 = not seen)
    bool licenseProblem = false;
    std::string lastError;                     // text of the first error messages
    bool parse(const std::string& line);       // true if progress changed
private:
    int collecting_ = 0;                       // 1 = continuation lines of an error message follow
};

// The first error / bug messages from an ECLIPSE print file (CASE.PRT), for the failure message. Empty if none.
std::string eclipsePrtErrors(const std::string& prtPath, int maxMessages = 4);
