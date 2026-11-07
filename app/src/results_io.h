// C14: RunResults <-> ResSim results / restart file (.rsres) and the ECLIPSE result export of a run
#pragma once
#include "sim_runner.h"
#include "ressim/results_file.h"
#include <string>

// whole run (incl. the C13 geothermal / EDFM analysis and the engine restart states) -> results file
ressim::ResultsFile resultsToFile(const RunResults& r);
// results file -> RunResults (grid rebuilt from the stored case; no live simulator: streamline tracing / restore of
// report steps need a re-run). Throws std::runtime_error.
RunResults resultsFromFile(const ressim::ResultsFile& f);
// convenience
void saveRunResults(const RunResults& r, const std::string& path);
RunResults loadRunResults(const std::string& path);
// ECLIPSE binary result case (EGRID / INIT / SMSPEC / UNSMRY / UNRST) of a ResSim run; returns "" or the error
std::string exportRunEclipse(const RunResults& r, const std::string& base, std::vector<std::string>* files = nullptr,
                             std::vector<std::string>* warnings = nullptr);
