// C14: ResSim results + restart file (.rsres) - one self-contained binary file per run:
//   CASE      the project text of the case that was run (projectToString), so the results reopen without the project
//   META      engine, unit system, wells, components, statistics, initial in-place
//   SUMMARY   every accepted time step (SummaryPoint incl. per-well rates)
//   SNAP      one chunk per report step (all per-cell arrays, float32)
//   BLOB      named extra data of the producer (e.g. the app's geothermal / EDFM analysis)
//   RESTART   engine restart states (Simulator::saveState) at chosen report steps -> exact continuation
// Layout: 8-byte magic "RSRES\x1a\0\1", uint32 version, then chunks [8-char tag][uint64 size][payload]; unknown chunks
// are skipped (forward compatible). All numbers little-endian (x86-64 / ARM64 hosts).
#pragma once
#include "ressim/simulator.h"
#include "ressim/components.h"
#include <string>
#include <vector>

namespace ressim {

struct ResultsFile {
    std::string engine = "ResSim";              // "ResSim", "ResSim compositional", "ResSim thermal", "OPM Flow", ...
    std::string caseText;                       // projectToString(case)
    std::string title, created;                 // free text, ISO date
    std::vector<SummaryPoint> summary;
    std::vector<Snapshot> snaps;
    std::vector<std::string> wellNames;
    std::vector<int> wellTypes;                 // WellType as int
    std::vector<ComponentInfo> components;
    std::string trackName;
    InPlace initial;
    RunStats stats;
    double totalTime = 0;
    std::vector<std::pair<std::string, std::string>> blobs;   // (key, binary payload)
    struct Restart {
        int report = 0;                         // report step index (number of completed report steps)
        double time = 0;                        // days (black oil) / case time unit
        int snap = -1;                          // index into snaps of the same time (-1 none)
        size_t summaryCount = 0;                // summary points up to and including this restart
        std::string kind;                       // "BO", "COMP", "THERMAL"
        std::string state;                      // Simulator::saveState()
    };
    std::vector<Restart> restarts;
    const std::string* blob(const std::string& key) const {
        for (const auto& b : blobs) if (b.first == key) return &b.second;
        return nullptr;
    }
    const Restart* restartAt(int report) const {
        for (const auto& r : restarts) if (r.report == report) return &r;
        return nullptr;
    }
};

// Throw std::runtime_error on I/O or format errors.
void saveResultsFile(const ResultsFile& r, const std::string& path);
ResultsFile loadResultsFile(const std::string& path, bool loadSnapshots = true);
// header-only peek (case text, meta, list of restart steps without their states)
ResultsFile peekResultsFile(const std::string& path);

} // namespace ressim
