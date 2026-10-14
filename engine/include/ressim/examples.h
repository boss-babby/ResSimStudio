// Built-in example library: ready-to-run cases from a 1D Buckley-Leverett check to faulted corner-point and
// channelised models. Every example is generated in code (deterministic), so it can also be exported as an
// ECLIPSE deck (see tools/export_examples or File > Examples > Export in the app).
#pragma once
#include "ressim/case.h"
#include <string>
#include <vector>

namespace ressim {

enum class ExampleLevel { Basic, Intermediate, Advanced };

struct ExampleInfo {
    std::string id;            // stable identifier, also the exported deck name (upper-cased)
    std::string title;
    ExampleLevel level = ExampleLevel::Basic;
    std::string summary;       // one line
    std::string details;       // what it demonstrates / what to look at
    std::string grid;          // e.g. "25 x 25 x 1 Cartesian"
    double runSeconds = 1;     // rough run time on a 4-core desktop
};

const char* toString(ExampleLevel l);
std::vector<ExampleInfo> listExamples();
// Builds the example. Throws std::runtime_error for an unknown id. The history-matching exercise runs its hidden
// "truth" model once to create the observed data (a few seconds).
Case makeExample(const std::string& id);
// SPE1 black-oil PVT + generic sandstone SWOF/SGOF (FIELD or METRIC): defaults for decks imported leniently.
void defaultBlackOil(Case& c, bool metric);

} // namespace ressim
