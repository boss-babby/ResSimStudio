// Project files (.rsim, human-readable text) and ECLIPSE deck export.
#pragma once
#include "ressim/case.h"
#include <string>
#include <vector>

namespace ressim {

// Throws std::runtime_error on I/O or parse errors.
void saveProject(const Case& c, const std::string& path);
Case loadProject(const std::string& path);
// complete case as text (project format 2: ECLIPSE deck + "-- RESSIM:" lines) and back (formats 1 and 2); used by the
// results files (.rsres) to store the case that was run
std::string projectText(const Case& c);
Case projectFromText(const std::string& text, const std::string& baseDir = "");
std::string projectToString(const Case& c);
Case projectFromString(const std::string& text);

// ---- ECLIPSE 100 decks (two-way).
// Export writes RUNSPEC/GRID/EDIT/PROPS/REGIONS/SOLUTION/SUMMARY/SCHEDULE; quantities that ECLIPSE cannot
// express (e.g. relative Newton pressure limit, thread count) go into "-- RESSIM:" comment lines,
// so export -> import -> export is lossless.
void exportEclipseDeck(const Case& c, const std::string& path);
std::string eclipseDeckToString(const Case& c);

// ---- ECLIPSE deck import.
// Two stages: deck_parser.h reads any deck generically (2900+ keywords from the OPM keyword database, INCLUDE /
// PATHS / SKIP, multi-record and table keywords), then the importer maps the supported keywords onto a Case:
// RUNSPEC dims/phases/units/START, GRID (Cartesian or corner-point arrays, BOX, EQUALS/MULTIPLY/ADD/COPY/
// MINVALUE/MAXVALUE, EQUALREG/MULTIREG/ADDREG/COPYREG, FAULTS/MULTFLT/MULTREGT/NNC, MINPV, MULTPV),
// PROPS (multi-region PVTW/ROCK/DENSITY/GRAVITY/PVDO/PVCDO/PVDG/PVTO, SWOF/SGOF and SWFN/SGFN/SOF3/SOF2/SLGOF),
// REGIONS (SATNUM/PVTNUM/EQLNUM/FIPNUM/MULTNUM/...), SOLUTION (multi-region EQUIL with RSVD/PBVD, or
// enumerated PRESSURE/SWAT/SGAS/RS, AQUFETP/AQUCT/AQUANCON) and SCHEDULE (WELSPECS/COMPDAT/WCONPROD/WCONHIST/
// WCONINJE/WCONINJH/WELTARG/WELOPEN/WPIMULT/TSTEP/DATES/TUNING). Everything else is skipped with a message.
// "-- RESSIM:" comment lines carry TUNING extras, AQVISC, PARAM and OBS (history data).
// Throws std::runtime_error ("KEYWORD (file:line): text") only when no runnable model can be built.
struct DeckMessage {
    int level = 1;                 // 0 info, 1 warning, 2 error
    std::string keyword, file;
    int line = 0;
    std::string text;
};
Case importEclipseDeck(const std::string& path, std::vector<std::string>* warnings = nullptr,
                       std::vector<DeckMessage>* report = nullptr);
// C31 ACTIONX: applies SCHEDULE keywords (deck text) to a case at run time - they act from report step 'step' on and win
// over the deck's own records of that step; keywords the engine cannot apply at run time are reported in 'messages'
bool applyScheduleText(Case& c, int step, const std::string& text, std::vector<std::string>* messages = nullptr);
Case eclipseDeckFromString(const std::string& text, const std::string& baseDir = "",
                           std::vector<std::string>* warnings = nullptr, std::vector<DeckMessage>* report = nullptr);
namespace deck { struct RawDeck; }
// E300 compositional (COMPS) / CO2STORE data of a parsed deck -> c.comp (called by the importer). False: not compositional.
bool compositionalFromDeck(const deck::RawDeck& D, Case& c, std::vector<DeckMessage>* messages = nullptr);
// Deck text of the compositional keywords of one section (0 RUNSPEC, 1 PROPS, 2 SOLUTION, 3 SCHEDULE start) written
// from c.comp (structured export; also used by the UI deck preview)
std::string compositionalKeywordsText(const Case& c, int section);

} // namespace ressim
