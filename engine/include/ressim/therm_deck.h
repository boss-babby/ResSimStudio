// Thermal deck hook-up (milestone C7): ECLIPSE 300 THERMAL decks <-> Case::thermal <-> thermal::Case + schedule.
//
// Import (thermalFromDeck, called by the deck importer when THERMAL is present): components (CNAMES, MW, TCRIT, PCRIT,
// CVTYPE, KVCR or Wilson K-values from PCRIT / TCRIT / ACF with KVWI), liquid density (DREF, CREF, THERMEX1), enthalpies
// (SPECHA / SPECHB / SPECHG / SPECHH, HEATVAPS, HEATVAP, HEATVAPE), gas Z-factor (ZFACTOR), viscosities (OILVISCT,
// GASVISCT or GASVISCF, WATVISCT), STCOND, rock (HEATCR, HEATCRT, THCONR, THCONSF), initial state (TEMPI or TEMPVD,
// ZI or XMF / YMF with EQUIL or PRESSURE / SWAT / SGAS) and steam injection (WINJTEMP). Geometry, rock compressibility,
// saturation functions, wells and the schedule come from the black-oil importer.
//
// Export (thermalKeywordsText): the same keywords written from Case::thermal. ECLIPSE runs the exported deck with
// THSTT97 (IAPWS-IF97 water / steam, as the ResSim core) and a WATVISCT table from IAPWS-2008 when none was given;
// heat loss to over / underburden is written as a '-- RESSIM: HEATLOSS' line (ECLIPSE ROCKCON / ROCKPROP: see ROADMAP).
#pragma once
#include "ressim/case.h"
#include "ressim/thermsim.h"
#include <functional>
#include <set>
#include <string>
#include <vector>

namespace ressim {
struct DeckRestart;   // binio.h
namespace deck { struct RawDeck; struct Record; }
struct DeckMessage;

// E300 THERMAL keywords read into Case::thermal (kept verbatim unless thermal.structured)
const std::set<std::string>& thermalKeywords();
// Parsed deck -> c.thermal (false: no THERMAL keyword / no usable fluid). Called by the importer.
bool thermalFromDeck(const deck::RawDeck& D, Case& c, std::vector<DeckMessage>* messages = nullptr);
// Deck text of the thermal keywords of one section: 0 RUNSPEC, 1 PROPS, 2 SOLUTION, 3 SCHEDULE (WINJTEMP of report step
// 'step'), 4 GRID (rock heat arrays; written in both export modes)
std::string thermalKeywordsText(const Case& c, int section, int step = 0);
// C25 shared '-- RESSIM: GEOMECH / MECHBC / GEOPLAST' records (thermal: Case::thermal.geo, black oil / compositional:
// Case::geoSeq): parse one record (false = not a mechanics key) / write the lines of an enabled input
bool parseGeoRecord(const std::string& key, const deck::Record& r, geomech::Input& G);
std::string geoKeywordsText(const geomech::Input& G);

namespace thermal {
struct DeckModel {
    Case base;                                      // geometry, rock, fluid, initial state (wells of step 0)
    std::vector<double> stepDays;                   // report step lengths (days)
    std::vector<std::vector<Well>> stepWells;       // all wells (open flag) and their controls per report step
    std::vector<std::string> wellNames;
    std::vector<int> cellGlobal;                    // simulator cell -> global grid cell
    bool field = false;
    std::vector<std::string> notes;
};
// Case with c.thermal.enabled (imported deck or built in the UI) -> thermal core model
bool deckModel(const ressim::Case& c, DeckModel& out, std::string* error = nullptr);
bool deckModelFromFile(const std::string& path, DeckModel& out, std::string* error = nullptr);

// summary of one report step (deck units: METRIC sm3, bar, C, kJ | FIELD stb, Mscf, psia, F, Btu)
struct Summary {
    double time = 0;                                // days
    double fopr = 0, fwpr = 0, fgpr = 0, fwir = 0;  // fwir: cold-water-equivalent steam / water injection
    double fopt = 0, fwpt = 0, fgpt = 0, fwit = 0;
    double fpr = 0, favgT = 0, fmaxT = 0, foip = 0; // pore-volume averages, max temperature
    double fhlt = 0, fenergyIn = 0, csor = 0;       // cumulative heat loss, net injected enthalpy, cum. steam-oil ratio
    double fsteamPV = 0;                            // steam (gas) pore volume fraction
    double fgir = 0, fgit = 0;                      // gas injection (air / O2 / CO2 / N2 streams), surface gas volume
    struct Well { double oil = 0, water = 0, gas = 0, winj = 0, bhp = 0, hrate = 0, tfluid = 0; bool open = false, onBhp = false;
                  std::vector<double> molar; };   // report-step average component rates (mol/s, nc + 1, water last; production > 0)
    std::vector<Well> wells;
};
bool runDeckModel(const DeckModel& m, std::vector<Summary>& out, const std::function<bool(const Summary&, Simulator&)>& cb = {},
                  std::string* error = nullptr, int maxSteps = -1, Stats* stats = nullptr, const DeckRestart* restart = nullptr);
}  // namespace thermal
}  // namespace ressim
