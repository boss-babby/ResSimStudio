// Compositional deck hook-up (milestone C3): ECLIPSE 300 decks (COMPS) -> CompCase + schedule, run and summary.
//
// Geometry, rock, wells and schedule come from the black-oil importer (corner point, faults, NNC, MULTREGT, PINCH, LGR,
// dual porosity: the processed Grid connection list is used as is, in SI). The fluid comes from eos::fluidFromDeck
// (EOS, CNAMES, TCRIT, ..., RTEMP, ZMFVD, ZI, STCOND). Saturation functions: SWFN / SGFN / SOF3 (or SWOF / SGOF),
// table 1. Separator: FIELDSEP (stages, last liquid to standard conditions). Injection gas: WELLSTRE + WINJGAS STREAM.
//
// Initialisation (E300 EQUIL with ZMFVD / ZI): hydrocarbon (oil-phase) pressure integrated from the datum with the EOS
// density of the local composition z(depth); water pressure integrated from the OWC with p_w(OWC) = p_o(OWC) - Pcow(OWC);
// Sw = Pcow^-1(p_o - p_w) at the cell centre (sharp contact without capillary pressure); gas caps follow from the
// composition (flash).
#pragma once
#include "ressim/case.h"
#include "ressim/compsim.h"
#include <functional>
#include <string>
#include <vector>

namespace ressim {
struct DeckRestart;   // binio.h
namespace deck { struct RawDeck; }
namespace comp {

struct DeckModel {
    CompCase base;                                  // geometry, rock, fluid, initial state (wells empty)
    std::vector<double> stepDays;                   // report step lengths (days)
    std::vector<std::vector<CompWell>> stepWells;   // wells and controls in force during each report step
    std::vector<std::vector<CompGroup>> stepGroups; // group targets (GCONPROD / GCONINJE RATE) of each report step
    std::vector<std::string> wellNames;
    std::vector<int> cellGlobal;                    // simulator cell -> global grid cell
    bool field = false;                             // deck unit system (summary output units)
    std::vector<std::string> notes;                 // simplifications / unsupported items
};

// Builds the model from a parsed deck and its imported case. Returns false (with *error) when the deck cannot be
// simulated by the compositional core.
bool deckModel(const deck::RawDeck& D, const Case& c, DeckModel& out, std::string* error = nullptr);
// Builds the model from a case with c.comp.enabled (imported deck or built in the UI): compositional EOS or CO2STORE
bool deckModel(const Case& c, DeckModel& out, std::string* error = nullptr);
bool equilibrateCo2(const Case& c, CompCase& cc, std::string* error = nullptr);
bool deckModelFromFile(const std::string& path, DeckModel& out, std::string* error = nullptr);
// hydrostatic initialisation of base (p0, z0, sw0) from EQUIL + ZMFVD / ZI; called by deckModel
bool equilibrate(const Case& c, const eos::CompFluid& fl, CompCase& cc, std::string* error = nullptr);

// summary of one report step (deck units: FIELD stb, Mscf, psia / METRIC sm3, bar)
struct CompSummary {
    double time = 0;                                // days
    double fopr = 0, fgpr = 0, fwpr = 0, fgir = 0, fwir = 0;
    double fopt = 0, fgpt = 0, fwpt = 0, fgit = 0, fwit = 0;
    double fpr = 0, foip = 0, fgip = 0, fwip = 0;
    // CO2STORE inventory (kg CO2): free gas not residually trapped (incl. the structurally held part co2Struct) /
    // residually trapped (Land / critical saturation) / dissolved in brine / mineralised; co2KrImmob = gas in cells with
    // krg = 0 (overlaps the others); saltSolid = precipitated NaCl (kg)
    double co2Mob = 0, co2Immob = 0, co2Diss = 0, co2Struct = 0, co2Min = 0, co2KrImmob = 0, saltSolid = 0;
    struct Well { double oil = 0, gas = 0, water = 0, ginj = 0, winj = 0, bhp = 0; bool open = false, onBhp = false;
                  std::vector<double> molar; };   // report-step average component rates (mol/s, production > 0, injection < 0)
    std::vector<Well> wells;
};
// Runs the schedule (maxSteps < 0: all). The callback gets every report step; returning false stops the run.
bool runDeckModel(const DeckModel& m, std::vector<CompSummary>& out, const std::function<bool(const CompSummary&, CompSimulator&)>& cb = {},
                  std::string* error = nullptr, int maxSteps = -1, CompStats* stats = nullptr, const DeckRestart* restart = nullptr);

}  // namespace comp
}  // namespace ressim
