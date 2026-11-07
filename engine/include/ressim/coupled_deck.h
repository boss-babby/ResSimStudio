// C28: headless coupled runs of the compositional / CO2STORE and thermal deck engines (CLI, scripts, tests) - the same
// sequence as the GUI runners: after every report step FE geomechanics (compositional GEOMECH: pore volume + permeability
// two-way) or the THM mechanics stress of the thermal engine -> fault mechanics (FAULTMECH, FAULTSLIP quasi-static slip
// with DDM stress transfer, thermo-elastic stress 'STRESSFIELD ... TH') -> two-way fault transmissibility; post-run RSF.
// Results: the engine summaries plus the scalar KPIs of coupledbo::Output (same names as the black-oil runner).
#pragma once
#include "ressim/case.h"
#include "ressim/coupled_bo.h"
#include "ressim/comp_deck.h"
#include "ressim/therm_deck.h"

namespace ressim {
namespace coupleddeck {

// true when a compositional / thermal case asks for coupled mechanics handled here (fault mechanics, compositional GEOMECH)
bool needed(const Case& c);

bool runComp(const Case& c, const comp::DeckModel& m, std::vector<comp::CompSummary>& out, coupledbo::Output& co,
             const std::function<void(const std::string&)>& log = {}, const std::function<bool(const comp::CompSummary&)>& progress = {},
             std::string* error = nullptr, int maxSteps = -1, comp::CompStats* stats = nullptr,
             const std::function<void(const comp::CompSimulator&, const comp::CompSummary&)>& onState = {});   // C48: per report step (CLI snapshots)
bool runThermal(const Case& c, const thermal::DeckModel& m, std::vector<thermal::Summary>& out, coupledbo::Output& co,
                const std::function<void(const std::string&)>& log = {}, const std::function<bool(const thermal::Summary&)>& progress = {},
                std::string* error = nullptr, int maxSteps = -1, thermal::Stats* stats = nullptr,
                const std::function<void(const thermal::Simulator&, const thermal::Summary&)>& onState = {});

}  // namespace coupleddeck
}  // namespace ressim
