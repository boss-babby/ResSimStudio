// C48: per-report-step cell maps (Snapshot) of the compositional / CO2STORE and thermal deck engines, shared by the GUI
// runner (app/src/sim_runner.cpp) and the headless CLI (ResSimCLI -o out.rsres), so both write the same 3-D arrays.
// Engine-only content: coupled FE-geomechanics maps of the compositional runner and fault-mechanics updates stay with
// the caller.
#pragma once
#include "ressim/simulator.h"
#include "ressim/compsim.h"
#include "ressim/thermsim.h"
#include "ressim/case.h"

namespace ressim {
namespace comp {
// PRESSURE (case units: pu = Pa per unit), SWAT / SGAS / SOIL, XCO2 slot (CO2STORE: dissolved CO2 of the brine; else z1
// or the overall CO2 fraction with CO2SOL), SGTRAP [SGMAX], SSALT, MINCO2, PORFAC / PERMFAC, ZMFk / XMFk / YMFk (nc <= 16),
// AMFk / XCO2AQ (CO2SOL)
Snapshot cellSnapshot(const CompSimulator& sim, double time, double pu);
}  // namespace comp
namespace thermal {
// PRESSURE, saturations, TEMP / DTEMP (output units: C or F; t0 = the step-0 temperature map, filled when empty), SSTEAM,
// ZMFk, XMFk / YMFk / AMFk, SGTRAP [SGMAX], XCO2 slot (tracked component), POROF / PERMMULT, SOLID_<name>, FRACAP / FRACP /
// FRACT (EDFM), THM mechanics (UPLIFT, DISPX, STRAINV, SIGEFF*, FAILURE, POROM, SIGMIN [EPSPLAS, EPSPV, KDAMAGE]).
// n = number of active cells of the grid (simulator cells are the active cells in order).
Snapshot cellSnapshot(const Simulator& sim, const ThermalInput& TI, double time, size_t n, bool field, std::vector<float>& t0);
}  // namespace thermal
}  // namespace ressim
