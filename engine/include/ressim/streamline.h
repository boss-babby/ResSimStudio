#pragma once
// Streamline tracing (Pollock's semi-analytical method) and time-of-flight.
// Tracing is done in each cell's reference unit cube using the total reservoir face fluxes, so it works
// unchanged on corner-point grids with faults and non-neighbour connections.
#include <string>
#include <vector>
#include "ressim/grid.h"

namespace ressim {
class Simulator;

struct Streamline {
    std::vector<Vec3> pts;     // physical points (x, y, depth)
    std::vector<float> tof;    // time-of-flight at each point [days]
    std::vector<int> cells;    // active cells visited
    int source = -1;           // index into StreamlineResult::sources
    int sink = -1;             // index into StreamlineResult::sinks, -1 = did not reach a producer
    double flux = 0;           // reservoir rate carried [rb/d or rm3/d]
};

struct StreamlineResult {
    std::vector<Streamline> lines;
    std::vector<float> cellTof;     // per active cell: flux-weighted mean TOF from injectors (-1 = not swept)
    std::vector<int> cellSink;      // per active cell: producer that drains it (-1 = none)
    std::vector<int> cellSource;    // per active cell: injector that sweeps it (-1 = none)
    std::vector<std::string> sources, sinks;    // injector / aquifer names, producer names
    std::vector<std::vector<double>> pairFlux;  // [source][sink] reservoir rate (well allocation factors)
    double maxTof = 0;
    int stalled = 0;                            // lines that stopped before a producer
};

struct StreamlineOptions {
    int lines = 400;          // target number of streamlines
    int pointsPerCell = 3;    // drawing resolution inside each cell
    int maxCells = 0;         // per line (0 = 3 x active cells)
    // optional total reservoir-volume flux per grid connection (+ from a to b), e.g. the inter-block flows written by
    // ECLIPSE / OPM / FrontSim (RPTRST FLOWS: FLROILI+, FLRWATI+, FLRGASI+ ...). NaN entries use the model flux.
    const std::vector<double>* connFlux = nullptr;
};

// Uses the simulator's current state (call Simulator::restoreReport() first to trace an earlier report step).
StreamlineResult traceStreamlines(const Simulator& sim, const StreamlineOptions& opt = {});

} // namespace ressim
