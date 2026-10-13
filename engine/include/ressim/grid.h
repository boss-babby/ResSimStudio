// Processed grid: active cells, pore volumes, depths and the connection (transmissibility) graph.
// The solver is connection-based, so corner-point grids, faults and NNCs plug in later
// by just adding connections.
#pragma once
#include "ressim/case.h"
#include <vector>
#include <array>

namespace ressim {

struct Vec3 { float x, y, z; };

struct Connection {
    int a, b;          // active cell indices
    double trans;      // transmissibility (cDarcy included)
    double dDepth;     // depth(b) - depth(a)
    char dir;          // 'X','Y','Z', 'N' explicit NNC, 'M' matrix-fracture (dual porosity)
};

class Grid {
public:
    void build(const Case& c);            // throws std::runtime_error on bad input

    int nx = 0, ny = 0, nz = 0;
    int numActive() const { return int(globalOf.size()); }
    int numCells() const { return nx * ny * nz; }

    std::vector<int> activeOf;            // global -> active (-1 inactive)
    std::vector<int> globalOf;            // active -> global
    std::vector<double> poreVolRef;       // per active cell, at reference pressure (reservoir volume units)
    std::vector<double> depth;            // cell-centre depth, per active cell
    std::vector<double> dzCell;           // thickness, per active cell
    std::vector<double> topDepth;         // per global cell (filled tops)
    std::vector<double> cellDx, cellDy, cellDz, bulkVol;   // per global cell (for well indices / statistics)
    std::vector<Vec3> center;             // per global cell centre
    std::vector<Connection> conns;

    // adjacency (CSR over active cells): for cell c, entries adjStart[c]..adjStart[c+1]
    std::vector<int> adjStart, adjConn, adjNbr;

    // geometry for visualisation: 8 corners per global cell (x right, y forward, z depth)
    // order: (i,j,k) top: 0:(0,0) 1:(1,0) 2:(0,1) 3:(1,1); bottom: 4..7 same pattern
    std::vector<std::array<Vec3, 8>> corners;
    Vec3 bbMin{}, bbMax{};

    // area of a cell face (dir 'X','Y','Z'; minus = low side), from the corner geometry
    double faceArea(int g, char dir, bool minus) const;
    bool cornerPoint = false;
    int numNNC = 0;
    // dual porosity: number of matrix layers (0 = single porosity), matrix-fracture connections ('M'), and the factor
    // applied to the input permeability of every global cell (fracture porosity unless NODPPM; empty = 1)
    int nzMatrix = 0, dualMode = 0, numMF = 0, numPinch = 0, numPinchedCells = 0;
    std::vector<double> permScale;
    bool isFracture(int g) const { return nzMatrix > 0 && g / (nx * ny) >= nzMatrix; }
    double effPerm(const std::vector<double>& k, int g) const { return k[size_t(g)] * (permScale.empty() ? 1.0 : permScale[size_t(g)]); }                       // connections that are not k-to-k neighbours (faults / explicit NNC)

    void ijk(int g, int& i, int& j, int& k) const { i = g % nx; j = (g / nx) % ny; k = g / (nx * ny); }
};

} // namespace ressim
