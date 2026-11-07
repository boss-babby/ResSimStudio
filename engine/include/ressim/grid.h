// Processed grid: active cells, pore volumes, depths and the connection (transmissibility) graph.
// The solver is connection-based, so corner-point grids, faults and NNCs plug in later
// by just adding connections.
#pragma once
#include "ressim/case.h"
#include <vector>
#include <array>
#include <string>

namespace ressim {

struct Vec3 { float x, y, z; };

struct Connection {
    int a, b;          // active cell indices
    double trans;      // transmissibility (cDarcy included)
    double dDepth;     // depth(b) - depth(a)
    // adjoint field derivatives (C11): dT / d ln k of cell a / b in the connection's direction (harmonic average:
    // dT/d ln t_a = T (1/t_a) / (1/t_a + 1/t_b)), and the global cell whose directional multiplier (MULTX/Y/Z) acts on T
    double dTa = 0, dTb = 0;
    // conduction geometry A / (d_a + d_b) in grid length units (T without CDARCY, permeability and NTG; multipliers kept),
    // 0 = unknown (LGR / explicit NNC: the thermal engine falls back to T / k_harmonic). C17: T / k_harmonic is only exact for
    // equal half-cell distances - it halved matrix-fracture conduction and mis-weighted layers of different DZ and k.
    double geo = 0;
    int mq = -1;       // (C45: mq and dir packed at the end - 56 instead of 64 bytes per connection)
    char dir;          // 'X','Y','Z', 'N' explicit NNC, 'M' matrix-fracture (dual porosity)
};
static_assert(sizeof(Connection) <= 56, "Connection packing");

class Grid {
public:
    int numNonnc = 0;                     // C38 fault-throw connections dropped by NONNC
    int numEditNnc = 0;                   // C34 connections changed by EDITNNC / EDITNNCR
    void build(const Case& c);            // throws std::runtime_error on bad input
    // MINC shells of a DUALPORO grid (called by build after the matrix-fracture pass; appends the inner shell cells)
    void buildMinc(const GridInput& g, double cDarcy, size_t n, size_t nl, int nzh, int nm, const std::vector<double>& kx,
                   const std::vector<double>& ky, const std::vector<double>& kz);

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
    int numMultRegT = 0;                  // connections modified by MULTREGT
    std::vector<std::pair<int, double>> numAquP0;   // C34 numerical aquifer cells (active index, initial pressure or < 0)
    int numAquConns = 0;                  // C34 reservoir connections of numerical aquifers
    std::vector<double> permScale;
    bool fieldDerivOk = true;             // Connection::dTa / dTb / mq valid (false with LGRs / dual porosity)
    bool isFracture(int g) const { return nzMatrix > 0 && g < nx * ny * nz && g / (nx * ny) >= nzMatrix; }
    double effPerm(const std::vector<double>& k, int g) const { return k[size_t(g)] * (permScale.empty() ? 1.0 : permScale[size_t(g)]); }                       // connections that are not k-to-k neighbours (faults / explicit NNC)

    // ---- local grid refinements (CARFIN): global cell indices >= nx*ny*nz are local cells, appended LGR by LGR in local
    // natural order. Host cells are inactive; per-cell arrays (corners, centre, cellDx, ...) cover the local cells too.
    struct LgrInfo { std::string name; int offset = 0, nx = 1, ny = 1, nz = 1, i1 = 0, j1 = 0, k1 = 0, i2 = 0, j2 = 0, k2 = 0; };
    std::vector<LgrInfo> lgr;
    std::vector<int> lgrHost;             // per local cell (g - nx*ny*nz): host global cell
    std::vector<int> hostLgr;             // per global cell: LGR index of a host cell, -1 otherwise (empty without LGRs)
    int numLgrCells() const { return int(lgrHost.size()); }
    int numAllCells() const { return numCells() + numLgrCells() + numMincCells() + numHfCells(); }
    bool isLocal(int g) const { return g >= numCells() && g < numCells() + numLgrCells(); }
    // ---- MINC shells (global index >= numCells() + numLgrCells()): host matrix cell and shell number (2 .. mincN)
    std::vector<int> mincHost, mincShell;
    int mincN = 0;
    std::vector<double> mincFrac;         // volume fraction of each shell (outer -> inner)
    int numMincCells() const { return int(mincHost.size()); }
    bool isMincShell(int g) const { return !mincHost.empty() && g >= numCells() + numLgrCells() && g < hfBase(); }
    // ---- C15 hydraulic fracture cells (EDFM; global index >= hfBase()), one per fracture piece in a host cell
    struct HfCell { int frac = 0, host = 0, global = 0; double kf = 0, w = 0, len = 0, h = 0, gamma = 0; };
    std::vector<HfCell> hfCells;
    int numHfCells() const { return int(hfCells.size()); }
    int hfBase() const { return numCells() + numLgrCells() + numMincCells(); }
    bool isHf(int g) const { return !hfCells.empty() && g >= hfBase(); }
    void buildHydraulicFractures(const Case& c, double cDarcy);
    bool isHost(int g) const { return !hostLgr.empty() && g < numCells() && hostLgr[size_t(g)] >= 0; }
    int localCell(int l, int i, int j, int k) const {
        const LgrInfo& L = lgr[size_t(l)];
        return numCells() + L.offset + i + L.nx * (j + L.ny * k);
    }
    // LGR index and local (i, j, k) of a local cell; false for a global cell
    bool localIJK(int g, int& l, int& i, int& j, int& k) const {
        if (!isLocal(g)) return false;
        const int o = g - numCells();
        for (l = int(lgr.size()) - 1; l > 0 && lgr[size_t(l)].offset > o; --l) {}
        const LgrInfo& L = lgr[size_t(l)];
        const int r = o - L.offset;
        i = r % L.nx; j = (r / L.nx) % L.ny; k = r / (L.nx * L.ny);
        return true;
    }
    // local cells of host cell 'host' on its face dir ('X','Y','Z'; minus = low side)
    std::vector<int> localCellsOnFace(int host, char dir, bool minus) const;
    // a global completion (COMPDAT) in a refined cell: the local cells of the host along the well direction through the
    // local column nearest the host centre (ECLIPSE expects COMPDATL there)
    std::vector<int> hostWellColumn(int host, char dir) const;
    // "(i, j, k)" (1-based) of a global cell, "NAME (i, j, k) in (I, J, K)" of a local cell
    std::string cellLabel(int g) const;
    // (i, j, k) of a global cell; a local cell reports its host cell
    void ijk(int g, int& i, int& j, int& k) const {
        if (isHf(g)) g = hfCells[size_t(g - hfBase())].host;
        if (isMincShell(g)) g = mincHost[size_t(g - numCells() - numLgrCells())];
        else if (g >= nx * ny * nz && !lgrHost.empty()) g = lgrHost[size_t(g - nx * ny * nz)];
        i = g % nx; j = (g / nx) % ny; k = g / (nx * ny);
    }
};

// Appends the local cells of every LGR to the per-cell arrays of the case (host values, or the LGR's own 'props'),
// so that index g >= nx*ny*nz addresses local cells everywhere. Transmissibility multipliers are not inherited
// (the host's face multipliers act on the LGR boundary through the redistributed host connections).
void expandLgrArrays(Case& c);
// DPNUM = 0 cells of a dual porosity model become single porosity: the cell keeps its matrix rock, regions and initial
// state but is represented by its fracture-layer cell (the matrix-layer cell is made inactive), so it connects to the
// neighbouring fracture system (ECLIPSE: single porosity matrix cells connect to the fracture cells of dual neighbours).
// Idempotent; returns the number of single porosity cells.
int applyDpnum(Case& c);
// MINC: appends the inner shells of every matrix cell (host values) to the per-cell arrays (see GridInput::minc).
// Idempotent; returns the number of shell cells added.
int expandMincArrays(Case& c);
// C15: hydraulic fractures (GridInput::hfracs): computes the fracture pieces from the host geometry and appends one entry
// per piece (the host's value) to every per-cell array; returns the number of fracture cells (idempotent)
int prepareHydraulicFractures(Case& c);

} // namespace ressim
