// Input generators used by the UI: PVT correlations, Corey relative permeability, heterogeneous grids.
#pragma once
#include "ressim/case.h"
#include <cstdint>
#include <string>

namespace ressim {

// ---- black-oil PVT from correlations (Standing Rs/Bo, Beggs-Robinson + Vasquez-Beggs viscosity,
//      Vasquez-Beggs undersaturated compressibility, Sutton pseudo-criticals + Hall-Yarborough z,
//      Lee-Gonzalez-Eakin gas viscosity). Produces PVTO + PVDG + surface densities in the case units.
struct CorrelationInput {
    double api = 35, gasGravity = 0.75, temperatureF = 200;
    double pbMax = 4500;            // highest bubble point in the table (psia)
    double pMax = 9000;             // highest pressure (psia)
    int nRows = 12;
};
void makeCorrelationPVT(const CorrelationInput& in, UnitSystem units, PropsInput& out, bool deadOil = false);
double hallYarboroughZ(double ppr, double tpr);

// ---- Corey curves
struct CoreyInput {
    double swc = 0.2, sorw = 0.2, sgc = 0.05, sorg = 0.2;
    double nw = 2, now = 2, ng = 2, nog = 2;
    double krwMax = 0.6, kroMax = 1.0, krgMax = 0.8;
    double pcEntry = 0, pcExponent = 2;   // Pcow = pcEntry * ((1-Se))^n  (pressure units)
    int rows = 16;
};
void makeCoreySWOF(const CoreyInput& c, std::vector<SatRow>& out);
void makeCoreySGOF(const CoreyInput& c, std::vector<SatRow>& out);

// ---- heterogeneous rock: lognormal permeability with Gaussian correlation (moving-average filter)
struct RockGenInput {
    std::vector<double> layerPerm, layerPoro, layerKvKh, layerNtg;   // per layer (size nz, or size 1 = all)
    double sigmaLn = 0.0;            // std dev of ln(k); 0 = homogeneous layers
    double corrX = 3, corrY = 3;     // correlation lengths in cells
    double poroPermExp = 0.12;       // poro = layerPoro * (k/kLayer)^exp  (correlated porosity)
    uint32_t seed = 12345;
};
void generateRock(const RockGenInput& in, GridInput& g);

// Uniform Cartesian geometry with optional dip; tops computed for top layer.
struct GeomInput {
    int nx = 20, ny = 20, nz = 4;
    double dx = 300, dy = 300;
    std::vector<double> dz{25, 25, 25, 25};   // per layer
    double top = 8325, dipX = 0, dipY = 0;     // dip in degrees along +x / +y
};
void generateGeometry(const GeomInput& in, GridInput& g);

// Convert DX/DY/DZ/TOPS to COORD/ZCORN with vertical pillars (cells keep their own flat tops/bottoms;
// stair-stepped dip becomes partial overlaps). DX/DY/DZ/TOPS are cleared afterwards.
void convertToCornerPoint(GridInput& g);
// Vertical fault on the pillar line between columns pos-1 and pos (along I when alongJ = false,
// i.e. an X-face fault; along J otherwise). Cells on the high side are moved down by throwDepth.
// Also appends a FAULTS entry with the given transmissibility multiplier (MULTFLT).
void addFaultThrow(GridInput& g, int pos, double throwDepth, double mult, const std::string& name, bool alongJ = false);

// Dual porosity conversion (ECLIPSE layout). makeDualPorosity doubles NZ: every per-cell array is copied to the new
// fracture layers NZ+1..2NZ, which get porosity fracPoro and permeability fracPerm (intrinsic: NODPPM is set), and
// SIGMA = 4 (2/l^2 + 1/h^2) from the matrix block size l (lateral) and h (height). With fractureRegion, the fracture
// cells get a new SATNUM region with straight-line relative permeability and zero capillary pressure.
// makeSinglePorosity drops the fracture layers again (keeps the matrix).
void makeDualPorosity(Case& c, int mode, double fracPoro, double fracPerm, double blockL, double blockH, bool fractureRegion);
void makeSinglePorosity(Case& c);

} // namespace ressim
