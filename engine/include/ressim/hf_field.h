// C30: one-way field output of the in-run hydraulic fractures for the 3-D viewers (display only - nothing is fed back to
// the flow model; the conductance coupling stays with DDMWIDTH / PHASEFIELD TWOWAY / G1-c, which own setHfConductance).
//
// After every report step:
//   1. net pressure of every piece p_net = p_piece - sigma_closure (FE minimum stress with FESTRESS, else sigmaMin +
//      sigmaGrad (z - z_mid)), opening of the OPEN pieces from the dense 3-D displacement-discontinuity system
//      K_o w = -p_net with contact (w >= 0) - the same solve as DDMWIDTH (ddm.h, Nikkhoo & Walter 2015 triangular
//      dislocations; half space when the fractures are shallower than twice their extent, else full space)
//   2. stress, displacement and volumetric strain change of the matrix cells around the fractures (superposition of the
//      opening of every open piece over its two triangles), cells within radius x the largest fracture extent, nearest
//      first, capped by the evaluation budget (cells x open triangles); cells whose centre lies in a fracture plane are
//      evaluated at +-1/4 of the cell width along the plane normal and averaged (the TDE field is singular in the plane)
// Arrays (Snapshot::extra, per active cell, 0 elsewhere):
//   HFWIDTH  mm      opening of the piece (fracture cells) / largest opening of the pieces in the cell (host cells)
//   HFPNET   p unit  net pressure (fracture cells; host cells: largest of the open pieces)
//   HFCOND   -       conductance multiplier in use (fracture cells)
//   HFDSXX HFDSYY HFDSZZ HFDSXY  p unit  stress change, compression + (x east, y north, z vertical)
//   HFDSNN   p unit  normal-stress change on planes parallel to the nearest fracture ("stress shadow")
//   HFDSM    p unit  mean stress change;  HFDTAU p unit  largest horizontal shear-stress change
//   HFDISPX HFDISPY HFDISPZ HFDISP  mm  displacement (z + up) and its magnitude;  HFSTRV  microstrain  volumetric (+ dilation)
// Deck: '-- RESSIM: HFFIELD ON|OFF [E_GPa nu] [radius] [budget]' (Case::hfField; default: on for HFGROW / DDMWIDTH fractures).
#pragma once
#include "ressim/ddm.h"
#include <string>
#include <vector>

namespace ressim {
class Simulator;
struct Snapshot;
namespace hffield {

struct Stats {
    int pieces = 0, cells = 0, cellsInRadius = 0, open = 0, updates = 0, evaluated = 0;
    bool halfSpace = false;
    double E = 0, nu = 0, radiusM = 0, setupSeconds = 0, lastSeconds = 0, maxWidth = 0, maxShadow = 0;
    std::string note;
};

class Field {
public:
    // false (with a note) when there is nothing to do: no eligible fracture, too many pieces, HFFIELD OFF
    bool setup(const Simulator& sim, std::string& note);
    void addArrays(const Simulator& sim, Snapshot& s);
    const Stats& stats() const { return st_; }
    // stress change (Pa, tension +, xx yy zz xy xz yz) and displacement (m) at a point (m, z up) for the openings of the
    // last update (tests / probes)
    void fieldAt(const double x[3], double S[6], double u[3]) const;
    const std::vector<double>& width() const { return wAll_; }   // per fracture cell (m), last update
    const std::vector<double>& pnet() const { return pnAll_; }   // per fracture cell (case pressure units), last update
private:
    struct Obs { int act = 0; ddm::V3 x{}, x2{}; bool twin = false; ddm::V3 n{}; };
    std::vector<int> pieces_;              // fracture cell of every element
    std::vector<ddm::Element> el_;         // metres, z up
    std::vector<std::array<ddm::V3, 2>> slipDir_;   // per element and triangle: unit opening in the triangle frame (s, d, n)
    std::vector<double> Ko_;               // dense opening matrix (Pa / m)
    std::vector<double> zMid_;             // piece mid depth (case units)
    std::vector<Obs> obs_;
    std::vector<double> w_, wPrev_, wAll_, pnAll_;
    std::vector<std::vector<float>> cache_;
    bool active_ = false, have_ = false;
    double G_ = 0, nu_ = 0.25, Mlen_ = 1, Pa_ = 1e5, sense_ = 1;
    Stats st_;
    void evalPoint(const ddm::V3& x, double S[6], double u[3]) const;
};

}  // namespace hffield
}  // namespace ressim
