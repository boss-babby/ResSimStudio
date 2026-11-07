// C25: finite-element geomechanics and fracture growth for the black-oil and compositional engines (sequential coupling)
// on Cartesian AND corner-point grids - the non-XFEM path of the geological / geomechanical toolchain.
//
// Coupler   the mechanics mesh is the box of the flow grid (nx ny nz, optional over / under / side-burden padding of
//           geomech::Input). Cartesian tensor grids use the brick elements of geomech::Solver::setup; corner-point grids
//           (COORD / ZCORN, dipping TOPS, irregular DX) use trilinear isoparametric hexahedra (Solver::setupHex) whose
//           nodes are the averages of the corner points of the cells sharing a node (fault throws are smeared over one
//           element: the mechanics mesh is conforming, the flow grid keeps its NNCs).
//           After every report step: dp = p - p0 per box cell -> FE solve (elastic or Drucker-Prager) -> full stress
//           tensor, strain, displacement, failure index per cell.
//           Feedback (GEOMECH coupling 1, explicit / one report step lag):
//             pore volume   PV_flow x m,  m = (phi_FE / phi0) exp(-c_r (p - p0))  (the FE porosity replaces the flow rock
//                           compressibility; the flow engines keep the moles / masses of the old time level, so a PV change
//                           is mass conservative and acts on the pressure of the next step)
//             permeability  the cell multiplier of geomech (damage exp(b eps_p) x stress law exp(-permStress dsigma'_m)),
//                           transmissibility x harmonic mean of the two cells (matrix-matrix connections only)
//           Products: sigma_min per cell for HFGROW ... FESTRESS and the total stress tensor for STRESSFIELD ... FE.
// HfSequential  LEFM growth of the HFRAC pieces between report steps for the compositional engine (the black-oil engine
//           grows them inside the time step): same rule as Simulator::growFractures (open the next column when
//           K_I = (p_tip - sigma) sqrt(pi c) >= K_IC, c = min(L_wing, H / 2)); closed pieces conduct like host rock
//           (fracture-fracture transmissibility x min(1, k_n / k_f) x closedMult).
#pragma once
#include "ressim/geomech.h"
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace ressim {
class Grid;
struct Case;
class Simulator;
struct Snapshot;
namespace geoseq {

// nodes (x, y, depth in m) of the box grid nx ny nz from the cell corners: ((nx+1)(ny+1)(nz+1)) x 3, empty without corners
std::vector<double> gridNodes(const Grid& g, double lenToM);
// true when the grid is a tensor-product Cartesian box (uniform top, DX per column, DY per row, DZ per layer)
bool tensorBox(const Grid& g, double tol = 1e-6);

class Coupler {
public:
    // cellGlobal: flow cell -> global cell of g (cells outside the box - LGR / MINC / HF - have no mechanics element);
    // pPa / phi per flow cell (initial state). generalHex forces the isoparametric mesh even for tensor boxes.
    bool setup(const geomech::Input& in, const Grid& g, const std::vector<int>& cellGlobal, double lenToM,
               const std::vector<double>& pPa, const std::vector<double>& phi, std::string& err, bool generalHex = false);
    bool ready() const { return ready_; }
    bool twoWay() const { return in_.coupling != 0; }
    bool hexMesh() const { return S_.generalHex(); }
    // one FE solve at the current flow state (Pa per flow cell; temperature change K per flow cell optional) + commit
    bool step(const std::vector<double>& pPa, const std::vector<double>* dTK = nullptr);
    // feedback factors per flow cell (1 for unmapped cells)
    std::vector<double> pvFactor(const std::vector<double>& pPa, double crPerPa) const;
    std::vector<double> permFactor() const;
    // transmissibility multipliers of the given flow connections (pairs of flow cells): harmonic cell multiplier when both
    // cells have mechanics elements, else 1
    std::vector<double> connFactor(const std::vector<std::pair<int, int>>& conns) const;
    std::vector<double> minStress() const;                 // Pa per flow cell (0 = unmapped)
    bool totalStress(int cell, double S[6]) const;         // flow cell, Pa, compression +
    int boxOf(int cell) const { return cell >= 0 && size_t(cell) < box_.size() ? box_[size_t(cell)] : -1; }
    const geomech::Solver& solver() const { return S_; }
    std::vector<float> map(const std::vector<double>& boxValues, double scale) const;   // per flow cell
    int steps() const { return steps_; }
    double maxFailure() const;
private:
    bool ready_ = false;
    geomech::Input in_;
    geomech::Solver S_;
    std::vector<int> box_;               // flow cell -> box cell (-1)
    std::vector<double> p0_;             // flow cell initial pressure (Pa)
    std::vector<double> phi0b_;          // box porosity at setup
    int nBox_ = 0, steps_ = 0;
};

// C26 (C25-4a, GEOMECH coupling 2): one black-oil report step with ITERATIVE flow <-> FE coupling. The step is run, the FE
// state solved, and while the pore-volume multipliers m(p) differ from the ones the flow step used by more than tol (max
// absolute change), the step is repeated from its saved start state (Simulator::saveState / loadState; the old-time masses
// keep the start multipliers, the new time level uses the latest m) - a fixed point of flow and mechanics per report step
// instead of the one-report-step lag of the explicit scheme. pvApplied: multipliers in the flow model (updated). onRewind
// is called before every repetition (callers drop summary points they collected for the discarded pass).
struct IterStep { bool ran = false, feOk = true, converged = false; int passes = 0; double diff = 0; };
IterStep iterativeReportStep(Simulator& sim, Coupler& gs, std::vector<double>& pvApplied, double paPerUnit, double crPerPa,
                             int maxIt, double tol, const std::function<void()>& onRewind = {});

class HfSequential {
public:
    struct Event { double t = 0; int frac = 0, side = 0; double length = 0, KI = 0, pnet = 0; };
    // cx: case with the hydraulic fractures prepared (prepareHydraulicFractures), gx: Grid built from cx (the flow cells of
    // the compositional model are the active cells of gx in order)
    bool setup(const Case& cx, const Grid& gx, std::string& err);
    bool active() const { return active_; }
    // growth test at the current state: pressure per flow cell (case units), sigma_min per flow cell (case units; empty =
    // HFGROW gradient law); returns the number of opened piece columns
    int grow(double tDays, const std::vector<double>& pCase, const std::vector<double>& sMinCase);
    // multipliers of the fracture-fracture connections touching closed pieces (cell pairs of the flow cells)
    void connMultipliers(std::map<std::pair<int, int>, double>& m) const;
    const std::vector<Event>& events() const { return ev_; }
    int openPieces() const { int n = 0; for (char o : open_) n += o; return n; }
    int pieces() const { return int(open_.size()); }
    double wingLength(int frac, int side) const;
private:
    bool active_ = false;
    int nx_ = 0, ny_ = 0, nz_ = 0;
    double toPa_ = 1e5, toM_ = 1;
    std::vector<int> frac_, act_, hostAct_;          // per piece: fracture, flow cell, host flow cell
    std::vector<double> s_, half_, h_, depth_, kf_, kn_;
    std::vector<char> open_;
    struct FF { int a, b, ka, kb; };                   // fracture-fracture connection (flow cells, pieces)
    std::vector<FF> ff_;
    struct Frac { double kic = 0, sigma = 0, sigmaGrad = 0, zTop = -1, zBot = -1, closedMult = 1; int fe = 0; std::string name; };
    std::vector<Frac> F_;
    std::vector<Event> ev_;
};

// C30: finite-element geomechanics maps of a report snapshot (Snapshot::extra, per snapshot cell q, box(q) = mechanics
// cell or -1): UPLIFT DISPX DISPY DISPMAG (m, uplift + up), STRAINV (-), SIGEFFM SIGEFFV SIGEFFH SIGMIN SQRTJ2 (pressure
// units, compression +), FAILURE (Mohr-Coulomb index), POROM, KDAMAGE, with plasticity EPSPLAS EPSPV. Existing arrays of the
// same name are replaced. pu = Pa per case pressure unit.
void mechMaps(Snapshot& sn, const geomech::Solver& g, const std::function<int(size_t)>& box, size_t n, double pu, bool plastic);

// C30: one-way FE geomechanics as a POST-PROCESSING step of finished results - any engine (black oil, compositional,
// thermal) and imported ECLIPSE / OPM restart results: the snapshots' pressure (and TEMP, when present: thermal stress)
// changes from the first snapshot load the mechanics mesh of the grid (Cartesian bricks or corner-point hexahedra, padding
// and boundary conditions of 'in'), one FE solve per snapshot in time order (path-dependent with plasticity), and the
// mechMaps arrays are added to every snapshot. snaps: per active cell of g. progress(i, n) returning false cancels.
struct PostResult { int steps = 0; double seconds = 0, maxUplift = 0, maxSubsidence = 0, maxFailure = 0; std::string error; };   // uplift / subsidence: reservoir cells (m)
PostResult postProcess(const geomech::Input& in, const Grid& g, bool fieldUnits, std::vector<Snapshot>& snaps,
                       const std::function<bool(int, int)>& progress = {});

}  // namespace geoseq
}  // namespace ressim
