// C24: fault mechanics and induced seismicity INSIDE the flow run (the faults of the model are the scene objects).
//
// The faults are the FAULTS of the deck (the same faces that carry MULTFLT). For every fault face of an active cell:
//   geometry  centre, unit normal and area from the corner-point corners of the cell face; the connections that cross it
//             (k-neighbours and the NNCs of displaced columns) from the grid adjacency
//   stress    total in-situ tensor from '-- RESSIM: STRESSFIELD': Sv = sv0 + svGrad z, SHmax = rH Sv along azimuth theta,
//             Shmin = rh Sv, plus the poroelastic horizontal stress path  dS_H = dS_h = A (p - p_init)
//             sigma_n = n.S.n, tau = |S n - sigma_n n|, sigma'_n = sigma_n - p (p = mean of the two face cells)
//   Coulomb   CFF = tau - mu sigma'_n - c  (>= 0: the face is reactivated), dCFS = CFF - CFF(t0), slip tendency tau / sigma'_n
//   seismicity Dieterich (1994) rate-and-state seismicity rate per face with the background stressing rate Sr = A sigma / t_a:
//             d gamma / dt = (1 - gamma dS/dt) / (A sigma),  R = r / (gamma Sr),  r = r0 x face area
//             (exact gamma update for a constant stressing rate in 16 sub-steps per report step, Simpson-integrated count)
//             + a stochastic catalogue realisation (Poisson counts, truncated Gutenberg-Richter magnitudes between Mmin and
//             the magnitude of a full-fault rupture with a 3 MPa stress drop, location uniform on the face) - fixed seed
//   two-way   '-- RESSIM: FAULTMECH ... permFactor': a reactivated face multiplies the transmissibility of its connections by
//             permFactor^x, x = clamp(CFF / (w sigma'_n), 0, 1) (w = item 10, default 0.05), irreversible (max over time) - fault-valve / slip-induced
//             permeability; applied after each report step (explicit, one report step lag), via Simulator::setConnectionMultipliers
//   RSF       '-- RESSIM: FAULTRSF name layer a b Dc': after the run a quasi-dynamic rate-and-state simulation (rsf.h, G3) on the
//             fault trace of that layer, loaded by the simulated pressure and stress changes of its faces -> deterministic events
// C25 engines and stress sources: black oil, compositional (CO2 storage) and thermal (EGS) runs. 'STRESSFIELD ... FE': the
// stress CHANGE is the full tensor change of the finite-element geomechanics at the two face cells since t0
// (setStressProvider: GEOMECH sequential coupling of black oil / compositional, the THM mechanics of thermal runs; includes
// poro- and thermo-elastic stress, shear and stress arching) added to the analytic in-situ field (whose SHmax azimuth the
// FE initial state - Sv, k0 Sv - cannot represent); without FE the change is the poroelastic path, plus for thermal runs the
// thermo-elastic horizontal stress dS_H = dS_h = E alpha / (1 - nu) dT ('TH E_GPa alpha nu').
// Two-way coupling is engine neutral (connMultipliers: connection index of the grid -> multiplier) and is composed with the
// geomechanics transmissibility multipliers by the caller.
// Limits: no stress transfer between faces / no slip-induced stress change (the FE stress is that of the continuum), faces
// of LGR cells skipped.
#pragma once
#include <array>
#include <atomic>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ressim {
struct Case;
class Grid;
class Simulator;
namespace faultmech {

struct Face {
    int fault = 0;                         // index into Output::names / Case::faultMech.faults
    int g = -1, a = -1, b = -1;            // global cell, its active index, dominant active neighbour across the face (-1: none)
    int k = 0;                             // 0-based layer
    std::array<float, 3> c{};              // centre (case units: x, y, depth)
    std::array<float, 3> n{};              // unit normal
    std::array<std::array<float, 3>, 4> quad{};
    double area = 0;                       // case length units^2
    std::vector<int> conns;                // connections across the face
    // initial state
    double p0 = 0, sn0 = 0, tau0 = 0, cff0 = 0, aSig = 0, T0 = 0;
    std::array<double, 6> fe0{};           // C25 FE total stress at t0 (case units; xx yy zz xy yz xz)
    bool feOk = false;
};

struct StepOut {
    double t = 0;                          // days (case time units)
    std::vector<float> cff, dcfs, st, sne, rate, cum, mult;   // per face: CFF, dCFS, slip tendency, sigma'_n, R / r, expected events, T multiplier
    std::vector<float> slip, dtauSlip;     // C26 (FAULTSLIP): cumulative slip (m), shear-stress change from slip of all faces (case units)
    double moment = 0, momentStep = 0;     // C26: cumulative / this step seismic moment of the DDM slip (N m)
    int slipping = 0;                      // C26: faces that slipped in this step
};

struct Series {                            // per fault and report step
    std::vector<double> t, maxDcfs, maxCff, maxSt, reactFrac, expected, events, maxMw;
    std::vector<double> maxSlip, moment;   // C26 (FAULTSLIP): maximum cumulative slip (m), cumulative moment (N m)
};

struct Event { double t = 0; float x = 0, y = 0, z = 0; float mw = 0; int fault = 0; int kind = 0; int faces = 0; };   // kind 0 Dieterich (stochastic), 1 RSF, 2 DDM slip event (faces = slipped faces)

struct RsfInfo { int fault = 0, layer = 0, faces = 0; double hStar = 0, Lb = 0, dx = 0, seconds = 0; int events = 0; std::string error; };

struct Output {
    bool ok = false; std::string error;
    std::vector<std::string> names;        // fault names
    std::vector<double> mMax;              // per fault: full-rupture magnitude cap
    std::vector<Face> faces;
    std::vector<StepOut> steps;
    std::vector<Series> series;
    std::vector<Event> events;
    std::vector<RsfInfo> rsf;
    double paPerUnit = 1, mPerUnit = 1;    // case pressure / length unit in Pa / m
    int criticalAtStart = 0;               // faces with CFF >= 0 at t0
    long twoWayUpdates = 0;
    bool feStress = false, thermalStress = false;   // C25 stress source flags of the run
    // C26 DDM stress transfer
    bool slip = false; bool halfSpace = true; double shearModulus = 0;
    long ddmPairs = 0, ddmDropped = 0; double ddmSeconds = 0; int maxSlipIter = 0; std::string slipNote;
    std::string engine = "black oil";
};

class Monitor {
public:
    // false + err when the case has no usable fault mechanics (no STRESSFIELD, no FAULTMECH matching a FAULTS name, ...)
    bool setup(const Case& c, const Grid& g, std::string& err);
    bool active() const { return active_; }
    // call once after initialisation (t = 0) and after every report step: pressure per active cell (case units)
    // C25: temperature per active cell (any unit, only changes are used) for the thermo-elastic stress of thermal runs
    void step(double tDays, const std::vector<float>& p, const std::vector<float>* temperature = nullptr);
    // C25 FE stress source: total stress tensor (Pa, compression +; xx yy zz xy yz xz) of an active cell; false = no element
    using StressFn = std::function<bool(int cell, double S[6])>;
    void setStressProvider(StressFn f) { fe_ = std::move(f); out_.feStress = bool(fe_); }
    bool wantsFeStress() const { return source_ == 1; }
    // C25 engine-neutral two-way coupling: grid connection -> transmissibility multiplier (max over the faces crossing it)
    std::map<int, double> connMultipliers() const;
    // two-way coupling: apply the current face multipliers to the simulator (no-op for one-way faults)
    void applyCoupling(Simulator& sim);
    bool twoWay() const { return twoWay_; }
    // after the run: quasi-dynamic RSF on the FAULTRSF traces (adds kind-1 events)
    void finish(std::atomic<bool>* cancel = nullptr, const std::function<void(const std::string&)>& log = {});
    const Output& out() const { return out_; }
    Output& out() { return out_; }
private:
    bool active_ = false, twoWay_ = false;
    Output out_;
    std::vector<double> gamma_, mult_;     // per face
    std::vector<double> mu_, coh_, r0_, ta_, perm_, bval_, mmin_, ramp_;   // per fault
    double svGrad_ = 0, sv0_ = 0, rh_ = 0.7, rH_ = 0.9, az_ = 0, path_ = 0.5;
    std::vector<int> rsfLayer_; std::vector<double> rsfA_, rsfB_, rsfDc_;
    unsigned seed_ = 20261021u;
    int source_ = 0; double thCoef_ = 0;   // C25: stress source, thermo-elastic coefficient (case pressure units per degree)
    StressFn fe_;
    // C26 DDM stress transfer / quasi-static slip
    std::shared_ptr<void> ddm_;            // ddm::ShearMatrix (opaque here)
    std::vector<double> slip_;             // cumulative slip per face (m)
    std::vector<std::array<double, 3>> rake_;   // shear direction at t0 (case frame, unit)
    bool slipOn_ = false; double slipNu_ = 0.25, slipCut_ = 25; int slipIterMax_ = 4000;
    double dMu_ = 0, slipRef_ = 0.01;
    std::vector<char> failed_;             // weakening: face has failed (dynamic friction from then on)
    void setupSlip(const Case& c);
    void slipStep(const std::vector<double>& pf, const std::vector<double>& snR, const std::vector<std::array<double, 3>>& shR,
                  std::vector<double>& sn, std::vector<double>& tau, StepOut& S, double tDays, double tPrev);
    void stress(const Face& f, double p, double dT, const double* dS6, double& sn, double& tau, double* shear = nullptr) const;
};

// C47: per-face history CSV (t, face, fault, k, x, y, depth, area, aSig, cff, dcfs, rate, cum); also written by finish() when
// the environment variable RESSIM_FAULT_CSV names a file (CLI: --fault-csv file)
bool writeFaultCsv(const Output& o, const std::string& path);
}  // namespace faultmech
}  // namespace ressim
