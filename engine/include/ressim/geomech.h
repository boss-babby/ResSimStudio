// Geomechanics (milestone C8): linear thermo-poroelastic finite elements for THM / THMC simulation.
//
// Model (quasi-static, small strain, tension positive):
//   div sigma = 0,   sigma = D eps - (alpha dp + 3 K_dr beta dT) m,   eps = sym grad u,   m = (1 1 1 0 0 0)
// with incremental pressure / temperature changes dp, dT from the initial (equilibrium) state, drained isotropic
// stiffness D(E, nu), Biot coefficient alpha, linear thermal expansion beta and K_dr = E / (3 (1 - 2 nu)).
// Discretisation: trilinear hexahedra (8 nodes, 2x2x2 Gauss) on the Cartesian flow grid, optionally padded with
// non-reservoir over / underburden layers and side-burden columns. Boundary conditions: rollers on the sides and the
// base (normal displacement 0), free or fixed top surface; 2-D sections (ny = 1 or nx = 1) are plane strain. The SPD
// system is solved by preconditioned conjugate gradients (symmetric Gauss-Seidel / SSOR), warm-started from the
// previous displacement field.
//
// Coupling to flow (Kim, Tchelepi & Juanes 2011 fixed-stress split; unconditionally stable): Lagrangian porosity
//   phi = phi0 + alpha eps_v + dp / N - beta_phi dT,   1 / N = (alpha - phi0) / K_s,   beta_phi = 3 beta (alpha - phi0)
// and, at a fixed mean total stress, the flow sees  d phi = c_p dp + c_T dT  with
//   c_p = alpha^2 / K_dr + 1 / N,   c_T = 3 alpha beta - beta_phi.
// Stresses (compression positive for reporting): sigma = sigma_0 - d sigma, sigma_0 from an overburden gradient and a
// horizontal / vertical effective-stress ratio; effective stress sigma' = sigma - p; Mohr-Coulomb shear failure index
// tau_m / (c cos phi_f + sigma'_m sin phi_f) and tensile failure (sigma'_3 < -T_0).
//
// C16 elastoplasticity (plasticity = 1): Drucker-Prager yield surface matched to Mohr-Coulomb (outer / inner cone or
// plane strain), F = sqrt(J2) - eta p' - xi c(eps_p) <= 0 (p' mean effective stress, compression +), non-associated flow
// potential G = sqrt(J2) - eta_psi p', linear cohesion hardening / softening c = max(c_res, c0 + H eps_p) and the
// closed-form return mapping to the cone or its apex (de Souza Neto, Peric & Owen 2008, ch. 8). Global equilibrium by
// initial-stiffness (modified Newton) iterations with the elastic PCG; the plastic state is path dependent and committed
// at accepted time steps (commit()). Damage-permeability coupling k / k0 = min(k_max, exp(b eps_p)) (shear-enhanced
// permeability of fault / damage zones). Initial stresses violating the yield surface are made admissible by a local
// cohesion increase (the initial state is in equilibrium by definition).
// Rigid plates (bcType 3): the normal displacement of a face is uniform, either prescribed (plateDispControl) or found
// so that the net normal force increment equals bcStress x face area (secant iteration) - Mandel's problem, platen tests.
#pragma once
#include <string>
#include <vector>

namespace ressim {
namespace geomech {

struct Input {
    bool enabled = false;
    int coupling = 1;                       // 0 one-way (mechanics after each step, porosity lagged), 1 iterative fixed-stress
    std::vector<double> E{5e9}, nu{0.25}, biot{1.0}, alphaT{1e-5};   // reservoir cells (size 1 = uniform): Pa, -, -, 1/K
    double Ks = 0;                          // grain bulk modulus (Pa); <= 0: incompressible grains (1 / N = 0)
    int overLayers = 0, underLayers = 0, sidePad = 0;      // padding elements above / below / around the reservoir
    double overThick = 300, underThick = 300, sideThick = 300;   // m (total thickness of each padding zone)
    double padE = 5e9, padNu = 0.25, padAlphaT = 1e-5;     // padding rock
    bool topFree = true;                    // free (traction-free) top surface; false: fixed vertical displacement
    // boundary conditions of the six faces of the mechanics grid (OPM BCMECH-like): 0 x-, 1 x+, 2 y-, 3 y+, 4 top, 5 base
    //   type 0 roller (normal displacement fixed), 1 free (optionally loaded by bcStress), 2 fixed (all components);
    //   -1 = default (sides / base roller, top from topFree). bcStress: normal stress increment on a free face (Pa,
    //   compression positive) applied as consistent nodal loads - tectonic / surface load changes
    int bcType[6] = {-1, -1, -1, -1, -1, -1};
    double bcStress[6] = {0, 0, 0, 0, 0, 0};
    double sigmaVGrad = 22.6e3;             // overburden stress gradient (Pa/m)
    double k0 = 0.8;                        // horizontal / vertical effective stress ratio of the initial state
    double cohesion = 2e6, frictionDeg = 30, tensile = 1e6;   // Mohr-Coulomb (Pa, degrees, Pa)
    // C16 plasticity and rigid plates (see the header comment)
    int plasticity = 0;                     // 0 linear elastic, 1 Drucker-Prager matched to Mohr-Coulomb (cohesion, frictionDeg)
    int dpMatch = 2;                        // 0 outer cone (triaxial compression), 1 inner cone (extension), 2 plane strain
    double dilationDeg = 0;                 // dilation angle psi (<= frictionDeg; = frictionDeg: associated flow)
    double hardening = 0;                   // cohesion hardening (+) / softening (-) modulus H (Pa per unit eq. plastic strain)
    double residualCohesion = 0;            // softening floor (Pa)
    double permPlastic = 0;                 // damage permeability k / k0 = exp(permPlastic eps_p) (0: off)
    double permPlasticMax = 100;            // cap of that multiplier
    int plasticIter = 80;                   // equilibrium iterations per mechanics solve
    double plasticTol = 1e-6;               // relative force residual
    bool plateDispControl[6] = {false, false, false, false, false, false};   // bcType 3: displacement instead of force control
    double bcDisp[6] = {0, 0, 0, 0, 0, 0};  // prescribed plate displacement into the body (m)
    int maxIter = 4;                        // fixed-stress iterations per time step
    double tol = 1e-4;                      // porosity change for convergence of the fixed-stress iteration
    // C25: stress-dependent permeability of the matrix k / k0 = exp(-permStress (sigma'_m - sigma'_m0)) (1/Pa, 0 = off),
    // multiplied with the damage multiplier; capped to [1e-4, permPlasticMax]
    double permStress = 0;
};

class Solver {
public:
    // Cartesian reservoir box: sizes (size 1 = uniform, else per cell), top depth, initial pressure (Pa) / porosity per cell
    bool setup(const Input& in, int nx, int ny, int nz, const std::vector<double>& dx, const std::vector<double>& dy,
               const std::vector<double>& dz, double top, const std::vector<double>& p0, const std::vector<double>& phi0,
               std::string* error = nullptr);
    // C25 corner-point / general hexahedral mesh: the reservoir nodes are given explicitly ((nx+1)(ny+1)(nz+1) x 3
    // coordinates x, y, depth in m, node (i, j, k) at index ((k (ny+1) + j) (nx+1) + i)); trilinear isoparametric elements
    // with the full Jacobian. Padding nodes are the boundary nodes shifted by the padding thicknesses (vertical / horizontal).
    bool setupHex(const Input& in, int nx, int ny, int nz, const std::vector<double>& nodes, const std::vector<double>& p0,
                  const std::vector<double>& phi0, std::string* error = nullptr);
    bool generalHex() const { return !X_.empty(); }
    bool ready() const { return ready_; }
    // displacements for pressure / temperature changes (per reservoir cell, from the initial state)
    bool solve(const std::vector<double>& dp, const std::vector<double>& dT);
    // per reservoir cell results of the last solve
    const std::vector<double>& volStrain() const { return ev_; }
    const std::vector<double>& porosity() const { return phi_; }          // Lagrangian porosity
    const std::vector<double>& dispZ() const { return uzc_; }             // cell-centre vertical displacement (m, + down)
    const std::vector<double>& dispX() const { return uxc_; }
    const std::vector<double>& dispY() const { return uyc_; }
    const std::vector<double>& sigEffMean() const { return sm_; }         // mean effective stress (Pa, compression +)
    const std::vector<double>& sigEffV() const { return szz_; }           // vertical effective stress (Pa, compression +)
    const std::vector<double>& sigEffH() const { return sxx_; }           // x effective stress (Pa, compression +)
    const std::vector<double>& failure() const { return fail_; }          // max(shear ratio, tensile ratio); >= 1 failed
    // C25: full effective stress tensor per reservoir cell (Pa, compression +; xx, yy, zz, xy, yz, xz) and the pore
    // pressure of the last solve; total stress = effective + p on the diagonal
    const std::vector<double>& sigEffFull() const { return sfull_; }
    void totalStress(int cell, double S[6]) const;
    double minTotalStress(int cell) const;                                // smallest principal total stress (Pa, compression +)
    double surfaceUplift() const { return uplift_; }                       // max upward top-surface displacement (m)
    double reservoirTopUplift() const { return upliftRes_; }               // max upward displacement of the reservoir top (m)
    // fixed-stress coefficients per bulk volume (1/Pa, 1/K) and reference porosity
    double cp(int cell) const { return cp_[size_t(cell)]; }
    double cT(int cell) const { return cT_[size_t(cell)]; }
    double phi0(int cell) const { return phi0_[size_t(cell)]; }
    int lastCg() const { return lastCg_; }
    // C16 plasticity / plates
    void commit();                                                         // accept the plastic state of the last solve
    void setFaceLoad(int face, double stress, double disp) { if (face >= 0 && face < 6) { in_.bcStress[face] = stress; in_.bcDisp[face] = disp; bcStress_[face] = stress; } }
    const std::vector<double>& plasticStrain() const { return epc_; }     // cell mean equivalent plastic strain (last solve)
    const std::vector<double>& plasticVolStrain() const { return evp_; }  // cell plastic volumetric strain (+ dilation)
    const std::vector<double>& permMultiplier() const { return kmult_; }  // damage permeability multiplier (last solve)
    int plasticPoints() const { return nPlast_; }                         // yielding Gauss points in the last solve
    int lastPlasticIter() const { return lastPIt_; }
    bool lastConverged() const { return lastConv_; }
    double plateDisp(int f) const { return plateU_[size_t(f)]; }          // m into the body
    double plateForce(int f) const { return plateF_[size_t(f)]; }         // net normal force increment (N, compression +)
    double faceArea(int f) const;
    // C14 restart: the state and its derived fields (C16: + committed plastic strains and plate displacements)
    template <class Ar> void stateIO(Ar& ar) { ar(u_, ev_, phi_, uxc_, uyc_, uzc_, sm_, szz_, sxx_, fail_, uplift_, upliftRes_, epsPc_, ebarC_, plateU_); }
    int dofs() const { return int(K_.n); }
private:
    struct Csr { size_t n = 0; std::vector<int> ptr, col; std::vector<double> val; std::vector<int> diag; };
    Input in_;
    bool ready_ = false;
    int nx_ = 0, ny_ = 0, nz_ = 0;          // reservoir
    int NX_ = 0, NY_ = 0, NZ_ = 0;          // mechanics grid (elements)
    int ox_ = 0, oy_ = 0, oz_ = 0;          // reservoir offset in the mechanics grid
    std::vector<double> ex_, ey_, ez_;      // element sizes along x, y, z (per column / row / layer)
    std::vector<double> eE_, eNu_, eA_, eB_;   // per element: E, nu, alpha, beta
    std::vector<int> resCell_;              // element -> reservoir cell (-1: padding)
    std::vector<double> sx0_, depthC_;      // initial stresses / centre depths per reservoir cell
    std::vector<double> p0c_, phi0_, cp_, cT_, invN_, betaPhi_;
    std::vector<char> fixed_;               // constrained dofs
    int bc_[6] = {0, 0, 0, 0, 1, 0};        // resolved face types
    double bcStress_[6] = {0, 0, 0, 0, 0, 0};
    Csr K_;
    std::vector<double> u_, f_;
    std::vector<double> ev_, phi_, uxc_, uyc_, uzc_, sm_, szz_, sxx_, fail_;
    std::vector<double> sfull_, pAbs_, sm0_;   // C25: full effective stress (6 per cell), pore pressure, initial mean effective stress
    std::vector<double> X_;                 // C25 general hexahedra: node coordinates (3 per mechanics node; empty = tensor box)
    bool setupImpl(const Input& in, int nx, int ny, int nz, const std::vector<double>& dx, const std::vector<double>& dy,
                   const std::vector<double>& dz, double top, const std::vector<double>& p0, const std::vector<double>& phi0,
                   std::string* error, const std::vector<double>* resNodes);
    double elemGrad(int i, int j, int k, double xi, double eta, double zeta, double G[8][3]) const;   // returns det J
    double uplift_ = 0, upliftRes_ = 0;
    int lastCg_ = 0;
    // C16 plasticity state per Gauss point (8 per element): committed / trial plastic strain (6, tension +, engineering
    // shear), equivalent plastic strain; per element cohesion adjustment for an admissible initial state
    std::vector<double> epsPc_, epsPt_, ebarC_, ebarT_, cAdj_, sgp_;
    std::vector<double> epc_, evp_, kmult_;
    std::vector<double> plateU_ = std::vector<double>(6, 0.0), plateF_ = std::vector<double>(6, 0.0);
    std::vector<std::vector<int>> plateDofs_ = std::vector<std::vector<int>>(6);
    int nPlast_ = 0, lastPIt_ = 0; bool lastConv_ = true;
    double eta_ = 0, xi_ = 0, etaB_ = 0;
    void internalForce(const std::vector<double>& dp, const std::vector<double>& dT, std::vector<double>& fint, bool storeGp);
    bool equilibrate(const std::vector<double>& dp, const std::vector<double>& dT, const std::vector<double>& fext);
    void cellResults(const std::vector<double>& dp, const std::vector<double>& dT, bool fromGp);
    int node(int i, int j, int k) const { return (k * (NY_ + 1) + j) * (NX_ + 1) + i; }
    int elem(int i, int j, int k) const { return (k * NY_ + j) * NX_ + i; }
    void elemNodes(int i, int j, int k, int nd[8]) const;
    int pcg(const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt);
};

// isotropic drained stiffness helpers
inline double bulkModulus(double E, double nu) { return E / (3 * (1 - 2 * nu)); }
inline double shearModulus(double E, double nu) { return E / (2 * (1 + nu)); }

} // namespace geomech
} // namespace ressim
