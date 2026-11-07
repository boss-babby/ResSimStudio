// C26: 3-D displacement discontinuity method (DDM) on triangulated surfaces - the elastic interaction between fault faces
// (stress transfer: slip on one face loads its neighbours) and the opening of pressurised fractures (mechanical apertures).
//
// Kernel   triangular dislocation elements (TDE) of Nikkhoo & Walter (2015), full space or homogeneous half space with a
//          traction-free surface at z = 0 (z up, z <= 0 in the rock), isotropic linear elasticity (G, nu). Constant
//          displacement discontinuity per element; an element is a planar or slightly warped polygon split into triangles
//          (corner-point cell faces -> 2 triangles). Collocation at the element centroid (receivers exactly on a triangle
//          edge line are moved by 1e-6 of the element size inside the element plane: the edge singularities of the two
//          triangles of one element cancel, so the result is the regular limit).
// Shear    each element slips along a fixed unit direction d_i in its plane (the rake = direction of the initial resolved
//          shear traction). K_tau(i,j) = shear-traction change along d_i at the centroid of i per metre of slip of j,
//          K_n(i,j) = normal-stress change (compression +) at i. The sense of d_j is chosen so that slip relieves the
//          element's own shear (K_tau(j,j) < 0), which makes the matrix physically oriented whatever the face orientation.
//          Pairs further apart than cutoff x (size_i + size_j) / 2 are dropped (stress decays as r^-3; default 25).
// Opening  K_o(i,j) = normal traction (tension +) at i per metre of opening of j (dense; fractures have few elements).
// Solvers  quasi-static frictional slip: find monotone slip increments ds >= 0 with
//              tau_i + sum_j (K_tau - mu_i K_n)(i,j) ds_j - mu_i sigma'_i - c_i = 0   where ds_i > 0,   <= 0 elsewhere
//          (linear complementarity, projected Gauss-Seidel with over-relaxation, exact for a single element);
//          fracture opening under net pressure: K_o w = -p_net with contact w >= 0 (same solver).
// Verification (tests -c26): strain vs finite differences of the displacement, Betti reciprocity, Okada-type symmetry,
// penny-shaped crack volume (Sneddon) and opening profile, single-element stress drop, frictional equilibrium after slip.
#pragma once
#include <array>
#include <string>
#include <vector>

namespace ressim {
namespace ddm {

using V3 = std::array<double, 3>;
struct Tri { V3 a, b, c; };

// stress (Pa, tension +; xx yy zz xy xz yz) at obs from a constant displacement discontinuity (m) on a triangle, given in the
// triangle's own frame (strike, dip, tensile) - see triFrame. halfSpace: free surface z = 0, all points must have z <= 0.
bool triStress(const V3& obs, const Tri& t, const V3& slipSDT, double G, double nu, bool halfSpace, double S[6]);
bool triDisp(const V3& obs, const Tri& t, const V3& slipSDT, double nu, bool halfSpace, double u[3]);
bool triStrain(const V3& obs, const Tri& t, const V3& slipSDT, double nu, bool halfSpace, double e[6]);
// the TDE frame of a triangle: unit normal, strike and dip vectors (a slip vector b has components (b.s, b.d, b.n))
void triFrame(const Tri& t, V3& n, V3& s, V3& d);

struct Element {
    std::vector<Tri> tris;   // metres, z up
    V3 c{}, n{};             // centroid, unit normal (all triangles are oriented along n)
    double area = 0, size = 0;
};
// build an element from a quadrilateral (corner loop order) - the two triangles get the orientation of the quad normal
Element quadElement(const V3& p0, const V3& p1, const V3& p2, const V3& p3);
// stress at obs from a uniform discontinuity vector b (global frame, m) on the element (sum over its triangles)
bool elementStress(const V3& obs, const Element& e, const V3& b, double G, double nu, bool halfSpace, double S[6]);

struct ShearMatrix {
    int n = 0;
    std::vector<int> start, col;     // CSR per receiver row
    std::vector<double> kt, kn;      // Pa per m
    std::vector<V3> dir;             // shear direction d of each element (in plane, unit); slip along d relieves it
    double sense = 1;                // TDE discontinuity b = sense * slip * d
    std::vector<double> selfKt;      // K_tau(i,i)
    long pairs = 0, dropped = 0;
    double seconds = 0;
};
// dir: rake per element (unit, in plane; any sign). cutoff <= 0: keep all pairs.
ShearMatrix buildShear(const std::vector<Element>& el, const std::vector<V3>& dir, double G, double nu, bool halfSpace,
                       double cutoff = 25.0, int threads = 0);

// dense opening matrix (row-major n x n): normal traction (tension +) at i per metre of opening of j
std::vector<double> buildOpening(const std::vector<Element>& el, double G, double nu, bool halfSpace);

struct SlipResult { int iterations = 0, slipping = 0; double residual = 0; bool converged = false; };
// excess[i] = tau_i - mu_i sigma'_i - c_i at the current slip (Pa); ds (in/out, m) = slip increments of this solve (>= 0);
// mu[i] = friction used for the K_n coupling and for the slipping elements (dynamic friction for weakening faults)
SlipResult solveSlip(const ShearMatrix& K, const std::vector<double>& excess, const std::vector<double>& mu, std::vector<double>& ds,
                     double tolPa = 1.0, int maxIter = 5000);
// opening w >= 0 with K_o w = -pNet (Pa); returns the residual of the open elements
SlipResult solveOpening(const std::vector<double>& Ko, int n, const std::vector<double>& pNet, std::vector<double>& w,
                        double tolPa = 1.0, int maxIter = 20000);

// analytic references used by the tests and the docs
double pennyVolume(double p, double a, double E, double nu);          // 16 (1 - nu^2) p a^3 / (3 E)
double pennyOpening(double p, double a, double r, double E, double nu);  // 8 (1 - nu^2) p sqrt(a^2 - r^2) / (pi E)

}  // namespace ddm
}  // namespace ressim
