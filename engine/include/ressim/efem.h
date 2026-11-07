// C20 geomechanics G1: embedded finite elements (EFEM) for fractures with frictional contact - 2-D plane strain core.
//
// DEFAULT (jumpModel 0): continuous jumps by nodal Heaviside enrichment (u = sum N_a u_a + sum N_a (H - H_a) b_a,
// Moes et al. 1999 / Dolbow et al. 2001 contact), integration on the two sub-polygons of every cut element (exact for
// the bilinear products), contact by penalty + Coulomb return mapping at 2 Gauss points per embedded segment.
// Research finding F9 (research/FINDINGS_LOG.md): the classic element-constant EFEM below is accurate on aligned
// meshes but its slip error decays only like O(h) when the fracture crosses the mesh obliquely (+36 % at a/h = 8).
// Element-constant variant (jumpModel 1; strong discontinuity embedded in bilinear quadrilaterals, one jump per cut element;
// Oliver 1996, Jirasek 2000, Borja 2000, Cusini et al. 2021 / GEOS):
//   u = sum_a N_a u_a + (H_G - phi) [[u]],   phi = sum_{a in Omega+} N_a,   [[u]] = w_n n + w_t m
// so the nodal values are the true displacements and, away from the fracture, eps = B u - sym(grad phi (x) [[u]]).
// Statically and kinematically optimal non-symmetric formulation (SKON, the most accurate variant in Jirasek's
// comparison): the element equation for the jump enforces traction continuity in the mean,
//   (1 / A_e) int_e R^T sigma n dA = t(w) - p e_n        (R = [n m], tension positive, p = fluid pressure in the crack)
// with the contact traction t(w): open (w_n >= 0) t = 0; closed: penalty t_n = K_n w_n, Coulomb friction with return
// mapping |t_t| <= c + mu |t_n| (stick: t_t = K_t (w_t - w_t^p); slip: t_t = sign (c + mu |t_n|), slip history w_t^p).
// The jumps are condensed element by element; the global displacement system is solved by a banded LU and a Newton /
// active-set iteration on the contact states.
// Boundary conditions: rollers on the left (u_x = 0) and bottom (u_y = 0), normal tractions sxx on the right and syy on
// the top (far-field total stress, tension positive) - the uncracked solution is the uniform stress (sxx, syy).
// Verified against Sneddon's pressurised crack and frictional slip on an inclined crack under compression
// (tests/c20_tests.inc).
#pragma once
#include <string>
#include <vector>

namespace ressim {
namespace efem {

struct Fracture {
    double x0 = 0, y0 = 0, x1 = 1, y1 = 0;   // straight segment (m); end points preferably on element edges
    double pressure = 0;                     // fluid pressure in the fracture (Pa)
    // C21 G1-c: optional pressure profile along the fracture (s = signed distance from the midpoint towards (x1, y1), m;
    // Pa), linear in between, constant beyond the ends; overrides 'pressure' when it has >= 2 points
    std::vector<double> profS, profP;
    double friction = 0.6, cohesion = 0;     // Coulomb (-, Pa)
};

struct Input {
    double Lx = 10, Ly = 10;                 // domain (m)
    int nx = 80, ny = 80;                    // bilinear elements
    double E = 10e9, nu = 0.25;              // Pa, -
    double sxx = 0, syy = 0;                 // far-field total normal stresses on the right / top faces (Pa, tension +)
    std::vector<Fracture> fracs;
    double penalty = 1e3;                    // contact stiffness K_n = K_t = penalty E / h
    // jump interpolation: 0 continuous (default): Heaviside enrichment of the nodes whose support the fracture cuts,
    //   [[u]] = sum_a N_a b_a along the fracture (piecewise linear, zero at the tips), contact at 2 Gauss points per cut
    //   element - converges like the regular field; 1 element-constant jumps (classic EFEM, SKON, condensed per element)
    //   - kept for comparison: slip / opening errors only decay like O(h) on meshes not aligned with the fracture
    int jumpModel = 0;
    int maxNewton = 60;
    double tol = 1e-9;                       // relative residual
};

struct CutElement {
    int elem = 0, frac = 0;
    double xc = 0, yc = 0, s = 0, len = 0;   // centre of the embedded segment, signed distance along the fracture from its midpoint, length
    double nx = 0, ny = 0;                   // unit normal n = (-m_y, m_x), m = unit tangent from (x0, y0) to (x1, y1)
    double wn = 0, wt = 0;                   // normal opening and tangential slip (m)
    double tn = 0, tt = 0;                   // effective contact tractions (Pa, tension +)
    int state = 0;                           // 0 open, 1 stick, 2 slip
    double wtp = 0;                          // committed plastic slip
};

struct Result {
    bool ok = false; std::string error;
    int newton = 0, dofs = 0;
    std::vector<double> u;                   // nodal displacements (2 per node, node = i + (nx + 1) j)
    std::vector<CutElement> cut;
    std::vector<double> sxx, syy, sxy;       // element-mean stresses (Pa, tension +)
    double seconds = 0;
};

// Solve the static problem (one load step from the unloaded state, or from 'prev' cut-element slip history)
bool solve(const Input& in, Result& out, const Result* prev = nullptr);

}  // namespace efem
}  // namespace ressim
