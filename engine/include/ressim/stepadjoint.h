#pragma once
// C18: engine-independent discrete adjoint of a fully implicit time-step sequence
//     R_n(x_n, x_{n-1}, theta) = 0,  n = 1..N,     F = F(x_1..x_N, theta)
// Backward recursion (one transposed linear solve per step, independent of the number of parameters):
//     J_n^T lam_n = -( dF/dx_n^T + (dR_{n+1}/dx_n)^T lam_{n+1} ),   J_n = dR_n/dx_n at the converged state
//     dF/dtheta  = dF/dtheta|_x + sum_n lam_n^T dR_n/dtheta + lam_1^T dR_1/dx_0 dx_0/dtheta
// The engines (compositional, thermal) supply the step Jacobians, the coupling to the previous state, the objective
// gradient and the parameter terms as callbacks.
#include "ressim/bsr.h"
#include <functional>
#include <string>
#include <vector>

namespace ressim {

struct StepAdjointProblem {
    int steps = 0;
    std::function<void(int n, BSR& J)> jacobian;                                                    // dR_n/dx_n
    std::function<void(int n, std::vector<double>& rhs)> addObjectiveGradient;                      // rhs += dF/dx_n
    std::function<void(int n, const std::vector<double>& lam, std::vector<double>& rhs)> addPrevCoupling; // rhs += (dR_n/dx_{n-1})^T lam
    std::function<void(int n, const std::vector<double>& lam, std::vector<double>& grad)> addParamTerm;   // grad += lam^T dR_n/dtheta
    std::function<void(const std::vector<double>& lam1, std::vector<double>& grad)> addInitialTerm;        // grad += lam_1^T dR_1/dx_0 dx_0/dtheta
};

struct StepAdjointStats { bool ok = true; int linear = 0; double worstRel = 0; std::string error; };

// grad must be pre-sized to the number of parameters (and may already hold dF/dtheta|_x)
StepAdjointStats solveStepAdjoint(const StepAdjointProblem& P, std::vector<double>& grad);

// x = J^-T b  (explicit block transpose + BILU0 + GMRES)
LinResult solveTransposed(const BSR& J, const std::vector<double>& b, std::vector<double>& x, double tol = 1e-11, int maxIt = 800);
BSR transposeBsr(const BSR& J);

} // namespace ressim
