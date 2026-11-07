// C18: generic discrete adjoint driver (see stepadjoint.h)
#include "ressim/stepadjoint.h"
#include <algorithm>
#include <cmath>
#include <string>

namespace ressim {

BSR transposeBsr(const BSR& J) {
    BSR T; T.n = J.n; T.nb = J.nb;
    const int n = J.n, nb = J.nb, b2 = nb * nb;
    std::vector<int> cnt(size_t(n) + 1, 0);
    for (int k = 0; k < J.nnzb(); ++k) cnt[size_t(J.col[size_t(k)]) + 1]++;
    for (int i = 0; i < n; ++i) cnt[size_t(i) + 1] += cnt[size_t(i)];
    T.rowPtr = cnt; T.col.assign(size_t(J.nnzb()), 0); T.val.assign(J.val.size(), 0.0);
    std::vector<int> pos(cnt.begin(), cnt.end() - 1);
    for (int r = 0; r < n; ++r)
        for (int k = J.rowPtr[size_t(r)]; k < J.rowPtr[size_t(r) + 1]; ++k) {
            const int c = J.col[size_t(k)], d = pos[size_t(c)]++;
            T.col[size_t(d)] = r;   // rows visited in increasing order -> columns of T sorted
            const double* B = J.blk(k); double* D = T.blk(d);
            for (int i = 0; i < nb; ++i) for (int j = 0; j < nb; ++j) D[j * nb + i] = B[i * nb + j];
        }
    (void)b2;
    T.finalizeDiag();
    return T;
}

LinResult solveTransposed(const BSR& J, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt) {
    const BSR T = transposeBsr(J);
    BILU0 M; M.factor(T);
    BsrOperator op(T);
    x.assign(b.size(), 0.0);
    double bn = 0; for (double v : b) bn += v * v;
    if (bn == 0) { LinResult r; r.converged = true; return r; }
    LinResult r = gmres(op, M, b, x, tol, maxIt, 80);
    if (!r.converged) {   // second chance from the current iterate with BiCGSTAB
        LinResult r2 = bicgstab(op, M, b, x, tol, maxIt);
        r2.iterations += r.iterations;
        return r2;
    }
    return r;
}

StepAdjointStats solveStepAdjoint(const StepAdjointProblem& P, std::vector<double>& grad) {
    StepAdjointStats st;
    std::vector<double> lamNext, lam, rhs, neg;
    BSR J;
    for (int n = P.steps; n >= 1; --n) {
        P.jacobian(n, J);
        rhs.assign(size_t(J.size()), 0.0);
        if (P.addObjectiveGradient) P.addObjectiveGradient(n, rhs);
        if (n < P.steps && P.addPrevCoupling) P.addPrevCoupling(n + 1, lamNext, rhs);
        neg.resize(rhs.size());
        for (size_t i = 0; i < rhs.size(); ++i) neg[i] = -rhs[i];
        const LinResult r = solveTransposed(J, neg, lam);
        st.linear += r.iterations;
        if (!r.converged) { st.worstRel = std::max(st.worstRel, r.relResidual); if (r.relResidual > 1e-6) { st.ok = false; st.error = "adjoint solve did not converge at step " + std::to_string(n); } }
        if (P.addParamTerm) P.addParamTerm(n, lam, grad);
        lamNext.swap(lam);
    }
    if (P.steps >= 1 && P.addInitialTerm) P.addInitialTerm(lamNext, grad);
    return st;
}

} // namespace ressim
