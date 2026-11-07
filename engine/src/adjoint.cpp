// Discrete adjoint and forward sensitivities of a history-matching misfit for the black-oil FIM simulator (C9).
//
// Forward model (one fully implicit time step n, scaled residual of simulator.cpp):
//     R_n(x_n, x_{n-1}, theta) = [acc(x_n) - acc(x_{n-1})] / PVref + dt flux(x_n) / PVref   (+ well equations) = 0
// Misfit F = sum_i ((y_i - d_i) / sigma_i)^2, y_i = observation operators on the per-step summary (time interpolation
// exactly as study.cpp simulatedValue()).
//
// Adjoint (reverse) mode: one transposed solve per time step, independent of the number of parameters
//     J_n^T lam_n = -dF/dx_n^T - (dR_{n+1}/dx_n)^T lam_{n+1},   dR_{n+1}/dx_n = -d acc(x_n)/dx_n / PVref (block diagonal)
//     dF/dtheta  = dF/dtheta|_x + sum_n lam_n^T dR_n/dtheta
// Forward (tangent) mode: one solve per parameter and step, gives the full sensitivity matrix dy/dtheta
//     J_n s_n = -dR_n/dtheta - dR_n/dx_{n-1} s_{n-1}
//
// J_n is the exact AD Jacobian of the Newton solver (assemble()) re-evaluated at the recorded converged state.
// dR_n/dtheta and dy/dtheta|_x are central differences of the residual / observation operator of two perturbed models
// evaluated at the SAME recorded states (no extra simulation; the residual is smooth in theta). The dependence of the
// initial equilibrium state on theta (contacts, datum pressure, Pc) enters through step 1: the perturbed models use their
// own initial state as x_0.
// Explicit well heads (computeHeads at the start of each step) are an augmented lagged state d_n = H(x_{n-1}, bhp_{n-1},
// d_{n-1}; theta) with multipliers mu_n:
//     mu_n = H_{n+1,d}^T mu_{n+1} - R_{n,d}^T lam_n - dF/dd_n,   J_n^T lam_n gains + H_{n+1,x}^T mu_{n+1},
//     dF/dtheta gains - sum_n mu_n^T H_{n,theta}       (H Jacobians: central differences per well, small)
// Other lagged data - group / THP targets, RESV conversion, hysteresis history, DRSDT limits, aquifer pressures - are
// frozen at the recorded values (frozen-coefficient discrete adjoint, cf. OPM / MRST adjoint).
#include "ressim/simulator.h"
#include "ressim/ad.h"
#include <cmath>
#include <chrono>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <thread>
#include <mutex>
#include <atomic>
#include <functional>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <cstdlib>

namespace ressim {

namespace {


bool matchName(const std::string& pattern, const std::string& name) {
    if (!pattern.empty() && pattern.back() == '*')
        return name.compare(0, pattern.size() - 1, pattern, 0, pattern.size() - 1) == 0;
    return pattern == name;
}

bool obsIsWellQ(ObsQty q) {
    return q == ObsQty::WOPR || q == ObsQty::WWPR || q == ObsQty::WGPR || q == ObsQty::WWIR || q == ObsQty::WBHP ||
           q == ObsQty::WWCT || q == ObsQty::WGOR;
}

// identical to study.cpp pointValue()
double pointVal(const SummaryPoint& s, ObsQty q, int w) {
    if (obsIsWellQ(q)) {
        if (w < 0 || w >= int(s.wells.size())) return NAN;
        const WellRates& r = s.wells[size_t(w)];
        switch (q) {
            case ObsQty::WOPR: return r.oil;
            case ObsQty::WWPR: return r.water;
            case ObsQty::WGPR: return r.gas;
            case ObsQty::WWIR: return r.waterInj;
            case ObsQty::WBHP: return r.bhp;
            case ObsQty::WWCT: return (r.oil + r.water) > 0 ? r.water / (r.oil + r.water) : 0.0;
            case ObsQty::WGOR: return r.oil > 0 ? r.gas / r.oil : 0.0;
            default: return NAN;
        }
    }
    switch (q) {
        case ObsQty::FOPR: return s.fopr;
        case ObsQty::FWPR: return s.fwpr;
        case ObsQty::FGPR: return s.fgpr;
        case ObsQty::FPR: return s.fpr;
        case ObsQty::FWCT: return s.fwct;
        case ObsQty::FGOR: return s.fgor;
        default: return NAN;
    }
}

void transposePattern(const CSR& A, CSR& T, std::vector<int>& map) {
    const int n = A.n;
    T.n = n;
    T.rowPtr.assign(size_t(n) + 1, 0);
    for (int q = 0; q < A.rowPtr[n]; ++q) T.rowPtr[size_t(A.col[size_t(q)]) + 1]++;
    for (int i = 0; i < n; ++i) T.rowPtr[size_t(i) + 1] += T.rowPtr[size_t(i)];
    T.col.assign(size_t(A.rowPtr[n]), 0);
    T.val.assign(size_t(A.rowPtr[n]), 0.0);
    map.assign(size_t(A.rowPtr[n]), 0);
    std::vector<int> fill(T.rowPtr.begin(), T.rowPtr.end() - 1);
    for (int i = 0; i < n; ++i)
        for (int q = A.rowPtr[size_t(i)]; q < A.rowPtr[size_t(i) + 1]; ++q) {
            const int c = A.col[size_t(q)];
            const int pos = fill[size_t(c)]++;
            T.col[size_t(pos)] = i;
            map[size_t(q)] = pos;
        }
    T.finalizeDiag();
}

double norm2(const std::vector<double>& v) { double s = 0; for (double x : v) s += x * x; return std::sqrt(s); }

using Sparse = std::vector<std::pair<int, double>>;

// runs fn(j) for j in [0, n) on up to hardware-concurrency std::threads (each with a 1-thread OpenMP team, so the
// replay models do not oversubscribe the cores); the first exception is re-thrown on the caller's thread
void parallelFor(int n, const std::function<void(int)>& fn) {
    int nt = int(std::thread::hardware_concurrency()); if (nt < 1) nt = 1;
    nt = std::min(nt, n);
    if (nt <= 1) { for (int j = 0; j < n; ++j) fn(j); return; }
    std::atomic<int> next{0};
    std::mutex em; std::string err;
    auto work = [&]() {
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
        for (;;) {
            const int j = next++;
            if (j >= n) break;
            try { fn(j); } catch (const std::exception& e) { std::lock_guard<std::mutex> l(em); if (err.empty()) err = e.what(); }
        }
    };
    std::vector<std::thread> th;
    for (int t = 0; t < nt; ++t) th.emplace_back(work);
    for (auto& t : th) t.join();
    if (!err.empty()) throw std::runtime_error(err);
}
} // namespace

// ------------------------------------------------------------------------------------------------ recording
void Simulator::adjCapture(double t0, double dt) {
    AdjRec r;
    r.report = step_; r.t0 = t0; r.dt = dt;
    r.p = p_; r.sw = sw_; r.x = x_; r.hg = hasGas_;
    r.hfo = hfOpen_; r.hff = hfFrac_;   // C20: growth is tested after the capture, so this is the topology the step was solved with
    r.w.resize(wells_.size());
    for (size_t i = 0; i < wells_.size(); ++i) {
        const WellRT& W = wells_[i];
        AdjWell& a = r.w[i];
        a.bhp = W.bhp; a.mode = W.mode; a.type = W.type; a.open = W.open;
        a.ctrl = W.ctrl; a.own = W.own; a.dpConn = W.dpConn; a.headRho = W.headRho;
        for (int k = 0; k < 3; ++k) a.resvC[k] = W.resvC[k];
        a.resvOk = W.resvOk;
        a.oilH = W.rates.oilH; a.waterH = W.rates.waterH; a.gasH = W.rates.gasH;
    }
    for (const AquiferRT& Q : aquifers_) { r.aqPa.push_back(Q.pa); r.aqWe.push_back(Q.We); }
    r.rsLim = rsLim_; r.sgHy = sgHy_; r.soHy = soHy_; r.dG = dG_; r.dO = dO_;
    r.soMax = soMaxV_; r.fRv = fRv_; r.fRs = fRs_;
    if (dt > 0 && adjHeadIn_.size() == wells_.size()) {
        r.hbhp.resize(wells_.size()); r.dOld.resize(wells_.size());
        for (size_t i = 0; i < wells_.size(); ++i) { r.hbhp[i] = adjHeadIn_[i].first; r.dOld[i] = adjHeadIn_[i].second; }
    }
    adj_.push_back(std::move(r));
}

void Simulator::adjSetState(const AdjRec& r) {
    p_ = r.p; sw_ = r.sw; x_ = r.x; hasGas_ = r.hg;
    rsLim_ = r.rsLim;
    if (!r.sgHy.empty()) sgHy_ = r.sgHy;
    if (!r.soHy.empty()) soHy_ = r.soHy;
    if (!r.dG.empty()) dG_ = r.dG;
    if (!r.dO.empty()) dO_ = r.dO;
    if (!r.soMax.empty()) soMaxV_ = r.soMax;
    if (!r.fRv.empty()) fRv_ = r.fRv;
    if (!r.fRs.empty()) fRs_ = r.fRs;
}

// residual + Jacobian of time step 'cur' (old level 'prev') with this model's parameters at the recorded states
void Simulator::adjAssembleStep(const AdjRec& prev, const AdjRec& cur) {
    const int na = grid_.numActive();
    updateConnections(cur.report);
    // C20 HFGROW: the fracture transmissibilities of the recorded step (frozen topology). The open set is piecewise
    // constant in the parameters, so its derivative is zero except at opening-time switches: the gradient is exact on
    // the current branch, like control / phase switching events.
    if (!cur.hfo.empty() && cur.hfo.size() == hfOpen_.size() && (cur.hfo != hfOpen_ || cur.hff != hfFrac_)) {
        hfOpen_ = cur.hfo;
        if (cur.hff.size() == hfFrac_.size()) hfFrac_ = cur.hff;   // C21 sub-grid tips: the step's open fractions
        applyHfTrans();
    }
    // old accumulation (prev state with its own lagged data, as storeOld() computed it)
    adjSetState(prev);
    evalAll();
    for (int a = 0; a < na; ++a)
        for (int p = 0; p < 3; ++p) accOld_[size_t(a) * 3 + p] = (eq_[p] >= 0) ? cf_[size_t(a)].acc[p].v : 0.0;
    // own well controls of this model at the report step (rate-multiplier parameters scale the recorded targets)
    std::vector<WellControl> own(wells_.size());
    std::vector<char> hasOwn(wells_.size(), 0);
    for (int s = 0; s <= cur.report && s < int(case_.schedule.size()); ++s)
        for (const WellControl& wc : case_.schedule[size_t(s)].controls)
            for (size_t w = 0; w < wells_.size(); ++w)
                if (matchName(wc.well, case_.wells[size_t(wells_[w].input)].name)) { own[w] = wc; hasOwn[w] = 1; }
    for (size_t w = 0; w < wells_.size() && w < cur.w.size(); ++w) {
        WellRT& W = wells_[w];
        const AdjWell& R = cur.w[w];
        W.bhp = R.bhp; W.mode = R.mode; W.type = R.type; W.open = R.open;
        W.ctrl = R.ctrl; W.own = R.own;
        if (R.dpConn.size() == W.cell.size()) W.dpConn = R.dpConn;
        W.headRho = R.headRho;
        for (int k = 0; k < 3; ++k) W.resvC[k] = R.resvC[k];
        W.resvOk = R.resvOk;
        W.rates.oilH = R.oilH; W.rates.waterH = R.waterH; W.rates.gasH = R.gasH;
        if (hasOwn[w]) {
            auto scale = [](double& tgt, double baseOwn, double mine) { if (tgt > 0 && baseOwn > 0 && mine > 0) tgt *= mine / baseOwn; };
            scale(W.ctrl.orat, R.own.orat, own[w].orat); scale(W.ctrl.wrat, R.own.wrat, own[w].wrat);
            scale(W.ctrl.grat, R.own.grat, own[w].grat); scale(W.ctrl.lrat, R.own.lrat, own[w].lrat);
            scale(W.ctrl.resv, R.own.resv, own[w].resv); scale(W.ctrl.rate, R.own.rate, own[w].rate);
        }
    }
    for (size_t q = 0; q < aquifers_.size() && q < cur.aqPa.size(); ++q) { aquifers_[q].pa = cur.aqPa[q]; aquifers_[q].We = cur.aqWe[q]; }
    time_ = cur.t0;
    adjSetState(cur);
    evalAll();
    ensureScalarJacobian();   // C45: the adjoint reads the scalar Jacobian J_
    assemble(cur.dt);
}

// state and well inputs of the explicit well-head update at the start of step 'cur' (computeHeads() not called):
// cell state of the previous step, lagged data of this step, BHP and connection heads before the update
void Simulator::adjHeads(const AdjRec& prev, const AdjRec& cur) {
    updateConnections(cur.report);
    adjSetState(prev);
    rsLim_ = cur.rsLim;
    if (!cur.sgHy.empty()) sgHy_ = cur.sgHy;
    if (!cur.soHy.empty()) soHy_ = cur.soHy;
    if (!cur.dG.empty()) dG_ = cur.dG;
    if (!cur.dO.empty()) dO_ = cur.dO;
    if (!cur.soMax.empty()) soMaxV_ = cur.soMax;
    if (!cur.fRv.empty()) fRv_ = cur.fRv;
    if (!cur.fRs.empty()) fRs_ = cur.fRs;
    for (size_t w = 0; w < wells_.size() && w < cur.w.size(); ++w) {
        WellRT& W = wells_[w];
        W.open = cur.w[w].open; W.type = cur.w[w].type; W.mode = cur.w[w].mode;
        if (w < cur.hbhp.size()) { W.bhp = cur.hbhp[w]; W.dpConn = cur.dOld[w]; }
    }
}

SummaryPoint Simulator::adjPoint() const {
    SummaryPoint s;
    for (const auto& w : wells_) {
        WellRates r = w.rates; r.bhp = w.bhp; r.mode = w.mode; r.open = w.open;
        s.wells.push_back(r);
        s.fopr += r.oil; s.fwpr += r.water; s.fgpr += r.gas; s.fwir += r.waterInj; s.fgir += r.gasInj;
    }
    const double liq = s.fopr + s.fwpr;
    s.fwct = liq > 0 ? s.fwpr / liq : 0; s.fgor = s.fopr > 0 ? s.fgpr / s.fopr : 0;
    double ppv = 0, pv = 0, phc = 0, hc = 0;
    const int na = grid_.numActive();
    for (int a = 0; a < na; ++a) {
        ppv += p_[size_t(a)] * grid_.poreVolRef[size_t(a)]; pv += grid_.poreVolRef[size_t(a)];
        const double h = grid_.poreVolRef[size_t(a)] * (iSw_ >= 0 ? std::max(0.0, 1.0 - sw_[size_t(a)]) : 1.0);
        phc += p_[size_t(a)] * h; hc += h;
    }
    s.fprp = ppv / pv;
    s.fpr = (eq_[PO] >= 0 || eq_[PG] >= 0) && hc > 1e-12 * pv ? phc / hc : s.fprp;
    return s;
}

// ------------------------------------------------------------------------------------------------ gradient
Simulator::AdjResult Simulator::adjointGradient(const std::vector<AdjObs>& obs, const std::vector<AdjParam>& params, bool forwardMode,
                                                 const AdjFieldReq* field) {
    AdjResult res;
    const auto tStart = std::chrono::steady_clock::now();
    try {
        if (!adjRec_ || adj_.size() < 2) throw std::runtime_error("the forward run was not recorded for the adjoint (enable recording before initialize)");
        if (poly_) throw std::runtime_error("the adjoint does not support the polymer option (sequential transport)");
        if (summary_.size() != adj_.size()) throw std::runtime_error("summary and adjoint records are out of step (restarted run?)");
        const int N = int(adj_.size()) - 1;
        const int na = grid_.numActive(), nv = nv_, nw = int(wells_.size());
        const int nx = J_.n;
        const int no = int(obs.size()), np = int(params.size());
        res.steps = N;

        // ---- observation operators: y_i = sum_terms c * value(summary point k)
        struct Term { int obs; double c; };
        std::vector<std::vector<Term>> termsAt(size_t(N) + 1);
        std::vector<double> dFdy(size_t(no), 0.0);
        res.y.assign(size_t(no), 0.0);
        for (int i = 0; i < no; ++i) {
            const AdjObs& o = obs[size_t(i)];
            if (obsIsWellQ(o.qty) && (o.well < 0 || o.well >= nw)) throw std::runtime_error("observation of an unknown well");
            const auto& S = summary_;
            size_t k = 0;
            while (k < S.size() && S[k].time < o.time - 1e-6) ++k;
            std::vector<std::pair<int, double>> tk;
            const bool smooth = o.qty == ObsQty::FPR || o.qty == ObsQty::WBHP;
            if (k >= S.size()) tk.push_back({int(S.size()) - 1, 1.0});
            else if (smooth && k > 0 && S[k].time > S[k - 1].time) {
                const double a = (o.time - S[k - 1].time) / (S[k].time - S[k - 1].time);
                tk.push_back({int(k) - 1, 1 - a}); tk.push_back({int(k), a});
            } else if (k == 0 && S.size() > 1 && !smooth) tk.push_back({1, 1.0});
            else tk.push_back({int(k), 1.0});
            double y = 0;
            for (auto& t : tk) { y += t.second * pointVal(S[size_t(t.first)], o.qty, o.well); termsAt[size_t(t.first)].push_back({i, t.second}); }
            if (!std::isfinite(y)) throw std::runtime_error("observation value is not finite");
            res.y[size_t(i)] = y;
            const double r = (y - o.value) / o.sigma;
            res.F += r * r;
            dFdy[size_t(i)] = 2.0 * r / o.sigma;
        }

        // ---- replay models (fresh: setup + own initial equilibrium)
        auto make = [&](const Case& c) {
            std::unique_ptr<Simulator> s(new Simulator);
            s->adjRec_ = true;
            s->setup(c);
            if (s->hfOpen_.size() != hfOpen_.size())
                throw std::runtime_error("a perturbed model changed the hydraulic-fracture discretisation (HFGROW adjoint needs the same pieces)");
            s->initialize();
            if (s->grid_.numActive() != na || int(s->wells_.size()) != nw || s->nv_ != nv || s->J_.n != nx)
                throw std::runtime_error("a perturbed model changed the grid / well / phase structure (adjoint needs the same structure)");
            return s;
        };
        auto prevOf = [&](Simulator& S, int n) -> const AdjRec& { return n == 1 ? S.adj_[0] : adj_[size_t(n) - 1]; };

        // d y / d x at the current state of S (after adjAssembleStep with adjCollect_), scaled by coef, appended to out
        // (outD: derivative w.r.t. the connection heads d_n of the step, index coff[well] + k)
        std::vector<int> coff(size_t(nw) + 1, 0);
        for (int w = 0; w < nw; ++w) coff[size_t(w) + 1] = coff[size_t(w)] + int(wells_[size_t(w)].cell.size());
        const int nct = coff[size_t(nw)];
        auto stateGrad = [&](Simulator& S, ObsQty q, int w, double coef, Sparse& out, Sparse& outD, std::vector<double>* rc = nullptr) {
            auto addRate = [&](int well, int qi, double c) {
                if (c == 0) return;
                if (rc) { (*rc)[size_t(well) * 5 + size_t(qi)] += c; return; }   // dF/d(rate) only (field gradients)
                for (const AdjRateD& rd : S.adjRd_) {
                    if (rd.well != well) continue;
                    for (int v = 0; v < nv; ++v) if (rd.d[qi][v] != 0) out.push_back({rd.cell * nv + v, c * rd.d[qi][v]});
                    if (rd.d[qi][3] != 0) {   // rate = f(cell, p_conn = bhp + d_k): same derivative w.r.t. the BHP and the head
                        out.push_back({na * nv + well, c * rd.d[qi][3]});
                        outD.push_back({coff[size_t(well)] + rd.k, c * rd.d[qi][3]});
                    }
                }
            };
            auto ratesOf = [&](int well, double& o, double& wa, double& g) {
                const WellRates& r = S.wells_[size_t(well)].rates; o = r.oil; wa = r.water; g = r.gas;
            };
            switch (q) {
                case ObsQty::WOPR: addRate(w, 0, coef); break;
                case ObsQty::WWPR: addRate(w, 1, coef); break;
                case ObsQty::WGPR: addRate(w, 2, coef); break;
                case ObsQty::WWIR: addRate(w, 3, coef); break;
                case ObsQty::WBHP: if (!rc) out.push_back({na * nv + w, coef}); break;
                case ObsQty::WWCT: {
                    double o, wa, g; ratesOf(w, o, wa, g);
                    const double L = o + wa;
                    if (L > 0) { addRate(w, 1, coef * o / (L * L)); addRate(w, 0, -coef * wa / (L * L)); }
                } break;
                case ObsQty::WGOR: {
                    double o, wa, g; ratesOf(w, o, wa, g);
                    if (o > 0) { addRate(w, 2, coef / o); addRate(w, 0, -coef * g / (o * o)); }
                } break;
                case ObsQty::FOPR: for (int k = 0; k < nw; ++k) addRate(k, 0, coef); break;
                case ObsQty::FWPR: for (int k = 0; k < nw; ++k) addRate(k, 1, coef); break;
                case ObsQty::FGPR: for (int k = 0; k < nw; ++k) addRate(k, 2, coef); break;
                case ObsQty::FWCT: case ObsQty::FGOR: {
                    double O = 0, Wt = 0, G = 0;
                    for (int k = 0; k < nw; ++k) { double o, wa, g; ratesOf(k, o, wa, g); O += o; Wt += wa; G += g; }
                    if (q == ObsQty::FWCT) {
                        const double L = O + Wt;
                        if (L > 0) for (int k = 0; k < nw; ++k) { addRate(k, 1, coef * O / (L * L)); addRate(k, 0, -coef * Wt / (L * L)); }
                    } else if (O > 0) for (int k = 0; k < nw; ++k) { addRate(k, 2, coef / O); addRate(k, 0, -coef * G / (O * O)); }
                } break;
                case ObsQty::FPR: {
                    if (rc) break;
                    double pv = 0, phc = 0, hc = 0;
                    for (int a = 0; a < na; ++a) {
                        const double V = S.grid_.poreVolRef[size_t(a)];
                        const double h = V * (iSw_ >= 0 ? std::max(0.0, 1.0 - S.sw_[size_t(a)]) : 1.0);
                        pv += V; phc += S.p_[size_t(a)] * h; hc += h;
                    }
                    const bool hcw = (eq_[PO] >= 0 || eq_[PG] >= 0) && hc > 1e-12 * pv;
                    const double fpr = hcw ? phc / hc : 0;
                    for (int a = 0; a < na; ++a) {
                        const double V = S.grid_.poreVolRef[size_t(a)];
                        if (!hcw) { out.push_back({a * nv, coef * V / pv}); continue; }
                        const double so = iSw_ >= 0 ? 1.0 - S.sw_[size_t(a)] : 1.0;
                        if (so <= 0) continue;
                        out.push_back({a * nv, coef * V * so / hc});
                        if (iSw_ >= 0) out.push_back({a * nv + iSw_, -coef * V * (S.p_[size_t(a)] - fpr) / hc});
                    }
                } break;
                default: break;
            }
        };

        std::mutex resM;
        // linear solve with ILU(0)-GMRES, tight tolerance; returns the true relative residual
        auto solve = [&](const CSR& A, ILU0& M, const std::vector<double>& b, std::vector<double>& x) {
            x.assign(b.size(), 0.0);
            const double bn = norm2(b);
            if (bn == 0) return 0.0;
            LinResult lr = gmres(A, M, b, x, 1e-12, 1500, 100);
            int its = lr.iterations;
            std::vector<double> Ax(b.size());
            A.spmv(x.data(), Ax.data());
            double rn = 0; for (size_t i = 0; i < b.size(); ++i) rn += (Ax[i] - b[i]) * (Ax[i] - b[i]);
            double rel = std::sqrt(rn) / bn;
            if (rel > 1e-8) {   // restart from the current iterate with BiCGSTAB as a second chance
                std::vector<double> r(b.size()), dx;
                for (size_t i = 0; i < b.size(); ++i) r[i] = b[i] - Ax[i];
                dx.assign(b.size(), 0.0);
                LinResult l2 = bicgstab(A, M, r, dx, 1e-10, 3000);
                its += l2.iterations;
                for (size_t i = 0; i < b.size(); ++i) x[i] += dx[i];
                A.spmv(x.data(), Ax.data());
                rn = 0; for (size_t i = 0; i < b.size(); ++i) rn += (Ax[i] - b[i]) * (Ax[i] - b[i]);
                rel = std::sqrt(rn) / bn;
            }
            { std::lock_guard<std::mutex> l(resM); res.maxLinRes = std::max(res.maxLinRes, rel); res.linearIts += its; }
            return rel;
        };

        std::unique_ptr<Simulator> B = make(case_);
        B->adjCollect_ = true;
        auto replayCheck = [&](int n) {   // the replay must reproduce the forward run (same rates to round-off)
            const SummaryPoint sp = B->adjPoint();
            const SummaryPoint& f = summary_[size_t(n)];
            auto dev = [](double a, double b) { return std::fabs(a - b) / std::max({std::fabs(a), std::fabs(b), 1e-3}); };
            for (int w = 0; w < nw && w < int(f.wells.size()); ++w) {
                res.replayDev = std::max(res.replayDev, dev(sp.wells[size_t(w)].oil, f.wells[size_t(w)].oil));
                res.replayDev = std::max(res.replayDev, dev(sp.wells[size_t(w)].water, f.wells[size_t(w)].water));
                res.replayDev = std::max(res.replayDev, dev(sp.wells[size_t(w)].gas, f.wells[size_t(w)].gas));
                res.replayDev = std::max(res.replayDev, dev(sp.wells[size_t(w)].waterInj, f.wells[size_t(w)].waterInj));
                res.replayDev = std::max(res.replayDev, dev(sp.wells[size_t(w)].gasInj, f.wells[size_t(w)].gasInj));
            }
        };

        // ---- explicit well heads d_n = H(x_{n-1}, bhp, d_{n-1}; theta) (computeHeads at the start of each step):
        // augmented lagged variables with their own multipliers mu_n
        auto headsOf = [&](Simulator& S, std::vector<double>& d) {
            d.assign(size_t(nct), 0.0);
            for (int w = 0; w < nw; ++w) {
                const auto& dp = S.wells_[size_t(w)].dpConn;
                for (size_t k = 0; k < dp.size() && int(k) < coff[size_t(w) + 1] - coff[size_t(w)]; ++k) d[size_t(coff[size_t(w)]) + k] = dp[k];
            }
        };
        // heads of step n at the recorded inputs with the parameters of S
        auto headsAt = [&](Simulator& S, int n, std::vector<double>& d) {
            S.adjHeads(prevOf(S, n), adj_[size_t(n)]);
            S.headsOnly_ = -1;
            S.computeHeads();
            headsOf(S, d);
        };
        struct HeadJac { int w = 0, nk = 0; std::vector<int> inX, inD; std::vector<double> Jx, Jd; };   // Jx: nk x inX, Jd: nk x inD
        // Jacobian of the head update of step n w.r.t. x_{n-1} (well cells, BHP when it is the previous step's unknown)
        // and d_{n-1} (when the old heads are the previous step's heads), central differences per well
        auto headJac = [&](Simulator& S, int n) {
            std::vector<HeadJac> out;
            const AdjRec& cur = adj_[size_t(n)];
            if (cur.hbhp.empty() || n < 2) return out;   // step 1: x_0 / d_0 depend on theta only (perturbed models)
            const AdjRec& prev = adj_[size_t(n) - 1];
            S.adjHeads(prev, cur);
            for (int w = 0; w < nw; ++w) {
                WellRT& W = S.wells_[size_t(w)];
                if (!W.open || W.cell.empty()) continue;
                HeadJac H; H.w = w; H.nk = int(W.cell.size());
                std::vector<int> cells = W.cell;
                std::sort(cells.begin(), cells.end()); cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
                for (int c : cells) for (int v = 0; v < nv; ++v) H.inX.push_back(c * nv + v);
                const bool chainB = cur.hbhp[size_t(w)] == prev.w[size_t(w)].bhp;
                if (chainB) H.inX.push_back(na * nv + w);
                const bool chainD = !cur.dOld[size_t(w)].empty() && cur.dOld[size_t(w)] == prev.w[size_t(w)].dpConn;
                if (chainD) for (int k = 0; k < H.nk; ++k) H.inD.push_back(k);
                const std::vector<double> d0 = cur.dOld[size_t(w)];
                const double b0 = cur.hbhp[size_t(w)];
                auto evalW = [&]() { W.bhp = b0; W.dpConn = d0; S.headsOnly_ = w; S.computeHeads(); S.headsOnly_ = -1; std::vector<double> r = W.dpConn; r.resize(size_t(H.nk), 0.0); return r; };
                H.Jx.assign(size_t(H.nk) * H.inX.size(), 0.0);
                H.Jd.assign(size_t(H.nk) * H.inD.size(), 0.0);
                for (size_t q = 0; q < H.inX.size(); ++q) {
                    const int gi = H.inX[q];
                    double* var; bool pres;
                    double bb = b0;
                    if (gi >= na * nv) { var = &bb; pres = true; }
                    else {
                        const int c = gi / nv, v = gi % nv;
                        pres = v == 0;
                        var = v == 0 ? &S.p_[size_t(c)] : (v == iSw_ ? &S.sw_[size_t(c)] : &S.x_[size_t(c)]);
                    }
                    const double v0 = *var;
                    const double h = pres ? 1e-7 * std::max(std::fabs(v0), 1.0) : 1e-6 * std::max(std::fabs(v0), 1e-2);
                    const bool central = pres || v0 - h >= 0;
                    auto evalAt = [&](double val) {
                        if (gi >= na * nv) { W.dpConn = d0; W.bhp = val; S.headsOnly_ = w; S.computeHeads(); S.headsOnly_ = -1; std::vector<double> r = W.dpConn; r.resize(size_t(H.nk), 0.0); return r; }
                        *var = val; std::vector<double> r = evalW(); *var = v0; return r;
                    };
                    const std::vector<double> a = evalAt(v0 + h), b = central ? evalAt(v0 - h) : evalAt(v0);
                    const double den = central ? 2 * h : h;
                    for (int k = 0; k < H.nk; ++k) H.Jx[size_t(k) * H.inX.size() + q] = (a[size_t(k)] - b[size_t(k)]) / den;
                }
                for (size_t q = 0; q < H.inD.size(); ++q) {
                    const int kk = H.inD[q];
                    const double v0 = d0[size_t(kk)], h = 1e-7 * std::max(std::fabs(v0), 1.0);
                    auto evalAt = [&](double val) { std::vector<double> dd = d0; dd[size_t(kk)] = val; W.bhp = b0; W.dpConn = dd; S.headsOnly_ = w; S.computeHeads(); S.headsOnly_ = -1; std::vector<double> r = W.dpConn; r.resize(size_t(H.nk), 0.0); return r; };
                    const std::vector<double> a = evalAt(v0 + h), b = evalAt(v0 - h);
                    for (int k = 0; k < H.nk; ++k) H.Jd[size_t(k) * H.inD.size() + q] = (a[size_t(k)] - b[size_t(k)]) / (2 * h);
                }
                W.bhp = b0; W.dpConn = d0;
                out.push_back(std::move(H));
            }
            return out;
        };
        // R_n,d : derivative of the step-n residual (B's last assembly) w.r.t. the heads, applied / transposed
        auto RdT = [&](const std::vector<double>& lam, double dt, std::vector<double>& out) {   // out = R_d^T lam
            out.assign(size_t(nct), 0.0);
            for (const AdjConnD& cd : B->adjCd_) {
                const double sc = dt / B->grid_.poreVolRef[size_t(cd.cell)];
                double v = 0;
                for (int e = 0; e < nv; ++e) v += cd.dc[e] * sc * lam[size_t(cd.cell) * nv + e];
                v += cd.dw * lam[size_t(na) * nv + cd.well];
                out[size_t(coff[size_t(cd.well)] + cd.k)] += v;
            }
        };
        auto RdApply = [&](const std::vector<double>& sd, double dt, std::vector<double>& rhs, double sign) {   // rhs += sign R_d sd
            for (const AdjConnD& cd : B->adjCd_) {
                const double x = sd[size_t(coff[size_t(cd.well)] + cd.k)];
                if (x == 0) continue;
                const double sc = dt / B->grid_.poreVolRef[size_t(cd.cell)];
                for (int e = 0; e < nv; ++e) rhs[size_t(cd.cell) * nv + e] += sign * cd.dc[e] * sc * x;
                rhs[size_t(na) * nv + cd.well] += sign * cd.dw * x;
            }
        };
        // head check: the replayed head update must reproduce the recorded heads
        auto headCheck = [&](int n) {
            if (adj_[size_t(n)].hbhp.empty()) return;
            std::vector<double> d; headsAt(*B, n, d);
            for (int w = 0; w < nw; ++w) {
                const auto& rec = adj_[size_t(n)].w[size_t(w)].dpConn;
                for (size_t k = 0; k < rec.size(); ++k)
                    res.replayDev = std::max(res.replayDev, std::fabs(d[size_t(coff[size_t(w)]) + k] - rec[k]) / std::max(std::fabs(rec[k]), 1.0));
            }
        };

        // ---- field parameters (analytic dR/d ln theta, accumulated during the reverse pass)
        const bool fT = field && (field->perm || field->trans), fW = field && (field->perm || field->wi), fP = field && field->poro;
        if (field && forwardMode) throw std::runtime_error("field gradients need the adjoint (reverse) mode");
        if (field && (field->perm || field->trans) && !grid_.fieldDerivOk)
            throw std::runtime_error("permeability / transmissibility field gradients are not available with local grid refinements or dual porosity");
        const size_t nconn = grid_.conns.size();
        std::vector<double> gT(fT ? nconn : 0, 0.0), gP(fP ? size_t(na) : 0, 0.0);
        std::vector<std::vector<double>> gW;
        if (fW) for (int w = 0; w < nw; ++w) gW.emplace_back(wells_[size_t(w)].cell.size(), 0.0);
        double secField = 0;
        // explicit FPR derivative w.r.t. ln(pore volume) of every cell at state (p, sw)
        auto fprPoro = [&](const std::vector<double>& P, const std::vector<double>& SW, double coef) {
            double pv = 0, ppv = 0, phc = 0, hc = 0;
            for (int a = 0; a < na; ++a) {
                const double V = grid_.poreVolRef[size_t(a)];
                const double h = V * (iSw_ >= 0 ? std::max(0.0, 1.0 - SW[size_t(a)]) : 1.0);
                pv += V; ppv += P[size_t(a)] * V; phc += P[size_t(a)] * h; hc += h;
            }
            const bool hcw = (eq_[PO] >= 0 || eq_[PG] >= 0) && hc > 1e-12 * pv;
            const double fpr = hcw ? phc / hc : ppv / pv;
            for (int a = 0; a < na; ++a) {
                const double V = grid_.poreVolRef[size_t(a)];
                const double h = hcw ? V * (iSw_ >= 0 ? std::max(0.0, 1.0 - SW[size_t(a)]) : 1.0) : V;
                gP[size_t(a)] += coef * h * (P[size_t(a)] - fpr) / (hcw ? hc : pv);
            }
        };
        if (fP)
            for (const Term& t : termsAt[0])
                if (obs[size_t(t.obs)].qty == ObsQty::FPR) fprPoro(adj_[0].p, adj_[0].sw, dFdy[size_t(t.obs)] * t.c);

        res.grad.assign(size_t(np), 0.0);
        if (!forwardMode) {
            // ================= reverse (adjoint) pass
            std::vector<std::vector<double>> lam(size_t(N) + 1), mu(size_t(N) + 1);
            std::vector<double> rhs(static_cast<size_t>(nx)), lamNext, muNext;
            CSR T; std::vector<int> tmap; ILU0 ilu;
            bool patt = false;
            for (int n = N; n >= 1; --n) {
                const bool needLam = !lamNext.empty() || !termsAt[size_t(n)].empty();
                if (!needLam) continue;
                std::vector<HeadJac> HJ;
                if (!muNext.empty() && n < N) HJ = headJac(*B, n + 1);
                if (n % 7 == 0 || n == N) headCheck(n);
                B->adjAssembleStep(prevOf(*B, n), adj_[size_t(n)]);
                replayCheck(n);
                std::fill(rhs.begin(), rhs.end(), 0.0);
                Sparse g, gD;   // -dF/dx_n, -dF/dd_n
                for (const Term& t : termsAt[size_t(n)]) stateGrad(*B, obs[size_t(t.obs)].qty, obs[size_t(t.obs)].well, -dFdy[size_t(t.obs)] * t.c, g, gD);
                for (auto& e : g) rhs[size_t(e.first)] += e.second;
                if (!lamNext.empty())
                    for (int a = 0; a < na; ++a) {
                        const double V = B->grid_.poreVolRef[size_t(a)];
                        for (int p = 0; p < 3; ++p) {
                            if (eq_[p] < 0) continue;
                            const double l = lamNext[size_t(a) * nv + eq_[p]] / V;
                            if (l == 0) continue;
                            const AD3& acc = B->cf_[size_t(a)].acc[p];
                            for (int k = 0; k < nv; ++k) rhs[size_t(a) * nv + k] += acc.d[k] * l;
                        }
                    }
                for (const HeadJac& H : HJ)   // + H_{n+1,x}^T mu_{n+1}
                    for (int k = 0; k < H.nk; ++k) {
                        const double m = muNext[size_t(coff[size_t(H.w)] + k)];
                        if (m == 0) continue;
                        for (size_t q = 0; q < H.inX.size(); ++q) rhs[size_t(H.inX[q])] += H.Jx[size_t(k) * H.inX.size() + q] * m;
                    }
                if (!patt) { transposePattern(B->J_, T, tmap); patt = true; }
                for (size_t q = 0; q < tmap.size(); ++q) T.val[size_t(tmap[q])] = B->J_.val[q];
                ilu.factor(T);
                std::vector<double> x;
                solve(T, ilu, rhs, x);
                // mu_n = H_{n+1,d}^T mu_{n+1} - R_{n,d}^T lam_n - dF/dd_n
                std::vector<double> m;
                RdT(x, adj_[size_t(n)].dt, m);
                for (double& v : m) v = -v;
                for (auto& e : gD) m[size_t(e.first)] += e.second;
                for (const HeadJac& H : HJ)
                    for (int k = 0; k < H.nk; ++k) {
                        const double mk = muNext[size_t(coff[size_t(H.w)] + k)];
                        if (mk == 0) continue;
                        for (size_t q = 0; q < H.inD.size(); ++q) m[size_t(coff[size_t(H.w)] + H.inD[q])] += H.Jd[size_t(k) * H.inD.size() + q] * mk;
                    }
                if (field) {   // lam_n^T dR_n/d ln theta for every cell / face / connection, + explicit dF/d theta|_x
                    const auto tf = std::chrono::steady_clock::now();
                    const double dtn = adj_[size_t(n)].dt;
                    const std::vector<double>& L = x;
                    const auto& PV = B->grid_.poreVolRef;
                    if (fT)
                        for (size_t f = 0; f < nconn; ++f) {
                            const Connection& cn = B->grid_.conns[f];
                            const double sa = dtn / PV[size_t(cn.a)], sb = dtn / PV[size_t(cn.b)];
                            const double* F = &B->adjFace_[f * size_t(nv)];
                            double v = 0;
                            for (int e = 0; e < nv; ++e) v += F[e] * (L[size_t(cn.a) * nv + e] * sa - L[size_t(cn.b) * nv + e] * sb);
                            gT[f] += v;
                        }
                    if (fW) {
                        for (const AdjConnD& cd : B->adjCd_) {
                            const double sc = dtn / PV[size_t(cd.cell)];
                            double v = cd.vw * L[size_t(na) * nv + cd.well];
                            for (int e = 0; e < nv; ++e) v += cd.vc[e] * sc * L[size_t(cd.cell) * nv + e];
                            gW[size_t(cd.well)][size_t(cd.k)] += v;
                        }
                        if (!termsAt[size_t(n)].empty()) {   // rates are linear in WI: d y / d ln WI_k = connection k's rate terms
                            std::vector<double> rc(size_t(nw) * 5, 0.0);
                            Sparse d1, d2;
                            for (const Term& t : termsAt[size_t(n)]) stateGrad(*B, obs[size_t(t.obs)].qty, obs[size_t(t.obs)].well, dFdy[size_t(t.obs)] * t.c, d1, d2, &rc);
                            for (const AdjRateD& rd : B->adjRd_)
                                for (int qi = 0; qi < 5; ++qi) gW[size_t(rd.well)][size_t(rd.k)] += rc[size_t(rd.well) * 5 + size_t(qi)] * rd.v[qi];
                        }
                        if (!adj_[size_t(n)].hbhp.empty()) {   // - mu_n^T dH_n/d ln WI_k (head update of the step, per well)
                            const AdjRec& cur = adj_[size_t(n)];
                            B->adjHeads(prevOf(*B, n), cur);
                            const double he = 1e-4;
                            for (int w = 0; w < nw; ++w) {
                                WellRT& W = B->wells_[size_t(w)];
                                if (!W.open || W.cell.empty() || size_t(w) >= cur.hbhp.size()) continue;
                                bool any = false;
                                for (int k = coff[size_t(w)]; k < coff[size_t(w) + 1]; ++k) if (m[size_t(k)] != 0) any = true;
                                if (!any) continue;
                                for (size_t k = 0; k < W.cell.size(); ++k) {
                                    const double wi0 = W.wi[k];
                                    if (!(wi0 > 0)) continue;
                                    auto ev = [&](double f) {
                                        W.wi[k] = wi0 * f; W.bhp = cur.hbhp[size_t(w)]; W.dpConn = cur.dOld[size_t(w)];
                                        B->headsOnly_ = w; B->computeHeads(); B->headsOnly_ = -1;
                                        std::vector<double> r = W.dpConn; r.resize(W.cell.size(), 0.0); return r;
                                    };
                                    const std::vector<double> hp = ev(std::exp(he)), hm = ev(std::exp(-he));
                                    W.wi[k] = wi0;
                                    double sum = 0;
                                    for (size_t j = 0; j < W.cell.size(); ++j) sum += m[size_t(coff[size_t(w)]) + j] * (hp[j] - hm[j]);
                                    gW[size_t(w)][k] -= sum / (2 * he);
                                }
                                W.bhp = cur.hbhp[size_t(w)]; W.dpConn = cur.dOld[size_t(w)];
                            }
                        }
                    }
                    if (fP) {
                        for (int a = 0; a < na; ++a)
                            for (int p = 0; p < 3; ++p) {
                                if (eq_[p] < 0) continue;
                                const int r = a * nv + eq_[p];
                                const double dacc = B->cf_[size_t(a)].acc[p].v - B->accOld_[size_t(a) * 3 + p];
                                gP[size_t(a)] += L[size_t(r)] * (dacc / PV[size_t(a)] - B->R_[size_t(r)]);
                            }
                        for (const Term& t : termsAt[size_t(n)])
                            if (obs[size_t(t.obs)].qty == ObsQty::FPR) fprPoro(adj_[size_t(n)].p, adj_[size_t(n)].sw, dFdy[size_t(t.obs)] * t.c);
                    }
                    secField += std::chrono::duration<double>(std::chrono::steady_clock::now() - tf).count();
                }
                lam[size_t(n)] = x; mu[size_t(n)] = m;
                lamNext = std::move(x); muNext = std::move(m);
            }
            // ================= parameter contractions
            parallelFor(np, [&](int j) {
                const AdjParam& P = params[size_t(j)];
                if (!P.plus || !P.minus || P.h <= 0) throw std::runtime_error("adjoint parameter without perturbed models");
                std::unique_ptr<Simulator> Sp = make(*P.plus), Sm = make(*P.minus);
                const double inv2h = 0.5 / P.h;
                double g = 0;
                for (const Term& t : termsAt[0]) {
                    const AdjObs& o = obs[size_t(t.obs)];
                    g += dFdy[size_t(t.obs)] * t.c * (pointVal(Sp->summary_[0], o.qty, o.well) - pointVal(Sm->summary_[0], o.qty, o.well)) * inv2h;
                }
                std::vector<double> dp, dm;
                for (int n = 1; n <= N; ++n) {
                    if (lam[size_t(n)].empty() && termsAt[size_t(n)].empty()) continue;
                    if (!mu[size_t(n)].empty() && !adj_[size_t(n)].hbhp.empty()) {   // - mu_n^T H_{n,theta}
                        headsAt(*Sp, n, dp); headsAt(*Sm, n, dm);
                        double s = 0;
                        for (int q = 0; q < nct; ++q) s += mu[size_t(n)][size_t(q)] * (dp[size_t(q)] - dm[size_t(q)]);
                        g -= s * inv2h;
                    }
                    Sp->adjAssembleStep(prevOf(*Sp, n), adj_[size_t(n)]);
                    Sm->adjAssembleStep(prevOf(*Sm, n), adj_[size_t(n)]);
                    if (!lam[size_t(n)].empty()) {
                        const std::vector<double>& L = lam[size_t(n)];
                        double s = 0;
                        for (int r = 0; r < nx; ++r) s += L[size_t(r)] * (Sp->R_[size_t(r)] - Sm->R_[size_t(r)]);
                        g += s * inv2h;
                    }
                    if (!termsAt[size_t(n)].empty()) {
                        const SummaryPoint a = Sp->adjPoint(), b = Sm->adjPoint();
                        for (const Term& t : termsAt[size_t(n)]) {
                            const AdjObs& o = obs[size_t(t.obs)];
                            g += dFdy[size_t(t.obs)] * t.c * (pointVal(a, o.qty, o.well) - pointVal(b, o.qty, o.well)) * inv2h;
                        }
                    }
                }
                res.grad[size_t(j)] = g;
            });
            if (field) {   // map to the field outputs
                if (fT) {
                    res.gLogTrans = gT;
                    if (field->perm) {
                        res.gLogKx.assign(size_t(na), 0.0); res.gLogKy.assign(size_t(na), 0.0); res.gLogKz.assign(size_t(na), 0.0);
                        for (size_t f = 0; f < nconn; ++f) {
                            const Connection& cn = grid_.conns[f];
                            std::vector<double>* G = cn.dir == 'X' ? &res.gLogKx : cn.dir == 'Y' ? &res.gLogKy : cn.dir == 'Z' ? &res.gLogKz : nullptr;
                            if (!G || !(cn.trans > 0)) continue;
                            (*G)[size_t(cn.a)] += gT[f] * cn.dTa / cn.trans;
                            (*G)[size_t(cn.b)] += gT[f] * cn.dTb / cn.trans;
                        }
                    }
                }
                if (fW) {
                    res.gLogWI = gW;
                    if (field->perm)
                        for (int w = 0; w < nw; ++w) {
                            const WellRT& W = wells_[size_t(w)];
                            for (size_t k = 0; k < W.cell.size() && k < W.dlwi.size(); ++k) {
                                const size_t c = size_t(W.cell[k]);
                                res.gLogKx[c] += gW[size_t(w)][k] * W.dlwi[k][0];
                                res.gLogKy[c] += gW[size_t(w)][k] * W.dlwi[k][1];
                                res.gLogKz[c] += gW[size_t(w)][k] * W.dlwi[k][2];
                            }
                        }
                }
                if (fP) res.gLogPoro = gP;
                res.secField = secField;
            }
        } else {
            // ================= forward (tangent) sensitivities, all parameters together
            std::vector<std::unique_ptr<Simulator>> Sp(static_cast<size_t>(np)), Sm(static_cast<size_t>(np));
            parallelFor(np, [&](int j) {
                const AdjParam& P = params[size_t(j)];
                if (!P.plus || !P.minus || P.h <= 0) throw std::runtime_error("sensitivity parameter without perturbed models");
                Sp[size_t(j)] = make(*P.plus); Sm[size_t(j)] = make(*P.minus);
            });
            res.dydp.assign(size_t(no), std::vector<double>(size_t(np), 0.0));
            for (const Term& t : termsAt[0])
                for (int j = 0; j < np; ++j) {
                    const AdjObs& o = obs[size_t(t.obs)];
                    res.dydp[size_t(t.obs)][size_t(j)] += t.c * (pointVal(Sp[size_t(j)]->summary_[0], o.qty, o.well) -
                                                                 pointVal(Sm[size_t(j)]->summary_[0], o.qty, o.well)) * 0.5 / params[size_t(j)].h;
                }
            std::vector<std::vector<double>> s(size_t(np), std::vector<double>(size_t(nx), 0.0));
            std::vector<std::vector<double>> sd(size_t(np), std::vector<double>(size_t(nct), 0.0));
            std::vector<double> accPrev;   // d acc / dx at x_{n-1}: na * 3 * nv
            ILU0 ilu;
            for (int n = 1; n <= N; ++n) {
                const bool heads = !adj_[size_t(n)].hbhp.empty();
                std::vector<HeadJac> HJ;
                if (heads) HJ = headJac(*B, n);
                if (heads) {   // head sensitivities sd_n = H_theta + H_x s_{n-1} + H_d sd_{n-1}
                    parallelFor(np, [&](int j) {
                        std::vector<double> dp, dm;
                        headsAt(*Sp[size_t(j)], n, dp); headsAt(*Sm[size_t(j)], n, dm);
                        std::vector<double> nd(static_cast<size_t>(nct));
                        for (int q = 0; q < nct; ++q) nd[size_t(q)] = (dp[size_t(q)] - dm[size_t(q)]) * 0.5 / params[size_t(j)].h;
                        for (const HeadJac& H : HJ)
                            for (int k = 0; k < H.nk; ++k) {
                                double v = 0;
                                for (size_t q = 0; q < H.inX.size(); ++q) v += H.Jx[size_t(k) * H.inX.size() + q] * s[size_t(j)][size_t(H.inX[q])];
                                for (size_t q = 0; q < H.inD.size(); ++q) v += H.Jd[size_t(k) * H.inD.size() + q] * sd[size_t(j)][size_t(coff[size_t(H.w)] + H.inD[q])];
                                nd[size_t(coff[size_t(H.w)] + k)] += v;
                            }
                        sd[size_t(j)] = std::move(nd);
                    });
                } else for (auto& v : sd) std::fill(v.begin(), v.end(), 0.0);
                if (n % 7 == 0 || n == N) headCheck(n);
                B->adjAssembleStep(prevOf(*B, n), adj_[size_t(n)]);
                replayCheck(n);
                ilu.factor(B->J_);
                std::vector<Sparse> tg(termsAt[size_t(n)].size()), tgD(termsAt[size_t(n)].size());
                for (size_t ti = 0; ti < termsAt[size_t(n)].size(); ++ti) {
                    const Term& t = termsAt[size_t(n)][ti];
                    stateGrad(*B, obs[size_t(t.obs)].qty, obs[size_t(t.obs)].well, 1.0, tg[ti], tgD[ti]);
                }
                parallelFor(np, [&](int j) {
                    std::vector<double> rhs(static_cast<size_t>(nx));
                    Simulator& A = *Sp[size_t(j)]; Simulator& C = *Sm[size_t(j)];
                    A.adjAssembleStep(prevOf(A, n), adj_[size_t(n)]);
                    C.adjAssembleStep(prevOf(C, n), adj_[size_t(n)]);
                    const double inv2h = 0.5 / params[size_t(j)].h;
                    for (int r = 0; r < nx; ++r) rhs[size_t(r)] = -(A.R_[size_t(r)] - C.R_[size_t(r)]) * inv2h;
                    if (heads) RdApply(sd[size_t(j)], adj_[size_t(n)].dt, rhs, -1.0);
                    if (!accPrev.empty()) {
                        const std::vector<double>& sp = s[size_t(j)];
                        for (int a = 0; a < na; ++a) {
                            const double V = B->grid_.poreVolRef[size_t(a)];
                            for (int p = 0; p < 3; ++p) {
                                if (eq_[p] < 0) continue;
                                double d = 0;
                                for (int k = 0; k < nv; ++k) d += accPrev[(size_t(a) * 3 + p) * nv + k] * sp[size_t(a) * nv + k];
                                rhs[size_t(a) * nv + eq_[p]] += d / V;
                            }
                        }
                    }
                    std::vector<double> x;
                    solve(B->J_, ilu, rhs, x);
                    s[size_t(j)] = std::move(x);
                    if (!termsAt[size_t(n)].empty()) {
                        const SummaryPoint pa = A.adjPoint(), pc = C.adjPoint();
                        for (size_t ti = 0; ti < termsAt[size_t(n)].size(); ++ti) {
                            const Term& t = termsAt[size_t(n)][ti];
                            const AdjObs& o = obs[size_t(t.obs)];
                            double d = (pointVal(pa, o.qty, o.well) - pointVal(pc, o.qty, o.well)) * inv2h;
                            for (auto& e : tg[ti]) d += e.second * s[size_t(j)][size_t(e.first)];
                            for (auto& e : tgD[ti]) d += e.second * sd[size_t(j)][size_t(e.first)];
                            res.dydp[size_t(t.obs)][size_t(j)] += t.c * d;
                        }
                    }
                });
                accPrev.assign(size_t(na) * 3 * nv, 0.0);
                for (int a = 0; a < na; ++a)
                    for (int p = 0; p < 3; ++p) if (eq_[p] >= 0)
                        for (int k = 0; k < nv; ++k) accPrev[(size_t(a) * 3 + p) * nv + k] = B->cf_[size_t(a)].acc[p].d[k];
            }
            for (int j = 0; j < np; ++j) {
                double g = 0;
                for (int i = 0; i < no; ++i) g += dFdy[size_t(i)] * res.dydp[size_t(i)][size_t(j)];
                res.grad[size_t(j)] = g;
            }
        }
        res.ok = true;
    } catch (const std::exception& e) { res.ok = false; res.error = e.what(); }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
    return res;
}

} // namespace ressim
