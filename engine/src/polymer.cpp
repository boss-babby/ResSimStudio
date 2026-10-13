// Polymer flooding (ECLIPSE / OPM Flow POLYMER option) for the black-oil simulator.
//
// Model (identical to OPM Flow's blackoilpolymermodules):
//   storage   M = PV(p) * Sw * bw * c * (1 - DPV) + (V - PV(p)) * rho_r * Ca(c)          [mass]
//   flux      F_c = F_w * c_up * (mu_w,eff / mu_p,eff)(c_up)                               (water upwind)
//   wells     injection: q_w * c_inj; production: q_w * c * (mu_w,eff / mu_p,eff)
//   water mobility  lambda_w <- lambda_w * (mu_w / mu_w,eff) / R_k,  R_k = 1 + (RRF - 1) * Ca / Ca,max
//   Todd-Longstaff: mu_m = M(c) mu_w, mu_p = M(cmax) mu_w, mu_p,eff = mu_m^w mu_p^(1-w), mu_w,eff' = mu_m^w mu_w^(1-w)
//                   1 / mu_w,eff = (1 - c/cmax) / mu_w,eff' + (c/cmax) / mu_p,eff
//   irreversible adsorption (PLYROCK item 4 = 2): Ca = max(historical maximum, Ca(c))
//
// Solution: sequential-implicit. The flow equations are solved with the water mobility multiplier of the current
// concentration iterate; the polymer equation is then solved implicitly (Newton, exact Jacobian, ILU(0)-BiCGSTAB)
// with the converged fluxes; flow and transport are iterated inside the time step until the multiplier no longer
// changes. At convergence the result equals the fully coupled solution.
#include "ressim/simulator.h"
#include "ressim/ad.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <stdexcept>

namespace ressim {

namespace {
// piecewise-linear table (x, y rows), linear extrapolation beyond the ends (OPM: extrapolate = true)
double tab(const std::vector<double>& t, double x, double* dydx) {
    const size_t n = t.size() / 2;
    if (n == 0) { if (dydx) *dydx = 0; return 0; }
    if (n == 1) { if (dydx) *dydx = 0; return t[1]; }
    size_t i = 0;
    if (x <= t[0]) i = 0;
    else if (x >= t[2 * (n - 1)]) i = n - 2;
    else { while (i + 2 < n && x > t[2 * (i + 1)]) ++i; }
    const double x0 = t[2 * i], x1 = t[2 * i + 2], y0 = t[2 * i + 1], y1 = t[2 * i + 3];
    const double s = x1 != x0 ? (y1 - y0) / (x1 - x0) : 0.0;
    if (dydx) *dydx = s;
    return y0 + s * (x - x0);
}
std::string f6(double v) { char b[64]; std::snprintf(b, sizeof b, "%.6g", v); return b; }
}

Simulator::PolyEval Simulator::polyEval(int a, double c) const {
    const PolymerInput& P = case_.polymer;
    PolyEval e{};
    const auto& vis = P.plyvisc[std::min(size_t(pvtReg_[size_t(a)]), P.plyvisc.size() - 1)];
    const size_t sr = std::min(size_t(satReg_[size_t(a)]), P.plyrock.size() - 1);
    const PolymerRock& R = P.plyrock[sr];
    const auto& ads = P.plyads[std::min(size_t(satReg_[size_t(a)]), P.plyads.size() - 1)];
    // viscosity corrections (AD in c)
    using A1 = AD<1>;
    const A1 C = A1::variable(c, 0);
    double dm = 0; const double m0 = tab(vis, c, &dm);
    const A1 m = chain(C, std::max(m0, 1e-12), m0 > 1e-12 ? dm : 0.0);
    const size_t mr = mixReg_.empty() ? 0 : size_t(mixReg_[size_t(a)]);
    const double w = P.mixing[std::min(mr, P.mixing.size() - 1)];
    const double cmax = P.cmax[std::min(mr, P.cmax.size() - 1)];
    const double mp = std::max(1e-12, tab(vis, cmax, nullptr));
    const A1 mw = pow(m, w);
    const double mp1 = std::pow(mp, 1.0 - w);
    const A1 cbar = cmax > 0 ? C / cmax : A1(0.0);
    const A1 wvc = (1.0 - cbar) / mw + cbar / (mw * mp1);
    const A1 pvc = 1.0 / (wvc * mw * mp1);
    e.wvc = wvc.v; e.dwvc = wvc.d[0]; e.pvc = pvc.v; e.dpvc = pvc.d[0];
    // adsorption
    double da = 0; double ca = tab(ads, std::max(c, 0.0), &da);
    if (c < 0) da = 0;
    if (R.adsIndex == 2 && !pcAdsMax_.empty() && pcAdsMax_[size_t(a)] > ca) { ca = pcAdsMax_[size_t(a)]; da = 0; }
    e.ads = ca; e.dads = da;
    e.rk = R.maxAdsorption > 0 ? 1.0 + (R.rrf - 1.0) * ca / R.maxAdsorption : 1.0;
    return e;
}

// polymer mass in a cell for concentration c, from the current flow state (cf_)
double Simulator::polyMass(int a, double c, double* dM) const {
    const PolymerInput& P = case_.polymer;
    const PolymerRock& R = P.plyrock[std::min(size_t(satReg_[size_t(a)]), P.plyrock.size() - 1)];
    const double wat = cf_[size_t(a)].acc[PW].v;                    // PV(p) * bw * Sw  (surface water volume)
    const double pv = fl(a).pvMult(p_[size_t(a)]) * grid_.poreVolRef[size_t(a)];
    const double rockVol = std::max(0.0, grid_.bulkVol[size_t(grid_.globalOf[size_t(a)])] - pv);
    const PolyEval e = polyEval(a, c);
    const double m = std::max(wat, 1e-10 * pv) * c * (1.0 - R.deadPoreVolume) + rockVol * R.rockDensity * e.ads;
    if (dM) *dM = std::max(wat, 1e-10 * pv) * (1.0 - R.deadPoreVolume) + rockVol * R.rockDensity * e.dads;
    return m;
}

void Simulator::polymerSetup() {
    poly_ = case_.phases.polymer && W_;
    pc_.clear(); pcAdsMax_.clear(); pcMassOld_.clear(); wmul_.clear(); pcInj_.clear();
    pcCum_[0] = pcCum_[1] = 0;
    if (!poly_) return;
    PolymerInput& P = case_.polymer;
    if (P.plyvisc.empty()) throw std::runtime_error("POLYMER: PLYVISC (polymer viscosity multiplier) is required");
    if (P.plyrock.empty()) throw std::runtime_error("POLYMER: PLYROCK (polymer-rock properties) is required");
    if (P.plyads.empty()) throw std::runtime_error("POLYMER: PLYADS (polymer adsorption) is required");
    if (P.cmax.empty() || P.mixing.empty()) throw std::runtime_error("POLYMER: PLYMAX and PLMIXPAR are required");
    for (double v : P.cmax) if (!(v > 0)) throw std::runtime_error("POLYMER: PLYMAX (maximum polymer concentration) must be positive");
    for (double v : P.mixing) if (v < 0 || v > 1) throw std::runtime_error("POLYMER: the PLMIXPAR mixing parameter must be in [0, 1]");
    for (auto& t : P.plyvisc) if (t.size() < 2 || t.size() % 2) throw std::runtime_error("POLYMER: invalid PLYVISC table");
    for (auto& t : P.plyads) if (t.size() < 2 || t.size() % 2) throw std::runtime_error("POLYMER: invalid PLYADS table");
    const int na = grid_.numActive();
    mixReg_.assign(size_t(na), 0);
    if (!case_.grid.plmixnum.empty()) {
        const int nmix = P.numMixRegions();
        for (int a = 0; a < na; ++a) {
            const size_t g = size_t(grid_.globalOf[size_t(a)]);
            const int r = g < case_.grid.plmixnum.size() ? case_.grid.plmixnum[g] - 1 : 0;
            if (r < 0 || r >= nmix) throw std::runtime_error("POLYMER: PLMIXNUM region out of range (PLYMAX / PLMIXPAR records)");
            mixReg_[size_t(a)] = r;
        }
    }
    // shear (PLYSHLOG): log velocities and reference multipliers converted as in OPM Flow
    pcShear_ = !P.plyshlog.empty();
    shLogV_.clear(); shM_.clear(); pcArea_.clear(); pcZw_.clear(); pcZp_.clear();
    if (pcShear_) {
        for (size_t r = 0; r < P.plyshlog.size(); ++r) {
            const auto& S = P.plyshlog[r];
            const auto& vis = P.plyvisc[std::min(r, P.plyvisc.size() - 1)];
            const double pref = tab(vis, S.refConc, nullptr);
            std::vector<double> lv, m;
            for (size_t i = 0; i + 1 < S.rows.size(); i += 2) {
                lv.push_back(std::log(S.rows[i]));
                m.push_back(std::fabs(pref - 1.0) > 1e-14 ? (S.rows[i + 1] * pref - 1.0) / (pref - 1.0) : 1.0);
            }
            shLogV_.push_back(lv); shM_.push_back(m);
        }
        const size_t nc = grid_.conns.size();
        pcArea_.assign(nc, 0.0); pcZw_.assign(nc, 1.0); pcZp_.assign(nc, 1.0);
        for (size_t ci = 0; ci < nc; ++ci) {
            const Connection& cn = grid_.conns[ci];
            const int ga = grid_.globalOf[size_t(cn.a)], gb = grid_.globalOf[size_t(cn.b)];
            const char d = (cn.dir == 'X' || cn.dir == 'Y' || cn.dir == 'Z') ? cn.dir : 'Z';
            const bool aLow = ga < gb;
            const double A1v = grid_.faceArea(ga, d, !aLow), A2v = grid_.faceArea(gb, d, aLow);
            pcArea_[ci] = std::max(1e-12, 0.5 * (A1v + A2v));
        }
    }
    for (auto& W : wells_) { W.zwConn.assign(W.cell.size(), 1.0); W.wmulConn.assign(W.cell.size(), 1.0); }
    pc_.assign(size_t(na), 0.0);
    if (!case_.grid.spoly.empty())
        for (int a = 0; a < na; ++a) {
            const size_t g = size_t(grid_.globalOf[size_t(a)]);
            if (g < case_.grid.spoly.size()) pc_[size_t(a)] = std::max(0.0, case_.grid.spoly[g]);
        }
    pcAdsMax_.assign(size_t(na), 0.0);
    for (int a = 0; a < na; ++a) pcAdsMax_[size_t(a)] = polyEval(a, pc_[size_t(a)]).ads;
    pcMassOld_.assign(size_t(na), 0.0);
    wmul_.assign(size_t(na), 1.0);
    pcInj_.assign(wells_.size(), 0.0);
    polymerMultipliers();
    // scalar Jacobian pattern: self + neighbours (sorted, unique)
    pcJ_ = CSR{}; pcJ_.n = na; pcJ_.rowPtr.assign(size_t(na) + 1, 0);
    pcSelf_.assign(size_t(na), 0); pcAdj_.assign(grid_.adjNbr.size(), 0);
    std::vector<int> cols;
    for (int a = 0; a < na; ++a) {
        cols.clear(); cols.push_back(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) cols.push_back(grid_.adjNbr[size_t(e)]);
        std::sort(cols.begin(), cols.end()); cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
        const int base = int(pcJ_.col.size());
        for (int q : cols) pcJ_.col.push_back(q);
        pcJ_.rowPtr[size_t(a) + 1] = int(pcJ_.col.size());
        auto pos = [&](int q) { return base + int(std::lower_bound(cols.begin(), cols.end(), q) - cols.begin()); };
        pcSelf_[size_t(a)] = pos(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) pcAdj_[size_t(e)] = pos(grid_.adjNbr[size_t(e)]);
    }
    pcJ_.val.assign(pcJ_.col.size(), 0.0);
    pcJ_.finalizeDiag();
    std::string info = "Polymer flooding: PLYMAX " + f6(P.cmax[0]) + ", Todd-Longstaff mixing " + f6(P.mixing[0]) +
                       (P.numMixRegions() > 1 ? " (" + std::to_string(P.numMixRegions()) + " PLMIXNUM regions)" : std::string());
    const PolymerRock& R = P.plyrock[0];
    info += ", dead pore volume " + f6(R.deadPoreVolume) + ", RRF " + f6(R.rrf) +
            (R.adsIndex == 2 ? ", irreversible adsorption" : ", reversible adsorption");
    if (pcShear_) info += ", shear thinning / thickening (PLYSHLOG)";
    log_(info);
}

// OPM Flow / ECLIPSE shear model: the sheared velocity v solves log(v) + log(Z(v)) = log(v0) with
// Z = (1 + (P - 1) M(v)) / P, P = viscosity multiplier at c; returns Z (1 = no effect)
double Simulator::shearFactor(int pvtReg, double c, double v0) const {
    if (!pcShear_) return 1.0;
    const size_t r = std::min(size_t(pvtReg), shLogV_.size() - 1);
    const auto& vis = case_.polymer.plyvisc[std::min(size_t(pvtReg), case_.polymer.plyvisc.size() - 1)];
    const double P = tab(vis, c, nullptr);
    if (std::fabs(P - 1.0) < 1e-14 || !(std::fabs(v0) > 0)) return 1.0;
    const auto& lv = shLogV_[r]; const auto& M = shM_[r];
    const double l0 = std::log(std::fabs(v0));
    if (l0 < lv[0]) return 1.0;
    std::vector<double> t(2 * lv.size());
    for (size_t i = 0; i < lv.size(); ++i) { t[2 * i] = lv[i]; t[2 * i + 1] = std::log(std::max(1e-300, (1.0 + (P - 1.0) * M[i]) / P)); }
    double u = l0;
    for (int it = 0; it < 30; ++it) {
        double d = 0; const double f = u + tab(t, u, &d) - l0;
        if (std::fabs(f) < 1e-12) break;
        u -= f / (1.0 + d);
    }
    return std::exp(tab(t, u, nullptr));
}

void Simulator::polymerControls(int step) {
    if (!poly_ || step < 0 || step >= int(case_.schedule.size())) return;
    for (const WellPolymer& wp : case_.schedule[size_t(step)].wpolymer)
        for (size_t w = 0; w < wells_.size(); ++w) {
            const std::string& name = case_.wells[size_t(wells_[w].input)].name;
            const std::string& pat = wp.well;
            bool hit = pat == name;
            if (!hit && !pat.empty() && pat.back() == '*') hit = name.compare(0, pat.size() - 1, pat, 0, pat.size() - 1) == 0;
            if (hit) pcInj_[w] = std::max(0.0, wp.conc);
        }
}

void Simulator::polymerMultipliers() {
    const int na = grid_.numActive();
    for (int a = 0; a < na; ++a) { const PolyEval e = polyEval(a, pc_[size_t(a)]); wmul_[size_t(a)] = e.wvc / e.rk; }
    for (auto& W : wells_) {
        W.wmulConn.assign(W.cell.size(), 1.0);
        if (W.zwConn.size() != W.cell.size()) W.zwConn.assign(W.cell.size(), 1.0);
        for (size_t k = 0; k < W.cell.size(); ++k) {
            double base = 1.0;
            if (W.type != WellType::Producer) {          // injector: fully mixed polymer viscosity in the wellbore
                const int a = W.cell[k];
                const auto& vis = case_.polymer.plyvisc[std::min(size_t(pvtReg_[size_t(a)]), case_.polymer.plyvisc.size() - 1)];
                base = 1.0 / (polyEval(a, pc_[size_t(a)]).wvc * std::max(1e-12, tab(vis, pc_[size_t(a)], nullptr)));
            }
            W.wmulConn[k] = base / W.zwConn[k];
        }
    }
}

double Simulator::polymerTransport(double dt) {
    const int na = grid_.numActive();
    const double g = U_.gravity;
    // water surface flux per connection (a -> b) at the converged flow state, without the shear factor
    const size_t nc = grid_.conns.size();
    pcFw_.assign(nc, 0.0);
    for (size_t ci = 0; ci < nc; ++ci) {
        const Connection& cn = grid_.conns[ci];
        const CellFlux& A = cf_[size_t(cn.a)]; const CellFlux& B = cf_[size_t(cn.b)];
        const double rho = 0.5 * (A.rho[PW].v + B.rho[PW].v);
        double dphi = A.P[PW].v - B.P[PW].v + rho * g * cn.dDepth;
        if (!thp_.empty() && thp_[ci] > 0) { const double t = thp_[ci]; dphi = dphi > t ? dphi - t : dphi < -t ? dphi + t : 0.0; }
        const double mob = dphi >= 0 ? A.mob[PW].v : B.mob[PW].v;
        pcFw_[ci] = mob * dphi * cn.trans;
    }
    double cScale = 0;
    for (double v : case_.polymer.cmax) cScale = std::max(cScale, v);
    for (double v : pcInj_) cScale = std::max(cScale, v);
    std::vector<double> c = pc_, R(static_cast<size_t>(na)), dx, rhs(static_cast<size_t>(na));
    std::vector<PolyEval> E(static_cast<size_t>(na));
    // flux used by the flow solve = unsheared flux / water shear factor (lagged); the polymer flux also divides by
    // the polymer shear factor
    auto zp = [&](size_t ci) { return pcZp_.empty() ? 1.0 : pcZp_[ci]; };
    bool ok = false;
    for (int it = 0; it < 25; ++it) {
        std::fill(pcJ_.val.begin(), pcJ_.val.end(), 0.0);
        for (int a = 0; a < na; ++a) {
            double dM = 0; const double M = polyMass(a, c[size_t(a)], &dM);
            R[size_t(a)] = (M - pcMassOld_[size_t(a)]) / dt;
            pcJ_.val[size_t(pcSelf_[size_t(a)])] += dM / dt;
            E[size_t(a)] = polyEval(a, c[size_t(a)]);
        }
        for (int a = 0; a < na; ++a)
            for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) {
                const size_t ci = size_t(grid_.adjConn[size_t(e)]);
                const Connection& cn = grid_.conns[ci];
                const double F0 = cn.a == a ? pcFw_[ci] : -pcFw_[ci];   // out of a
                const double F = F0 / zp(ci);                            // OPM: unsheared water flux / polymer shear factor
                const int b = grid_.adjNbr[size_t(e)];
                const int up = F0 >= 0 ? a : b;
                const PolyEval& U = E[size_t(up)];
                const double cu = c[size_t(up)];
                R[size_t(a)] += F * U.pvc * cu;
                pcJ_.val[size_t(up == a ? pcSelf_[size_t(a)] : pcAdj_[size_t(e)])] += F * (U.pvc + cu * U.dpvc);
            }
        for (size_t w = 0; w < wells_.size(); ++w) {
            const WellRT& W = wells_[w];
            if (!W.open) continue;
            for (size_t k = 0; k < W.cell.size() && k < W.qwConn.size(); ++k) {
                const int a = W.cell[k]; const double q = W.qwConn[k];
                if (q > 0) {
                    const PolyEval& U = E[size_t(a)];
                    R[size_t(a)] += q * U.pvc * c[size_t(a)];
                    pcJ_.val[size_t(pcSelf_[size_t(a)])] += q * (U.pvc + c[size_t(a)] * U.dpvc);
                } else if (q < 0) R[size_t(a)] += q * pcInj_[w];
            }
        }
        // convergence: residual as a concentration change (R dt / (PV bw)) relative to the concentration scale
        double rmax = 0;
        for (int a = 0; a < na; ++a) {
            const double s = dt / std::max(cf_[size_t(a)].acc[PW].v, 1e-10 * grid_.poreVolRef[size_t(a)]);
            rmax = std::max(rmax, std::fabs(R[size_t(a)]) * s);
        }
        if (rmax <= 1e-7 * std::max(cScale, 1e-12)) { ok = true; break; }
        for (int a = 0; a < na; ++a) {
            const double s = dt / grid_.poreVolRef[size_t(a)];
            rhs[size_t(a)] = -R[size_t(a)] * s;
            for (int q = pcJ_.rowPtr[size_t(a)]; q < pcJ_.rowPtr[size_t(a) + 1]; ++q) pcJ_.val[size_t(q)] *= s;
        }
        pcIlu_.factor(pcJ_);
        dx.assign(size_t(na), 0.0);
        LinResult lr = bicgstab(pcJ_, pcIlu_, rhs, dx, 1e-9, 500);
        if (!lr.converged) { LinResult lg = gmres(pcJ_, pcIlu_, rhs, dx, 1e-9, 500, 40); if (!lg.converged && lg.relResidual > 1e-4) break; }
        double dmax = 0;
        for (int a = 0; a < na; ++a) {
            double d = dx[size_t(a)];
            if (!std::isfinite(d)) { dmax = -1; break; }
            const double lim = 0.5 * cScale;                         // Appleyard-type chop of large updates
            if (d > lim) d = lim; else if (d < -lim) d = -lim;
            c[size_t(a)] = std::max(0.0, c[size_t(a)] + d);
            dmax = std::max(dmax, std::fabs(d));
        }
        if (dmax < 0) break;
        if (dmax <= 1e-10 * std::max(cScale, 1e-12)) { ok = true; break; }
    }
    if (!ok) return -1.0;
    pc_.swap(c);
    // well polymer rates
    for (size_t w = 0; w < wells_.size(); ++w) {
        WellRT& W = wells_[w];
        double inj = 0, prod = 0;
        if (W.open)
            for (size_t k = 0; k < W.cell.size() && k < W.qwConn.size(); ++k) {
                const int a = W.cell[k]; const double q = W.qwConn[k];
                if (q > 0) prod += q * pc_[size_t(a)] * polyEval(a, pc_[size_t(a)]).pvc;
                else inj -= q * pcInj_[w];
            }
        W.rates.polyInj = inj; W.rates.polyProd = prod;
    }
    double ch = 0;
    auto upd = [&](double& old, double v) { ch = std::max(ch, std::fabs(v - old) / std::max(std::fabs(old), 1e-12)); old = v; };
    // shear factors (OPM Flow definitions) from the unsheared water velocity and the new concentrations
    if (pcShear_) {
        auto poro = [&](int a) {
            const double pv = fl(a).pvMult(p_[size_t(a)]) * grid_.poreVolRef[size_t(a)];
            return pv / std::max(grid_.bulkVol[size_t(grid_.globalOf[size_t(a)])], 1e-30);
        };
        auto swOf = [&](int a) { return iSw_ >= 0 ? sw_[size_t(a)] : 1.0; };
        auto swcrOf = [&](int a) { const auto* e = ep(a); return e ? e->swcr : st(a).tableEndPoints().swcr; };
        for (size_t ci = 0; ci < nc; ++ci) {
            const Connection& cn = grid_.conns[ci];
            const int up = pcFw_[ci] >= 0 ? cn.a : cn.b;
            const double bw = std::max(1e-12, cf_[size_t(up)].rho[PW].v / fl(up).rhoSC().water);
            const double den = std::max(0.5 * (poro(cn.a) + poro(cn.b)) * (swOf(up) - swcrOf(up)), 1e-12);
            const double v0 = pcFw_[ci] / bw / pcArea_[ci] / den;
            const double cu = pc_[size_t(up)];
            const int pr = pvtReg_[size_t(up)];
            upd(pcZw_[ci], shearFactor(pr, cu, v0));
            upd(pcZp_[ci], shearFactor(pr, cu, v0 * polyEval(up, cu).pvc));
        }
        for (size_t w = 0; w < wells_.size(); ++w) {
            WellRT& W = wells_[w];
            if (!W.open) continue;
            const bool inj = W.type != WellType::Producer;
            for (size_t k = 0; k < W.cell.size() && k < W.qwConn.size(); ++k) {
                const int a = W.cell[k];
                double z = 1.0;
                if (!(inj && pcInj_[w] == 0.0) && W.qwConn[k] != 0.0) {
                    const double q0 = W.qwConn[k] * W.zwConn[k];            // rate without the shear factor
                    const double bw = std::max(1e-12, cf_[size_t(a)].rho[PW].v / fl(a).rhoSC().water);
                    const double area = 2.0 * M_PI * std::sqrt(W.re[k] * W.rw[k]) * W.plen[k];
                    const double den = std::max(area * poro(a) * (swOf(a) - swcrOf(a)), 1e-12);
                    z = shearFactor(pvtReg_[size_t(a)], pc_[size_t(a)], q0 / bw / den);
                }
                upd(W.zwConn[k], z);
            }
        }
    }
    // new water mobility multipliers (cells and well connections) and their largest relative change
    std::vector<double> oldConn;
    for (auto& W : wells_) oldConn.insert(oldConn.end(), W.wmulConn.begin(), W.wmulConn.end());
    for (int a = 0; a < na; ++a) {
        const PolyEval e = polyEval(a, pc_[size_t(a)]);
        upd(wmul_[size_t(a)], e.wvc / e.rk);
    }
    polymerMultipliers();
    size_t q = 0;
    for (auto& W : wells_) for (double v : W.wmulConn) { double o = oldConn[q++]; upd(o, v); }
    return ch;
}

void Simulator::polymerAccept() {
    if (!poly_) return;
    const int na = grid_.numActive();
    const PolymerInput& P = case_.polymer;
    for (int a = 0; a < na; ++a) {
        const PolymerRock& R = P.plyrock[std::min(size_t(satReg_[size_t(a)]), P.plyrock.size() - 1)];
        if (R.adsIndex == 2) pcAdsMax_[size_t(a)] = std::max(pcAdsMax_[size_t(a)], polyEval(a, pc_[size_t(a)]).ads);
        pcMassOld_[size_t(a)] = polyMass(a, pc_[size_t(a)], nullptr);
    }
}

} // namespace ressim
