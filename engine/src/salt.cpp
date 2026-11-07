// ResSim - C29 salt precipitation and dry-out in the black-oil engine (TOUGH3 EWASG-type)
//
// Input: '-- RESSIM: SALT c0 [cSol [rhoSalt [phiR [Gamma [T_C [cInj]]]]]]' (Case::salt).
//
// Model (operator split after every converged flow step, the order the TOUGH EWASG equations are lagged in sequential
// reactive-transport codes):
//  1. dissolved NaCl is the internal water-phase tracer NACL (kg per surface m3 of water) of the tracer module: implicit
//     upwind transport with the converged water fluxes, injection concentration cInj (or WTRACER NACL ...);
//  2. local solubility equilibrium per cell: excess above cSol precipitates as solid, solid redissolves into brine
//     undersaturated in salt; S_s = solid volume / reference pore volume;
//  3. dry-out: gas entering a water-bearing cell undersaturated in water vapour picks up water to saturation (local
//     equilibrium; humidity y of the gas leaving dry cells from the inflow mixing, solved by Gauss-Seidel sweeps over the
//     converged gas fluxes; injected gas is dry, y = 0). Evaporation rate E = omega * sum_in q_g,res (1 - y_up),
//     omega = p_sat(T) M_w / (R T rho_w) (ideal gas, vapour volume of liquid water per reservoir volume of gas). E is the
//     water sink of the next flow step, applied as E * S_w / (S_w + 0.005) (implicit in S_w: no negative saturations);
//     the salt stays behind, so the brine concentrates and precipitates (the EWASG dry-out mechanism);
//  4. porosity phi (1 - S_s) through the pore-volume reference, permeability by the Verma & Pruess (1988) tubes-in-series
//     model on the connection transmissibilities (harmonic of the cell factors) and the well connection factors.
// Limitations: the vapour is not carried as a component (it leaves with the gas; no condensation), the PV / T factors are
// applied relative to the values at the start of the run (a two-way GEOMECH / FAULTMECH permeability coupling of the same
// cells is last-writer-wins), halite is the only salt.
#include "ressim/simulator.h"
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace ressim {

namespace {
std::string fmt(const char* f, double a = 0, double b = 0, double c = 0, double d = 0) {
    char buf[512]; std::snprintf(buf, sizeof buf, f, a, b, c, d); return buf;
}
// saturated vapour pressure of water, Pa (Wagner & Pruss 2002 / IAPWS-IF97 region 4 auxiliary equation)
double psatWater(double TC) {
    const double T = TC + 273.15, Tc = 647.096, pc = 22.064e6;
    const double t = 1 - T / Tc;
    if (t <= 0) return pc;
    const double a1 = -7.85951783, a2 = 1.84408259, a3 = -11.7866497, a4 = 22.6807411, a5 = -15.9618719, a6 = 1.80122502;
    const double s = (Tc / T) * (a1 * t + a2 * std::pow(t, 1.5) + a3 * t * t * t + a4 * std::pow(t, 3.5) + a5 * std::pow(t, 4) + a6 * std::pow(t, 7.5));
    return pc * std::exp(s);
}
}  // namespace

// Verma & Pruess (1988) permeability factor for a solid fraction Ss of the pore space
static double vermaPruess(double Ss, double phiR, double G) {
    const double phiRel = 1.0 - Ss;
    if (phiRel <= phiR) return 0.0;
    const double th = (phiRel - phiR) / (1.0 - phiR);
    const double w = 1.0 + (1.0 / G) / (1.0 / phiR - 1.0);
    const double r = th / (th + w - 1.0);
    return th * th * (1.0 - G + G / (w * w)) / (1.0 - G + G * r * r);
}

void Simulator::saltSetup() {
    if (!saltPv0_.empty() && saltPv0_.size() <= grid_.poreVolRef.size() && saltT0_.size() == grid_.conns.size()) {   // re-initialisation
        for (size_t a = 0; a < saltPv0_.size(); ++a) grid_.poreVolRef[a] = saltPv0_[a];
        for (size_t c = 0; c < saltT0_.size(); ++c) grid_.conns[c].trans = saltT0_[c];
        for (WellRT& w : wells_)
            for (size_t k = 0; k < w.cell.size(); ++k) {
                const size_t a = size_t(w.cell[k]);
                if (a >= saltK_.size() || !(saltK_[a] > 0)) continue;
                if (k < w.wi.size()) w.wi[k] /= saltK_[a];
                if (k < w.wiBase.size()) w.wiBase[k] /= saltK_[a];
            }
    }
    saltS_.clear(); saltK_.clear(); saltEvap_.clear(); saltPvF_.clear(); saltTF_.clear();
    saltEvapCum_ = 0; saltOmega_ = 0; saltSolidPrev_ = 0;
    saltTr_ = -1;
    const auto& S = case_.salt;
    if (!S.on) return;
    if (eq_[PW] < 0) { log_("NOTE: SALT needs the water phase - ignored"); return; }
    for (size_t t = 0; t < trInfo_.size(); ++t) if (trInfo_[t].name == "NACL") saltTr_ = int(t);
    if (saltTr_ < 0) { log_("NOTE: SALT: internal NACL tracer missing - ignored"); return; }
    const size_t na = size_t(grid_.numActive());
    saltS_.assign(na, 0.0); saltK_.assign(na, 1.0); saltEvap_.assign(na, 0.0); saltPvF_.assign(na, 1.0);
    saltTF_.assign(grid_.conns.size(), 1.0);
    saltPv0_.assign(grid_.poreVolRef.begin(), grid_.poreVolRef.begin() + std::min(na, grid_.poreVolRef.size()));
    saltT0_.resize(grid_.conns.size());
    for (size_t c = 0; c < grid_.conns.size(); ++c) saltT0_[c] = grid_.conns[c].trans;
    const double TK = S.tempC + 273.15;
    saltOmega_ = G_ ? psatWater(S.tempC) * 0.018015 / (8.314462618 * TK * 1000.0) : 0.0;
    log_(fmt("Salt (EWASG-type): initial NaCl %.4g kg/m3 water, solubility %.4g kg/m3, halite %.4g kg/m3, ", S.c0, S.cSol, S.rhoSalt) +
         fmt("Verma-Pruess phiR %.3g Gamma %.3g, ", S.phiR, S.gamma) +
         fmt("T %.4g C: vapour capacity %.4g m3 water / m3 gas", S.tempC, saltOmega_) + (G_ ? "" : " (no gas phase: no dry-out)"));
    saltEquilibrium();   // an initial brine above the solubility precipitates at t = 0
    saltApply();
}

double Simulator::saltSolidMass() const {
    if (saltS_.empty()) return 0.0;
    const double vf = case_.units == UnitSystem::Field ? 0.158987294928 : 1.0;
    double m = 0;
    for (size_t a = 0; a < saltS_.size() && a < saltPv0_.size(); ++a) m += saltS_[a] * case_.salt.rhoSalt * saltPv0_[a] * vf;
    return m;
}

// solubility equilibrium of the NACL tracer and the solid per cell
void Simulator::saltEquilibrium() {
    if (saltTr_ < 0 || saltS_.empty()) return;
    const auto& S = case_.salt;
    const size_t na = saltS_.size();
    const double vf = case_.units == UnitSystem::Field ? 0.158987294928 : 1.0;   // bbl -> m3
    double* c = &trC_[size_t(saltTr_) * na];
    if (trVolOld_.size() < (size_t(saltTr_) + 1) * na) tracerVolumes();
    const double* cap = &trVolOld_[size_t(saltTr_) * na];
    for (size_t a = 0; a < na; ++a) {
        const double pvm = saltPv0_[a] * vf;
        if (!(pvm > 0)) continue;
        const double capm = cap[a] * vf;                       // surface water volume, m3
        double Ms = saltS_[a] * S.rhoSalt * pvm;               // solid salt, kg
        if (capm < 1e-9 * pvm) { Ms += std::max(0.0, c[a]) * capm; c[a] = 0.0; }   // dried out: everything is solid
        else if (c[a] > S.cSol) { Ms += (c[a] - S.cSol) * capm; c[a] = S.cSol; }
        else if (Ms > 0) { const double d = std::min(Ms, (S.cSol - c[a]) * capm); Ms -= d; c[a] += d / capm; }
        saltS_[a] = std::clamp(Ms / (S.rhoSalt * pvm), 0.0, 0.999);
    }
}

// porosity / permeability factors from the solid fractions (relative to the start-of-run values: restart-exact)
void Simulator::saltApply() {
    if (saltS_.empty()) return;
    const auto& S = case_.salt;
    const size_t na = saltS_.size();
    std::vector<double> kOld = saltK_;
    for (size_t a = 0; a < na; ++a) {
        saltPvF_[a] = std::max(1e-3, 1.0 - saltS_[a]);
        saltK_[a] = std::max(1e-8, vermaPruess(saltS_[a], S.phiR, S.gamma));
        if (a < grid_.poreVolRef.size()) grid_.poreVolRef[a] = saltPv0_[a] * saltPvF_[a];
    }
    for (size_t ci = 0; ci < grid_.conns.size() && ci < saltT0_.size(); ++ci) {
        const Connection& cn = grid_.conns[ci];
        const double ka = saltK_[size_t(cn.a)], kb = saltK_[size_t(cn.b)];
        saltTF_[ci] = 2 * ka * kb / (ka + kb);
        grid_.conns[ci].trans = saltT0_[ci] * saltTF_[ci];
    }
    for (WellRT& w : wells_)
        for (size_t k = 0; k < w.cell.size(); ++k) {
            const size_t a = size_t(w.cell[k]);
            if (a >= na || kOld[a] == saltK_[a]) continue;
            const double r = saltK_[a] / kOld[a];
            if (k < w.wi.size()) w.wi[k] *= r;
            if (k < w.wiBase.size()) w.wiBase[k] *= r;
        }
}

void Simulator::saltStep(double dt) {
    if (saltTr_ < 0 || saltS_.empty() || dt <= 0) return;
    const size_t na = saltS_.size();
    const double g = U_.gravity;
    // evaporated water of the step just finished (the sink as applied in the residual)
    std::vector<double> sw(na, 0.0), bw(na, 1.0), bg(na, 1.0);
    for (size_t a = 0; a < na; ++a) {
        const CellFull C = evalCell(int(a));
        sw[a] = C.S[PW].v; bw[a] = std::max(C.b[PW].v, 1e-30);
        if (G_ && eq_[PG] >= 0) bg[a] = std::max(C.b[PG].v, 1e-30);
    }
    for (size_t a = 0; a < na; ++a) if (saltEvap_[a] > 0) saltEvapCum_ += saltEvap_[a] * sw[a] / (sw[a] + 0.005) * dt;
    saltEquilibrium();
    // dry-out: humidity sweeps over the converged gas fluxes
    std::fill(saltEvap_.begin(), saltEvap_.end(), 0.0);
    if (G_ && eq_[PG] >= 0 && saltOmega_ > 0) {
        const size_t nc = grid_.conns.size();
        std::vector<double> Fr(nc, 0.0), qinj(na, 0.0);
        for (size_t ci = 0; ci < nc; ++ci) {
            const Connection& cn = grid_.conns[ci];
            const CellFlux& A = cf_[size_t(cn.a)]; const CellFlux& B = cf_[size_t(cn.b)];
            const double rho = 0.5 * (A.rho[PG].v + B.rho[PG].v);
            const double dphi = A.P[PG].v - B.P[PG].v + rho * g * cn.dDepth;
            const double F = (dphi >= 0 ? A.mob[PG].v : B.mob[PG].v) * dphi * cn.trans;    // surface gas a -> b
            Fr[ci] = F / bg[size_t(dphi >= 0 ? cn.a : cn.b)];                            // reservoir volume
        }
        for (const WellRT& W : wells_) {
            if (!W.open || W.type != WellType::GasInjector) continue;
            for (size_t k = 0; k < W.cell.size(); ++k) {
                if (W.wi[k] <= 0) continue;
                const int a = W.cell[k];
                const CellFull C = evalCell(a);
                const double dd = W.pConn(k, g) - C.P[PG].v;
                if (dd <= 0) continue;
                double lt = 0;
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) lt += C.mob[p].v / std::max(C.b[p].v, 1e-30);
                qinj[size_t(a)] += W.wi[k] * lt * dd;                                         // reservoir volume of dry gas
            }
        }
        std::vector<double> y(na, 1.0), qin(na, 0.0), dry(na, 0.0);
        const double wetEps = 1e-4;
        for (int sweep = 0; sweep < 400; ++sweep) {
            for (size_t a = 0; a < na; ++a) { qin[a] = qinj[a]; dry[a] = qinj[a]; }
            for (size_t ci = 0; ci < nc; ++ci) {
                const Connection& cn = grid_.conns[ci];
                if (Fr[ci] > 0) { qin[size_t(cn.b)] += Fr[ci]; dry[size_t(cn.b)] += Fr[ci] * (1 - y[size_t(cn.a)]); }
                else if (Fr[ci] < 0) { qin[size_t(cn.a)] -= Fr[ci]; dry[size_t(cn.a)] -= Fr[ci] * (1 - y[size_t(cn.b)]); }
            }
            double change = 0;
            for (size_t a = 0; a < na; ++a) {
                const double yn = sw[a] > wetEps ? 1.0 : (qin[a] > 0 ? std::clamp(1.0 - dry[a] / qin[a], 0.0, 1.0) : 1.0);
                change = std::max(change, std::fabs(yn - y[a]));
                y[a] = yn;
            }
            if (change < 1e-9) break;
        }
        const double* cap = &trVolOld_[size_t(saltTr_) * na];
        for (size_t a = 0; a < na; ++a) {
            if (sw[a] <= wetEps || dry[a] <= 0) continue;
            const double E = saltOmega_ * dry[a] * bw[a];             // surface water volume per day
            saltEvap_[a] = std::min(E, 0.9 * cap[a] / dt);            // at most the water of the cell within one step
        }
    }
    saltApply();
}

}  // namespace ressim
