// ResSim - C35 SOLVENT: the four-component black-oil solvent model (ECLIPSE SOLVENT / OPM Flow blackoilsolventmodules).
//
// Input (kept verbatim, exact deck round trip): RUNSPEC SOLVENT, PROPS PVDS SDENSITY SSFN, SOLUTION SSOL,
// SCHEDULE WSOLVENT (solvent fraction of the injected gas). MISC / PMISC / MSFN / SORWMIS / SGCWMIS / TLMIXPAR /
// TLPMIXPA (miscible solvent, Todd-Longstaff mixing) are kept for OPM Flow; the engine uses the immiscible model.
//
// Formulation (docs/dev/C35_SOLVENT_DESIGN.md, ENGINES_AND_MATH 6.29):
//   flow (fully implicit): water, oil and the total gas-like component G+S (free + dissolved gas + solvent, surface
//   volumes); gas-phase variable SgT = Sg + Ss (state 1) - the saturation functions see SgT as in OPM.
//   lagged composition z = Ns / (Ns + Ng) per cell; the flash in evalCell splits SgT:
//       Ss = z (SgT bg + So bo Rs) / (bs (1 - z) + z bg)      (free gas present, Rs = Rs_sat)
//       Ss = SgT, Sg = 0, Rs = Ss bs (1 - z) / (z So bo)      (gas phase is pure solvent; undersaturated oil)
//   immiscible relative permeabilities: krg = krgT(SgT) SSFN_g(Sg / SgT), krs = krgT(SgT) SSFN_s(Ss / SgT)
//   solvent flux with the gas-phase pressure and the solvent density (SDENSITY x 1/Bs from PVDS).
//   composition (after the flow solve, implicit in z, flow state and potentials frozen):
//       [z N_T - (z N_T)^n] / dt + sum_out T mob_s(z_up) dPhi_s + q_prod(z) - fS q_inj = 0
//   iterated with the flow solve until z stops changing (outer coupling shared with polymer / chemicals).
#include "ressim/simulator.h"
#include "ressim/linsolve.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>

namespace ressim {
namespace {
std::string svUp(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
bool svMatch(const std::string& pat, const std::string& name) {
    const std::string a = svUp(pat), b = svUp(name);
    if (!a.empty() && a.back() == '*') return b.compare(0, a.size() - 1, a, 0, a.size() - 1) == 0;
    return a == b;
}
char svBuf[512];
template <class... A> std::string sfmt(const char* f, A... a) { std::snprintf(svBuf, sizeof svBuf, f, a...); return svBuf; }
}  // namespace

void Simulator::solventSetup() {
    solv_ = solvMisc_ = false;
    zs_.clear(); zsBak_.clear(); ntOld_.clear(); svInjF_.clear(); svB_.clear(); svBMu_.clear(); svRho_.clear();
    ssfnG_.clear(); ssfnS_.clear(); svF_.clear(); svWq_.clear(); svCum_[0] = svCum_[1] = 0; svOuter_ = 0;
    bool on = false;
    for (const ExtraKeyword& k : case_.extraKeywords)
        if (k.section == DeckSection::Runspec && k.name == "SOLVENT") on = true;
    if (!on) return;
    if (!G_) throw std::runtime_error("SOLVENT needs the GAS phase");
    if (daOn_) throw std::runtime_error("SOLVENT: the MPI domain-decomposed black-oil run does not support the solvent model - run without -mpi");
    auto tables = [&](const char* name, size_t ncol) {
        std::vector<std::vector<double>> out;
        for (const ExtraKeyword& k : case_.extraKeywords) {
            if (k.name != name || k.section == DeckSection::Schedule) continue;
            out.clear();
            for (const auto& r : verbatimRecords(k.body)) {
                std::vector<double> v;
                for (const auto& s : r) if (!s.empty()) v.push_back(std::atof(s.c_str()));
                if (v.empty()) { if (!out.empty()) out.push_back(out.back()); continue; }
                if (v.size() % ncol) throw std::runtime_error(sfmt("%s: every table row needs %d values", name, int(ncol)));
                for (size_t i = ncol; i < v.size(); i += ncol)
                    if (!(v[i] > v[i - ncol])) throw std::runtime_error(std::string(name) + ": the first column must increase");
                out.push_back(v);
            }
        }
        return out;
    };
    const auto pvds = tables("PVDS", 3), ssfn = tables("SSFN", 3);
    if (pvds.empty()) throw std::runtime_error("SOLVENT: PVDS (solvent formation volume factor and viscosity) is required");
    if (ssfn.empty()) throw std::runtime_error("SOLVENT: SSFN (gas / solvent relative permeability multipliers) is required");
    for (const auto& t : pvds) {
        Table1D b, bm;
        for (size_t i = 0; i + 2 < t.size(); i += 3) {
            if (!(t[i + 1] > 0) || !(t[i + 2] > 0)) throw std::runtime_error("SOLVENT: PVDS formation volume factors and viscosities must be positive");
            b.x.push_back(t[i]); b.y.push_back(1.0 / t[i + 1]);
            bm.x.push_back(t[i]); bm.y.push_back(1.0 / (t[i + 1] * t[i + 2]));
        }
        svB_.push_back(b); svBMu_.push_back(bm);
    }
    for (const auto& t : ssfn) {
        Table1D g, s;
        g.clampEnds = s.clampEnds = true;
        for (size_t i = 0; i + 2 < t.size(); i += 3) { g.x.push_back(t[i]); g.y.push_back(t[i + 1]); s.x.push_back(t[i]); s.y.push_back(t[i + 2]); }
        ssfnG_.push_back(g); ssfnS_.push_back(s);
    }
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.name != "SDENSITY" || k.section == DeckSection::Schedule) continue;
        svRho_.clear();
        for (const auto& r : verbatimRecords(k.body)) {
            if (r.empty()) { if (!svRho_.empty()) svRho_.push_back(svRho_.back()); continue; }
            svRho_.push_back(verbatimNum(r, 0, 0));
        }
    }
    if (svRho_.empty()) throw std::runtime_error("SOLVENT: SDENSITY (solvent surface density) is required");
    for (const ExtraKeyword& k : case_.extraKeywords)
        if (k.section != DeckSection::Schedule && (k.name == "MISC" || k.name == "TLMIXPAR")) solvMisc_ = true;
    if (solvMisc_) log_("NOTE: SOLVENT with MISCIBLE data (MISC / TLMIXPAR) - the miscible solvent model (Todd-Longstaff mixing, "
                        "miscible relative permeabilities) is not modelled by the ResSim engine; the immiscible SSFN model is used "
                        "(the data are kept verbatim for OPM Flow / ECLIPSE)");
    for (const ExtraKeyword& k : case_.extraKeywords)
        if (k.section == DeckSection::Runspec && (k.name == "CO2SOL" || k.name == "H2SOL"))
            throw std::runtime_error("SOLVENT: " + k.name + " (solvent dissolving in water) is not modelled by the ResSim engine - run the deck with OPM Flow");
    const int na = grid_.numActive();
    solv_ = true;
    zs_.assign(size_t(na), 0.0); ntOld_.assign(size_t(na), 0.0);
    svInjF_.assign(wells_.size(), 0.0);
    svF_.assign(grid_.conns.size(), 0.0);
    svWq_.assign(wells_.size(), {});
    // scalar Jacobian pattern: self + neighbours
    svJ_ = CSR{}; svJ_.n = na; svJ_.rowPtr.assign(size_t(na) + 1, 0);
    svSelf_.assign(size_t(na), 0); svAdj_.assign(grid_.adjNbr.size(), 0);
    std::vector<int> cols;
    for (int a = 0; a < na; ++a) {
        cols.clear(); cols.push_back(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) cols.push_back(grid_.adjNbr[size_t(e)]);
        std::sort(cols.begin(), cols.end()); cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
        const int base = int(svJ_.col.size());
        for (int q : cols) svJ_.col.push_back(q);
        svJ_.rowPtr[size_t(a) + 1] = int(svJ_.col.size());
        auto pos = [&](int q) { return base + int(std::lower_bound(cols.begin(), cols.end(), q) - cols.begin()); };
        svSelf_[size_t(a)] = pos(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) svAdj_[size_t(e)] = pos(grid_.adjNbr[size_t(e)]);
    }
    svJ_.val.assign(svJ_.col.size(), 0.0);
    svJ_.finalizeDiag();
    // SSOL: initial solvent saturation - part of the gas phase (SgT = Sg + Ss, oil reduced), z from the solvent amount
    int nInit = 0;
    if (!case_.grid.ssol.empty()) {
        for (int a = 0; a < na; ++a) {
            const size_t g = size_t(grid_.globalOf[size_t(a)]);
            const double ss = g < case_.grid.ssol.size() ? std::clamp(case_.grid.ssol[g], 0.0, 1.0) : 0.0;
            if (!(ss > 0)) continue;
            const CellFull C = evalCell(a);                 // z = 0: the gas-phase state of the equilibration
            const double ng = std::max(0.0, C.acc[PG].v);
            const double pv = fl(a).pvMult(AD3(p_[size_t(a)])).v * grid_.poreVolRef[size_t(a)];
            const double ns = pv * C.bs.v * ss;
            if (O_ && hasGas_[size_t(a)] == 0) { hasGas_[size_t(a)] = 1; x_[size_t(a)] = ss; }
            else if (hasGas_[size_t(a)] == 2) {}            // no oil: SgT = 1 - Sw already holds the solvent
            else x_[size_t(a)] += ss;
            zs_[size_t(a)] = ns / std::max(ns + ng, 1e-300);
            ++nInit;
        }
        evalAll();
    }
    log_("Solvent model (SOLVENT): " + sfmt("%d PVDS, %d SSFN table(s)", int(svB_.size()), int(ssfnG_.size())) +
         (nInit ? sfmt(", initial solvent (SSOL) in %d cells", nInit) : std::string()) +
         (solvMisc_ ? ", immiscible (miscible data ignored)" : ", immiscible"));
}

void Simulator::solventControls(int step) {
    if (!solv_) return;
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.section != DeckSection::Schedule || k.step != step || k.name != "WSOLVENT") continue;
        for (const auto& r : verbatimRecords(k.body)) {
            if (r.empty() || r[0].empty()) continue;
            const double f = std::clamp(verbatimNum(r, 1, 0), 0.0, 1.0);
            for (size_t w = 0; w < wells_.size(); ++w)
                if (svMatch(r[0], case_.wells[size_t(wells_[w].input)].name)) {
                    if (f > 0 && wells_[w].type != WellType::GasInjector)
                        log_("NOTE: WSOLVENT " + case_.wells[size_t(wells_[w].input)].name + " - only gas injectors inject solvent, ignored");
                    svInjF_[w] = wells_[w].type == WellType::GasInjector ? f : 0.0;
                }
        }
    }
}

double Simulator::solventTotalNT(int a) const {
    double nt = cf_[size_t(a)].acc[PG].v;
    if (!langK_.empty() && langK_[size_t(a)] > 0) { const double P = p_[size_t(a)]; nt -= P * (langK_[size_t(a)] / langPL_) / (P / langPL_ + 1.0); }
    return std::max(nt, 0.0);
}

double Simulator::solventInPlace() const {
    if (!solv_) return 0;
    double s = 0;   // the solvent amount of the last accepted step (exactly what the solvent transport conserves)
    for (double v : ntOld_) s += v;
    return s;
}

void Simulator::solventAccept() {
    if (!solv_) return;
    const int na = grid_.numActive();
    for (int a = 0; a < na; ++a) ntOld_[size_t(a)] = zs_[size_t(a)] * solventTotalNT(a);
}

void Simulator::solventFixStates() {}   // (state appearance handled inside solventTransport)

double Simulator::solventTransport(double dt) {
    if (!solv_ || dt <= 0) return 0;
    const int na = grid_.numActive();
    const size_t nau = size_t(na), nc = grid_.conns.size();
    const double g = U_.gravity;
    const int refP = O_ ? PO : PW;
    // frozen flow state: G+S amount per cell and solvent potential differences per connection
    std::vector<double> nt(nau), dphi(nc, 0.0);
    for (int a = 0; a < na; ++a) nt[size_t(a)] = solventTotalNT(a);
    for (size_t ci = 0; ci < nc; ++ci) {
        const Connection& cn = grid_.conns[ci];
        const CellFlux& A = cf_[size_t(cn.a)]; const CellFlux& B = cf_[size_t(cn.b)];
        double d = A.P[PG].v - B.P[PG].v + 0.5 * (cfS_[size_t(cn.a)][1].v + cfS_[size_t(cn.b)][1].v) * g * cn.dDepth;
        if (!thp_.empty() && thp_[ci] > 0) { const double th = thp_[ci]; d = d > th ? d - th : d < -th ? d + th : 0.0; }
        dphi[ci] = d;
    }
    // solvent mobility of a cell at composition z (flow state frozen) and its derivative
    std::vector<double> mob(nau), dmob(nau);
    auto mobAt = [&](int a, double z) {
        const double keep = zs_[size_t(a)];
        zs_[size_t(a)] = z;
        const double m = evalCell(a).mobS.v;
        zs_[size_t(a)] = keep;
        return m;
    };
    auto evalMob = [&](const std::vector<double>& z) {
        for (int a = 0; a < na; ++a) {
            const size_t u = size_t(a);
            if (G_ && O_ && hasGas_[u] == 0) { mob[u] = dmob[u] = 0; continue; }   // no gas phase: solvent cannot flow out
            const double h = 1e-6;
            const double z0 = z[u], m0 = mobAt(a, z0);
            const double z1 = z0 + h <= 1.0 ? z0 + h : z0 - h;
            mob[u] = m0; dmob[u] = (mobAt(a, z1) - m0) / (z1 - z0);
        }
    };
    // well connections: producer drawdowns (frozen) and injected solvent (fixed)
    struct WC { int cell; double wiDd; double inj; };   // outflow = wiDd * mob(z) (producers), inj = injected solvent (+)
    std::vector<std::vector<WC>> wc(wells_.size());
    for (size_t w = 0; w < wells_.size(); ++w) {
        const WellRT& W = wells_[w];
        if (!W.open) continue;
        const bool prod = W.type == WellType::Producer;
        std::vector<std::pair<int, double>> qi;
        double sq = 0;
        for (size_t k = 0; k < W.cell.size(); ++k) {
            if (W.wi[k] <= 0) continue;
            const int a = W.cell[k];
            const CellFull C = evalCell(a);
            if (prod) {
                const double dd = C.P[refP].v - W.pConn(k, g);
                if (dd > 0) wc[w].push_back({a, W.wi[k] * dd * W.eff, 0.0});
            } else if (W.type == WellType::GasInjector && svInjF_[w] > 0) {
                const double dd = W.pConn(k, g) - C.P[PG].v;
                if (dd <= 0) continue;
                double lt = 0;
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && C.b[p].v > 0) lt += C.mob[p].v / C.b[p].v;
                if (C.bs.v > 0) lt += C.mobS.v / C.bs.v;
                const double f = svInjF_[w];
                const double binj = C.bs.v > 0 ? 1.0 / ((1.0 - f) / C.b[PG].v + f / C.bs.v) : C.b[PG].v;
                const double q = W.wi[k] * lt * binj * dd;
                qi.push_back({a, q}); sq += q;
            }
        }
        if (!qi.empty() && sq > 0) {   // honour the well model's injected rate exactly (last assembly)
            const double tot = W.rates.gasInj > 0 ? W.rates.gasInj : sq;
            for (auto& e : qi) wc[w].push_back({e.first, 0.0, e.second / sq * tot * svInjF_[w] * W.eff});
        }
    }
    std::vector<double> z = zs_, z0 = zs_, R(nau), rhs(nau), dx(nau);
    double scale = 1e-12;
    for (double v : nt) scale = std::max(scale, v);
    bool ok = false;
    for (int it = 0; it < 40; ++it) {
        evalMob(z);
        std::fill(svJ_.val.begin(), svJ_.val.end(), 0.0);
        for (int a = 0; a < na; ++a) {
            const size_t u = size_t(a);
            const bool reg = nt[u] < 1e-12 * scale;   // no G+S in the cell: z is undefined, keep it
            R[u] = reg ? 0.0 : (z[u] * nt[u] - ntOld_[u]) / dt;
            svJ_.val[size_t(svSelf_[u])] += reg ? 1.0 / dt * std::max(1e-6 * scale, 1e-300) : nt[u] / dt;
            if (reg) continue;
            for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) {
                const size_t ci = size_t(grid_.adjConn[size_t(e)]);
                const double d = grid_.conns[ci].a == a ? dphi[ci] : -dphi[ci];   // out of a
                const int b = grid_.adjNbr[size_t(e)];
                const double T = grid_.conns[ci].trans;
                if (d >= 0) { R[u] += T * mob[u] * d; svJ_.val[size_t(svSelf_[u])] += T * dmob[u] * d; }
                else { R[u] += T * mob[size_t(b)] * d; svJ_.val[size_t(svAdj_[size_t(e)])] += T * dmob[size_t(b)] * d; }
            }
        }
        for (size_t w = 0; w < wc.size(); ++w)
            for (const WC& c : wc[w]) {
                const size_t u = size_t(c.cell);
                if (c.wiDd > 0) { R[u] += c.wiDd * mob[u]; svJ_.val[size_t(svSelf_[u])] += c.wiDd * dmob[u]; }
                R[u] -= c.inj;
            }
        double rmax = 0;
        for (int a = 0; a < na; ++a) rmax = std::max(rmax, std::fabs(R[size_t(a)]) * dt / std::max(nt[size_t(a)], 1e-6 * scale));
        if (std::getenv("RESSIM_DEBUG")) std::fprintf(stderr, "      solvent it %d rmax %.3e\n", it, rmax);
        if (rmax <= 1e-9) { ok = true; break; }
        for (int a = 0; a < na; ++a) {
            const size_t u = size_t(a);
            const double s = dt / std::max(nt[u], 1e-6 * scale);
            rhs[u] = -R[u] * s;
            for (int q = svJ_.rowPtr[u]; q < svJ_.rowPtr[u + 1]; ++q) svJ_.val[size_t(q)] *= s;
        }
        svIlu_.factor(svJ_);
        std::fill(dx.begin(), dx.end(), 0.0);
        LinResult lr = bicgstab(svJ_, svIlu_, rhs, dx, 1e-12, 500);
        if (!lr.converged) { LinResult lg = gmres(svJ_, svIlu_, rhs, dx, 1e-12, 500, 40); if (!lg.converged && lg.relResidual > 1e-4) break; }
        double dmax = 0;
        bool bad = false;
        for (int a = 0; a < na; ++a) {
            double d = dx[size_t(a)];
            if (!std::isfinite(d)) { bad = true; break; }
            d = std::clamp(d, -0.2, 0.2);
            const double zn = std::clamp(z[size_t(a)] + d, 0.0, 1.0);
            dmax = std::max(dmax, std::fabs(zn - z[size_t(a)]));
            z[size_t(a)] = zn;
        }
        if (bad) break;
        if (dmax <= 1e-12) { ok = true; break; }
    }
    static const bool dbg = std::getenv("RESSIM_DEBUG") != nullptr;
    if (!ok) { if (dbg) std::fprintf(stderr, "    [solvent transport failed]\n"); return -1.0; }
    // converged solvent fluxes / well rates (foam carried by the solvent, reporting)
    evalMob(z);
    for (size_t ci = 0; ci < nc; ++ci) {
        const Connection& cn = grid_.conns[ci];
        svF_[ci] = (dphi[ci] >= 0 ? mob[size_t(cn.a)] : mob[size_t(cn.b)]) * dphi[ci] * cn.trans;
    }
    for (size_t w = 0; w < wc.size(); ++w) {
        svWq_[w].clear();
        for (const WC& c : wc[w]) {
            const double q = c.wiDd > 0 ? c.wiDd * mob[size_t(c.cell)] : -c.inj;
            if (q != 0) svWq_[w].push_back({double(c.cell), q});
        }
    }
    double ch = 0;
    for (int a = 0; a < na; ++a) ch = std::max(ch, std::fabs(z[size_t(a)] - z0[size_t(a)]));
    zs_.swap(z);
    // the gas phase appears where solvent entered an undersaturated cell: SgT = Ns / (PV bs) (pure-solvent gas phase,
    // the flash then returns the cell's own Rs), and the flow is solved again
    int nApp = 0;
    if (O_)
        for (int a = 0; a < na; ++a) {
            const size_t u = size_t(a);
            if (hasGas_[u] != 0 || zs_[u] <= 1e-10) continue;
            const CellFull C = evalCell(a);
            const double pv = fl(a).pvMult(AD3(p_[u])).v * grid_.poreVolRef[u];
            const double bs = C.bs.v > 0 ? C.bs.v : svB_[std::min(size_t(pvtReg_[u]), svB_.size() - 1)](AD3(p_[u])).v;
            hasGas_[u] = 1; x_[u] = std::min(zs_[u] * nt[u] / (pv * bs), 0.999);
            ++nApp;
        }
    if (nApp) { evalAll(); ch = std::max(ch, 1.0); }
    if (dbg) {
        double zmax = 0; for (double v : zs_) zmax = std::max(zmax, v);
        std::fprintf(stderr, "    [solvent: dz %.3e, z max %.4f, %d cells gain a gas phase]\n", ch, zmax, nApp);
    }
    return ch;
}

} // namespace ressim
