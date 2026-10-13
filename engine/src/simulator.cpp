#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include "ressim/simulator.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <map>

namespace ressim {

using AD4 = AD<4>;

static std::string fmt(const char* f, double a = 0, double b = 0, double c = 0, double d = 0) {
    char buf[512]; std::snprintf(buf, sizeof(buf), f, a, b, c, d); return buf;
}

// ======================================================================== setup
void Simulator::setup(const Case& c) {
    const std::string err = validateCase(c);
    if (!err.empty()) throw std::runtime_error(err);
    case_ = c;
    U_ = Units::of(c.units);
    W_ = c.phases.water; O_ = c.phases.oil; G_ = c.phases.gas;
    if (W_ && G_ && !O_) throw std::runtime_error("Water-gas systems without oil are not supported yet");
    if (!W_ && !O_) throw std::runtime_error("At least water or oil must be active");
    live_ = O_ && G_ && c.phases.disgas;
    vap_ = O_ && G_ && c.phases.vapoil;
    nv_ = 0;
    for (int ph = 0; ph < 3; ++ph) eq_[ph] = -1;
    if (W_) eq_[PW] = nv_++;
    if (O_) eq_[PO] = nv_++;
    if (G_) eq_[PG] = nv_++;
    iSw_ = (W_ && O_) ? 1 : -1;
    iX_ = G_ ? (iSw_ >= 0 ? 2 : 1) : -1;

    grid_.build(c);
    {   // ---- PVTNUM / SATNUM / EQLNUM regions
        const int npvt = c.props.numPvtRegions(), nsat = c.props.numSatRegions();
        fluids_.assign(size_t(npvt), FluidModel{});
        sats_.assign(size_t(nsat), SatFunctions{});
        for (int r = 0; r < npvt; ++r) {
            try { fluids_[size_t(r)].setup(c, r); }
            catch (const std::exception& e) { throw std::runtime_error(fmt("PVT region %g: ", r + 1) + e.what()); }
        }
        for (int r = 0; r < nsat; ++r) {
            try { sats_[size_t(r)].setup(c, r); }
            catch (const std::exception& e) { throw std::runtime_error(fmt("Saturation region %g: ", r + 1) + e.what()); }
        }
        const int na0 = grid_.numActive();
        const int neql = 1 + int(c.equilExtra.size());
        auto mapRegion = [&](const std::vector<int>& arr, int nreg, std::vector<int>& out) {
            out.assign(size_t(na0), 0);
            if (arr.empty()) return;
            for (int a = 0; a < na0; ++a) {
                const size_t q = size_t(grid_.globalOf[size_t(a)]);
                const int r = q < arr.size() ? arr[q] - 1 : 0;
                out[size_t(a)] = (r >= 0 && r < nreg) ? r : 0;
            }
        };
        mapRegion(c.grid.pvtnum, npvt, pvtReg_);
        mapRegion(c.grid.satnum, nsat, satReg_);
        mapRegion(c.grid.eqlnum, neql, eqlReg_);
    }
    setupEndPointScaling(c);
    setNumThreads(c.tuning.threads);

    const int na = grid_.numActive();
    p_.assign(na, 0); sw_.assign(na, 0); x_.assign(na, 0); hasGas_.assign(na, 0);
    accOld_.assign(size_t(na) * 3, 0.0);
    cf_.assign(na, CellFlux{});
    buildWells();
    wellNames_.clear(); for (auto& w : c.wells) wellNames_.push_back(w.name);
    buildAquifers();
    buildPattern();
    polymerSetup();
    useCPR_ = c.tuning.linearSolver == 1 || (c.tuning.linearSolver == 2 && na >= 500);
    time_ = 0; step_ = 0; dtNext_ = c.tuning.initialDt;
    for (double& v : cum_) v = 0;
    summary_.clear(); snaps_.clear(); stats_ = {}; stats_.cpr = useCPR_;
    log_(fmt("Grid %g x %g x %g, active cells %g", grid_.nx, grid_.ny, grid_.nz, na) +
         fmt(", unknowns %g, threads %g", double(J_.n), numThreads()) +
         (useCPR_ ? ", linear solver CPR-AMG" : ", linear solver ILU(0)-BiCGSTAB"));
}

void Simulator::buildWells() {
    const GridInput& g = case_.grid;
    const int na = grid_.numActive();
    cellWells_.assign(na, {});
    wells_.clear();
    auto at = [&](const std::vector<double>& v, int q, double def) { return v.empty() ? def : v[q]; };
    for (int wi = 0; wi < int(case_.wells.size()); ++wi) {
        const WellInput& in = case_.wells[wi];
        WellRT w; w.input = wi; w.type = in.type;
        const int widx = int(wells_.size());
        for (const Completion& cp : in.completions) {
            for (int k = std::min(cp.k1, cp.k2); k <= std::max(cp.k1, cp.k2); ++k) {
                if (cp.i < 0 || cp.i >= g.nx || cp.j < 0 || cp.j >= g.ny || k < 0 || k >= g.nz)
                    throw std::runtime_error("Well " + in.name + ": completion outside the grid");
                // DUALPORO: wells connect to the fracture cell of the block (a matrix-layer completion is moved)
                const int kk = (g.dual == 1 && k < g.nzMatrix()) ? k + g.nzMatrix() : k;
                const int q = g.index(cp.i, cp.j, kk);
                const int a = grid_.activeOf[q];
                if (a < 0 || !cp.open) continue;
                double wi_ = cp.cf;
                const double kx = grid_.effPerm(g.permx, q), ky = g.permy.empty() ? kx : grid_.effPerm(g.permy, q),
                             kz = g.permz.empty() ? kx : grid_.effPerm(g.permz, q);
                const double ntg = at(g.ntg, q, 1.0);
                double k1, k2, d1, d2, h;
                const double DX = grid_.cellDx[q], DY = grid_.cellDy[q], DZ = grid_.cellDz[q];
                if (cp.dir == 'X' || cp.dir == 'x') { k1 = ky; k2 = kz; d1 = DY; d2 = DZ; h = DX; }
                else if (cp.dir == 'Y' || cp.dir == 'y') { k1 = kx; k2 = kz; d1 = DX; d2 = DZ; h = DY; }
                else { k1 = kx; k2 = ky; d1 = DX; d2 = DY; h = DZ * ntg; }
                w.rw.push_back(0.5 * cp.diameter); w.re.push_back(std::sqrt(d1 * d2 / M_PI));
                w.plen.push_back(cp.kh > 0 && k1 > 0 && k2 > 0 ? cp.kh / std::sqrt(k1 * k2) : h);
                if (wi_ <= 0) {
                    const double kh = cp.kh > 0 ? cp.kh : std::sqrt(k1 * k2) * h;
                    const double r12 = std::sqrt(k2 / k1), r21 = std::sqrt(k1 / k2);
                    const double ro = 0.28 * std::sqrt(r12 * d1 * d1 + r21 * d2 * d2) /
                                      (std::pow(k2 / k1, 0.25) + std::pow(k1 / k2, 0.25));
                    const double rw = 0.5 * cp.diameter;
                    wi_ = U_.cDarcy * 2.0 * M_PI * kh / (std::log(ro / rw) + cp.skin);
                    if (!(wi_ > 0)) throw std::runtime_error("Well " + in.name + ": non-positive connection factor");
                }
                wi_ *= in.piMult;
                w.cell.push_back(a); w.wi.push_back(wi_); w.cdepth.push_back(grid_.depth[a]);
                w.wiBase.push_back(wi_); w.cFrom.push_back(std::max(0, cp.fromStep)); w.cUntil.push_back(cp.untilStep);
                auto& cw = cellWells_[a];
                int idx = -1;
                for (int m = 0; m < int(cw.size()); ++m) if (cw[m] == widx) idx = m;
                if (idx < 0) { idx = int(cw.size()); cw.push_back(widx); }
                w.cellWellIdx.push_back(idx);
            }
        }
        w.rowCells = w.cell;
        std::sort(w.rowCells.begin(), w.rowCells.end());
        w.rowCells.erase(std::unique(w.rowCells.begin(), w.rowCells.end()), w.rowCells.end());
        for (int a : w.cell) w.cpos.push_back(int(std::lower_bound(w.rowCells.begin(), w.rowCells.end(), a) - w.rowCells.begin()));
        if (in.refDepth > 0) w.refDepth = in.refDepth;
        else { w.refDepth = 1e30; for (double d : w.cdepth) w.refDepth = std::min(w.refDepth, d); if (w.cdepth.empty()) w.refDepth = 0; }
        w.ctrl.well = in.name; w.ctrl.status = WellStatus::Shut;
        w.open = false;
        wells_.push_back(w);
    }
}

// ======================================================================== aquifers
namespace {
// Infinite-acting radial aquifer dimensionless pressure (Edwardson et al. fit) and its derivative
void carterTracyPD(double tD, double& pd, double& dpd) {
    auto f = [](double t) {
        if (t <= 0) return 0.0;
        if (t < 0.01) return 2.0 * std::sqrt(t / M_PI);
        if (t < 500.0) {
            const double st = std::sqrt(t);
            return (370.529 * st + 137.582 * t + 5.69549 * t * st) / (328.834 + 265.488 * st + 45.2157 * t + t * st);
        }
        return 0.5 * (std::log(t) + 0.80907);
    };
    pd = f(tD);
    if (tD < 0.01) dpd = 1.0 / std::sqrt(M_PI * std::max(tD, 1e-30));
    else if (tD >= 500.0) dpd = 0.5 / tD;
    else { const double h = 1e-5 * tD; dpd = (f(tD + h) - f(tD - h)) / (2 * h); }
}
}

void Simulator::buildAquifers() {
    aquifers_.clear();
    const GridInput& g = case_.grid;
    for (const AquiferInput& in : case_.aquifers) {
        AquiferRT Q; Q.id = in.id; Q.type = in.type;
        std::vector<double> area;
        for (const AquiferConnection& c : in.conns)
            for (int k = c.k1; k <= c.k2; ++k)
                for (int j = c.j1; j <= c.j2; ++j)
                    for (int i = c.i1; i <= c.i2; ++i) {
                        const int q = g.index(i, j, k);
                        const int a = grid_.activeOf[q];
                        if (a < 0) continue;
                        int ni = i, nj = j, nk = k;
                        const int s = c.minus ? -1 : 1;
                        if (c.face == 'X') ni += s; else if (c.face == 'Y') nj += s; else nk += s;
                        if (ni >= 0 && ni < g.nx && nj >= 0 && nj < g.ny && nk >= 0 && nk < g.nz &&
                            grid_.activeOf[g.index(ni, nj, nk)] >= 0) continue;     // internal face: not connected
                        const double A = grid_.faceArea(q, c.face, c.minus) * c.mult;
                        if (A <= 0) continue;
                        Q.cell.push_back(a); area.push_back(A); Q.cdepth.push_back(grid_.depth[a]);
                    }
        if (Q.cell.empty()) { log_(fmt("Warning: aquifer %g has no valid connections - ignored", in.id)); continue; }
        double sa = 0; for (double v : area) sa += v;
        for (double v : area) Q.alpha.push_back(v / sa);
        if (in.depth > 0) Q.datum = in.depth;
        else { double d = 0; for (size_t k = 0; k < Q.cell.size(); ++k) d += Q.alpha[k] * Q.cdepth[k]; Q.datum = d; }
        Q.autoP0 = !(in.p0 > 0);
        Q.muw = in.visc > 0 ? in.visc : case_.props.pvtw.visc;
        Q.p0 = in.p0;
        if (in.type == AquiferType::Fetkovich) { Q.ctV = in.ct * in.volume; Q.J = in.pi; }
        else {
            const double mu = in.visc > 0 ? in.visc : case_.props.pvtw.visc;
            const double c1 = U_.cDarcy / U_.poreVolume;
            Q.Tc = mu * in.poro * in.ct * in.rInner * in.rInner / (in.perm * c1);
            Q.beta = 2.0 * M_PI * U_.poreVolume * in.thickness * (in.angle / 360.0) * in.poro * in.ct * in.rInner * in.rInner;
        }
        aquifers_.push_back(Q);
    }
    if (!aquifers_.empty()) log_(fmt("Analytic aquifers: %g", double(aquifers_.size())));
}

void Simulator::initAquifers() {
    const double g = U_.gravity;
    for (AquiferRT& Q : aquifers_) {
        if (Q.autoP0) {
            double p = 0;
            for (size_t k = 0; k < Q.cell.size(); ++k) {
                const CellFull C = evalCell(Q.cell[k]);
                p += Q.alpha[k] * (C.P[PW].v - C.rho[PW].v * g * (Q.cdepth[k] - Q.datum));
            }
            Q.p0 = p;
        }
        Q.pa = Q.p0; Q.We = 0; Q.rate = 0; Q.WeS = 0; Q.rateS = 0;
        log_(fmt("  aquifer %g: %g connections, initial pressure %.6g at datum %.6g", Q.id, double(Q.cell.size()), Q.p0, Q.datum));
    }
}

std::vector<AquiferState> Simulator::aquiferStates() const {
    std::vector<AquiferState> out;
    for (const AquiferRT& Q : aquifers_) {
        AquiferState s; s.id = Q.id; s.fetkovich = Q.type == AquiferType::Fetkovich; s.cells = int(Q.cell.size());
        s.pressure = s.fetkovich ? Q.pa : Q.p0; s.rate = Q.rate; s.cumulative = Q.We;
        out.push_back(s);
    }
    return out;
}

void Simulator::buildPattern() {
    const int na = grid_.numActive(), nw = int(wells_.size()), nv = nv_;
    blkStart_.assign(na + 1, 0); blkCells_.clear();
    adjPos_.assign(grid_.adjNbr.size(), 0); selfPos_.assign(na, 0);
    std::vector<int> tmp;
    for (int a = 0; a < na; ++a) {
        tmp.clear(); tmp.push_back(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) tmp.push_back(grid_.adjNbr[e]);
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        blkStart_[a] = int(blkCells_.size());
        blkCells_.insert(blkCells_.end(), tmp.begin(), tmp.end());
        selfPos_[a] = int(std::lower_bound(tmp.begin(), tmp.end(), a) - tmp.begin());
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e)
            adjPos_[e] = int(std::lower_bound(tmp.begin(), tmp.end(), grid_.adjNbr[e]) - tmp.begin());
    }
    blkStart_[na] = int(blkCells_.size());

    const int N = na * nv + nw;
    J_ = CSR{}; J_.n = N; J_.rowPtr.assign(N + 1, 0);
    for (int a = 0; a < na; ++a) {
        const int len = (blkStart_[a + 1] - blkStart_[a]) * nv + int(cellWells_[a].size());
        for (int e = 0; e < nv; ++e) J_.rowPtr[a * nv + e + 1] = len;
    }
    for (int w = 0; w < nw; ++w) J_.rowPtr[na * nv + w + 1] = int(wells_[w].rowCells.size()) * nv + 1;
    for (int r = 0; r < N; ++r) J_.rowPtr[r + 1] += J_.rowPtr[r];
    J_.col.resize(J_.rowPtr[N]); J_.val.assign(J_.rowPtr[N], 0.0);
    for (int a = 0; a < na; ++a)
        for (int e = 0; e < nv; ++e) {
            int p = J_.rowPtr[a * nv + e];
            for (int b = blkStart_[a]; b < blkStart_[a + 1]; ++b)
                for (int v = 0; v < nv; ++v) J_.col[p++] = blkCells_[b] * nv + v;
            for (int w : cellWells_[a]) J_.col[p++] = na * nv + w;
        }
    for (int w = 0; w < nw; ++w) {
        int p = J_.rowPtr[na * nv + w];
        for (int a : wells_[w].rowCells) for (int v = 0; v < nv; ++v) J_.col[p++] = a * nv + v;
        J_.col[p++] = na * nv + w;
    }
    J_.finalizeDiag();
    R_.assign(N, 0.0);
}

// ======================================================================== controls
static bool matchName(const std::string& pattern, const std::string& name) {
    if (!pattern.empty() && pattern.back() == '*')
        return name.compare(0, pattern.size() - 1, pattern, 0, pattern.size() - 1) == 0;
    return pattern == name;
}

void Simulator::updateConnections(int step) {
    for (auto& w : wells_) {
        w.nActive = 0;
        for (size_t c = 0; c < w.cell.size(); ++c) {
            const bool on = step >= w.cFrom[c] && (w.cUntil[c] < 0 || step < w.cUntil[c]);
            w.wi[c] = on ? w.wiBase[c] : 0.0;
            w.nActive += on ? 1 : 0;
        }
    }
}

void Simulator::applyControls(int step) {
    // connections opening / closing at this step: wells without open connections are shut, wells whose
    // first connections open now start producing if their control is OPEN
    std::vector<int> before(wells_.size());
    for (size_t w = 0; w < wells_.size(); ++w) before[w] = wells_[w].nActive;
    updateConnections(step);
    for (size_t wi = 0; wi < wells_.size(); ++wi) {
        WellRT& w = wells_[wi];
        if (w.nActive == before[wi] && step > 0) continue;
        const bool want = w.ctrl.status == WellStatus::Open && w.nActive > 0;
        if (want == w.open) continue;
        if (want) {
            int a0 = w.cell.front();
            for (size_t c = 0; c < w.cell.size(); ++c) if (w.wi[c] > 0) { a0 = w.cell[c]; break; }
            const double pc = p_[a0];
            w.bhp = (w.type == WellType::Producer) ? std::max(pc * 0.9, U_.atm) : pc * 1.05;
            if (w.mode == Control::BHP && w.ctrl.bhp > 0) w.bhp = w.ctrl.bhp;
            log_(fmt("Report step %g: well ", step + 1) + case_.wells[size_t(w.input)].name + " opens (connections added)");
        } else if (w.open) log_(fmt("Report step %g: well ", step + 1) + case_.wells[size_t(w.input)].name + " shut - no open connections");
        w.open = want;
    }
    for (const WellControl& wc : case_.schedule[step].controls) {
        bool found = false;
        for (auto& w : wells_) {
            const WellInput& in = case_.wells[w.input];
            if (!matchName(wc.well, in.name)) continue;
            found = true;
            const bool wasOpen = w.open;
            w.ctrl = wc; w.ctrl.well = in.name;
            w.open = wc.status == WellStatus::Open && w.nActive > 0;
            w.mode = wc.mode;
            if (w.type != WellType::Producer && wc.phase) w.type = wc.phase == 2 ? WellType::GasInjector : WellType::WaterInjector;
            if (w.open && !wasOpen) {
                int a0 = w.cell.front();
                for (size_t c = 0; c < w.cell.size(); ++c) if (w.wi[c] > 0) { a0 = w.cell[c]; break; }
                const double pc = p_[a0];
                w.bhp = (w.type == WellType::Producer) ? std::max(pc * 0.9, U_.atm) : pc * 1.05;
            }
            if (w.mode == Control::BHP && wc.bhp > 0) w.bhp = wc.bhp;
            if (wc.history && wc.mode == Control::RESV && wc.resv < 0 && w.type == WellType::Producer && !w.cell.empty()) {
                // WCONHIST RESV: reservoir-voidage target equivalent to the observed surface rates (ECLIPSE practice),
                // evaluated at the mean pressure of the connected cells
                double pm = 0; for (int a : w.cell) pm += p_[a]; pm /= double(w.cell.size());
                const FluidModel& F = fl(w.cell.front());
                const double qo = std::max(0.0, wc.orat), qw = std::max(0.0, wc.wrat), qg = std::max(0.0, wc.grat);
                double rs = 0, rv = 0;
                if (live_) rs = qo > 0 ? std::min(F.rsSat(pm), qg / qo) : 0.0;
                if (vap_) rv = qg > 0 ? std::min(F.rvSat(pm), qo / qg) : 0.0;
                // split the surface rates into reservoir oil / free gas: O + rv G = qo, G + rs O = qg
                const double den = std::max(1e-12, 1.0 - rs * rv);
                const double O = std::max(0.0, (qo - rv * qg) / den), Gs = std::max(0.0, (qg - rs * qo) / den);
                double resv = 0;
                if (O_) resv += O / F.boS(pm, rs);
                if (W_) resv += qw / F.bwS(pm);
                if (G_) resv += Gs / (vap_ ? F.bg(AD3(pm), AD3(rv)).v : F.bgS(pm));
                w.ctrl.resv = resv;
            }
            w.rates.oilH = wc.history ? std::max(0.0, wc.orat) : -1;
            w.rates.waterH = wc.history ? std::max(0.0, wc.wrat) : -1;
            w.rates.gasH = wc.history ? std::max(0.0, wc.grat) : -1;
        }
        if (!found) log_("Warning: control for unknown well '" + wc.well + "'");
    }
}

double Simulator::wellRateOf(const WellRT& w, Control m) const {
    switch (m) {
        case Control::ORAT: return w.rates.oil;
        case Control::WRAT: return w.rates.water;
        case Control::GRAT: return w.rates.gas;
        case Control::LRAT: return w.rates.oil + w.rates.water;
        case Control::RESV: return w.rates.resv;
        case Control::RATE: return w.type == WellType::GasInjector ? w.rates.gasInj : w.rates.waterInj;
        default: return 0;
    }
}

static double targetOf(const WellControl& c, Control m) {
    switch (m) {
        case Control::ORAT: return c.orat; case Control::WRAT: return c.wrat; case Control::GRAT: return c.grat;
        case Control::LRAT: return c.lrat; case Control::RESV: return c.resv; case Control::RATE: return c.rate;
        default: return c.bhp;
    }
}

bool Simulator::checkWellSwitching() {
    bool sw = false;
    for (auto& w : wells_) {
        if (!w.open) continue;
        const bool prod = w.type == WellType::Producer;
        const double lim = w.ctrl.bhp > 0 ? w.ctrl.bhp : (prod ? U_.atm : 1e12);
        if (w.mode != Control::BHP) {
            if ((prod && w.bhp < lim) || (!prod && w.bhp > lim)) {
                w.mode = Control::BHP; w.bhp = lim; w.ctrl.bhp = lim; sw = true;
                log_("  well " + case_.wells[w.input].name + " switched to BHP control");
            }
        } else {
            Control best = Control::BHP; double worst = 1.0 + 1e-6;
            const Control cand[] = {Control::ORAT, Control::WRAT, Control::GRAT, Control::LRAT, Control::RESV, Control::RATE};
            for (Control m : cand) {
                if (prod == (m == Control::RATE)) continue;
                const double t = targetOf(w.ctrl, m);
                if (t <= 0) continue;
                const double ratio = wellRateOf(w, m) / t;
                if (ratio > worst) { worst = ratio; best = m; }
            }
            if (best != Control::BHP) {
                w.mode = best; sw = true;
                log_(std::string("  well ") + case_.wells[w.input].name + " switched to " + toString(best) + " control");
            }
        }
    }
    return sw;
}

// ======================================================================== end-point scaling
void Simulator::setupEndPointScaling(const Case& c) {
    eps_.clear();
    const GridInput& g = c.grid;
    if (!g.hasEndPointScaling()) return;
    const int na = grid_.numActive();
    eps_.resize(size_t(na));
    int nClamped = 0;
    for (int a = 0; a < na; ++a) {
        const size_t q = size_t(grid_.globalOf[size_t(a)]);
        SatFunctions::EndPoints e = st(a).tableEndPoints();
        auto get = [&](const std::vector<double>& arr, double& v) {
            if (q < arr.size() && std::isfinite(arr[q]) && arr[q] >= 0) v = std::min(1.0, arr[q]);
        };
        get(g.swl, e.swl); get(g.swcr, e.swcr); get(g.swu, e.swu);
        get(g.sgl, e.sgl); get(g.sgcr, e.sgcr); get(g.sgu, e.sgu);
        get(g.sowcr, e.sowcr); get(g.sogcr, e.sogcr);
        // physical consistency: SWL <= SWCR < SWU, SGL <= SGCR < SGU, critical oil below the mobile range
        const SatFunctions::EndPoints o = e;
        e.swcr = std::max(e.swcr, e.swl);
        e.swu = std::max(e.swu, e.swcr + 1e-6);
        e.sgcr = std::max(e.sgcr, e.sgl);
        e.sgu = std::max(e.sgu, e.sgcr + 1e-6);
        e.sowcr = std::min(e.sowcr, std::max(0.0, 1.0 - e.swl - 1e-6));
        e.sogcr = std::min(e.sogcr, std::max(0.0, 1.0 - e.swl - e.sgl - 1e-6));
        if (o.swcr != e.swcr || o.swu != e.swu || o.sgcr != e.sgcr || o.sgu != e.sgu || o.sowcr != e.sowcr || o.sogcr != e.sogcr) ++nClamped;
        eps_[size_t(a)] = e;
    }
    std::string used;
    const std::pair<const char*, const std::vector<double>*> arrs[] = {{"SWL", &g.swl}, {"SWCR", &g.swcr}, {"SWU", &g.swu}, {"SGL", &g.sgl},
                                                                        {"SGCR", &g.sgcr}, {"SGU", &g.sgu}, {"SOWCR", &g.sowcr}, {"SOGCR", &g.sogcr}};
    for (auto& p : arrs) if (!p.second->empty()) used += std::string(used.empty() ? "" : " ") + p.first;
    log_("End-point scaling (two-point): " + used + (nClamped ? fmt(" - %g cells adjusted to consistent end points", double(nClamped)) : std::string()));
}

// ======================================================================== cell properties
Simulator::CellFull Simulator::evalCell(int c) const {
    CellFull f;
    const FluidModel& fluid_ = fl(c);
    const SatFunctions& sat_ = st(c);
    const SatFunctions::EndPoints* E = ep(c);
    const DensitySC& rhoSC = fluid_.rhoSC();
    const AD3 P = AD3::variable(p_[c], 0);
    AD3 Sw(W_ && !O_ ? 1.0 : 0.0), Sg(0.0), rs(0.0), rv(0.0);
    if (iSw_ >= 0) Sw = AD3::variable(sw_[c], iSw_);
    if (G_) {
        const int state = hasGas_[c];
        if (state == 1) {
            Sg = AD3::variable(x_[c], iX_);
            if (live_) { rs = fluid_.rsSat(P); if (!rsLim_.empty() && rs.v > rsLim_[size_t(c)]) rs = AD3(rsLim_[size_t(c)]); }
            if (vap_) rv = fluid_.rvSat(P);
        }
        else if (state == 0) { rs = AD3::variable(x_[c], iX_); if (vap_) rv = fluid_.rvSat(P); }
        else { Sg = 1.0 - Sw; rv = AD3::variable(x_[c], iX_); if (live_) rs = fluid_.rsSat(P); }
    }
    const AD3 So = (O_ && !(G_ && hasGas_[c] == 2)) ? (1.0 - Sw - Sg) : AD3(0.0);
    const AD3 pv = fluid_.pvMult(P) * grid_.poreVolRef[c];
    const double g = U_.gravity;
    (void)g;
    if (W_) {
        f.P[PW] = (O_ ? P - sat_.pcow(Sw, E) * (pcScale_.empty() ? 1.0 : pcScale_[size_t(c)]) : P);
        f.b[PW] = fluid_.bw(f.P[PW]);
        const AD3 kr = O_ ? sat_.krw(Sw, E) : AD3(1.0);
        f.mob[PW] = fluid_.bwMu(f.P[PW]) * kr;
        if (poly_ && !wmul_.empty()) f.mob[PW] *= wmul_[size_t(c)];   // polymer viscosity correction / RRF
        f.rho[PW] = f.b[PW] * rhoSC.water;
        f.acc[PW] = pv * f.b[PW] * Sw;
        f.S[PW] = Sw;
    }
    if (O_) {
        f.P[PO] = P;
        f.b[PO] = fluid_.bo(P, rs);
        const AD3 kr = (W_ || G_) ? sat_.kro(Sw, Sg, E) : AD3(1.0);
        f.mob[PO] = fluid_.boMu(P, rs) * kr;
        f.rho[PO] = f.b[PO] * (rs * (rhoSC.gas * fluid_.gasFactor()) + rhoSC.oil);
        f.acc[PO] = pv * f.b[PO] * So;
        f.mobRs = f.mob[PO] * rs;
        f.S[PO] = So;
    }
    if (G_) {
        f.P[PG] = P + sat_.pcog(Sg, E);
        if (vap_) {
            f.b[PG] = fluid_.bg(f.P[PG], rv);
            f.mob[PG] = fluid_.bgMu(f.P[PG], rv) * sat_.krg(Sg, E);
            f.rho[PG] = f.b[PG] * (rhoSC.gas * fluid_.gasFactor() + rv * rhoSC.oil);
            f.mobRv = f.mob[PG] * rv;
            if (O_) f.acc[PO] += pv * rv * f.b[PG] * Sg;
        } else {
            f.b[PG] = fluid_.bg(f.P[PG]);
            f.mob[PG] = fluid_.bgMu(f.P[PG]) * sat_.krg(Sg, E);
            f.rho[PG] = f.b[PG] * (rhoSC.gas * fluid_.gasFactor());
        }
        f.acc[PG] = pv * (f.b[PG] * Sg + (O_ ? rs * f.b[PO] * So : AD3(0.0)));
        f.S[PG] = Sg;
    }
    f.rs = rs; f.rv = rv;
    return f;
}

void Simulator::evalAll() {
    const int na = grid_.numActive();
    if (useParallel(na > 2000)) {
        #pragma omp parallel for schedule(static)
        for (int c = 0; c < na; ++c) cf_[c] = static_cast<const CellFlux&>(evalCell(c));
    } else {
        for (int c = 0; c < na; ++c) cf_[c] = static_cast<const CellFlux&>(evalCell(c));
    }
}

// ======================================================================== assembly
void Simulator::assemble(double dt) {
    const int na = grid_.numActive(), nv = nv_;
    std::fill(J_.val.begin(), J_.val.end(), 0.0);
    std::fill(R_.begin(), R_.end(), 0.0);
    const double g = U_.gravity;
    const double invDt = 1.0 / dt;
    bool ph[3] = {W_, O_, G_};

    if (useParallel(na > 2000)) {
        #pragma omp parallel for schedule(dynamic, 256)
        for (int a = 0; a < na; ++a) {
            const CellFlux& A = cf_[a];
            const int self = selfPos_[a];
            auto addSelf = [&](int e, const AD3& v, double scale) {
                const int r = a * nv + e;
                R_[r] += v.v * scale;
                double* row = &J_.val[J_.rowPtr[r] + self * nv];
                for (int k = 0; k < nv; ++k) row[k] += v.d[k] * scale;
            };
            for (int p = 0; p < 3; ++p) if (ph[p]) {
                AD3 acc = A.acc[p]; acc -= accOld_[size_t(a) * 3 + p];
                addSelf(eq_[p], acc, invDt);
            }
            for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) {
                const Connection& cn = grid_.conns[grid_.adjConn[e]];
                const int b = grid_.adjNbr[e];
                const double dD = (cn.a == a) ? cn.dDepth : -cn.dDepth;
                const CellFlux& B = cf_[b];
                const int nbp = adjPos_[e];
                auto add = [&](int eqi, const AD6& F) {
                    const int r = a * nv + eqi;
                    R_[r] += F.v;
                    double* rs = &J_.val[J_.rowPtr[r] + self * nv];
                    double* rn = &J_.val[J_.rowPtr[r] + nbp * nv];
                    for (int k = 0; k < nv; ++k) { rs[k] += F.d[k]; rn[k] += F.d[3 + k]; }
                };
                for (int p = 0; p < 3; ++p) if (ph[p]) {
                    const AD6 rho = (lift<6>(A.rho[p], 0) + lift<6>(B.rho[p], 3)) * 0.5;
                    AD6 dphi = lift<6>(A.P[p], 0) - lift<6>(B.P[p], 3) + rho * (g * dD);
                    if (!thp_.empty()) {
                        const double t = thp_[size_t(grid_.adjConn[e])];
                        if (t > 0) {
                            if (dphi.v > t) dphi -= t; else if (dphi.v < -t) dphi += t; else dphi = AD6(0.0);
                        }
                    }
                    const bool upA = dphi.v >= 0;
                    const AD6 mob = upA ? lift<6>(A.mob[p], 0) : lift<6>(B.mob[p], 3);
                    const double tr = (p == PW && !pcZw_.empty()) ? cn.trans / pcZw_[size_t(grid_.adjConn[e])] : cn.trans;
                    const AD6 F = mob * dphi * tr;
                    add(eq_[p], F);
                    if (p == PO && live_) {
                        const AD6 mr = upA ? lift<6>(A.mobRs, 0) : lift<6>(B.mobRs, 3);
                        add(eq_[PG], mr * dphi * cn.trans);
                    }
                    if (p == PG && vap_) {
                        const AD6 mv = upA ? lift<6>(A.mobRv, 0) : lift<6>(B.mobRv, 3);
                        add(eq_[PO], mv * dphi * cn.trans);
                    }
                }
            }
        }
    } else {
        for (int a = 0; a < na; ++a) {
            const CellFlux& A = cf_[a];
            const int self = selfPos_[a];
            auto addSelf = [&](int e, const AD3& v, double scale) {
                const int r = a * nv + e;
                R_[r] += v.v * scale;
                double* row = &J_.val[J_.rowPtr[r] + self * nv];
                for (int k = 0; k < nv; ++k) row[k] += v.d[k] * scale;
            };
            for (int p = 0; p < 3; ++p) if (ph[p]) {
                AD3 acc = A.acc[p]; acc -= accOld_[size_t(a) * 3 + p];
                addSelf(eq_[p], acc, invDt);
            }
            for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) {
                const Connection& cn = grid_.conns[grid_.adjConn[e]];
                const int b = grid_.adjNbr[e];
                const double dD = (cn.a == a) ? cn.dDepth : -cn.dDepth;
                const CellFlux& B = cf_[b];
                const int nbp = adjPos_[e];
                auto add = [&](int eqi, const AD6& F) {
                    const int r = a * nv + eqi;
                    R_[r] += F.v;
                    double* rs = &J_.val[J_.rowPtr[r] + self * nv];
                    double* rn = &J_.val[J_.rowPtr[r] + nbp * nv];
                    for (int k = 0; k < nv; ++k) { rs[k] += F.d[k]; rn[k] += F.d[3 + k]; }
                };
                for (int p = 0; p < 3; ++p) if (ph[p]) {
                    const AD6 rho = (lift<6>(A.rho[p], 0) + lift<6>(B.rho[p], 3)) * 0.5;
                    AD6 dphi = lift<6>(A.P[p], 0) - lift<6>(B.P[p], 3) + rho * (g * dD);
                    if (!thp_.empty()) {
                        const double t = thp_[size_t(grid_.adjConn[e])];
                        if (t > 0) {
                            if (dphi.v > t) dphi -= t; else if (dphi.v < -t) dphi += t; else dphi = AD6(0.0);
                        }
                    }
                    const bool upA = dphi.v >= 0;
                    const AD6 mob = upA ? lift<6>(A.mob[p], 0) : lift<6>(B.mob[p], 3);
                    const double tr = (p == PW && !pcZw_.empty()) ? cn.trans / pcZw_[size_t(grid_.adjConn[e])] : cn.trans;
                    const AD6 F = mob * dphi * tr;
                    add(eq_[p], F);
                    if (p == PO && live_) {
                        const AD6 mr = upA ? lift<6>(A.mobRs, 0) : lift<6>(B.mobRs, 3);
                        add(eq_[PG], mr * dphi * cn.trans);
                    }
                    if (p == PG && vap_) {
                        const AD6 mv = upA ? lift<6>(A.mobRv, 0) : lift<6>(B.mobRv, 3);
                        add(eq_[PO], mv * dphi * cn.trans);
                    }
                }
            }
        }
    }

    // ---- wells (serial; few connections)
    const int refP = O_ ? PO : PW;
    for (int w = 0; w < int(wells_.size()); ++w) {
        WellRT& W = wells_[w];
        const int wr = na * nv + w;
        const int wbase = J_.rowPtr[wr];
        const int bhpPos = J_.rowPtr[wr + 1] - 1;
        WellRates rt; rt.bhp = W.bhp; rt.mode = W.mode; rt.open = W.open;
        rt.oilH = W.rates.oilH; rt.waterH = W.rates.waterH; rt.gasH = W.rates.gasH;
        const double target = W.open ? targetOf(W.ctrl, W.mode) : 0.0;
        W.qwConn.assign(W.cell.size(), 0.0);
        // a rate target of zero stops the well (ECLIPSE behaviour) - no cross-flow through a zero-rate well
        if (!W.open || (W.mode != Control::BHP && target == 0.0)) { J_.val[bhpPos] = 1.0; R_[wr] = 0.0; W.rates = rt; continue; }
        const bool prod = W.type == WellType::Producer;
        const double scale = (W.mode == Control::BHP) ? 1.0 : 1.0 / std::max(std::fabs(target), 1.0);
        double qsel = 0;
        for (size_t k = 0; k < W.cell.size(); ++k) {
            const int a = W.cell[k];
            const CellFull C = evalCell(a);
            const AD4 bhp = AD4::variable(W.bhp, 3);
            const AD4 pconn = bhp + (W.pConn(k, g) - W.bhp);
            const AD4 Pc = lift<4>(C.P[refP], 0);
            const double wi = W.wi[k];
            const int nb = blkStart_[a + 1] - blkStart_[a];
            auto addCell = [&](int eqi, const AD4& q) {      // q: outflow (+) from the cell
                const int r = a * nv + eqi;
                R_[r] += q.v;
                double* row = &J_.val[J_.rowPtr[r] + selfPos_[a] * nv];
                for (int v = 0; v < nv; ++v) row[v] += q.d[v];
                J_.val[J_.rowPtr[r] + nb * nv + W.cellWellIdx[k]] += q.d[3];
            };
            AD4 sel(0.0);
            if (prod) {
                const AD4 dd = Pc - pconn;
                AD4 q[3]; AD4 qd(0.0), qv(0.0);   // qd: dissolved gas (in oil), qv: vaporised oil (in gas)
                for (int p = 0; p < 3; ++p) if (ph[p]) q[p] = lift<4>(C.mob[p], 0) * dd * wi;
                if (ph[PW] && k < W.wmulConn.size()) q[PW] *= W.wmulConn[k];
                if (live_) qd = lift<4>(C.mobRs, 0) * dd * wi;
                if (vap_) qv = lift<4>(C.mobRv, 0) * dd * wi;
                if (dd.v > 0) {
                    for (int p = 0; p < 3; ++p) if (ph[p]) addCell(eq_[p], p == PG ? q[p] + qd : p == PO ? q[p] + qv : q[p]);
                    rt.oil += (ph[PO] ? q[PO].v : 0) + qv.v; rt.water += ph[PW] ? q[PW].v : 0;
                    if (ph[PW]) W.qwConn[k] = q[PW].v;
                    rt.gas += (ph[PG] ? q[PG].v : 0) + qd.v;
                    double resv = 0;
                    if (ph[PO]) resv += q[PO].v / C.b[PO].v;
                    if (ph[PW]) resv += q[PW].v / C.b[PW].v;
                    if (ph[PG]) resv += q[PG].v / C.b[PG].v;
                    rt.resv += resv;
                }
                switch (W.mode) {
                    case Control::ORAT: sel = q[PO] + qv; break;
                    case Control::WRAT: sel = q[PW]; break;
                    case Control::GRAT: sel = q[PG] + qd; break;
                    case Control::LRAT: sel = q[PO] + qv + q[PW]; break;
                    case Control::RESV: {
                        if (ph[PO]) sel += q[PO] / lift<4>(C.b[PO], 0);
                        if (ph[PW]) sel += q[PW] / lift<4>(C.b[PW], 0);
                        if (ph[PG]) sel += q[PG] / lift<4>(C.b[PG], 0);
                    } break;
                    default: break;
                }
            } else {
                const AD4 dd = pconn - Pc;
                AD4 lt(0.0);
                for (int p = 0; p < 3; ++p) if (ph[p]) lt += lift<4>(C.mob[p], 0) / lift<4>(C.b[p], 0) * (p == PW && k < W.wmulConn.size() ? W.wmulConn[k] : 1.0);
                const int ip = (W.type == WellType::GasInjector) ? PG : PW;
                if (!ph[ip]) throw std::runtime_error("Injector phase is not active in the model");
                const AD4 q = lt * lift<4>(C.b[ip], 0) * dd * wi;
                if (dd.v > 0) {
                    addCell(eq_[ip], -q);
                    if (ip == PW) { rt.waterInj += q.v; W.qwConn[k] = -q.v; } else rt.gasInj += q.v;
                    rt.resv -= q.v / C.b[ip].v;
                }
                sel = q;
            }
            if (W.mode != Control::BHP) {
                qsel += sel.v;
                const int cp = W.cpos[k];
                for (int v = 0; v < nv; ++v) J_.val[wbase + cp * nv + v] += sel.d[v] * scale;
                J_.val[bhpPos] += sel.d[3] * scale;
            }
        }
        if (W.mode == Control::BHP) { R_[wr] = W.bhp - target; J_.val[bhpPos] = 1.0; }
        else {
            R_[wr] = (qsel - target) * scale;
            if (std::fabs(J_.val[bhpPos]) < 1e-14) J_.val[bhpPos] = prod ? -1e-8 : 1e-8;
        }
        W.rates = rt;
    }

    // ---- analytic aquifers: water source in the connected cells (reservoir volume rate * bw)
    for (AquiferRT& Q : aquifers_) {
        double K1 = 0, K0 = 0, base = Q.p0;
        if (Q.type == AquiferType::Fetkovich) {
            K1 = Q.ctV / dt * (1.0 - std::exp(-Q.J * dt / Q.ctV));
            base = Q.pa;
        } else {
            const double tn = time_ / Q.Tc, tn1 = (time_ + dt) / Q.Tc;
            double pd, dpd; carterTracyPD(tn1, pd, dpd);
            const double den = std::max(pd - tn * dpd, 1e-12);
            K1 = Q.beta * (tn1 - tn) / (dt * den);
            K0 = Q.We * dpd * (tn1 - tn) / (dt * den);
        }
        Q.rate = 0; Q.rateS = 0;
        for (size_t k = 0; k < Q.cell.size(); ++k) {
            const int a = Q.cell[k];
            const CellFull C = evalCell(a);
            const double head = base + C.rho[PW].v * g * (Q.cdepth[k] - Q.datum);
            const AD3 q = (head - C.P[PW]) * (Q.alpha[k] * K1) - Q.alpha[k] * K0;   // into the cell (res. volume)
            // inflow carries aquifer water; outflow is upwinded with the cell's water mobility (krw * b / mu)
            const AD3 src = q.v >= 0 ? C.b[PW] * q : C.mob[PW] * q * Q.muw;
            Q.rate += q.v >= 0 ? q.v : src.v / C.b[PW].v;
            Q.rateS += src.v;
            const int r = a * nv + eq_[PW];
            R_[r] -= src.v;
            double* row = &J_.val[J_.rowPtr[r] + selfPos_[a] * nv];
            for (int v = 0; v < nv; ++v) row[v] -= src.d[v];
        }
    }

    // ---- row scaling of cell equations: R * dt / PVref (dimensionless, saturation-like)
    if (useParallel(na > 2000)) {
        #pragma omp parallel for schedule(static)
        for (int a = 0; a < na; ++a) {
            const double s = dt / grid_.poreVolRef[a];
            for (int e = 0; e < nv; ++e) {
                const int r = a * nv + e;
                R_[r] *= s;
                for (int q = J_.rowPtr[r]; q < J_.rowPtr[r + 1]; ++q) J_.val[q] *= s;
            }
        }
    } else {
        for (int a = 0; a < na; ++a) {
            const double s = dt / grid_.poreVolRef[a];
            for (int e = 0; e < nv; ++e) {
                const int r = a * nv + e;
                R_[r] *= s;
                for (int q = J_.rowPtr[r]; q < J_.rowPtr[r + 1]; ++q) J_.val[q] *= s;
            }
        }
    }
}

bool Simulator::converged(double dt, double& cnvMax, double& mbMax) const {
    (void)dt;
    const int na = grid_.numActive(), nv = nv_;
    double cnv[3] = {0, 0, 0}, sum[3] = {0, 0, 0}, pvb[3] = {0, 0, 0};
    for (int a = 0; a < na; ++a) {
        const double pv = grid_.poreVolRef[a];
        for (int p = 0; p < 3; ++p) {
            const int e = eq_[p];
            if (e < 0) continue;
            const double r = R_[a * nv + e];
            cnv[p] = std::max(cnv[p], std::fabs(r) / bRef_[p]);
            sum[p] += r * pv;
            pvb[p] += pv * bRef_[p];
        }
    }
    cnvMax = 0; mbMax = 0;
    for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) {
        cnvMax = std::max(cnvMax, cnv[p]);
        mbMax = std::max(mbMax, std::fabs(sum[p]) / pvb[p]);
    }
    bool wellsOk = true;
    for (int w = 0; w < int(wells_.size()); ++w) {
        const double r = std::fabs(R_[na * nv + w]);
        const double tol = (wells_[w].mode == Control::BHP) ? 1e-4 * U_.atm : 1e-5;
        if (wells_[w].open && r > tol) wellsOk = false;
    }
    return cnvMax < case_.tuning.cnvTol && mbMax < case_.tuning.mbTol && wellsOk;
}

// ======================================================================== Newton update
void Simulator::updateState(const std::vector<double>& dx) {
    const int na = grid_.numActive(), nv = nv_;
    const Tuning& T = case_.tuning;
    // primary variable update with Appleyard-type chopping and the three-state variable switch:
    //   state 0 (no free gas): x = Rs     state 1 (oil + gas): x = Sg     state 2 (no oil, wet gas): x = Rv
    auto update = [&](int a) {
        double dp = dx[a * nv];
        const double dpMax = T.dpMaxRel * std::max(p_[a], 10.0 * U_.atm);
        dp = std::max(-dpMax, std::min(dpMax, dp));
        const double dsw = iSw_ >= 0 ? dx[a * nv + iSw_] : 0.0;
        const double dX = iX_ >= 0 ? dx[a * nv + iX_] : 0.0;
        const int st = G_ ? hasGas_[a] : 1;
        const double dsg = !G_ ? 0.0 : st == 1 ? dX : st == 2 ? -dsw : 0.0;
        const double dso = st == 2 ? 0.0 : -(dsw + dsg);
        const double m = std::max(std::fabs(dsw), std::max(std::fabs(dsg), O_ ? std::fabs(dso) : 0.0));
        const double f = m > T.dsMax ? T.dsMax / m : 1.0;
        sw_[a] += f * dsw;
        p_[a] = std::max(p_[a] + dp, 0.05 * U_.atm);
        if (iSw_ >= 0) sw_[a] = std::max(0.0, std::min(1.0, sw_[a]));
        if (!G_) return;
        const FluidModel& F = fl(a);
        if (st == 1) x_[a] += f * dX;
        else if (st == 0) {
            const double lim = 0.25 * std::max(F.rsSat(p_[a]), std::max(x_[a], 1e-6));
            x_[a] += std::max(-lim, std::min(lim, dX));
        } else {
            const double lim = 0.25 * std::max(F.rvSat(p_[a]), std::max(x_[a], 1e-9));
            x_[a] += std::max(-lim, std::min(lim, dX));
        }
        if (!O_) { if (x_[a] < 0) x_[a] = 0; return; }
        const double sw = iSw_ >= 0 ? sw_[a] : 0.0;
        // without oil and gas (water only) the Rs / Rv variables are undefined: such cells use state 1 with Sg = 0
        const bool waterOnly = iSw_ >= 0 && sw > 1.0 - 1e-9;
        if (st == 1) {
            if (x_[a] < 0) {
                if (live_ && !waterOnly) { hasGas_[a] = 0; x_[a] = F.rsSat(p_[a]); if (!rsLim_.empty()) x_[a] = std::min(x_[a], rsLim_[size_t(a)]); }
                else x_[a] = 0;
            } else if (sw + x_[a] > 1.0) {
                if (vap_ && !waterOnly) { hasGas_[a] = 2; x_[a] = F.rvSat(p_[a]); }   // oil vanished: gas carries it
                else { const double s = sw + x_[a]; if (iSw_ >= 0) sw_[a] /= s; x_[a] /= s; }
            }
        } else if (st == 0) {
            if (x_[a] < 0) x_[a] = 0;
            const double rsCap = rsLim_.empty() ? F.rsSat(p_[a]) : std::min(F.rsSat(p_[a]), rsLim_[size_t(a)] * (1 + 1e-9) + 1e-12);
            if (waterOnly || x_[a] > rsCap) { hasGas_[a] = 1; x_[a] = 0.0; }
        } else {
            if (x_[a] < 0) x_[a] = 0;
            if (waterOnly) { hasGas_[a] = 1; x_[a] = 0.0; }
            else if (x_[a] > F.rvSat(p_[a])) { hasGas_[a] = 1; x_[a] = std::max(0.0, 1.0 - sw); }   // condensate drops out
        }
    };
    if (useParallel(na > 4000)) {
        #pragma omp parallel for schedule(static)
        for (int a = 0; a < na; ++a) update(a);
    } else {
        for (int a = 0; a < na; ++a) update(a);
    }
    for (int w = 0; w < int(wells_.size()); ++w) {
        if (!wells_[w].open) continue;
        double d = dx[na * nv + w];
        const double lim = T.dpMaxRel * std::max(wells_[w].bhp, 10.0 * U_.atm);
        d = std::max(-lim, std::min(lim, d));
        wells_[w].bhp = std::max(wells_[w].bhp + d, 0.05 * U_.atm);
    }
}

bool Simulator::newton(double dt, int& its, int& lin) {
    const Tuning& T = case_.tuning;
    std::vector<double> rhs(J_.n), dx;
    lin = 0;
    std::vector<double> prevWellR(wells_.size(), 0.0);
    for (its = 0; its <= T.maxNewton; ++its) {
        evalAll();
        assemble(dt);
        double cnv, mb;
        const bool conv = converged(dt, cnv, mb);
        {   // a rate target the well cannot deliver (its residual does not move) -> BHP control at the limit
            const int na = grid_.numActive();
            for (int w = 0; w < int(wells_.size()); ++w) {
                WellRT& W = wells_[size_t(w)];
                const double r = std::fabs(R_[size_t(na) * nv_ + w]);
                if (its >= 3 && W.open && W.mode != Control::BHP && r > 0.5 && r > 0.95 * prevWellR[size_t(w)]) {
                    const bool prod = W.type == WellType::Producer;
                    W.mode = Control::BHP; W.bhp = W.ctrl.bhp > 0 ? W.ctrl.bhp : (prod ? U_.atm : 1e12); W.ctrl.bhp = W.bhp;
                    log_("  well " + case_.wells[W.input].name + " cannot meet its rate target - switched to BHP control");
                }
                prevWellR[size_t(w)] = r;
            }
        }
        if (std::getenv("RESSIM_DEBUG")) {
            const int na = grid_.numActive(); int worst = 0; double wv = 0;
            for (int r = 0; r < na * nv_; ++r) if (std::fabs(R_[r]) > wv) { wv = std::fabs(R_[r]); worst = r; }
            int nsat = 0; for (int a = 0; a < na; ++a) nsat += hasGas_[a] == 1;
            std::fprintf(stderr, "    it %d cnv %.3e mb %.3e worst cell %d eq %d (p %.1f sw %.3f x %.4f gas %d) ngas %d wellR",
                         its, cnv, mb, worst / nv_, worst % nv_, p_[worst / nv_], sw_[worst / nv_], x_[worst / nv_],
                         int(hasGas_[worst / nv_]), nsat);
            for (int w = 0; w < int(wells_.size()); ++w) std::fprintf(stderr, " %.2e", R_[na * nv_ + w]);
            std::fprintf(stderr, "\n");
        }
        if (conv) return true;
        if (its == T.maxNewton) break;
        for (int i = 0; i < J_.n; ++i) rhs[i] = -R_[i];
        if (std::getenv("RESSIM_DEBUG") && its == 0) {
            int zeroDiag = 0, nonFinite = 0, zeroRow = 0; double maxd = 0, mind = 1e300;
            for (int i = 0; i < J_.n; ++i) {
                double rs = 0; bool fin = true;
                for (int k = J_.rowPtr[i]; k < J_.rowPtr[i + 1]; ++k) { rs += std::fabs(J_.val[k]); if (!std::isfinite(J_.val[k])) fin = false; }
                const double d = std::fabs(J_.val[J_.diag[i]]);
                if (!fin) { if (nonFinite < 5) std::fprintf(stderr, "      non-finite row %d\n", i); ++nonFinite; }
                if (rs == 0) { if (zeroRow < 5) std::fprintf(stderr, "      zero row %d (cell %d eq %d)\n", i, i / nv_, i % nv_); ++zeroRow; }
                else if (d < 1e-14 * rs) { if (zeroDiag < 5) std::fprintf(stderr, "      tiny diag row %d (cell %d eq %d) d=%g rowsum=%g\n", i, i / nv_, i % nv_, d, rs); ++zeroDiag; }
                if (d > 0) { maxd = std::max(maxd, d); mind = std::min(mind, d); }
            }
            std::fprintf(stderr, "    matrix n=%d zeroRows=%d tinyDiag=%d nonFinite=%d diag range %.3e..%.3e\n", J_.n, zeroRow, zeroDiag, nonFinite, mind, maxd);
        }
        const Preconditioner* M = &ilu_;
        if (useCPR_) {
            const double* acc = nullptr;
            if (case_.tuning.cprWeights == 1) {
                const int na = grid_.numActive(), nv = nv_;
                accD_.assign(size_t(na) * nv * nv, 0.0);
                const bool ph[3] = {W_, O_, G_};
                // same per-cell row scaling as the assembled system (dt / PV; the constant dt does not matter)
                for (int a = 0; a < na; ++a) {
                    const double sc = 1.0 / grid_.poreVolRef[a];
                    for (int p = 0; p < 3; ++p) if (ph[p])
                        for (int k = 0; k < nv; ++k) accD_[(size_t(a) * nv + eq_[p]) * nv + k] = cf_[a].acc[p].d[k] * sc;
                }
                acc = accD_.data();
            }
            cpr_.setup(J_, grid_.numActive(), nv_, int(wells_.size()), acc);
            M = &cpr_; stats_.amgLevels = cpr_.amg().levels();
        }
        else ilu_.factor(J_);
        LinResult lr = bicgstab(J_, *M, rhs, dx, T.linTol, T.maxLinear);
        lin += lr.iterations;
        if (!lr.converged) {
            LinResult lg = gmres(J_, *M, rhs, dx, T.linTol, T.maxLinear, 40);
            lin += lg.iterations;
            if (!lg.converged && lg.relResidual > 1e-2) { log_(fmt("  linear solver failed (res %.2e)", lg.relResidual)); return false; }
        }
        for (double v : dx) if (!std::isfinite(v)) { log_("  non-finite Newton update"); return false; }
        updateState(dx);
        // rate-controlled wells whose BHP crosses the limit switch to BHP control inside the Newton loop
        // (otherwise a target that cannot be met makes Newton stall and the step is cut forever)
        for (auto& w : wells_) {
            if (!w.open || w.mode == Control::BHP) continue;
            const bool prod = w.type == WellType::Producer;
            const double blim = w.ctrl.bhp > 0 ? w.ctrl.bhp : (prod ? U_.atm : 1e12);
            if ((prod && w.bhp < blim) || (!prod && w.bhp > blim)) {
                w.mode = Control::BHP; w.bhp = blim; w.ctrl.bhp = blim;
                log_("  well " + case_.wells[w.input].name + " switched to BHP control");
            }
        }
    }
    return false;
}

// ======================================================================== time stepping
void Simulator::computeHeads() {
    // Connection pressures of a standard (non-segmented) well, OPM Flow / ECLIPSE style, evaluated explicitly
    // at the start of each time step:
    //  1. surface component rates of every connection at the current state (old bhp and connection heads);
    //  2. walking from the bottom connection upwards, the flow leaving each connection up the wellbore is the
    //     sum of the inflows at and below it; its component mixture is split into phases with Rs/Rv capped by
    //     the saturated values at the average wellbore pressure of the segment;
    //  3. segment density = sum(component * rho_sc) / sum(phase / b), and the connection pressure differences are
    //     accumulated from the reference depth downwards.
    // Stopped / zero-rate producers use the mobility-weighted cell densities of each connection.
    const double g = U_.gravity;
    const int refP = O_ ? PO : PW;
    for (auto& w : wells_) {
        if (!w.open || w.cell.empty()) continue;
        const size_t n = w.cell.size();
        const bool prod = w.type == WellType::Producer;
        std::vector<size_t> ord(n);
        for (size_t k = 0; k < n; ++k) ord[k] = k;
        std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return w.cdepth[a] < w.cdepth[b]; });
        if (w.dpConn.size() != n) {   // first step / changed connections: start from the cell densities
            w.dpConn.assign(n, 0.0);
            double dp = 0, zAbove = w.refDepth;
            for (size_t r = 0; r < n; ++r) {
                const size_t k = ord[r];
                const CellFull C = evalCell(w.cell[k]);
                double num = 0, den = 0;
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) { const double m = C.mob[p].v / std::max(C.b[p].v, 1e-30); num += m * C.rho[p].v; den += m; }
                const double rc = den > 0 ? num / den : (O_ ? C.rho[PO].v : C.rho[PW].v);
                dp += (w.cdepth[k] - zAbove) * rc * g; w.dpConn[k] = dp; zAbove = w.cdepth[k];
            }
        }
        std::vector<double> pOld(n);
        for (size_t k = 0; k < n; ++k) pOld[k] = w.pConn(k, g);
        // 1. connection component rates (oil, water, gas surface volumes; + = out of the reservoir)
        std::vector<double> qc(3 * n, 0.0);
        std::vector<double> rhoCell(n, 0.0);
        double tot[3] = {0, 0, 0};
        for (size_t k = 0; k < n; ++k) {
            const CellFull C = evalCell(w.cell[k]);
            double num = 0, den = 0;
            for (int p = 0; p < 3; ++p) {
                if (eq_[p] < 0) continue;
                const double m = C.mob[p].v / std::max(C.b[p].v, 1e-30);
                num += m * C.rho[p].v; den += m;
            }
            rhoCell[k] = den > 0 ? num / den : (O_ ? C.rho[PO].v : C.rho[PW].v);
            if (w.wi[k] <= 0) continue;
            const double dd = C.P[refP].v - pOld[k];
            if (prod && dd > 0) {
                const double wm = (k < w.wmulConn.size()) ? w.wmulConn[k] : 1.0;
                if (eq_[PO] >= 0) qc[3 * k + 0] = C.mob[PO].v * dd * w.wi[k];
                if (eq_[PW] >= 0) qc[3 * k + 1] = C.mob[PW].v * dd * w.wi[k] * wm;
                if (eq_[PG] >= 0) qc[3 * k + 2] = C.mob[PG].v * dd * w.wi[k];
                if (live_) qc[3 * k + 2] += C.mobRs.v * dd * w.wi[k];
                if (vap_) qc[3 * k + 0] += C.mobRv.v * dd * w.wi[k];
            } else if (!prod && dd < 0) {
                const int ip = w.type == WellType::GasInjector ? PG : PW;
                double lt = 0;
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) lt += C.mob[p].v / std::max(C.b[p].v, 1e-30);
                if (eq_[ip] >= 0) qc[3 * k + (ip == PG ? 2 : 1)] = lt * C.b[ip].v * dd * w.wi[k];   // < 0: injection
            }
            for (int c = 0; c < 3; ++c) tot[c] += qc[3 * k + c];
        }
        std::vector<double> rhoSeg(n, 0.0);
        const bool flowing = (std::fabs(tot[0]) + std::fabs(tot[1]) + std::fabs(tot[2])) > 0;
        if (prod && !flowing) {
            for (size_t k = 0; k < n; ++k) rhoSeg[k] = rhoCell[k];
        } else {
            // 2. outflow up the wellbore at each connection, bottom to top
            std::vector<double> qout(3 * n, 0.0);
            double below[3] = {0, 0, 0};
            for (size_t r = n; r-- > 0;) {
                const size_t k = ord[r];
                for (int c = 0; c < 3; ++c) { below[c] += qc[3 * k + c]; qout[3 * k + c] = below[c]; }
            }
            const double aO = std::fabs(tot[0]), aG = std::fabs(tot[2]);
            double prevMix[3] = {0, 0, 0};
            for (size_t r = 0; r < n; ++r) {
                const size_t k = ord[r];
                const FluidModel& F = fl(w.cell[k]);
                const DensitySC& rs0 = F.rhoSC();
                const double gf = F.gasFactor();
                const double pAbove = (r == 0) ? w.bhp : pOld[ord[r - 1]];
                const double pAvg = 0.5 * (pOld[k] + pAbove);
                double cm[3] = {0, 0, 0};
                const double sum = std::fabs(qout[3 * k]) + std::fabs(qout[3 * k + 1]) + std::fabs(qout[3 * k + 2]);
                if (sum > 0) for (int c = 0; c < 3; ++c) cm[c] = std::fabs(qout[3 * k + c]) / sum;
                else if (!prod) cm[w.type == WellType::GasInjector ? 2 : 1] = 1.0;
                else if (r == 0) cm[O_ ? 0 : 1] = 1.0;
                else for (int c = 0; c < 3; ++c) cm[c] = prevMix[c];
                double pm[3] = {cm[0], cm[1], cm[2]};
                const double rsMax = live_ ? F.rsSat(pAvg) : 0.0;
                const double rvMax = vap_ ? F.rvSat(pAvg) : 0.0;
                if (O_ && G_) {
                    double rs = 0, rv = 0;
                    if (live_ && cm[0] > 1e-12) rs = std::min(cm[2] / cm[0], rsMax);
                    if (vap_ && cm[2] > 1e-12) rv = std::min(cm[0] / cm[2], rvMax);
                    const double d = 1.0 - rs * rv;
                    if (d > 1e-12) {
                        if (rs > 0) pm[2] = (cm[2] - cm[0] * rs) / d;
                        if (rv > 0) pm[0] = (cm[0] - cm[2] * rv) / d;
                    }
                }
                for (int c = 0; c < 3; ++c) prevMix[c] = pm[c];
                // formation volume factors at the segment pressure, Rs/Rv from the well's surface rates
                double num = 0, vol = 0;
                if (O_) {
                    double rs = 0;
                    if (live_ && aG > 0) { rs = aO > 0 ? aG / aO : 0.0; rs = std::min(rs, rsMax); }
                    num += cm[0] * rs0.oil; vol += pm[0] / std::max(F.boS(pAvg, rs), 1e-30);
                }
                if (W_) { num += cm[1] * rs0.water; vol += pm[1] / std::max(F.bwS(pAvg), 1e-30); }
                if (G_) {
                    double bgv;
                    if (vap_) {
                        double rv = 0;
                        if (aO > 0) { rv = aG > 0 ? aO / aG : 0.0; rv = std::min(rv, rvMax); }
                        bgv = F.bg(AD3(pAvg), AD3(rv)).v;
                    } else bgv = F.bgS(pAvg);
                    num += cm[2] * rs0.gas * gf; vol += pm[2] / std::max(bgv, 1e-30);
                }
                rhoSeg[k] = vol > 0 ? num / vol : rhoCell[k];
            }
        }
        // 3. accumulate the pressure differences from the reference depth
        w.dpConn.assign(n, 0.0);
        double dp = 0, zAbove = w.refDepth;
        for (size_t r = 0; r < n; ++r) {
            const size_t k = ord[r];
            dp += (w.cdepth[k] - zAbove) * rhoSeg[k] * g;
            w.dpConn[k] = dp; zAbove = w.cdepth[k];
        }
        w.headRho = rhoSeg[ord[0]];
    }
}

void Simulator::storeOld() {
    const int na = grid_.numActive();
    for (int a = 0; a < na; ++a)
        for (int p = 0; p < 3; ++p) accOld_[size_t(a) * 3 + p] = (eq_[p] >= 0) ? cf_[a].acc[p].v : 0.0;
    // reference b per phase for convergence norms
    bRef_.assign(3, 1.0);
    double s[3] = {0, 0, 0}; int n = 0;
    for (int a = 0; a < na; a += std::max(1, na / 2000)) {
        const CellFull C = evalCell(a);
        for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) s[p] += C.b[p].v;
        ++n;
    }
    for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && n > 0) bRef_[p] = std::max(1e-12, s[p] / n);
    polymerAccept();
    drsdtAccept();
}

void Simulator::drsdtAccept() {
    if (drsdt_ < 0 || !live_ || !G_ || !O_) { lastRs_.clear(); return; }
    const int na = grid_.numActive();
    lastRs_.assign(size_t(na), 1e300);
    for (int a = 0; a < na; ++a) {
        const int st = hasGas_[a];
        if (st == 2) continue;
        const bool freeGas = st == 1 && x_[a] > 1e-7;          // OPM freeGasMinSaturation
        if (drsdtFree_ && !freeGas) continue;
        double rs = st == 0 ? x_[a] : fl(a).rsSat(p_[a]);
        if (st == 1 && !rsLim_.empty()) rs = std::min(rs, rsLim_[size_t(a)]);
        lastRs_[size_t(a)] = rs;
    }
}

void Simulator::drsdtControls(int step) {
    const ReportStep& S = case_.schedule[size_t(step)];
    if (S.drsdt < -1.5) return;
    drsdt_ = S.drsdt; drsdtFree_ = S.drsdtFree;
    if (drsdt_ >= 0) log_(fmt("  DRSDT %.4g /day", drsdt_) + (drsdtFree_ ? " (FREE)" : ""));
    drsdtAccept();
    if (drsdt_ < 0) rsLim_.clear();
}

bool Simulator::runReportStep() {
    if (finished()) return false;
    auto t0 = std::chrono::steady_clock::now();
    applyControls(step_);
    polymerControls(step_);
    drsdtControls(step_);
    const Tuning& T = case_.tuning;
    const double tEnd = time_ + case_.schedule[step_].length;
    struct WB { double bhp; Control mode; };
    std::vector<WB> wb(wells_.size());
    while (time_ < tEnd - 1e-9) {
        if (cancel && cancel->load()) { log_("Run cancelled"); return false; }
        double dt = std::min(dtNext_, T.maxDt);
        const double remaining = tEnd - time_;
        bool clipped = false;
        if (dt >= remaining) { dt = remaining; clipped = true; }
        else if (dt > 0.5 * remaining) { dt = 0.5 * remaining; clipped = true; }
        // backup
        const auto p0 = p_, s0 = sw_, x0 = x_; const auto h0 = hasGas_;
        const auto pc0 = pc_, wm0 = wmul_;
        for (size_t w = 0; w < wells_.size(); ++w) wb[w] = {wells_[w].bhp, wells_[w].mode};
        if (!lastRs_.empty()) {   // DRSDT: Rs may grow by at most DRSDT * dt in this time step
            rsLim_.resize(lastRs_.size());
            for (size_t a = 0; a < lastRs_.size(); ++a) rsLim_[a] = lastRs_[a] + drsdt_ * dt;
        } else rsLim_.clear();
        computeHeads();
        int its = 0, lin = 0, totIts = 0, totLin = 0;
        bool ok = newton(dt, its, lin);
        totIts += its; totLin += lin;
        for (int sw = 0; ok && sw < 4 && checkWellSwitching(); ++sw) {
            ok = newton(dt, its, lin); totIts += its; totLin += lin;
        }
        if (ok && poly_) {   // flow <-> polymer outer coupling iterations
            const double tol = 1e-4;
            for (int outer = 0;; ++outer) {
                const double ch = polymerTransport(dt);
                stats_.polyOuter++;
                if (ch < 0) { ok = false; log_("  polymer transport did not converge"); break; }
                if (ch <= tol) break;
                if (outer >= 9) { ok = false; log_(fmt("  polymer coupling not converged (mobility change %.2e)", ch)); break; }
                int its2 = 0, lin2 = 0;
                ok = newton(dt, its2, lin2); totIts += its2; totLin += lin2;
                its = std::max(its, its2);
                if (!ok) break;
            }
        }
        stats_.newton += totIts; stats_.linear += totLin;
        if (!ok) {
            p_ = p0; sw_ = s0; x_ = x0; hasGas_ = h0; pc_ = pc0; wmul_ = wm0;
            for (size_t w = 0; w < wells_.size(); ++w) { wells_[w].bhp = wb[w].bhp; wells_[w].mode = wb[w].mode; }
            dtNext_ = dt * T.cutFactor; stats_.cuts++;
            log_(fmt("  timestep cut: dt %.4g -> %.4g days", dt, dtNext_));
            if (dtNext_ < T.minDt) throw std::runtime_error(fmt("Time step below minimum at t = %.4f days - simulation stopped", time_));
            continue;
        }
        for (AquiferRT& Q : aquifers_) {
            Q.We += Q.rate * dt; Q.WeS += Q.rateS * dt;
            if (Q.type == AquiferType::Fetkovich) Q.pa -= Q.rate * dt / Q.ctV;
        }
        storeOld();
        time_ += dt;
        stats_.timesteps++;
        record(dt, totIts, totLin);
        double next = dt;
        if (its <= 5) next = dt * T.growth; else if (its > 9) next = dt * 0.6;
        if (clipped) next = std::max(next, dtNext_);
        dtNext_ = std::min(next, T.maxDt);
    }
    ++step_;
    snapshot();
    stats_.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

void Simulator::record(double dt, int its, int lin) {
    SummaryPoint s;
    s.time = time_; s.dt = dt; s.newton = its; s.linear = lin;
    for (auto& w : wells_) {
        WellRates r = w.rates; r.bhp = w.bhp; r.mode = w.mode; r.open = w.open;
        s.wells.push_back(r);
        s.fopr += r.oil; s.fwpr += r.water; s.fgpr += r.gas; s.fwir += r.waterInj; s.fgir += r.gasInj;
    }
    cum_[0] += s.fopr * dt; cum_[1] += s.fwpr * dt; cum_[2] += s.fgpr * dt; cum_[3] += s.fwir * dt; cum_[4] += s.fgir * dt;
    s.fopt = cum_[0]; s.fwpt = cum_[1]; s.fgpt = cum_[2]; s.fwit = cum_[3]; s.fgit = cum_[4];
    const double liq = s.fopr + s.fwpr;
    s.fwct = liq > 0 ? s.fwpr / liq : 0; s.fgor = s.fopr > 0 ? s.fgpr / s.fopr : 0;
    double ppv = 0, pv = 0, phc = 0, hc = 0;
    const int na = grid_.numActive();
    for (int a = 0; a < na; ++a) {
        ppv += p_[a] * grid_.poreVolRef[a]; pv += grid_.poreVolRef[a];
        const double h = grid_.poreVolRef[a] * (iSw_ >= 0 ? std::max(0.0, 1.0 - sw_[a]) : 1.0);
        phc += p_[a] * h; hc += h;
        if (eq_[PO] >= 0) s.oip += cf_[a].acc[PO].v;
        if (eq_[PW] >= 0) s.wip += cf_[a].acc[PW].v;
        if (eq_[PG] >= 0) s.gip += cf_[a].acc[PG].v;
    }
    // ECLIPSE definitions: FPR = hydrocarbon-pore-volume weighted (PV weighted when there is no hydrocarbon), FPRP = PV weighted
    s.fprp = ppv / pv;
    s.fpr = (eq_[PO] >= 0 || eq_[PG] >= 0) && hc > 1e-12 * pv ? phc / hc : s.fprp;
    for (const AquiferRT& Q : aquifers_) { s.faqr += Q.rate; s.faqt += Q.We; s.faqrS += Q.rateS; s.faqtS += Q.WeS; }
    if (poly_) {
        for (const WellRates& r : s.wells) { s.fcir += r.polyInj; s.fcpr += r.polyProd; }
        pcCum_[0] += s.fcir * dt; pcCum_[1] += s.fcpr * dt;
        s.fcit = pcCum_[0]; s.fcpt = pcCum_[1];
        const PolymerInput& P = case_.polymer;
        for (int a = 0; a < na; ++a) {
            const PolymerRock& R = P.plyrock[std::min(size_t(satReg_[size_t(a)]), P.plyrock.size() - 1)];
            const double pv = fl(a).pvMult(p_[a]) * grid_.poreVolRef[a];
            s.fcip += cf_[a].acc[PW].v * pc_[a] * (1.0 - R.deadPoreVolume);
            s.fcad += std::max(0.0, grid_.bulkVol[size_t(grid_.globalOf[a])] - pv) * R.rockDensity * polyEval(a, pc_[a]).ads;
        }
    }
    summary_.push_back(s);
    if (onTimestep) onTimestep(s);
}

void Simulator::snapshot() {
    const int na = grid_.numActive();
    Snapshot s; s.time = time_;
    s.pressure.resize(na); s.sw.resize(na); s.sg.resize(na); s.so.resize(na); s.rs.resize(na);
    if (vap_) s.rv.resize(na);
    for (int a = 0; a < na; ++a) {
        const double sw = iSw_ >= 0 ? sw_[a] : (W_ && !O_ ? 1.0 : 0.0);
        const int st = G_ ? hasGas_[a] : 1;
        const double sg = !G_ ? 0.0 : st == 1 ? x_[a] : st == 2 ? 1.0 - sw : 0.0;
        double rs = live_ ? (st == 0 ? x_[a] : fl(a).rsSat(p_[a])) : 0.0;
        if (live_ && st == 1 && !rsLim_.empty()) rs = std::min(rs, rsLim_[size_t(a)]);
        const double rv = vap_ ? (st == 2 ? x_[a] : fl(a).rvSat(p_[a])) : 0.0;
        s.pressure[a] = float(p_[a]); s.sw[a] = float(sw); s.sg[a] = float(sg);
        s.so[a] = float(O_ && st != 2 ? std::max(0.0, 1.0 - sw - sg) : 0.0); s.rs[a] = float(rs);
        if (vap_) s.rv[a] = float(rv);
    }
    if (poly_) {
        s.poly.resize(na); s.polyAds.resize(na);
        for (int a = 0; a < na; ++a) { s.poly[a] = float(pc_[a]); s.polyAds[a] = float(polyEval(a, pc_[a]).ads); }
    }
    snaps_.push_back(std::move(s));
}

double Simulator::totalSimTime() const {
    double t = 0; for (auto& r : case_.schedule) t += r.length; return t;
}

InPlace Simulator::inPlace() const {
    InPlace ip;
    const int na = grid_.numActive();
    for (int a = 0; a < na; ++a) {
        const CellFull C = evalCell(a);
        const double pv = fl(a).pvMult(p_[a]) * grid_.poreVolRef[a];
        if (O_) { ip.oil += pv * C.b[PO].v * C.S[PO].v; ip.gasDissolved += pv * C.b[PO].v * C.S[PO].v * C.rs.v; }
        if (W_) ip.water += pv * C.b[PW].v * C.S[PW].v;
        if (G_) ip.gasFree += pv * C.b[PG].v * C.S[PG].v;
        if (vap_) ip.oilVaporized += pv * C.b[PG].v * C.S[PG].v * C.rv.v;
    }
    return ip;
}

// ======================================================================== flux field (streamlines)
void Simulator::connectionFluxes(std::vector<double>& qw, std::vector<double>& qo, std::vector<double>& qg) const {
    const int na = grid_.numActive();
    const size_t nc = grid_.conns.size();
    qw.assign(nc, 0.0); qo.assign(nc, 0.0); qg.assign(nc, 0.0);
    std::vector<CellFull> C(na);
    if (useParallel(na > 2000)) {
        #pragma omp parallel for schedule(static)
        for (int a = 0; a < na; ++a) C[a] = evalCell(a);
    } else {
        for (int a = 0; a < na; ++a) C[a] = evalCell(a);
    }
    const double g = U_.gravity;
    if (useParallel(nc > 4000)) {
        #pragma omp parallel for schedule(static)
        for (long long ci = 0; ci < (long long)nc; ++ci) {
            const Connection& cn = grid_.conns[ci];
            const CellFull& A = C[cn.a]; const CellFull& B = C[cn.b];
            for (int p = 0; p < 3; ++p) {
                if (eq_[p] < 0) continue;
                const double rho = 0.5 * (A.rho[p].v + B.rho[p].v);
                double dphi = A.P[p].v - B.P[p].v + rho * g * cn.dDepth;
                if (!thp_.empty() && thp_[size_t(ci)] > 0) { const double t = thp_[size_t(ci)]; dphi = dphi > t ? dphi - t : dphi < -t ? dphi + t : 0.0; }
                const CellFull& U = dphi >= 0 ? A : B;
                const double tr = (p == PW && !pcZw_.empty()) ? cn.trans / pcZw_[size_t(ci)] : cn.trans;   // polymer shear factor
                const double q = U.mob[p].v / std::max(U.b[p].v, 1e-12) * tr * dphi;   // reservoir volume / day
                (p == PW ? qw : p == PO ? qo : qg)[ci] = q;
            }
        }
    } else {
        for (long long ci = 0; ci < (long long)nc; ++ci) {
            const Connection& cn = grid_.conns[ci];
            const CellFull& A = C[cn.a]; const CellFull& B = C[cn.b];
            for (int p = 0; p < 3; ++p) {
                if (eq_[p] < 0) continue;
                const double rho = 0.5 * (A.rho[p].v + B.rho[p].v);
                double dphi = A.P[p].v - B.P[p].v + rho * g * cn.dDepth;
                if (!thp_.empty() && thp_[size_t(ci)] > 0) { const double t = thp_[size_t(ci)]; dphi = dphi > t ? dphi - t : dphi < -t ? dphi + t : 0.0; }
                const CellFull& U = dphi >= 0 ? A : B;
                const double tr = (p == PW && !pcZw_.empty()) ? cn.trans / pcZw_[size_t(ci)] : cn.trans;   // polymer shear factor
                const double q = U.mob[p].v / std::max(U.b[p].v, 1e-12) * tr * dphi;   // reservoir volume / day
                (p == PW ? qw : p == PO ? qo : qg)[ci] = q;
            }
        }
    }
}

std::vector<Simulator::WellCellRate> Simulator::wellCellRates() const {
    std::vector<WellCellRate> out;
    const double g = U_.gravity;
    const int refP = O_ ? PO : PW;
    for (int w = 0; w < int(wells_.size()); ++w) {
        const WellRT& W = wells_[w];
        if (!W.open) continue;
        for (size_t k = 0; k < W.cell.size(); ++k) {
            const CellFull C = evalCell(W.cell[k]);
            const double pconn = W.pConn(k, g);
            const double dd = C.P[refP].v - pconn;
            double q = 0;
            if (W.type == WellType::Producer) {
                const double wm = k < W.wmulConn.size() ? W.wmulConn[k] : 1.0;   // polymer water mobility factor
                if (dd > 0) for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) q += C.mob[p].v / C.b[p].v * dd * W.wi[k] * (p == PW ? wm : 1.0);
            } else if (dd < 0) {
                const double wm = k < W.wmulConn.size() ? W.wmulConn[k] : 1.0;
                double lt = 0; for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) lt += C.mob[p].v / C.b[p].v * (p == PW ? wm : 1.0);
                q = lt * dd * W.wi[k];
            }
            out.push_back({W.cell[k], q, w});
        }
    }
    for (const AquiferRT& Q : aquifers_)
        for (size_t k = 0; k < Q.cell.size(); ++k) out.push_back({Q.cell[k], -Q.alpha[k] * Q.rate, -1 - Q.id});   // influx = source
    return out;
}

void Simulator::restoreExternal(const Snapshot& sn, int reportStep, const std::vector<double>& bhp) {
    const int na = grid_.numActive();
    if (int(sn.pressure.size()) != na) throw std::runtime_error("restoreExternal: snapshot does not match the active cells of the grid");
    Snapshot s2 = sn;
    auto fill = [&](std::vector<float>& v, float d) { if (int(v.size()) != na) v.assign(size_t(na), d); };
    fill(s2.sw, 0.f); fill(s2.sg, 0.f); fill(s2.rs, 0.f);
    if (int(s2.so.size()) != na) { s2.so.resize(size_t(na)); for (int a = 0; a < na; ++a) s2.so[size_t(a)] = 1.f - s2.sw[size_t(a)] - s2.sg[size_t(a)]; }
    restoreCells(s2);
    updateConnections(std::max(0, reportStep));
    for (size_t w = 0; w < wells_.size(); ++w) {
        const double b = w < bhp.size() ? bhp[w] : 0.0;
        wells_[w].bhp = b;
        wells_[w].open = std::isfinite(b) && b > 0 && wells_[w].nActive > 0;
        wells_[w].mode = Control::BHP;
    }
    if (poly_ && s2.poly.size() == size_t(na)) { for (int a = 0; a < na; ++a) pc_[a] = s2.poly[a]; polymerMultipliers(); }
    evalAll();
    computeHeads();
}

void Simulator::restoreCells(const Snapshot& sn) {
    const int na = grid_.numActive();
    for (int a = 0; a < na; ++a) {
        p_[a] = sn.pressure[a];
        if (iSw_ >= 0) sw_[a] = sn.sw[a];
        if (G_) {
            if (vap_ && a < int(sn.rv.size()) && sn.so[a] < 1e-7 && sn.sg[a] > 1e-7 && sn.rv[a] < fl(a).rvSat(p_[a]) * (1 - 1e-6)) { hasGas_[a] = 2; x_[a] = sn.rv[a]; }
            else if (!live_) { hasGas_[a] = 1; x_[a] = sn.sg[a]; }
            else if (sn.sg[a] > 1e-7 || sn.rs[a] >= fl(a).rsSat(p_[a]) * (1 - 1e-6)) { hasGas_[a] = 1; x_[a] = sn.sg[a]; }
            else { hasGas_[a] = 0; x_[a] = sn.rs[a]; }
        }
    }
}

void Simulator::restoreReport(int snapIndex) {
    if (snapIndex < 0 || snapIndex >= int(snaps_.size())) return;
    const Snapshot& sn = snaps_[snapIndex];
    const int na = grid_.numActive();
    restoreCells(sn);
    const SummaryPoint* sp = nullptr;
    for (const SummaryPoint& s : summary_) if (s.time <= sn.time + 1e-9) sp = &s;
    updateConnections(std::max(0, snapIndex - 1));    // connections of the report step that ended at this snapshot
    if (sp)
        for (size_t w = 0; w < wells_.size() && w < sp->wells.size(); ++w) {
            wells_[w].bhp = sp->wells[w].bhp; wells_[w].open = sp->wells[w].open && wells_[w].nActive > 0; wells_[w].mode = sp->wells[w].mode;
        }
    if (poly_ && sn.poly.size() == size_t(na)) { for (int a = 0; a < na; ++a) pc_[a] = sn.poly[a]; polymerMultipliers(); }
    // aquifer rates: use the summary total split by the stored fractions (approximate for restored states)
    evalAll();
    computeHeads();
}

// ======================================================================== equilibration
namespace {
struct Curve {
    std::vector<double> z, p;
    double at(double zz) const {
        if (zz <= z.front()) return p.front() + (zz - z.front()) * (p[1] - p[0]) / (z[1] - z[0]);
        if (zz >= z.back()) { size_t n = z.size(); return p[n - 1] + (zz - z[n - 1]) * (p[n - 1] - p[n - 2]) / (z[n - 1] - z[n - 2]); }
        size_t i = size_t(std::upper_bound(z.begin(), z.end(), zz) - z.begin()) - 1;
        return p[i] + (zz - z[i]) * (p[i + 1] - p[i]) / (z[i + 1] - z[i]);
    }
};
template <class F>
Curve integrate(const std::vector<double>& zs, double zref, double pref, double g, F rho) {
    Curve c; c.z = zs; c.p.assign(zs.size(), 0);
    auto stepTo = [&](double z0, double p0, double z1) {
        const int nsub = 4; double h = (z1 - z0) / nsub, p = p0, z = z0;
        for (int s = 0; s < nsub; ++s) {
            const double k1 = rho(p, z) * g, k2 = rho(p + 0.5 * h * k1, z + 0.5 * h) * g;
            p += h * k2; z += h;
        }
        return p;
    };
    size_t i0 = size_t(std::lower_bound(zs.begin(), zs.end(), zref) - zs.begin());
    double z = zref, p = pref;
    for (size_t i = i0; i < zs.size(); ++i) { p = stepTo(z, p, zs[i]); z = zs[i]; c.p[i] = p; }
    z = zref; p = pref;
    for (size_t i = i0; i-- > 0;) { p = stepTo(z, p, zs[i]); z = zs[i]; c.p[i] = p; }
    return c;
}
}

void Simulator::initialize() {
    const double g = U_.gravity;
    const int na = grid_.numActive();
    std::string rsInfo;
    pcScale_.clear();
    std::vector<double> swi;
    if (!case_.init.active() && W_ && O_ && case_.grid.swatinit.size() == size_t(case_.grid.size())) {
        swi = case_.grid.swatinit;
        pcScale_.assign(size_t(na), 1.0);
    }
    int nSwatScaled = 0, nSwatIgnored = 0;
    if (case_.init.active()) {
        // ---- enumerated initial conditions (PRESSURE / SWAT / SGAS / RS)
        const InitialState& I = case_.init;
        auto val = [&](const std::vector<double>& v, int a, double def) {
            const size_t q = size_t(grid_.globalOf[size_t(a)]);
            return q < v.size() ? v[q] : def;
        };
        for (int a = 0; a < na; ++a) {
            p_[a] = val(I.pressure, a, 0.0);
            const double sw = (W_ && !O_) ? 1.0 : (W_ ? val(I.swat, a, eps_.empty() ? st(a).swco() : eps_[size_t(a)].swl) : 0.0);
            const double sg = G_ ? val(I.sgas, a, 0.0) : 0.0;
            sw_[a] = sw;
            if (G_) {
                const double rsat = live_ ? fl(a).rsSat(p_[a]) : 0.0;
                const double rvsat = vap_ ? fl(a).rvSat(p_[a]) : 0.0;
                if (vap_ && O_ && sg > 1e-9 && sw + sg >= 1.0 - 1e-9 && sw < 1.0 - 1e-9) { hasGas_[a] = 2; x_[a] = std::min(val(I.rv, a, rvsat), rvsat); }
                else if (!live_ || sg > 1e-9 || (iSw_ >= 0 && sw > 1.0 - 1e-9)) { hasGas_[a] = 1; x_[a] = sg; }
                else { hasGas_[a] = 0; x_[a] = std::min(val(I.rs, a, rsat), rsat); }
            }
        }
        rsInfo = "enumerated";
    } else {
        // ---- hydrostatic equilibration, one pass per EQLNUM region
        const int neql = 1 + int(case_.equilExtra.size());
        std::vector<std::vector<int>> cellsOf(static_cast<size_t>(neql));
        for (int a = 0; a < na; ++a) cellsOf[size_t(eqlReg_[size_t(a)])].push_back(a);
        initRs_.assign(size_t(neql), -1.0);
        initRv_.assign(size_t(neql), -1.0);
        for (int e = 0; e < neql; ++e) {
            const std::vector<int>& cells = cellsOf[size_t(e)];
            if (cells.empty()) continue;
            const EquilInput& E = e == 0 ? case_.equil : case_.equilExtra[size_t(e - 1)];
            // fluid of the region: the most common PVTNUM among its cells
            std::vector<int> cnt(fluids_.size(), 0);
            for (int a : cells) ++cnt[size_t(pvtReg_[size_t(a)])];
            const FluidModel& F = fluids_[size_t(std::max_element(cnt.begin(), cnt.end()) - cnt.begin())];
            const DensitySC& rho = F.rhoSC();
            const double gf = F.gasFactor();
            double zmin = 1e300, zmax = -1e300;
            for (int a : cells) {
                zmin = std::min(zmin, grid_.depth[size_t(a)] - grid_.dzCell[size_t(a)]);
                zmax = std::max(zmax, grid_.depth[size_t(a)] + grid_.dzCell[size_t(a)]);
            }
            double lo = std::min({zmin, E.datumDepth, E.owc, E.goc});
            double hi = std::max({zmax, E.datumDepth, E.owc, E.goc});
            lo -= 1.0; hi += 1.0;
            const int N = 400;
            std::vector<double> zs(N + 1);
            for (int i = 0; i <= N; ++i) zs[i] = lo + (hi - lo) * i / N;

            const bool rsvd = live_ && E.rsvdDepth.size() >= 1 && E.rsvdDepth.size() == E.rsvdRs.size();
            Table1D rsTab; rsTab.clampEnds = true;
            if (rsvd) { rsTab.x = E.rsvdDepth; rsTab.y = E.rsvdRs; if (rsTab.x.size() == 1) { rsTab.x.push_back(rsTab.x[0] + 1); rsTab.y.push_back(rsTab.y[0]); } }
            double rsConst = E.initialRs;
            auto rsAt = [&](double p, double z) {
                if (!live_) return 0.0;
                const double cap = F.rsSat(p);
                return std::min(rsvd ? rsTab.eval(z) : rsConst, cap);
            };
            // wet gas: Rv from RVVD, a constant, or saturated at the GOC pressure (ECLIPSE default), capped by the dew point
            const bool rvvd = vap_ && E.rvvdDepth.size() >= 1 && E.rvvdDepth.size() == E.rvvdRv.size();
            Table1D rvTab; rvTab.clampEnds = true;
            if (rvvd) { rvTab.x = E.rvvdDepth; rvTab.y = E.rvvdRv; if (rvTab.x.size() == 1) { rvTab.x.push_back(rvTab.x[0] + 1); rvTab.y.push_back(rvTab.y[0]); } }
            double rvConst = E.initialRv;
            auto rvAt = [&](double p, double z) {
                if (!vap_) return 0.0;
                return std::min(rvvd ? rvTab.eval(z) : rvConst, F.rvSat(p));
            };
            auto rhoW = [&](double p, double) { return rho.water * F.bwS(p); };
            auto rhoG = [&](double p, double z) {
                if (!vap_) return rho.gas * gf * F.bgS(p);
                const double rv = rvAt(p, z);
                return (rho.gas * gf + rv * rho.oil) * F.bg(AD3(p), AD3(rv)).v;
            };
            auto rhoO = [&](double p, double z) { const double rs = rsAt(p, z); return (rho.oil + rs * rho.gas * gf) * F.boS(p, rs); };

            enum { DW, DO, DG } dph = DO;
            if (!O_) dph = DW;
            else if (W_ && E.datumDepth > E.owc) dph = DW;
            else if (G_ && E.datumDepth < E.goc) dph = DG;

            Curve cw, co, cg;
            const bool iterRs = live_ && !rsvd && E.initialRs < 0;
            const bool iterRv = vap_ && !rvvd && E.initialRv < 0;
            for (int iter = 0; iter < ((iterRs || iterRv) ? 3 : 1); ++iter) {
                if (iterRs && iter == 0) rsConst = 1e30;   // first pass: saturated everywhere
                if (iterRv && iter == 0) rvConst = 1e30;
                if (dph == DO) {
                    co = integrate(zs, E.datumDepth, E.datumPressure, g, rhoO);
                    if (W_) cw = integrate(zs, E.owc, co.at(E.owc) - E.pcowOWC, g, rhoW);
                    if (G_) cg = integrate(zs, E.goc, co.at(E.goc) + E.pcogGOC, g, rhoG);
                } else if (dph == DW) {
                    cw = integrate(zs, E.datumDepth, E.datumPressure, g, rhoW);
                    if (O_) co = integrate(zs, E.owc, cw.at(E.owc) + E.pcowOWC, g, rhoO);
                    if (G_) cg = integrate(zs, E.goc, co.at(E.goc) + E.pcogGOC, g, rhoG);
                } else {
                    cg = integrate(zs, E.datumDepth, E.datumPressure, g, rhoG);
                    co = integrate(zs, E.goc, cg.at(E.goc) - E.pcogGOC, g, rhoO);
                    if (W_) cw = integrate(zs, E.owc, co.at(E.owc) - E.pcowOWC, g, rhoW);
                }
                if (iterRs) rsConst = F.rsSat(co.at(E.goc));
                if (iterRv) rvConst = F.rvSat(G_ ? cg.at(E.goc) : co.at(E.goc));
            }
            if (vap_ && !rvvd) initRv_[size_t(e)] = rvConst;
            if (vap_) {
                if (!rsInfo.empty()) rsInfo += ", ";
                rsInfo += std::string("Rv ") + (rvvd ? std::string("RVVD") : fmt("%.5g", rvConst));
            }
            if (live_ && !rsvd) initRs_[size_t(e)] = rsConst;
            if (live_) {
                if (!rsInfo.empty()) rsInfo += ", ";
                rsInfo += rsvd ? std::string("RSVD") : fmt("%.5g", rsConst);
            }

            const int ns = E.nacc == 0 ? 1 : 2 * std::abs(E.nacc);
            for (int a : cells) {
                const SatFunctions& S = st(a);
                const double top = grid_.depth[a] - 0.5 * grid_.dzCell[a];
                double swa = 0, sga = 0, poa = 0;
                // OPM / ECLIPSE: the oil pressure of a cell is the oil curve, except where the oil phase is absent
                // (water zone at Swu, gas cap at Sgu or an overlapping G/W transition): there it follows the
                // pressure of the present phase, Po = Pw + Pcow(Sw) or Po = Pg - Pcgo(Sg).
                const auto& epc = ep(a) ? *ep(a) : S.tableEndPoints();
                for (int s = 0; s < ns; ++s) {
                    const double z = (ns == 1) ? grid_.depth[a] : top + grid_.dzCell[a] * (s + 0.5) / ns;
                    double sw = (W_ && !O_) ? 1.0 : 0.0, sg = 0.0;
                    double po = O_ ? co.at(z) : cw.at(z);
                    bool overlap = false;
                    if (W_ && O_) sw = S.swInverseFromPc(co.at(z) - cw.at(z), ep(a));
                    if (G_ && O_) {
                        sg = S.sgInverseFromPc(cg.at(z) - co.at(z), ep(a));
                        if (W_ && sg > 0 && sw + sg > 1.0) { sw = S.swInverseFromPcGW(cg.at(z) - cw.at(z), ep(a)); sg = 1.0 - sw; overlap = true; }
                        sg = std::min(sg, 1.0 - sw);
                    }
                    if (O_) {
                        if (overlap) po = cg.at(z) - S.pcog(AD3(sg), ep(a)).v;
                        if (W_ && sw + 1e-6 > epc.swu) po = cw.at(z) + S.pcow(AD3(epc.swu), ep(a)).v;
                        if (G_ && sg + 1e-6 > epc.sgu) po = cg.at(z) - S.pcog(AD3(epc.sgu), ep(a)).v;
                    }
                    swa += sw; sga += sg; poa += po;
                }
                swa /= ns; sga /= ns; poa /= ns;
                if (!swi.empty() && W_ && O_) {
                    // SWATINIT: honour the given water saturation by scaling Pcow so that the cell is in equilibrium
                    // (ECLIPSE / OPM Flow behaviour). Cells in the water zone or with a zero Pc curve keep the
                    // equilibrium saturation.
                    const double s0 = swi[size_t(grid_.globalOf[size_t(a)])];
                    if (s0 >= 0) {
                        const double pcEq = co.at(grid_.depth[a]) - cw.at(grid_.depth[a]);
                        const double pcT = S.pcow(AD3(std::min(1.0, s0)), ep(a)).v;
                        const double pcAbs = std::max(std::fabs(S.pcow(AD3(0.0), ep(a)).v), std::fabs(pcT));
                        if (pcEq > 1e-9 && pcT > 1e-9 * std::max(1.0, pcAbs)) {
                            pcScale_[size_t(a)] = pcEq / pcT; swa = std::min(1.0, s0); ++nSwatScaled;
                            if (G_) sga = std::min(sga, 1.0 - swa);
                        } else if (pcEq > 1e-9) ++nSwatIgnored;
                    }
                }
                p_[a] = poa;
                sw_[a] = swa;
                if (G_) {
                    if (vap_ && O_ && sga > 1e-9 && swa + sga >= 1.0 - 1e-9 && !(iSw_ >= 0 && swa > 1.0 - 1e-9)) {
                        const double pg = p_[a] + S.pcog(AD3(sga), ep(a)).v;   // Rv at the gas pressure (OPM)
                        hasGas_[a] = 2; x_[a] = std::min(rvAt(pg, grid_.depth[a]), fl(a).rvSat(pg));   // gas zone: no oil
                    }
                    else if (sga > 1e-9 || !live_ || (iSw_ >= 0 && swa > 1.0 - 1e-9)) { hasGas_[a] = 1; x_[a] = sga; }
                    else { hasGas_[a] = 0; x_[a] = std::min(rsAt(p_[a], grid_.depth[a]), fl(a).rsSat(p_[a])); }
                }
            }
        }
    }
    if (!swi.empty()) {
        double smin = 1e300, smax = 0;
        for (double s : pcScale_) if (s != 1.0) { smin = std::min(smin, s); smax = std::max(smax, s); }
        log_(fmt("SWATINIT: Pcow scaled in %g cells", double(nSwatScaled)) +
             (nSwatScaled ? fmt(" (multiplier %.3g", smin) + fmt("..%.3g)", smax) : std::string()) +
             (nSwatIgnored ? fmt(", %g cells above the contact keep the equilibrium Sw (Pc = 0 at SWATINIT)", double(nSwatIgnored)) : std::string()));
    }
    polymerSetup();
    evalAll();
    setupThresholds();
    storeOld();
    time_ = 0;
    initAquifers();
    summary_.clear(); snaps_.clear();
    for (double& v : cum_) v = 0;
    pcCum_[0] = pcCum_[1] = 0;
    for (auto& w : wells_) { w.rates = WellRates{}; }
    record(0, 0, 0);
    snapshot();
    const InPlace ip = inPlace();
    log_(std::string("Initial fluids in place (") + U_.liq + " / " + U_.gas + "):");
    log_(fmt("  oil %.6g   water %.6g   free gas %.6g   dissolved gas %.6g", ip.oil, ip.water, ip.gasFree, ip.gasDissolved));
    if (vap_) log_(fmt("  vaporised oil (in gas) %.6g   total oil %.6g", ip.oilVaporized, ip.oil + ip.oilVaporized));
    if (live_) log_("  initial Rs: " + rsInfo + " " + U_.rs);
    if (fluids_.size() > 1 || sats_.size() > 1 || !case_.equilExtra.empty())
        log_(fmt("  regions: PVT %g, saturation %g, ", double(fluids_.size()), double(sats_.size())) +
             fmt("equilibration %g", double(1 + case_.equilExtra.size())));
}

// THPRES: per-connection threshold pressures between EQLNUM regions; defaulted values = largest initial potential
// difference (over the mobile phases) across the region boundary, so the initial equilibrium is preserved
void Simulator::setupThresholds() {
    thp_.clear();
    const auto& eqn = case_.grid.eqlnum;
    if (case_.thpres.empty() || eqn.size() != size_t(case_.grid.size())) return;
    const size_t nc = grid_.conns.size();
    thp_.assign(nc, 0.0);
    std::map<std::pair<int, int>, double> computed;
    auto regOf = [&](int a) { return eqn[size_t(grid_.globalOf[size_t(a)])]; };
    const double g = U_.gravity;
    // pass 1: largest initial potential difference per region pair
    for (size_t ci = 0; ci < nc; ++ci) {
        const Connection& cn = grid_.conns[ci];
        int ra = regOf(cn.a), rb = regOf(cn.b);
        if (ra == rb) continue;
        const auto key = std::make_pair(std::min(ra, rb), std::max(ra, rb));
        const CellFlux& A = cf_[cn.a]; const CellFlux& B = cf_[cn.b];
        double m = 0;
        for (int p = 0; p < 3; ++p) {
            if (eq_[p] < 0) continue;
            const double dphi = A.P[p].v - B.P[p].v + 0.5 * (A.rho[p].v + B.rho[p].v) * g * cn.dDepth;
            const CellFlux& Up = dphi >= 0 ? A : B;
            if (Up.mob[p].v > 0) m = std::max(m, std::fabs(dphi));
        }
        double& v = computed[key];
        v = std::max(v, m);
    }
    int nset = 0;
    for (size_t ci = 0; ci < nc; ++ci) {
        const Connection& cn = grid_.conns[ci];
        int ra = regOf(cn.a), rb = regOf(cn.b);
        if (ra == rb) continue;
        for (const ThresholdPressure& t : case_.thpres)
            if ((t.r1 == ra && t.r2 == rb) || (t.r1 == rb && t.r2 == ra)) {
                thp_[ci] = t.value >= 0 ? t.value : computed[std::make_pair(std::min(ra, rb), std::max(ra, rb))];
                ++nset;
            }
    }
    std::string info;
    for (const ThresholdPressure& t : case_.thpres) {
        const double v = t.value >= 0 ? t.value : computed[std::make_pair(std::min(t.r1, t.r2), std::max(t.r1, t.r2))];
        info += fmt(" %g-", double(t.r1)) + fmt("%g:", double(t.r2)) + fmt(" %.4g", v) + (t.value < 0 ? "(init)" : "");
    }
    log_(fmt("THPRES: threshold pressures on %g connections;", double(nset)) + info + " " + U_.p);
    if (nset == 0) thp_.clear();
}

std::vector<Simulator::ConnInfo> Simulator::connections() const {
    std::vector<ConnInfo> out;
    for (const WellRT& w : wells_) {
        const std::string& name = case_.wells[size_t(w.input)].name;
        for (size_t c = 0; c < w.cell.size(); ++c) {
            ConnInfo ci; ci.well = name;
            grid_.ijk(grid_.globalOf[size_t(w.cell[c])], ci.i, ci.j, ci.k);
            ci.wi = c < w.wi.size() ? w.wi[c] : 0; ci.depth = c < w.cdepth.size() ? w.cdepth[c] : 0;
            out.push_back(ci);
        }
    }
    return out;
}

} // namespace ressim
