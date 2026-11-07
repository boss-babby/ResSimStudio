#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <climits>
#include <cstring>
#include <cstddef>
#include <cctype>
#include "ressim/simulator.h"
#include "ressim/pyaction.h"
#include "ressim/nonlinear.h"
#include "ressim/binio.h"
#include "ressim/components.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <chrono>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "ressim/phasefield.h"
#include "ressim/hf_field.h"
#include "ressim/udq.h"

namespace ressim {

// C28 in-run phase-field link (see pfReportStep)
struct Simulator::PfLink {
    int frac = 0;
    phasefield::Input in;
    std::thread th; std::mutex m; std::condition_variable cv;
    bool want = false, stop = false, done = false, ready = false, started = false;
    double target = 0; phasefield::Step last; double tipAng = 0; std::string error;
    double vinj = 0; size_t sumIdx = 0;                 // cumulative injected reservoir volume (m3), summary points used
    double C = 0, Vlast = 0, wRef = 0;                  // C29 TWOWAY: compliance V / p_net (m2/m/Pa), last volume, reference aperture (m)
    double vinjLast = 0;                                // C29 TWOWAY: injected volume at the last phase-field step (mass cap)
    double Hm = 1, toM = 1, toPa = 1, sigmaCase = 0;    // fracture height (m), unit factors, closure stress (case units)
    std::shared_ptr<phasefield::Field> field;           // C30: field of the latest phase-field step (under m)
    std::vector<int> wells;                             // indices into wells_ completed in the fracture's host cells
    ~PfLink() {
        { std::lock_guard<std::mutex> lk(m); stop = true; }
        cv.notify_all();
        if (th.joinable()) th.join();
    }
    void start() {
        started = true;
        in.nextVolume = [this](int, double& V) {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [&] { return want || stop; });
            if (stop) return false;
            V = target; want = false;
            return true;
        };
        in.onField = [this](const phasefield::Field& F) {
            auto f = std::make_shared<phasefield::Field>(F);
            std::lock_guard<std::mutex> lk(m);
            field = std::move(f);
        };
        in.onStep = [this](const phasefield::Step& S, double ang) {
            std::lock_guard<std::mutex> lk(m);
            last = S; tipAng = ang; ready = true;
            cv.notify_all();
        };
        th = std::thread([this] {
            phasefield::Output out;
            phasefield::run(in, out);
            std::lock_guard<std::mutex> lk(m);
            done = true; error = out.error.empty() ? "the crack reached the edge of the phase-field core (maximum extent)" : out.error;
            cv.notify_all();
        });
    }
    bool advance(double V, phasefield::Step& S) {
        if (!started) start();
        {
            std::lock_guard<std::mutex> lk(m);
            if (done) return false;
            target = V; want = true; ready = false;
        }
        cv.notify_all();
        std::unique_lock<std::mutex> lk(m);
        cv.wait(lk, [&] { return ready || done; });
        if (!ready) return false;
        S = last;
        return true;
    }
};


using AD4 = AD<4>;

static std::string fmt(const char* f, double a = 0, double b = 0, double c = 0, double d = 0) {
    char buf[512]; std::snprintf(buf, sizeof(buf), f, a, b, c, d); return buf;
}

// ======================================================================== setup
void Simulator::setup(const Case& c0) {
    const std::string err = validateCase(c0);
    if (!err.empty()) throw std::runtime_error(err);
    case_ = c0;
    saltPv0_.clear(); saltT0_.clear();   // C29: salt references belong to the grid built below
    expandLgrArrays(case_);          // local grid refinements: per-cell arrays also cover the local cells
    expandMincArrays(case_);         // MINC: per-cell arrays also cover the inner matrix shells
    applyDpnum(case_);               // DPNUM = 0: single porosity blocks
    prepareHydraulicFractures(case_);   // C15 hydraulic fractures (EDFM): per-cell arrays also cover the fracture cells
    const Case& c = case_;
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
        // C34 SATOPTS DIRECT: KRNUMX / KRNUMY / KRNUMZ tables for the flux terms in each direction (wells use SATNUM)
        dirKr_ = false;
        const std::vector<int>* kn[3] = {&c.grid.krnumx, &c.grid.krnumy, &c.grid.krnumz};
        int nDiff = 0;
        for (int d = 0; d < 3; ++d) {
            krReg_[d].clear(); cfDir_[d].clear();
            if (kn[d]->empty()) continue;
            krReg_[d] = satReg_;
            for (int a = 0; a < na0; ++a) {
                const size_t q = size_t(grid_.globalOf[size_t(a)]);
                const int r = q < kn[d]->size() ? (*kn[d])[q] - 1 : -1;
                if (r >= 0 && r < nsat) krReg_[d][size_t(a)] = r;
                if (krReg_[d][size_t(a)] != satReg_[size_t(a)]) ++nDiff;
            }
            dirKr_ = true;
        }
        if (dirKr_) log_(fmt("Directional relative permeability (KRNUMX/Y/Z): %g cell-direction pairs differ from SATNUM", double(nDiff)));
    }
    {   // ---- C15 Langmuir adsorption and stress-dependent permeability
        const int na0 = grid_.numActive();
        langK_.clear(); sGam_.clear(); pI_.clear(); langPL_ = 0;
        const GridInput& gi = c.grid;
        if (G_ && gi.langVL > 0 && gi.langPL > 0) {
            langK_.assign(size_t(na0), 0.0); langPL_ = gi.langPL;
            // surface gas per bulk length^3: METRIC m3 * t/m3 * sm3/t; FIELD ft3 -> m3 -> short ton * scf/ton -> Mscf
            const double conv = c.units == UnitSystem::Field ? 0.0283168 * 1.10231 / 1000.0 : 1.0;
            for (int a = 0; a < na0; ++a) {
                const int q = grid_.globalOf[size_t(a)];
                if (grid_.isHf(q) || grid_.isFracture(q)) continue;   // adsorbed gas lives in the matrix
                const double phi = q < int(gi.poro.size()) ? gi.poro[size_t(q)] : 0.0, ntg = q < int(gi.ntg.size()) ? gi.ntg[size_t(q)] : 1.0;
                langK_[size_t(a)] = grid_.bulkVol[size_t(q)] * ntg * std::max(0.0, 1 - phi) * gi.langRho * gi.langVL * conv;
            }
        }
        bool anyG = gi.permGamma > 0;
        for (const auto& h : grid_.hfCells) if (h.gamma > 0) anyG = true;
        if (anyG) {
            sGam_.assign(size_t(na0), 0.0);
            for (int a = 0; a < na0; ++a) {
                const int q = grid_.globalOf[size_t(a)];
                sGam_[size_t(a)] = grid_.isHf(q) ? grid_.hfCells[size_t(q - grid_.hfBase())].gamma : gi.permGamma;
            }
        }
    }
    setupEndPointScaling(c);
    setupHysteresis(c);
    vap1_ = std::max(0.0, c.props.vap1); vap2_ = std::max(0.0, c.props.vap2);
    setNumThreads(c.tuning.threads);

    const int na = grid_.numActive();
    p_.assign(na, 0); sw_.assign(na, 0); x_.assign(na, 0); hasGas_.assign(na, 0);
    accOld_.assign(size_t(na) * 3, 0.0);
    cf_.assign(na, CellFlux{});
    setupGravityDrainage(c);
    buildWells();
    wellNames_.clear(); for (auto& w : c.wells) wellNames_.push_back(w.name);
    buildAquifers();
    buildPattern();
    setupHfGrowth();
    {   // C30 fracture stress / opening field output (one-way, display only)
        hff_.reset();
        auto F = std::make_shared<hffield::Field>();
        std::string note;
        if (F->setup(*this, note)) {
            hff_ = F;
            const auto& st = F->stats();
            log_(fmt("Fracture field output: %g pieces, ", st.pieces) + fmt("%g matrix cells", st.cells) + fmt(" within %.4g m", st.radiusM) +
                 (st.halfSpace ? std::string(", half space") : std::string(", full space")) + fmt(", E %.4g GPa", st.E * 1e-9) + fmt(" nu %.3g", st.nu) +
                 fmt(" (setup %.3g s)", st.setupSeconds) + (note.empty() ? std::string() : " - " + note));
        } else if (c.hfField.on == 1) log_("NOTE: HFFIELD: " + note);
    }
    polymerSetup();
    useCPR_ = c.tuning.linearSolver == 1 || (c.tuning.linearSolver == 2 && na >= 500);
    time_ = 0; step_ = 0; dtNext_ = c.tuning.initialDt;
    for (double& v : cum_) v = 0;
    summary_.clear(); snaps_.clear(); stats_ = {}; stats_.cpr = useCPR_;
    log_(fmt("Grid %g x %g x %g, active cells %g", grid_.nx, grid_.ny, grid_.nz, na) +
         fmt(", unknowns %g, threads %g", double(J_.n), numThreads()) +
         (useCPR_ ? ", linear solver CPR-AMG" : ", linear solver ILU(0)-BiCGSTAB"));
    udqSetup();                               // C31 UDQ / UDA / ACTIONX
    pyActionSetup();                          // C37 PYACTION
}

void Simulator::buildWells() {
    const GridInput& g = case_.grid;
    const int na = grid_.numActive();
    cellWells_.assign(na, {});
    wells_.clear();
    std::vector<std::vector<int>> hfOfHost;   // C15: fracture cells per host (wells completed in a host feed its fractures)
    if (grid_.numHfCells()) {
        hfOfHost.resize(size_t(grid_.hfBase()));
        for (int k = 0; k < grid_.numHfCells(); ++k) hfOfHost[size_t(grid_.hfCells[size_t(k)].host)].push_back(k);
    }
    auto at = [&](const std::vector<double>& v, int q, double def) { return v.empty() ? def : v[q]; };
    for (int wi = 0; wi < int(case_.wells.size()); ++wi) {
        const WellInput& in = case_.wells[wi];
        WellRT w; w.input = wi; w.type = in.type;
        const int widx = int(wells_.size());
        // C48: defaulted BHP reference depth = the first connections the well ever had, SHUT ones included (OPM
        // Well::updateRefDepth: depth of connection 0, fixed once set; e.g. COMPDAT ... SHUT + WELOPEN of a later COMPLUMP)
        int rdStep = std::numeric_limits<int>::max(); double rdDepth = 1e30;
        for (const Completion& cp : in.completions) {
            const int lg = cp.lgr;
            if (lg >= int(grid_.lgr.size())) throw std::runtime_error("Well " + in.name + ": completion in an unknown local grid");
            const int NX = lg >= 0 ? grid_.lgr[size_t(lg)].nx : g.nx, NY = lg >= 0 ? grid_.lgr[size_t(lg)].ny : g.ny,
                      NZ = lg >= 0 ? grid_.lgr[size_t(lg)].nz : g.nz;
            for (int k = std::min(cp.k1, cp.k2); k <= std::max(cp.k1, cp.k2); ++k) {
                if (cp.i < 0 || cp.i >= NX || cp.j < 0 || cp.j >= NY || k < 0 || k >= NZ)
                    throw std::runtime_error("Well " + in.name + ": completion outside the " + (lg >= 0 ? "local grid " + grid_.lgr[size_t(lg)].name : std::string("grid")));
                // DUALPORO: wells connect to the fracture cell of the block (a matrix-layer completion is moved)
                const int kk = (g.dual == 1 && k < g.nzMatrix()) ? k + g.nzMatrix() : k;
                std::vector<int> targets{lg >= 0 ? grid_.localCell(lg, cp.i, cp.j, k) : g.index(cp.i, cp.j, kk)};
                if (lg < 0 && grid_.isHost(targets[0])) targets = grid_.hostWellColumn(targets[0], cp.dir);   // global COMPDAT in a refined cell
                for (const int q : targets) {
                const int a = grid_.activeOf[q];
                if (a >= 0) {
                    const int st = cp.defStep >= 0 ? std::min(cp.defStep, std::max(0, cp.fromStep)) : std::max(0, cp.fromStep);
                    if (st < rdStep) { rdStep = st; rdDepth = grid_.depth[size_t(a)]; }
                    else if (st == rdStep) rdDepth = std::min(rdDepth, grid_.depth[size_t(a)]);
                }
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
                // Peaceman connection factor of permeabilities (kx, ky, kz)
                auto peaceman = [&](double KX, double KY, double KZ) {
                    double a1, a2;
                    if (cp.dir == 'X' || cp.dir == 'x') { a1 = KY; a2 = KZ; }
                    else if (cp.dir == 'Y' || cp.dir == 'y') { a1 = KX; a2 = KZ; }
                    else { a1 = KX; a2 = KY; }
                    const double kh = cp.kh > 0 ? cp.kh : std::sqrt(a1 * a2) * h;
                    const double r12 = std::sqrt(a2 / a1), r21 = std::sqrt(a1 / a2);
                    const double ro = 0.28 * std::sqrt(r12 * d1 * d1 + r21 * d2 * d2) /
                                      (std::pow(a2 / a1, 0.25) + std::pow(a1 / a2, 0.25));
                    return U_.cDarcy * 2.0 * M_PI * kh / (std::log(ro / (0.5 * cp.diameter)) + cp.skin);
                };
                std::array<double, 3> dl{0.0, 0.0, 0.0};
                if (wi_ <= 0) {
                    wi_ = peaceman(kx, ky, kz);
                    if (!(wi_ > 0)) throw std::runtime_error("Well " + in.name + ": non-positive connection factor");
                    // d ln WI / d ln k (adjoint field gradients): central differences of the formula in ln k
                    const double e = 1e-4, ep = std::exp(e), em = std::exp(-e);
                    dl[0] = std::log(peaceman(kx * ep, ky, kz) / peaceman(kx * em, ky, kz)) / (2 * e);
                    dl[1] = std::log(peaceman(kx, ky * ep, kz) / peaceman(kx, ky * em, kz)) / (2 * e);
                    dl[2] = std::log(peaceman(kx, ky, kz * ep) / peaceman(kx, ky, kz * em)) / (2 * e);
                    for (double& v : dl) if (!std::isfinite(v)) v = 0;
                }
                w.dlwi.push_back(dl);
                wi_ *= in.piMult * cp.mult;
                w.cell.push_back(a); w.wi.push_back(wi_); w.cdepth.push_back(grid_.depth[a]);
                w.wiBase.push_back(wi_); w.cFrom.push_back(std::max(0, cp.fromStep)); w.cUntil.push_back(cp.untilStep); w.cComp.push_back(int(&cp - in.completions.data()));
                w.cDef.push_back(cp.defStep >= 0 ? std::min(cp.defStep, std::max(0, cp.fromStep)) : std::max(0, cp.fromStep));
                w.cMult.push_back(cp.mult); w.cWelpi.push_back(1.0);
                auto& cw = cellWells_[a];
                int idx = -1;
                for (int m = 0; m < int(cw.size()); ++m) if (cw[m] == widx) idx = m;
                if (idx < 0) { idx = int(cw.size()); cw.push_back(widx); }
                w.cellWellIdx.push_back(idx);
                if (!hfOfHost.empty() && q < int(hfOfHost.size()))
                    for (const int k : hfOfHost[size_t(q)]) {   // EDFM well index of the fracture piece: 2 pi k_f w / ln(r_e / r_w)
                        const Grid::HfCell& H = grid_.hfCells[size_t(k)];
                        const int fa = grid_.activeOf[size_t(H.global)];
                        if (fa < 0) continue;
                        const double rw = 0.5 * cp.diameter, re = 0.14 * std::sqrt(H.len * H.len + H.h * H.h);
                        const double wf = U_.cDarcy * 2 * M_PI * H.kf * H.w / std::max(std::log(std::max(re, 1.01 * rw) / rw), 1.0) * in.piMult * cp.mult;
                        w.rw.push_back(rw); w.re.push_back(re); w.plen.push_back(H.h); w.dlwi.push_back({0.0, 0.0, 0.0});
                        w.cell.push_back(fa); w.wi.push_back(wf); w.cdepth.push_back(grid_.depth[size_t(fa)]);
                        w.wiBase.push_back(wf); w.cFrom.push_back(std::max(0, cp.fromStep)); w.cUntil.push_back(cp.untilStep); w.cComp.push_back(int(&cp - in.completions.data()));
                w.cDef.push_back(cp.defStep >= 0 ? std::min(cp.defStep, std::max(0, cp.fromStep)) : std::max(0, cp.fromStep));
                w.cMult.push_back(cp.mult); w.cWelpi.push_back(1.0);
                        auto& cwf = cellWells_[size_t(fa)];
                        int ix = -1;
                        for (int m = 0; m < int(cwf.size()); ++m) if (cwf[m] == widx) ix = m;
                        if (ix < 0) { ix = int(cwf.size()); cwf.push_back(widx); }
                        w.cellWellIdx.push_back(ix);
                    }
                }
            }
        }
        w.rowCells = w.cell;
        std::sort(w.rowCells.begin(), w.rowCells.end());
        w.rowCells.erase(std::unique(w.rowCells.begin(), w.rowCells.end()), w.rowCells.end());
        for (int a : w.cell) w.cpos.push_back(int(std::lower_bound(w.rowCells.begin(), w.rowCells.end(), a) - w.rowCells.begin()));
        if (in.refDepth > 0) w.refDepth = in.refDepth;
        else if (rdDepth < 1e29) w.refDepth = rdDepth;
        else { w.refDepth = 1e30; for (double d : w.cdepth) w.refDepth = std::min(w.refDepth, d); if (w.cdepth.empty()) w.refDepth = 0; }
        {   // C41: crossflow flag = WELSPECS item 10 (4th item of the tail 7..; N* repeats defaults)
            std::istringstream is(in.specsTail); std::string t; std::vector<std::string> it;
            while (is >> t && it.size() < 8) {
                if (t == "/") break;
                const size_t st = t.find('*');
                if (st != std::string::npos && st > 0 && t.find_first_not_of("0123456789") == st) { const int nrep = std::atoi(t.substr(0, st).c_str()); for (int r = 0; r < nrep; ++r) it.push_back("*"); }
                else it.push_back(t);
            }
            std::string x = it.size() > 3 ? it[3] : "*";
            for (char& ch : x) ch = char(std::toupper((unsigned char)ch));
            x.erase(std::remove(x.begin(), x.end(), '\''), x.end());
            w.xflow = !(x == "NO" || x == "N");
            // C44: item 9 = automatic shut-in instruction for wells closed by the simulator (SHUT default, STOP)
            std::string a9 = it.size() > 2 ? it[2] : "*";
            for (char& ch : a9) ch = char(std::toupper((unsigned char)ch));
            a9.erase(std::remove(a9.begin(), a9.end(), '\''), a9.end());
            w.autoStop = a9 == "STOP";
            if (std::getenv("RESSIM_NO_XFLOW")) w.xflow = false;
        }
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
                        const int q0 = g.index(i, j, k);
                        if (grid_.activeOf[q0] < 0 && !grid_.isHost(q0)) continue;
                        int ni = i, nj = j, nk = k;
                        const int s = c.minus ? -1 : 1;
                        if (c.face == 'X') ni += s; else if (c.face == 'Y') nj += s; else nk += s;
                        if (ni >= 0 && ni < g.nx && nj >= 0 && nj < g.ny && nk >= 0 && nk < g.nz &&
                            (grid_.activeOf[g.index(ni, nj, nk)] >= 0 || grid_.isHost(g.index(ni, nj, nk)))) continue;     // internal face: not connected
                        // a refined cell connects through the local cells on that face
                        const std::vector<int> cells = grid_.isHost(q0) ? grid_.localCellsOnFace(q0, c.face, c.minus) : std::vector<int>{q0};
                        for (const int q : cells) {
                            const int a = grid_.activeOf[q];
                            if (a < 0) continue;
                            const double A = (c.coeff > 0 ? c.coeff : grid_.faceArea(q, c.face, c.minus)) * c.mult;   // C34: explicit influx coefficient
                            if (A <= 0) continue;
                            Q.cell.push_back(a); area.push_back(A); Q.cdepth.push_back(grid_.depth[a]);
                        }
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
            Q.tdTab = in.tdTab; Q.pdTab = in.pdTab;
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
    static const bool scalarEnv = std::getenv("RESSIM_SCALAR_JAC") || std::getenv("RESSIM_DUMP") || std::getenv("RESSIM_DUMP_FAIL");
    blkAsm_ = case_.tuning.blockSolver && !adjRec_ && !forceScalarJ_ && !scalarEnv;
    if (blkAsm_) schur_.bindCells(na, nv, blkStart_, blkCells_);
    J_ = CSR{}; J_.n = N; J_.rowPtr.assign(N + 1, 0);
    for (int a = 0; a < na; ++a) {
        const int len = (blkAsm_ ? 0 : (blkStart_[a + 1] - blkStart_[a]) * nv) + int(cellWells_[a].size());
        for (int e = 0; e < nv; ++e) J_.rowPtr[a * nv + e + 1] = len;
    }
    for (int w = 0; w < nw; ++w) J_.rowPtr[na * nv + w + 1] = int(wells_[w].rowCells.size()) * nv + 1;
    for (int r = 0; r < N; ++r) J_.rowPtr[r + 1] += J_.rowPtr[r];
    J_.col.resize(J_.rowPtr[N]); J_.val.assign(J_.rowPtr[N], 0.0);
    for (int a = 0; a < na; ++a)
        for (int e = 0; e < nv; ++e) {
            int p = J_.rowPtr[a * nv + e];
            if (!blkAsm_)
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

CSR Simulator::fullJacobian() const {
    const int na = grid_.numActive(), nw = int(wells_.size()), nv = nv_, N = na * nv + nw;
    if (!blkAsm_) return J_;
    CSR F; F.n = N; F.rowPtr.assign(size_t(N) + 1, 0);
    for (int a = 0; a < na; ++a) {
        const int len = (blkStart_[a + 1] - blkStart_[a]) * nv + int(cellWells_[a].size());
        for (int e = 0; e < nv; ++e) F.rowPtr[size_t(a) * nv + e + 1] = len;
    }
    for (int w = 0; w < nw; ++w) F.rowPtr[size_t(na) * nv + w + 1] = J_.rowPtr[size_t(na) * nv + w + 1] - J_.rowPtr[size_t(na) * nv + w];
    for (int r = 0; r < N; ++r) F.rowPtr[size_t(r) + 1] += F.rowPtr[size_t(r)];
    F.col.resize(size_t(F.rowPtr[size_t(N)])); F.val.resize(size_t(F.rowPtr[size_t(N)]));
    for (int a = 0; a < na; ++a)
        for (int e = 0; e < nv; ++e) {
            int p = F.rowPtr[size_t(a) * nv + e];
            for (int k = blkStart_[a]; k < blkStart_[a + 1]; ++k) {
                const double* B = schur_.cellBlock(size_t(k));
                for (int v = 0; v < nv; ++v) { F.col[size_t(p)] = blkCells_[size_t(k)] * nv + v; F.val[size_t(p++)] = B[e * nv + v]; }
            }
            for (int q = J_.rowPtr[size_t(a) * nv + e]; q < J_.rowPtr[size_t(a) * nv + e + 1]; ++q) { F.col[size_t(p)] = J_.col[size_t(q)]; F.val[size_t(p++)] = J_.val[size_t(q)]; }
        }
    for (int r = na * nv; r < N; ++r) {
        int p = F.rowPtr[size_t(r)];
        for (int q = J_.rowPtr[size_t(r)]; q < J_.rowPtr[size_t(r) + 1]; ++q) { F.col[size_t(p)] = J_.col[size_t(q)]; F.val[size_t(p++)] = J_.val[size_t(q)]; }
    }
    F.finalizeDiag();
    return F;
}

// ======================================================================== controls
static bool matchName(const std::string& pat, const std::string& name) {   // '*' / '?' anywhere (C33, as deck.cpp)
    if (pat.find_first_of("*?") == std::string::npos) return pat == name;
    size_t p = 0, n = 0, star = std::string::npos, mark = 0;
    while (n < name.size()) {
        if (p < pat.size() && (pat[p] == '?' || pat[p] == name[n])) { ++p; ++n; }
        else if (p < pat.size() && pat[p] == '*') { star = p++; mark = n; }
        else if (star != std::string::npos) { p = star + 1; n = ++mark; }
        else return false;
    }
    while (p < pat.size() && pat[p] == '*') ++p;
    return p == pat.size();
}

bool mswConnectionShut(const msw::Well& M, int seg);   // sim_msw.cpp

void Simulator::updateConnections(int step) {
    for (size_t wi = 0; wi < wells_.size(); ++wi) {
        WellRT& w = wells_[wi];
        w.nActive = 0;
        const int m = mswOfWell(int(wi));   // C32: connections below a SHUT segment device are closed
        for (size_t c = 0; c < w.cell.size(); ++c) {
            bool on = step >= w.cFrom[c] && (w.cUntil[c] < 0 || step < w.cUntil[c]);
            if (on && m >= 0 && c < mswSegOf_[wi].size() && mswConnectionShut(msw_[size_t(m)], mswSegOf_[wi][c])) on = false;
            w.wi[c] = on ? w.wiBase[c] : 0.0;
            w.nActive += on ? 1 : 0;
        }
    }
}

static double targetOf(const WellControl& c, Control m);   // (defined with the well-switching logic)

// C34 WELPI (ECLIPSE / OPM): the connection factors of the connections present now are multiplied by
// PI_target / PI_current, PI = sum_k WI_k lambda_k of the preferred phase at the current reservoir state (surface rate
// per pressure drop, gas incl. Rs x oil, oil incl. Rv x gas as OPM; injectors: total mobility times the b of the phase injected under this step's controls). Later COMPDAT
// connections are not scaled.
// C46 surface-volume PI per phase of well wi at the current state (OPM computeConnLevelProdInd / InjInd): sum_k WI_k x
// mobility x b over the open connections; producers: oil + Rv x gas (vaporised oil), gas + Rs x oil (dissolved gas);
// injectors: total mobility x b of the injected phase, reported under that phase (water slot for water injectors)
std::array<double, 3> Simulator::wellPIs(size_t wi) const {
    std::array<double, 3> out{0, 0, 0};
    const WellRT& w = wells_[wi];
    const bool prod = w.type == WellType::Producer;
    const int ip = w.type == WellType::GasInjector ? PG : PW;
    for (size_t k = 0; k < w.cell.size(); ++k) {
        if (!(w.wi[k] > 0)) continue;
        const CellFull C = evalCell(w.cell[k]);
        if (!prod) {
            double lt = 0; for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && C.b[p].v > 0) lt += C.mob[p].v / C.b[p].v;
            out[ip == PG ? 2 : 1] += w.wi[k] * lt * (eq_[ip] >= 0 ? C.b[ip].v : 1.0);
            continue;
        }
        const double mo = eq_[PO] >= 0 ? C.mob[PO].v : 0, mg = eq_[PG] >= 0 ? C.mob[PG].v : 0, mw = eq_[PW] >= 0 ? C.mob[PW].v : 0;
        const bool og = eq_[PO] >= 0 && eq_[PG] >= 0;
        const double rs = og ? C.rs.v : 0, rv = og ? C.rv.v : 0;
        out[0] += w.wi[k] * (mo + rv * mg); out[1] += w.wi[k] * mw; out[2] += w.wi[k] * (mg + rs * mo);
    }
    return out;
}

void Simulator::applyWelpi(int step) {
    for (const auto& [pat, target] : case_.schedule[size_t(step)].welpi) {
        if (!(target > 0)) continue;
        for (WellRT& w : wells_) {
            const WellInput& in = case_.wells[size_t(w.input)];
            if (!matchName(pat, in.name)) continue;
            std::string ph = in.prefPhase;
            if (ph.empty()) ph = w.type == WellType::Producer ? "OIL" : w.type == WellType::GasInjector ? "GAS" : "WATER";
            const auto P = wellPIs(size_t(&w - wells_.data()));
            const double pi = w.type != WellType::Producer ? P[w.type == WellType::GasInjector ? 2 : 0] + P[1]
                            : ph == "GAS" ? P[2] : ph == "WATER" ? P[1] : ph == "LIQ" ? P[0] + P[1] : P[0];
            if (!(pi > 0)) { log_("  WELPI " + in.name + ": no flowing connection - ignored"); continue; }
            const double f = target / pi;
            for (size_t k = 0; k < w.cell.size(); ++k) if (w.cDef[k] <= step && k < w.wiBase.size()) { w.wiBase[k] *= f; w.cWelpi[k] *= f; if (w.wi[k] > 0) w.wi[k] = w.wiBase[k]; }
            log_(std::string("  WELPI ") + in.name + fmt(": PI %.6g", pi) + fmt(" -> %.6g", target) + fmt(" (connection factors x %.4g)", f));
        }
    }
}

void Simulator::applySchedMultflt(int step) {
    // C46 MULTFLT in SCHEDULE: every event with event.step <= step not applied yet (restarts start from a fresh grid and
    // re-apply the earlier events). The factor multiplies the current transmissibility of every connection across the
    // fault faces (OPM: the new multiplier is the old one times the factor), incl. NNCs of displaced columns.
    const auto& E = case_.schedMultflt;
    if (E.empty()) return;
    if (multfltDone_.size() != E.size()) multfltDone_.assign(E.size(), 0);
    const int nx = grid_.nx, ny = grid_.ny, nz = grid_.nz, nCells = grid_.numCells();
    auto nameMatch = [](const std::string& pat, const std::string& n) {
        auto up = [](std::string s) { for (char& ch : s) ch = char(std::toupper((unsigned char)ch)); return s; };
        const std::string P = up(pat), N = up(n);
        if (!P.empty() && P.back() == '*') return N.compare(0, P.size() - 1, P, 0, P.size() - 1) == 0;
        return P == N;
    };
    for (size_t e = 0; e < E.size(); ++e) {
        if (multfltDone_[e] || E[e].step > step) continue;
        multfltDone_[e] = 1;
        const double m = std::isfinite(E[e].mult) ? std::max(0.0, E[e].mult) : 1.0;
        std::set<int> conns;
        for (const FaultInput& F : case_.grid.faults) {
            if (!nameMatch(E[e].fault, F.name)) continue;
            for (const FaultFace& ff : F.faces) {
                const int di = ff.face == 'X' ? (ff.minus ? -1 : 1) : 0, dj = ff.face == 'Y' ? (ff.minus ? -1 : 1) : 0, dk = ff.face == 'Z' ? (ff.minus ? -1 : 1) : 0;
                for (int k = std::max(0, ff.k1); k <= std::min(nz - 1, ff.k2); ++k)
                    for (int j = std::max(0, ff.j1); j <= std::min(ny - 1, ff.j2); ++j)
                        for (int i = std::max(0, ff.i1); i <= std::min(nx - 1, ff.i2); ++i) {
                            const int gl = i + nx * (j + ny * k);
                            if (gl >= nCells || size_t(gl) >= grid_.activeOf.size()) continue;
                            const int a = grid_.activeOf[size_t(gl)];
                            if (a < 0) continue;
                            for (int q = grid_.adjStart[size_t(a)]; q < grid_.adjStart[size_t(a) + 1]; ++q) {
                                const int gb = grid_.globalOf[size_t(grid_.adjNbr[size_t(q)])];
                                if (gb >= nCells) continue;
                                const int ib = gb % nx, jb = (gb / nx) % ny, kb = gb / (nx * ny);
                                const bool across = dk != 0 ? (ib == i && jb == j && kb == k + dk) : (ib == i + di && jb == j + dj);
                                if (across) conns.insert(grid_.adjConn[size_t(q)]);
                            }
                        }
            }
        }
        for (int c : conns) {
            grid_.conns[size_t(c)].trans *= m;
            auto it = connT0_.find(c); if (it != connT0_.end()) it->second *= m;
            if (size_t(c) < hfT0_.size() && hfT0_[size_t(c)] >= 0) hfT0_[size_t(c)] *= m;
        }
        log_(fmt("MULTFLT (SCHEDULE) at step %g: x%g on %g connections", double(step), m, double(conns.size())) + " of fault " + E[e].fault);
    }
}

void Simulator::applyControls(int step) {
    applySchedMultflt(step);
    // connections opening / closing at this step: wells without open connections are shut, wells whose
    // first connections open now start producing if their control is OPEN
    mswSetup(step);   // C32 multisegment wells in force at this step
    // C37 WCYCLE: a well opened by the schedule at this report step (status change, new well, injector <-> producer)
    // restarts its on-period (OPM Flow: well_open_times_ at REQUEST_OPEN / WELL_STATUS_CHANGE events)
    std::vector<char> schedOpen0(wells_.size()), reopenEv(wells_.size(), 0); std::vector<WellType> type0(wells_.size());
    for (size_t w = 0; w < wells_.size(); ++w) { schedOpen0[w] = wells_[w].ctrl.status == WellStatus::Open && wells_[w].nActive > 0 && step > 0; type0[w] = wells_[w].type; }
    std::vector<int> before(wells_.size());
    for (size_t w = 0; w < wells_.size(); ++w) before[w] = wells_[w].nActive;
    updateConnections(step);
    for (size_t wi = 0; wi < wells_.size(); ++wi) {
        WellRT& w = wells_[wi];
        if (w.nActive == before[wi] && step > 0) continue;
        const bool want = w.ctrl.status == WellStatus::Open && w.nActive > 0 && w.econ != 1;
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
            if (wc.reopen && wc.status == WellStatus::Open) reopenEv[size_t(&w - wells_.data())] = 1;   // C37 WCYCLE
            // a well stopped by a zero target (e.g. WCONHIST ORAT 0 0 0) keeps a stale BHP: restart it like a newly opened well
            const bool wasOpen = w.open && !(w.mode != Control::BHP && targetOf(w.ctrl, w.mode) <= 0.0 && !w.resvHist);
            w.ctrl = wc; w.ctrl.well = in.name;
            w.open = wc.status == WellStatus::Open && w.nActive > 0;
            if (w.econ == 1) w.open = false;   // C27: shut by an economic limit - only WTEST re-opens it
            w.mode = wc.mode;
            if (w.type != WellType::Producer && wc.phase) w.type = wc.phase == 2 ? WellType::GasInjector : WellType::WaterInjector;
            // C33: a producer record (WCONPROD / WCONHIST rate modes) turns an injector into a producer (OPM / ECLIPSE:
            // the well type follows the last control keyword, e.g. WELSPECS + WCONHIST in an ACTIONX body)
            if (w.type != WellType::Producer && !wc.phase && (wc.history || wc.mode == Control::ORAT || wc.mode == Control::WRAT ||
                                                              wc.mode == Control::GRAT || wc.mode == Control::LRAT || wc.mode == Control::RESV)) {
                w.type = WellType::Producer;
                log_(fmt("Report step %g: well ", step + 1) + in.name + " becomes a producer");
            }
            if (w.type == WellType::Producer && wc.phase)
                w.type = wc.phase == 2 ? WellType::GasInjector : WellType::WaterInjector;   // producer -> injector (WCONINJE)
            if (w.open && !wasOpen) {
                int a0 = w.cell.front();
                for (size_t c = 0; c < w.cell.size(); ++c) if (w.wi[c] > 0) { a0 = w.cell[c]; break; }
                const double pc = p_[a0];
                w.bhp = (w.type == WellType::Producer) ? std::max(pc * 0.9, U_.atm) : pc * 1.05;
            }
            if (w.mode == Control::BHP && wc.bhp > 0) w.bhp = wc.bhp;
            // WCONHIST RESV: the reservoir-voidage target is the observed surface rates converted at field-average
            // conditions (updateResvConversion, every time step - ECLIPSE / OPM practice)
            w.resvHist = wc.history && wc.mode == Control::RESV && wc.resv < 0 && w.type == WellType::Producer && !w.cell.empty();
            if (w.resvHist) w.ctrl.resv = 0;
            w.own = w.ctrl; w.bhpDeck = w.ctrl.bhp; w.grpCtl = false;
            w.rates.oilH = wc.history ? std::max(0.0, wc.orat) : -1;
            w.rates.waterH = wc.history ? std::max(0.0, wc.wrat) : -1;
            w.rates.gasH = wc.history ? std::max(0.0, wc.grat) : -1;
        }
        if (!found) log_("Warning: control for unknown well '" + wc.well + "'");
    }
    // C44: WELPI after the controls of this step (the injected phase / well type of this step decides the injector PI)
    if (step >= 0 && step < int(case_.schedule.size()) && !case_.schedule[size_t(step)].welpi.empty()) applyWelpi(step);
    econControls(step);
    tracerControls(step);
    facilityControls(step);                   // C35 GSATPROD / GCONSUMP
    networkSetup(step);                       // C29 NODEPROP / BRANPROP network of this report step
    applyWefac(step);
    applyWcycle(step);
    for (size_t w = 0; w < wells_.size(); ++w) {
        WellRT& W = wells_[w];
        const bool schedOpen = W.ctrl.status == WellStatus::Open && W.nActive > 0 && W.econ != 1;
        if (!schedOpen) { W.cycOpen = W.cycClose = -1e300; continue; }
        if ((W.open || reopenEv[w]) && (!schedOpen0[w] || type0[w] != W.type || reopenEv[w])) { W.cycOpen = time_; W.cycClose = -1e300; W.cycSched = true; }
    }
    updateResvConversion();
}

// C37 WCYCLE 'well  on-time  off-time  start-up-time  max-step-on  YES|NO' (days). The records in force are re-derived from
// all WCYCLE extra keywords up to this step in deck order (a later record for the same well name replaces the earlier one).
void Simulator::applyWcycle(int step) {
    wcyc_.clear();
    std::vector<const ExtraKeyword*> ks;
    for (const ExtraKeyword& e : case_.extraKeywords)
        if (e.name == "WCYCLE" && e.section == DeckSection::Schedule && e.step <= step) ks.push_back(&e);
    std::stable_sort(ks.begin(), ks.end(), [](const ExtraKeyword* a, const ExtraKeyword* b) { return a->step < b->step; });
    for (const ExtraKeyword* e : ks) {
        std::string clean; bool com = false;
        for (size_t i = 0; i < e->body.size(); ++i) {
            const char ch = e->body[i];
            if (ch == '\n') { com = false; clean += ' '; continue; }
            if (!com && ch == '-' && i + 1 < e->body.size() && e->body[i + 1] == '-') com = true;
            if (!com) clean += ch;
        }
        std::stringstream ss(clean); std::string rec;
        while (std::getline(ss, rec, '/')) {
            // tokens with 'N*' default repeats expanded
            std::vector<std::string> tok; std::stringstream rs(rec); std::string t;
            while (rs >> t) {
                const size_t st = t.find('*');
                if (st != std::string::npos && t.front() != '\'' && t.front() != '"') {
                    int n = 1; if (st > 0) { try { n = std::stoi(t.substr(0, st)); } catch (...) { n = 1; } }
                    for (int q = 0; q < n; ++q) tok.push_back("");
                } else {
                    if (!t.empty() && (t.front() == '\'' || t.front() == '"')) t = t.substr(1, t.size() >= 2 ? t.size() - 2 : 0);
                    tok.push_back(t);
                }
            }
            if (tok.empty() || tok[0].empty()) continue;
            auto num = [&](size_t i) { if (i >= tok.size() || tok[i].empty()) return 0.0; try { return std::stod(tok[i]); } catch (...) { return 0.0; } };
            WCycleRec r; r.well = tok[0]; r.on = num(1); r.off = num(2); r.start = num(3); r.maxDt = num(4);
            std::string c = tok.size() > 5 ? tok[5] : std::string();
            for (char& ch : c) ch = char(std::toupper((unsigned char)ch));
            r.ctrl = c.rfind("Y", 0) == 0;
            bool repl = false;
            for (auto& q : wcyc_) if (q.well == r.well) { q = r; repl = true; }
            if (!repl) wcyc_.push_back(r);
        }
    }
}

// C37 WCYCLE at the start of every time step (OPM Flow WCYCLE::nextTimeStep / wellStatus / efficiencyScale): the step is cut
// to land on the next switch (CONTROLLED_TIMESTEP YES) and limited to MAX_TIMESTEP when a well re-opens; a well open for
// ON_TIME is shut for OFF_TIME and then re-opened; during START_TIME after a cycle re-opening its efficiency factor ramps
// linearly, (t - t_open + dt) / START_TIME. A cycle only switches wells that the schedule keeps open.
void Simulator::wcycleStep(double& dt, bool& clipped) {
    if (wcyc_.empty()) return;
    const double t = time_;
    for (WellRT& W : wells_) W.eff = W.effBase;
    double dtNew = dt;
    for (const WCycleRec& r : wcyc_)
        for (WellRT& W : wells_) {
            if (!matchName(r.well, case_.wells[size_t(W.input)].name)) continue;
            if (r.off > 0 && r.ctrl && W.cycOpen > -1e299) {
                const double tt = W.cycOpen + r.on;
                if (tt > t + 1e-9 && tt < t + dtNew - 1e-9) dtNew = tt - t;
            }
            if (r.on > 0 && W.cycClose > -1e299) {
                const double tt = W.cycClose + r.off;
                if (r.ctrl && tt > t + 1e-9 && tt < t + dtNew - 1e-9) dtNew = tt - t;
                if (t >= tt - 1e-9 && r.maxDt > 0) dtNew = std::min(dtNew, r.maxDt);
            }
        }
    if (dtNew < dt) { dt = dtNew; clipped = true; }
    for (const WCycleRec& r : wcyc_)
        for (WellRT& W : wells_) {
            const std::string& nm = case_.wells[size_t(W.input)].name;
            if (!matchName(r.well, nm)) continue;
            const bool schedOpen = W.ctrl.status == WellStatus::Open && W.nActive > 0 && W.econ != 1;
            int want = -1;   // -1 = no decision, 0 = shut, 1 = open
            bool reopened = false;
            if (r.on > 0 && W.cycClose > -1e299) {
                if (t < W.cycClose + r.off - 1e-9) want = 0;
                else { want = 1; if (r.off > 0) { W.cycClose = -1e300; W.cycOpen = t; W.cycSched = false; reopened = true; } }
            }
            if (r.off > 0) {
                if (W.cycOpen > -1e299) {
                    if (t < W.cycOpen + r.on - 1e-9) { if (want < 0) want = 1; }
                    else { if (want < 0) want = 0; W.cycOpen = -1e300; if (r.on > 0) W.cycClose = t; }
                } else if (r.on < 0 && want < 0) want = 0;
            }
            if (!schedOpen || want < 0) continue;
            if (want == 1 && !W.open) {
                int a0 = W.cell.front();
                for (size_t c = 0; c < W.cell.size(); ++c) if (W.wi[c] > 0) { a0 = W.cell[c]; break; }
                const double pc = p_[size_t(a0)];
                W.bhp = W.type == WellType::Producer ? std::max(pc * 0.9, U_.atm) : pc * 1.05;
                if (W.mode == Control::BHP && W.ctrl.bhp > 0) W.bhp = W.ctrl.bhp;
                W.open = true;
                log_(fmt("t = %.4g d: WCYCLE opens well ", t) + nm);
            } else if (want == 0 && W.open) {
                W.open = false;
                log_(fmt("t = %.4g d: WCYCLE shuts well ", t) + nm);
            }
            (void)reopened;
            if (W.open && r.on > 0 && r.start > 0 && W.cycOpen > -1e299) {
                const double ts = t - W.cycOpen;
                if (ts < r.start - 1e-12 && !W.cycSched)
                    W.eff = W.effBase * (ts + dt > r.start ? 1.0 : (ts + dt) / r.start);
            }
        }
}

// C33 WEFAC (well efficiency factor): 'well-or-pattern  factor  [network flag]'. The well equations work on the flowing
// rates; the reservoir source, group / field rates and all cumulatives use factor x rate (ECLIPSE / OPM Flow). The keyword
// stays an extra keyword of the case (exported verbatim); the factors in force are re-derived from all WEFAC records
// up to this step, in deck order (also the ones an ACTIONX / script added at run time).
void Simulator::applyWefac(int step) {
    for (auto& w : wells_) w.eff = 1.0;
    std::vector<const ExtraKeyword*> ks;
    for (const ExtraKeyword& e : case_.extraKeywords)
        if (e.name == "WEFAC" && e.section == DeckSection::Schedule && e.step <= step) ks.push_back(&e);
    std::stable_sort(ks.begin(), ks.end(), [](const ExtraKeyword* a, const ExtraKeyword* b) { return a->step < b->step; });
    for (const ExtraKeyword* e : ks) {
        std::string body = e->body, rec;
        // records are '/'-terminated; '--' comments to the end of the line
        std::string clean; bool com = false;
        for (size_t i = 0; i < body.size(); ++i) {
            if (body[i] == '\n') { com = false; clean += ' '; continue; }
            if (!com && body[i] == '-' && i + 1 < body.size() && body[i + 1] == '-') com = true;
            if (!com) clean += body[i];
        }
        std::stringstream ss(clean);
        while (std::getline(ss, rec, '/')) {
            std::stringstream rs(rec); std::string nm, fs;
            if (!(rs >> nm)) continue;
            if (!nm.empty() && (nm.front() == '\'' || nm.front() == '"')) nm = nm.substr(1, nm.size() >= 2 ? nm.size() - 2 : 0);
            double f = 1.0;
            if (rs >> fs) {
                if (fs.size() >= 2 && fs.substr(fs.size() - 1) == "*" ) f = 1.0;   // 1* default
                else { try { f = std::stod(fs); } catch (...) { f = 1.0; } }
            }
            f = std::clamp(f, 0.0, 1.0);
            for (auto& w : wells_) if (matchName(nm, case_.wells[size_t(w.input)].name)) w.eff = f;
        }
    }
    for (auto& w : wells_) w.effBase = w.eff;   // C37 WCYCLE ramps on top of WEFAC
}

void Simulator::updateResvConversion() {
    bool any = false;
    for (auto& w : wells_) if (w.type == WellType::Producer && !w.cell.empty()) any = true;
    if (!any) return;
    // field averages weighted by hydrocarbon pore volume (OPM RateConverter::defineState, one region)
    const int na = grid_.numActive();
    double sp = 0, srs = 0, srv = 0, sw8 = 0, spv = 0, spp = 0;
    #pragma omp parallel for reduction(+ : sp, srs, srv, sw8, spv, spp) schedule(static)
    for (int a = 0; a < na; ++a) {
        const double pv = grid_.poreVolRef[size_t(a)];
        const double hc = pv * (1.0 - (iSw_ >= 0 ? sw_[size_t(a)] : 0.0));
        spv += pv; spp += pv * p_[size_t(a)];
        if (hc <= 0) continue;
        sw8 += hc; sp += hc * p_[size_t(a)];
        if (live_ || vap_) { const CellFull C = evalCell(a); srs += hc * C.rs.v; srv += hc * C.rv.v; }
    }
    const double pa = sw8 > 0 ? sp / sw8 : (spv > 0 ? spp / spv : U_.atm);
    const double rsA = sw8 > 0 ? srs / sw8 : 0.0, rvA = sw8 > 0 ? srv / sw8 : 0.0;
    for (auto& w : wells_) {
        w.resvOk = false;
        if (w.type != WellType::Producer || w.cell.empty()) continue;
        const FluidModel& F = fl(w.cell.front());
        // producing ratios: observed rates of a history well, else the last computed rates
        const double qo = w.resvHist ? std::max(0.0, w.rates.oilH) : w.rates.oil, qg = w.resvHist ? std::max(0.0, w.rates.gasH) : w.rates.gas;
        double rs = live_ ? rsA : 0.0, rv = vap_ ? rvA : 0.0;
        if (live_ && qo > 0 && qg >= 0) rs = std::min(rs, qg / qo);
        if (vap_ && qg > 0 && qo >= 0) rv = std::min(rv, qo / qg);
        const double den = std::max(1e-12, 1.0 - rs * rv);
        const double bo = O_ ? F.boS(pa, rs) : 1.0, bw = W_ ? F.bwS(pa) : 1.0;
        const double bg = G_ ? (vap_ ? F.bg(AD3(pa), AD3(rv)).v : F.bgS(pa)) : 1.0;
        w.resvC[0] = W_ ? 1.0 / bw : 0.0;
        w.resvC[1] = (O_ ? 1.0 / (den * bo) : 0.0) - (G_ && live_ ? rs / (den * bg) : 0.0);
        w.resvC[2] = (G_ ? 1.0 / (den * bg) : 0.0) - (O_ && vap_ ? rv / (den * bo) : 0.0);
        if (w.resvHist) {   // the observed ratios cap Rs / Rv, so every reservoir-phase term is >= 0
            const double t = w.resvC[1] * std::max(0.0, w.rates.oilH) + w.resvC[0] * std::max(0.0, w.rates.waterH) + w.resvC[2] * std::max(0.0, w.rates.gasH);
            w.ctrl.resv = std::max(0.0, t); w.own.resv = w.ctrl.resv;
        }
        w.resvOk = w.resvC[0] >= 0 && w.resvC[1] >= 0 && w.resvC[2] >= 0;   // else: in-situ voidage at connection conditions
    }
}

// ======================================================================== groups and THP limits
// ---------------------------------------------------------------- C27 economic limits (WECON) and well tests (WTEST)
// The records stay verbatim extra keywords (exact deck round trip); the engine reads them at their report step.
namespace {
std::vector<std::vector<std::string>> econRecords(const std::string& body) {
    // '/'-terminated records -> tokens; quotes stripped, 'n*' expanded to n empty items, '--' comments removed
    std::vector<std::vector<std::string>> out;
    std::vector<std::string> cur;
    std::string text;
    {   std::istringstream ls(body); std::string line;
        while (std::getline(ls, line)) { const size_t c = line.find("--"); if (c != std::string::npos) line.resize(c); text += line + '\n'; } }
    size_t i = 0;
    while (i < text.size()) {
        const char ch = text[i];
        if (std::isspace(static_cast<unsigned char>(ch))) { ++i; continue; }
        if (ch == '/') { out.push_back(cur); cur.clear(); ++i; continue; }
        std::string tok;
        if (ch == '\'' || ch == '"') { const size_t e = text.find(ch, i + 1); tok = text.substr(i + 1, (e == std::string::npos ? text.size() : e) - i - 1); i = e == std::string::npos ? text.size() : e + 1; cur.push_back(tok); continue; }
        while (i < text.size() && !std::isspace(static_cast<unsigned char>(text[i])) && text[i] != '/') tok += text[i++];
        const size_t st = tok.find('*');
        if (st != std::string::npos && tok.find_first_not_of("0123456789") == st) {
            const int n = st == 0 ? 1 : std::max(1, std::atoi(tok.substr(0, st).c_str()));
            const std::string v = tok.substr(st + 1);
            for (int k = 0; k < n; ++k) cur.push_back(v);
        } else cur.push_back(tok);
    }
    return out;
}
double econNum(const std::vector<std::string>& r, size_t i, double def) {
    if (i >= r.size() || r[i].empty()) return def;
    char* e = nullptr; const double v = std::strtod(r[i].c_str(), &e);
    return (e && e != r[i].c_str()) ? v : def;
}
std::string econStr(const std::vector<std::string>& r, size_t i, const std::string& def) {
    if (i >= r.size() || r[i].empty()) return def;
    std::string s = r[i]; for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c)));
    return s;
}
}  // namespace

std::vector<std::vector<std::string>> verbatimRecords(const std::string& body) { return econRecords(body); }
double verbatimNum(const std::vector<std::string>& r, size_t i, double def) { return econNum(r, i, def); }
std::string verbatimStr(const std::vector<std::string>& r, size_t i, const std::string& def) { return econStr(r, i, def); }

void Simulator::econControls(int step) {
    econLim_.resize(wells_.size()); wtest_.resize(wells_.size());
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.section != DeckSection::Schedule || k.step != step || (k.name != "WECON" && k.name != "WTEST" && k.name != "CECON")) continue;
        if (k.name == "CECON") {   // C41: limits attach to the matching connections of the well (later records override)
            cecon_.resize(wells_.size());
            for (const auto& r : econRecords(k.body)) {
                if (r.empty()) continue;
                for (size_t wi = 0; wi < wells_.size(); ++wi) {
                    WellRT& w = wells_[wi];
                    if (!matchName(r[0], case_.wells[size_t(w.input)].name)) continue;
                    const int I = int(econNum(r, 1, 0)), J = int(econNum(r, 2, 0)), K1 = int(econNum(r, 3, 0)), K2 = int(econNum(r, 4, 0));
                    ConnEcon L;
                    L.maxWct = econNum(r, 5, 0); L.maxGor = econNum(r, 6, 0); L.maxWgr = econNum(r, 7, 0);
                    L.workover = econStr(r, 8, "CON"); L.checkStopped = econStr(r, 9, "NO") == "YES";
                    L.minOil = econNum(r, 10, 0); L.minGas = econNum(r, 11, 0);
                    L.on = L.maxWct > 0 || L.maxGor > 0 || L.maxWgr > 0 || L.minOil > 0 || L.minGas > 0;
                    auto& V = cecon_[wi]; if (V.size() < w.cell.size()) V.resize(w.cell.size());
                    for (size_t c = 0; c < w.cell.size(); ++c) {
                        int i, j, kk; grid_.ijk(grid_.globalOf[size_t(w.cell[c])], i, j, kk);
                        if ((I <= 0 || i + 1 == I) && (J <= 0 || j + 1 == J) && (K1 <= 0 || kk + 1 >= K1) && (K2 <= 0 || kk + 1 <= K2)) V[c] = L;
                    }
                    if (L.on) econAny_ = true;
                }
            }
            continue;
        }
        for (const auto& r : econRecords(k.body)) {
            if (r.empty()) continue;
            for (size_t wi = 0; wi < wells_.size(); ++wi) {
                const std::string& nm = case_.wells[size_t(wells_[wi].input)].name;
                if (!matchName(r[0], nm)) continue;
                if (k.name == "WECON") {
                    EconLimit L;
                    L.minOil = econNum(r, 1, 0); L.minGas = econNum(r, 2, 0); L.maxWct = econNum(r, 3, 0);
                    L.maxGor = econNum(r, 4, 0); L.maxWgr = econNum(r, 5, 0);
                    L.action = econStr(r, 6, "NONE"); L.endRun = econStr(r, 7, "NO") == "YES";
                    L.minLiq = econNum(r, 13, 0);
                    L.maxWct2 = econNum(r, 10, 0); L.action2 = econStr(r, 11, L.action);   // C34 secondary water cut + workover
                    L.potn = econStr(r, 9, "RATE") == "POTN";                               // C34 limits on potentials
                    L.on = L.minOil > 0 || L.minGas > 0 || L.maxWct > 0 || L.maxGor > 0 || L.maxWgr > 0 || L.minLiq > 0;
                    econLim_[wi] = L;
                    if (L.on) econAny_ = true;
                } else {
                    WellTest T; T.interval = econNum(r, 1, 0); T.reasons = econStr(r, 2, ""); T.maxTests = int(econNum(r, 3, 0));
                    wtest_[wi] = T;
                    WellRT& w = wells_[wi];
                    if (w.econ == 1 && T.interval > 0 && T.reasons.find('E') != std::string::npos && w.econNextTest < 0) w.econNextTest = time_ + T.interval;
                    if (T.interval > 0 && T.reasons.find('P') != std::string::npos) econAny_ = true;   // C37 physical closures
                }
            }
        }
    }
}

void Simulator::econCheck() {
    if (!econAny_) return;
    // C37 WTEST 'P' (physical reasons): an open well that cannot flow at its limits (an injector whose BHP limit is below
    // the reservoir pressure, a producer whose BHP limit is above it) is closed and re-tested every WTEST interval
    // (OPM Flow: non-operable wells are shut, WellTestConfig::Reason::PHYSICAL). Without WTEST 'P' it stays open at zero rate.
    for (size_t wi = 0; wi < wells_.size() && wi < wtest_.size(); ++wi) {
        WellRT& w = wells_[wi];
        const WellTest& T = wtest_[wi];
        if (T.interval <= 0 || T.reasons.find('P') == std::string::npos || !w.open || w.econ == 1) continue;
        const bool inj = w.type != WellType::Producer;
        if (!inj && w.mode != Control::BHP && targetOf(w.ctrl, w.mode) <= 0.0) continue;   // stopped by a zero target
        const double q = inj ? std::fabs(w.rates.waterInj) + std::fabs(w.rates.gasInj) : std::fabs(w.rates.oil) + std::fabs(w.rates.water) + std::fabs(w.rates.gas);
        if (q > 1e-6) continue;
        const bool canTest = T.maxTests <= 0 || w.econTests < T.maxTests;
        w.open = false; w.econ = 1; w.econNextTest = canTest ? time_ + T.interval : -1; stats_.econShut++; w.pweSim = w.autoStop ? 4 : 7;
        log_(fmt("t = %.4g d: well ", time_) + case_.wells[size_t(w.input)].name + " cannot flow at its limits - shut (WTEST P" +
             (canTest ? fmt(", next test at %.4g d)", w.econNextTest) : std::string(", no tests left)")));
    }
    for (size_t wi = 0; wi < wells_.size() && wi < econLim_.size(); ++wi) {
        WellRT& w = wells_[wi];
        const EconLimit& L = econLim_[wi];
        if (!L.on || !w.open || w.type != WellType::Producer) continue;
        double qo = w.rates.oil, qw = w.rates.water, qg = w.rates.gas, ql = qo + qw;
        std::string why;
        bool rateLimit = false;
        char b[160];
        int ratio = 0;   // violated ratio limit: 1 WCT, 2 GOR, 3 WGR, 4 secondary WCT
        bool noConns = false;   // C44: connection workover left no open connection (the well is then shut, not stopped)
        const double qoP = L.potn ? wellPotential(w, Control::ORAT) : qo, qgP = L.potn ? wellPotential(w, Control::GRAT) : qg,
                     qlP = L.potn ? wellPotential(w, Control::LRAT) : ql;
        if (L.minOil > 0 && qoP < L.minOil) { std::snprintf(b, sizeof b, "oil rate %.4g < %.4g", qoP, L.minOil); why = b; rateLimit = true; }
        else if (L.minGas > 0 && qgP < L.minGas) { std::snprintf(b, sizeof b, "gas rate %.4g < %.4g", qgP, L.minGas); why = b; rateLimit = true; }
        else if (L.minLiq > 0 && qlP < L.minLiq) { std::snprintf(b, sizeof b, "liquid rate %.4g < %.4g", qlP, L.minLiq); why = b; rateLimit = true; }
        else if (L.maxWct2 > 0 && ql > 0 && qw / ql > L.maxWct2) { std::snprintf(b, sizeof b, "water cut %.4g > secondary %.4g", qw / ql, L.maxWct2); why = b; ratio = 4; }
        else if (L.maxWct > 0 && ql > 0 && qw / ql > L.maxWct) { std::snprintf(b, sizeof b, "water cut %.4g > %.4g", qw / ql, L.maxWct); why = b; ratio = 1; }
        else if (L.maxGor > 0 && qo > 0 && qg / qo > L.maxGor) { std::snprintf(b, sizeof b, "GOR %.4g > %.4g", qg / qo, L.maxGor); why = b; ratio = 2; }
        else if (L.maxWgr > 0 && qg > 0 && qw / qg > L.maxWgr) { std::snprintf(b, sizeof b, "WGR %.4g > %.4g", qw / qg, L.maxWgr); why = b; ratio = 3; }
        if (why.empty()) { if (w.econ == 2) w.econ = 0; continue; }   // a tested well that stays inside its limits is back
        const std::string& nm = case_.wells[size_t(w.input)].name;
        const std::string action = ratio == 4 ? L.action2 : L.action;
        // C34 connection workovers (ECLIPSE WECON item 7): CON closes the worst-offending connection (or COMPLUMP completion),
        // +CON also every connection below it, PLUG plugs back (the worst and all deeper); the well stays open while it
        // has open connections. Ratios per completion from the connection rates of the converged state.
        if (!rateLimit && (action == "CON" || action == "+CON" || action == "PLUG")) {
            const std::vector<int> lump = completionLumps(wi);
            std::map<int, std::array<double, 4>> agg;   // lump -> oil, water, gas, depth (max)
            const auto cr = connRatesOf(wi);
            for (size_t k = 0; k < w.cell.size(); ++k) {
                if (!(w.wi[k] > 0) || k >= cr.size()) continue;
                auto& A = agg.emplace(lump[k], std::array<double, 4>{0, 0, 0, -1e30}).first->second;
                A[0] += cr[k][0]; A[1] += cr[k][1]; A[2] += cr[k][2]; A[3] = std::max(A[3], w.cdepth[k]);
            }
            int worst = -1; double wv = -1, wdepth = 0;
            for (auto& [id, A] : agg) {
                const double r = ratio == 2 ? (A[0] > 0 ? A[2] / A[0] : (A[2] > 0 ? 1e30 : 0)) : ratio == 3 ? (A[2] > 0 ? A[1] / A[2] : (A[1] > 0 ? 1e30 : 0))
                                                                     : (A[0] + A[1] > 0 ? A[1] / (A[0] + A[1]) : 0);
                if (r > wv) { wv = r; worst = id; wdepth = A[3]; }
            }
            int closed = 0;
            for (size_t k = 0; k < w.cell.size() && worst >= 0; ++k) {
                if (!(w.wi[k] > 0)) continue;
                const bool hit = lump[k] == worst || (action != "CON" && w.cdepth[k] >= wdepth - 1e-9);
                if (!hit) continue;
                if (lump[k] != worst) { w.pweSwept.insert(w.cell[k]); w.pweReach = true; }   // C43 WPWE: reach, not an offender
                w.econConns.emplace(w.cell[k], w.cUntil[k]);
                w.wi[k] = 0; w.cUntil[k] = w.cUntil[k] < 0 ? step_ : std::min(w.cUntil[k], step_); --w.nActive; ++closed;
            }
            log_(fmt("t = %.4g d: well ", time_) + nm + " economic limit (" + why + "): workover " + action + fmt(" closes %g connection(s)", closed));
            if (w.nActive > 0) { if (L.endRun) econStop_ = true; continue; }
            log_("  no open connections left - the well is shut");
            noConns = true;
        }
        const bool act = rateLimit || action != "NONE";
        if (act) {
            w.open = false; w.econ = 1; stats_.econShut++; w.pweSim = w.autoStop && !noConns ? 4 : 7;
            const WellTest& T = wi < wtest_.size() ? wtest_[wi] : WellTest{};
            const bool canTest = T.interval > 0 && T.reasons.find('E') != std::string::npos && (T.maxTests <= 0 || w.econTests < T.maxTests);
            w.econNextTest = canTest ? time_ + T.interval : -1;
            log_(fmt("t = %.4g d: well ", time_) + nm + " shut by economic limit (" + why + ")" +
                 (rateLimit || action == "WELL" || action == "CON" || action == "+CON" || action == "PLUG" ? "" : " - workover " + action) +
                 (canTest ? fmt(", WTEST in %.4g d", T.interval) : ""));
        }
        if (L.endRun) { econStop_ = true; log_(fmt("t = %.4g d: economic limit of well ", time_) + nm + " ends the run (WECON item 8 = YES)"); }
    }
    ceconCheck();
}

// C41 CECON (OPM WellTest::updateWellTestStateCECON): after the WECON checks, every open completion (COMPLUMP lump or
// single connection) of a producer that carries limits is tested on its own summed surface rates. Water cut
// qw / (qo + qw) (1 when the oil rate is negative), GOR (infinite without oil), WGR (infinite without gas). CON closes the
// completion, +CON also every completion after it in the connection order, WELL shuts the well; a well whose completions
// are all closed is shut. Stopped wells (zero target) are skipped unless item 10 = YES.
// C42 extensions beyond OPM (which ignores them): items 11 / 12 minimum oil / gas rate of the completion - a producing
// completion below the limit offends and is closed (CON, whatever item 9 says: a rate limit is not a ratio violation that
// spreads with depth). PLUG = plug-back by depth: a water violation (WCT / WGR) closes the offending completion and every
// deeper one, a gas violation (GOR) the offending one and every shallower one. The WPLUG plug-back lengths are not used.
void Simulator::ceconCheck() {
    for (size_t wi = 0; wi < wells_.size() && wi < cecon_.size(); ++wi) {
        WellRT& w = wells_[wi];
        auto& V = cecon_[wi];
        if (V.empty() || !w.open || w.type != WellType::Producer) continue;
        const bool stopped = w.mode != Control::BHP && targetOf(w.ctrl, w.mode) <= 0.0;
        const std::vector<int> lump = completionLumps(wi);
        const auto cr = connRatesOf(wi);
        const std::string& nm = case_.wells[size_t(w.input)].name;
        std::set<int> done;
        for (size_t k = 0; k < w.cell.size() && k < V.size(); ++k) {
            const ConnEcon& L = V[k];
            if (!L.on || !(w.wi[k] > 0) || done.count(lump[k])) continue;
            if (stopped && !L.checkStopped) continue;
            done.insert(lump[k]);
            double qo = 0, qw = 0, qg = 0;
            for (size_t m = 0; m < w.cell.size() && m < cr.size(); ++m)
                if (lump[m] == lump[k] && w.wi[m] > 0) { qo += cr[m][0]; qw += cr[m][1]; qg += cr[m][2]; }
            char b[160]; b[0] = 0;
            if (L.maxWct > 0 && qo + qw > 0) { const double wc = qo >= 0 ? qw / (qo + qw) : 1.0; if (wc > L.maxWct) std::snprintf(b, sizeof b, "water cut %.4g > %.4g", wc, L.maxWct); }
            if (!b[0] && L.maxGor > 0 && qg > 0 && (qo <= 0 || qg / qo > L.maxGor)) std::snprintf(b, sizeof b, "GOR %.4g > %.4g", qo > 0 ? qg / qo : 1e30, L.maxGor);
            if (!b[0] && L.maxWgr > 0 && qw > 0 && (qg <= 0 || qw / qg > L.maxWgr)) std::snprintf(b, sizeof b, "WGR %.4g > %.4g", qg > 0 ? qw / qg : 1e30, L.maxWgr);
            int cause = b[0] ? (b[0] == 'G' ? 2 : 1) : 0;   // 1 water, 2 gas, 3 minimum rate
            if (!b[0] && L.minOil > 0 && qo + qw + qg > 0 && qo < L.minOil) { std::snprintf(b, sizeof b, "oil rate %.4g < %.4g", qo, L.minOil); cause = 3; }
            if (!b[0] && L.minGas > 0 && qo + qw + qg > 0 && qg < L.minGas) { std::snprintf(b, sizeof b, "gas rate %.4g < %.4g", qg, L.minGas); cause = 3; }
            if (!b[0]) continue;
            const std::string wo = cause == 3 ? std::string("CON") : L.workover;
            if (wo == "WELL") {
                w.open = false; w.econ = 1; stats_.econShut++; w.pweSim = w.autoStop ? 4 : 7;
                const WellTest& T = wi < wtest_.size() ? wtest_[wi] : WellTest{};
                const bool canTest = T.interval > 0 && T.reasons.find('E') != std::string::npos && (T.maxTests <= 0 || w.econTests < T.maxTests);
                w.econNextTest = canTest ? time_ + T.interval : -1;
                log_(fmt("t = %.4g d: well ", time_) + nm + " shut by the CECON limit of a completion (" + b + ")" + (canTest ? fmt(", WTEST in %.4g d", T.interval) : std::string()));
                break;
            }
            std::set<int> shut = {lump[k]};
            if (wo == "+CON") {
                bool below = false;
                for (size_t m = 0; m < w.cell.size(); ++m) { if (lump[m] == lump[k]) { below = true; continue; } if (below) shut.insert(lump[m]); }
            } else if (wo == "PLUG" && k < w.cdepth.size()) {
                const double d0 = w.cdepth[k];
                for (size_t m = 0; m < w.cell.size() && m < w.cdepth.size(); ++m)
                    if (cause == 2 ? w.cdepth[m] <= d0 : w.cdepth[m] >= d0) shut.insert(lump[m]);
            }
            int closed = 0;
            for (size_t m = 0; m < w.cell.size(); ++m) {
                if (!shut.count(lump[m]) || !(w.wi[m] > 0)) continue;
                if (lump[m] != lump[k]) { w.pweSwept.insert(w.cell[m]); w.pweReach = true; }   // C43 WPWE: '+CON' / PLUG reach
                w.econConns.emplace(w.cell[m], w.cUntil[m]);
                w.wi[m] = 0; w.cUntil[m] = w.cUntil[m] < 0 ? step_ : std::min(w.cUntil[m], step_); --w.nActive; ++closed;
            }
            log_(fmt("t = %.4g d: well ", time_) + nm + " CECON limit (" + b + "): workover " + wo + fmt(" closes %g connection(s)", closed));
            if (w.nActive <= 0) {   // no open connection left: shut whatever WELSPECS item 9 says (OPM)
                w.open = false; w.econ = 1; stats_.econShut++; w.pweSim = 7;
                const WellTest& T = wi < wtest_.size() ? wtest_[wi] : WellTest{};
                const bool canTest = T.interval > 0 && T.reasons.find('E') != std::string::npos && (T.maxTests <= 0 || w.econTests < T.maxTests);
                w.econNextTest = canTest ? time_ + T.interval : -1;
                log_("  no open connections left - the well is shut");
                break;
            }
        }
    }
}

// C34 COMPLUMP: completion number per connection of well wi (the latest records up to this report step; unlumped
// connections are their own completion, numbered after the lumps)
std::vector<int> Simulator::completionLumps(size_t wi) const {
    const WellRT& w = wells_[wi];
    const std::string& nm = case_.wells[size_t(w.input)].name;
    std::vector<int> lump(w.cell.size(), -1);
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.section != DeckSection::Schedule || k.step > step_ || k.name != "COMPLUMP") continue;
        for (const auto& r : verbatimRecords(k.body)) {
            if (r.size() < 6 || !matchName(r[0], nm)) continue;
            auto iv = [&](size_t i) { return r[i].empty() || r[i] == "*" ? 0 : std::atoi(r[i].c_str()); };
            const int I = iv(1), J = iv(2), K1 = iv(3), K2 = iv(4), N = iv(5);
            for (size_t c = 0; c < w.cell.size(); ++c) {
                int i, j, kk; grid_.ijk(grid_.globalOf[size_t(w.cell[c])], i, j, kk);
                if ((I <= 0 || i + 1 == I) && (J <= 0 || j + 1 == J) && (K1 <= 0 || kk + 1 >= K1) && (K2 <= 0 || kk + 1 <= K2)) lump[c] = N;
            }
        }
    }
    int next = 100000;
    for (int& l : lump) if (l < 0) l = next++;
    return lump;
}

// C34 surface rates (oil incl. vaporised, water, gas incl. dissolved) of each connection of producer wi at the current state
std::vector<std::array<double, 3>> Simulator::connRatesOf(size_t wi) const {
    const WellRT& w = wells_[wi];
    std::vector<std::array<double, 3>> out(w.cell.size(), std::array<double, 3>{0, 0, 0});
    const double g = U_.gravity;
    const int refP = O_ ? PO : PW;
    for (size_t k = 0; k < w.cell.size(); ++k) {
        if (!(w.wi[k] > 0)) continue;
        const CellFull C = evalCell(w.cell[k]);
        const double dd = C.P[refP].v - w.pConn(k, g);
        if (dd <= 0) continue;
        double q[3] = {0, 0, 0};
        for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) q[p] = w.wi[k] * C.mob[p].v * dd;
        const double qd = live_ ? w.wi[k] * C.mobRs.v * dd : 0.0, qv = vap_ ? w.wi[k] * C.mobRv.v * dd : 0.0;
        out[k] = {q[PO] + qv, q[PW], q[PG] + qd};
    }
    return out;
}

// C43 WPWE0-7 (OPM data::WellPerformanceEvents, opm-simulators PR 7395): a snapshot comparison between consecutive
// accepted time steps. Connections are identified by their cell (LGR cells have their own active index). Opened
// connections count only while the well is open (neither shut nor stopped); closed connections count whatever closed them,
// except those closed only by the reach of a '+CON' (or ResSim PLUG) workover - just the offender counts. Closed to bottom:
// a reach past the offender, or this step's closures leave the well flowing over its topmost completion alone or not at
// all. Status events only for simulator decisions (economic limits, failed WTEST); deck SHUT / STOP and WCYCLE are not
// events. Economic closures run after the step is accepted, so retried attempts cannot report or lose events.
void Simulator::pweBaseline(WellRT& w) {
    w.pweOpenCells.clear();
    for (size_t k = 0; k < w.cell.size() && k < w.wi.size(); ++k) if (w.wi[k] > 0) w.pweOpenCells.push_back(w.cell[k]);
    std::sort(w.pweOpenCells.begin(), w.pweOpenCells.end());
    w.pweInj = w.type != WellType::Producer;
    w.pweSwept.clear(); w.pweReach = false; w.pweSim = 0; w.pweBase = true;
}

bool Simulator::pweExists(const WellRT& w) const {   // the well has a connection in the schedule up to this report step
    for (size_t k = 0; k < w.cFrom.size(); ++k) if (w.cFrom[k] <= step_) return true;
    return false;
}

void Simulator::pweTrack() {
    for (size_t wi = 0; wi < wells_.size(); ++wi) {
        WellRT& w = wells_[wi];
        int ev[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        if (w.pweBase) {
            std::vector<int> now;
            for (size_t k = 0; k < w.cell.size() && k < w.wi.size(); ++k) if (w.wi[k] > 0) now.push_back(w.cell[k]);
            std::sort(now.begin(), now.end());
            const bool flowing = w.open && w.ctrl.status == WellStatus::Open;
            int opened = 0, closed = 0, closedAny = 0;
            for (int c : now) if (!std::binary_search(w.pweOpenCells.begin(), w.pweOpenCells.end(), c)) ++opened;
            for (int c : w.pweOpenCells)
                if (!std::binary_search(now.begin(), now.end(), c)) { ++closedAny; if (!w.pweSwept.count(c)) ++closed; }
            ev[1] = flowing ? opened : 0;
            ev[2] = closed;
            bool bottom = w.pweReach;
            if (!bottom && closedAny > 0) {   // flowing over the topmost completion alone, or not at all
                // topmost completion = the well's first deck connection (also when it is shut and so has no open
                // connection here); with COMPLUMP its whole lump. Without completion indices: the shallowest connection.
                const std::vector<int> lump = completionLumps(wi);
                int top = -1;
                if (w.cComp.size() == w.cell.size() && !w.cComp.empty()) {
                    int cmin = INT_MAX; for (int cc : w.cComp) cmin = std::min(cmin, cc);
                    // (a first connection in an inactive cell does not exist for the run: then the first one present counts)
                    const auto& comps = case_.wells[size_t(w.input)].completions;
                    const bool firstPresent = cmin == 0 || (!comps.empty() && comps[0].open);
                    for (size_t k = 0; k < w.cell.size() && k < lump.size(); ++k) if (firstPresent && w.cComp[k] == cmin) { top = lump[k]; break; }
                    if (!firstPresent) top = INT_MIN;   // the first connection is not in the run: only "no open connection" counts
                } else {
                    double dmin = 1e300;
                    for (size_t k = 0; k < w.cell.size() && k < lump.size(); ++k) {
                        const double d = k < w.cdepth.size() ? w.cdepth[k] : double(k);
                        if (d < dmin) { dmin = d; top = lump[k]; }
                    }
                }
                bool onlyTop = true;
                for (size_t k = 0; k < w.cell.size() && k < w.wi.size() && k < lump.size(); ++k) if (w.wi[k] > 0 && lump[k] != top) onlyTop = false;
                bottom = onlyTop;
            }
            ev[3] = bottom ? 1 : 0;
            if (w.pweSim == 4) ev[4] = 1;
            if (w.pweSim == 7) ev[7] = 1;
            const bool inj = w.type != WellType::Producer;
            ev[5] = w.pweInj && !inj ? 1 : 0;
            ev[6] = !w.pweInj && inj ? 1 : 0;
        }
        for (int e = 0; e < 8; ++e) w.pweLast[e] = ev[e];
        if (!summary_.empty() && wi < summary_.back().wells.size())
            for (int e = 0; e < 8; ++e) summary_.back().wells[wi].pwe[e] = ev[e];
        if (w.pweBase || pweExists(w)) pweBaseline(w);
    }
}

void Simulator::econTest() {
    if (!econAny_) return;
    for (WellRT& w : wells_) {
        if (w.econ != 1 || w.econNextTest < 0 || time_ < w.econNextTest - 1e-9) continue;
        w.econNextTest = -1;
        // C43 (OPM WTEST 'E'): a well left without open connections by CON / +CON workovers is retested with those
        // connections re-opened (deck closures that fell due meanwhile stay closed)
        if (w.ctrl.status == WellStatus::Open && w.nActive <= 0 && !w.econConns.empty()) {
            int re = 0;
            for (size_t k = 0; k < w.cell.size(); ++k) {
                auto it = w.econConns.find(w.cell[k]);
                if (it == w.econConns.end() || w.wi[k] > 0) continue;
                if (step_ < w.cFrom[k] || (it->second >= 0 && step_ >= it->second)) continue;
                w.cUntil[k] = it->second; w.wi[k] = w.wiBase[k]; ++w.nActive; ++re;
            }
            w.econConns.clear();
            if (re) log_(fmt("t = %.4g d: WTEST re-opens %g connection(s) closed by economic workovers of well ", time_, re) + case_.wells[size_t(w.input)].name);
        }
        if (w.ctrl.status != WellStatus::Open || w.nActive <= 0) continue;
        int a0 = w.cell.front();
        for (size_t c = 0; c < w.cell.size(); ++c) if (w.wi[c] > 0) { a0 = w.cell[c]; break; }
        w.bhp = std::max(p_[size_t(a0)] * 0.9, U_.atm);
        if (w.mode == Control::BHP && w.ctrl.bhp > 0) w.bhp = w.ctrl.bhp;
        w.open = true; w.econ = 2; ++w.econTests; stats_.econTests++;
        w.cycOpen = time_; w.cycClose = -1e300; w.cycSched = true;   // C37 WCYCLE: a WTEST re-opening starts an on-period
        log_(fmt("t = %.4g d: WTEST re-opens well ", time_) + case_.wells[size_t(w.input)].name + fmt(" (test %g)", w.econTests));
    }
}

void Simulator::groupControls(int step) {
    if (step == 0 || wellChain_.size() != wells_.size()) {
        gParent_.clear();
        for (auto& [ch, pa] : case_.groupTree) gParent_[ch] = pa;
        wellChain_.assign(wells_.size(), {});
        for (size_t i = 0; i < wells_.size(); ++i) {
            std::string g = case_.wells[size_t(wells_[i].input)].group;
            for (int d = 0; d < 32 && !g.empty(); ++d) {
                wellChain_[i].push_back(g);
                if (g == "FIELD") break;
                auto it = gParent_.find(g);
                g = it == gParent_.end() ? std::string("FIELD") : it->second;
            }
        }
    }
    const ReportStep& S = case_.schedule[size_t(step)];
    // C34 GUIDERAT in force: the latest record up to this step (also after a restart)
    guide_ = GuideRateSpec();
    for (int s = 0; s <= step && s < int(case_.schedule.size()); ++s) if (case_.schedule[size_t(s)].guiderat.set) guide_ = case_.schedule[size_t(s)].guiderat;
    if (S.guiderat.set) { for (auto& w : wells_) { w.guide = -1; w.guideT = -1e30; } log_("  guide-rate formula GUIDERAT phase " + guide_.phase); }
    for (const GroupProdControl& g : S.gprod) {
        const bool any = g.mode != Control::BHP || g.orat > 0 || g.wrat > 0 || g.grat > 0 || g.lrat > 0 || g.resv > 0 || !g.uda.empty();
        if (any) gprod_[g.group] = g; else gprod_.erase(g.group);
        log_(std::string("  group ") + g.group + (any ? std::string(" production control ") + (g.mode == Control::BHP ? "NONE (limits)" : toString(g.mode)) : std::string(" production control removed")));
    }
    for (const GroupInjControl& g : S.ginj) {
        auto& m = ginj_[g.phase == 2 ? 1 : 0];
        if (g.mode == "NONE" && !(g.rate > 0) && !(g.resv > 0) && g.uda.empty()) m.erase(g.group); else m[g.group] = g;
        log_(std::string("  group ") + g.group + (g.phase == 2 ? " gas" : " water") + " injection control " + g.mode);
    }
    groupsOn_ = !gprod_.empty() || !ginj_[0].empty() || !ginj_[1].empty();
}

std::vector<Simulator::ConnRate> Simulator::connectionRates() const {
    std::vector<ConnRate> out;
    const double g = U_.gravity;
    const int refP = O_ ? PO : PW;
    for (const WellRT& w : wells_) {
        if (!w.open) continue;
        const bool prod = w.type == WellType::Producer;
        for (size_t k = 0; k < w.cell.size(); ++k) {
            if (w.wi[k] <= 0) continue;
            const int a = w.cell[k];
            const CellFull C = evalCell(a);
            ConnRate r; r.well = case_.wells[size_t(w.input)].name;
            grid_.ijk(grid_.globalOf[size_t(a)], r.i, r.j, r.k);
            r.pBlock = C.P[refP].v; r.pConn = w.pConn(k, g); r.wi = w.wi[k];
            r.sw = iSw_ >= 0 ? sw_[size_t(a)] : 0.0; r.sg = C.S[PG].v;
            if (prod) {
                const double dd = C.P[refP].v - r.pConn;
                if (dd > 0) {
                    double q[3] = {0, 0, 0};
                    for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) q[p] = w.wi[k] * C.mob[p].v * dd;
                    if (eq_[PW] >= 0 && k < w.wmulConn.size()) q[PW] *= w.wmulConn[k];
                    const double qd = live_ ? w.wi[k] * C.mobRs.v * dd : 0.0, qv = vap_ ? w.wi[k] * C.mobRv.v * dd : 0.0;
                    r.oil = q[PO] + qv; r.water = q[PW]; r.gas = q[PG] + qd;
                    for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && C.b[p].v > 0) r.resv += q[p] / C.b[p].v;
                }
            } else {
                const int ip = w.type == WellType::GasInjector ? PG : PW;
                const double dd = r.pConn - C.P[ip].v;
                if (dd > 0 && eq_[ip] >= 0) {
                    double lt = 0;
                    for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && C.b[p].v > 0) lt += C.mob[p].v / C.b[p].v * (p == PW && k < w.wmulConn.size() ? w.wmulConn[k] : 1.0);
                    const double q = w.wi[k] * lt * C.b[ip].v * dd;
                    (ip == PW ? r.water : r.gas) = -q; r.resv = -q / C.b[ip].v;
                } else if (dd < 0 && w.xflow && (mswOfWell(int(&w - wells_.data())) < 0 || mswXflow())) {   // C41 injector crossflow: the connection produces (as in the assembly; C45 also MSW)
                    double q[3] = {0, 0, 0};
                    for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) q[p] = -w.wi[k] * C.mob[p].v * dd;
                    if (eq_[PW] >= 0 && k < w.wmulConn.size()) q[PW] *= w.wmulConn[k];
                    const double qd = live_ ? -w.wi[k] * C.mobRs.v * dd : 0.0, qv = vap_ ? -w.wi[k] * C.mobRv.v * dd : 0.0;
                    r.oil = q[PO] + qv; r.water = q[PW]; r.gas = q[PG] + qd;
                    for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && C.b[p].v > 0) r.resv += q[p] / C.b[p].v;
                }
            }
            out.push_back(r);
        }
    }
    return out;
}

double Simulator::wellPotential(const WellRT& w, Control m) const {
    const bool prod = w.type == WellType::Producer;
    double lim = w.ctrl.bhp > 0 ? w.ctrl.bhp : (prod ? U_.atm : -1);
    const double g = U_.gravity;
    double qo = 0, qw = 0, qg = 0, resv = 0, qi = 0;
    for (size_t k = 0; k < w.cell.size(); ++k) {
        if (w.wi[k] <= 0) continue;
        const int a = w.cell[k];
        const CellFull C = evalCell(a);
        if (!prod && lim <= 0) lim = 1.5 * p_[a];
        const double pc = lim + (k < w.dpConn.size() ? w.dpConn[k] : w.headRho * g * (w.cdepth[k] - w.refDepth));
        if (prod) {
            double q[3] = {0, 0, 0};
            for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) q[p] = w.wi[k] * C.mob[p].v * std::max(0.0, C.P[p].v - pc);
            qw += q[PW]; qo += q[PO] + C.rv.v * q[PG]; qg += q[PG] + C.rs.v * q[PO];
            if (w.resvOk) resv += w.resvC[0] * q[PW] + w.resvC[1] * (q[PO] + C.rv.v * q[PG]) + w.resvC[2] * (q[PG] + C.rs.v * q[PO]);
            else for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && C.b[p].v > 0) resv += q[p] / C.b[p].v;
        } else {
            const int ip = w.type == WellType::GasInjector ? PG : PW;
            if (eq_[ip] < 0) continue;
            double lt = 0;
            for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && C.b[p].v > 0) lt += C.mob[p].v / C.b[p].v;
            qi += w.wi[k] * lt * C.b[ip].v * std::max(0.0, pc - C.P[ip].v);
        }
    }
    switch (m) {
        case Control::ORAT: return qo; case Control::WRAT: return qw; case Control::GRAT: return qg;
        case Control::LRAT: return qo + qw; case Control::RESV: return resv;
        case Control::BHP: return qo + qw + qg + resv;   // (unused)
        default: return qi;
    }
}

// C34 GUIDERAT (ECLIPSE / OPM): GR = POT^A / (B + C R1^D + E R2^F) for the nominated phase, R1 / R2 = water cut and GOR
// (OIL, LIQ, RES, WAT) or water-gas and oil-gas ratios (GAS) of the well's current production; recomputed when the
// minimum interval has passed, damped GR = d GR_new + (1 - d) GR_old, optionally not allowed to increase. For a group
// controlled on another quantity the guide rate is converted with the well's production ratio q_m / q_nominated.
double Simulator::guideRateOf(WellRT& w, Control m) {
    const std::string& ph = guide_.phase;
    const Control nom = ph == "OIL" ? Control::ORAT : ph == "GAS" ? Control::GRAT : ph == "WAT" ? Control::WRAT : ph == "RES" ? Control::RESV : Control::LRAT;
    const WellRates& r = w.rates;
    auto q = [&](Control c) {
        switch (c) { case Control::ORAT: return r.oil; case Control::WRAT: return r.water; case Control::GRAT: return r.gas;
                     case Control::LRAT: return r.oil + r.water; case Control::RESV: return std::max(0.0, r.resv); default: return 0.0; }
    };
    const bool flowing = r.oil + r.water + r.gas > 0;
    // C44 OPM GuideRate::get: a well without a guide-rate value yet (OPM evaluates the formula after the step, from that
    // step's potentials) gets its POTENTIAL of the controlled quantity - the first step of the run / of a new well
    // (opm-tests model2 GUIDERAT_DEFAULT_WGRUPCON: day-1 shares follow the potentials, not POT^3)
    if (w.guide < 0 && !flowing) return std::max(0.0, wellPotential(w, m));
    if (w.guide < 0 || time_ - w.guideT >= std::max(guide_.minInterval, 1e-6)) {
        // OPM GuideRateModel::eval: the potential and both ratios come from the phase POTENTIALS
        // (OIL / RES: R1 = WOR, R2 = GOR; LIQ: R1 = WCT, R2 = gas / liquid; GAS: R1 = WGR, R2 = OGR)
        const double po = std::max(0.0, wellPotential(w, Control::ORAT)), pwat = std::max(0.0, wellPotential(w, Control::WRAT)),
                     pg = std::max(0.0, wellPotential(w, Control::GRAT));
        double pot = ph == "OIL" ? po : ph == "GAS" ? pg : ph == "WAT" ? pwat : ph == "RES" ? po + pwat + pg : po + pwat;
        double r1 = 0, r2 = 0;
        if (ph == "GAS") { r1 = pg > 0 ? pwat / pg : 0; r2 = pg > 0 ? po / pg : 0; }
        else if (ph == "LIQ") { r1 = po + pwat > 0 ? pwat / (po + pwat) : 0; r2 = po + pwat > 0 ? pg / (po + pwat) : 0; }
        else { r1 = po > 0 ? pwat / po : 0; r2 = po > 0 ? pg / po : 0; }
        auto pw = [](double x, double e) { return e == 0 ? 1.0 : (x > 0 ? std::pow(x, e) : 0.0); };
        const double den = guide_.B + guide_.C * pw(r1, guide_.D) + guide_.E * pw(r2, guide_.F);
        double gr = pot <= 0 || den <= 0 ? 0.0 : pw(pot, guide_.A) / den;
        if (w.guide >= 0) {
            if (!guide_.allowIncrease) gr = std::min(gr, w.guide);
            const double d = std::clamp(guide_.damping, 0.0, 1.0);
            gr = d * gr + (1 - d) * w.guide;
        }
        w.guide = gr; w.guideT = time_;
        if (std::getenv("RESSIM_THPDBG")) std::fprintf(stderr, "GR %s t %.2f pot %.2f r1 %.4f r2 %.4f gr %.3f\n", case_.wells[size_t(w.input)].name.c_str(), time_, pot, r1, r2, gr);
    }
    if (m == nom) return w.guide;
    // conversion to the controlled quantity with the well's current production ratio (OPM GuideRate::get): when the
    // well has no rate in the guide-rate phase yet (first step) OPM returns the guide rate unscaled - C41: the earlier
    // fallback to the potential ratio gave high-water-cut wells a 2-3x larger first share (model2 GRCTRL STW 15-22 %)
    const double qn = ph == "RES" ? q(Control::ORAT) + q(Control::WRAT) + q(Control::GRAT) : q(nom), qm = q(m);
    if (!flowing || !(qn > 1e-9)) return w.guide;
    return w.guide * qm / qn;
}

void Simulator::applyGroupTargets() {
    if (!groupsOn_) return;
    // wells under group control in the previous time step return to their own records
    for (auto& w : wells_) {
        if (!w.grpCtl) continue;
        const double b = w.ctrl.bhp;
        w.ctrl = w.own; w.ctrl.bhp = b;
        if (w.mode != Control::BHP) w.mode = w.own.mode;
        w.grpCtl = false;
    }
    // well potentials (rates at the BHP / THP limit at the current reservoir state), cached for this time step
    std::vector<std::map<int, double>> potC(wells_.size());
    auto pot = [&](size_t i, Control m) {
        auto it = potC[i].find(int(m));
        if (it != potC[i].end()) return it->second;
        const double v = wellPotential(wells_[i], m);
        potC[i][int(m)] = v;
        return v;
    };
    auto depth = [&](const std::string& g) { int d = 0; std::string x = g; while (x != "FIELD" && d < 32) { auto it = gParent_.find(x); x = it == gParent_.end() ? "FIELD" : it->second; ++d; } return d; };
    auto inGroup = [&](size_t i, const std::string& g) { for (auto& x : wellChain_[i]) if (x == g) return true; return false; };
    // C33 WEFAC: group quantities are efficiency-weighted (eff x flowing rate); shares are converted back to flowing targets
    auto prodQ = [&](const WellRT& w, Control m) {
        switch (m) {
            case Control::ORAT: return w.eff * w.rates.oil; case Control::WRAT: return w.eff * w.rates.water; case Control::GRAT: return w.eff * w.rates.gas;
            case Control::LRAT: return w.eff * (w.rates.oil + w.rates.water); case Control::RESV: return w.eff * std::max(0.0, w.rates.resv);
            default: return 0.0;
        }
    };
    auto field = [](WellControl& c, Control m) -> double& {
        switch (m) { case Control::ORAT: return c.orat; case Control::WRAT: return c.wrat; case Control::GRAT: return c.grat;
                     case Control::LRAT: return c.lrat; case Control::RESV: return c.resv; default: return c.rate; }
    };
    // C35 GSATPROD: satellite production below a group counts towards its rates as a fixed contribution
    auto satQ = [&](const std::string& g, Control m) {
        if (gsat_.empty()) return 0.0;
        switch (m) {
            case Control::ORAT: return satRate(g, 0); case Control::WRAT: return satRate(g, 1); case Control::GRAT: return satRate(g, 2);
            case Control::LRAT: return satRate(g, 0) + satRate(g, 1); case Control::RESV: return satRate(g, 3); default: return 0.0;
        }
    };
    // allocate a group target among wells: P = wells sharing it, the others keep their rates
    auto allocate = [&](const std::vector<size_t>& W, std::vector<size_t> P, Control m, double T, bool inj, const char* what, const std::string& gname) {
        double fixed = 0, Q = 0;
        for (size_t i : W) {
            const double q = inj ? wells_[i].eff * (wells_[i].type == WellType::GasInjector ? wells_[i].rates.gasInj : wells_[i].rates.waterInj) : prodQ(wells_[i], m);
            Q += q;
            if (std::find(P.begin(), P.end(), i) == P.end()) fixed += q;
        }
        if (!inj) { const double sq = satQ(gname, m); fixed += sq; Q += sq; }
        const double avail = std::max(0.0, T - fixed);
        std::vector<double> potW(P.size());
        double sum = 0;
        for (size_t j = 0; j < P.size(); ++j) {
            potW[j] = std::max(1e-9, wells_[P[j]].eff * pot(P[j], inj ? Control::RATE : m));
            // C34 GUIDERAT: producers share the target in proportion to their guide rates (converted to the controlled
            // quantity); the cap below stays the potential
            if (!inj && guide_.phase != "NONE" && guide_.phase != "") { const double gr = guideRateOf(wells_[P[j]], m); if (gr > 0) potW[j] = std::max(1e-12, wells_[P[j]].eff * gr); }
            sum += potW[j];
        }
        // C34 hierarchical guide rates: a subgroup with its own guide rate (GCONPROD items 9 / 10: a value, FORM =
        // GUIDERAT on the group potentials, POTN = potentials) takes the share GR_sub / sum(GR_siblings) of the parent
        // target, which its wells then share by their own guide rates; default subgroups pass the sum of their wells
        if (!inj && P.size() > 1) {
            bool anySub = false;
            for (auto& [gn, gc] : gprod_) if (gn != gname && (gc.guideRate > 0 || gc.guideDef == "FORM" || gc.guideDef == "POTN")) anySub = true;
            if (anySub) {
                std::map<size_t, double> wOf; for (size_t j = 0; j < P.size(); ++j) wOf[P[j]] = potW[j];
                std::vector<double> frac(P.size(), 0.0);
                // child of 'g' on the chain of well i (the group directly below g, or "" when the well itself is the child)
                auto childOf = [&](size_t i, const std::string& g) -> std::string {
                    const auto& ch = wellChain_[i];
                    for (size_t k = 0; k < ch.size(); ++k) if (ch[k] == g) return k == 0 ? std::string() : ch[k - 1];
                    return std::string();
                };
                std::function<void(const std::string&, const std::vector<size_t>&, double)> split = [&](const std::string& g, const std::vector<size_t>& idx, double f) {
                    std::map<std::string, std::vector<size_t>> kids; std::vector<size_t> direct;
                    for (size_t j : idx) { const std::string c = childOf(P[j], g); if (c.empty()) direct.push_back(j); else kids[c].push_back(j); }
                    std::vector<std::pair<double, std::vector<size_t>>> nodes; std::vector<std::string> names;
                    for (size_t j : direct) { nodes.push_back({potW[j], {j}}); names.push_back(""); }
                    for (auto& [c, js] : kids) {
                        double gr = 0; for (size_t j : js) gr += potW[j];
                        auto it = gprod_.find(c);
                        if (it != gprod_.end() && (it->second.guideDef == "FORM" || it->second.guideDef == "POTN" || it->second.guideRate > 0)) {
                            double po = 0, pwat = 0, pg = 0;
                            for (size_t j : js) { po += wells_[P[j]].eff * std::max(0.0, pot(P[j], Control::ORAT)); pwat += wells_[P[j]].eff * std::max(0.0, pot(P[j], Control::WRAT)); pg += wells_[P[j]].eff * std::max(0.0, pot(P[j], Control::GRAT)); }
                            const std::string ph = guide_.phase;
                            if (it->second.guideDef == "FORM" && ph != "NONE" && !ph.empty()) {
                                double pt = ph == "OIL" ? po : ph == "GAS" ? pg : ph == "WAT" ? pwat : ph == "RES" ? po + pwat + pg : po + pwat, r1, r2;
                                if (ph == "GAS") { r1 = pg > 0 ? pwat / pg : 0; r2 = pg > 0 ? po / pg : 0; }
                                else if (ph == "LIQ") { r1 = po + pwat > 0 ? pwat / (po + pwat) : 0; r2 = po + pwat > 0 ? pg / (po + pwat) : 0; }
                                else { r1 = po > 0 ? pwat / po : 0; r2 = po > 0 ? pg / po : 0; }
                                auto pw = [](double x, double e) { return e == 0 ? 1.0 : (x > 0 ? std::pow(x, e) : 0.0); };
                                const double den = guide_.B + guide_.C * pw(r1, guide_.D) + guide_.E * pw(r2, guide_.F);
                                gr = pt <= 0 || den <= 0 ? 0.0 : pw(pt, guide_.A) / den;
                                // converted to the controlled quantity with the group's production ratio
                                double qo = 0, qw = 0, qg = 0; for (size_t j : js) { qo += wells_[P[j]].rates.oil; qw += wells_[P[j]].rates.water; qg += wells_[P[j]].rates.gas; }
                                if (qo + qw + qg <= 0) { qo = po; qw = pwat; qg = pg; }
                                const double qn = ph == "OIL" ? qo : ph == "GAS" ? qg : ph == "WAT" ? qw : ph == "RES" ? qo + qw + qg : qo + qw;
                                const double qm = m == Control::ORAT ? qo : m == Control::WRAT ? qw : m == Control::GRAT ? qg : m == Control::LRAT ? qo + qw : qn;
                                if (qn > 0) gr *= qm / qn;
                            } else if (it->second.guideDef == "POTN") gr = m == Control::ORAT ? po : m == Control::WRAT ? pwat : m == Control::GRAT ? pg : po + pwat;
                            else if (it->second.guideRate > 0) gr = it->second.guideRate;
                        }
                        nodes.push_back({gr, js}); names.push_back(c);
                    }
                    double tot = 0; for (auto& nd : nodes) tot += nd.first;
                    if (tot <= 0) return;
                    for (size_t k = 0; k < nodes.size(); ++k) {
                        const double fk = f * nodes[k].first / tot;
                        if (names[k].empty()) frac[nodes[k].second.front()] = fk; else split(names[k], nodes[k].second, fk);
                    }
                };
                std::vector<size_t> all(P.size()); for (size_t j = 0; j < P.size(); ++j) all[j] = j;
                split(gname, all, 1.0);
                double fs = 0; for (double v : frac) fs += v;
                if (fs > 0) for (size_t j = 0; j < P.size(); ++j) potW[j] = std::max(1e-12, frac[j] * sum);
                (void)wOf;
            }
        }
        // water filling: shares in proportion to the weights (potentials / guide rates), capped by each well's own target and potential;
        // what a capped well cannot take goes to the others
        std::vector<double> cap(P.size()), shareV(P.size(), 0.0);
        std::vector<char> fixedW(P.size(), 0);
        for (size_t j = 0; j < P.size(); ++j) {
            const WellRT& w = wells_[P[j]];
            WellControl o = w.own;
            const double own = inj ? w.own.rate : field(o, m);
            cap[j] = std::max(1e-9, potW[j]);
            // an injector on its own RATE control with an explicit target of 0 is stopped: it takes no share of the
            // group target (ECLIPSE / OPM: group targets are shared only by wells under / available for group control)
            if (own > 0 || (inj && own == 0.0)) cap[j] = std::min(cap[j], std::max(0.0, own) * w.eff);
        }
        double left = avail;
        for (int pass = 0; pass < int(P.size()) + 1 && left > 1e-9; ++pass) {
            double sw = 0;
            for (size_t j = 0; j < P.size(); ++j) if (!fixedW[j]) sw += potW[j];
            if (sw <= 0) break;
            bool capped = false;
            for (size_t j = 0; j < P.size(); ++j) if (!fixedW[j] && left * potW[j] / sw >= cap[j]) { shareV[j] = cap[j]; fixedW[j] = 1; capped = true; }
            if (!capped) { for (size_t j = 0; j < P.size(); ++j) if (!fixedW[j]) shareV[j] = left * potW[j] / sw; break; }
            left = avail; for (size_t j = 0; j < P.size(); ++j) if (fixedW[j]) left -= shareV[j];
            left = std::max(0.0, left);
        }
        (void)sum;
        for (size_t j = 0; j < P.size(); ++j) {
            WellRT& w = wells_[P[j]];
            double share = shareV[j] / std::max(w.eff, 1e-9);   // flowing-rate target
            double& f = field(w.ctrl, inj ? Control::RATE : m);
            f = w.grpCtl && f > 0 ? std::min(f, share) : share;
            if (std::getenv("RESSIM_THPDBG")) std::fprintf(stderr, "GRP %s %s share %.1f avail %.1f fixed %.1f T %.1f\n", gname.c_str(), case_.wells[size_t(w.input)].name.c_str(), share, avail, fixed, T);
            if (w.mode == Control::BHP || !w.grpCtl) w.mode = inj ? Control::RATE : m;
            w.grpCtl = true; w.lockBhp = false;
        }
        (void)Q; (void)what; (void)gname;
    };
    // ---- production: top-down
    std::vector<std::pair<int, std::string>> order;
    for (auto& [n, g] : gprod_) order.push_back({depth(n), n});
    std::sort(order.begin(), order.end());
    for (auto& [d, name] : order) {
        const GroupProdControl& g = gprod_[name];
        std::vector<size_t> W, P;
        for (size_t i = 0; i < wells_.size(); ++i) if (wells_[i].open && wells_[i].type == WellType::Producer && inGroup(i, name)) W.push_back(i);
        if (W.empty()) continue;
        auto lim = [&](Control m) { switch (m) { case Control::ORAT: return g.orat; case Control::WRAT: return g.wrat; case Control::GRAT: return g.grat;
                                                 case Control::LRAT: return g.lrat; case Control::RESV: return g.resv; default: return -1.0; } };
        // rate of a well without group control: its potential capped by its own target (history wells: observed)
        auto unc = [&](size_t i, Control q) {
            const WellRT& w = wells_[i];
            if (w.own.history) return prodQ(w, q);
            double v = pot(i, q);
            WellControl o = w.own;
            const double own = field(o, q);
            if (own > 0) v = std::min(v, own);
            // the well's own controlling target also limits the other quantities (proportionally)
            if (w.own.mode != q && w.own.mode != Control::BHP) {
                const double t = field(o, w.own.mode), pm = pot(i, w.own.mode);
                if (t > 0 && pm > t) v *= t / pm;
            }
            return w.eff * v;
        };
        Control m = g.mode; double T = lim(m);
        if (T <= 0) m = Control::BHP;
        double worst = 1.0 + 1e-3;
        // C34: the group's own target competes with the other limits - the most violated constraint (largest
        // uncontrolled rate / limit) controls (e.g. LRAT 900 with an ORAT limit 1000 stays on LRAT)
        if (m != Control::BHP) { double Qm = satQ(name, m); for (size_t i : W) Qm += unc(i, m); worst = std::max(worst, Qm / T); }
        const Control qs[] = {Control::ORAT, Control::WRAT, Control::GRAT, Control::LRAT, Control::RESV};
        for (Control q : qs) {
            const double L = lim(q);
            if (L <= 0 || q == g.mode) continue;
            if (g.actionFor(q) == "NONE") continue;   // C35 (OPM / ECLIPSE): procedure NONE - the limit is only reported
            double Q = satQ(name, q); for (size_t i : W) Q += unc(i, q);
            if (Q / L > worst) { worst = Q / L; m = q; T = L; }   // a violated limit takes control (procedure RATE)
        }
        if (m == Control::BHP) continue;
        if (m == g.mode && chokeGroup(name)) continue;   // C35: the automatic network choke meets this target (common THP)
        double Q = satQ(name, m); for (size_t i : W) Q += unc(i, m);
        for (size_t i : W) if (!wells_[i].own.history && wells_[i].own.deckMode == "GRUP") P.push_back(i);
        const bool exceed = Q > T * (1 + 1e-3);
        if (P.empty() || (exceed && m != g.mode)) { if (!exceed) continue; P.clear(); for (size_t i : W) if (!wells_[i].own.history) P.push_back(i); }
        if (P.empty()) continue;
        allocate(W, P, m, T, false, "production", name);
    }
    // ---- injection
    for (int ph = 0; ph < 2; ++ph) {
        order.clear();
        for (auto& [n, g] : ginj_[ph]) order.push_back({depth(n), n});
        std::sort(order.begin(), order.end());
        const WellType wt = ph == 1 ? WellType::GasInjector : WellType::WaterInjector;
        for (auto& [d, name] : order) {
            const GroupInjControl& g = ginj_[ph][name];
            std::vector<size_t> W, P;
            for (size_t i = 0; i < wells_.size(); ++i) if (wells_[i].open && wells_[i].type == wt && inGroup(i, name)) W.push_back(i);
            if (W.empty()) continue;
            // formation volume factor of the injected phase at the mean pressure of the injection cells
            double pm = 0; int nc = 0;
            for (size_t i : W) for (size_t k = 0; k < wells_[i].cell.size(); ++k) if (wells_[i].wi[k] > 0) { pm += p_[wells_[i].cell[k]]; ++nc; }
            if (!nc) continue;
            pm /= nc;
            const FluidModel& F = fl(wells_[W.front()].cell.front());
            const double bInj = ph == 1 ? F.bgS(pm) : F.bwS(pm);
            const double bOther = ph == 1 ? F.bwS(pm) : F.bgS(pm);
            auto sumProd = [&](const std::string& grp, int what) {   // 0 gas, 1 water, 2 voidage
                double q = 0;
                for (size_t i = 0; i < wells_.size(); ++i) {
                    if (!wells_[i].open || wells_[i].type != WellType::Producer || !inGroup(i, grp)) continue;
                    q += wells_[i].eff * (what == 0 ? wells_[i].rates.gas : what == 1 ? wells_[i].rates.water : std::max(0.0, wells_[i].rates.resv));
                }
                return q;
            };
            double T = -1;
            if (g.mode == "RATE") T = g.rate;
            else if (g.mode == "RESV" && g.resv > 0) T = g.resv * bInj;
            else if (g.mode == "REIN" && g.rein >= 0) {
                const std::string rg = g.reinGroup.empty() ? name : g.reinGroup;
                // C35: gas re-injection works on the gas available after consumption / import (GCONSUMP) and incl. satellites
                T = ph == 1 ? g.rein * std::max(0.0, groupGasAvail(rg)) : g.rein * sumProd(rg, 1);
            }
            else if (g.mode == "VREP" && g.vrep >= 0) {
                double other = 0;   // reservoir volume injected by the other phase in this group
                for (size_t i = 0; i < wells_.size(); ++i) {
                    const WellRT& w = wells_[i];
                    if (!w.open || w.type == WellType::Producer || w.type == wt || !inGroup(i, name)) continue;
                    other += w.eff * (w.type == WellType::GasInjector ? w.rates.gasInj : w.rates.waterInj) / std::max(1e-12, bOther);
                }
                T = std::max(0.0, g.vrep * sumProd(g.vrepGroup.empty() ? name : g.vrepGroup, 2) - other) * bInj;
            }
            if (g.mode != "RATE" && g.rate > 0) T = T < 0 ? g.rate : std::min(T, g.rate);
            if (g.mode != "RESV" && g.resv > 0) T = T < 0 ? g.resv * bInj : std::min(T, g.resv * bInj);
            if (ph == 1) {   // C35 GCONSALE: the group's gas injection is limited so that its sales target is met
                std::string un = name; for (char& c : un) c = char(std::toupper(static_cast<unsigned char>(c)));
                auto it = gsale_.find(un);
                if (it != gsale_.end() && it->second.target >= 0) {
                    const double ts = std::max(0.0, groupGasAvail(name) - it->second.target);
                    T = T < 0 ? ts : std::min(T, ts);
                }
            }
            if (T < 0) continue;
            double Q = 0;
            for (size_t i : W) {
                const WellRT& w = wells_[i];
                if (w.own.history) { Q += w.eff * (ph == 1 ? w.rates.gasInj : w.rates.waterInj); continue; }
                double v = pot(i, Control::RATE);
                if (w.own.rate >= 0) v = std::min(v, w.own.rate);   // an explicit target of 0 = stopped (rate < 0: no target)
                Q += w.eff * v;
            }
            for (size_t i : W) if (!wells_[i].own.history && wells_[i].own.deckMode == "GRUP") P.push_back(i);
            const bool exceed = Q > T * (1 + 1e-3);
            if (P.empty()) { if (!exceed) continue; for (size_t i : W) if (!wells_[i].own.history) P.push_back(i); }
            if (P.empty()) continue;
            allocate(W, P, Control::RATE, T, true, "injection", name);
        }
    }
}

// ALQ value of a VFPPROD lookup. A well on gas-lift optimisation (C35, liftRate) carries its lift-gas RATE, converted to the
// table's ALQ quantity: GRAT / blank = lift gas rate, IGLR = lift gas / liquid, TGLR = (produced + lift gas) / liquid.
// Otherwise the WCONPROD / WELTARG / WLIFTOPT value is the table quantity itself.
static double alqLookup(const VfpTable& t, double alq, bool liftRate, double qo, double qw, double qg) {
    if (!liftRate) return std::max(0.0, alq);
    std::string ty = t.alq;
    ty.erase(std::remove(ty.begin(), ty.end(), '\''), ty.end());
    const double ql = qo + qw;
    if (ty == "IGLR") return ql > 0 ? alq / ql : 0.0;
    if (ty == "TGLR") return ql > 0 ? (qg + alq) / ql : 0.0;
    return std::max(0.0, alq);
}

const VfpTable* Simulator::wellVfp(const WellRT& w) const {
    if (w.own.vfp <= 0) return nullptr;
    const bool prod = w.type == WellType::Producer;
    const VfpTable* t = nullptr;
    for (const VfpTable& v : case_.vfp) if (v.id == w.own.vfp && v.inj == !prod) t = &v;
    return t;
}

// BHP at the THP limit = operating point of the well's inflow (IPR at the current reservoir state) and the lift
// curve: bhp = VFP(q(bhp)), found by bisection (OPM Flow computeBhpAtThpLimit practice). Stable where the lagged
// fixed point (bhp from the last rates) oscillates or locks a well at an unreachable BHP. q (optional): producer surface
// rates [W, O, G] at max(bhp, BHP limit) - used by the gas-lift optimisation (C35) for the lift performance per ALQ.
bool Simulator::thpOperatingPoint(const WellRT& w, double alq, bool liftRate, double& bhpOut, double* qOut) const {
    const double g = U_.gravity;
    const double thpE = std::max(w.own.thp, w.thpNet);   // C29: deck THP limit or the network node pressure
    if (!w.open || !(thpE > 0) || w.own.vfp <= 0) return false;
    const bool prod = w.type == WellType::Producer;
    const VfpTable* t = wellVfp(w);
    if (!t) return false;
    // inflow: per connection and phase q = a * (P - dp - bhp)+ (surface rates), with Rs / Rv mixing
    struct Conn { double a[3], P[3], rs, rv, lt, bInj, wi; };
    std::vector<Conn> cs;
    const int ip = w.type == WellType::GasInjector ? PG : PW;
    double pMax = U_.atm;
    for (size_t k = 0; k < w.cell.size(); ++k) {
        if (w.wi[k] <= 0) continue;
        const CellFull C = evalCell(w.cell[k]);
        const double dp = k < w.dpConn.size() ? w.dpConn[k] : w.headRho * g * (w.cdepth[k] - w.refDepth);
        Conn c{};
        c.lt = 0; c.wi = w.wi[k];
        for (int p = 0; p < 3; ++p) {
            c.a[p] = eq_[p] >= 0 ? w.wi[k] * C.mob[p].v : 0.0;
            c.P[p] = eq_[p] >= 0 ? C.P[p].v - dp : 0.0;
            if (eq_[p] >= 0 && C.b[p].v > 0) c.lt += C.mob[p].v / C.b[p].v;
            if (eq_[p] >= 0) pMax = std::max(pMax, C.P[p].v - dp);
        }
        c.rs = C.rs.v; c.rv = C.rv.v; c.bInj = eq_[ip] >= 0 ? C.b[ip].v : 1.0;
        cs.push_back(c);
    }
    if (cs.empty()) return false;
    auto rates = [&](double bh, double& o, double& wv, double& gv) {
        o = wv = gv = 0;
        for (const Conn& c : cs) {
            double q[3];
            for (int p = 0; p < 3; ++p) q[p] = c.a[p] * std::max(0.0, c.P[p] - bh);
            wv += q[PW]; o += q[PO] + c.rv * q[PG]; gv += q[PG] + c.rs * q[PO];
        }
    };
    auto vfpAt = [&](double bhp) {
        if (prod) {
            double qo = 0, qw = 0, qg = 0, ro = 0, rw = 0, rg = 0;
            rates(bhp, qo, qw, qg);
            ro = qo; rw = qw; rg = qg;
            // C30: at (or above) the shut-in pressure the ratios are their limits at a vanishing drawdown, not 0 (a GOR of 0
            // made the lift curve heaviest exactly at q -> 0, so a gassy well was declared unable to flow)
            if (!(qo + qw + qg > 1e-12)) rates(bhp - 1e-3 * std::max(pMax, U_.atm), ro, rw, rg);
            const double flo = t->flo == "LIQ" ? qo + qw : t->flo == "GAS" ? qg : qo;
            const double wfr = t->wfr == "WOR" ? (ro > 0 ? rw / ro : 0) : t->wfr == "WGR" ? (rg > 0 ? rw / rg : 0) : (ro + rw > 0 ? rw / (ro + rw) : 0);
            const double gfr = t->gfr == "GLR" ? (ro + rw > 0 ? rg / (ro + rw) : 0) : t->gfr == "OGR" ? (rg > 0 ? ro / rg : 0) : (ro > 0 ? rg / ro : 0);
            return t->eval(flo, thpE, wfr, gfr, alqLookup(*t, alq, liftRate, ro, rw, rg));
        }
        double q = 0;
        for (const Conn& c : cs) q += c.wi * c.lt * c.bInj * std::max(0.0, bhp - c.P[ip]);   // total mobility injection
        return t->eval(q, thpE, 0, 0, 0);
    };
    const double datumCorr = w.headRho * g * (w.refDepth - t->datum);   // table datum -> well reference depth
    auto f = [&](double bhp) { return bhp - (vfpAt(bhp) + datumCorr); };
    double lo = U_.atm, hi = prod ? pMax : std::max(pMax * 3.0, w.bhpDeck > 0 ? std::min(w.bhpDeck, 1e6) : pMax * 3.0);
    double b;
    if (prod) {
        // C35: a J-shaped lift curve (liquid loading at low rates) gives two operating points; the stable one is the
        // lowest BHP / highest rate root, where f crosses from - to + with increasing BHP (OPM computeBhpAtThpLimitProd
        // takes the stable solution). Bisection over [atm, p_res] could land on either root, and the BHP limit then
        // flipped between them with the lagged state (well switching ORAT <-> BHP every iteration, rate limit exceeded).
        const double flo = f(lo), fhi = f(hi);
        if (flo >= 0) b = lo;
        else {
            const int ns = 48;
            double a0 = lo, fa = flo, b0 = -1;
            double xBest = lo, fBest = flo;
            for (int s2 = 1; s2 <= ns; ++s2) {
                const double x = lo + (hi - lo) * double(s2) / ns, fx = s2 == ns ? fhi : f(x);
                if (fa < 0 && fx >= 0) { b0 = x; break; }
                if (fx > fBest) { fBest = fx; xBest = x; }
                a0 = x; fa = fx;
            }
            if (b0 < 0) {
                // C35: the window where the well can flow (f > 0 between the stable and the unstable root) may be narrower
                // than the sample spacing (low THP, rates close to the shut-in pressure): maximise f around the best sample
                // (golden section) before declaring that the well cannot lift. Missing the window made the limit flip
                // between the stable root and 'no flow' with tiny state changes (ORAT <-> BHP switching).
                const double dx = (hi - lo) / ns;
                double gl = std::max(lo, xBest - dx), gh = std::min(hi, xBest + dx);
                const double r = 0.6180339887498949;
                double x1 = gh - r * (gh - gl), x2 = gl + r * (gh - gl), f1 = f(x1), f2 = f(x2);
                for (int it = 0; it < 40 && std::max(f1, f2) < 0; ++it) {
                    if (f1 > f2) { gh = x2; x2 = x1; f2 = f1; x1 = gh - r * (gh - gl); f1 = f(x1); }
                    else { gl = x1; x1 = x2; f1 = f2; x2 = gl + r * (gh - gl); f2 = f(x2); }
                }
                const double xm = f1 >= f2 ? x1 : x2;
                if (std::max(f1, f2) >= 0 && f(std::max(lo, xBest - dx)) < 0) { a0 = std::max(lo, xBest - dx); b0 = xm; }
            }
            if (b0 < 0) b = hi;               // the well cannot lift any fluid: limit at the reservoir pressure (no flow)
            else { double l2 = a0, h2 = b0; for (int it = 0; it < 50; ++it) { const double m = 0.5 * (l2 + h2); (f(m) > 0 ? h2 : l2) = m; } b = 0.5 * (l2 + h2); }
        }
    } else {
        lo = std::max(U_.atm, pMax * 0.5);
        if (f(lo) >= 0) b = lo;
        else if (f(hi) <= 0) b = hi;
        else { for (int it = 0; it < 60; ++it) { const double m = 0.5 * (lo + hi); (f(m) > 0 ? hi : lo) = m; } b = 0.5 * (lo + hi); }
    }
    bhpOut = b;
    if (qOut) {
        double qo = 0, qw = 0, qg = 0;
        if (prod) rates(std::max(b, w.bhpDeck > 0 ? w.bhpDeck : U_.atm), qo, qw, qg);
        qOut[PW] = qw; qOut[PO] = qo; qOut[PG] = qg;
    }
    return true;
}

void Simulator::updateThp() {
    for (auto& w : wells_) {
        const double thpE = std::max(w.own.thp, w.thpNet);
        if (!w.open || !(thpE > 0) || w.own.vfp <= 0) continue;
        const bool prod = w.type == WellType::Producer;
        const std::string& nm = case_.wells[size_t(w.input)].name;
        if (!wellVfp(w)) {
            if (vfpMissing_.insert(nm).second) log_("Warning: well " + nm + fmt(": VFP table %g not found - THP limit ignored", double(w.own.vfp)));
            continue;
        }
        double b = 0;
        if (!thpOperatingPoint(w, prod ? w.alq : 0.0, w.alqLift, b, nullptr)) continue;
        const double deck = w.bhpDeck > 0 ? w.bhpDeck : (prod ? U_.atm : 1e12);
        w.ctrl.bhp = prod ? std::max(deck, b) : std::min(deck, b);
        if (std::getenv("RESSIM_THPDBG")) std::fprintf(stderr, "THP %s q %.1f/%.1f/%.1f/%.1f -> bhp(thp) %.1f lim %.1f bhp %.1f mode %d alq %.4g\n", nm.c_str(),
            w.rates.oil, w.rates.water, w.rates.gas, w.rates.waterInj, b, w.ctrl.bhp, w.bhp, int(w.mode), w.alq);
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
        } else if (!w.lockBhp) {   // (a well locked by the in-Newton rule retries its rate target at the next report step)
            Control best = Control::BHP; double worst = 1.0 + 1e-6;
            const Control cand[] = {Control::ORAT, Control::WRAT, Control::GRAT, Control::LRAT, Control::RESV, Control::RATE};
            for (Control m : cand) {
                if (prod == (m == Control::RATE)) continue;
                // WCONHIST / WCONINJH: the observed rates are not limits - only the history control mode returns
                if (w.ctrl.history && m != w.own.mode) continue;
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
        // consistency: SWCR < SWU, SGCR < SGU, critical oil below the mobile range. C44: SWCR / SGCR are NOT raised to SWL / SGL
        // - OPM (EclEpsScalingPointsInfo::extractScaled) keeps a cell SWL above the table SWCR, so water is mobile at connate
        // saturation there (opm-tests model2 9_* / 6_ENDSCALE: SWL array up to 0.19, table SWCR 0.10)
        const SatFunctions::EndPoints o = e;
        e.swu = std::max(e.swu, e.swcr + 1e-6);
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

// ======================================================================== hysteresis / VAPPARS
void Simulator::setupHysteresis(const Case& c) {
    hyst_ = c.props.hyst.on && O_ && (W_ || G_);
    hystWetImb_ = hyst_ && (c.props.hyst.model == 1 || c.props.hyst.model == 3);
    imbReg_.clear(); ieps_.clear(); sgHy_.clear(); soHy_.clear(); dG_.clear(); dO_.clear();
    if (!hyst_) return;
    const int na = grid_.numActive(), nsat = int(sats_.size());
    imbReg_ = satReg_;
    if (!c.grid.imbnum.empty())
        for (int a = 0; a < na; ++a) {
            const size_t q = size_t(grid_.globalOf[size_t(a)]);
            const int r = q < c.grid.imbnum.size() ? c.grid.imbnum[q] - 1 : satReg_[size_t(a)];
            imbReg_[size_t(a)] = (r >= 0 && r < nsat) ? r : satReg_[size_t(a)];
        }
    const GridInput& g = c.grid;
    const bool iArrays = !(g.iswl.empty() && g.iswcr.empty() && g.iswu.empty() && g.isgl.empty() && g.isgcr.empty() && g.isgu.empty() &&
                           g.isowcr.empty() && g.isogcr.empty());
    if (!eps_.empty() || iArrays) {
        ieps_.resize(size_t(na));
        for (int a = 0; a < na; ++a) {
            const size_t q = size_t(grid_.globalOf[size_t(a)]);
            SatFunctions::EndPoints e = sti(a).tableEndPoints();
            auto get = [&](const std::vector<double>& arr, double& v) { if (q < arr.size() && std::isfinite(arr[q]) && arr[q] >= 0) v = std::min(1.0, arr[q]); };
            get(g.iswl, e.swl); get(g.iswcr, e.swcr); get(g.iswu, e.swu); get(g.isgl, e.sgl);
            get(g.isgcr, e.sgcr); get(g.isgu, e.sgu); get(g.isowcr, e.sowcr); get(g.isogcr, e.sogcr);
            e.swu = std::max(e.swu, e.swcr + 1e-6);   // C44: no SWCR >= SWL / SGCR >= SGL clamp (as OPM)
            e.sgu = std::max(e.sgu, e.sgcr + 1e-6);
            e.sowcr = std::min(e.sowcr, std::max(0.0, 1.0 - e.swl - 1e-6));
            e.sogcr = std::min(e.sogcr, std::max(0.0, 1.0 - e.swl - e.sgl - 1e-6));
            ieps_[size_t(a)] = e;
        }
    }
    int nImb = 0;
    for (int a = 0; a < na; ++a) if (imbReg_[size_t(a)] != satReg_[size_t(a)]) ++nImb;
    const int m = c.props.hyst.model;
    log_(std::string("Relative permeability hysteresis: Carlson (non-wetting gas and oil), ") +
         (m == 1 || m == 3 ? "wetting phases (krw, krog) on the imbibition curve" : "wetting phases on drainage") +
         (m >= 2 ? " - EHYSTR Killough model approximated by Carlson" : "") +
         fmt(" - %g cells with a separate imbibition table", double(nImb)) + (iArrays ? ", imbibition end points" : ""));
}

double Simulator::soOW(int a, double sw, double sg) const {
    return st(a).kroModel() != 0 ? 1.0 - sw : 1.0 - sw - sg;
}

void Simulator::updateHysteresis() {
    if (!hyst_) return;
    const int na = grid_.numActive();
    const bool init = sgHy_.empty();
    if (init) { sgHy_.assign(size_t(na), 0.0); soHy_.assign(size_t(na), 0.0); dG_.assign(size_t(na), NAN); dO_.assign(size_t(na), NAN); }
    // inverse of a monotone (non-decreasing) kr curve on [lo, hi]
    auto inv = [](auto&& kr, double target, double lo, double hi) {
        if (kr(hi) <= target) return hi;
        for (int it = 0; it < 50; ++it) { const double m = 0.5 * (lo + hi); (kr(m) < target ? lo : hi) = m; }
        return 0.5 * (lo + hi);
    };
    auto body = [&](int a) {
        const double sw = iSw_ >= 0 ? sw_[a] : 0.0;
        const int stt = G_ ? hasGas_[a] : 1;
        const double sg = !G_ ? 0.0 : stt == 1 ? x_[a] : stt == 2 ? 1.0 - sw : 0.0;
        const SatFunctions& D = st(a); const SatFunctions& I = sti(a);
        const SatFunctions::EndPoints* E = ep(a); const SatFunctions::EndPoints* IE = iep(a);
        const size_t u = size_t(a);
        if (G_) {
            if (init || sg > sgHy_[u]) {
                sgHy_[u] = sg;
                const double kd = D.krg(AD3(sg), E).v;
                if (kd <= 1e-12) dG_[u] = NAN;
                else {
                    const double lo = IE ? IE->sgcr : I.tableEndPoints().sgcr, hi = IE ? IE->sgu : I.tableEndPoints().sgu;
                    const double x = inv([&](double s) { return I.krg(AD3(s), IE).v; }, kd, lo, hi);
                    dG_[u] = x - sg;
                }
            }
        }
        if (W_) {
            const double so = soOW(a, sw, sg);
            if (init || so > soHy_[u]) {
                soHy_[u] = so;
                const double kd = D.krowAtSo(AD3(so), E).v;
                if (kd <= 1e-12) dO_[u] = NAN;
                else {
                    const SatFunctions::EndPoints& T = IE ? *IE : I.tableEndPoints();
                    const double x = inv([&](double s) { return I.krowAtSo(AD3(s), IE).v; }, kd, T.sowcr, 1.0 - T.swl);
                    dO_[u] = x - so;
                }
            }
        }
    };
    if (useParallel(na > 4000)) {
#pragma omp parallel for schedule(static)
        for (int a = 0; a < na; ++a) body(a);
    } else for (int a = 0; a < na; ++a) body(a);
}

void Simulator::updateVappars() {
    if (!(O_ && G_) || (vap1_ <= 0 && vap2_ <= 0)) { fRv_.clear(); fRs_.clear(); return; }
    const int na = grid_.numActive();
    if (soMaxV_.size() != size_t(na)) soMaxV_.assign(size_t(na), 0.0);
    if (vap1_ > 0 && vap_) fRv_.assign(size_t(na), 1.0); else fRv_.clear();
    if (vap2_ > 0 && live_) fRs_.assign(size_t(na), 1.0); else fRs_.clear();
    for (int a = 0; a < na; ++a) {
        const double sw = iSw_ >= 0 ? sw_[a] : 0.0;
        const int stt = hasGas_[a];
        const double sg = stt == 1 ? x_[a] : stt == 2 ? 1.0 - sw : 0.0;
        const double so = stt == 2 ? 0.0 : std::max(0.0, 1.0 - sw - sg);
        double& smax = soMaxV_[size_t(a)];
        smax = std::min(1.0, std::max(smax, so));
        // OPM WetGasPvt / LiveOilPvt: factor (So / SoMax)^vap, floored at 1e-3, only below the historical maximum
        if (smax > 0.01 && so < smax) {
            const double r = std::max(so, 1e-3) / smax;
            if (!fRv_.empty()) fRv_[size_t(a)] = std::max(1e-3, std::pow(r, vap1_));
            if (!fRs_.empty()) fRs_[size_t(a)] = std::max(1e-3, std::pow(r, vap2_));
        }
    }
}

// ======================================================================== cell properties
Simulator::CellFull Simulator::evalCell(int c, int satOverride) const {
    CellFull f;
    const FluidModel& fluid_ = fl(c);
    const SatFunctions& sat_ = satOverride >= 0 ? sats_[size_t(satOverride)] : st(c);
    const SatFunctions::EndPoints* E = ep(c);
    const DensitySC& rhoSC = fluid_.rhoSC();
    const AD3 P = AD3::variable(p_[c], 0);
    AD3 Sw(W_ && !O_ ? 1.0 : 0.0), Sg(0.0), rs(0.0), rv(0.0);
    const double fRv = fRv_.empty() ? 1.0 : fRv_[size_t(c)], fRs = fRs_.empty() ? 1.0 : fRs_[size_t(c)];
    if (iSw_ >= 0) Sw = AD3::variable(sw_[c], iSw_);
    if (G_) {
        const int state = hasGas_[c];
        if (state == 1) {
            Sg = AD3::variable(x_[c], iX_);
            if (live_) { rs = fluid_.rsSat(P) * fRs; if (!rsLim_.empty() && rs.v > rsLim_[size_t(c)]) rs = AD3(rsLim_[size_t(c)]); }
            if (vap_) rv = fluid_.rvSat(P) * fRv;
        }
        else if (state == 0) { rs = AD3::variable(x_[c], iX_); if (vap_) rv = fluid_.rvSat(P) * fRv; }
        else { Sg = 1.0 - Sw; rv = AD3::variable(x_[c], iX_); if (live_) rs = fluid_.rsSat(P) * fRs; }
    }
    const AD3 So = (O_ && !(G_ && hasGas_[c] == 2)) ? (1.0 - Sw - Sg) : AD3(0.0);
    // C35 SOLVENT: with a solvent the gas-phase variable is the total Sg + Ss (saturation functions see the total, OPM);
    // the split follows from the lagged solvent fraction z of the G+S component (flash, solvent.cpp header)
    AD3 SgF = Sg, Ss(0.0), bsv(0.0), pgS(0.0);
    const double zS = solv_ ? zs_[size_t(c)] : 0.0;
    if (solv_ && G_) {
        pgS = P + sat_.pcog(Sg, E);
        bsv = svB_[std::min(size_t(pvtReg_[size_t(c)]), svB_.size() - 1)](pgS);
        if (zS > 0 && (hasGas_[c] != 0 || !O_)) {
            const AD3 bgv = vap_ ? fluid_.bg(pgS, rv) : fluid_.bg(pgS);
            if (O_ && hasGas_[c] == 1) {
                Ss = zS * (Sg * bgv + So * fluid_.bo(P, rs) * rs) / (bsv * (1.0 - zS) + zS * bgv);
                if (Ss.v > Sg.v) {                       // the gas phase is pure solvent: Rs from the composition
                    Ss = Sg;
                    if (live_ && So.v > 1e-12) {
                        AD3 r = rs;
                        for (int it = 0; it < 8; ++it) {
                            AD3 rn = Ss * bsv * (1.0 - zS) / (zS * So * fluid_.bo(P, r));
                            if (rn.v > rs.v) rn = rs;
                            const double d = std::fabs(rn.v - r.v);
                            r = rn;
                            if (d < 1e-12 * std::max(1.0, rs.v)) break;
                        }
                        rs = r;
                    }
                }
            } else Ss = zS * Sg * bgv / (bsv * (1.0 - zS) + zS * bgv);
            SgF = Sg - Ss;
            if (SgF.v < 0) SgF = AD3(0.0);
        }
    }
    AD3 pv = fluid_.pvMult(P) * grid_.poreVolRef[c];
    if (rockOn_) pv *= rockMult(c, P, Sw, false);   // C34 ROCKTAB / ROCK2D pore volume multiplier
    const double g = U_.gravity;
    (void)g;
    // C35 surfactant: miscibility F from the capillary number - water-oil end points interpolated between the immiscible
    // table and the miscible SURFNUM table, both curves scaled to them and blended with F (lagged F, ECLIPSE SURFACT)
    const double Fm = chMisc_.empty() ? 0.0 : chMisc_[size_t(c)];
    const SatFunctions* Sm = nullptr;
    SatFunctions::EndPoints Eb{};
    if (Fm > 0 && !chMisSat_.empty() && chMisSat_[size_t(c)] >= 0) {
        Sm = &sats_[size_t(chMisSat_[size_t(c)])];
        const SatFunctions::EndPoints& Ei = E ? *E : sat_.tableEndPoints();
        const SatFunctions::EndPoints& Em = Sm->tableEndPoints();
        Eb = Ei;
        Eb.swcr = (1 - Fm) * Ei.swcr + Fm * Em.swcr; Eb.sowcr = (1 - Fm) * Ei.sowcr + Fm * Em.sowcr; Eb.sogcr = (1 - Fm) * Ei.sogcr + Fm * Em.sogcr;
    }
    if (W_) {
        f.P[PW] = (O_ ? P - sat_.pcow(Sw, E) * (pcScale_.empty() ? 1.0 : pcScale_[size_t(c)]) * (chPcMul_.empty() ? 1.0 : chPcMul_[size_t(c)]) : P);
        f.b[PW] = fluid_.bw(f.P[PW]);
        const AD3 kr = !O_ ? AD3(1.0) : Sm ? (1 - Fm) * sat_.krw(Sw, &Eb) + Fm * Sm->krw(Sw, &Eb)
                     : hystWetImb_ ? sti(c).krw(Sw, iep(c)) : sat_.krw(Sw, E);   // C44 EHYSTR 1/3: wetting water on imbibition
        f.mob[PW] = fluid_.bwMu(f.P[PW]) * kr;
        if (poly_ && !wmul_.empty()) f.mob[PW] *= wmul_[size_t(c)];   // polymer viscosity correction / RRF
        if (!chWmul_.empty()) f.mob[PW] *= chWmul_[size_t(c)];         // C35 surfactant water viscosity
        f.rho[PW] = f.b[PW] * rhoSC.water;
        f.acc[PW] = pv * f.b[PW] * Sw;
        f.S[PW] = Sw;
    }
    if (O_) {
        f.P[PO] = P;
        f.b[PO] = fluid_.bo(P, rs);
        AD3 kr(1.0);
        const bool carlsonO = hyst_ && W_ && !dO_.empty() && std::isfinite(dO_[size_t(c)]);
        if (!Sm && (carlsonO || (hystWetImb_ && G_))) {
            // Carlson scanning curve for oil (non-wetting in oil-water): imbibition curve shifted to meet drainage at the turning
            // point. C44 EHYSTR 1/3: oil is the wetting phase of the gas-oil system -> krog on the imbibition curve (OPM
            // EclHysteresisTwoPhaseLaw::twoPhaseSatKrw)
            const bool stone = sat_.kroModel() != 0;
            const double swl = E ? E->swl : sat_.swco();
            const AD3 soW = stone ? 1.0 - Sw : So;
            const AD3 krow = !W_ ? AD3(0.0) : carlsonO && soW.v < soHy_[size_t(c)] - 1e-9 ? sti(c).krowAtSo(soW + dO_[size_t(c)], iep(c)) : sat_.krowAtSo(soW, E);
            const AD3 krog = !G_ ? AD3(0.0) : hystWetImb_ ? sti(c).krogAtSo(stone ? 1.0 - swl - Sg : So, iep(c)) : sat_.krogAtSo(stone ? 1.0 - swl - Sg : So, E);
            kr = sat_.kroCombine(Sw, Sg, krow, krog, E);
        } else if (Sm) kr = (1 - Fm) * sat_.kro(Sw, Sg, &Eb) + Fm * Sm->kro(Sw, Sg, &Eb);   // C35 surfactant
        else if (W_ || G_) kr = sat_.kro(Sw, Sg, E);
        f.mob[PO] = fluid_.boMu(P, rs) * kr;
        f.rho[PO] = f.b[PO] * (rs * (rhoSC.gas * fluid_.gasFactor()) + rhoSC.oil);
        f.acc[PO] = pv * f.b[PO] * So;
        f.mobRs = f.mob[PO] * rs;
        f.S[PO] = So;
    }
    if (G_) {
        f.P[PG] = P + sat_.pcog(Sg, E);
        AD3 krgC;
        if (hyst_ && O_ && !dG_.empty() && std::isfinite(dG_[size_t(c)]) && Sg.v < sgHy_[size_t(c)] - 1e-9)
            krgC = sti(c).krg(Sg + dG_[size_t(c)], iep(c));     // Carlson: trapped gas on imbibition
        else krgC = sat_.krg(Sg, E);
        if (vap_) {
            f.b[PG] = fluid_.bg(f.P[PG], rv);
            f.mob[PG] = fluid_.bgMu(f.P[PG], rv) * krgC;
            f.rho[PG] = f.b[PG] * (rhoSC.gas * fluid_.gasFactor() + rv * rhoSC.oil);
            if (O_) f.acc[PO] += pv * rv * f.b[PG] * SgF;
        } else {
            f.b[PG] = fluid_.bg(f.P[PG]);
            f.mob[PG] = fluid_.bgMu(f.P[PG]) * krgC;
            f.rho[PG] = f.b[PG] * (rhoSC.gas * fluid_.gasFactor());
        }
        bool foamOnSolvent = false;
        if (solv_) {   // C35 SOLVENT: immiscible split of the total gas kr by the SSFN multipliers (OPM blackoilsolventmodules)
            const size_t sr = std::min(size_t(satReg_[size_t(c)]), ssfnG_.size() - 1);
            const AD3 Fs = Sg.v > 1e-12 ? Ss / Sg : AD3(0.0);
            f.mob[PG] *= ssfnG_[sr](1.0 - Fs);
            f.mobS = svBMu_[std::min(size_t(pvtReg_[size_t(c)]), svBMu_.size() - 1)](pgS) * krgC * ssfnS_[sr](Fs);
            f.rhoS = bsv * (svRho_[std::min(size_t(pvtReg_[size_t(c)]), svRho_.size() - 1)] * fluid_.gasFactor());
            foamOnSolvent = chFoam_ >= 0 && chem_[size_t(chFoam_)].phase == 3;
            if (foamOnSolvent && !chGmul_.empty()) f.mobS *= chGmul_[size_t(c)];   // FOAMOPTS SOLVENT
        }
        if (vap_) f.mobRv = f.mob[PG] * rv;
        if (!chGmul_.empty() && !foamOnSolvent) f.mob[PG] *= chGmul_[size_t(c)];   // C35 foam gas mobility reduction
        f.acc[PG] = pv * (f.b[PG] * SgF + (O_ ? rs * f.b[PO] * So : AD3(0.0)));
        if (solv_) f.acc[PG] += pv * bsv * Ss;
        if (!langK_.empty() && langK_[size_t(c)] > 0) f.acc[PG] += P * (langK_[size_t(c)] / langPL_) / (P / langPL_ + 1.0);   // Langmuir
        f.S[PG] = SgF;
    }
    if (rockTm_) {   // C34 ROCKTAB / ROCK2DTR transmissibility multiplier (upstream cell; wells use their cell's value as OPM)
        const AD3 tm = rockMult(c, P, Sw, true);
        for (int p = 0; p < 3; ++p) f.mob[p] *= tm;
        f.mobRs *= tm; f.mobRv *= tm; f.mobS *= tm;
    }
    f.rs = rs; f.rv = rv; f.ss = Ss; f.bs = bsv;
    return f;
}

// ======================================================================== C34 rock compaction
void Simulator::setupRockComp() {
    const RockCompInput& R = case_.props.rockComp;
    rockOn_ = rockTm_ = false;
    rockReg_.clear(); rockPob_.clear(); rockPmin_.clear(); rockSwMax_.clear(); rockSw0_.clear(); rockP0_.clear();
    if (!R.active()) return;
    const int na = grid_.numActive();
    rockReg_.assign(size_t(na), 0);
    for (int a = 0; a < na; ++a) {
        const size_t q = size_t(grid_.globalOf[size_t(a)]);
        int r = 0;
        if (R.regionArray == "ROCKNUM") r = q < case_.grid.rocknum.size() ? case_.grid.rocknum[q] - 1 : 0;
        else if (R.regionArray == "SATNUM") r = satReg_[size_t(a)];
        else r = pvtReg_[size_t(a)];
        const int nt = int(R.tab.empty() ? R.tab2.size() : R.tab.size());
        rockReg_[size_t(a)] = std::clamp(r, 0, std::max(0, nt - 1));
    }
    auto interp = [](const std::vector<double>& X, const std::vector<double>& Y, double x) {
        if (X.empty() || Y.size() != X.size()) return 0.0;
        if (X.size() == 1) return Y[0];
        size_t i = size_t(std::upper_bound(X.begin(), X.end(), x) - X.begin()); i = std::clamp<size_t>(i, 1, X.size() - 1);
        return Y[i - 1] + (Y[i] - Y[i - 1]) * (x - X[i - 1]) / std::max(X[i] - X[i - 1], 1e-30);
    };
    rockPob_.assign(size_t(na), 0.0);
    if (!R.overburd.empty() && !R.store)
        for (int a = 0; a < na; ++a) {
            const auto& v = R.overburd[std::min(size_t(rockReg_[size_t(a)]), R.overburd.size() - 1)];
            std::vector<double> d, p; for (size_t m = 0; m + 1 < v.size(); m += 2) { d.push_back(v[m]); p.push_back(v[m + 1]); }
            rockPob_[size_t(a)] = interp(d, p, grid_.depth[size_t(a)]);
        }
    if (R.irreversible) rockPmin_.assign(size_t(na), 1e99);
    if (!R.tab2.empty()) { rockSwMax_.assign(size_t(na), 0.0); rockSw0_.assign(size_t(na), 0.0); for (int a = 0; a < na; ++a) rockSw0_[size_t(a)] = iSw_ >= 0 || W_ ? sw_[size_t(a)] : 0.0; }
    if (R.store) rockP0_ = p_;
    rockOn_ = true;
    if (!R.tab.empty()) for (const auto& t : R.tab) if (!t.tm.empty()) rockTm_ = true;
    if (R.tab.empty()) for (const auto& t : R.tab2) if (!t.tm.empty()) rockTm_ = true;
    log_(std::string("Rock compaction: ") + (R.tab.empty() ? "ROCK2D water-induced" : "ROCKTAB") + (rockTm_ ? " with transmissibility multipliers" : "") +
         (R.hyster ? ", hysteresis (ROCKTABH elastic curves)" : R.irreversible ? ", irreversible" : ", reversible") + (R.overburd.empty() ? "" : ", OVERBURD") + (R.store ? ", ROCKOPTS STORE" : ""));
}

void Simulator::updateRockHistory() {   // after each converged time step (OPM FlowProblem::updateMaxWaterSaturation / updateMinPressure)
    if (!rockOn_) return;
    const int na = grid_.numActive();
    for (int a = 0; a < na; ++a) {
        if (!rockPmin_.empty()) rockPmin_[size_t(a)] = std::min(rockPmin_[size_t(a)], p_[size_t(a)]);
        if (!rockSwMax_.empty() && (iSw_ >= 0 || W_)) rockSwMax_[size_t(a)] = std::max(rockSwMax_[size_t(a)], sw_[size_t(a)]);
    }
}

AD3 Simulator::rockMult(int c, const AD3& P, const AD3& Sw, bool trans) const {
    const RockCompInput& R = case_.props.rockComp;
    const int r = rockReg_[size_t(c)];
    const bool above = !rockPmin_.empty() && P.v > rockPmin_[size_t(c)];
    const double off = rockPob_[size_t(c)] + (rockP0_.empty() ? 0.0 : rockP0_[size_t(c)]);
    AD3 pe = (above && !R.hyster) ? AD3(rockPmin_[size_t(c)]) : P;
    pe -= off;
    auto lin = [](const std::vector<double>& X, const std::vector<double>& Y, const AD3& x) -> AD3 {
        if (X.empty() || Y.size() != X.size()) return AD3(1.0);
        if (X.size() == 1) return AD3(Y[0]);
        size_t i = size_t(std::upper_bound(X.begin(), X.end(), x.v) - X.begin()); i = std::clamp<size_t>(i, 1, X.size() - 1);
        const double s = (Y[i] - Y[i - 1]) / std::max(X[i] - X[i - 1], 1e-30);
        return (x - X[i - 1]) * s + Y[i - 1];   // linear extrapolation beyond the ends (OPM)
    };
    if (above && R.hyster && !R.tabH.empty()) {
        // C41 ROCKCOMP HYSTER: above the turning pressure tp (lowest pressure reached) the elastic reload curves,
        // interpolated as OPM UniformXTabulated2DFunction LeftExtreme: the two curves i, i+1 whose start pressures x_i,
        // x_i+1 bracket tp, alpha = (tp - x_i) / (x_i+1 - x_i), each curve read at p shifted parallel to the start
        // pressures (p - alpha dx on curve i, p + (1 - alpha) dx on curve i+1) - continuous with the deflation curve at p = tp
        const auto& curves = R.tabH[std::min(size_t(r), R.tabH.size() - 1)];
        const double tp = rockPmin_[size_t(c)] - off;
        if (curves.size() == 1) return lin(curves[0].p, trans ? curves[0].tm : curves[0].pv, pe);
        size_t i = 0;
        while (i + 2 < curves.size() && curves[i + 1].p[0] <= tp) ++i;
        const double dx = curves[i + 1].p[0] - curves[i].p[0];
        const double al = dx > 1e-30 ? (tp - curves[i].p[0]) / dx : 0.0;
        const AD3 s1 = lin(curves[i].p, trans ? curves[i].tm : curves[i].pv, pe - al * dx);
        const AD3 s2 = lin(curves[i + 1].p, trans ? curves[i + 1].tm : curves[i + 1].pv, pe + (1.0 - al) * dx);
        return s1 * (1.0 - al) + s2 * al;
    }
    if (!R.tab.empty()) {
        const auto& t = R.tab[std::min(size_t(r), R.tab.size() - 1)];
        return lin(t.p, trans ? t.tm : t.pv, pe);
    }
    const auto& t = R.tab2[std::min(size_t(r), R.tab2.size() - 1)];
    const auto& rows = trans ? t.tm : t.pv;
    if (rows.empty() || t.p.empty() || R.swNodes.empty()) return AD3(1.0);
    const auto& S = R.swNodes[std::min(size_t(r), R.swNodes.size() - 1)];
    AD3 ds = (!rockSwMax_.empty() && Sw.v < rockSwMax_[size_t(c)]) ? AD3(rockSwMax_[size_t(c)]) : Sw;
    if (!rockSw0_.empty()) ds -= rockSw0_[size_t(c)];
    // bilinear (OPM TabulatedTwoDFunction, vertical policy): along the ROCKWNOD axis in the two bracketing pressure rows
    const size_t np = std::min(t.p.size(), rows.size());
    auto rowVal = [&](size_t i) { std::vector<double> y = rows[i]; y.resize(S.size(), y.empty() ? 1.0 : y.back()); return lin(S, y, ds); };
    if (np == 1) return rowVal(0);
    size_t i = size_t(std::upper_bound(t.p.begin(), t.p.begin() + long(np), pe.v) - t.p.begin()); i = std::clamp<size_t>(i, 1, np - 1);
    const AD3 y0 = rowVal(i - 1), y1 = rowVal(i);
    const double dp = std::max(t.p[i] - t.p[i - 1], 1e-30);
    return y0 + (y1 - y0) * ((pe - t.p[i - 1]) / dp);
}

// C28 distributed assembly: cells whose properties this rank needs = owned cells, their neighbours (flux terms of the
// owned rows), every well-completion cell (the well rows are assembled on all ranks) and the aquifer cells
void Simulator::daNeedCells() {
    const int na = grid_.numActive();
    daNeed_.assign(size_t(na), 0);
    const int a0 = std::max(0, da_.a0), a1 = std::min(na, da_.a1);
    for (int a = a0; a < a1; ++a) {
        daNeed_[size_t(a)] = 1;
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) daNeed_[size_t(grid_.adjNbr[e])] = 1;
    }
    for (const auto& w : wells_) for (int c : w.cell) if (c >= 0 && c < na) daNeed_[size_t(c)] = 1;
    for (const auto& q : aquifers_) for (int c : q.cell) if (c >= 0 && c < na) daNeed_[size_t(c)] = 1;
}
std::vector<int> Simulator::cellPartition(int parts) const {
    const int na = grid_.numActive();
    parts = std::max(1, std::min(parts, std::max(1, na)));
    std::vector<int> st(size_t(parts) + 1, na);
    st[0] = 0;
    const double total = double(grid_.adjStart[na]) + 2.0 * na;    // block-row nnz (+ diagonal) + vector work per row
    int p = 1;
    for (int a = 0; a < na && p < parts; ++a) {
        const double done = double(grid_.adjStart[a + 1]) + 2.0 * (a + 1);
        if (done >= total * p / parts) st[size_t(p++)] = a + 1;
    }
    for (int q = 1; q <= parts; ++q) st[size_t(q)] = std::max(st[size_t(q)], st[size_t(q) - 1]);
    st[size_t(parts)] = na;
    return st;
}

// C34 directional relative permeability: mobilities of the cells whose KRNUMX / Y / Z table differs from SATNUM
void Simulator::evalDirKr(bool restricted) {
    const int na = grid_.numActive();
    const bool useDa = restricted && daOn_ && daNeed_.size() == size_t(na);
    for (int d = 0; d < 3; ++d) {
        if (krReg_[d].empty()) continue;
        cfDir_[d].resize(size_t(na));
        if (solv_) cfDirS_[d].resize(size_t(na)); else cfDirS_[d].clear();
        #pragma omp parallel for schedule(static) if (na > 2000 && useParallel(true))
        for (int c = 0; c < na; ++c) {
            if (useDa && !daNeed_[size_t(c)]) continue;
            const int r = krReg_[d][size_t(c)];
            if (r == satReg_[size_t(c)]) {
                cfDir_[d][size_t(c)] = cf_[size_t(c)];
                if (solv_) cfDirS_[d][size_t(c)] = cfS_[size_t(c)][0];
            } else {
                const CellFull f = evalCell(c, r);
                cfDir_[d][size_t(c)] = static_cast<const CellFlux&>(f);
                if (solv_) cfDirS_[d][size_t(c)] = f.mobS;
            }
        }
    }
}

void Simulator::evalAll(bool restricted) {
    const int na = grid_.numActive();
    if (solv_ && cfS_.size() != size_t(na)) cfS_.assign(size_t(na), std::array<AD3, 2>{});
    else if (!solv_ && !cfS_.empty()) std::vector<std::array<AD3, 2>>().swap(cfS_);
    if (restricted && daOn_ && daNeed_.size() == size_t(na)) {   // C28 distributed assembly: owned + halo + well cells
        if (useParallel(na > 2000)) {
            #pragma omp parallel for schedule(static)
            for (int c = 0; c < na; ++c) if (daNeed_[size_t(c)]) { const CellFull f = evalCell(c); storeCell(c, f); if (gravdr_) evalGravityHeights(c, f); }
        } else {
            for (int c = 0; c < na; ++c) if (daNeed_[size_t(c)]) { const CellFull f = evalCell(c); storeCell(c, f); if (gravdr_) evalGravityHeights(c, f); }
        }
        if (dirKr_) evalDirKr(true);
        return;
    }
    if (gravdr_) {
        if (useParallel(na > 2000)) {
            #pragma omp parallel for schedule(static)
            for (int c = 0; c < na; ++c) { const CellFull f = evalCell(c); storeCell(c, f); evalGravityHeights(c, f); }
        } else {
            for (int c = 0; c < na; ++c) { const CellFull f = evalCell(c); storeCell(c, f); evalGravityHeights(c, f); }
        }
        if (dirKr_) evalDirKr(false);
        return;
    }
    if (useParallel(na > 2000)) {
        #pragma omp parallel for schedule(static)
        for (int c = 0; c < na; ++c) storeCell(c, evalCell(c));
    } else {
        for (int c = 0; c < na; ++c) storeCell(c, evalCell(c));
    }
    if (dirKr_) evalDirKr(false);
}

// ---- GRAVDR: ECLIPSE dual porosity gravity drainage. The fluids in a matrix block of height Z and in the fracture around
// it are taken as vertically segregated; a column of water of height Hw = Z * Swn (Swn = (Sw - Swcr) / (1 - Swcr - Sowcr))
// and of gas Hg = Z * Sgn (Sgn = (Sg - Sgcr) / (1 - Swl - Sgcr - Sogcr)), clamped to [0, 1], adds a gravity head to the
// matrix-fracture potential difference of that phase:
//     dPhi_w(a->b) += (rho_w - rho_o) g Z (Hwn_a - Hwn_b),   dPhi_g(a->b) += (rho_o - rho_g) g Z (Hgn_a - Hgn_b)
// so water imbibes from below and gas enters from above when the fracture column is higher than the matrix one.
void Simulator::setupGravityDrainage(const Case& c) {
    gravdr_ = false; gdZ_.clear(); gdH_.clear();
    const GridInput& g = c.grid;
    if (!g.dual || !g.gravdr || !O_ || (!W_ && !G_)) return;
    gdZ_.assign(grid_.conns.size(), 0.0);
    int nset = 0, nDefault = 0;
    for (size_t n = 0; n < grid_.conns.size(); ++n) {
        const Connection& cn = grid_.conns[n];
        if (cn.dir != 'M') continue;
        const int qm = grid_.globalOf[size_t(cn.a)];
        double Z = size_t(qm) < g.dzmtrxv.size() ? g.dzmtrxv[size_t(qm)] : 0.0;
        if (!(Z > 0)) { Z = grid_.cellDz[size_t(qm)]; ++nDefault; }
        gdZ_[n] = Z; ++nset;
    }
    if (!nset) { gdZ_.clear(); return; }
    gdH_.assign(size_t(grid_.numActive()) * 2, AD3(0.0));
    gravdr_ = true;
    log_(fmt("GRAVDR: gravity drainage on %g matrix-fracture connections", nset) +
         (nDefault ? fmt(" (%g without DZMTRX/DZMTRXV use the cell thickness)", nDefault) : std::string()));
}

void Simulator::evalGravityHeights(int c, const CellFull& f) {
    const SatFunctions::EndPoints* E = ep(c);
    const SatFunctions::EndPoints& e = E ? *E : st(c).tableEndPoints();
    auto clamp01 = [](AD3 v) { if (v.v < 0) return AD3(0.0); if (v.v > 1) return AD3(1.0); return v; };
    AD3 hw(0.0), hg(0.0);
    if (W_) { const double d = 1.0 - e.swcr - e.sowcr; if (d > 1e-9) hw = clamp01((f.S[PW] - e.swcr) * (1.0 / d)); }
    if (G_) { const double d = 1.0 - e.swl - e.sgcr - e.sogcr; if (d > 1e-9) hg = clamp01((f.S[PG] - e.sgcr) * (1.0 / d)); }
    gdH_[size_t(c) * 2] = hw; gdH_[size_t(c) * 2 + 1] = hg;
}

// ======================================================================== assembly
void Simulator::assemble(double dt) {
    const int na = grid_.numActive(), nv = nv_;
    std::fill(J_.val.begin(), J_.val.end(), 0.0);
    if (blkAsm_) schur_.zeroValues();
    std::fill(R_.begin(), R_.end(), 0.0);
    const double g = U_.gravity;
    const double invDt = 1.0 / dt;
    bool ph[3] = {W_, O_, G_};
    if (adjCollect_) adjFace_.assign(grid_.conns.size() * size_t(nv), 0.0);

    const int aLo = daOn_ ? std::max(0, da_.a0) : 0, aHi = daOn_ ? std::min(na, da_.a1) : na;   // C28: owned rows only
    if (useParallel(na > 2000)) {
        #pragma omp parallel for schedule(dynamic, 256)
        for (int a = aLo; a < aHi; ++a) {
            const CellFlux& A = cf_[a];
            const int self = selfPos_[a];
            auto addSelf = [&](int e, const AD3& v, double scale) {
                const int r = a * nv + e;
                R_[r] += v.v * scale;
                double* row = jrow(a, e, self);
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
                const int kd = dirKr_ ? (cn.dir == 'X' ? 0 : cn.dir == 'Y' ? 1 : cn.dir == 'Z' ? 2 : -1) : -1;   // C34 KRNUMX/Y/Z
                const CellFlux& MA = kd >= 0 && !cfDir_[kd].empty() ? cfDir_[kd][size_t(a)] : A;
                const CellFlux& MB = kd >= 0 && !cfDir_[kd].empty() ? cfDir_[kd][size_t(b)] : B;
                const int nbp = adjPos_[e];
                // C15 stress-dependent permeability: T * exp(-(gamma_a (p_ia - p_a) + gamma_b (p_ib - p_b)) / 2)
                AD6 smul(1.0);
                const bool sOn = !sGam_.empty() && !pI_.empty() && (sGam_[size_t(a)] > 0 || sGam_[size_t(b)] > 0);
                if (sOn) {
                    AD6 ex = (lift<6>(AD3::variable(p_[size_t(a)], 0), 0) - pI_[size_t(a)]) * (0.5 * sGam_[size_t(a)]) +
                             (lift<6>(AD3::variable(p_[size_t(b)], 0), 3) - pI_[size_t(b)]) * (0.5 * sGam_[size_t(b)]);
                    if (ex.v > 3.0) ex = AD6(3.0); else if (ex.v < -20.0) ex = AD6(-20.0);
                    smul = exp(ex);
                }
                auto add = [&](int eqi, const AD6& F) {
                    const int r = a * nv + eqi;
                    R_[r] += F.v;
                    double* rs = jrow(a, eqi, self);
                    double* rn = jrow(a, eqi, nbp);
                    for (int k = 0; k < nv; ++k) { rs[k] += F.d[k]; rn[k] += F.d[3 + k]; }
                };
                for (int p = 0; p < 3; ++p) if (ph[p]) {
                    const AD6 rho = (lift<6>(A.rho[p], 0) + lift<6>(B.rho[p], 3)) * 0.5;
                    AD6 dphi = lift<6>(A.P[p], 0) - lift<6>(B.P[p], 3) + rho * (g * dD);
                    if (gravdr_ && p != PO && gdZ_[size_t(grid_.adjConn[e])] > 0) {   // GRAVDR column heads (matrix <-> fracture)
                        const double Z = gdZ_[size_t(grid_.adjConn[e])];
                        const int hp = p == PW ? 0 : 1;
                        const AD6 rhoO = (lift<6>(A.rho[PO], 0) + lift<6>(B.rho[PO], 3)) * 0.5;
                        const AD6 drho = p == PW ? rho - rhoO : rhoO - rho;
                        dphi += drho * (lift<6>(gdH_[size_t(a) * 2 + hp], 0) - lift<6>(gdH_[size_t(b) * 2 + hp], 3)) * (g * Z);
                    }
                    if (!thp_.empty()) {
                        const double t = thp_[size_t(grid_.adjConn[e])];
                        if (t > 0) {
                            if (dphi.v > t) dphi -= t; else if (dphi.v < -t) dphi += t; else dphi = AD6(0.0);
                        }
                    }
                    const bool upA = dphi.v >= 0;
                    const AD6 mob = upA ? lift<6>(MA.mob[p], 0) : lift<6>(MB.mob[p], 3);
                    const double tr = (p == PW && !pcZw_.empty()) ? cn.trans / pcZw_[size_t(grid_.adjConn[e])] : cn.trans;
                    const AD6 F = sOn ? mob * dphi * smul * tr : mob * dphi * tr;
                    add(eq_[p], F);
                    double* fr = (adjCollect_ && cn.a == a) ? &adjFace_[size_t(grid_.adjConn[e]) * nv] : nullptr;
                    if (fr) fr[eq_[p]] += F.v;
                    if (p == PO && live_) {
                        const AD6 mr = upA ? lift<6>(MA.mobRs, 0) : lift<6>(MB.mobRs, 3);
                        const AD6 Fd = sOn ? mr * dphi * smul * cn.trans : mr * dphi * cn.trans;
                        add(eq_[PG], Fd);
                        if (fr) fr[eq_[PG]] += Fd.v;
                    }
                    if (p == PG && vap_) {
                        const AD6 mv = upA ? lift<6>(MA.mobRv, 0) : lift<6>(MB.mobRv, 3);
                        const AD6 Fv = sOn ? mv * dphi * smul * cn.trans : mv * dphi * cn.trans;
                        add(eq_[PO], Fv);
                        if (fr) fr[eq_[PO]] += Fv.v;
                    }
                }
                if (solv_) {   // C35 SOLVENT flux in the G+S equation: gas-phase pressure, solvent density, upstream mobility
                    const AD6 rhoS = (lift<6>(cfS_[size_t(a)][1], 0) + lift<6>(cfS_[size_t(b)][1], 3)) * 0.5;
                    AD6 dphi = lift<6>(A.P[PG], 0) - lift<6>(B.P[PG], 3) + rhoS * (g * dD);
                    if (!thp_.empty()) {
                        const double t = thp_[size_t(grid_.adjConn[e])];
                        if (t > 0) { if (dphi.v > t) dphi -= t; else if (dphi.v < -t) dphi += t; else dphi = AD6(0.0); }
                    }
                    const AD6 ms = dphi.v >= 0 ? lift<6>(mobSDir(kd, a), 0) : lift<6>(mobSDir(kd, b), 3);
                    const AD6 F = sOn ? ms * dphi * smul * cn.trans : ms * dphi * cn.trans;
                    add(eq_[PG], F);
                }
            }
        }
    } else {
        for (int a = aLo; a < aHi; ++a) {
            const CellFlux& A = cf_[a];
            const int self = selfPos_[a];
            auto addSelf = [&](int e, const AD3& v, double scale) {
                const int r = a * nv + e;
                R_[r] += v.v * scale;
                double* row = jrow(a, e, self);
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
                const int kd = dirKr_ ? (cn.dir == 'X' ? 0 : cn.dir == 'Y' ? 1 : cn.dir == 'Z' ? 2 : -1) : -1;   // C34 KRNUMX/Y/Z
                const CellFlux& MA = kd >= 0 && !cfDir_[kd].empty() ? cfDir_[kd][size_t(a)] : A;
                const CellFlux& MB = kd >= 0 && !cfDir_[kd].empty() ? cfDir_[kd][size_t(b)] : B;
                const int nbp = adjPos_[e];
                // C15 stress-dependent permeability: T * exp(-(gamma_a (p_ia - p_a) + gamma_b (p_ib - p_b)) / 2)
                AD6 smul(1.0);
                const bool sOn = !sGam_.empty() && !pI_.empty() && (sGam_[size_t(a)] > 0 || sGam_[size_t(b)] > 0);
                if (sOn) {
                    AD6 ex = (lift<6>(AD3::variable(p_[size_t(a)], 0), 0) - pI_[size_t(a)]) * (0.5 * sGam_[size_t(a)]) +
                             (lift<6>(AD3::variable(p_[size_t(b)], 0), 3) - pI_[size_t(b)]) * (0.5 * sGam_[size_t(b)]);
                    if (ex.v > 3.0) ex = AD6(3.0); else if (ex.v < -20.0) ex = AD6(-20.0);
                    smul = exp(ex);
                }
                auto add = [&](int eqi, const AD6& F) {
                    const int r = a * nv + eqi;
                    R_[r] += F.v;
                    double* rs = jrow(a, eqi, self);
                    double* rn = jrow(a, eqi, nbp);
                    for (int k = 0; k < nv; ++k) { rs[k] += F.d[k]; rn[k] += F.d[3 + k]; }
                };
                for (int p = 0; p < 3; ++p) if (ph[p]) {
                    const AD6 rho = (lift<6>(A.rho[p], 0) + lift<6>(B.rho[p], 3)) * 0.5;
                    AD6 dphi = lift<6>(A.P[p], 0) - lift<6>(B.P[p], 3) + rho * (g * dD);
                    if (gravdr_ && p != PO && gdZ_[size_t(grid_.adjConn[e])] > 0) {   // GRAVDR column heads (matrix <-> fracture)
                        const double Z = gdZ_[size_t(grid_.adjConn[e])];
                        const int hp = p == PW ? 0 : 1;
                        const AD6 rhoO = (lift<6>(A.rho[PO], 0) + lift<6>(B.rho[PO], 3)) * 0.5;
                        const AD6 drho = p == PW ? rho - rhoO : rhoO - rho;
                        dphi += drho * (lift<6>(gdH_[size_t(a) * 2 + hp], 0) - lift<6>(gdH_[size_t(b) * 2 + hp], 3)) * (g * Z);
                    }
                    if (!thp_.empty()) {
                        const double t = thp_[size_t(grid_.adjConn[e])];
                        if (t > 0) {
                            if (dphi.v > t) dphi -= t; else if (dphi.v < -t) dphi += t; else dphi = AD6(0.0);
                        }
                    }
                    const bool upA = dphi.v >= 0;
                    const AD6 mob = upA ? lift<6>(MA.mob[p], 0) : lift<6>(MB.mob[p], 3);
                    const double tr = (p == PW && !pcZw_.empty()) ? cn.trans / pcZw_[size_t(grid_.adjConn[e])] : cn.trans;
                    const AD6 F = sOn ? mob * dphi * smul * tr : mob * dphi * tr;
                    add(eq_[p], F);
                    double* fr = (adjCollect_ && cn.a == a) ? &adjFace_[size_t(grid_.adjConn[e]) * nv] : nullptr;
                    if (fr) fr[eq_[p]] += F.v;
                    if (p == PO && live_) {
                        const AD6 mr = upA ? lift<6>(MA.mobRs, 0) : lift<6>(MB.mobRs, 3);
                        const AD6 Fd = sOn ? mr * dphi * smul * cn.trans : mr * dphi * cn.trans;
                        add(eq_[PG], Fd);
                        if (fr) fr[eq_[PG]] += Fd.v;
                    }
                    if (p == PG && vap_) {
                        const AD6 mv = upA ? lift<6>(MA.mobRv, 0) : lift<6>(MB.mobRv, 3);
                        const AD6 Fv = sOn ? mv * dphi * smul * cn.trans : mv * dphi * cn.trans;
                        add(eq_[PO], Fv);
                        if (fr) fr[eq_[PO]] += Fv.v;
                    }
                }
                if (solv_) {   // C35 SOLVENT flux in the G+S equation: gas-phase pressure, solvent density, upstream mobility
                    const AD6 rhoS = (lift<6>(cfS_[size_t(a)][1], 0) + lift<6>(cfS_[size_t(b)][1], 3)) * 0.5;
                    AD6 dphi = lift<6>(A.P[PG], 0) - lift<6>(B.P[PG], 3) + rhoS * (g * dD);
                    if (!thp_.empty()) {
                        const double t = thp_[size_t(grid_.adjConn[e])];
                        if (t > 0) { if (dphi.v > t) dphi -= t; else if (dphi.v < -t) dphi += t; else dphi = AD6(0.0); }
                    }
                    const AD6 ms = dphi.v >= 0 ? lift<6>(mobSDir(kd, a), 0) : lift<6>(mobSDir(kd, b), 3);
                    const AD6 F = sOn ? ms * dphi * smul * cn.trans : ms * dphi * cn.trans;
                    add(eq_[PG], F);
                }
            }
        }
    }

    // ---- wells (serial; few connections)
    const int refP = O_ ? PO : PW;
    if (adjCollect_) { adjRd_.clear(); adjCd_.clear(); }
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
        if (!W.open || (W.mode != Control::BHP && target == 0.0)) {
            const bool stopSt = !W.open && W.ctrl.status == WellStatus::Stop && W.nActive > 0;   // C49 STOP status
            if ((W.open || stopSt) && !W.cell.empty()) {   // C33 stopped well (open, zero target): report the zero-inflow BHP like OPM / ECLIPSE
                double sw = 0, sp = 0;
                for (size_t k = 0; k < W.cell.size(); ++k) {
                    if (W.wi[k] <= 0) continue;
                    sw += W.wi[k]; sp += W.wi[k] * (p_[W.cell[k]] - (W.pConn(k, g) - W.bhp));
                }
                if (sw > 0) { W.bhp = sp / sw; rt.bhp = W.bhp; }
            }
            J_.val[bhpPos] = 1.0; R_[wr] = 0.0; W.rates = rt; continue;
        }
        const bool prod = W.type == WellType::Producer;
        const double scale = (W.mode == Control::BHP) ? 1.0 : 1.0 / std::max(std::fabs(target), 1.0);
        double qsel = 0;
        double xfTot[3] = {0, 0, 0};   // C41 crossflow (surface oil, water, gas into the reservoir, < 0)
        double recNew[3] = {0, 0, 0}, injTotNew = 0;   // C44 injector crossflow recirculation (next assembly)
        const bool isMswW = mswOfWell(w) >= 0;
        const bool xfOk = W.xflow && (!isMswW || mswXflow());   // C45: multisegment wells cross-flow too when enabled (segment mixtures)
        std::vector<double> gIn(isMswW && prod && xfOk ? 3 * W.cell.size() : 0, 0.0);   // C45 gross inflow per connection
        for (size_t k = 0; k < W.cell.size(); ++k) {
            const int a = W.cell[k];
            const CellFull C = evalCell(a);
            const AD4 bhp = AD4::variable(W.bhp, 3);
            const AD4 pconn = bhp + (W.pConn(k, g) - W.bhp);
            const AD4 Pc = lift<4>(C.P[refP], 0);
            const double wi = W.wi[k];
            const int nb = blkStart_[a + 1] - blkStart_[a];
            const double ef = W.eff;                          // C33 WEFAC: the reservoir sees eff x the flowing rate
            auto addCell = [&](int eqi, const AD4& q) {      // q: outflow (+) from the cell
                const int r = a * nv + eqi;
                R_[r] += ef * q.v;
                double* row = jrow(a, eqi, selfPos_[a]);
                for (int v = 0; v < nv; ++v) row[v] += ef * q.d[v];
                J_.val[J_.rowPtr[r] + (blkAsm_ ? 0 : nb * nv) + W.cellWellIdx[k]] += ef * q.d[3];
            };
            AD4 sel(0.0);
            if (prod) {
                const AD4 dd = Pc - pconn;
                AD4 q[3]; AD4 qd(0.0), qv(0.0);   // qd: dissolved gas (in oil), qv: vaporised oil (in gas)
                for (int p = 0; p < 3; ++p) if (ph[p]) q[p] = lift<4>(C.mob[p], 0) * dd * wi;
                if (ph[PW] && k < W.wmulConn.size()) q[PW] *= W.wmulConn[k];
                if (live_) qd = lift<4>(C.mobRs, 0) * dd * wi;
                if (vap_) qv = lift<4>(C.mobRv, 0) * dd * wi;
                const AD4 qs = solv_ ? lift<4>(C.mobS, 0) * dd * wi : AD4(0.0);   // C35 solvent (G+S equation)
                if (dd.v > 0) {
                    for (int p = 0; p < 3; ++p) if (ph[p]) addCell(eq_[p], p == PG ? q[p] + qd + qs : p == PO ? q[p] + qv : q[p]);
                    if (solv_) { rt.solvProd += qs.v; if (C.bs.v > 0) rt.resv += qs.v / C.bs.v; }
                    if (adjCollect_) {   // adjoint: rate derivatives w.r.t. the cell variables and the BHP
                        AdjRateD rd; rd.well = w; rd.cell = a; rd.k = int(k);
                        for (int v = 0; v < 4; ++v) {
                            rd.d[0][v] = (ph[PO] ? q[PO].d[v] : 0.0) + qv.d[v];
                            rd.d[1][v] = ph[PW] ? q[PW].d[v] : 0.0;
                            rd.d[2][v] = (ph[PG] ? q[PG].d[v] : 0.0) + qd.d[v];
                            rd.d[3][v] = rd.d[4][v] = 0.0;
                        }
                        rd.v[0] = (ph[PO] ? q[PO].v : 0.0) + qv.v; rd.v[1] = ph[PW] ? q[PW].v : 0.0;
                        rd.v[2] = (ph[PG] ? q[PG].v : 0.0) + qd.v;
                        adjRd_.push_back(rd);
                    }
                    rt.oil += (ph[PO] ? q[PO].v : 0) + qv.v; rt.water += ph[PW] ? q[PW].v : 0;
                    if (!gIn.empty()) { gIn[3 * k] = (ph[PO] ? q[PO].v : 0) + qv.v; gIn[3 * k + 1] = ph[PW] ? q[PW].v : 0; gIn[3 * k + 2] = (ph[PG] ? q[PG].v : 0) + qd.v; }
                    if (ph[PW]) W.qwConn[k] = q[PW].v;
                    rt.gas += (ph[PG] ? q[PG].v : 0) + qd.v;
                    double resv = 0;
                    if (W.resvOk) resv = W.resvC[0] * (ph[PW] ? q[PW].v : 0) + W.resvC[1] * ((ph[PO] ? q[PO].v : 0) + qv.v) + W.resvC[2] * ((ph[PG] ? q[PG].v : 0) + qd.v);
                    else {
                        if (ph[PO]) resv += q[PO].v / C.b[PO].v;
                        if (ph[PW]) resv += q[PW].v / C.b[PW].v;
                        if (ph[PG]) resv += q[PG].v / C.b[PG].v;
                    }
                    rt.resv += resv;
                }
                // C41 crossflow (WELSPECS item 10, default YES; OPM StandardWell::computePerfRate): a connection whose
                // cell pressure is below the connection pressure takes the wellbore mixture (surface fractions x_c of the
                // well's production at the last assembly): q_t,res = lambda_t WI dd, q_c = x_c q_t,res / V with the volume
                // ratio V = x_w/b_w + (x_o - Rv x_g)/(d b_o) + (x_g - Rs x_o)/(d b_g), d = 1 - Rs Rv.
                // C45 multisegment wells: the mixture of the connection's segment (OPM MultisegmentWell: segment fractions)
                const double* xmk = isMswW && W.xmixConn.size() == 3 * W.cell.size() && (W.xmixConn[3 * k] + W.xmixConn[3 * k + 1] + W.xmixConn[3 * k + 2]) > 0
                                        ? &W.xmixConn[3 * k] : W.xmix;
                const bool xf = xfOk && dd.v < 0 && (xmk[0] + xmk[1] + xmk[2]) > 0;
                AD4 xc[3] = {AD4(0.0), AD4(0.0), AD4(0.0)};
                if (xf) {
                    AD4 lt(0.0);
                    for (int p = 0; p < 3; ++p) if (ph[p]) lt += lift<4>(C.mob[p], 0) / lift<4>(C.b[p], 0) * (p == PW && k < W.wmulConn.size() ? W.wmulConn[k] : 1.0);
                    const double xo = ph[PO] ? xmk[0] : 0.0, xw = ph[PW] ? xmk[1] : 0.0, xg = ph[PG] ? xmk[2] : 0.0;
                    AD4 vr(0.0);
                    if (ph[PW]) vr += xw / lift<4>(C.b[PW], 0);
                    const AD4 rs = live_ ? lift<4>(C.rs, 0) : AD4(0.0), rv = vap_ ? lift<4>(C.rv, 0) : AD4(0.0);
                    if (ph[PO] && ph[PG]) {
                        const AD4 d = 1.0 - rs * rv;
                        const bool dok = d.v > 1e-12;
                        vr += (dok ? (xo - rv * xg) / d : AD4(xo)) / lift<4>(C.b[PO], 0);
                        vr += (dok ? (xg - rs * xo) / d : AD4(xg)) / lift<4>(C.b[PG], 0);
                    } else {
                        if (ph[PO]) vr += xo / lift<4>(C.b[PO], 0);
                        if (ph[PG]) vr += xg / lift<4>(C.b[PG], 0);
                    }
                    if (vr.v > 1e-30) {
                        const AD4 qts = lt * dd * wi / vr;   // < 0: into the cell, surface volume of mixture
                        xc[0] = xo * qts; xc[1] = xw * qts; xc[2] = xg * qts;
                        if (ph[PO]) addCell(eq_[PO], xc[0]);
                        if (ph[PW]) addCell(eq_[PW], xc[1]);
                        if (ph[PG]) addCell(eq_[PG], xc[2]);
                        rt.oil += xc[0].v; rt.water += xc[1].v; rt.gas += xc[2].v;
                        xfTot[0] += xc[0].v; xfTot[1] += xc[1].v; xfTot[2] += xc[2].v;
                        if (ph[PW]) W.qwConn[k] = xc[1].v;
                        double resv = 0;
                        if (W.resvOk) resv = W.resvC[0] * xc[1].v + W.resvC[1] * xc[0].v + W.resvC[2] * xc[2].v;
                        else resv = lt.v * dd.v * wi;
                        rt.resv += resv;
                    }
                }
                // the control equation counts the connections that flow (dd > 0, as in addCell and the reported rates) plus
                // the crossflow connections; counting blocked back-flow connections made rate targets look unreachable
                // (switch to BHP) and then over-met (switch back) - an endless control flip-flop on Norne D-2H
                if (xf) switch (W.mode) {
                    case Control::ORAT: sel = xc[0]; break;
                    case Control::WRAT: sel = xc[1]; break;
                    case Control::GRAT: sel = xc[2]; break;
                    case Control::LRAT: sel = xc[0] + xc[1]; break;
                    case Control::RESV: if (W.resvOk) sel = xc[1] * W.resvC[0] + xc[0] * W.resvC[1] + xc[2] * W.resvC[2]; break;
                    default: break;
                }
                if (dd.v > 0) switch (W.mode) {
                    case Control::ORAT: sel = q[PO] + qv; break;
                    case Control::WRAT: sel = q[PW]; break;
                    case Control::GRAT: sel = q[PG] + qd; break;
                    case Control::LRAT: sel = q[PO] + qv + q[PW]; break;
                    case Control::RESV: if (W.resvOk) {
                        if (ph[PW]) sel += q[PW] * W.resvC[0];
                        if (ph[PO]) sel += (q[PO] + qv) * W.resvC[1]; else sel += qv * W.resvC[1];
                        if (ph[PG]) sel += (q[PG] + qd) * W.resvC[2]; else sel += qd * W.resvC[2];
                    } else {
                        if (ph[PO]) sel += q[PO] / lift<4>(C.b[PO], 0);
                        if (ph[PW]) sel += q[PW] / lift<4>(C.b[PW], 0);
                        if (ph[PG]) sel += q[PG] / lift<4>(C.b[PG], 0);
                        if (solv_ && C.bs.v > 0) sel += qs / lift<4>(C.bs, 0);
                    } break;
                    default: break;
                }
            } else {
                const AD4 dd = pconn - Pc;
                AD4 lt(0.0);
                for (int p = 0; p < 3; ++p) if (ph[p]) lt += lift<4>(C.mob[p], 0) / lift<4>(C.b[p], 0) * (p == PW && k < W.wmulConn.size() ? W.wmulConn[k] : 1.0);
                const int ip = (W.type == WellType::GasInjector) ? PG : PW;
                if (!ph[ip]) throw std::runtime_error("Injector phase is not active in the model");
                AD4 binj = lift<4>(C.b[ip], 0);
                double fS = 0;
                if (solv_) {   // C35: solvent mobility in the injectivity; injected gas with the WSOLVENT fraction fS
                    if (C.bs.v > 0) lt += lift<4>(C.mobS, 0) / lift<4>(C.bs, 0);
                    fS = ip == PG && size_t(w) < svInjF_.size() ? svInjF_[size_t(w)] : 0.0;
                    if (fS > 0 && C.bs.v > 0) binj = 1.0 / ((1.0 - fS) / lift<4>(C.b[PG], 0) + fS / lift<4>(C.bs, 0));
                }
                const AD4 q = lt * binj * dd * wi;
                if (dd.v > 0) {
                    addCell(eq_[ip], -q);
                    // C44: fluid back-produced at crossflow connections stays in the wellbore and leaves with the injection
                    // (ECLIPSE / OPM wellbore mixture): the non-injected phases re-enter through the injecting connections in
                    // proportion to their injection (lagged by one assembly; the injected phase is already netted)
                    if (W.xflow && W.xinjTot > 0) {
                        const double sh = q.v / W.xinjTot;
                        const int po[3] = {PO, PW, PG};
                        for (int i = 0; i < 3; ++i) if (po[i] != ip && ph[po[i]] && W.xrec[i] > 0) addCell(eq_[po[i]], AD4(-W.xrec[i] * sh));
                    }
                    injTotNew += q.v;
                    if (adjCollect_) {
                        AdjRateD rd; rd.well = w; rd.cell = a; rd.k = int(k);
                        for (int v = 0; v < 4; ++v) { for (int qi = 0; qi < 5; ++qi) rd.d[qi][v] = 0.0; rd.d[ip == PW ? 3 : 4][v] = q.d[v]; }
                        rd.v[ip == PW ? 3 : 4] = q.v;
                        adjRd_.push_back(rd);
                    }
                    if (ip == PW) { rt.waterInj += q.v; W.qwConn[k] = -q.v; } else { rt.gasInj += q.v; rt.solvInj += fS * q.v; }
                    rt.resv -= q.v / binj.v;
                    sel = q;
                } else if (xfOk && dd.v < 0) {
                    // C41 injector crossflow (WELSPECS item 10, OPM StandardWell::computePerfRate): a connection whose cell
                    // pressure exceeds the connection pressure produces the reservoir fluid at the cell mobilities; the
                    // rate control counts the net injected-phase surface rate (injection minus back-production)
                    const AD4 dp = -dd;
                    AD4 qp[3] = {AD4(0.0), AD4(0.0), AD4(0.0)};
                    for (int p = 0; p < 3; ++p) if (ph[p]) qp[p] = lift<4>(C.mob[p], 0) * dp * wi;
                    if (ph[PW] && k < W.wmulConn.size()) qp[PW] *= W.wmulConn[k];
                    const AD4 qd = live_ ? lift<4>(C.mobRs, 0) * dp * wi : AD4(0.0);
                    const AD4 qv = vap_ ? lift<4>(C.mobRv, 0) * dp * wi : AD4(0.0);
                    if (ph[PO]) addCell(eq_[PO], qp[PO] + qv);
                    if (ph[PW]) addCell(eq_[PW], qp[PW]);
                    if (ph[PG]) addCell(eq_[PG], qp[PG] + qd);
                    const AD4 qInj = ip == PW ? qp[PW] : qp[PG] + qd;
                    // C44: the non-injected phases are recirculated (not surface production) while the well has injecting
                    // connections; only a well that injects nowhere would lift them to surface
                    const double bo = (ph[PO] ? qp[PO].v : 0.0) + qv.v, bw = ph[PW] ? qp[PW].v : 0.0, bg = (ph[PG] ? qp[PG].v : 0.0) + qd.v;
                    const bool recirc = W.xinjTot > 0;
                    if (ip == PW) { rt.waterInj -= qp[PW].v; W.qwConn[k] = qp[PW].v; recNew[0] += bo; recNew[2] += bg; if (!recirc) { rt.oil += bo; rt.gas += bg; } }
                    else { rt.gasInj -= qInj.v; recNew[0] += bo; recNew[1] += bw; if (!recirc) { rt.oil += bo; rt.water += bw; } }
                    AD4 rv(0.0);
                    for (int p = 0; p < 3; ++p) if (ph[p]) rv += qp[p] / lift<4>(C.b[p], 0);
                    rt.resv += rv.v;
                    sel = -qInj;
                }
            }
            if (adjCollect_) {   // adjoint: d(residual)/d(connection head) = the connection's BHP derivatives
                AdjConnD cd; cd.well = w; cd.k = int(k); cd.cell = a;
                for (int e = 0; e < 3; ++e) cd.dc[e] = 0.0;
                if (prod) {
                    const AD4 dd = Pc - pconn;
                    if (dd.v > 0) {
                        AD4 qq[3]; AD4 qd(0.0), qv(0.0);
                        for (int p = 0; p < 3; ++p) if (ph[p]) qq[p] = lift<4>(C.mob[p], 0) * dd * wi;
                        if (ph[PW] && k < W.wmulConn.size()) qq[PW] *= W.wmulConn[k];
                        if (live_) qd = lift<4>(C.mobRs, 0) * dd * wi;
                        if (vap_) qv = lift<4>(C.mobRv, 0) * dd * wi;
                        for (int p = 0; p < 3; ++p) if (ph[p]) {
                            const AD4 t = p == PG ? qq[p] + qd : p == PO ? qq[p] + qv : qq[p];
                            cd.dc[eq_[p]] = t.d[3]; cd.vc[eq_[p]] = t.v;
                        }
                    }
                } else {
                    const AD4 dd = pconn - Pc;
                    if (dd.v > 0) {
                        AD4 lt(0.0);
                        for (int p = 0; p < 3; ++p) if (ph[p]) lt += lift<4>(C.mob[p], 0) / lift<4>(C.b[p], 0) * (p == PW && k < W.wmulConn.size() ? W.wmulConn[k] : 1.0);
                        const int ip = (W.type == WellType::GasInjector) ? PG : PW;
                        const AD4 t = lt * lift<4>(C.b[ip], 0) * dd * wi;
                        cd.dc[eq_[ip]] = -t.d[3]; cd.vc[eq_[ip]] = -t.v;
                    }
                }
                for (int e = 0; e < 3; ++e) { cd.dc[e] *= W.eff; cd.vc[e] *= W.eff; }   // C33 WEFAC
                cd.dw = W.mode != Control::BHP ? sel.d[3] * scale : 0.0;
                cd.vw = W.mode != Control::BHP ? sel.v * scale : 0.0;
                adjCd_.push_back(cd);
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
        if (!prod) { for (int i = 0; i < 3; ++i) W.xrec[i] = recNew[i]; W.xinjTot = injTotNew; }   // C44 recirculation state
        if (prod) {   // C41: wellbore mixture for the crossflow connections of the next assembly (gross production)
            const double go = std::max(0.0, rt.oil - xfTot[0]), gw = std::max(0.0, rt.water - xfTot[1]), gg = std::max(0.0, rt.gas - xfTot[2]);
            const double st = go + gw + gg;
            if (st > 0) { W.xmix[0] = go / st; W.xmix[1] = gw / st; W.xmix[2] = gg / st; }
            if (!gIn.empty()) mswMixture(W, gIn);
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
            double pd, dpd;
            if (Q.tdTab.size() >= 2) {   // C34 AQUTAB: piecewise linear pD(tD), linear extrapolation beyond the ends (OPM)
                const auto& X = Q.tdTab; const auto& Y = Q.pdTab;
                size_t i = size_t(std::upper_bound(X.begin(), X.end(), tn1) - X.begin());
                i = std::clamp<size_t>(i, 1, X.size() - 1);
                dpd = (Y[i] - Y[i - 1]) / std::max(X[i] - X[i - 1], 1e-30);
                pd = Y[i - 1] + dpd * (tn1 - X[i - 1]);
            } else carterTracyPD(tn1, pd, dpd);
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
            double* row = jrow(a, eq_[PW], selfPos_[a]);
            for (int v = 0; v < nv; ++v) row[v] -= src.d[v];
        }
    }

    // ---- C29 SALT dry-out: lagged evaporation sink E (surface water / day), implicit fade-out E S_w / (S_w + 0.005)
    if (!saltEvap_.empty() && eq_[PW] >= 0)
        for (int a = aLo; a < aHi; ++a) {
            const double E = saltEvap_[size_t(a)];
            if (!(E > 0)) continue;
            const CellFull C = evalCell(a);
            if (!(C.S[PW].v > 0)) continue;
            const AD3 src = C.S[PW] / (C.S[PW] + 0.005) * E;
            const int r = a * nv + eq_[PW];
            R_[r] += src.v;
            double* row = jrow(a, eq_[PW], selfPos_[a]);
            for (int v = 0; v < nv; ++v) row[v] += src.d[v];
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
                if (blkAsm_)
                    for (int k = blkStart_[a]; k < blkStart_[a + 1]; ++k) { double* bl = schur_.cellRow(size_t(k), e); for (int v = 0; v < nv; ++v) bl[v] *= s; }
            }
        }
    } else {
        for (int a = 0; a < na; ++a) {
            const double s = dt / grid_.poreVolRef[a];
            for (int e = 0; e < nv; ++e) {
                const int r = a * nv + e;
                R_[r] *= s;
                for (int q = J_.rowPtr[r]; q < J_.rowPtr[r + 1]; ++q) J_.val[q] *= s;
                if (blkAsm_)
                    for (int k = blkStart_[a]; k < blkStart_[a + 1]; ++k) { double* bl = schur_.cellRow(size_t(k), e); for (int v = 0; v < nv; ++v) bl[v] *= s; }
            }
        }
    }
}

bool Simulator::converged(double dt, double& cnvMax, double& mbMax) const {
    (void)dt;
    const int na = grid_.numActive(), nv = nv_;
    double cnv[3] = {0, 0, 0}, sum[3] = {0, 0, 0}, pvb[3] = {0, 0, 0};
    const int aLo = daOn_ ? std::max(0, da_.a0) : 0, aHi = daOn_ ? std::min(na, da_.a1) : na;
    for (int a = aLo; a < aHi; ++a) {
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
    if (daOn_) {   // C28: global norms (identical on every rank -> identical Newton decisions)
        double sv[6] = {sum[0], sum[1], sum[2], pvb[0], pvb[1], pvb[2]};
        if (da_.sumN) da_.sumN(sv, 6);
        if (da_.maxN) da_.maxN(cnv, 3);
        for (int p = 0; p < 3; ++p) { sum[p] = sv[p]; pvb[p] = sv[3 + p]; }
    }
    cnvMax = 0; mbMax = 0;
    for (int p = 0; p < 3; ++p) cnvPh_[p] = cnv[p];
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
        if (st == 1) x_[a] += f * dX;
        else if (st == 0) {
            const double lim = 0.25 * std::max(rsSatC(a, p_[a]), std::max(x_[a], 1e-6));
            x_[a] += std::max(-lim, std::min(lim, dX));
        } else {
            const double lim = 0.25 * std::max(rvSatC(a, p_[a]), std::max(x_[a], 1e-9));
            x_[a] += std::max(-lim, std::min(lim, dX));
        }
        if (!O_) { if (x_[a] < 0) x_[a] = 0; return; }
        const double sw = iSw_ >= 0 ? sw_[a] : 0.0;
        // without oil and gas (water only) the Rs / Rv variables are undefined: such cells use state 1 with Sg = 0
        const bool waterOnly = iSw_ >= 0 && sw > 1.0 - 1e-9;
        if (st == 1) {
            if (x_[a] < 0) {
                if (live_ && !waterOnly) { hasGas_[a] = 0; x_[a] = rsSatC(a, p_[a]); if (!rsLim_.empty()) x_[a] = std::min(x_[a], rsLim_[size_t(a)]); }
                else x_[a] = 0;
            } else if (sw + x_[a] > 1.0) {
                if (vap_ && !waterOnly) { hasGas_[a] = 2; x_[a] = rvSatC(a, p_[a]); }   // oil vanished: gas carries it
                else { const double s = sw + x_[a]; if (iSw_ >= 0) sw_[a] /= s; x_[a] /= s; }
            }
        } else if (st == 0) {
            if (x_[a] < 0) x_[a] = 0;
            const double rsCap = rsLim_.empty() ? rsSatC(a, p_[a]) : std::min(rsSatC(a, p_[a]), rsLim_[size_t(a)] * (1 + 1e-9) + 1e-12);
            if (waterOnly || x_[a] > rsCap) { hasGas_[a] = 1; x_[a] = 0.0; }
        } else {
            if (x_[a] < 0) x_[a] = 0;
            if (waterOnly) { hasGas_[a] = 1; x_[a] = 0.0; }
            else if (x_[a] > rvSatC(a, p_[a])) { hasGas_[a] = 1; x_[a] = std::max(0.0, 1.0 - sw); }   // condensate drops out
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
        // relative to the reservoir pressure scale as well: a BHP far below the cell pressure (1 atm limit) must be able to
        // return to a rate-controlled value in a few iterations
        const double pref = wells_[w].cell.empty() ? 0.0 : p_[size_t(wells_[w].cell.front())];
        const double lim = T.dpMaxRel * std::max({wells_[w].bhp, pref, 10.0 * U_.atm});
        d = std::max(-lim, std::min(lim, d));
        const double b0 = wells_[w].bhp;
        wells_[w].bhp = std::max(wells_[w].bhp + d, 0.05 * U_.atm);
        // C29 well-solution repair: a rate-controlled well whose BHP left the flowing window (no connection flows: zero rate
        // and zero derivative -> the next update explodes and the iteration flip-flops, seen with gas injected into a tiny
        // EDFM fracture cell) is put back just inside it
        WellRT& W = wells_[w];
        if (W.mode == Control::BHP || W.cell.empty()) continue;
        const bool prod = W.type == WellType::Producer;
        // flowing window: an injector flows while its BHP exceeds the lowest connection pressure (minus head), a producer
        // while its BHP is below the highest one
        double edge = prod ? -1e300 : 1e300;
        for (size_t k = 0; k < W.cell.size(); ++k) {
            if (k < W.wi.size() && W.wi[k] <= 0) continue;
            const double head = W.pConn(k, U_.gravity) - W.bhp;
            const double pc = p_[size_t(W.cell[k])] - head;
            edge = prod ? std::max(edge, pc) : std::min(edge, pc);
        }
        if (std::fabs(edge) > 1e299) continue;
        const double margin = 0.002 * std::max(std::fabs(edge), 10.0 * U_.atm);
        // back inside: half-way between the edge and the previous (flowing) BHP, or a small margin if that was outside too
        if (!prod && W.bhp <= edge) W.bhp = b0 > edge ? edge + 0.5 * (b0 - edge) : edge + margin;
        else if (prod && W.bhp >= edge) W.bhp = std::max(b0 < edge ? edge - 0.5 * (edge - b0) : edge - margin, 0.05 * U_.atm);
    }
}

bool Simulator::newton(double dt, int& its, int& lin) {
    const Tuning& T = case_.tuning;
    std::vector<double> rhs(J_.n), dx;
    lin = 0;
    std::vector<double> prevWellR(wells_.size(), 0.0);
    for (auto& w : wells_) w.stallIts = 0;
    const bool dbg = std::getenv("RESSIM_DEBUG") != nullptr;
    int nEqAct = 0; for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) ++nEqAct;
    OscillationDetector osc(nEqAct, 0.2, 0.5, 0.1, std::min(2, nEqAct));
    auto now = [] { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
    for (its = 0; its <= T.maxNewton; ++its) {
        const double tA0 = now();
        if (its < 4) updateThp();      // THP limits from the VFP tables at the latest rates (frozen later for convergence)
        if (daOn_) daNeedCells();      // C28: completions may change between iterations (well shut-in / reopen)
        evalAll(daOn_);
        assemble(dt);
        stats_.secAssemble += now() - tA0;
        if (dbg) std::fprintf(stderr, "    [assemble %.2fs]", now() - tA0);
        double cnv, mb;
        const bool conv = converged(dt, cnv, mb);
        if (case_.tuning.relaxOscillations) {
            double rr[3]; int k = 0;
            for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) rr[k++] = cnvPh_[p];
            osc.push(rr);
        }
        {   // a rate target the well cannot deliver (its residual does not move) -> BHP control at the limit
            const int na = grid_.numActive();
            for (int w = 0; w < int(wells_.size()); ++w) {
                WellRT& W = wells_[size_t(w)];
                const double r = std::fabs(R_[size_t(na) * nv_ + w]);
                // stalled = the residual did not drop by 5 % for two consecutive iterations (after a few iterations)
                W.stallIts = (its >= 2 && W.open && W.mode != Control::BHP && r > 0.5 && r > 0.95 * prevWellR[size_t(w)]) ? W.stallIts + 1 : 0;
                // a well producing / injecting MORE than its target can deliver it (slow convergence, not a deliverability
                // limit): it is never locked at the BHP limit
                const double tgt = W.open && W.mode != Control::BHP ? targetOf(W.ctrl, W.mode) : 0.0;
                const bool deliverable = tgt > 0 && wellRateOf(W, W.mode) > tgt;
                if (its >= 5 && W.stallIts >= 2 && !deliverable) {
                    const bool prod = W.type == WellType::Producer;
                    W.mode = Control::BHP; W.bhp = W.ctrl.bhp > 0 ? W.ctrl.bhp : (prod ? U_.atm : 1e12); W.ctrl.bhp = W.bhp;
                    W.lockBhp = true;
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
            if (std::getenv("RESSIM_DEBUG_WELLS"))   // C47: per-well state (open, mode, target, rate, stall counter, BHP)
                for (int w = 0; w < int(wells_.size()); ++w) {
                    const WellRT& W = wells_[size_t(w)];
                    const double tg = W.open && W.mode != Control::BHP ? targetOf(W.ctrl, W.mode) : 0.0;
                    int nOpen = 0; for (double v : W.wi) nOpen += v > 0 ? 1 : 0;
                    std::fprintf(stderr, "      well %s open %d mode %d tgt %.6g rate %.6g stall %d bhp %.4g conns %d/%d\n",
                                 case_.wells[size_t(W.input)].name.c_str(), int(W.open), int(W.mode), tg, W.open ? wellRateOf(W, W.mode) : 0.0,
                                 W.stallIts, W.bhp, nOpen, int(W.wi.size()));
                    for (size_t k = 0; k < W.cell.size() && k < 4; ++k)
                        std::fprintf(stderr, "        conn %zu cell %d z %.6g (ref %.6g) p %.5g dp %.5g pconn %.5g wi %.4g\n", k, W.cell[k], W.cdepth[k], W.refDepth,
                                     p_[size_t(W.cell[k])], k < W.dpConn.size() ? W.dpConn[k] : -1.0, W.pConn(k, U_.gravity), k < W.wi.size() ? W.wi[k] : 0.0);
                }
        }
        // C38: at least one Newton update per time step (OPM Flow NewtonMinIterations): a near-steady state whose residual is
        // already below the tolerances at iteration 0 would otherwise keep the previous well rates / pressures frozen
        // (opm-tests UDT-1D: the rates of wells at their BHP limit stopped declining). Opt out: RESSIM_NEWTON_MIN0=1
        static const bool min0 = std::getenv("RESSIM_NEWTON_MIN0") != nullptr;
        if (conv && (its >= 1 || min0)) {
            if (daOn_) evalAll(false);   // C28: full properties for the step bookkeeping (storeOld, in-place, rates, output)
            return true;
        }
        if (its == T.maxNewton) break;
        for (int i = 0; i < J_.n; ++i) rhs[i] = -R_[i];
        if (std::getenv("RESSIM_DEBUG") && its == 0 && !blkAsm_) {   // (block assembly: no scalar cell rows - RESSIM_SCALAR_JAC=1)
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
        if (const char* dp = std::getenv("RESSIM_DUMP"); dp && its == 0 && !dumped_) {   // Jacobian + rhs for offline solver work
            dumped_ = true;
            if (FILE* f = std::fopen(dp, "wb")) {
                const int na = grid_.numActive(), nw = int(wells_.size());
                std::fwrite(&J_.n, 4, 1, f); std::fwrite(&na, 4, 1, f); std::fwrite(&nv_, 4, 1, f); std::fwrite(&nw, 4, 1, f);
                const int nnz = J_.rowPtr[J_.n];
                std::fwrite(&nnz, 4, 1, f);
                std::fwrite(J_.rowPtr.data(), 4, size_t(J_.n) + 1, f); std::fwrite(J_.col.data(), 4, size_t(nnz), f);
                std::fwrite(J_.val.data(), 8, size_t(nnz), f); std::fwrite(rhs.data(), 8, size_t(J_.n), f);
                std::fclose(f);
            }
        }
        auto dumpFail = [&] {
            if (const char* dp = std::getenv("RESSIM_DUMP_FAIL"); dp && !dumpedFail_) {
                    dumpedFail_ = true;
                    if (FILE* f = std::fopen(dp, "wb")) {
                        const int na = grid_.numActive(), nw = int(wells_.size()), nnz = J_.rowPtr[J_.n];
                        std::fwrite(&J_.n, 4, 1, f); std::fwrite(&na, 4, 1, f); std::fwrite(&nv_, 4, 1, f); std::fwrite(&nw, 4, 1, f); std::fwrite(&nnz, 4, 1, f);
                        std::fwrite(J_.rowPtr.data(), 4, size_t(J_.n) + 1, f); std::fwrite(J_.col.data(), 4, size_t(nnz), f);
                        std::fwrite(J_.val.data(), 8, size_t(nnz), f); std::fwrite(rhs.data(), 8, size_t(J_.n), f);
                        std::fclose(f);
                    }
                }
        };
        const double tP0 = now();
        double tSol0 = tP0;
        // ---- primary path: Schur-eliminate the well BHPs, solve the pure block system with block CPR (AMG + BILU0)
        // or BILU0. Falls back to the scalar CSR ladder below if the elimination is not possible.
        const int bicgMax = std::max(20, std::min(T.maxLinear, 120));
        bool blockPath = case_.tuning.blockSolver &&
                         (blkAsm_ ? schur_.buildInPlace(J_, int(wells_.size()), blkStart_, blkCells_)
                                  : schur_.build(J_, grid_.numActive(), nv_, int(wells_.size())));
        // C45 block assembly: the scalar ladder below needs the whole system - built from the (unreduced) blocks for
        // this iteration only, the well-part J_ is restored at the end of the iteration
        CSR jwKeep;
        struct JRestore { CSR& a; CSR& b; bool on = false; ~JRestore() { if (on) std::swap(a, b); } } jRestore{J_, jwKeep};
        if (!blockPath && blkAsm_) { jwKeep = fullJacobian(); std::swap(J_, jwKeep); jRestore.on = true; }
        if (blockPath) {
            const BSR& S = schur_.reduced();
            bc_.assign(size_t(S.size()), 0.0);
            schur_.reduceRhs(rhs.data(), bc_.data());
            const Preconditioner* MB = &bilu_;
            // C21: when the MPI DD or the device solver will do the solve, the CPU preconditioner is built lazily (only if
            // that solve fails) - it was previously always set up and thrown away (~20% of the run on a laptop)
            const bool offload = daOn_ || distributedSolver() != nullptr || (case_.tuning.gpu && gpuBlockSolver() != nullptr);
            const bool allowReuse = !offload;          // a lazily built AMG must not reuse a hierarchy from an older matrix
            bool mbReady = false;
            auto setupBlock = [&](bool trueImpes) {
                if (useCPR_) {
            const double* acc = nullptr;
                            if (trueImpes) {
                        const int na = grid_.numActive(), nv = nv_;
                        accD_.assign(size_t(na) * nv * nv, 0.0);
                        const bool ph[3] = {W_, O_, G_};
                        // same per-cell row scaling as the assembled system (dt / PV; the constant dt does not matter)
                        #pragma omp parallel for schedule(static) if (useParallel(na > 4000))
                        for (int a = 0; a < na; ++a) {
                            const double sc = 1.0 / grid_.poreVolRef[a];
                            for (int p = 0; p < 3; ++p) if (ph[p])
                                for (int k = 0; k < nv; ++k) accD_[(size_t(a) * nv + eq_[p]) * nv + k] = cf_[a].acc[p].d[k] * sc;
                        }
                        acc = accD_.data();
                        if (const char* dp = std::getenv("RESSIM_DUMP"); dp && its == 0)
                                if (FILE* f = std::fopen((std::string(dp) + ".acc").c_str(), "wb")) { std::fwrite(accD_.data(), 8, accD_.size(), f); std::fclose(f); }
                    }
                    bcpr_.setMixed(case_.tuning.mixedPrecision);
                    bcpr_.setLean(leanLinear());
                    // C26: across time steps too (Tuning amgReuse > 0) while the hierarchy is young and still effective
                    const bool crossStep = its == 0 && case_.tuning.amgReuse > 0 && !amgStale_ && amgAge_ < case_.tuning.amgReuse;
                    bcpr_.setup(S, acc, (its > 0 || crossStep) && allowReuse, acc ? &accD_ : nullptr);
                    if (its == 0) { if (bcpr_.amgReused()) { ++amgAge_; ++amgCrossReuse_; } else { amgAge_ = 0; amgBaseIts_ = -1; amgStale_ = false; } }
                    MB = &bcpr_; stats_.amgLevels = bcpr_.amg().levels();
                } else if (leanLinear()) { bilu_.setLean(true); bilu_.factor(S); MB = &bilu_; }   // C45 low-memory BILU0
                else if (case_.tuning.mixedPrecision) {   // C38: level-parallel BILU0 with float32 factors where possible
                    bilu_.setFloatFactors(true); bilu_.factor(S);
                    if (bilu_.contiguousLevels()) MB = &bilu_; else { biluf_.factor(S); MB = &biluf_; }
                }
                else { bilu_.setLean(false); bilu_.setFloatFactors(false); bilu_.factor(S); MB = &bilu_; }
            };
            auto ensureMB = [&] {
                if (mbReady) return;
                const double t0 = now();
                setupBlock(case_.tuning.cprWeights == 1 && !quasiFallback_);
                mbReady = true; stats_.secPrecond += now() - t0;
            };
            stats_.secPrecond += now() - tP0;          // Schur reduction
            if (!offload) ensureMB();
            const double tS0 = now();
            tSol0 = tS0;
            if (dbg) std::fprintf(stderr, "    [block precond setup %.2fs, Schur fill %d blocks, AMG %s]", tS0 - tP0, schur_.fillBlocks(),
                                  useCPR_ && bcpr_.amgReused() ? "reused" : "new");
            BsrOperator op(S);
            std::vector<double> xc;
            LinResult lr;
            bool deviceDone = false;
            if (daOn_) {
                // C28 distributed assembly: only the owned rows exist on this rank - the collective solve is the only
                // possible path; a failure chops the step on every rank alike (the outcome is a global quantity)
                int itsD = 0;
                const double tD0 = now();
                const bool okD = da_.solve(S, bc_, xc, T.linTol, bicgMax, itsD);
                stats_.ddSec += now() - tD0; stats_.ddSolves++;
                lin += itsD;
                if (!okD) { stats_.ddFallbacks++; log_("  distributed linear solve failed - time step cut"); return false; }
                deviceDone = true; lr.converged = true; lr.iterations = itsD; lr.relResidual = T.linTol;
            }
            if (!deviceDone && distributedSolver()) {
                // C18: MPI domain-decomposed solve (double precision, all ranks); CPU path on failure
                int itsD = 0;
                const double tD0 = now();
                const bool okD = distributedSolver()->solve(S, bc_, xc, T.linTol, bicgMax, itsD);
                stats_.ddSec += now() - tD0; stats_.ddSolves++;
                lin += itsD;
                if (okD) { deviceDone = true; lr.converged = true; lr.iterations = itsD; lr.relResidual = T.linTol; }
                else { stats_.ddFallbacks++; xc.clear(); }
            }
            if (!deviceDone && case_.tuning.gpu && gpuBlockSolver()) {
                // C16: float32 device Krylov solve inside double iterative refinement; CPU path on failure
                const double dev0 = gpuBlockSolver()->deviceSeconds(), tg0 = now();
                lr = refineWithDevice(*gpuBlockSolver(), S, bc_, xc, T.linTol, 6, bicgMax);
                stats_.gpuWallSec += now() - tg0;
                const double dev1 = gpuBlockSolver()->deviceSeconds();
                if (dev0 >= 0 && dev1 >= dev0) stats_.gpuDeviceSec = std::max(0.0, stats_.gpuDeviceSec) + (dev1 - dev0);
                lin += lr.iterations; stats_.gpuSolves++;
                if (lr.converged) deviceDone = true;
                else { stats_.gpuFallbacks++; xc.clear(); }
            }
            if (!deviceDone) {
                ensureMB(); lr = bicgstab(op, *MB, bc_, xc, T.linTol, bicgMax);
                if (useCPR_ && case_.tuning.amgReuse > 0) {   // C26 quality monitor of the reused hierarchy
                    if (amgBaseIts_ < 0 && !bcpr_.amgReused()) amgBaseIts_ = lr.iterations;
                    else if (!lr.converged || (amgBaseIts_ >= 0 && lr.iterations > 1.5 * amgBaseIts_ + 3)) amgStale_ = true;
                }
            }
            if (dbg) std::fprintf(stderr, " [bicgstab %d its, %.2fs, res %.1e]\n", lr.iterations, now() - tS0, lr.relResidual);
            if (!deviceDone) lin += lr.iterations;
            if (!lr.converged) {
                if (useCPR_ && case_.tuning.cprWeights == 1 && !quasiFallback_) {
                    quasiFallback_ = true;
                    log_("  CPR: true-IMPES decoupling did not converge - using quasi-IMPES weights for the rest of the run");
                    setupBlock(false);
                }
                std::vector<double> xg;
                LinResult lg = gmres(op, *MB, bc_, xg, T.linTol, T.maxLinear, 60);
                lin += lg.iterations;
                if (dbg) std::fprintf(stderr, "      [gmres %d its, res %.1e]\n", lg.iterations, lg.relResidual);
                if (lg.relResidual < lr.relResidual || xc.empty()) { xc.swap(xg); lr = lg; }
                if (!lr.converged && lr.relResidual > 1e-3 && (useCPR_ || case_.tuning.mixedPrecision || leanLinear())) {   // double BILU0 rescue
                    bilu_.setLean(false); bilu_.factor(S);
                    std::vector<double> xi;
                    LinResult li = gmres(op, bilu_, bc_, xi, T.linTol, T.maxLinear, 60);
                    lin += li.iterations;
                    if (li.relResidual < lr.relResidual) { xc.swap(xi); lr = li; }
                }
                if (!lr.converged && lr.relResidual > 1e-2) {
                    dumpFail();
                    log_(fmt("  linear solver failed (res %.2e)", lr.relResidual)); return false;
                }
            }
            dx.assign(size_t(J_.n), 0.0);
            schur_.recover(rhs.data(), xc.data(), dx.data());
        }
        if (!blockPath && daOn_) { log_("  distributed assembly needs the block (Schur) solver path - time step cut"); return false; }
        if (!blockPath) {
            const Preconditioner* M = &ilu_;
            if (useCPR_) {
                const double* acc = nullptr;
                if (case_.tuning.cprWeights == 1 && !quasiFallback_) {
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
                    if (const char* dp = std::getenv("RESSIM_DUMP"); dp && its == 0)
                        if (FILE* f = std::fopen((std::string(dp) + ".acc").c_str(), "wb")) { std::fwrite(accD_.data(), 8, accD_.size(), f); std::fclose(f); }
                }
                cpr_.setup(J_, grid_.numActive(), nv_, int(wells_.size()), acc);
                M = &cpr_; stats_.amgLevels = cpr_.amg().levels();
            }
            else ilu_.factor(J_);
            const double tS0 = now();
            stats_.secPrecond += tS0 - tP0; tSol0 = tS0;
            if (dbg) std::fprintf(stderr, "    [precond setup %.2fs]", tS0 - tP0);
            // Krylov with a robustness ladder. BiCGSTAB stagnates long before maxLinear when the preconditioner is
            // poor, so its budget is capped; then quasi-IMPES CPR + GMRES (true-IMPES decoupling can fail in gas caps /
            // phase-absent cells at large time steps, and is then switched off for the rest of the run); then ILU0 GMRES.
            LinResult lr = bicgstab(J_, *M, rhs, dx, T.linTol, bicgMax);
            if (dbg) std::fprintf(stderr, " [bicgstab %d its, %.2fs, res %.1e]\n", lr.iterations, now() - tS0, lr.relResidual);
            lin += lr.iterations;
            if (!lr.converged) {
                if (useCPR_ && case_.tuning.cprWeights == 1 && !quasiFallback_) {
                    quasiFallback_ = true;
                    log_("  CPR: true-IMPES decoupling did not converge - using quasi-IMPES weights for the rest of the run");
                    cpr_.setup(J_, grid_.numActive(), nv_, int(wells_.size()), nullptr);
                    M = &cpr_;
                }
                LinResult lg = gmres(J_, *M, rhs, dx, T.linTol, T.maxLinear, 60);
                lin += lg.iterations;
                if (dbg) std::fprintf(stderr, "      [gmres %d its, res %.1e]\n", lg.iterations, lg.relResidual);
                if (!lg.converged && lg.relResidual > 1e-3) {
                    ilu_.factor(J_);
                    std::vector<double> dx2;
                    LinResult li = gmres(J_, ilu_, rhs, dx2, T.linTol, T.maxLinear, 60);
                    lin += li.iterations;
                    if (li.relResidual < lg.relResidual) { dx.swap(dx2); lg = li; }
                }
                if (!lg.converged && lg.relResidual > 1e-2) {
                    dumpFail();
                    log_(fmt("  linear solver failed (res %.2e)", lg.relResidual)); return false;
                }
            }
        }
        const double tU0 = now();
        stats_.secSolve += tU0 - tSol0;   // Krylov iterations (incl. fallbacks)
        for (double v : dx) if (!std::isfinite(v)) { log_("  non-finite Newton update"); return false; }
        if (osc.omega() < 1.0) {           // oscillating / stagnating residuals: relaxed (damped) Newton update
            for (double& v : dx) v *= osc.omega();
            ++stats_.relaxedIts;
            if (dbg) std::fprintf(stderr, "      [relaxation omega %.2f%s]\n", osc.omega(), osc.oscillating() ? " (oscillation)" : " (stagnation)");
        }
        updateState(dx);
        stats_.secUpdate += now() - tU0;
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
        if (headsOnly_ >= 0 && &w != &wells_[size_t(headsOnly_)]) continue;
        if (headsMswOnly_ && mswOfWell(int(&w - wells_.data())) < 0) continue;
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
        if (mswOfWell(int(&w - wells_.data())) >= 0) {   // C32: device conductance factors in force (restored after a report-step reset / restart)
            const size_t ns = msw_[size_t(mswOfWell(int(&w - wells_.data())))].segs.size();
            if (w.dpF.size() >= ns + n)
                for (size_t k = 0; k < n; ++k) if (w.wi[k] > 0 && k < w.wiBase.size()) w.wi[k] = w.wiBase[k] * w.dpF[ns + k];
        }
        std::vector<double> pOld(n);
        for (size_t k = 0; k < n; ++k) pOld[k] = w.pConn(k, g);
        // 1. connection component rates (oil, water, gas surface volumes; + = out of the reservoir)
        std::vector<double> qc(3 * n, 0.0);
        std::vector<double> rhoCell(n, 0.0), qresC(n, 0.0), pCellC(n, 0.0), aResC(n, 0.0);   // C29: connection reservoir rate, cell pressure, reservoir PI
        const bool isMsw = mswOfWell(int(&w - wells_.data())) >= 0;   // C32
        std::vector<double> muC(isMsw ? n : 0, 1.0);
        std::vector<double> fgC(isMsw ? n : 0, 0.0), rgC(isMsw ? n : 0, 0.0), rlC(isMsw ? n : 0, 0.0);   // C36 drift flux inputs
        const bool xfH = isMsw && prod && w.xflow && mswXflow();   // C45 MSW cross-flow volume for the tubing storage sweep
        w.xfRes.assign(xfH ? n : 0, 0.0);
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
            pCellC[k] = C.P[refP].v;
            if (isMsw) {   // C36: free-gas fraction of the reservoir inflow, gas / liquid densities
                double ml = 0, rl = 0;
                for (int p : {PO, PW}) if (eq_[p] >= 0) { const double m = C.mob[p].v / std::max(C.b[p].v, 1e-30); ml += m; rl += m * C.rho[p].v; }
                rlC[k] = ml > 0 ? rl / ml : (O_ ? C.rho[PO].v : C.rho[PW].v);
                if (G_ && eq_[PG] >= 0) { rgC[k] = C.rho[PG].v; fgC[k] = den > 0 ? C.mob[PG].v / std::max(C.b[PG].v, 1e-30) / den : 0.0; }
            }
            if (w.wi[k] <= 0) continue;
            if (isMsw) {   // C32: inflow mixture viscosity (reservoir-rate weighted), cP
                const FluidModel& F = fl(w.cell[k]);
                double sw = 0, sm = 0;
                for (int p = 0; p < 3; ++p) {
                    if (eq_[p] < 0) continue;
                    const double wt = C.mob[p].v / std::max(C.b[p].v, 1e-30);
                    double mup = 1.0;
                    const AD3 P(C.P[refP].v);
                    if (p == PW) mup = F.bw(P).v / std::max(F.bwMu(P).v, 1e-30);
                    else if (p == PO) mup = F.bo(P, AD3(C.rs.v)).v / std::max(F.boMu(P, AD3(C.rs.v)).v, 1e-30);
                    else mup = vap_ ? F.bg(P, AD3(C.rv.v)).v / std::max(F.bgMu(P, AD3(C.rv.v)).v, 1e-30) : F.bg(P).v / std::max(F.bgMu(P).v, 1e-30);
                    sw += wt; sm += wt * mup;
                }
                muC[k] = sw > 0 ? sm / sw : 1.0;
            }
            if (!case_.wellbores.empty() || isMsw) {
                const double wm0 = (k < w.wmulConn.size()) ? w.wmulConn[k] : 1.0;
                const double wiA = isMsw && k < w.wiBase.size() ? w.wiBase[k] : w.wi[k];   // C32: device-free reservoir PI
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) aResC[k] += C.mob[p].v / std::max(C.b[p].v, 1e-30) * wiA * (p == PW && prod ? wm0 : 1.0);
            }
            const double dd = C.P[refP].v - pOld[k];
            if (prod && dd > 0) {
                const double wm = (k < w.wmulConn.size()) ? w.wmulConn[k] : 1.0;
                if (eq_[PO] >= 0) qc[3 * k + 0] = C.mob[PO].v * dd * w.wi[k];
                if (eq_[PW] >= 0) qc[3 * k + 1] = C.mob[PW].v * dd * w.wi[k] * wm;
                if (eq_[PG] >= 0) qc[3 * k + 2] = C.mob[PG].v * dd * w.wi[k];
                if (live_) qc[3 * k + 2] += C.mobRs.v * dd * w.wi[k];
                if (vap_) qc[3 * k + 0] += C.mobRv.v * dd * w.wi[k];
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) qresC[k] += C.mob[p].v / std::max(C.b[p].v, 1e-30) * dd * w.wi[k] * (p == PW ? wm : 1.0);
            } else if (prod && dd < 0 && xfH) {
                // C45: cross-flowing producer connection of a multisegment well (total-mobility injection of the tubing fluid)
                double lt = 0;
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) lt += C.mob[p].v / std::max(C.b[p].v, 1e-30);
                w.xfRes[k] = -lt * dd * w.wi[k];
            } else if (!prod && dd < 0) {
                const int ip = w.type == WellType::GasInjector ? PG : PW;
                double lt = 0;
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) lt += C.mob[p].v / std::max(C.b[p].v, 1e-30);
                if (solv_ && C.bs.v > 0) lt += C.mobS.v / C.bs.v;   // C35 solvent
                if (eq_[ip] >= 0) qc[3 * k + (ip == PG ? 2 : 1)] = lt * C.b[ip].v * dd * w.wi[k];   // < 0: injection
                qresC[k] = -lt * dd * w.wi[k];
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
        // 4. C29 WELLBORE: lagged friction + ICD drops (FlexWell-lite)
        if (isMsw) mswHeads(w, qresC, rhoCell, muC, pCellC, aResC, &fgC, &rgC, &rlC, &qc);   // C32 segment network replaces the standard heads
        else if (!case_.wellbores.empty()) wellboreFriction(w, qresC, rhoCell, pCellC, aResC);
    }
}

// C29 discretised-wellbore friction and inflow-control devices (Case::wellbores). The tubing is the polyline datum ->
// connection 0 -> 1 -> ... in completion (COMPDAT) order (heel to toe / top to bottom); the segment between two connections
// carries the reservoir rates of all connections beyond it (toe side) at their mobility-weighted in-situ density. Darcy
// friction factor: 64 / Re (Re < 2300), Haaland (1983) otherwise, mixture viscosity 1 cP. Drop per segment
// f rho v |v| L / (2 D). ICD per connection: K rho_mix / rho_cal q^2. Producers: the connection pressure rises by the
// accumulated friction + its ICD drop (more drawdown for the same BHP); injectors: falls. Under-relaxed 1/2 with the
// previous step (explicit lag, as the hydrostatic heads).
void Simulator::wellboreFriction(WellRT& w, const std::vector<double>& qres, const std::vector<double>& rhoMix, const std::vector<double>& pCell,
                                 const std::vector<double>& aRes) {
    const std::string& nm = case_.wells[size_t(w.input)].name;
    const Case::WellboreSpec* B = nullptr;
    for (const auto& b : case_.wellbores) if (b.well == nm) B = &b;
    const size_t n = w.cell.size();
    if (!B || n == 0) { w.dpFric.clear(); w.dpF.clear(); return; }
    const bool field = case_.units == UnitSystem::Field;
    const double toM = field ? 0.3048 : 1.0, toRho = field ? 16.01846337 : 1.0;
    const double toQ = (field ? 0.158987294928 : 1.0) / 86400.0;            // reservoir rate -> m3/s
    const double toPa = 101325.0 / U_.atm;
    const double D = std::max(1e-3, B->id * toM), eps = std::max(0.0, B->rough * toM), A = 0.25 * 3.14159265358979323846 * D * D;
    const double mu = 1e-3;
    const bool prod = w.type == WellType::Producer;
    std::vector<Vec3> pos(n);
    for (size_t k = 0; k < n; ++k) pos[k] = grid_.center[size_t(grid_.globalOf[size_t(w.cell[k])])];
    // toe-side cumulative rates and densities
    std::vector<double> qSeg(n, 0.0), rSeg(n, 0.0);
    double qs = 0, rq = 0;
    for (size_t k = n; k-- > 0;) {
        const double q = std::fabs(qres[k]);
        qs += q; rq += q * rhoMix[k];
        qSeg[k] = qs; rSeg[k] = qs > 0 ? rq / qs : rhoMix[k];
    }
    std::vector<double> Fk(n, 0.0);
    double F = 0;
    for (size_t k = 0; k < n; ++k) {
        double L;
        if (k == 0) L = std::fabs(w.cdepth[0] - w.refDepth);
        else { const double dx = pos[k].x - pos[k - 1].x, dy = pos[k].y - pos[k - 1].y, dz = w.cdepth[k] - w.cdepth[k - 1]; L = std::sqrt(dx * dx + dy * dy + dz * dz); }
        const double v = qSeg[k] * toQ / A, rho = rSeg[k] * toRho;
        if (v > 0 && L > 0) {
            const double Re = rho * v * D / mu;
            double f;
            if (Re < 2300) f = 64.0 / std::max(Re, 1e-12);
            else { const double t = -1.8 * std::log10(std::pow(eps / D / 3.7, 1.11) + 6.9 / Re); f = 1.0 / (t * t); }
            F += f * rho * v * v * L * toM / (2 * D) / toPa;
        }
        Fk[k] = F;
    }
    // friction: explicit, under-relaxed 1/2 with the previous step
    if (w.dpF.size() == n) for (size_t k = 0; k < n; ++k) w.dpF[k] = 0.5 * (w.dpF[k] + Fk[k]);
    else w.dpF = Fk;
    // ICD: local balance dp = c q^2 with q = a (D0 - dp) (reservoir PI a, drawdown D0 without the ICD at the lagged BHP):
    // u = D0 - dp solves c a^2 u^2 + u - D0 = 0 - exact for the connection, no explicit-feedback instability
    w.dpFric.assign(n, 0.0);
    const double sg = prod ? 1.0 : -1.0;
    for (size_t k = 0; k < n; ++k) {
        double icd = 0;
        if (B->icdK > 0 && k < aRes.size() && aRes[k] > 0 && k < w.dpConn.size()) {
            const double rc = B->rhoCal > 0 ? B->rhoCal : fl(w.cell[k]).rhoSC().water;
            const double c = B->icdK * rhoMix[k] / std::max(rc, 1e-12);
            const double pc = w.bhp + w.dpConn[k] + sg * w.dpF[k];
            const double D0 = prod ? pCell[k] - pc : pc - pCell[k];
            if (D0 > 0 && c > 0) {
                const double ca2 = c * aRes[k] * aRes[k];
                const double u = 2 * D0 / (1 + std::sqrt(1 + 4 * ca2 * D0));   // stable root form
                icd = D0 - u;
            }
        }
        w.dpFric[k] = sg * (w.dpF[k] + icd);
    }
    for (size_t k = 0; k < n && k < w.dpConn.size(); ++k) w.dpConn[k] += w.dpFric[k];
}

void Simulator::storeOld() {
    const int na = grid_.numActive();
    // C35: the accepted state with the final solvent fraction (the last transport follows the flow solve); the solvent
    // amount of the transport is kept exactly (old-time solvent mass), the flow accumulations use the final split
    std::vector<double> ntKeep;
    if (solv_) { ntKeep.resize(size_t(na)); for (int a = 0; a < na; ++a) ntKeep[size_t(a)] = zs_[size_t(a)] * solventTotalNT(a); evalAll(); }
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
    solventAccept();
    if (!ntKeep.empty()) ntOld_ = ntKeep;
    chemAccept();
    updateHysteresis();
    updateVappars();
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
        double rs = st == 0 ? x_[a] : rsSatC(a, p_[a]);
        if (st == 1 && !rsLim_.empty()) rs = std::min(rs, rsLim_[size_t(a)]);
        lastRs_[size_t(a)] = rs;
    }
}

void Simulator::drsdtControls(int step) {
    const ReportStep& S = case_.schedule[size_t(step)];
    // DRSDT and VAPPARS exclude each other (OPM OilVaporizationProperties / ECLIPSE): the later keyword wins
    auto applyVap = [&] {
        if (!(O_ && G_)) return;
        vap1_ = S.vap1; vap2_ = std::max(0.0, S.vap2);
        log_(fmt("  VAPPARS %.4g %.4g", vap1_, vap2_) + (drsdt_ >= 0 ? " (switches DRSDT off)" : ""));
        drsdt_ = -1; drsdtFree_ = false; rsLim_.clear(); lastRs_.clear();
        updateVappars();
    };
    auto applyDrsdt = [&] {
        drsdt_ = S.drsdt; drsdtFree_ = S.drsdtFree;
        if (drsdt_ >= 0) log_(fmt("  DRSDT %.4g /day", drsdt_) + (drsdtFree_ ? " (FREE)" : ""));
        if (vap1_ > 0 || vap2_ > 0) { log_("  DRSDT switches VAPPARS off"); vap1_ = vap2_ = 0; updateVappars(); }
        drsdtAccept();
        if (drsdt_ < 0) rsLim_.clear();
    };
    const bool hv = S.vap1 >= 0, hd = S.drsdt > -1.5;
    if (hv && hd && S.vapAfterDrsdt) { applyDrsdt(); applyVap(); }
    else { if (hv) applyVap(); if (hd) applyDrsdt(); }
}

bool Simulator::runReportStep() {
    if (finished()) return false;
    auto t0 = std::chrono::steady_clock::now();
    applyControls(step_);
    // C26: the CPR-AMG hierarchy is reused across the time steps of one report step only - a fresh build at every report
    // step keeps restarts (and iterative-coupling rewinds) from a report-step state bitwise identical to the continuous run
    amgAge_ = std::max(amgAge_, case_.tuning.amgReuse);
    for (auto& w : wells_) { w.lockBhp = false; w.stallIts = 0; }   // new report step: rate targets are retried
    polymerControls(step_);
    solventControls(step_);
    chemControls(step_);
    drsdtControls(step_);
    groupControls(step_);
    gasLiftControls(step_);                   // C35 LIFTOPT / WLIFTOPT / GLIFTOPT
    // C31: ASSIGN records due at this report step (DEFINEs keep the values of the last accepted step - a restart restores
    // them from the state); the very first report step evaluates everything so that UDAs have values for its first step
    if (udq_) udqEvaluate(nullptr, nullptr, !(step_ == 0 && summary_.empty()));
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
        wcycleStep(dt, clipped);               // C37 WCYCLE: switches, start-up ramp, step landing on the switch times
        headDt_ = dt;                          // C37 segment storage
        // backup
        const auto p0 = p_, s0 = sw_, x0 = x_; const auto h0 = hasGas_;
        const auto pc0 = pc_, wm0 = wmul_;
        chemBackup();
        if (solv_) zsBak_ = zs_;
        const auto tip0 = hfTip_;                  // C21 sub-grid fracture tips at the start of the step
        for (size_t w = 0; w < wells_.size(); ++w) wb[w] = {wells_[w].bhp, wells_[w].mode};
        if (!lastRs_.empty()) {   // DRSDT: Rs may grow by at most DRSDT * dt in this time step
            rsLim_.resize(lastRs_.size());
            for (size_t a = 0; a < lastRs_.size(); ++a) rsLim_[a] = lastRs_[a] + drsdt_ * dt;
        } else rsLim_.clear();
        // C43 WPWE: the initial configuration is no event; C44: per well - a well that enters the schedule later (its
        // first connections appear) starts its baseline then, so its arrival is not counted as opened connections
        for (auto& w : wells_) if (!w.pweBase && pweExists(w)) pweBaseline(w);
        econTest();
        updateResvConversion();
        applyUda();                           // C31 UDA targets from the latest UDQ values
        if (gasLiftOptimize()) updateThp();   // C35 lift gas allocation -> THP-limit BHPs before the group targets (potentials)
        applyGroupTargets();
        // C29 surface network: node pressures -> well THP / BHP limits; C36 NETBALAN interval > 0: balanced at the first time
        // step of each interval, the node pressures held in between; NETCOMPA compressors switched before the balance
        const bool netDue = netBal_.interval <= 0 || time_ - netLastBal_ >= netBal_.interval * (1 - 1e-9);
        if (netDue && !netNodes_.empty()) { networkCompressors(); networkBalance(); }
        if (adjRec_) {   // adjoint: inputs of the explicit well-head update (old BHP and connection heads)
            adjHeadIn_.resize(wells_.size());
            for (size_t w = 0; w < wells_.size(); ++w) adjHeadIn_[w] = {wells_[w].bhp, wells_[w].dpConn};
        }
        computeHeads();
        int its = 0, lin = 0, totIts = 0, totLin = 0;
        bool ok = newton(dt, its, lin);
        totIts += its; totLin += lin;
        for (int sw = 0; ok && sw < 4 && checkWellSwitching(); ++sw) {
            ok = newton(dt, its, lin); totIts += its; totLin += lin;
        }
        // C49 converged multisegment heads (Picard coupling): the segment network is re-solved from the end-of-step state and,
        // while any connection head moved by more than the tolerance, the step is re-solved from the converged state (warm
        // start, as the VREP loop). From the third pass on the new heads are averaged with the previous ones (the explicit
        // heads flip-flopped between a gas- and a liquid-filled tubing, ENGINES_AND_MATH §9.v). Not used while an adjoint
        // records the forward run (the adjoint differentiates the explicit start-of-step heads).
        if (ok && !msw_.empty() && !adjRec_ && mswConvergedHeads()) {
            // fixed point h = G(h) of the connection heads h (all multisegment connections, one vector): G = heads of the
            // end-of-step state solved with h. G is decreasing in each h_j (a heavier head -> less inflow -> the column holds less
            // of the liquid below / more of the swept fluid) and very steep where a connection changes between producing and
            // cross-flow: plain substitution and Aitken relaxation both ran into a 2-cycle (MSW-2D-VERT-02: gas- / liquid-filled
            // tubing, 26 bar). Each head is therefore bracketed separately: regula falsi with the Illinois modification on
            // r_j = G_j - h_j (as the network balance), the first trial being the substitution h = G.
            const double tol = 0.02 * U_.atm;
            const int maxPass = 8;
            bool stepConv = false; int passes = 0;
            std::vector<double> lo, hi, rlo, rhi; std::vector<int> last, hasLo, hasHi;
            for (int it = 0; ok && it <= maxPass; ++it) {
                std::vector<std::vector<double>> before(wells_.size());
                for (size_t i = 0; i < wells_.size(); ++i) if (mswOfWell(int(i)) >= 0) before[i] = wells_[i].dpConn;
                headsMswOnly_ = true; computeHeads(); headsMswOnly_ = false;
                std::vector<double> h, gv;
                for (size_t i = 0; i < wells_.size(); ++i) {
                    WellRT& w = wells_[i];
                    if (mswOfWell(int(i)) < 0 || !w.open || before[i].size() != w.dpConn.size()) continue;
                    for (size_t k = 0; k < w.dpConn.size(); ++k) { h.push_back(before[i][k]); gv.push_back(w.dpConn[k]); }
                }
                if (lo.size() != h.size()) { lo.assign(h.size(), 0); hi = lo; rlo = lo; rhi = lo; last.assign(h.size(), 0); hasLo = last; hasHi = last; }
                double ch = 0, chB = 0;   // largest residual; largest residual of heads not yet bracketed within tol
                std::vector<double> next(gv);
                for (size_t j = 0; j < h.size(); ++j) {
                    const double r = gv[j] - h[j];
                    ch = std::max(ch, std::fabs(r));
                    if (std::fabs(r) < tol) { next[j] = h[j] + r; continue; }
                    // a head whose root is bracketed within tol (G jumps there: a connection changing between producing and
                    // cross-flow) counts as located; it keeps the bracketed value
                    if (r > 0) { lo[j] = h[j]; rlo[j] = r; hasLo[j] = 1; if (last[j] == 1 && hasHi[j]) rhi[j] *= 0.5; last[j] = 1; }
                    else { hi[j] = h[j]; rhi[j] = r; hasHi[j] = 1; if (last[j] == -1 && hasLo[j]) rlo[j] *= 0.5; last[j] = -1; }
                    const bool located = hasLo[j] && hasHi[j] && std::fabs(hi[j] - lo[j]) < tol;   // h is one end of the bracket
                    if (!located) chB = std::max(chB, std::fabs(r));
                    if (hasLo[j] && hasHi[j]) next[j] = lo[j] + (hi[j] - lo[j]) * rlo[j] / std::max(rlo[j] - rhi[j], 1e-300);
                }
                if (std::getenv("RESSIM_DEBUG_MSWHEADS")) {
                    std::fprintf(stderr, "t=%.3f dt=%.3f pass %d max dhead %.4g |", time_, dt, it, ch / U_.atm);
                    for (size_t j = 0; j < h.size(); ++j) std::fprintf(stderr, " h%zu %.2f->G %.2f next %.2f", j, h[j], gv[j], next[j]);
                    std::fprintf(stderr, "\n");
                }
                if (ch < tol || chB < tol) {   // converged, or every remaining head located at a jump of G: keep the solved heads
                    stepConv = true;
                    if (ch >= tol)
                        for (size_t i = 0; i < wells_.size(); ++i)
                            if (mswOfWell(int(i)) >= 0 && wells_[i].open && before[i].size() == wells_[i].dpConn.size()) wells_[i].dpConn = before[i];
                    break;
                }
                if (it == maxPass) {   // not converged: keep the heads the accepted solution was computed with
                    for (size_t i = 0; i < wells_.size(); ++i)
                        if (mswOfWell(int(i)) >= 0 && wells_[i].open && before[i].size() == wells_[i].dpConn.size()) wells_[i].dpConn = before[i];
                    break;
                }
                for (size_t i = 0, j = 0; i < wells_.size(); ++i) {
                    WellRT& w = wells_[i];
                    if (mswOfWell(int(i)) < 0 || !w.open || before[i].size() != w.dpConn.size()) continue;
                    for (size_t k = 0; k < w.dpConn.size(); ++k, ++j) w.dpConn[k] = next[j];
                }
                ++passes;
                ok = newton(dt, its, lin); totIts += its; totLin += lin;
                for (int sw = 0; ok && sw < 2 && checkWellSwitching(); ++sw) { ok = newton(dt, its, lin); totIts += its; totLin += lin; }
            }
            mswPicardIts_ += passes; ++mswPicardSteps_;
            if (ok && !stepConv) ++mswPicardUnconv_;
        }
        // C34 NUPCOL-like coupling of VREP / REIN injection targets to the production of THIS time step (OPM updates group
        // targets inside the first NUPCOL Newton iterations; the lagged target injected nothing in the first step and
        // under-replaced voidage while production rose). Producers keep their allocation; injection targets that moved
        // by more than 1 % re-solve the step from the converged state.
        bool voidTargets = false;
        for (int ph = 0; ph < 2; ++ph) for (auto& [n, g] : ginj_[ph]) if (g.mode == "VREP" || g.mode == "REIN") voidTargets = true;
        for (int it = 0; ok && voidTargets && it < 2; ++it) {
            struct Keep { WellControl ctrl; Control mode; bool grp; };
            std::vector<Keep> keep(wells_.size());
            std::vector<double> before(wells_.size());
            for (size_t i = 0; i < wells_.size(); ++i) { keep[i] = {wells_[i].ctrl, wells_[i].mode, wells_[i].grpCtl}; before[i] = wells_[i].ctrl.rate; }
            applyGroupTargets();
            double ch = 0;
            for (size_t i = 0; i < wells_.size(); ++i) {
                WellRT& w = wells_[i];
                if (w.type == WellType::Producer) { w.ctrl = keep[i].ctrl; w.mode = keep[i].mode; w.grpCtl = keep[i].grp; continue; }
                if (!w.open || !w.grpCtl) continue;
                const double a = before[i], b = w.ctrl.rate;
                if (keep[i].mode == Control::BHP && w.mode != Control::BHP && w.lockBhp) w.mode = Control::BHP;
                ch = std::max(ch, std::fabs(b - a) / std::max({std::fabs(a), std::fabs(b), 1e-9}));
            }
            if (ch < 0.01) break;
            computeHeads();
            ok = newton(dt, its, lin); totIts += its; totLin += lin;
            for (int sw = 0; ok && sw < 2 && checkWellSwitching(); ++sw) { ok = newton(dt, its, lin); totIts += its; totLin += lin; }
        }
        if (ok && !netNodes_.empty()) {
            // C30 implicit network balance: node pressures p must equal the branch pressures G(p) computed from the rates the
            // wells deliver at p. G is decreasing in p (higher back-pressure -> less flow -> smaller branch loss), so the
            // residual r = G(p) - p is decreasing and has one root per node: regula falsi with the Illinois modification,
            // first trial = the fixed-point value. Each trial re-solves the step from its start. The explicit C29 balance
            // (lagged rates) produced a limit cycle when the loop gain PI x dG/dq approached 1 (well dies / revives).
            std::vector<double> pU, pN; std::vector<int> wn;
            networkPressures(&pU, &wn);
            for (size_t i = 0; i < pU.size(); ++i) pU[i] = (pU[i] >= 0 && !netNodes_[i].fixed) ? netNodes_[i].p : pU[i];
            const size_t nn = pU.size();
            std::vector<double> lo(nn, -1), hi(nn, -1), rlo(nn, 0), rhi(nn, 0); std::vector<int> lastSide(nn, 0);
            double worst = 0;
            const int netMax = std::getenv("RESSIM_NET_EXPLICIT") || !netDue ? 0 : netBal_.given ? netBal_.maxIter : netIterMax_;   // diagnostics: C29 explicit balance
            double worstAbs = 0; bool netConv = true;
            for (int it = 0; ok && it < netMax; ++it) {
                networkPressures(&pN, &wn);
                worst = 0; worstAbs = 0; bool conv = true;
                std::vector<double> next = pU;
                for (size_t i = 0; i < nn; ++i) {
                    if (pN[i] < 0 || netNodes_[i].fixed || pU[i] < 0) continue;
                    const double r = pN[i] - pU[i], tol = netBal_.given ? netBal_.ptol : 2e-3 * std::max(pN[i], U_.atm) + 1e-3 * U_.atm;
                    worst = std::max(worst, std::fabs(r) / std::max(pN[i], U_.atm));
                    worstAbs = std::max(worstAbs, std::fabs(r));
                    if (std::fabs(r) <= tol) continue;
                    conv = false;
                    if (r > 0) { lo[i] = pU[i]; rlo[i] = r; if (lastSide[i] == 1 && hi[i] >= 0) rhi[i] *= 0.5; lastSide[i] = 1; }
                    else { hi[i] = pU[i]; rhi[i] = r; if (lastSide[i] == -1 && lo[i] >= 0) rlo[i] *= 0.5; lastSide[i] = -1; }
                    if (lo[i] >= 0 && hi[i] >= 0) next[i] = lo[i] + (hi[i] - lo[i]) * rlo[i] / std::max(rlo[i] - rhi[i], 1e-300);
                    else next[i] = pN[i];
                }
                netConv = conv;
                if (conv) break;
                ++netIters_;
                p_ = p0; sw_ = s0; x_ = x0; hasGas_ = h0; pc_ = pc0; wmul_ = wm0; chemRestore(); if (solv_) zs_ = zsBak_;
                for (size_t w = 0; w < wells_.size(); ++w) { wells_[w].bhp = wb[w].bhp; wells_[w].mode = wb[w].mode; }
                pU = next;
                networkApply(pU, wn);
                computeHeads();
                ok = newton(dt, its, lin); totIts += its; totLin += lin;
                for (int sw = 0; ok && sw < 4 && checkWellSwitching(); ++sw) { ok = newton(dt, its, lin); totIts += its; totLin += lin; }
                ++netResolves_;
                if (it == netMax - 1 && ok && netWarn_.insert(step_).second)
                    log_(fmt("  network balance: not converged in %g iterations at t = %.2f d (largest node mismatch %.2g %%) - step accepted", double(netMax), time_ + dt, 100 * worst));
            }
            if (netDue && netMax > 0) netLastBal_ = time_;
            if (ok && netMax > 0 && !netConv && netBal_.maxErr > 0 && worstAbs > netBal_.maxErr) {   // C36 NETBALAN item 7
                ok = false;
                log_(fmt("  network balance error %.3g above NETBALAN MAX_BALANCE_ERROR %.3g - time step cut", worstAbs, netBal_.maxErr));
            }
        }
        if (ok && (poly_ || !chem_.empty() || solv_)) {   // flow <-> polymer / C35 chemical / solvent outer coupling iterations
            const double tol = 1e-4, tolC = 1e-3, tolS = 1e-4;
            for (int outer = 0;; ++outer) {
                const double chS = solv_ ? solventTransport(dt) : 0.0;
                if (chS < 0) { ok = false; log_("  solvent transport did not converge"); break; }
                if (solv_) svOuter_ = std::max(svOuter_, outer + 1);
                if (chS >= 1.0 && outer < 9) {   // a gas phase appeared (solvent entered undersaturated oil): flow first
                    int its2 = 0, lin2 = 0;
                    ok = newton(dt, its2, lin2); totIts += its2; totLin += lin2;
                    its = std::max(its, its2);
                    if (!ok) break;
                    continue;
                }
                const double ch = poly_ ? polymerTransport(dt) : 0.0;
                stats_.polyOuter++;
                if (ch < 0) { ok = false; log_("  polymer transport did not converge"); break; }
                const double chC = chem_.empty() ? 0.0 : chemTransport(dt);
                if (chC < 0) { ok = false; break; }
                chOuterMax_ = std::max(chOuterMax_, outer + 1);
                if (ch <= tol && chC <= tolC && chS <= tolS) break;
                if (outer >= 9 && ch <= tol && chC <= 0.02 && chS <= 1e-2) { log_(fmt("  chemical / solvent coupling: change %.2e after 10 outer iterations - step accepted", std::max(chC, chS))); break; }
                if (outer >= 9) { ok = false; log_(fmt("  polymer / chemical / solvent coupling not converged (change %.2e)", std::max(ch, std::max(chC, chS)))); break; }
                int its2 = 0, lin2 = 0;
                ok = newton(dt, its2, lin2); totIts += its2; totLin += lin2;
                its = std::max(its, its2);
                if (!ok) break;
            }
        }
        if (ok && hfSubgrid_) {
            // C21 (H2) implicit sub-grid tip: find the wing lengths L >= L_n with K_I(L, p(L)) = K_IC at the end of the step
            // (bracket expansion, then safeguarded regula falsi); every trial re-solves the step from its start
            const GridInput& g = case_.grid;
            const double toM = case_.units == UnitSystem::Field ? 0.3048 : 1.0;
            const int nw = int(hfTip_.size());
            std::vector<double> lo(tip0), hi(size_t(nw), 1e300), rlo(size_t(nw), 0.0), rhi(size_t(nw), 0.0);
            const double tolR = 0.01, tolL = 0.05 / toM;   // 1 % in K_I or a 0.05 m bracket
            for (int trial = 0; trial < 20 && ok; ++trial) {
                std::vector<double> next = hfTip_;
                bool moved = false;
                for (int w = 0; w < nw; ++w) {
                    const int f = w / 2, side = (w % 2) ? 1 : -1;
                    const auto& F = g.hfracs[size_t(f)];
                    if (!F.subgrid || !(F.kic > 0)) continue;
                    double KI = 0, pnet = 0;
                    if (!hfTipKI(f, side, KI, pnet)) continue;
                    const double r = KI / (F.kic * 1e6) - 1.0, L = hfTip_[size_t(w)];
                    if (std::getenv("RESSIM_HFDBG")) std::fprintf(stderr, "  [tip t=%.3f dt=%.3f w=%d L=%.3f r=%.3f pnet=%.3f]\n", time_, dt, w, L, r, pnet);
                    if (r > tolR) { lo[size_t(w)] = L; rlo[size_t(w)] = r; }
                    else if (r < -tolR && L > tip0[size_t(w)] + 1e-12) { hi[size_t(w)] = L; rhi[size_t(w)] = r; }
                    else continue;                          // converged, or no growth this step
                    if (L >= hfTipMax_[size_t(w)] - 1e-9 && r > 0) continue;   // wing at its maximum extent
                    double Ln;
                    if (hi[size_t(w)] < 1e299) {
                        const double a0 = lo[size_t(w)], b0 = hi[size_t(w)];
                        if (b0 - a0 < tolL) continue;
                        Ln = a0 + (b0 - a0) * rlo[size_t(w)] / std::max(rlo[size_t(w)] - rhi[size_t(w)], 1e-12);
                        Ln = std::clamp(Ln, a0 + 0.1 * (b0 - a0), b0 - 0.1 * (b0 - a0));
                    } else {
                        double piece = 0;   // expansion: one tip piece, then doubling
                        for (size_t k = 0; k < hfFrac_.size(); ++k) if (g.hfPieces[k].frac == f) piece = std::max(piece, 2 * hfHalf_[k]);
                        Ln = L + std::max(piece, L - tip0[size_t(w)]);
                    }
                    Ln = std::min(Ln, hfTipMax_[size_t(w)]);
                    if (std::fabs(Ln - L) > 1e-9) { next[size_t(w)] = Ln; moved = true; }
                }
                if (!moved) break;
                p_ = p0; sw_ = s0; x_ = x0; hasGas_ = h0; pc_ = pc0; wmul_ = wm0; chemRestore(); if (solv_) zs_ = zsBak_;
                for (size_t w = 0; w < wells_.size(); ++w) { wells_[w].bhp = wb[w].bhp; wells_[w].mode = wb[w].mode; }
                hfTip_ = next; setHfTips();
                computeHeads();
                ok = newton(dt, its, lin); totIts += its; totLin += lin;
                for (int sw = 0; ok && sw < 4 && checkWellSwitching(); ++sw) { ok = newton(dt, its, lin); totIts += its; totLin += lin; }
            }
            if (ok)
                for (int w = 0; w < nw; ++w)
                    if (hfTip_[size_t(w)] > tip0[size_t(w)] + 1e-9) {
                        const int f = w / 2, side = (w % 2) ? 1 : -1;
                        double KI = 0, pnet = 0; hfTipKI(f, side, KI, pnet);
                        HfEvent ev; ev.t = time_ + dt; ev.frac = f; ev.side = side; ev.length = hfTip_[size_t(w)]; ev.KI = KI; ev.pnet = pnet;
                        hfEvents_.push_back(ev);
                    }
        }
        stats_.newton += totIts; stats_.linear += totLin;
        if (!ok) {
            p_ = p0; sw_ = s0; x_ = x0; hasGas_ = h0; pc_ = pc0; wmul_ = wm0; chemRestore(); if (solv_) zs_ = zsBak_;
            if (hfSubgrid_ && hfTip_ != tip0) { hfTip_ = tip0; setHfTips(); }
            for (size_t w = 0; w < wells_.size(); ++w) { wells_[w].bhp = wb[w].bhp; wells_[w].mode = wb[w].mode; }
            dtNext_ = dt * T.cutFactor; stats_.cuts++;
            log_(fmt("  timestep cut: dt %.4g -> %.4g days", dt, dtNext_));
            if (dtNext_ < T.minDt) throw std::runtime_error(fmt("Time step below minimum at t = %.4f days - simulation stopped", time_));
            continue;
        }
        if (adjRec_) adjCapture(time_, dt);   // converged state + the lagged data this step used (adjoint replay)
        for (AquiferRT& Q : aquifers_) {
            Q.We += Q.rate * dt; Q.WeS += Q.rateS * dt;
            if (Q.type == AquiferType::Fetkovich) Q.pa -= Q.rate * dt / Q.ctV;
        }
        updateRockHistory();                  // C34 irreversible rock compaction history
        tracerTransport(dt);                  // C27 passive tracers on the converged fluxes
        saltStep(dt);                         // C29 salt precipitation / dry-out (after the NACL transport)
        storeOld();
        time_ += dt;
        stats_.timesteps++;
        record(dt, totIts, totLin);
        econCheck();
        pweTrack();
        if (econStop_) break;
        double next = dt;
        if (its <= 5) next = dt * T.growth; else if (its > 9) next = dt * 0.6;
        if (T.dsTarget > 0 || T.dpTarget > 0) {
            // C27: target-change control (CMG NORM / Aziz-Settari): f = (1 + w) eta / (d + w eta), w = 1, per quantity
            const int na = grid_.numActive();
            double dP = 0, dS = 0;
            auto sgOf = [&](const std::vector<double>& sw, const std::vector<double>& x, const std::vector<uint8_t>& h, int a) {
                const uint8_t st = size_t(a) < h.size() ? h[size_t(a)] : 0;
                return st == 1 ? x[size_t(a)] : (st == 2 ? 1.0 - (sw.empty() ? 0.0 : sw[size_t(a)]) : 0.0);
            };
            for (int a = 0; a < na; ++a) {
                dP = std::max(dP, std::fabs(p_[size_t(a)] - p0[size_t(a)]));
                if (!sw_.empty()) dS = std::max(dS, std::fabs(sw_[size_t(a)] - s0[size_t(a)]));
                if (!x_.empty() && !hasGas_.empty()) dS = std::max(dS, std::fabs(sgOf(sw_, x_, hasGas_, a) - sgOf(s0, x0, h0, a)));
            }
            double f = T.growth;
            if (T.dpTarget > 0) f = std::min(f, 2.0 * T.dpTarget / (dP + T.dpTarget));
            if (T.dsTarget > 0) f = std::min(f, 2.0 * T.dsTarget / (dS + T.dsTarget));
            next = std::min(next, dt * std::max(f, T.cutFactor));
            stats_.lastDP = dP; stats_.lastDS = dS;
        }
        if (clipped) next = std::max(next, dtNext_);
        dtNext_ = std::min(next, T.maxDt);
        // C18: LEFM growth test at the converged state; a newly opened piece fills from the fracture in the next steps, so
        // the step is not allowed to grow while the fracture is moving
        if (!hfOpen_.empty() && growFractures()) dtNext_ = std::min(dtNext_, dt);
        if (hfSubgrid_ && hfTip_ != tip0) dtNext_ = std::min(dtNext_, dt);   // C21: no step growth while a tip moves
    }
    if (!pfLinks_.empty()) pfReportStep();   // C28 in-run phase-field growth at the end of the report step
    ++step_;
    actionsEndOfStep();                       // C31 ACTIONX: conditions at the end of the report step, body acts from the next
    pyActionsEndOfStep();                     // C37 PYACTION scripts (same boundary, after the ACTIONX blocks)
    snapshot();
    if (finished() && mswPicardSteps_ > 0) {   // C49: converged multisegment heads - one summary line at the end of the run
        char b[200];
        std::snprintf(b, sizeof b, "Multisegment heads: %d time steps checked, %d extra solves; %d steps ended at the pass limit (bracket still open)",
                      mswPicardSteps_, mswPicardIts_, mswPicardUnconv_);
        log_(b);
    }
    stats_.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

void Simulator::record(double dt, int its, int lin) {
    SummaryPoint s;
    s.time = time_; s.dt = dt; s.newton = its; s.linear = lin;
    for (auto& w : wells_) {
        WellRates r = w.rates; r.bhp = w.bhp; r.mode = w.mode; r.open = w.open; r.eff = w.eff;
        r.stopped = !w.open && w.ctrl.status == WellStatus::Stop && w.nActive > 0 && !w.cell.empty();
        if (w.open) {   // C46 WPI vectors
            const auto P = wellPIs(size_t(&w - wells_.data()));
            for (int p = 0; p < 3; ++p) r.pi[p] = P[size_t(p)];
            const std::string& ph = case_.wells[size_t(w.input)].prefPhase;
            r.piPref = w.type != WellType::Producer ? P[1] + P[2] : ph == "GAS" ? P[2] : ph == "WATER" ? P[1] : ph == "LIQ" ? P[0] + P[1] : P[0];
        }
        r.thpNet = w.thpNet > 0 ? w.thpNet : 0.0;
        r.dpWellbore = 0; for (double d : w.dpFric) if (std::fabs(d) > std::fabs(r.dpWellbore)) r.dpWellbore = d;
        s.wells.push_back(r);
        s.fopr += r.eff * r.oil; s.fwpr += r.eff * r.water; s.fgpr += r.eff * r.gas; s.fwir += r.eff * r.waterInj; s.fgir += r.eff * r.gasInj;
        s.fnir += r.eff * r.solvInj; s.fnpr += r.eff * r.solvProd;
        s.wells.back().lift = liftGasRate(w); s.fglir += r.eff * s.wells.back().lift;
    }
    gloCum_ += s.fglir * dt; s.fglit = gloCum_;
    if (!gsat_.empty()) { s.fopr += satRate("FIELD", 0); s.fwpr += satRate("FIELD", 1); s.fgpr += satRate("FIELD", 2); }   // C35 GSATPROD
    for (const auto& [g, C] : gcons_) { s.fgcr += C.cons; s.fgimr += C.imp; }
    s.fgcr += compressorGas();   // C36 NETCOMPA fuel gas                                            // C35 GCONSUMP
    s.fgsr = s.fgpr - s.fgir - s.fgcr + s.fgimr;   // gas sales (ECLIPSE FGSR)
    gasBalCum_[0] += s.fgcr * dt; gasBalCum_[1] += s.fgimr * dt; gasBalCum_[2] += s.fgsr * dt;
    s.fgct = gasBalCum_[0]; s.fgimt = gasBalCum_[1]; s.fgst = gasBalCum_[2];
    if (solv_) { svCum_[0] += s.fnir * dt; svCum_[1] += s.fnpr * dt; s.fnit = svCum_[0]; s.fnpt = svCum_[1]; s.fnip = solventInPlace(); }
    if (!saltS_.empty()) s.saltSolid = saltSolidMass();   // C29 SALT: precipitated NaCl, kg
    if (!netNodes_.empty()) { NetRecord R; R.t = time_; for (const auto& N : netNodes_) R.p.push_back(N.p); netHist_.push_back(std::move(R)); }
    if (!trInfo_.empty() || !chem_.empty()) {   // C27 tracers (+ C35 chemicals) ride on the component vectors: rates amount/s (+ produced, - injected per well)
        const size_t ncm = trInfo_.size() + chem_.size();
        s.fcompP.assign(ncm, 0.0); s.fcompI.assign(ncm, 0.0);
        for (WellRates& r : s.wells) {
            if (r.comp.size() < ncm) r.comp.resize(ncm, 0.0);
            for (size_t t = 0; t < r.comp.size() && t < ncm; ++t) { s.fcompP[t] += r.eff * std::max(0.0, r.comp[t]); s.fcompI[t] += r.eff * std::max(0.0, -r.comp[t]); }
        }
        accumulateComponentPoint(summary_.empty() ? nullptr : &summary_.back(), s);
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
        if (eq_[PG] >= 0) s.gip += cf_[a].acc[PG].v * (solv_ ? 1.0 - zs_[size_t(a)] : 1.0);   // C35: FGIP without the solvent
    }
    // ECLIPSE definitions: FPR = hydrocarbon-pore-volume weighted (PV weighted when there is no hydrocarbon), FPRP = PV weighted
    s.fprp = ppv / pv;
    s.fpr = (eq_[PO] >= 0 || eq_[PG] >= 0) && hc > 1e-12 * pv ? phc / hc : s.fprp;
    for (const AquiferRT& Q : aquifers_) { s.faqr += Q.rate; s.faqt += Q.We; s.faqrS += Q.rateS; s.faqtS += Q.WeS; }
    if (poly_) {
        for (const WellRates& r : s.wells) { s.fcir += r.eff * r.polyInj; s.fcpr += r.eff * r.polyProd; }   // C33 WEFAC
        pcCum_[0] += s.fcir * dt; pcCum_[1] += s.fcpr * dt;
        s.fcit = pcCum_[0]; s.fcpt = pcCum_[1];
        const PolymerInput& P = case_.polymer;
        for (int a = 0; a < na; ++a) {
            const PolymerRock& R = P.plyrock[std::min(size_t(satReg_[size_t(a)]), P.plyrock.size() - 1)];
            const double pv = fl(a).pvMult(p_[a]) * grid_.poreVolRef[a];
            s.fcip += cf_[a].acc[PW].v * pc_[a] * (1.0 - R.deadPoreVolume);
            s.fcad += std::max(0.0, grid_.bulkVol[size_t(grid_.globalOf[a])] - pv / U_.poreVolume) * R.rockDensity * polyEval(a, pc_[a]).ads;
        }
    }
    if (udq_) {   // C31: well cumulatives and UDQ values after the accepted step
        if (wCum_.size() != wells_.size() * 5) wCum_.assign(wells_.size() * 5, 0.0);
        for (size_t i = 0; i < wells_.size(); ++i) {
            const WellRates& r = wells_[i].rates; const double ef = wells_[i].eff;
            wCum_[i * 5] += ef * r.oil * dt; wCum_[i * 5 + 1] += ef * r.water * dt; wCum_[i * 5 + 2] += ef * r.gas * dt;
            wCum_[i * 5 + 3] += ef * r.waterInj * dt; wCum_[i * 5 + 4] += ef * r.gasInj * dt;
        }
        udqEvaluate(&s, &s);
    }
    if (!msw_.empty()) mswRecord();   // C32
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
        double rs = live_ ? (st == 0 ? x_[a] : rsSatC(a, p_[a])) : 0.0;
        if (live_ && st == 1 && !rsLim_.empty()) rs = std::min(rs, rsLim_[size_t(a)]);
        const double rv = vap_ ? (st == 2 ? x_[a] : rvSatC(a, p_[a])) : 0.0;
        s.pressure[a] = float(p_[a]); s.sw[a] = float(sw); s.sg[a] = float(sg);
        s.so[a] = float(O_ && st != 2 ? std::max(0.0, 1.0 - sw - sg) : 0.0); s.rs[a] = float(rs);
        if (vap_) s.rv[a] = float(rv);
    }
    if (solv_ && G_) {   // C35 SOLVENT: SGAS = free gas, SSOLVENT, Rs from the flash (SOIL = 1 - SWAT - SGAS - SSOLVENT unchanged)
        std::vector<float> ss(static_cast<size_t>(na)), zf(static_cast<size_t>(na));
        for (int a = 0; a < na; ++a) {
            const CellFull C = evalCell(a);
            s.sg[a] = float(C.S[PG].v); ss[size_t(a)] = float(C.ss.v); zf[size_t(a)] = float(zs_[size_t(a)]);
            if (live_) s.rs[a] = float(C.rs.v);
        }
        s.extra.emplace_back("SSOLVENT", std::move(ss));
        s.extra.emplace_back("SOLVFRAC", std::move(zf));   // solvent fraction of the G+S component
    }
    if (poly_) {
        s.poly.resize(na); s.polyAds.resize(na);
        for (int a = 0; a < na; ++a) { s.poly[a] = float(pc_[a]); s.polyAds[a] = float(polyEval(a, pc_[a]).ads); }
    }
    if (!hfOpen_.empty()) {   // C18: 1 = open hydraulic-fracture piece (growth runs), 0 = closed piece or matrix cell
        std::vector<float> o(size_t(na), 0.0f);
        for (size_t k = 0; k < hfOpen_.size(); ++k) {
            const int a = grid_.activeOf[size_t(grid_.hfCells[k].global)];
            if (a >= 0 && hfOpen_[k]) o[size_t(a)] = hfFrac_.empty() ? 1.0f : float(hfFrac_[k]);
        }
        s.extra.emplace_back("HFOPEN", std::move(o));
    }
    for (size_t t = 0; t < trInfo_.size(); ++t) {   // C27 tracer concentrations
        std::vector<float> v(static_cast<size_t>(na));
        for (int a = 0; a < na; ++a) v[size_t(a)] = float(trC_[t * size_t(na) + size_t(a)]);
        s.extra.emplace_back("TR_" + trInfo_[t].name, std::move(v));
    }
    for (const ChemSp& S : chem_) {                  // C35 chemical concentrations (+ surfactant capillary number / miscibility, foam factor)
        std::vector<float> v(static_cast<size_t>(na));
        for (int a = 0; a < na; ++a) v[size_t(a)] = float(S.c[size_t(a)]);
        s.extra.emplace_back(S.kind == 0 ? "SURFACT" : S.kind == 1 ? "ALKALINE" : "FOAM", std::move(v));
    }
    if (chSurf_ >= 0) {
        std::vector<float> v(static_cast<size_t>(na)), m(static_cast<size_t>(na)), pc(static_cast<size_t>(na));
        for (int a = 0; a < na; ++a) { v[size_t(a)] = float(chLogNc_[size_t(a)]); m[size_t(a)] = float(chMisc_[size_t(a)]); pc[size_t(a)] = float(chPcMul_[size_t(a)]); }
        s.extra.emplace_back("SURFCNM", std::move(v));   // log10 capillary number
        s.extra.emplace_back("SURFMISC", std::move(m));  // miscibility fraction (SURFCAPD)
        s.extra.emplace_back("SURFSTM", std::move(pc));  // surface tension / SURFST(0)
    }
    if (chFoam_ >= 0) {
        std::vector<float> v(static_cast<size_t>(na));
        for (int a = 0; a < na; ++a) v[size_t(a)] = float(chGmul_[size_t(a)]);
        s.extra.emplace_back("FOAMMRF", std::move(v));   // gas mobility multiplier
    }
    if (!saltS_.empty()) {   // C29 SALT: solid fraction of the pore space and permeability factor
        std::vector<float> v(static_cast<size_t>(na)), k(static_cast<size_t>(na));
        for (int a = 0; a < na; ++a) { v[size_t(a)] = float(saltS_[size_t(a)]); k[size_t(a)] = float(saltK_[size_t(a)]); }
        s.extra.emplace_back("SSALT", std::move(v));
        s.extra.emplace_back("KSALT", std::move(k));
    }
    if (!pfLinks_.empty()) pfFieldArrays(s);   // C30 phase-field damage / stress / opening
    if (hff_) hff_->addArrays(*this, s);      // C30 DDM opening + stress / displacement field of the fractures
    snaps_.push_back(std::move(s));
}

const hffield::Stats* Simulator::hfFieldStats() const { return hff_ ? &hff_->stats() : nullptr; }
bool Simulator::hfFieldAt(const double x[3], double S[6], double u[3]) const {
    if (!hff_) return false;
    hff_->fieldAt(x, S, u);
    return true;
}

bool Simulator::phaseFieldFrame(int frac, double& ox, double& oy, double& tx, double& ty) const {
    const GridInput& g = case_.grid;
    if (frac < 0 || size_t(frac) >= g.hfracs.size() || hfS_.size() != g.hfPieces.size()) return false;
    const auto& F = g.hfracs[size_t(frac)];
    const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0);
    if (!(L > 0)) return false;
    tx = (F.x1 - F.x0) / L; ty = (F.y1 - F.y0) / L;
    for (size_t k = 0; k < g.hfPieces.size(); ++k) {
        const auto& P = g.hfPieces[k];
        if (P.frac != frac) continue;
        ox = 0.5 * (P.xa + P.xb) - hfS_[k] * tx; oy = 0.5 * (P.ya + P.yb) - hfS_[k] * ty;
        return true;
    }
    return false;
}

// C30: the latest phase-field step of every PHASEFIELD fracture mapped onto the flow grid. The 2-D plane-strain model is
// the horizontal section through the fracture (local x along the strike from the injection point, y the horizontal
// normal); every node / element is assigned to the grid column with the nearest centre and the values are spread over the
// layers of the fracture height:
//   PFDAMAGE  -       largest damage 1 - v of the column (the regularised crack is narrower than a cell)
//   PFDSXX PFDSYY PFDSXY  p unit  mean stress change of the column, global axes, compression +
//   PFDSNN    p unit  normal-stress change on planes parallel to the crack (stress shadow)
//   PFDISP    mm      largest displacement of the intact material (v > 0.5) in the column
//   PFWIDTH   mm      crack opening w(s) = int u . grad v dy at the piece position (fracture and host cells)
void Simulator::pfFieldArrays(Snapshot& s) const {
    const size_t na = s.pressure.size();
    const GridInput& g = case_.grid;
    const int nx = grid_.nx, ny = grid_.ny, nz = grid_.nz;
    if (grid_.corners.size() < size_t(grid_.numCells())) return;
    std::vector<float> D(na, 0.0f), SX(na, 0.0f), SY(na, 0.0f), SXY(na, 0.0f), SN(na, 0.0f), U(na, 0.0f), W(na, 0.0f);
    bool any = false;
    const double pu = 101325.0 / Units::of(case_.units).atm;
    for (const auto& lp : pfLinks_) {
        const PfLink& P = *lp;
        std::shared_ptr<phasefield::Field> Fp;
        { std::lock_guard<std::mutex> lk(const_cast<std::mutex&>(P.m)); Fp = P.field; }
        if (!Fp || Fp->v.empty()) continue;
        const phasefield::Field& F = *Fp;
        double ox, oy, tx, ty;
        if (!phaseFieldFrame(P.frac, ox, oy, tx, ty)) continue;
        any = true;
        const double toM = P.toM;
        // layers of the fracture: depth range of its pieces
        double zt = 1e300, zb = -1e300;
        for (size_t k = 0; k < g.hfPieces.size(); ++k) if (g.hfPieces[k].frac == P.frac) { zt = std::min(zt, g.hfPieces[k].zt); zb = std::max(zb, g.hfPieces[k].zb); }
        // column centres (plan) and the layers inside [zt, zb]
        std::vector<double> cx(size_t(nx) * size_t(ny), 0.0), cy(cx.size(), 0.0), csz(cx.size(), 0.0);
        std::vector<std::vector<int>> colCells(cx.size());
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const size_t col = size_t(i) + size_t(nx) * size_t(j);
            int nk = 0;
            for (int k = 0; k < nz; ++k) {
                const int gq = i + nx * (j + ny * k);
                const auto& C = grid_.corners[size_t(gq)];
                double x = 0, y = 0, z = 0;
                for (int v = 0; v < 8; ++v) { x += C[size_t(v)].x; y += C[size_t(v)].y; z += C[size_t(v)].z; }
                x /= 8; y /= 8; z /= 8;
                if (nk == 0 || (z >= zt && z <= zb)) { if (nk == 0) { cx[col] = x; cy[col] = y; } }
                if (z >= zt - 1e-9 && z <= zb + 1e-9) { colCells[col].push_back(gq); cx[col] = x; cy[col] = y; ++nk; }
            }
            const int g0 = i + nx * j;
            csz[col] = std::max(grid_.cellDx.size() > size_t(g0) ? grid_.cellDx[size_t(g0)] : 1.0, grid_.cellDy.size() > size_t(g0) ? grid_.cellDy[size_t(g0)] : 1.0);
        }
        // uniform bins of the column centres
        double bx0 = 1e300, by0 = 1e300, bx1 = -1e300, by1 = -1e300, hs = 0;
        for (size_t c = 0; c < cx.size(); ++c) { bx0 = std::min(bx0, cx[c]); bx1 = std::max(bx1, cx[c]); by0 = std::min(by0, cy[c]); by1 = std::max(by1, cy[c]); hs = std::max(hs, csz[c]); }
        const double bs = std::max(hs, 1e-9);
        const int nbx = std::max(1, int((bx1 - bx0) / bs) + 1), nby = std::max(1, int((by1 - by0) / bs) + 1);
        std::vector<std::vector<int>> bins(size_t(nbx) * size_t(nby));
        for (size_t c = 0; c < cx.size(); ++c) {
            if (colCells[c].empty()) continue;
            const int bi = std::clamp(int((cx[c] - bx0) / bs), 0, nbx - 1), bj = std::clamp(int((cy[c] - by0) / bs), 0, nby - 1);
            bins[size_t(bi) + size_t(nbx) * size_t(bj)].push_back(int(c));
        }
        auto columnOf = [&](double X, double Y) {   // local metres -> column (-1 outside)
            const double x = ox + (X * tx - Y * ty) / toM, y = oy + (X * ty + Y * tx) / toM;
            const int bi = int(std::floor((x - bx0) / bs)), bj = int(std::floor((y - by0) / bs));
            int best = -1; double bd = 1e300;
            for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di) {
                const int a = bi + di, b = bj + dj;
                if (a < 0 || b < 0 || a >= nbx || b >= nby) continue;
                for (int c : bins[size_t(a) + size_t(nbx) * size_t(b)]) {
                    const double d = std::hypot(cx[size_t(c)] - x, cy[size_t(c)] - y);
                    if (d < bd) { bd = d; best = c; }
                }
            }
            return best >= 0 && bd <= 0.75 * csz[size_t(best)] ? best : -1;
        };
        const int nx1 = int(F.x.size()), ny1 = int(F.y.size());
        std::vector<float> cd(cx.size(), 0.0f), cu(cx.size(), 0.0f);
        std::vector<double> sxx(cx.size(), 0.0), syy(cx.size(), 0.0), sxy(cx.size(), 0.0), snn(cx.size(), 0.0), cnt(cx.size(), 0.0);
        for (int j = 0; j < ny1; ++j) for (int i = 0; i < nx1; ++i) {
            const int col = columnOf(F.x[size_t(i)], F.y[size_t(j)]);
            if (col < 0) continue;
            const size_t q = size_t(i) + size_t(nx1) * size_t(j);
            cd[size_t(col)] = std::max(cd[size_t(col)], 1.0f - F.v[q]);
            if (F.v[q] > 0.5f)   // inside the broken band the displacement of the (stiffness k) material is not defined
                cu[size_t(col)] = std::max(cu[size_t(col)], float(std::hypot(F.ux[q], F.uy[q]) * 1e3));
        }
        const double c2 = tx * tx, s2 = ty * ty, cs = tx * ty;
        for (int j = 0; j + 1 < ny1; ++j) for (int i = 0; i + 1 < nx1; ++i) {
            const size_t e = size_t(i) + size_t(nx1 - 1) * size_t(j);
            if (e >= F.sxx.size()) continue;
            const int col = columnOf(0.5 * (F.x[size_t(i)] + F.x[size_t(i) + 1]), 0.5 * (F.y[size_t(j)] + F.y[size_t(j) + 1]));
            if (col < 0) continue;
            const double a = F.sxx[e], b = F.syy[e], t = F.sxy[e];   // local frame, tension +
            sxx[size_t(col)] += a * c2 + b * s2 - 2 * t * cs;
            syy[size_t(col)] += a * s2 + b * c2 + 2 * t * cs;
            sxy[size_t(col)] += (a - b) * cs + t * (c2 - s2);
            snn[size_t(col)] += b; cnt[size_t(col)] += 1;
        }
        for (size_t c = 0; c < cx.size(); ++c)
            for (int gq : colCells[c]) {
                const int a = grid_.activeOf[size_t(gq)];
                if (a < 0 || size_t(a) >= na) continue;
                D[size_t(a)] = std::max(D[size_t(a)], cd[c]); U[size_t(a)] = std::max(U[size_t(a)], cu[c]);
                if (cnt[c] > 0) {
                    SX[size_t(a)] = float(-sxx[c] / cnt[c] / pu); SY[size_t(a)] = float(-syy[c] / cnt[c] / pu);
                    SXY[size_t(a)] = float(-sxy[c] / cnt[c] / pu); SN[size_t(a)] = float(-snn[c] / cnt[c] / pu);
                }
            }
        // opening profile on the pieces (fracture cells and their hosts)
        if (!F.cod.empty())
            for (size_t k = 0; k < g.hfPieces.size() && k < hfS_.size(); ++k) {
                if (g.hfPieces[k].frac != P.frac) continue;
                const double X = hfS_[k] * toM;
                const auto it = std::lower_bound(F.codX.begin(), F.codX.end(), float(X));
                double w = 0;
                if (it == F.codX.begin()) w = F.cod.front();
                else if (it == F.codX.end()) w = F.cod.back();
                else {
                    const size_t r = size_t(it - F.codX.begin());
                    const double x0 = F.codX[r - 1], x1 = F.codX[r], f = (X - x0) / std::max(1e-12, x1 - x0);
                    w = F.cod[r - 1] + f * (F.cod[r] - F.cod[r - 1]);
                }
                const float wmm = float(std::max(0.0, w) * 1e3);
                const int a = grid_.activeOf[size_t(grid_.hfCells[k].global)];
                if (a >= 0 && size_t(a) < na) W[size_t(a)] = wmm;
                const int h = grid_.hfCells[k].host;
                const int ah = h >= 0 && size_t(h) < grid_.activeOf.size() ? grid_.activeOf[size_t(h)] : -1;
                if (ah >= 0 && size_t(ah) < na) W[size_t(ah)] = std::max(W[size_t(ah)], wmm);
            }
    }
    (void)any;
    s.extra.emplace_back("PFDAMAGE", std::move(D));
    s.extra.emplace_back("PFDSXX", std::move(SX));
    s.extra.emplace_back("PFDSYY", std::move(SY));
    s.extra.emplace_back("PFDSXY", std::move(SXY));
    s.extra.emplace_back("PFDSNN", std::move(SN));
    s.extra.emplace_back("PFDISP", std::move(U));
    s.extra.emplace_back("PFWIDTH", std::move(W));
}

double Simulator::totalSimTime() const {
    double t = 0; for (auto& r : case_.schedule) t += r.length; return t;
}

double Simulator::adsorbedGas() const {
    double v = 0;
    if (langK_.empty()) return 0;
    for (int a = 0; a < grid_.numActive(); ++a) v += langK_[size_t(a)] * p_[size_t(a)] / (langPL_ + p_[size_t(a)]);
    return v;
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
    ip.gasAdsorbed = adsorbedGas();
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
                if (gravdr_ && p != PO && gdZ_[size_t(ci)] > 0) {   // GRAVDR column heads, as in the assembly
                    const double rhoO = 0.5 * (A.rho[PO].v + B.rho[PO].v), drho = p == PW ? rho - rhoO : rhoO - rho;
                    const int hp = p == PW ? 0 : 1;
                    dphi += drho * (gdH_[size_t(cn.a) * 2 + hp].v - gdH_[size_t(cn.b) * 2 + hp].v) * (g * gdZ_[size_t(ci)]);
                }
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
                if (gravdr_ && p != PO && gdZ_[size_t(ci)] > 0) {   // GRAVDR column heads, as in the assembly
                    const double rhoO = 0.5 * (A.rho[PO].v + B.rho[PO].v), drho = p == PW ? rho - rhoO : rhoO - rho;
                    const int hp = p == PW ? 0 : 1;
                    dphi += drho * (gdH_[size_t(cn.a) * 2 + hp].v - gdH_[size_t(cn.b) * 2 + hp].v) * (g * gdZ_[size_t(ci)]);
                }
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
            out.push_back({W.cell[k], q * W.eff, w});   // C33 WEFAC: time-averaged reservoir flux
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
                    // C38: SWATINIT below the (scaled) connate saturation is raised to SWL, SWL itself is not changed
                    // (OPM EclMaterialLawManager::applySwatinit)
                    const auto& ec = ep(a) ? *ep(a) : S.tableEndPoints();
                    const double sIn = swi[size_t(grid_.globalOf[size_t(a)])];
                    const double s0 = sIn >= 0 ? std::max(ec.swl, sIn) : sIn;
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
    // C34 numerical aquifer cells with an AQUNUM initial pressure (default: equilibrium at the aquifer depth)
    for (auto& [a, p0] : grid_.numAquP0) if (p0 > 0 && a >= 0 && a < na) p_[size_t(a)] = p0;
    setupRockComp();   // C34 (needs the initial pressure / saturation)
    if (!grid_.numAquP0.empty()) log_(fmt("numerical aquifers: %g cells", double(grid_.numAquP0.size())) + fmt(", %g reservoir connections", double(grid_.numAquConns)));
    polymerSetup();
    solventSetup();                     // C35 SOLVENT
    chemSetup();                        // C35 surfactant / alkaline / foam
    if (!sGam_.empty()) pI_ = p_;   // C15: stress-dependent permeability is measured from the initial pressure
    evalAll();
    setupThresholds();
    storeOld();
    time_ = 0;
    initAquifers();
    summary_.clear(); snaps_.clear(); netHist_.clear(); netStep_ = -1;
    for (double& v : cum_) v = 0;
    pcCum_[0] = pcCum_[1] = 0;
    for (auto& w : wells_) { w.rates = WellRates{}; }
    tracerSetup();
    saltSetup();                       // C29 SALT (after the tracers: dissolved NaCl is the internal NACL tracer)
    record(0, 0, 0);
    if (adjRec_) { adj_.clear(); adjCapture(0.0, 0.0); }
    snapshot();
    const InPlace ip = inPlace();
    log_(std::string("Initial fluids in place (") + U_.liq + " / " + U_.gas + "):");
    log_(fmt("  oil %.6g   water %.6g   free gas %.6g   dissolved gas %.6g", ip.oil, ip.water, ip.gasFree, ip.gasDissolved));
    if (vap_) log_(fmt("  vaporised oil (in gas) %.6g   total oil %.6g", ip.oilVaporized, ip.oil + ip.oilVaporized));
    if (!langK_.empty()) log_(fmt("  adsorbed gas (Langmuir) %.6g   free + dissolved + adsorbed %.6g", ip.gasAdsorbed, ip.gasFree + ip.gasDissolved + ip.gasAdsorbed));
    if (grid_.numHfCells()) log_(fmt("  hydraulic fractures: %g fracture cells (EDFM)", double(grid_.numHfCells())));
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


// ---------------------------------------------------------------- C14 restart state
template <class Ar> void Simulator::stateIO(Ar& ar) {
    int na = grid_.numActive(); uint64_t nw = wells_.size(), naq = aquifers_.size();
    const int na0 = na; const uint64_t nw0 = nw, naq0 = naq;
    ar(na, nw, naq);
    if (Ar::reading && (na != na0 || nw != nw0 || naq != naq0))
        throw std::runtime_error(fmt("restart state does not match the case: %g active cells", double(na)) +
                                 fmt(" / %g wells", double(nw)) + fmt(" / %g aquifers in the file vs ", double(naq)) +
                                 fmt("%g / ", double(na0)) + fmt("%g / ", double(nw0)) + fmt("%g in the case", double(naq0)));
    ar(time_, dtNext_, step_);
    for (double& v : cum_) ar(v);
    ar(pcCum_[0], pcCum_[1]);
    ar(p_, sw_, x_, hasGas_);
    for (WellRT& w : wells_) {
        ar(w.bhp, w.mode, w.open, w.dpConn, w.headRho, w.qwConn, w.wmulConn, w.zwConn);
        WellRates& r = w.rates;
        ar(r.oil, r.water, r.gas, r.waterInj, r.gasInj, r.resv, r.polyInj, r.polyProd, r.bhp, r.oilH, r.waterH, r.gasH, r.mode, r.open);
        ar(r.comp, r.compCumP, r.compCumI);
    }
    for (AquiferRT& Q : aquifers_) ar(Q.p0, Q.pa, Q.We, Q.rate, Q.WeS, Q.rateS);
    ar(sgHy_, soHy_, dG_, dO_, soMaxV_, fRv_, fRs_);
    ar(pc_, pcAdsMax_, pcMassOld_, wmul_, pcZw_, pcZp_);
    ar(lastRs_, rsLim_, quasiFallback_);
    // RunStats: raw POD in state version 1 (C14-C17 layout, ends before the C18 GPU fields); size-prefixed from version 2
    if constexpr (Ar::reading) {
        if (stateVer_ < 2) ar.r.raw(&stats_, offsetof(RunStats, gpuWallSec));
        else {
            const uint64_t n = ar.r.u64();
            ar.r.need(size_t(n));
            std::memcpy(&stats_, ar.r.p, std::min<size_t>(size_t(n), sizeof(RunStats))); ar.r.p += n;
        }
    } else { ar.w.u64(sizeof(RunStats)); ar.w.raw(&stats_, sizeof(RunStats)); }
    if (!sGam_.empty()) ar(pI_);   // C15 stress-dependent permeability reference (only cases that use it)
    if (stateVer_ >= 2) ar(hfOpen_, hfEvents_);   // C18 fracture growth state (empty vectors without growth)
    if (stateVer_ >= 3) ar(hfTip_, hfFrac_);      // C21 sub-grid tips (empty without)
    if (stateVer_ >= 4) {                         // C27 economic-limit state
        for (WellRT& w : wells_) ar(w.econ, w.econNextTest, w.econTests);
        ar(econStop_);
        ar(trC_, trDecayCum_);                    // C27 tracer concentrations + decayed amounts
    }
    if (stateVer_ >= 5) {                         // C29 salt state, wellbore friction, network THP
        ar(saltS_, saltEvap_, saltEvapCum_);
        for (WellRT& w : wells_) ar(w.dpFric, w.dpF, w.thpNet);
    }
    if (stateVer_ >= 6) {                         // C31 UDQ values, well cumulatives, action runs and applied action bodies
        ar(wCum_);
        std::string ub = (!Ar::reading && udq_) ? udq_->save() : std::string();
        ar(ub);
        if (Ar::reading && udq_ && !ub.empty()) udq_->load(ub);
        std::vector<int> runs; std::vector<double> last;
        if (!Ar::reading && actions_) for (const auto& a : *actions_) { runs.push_back(a.runs); last.push_back(a.lastRun); }
        ar(runs, last);
        if (Ar::reading && actions_)
            for (size_t i = 0; i < actions_->size() && i < runs.size() && i < last.size(); ++i) { (*actions_)[i].runs = runs[i]; (*actions_)[i].lastRun = last[i]; }
        ar(actStep_, actText_);
    }
    if (stateVer_ >= 7) {                         // C34 WEFAC, guide rates, WELPI / econ-closed connections, rock-compaction history
        for (WellRT& w : wells_) ar(w.rates.eff, w.eff, w.guide, w.guideT, w.wiBase, w.cUntil);
        ar(rockPmin_, rockSwMax_);
    }
    if (stateVer_ >= 8) {                         // C35 chemical EOR: concentrations, adsorption history, coupling multipliers
        for (ChemSp& S : chem_) { ar(S.c, S.adsMax, S.adsMul); ar(S.cumDecay); }
        ar(chWmul_, chGmul_, chMisc_, chPcMul_, polyAdsMul_, chVg_, chLogNc_);
    }
    if (stateVer_ >= 9) {                         // C35 SOLVENT: solvent fraction of G+S, field totals
        ar(zs_, ntOld_); ar(svCum_[0]); ar(svCum_[1]);
    }
    if (stateVer_ >= 10) {                        // C35 gas lift: ALQ in force, last optimisation, lift gas total
        for (WellRT& w : wells_) ar(w.alq, w.alqLift);
        ar(gloLast_, gloCum_, gloRuns_);
        // the BHP limit of THP-controlled wells (VFP operating point, updated in the Newton iterations) and the group-allocated
        // rate target carry over from the last time step: the replay of the schedule resets them to the deck values
        for (WellRT& w : wells_) { ar(w.ctrl.bhp, w.ctrl.rate); bool g = w.grpCtl; ar(g); w.grpCtl = g; }
        // C35 facilities: gas consumption / import / sales totals, automatic choke pressure drops of the network nodes
        ar(gasBalCum_[0]); ar(gasBalCum_[1]); ar(gasBalCum_[2]);
        std::vector<double> dc; for (const NetNode& N : netNodes_) dc.push_back(N.dpChoke);
        ar(dc);
        for (size_t i = 0; i < dc.size() && i < netNodes_.size(); ++i) netNodes_[i].dpChoke = dc[i];
    }
    if (stateVer_ >= 11) {                        // C36 network: last balance time (NETBALAN interval), compressor levels
        ar(netLastBal_);
        std::vector<int> lv; for (const NetNode& N : netNodes_) lv.push_back(N.compLevel);
        ar(lv);
        for (size_t i = 0; i < lv.size() && i < netNodes_.size(); ++i) netNodes_[i].compLevel = lv[i];
    }
    if (stateVer_ >= 12) {                        // C37 PYACTION call counts / finished flags (the Python globals are not saved)
        std::vector<int> calls, done;
        if (!Ar::reading && pyActions_) for (const auto& a : *pyActions_) { calls.push_back(a.calls); done.push_back(a.done ? 1 : 0); }
        ar(calls, done);
        if (Ar::reading && pyActions_)
            for (size_t i = 0; i < pyActions_->size() && i < calls.size() && i < done.size(); ++i) { (*pyActions_)[i].calls = calls[i]; (*pyActions_)[i].done = done[i] != 0; }
    }
    if (stateVer_ >= 13) {                        // C37 WCYCLE open / shut times
        for (WellRT& w : wells_) { ar(w.cycOpen, w.cycClose); bool sc = w.cycSched; ar(sc); w.cycSched = sc; }
    }
    if (stateVer_ >= 14) {                        // C41 cross-flow wellbore mixture of the last assembly
        for (WellRT& w : wells_) ar(w.xmix[0], w.xmix[1], w.xmix[2]);
    }
    if (stateVer_ >= 15) {                        // C44: economic connection closures (WTEST 'E' re-opening), WPWE baseline, injector recirculation
        for (WellRT& w : wells_) {
            std::vector<int> ec, eu;
            if (!Ar::reading) for (const auto& [cell, until] : w.econConns) { ec.push_back(cell); eu.push_back(until); }
            ar(ec, eu);
            if (Ar::reading) { w.econConns.clear(); for (size_t i = 0; i < ec.size() && i < eu.size(); ++i) w.econConns[ec[i]] = eu[i]; }
            bool b = w.pweBase, inj = w.pweInj; ar(b, inj); w.pweBase = b; w.pweInj = inj;
            ar(w.pweOpenCells);
            std::vector<int> last(w.pweLast, w.pweLast + 8); ar(last);
            for (size_t i = 0; i < 8 && i < last.size(); ++i) w.pweLast[i] = last[i];
            ar(w.xrec[0], w.xrec[1], w.xrec[2]); ar(w.xinjTot);
            if (stateVer_ >= 16) ar(w.xmixConn);   // C45 MSW segment mixtures
        }
    }
}

std::string Simulator::saveState() const {
    StateOut o;
    o.w.str("RESSIM-BO-STATE-16");
    const_cast<Simulator*>(this)->stateVer_ = 16;
    const_cast<Simulator*>(this)->stateIO(o);
    return o.w.buf;
}

void Simulator::loadState(const std::string& blob) {
    StateIn in(blob);
    {
        const std::string tag = in.r.str();
        if (tag == "RESSIM-BO-STATE-1") stateVer_ = 1;
        else if (tag == "RESSIM-BO-STATE-2") stateVer_ = 2;
        else if (tag == "RESSIM-BO-STATE-3") stateVer_ = 3;
        else if (tag == "RESSIM-BO-STATE-4") stateVer_ = 4;
        else if (tag == "RESSIM-BO-STATE-5") stateVer_ = 5;
        else if (tag == "RESSIM-BO-STATE-6") stateVer_ = 6;
        else if (tag == "RESSIM-BO-STATE-7") stateVer_ = 7;
        else if (tag == "RESSIM-BO-STATE-8") stateVer_ = 8;
        else if (tag == "RESSIM-BO-STATE-9") stateVer_ = 9;
        else if (tag == "RESSIM-BO-STATE-10") stateVer_ = 10;
        else if (tag == "RESSIM-BO-STATE-11") stateVer_ = 11;
        else if (tag == "RESSIM-BO-STATE-12") stateVer_ = 12;
        else if (tag == "RESSIM-BO-STATE-13") stateVer_ = 13;
        else if (tag == "RESSIM-BO-STATE-14") stateVer_ = 14;
        else if (tag == "RESSIM-BO-STATE-15") stateVer_ = 15;
        else if (tag == "RESSIM-BO-STATE-16") stateVer_ = 16;
        else throw std::runtime_error("not a black-oil restart state");
    }
    const auto hfOpenCase = hfOpen_;
    const auto hfTipCase = hfTip_, hfFracCase = hfFrac_;
    stateIO(in);
    if (actions_ && (!baseSchedule_.empty() || !actStep_.empty())) {   // C31: the schedule as modified by the actions that fired before the restart (C33: or by script input)
        if (baseSchedule_.empty()) { baseSchedule_ = case_.schedule; baseWells_ = case_.wells; }
        case_.schedule = baseSchedule_; case_.wells = baseWells_;
        auto fresh = std::make_shared<udq::Set>();
        for (const ExtraKeyword& e : case_.extraKeywords) if (e.name == "UDQ") { std::string err; fresh->add(e.body, e.step, err); }
        udq_ = fresh;
        const auto st = actStep_; const auto tx = actText_;
        for (size_t k = 0; k < st.size() && k < tx.size(); ++k) actionApply(st[k], tx[k], true);
        refreshConnectionWindows(true);
    }
    // the schedule events of the report steps before the restart (well controls, connection windows, polymer / DRSDT
    // settings, group controls) define persistent run-time data: replay them silently, then restore the dynamic state
    // again (the replay re-initialises BHPs / modes of opening wells)
    if (step_ > int(case_.schedule.size())) throw std::runtime_error("restart step beyond the schedule of this case");
    {
        auto keepLog = log; log = nullptr;
        for (int s = 0; s < step_; ++s) { applyControls(s); polymerControls(s); solventControls(s); chemControls(s); drsdtControls(s); groupControls(s); gasLiftControls(s); }
        log = keepLog;
        StateIn again(blob);
        again.r.str();
        stateIO(again);
    }
    if (hfOpen_.size() != hfOpenCase.size()) hfOpen_ = hfOpenCase;
    if (stateVer_ < 6 && udq_) wCum_.assign(wells_.size() * 5, 0.0);   // v1 state / different fracture set: initial openings
    if (!pfLinks_.empty()) log_("  NOTE: restart - the phase-field fracture models (PHASEFIELD) are not part of the saved state; they restart from the "
                                "initial crack and integrate the injection after the restart (the open pieces of the base run are kept)");
    if (hfTip_.size() != hfTipCase.size() || hfFrac_.size() != hfFracCase.size()) { hfTip_ = hfTipCase; hfFrac_ = hfFracCase; }
    if (hfSubgrid_) setHfTips(); else applyHfTrans();
    stateVer_ = 6;
    if (case_.salt.on && saltTr_ >= 0) {         // C29: re-apply the porosity / permeability factors of the saved solid
        const size_t na0 = size_t(grid_.numActive());
        if (saltS_.size() != na0) { saltS_.assign(na0, 0.0); saltEvap_.assign(na0, 0.0); }
        saltApply();
    }
    evalAll();
    if (!trInfo_.empty()) {                       // C27 tracers: a state without (other) tracers restarts from the case values
        const size_t nt = trInfo_.size(), na0 = size_t(grid_.numActive());
        if (trC_.size() != nt * na0 || trDecayCum_.size() != nt) tracerSetup();
        tracerVolumes();
    } else { trC_.clear(); trDecayCum_.clear(); }
    if (ntOld_.size() != zs_.size()) solventAccept();   // (state version 9 restores the old-time solvent amount)
    chemAccept();                                 // C35: old-time chemical mass of the restored state
    {   // old-time accumulation and convergence reference only: the step-acceptance updates (polymer, hysteresis,
        // VAPPARS, DRSDT histories) already ran before the state was saved and are part of it
        const int na = grid_.numActive();
        for (int a = 0; a < na; ++a)
            for (int p = 0; p < 3; ++p) accOld_[size_t(a) * 3 + p] = (eq_[p] >= 0) ? cf_[a].acc[p].v : 0.0;
        bRef_.assign(3, 1.0);
        double sb[3] = {0, 0, 0}; int n = 0;
        for (int a = 0; a < na; a += std::max(1, na / 2000)) {
            const CellFull C = evalCell(a);
            for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) sb[p] += C.b[p].v;
            ++n;
        }
        for (int p = 0; p < 3; ++p) if (eq_[p] >= 0 && n > 0) bRef_[p] = std::max(1e-12, sb[p] / n);
    }
    // C26: keep the history up to the restored time (a fresh simulator has none; a rewound one - iterative geomechanics
    // coupling repeats a report step from its start state - keeps the summary and snapshots of the accepted steps)
    {
        const double eps = 1e-9 * std::max(1.0, std::fabs(time_)), tk = time_ + eps;
        const bool rewind = time_ > 0 && !summary_.empty() && summary_.back().time >= time_ - eps;   // this run passed the restored time
        if (!rewind) { summary_.clear(); snaps_.clear(); }
        while (!summary_.empty() && summary_.back().time > tk) summary_.pop_back();
        while (!snaps_.empty() && snaps_.back().time > tk) snaps_.pop_back();
    }
    adj_.clear(); adjRec_ = false;
}

} // namespace ressim

// ====================================================================== C18 in-run hydraulic-fracture growth (LEFM)
namespace ressim {
void Simulator::setupHfGrowth() {
    hfOpen_.clear(); hfS_.clear(); hfHalf_.clear(); hfConnsOf_.clear(); hfT0_.clear(); hfEvents_.clear();
    hfTip_.clear(); hfTipMax_.clear(); hfFrac_.clear(); hfHostConnsOf_.clear(); hfSubgrid_ = false; hfCond_.clear(); hfFF_.clear();
    pfLinks_.clear(); pfRec_.clear(); hfPrest_.clear(); hfLeakM_.clear(); hfSg0_.clear();
    hfSubLeak_ = std::getenv("RESSIM_SUBLEAK") && std::atoi(std::getenv("RESSIM_SUBLEAK")) > 0;
    for (const auto& F : case_.grid.hfracs) if (F.subgrid && F.subLeak) hfSubLeak_ = true;
    const GridInput& g = case_.grid;
    bool any = false;
    for (const auto& F : g.hfracs) if (F.kic > 0) any = true;
    const int nH = grid_.numHfCells();
    if (!any || nH == 0 || int(g.hfPieces.size()) != nH) return;
    hfOpen_.assign(size_t(nH), 1); hfS_.assign(size_t(nH), 0.0); hfHalf_.assign(size_t(nH), 0.0);
    // injection point per fracture: pieces whose host carries a well completion
    std::vector<char> wellHost(size_t(grid_.hfBase()), 0);
    for (const auto& w : case_.wells)
        for (const auto& cp : w.completions)
            for (int k = cp.k1; k <= cp.k2; ++k) {
                if (cp.i < 0 || cp.j < 0 || k < 0 || cp.i >= g.nx || cp.j >= g.ny || k >= g.nz) continue;
                const int q = g.index(cp.i, cp.j, k);
                if (q < int(wellHost.size())) wellHost[size_t(q)] = 1;
            }
    for (size_t f = 0; f < g.hfracs.size(); ++f) {
        const auto& F = g.hfracs[f];
        const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0), tx = (F.x1 - F.x0) / std::max(L, 1e-12), ty = (F.y1 - F.y0) / std::max(L, 1e-12);
        double sInj = 0; int nInj = 0;
        for (int k = 0; k < nH; ++k) {
            const auto& P = g.hfPieces[size_t(k)];
            if (P.frac != int(f)) continue;
            const double s = (0.5 * (P.xa + P.xb) - F.x0) * tx + (0.5 * (P.ya + P.yb) - F.y0) * ty;
            hfS_[size_t(k)] = s;
            hfHalf_[size_t(k)] = 0.5 * std::hypot(P.xb - P.xa, P.yb - P.ya);
            if (P.host < int(wellHost.size()) && wellHost[size_t(P.host)]) { sInj += s; ++nInj; }
        }
        sInj = nInj ? sInj / nInj : 0.5 * L;
        for (int k = 0; k < nH; ++k) {
            if (g.hfPieces[size_t(k)].frac != int(f)) continue;
            hfS_[size_t(k)] -= sInj;
            if (F.kic > 0) hfOpen_[size_t(k)] = (std::fabs(hfS_[size_t(k)]) - hfHalf_[size_t(k)] <= F.initLen + 1e-9) ? 1 : 0;
        }
    }
    buildHfConns();
    for (const auto& F : g.hfracs) if (F.kic > 0 && F.subgrid) hfSubgrid_ = true;
    // C28: persistent phase-field models of the PHASEFIELD fractures (threads start at the first report end)
    for (size_t f = 0; f < g.hfracs.size(); ++f) {
        const auto& F = g.hfracs[f];
        if (!F.phaseField || !(F.kic > 0)) continue;
        auto L = std::make_shared<PfLink>();
        PfLink& P = *L;
        P.frac = int(f);
        const Units Uu = Units::of(case_.units);
        P.toM = case_.units == UnitSystem::Field ? 0.3048 : 1.0;
        P.toPa = 101325.0 / Uu.atm;
        double Lmax = 0, aInj = 0, Hc = 0;
        for (int k = 0; k < nH; ++k) {
            if (g.hfPieces[size_t(k)].frac != int(f)) continue;
            const double s = std::fabs(hfS_[size_t(k)]), h = hfHalf_[size_t(k)];
            Lmax = std::max(Lmax, s + h);
            if (s <= h + 1e-9) { aInj = std::max(aInj, h); Hc += grid_.hfCells[size_t(k)].h; }
        }
        P.Hm = (F.zTop >= 0 && F.zBot > F.zTop ? F.zBot - F.zTop : std::max(Hc, 1e-6)) * P.toM;
        const double zMid = F.zTop >= 0 && F.zBot > F.zTop ? 0.5 * (F.zTop + F.zBot) : 0.0;
        (void)zMid;
        P.sigmaCase = F.sigma;
        const double a = std::max(F.initLen, aInj) * P.toM, LmaxM = std::max(Lmax * P.toM, 2 * a);
        phasefield::Input& in = P.in;
        in.crackA = std::max(a, 1e-3);
        in.h = F.pfH > 0 ? F.pfH : std::max(std::min(in.crackA / 8, LmaxM / 60), LmaxM / 400);
        in.coreX = 1.08 * LmaxM / in.crackA; in.coreY = std::max(0.5, 6 * in.h / in.crackA);
        in.stopLength = 0.93;
        in.E = F.pfE > 0 ? F.pfE * 1e9 : F.ddmE > 0 ? F.ddmE * 1e9 : (case_.geoSeq.enabled && !case_.geoSeq.E.empty() && case_.geoSeq.E[0] > 0 ? case_.geoSeq.E[0] : 20e9);
        in.nu = F.pfE > 0 ? F.pfNu : F.ddmE > 0 ? F.ddmNu : (case_.geoSeq.enabled && !case_.geoSeq.E.empty() && case_.geoSeq.E[0] > 0 && !case_.geoSeq.nu.empty() ? case_.geoSeq.nu[0] : 0.25);
        const double Ep = in.E / (1 - in.nu * in.nu), K = F.kic * 1e6;
        in.Gc = K * K / Ep;
        // net-pressure formulation: an isotropic in-plane closure stress sigma adds exactly sigma to the crack pressure of a
        // straight crack (superposition), so the model runs at sigma = 0 and the reported pressure is sigma + p_net. Running
        // with sigma ~ 100 p_net directly would let O(h) errors of the compressed regularised crack swamp the net pressure.
        in.sxx = in.syy = 0;
        in.control = 1; in.steps = 1; in.maxStagger = 200; in.tol = 1e-3;
        // wells completed in the hosts of this fracture's pieces
        std::vector<char> host(size_t(grid_.hfBase()), 0);
        for (int k = 0; k < nH; ++k) if (g.hfPieces[size_t(k)].frac == int(f) && g.hfPieces[size_t(k)].host < int(host.size())) host[size_t(g.hfPieces[size_t(k)].host)] = 1;
        for (int w = 0; w < int(wells_.size()); ++w) {
            const auto& WI = case_.wells[size_t(wells_[size_t(w)].input)];
            bool hit = false;
            for (const auto& cp : WI.completions)
                for (int kk = cp.k1; kk <= cp.k2 && !hit; ++kk)
                    if (cp.i >= 0 && cp.j >= 0 && kk >= 0 && cp.i < g.nx && cp.j < g.ny && kk < g.nz && host[size_t(g.index(cp.i, cp.j, kk))]) hit = true;
            if (hit) P.wells.push_back(w);
        }
        P.sumIdx = summary_.size();
        pfLinks_.push_back(L);
    }
    if (hfSubgrid_) {
        const size_t nf = g.hfracs.size();
        hfTip_.assign(2 * nf, 0.0); hfTipMax_.assign(2 * nf, 0.0);
        hfFrac_.assign(size_t(nH), 0.0);
        for (int k = 0; k < nH; ++k) {
            const int f = g.hfPieces[size_t(k)].frac;
            for (int side : {-1, 1}) {
                const size_t w = size_t(2 * f + (side > 0));
                hfTipMax_[w] = std::max(hfTipMax_[w], side * hfS_[size_t(k)] + hfHalf_[size_t(k)]);
            }
        }
        for (size_t f = 0; f < nf; ++f)
            for (int side : {-1, 1}) hfTip_[2 * f + (side > 0)] = std::min(g.hfracs[f].initLen, hfTipMax_[2 * f + (side > 0)]);
        setHfTips();
    } else applyHfTrans();
    int open = 0; for (auto o : hfOpen_) open += o;
    log_(fmt("Hydraulic-fracture growth (LEFM): %g fracture pieces, ", nH) + fmt("%g open at the start", open) +
         (hfSubgrid_ ? std::string(", sub-grid tips (implicit K_I = K_IC)") : std::string()));
}

void Simulator::buildHfConns() {
    // connections touching each fracture cell: fracture-fracture (scaled when closed / by the contact conductance) and
    // host-fracture (leak-off area of a sub-grid tip piece); original transmissibilities kept in hfT0_
    const int nH = grid_.numHfCells();
    std::vector<int> hfOfActive(size_t(grid_.numActive()), -1);
    for (int k = 0; k < nH; ++k) {
        const int a = grid_.activeOf[size_t(grid_.hfCells[size_t(k)].global)];
        if (a >= 0) hfOfActive[size_t(a)] = k;
    }
    hfConnsOf_.assign(size_t(nH), {}); hfHostConnsOf_.assign(size_t(nH), {}); hfFF_.clear();
    hfT0_.assign(grid_.conns.size(), -1.0);
    for (size_t c = 0; c < grid_.conns.size(); ++c) {
        const auto& cn = grid_.conns[c];
        const int ka = hfOfActive[size_t(cn.a)], kb = hfOfActive[size_t(cn.b)];
        if (ka < 0 && kb < 0) continue;
        hfT0_[c] = cn.trans;
        if (ka >= 0 && kb >= 0) {
            hfConnsOf_[size_t(ka)].push_back(int(c)); hfConnsOf_[size_t(kb)].push_back(int(c));
            hfFF_.push_back({int(c), ka, kb});
        } else hfHostConnsOf_[size_t(ka >= 0 ? ka : kb)].push_back(int(c));
    }
}

void Simulator::setConnectionMultipliers(const std::vector<int>& conns, const std::vector<double>& mult) {
    const size_t nc = grid_.conns.size();
    for (size_t q = 0; q < conns.size() && q < mult.size(); ++q) {
        const int c = conns[q];
        if (c < 0 || size_t(c) >= nc) continue;
        auto it = connT0_.find(c);
        if (it == connT0_.end()) it = connT0_.emplace(c, grid_.conns[size_t(c)].trans).first;
        const double T0 = it->second;
        const double m = std::isfinite(mult[q]) ? std::clamp(mult[q], 1e-8, 1e8) : 1.0;
        grid_.conns[size_t(c)].trans = T0 * m;
    }
}

void Simulator::setPoreVolumeMultipliers(const std::vector<double>& mult) {
    const size_t na = size_t(grid_.numActive());
    if (pvRef0_.size() != na) pvRef0_.assign(grid_.poreVolRef.begin(), grid_.poreVolRef.begin() + std::min(na, grid_.poreVolRef.size()));
    for (size_t a = 0; a < pvRef0_.size(); ++a) {
        const double m = a < mult.size() && std::isfinite(mult[a]) ? std::clamp(mult[a], 0.2, 5.0) : 1.0;
        grid_.poreVolRef[a] = pvRef0_[a] * m;
    }
}

// C25: closure stress of HF piece k of fracture f at depth z (case units): FE minimum principal stress of the host cell, or
// the HFGROW gradient law
double Simulator::hfPieceSigma(int f, size_t k, double z) const {
    const auto& F = case_.grid.hfracs[size_t(f)];
    if (F.feStress && !sMin_.empty() && k < grid_.hfCells.size()) {
        const int host = grid_.hfCells[k].host;
        const int a = host >= 0 && size_t(host) < grid_.activeOf.size() ? grid_.activeOf[size_t(host)] : -1;
        if (a >= 0 && size_t(a) < sMin_.size() && std::isfinite(sMin_[size_t(a)]) && sMin_[size_t(a)] > 0) return sMin_[size_t(a)];
    }
    const double zMid = 0.5 * (F.zTop + F.zBot) > 0 ? 0.5 * (F.zTop + F.zBot) : z;
    return F.sigma + F.sigmaGrad * (z - zMid);
}

void Simulator::setHfConductance(const std::vector<double>& mult) {
    const int nH = grid_.numHfCells();
    if (nH == 0) return;
    if (hfT0_.empty()) buildHfConns();
    hfCond_.clear();
    if (int(mult.size()) == nH) { hfCond_ = mult; for (double& m : hfCond_) m = std::clamp(std::isfinite(m) ? m : 1.0, 1e-8, 1e8); }
    applyHfTrans();
}

double Simulator::hfCellPressure(int k) const {
    if (k < 0 || k >= grid_.numHfCells()) return 0;
    const int a = grid_.activeOf[size_t(grid_.hfCells[size_t(k)].global)];
    return a >= 0 ? p_[size_t(a)] : 0;
}

void Simulator::applyHfTrans() {
    if (hfT0_.empty()) return;
    const GridInput& g = case_.grid;
    for (size_t c = 0; c < grid_.conns.size(); ++c) if (hfT0_[c] >= 0) grid_.conns[c].trans = hfT0_[c];
    if (!hfCond_.empty())   // C21 G1-c: aperture-dependent conductance, harmonic mean of the two pieces
        for (const HfFF& e : hfFF_) {
            const double ma = hfCond_[size_t(e.a)], mb = hfCond_[size_t(e.b)];
            grid_.conns[size_t(e.c)].trans = hfT0_[size_t(e.c)] * 2.0 * ma * mb / (ma + mb);
        }
    if (hfOpen_.empty()) return;
    // a closed piece conducts like the host rock: fracture-fracture connections touching it are scaled by k_host,n / k_f
    // (times closedMult); its host-fracture connection is a matrix transmissibility already and stays (the piece is
    // then just a sliver of host rock). Closed-closed pairs are scaled once.
    for (size_t k = 0; k < hfOpen_.size(); ++k) {
        const auto& P = g.hfPieces[k];
        const bool sub = !hfFrac_.empty() && g.hfracs[size_t(P.frac)].subgrid;
        const double phi = sub ? hfFrac_[k] : (hfOpen_[k] ? 1.0 : 0.0);
        if (phi >= 1.0) continue;
        if (phi > 0.0) {
            // C21 sub-grid tip piece: hydraulically part of the fracture (its conductance is effectively infinite either way),
            // leak-off area = open fraction -> host-fracture transmissibility x phi. Continuous in the tip position.
            if (k < hfHostConnsOf_.size())
                for (int c : hfHostConnsOf_[k]) grid_.conns[size_t(c)].trans = std::min(grid_.conns[size_t(c)].trans, hfT0_[size_t(c)] * phi *
                                                                                         (k < hfLeakM_.size() ? std::max(hfLeakM_[k], 1e-3) : 1.0));
            continue;
        }
        const double kf = std::max(grid_.hfCells[k].kf, 1e-30);
        const double m = std::min(1.0, std::max(P.kn, 0.0) / kf) * g.hfracs[size_t(P.frac)].closedMult;
        for (int c : hfConnsOf_[k]) {
            grid_.conns[size_t(c)].trans = std::min(grid_.conns[size_t(c)].trans, hfT0_[size_t(c)] * m);
        }
    }
}

double Simulator::hfWingLength(int frac, int side) const {
    double L = 0;
    if (hfOpen_.empty()) return 0;
    if (hfSubgrid_ && frac >= 0 && size_t(frac) < case_.grid.hfracs.size() && case_.grid.hfracs[size_t(frac)].subgrid &&
        case_.grid.hfracs[size_t(frac)].kic > 0)
        return hfTip_[size_t(2 * frac + (side > 0))];   // C21 continuous tip
    for (size_t k = 0; k < hfOpen_.size(); ++k)
        if (hfOpen_[k] && case_.grid.hfPieces[k].frac == frac && (side < 0 ? hfS_[k] <= 0 : hfS_[k] >= 0))
            L = std::max(L, std::fabs(hfS_[k]) + hfHalf_[k]);
    return L;
}

void Simulator::setHfTips() {
    const GridInput& g = case_.grid;
    for (size_t k = 0; k < hfFrac_.size(); ++k) {
        const int f = g.hfPieces[k].frac;
        if (!g.hfracs[size_t(f)].subgrid || !(g.hfracs[size_t(f)].kic > 0)) { hfFrac_[k] = hfOpen_[k] ? 1.0 : 0.0; continue; }
        const double s = hfS_[k], h = std::max(hfHalf_[k], 1e-12);
        double phi;
        if (std::fabs(s) - h <= 1e-9) phi = 1.0;                          // the piece at the injection point
        else {
            const int side = s > 0 ? 1 : -1;
            const double L = hfTip_[size_t(2 * f + (side > 0))];
            phi = std::clamp((L - (std::fabs(s) - h)) / (2 * h), 0.0, 1.0);
        }
        hfFrac_[k] = phi;
        hfOpen_[k] = phi > 1e-12 ? 1 : 0;
    }
    // C28 H2b (research, RESSIM_SUBLEAK=1 or HFGROW ... SUBGRID LEAK): sub-grid leak-off. The exposed part phi of the tip host
    // cell is treated as its own storage: everything the cell absorbed since the tip entered it sits in that part, so its
    // pressure is p_sub = p_rest + (p_cell - p_rest) / phi (p_rest = host pressure at first exposure), and the leak-off of the
    // tip piece is phi T (p_f - p_sub) instead of phi T (p_f - p_cell): multiplier m = (p_f - p_sub) / (p_f - p_cell) in [0, 1],
    // lagged (start-of-step pressures). m -> 1 continuously as phi -> 1.
    if (hfSubLeak_) {
        if (hfPrest_.size() != hfFrac_.size()) hfPrest_.assign(hfFrac_.size(), std::nan(""));
        hfLeakM_.assign(hfFrac_.size(), 1.0);
        for (size_t k = 0; k < hfFrac_.size(); ++k) {
            const double phi = hfFrac_[k];
            if (!(phi > 1e-12)) { hfPrest_[k] = std::nan(""); continue; }
            const int ah = grid_.activeOf[size_t(g.hfPieces[k].host)];
            if (ah < 0) continue;
            const double pc = p_[size_t(ah)], pf = hfCellPressure(int(k));
            const bool gasAware = g.hfracs[size_t(g.hfPieces[k].frac)].subLeak == 2 && G_;
            const double sgc = !G_ ? 0.0 : hasGas_[size_t(ah)] == 1 ? x_[size_t(ah)] : hasGas_[size_t(ah)] == 2 ? 1.0 - sw_[size_t(ah)] : 0.0;
            if (hfSg0_.size() != hfFrac_.size()) hfSg0_.assign(hfFrac_.size(), 0.0);
            if (std::isnan(hfPrest_[k])) { hfPrest_[k] = pc; hfSg0_[k] = sgc; }
            if (phi >= 1.0 || pf <= pc + 1e-12) continue;
            // C48 H2e (LEAKC): what the cell absorbed since exposure sits in the exposed part, gas included. The stored volume
            // phi V ct_sub dp_sub = V ct_cell dp_cell, with Sg_sub = Sg0 + (Sg_cell - Sg0) / phi (capped at 1 - Sw), so an
            // injected gas that concentrates in the sub-cell raises its compressibility and lowers p_sub. ct ratio = 1 for a
            // single-phase or uniform-saturation cell (reduces to LEAK).
            double ratio = 1.0;
            if (gasAware) {
                const double swc = iSw_ >= 0 ? sw_[size_t(ah)] : 0.0;
                const double sgs = std::clamp(hfSg0_[k] + (sgc - hfSg0_[k]) / phi, 0.0, std::max(0.0, 1.0 - swc));
                const double cc = cellCt(ah, sgc, pc), cs = cellCt(ah, sgs, pc);
                if (cc > 0 && cs > 0) ratio = cc / cs;
            }
            const double psub = std::min(pf, hfPrest_[k] + (pc - hfPrest_[k]) * ratio / phi);
            hfLeakM_[k] = std::clamp((pf - psub) / (pf - pc), 0.0, 1.0);
        }
    }
    applyHfTrans();
}

double Simulator::cellCt(int c, double sg, double p) const {
    // C48 H2e: c_t = c_r + Sw c_w + So c_o + Sg c_g from the cell's PVT (c_x = d ln b_x / dp, b = 1/B), Rs/Rv frozen
    const FluidModel& F = fl(c);
    const AD3 P = AD3::variable(p, 0);
    const double sw = iSw_ >= 0 ? sw_[size_t(c)] : (W_ && !O_ ? 1.0 : 0.0);
    const double so = O_ ? std::max(0.0, 1.0 - sw - sg) : 0.0;
    double ct = 0;
    { const AD3 m = F.pvMult(P); if (m.v > 0) ct += m.d[0] / m.v; }
    if (W_ && sw > 0) { const AD3 b = F.bw(P); if (b.v > 0) ct += sw * b.d[0] / b.v; }
    const double rs = G_ && live_ && hasGas_[size_t(c)] == 0 ? x_[size_t(c)] : (G_ && live_ ? F.rsSat(P).v : 0.0);
    if (O_ && so > 0) { const AD3 b = F.bo(P, AD3(rs)); if (b.v > 0) ct += so * b.d[0] / b.v; }
    if (G_ && sg > 0) { const AD3 b = vap_ ? F.bg(P, AD3(F.rvSat(P).v)) : F.bg(P); if (b.v > 0) ct += sg * b.d[0] / b.v; }
    return ct;
}

bool Simulator::hfTipKI(int f, int side, double& KI, double& pnet) const {
    // net pressure = open-length weighted mean over the wing (weights phi_k * length_k * h_k). The uniform-pressure LEFM
    // criterion (KGD / PKN) and, unlike the pressure of one piece, continuous in L: a newly entered tip piece starts with
    // zero weight. (C21 H2: the tip-column pressure made K_I jump at every piece boundary and re-quantised the growth.)
    const GridInput& g = case_.grid;
    const auto& F = g.hfracs[size_t(f)];
    const Units U = Units::of(case_.units);
    const double toPa = 101325.0 / U.atm, toM = case_.units == UnitSystem::Field ? 0.3048 : 1.0;
    const double L = hfTip_[size_t(2 * f + (side > 0))];
    double sw = 0, sp = 0, sz = 0, Htip = 0, Hc = 0, ss = 0;
    for (size_t k = 0; k < hfFrac_.size(); ++k) {
        if (g.hfPieces[k].frac != f || hfFrac_[k] <= 0) continue;
        const double sk = side * hfS_[k], h = hfHalf_[k];
        if (sk < -h) continue;                                       // other wing
        const int a = grid_.activeOf[size_t(grid_.hfCells[k].global)];
        if (a < 0) continue;
        const double w = hfFrac_[k] * 2 * h * grid_.hfCells[k].h;
        sw += w; sp += w * p_[size_t(a)]; sz += w * grid_.depth[size_t(a)];
        if (F.feStress) ss += w * hfPieceSigma(f, k, grid_.depth[size_t(a)]);
        if (L > sk - h - 1e-9 && L <= sk + h + 1e-9) Htip += grid_.hfCells[k].h;   // column containing the tip
        if (sk <= h) Hc += grid_.hfCells[k].h;                                     // injection column
    }
    if (!(sw > 0)) return false;
    const double pAvg = sp / sw, zAvg = sz / sw, H = Htip > 0 ? Htip : Hc;
    const double zMid = 0.5 * (F.zTop + F.zBot) > 0 ? 0.5 * (F.zTop + F.zBot) : zAvg;
    const double sigma = F.feStress && !sMin_.empty() ? ss / sw : F.sigma + F.sigmaGrad * (zAvg - zMid);
    pnet = pAvg - sigma;
    const double c = std::min(L, 0.5 * H) * toM;
    KI = pnet > 0 ? pnet * toPa * std::sqrt(3.14159265358979323846 * c) : 0.0;
    return true;
}

// ---------------------------------------------------------------------------------------------------------------------
// C28 in-run phase-field growth. The phase-field model (phasefield::run, volume control) runs in a worker thread and keeps
// its state (displacement, phase field, history) between report steps: at every report end the simulator posts the next
// crack volume, the worker performs that load step (staggered momentum / phase-field iterations) and reports the crack.
bool Simulator::setHfWingLength(int frac, double L) {
    if (hfOpen_.empty() || frac < 0 || size_t(frac) >= case_.grid.hfracs.size()) return false;
    const double before[2] = {hfWingLength(frac, -1), hfWingLength(frac, 1)};
    bool changed = false;
    for (size_t k = 0; k < hfOpen_.size(); ++k)
        if (case_.grid.hfPieces[k].frac == frac && !hfOpen_[k] && std::fabs(hfS_[k]) <= L + 1e-9) { hfOpen_[k] = 1; changed = true; }
    if (!changed) return false;
    if (!hfFrac_.empty()) for (size_t k = 0; k < hfFrac_.size(); ++k) if (case_.grid.hfPieces[k].frac == frac) hfFrac_[k] = hfOpen_[k] ? 1.0 : 0.0;
    applyHfTrans();
    for (int side : {-1, 1}) {
        const double Lw = hfWingLength(frac, side);
        if (Lw > before[side > 0] + 1e-12) { HfEvent ev; ev.t = time_; ev.frac = frac; ev.side = side; ev.length = Lw; hfEvents_.push_back(ev); }
    }
    return true;
}

void Simulator::pfReportStep() {
    const GridInput& g = case_.grid;
    for (auto& lp : pfLinks_) {
        PfLink& P = *lp;
        const auto& F = g.hfracs[size_t(P.frac)];
        // injected volume since the last report end (surface water ~ reservoir volume; injected gas at its surface volume
        // is not counted - water / slickwater fracturing)
        for (; P.sumIdx < summary_.size(); ++P.sumIdx) {
            const SummaryPoint& s = summary_[P.sumIdx];
            for (int w : P.wells) if (size_t(w) < s.wells.size()) P.vinj += s.wells[size_t(w)].waterInj * s.dt * (case_.units == UnitSystem::Field ? 0.158987 : 1.0);
        }
        PfRecord R; R.t = time_; R.frac = P.frac; R.sigma = P.sigmaCase;
        double V = F.pfEta * P.vinj / std::max(P.Hm, 1e-6);
        double pw = 0, ww = 0;
        for (size_t k = 0; k < hfOpen_.size(); ++k)
            if (g.hfPieces[k].frac == P.frac && hfOpen_[k]) { const double w = hfHalf_[k] * grid_.hfCells[k].h; pw += w * hfCellPressure(int(k)); ww += w; }
        R.pFlow = ww > 0 ? pw / ww : 0.0;
        // C29 TWOWAY: pressure-driven volume V = C p_net(flow) (C from the previous phase-field step; irreversible crack:
        // the volume never decreases); the first step uses the eta volume
        //   mass cap: the crack cannot gain more volume within a step than the fluid injected in it (eta = 1), which keeps an
        //   unstably propagating crack (p_flow above the LEFM pressure) volume-limited as under rate-controlled injection
        if (F.pfTwoWay && P.C > 0) {
            const double pnet = (R.pFlow - P.sigmaCase) * P.toPa;
            const double vCap = P.Vlast + std::max(0.0, P.vinj - P.vinjLast) / std::max(P.Hm, 1e-6);
            V = std::max(P.Vlast, std::min(vCap, P.C * std::max(0.0, pnet)));
            R.note = "two-way";
        }
        phasefield::Step S;
        if (!(V > 0)) { R.note = "no injection yet"; R.halfLength = P.in.crackA; pfRec_.push_back(R); continue; }
        // two-way: the volume increment is applied in sub-steps (x1.5 at most each, <= 16) so that the staggered phase-field
        // iterations follow a fast-growing crack (one large volume jump leaves the solver unconverged)
        bool advOk = true;
        if (F.pfTwoWay && P.Vlast > 0 && V > 1.5 * P.Vlast) {
            double v = P.Vlast;
            for (int sub = 0; sub < 15 && v * 1.5 < V && advOk; ++sub) { v *= 1.5; advOk = P.advance(v, S); }
        }
        if (!advOk || !P.advance(V, S)) {
            R.active = false; R.note = P.error; R.halfLength = hfWingLength(P.frac, 1) * P.toM;
            if (pfRec_.empty() || pfRec_.back().active || pfRec_.back().frac != P.frac) log_("  phase field " + F.name + ": " + P.error);
            pfRec_.push_back(R);
            continue;
        }
        R.volume = S.volume; R.pressure = S.p + P.sigmaCase * P.toPa; R.halfLength = S.lengthHalf; R.stagger = S.stagger;
        { std::lock_guard<std::mutex> lk(P.m); R.field = P.field; }
        pfRec_.push_back(R);
        {   // memory cap of the C30 field history: beyond 48 stored fields drop every second older one (the last is kept)
            int nf = 0;
            for (const auto& q : pfRec_) nf += q.field ? 1 : 0;
            if (nf > 48) { int k = 0; for (size_t q = 0; q + 1 < pfRec_.size(); ++q) if (pfRec_[q].field && (k++ % 2)) pfRec_[q].field.reset(); }
        }
        if (F.pfTwoWay && S.p > 0 && S.volume > 0) {
            P.C = S.volume / S.p; P.Vlast = S.volume; P.vinjLast = P.vinj;
            if (!S.converged && pfRec_.back().note == "two-way") pfRec_.back().note = "two-way (unconverged)";
            // piece conductance from the Sneddon aperture of the phase-field crack w(x) = 4 p_net / E' sqrt(a^2 - x^2)
            const double Ep = P.in.E / (1 - P.in.nu * P.in.nu), a = S.lengthHalf;
            const double w0 = 4 * S.p / Ep * a;
            if (!(P.wRef > 0)) P.wRef = w0;
            const int nH = grid_.numHfCells();
            std::vector<double> m = int(hfCond_.size()) == nH ? hfCond_ : std::vector<double>(size_t(nH), 1.0);
            for (int k = 0; k < nH; ++k) {
                if (g.hfPieces[size_t(k)].frac != P.frac) continue;
                const double x = std::fabs(hfS_[size_t(k)]) * P.toM;
                const double wx = x < a ? 4 * S.p / Ep * std::sqrt(a * a - x * x) : 0.0;
                const double r = wx / P.wRef;
                m[size_t(k)] = std::clamp(r * r * r, 1e-3, 1e3);
            }
            setHfConductance(m);
            log_("  phase field " + F.name + fmt(" (two-way): compliance %.4g m2/m/MPa, aperture at the well %.4g mm", P.C * 1e6, w0 * 1e3) +
                 fmt(" (x%.4g of the first step)", w0 / P.wRef));
        }
        const double Lcase = S.lengthHalf / P.toM;
        if (setHfWingLength(P.frac, Lcase))
            log_("  phase field " + F.name + fmt(": crack volume %.4g m2/m, pressure %.4g MPa", S.volume, R.pressure * 1e-6) +
                 fmt(" (flow model %.4g MPa)", R.pFlow * P.toPa * 1e-6) + fmt(", half-length %.4g m -> wings open to ", S.lengthHalf) + fmt("%.4g ", Lcase) +
                 Units::of(case_.units).len);
    }
}

bool Simulator::growFractures() {
    const GridInput& g = case_.grid;
    const Units U = Units::of(case_.units);
    const double toPa = 101325.0 / U.atm, toM = case_.units == UnitSystem::Field ? 0.3048 : 1.0;
    const double pi = 3.14159265358979323846;
    bool changed = false;
    for (size_t f = 0; f < g.hfracs.size(); ++f) {
        const auto& F = g.hfracs[f];
        if (!(F.kic > 0) || F.subgrid || F.phaseField) continue;   // C21: sub-grid tips inside the step; C28: phase field per report step
        for (int side : {-1, 1}) {
            // tip = outermost open piece(s) of this wing (all layers at that distance), next = nearest closed beyond it
            double sTip = -1; int tip = -1;
            for (size_t k = 0; k < hfOpen_.size(); ++k) {
                if (g.hfPieces[k].frac != int(f) || !hfOpen_[k]) continue;
                const double s = side * hfS_[k];
                if (s >= -hfHalf_[k] && s > sTip) { sTip = s; tip = int(k); }
            }
            if (tip < 0) continue;
            double sNext = 1e300;
            for (size_t k = 0; k < hfOpen_.size(); ++k)
                if (g.hfPieces[k].frac == int(f) && !hfOpen_[k] && side * hfS_[k] > sTip && side * hfS_[k] < sNext) sNext = side * hfS_[k];
            if (sNext > 1e299) continue;   // fully open wing
            // fracture pressure at the tip: highest pressure among the tip-distance pieces (the open column)
            double pTip = -1e300, zTip = 0, H = 0; size_t kTip = 0;
            for (size_t k = 0; k < hfOpen_.size(); ++k) {
                if (g.hfPieces[k].frac != int(f) || !hfOpen_[k] || std::fabs(side * hfS_[k] - sTip) > 1e-6 * (1 + std::fabs(sTip))) continue;
                const int a = grid_.activeOf[size_t(grid_.hfCells[k].global)];
                if (a < 0) continue;
                if (p_[size_t(a)] > pTip) { pTip = p_[size_t(a)]; zTip = grid_.depth[size_t(a)]; kTip = k; }
                H += grid_.hfCells[k].h;
            }
            if (pTip < -1e299) continue;
            double sigma = hfPieceSigma(int(f), kTip, zTip);
            double pnet = pTip - sigma;
            if (F.feStress && !sMin_.empty()) {
                // C25: with the FE closure stress the stress varies from layer to layer (stiffness contrasts, poro-elastic
                // stressing): the tip column is driven by its LARGEST net pressure p_k - sigma_min,k, not by the net pressure
                // of the highest-pressure piece (identical for a linear sigma(z) parallel to the fluid gradient)
                for (size_t k = 0; k < hfOpen_.size(); ++k) {
                    if (g.hfPieces[k].frac != int(f) || !hfOpen_[k] || std::fabs(side * hfS_[k] - sTip) > 1e-6 * (1 + std::fabs(sTip))) continue;
                    const int a = grid_.activeOf[size_t(grid_.hfCells[k].global)];
                    if (a < 0) continue;
                    const double sk = hfPieceSigma(int(f), k, grid_.depth[size_t(a)]);
                    if (p_[size_t(a)] - sk > pnet) { pnet = p_[size_t(a)] - sk; sigma = sk; }
                }
            }
            const double Lw = sTip + hfHalf_[size_t(tip)];
            const double c = std::min(Lw, 0.5 * H) * toM;
            const double KI = pnet > 0 ? pnet * toPa * std::sqrt(pi * c) : 0.0;
            if (KI < F.kic * 1e6) continue;
            // open the next column of pieces (all layers at distance sNext)
            for (size_t k = 0; k < hfOpen_.size(); ++k)
                if (g.hfPieces[k].frac == int(f) && !hfOpen_[k] && std::fabs(side * hfS_[k] - sNext) <= 1e-6 * (1 + std::fabs(sNext))) hfOpen_[k] = 1;
            HfEvent ev; ev.t = time_; ev.frac = int(f); ev.side = side; ev.length = hfWingLength(int(f), side); ev.KI = KI; ev.pnet = pnet;
            hfEvents_.push_back(ev);
            log_("  fracture " + F.name + (side < 0 ? " wing - " : " wing + ") +
                 fmt("grows to %.4g ", ev.length) + U.len + fmt(" at t = %.5g d", time_) + fmt(" (K_I %.3g MPa m^0.5", KI * 1e-6) + fmt(" >= K_IC %.3g)", F.kic));
            changed = true;
        }
    }
    if (changed) applyHfTrans();
    return changed;
}
} // namespace ressim
