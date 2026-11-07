#include <map>
#include "sim_runner.h"
#include "ressim/binio.h"
#include "ressim/geomech.h"
#include "ressim/cell_snapshot.h"
#include "ressim/geoseq.h"
#include "ressim/faultmech.h"
#include <functional>

namespace {
// C25: finite-element geomechanics maps of a report snapshot; box(q) = mechanics cell of flow cell q (-1 none), pu = Pa per
// case pressure unit
void addMechMaps(ressim::Snapshot& sn, const ressim::geomech::Solver& g, const std::function<int(size_t)>& box, size_t n, double pu, bool plastic) {
    ressim::geoseq::mechMaps(sn, g, box, n, pu, plastic);   // C30: shared with the post-processing (geoseq.h)
}
}  // namespace
#include "ecl_engine.h"
#include <chrono>
#include <cstdio>
#include "ressim/ecl.h"
#include <algorithm>
#include <cmath>

// wells of a run whose case lists none (decks ResSim cannot model): taken from the W* summary vectors
// field component rates from one well's net component rates (production > 0, injection < 0)
static void addCompRates(ressim::SummaryPoint& p, const std::vector<double>& q) {
    if (p.fcompP.size() < q.size()) { p.fcompP.resize(q.size(), 0.0); p.fcompI.resize(q.size(), 0.0); }
    for (size_t i = 0; i < q.size(); ++i) { if (q[i] > 0) p.fcompP[i] += q[i]; else p.fcompI[i] -= q[i]; }
}

static void wellsFromSummary(const ressim::ecl::Summary& s, std::vector<std::string>& names, std::vector<ressim::WellType>& types) {
    if (!names.empty()) return;
    std::map<std::string, ressim::WellType> m;
    for (size_t v = 0; v < s.vars.size(); ++v) {
        const auto& x = s.vars[v];
        if (x.category != 'W' || x.wgname.empty() || x.wgname[0] == ':') continue;
        auto it = m.emplace(x.wgname, ressim::WellType::Producer).first;
        const bool wi = x.keyword == "WWIR" || x.keyword == "WWIT", gi = x.keyword == "WGIR" || x.keyword == "WGIT";
        if ((wi || gi) && v < s.values.size()) {
            float mx = 0; for (float f : s.values[v]) mx = std::max(mx, f);
            if (mx > 0) it->second = wi ? ressim::WellType::WaterInjector : ressim::WellType::GasInjector;
        }
    }
    for (auto& [n, t] : m) { names.push_back(n); types.push_back(t); }
}

// C14: base results up to the restart step (caller holds m_)
void SimRunner::prefillFromBase(const RestartRequest& rq) {
    if (!rq.base) throw std::runtime_error("restart: base results missing");
    const RunResults& B = *rq.base;
    const size_t ns = std::min(B.snaps.size(), size_t(rq.step) + 1);
    if (ns == 0 || size_t(rq.step) >= B.snaps.size()) throw std::runtime_error("restart: the base results have no report step " + std::to_string(rq.step));
    const double tEnd = B.snaps[ns - 1].time;
    r_.snaps.assign(B.snaps.begin(), B.snaps.begin() + std::ptrdiff_t(ns));
    r_.summary.clear();
    for (const auto& p : B.summary) if (p.time <= tEnd + 1e-9 * std::max(1.0, std::fabs(tEnd))) r_.summary.push_back(p);
    r_.initial = B.initial;
    r_.components = B.components; r_.trackName = B.trackName;
    r_.fracNames = B.fracNames; r_.fracSegs = B.fracSegs; r_.wellFracLinks = B.wellFracLinks; r_.heatInPlace = B.heatInPlace; r_.tInjRef = B.tInjRef;
    auto cut = [ns](auto& dst, const auto& src) { dst.assign(src.begin(), src.begin() + std::ptrdiff_t(std::min(ns, src.size()))); };
    cut(r_.fracSnap, B.fracSnap); cut(r_.wellT, B.wellT); cut(r_.wellMW, B.wellMW); cut(r_.cooledFrac, B.cooledFrac);
    r_.restarts.clear();
    for (const auto& R : B.restarts) if (R.report <= rq.step) r_.restarts.push_back(R);
    r_.restartedFrom = rq.step;
    r_.valid = true;
}

void SimRunner::start(const ressim::Case& c) {
    if (running_) return;
    join();
    cancel_ = false; running_ = true; progress_ = 0; wall_ = 0;
    {
        std::lock_guard<std::mutex> l(m_);
        error_.clear();
        r_ = RunResults{};
    }
    rqActive_ = std::move(rq_); rq_.reset();
    ++version_;
    th_ = std::thread(&SimRunner::work, this, c);
}

void SimRunner::work(ressim::Case c) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    try {
        auto sim = std::make_shared<ressim::Simulator>();
        sim->log = [this](const std::string& s) { std::lock_guard<std::mutex> l(m_); log_.push_back(s); };
        sim->cancel = &cancel_;
        addLog("=== Run started: " + c.title);
        if (!c.extraKeywords.empty()) {
            std::string names; int shown = 0;
            for (auto& e : c.extraKeywords) {
                if (e.section == ressim::DeckSection::Summary || names.find(e.name) != std::string::npos ||
                    e.note.rfind("minimum dimensions", 0) == 0) continue;
                // C27: simulated by the ResSim engine from the verbatim records (economic limits, well tests, tracers)
                if (e.name == "WECON" || e.name == "WTEST" || e.name == "TRACER" || e.name == "WTRACER" || e.name.rfind("TBLK", 0) == 0) continue;
                if (shown++ < 12) names += (names.empty() ? "" : " ") + e.name;
            }
            if (!names.empty())
                addLog("NOTE: advanced keywords used only by OPM Flow are not simulated by the ResSim engine: " + names +
                       (shown > 12 ? " ..." : "") + " (switch the engine to OPM Flow to include them)");
        }
        if (!c.grid.disperc.empty()) addLog("NOTE: DISPERC (mechanical dispersion) is modelled by the compositional / CO2 engine only - ignored by the black-oil engine");
        sim->setup(c);
        const double total = sim->totalSimTime();
        sim->onTimestep = [this, total](const ressim::SummaryPoint& s) {
            std::lock_guard<std::mutex> l(m_);
            r_.summary.push_back(s);
            progress_ = total > 0 ? s.time / total : 0;
            ++version_;
        };
        sim->initialize();
        if (!sim->tracers().empty() || sim->numChem() > 0) {   // C27 passive tracers (+ C35 chemicals) ride on the component vectors (kind 4, amount units, mw 1)
            std::vector<ressim::ComponentInfo> ci;
            for (const auto& T : sim->tracers()) { ressim::ComponentInfo q; q.name = T.name; q.mw = 1.0; q.kind = 4; q.species = T.phase; ci.push_back(q); }
            for (size_t k = 0; k < sim->numChem(); ++k) { ressim::ComponentInfo q; q.name = sim->chemName(k); q.mw = 1.0; q.kind = 4; q.species = 0; ci.push_back(q); }
            std::lock_guard<std::mutex> l(m_);
            r_.components = std::move(ci);
        }
        const std::shared_ptr<RestartRequest> rq = rqActive_;
        // ---- C24: coupled geomechanics of the main run. Fault mechanics (STRESSFIELD + FAULTMECH): Coulomb stress,
        // slip tendency, Dieterich seismicity and the reactivation -> fault transmissibility feedback on every FAULTS
        // face after every report step; G1-c (GEOCOUPLE): EFEM contact solve of the hydraulic fractures -> fracture
        // transmissibilities. Both act between report steps (sequential, one report step lag).
        std::unique_ptr<ressim::faultmech::Monitor> fm;
        std::unique_ptr<ressim::geocouple::Output> gco;
        ressim::geocouple::Input gin; ressim::efem::Result gprev;
        const double paU = c.units == ressim::UnitSystem::Field ? 6894.757 : 1e5;
        auto geoStep = [&]() {
            if (!gco) return;
            ressim::geocouple::StepOut st; std::string err;
            if (!ressim::geocouple::coupleOnce(*sim, gin, gprev, st, gco.get(), err)) { addLog("NOTE: G1-c contact coupling stopped: " + err); gco->error = err; gco.reset(); return; }
            char b[256];
            std::snprintf(b, sizeof b, "  G1-c t = %.4g d: %d open / %d slipping pieces, max opening %.3g mm, max slip %.3g mm, conductance x %.3g .. %.3g",
                          st.t, st.openPieces, st.slipPieces, st.maxOpening * 1e3, st.maxSlip * 1e3, st.minMult, st.maxMult);
            addLog(b);
            gco->steps.push_back(std::move(st)); gco->ok = true;
        };
        std::function<void()> connStepFn;   // C25: composed fault x geomechanics transmissibility update (set below)
        auto faultStep = [&]() {
            if (fm) fm->step(sim->time(), sim->snapshots().back().pressure);
            if (connStepFn) connStepFn();
            else if (fm) fm->applyCoupling(*sim);
        };
        if (c.faultMech.on()) {
            fm = std::make_unique<ressim::faultmech::Monitor>();
            std::string err;
            if (!fm->setup(c, sim->grid(), err)) { addLog("NOTE: fault mechanics off - " + err); fm.reset(); }
            else {
                char b[200]; std::snprintf(b, sizeof b, "=== Fault mechanics: %zu fault faces on %zu fault(s), %s coupling", fm->out().faces.size(), fm->out().names.size(), fm->twoWay() ? "two-way (reactivation -> transmissibility)" : "one-way");
                addLog(b);
            }
        }
        if (c.geoCouple.on) {
            if (c.grid.hfracs.empty()) addLog("NOTE: GEOCOUPLE needs hydraulic fractures (HFRAC) - contact coupling off");
            else {
                gin.layer = std::max(0, c.geoCouple.layer - 1); gin.E = c.geoCouple.eGPa * 1e9; gin.nu = c.geoCouple.nu;
                gin.sxx = -c.geoCouple.sxx * paU; gin.syy = -c.geoCouple.syy * paU; gin.friction = c.geoCouple.friction;
                gco = std::make_unique<ressim::geocouple::Output>();
                addLog("=== G1-c fracture contact coupling in the run (EFEM slice of layer " + std::to_string(c.geoCouple.layer) + ")");
            }
        }
        // ---- C26: DDM mechanical apertures of the in-run hydraulic fractures (HFGROW ... DDMWIDTH)
        std::unique_ptr<ressim::hfddm::Aperture> hfd;
        {
            bool want = false;
            for (const auto& F : c.grid.hfracs) if (F.ddmWidth) want = true;
            if (want) {
                hfd = std::make_unique<ressim::hfddm::Aperture>();
                std::string err;
                if (!hfd->setup(*sim, err)) { addLog("NOTE: DDMWIDTH apertures off - " + err); hfd.reset(); }
                else {
                    addLog("=== Hydraulic-fracture apertures from a 3-D displacement-discontinuity solve (DDMWIDTH): " + std::to_string(hfd->out().elements) +
                           " pieces, conductance from the cubic law after every report step");
                    if (gco) addLog("NOTE: GEOCOUPLE (G1-c) and DDMWIDTH both set fracture conductances - the DDM multipliers are applied last on DDMWIDTH pieces");
                }
            }
        }
        auto ddmStep = [&]() {
            if (!hfd) return;
            hfd->update(*sim);
            const auto& o = hfd->out();
            char b[220];
            std::snprintf(b, sizeof b, "  DDM apertures t = %.4g d: %d / %d pieces open, max aperture %.3g mm, opening volume %.4g m3", sim->time(), o.open, o.elements, o.maxWidth * 1e3, o.volume);
            addLog(b);
        };
        // ---- C25: sequential FE geomechanics (GEOMECH of a black-oil deck) on the Cartesian / corner-point grid
        std::unique_ptr<ressim::geoseq::Coupler> gs;
        bool anyFeHf = false;
        for (const auto& F : c.grid.hfracs) if (F.kic > 0 && F.feStress) anyFeHf = true;
        const bool geoPerm = (c.geoSeq.plasticity == 1 && c.geoSeq.permPlastic > 0) || c.geoSeq.permStress > 0;
        if (c.geoSeq.enabled) {
            const ressim::Grid& G = sim->grid();
            const size_t na = size_t(G.numActive());
            const auto& P0 = sim->snapshots().back().pressure;
            std::vector<double> pPa(na, 0.0), phi(na, 0.2);
            for (size_t a = 0; a < na; ++a) {
                pPa[a] = a < P0.size() ? double(P0[a]) * paU : 0.0;
                const int gq = G.globalOf[a];
                if (gq >= 0 && size_t(gq) < G.bulkVol.size() && G.bulkVol[size_t(gq)] > 0) phi[a] = std::clamp(G.poreVolRef[a] / G.bulkVol[size_t(gq)], 0.01, 0.6);
            }
            gs = std::make_unique<ressim::geoseq::Coupler>();
            std::string err;
            if (!gs->setup(c.geoSeq, G, G.globalOf, c.units == ressim::UnitSystem::Field ? 0.3048 : 1.0, pPa, phi, err)) { addLog("NOTE: FE geomechanics off - " + err); gs.reset(); }
            else {
                char b[300];
                std::snprintf(b, sizeof b, "=== FE geomechanics (sequential): %d x %d x %d box, %s, %s coupling%s%s", G.nx, G.ny, G.nz,
                              gs->hexMesh() ? "isoparametric hexahedra on the corner-point nodes" : "Cartesian bricks",
                              gs->twoWay() ? "two-way (pore volume + permeability after every report step)" : "one-way",
                              c.geoSeq.plasticity == 1 ? ", Drucker-Prager plasticity" : "", anyFeHf ? ", HFGROW closure stress from the FE stress" : "");
                addLog(b);
            }
        } else if (anyFeHf) addLog("NOTE: HFGROW ... FESTRESS needs GEOMECH - the HFGROW stress gradient is used");
        if (fm && fm->wantsFeStress()) {
            if (gs) { ressim::geoseq::Coupler* gp = gs.get(); fm->setStressProvider([gp](int cell, double S[6]) { return gp->totalStress(cell, S); }); }
            else addLog("NOTE: STRESSFIELD ... FE needs GEOMECH - the analytic stress field is used");
        }
        std::vector<int> geoTouched;   // connections carrying a composed (fault x geomechanics) multiplier
        auto connStep = [&]() {
            // C25: one absolute multiplier per connection = fault reactivation x geomechanics permeability (setConnectionMultipliers
            // is last-call-wins per connection, so the two are composed here); hydraulic-fracture connections are left to the
            // fracture model (applyHfTrans owns them)
            std::map<int, double> M;
            const ressim::Grid& G = sim->grid();
            if (gs && gs->twoWay() && geoPerm) {
                std::vector<std::pair<int, int>> pr; std::vector<int> idx;
                for (size_t q = 0; q < G.conns.size(); ++q) {
                    const auto& cn = G.conns[q];
                    if (G.globalOf[size_t(cn.a)] >= G.numCells() || G.globalOf[size_t(cn.b)] >= G.numCells()) continue;
                    pr.push_back({cn.a, cn.b}); idx.push_back(int(q));
                }
                const std::vector<double> f = gs->connFactor(pr);
                for (size_t q = 0; q < f.size(); ++q) if (f[q] != 1.0) M[idx[q]] = f[q];
            }
            if (fm) {
                const std::map<int, double> F = fm->connMultipliers();
                for (const auto& e : F) { auto it = M.find(e.first); if (it == M.end()) M[e.first] = e.second; else it->second *= e.second; }
                if (!F.empty()) ++fm->out().twoWayUpdates;
            }
            for (int q : geoTouched) if (!M.count(q)) M[q] = 1.0;
            if (M.empty()) return;
            std::vector<int> cs; std::vector<double> ms;
            for (const auto& e : M) { cs.push_back(e.first); ms.push_back(e.second); }
            sim->setConnectionMultipliers(cs, ms);
            geoTouched.clear();
            for (const auto& e : M) if (e.second != 1.0) geoTouched.push_back(e.first);
        };
        // C26 (C25-4a): iterative two-way coupling (GEOMECH coupling 2) - the report step is repeated from its saved start
        // state with the pore-volume multipliers of the new mechanics state until they change by less than tol (fixed point
        // of flow <-> mechanics per report step, at most maxIter passes); the explicit scheme lags them one report step
        const bool gsIter = gs && c.geoSeq.coupling == 2;
        std::vector<double> pvApplied;          // multipliers in the flow model now (iterative scheme)
        bool gsStepped = false;                 // the iterative loop already ran the FE solve of this report step
        int gsIterTotal = 0, gsIterSteps = 0, gsIterMaxed = 0;
        auto mechStep = [&]() {
            if (!gs) return;
            const auto& P = sim->snapshots().back().pressure;
            std::vector<double> pPa(P.size());
            for (size_t a = 0; a < P.size(); ++a) pPa[a] = double(P[a]) * paU;
            if (gsStepped) {   // iterative: FE state and pore volume already consistent - only the products
                gsStepped = false;
                if (anyFeHf) { std::vector<double> sm = gs->minStress(); for (double& v : sm) v /= paU; sim->setMinStress(sm); }
                return;
            }
            if (!gs->step(pPa)) {
                addLog("NOTE: FE geomechanics: the mechanics solve did not converge - coupling stopped");
                if (fm) fm->setStressProvider({});
                gs.reset(); return;
            }
            double mlo = 1, mhi = 1;
            if (gs->twoWay()) {
                const std::vector<double> m = gs->pvFactor(pPa, c.props.rock.cr / paU);
                for (double v : m) { mlo = std::min(mlo, v); mhi = std::max(mhi, v); }
                sim->setPoreVolumeMultipliers(m);
                pvApplied = m;
            }
            if (anyFeHf) { std::vector<double> sm = gs->minStress(); for (double& v : sm) v /= paU; sim->setMinStress(sm); }
            const auto& g = gs->solver();
            char b[256];
            std::snprintf(b, sizeof b, "  geomechanics t = %.4g d: surface uplift %.3g mm, max failure index %.3f%s, PV factor %.4f .. %.4f, %d CG its",
                          sim->time(), g.surfaceUplift() * 1e3, gs->maxFailure(), g.plasticPoints() ? (", " + std::to_string(g.plasticPoints()) + " yielding Gauss points").c_str() : "",
                          mlo, mhi, g.lastCg());
            addLog(b);
        };
        auto mechMaps = [&](ressim::Snapshot& sn) {
            if (!gs) return;
            const ressim::geoseq::Coupler* gp = gs.get();
            addMechMaps(sn, gp->solver(), [gp](size_t q) { return gp->boxOf(int(q)); }, sn.pressure.size(), paU, c.geoSeq.plasticity == 1);
        };
        auto mechSummary = [&]() {
            if (!gs) return;
            std::lock_guard<std::mutex> l(m_);
            if (r_.summary.empty()) return;
            r_.summary.back().fuplift = gs->solver().surfaceUplift(); r_.summary.back().fupliftRes = gs->solver().reservoirTopUplift();
            r_.summary.back().fmaxFail = gs->maxFailure();
        };
        if (rq) {
            if (rq->kind != "BO") throw std::runtime_error("restart: the saved state is from the " + rq->kind + " engine, not black oil");
            sim->loadState(rq->state);
            char rb[160]; std::snprintf(rb, sizeof rb, "=== Restart from report step %d (time %.6g) - continuing the base run exactly", rq->step, sim->time());
            addLog(rb);
        }
        {
            std::lock_guard<std::mutex> l(m_);
            r_.grid = std::make_shared<ressim::Grid>(sim->grid());
            r_.input = c;
            if (rq) prefillFromBase(*rq);
            else {
                r_.summary = sim->summary();
                r_.snaps = sim->snapshots();
                r_.initial = sim->inPlace();
            }
            r_.totalTime = total;
            for (auto& w : c.wells) { r_.wellNames.push_back(w.name); r_.wellTypes.push_back(w.type); }
            r_.valid = true;
            ++version_;
        }
        if (rq && (fm || gco || gs || hfd)) addLog("NOTE: restarted run - geomechanics / fault mechanics / contact coupling start from the restart state");
        connStepFn = connStep;
        mechStep(); faultStep(); geoStep(); ddmStep(); mechSummary();
        if (gs) { std::lock_guard<std::mutex> l(m_); if (!r_.snaps.empty()) mechMaps(r_.snaps.back()); }
        auto publish = [&]() {
            std::shared_ptr<const ressim::faultmech::Output> a; std::shared_ptr<const ressim::geocouple::Output> b;
            std::shared_ptr<const ressim::hfddm::Output> d;
            if (fm) a = std::make_shared<const ressim::faultmech::Output>(fm->out());
            if (gco) b = std::make_shared<const ressim::geocouple::Output>(*gco);
            if (hfd) d = std::make_shared<const ressim::hfddm::Output>(hfd->out());
            std::lock_guard<std::mutex> l(m_);
            if (a) r_.faultMech = a;
            if (b) r_.geoCouple = b;
            if (d) r_.hfDdm = d;
            if (!sim->phaseFieldRecords().empty()) r_.hfPf = std::make_shared<const std::vector<ressim::Simulator::PfRecord>>(sim->phaseFieldRecords());
            r_.actions = sim->actionEvents();
            if (!sim->mswRecords().empty()) { r_.msw = sim->mswWells(); r_.mswRec = std::make_shared<const std::vector<ressim::Simulator::MswRecord>>(sim->mswRecords()); }
            if (!sim->networkHistory().empty()) { r_.netNodes = sim->networkNodes(); r_.netRec = std::make_shared<const std::vector<ressim::Simulator::NetRecord>>(sim->networkHistory()); }
            ++version_;
        };
        publish();
        // one report step; iterative geomechanics: repeat it until the pore-volume multipliers agree with the mechanics
        auto reportStep = [&]() -> bool {
            if (!gsIter || !gs) return sim->runReportStep();
            size_t nSum0 = 0; { std::lock_guard<std::mutex> l(m_); nSum0 = r_.summary.size(); }
            const ressim::geoseq::IterStep R = ressim::geoseq::iterativeReportStep(*sim, *gs, pvApplied, paU, c.props.rock.cr / paU, cancel_ ? 1 : c.geoSeq.maxIter, c.geoSeq.tol,
                [&] { std::lock_guard<std::mutex> l(m_); if (r_.summary.size() > nSum0) r_.summary.resize(nSum0); });
            if (!R.ran) return false;
            if (!R.feOk) return true;   // mechStep repeats the FE solve, reports the failure and stops the coupling
            gsStepped = true;
            gsIterTotal += R.passes; ++gsIterSteps; if (!R.converged) ++gsIterMaxed;
            char b[200]; std::snprintf(b, sizeof b, "  iterative geomechanics t = %.4g d: %d pass%s, max pore-volume change %.2e%s", sim->time(), R.passes, R.passes > 1 ? "es" : "", R.diff,
                                       R.converged ? "" : " (not converged - accepted after maxIter)");
            addLog(b);
            return true;
        };
        while (!cancel_ && reportStep()) {
            mechStep(); faultStep(); geoStep(); ddmStep(); mechSummary();
            publish();
            std::lock_guard<std::mutex> l(m_);
            r_.snaps.push_back(sim->snapshots().back());
            mechMaps(r_.snaps.back());
            r_.stats = sim->stats();
            if (wantsRestart(sim->reportStep(), sim->finished())) {   // C14: engine state for exact restarts
                ressim::ResultsFile::Restart R;
                R.report = sim->reportStep(); R.time = sim->time(); R.kind = "BO"; R.state = sim->saveState();
                R.summaryCount = r_.summary.size(); R.snap = int(r_.snaps.size()) - 1;
                r_.restarts.push_back(std::move(R));
            }
            wall_ = std::chrono::duration<double>(clock::now() - t0).count();
            ++version_;
        }
        if (gsIterSteps > 0) {
            char b[200]; std::snprintf(b, sizeof b, "=== Iterative geomechanics: %d report steps, %.2f flow passes per step on average, %d accepted at maxIter", gsIterSteps,
                                       double(gsIterTotal) / gsIterSteps, gsIterMaxed);
            addLog(b);
        }
        if (fm && !cancel_) {
            fm->finish(&cancel_, [this](const std::string& s) { addLog(s); });
            const auto& O = fm->out();
            for (size_t f = 0; f < O.series.size(); ++f) if (!O.series[f].t.empty()) {
                const auto& S = O.series[f];
                char b[256]; std::snprintf(b, sizeof b, "    fault %s: max dCFS %.3g, max slip tendency %.3f, %.0f %% of the faces reactivated at the end, expected events %.4g, catalogue %.0f (Mw max %.2f, cap %.2f)",
                                           O.names[f].c_str(), *std::max_element(S.maxDcfs.begin(), S.maxDcfs.end()), *std::max_element(S.maxSt.begin(), S.maxSt.end()),
                                           100 * S.reactFrac.back(), S.expected.back(), S.events.back(), S.maxMw.back(), f < O.mMax.size() ? O.mMax[f] : 0.0);
                addLog(b);
            }
            publish();
        }
        const auto& st = sim->stats();
        char buf[256];
        std::snprintf(buf, sizeof(buf), "=== %s: %d timesteps, %d cuts, %d Newton, %d linear iterations, %.2f s",
                      cancel_ ? "Run stopped" : "Run finished", st.timesteps, st.cuts, st.newton, st.linear,
                      std::chrono::duration<double>(clock::now() - t0).count());
        addLog(buf);
        if (st.gpuSolves > 0) {
            std::snprintf(buf, sizeof(buf), "    GPU linear solver: %d device solves (%d CPU fallbacks), %.2f s wall in device solves, GPU busy %s", st.gpuSolves,
                          st.gpuFallbacks, st.gpuWallSec, st.gpuDeviceSec >= 0 ? (std::to_string(st.gpuDeviceSec).substr(0, 6) + " s (timestamp queries)").c_str() : "n/a");
            addLog(buf);
        }
        std::lock_guard<std::mutex> l(m_);
        r_.stats = st;
        sim->onTimestep = nullptr; sim->log = nullptr; sim->cancel = nullptr;
        r_.sim = sim;
    } catch (const std::exception& e) {
        std::lock_guard<std::mutex> l(m_);
        error_ = e.what();
        log_.push_back(std::string("ERROR: ") + e.what());
    }
    wall_ = std::chrono::duration<double>(clock::now() - t0).count();
    running_ = false;
    ++version_;
}

// ================================================================================================ compositional core (E300 decks)
#include "ressim/comp_deck.h"
#include "ressim/deck_parser.h"

void SimRunner::startComp(const ressim::Case& c, const std::string& deckPath) {
    if (running_) return;
    join();
    cancel_ = false; running_ = true; progress_ = 0; wall_ = 0;
    {
        std::lock_guard<std::mutex> l(m_);
        error_.clear();
        r_ = RunResults{};
    }
    rqActive_ = std::move(rq_); rq_.reset();
    ++version_;
    th_ = std::thread(&SimRunner::workComp, this, c, deckPath);
}

void SimRunner::workComp(ressim::Case c, std::string deckPath) {
    using clock = std::chrono::steady_clock;
    using namespace ressim;
    const auto t0 = clock::now();
    try {
        const bool co2 = c.comp.co2store;
        addLog("=== " + std::string(co2 ? "CO2 storage" : "Compositional") + " run started: " + c.title + (co2 ? " (ResSim CO2STORE core: Spycher-Pruess / Duan-Sun CO2-brine, FIM)" : c.comp.fluid.co2sol ? " (ResSim EOS core, FIM, CO2SOL: CO2 dissolves in the brine)" : " (ResSim EOS core, FIM)"));
        // fluid, separator, injection streams, grid, wells and schedule: all from the (edited) case
        comp::DeckModel m; std::string err;
        if (!c.comp.enabled && !deckPath.empty()) {   // legacy: fluid from the original deck
            const deck::RawDeck D = deck::parseDeckFile(deckPath);
            if (!comp::deckModel(D, c, m, &err)) throw std::runtime_error("compositional model: " + err);
        } else if (!comp::deckModel(c, m, &err)) throw std::runtime_error("compositional model: " + err);
        char buf[320];
        std::snprintf(buf, sizeof(buf), "Compositional model: %d cells, %zu connections, %d components (%s EOS), %zu wells, %zu report steps, T %.1f K",
                      m.base.size(), m.base.conns.size(), m.base.fluid.nc(), m.base.fluid.type == eos::EosType::PR ? "Peng-Robinson" : "SRK",
                      m.wellNames.size(), m.stepDays.size(), m.base.T);
        addLog(buf);
        if (co2) {
            const co2::Fluid cf = co2::makeFluid(m.base.T, m.base.salinity);
            addLog("CO2-brine fluid: " + co2::describe(cf));
        }
        for (const auto& n : m.notes) addLog("NOTE: " + n);
        double total = 0; for (double d : m.stepDays) total += d;
        auto grid = std::make_shared<Grid>(); grid->build(c);
        {
            std::lock_guard<std::mutex> l(m_);
            r_.grid = grid; r_.input = c; r_.totalTime = total; r_.engine = co2 ? "ResSim CO2STORE" : "ResSim compositional";
            for (auto& w : c.wells) { r_.wellNames.push_back(w.name); r_.wellTypes.push_back(w.type); }
        }
        const double pu = c.units == UnitSystem::Field ? 6894.757293168 : 1e5;
        size_t step = 0;
        std::vector<comp::CompSummary> out; comp::CompStats st;
        double lastT = 0;
        // C14 restart: save / continue
        const std::shared_ptr<RestartRequest> rq = rqActive_;
        DeckRestart dr; dr.saveEvery = rstEvery_; dr.saveLast = rstLast_;
        dr.onSave = [&](int stp, double days, const std::string& blob) {
            std::lock_guard<std::mutex> l(m_);
            ressim::ResultsFile::Restart R; R.report = stp; R.time = days; R.kind = "COMP"; R.state = blob;
            R.summaryCount = r_.summary.size() + 1; R.snap = int(r_.snaps.size());
            r_.restarts.push_back(std::move(R));
        };
        if (rq) {
            if (rq->kind != "COMP") throw std::runtime_error("restart: the saved state is from the " + rq->kind + " engine, not compositional");
            dr.startStep = rq->step; dr.state = rq->state;
            std::lock_guard<std::mutex> l(m_);
            prefillFromBase(*rq);
            step = size_t(rq->step) + 1;
            lastT = r_.summary.empty() ? 0.0 : r_.summary.back().time;
            addLog("=== Restart from report step " + std::to_string(rq->step) + " - continuing the base run exactly");
        }
        auto toPoint = [&](const comp::CompSummary& S) {
            SummaryPoint p;
            p.time = S.time; p.dt = S.time - lastT; lastT = S.time;
            p.fopr = S.fopr; p.fgpr = S.fgpr; p.fwpr = S.fwpr; p.fgir = S.fgir; p.fwir = S.fwir;
            p.fopt = S.fopt; p.fgpt = S.fgpt; p.fwpt = S.fwpt; p.fgit = S.fgit; p.fwit = S.fwit;
            p.fpr = S.fpr; p.fprp = S.fpr; p.oip = S.foip; p.gip = S.fgip; p.wip = S.fwip;
            p.fwct = S.fopr + S.fwpr > 0 ? S.fwpr / (S.fopr + S.fwpr) : 0.0;
            p.fgor = S.fopr > 0 ? S.fgpr / S.fopr : 0.0;
            p.co2Mob = S.co2Mob; p.co2Trap = S.co2Immob; p.co2Diss = S.co2Diss;
            p.co2Struct = S.co2Struct; p.co2Min = S.co2Min; p.co2KrImmob = S.co2KrImmob; p.saltSolid = S.saltSolid;
            const std::vector<comp::CompWell>* ws = step > 0 && step - 1 < m.stepWells.size() ? &m.stepWells[step - 1] : nullptr;
            for (size_t w = 0; w < S.wells.size(); ++w) {
                WellRates r; const auto& W = S.wells[w];
                r.oil = W.oil; r.gas = W.gas; r.water = W.water; r.waterInj = W.winj; r.gasInj = W.ginj; r.bhp = W.bhp; r.open = W.open;
                r.mode = Control::BHP;
                if (ws && w < ws->size() && !W.onBhp) {
                    const comp::CompWell& cw = (*ws)[w];
                    if (cw.injector) r.mode = Control::RATE;
                    else switch (cw.control) { case 2: r.mode = Control::RESV; break; case 3: r.mode = Control::ORAT; break; case 4: r.mode = Control::GRAT; break;
                                               case 5: r.mode = Control::WRAT; break; case 6: r.mode = Control::LRAT; break; default: break; }
                }
                r.comp = W.molar;
                addCompRates(p, r.comp);
                p.wells.push_back(r);
            }
            return p;
        };
        // ---- C25 coupled geomechanics of the compositional / CO2 run (sequential, after every report step): FE geomechanics
        // (GEOMECH), fault mechanics + seismicity (STRESSFIELD / FAULTMECH) and report-step LEFM growth of HFRAC (HFGROW)
        std::unique_ptr<geoseq::Coupler> gsC; std::unique_ptr<faultmech::Monitor> fmC; std::unique_ptr<geoseq::HfSequential> hfC;
        std::shared_ptr<Grid> gxC;   // grid of the simulator cells (LGR / MINC / HF expanded as in comp::deckModel)
        const bool geoPermC = (c.geoSeq.plasticity == 1 && c.geoSeq.permPlastic > 0) || c.geoSeq.permStress > 0;
        bool anyFeHfC = false, anyHfGrow = false;
        for (const auto& F : c.grid.hfracs) { if (F.kic > 0) anyHfGrow = true; if (F.kic > 0 && F.feStress) anyFeHfC = true; }
        const bool wantGeo = c.geoSeq.enabled || c.faultMech.on() || anyHfGrow;
        auto coupleComp = [&](comp::CompSimulator& sim, double tDays) {
            if (!wantGeo || rq) return;
            const size_t ns = sim.pressure().size();
            if (step == 0) {
                Case cx = c;
                try { expandLgrArrays(cx); expandMincArrays(cx); applyDpnum(cx); prepareHydraulicFractures(cx); gxC = std::make_shared<Grid>(); gxC->build(cx); }
                catch (const std::exception& e) { addLog(std::string("NOTE: coupled geomechanics off - grid: ") + e.what()); gxC.reset(); }
                if (gxC && size_t(gxC->numActive()) != ns) { addLog("NOTE: coupled geomechanics off - simulator cells do not match the grid"); gxC.reset(); }
                if (!gxC) return;
                const Grid& G = *gxC;
                const double lm = c.units == UnitSystem::Field ? 0.3048 : 1.0;
                if (c.geoSeq.enabled) {
                    std::vector<double> phi(ns, 0.2);
                    for (size_t q = 0; q < ns && q < m.base.cellPV.size(); ++q) {
                        const int gq = G.globalOf[q];
                        if (gq >= 0 && size_t(gq) < G.bulkVol.size() && G.bulkVol[size_t(gq)] > 0) phi[q] = std::clamp(m.base.cellPV[q] / (G.bulkVol[size_t(gq)] * lm * lm * lm), 0.01, 0.6);
                    }
                    gsC = std::make_unique<geoseq::Coupler>(); std::string e;
                    if (!gsC->setup(c.geoSeq, G, G.globalOf, lm, sim.pressure(), phi, e)) { addLog("NOTE: FE geomechanics off - " + e); gsC.reset(); }
                    else addLog(std::string("=== FE geomechanics (sequential, compositional): ") + (gsC->hexMesh() ? "isoparametric hexahedra on the corner-point nodes" : "Cartesian bricks") +
                                (gsC->twoWay() ? ", two-way (pore volume + permeability after every report step)" : ", one-way"));
                    if (gsC && c.geoSeq.coupling == 2) addLog("NOTE: iterative geomechanics coupling (2) is implemented for the black-oil engine - the compositional run uses the explicit two-way scheme");
                }
                if (c.faultMech.on()) {
                    fmC = std::make_unique<faultmech::Monitor>(); std::string e;
                    if (!fmC->setup(c, G, e)) { addLog("NOTE: fault mechanics off - " + e); fmC.reset(); }
                    else {
                        addLog("=== Fault mechanics (" + fmC->out().engine + "): " + std::to_string(fmC->out().faces.size()) + " fault faces, " + (fmC->twoWay() ? "two-way" : "one-way"));
                        if (fmC->wantsFeStress()) {
                            if (gsC) { geoseq::Coupler* gp = gsC.get(); fmC->setStressProvider([gp](int cell, double S6[6]) { return gp->totalStress(cell, S6); }); }
                            else addLog("NOTE: STRESSFIELD ... FE needs GEOMECH - the analytic stress field is used");
                        }
                    }
                }
                if (anyHfGrow) {
                    hfC = std::make_unique<geoseq::HfSequential>(); std::string e;
                    if (!hfC->setup(cx, G, e)) { addLog("NOTE: HFGROW off - " + e); hfC.reset(); }
                    else addLog("=== Hydraulic-fracture growth (LEFM, report-step sequential in the compositional engine): " + std::to_string(hfC->pieces()) + " pieces, " +
                                std::to_string(hfC->openPieces()) + " open at the start" + (anyFeHfC && !gsC ? " (FESTRESS needs GEOMECH: stress gradient used)" : ""));
                }
            }
            if (!gxC) return;
            std::vector<double> sMinCase;
            if (gsC) {
                if (!gsC->step(sim.pressure())) {
                    addLog("NOTE: FE geomechanics: the mechanics solve did not converge - coupling stopped");
                    if (fmC) fmC->setStressProvider({});
                    gsC.reset();
                } else {
                    if (gsC->twoWay()) sim.setGeoFactors(gsC->pvFactor(sim.pressure(), m.base.cr), geoPermC ? gsC->permFactor() : std::vector<double>());
                    if (hfC) { sMinCase = gsC->minStress(); for (double& v : sMinCase) v /= pu; }
                    char b[200];
                    std::snprintf(b, sizeof b, "  geomechanics t = %.4g d: surface uplift %.3g mm, max failure index %.3f, %d CG its", tDays,
                                  gsC->solver().surfaceUplift() * 1e3, gsC->maxFailure(), gsC->solver().lastCg());
                    addLog(b);
                }
            }
            std::vector<float> pc(ns);
            for (size_t q = 0; q < ns; ++q) pc[q] = float(sim.pressure()[q] / pu);
            if (fmC) fmC->step(tDays, pc);
            if (hfC && step > 0) {
                std::vector<double> pcd(pc.begin(), pc.end());
                const size_t ne = hfC->events().size();
                if (hfC->grow(tDays, pcd, sMinCase) > 0)
                    for (size_t e = ne; e < hfC->events().size(); ++e) {
                        const auto& E = hfC->events()[e]; char b[200];
                        std::snprintf(b, sizeof b, "  fracture %s wing %c grows to %.4g (K_I %.3g MPa m^0.5, net pressure %.4g)", c.grid.hfracs[size_t(E.frac)].name.c_str(),
                                      E.side < 0 ? '-' : '+', E.length, E.KI * 1e-6, E.pnet);
                        addLog(b);
                    }
            }
            // transmissibility multipliers: fault reactivation x closed fracture pieces (absolute, per cell pair)
            std::map<std::pair<int, int>, double> M;
            if (hfC) hfC->connMultipliers(M);
            if (fmC) {
                for (const auto& e : fmC->connMultipliers()) {
                    const auto& cn = gxC->conns[size_t(e.first)];
                    const std::pair<int, int> key{std::min(cn.a, cn.b), std::max(cn.a, cn.b)};
                    auto it = M.find(key); if (it == M.end()) M[key] = e.second; else it->second *= e.second;
                }
            }
            if (!M.empty()) {
                std::vector<std::pair<int, int>> pr; std::vector<double> mv;
                for (const auto& e : M) { pr.push_back(e.first); mv.push_back(e.second); }
                sim.setConnectionMultipliersByCells(pr, mv);
            }
            std::shared_ptr<const faultmech::Output> fo;
            if (fmC) fo = std::make_shared<const faultmech::Output>(fmC->out());
            if (fo) { std::lock_guard<std::mutex> l(m_); r_.faultMech = fo; }
        };
        if (rq && wantGeo) addLog("NOTE: restarted compositional run - the coupled geomechanics is not continued after a restart");
        const bool ok = comp::runDeckModel(m, out, [&](const comp::CompSummary& S, comp::CompSimulator& sim) {
            coupleComp(sim, S.time);
            if (step == 0) {   // component list of the run (EOS components + mobile water)
                std::vector<ComponentInfo> ci;
                const auto& F = sim.model().fluid;
                for (int i = 0; i < sim.nc(); ++i) {
                    ComponentInfo c; c.name = i < int(F.comp.size()) ? F.comp[size_t(i)].name : "C" + std::to_string(i + 1);
                    c.mw = i < int(F.comp.size()) ? F.comp[size_t(i)].mw : 0.0;
                    if (co2 && (c.name.empty() || i == 1)) c.name = i == 0 ? "CO2" : "H2O";
                    ci.push_back(c);
                }
                if (sim.ncw() > sim.nc()) { ComponentInfo w; w.name = "H2O"; w.mw = 0.0180153; w.kind = 3; ci.push_back(w); }
                assignHeatingValues(ci);
                std::lock_guard<std::mutex> l(m_);
                r_.components = std::move(ci);
                const comp::CompCase& SM0 = sim.model();
                r_.trackName = co2 ? std::string("CO2") : (SM0.co2sol && SM0.co2Comp >= 0 && SM0.co2Comp < int(r_.components.size())) ? r_.components[size_t(SM0.co2Comp)].name : std::string();
            }
            Snapshot sn = comp::cellSnapshot(sim, S.time, pu);   // C48: shared with the CLI (engine/src/cell_snapshot.cpp)
            if (gsC) {   // C25 FE geomechanics maps / summary
                const geoseq::Coupler* gp = gsC.get();
                addMechMaps(sn, gp->solver(), [gp](size_t q) { return gp->boxOf(int(q)); }, sn.pressure.size(), pu, c.geoSeq.plasticity == 1);
                if (gsC->twoWay()) { std::vector<float> f(sn.pressure.size(), 1.0f); const auto& g = sim.geoPvFactor(); for (size_t q = 0; q < f.size() && q < g.size(); ++q) f[q] = float(g[q]); sn.extra.push_back({"PVGEO", std::move(f)}); }
            }
            SummaryPoint pt = toPoint(S);
            if (gsC) { pt.fuplift = gsC->solver().surfaceUplift(); pt.fupliftRes = gsC->solver().reservoirTopUplift(); pt.fmaxFail = gsC->maxFailure(); }
            {
                std::lock_guard<std::mutex> l(m_);
                accumulateComponentPoint(r_.summary.empty() ? nullptr : &r_.summary.back(), pt);
                r_.summary.push_back(pt);
                r_.snaps.push_back(std::move(sn));
                if (step == 0) { r_.initial.oil = S.foip; r_.initial.water = S.fwip; r_.initial.gasFree = S.fgip; r_.valid = true; }
                const comp::CompStats& cs = sim.stats();
                r_.stats.timesteps = cs.steps; r_.stats.newton = cs.newton; r_.stats.linear = cs.linear; r_.stats.cuts = cs.cuts; r_.stats.cpr = true;
                progress_ = total > 0 ? S.time / total : 0;
                wall_ = std::chrono::duration<double>(clock::now() - t0).count();
                ++version_;
            }
            ++step;
            return !cancel_.load();
        }, &err, -1, &st, &dr);
        if (fmC && !cancel_) {   // C25: post-run RSF traces of the compositional run
            fmC->finish(&cancel_, [this](const std::string& s2) { addLog(s2); });
            auto fo = std::make_shared<const faultmech::Output>(fmC->out());
            std::lock_guard<std::mutex> l(m_); r_.faultMech = fo; ++version_;
        }
        std::snprintf(buf, sizeof(buf), "=== %s: %d timesteps, %d cuts, %d Newton, %d linear iterations, %d flashes, %.2f s",
                      !ok ? "Run failed" : cancel_ ? "Run stopped" : "Run finished", st.steps, st.cuts, st.newton, st.linear, st.flashes,
                      std::chrono::duration<double>(clock::now() - t0).count());
        addLog(buf);
        if (!ok) throw std::runtime_error(err);
    } catch (const std::exception& e) {
        std::lock_guard<std::mutex> l(m_);
        error_ = e.what();
        log_.push_back(std::string("ERROR: ") + e.what());
    }
    wall_ = std::chrono::duration<double>(clock::now() - t0).count();
    running_ = false;
    ++version_;
}

// ================================================================================================ thermal core (E300 THERMAL)
#include "ressim/therm_deck.h"
void SimRunner::startThermal(const ressim::Case& c) {
    if (running_) return;
    join();
    cancel_ = false; running_ = true; progress_ = 0; wall_ = 0;
    {
        std::lock_guard<std::mutex> l(m_);
        error_.clear();
        r_ = RunResults{};
    }
    rqActive_ = std::move(rq_); rq_.reset();
    ++version_;
    th_ = std::thread(&SimRunner::workThermal, this, c);
}

void SimRunner::workThermal(ressim::Case c) {
    using clock = std::chrono::steady_clock;
    using namespace ressim;
    const auto t0 = clock::now();
    try {
        addLog("=== Thermal run started: " + c.title + " (ResSim thermal core: FIM p-N-U, IAPWS-IF97 steam, K-value / free-water flash, "
               "Vinsome-Westerveld heat loss)");
        thermal::DeckModel m; std::string err;
        if (!thermal::deckModel(c, m, &err)) throw std::runtime_error("thermal model: " + err);
        char buf[320];
        std::snprintf(buf, sizeof(buf), "Thermal model: %d cells, %zu connections, %zu components + water, %zu wells, %zu report steps, heat loss %s",
                      m.base.size(), m.base.conns.size(), m.base.comps.size(), m.wellNames.size(), m.stepDays.size(), m.base.heatLoss ? "on" : "off");
        addLog(buf);
        for (const auto& n : m.notes) addLog("NOTE: " + n);
        double total = 0; for (double d : m.stepDays) total += d;
        auto grid = std::make_shared<Grid>(); grid->build(c);
        {
            std::lock_guard<std::mutex> l(m_);
            r_.grid = grid; r_.input = c; r_.totalTime = total; r_.engine = "ResSim thermal";
            for (auto& w : c.wells) { r_.wellNames.push_back(w.name); r_.wellTypes.push_back(w.type); }
        }
        const bool fld = c.units == UnitSystem::Field;
        const double pu = fld ? 6894.757293168 : 1e5;
        auto tOut = [&](double K) { return fld ? K * 1.8 - 459.67 : K - 273.15; };
        const int ng = grid->numActive();
        size_t step = 0;
        std::vector<thermal::Summary> out; thermal::Stats st;
        double lastT = 0;
        std::vector<float> T0map;      // temperature map at t = 0 (cooling map DTEMP)
        std::vector<double> Tstart;    // initial cell temperatures (K) for the cooled-volume fraction
        const std::shared_ptr<RestartRequest> rq = rqActive_;
        DeckRestart dr; dr.saveEvery = rstEvery_; dr.saveLast = rstLast_;
        dr.onSave = [&](int stp, double days, const std::string& blob) {
            std::lock_guard<std::mutex> l(m_);
            ressim::ResultsFile::Restart R; R.report = stp; R.time = days; R.kind = "THERMAL"; R.state = blob;
            R.summaryCount = r_.summary.size() + 1; R.snap = int(r_.snaps.size());
            r_.restarts.push_back(std::move(R));
        };
        if (rq) {
            if (rq->kind != "THERMAL") throw std::runtime_error("restart: the saved state is from the " + rq->kind + " engine, not thermal");
            dr.startStep = rq->step; dr.state = rq->state;
            std::lock_guard<std::mutex> l(m_);
            prefillFromBase(*rq);
            step = size_t(rq->step) + 1;
            lastT = r_.summary.empty() ? 0.0 : r_.summary.back().time;
            if (!r_.snaps.empty())
                for (const auto& e : r_.snaps[0].extra)
                    if (e.first == "TEMP") {
                        T0map = e.second;
                        Tstart.resize(e.second.size());
                        for (size_t q = 0; q < e.second.size(); ++q) Tstart[q] = fld ? (e.second[q] + 459.67) / 1.8 : e.second[q] + 273.15;
                    }
            addLog("=== Restart from report step " + std::to_string(rq->step) + " - continuing the base run exactly");
        }
        // ---- C25 fault mechanics + seismicity of the thermal (EGS) run: thermo-poro-elastic Coulomb stress on the FAULTS faces
        // (analytic field + poroelastic path + thermo-elastic dS_h, or the THM mechanics stress with STRESSFIELD ... FE),
        // two-way fault transmissibility, post-run RSF traces
        std::unique_ptr<faultmech::Monitor> fmT;
        const thermal::Simulator* simT = nullptr;
        if (c.faultMech.on()) {
            fmT = std::make_unique<faultmech::Monitor>(); std::string e;
            if (!fmT->setup(c, *grid, e)) { addLog("NOTE: fault mechanics off - " + e); fmT.reset(); }
            else {
                addLog("=== Fault mechanics (thermal): " + std::to_string(fmT->out().faces.size()) + " fault faces, " + (fmT->twoWay() ? "two-way" : "one-way") +
                       (fmT->out().thermalStress ? ", thermo-elastic stress" : ""));
                if (fmT->wantsFeStress()) {
                    if (c.thermal.geo.enabled)
                        fmT->setStressProvider([&simT](int cell, double S6[6]) {
                            const geomech::Solver* g = simT ? simT->mechanics() : nullptr;
                            if (!g) return false;
                            const auto& gm = simT->mechCellMap();
                            if (cell < 0 || size_t(cell) >= gm.size() || gm[size_t(cell)] < 0) return false;
                            g->totalStress(gm[size_t(cell)], S6);
                            return true;
                        });
                    else addLog("NOTE: STRESSFIELD ... FE needs the THM mechanics (GEOMECH) - the analytic stress field is used");
                }
            }
        }
        if (rq && fmT) addLog("NOTE: restarted thermal run - fault mechanics starts from the restart state");
        const bool ok = thermal::runDeckModel(m, out, [&](const thermal::Summary& S, thermal::Simulator& sim) {
            simT = &sim;
            if (step == 0) {   // component list of the run (K-value / gas / solid components + water)
                std::vector<ComponentInfo> ci;
                for (const auto& q : sim.model().comps) {
                    ComponentInfo c; c.name = q.name; c.mw = q.mw; c.kind = q.type == 2 ? 2 : q.type == 1 ? 1 : 0;
                    ci.push_back(c);
                }
                ComponentInfo w; w.name = "H2O"; w.mw = 0.0180153; w.kind = 3; ci.push_back(w);
                assignHeatingValues(ci);
                std::lock_guard<std::mutex> l(m_);
                r_.components = std::move(ci);
                const int tk = sim.model().track;
                r_.trackName = tk >= 0 ? sim.model().comps[size_t(tk)].name : std::string();
            }
            SummaryPoint p;
            p.time = S.time; p.dt = S.time - lastT; lastT = S.time;
            p.fopr = S.fopr; p.fwpr = S.fwpr; p.fgpr = S.fgpr; p.fwir = S.fwir;
            p.fopt = S.fopt; p.fwpt = S.fwpt; p.fgpt = S.fgpt; p.fwit = S.fwit;
            p.fgir = S.fgir; p.fgit = S.fgit;
            p.fpr = S.fpr; p.fprp = S.fpr; p.oip = S.foip;
            p.fwct = S.fopr + S.fwpr > 0 ? S.fwpr / (S.fopr + S.fwpr) : 0.0;
            p.fgor = S.fopr > 0 ? S.fgpr / S.fopr : 0.0;
            p.favgT = S.favgT; p.fmaxT = S.fmaxT; p.fsteamPV = S.fsteamPV; p.fhlt = S.fhlt; p.fenergyIn = S.fenergyIn; p.csor = S.csor;
            p.freactRate = sim.reactionHeatRate(); p.freactCum = sim.cumReactionHeat();
            if (sim.model().track >= 0) {   // storage / trapping inventory of the tracked component (same classes as the CO2 cores)
                const thermal::Inventory I = sim.trackedInventory();
                p.co2Mob = I.mobile + I.structural; p.co2Struct = I.structural; p.co2Trap = I.residual; p.co2Diss = I.dissolvedAq;
                p.co2Oil = I.dissolvedOil; p.co2Min = I.mineral; p.co2KrImmob = I.immobileKr;
            }
            if (const geomech::Solver* g = sim.mechanics()) {
                p.fuplift = g->surfaceUplift(); p.fupliftRes = g->reservoirTopUplift();
                for (double f : g->failure()) p.fmaxFail = std::max(p.fmaxFail, f);
            }
            const std::vector<thermal::Well>* ws = step > 0 && step - 1 < m.stepWells.size() ? &m.stepWells[step - 1] : nullptr;
            for (size_t w = 0; w < S.wells.size(); ++w) {
                WellRates r; const auto& W = S.wells[w];
                r.oil = W.oil; r.gas = W.gas; r.water = W.water; r.waterInj = W.winj; r.bhp = W.bhp; r.open = W.open;
                r.mode = Control::BHP;
                if (ws && w < ws->size() && !W.onBhp) {
                    const thermal::Well& tw = (*ws)[w];
                    if (tw.injector) r.mode = Control::RATE;
                    else switch (tw.control) { case 1: r.mode = Control::ORAT; break; case 2: r.mode = Control::WRAT; break; case 3: r.mode = Control::LRAT; break; default: break; }
                }
                r.comp = W.molar;
                addCompRates(p, r.comp);
                p.wells.push_back(r);
            }
            // cell maps on the full active grid (simulator cells are the active cells in order)
            const size_t n = size_t(ng);
            if (step == 0) T0map.clear();
            Snapshot sn = thermal::cellSnapshot(sim, c.thermal, S.time, n, fld, T0map);   // C48: shared with the CLI
            if (fmT) {   // C25 fault mechanics after the report step (temperature in K: only changes are used)
                std::vector<float> Tc(n, 0.0f);
                for (size_t q = 0; q < n && q < sim.temperature().size(); ++q) Tc[q] = float(sim.temperature()[q]);
                fmT->step(S.time, sn.pressure, &Tc);
                const std::map<int, double> M = fmT->connMultipliers();
                if (!M.empty()) {
                    std::vector<std::pair<int, int>> pr; std::vector<double> mv;
                    for (const auto& e : M) { const auto& cn = grid->conns[size_t(e.first)]; pr.push_back({cn.a, cn.b}); mv.push_back(e.second); }
                    sim.setConnectionMultipliersByCells(pr, mv);
                    ++fmT->out().twoWayUpdates;
                }
            }
            // ---- C13: geothermal / EDFM analysis data
            {
                const auto& FC = sim.model().fracCells;
                const double lenU = fld ? 0.3048 : 1.0;
                const int nbox = sim.model().nx * sim.model().ny * sim.model().nz;
                RunResults::FracState fs;
                if (!FC.empty()) {
                    const std::vector<double> wf = sim.fractureApertures();
                    fs.ap.resize(FC.size()); fs.p.resize(FC.size()); fs.T.resize(FC.size()); fs.cond.resize(FC.size());
                    for (size_t f = 0; f < FC.size(); ++f) {
                        const double w = f < wf.size() ? wf[f] : FC[f].w0;
                        const bool cubic = std::fabs(FC[f].kf - FC[f].w0 * FC[f].w0 / 12.0) <= 1e-6 * FC[f].kf;
                        const double kf = cubic ? w * w / 12.0 : FC[f].kf;
                        const size_t cc = size_t(FC[f].cell);
                        fs.ap[f] = float(w * 1e3);
                        fs.p[f] = cc < sim.pressure().size() ? float(sim.pressure()[cc] / pu) : 0.0f;
                        fs.T[f] = cc < sim.temperature().size() ? float(tOut(sim.temperature()[cc])) : 0.0f;
                        fs.cond[f] = float(kf * w / 9.869233e-16);
                    }
                }
                std::vector<float> wt(S.wells.size()), wmw(S.wells.size());
                const double eJ = fld ? 1055.05585 : 1000.0;
                for (size_t w = 0; w < S.wells.size(); ++w) { wt[w] = float(S.wells[w].tfluid); wmw[w] = float(S.wells[w].hrate * eJ / 86400.0 * 1e-6); }
                // matrix cooling (> 10 K) and the accessible heat in place at t = 0
                const auto& Tk = sim.temperature();
                double vb = 0, vc = 0;
                if (step == 0) Tstart.assign(Tk.begin(), Tk.begin() + std::min<size_t>(Tk.size(), size_t(nbox)));
                for (int q = 0; q < nbox && size_t(q) < Tk.size() && size_t(q) < Tstart.size(); ++q) {
                    const double b = sim.bulkVolume(q);
                    vb += b; if (Tstart[size_t(q)] - Tk[size_t(q)] > 10.0) vc += b;
                }
                std::lock_guard<std::mutex> l(m_);
                if (step == 0) {
                    r_.fracNames.clear(); for (const auto& F : c.thermal.fractures) r_.fracNames.push_back(F.name);
                    r_.fracSegs.clear();
                    for (const auto& f : FC) {
                        RunResults::FracSeg g; g.frac = f.frac; g.host = f.host;
                        g.xa = float(f.xa / lenU); g.ya = float(f.ya / lenU); g.xb = float(f.xb / lenU); g.yb = float(f.yb / lenU);
                        g.zt = float(f.zt / lenU); g.zb = float(f.zb / lenU);
                        r_.fracSegs.push_back(g);
                    }
                    r_.wellFracLinks.clear();
                    for (size_t w = 0; w < sim.model().wells.size(); ++w)
                        for (int q : sim.model().wells[w].cells)
                            for (const auto& f : FC)
                                if (f.cell == q && f.frac >= 0 && size_t(f.frac) < r_.fracNames.size()) {
                                    std::pair<std::string, std::string> L{m.wellNames[w], r_.fracNames[size_t(f.frac)]};
                                    if (std::find(r_.wellFracLinks.begin(), r_.wellFracLinks.end(), L) == r_.wellFracLinks.end()) r_.wellFracLinks.push_back(L);
                                }
                    // heat in place above the (lowest) injection temperature: rock (per rock volume) + water
                    double tref = 1e9;
                    for (const auto& w : sim.model().wells) if (w.injector && w.tinj > 0) tref = std::min(tref, w.tinj);
                    for (const auto& wj : c.thermal.winj) if (wj.temp > -200) tref = std::min(tref, fld ? (wj.temp + 459.67) * 5.0 / 9.0 : wj.temp + 273.15);
                    if (tref > 1e8) tref = 288.15;
                    double hip = 0;
                    const auto& M = sim.model();
                    for (int q = 0; q < nbox && size_t(q) < Tk.size(); ++q) {
                        const double cv = M.rockCvCell.size() > size_t(q) ? M.rockCvCell[size_t(q)] : M.rockCv;
                        const double pv = size_t(q) < M.cellPV.size() ? M.cellPV[size_t(q)] : 0.0, bk = sim.bulkVolume(q);
                        hip += ((bk - pv) * cv + pv * 4.18e6) * std::max(0.0, Tk[size_t(q)] - tref);
                    }
                    r_.heatInPlace = hip; r_.tInjRef = tOut(tref);
                    r_.fracSnap.clear(); r_.wellT.clear(); r_.wellMW.clear(); r_.cooledFrac.clear();
                }
                r_.fracSnap.push_back(std::move(fs)); r_.wellT.push_back(std::move(wt)); r_.wellMW.push_back(std::move(wmw));
                r_.cooledFrac.push_back(float(vb > 0 ? vc / vb : 0.0));
            }
            {
                std::lock_guard<std::mutex> l(m_);
                accumulateComponentPoint(r_.summary.empty() ? nullptr : &r_.summary.back(), p);
                if (fmT) r_.faultMech = std::make_shared<const faultmech::Output>(fmT->out());
                r_.summary.push_back(p);
                r_.snaps.push_back(std::move(sn));
                if (step == 0) { r_.initial.oil = S.foip; r_.valid = true; }
                const thermal::Stats& cs = sim.stats();
                r_.stats.timesteps = cs.steps; r_.stats.newton = cs.newton; r_.stats.linear = cs.linear; r_.stats.cuts = cs.cuts; r_.stats.cpr = true;
                progress_ = total > 0 ? S.time / total : 0;
                wall_ = std::chrono::duration<double>(clock::now() - t0).count();
                ++version_;
            }
            ++step;
            return !cancel_.load();
        }, &err, -1, &st, &dr);
        simT = nullptr;
        if (fmT && !cancel_) {   // C25: post-run RSF traces of the thermal run
            fmT->finish(&cancel_, [this](const std::string& s2) { addLog(s2); });
            auto fo = std::make_shared<const faultmech::Output>(fmT->out());
            std::lock_guard<std::mutex> l(m_); r_.faultMech = fo; ++version_;
        }
        std::snprintf(buf, sizeof(buf), "=== %s: %d timesteps, %d cuts, %d Newton, %d linear iterations, %d flashes, %.2f s",
                      !ok ? "Run failed" : cancel_ ? "Run stopped" : "Run finished", st.steps, st.cuts, st.newton, st.linear, st.flashes,
                      std::chrono::duration<double>(clock::now() - t0).count());
        addLog(buf);
        if (!ok) throw std::runtime_error(err);
    } catch (const std::exception& e) {
        std::lock_guard<std::mutex> l(m_);
        error_ = e.what();
        log_.push_back(std::string("ERROR: ") + e.what());
    }
    wall_ = std::chrono::duration<double>(clock::now() - t0).count();
    running_ = false;
    ++version_;
}

// ================================================================================================ OPM Flow engine
#include "ressim/ecl.h"
#include "ressim/io.h"
#include <filesystem>
#include <sstream>

void SimRunner::startOpm(const ressim::Case& c, const OpmOptions& o) {
    if (running_) return;
    join();
    cancel_ = false; running_ = true; progress_ = 0; wall_ = 0;
    {
        std::lock_guard<std::mutex> l(m_);
        error_.clear();
        r_ = RunResults{};
    }
    ++version_;
    th_ = std::thread(&SimRunner::workOpm, this, c, o);
}

void SimRunner::workOpm(ressim::Case c, OpmOptions o) {
    namespace fs = std::filesystem;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    auto fail = [&](const std::string& msg) {
        std::lock_guard<std::mutex> l(m_);
        error_ = msg;
        log_.push_back("ERROR: " + msg);
    };
    const bool ecl = o.kind == 1;
    const std::string engName = ecl ? "ECLIPSE" : "OPM Flow";
    try {
        if (ecl) {
            if (o.flowPath.empty()) o.flowPath = detectEclipseExecutable();
            if (o.flowPath.empty())
                throw std::runtime_error("ECLIPSE (eclrun.exe / eclipse.exe) was not found. Choose the executable in the simulation "
                                         "engine settings (Run menu).");
            if (!fs::is_regular_file(fs::u8path(o.flowPath)))
                throw std::runtime_error("ECLIPSE executable not found: " + o.flowPath + " (simulation engine settings)");
            // direct executable (no eclrun): use the selected program (e300 / frontsim) from the same bin folder
            if (!isEclrunLauncher(o.flowPath) && !o.eclProgram.empty()) {
                const fs::path cur = fs::u8path(o.flowPath);
                const std::string want = o.eclProgram + cur.extension().u8string();
                std::string stem = cur.stem().u8string();
                for (char& ch : stem) ch = char(std::tolower((unsigned char)ch));
                if (stem != o.eclProgram && o.eclProgram != "eclipse") {
                    std::error_code fe;
                    const fs::path alt = cur.parent_path() / want;
                    if (fs::is_regular_file(alt, fe)) o.flowPath = alt.u8string();
                    else throw std::runtime_error(o.eclProgram + " was selected but " + want + " is not next to " + o.flowPath +
                                                  " - choose eclrun.exe (recommended) or the " + o.eclProgram + " executable in the engine settings");
                }
            }
        } else {
            if (o.flowPath.empty()) o.flowPath = detectFlowExecutable();
            if (o.flowPath.empty())
                throw std::runtime_error("OPM Flow (flow.exe) was not found. Install OPM Flow for Windows or choose the executable in "
                                         "Run > OPM Flow settings.");
            if (!fs::is_regular_file(fs::u8path(o.flowPath)))
                throw std::runtime_error("OPM Flow executable not found: " + o.flowPath + " (Run > OPM Flow settings)");
            const std::string pre = flowPrerequisiteProblem(o.flowPath);
            if (!pre.empty()) addLog((pre.rfind("Note:", 0) == 0 ? std::string() : std::string("WARNING: ")) + pre);
            for (auto& n : opmCompatibilityNotes(c)) addLog("NOTE: " + n);
        }
        const fs::path dir = fs::u8path(o.runDir);
        std::error_code ec;
        fs::create_directories(dir, ec);
        if (!fs::is_directory(dir)) throw std::runtime_error("Cannot create the " + engName + " run folder " + o.runDir);
        // Original-deck mode (imported deck not edited in the UI): the source deck and its INCLUDE files are run unchanged,
        // so everything the external simulator supports works even where ResSim's model cannot represent it.
        const bool passThrough = !o.sourceDeck.empty() && fs::is_regular_file(fs::u8path(o.sourceDeck));
        const fs::path src = passThrough ? fs::absolute(fs::u8path(o.sourceDeck), ec) : fs::path();
        std::string base = passThrough ? src.stem().u8string() : caseBaseName(c);
        // ECLIPSE writes next to the dataset (run in the deck folder); OPM writes to --output-dir
        const fs::path outDir = passThrough && ecl ? src.parent_path() : dir;
        // stale output of an earlier run must not be read back
        for (const char* ext : {".SMSPEC", ".UNSMRY", ".UNRST", ".EGRID", ".INIT", ".PRT", ".DBG", ".RSM", ".ESMRY", ".MSG", ".ECLEND", ".LOG", ".h5"})
            fs::remove(outDir / (base + ext), ec);
        const fs::path deck = passThrough ? (ecl ? src : dir / (base + ".DATA")) : dir / (base + ".DATA");
        if (!passThrough) {
            if (ecl && !c.grid.disperc.empty()) {   // C17: DISPERC is an OPM Flow keyword - ECLIPSE would stop on it
                ressim::Case ce = c; ce.grid.disperc.clear();
                ressim::exportEclipseDeck(ce, deck.u8string());
                addLog("NOTE: DISPERC (mechanical dispersion, OPM Flow keyword) is not written for ECLIPSE - the ECLIPSE run has no dispersion");
            } else ressim::exportEclipseDeck(c, deck.u8string());
        }
        addLog("=== " + engName + " run started: " + c.title);
        addLog("    engine  " + o.flowPath);
        if (passThrough) {
            addLog("    deck    " + src.u8string() + "  (original imported deck, run unchanged with its INCLUDE files)");
            if (!ecl) addLog("    output  " + dir.u8string());
        } else addLog("    deck    " + deck.u8string() + "  (generated from the model)");
        double total = 0;
        for (auto& s : c.schedule) total += s.length;
        {
            std::lock_guard<std::mutex> l(m_);
            r_.input = c;
            r_.totalTime = total;
            r_.engine = engName;
            r_.opmDeck = deck.u8string();
            for (auto& w : c.wells) { r_.wellNames.push_back(w.name); r_.wellTypes.push_back(w.type); }
        }
        std::vector<std::string> args;
        std::string stdinText;
        if (ecl) {
            // eclrun [options] [-v VERSION] eclipse|e300 CASE   - or the simulator itself, which asks for the case name
            if (isEclrunLauncher(o.flowPath)) {
                { std::istringstream ss(o.extraArgs); std::string a; while (ss >> a) args.push_back(a); }
                if (!o.eclVersion.empty()) { args.push_back("-v"); args.push_back(o.eclVersion); }
                args.push_back(o.eclProgram.empty() ? std::string("eclipse") : o.eclProgram);
                args.push_back(base);
            } else {
                { std::istringstream ss(o.extraArgs); std::string a; while (ss >> a) args.push_back(a); }
                stdinText = base + "\n" + base + "\n";
                addLog("    input   case name '" + base + "' on standard input");
            }
        } else {
            args = {passThrough ? src.u8string() : base + ".DATA", "--output-dir=."};
            int threads = o.threads > 0 ? o.threads : std::min(8, int(std::max(1u, std::thread::hardware_concurrency())));
            args.push_back("--threads-per-process=" + std::to_string(threads));
            if (o.maxDtDays > 0) { std::ostringstream s; s << "--solver-max-time-step-in-days=" << o.maxDtDays; args.push_back(s.str()); }
            if (o.enableTuning) args.push_back("--enable-tuning=true");
            if (o.keepWellsOpen) {   // OPM shuts wells whose equations fail to converge - silently changes the physics
                args.push_back("--shut-unsolvable-wells=false");
                args.push_back("--min-time-step-before-shutting-problematic-wells-in-days=1e-12");
            }
            { std::istringstream ss(o.extraArgs); std::string a; while (ss >> a) args.push_back(a); }
        }
        std::string cmdLine = "    command " + o.flowPath;
        for (auto& a : args) cmdLine += " " + a;
        addLog(cmdLine);

        ChildProcess proc;
        std::string err;
        if (!proc.start(o.flowPath, args, outDir.u8string(), &err, stdinText)) throw std::runtime_error(engName + " " + err + ": " + o.flowPath);
        { std::lock_guard<std::mutex> l(procM_); proc_ = &proc; }
        if (cancel_) proc.kill();
        OpmProgress prog;
        EclProgress eprog;
        const std::string smspec = (outDir / (base + ".SMSPEC")).u8string();
        auto lastLive = clock::now();
        int lastReport = -1;
        auto liveSummary = [&]() {
            // summary files are rewritten at report steps; a half-written file is simply skipped
            ressim::ecl::Summary s; std::string e;
            try {
                if (!ressim::ecl::readSummary(smspec, s, &e) || s.time.empty()) return;
                std::lock_guard<std::mutex> l(m_);
                wellsFromSummary(s, r_.wellNames, r_.wellTypes);
                auto pts = ressim::ecl::runFromSummary(s, r_.wellNames);
                r_.summary = std::move(pts);
                ++version_;
            } catch (...) {}
        };
        const int rc = proc.pump([&](const std::string& line) {
            {
                std::lock_guard<std::mutex> l(m_);
                log_.push_back(line);
                if (!ecl && (line.find("will be shut because it fails to converge") != std::string::npos ||
                    line.find("Problematic well(s) were shut") != std::string::npos))
                    log_.push_back("WARNING: OPM Flow shut a well because its equations did not converge - the results no longer "
                                   "represent the deck. Enable 'Keep problematic wells open' in the OPM engine settings.");
            }
            const bool changed = ecl ? eprog.parse(line) : prog.parse(line);
            if (changed) {
                if (ecl) { prog.day = eprog.day; prog.timesteps = eprog.timesteps; prog.newton = eprog.newton; prog.reportStep = eprog.timesteps; }
                const double tot = !ecl && prog.totalDays > 0 ? prog.totalDays : total;
                if (tot > 0) progress_ = std::min(1.0, prog.day / tot);
                wall_ = std::chrono::duration<double>(clock::now() - t0).count();
                ++version_;
                if (prog.reportStep != lastReport && clock::now() - lastLive > std::chrono::milliseconds(1000)) {
                    lastReport = prog.reportStep; lastLive = clock::now();
                    liveSummary();
                }
            }
        });
        { std::lock_guard<std::mutex> l(procM_); proc_ = nullptr; }
        const double secs = std::chrono::duration<double>(clock::now() - t0).count();
        if (cancel_) addLog("=== " + engName + " run stopped by the user");
        else if (ecl) {
            // eclrun may return 0 although the simulator stopped on input errors: judge by the messages
            const int nErr = std::max(eprog.errors + eprog.bugs, eprog.reportedErrors);
            std::string why = eprog.lastError;
            if (why.empty() && (nErr > 0 || rc != 0)) why = eclipsePrtErrors((outDir / (base + ".PRT")).u8string());
            if (rc != 0 || nErr > 0 || eprog.licenseProblem) {
                std::string msg = "ECLIPSE failed";
                if (rc != 0) msg += " (exit code " + std::to_string(rc) + ")";
                if (nErr > 0) msg += ", " + std::to_string(nErr) + " error(s)";
                if (eprog.licenseProblem) msg += " - no ECLIPSE licence could be checked out (check the licence server / SLBSLS_LICENSE_FILE)";
                if (!why.empty()) msg += ": " + why;
                while (!msg.empty() && msg.back() == '.') msg.pop_back();
                msg += ". See the log and " + base + ".PRT in " + outDir.u8string();
                fail(msg);
            }
        } else if (rc != 0) {
            std::string msg = "OPM Flow failed (exit code " + std::to_string(rc) + ")";
            if (rc == -1073741515 || unsigned(rc) == 0xC0000135u)
                msg += ": a required DLL is missing - most likely the Microsoft MPI runtime (msmpi.dll). Install it with "
                       "'winget install Microsoft.msmpi' or redist\\msmpisetup.exe from the OPM Flow folder";
            else if (!prog.lastError.empty()) msg += ": " + prog.lastError;
            while (!msg.empty() && msg.back() == '.') msg.pop_back();
            msg += ". See the log and " + base + ".PRT in " + outDir.u8string();
            fail(msg);
        }
        // read back whatever was written (complete run, or partial after stop / failure)
        ressim::ecl::ResultCase rcs;
        std::string e;
        if (fs::exists(fs::u8path(smspec)) && ressim::ecl::openResultCase(deck.u8string(), rcs, &e)) {
            std::vector<std::string> wn; std::vector<ressim::WellType> wt;
            { std::lock_guard<std::mutex> l(m_); wn = r_.wellNames; wt = r_.wellTypes; }
            wellsFromSummary(rcs.summary, wn, wt);
            auto pts = ressim::ecl::runFromSummary(rcs.summary, wn);
            auto snaps = ressim::ecl::snapshotsFromRestart(rcs);
            for (auto& w : rcs.warnings) addLog("WARNING: " + w);
            std::lock_guard<std::mutex> l(m_);
            r_.wellNames = wn; r_.wellTypes = wt;
            r_.summary = std::move(pts);
            r_.snaps = std::move(snaps);
            r_.grid = rcs.grid;
            if (!r_.summary.empty()) {
                const auto& s0 = r_.summary.front();
                r_.initial.oil = s0.oip; r_.initial.water = s0.wip; r_.initial.gasFree = s0.gip; r_.initial.gasDissolved = 0;
                const int ol = rcs.summary.find("FOIPL"), og = rcs.summary.find("FOIPG");
                if (ol >= 0 && og >= 0 && !rcs.summary.values[size_t(og)].empty()) {
                    r_.initial.oil = rcs.summary.values[size_t(ol)][0]; r_.initial.oilVaporized = rcs.summary.values[size_t(og)][0];
                }
                const int a = rcs.summary.find("FGIPG"), b = rcs.summary.find("FGIPL");
                if (a >= 0 && b >= 0 && !rcs.summary.values[size_t(a)].empty()) {
                    r_.initial.gasFree = rcs.summary.values[size_t(a)][0]; r_.initial.gasDissolved = rcs.summary.values[size_t(b)][0];
                }
            }
            r_.stats.timesteps = prog.timesteps ? prog.timesteps : int(r_.summary.size());
            r_.stats.newton = prog.newton; r_.stats.linear = prog.linear; r_.stats.seconds = secs;
            r_.valid = r_.grid && !r_.snaps.empty();
            if (!r_.valid && !r_.summary.empty() && r_.grid) r_.valid = true;
            char buf[320];
            std::snprintf(buf, sizeof(buf), "=== %s %s: %d timesteps, %d Newton, %d linear iterations, %.2f s  (%zu report steps read)",
                          engName.c_str(), cancel_ ? "stopped" : error_.empty() ? "finished" : "failed", r_.stats.timesteps, prog.newton, prog.linear, secs, r_.snaps.size());
            log_.push_back(buf);
        } else if (rc == 0 && !cancel_ && error().empty()) {
            fail(engName + " finished but its results could not be read: " + (e.empty() ? smspec : e) +
                 (ecl ? " (ECLIPSE writes them only when the run gets past the input checks - see " + base + ".PRT)" : std::string()));
        }
    } catch (const std::exception& ex) {
        fail(ex.what());
    }
    { std::lock_guard<std::mutex> l(procM_); proc_ = nullptr; }
    wall_ = std::chrono::duration<double>(clock::now() - t0).count();
    progress_ = 1.0;
    running_ = false;
    ++version_;
}

// ---------------------------------------------------------------- C18 external ResSim process (MPI domain decomposition)
#include "results_io.h"
void SimRunner::startCli(const ressim::Case& c, const CliRunSpec& s) {
    if (running_) return;
    join();
    cancel_ = false; running_ = true; progress_ = 0; wall_ = 0;
    {
        std::lock_guard<std::mutex> l(m_);
        error_.clear();
        r_ = RunResults{};
    }
    ++version_;
    th_ = std::thread(&SimRunner::workCli, this, c, s);
}

void SimRunner::workCli(ressim::Case c, CliRunSpec s) {
    namespace fs = std::filesystem;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    try {
        double total = 0; for (auto& st : c.schedule) total += st.length;
        std::string cmd = "    command " + s.exe;
        for (auto& a : s.args) cmd += " " + a;
        addLog("=== " + s.engineName + " started");
        addLog(cmd);
        std::error_code ec;
        fs::remove(fs::u8path(s.resultsFile), ec);
        ChildProcess proc;
        std::string err;
        if (!proc.start(s.exe, s.args, s.workDir, &err)) throw std::runtime_error(s.engineName + ": " + err + ": " + s.exe);
        { std::lock_guard<std::mutex> l(procM_); proc_ = &proc; }
        if (cancel_) proc.kill();
        std::string lastError;
        const int rc = proc.pump([&](const std::string& line) {
            { std::lock_guard<std::mutex> l(m_); log_.push_back(line); }
            double t = 0;
            if (std::sscanf(line.c_str(), "  t = %lf d", &t) == 1 && total > 0) { progress_ = std::min(1.0, t / total); wall_ = std::chrono::duration<double>(clock::now() - t0).count(); ++version_; }
            if (line.rfind("ERROR", 0) == 0 || line.find("MPI:") != std::string::npos) lastError = line;
        });
        { std::lock_guard<std::mutex> l(procM_); proc_ = nullptr; }
        if (cancel_) addLog("=== " + s.engineName + " run stopped by the user");
        else if (rc != 0 || !fs::exists(fs::u8path(s.resultsFile), ec)) {
            std::string msg = s.engineName + " failed (exit code " + std::to_string(rc) + ")";
            if (!lastError.empty()) msg += ": " + lastError;
            msg += ". Check MPI with Settings -> Parallel computing -> Test MPI.";
            std::lock_guard<std::mutex> l(m_);
            error_ = msg; log_.push_back("ERROR: " + msg);
        } else {
            RunResults r = loadRunResults(s.resultsFile);
            r.engine = s.engineName;
            r.input = c;
            std::lock_guard<std::mutex> l(m_);
            r_ = std::move(r);
            log_.push_back("=== " + s.engineName + " results loaded from " + s.resultsFile);
        }
    } catch (const std::exception& ex) {
        std::lock_guard<std::mutex> l(m_);
        error_ = ex.what(); log_.push_back(std::string("ERROR: ") + ex.what());
    }
    { std::lock_guard<std::mutex> l(procM_); proc_ = nullptr; }
    wall_ = std::chrono::duration<double>(clock::now() - t0).count();
    progress_ = 1.0;
    running_ = false;
    ++version_;
}
