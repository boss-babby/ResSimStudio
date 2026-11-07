// C28: headless coupled compositional / thermal runs (see coupled_deck.h). Closes the C24-a / C26-1 open item "FAULTSLIP
// for the compositional and thermal engines" outside the GUI: the CLI, scripts and tests run the same fault mechanics
// (incl. DDM slip + stress transfer) as the GUI runners, with the thermo-elastic or THM-mechanics stress of thermal runs.
#include "ressim/coupled_deck.h"
#include "ressim/faultmech.h"
#include "ressim/geoseq.h"
#include "ressim/geomech.h"
#include "ressim/grid.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>

namespace ressim {
namespace coupleddeck {

bool needed(const Case& c) {
    if (!c.thermal.enabled && !c.comp.enabled) return false;
    return c.faultMech.on() || (c.comp.enabled && !c.thermal.enabled && c.geoSeq.enabled);
}

namespace {
double pressureUnitPa(const Case& c) { return c.units == UnitSystem::Field ? 6894.757 : 1e5; }

// grid of the simulator cells: LGR / MINC / DPNUM / HF expanded exactly as the deck engines do
std::shared_ptr<Grid> engineGrid(const Case& c, size_t cells, std::string& err) {
    Case cx = c;
    auto g = std::make_shared<Grid>();
    try { expandLgrArrays(cx); expandMincArrays(cx); applyDpnum(cx); prepareHydraulicFractures(cx); g->build(cx); }
    catch (const std::exception& e) { err = std::string("grid: ") + e.what(); return nullptr; }
    if (size_t(g->numActive()) != cells) { err = "simulator cells do not match the grid"; return nullptr; }
    return g;
}

template <class Sim>
void applyFaultMultipliers(const faultmech::Monitor& fm, const Grid& G, Sim& sim, faultmech::Monitor& fmw) {
    const std::map<int, double> M = fm.connMultipliers();
    if (M.empty()) return;
    std::vector<std::pair<int, int>> pr; std::vector<double> mv;
    for (const auto& e : M) { const auto& cn = G.conns[size_t(e.first)]; pr.push_back({std::min(cn.a, cn.b), std::max(cn.a, cn.b)}); mv.push_back(e.second); }
    sim.setConnectionMultipliersByCells(pr, mv);
    ++fmw.out().twoWayUpdates;
}
}  // namespace

bool runComp(const Case& c, const comp::DeckModel& m, std::vector<comp::CompSummary>& out, coupledbo::Output& co,
             const std::function<void(const std::string&)>& logFn, const std::function<bool(const comp::CompSummary&)>& progress,
             std::string* error, int maxSteps, comp::CompStats* stats,
             const std::function<void(const comp::CompSimulator&, const comp::CompSummary&)>& onState) {
    auto log = [&](const std::string& s) { if (logFn) logFn(s); };
    const double pu = pressureUnitPa(c);
    std::shared_ptr<Grid> G; std::unique_ptr<geoseq::Coupler> gs; std::unique_ptr<faultmech::Monitor> fm;
    const bool geoPerm = (c.geoSeq.plasticity == 1 && c.geoSeq.permPlastic > 0) || c.geoSeq.permStress > 0;
    int step = 0; bool setupDone = false;
    const bool ok = comp::runDeckModel(m, out, [&](const comp::CompSummary& S, comp::CompSimulator& sim) {
        const size_t ns = sim.pressure().size();
        if (!setupDone) {
            setupDone = true;
            std::string err;
            G = engineGrid(c, ns, err);
            if (!G) log("NOTE: coupled mechanics off - " + err);
            else {
                if (c.geoSeq.enabled) {
                    std::vector<double> phi(ns, 0.2);
                    for (size_t q = 0; q < ns && q < m.base.cellPV.size(); ++q) {
                        const int g = G->globalOf[q];
                        if (g >= 0 && size_t(g) < G->bulkVol.size() && G->bulkVol[size_t(g)] > 0) phi[q] = std::clamp(m.base.cellPV[q] / (G->bulkVol[size_t(g)] * (c.units == UnitSystem::Field ? 0.028316846592 : 1.0)), 0.01, 0.6);   // m3 vs ft3 (C35 fix)
                    }
                    gs = std::make_unique<geoseq::Coupler>();
                    if (!gs->setup(c.geoSeq, *G, G->globalOf, c.units == UnitSystem::Field ? 0.3048 : 1.0, sim.pressure(), phi, err)) { log("NOTE: FE geomechanics off - " + err); gs.reset(); }
                    else log(std::string("=== FE geomechanics (sequential, compositional)") + (gs->twoWay() ? ", two-way" : ", one-way"));
                }
                if (c.faultMech.on()) {
                    fm = std::make_unique<faultmech::Monitor>();
                    if (!fm->setup(c, *G, err)) { log("NOTE: fault mechanics off - " + err); fm.reset(); }
                    else {
                        log("=== fault mechanics (" + fm->out().engine + "): " + std::to_string(fm->out().faces.size()) + " faces" + (c.faultMech.slipOn ? ", quasi-static slip (FAULTSLIP)" : ""));
                        if (fm->wantsFeStress() && gs) { geoseq::Coupler* gp = gs.get(); fm->setStressProvider([gp](int cell, double S6[6]) { return gp->totalStress(cell, S6); }); }
                        else if (fm->wantsFeStress()) log("NOTE: STRESSFIELD ... FE needs GEOMECH - the analytic stress field is used");
                    }
                }
            }
        }
        if (gs) {
            if (!gs->step(sim.pressure())) { log("NOTE: FE geomechanics did not converge - coupling stopped"); if (fm) fm->setStressProvider({}); gs.reset(); }
            else if (gs->twoWay()) sim.setGeoFactors(gs->pvFactor(sim.pressure(), m.base.cr), geoPerm ? gs->permFactor() : std::vector<double>());
        }
        if (fm) {
            std::vector<float> pc(ns);
            for (size_t q = 0; q < ns; ++q) pc[q] = float(sim.pressure()[q] / pu);
            fm->step(S.time, pc);
            if (G) applyFaultMultipliers(*fm, *G, sim, *fm);
        }
        ++step;
        if (onState) onState(sim, S);
        return progress ? progress(S) : true;
    }, error, maxSteps, stats);
    co.finished = ok;
    if (fm) { fm->finish(nullptr, log); coupledbo::faultReport(*fm, co); }
    if (gs) {
        char b[200];
        std::snprintf(b, sizeof b, "  geomechanics: surface uplift %.4g mm, max failure index %.3f", gs->solver().surfaceUplift() * 1e3, gs->maxFailure());
        co.report.push_back(b);
        co.kpi.push_back({"surface_uplift_m", gs->solver().surfaceUplift()});
        co.kpi.push_back({"max_failure_index", gs->maxFailure()});
    }
    return ok;
}

bool runThermal(const Case& c, const thermal::DeckModel& m, std::vector<thermal::Summary>& out, coupledbo::Output& co,
                const std::function<void(const std::string&)>& logFn, const std::function<bool(const thermal::Summary&)>& progress,
                std::string* error, int maxSteps, thermal::Stats* stats,
                const std::function<void(const thermal::Simulator&, const thermal::Summary&)>& onState) {
    auto log = [&](const std::string& s) { if (logFn) logFn(s); };
    const double pu = pressureUnitPa(c);
    std::shared_ptr<Grid> G; std::unique_ptr<faultmech::Monitor> fm;
    const thermal::Simulator* simT = nullptr;
    bool setupDone = false;
    const bool ok = thermal::runDeckModel(m, out, [&](const thermal::Summary& S, thermal::Simulator& sim) {
        simT = &sim;
        const size_t ns = sim.pressure().size();
        if (!setupDone) {
            setupDone = true;
            std::string err;
            G = std::make_shared<Grid>();
            try { G->build(c); } catch (const std::exception& e) { err = e.what(); G.reset(); }
            // the thermal engine's cells are the active cells of the deck grid (+ EDFM / MINC cells appended at the end)
            if (G && size_t(G->numActive()) > ns) { err = "simulator cells do not match the grid"; G.reset(); }
            if (!G) log("NOTE: fault mechanics off - " + err);
            else if (c.faultMech.on()) {
                fm = std::make_unique<faultmech::Monitor>();
                if (!fm->setup(c, *G, err)) { log("NOTE: fault mechanics off - " + err); fm.reset(); }
                else {
                    log("=== fault mechanics (thermal): " + std::to_string(fm->out().faces.size()) + " faces" + (fm->out().thermalStress ? ", thermo-elastic stress" : "") +
                        (c.faultMech.slipOn ? ", quasi-static slip (FAULTSLIP)" : ""));
                    if (fm->wantsFeStress()) {
                        if (c.thermal.geo.enabled)
                            fm->setStressProvider([&simT](int cell, double S6[6]) {
                                const geomech::Solver* g = simT ? simT->mechanics() : nullptr;
                                if (!g) return false;
                                const auto& gm = simT->mechCellMap();
                                if (cell < 0 || size_t(cell) >= gm.size() || gm[size_t(cell)] < 0) return false;
                                g->totalStress(gm[size_t(cell)], S6);
                                return true;
                            });
                        else log("NOTE: STRESSFIELD ... FE needs the THM mechanics (GEOMECH) - the analytic stress field is used");
                    }
                }
            }
        }
        if (fm && G) {
            const size_t na = size_t(G->numActive());
            std::vector<float> p(na), T(na);
            for (size_t q = 0; q < na; ++q) { p[q] = float(sim.pressure()[q] / pu); T[q] = float(sim.temperature()[q]); }
            fm->step(S.time, p, &T);
            applyFaultMultipliers(*fm, *G, sim, *fm);
        }
        if (onState) onState(sim, S);
        return progress ? progress(S) : true;
    }, error, maxSteps, stats);
    simT = nullptr;
    co.finished = ok;
    if (fm) { fm->setStressProvider({}); fm->finish(nullptr, log); coupledbo::faultReport(*fm, co); }
    return ok;
}

}  // namespace coupleddeck
}  // namespace ressim
