// C27: headless coupled black-oil runner (see coupled_bo.h) - formerly runCoupledBO in tools/ressim_cli.cpp (C26)
#include "ressim/coupled_bo.h"
#include "ressim/faultmech.h"
#include "ressim/geoseq.h"
#include "ressim/hf_ddm.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <memory>

namespace ressim {
namespace coupledbo {

bool needed(const Case& c) {
    if (c.thermal.enabled || c.comp.enabled) return false;   // the own comp / thermal engines have their own couplings
    bool coupled = c.faultMech.on() || c.geoSeq.enabled;
    for (const auto& F : c.grid.hfracs) coupled |= F.ddmWidth != 0;
    return coupled;
}

double kpiSum(const std::vector<std::pair<std::string, double>>& kpi, const std::string& suffix) {
    double s = 0; bool any = false;
    for (const auto& k : kpi)
        if (k.first.size() >= suffix.size() && k.first.compare(k.first.size() - suffix.size(), suffix.size(), suffix) == 0) { s += k.second; any = true; }
    return any ? s : std::numeric_limits<double>::quiet_NaN();
}
double kpiMax(const std::vector<std::pair<std::string, double>>& kpi, const std::string& suffix) {
    double s = -1e300; bool any = false;
    for (const auto& k : kpi)
        if (k.first.size() >= suffix.size() && k.first.compare(k.first.size() - suffix.size(), suffix.size(), suffix) == 0) { s = std::max(s, k.second); any = true; }
    return any ? s : std::numeric_limits<double>::quiet_NaN();
}

void faultReport(const faultmech::Monitor& fm, Output& out) {
    char b[400];
    {
        const auto& O = fm.out();
        for (size_t f = 0; f < O.series.size(); ++f) if (!O.series[f].t.empty()) {
            const auto& S = O.series[f];
            const double mx = *std::max_element(S.maxDcfs.begin(), S.maxDcfs.end());
            std::string line;
            std::snprintf(b, sizeof b, "  fault %s: max dCFS %.8g, expected events %.10g, catalogue %.0f, Mw max %.2f", O.names[f].c_str(), mx, S.expected.back(), S.events.back(), S.maxMw.back());
            line = b;
            if (!S.maxSlip.empty()) { std::snprintf(b, sizeof b, ", max slip %.4g mm, M0 %.4g N m", S.maxSlip.back() * 1e3, S.moment.empty() ? 0.0 : S.moment.back()); line += b; }
            out.report.push_back(line);
            const std::string p = "fault_" + O.names[f];
            out.kpi.push_back({p + "_expected_events", S.expected.back()});
            out.kpi.push_back({p + "_max_dcfs", mx});
            out.kpi.push_back({p + "_catalogue_events", S.events.back()});
            if (!S.maxMw.empty()) out.kpi.push_back({p + "_max_mw", S.maxMw.back()});
            if (!S.maxSlip.empty()) { out.kpi.push_back({p + "_max_slip_m", S.maxSlip.back()}); out.kpi.push_back({p + "_moment_Nm", S.moment.empty() ? 0.0 : S.moment.back()}); }
        }
    }
}

bool run(Simulator& sim, const Case& c, Output& out, const std::function<void(const std::string&)>& logFn,
         std::atomic<bool>* cancel, const std::function<void(double)>& progress) {
    auto log = [&](const std::string& m) { if (logFn) logFn(m); };
    char b[400];
    sim.cancel = cancel;
    const double paU = c.units == UnitSystem::Field ? 6894.757 : 1e5;
    std::unique_ptr<faultmech::Monitor> fm; std::unique_ptr<geoseq::Coupler> gs; std::unique_ptr<hfddm::Aperture> hfd;
    std::string err;
    if (c.faultMech.on()) {
        fm = std::make_unique<faultmech::Monitor>();
        if (!fm->setup(c, sim.grid(), err)) { log("NOTE: fault mechanics off - " + err); fm.reset(); }
        else log("=== fault mechanics: " + std::to_string(fm->out().faces.size()) + " faces" + (c.faultMech.slipOn ? ", quasi-static slip (FAULTSLIP)" : ""));
    }
    bool anyFeHf = false; for (const auto& F : c.grid.hfracs) if (F.kic > 0 && F.feStress) anyFeHf = true;
    const bool geoPerm = (c.geoSeq.plasticity == 1 && c.geoSeq.permPlastic > 0) || c.geoSeq.permStress > 0;
    if (c.geoSeq.enabled) {
        const Grid& G = sim.grid(); const size_t na = size_t(G.numActive());
        const auto& P0 = sim.snapshots().back().pressure;
        std::vector<double> pPa(na, 0.0), phi(na, 0.2);
        for (size_t a = 0; a < na; ++a) {
            pPa[a] = a < P0.size() ? double(P0[a]) * paU : 0.0;
            const int gq = G.globalOf[a];
            if (gq >= 0 && size_t(gq) < G.bulkVol.size() && G.bulkVol[size_t(gq)] > 0) phi[a] = std::clamp(G.poreVolRef[a] * (c.units == UnitSystem::Field ? 5.614583 : 1.0) / G.bulkVol[size_t(gq)], 0.01, 0.6);   // rb vs ft3 (C35 fix)
        }
        gs = std::make_unique<geoseq::Coupler>();
        if (!gs->setup(c.geoSeq, G, G.globalOf, c.units == UnitSystem::Field ? 0.3048 : 1.0, pPa, phi, err)) { log("NOTE: FE geomechanics off - " + err); gs.reset(); }
        else log("=== FE geomechanics (sequential)");
    }
    if (fm && fm->wantsFeStress() && gs) { geoseq::Coupler* gp = gs.get(); fm->setStressProvider([gp](int cell, double S[6]) { return gp->totalStress(cell, S); }); }
    bool wantDdm = false; for (const auto& F : c.grid.hfracs) if (F.ddmWidth) wantDdm = true;
    if (wantDdm) {
        hfd = std::make_unique<hfddm::Aperture>();
        if (!hfd->setup(sim, err)) { log("NOTE: DDMWIDTH off - " + err); hfd.reset(); }
        else log("=== DDM apertures: " + std::to_string(hfd->out().elements) + " pieces");
    }
    std::vector<int> touched;
    const bool gsIter = gs && c.geoSeq.coupling == 2;   // C26: iterative two-way coupling (report step repeated to a fixed point)
    std::vector<double> pvApplied; bool gsStepped = false;
    auto coupleStep = [&]() {
        const auto& P = sim.snapshots().back().pressure;
        if (gs && gsStepped) {
            gsStepped = false;
            if (anyFeHf) { std::vector<double> sm = gs->minStress(); for (double& v : sm) v /= paU; sim.setMinStress(sm); }
        } else if (gs) {
            std::vector<double> pPa(P.size()); for (size_t a = 0; a < P.size(); ++a) pPa[a] = double(P[a]) * paU;
            if (!gs->step(pPa)) { log("NOTE: FE geomechanics did not converge - coupling stopped"); if (fm) fm->setStressProvider({}); gs.reset(); }
            else {
                if (gs->twoWay()) { pvApplied = gs->pvFactor(pPa, c.props.rock.cr / paU); sim.setPoreVolumeMultipliers(pvApplied); }
                if (anyFeHf) { std::vector<double> sm = gs->minStress(); for (double& v : sm) v /= paU; sim.setMinStress(sm); }
            }
        }
        if (fm) fm->step(sim.time(), P);
        // composed connection multipliers (fault reactivation x geomechanics permeability)
        std::map<int, double> M;
        const Grid& G = sim.grid();
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
        if (fm) for (const auto& e : fm->connMultipliers()) { auto it = M.find(e.first); if (it == M.end()) M[e.first] = e.second; else it->second *= e.second; }
        for (int q : touched) if (!M.count(q)) M[q] = 1.0;
        if (!M.empty()) {
            std::vector<int> cs; std::vector<double> ms; for (const auto& e : M) { cs.push_back(e.first); ms.push_back(e.second); }
            sim.setConnectionMultipliers(cs, ms);
            touched.clear(); for (const auto& e : M) if (e.second != 1.0) touched.push_back(e.first);
        }
        if (hfd) hfd->update(sim);
    };
    coupleStep();
    auto reportStep = [&]() -> bool {
        if (!gsIter || !gs) return sim.runReportStep();
        const geoseq::IterStep R = geoseq::iterativeReportStep(sim, *gs, pvApplied, paU, c.props.rock.cr / paU, c.geoSeq.maxIter, c.geoSeq.tol);
        if (!R.ran) return false;
        if (R.feOk) { gsStepped = true; ++out.iterSteps; out.iterPasses += R.passes; }
        std::snprintf(b, sizeof b, "  iterative geomechanics: %d pass(es), max pore-volume change %.2e%s", R.passes, R.diff, R.converged ? "" : " (maxIter)");
        log(b);
        return true;
    };
    while (reportStep()) {
        coupleStep();
        if (progress && !sim.summary().empty()) progress(sim.summary().back().time);
    }
    out.finished = sim.finished();
    if (out.iterSteps > 0) {
        std::snprintf(b, sizeof b, "=== Iterative geomechanics: %d report steps, %.2f flow passes per step", out.iterSteps, double(out.iterPasses) / out.iterSteps);
        log(b);
    }
    if (fm) fm->finish(nullptr, log);
    if (fm) faultReport(*fm, out);
    if (gs) {
        std::snprintf(b, sizeof b, "  geomechanics: surface uplift %.4g mm, max failure index %.3f", gs->solver().surfaceUplift() * 1e3, gs->maxFailure());
        out.report.push_back(b);
        out.kpi.push_back({"surface_uplift_m", gs->solver().surfaceUplift()});
        out.kpi.push_back({"max_failure_index", gs->maxFailure()});
    }
    if (hfd) {
        const auto& o = hfd->out();
        std::snprintf(b, sizeof b, "  DDM apertures: %d / %d open, max %.4g mm, volume %.4g m3", o.open, o.elements, o.maxWidth * 1e3, o.volume);
        out.report.push_back(b);
        out.kpi.push_back({"hf_max_aperture_m", o.maxWidth}); out.kpi.push_back({"hf_opening_volume_m3", o.volume});
    }
    return out.finished;
}

}  // namespace coupledbo
}  // namespace ressim
