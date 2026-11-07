// C48: shared per-report-step cell maps of the compositional / CO2STORE and thermal engines (see cell_snapshot.h).
// Moved verbatim from the GUI runner (app/src/sim_runner.cpp, C13-C25) so that the CLI writes the same arrays.
#include "ressim/cell_snapshot.h"
#include "ressim/geomech.h"
#include <algorithm>
#include <string>

namespace ressim {
namespace comp {

Snapshot cellSnapshot(const CompSimulator& sim, double time, double pu) {
    const CompCase& SM = sim.model();
    const bool co2 = SM.co2store;
    Snapshot sn; sn.time = time;
    const auto& P = sim.pressure();
    const std::vector<double> sg = sim.gasSaturation(), sw = sim.waterSaturation();
    const size_t n = P.size();
    sn.pressure.resize(n); sn.sw.resize(n); sn.sg.resize(n); sn.so.resize(n); sn.rs.assign(n, 0.0f);
    for (size_t q = 0; q < n; ++q) {
        sn.pressure[q] = float(P[q] / pu); sn.sw[q] = float(sw[q]); sn.sg[q] = float(sg[q]); sn.so[q] = float(std::max(0.0, 1.0 - sw[q] - sg[q]));
    }
    if (co2) {   // brine saturation = 1 - Sg; dissolved CO2 mole fraction of the brine
        const std::vector<double> x = sim.dissolvedCO2();
        sn.xco2.resize(n);
        for (size_t q = 0; q < n; ++q) { sn.sw[q] = float(1.0 - sg[q]); sn.so[q] = 0.0f; sn.xco2[q] = float(q < x.size() ? x[q] : 0.0); }
    } else {     // overall mole fraction of component 1
        const auto& N = sim.moles(); const int ncw = sim.ncw(), ncc = sim.nc();
        sn.xco2.resize(n);
        for (size_t q = 0; q < n; ++q) { double t = 0; for (int i = 0; i < ncc; ++i) t += N[q * size_t(ncw) + size_t(i)]; sn.xco2[q] = float(t > 0 ? N[q * size_t(ncw)] / t : 0.0); }
    }
    auto put = [&](const char* name, const std::vector<double>& v) {
        std::vector<float> f(n, 0.0f);
        for (size_t q = 0; q < n && q < v.size(); ++q) f[q] = float(v[q]);
        sn.extra.push_back({name, std::move(f)});
    };
    if (SM.hyst.on) { put("SGMAX", sim.sgMax()); put("SGTRAP", sim.trappedGasSaturation()); }
    else put("SGTRAP", sim.trappedGasSaturation());
    if (co2 && SM.trap.saltPrecip) put("SSALT", sim.saltSaturation());
    if (co2 && !SM.trap.minerals.empty()) put("MINCO2", sim.mineralCO2());
    if (co2 && (SM.trap.saltPrecip || !SM.trap.minerals.empty())) { put("PORFAC", sim.porosityFactor()); put("PERMFAC", sim.permFactorCells()); }
    const int ncc = sim.nc();
    if (ncc <= 16) {
        std::vector<double> X, Y, Z; sim.compositions(X, Y, Z);
        const char* nm[3] = {"ZMF", "XMF", "YMF"};
        const std::vector<double>* src[3] = {&Z, &X, &Y};
        for (int a = 0; a < 3; ++a)
            for (int i = 0; i < ncc; ++i) {
                std::vector<float> f(n, 0.0f);
                for (size_t q = 0; q < n && q * size_t(ncc) + size_t(i) < src[a]->size(); ++q) f[q] = float((*src[a])[q * size_t(ncc) + size_t(i)]);
                sn.extra.push_back({std::string(nm[a]) + std::to_string(i + 1), std::move(f)});
            }
    }
    // CO2SOL: CO2 as an EOS component that also dissolves in the brine (E300 AMF / RSW style arrays)
    if (!co2 && SM.co2sol && SM.co2Comp >= 0 && SM.co2Comp < ncc) {
        const std::vector<double> xa = sim.dissolvedCO2();
        for (int i = 0; i < ncc; ++i) {
            std::vector<float> f(n, 0.0f);
            if (i == SM.co2Comp) for (size_t q = 0; q < n && q < xa.size(); ++q) f[q] = float(xa[q]);
            sn.extra.push_back({"AMF" + std::to_string(i + 1), std::move(f)});
        }
        put("XCO2AQ", xa);
        const auto& N = sim.moles(); const int ncw = sim.ncw();   // overall CO2 mole fraction (hydrocarbon + aqueous) in the XCO2 slot
        for (size_t q = 0; q < n; ++q) { double t = 0; for (int i = 0; i < ncc; ++i) t += N[q * size_t(ncw) + size_t(i)];
                                         sn.xco2[q] = float(t > 0 ? N[q * size_t(ncw) + size_t(SM.co2Comp)] / t : 0.0); }
    }
    return sn;
}

}  // namespace comp

namespace thermal {

Snapshot cellSnapshot(const Simulator& sim, const ThermalInput& TI, double time, size_t n, bool fld, std::vector<float>& T0map) {
    const double pu = fld ? 6894.757293168 : 1e5;
    auto tOut = [&](double K) { return fld ? K * 1.8 - 459.67 : K - 273.15; };
    Snapshot sn; sn.time = time;
    sn.pressure.assign(n, 0.0f); sn.sw.assign(n, 0.0f); sn.sg.assign(n, 0.0f); sn.so.assign(n, 0.0f); sn.rs.assign(n, 0.0f);
    std::vector<float> T(n, 0.0f);
    const int nc = sim.nc(), nb = nc + 3;
    std::vector<std::vector<float>> z(size_t(nc), std::vector<float>(n, 0.0f));
    const auto& X = sim.moles();
    for (int q = 0; q < sim.size() && size_t(q) < n; ++q) {
        sn.pressure[size_t(q)] = float(sim.pressure()[size_t(q)] / pu);
        sn.sw[size_t(q)] = float(sim.sw()[size_t(q)]); sn.so[size_t(q)] = float(sim.so()[size_t(q)]); sn.sg[size_t(q)] = float(sim.sg()[size_t(q)]);
        T[size_t(q)] = float(tOut(sim.temperature()[size_t(q)]));
        double tot = 0; for (int i = 0; i < nc; ++i) if (sim.model().comps[size_t(i)].type != 2) tot += X[size_t(q) * size_t(nb) + 1 + size_t(i)];
        for (int i = 0; i < nc; ++i) if (sim.model().comps[size_t(i)].type != 2) z[size_t(i)][size_t(q)] = float(tot > 0 ? X[size_t(q) * size_t(nb) + 1 + size_t(i)] / tot : 0.0);
    }
    {   // C13: cooling map T - T(t = 0)
        if (T0map.empty()) T0map = T;
        std::vector<float> dT(n, 0.0f);
        for (size_t q = 0; q < n && q < T0map.size(); ++q) dT[q] = T[q] - T0map[q];
        sn.extra.push_back({"DTEMP", std::move(dT)});
    }
    sn.extra.push_back({"TEMP", std::move(T)});
    {   // steam zone indicator: Sg where the cell is at the water saturation temperature
        std::vector<float> st(n, 0.0f);
        for (size_t q = 0; q < n && q < sim.sg().size(); ++q) st[q] = sn.sg[q];
        sn.extra.push_back({"SSTEAM", std::move(st)});
    }
    for (int i = 0; i < nc && nc <= 16; ++i) sn.extra.push_back({"ZMF" + std::to_string(i + 1), std::move(z[size_t(i)])});
    if (nc <= 16) {   // phase compositions (index k = component k, nc + 1 = water): oil XMF, gas YMF, aqueous AMF
        std::vector<double> xo, yg, aq; sim.phaseCompositions(xo, yg, aq);
        const size_t nk = size_t(nc) + 1, ns = std::min(n, size_t(sim.size()));
        const bool anyAq = [&]() { for (const auto& q : sim.model().comps) if (q.soluble()) return true; return false; }();
        const char* nm[3] = {"XMF", "YMF", "AMF"};
        const std::vector<double>* src[3] = {&xo, &yg, &aq};
        for (int a = 0; a < 3; ++a) {
            if (a == 2 && !anyAq) continue;
            for (size_t k = 0; k < nk; ++k) {
                if (a == 0 && k == nk - 1) continue;   // no water in the oil phase
                std::vector<float> f(n, 0.0f);
                for (size_t q = 0; q < ns && q * nk + k < src[a]->size(); ++q) f[q] = float((*src[a])[q * nk + k]);
                sn.extra.push_back({std::string(nm[a]) + std::to_string(k + 1), std::move(f)});
            }
        }
        const std::vector<double> st = sim.trappedGasSaturation();
        std::vector<float> f(n, 0.0f); for (size_t q = 0; q < ns && q < st.size(); ++q) f[q] = float(st[q]);
        sn.extra.push_back({"SGTRAP", std::move(f)});
        if (sim.model().hyst.on) {
            const std::vector<double>& sm = sim.sgMax();
            std::vector<float> g(n, 0.0f); for (size_t q = 0; q < ns && q < sm.size(); ++q) g[q] = float(sm[q]);
            sn.extra.push_back({"SGMAX", std::move(g)});
        }
        const int tk = sim.model().track;
        if (tk >= 0) {   // XCO2 slot: overall mole fraction of the tracked component (incl. water in the total)
            sn.xco2.assign(n, 0.0f);
            for (size_t q = 0; q < ns; ++q) {
                double t = 0; for (int i = 0; i <= nc; ++i) if (i == nc || sim.model().comps[size_t(i)].type != 2) t += X[q * size_t(nb) + 1 + size_t(i)];
                sn.xco2[q] = float(t > 0 ? X[q * size_t(nb) + 1 + size_t(tk)] / t : 0.0);
            }
        }
    }
    {   // THMC maps: fluid porosity, permeability multiplier, solids, mechanics
        const size_t ns = std::min(n, size_t(sim.size()));
        if (!TI.reactions.empty() || TI.permLaw != 0 || TI.geo.enabled) {
            const std::vector<double> pf = sim.fluidPorosity();
            std::vector<float> a(n, 0.0f), km(n, 1.0f);
            for (size_t q = 0; q < ns && q < pf.size(); ++q) a[q] = float(pf[q]);
            for (size_t q = 0; q < ns && q < sim.permMultiplier().size(); ++q) km[q] = float(sim.permMultiplier()[q]);
            sn.extra.push_back({"POROF", std::move(a)});
            sn.extra.push_back({"PERMMULT", std::move(km)});
        }
        for (int i = 0; i < nc && i < int(TI.comps.size()); ++i) {
            if (TI.comps[size_t(i)].type != 2) continue;   // solid concentration, kg per m3 of bulk volume
            std::vector<float> a(n, 0.0f);
            for (size_t q = 0; q < ns; ++q) a[q] = float(X[q * size_t(nb) + 1 + size_t(i)] * TI.comps[size_t(i)].mw * 1e-3 / std::max(sim.bulkVolume(int(q)), 1e-12));
            sn.extra.push_back({"SOLID_" + TI.comps[size_t(i)].name, std::move(a)});
        }
        if (!sim.model().fracCells.empty()) {   // C12 EDFM: fracture aperture (mm) / pressure / temperature on the host cells
            const std::vector<double> wf = sim.fractureApertures();
            std::vector<float> ap(n, 0.0f), fp(n, 0.0f), ft(n, 0.0f);
            const auto& FC = sim.model().fracCells;
            for (size_t f = 0; f < FC.size() && f < wf.size(); ++f) {
                const size_t h = size_t(FC[f].host), c = size_t(FC[f].cell);
                if (FC[f].host < 0 || h >= n) continue;
                if (float(wf[f] * 1e3) > ap[h]) {
                    ap[h] = float(wf[f] * 1e3); fp[h] = float(sim.pressure()[c] / pu); ft[h] = float(tOut(sim.temperature()[c]));
                }
            }
            sn.extra.push_back({"FRACAP", std::move(ap)});
            sn.extra.push_back({"FRACP", std::move(fp)});
            sn.extra.push_back({"FRACT", std::move(ft)});
        }
        if (const geomech::Solver* g = sim.mechanics()) {
            const std::vector<int>& gm = sim.mechCellMap();
            auto mapv = [&](const std::vector<double>& v, double sc) {
                std::vector<float> a(n, 0.0f);
                for (size_t q = 0; q < ns && q < gm.size(); ++q) { const size_t k = size_t(gm[q]); if (k < v.size()) a[q] = float(v[k] * sc); }
                return a;
            };
            sn.extra.push_back({"UPLIFT", mapv(g->dispZ(), -1.0)});            // m, + up
            sn.extra.push_back({"DISPX", mapv(g->dispX(), 1.0)});
            sn.extra.push_back({"STRAINV", mapv(g->volStrain(), 1.0)});
            sn.extra.push_back({"SIGEFFM", mapv(g->sigEffMean(), 1.0 / pu)});   // pressure units
            sn.extra.push_back({"SIGEFFV", mapv(g->sigEffV(), 1.0 / pu)});
            sn.extra.push_back({"SIGEFFH", mapv(g->sigEffH(), 1.0 / pu)});
            sn.extra.push_back({"FAILURE", mapv(g->failure(), 1.0)});
            sn.extra.push_back({"POROM", mapv(g->porosity(), 1.0)});
            {   // C25 smallest principal total stress (pressure units)
                std::vector<float> smin(n, 0.0f);
                for (size_t q = 0; q < ns && q < gm.size(); ++q) if (gm[q] >= 0) smin[q] = float(g->minTotalStress(gm[q]) / pu);
                sn.extra.push_back({"SIGMIN", std::move(smin)});
            }
            if (sim.model().geo.plasticity == 1) {   // C16 elasto-plasticity
                sn.extra.push_back({"EPSPLAS", mapv(g->plasticStrain(), 1.0)});
                sn.extra.push_back({"EPSPV", mapv(g->plasticVolStrain(), 1.0)});
                sn.extra.push_back({"KDAMAGE", mapv(g->permMultiplier(), 1.0)});
            }
        }
    }
    return sn;
}

}  // namespace thermal
}  // namespace ressim
