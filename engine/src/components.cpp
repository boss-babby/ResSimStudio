// Multi-component product analysis: species recognition, standard heating values, produced-stream composition and
// product metrics (H2, syngas CO + H2, H2/CO ratio, LHV / HHV, CO2 / O2 slip), component cumulative totals.
#include "ressim/components.h"
#include "ressim/simulator.h"

#include <algorithm>
#include <cctype>
#include <cmath>

namespace ressim {

namespace {
struct SpeciesData { const char* formula; double mw; double hhv, lhv; };   // kg/mol, kJ/mol at 25 degC
// standard heats of combustion (NIST / GPSA 2172): HHV = products incl. liquid water, LHV = water vapour
const SpeciesData kSp[SP_COUNT] = {
    {"H2", 2.01588e-3, 285.83, 241.82},   {"CO", 28.0101e-3, 282.98, 282.98},  {"CO2", 44.0095e-3, 0, 0},
    {"CH4", 16.0425e-3, 890.36, 802.34},  {"C2H6", 30.069e-3, 1560.7, 1428.6}, {"C3H8", 44.0956e-3, 2219.2, 2043.1},
    {"C4H10", 58.1222e-3, 2877.5, 2657.3}, {"C5H12", 72.1488e-3, 3535.8, 3271.7}, {"C6H14", 86.1754e-3, 4194.8, 3886.7},
    {"N2", 28.0134e-3, 0, 0},            {"O2", 31.9988e-3, 0, 0},            {"H2S", 34.081e-3, 562.0, 518.0},
    {"H2O", 18.0153e-3, 0, 0},           {"NH3", 17.0305e-3, 382.8, 316.8},   {"AR", 39.948e-3, 0, 0}};

std::string norm(const std::string& s) {
    std::string r;
    for (char c : s) if (!std::isspace(static_cast<unsigned char>(c)) && c != '-' && c != '_' && c != '\'') r += char(std::toupper(static_cast<unsigned char>(c)));
    return r;
}
bool isLight(int sp) { return sp >= 0 && sp != SP_H2O && sp != SP_C5H12 && sp != SP_C6H14; }
} // namespace

const char* speciesFormula(int sp) { return sp >= 0 && sp < SP_COUNT ? kSp[sp].formula : ""; }

int recogniseSpecies(const std::string& name) {
    const std::string n = norm(name);
    if (n.empty()) return -1;
    static const std::pair<const char*, int> m[] = {
        {"H2", SP_H2}, {"HYDROGEN", SP_H2}, {"CO", SP_CO}, {"CARBONMONOXIDE", SP_CO}, {"CO2", SP_CO2}, {"CARBONDIOXIDE", SP_CO2},
        {"C1", SP_CH4}, {"CH4", SP_CH4}, {"METHANE", SP_CH4}, {"C2", SP_C2H6}, {"C2H6", SP_C2H6}, {"ETHANE", SP_C2H6},
        {"C3", SP_C3H8}, {"C3H8", SP_C3H8}, {"PROPANE", SP_C3H8}, {"C4", SP_C4H10}, {"IC4", SP_C4H10}, {"NC4", SP_C4H10},
        {"C4H10", SP_C4H10}, {"BUTANE", SP_C4H10}, {"C5", SP_C5H12}, {"IC5", SP_C5H12}, {"NC5", SP_C5H12}, {"C5H12", SP_C5H12},
        {"PENTANE", SP_C5H12}, {"C6", SP_C6H14}, {"NC6", SP_C6H14}, {"C6H14", SP_C6H14}, {"HEXANE", SP_C6H14},
        {"N2", SP_N2}, {"NITROGEN", SP_N2}, {"O2", SP_O2}, {"OXYGEN", SP_O2}, {"H2S", SP_H2S}, {"H2O", SP_H2O},
        {"WATER", SP_H2O}, {"STEAM", SP_H2O}, {"WAT", SP_H2O}, {"NH3", SP_NH3}, {"AMMONIA", SP_NH3}, {"AR", SP_AR}, {"ARGON", SP_AR}};
    for (const auto& e : m) if (n == e.first) return e.second;
    return -1;
}

void assignHeatingValues(std::vector<ComponentInfo>& comps) {
    for (ComponentInfo& c : comps) {
        c.hhv = c.lhv = 0;
        if (c.kind == 2) { c.species = -1; continue; }          // solids (coke, coal, fuel) are never produced
        c.species = c.kind == 3 ? SP_H2O : recogniseSpecies(c.name);
        if (c.mw <= 0 && c.species >= 0) c.mw = kSp[c.species].mw;
        if (c.species >= 0) { c.hhv = kSp[c.species].hhv * 1e3; c.lhv = kSp[c.species].lhv * 1e3; continue; }
        if (c.mw >= 0.072) { c.lhv = 44.0e6 * c.mw; c.hhv = 47.3e6 * c.mw; }   // C5+ pseudo components (paraffinic oil)
        else if (c.mw > 0 && c.kind == 0) {                                     // unknown light hydrocarbon: interpolate n-alkanes
            const double n = std::max(1.0, (c.mw * 1e3 - 2.016) / 14.027);     // CnH2n+2
            c.lhv = (802.34 + (n - 1) * 613.1) * 1e3; c.hhv = (890.36 + (n - 1) * 658.8) * 1e3;
        }
    }
}

StreamAnalysis analyseStream(const std::vector<ComponentInfo>& comps, const std::vector<double>& rates,
                             const std::vector<double>* inj) {
    StreamAnalysis A;
    const size_t nc = std::min(comps.size(), rates.size());
    A.yDry.assign(comps.size(), 0.0);
    A.yGas.assign(comps.size(), 0.0);
    std::vector<char> gasLike(comps.size(), 0);
    for (size_t i = 0; i < nc; ++i) {
        const ComponentInfo& c = comps[i];
        const double q = std::max(0.0, rates[i]);
        if (c.kind == 2 || q <= 0) continue;
        A.total += q;
        A.mass += q * c.mw;
        if (c.kind == 3 || c.species == SP_H2O) continue;
        A.dry += q;
        gasLike[i] = c.kind == 1 || isLight(c.species);
        if (gasLike[i]) A.gas += q;
        A.lhvRate += q * c.lhv; A.hhvRate += q * c.hhv;
        switch (c.species) {
            case SP_H2: A.h2 += q; break;
            case SP_CO: A.co += q; break;
            case SP_CO2: A.co2 += q; break;
            case SP_CH4: A.ch4 += q; break;
            case SP_O2: A.o2 += q; break;
            case SP_N2: A.n2 += q; break;
            case SP_H2S: A.h2s += q; break;
            default: break;
        }
    }
    for (size_t i = 0; i < nc; ++i) {
        const ComponentInfo& c = comps[i];
        const double q = std::max(0.0, rates[i]);
        if (c.kind == 2 || c.kind == 3 || c.species == SP_H2O) continue;
        if (A.dry > 0) A.yDry[i] = q / A.dry;
        if (A.gas > 0 && gasLike[i]) A.yGas[i] = q / A.gas;
    }
    A.syngas = A.h2 + A.co;
    A.h2co = A.co > 0 ? A.h2 / A.co : 0;
    A.h2Mass = A.h2 * kSp[SP_H2].mw;
    if (A.gas > 0) {   // heating value of the produced gas per standard m3 (gas-like components only)
        double l = 0, h = 0;
        for (size_t i = 0; i < nc; ++i) if (gasLike[i]) { const double q = std::max(0.0, rates[i]); l += q * comps[i].lhv; h += q * comps[i].hhv; }
        A.lhvGas = l / A.gas * kMolPerSm3 * 1e-6;
        A.hhvGas = h / A.gas * kMolPerSm3 * 1e-6;
    }
    if (inj) {
        double iCO2 = 0, iO2 = 0;
        for (size_t i = 0; i < comps.size() && i < inj->size(); ++i) {
            const double q = std::fabs((*inj)[i]);
            if (comps[i].species == SP_CO2) iCO2 += q;
            if (comps[i].species == SP_O2) iO2 += q;
        }
        A.co2Slip = iCO2 > 0 ? A.co2 / iCO2 : 0;
        A.o2Slip = iO2 > 0 ? A.o2 / iO2 : 0;
    }
    return A;
}

void accumulateComponentPoint(const SummaryPoint* prev, SummaryPoint& p) {
    const double dt = p.dt * 86400.0;
    auto acc = [&](const std::vector<double>* cum0, const std::vector<double>& rate, double sign, std::vector<double>& out) {
        out.assign(rate.size(), 0.0);
        for (size_t i = 0; i < rate.size(); ++i)
            out[i] = (cum0 && i < cum0->size() ? (*cum0)[i] : 0.0) + std::max(0.0, sign * rate[i]) * dt;
    };
    acc(prev ? &prev->fcompPT : nullptr, p.fcompP, 1.0, p.fcompPT);
    acc(prev ? &prev->fcompIT : nullptr, p.fcompI, 1.0, p.fcompIT);
    for (size_t w = 0; w < p.wells.size(); ++w) {
        const WellRates* q = prev && w < prev->wells.size() ? &prev->wells[w] : nullptr;
        acc(q ? &q->compCumP : nullptr, p.wells[w].comp, 1.0, p.wells[w].compCumP);
        acc(q ? &q->compCumI : nullptr, p.wells[w].comp, -1.0, p.wells[w].compCumI);
    }
}

void accumulateComponentTotals(std::vector<SummaryPoint>& summary) {
    std::vector<double> cp, ci;
    std::vector<std::vector<double>> wp, wi;
    for (SummaryPoint& s : summary) {
        const double dt = s.dt * 86400.0;   // summary dt in days, rates in mol/s
        auto acc = [&](std::vector<double>& cum, const std::vector<double>& rate, std::vector<double>& out) {
            if (cum.size() < rate.size()) cum.resize(rate.size(), 0.0);
            for (size_t i = 0; i < rate.size(); ++i) cum[i] += std::max(0.0, rate[i]) * dt;
            out = cum;
        };
        acc(cp, s.fcompP, s.fcompPT);
        acc(ci, s.fcompI, s.fcompIT);
        if (wp.size() < s.wells.size()) { wp.resize(s.wells.size()); wi.resize(s.wells.size()); }
        for (size_t w = 0; w < s.wells.size(); ++w) {
            WellRates& r = s.wells[w];
            std::vector<double> pr(r.comp.size()), ir(r.comp.size());
            for (size_t i = 0; i < r.comp.size(); ++i) { pr[i] = std::max(0.0, r.comp[i]); ir[i] = std::max(0.0, -r.comp[i]); }
            acc(wp[w], pr, r.compCumP);
            acc(wi[w], ir, r.compCumI);
        }
    }
}

} // namespace ressim
