// E300-style compositional PROPS keywords -> eos::CompFluid (SI units).
#include "ressim/deck_parser.h"
#include "ressim/eos.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>

namespace ressim {
namespace eos {

namespace {
std::string up(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
// all numeric items of record r (defaulted items -> NaN)
std::vector<double> nums(const deck::Record& r) {
    std::vector<double> v;
    for (size_t i = 0; i < r.items.size(); ++i) {
        if (r.items[i].def) { v.push_back(std::nan("")); continue; }
        try { v.push_back(r.num(i, 0.0)); } catch (...) { v.push_back(std::nan("")); }
    }
    return v;
}
}  // namespace

Component standardComponent(const std::string& name) {
    struct D { const char* n; double tc, pcBar, w, mw, vcCm3; };
    static const D t[] = {
        {"C1", 190.56, 45.99, 0.0115, 16.043, 98.6},  {"C2", 305.32, 48.72, 0.0995, 30.070, 145.5},
        {"C3", 369.83, 42.48, 0.1523, 44.097, 200.0}, {"IC4", 407.80, 36.40, 0.1844, 58.123, 262.7},
        {"NC4", 425.12, 37.96, 0.2002, 58.123, 255.0}, {"IC5", 460.40, 33.80, 0.2275, 72.150, 306.0},
        {"NC5", 469.70, 33.70, 0.2515, 72.150, 313.0}, {"C6", 507.60, 30.25, 0.3013, 86.177, 371.0},
        {"N2", 126.20, 33.98, 0.0377, 28.014, 90.1},   {"CO2", 304.13, 73.77, 0.2239, 44.010, 94.07},
        {"H2S", 373.10, 89.63, 0.0942, 34.080, 98.6},  {"H2O", 647.14, 220.64, 0.3443, 18.015, 55.95},
        {"H2", 33.19, 13.13, -0.216, 2.016, 64.1},
    };
    const std::string u = up(name);
    for (const D& d : t)
        if (u == d.n) {
            Component c; c.name = d.n; c.tc = d.tc; c.pc = d.pcBar * 1e5; c.omega = d.w; c.mw = d.mw * 1e-3; c.vc = d.vcCm3 * 1e-6;
            c.zc = c.pc * c.vc / (Rgas * c.tc);
            return c;
        }
    Component c; c.name = name;
    return c;
}

bool fluidFromDeck(const deck::RawDeck& d, bool field, CompFluid& out, std::vector<std::string>* msgs) {
    out = CompFluid();
    auto msg = [&](const std::string& s) { if (msgs) msgs->push_back(s); };
    int nc = 0;
    std::map<std::string, const deck::Keyword*> kw;
    for (const deck::Keyword& k : d.keywords) {
        kw[k.name] = &k;   // the last occurrence wins (E300 semantics for PROPS data)
        if ((k.name == "COMPS" || k.name == "NCOMPS") && !k.records.empty()) nc = std::max(nc, k.records[0].integer(0, 0));
        if (k.name == "CO2STORE") out.co2store = true;
        if (k.name == "CO2SOL") out.co2sol = true;
        if (k.name == "WATER") out.water = true;
        if (k.name == "PRCORR") out.pr78 = true;
    }
    if (nc <= 0 && out.co2store) {   // OPM-style CO2STORE deck (OIL/WATER + GAS, no COMPS): CO2 + H2O
        nc = 2;
        out.names = {"CO2", "H2O"};
    }
    if (nc <= 0) return false;
    out.present = true;
    if (kw.count("SALINITY") && !kw["SALINITY"]->records.empty() && kw["SALINITY"]->records[0].has(0))
        out.salinity = std::max(0.0, kw["SALINITY"]->records[0].num(0, 0));
    // unit factors to SI
    const double fT = field ? 5.0 / 9.0 : 1.0;                               // TCRIT
    const double fP = field ? 6894.757293168 : 1e5;                          // PCRIT
    const double fV = field ? 0.028316846592 / 453.59237 : 1e-3;             // VCRIT (ft3/lb-mol or m3/kg-mol)
    const double fD = field ? 0.3048 : 1.0;                                  // depth
    auto tempToK = [&](double t) { return field ? (t + 459.67) * 5.0 / 9.0 : t + 273.15; };
    // per-region arrays of nc values
    auto regionArrays = [&](const char* name) {
        std::vector<std::vector<double>> r;
        auto it = kw.find(name);
        if (it == kw.end()) return r;
        for (const deck::Record& rec : it->second->records) {
            if (rec.terminator) continue;
            std::vector<double> v = nums(rec);
            if (v.empty()) { if (!r.empty()) r.push_back(r.back()); continue; }   // empty record: copy the previous region
            if (int(v.size()) != nc && std::string(name) != "BIC")
                msg(std::string(name) + ": " + std::to_string(v.size()) + " values for " + std::to_string(nc) + " components");
            r.push_back(v);
        }
        return r;
    };
    // number of EOS regions = the largest record count among the property keywords
    const char* props[] = {"TCRIT", "PCRIT", "ACF", "MW", "VCRIT", "ZCRIT", "SSHIFT", "PARACHOR", "VCRITVIS", "ZCRITVIS", "OMEGAA", "OMEGAB", "BIC"};
    std::map<std::string, std::vector<std::vector<double>>> A;
    size_t nreg = 1;
    for (const char* p : props) { A[p] = regionArrays(p); nreg = std::max(nreg, A[p].size()); }
    std::vector<std::string> eosType;
    if (kw.count("EOS"))
        for (const deck::Record& r : kw["EOS"]->records) if (!r.terminator) eosType.push_back(up(r.str(0, "PR")));
    nreg = std::max(nreg, eosType.size());
    if (kw.count("CNAMES") && !kw["CNAMES"]->records.empty())
        for (const deck::Item& it : kw["CNAMES"]->records[0].items) out.names.push_back(it.s);
    while (int(out.names.size()) < nc) out.names.push_back("C" + std::to_string(out.names.size() + 1));
    out.names.resize(size_t(nc));
    auto pick = [&](const char* name, size_t reg, int i, double def) {
        const auto& v = A[name];
        if (v.empty()) return def;
        const auto& r = v[std::min(reg, v.size() - 1)];
        return (i < int(r.size()) && !std::isnan(r[size_t(i)])) ? r[size_t(i)] : def;
    };
    for (size_t reg = 0; reg < nreg; ++reg) {
        Model m;
        const std::string et = reg < eosType.size() ? eosType[reg] : eosType.empty() ? "PR" : eosType.back();
        if (et == "SRK") m.type = EosType::SRK;
        else if (et == "PR" || et.empty()) m.type = EosType::PR;
        else { msg("EOS " + et + " is not supported - Peng-Robinson is used"); m.type = EosType::PR; }
        m.pr78 = out.pr78;
        for (int i = 0; i < nc; ++i) {
            Component c;
            c.name = out.names[size_t(i)];
            const Component lib = standardComponent(c.name);
            c.tc = pick("TCRIT", reg, i, lib.tc / fT) * fT;
            c.pc = pick("PCRIT", reg, i, lib.pc / fP) * fP;
            c.omega = pick("ACF", reg, i, lib.omega);
            c.mw = pick("MW", reg, i, lib.mw * 1e3) * 1e-3;
            c.vc = pick("VCRIT", reg, i, 0.0) * fV;
            c.zc = pick("ZCRIT", reg, i, 0.0);
            if (c.vc <= 0 && c.zc > 0 && c.tc > 0) c.vc = c.zc * Rgas * c.tc / c.pc;
            if (c.zc <= 0 && c.vc > 0) c.zc = c.pc * c.vc / (Rgas * c.tc);
            if (c.vc <= 0 && lib.vc > 0) { c.vc = lib.vc; c.zc = c.pc * c.vc / (Rgas * c.tc); }
            c.sshift = pick("SSHIFT", reg, i, 0.0);
            c.parachor = pick("PARACHOR", reg, i, 0.0);
            c.vcVisc = pick("VCRITVIS", reg, i, 0.0) * fV;
            const double zcv = pick("ZCRITVIS", reg, i, 0.0);
            if (c.vcVisc <= 0 && zcv > 0) c.vcVisc = zcv * Rgas * c.tc / c.pc;
            c.omegaA = pick("OMEGAA", reg, i, -1.0);
            c.omegaB = pick("OMEGAB", reg, i, -1.0);
            if (!(c.tc > 0) || !(c.pc > 0) || !(c.mw > 0)) msg("component " + c.name + ": TCRIT / PCRIT / MW missing");
            m.comp.push_back(c);
        }
        // BIC: lower triangle by rows (k21; k31 k32; k41 k42 k43; ...)
        if (!A["BIC"].empty()) {
            const auto& b = A["BIC"][std::min(reg, A["BIC"].size() - 1)];
            size_t q = 0;
            for (int i = 1; i < nc; ++i)
                for (int j = 0; j < i; ++j, ++q)
                    if (q < b.size() && !std::isnan(b[q]) && b[q] != 0) m.setK(i, j, b[q]);
            if (b.size() != size_t(nc) * size_t(nc - 1) / 2) msg("BIC: " + std::to_string(b.size()) + " values, expected " + std::to_string(nc * (nc - 1) / 2));
        }
        if (kw.count("LBCCOEF") && !kw["LBCCOEF"]->records.empty()) {
            const std::vector<double> v = nums(kw["LBCCOEF"]->records[0]);
            for (size_t i = 0; i < 5 && i < v.size(); ++i) if (!std::isnan(v[i])) m.lbc[i] = v[i];
        }
        out.eosRegions.push_back(m);
    }
    // temperature
    for (const char* t : {"RTEMP", "RTEMPA"})
        if (kw.count(t) && !kw[t]->records.empty() && kw[t]->records[0].has(0)) out.resT = tempToK(kw[t]->records[0].num(0, 0));
    for (const char* t : {"TEMPVD", "RTEMPVD"})
        if (kw.count(t) && !kw[t]->records.empty()) {
            const std::vector<double> v = nums(kw[t]->records[0]);
            for (size_t i = 0; i + 1 < v.size(); i += 2) out.tempvd.push_back({v[i] * fD, tempToK(v[i + 1])});
            if (out.resT <= 0 && !out.tempvd.empty()) out.resT = out.tempvd.front().second;
        }
    // C33: OPM Flow's default reservoir temperature for CO2STORE / H2STORE decks without RTEMP / RTEMPVD is the standard
    // temperature 288.71 K (60 F)
    if (out.co2store && !(out.resT > 0)) out.resT = 288.71;
    // CO2SOL solubility table (SOLUBILI, table 1): P, RSW, BW, CW, VISC (5 columns; 4-column tables without CW accepted)
    if (kw.count("SOLUBILI") && !kw["SOLUBILI"]->records.empty()) {
        const std::vector<double> v = nums(kw["SOLUBILI"]->records[0]);
        const size_t nc5 = v.size() % 5 == 0 ? 5 : v.size() % 4 == 0 ? 4 : 5;
        const double fp = field ? 6894.757293168 : 1e5, frs = field ? 1000.0 / 5.614583333 : 1.0, fmu = 1e-3;
        for (size_t i = 0; i + nc5 <= v.size(); i += nc5) {
            CompFluid::SolRow r;
            r.p = v[i] * fp; r.rsw = v[i + 1] * frs; r.bw = v[i + 2];
            if (nc5 == 5) { r.cw = v[i + 3] / fp; r.muw = v[i + 4] * fmu; } else r.muw = v[i + 3] * fmu;
            out.solubili.push_back(r);
        }
        if (nc5 == 4) msg("SOLUBILI read with 4 columns (P, RSW, BW, VISC)");
    }
    // initial composition
    if (kw.count("ZMFVD"))
        for (const deck::Record& r : kw["ZMFVD"]->records) {
            if (r.terminator) continue;
            const std::vector<double> v = nums(r);
            ZmfTable t;
            for (size_t i = 0; i + size_t(nc) < v.size(); i += size_t(nc) + 1) {
                t.depth.push_back(v[i] * fD);
                std::vector<double> z(v.begin() + long(i) + 1, v.begin() + long(i) + 1 + nc);
                double s = 0; for (double q : z) s += q;
                if (s > 0) for (double& q : z) q /= s;
                t.z.push_back(z);
            }
            out.zmfvd.push_back(t);
        }
    if (kw.count("ZI") && !kw["ZI"]->records.empty()) {
        out.zi = nums(kw["ZI"]->records[0]);
        out.zi.resize(size_t(nc), 0.0);
        double s = 0; for (double q : out.zi) s += q;
        if (s > 0) for (double& q : out.zi) q /= s;
    }
    if (kw.count("STCOND") && !kw["STCOND"]->records.empty()) {
        const deck::Record& r = kw["STCOND"]->records[0];
        if (r.has(0)) out.stdT = tempToK(r.num(0, 0));
        if (r.has(1)) out.stdP = r.num(1, 0) * (field ? 6894.757293168 : 1e5);
    } else {
        out.stdT = field ? (60.0 + 459.67) * 5.0 / 9.0 : 288.71;
        out.stdP = field ? 14.6959 * 6894.757293168 : 1.01325e5;
    }
    return true;
}

}  // namespace eos
}  // namespace ressim
