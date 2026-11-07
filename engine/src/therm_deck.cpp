// Thermal deck hook-up (milestone C7) - see therm_deck.h.
#include "ressim/binio.h"
#include "ressim/gasdata.h"
#include "ressim/therm_deck.h"
#include "ressim/geoseq.h"
#include "ressim/deck_parser.h"
#include "ressim/grid.h"
#include "ressim/io.h"
#include "ressim/steam.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace ressim {

namespace {
constexpr double GRAV = 9.80665;
const double NaN = std::numeric_limits<double>::quiet_NaN();

// deck unit system -> SI
struct TU {
    bool field;
    double len, p, liq, gas, dens, trans, comp, pv, absT, spech, hv, heatcr, thcon;
    explicit TU(bool f) : field(f) {
        if (f) {
            len = 0.3048; p = 6894.757293168; liq = 0.158987294928; gas = 28.316846592; dens = 16.01846337;
            trans = 0.158987294928 * 1e-3 / (86400.0 * 6894.757293168); comp = 1.0 / 6894.757293168; pv = 0.158987294928;
            absT = 5.0 / 9.0; spech = 4186.8; hv = 2326.0; heatcr = 67066.1; thcon = 1055.05585 / (0.3048 * 86400.0 * 5.0 / 9.0);
        } else {
            len = 1; p = 1e5; liq = 1; gas = 1; dens = 1; trans = 1e-3 / (86400.0 * 1e5); comp = 1e-5; pv = 1;
            absT = 1; spech = 1000; hv = 1000; heatcr = 1000; thcon = 1000.0 / 86400.0;
        }
    }
    double tempK(double t) const { return field ? (t + 459.67) * 5.0 / 9.0 : t + 273.15; }
    double tempDeck(double T) const { return field ? T * 1.8 - 459.67 : T - 273.15; }
    double energy() const { return field ? 1055.05585 : 1000.0; }   // J per kJ | Btu
    double stdT() const { return field ? 60.0 : 15.56; }
    double stdP() const { return field ? 14.696 : 1.01325; }
};

std::string up(std::string s) { for (char& ch : s) ch = char(std::toupper(static_cast<unsigned char>(ch))); return s; }
double lin(const std::vector<double>& x, const std::vector<double>& y, double v) {
    if (x.empty()) return 0;
    if (v <= x.front()) return y.front();
    if (v >= x.back()) return y.back();
    const size_t i = size_t(std::upper_bound(x.begin(), x.end(), v) - x.begin());
    const double t = (v - x[i - 1]) / std::max(1e-300, x[i] - x[i - 1]);
    return y[i - 1] + t * (y[i] - y[i - 1]);
}
std::string num(double v) {
    char b[40];
    for (int p = 8; p <= 17; ++p) { std::snprintf(b, sizeof b, "%.*g", p, v); if (std::strtod(b, nullptr) == v) break; }
    return b;
}
bool matchWell(const std::string& pat, const std::string& well) {
    const std::string P = up(pat), W = up(well);
    if (!P.empty() && P.back() == '*') return W.compare(0, P.size() - 1, P, 0, P.size() - 1) == 0;
    return P == W;
}
// all items of the first data record (NaN = defaulted)
std::vector<double> recordValues(const deck::Keyword& k) {
    std::vector<double> v;
    for (const deck::Record& r : k.records) {
        if (r.terminator) continue;
        for (size_t i = 0; i < r.items.size(); ++i) v.push_back(r.has(i) && r.isNumber(i) ? r.num(i, 0.0) : NaN);
        break;
    }
    return v;
}
// least-squares fit mu = A T^B on a table (T > 0, mu > 0): returns false when fewer than 2 usable rows
bool fitPower(const std::vector<double>& T, const std::vector<double>& mu, double& A, double& B) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0; int n = 0;
    for (size_t i = 0; i < T.size(); ++i) {
        if (!(T[i] > 0 && mu[i] > 0)) continue;
        const double x = std::log(T[i]), y = std::log(mu[i]);
        sx += x; sy += y; sxx += x * x; sxy += x * y; ++n;
    }
    if (n < 2 || std::fabs(n * sxx - sx * sx) < 1e-300) return false;
    B = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    A = std::exp((sy - B * sx) / n);
    return true;
}
}  // namespace

const std::set<std::string>& thermalKeywords() {
    static const std::set<std::string> s = {"THERMAL", "DEADOIL", "COMPS", "CNAMES", "MW", "TCRIT", "PCRIT", "ACF", "CVTYPE", "KVCR", "KVWI",
        "DREF", "CREF", "THERMEX1", "SPECHA", "SPECHB", "SPECHG", "SPECHH", "HEATVAPS", "HEATVAP", "HEATVAPE", "ZFACTOR", "ZFACT1",
        "OILVISCT", "GASVISCT", "GASVISCF", "WATVISCT", "STCOND", "THSTT97", "TEMPVD", "TEMPI", "ZI", "XMF", "YMF", "WINJTEMP",
        "HEATCR", "HEATCRT", "THCONR", "THCONSF", "PREF", "TREF"};
    return s;
}

// ======================================================================== deck -> Case::thermal
bool parseGeoRecord(const std::string& key, const deck::Record& r, geomech::Input& G) {
    auto U = [](std::string v) { for (char& ch : v) ch = char(std::toupper((unsigned char)ch)); return v; };
    if (key == "GEOMECH") {
        G.enabled = U(r.str(1, "ON")) != "OFF";
        G.coupling = int(r.num(2, 1.0));
        G.E = {r.num(3, 5e9)}; G.nu = {r.num(4, 0.25)}; G.biot = {r.num(5, 1.0)}; G.alphaT = {r.num(6, 1e-5)};
        G.Ks = r.num(7, 0.0);
        G.overLayers = int(r.num(8, 0.0)); G.underLayers = int(r.num(9, 0.0)); G.sidePad = int(r.num(10, 0.0));
        G.overThick = r.num(11, 300.0); G.underThick = r.num(12, 300.0); G.sideThick = r.num(13, 300.0);
        G.padE = r.num(14, 5e9); G.padNu = r.num(15, 0.25); G.padAlphaT = r.num(16, 1e-5);
        G.topFree = r.num(17, 1.0) != 0; G.sigmaVGrad = r.num(18, 22.6e3); G.k0 = r.num(19, 0.8);
        G.cohesion = r.num(20, 2e6); G.frictionDeg = r.num(21, 30.0); G.tensile = r.num(22, 1e6);
        G.maxIter = int(r.num(23, 4.0)); G.tol = r.num(24, 1e-4);
        G.permStress = r.num(25, 0.0);   // C25 stress-dependent permeability (1/Pa)
        return true;
    }
    if (key == "MECHBC") {   // MECHBC face(0 x- .. 5 base) type(-1 default, 0 roller, 1 free, 2 fixed, 3 plate) stress(Pa) [disp dispControl]
        const int f = int(r.num(1, -1.0));
        if (f >= 0 && f < 6) {
            G.bcType[f] = int(r.num(2, -1.0)); G.bcStress[f] = r.num(3, 0.0);
            G.bcDisp[f] = r.num(4, 0.0); G.plateDispControl[f] = r.num(5, 0.0) != 0;
        }
        return true;
    }
    if (key == "GEOPLAST") {   // C16: GEOPLAST model(0/1) match(0/1/2) dilation H c_res permCoef permMax [iter tol]
        G.plasticity = int(r.num(1, 1.0)); G.dpMatch = int(r.num(2, 2.0)); G.dilationDeg = r.num(3, 0.0);
        G.hardening = r.num(4, 0.0); G.residualCohesion = r.num(5, 0.0); G.permPlastic = r.num(6, 0.0);
        G.permPlasticMax = r.num(7, 100.0); G.plasticIter = int(r.num(8, 80.0)); G.plasticTol = r.num(9, 1e-6);
        return true;
    }
    return false;
}

std::string geoKeywordsText(const geomech::Input& G) {
    if (!G.enabled) return {};
    std::ostringstream o;
    auto f0 = [](const std::vector<double>& v, double d) { return v.empty() ? d : v[0]; };
    o << "-- RESSIM: GEOMECH ON " << G.coupling << ' ' << num(f0(G.E, 5e9)) << ' ' << num(f0(G.nu, 0.25)) << ' ' << num(f0(G.biot, 1.0)) << ' '
      << num(f0(G.alphaT, 1e-5)) << ' ' << num(G.Ks) << ' ' << G.overLayers << ' ' << G.underLayers << ' ' << G.sidePad << ' '
      << num(G.overThick) << ' ' << num(G.underThick) << ' ' << num(G.sideThick) << ' ' << num(G.padE) << ' ' << num(G.padNu) << ' '
      << num(G.padAlphaT) << ' ' << (G.topFree ? 1 : 0) << ' ' << num(G.sigmaVGrad) << ' ' << num(G.k0) << ' ' << num(G.cohesion) << ' '
      << num(G.frictionDeg) << ' ' << num(G.tensile) << ' ' << G.maxIter << ' ' << num(G.tol);
    if (G.permStress > 0) o << ' ' << num(G.permStress);
    o << "\n";
    for (int f = 0; f < 6; ++f)
        if (G.bcType[f] >= 0 || G.bcStress[f] != 0) {
            o << "-- RESSIM: MECHBC " << f << ' ' << G.bcType[f] << ' ' << num(G.bcStress[f]);
            if (G.bcType[f] == 3) o << ' ' << num(G.bcDisp[f]) << ' ' << (G.plateDispControl[f] ? 1 : 0);
            o << "\n";
        }
    if (G.plasticity)
        o << "-- RESSIM: GEOPLAST " << G.plasticity << ' ' << G.dpMatch << ' ' << num(G.dilationDeg) << ' ' << num(G.hardening) << ' '
          << num(G.residualCohesion) << ' ' << num(G.permPlastic) << ' ' << num(G.permPlasticMax) << ' ' << G.plasticIter << ' ' << num(G.plasticTol) << "\n";
    return o.str();
}

bool thermalFromDeck(const deck::RawDeck& D, Case& c, std::vector<DeckMessage>* msgs) {
    bool th = false;
    for (const deck::Keyword& k : D.keywords) if (k.name == "THERMAL") th = true;
    if (!th) return false;
    const bool field = c.units == UnitSystem::Field;
    const TU U(field);
    ThermalInput& T = c.thermal;
    // GRID / SOLUTION arrays were read by the importer (BOX / EQUALS / ... applied) - keep them
    const std::vector<double> heatcr = T.heatcr, heatcrt = T.heatcrt, thconr = T.thconr, thconsf = T.thconsf, tempi = T.tempi;
    T = ThermalInput();
    T.heatcr = heatcr; T.heatcrt = heatcrt; T.thconr = thconr; T.thconsf = thconsf; T.tempi = tempi;
    T.enabled = true;
    std::vector<std::string> notes;
    std::map<std::string, const deck::Keyword*> first;
    std::string section; int step = 0, nc = 0;
    static const char* kSections[] = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};
    std::vector<std::string> names;
    std::map<std::string, std::string> ctypes;   // '-- RESSIM: CTYPE name SOLID|GAS'
    std::map<std::string, std::vector<std::string>> aqsol, crit, formula;   // C10 per-component RESSIM lines (raw items)
    bool mechKw = false;
    for (const deck::Keyword& k : D.keywords) {
        bool isSec = false;
        for (const char* q : kSections) if (k.name == q) { section = q; isSec = true; }
        if (isSec) continue;
        if (section != "SCHEDULE" && !first.count(k.name)) first[k.name] = &k;
        if (k.name == "COMPS" && !k.records.empty() && !k.records[0].items.empty()) nc = k.records[0].integer(0, 0);
        else if (k.name == "DEADOIL") T.deadOil = true;
        else if (k.name == "CNAMES") { for (const deck::Record& r : k.records) for (size_t i = 0; i < r.items.size(); ++i) names.push_back(r.str(i)); }
        else if (k.name == "WINJTEMP" && section == "SCHEDULE") {
            for (const deck::Record& r : k.records) {
                if (r.terminator || r.items.empty()) continue;
                ThermalWinj w; w.well = r.str(0); w.step = step;
                w.quality = r.num(1, 1.0);
                w.temp = r.has(2) ? r.num(2, 0.0) : -1e9;
                w.pres = r.num(3, 0.0); w.enth = r.num(4, 0.0);
                T.winj.push_back(w);
            }
        } else if (k.name == "__RESSIM") {
            for (const deck::Record& r : k.records) {
                if (r.items.empty()) continue;
                const std::string key = up(r.str(0));
                if (key == "HEATLOSS") {
                    T.heatLoss = up(r.str(1, "ON")) != "OFF";
                    T.obCv = r.num(2, 0.0); T.obK = r.num(3, 0.0);
                } else if (key == "CTYPE") {
                    ctypes[up(r.str(1))] = up(r.str(2, "KVAL"));
                } else if (key == "REACTION") {
                    // REACTION name A Ea dH tMin o2pp reversible keqA keqB STO .. PROD .. [ORDER ..] [PHASE ..]
                    ThermalReaction rx; rx.name = r.str(1, "R");
                    rx.A = r.num(2, 0.0); rx.Ea = r.num(3, 0.0); rx.dH = r.num(4, 0.0); rx.tMin = r.num(5, 0.0);
                    rx.o2pp = int(r.num(6, -1.0)); rx.reversible = r.num(7, 0.0) != 0; rx.keqA = r.num(8, 0.0); rx.keqB = r.num(9, 0.0);
                    int mode = -1;
                    for (size_t i = 10; i < r.items.size(); ++i) {
                        const std::string t = up(r.str(i));
                        if (t == "STO") { mode = 0; continue; } if (t == "PROD") { mode = 1; continue; }
                        if (t == "ORDER") { mode = 2; continue; } if (t == "PHASE") { mode = 3; continue; }
                        const double v = r.num(i, 0.0);
                        if (mode == 0) rx.sto.push_back(v); else if (mode == 1) rx.prod.push_back(v);
                        else if (mode == 2) rx.order.push_back(v); else if (mode == 3) rx.phase.push_back(int(v));
                    }
                    T.reactions.push_back(rx);
                } else if (key == "AQSOL" || key == "CRIT" || key == "FORMULA") {
                    // AQSOL name LIB | AQSOL name model a b c vInf ksalt ; CRIT name tc pc omega ; FORMULA name formula [trackEq]
                    std::vector<std::string> v; for (size_t i = 2; i < r.items.size(); ++i) v.push_back(r.str(i));
                    (key == "AQSOL" ? aqsol : key == "CRIT" ? crit : formula)[up(r.str(1))] = v;
                } else if (key == "GASEOS") {
                    T.gasPR = up(r.str(1, "IDEAL")) == "PR";
                } else if (key == "SALINITY") {
                    T.salinity = r.num(1, 0.0);
                } else if (key == "HYSTGAS") {   // HYSTGAS sgtMax | HYSTGAS OFF
                    const std::string a = up(r.str(1, ""));
                    T.hystGas = a != "OFF"; if (T.hystGas) T.sgtMax = r.num(1, 0.3);
                } else if (key == "TRACK") {     // TRACK name [structTol]
                    T.track = up(r.str(1, "")); T.structTol = r.num(2, 0.01);
                } else if (key == "FRACTURE") {   // FRACTURE name x0 y0 x1 y1 zTop zBot aperture perm [MECH kn0 vm mu ks dil]
                    ThermalFracture f; f.name = r.str(1, "F");
                    f.x0 = r.num(2, 0.0); f.y0 = r.num(3, 0.0); f.x1 = r.num(4, 0.0); f.y1 = r.num(5, 0.0);
                    f.zTop = r.num(6, -1.0); f.zBot = r.num(7, -1.0); f.aperture = r.num(8, 1e-3); f.perm = r.num(9, 0.0);
                    if (up(r.str(10, "")) == "MECH") {
                        f.mech = true; f.kn0 = r.num(11, 20e9); f.vm = r.num(12, 0.4e-3); f.mu = r.num(13, 0.6); f.ks = r.num(14, 0.0); f.dilationDeg = r.num(15, 5.0);
                    }
                    T.fractures.push_back(f);
                } else if (key == "FRACSTRESS") {   // FRACSTRESS sHmaxGrad sHminGrad azimuth E nu alphaT
                    ThermalFracStress& F = T.fracStress;
                    F.sHmaxGrad = r.num(1, 0.0); F.sHminGrad = r.num(2, 0.0); F.azimuthDeg = r.num(3, 0.0);
                    F.E = r.num(4, 30e9); F.nu = r.num(5, 0.25); F.alphaT = r.num(6, 0.0);
                } else if (key == "SORB") {   // SORB comp capacity b tRef dH tau [AUTO|GAS|AQ] [adsRho] [site] [EQ|NOEQ]
                    ThermalSorb so; so.comp = up(r.str(1, ""));
                    so.capacity = r.num(2, 0.0); so.b = r.num(3, 0.0); so.tRef = r.num(4, 298.15); so.dH = r.num(5, 0.0); so.tau = r.num(6, 3600.0);
                    const std::string ph = up(r.str(7, "AUTO"));
                    so.phase = ph == "GAS" ? 2 : ph == "AQ" ? 0 : -1;
                    so.adsRho = r.num(8, 0.0);
                    if (r.has(9)) so.site = up(r.str(9, "SITES"));
                    so.initEq = up(r.str(10, "EQ")) != "NOEQ";
                    if (!so.comp.empty()) T.sorb.push_back(so);
                } else if (key == "STEAMTRAP") {   // STEAMTRAP well dT (SCHEDULE: from the current report step)
                    ThermalTrap tp; tp.well = r.str(1); tp.step = section == "SCHEDULE" ? step : 0; tp.subcool = r.num(2, 0.0); tp.steamMax = r.num(3, 0.0);
                    if (!tp.well.empty()) T.traps.push_back(tp);
                } else if (key == "SOLIDCONC") {
                    T.solidConc.clear(); for (size_t i = 1; i < r.items.size(); ++i) T.solidConc.push_back(r.num(i, 0.0));
                } else if (key == "PERMLAW") {
                    T.permLaw = int(r.num(1, 0.0)); T.permC = r.num(2, 20.0); T.permMultMax = r.num(3, 1e5);
                } else if (parseGeoRecord(key, r, T.geo)) {                }
            }
        } else if (k.name == "MECH") {
            T.geo.enabled = true; mechKw = true;
        } else if (k.name == "WELLSTRE") {
            for (const deck::Record& r : k.records) {
                if (r.terminator || r.items.empty()) continue;
                CompStream st; st.name = up(r.str(0));
                for (size_t i = 1; i < r.items.size(); ++i) st.z.push_back(r.num(i, 0.0));
                bool found = false;
                for (auto& x : T.streams) if (x.name == st.name) { x = st; found = true; }
                if (!found) T.streams.push_back(st);
            }
        } else if (k.name == "WINJGAS" && section == "SCHEDULE") {
            for (const deck::Record& r : k.records) {
                if (r.terminator || r.items.empty()) continue;
                const std::string how = up(r.str(1, "STREAM"));
                if ((how == "STREAM" || how == "WELLSTRE") && size_t(step) < c.schedule.size()) {
                    auto& v = c.schedule[size_t(step)].injStreams;
                    const std::pair<std::string, std::string> e{r.str(0), up(r.str(2))};
                    if (std::find(v.begin(), v.end(), e) == v.end()) v.push_back(e);
                }
            }
        } else if (k.name == "TSTEP" && section == "SCHEDULE") {
            int nv = 0; for (const deck::Record& r : k.records) nv += int(r.items.size());
            step += nv;
        } else if (k.name == "DATES" && section == "SCHEDULE") {
            for (const deck::Record& r : k.records) if (!r.terminator && !r.items.empty()) ++step;
        }
    }
    if (nc <= 0) nc = int(names.size());
    auto vals = [&](const char* name) {
        std::vector<double> v;
        auto it = first.find(name);
        if (it != first.end()) v = recordValues(*it->second);
        return v;
    };
    if (nc <= 0) { const std::vector<double> mw = vals("MW"); nc = int(mw.size()); }
    if (nc <= 0) {
        if (msgs) { DeckMessage m; m.level = 1; m.keyword = "THERMAL"; m.text = "no hydrocarbon components (COMPS / CNAMES) - the thermal model is not built"; msgs->push_back(m); }
        T = ThermalInput();
        return false;
    }
    auto at = [](const std::vector<double>& v, size_t i, double d) { return i < v.size() && std::isfinite(v[i]) ? v[i] : d; };
    const std::vector<double> mw = vals("MW"), tc = vals("TCRIT"), pc = vals("PCRIT"), acf = vals("ACF"), kv = vals("KVCR"),
        dref = vals("DREF"), cref = vals("CREF"), tx1 = vals("THERMEX1"), sa = vals("SPECHA"), sb = vals("SPECHB"), sg = vals("SPECHG"),
        sh = vals("SPECHH"), hvs = vals("HEATVAPS"), hv = vals("HEATVAP"), hve = vals("HEATVAPE"), zf = vals("ZFACTOR"), zf1 = vals("ZFACT1"),
        gvf = vals("GASVISCF");
    std::vector<std::string> cvt;
    if (first.count("CVTYPE")) for (const deck::Record& r : first["CVTYPE"]->records) { if (r.terminator) continue; for (size_t i = 0; i < r.items.size(); ++i) cvt.push_back(up(r.str(i))); break; }
    const bool kvwi = first.count("KVWI") > 0;
    const size_t n = size_t(nc);
    bool wilson = false, noK = false;
    for (size_t i = 0; i < n; ++i) {
        ThermalComponent q;
        q.name = i < names.size() ? names[i] : "C" + std::to_string(i + 1);
        q.mw = at(mw, i, 100.0);
        if (!(at(mw, i, 0.0) > 0)) notes.push_back("MW missing for component " + q.name + " - 100 used");
        q.tcrit = at(tc, i, 1e5); q.pcrit = at(pc, i, 1.0);
        const std::string ty = i < cvt.size() ? cvt[i] : "";
        bool anyK = false;
        for (int j = 0; j < 5; ++j) { q.kvcr[j] = at(kv, size_t(j) * n + i, 0.0); anyK = anyK || q.kvcr[j] != 0; }
        q.dead = ty == "DEAD" || (ty.empty() && q.tcrit >= 1e4 && !anyK);
        if (!q.dead && !anyK) {
            if (q.tcrit <= 1e-6 || ty == "GAS") { q.kvcr[0] = 1e4; }   // non-condensable: stays in the gas
            else if (at(pc, i, 0.0) > 0 && at(tc, i, 0.0) > 0) {        // Wilson K-values (E300 KVWI form)
                const double a = 5.372697 * (1 + at(acf, i, 0.0));
                q.kvcr[1] = q.pcrit * std::exp(a); q.kvcr[3] = a * q.tcrit;
                wilson = true;
            } else { q.dead = true; noK = true; }
        }
        q.dref = at(dref, i, 950.0 / U.dens);
        q.cref = at(cref, i, 0.0); q.thermex1 = at(tx1, i, 0.0);
        q.specha = at(sa, i, 0.0); q.spechb = at(sb, i, 0.0); q.spechg = at(sg, i, 0.0); q.spechh = at(sh, i, 0.0);
        q.heatvaps = at(hvs, i, 0.0); q.heatvap = at(hv, i, 0.0); q.heatvape = at(hve, i, 0.38);
        q.zfactor = at(zf, i, 0.96); q.zfact1 = at(zf1, i, 0.0);
        q.gvA = at(gvf, i, 0.0); q.gvB = at(gvf, n + i, 0.0);
        if (ty == "GAS") q.type = 1;
        auto ct = ctypes.find(up(q.name));
        if (ct != ctypes.end()) { if (ct->second == "SOLID") { q.type = 2; q.dead = true; } else if (ct->second == "GAS") q.type = 1; else q.type = 0; }
        auto dnum = [](const std::vector<std::string>& v, size_t i, double d) {
            if (i >= v.size()) return d;
            char* e = nullptr; const double x = std::strtod(v[i].c_str(), &e);
            return (e && e != v[i].c_str()) ? x : d;
        };
        auto aq = aqsol.find(up(q.name));
        if (aq != aqsol.end() && !aq->second.empty()) {
            const std::vector<std::string>& v = aq->second;
            if (up(v[0]) == "LIB") {
                const gasdata::GasInfo* gi = gasdata::find(q.name);
                if (gi && gi->hA != 0) { q.aqModel = 1; q.hA = gi->hA; q.hB = gi->hB; q.hC = gi->hC; q.vInf = gi->vInf; q.ksalt = gi->ksalt; }
                else notes.push_back("AQSOL " + q.name + " LIB: no library Henry data for this name (component stays insoluble)");
            } else {
                q.aqModel = int(dnum(v, 0, 0.0));
                if (q.aqModel == 2) { q.kh298 = dnum(v, 1, 0.0); q.hdT = dnum(v, 2, 0.0); }
                else { q.hA = dnum(v, 1, 0.0); q.hB = dnum(v, 2, 0.0); q.hC = dnum(v, 3, 0.0); }
                q.vInf = dnum(v, 4, 0.0); q.ksalt = dnum(v, 5, 0.0);
            }
        }
        auto cr = crit.find(up(q.name));
        if (cr != crit.end()) { q.tcG = dnum(cr->second, 0, 0.0); q.pcG = dnum(cr->second, 1, 0.0); q.omegaG = dnum(cr->second, 2, 0.0); }
        auto fo = formula.find(up(q.name));
        if (fo != formula.end() && !fo->second.empty()) { q.formula = fo->second[0]; q.trackEq = dnum(fo->second, 1, 0.0); }
        T.comps.push_back(q);
    }
    if (wilson) notes.push_back(std::string("volatile component(s) without KVCR: Wilson K-values from PCRIT / TCRIT / ACF") + (kvwi ? " (KVWI)" : " (E300 default)"));
    if (noK) notes.push_back("component(s) without CVTYPE / KVCR and without PCRIT / TCRIT: treated as non-volatile (dead)");
    for (const auto& q : T.comps) if (q.zfact1 != 0) { notes.push_back("ZFACT1 is kept for ECLIPSE but not used by the ResSim thermal core"); break; }
    auto rows = [&](const char* name, size_t width) {
        std::vector<std::vector<double>> out;
        const std::vector<double> v = vals(name);
        for (size_t r = 0; r + width <= v.size(); r += width) out.push_back(std::vector<double>(v.begin() + long(r), v.begin() + long(r + width)));
        return out;
    };
    T.oilViscT = rows("OILVISCT", n + 1);
    T.gasViscT = rows("GASVISCT", n + 1);
    for (auto& r : rows("WATVISCT", 2)) T.watViscT.push_back({r[0], r[1]});
    if (first.count("STCOND")) { const std::vector<double> v = vals("STCOND"); T.stdT = at(v, 0, U.stdT()); T.stdP = at(v, 1, U.stdP()); }
    for (auto& r : rows("TEMPVD", 2)) T.tempvd.push_back({r[0], r[1]});
    T.zi = vals("ZI"); if (!T.zi.empty()) { T.zi.resize(n, 0.0); for (double& x : T.zi) if (!std::isfinite(x)) x = 0; }
    const size_t ncell = size_t(c.grid.nx) * size_t(c.grid.ny) * size_t(c.grid.nz);
    auto cellComp = [&](const char* name, std::vector<double>& dst) {
        auto it = first.find(name);
        if (it == first.end()) return;
        dst = it->second->values;
        if (dst.size() != ncell * n) { notes.push_back(std::string(name) + ": expected NX*NY*NZ*NCOMPS values - ignored"); dst.clear(); }
        for (double& x : dst) if (!std::isfinite(x)) x = 0;
    };
    cellComp("XMF", T.xmf); cellComp("YMF", T.ymf);
    // C8 data sized to the component count
    for (auto& rx : T.reactions) {
        rx.sto.resize(n + 1, 0.0); rx.prod.resize(n + 1, 0.0);
        if (!rx.order.empty()) rx.order.resize(n + 1, 0.0);
        if (!rx.phase.empty()) rx.phase.resize(n + 1, -1);
    }
    if (!T.solidConc.empty()) T.solidConc.resize(n, 0.0);
    for (auto& st : T.streams) st.z.resize(n, 0.0);
    // OPM mechanics keywords (MECH + YMODULE GPa / PRATIO / BIOTCOEF / THERMEXR 1/K, GRID arrays)
    auto cellArr = [&](const char* name, double f, std::vector<double>& dst) {
        auto it = first.find(name);
        if (it == first.end() || it->second->values.empty()) return;
        dst = it->second->values;
        for (double& x : dst) x = std::isfinite(x) ? x * f : 0.0;
        bool uni = true; for (double x : dst) if (x != dst[0]) uni = false;
        if (uni) dst.resize(1);
    };
    if (mechKw || first.count("YMODULE")) {
        cellArr("YMODULE", 1e9, T.geo.E); cellArr("PRATIO", 1.0, T.geo.nu); cellArr("BIOTCOEF", 1.0, T.geo.biot);
        cellArr("THERMEXR", field ? 1.8 : 1.0, T.geo.alphaT);
        if (mechKw) notes.push_back("MECH: geomechanics (linear thermo-poroelastic FEM, fixed-stress coupling) from YMODULE / PRATIO / BIOTCOEF / THERMEXR");
    }
    if (!T.reactions.empty()) notes.push_back(std::to_string(T.reactions.size()) + " reaction(s) ('-- RESSIM: REACTION'): run by the ResSim thermal core only");
    if (first.count("PREF") || first.count("TREF")) notes.push_back("PREF / TREF: the ResSim thermal core uses the standard conditions (STCOND) as reference");
    if (first.count("OILVISCC") || first.count("OILVINDX")) notes.push_back("OILVISCC / OILVINDX correlations are not used - give OILVISCT");
    if (first.count("THSTT97") == 0) notes.push_back("water / steam: the ResSim thermal core always uses IAPWS-IF97 (ECLIPSE: THSTT97)");
    if (T.oilViscT.empty()) notes.push_back("no OILVISCT: oil-phase viscosity 1 cP");
    if (T.tempi.empty() && T.tempvd.empty()) notes.push_back("no TEMPI / TEMPVD: initial temperature = standard temperature");
    T.structured = false;
    // the case is now represented: remove the "black-oil engine cannot model it" markers
    auto& Un = c.engineUnsupported;
    Un.erase(std::remove_if(Un.begin(), Un.end(), [&](const std::string& r) {
        return r.rfind("thermal model", 0) == 0 || r.rfind("compositional model", 0) == 0;
    }), Un.end());
    c.comp = CompositionalInput();
    if (msgs) {
        DeckMessage m; m.level = 0; m.keyword = "THERMAL";
        std::string nm; for (auto& q : T.comps) nm += (nm.empty() ? "" : " ") + q.name + (q.dead ? "(dead)" : "");
        m.text = std::to_string(nc) + "-component thermal model (" + nm + "), " + std::to_string(T.winj.size()) +
                 " WINJTEMP record(s) - runs on the ResSim thermal core (steam tables IAPWS-IF97)";
        msgs->push_back(m);
        for (const auto& s : notes) { DeckMessage w; w.level = 1; w.keyword = "THERMAL"; w.text = s; msgs->push_back(w); }
    }
    return true;
}

// ======================================================================== Case::thermal -> deck text
std::string thermalKeywordsText(const Case& c, int section, int step) {
    const ThermalInput& T = c.thermal;
    if (!T.enabled) return "";
    const bool field = c.units == UnitSystem::Field;
    const TU U(field);
    const size_t n = T.comps.size();
    std::ostringstream o;
    auto perComp = [&](const char* name, auto get, const char* comment = nullptr) {
        o << name << "\n";
        if (comment) o << "-- " << comment << "\n";
        for (size_t i = 0; i < n; ++i) o << ' ' << num(get(T.comps[i]));
        o << " /\n\n";
    };
    auto cellArray = [&](const char* name, const std::vector<double>& v) {
        if (v.empty()) return;
        o << name << "\n";
        const size_t ncell = size_t(c.grid.nx) * size_t(c.grid.ny) * size_t(c.grid.nz);
        if (v.size() == 1) o << ' ' << ncell << '*' << num(v[0]) << " /\n\n";
        else {
            size_t i = 0, col = 0;
            while (i < v.size()) {
                size_t j = i; while (j < v.size() && v[j] == v[i]) ++j;
                o << ' ' << (j - i > 1 ? std::to_string(j - i) + "*" : std::string()) << num(v[i]);
                if (++col % 8 == 0) o << "\n";
                i = j;
            }
            o << " /\n\n";
        }
    };
    // RUNSPEC / PROPS / SCHEDULE keywords and the initial compositions: written from the data only in structured mode
    // (imported decks keep them verbatim); the rock heat arrays and TEMPI are always written from the data (the importer
    // reads them as arrays with BOX / EQUALS / ... applied)
    if (!T.structured && (section == 0 || section == 1 || section == 3)) return "";
    if (section == 0) {
        o << "THERMAL\n";
        if (T.deadOil) o << "DEADOIL\n";
        o << "COMPS\n " << n << " /\n";
        return o.str();
    }
    if (section == 4) {
        cellArray("HEATCR", T.heatcr); cellArray("HEATCRT", T.heatcrt); cellArray("THCONR", T.thconr); cellArray("THCONSF", T.thconsf);
        return o.str();
    }
    if (section == 1) {
        o << "-- ---------------- thermal fluid (E300 THERMAL, K-value components; water / steam IAPWS-IF97)\n";
        o << "CNAMES\n"; for (auto& q : T.comps) o << " '" << q.name << "'"; o << " /\n\n";
        perComp("MW", [](const ThermalComponent& q) { return q.mw; });
        perComp("TCRIT", [](const ThermalComponent& q) { return q.tcrit; });
        perComp("PCRIT", [](const ThermalComponent& q) { return q.pcrit; });
        o << "CVTYPE\n"; for (auto& q : T.comps) o << ' ' << (q.type == 1 ? "GAS" : (q.dead || q.type == 2) ? "DEAD" : "LIVE"); o << " /\n\n";
        bool anyK = false; for (auto& q : T.comps) if (!q.dead && q.type != 2) anyK = true;
        if (anyK) {
            o << "KVCR\n-- A (one value per component), B, C, D, E: K = (A + B/p + C p) exp(-D / (T - E))\n";
            for (int j = 0; j < 5; ++j) {
                for (auto& q : T.comps) o << ' ' << num((q.dead || q.type == 2) ? 0.0 : (q.type == 1 && j == 0 && q.kvcr[0] == 0 ? 1e4 : q.kvcr[j]));
                o << "\n";
            }
            o << "/\n\n";
        }
        perComp("DREF", [](const ThermalComponent& q) { return q.dref; });
        perComp("CREF", [](const ThermalComponent& q) { return q.cref; });
        perComp("THERMEX1", [](const ThermalComponent& q) { return q.thermex1; });
        perComp("SPECHA", [](const ThermalComponent& q) { return q.specha; });
        perComp("SPECHB", [](const ThermalComponent& q) { return q.spechb; });
        perComp("SPECHG", [](const ThermalComponent& q) { return q.spechg; });
        perComp("SPECHH", [](const ThermalComponent& q) { return q.spechh; });
        perComp("HEATVAPS", [](const ThermalComponent& q) { return q.heatvaps; });
        perComp("HEATVAP", [](const ThermalComponent& q) { return q.heatvap; });
        perComp("HEATVAPE", [](const ThermalComponent& q) { return q.heatvape; });
        perComp("ZFACTOR", [](const ThermalComponent& q) { return q.zfactor; });
        bool z1 = false; for (auto& q : T.comps) if (q.zfact1 != 0) z1 = true;
        if (z1) perComp("ZFACT1", [](const ThermalComponent& q) { return q.zfact1; });
        auto table = [&](const char* name, const std::vector<std::vector<double>>& t) {
            if (t.empty()) return;
            o << name << "\n";
            const size_t nw = size_t(T.nc()) + 1;   // exactly T + one column per component (the importer reads nc + 1 values per row)
            for (auto& r : t) {
                if (r.empty()) continue;
                for (size_t i = 0; i < nw; ++i) o << ' ' << num(i < r.size() ? r[i] : r.back());
                o << "\n";
            }
            o << "/\n\n";
        };
        table("OILVISCT", T.oilViscT);
        bool gvf = false; for (auto& q : T.comps) if (q.gvA != 0) gvf = true;
        if (gvf) {
            o << "GASVISCF\n";
            for (auto& q : T.comps) { o << ' ' << num(q.gvA); }
            o << "\n";
            for (auto& q : T.comps) { o << ' ' << num(q.gvB); }
            o << "\n/\n\n";
        } else table("GASVISCT", T.gasViscT);
        o << "WATVISCT\n";
        if (!T.watViscT.empty()) for (auto& r : T.watViscT) o << ' ' << num(r.first) << ' ' << num(r.second) << "\n";
        else {   // IAPWS-2008 liquid viscosity on the saturation line (compressed liquid at 1.05 Psat, >= 1 bar)
            o << "-- IAPWS-2008 (ResSim water model)\n";
            for (int i = 0; i <= 17; ++i) {
                const double TK = 283.15 + 20.0 * i;
                const double p = std::max(1e5, 1.05 * steam::psat(TK));
                const double rho = 1.0 / steam::liquid(TK, p).v;
                char b[80]; std::snprintf(b, sizeof b, " %.6g %.6g\n", U.tempDeck(TK), steam::viscosity(TK, rho) * 1e3);
                o << b;
            }
        }
        o << "/\n\n";
        o << "THSTT97\n\n";
        o << "STCOND\n " << num(T.stdT > -1e8 ? T.stdT : U.stdT()) << ' ' << num(T.stdP > 0 ? T.stdP : U.stdP()) << " /\n\n";
        if (T.xmf.empty() && !T.zi.empty()) { o << "ZI\n"; for (double x : T.zi) o << ' ' << num(x); o << " /\n\n"; }
        if (T.heatLoss) o << "-- RESSIM: HEATLOSS ON " << num(T.obCv) << ' ' << num(T.obK) << "\n\n";
        // C8: chemistry / solids / geomechanics (ResSim thermal core; ECLIPSE ignores the comment lines)
        bool c8 = false;
        for (auto& q : T.comps) if (q.type == 2) { o << "-- RESSIM: CTYPE " << q.name << " SOLID\n"; c8 = true; }
        for (const ThermalReaction& rx : T.reactions) {
            o << "-- RESSIM: REACTION " << rx.name << ' ' << num(rx.A) << ' ' << num(rx.Ea) << ' ' << num(rx.dH) << ' ' << num(rx.tMin) << ' '
              << rx.o2pp << ' ' << (rx.reversible ? 1 : 0) << ' ' << num(rx.keqA) << ' ' << num(rx.keqB) << " STO";
            for (double v : rx.sto) o << ' ' << num(v);
            o << " PROD"; for (double v : rx.prod) o << ' ' << num(v);
            if (!rx.order.empty()) { o << " ORDER"; for (double v : rx.order) o << ' ' << num(v); }
            if (!rx.phase.empty()) { o << " PHASE"; for (int v : rx.phase) o << ' ' << v; }
            o << "\n"; c8 = true;
        }
        if (!T.solidConc.empty()) { o << "-- RESSIM: SOLIDCONC"; for (double v : T.solidConc) o << ' ' << num(v); o << "\n"; c8 = true; }
        if (T.permLaw) { o << "-- RESSIM: PERMLAW " << T.permLaw << ' ' << num(T.permC) << ' ' << num(T.permMultMax) << "\n"; c8 = true; }
        if (T.geo.enabled) {
            o << geoKeywordsText(T.geo);
            c8 = true;
        }
        // C10: dissolved gases / PR gas / hysteresis / tracking
        for (auto& q : T.comps) {
            if (q.aqModel == 1) { o << "-- RESSIM: AQSOL " << q.name << " 1 " << num(q.hA) << ' ' << num(q.hB) << ' ' << num(q.hC) << ' ' << num(q.vInf) << ' ' << num(q.ksalt) << "\n"; c8 = true; }
            else if (q.aqModel == 2) { o << "-- RESSIM: AQSOL " << q.name << " 2 " << num(q.kh298) << ' ' << num(q.hdT) << " 0 " << num(q.vInf) << ' ' << num(q.ksalt) << "\n"; c8 = true; }
            if (q.tcG > 0) { o << "-- RESSIM: CRIT " << q.name << ' ' << num(q.tcG) << ' ' << num(q.pcG) << ' ' << num(q.omegaG) << "\n"; c8 = true; }
            if (!q.formula.empty()) { o << "-- RESSIM: FORMULA " << q.name << ' ' << q.formula; if (q.trackEq != 0) o << ' ' << num(q.trackEq); o << "\n"; c8 = true; }
        }
        if (T.gasPR) { o << "-- RESSIM: GASEOS PR\n"; c8 = true; }
        if (T.salinity > 0) { o << "-- RESSIM: SALINITY " << num(T.salinity) << "\n"; c8 = true; }
        if (T.hystGas) { o << "-- RESSIM: HYSTGAS " << num(T.sgtMax) << "\n"; c8 = true; }
        if (!T.track.empty()) { o << "-- RESSIM: TRACK " << T.track << ' ' << num(T.structTol) << "\n"; c8 = true; }
        for (const ThermalFracture& f : T.fractures) {
            o << "-- RESSIM: FRACTURE " << f.name << ' ' << num(f.x0) << ' ' << num(f.y0) << ' ' << num(f.x1) << ' ' << num(f.y1) << ' ' << num(f.zTop) << ' '
              << num(f.zBot) << ' ' << num(f.aperture) << ' ' << num(f.perm);
            if (f.mech) o << " MECH " << num(f.kn0) << ' ' << num(f.vm) << ' ' << num(f.mu) << ' ' << num(f.ks) << ' ' << num(f.dilationDeg);
            o << "\n"; c8 = true;
        }
        if (!T.fractures.empty()) {
            const ThermalFracStress& F = T.fracStress;
            o << "-- RESSIM: FRACSTRESS " << num(F.sHmaxGrad) << ' ' << num(F.sHminGrad) << ' ' << num(F.azimuthDeg) << ' ' << num(F.E) << ' ' << num(F.nu) << ' ' << num(F.alphaT) << "\n";
        }
        for (const ThermalSorb& so : T.sorb) {
            o << "-- RESSIM: SORB " << so.comp << ' ' << num(so.capacity) << ' ' << num(so.b) << ' ' << num(so.tRef) << ' ' << num(so.dH) << ' ' << num(so.tau) << ' '
              << (so.phase == 2 ? "GAS" : so.phase == 0 ? "AQ" : "AUTO") << ' ' << num(so.adsRho) << ' ' << so.site << ' ' << (so.initEq ? "EQ" : "NOEQ") << "\n";
            c8 = true;
        }
        if (c8) o << "\n";
        return o.str();
    }
    if (section == 2) {
        if (!T.tempi.empty()) cellArray("TEMPI", T.tempi);
        if (!T.structured) return o.str();
        if (!T.tempi.empty()) {}
        else if (!T.tempvd.empty()) {
            o << "TEMPVD\n";
            for (auto& r : T.tempvd) o << ' ' << num(r.first) << ' ' << num(r.second) << "\n";
            o << "/\n\n";
        }
        if (!T.xmf.empty()) { o << "XMF\n"; size_t col = 0; for (double x : T.xmf) { o << ' ' << num(x); if (++col % 8 == 0) o << "\n"; } o << " /\n\n"; }
        if (!T.ymf.empty()) { o << "YMF\n"; size_t col = 0; for (double x : T.ymf) { o << ' ' << num(x); if (++col % 8 == 0) o << "\n"; } o << " /\n\n"; }
        return o.str();
    }
    if (section == 3) {
        if (step == 0 && !T.streams.empty()) {
            o << "WELLSTRE\n";
            for (const auto& st : T.streams) { o << " '" << st.name << "'"; for (double z : st.z) o << ' ' << num(z); o << " /\n"; }
            o << "/\n\n";
        }
        if (size_t(step) < c.schedule.size() && !c.schedule[size_t(step)].injStreams.empty()) {
            o << "WINJGAS\n";
            for (auto& ws : c.schedule[size_t(step)].injStreams) o << " '" << ws.first << "' STREAM '" << ws.second << "' /\n";
            o << "/\n\n";
        }
        bool any = false;
        for (const ThermalWinj& w : T.winj) {
            if (w.step != step) continue;
            if (!any) o << "WINJTEMP\n";
            any = true;
            o << " '" << w.well << "' " << num(w.quality) << ' ' << (w.temp > -1e8 ? num(w.temp) : std::string("1*")) << ' '
              << (w.pres > 0 ? num(w.pres) : std::string("1*")) << ' ' << (w.enth > 0 ? num(w.enth) : std::string("1*")) << " /\n";
        }
        if (any) o << "/\n\n";
        for (const ThermalTrap& tp : T.traps)
            if (tp.step == step) o << "-- RESSIM: STEAMTRAP '" << tp.well << "' " << num(tp.subcool) << (tp.steamMax > 0 ? " " + num(tp.steamMax) : std::string()) << "\n";
        return o.str();
    }
    return "";
}

namespace thermal {

// ======================================================================== Case -> thermal core model
bool deckModel(const ressim::Case& c, DeckModel& out, std::string* error) {
    auto fail = [&](const std::string& m) { if (error) *error = m; return false; };
    {   // C17 review fix: grid arrays cover LGR cells (own local values), MINC shells and DPNUM = 0 folding, as in Simulator::setup
        const GridInput& g0 = c.grid;
        if ((!g0.lgrs.empty() && !g0.lgrExpanded) || (g0.minc > 1 && g0.dual == 1 && !g0.mincExpanded) || (g0.dual && !g0.dpnum.empty() && !g0.dpnumApplied)) {
            ressim::Case cx = c;
            try { expandLgrArrays(cx); expandMincArrays(cx); applyDpnum(cx); }
            catch (const std::exception& e) { return fail(std::string("grid: ") + e.what()); }
            return deckModel(cx, out, error);
        }
    }
    out = DeckModel();
    const bool field = c.units == UnitSystem::Field;
    out.field = field;
    const TU U(field);
    const ThermalInput& TI = c.thermal;
    if (!TI.enabled || TI.comps.empty()) return fail("the case has no thermal fluid (THERMAL with COMPS / CNAMES)");
    Case& cc = out.base;
    const size_t nc = TI.comps.size();
    cc.stdT = U.tempK(TI.stdT > -1e8 ? TI.stdT : U.stdT());
    cc.stdP = (TI.stdP > 0 ? TI.stdP : U.stdP()) * U.p;
    cc.tRef = cc.stdT; cc.pRef = cc.stdP;
    // ---- components (deck units -> SI, molar)
    for (size_t i = 0; i < nc; ++i) {
        const ThermalComponent& q = TI.comps[i];
        Component k;
        k.name = q.name;
        k.mw = q.mw * 1e-3;
        k.dead = q.dead || q.type == 2;
        k.type = q.type;
        if (q.type == 1) k.dead = false;
        if (!k.dead && q.type == 0) {
            k.kv1 = q.kvcr[0]; k.kv2 = q.kvcr[1] * U.p; k.kv3 = q.kvcr[2] / U.p;
            k.kv4 = -q.kvcr[3] * U.absT; k.kv5 = q.kvcr[4] * U.absT;
        }
        k.rhoRef = q.dref * U.dens;
        k.cpL = q.cref * U.comp;
        k.ct1 = q.thermex1 / U.absT;
        k.cpl1 = q.specha * U.spech * k.mw; k.cpl2 = q.spechb * U.spech / U.absT * k.mw;
        k.cpg1 = q.spechg * U.spech * k.mw; k.cpg2 = q.spechh * U.spech / U.absT * k.mw;
        k.hvaps = q.heatvaps * U.hv * k.mw; k.hvapA = q.heatvap * U.hv * k.mw; k.hvapB = q.heatvape;
        k.tcrit = q.tcrit * U.absT;
        k.z0 = q.zfactor > 0 ? q.zfactor : 0.96;
        // oil-phase viscosity: OILVISCT column (log-linear in T); default 1 cP
        if (!TI.oilViscT.empty()) {
            for (const auto& r : TI.oilViscT) if (r.size() > i + 1 && r[i + 1] > 0) { k.viscT.push_back(U.tempK(r[0])); k.viscV.push_back(r[i + 1] * 1e-3); }
            if (k.viscT.size() == 1) { k.avisc = k.viscV[0]; k.bvisc = 0; k.viscT.clear(); k.viscV.clear(); }
        } else { k.avisc = 1e-3; k.bvisc = 0; }
        // gas-phase viscosity: GASVISCF mu = A T^B (cP, T in K | R) or a power-law fit of GASVISCT
        if (q.gvA > 0) { k.avg = q.gvA * 1e-3 * std::pow(field ? 1.8 : 1.0, q.gvB); k.bvg = q.gvB; }
        else if (!TI.gasViscT.empty()) {
            std::vector<double> Tk, mu;
            for (const auto& r : TI.gasViscT) if (r.size() > i + 1) { Tk.push_back(U.tempK(r[0])); mu.push_back(r[i + 1] * 1e-3); }
            double A, B;
            if (fitPower(Tk, mu, A, B)) { k.avg = A; k.bvg = B; }
        }
        // C10: solubility (SI already), PR critical constants (library by name when GASEOS PR and no CRIT line), formula
        k.aqModel = q.aqModel; k.hA = q.hA; k.hB = q.hB; k.hC = q.hC; k.kh298 = q.kh298; k.hdT = q.hdT; k.vInf = q.vInf; k.ksalt = q.ksalt;
        k.tcG = q.tcG; k.pcG = q.pcG; k.omegaG = q.omegaG;
        const gasdata::GasInfo* gi = gasdata::find(q.name);
        if (TI.gasPR && !(k.tcG > 0) && gi && q.type != 2) { k.tcG = gi->tc; k.pcG = gi->pc; k.omegaG = gi->omega; }
        k.formula = !q.formula.empty() ? q.formula : (gi ? std::string(gi->formula) : std::string());
        k.trackEq = q.trackEq;
        cc.comps.push_back(k);
    }
    if (!TI.gasViscT.empty()) out.notes.push_back("GASVISCT: fitted by mu = A T^B per component");
    for (const auto& r : TI.watViscT) { cc.watViscT.push_back(U.tempK(r.first)); cc.watViscV.push_back(r.second * 1e-3); }
    // ---- C8: reactions, solids, permeability alteration (SI in the case data)
    for (const ThermalReaction& rx : TI.reactions) {
        Reaction r; r.name = rx.name; r.sto = rx.sto; r.prod = rx.prod; r.order = rx.order; r.phase = rx.phase;
        r.sto.resize(nc + 1, 0.0); r.prod.resize(nc + 1, 0.0);
        if (!r.order.empty()) r.order.resize(nc + 1, 0.0);
        if (!r.phase.empty()) r.phase.resize(nc + 1, -1);
        r.A = rx.A; r.Ea = rx.Ea; r.dH = rx.dH; r.tMin = rx.tMin; r.o2pp = rx.o2pp; r.reversible = rx.reversible; r.keqA = rx.keqA; r.keqB = rx.keqB;
        cc.reactions.push_back(r);
    }
    if (!TI.solidConc.empty()) { std::vector<double> s0 = TI.solidConc; s0.resize(nc, 0.0); cc.solid0.assign(1, s0); }
    cc.permLaw = TI.permLaw; cc.permC = TI.permC; cc.permMultMax = TI.permMultMax;
    cc.salinity = TI.salinity; cc.gasPR = TI.gasPR; cc.hyst.on = TI.hystGas; cc.hyst.sgtMax = TI.sgtMax; cc.structTol = TI.structTol;
    cc.track = -1;
    for (size_t i = 0; i < nc; ++i) {
        const std::string nm = up(TI.comps[i].name);
        const gasdata::GasInfo* gi = gasdata::find(nm);
        if (!TI.track.empty() ? nm == up(TI.track) : (gi && std::string(gi->name) == "CO2")) { cc.track = int(i); break; }
    }
    if (!TI.track.empty() && cc.track < 0) out.notes.push_back("TRACK " + TI.track + ": no such component (no storage inventory)");
    // ---- grid (processed connection list, SI)
    Grid g;
    if (!c.grid.hfracs.empty())
        return fail("HFRAC hydraulic fractures are for black-oil / compositional cases - thermal cases use '-- RESSIM: FRACTURE' "
                    "(EDFM with aperture mechanics, stimulation and propagation, SI units)");
    if (c.grid.langVL > 0) out.notes.push_back("LANGMUIR is ignored by the thermal engine - use SORB (temperature-dependent Langmuir sorption)");
    if (c.grid.permGamma > 0) out.notes.push_back("PERMSTRESS is ignored by the thermal engine - use the fracture mechanics (FRACTURE ... MECH) / ROCKCOMP options");
    try { g.build(c); } catch (const std::exception& e) { return fail(std::string("grid: ") + e.what()); }
    const int na = g.numActive();
    const GridInput& gr = c.grid;
    // C17 MINC: inner matrix shells (global index past the grid) take the per-cell input of their host matrix cell
    // (C17 review fix: LGR cells - also past the grid - read HEATCR / THCONR / TEMPI / XMF of their host cell instead of 0 / out of range)
    auto hostQ = [&](int gq) {
        if (g.isMincShell(gq)) return g.mincHost[size_t(gq - g.numCells() - g.numLgrCells())];
        if (g.isLocal(gq)) return g.lgrHost[size_t(gq - g.numCells())];
        return gq;
    };
    if (g.numMincCells() > 0)
        out.notes.push_back("MINC: " + std::to_string(g.mincN) + " nested matrix shells per block (" + std::to_string(g.numMincCells()) +
                            " shell cells) - heat conducts between the shells over the same proximity geometry as the flow");
    const double l3 = U.len * U.len * U.len;
    cc.cellPV.resize(size_t(na)); cc.cellDepth.resize(size_t(na)); cc.cellBulk.resize(size_t(na));
    for (int a = 0; a < na; ++a) {
        const int gq = g.globalOf[size_t(a)];
        cc.cellPV[size_t(a)] = g.poreVolRef[size_t(a)] * U.pv;
        cc.cellDepth[size_t(a)] = g.depth[size_t(a)] * U.len;
        double bulk = gq < int(g.bulkVol.size()) ? g.bulkVol[size_t(gq)] * l3 : 0.0;
        if (!(bulk > cc.cellPV[size_t(a)])) {
            const double phi = gq < int(gr.poro.size()) ? gr.poro[size_t(gq)] : 0.3;   // expanded grid array
            bulk = cc.cellPV[size_t(a)] / std::max(1e-3, phi);
        }
        cc.cellBulk[size_t(a)] = bulk;
    }
    const double mD = 9.869233e-16;
    auto permDir = [&](int gq, char d) {
        const std::vector<double>& k = (d == 'Y' && !gr.permy.empty()) ? gr.permy : (d == 'Z' && !gr.permz.empty()) ? gr.permz : gr.permx;
        return gq < int(k.size()) ? g.effPerm(k, gq) : 0.0;   // expanded grid array (LGR cells: their own PERM)
    };
    for (const Connection& cn : g.conns) {
        const double Tr = cn.trans * U.trans;
        const int ga = g.globalOf[size_t(cn.a)], gb = g.globalOf[size_t(cn.b)];
        const char d = cn.dir == 'M' ? 'X' : cn.dir;
        const double ka = permDir(ga, d), kb = permDir(gb, d);
        double G = 0;   // conduction geometry A / d = T / k_harmonic (multipliers act on conduction too)
        if (cn.geo > 0) G = cn.geo * U.len;   // C17: geometric A / d (matrix-fracture: sigma V, MINC shells: shell area / distance)
        else if (ka > 0 && kb > 0) G = Tr / (2 * ka * kb / (ka + kb) * mD);
        cc.conns.push_back({cn.a, cn.b, Tr, G});
    }
    out.cellGlobal = g.globalOf;
    // geomechanics: the mechanics mesh is the (tensor-product) Cartesian box of the grid
    if (TI.geo.enabled) {
        cc.geo = TI.geo;
        cc.geoNx = gr.nx; cc.geoNy = gr.ny; cc.geoNz = gr.nz;
        cc.geoDx.resize(size_t(gr.size())); cc.geoDy.resize(size_t(gr.size())); cc.geoDz.resize(size_t(gr.size()));
        for (int q = 0; q < gr.size(); ++q) {
            cc.geoDx[size_t(q)] = g.cellDx[size_t(q)] * U.len; cc.geoDy[size_t(q)] = g.cellDy[size_t(q)] * U.len; cc.geoDz[size_t(q)] = g.cellDz[size_t(q)] * U.len;
        }
        cc.geoTop = gr.tops.empty() ? 0.0 : gr.tops[0] * U.len;
        if (!g.corners.empty()) cc.geoTop = double(g.corners[0][0].z) * U.len;
        // C25: corner-point / dipping / irregular grids -> isoparametric hexahedra on the corner-point nodes
        if (!geoseq::tensorBox(g)) {
            cc.geoNodes = geoseq::gridNodes(g, U.len);
            out.notes.push_back("geomechanics on a corner-point / irregular grid: isoparametric hexahedra on the averaged corner-point nodes");
        }
        cc.geoCell = g.globalOf;
        for (int& q : cc.geoCell) if (g.isMincShell(q)) q = -1;   // MINC inner shells: no own mechanics cell (the box holds the outer shell)
        if (gr.dual) out.notes.push_back("geomechanics on a dual-porosity grid: the mechanics box spans both continua");
        // OPM arrays are per global cell: keep them aligned with the box
        auto fix = [&](std::vector<double>& v) { if (v.size() > 1 && v.size() != size_t(gr.size())) v.resize(1); };
        fix(cc.geo.E); fix(cc.geo.nu); fix(cc.geo.biot); fix(cc.geo.alphaT);
    }
    // heat loss: top face of the first and bottom face of the last active cell of every column
    if (TI.heatLoss) {
        cc.heatLoss = true;
        if (TI.obCv > 0) cc.obCv = TI.obCv * U.heatcr;
        if (TI.obK > 0) cc.obK = TI.obK * U.thcon;
        for (int j = 0; j < gr.ny; ++j)
            for (int i = 0; i < gr.nx; ++i) {
                int top = -1, bot = -1, gt = -1, gb = -1;
                for (int k = 0; k < gr.nz; ++k) {
                    const int q = gr.index(i, j, k);
                    const int a = q < int(g.activeOf.size()) ? g.activeOf[size_t(q)] : -1;
                    if (a < 0) continue;
                    if (top < 0) { top = a; gt = q; }
                    bot = a; gb = q;
                }
                if (top < 0) continue;
                cc.lossFaces.push_back({top, g.cellDx[size_t(gt)] * g.cellDy[size_t(gt)] * U.len * U.len});
                if (bot != top || gr.nz == 1) cc.lossFaces.push_back({bot, g.cellDx[size_t(gb)] * g.cellDy[size_t(gb)] * U.len * U.len});
            }
    }
    // ---- rock
    cc.cr = c.props.rock.cr * U.comp; cc.pref = c.props.rock.pref * U.p; cc.ct = 0;
    auto cellVals = [&](const std::vector<double>& v, double f, std::vector<double>& dst) {
        dst.clear();
        if (v.empty()) return false;
        if (v.size() == 1) { dst.assign(1, v[0] * f); return true; }
        dst.resize(size_t(na));
        for (int a = 0; a < na; ++a) { const int gq = hostQ(g.globalOf[size_t(a)]); const double x = gq < int(v.size()) ? v[size_t(gq)] : NaN; dst[size_t(a)] = std::isfinite(x) ? x * f : 0.0; }
        return true;
    };
    if (!cellVals(TI.heatcr, U.heatcr, cc.rockCvCell)) out.notes.push_back("HEATCR missing: rock heat capacity 2350 kJ/m3/K");
    cellVals(TI.heatcrt, U.heatcr / U.absT, cc.rockCvTCell);
    if (!cellVals(TI.thconr, U.thcon, cc.thconr)) out.notes.push_back("THCONR missing: conductivity from phase / rock defaults");
    cellVals(TI.thconsf, 1.0, cc.thconsfCell);
    if (cc.thconsfCell.empty()) cc.thconsf = 0;   // E300: THCONSF defaults to 0
    // ---- saturation functions (SWOF / SGOF of SATNUM region 1)
    {
        const SatRegion R = c.props.satRegion(0);
        for (const SatRow& r : R.swof) { cc.swT.push_back(r.s); cc.krwT.push_back(r.krw); cc.krowT.push_back(r.kro); }
        cc.swc = R.swof.empty() ? 0.0 : R.swof.front().s;
        for (size_t i = R.sgof.size(); i-- > 0;) { cc.slT.push_back(1.0 - R.sgof[i].s); cc.krgT.push_back(R.sgof[i].krw); cc.krogT.push_back(R.sgof[i].kro); }
        // capillary pressure (C8): Pcow(Sw) from SWOF, Pcgo(Sl) from SGOF (aligned with slT)
        bool pc = false;
        for (const SatRow& r : R.swof) if (r.pc != 0) pc = true;
        for (const SatRow& r : R.sgof) if (r.pc != 0) pc = true;
        if (pc) {
            for (const SatRow& r : R.swof) cc.pcowT.push_back(r.pc * U.p);
            for (size_t i = R.sgof.size(); i-- > 0;) cc.pcgoT.push_back(R.sgof[i].pc * U.p);
            if (cc.pcowT.size() != cc.swT.size()) cc.pcowT.clear();
            if (cc.pcgoT.size() != cc.slT.size()) cc.pcgoT.clear();
        }
        if (c.props.numSatRegions() > 1) out.notes.push_back("several SATNUM regions: table 1 is used by the thermal core");
        if (cc.swT.empty() && cc.slT.empty()) out.notes.push_back("no SWOF / SGOF: Corey defaults");
    }
    // ---- initial state
    const size_t ng = size_t(gr.nx) * size_t(gr.ny) * size_t(gr.nz);
    cc.T0.resize(size_t(na)); cc.p0.resize(size_t(na)); cc.sw0.resize(size_t(na)); cc.sg0.assign(size_t(na), 0.0);
    std::vector<double> tvdD, tvdT;
    for (auto& r : TI.tempvd) { tvdD.push_back(r.first * U.len); tvdT.push_back(U.tempK(r.second)); }
    for (int a = 0; a < na; ++a) {
        const int gq = hostQ(g.globalOf[size_t(a)]);
        double Tk = cc.stdT;
        if (!TI.tempi.empty()) { const double t = TI.tempi.size() == 1 ? TI.tempi[0] : (size_t(gq) < TI.tempi.size() ? TI.tempi[size_t(gq)] : NaN); if (std::isfinite(t)) Tk = U.tempK(t); }
        else if (!tvdD.empty()) Tk = lin(tvdD, tvdT, cc.cellDepth[size_t(a)]);
        cc.T0[size_t(a)] = Tk;
    }
    // compositions
    if (!TI.xmf.empty() && TI.xmf.size() == ng * nc) {
        cc.x0.resize(size_t(na));
        for (int a = 0; a < na; ++a) { std::vector<double> x(nc); for (size_t i = 0; i < nc; ++i) x[i] = TI.xmf[i * ng + size_t(hostQ(g.globalOf[size_t(a)]))]; cc.x0[size_t(a)] = x; }
        if (!TI.ymf.empty() && TI.ymf.size() == ng * nc) {
            cc.y0.resize(size_t(na));
            for (int a = 0; a < na; ++a) { std::vector<double> y(nc); for (size_t i = 0; i < nc; ++i) y[i] = TI.ymf[i * ng + size_t(hostQ(g.globalOf[size_t(a)]))]; cc.y0[size_t(a)] = y; }
        }
    } else if (TI.zi.size() == nc) cc.x0.assign(1, TI.zi);
    else if (nc == 1) cc.x0.assign(1, std::vector<double>{1.0});
    else return fail("initial oil composition is missing (ZI or XMF)");
    for (auto& x : cc.x0) { double s = 0; for (double v : x) s += v; if (s > 0) for (double& v : x) v /= s; else return fail("initial oil composition sums to zero"); }
    for (auto& y : cc.y0) { double s = 0; for (double v : y) s += v; if (s > 0) for (double& v : y) v /= s; else y = cc.x0.size() == 1 ? cc.x0[0] : y; }
    if (c.init.active()) {
        for (int a = 0; a < na; ++a) {
            const size_t q = size_t(hostQ(g.globalOf[size_t(a)]));
            auto pick = [&](const std::vector<double>& v, double d) { return v.empty() ? d : v[std::min(q, v.size() - 1)]; };
            cc.p0[size_t(a)] = pick(c.init.pressure, 0.0) * U.p;
            cc.sw0[size_t(a)] = pick(c.init.swat, cc.swc);
            cc.sg0[size_t(a)] = pick(c.init.sgas, 0.0);
        }
    } else {   // EQUIL: hydrostatic oil / water columns (sharp contact; Pcow = 0)
        const EquilInput& E = c.equil;
        const double zD = E.datumDepth * U.len, pD = E.datumPressure * U.p, owc = E.owc * U.len;
        const std::vector<double>& x = cc.x0.front();
        const double Td = tvdD.empty() ? cc.T0.front() : lin(tvdD, tvdT, zD);
        double m = 0, v = 0;
        for (size_t i = 0; i < nc; ++i) {
            const Component& k = cc.comps[i];
            const double dT = Td - cc.tRef;
            const double rho = k.rhoRef / ((1 - k.cpL * (pD - cc.pRef)) * (1 + k.ct1 * dT));
            m += x[i] * k.mw; v += x[i] * k.mw / rho;
        }
        const double rhoO = v > 0 ? m / v : 900.0;
        const double pOwc = pD + rhoO * GRAV * (owc - zD);
        const double rhoW = 1.0 / steam::liquid(Td, std::max(pOwc, 1e5)).v;
        for (int a = 0; a < na; ++a) {
            const double z = cc.cellDepth[size_t(a)];
            if (z > owc) { cc.p0[size_t(a)] = pOwc + rhoW * GRAV * (z - owc); cc.sw0[size_t(a)] = 1.0; }
            else { cc.p0[size_t(a)] = pD + rhoO * GRAV * (z - zD); cc.sw0[size_t(a)] = cc.swc; }
        }
        out.notes.push_back("EQUIL: hydrostatic oil / water columns with a sharp contact (no capillary transition)");
    }
    // ---- wells: connections (Peaceman) and controls per report step
    const double rhoWsc = 1.0 / steam::liquid(cc.stdT, std::max(cc.stdP, 1e3)).v;
    const double molPerSm3W = rhoWsc / steam::kMw;
    const size_t nw = c.wells.size();
    for (const auto& w : c.wells) out.wellNames.push_back(w.name);
    std::vector<WellControl> ctrl(nw);
    std::vector<bool> haveCtrl(nw, false);
    bool notedGas = false, notedCtl = false, notedCold = false;
    for (size_t s = 0; s < c.schedule.size(); ++s) {
        const ReportStep& rs = c.schedule[s];
        out.stepDays.push_back(rs.length);
        for (const WellControl& wc : rs.controls)
            for (size_t w = 0; w < nw; ++w) if (c.wells[w].name == wc.well) { ctrl[w] = wc; haveCtrl[w] = true; }
        std::vector<Well> ws;
        for (size_t w = 0; w < nw; ++w) {
            const WellInput& in = c.wells[w];
            Well tw; tw.name = in.name;
            tw.refDepth = in.refDepth > 0 ? in.refDepth * U.len : -1;
            for (const Completion& cp : in.completions) {
                if (int(s) < std::max(0, cp.fromStep) || (cp.untilStep >= 0 && int(s) >= cp.untilStep) || !cp.open) continue;
                for (int k = std::min(cp.k1, cp.k2); k <= std::max(cp.k1, cp.k2); ++k) {
                    const int kk = (gr.dual == 1 && k < gr.nzMatrix()) ? k + gr.nzMatrix() : k;
                    std::vector<int> targets{cp.lgr >= 0 ? g.localCell(cp.lgr, cp.i, cp.j, k) : gr.index(cp.i, cp.j, kk)};
                    if (cp.lgr < 0 && g.isHost(targets[0])) targets = g.hostWellColumn(targets[0], cp.dir);
                    for (int q : targets) {
                        const int a = q < int(g.activeOf.size()) ? g.activeOf[size_t(q)] : -1;
                        if (a < 0) continue;
                        double wi = cp.cf;
                        if (wi <= 0) {
                            const double kx = q < int(gr.permx.size()) ? g.effPerm(gr.permx, q) : 0.0;
                            const double ky = q < int(gr.permy.size()) && !std::isnan(gr.permy[size_t(q)]) ? g.effPerm(gr.permy, q) : kx,
                                         kz = q < int(gr.permz.size()) && !std::isnan(gr.permz[size_t(q)]) ? g.effPerm(gr.permz, q) : kx;
                            const double DX = g.cellDx[size_t(q)], DY = g.cellDy[size_t(q)], DZ = g.cellDz[size_t(q)];
                            double k1, k2, d1, d2, hh;
                            if (cp.dir == 'X' || cp.dir == 'x') { k1 = ky; k2 = kz; d1 = DY; d2 = DZ; hh = DX; }
                            else if (cp.dir == 'Y' || cp.dir == 'y') { k1 = kx; k2 = kz; d1 = DX; d2 = DZ; hh = DY; }
                            else { k1 = kx; k2 = ky; d1 = DX; d2 = DY; hh = DZ * (q < int(gr.ntg.size()) ? gr.ntg[size_t(q)] : 1.0); }
                            const double kh = cp.kh > 0 ? cp.kh : std::sqrt(k1 * k2) * hh;
                            const double r12 = std::sqrt(k2 / k1), r21 = std::sqrt(k1 / k2);
                            const double ro = 0.28 * std::sqrt(r12 * d1 * d1 + r21 * d2 * d2) / (std::pow(k2 / k1, 0.25) + std::pow(k1 / k2, 0.25));
                            wi = Units::of(c.units).cDarcy * 2.0 * M_PI * kh / (std::log(ro / (0.5 * cp.diameter)) + cp.skin);
                        }
                        wi *= in.piMult * cp.mult;
                        tw.cells.push_back(a); tw.wi.push_back(wi * U.trans);
                    }
                }
            }
            const WellControl& wc = ctrl[w];
            tw.open = haveCtrl[w] && wc.status == WellStatus::Open && !tw.cells.empty();
            tw.injector = in.type != WellType::Producer;
            if (tw.injector) {
                tw.zinj.assign(nc + 1, 0.0); tw.zinj[nc] = 1.0;
                const int ph = wc.phase ? wc.phase : (in.type == WellType::WaterInjector ? 1 : 2);
                bool gasInj = false;
                if (ph != 1) {
                    // gas injector (air / O2 / N2 / CO2 ...): WINJGAS stream in force (WELLSTRE composition of the nc components)
                    std::string sname;
                    for (size_t s2 = 0; s2 <= s; ++s2)
                        for (auto& ws : c.schedule[s2].injStreams) if (matchWell(ws.first, in.name)) sname = up(ws.second);
                    const CompStream* st = nullptr;
                    for (auto& x : TI.streams) if (up(x.name) == sname) st = &x;
                    double zs = 0;
                    // solids cannot be injected: their stream fractions are dropped (and noted)
                    auto isSolid = [&](size_t i2) { return i2 < TI.comps.size() && TI.comps[i2].type == 2; };
                    if (st) for (size_t i2 = 0; i2 < nc && i2 < st->z.size(); ++i2) if (!isSolid(i2)) zs += std::max(0.0, st->z[i2]);
                    if (st) for (size_t i2 = 0; i2 < nc && i2 < st->z.size(); ++i2)
                        if (isSolid(i2) && st->z[i2] > 0) { out.notes.push_back("stream " + st->name + ": solid component " + TI.comps[i2].name + " ignored"); break; }
                    if (st && zs > 0) {
                        for (size_t i2 = 0; i2 < nc; ++i2) tw.zinj[i2] = i2 < st->z.size() && !isSolid(i2) ? std::max(0.0, st->z[i2]) / zs : 0.0;
                        tw.zinj[nc] = 0.0;
                        gasInj = true;
                    } else {
                        tw.open = false;
                        if (!notedGas) { out.notes.push_back("gas injector " + in.name + " without a WINJGAS / WELLSTRE stream: shut (steam / hot water: WATER injectors + WINJTEMP)"); notedGas = true; }
                    }
                }
                tw.bhpMax = wc.bhp > 0 ? wc.bhp * U.p : 0.0;
                if (wc.mode == Control::BHP) { tw.bhpControl = true; tw.bhpMax = wc.bhp * U.p; }
                else if (gasInj) tw.rate = std::max(0.0, wc.rate) * U.gas * cc.stdP / (8.314462618 * cc.stdT) / 86400.0;   // surface gas -> mol/s
                else tw.rate = std::max(0.0, wc.rate) * U.liq * molPerSm3W / 86400.0;   // cold-water equivalent
                // WINJTEMP in force (last record of a step <= s matching the well)
                const ThermalWinj* wt = nullptr;
                for (const ThermalWinj& x : TI.winj) if (x.step <= int(s) && matchWell(x.well, in.name)) wt = &x;
                if (wt && gasInj) {
                    tw.quality = -1;
                    tw.tinj = wt->temp > -1e8 ? U.tempK(wt->temp) : cc.stdT;
                } else if (gasInj) {
                    tw.quality = -1; tw.tinj = cc.stdT;
                } else if (wt) {
                    tw.quality = wt->quality;
                    tw.tinj = wt->temp > -1e8 ? U.tempK(wt->temp) : 0.0;
                    tw.pinj = wt->pres > 0 ? wt->pres * U.p : 0.0;
                    tw.hinj = wt->enth > 0 ? wt->enth * U.hv : 0.0;
                    if (tw.quality < 0 && tw.tinj <= 0) tw.tinj = cc.stdT;
                } else {
                    tw.quality = -1; tw.tinj = cc.stdT;
                    if (tw.open && !notedCold) { out.notes.push_back("injector " + in.name + " without WINJTEMP: cold water at the standard temperature"); notedCold = true; }
                }
            } else {
                tw.bhp = (wc.bhp > 0 ? wc.bhp : U.stdP()) * U.p;
                switch (wc.mode) {
                    case Control::ORAT: tw.control = 1; tw.target = wc.orat * U.liq / 86400.0; break;
                    case Control::WRAT: tw.control = 2; tw.target = wc.wrat * U.liq / 86400.0; break;
                    case Control::LRAT: tw.control = 3; tw.target = wc.lrat * U.liq / 86400.0; break;
                    case Control::BHP: tw.control = 0; break;
                    default:
                        tw.control = 0;
                        if (!notedCtl) { out.notes.push_back(std::string("producer control ") + toString(wc.mode) + " is not modelled by the thermal core - BHP control"); notedCtl = true; }
                        break;
                }
                if (tw.control != 0 && !(tw.target > 0)) tw.control = 0;
                for (const ThermalTrap& tp : c.thermal.traps)   // last STEAMTRAP record in force for this well
                    if (tp.step <= int(s) && (tp.well == in.name || tp.well == "*")) { tw.subcool = tp.subcool > 0 ? tp.subcool * (field ? 5.0 / 9.0 : 1.0) : 0.0; tw.steamMax = tp.steamMax > 0 ? tp.steamMax * U.liq / 86400.0 : 0.0; }
            }
            ws.push_back(tw);
        }
        out.stepWells.push_back(ws);
    }
    if (out.stepWells.empty()) return fail("the schedule has no report steps");
    cc.wells = out.stepWells.front();
    // C12: Langmuir sorption appends site / sorbed solids before water (every schedule step's injection streams too)
    for (const ThermalSorb& so : TI.sorb) {
        Sorption S; S.comp = so.comp; S.capacity = so.capacity; S.bRef = so.b; S.tRef = so.tRef; S.dH = so.dH; S.tau = so.tau;
        S.phase = so.phase; S.adsRho = so.adsRho; S.site = so.site; S.initEquilibrium = so.initEq;
        const size_t before = cc.comps.size();
        std::string e;
        if (!addSorption(cc, S, &e)) return fail(e);
        const size_t added = cc.comps.size() - before;
        for (auto& ws : out.stepWells)
            for (auto& w : ws)
                if (w.zinj.size() >= before) w.zinj.insert(w.zinj.begin() + long(before), added, 0.0);
    }
    if (!TI.sorb.empty()) cc.wells = out.stepWells.front();
    // C12: embedded discrete fractures (all-active Cartesian box: the cell order is the box order)
    if (!TI.fractures.empty()) {
        if (!gr.zcorn.empty() || na != gr.size()) return fail("FRACTURE (EDFM) needs a Cartesian grid with all cells active");
        cc.nx = gr.nx; cc.ny = gr.ny; cc.nz = gr.nz;
        cc.dx.resize(size_t(na)); cc.dy.resize(size_t(na)); cc.dz.resize(size_t(na));
        cc.permx.resize(size_t(na)); cc.permy.resize(size_t(na));
        for (int q = 0; q < na; ++q) {
            cc.dx[size_t(q)] = g.cellDx[size_t(q)] * U.len; cc.dy[size_t(q)] = g.cellDy[size_t(q)] * U.len; cc.dz[size_t(q)] = g.cellDz[size_t(q)] * U.len;
            cc.permx[size_t(q)] = permDir(q, 'X') * mD; cc.permy[size_t(q)] = permDir(q, 'Y') * mD;
        }
        cc.top = gr.tops.empty() ? 0.0 : gr.tops[0] * U.len;
        std::vector<Fracture> fr;
        for (const ThermalFracture& f : TI.fractures) {
            Fracture F; F.name = f.name; F.x0 = f.x0; F.y0 = f.y0; F.x1 = f.x1; F.y1 = f.y1; F.zTop = f.zTop; F.zBot = f.zBot;
            F.aperture = f.aperture; F.perm = f.perm; F.mech = f.mech; F.kn0 = f.kn0; F.vm = f.vm; F.mu = f.mu; F.ks = f.ks; F.dilationDeg = f.dilationDeg;
            fr.push_back(F);
        }
        FractureStress st; st.sHmaxGrad = TI.fracStress.sHmaxGrad; st.sHminGrad = TI.fracStress.sHminGrad; st.azimuthDeg = TI.fracStress.azimuthDeg;
        st.E = TI.fracStress.E; st.nu = TI.fracStress.nu; st.alphaT = TI.fracStress.alphaT;
        std::string e;
        if (!addFractures(cc, fr, st, &e)) return fail(e);
        int nwf = 0;
        for (auto& ws : out.stepWells) nwf = connectWellsToFractures(cc, ws);
        cc.wells = out.stepWells.front();
        if (nwf > 0) out.notes.push_back("EDFM: " + std::to_string(nwf) + " well - fracture connections");
        out.notes.push_back("EDFM: " + std::to_string(cc.fracCells.size()) + " fracture cells appended after the matrix cells");
    }
    if (!c.aquifers.empty()) out.notes.push_back("analytic aquifers are not modelled by the thermal core");
    return true;
}

bool deckModelFromFile(const std::string& path, DeckModel& out, std::string* error) {
    try {
        const ressim::Case c = importEclipseDeck(path);
        return deckModel(c, out, error);
    } catch (const std::exception& e) {
        if (error) *error = e.what();
        return false;
    }
}

// ======================================================================== run
bool runDeckModel(const DeckModel& m, std::vector<Summary>& out, const std::function<bool(const Summary&, Simulator&)>& cb,
                  std::string* error, int maxSteps, Stats* stats, const DeckRestart* rst) {
    const TU U(m.field);
    Simulator sim;
    std::string err;
    if (!sim.setup(m.base, &err)) { if (error) *error = err; return false; }
    out.clear();
    const size_t nw = m.wellNames.size();
    std::vector<WellReport> prev(nw);
    double t = 0, lastDays = 0;
    auto summary = [&](double days) {
        Summary S; S.time = days;
        const auto& R = sim.wellReports();
        const Case& b = sim.model();
        double pv = 0, pp = 0, tt = 0, sgv = 0, tmax = 0;
        for (int q = 0; q < sim.size(); ++q) {
            const double v = b.cellPV[size_t(q)];
            pv += v; pp += v * sim.pressure()[size_t(q)]; tt += v * sim.temperature()[size_t(q)]; sgv += v * sim.sg()[size_t(q)];
            tmax = std::max(tmax, sim.temperature()[size_t(q)]);
        }
        S.fpr = pp / std::max(pv, 1e-30) / U.p;
        S.favgT = U.tempDeck(tt / std::max(pv, 1e-30));
        S.fmaxT = U.tempDeck(tmax);
        S.fsteamPV = sgv / std::max(pv, 1e-30);
        S.foip = sim.stockTankOilInPlace() / U.liq;
        S.fhlt = sim.cumHeatLoss() / U.energy();
        S.fenergyIn = sim.cumEnthalpyIn() / U.energy();
        S.wells.resize(nw);
        for (size_t w = 0; w < nw && w < R.size(); ++w) {
            const WellReport& r = R[w];
            const Well& wl = sim.model().wells[w];
            Summary::Well& o = S.wells[w];
            o.open = wl.open; o.onBhp = r.onBhp; o.bhp = r.bhp / U.p;
            o.hrate = r.hRate * 86400.0 / U.energy();
            o.tfluid = r.tFluid > 0 ? U.tempDeck(r.tFluid) : 0.0;
            o.molar.assign(r.cumMolP.size(), 0.0);
            if (days > lastDays)
                for (size_t i = 0; i < r.cumMolP.size(); ++i) {
                    const double p0 = i < prev[w].cumMolP.size() ? prev[w].cumMolP[i] : 0.0, i0 = i < prev[w].cumMolI.size() ? prev[w].cumMolI[i] : 0.0;
                    o.molar[i] = ((r.cumMolP[i] - p0) - (r.cumMolI[i] - i0)) / ((days - lastDays) * 86400.0);
                }
            if (wl.injector) {
                o.winj = r.qcwe * 86400.0 / U.liq; S.fwir += o.winj; S.fwit += r.cumCwe / U.liq;
                // injected non-water moles -> surface gas volume (ideal gas at standard conditions)
                const double vStd = 8.314462618 * 288.71 / 101325.0;   // m3/mol at 15.56 C (60 F), 1 atm
                const size_t nci = r.cumMolI.empty() ? 0 : r.cumMolI.size() - 1;
                double cum = 0, rate = 0;
                for (size_t i = 0; i < nci; ++i) {
                    cum += r.cumMolI[i];
                    const double i0 = i < prev[w].cumMolI.size() ? prev[w].cumMolI[i] : 0.0;
                    if (days > lastDays) rate += (r.cumMolI[i] - i0) / ((days - lastDays) * 86400.0);
                }
                S.fgir += rate * vStd * 86400.0 / U.gas; S.fgit += cum * vStd / U.gas;
            }
            else {
                o.oil = r.qo * 86400.0 / U.liq; o.water = r.qw * 86400.0 / U.liq; o.gas = r.qg * 86400.0 / U.gas;
                S.fopr += o.oil; S.fwpr += o.water; S.fgpr += o.gas;
                S.fopt += r.cumO / U.liq; S.fwpt += r.cumW / U.liq; S.fgpt += r.cumG / U.gas;
            }
            prev[w].cumMolP = r.cumMolP; prev[w].cumMolI = r.cumMolI;
        }
        S.csor = S.fopt > 0 ? S.fwit / S.fopt : 0.0;
        lastDays = days;
        return S;
    };
    size_t s0 = 0;
    if (rst && rst->startStep > 0) {
        try {
            StateIn in(rst->state);
            if (in.r.str() != "RESSIM-THERMAL-RUN-1") throw std::runtime_error("not a thermal run restart");
            in(t);
            s0 = size_t(rst->startStep);
            if (s0 > m.stepDays.size()) throw std::runtime_error("restart step beyond the schedule");
            if (s0 - 1 < m.stepWells.size()) sim.wells() = m.stepWells[s0 - 1];   // controls of the restart step (change detection)
            sim.loadState(in.r.str());
            const auto& R0 = sim.wellReports();
            for (size_t w = 0; w < nw && w < R0.size(); ++w) { prev[w].cumMolP = R0[w].cumMolP; prev[w].cumMolI = R0[w].cumMolI; }
            lastDays = t / 86400.0;
            Summary S = summary(t / 86400.0);
            out.push_back(S);
        } catch (const std::exception& e) { if (error) *error = std::string("restart: ") + e.what(); return false; }
    } else {
        Summary S0 = summary(0);
        out.push_back(S0);
        if (cb && !cb(S0, sim)) return true;
    }
    const size_t ns = m.stepDays.size();
    for (size_t s = s0; s < ns && (maxSteps < 0 || int(s) < maxSteps); ++s) {
        if (s < m.stepWells.size()) {
            std::vector<Well>& ws = sim.wells();
            bool change = false;
            for (size_t w = 0; w < ws.size() && w < m.stepWells[s].size(); ++w) {
                const Well& a = ws[w]; const Well& b = m.stepWells[s][w];
                if (a.open != b.open || a.control != b.control || a.subcool != b.subcool || a.steamMax != b.steamMax || a.bhpControl != b.bhpControl || a.rate != b.rate || a.target != b.target ||
                    a.bhp != b.bhp || a.bhpMax != b.bhpMax || a.tinj != b.tinj || a.quality != b.quality) change = true;
                ws[w] = m.stepWells[s][w];
            }
            if (change && s > 0) sim.limitNextDt(0.05 * 86400.0);   // well opened / shut / re-targeted: restart with a small step
        }
        const double len = m.stepDays[s] * 86400.0;
        t += len;
        if (!sim.advanceTo(t, std::min(len, 30.0 * 86400.0))) {
            if (error) *error = "thermal core: time step failed at day " + std::to_string(sim.time() / 86400.0);
            if (stats) *stats = sim.stats();
            return false;
        }
        Summary S = summary(t / 86400.0);
        out.push_back(S);
        if (rst && rst->wants(int(s + 1), s + 1 == ns)) {
            StateOut o; o.w.str("RESSIM-THERMAL-RUN-1"); o(t); o.w.str(sim.saveState());
            rst->onSave(int(s + 1), t / 86400.0, o.w.buf);
        }
        if (cb && !cb(S, sim)) break;
    }
    if (stats) *stats = sim.stats();
    return true;
}

}  // namespace thermal
}  // namespace ressim
