// C14: user-defined result plots (custom plot library). See custom_plots.h for the curve key grammar.
#include "app.h"
#include "icons.h"
#include "widgets.h"
#include "implot.h"
#include "plotkit.h"
#include "custom_plots.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>

using namespace ressim;
namespace fs = std::filesystem;

// ============================================================================ library file
std::string customPlotsToText(const std::vector<CustomPlot>& plots) {
    std::ostringstream o;
    o << "# ResSim Studio custom plots 1\n";
    auto clean = [](std::string s) { for (char& c : s) if (c == '\t' || c == '\n' || c == '\r') c = ' '; return s; };
    for (const CustomPlot& p : plots) {
        o << "PLOT\t" << clean(p.name) << "\n";
        o << "TAB\t" << (p.ownTab ? 1 : 0) << "\nLOGY\t" << (p.logY ? 1 : 0) << "\n";
        if (!p.ylabel.empty()) o << "YLABEL\t" << clean(p.ylabel) << "\n";
        if (!p.y2label.empty()) o << "Y2LABEL\t" << clean(p.y2label) << "\n";
        for (const CustomCurve& c : p.curves) {
            char sc[40]; std::snprintf(sc, sizeof(sc), "%.10g", c.scale);
            o << "CURVE\t" << c.axis << "\t" << c.style << "\t" << sc << "\t" << clean(c.key) << "\t" << clean(c.label) << "\n";
        }
        o << "END\n";
    }
    return o.str();
}

std::vector<CustomPlot> customPlotsFromText(const std::string& text) {
    std::vector<CustomPlot> out;
    std::istringstream in(text);
    std::string line;
    CustomPlot* cur = nullptr;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> f;
        size_t a = 0;
        for (;;) { const size_t b = line.find('\t', a); f.push_back(line.substr(a, b == std::string::npos ? std::string::npos : b - a)); if (b == std::string::npos) break; a = b + 1; }
        const std::string& k = f[0];
        if (k == "PLOT") { out.emplace_back(); cur = &out.back(); cur->name = f.size() > 1 ? f[1] : "Plot"; continue; }
        if (!cur) continue;
        try {
            if (k == "TAB" && f.size() > 1) cur->ownTab = std::stoi(f[1]) != 0;
            else if (k == "LOGY" && f.size() > 1) cur->logY = std::stoi(f[1]) != 0;
            else if (k == "YLABEL" && f.size() > 1) cur->ylabel = f[1];
            else if (k == "Y2LABEL" && f.size() > 1) cur->y2label = f[1];
            else if (k == "CURVE" && f.size() >= 5) {
                CustomCurve c;
                c.axis = std::clamp(std::stoi(f[1]), 0, 1);
                c.style = std::clamp(std::stoi(f[2]), 0, 3);
                c.scale = std::stod(f[3]);
                c.key = f[4];
                if (f.size() > 5) c.label = f[5];
                cur->curves.push_back(c);
            } else if (k == "END") cur = nullptr;
        } catch (...) {}
    }
    return out;
}

bool writeCustomPlots(const std::string& path, const std::vector<CustomPlot>& plots) {
    std::ofstream f(fs::u8path(path), std::ios::binary);
    if (!f) return false;
    f << customPlotsToText(plots);
    return bool(f);
}

bool readCustomPlots(const std::string& path, std::vector<CustomPlot>& plots) {
    std::ifstream f(fs::u8path(path), std::ios::binary);
    if (!f) return false;
    std::stringstream ss; ss << f.rdbuf();
    plots = customPlotsFromText(ss.str());
    return true;
}

// ============================================================================ curve catalog / resolution
namespace {
enum UnitKind { U_LR, U_GR, U_LV, U_GV, U_P, U_FRAC, U_GOR, U_KG, U_T, U_EN, U_MW, U_GJ, U_MM, U_DAY, U_N, U_PM, U_PCT, U_MOLS, U_MOL, U_KGD, U_RV };

std::string unitText(int u, const RunResults& R) {
    const Units U = Units::of(R.input.units);
    const bool fld = R.input.units == UnitSystem::Field;
    switch (u) {
        case U_LR: return std::string(U.liq) + "/d";
        case U_GR: return std::string(U.gas) + "/d";
        case U_LV: return U.liq;
        case U_GV: return U.gas;
        case U_P: return U.p;
        case U_FRAC: return "-";
        case U_GOR: return U.rs;
        case U_KG: return "kg";
        case U_T: return fld ? "degF" : "degC";
        case U_EN: return fld ? "Btu" : "kJ";
        case U_MW: return "MW";
        case U_GJ: return "GJ";
        case U_MM: return "mm";
        case U_DAY: return "days";
        case U_PM: return fld ? "lb" : "kg";
        case U_PCT: return "%";
        case U_MOLS: return "mol/s";
        case U_MOL: return "mol";
        case U_KGD: return "kg/d";
        case U_RV: return fld ? "rb/d" : "rm3/d";
        default: return "";
    }
}

// C27: tracer runs (black-oil tracer module) reuse the component vectors - amounts instead of mol / kg
bool tracerRun(const RunResults& R) { return !R.components.empty() && R.components[0].kind == 4; }
std::string compUnitText(int u, const RunResults& R) {
    if (tracerRun(R)) {
        if (u == U_MOLS) return "amount/s";
        if (u == U_MOL || u == U_KG) return "amount";
        if (u == U_KGD) return "amount/d";
    }
    return unitText(u, R);
}

using SP = SummaryPoint;
struct FDef { const char* name; const char* desc; int unit; bool rate; double (*get)(const SP&); bool always; };
const FDef kField[] = {
    {"FOPR", "oil production rate", U_LR, true, [](const SP& s) { return s.fopr; }, true},
    {"FWPR", "water production rate", U_LR, true, [](const SP& s) { return s.fwpr; }, true},
    {"FGPR", "gas production rate", U_GR, true, [](const SP& s) { return s.fgpr; }, true},
    {"FLPR", "liquid production rate", U_LR, true, [](const SP& s) { return s.fopr + s.fwpr; }, true},
    {"FWIR", "water injection rate", U_LR, true, [](const SP& s) { return s.fwir; }, true},
    {"FGIR", "gas injection rate", U_GR, true, [](const SP& s) { return s.fgir; }, true},
    {"FOPT", "cumulative oil production", U_LV, false, [](const SP& s) { return s.fopt; }, true},
    {"FWPT", "cumulative water production", U_LV, false, [](const SP& s) { return s.fwpt; }, true},
    {"FGPT", "cumulative gas production", U_GV, false, [](const SP& s) { return s.fgpt; }, true},
    {"FWIT", "cumulative water injection", U_LV, false, [](const SP& s) { return s.fwit; }, true},
    {"FGIT", "cumulative gas injection", U_GV, false, [](const SP& s) { return s.fgit; }, true},
    {"FPR", "average pressure (HCPV weighted)", U_P, false, [](const SP& s) { return s.fpr; }, true},
    {"FPRP", "average pressure (PV weighted)", U_P, false, [](const SP& s) { return s.fprp; }, false},
    {"FWCT", "water cut", U_FRAC, true, [](const SP& s) { return s.fwct; }, true},
    {"FGOR", "gas-oil ratio", U_GOR, true, [](const SP& s) { return s.fgor; }, true},
    {"FOIP", "oil in place", U_LV, false, [](const SP& s) { return s.oip; }, true},
    {"FWIP", "water in place", U_LV, false, [](const SP& s) { return s.wip; }, true},
    {"FGIP", "gas in place", U_GV, false, [](const SP& s) { return s.gip; }, true},
    {"FAQR", "aquifer influx rate (surface)", U_LR, true, [](const SP& s) { return s.faqrS; }, false},
    {"FAQT", "cumulative aquifer influx (surface)", U_LV, false, [](const SP& s) { return s.faqtS; }, false},
    {"FCIR", "polymer injection rate", U_PM, true, [](const SP& s) { return s.fcir; }, false},
    {"FCPR", "polymer production rate", U_PM, true, [](const SP& s) { return s.fcpr; }, false},
    {"FCIT", "cumulative polymer injected", U_PM, false, [](const SP& s) { return s.fcit; }, false},
    {"FCPT", "cumulative polymer produced", U_PM, false, [](const SP& s) { return s.fcpt; }, false},
    {"FCIP", "polymer in solution", U_PM, false, [](const SP& s) { return s.fcip; }, false},
    {"FCAD", "polymer adsorbed", U_PM, false, [](const SP& s) { return s.fcad; }, false},
    {"FGLIR", "lift gas injection rate (gas lift)", U_GR, true, [](const SP& s) { return s.fglir; }, false},
    {"FGLIT", "cumulative lift gas injected", U_GV, false, [](const SP& s) { return s.fglit; }, false},
    {"FGCR", "gas consumption rate (GCONSUMP)", U_GR, true, [](const SP& s) { return s.fgcr; }, false},
    {"FGCT", "cumulative gas consumed", U_GV, false, [](const SP& s) { return s.fgct; }, false},
    {"FGIMR", "gas import rate (GCONSUMP)", U_GR, true, [](const SP& s) { return s.fgimr; }, false},
    {"FGIMT", "cumulative gas imported", U_GV, false, [](const SP& s) { return s.fgimt; }, false},
    {"FGSR", "gas sales rate (produced - injected - consumed + imported)", U_GR, true, [](const SP& s) { return s.fgsr; }, false},
    {"FGST", "cumulative gas sales", U_GV, false, [](const SP& s) { return s.fgst; }, false},
    {"FNIR", "solvent injection rate", U_GR, true, [](const SP& s) { return s.fnir; }, false},
    {"FNPR", "solvent production rate", U_GR, true, [](const SP& s) { return s.fnpr; }, false},
    {"FNIT", "cumulative solvent injected", U_GV, false, [](const SP& s) { return s.fnit; }, false},
    {"FNPT", "cumulative solvent produced", U_GV, false, [](const SP& s) { return s.fnpt; }, false},
    {"FNIP", "solvent in place", U_GV, false, [](const SP& s) { return s.fnip; }, false},
    {"FGMIP", "stored CO2 / tracked component, total", U_KG, false, [](const SP& s) { return s.co2Mob + s.co2Trap + s.co2Diss + s.co2Min + s.co2Oil; }, false},
    {"FGMMO", "CO2 mobile", U_KG, false, [](const SP& s) { return s.co2Mob; }, false},
    {"FGMST", "CO2 structurally held", U_KG, false, [](const SP& s) { return s.co2Struct; }, false},
    {"FGMTR", "CO2 residually trapped", U_KG, false, [](const SP& s) { return s.co2Trap; }, false},
    {"FGMDS", "CO2 dissolved in brine", U_KG, false, [](const SP& s) { return s.co2Diss; }, false},
    {"FGMOI", "CO2 dissolved in oil", U_KG, false, [](const SP& s) { return s.co2Oil; }, false},
    {"FGMMN", "CO2 mineralised", U_KG, false, [](const SP& s) { return s.co2Min; }, false},
    {"FGMIM", "CO2 immobile (krg = 0)", U_KG, false, [](const SP& s) { return s.co2KrImmob; }, false},
    {"FSALT", "precipitated salt", U_KG, false, [](const SP& s) { return s.saltSolid; }, false},
    {"FAVGT", "average temperature", U_T, false, [](const SP& s) { return s.favgT; }, false},
    {"FMAXT", "maximum temperature", U_T, false, [](const SP& s) { return s.fmaxT; }, false},
    {"FSTPV", "steam pore-volume fraction", U_FRAC, false, [](const SP& s) { return s.fsteamPV; }, false},
    {"FHLT", "cumulative heat loss", U_EN, false, [](const SP& s) { return s.fhlt; }, false},
    {"FHIN", "net injected enthalpy", U_EN, false, [](const SP& s) { return s.fenergyIn; }, false},
    {"FCSOR", "cumulative steam-oil ratio", U_FRAC, false, [](const SP& s) { return s.csor; }, false},
    {"FRHR", "reaction heat rate", U_MW, true, [](const SP& s) { return s.freactRate * 1e-6; }, false},
    {"FRHT", "cumulative reaction heat", U_GJ, false, [](const SP& s) { return s.freactCum * 1e-9; }, false},
    {"FUPL", "surface uplift", U_MM, false, [](const SP& s) { return s.fuplift * 1e3; }, false},
    {"FUPR", "reservoir-top uplift", U_MM, false, [](const SP& s) { return s.fupliftRes * 1e3; }, false},
    {"FMXF", "maximum failure index", U_FRAC, false, [](const SP& s) { return s.fmaxFail; }, false},
    {"DT", "time-step size", U_DAY, true, [](const SP& s) { return s.dt; }, true},
    {"NEWTON", "Newton iterations per step", U_N, true, [](const SP& s) { return double(s.newton); }, true},
    {"LINEAR", "linear iterations per step", U_N, true, [](const SP& s) { return double(s.linear); }, true},
};
const int kNumField = int(sizeof(kField) / sizeof(kField[0]));

using WR = WellRates;
struct WDef { const char* name; const char* desc; int unit; bool cum; double (*get)(const WR&); };
const WDef kWell[] = {
    {"WOPR", "oil production rate", U_LR, false, [](const WR& w) { return w.oil; }},
    {"WWPR", "water production rate", U_LR, false, [](const WR& w) { return w.water; }},
    {"WGPR", "gas production rate", U_GR, false, [](const WR& w) { return w.gas; }},
    {"WLPR", "liquid production rate", U_LR, false, [](const WR& w) { return w.oil + w.water; }},
    {"WWIR", "water injection rate", U_LR, false, [](const WR& w) { return w.waterInj; }},
    {"WGIR", "gas injection rate", U_GR, false, [](const WR& w) { return w.gasInj; }},
    {"WVPR", "reservoir voidage rate", U_RV, false, [](const WR& w) { return w.resv; }},
    {"WBHP", "bottom-hole pressure", U_P, false, [](const WR& w) { return w.open ? w.bhp : double(NAN); }},
    {"WTHPN", "network THP limit (node pressure)", U_P, false, [](const WR& w) { return w.thpNet >= 0 ? w.thpNet : double(NAN); }},
    {"WDPWB", "wellbore friction + ICD drop (deepest connection)", U_P, false, [](const WR& w) { return w.dpWellbore; }},
    {"WWCT", "water cut", U_FRAC, false, [](const WR& w) { const double l = w.oil + w.water; return l > 0 ? w.water / l : 0.0; }},
    {"WGOR", "gas-oil ratio", U_GOR, false, [](const WR& w) { return w.oil > 0 ? w.gas / w.oil : 0.0; }},
    {"WOPT", "cumulative oil production", U_LV, true, [](const WR& w) { return w.oil; }},
    {"WWPT", "cumulative water production", U_LV, true, [](const WR& w) { return w.water; }},
    {"WGPT", "cumulative gas production", U_GV, true, [](const WR& w) { return w.gas; }},
    {"WWIT", "cumulative water injection", U_LV, true, [](const WR& w) { return w.waterInj; }},
    {"WGIT", "cumulative gas injection", U_GV, true, [](const WR& w) { return w.gasInj; }},
    {"WCIR", "polymer injection rate", U_PM, false, [](const WR& w) { return w.polyInj; }},
    {"WCPR", "polymer production rate", U_PM, false, [](const WR& w) { return w.polyProd; }},
    {"WGLIR", "lift gas injection rate (gas lift)", U_GR, false, [](const WR& w) { return w.lift; }},
    {"WGLIT", "cumulative lift gas injected", U_GV, true, [](const WR& w) { return w.lift; }},
    {"WNIR", "solvent injection rate", U_GR, false, [](const WR& w) { return w.solvInj; }},
    {"WNPR", "solvent production rate", U_GR, false, [](const WR& w) { return w.solvProd; }},
    {"WOPRH", "observed oil rate (history)", U_LR, false, [](const WR& w) { return w.oilH >= 0 ? w.oilH : double(NAN); }},
    {"WWPRH", "observed water rate (history)", U_LR, false, [](const WR& w) { return w.waterH >= 0 ? w.waterH : double(NAN); }},
    {"WGPRH", "observed gas rate (history)", U_GR, false, [](const WR& w) { return w.gasH >= 0 ? w.gasH : double(NAN); }},
};
const int kNumWell = int(sizeof(kWell) / sizeof(kWell[0]));

struct CDef { const char* name; const char* desc; int unit; };
const CDef kCompField[] = {
    {"FCPR", "component production rate", U_MOLS}, {"FCIR", "component injection rate", U_MOLS},
    {"FCPT", "cumulative component produced", U_MOL}, {"FCIT", "cumulative component injected", U_MOL},
    {"FCMPR", "component mass production rate", U_KGD}, {"FCMPT", "cumulative component mass produced", U_KG},
};
const CDef kCompWell[] = {
    {"WCPR", "well component production rate", U_MOLS}, {"WCIR", "well component injection rate", U_MOLS},
    {"WCPT", "well cumulative component produced", U_MOL}, {"WCIT", "well cumulative component injected", U_MOL},
};
const char* kObsNames[] = {"WOPR", "WWPR", "WGPR", "WWIR", "WBHP", "WWCT", "WGOR", "FOPR", "FWPR", "FGPR", "FPR", "FWCT", "FGOR"};

std::vector<std::string> splitKey(const std::string& k) {
    std::vector<std::string> f;
    size_t a = 0;
    for (;;) { const size_t b = k.find(':', a); f.push_back(k.substr(a, b == std::string::npos ? std::string::npos : b - a)); if (b == std::string::npos) break; a = b + 1; }
    return f;
}

std::vector<int> wellsFor(const RunResults& R, const std::string& t) {
    std::vector<int> w;
    for (int i = 0; i < int(R.wellNames.size()); ++i) {
        const bool prod = i < int(R.wellTypes.size()) && R.wellTypes[size_t(i)] == WellType::Producer;
        if (t == "*" || (t == "*P" && prod) || (t == "*I" && !prod) || R.wellNames[size_t(i)] == t) w.push_back(i);
    }
    return w;
}

int compIndex(const RunResults& R, const std::string& n) {
    for (int i = 0; i < int(R.components.size()); ++i) if (R.components[size_t(i)].name == n) return i;
    return -1;
}

std::string arrayUnitOf(const std::string& n, const RunResults& R) {
    const Units U = Units::of(R.input.units);
    if (n == "PRESSURE") return U.p;
    if (n == "TEMP") return R.input.units == UnitSystem::Field ? "degF" : "degC";
    if (n == "RS") return U.rs;
    if (n == "SWAT" || n == "SGAS" || n == "SOIL") return "-";
    return "";
}

struct Resolved { std::string label, unit; std::vector<double> x, y; bool markers = false; };

// all curves a key stands for (wildcard well keys give one curve per well)
std::vector<Resolved> resolveKey(const RunResults& R, const std::string& key, double xs, std::string* err = nullptr) {
    std::vector<Resolved> out;
    const auto& S = R.summary;
    const std::vector<std::string> f = splitKey(key);
    auto fail = [&](const std::string& m) { if (err) *err = m; return std::vector<Resolved>{}; };
    if (f.size() < 2) return fail("malformed key");
    const std::string& kind = f[0];
    if (kind == "F") {
        for (const FDef& d : kField) {
            if (f[1] != d.name) continue;
            Resolved r; r.label = d.name; r.unit = unitText(d.unit, R);
            for (size_t k = (d.rate && S.size() > 1) ? 1 : 0; k < S.size(); ++k) { r.x.push_back(S[k].time * xs); r.y.push_back(d.get(S[k])); }
            out.push_back(std::move(r));
            return out;
        }
        return fail("unknown field vector " + f[1]);
    }
    if (kind == "W") {
        if (f.size() < 3) return fail("well missing");
        const WDef* d = nullptr;
        for (const WDef& w : kWell) if (f[1] == w.name) d = &w;
        if (!d) return fail("unknown well vector " + f[1]);
        const std::vector<int> ws = wellsFor(R, f[2]);
        if (ws.empty()) return fail("no well matches " + f[2]);
        for (int w : ws) {
            Resolved r; r.label = std::string(d->name) + " " + R.wellNames[size_t(w)]; r.unit = unitText(d->unit, R);
            double cum = 0;
            for (size_t k = 0; k < S.size(); ++k) {
                const bool has = size_t(w) < S[k].wells.size();
                const double v = has ? d->get(S[k].wells[size_t(w)]) : double(NAN);
                if (d->cum) {
                    if (k > 0 && has && std::isfinite(v)) cum += v * (S[k].time - S[k - 1].time);
                    r.x.push_back(S[k].time * xs); r.y.push_back(cum);
                } else if (k > 0 || S.size() == 1) { r.x.push_back(S[k].time * xs); r.y.push_back(v); }
            }
            out.push_back(std::move(r));
        }
        return out;
    }
    if (kind == "C") {
        if (f.size() < 3) return fail("component missing");
        const int c = compIndex(R, f[2]);
        if (c < 0) return fail("component " + f[2] + " not in this run");
        const double mw = R.components[size_t(c)].mw;
        int which = -1;
        for (int i = 0; i < 6; ++i) if (f[1] == kCompField[i].name) which = i;
        if (which < 0) return fail("unknown component vector " + f[1]);
        Resolved r; r.label = f[1] + " " + f[2]; r.unit = compUnitText(kCompField[which].unit, R);
        const bool rate = which == 0 || which == 1 || which == 4;
        auto at = [&](const std::vector<double>& v) { return size_t(c) < v.size() ? v[size_t(c)] : double(NAN); };
        for (size_t k = (rate && S.size() > 1) ? 1 : 0; k < S.size(); ++k) {
            double v = NAN;
            switch (which) {
                case 0: v = at(S[k].fcompP); break;
                case 1: v = at(S[k].fcompI); break;
                case 2: v = at(S[k].fcompPT); break;
                case 3: v = at(S[k].fcompIT); break;
                case 4: v = at(S[k].fcompP) * mw * 86400.0; break;
                case 5: v = at(S[k].fcompPT) * mw; break;
            }
            r.x.push_back(S[k].time * xs); r.y.push_back(v);
        }
        out.push_back(std::move(r));
        return out;
    }
    if (kind == "WC") {
        if (f.size() < 4) return fail("component / well missing");
        const int c = compIndex(R, f[2]);
        if (c < 0) return fail("component " + f[2] + " not in this run");
        int which = -1;
        for (int i = 0; i < 4; ++i) if (f[1] == kCompWell[i].name) which = i;
        if (which < 0) return fail("unknown well component vector " + f[1]);
        const std::vector<int> ws = wellsFor(R, f[3]);
        if (ws.empty()) return fail("no well matches " + f[3]);
        for (int w : ws) {
            Resolved r; r.label = f[1] + " " + f[2] + " " + R.wellNames[size_t(w)]; r.unit = compUnitText(kCompWell[which].unit, R);
            for (size_t k = (which < 2 && S.size() > 1) ? 1 : 0; k < S.size(); ++k) {
                double v = NAN;
                if (size_t(w) < S[k].wells.size()) {
                    const WellRates& W = S[k].wells[size_t(w)];
                    auto at = [&](const std::vector<double>& a) { return size_t(c) < a.size() ? a[size_t(c)] : double(NAN); };
                    switch (which) {
                        case 0: v = std::max(0.0, at(W.comp)); break;
                        case 1: v = std::max(0.0, -at(W.comp)); break;
                        case 2: v = at(W.compCumP); break;
                        case 3: v = at(W.compCumI); break;
                    }
                }
                r.x.push_back(S[k].time * xs); r.y.push_back(v);
            }
            out.push_back(std::move(r));
        }
        return out;
    }
    if (kind == "G") {
        auto timeOf = [&](size_t k, size_t n) {
            if (n == S.size() && k < S.size()) return S[k].time * xs;
            if (n == R.snaps.size() && k < R.snaps.size()) return R.snaps[k].time * xs;
            return double(k);
        };
        if (f[1] == "COOLED") {
            Resolved r; r.label = "Cooled volume (>10 K)"; r.unit = "%";
            for (size_t k = 0; k < R.cooledFrac.size(); ++k) { r.x.push_back(timeOf(k, R.cooledFrac.size())); r.y.push_back(100.0 * R.cooledFrac[k]); }
            if (r.x.empty()) return fail("no geothermal analysis in this run");
            out.push_back(std::move(r));
            return out;
        }
        if (f.size() < 3) return fail("well missing");
        const bool temp = f[1] == "WTEMP";
        if (!temp && f[1] != "WMW") return fail("unknown geothermal vector " + f[1]);
        const auto& A = temp ? R.wellT : R.wellMW;
        if (A.empty()) return fail("no geothermal analysis in this run");
        for (int w : wellsFor(R, f[2])) {
            Resolved r; r.label = f[1] + " " + R.wellNames[size_t(w)]; r.unit = temp ? unitText(U_T, R) : "MW";
            for (size_t k = 0; k < A.size(); ++k) {
                double v = size_t(w) < A[k].size() ? double(A[k][size_t(w)]) : double(NAN);
                if (temp && v == 0) v = NAN;   // no flow
                r.x.push_back(timeOf(k, A.size())); r.y.push_back(v);
            }
            out.push_back(std::move(r));
        }
        if (out.empty()) return fail("no well matches " + f[2]);
        return out;
    }
    if (kind == "B") {
        if (f.size() < 3) return fail("cell missing");
        int i = 0, j = 0, kk = 0;
        if (std::sscanf(f[2].c_str(), "%d,%d,%d", &i, &j, &kk) != 3) return fail("cell must be i,j,k");
        const Grid* g = R.grid.get();
        if (!g) return fail("no grid");
        if (i < 1 || j < 1 || kk < 1 || i > g->nx || j > g->ny || kk > g->nz) return fail("cell outside the grid");
        const size_t gi = size_t(i - 1) + size_t(g->nx) * (size_t(j - 1) + size_t(g->ny) * size_t(kk - 1));
        const int a = gi < g->activeOf.size() ? g->activeOf[gi] : -1;
        if (a < 0) return fail("cell is inactive");
        Resolved r; r.label = f[1] + " (" + f[2] + ")"; r.unit = arrayUnitOf(f[1], R); r.markers = R.snaps.size() < 40;
        for (const Snapshot& sn : R.snaps) {
            const std::vector<float>& v = sn.field(f[1]);
            if (size_t(a) < v.size()) { r.x.push_back(sn.time * xs); r.y.push_back(v[size_t(a)]); }
        }
        if (r.x.empty()) return fail("array " + f[1] + " not in the report steps");
        out.push_back(std::move(r));
        return out;
    }
    if (kind == "O") {
        int q = -1;
        for (int i = 0; i < 13; ++i) if (f[1] == kObsNames[i]) q = i;
        if (q < 0) return fail("unknown observed quantity " + f[1]);
        const std::string well = f.size() > 2 ? f[2] : std::string();
        Resolved r; r.label = f[1] + (well.empty() ? "" : " " + well) + " (observed)"; r.markers = true;
        std::vector<std::pair<double, double>> pts;
        for (const Observation& o : R.input.observations)
            if (int(o.qty) == q && o.well == well) pts.push_back({o.time * xs, o.value});
        std::sort(pts.begin(), pts.end());
        for (const auto& p : pts) { r.x.push_back(p.first); r.y.push_back(p.second); }
        if (r.x.empty()) return fail("no observations for " + key);
        out.push_back(std::move(r));
        return out;
    }
    return fail("unknown key type " + kind);
}

bool anyNonZero(const RunResults& R, const FDef& d) {
    if (d.always) return true;
    for (const SP& s : R.summary) { const double v = d.get(s); if (std::isfinite(v) && v != 0) return true; }
    return false;
}
bool anyNonZeroWell(const RunResults& R, const WDef& d) {
    for (const SP& s : R.summary) for (const WR& w : s.wells) { const double v = d.get(w); if (std::isfinite(v) && v != 0) return true; }
    return false;
}

bool matches(const char* filter, const char* a, const char* b) {
    if (!filter[0]) return true;
    auto has = [&](const char* s) {
        std::string h = s, n = filter;
        std::transform(h.begin(), h.end(), h.begin(), ::tolower);
        std::transform(n.begin(), n.end(), n.begin(), ::tolower);
        return h.find(n) != std::string::npos;
    };
    return has(a) || has(b);
}

std::string defaultPathFrom(const std::string& settingsPath) {
    if (settingsPath.empty()) return "custom_plots.txt";
    return (fs::u8path(settingsPath).parent_path() / "custom_plots.txt").u8string();
}
}  // namespace

// ============================================================================ App
void App::loadCustomPlots() {
    cplotsPath_ = defaultPathFrom(settingsPath_);
    if (readCustomPlots(cplotsPath_, cplots_)) return;
    // first start: a small example library (edit / delete freely)
    cplots_.clear();
    CustomPlot a; a.name = "Production overview";
    a.curves = {{"F:FOPR", "", 0, 0, 1}, {"F:FWPR", "", 0, 0, 1}, {"F:FWIR", "", 0, 1, 1}, {"F:FGPR", "", 1, 0, 1}};
    CustomPlot b; b.name = "Pressure and BHP";
    b.curves = {{"F:FPR", "", 0, 0, 1}, {"W:WBHP:*", "", 0, 1, 1}};
    CustomPlot c; c.name = "Producer water cut"; c.ownTab = false;
    c.curves = {{"W:WWCT:*P", "", 0, 0, 1}};
    cplots_ = {a, b, c};
    saveCustomPlots();
}

void App::saveCustomPlots() {
    if (cplotsPath_.empty()) cplotsPath_ = defaultPathFrom(settingsPath_);
    if (!writeCustomPlots(cplotsPath_, cplots_)) notify("Could not save the custom plots to " + cplotsPath_);
}

void App::drawCustomPlot(const CustomPlot& P, const RunResults& R, double xs, const char* xl, const std::function<void()>& cursor,
                         ImVec2 size, int id) {
    struct Item { Resolved r; const CustomCurve* c; };
    std::vector<Item> items;
    std::string firstErr;
    for (const CustomCurve& c : P.curves) {
        std::string err;
        for (Resolved& r : resolveKey(R, c.key, xs, &err)) {
            if (!c.label.empty()) r.label = c.label + (r.label.find(' ') != std::string::npos && c.key.find('*') != std::string::npos ? r.label.substr(r.label.find(' ')) : "");
            if (c.scale != 1) for (double& v : r.y) v *= c.scale;
            items.push_back({std::move(r), &c});
        }
        if (!err.empty() && firstErr.empty()) firstErr = c.key + ": " + err;
    }
    bool y2 = false;
    std::string yl = P.ylabel, y2l = P.y2label;
    for (const Item& it : items) {
        if (it.c->axis == 1) { y2 = true; if (P.y2label.empty() && y2l.empty()) y2l = it.r.unit; }
        else if (P.ylabel.empty() && yl.empty()) yl = it.r.unit;
    }
    char pid[300];
    std::snprintf(pid, sizeof(pid), "%s##cplot%d", P.name.c_str(), id);
    if (!ImPlot::BeginPlot(pid, size)) return;
    recPanel(P.name.c_str(), xl, yl.c_str(), y2 ? y2l.c_str() : "");
    const ImPlotAxisFlags af = ImPlotAxisFlags_AutoFit;
    ImPlot::SetupAxes(xl, yl.c_str(), af, af);
    if (P.logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10); else setupSI(ImAxis_Y1);
    if (y2) {
        ImPlot::SetupAxis(ImAxis_Y2, y2l.c_str(), af | ImPlotAxisFlags_AuxDefault);
        if (P.logY) ImPlot::SetupAxisScale(ImAxis_Y2, ImPlotScale_Log10); else setupSI(ImAxis_Y2);
    }
    if (items.empty()) {
        ImPlot::PlotText(firstErr.empty() ? "No curves - press Edit and add curves from the catalog" : firstErr.c_str(), 0.5, 0.5);
    }
    int ci = 0;
    for (const Item& it : items) {
        const ImU32 col = seriesColor(ci++);
        useY2(it.c->axis == 1);
        const int style = it.r.markers && it.c->style == 0 ? 2 : it.c->style;
        ImPlotSpec sp;
        sp.LineColor = col4(col);
        sp.LineWeight = 2.0f * g_plotLineScale;
        if (style == 1) sp.LineWeight = 1.0f * g_plotLineScale;   // thin
        if (style >= 2) { sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 3.5f; sp.MarkerFillColor = col4(col); sp.MarkerLineColor = col4(col); }
        // finite (and, on a log axis, positive) runs only: ImPlot gets no NaN; gaps split the curve (one legend entry)
        std::vector<double> X, Y, RX, RY;
        auto flush = [&]() {
            if (X.empty()) return;
            if (style == 2) ImPlot::PlotScatter(it.r.label.c_str(), X.data(), Y.data(), int(X.size()), sp);
            else ImPlot::PlotLine(it.r.label.c_str(), X.data(), Y.data(), int(X.size()), sp);
            X.clear(); Y.clear();
        };
        for (size_t i = 0; i < it.r.x.size(); ++i) {
            const double v = it.r.y[i];
            if (!std::isfinite(v) || (P.logY && v <= 0)) { flush(); continue; }
            X.push_back(it.r.x[i]); Y.push_back(v); RX.push_back(it.r.x[i]); RY.push_back(v);
        }
        flush();
        if (!RX.empty()) recSeries(it.r.label.c_str(), RX.data(), RY.data(), int(RX.size()), col, 2.0f, style == 2, style == 3);
    }
    useY2(false);
    cursor();
    recEnd();
    ImPlot::EndPlot();
}

void App::customPlotsTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor) {
    (void)tSnap;
    static bool editing = false, grid = false;
    static char filter[64] = "";
    static int wellPick = 0, compPick = 0, arrPick = 0, cell[3] = {1, 1, 1};
    const float fs = ImGui::GetFontSize();
    if (cplots_.empty()) cplotSel_ = -1;
    cplotSel_ = std::clamp(cplotSel_, cplots_.empty() ? -1 : 0, int(cplots_.size()) - 1);
    bool changed = false;

    // ---------------------------------------------------------------- library list
    ImGui::BeginChild("cplib", ImVec2(fs * 15, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
    ImGui::TextDisabled("My plots (saved in your profile)");
    for (int i = 0; i < int(cplots_.size()); ++i) {
        ImGui::PushID(i);
        char lab[300]; std::snprintf(lab, sizeof(lab), "%s%s", cplots_[size_t(i)].ownTab ? ICON_STAR " " : "   ", cplots_[size_t(i)].name.c_str());
        if (ImGui::Selectable(lab, cplotSel_ == i && !grid)) { cplotSel_ = i; grid = false; }
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%d curve(s)%s", int(cplots_[size_t(i)].curves.size()), cplots_[size_t(i)].ownTab ? ", shown as its own tab" : "");
        ImGui::PopID();
    }
    ImGui::Separator();
    if (ImGui::Button("New")) {
        CustomPlot p; p.name = "Plot " + std::to_string(cplots_.size() + 1);
        cplots_.push_back(p); cplotSel_ = int(cplots_.size()) - 1; editing = true; grid = false; changed = true;
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(cplotSel_ < 0);
    if (ImGui::Button("Duplicate")) { CustomPlot p = cplots_[size_t(cplotSel_)]; p.name += " (copy)"; cplots_.insert(cplots_.begin() + cplotSel_ + 1, p); ++cplotSel_; changed = true; }
    ImGui::SameLine();
    if (ImGui::Button("Delete")) { cplots_.erase(cplots_.begin() + cplotSel_); cplotSel_ = std::min(cplotSel_, int(cplots_.size()) - 1); changed = true; }
    if (cplotSel_ >= 0) {
        if (ImGui::ArrowButton("##up", ImGuiDir_Up) && cplotSel_ > 0) { std::swap(cplots_[size_t(cplotSel_)], cplots_[size_t(cplotSel_ - 1)]); --cplotSel_; changed = true; }
        ImGui::SameLine();
        if (ImGui::ArrowButton("##dn", ImGuiDir_Down) && cplotSel_ + 1 < int(cplots_.size())) { std::swap(cplots_[size_t(cplotSel_)], cplots_[size_t(cplotSel_ + 1)]); ++cplotSel_; changed = true; }
        ImGui::SameLine();
    }
    ImGui::EndDisabled();
    if (ImGui::Checkbox("Edit", &editing) && editing) grid = false;
    if (ImGui::Checkbox("Show all plots together", &grid) && grid) editing = false;
    ImGui::Separator();
    if (ImGui::Button("Import...")) {
        std::string p;
        if (fileDialog(false, "ResSim plot library (*.rsplots)\0*.rsplots\0All files (*.*)\0*.*\0", "rsplots", p)) {
            std::vector<CustomPlot> in;
            if (readCustomPlots(p, in) && !in.empty()) { cplots_.insert(cplots_.end(), in.begin(), in.end()); changed = true; notify("Imported " + std::to_string(in.size()) + " plot(s)"); }
            else notify("No plots found in " + p);
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Export...")) {
        std::string p = "my_plots.rsplots";
        if (fileDialog(true, "ResSim plot library (*.rsplots)\0*.rsplots\0All files (*.*)\0*.*\0", "rsplots", p)) {
            if (writeCustomPlots(p, cplots_)) notify("Saved " + std::to_string(cplots_.size()) + " plot(s) to " + p);
        }
    }
    ImGui::TextWrapped(ICON_STAR " = own tab in the Plots window. Right-click any plot and choose \"Export this plot\" to export just that plot.");
    ImGui::EndChild();
    ImGui::SameLine();

    // ---------------------------------------------------------------- plot area
    ImGui::BeginChild("cpmain", ImVec2(0, 0));
    if (grid && !cplots_.empty()) {
        const ImVec2 av = ImGui::GetContentRegionAvail();
        const PlotGrid G = fitGrid(int(cplots_.size()), av, fs * 22, fs * 12, 1.6f);
        fig::recorder().beginFigure("run.custom", G.rows, G.cols);
        if (G.scroll) ImGui::BeginChild("##cpscroll", av);
        const ImVec2 cellSz((G.size.x - ImGui::GetStyle().ItemSpacing.x * float(G.cols - 1)) / float(G.cols),
                            (G.size.y - ImGui::GetStyle().ItemSpacing.y * float(G.rows - 1)) / float(G.rows));
        for (int i = 0; i < int(cplots_.size()); ++i) {
            if (i % G.cols) ImGui::SameLine();
            drawCustomPlot(cplots_[size_t(i)], R, xs, xl, cursor, cellSz, i);
        }
        if (G.scroll) ImGui::EndChild();
        fig::recorder().endFigure();
    } else if (cplotSel_ >= 0) {
        CustomPlot& P = cplots_[size_t(cplotSel_)];
        if (editing) {
            ImGui::BeginChild("cpedit", ImVec2(0, fs * 15), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeY);
            // left: plot settings + curves
            ImGui::BeginChild("cpcurves", ImVec2(ImGui::GetContentRegionAvail().x * 0.55f, 0), ImGuiChildFlags_ResizeX);
            char buf[256];
            std::snprintf(buf, sizeof(buf), "%s", P.name.c_str());
            ImGui::SetNextItemWidth(fs * 14);
            if (ImGui::InputText("Name", buf, sizeof(buf))) { P.name = buf; changed = true; }
            ImGui::SameLine();
            if (ImGui::Checkbox("Own tab", &P.ownTab)) changed = true;
            ImGui::SameLine();
            if (ImGui::Checkbox("Log Y", &P.logY)) changed = true;
            std::snprintf(buf, sizeof(buf), "%s", P.ylabel.c_str());
            ImGui::SetNextItemWidth(fs * 9);
            if (ImGui::InputTextWithHint("##yl", "left axis label (auto)", buf, sizeof(buf))) { P.ylabel = buf; changed = true; }
            ImGui::SameLine();
            std::snprintf(buf, sizeof(buf), "%s", P.y2label.c_str());
            ImGui::SetNextItemWidth(fs * 9);
            if (ImGui::InputTextWithHint("##y2l", "right axis label (auto)", buf, sizeof(buf))) { P.y2label = buf; changed = true; }
            if (ImGui::BeginTable("cpc", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit, ImVec2(0, -FLT_MIN))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Curve key", ImGuiTableColumnFlags_WidthStretch);
                ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthFixed, fs * 7);
                ImGui::TableSetupColumn("Axis", ImGuiTableColumnFlags_WidthFixed, fs * 4.5f);
                ImGui::TableSetupColumn("Style", ImGuiTableColumnFlags_WidthFixed, fs * 6.5f);
                ImGui::TableSetupColumn("Scale", ImGuiTableColumnFlags_WidthFixed, fs * 4.5f);
                ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 1.6f);
                ImGui::TableHeadersRow();
                int del = -1;
                for (int i = 0; i < int(P.curves.size()); ++i) {
                    CustomCurve& c = P.curves[size_t(i)];
                    ImGui::PushID(i);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    std::snprintf(buf, sizeof(buf), "%s", c.key.c_str());
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::InputText("##k", buf, sizeof(buf))) { c.key = buf; changed = true; }
                    std::string err;
                    if (resolveKey(R, c.key, xs, &err).empty() && ImGui::IsItemHovered()) ImGui::SetTooltip("%s", err.c_str());
                    ImGui::TableNextColumn();
                    std::snprintf(buf, sizeof(buf), "%s", c.label.c_str());
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::InputTextWithHint("##l", "auto", buf, sizeof(buf))) { c.label = buf; changed = true; }
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::Combo("##a", &c.axis, "left\0right\0")) changed = true;
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::Combo("##s", &c.style, "line\0thin line\0markers\0line+markers\0")) changed = true;
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::InputDouble("##sc", &c.scale, 0, 0, "%.4g")) changed = true;
                    ImGui::TableNextColumn();
                    if (ImGui::SmallButton("x")) del = i;
                    ImGui::PopID();
                }
                if (del >= 0) { P.curves.erase(P.curves.begin() + del); changed = true; }
                ImGui::EndTable();
            }
            ImGui::EndChild();
            ImGui::SameLine();
            // right: catalog
            ImGui::BeginChild("cpcat", ImVec2(0, 0));
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputTextWithHint("##flt", "filter vectors (e.g. rate, CO2, BHP)", filter, sizeof(filter));
            std::vector<std::string> wellOpts = {"*", "*P", "*I"};
            for (const auto& n : R.wellNames) wellOpts.push_back(n);
            wellPick = std::clamp(wellPick, 0, int(wellOpts.size()) - 1);
            auto wellLabel = [&](int i) { return i == 0 ? std::string("all wells") : i == 1 ? std::string("all producers") : i == 2 ? std::string("all injectors") : wellOpts[size_t(i)]; };
            auto add = [&](const std::string& key) { P.curves.push_back({key, "", 0, 0, 1}); changed = true; };
            ImGui::BeginChild("cpcatl", ImVec2(0, 0));
            if (ImGui::TreeNodeEx("Field", ImGuiTreeNodeFlags_DefaultOpen)) {
                for (const FDef& d : kField) {
                    if (!matches(filter, d.name, d.desc) || !anyNonZero(R, d)) continue;
                    char lab[160]; std::snprintf(lab, sizeof(lab), "%-7s %s", d.name, d.desc);
                    if (ImGui::Selectable(lab)) add(std::string("F:") + d.name);
                }
                ImGui::TreePop();
            }
            if (!R.wellNames.empty() && ImGui::TreeNodeEx("Wells", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::SetNextItemWidth(fs * 10);
                if (ImGui::BeginCombo("Well(s)", wellLabel(wellPick).c_str())) {
                    for (int i = 0; i < int(wellOpts.size()); ++i) if (ImGui::Selectable(wellLabel(i).c_str(), i == wellPick)) wellPick = i;
                    ImGui::EndCombo();
                }
                for (const WDef& d : kWell) {
                    if (!matches(filter, d.name, d.desc) || !anyNonZeroWell(R, d)) continue;
                    char lab[160]; std::snprintf(lab, sizeof(lab), "%-7s %s", d.name, d.desc);
                    if (ImGui::Selectable(lab)) add(std::string("W:") + d.name + ":" + wellOpts[size_t(wellPick)]);
                }
                ImGui::TreePop();
            }
            if (!R.components.empty() && ImGui::TreeNodeEx("Components", ImGuiTreeNodeFlags_DefaultOpen)) {
                compPick = std::clamp(compPick, 0, int(R.components.size()) - 1);
                ImGui::SetNextItemWidth(fs * 10);
                if (ImGui::BeginCombo("Component", R.components[size_t(compPick)].name.c_str())) {
                    for (int i = 0; i < int(R.components.size()); ++i) if (ImGui::Selectable(R.components[size_t(i)].name.c_str(), i == compPick)) compPick = i;
                    ImGui::EndCombo();
                }
                const std::string cn = R.components[size_t(compPick)].name;
                for (const CDef& d : kCompField) {
                    if (!matches(filter, d.name, d.desc)) continue;
                    char lab[160]; std::snprintf(lab, sizeof(lab), "%-7s %s (field)", d.name, d.desc);
                    if (ImGui::Selectable(lab)) add(std::string("C:") + d.name + ":" + cn);
                }
                for (const CDef& d : kCompWell) {
                    if (!matches(filter, d.name, d.desc) || R.wellNames.empty()) continue;
                    char lab[200]; std::snprintf(lab, sizeof(lab), "%-7s %s (%s)", d.name, d.desc, wellLabel(wellPick).c_str());
                    if (ImGui::Selectable(lab)) add(std::string("WC:") + d.name + ":" + cn + ":" + wellOpts[size_t(wellPick)]);
                }
                ImGui::TreePop();
            }
            if (!R.wellT.empty() && ImGui::TreeNodeEx("Geothermal", ImGuiTreeNodeFlags_DefaultOpen)) {
                char lab[200];
                std::snprintf(lab, sizeof(lab), "WTEMP   well-stream temperature (%s)", wellLabel(wellPick).c_str());
                if (matches(filter, "WTEMP", "temperature") && ImGui::Selectable(lab)) add("G:WTEMP:" + wellOpts[size_t(wellPick)]);
                std::snprintf(lab, sizeof(lab), "WMW     enthalpy rate (%s)", wellLabel(wellPick).c_str());
                if (matches(filter, "WMW", "enthalpy power") && ImGui::Selectable(lab)) add("G:WMW:" + wellOpts[size_t(wellPick)]);
                if (matches(filter, "COOLED", "cooled volume") && ImGui::Selectable("COOLED  bulk volume cooled > 10 K")) add("G:COOLED");
                ImGui::TreePop();
            }
            if (!R.snaps.empty() && R.grid && ImGui::TreeNodeEx("Cell (block) values", ImGuiTreeNodeFlags_DefaultOpen)) {
                std::vector<std::string> arrs;
                const Snapshot& s0 = R.snaps.back();
                for (const char* n : {"PRESSURE", "SWAT", "SOIL", "SGAS", "RS", "RV", "POLYMER", "PADS", "XCO2"}) if (!s0.field(n).empty()) arrs.push_back(n);
                for (const auto& e : s0.extra) if (!e.second.empty()) arrs.push_back(e.first);
                if (!arrs.empty()) {
                    arrPick = std::clamp(arrPick, 0, int(arrs.size()) - 1);
                    ImGui::SetNextItemWidth(fs * 8);
                    if (ImGui::BeginCombo("Array", arrs[size_t(arrPick)].c_str())) {
                        for (int i = 0; i < int(arrs.size()); ++i) if (ImGui::Selectable(arrs[size_t(i)].c_str(), i == arrPick)) arrPick = i;
                        ImGui::EndCombo();
                    }
                    ImGui::SetNextItemWidth(fs * 8);
                    ImGui::InputInt3("i, j, k", cell);
                    ImGui::SameLine();
                    if (ImGui::SmallButton("Add cell curve")) {
                        char key[160]; std::snprintf(key, sizeof(key), "B:%s:%d,%d,%d", arrs[size_t(arrPick)].c_str(), cell[0], cell[1], cell[2]);
                        std::string err;
                        if (resolveKey(R, key, xs, &err).empty()) notify(std::string(key) + ": " + err); else add(key);
                    }
                }
                ImGui::TreePop();
            }
            if (!R.input.observations.empty() && ImGui::TreeNodeEx("Observed data", ImGuiTreeNodeFlags_DefaultOpen)) {
                std::vector<std::pair<int, std::string>> seen;
                for (const Observation& o : R.input.observations) {
                    const std::pair<int, std::string> k{int(o.qty), o.well};
                    if (std::find(seen.begin(), seen.end(), k) != seen.end()) continue;
                    seen.push_back(k);
                }
                for (const auto& k : seen) {
                    char lab[160]; std::snprintf(lab, sizeof(lab), "%s %s (observed)", kObsNames[k.first], k.second.c_str());
                    if (!matches(filter, lab, "observed history")) continue;
                    if (ImGui::Selectable(lab)) add(std::string("O:") + kObsNames[k.first] + (k.second.empty() ? "" : ":" + k.second));
                }
                ImGui::TreePop();
            }
            ImGui::EndChild();
            ImGui::EndChild();
            ImGui::EndChild();
        }
        fig::recorder().beginFigure("run.custom", 1, 1);
        drawCustomPlot(P, R, xs, xl, cursor, ImVec2(-1, -1), cplotSel_);
        fig::recorder().endFigure();
    } else {
        ui::hint("No custom plots yet. Press New, then pick curves from the catalog.");
    }
    ImGui::EndChild();
    if (changed) saveCustomPlots();
}
