// C35 EOR section, part 2: chemical flooding with the own engine - surfactant (SURFACT), alkaline (ALKALINE) and foam
// (FOAM). The ECLIPSE keywords are kept verbatim in Case::extraKeywords (that is what the engine reads and what the deck
// writer exports), so this editor parses / re-writes those keyword bodies. A keyword is only re-written when it is
// edited here, so imported decks keep their layout and comments otherwise.
#include "app.h"
#include "widgets.h"
#include "implot.h"
#include "ressim/ecl.h"
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <sstream>
#include <string>

using namespace ressim;

namespace {
using Row2 = std::array<double, 2>;
using Row3 = std::array<double, 3>;

// ---- verbatim keyword access
ExtraKeyword* findKw(Case& c, const std::string& n) {
    for (auto& k : c.extraKeywords) if (k.name == n) return &k;
    return nullptr;
}
void setKw(Case& c, const std::string& n, DeckSection s, const std::string& body) {
    if (ExtraKeyword* k = findKw(c, n)) { k->body = body; return; }
    ExtraKeyword k; k.name = n; k.section = s; k.body = body; k.note = "chemical EOR (C35 editor)";
    c.extraKeywords.push_back(k);
}
void removeKw(Case& c, const std::string& n) {
    c.extraKeywords.erase(std::remove_if(c.extraKeywords.begin(), c.extraKeywords.end(), [&](const ExtraKeyword& k) { return k.name == n; }),
                          c.extraKeywords.end());
}
// records of a keyword body ('/'-terminated, '--' comments removed); each record = whitespace tokens with N*v expanded
std::vector<std::vector<std::string>> records(const std::string& body) {
    std::vector<std::vector<std::string>> out;
    std::vector<std::string> cur;
    std::istringstream is(body);
    std::string line;
    while (std::getline(is, line)) {
        const size_t cm = line.find("--");
        if (cm != std::string::npos) line.erase(cm);
        std::string tok;
        bool q = false;
        auto flush = [&] {
            if (tok.empty()) return;
            const size_t st = tok.find('*');
            if (st != std::string::npos && st > 0 && tok[0] != '\'') {
                const int n = std::atoi(tok.substr(0, st).c_str());
                const std::string v = tok.substr(st + 1);
                for (int k = 0; k < n; ++k) cur.push_back(v.empty() ? "*" : v);
            } else cur.push_back(tok);
            tok.clear();
        };
        for (char ch : line) {
            if (ch == '\'') { q = !q; tok += ch; continue; }
            if (!q && ch == '/') { flush(); out.push_back(cur); cur.clear(); continue; }
            if (!q && std::isspace(static_cast<unsigned char>(ch))) { flush(); continue; }
            tok += ch;
        }
        flush();
    }
    return out;
}
bool num(const std::string& t, double& v) {
    char* e = nullptr; v = std::strtod(t.c_str(), &e);
    return e && e != t.c_str() && *e == 0 && std::isfinite(v);
}
// tables of ncol columns; an empty record copies the previous table (ECLIPSE); false = not purely numeric
bool parseTables(const std::string& body, int ncol, std::vector<std::vector<double>>& tabs) {
    tabs.clear();
    for (const auto& r : records(body)) {
        std::vector<double> t;
        if (r.empty()) { if (tabs.empty()) return false; tabs.push_back(tabs.back()); continue; }
        for (const auto& s : r) { double v; if (!num(s, v)) return false; t.push_back(v); }
        if (int(t.size()) % ncol) return false;
        tabs.push_back(t);
    }
    return !tabs.empty();
}
std::string writeTables(const std::vector<std::vector<double>>& tabs, int ncol) {
    std::string b;
    char buf[64];
    for (const auto& t : tabs) {
        for (size_t i = 0; i < t.size(); ++i) {
            std::snprintf(buf, sizeof buf, "%s%.6g", i % size_t(ncol) == 0 ? "  " : "   ", t[i]);
            b += buf;
            if (i % size_t(ncol) == size_t(ncol - 1)) b += i + 1 == t.size() ? " /\n" : "\n";
        }
        if (t.empty()) b += " /\n";
    }
    return b;
}
// unit helpers: seeds are FIELD values from the ECLIPSE ASP / FOAM examples
struct ChemUnits { bool field = true; double conc = 1, st = 1, p = 1, vel = 1; };
ChemUnits unitsOf(UnitSystem u) {
    ChemUnits r; r.field = u == UnitSystem::Field;
    if (!r.field) { r.conc = 2.853010; r.st = 175.1268; r.p = 0.0689476; r.vel = 0.3048; }
    return r;
}
// one table keyword description
struct TabKw {
    const char* name; const char* x; const char* y; bool perPvt;   // per PVTNUM (else per SATNUM)
    const char* tip;
};
const TabKw kSurf[] = {
    {"SURFVISC", "c surf", "water viscosity (cP)", true, "Water viscosity with dissolved surfactant (the water mobility is multiplied by SURFVISC(0) / SURFVISC(c))."},
    {"SURFST", "c surf", "oil-water IFT", true, "Oil-water interfacial tension vs surfactant concentration (field lb/in, metric N/m). Pcow is scaled by IFT(c) / IFT(0)."},
    {"SURFCAPD", "log10 Nc", "miscibility F (0-1)", false, "Miscibility function of the capillary number: 0 = immiscible (SATNUM curves), 1 = miscible (SURFNUM curves)."},
    {"SURFADS", "c surf", "adsorbed (mass/mass rock)", false, "Surfactant adsorption isotherm (reversible or irreversible from SURFROCK)."},
};
const TabKw kAlk[] = {
    {"ALSURFST", "c alk", "IFT multiplier", true, "Multiplier on the surfactant IFT (SURFST) vs alkaline concentration."},
    {"ALSURFAD", "c alk", "surfactant ads. multiplier", false, "Multiplier on the surfactant adsorption vs alkaline concentration."},
    {"ALPOLADS", "c alk", "polymer ads. multiplier", false, "Multiplier on the polymer adsorption (PLYADS) vs alkaline concentration."},
    {"ALKADS", "c alk", "adsorbed (mass/mass rock)", false, "Alkaline adsorption isotherm (needs ALKROCK)."},
};
const TabKw kFoam[] = {
    {"FOAMMOB", "c foam", "gas mobility reduction", true, "Gas mobility reduction factor vs foam concentration (1 = none)."},
    {"FOAMMOBP", "oil pressure", "weakening", true, "Pressure weakening of the foam (functional reduction 1 - (1-Mc)(1-Mp)(1-Ms))."},
    {"FOAMMOBS", "gas velocity", "weakening", true, "Shear weakening: foam effect lost at high gas Darcy velocity."},
    {"FOAMADS", "c foam", "adsorbed (mass/mass rock)", false, "Foam adsorption isotherm (needs FOAMROCK)."},
    {"FOAMDCYW", "Sw", "half-life (days)", false, "Foam half-life vs water saturation (first-order decay)."},
    {"FOAMDCYO", "So", "half-life (days)", false, "Foam half-life vs oil saturation; the decay rates of FOAMDCYW and FOAMDCYO add up."},
};
std::vector<double> seedTable(const std::string& n, const ChemUnits& u) {
    const double c = u.conc;
    if (n == "SURFVISC") return {0, 0.5, 30 * c, 5.0};
    if (n == "SURFST") return {0, 1e-4 * u.st, 1 * c, 1e-8 * u.st, 30 * c, 1e-8 * u.st};
    if (n == "SURFCAPD") return {-9, 0, -4.5, 0, -2, 1, 10, 1};
    if (n == "SURFADS") return {0, 0, 1 * c, 0.0005, 30 * c, 0.0005};
    if (n == "ALSURFST" || n == "ALPOLADS") return {0, 1, 0.1 * c, 0.9, 0.2 * c, 0.8, 0.5 * c, 0.7};
    if (n == "ALSURFAD") return {0, 1, 0.1 * c, 0.99995, 0.2 * c, 0.995, 0.5 * c, 0.5};
    if (n == "ALKADS") return {0, 0, 1 * c, 0.0001, 50 * c, 0.0001};
    if (n == "FOAMMOB") return {0, 1, 0.001 * c, 0.4, 0.1 * c, 0.1, 1.2 * c, 0.05};
    if (n == "FOAMMOBP") return {3000 * u.p, 0, 6000 * u.p, 0.2};
    if (n == "FOAMMOBS") return {0, 0, 4 * u.vel, 0.1};
    if (n == "FOAMADS") return {0, 0, 1 * c, 0.00005, 30 * c, 0.00005};
    if (n == "FOAMDCYW") return {0, 3000, 1, 2000};
    if (n == "FOAMDCYO") return {0, 3000, 1, 2500};
    return {};
}
void seedKw(Case& c, const std::string& n, int ntab, const ChemUnits& u) {
    if (findKw(c, n)) return;
    std::vector<std::vector<double>> t(size_t(std::max(1, ntab)), seedTable(n, u));
    setKw(c, n, DeckSection::Props, writeTables(t, 2));
}
// per-region scalar records (SURFROCK / FOAMROCK: adsorption index, rock density; ALKROCK: max adsorption... ) as text
void seedRock(Case& c, const char* n, int ntab, const ChemUnits& u) {
    if (findKw(c, n)) return;
    std::string b;
    char buf[64];
    for (int k = 0; k < std::max(1, ntab); ++k) {
        if (std::string(n) == "ALKROCK") std::snprintf(buf, sizeof buf, "  0.0001 /\n");
        else if (std::string(n) == "FOAMROCK")   // FOAMROCK (as PLYROCK) in lb/ft3 (field), SURFROCK in lb/rb - as OPM Flow
            std::snprintf(buf, sizeof buf, "  1  %.6g /\n", u.field ? 2650 * 0.0624280 : 2650.0);
        else std::snprintf(buf, sizeof buf, "  1  %.6g /\n", u.field ? 2650 * 0.350507 : 2650.0);   // kg/rm3 or lb/rb
        b += buf;
    }
    setKw(c, n, DeckSection::Props, b);
}
// table editor of one keyword (2 columns, one table per region)
bool tableEditor(Case& c, const TabKw& K, int ntab, int& reg) {
    ExtraKeyword* k = findKw(c, K.name);
    ImGui::PushID(K.name);
    bool ch = false;
    ImGui::TextUnformatted(K.name);
    ImGui::SameLine(); ui::hint(K.tip);
    if (!k) {
        ImGui::TextDisabled("not given");
        ImGui::SameLine();
        if (ImGui::SmallButton("Add with typical data")) { seedKw(c, K.name, ntab, unitsOf(c.units)); ch = true; }
        ImGui::PopID();
        return ch;
    }
    std::vector<std::vector<double>> tabs;
    if (!parseTables(k->body, 2, tabs)) {
        ImGui::TextDisabled("non-numeric records - edit the text");
        static std::map<std::string, std::string> buf;
        std::string& s = buf[K.name];
        if (s.empty()) s = k->body;
        char tmp[4096]; std::snprintf(tmp, sizeof tmp, "%s", s.c_str());
        if (ImGui::InputTextMultiline("##raw", tmp, sizeof tmp, ImVec2(-1, ImGui::GetFontSize() * 6))) { s = tmp; k->body = s; ch = true; }
        ImGui::PopID();
        return ch;
    }
    while (int(tabs.size()) < ntab) tabs.push_back(tabs.back());
    if (ntab > 1) {
        reg = std::clamp(reg, 0, ntab - 1);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::SliderInt("##reg", &reg, 0, ntab - 1, K.perPvt ? "PVTNUM %d" : "SATNUM %d");
        ImGui::SameLine();
    } else reg = 0;
    if (ImGui::SmallButton("Remove")) { removeKw(c, K.name); ImGui::PopID(); return true; }
    std::vector<double>& t = tabs[size_t(std::min(reg, int(tabs.size()) - 1))];
    std::vector<Row2> rows;
    for (size_t i = 0; i + 1 < t.size(); i += 2) rows.push_back({t[i], t[i + 1]});
    const float h = ImGui::GetFontSize() * 7.5f;
    if (ImGui::BeginTable("tp", 2, ImGuiTableFlags_SizingStretchSame)) {
        ImGui::TableNextColumn();
        if (ui::editTable<Row2>("t", rows, {K.x, K.y}, [](Row2& r, int col) -> double* { return &r[size_t(col)]; }, h)) {
            std::sort(rows.begin(), rows.end(), [](const Row2& a, const Row2& b) { return a[0] < b[0]; });
            t.clear();
            for (auto& r : rows) { t.push_back(r[0]); t.push_back(r[1]); }
            k->body = writeTables(tabs, 2);
            ch = true;
        }
        ImGui::TableNextColumn();
        std::vector<double> x, y;
        for (auto& r : rows) { x.push_back(r[0]); y.push_back(r[1]); }
        if (ImPlot::BeginPlot("##plot", ImVec2(-1, h), ImPlotFlags_NoLegend | ImPlotFlags_NoMenus)) {
            ImPlot::SetupAxes(K.x, K.y, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlotSpec s; s.LineColor = col4(Pal::Accent); s.LineWeight = 2; s.Marker = ImPlotMarker_Circle;
            if (!x.empty()) ImPlot::PlotLine(K.name, x.data(), y.data(), int(x.size()), s);
            ImPlot::EndPlot();
        }
        ImGui::EndTable();
    }
    ImGui::PopID();
    return ch;
}
// SURFROCK / FOAMROCK (adsorption index 1 = desorption allowed, 2 = irreversible; rock density) and ALKROCK
bool rockEditor(Case& c, const char* n, int ntab) {
    ExtraKeyword* k = findKw(c, n);
    ImGui::PushID(n);
    bool ch = false;
    const bool alk = std::string(n) == "ALKROCK";
    ImGui::TextUnformatted(n);
    ImGui::SameLine();
    ui::hint(alk ? "Maximum alkaline adsorption per SATNUM region."
                 : std::string(n) == "FOAMROCK" ? "Per SATNUM region: adsorption index (1 = reversible, 2 = irreversible) and rock mass density (field lb/ft3, metric kg/rm3)."
                                                : "Per SATNUM region: adsorption index (1 = reversible, 2 = irreversible) and rock mass density (field lb/rb, metric kg/rm3).");
    if (!k) {
        ImGui::TextDisabled("not given"); ImGui::SameLine();
        if (ImGui::SmallButton("Add")) { seedRock(c, n, ntab, unitsOf(c.units)); ch = true; }
        ImGui::PopID(); return ch;
    }
    struct R { int idx = 1; double rho = 2650; double mx = 0; };
    std::vector<R> rows;
    for (const auto& r : records(k->body)) {
        R x; double v;
        if (alk) { if (!r.empty() && num(r[0], v)) x.mx = v; }
        else { if (r.size() > 0 && num(r[0], v)) x.idx = int(v); if (r.size() > 1 && num(r[1], v)) x.rho = v; }
        if (r.empty() && !rows.empty()) x = rows.back();
        rows.push_back(x);
    }
    while (int(rows.size()) < ntab) rows.push_back(rows.empty() ? R{} : rows.back());
    ui::RecordTable<R> t;
    if (alk) t.num("max adsorption", [](R& x) -> double& { return x.mx; });
    else t.choice("adsorption", [](R& x) -> int& { return x.idx; }, {"reversible", "irreversible"}, 1)
          .num("rock density", [](R& x) -> double& { return x.rho; });
    t.fixedRows();
    if (t.draw("rk", rows, ImGui::GetFontSize() * (2.5f + float(rows.size())))) {
        std::string b; char buf[64];
        for (auto& x : rows) { if (alk) std::snprintf(buf, sizeof buf, "  %.6g /\n", x.mx); else std::snprintf(buf, sizeof buf, "  %d  %.6g /\n", x.idx, x.rho); b += buf; }
        k->body = b; ch = true;
    }
    if (ImGui::SmallButton("Remove")) { removeKw(c, n); ch = true; }
    ImGui::PopID();
    return ch;
}
// injection records: WSURFACT / WALKALIN / WFOAM per report step
struct InjRow { int step = 0; int kind = 0; std::string well; double c = 0; };
const char* kInjKw[4] = {"WSURFACT", "WALKALIN", "WFOAM", "WSOLVENT"};
std::vector<InjRow> readInj(const Case& c) {
    std::vector<InjRow> v;
    for (const auto& k : c.extraKeywords)
        for (int t = 0; t < 4; ++t)
            if (k.name == kInjKw[t])
                for (const auto& r : records(k.body)) {
                    if (r.empty()) continue;
                    InjRow x; x.step = k.step; x.kind = t; x.well = r[0];
                    if (x.well.size() >= 2 && x.well.front() == '\'') x.well = x.well.substr(1, x.well.size() - 2);
                    if (r.size() > 1) num(r[1], x.c);
                    v.push_back(x);
                }
    std::stable_sort(v.begin(), v.end(), [](const InjRow& a, const InjRow& b) { return a.step < b.step; });
    return v;
}
void writeInj(Case& c, const std::vector<InjRow>& v) {
    c.extraKeywords.erase(std::remove_if(c.extraKeywords.begin(), c.extraKeywords.end(), [](const ExtraKeyword& k) {
                              return k.name == "WSURFACT" || k.name == "WALKALIN" || k.name == "WFOAM" || k.name == "WSOLVENT"; }), c.extraKeywords.end());
    std::map<std::pair<int, int>, std::string> body;
    for (const auto& x : v) {
        char b[160]; std::snprintf(b, sizeof b, " '%s'  %.6g /\n", x.well.c_str(), x.c);
        body[{x.step, x.kind}] += b;
    }
    for (const auto& [key, b] : body) {
        ExtraKeyword k; k.name = kInjKw[key.second]; k.section = DeckSection::Schedule; k.step = key.first; k.body = b;
        k.beforeControls = false; k.note = key.second == 3 ? "solvent (C35 editor)" : "chemical EOR (C35 editor)";
        c.extraKeywords.push_back(k);
    }
}
// SOLVENT: PVDS (p, Bs, mu_s per PVTNUM), SDENSITY (surface density per PVTNUM), SSFN (Fs, krg mult, krs mult per SATNUM)
std::vector<double> seedSolv(const std::string& n, const ChemUnits& u) {
    if (n == "PVDS") {   // a CO2-like solvent (field: psia, rb/Mscf, cP)
        const double pp[5] = {400, 1000, 2000, 4000, 8000}, bb[5] = {7.5, 2.8, 1.2, 0.55, 0.38}, mm[5] = {0.016, 0.019, 0.03, 0.06, 0.09};
        std::vector<double> t;
        for (int i = 0; i < 5; ++i) { t.push_back(pp[i] * u.p); t.push_back(bb[i] * (u.field ? 1.0 : 0.1781076)); t.push_back(mm[i]); }
        return t;
    }
    if (n == "SSFN") return {0, 0, 0, 1, 1, 1};
    return {};
}
bool table3Editor(Case& c, const char* name, const char* const cols[3], const char* tip, int ntab, int& reg, bool perPvt) {
    ExtraKeyword* k = findKw(c, name);
    ImGui::PushID(name);
    bool ch = false;
    ImGui::TextUnformatted(name);
    ImGui::SameLine(); ui::hint(tip);
    if (!k) {
        ImGui::TextDisabled("not given"); ImGui::SameLine();
        if (ImGui::SmallButton("Add with typical data")) {
            std::vector<std::vector<double>> t(size_t(std::max(1, ntab)), seedSolv(name, unitsOf(c.units)));
            setKw(c, name, DeckSection::Props, writeTables(t, 3)); ch = true;
        }
        ImGui::PopID(); return ch;
    }
    std::vector<std::vector<double>> tabs;
    if (!parseTables(k->body, 3, tabs) || tabs.empty()) {
        ImGui::TextDisabled("non-numeric records - edit the deck text (Deck tab)");
        ImGui::PopID(); return ch;
    }
    while (int(tabs.size()) < ntab) tabs.push_back(tabs.back());
    if (ntab > 1) {
        reg = std::clamp(reg, 0, ntab - 1);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::SliderInt("##reg", &reg, 0, ntab - 1, perPvt ? "PVTNUM %d" : "SATNUM %d");
        ImGui::SameLine();
    } else reg = 0;
    if (ImGui::SmallButton("Remove")) { removeKw(c, name); ImGui::PopID(); return true; }
    std::vector<double>& t = tabs[size_t(std::min(reg, int(tabs.size()) - 1))];
    std::vector<Row3> rows;
    for (size_t i = 0; i + 2 < t.size(); i += 3) rows.push_back({t[i], t[i + 1], t[i + 2]});
    const float h = ImGui::GetFontSize() * 7.5f;
    if (ImGui::BeginTable("tp", 2, ImGuiTableFlags_SizingStretchSame)) {
        ImGui::TableNextColumn();
        if (ui::editTable<Row3>("t", rows, {cols[0], cols[1], cols[2]}, [](Row3& r, int col) -> double* { return &r[size_t(col)]; }, h)) {
            std::sort(rows.begin(), rows.end(), [](const Row3& a, const Row3& b) { return a[0] < b[0]; });
            t.clear();
            for (auto& r : rows) for (double v : r) t.push_back(v);
            k->body = writeTables(tabs, 3);
            ch = true;
        }
        ImGui::TableNextColumn();
        std::vector<double> x, y1, y2;
        for (auto& r : rows) { x.push_back(r[0]); y1.push_back(r[1]); y2.push_back(r[2]); }
        if (ImPlot::BeginPlot("##plot", ImVec2(-1, h), ImPlotFlags_NoMenus)) {
            ImPlot::SetupAxes(cols[0], nullptr, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, nullptr, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Opposite);
            ImPlotSpec a; a.LineColor = col4(Pal::Accent); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle;
            ImPlotSpec b = a; b.LineColor = col4(Pal::Warn);
            if (!x.empty()) {
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1); ImPlot::PlotLine(cols[1], x.data(), y1.data(), int(x.size()), a);
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2); ImPlot::PlotLine(cols[2], x.data(), y2.data(), int(x.size()), b);
            }
            ImPlot::EndPlot();
        }
        ImGui::EndTable();
    }
    ImGui::PopID();
    return ch;
}
// SDENSITY: one surface density per PVTNUM region
bool sdensityEditor(Case& c, int npvt) {
    ExtraKeyword* k = findKw(c, "SDENSITY");
    ImGui::PushID("SDENSITY");
    bool ch = false;
    const bool field = c.units == UnitSystem::Field;
    ImGui::TextUnformatted("SDENSITY");
    ImGui::SameLine(); ui::hint("Solvent density at surface conditions per PVTNUM region (field lb/ft3, metric kg/m3).");
    if (!k) {
        ImGui::TextDisabled("not given"); ImGui::SameLine();
        if (ImGui::SmallButton("Add")) {
            std::string b; char buf[48];
            for (int r = 0; r < npvt; ++r) { std::snprintf(buf, sizeof buf, "  %.6g /\n", field ? 0.1167 : 1.87); b += buf; }   // CO2
            setKw(c, "SDENSITY", DeckSection::Props, b); ch = true;
        }
        ImGui::PopID(); return ch;
    }
    struct R { double rho = 0.1167; };
    std::vector<R> rows;
    for (const auto& r : records(k->body)) { R x; double v; if (!r.empty() && num(r[0], v)) x.rho = v; else if (!rows.empty()) x = rows.back(); rows.push_back(x); }
    while (int(rows.size()) < npvt) rows.push_back(rows.empty() ? R{} : rows.back());
    ui::RecordTable<R> t;
    t.num(field ? "surface density (lb/ft3)" : "surface density (kg/m3)", [](R& x) -> double& { return x.rho; });
    t.fixedRows();
    if (t.draw("sd", rows, ImGui::GetFontSize() * (2.5f + float(rows.size())))) {
        std::string b; char buf[48];
        for (auto& x : rows) { std::snprintf(buf, sizeof buf, "  %.6g /\n", x.rho); b += buf; }
        k->body = b; ch = true;
    }
    ImGui::PopID();
    return ch;
}
} // namespace

void App::editChemEOR() {
    Case& C = case_;
    bool ch = false;
    const ChemUnits u = unitsOf(C.units);
    const int npvt = std::max(1, C.props.numPvtRegions()), nsat = std::max(1, C.props.numSatRegions());
    ui::header("Chemical / solvent EOR: surfactant, alkaline, foam, solvent");
    ui::hint("ECLIPSE SURFACT / ALKALINE / FOAM options in the own engine: transport of each chemical in its carrier phase with "
             "adsorption, decay (foam) and the coupled effects on the flow - surfactant viscosity, IFT -> capillary number -> "
             "miscible relative permeability blend (SURFCAPD, SURFNUM curves) and Pc scaling, alkaline multipliers on IFT and "
             "adsorption, foam gas-mobility reduction. Solved sequentially-implicitly with outer coupling, mass-conservative "
             "(see ENGINES_AND_MATH 6.28). Solvent: a fourth (gas-like) pseudo-component sharing the gas saturation, "
             "immiscible SSFN model (see 6.29).");
    bool on[4] = {findKw(C, "SURFACT") || findKw(C, "SURFACTW"), findKw(C, "ALKALINE") != nullptr, findKw(C, "FOAM") != nullptr,
                  findKw(C, "SOLVENT") != nullptr};
    const char* opt[4] = {"Surfactant in water (SURFACT)", "Alkaline in water (ALKALINE)", "Foam (FOAM)", "Solvent (SOLVENT)"};
    for (int t = 0; t < 4; ++t) {
        if (t) ImGui::SameLine();
        bool v = on[t];
        if (ImGui::Checkbox(opt[t], &v)) {
            const char* kw = t == 0 ? "SURFACT" : t == 1 ? "ALKALINE" : t == 2 ? "FOAM" : "SOLVENT";
            if (v) {
                setKw(C, kw, DeckSection::Runspec, "");
                if (t == 0) {
                    for (auto& K : kSurf) seedKw(C, K.name, K.perPvt ? npvt : nsat, u);
                    seedRock(C, "SURFROCK", nsat, u);
                    if (C.grid.surfnum.empty()) notify("Surfactant: choose the miscible saturation table (SURFNUM) in the Surfactant tab");
                } else if (t == 1) { seedKw(C, "ALSURFST", npvt, u); seedKw(C, "ALSURFAD", nsat, u); }
                else if (t == 2) { seedKw(C, "FOAMMOB", npvt, u); if (!findKw(C, "FOAMOPTS")) setKw(C, "FOAMOPTS", DeckSection::Props, " GAS  TAB /\n"); }
                else {
                    for (const char* n : {"PVDS", "SSFN"})
                        if (!findKw(C, n)) {
                            std::vector<std::vector<double>> tb(size_t(std::string(n) == "PVDS" ? npvt : nsat), seedSolv(n, u));
                            setKw(C, n, DeckSection::Props, writeTables(tb, 3));
                        }
                    if (!findKw(C, "SDENSITY")) {
                        std::string b; char buf[48];
                        for (int r = 0; r < npvt; ++r) { std::snprintf(buf, sizeof buf, "  %.6g /\n", u.field ? 0.1167 : 1.87); b += buf; }
                        setKw(C, "SDENSITY", DeckSection::Props, b);
                    }
                    if (!C.phases.gas || !C.phases.oil) notify("Solvent needs oil and gas phases (gas-water solvent runs with the OPM Flow engine)");
                }
            } else {
                removeKw(C, kw);
                if (t == 0) removeKw(C, "SURFACTW");
            }
            on[t] = v; ch = true;
        }
    }
    if (!on[0] && !on[1] && !on[2] && !on[3]) {
        ui::hint("Brine / low salinity and the miscible solvent model (MISC, TLMIXPAR) are on the roadmap; decks using them are "
                 "kept verbatim and run with the OPM Flow engine (Run > Engine).");
        if (ch) markDirty();
        return;
    }
    static int reg[16] = {0};
    if (ImGui::BeginTabBar("chem")) {
        if (on[0] && ui::tabItem("Surfactant")) {
            int i = 0;
            for (auto& K : kSurf) { ch |= tableEditor(C, K, K.perPvt ? npvt : nsat, reg[i++]); ImGui::Separator(); }
            ch |= rockEditor(C, "SURFROCK", nsat);
            ImGui::Separator();
            ImGui::TextUnformatted("SURFNUM - miscible saturation table");
            ImGui::SameLine(); ui::hint("Saturation table used for the miscible (high capillary number) curves, per cell. "
                                        "Set one table for the whole model or assign boxes of cells.");
            const int n = C.grid.size();
            auto& sn = C.grid.surfnum;
            {
                std::vector<int> cnt(size_t(nsat) + 1, 0);
                for (int v : sn) if (v >= 1 && v <= nsat) ++cnt[size_t(v)]; else ++cnt[0];
                if (sn.empty()) ImGui::TextDisabled("not set: no miscible curves (surfactant only changes viscosity and Pc)");
                else for (int k = 1; k <= nsat; ++k) { ImGui::Text("table %d: %d cells", k, cnt[size_t(k)]); if (k < nsat) ImGui::SameLine(0, 18); }
                if (!sn.empty() && cnt[0]) ImGui::TextColored(col4(Pal::Err), "%d cells refer to a missing table", cnt[0]);
            }
            static int all = 2;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6); ImGui::InputInt("##snall", &all, 0); all = std::clamp(all, 1, nsat);
            ImGui::SameLine();
            if (ImGui::Button("All cells")) { sn.assign(size_t(n), all); ch = true; }
            ImGui::SameLine();
            if (ImGui::Button("Clear SURFNUM")) { sn.clear(); ch = true; }
            static int box[7] = {1, 1, 1, 1, 1, 1, 2};
            const GridInput& g = C.grid;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9); ImGui::InputInt2("I1, I2##sn", &box[0]);
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9); ImGui::InputInt2("J1, J2##sn", &box[2]);
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9); ImGui::InputInt2("K1, K2##sn", &box[4]);
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5); ImGui::InputInt("table##snb", &box[6], 0);
            box[6] = std::clamp(box[6], 1, nsat);
            if (ImGui::Button("Assign box##sn")) {
                if (sn.size() != size_t(n)) sn.assign(size_t(n), 1);
                const int i1 = std::clamp(std::min(box[0], box[1]), 1, g.nx), i2 = std::clamp(std::max(box[0], box[1]), 1, g.nx);
                const int j1 = std::clamp(std::min(box[2], box[3]), 1, g.ny), j2 = std::clamp(std::max(box[2], box[3]), 1, g.ny);
                const int k1 = std::clamp(std::min(box[4], box[5]), 1, g.nz), k2 = std::clamp(std::max(box[4], box[5]), 1, g.nz);
                for (int k = k1; k <= k2; ++k) for (int j = j1; j <= j2; ++j) for (int i = i1; i <= i2; ++i)
                    sn[size_t((i - 1) + g.nx * ((j - 1) + g.ny * (k - 1)))] = box[6];
                ch = true;
            }
            if (nsat < 2) ui::hint("Add a second saturation table (Properties > Saturation functions) for the miscible curves.");
            ImGui::EndTabItem();
        }
        if (on[1] && ui::tabItem("Alkaline")) {
            int i = 4;
            for (auto& K : kAlk) { ch |= tableEditor(C, K, K.perPvt ? npvt : nsat, reg[i++]); ImGui::Separator(); }
            ch |= rockEditor(C, "ALKROCK", nsat);
            ImGui::EndTabItem();
        }
        if (on[2] && ui::tabItem("Foam")) {
            ExtraKeyword* fo = findKw(C, "FOAMOPTS");
            std::string ph = "GAS", md = "TAB";
            if (fo) { auto r = records(fo->body); if (!r.empty() && !r[0].empty()) ph = r[0][0]; if (!r.empty() && r[0].size() > 1) md = r[0][1]; }
            int iph = ph.find("WAT") != std::string::npos ? 1 : 0;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
            if (ph.find("SOLV") != std::string::npos) iph = 2;
            if (ImGui::Combo("Foam transport phase (FOAMOPTS)", &iph, "gas\0water\0solvent\0")) {
                setKw(C, "FOAMOPTS", DeckSection::Props, std::string(iph == 1 ? " WATER" : iph == 2 ? " SOLVENT" : " GAS") + "  " + md + " /\n"); ch = true;
            }
            if (iph == 2 && !on[3]) ui::hint("FOAMOPTS SOLVENT needs the SOLVENT option.");
            if (md.find("FUNC") != std::string::npos) ui::hint("FOAMOPTS FUNC (functional model) is read as the tabular model in the own engine.");
            int i = 8;
            for (auto& K : kFoam) { ch |= tableEditor(C, K, K.perPvt ? npvt : nsat, reg[i++]); ImGui::Separator(); }
            ch |= rockEditor(C, "FOAMROCK", nsat);
            ImGui::Separator();
            ImGui::TextUnformatted("SFOAM - initial foam concentration");
            ImGui::SameLine();
            double s0 = C.grid.sfoam.empty() ? 0.0 : C.grid.sfoam[0];
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            if (ImGui::InputDouble("##sf", &s0, 0, 0, "%.4g")) {
                if (s0 <= 0) C.grid.sfoam.clear(); else C.grid.sfoam.assign(size_t(C.grid.size()), s0);
                ch = true;
            }
            ImGui::EndTabItem();
        }
        if (on[3] && ui::tabItem("Solvent")) {
            const bool fd = C.units == UnitSystem::Field;
            static const char* pvdsF[3] = {"pressure (psia)", "Bs (rb/Mscf)", "viscosity (cP)"};
            static const char* pvdsM[3] = {"pressure (bar)", "Bs (rm3/sm3)", "viscosity (cP)"};
            static const char* ssfn[3] = {"solvent fraction Fs", "gas kr multiplier", "solvent kr multiplier"};
            ch |= table3Editor(C, "PVDS", fd ? pvdsF : pvdsM, "Solvent formation volume factor and viscosity vs pressure (dry-gas style).", npvt, reg[14], true);
            ImGui::Separator();
            ch |= sdensityEditor(C, npvt);
            ImGui::Separator();
            ch |= table3Editor(C, "SSFN", ssfn,
                               "Immiscible split of the gas-phase relative permeability: krg(Sg+Ss) x SSFN_g(Sg/(Sg+Ss)) for gas and "
                               "krg(Sg+Ss) x SSFN_s(Ss/(Sg+Ss)) for solvent. A straight line (0 0 0 / 1 1 1) is the classical mixed-gas case.",
                               nsat, reg[15], false);
            ImGui::Separator();
            ImGui::TextUnformatted("SSOL - initial solvent saturation");
            ImGui::SameLine();
            double s0 = C.grid.ssol.empty() ? 0.0 : C.grid.ssol[0];
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            if (ImGui::InputDouble("##ss", &s0, 0, 0, "%.4g")) {
                if (s0 <= 0) C.grid.ssol.clear(); else C.grid.ssol.assign(size_t(C.grid.size()), std::min(1.0, s0));
                ch = true;
            }
            for (const char* n : {"MISCIBLE", "MISC", "TLMIXPAR", "PMISC", "SORWMIS", "SGCWMIS", "MSFN", "SOF2"})
                if (findKw(C, n)) {
                    ImGui::TextColored(col4(Pal::Warn), "%s present: the miscible (Todd-Longstaff) model is not in the own engine - it runs the "
                                       "immiscible SSFN model (the OPM Flow engine runs the miscible model).", n);
                    break;
                }
            ui::hint("Inject solvent with WSOLVENT (Injection tab): the fraction of the gas injection rate that is solvent.");
            ImGui::EndTabItem();
        }
        if (ui::tabItem("Injection")) {
            std::vector<InjRow> rows = readInj(C);
            std::vector<std::string> wells, kinds;
            for (auto& w : C.wells) wells.push_back(w.name);
            for (int t = 0; t < 4; ++t) kinds.push_back(kInjKw[t]);
            ui::hint("Injected concentration per well and report step (WSURFACT / WALKALIN / WFOAM), or the solvent fraction of the "
                     "gas injection rate (WSOLVENT, 0-1). A value holds until it is changed; set 0 to stop (chase water / gas). "
                     "Step 1 = before the first TSTEP / DATES.");
            ui::RecordTable<InjRow> t;
            t.integer("Step", [](InjRow& x) -> int& { return x.step; }, 1)
             .choice("Keyword", [](InjRow& x) -> int& { return x.kind; }, kinds)
             .choiceStr("Well", [](InjRow& x) -> std::string& { return x.well; }, wells)
             .num(C.units == UnitSystem::Field ? "Concentration (lb/stb; foam: lb/stb or lb/Mscf)" : "Concentration (kg/sm3)", [](InjRow& x) -> double& { return x.c; });
            t.onChange = [&](InjRow& x) {
                x.step = std::clamp(x.step, 0, std::max(0, int(C.schedule.size())));
                x.c = x.kind == 3 ? std::clamp(x.c, 0.0, 1.0) : std::max(0.0, x.c);
            };
            t.newRow = [&](int at) {
                InjRow x = rows.empty() ? InjRow{} : rows[size_t(std::min(at, int(rows.size()) - 1))];
                if (x.well.empty() && !wells.empty()) x.well = wells[0];
                if (!on[x.kind]) x.kind = on[0] ? 0 : on[1] ? 1 : on[2] ? 2 : 3;
                return x;
            };
            if (t.draw("inj", rows, ImGui::GetFontSize() * 14)) { writeInj(C, rows); ch = true; }
            if (rows.empty() && ui::accentButton("Add an injection record")) {
                InjRow x; x.kind = on[0] ? 0 : on[1] ? 1 : on[2] ? 2 : 3;
                for (auto& w : C.wells) if (w.type != WellType::Producer) { x.well = w.name; break; }
                rows.push_back(x); writeInj(C, rows); ch = true;
            }
            ImGui::EndTabItem();
        }
        if (ui::tabItem("Results")) {
            const std::vector<std::string> names = ecl::componentNames(C);
            std::vector<double> tt;
            std::vector<std::vector<double>> it, pt;
            std::vector<std::string> nm;
            std::vector<double> sIT, sPT, sIP;   // solvent totals / in place
            {
                std::lock_guard<std::mutex> l(runner_.mutex());
                const RunResults& R = runner_.results();
                if (R.valid && R.summary.size() > 1) {
                    for (size_t k = 0; k < names.size(); ++k)
                        if (names[k] == "SURFACT" || names[k] == "ALKALINE" || names[k] == "FOAM") { nm.push_back(names[k]); it.emplace_back(); pt.emplace_back(); }
                    const size_t first = names.size() - nm.size();
                    for (auto& sp : R.summary) {
                        tt.push_back(sp.time / 365.25);
                        if (on[3]) { sIT.push_back(sp.fnit); sPT.push_back(sp.fnpt); sIP.push_back(sp.fnip); }
                        for (size_t j = 0; j < nm.size(); ++j) {
                            const size_t q = first + j;
                            it[j].push_back(q < sp.fcompIT.size() ? sp.fcompIT[q] : 0.0);
                            pt[j].push_back(q < sp.fcompPT.size() ? sp.fcompPT[q] : 0.0);
                        }
                    }
                }
            }
            if (!tt.empty() && sIT.size() == tt.size()) {
                const char* gu = C.units == UnitSystem::Field ? "solvent (Mscf)" : "solvent (sm3)";
                if (ImPlot::BeginPlot("Solvent", ImVec2(-1, ImGui::GetFontSize() * 14))) {
                    ImPlot::SetupAxes("years", gu, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec a; a.LineColor = col4(Pal::Accent); a.LineWeight = 2;
                    ImPlotSpec b = a; b.LineColor = col4(Pal::Warn);
                    ImPlotSpec d = a; d.LineColor = col4(Pal::Gas); d.LineWeight = 1.5f;
                    ImPlot::PlotLine("injected (FNIT)", tt.data(), sIT.data(), int(tt.size()), a);
                    ImPlot::PlotLine("produced (FNPT)", tt.data(), sPT.data(), int(tt.size()), b);
                    ImPlot::PlotLine("in place (FNIP)", tt.data(), sIP.data(), int(tt.size()), d);
                    ImPlot::EndPlot();
                }
                ui::hint("Solvent: 3-D maps SSOLVENT (saturation) and SOLVFRAC (fraction of the gas-like phase); summary vectors "
                         "FNIR/FNPR/FNIT/FNPT/FNIP and WNIR/WNPR/WNIT/WNPT per well.");
            }
            if (tt.empty() || (nm.empty() && sIT.empty())) ui::hint("Run the case to see the injected / produced chemical totals.");
            else if (nm.empty()) {}
            else {
                const char* mu = C.units == UnitSystem::Field ? "lb" : "kg";
                char ya[48]; std::snprintf(ya, sizeof ya, "cumulative (%s)", mu);
                if (ImPlot::BeginPlot("Chemical totals", ImVec2(-1, ImGui::GetFontSize() * 14))) {
                    ImPlot::SetupAxes("years", ya, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    const ImU32 cols[3] = {Pal::Accent, Pal::Warn, Pal::Gas};
                    for (size_t j = 0; j < nm.size(); ++j) {
                        ImPlotSpec a; a.LineColor = col4(cols[j % 3]); a.LineWeight = 2;
                        ImPlotSpec b = a; b.LineWeight = 1.0f; b.Marker = ImPlotMarker_Square; b.MarkerSize = 2;
                        ImPlot::PlotLine((nm[j] + " injected").c_str(), tt.data(), it[j].data(), int(tt.size()), a);
                        ImPlot::PlotLine((nm[j] + " produced").c_str(), tt.data(), pt[j].data(), int(tt.size()), b);
                    }
                    ImPlot::EndPlot();
                }
                ui::hint("3-D maps (Results / viewer property list): SURFACT, ALKALINE, FOAM concentrations, SURFCNM (log10 capillary "
                         "number), SURFMISC (miscibility F), SURFSTM (IFT ratio), FOAMMRF (gas mobility factor). Summary vectors: "
                         "FTPR/FTPT/FTIR/FTIT + SUR/ALK/FOA, FSPR/FSPT/FSIR/FSIT, WTPR/WTIR per well.");
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    if (ch) markDirty();
}
