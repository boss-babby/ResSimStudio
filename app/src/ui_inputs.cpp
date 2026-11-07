// Inspector editors for every section of the case (no deck editing needed).
#include "plotkit.h"
#include "app.h"
#include <functional>
#include <map>
#include <set>
#include "widgets.h"
#include "i18n.h"
#include "implot.h"
#include "ressim/fluid.h"
#include "ressim/tools.h"
#include "ressim/geostat.h"
#include "ressim/geo_analysis.h"
#include "icons.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <thread>
#include <random>

using namespace ressim;

static void clampWellsToGrid(Case& c) {
    const GridInput& g = c.grid;
    auto dims = [&](int l, int& nx, int& ny, int& nz) {
        if (l >= 0 && l < int(g.lgrs.size())) { nx = g.lgrs[size_t(l)].nx; ny = g.lgrs[size_t(l)].ny; nz = g.lgrs[size_t(l)].nz; }
        else { nx = g.nx; ny = g.ny; nz = g.nz; }
    };
    for (auto& w : c.wells) {
        int nx, ny, nz; dims(w.lgr, nx, ny, nz);
        w.i = std::max(0, std::min(nx - 1, w.i)); w.j = std::max(0, std::min(ny - 1, w.j));
        for (auto& cp : w.completions) {
            dims(cp.lgr, nx, ny, nz);
            cp.i = std::max(0, std::min(nx - 1, cp.i)); cp.j = std::max(0, std::min(ny - 1, cp.j));
            cp.k1 = std::max(0, std::min(nz - 1, cp.k1)); cp.k2 = std::max(cp.k1, std::min(nz - 1, cp.k2));
        }
    }
}


// ============================================================================ GRID: local grid refinements (CARFIN)
namespace {
std::string intsToText(const std::vector<int>& v) { std::string s; for (int x : v) { if (!s.empty()) s += ' '; s += std::to_string(x); } return s; }
std::string dblsToText(const std::vector<double>& v) {
    std::string s; char b[32];
    for (double x : v) { std::snprintf(b, sizeof b, "%.6g", x); if (!s.empty()) s += ' '; s += b; }
    return s;
}
// "1 2 3", also accepts ECLIPSE repeat counts "3*2"
template <class T> std::vector<T> textToNums(const std::string& s) {
    std::vector<T> out;
    size_t p = 0;
    while (p < s.size()) {
        while (p < s.size() && (s[p] == ' ' || s[p] == ',' || s[p] == '\t')) ++p;
        if (p >= s.size()) break;
        size_t e = p; while (e < s.size() && s[e] != ' ' && s[e] != ',' && s[e] != '\t') ++e;
        const std::string tok = s.substr(p, e - p);
        const size_t star = tok.find('*');
        int rep = 1; std::string val = tok;
        if (star != std::string::npos) { rep = std::atoi(tok.substr(0, star).c_str()); val = tok.substr(star + 1); }
        const double v = std::atof(val.c_str());
        for (int r = 0; r < std::max(rep, 0); ++r) out.push_back(T(v));
        p = e;
    }
    return out;
}
}  // namespace

void App::editLgrs() {
    GridInput& g = case_.grid;
    ui::header("Local grid refinement (LGR)");
    ui::hint("Cartesian refinements (ECLIPSE CARFIN / ENDFIN, OPM Flow compatible): every host cell in the box is replaced by "
             "its children. Children inherit the host rock and regions unless overridden below; host-child and child-child "
             "transmissibilities are computed from the refined geometry, host connections (faults, MULTFLT, NNCs, pinch-outs) are "
             "redistributed over the children so the coarse flux capacity is preserved. Wells inside use WELSPECL / COMPDATL.");
    if (g.dual) { ui::hint("Dual porosity model: local grid refinements are not available in ResSim (ECLIPSE only)."); return; }
    static int sel = 0;
    sel = std::max(0, std::min(int(g.lgrs.size()) - 1, sel));
    int del = -1;
    if (!g.lgrs.empty() &&
        ui::beginTable("lgrs", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ui::tableColumn("Name"); ui::tableColumn("Host box I J K"); ui::tableColumn("Local NX NY NZ");
        ui::tableColumn("Cells"); ui::tableColumn("Wells"); ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
        ImGui::TableHeadersRow();
        for (int l = 0; l < int(g.lgrs.size()); ++l) {
            const LgrInput& L = g.lgrs[size_t(l)];
            ImGui::PushID(l);
            ImGui::TableNextRow(); ImGui::TableNextColumn();
            if (ImGui::Selectable(L.name.c_str(), sel == l, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) sel = l;
            ImGui::TableNextColumn(); ImGui::Text("%d-%d  %d-%d  %d-%d", L.i1 + 1, L.i2 + 1, L.j1 + 1, L.j2 + 1, L.k1 + 1, L.k2 + 1);
            ImGui::TableNextColumn(); ImGui::Text("%d x %d x %d", L.nx, L.ny, L.nz);
            ImGui::TableNextColumn(); ImGui::Text("%d", L.size());
            ImGui::TableNextColumn();
            std::string wl; for (auto& w : case_.wells) if (w.lgr == l) { if (!wl.empty()) wl += ", "; wl += w.name; }
            ImGui::TextUnformatted(wl.empty() ? "-" : wl.c_str());
            ImGui::TableNextColumn(); if (ImGui::SmallButton(TRL("\xC3\x97"))) del = l;
            ImGui::PopID();
        }
        ui::endTable();
    }
    if (del >= 0) { removeLgr(case_, del); markDirty(true); notify("Local grid removed (its wells moved back to the global grid)"); return; }

    // ---- add
    if (ImGui::TreeNode(TRL("Add a refinement"))) {
        static int mode = 0, wsel = 0, half = 1, rxy = 3, rz = 1; static bool compLayers = true, moveWell = true;
        static int box[6] = {1, 1, 1, 1, 1, 1}, ref[3] = {3, 3, 1};
        static const char* modes[] = {"Around a well", "Host box"};
        if (ui::beginProps("lgradd")) {
            ui::propCombo("Placement", &mode, modes, 2);
            if (mode == 0) {
                std::vector<const char*> names; for (auto& w : case_.wells) if (w.lgr < 0) names.push_back(w.name.c_str());
                if (names.empty()) ui::propText("Well", "no global wells");
                else {
                    wsel = std::max(0, std::min(int(names.size()) - 1, wsel));
                    ui::propCombo("Well", &wsel, names.data(), int(names.size()));
                }
                ui::propI("Host cells around the well", &half, 0, 20, "Half width of the host box in I and J (1 = 3 x 3 host columns).");
                ui::propI("Refinement in I and J", &rxy, 1, 25, "Children per host cell in I and in J (odd numbers keep the well in the centre child).");
                ui::propI("Refinement in K", &rz, 1, 25);
                ui::propBool("Completed layers only", &compLayers, "Host box spans the completed layers; off = all layers.");
                ui::propBool("Move the well into the LGR", &moveWell, "Converts its WELSPECS / COMPDAT to WELSPECL / COMPDATL.");
            } else {
                ui::label("I1 I2 (1-based)"); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputInt2("##bi", box);
                ui::label("J1 J2"); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputInt2("##bj", box + 2);
                ui::label("K1 K2"); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputInt2("##bk", box + 4);
                ui::label("Children per host cell I J K"); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputInt3("##br", ref);
            }
            ui::endProps();
        }
        if (ui::accentButton("Add LGR")) {
            try {
                if (mode == 0) {
                    int wi = -1, c = 0;
                    for (int w = 0; w < int(case_.wells.size()); ++w) if (case_.wells[size_t(w)].lgr < 0) { if (c == wsel) wi = w; ++c; }
                    if (wi < 0) throw std::runtime_error("no well selected");
                    LgrInput L = makeLgrAroundWell(case_, wi, half, rxy, rz, compLayers);
                    g.lgrs.push_back(L);
                    if (moveWell) moveWellToLgr(case_, wi, int(g.lgrs.size()) - 1);
                } else {
                    LgrInput L;
                    L.i1 = std::max(0, std::min(box[0], box[1]) - 1); L.i2 = std::min(g.nx - 1, std::max(box[0], box[1]) - 1);
                    L.j1 = std::max(0, std::min(box[2], box[3]) - 1); L.j2 = std::min(g.ny - 1, std::max(box[2], box[3]) - 1);
                    L.k1 = std::max(0, std::min(box[4], box[5]) - 1); L.k2 = std::min(g.nz - 1, std::max(box[4], box[5]) - 1);
                    L.nx = (L.i2 - L.i1 + 1) * std::max(1, ref[0]); L.ny = (L.j2 - L.j1 + 1) * std::max(1, ref[1]); L.nz = (L.k2 - L.k1 + 1) * std::max(1, ref[2]);
                    for (int s = int(g.lgrs.size()) + 1; ; ++s) {
                        L.name = "LGR" + std::to_string(s);
                        bool clash = false; for (auto& o : g.lgrs) if (o.name == L.name) clash = true;
                        if (!clash) break;
                    }
                    for (auto& o : g.lgrs)
                        if (L.i1 <= o.i2 && o.i1 <= L.i2 && L.j1 <= o.j2 && o.j1 <= L.j2 && L.k1 <= o.k2 && o.k1 <= L.k2)
                            throw std::runtime_error("the box overlaps " + o.name);
                    g.lgrs.push_back(L);
                }
                sel = int(g.lgrs.size()) - 1;
                markDirty(true);
                notify("Local grid " + g.lgrs.back().name + " added");
            } catch (const std::exception& e) { notify(std::string("LGR: ") + e.what()); }
        }
        ImGui::TreePop();
    }
    if (g.lgrs.empty()) return;

    // ---- selected refinement
    LgrInput& L = g.lgrs[size_t(sel)];
    ImGui::PushID("lgrsel");
    ImGui::SeparatorText(("Refinement " + L.name).c_str());
    bool ch = false;
    if (ui::beginProps("lgred")) {
        ui::label("Name", "Up to 8 characters in ECLIPSE.");
        std::string nm = L.name;
        if (ui::inputString("##lname", &nm, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsNoBlank)) {
            bool clash = false; for (auto& o : g.lgrs) if (&o != &L && o.name == nm) clash = true;
            if (!nm.empty() && !clash) { L.name = nm; ch = true; }
        }
        int bi[2] = {L.i1 + 1, L.i2 + 1}, bj[2] = {L.j1 + 1, L.j2 + 1}, bk[2] = {L.k1 + 1, L.k2 + 1};
        int hn[3] = {L.i2 - L.i1 + 1, L.j2 - L.j1 + 1, L.k2 - L.k1 + 1};
        auto box = [&](const char* lab, int* v, int& a, int& b, int nmax) {
            ui::label(lab); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputInt2((std::string("##") + lab).c_str(), v, ImGuiInputTextFlags_EnterReturnsTrue)) {
                a = std::max(0, std::min(nmax - 1, std::min(v[0], v[1]) - 1)); b = std::max(a, std::min(nmax - 1, std::max(v[0], v[1]) - 1)); return true;
            }
            return false;
        };
        bool bch = box("Host I1 I2", bi, L.i1, L.i2, g.nx);
        bch |= box("Host J1 J2", bj, L.j1, L.j2, g.ny);
        bch |= box("Host K1 K2", bk, L.k1, L.k2, g.nz);
        if (bch) {   // keep the refinement ratio, drop per-host splits that no longer fit
            const int nh[3] = {L.i2 - L.i1 + 1, L.j2 - L.j1 + 1, L.k2 - L.k1 + 1};
            L.nx = std::max(1, L.nx / std::max(1, hn[0])) * nh[0]; L.ny = std::max(1, L.ny / std::max(1, hn[1])) * nh[1]; L.nz = std::max(1, L.nz / std::max(1, hn[2])) * nh[2];
            L.nxfin.clear(); L.nyfin.clear(); L.nzfin.clear(); L.hxfin.clear(); L.hyfin.clear(); L.hzfin.clear();
            L.props.clear();
            ch = true;
        }
        int nn[3] = {L.nx, L.ny, L.nz};
        ui::label("Local NX NY NZ", "Total number of local cells in each direction (CARFIN items 8-10). Without NXFIN they are spread "
                                    "evenly over the host cells, so NX must be a multiple of the number of host columns.");
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::InputInt3("##lnn", nn, ImGuiInputTextFlags_EnterReturnsTrue)) {
            L.nx = std::max(1, nn[0]); L.ny = std::max(1, nn[1]); L.nz = std::max(1, nn[2]);
            if (int(L.hxfin.size()) != L.nx) L.hxfin.clear();
            if (int(L.hyfin.size()) != L.ny) L.hyfin.clear();
            if (int(L.hzfin.size()) != L.nz) L.hzfin.clear();
            L.props.clear();
            ch = true;
        }
        auto finI = [&](const char* lab, const char* tipText, std::vector<int>& v, int& N) {
            ui::label(lab, tipText);
            std::string t = intsToText(v);
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ui::inputString((std::string("##") + lab).c_str(), &t, ImGuiInputTextFlags_EnterReturnsTrue)) {
                v = textToNums<int>(t);
                if (!v.empty()) { int s = 0; for (int x : v) s += x; N = s; L.props.clear(); }
                ch = true;
            }
        };
        auto finD = [&](const char* lab, const char* tipText, std::vector<double>& v) {
            ui::label(lab, tipText);
            std::string t = dblsToText(v);
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ui::inputString((std::string("##") + lab).c_str(), &t, ImGuiInputTextFlags_EnterReturnsTrue)) { v = textToNums<double>(t); ch = true; }
        };
        finI("NXFIN", "Local cells in each host column (one value per host column, blank = even). Sets NX.", L.nxfin, L.nx);
        finI("NYFIN", "Local cells in each host row (blank = even). Sets NY.", L.nyfin, L.ny);
        finI("NZFIN", "Local layers in each host layer (blank = even). Sets NZ.", L.nzfin, L.nz);
        finD("HXFIN", "Relative widths of the local cells in I (NX values, blank = equal), e.g. log spacing towards a well.", L.hxfin);
        finD("HYFIN", "Relative widths in J (NY values).", L.hyfin);
        finD("HZFIN", "Relative thicknesses in K (NZ values).", L.hzfin);
        std::string why;
        if (!L.valid(&why)) { ui::label("Status"); ImGui::TextColored(col4(Pal::Warn), "%s", why.c_str()); }
        else {
            char b[96]; std::snprintf(b, sizeof b, "%d local cells replace %d host cells", L.size(), hn[0] * hn[1] * hn[2]);
            ui::propText("Cells", b);
        }
        ui::endProps();
    }
    if (ImGui::SmallButton(TRL("Logarithmic I/J spacing towards the centre"))) {
        auto logw = [](int n) {
            std::vector<double> w(static_cast<size_t>(n));
            for (int m = 0; m < n; ++m) { const double d = std::abs(m - 0.5 * (n - 1)); w[size_t(m)] = std::pow(1.6, d); }
            return w;
        };
        L.nxfin.clear(); L.nyfin.clear();
        L.hxfin = logw(L.nx); L.hyfin = logw(L.ny);
        ch = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("HXFIN / HYFIN: geometric (factor 1.6) relative widths, smallest at the centre of the LGR. "
                                                  "Widths are normalised inside each host cell, so host faces stay where they are."));

    // ---- local property overrides (CARFIN ... ENDFIN arrays, REFINE in PROPS / REGIONS / SOLUTION)
    if (ImGui::TreeNode(TRL("Local properties"))) {
        ui::hint("Arrays defined inside the refinement. Cells without a local value inherit the host cell. Written as "
                 "CARFIN ... PORO / PERMX ... ENDFIN (GRID) and REFINE ... ENDFIN (PROPS, REGIONS, SOLUTION).");
        std::string drop;
        for (auto& [name, v] : L.props) {
            double lo = 1e300, hi = -1e300; int nset = 0;
            for (double x : v) if (!std::isnan(x)) { lo = std::min(lo, x); hi = std::max(hi, x); ++nset; }
            ImGui::BulletText("%s: %d of %d cells set, %.4g - %.4g", name.c_str(), nset, int(v.size()), nset ? lo : 0.0, nset ? hi : 0.0);
            ImGui::SameLine(); ImGui::PushID(name.c_str()); if (ImGui::SmallButton(TRL("Clear"))) drop = name; ImGui::PopID();
        }
        if (!drop.empty()) { L.props.erase(drop); ch = true; }
        static const char* names[] = {"PORO", "PERMX", "PERMY", "PERMZ", "NTG", "MULTX", "MULTY", "MULTZ", "MULTPV", "ACTNUM",
                                      "SATNUM", "PVTNUM", "EQLNUM", "FIPNUM", "SWATINIT", "SWL", "SWCR", "SGCR", "SOWCR"};
        static int pn = 0, lay0 = 1, lay1 = 1; static double val = 1.0;
        lay1 = std::max(lay1, lay0);
        if (ui::beginProps("lgrprop")) {
            ui::propCombo("Array", &pn, names, IM_ARRAYSIZE(names));
            ui::propD("Value / factor", &val, "%.6g");
            ui::label("Local layers K1 K2"); int kk[2] = {lay0, lay1}; ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputInt2("##lk", kk)) { lay0 = std::max(1, std::min(L.nz, kk[0])); lay1 = std::max(lay0, std::min(L.nz, kk[1])); }
            ui::endProps();
        }
        lay0 = std::min(lay0, L.nz); lay1 = std::min(std::max(lay1, lay0), L.nz);
        auto apply = [&](bool mult) {
            if (!L.valid()) return;
            std::vector<double>& v = L.props[names[pn]];
            v.resize(size_t(L.size()), NAN);
            const std::vector<int> host = lgrHostCells(g, L);
            const std::string nmp = names[pn];
            const std::vector<double>* hv = nmp == "PORO" ? &g.poro : nmp == "PERMX" ? &g.permx : nmp == "PERMY" ? (g.permy.empty() ? &g.permx : &g.permy)
                                          : nmp == "PERMZ" ? (g.permz.empty() ? &g.permx : &g.permz) : nmp == "NTG" ? &g.ntg : nmp == "MULTPV" ? &g.multpv : nullptr;
            for (int k = lay0 - 1; k <= lay1 - 1; ++k)
                for (int j = 0; j < L.ny; ++j)
                    for (int i = 0; i < L.nx; ++i) {
                        const size_t m = size_t(L.index(i, j, k));
                        if (!mult) { v[m] = val; continue; }
                        double base = v[m];
                        if (std::isnan(base)) base = (hv && hv->size() == size_t(g.size())) ? (*hv)[size_t(host[m])] : 1.0;
                        v[m] = base * val;
                    }
            ch = true;
        };
        if (ImGui::Button(TRL("Set (EQUALS)"))) apply(false);
        ImGui::SameLine();
        if (ImGui::Button(TRL("Multiply (MULTIPLY)"))) apply(true);
        ImGui::TreePop();
    }
    ImGui::PopID();
    if (ch) { clampWellsToGrid(case_); markDirty(true); }
}

// ============================================================================ RUNSPEC
void App::editRunspec() {
    ui::header("Case");
    if (ui::beginProps("rs")) {
        ui::label("Title");
        if (ui::inputString("##title", &case_.title)) markDirty();
        static const char* units[] = {"METRIC (bar, m, sm3)", "FIELD (psia, ft, stb, Mscf)"};
        int u = case_.units == UnitSystem::Field ? 1 : 0;
        if (ui::propCombo("Unit system", &u, units, 2, "Values are NOT converted when switching; re-enter data in the new units.")) {
            case_.units = u ? UnitSystem::Field : UnitSystem::Metric; markDirty(true);
        }
        ui::endProps();
    }
    ui::header("Active phases");
    bool ch = false;
    ch |= ImGui::Checkbox(TRL("Oil"), &case_.phases.oil); ImGui::SameLine();
    ch |= ImGui::Checkbox(TRL("Water"), &case_.phases.water); ImGui::SameLine();
    ch |= ImGui::Checkbox(TRL("Gas"), &case_.phases.gas); ImGui::SameLine();
    ImGui::BeginDisabled(!(case_.phases.gas && case_.phases.oil));
    ch |= ImGui::Checkbox(TRL("Dissolved gas (DISGAS)"), &case_.phases.disgas); ImGui::SameLine();
    const bool vapWas = case_.phases.vapoil;
    ch |= ImGui::Checkbox(TRL("Vaporised oil (VAPOIL)"), &case_.phases.vapoil);
    ImGui::EndDisabled();
    if (case_.phases.vapoil && !(case_.phases.gas && case_.phases.oil)) case_.phases.vapoil = false;
    if (case_.phases.vapoil && !vapWas && case_.props.pvtg.size() < 2) {
        // seed PVTG from PVDG: Rv = 0 everywhere (a dry gas expressed as wet gas) - the user edits Rv afterwards
        auto& P = case_.props; P.pvtg.clear();
        for (auto& r : P.pvdg) { PVTGRecord g; g.pg = r.p; g.rows.push_back({0.0, r.B, r.mu}); P.pvtg.push_back(g); }
        if (!P.pvtg.empty()) { auto& L = P.pvtg.back(); L.rows.front().rv = 1e-4; L.rows.push_back({0.0, L.rows.front().B * 0.99, L.rows.front().mu * 0.98}); }
    }
    if (ch) markDirty();
    ui::hint("Black-oil formulation: water, oil and gas components; gas may dissolve in oil (live oil, PVTO) and oil may vaporise "
             "in gas (wet gas / gas condensate, PVTG). With DISGAS off, oil uses PVDO (dead oil); with VAPOIL off, gas uses PVDG "
             "(dry gas). Supported systems: oil-water, gas-oil, 3-phase, single-phase water or oil.");
    ui::header("Model summary");
    if (ui::beginProps("sum")) {
        char b[128];
        std::snprintf(b, sizeof(b), "%d x %d x %d = %d cells", case_.grid.nx, case_.grid.ny, case_.grid.nz, case_.grid.size());
        ui::propText("Grid", b);
        if (preview_) { std::snprintf(b, sizeof(b), "%d", preview_->numActive()); ui::propText("Active cells", b); }
        std::snprintf(b, sizeof(b), "%d", int(case_.wells.size())); ui::propText("Wells", b);
        double T = 0; for (auto& s : case_.schedule) T += s.length;
        std::snprintf(b, sizeof(b), "%d report steps, %.0f days (%.2f years)", int(case_.schedule.size()), T, T / 365.25);
        ui::propText("Schedule", b);
        ui::endProps();
    }
}

// ============================================================================ GRID
void App::editGrid() {
    bool geomCh = false, rockCh = false;
    const bool cp = case_.grid.cornerPoint();
    if (cp) {
        ui::header("Geometry (corner-point COORD / ZCORN)");
        if (ui::beginProps("cpgeo")) {
            char b[128];
            std::snprintf(b, sizeof(b), "%d x %d x %d  (%d cells)", case_.grid.nx, case_.grid.ny, case_.grid.nz, case_.grid.size());
            ui::propText("Dimensions", b);
            std::snprintf(b, sizeof(b), "%d pillars, %d corner depths", (case_.grid.nx + 1) * (case_.grid.ny + 1), int(case_.grid.zcorn.size()));
            ui::propText("Geometry", b);
            if (preview_) { std::snprintf(b, sizeof(b), "%d (faults / explicit)", preview_->numNNC); ui::propText("Non-neighbour connections", b); }
            ui::endProps();
        }
        ui::hint("Corner-point geometry (imported deck, converted grid or fault throws). Layer DZ and dips come from ZCORN; "
                 "rock properties per layer can still be regenerated below.");
        if (ImGui::Button(TRL("Discard corner-point geometry (rebuild Cartesian)"))) {
            try { generateGeometry(geom_, case_.grid); case_.grid.coord.clear(); case_.grid.zcorn.clear(); clampWellsToGrid(case_); markDirty(true); }
            catch (const std::exception& e) { gridError_ = e.what(); }
        }
    } else {
    ui::header("Geometry (Cartesian block-centred)");
    {   // C13: refined / variable-spacing grids (EGS designer, imported DX/DY arrays) would silently become uniform on edit
        const GridInput& g0 = case_.grid;
        auto span = [](const std::vector<double>& v, double& lo, double& hi) {
            lo = 1e300; hi = -1e300; for (double x : v) { lo = std::min(lo, x); hi = std::max(hi, x); }
            return !v.empty() && hi > lo * (1 + 1e-9);
        };
        double xl = 0, xh = 0, yl = 0, yh = 0;
        const bool vx = span(g0.dx, xl, xh), vy = span(g0.dy, yl, yh);
        if (vx || vy) {
            char b[400];
            std::snprintf(b, sizeof(b), "Variable cell sizes: DX %.3g - %.3g, DY %.3g - %.3g (refined grid). Editing the geometry "
                          "or the layer rock table below (or Build grid) rebuilds a uniform grid and removes the refinement - edit DX / DY in the Deck Editor instead.",
                          vx ? xl : g0.dx.empty() ? 0.0 : g0.dx[0], vx ? xh : g0.dx.empty() ? 0.0 : g0.dx[0],
                          vy ? yl : g0.dy.empty() ? 0.0 : g0.dy[0], vy ? yh : g0.dy.empty() ? 0.0 : g0.dy[0]);
            ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Warn)); ImGui::TextWrapped("%s", b); ImGui::PopStyleColor();
        }
    }
    if (ui::beginProps("geo")) {
        geomCh |= ui::propI("NX (cells in I)", &geom_.nx, 1, 2000);
        geomCh |= ui::propI("NY (cells in J)", &geom_.ny, 1, 2000);
        geomCh |= ui::propI("NZ (layers)", &geom_.nz, 1, 500);
        geomCh |= ui::propD("DX (cell length I)", &geom_.dx, "%.4g");
        geomCh |= ui::propD("DY (cell length J)", &geom_.dy, "%.4g");
        geomCh |= ui::propD("Top depth", &geom_.top, "%.6g", "Depth of the top face of layer 1 at the origin corner (TOPS).");
        geomCh |= ui::propSliderD("Dip along I (deg)", &geom_.dipX, -30, 30, "%.1f");
        geomCh |= ui::propSliderD("Dip along J (deg)", &geom_.dipY, -30, 30, "%.1f");
        if (ui::propBool("Geometric transmissibility", &case_.grid.newtran,
                         "NEWTRAN: neighbouring columns at different depths connect through their actual face overlap "
                         "(non-neighbour connections where layers are offset). This is what OPM Flow always does.\n"
                         "Off = OLDTRAN (ECLIPSE default for block-centred grids): cell (i,j,k) connects to (i+1,j,k) "
                         "whatever the depth offset. Matters for dipping grids when the column offset approaches the layer thickness."))
            markDirty(true);
        ui::endProps();
    }
    if (ImGui::SmallButton(TRL("Convert to corner-point"))) { convertToCornerPoint(case_.grid); markDirty(true); notify("Converted to corner-point geometry (COORD/ZCORN)"); }
    ImGui::SameLine(); ImGui::TextDisabled("%s", TR("enables fault throws below"));
    }
    // ---- dual porosity / dual permeability (ECLIPSE DUALPORO / DUALPERM)
    ui::header("Dual porosity (naturally fractured reservoirs)");
    {
        GridInput& g = case_.grid;
        static double fPoro = 0.005, fPerm = 500, bL = 10, bH = 10; static bool fRegion = true;
        if (!g.dual) {
            ui::hint("Single porosity. Converting doubles NZ: layers NZ+1..2NZ become the fracture system of the same blocks "
                     "(ECLIPSE layout). Each matrix block couples to its fracture cell with T = CDARCY V sigma K.");
            if (ui::beginProps("dpnew")) {
                ui::propD("Fracture porosity", &fPoro, "%.4g", "Fracture pore volume per bulk volume (typically 0.1 - 1 %).");
                ui::propD("Fracture permeability (mD)", &fPerm, "%.5g", "Intrinsic (bulk) fracture-system permeability; NODPPM is set.");
                ui::propD("Matrix block size L (lateral)", &bL, "%.4g", "Fracture spacing in I and J (length units of the case).");
                ui::propD("Matrix block height H", &bH, "%.4g", "Fracture spacing in the vertical direction.");
                ui::propBool("Straight-line kr in the fractures", &fRegion, "Adds a SATNUM region for the fracture cells: krw = Sw, kro = So, Pc = 0.");
                ui::endProps();
            }
            const bool hasLgr = !g.lgrs.empty();
            if (hasLgr) ui::hint("Remove the local grid refinements first: LGRs inside dual porosity models are not modelled by ResSim.");
            ImGui::BeginDisabled(hasLgr);
            if (ImGui::Button(TRL("Convert to dual porosity (DUALPORO)"))) {
                makeDualPorosity(case_, 1, fPoro, fPerm, std::max(1e-3, bL), std::max(1e-3, bH), fRegion);
                syncBuildersFromCase(); markDirty(true);
                notify("Dual porosity: NZ doubled, fracture layers added (SIGMA from the block size)");
            }
            ImGui::SameLine();
            if (ImGui::Button(TRL("... dual permeability (DUALPERM)"))) {
                makeDualPorosity(case_, 2, fPoro, fPerm, std::max(1e-3, bL), std::max(1e-3, bH), fRegion);
                syncBuildersFromCase(); markDirty(true);
                notify("Dual porosity / dual permeability model created");
            }
            ImGui::EndDisabled();
        } else {
            bool ch = false;
            if (ui::beginProps("dp")) {
                int mode = g.dual - 1;
                static const char* modes[] = {"Dual porosity (DUALPORO)", "Dual porosity + dual permeability (DUALPERM)"};
                if (ui::propCombo("Model", &mode, modes, 2, "DUALPORO: the matrix blocks only exchange fluid with their fracture cell, wells "
                                  "connect to the fractures.\nDUALPERM: matrix blocks also flow to each other and wells connect to both media.")) {
                    g.dual = mode + 1; ch = true;
                }
                char b[96]; std::snprintf(b, sizeof b, "%d matrix + %d fracture", g.nz / 2, g.nz / 2);
                ui::propText("Layers", b);
                int sm = g.ltosigma ? 1 : 0;
                static const char* sms[] = {"SIGMA (shape factor)", "From block sizes (LTOSIGMA, LX/LY/LZ)"};
                if (ui::propCombo("Coupling", &sm, sms, 2)) { g.ltosigma = sm == 1; ch = true; }
                if (!g.ltosigma) {
                    ch |= ui::propD("SIGMA (1/length^2)", &g.sigma, "%.5g",
                                    "Matrix-fracture interface area per unit volume. Kazemi: sigma = 4 (1/lx^2 + 1/ly^2 + 1/lz^2); "
                                    "Warren-Root / Gilman-Kazemi variants differ by a constant factor.");
                    if (!g.sigmav.empty()) ui::propText("SIGMAV", "per-cell values from the deck (override SIGMA)");
                } else {
                    double lx = g.lx.empty() ? 10.0 : g.lx[0], ly = g.ly.empty() ? 10.0 : g.ly[0], lz = g.lz.empty() ? 10.0 : g.lz[0];
                    bool lch = ui::propD("LX (block length I)", &lx, "%.4g");
                    lch |= ui::propD("LY (block length J)", &ly, "%.4g");
                    lch |= ui::propD("LZ (block height)", &lz, "%.4g");
                    if (lch) { g.lx.assign(size_t(g.size()), lx); g.ly.assign(size_t(g.size()), ly); g.lz.assign(size_t(g.size()), lz); ch = true; }
                    ch |= ui::propD("Factor fx", &g.ltoF[0], "%.3g"); ch |= ui::propD("Factor fy", &g.ltoF[1], "%.3g"); ch |= ui::propD("Factor fz", &g.ltoF[2], "%.3g");
                }
                ch |= ui::propBool("Intrinsic fracture permeability (NODPPM)", &g.nodppm,
                                   "Off (ECLIPSE default): the fracture PERMX/Y/Z are multiplied by the fracture porosity.");
                if (preview_) { std::snprintf(b, sizeof b, "%d", preview_->numMF); ui::propText("Matrix-fracture connections", b); }
                // C14 MINC: nested matrix shells (Pruess & Narasimhan) - transient matrix-fracture exchange
                if (g.dual == 1) {
                    bool on = g.minc > 1;
                    if (ui::propBool("MINC nested shells", &on,
                                     "Multiple interacting continua (Pruess & Narasimhan 1985): each matrix block is split into nested "
                                     "shells, the outermost touching the fractures. Resolves the transient (early-time) matrix drainage / "
                                     "imbibition that the pseudo-steady dual porosity transfer smears out. ResSim engine (black oil) only; "
                                     "ECLIPSE / OPM runs use the dual porosity model.")) { g.minc = on ? 6 : 0; ch = true; }
                    if (g.minc > 1) {
                        int nsh = g.minc;
                        if (ui::propI("Shells", &nsh, 2, 50)) { g.minc = nsh; if (!g.mincVol.empty()) g.mincVol.clear(); ch = true; }
                        int dims = g.mincDims - 1;
                        static const char* dm[] = {"1 set (horizontal fractures, Z)", "2 sets (vertical, X + Y)", "3 sets (X + Y + Z)"};
                        if (ui::propCombo("Fracture sets", &dims, dm, 3, "Used when LX / LY / LZ are not given: the block size follows from "
                                          "SIGMA (L = sqrt(4 sets / sigma)). With LX / LY / LZ the non-zero sizes define the sets.")) { g.mincDims = dims + 1; ch = true; }
                        std::string fr;
                        if (g.mincVol.empty()) fr = "geometric 1 : 2 : 4 ... (outer shells thin)";
                        else for (double f : g.mincVol) { char t[24]; std::snprintf(t, sizeof t, "%s%.3g", fr.empty() ? "" : " ", f); fr += t; }
                        ui::propText("Volume fractions", fr.c_str());
                        if (preview_ && preview_->mincN > 1) { std::snprintf(b, sizeof b, "%d", preview_->numMincCells()); ui::propText("Shell cells added", b); }
                    }
                }
                ch |= ui::propBool("Gravity drainage (GRAVDR)", &g.gravdr,
                                   "Adds the gravity head of segregated fluid columns to the matrix-fracture flow (ECLIPSE GRAVDR): "
                                   "water imbibes from below and gas drains oil from above when the fracture column is higher than "
                                   "the matrix one. dPhi += drho g Z (H_matrix - H_fracture), H = normalised saturation.");
                if (g.gravdr) {
                    bool uni = true; for (double v : g.dzmtrxv) if (v != g.dzmtrxv.front()) uni = false;
                    if (uni) {
                        double z = g.dzmtrxv.empty() ? 0.0 : g.dzmtrxv.front();
                        if (ui::propD("Matrix block height (DZMTRX)", &z, "%.4g", "Vertical size of the matrix blocks; 0 = the cell thickness.")) {
                            if (z > 0) g.dzmtrxv.assign(size_t(g.size()), z); else g.dzmtrxv.clear();
                            ch = true;
                        }
                    } else ui::propText("Matrix block height", "DZMTRXV per cell (from the deck)");
                }
                if (!g.dpnum.empty()) {
                    int nsp = 0; for (int v : g.dpnum) nsp += v == 0;
                    std::snprintf(b, sizeof b, "%d single porosity blocks (DPNUM = 0)", nsp);
                    ui::propText("DPNUM", b);
                }
                ui::endProps();
            }
            if (ch) markDirty(true);
            ui::hint("Rock properties of the fracture system are the lower half of the layer table. OPM Flow does not support "
                     "dual porosity - run it with ResSim or ECLIPSE.");
            if (ImGui::SmallButton(TRL("Back to single porosity (drop the fracture layers)"))) {
                makeSinglePorosity(case_); syncBuildersFromCase(); markDirty(true);
            }
        }
    }

    editLgrs();
    editHydraulicFractures();

    const int nz = geom_.nz;
    auto fit = [nz](std::vector<double>& v, double def) { if (int(v.size()) != nz) { const double last = v.empty() ? def : v.back(); v.resize(nz, last); } };
    fit(geom_.dz, 10); fit(rock_.layerPerm, 100); fit(rock_.layerPoro, 0.2); fit(rock_.layerKvKh, 0.1); fit(rock_.layerNtg, 1.0);

    ui::header("Layers");
    {   // spreadsheet grid: one row per layer (paste whole columns from Excel / a log table)
        ui::GridAccess ga;
        ga.count = [nz] { return nz; };
        ga.cols = 5;
        ga.num = [&](int k, int c) -> double* {
            if (cp && c == 0) return nullptr;   // corner-point: thickness comes from ZCORN
            double* cols[5] = {&geom_.dz[size_t(k)], &rock_.layerPerm[size_t(k)], &rock_.layerPoro[size_t(k)], &rock_.layerKvKh[size_t(k)], &rock_.layerNtg[size_t(k)]};
            return cols[c];
        };
        ga.getText = [&](int k, int) { char b[32]; std::snprintf(b, sizeof b, "%.4g (ZCORN)", geom_.dz[size_t(k)]); return std::string(b); };
        ga.rowLabel = [](int k) { return "K" + std::to_string(k + 1); };
        std::vector<double> before = geom_.dz;
        if (ui::gridTable("layers", ga, {"DZ", "Perm (mD)", "Porosity", "Kv/Kh", "NTG"}, std::min(ImGui::GetFontSize() * 18, ImGui::GetFrameHeight() * (nz + 1.6f) + 4), nullptr, "%.4g", false)) {
            if (before != geom_.dz) geomCh = true;
            rockCh = true;
        }
    }
    ui::header("Heterogeneity (geostatistics)");
    if (ui::beginProps("het", 0.56f)) {
        rockCh |= ui::propSliderD("Perm variability  \xCF\x83(ln k)", &rock_.sigmaLn, 0.0, 2.5, "%.2f",
                                  "Standard deviation of ln(k). 0 = homogeneous layers; 1 = strongly heterogeneous. Mean permeability per layer is preserved.");
        rockCh |= ui::propSliderD("Correlation length I (cells)", &rock_.corrX, 0.0, 30.0, "%.1f");
        rockCh |= ui::propSliderD("Correlation length J (cells)", &rock_.corrY, 0.0, 30.0, "%.1f");
        rockCh |= ui::propSliderD("Poro-perm exponent", &rock_.poroPermExp, 0.0, 0.5, "%.2f", "poro = layer poro * (k / k_layer)^exponent");
        ui::label("Random seed");
        int seed = int(rock_.seed);
        ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - ImGui::GetFontSize() * 5);
        if (ImGui::InputInt("##seed", &seed)) { rock_.seed = uint32_t(seed); rockCh = true; }
        ImGui::SameLine();
        if (ImGui::Button(TRL("New seed"))) { rock_.seed = std::random_device{}(); rockCh = true; }
        ui::endProps();
    }
    // ---- C28 kriging / SGS conditioned to the wells (CMG Builder / Petrel-style property modelling, GSLIB conventions)
    if (ImGui::CollapsingHeader(TRL("Kriging / sequential Gaussian simulation (conditioned to wells)"))) {
        static geostat::ApplyInput gi;
        static int vm = 0;
        static double lagW = 0; static int nLag = 15, dirE = 0;
        static geostat::ExpVario ev; static int evProp = -1;
        static std::string gmsg;
        ui::hint("Fills PORO, PERMX (PERMY = PERMX, PERMZ = kv/kh PERMX; modelled as log10 k) or NTG of every active cell from hard "
                 "data = the current values at the completed cells of all wells. Kriging gives the smooth best estimate; SGS one "
                 "equiprobable realisation with the variogram's variability (use different seeds in a sensitivity study). Check the "
                 "model variogram against the experimental one of the current property below. Applying switches Auto-apply off so "
                 "the layer generator does not overwrite the result.");
        if (ui::beginProps("gsp", 0.56f)) {
            static const char* kGP[] = {"PORO", "PERMX (log10)", "NTG"};
            static const char* kGM[] = {"Ordinary kriging", "Sequential Gaussian simulation"};
            static const char* kVM[] = {"Spherical", "Exponential", "Gaussian"};
            ui::propCombo("Property", &gi.property, kGP, 3);
            ui::propCombo("Method", &gi.method, kGM, 2);
            if (ui::propCombo("Variogram", &vm, kVM, 3)) gi.opt.vg.model = geostat::VarioModel(vm);
            ui::propD("Nugget (fraction of sill)", &gi.opt.vg.nugget, "%.3g", "Nugget / total sill; the structured sill is 1 - nugget (normal scores)");
            gi.opt.vg.nugget = std::clamp(gi.opt.vg.nugget, 0.0, 0.99); gi.opt.vg.sill = 1.0 - gi.opt.vg.nugget;
            ui::propD("Range major (length)", &gi.opt.vg.rangeMax);
            ui::propD("Range minor (length)", &gi.opt.vg.rangeMin);
            ui::propD("Range vertical (length)", &gi.opt.vg.rangeVert);
            ui::propD("Azimuth of major (deg)", &gi.opt.vg.azimuthDeg, "%.4g", "Clockwise from north (+y)");
            ui::propI("Max conditioning points", &gi.opt.maxNbr, 4, 64);
            if (gi.property == 1) ui::propD("kv/kh", &gi.kvkh, "%.3g");
            if (gi.method == 1) { int sd = int(gi.opt.seed); if (ui::propI("Seed", &sd, 1, 1 << 30)) gi.opt.seed = unsigned(sd); }
            ui::endProps();
        }
        gi.opt.vg.rangeMax = std::max(gi.opt.vg.rangeMax, 1e-6); gi.opt.vg.rangeMin = std::max(gi.opt.vg.rangeMin, 1e-6);
        gi.opt.vg.rangeVert = std::max(gi.opt.vg.rangeVert, 1e-6);
        if (gi.property == 0) { gi.minV = 0.0; gi.maxV = 0.6; } else if (gi.property == 2) { gi.minV = 0.0; gi.maxV = 1.0; } else { gi.minV = 1e-4; gi.maxV = 1e5; }
        if (ui::accentButton("Model the property from the wells")) {
            const geostat::ApplyReport rep = geostat::applyToCase(case_, gi);
            char b[256];
            if (rep.ok) {
                std::snprintf(b, sizeof b, "%d cells from %d data points: mean %.4g, sd %.4g, range [%.4g, %.4g]", rep.cells, rep.dataPoints, rep.mean, rep.stdDev, rep.minV, rep.maxV);
                autoBuild_ = false; evProp = -1;
                {   // C40: keep the recipe with the case, so Studies can re-simulate it (variogram range / realisation parameters)
                    PropModelInput m;
                    m.property = gi.property; m.method = gi.method; m.model = int(gi.opt.vg.model);
                    m.nugget = gi.opt.vg.nugget; m.sill = gi.opt.vg.sill; m.rangeMax = gi.opt.vg.rangeMax; m.rangeMin = gi.opt.vg.rangeMin;
                    m.rangeVert = gi.opt.vg.rangeVert; m.azimuthDeg = gi.opt.vg.azimuthDeg; m.seed = gi.opt.seed; m.maxNbr = gi.opt.maxNbr;
                    m.fromWells = gi.fromWells; m.minV = gi.minV; m.maxV = gi.maxV; m.kvkh = gi.kvkh; m.source = "Model geostatistics";
                    for (const auto& d : gi.data) m.data.push_back({d.x, d.y, d.z, d.v});
                    auto& PM = case_.propModels;
                    PM.erase(std::remove_if(PM.begin(), PM.end(), [&](const PropModelInput& q) { return q.property == m.property && q.source == m.source; }), PM.end());
                    PM.push_back(m);
                }
                markDirty(true);
            } else std::snprintf(b, sizeof b, "%s", rep.message.c_str());
            gmsg = b;
        }
        if (!gmsg.empty()) { ImGui::SameLine(); ImGui::TextDisabled("%s", gmsg.c_str()); }
        // experimental variogram of the current property (normal scores are not used here: raw / log10 values, sill = variance)
        if (preview_ && preview_->numCells() == case_.grid.size()) {
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
            ImGui::InputInt("Lags", &nLag); nLag = std::clamp(nLag, 3, 60);
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
            ImGui::InputDouble("Lag width (0 = auto)", &lagW, 0, 0, "%.4g");
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9);
            static const char* kDir[] = {"Horizontal (omni)", "Vertical", "Along azimuth"};
            ImGui::Combo("##evdir", &dirE, kDir, 3);
            ImGui::SameLine();
            const std::vector<double>* P = gi.property == 0 ? &case_.grid.poro : gi.property == 1 ? &case_.grid.permx : &case_.grid.ntg;
            if (ImGui::Button(TRL("Experimental variogram")) && P->size() == size_t(case_.grid.size())) {
                const Grid& g = *preview_;
                std::vector<geostat::Point> pts;
                double ext = 0;
                for (int a = 0; a < g.numActive(); ++a) {
                    const int q = g.globalOf[size_t(a)];
                    const double v = (*P)[size_t(q)];
                    if (gi.property == 1 && v <= 0) continue;
                    pts.push_back({g.center[size_t(q)].x, g.center[size_t(q)].y, g.center[size_t(q)].z, gi.property == 1 ? std::log10(v) : v});
                }
                double xmn = 1e300, xmx = -1e300, ymn = 1e300, ymx = -1e300, zmn = 1e300, zmx = -1e300;
                for (const auto& p : pts) { xmn = std::min(xmn, p.x); xmx = std::max(xmx, p.x); ymn = std::min(ymn, p.y); ymx = std::max(ymx, p.y); zmn = std::min(zmn, p.z); zmx = std::max(zmx, p.z); }
                ext = dirE == 1 ? zmx - zmn : std::hypot(xmx - xmn, ymx - ymn);
                const double lw = lagW > 0 ? lagW : std::max(ext, 1e-6) / (2.0 * nLag);
                ev = geostat::experimental(pts, lw, nLag, dirE, gi.opt.vg.azimuthDeg);
                double m = 0, v2 = 0;
                for (const auto& p : pts) m += p.v;
                m /= std::max<size_t>(1, pts.size());
                for (const auto& p : pts) v2 += (p.v - m) * (p.v - m);
                v2 /= std::max<size_t>(1, pts.size());
                for (auto& gm : ev.gamma) gm = v2 > 0 ? gm / v2 : gm;   // standardised: comparable with the normal-score model
                evProp = gi.property;
            }
            if (evProp >= 0 && !ev.lag.empty() && ImPlot::BeginPlot("Variogram (standardised by the variance)", ImVec2(-1, ImGui::GetFontSize() * 14))) {
                recPanel("Variogram (standardised by the variance)", "lag distance", "gamma / variance", "");
                ImPlot::SetupAxes("lag distance", "gamma / variance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlotSpec sp; sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 4; sp.LineColor = col4(Pal::Accent);
                sp.MarkerFillColor = col4(Pal::Accent); sp.MarkerLineColor = col4(Pal::Accent);
                std::vector<double> lx, ly;
                for (size_t i = 0; i < ev.lag.size(); ++i) if (ev.pairs[i] > 0) { lx.push_back(ev.lag[i]); ly.push_back(ev.gamma[i]); }
                if (fig::recorder().recording()) recSeries("experimental", lx.data(), ly.data(), int(lx.size()), Pal::Accent, 2.0f, true);
                ImPlot::PlotScatter("experimental", lx.data(), ly.data(), int(lx.size()), sp);
                std::vector<double> mx, my;
                const double hMax = ev.lag.back() * 1.05;
                const double az = gi.opt.vg.azimuthDeg * M_PI / 180.0;
                for (int i = 0; i <= 100; ++i) {
                    const double h = hMax * i / 100.0;
                    double dx = h, dy = 0, dz = 0;
                    if (dirE == 1) { dx = 0; dz = h; } else if (dirE == 2) { dx = h * std::sin(az); dy = h * std::cos(az); }
                    mx.push_back(h); my.push_back(geostat::gamma(gi.opt.vg, dx, dy, dz));
                }
                ImPlotSpec sm; sm.LineColor = col4(Pal::Warn); sm.LineWeight = 2.0f * g_plotLineScale;
                if (fig::recorder().recording()) recSeries("model", mx.data(), my.data(), int(mx.size()), Pal::Warn, 2.0f, false);
                ImPlot::PlotLine("model", mx.data(), my.data(), int(mx.size()), sm);
                recEnd();
                ImPlot::EndPlot();
            }
            if (evProp >= 0 && !ev.lag.empty()) {   // C40: weighted least-squares fit of the current model type to the points
                if (ImGui::Button(TRL(ICON_WAND_SPARKLES " Fit the model to the points"))) {
                    const auto f = geo::fitVariogram(ev, gi.opt.vg.model);
                    if (f.ok) {
                        gi.opt.vg.nugget = f.nugget; gi.opt.vg.sill = f.sill;
                        if (dirE == 1) gi.opt.vg.rangeVert = f.range; else if (dirE == 0) gi.opt.vg.rangeMax = gi.opt.vg.rangeMin = f.range; else gi.opt.vg.rangeMax = f.range;
                        char b[160]; std::snprintf(b, sizeof b, "fitted: nugget %.3f, sill %.3f, range %.4g (rms %.3g)", f.nugget, f.sill, f.range, f.rmse);
                        gmsg = b;
                    } else gmsg = "variogram fit failed (flat or too few lags)";
                }
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Cressie-weighted least squares; sets nugget, sill and the range of the selected direction."));
            }
        }
    }
    // ---- C41 facies by truncated Gaussian simulation (stored as a property model: Studies re-simulates it)
    if (ImGui::CollapsingHeader(TRL("Facies (truncated Gaussian simulation)"))) {
        static PropModelInput fm = [] {
            PropModelInput m; m.property = 3; m.method = 1; m.rangeMax = m.rangeMin = 1000; m.rangeVert = 20;
            m.facProp = {0.6, 0.4}; m.facPoro = {0.25, 0.08}; m.facPerm = {500, 5}; m.facSatnum = false; m.fromWells = false; m.kvkh = 0.1;
            m.source = "Model facies"; return m;
        }();
        static std::string fmsg;
        ui::hint("Facies by truncated Gaussian simulation: one Gaussian field with the variogram below, cut at the thresholds of "
                 "the cumulative target proportions (ordered facies: facies 1 touches 2, 2 touches 3, ...). Each cell gets its "
                 "facies' porosity and permeability (PERMY = PERMX, PERMZ = kv/kh PERMX; a negative value keeps the current one) "
                 "and optionally SATNUM = facies number, so every facies uses its own SWOF / SGOF table. Conditioning: with "
                 "'Condition on SATNUM at wells' the current SATNUM of the completed cells is honoured. The recipe is kept with "
                 "the case: Studies' geomodel parameters (variogram range, realisation) re-simulate it.");
        int K = int(fm.facProp.size());
        if (ui::beginProps("fac", 0.56f)) {
            static const char* kVM[] = {"Spherical", "Exponential", "Gaussian"};
            if (ui::propI("Number of facies", &K, 2, 6)) { fm.facProp.resize(size_t(K), 0.1); fm.facPoro.resize(size_t(K), 0.15); fm.facPerm.resize(size_t(K), 50); }
            ui::propCombo("Variogram##fac", &fm.model, kVM, 3);
            ui::propD("Range major (length)##fac", &fm.rangeMax);
            ui::propD("Range minor (length)##fac", &fm.rangeMin);
            ui::propD("Range vertical (length)##fac", &fm.rangeVert);
            ui::propD("Azimuth of major (deg)##fac", &fm.azimuthDeg, "%.4g");
            ui::propD("kv/kh##fac", &fm.kvkh, "%.3g");
            ui::propBool("Set SATNUM = facies", &fm.facSatnum);
            if (fm.facSatnum) ui::propBool("Condition on SATNUM at wells", &fm.fromWells);
            int sd = int(fm.seed); if (ui::propI("Seed##fac", &sd, 1, 1 << 30)) fm.seed = unsigned(sd);
            ui::endProps();
        }
        fm.rangeMax = std::max(fm.rangeMax, 1e-6); fm.rangeMin = std::max(fm.rangeMin, 1e-6); fm.rangeVert = std::max(fm.rangeVert, 1e-6);
        if (ImGui::BeginTable("##factab", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_SizingStretchSame)) {
            ImGui::TableSetupColumn("Facies"); ImGui::TableSetupColumn("Proportion"); ImGui::TableSetupColumn("PORO"); ImGui::TableSetupColumn(case_.units == UnitSystem::Field ? "PERMX (mD)" : "PERMX (mD)");
            ImGui::TableHeadersRow();
            for (int k = 0; k < K; ++k) {
                ImGui::PushID(k);
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("%d", k + 1);
                ImGui::TableSetColumnIndex(1); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##p", &fm.facProp[size_t(k)], 0, 0, "%.3f");
                ImGui::TableSetColumnIndex(2); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##phi", &fm.facPoro[size_t(k)], 0, 0, "%.3f");
                ImGui::TableSetColumnIndex(3); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##k", &fm.facPerm[size_t(k)], 0, 0, "%.4g");
                ImGui::PopID();
            }
            ImGui::EndTable();
        }
        if (ui::accentButton("Simulate facies")) {
            const geostat::ApplyReport rep = geostat::applyFaciesModel(case_, fm);
            fmsg = rep.message;
            if (rep.ok) {
                autoBuild_ = false;
                auto& PM = case_.propModels;
                PM.erase(std::remove_if(PM.begin(), PM.end(), [&](const PropModelInput& q) { return q.property == 3 && q.source == fm.source; }), PM.end());
                PM.push_back(fm);
                markDirty(true);
            }
        }
        if (!fmsg.empty()) { ImGui::PushTextWrapPos(0); ImGui::TextDisabled("%s", fmsg.c_str()); ImGui::PopTextWrapPos(); }
    }
    ImGui::Spacing();
    ImGui::Checkbox(TRL("Auto-apply changes"), &autoBuild_);
    ImGui::SameLine();
    const bool apply = ui::accentButton(cp ? "Regenerate rock properties" : "Build grid", ImVec2(-FLT_MIN, 0));
    if (apply || (autoBuild_ && (geomCh || rockCh))) {
        try {
            if (!cp) generateGeometry(geom_, case_.grid);
            generateRock(rock_, case_.grid);
            clampWellsToGrid(case_);
            markDirty(true);
        } catch (const std::exception& e) { gridError_ = e.what(); }
    }
    if (!gridError_.empty()) ImGui::TextColored(col4(Pal::Err), "%s", gridError_.c_str());

    if (preview_ && preview_->numCells() == case_.grid.size() && case_.grid.permx.size() == size_t(case_.grid.size())) {
        ui::header("Statistics");
        const Grid& g = *preview_;
        double pv = 0, bulk = 0, kavg = 0, kmin = 1e300, kmax = 0;
        const GridInput& gi = case_.grid;
        std::vector<double> logk;
        int nExtra = 0;
        logk.reserve(g.numActive());
        for (int a = 0; a < g.numActive(); ++a) {
            const int q = g.globalOf[a];
            pv += g.poreVolRef[a];
            bulk += q < int(g.bulkVol.size()) ? g.bulkVol[q] : 0.0;
            if (q >= int(gi.permx.size()) || q >= g.numCells()) { ++nExtra; continue; }   // LGR / MINC / hydraulic fracture cells
            kavg += gi.permx[q]; kmin = std::min(kmin, gi.permx[q]); kmax = std::max(kmax, gi.permx[q]);
            logk.push_back(std::log10(std::max(gi.permx[q], 1e-6)));
        }
        if (!logk.empty()) kavg /= double(logk.size()); else kmin = 0;
        const Units U = Units::of(case_.units);
        if (ui::beginProps("stats")) {
            char b[96];
            if (nExtra) std::snprintf(b, sizeof(b), "%d of %d + %d local / fracture", g.numActive() - nExtra, g.numCells(), nExtra);
            else std::snprintf(b, sizeof(b), "%d of %d", g.numActive(), g.numCells());
            ui::propText("Active cells", b);
            std::snprintf(b, sizeof(b), "%.4g %s3", bulk, U.len); ui::propText("Bulk volume", b);
            std::snprintf(b, sizeof(b), "%.4g %s", pv, case_.units == UnitSystem::Field ? "rb" : "rm3"); ui::propText("Pore volume", b);
            std::snprintf(b, sizeof(b), "%.1f  (min %.2f, max %.1f) mD", kavg, kmin, kmax); ui::propText("Permeability X", b);
            std::snprintf(b, sizeof(b), "%d", int(g.conns.size())); ui::propText("Connections", b);
            std::snprintf(b, sizeof(b), "%.4g to %.4g %s", g.bbMin.z, g.bbMax.z, U.len); ui::propText("Depth range", b);
            ui::endProps();
        }
        if (ImPlot::BeginPlot("Permeability distribution", ImVec2(-1, ImGui::GetFontSize() * 11))) {
            ImPlot::SetupAxes("log10(PERMX, mD)", "cells", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlotSpec sp; sp.FillColor = col4(Pal::Accent); sp.FillAlpha = 0.75f; sp.LineColor = col4(Pal::Accent);
            if (!logk.empty()) ImPlot::PlotHistogram("##h", logk.data(), int(logk.size()), 30, 1.0, ImPlotRange(), sp);
            ImPlot::EndPlot();
        }
    }

    // ---- faults (FAULTS + MULTFLT) and fault-throw tool
    GridInput& gi = case_.grid;
    ui::header("Faults (FAULTS / MULTFLT)");
    {
        static int pos = 0; static double thr = 20; static double mult = 1.0; static bool alongJ = false; static char name[32] = "F1";
        if (pos <= 0) pos = std::max(1, gi.nx / 2);
        if (ImGui::TreeNodeEx(TRL("Add a fault with throw"), ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ui::beginProps("fthrow")) {
                static const char* dirs[] = {"Constant I (fault plane along J)", "Constant J (fault plane along I)"};
                int d = alongJ ? 1 : 0; if (ui::propCombo("Orientation", &d, dirs, 2)) alongJ = d == 1;
                ui::propI(alongJ ? "Between rows J and J+1" : "Between columns I and I+1", &pos, 1, (alongJ ? gi.ny : gi.nx) - 1);
                ui::propD("Throw (down on the far side)", &thr, "%.4g", "Vertical displacement of all cells beyond the fault plane. 0 = vertical barrier without offset.");
                ui::propSliderD("Transmissibility multiplier", &mult, 0.0, 1.0, "%.3f", "MULTFLT: 0 = sealing, 1 = open");
                ui::label("Name"); ImGui::InputText("##fn", name, sizeof(name));
                ui::endProps();
            }
            const int nmax = alongJ ? gi.ny : gi.nx;
            ImGui::BeginDisabled(nmax < 2);
            if (ui::accentButton("Add fault")) {
                try {
                    addFaultThrow(gi, std::min(pos, nmax - 1), thr, mult, name, alongJ);
                    syncBuildersFromCase();
                    markDirty(true);
                    notify(std::string("Fault ") + name + " added (corner-point grid)");
                    std::snprintf(name, sizeof(name), "F%d", int(gi.faults.size()) + 1);
                } catch (const std::exception& e) { gridError_ = e.what(); }
            }
            ImGui::EndDisabled();
            ImGui::TreePop();
        }
    }
    if (gi.faults.empty()) ui::hint("No faults defined.");
    else {   // C15 spreadsheet table: names and MULTFLT (paste a column of multipliers from Excel)
        ui::GridAccess ga;
        ga.count = [&] { return int(gi.faults.size()); };
        ga.cols = 3;
        ga.getText = [&](int r, int c) -> std::string {
            const FaultInput& F = gi.faults[size_t(r)];
            if (c == 0) return F.name;
            if (c == 2) { int cells = 0; for (auto& fc : F.faces) cells += (fc.i2 - fc.i1 + 1) * (fc.j2 - fc.j1 + 1) * (fc.k2 - fc.k1 + 1);
                          return std::to_string(F.faces.size()) + " (" + std::to_string(cells) + " cells)"; }
            return std::string();
        };
        ga.setText = [&](int r, int c, const std::string& v) { if (c == 0 && !v.empty()) gi.faults[size_t(r)].name = v; };
        ga.numGet = [&](int r, int c, double& v) -> bool { if (c != 1) return false; v = gi.faults[size_t(r)].mult; return true; };
        ga.numSet = [&](int r, int c, double v) { if (c == 1) gi.faults[size_t(r)].mult = std::max(0.0, v); };
        ga.erase = [&](int r0, int r1) { r1 = std::min(r1, int(gi.faults.size()) - 1); if (r0 <= r1) gi.faults.erase(gi.faults.begin() + r0, gi.faults.begin() + r1 + 1); };
        const float h = std::min(ImGui::GetFontSize() * 14, ImGui::GetFrameHeight() * (float(gi.faults.size()) + 2.6f) + 4);
        if (ui::gridTable("faults", ga, {"Name", "MULTFLT", "Faces (cells)"}, h, nullptr, "%.4g", true)) markDirty(true);
        ui::hint("Removing a fault keeps the geometric throw (ZCORN) but removes its MULTFLT barrier.");
    }
    editReservoirGeomech(); // C25: FE geomechanics of black-oil / compositional runs
    editFaultMechanics();   // C24: stress field, fault mechanics, seismicity, G1-c coupling (run with the model)

    // ---- explicit non-neighbour connections
    ui::header("Explicit non-neighbour connections (NNC)");
    {
        ui::GridAccess ga;
        ga.count = [&] { return int(gi.nnc.size()); };
        ga.cols = 7;
        ga.numGet = [&](int r, int c, double& v) -> bool {
            const NNCInput& q = gi.nnc[size_t(r)];
            const int* idx[6] = {&q.i1, &q.j1, &q.k1, &q.i2, &q.j2, &q.k2};
            if (c < 6) v = *idx[c] + 1; else v = q.trans;
            return true;
        };
        ga.numSet = [&](int r, int c, double v) {
            NNCInput& q = gi.nnc[size_t(r)];
            int* idx[6] = {&q.i1, &q.j1, &q.k1, &q.i2, &q.j2, &q.k2};
            const int lim[6] = {gi.nx, gi.ny, gi.nz, gi.nx, gi.ny, gi.nz};
            if (c < 6) *idx[c] = std::clamp(int(std::lround(v)), 1, lim[c]) - 1; else q.trans = std::max(0.0, v);
        };
        ga.insert = [&](int at, int n) {
            at = std::clamp(at, 0, int(gi.nnc.size()));
            NNCInput q; q.i2 = std::min(1, gi.nx - 1); q.trans = 10;
            if (!gi.nnc.empty()) q = gi.nnc[size_t(std::min(at, int(gi.nnc.size()) - 1))];
            gi.nnc.insert(gi.nnc.begin() + at, size_t(n), q);
        };
        ga.erase = [&](int r0, int r1) { r1 = std::min(r1, int(gi.nnc.size()) - 1); if (r0 <= r1) gi.nnc.erase(gi.nnc.begin() + r0, gi.nnc.begin() + r1 + 1); };
        if (!gi.nnc.empty()) {
            const float h = std::min(ImGui::GetFontSize() * 16, ImGui::GetFrameHeight() * (float(gi.nnc.size()) + 2.6f) + 4);
            if (ui::gridTable("nnc", ga, {"I1", "J1", "K1", "I2", "J2", "K2", "Trans"}, h, nullptr, "%.5g", true)) markDirty(true);
        }
        if (ImGui::Button(TRL("+ NNC"))) { NNCInput c; c.i2 = std::min(1, gi.nx - 1); c.trans = 10; gi.nnc.push_back(c); markDirty(true); }
        ImGui::SameLine(); ImGui::TextDisabled("%s", TR("transmissibility in simulator units (cDarcy included), 1-based indices"));
    }
}

// ============================================================================ AQUIFERS
void App::editAquifers() {
    auto& aq = case_.aquifers;
    const GridInput& g = case_.grid;
    bool ch = false;
    if (ImGui::Button(TRL("+ Fetkovich aquifer"))) {
        AquiferInput a; a.id = int(aq.size()) + 1; a.type = AquiferType::Fetkovich;
        a.depth = case_.equil.owc; a.volume = case_.units == UnitSystem::Field ? 2e9 : 3e8; a.pi = case_.units == UnitSystem::Field ? 50 : 5;
        a.ct = case_.units == UnitSystem::Field ? 1e-5 : 1.5e-4;
        AquiferConnection c; c.i1 = 0; c.i2 = g.nx - 1; c.j1 = c.j2 = 0; c.k1 = 0; c.k2 = g.nz - 1; c.face = 'Y'; c.minus = true;
        a.conns.push_back(c); aq.push_back(a); ch = true;
    }
    ImGui::SameLine();
    if (ImGui::Button(TRL("+ Carter-Tracy aquifer"))) {
        AquiferInput a; a.id = int(aq.size()) + 1; a.type = AquiferType::CarterTracy;
        a.depth = case_.equil.owc; a.perm = 200; a.poro = 0.22; a.thickness = case_.units == UnitSystem::Field ? 100 : 30;
        a.rInner = case_.units == UnitSystem::Field ? 3000 : 900; a.angle = 90;
        a.ct = case_.units == UnitSystem::Field ? 1e-5 : 1.5e-4;
        AquiferConnection c; c.i1 = 0; c.i2 = g.nx - 1; c.j1 = c.j2 = g.ny - 1; c.k1 = 0; c.k2 = g.nz - 1; c.face = 'Y'; c.minus = false;
        a.conns.push_back(c); aq.push_back(a); ch = true;
    }
    if (aq.empty()) {
        ui::hint("No aquifers. A Fetkovich aquifer is a finite tank with a productivity index (pseudo-steady state); Carter-Tracy "
                 "approximates the van Everdingen-Hurst infinite-acting radial aquifer. Both connect to cell faces (AQUANCON).");
    }
    // C16: Excel-style tables - all aquifers in one table (AQUFETP / AQUCT), connections of the selected aquifer below
    static int aqSel = 0;
    if (!aq.empty()) {
        ui::header("Aquifers (AQUFETP / AQUCT)");
        ui::RecordTable<AquiferInput> t;
        t.integer("ID", [](AquiferInput& a) -> int& { return a.id; })
         .choiceFn("Model", [](AquiferInput& a) { return int(a.type); }, [](AquiferInput& a, int k) { a.type = AquiferType(k); }, {"Fetkovich", "Carter-Tracy"})
         .num("Datum depth", [](AquiferInput& a) -> double& { return a.depth; })
         .num("Initial p (-1 = equil.)", [](AquiferInput& a) -> double& { return a.p0; })
         .num("Total compr.", [](AquiferInput& a) -> double& { return a.ct; })
         .num("Volume (Fetk.)", [](AquiferInput& a) -> double& { return a.volume; })
         .num("PI (Fetk.)", [](AquiferInput& a) -> double& { return a.pi; })
         .num("Perm mD (C-T)", [](AquiferInput& a) -> double& { return a.perm; })
         .num("Porosity (C-T)", [](AquiferInput& a) -> double& { return a.poro; })
         .num("Inner radius (C-T)", [](AquiferInput& a) -> double& { return a.rInner; })
         .num("Thickness (C-T)", [](AquiferInput& a) -> double& { return a.thickness; })
         .num("Angle deg (C-T)", [](AquiferInput& a) -> double& { return a.angle; })
         .num("Water visc (-1 = PVTW)", [](AquiferInput& a) -> double& { return a.visc; });
        t.canInsert = false;   // new aquifers through the buttons above (sensible defaults per model)
        t.onChange = [](AquiferInput& a) { a.angle = std::clamp(a.angle, 1.0, 360.0); a.id = std::clamp(a.id, 1, 99); if (a.p0 < 0) a.p0 = -1; if (a.visc < 0) a.visc = -1; };
        const float h = std::min(ImGui::GetFontSize() * 12, ImGui::GetFrameHeight() * (float(aq.size()) + 1.8f) + 4);
        ch |= t.draw("aquifers", aq, h, &aqSel, "%.5g", true);
        ui::hint("Fetkovich: finite tank with a productivity index (volume in reservoir volume units, PI per day per pressure). "
                 "Carter-Tracy: infinite-acting radial aquifer (van Everdingen-Hurst). Initial pressure -1: in equilibrium with the "
                 "reservoir at the datum; water viscosity -1: from PVTW. Select a row to edit its connections.");
    }
    aqSel = aq.empty() ? 0 : std::clamp(aqSel, 0, int(aq.size()) - 1);
    if (!aq.empty()) {
        AquiferInput& a = aq[size_t(aqSel)];
        char hdr[96]; std::snprintf(hdr, sizeof(hdr), "Connections of aquifer %d (AQUANCON, 1-based cell ranges)", a.id);
        ui::header(hdr);
        const int nx = g.nx, ny = g.ny, nz = g.nz;
        ui::RecordTable<AquiferConnection> ct;
        ct.integer("I1", [](AquiferConnection& c) -> int& { return c.i1; }, 1).integer("I2", [](AquiferConnection& c) -> int& { return c.i2; }, 1)
          .integer("J1", [](AquiferConnection& c) -> int& { return c.j1; }, 1).integer("J2", [](AquiferConnection& c) -> int& { return c.j2; }, 1)
          .integer("K1", [](AquiferConnection& c) -> int& { return c.k1; }, 1).integer("K2", [](AquiferConnection& c) -> int& { return c.k2; }, 1)
          .choiceFn("Face", [](AquiferConnection& c) { return (c.face == 'X' ? 0 : c.face == 'Y' ? 2 : 4) + (c.minus ? 0 : 1); },
                    [](AquiferConnection& c, int f) { c.face = "XYZ"[f / 2]; c.minus = (f % 2) == 0; },
                    {"I- (west)", "I+ (east)", "J- (south)", "J+ (north)", "K- (top)", "K+ (bottom)"})
          .num("Mult", [](AquiferConnection& c) -> double& { return c.mult; });
        ct.onChange = [nx, ny, nz](AquiferConnection& c) {
            c.i1 = std::clamp(c.i1, 0, nx - 1); c.i2 = std::clamp(c.i2, 0, nx - 1); c.j1 = std::clamp(c.j1, 0, ny - 1);
            c.j2 = std::clamp(c.j2, 0, ny - 1); c.k1 = std::clamp(c.k1, 0, nz - 1); c.k2 = std::clamp(c.k2, 0, nz - 1);
        };
        ct.newRow = [&](int) { AquiferConnection c; c.i2 = nx - 1; c.k2 = nz - 1; c.face = 'Y'; c.minus = true; return c; };
        const float h = std::min(ImGui::GetFontSize() * 14, ImGui::GetFrameHeight() * (float(a.conns.size()) + 1.8f) + 4);
        ch |= ct.draw("aqconn", a.conns, h, nullptr, "%.4g", true);
        if (ImGui::SmallButton(TRL("+ connection"))) { a.conns.push_back(ct.newRow(0)); ch = true; }
    }
    if (ch) markDirty();
    // influx results from the last run
    std::lock_guard<std::mutex> l(runner_.mutex());
    const RunResults& R = runner_.results();
    if (R.valid && R.summary.size() > 1 && !aq.empty()) {
        ui::header("Aquifer influx (last run)");
        std::vector<double> t, we, qa;
        const double xs = 1.0 / 365.25;
        for (auto& sp : R.summary) { t.push_back(sp.time * xs); we.push_back(sp.faqt); qa.push_back(sp.faqr); }
        if (ImPlot::BeginPlot("##aqinf", ImVec2(-1, ImGui::GetFontSize() * 13))) {
            ImPlot::SetupAxes("years", "cumulative influx", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "rate / day", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
            ImPlotSpec a; a.LineColor = col4(Pal::Water); a.LineWeight = 2;
            ImPlot::PlotLine("We (cum)", t.data(), we.data(), int(t.size()), a);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
            ImPlotSpec b; b.LineColor = col4(Pal::Accent); b.LineWeight = 1.5f;
            ImPlot::PlotLine("influx rate", t.data(), qa.data(), int(t.size()), b);
            ImPlot::EndPlot();
        }
    }
}

// ============================================================================ PVT
void App::editPVT() {
    if (!editFluidModel()) return;   // compositional / CO2 storage editors
    const Units U = Units::of(case_.units);
    PropsInput& p = case_.props;
    const bool live = case_.phases.oil && case_.phases.gas && case_.phases.disgas;
    if (live && !case_.schedule.empty()) {   // DRSDT from the start of the run (later changes: deck editor)
        ReportStep& s0 = case_.schedule[0];
        bool on = s0.drsdt >= 0;
        if (ImGui::Checkbox(TRL("Limit gas re-dissolution (DRSDT)"), &on)) { s0.drsdt = on ? 0.0 : -2; markDirty(); }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Maximum rate at which Rs may increase in a cell (%s per day).\n"
                              "DRSDT 0: free gas cannot re-dissolve, the GOR of undersaturated oil cannot rise.\n"
                              "Common in field models to suppress unphysical re-dissolution in coarse cells.", U.rs);
        if (on) {
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
            double v = s0.drsdt;
            if (ImGui::InputDouble("per day##drsdt", &v, 0, 0, "%.4g") && v >= 0) { s0.drsdt = v; markDirty(); }
            ImGui::SameLine();
            if (ImGui::Checkbox(TRL("only cells with free gas (FREE)"), &s0.drsdtFree)) markDirty();
        }
    }
    if (!ImGui::BeginTabBar("pvt")) return;
    if (case_.phases.oil && ui::tabItem(live ? "Oil (PVTO)" : "Oil (PVDO)")) {
        static int selRec = 0;
        if (live) {
            ui::hint("Saturated live-oil table: each row is a bubble point (Rs, Pb). Select a row to edit its undersaturated branch (p > Pb).");
            char c0[48], c1[48];
            std::snprintf(c0, sizeof(c0), "Rs (%s)", U.rs); std::snprintf(c1, sizeof(c1), "Pb (%s)", U.p);
            if (ui::editTable<PVTORecord>("pvto", p.pvto, {c0, c1, "Bo (rb/stb)", "\xCE\xBCo (cP)"},
                                          [](PVTORecord& r, int c) -> double* { return c == 0 ? &r.rs : c == 1 ? &r.pb : c == 2 ? &r.bo : &r.mu; },
                                          ImGui::GetFontSize() * 12, &selRec))
                markDirty();
            if (selRec >= 0 && selRec < int(p.pvto.size())) {
                ImGui::TextDisabled("Undersaturated branch of record %d (Rs = %.4g)", selRec + 1, p.pvto[selRec].rs);
                if (ui::editTable<PVDRow>("under", p.pvto[selRec].under, {"p", "Bo", "\xCE\xBCo"},
                                          [](PVDRow& r, int c) -> double* { return c == 0 ? &r.p : c == 1 ? &r.B : &r.mu; },
                                          ImGui::GetFontSize() * 6))
                    markDirty();
            }
        } else {
            if (ui::editTable<PVDRow>("pvdo", p.pvdo, {"p", "Bo", "\xCE\xBCo (cP)"},
                                      [](PVDRow& r, int c) -> double* { return c == 0 ? &r.p : c == 1 ? &r.B : &r.mu; }))
                markDirty();
        }
        // plots
        if (ImPlot::BeginSubplots("##oilplots", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 14))) {
            if (ImPlot::BeginPlot("Bo")) {
                ImPlot::SetupAxes(U.p, "Bo", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                if (live) {
                    std::vector<double> x, y, rs;
                    for (auto& r : p.pvto) { x.push_back(r.pb); y.push_back(r.bo); rs.push_back(r.rs); }
                    ImPlotSpec sp; sp.LineColor = col4(Pal::Oil); sp.LineWeight = 2; sp.Marker = ImPlotMarker_Circle;
                    ImPlot::PlotLine("Saturated", x.data(), y.data(), int(x.size()), sp);
                    for (auto& r : p.pvto) if (!r.under.empty()) {
                        std::vector<double> ux{r.pb}, uy{r.bo};
                        for (auto& u : r.under) { ux.push_back(u.p); uy.push_back(u.B); }
                        ImPlotSpec s2; s2.LineColor = ImVec4(0.6f, 0.9f, 0.6f, 0.7f); s2.Flags = ImPlotItemFlags_NoLegend;
                        ImPlot::PlotLine("##u", ux.data(), uy.data(), int(ux.size()), s2);
                    }
                } else {
                    std::vector<double> x, y; for (auto& r : p.pvdo) { x.push_back(r.p); y.push_back(r.B); }
                    ImPlotSpec sp; sp.LineColor = col4(Pal::Oil); sp.LineWeight = 2;
                    ImPlot::PlotLine("Bo", x.data(), y.data(), int(x.size()), sp);
                }
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot(live ? "Rs and viscosity" : "Viscosity")) {
                ImPlot::SetupAxes(U.p, live ? U.rs : "cP", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                if (live) ImPlot::SetupAxis(ImAxis_Y2, "cP", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
                std::vector<double> x, y, m;
                if (live) for (auto& r : p.pvto) { x.push_back(r.pb); y.push_back(r.rs); m.push_back(r.mu); }
                else for (auto& r : p.pvdo) { x.push_back(r.p); m.push_back(r.mu); }
                ImPlotSpec a; a.LineColor = col4(Pal::Gas); a.LineWeight = 2;
                ImPlotSpec b; b.LineColor = col4(Pal::Warn); b.LineWeight = 2;
                if (live) { ImPlot::PlotLine("Rs", x.data(), y.data(), int(x.size()), a); ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2); }
                ImPlot::PlotLine("\xCE\xBCo", x.data(), m.data(), int(x.size()), b);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        ImGui::EndTabItem();
    }
    const bool vapo = case_.phases.gas && case_.phases.oil && case_.phases.vapoil;
    if (vapo && ui::tabItem("Gas (PVTG)")) {
        static int selG = 0;
        char bc[48]; std::snprintf(bc, sizeof(bc), "Bg (%s)", case_.units == UnitSystem::Field ? "rb/Mscf" : "rm3/sm3");
        char rvc[48]; std::snprintf(rvc, sizeof(rvc), "Rv sat (%s)", U.rv);
        ui::hint("Wet-gas table: each row is a gas pressure with its saturated (dew-point) Rv, Bg and \xCE\xBCg. Select a row to edit its "
                 "undersaturated branch (Rv below saturation, down to 0). Records without a branch borrow the shape of the next "
                 "higher-pressure record that has one (as OPM Flow / ECLIPSE).");
        auto fix = [](PVTGRecord& r) { if (r.rows.empty()) r.rows.push_back({0.0, 1.0, 0.02}); };
        for (auto& r : p.pvtg) fix(r);
        bool chg = ui::editTable<PVTGRecord>("pvtg", p.pvtg, {U.p, rvc, bc, "\xCE\xBCg (cP)"},
            [&](PVTGRecord& r, int c) -> double* { fix(r); return c == 0 ? &r.pg : c == 1 ? &r.rows.front().rv : c == 2 ? &r.rows.front().B : &r.rows.front().mu; },
            ImGui::GetFontSize() * 12, &selG);
        if (selG >= 0 && selG < int(p.pvtg.size())) {
            PVTGRecord& R = p.pvtg[size_t(selG)];
            ImGui::TextDisabled("Undersaturated branch of record %d (p = %.5g, Rv sat = %.4g)", selG + 1, R.pg, R.rows.front().rv);
            std::vector<PVTGRow> und(R.rows.begin() + 1, R.rows.end());
            char rvu[48]; std::snprintf(rvu, sizeof(rvu), "Rv (%s)", U.rv);
            if (ui::editTable<PVTGRow>("pvtgu", und, {rvu, "Bg", "\xCE\xBCg"},
                                       [](PVTGRow& r, int c) -> double* { return c == 0 ? &r.rv : c == 1 ? &r.B : &r.mu; },
                                       ImGui::GetFontSize() * 6)) {
                R.rows.resize(1); R.rows.insert(R.rows.end(), und.begin(), und.end()); chg = true;
            }
        }
        if (chg) {
            p.pvdg.clear(); for (auto& g : p.pvtg) if (!g.rows.empty()) p.pvdg.push_back({g.pg, g.rows.front().B, g.rows.front().mu});
            markDirty();
        }
        if (ImPlot::BeginSubplots("##gasplots", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 14))) {
            std::vector<double> x, bg, rv, mu;
            for (auto& g : p.pvtg) if (!g.rows.empty()) { x.push_back(g.pg); bg.push_back(g.rows.front().B); rv.push_back(g.rows.front().rv); mu.push_back(g.rows.front().mu); }
            if (ImPlot::BeginPlot("Bg (saturated)")) {
                ImPlot::SetupAxes(U.p, bc, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                ImPlotSpec a; a.LineColor = col4(Pal::Gas); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle;
                ImPlot::PlotLine("Bg", x.data(), bg.data(), int(x.size()), a);
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Dew-point Rv and viscosity")) {
                ImPlot::SetupAxes(U.p, U.rv, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxis(ImAxis_Y2, "cP", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
                ImPlotSpec a; a.LineColor = col4(Pal::Oil); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle;
                ImPlotSpec b; b.LineColor = col4(Pal::Warn); b.LineWeight = 2;
                ImPlot::PlotLine("Rv sat", x.data(), rv.data(), int(x.size()), a);
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                ImPlot::PlotLine("\xCE\xBCg", x.data(), mu.data(), int(x.size()), b);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        ImGui::EndTabItem();
    }
    if (case_.phases.gas && !vapo && ui::tabItem("Gas (PVDG)")) {
        char bc[48]; std::snprintf(bc, sizeof(bc), "Bg (%s)", case_.units == UnitSystem::Field ? "rb/Mscf" : "rm3/sm3");
        if (ui::editTable<PVDRow>("pvdg", p.pvdg, {U.p, bc, "\xCE\xBCg (cP)"},
                                  [](PVDRow& r, int c) -> double* { return c == 0 ? &r.p : c == 1 ? &r.B : &r.mu; }))
            markDirty();
        if (ImPlot::BeginPlot("Dry gas", ImVec2(-1, ImGui::GetFontSize() * 14))) {
            ImPlot::SetupAxes(U.p, bc, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "cP", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            std::vector<double> x, y, m; for (auto& r : p.pvdg) { x.push_back(r.p); y.push_back(r.B); m.push_back(r.mu); }
            ImPlotSpec a; a.LineColor = col4(Pal::Gas); a.LineWeight = 2;
            ImPlotSpec b; b.LineColor = col4(Pal::Warn); b.LineWeight = 2;
            ImPlot::PlotLine("Bg", x.data(), y.data(), int(x.size()), a);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
            ImPlot::PlotLine("\xCE\xBCg", x.data(), m.data(), int(x.size()), b);
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }
    if (ui::tabItem("Water, rock, densities")) {
        bool ch = false;
        ui::header("Water (PVTW)");
        if (ui::beginProps("pvtw")) {
            ch |= ui::propD("Reference pressure", &p.pvtw.pref);
            ch |= ui::propD("Bw at reference", &p.pvtw.bw);
            ch |= ui::propD("Compressibility cw (1/p)", &p.pvtw.cw, "%.4g");
            ch |= ui::propD("Viscosity \xCE\xBCw (cP)", &p.pvtw.visc);
            ch |= ui::propD("Viscosibility (1/p)", &p.pvtw.cv, "%.4g");
            ui::endProps();
        }
        ui::header("Rock (ROCK)");
        if (ui::beginProps("rock")) {
            ch |= ui::propD("Reference pressure", &p.rock.pref);
            ch |= ui::propD("Pore compressibility cr (1/p)", &p.rock.cr, "%.4g");
            ui::endProps();
        }
        ui::header("Surface densities (DENSITY)");
        if (ui::beginProps("dens")) {
            ch |= ui::propD("Oil", &p.density.oil); ch |= ui::propD("Water", &p.density.water); ch |= ui::propD("Gas", &p.density.gas, "%.5g");
            char b[64];
            const double rhoW = case_.units == UnitSystem::Field ? 62.37 : 999.0;
            const double sg = p.density.oil / rhoW;
            std::snprintf(b, sizeof(b), "%.1f \xC2\xB0" "API", 141.5 / sg - 131.5); ui::propText("Oil gravity", b);
            std::snprintf(b, sizeof(b), "%.3f (air = 1)", p.density.gas / (case_.units == UnitSystem::Field ? 0.0764 : 1.225)); ui::propText("Gas gravity", b);
            ui::endProps();
        }
        if (ch) markDirty();
        ImGui::EndTabItem();
    }
    if (case_.phases.oil && ui::tabItem("Correlations")) {
        ui::hint("Generate a complete black-oil PVT description from a few field measurements. "
                 "Standing (Rs, Bo), Beggs-Robinson + Vasquez-Beggs (oil viscosity), Vasquez-Beggs (undersaturated co), "
                 "Sutton pseudo-criticals + Hall-Yarborough (z-factor, Bg), Lee-Gonzalez-Eakin (gas viscosity).");
        if (ui::beginProps("corr")) {
            ui::propD("Oil gravity (\xC2\xB0" "API)", &corr_.api, "%.1f");
            ui::propD("Gas specific gravity (air=1)", &corr_.gasGravity, "%.3f");
            ui::propD("Reservoir temperature (\xC2\xB0" "F)", &corr_.temperatureF, "%.1f");
            ui::propD("Max bubble point (psia)", &corr_.pbMax, "%.0f");
            ui::propD("Max table pressure (psia)", &corr_.pMax, "%.0f");
            ui::propI("Table rows", &corr_.nRows, 4, 40);
            ui::endProps();
        }
        if (ui::accentButton("Generate PVT tables", ImVec2(-FLT_MIN, 0))) {
            try {
                makeCorrelationPVT(corr_, case_.units, case_.props, !live);
                markDirty();
                notify("PVT tables generated from correlations");
            } catch (const std::exception& e) { notify(std::string("Correlation error: ") + e.what()); }
        }
        ImGui::TextDisabled("%s", TR("Pressures are entered in psia; tables are converted to the case units."));
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

// ============================================================================ SCAL
void App::editSCAL() {
    PropsInput& p = case_.props;
    const bool wo = case_.phases.water && case_.phases.oil, go = case_.phases.gas && case_.phases.oil;
    if (!wo && !go) { ui::hint("Relative permeability is only needed with two or more phases."); return; }
    if (ImGui::CollapsingHeader(TRL("Corey generator"), ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ui::beginProps("corey", 0.56f)) {
            ui::propSliderD("Swc (connate water)", &corey_.swc, 0.0, 0.5, "%.3f");
            if (wo) {
                ui::propSliderD("Sorw (residual oil to water)", &corey_.sorw, 0.0, 0.5, "%.3f");
                ui::propSliderD("nw (water exponent)", &corey_.nw, 1.0, 6.0, "%.2f");
                ui::propSliderD("now (oil exponent, o/w)", &corey_.now, 1.0, 6.0, "%.2f");
                ui::propSliderD("krw at Sorw", &corey_.krwMax, 0.05, 1.0, "%.2f");
            }
            ui::propSliderD("kro at Swc", &corey_.kroMax, 0.05, 1.0, "%.2f");
            if (go) {
                ui::propSliderD("Sgc (critical gas)", &corey_.sgc, 0.0, 0.3, "%.3f");
                ui::propSliderD("Sorg (residual oil to gas)", &corey_.sorg, 0.0, 0.5, "%.3f");
                ui::propSliderD("ng (gas exponent)", &corey_.ng, 1.0, 6.0, "%.2f");
                ui::propSliderD("nog (oil exponent, g/o)", &corey_.nog, 1.0, 6.0, "%.2f");
                ui::propSliderD("krg at max Sg", &corey_.krgMax, 0.05, 1.0, "%.2f");
            }
            if (wo) {
                ui::propD("Pcow entry pressure", &corey_.pcEntry, "%.4g", "Capillary pressure at Swc (0 = no capillary pressure).");
                ui::propSliderD("Pcow exponent", &corey_.pcExponent, 0.5, 6.0, "%.2f");
            }
            ui::endProps();
        }
        if (ui::accentButton("Generate SWOF / SGOF", ImVec2(-FLT_MIN, 0))) {
            try {
                if (wo) makeCoreySWOF(corey_, p.swof);
                if (go) makeCoreySGOF(corey_, p.sgof);
                markDirty();
            } catch (const std::exception& e) { notify(e.what()); }
        }
    }
    if (ImGui::CollapsingHeader(TRL("van Genuchten - Mualem generator (TOUGH2 / TOUGH3 style)"))) {
        // C27: the standard curves of the TOUGH family (RPCAP / ICP type 7), e.g. for CO2-brine, vadose-zone or geothermal work
        static const char* kTargets[] = {"SWOF (water wetting, oil non-wetting)", "SGOF (liquid wetting, gas non-wetting)"};
        if (!wo) vgTarget_ = 1;
        if (!go) vgTarget_ = 0;
        if (ui::beginProps("vgm", 0.56f)) {
            if (wo && go) ui::propCombo("Table", &vgTarget_, kTargets, 2);
            ui::propSliderD(vgTarget_ == 0 ? "Swr (residual water)" : "Slr (residual liquid, >= Swc)", &vg_.swr, 0.0, 0.6, "%.3f");
            ui::propSliderD(vgTarget_ == 0 ? "Snr (residual oil)" : "Sgr (residual gas)", &vg_.snr, 0.0, 0.4, "%.3f");
            ui::propSliderD("m = 1 - 1/n (TOUGH lambda)", &vg_.m, 0.1, 0.9, "%.3f");
            ui::propD("P0 = 1/alpha (pressure units, 0 = no Pc)", &vg_.p0, "%.4g", "Capillary strength: Pc = P0 (Sc^(-1/m) - 1)^(1-m), Sc = (S - Sr)/(1 - Sr).");
            ui::propD("Pc cap (0 = 100 P0)", &vg_.pcMax, "%.4g");
            ui::propSliderD("kr wetting end point", &vg_.krwMax, 0.05, 1.0, "%.2f");
            ui::propSliderD("kr non-wetting end point", &vg_.krnMax, 0.05, 1.0, "%.2f");
            static const char* kNw[] = {"Parker / VGM (1 - Se)^1/2 (1 - Se^(1/m))^2m", "Corey (TOUGH RPCAP 7): (1 - S)^2 (1 - S^2)"};
            ui::propCombo("Non-wetting kr", &vg_.nonWetting, kNw, 2);
            ui::propI("Rows", &vg_.rows, 6, 60);
            ui::endProps();
        }
        if (ui::accentButton(vgTarget_ == 0 ? "Generate SWOF (van Genuchten - Mualem)" : "Generate SGOF (van Genuchten - Mualem)", ImVec2(-FLT_MIN, 0))) {
            try {
                if (vgTarget_ == 0) makeVanGenuchtenSWOF(vg_, p.swof);
                else makeVanGenuchtenSGOF(vg_, p.swof.empty() ? 0.0 : p.swof.front().s, p.sgof);
                markDirty();
                notify("van Genuchten - Mualem table generated");
            } catch (const std::exception& e) { notify(e.what()); }
        }
        ImGui::TextDisabled("%s", TR("Typical sandstone: m 0.4 - 0.6, P0 0.1 - 0.5 bar; CO2-brine (Pruess benchmark): m 0.457, Swr 0.3, P0 0.196 bar."));
    }
    if (wo && go) {
        static const char* kModels[] = {"ECLIPSE default (saturation-weighted)", "Stone I (STONE1)", "Stone II (STONE2)"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 16);
        if (ImGui::Combo("Three-phase oil kr", &p.kroModel, kModels, 3)) markDirty();
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", TR("How kro is built from krow (SWOF) and krog (SGOF) when water, oil and gas flow together.\n"
                              "Default: kro = (Sg krog + (Sw - Swc) krow) / (Sg + Sw - Swc).\n"
                              "Stone II: kro = krocw ((krow/krocw + krw)(krog/krocw + krg) - krw - krg).\n"
                              "Stone I: kro = beta krow krog / krocw (normalised saturations, exponent STONE1EX)."));
        if (p.kroModel == 1) {
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
            double e = p.stone1Exp;
            if (ImGui::InputDouble("STONE1EX", &e, 0, 0, "%.3g") && e > 0) { p.stone1Exp = e; markDirty(); }
        }
    }
    if (ImGui::BeginTabBar("scal")) {
        if (wo && ui::tabItem("Water-oil (SWOF)")) {
            if (ui::editTable<SatRow>("swof", p.swof, {"Sw", "krw", "krow", "Pcow"},
                                      [](SatRow& r, int c) -> double* { return c == 0 ? &r.s : c == 1 ? &r.krw : c == 2 ? &r.kro : &r.pc; },
                                      ImGui::GetFontSize() * 10))
                markDirty();
            if (ImPlot::BeginSubplots("##wo", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 15))) {
                std::vector<double> s, a, b, pc;
                for (auto& r : p.swof) { s.push_back(r.s); a.push_back(r.krw); b.push_back(r.kro); pc.push_back(r.pc); }
                if (ImPlot::BeginPlot("Relative permeability")) {
                    ImPlot::SetupAxes("Sw", "kr", 0, 0);
                    ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
                    ImPlotSpec sw; sw.LineColor = col4(Pal::Water); sw.LineWeight = 2.5f;
                    ImPlotSpec so; so.LineColor = col4(Pal::Oil); so.LineWeight = 2.5f;
                    ImPlot::PlotLine("krw", s.data(), a.data(), int(s.size()), sw);
                    ImPlot::PlotLine("krow", s.data(), b.data(), int(s.size()), so);
                    ImPlot::EndPlot();
                }
                // fractional flow + Welge tangent (Buckley-Leverett)
                if (ImPlot::BeginPlot("Fractional flow (Buckley-Leverett)")) {
                    ImPlot::SetupAxes("Sw", "fw", 0, 0);
                    ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
                    double muw = p.pvtw.visc, muo = 1.0;
                    if (case_.phases.gas && case_.phases.disgas && !p.pvto.empty()) muo = p.pvto.back().mu;
                    else if (!p.pvdo.empty()) muo = p.pvdo[p.pvdo.size() / 2].mu;
                    std::vector<double> fs, fw;
                    double swc = s.empty() ? 0 : s.front();
                    double bestSlope = 0, sf = swc, fwf = 0;
                    for (int i = 0; i <= 200 && !s.empty(); ++i) {
                        const double sv = swc + (s.back() - swc) * i / 200.0;
                        auto interp = [&](const std::vector<double>& y) {
                            size_t k = std::upper_bound(s.begin(), s.end(), sv) - s.begin();
                            if (k == 0) return y.front();
                            if (k >= s.size()) return y.back();
                            const double t = (sv - s[k - 1]) / (s[k] - s[k - 1]); return y[k - 1] + t * (y[k] - y[k - 1]);
                        };
                        const double lw = interp(a) / muw, lo = interp(b) / muo;
                        const double f = (lw + lo) > 0 ? lw / (lw + lo) : 1.0;
                        fs.push_back(sv); fw.push_back(f);
                        if (sv > swc + 1e-6) { const double sl = f / (sv - swc); if (sl > bestSlope) { bestSlope = sl; sf = sv; fwf = f; } }
                    }
                    ImPlotSpec sp; sp.LineColor = col4(Pal::Water); sp.LineWeight = 2.5f;
                    ImPlot::PlotLine("fw", fs.data(), fw.data(), int(fs.size()), sp);
                    if (bestSlope > 0) {
                        const double tx[2] = {swc, std::min(1.0, swc + 1.0 / bestSlope)}, ty[2] = {0, std::min(1.0, bestSlope * (tx[1] - swc))};
                        ImPlotSpec ts; ts.LineColor = col4(Pal::Warn); ts.LineWeight = 1.5f;
                        ImPlot::PlotLine("Welge tangent", tx, ty, 2, ts);
                        ImPlotSpec ms; ms.Marker = ImPlotMarker_Circle; ms.MarkerSize = 5; ms.MarkerFillColor = col4(Pal::Warn);
                        ImPlot::PlotScatter("Shock front", &sf, &fwf, 1, ms);
                        char b2[96]; std::snprintf(b2, sizeof(b2), "Swf = %.3f, fwf = %.3f", sf, fwf);
                        ImPlot::PlotText(b2, 0.3, 0.08);
                    }
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            ImGui::EndTabItem();
        }
        if (go && ui::tabItem("Gas-oil (SGOF)")) {
            if (ui::editTable<SatRow>("sgof", p.sgof, {"Sg", "krg", "krog", "Pcog"},
                                      [](SatRow& r, int c) -> double* { return c == 0 ? &r.s : c == 1 ? &r.krw : c == 2 ? &r.kro : &r.pc; },
                                      ImGui::GetFontSize() * 10))
                markDirty();
            if (ImPlot::BeginPlot("Gas-oil relative permeability", ImVec2(-1, ImGui::GetFontSize() * 15))) {
                ImPlot::SetupAxes("Sg", "kr", 0, 0);
                ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
                std::vector<double> s, a, b;
                for (auto& r : p.sgof) { s.push_back(r.s); a.push_back(r.krw); b.push_back(r.kro); }
                ImPlotSpec sg; sg.LineColor = col4(Pal::Gas); sg.LineWeight = 2.5f;
                ImPlotSpec so; so.LineColor = col4(Pal::Oil); so.LineWeight = 2.5f;
                ImPlot::PlotLine("krg", s.data(), a.data(), int(s.size()), sg);
                ImPlot::PlotLine("krog", s.data(), b.data(), int(s.size()), so);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if (go && ui::tabItem("Hysteresis")) {
            editHysteresis(case_.comp.enabled);
            ImGui::EndTabItem();
        }
        if (ui::tabItem("End-point scaling")) {
            editEndPointScaling();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
}

// End-point scaling: ENDSCALE + per-cell saturation end points (constant or per layer), live preview of the scaled
// curves against the tables of SATNUM region 1.
void App::editEndPointScaling() {
    GridInput& g = case_.grid;
    const bool wo = case_.phases.water && case_.phases.oil, go = case_.phases.gas && case_.phases.oil;
    const float fs = ImGui::GetFontSize();
    ImGui::TextWrapped("%s", TR("Scales the relative permeability and capillary pressure curves per cell to its own saturation end points "
                       "(ECLIPSE two-point scaling: kr between the critical and maximum saturation, Pc between SWL and SWU). "
                       "Cells without a value use the table end points of their SATNUM region. Written as ENDSCALE + SWL ... in the deck."));
    if (ImGui::Checkbox(TRL("Enable end-point scaling (ENDSCALE)"), &g.endscale)) markDirty();
    SatFunctions sf;
    bool sfOk = true;
    try { sf.setup(case_, 0); } catch (...) { sfOk = false; }
    const SatFunctions::EndPoints T = sfOk ? sf.tableEndPoints() : SatFunctions::EndPoints{};
    struct Arr { const char* name; const char* meaning; std::vector<double>* v; double T; bool show; };
    Arr arrs[] = {
        {"SWL", "connate water (Pc and kro end)", &g.swl, T.swl, wo}, {"SWCR", "critical water (krw = 0 below)", &g.swcr, T.swcr, wo},
        {"SWU", "maximum water", &g.swu, T.swu, wo}, {"SOWCR", "residual oil to water", &g.sowcr, T.sowcr, wo},
        {"SGL", "connate gas", &g.sgl, T.sgl, go}, {"SGCR", "critical gas (krg = 0 below)", &g.sgcr, T.sgcr, go},
        {"SGU", "maximum gas", &g.sgu, T.sgu, go}, {"SOGCR", "residual oil to gas", &g.sogcr, T.sogcr, go},
    };
    const int n = g.size(), nz = std::max(1, g.nz), perLayer = std::max(1, n / nz);
    static int sel = 0;
    static double constVal = 0.2;
    if (ui::beginTable("##eps", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
        ui::tableColumn("Array", ImGuiTableColumnFlags_WidthFixed, fs * 4);
        ui::tableColumn("Meaning", ImGuiTableColumnFlags_WidthFixed, fs * 12);
        ui::tableColumn("Table", ImGuiTableColumnFlags_WidthFixed, fs * 3.5f);
        ui::tableColumn("Cells", ImGuiTableColumnFlags_WidthFixed, fs * 9);
        ui::tableColumn("", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();
        for (int i = 0; i < 8; ++i) {
            Arr& a = arrs[i];
            if (!a.show) continue;
            ImGui::PushID(i);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            if (ImGui::Selectable(a.name, sel == i, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) sel = i;
            ImGui::TableSetColumnIndex(1); ImGui::TextDisabled("%s", a.meaning);
            ImGui::TableSetColumnIndex(2); ImGui::Text("%.3f", a.T);
            ImGui::TableSetColumnIndex(3);
            if (a.v->empty()) ImGui::TextDisabled("%s", TR("table"));
            else {
                double lo = 1e9, hi = -1e9; int set = 0;
                for (double x : *a.v) if (x >= 0) { lo = std::min(lo, x); hi = std::max(hi, x); ++set; }
                if (set == 0) ImGui::TextDisabled("%s", TR("table"));
                else if (hi - lo < 1e-12) ImGui::Text("%.3f%s", lo, set < n ? " (part)" : "");
                else ImGui::Text("%.3f - %.3f", lo, hi);
            }
            ImGui::TableSetColumnIndex(4);
            if (ImGui::SmallButton(TRL("Constant"))) { a.v->assign(size_t(n), std::clamp(constVal, 0.0, 1.0)); g.endscale = true; markDirty(); }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Set every cell to the value below"));
            ImGui::SameLine();
            if (ImGui::SmallButton(TRL("Per layer"))) { sel = i; if (a.v->empty()) a.v->assign(size_t(n), a.T); g.endscale = true; markDirty(); }
            ImGui::SameLine();
            ImGui::BeginDisabled(a.v->empty());
            if (ImGui::SmallButton(TRL("Clear"))) { a.v->clear(); markDirty(); }
            ImGui::EndDisabled();
            ImGui::PopID();
        }
        ui::endTable();
    }
    ImGui::SetNextItemWidth(fs * 8);
    ImGui::InputDouble("Value for 'Constant'", &constVal, 0.01, 0.05, "%.3f");
    Arr* Ap = &arrs[std::clamp(sel, 0, 7)];
    if (!Ap->show) for (auto& a : arrs) if (a.show) { Ap = &a; break; }
    Arr& A = *Ap;
    // per-layer editor of the selected array
    if (!A.v->empty() && int(A.v->size()) == n && ImGui::TreeNodeEx(TRL("##layers"), ImGuiTreeNodeFlags_DefaultOpen, "%s per layer", A.name)) {
        if (ui::beginTable("##epsl", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit,
                              ImVec2(fs * 14, std::min(fs * 12, fs * 1.6f * (nz + 1))))) {
            ui::tableColumn("Layer", ImGuiTableColumnFlags_WidthFixed, fs * 3.5f);
            ui::tableColumn(A.name, ImGuiTableColumnFlags_WidthFixed, fs * 9);
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableHeadersRow();
            for (int k = 0; k < nz; ++k) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("%d", k + 1);
                ImGui::TableSetColumnIndex(1);
                double v = (*A.v)[size_t(k) * size_t(perLayer)];
                ImGui::PushID(k);
                ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::InputDouble("##v", &v, 0, 0, v < 0 ? "table" : "%.4f")) {
                    for (int q = 0; q < perLayer; ++q) (*A.v)[size_t(k) * size_t(perLayer) + size_t(q)] = std::min(1.0, v);
                    markDirty();
                }
                ImGui::PopID();
            }
            ui::endTable();
        }
        ImGui::TextDisabled("%s", TR("Values < 0 = table end point. Cell-by-cell values come from imported decks (SWL ... / EQUALS / BOX)."));
        ImGui::TreePop();
    }
    // preview: tables (dashed) against the scaled curves of the average cell end points
    if (!sfOk) { ui::hint("Saturation tables are incomplete - no preview."); return; }
    auto mean = [&](const std::vector<double>& v, double def) {
        double s = 0; int c = 0;
        for (double x : v) if (x >= 0) { s += x; ++c; }
        return c ? s / c : def;
    };
    SatFunctions::EndPoints E;
    E.swl = mean(g.swl, T.swl); E.swcr = std::max(E.swl, mean(g.swcr, T.swcr)); E.swu = std::max(E.swcr + 1e-6, mean(g.swu, T.swu));
    E.sgl = mean(g.sgl, T.sgl); E.sgcr = std::max(E.sgl, mean(g.sgcr, T.sgcr)); E.sgu = std::max(E.sgcr + 1e-6, mean(g.sgu, T.sgu));
    E.sowcr = mean(g.sowcr, T.sowcr); E.sogcr = mean(g.sogcr, T.sogcr);
    const int N = 120;
    if (ImPlot::BeginPlot("Table (dashed) vs scaled with the average cell end points", ImVec2(-1, std::max(fs * 14, ImGui::GetContentRegionAvail().y)))) {
        ImPlot::SetupAxes("Saturation", "kr", 0, 0);
        ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
        auto line = [&](const char* lab, ImU32 col, bool dashed, auto f) {
            std::vector<double> x(N + 1), y(N + 1);
            for (int i = 0; i <= N; ++i) { x[size_t(i)] = double(i) / N; y[size_t(i)] = f(x[size_t(i)]); }
            ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = dashed ? 1.2f : 2.5f;
            if (dashed) sp.LineColor.w = 0.55f;
            ImPlot::PlotLine(lab, x.data(), y.data(), N + 1, sp);
        };
        if (wo) {
            line("krw table", Pal::Water, true, [&](double s) { return sf.krw(AD3(s), nullptr).v; });
            line("krow table", Pal::Oil, true, [&](double s) { return sf.kro(AD3(s), AD3(0.0), nullptr).v; });
            line("krw scaled", Pal::Water, false, [&](double s) { return sf.krw(AD3(s), &E).v; });
            line("krow scaled", Pal::Oil, false, [&](double s) { return sf.kro(AD3(s), AD3(0.0), &E).v; });
        }
        if (go) {
            line("krg table (vs Sg)", Pal::Gas, true, [&](double s) { return sf.krg(AD3(s), nullptr).v; });
            line("krg scaled (vs Sg)", Pal::Gas, false, [&](double s) { return sf.krg(AD3(s), &E).v; });
        }
        ImPlot::EndPlot();
    }
}

// ============================================================================ INITIALISATION
void App::editInit() {
    EquilInput& e = case_.equil;
    const Units U = Units::of(case_.units);
    bool ch = false;
    ui::header("Hydrostatic equilibrium (EQUIL)");
    if (ui::beginProps("eq")) {
        ch |= ui::propD("Datum depth", &e.datumDepth, "%.6g");
        ch |= ui::propD("Pressure at datum", &e.datumPressure, "%.6g");
        if (case_.phases.water && case_.phases.oil) {
            ch |= ui::propD("Oil-water contact (OWC)", &e.owc, "%.6g");
            ch |= ui::propD("Pcow at OWC", &e.pcowOWC, "%.4g");
        }
        if (case_.phases.gas && case_.phases.oil) {
            ch |= ui::propD("Gas-oil contact (GOC)", &e.goc, "%.6g");
            ch |= ui::propD("Pcog at GOC", &e.pcogGOC, "%.4g");
        }
        ch |= ui::propI("Accuracy (NACC)", &e.nacc, 0, 20, "0: saturations at cell centre. N > 0: averaged over 2N horizontal slices per cell (smoother transition zones).");
        if (case_.phases.gas && case_.phases.oil && case_.phases.disgas) {
            bool sat = e.initialRs < 0;
            ui::label("Initial Rs", "Saturated at the GOC pressure is the ECLIPSE default; otherwise use a constant Rs (RSVD).");
            if (ImGui::Checkbox(TRL("Saturated at GOC"), &sat)) { e.initialRs = sat ? -1 : 1.0; ch = true; }
            if (!sat) { ui::label("Constant Rs"); ch |= ImGui::InputDouble("##rs", &e.initialRs, 0, 0, "%.5g"); }
        }
        if (case_.phases.gas && case_.phases.oil && case_.phases.vapoil) {
            bool sat = e.initialRv < 0;
            ui::label("Initial Rv", "Saturated at the GOC pressure (capped by the local dew point) is the ECLIPSE default; "
                                    "otherwise a constant Rv. An RVVD table below overrides both.");
            if (ImGui::Checkbox(TRL("Saturated at GOC##rv"), &sat)) { e.initialRv = sat ? -1 : 0.0; ch = true; }
            if (!sat) { ui::label("Constant Rv"); ch |= ImGui::InputDouble("##rvc", &e.initialRv, 0, 0, "%.5g"); }
        }
        ui::endProps();
    }
    // depth tables (RSVD / RVVD) - optional, override the constant / saturated choice
    struct DV { double d, v; };
    auto depthTable = [&](const char* id, const char* title, const char* vcol, std::vector<double>& dep, std::vector<double>& val) {
        if (!ImGui::TreeNode(id, "%s (%d rows)", title, int(dep.size()))) return;
        std::vector<DV> t; for (size_t i = 0; i < dep.size() && i < val.size(); ++i) t.push_back({dep[i], val[i]});
        if (ui::editTable<DV>(id, t, {U.len, vcol}, [](DV& r, int c) -> double* { return c == 0 ? &r.d : &r.v; })) {
            dep.clear(); val.clear(); for (auto& r : t) { dep.push_back(r.d); val.push_back(r.v); } ch = true;
        }
        if (!dep.empty() && ImGui::SmallButton(TRL("Clear table"))) { dep.clear(); val.clear(); ch = true; }
        ImGui::TreePop();
    };
    if (case_.phases.gas && case_.phases.oil && case_.phases.disgas) {
        char c[48]; std::snprintf(c, sizeof(c), "Rs (%s)", U.rs);
        depthTable("rsvd", "Rs vs depth (RSVD)", c, e.rsvdDepth, e.rsvdRs);
    }
    if (case_.phases.gas && case_.phases.oil && case_.phases.vapoil) {
        char c[48]; std::snprintf(c, sizeof(c), "Rv (%s)", U.rv);
        depthTable("rvvd", "Rv vs depth (RVVD)", c, e.rvvdDepth, e.rvvdRv);
    }
    if (ch) markDirty();
    ui::hint("Phase pressures are integrated from the datum with dP/dz = \xCF\x81(P) g; saturations follow from the capillary pressure "
             "curves. Contacts above the reservoir top give undersaturated oil (no free gas).");
    if (ui::accentButton("Compute initial state", ImVec2(-FLT_MIN, 0))) computeInit();
    if (!init_.error.empty()) { ImGui::PushTextWrapPos(0); ImGui::TextColored(col4(Pal::Err), "%s", init_.error.c_str()); ImGui::PopTextWrapPos(); }
    if (init_.valid) {
        if (init_.caseVersion != caseVersion_) ImGui::TextColored(col4(Pal::Warn), "Inputs changed since this initialisation.");
        ui::header("Fluids in place");
        if (ui::beginProps("fip")) {
            char b[96];
            auto big = [&](double v, const char* u) { if (std::fabs(v) >= 1e6) std::snprintf(b, sizeof(b), "%.4f MM %s", v / 1e6, u); else std::snprintf(b, sizeof(b), "%.6g %s", v, u); return b; };
            if (case_.phases.oil) ui::propText("Oil (STOIIP)", big(init_.ip.oil, U.liq));
            if (case_.phases.water) ui::propText("Water", big(init_.ip.water, U.liq));
            if (case_.phases.gas) ui::propText("Free gas", big(init_.ip.gasFree, U.gas));
            if (case_.phases.disgas && case_.phases.gas) ui::propText("Dissolved gas", big(init_.ip.gasDissolved, U.gas));
            if (case_.phases.vapoil && case_.phases.gas) {
                ui::propText("Vaporised oil (in gas)", big(init_.ip.oilVaporized, U.liq));
                ui::propText("Total oil (liquid + vapour)", big(init_.ip.oil + init_.ip.oilVaporized, U.liq));
            }
            ui::endProps();
        }
        if (ImPlot::BeginSubplots("##initplots", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 18), ImPlotSubplotFlags_LinkAllY)) {
            const int n = int(init_.depth.size());
            if (ImPlot::BeginPlot("Pressure vs depth")) {
                ImPlot::SetupAxes(U.p, U.len, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
                std::vector<double> pp(init_.p.begin(), init_.p.end());
                ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 2.5f; s.MarkerFillColor = col4(Pal::Pressure); s.MarkerLineColor = col4(Pal::Pressure);
                ImPlot::PlotScatter("cells", pp.data(), init_.depth.data(), n, s);
                ImPlotSpec h; h.Flags = ImPlotInfLinesFlags_Horizontal; h.LineColor = col4(Pal::Water);
                if (case_.phases.water) ImPlot::PlotInfLines("OWC", &e.owc, 1, h);
                h.LineColor = col4(Pal::Gas);
                if (case_.phases.gas) ImPlot::PlotInfLines("GOC", &e.goc, 1, h);
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Saturation vs depth")) {
                ImPlot::SetupAxes("S", U.len, 0, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
                ImPlot::SetupAxisLimits(ImAxis_X1, 0, 1, ImPlotCond_Always);
                std::vector<double> sw(init_.sw.begin(), init_.sw.end()), sg(init_.sg.begin(), init_.sg.end());
                ImPlotSpec a; a.Marker = ImPlotMarker_Circle; a.MarkerSize = 2.5f; a.MarkerFillColor = col4(Pal::Water); a.MarkerLineColor = col4(Pal::Water);
                ImPlot::PlotScatter("Sw", sw.data(), init_.depth.data(), n, a);
                if (case_.phases.gas) {
                    ImPlotSpec b2 = a; b2.MarkerFillColor = col4(Pal::Gas); b2.MarkerLineColor = col4(Pal::Gas);
                    ImPlot::PlotScatter("Sg", sg.data(), init_.depth.data(), n, b2);
                }
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
    }
}

// ============================================================================ WELLS
// ============================================================================ GROUPS
// Group tree (WELSPECS groups + GRUPTREE) and group controls (GCONPROD / GCONINJE / ... kept as deck keywords for
// OPM Flow; edited with schema forms on the Advanced page).
void App::editGroups() {
    namespace di = ressim::deckintel;
    auto& X = case_.extraKeywords;
    // GRUPTREE (first one in the schedule is the editable tree)
    int treeIdx = -1;
    for (int i = 0; i < int(X.size()); ++i) if (X[size_t(i)].name == "GRUPTREE") { treeIdx = i; break; }
    std::map<std::string, std::string> parent;
    di::Form tree;
    if (treeIdx >= 0) {
        tree = di::parseForm("GRUPTREE", di::splitLines(X[size_t(treeIdx)].body));
        for (auto& r : tree.records) if (r.size() >= 2) parent[di::unquote(r[0])] = di::unquote(r[1]) == "1*" ? "FIELD" : di::unquote(r[1]);
    }
    std::set<std::string> groups;
    for (auto& w : case_.wells) if (!w.group.empty()) groups.insert(w.group);
    for (auto& [c, p] : parent) { groups.insert(c); if (p != "FIELD") groups.insert(p); }
    groups.erase("FIELD");
    auto parentOf = [&](const std::string& g) { auto it = parent.find(g); return it == parent.end() ? std::string("FIELD") : it->second; };
    auto isDescendant = [&](std::string a, const std::string& of) {   // a below 'of'?
        for (int d = 0; d < 64 && a != "FIELD"; ++d) { a = parentOf(a); if (a == of) return true; }
        return false;
    };
    const ImGuiTreeNodeFlags tf = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_SpanAvailWidth;
    std::function<void(const std::string&, int)> node = [&](const std::string& g, int depth) {
        if (depth > 32) return;
        ImGuiTreeNodeFlags f = tf | (selGroup_ == g ? ImGuiTreeNodeFlags_Selected : 0);
        const bool open = ImGui::TreeNodeEx(g.c_str(), f, "%s", g.c_str());
        if (ImGui::IsItemClicked() && !ImGui::IsItemToggledOpen()) selGroup_ = g;
        if (!open) return;
        for (auto& c : groups) if (parentOf(c) == g) node(c, depth + 1);
        for (int i = 0; i < int(case_.wells.size()); ++i) {
            const auto& w = case_.wells[size_t(i)];
            if (w.group != g) continue;
            ImGui::PushID(i);
            const ImU32 col = w.type == WellType::Producer ? Pal::Oil : w.type == WellType::WaterInjector ? Pal::Water : Pal::Gas;
            ImGui::TextColored(col4(col), "%s", w.type == WellType::Producer ? "P " : "I ");
            ImGui::SameLine(0, 0);
            if (ImGui::Selectable(w.name.c_str(), selWell_ == i)) selWell_ = i;
            ImGui::PopID();
        }
        ImGui::TreePop();
    };
    const float fs = ImGui::GetFontSize();
    ImGui::BeginChild("##gtree", ImVec2(fs * 13, fs * 14), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
    node("FIELD", 0);
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("##gedit", ImVec2(0, fs * 14), ImGuiChildFlags_Borders);
    if (selGroup_.empty() || (selGroup_ != "FIELD" && !groups.count(selGroup_))) selGroup_ = "FIELD";
    ImGui::Text("Group %s", selGroup_.c_str());
    if (selGroup_ != "FIELD") {
        ImGui::SameLine(0, fs * 2);
        ImGui::SetNextItemWidth(fs * 9);
        const std::string cur = parentOf(selGroup_);
        if (ImGui::BeginCombo("Parent", cur.c_str())) {
            std::vector<std::string> cands = {"FIELD"};
            for (auto& g : groups) if (g != selGroup_ && !isDescendant(g, selGroup_)) cands.push_back(g);
            for (auto& c : cands) {
                if (!ImGui::Selectable(c.c_str(), c == cur) || c == cur) continue;
                // update / create GRUPTREE
                if (treeIdx < 0) {
                    ExtraKeyword e; e.name = "GRUPTREE"; e.section = DeckSection::Schedule; e.step = 0; e.beforeControls = true; e.note = "group tree (ResSim Studio)";
                    X.insert(X.begin(), e); treeIdx = 0;
                    tree = di::parseForm("GRUPTREE", {});
                }
                bool found = false;
                for (auto& r : tree.records) if (r.size() >= 2 && di::unquote(r[0]) == selGroup_) { r[1] = "'" + c + "'"; found = true; }
                if (!found) tree.records.push_back({"'" + selGroup_ + "'", "'" + c + "'"});
                X[size_t(treeIdx)].body = di::joinLines(di::writeForm(tree)) + "\n";
                markDirty();
            }
            ImGui::EndCombo();
        }
    }
    // wells of the group: assignment
    int nw = 0;
    for (auto& w : case_.wells) nw += w.group == selGroup_;
    ImGui::TextDisabled("%d well(s). Move the selected well into this group:", nw);
    ImGui::SameLine();
    const bool canMove = selWell_ >= 0 && selWell_ < int(case_.wells.size()) && case_.wells[size_t(selWell_)].group != selGroup_ && selGroup_ != "FIELD";
    ImGui::BeginDisabled(!canMove);
    if (ImGui::SmallButton(canMove ? ("Move " + case_.wells[size_t(selWell_)].name).c_str() : "Move well")) { case_.wells[size_t(selWell_)].group = selGroup_; markDirty(); }
    ImGui::EndDisabled();
    static char newGroup[32] = "";
    ImGui::SetNextItemWidth(fs * 7);
    ImGui::InputTextWithHint("##ng", "new group", newGroup, sizeof(newGroup), ImGuiInputTextFlags_CharsUppercase);
    ImGui::SameLine();
    if (ImGui::SmallButton(TRL("+ Group under this")) && newGroup[0]) {
        const std::string ng = newGroup;
        if (!groups.count(ng) && ng != "FIELD") {
            if (treeIdx < 0) {
                ExtraKeyword e; e.name = "GRUPTREE"; e.section = DeckSection::Schedule; e.step = 0; e.beforeControls = true; e.note = "group tree (ResSim Studio)";
                X.insert(X.begin(), e); treeIdx = 0;
                tree = di::parseForm("GRUPTREE", {});
            }
            tree.records.push_back({"'" + ng + "'", "'" + selGroup_ + "'"});
            X[size_t(treeIdx)].body = di::joinLines(di::writeForm(tree)) + "\n";
            selGroup_ = ng;
            newGroup[0] = 0;
            markDirty();
        }
    }
    // group controls of this group
    ImGui::SeparatorText(TR("Group controls"));
    static const char* kGroupKw[] = {"GCONPROD", "GCONINJE", "GCONSUMP", "GCONSALE", "GEFAC", "GRUPTARG", "GCONPRI", "GUIDERAT"};
    int shown = 0;
    for (int i = 0; i < int(X.size()); ++i) {
        const ExtraKeyword& e = X[size_t(i)];
        bool isG = false;
        for (auto* k : kGroupKw) if (e.name == k) isG = true;
        if (!isG) continue;
        const di::Form f = di::parseForm(e.name, di::splitLines(e.body));
        for (auto& r : f.records) {
            if (r.empty()) continue;
            const std::string g = di::unquote(r[0]);
            if (g != selGroup_ && !(g.find('*') != std::string::npos && selGroup_.rfind(g.substr(0, g.find('*')), 0) == 0)) continue;
            std::string txt;
            for (size_t k = 1; k < r.size() && k < 8; ++k) txt += r[k] + " ";
            ImGui::PushID(i * 1000 + shown);
            ImGui::Text("%s  step %d:", e.name.c_str(), e.step);
            ImGui::SameLine();
            ImGui::TextDisabled("%s", txt.c_str());
            ImGui::SameLine();
            if (ImGui::SmallButton(TRL("Edit"))) { advSel_ = i; advFormView_ = true; section_ = Section::Advanced; }
            ImGui::PopID();
            ++shown;
        }
    }
    if (!shown) ImGui::TextDisabled("%s", TR("none"));
    auto addCtl = [&](const char* kw, const std::string& rec) {
        ExtraKeyword e; e.name = kw; e.section = DeckSection::Schedule; e.step = 0; e.beforeControls = true; e.note = "group control (ResSim Studio)";
        e.body = rec + "\n/\n";
        X.push_back(e);
        advSel_ = int(X.size()) - 1; advFormView_ = true; section_ = Section::Advanced;
        markDirty();
    };
    const bool fieldU = case_.units == UnitSystem::Field;
    if (ImGui::SmallButton(TRL("+ Production target (GCONPROD)"))) addCtl("GCONPROD", " '" + selGroup_ + "' ORAT " + (fieldU ? "5000" : "800") + " 1* 1* 1* RATE /");
    ImGui::SameLine();
    if (ImGui::SmallButton(TRL("+ Water injection (GCONINJE)"))) addCtl("GCONINJE", " '" + selGroup_ + "' WATER VREP 1* 1* 1* 1.0 /");
    if (case_.phases.gas) { ImGui::SameLine(); if (ImGui::SmallButton(TRL("+ Gas injection"))) addCtl("GCONINJE", " '" + selGroup_ + "' GAS REIN 1* 1* 0.8 /"); }
    ImGui::SameLine();
    if (ImGui::SmallButton(TRL("+ Guide-rate formula (GUIDERAT)"))) {
        ExtraKeyword e; e.name = "GUIDERAT"; e.section = DeckSection::Schedule; e.step = 0; e.beforeControls = true; e.note = "guide-rate formula (ResSim Studio)";
        e.body = " 0 OIL 1.0 1.0 1* 1* 1* 1* YES 1.0 /\n";   // GR = POT^A / (B + C R1^D + E R2^F)
        X.push_back(e); advSel_ = int(X.size()) - 1; advFormView_ = true; section_ = Section::Advanced; markDirty();
    }
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", TR("GUIDERAT: guide rate of each well = POT^A / (B + C R1^D + E R2^F) of the nominated phase\n"
                          "OIL: R1 = WOR, R2 = GOR; LIQ: WCT, GLR; GAS: WGR, OGR. Shares group targets among GRUP wells.\n"
                          "Simulated by the compositional core, OPM Flow and ECLIPSE."));
    ImGui::TextColored(col4(Pal::Warn), case_.comp.enabled
        ? "Compositional core: GCONPROD targets (ORAT/GRAT/WRAT/LRAT/RESV), GCONINJE RATE and GUIDERAT are simulated; other options run in OPM Flow / ECLIPSE."
        : "Black-oil engine: GCONPROD / GCONINJE (RATE, RESV, REIN, VREP), GUIDERAT, GLIFTOPT and the surface network are simulated; "
          "other group options run in OPM Flow / ECLIPSE.");
    ImGui::EndChild();
}

void App::editWells() {
    const GridInput& g = case_.grid;
    if (ImGui::CollapsingHeader(TRL("Groups & group controls"))) editGroups();
    auto addWell = [&](WellType t) {
        WellInput w;
        int n = 1;
        const char* prefix = t == WellType::Producer ? "PROD" : t == WellType::WaterInjector ? "WINJ" : "GINJ";
        for (;;) {
            std::string nm = prefix + std::to_string(n);
            bool used = false; for (auto& x : case_.wells) if (x.name == nm) used = true;
            if (!used) { w.name = nm; break; }
            ++n;
        }
        w.type = t; w.i = g.nx / 2; w.j = g.ny / 2;
        Completion cp; cp.i = w.i; cp.j = w.j; cp.k1 = 0; cp.k2 = g.nz - 1;
        w.completions.push_back(cp);
        case_.wells.push_back(w);
        if (case_.schedule.empty()) case_.schedule.push_back(ReportStep{});
        WellControl wc; wc.well = w.name;
        if (t == WellType::Producer) { wc.mode = Control::ORAT; wc.orat = case_.units == UnitSystem::Field ? 1000 : 150; wc.bhp = case_.units == UnitSystem::Field ? 1000 : 70; }
        else { wc.mode = Control::RATE; wc.rate = case_.units == UnitSystem::Field ? 1000 : 150; wc.bhp = case_.units == UnitSystem::Field ? 6000 : 400; }
        case_.schedule[0].controls.push_back(wc);
        selWell_ = int(case_.wells.size()) - 1;
        markDirty();
    };
    if (ImGui::Button(TRL("+ Producer"))) addWell(WellType::Producer);
    ImGui::SameLine();
    if (ImGui::Button(TRL("+ Water injector"))) addWell(WellType::WaterInjector);
    ImGui::SameLine();
    ImGui::BeginDisabled(!case_.phases.gas);
    if (ImGui::Button(TRL("+ Gas injector"))) addWell(WellType::GasInjector);
    ImGui::EndDisabled();
    if (case_.wells.empty()) { ui::hint("No wells yet. Add a producer and an injector, then set their controls in the Schedule."); return; }
    selWell_ = std::max(0, std::min(selWell_, int(case_.wells.size()) - 1));

    // ---- areal map with wells (click to place)
    static bool placeMode = false;
    ui::header("Well locations (top view)");
    ImGui::Checkbox(TRL("Click on the map to move the selected well"), &placeMode);
    if (preview_ && preview_->numCells() == g.size() && ImPlot::BeginPlot("##map", ImVec2(-1, ImGui::GetFontSize() * 16), ImPlotFlags_Equal | ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxes("I", "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
        std::vector<double> m(size_t(g.nx) * g.ny);
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                double s = 0; int cnt = 0;
                for (int k = 0; k < g.nz; ++k) { const int q = g.index(i, j, k); if (preview_->activeOf[q] >= 0) { s += std::log10(std::max(g.permx[q], 1e-3)); ++cnt; } }
                m[size_t(j) * g.nx + i] = cnt ? s / cnt : NAN;
            }
        // inactive columns: NaN must never reach the colormap (it produced a negative table index -> assertion);
        // they get the lowest value here and are painted grey on top
        double lo = 1e300, hi = -1e300;
        for (double v : m) if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
        if (lo > hi) { lo = 0; hi = 1; } else if (hi - lo < 1e-9) hi = lo + 1;
        std::vector<char> inactive(m.size(), 0);
        for (size_t q = 0; q < m.size(); ++q) if (!std::isfinite(m[q])) { m[q] = lo; inactive[q] = 1; }
        ImPlot::PushColormap(ImPlotColormap_Viridis);
        ImPlot::PlotHeatmap("logk", m.data(), g.ny, g.nx, lo, hi, nullptr, ImPlotPoint(0.5, g.ny + 0.5), ImPlotPoint(g.nx + 0.5, 0.5));
        ImPlot::PopColormap();
        {
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            ImPlot::PushPlotClipRect();
            for (int j = 0; j < g.ny; ++j)
                for (int i = 0; i < g.nx; ++i)
                    if (inactive[size_t(j) * g.nx + i]) {
                        pk::rect(i + 0.5, j + 0.5, i + 1.5, j + 1.5, IM_COL32(48, 50, 56, 255));
                    }
            for (const LgrInput& L : g.lgrs) {   // local grid refinements: outline + local cell lines
                const ImVec2 a = ImPlot::PlotToPixels(L.i1 + 0.5, L.j1 + 0.5), b = ImPlot::PlotToPixels(L.i2 + 1.5, L.j2 + 1.5);
                const ImVec2 lo(std::min(a.x, b.x), std::min(a.y, b.y)), hi(std::max(a.x, b.x), std::max(a.y, b.y));
                std::vector<int> h; std::vector<double> f0, f1;
                L.split(0, h, f0, f1);
                for (int m = 1; m < L.nx && m < int(h.size()); ++m) {
                    const double x = h[size_t(m)] + 0.5 + f0[size_t(m)];
                    pk::line(x, L.j1 + 0.5, x, L.j2 + 1.5, IM_COL32(255, 255, 255, 70));
                }
                L.split(1, h, f0, f1);
                for (int m = 1; m < L.ny && m < int(h.size()); ++m) {
                    const double y = h[size_t(m)] + 0.5 + f0[size_t(m)];
                    pk::line(L.i1 + 0.5, y, L.i2 + 1.5, y, IM_COL32(255, 255, 255, 70));
                }
                pk::rect(L.i1 + 0.5, L.j1 + 0.5, L.i2 + 1.5, L.j2 + 1.5, 0, IM_COL32(255, 210, 90, 255), 1.5f);
                pk::text(std::min(a.x, b.x) == a.x ? L.i1 + 0.5 : L.i2 + 1.5, std::min(a.y, b.y) == a.y ? L.j1 + 0.5 : L.j2 + 1.5, L.name.c_str(),
                         IM_COL32(255, 210, 90, 255), ImVec2(2, ImGui::GetFontSize() * 0.5f + 1), 0);
                (void)lo; (void)hi; (void)dl;
            }
            ImPlot::PopPlotClipRect();
        }
        for (int w = 0; w < int(case_.wells.size()); ++w) {
            const auto& wl = case_.wells[w];
            double x = wl.i + 1, y = wl.j + 1;
            if (wl.lgr >= 0 && wl.lgr < int(g.lgrs.size())) {   // local location -> position inside the host box
                const LgrInput& L = g.lgrs[size_t(wl.lgr)];
                std::vector<int> h; std::vector<double> f0, f1;
                L.split(0, h, f0, f1); if (wl.i < int(h.size())) x = h[size_t(wl.i)] + 0.5 + 0.5 * (f0[size_t(wl.i)] + f1[size_t(wl.i)]);
                L.split(1, h, f0, f1); if (wl.j < int(h.size())) y = h[size_t(wl.j)] + 0.5 + 0.5 * (f0[size_t(wl.j)] + f1[size_t(wl.j)]);
            }
            ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = w == selWell_ ? 7.0f : 5.0f;
            s.MarkerFillColor = col4(wl.type == WellType::Producer ? Pal::Oil : wl.type == WellType::WaterInjector ? Pal::Water : Pal::Gas);
            s.MarkerLineColor = w == selWell_ ? ImVec4(1, 1, 1, 1) : ImVec4(0, 0, 0, 1);
            s.LineWeight = 2.0f;
            ImPlot::PlotScatter(wl.name.c_str(), &x, &y, 1, s);
            ImPlot::PlotText(wl.name.c_str(), x, y, ImVec2(0, -14));
        }
        if (placeMode && ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) {
            const ImPlotPoint mp = ImPlot::GetPlotMousePos();
            const int i = std::max(0, std::min(g.nx - 1, int(std::lround(mp.x)) - 1));
            const int j = std::max(0, std::min(g.ny - 1, int(std::lround(mp.y)) - 1));
            WellInput& w = case_.wells[selWell_];
            if (w.lgr >= 0) notify("Well " + w.name + " is in a local grid: edit its local location, or move it to the global grid first");
            else {
                w.i = i; w.j = j;
                for (auto& cp : w.completions) { cp.i = i; cp.j = j; }
                markDirty();
            }
        }
        ImPlot::EndPlot();
    }

    // ---- list + details
    ImGui::BeginChild("wlist", ImVec2(ImGui::GetFontSize() * 8, ImGui::GetFontSize() * 11.5f), ImGuiChildFlags_Borders);
    for (int w = 0; w < int(case_.wells.size()); ++w) {
        const auto& wl = case_.wells[w];
        const ImU32 c = wl.type == WellType::Producer ? Pal::Oil : wl.type == WellType::WaterInjector ? Pal::Water : Pal::Gas;
        ImGui::PushStyleColor(ImGuiCol_Text, col4(c));
        ImGui::TextUnformatted("\xE2\x97\x8F"); ImGui::PopStyleColor(); ImGui::SameLine();
        if (ImGui::Selectable(wl.name.c_str(), selWell_ == w)) selWell_ = w;
    }
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("wdet", ImVec2(0, ImGui::GetFontSize() * 11.5f));
    WellInput& w = case_.wells[selWell_];
    if (ui::beginProps("well")) {
        ui::label("Name");
        std::string nm = w.name;
        if (ui::inputString("##name", &nm, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsNoBlank)) {
            bool clash = false; for (auto& x : case_.wells) if (&x != &w && x.name == nm) clash = true;
            if (!nm.empty() && !clash) {
                for (auto& s : case_.schedule) for (auto& c : s.controls) if (c.well == w.name) c.well = nm;
                w.name = nm; markDirty();
            }
        }
        static const char* types[] = {"Producer", "Water injector", "Gas injector"};
        int t = int(w.type);
        if (ui::propCombo("Type", &t, types, 3)) {
            w.type = WellType(t);
            for (auto& s : case_.schedule) for (auto& c : s.controls) if (c.well == w.name) {
                if (w.type == WellType::Producer && c.mode == Control::RATE) { c.mode = Control::ORAT; c.orat = std::max(c.rate, 0.0); }
                if (w.type != WellType::Producer && c.mode != Control::RATE && c.mode != Control::BHP) { c.mode = Control::RATE; c.rate = std::max(c.orat, 0.0); }
            }
            markDirty();
        }
        if (!g.lgrs.empty() || w.lgr >= 0) {
            std::vector<std::string> gn{"Global grid (WELSPECS / COMPDAT)"};
            for (auto& L : g.lgrs) gn.push_back("LGR " + L.name + " (WELSPECL / COMPDATL)");
            std::vector<const char*> gp; for (auto& x : gn) gp.push_back(x.c_str());
            int gl = std::max(-1, std::min(int(g.lgrs.size()) - 1, w.lgr)) + 1;
            if (ui::propCombo("Grid", &gl, gp.data(), int(gp.size()),
                              "Grid the well and its connections live in. Moving a well into an LGR places each connection in the middle "
                              "child column of its host column and in every refined layer of its host layer (explicit CF split); "
                              "moving it back maps local cells to their host cells.")) {
                try { moveWellToLgr(case_, selWell_, gl - 1); markDirty(); }
                catch (const std::exception& e) { notify(std::string("LGR: ") + e.what()); }
            }
        }
        const LgrInput* WL = (w.lgr >= 0 && w.lgr < int(g.lgrs.size())) ? &g.lgrs[size_t(w.lgr)] : nullptr;
        const int wnx = WL ? WL->nx : g.nx, wny = WL ? WL->ny : g.ny;
        int ij[2] = {w.i + 1, w.j + 1};
        ui::label(WL ? "Local location I, J" : "Location I, J");
        if (ImGui::InputInt2("##ij", ij)) {
            w.i = std::max(0, std::min(wnx - 1, ij[0] - 1)); w.j = std::max(0, std::min(wny - 1, ij[1] - 1));
            for (auto& cp : w.completions) { cp.i = w.i; cp.j = w.j; }
            markDirty();
        }
        if (ui::propD("PI multiplier (WPIMULT)", &w.piMult, "%.4g", "Multiplies all connection factors of the well for the whole run (stimulation / damage, history matching). Deck WPIMULT records are applied per connection and step: column 'PI mult' of the completions.")) {
            w.piMult = std::max(0.0, w.piMult); markDirty();
        }
        bool autoRef = w.refDepth <= 0;
        ui::label("BHP reference depth");
        if (ImGui::Checkbox(TRL("Top perforation"), &autoRef)) { w.refDepth = autoRef ? -1 : case_.equil.datumDepth; markDirty(); }
        if (!autoRef) { ImGui::SameLine(); ImGui::SetNextItemWidth(-FLT_MIN); if (ImGui::InputDouble("##ref", &w.refDepth, 0, 0, "%.6g")) markDirty(); }
        ui::label("WELSPECS items 7+ (OPM)", "Written after the preferred phase in WELSPECS, for OPM Flow:\n"
                  "drainage radius, inflow equation, auto shut-in (SHUT/STOP), crossflow (YES/NO), PVT table, density calc (SEG/AVG), ...\n"
                  "e.g.  0.0 STD SHUT NO 1* SEG   - not used by the ResSim engine.");
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ui::inputString("##wtail", w.specsTail)) markDirty();
        ui::endProps();
    }
    ImGui::EndChild();
    ui::header("Completions");
    ImGui::PushTextWrapPos(0);
    ImGui::TextDisabled("%s", TR("Layers K1..K2 (1-based). CF <= 0: Peaceman well index from diameter, skin and cell properties. "
                        "From / Until: report step (1-based) at which the connection is perforated / shut (late COMPDAT); "
                        "Until 0 = never."));
    ImGui::PopTextWrapPos();
    const int nsteps = std::max(1, int(case_.schedule.size()));
    static bool opmCols = false;
    bool anyOpm = false;
    for (auto& q : w.completions) if (q.satTable > 0 || q.dfactor >= 0 || q.r0 > 0) anyOpm = true;
    ImGui::Checkbox(TRL("OPM Flow connection items (sat. table, D-factor, r0)"), &opmCols);
    if (anyOpm && !opmCols) { ImGui::SameLine(); ImGui::TextColored(col4(Pal::Warn), "(set for this well)"); }
    {   // C15: spreadsheet table (multi-cell selection, copy / paste from Excel or COMPDAT columns, row insert / delete)
        const float fsz = ImGui::GetFontSize();
        auto dimsOf = [&](const Completion& cp, int& nxl, int& nyl, int& nzl) {
            if (cp.lgr >= 0 && cp.lgr < int(g.lgrs.size())) { const LgrInput& L = g.lgrs[size_t(cp.lgr)]; nxl = L.nx; nyl = L.ny; nzl = L.nz; }
            else { nxl = g.nx; nyl = g.ny; nzl = g.nz; }
        };
        // columns: 0 I, 1 J, 2 K1, 3 K2, 4 skin, 5 diameter, 6 CF, 7 Dir, 8 Open, 9 From, 10 Until, [11 Sat tbl, 12 D-factor, 13 r0]
        ui::GridAccess ga;
        ga.count = [&] { return int(w.completions.size()); };
        ga.cols = opmCols ? 15 : 12;
        ga.getText = [&](int r, int c) -> std::string {
            const Completion& cp = w.completions[size_t(r)];
            if (c == 7) return std::string(1, cp.dir);
            if (c == 8) return cp.open ? "OPEN" : "SHUT";
            return std::string();
        };
        ga.setText = [&](int r, int c, const std::string& v) {
            Completion& cp = w.completions[size_t(r)];
            if (v.empty()) return;
            const char ch0 = char(std::toupper(static_cast<unsigned char>(v[0])));
            if (c == 7 && (ch0 == 'X' || ch0 == 'Y' || ch0 == 'Z')) cp.dir = ch0;
            if (c == 8) cp.open = !(ch0 == 'S' || ch0 == 'N' || ch0 == '0' || ch0 == 'F');
        };
        ga.numGet = [&](int r, int c, double& v) -> bool {
            const Completion& cp = w.completions[size_t(r)];
            switch (c) {
                case 0: v = cp.i + 1; return true; case 1: v = cp.j + 1; return true;
                case 2: v = cp.k1 + 1; return true; case 3: v = cp.k2 + 1; return true;
                case 4: v = cp.skin; return true; case 5: v = cp.diameter; return true; case 6: v = cp.cf; return true;
                case 9: v = cp.fromStep + 1; return true; case 10: v = cp.untilStep < 0 ? 0 : cp.untilStep + 1; return true;
                case 11: v = cp.mult; return true;
                case 12: v = cp.satTable; return true; case 13: v = cp.dfactor; return true; case 14: v = cp.r0; return true;
                default: return false;
            }
        };
        ga.numSet = [&](int r, int c, double v) {
            Completion& cp = w.completions[size_t(r)];
            int nxl, nyl, nzl; dimsOf(cp, nxl, nyl, nzl);
            const int iv = int(std::lround(v));
            switch (c) {
                case 0: cp.i = std::clamp(iv, 1, nxl) - 1; break;
                case 1: cp.j = std::clamp(iv, 1, nyl) - 1; break;
                case 2: cp.k1 = std::clamp(iv, 1, nzl) - 1; cp.k2 = std::max(cp.k1, cp.k2); break;
                case 3: cp.k2 = std::max(cp.k1, std::clamp(iv, 1, nzl) - 1); break;
                case 4: cp.skin = v; break;
                case 5: if (v > 0) cp.diameter = v; break;
                case 6: cp.cf = v; break;
                case 9: cp.fromStep = std::clamp(iv, 1, nsteps) - 1; break;
                case 10: cp.untilStep = iv <= 0 ? -1 : std::max(cp.fromStep + 1, std::min(nsteps, iv - 1)); break;
                case 11: if (v >= 0) cp.mult = v; break;
                case 12: cp.satTable = std::max(0, iv); break;
                case 13: cp.dfactor = v; break;
                case 14: cp.r0 = v; break;
            }
        };
        ga.insert = [&](int at, int n) {
            at = std::clamp(at, 0, int(w.completions.size()));
            Completion cp; cp.i = w.i; cp.j = w.j; cp.lgr = w.lgr;
            cp.k2 = (w.lgr >= 0 && w.lgr < int(g.lgrs.size())) ? g.lgrs[size_t(w.lgr)].nz - 1 : g.nz - 1;
            if (!w.completions.empty()) cp = w.completions[size_t(std::min(at, int(w.completions.size()) - 1))];
            w.completions.insert(w.completions.begin() + at, size_t(n), cp);
        };
        ga.erase = [&](int r0, int r1) {
            r1 = std::min(r1, int(w.completions.size()) - 1);
            if (r0 <= r1) w.completions.erase(w.completions.begin() + r0, w.completions.begin() + r1 + 1);
        };
        std::vector<const char*> hd = {"I", "J", "K1", "K2", "Skin", "Diameter", "CF (<=0 Peaceman)", "Dir (X/Y/Z)", "Open / Shut",
                                       "From step", "Until step (0 never)", "PI mult (WPIMULT)"};
        if (opmCols) { hd.push_back("Sat tbl"); hd.push_back("D-factor"); hd.push_back("r0"); }
        const float h = std::min(fsz * 16, ImGui::GetFrameHeight() * (float(w.completions.size()) + 2.6f) + 4);
        if (ui::gridTable("comp", ga, hd, h, nullptr, "%.4g", true)) markDirty();
        if (w.lgr >= 0) ui::hint("I / J / K are local cell indices of the LGR.");
    }
    if (ImGui::SmallButton(TRL("+ Add completion"))) {
        Completion cp; cp.i = w.i; cp.j = w.j; cp.k1 = 0; cp.lgr = w.lgr;
        cp.k2 = (w.lgr >= 0 && w.lgr < int(g.lgrs.size())) ? g.lgrs[size_t(w.lgr)].nz - 1 : g.nz - 1;
        w.completions.push_back(cp); markDirty();
    }
    ImGui::SameLine();
    if (ui::dangerButton("Delete well")) {
        const std::string nm = w.name;
        for (auto& s : case_.schedule)
            s.controls.erase(std::remove_if(s.controls.begin(), s.controls.end(), [&](const WellControl& c) { return c.well == nm; }), s.controls.end());
        case_.wells.erase(case_.wells.begin() + selWell_);
        selWell_ = std::max(0, selWell_ - 1);
        markDirty();
        return;
    }
    // C32: multisegment well definition from the completions (WELSEGS / COMPSEGS / devices)
    if (ImGui::CollapsingHeader(TRL("Multisegment well (WELSEGS / COMPSEGS / ICDs)"))) mswGenerator(selWell_);
    if (ImGui::CollapsingHeader(TRL("Gas-lift optimisation (LIFTOPT / WLIFTOPT / GLIFTOPT)"))) editGasLift();   // C35
}

// ============================================================================ SCHEDULE
void App::editSchedule() {
    static double totalYears = 10, stepDays = 91.3125;
    if (ImGui::CollapsingHeader(TRL("Run-time logic: UDQ, UDA, ACTIONX"))) editRunLogic();   // C31
    ui::header("Report steps");
    double T = 0; for (auto& s : case_.schedule) T += s.length;
    ImGui::Text("%d report steps, total %.0f days (%.2f years)", int(case_.schedule.size()), T, T / 365.25);
    if (ImGui::TreeNode(TRL("Rebuild uniform steps"))) {
        if (ui::beginProps("rb")) {
            ui::propD("Total time (years)", &totalYears, "%.4g");
            ui::propD("Step length (days)", &stepDays, "%.4g");
            ui::endProps();
        }
        if (ImGui::Button(TRL("Rebuild (keeps well controls at their times)"))) {
            const int n = std::max(1, int(std::ceil(totalYears * 365.25 / std::max(1e-3, stepDays) - 1e-9)));
            std::vector<std::pair<double, WellControl>> ctrls;
            double t = 0;
            for (auto& s : case_.schedule) { for (auto& c : s.controls) ctrls.push_back({t, c}); t += s.length; }
            std::vector<ReportStep> ns(n);
            for (auto& s : ns) s.length = stepDays;
            for (auto& [tc, c] : ctrls) { int k = std::min(n - 1, int(std::floor(tc / stepDays + 1e-6))); ns[k].controls.push_back(c); }
            case_.schedule = ns;
            selStep_ = 0;
            markDirty();
        }
        ImGui::TreePop();
    }
    if (case_.schedule.empty()) case_.schedule.push_back(ReportStep{});
    selStep_ = std::max(0, std::min(selStep_, int(case_.schedule.size()) - 1));

    // timeline table
    {   // spreadsheet grid of the report steps: paste a column of step lengths (TSTEP) from Excel
        std::vector<double> start(case_.schedule.size() + 1, 0.0);
        for (size_t k = 0; k < case_.schedule.size(); ++k) start[k + 1] = start[k] + case_.schedule[k].length;
        ui::GridAccess ga;
        ga.count = [&] { return int(case_.schedule.size()); };
        ga.cols = 3;
        ga.num = [&](int r, int c) -> double* { return c == 1 ? &case_.schedule[size_t(r)].length : nullptr; };
        ga.getText = [&](int r, int c) -> std::string {
            if (c == 0) { char b[32]; std::snprintf(b, sizeof b, "%.2f", r < int(start.size()) ? start[size_t(r)] : 0.0); return b; }
            std::string l;
            for (auto& cc : case_.schedule[size_t(r)].controls) { if (!l.empty()) l += ", "; l += cc.well; }
            return l;
        };
        ga.insert = [&](int at, int n) {
            at = std::max(0, std::min(at, int(case_.schedule.size())));
            ReportStep r; r.length = case_.schedule.empty() ? 30.0 : case_.schedule[size_t(std::min(at, int(case_.schedule.size()) - 1))].length;
            case_.schedule.insert(case_.schedule.begin() + at, size_t(n), r);
        };
        ga.erase = [&](int r0, int r1) {
            r1 = std::min(r1, int(case_.schedule.size()) - 1);
            if (r0 > r1 || int(case_.schedule.size()) - (r1 - r0 + 1) < 1) return;   // keep at least one step
            std::vector<WellControl> moved;
            for (int k = r0; k <= r1; ++k) for (auto& c : case_.schedule[size_t(k)].controls) moved.push_back(c);
            case_.schedule.erase(case_.schedule.begin() + r0, case_.schedule.begin() + r1 + 1);
            const int dst = std::max(0, std::min(r0, int(case_.schedule.size()) - 1));
            for (auto& c : moved) case_.schedule[size_t(dst)].controls.push_back(c);   // control changes are kept
        };
        if (ui::gridTable("steps", ga, {"Start (days)", "Length (days)", "Control changes"}, ImGui::GetFontSize() * 13, &selStep_, "%.4g", false)) {
            for (auto& r : case_.schedule) r.length = std::max(1e-3, r.length);
            markDirty();
        }
        selStep_ = std::max(0, std::min(selStep_, int(case_.schedule.size()) - 1));
    }
    ReportStep& st = case_.schedule[selStep_];
    ui::header(("Step " + std::to_string(selStep_ + 1)).c_str());
    if (ui::beginProps("st")) {
        if (ui::propD("Length (days)", &st.length, "%.6g")) { st.length = std::max(1e-3, st.length); markDirty(); }
        ui::endProps();
    }
    ImGui::SameLine();
    ImGui::TextDisabled("%s", TR("Controls set here stay in force until changed in a later step."));
    if (ImGui::SmallButton(TRL("Insert step after"))) { case_.schedule.insert(case_.schedule.begin() + selStep_ + 1, [&] { ReportStep r; r.length = st.length; return r; }()); markDirty(); }
    ImGui::SameLine();
    if (case_.schedule.size() > 1 && ImGui::SmallButton(TRL("Delete step"))) {
        auto ctrls = st.controls;
        case_.schedule.erase(case_.schedule.begin() + selStep_);
        selStep_ = std::max(0, selStep_ - 1);
        for (auto& c : ctrls) case_.schedule[selStep_].controls.push_back(c);
        markDirty();
        return;
    }

    ui::header("Well controls in this step");
    int del = -1;
    for (int ci = 0; ci < int(st.controls.size()); ++ci) {
        WellControl& c = st.controls[ci];
        ImGui::PushID(ci);
        const WellInput* w = nullptr; for (auto& x : case_.wells) if (x.name == c.well) w = &x;
        const bool prod = w ? w->type == WellType::Producer : c.mode != Control::RATE;
        ImGui::BeginChild("ctl", ImVec2(0, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_AutoResizeY);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
        if (ImGui::BeginCombo("##well", c.well.c_str())) {
            for (auto& x : case_.wells) if (ImGui::Selectable(x.name.c_str(), x.name == c.well)) { c.well = x.name; markDirty(); }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        static const char* stats[] = {"OPEN", "SHUT", "STOP"};
        int sidx = int(c.status);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5);
        if (ImGui::Combo("##status", &sidx, stats, 3)) { c.status = WellStatus(sidx); markDirty(); }
        ImGui::SameLine();
        static const char* pm[] = {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP"};
        static const char* im[] = {"RATE", "BHP"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5);
        if (prod) {
            int m = std::min(int(c.mode), 5);
            if (ImGui::Combo("##mode", &m, pm, 6)) { c.mode = Control(m); markDirty(); }
        } else {
            int m = c.mode == Control::BHP ? 1 : 0;
            if (ImGui::Combo("##mode", &m, im, 2)) { c.mode = m ? Control::BHP : Control::RATE; markDirty(); }
        }
        ImGui::SameLine();
        if (prod) { if (ImGui::Checkbox(TRL("Hist"), &c.history)) markDirty(); if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("History mode (WCONHIST): rates are observed data")); ImGui::SameLine(); }
        if (ImGui::SmallButton(TRL("\xC3\x97"))) del = ci;
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Remove this control"));
        // targets
        const Units U = Units::of(case_.units);
        char lb[64];
        // targets/limits in a 2 x (label, value) grid that fits the panel width
        auto fieldD = [&](const char* label, const char* id, double* v, bool active) {
            ImGui::TableNextColumn();
            if (active) ImGui::TextColored(col4(Pal::Accent), "%s", label); else ImGui::TextDisabled("%s", label);
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (active) ImGui::PushStyleColor(ImGuiCol_FrameBg, themeSelectedField());
            const bool chg = ImGui::InputDouble(id, v, 0, 0, "%.6g");
            if (active) ImGui::PopStyleColor();
            if (chg) markDirty();
            if (ImGui::IsItemHovered()) ImGui::SetTooltip(active ? "Target of the active control mode" : "Limit (negative = none)");
        };
        if (ui::beginTable("targets", 4, ImGuiTableFlags_SizingStretchProp)) {
            ui::tableColumn("l1", ImGuiTableColumnFlags_WidthStretch, 0.8f);
            ui::tableColumn("v1", ImGuiTableColumnFlags_WidthStretch, 1.0f);
            ui::tableColumn("l2", ImGuiTableColumnFlags_WidthStretch, 0.8f);
            ui::tableColumn("v2", ImGuiTableColumnFlags_WidthStretch, 1.0f);
            char lb2[64];
            if (prod) {
                std::snprintf(lb, sizeof(lb), "Oil (%s/d)", U.liq);   fieldD(lb, "##o", &c.orat, c.mode == Control::ORAT);
                std::snprintf(lb2, sizeof(lb2), "Water (%s/d)", U.liq); fieldD(lb2, "##w", &c.wrat, c.mode == Control::WRAT);
                std::snprintf(lb, sizeof(lb), "Gas (%s/d)", U.gas);   fieldD(lb, "##g", &c.grat, c.mode == Control::GRAT);
                std::snprintf(lb2, sizeof(lb2), "Liquid (%s/d)", U.liq); fieldD(lb2, "##l", &c.lrat, c.mode == Control::LRAT);
                fieldD("Voidage (RESV)", "##r", &c.resv, c.mode == Control::RESV);
                std::snprintf(lb, sizeof(lb), "BHP (%s)", U.p);       fieldD(lb, "##b", &c.bhp, c.mode == Control::BHP);
            } else {
                std::snprintf(lb, sizeof(lb), "Rate (%s/d)", w && w->type == WellType::GasInjector ? U.gas : U.liq);
                fieldD(lb, "##rate", &c.rate, c.mode == Control::RATE);
                std::snprintf(lb2, sizeof(lb2), "BHP (%s)", U.p);    fieldD(lb2, "##b", &c.bhp, c.mode == Control::BHP);
            }
            ui::endTable();
        }
        // ---- OPM Flow-only control: group / THP control, VFP table, artificial lift
        {
            const bool opmSet = !c.deckMode.empty() || c.thp > 0 || c.vfp > 0 || c.alq > 0;
            if (opmSet) ImGui::SetNextItemOpen(true, ImGuiCond_Once);
            if (ImGui::TreeNodeEx(TRL("OPM Flow control (group / THP / VFP)"), ImGuiTreeNodeFlags_SpanAvailWidth)) {
                static const char* prodModes[] = {"(same as above)", "GRUP", "THP", "NONE", "CRAT"};
                static const char* injModes[]  = {"(same as above)", "GRUP", "THP", "RESV"};
                const char** modes = prod ? prodModes : injModes;
                const int nModes = prod ? 5 : 4;
                const std::string dm = c.deckMode == "RESV*" ? std::string("RESV") : c.deckMode;
                int cur = 0; for (int m = 1; m < nModes; ++m) if (dm == modes[m]) cur = m;
                ImGui::SetNextItemWidth(170);
                if (ImGui::Combo("Mode sent to OPM Flow", &cur, modes, nModes)) {
                    c.deckMode = cur == 0 ? std::string() : std::string(modes[cur]);
                    if (c.deckMode == "RESV" && c.rate <= 0) c.deckMode = "RESV*";
                    markDirty();
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s", TR("OPM Flow runs this control mode (GRUP = under its group's GCONPROD/GCONINJE, THP = tubing-head\n"
                                      "pressure through the VFP table). The ResSim engine does not model groups or VFP and runs\n"
                                      "the mode chosen above, using the targets/limits you give."));
                if (ui::beginTable("opmctl", 4, ImGuiTableFlags_SizingStretchProp)) {
                    ui::tableColumn("l1", ImGuiTableColumnFlags_WidthStretch, 0.8f);
                    ui::tableColumn("v1", ImGuiTableColumnFlags_WidthStretch, 1.0f);
                    ui::tableColumn("l2", ImGuiTableColumnFlags_WidthStretch, 0.8f);
                    ui::tableColumn("v2", ImGuiTableColumnFlags_WidthStretch, 1.0f);
                    char l3[64]; std::snprintf(l3, sizeof(l3), "THP (%s)", U.p);
                    fieldD(l3, "##thp", &c.thp, dm == "THP");
                    ImGui::TableNextColumn(); ImGui::TextDisabled("%s", TR("VFP table"));
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::InputInt("##vfp", &c.vfp, 0)) { c.vfp = std::max(0, c.vfp); markDirty(); }
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("VFPPROD / VFPINJ table number (0 = none). Define the tables on the Advanced keywords page."));
                    if (prod) fieldD("Artificial lift (ALQ)", "##alq", &c.alq, false);
                    ui::endTable();
                }
                if ((dm == "THP" || c.thp > 0) && c.vfp <= 0)
                    ImGui::TextColored(col4(Pal::Warn), "A THP control/limit needs a VFP table number.");
                ImGui::TreePop();
            }
        }
        ImGui::EndChild();
        ImGui::PopID();
    }
    if (del >= 0) { st.controls.erase(st.controls.begin() + del); markDirty(); }
    ImGui::BeginDisabled(case_.wells.empty());
    if (ImGui::Button(TRL("+ Add well control"))) {
        WellControl c; c.well = case_.wells[std::min(selWell_, int(case_.wells.size()) - 1)].name;
        const WellInput* w = nullptr; for (auto& x : case_.wells) if (x.name == c.well) w = &x;
        if (w && w->type != WellType::Producer) { c.mode = Control::RATE; c.rate = 0; } else { c.mode = Control::ORAT; c.orat = 0; }
        st.controls.push_back(c);
        markDirty();
    }
    ImGui::EndDisabled();
    ImGui::Separator();
    editProductionData();   // C38
}

// ============================================================================ NUMERICS
void App::editNumerics() {
    Tuning& t = case_.tuning;
    bool ch = false;
    ui::header("Time stepping");
    if (ui::beginProps("ts")) {
        ch |= ui::propD("Initial step (days)", &t.initialDt, "%.4g");
        ch |= ui::propD("Maximum step (days)", &t.maxDt, "%.4g");
        ch |= ui::propD("Minimum step (days)", &t.minDt, "%.3g");
        ch |= ui::propD("Growth factor", &t.growth, "%.3g");
        ch |= ui::propD("Cut factor after failure", &t.cutFactor, "%.3g");
        ch |= ui::propD("Max saturation change / iter", &t.dsMax, "%.3g", "Appleyard-type chop of Newton updates");
        ch |= ui::propD("Max relative pressure change / iter", &t.dpMaxRel, "%.3g");
        ui::endProps();
    }
    ui::header("Newton and linear solver");
    if (ui::beginProps("nl")) {
        ch |= ui::propI("Max Newton iterations", &t.maxNewton, 2, 50);
        ch |= ui::propD("CNV tolerance (max local)", &t.cnvTol, "%.3g");
        ch |= ui::propD("Material balance tolerance", &t.mbTol, "%.3g");
        ch |= ui::propI("Max linear iterations", &t.maxLinear, 10, 5000);
        ch |= ui::propD("Linear relative tolerance", &t.linTol, "%.3g");
        static const char* ls[] = {"ILU(0) - BiCGStab", "CPR-AMG (two-stage)", "Automatic (CPR for >= 500 cells)"};
        ch |= ui::propCombo("Preconditioner", &t.linearSolver, ls, 3,
                            "CPR: constrained pressure residual - algebraic multigrid on the pressure system, then ILU(0) on the full system. "
                            "Linear iterations become nearly independent of grid size.");
        static const char* cw[] = {"Quasi-IMPES", "True-IMPES"};
        ch |= ui::propCombo("CPR pressure weights", &t.cprWeights, cw, 2,
                            "How CPR forms the pressure equation. True-IMPES (default, as in OPM Flow) combines the mass balances with "
                            "weights from the accumulation terms only - more robust for strongly coupled or gravity-dominated flow. "
                            "Quasi-IMPES uses the full diagonal block of the Jacobian.");
        const int hw = int(std::max(1u, std::thread::hardware_concurrency()));
        ch |= ui::propI("CPU threads (0 = all)", &t.threads, 0, hw);
        ch |= ui::propBool("Mixed precision (float32 BILU0)", &t.mixedPrecision,
                           "The ILU(0) factors (standalone or CPR second stage) are stored and applied in single precision: half the "
                           "memory traffic of the triangular sweeps. The Krylov method stays in double precision, so the converged "
                           "solution is the same; for extremely ill-conditioned systems the double ILU(0) takes over automatically.");
        {   // C45 low-memory linear algebra (TUNING item 10 / SOLVER 4th token)
            static const char* lm[] = {"Auto (500k cells and up)", "Off", "On"};
            int li = std::clamp(t.leanMemory + 1, 0, 2);
            if (ui::propCombo("Low-memory linear solver", &li, lm, 3,
                              "C45: the BILU0 second stage is factorised row by row straight into float32 storage (no double copy), the "
                              "AMG takes over the fine matrix instead of copying it and the accumulation scratch is released. About 30 % "
                              "less peak memory on large models (1 M cells fit in ~2.8 GB). The Krylov method stays double; if the float "
                              "factors stall, a double factorisation takes over automatically. Auto switches it on from 500k active cells.")) {
                t.leanMemory = li - 1; ch = true;
            }
        }
        ch |= ui::propI("AMG reuse across time steps (0 = off)", &t.amgReuse, 0, 50);
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
            ImGui::SetTooltip("C26: keep the CPR-AMG hierarchy (aggregates, transfer operators) for up to N time steps and only refresh the\n"
                              "coarse-level values (Galerkin product) - saves most of the AMG set-up. A quality monitor rebuilds it early\n"
                              "when the linear iterations rise 50 %% above the count right after a fresh build. Converged solutions are\n"
                              "unchanged (same tolerances). Since the C26 restart fix it is limited to one report step (Egg -3 %%).");
        ch |= ui::propD("Target saturation change / step (0 = off)", &t.dsTarget, "%.3g",
                        "C27 (CMG NORM SATUR / Aziz-Settari): the next time step is scaled by 2*target / (change + target), where change is\n"
                        "the largest saturation change of the last step - steps grow while the fronts move slowly and shrink before\n"
                        "Newton has to cut them. Typical 0.05 - 0.1. 0 = iteration-count rule only.");
        ch |= ui::propD("Target pressure change / step (0 = off)", &t.dpTarget, "%.4g",
                        "C27 (CMG NORM PRESS): same control for the largest pressure change of a step, in case pressure units\n"
                        "(e.g. 50 bar / 700 psi). Both targets may be combined; the stricter one wins.");
        if (t.dsTarget < 0) t.dsTarget = 0;
        if (t.dpTarget < 0) t.dpTarget = 0;
        bool g = t.gpu != 0;
        if (ui::propBool("GPU linear solver", &g, "Solve the reduced block system on the GPU (Settings -> GPU acceleration must be on): "
                         "float32 BiCGSTAB with block-Jacobi on the device inside double-precision iterative refinement. Falls back "
                         "to the CPU solver per system if the device is unavailable or does not converge.")) { t.gpu = g ? 1 : 0; ch = true; }
        ui::endProps();
    }
    if (ch) markDirty();
    if (ImGui::Button(TRL("Restore defaults"))) { t = Tuning{}; markDirty(); }
    ui::hint("Fully-implicit black-oil: Newton-Raphson with automatic-differentiation Jacobians, BiCGSTAB/GMRES with CPR-AMG or ILU(0), "
             "OpenMP-parallel assembly and sparse kernels. Convergence: CNV (max scaled residual) and material balance per component.");
}
