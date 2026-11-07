// Advanced keywords editor: every ECLIPSE / OPM Flow keyword the ResSim engine does not model natively (groups,
// VFP tables, multisegment wells, ACTIONX / UDQ, tracers, extra summary requests, ...) is kept verbatim in
// Case::extraKeywords. This page lists, edits, adds, re-orders and checks them; the deck exporter writes them in
// place, so OPM Flow runs the full case.
#include "app.h"
#include "widgets.h"
#include "theme.h"
#include "ressim/deck_parser.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <vector>

using namespace ressim;

namespace {

const char* kSecNames[] = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};

std::string upperTrim(const std::string& s) {
    std::string r;
    for (char ch : s) if (!std::isspace((unsigned char)ch)) r += char(std::toupper((unsigned char)ch));
    return r;
}

// Parses "SECTION / NAME / body" with the generic deck reader and returns its messages (empty = clean)
std::vector<std::string> checkKeyword(const ExtraKeyword& e) {
    std::vector<std::string> out;
    if (e.name.empty()) { out.push_back("Error: keyword name is empty"); return out; }
    if (!deck::findKeyword(e.name)) out.push_back("Warning: " + e.name + " is not in the OPM Flow / ECLIPSE keyword index (typing error?)");
    std::string text = std::string(kSecNames[int(e.section)]) + "\n" + e.name + "\n" + e.body;
    if (!text.empty() && text.back() != '\n') text += '\n';
    try {
        const deck::RawDeck D = deck::parseDeckText(text, "", "keyword");
        for (auto& m : D.messages) {
            if (m.level == deck::Level::Info) continue;
            out.push_back(std::string(m.level == deck::Level::Error ? "Error" : "Warning") + " (line " + std::to_string(std::max(0, m.line - 2)) + "): " + m.text);
        }
        int n = 0; for (auto& k : D.keywords) if (k.name == e.name) ++n;
        if (n == 0) out.push_back("Error: the keyword was not recognised");
        else if (D.keywords.size() > 2) out.push_back("Warning: the text contains more than one keyword - each one is written as is");
    } catch (const std::exception& ex) { out.push_back(std::string("Error: ") + ex.what()); }
    return out;
}

}  // namespace

void App::editAdvanced() {
    auto& X = case_.extraKeywords;
    const float fs = ImGui::GetFontSize();
    ui::hint("Run with OPM Flow to simulate them. SCHEDULE keywords are written at their report step. "
             "The ResSim black-oil engine also simulates WECON / WTEST (economic limits, well tests) and TRACER / WTRACER / "
             "TBLKF<name> (passive tracers, see Tracers below) from these records.");

    // ---- C27 passive tracers: TRACER / WTRACER stay verbatim records below; decay, adsorption and the initial value are
    // the ResSim extension RESSIM TRACERPROP (case_.tracerProps)
    if (ImGui::CollapsingHeader("Tracers (ResSim engine + OPM Flow)")) {
        std::vector<std::pair<std::string, std::string>> trs;   // name, phase
        for (const auto& k : X) if (k.name == "TRACER")
            for (const auto& r : verbatimRecords(k.body)) if (r.size() >= 2 && !r[0].empty()) trs.push_back({upperTrim(r[0]), upperTrim(r[1])});
        ui::hint("Interwell tracer tests / water-source tracking: a tracer rides with the free phase it is declared in. TRACER (PROPS) "
                 "declares it, WTRACER (SCHEDULE) sets the injected concentration per well, TBLKF<name> (SOLUTION) optional per-cell "
                 "initial values. Half-life (first-order decay, TOUGH3 EOS7R-like) and linear adsorption kd (adsorbed amount per pore "
                 "volume = kd c, retardation 1 + kd / (b S)) are ResSim extensions. C28: partition coefficient K into the partner phase "
                 "(oil for WAT / GAS tracers - partitioning interwell tracer test: retardation 1 + K So / Sw gives the residual oil), "
                 "decay-chain parent (radionuclide / reactive daughters, exact Bateman for decay alone), longitudinal dispersivity and "
                 "molecular diffusion (upwind transport adds ~dx/2 of numerical dispersivity on top). Results: Tracers tab, Custom plots "
                 "(FCPR / WCPR of the tracer), 3D array TR_<name>.");
        if (trs.empty()) ImGui::TextDisabled("No TRACER keyword in the case.");
        else if (ui::beginTable("##trp", 9, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
            ui::tableColumn("Tracer"); ui::tableColumn("Phase"); ui::tableColumn("Half-life (d, 0 = stable)");
            ui::tableColumn("Adsorption kd (-)"); ui::tableColumn("Initial concentration");
            ui::tableColumn("Partition K"); ui::tableColumn("Decay parent"); ui::tableColumn("Dispersivity aL"); ui::tableColumn("Diffusion Dm (L2/d)");
            ImGui::TableHeadersRow();
            for (const auto& [nm, ph] : trs) {
                ImGui::PushID(nm.c_str());
                Case::TracerProp* P = nullptr;
                for (auto& q : case_.tracerProps) if (q.name == nm) P = &q;
                Case::TracerProp tmp; tmp.name = nm;
                Case::TracerProp& E = P ? *P : tmp;
                bool ch2 = false;
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted(nm.c_str());
                ImGui::TableNextColumn(); ImGui::TextUnformatted(ph.c_str());
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ch2 |= ImGui::InputDouble("##hl", &E.halfLife, 0, 0, "%.6g");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ch2 |= ImGui::InputDouble("##kd", &E.kd, 0, 0, "%.6g");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ch2 |= ImGui::InputDouble("##in", &E.init, 0, 0, "%.6g");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ch2 |= ImGui::InputDouble("##kp", &E.kPart, 0, 0, "%.6g");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::BeginCombo("##par", E.parent.empty() ? "(none)" : E.parent.c_str())) {
                    if (ImGui::Selectable("(none)", E.parent.empty())) { E.parent.clear(); ch2 = true; }
                    for (const auto& q : trs) if (q.first != nm && ImGui::Selectable(q.first.c_str(), E.parent == q.first)) { E.parent = q.first; ch2 = true; }
                    ImGui::EndCombo();
                }
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ch2 |= ImGui::InputDouble("##al", &E.alphaL, 0, 0, "%.6g");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ch2 |= ImGui::InputDouble("##dm", &E.dm, 0, 0, "%.6g");
                if (ch2) {
                    E.halfLife = std::max(0.0, E.halfLife); E.kd = std::max(0.0, E.kd); E.init = std::max(0.0, E.init);
                    E.kPart = std::max(0.0, E.kPart); E.alphaL = std::max(0.0, E.alphaL); E.dm = std::max(0.0, E.dm);
                    if (!P) case_.tracerProps.push_back(E);
                    markDirty();
                }
                ImGui::PopID();
            }
            ui::endTable();
        }
        ImGui::SetNextItemWidth(fs * 6);
        ui::inputString("##trnew", advTracerName_, ImVec2(0, 0), false, ImGuiInputTextFlags_CharsUppercase);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5);
        ImGui::Combo("##trph", &advTracerPhase_, "WAT\0OIL\0GAS\0");
        ImGui::SameLine();
        const std::string tn = upperTrim(advTracerName_);
        ImGui::BeginDisabled(tn.empty());
        if (ImGui::Button("+ Add tracer")) {
            static const char* kP[] = {"WAT", "OIL", "GAS"};
            const std::string rec = " '" + tn + "' '" + kP[std::clamp(advTracerPhase_, 0, 2)] + "' /\n";
            ExtraKeyword* T = nullptr;
            for (auto& k : X) if (k.name == "TRACER") T = &k;
            if (T) { std::string& b = T->body; const size_t e = b.rfind('/'); b.insert(e == std::string::npos ? b.size() : e, rec); }
            else {
                ExtraKeyword e; e.name = "TRACER"; e.section = DeckSection::Props; e.body = rec + "/\n"; e.note = "added in ResSim Studio (tracer)";
                int at = int(X.size());
                for (int i = int(X.size()) - 1; i >= 0; --i) if (int(X[size_t(i)].section) <= int(DeckSection::Props)) { at = i + 1; break; }
                if (X.empty() || int(X.front().section) > int(DeckSection::Props)) at = 0;
                X.insert(X.begin() + at, e);
            }
            // injected at concentration 1 by every injector of the matching phase from the first report step
            std::string wt;
            for (const auto& w : case_.wells)
                if ((advTracerPhase_ == 0 && w.type == WellType::WaterInjector) || (advTracerPhase_ == 2 && w.type == WellType::GasInjector))
                    wt += " '" + w.name + "' '" + tn + "' 1.0 /\n";
            if (!wt.empty()) {
                ExtraKeyword e; e.name = "WTRACER"; e.section = DeckSection::Schedule; e.step = 0; e.body = wt + "/\n"; e.note = "added in ResSim Studio (tracer)";
                X.push_back(e);
            }
            notify("Tracer " + tn + " added (TRACER" + std::string(wt.empty() ? "" : " + WTRACER at step 0") + ") - edit the records below");
            advTracerName_.clear(); markDirty();
        }
        ImGui::EndDisabled();
    }


    // ---- C29 salt precipitation / dry-out (RESSIM SALT, TOUGH3 EWASG-type)
    if (ImGui::CollapsingHeader("Salt precipitation and dry-out (ResSim engine, EWASG-type)")) {
        Case::SaltInput& S = case_.salt;
        ui::hint("Dissolved NaCl rides with the water (internal tracer NACL, kg per m3 of water). Above the solubility it precipitates "
                 "as halite; dry gas entering a water-bearing cell evaporates water (local equilibrium, vapour capacity from the "
                 "saturation pressure at the reservoir temperature) and the brine concentrates - near-well dry-out of gas producers, "
                 "CO2 / gas storage injectors. Porosity phi (1 - S_s); permeability by the Verma & Pruess (1988) tubes-in-series model "
                 "(phiR = fraction of the porosity at which k vanishes, Gamma = pore-body fraction). Results: 3D arrays SSALT (solid "
                 "fraction of the pore volume) and KSALT (permeability factor). Deck: '-- RESSIM: SALT' (OPM Flow ignores it).");
        bool c2 = ImGui::Checkbox("Salt model on", &S.on);
        ImGui::BeginDisabled(!S.on);
        if (ui::beginProps("salt")) {
            c2 |= ui::propD("Initial NaCl (kg/m3 water)", &S.c0, "%.5g");
            c2 |= ui::propD("Solubility (kg/m3 water)", &S.cSol, "%.5g", "Halite in water: ~317 kg/m3 at 25-150 C (26 wt%)");
            c2 |= ui::propD("Halite density (kg/m3)", &S.rhoSalt, "%.5g");
            c2 |= ui::propD("Verma-Pruess phiR (-)", &S.phiR, "%.3g", "permeability is zero once the porosity has fallen to phiR x phi0 (0.8 - 0.95)");
            c2 |= ui::propD("Verma-Pruess Gamma (-)", &S.gamma, "%.3g", "fractional length of the pore bodies (0.8 typical)");
            c2 |= ui::propD("Reservoir temperature (C)", &S.tempC, "%.4g", "sets the water-vapour capacity of the gas");
            c2 |= ui::propD("Injected-water NaCl (kg/m3)", &S.cInj, "%.5g");
            ui::endProps();
        }
        ImGui::EndDisabled();
        if (c2) {
            S.c0 = std::max(0.0, S.c0); S.cSol = std::max(1e-6, S.cSol); S.rhoSalt = std::max(100.0, S.rhoSalt);
            S.phiR = std::clamp(S.phiR, 0.0, 0.999); S.gamma = std::clamp(S.gamma, 0.01, 0.999); S.cInj = std::max(0.0, S.cInj);
            markDirty();
        }
    }

    // ---- C29 wellbore friction + inflow-control devices (RESSIM WELLBORE)
    if (ImGui::CollapsingHeader("Wellbore friction and inflow-control devices (ResSim engine)")) {
        ui::hint("Tubing friction between the datum and the connections (Darcy-Weisbach, Haaland friction factor, mixture density "
                 "of the toe-side inflow) and an ICD per connection: dp = K (rho_mix / rho_cal) q^2 with q the connection's "
                 "reservoir rate. The connection order is the COMPDAT order (heel to toe). Long horizontal wells: friction shifts "
                 "inflow to the heel; ICDs equalise it. Deck: '-- RESSIM: WELLBORE well ID rough [K [rho_cal]]'.");
        const bool field = case_.units == UnitSystem::Field;
        if (case_.wells.empty()) ImGui::TextDisabled("No wells.");
        else if (ui::beginTable("##wbore", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
            ui::tableColumn("Well"); ui::tableColumn("On"); ui::tableColumn(field ? "Tubing ID (ft)" : "Tubing ID (m)");
            ui::tableColumn(field ? "Roughness (ft)" : "Roughness (m)");
            ui::tableColumn(field ? "ICD K (psi/(rb/d)^2, 0 = none)" : "ICD K (bar/(rm3/d)^2, 0 = none)");
            ui::tableColumn(field ? "ICD rho_cal (lb/ft3)" : "ICD rho_cal (kg/m3)");
            ImGui::TableHeadersRow();
            for (const auto& w : case_.wells) {
                ImGui::PushID(w.name.c_str());
                int idx = -1;
                for (size_t i = 0; i < case_.wellbores.size(); ++i) if (case_.wellbores[i].well == w.name) idx = int(i);
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted(w.name.c_str());
                ImGui::TableNextColumn();
                bool on = idx >= 0;
                if (ImGui::Checkbox("##on", &on)) {
                    if (on) { Case::WellboreSpec B; B.well = w.name; if (field) { B.id = 0.25; B.rough = 5e-5; B.rhoCal = 62.4; } else { B.id = 0.076; B.rough = 1.5e-5; } case_.wellbores.push_back(B); }
                    else case_.wellbores.erase(case_.wellbores.begin() + idx);
                    markDirty();
                    idx = -1; for (size_t i = 0; i < case_.wellbores.size(); ++i) if (case_.wellbores[i].well == w.name) idx = int(i);
                }
                if (idx >= 0) {
                    Case::WellboreSpec& B = case_.wellbores[size_t(idx)];
                    bool c3 = false;
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); c3 |= ImGui::InputDouble("##id", &B.id, 0, 0, "%.5g");
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); c3 |= ImGui::InputDouble("##ro", &B.rough, 0, 0, "%.3g");
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); c3 |= ImGui::InputDouble("##k", &B.icdK, 0, 0, "%.4g");
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); c3 |= ImGui::InputDouble("##rc", &B.rhoCal, 0, 0, "%.5g");
                    if (c3) { B.id = std::max(1e-3, B.id); B.rough = std::max(0.0, B.rough); B.icdK = std::max(0.0, B.icdK); B.rhoCal = std::max(1e-6, B.rhoCal); markDirty(); }
                } else { for (int k = 0; k < 4; ++k) { ImGui::TableNextColumn(); ImGui::TextDisabled("-"); } }
                ImGui::PopID();
            }
            ui::endTable();
        }
    }

    // ---- C29 surface network (NODEPROP / BRANPROP)
    if (ImGui::CollapsingHeader("Surface network & facilities (NODEPROP / BRANPROP / GRUPNET, GSATPROD / GCONSUMP / GCONSALE)")) {
        ui::hint("Extended network: NODEPROP gives the fixed pressure of a terminal node (item 2), AS_CHOKE (item 3: automatic choke "
                 "meeting the group rate target with a common THP), ADD_GAS_LIFT_GAS (item 4); BRANPROP the pipeline of each branch "
                 "(downtree, uptree, VFPPROD table: p_in = VFP(q, p_out, ALQ), 9999 = flow line without pressure loss, 0 = branch "
                 "removed). Standard network: GRUPNET (terminal pressure / VFP per group, follows GRUPTREE moves). NEFAC scales a "
                 "node's flow. The ResSim engine balances the network implicitly every time step and gives every producer of a node "
                 "the node pressure as its THP limit. Facilities: GSATPROD satellite production (in group / field rates), GCONSUMP "
                 "gas consumption / import (taken out of the network flow and the re-injected gas), GCONSALE sales target (met by the "
                 "group gas injection) -> FGCR / FGIMR / FGSR. Results: Custom plots. The records are edited in the keyword list below.");
        int nNode = 0, nBr = 0;
        for (const auto& k : X) { if (k.name == "NODEPROP") ++nNode; if (k.name == "BRANPROP") ++nBr; }
        ImGui::Text("NODEPROP keywords: %d   BRANPROP keywords: %d   VFP tables: %d", nNode, nBr, int(case_.vfp.size()));
        static double pTerm = 0; static int vfpId = 1;
        if (pTerm <= 0) pTerm = case_.units == UnitSystem::Field ? 300.0 : 20.0;
        ImGui::SetNextItemWidth(fs * 6); ImGui::InputDouble("terminal pressure##npT", &pTerm, 0, 0, "%.4g");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4); ImGui::InputInt("branch VFP##npV", &vfpId, 0);
        ImGui::SameLine();
        ImGui::BeginDisabled(nNode > 0 || case_.wells.empty());
        if (ImGui::Button("+ Create network (FIELD + one branch per well group)")) {
            std::vector<std::string> grps;
            for (const auto& w : case_.wells) if (w.type == WellType::Producer && std::find(grps.begin(), grps.end(), w.group) == grps.end()) grps.push_back(w.group);
            char nb[96]; std::snprintf(nb, sizeof nb, " 'FIELD' %.6g /\n/\n", pTerm);
            ExtraKeyword np; np.name = "NODEPROP"; np.section = DeckSection::Schedule; np.step = 0; np.body = nb; np.note = "added in ResSim Studio (network)";
            std::string bb;
            for (const auto& g : grps) if (g != "FIELD") bb += " '" + g + "' 'FIELD' " + std::to_string(vfpId > 0 ? vfpId : 9999) + " /\n";   // 9999 = no loss
            ExtraKeyword bp; bp.name = "BRANPROP"; bp.section = DeckSection::Schedule; bp.step = 0; bp.body = bb + "/\n"; bp.note = np.note;
            X.push_back(np); X.push_back(bp);
            notify("Network added: NODEPROP FIELD + BRANPROP for " + std::to_string(grps.size()) + " group(s) - define the VFPPROD table(s)");
            markDirty();
        }
        ImGui::EndDisabled();
        if (nBr > 0 && case_.vfp.empty()) ImGui::TextColored(col4(Pal::Warn), "BRANPROP refers to VFP tables, but the case has none (VFPPROD).");
        // ---- C35 facilities: quick records (appended at the first report step; edit / move them in the list below)
        int nSat = 0, nCons = 0, nSale = 0, nStd = 0, nChk = 0;
        for (const auto& k : X) {
            nSat += k.name == "GSATPROD"; nCons += k.name == "GCONSUMP"; nSale += k.name == "GCONSALE"; nStd += k.name == "GRUPNET";
            if (k.name == "NODEPROP" && (k.body.find("YES") != std::string::npos || k.body.find("yes") != std::string::npos)) ++nChk;
        }
        ImGui::Text("GRUPNET: %d   GSATPROD: %d   GCONSUMP: %d   GCONSALE: %d   NODEPROP with YES options: %d", nStd, nSat, nCons, nSale, nChk);
        static char fGrp[32] = "SAT"; static double fA = 0, fB = 0, fC = 0, fD = 0; static char fNode[32] = "";
        ImGui::SetNextItemWidth(fs * 5); ImGui::InputText("group / node##fg", fGrp, sizeof fGrp, ImGuiInputTextFlags_CharsUppercase);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4.5f); ImGui::InputDouble("a##fa", &fA, 0, 0, "%.4g");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4.5f); ImGui::InputDouble("b##fb", &fB, 0, 0, "%.4g");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4.5f); ImGui::InputDouble("c##fc", &fC, 0, 0, "%.4g");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4.5f); ImGui::InputDouble("d##fd", &fD, 0, 0, "%.4g");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::InputText("node##fn", fNode, sizeof fNode, ImGuiInputTextFlags_CharsUppercase);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("GCONSUMP item 4: network node where the gas leaves the flow ('NO-NODES' = none; empty = default)");
        auto addFac = [&](const char* kw, const std::string& rec) {
            ExtraKeyword e; e.name = kw; e.section = DeckSection::Schedule; e.step = 0; e.body = rec + "/\n"; e.note = "added in ResSim Studio (facilities)";
            X.push_back(e); markDirty(); notify(std::string(kw) + " record added at the first report step - move / edit it in the keyword list");
        };
        auto g6 = [](double v) { char b[32]; std::snprintf(b, sizeof b, "%.6g", v); return std::string(b); };
        const std::string gq = std::string(" '") + fGrp + "' ";
        ImGui::BeginDisabled(fGrp[0] == 0);
        if (ImGui::Button("+ GSATPROD (a oil, b water, c gas, d lift gas)")) addFac("GSATPROD", gq + g6(fA) + " " + g6(fB) + " " + g6(fC) + " 1* " + g6(fD) + " /\n");
        ImGui::SameLine();
        if (ImGui::Button("+ GCONSUMP (a consumption, b import)")) addFac("GCONSUMP", gq + g6(fA) + " " + g6(fB) + (fNode[0] ? std::string(" '") + fNode + "'" : std::string()) + " /\n");
        ImGui::SameLine();
        if (ImGui::Button("+ GCONSALE (a target, b max, c min)")) addFac("GCONSALE", gq + g6(fA) + " " + g6(fB) + " " + g6(fC) + " RATE /\n");
        ImGui::SameLine();
        if (ImGui::Button("+ automatic choke node")) addFac("NODEPROP", gq + "1* YES NO /\n");
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("NODEPROP <group> 1* YES: the node pressure is raised until the group's GCONPROD rate target is met\n"
                                                      "(join the node to its parent with a BRANPROP 9999 branch)");
        ImGui::EndDisabled();
    }

    // ---- toolbar: filter + add
    ImGui::SetNextItemWidth(fs * 10);
    ui::inputString("##flt", advFilter_, ImVec2(0, 0), false, ImGuiInputTextFlags_CharsUppercase);
    if (advFilter_.empty()) { ImGui::SameLine(fs * 0.6f); ImGui::TextDisabled("filter"); }
    ImGui::SameLine(fs * 11.5f);
    ImGui::SetNextItemWidth(fs * 7);
    ui::inputString("##newkw", advNewName_, ImVec2(0, 0), false, ImGuiInputTextFlags_CharsUppercase);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7);
    ImGui::Combo("##newsec", &advNewSection_, kSecNames, 8);
    ImGui::SameLine();
    const std::string nn = upperTrim(advNewName_);
    ImGui::BeginDisabled(nn.empty());
    if (ImGui::Button("+ Add keyword")) {
        ExtraKeyword e; e.name = nn; e.section = DeckSection(advNewSection_); e.note = "added in ResSim Studio";
        const deck::KeywordSpec* sp = deck::findKeyword(nn);
        if (sp && sp->kind == deck::KwSize::Slash) e.body = " /\n/\n";
        else if (sp && sp->kind == deck::KwSize::None) e.body = "";
        else e.body = " /\n";
        if (e.section == DeckSection::Schedule) e.step = 0;
        // insert after the last keyword of the same section (keeps the deck order readable)
        int at = int(X.size());
        for (int i = int(X.size()) - 1; i >= 0; --i) if (int(X[size_t(i)].section) <= advNewSection_) { at = i + 1; break; }
        if (X.empty() || int(X.front().section) > advNewSection_) at = 0;
        X.insert(X.begin() + at, e);
        advSel_ = at; advNewName_.clear(); markDirty();
    }
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Type a keyword name (e.g. GRUPTREE, GCONPROD, VFPPROD, WELSEGS, ACTIONX, WECON, TRACER),\n"
                                                  "pick its section and add it. The keyword manual is shown below the editor.");

    // ---- list
    double day = 0; std::vector<double> stepDay(case_.schedule.size() + 1, 0.0);
    for (size_t s = 0; s < case_.schedule.size(); ++s) { stepDay[s] = day; day += case_.schedule[s].length; }
    stepDay.back() = day;
    // list height: the rows it needs, at most ~40 % of the panel (the editor goes below)
    const float listMax = std::max(fs * 10, ImGui::GetContentRegionAvail().y * 0.4f);
    const float listH = std::min(listMax, ImGui::GetTextLineHeightWithSpacing() * 1.25f * (float(X.size()) + 1.5f) + fs);
    const bool wide = ImGui::GetContentRegionAvail().x > fs * 30;   // room for the Note column
    int moveFrom = -1, moveTo = -1, del = -1;
    if (ui::beginTable("xkw", wide ? 6 : 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit,
                          ImVec2(0, std::max(listH, fs * 5)))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ui::tableColumn("#", 0, fs * 2);
        ui::tableColumn("Section", 0, fs * 5);
        ui::tableColumn("Keyword", 0, fs * 5.5f);
        ui::tableColumn("When", wide ? 0 : ImGuiTableColumnFlags_WidthStretch, fs * 6.5f);
        if (wide) ui::tableColumn("Note", ImGuiTableColumnFlags_WidthStretch);
        ui::tableColumn("", 0, fs * 4.2f);
        ImGui::TableHeadersRow();
        for (int i = 0; i < int(X.size()); ++i) {
            const ExtraKeyword& e = X[size_t(i)];
            if (!advFilter_.empty() && e.name.find(upperTrim(advFilter_)) == std::string::npos &&
                std::string(kSecNames[int(e.section)]).find(upperTrim(advFilter_)) == std::string::npos) continue;
            ImGui::PushID(i);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            char b[16]; std::snprintf(b, sizeof(b), "%d", i + 1);
            if (ImGui::Selectable(b, advSel_ == i, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) advSel_ = i;
            if (!wide && !e.note.empty() && ImGui::IsItemHovered()) ImGui::SetTooltip("%s", e.note.c_str());
            ImGui::TableNextColumn(); ImGui::TextDisabled("%s", kSecNames[int(e.section)]);
            ImGui::TableNextColumn(); ImGui::TextUnformatted(e.name.c_str());
            ImGui::TableNextColumn();
            if (e.section == DeckSection::Schedule) {
                const int s = std::clamp(e.step, 0, int(case_.schedule.size()));
                ImGui::Text("step %d, day %.0f", s, stepDay[size_t(s)]);
            }
            if (wide) { ImGui::TableNextColumn(); ImGui::TextDisabled("%s", e.note.c_str()); }
            ImGui::TableNextColumn();
            if (ImGui::SmallButton("\xE2\x86\x91") && i > 0) { moveFrom = i; moveTo = i - 1; }
            ImGui::SameLine(0, 2);
            if (ImGui::SmallButton("\xE2\x86\x93") && i + 1 < int(X.size())) { moveFrom = i; moveTo = i + 1; }
            ImGui::SameLine(0, 2);
            if (ImGui::SmallButton("\xC3\x97")) del = i;
            ImGui::PopID();
        }
        ui::endTable();
    }
    if (moveFrom >= 0) { std::swap(X[size_t(moveFrom)], X[size_t(moveTo)]); if (advSel_ == moveFrom) advSel_ = moveTo; else if (advSel_ == moveTo) advSel_ = moveFrom; markDirty(); }
    if (del >= 0) { X.erase(X.begin() + del); if (advSel_ >= int(X.size())) advSel_ = int(X.size()) - 1; markDirty(); }
    if (X.empty()) ImGui::TextDisabled("No extra keywords. Imported decks keep every keyword ResSim does not model here; or add one above.");

    // ---- editor of the selected keyword
    if (advSel_ < 0 || advSel_ >= int(X.size())) return;
    ExtraKeyword& e = X[size_t(advSel_)];
    ui::header(("Keyword " + std::to_string(advSel_ + 1) + ": " + e.name).c_str());
    bool ch = false;
    if (ui::beginProps("xkwp")) {
        ui::label("Name");
        if (ui::inputString("##nm", e.name, ImVec2(0, 0), false, ImGuiInputTextFlags_CharsUppercase)) { e.name = upperTrim(e.name); ch = true; }
        ui::label("Section");
        int sec = int(e.section);
        if (ImGui::Combo("##sec", &sec, kSecNames, 8)) { e.section = DeckSection(sec); ch = true; }
        if (e.section == DeckSection::Schedule) {
            ui::label("Report step", "The keyword is written after the well controls of this report step, before its TSTEP. "
                                     "Step 0 = start of the schedule; the last value = after the final step.");
            int st = e.step;
            if (ImGui::SliderInt("##step", &st, 0, int(case_.schedule.size()), "%d")) { e.step = st; ch = true; }
            const int s = std::clamp(e.step, 0, int(case_.schedule.size()));
            ui::label("Time"); ImGui::Text("day %.1f", stepDay[size_t(s)]);
            ui::label("Order in the step", "Before: written ahead of the step's WCONPROD / WCONINJE records (needed for VFP tables, "
                                           "group trees, UDQs they refer to). After: written behind them (e.g. records that modify them).");
            int bc = e.beforeControls ? 0 : 1;
            if (ImGui::Combo("##bc", &bc, "Before the well controls\0After the well controls\0")) { e.beforeControls = bc == 0; ch = true; }
        }
        std::string sn; const deck::Support sup = deck::ressimSupport(e.name, &sn);
        ui::label("ResSim engine");
        ImGui::TextWrapped("%s%s%s", deck::toString(sup), sn.empty() ? "" : " - ", sn.c_str());
        if (sup == deck::Support::Full)
            ImGui::TextColored(col4(Pal::Warn), "ResSim edits this keyword natively - the native data is also written (possible duplicate).");
        ui::endProps();
    }
    const float bodyH = std::max(fs * 10, ImGui::GetContentRegionAvail().y - fs * 5);
    // schema form view (parameter names, units, defaults, option lists from the OPM keyword index)
    namespace di = ressim::deckintel;
    const bool hugeBody = e.body.size() > 256 * 1024;
    if (!hugeBody) {
        if (ImGui::RadioButton("Form", advFormView_)) advFormView_ = true;
        ImGui::SameLine();
        if (ImGui::RadioButton("Text", !advFormView_)) advFormView_ = false;
        ImGui::SameLine();
    }
    bool formShown = false;
    if (advFormView_ && !hugeBody) {
        if (advFormKey_ != advSel_ || advFormSrc_ != e.body || advForm_.keyword != e.name) {
            advForm_ = di::parseForm(e.name, di::splitLines(e.body));
            advFormKey_ = advSel_;
            advFormSrc_ = e.body;
        }
        if (advForm_.table || advForm_.editable()) {
            formShown = true;
            ImGui::TextDisabled("Records of %s - hover a column for its meaning, units and default", e.name.c_str());
            ImGui::BeginChild("##advform", ImVec2(0, bodyH), ImGuiChildFlags_Borders, ImGuiWindowFlags_HorizontalScrollbar);
            if (deckFormRecords(advForm_, case_.units == UnitSystem::Metric ? "METRIC" : "FIELD")) {
                e.body = di::joinLines(di::writeForm(advForm_)) + "\n";
                advFormSrc_ = e.body;
                ch = true;
            }
            ImGui::EndChild();
        } else {
            ImGui::TextDisabled("(no record schema for a form - text view)");
        }
    }
    if (!formShown) ImGui::TextDisabled("Records (ECLIPSE syntax: items separated by spaces, 1* = default, each record ends with /)");
    constexpr size_t kBig = 256 * 1024;      // larger bodies (big VFP tables) are shown read-only unless asked
    if (formShown) {
    } else if (g_fonts.mono) ImGui::PushFont(g_fonts.mono, 0.0f);
    if (formShown) {
    } else if (e.body.size() > kBig && advBigEdit_ != advSel_) {
        size_t lines = 0; for (char ch2 : e.body) lines += ch2 == '\n';
        size_t cut = 0; for (int n = 0; n < 300 && cut != std::string::npos; ++n) cut = e.body.find('\n', cut + 1);
        std::string head = e.body.substr(0, cut == std::string::npos ? e.body.size() : cut);
        head += "\n...";
        ImGui::InputTextMultiline("##bodyro", head.data(), head.size() + 1, ImVec2(-FLT_MIN, bodyH), ImGuiInputTextFlags_ReadOnly);
        if (g_fonts.mono) ImGui::PopFont();
        ImGui::TextColored(col4(Pal::Warn), "Large keyword: %.1f MB, %zu lines - the first 300 lines are shown.", double(e.body.size()) / 1048576.0, lines);
        ImGui::SameLine();
        if (ImGui::SmallButton("Edit all")) advBigEdit_ = advSel_;
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Edit the full text (can be slow for very large tables)");
    } else {
        if (ui::inputString("##body", e.body, ImVec2(-FLT_MIN, bodyH), true, ImGuiInputTextFlags_AllowTabInput)) ch = true;
        if (g_fonts.mono) ImGui::PopFont();
    }
    if (ch) markDirty();
    // syntax check (re-parsed only when the keyword changes)
    if (advCheckSel_ != advSel_ || advCheckBody_ != e.body || advCheckName_ != e.name + char('0' + int(e.section))) {
        advCheckSel_ = advSel_; advCheckBody_ = e.body; advCheckName_ = e.name + char('0' + int(e.section));
        advCheckMsgs_ = checkKeyword(e);
    }
    if (advCheckMsgs_.empty()) ImGui::TextColored(col4(Pal::Ok), "\xE2\x9C\x93 Syntax OK (record structure as OPM Flow reads it)");
    for (auto& m : advCheckMsgs_) ImGui::TextColored(col4(m.rfind("Error", 0) == 0 ? Pal::Err : Pal::Warn), "%s", m.c_str());
    if (ImGui::CollapsingHeader("Keyword manual")) keywordDocView(e.name);
}
