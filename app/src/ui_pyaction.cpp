// C37 PYACTION editor (Model mode, Schedule -> Run-time logic):
//   records   - PYACTION name / run count (SINGLE, UNLIMITED, FIRST_TRUE) / script file / report step (Case::extraKeywords,
//               written verbatim, so OPM Flow runs the same deck)
//   script    - the .py file next to the deck: editor, templates (OPM embedded-Python style and the legacy run() form),
//               syntax check with the configured interpreter, save / reload
//   Python    - interpreter path (settings), test button
// The ResSim engine runs the scripts in a Python child process (engine/src/pyaction.cpp) at the end of every report step.
#include "app.h"
#include "widgets.h"
#include "theme.h"
#include "ressim/pyaction.h"
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>

using namespace ressim;
namespace fs = std::filesystem;

namespace {
const char* kRunCounts[] = {"SINGLE", "UNLIMITED", "FIRST_TRUE"};

struct PyTemplate { const char* label; const char* file; const char* text; };
const PyTemplate kPyTemplates[] = {
    {"Well controller: open wells one by one when the field oil rate falls (OPM example)", "wellcontroller.py",
     "import opm_embedded\n"
     "# OPM Flow embedded-Python style: the whole script runs at every call, module globals persist between calls\n"
     "summary_state = opm_embedded.current_summary_state\n"
     "schedule = opm_embedded.current_schedule\n"
     "\n"
     "if 'setup_done' not in locals():\n"
     "    setup_done = True\n"
     "    wells_to_open = ['PROD2']      # opened in this order\n"
     "    target_fopr = 8000.0           # field oil rate below which the next well opens\n"
     "    min_wait_days = 50.0           # minimum time between two openings\n"
     "    last_open = -1e30\n"
     "\n"
     "days = summary_state.elapsed() / 86400.0\n"
     "if wells_to_open and summary_state['FOPR'] < target_fopr and days - last_open >= min_wait_days:\n"
     "    well = wells_to_open.pop(0)\n"
     "    schedule.open_well(well)\n"
     "    last_open = days\n"
     "    opm_embedded.OpmLog.info('opened %s at day %.1f (FOPR %.1f)' % (well, days, summary_state['FOPR']))\n"},
    {"Shut producers above a water cut", "wct_shut.py",
     "import opm_embedded\n"
     "summary_state = opm_embedded.current_summary_state\n"
     "schedule = opm_embedded.current_schedule\n"
     "\n"
     "if 'shut' not in locals():\n"
     "    shut = set()\n"
     "    max_wct = 0.9\n"
     "\n"
     "for well in summary_state.wells():\n"
     "    if well in shut or not summary_state.has_well_var(well, 'WWCT'):\n"
     "        continue\n"
     "    if summary_state.well_var(well, 'WWCT') > max_wct and summary_state.well_var(well, 'WOPR') > 0:\n"
     "        schedule.shut_well(well)\n"
     "        shut.add(well)\n"
     "        print('shut', well, 'water cut', summary_state.well_var(well, 'WWCT'))\n"},
    {"Cut the producer rate once the cumulative oil exceeds a value (insert_keywords)", "rate_cut.py",
     "import opm_embedded\n"
     "summary_state = opm_embedded.current_summary_state\n"
     "schedule = opm_embedded.current_schedule\n"
     "\n"
     "if 'done' not in locals():\n"
     "    done = False\n"
     "\n"
     "if not done and summary_state['FOPT'] > 1.0e7:\n"
     "    done = True\n"
     "    # any SCHEDULE keywords; they act from the next report step (like an ACTIONX body)\n"
     "    schedule.insert_keywords('''\n"
     "WELTARG\n"
     " 'PROD' ORAT 5000 /\n"
     "/\n"
     "''')\n"
     "    opm_embedded.OpmLog.info('rate cut at day %.1f' % (summary_state.elapsed() / 86400.0))\n"},
    {"Legacy form: run() with the ACTIONX callback", "legacy_action.py",
     "# legacy PYACTION form (OPM 2020-2023): the module body runs once, run() at every call.\n"
     "# Return True / False (FIRST_TRUE stops after the first True). actionx_callback(name, wells) applies the\n"
     "# keywords of the ACTIONX block 'name' of the deck ('?' in its body = the wells given here).\n"
     "state = {'calls': 0}\n"
     "\n"
     "def run(ecl_state, schedule, report_step, summary_state, actionx_callback):\n"
     "    state['calls'] += 1\n"
     "    wet = [w for w in summary_state.wells() if summary_state.has_well_var(w, 'WWCT') and summary_state.well_var(w, 'WWCT') > 0.8]\n"
     "    if wet:\n"
     "        actionx_callback('ACT1', wet)\n"
     "        return True\n"
     "    return False\n"},
};

std::string readText(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return std::string();
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}
std::vector<std::string> recTokens(const std::string& body) {
    std::vector<pyact::Action> a;
    std::string err;
    a = pyact::parseActions({{0, body}}, err);
    if (a.empty()) return {};
    return {a[0].name, a[0].runCount, a[0].file};
}
std::string writeRec(const std::string& name, const std::string& rc, const std::string& file) {
    return " '" + name + "' '" + rc + "' /\n '" + file + "' /\n";
}
}  // namespace

void App::editPyActions() {
    auto& X = case_.extraKeywords;
    const float fs = ImGui::GetFontSize();
    const int nSteps = std::max(1, int(case_.schedule.size()));
    static int sel = -1;
    static std::string script, loadedFor, checkMsg;
    static bool scriptDirty = false, checkOk = false;
    static std::shared_ptr<pyact::Host> host;
    static std::string testMsg;
    ui::header("PYACTION - Python actions");
    ui::hint("Python scripts called at the end of every report step (OPM Flow keyword PYACTION; ResSim runs them too, in a Python 3 "
             "child process). A script reads summary values (opm_embedded.current_summary_state: ss['FOPR'], ss.well_var('P1', "
             "'WWCT'), ss.elapsed()) and changes the schedule (current_schedule.open_well / shut_well / insert_keywords), which "
             "acts from the next report step like an ACTIONX body. Globals persist between calls. OPM Flow calls the scripts every "
             "time step and needs a build with embedded Python.");
    const fs::path base = !case_.deckDir.empty() ? fs::u8path(case_.deckDir) : (!path_.empty() ? fs::u8path(path_).parent_path() : fs::path());
    auto resolve = [&](const std::string& f) {
        fs::path p = fs::u8path(f);
        if (p.is_relative() && !base.empty()) p = base / p;
        return p;
    };
    // ---- interpreter
    if (ui::beginProps("##pyexe", 0.32f)) {
        ui::label("Python interpreter", "python.exe / python3 path. Empty: RESSIM_PYTHON, then python / python3 / py on the PATH.");
        ImGui::SetNextItemWidth(-fs * 9);
        bool ch = ui::inputString("##pyx", &settings_.pythonExe);
        ImGui::SameLine();
        if (ImGui::SmallButton("Browse##pyx")) {
            std::string p = settings_.pythonExe;
            if (fileDialog(false, "Python (python*.exe)\0python*.exe;py.exe\0All files (*.*)\0*.*\0", "exe", p)) { settings_.pythonExe = p; ch = true; }
        }
        ImGui::SameLine();
        if (ImGui::SmallButton("Test##pyx")) {
            pyact::Host h; std::string err;
            testMsg = h.start(settings_.pythonExe, err) ? "OK: Python " + h.version() : "Error: " + err;
            h.stop();
        }
        if (ch) { saveSettings(); host.reset(); }
        ui::endProps();
    }
    if (!testMsg.empty()) ImGui::TextColored(col4(testMsg.rfind("OK", 0) == 0 ? Pal::Ok : Pal::Warn), "%s", testMsg.c_str());
    // ---- records
    std::vector<int> items;
    for (int i = 0; i < int(X.size()); ++i) if (X[size_t(i)].name == "PYACTION") items.push_back(i);
    if (sel >= int(X.size()) || (sel >= 0 && X[size_t(sel)].name != "PYACTION")) sel = -1;
    if (ui::beginTable("##pyacts", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY,
                       ImVec2(-1, fs * std::clamp(1.6f + 1.35f * float(items.size()), 3.0f, 8.0f)))) {
        ui::tableColumn("Name", ImGuiTableColumnFlags_WidthFixed, fs * 9);
        ui::tableColumn("Run count", ImGuiTableColumnFlags_WidthFixed, fs * 6);
        ui::tableColumn("Script");
        ui::tableColumn("From step", ImGuiTableColumnFlags_WidthFixed, fs * 5);
        ui::tableColumn("Status", ImGuiTableColumnFlags_WidthFixed, fs * 9);
        ImGui::TableHeadersRow();
        for (int i : items) {
            const ExtraKeyword& e = X[size_t(i)];
            const auto t = recTokens(e.body);
            ImGui::TableNextRow(); ImGui::PushID(i);
            ImGui::TableNextColumn();
            if (ImGui::Selectable(t.empty() ? "(invalid)" : t[0].c_str(), sel == i, ImGuiSelectableFlags_SpanAllColumns)) { sel = i; loadedFor.clear(); }
            ImGui::TableNextColumn(); ImGui::TextUnformatted(t.empty() ? "" : t[1].c_str());
            ImGui::TableNextColumn(); ImGui::TextUnformatted(t.empty() ? "" : t[2].c_str());
            ImGui::TableNextColumn(); ImGui::Text("%d", e.step + 1);
            ImGui::TableNextColumn();
            std::error_code ec;
            if (t.empty()) ImGui::TextColored(col4(Pal::Warn), "bad records");
            else if (base.empty() && fs::u8path(t[2]).is_relative()) ImGui::TextColored(col4(Pal::Warn), "save the deck first");
            else if (fs::exists(resolve(t[2]), ec)) ImGui::TextColored(col4(Pal::Ok), "script found");
            else ImGui::TextColored(col4(Pal::Warn), "script missing");
            ImGui::PopID();
        }
        ui::endTable();
    }
    static int newStep = 1, tplNew = 0;
    newStep = std::clamp(newStep, 1, nSteps);
    ImGui::SetNextItemWidth(fs * 6); ImGui::InputInt("at report step##py", &newStep); newStep = std::clamp(newStep, 1, nSteps);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 16);
    {
        std::vector<std::string> labels;
        for (const auto& t : kPyTemplates) labels.push_back(t.label);
        if (ImGui::BeginCombo("##pytplnew", labels[size_t(tplNew)].c_str(), ImGuiComboFlags_HeightLarge)) {
            for (int k = 0; k < int(labels.size()); ++k) if (ImGui::Selectable(labels[size_t(k)].c_str(), k == tplNew)) tplNew = k;
            ImGui::EndCombo();
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("New PYACTION")) {
        const PyTemplate& T = kPyTemplates[tplNew];
        std::string name = "PYACT" + std::to_string(items.size() + 1);
        ExtraKeyword e; e.name = "PYACTION"; e.section = DeckSection::Schedule; e.step = newStep - 1;
        e.body = writeRec(name, tplNew == 3 ? "FIRST_TRUE" : "UNLIMITED", T.file);
        e.note = "Python action (C37 editor)";
        X.push_back(e); sel = int(X.size()) - 1;
        const fs::path p = resolve(T.file);
        std::error_code ec;
        if (!base.empty() && !fs::exists(p, ec)) { std::ofstream f(p, std::ios::binary); f << T.text; }
        loadedFor.clear();
        if (base.empty()) { script = T.text; loadedFor = std::to_string(sel) + "\x1f" + T.file; scriptDirty = true; }
        markDirty();
    }
    if (sel >= 0) {
        ImGui::SameLine();
        if (ui::dangerButton("Delete##py")) { X.erase(X.begin() + sel); sel = -1; loadedFor.clear(); markDirty(); return; }
    }
    if (sel < 0) { if (!items.empty()) ImGui::TextDisabled("Select a PYACTION to edit it and its script."); return; }
    // ---- selected record
    ExtraKeyword& e = X[size_t(sel)];
    auto t = recTokens(e.body);
    if (t.empty()) t = {"PYACT", "SINGLE", "script.py"};
    bool ch = false;
    int rc = 0; for (int k = 0; k < 3; ++k) if (t[1] == kRunCounts[k]) rc = k;
    int st = e.step + 1;
    if (ui::beginProps("##pyrec", 0.32f)) {
        ch |= ui::propString("Name", &t[0], "action name (unique; a later PYACTION with the same name replaces it)");
        if (ui::propCombo("Run count", &rc, kRunCounts, 3, "SINGLE: called once. UNLIMITED: every report step. FIRST_TRUE: until the "
                                                           "legacy run() returns True.")) { t[1] = kRunCounts[rc]; ch = true; }
        ch |= ui::propString("Script file", &t[2], "relative to the deck folder (as OPM Flow) or absolute");
        if (ui::propI("From report step", &st, 1, nSteps)) { e.step = st - 1; markDirty(); }
        ui::endProps();
    }
    if (ImGui::SmallButton("Choose script...")) {
        std::string p = resolve(t[2]).u8string();
        if (fileDialog(false, "Python script (*.py)\0*.py\0All files (*.*)\0*.*\0", "py", p)) {
            std::error_code ec;
            const fs::path rel = base.empty() ? fs::u8path(p) : fs::relative(fs::u8path(p), base, ec);
            t[2] = (ec || rel.empty() ? fs::u8path(p) : rel).generic_u8string();
            ch = true;
        }
    }
    if (ch) { e.body = writeRec(t[0], t[1], t[2]); markDirty(); }
    // ---- script editor
    const fs::path sp = resolve(t[2]);
    const std::string key = std::to_string(sel) + "\x1f" + t[2];
    if (loadedFor != key) { script = readText(sp); loadedFor = key; scriptDirty = false; checkMsg.clear(); }
    ImGui::TextDisabled("Script: %s%s", sp.u8string().c_str(), scriptDirty ? "  (modified)" : "");
    if (ui::inputString("##pyscript", script, ImVec2(-1, fs * 16), true, ImGuiInputTextFlags_AllowTabInput)) { scriptDirty = true; checkMsg.clear(); }
    ImGui::BeginDisabled(base.empty() && sp.is_relative());
    if (ImGui::Button("Save script")) {
        std::error_code ec;
        fs::create_directories(sp.parent_path(), ec);
        std::ofstream f(sp, std::ios::binary);
        if (f) { f << script; scriptDirty = false; notify("Saved " + sp.u8string()); }
        else notify("Cannot write " + sp.u8string());
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Reload")) { loadedFor.clear(); }
    ImGui::SameLine();
    if (ImGui::Button("Check syntax")) {
        if (!host || !host->running()) {
            host = std::make_shared<pyact::Host>();
            std::string err;
            if (!host->start(settings_.pythonExe, err)) { checkOk = false; checkMsg = err; host.reset(); }
        }
        if (host) {
            pyact::Json req = pyact::Json::obj(), resp;
            req.set("cmd", pyact::Json::str("check")).set("source", pyact::Json::str(script)).set("file", pyact::Json::str(t[2]));
            std::string err;
            if (!host->call(req, nullptr, resp, err, 20)) { checkOk = false; checkMsg = err; host.reset(); }
            else {
                const pyact::Json* ok = resp.get("ok");
                checkOk = ok && ok->b;
                const pyact::Json* em = resp.get("error");
                checkMsg = checkOk ? "Syntax OK (Python " + host->version() + ")" : (em ? em->s : std::string("error"));
            }
        }
    }
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 14);
    if (ImGui::BeginCombo("##pytplrep", "Replace with a template...", ImGuiComboFlags_HeightLarge)) {
        for (const auto& T : kPyTemplates) if (ImGui::Selectable(T.label)) { script = T.text; scriptDirty = true; checkMsg.clear(); }
        ImGui::EndCombo();
    }
    if (!checkMsg.empty()) ImGui::TextColored(col4(checkOk ? Pal::Ok : Pal::Warn), "%s", checkMsg.c_str());
    ui::help("opm_embedded API in ResSim: current_summary_state (ss['FOPR'], 'FOPR' in ss, ss.get(k, default), ss.well_var(w, v), "
             "ss.group_var(g, v), has_well_var / has_group_var, ss.wells(pattern), ss.groups(), ss.elapsed() in seconds, "
             "ss.update(k, v)), current_schedule (open_well / shut_well / stop_well(name), insert_keywords(text), well_names(pattern), "
             "get_well(name).status(), .start), current_report_step, OpmLog.info / warning / error. Any F* / W* / G* vector the UDQ "
             "engine knows (and UDQ values) can be read. print() goes to the run log. The Python state is not stored in restart files.");
}
