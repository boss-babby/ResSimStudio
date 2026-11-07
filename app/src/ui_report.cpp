// C33 Report generator (Word / PowerPoint / HTML / Markdown from the current run) and deck regression runner
// (manifest of decks + reference summaries, e.g. OPM/opm-tests fetched by tools/fetch_opm_tests.py).
#include "app.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "ressim/api.h"
#include "ressim/report.h"
#include "ressim/regress.h"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include "platform.h"

namespace fs = std::filesystem;

namespace {
inline ImVec4 C4(ImU32 c) { return ImGui::ColorConvertU32ToFloat4(c); }
void openPath(const std::string& p) { plat::openPath(p); }
const char* kFormats[] = {"Word (.docx)", "PowerPoint (.pptx)", "HTML (self-contained)", "Markdown (.md)"};
const char* kExt[] = {"docx", "pptx", "html", "md"};
const char* kFilter[] = {"Word document (*.docx)\0*.docx\0All files (*.*)\0*.*\0", "PowerPoint (*.pptx)\0*.pptx\0All files (*.*)\0*.*\0",
                         "HTML (*.html)\0*.html\0All files (*.*)\0*.*\0", "Markdown (*.md)\0*.md\0All files (*.*)\0*.*\0"};
}  // namespace

// ---------------------------------------------------------------- report generator
void App::reportWindow() {
    if (!showReport_) return;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 40, ImGui::GetFontSize() * 31), ImGuiCond_Appearing);
    if (!ImGui::Begin(ICON_FILE_TEXT "  Report generator###reportwin", &showReport_, ImGuiWindowFlags_NoDocking)) { ImGui::End(); return; }
    bool has = false; std::string engine;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        has = R.valid && R.summary.size() >= 2 && R.engine.rfind("ResSim", 0) == 0;
        engine = R.engine;
    }
    ImGui::TextWrapped("Writes a report of the current ResSim run: model description, in-place volumes, field and well "
                       "charts (native, editable charts in Word / PowerPoint), pressure / saturation maps, engine messages and "
                       "an optional comparison with a reference run (OPM Flow, ECLIPSE or history).");
    ImGui::Separator();
    if (!has) {
        ImGui::TextColored(C4(Pal::Warn), "No ResSim results yet%s - run the model first (F5).", engine.empty() || engine.rfind("ResSim", 0) == 0 ? "" : " (the last run was external)");
        ImGui::End(); return;
    }
    auto& S = rep_;
    if (ui::beginProps("##repprops")) {
        ui::label("Format");
        ImGui::SetNextItemWidth(-1); ImGui::Combo("##fmt", &S.format, kFormats, 4);
        ui::label("Title");
        ImGui::SetNextItemWidth(-1); ImGui::InputTextWithHint("##title", "(case title)", S.title, sizeof S.title);
        ui::label("Author");
        ImGui::SetNextItemWidth(-1); ImGui::InputText("##author", S.author, sizeof S.author);
        ui::label("Model tables");
        ImGui::Checkbox("##inputs", &S.inputs);
        ui::label("Well charts");
        ImGui::Checkbox("##wells", &S.wellCharts); ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
        ImGui::InputInt("max wells##mw", &S.maxWells); S.maxWells = std::clamp(S.maxWells, 1, 64);
        ui::label("Maps");
        ImGui::Checkbox("##maps", &S.maps); ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
        ImGui::InputInt("layer (0 = auto)##ly", &S.layer); S.layer = std::max(0, S.layer);
        ui::label("Engine messages");
        ImGui::Checkbox("##log", &S.log);
        ui::label("Reference run");
        ImGui::SetNextItemWidth(-ImGui::GetFontSize() * 5); ImGui::InputTextWithHint("##ref", "optional .SMSPEC (OPM Flow / ECLIPSE)", S.ref, sizeof S.ref);
        ImGui::SameLine(); if (ImGui::SmallButton("Browse##ref")) { std::string p; if (fileDialog(false, "Summary (*.SMSPEC)\0*.SMSPEC\0All files (*.*)\0*.*\0", "SMSPEC", p)) std::snprintf(S.ref, sizeof S.ref, "%s", p.c_str()); }
        ui::label("Reference label");
        ImGui::SetNextItemWidth(-1); ImGui::InputText("##rlabel", S.refLabel, sizeof S.refLabel);
        ui::endProps();
    }
    ImGui::TextUnformatted("Notes (paragraph after the summary)");
    ImGui::InputTextMultiline("##notes", S.notes, sizeof S.notes, ImVec2(-1, ImGui::GetFontSize() * 4));
    ImGui::Checkbox("Open the report when written", &S.openAfter);
    ImGui::Separator();
    if (ImGui::Button(ICON_DOWNLOAD "  Generate report...")) {
        std::string p = (fs::u8path(runDirFor(path_, case_.title, "REPORTS")) / (caseBaseName(case_) + "_report." + kExt[S.format])).u8string();
        { std::error_code ec; fs::create_directories(fs::u8path(p).parent_path(), ec); }
        if (fileDialog(true, kFilter[S.format], kExt[S.format], p)) {
            ressim::api::RunOutput run;
            ressim::Case input;
            {
                std::lock_guard<std::mutex> l(runner_.mutex());
                const RunResults& R = runner_.results();
                input = R.input;
                run.ok = true; run.engine = R.engine; run.summary = R.summary; run.snaps = R.snaps;
                if (R.grid && R.engine == "ResSim") { run.globalOf = R.grid->globalOf; run.depth = R.grid->depth; run.porv = R.grid->poreVolRef; }
                else run.snaps.clear();            // maps only for black-oil runs (per-active-cell arrays)
                run.stats = R.stats; run.seconds = R.totalTime; run.wellNames = R.wellNames; run.initial = R.initial;
                const auto& L = runner_.log();
                for (size_t i = L.size() > 4000 ? L.size() - 4000 : 0; i < L.size(); ++i) run.log.push_back(L[i]);
            }
            ressim::report::Options o;
            o.title = S.title; o.author = S.author; o.notes = S.notes; o.inputs = S.inputs; o.wellCharts = S.wellCharts;
            o.maxWells = S.maxWells; o.maps = S.maps; o.mapLayer = S.layer - 1; o.log = S.log;
            const ressim::ecl::Summary sum = ressim::ecl::summaryFromRun(run.summary, run.wellNames, input.units);
            ressim::ecl::Summary refS; const ressim::ecl::Summary* rp = nullptr;
            std::string err;
            if (S.ref[0]) { if (ressim::ecl::readSummary(S.ref, refS, &err)) rp = &refS; else notify("Reference not read: " + err); }
            err.clear();
            const ressim::report::Doc d = ressim::report::simulationReport(input, run, sum, o, rp, S.refLabel[0] ? S.refLabel : "reference");
            if (!ressim::report::writeDoc(p, d, err)) notify("Report not written: " + err);
            else {
                std::error_code ec; const auto sz = fs::file_size(fs::u8path(p), ec);
                char b[96]; std::snprintf(b, sizeof b, " (%.0f KB)", ec ? 0.0 : double(sz) / 1024.0);
                notify("Wrote " + fs::u8path(p).filename().u8string() + b);
                runner_.addLog("=== Report written: " + p);
                S.last = p;
                if (S.openAfter) openPath(p);
            }
        }
    }
    if (!S.last.empty()) { ImGui::SameLine(); if (ImGui::Button(ICON_FOLDER "  Open last")) openPath(S.last); ImGui::SameLine(); ImGui::TextDisabled("%s", fs::u8path(S.last).filename().u8string().c_str()); }
    ImGui::End();
}

// ---------------------------------------------------------------- deck regression
struct App::RegressJob {
    std::mutex m;
    std::vector<std::string> log;
    struct Row { std::string name, group, status, key, msg; double worst = 0, sec = 0, prev = -1; bool slow = false; };
    std::vector<Row> rows;
    std::atomic<bool> running{false}, cancel{false};
    std::string outDir, error;
    int bad = 0;
};

void App::regressionWindow() {
    if (!showRegress_) return;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 56, ImGui::GetFontSize() * 36), ImGuiCond_Appearing);
    if (!ImGui::Begin(ICON_CHECK "  Deck regression (OPM / ECLIPSE references)###regresswin", &showRegress_, ImGuiWindowFlags_NoDocking)) { ImGui::End(); return; }
    auto& S = reg_;
    ImGui::TextWrapped("Runs every deck of a manifest with the ResSim engine and compares the summary vectors with reference "
                       "results (relative RMS per vector). tools/fetch_opm_tests.py downloads OPM/opm-tests with the OPM Flow "
                       "reference runs and writes the manifest; the nightly scripts run the same harness. Writes report.html, "
                       "report.docx, junit.xml and a run-time history (performance regressions are flagged).");
    ImGui::Separator();
    const bool busy = regJob_ && regJob_->running;
    if (ui::beginProps("##regprops")) {
        ui::label("Manifest");
        ImGui::SetNextItemWidth(-ImGui::GetFontSize() * 5); ImGui::InputText("##man", S.manifest, sizeof S.manifest);
        ImGui::SameLine(); if (ImGui::SmallButton("Browse##man")) { std::string p; if (fileDialog(false, "Manifest (*.txt)\0*.txt\0All files (*.*)\0*.*\0", "txt", p)) std::snprintf(S.manifest, sizeof S.manifest, "%s", p.c_str()); }
        ui::label("Output folder");
        ImGui::SetNextItemWidth(-1); ImGui::InputTextWithHint("##out", "(next to the manifest: regression_out)", S.outDir, sizeof S.outDir);
        ui::label("Filter (name contains)");
        ImGui::SetNextItemWidth(-1); ImGui::InputText("##filt", S.filter, sizeof S.filter);
        ui::label("Tolerance override (%)");
        ImGui::SetNextItemWidth(-1); ImGui::InputDouble("##tol", &S.tolPct, 0, 0, "%.2f"); ImGui::SetItemTooltip("0 = the manifest tolerance");
        ui::label("Max seconds per case");
        ImGui::SetNextItemWidth(-1); ImGui::InputDouble("##maxs", &S.maxSec, 0, 0, "%.0f");
        ui::endProps();
    }
    ImGui::BeginDisabled(busy || !S.manifest[0]);
    if (ImGui::Button(ICON_PLAY "  Run regression")) {
        auto job = std::make_shared<RegressJob>();
        ressim::regress::Manifest m; std::string err;
        if (!ressim::regress::readManifest(S.manifest, m, err)) notify("Manifest: " + err);
        else {
            ressim::regress::Options o;
            o.outDir = S.outDir[0] ? std::string(S.outDir) : (fs::u8path(S.manifest).parent_path() / "regression_out").u8string();
            o.filter = S.filter; o.tolOverride = S.tolPct > 0 ? S.tolPct / 100.0 : -1; o.maxSeconds = S.maxSec;
            job->outDir = o.outDir; job->running = true;
            std::weak_ptr<RegressJob> wj = job;
            o.log = [wj](const std::string& s) { if (auto j = wj.lock()) { std::lock_guard<std::mutex> l(j->m); j->log.push_back(s); } };
            o.cancel = &job->cancel;   // the thread holds the job (shared_ptr) for its whole lifetime
            regJob_ = job;
            std::thread([job, m, o]() {
                try {
                    const auto R = ressim::regress::runAll(m, o);
                    std::string e;
                    const int bad = ressim::regress::writeOutputs(R, o, e);
                    std::lock_guard<std::mutex> l(job->m);
                    for (const auto& r : R) job->rows.push_back({r.name, r.group, r.status, r.worstKey, r.message, r.worst, r.seconds, r.prevSeconds, r.slow});
                    job->bad = bad; job->error = e;
                } catch (const std::exception& e) { std::lock_guard<std::mutex> l(job->m); job->error = e.what(); }
                job->running = false;
            }).detach();
        }
    }
    ImGui::EndDisabled();
    if (busy) {
        ImGui::SameLine();
        ImGui::BeginDisabled(regJob_->cancel.load());
        if (ImGui::Button(ICON_X "  Cancel")) { regJob_->cancel = true; notify("Deck regression: cancelling - the running case is stopped, the report covers the finished cases"); }
        ImGui::EndDisabled();
    }
    if (regJob_) {
        std::lock_guard<std::mutex> l(regJob_->m);
        const auto& J = *regJob_;
        ImGui::SameLine();
        if (J.running) ImGui::TextColored(C4(Pal::Warn), J.cancel.load() ? "cancelling ... %d log lines" : "running ... %d log lines", int(J.log.size()));
        else {
            int np = 0, nf = 0, ne = 0, ns = 0;
            for (const auto& r : J.rows) { np += r.status == "PASS"; nf += r.status == "FAIL"; ne += r.status == "ERROR"; ns += r.status == "SKIP"; }
            ImGui::TextColored(nf + ne ? C4(Pal::Warn) : C4(Pal::Ok), "%d PASS  %d FAIL  %d ERROR  %d SKIP", np, nf, ne, ns);
            ImGui::SameLine(); if (ImGui::SmallButton("Open report.html")) openPath(J.outDir + "/report.html");
            ImGui::SameLine(); if (ImGui::SmallButton("Open folder")) openPath(J.outDir);
            if (!J.error.empty()) ImGui::TextColored(C4(Pal::Err), "%s", J.error.c_str());
        }
        if (ImGui::BeginTabBar("##regtabs")) {
            if (ui::tabItem("Results")) {
                if (ImGui::BeginTable("##regt", 6, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable | ImGuiTableFlags_Sortable, ImVec2(0, -1))) {
                    ImGui::TableSetupScrollFreeze(0, 1);
                    ImGui::TableSetupColumn("Case"); ImGui::TableSetupColumn("Group"); ImGui::TableSetupColumn("Status");
                    ImGui::TableSetupColumn("Worst vector"); ImGui::TableSetupColumn("relRMS %"); ImGui::TableSetupColumn("Time s");
                    ImGui::TableHeadersRow();
                    std::vector<int> idx(J.rows.size()); for (size_t i = 0; i < idx.size(); ++i) idx[i] = int(i);
                    if (ImGuiTableSortSpecs* ss = ImGui::TableGetSortSpecs(); ss && ss->SpecsCount > 0) {
                        const int c = ss->Specs[0].ColumnIndex; const bool asc = ss->Specs[0].SortDirection == ImGuiSortDirection_Ascending;
                        std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) {
                            const auto& A = J.rows[size_t(a)]; const auto& B = J.rows[size_t(b)];
                            bool lt = c == 0 ? A.name < B.name : c == 1 ? A.group < B.group : c == 2 ? A.status < B.status : c == 3 ? A.key < B.key : c == 4 ? A.worst < B.worst : A.sec < B.sec;
                            bool gt = c == 0 ? B.name < A.name : c == 1 ? B.group < A.group : c == 2 ? B.status < A.status : c == 3 ? B.key < A.key : c == 4 ? B.worst < A.worst : B.sec < A.sec;
                            return asc ? lt : gt; });
                    }
                    for (int i : idx) {
                        const auto& r = J.rows[size_t(i)];
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(r.name.c_str()); if (!r.msg.empty()) ImGui::SetItemTooltip("%s", r.msg.c_str());
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(r.group.c_str());
                        ImGui::TableNextColumn(); ImGui::TextColored(r.status == "PASS" ? C4(Pal::Ok) : r.status == "SKIP" ? ImVec4(0.6f, 0.6f, 0.6f, 1) : r.status == "FAIL" ? C4(Pal::Warn) : C4(Pal::Err), "%s", r.status.c_str());
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(r.key.c_str());
                        ImGui::TableNextColumn(); if (r.status == "PASS" || r.status == "FAIL") ImGui::Text("%.2f", 100 * r.worst);
                        ImGui::TableNextColumn(); ImGui::Text("%.1f%s", r.sec, r.slow ? "  (slower)" : "");
                    }
                    ImGui::EndTable();
                }
                ImGui::EndTabItem();
            }
            if (ui::tabItem("Log")) {
                ImGui::BeginChild("##reglog", ImVec2(0, -1), ImGuiChildFlags_Borders);
                for (const auto& s : J.log) ImGui::TextUnformatted(s.c_str());
                if (J.running && ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 4) ImGui::SetScrollHereY(1.0f);
                ImGui::EndChild();
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
    }
    ImGui::End();
}
