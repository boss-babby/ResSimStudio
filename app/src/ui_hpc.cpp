// C16: high-performance computing settings and back ends
//  * MPI cluster batches: Studies (sensitivity, Monte-Carlo, history-matching populations, ensembles) hand their runs to
//    "mpiexec -n N ResSimCLI.exe --mpi-batch <jobdir>". Rank 0 dispatches cases dynamically, every other rank simulates;
//    with a shared (UNC) job folder and mpiexec -hosts the runs spread over several machines. Results come back as
//    .rsres files and enter the study exactly like local runs. If mpiexec cannot start, the local thread pool is used.
//  * GPU linear solver: the sokol_gfx compute back end (gpu_compute.cpp) is registered with the engine; cases with
//    Numerics -> "GPU linear solver" solve their block systems on the device (float32) inside double refinement.
#include <sstream>
#include "platform.h"
#include "app.h"
#include "widgets.h"
#include "theme.h"
#include "opm_engine.h"
#include "gpu_compute.h"
#include "ressim/io.h"
#include "ressim/study.h"
#include "ressim/results_file.h"
#include "ressim/bsr.h"
#include "ressim/linsolve.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <thread>
#include <mutex>

namespace fs = std::filesystem;

namespace {

std::string exeFolder() { return fs::u8path(plat::exePath()).parent_path().u8string(); }

std::string findMpiexec(const std::string& configured) {
    std::error_code ec;
    if (!configured.empty()) {
        fs::path p = fs::u8path(configured);
        if (fs::is_directory(p, ec)) p /= std::string("mpiexec") + plat::exeSuffix();
        if (fs::exists(p, ec)) return p.u8string();
    }
#ifdef _WIN32
    wchar_t buf[MAX_PATH] = {0};
    if (GetEnvironmentVariableW(L"MSMPI_BIN", buf, MAX_PATH)) {
        fs::path p = fs::path(buf) / "mpiexec.exe";
        if (fs::exists(p, ec)) return p.u8string();
    }
    for (const char* c : {"C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe", "C:\\Program Files (x86)\\Microsoft MPI\\Bin\\mpiexec.exe"})
        if (fs::exists(fs::u8path(c), ec)) return c;
    if (SearchPathW(nullptr, L"mpiexec.exe", nullptr, MAX_PATH, buf, nullptr)) return fs::path(buf).u8string();
#else
    // POSIX: mpiexec on PATH (MPICH / Open MPI / Intel MPI)
    const std::string path = plat::getEnv("PATH");
    size_t a = 0;
    while (a <= path.size()) {
        size_t b = path.find(':', a); if (b == std::string::npos) b = path.size();
        const fs::path p = fs::u8path(path.substr(a, b - a)) / "mpiexec";
        if (b > a && fs::exists(p, ec)) return p.u8string();
        a = b + 1;
    }
#endif
    return std::string();
}

// C17: which MPI family an mpiexec belongs to, and the matching runtime DLL for ResSimCLI. Every rank must load the library
// of the launcher that started it (MS-MPI msmpi.dll <-> MS-MPI mpiexec, Intel MPI impi.dll <-> Intel hydra mpiexec), otherwise
// each rank initialises a separate one-process world. Both share the MPICH ABI that ResSimCLI binds at run time.
struct MpiFamily { std::string name, dll; };
MpiFamily mpiFamily(const std::string& mpiexec) {
    MpiFamily f;
    if (mpiexec.empty()) return f;
    std::error_code ec;
    const fs::path b = fs::u8path(mpiexec).parent_path();
    std::string low = mpiexec; for (char& c : low) c = char(std::tolower((unsigned char)c));
    for (const fs::path& d : {b / "impi.dll", b / "release" / "impi.dll", b.parent_path() / "bin" / "release" / "impi.dll",
                              b.parent_path() / "redist" / "intel64" / "impi.dll", b.parent_path() / "lib" / "release" / "impi.dll"})
        if (fs::exists(d, ec)) { f.name = "Intel MPI"; f.dll = d.u8string(); return f; }
    if (low.find("intel") != std::string::npos) { f.name = "Intel MPI"; f.dll = "impi.dll"; return f; }
    f.name = "Microsoft MPI"; f.dll = "msmpi.dll";
    return f;
}
// set before every mpiexec start (the child environment is inherited from this process)
void prepareMpiEnv(const std::string& mpiexec) {
    const MpiFamily f = mpiFamily(mpiexec);
    if (!f.dll.empty()) plat::setEnv("RESSIM_MPI_LIB", f.dll);
    // C21: rank 0 runs the Newton assembly with OpenMP threads while the other ranks sleep; its threads must not spin
    // during the distributed solve (active OpenMP waiting there made a 2-rank run 3.5x slower in the H3 tests)
    if (plat::getEnv("OMP_WAIT_POLICY").empty()) plat::setEnv("OMP_WAIT_POLICY", "PASSIVE");
}
// readable diagnosis of a failed "Test MPI"
std::string mpiDiagnosis(const std::string& log, const MpiFamily& fam, int np) {
    std::string d;
    if (log.find("error code 5") != std::string::npos || log.find("Access is denied") != std::string::npos)
        d += "Diagnosis: Windows refused to start ResSimCLI.exe from the MPI launcher (error 5). Checks: (1) ResSimCLI.exe -> Properties -> "
             "'Unblock' (files downloaded from the internet), or copy the release folder out of Downloads (e.g. C:\\ResSim); "
             "(2) antivirus / Controlled Folder Access may block the launcher; (3) Intel MPI older than 2021.6 has a known bug when "
             "mpiexec is started by another program - add '-localonly' to Extra mpiexec arguments or install Microsoft MPI "
             "(recommended, msmpisetup.exe) and leave mpiexec empty.\n";
    int ones = 0, pos = 0;
    while ((pos = int(log.find(" of 1 ", size_t(pos)))) != int(std::string::npos)) { ++ones; ++pos; }
    if (np > 1 && ones >= 2)
        d += "Diagnosis: every rank started its own one-process world - the MPI library loaded by ResSimCLI does not match the launcher (" +
             fam.name + ", expected " + fam.dll + "). Install that runtime or set the mpiexec of the installed MPI.\n";
    if (log.find("libfabric") != std::string::npos || log.find("fi_getinfo") != std::string::npos)
        d += "Diagnosis: the Intel MPI fabric library was not found - install the full Intel MPI runtime or Microsoft MPI.\n";
    return d;
}

struct MpiConfig { std::string mpiexec, cli, args, jobRoot; int procs = 0; };

std::vector<std::string> splitArgs(const std::string& s) {
    std::vector<std::string> out; std::string cur; bool q = false;
    for (char c : s) {
        if (c == '"') { q = !q; continue; }
        if (!q && (c == ' ' || c == '\t')) { if (!cur.empty()) { out.push_back(cur); cur.clear(); } continue; }
        cur += c;
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

// the batch runner registered with the engine's study module
std::vector<ressim::RunResult> mpiRunBatch(const MpiConfig& cfg, const std::vector<ressim::Case>& cases, std::atomic<bool>* cancel,
                                           const std::function<void(int, int)>& progress) {
    std::vector<ressim::RunResult> none;
    if (cfg.mpiexec.empty() || cfg.cli.empty()) return none;
    std::error_code ec;
    fs::path root = cfg.jobRoot.empty() ? fs::temp_directory_path(ec) / "ResSimMPI" : fs::u8path(cfg.jobRoot);
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    fs::path jd = root / ("job_" + std::to_string(plat::processId()) + "_" + std::to_string(stamp % 100000000));
    fs::create_directories(jd, ec);
    if (ec) return none;
    const int n = int(cases.size());
    for (int i = 0; i < n; ++i) {
        char nm[32]; std::snprintf(nm, sizeof(nm), "case_%05d.rsim", i);
        std::ofstream f(jd / nm, std::ios::binary);
        f << ressim::projectText(cases[size_t(i)]);
    }
    { std::ofstream f(jd / "job.txt"); f << "threads 1\n"; }
    int np = cfg.procs > 0 ? cfg.procs : int(std::thread::hardware_concurrency()) + 1;
    np = std::max(1, std::min(np, n + 1));
    std::vector<std::string> args = {"-n", std::to_string(np)};
    for (auto& a : splitArgs(cfg.args)) args.push_back(a);
    args.push_back(cfg.cli); args.push_back("--mpi-batch"); args.push_back(jd.u8string()); args.push_back("-v");
    ChildProcess proc;
    std::string err;
    prepareMpiEnv(cfg.mpiexec);
    if (!proc.start(cfg.mpiexec, args, jd.u8string(), &err)) { fs::remove_all(jd, ec); return none; }
    std::atomic<bool> finished{false};
    std::mutex progMx;   // progress is reported from the watcher and the output pump: serialise (as the local runner does)
    int lastDone = -1;
    auto report = [&](int d, int t) { std::lock_guard<std::mutex> lk(progMx); if (progress && d > lastDone) { lastDone = d; progress(d, t); } };
    std::thread watch([&]() {
        while (!finished.load()) {
            if (cancel && cancel->load()) { std::ofstream(jd / "cancel.txt") << "1"; std::this_thread::sleep_for(std::chrono::seconds(2)); proc.kill(); break; }
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            // progress.txt ("done total np") is written by rank 0 after every case: robust against pipe buffering of -v output
            std::ifstream pf(jd / "progress.txt");
            int d = 0, t = 0;
            if (pf >> d >> t && t > 0) report(d, t);
        }
    });
    std::string log;
    proc.pump([&](const std::string& line) {
        log += line + "\n";
        int idx = 0, rk = 0, d = 0, t = 0;
        if (std::sscanf(line.c_str(), "[mpi] case %d done on rank %d (%d/%d)", &idx, &rk, &d, &t) == 4) report(d, t);
    });
    finished = true;
    watch.join();
    if (!fs::exists(jd / "done.txt", ec) && !(cancel && cancel->load())) {
        // mpiexec / MPI failed as a whole: report through an empty result -> local fallback
        std::ofstream(jd / "mpiexec_log.txt") << log;
        bool any = false;
        for (int i = 0; i < n && !any; ++i) { char nm[32]; std::snprintf(nm, sizeof(nm), "case_%05d.rsres", i); any = fs::exists(jd / nm, ec); }
        if (!any) return none;
    }
    std::vector<ressim::RunResult> out(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        char nm[32]; std::snprintf(nm, sizeof(nm), "case_%05d.rsres", i);
        ressim::RunResult& r = out[size_t(i)];
        try {
            const ressim::ResultsFile rf = ressim::loadResultsFile((jd / nm).u8string(), false);
            r.ok = rf.title == "OK";
            if (!r.ok) r.error = rf.title.rfind("ERROR: ", 0) == 0 ? rf.title.substr(7) : rf.title;
            r.summary = rf.summary; r.stats = rf.stats;
            r.oip0 = r.summary.empty() ? 0 : r.summary.front().oip;
            if (const std::string* kb = rf.blob("KPI")) {   // C28 coupled-mechanics KPIs of the rank's run
                std::istringstream ks(*kb); std::string nmK; double v;
                while (ks >> nmK >> v) r.kpi.push_back({nmK, v});
            }
        } catch (const std::exception& e) { r.ok = false; r.error = cancel && cancel->load() ? "cancelled" : std::string("MPI: ") + e.what(); }
    }
    fs::remove_all(jd, ec);
    return out;
}

} // namespace

void App::applyHpcSettings() {
    if (settings_.mpiEnabled) {
        MpiConfig cfg;
        cfg.mpiexec = findMpiexec(settings_.mpiExec);
        cfg.cli = (fs::u8path(exeFolder()) / (std::string("ResSimCLI") + plat::exeSuffix())).u8string();
        cfg.args = settings_.mpiArgs; cfg.jobRoot = settings_.mpiJobDir; cfg.procs = settings_.mpiProcs;
        std::error_code ec;
        if (!cfg.mpiexec.empty() && fs::exists(fs::u8path(cfg.cli), ec))
            ressim::setBatchRunner([cfg](const std::vector<ressim::Case>& cs, std::atomic<bool>* cancel, const std::function<void(int, int)>& prog) {
                return mpiRunBatch(cfg, cs, cancel, prog);
            });
        else ressim::setBatchRunner(nullptr);
    } else ressim::setBatchRunner(nullptr);
    if (settings_.gpuEnabled) gpu::setPreconditioner(settings_.gpuPrec);
    ressim::setGpuBlockSolver(settings_.gpuEnabled ? gpu::solver() : nullptr);
}

void App::hpcSettingsUi(bool& ch) {
    bool mch = false;
    // ---- C21 performance advisor: what actually speeds up a run on THIS machine (research/CYCLE3_REPORT.md, H3)
    ui::header("Performance advisor (this PC)");
    {
        const int phys = ressim::physicalCores(), logi = int(std::max(1u, std::thread::hardware_concurrency()));
        ImGui::Text("%d physical cores, %d logical processors", phys, logi);
        ui::hint("On one PC, OpenMP threads (Numerics -> Threads, 0 = one per physical core) are the fastest option. "
                 "MPI domain decomposition pays off on several machines or very large models: the Newton assembly stays on rank 0 "
                 "and the slab preconditioner needs ~40 % more iterations (54k cells on 2 cores: threads 33.5 s, DD 36.0 s, "
                 "serial 48.6 s). The GPU solver has no AMG on the device, so its iterations grow with the model size, and an "
                 "integrated laptop GPU shares the memory bandwidth of the CPU - enable it only if 'Benchmark CPU vs GPU' "
                 "shows a gain for your model.");
        if (settings_.mpiEnabled && settings_.mpiDD && settings_.mpiDDProcs > phys)
            ImGui::TextColored(col4(Pal::Warn), "Domain decomposition uses %d ranks on %d physical cores - oversubscribed, expect a slow run.",
                               settings_.mpiDDProcs, phys);
        if (settings_.mpiEnabled && settings_.mpiDD && settings_.mpiArgs.find("-host") == std::string::npos)
            ImGui::TextColored(col4(Pal::Warn), "Domain decomposition on a single PC is usually slower than threads.");
        if (settings_.gpuEnabled) ImGui::TextColored(col4(Pal::Warn), "GPU solver enabled - check 'Benchmark CPU vs GPU' first.");
        ui::hint("High CPU % in the Task Manager while running is partly idle OpenMP threads busy-waiting between parallel "
                 "regions (not extra work). To trade ~10 % speed for lower CPU load and heat, set the Windows environment "
                 "variable OMP_WAIT_POLICY=PASSIVE and restart ResSim. Idle MPI ranks now sleep instead of spinning.");
        if (ImGui::Button("Apply recommended settings for this PC")) {
            settings_.mpiDD = false; settings_.gpuEnabled = false;
            if (case_.tuning.threads > phys) { case_.tuning.threads = 0; markDirty(); }
            mch = true;
        }
    }
    ui::header("Parallel computing (MPI cluster)");
    mch |= ImGui::Checkbox("Run studies / ensembles with MPI (mpiexec + ResSimCLI.exe)", &settings_.mpiEnabled);
    const std::string mpiexec = findMpiexec(settings_.mpiExec);
    std::error_code ec;
    const std::string cli = (fs::u8path(exeFolder()) / (std::string("ResSimCLI") + plat::exeSuffix())).u8string();
    const bool haveCli = fs::exists(fs::u8path(cli), ec);
    if (settings_.mpiEnabled) {
        if (ui::beginProps("mpiopt", 0.36f)) {
            ui::label("mpiexec", "Microsoft MPI (msmpisetup.exe) installs C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe; Intel MPI works too. Empty = automatic.");
            mch |= ui::inputString("##mpiexec", &settings_.mpiExec);
            ui::label("Processes (ranks)", "0 = cores + 1. Rank 0 dispatches the cases, every other rank runs one case at a time.");
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            if (ImGui::InputInt("##mpin", &settings_.mpiProcs)) { settings_.mpiProcs = std::clamp(settings_.mpiProcs, 0, 4096); mch = true; }
            ui::label("Extra mpiexec arguments", "Several machines (MS-MPI): -hosts 2 node1 8 node2 8 (smpd / MS-MPI Launch Service must run on "
                      "every node, ResSimCLI.exe in the same folder on each). Intel MPI: -hosts node1,node2 -ppn 8.");
            mch |= ui::inputString("##mpiargs", &settings_.mpiArgs);
            ui::label("Shared job folder", "Folder for the case files and results; for several machines a share every node can reach "
                      "(\\\\server\\share\\ressim). Empty = local %TEMP%\\ResSimMPI.");
            mch |= ui::inputString("##mpijob", &settings_.mpiJobDir);
            ui::label("Single runs: domain decomposition", "C18: one black-oil simulation over several ranks - the Newton assembly runs on "
                      "rank 0, every linear solve is split into contiguous cell slabs (one per rank) and solved by a distributed "
                      "BiCGSTAB with a CPR preconditioner (global AMG pressure stage + block-Jacobi BILU0). Use at most one rank per "
                      "physical core. Compositional / thermal runs stay local.");
            mch |= ImGui::Checkbox("##mpidd", &settings_.mpiDD);
            if (settings_.mpiDD) {
                ImGui::SameLine();
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
                if (ImGui::InputInt("ranks##mpiddn", &settings_.mpiDDProcs)) { settings_.mpiDDProcs = std::clamp(settings_.mpiDDProcs, 1, 4096); mch = true; }
            }
            if (settings_.mpiDD) {
                ui::label("Assembly", "C28: 'Distributed' = every rank evaluates the cell properties and assembles the Jacobian rows of its "
                          "own cell slab (replicated state, reduced convergence norms, no matrix scatter, all-gathered update); "
                          "SPE9 on 2 ranks: 5.8 s vs 7.3 s with rank-0 assembly and 7.5 s serial. 'Rank 0' = the C18 mode (rank 0 "
                          "assembles, OpenMP threads there). Memory is still replicated per rank (fine up to a few million cells).");
                int m = settings_.mpiDistAsm ? 1 : 0;
                static const char* kAsm[] = {"Rank 0 (C18)", "Distributed (C28)"};
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
                if (ImGui::Combo("##mpida", &m, kAsm, 2)) { settings_.mpiDistAsm = m == 1; mch = true; }
            }
            ui::endProps();
        }
        if (mpiexec.empty()) ImGui::TextColored(col4(Pal::Warn), "mpiexec.exe not found - install Microsoft MPI or set the path. Studies run on the local thread pool.");
        else {
            const MpiFamily fam = mpiFamily(mpiexec);
            ImGui::TextColored(col4(Pal::Ok), "\xE2\x97\x8F mpiexec: %s  (%s)", mpiexec.c_str(), fam.name.c_str());
            if (fam.name == "Intel MPI" && mpiexec.find("Common Files") != std::string::npos)
                ImGui::TextColored(col4(Pal::Warn), "This is the Intel runtime redistributable installed by another program; Microsoft MPI is the tested choice.");
        }
        if (!haveCli) ImGui::TextColored(col4(Pal::Warn), "ResSimCLI.exe is missing next to ResSimStudio.exe.");
        const bool busy = hpcBusy_->load() != 0;
        ImGui::BeginDisabled(busy || mpiexec.empty() || !haveCli);
        if (ImGui::Button("Test MPI")) {
            hpcBusy_->store(1);
            auto busyFlag = hpcBusy_; auto mx = hpcMx_;
            const std::string me = mpiexec, c = cli, extra = settings_.mpiArgs;
            const int np = settings_.mpiProcs > 0 ? settings_.mpiProcs : std::min(4, int(std::thread::hardware_concurrency()) + 1);
            std::string* out = &hpcMpiReport_;
            std::thread([busyFlag, mx, me, c, extra, np, out]() {
                ChildProcess p; std::string err, log;
                std::vector<std::string> args = {"-n", std::to_string(np)};
                for (auto& a : splitArgs(extra)) args.push_back(a);
                args.push_back(c); args.push_back("--mpi-info");
                int rc = -1;
                const MpiFamily fam = mpiFamily(me);
                prepareMpiEnv(me);
                if (p.start(me, args, "", &err)) rc = p.pump([&](const std::string& l) { log += l + "\n"; });
                else log = "could not start mpiexec: " + err;
                const std::string diag = mpiDiagnosis(log, fam, np);
                const bool ok = rc == 0 && diag.empty();
                { std::lock_guard<std::mutex> lk(*mx); *out = (ok ? "OK  (" : "FAILED (exit " + std::to_string(rc) + ", ") + fam.name + ", " + fam.dll + ")\n" + log + diag; }
                busyFlag->store(0);
            }).detach();
        }
        ImGui::EndDisabled();
        if (busy) { ImGui::SameLine(); ImGui::TextDisabled("running..."); }
        {
            std::lock_guard<std::mutex> lk(*hpcMx_);
            if (!hpcMpiReport_.empty()) { ImGui::PushTextWrapPos(ImGui::GetFontSize() * 42); ImGui::TextUnformatted(hpcMpiReport_.c_str()); ImGui::PopTextWrapPos(); }
        }
        ui::hint("MPI distributes many runs (studies, uncertainty, history-matching populations) over cores and machines; with "
                 "'Single runs: domain decomposition' one black-oil run also splits its linear solve over the ranks. Otherwise "
                 "single runs use OpenMP threads on one machine. Results are identical to local runs.");
    }
    ui::header("GPU acceleration");
    bool gch = ImGui::Checkbox("Enable the GPU linear solver (sokol_gfx compute, float32 + double refinement)", &settings_.gpuEnabled);
    if (settings_.gpuEnabled) {
        const char* precs[] = {"Block Jacobi (one dispatch, weakest)", "Multicolour block ILU(0) (default)", "CPR-AMG (device V-cycle, C29)"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 18);
        if (ImGui::Combo("Device preconditioner", &settings_.gpuPrec, precs, 3)) gch = true;
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("CPR-AMG: quasi-IMPES pressure restriction, aggregation AMG V-cycle with damped-Jacobi smoothing\n"
                                                      "and a dense coarsest inverse, block-Jacobi second stage - all as device SpMV / axpy kernels");
        const bool busy = hpcBusy_->load() != 0;
        ImGui::BeginDisabled(busy);
        if (ImGui::Button("Test GPU")) {
            hpcBusy_->store(2);
            auto busyFlag = hpcBusy_; auto mx = hpcMx_;
            std::string* out = &hpcGpuReport_;
            std::thread([busyFlag, mx, out]() {
                bool ok = false;
                const std::string r = gpu::selfTest(ok);
                { std::lock_guard<std::mutex> lk(*mx); *out = (ok ? "PASS  " : "FAIL  ") + r; }
                busyFlag->store(0);
            }).detach();
        }
        ImGui::SameLine();
        if (ImGui::Button("Benchmark CPU vs GPU")) {
            hpcBusy_->store(2);
            auto busyFlag = hpcBusy_; auto mx = hpcMx_;
            std::string* out = &hpcGpuReport_;
            const std::string file = (fs::u8path(exeFolder()) / "gpu_benchmark.txt").u8string();
            std::thread([busyFlag, mx, out, file]() {
                const std::string r = gpu::benchmark(1000000, 0);
                std::ofstream f(fs::u8path(file)); f << r;
                { std::lock_guard<std::mutex> lk(*mx); *out = r + "saved: " + file; }
                busyFlag->store(0);
            }).detach();
        }
        if (busy) { ImGui::SameLine(); ImGui::TextDisabled("running..."); }
        ImGui::EndDisabled();
        {
            std::lock_guard<std::mutex> lk(*hpcMx_);
            if (hpcGpuReport_.find('\n') != std::string::npos) {   // benchmark table: unwrapped, horizontally scrollable
                const float h = ImGui::GetTextLineHeightWithSpacing() * std::min(24.0f, 1.0f + float(std::count(hpcGpuReport_.begin(), hpcGpuReport_.end(), '\n')));
                ImGui::BeginChild("##gpubench", ImVec2(0, h), ImGuiChildFlags_Borders, ImGuiWindowFlags_HorizontalScrollbar);
                ImGui::TextUnformatted(hpcGpuReport_.c_str());
                ImGui::EndChild();
            } else if (!hpcGpuReport_.empty()) { ImGui::PushTextWrapPos(ImGui::GetFontSize() * 42); ImGui::TextUnformatted(hpcGpuReport_.c_str()); ImGui::PopTextWrapPos(); }
        }
        ui::hint("Per case: Numerics -> GPU linear solver. The device runs BiCGSTAB in float32 with a multicolour block ILU(0) "
                 "(cells reordered by a greedy colouring so every colour is one parallel GPU dispatch) or block Jacobi; the "
                 "engine refines the result in double precision, so the Newton solution is unchanged. If the device does not "
                 "converge the CPU solver takes over for that system (counted in the run statistics). Best for large models; "
                 "CPR-AMG stays on the CPU.");
    }
    if (mch || gch) { applyHpcSettings(); ch = true; }
}

// C18: Run -> current black-oil case as ONE domain-decomposed simulation: mpiexec -n N ResSimCLI case.rsim --mpi-dd -o out.rsres
bool App::startMpiDD() {
    if (!settings_.mpiEnabled || !settings_.mpiDD || settings_.mpiDDProcs < 2) return false;
    if (case_.comp.enabled || case_.thermal.enabled) { runner_.addLog("NOTE: MPI domain decomposition from the app applies to black-oil runs - running locally. Compositional / thermal decks: ResSimCLI case.DATA --mpi-dd under mpiexec (summary results)"); return false; }
    const std::string mpiexec = findMpiexec(settings_.mpiExec);
    std::error_code ec;
    const std::string cli = (fs::u8path(exeFolder()) / (std::string("ResSimCLI") + plat::exeSuffix())).u8string();
    if (mpiexec.empty() || !fs::exists(fs::u8path(cli), ec)) {
        runner_.addLog("WARNING: MPI domain decomposition needs mpiexec and ResSimCLI.exe next to ResSimStudio.exe - running locally");
        return false;
    }
    const fs::path dir = fs::u8path(runDirFor(path_, case_.title, "MPI"));
    fs::create_directories(dir, ec);
    const fs::path cs = dir / "case.rsim", out = dir / "case.rsres";
    { std::ofstream f(cs, std::ios::binary); f << ressim::projectText(case_); }
    CliRunSpec s;
    s.exe = mpiexec;
    s.args = {"-n", std::to_string(settings_.mpiDDProcs)};
    for (auto& a : splitArgs(settings_.mpiArgs)) s.args.push_back(a);
    s.args.push_back(cli); s.args.push_back(cs.u8string()); s.args.push_back(settings_.mpiDistAsm ? "--mpi-da" : "--mpi-dd");
    s.args.push_back("-o"); s.args.push_back(out.u8string()); s.args.push_back("-threads");
    s.args.push_back(settings_.mpiDistAsm ? std::string("1")      // C28: one rank per core, each assembles its slab
                                          : std::to_string(std::max(1, std::min(settings_.mpiDDProcs, ressim::physicalCores()))));   // C21 hybrid rank 0
    s.workDir = dir.u8string(); s.resultsFile = out.u8string();
    s.engineName = std::string(settings_.mpiDistAsm ? "ResSim (MPI distributed assembly, " : "ResSim (MPI domain decomposition, ") + std::to_string(settings_.mpiDDProcs) + " ranks)";
    prepareMpiEnv(mpiexec);
    runner_.startCli(case_, s);
    return true;
}
