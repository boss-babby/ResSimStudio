#include "sim_runner.h"
#include "ecl_engine.h"
#include <chrono>
#include <cstdio>

void SimRunner::start(const ressim::Case& c) {
    if (running_) return;
    join();
    cancel_ = false; running_ = true; progress_ = 0; wall_ = 0;
    {
        std::lock_guard<std::mutex> l(m_);
        error_.clear();
        r_ = RunResults{};
    }
    ++version_;
    th_ = std::thread(&SimRunner::work, this, c);
}

void SimRunner::work(ressim::Case c) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    try {
        auto sim = std::make_shared<ressim::Simulator>();
        sim->log = [this](const std::string& s) { std::lock_guard<std::mutex> l(m_); log_.push_back(s); };
        sim->cancel = &cancel_;
        addLog("=== Run started: " + c.title);
        if (!c.extraKeywords.empty()) {
            std::string names; int shown = 0;
            for (auto& e : c.extraKeywords) {
                if (e.section == ressim::DeckSection::Summary || names.find(e.name) != std::string::npos ||
                    e.note.rfind("minimum dimensions", 0) == 0) continue;
                if (shown++ < 12) names += (names.empty() ? "" : " ") + e.name;
            }
            if (!names.empty())
                addLog("NOTE: advanced keywords used only by OPM Flow are not simulated by the ResSim engine: " + names +
                       (shown > 12 ? " ..." : "") + " (switch the engine to OPM Flow to include them)");
        }
        sim->setup(c);
        const double total = sim->totalSimTime();
        sim->onTimestep = [this, total](const ressim::SummaryPoint& s) {
            std::lock_guard<std::mutex> l(m_);
            r_.summary.push_back(s);
            progress_ = total > 0 ? s.time / total : 0;
            ++version_;
        };
        sim->initialize();
        {
            std::lock_guard<std::mutex> l(m_);
            r_.grid = std::make_shared<ressim::Grid>(sim->grid());
            r_.input = c;
            r_.summary = sim->summary();
            r_.snaps = sim->snapshots();
            r_.initial = sim->inPlace();
            r_.totalTime = total;
            for (auto& w : c.wells) { r_.wellNames.push_back(w.name); r_.wellTypes.push_back(w.type); }
            r_.valid = true;
            ++version_;
        }
        while (!cancel_ && sim->runReportStep()) {
            std::lock_guard<std::mutex> l(m_);
            r_.snaps.push_back(sim->snapshots().back());
            r_.stats = sim->stats();
            wall_ = std::chrono::duration<double>(clock::now() - t0).count();
            ++version_;
        }
        const auto& st = sim->stats();
        char buf[256];
        std::snprintf(buf, sizeof(buf), "=== %s: %d timesteps, %d cuts, %d Newton, %d linear iterations, %.2f s",
                      cancel_ ? "Run stopped" : "Run finished", st.timesteps, st.cuts, st.newton, st.linear,
                      std::chrono::duration<double>(clock::now() - t0).count());
        addLog(buf);
        std::lock_guard<std::mutex> l(m_);
        r_.stats = st;
        sim->onTimestep = nullptr; sim->log = nullptr; sim->cancel = nullptr;
        r_.sim = sim;
    } catch (const std::exception& e) {
        std::lock_guard<std::mutex> l(m_);
        error_ = e.what();
        log_.push_back(std::string("ERROR: ") + e.what());
    }
    wall_ = std::chrono::duration<double>(clock::now() - t0).count();
    running_ = false;
    ++version_;
}

// ================================================================================================ OPM Flow engine
#include "ressim/ecl.h"
#include "ressim/io.h"
#include <filesystem>
#include <sstream>

void SimRunner::startOpm(const ressim::Case& c, const OpmOptions& o) {
    if (running_) return;
    join();
    cancel_ = false; running_ = true; progress_ = 0; wall_ = 0;
    {
        std::lock_guard<std::mutex> l(m_);
        error_.clear();
        r_ = RunResults{};
    }
    ++version_;
    th_ = std::thread(&SimRunner::workOpm, this, c, o);
}

void SimRunner::workOpm(ressim::Case c, OpmOptions o) {
    namespace fs = std::filesystem;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    auto fail = [&](const std::string& msg) {
        std::lock_guard<std::mutex> l(m_);
        error_ = msg;
        log_.push_back("ERROR: " + msg);
    };
    const bool ecl = o.kind == 1;
    const std::string engName = ecl ? "ECLIPSE" : "OPM Flow";
    try {
        if (ecl) {
            if (o.flowPath.empty()) o.flowPath = detectEclipseExecutable();
            if (o.flowPath.empty())
                throw std::runtime_error("ECLIPSE (eclrun.exe / eclipse.exe) was not found. Choose the executable in the simulation "
                                         "engine settings (Run menu).");
            if (!fs::is_regular_file(fs::u8path(o.flowPath)))
                throw std::runtime_error("ECLIPSE executable not found: " + o.flowPath + " (simulation engine settings)");
            // direct executable (no eclrun): use the selected program (e300 / frontsim) from the same bin folder
            if (!isEclrunLauncher(o.flowPath) && !o.eclProgram.empty()) {
                const fs::path cur = fs::u8path(o.flowPath);
                const std::string want = o.eclProgram + cur.extension().u8string();
                std::string stem = cur.stem().u8string();
                for (char& ch : stem) ch = char(std::tolower((unsigned char)ch));
                if (stem != o.eclProgram && o.eclProgram != "eclipse") {
                    std::error_code fe;
                    const fs::path alt = cur.parent_path() / want;
                    if (fs::is_regular_file(alt, fe)) o.flowPath = alt.u8string();
                    else throw std::runtime_error(o.eclProgram + " was selected but " + want + " is not next to " + o.flowPath +
                                                  " - choose eclrun.exe (recommended) or the " + o.eclProgram + " executable in the engine settings");
                }
            }
        } else {
            if (o.flowPath.empty()) o.flowPath = detectFlowExecutable();
            if (o.flowPath.empty())
                throw std::runtime_error("OPM Flow (flow.exe) was not found. Install OPM Flow for Windows or choose the executable in "
                                         "Run > OPM Flow settings.");
            if (!fs::is_regular_file(fs::u8path(o.flowPath)))
                throw std::runtime_error("OPM Flow executable not found: " + o.flowPath + " (Run > OPM Flow settings)");
            const std::string pre = flowPrerequisiteProblem(o.flowPath);
            if (!pre.empty()) addLog((pre.rfind("Note:", 0) == 0 ? std::string() : std::string("WARNING: ")) + pre);
            for (auto& n : opmCompatibilityNotes(c)) addLog("NOTE: " + n);
        }
        const fs::path dir = fs::u8path(o.runDir);
        std::error_code ec;
        fs::create_directories(dir, ec);
        if (!fs::is_directory(dir)) throw std::runtime_error("Cannot create the " + engName + " run folder " + o.runDir);
        // Original-deck mode (imported deck not edited in the UI): the source deck and its INCLUDE files are run unchanged,
        // so everything the external simulator supports works even where ResSim's model cannot represent it.
        const bool passThrough = !o.sourceDeck.empty() && fs::is_regular_file(fs::u8path(o.sourceDeck));
        const fs::path src = passThrough ? fs::absolute(fs::u8path(o.sourceDeck), ec) : fs::path();
        std::string base = passThrough ? src.stem().u8string() : caseBaseName(c);
        // ECLIPSE writes next to the dataset (run in the deck folder); OPM writes to --output-dir
        const fs::path outDir = passThrough && ecl ? src.parent_path() : dir;
        // stale output of an earlier run must not be read back
        for (const char* ext : {".SMSPEC", ".UNSMRY", ".UNRST", ".EGRID", ".INIT", ".PRT", ".DBG", ".RSM", ".ESMRY", ".MSG", ".ECLEND", ".LOG", ".h5"})
            fs::remove(outDir / (base + ext), ec);
        const fs::path deck = passThrough ? (ecl ? src : dir / (base + ".DATA")) : dir / (base + ".DATA");
        if (!passThrough) ressim::exportEclipseDeck(c, deck.u8string());
        addLog("=== " + engName + " run started: " + c.title);
        addLog("    engine  " + o.flowPath);
        if (passThrough) {
            addLog("    deck    " + src.u8string() + "  (original imported deck, run unchanged with its INCLUDE files)");
            if (!ecl) addLog("    output  " + dir.u8string());
        } else addLog("    deck    " + deck.u8string() + "  (generated from the model)");
        double total = 0;
        for (auto& s : c.schedule) total += s.length;
        {
            std::lock_guard<std::mutex> l(m_);
            r_.input = c;
            r_.totalTime = total;
            r_.engine = engName;
            r_.opmDeck = deck.u8string();
            for (auto& w : c.wells) { r_.wellNames.push_back(w.name); r_.wellTypes.push_back(w.type); }
        }
        std::vector<std::string> args;
        std::string stdinText;
        if (ecl) {
            // eclrun [options] [-v VERSION] eclipse|e300 CASE   - or the simulator itself, which asks for the case name
            if (isEclrunLauncher(o.flowPath)) {
                { std::istringstream ss(o.extraArgs); std::string a; while (ss >> a) args.push_back(a); }
                if (!o.eclVersion.empty()) { args.push_back("-v"); args.push_back(o.eclVersion); }
                args.push_back(o.eclProgram.empty() ? std::string("eclipse") : o.eclProgram);
                args.push_back(base);
            } else {
                { std::istringstream ss(o.extraArgs); std::string a; while (ss >> a) args.push_back(a); }
                stdinText = base + "\n" + base + "\n";
                addLog("    input   case name '" + base + "' on standard input");
            }
        } else {
            args = {passThrough ? src.u8string() : base + ".DATA", "--output-dir=."};
            int threads = o.threads > 0 ? o.threads : std::min(8, int(std::max(1u, std::thread::hardware_concurrency())));
            args.push_back("--threads-per-process=" + std::to_string(threads));
            if (o.maxDtDays > 0) { std::ostringstream s; s << "--solver-max-time-step-in-days=" << o.maxDtDays; args.push_back(s.str()); }
            if (o.enableTuning) args.push_back("--enable-tuning=true");
            if (o.keepWellsOpen) {   // OPM shuts wells whose equations fail to converge - silently changes the physics
                args.push_back("--shut-unsolvable-wells=false");
                args.push_back("--min-time-step-before-shutting-problematic-wells-in-days=1e-12");
            }
            { std::istringstream ss(o.extraArgs); std::string a; while (ss >> a) args.push_back(a); }
        }
        std::string cmdLine = "    command " + o.flowPath;
        for (auto& a : args) cmdLine += " " + a;
        addLog(cmdLine);

        ChildProcess proc;
        std::string err;
        if (!proc.start(o.flowPath, args, outDir.u8string(), &err, stdinText)) throw std::runtime_error(engName + " " + err + ": " + o.flowPath);
        { std::lock_guard<std::mutex> l(procM_); proc_ = &proc; }
        if (cancel_) proc.kill();
        OpmProgress prog;
        EclProgress eprog;
        const std::string smspec = (outDir / (base + ".SMSPEC")).u8string();
        auto lastLive = clock::now();
        int lastReport = -1;
        auto liveSummary = [&]() {
            // summary files are rewritten at report steps; a half-written file is simply skipped
            ressim::ecl::Summary s; std::string e;
            try {
                if (!ressim::ecl::readSummary(smspec, s, &e) || s.time.empty()) return;
                auto pts = ressim::ecl::runFromSummary(s, r_.wellNames);
                std::lock_guard<std::mutex> l(m_);
                r_.summary = std::move(pts);
                ++version_;
            } catch (...) {}
        };
        const int rc = proc.pump([&](const std::string& line) {
            {
                std::lock_guard<std::mutex> l(m_);
                log_.push_back(line);
                if (!ecl && (line.find("will be shut because it fails to converge") != std::string::npos ||
                    line.find("Problematic well(s) were shut") != std::string::npos))
                    log_.push_back("WARNING: OPM Flow shut a well because its equations did not converge - the results no longer "
                                   "represent the deck. Enable 'Keep problematic wells open' in the OPM engine settings.");
            }
            const bool changed = ecl ? eprog.parse(line) : prog.parse(line);
            if (changed) {
                if (ecl) { prog.day = eprog.day; prog.timesteps = eprog.timesteps; prog.newton = eprog.newton; prog.reportStep = eprog.timesteps; }
                const double tot = !ecl && prog.totalDays > 0 ? prog.totalDays : total;
                if (tot > 0) progress_ = std::min(1.0, prog.day / tot);
                wall_ = std::chrono::duration<double>(clock::now() - t0).count();
                ++version_;
                if (prog.reportStep != lastReport && clock::now() - lastLive > std::chrono::milliseconds(1000)) {
                    lastReport = prog.reportStep; lastLive = clock::now();
                    liveSummary();
                }
            }
        });
        { std::lock_guard<std::mutex> l(procM_); proc_ = nullptr; }
        const double secs = std::chrono::duration<double>(clock::now() - t0).count();
        if (cancel_) addLog("=== " + engName + " run stopped by the user");
        else if (ecl) {
            // eclrun may return 0 although the simulator stopped on input errors: judge by the messages
            const int nErr = std::max(eprog.errors + eprog.bugs, eprog.reportedErrors);
            std::string why = eprog.lastError;
            if (why.empty() && (nErr > 0 || rc != 0)) why = eclipsePrtErrors((outDir / (base + ".PRT")).u8string());
            if (rc != 0 || nErr > 0 || eprog.licenseProblem) {
                std::string msg = "ECLIPSE failed";
                if (rc != 0) msg += " (exit code " + std::to_string(rc) + ")";
                if (nErr > 0) msg += ", " + std::to_string(nErr) + " error(s)";
                if (eprog.licenseProblem) msg += " - no ECLIPSE licence could be checked out (check the licence server / SLBSLS_LICENSE_FILE)";
                if (!why.empty()) msg += ": " + why;
                while (!msg.empty() && msg.back() == '.') msg.pop_back();
                msg += ". See the log and " + base + ".PRT in " + outDir.u8string();
                fail(msg);
            }
        } else if (rc != 0) {
            std::string msg = "OPM Flow failed (exit code " + std::to_string(rc) + ")";
            if (rc == -1073741515 || unsigned(rc) == 0xC0000135u)
                msg += ": a required DLL is missing - most likely the Microsoft MPI runtime (msmpi.dll). Install it with "
                       "'winget install Microsoft.msmpi' or redist\\msmpisetup.exe from the OPM Flow folder";
            else if (!prog.lastError.empty()) msg += ": " + prog.lastError;
            while (!msg.empty() && msg.back() == '.') msg.pop_back();
            msg += ". See the log and " + base + ".PRT in " + outDir.u8string();
            fail(msg);
        }
        // read back whatever was written (complete run, or partial after stop / failure)
        ressim::ecl::ResultCase rcs;
        std::string e;
        if (fs::exists(fs::u8path(smspec)) && ressim::ecl::openResultCase(deck.u8string(), rcs, &e)) {
            auto pts = ressim::ecl::runFromSummary(rcs.summary, r_.wellNames);
            auto snaps = ressim::ecl::snapshotsFromRestart(rcs);
            for (auto& w : rcs.warnings) addLog("WARNING: " + w);
            std::lock_guard<std::mutex> l(m_);
            r_.summary = std::move(pts);
            r_.snaps = std::move(snaps);
            r_.grid = rcs.grid;
            if (!r_.summary.empty()) {
                const auto& s0 = r_.summary.front();
                r_.initial.oil = s0.oip; r_.initial.water = s0.wip; r_.initial.gasFree = s0.gip; r_.initial.gasDissolved = 0;
                const int ol = rcs.summary.find("FOIPL"), og = rcs.summary.find("FOIPG");
                if (ol >= 0 && og >= 0 && !rcs.summary.values[size_t(og)].empty()) {
                    r_.initial.oil = rcs.summary.values[size_t(ol)][0]; r_.initial.oilVaporized = rcs.summary.values[size_t(og)][0];
                }
                const int a = rcs.summary.find("FGIPG"), b = rcs.summary.find("FGIPL");
                if (a >= 0 && b >= 0 && !rcs.summary.values[size_t(a)].empty()) {
                    r_.initial.gasFree = rcs.summary.values[size_t(a)][0]; r_.initial.gasDissolved = rcs.summary.values[size_t(b)][0];
                }
            }
            r_.stats.timesteps = prog.timesteps ? prog.timesteps : int(r_.summary.size());
            r_.stats.newton = prog.newton; r_.stats.linear = prog.linear; r_.stats.seconds = secs;
            r_.valid = r_.grid && !r_.snaps.empty();
            if (!r_.valid && !r_.summary.empty() && r_.grid) r_.valid = true;
            char buf[320];
            std::snprintf(buf, sizeof(buf), "=== %s %s: %d timesteps, %d Newton, %d linear iterations, %.2f s  (%zu report steps read)",
                          engName.c_str(), cancel_ ? "stopped" : error_.empty() ? "finished" : "failed", r_.stats.timesteps, prog.newton, prog.linear, secs, r_.snaps.size());
            log_.push_back(buf);
        } else if (rc == 0 && !cancel_ && error().empty()) {
            fail(engName + " finished but its results could not be read: " + (e.empty() ? smspec : e) +
                 (ecl ? " (ECLIPSE writes them only when the run gets past the input checks - see " + base + ".PRT)" : std::string()));
        }
    } catch (const std::exception& ex) {
        fail(ex.what());
    }
    { std::lock_guard<std::mutex> l(procM_); proc_ = nullptr; }
    wall_ = std::chrono::duration<double>(clock::now() - t0).count();
    progress_ = 1.0;
    running_ = false;
    ++version_;
}
