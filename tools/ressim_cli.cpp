// ResSimCLI - console front end of the ResSim engines (batch / cluster use).
//   ResSimCLI deck.DATA|case.rsim [-o out.rsres] [-threads n] [-q]     run one black-oil case, write a results file
//   ResSimCLI --mpi-info                                             (under mpiexec) print rank / size / host / library
//   mpiexec -n N ResSimCLI --mpi-batch <jobdir> [-v]                 distributed study: case_XXXX.DATA|.rsim -> case_XXXX.rsres
//   mpiexec -n N ResSimCLI deck.DATA --mpi-dd [-o out.rsres]         C18: ONE simulation, linear solves domain-decomposed over
//                                                                    the N ranks (results incl. 3D snapshots in out.rsres)
//   mpiexec -n N ResSimCLI deck.DATA --mpi-da [-o out.rsres]         C28: distributed assembly - every rank evaluates and
//                                                                    assembles its own cell rows (replicated state), collective
//                                                                    solve without matrix scatter (black oil)
// C26 (usability / automation):
//   --json out.json   summary (field vectors per time step, final values, run statistics) as JSON for scripts / notebooks
//   --csv out.csv     field summary table (+ per-well rates and BHP for black-oil runs, C48)
//   --set KEY=VALUE   edit the deck before the run (repeatable): KEY a single-record ECLIPSE keyword (data replaced up to
//                     its '/'), e.g. --set "PERMX=3000*80", or RESSIM:KEY for a '-- RESSIM:' line, e.g.
//                     --set "RESSIM:FAULTSLIP=20 0.25 0.1" - parameter sweeps from a shell script without editing files
//   --check           import + validate only (exit 0 valid, 1 problems), prints every warning
//   --version
// C27:
//   --list-examples              the built-in example library (id, level, grid, summary)
//   --example ID [out.DATA]      write a built-in example as an ECLIPSE deck (default ID.DATA)
//   Coupled black-oil physics run headless like in the GUI: fault mechanics (FAULTMECH / FAULTSLIP), FE geomechanics
//   (GEOMECH), DDM fracture apertures (HFGROW ... DDMWIDTH).
// The MPI library is loaded at run time (MS-MPI msmpi.dll on Windows, MPICH libmpi.so.12 on Linux), so the same
// executable also runs without MPI installed.
#include "ressim/cell_snapshot.h"
#include "ressim/geomodel.h"
#include "ressim/api.h"
#include "ressim/regress.h"
#include "ressim/report.h"
#include "ressim/geocouple.h"
#include "ressim/faultmech.h"
#include "ressim/coupled_deck.h"
#include "ressim/geoseq.h"
#include "ressim/hf_ddm.h"
#include "ressim/coupled_bo.h"
#include "ressim/examples.h"
#include "ressim/case.h"
#include <map>
#include <memory>
#include "ressim/io.h"
#include "ressim/proddata.h"
#include "ressim/study.h"
#include "ressim/results_file.h"
#include "ressim/mpi_rt.h"
#include "ressim/mpi_dd.h"
#include "ressim/simulator.h"
#include "ressim/comp_deck.h"
#include "ressim/therm_deck.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
using namespace ressim;

static int usage() {
    std::printf("ResSimCLI  deck.DATA|case.rsim [-o out.rsres] [-threads n] [-q]\n"
                "           --mpi-info\n"
                "           --mpi-batch <jobdir> [-v]     (start with mpiexec -n N)\n"
                "           deck.DATA|case.rsim --mpi-dd [-o out.rsres]   (start with mpiexec -n N: one domain-decomposed run)\n"
                "           deck.DATA|case.rsim --mpi-da [-o out.rsres]   (mpiexec -n N: distributed assembly + solve, black oil)\n"
                "           deck.DATA|case.rsim --geocouple sxx,syy[,E_GPa,nu,friction,layer]   (C21 G1-c: flow <-> fracture contact coupling; MPa, tension +)\n"
                "           --list-examples | --example ID [out.DATA]   (C27: built-in example library)\n"
                "           --geo-demo [dir] [--field]   (C32 ResGeo: demo geomodelling project + headless workflow -> GRDECL)\n"
                "           --report deck.DATA out.docx|.pptx|.html|.md [--ref REF.SMSPEC] [--label NAME] [--set K=V] [--opt \"title=..;maps=0\"]   (C33)\n"
                "           --regress manifest.txt [-o outdir] [--filter TEXT] [--tol 0.05] [--history file.tsv] [--max-seconds S]   (C33 deck regression vs reference summaries)\n"
                "           --compare REF.SMSPEC OTHER.SMSPEC|deck.DATA [KEY,KEY,...]   (C33 per-vector deviations)\n"
                "           --prod-import table.csv deck.DATA -o out.DATA [--rates] [--ends] [--mdy] [--ctrl RESV|ORAT|LRAT..] [--keep-controls] [--no-wefac]\n"
                "                         [--liq F] [--gas F] [--include out.INC] [--observed obs.csv]   (C38 production data -> QC -> WCONHIST / WCONINJH / WEFAC)\n"
                "  options: --json out.json  --csv out.csv  --fault-csv faults.csv  --set KEY=VALUE (repeatable; KEY = ECLIPSE keyword or RESSIM:KEY)  --check  --version\n");
    return 1;
}

// ---- C26 helpers
static std::string jsonEsc(const std::string& s) {
    std::string o;
    for (char c : s) { if (c == '"' || c == '\\') { o += '\\'; o += c; } else if (c == '\n') o += "\\n"; else if (static_cast<unsigned char>(c) < 32) o += ' '; else o += c; }
    return o;
}
struct Extra { std::string json, csv, engine, title; std::vector<std::pair<std::string, double>> kpi; std::vector<std::string> wells; };   // wells: C48 per-well CSV columns
static Extra g_extra;
static void writeExtras(const std::vector<SummaryPoint>& S, const RunStats* st, double sec) {
    struct Col { const char* name; double SummaryPoint::*m; };
    static const Col cols[] = {{"FOPR", &SummaryPoint::fopr}, {"FWPR", &SummaryPoint::fwpr}, {"FGPR", &SummaryPoint::fgpr}, {"FWIR", &SummaryPoint::fwir},
                               {"FGIR", &SummaryPoint::fgir}, {"FOPT", &SummaryPoint::fopt}, {"FWPT", &SummaryPoint::fwpt}, {"FGPT", &SummaryPoint::fgpt},
                               {"FWIT", &SummaryPoint::fwit}, {"FGIT", &SummaryPoint::fgit}, {"FPR", &SummaryPoint::fpr}, {"FWCT", &SummaryPoint::fwct},
                               {"FGOR", &SummaryPoint::fgor}};
    if (!g_extra.csv.empty()) {
        std::ofstream f(g_extra.csv);
        // C48: + per-well WOPR / WWPR / WGPR / WWIR / WGIR / WBHP columns (black-oil runs; ECLIPSE vector names)
        const auto& W = g_extra.wells;
        f << "TIME"; for (const Col& c : cols) f << ',' << c.name;
        for (const auto& w : W) for (const char* v : {"WOPR", "WWPR", "WGPR", "WWIR", "WGIR", "WBHP"}) f << ',' << v << ':' << w;
        f << '\n';
        for (const auto& p : S) {
            f << p.time; for (const Col& c : cols) f << ',' << p.*(c.m);
            for (size_t k = 0; k < W.size(); ++k) {
                if (k < p.wells.size()) { const WellRates& r = p.wells[k]; f << ',' << r.oil << ',' << r.water << ',' << r.gas << ',' << r.waterInj << ',' << r.gasInj << ',' << r.bhp; }
                else f << ",0,0,0,0,0,0";
            }
            f << '\n';
        }
        std::printf("wrote %s\n", g_extra.csv.c_str());
    }
    if (!g_extra.json.empty()) {
        std::ofstream f(g_extra.json);
        char b[64];
        auto num = [&](double v) { if (!std::isfinite(v)) return std::string("null"); std::snprintf(b, sizeof b, "%.10g", v); return std::string(b); };
        f << "{\n  \"engine\": \"" << jsonEsc(g_extra.engine) << "\",\n  \"title\": \"" << jsonEsc(g_extra.title) << "\",\n  \"seconds\": " << num(sec) << ",\n";
        if (st) f << "  \"stats\": {\"timesteps\": " << st->timesteps << ", \"newton\": " << st->newton << ", \"linear\": " << st->linear << "},\n";
        f << "  \"time\": ["; for (size_t i = 0; i < S.size(); ++i) f << (i ? ", " : "") << num(S[i].time); f << "],\n";
        for (const Col& c : cols) { f << "  \"" << c.name << "\": ["; for (size_t i = 0; i < S.size(); ++i) f << (i ? ", " : "") << num(S[i].*(c.m)); f << "],\n"; }
        f << "  \"final\": {";
        bool first = true;
        if (!S.empty()) for (const Col& c : cols) { f << (first ? "" : ", ") << "\"" << c.name << "\": " << num(S.back().*(c.m)); first = false; }
        for (const auto& k : g_extra.kpi) { f << (first ? "" : ", ") << "\"" << jsonEsc(k.first) << "\": " << num(k.second); first = false; }
        f << "}\n}\n";
        std::printf("wrote %s\n", g_extra.json.c_str());
    }
}
// --set KEY=VALUE on the deck text (ECLIPSE keyword: data up to its terminating '/'; RESSIM:KEY: the '-- RESSIM: KEY' line)
static bool applySet(std::string& text, const std::string& kv, std::string& err) {
    const size_t eq = kv.find('=');
    if (eq == std::string::npos || eq == 0) { err = "--set needs KEY=VALUE: " + kv; return false; }
    std::string key = kv.substr(0, eq), val = kv.substr(eq + 1);
    for (auto& ch : key) ch = char(std::toupper(static_cast<unsigned char>(ch)));
    std::vector<std::string> L;
    { size_t p = 0; while (p <= text.size()) { size_t e = text.find('\n', p); if (e == std::string::npos) e = text.size(); L.push_back(text.substr(p, e - p)); p = e + 1; } }
    auto firstTok = [](const std::string& l) { size_t a = l.find_first_not_of(" \t"); if (a == std::string::npos) return std::string(); size_t b = l.find_first_of(" \t\r/", a); std::string t = l.substr(a, b == std::string::npos ? std::string::npos : b - a); for (auto& ch : t) ch = char(std::toupper(static_cast<unsigned char>(ch))); return t; };
    if (key.rfind("RESSIM:", 0) == 0) {
        const std::string k = key.substr(7);
        for (auto& l : L) {
            if (l.rfind("-- RESSIM:", 0) != 0) continue;
            if (firstTok(l.substr(10)) == k) { l = "-- RESSIM: " + k + " " + val; goto done; }
        }
        {   // not present: add it after the last '-- RESSIM:' line (these lines are read in any section), else before SCHEDULE
            size_t at = L.size();
            for (size_t i = 0; i < L.size(); ++i) if (L[i].rfind("-- RESSIM:", 0) == 0) at = i + 1;
            if (at == L.size()) for (size_t i = 0; i < L.size(); ++i) if (firstTok(L[i]) == "SCHEDULE") { at = i; break; }
            L.insert(L.begin() + long(at), "-- RESSIM: " + k + " " + val);
            goto done;
        }
    }
    for (size_t i = 0; i < L.size(); ++i) {
        if (L[i].rfind("--", 0) == 0 || firstTok(L[i]) != key) continue;
        // the keyword line itself must hold only the keyword
        size_t j = i + 1;
        for (; j < L.size(); ++j) {
            std::string l = L[j]; const size_t c = l.find("--"); if (c != std::string::npos) l = l.substr(0, c);
            if (l.find('/') != std::string::npos) break;
        }
        if (j >= L.size()) { err = "--set: keyword " + key + " has no terminating '/'"; return false; }
        L.erase(L.begin() + long(i) + 1, L.begin() + long(j) + 1);
        L.insert(L.begin() + long(i) + 1, " " + val + " /");
        goto done;
    }
    err = "--set: keyword " + key + " not found in the deck"; return false;
done:
    text.clear();
    for (size_t i = 0; i < L.size(); ++i) { text += L[i]; if (i + 1 < L.size()) text += '\n'; }
    return true;
}

// C26/C27: black-oil run with the coupled mechanics of the GUI runner (fault mechanics incl. FAULTSLIP, FE geomechanics,
// DDM fracture apertures) - the runner itself lives in the engine (coupled_bo.h) and is shared with the studies
static int runCoupledBO(const Case& c, bool quiet, const std::string& out) {
    const auto t0 = std::chrono::steady_clock::now();
    auto log = [&](const std::string& m) { if (!quiet) { std::printf("%s\n", m.c_str()); std::fflush(stdout); } };
    Simulator sim; if (quiet) sim.log = [](const std::string&) {};
    sim.setup(c); sim.initialize();
    coupledbo::Output co;
    coupledbo::run(sim, c, co, log, nullptr, [&](double t) { if (!quiet) { std::printf("  t = %.1f d\n", t); std::fflush(stdout); } });
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const auto& e = sim.summary().back(); const auto& st = sim.stats();
    std::printf("DONE (coupled black oil) t=%.1f d  FOPT %.6g  FWPT %.6g  FWIT %.6g  FPR %.6g  steps %d newton %d linear %d  %.2f s\n", e.time, e.fopt, e.fwpt, e.fwit, e.fpr, st.timesteps, st.newton, st.linear, sec);
    for (const auto& l : co.report) std::printf("%s\n", l.c_str());
    for (const auto& k : co.kpi) g_extra.kpi.push_back(k);
    g_extra.engine = "ResSim coupled black oil"; g_extra.title = c.title;
    writeExtras(sim.summary(), &st, sec);
    if (!out.empty()) {
        ResultsFile rf; rf.engine = "ResSim"; rf.title = c.title.empty() ? "coupled" : c.title; rf.caseText = projectText(c);
        rf.summary = sim.summary(); rf.snaps = sim.snapshots(); rf.stats = st; rf.totalTime = sec; rf.initial = sim.inPlace();
        for (auto& w : c.wells) { rf.wellNames.push_back(w.name); rf.wellTypes.push_back(int(w.type)); }
        saveResultsFile(rf, out);
        std::printf("wrote %s\n", out.c_str());
    }
    return 0;
}

// C20: compositional / CO2STORE and thermal cases (own deck runners), serial or with the MPI domain-decomposed linear
// solve (the DD solver is installed by the caller). Summary-only results file.
static int runDeckEngine(const Case& c, bool quiet, const std::string& out, const char* tag) {
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<SummaryPoint> pts;
    std::string err, eng;
    int steps = 0, newton = 0, linear = 0;
    bool ok = false;
    coupledbo::Output co;
    std::vector<Snapshot> snaps;   // C48: comp / thermal 3-D maps (written with -o)
    auto clog = [&](const std::string& m) { if (!quiet) { std::printf("%s\n", m.c_str()); std::fflush(stdout); } };
    if (c.thermal.enabled) {
        eng = "ResSim thermal";
        thermal::DeckModel m;
        if (!thermal::deckModel(c, m, &err)) { std::printf("ERROR: thermal model: %s\n", err.c_str()); return 3; }
        std::vector<thermal::Summary> o; thermal::Stats st;
        const size_t nAct = m.cellGlobal.size();   // active cells of the grid (= GUI runner numActive)
        std::vector<float> t0;
        auto snapT = [&](const thermal::Simulator& sim, const thermal::Summary& S) {   // C48: 3-D maps for -o (same arrays as the GUI)
            if (!out.empty()) snaps.push_back(thermal::cellSnapshot(sim, c.thermal, S.time, nAct ? nAct : size_t(sim.size()), c.units == UnitSystem::Field, t0));
        };
        if (coupleddeck::needed(c))   // C28: fault mechanics (FAULTMECH / FAULTSLIP, thermo-elastic / THM stress) in the CLI thermal run
            ok = coupleddeck::runThermal(c, m, o, co, clog, [&](const thermal::Summary& S) {
                if (!quiet) { std::printf("  t = %.1f d\n", S.time); std::fflush(stdout); }
                return true; }, &err, -1, &st, snapT);
        else
        ok = thermal::runDeckModel(m, o, [&](const thermal::Summary& S, const thermal::Simulator& sim) {
            snapT(sim, S);
            if (!quiet) { std::printf("  t = %.1f d\n", S.time); std::fflush(stdout); }
            return true; }, &err, -1, &st);
        for (const auto& S : o) {
            SummaryPoint p; p.time = S.time; p.dt = pts.empty() ? S.time : S.time - pts.back().time;
            p.fopr = S.fopr; p.fwpr = S.fwpr; p.fgpr = S.fgpr; p.fwir = S.fwir; p.fgir = S.fgir;
            p.fopt = S.fopt; p.fwpt = S.fwpt; p.fgpt = S.fgpt; p.fwit = S.fwit; p.fgit = S.fgit; p.fpr = S.fpr; p.fprp = S.fpr; p.oip = S.foip;
            pts.push_back(p);
        }
        steps = st.steps; newton = st.newton; linear = st.linear;
    } else {
        eng = c.comp.co2store ? "ResSim CO2STORE" : "ResSim compositional";
        comp::DeckModel m;
        if (!comp::deckModel(c, m, &err)) { std::printf("ERROR: compositional model: %s\n", err.c_str()); return 3; }
        std::vector<comp::CompSummary> o; comp::CompStats st;
        const double pu = c.units == UnitSystem::Field ? 6894.757293168 : 1e5;
        auto snapC = [&](const comp::CompSimulator& sim, const comp::CompSummary& S) {   // C48: 3-D maps for -o
            if (!out.empty()) snaps.push_back(comp::cellSnapshot(sim, S.time, pu));
        };
        if (coupleddeck::needed(c))   // C28: FE geomechanics + fault mechanics (FAULTSLIP) in the CLI compositional run
            ok = coupleddeck::runComp(c, m, o, co, clog, [&](const comp::CompSummary& S) {
                if (!quiet) { std::printf("  t = %.1f d\n", S.time); std::fflush(stdout); }
                return true; }, &err, -1, &st, snapC);
        else
        ok = comp::runDeckModel(m, o, [&](const comp::CompSummary& S, const comp::CompSimulator& sim) {
            snapC(sim, S);
            if (!quiet) { std::printf("  t = %.1f d\n", S.time); std::fflush(stdout); }
            return true; }, &err, -1, &st);
        for (const auto& S : o) {
            SummaryPoint p; p.time = S.time; p.dt = pts.empty() ? S.time : S.time - pts.back().time;
            p.fopr = S.fopr; p.fwpr = S.fwpr; p.fgpr = S.fgpr; p.fwir = S.fwir; p.fgir = S.fgir;
            p.fopt = S.fopt; p.fwpt = S.fwpt; p.fgpt = S.fgpt; p.fwit = S.fwit; p.fgit = S.fgit; p.fpr = S.fpr; p.fprp = S.fpr;
            p.oip = S.foip; p.gip = S.fgip; p.wip = S.fwip;
            pts.push_back(p);
        }
        steps = st.steps; newton = st.newton; linear = st.linear;
    }
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!ok || pts.empty()) { std::printf("ERROR: %s\n", err.empty() ? "run failed" : err.c_str()); return 3; }
    const SummaryPoint& e = pts.back();
    std::printf("DONE (%s%s) t=%.1f d  FOPT %.6g  FWPT %.6g  FGPT %.6g  FPR %.6g  steps %d newton %d linear %d  %.2f s\n", eng.c_str(), tag,
                e.time, e.fopt, e.fwpt, e.fgpt, e.fpr, steps, newton, linear, sec);
    for (const auto& l : co.report) std::printf("%s\n", l.c_str());
    for (const auto& k : co.kpi) g_extra.kpi.push_back(k);
    g_extra.engine = eng; g_extra.title = c.title;
    writeExtras(pts, nullptr, sec);
    if (!out.empty()) {
        ResultsFile rf; rf.engine = eng; rf.title = c.title; rf.caseText = projectText(c);
        rf.summary = pts; rf.totalTime = sec; rf.snaps = std::move(snaps);
        for (auto& w : c.wells) { rf.wellNames.push_back(w.name); rf.wellTypes.push_back(int(w.type)); }
        saveResultsFile(rf, out);
        std::printf("wrote %s\n", out.c_str());
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) return usage();
    const std::string a1 = argv[1];
    if (a1 == "--version") { std::printf("ResSimCLI (ResSim Studio engine) C38\n"); return 0; }
    if (a1 == "--list-examples") {   // C27: the built-in example library (same as File > Examples in the app)
        const auto L = listExamples();
        std::printf("%-22s %-13s %-28s %s\n", "ID", "LEVEL", "GRID", "SUMMARY");
        for (const auto& e : L) std::printf("%-22s %-13s %-28s %s\n", e.id.c_str(), toString(e.level), e.grid.c_str(), e.summary.c_str());
        std::printf("%zu examples - write one as an ECLIPSE deck with: ResSimCLI --example ID [out.DATA]\n", L.size());
        return 0;
    }
    if (a1 == "--geo-demo") {   // C32 ResGeo: write the demo geomodelling project and run the whole workflow headless
        const std::string dir = argc > 2 && argv[2][0] != '-' ? argv[2] : std::string("resgeo_demo");
        bool field = false;
        for (int i = 2; i < argc; ++i) if (!std::strcmp(argv[i], "--field")) field = true;
        std::string err;
        const std::string rep = ressim::geo::runDemoWorkflow(dir, field, err);
        if (rep.empty()) { std::fprintf(stderr, "geo demo failed: %s\n", err.c_str()); return 1; }
        std::printf("%s", rep.c_str());
        return 0;
    }
    if (a1 == "--report") {   // C33: run a deck and write a Word / PowerPoint / HTML / Markdown report
        if (argc < 4) return usage();
        std::string ref, label = "reference", opts;
        std::vector<std::string> sets;
        for (int i = 4; i < argc; ++i) {
            if (!std::strcmp(argv[i], "--ref") && i + 1 < argc) ref = argv[++i];
            else if (!std::strcmp(argv[i], "--label") && i + 1 < argc) label = argv[++i];
            else if (!std::strcmp(argv[i], "--set") && i + 1 < argc) sets.push_back(argv[++i]);
            else if (!std::strcmp(argv[i], "--opt") && i + 1 < argc) opts = argv[++i];
        }
        api::Session S;
        if (!S.open(argv[2])) { std::printf("ERROR: %s\n", S.error().c_str()); return 2; }
        for (const auto& kv : sets) if (!S.setKeyword(kv)) { std::printf("ERROR: %s\n", S.error().c_str()); return 2; }
        std::printf("running %s (%s) ...\n", argv[2], api::engineName(S.caseData()).c_str()); std::fflush(stdout);
        if (!S.run([](double t) { std::printf("  t = %.1f d\n", t); std::fflush(stdout); })) { std::printf("ERROR: %s\n", S.error().c_str()); return 3; }
        ecl::Summary R; const ecl::Summary* rp = nullptr;
        if (!ref.empty()) { std::string err; if (!ecl::readSummary(ref, R, &err)) { std::printf("ERROR: reference %s: %s\n", ref.c_str(), err.c_str()); return 2; } rp = &R; }
        const report::Doc d = report::simulationReport(S.caseData(), S.output(), S.summary(), report::parseOptions(opts), rp, label);
        std::string err;
        if (!report::writeDoc(argv[3], d, err)) { std::printf("ERROR: %s\n", err.c_str()); return 2; }
        std::printf("DONE %.2f s; wrote %s\n", S.output().seconds, argv[3]);
        return 0;
    }
    if (a1 == "--prod-import") {   // C38: production table -> QC -> WCONHIST / WCONINJH / WEFAC merged into a deck
        // ResSimCLI --prod-import table.csv deck.DATA -o out.DATA [--rates] [--ends] [--mdy] [--ctrl RESV] [--keep-controls]
        //           [--no-wefac] [--liq F] [--gas F] [--include out.INC] [--observed obs.csv]
        if (argc < 4) return usage();
        const std::string tab = argv[2], deckPath = argv[3];
        std::string out, inc, obs;
        pdata::Options o;
        for (int i = 4; i < argc; ++i) {
            if (!std::strcmp(argv[i], "-o") && i + 1 < argc) out = argv[++i];
            else if (!std::strcmp(argv[i], "--rates")) o.volumes = false;
            else if (!std::strcmp(argv[i], "--ends")) o.datesAreEnds = true;
            else if (!std::strcmp(argv[i], "--mdy")) o.monthFirst = true;
            else if (!std::strcmp(argv[i], "--ctrl") && i + 1 < argc) o.prodControl = argv[++i];
            else if (!std::strcmp(argv[i], "--keep-controls")) o.replaceControls = false;
            else if (!std::strcmp(argv[i], "--no-wefac")) o.useWefac = false;
            else if (!std::strcmp(argv[i], "--liq") && i + 1 < argc) o.liqFactor = std::atof(argv[++i]);
            else if (!std::strcmp(argv[i], "--gas") && i + 1 < argc) o.gasFactor = std::atof(argv[++i]);
            else if (!std::strcmp(argv[i], "--include") && i + 1 < argc) inc = argv[++i];
            else if (!std::strcmp(argv[i], "--observed") && i + 1 < argc) obs = argv[++i];
        }
        auto slurp = [](const std::string& p) { std::ifstream f(p, std::ios::binary); std::ostringstream ss; ss << f.rdbuf(); return ss.str(); };
        const std::string text = slurp(tab), deck = slurp(deckPath);
        if (text.empty()) { std::fprintf(stderr, "cannot read %s\n", tab.c_str()); return 2; }
        std::vector<std::string> wells;
        try { Case c = importEclipseDeck(deckPath); o.startYear = c.startYear; o.startMonth = c.startMonth; o.startDay = c.startDay; for (auto& w : c.wells) wells.push_back(w.name); }
        catch (const std::exception& e) { std::fprintf(stderr, "deck: %s\n", e.what()); return 2; }
        pdata::Table t;
        const bool ok = pdata::parseTable(text, o, t);
        static const char* sev[] = {"info", "WARNING", "ERROR"};
        for (const auto& i : t.issues) std::printf("%-7s %5s %-8s %-10s %s\n", sev[i.severity], i.row ? std::to_string(i.row).c_str() : "", i.well.c_str(), i.date.c_str(), i.msg.c_str());
        if (!ok) return 1;
        if (!inc.empty()) { std::ofstream(inc, std::ios::binary) << pdata::historyInclude(t, o); std::printf("wrote %s\n", inc.c_str()); }
        if (!obs.empty()) { std::ofstream(obs, std::ios::binary) << pdata::observedCsv(t, o); std::printf("wrote %s\n", obs.c_str()); }
        if (!out.empty()) {
            std::vector<pdata::Issue> mi;
            const std::string merged = pdata::mergeIntoDeck(deck, t, o, mi, wells);
            bool fatal = false;
            for (const auto& i : mi) { std::printf("%-7s merge %-8s %s\n", sev[i.severity], i.well.c_str(), i.msg.c_str()); fatal |= i.severity == 2; }
            if (fatal) return 1;
            std::ofstream(out, std::ios::binary) << merged;
            std::printf("wrote %s\n", out.c_str());
        }
        return 0;
    }
    if (a1 == "--regress") {   // C33: deck regression against reference summaries (nightly opm-tests run)
        if (argc < 3) return usage();
        regress::Manifest m; regress::Options o;
        for (int i = 3; i < argc; ++i) {
            if (!std::strcmp(argv[i], "-o") && i + 1 < argc) o.outDir = argv[++i];
            else if (!std::strcmp(argv[i], "--filter") && i + 1 < argc) o.filter = argv[++i];
            else if (!std::strcmp(argv[i], "--tol") && i + 1 < argc) o.tolOverride = std::atof(argv[++i]);
            else if (!std::strcmp(argv[i], "--history") && i + 1 < argc) o.history = argv[++i];
            else if (!std::strcmp(argv[i], "--max-seconds") && i + 1 < argc) o.maxSeconds = std::atof(argv[++i]);
        }
        std::string err;
        if (!regress::readManifest(argv[2], m, err)) { std::printf("ERROR: %s\n", err.c_str()); return 2; }
#ifdef _WIN32
        std::system(("mkdir \"" + o.outDir + "\" 2>nul").c_str());
#else
        std::system(("mkdir -p '" + o.outDir + "'").c_str());
#endif
        o.log = [](const std::string& l) { std::printf("%s\n", l.c_str()); std::fflush(stdout); };
        const auto R = regress::runAll(m, o);
        const int bad = regress::writeOutputs(R, o, err);
        if (bad < 0) { std::printf("ERROR: %s\n", err.c_str()); return 2; }
        int np = 0, ns = 0; for (const auto& r : R) { np += r.status == "PASS"; ns += r.status == "SKIP"; }
        std::printf("REGRESSION %zu cases: %d pass, %d fail/error, %d skipped -> %s/report.html, report.docx, junit.xml\n", R.size(), np, bad, ns, o.outDir.c_str());
        return bad > 0 ? 1 : 0;
    }
    if (a1 == "--compare") {   // C33: deviations of a summary (or a ResSim run of a deck) from a reference summary
        if (argc < 4) return usage();
        ecl::Summary A, B; std::string err;
        if (!ecl::readSummary(argv[2], A, &err)) { std::printf("ERROR: %s\n", err.c_str()); return 2; }
        std::string b = argv[3], bu = b; for (auto& ch : bu) ch = char(std::toupper((unsigned char)ch));
        if (bu.size() > 7 && bu.substr(bu.size() - 7) == ".SMSPEC") { if (!ecl::readSummary(b, B, &err)) { std::printf("ERROR: %s\n", err.c_str()); return 2; } }
        else {
            api::Session S;
            if (!S.open(b) || !S.run()) { std::printf("ERROR: %s\n", S.error().c_str()); return 3; }
            B = S.summary();
        }
        std::vector<std::string> keys;
        if (argc > 4) { std::string k = argv[4]; size_t a = 0; while (a <= k.size()) { size_t e = k.find(',', a); if (e == std::string::npos) e = k.size(); if (e > a) keys.push_back(k.substr(a, e - a)); a = e + 1; } }
        const auto D = ecl::compareSummaries(A, B, keys);
        std::printf("%-22s %14s %14s %10s\n", "VECTOR", "REF FINAL", "OTHER FINAL", "relRMS %");
        for (const auto& d : D) std::printf("%-22s %14.6g %14.6g %10.4f\n", d.key.c_str(), d.refFinal, d.cmpFinal, 100 * d.relRms);
        return 0;
    }
    if (a1 == "--example") {
        if (argc < 3) return usage();
        try {
            const Case c = makeExample(argv[2]);
            std::string path = argc > 3 ? argv[3] : std::string();
            if (path.empty()) { path = argv[2]; for (auto& ch : path) ch = char(std::toupper((unsigned char)ch)); path += ".DATA"; }
            exportEclipseDeck(c, path);
            std::printf("wrote %s\n", path.c_str());
            return 0;
        } catch (const std::exception& ex) { std::printf("ERROR: %s\n", ex.what()); return 2; }
    }
    bool verbose = false;
    for (int i = 2; i < argc; ++i) if (!std::strcmp(argv[i], "-v")) verbose = true;
    if (a1 == "--mpi-info" || a1 == "--mpi-batch") {
        std::string err;
        if (!mpi::load(&err)) { std::fprintf(stderr, "MPI: %s\n", err.c_str()); return 2; }
        if (!mpi::init()) { std::fprintf(stderr, "MPI: MPI_Init failed (%s)\n", mpi::libraryName().c_str()); return 2; }
        int rc = 0;
        if (a1 == "--mpi-info") {
            const double s = mpi::allreduceSum(1.0);
            std::printf("rank %d of %d on %s (%s), allreduce check %g\n", mpi::rank(), mpi::size(), mpi::processorName().c_str(),
                        mpi::libraryName().c_str(), s);
            std::fflush(stdout);
            rc = (int(s + 0.5) == mpi::size()) ? 0 : 4;
            // C17: launched by an MPI whose process manager reports more ranks than this library sees -> wrong runtime DLL
            if (const char* ps = std::getenv("PMI_SIZE"); ps && std::atoi(ps) > 1 && mpi::size() == 1) {
                std::printf("MISMATCH: the launcher started %s ranks but %s sees a world of 1 - set RESSIM_MPI_LIB to the launcher's MPI library\n",
                            ps, mpi::libraryName().c_str());
                rc = 5;
            }
        } else {
            if (argc < 3) { mpi::finalize(); return usage(); }
            rc = mpi::runBatchJob(argv[2], verbose);
        }
        mpi::barrier();
        mpi::finalize();
        return rc;
    }
    // single run
    std::string out, geoArg; int threads = 0; bool quiet = false, dd = false, checkOnly = false;
    std::vector<std::string> sets;
    bool da = false;
    for (int i = 2; i < argc; ++i) { if (!std::strcmp(argv[i], "--mpi-dd")) dd = true; if (!std::strcmp(argv[i], "--mpi-da")) da = true; }
    if (da) dd = false;
    if (da) {   // C28: distributed assembly - every rank runs the Newton loop on its own cell rows (no serving ranks)
        std::string err;
        if (!mpi::load(&err)) { std::fprintf(stderr, "MPI: %s\n", err.c_str()); return 2; }
        if (!mpi::init()) { std::fprintf(stderr, "MPI: MPI_Init failed (%s)\n", mpi::libraryName().c_str()); return 2; }
    }
    if (dd) {
        std::string err;
        if (!mpi::load(&err)) { std::fprintf(stderr, "MPI: %s\n", err.c_str()); return 2; }
        if (!mpi::init()) { std::fprintf(stderr, "MPI: MPI_Init failed (%s)\n", mpi::libraryName().c_str()); return 2; }
        if (mpi::rank() > 0) { const int rc = mpi::serveDD(); mpi::finalize(); return rc; }
    }
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "-o") && i + 1 < argc) out = argv[++i];
        else if (!std::strcmp(argv[i], "-threads") && i + 1 < argc) threads = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-q")) quiet = true;
        else if (!std::strcmp(argv[i], "--geocouple") && i + 1 < argc) geoArg = argv[++i];
        else if (!std::strcmp(argv[i], "--json") && i + 1 < argc) g_extra.json = argv[++i];
        else if (!std::strcmp(argv[i], "--csv") && i + 1 < argc) g_extra.csv = argv[++i];
        else if (!std::strcmp(argv[i], "--set") && i + 1 < argc) sets.push_back(argv[++i]);
        else if (!std::strcmp(argv[i], "--check")) checkOnly = true;
        else if (!std::strcmp(argv[i], "--fault-csv") && i + 1 < argc) {   // C47: per-face fault history (FAULTMECH runs)
#ifdef _WIN32
            _putenv_s("RESSIM_FAULT_CSV", argv[++i]);
#else
            setenv("RESSIM_FAULT_CSV", argv[++i], 1);
#endif
        }
    }
    try {
        std::string ext = a1.size() > 5 ? a1.substr(a1.size() - 5) : "";
        for (auto& ch : ext) ch = char(std::tolower((unsigned char)ch));
        std::vector<std::string> warn;
        Case c = (ext == ".rsim") ? loadProject(a1) : importEclipseDeck(a1, &warn);
        if (!quiet || checkOnly) for (auto& w : warn) std::printf("WARN: %s\n", w.c_str());
        if (!sets.empty()) {   // C26: deck edits from the command line
            std::string text = projectText(c), err;
            for (const auto& kv : sets) { if (!applySet(text, kv, err)) { std::printf("ERROR: %s\n", err.c_str()); return 2; } if (!quiet) std::printf("set %s\n", kv.c_str()); }
            { const size_t sl = a1.find_last_of("\\/"); c = projectFromText(text, sl == std::string::npos ? std::string() : a1.substr(0, sl + 1)); }
        }
        if (checkOnly) {
            const std::string v = validateCase(c);
            std::printf("%s: %d x %d x %d%s, %zu wells, %zu report steps, %zu warnings\n", a1.c_str(), c.grid.nx, c.grid.ny, c.grid.nz, c.grid.cornerPoint() ? " corner-point" : "",
                        c.wells.size(), c.schedule.size(), warn.size());
            if (!c.engineUnsupported.empty()) for (const auto& u : c.engineUnsupported) std::printf("OPM/ECLIPSE only: %s\n", u.c_str());
            std::printf("%s\n", v.empty() ? "VALID" : ("PROBLEMS: " + v).c_str());
            return v.empty() ? 0 : 1;
        }
        if (threads > 0) c.tuning.threads = threads;
        if (da) {
            const int me = mpi::rank(), P = mpi::size();
            int rc = 0;
            if (c.thermal.enabled || c.comp.enabled) {
                if (me == 0) std::printf("ERROR: --mpi-da is for the black-oil engine (compositional / thermal: use --mpi-dd)\n");
                mpi::finalize(); return 2;
            }
            if (threads <= 0) c.tuning.threads = 1;          // one rank per core
            c.tuning.blockSolver = true;
            mpi::DDSolver solver;
            const auto t0 = std::chrono::steady_clock::now();
            try {
                Simulator sim;
                if (me > 0) sim.log = [](const std::string&) {};
                sim.setup(c);
                DistAssembly d;
                d.starts = sim.cellPartition(P);
                d.a0 = d.starts[size_t(me)]; d.a1 = d.starts[size_t(me) + 1];
                d.sumN = [](double* v, int n) { mpi::allreduceSumN(v, n); };
                d.maxN = [](double* v, int n) { for (int i = 0; i < n; ++i) v[i] = mpi::allreduceMax(v[i]); };
                const std::vector<int> st = d.starts;
                d.solve = [&solver, st](const BSR& S, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt, int& its) {
                    return solver.solveCollective(S, b, x, st, tol, maxIt, its);
                };
                sim.setDistributedAssembly(d);
                sim.initialize();
                while (sim.runReportStep()) {
                    if (me == 0 && !quiet && !sim.summary().empty()) { std::printf("  t = %.1f d\n", sim.summary().back().time); std::fflush(stdout); }
                }
                const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                const auto& st2 = sim.stats();
                const double asmMax = mpi::allreduceMax(st2.secAssemble), asmMin = -mpi::allreduceMax(-st2.secAssemble);
                if (me == 0) {
                    const auto& e = sim.summary().back();
                    std::printf("DONE (MPI distributed assembly, %d ranks, own rows %d of %d) t=%.1f d  FOPT %.6g  FWPT %.6g  FPR %.6g  steps %d newton %d linear %d  %.2f s\n",
                                P, d.a1 - d.a0, d.starts.back(), e.time, e.fopt, e.fwpt, e.fpr, st2.timesteps, st2.newton, st2.linear, sec);
                    std::printf("  assembly per rank %.2f - %.2f s; collective solves %d: set-up %.2f s, Krylov %.2f s, gather %.2f s, halo %lld doubles\n",
                                asmMin, asmMax, solver.stats.solves, solver.stats.scatterSec, solver.stats.solveSec, solver.stats.haloSec, solver.stats.haloDoubles);
                    if (!out.empty()) {
                        ResultsFile rf; rf.engine = "ResSim"; rf.title = c.title.empty() ? "MPI distributed assembly" : c.title; rf.caseText = projectText(c);
                        rf.summary = sim.summary(); rf.snaps = sim.snapshots(); rf.stats = st2; rf.totalTime = sec; rf.initial = sim.inPlace();
                        for (auto& w : c.wells) { rf.wellNames.push_back(w.name); rf.wellTypes.push_back(int(w.type)); }
                        saveResultsFile(rf, out);
                        std::printf("wrote %s\n", out.c_str());
                    }
                }
            } catch (const std::exception& ex) { std::printf("ERROR (rank %d): %s\n", me, ex.what()); rc = 3; }
            mpi::finalize();
            return rc;
        }
        if (dd) {
            // rank 0: Newton loop here, every linear solve shared with the serving ranks
            mpi::DDSolver solver;
            setDistributedSolver(&solver);
            int rc = 0;
            // C21 hybrid: the serving ranks sleep during the Newton assembly, so rank 0 may use one thread per rank for it
            if (threads <= 0) c.tuning.threads = std::max(1, std::min(mpi::size(), physicalCores()));
            if (c.tuning.threads > 1 && !quiet) {
                const char* wp = std::getenv("OMP_WAIT_POLICY");
                if (!wp || (std::strcmp(wp, "PASSIVE") && std::strcmp(wp, "passive")))
                    std::printf("NOTE: set OMP_WAIT_POLICY=PASSIVE for --mpi-dd with -threads > 1 (spinning OpenMP threads slow the distributed solve)\n");
            }
            const auto t0 = std::chrono::steady_clock::now();
            if (c.thermal.enabled || c.comp.enabled) {   // C20: comp / thermal engines share the same distributed solver
                try { rc = runDeckEngine(c, quiet, out, ", MPI domain decomposition"); }
                catch (const std::exception& ex) { std::printf("ERROR: %s\n", ex.what()); rc = 3; }
                std::printf("  distributed solves %d (%d iterations): %.2f s, scatter %.2f s, halo %.2f s, %d pattern set-ups, %d ranks\n",
                            solver.stats.solves, solver.stats.iterations, solver.stats.solveSec, solver.stats.scatterSec, solver.stats.haloSec,
                            solver.stats.patternSetups, solver.stats.ranks);
                setDistributedSolver(nullptr);
                solver.shutdown();
                mpi::finalize();
                return rc;
            }
            try {
                Simulator sim; sim.setup(c); sim.initialize();
                while (sim.runReportStep()) {
                    if (!quiet && !sim.summary().empty()) { std::printf("  t = %.1f d\n", sim.summary().back().time); std::fflush(stdout); }
                }
                const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                const auto& e = sim.summary().back();
                const auto& st = sim.stats();
                std::printf("DONE (MPI domain decomposition, %d ranks) t=%.1f d  FOPT %.6g  FWPT %.6g  FPR %.6g  steps %d newton %d linear %d  %.2f s\n",
                            solver.stats.ranks, e.time, e.fopt, e.fwpt, e.fpr, st.timesteps, st.newton, st.linear, sec);
                std::printf("  distributed solves %d (CPU fallbacks %d): %.2f s, scatter %.2f s, halo %.2f s (%lld doubles), %d pattern set-ups\n", st.ddSolves,
                            st.ddFallbacks, solver.stats.solveSec, solver.stats.scatterSec, solver.stats.haloSec, solver.stats.haloDoubles, solver.stats.patternSetups);
                std::printf("  profile: assemble %.2f s, precond %.2f s, solve %.2f s, update %.2f s, distributed %.2f s, device %.2f s\n", st.secAssemble, st.secPrecond, st.secSolve, st.secUpdate, st.ddSec, st.gpuWallSec);
                if (!out.empty()) {
                    ResultsFile rf; rf.engine = "ResSim"; rf.title = c.title.empty() ? "MPI domain decomposition" : c.title; rf.caseText = projectText(c);
                    rf.summary = sim.summary(); rf.snaps = sim.snapshots(); rf.stats = st; rf.totalTime = sec; rf.initial = sim.inPlace();
                    for (auto& w : c.wells) { rf.wellNames.push_back(w.name); rf.wellTypes.push_back(int(w.type)); }
                    saveResultsFile(rf, out);
                    std::printf("wrote %s\n", out.c_str());
                }
            } catch (const std::exception& ex) { std::printf("ERROR: %s\n", ex.what()); rc = 3; }
            setDistributedSolver(nullptr);
            solver.shutdown();
            mpi::finalize();
            return rc;
        }
        if (!geoArg.empty()) {
            // C21 G1-c: staggered EDFM flow <-> EFEM contact coupling (sxx,syy MPa [,E GPa, nu, friction, layer 1-based])
            double v[6] = {-40, -30, 20, 0.25, 0.6, 1}; int nv = 0;
            for (size_t a = 0; a <= geoArg.size() && nv < 6;) { size_t b = geoArg.find(',', a); if (b == std::string::npos) b = geoArg.size(); v[nv++] = std::atof(geoArg.substr(a, b - a).c_str()); a = b + 1; }
            geocouple::Input gi; gi.sxx = v[0] * 1e6; gi.syy = v[1] * 1e6; gi.E = v[2] * 1e9; gi.nu = v[3]; gi.friction = v[4]; gi.layer = std::max(0, int(v[5]) - 1);
            const auto t0 = std::chrono::steady_clock::now();
            Simulator sim; if (quiet) sim.log = [](const std::string&) {};
            sim.setup(c); sim.initialize();
            geocouple::Output go;
            const bool ok = geocouple::runCoupled(sim, gi, go, nullptr, [&](const std::string& m) { if (!quiet) { std::printf("%s\n", m.c_str()); std::fflush(stdout); } });
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (!ok) { std::printf("ERROR: %s\n", go.error.c_str()); return 3; }
            const auto& e = sim.summary().back(); const auto& st = sim.stats(); const auto& L = go.steps.back();
            std::printf("DONE (G1-c coupled) t=%.1f d  FOPT %.6g  FWPT %.6g  FWIT %.6g  FPR %.6g  steps %d newton %d  %.2f s (contact %.2f s)\n", e.time, e.fopt, e.fwpt, e.fwit, e.fpr, st.timesteps, st.newton, sec, go.efemSeconds);
            std::printf("  last coupling: %d open / %d slipping pieces, max opening %.4g mm, max slip %.4g mm, conductance x %.4g .. %.4g\n", L.openPieces, L.slipPieces, L.maxOpening * 1e3, L.maxSlip * 1e3, L.minMult, L.maxMult);
            if (!out.empty()) {
                ResultsFile rf; rf.engine = "ResSim + G1-c"; rf.title = c.title.empty() ? "G1-c coupled" : c.title; rf.caseText = projectText(c);
                rf.summary = sim.summary(); rf.snaps = sim.snapshots(); rf.stats = st; rf.totalTime = sec; rf.initial = sim.inPlace();
                for (auto& w : c.wells) { rf.wellNames.push_back(w.name); rf.wellTypes.push_back(int(w.type)); }
                saveResultsFile(rf, out);
                std::printf("wrote %s\n", out.c_str());
            }
            return 0;
        }
        if (c.thermal.enabled || c.comp.enabled) return runDeckEngine(c, quiet, out, "");
        if (coupledbo::needed(c)) return runCoupledBO(c, quiet, out);
        const auto t0 = std::chrono::steady_clock::now();
        RunResult r = runSingle(c, nullptr);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (!r.ok) { std::printf("ERROR: %s\n", r.error.c_str()); return 3; }
        const auto& e = r.summary.back();
        std::printf("DONE t=%.1f d  FOPT %.6g  FWPT %.6g  FGPT %.6g  FPR %.6g  steps %d newton %d linear %d  %.2f s\n", e.time, e.fopt,
                    e.fwpt, e.fgpt, e.fpr, r.stats.timesteps, r.stats.newton, r.stats.linear, sec);
        std::printf("  profile: assemble %.2f s, precond %.2f s, solve %.2f s, update %.2f s, distributed %.2f s, device %.2f s\n", r.stats.secAssemble, r.stats.secPrecond, r.stats.secSolve, r.stats.secUpdate, r.stats.ddSec, r.stats.gpuWallSec);
        g_extra.engine = "ResSim"; g_extra.title = c.title;
        g_extra.wells.clear(); for (const auto& w : c.wells) g_extra.wells.push_back(w.name);
        writeExtras(r.summary, &r.stats, sec);
        if (!out.empty()) {
            ResultsFile rf; rf.engine = "ResSim"; rf.title = "OK"; rf.caseText = projectText(c);
            rf.summary = r.summary; rf.stats = r.stats; rf.totalTime = sec;
            saveResultsFile(rf, out);
            std::printf("wrote %s\n", out.c_str());
        }
    } catch (const std::exception& e) {
        std::printf("ERROR: %s\n", e.what());
        if (dd) { mpi::DDSolver().shutdown(); mpi::finalize(); }
        return 2;
    }
    return 0;
}
