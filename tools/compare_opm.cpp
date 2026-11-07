// compare_opm: run a deck with ResSim and compare its summary with an OPM Flow / ECLIPSE result case.
//   compare_opm <deck.DATA> <opm result file or base> [-t seconds] [-csv out.csv] [-q]
//   compare_opm -dump <file>          list the arrays of an ECLIPSE file (reader test)
//   compare_opm -info <result file>   summary of a result case (grid, steps, vectors)
#include "ressim/ecl.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include "ressim/grid.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace ressim;

static int dump(const char* path) {
    ecl::File f; std::string err;
    if (!f.open(path, &err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 1; }
    std::printf("%s (%s), %zu arrays\n", path, f.formatted() ? "formatted" : "binary", f.arrays().size());
    int shown = 0;
    for (size_t i = 0; i < f.arrays().size(); ++i) {
        const auto& a = f.arrays()[i];
        if (++shown > 60) { std::printf("  ...\n"); break; }
        std::printf("  %-8s %-4s %10lld", a.name.c_str(), a.typeName.c_str(), (long long)a.count);
        if (a.count > 0 && a.count <= 6) {
            if (a.type == ecl::Type::CHAR) for (auto& s : f.strings(int(i))) std::printf(" '%s'", s.c_str());
            else for (double v : f.numeric(int(i))) std::printf(" %g", v);
        }
        std::printf("\n");
    }
    return 0;
}

static int info(const char* path) {
    ecl::ResultCase rc; std::string err;
    if (!ecl::openResultCase(path, rc, &err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 1; }
    std::printf("case %s\n", rc.name.c_str());
    if (rc.hasGrid) std::printf("  grid %dx%dx%d, %d active, bbox x %.1f..%.1f y %.1f..%.1f z %.1f..%.1f\n", rc.egrid.nx, rc.egrid.ny, rc.egrid.nz, rc.grid->numActive(),
                                rc.grid->bbMin.x, rc.grid->bbMax.x, rc.grid->bbMin.y, rc.grid->bbMax.y, rc.grid->bbMin.z, rc.grid->bbMax.z);
    if (rc.hasInit) { std::printf("  INIT arrays:"); for (auto& s : rc.initArrays) std::printf(" %s", s.c_str()); std::printf("\n"); }
    if (rc.hasRestart) {
        const auto& st = rc.restart.steps();
        std::printf("  restart: %zu steps (report %d..%d, day %.1f..%.1f)\n", st.size(), st.front().report, st.back().report, st.front().days, st.back().days);
        std::printf("  arrays at last step:"); for (auto& s : rc.restart.arrayNames(int(st.size()) - 1)) std::printf(" %s", s.c_str()); std::printf("\n");
        auto p = rc.restart.load(int(st.size()) - 1, "PRESSURE");
        if (!p.empty()) { double s = 0; for (float x : p) s += x; std::printf("  mean PRESSURE at last step %.3f\n", s / p.size()); }
    }
    if (rc.hasSummary) {
        const auto& S = rc.summary;
        std::printf("  summary: %zu vectors, %zu ministeps, %zu report steps, end %.1f d, start %02d-%02d-%04d, units %s\n", S.vars.size(), S.time.size(),
                    S.reportIndex.size(), S.time.empty() ? 0.0 : S.time.back(), S.startDay, S.startMonth, S.startYear, S.unitSystem.c_str());
        for (const char* k : {"FOPR", "FWCT", "FGOR", "FPR", "FOPT"}) {
            const int v = S.find(k);
            if (v >= 0 && !S.time.empty()) std::printf("    %-5s end = %.6g %s\n", k, S.values[size_t(v)].back(), S.vars[size_t(v)].unit.c_str());
        }
    }
    for (auto& w : rc.warnings) std::printf("  warning: %s\n", w.c_str());
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 3 && !std::strcmp(argv[1], "-dump")) return dump(argv[2]);
    if (argc >= 3 && !std::strcmp(argv[1], "-info")) return info(argv[2]);
    if (argc < 3) {
        std::fprintf(stderr, "usage: compare_opm <deck.DATA> <opm result file> [-t sec] [-csv out.csv] [-q]\n       compare_opm -dump <file> | -info <result file>\n");
        return 2;
    }
    double tmax = 600; std::string csv; bool quiet = false; std::vector<std::string> series;
    for (int i = 3; i < argc; ++i) {
        if (!std::strcmp(argv[i], "-t") && i + 1 < argc) tmax = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-csv") && i + 1 < argc) csv = argv[++i];
        else if (!std::strcmp(argv[i], "-q")) quiet = true;
        else if (!std::strcmp(argv[i], "-series") && i + 1 < argc) series.push_back(argv[++i]);
    }
    ecl::ResultCase rc; std::string err;
    if (!ecl::openResultCase(argv[2], rc, &err) || !rc.hasSummary) { std::fprintf(stderr, "OPM case: %s\n", err.empty() ? "no summary" : err.c_str()); return 1; }

    Case c;
    try { c = importEclipseDeck(argv[1]); }
    catch (const std::exception& e) { std::fprintf(stderr, "import: %s\n", e.what()); return 1; }
    Simulator s;
    s.log = [&](const std::string& m) { if (!quiet) std::fprintf(stderr, "  %s\n", m.c_str()); };
    const auto t0 = std::chrono::steady_clock::now();
    try {
        s.setup(c); s.initialize();
        while (s.runReportStep()) {
            if (std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() > tmax) { std::fprintf(stderr, "  time limit reached\n"); break; }
        }
    } catch (const std::exception& e) { std::fprintf(stderr, "run: %s\n", e.what()); }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::vector<std::string> names; for (auto& w : c.wells) names.push_back(w.name);
    const ecl::Summary mine = ecl::summaryFromRun(s.summary(), names, c.units, ecl::componentNames(c));
    const auto diffs = ecl::compareSummaries(rc.summary, mine);

    std::printf("deck %s: ResSim %.1f d in %.2f s (%d steps, %d Newton) vs OPM %.1f d\n", argv[1], mine.time.empty() ? 0.0 : mine.time.back(), secs,
                s.stats().timesteps, s.stats().newton, rc.summary.time.empty() ? 0.0 : rc.summary.time.back());
    std::printf("  %-16s %14s %14s %10s %8s\n", "vector", "OPM final", "ResSim final", "RMS", "relRMS%");
    double worstField = 0;
    for (const auto& d : diffs) {
        const bool fieldVec = d.key.find(':') == std::string::npos;
        if (fieldVec && d.key != "FOPT" && d.key != "FWPT" && d.key != "FGPT" && d.key != "FWIT" && d.key != "FGIT" && d.key[1] != 'C') worstField = std::max(worstField, d.relRms);
        if (!fieldVec && quiet) continue;
        std::printf("  %-16s %14.6g %14.6g %10.4g %8.2f\n", d.key.c_str(), d.refFinal, d.cmpFinal, d.rms, 100 * d.relRms);
    }
    std::printf("  worst field-rate/pressure relRMS: %.2f %%\n", 100 * worstField);
    for (const auto& key : series) {
        if (key.size() > 2 && key.substr(key.size() - 2) == ":*") {          // all wells, at the first report step
            const std::string kw = key.substr(0, key.size() - 2);
            const auto& R = rc.summary;
            if (R.reportIndex.empty()) continue;
            const int i = R.reportIndex[0];
            std::printf("%s at day %.2f\n", kw.c_str(), R.time[size_t(i)]);
            for (size_t v = 0; v < R.vars.size(); ++v) {
                if (R.vars[v].keyword != kw || R.vars[v].key.empty()) continue;
                const int b = mine.find(R.vars[v].key);
                std::printf("  %-12s %14.6g %14.6g\n", R.vars[v].key.c_str(), R.values[v][size_t(i)], b < 0 ? NAN : mine.at(b, R.time[size_t(i)]));
            }
            continue;
        }
        const int a = rc.summary.find(key), b = mine.find(key);
        if (a < 0 || b < 0) { std::printf("series %s: missing (%s)\n", key.c_str(), a < 0 ? "OPM" : "ResSim"); continue; }
        std::printf("series %s\n  %10s %14s %14s\n", key.c_str(), "day", "OPM", "ResSim");
        const auto& R = rc.summary;
        const size_t n = R.reportIndex.size(), stride = std::max<size_t>(1, n / 25);
        for (size_t k = 0; k < n; k += stride) {
            const int i = R.reportIndex[k];
            std::printf("  %10.2f %14.6g %14.6g\n", R.time[size_t(i)], R.values[size_t(a)][size_t(i)], mine.at(b, R.time[size_t(i)]));
        }
    }
    if (!csv.empty()) {
        std::ofstream o(csv);
        o << "vector,opm_final,ressim_final,rms,rel_rms\n";
        for (const auto& d : diffs) o << d.key << ',' << d.refFinal << ',' << d.cmpFinal << ',' << d.rms << ',' << d.relRms << "\n";
    }
    return 0;
}
