// norne_conn: compares per-connection producer rates / pressures with an ECLIPSE PRT production report.
// Usage: norne_conn DECK PRT T [WELL ...]   (T = report time in days; no WELL -> well-level table only)
#include "ressim/simulator.h"
#include "ressim/io.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace ressim;
struct Row { std::string name; int i = 0, j = 0, k = 0; std::string mode; double o = 0, w = 0, g = 0, v = 0, p1 = 0, p2 = 0; };
static std::string trim(std::string s) { while (!s.empty() && s.back() == ' ') s.pop_back(); size_t b = s.find_first_not_of(' '); return b == std::string::npos ? "" : s.substr(b); }
int main(int argc, char** argv) {
    if (argc < 4) { std::printf("usage: norne_conn DECK PRT T [WELL...]\n"); return 2; }
    const double T = std::atof(argv[3]);
    std::vector<std::string> wantW; for (int i = 4; i < argc; ++i) wantW.push_back(argv[i]);
    std::ifstream f(argv[2]); std::string line; bool inT = false, inProd = false;
    std::vector<Row> wells; std::map<std::string, std::vector<Row>> blocks; std::string cur;
    while (std::getline(f, line)) {
        size_t at = line.find("WELLS    AT");
        if (at != std::string::npos) { const double t = std::atof(line.c_str() + at + 11); if (inT && t != T) break; inT = std::fabs(t - T) < 1e-6; inProd = false; continue; }
        if (!inT) continue;
        if (line.find("PRODUCTION REPORT") != std::string::npos) { inProd = true; continue; }
        if (line.find("INJECTION REPORT") != std::string::npos || line.find("CUMULATIVE PRODUCTION") != std::string::npos) { inProd = false; continue; }
        if (!inProd || line.size() < 2 || line[0] != ':') continue;
        std::vector<std::string> c; { std::stringstream ss(line); std::string x; while (std::getline(ss, x, ':')) c.push_back(trim(x)); }
        if (c.size() < 13 || c[1].empty() || c[1] == "WELL" || c[1] == "NAME" || c[1] == "OR GRID" || c[1].find("FIELD") == 0 || c[1].find("GRP") != std::string::npos || c[2].find("GRP") != std::string::npos) continue;
        Row r; r.name = c[1]; r.mode = c[3];
        std::sscanf(c[2].c_str(), "%d, %d, %d", &r.i, &r.j, &r.k);
        r.o = std::atof(c[4].c_str()); r.w = std::atof(c[5].c_str()); r.g = std::atof(c[6].c_str()); r.v = std::atof(c[7].c_str());
        r.p1 = std::atof(c[11].c_str()); r.p2 = std::atof(c[12].c_str());
        if (r.name == "BLOCK") { if (!cur.empty()) blocks[cur].push_back(r); }
        else { cur = r.name; wells.push_back(r); }
    }
    std::printf("ECLIPSE report at %.1f d: %zu producers\n", T, wells.size());
    std::vector<std::string> warn; Case cs = importEclipseDeck(argv[1], &warn);
    Simulator s; s.log = [](const std::string&) {}; s.setup(cs); s.initialize();
    while (s.summary().back().time < T - 1e-6) if (!s.runReportStep()) break;
    std::printf("ours at %.1f d\n", s.summary().back().time);
    const auto cr = s.connectionRates();
    std::map<std::string, std::vector<Simulator::ConnRate>> ours;
    for (auto& c : cr) ours[c.well].push_back(c);
    const auto& names = s.wellNames(); const auto& sp = s.summary().back();
    std::printf("%-8s %5s %21s %21s %23s %17s\n", "well", "mode", "oil ours/ecl", "water ours/ecl", "gas ours/ecl", "BHP ours/ecl");
    for (auto& w : wells) {
        int k = -1; for (size_t q = 0; q < names.size(); ++q) if (names[q] == w.name) k = int(q);
        if (k < 0 || size_t(k) >= sp.wells.size()) continue;
        const auto& r = sp.wells[size_t(k)];
        std::printf("%-8s %5s %10.1f/%-10.1f %10.1f/%-10.1f %11.4g/%-11.4g %8.1f/%-8.1f\n", w.name.c_str(), w.mode.c_str(), r.oil, w.o, r.water, w.w, r.gas, w.g, r.bhp, w.p1);
    }
    for (auto& wn : wantW) {
        std::printf("\n%s connections (I J K): oil, water, gas ours/ecl | block P ours/ecl | conn P ours/ecl | ours Sw Sg\n", wn.c_str());
        for (auto& b : blocks[wn]) {
            const Simulator::ConnRate* m = nullptr;
            for (auto& c : ours[wn]) if (c.i + 1 == b.i && c.j + 1 == b.j && c.k + 1 == b.k) m = &c;
            if (!m) { std::printf("  %3d %3d %3d %-4s  ecl %.1f %.1f %.4g  (no open connection in ours)\n", b.i, b.j, b.k, b.mode.c_str(), b.o, b.w, b.g); continue; }
            std::printf("  %3d %3d %3d %-4s %8.1f/%-8.1f %7.2f/%-7.2f %10.4g/%-10.4g | %6.1f/%-6.1f | %6.1f/%-6.1f | %.3f %.3f\n", b.i, b.j, b.k, b.mode.c_str(),
                        m->oil, b.o, m->water, b.w, m->gas, b.g, m->pBlock, b.p2, m->pConn, b.p1, m->sw, m->sg);
        }
    }
}
