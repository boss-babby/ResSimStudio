// norne_inj OURS_WELLS.csv REF.SMSPEC [t1 t2 ...] : injector / producer comparison (rates, BHP) against an ECLIPSE summary
#include "ressim/ecl.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace ressim;
using namespace ressim::ecl;
struct Row { double t; int mode, open; double oil, water, gas, winj, ginj, bhp; };
int main(int argc, char** argv) {
    if (argc < 3) { std::printf("usage: norne_inj OURS_WELLS.csv REF.SMSPEC [times]\n"); return 2; }
    Summary es; std::string err;
    if (!readSummary(argv[2], es, &err)) { std::printf("error: %s\n", err.c_str()); return 2; }
    std::vector<double> times;
    for (int i = 3; i < argc; ++i) times.push_back(std::atof(argv[i]));
    if (times.empty()) times = {300, 600, 1000, 1500, 2000, 2500, 3000, 3300};
    std::map<std::string, std::vector<Row>> w;
    std::ifstream f(argv[1]); std::string line; std::getline(f, line);
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string tok; std::vector<std::string> v;
        while (std::getline(ss, tok, ',')) v.push_back(tok);
        if (v.size() < 10) continue;
        Row r{std::atof(v[0].c_str()), std::atoi(v[2].c_str()), std::atoi(v[3].c_str()), std::atof(v[4].c_str()), std::atof(v[5].c_str()),
              std::atof(v[6].c_str()), std::atof(v[7].c_str()), std::atof(v[8].c_str()), std::atof(v[9].c_str())};
        w[v[1]].push_back(r);
    }
    auto ours = [&](const std::vector<Row>& rs, double t) -> const Row* {
        const Row* best = nullptr;
        for (const Row& r : rs) if (r.t <= t + 1e-6) best = &r;
        return best;
    };
    auto ecl = [&](const std::string& k, double t) { const int i = es.find(k); return i < 0 ? NAN : es.at(i, t); };
    for (auto& [name, rs] : w) {
        double totInj = 0; for (const Row& r : rs) totInj += r.winj + r.ginj;
        const bool inj = totInj > 0 || ecl("WWIT:" + name, 1e9) > 0 || ecl("WGIT:" + name, 1e9) > 0;
        std::printf("%s %s\n", name.c_str(), inj ? "(injector)" : "(producer)");
        // cumulative injection (rate x step length) for a timing-independent comparison
        std::vector<double> cumW(rs.size()), cumG(rs.size());
        for (size_t k = 0; k < rs.size(); ++k) {
            const double dt = rs[k].t - (k ? rs[k - 1].t : 0.0);
            cumW[k] = (k ? cumW[k - 1] : 0.0) + rs[k].winj * dt; cumG[k] = (k ? cumG[k - 1] : 0.0) + rs[k].ginj * dt;
        }
        for (double t0 : times) {
            const Row* r = ours(rs, t0);
            if (!r) continue;
            const double t = r->t - 1e-3;   // ECLIPSE ministep values are step averages ending at the report time
            const size_t k = size_t(r - rs.data());
            if (inj)
                std::printf("  t %6.0f  WWIT %10.4g / %10.4g  WGIT %10.4g / %10.4g\n", r->t, cumW[k], ecl("WWIT:" + name, r->t), cumG[k], ecl("WGIT:" + name, r->t));
            if (inj)
                std::printf("  t %6.0f  WWIR %9.1f / %9.1f  WGIR %10.0f / %10.0f  BHP %6.1f / %6.1f  mode %d open %d\n", t, r->winj, ecl("WWIR:" + name, t), r->ginj,
                            ecl("WGIR:" + name, t), r->bhp, ecl("WBHP:" + name, t), r->mode, r->open);
            else
                std::printf("  t %6.0f  WOPR %8.1f / %8.1f  WWPR %8.1f / %8.1f  WGPR %10.0f / %10.0f  BHP %6.1f / %6.1f  mode %d\n", t, r->oil, ecl("WOPR:" + name, t),
                            r->water, ecl("WWPR:" + name, t), r->gas, ecl("WGPR:" + name, t), r->bhp, ecl("WBHP:" + name, t), r->mode);
        }
    }
}
