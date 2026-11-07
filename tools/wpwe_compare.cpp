// C44: WPWE event timelines of a ResSim run of a deck against an OPM Flow reference summary.
//   wpwe_compare DECK.DATA REF.SMSPEC      -> one line per (vector, well): reference events / ResSim events, MATCH / DIFF
// Event times are compared at the report steps of the reference (events are summed over the ministeps of each report step,
// because the two simulators take different time steps).
#include "ressim/api.h"
#include "ressim/ecl.h"
#include <cstdio>
#include <map>
#include <set>
#include <string>
using namespace ressim;
static std::map<std::string, std::map<int, double>> events(const ecl::Summary& S, const std::vector<double>& reportT) {
    std::map<std::string, std::map<int, double>> out;
    for (size_t v = 0; v < S.vars.size(); ++v) {
        const auto& V = S.vars[v];
        if (V.keyword.rfind("WPWE", 0) != 0 || V.keyword == "WPWEM") continue;
        for (size_t m = 0; m < S.time.size() && m < S.values[v].size(); ++m) {
            const double x = S.values[v][m];
            if (x == 0) continue;
            int r = 0; while (r < int(reportT.size()) && S.time[m] > reportT[size_t(r)] + 1e-6) ++r;
            out[V.keyword + ":" + V.wgname][r] += x;
        }
    }
    return out;
}
int main(int argc, char** argv) {
    if (argc < 3) { std::printf("usage: wpwe_compare DECK.DATA REF.SMSPEC\n"); return 1; }
    ecl::Summary R; std::string err;
    if (!ecl::readSummary(argv[2], R, &err)) { std::printf("ERROR %s\n", err.c_str()); return 2; }
    api::Session S;
    if (!S.open(argv[1]) || !S.run()) { std::printf("ERROR %s\n", S.error().c_str()); return 3; }
    const ecl::Summary& O = S.summary();
    std::vector<double> rt; for (int i : R.reportIndex) if (i >= 0 && size_t(i) < R.time.size()) rt.push_back(R.time[size_t(i)]);
    const auto er = events(R, rt), eo = events(O, rt);
    std::set<std::string> keys; for (auto& [k, v] : er) keys.insert(k); for (auto& [k, v] : eo) keys.insert(k);
    int diff = 0;
    for (const auto& k : keys) {
        auto fmtE = [&](const std::map<std::string, std::map<int, double>>& E) {
            std::string s; auto it = E.find(k); if (it == E.end()) return std::string("-");
            for (auto& [r, x] : it->second) { char b[64]; std::snprintf(b, sizeof b, "%s%.4gd=%g", s.empty() ? "" : " ", r < int(rt.size()) ? rt[size_t(r)] : -1.0, x); s += b; }
            return s;
        };
        const std::string a = fmtE(er), b = fmtE(eo);
        const bool same = a == b; diff += !same;
        std::printf("%-6s %-16s ref [%s]  ressim [%s]\n", same ? "MATCH" : "DIFF", k.c_str(), a.c_str(), b.c_str());
    }
    std::printf("%s: %zu event series, %d differ\n", argv[1], keys.size(), diff);
    return diff ? 4 : 0;
}
