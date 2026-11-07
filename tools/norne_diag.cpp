#include "ressim/ecl.h"
#include "ressim/simulator.h"
#include "ressim/io.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
using namespace ressim; using namespace ressim::ecl;
int main(int c, char** v) {
  Summary es; std::string e; readSummary(v[2], es, &e);
  std::vector<double> times; { std::string a = v[3]; size_t p0 = 0; while (p0 < a.size()) { size_t q = a.find(',', p0); if (q == std::string::npos) q = a.size(); times.push_back(std::atof(a.substr(p0, q - p0).c_str())); p0 = q + 1; } }
  std::vector<std::string> warn; Case cs = importEclipseDeck(v[1], &warn);
  Simulator s; s.log = [](const std::string&) {}; s.setup(cs); s.initialize();
  for (double tEnd : times) {
  while (s.summary().empty() || s.summary().back().time < tEnd - 1e-6) if (!s.runReportStep()) break;
  const double t = s.summary().back().time;
  const Snapshot& sn = s.snapshots().back();
  const Grid& G = s.grid();
  std::map<int, std::pair<double,double>> reg; std::map<int, double> roip;
  { const FluidModel& F = s.fluid(0);
    for (int a = 0; a < G.numActive(); ++a) { int g = G.globalOf[size_t(a)]; int r = cs.grid.fipnum.empty() ? 1 : cs.grid.fipnum[size_t(g)];
      const double p = sn.pressure[size_t(a)], so = sn.so[size_t(a)], sg = sn.sg[size_t(a)], rs = sn.rs.empty() ? 0.0 : sn.rs[size_t(a)], rv = sn.rv.empty() ? 0.0 : sn.rv[size_t(a)];
      roip[r] += G.poreVolRef[size_t(a)] * (so * F.boS(p, rs) + sg * F.bg(AD3(p), AD3(rv)).v * rv); } }
  for (int a = 0; a < G.numActive(); ++a) { int g = G.globalOf[size_t(a)]; int r = cs.grid.fipnum.empty() ? 1 : cs.grid.fipnum[size_t(g)];
    double w = G.poreVolRef[size_t(a)] * (1 - sn.sw[size_t(a)]); reg[r].first += w * sn.pressure[size_t(a)]; reg[r].second += w; }
  std::printf("t = %.1f (snapshot %.1f)\nregion RPR ours/ecl\n", t, sn.time);
  for (auto& [r, pw] : reg) { int id = es.find("RPR:" + std::to_string(r)); if (pw.second > 0) { int io = es.find("ROIP:" + std::to_string(r)); std::printf("  %2d  %7.2f / %7.2f   ROIP %.5g / %.5g\n", r, pw.first / pw.second, id < 0 ? NAN : es.at(id, t), roip[r], io < 0 ? NAN : es.at(io, t)); } }
  auto conns = s.connections();
  const auto& p = s.summary().back(); const auto& names = s.wellNames();
  std::printf("%-8s %17s %17s %8s  %s\n", "well", "BHP ours/ecl", "WBP ours/ecl", "ECL WPI", "rates ours (o w g / wi gi)");
  for (size_t k = 0; k < names.size(); ++k) {
    if (k >= p.wells.size() || !p.wells[k].open) continue;
    double sp = 0; int n = 0;
    for (auto& cn : conns) if (cn.well == names[k]) { int g = cn.i + G.nx * (cn.j + G.ny * cn.k); int a = G.activeOf[size_t(g)]; if (a >= 0) { sp += sn.pressure[size_t(a)]; ++n; } }
    auto ev = [&](const char* kw) { int id = es.find(std::string(kw) + ":" + names[k]); return id < 0 ? NAN : es.at(id, t); };
    const auto& w = p.wells[k];
    std::printf("%-8s %7.1f/%-7.1f %7.1f/%-7.1f %8.3g  %.4g %.4g %.4g / %.4g %.4g  (ecl %.4g %.4g %.4g / %.4g %.4g)\n", names[k].c_str(), w.bhp, ev("WBHP"), n ? sp / n : NAN, ev("WBP"), ev("WPI"),
      w.oil, w.water, w.gas, w.waterInj, w.gasInj, ev("WOPR"), ev("WWPR"), ev("WGPR"), ev("WWIR"), ev("WGIR"));
  }
  std::fflush(stdout);
  }
}
