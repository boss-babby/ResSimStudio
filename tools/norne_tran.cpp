#include "ressim/ecl.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include <cstdio>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>
using namespace ressim; using namespace ressim::ecl;
int main(int c, char** v) {
  ResultCase rc; std::string e; openResultCase(v[2], rc, &e);
  std::vector<std::string> warn; Case cs = importEclipseDeck(v[1], &warn);
  Simulator sim; sim.log = [](const std::string&) {}; sim.setup(cs);
  const Grid& G = sim.grid(); const int nx = G.nx, ny = G.ny;
  auto TX = rc.initArray("TRANX"), TY = rc.initArray("TRANY");
  std::vector<int> fl; { std::ifstream f(v[4]); std::string tok; bool on = false;
    while (f >> tok) { if (tok.rfind("--", 0) == 0) { std::string rest; std::getline(f, rest); continue; } if (tok == "FLUXNUM") { on = true; continue; } if (!on) continue; if (tok == "/") break;
      size_t st = tok.find('*'); if (st == std::string::npos) fl.push_back(std::atoi(tok.c_str())); else { int n = std::atoi(tok.substr(0, st).c_str()), x = std::atoi(tok.substr(st + 1).c_str()); for (int i = 0; i < n; ++i) fl.push_back(x); } } }
  std::printf("fluxnum size %zu\n", fl.size());
  { // NNCs across FLUXNUM boundaries
    auto tn = rc.init.numeric("TRANNNC"); File eg; eg.open(v[3]); auto n1 = eg.ints("NNC1"), n2 = eg.ints("NNC2");
    std::map<std::pair<int,int>, double> ours;
    for (const auto& k : G.conns) { int ga = G.globalOf[size_t(k.a)], gb = G.globalOf[size_t(k.b)]; if (ga > gb) std::swap(ga, gb); ours[{ga, gb}] += k.trans; }
    std::map<std::pair<int,int>, double> mr; { std::ifstream f(v[5]); std::string line; bool on = false;
      while (std::getline(f, line)) { if (line.rfind("--", 0) == 0) continue; if (line.rfind("MULTREGT", 0) == 0) { on = true; continue; } if (!on) continue;
        std::stringstream ss(line); int a, b; double m; if (ss >> a >> b >> m) { mr[{std::min(a,b), std::max(a,b)}] = m; } } }
    std::printf("MULTREGT pairs %zu\n", mr.size());
    double sE0 = 0, sO0 = 0, sEm = 0, sOm = 0; int n0 = 0, nm = 0;
    for (size_t i = 0; i < tn.size(); ++i) { int a = n1[i] - 1, b = n2[i] - 1; if (a > b) std::swap(a, b); int ra = fl[size_t(a)], rb = fl[size_t(b)]; if (ra == rb) continue;
      auto im = mr.find({std::min(ra,rb), std::max(ra,rb)}); double m = im == mr.end() ? 1.0 : im->second;
      auto it = ours.find({a, b}); double to = it != ours.end() ? it->second : 0;
      if (m == 0) { ++n0; sE0 += tn[i]; sO0 += to; } else if (m < 1) { ++nm; sEm += tn[i]; sOm += to; } }
    std::printf("NNCs where MULTREGT = 0: %d, ecl sum %.6g, ours %.6g ; 0<m<1: %d ecl %.6g ours %.6g\n", n0, sE0, sO0, nm, sEm, sOm);
    // neighbour connections across boundaries with m = 0
    double nE = 0, nO = 0; int nn = 0;
    for (const auto& k : G.conns) { int ga = G.globalOf[size_t(k.a)], gb = G.globalOf[size_t(k.b)]; if (ga > gb) std::swap(ga, gb); int ra = fl[size_t(ga)], rb = fl[size_t(gb)]; if (ra == rb) continue;
      auto im = mr.find({std::min(ra,rb), std::max(ra,rb)}); if (im == mr.end() || im->second != 0) continue; ++nn; nO += k.trans; }
    std::printf("ALL our connections across m=0 boundaries: %d, sum %.6g\n", nn, nO); (void)nE;
    int nb = 0; double sO = 0, sE = 0; int nz = 0;
    for (size_t i = 0; i < tn.size(); ++i) { int a = n1[i] - 1, b = n2[i] - 1; if (a > b) std::swap(a, b); if (fl[size_t(a)] == fl[size_t(b)]) continue;
      ++nb; sE += tn[i]; auto it = ours.find({a, b}); if (it != ours.end()) sO += it->second; if (tn[i] < 1e-6) ++nz; }
    std::printf("ECL NNCs across FLUXNUM boundaries: %d, sum ecl %.6g ours %.6g, ecl ~zero %d\n", nb, sE, sO, nz);
  }
  const char* bins[] = {"<0.01", "0.01-0.5", "0.5-0.9", "0.9-0.99", "0.99-1.01", "1.01-1.1", "1.1-2", "2-100", ">100"};
  const double edges[] = {0.01, 0.5, 0.9, 0.99, 1.01, 1.1, 2, 100};
  for (int d = 0; d < 2; ++d) {
    int cnt[9] = {0}, bnd[9] = {0}; double so[9] = {0}, se[9] = {0};
    for (const auto& k : G.conns) {
      int ga = G.globalOf[size_t(k.a)], gb = G.globalOf[size_t(k.b)]; if (ga > gb) std::swap(ga, gb);
      const int step = d == 0 ? 1 : nx; if (gb - ga != step || k.dir != "XY"[d]) continue;
      const double te = (d == 0 ? TX : TY)[size_t(G.activeOf[size_t(ga)])];
      const double r = te > 0 ? k.trans / te : 1e9;
      int b = 0; while (b < 8 && r >= edges[b]) ++b;
      cnt[b]++; so[b] += k.trans; se[b] += te; if (!fl.empty() && fl[size_t(ga)] != fl[size_t(gb)]) bnd[b]++;
    }
    std::printf("%c ratio ours/ecl:\n", "XY"[d]);
    for (int b = 0; b < 9; ++b) if (cnt[b]) std::printf("  %-10s n %6d (fluxnum boundary %5d)  sum ours %.4g ecl %.4g\n", bins[b], cnt[b], bnd[b], so[b], se[b]);
  }
}
