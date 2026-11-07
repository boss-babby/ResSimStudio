// norne_resv: checks the WCONHIST RESV conversion against ECLIPSE: history surface rates converted at a pressure
// (field FPR) vs the reported FVPR. Usage: norne_resv DECK SMSPEC
#include "ressim/ecl.h"
#include "ressim/simulator.h"
#include "ressim/io.h"
#include <cstdio>
#include <cmath>
using namespace ressim; using namespace ressim::ecl;
int main(int c, char** v) {
  Summary es; std::string e; readSummary(v[2], es, &e);
  std::vector<std::string> warn; Case cs = importEclipseDeck(v[1], &warn);
  Simulator s; s.log = [](const std::string&) {}; s.setup(cs);
  const FluidModel& F = s.fluid(0);
  const int iO = es.find("FOPRH"), iW = es.find("FWPRH"), iG = es.find("FGPRH"), iP = es.find("FPR"), iV = es.find("FVPR"),
            iWi = es.find("FWIR"), iVi = es.find("FVIR");
  std::printf("%8s %9s %12s %12s %12s  %9s %9s\n", "t", "FPR", "FVPR ecl", "conv@FPR", "ratio", "FWIR", "FVIR");
  for (double t : {8.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0, 3000.0}) {
    const double qo = es.at(iO, t), qw = es.at(iW, t), qg = es.at(iG, t), p = es.at(iP, t);
    const double rs = qo > 0 ? std::min(F.rsSat(p), qg / qo) : 0.0, rv = qg > 0 ? std::min(F.rvSat(p), qo / qg) : 0.0;
    const double den = 1 - rs * rv, O = (qo - rv * qg) / den, G = std::max(0.0, (qg - rs * qo) / den);
    const double resv = O / F.boS(p, rs) + qw / F.bwS(p) + G / F.bgS(p);
    std::printf("%8.0f %9.2f %12.6g %12.6g %12.5f  %9.4g %9.4g  rs %.4g rv %.4g\n", t, p, es.at(iV, t), resv, resv / es.at(iV, t), es.at(iWi, t), es.at(iVi, t), rs, rv);
  }
}
