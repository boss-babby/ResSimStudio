// compare_smry: compares two ECLIPSE / OPM summary results (e.g. an original deck and its ResSim re-export, both
// run by OPM Flow). Prints the worst vectors by relative RMS. Usage: compare_smry REF.SMSPEC CMP.SMSPEC [N]
#include "ressim/ecl.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
using namespace ressim::ecl;
int main(int argc, char** argv) {
    if (argc < 3) { std::printf("usage: compare_smry REF CMP [N]\n"); return 2; }
    ResultCase a, b; std::string err;
    if (!openResultCase(argv[1], a, &err) || !openResultCase(argv[2], b, &err)) { std::printf("error: %s\n", err.c_str()); return 2; }
    auto d = compareSummaries(a.summary, b.summary);
    std::sort(d.begin(), d.end(), [](const VectorDiff& x, const VectorDiff& y) { return x.relRms > y.relRms; });
    const int n = argc > 3 ? std::atoi(argv[3]) : 8;
    std::printf("%zu common vectors; worst relRMS %.4g %%\n", d.size(), d.empty() ? 0.0 : 100 * d[0].relRms);
    for (int i = 0; i < n && i < int(d.size()); ++i)
        std::printf("  %-22s ref %-12.6g cmp %-12.6g relRMS %.4g %%\n", d[size_t(i)].key.c_str(), d[size_t(i)].refFinal, d[size_t(i)].cmpFinal, 100 * d[size_t(i)].relRms);
    return 0;
}
