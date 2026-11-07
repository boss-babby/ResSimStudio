// Developer tool (C47): regression metric between two reference summaries (e.g. two OPM versions of the same opm-tests case),
// exactly as `ResSimCLI --regress` scores a run: relative RMS deviation (RMS / max|ref|) per vector + the C41 material worst.
//   refdrift NEW_REF.SMSPEC OLD_REF.SMSPEC [tol]
// Build: g++ -std=c++17 -O1 -fopenmp -Iengine/include tools/dev/refdrift.cpp build/linux/libressim_engine.a -o refdrift -lpthread -ldl
#include "ressim/ecl.h"
#include "ressim/regress.h"
#include <cstdio>
#include <cstdlib>
using namespace ressim;
int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: refdrift NEW_REF.SMSPEC OLD_REF.SMSPEC [tol]\n"); return 1; }
    regress::CaseResult r; std::string err;
    if (!ecl::readSummary(argv[1], r.ref, &err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 1; }
    if (!ecl::readSummary(argv[2], r.run, &err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 1; }
    regress::Manifest m;
    regress::compare(r, m.keys, argc > 3 ? std::atof(argv[3]) : 0.05);
    std::printf("worst %s %.3f %%  %s\n", r.worstKey.c_str(), 100 * r.worst, r.message.c_str());
    return 0;
}
