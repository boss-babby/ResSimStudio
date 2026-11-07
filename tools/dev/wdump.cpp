#include "ressim/io.h"
#include <cstdio>
using namespace ressim;
int main(int argc, char** argv) {
    std::vector<std::string> w; Case c = importEclipseDeck(argv[1], &w);
    for (const auto& W : c.wells) if (W.name == argv[2]) {
        for (const auto& q : W.completions) std::printf("i%d j%d k%d-%d open %d from %d until %d def %d cf %g\n", q.i+1, q.j+1, q.k1+1, q.k2+1, int(q.open), q.fromStep, q.untilStep, q.defStep, q.cf);
    }
}
