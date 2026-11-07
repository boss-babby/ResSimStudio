// compare_csv: compares a rundeck -csv field summary with an ECLIPSE / OPM summary (SMSPEC/UNSMRY).
// Usage: compare_csv OURS.csv REF.SMSPEC [every_n_rows] [VEC1,VEC2,...]
#include "ressim/ecl.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace ressim::ecl;
int main(int argc, char** argv) {
    if (argc < 3) { std::printf("usage: compare_csv OURS.csv REF.SMSPEC [N]\n"); return 2; }
    Summary es; std::string err;
    if (!readSummary(argv[2], es, &err)) { std::printf("error: %s\n", err.c_str()); return 2; }
    const int every = argc > 3 ? std::atoi(argv[3]) : 10;
    std::ifstream f(argv[1]); std::string line; std::getline(f, line);
    std::vector<std::string> hdr; { std::stringstream ss(line); std::string h; while (std::getline(ss, h, ',')) hdr.push_back(h); }
    std::vector<std::string> show = {"FPR", "FOPT", "FGPT", "FWPT", "FOPR", "FWPR", "FGPR"};
    if (argc > 4) { show.clear(); std::stringstream ss(argv[4]); std::string v; while (std::getline(ss, v, ',')) show.push_back(v); }
    std::printf("%8s", "time");
    for (auto s : show) std::printf(" %23s", (s + " ours/ref").c_str());
    std::printf("\n");
    std::vector<std::vector<double>> rows;
    while (std::getline(f, line)) { std::vector<double> r; std::stringstream ss(line); std::string v; while (std::getline(ss, v, ',')) r.push_back(std::atof(v.c_str())); rows.push_back(r); }
    for (size_t i = 0; i < rows.size(); ++i) {
        if (i % size_t(every) != 0 && i + 1 != rows.size()) continue;
        const double t = rows[i][0];
        std::printf("%8.1f", t);
        for (auto s : show) {
            int col = -1; for (size_t c = 0; c < hdr.size(); ++c) if (hdr[c] == s) col = int(c);
            const int id = es.find(s);
            const double a = col < 0 ? NAN : rows[i][size_t(col)], b = id < 0 ? NAN : es.at(id, t);
            std::printf(" %11.5g/%-11.5g", a, b);
        }
        std::printf("\n");
    }
    return 0;
}
