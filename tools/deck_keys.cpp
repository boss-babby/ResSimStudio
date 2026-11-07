#include "ressim/ecl.h"
#include <cstdio>
#include <map>
using namespace ressim::ecl;
int main(int c, char** v) { Summary s; std::string e; readSummary(v[1], s, &e); std::map<std::string,int> k; for (auto& x : s.vars) k[x.keyword]++; for (auto& [a,b] : k) std::printf("%s(%d) ", a.c_str(), b); std::printf("\n"); }
