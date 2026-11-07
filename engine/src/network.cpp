// ResSim - C29 surface production network of the black-oil engine (ECLIPSE / OPM "extended network")
//
// Input: the ECLIPSE keywords NODEPROP (node, fixed pressure of the terminal node, ...) and BRANPROP (downtree node,
// uptree node, VFPPROD table of the branch, ALQ) stay verbatim extra keywords (exact round trip, OPM Flow runs them too).
// Nodes are groups; a well belongs to the first network node met walking up its group chain (WELSPECS group, GRUPTREE).
//
// Balance (C29 explicit start-of-step estimate, C30 implicit: after the step converges the node pressures are recomputed
// from the new rates and the step is re-solved until they agree - Illinois / regula-falsi per node, see runReportStep): node rates = sum of the latest converged surface rates of the node's wells and of its downtree nodes; node
// pressures top-down from the fixed-pressure terminal node, p_down = VFP_branch(q_down, p_up, WCT / WGR, GOR / OGR /
// GLR, ALQ) (VFP 0: no pressure loss, 9999: branch removed). Every producer of a node gets the node pressure as its THP
// limit; with a well VFP table (WCONPROD item 13) updateThp converts it to the BHP limit at the current inflow (bisection
// on IPR = VFP), without one the BHP limit is the node pressure plus the mixture head to the well's datum.
#include <climits>
#include "ressim/simulator.h"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <map>
#include <functional>

namespace ressim {
namespace {
std::string fmt(const char* f, double a = 0, double b = 0, double c = 0, double d = 0) {
    char buf[512]; std::snprintf(buf, sizeof buf, f, a, b, c, d); return buf;
}
std::string nUpper(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
bool nDefault(const std::vector<std::string>& r, size_t i) { return i >= r.size() || r[i].empty() || r[i] == "1*" || r[i] == "*"; }
}  // namespace

namespace {
bool nMatch(const std::string& pat, const std::string& name) {   // '*' / '?' wildcards, upper case
    const std::string a = nUpper(pat), b = nUpper(name);
    size_t i = 0, j = 0, star = std::string::npos, mark = 0;
    while (j < b.size()) {
        if (i < a.size() && (a[i] == '?' || a[i] == b[j])) { ++i; ++j; }
        else if (i < a.size() && a[i] == '*') { star = i++; mark = j; }
        else if (star != std::string::npos) { i = star + 1; j = ++mark; }
        else return false;
    }
    while (i < a.size() && a[i] == '*') ++i;
    return i == a.size();
}
std::string nUnq(std::string s) { if (s.size() >= 2 && s.front() == '\'') s = s.substr(1, s.size() - 2); return s; }
bool nYes(const std::vector<std::string>& r, size_t i, bool def) {
    if (nDefault(r, i)) return def;
    const std::string v = nUpper(nUnq(r[i]));
    return v.rfind("Y", 0) == 0 || v == "FLO";
}
int nAlqEq(const std::vector<std::string>& r, size_t i) {
    if (nDefault(r, i)) return 0;
    const std::string v = nUpper(nUnq(r[i]));
    return v == "OIL" ? 1 : v == "GAS" ? 2 : 0;
}
}  // namespace

bool Simulator::groupBelow(const std::string& g, const std::string& top) const {
    std::string x = nUpper(g);
    const std::string t = nUpper(top);
    for (int guard = 0; guard < 64 && !x.empty(); ++guard) {
        if (x == t) return true;
        if (x == "FIELD") return false;
        std::string up = "FIELD";
        for (const auto& [c, p] : gParent_) if (nUpper(c) == x) { up = nUpper(p); break; }
        x = up;
    }
    return false;
}

double Simulator::satRate(const std::string& top, int what) const {
    double q = 0;
    for (const auto& [g, S] : gsat_) {
        if (!groupBelow(g, top)) continue;
        q += what == 0 ? S.o : what == 1 ? S.w : what == 2 ? S.g : what == 3 ? S.resv : S.lift;
    }
    return q;
}

// C35 GSATPROD (satellite group production: oil, water, gas, RESV, lift gas - added to the group / field rates, fixed in
// group control) and GCONSUMP (gas consumption and import of a group; item 4 = network node where it leaves the flow)
void Simulator::facilityControls(int step) {
    if (step <= 0) { gsat_.clear(); gcons_.clear(); gsale_.clear(); }
    if (step < 0) return;
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.section != DeckSection::Schedule || k.step != step) continue;
        if (k.name == "GSATPROD") {
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.empty() || r[0].empty()) continue;
                SatProd S;
                S.o = verbatimNum(r, 1, 0); S.w = verbatimNum(r, 2, 0); S.g = verbatimNum(r, 3, 0); S.resv = verbatimNum(r, 4, 0); S.lift = verbatimNum(r, 5, 0);
                gsat_[nUpper(nUnq(r[0]))] = S;
                log_("  satellite group " + nUpper(nUnq(r[0])) + fmt(": oil %g, water %g, gas %g", S.o, S.w, S.g) + fmt(", RESV %g, lift gas %g", S.resv, S.lift));
            }
        } else if (k.name == "GCONSALE") {
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.empty() || r[0].empty()) continue;
                GasSale S;
                S.target = verbatimNum(r, 1, -1); S.maxRate = verbatimNum(r, 2, -1); S.minRate = verbatimNum(r, 3, -1);
                S.proc = nDefault(r, 4) ? std::string("NONE") : nUpper(nUnq(r[4]));
                const std::string g = nUpper(nUnq(r[0]));
                gsale_[g] = S;
                log_("  gas sales of group " + g + fmt(": target %g, maximum %g, minimum %g", S.target, S.maxRate, S.minRate) + " (" + S.proc + ")");
            }
        } else if (k.name == "GCONSUMP") {
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.empty() || r[0].empty()) continue;
                GasConsump G;
                G.cons = verbatimNum(r, 1, 0); G.imp = verbatimNum(r, 2, 0);
                G.node = nDefault(r, 3) ? std::string() : nUpper(nUnq(r[3]));
                const std::string g = nUpper(nUnq(r[0]));
                if (G.cons == 0 && G.imp == 0) gcons_.erase(g); else gcons_[g] = G;
                log_("  gas consumption of group " + g + fmt(": %g, import %g", G.cons, G.imp) + (G.node.empty() ? std::string() : " at network node " + G.node));
            }
        }
    }
}

// gas available below a group for re-injection / sales: production (efficiency-weighted, + satellite) - consumption + import
double Simulator::groupGasAvail(const std::string& g) const {
    double q = satRate(g, 2);
    for (const WellRT& W : wells_)
        if (W.open && W.type == WellType::Producer && groupBelow(case_.wells[size_t(W.input)].group, g)) q += W.eff * W.rates.gas;
    for (const auto& [n, C] : gcons_) if (groupBelow(n, g)) q -= C.cons - C.imp;
    return q;
}

void Simulator::networkSetup(int step) {
    if (netStep_ == step) return;
    netStep_ = step;
    std::map<std::string, NetNode> nodes;
    std::vector<std::string> order;
    auto node = [&](const std::string& n) -> NetNode& {
        auto it = nodes.find(n);
        if (it == nodes.end()) { NetNode N; N.name = n; it = nodes.emplace(n, N).first; order.push_back(n); }
        return it->second;
    };
    // groups known to the case (GRUPNET name patterns, GRUPTREE parents)
    std::vector<std::string> groups = {"FIELD"};
    auto addG = [&](const std::string& g) { const std::string u = nUpper(g); if (!u.empty() && std::find(groups.begin(), groups.end(), u) == groups.end()) groups.push_back(u); };
    for (const auto& [c, p] : case_.groupTree) { addG(c); addG(p); }
    for (const auto& w : case_.wells) addG(w.group);
    // GRUPTREE as of the keyword being processed (GRUPNET takes the parent at that time; a later GRUPTREE move of a group
    // re-routes its standard-network branch, VFP / ALQ kept - OPM Schedule::addGroupToGroup)
    std::map<std::string, std::string> gpar;
    bool treeKw = false;
    for (const ExtraKeyword& k : case_.extraKeywords) if (k.name == "GRUPTREE" && k.section == DeckSection::Schedule) { treeKw = true; break; }
    if (!treeKw) for (const auto& [c, p] : case_.groupTree) gpar[nUpper(c)] = nUpper(p);
    auto parentOf = [&](const std::string& g) -> std::string {
        if (g == "FIELD") return std::string();
        auto it = gpar.find(g);
        return it == gpar.end() ? std::string("FIELD") : it->second;
    };
    bool any = false;
    netStd_ = false;
    bool compOff = false;
    netBal_ = NetBalan();
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.name == "NETBALAN" && k.step <= step) {   // C36 (SPECIAL or SCHEDULE; the latest in force)
            const auto R = verbatimRecords(k.body);
            netBal_ = NetBalan(); netBal_.given = true;
            if (!R.empty()) {
                const auto& r = R[0];
                const double bar = case_.units == UnitSystem::Field ? 14.503773773 : 1.0;
                netBal_.interval = nDefault(r, 0) ? 0.0 : std::atof(r[0].c_str());
                netBal_.ptol = nDefault(r, 1) ? 0.1 * bar : std::atof(r[1].c_str());
                netBal_.maxIter = nDefault(r, 2) ? 10 : std::max(0, std::atoi(r[2].c_str()));
                netBal_.targetErr = nDefault(r, 5) ? -1.0 : std::atof(r[5].c_str());
                netBal_.maxErr = nDefault(r, 6) ? -1.0 : std::atof(r[6].c_str());
            }
            continue;
        }
        if (k.section != DeckSection::Schedule || k.step > step) continue;
        if (k.name == "COMPOFF") { if (k.step == step) compOff = true; continue; }
        if (k.name == "GRUPTREE") {
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.size() < 2 || r[0].empty()) continue;
                const std::string c = nUpper(nUnq(r[0])), up = nDefault(r, 1) ? std::string("FIELD") : nUpper(nUnq(r[1]));
                gpar[c] = up; addG(c); addG(up);
                if (netStd_) { auto it = nodes.find(c); if (it != nodes.end() && !it->second.parent.empty()) { it->second.parent = up; node(up); } }
            }
            continue;
        }
        if (k.name == "BRANPROP") {
            // OPM / ECLIPSE: VFP table 0 removes the branch, 9999 = flow line without pressure loss; item 4 ALQ (e.g.
            // compressor power / pump), item 5 ALQ = surface density of oil / gas instead
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.size() < 3 || r[0].empty()) continue;
                const std::string d = nUpper(nUnq(r[0])), u = nUpper(nUnq(r[1]));
                const int vfp = nDefault(r, 2) ? 0 : std::atoi(r[2].c_str());
                NetNode& D = node(d);
                if (vfp == 0) { if (D.parent == u) D.parent.clear(); continue; }
                node(u);
                D.parent = u; D.noLoss = vfp == 9999; D.vfp = D.noLoss ? 0 : vfp;
                D.alq = nDefault(r, 3) ? 0.0 : std::atof(r[3].c_str());
                D.alqEq = nAlqEq(r, 4);
                any = true;
            }
        } else if (k.name == "NODEPROP") {
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.empty() || r[0].empty()) continue;
                NetNode& N = node(nUpper(nUnq(r[0])));
                N.fixed = !nDefault(r, 1);
                if (N.fixed) N.p = std::atof(r[1].c_str());
                N.choke = nYes(r, 2, false);
                N.addLift = nYes(r, 3, false);
                N.chokeGroup = nDefault(r, 4) ? std::string() : nUpper(nUnq(r[4]));
                any = true;
            }
        } else if (k.name == "GRUPNET") {
            // standard network (OPM handleGRUPNET): a terminal pressure >= 0 makes a fixed node (joined to its GRUPTREE
            // parent by a no-loss branch), a VFP table > 0 a branch to the parent, a defaulted / 0 table takes the group out
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.empty() || r[0].empty()) continue;
                const std::string pat = nUnq(r[0]);
                const bool term = !nDefault(r, 1) && std::atof(r[1].c_str()) >= 0;
                const int vfp = nDefault(r, 2) ? 0 : std::atoi(r[2].c_str());
                int hits = 0;
                for (const std::string& g : groups) {
                    if (!nMatch(pat, g)) continue;
                    ++hits;
                    const std::string up = parentOf(g);
                    if (term) {
                        NetNode& N = node(g);
                        N.fixed = true; N.p = std::atof(r[1].c_str());
                        N.addLift = nYes(r, 5, false);
                        if (!up.empty()) { node(up); N.parent = up; N.noLoss = true; N.vfp = 0; }
                        any = true;
                    } else if (vfp <= 0) {
                        auto it = nodes.find(g);
                        if (it != nodes.end() && it->second.parent == up) it->second.parent.clear();
                    } else if (!up.empty()) {
                        NetNode& N = node(g);
                        node(up);
                        N.parent = up; N.vfp = vfp == 9999 ? 0 : vfp; N.noLoss = vfp == 9999; N.fixed = false;
                        N.alq = nDefault(r, 3) ? 0.0 : std::atof(r[3].c_str());
                        N.addLift = nYes(r, 5, false);
                        N.alqEq = nAlqEq(r, 6);
                        any = true;
                    }
                }
                if (!hits && step == k.step) log_("Warning: GRUPNET - no group matches '" + pat + "'");
            }
            netStd_ = true;
        } else if (k.name == "NETCOMPA") {   // C36 automatic compressors (ECLIPSE NETCOMPA items 1-12)
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.size() < 2 || r[0].empty()) continue;
                const std::string in = nUpper(nUnq(r[0])), out = nUpper(nUnq(r[1]));
                NetNode& X = node(in);
                if (!X.parent.empty() && X.parent != out && step == k.step) log_("Warning: NETCOMPA " + in + " -> " + out + ": no such branch (the branch of " + in + " goes to " + X.parent + ")");
                X.compGroup = nDefault(r, 2) ? in : nUpper(nUnq(r[2]));
                X.compPhase = nDefault(r, 3) ? std::string("GAS") : nUpper(nUnq(r[3]));
                X.compVfp = nDefault(r, 4) ? 0 : std::max(0, std::atoi(r[4].c_str()));
                X.compAlq = nDefault(r, 5) ? 0.0 : std::atof(r[5].c_str());
                X.compGas = nDefault(r, 6) ? 0.0 : std::atof(r[6].c_str());
                X.compExtract = nDefault(r, 7) ? in : nUpper(nUnq(r[7]));
                const bool multi = !nDefault(r, 8) && nUpper(nUnq(r[8])).rfind("MULT", 0) == 0;
                X.compLevels = multi && !nDefault(r, 9) ? std::max(1, std::atoi(r[9].c_str())) : 1;
                X.compAlq1 = multi && !nDefault(r, 10) ? std::atof(r[10].c_str()) : X.compAlq;
                X.compSeq = nDefault(r, 11) ? 0 : std::atoi(r[11].c_str());
            }
        } else if (k.name == "NEFAC") {
            for (const auto& r : verbatimRecords(k.body)) {
                if (r.empty() || r[0].empty()) continue;
                node(nUpper(nUnq(r[0]))).nefac = nDefault(r, 1) ? 1.0 : std::clamp(std::atof(r[1].c_str()), 0.0, 1.0);
            }
        }
    }
    const bool had = !netNodes_.empty();
    std::map<std::string, double> keepChoke;
    std::map<std::string, int> keepLevel;
    std::map<std::string, double> keepP;   // C36: node pressures held between NETBALAN balances
    for (const auto& N : netNodes_) { keepChoke[N.name] = N.dpChoke; keepLevel[N.name] = N.compLevel; keepP[N.name] = N.p; }
    netNodes_.clear();
    if (!any) {
        if (had) { for (auto& w : wells_) w.thpNet = -1; log_("Surface network removed"); }
        return;
    }
    for (const std::string& n : order) netNodes_.push_back(nodes[n]);
    for (auto& N : netNodes_) { N.qo = N.qw = N.qg = 0; auto it = keepChoke.find(N.name); if (it != keepChoke.end()) N.dpChoke = it->second;
                                   auto ip = keepP.find(N.name); if (!N.fixed && ip != keepP.end()) N.p = ip->second; }
    int nComp = 0;
    for (auto& N : netNodes_) {
        if (N.compVfp < 0) continue;
        ++nComp;
        auto it = keepLevel.find(N.name);
        N.compLevel = compOff || it == keepLevel.end() ? 0 : std::min(it->second, N.compLevels);
    }
    if (compOff && nComp) log_(fmt("COMPOFF: %g automatic compressor(s) switched off", double(nComp)));
    int terminals = 0;
    for (const auto& N : netNodes_) if (N.parent.empty() && N.fixed) ++terminals;
    if (!had || step == 0)
        log_(std::string("Surface network (") + (netStd_ ? "standard, GRUPNET" : "extended, NODEPROP / BRANPROP") +
             fmt("): %g nodes, %g fixed-pressure terminal node(s), balanced every time step", double(netNodes_.size()), double(terminals)));
    if (terminals == 0) log_("NOTE: network without a fixed-pressure terminal node (NODEPROP item 2) - not balanced");
}

// node pressures from the current (latest converged) well rates; returns false when no node is known. Rates of a leaf
// (source) node = efficiency-weighted production of all wells below its group + satellite production (+ lift gas of
// its own wells / satellite with ADD_GAS_LIFT_GAS); a node passes its inflow minus the GCONSUMP net consumption taken at
// it up its branch (times NEFAC). Pressures top-down from the fixed-pressure nodes (OPM NetworkPressureComputation).
bool Simulator::networkPressures(std::vector<double>* pOut, std::vector<int>* wnodeOut) {
    if (netNodes_.empty()) return false;
    const size_t nn = netNodes_.size();
    std::map<std::string, int> idx;
    for (size_t i = 0; i < nn; ++i) idx[netNodes_[i].name] = int(i);
    std::vector<int> par(nn, -1);
    std::vector<char> hasKids(nn, 0);
    for (size_t i = 0; i < nn; ++i)
        if (!netNodes_[i].parent.empty()) { auto it = idx.find(netNodes_[i].parent); if (it != idx.end()) { par[i] = it->second; hasKids[size_t(it->second)] = 1; } }
    auto connected = [&](size_t i) { return par[i] >= 0 || hasKids[i]; };
    // network node of every producer (THP limit): the first connected node walking up its group chain
    std::vector<int> wnode(wells_.size(), -1);
    for (size_t w = 0; w < wells_.size(); ++w) {
        if (wells_[w].type != WellType::Producer) continue;
        std::string g = nUpper(case_.wells[size_t(wells_[w].input)].group);
        for (int guard = 0; guard < 64 && !g.empty(); ++guard) {
            auto it = idx.find(g);
            if (it != idx.end() && connected(size_t(it->second))) { wnode[w] = it->second; break; }
            if (g == "FIELD") break;
            std::string up = "FIELD";
            for (const auto& [c, p] : gParent_) if (nUpper(c) == g) { up = nUpper(p); break; }
            g = up;
        }
    }
    std::vector<NetNode> N = netNodes_;
    for (auto& X : N) { X.qo = X.qw = X.qg = 0; }
    for (size_t i = 0; i < nn; ++i) {
        if (hasKids[i] || !connected(i)) continue;   // leaf (source) node: the group's production
        NetNode& X = N[i];
        for (size_t w = 0; w < wells_.size(); ++w) {
            const WellRT& W = wells_[w];
            if (W.type != WellType::Producer || !W.open) continue;
            const std::string& wg = case_.wells[size_t(W.input)].group;
            if (!groupBelow(wg, X.name)) continue;
            X.qo += W.eff * W.rates.oil; X.qw += W.eff * W.rates.water; X.qg += W.eff * W.rates.gas;
            if (X.addLift && nUpper(wg) == X.name) X.qg += W.eff * liftGasRate(W);
        }
        X.qo += satRate(X.name, 0); X.qw += satRate(X.name, 1); X.qg += satRate(X.name, 2);
        if (X.addLift) { auto it = gsat_.find(X.name); if (it != gsat_.end()) X.qg += it->second.lift; }
    }
    // GCONSUMP net consumption per node: item 4 node, 'NO-NODES' none, default: the nearest connected node at or above
    // the group if that is a source node
    std::vector<double> cons(nn, 0.0);
    for (const auto& [g, C] : gcons_) {
        const double net = C.cons - C.imp;
        if (net == 0 || C.node == "NO-NODES") continue;
        if (!C.node.empty()) { auto it = idx.find(C.node); if (it != idx.end() && connected(size_t(it->second))) cons[size_t(it->second)] += net; continue; }
        std::string x = g;
        for (int guard = 0; guard < 64 && !x.empty(); ++guard) {
            auto it = idx.find(x);
            if (it != idx.end() && connected(size_t(it->second))) { if (!hasKids[size_t(it->second)]) cons[size_t(it->second)] += net; break; }
            if (x == "FIELD") break;
            std::string up = "FIELD";
            for (const auto& [c, p] : gParent_) if (nUpper(c) == x) { up = nUpper(p); break; }
            x = up;
        }
    }
    // C36: compressor fuel gas taken at the node of the extraction group (nearest connected node at or above it)
    for (size_t c = 0; c < nn; ++c) {
        const NetNode& X = netNodes_[c];
        if (X.compVfp < 0 || X.compLevel <= 0 || X.compGas <= 0) continue;
        const double gas = X.compGas * double(X.compLevel) / double(X.compLevels);
        std::string x = X.compExtract;
        for (int guard = 0; guard < 64 && !x.empty(); ++guard) {
            auto it = idx.find(x);
            if (it != idx.end() && connected(size_t(it->second))) { cons[size_t(it->second)] += gas; break; }
            if (x == "FIELD") break;
            std::string up = "FIELD";
            for (const auto& [ch, pp] : gParent_) if (nUpper(ch) == x) { up = nUpper(pp); break; }
            x = up;
        }
    }
    // roll the rates up the tree (leaves first), the consumption leaves at its node
    std::vector<int> depth(nn, 0);
    for (size_t i = 0; i < nn; ++i) { int d = 0; for (int j = par[i]; j >= 0 && d < int(nn); j = par[size_t(j)]) ++d; depth[i] = d; }
    std::vector<size_t> byDepth(nn);
    for (size_t i = 0; i < nn; ++i) byDepth[i] = i;
    std::sort(byDepth.begin(), byDepth.end(), [&](size_t a, size_t b) { return depth[a] > depth[b]; });
    for (size_t i : byDepth) {
        N[i].qg = std::max(0.0, N[i].qg - cons[i]);
        if (par[i] >= 0) { NetNode& P = N[size_t(par[i])]; const double e = N[i].nefac; P.qo += e * N[i].qo; P.qw += e * N[i].qw; P.qg += e * N[i].qg; }
    }
    // pressures top-down; every fixed-pressure node is known (also one joined to a parent by a no-loss branch)
    std::vector<char> known(nn, 0);
    for (size_t i = 0; i < nn; ++i) if (N[i].fixed && (connected(i) || par[i] < 0)) known[i] = 1;
    std::vector<size_t> topDown(byDepth.rbegin(), byDepth.rend());
    const DensitySC& rho = case_.props.density;
    for (size_t i : topDown) {
        if (known[i] || par[i] < 0 || !known[size_t(par[i])]) continue;
        NetNode& X = N[i];
        const double pu = N[size_t(par[i])].p;
        double pd = pu;
        const bool compOn = X.compVfp >= 0 && X.compLevel > 0;   // C36 NETCOMPA
        const int vfpUse = compOn && X.compVfp > 0 ? X.compVfp : X.vfp;
        if (vfpUse > 0 && !(X.noLoss && !compOn)) {
            const VfpTable* t = nullptr;
            for (const VfpTable& v : case_.vfp) if (v.id == vfpUse && !v.inj) t = &v;
            if (!t) { if (vfpMissing_.insert("NET:" + X.name + std::to_string(vfpUse)).second) log_(fmt("Warning: network branch VFP table %g not found - no pressure loss", double(vfpUse))); }
            else {
                const double qo = X.qo, qw = X.qw, qg = X.qg;
                const double flo = t->flo == "LIQ" ? qo + qw : t->flo == "GAS" ? qg : qo;
                const double wfr = t->wfr == "WOR" ? (qo > 0 ? qw / qo : 0) : t->wfr == "WGR" ? (qg > 0 ? qw / qg : 0) : (qo + qw > 0 ? qw / (qo + qw) : 0);
                const double gfr = t->gfr == "GLR" ? (qo + qw > 0 ? qg / (qo + qw) : 0) : t->gfr == "OGR" ? (qg > 0 ? qo / qg : 0) : (qo > 0 ? qg / qo : 0);
                double alq = X.alqEq == 1 ? rho.oil : X.alqEq == 2 ? rho.gas : X.alq;
                if (compOn) alq = X.compLevels > 1 ? X.compAlq1 + (X.compAlq - X.compAlq1) * double(X.compLevel - 1) / double(X.compLevels - 1) : X.compAlq;
                pd = t->eval(flo, pu, wfr, gfr, alq);
                if (!std::isfinite(pd) || pd <= U_.atm) pd = pu;   // invalid lookup: the upstream pressure (OPM)
            }
        }
        if (X.choke) {   // C35 automatic choke (NODEPROP AS_CHOKE): the node pressure that meets the group's rate target
            const double pc = chokePressure(i, pu, wnode);
            if (pc > pd) pd = pc;
            netNodes_[i].dpChoke = pd - pu;
        }
        X.p = pd; known[i] = 1;
    }
    std::vector<double> p(nn, -1.0);
    for (size_t i = 0; i < nn; ++i) {
        if (known[i]) p[i] = N[i].p;
        netNodes_[i].qo = N[i].qo; netNodes_[i].qw = N[i].qw; netNodes_[i].qg = N[i].qg;   // rates of the last evaluation
    }
    if (pOut) *pOut = p;
    if (wnodeOut) *wnodeOut = wnode;
    return true;
}

// apply node pressures (p[i] < 0: unknown) as THP limits of the producers of each node
void Simulator::networkApply(const std::vector<double>& p, const std::vector<int>& wnode) {
    const double g = U_.gravity;
    for (size_t i = 0; i < netNodes_.size() && i < p.size(); ++i) if (p[i] >= 0) netNodes_[i].p = p[i];
    for (size_t w = 0; w < wells_.size(); ++w) {
        WellRT& W = wells_[w];
        const int nd = w < wnode.size() ? wnode[w] : -1;
        if (nd < 0 || size_t(nd) >= p.size() || p[size_t(nd)] < 0) { W.thpNet = -1; continue; }
        W.thpNet = p[size_t(nd)];
        if (W.own.vfp <= 0) {   // no lift table: node pressure + mixture head from the surface (depth 0) to the datum
            const double lim = W.thpNet + std::max(0.0, W.headRho) * g * std::max(0.0, W.refDepth);
            const double deck = W.bhpDeck > 0 ? W.bhpDeck : U_.atm;
            W.ctrl.bhp = std::max(deck, lim);
        }
    }
}

// C35 automatic choke: the common tubing-head pressure p of the producers below node i at which their (efficiency-weighted)
// rate in the controlled phase equals the group target (GCONPROD of CHOKE_GROUP, default the node's group); the THP
// operating points at the current reservoir state give q(p), decreasing in p -> bisection. Returns < 0 when the target
// is not binding (q(pu) <= target: choke fully open) or the group has no rate target.
double Simulator::chokePressure(size_t i, double pu, const std::vector<int>& wnode) const {
    const NetNode& X = netNodes_[i];
    const std::string grp = X.chokeGroup.empty() ? X.name : X.chokeGroup;
    const GroupProdControl* G = nullptr;
    for (const auto& [n, g] : gprod_) if (nUpper(n) == grp) G = &g;
    if (!G) return -1;
    const Control m = G->mode;
    const double T = m == Control::ORAT ? G->orat : m == Control::WRAT ? G->wrat : m == Control::GRAT ? G->grat : m == Control::LRAT ? G->lrat : -1.0;
    if (!(T > 0)) return -1;
    std::vector<size_t> ws;
    for (size_t w = 0; w < wells_.size(); ++w) {
        const WellRT& W = wells_[w];
        if (W.type != WellType::Producer || !W.open || W.own.vfp <= 0) continue;
        if (w < wnode.size() && wnode[w] == int(i)) ws.push_back(w);
        else if (groupBelow(case_.wells[size_t(W.input)].group, X.name)) ws.push_back(w);
    }
    if (ws.empty()) return -1;
    auto rate = [&](double p) {
        double q = 0;
        for (size_t w : ws) {
            WellRT W = wells_[w];
            W.thpNet = p; W.own.thp = 0;
            double b = 0, r[3] = {0, 0, 0};
            if (!thpOperatingPoint(W, W.alq, W.alqLift, b, r)) continue;
            double f = 1.0;   // the well's own rate limits cap the deliverable rates
            if (W.own.orat > 0 && r[PO] > 0) f = std::min(f, W.own.orat / r[PO]);
            if (W.own.wrat > 0 && r[PW] > 0) f = std::min(f, W.own.wrat / r[PW]);
            if (W.own.grat > 0 && r[PG] > 0) f = std::min(f, W.own.grat / r[PG]);
            if (W.own.lrat > 0 && r[PO] + r[PW] > 0) f = std::min(f, W.own.lrat / (r[PO] + r[PW]));
            const double v = f * (m == Control::ORAT ? r[PO] : m == Control::WRAT ? r[PW] : m == Control::GRAT ? r[PG] : r[PO] + r[PW]);
            q += W.eff * std::max(0.0, v);
        }
        return q;
    };
    if (rate(pu) <= T) return -1;
    double lo = pu, hi = std::max(2.0 * pu, pu + 10.0 * U_.atm);
    for (int k = 0; k < 40 && rate(hi) > T; ++k) { lo = hi; hi *= 1.5; }
    for (int k = 0; k < 50 && hi - lo > 1e-6 * hi; ++k) { const double mid = 0.5 * (lo + hi); (rate(mid) > T ? lo : hi) = mid; }
    if (std::getenv("RESSIM_NETDBG")) {
        double qa = 0; for (size_t w : ws) qa += wells_[w].eff * wellRateOf(wells_[w], m);
        std::fprintf(stderr, "CHOKE %s t %.2f pu %.3f q(pu) %.1f T %.1f -> p %.3f q(p) %.1f  actual %.1f\n", X.name.c_str(), time_, pu, rate(pu), T, hi, rate(hi), qa);
    }
    return hi;
}

bool Simulator::chokeGroup(const std::string& g) const {
    const std::string u = nUpper(g);
    for (const NetNode& N : netNodes_) {
        if (!N.choke || N.parent.empty()) continue;
        if ((N.chokeGroup.empty() ? N.name : N.chokeGroup) == u) return true;
    }
    return false;
}

// C36 NETCOMPA switching at the start of a time step: the compressors of the lowest sequence number whose group misses its
// production target (GCONPROD rate of the compressor's phase, efficiency-weighted group rate below 99.9 % of it) go up one
// level (off -> on for single-level ones); the next sequence follows on later time steps if the target is still missed
void Simulator::networkCompressors() {
    if (time_ <= 0) return;   // no rates before the first time step
    int best = INT_MAX;
    std::vector<size_t> cand;
    for (size_t i = 0; i < netNodes_.size(); ++i) {
        const NetNode& X = netNodes_[i];
        if (X.compVfp < 0 || X.compLevel >= X.compLevels) continue;
        const GroupProdControl* G = nullptr;
        for (const auto& [n, g] : gprod_) if (nUpper(n) == X.compGroup) G = &g;
        if (!G) continue;
        const std::string& ph = X.compPhase;
        const double T = ph == "OIL" ? G->orat : ph == "WAT" ? G->wrat : ph == "LIQ" ? G->lrat : G->grat;
        if (!(T > 0)) continue;
        double q = 0;
        for (const WellRT& W : wells_) {
            if (W.type != WellType::Producer || !W.open || !groupBelow(case_.wells[size_t(W.input)].group, X.compGroup)) continue;
            q += W.eff * (ph == "OIL" ? W.rates.oil : ph == "WAT" ? W.rates.water : ph == "LIQ" ? W.rates.oil + W.rates.water : W.rates.gas);
        }
        if (q >= 0.999 * T) continue;
        if (X.compSeq < best) { best = X.compSeq; cand.clear(); }
        if (X.compSeq == best) cand.push_back(i);
    }
    for (size_t i : cand) {
        NetNode& X = netNodes_[i];
        ++X.compLevel;
        log_(fmt("  t = %.2f d: compressor on branch ", time_) + X.name + " -> " + X.parent +
             (X.compLevels > 1 ? fmt(" up to level %g", double(X.compLevel)) + fmt(" of %g", double(X.compLevels)) : std::string(" switched on")) +
             " (group " + X.compGroup + " below its " + X.compPhase + " target)");
    }
}

double Simulator::compressorGas() const {
    double g = 0;
    for (const NetNode& X : netNodes_) if (X.compVfp >= 0 && X.compLevel > 0) g += X.compGas * double(X.compLevel) / double(X.compLevels);
    return g;
}

void Simulator::networkBalance() {
    std::vector<double> p; std::vector<int> wnode;
    if (!networkPressures(&p, &wnode)) return;
    networkApply(p, wnode);
}

}  // namespace ressim
